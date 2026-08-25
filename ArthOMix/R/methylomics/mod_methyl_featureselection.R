## R/methylomics/mod_methyl_featureselection.R
## Methylomics sub-module: ML Feature Selection - a comprehensive, leakage-
## aware CpG/DNA-methylation feature-selection workflow. Upload your own
## beta/M-value matrix (+ optional phenotype/annotation files) or use the
## preloaded whole-blood cohort; filter probes on methylation-specific QC
## (missingness, variance, SNP/sex-chromosome/cross-reactive/MAF); run up to
## five independent selection methods (Univariate, Regularization, Tree-
## Based, RFE, Stability Selection), each behind its own explicit "Run"
## button; combine them into an overlap-based Consensus panel with a
## configurable count/weighted threshold and optional correlation-based
## redundancy reduction; inspect the final CpGs interactively; validate the
## panel with a leakage-safe or frozen-panel cross-validation; and save the
## whole configuration (not just the CpG list) as a reproducible RDS.
##
## Isolated strictly to this file/module: id = "featureselection" stays as
## already registered in MX_MODULES (submodules_registry.R); every filter/
## parsing/annotation/plotting helper below is reused, never modified, from
## qc.R / parse_upload.R / annotation.R / global.R. Nothing here touches the
## Transcriptomics Feature Selection module or any other Methylomics tab.

mod_methyl_featureselection_config <- list(
  id = "featureselection", title = "ML Feature Selection", icon = "sliders", group = "Biomarker modeling",
  description = "CpG feature selection using univariate, regularization, tree-based, RFE, and stability selection, each run on its own step. Combines results into a consensus set on the preloaded cohort or your own uploaded matrix."
)

## =============================================================================
## Helpers - data source / scale / grouping
## =============================================================================

## .rds companion to qc.R's own methyl_parse_matrix() (CSV/TSV only) - same
## list(ok, mat, error) shape, so callers don't need to branch on file type.
methyl_fs_parse_matrix_rds <- function(datapath) {
  x <- tryCatch(readRDS(datapath), error = function(e) e)
  if (inherits(x, "error")) {
    return(list(ok = FALSE, error = paste("Could not read this .rds file:", conditionMessage(x))))
  }
  if (is.matrix(x) && is.numeric(x)) return(list(ok = TRUE, mat = x))
  if (is.data.frame(x) && ncol(x) >= 2) {
    ids <- as.character(x[[1]])
    if (any(duplicated(ids))) {
      return(list(ok = FALSE, error = sprintf("%d duplicated probe ID(s) in the first column.", sum(duplicated(ids)))))
    }
    m <- tryCatch({
      mm <- as.matrix(x[, -1, drop = FALSE]); storage.mode(mm) <- "double"; rownames(mm) <- ids; mm
    }, error = function(e) NULL)
    if (is.null(m)) return(list(ok = FALSE, error = "Every column after the first must be numeric."))
    return(list(ok = TRUE, mat = m))
  }
  list(ok = FALSE, error = "The .rds file must contain either a numeric probe x sample matrix, or a data.frame with probe IDs in the first column.")
}

methyl_fs_parse_matrix_any <- function(datapath, filename) {
  if (grepl("\\.rds$", filename, ignore.case = TRUE)) methyl_fs_parse_matrix_rds(datapath)
  else methyl_parse_matrix(datapath, filename)
}

## Beta-vs-M-value heuristic: samples up to 20,000 finite values and checks
## whether the 0.1-99.9% quantile range sits inside [-0.05, 1.05]. Only used
## for uploads - the preloaded dataset's own input_scale is authoritative.
methyl_fs_detect_scale <- function(mat) {
  v <- mat[is.finite(mat)]
  if (length(v) == 0) return(list(scale = "beta", frac_out_of_range = 0))
  if (length(v) > 20000) v <- v[sample.int(length(v), 20000)]
  q <- stats::quantile(v, c(0.001, 0.999), na.rm = TRUE)
  scale <- if (q[1] >= -0.05 && q[2] <= 1.05) "beta" else "m"
  frac_out <- mean(v < -0.05 | v > 1.05)
  list(scale = scale, frac_out_of_range = frac_out)
}

methyl_fs_guess_group_col <- function(cols) {
  hit <- intersect(c("group", "Group", "disease", "Disease", "status", "Status", "phenotype", "Phenotype"), cols)
  if (length(hit) > 0) hit[1] else cols[1]
}

methyl_fs_cap_top_variance <- function(mat, n) {
  if (nrow(mat) <= n) return(mat)
  v <- methyl_row_vars(mat)
  top <- order(v, decreasing = TRUE)[seq_len(n)]
  mat[sort(top), , drop = FALSE]
}

## IQR companion to qc.R's methyl_filter_variance()/methyl_filter_sd() -
## same list(keep, note) shape, kept in this module since IQR-based
## thresholding is specific to this workflow's filter set.
methyl_fs_filter_iqr <- function(mat, min_iqr = 0) {
  q <- matrixStats::rowQuantiles(mat, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[, 2] - q[, 1]
  iqr[!is.finite(iqr)] <- 0
  keep <- iqr >= min_iqr
  list(keep = keep, note = sprintf("%d probe(s) below the IQR threshold (%.4g).", sum(!keep), min_iqr))
}

## Sample-axis missingness (no equivalent in qc.R, which is probe-axis only).
methyl_fs_sample_missing_ok <- function(mat, max_na_frac = 1) {
  frac <- colMeans(is.na(mat))
  keep <- frac <= max_na_frac
  list(keep = keep, note = sprintf("%d sample(s) exceed %.0f%% missing probes.", sum(!keep), max_na_frac * 100))
}

## Simple per-probe median imputation - the always-available default; k-NN
## (impute::impute.knn) is used instead when that optional package happens
## to be installed, degrading silently back to median otherwise (no new
## dependency is required for imputation to work).
methyl_fs_impute <- function(mat, method = "median") {
  if (identical(method, "knn") && requireNamespace("impute", quietly = TRUE)) {
    out <- tryCatch(impute::impute.knn(mat)$data, error = function(e) NULL)
    if (!is.null(out)) return(list(mat = out, method_used = "k-NN (impute::impute.knn)"))
  }
  med <- matrixStats::rowMedians(mat, na.rm = TRUE)
  out <- mat
  na_idx <- which(is.na(out), arr.ind = TRUE)
  if (nrow(na_idx) > 0) out[na_idx] <- med[na_idx[, 1]]
  list(mat = out, method_used = "per-probe median")
}

## =============================================================================
## Helpers - Univariate Selection
## =============================================================================

## Linear-model engine (t-test / moderated t-test / ANOVA / linear
## regression / Pearson) - all are just different design/contrast choices
## on the same limma::lmFit()/eBayes() fit, run on M-values, optionally with
## covariate adjustment (age/sex/batch/other metadata columns).
methyl_fs_univariate_linear <- function(m_mat, y, covariates = NULL, mode = c("moderated_t", "t_test", "anova", "linear_regression", "pearson")) {
  mode <- match.arg(mode)
  continuous <- mode %in% c("linear_regression", "pearson")
  if (continuous) {
    yv <- suppressWarnings(as.numeric(y))
    validate(need(!anyNA(yv), "The phenotype column must be numeric for linear regression / Pearson correlation."))
    df <- if (is.null(covariates)) data.frame(y = yv) else cbind(data.frame(y = yv), covariates)
    design <- stats::model.matrix(~., data = df)
    coef_name <- "y"
  } else {
    yf <- factor(as.character(y))
    validate(need(nlevels(yf) >= 2, "The group column needs at least two levels."))
    df <- if (is.null(covariates)) data.frame(grp = yf) else cbind(data.frame(grp = yf), covariates)
    design <- stats::model.matrix(~grp + ., data = df)
    coef_name <- grep("^grp", colnames(design), value = TRUE)
  }
  validate(need(qr(design)$rank == ncol(design), "The selected group/covariate combination is rank-deficient - remove a covariate or use fewer levels."))
  fit <- limma::lmFit(m_mat, design)
  fit <- limma::eBayes(fit, trend = identical(mode, "moderated_t"))
  tt <- if (identical(mode, "anova") && length(coef_name) > 1) {
    limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "F")
  } else {
    limma::topTable(fit, coef = coef_name[1], number = Inf, sort.by = "P")
  }
  data.frame(cpg = rownames(tt), statistic = if ("t" %in% colnames(tt)) tt$t else tt$F,
             p = tt$P.Value, fdr = tt$adj.P.Val, logfc = tt$logFC %||% NA_real_,
             row.names = NULL, stringsAsFactors = FALSE)
}

## Rank-based engine (Wilcoxon / Kruskal-Wallis / Spearman) - no covariate
## adjustment, base R at candidate-panel scale (already <= the Data & Filters
## cap, so a per-probe apply() is fast here unlike a genome-wide DMP scan).
methyl_fs_univariate_rank <- function(mat, y, mode = c("wilcoxon", "kruskal", "spearman")) {
  mode <- match.arg(mode)
  if (identical(mode, "spearman")) {
    yv <- suppressWarnings(as.numeric(y))
    validate(need(!anyNA(yv), "The phenotype column must be numeric for Spearman correlation."))
    res <- apply(mat, 1, function(r) {
      ct <- tryCatch(stats::cor.test(r, yv, method = "spearman"), error = function(e) NULL)
      if (is.null(ct)) c(rho = NA_real_, p = NA_real_) else c(rho = unname(ct$estimate), p = ct$p.value)
    })
    df <- data.frame(cpg = rownames(mat), statistic = res["rho", ], p = res["p", ], row.names = NULL)
  } else {
    yf <- factor(as.character(y))
    if (identical(mode, "wilcoxon")) {
      validate(need(nlevels(yf) == 2, "Wilcoxon rank-sum needs exactly two groups - use Kruskal-Wallis for more than two."))
      res <- apply(mat, 1, function(r) {
        wt <- tryCatch(stats::wilcox.test(r ~ yf), error = function(e) NULL)
        if (is.null(wt)) c(stat = NA_real_, p = NA_real_) else c(stat = unname(wt$statistic), p = wt$p.value)
      })
    } else {
      validate(need(nlevels(yf) >= 2, "Kruskal-Wallis needs at least two groups."))
      res <- apply(mat, 1, function(r) {
        kt <- tryCatch(stats::kruskal.test(r, yf), error = function(e) NULL)
        if (is.null(kt)) c(stat = NA_real_, p = NA_real_) else c(stat = unname(kt$statistic), p = kt$p.value)
      })
    }
    df <- data.frame(cpg = rownames(mat), statistic = res["stat", ], p = res["p", ], row.names = NULL)
  }
  df$fdr <- stats::p.adjust(df$p, method = "BH")
  df
}

## Applies the chosen selection rule (FDR / p-value / |dbeta| / top-N) on top
## of either engine's ranked table, and attaches mean-Beta/dbeta/direction
## effect-size columns computed directly from the beta-scale matrix.
methyl_fs_univariate_select <- function(ranked_df, beta_mat, y, rule = "fdr", threshold = 0.05, dbeta_min = 0, top_n = 50) {
  df <- ranked_df[order(ranked_df$p), , drop = FALSE]
  yf <- suppressWarnings(as.numeric(y))
  is_continuous <- !anyNA(yf) && length(unique(as.character(y))) > 6
  if (!is_continuous) {
    grp <- factor(as.character(y))
    lv <- levels(grp)
    mean1 <- rowMeans(beta_mat[df$cpg, grp == lv[1], drop = FALSE], na.rm = TRUE)
    mean2 <- rowMeans(beta_mat[df$cpg, grp == lv[length(lv)], drop = FALSE], na.rm = TRUE)
    df$mean_beta_ref <- mean1; df$mean_beta_comp <- mean2
    df$dbeta <- mean2 - mean1
    df$direction <- ifelse(df$dbeta > 0, "hyper", "hypo")
  } else {
    df$mean_beta_ref <- NA_real_; df$mean_beta_comp <- NA_real_
    df$dbeta <- df$logfc %||% NA_real_
    df$direction <- NA_character_
  }
  keep <- switch(rule,
    fdr    = !is.na(df$fdr) & df$fdr <= threshold,
    pvalue = !is.na(df$p) & df$p <= threshold,
    top_n  = seq_len(nrow(df)) <= top_n,
    rep(TRUE, nrow(df))
  )
  if (!is.na(dbeta_min) && dbeta_min > 0 && !is_continuous) keep <- keep & !is.na(df$dbeta) & abs(df$dbeta) >= dbeta_min
  list(ranked = df, selected_ids = df$cpg[keep])
}

## =============================================================================
## Helpers - Regularization (LASSO / Elastic Net)
## =============================================================================

methyl_fs_lasso_fit <- function(X, y, alpha = 1, nfolds = 10, lambda_choice = "lambda.min", nlambda = 100,
                                 standardize = TRUE, weights = NULL, seed = 1234, max_selected = NULL, coef_threshold = 0) {
  set.seed(seed)
  cv <- tryCatch(
    glmnet::cv.glmnet(X, y, family = "binomial", alpha = alpha, nfolds = max(3, min(nfolds, min(table(y)))),
                       nlambda = nlambda, standardize = standardize, weights = weights),
    error = function(e) validate(need(FALSE, paste("glmnet failed to fit:", conditionMessage(e))))
  )
  lambda_s <- if (identical(lambda_choice, "lambda.1se")) "lambda.1se" else "lambda.min"
  co <- stats::coef(cv, s = lambda_s)[-1, 1, drop = TRUE]
  co <- co[abs(co) > coef_threshold]
  co <- sort(co[co != 0], decreasing = TRUE)
  ids <- names(co)
  if (!is.null(max_selected) && length(ids) > max_selected) {
    ord <- order(abs(co), decreasing = TRUE)[seq_len(max_selected)]
    ids <- names(co)[ord]; co <- co[ord]
  }
  list(cv = cv, lambda_used = unname(cv[[lambda_s]]), coefficients = co, selected_ids = ids,
       ranked = data.frame(cpg = names(co), coefficient = as.numeric(co), abs_coefficient = abs(as.numeric(co)),
                            rank = seq_along(co), row.names = NULL))
}

## =============================================================================
## Helpers - Tree-Based Selection (Random Forest)
## =============================================================================

methyl_fs_rf_fit <- function(X, y, ntree = 1000, mtry_mode = "auto", mtry_manual = NULL, nodesize = 1,
                              maxnodes = NULL, sample_fraction = NULL, replace = TRUE, classwt = NULL,
                              seed = 1234, importance_type = "gini", rule = "top_n", top_n = 50, threshold = NULL) {
  p <- ncol(X)
  if (identical(mtry_mode, "manual") && !is.null(mtry_manual)) {
    mtry <- min(p, max(1, round(mtry_manual)))
  } else {
    mtry <- max(1, floor(sqrt(p)))
  }
  ## randomForest::randomForest() has non-NULL internal defaults for
  ## sampsize/maxnodes/classwt (missing(arg)-style checks) - passing NULL
  ## explicitly rather than omitting the argument breaks it ("missing value
  ## where TRUE/FALSE needed"), so unset optional args are only added to the
  ## call when actually provided.
  rf_args <- list(x = X, y = y, importance = TRUE, ntree = max(50, round(ntree)), mtry = mtry,
                   nodesize = max(1, round(nodesize)), replace = replace)
  if (!is.null(maxnodes)) rf_args$maxnodes <- maxnodes
  if (!is.null(classwt)) rf_args$classwt <- classwt
  if (!is.null(sample_fraction)) rf_args$sampsize <- pmax(1, round(sample_fraction * table(y)))
  set.seed(seed)
  rf <- tryCatch(
    do.call(randomForest::randomForest, rf_args),
    error = function(e) validate(need(FALSE, paste("randomForest failed to fit:", conditionMessage(e))))
  )
  imp_col <- if (identical(importance_type, "accuracy")) "MeanDecreaseAccuracy" else "MeanDecreaseGini"
  imp <- sort(rf$importance[, imp_col], decreasing = TRUE)
  if (identical(rule, "threshold") && !is.null(threshold)) {
    ids <- names(imp)[imp >= threshold]
  } else if (identical(rule, "percentile") && !is.null(threshold)) {
    ids <- names(imp)[imp >= stats::quantile(imp, threshold, na.rm = TRUE)]
  } else {
    ids <- head(names(imp), min(p, max(1, round(top_n))))
  }
  list(rf = rf, mtry = mtry, importance_type = imp_col, importance = imp, selected_ids = ids,
       ranked = data.frame(cpg = names(imp), importance = as.numeric(imp), rank = seq_along(imp), row.names = NULL))
}

## =============================================================================
## Helpers - RFE / Wrapper Selection
## =============================================================================

methyl_fs_rfe_sizes <- function(sizes_str, p) {
  sizes <- suppressWarnings(as.integer(trimws(strsplit(sizes_str %||% "", ",")[[1]])))
  sizes <- sort(unique(sizes[!is.na(sizes) & sizes > 0 & sizes <= p]))
  if (length(sizes) == 0) sizes <- sort(unique(pmin(p, c(10, 25, 50, 100, 250, 500, 1000, 2000, 5000))))
  sizes
}

methyl_fs_rfe_run <- function(X, y, sizes, k = 5, flavor = c("rf", "logistic"), seed = 1234) {
  flavor <- match.arg(flavor)
  fn <- if (identical(flavor, "rf")) caret::rfFuncs else caret::lrFuncs
  ctrl <- caret::rfeControl(functions = fn, method = "cv", number = max(3, k), verbose = FALSE)
  set.seed(seed)
  fit <- tryCatch(caret::rfe(x = as.data.frame(X), y = y, sizes = sizes, rfeControl = ctrl),
                   error = function(e) validate(need(FALSE, paste("caret::rfe failed:", conditionMessage(e)))))
  perf <- fit$results[, c("Variables", if ("Accuracy" %in% colnames(fit$results)) "Accuracy" else colnames(fit$results)[2])]
  colnames(perf) <- c("size", "performance")
  list(fit = fit, sizes = sizes, optimal_size = fit$optsize, selected_ids = caret::predictors(fit),
       perf = perf, ranked = data.frame(cpg = fit$optVariables, rank = seq_along(fit$optVariables), row.names = NULL))
}

## Ported near-verbatim from the transcriptomics module's own SVM-RFE
## (mod_featureselection.R fs_svm_rfe_rank()/fs_svm_rfe_curve()) - backward
## elimination by smallest squared linear-SVM weight (one feature per round,
## most-to-least-important ranking), then a k-fold CV error curve over the
## top-k ranked features to pick the panel size objectively.
methyl_fs_svm_rfe_rank <- function(X, y, cost = 1, tolerance = 0.001) {
  feats <- colnames(X)
  ranking <- character(0)
  while (length(feats) > 1) {
    m <- e1071::svm(X[, feats, drop = FALSE], y, kernel = "linear", scale = TRUE, cost = cost, tolerance = tolerance)
    w2 <- ((t(m$coefs) %*% m$SV)[1, ])^2
    drop <- names(sort(w2))[1]
    ranking <- c(drop, ranking)
    feats <- setdiff(feats, drop)
  }
  c(feats, ranking)
}

methyl_fs_svm_rfe_curve <- function(X, y, rank, cost = 1, seed = 1234, folds = 5, tolerance = 0.001) {
  ks <- seq_along(rank)
  err <- vapply(ks, function(k) {
    set.seed(seed)
    acc <- e1071::svm(X[, rank[seq_len(k)], drop = FALSE], y, kernel = "linear", scale = TRUE,
                       cost = cost, cross = folds, tolerance = tolerance)$tot.accuracy
    1 - acc / 100
  }, numeric(1))
  list(k = ks, err = err, best = ks[which.min(err)], besterr = min(err))
}

methyl_fs_svm_rfe_run <- function(X, y, cost = 1, tolerance = 0.001, k_folds = 5, seed = 1234,
                                   panel_mode = "auto", manual_k = 10) {
  set.seed(seed)
  rank <- methyl_fs_svm_rfe_rank(X, y, cost = cost, tolerance = tolerance)
  curve <- methyl_fs_svm_rfe_curve(X, y, rank, cost = cost, seed = seed, folds = k_folds, tolerance = tolerance)
  k_sel <- if (identical(panel_mode, "manual")) min(length(rank), max(1, round(manual_k))) else curve$best
  ids <- rank[seq_len(k_sel)]
  list(rank = rank, curve = curve, optimal_size = k_sel, selected_ids = ids,
       ranked = data.frame(cpg = rank, rank = seq_along(rank), row.names = NULL))
}

## =============================================================================
## Helpers - Stability Selection
## =============================================================================

## One generalized resampling loop (bootstrap / repeated k-fold /
## subsampling), all driving a fixed LASSO base selector at a shared
## reference lambda (Meinshausen-Buhlmann stability selection) - only how
## the resample index sets are drawn differs between the three modes.
methyl_fs_stability_resample_indices <- function(n, y, type = c("bootstrap", "repeated_kfold", "subsampling"),
                                                   n_resamples = 50, fraction = 0.8, k = 5, repeats = 5, seed = 1234) {
  type <- match.arg(type)
  set.seed(seed)
  if (identical(type, "bootstrap")) {
    lapply(seq_len(n_resamples), function(i) sample.int(n, n, replace = TRUE))
  } else if (identical(type, "subsampling")) {
    sz <- max(2, round(fraction * n))
    lapply(seq_len(n_resamples), function(i) sample.int(n, sz, replace = FALSE))
  } else {
    folds <- caret::createMultiFolds(y, k = max(2, k), times = max(1, repeats))
    folds
  }
}

methyl_fs_stability_run <- function(X, y, type = "bootstrap", n_resamples = 50, fraction = 0.8, k = 5, repeats = 5,
                                     seed = 1234, freq_threshold = 0.7, alpha = 1) {
  set.seed(seed)
  ref_cv <- tryCatch(glmnet::cv.glmnet(X, y, family = "binomial", alpha = alpha, nfolds = max(3, min(10, min(table(y))))),
                      error = function(e) validate(need(FALSE, paste("Reference LASSO fit failed:", conditionMessage(e)))))
  ref_lambda <- ref_cv$lambda.min
  idx_sets <- methyl_fs_stability_resample_indices(nrow(X), y, type = type, n_resamples = n_resamples,
                                                     fraction = fraction, k = k, repeats = repeats, seed = seed)
  sel_mat <- matrix(0L, nrow = ncol(X), ncol = length(idx_sets), dimnames = list(colnames(X), NULL))
  for (i in seq_along(idx_sets)) {
    idx <- idx_sets[[i]]
    if (length(unique(as.character(y[idx]))) < 2) next
    fit <- tryCatch(glmnet::glmnet(X[idx, , drop = FALSE], y[idx], family = "binomial", alpha = alpha, lambda = ref_lambda),
                     error = function(e) NULL)
    if (is.null(fit)) next
    co <- stats::coef(fit)[-1, 1]
    sel_mat[names(co)[co != 0], i] <- 1L
  }
  freq <- rowMeans(sel_mat)
  mean_rank <- apply(sel_mat, 1, function(r) if (sum(r) == 0) NA_real_ else mean(which(r == 1)))
  stability_score <- freq
  ranked <- data.frame(cpg = names(freq), selection_frequency = as.numeric(freq),
                        stability_score = as.numeric(stability_score), row.names = NULL)
  ranked <- ranked[order(-ranked$selection_frequency), , drop = FALSE]
  list(n_resamples = length(idx_sets), ref_lambda = ref_lambda, ranked = ranked,
       selected_ids = ranked$cpg[ranked$selection_frequency >= freq_threshold])
}

## =============================================================================
## Helpers - Consensus / Overlap / Correlation reduction
## =============================================================================

methyl_fs_consensus_table <- function(id_lists, weights = NULL) {
  all_ids <- Reduce(union, id_lists)
  methods <- names(id_lists)
  if (length(all_ids) == 0) return(list(table = data.frame(cpg = character(0)), all_ids = character(0)))
  w <- weights %||% stats::setNames(rep(1, length(methods)), methods)
  w <- w[methods]
  tbl <- data.frame(cpg = all_ids, stringsAsFactors = FALSE)
  for (m in methods) tbl[[m]] <- as.integer(tbl$cpg %in% id_lists[[m]])
  flag_mat <- as.matrix(tbl[, methods, drop = FALSE])
  tbl$n_methods <- rowSums(flag_mat)
  tbl$weighted_score <- as.numeric(flag_mat %*% matrix(w, ncol = 1)) / sum(w)
  tbl <- tbl[order(-tbl$weighted_score, -tbl$n_methods), , drop = FALSE]
  list(table = tbl, all_ids = all_ids, methods = methods, weights = w)
}

methyl_fs_consensus_select <- function(consensus_table, min_methods = 2, min_weighted = NULL) {
  keep <- consensus_table$n_methods >= min_methods
  if (!is.null(min_weighted)) keep <- keep | consensus_table$weighted_score >= min_weighted
  consensus_table$cpg[keep]
}

## Manual UpSet-style plot (dot-matrix intersection panel + intersection-
## size bar chart, combined via patchwork) - no UpSetR/ComplexUpset
## dependency needed for up to ~5 sets.
methyl_fs_upset_plot <- function(id_lists) {
  id_lists <- id_lists[lengths(id_lists) > 0]
  validate(need(length(id_lists) >= 2, "Need at least two non-empty methods for an UpSet plot."))
  methods <- names(id_lists)
  all_ids <- Reduce(union, id_lists)
  mem <- vapply(id_lists, function(s) all_ids %in% s, logical(length(all_ids)))
  combo_key <- apply(mem, 1, function(r) paste(as.integer(r), collapse = ""))
  combo_tbl <- as.data.frame(table(combo_key), stringsAsFactors = FALSE)
  combo_tbl <- combo_tbl[order(-combo_tbl$Freq), , drop = FALSE]
  combo_tbl <- utils::head(combo_tbl, 20)
  combo_tbl$combo_id <- factor(combo_tbl$combo_key, levels = combo_tbl$combo_key)

  bar <- ggplot(combo_tbl, aes(x = combo_id, y = Freq)) +
    geom_col(fill = ARTHOMIX_COLORS$blue) +
    geom_text(aes(label = Freq), vjust = -0.3, size = 3) +
    labs(x = NULL, y = "Intersection size") + theme_arthomix() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  dot_df <- do.call(rbind, lapply(seq_len(nrow(combo_tbl)), function(i) {
    bits <- as.integer(strsplit(combo_tbl$combo_key[i], "")[[1]])
    data.frame(combo_id = combo_tbl$combo_id[i], method = methods, member = bits == 1)
  }))
  dot_df$method <- factor(dot_df$method, levels = rev(methods))
  dots <- ggplot(dot_df, aes(x = combo_id, y = method)) +
    geom_point(aes(color = member, size = member)) +
    scale_color_manual(values = c("TRUE" = ARTHOMIX_COLORS$ink, "FALSE" = ARTHOMIX_COLORS$grid), guide = "none") +
    scale_size_manual(values = c("TRUE" = 3.2, "FALSE" = 1.6), guide = "none") +
    labs(x = "Intersection", y = NULL) + theme_arthomix() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank())

  patchwork::wrap_plots(bar, dots, ncol = 1, heights = c(2, 1.2))
}

methyl_fs_consensus_rank_plot <- function(consensus_table, rank_by = "weighted_score", n = 30) {
  df <- utils::head(consensus_table[order(-consensus_table[[rank_by]]), , drop = FALSE], n)
  df$cpg <- factor(df$cpg, levels = rev(df$cpg))
  ggplot(df, aes(x = .data[[rank_by]], y = cpg)) +
    geom_col(fill = ARTHOMIX_COLORS$aqua) +
    labs(x = tools::toTitleCase(gsub("_", " ", rank_by)), y = NULL) + theme_arthomix()
}

methyl_fs_method_heatmap <- function(consensus_table, methods) {
  df <- consensus_table[, c("cpg", methods), drop = FALSE]
  df <- utils::head(df[order(-rowSums(as.matrix(df[, methods, drop = FALSE]))), , drop = FALSE], 60)
  long <- stats::reshape(df, direction = "long", varying = methods, v.names = "selected",
                          times = methods, timevar = "method", idvar = "cpg")
  long$cpg <- factor(long$cpg, levels = rev(unique(df$cpg)))
  ggplot(long, aes(x = method, y = cpg, fill = factor(selected))) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("0" = ARTHOMIX_COLORS$grid, "1" = ARTHOMIX_COLORS$blue), guide = "none") +
    labs(x = NULL, y = NULL) + theme_arthomix() +
    theme(axis.text.y = element_text(size = 6))
}

## Greedy pairwise-correlation pruning of a (small) candidate panel: repeatedly
## finds the highest-|r| pair above `r_threshold` and drops whichever member
## has the lower consensus score, until no pair exceeds the threshold.
methyl_fs_correlation_reduce <- function(mat, ids, score, method = "pearson", r_threshold = 0.8) {
  ids <- intersect(ids, rownames(mat))
  if (length(ids) < 2) return(list(reduced_ids = ids, dropped = data.frame(kept = character(0), dropped = character(0), r = numeric(0))))
  sub <- t(mat[ids, , drop = FALSE])
  cmat <- stats::cor(sub, method = method, use = "pairwise.complete.obs")
  diag(cmat) <- 0
  scores <- stats::setNames(score[ids] %||% rep(0, length(ids)), ids)
  active <- ids
  dropped <- list()
  repeat {
    sub_c <- abs(cmat[active, active, drop = FALSE])
    if (length(active) < 2 || max(sub_c) <= r_threshold) break
    idx <- which(sub_c == max(sub_c), arr.ind = TRUE)[1, ]
    a <- active[idx[1]]; b <- active[idx[2]]
    loser <- if (scores[a] >= scores[b]) b else a
    keeper <- setdiff(c(a, b), loser)
    dropped[[length(dropped) + 1]] <- data.frame(kept = keeper, dropped = loser, r = cmat[a, b])
    active <- setdiff(active, loser)
  }
  list(reduced_ids = active, dropped = if (length(dropped) > 0) do.call(rbind, dropped) else data.frame(kept = character(0), dropped = character(0), r = numeric(0)))
}

## =============================================================================
## Helpers - Selected Features / annotation
## =============================================================================

## Same chr/pos/gene join pattern as mod_methyl_dmp.R:648-655. CpG-island /
## promoter-enhancer / genomic-region columns are intentionally left NA -
## the shared annotation.R (used by other methylomics modules, not modified
## here) only exposes chr/pos/Type/SNP-rs/gene.
methyl_fs_annotate_panel <- function(ids, anno_result) {
  df <- data.frame(cpg = ids, chr = NA_character_, pos = NA_real_, gene = NA_character_,
                    cpg_island = NA_character_, genomic_region = NA_character_, stringsAsFactors = FALSE)
  if (isTRUE(anno_result$ok)) {
    a <- anno_result$anno
    hit <- df$cpg %in% rownames(a)
    df$chr[hit] <- a[df$cpg[hit], "chr"]
    df$pos[hit] <- a[df$cpg[hit], "pos"]
    df$gene[hit] <- a[df$cpg[hit], "gene"]
  }
  df
}

## =============================================================================
## Helpers - Validation
## =============================================================================

methyl_fs_validate_frozen <- function(X, y, classifier = "glm", k = 5, repeats = 1, seed = 1234) {
  method <- switch(classifier, rf = "rf", svm = "svmLinear", "glm")
  ctrl <- caret::trainControl(method = if (repeats > 1) "repeatedcv" else "cv", number = max(3, k),
                               repeats = max(1, repeats), classProbs = TRUE, summaryFunction = caret::twoClassSummary,
                               savePredictions = "final")
  set.seed(seed)
  fit <- tryCatch(caret::train(x = as.data.frame(X), y = y, method = method, trControl = ctrl, metric = "ROC"),
                   error = function(e) validate(need(FALSE, paste("Validation model failed to fit:", conditionMessage(e)))))
  preds <- fit$pred
  roc_obj <- tryCatch(pROC::roc(preds$obs, preds[[levels(y)[2]]], quiet = TRUE), error = function(e) NULL)
  list(fit = fit, mean_auc = suppressMessages(tryCatch(as.numeric(pROC::auc(roc_obj)), error = function(e) NA_real_)),
       roc = roc_obj, resample_results = fit$resample)
}

## Manual outer k-fold loop reselecting the panel (Univariate + Regularization
## only, per the plan) inside each training fold - leakage-safe, mirroring
## the outer/inner CV shape used elsewhere in this codebase's diagnostic
## module (caret::createFolds() + a fixed-hyperparameter per-fold refit).
methyl_fs_validate_nested <- function(beta_mat, m_mat, y, uni_params, lasso_params, classifier = "glm",
                                       outer_k = 5, repeats = 1, seed = 1234) {
  set.seed(seed)
  yf <- factor(as.character(y))
  per_fold <- list()
  for (rep_i in seq_len(max(1, repeats))) {
    folds <- caret::createFolds(yf, k = max(2, outer_k), list = TRUE)
    for (fi in seq_along(folds)) {
      test_idx <- folds[[fi]]; train_idx <- setdiff(seq_along(yf), test_idx)
      y_tr <- yf[train_idx]; y_te <- yf[test_idx]
      if (length(unique(as.character(y_tr))) < 2 || length(unique(as.character(y_te))) < 2) next
      uni <- tryCatch(methyl_fs_univariate_linear(m_mat[, train_idx, drop = FALSE], y_tr, mode = "moderated_t"), error = function(e) NULL)
      if (is.null(uni)) next
      sel <- methyl_fs_univariate_select(uni, beta_mat[, train_idx, drop = FALSE], y_tr,
                                          rule = uni_params$rule %||% "top_n", threshold = uni_params$threshold %||% 0.05,
                                          top_n = uni_params$top_n %||% 100)
      ids <- sel$selected_ids
      if (length(ids) < 2) ids <- utils::head(uni[order(uni$p), "cpg"], 20)
      Xtr <- t(m_mat[ids, train_idx, drop = FALSE]); Xte <- t(m_mat[ids, test_idx, drop = FALSE])
      lasso <- tryCatch(methyl_fs_lasso_fit(Xtr, y_tr, alpha = lasso_params$alpha %||% 1, nfolds = 3, seed = seed), error = function(e) NULL)
      panel_ids <- if (!is.null(lasso) && length(lasso$selected_ids) >= 2) lasso$selected_ids else ids
      method <- switch(classifier, rf = "rf", svm = "svmLinear", "glm")
      ctrl0 <- caret::trainControl(method = "none", classProbs = TRUE)
      set.seed(seed)
      fit <- tryCatch(caret::train(x = as.data.frame(Xtr[, panel_ids, drop = FALSE]), y = y_tr, method = method, trControl = ctrl0), error = function(e) NULL)
      if (is.null(fit)) next
      prob <- tryCatch(stats::predict(fit, newdata = as.data.frame(Xte[, panel_ids, drop = FALSE]), type = "prob"), error = function(e) NULL)
      if (is.null(prob)) next
      roc_obj <- tryCatch(pROC::roc(y_te, prob[[levels(yf)[2]]], quiet = TRUE), error = function(e) NULL)
      auc_i <- if (!is.null(roc_obj)) as.numeric(pROC::auc(roc_obj)) else NA_real_
      per_fold[[length(per_fold) + 1]] <- data.frame(repeat_i = rep_i, fold = fi, n_panel = length(panel_ids), auc = auc_i)
    }
  }
  per_fold_df <- if (length(per_fold) > 0) do.call(rbind, per_fold) else data.frame(repeat_i = integer(0), fold = integer(0), n_panel = integer(0), auc = numeric(0))
  list(per_fold = per_fold_df, mean_auc = mean(per_fold_df$auc, na.rm = TRUE), sd_auc = stats::sd(per_fold_df$auc, na.rm = TRUE))
}

## =============================================================================
## UI
## =============================================================================

mod_methyl_featureselection_ui <- function(id) {
  ns <- NS(id)
  div(
    class = "tx-menu-wrap",
    tabsetPanel(
      id = ns("fs_subtabs"), type = "tabs",
      tabPanel("Data & Filters", br(), withSpinner(uiOutput(ns("filters_ui")), color = "#2563EB", type = 6)),
      tabPanel("Univariate Selection", br(), withSpinner(uiOutput(ns("uni_ui")), color = "#2563EB", type = 6)),
      tabPanel("LASSO", br(), withSpinner(uiOutput(ns("reg_ui")), color = "#2563EB", type = 6)),
      tabPanel("Tree-Based Selection", br(), withSpinner(uiOutput(ns("rf_ui")), color = "#2563EB", type = 6)),
      tabPanel("RFE / Wrapper Selection", br(), withSpinner(uiOutput(ns("rfe_ui")), color = "#2563EB", type = 6)),
      tabPanel("Stability Selection", br(), withSpinner(uiOutput(ns("stab_ui")), color = "#2563EB", type = 6)),
      tabPanel("Consensus / Overlap", br(), withSpinner(uiOutput(ns("consensus_ui")), color = "#2563EB", type = 6)),
      tabPanel("Selected Features", br(), withSpinner(uiOutput(ns("selected_ui")), color = "#2563EB", type = 6)),
      tabPanel("Model & Export", br(), withSpinner(uiOutput(ns("export_ui")), color = "#2563EB", type = 6))
    )
  )
}

mod_methyl_fs_empty_note <- function(msg = "Configure the parameters and click Run to generate results.") {
  div(class = "empty-note", icon("circle-info"), msg)
}

mod_methyl_fs_need_filters_note <- function() {
  div(class = "empty-note", icon("circle-info"), "Run Data & Filters first (see the Data & Filters tab).")
}

## =============================================================================
## Server
## =============================================================================

mod_methyl_featureselection_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## ---------------------------------------------------------------------
    ## Stage 0: source resolution (always live, never Run-gated - a plain
    ## dataset/filter summary is fine to show immediately; only actual
    ## computation is gated behind explicit Run buttons).
    ## ---------------------------------------------------------------------

    fs_own_matrix <- reactive({
      req(input$fs_matrix_file)
      p <- methyl_fs_parse_matrix_any(input$fs_matrix_file$datapath, input$fs_matrix_file$name)
      validate(need(isTRUE(p$ok), p$error %||% "Could not parse the uploaded matrix."))
      p$mat
    })

    fs_own_sheet <- reactive({
      if (is.null(input$fs_sheet_file)) return(NULL)
      p <- methyl_parse_sample_sheet(input$fs_sheet_file$datapath, input$fs_sheet_file$name)
      validate(need(isTRUE(p$ok), p$error %||% "Could not parse the uploaded phenotype file."))
      p$df
    })

    fs_exclusion_ids <- reactive({
      if (is.null(input$fs_exclusion_file)) return(NULL)
      p <- methyl_parse_probe_list(input$fs_exclusion_file$datapath, input$fs_exclusion_file$name)
      if (!isTRUE(p$ok)) return(NULL)
      p$ids
    })

    fs_maf_table <- reactive({
      if (is.null(input$fs_maf_file)) return(NULL)
      p <- methyl_parse_maf_list(input$fs_maf_file$datapath, input$fs_maf_file$name)
      if (!isTRUE(p$ok)) return(NULL)
      p$maf
    })

    fs_active_source <- reactive({
      if (identical(input$fs_source %||% "preloaded", "upload")) {
        mat <- fs_own_matrix()
        det <- if (identical(input$fs_scale_mode %||% "auto", "auto")) methyl_fs_detect_scale(mat) else list(scale = if (identical(input$fs_scale_mode, "beta")) "beta" else "m", frac_out_of_range = mean(mat[is.finite(mat)] < -0.05 | mat[is.finite(mat)] > 1.05))
        list(mat = mat, sheet = fs_own_sheet(), scale = det$scale, frac_out_of_range = det$frac_out_of_range,
             array_type = input$fs_upload_array_type %||% "450K", source_label = "Uploaded matrix", preloaded = FALSE)
      } else {
        validate(need(!is.null(dataset$beta), "No preloaded methylation dataset is available in this deployment - switch to \"Upload my own\" above."))
        v <- dataset$beta[is.finite(dataset$beta)]
        frac_out <- if (length(v) > 0) mean((if (length(v) > 20000) v[sample.int(length(v), 20000)] else v) < -0.05 | (if (length(v) > 20000) v[sample.int(length(v), 20000)] else v) > 1.05) else 0
        list(mat = dataset$beta, sheet = dataset$sample_sheet, scale = dataset$input_scale %||% "beta",
             frac_out_of_range = frac_out, array_type = dataset$array_type %||% "450K",
             source_label = dataset$source %||% "Preloaded", preloaded = TRUE)
      }
    })

    anno_result <- reactive({ methyl_get_annotation(fs_active_source()$array_type) })

    fs_group_choices <- reactive({
      src <- fs_active_source()
      if (is.null(src$sheet)) return(character(0))
      colnames(src$sheet)
    })

    ## ---------------------------------------------------------------------
    ## Stage 1: Data & Filters
    ## ---------------------------------------------------------------------

    output$filters_ui <- renderUI({
      tagList(
        div(class = "card",
          div(class = "card-title", icon("filter"), "Data source"),
          radioButtons(ns("fs_source"), NULL,
                       choices = if (METH_DATA_AVAILABLE) c("Preloaded whole-blood cohort" = "preloaded", "Upload my own" = "upload") else c("Upload my own" = "upload"),
                       selected = if (METH_DATA_AVAILABLE) "preloaded" else "upload", inline = TRUE),
          conditionalPanel(condition = sprintf("input['%s'] == 'upload'", ns("fs_source")),
            fluidRow(
              column(6,
                p(strong("Methylation matrix"), " - CSV/TSV/TXT/RDS, CpG/probe IDs in the first column (or rownames), samples as columns."),
                fileInput(ns("fs_matrix_file"), "Beta/M-value matrix", accept = c(".csv", ".tsv", ".txt", ".rds", ".Rds")),
                selectInput(ns("fs_upload_array_type"), "Array type (for annotation-dependent filters)", choices = METHYL_ARRAY_TYPES, selected = "450K"),
                radioButtons(ns("fs_scale_mode"), "Value scale", choices = c("Auto-detect" = "auto", "Force beta [0,1]" = "beta", "Force M-value" = "m"), selected = "auto", inline = TRUE)
              ),
              column(6,
                p(strong("Phenotype / sample metadata"), " (optional but needed for group-based methods) - CSV/TSV, one row per sample."),
                fileInput(ns("fs_sheet_file"), "Phenotype file", accept = c(".csv", ".tsv", ".txt")),
                p(strong("Optional: cross-reactive probe exclusion list"), " - one probe ID per line."),
                fileInput(ns("fs_exclusion_file"), "Exclusion list", accept = c(".csv", ".tsv", ".txt")),
                p(strong("Optional: MAF table"), " - probe_id,maf CSV/TSV."),
                fileInput(ns("fs_maf_file"), "MAF table", accept = c(".csv", ".tsv", ".txt"))
              )
            )
          ),
          uiOutput(ns("fs_source_summary_ui"))
        ),
        uiOutput(ns("fs_filter_params_ui")),
        uiOutput(ns("fs_filter_result_ui"))
      )
    })

    output$fs_source_summary_ui <- renderUI({
      src <- tryCatch(fs_active_source(), error = function(e) NULL)
      if (is.null(src)) return(mod_methyl_fs_empty_note("Load a dataset above."))
      tagList(
        p(class = "empty-note", icon("circle-info"),
          sprintf("%s: %s probes x %s samples, detected scale = %s.", src$source_label,
                  format(nrow(src$mat), big.mark = ","), ncol(src$mat), toupper(src$scale))),
        if (src$frac_out_of_range > 0.001)
          p(class = "empty-note", icon("triangle-exclamation"),
            sprintf("%.1f%% of sampled values fall outside the expected beta range [0, 1] (+/-0.05 tolerance) - confirm this is really a beta-value matrix, or set Value scale to \"Force M-value\" above.", src$frac_out_of_range * 100))
        else NULL
      )
    })

    output$fs_filter_params_ui <- renderUI({
      src <- tryCatch(fs_active_source(), error = function(e) NULL)
      req(src)
      cols <- fs_group_choices()
      arr_ok <- isTRUE(anno_result()$ok)
      div(class = "card",
        div(class = "card-title", icon("sliders"), "Filters"),
        fluidRow(
          column(4,
            tags$h5("Grouping"),
            if (length(cols) == 0) p(class = "empty-note", icon("triangle-exclamation"), "No phenotype/sample metadata loaded - group-dependent methods will be unavailable.")
            else selectInput(ns("fs_group_col"), "Group / phenotype column", choices = cols, selected = methyl_fs_guess_group_col(cols)),
            uiOutput(ns("fs_group_levels_ui"))
          ),
          column(4,
            tags$h5("Missingness"),
            numericInput(ns("fs_max_na_probe"), "Max missingness per CpG (%)", value = 5, min = 0, max = 100, step = 5),
            numericInput(ns("fs_max_na_sample"), "Max missingness per sample (%)", value = 20, min = 0, max = 100, step = 5),
            checkboxInput(ns("fs_impute"), "Impute remaining missing values", value = TRUE),
            radioButtons(ns("fs_impute_method"), NULL, choices = c("Median (default)" = "median", "k-NN (if available)" = "knn"), selected = "median", inline = TRUE)
          ),
          column(4,
            tags$h5("Variance filtering"),
            radioButtons(ns("fs_var_metric"), "Metric", choices = c("Variance" = "variance", "SD" = "sd", "IQR" = "iqr"), selected = "variance", inline = TRUE),
            numericInput(ns("fs_var_min"), "Minimum threshold", value = 0, min = 0, step = 0.001),
            selectInput(ns("fs_max_probes"), "Maximum candidate CpGs to carry forward (top by variance)",
                        choices = c("500" = 500, "1,000" = 1000, "5,000" = 5000, "10,000" = 10000, "Custom" = "custom"), selected = 5000),
            conditionalPanel(condition = sprintf("input['%s'] == 'custom'", ns("fs_max_probes")),
                              numericInput(ns("fs_max_probes_custom"), NULL, value = 5000, min = 500, max = 50000, step = 500))
          )
        ),
        fluidRow(
          column(4,
            tags$h5("Methylation-specific filters"),
            checkboxInput(ns("fs_filter_snp"), "Remove SNP-associated probes", value = FALSE),
            checkboxInput(ns("fs_filter_sexchr"), "Remove sex-chromosome probes", value = FALSE),
            if (!arr_ok) p(class = "empty-note", icon("triangle-exclamation"), "This filter requires compatible CpG/probe annotation.") else NULL
          ),
          column(4,
            tags$h5("Cross-reactive / MAF"),
            checkboxInput(ns("fs_filter_crossreactive"), "Apply uploaded cross-reactive exclusion list", value = FALSE),
            checkboxInput(ns("fs_filter_maf"), "Apply uploaded MAF filter", value = FALSE),
            numericInput(ns("fs_max_maf"), "Max MAF", value = 0.05, min = 0, max = 0.5, step = 0.01)
          ),
          column(4,
            tags$h5("Beta <-> M-value"),
            p(class = "submodule-desc", "Beta values are convenient for interpretation (0-1 scale); M-values are the better-behaved scale for statistical modeling. This module derives M-values internally wherever a method needs them - the displayed matrix stays on whichever scale was detected/selected above.")
          )
        ),
        actionButton(ns("filters_run_btn"), "Run Filters", icon = icon("play"), class = "btn-primary")
      )
    })

    output$fs_group_levels_ui <- renderUI({
      src <- tryCatch(fs_active_source(), error = function(e) NULL)
      req(src, input$fs_group_col, src$sheet)
      req(input$fs_group_col %in% colnames(src$sheet))
      lv <- sort(unique(as.character(stats::na.omit(src$sheet[[input$fs_group_col]]))))
      if (length(lv) < 2) return(p(class = "empty-note", icon("triangle-exclamation"), "This column has fewer than two distinct values."))
      tagList(
        selectInput(ns("fs_ref_group"), "Reference group", choices = lv, selected = lv[1]),
        selectInput(ns("fs_comp_group"), "Comparison group", choices = lv, selected = lv[min(2, length(lv))])
      )
    })

    filters_has_run <- reactiveVal(FALSE)
    observeEvent(dataset$beta, filters_has_run(FALSE), ignoreNULL = TRUE)
    observeEvent(input$fs_source, filters_has_run(FALSE))

    fs_filter_result <- eventReactive(input$filters_run_btn, {
      src <- fs_active_source()
      sheet <- src$sheet
      validate(need(!is.null(sheet) && !is.null(input$fs_group_col), "Load phenotype/sample metadata and pick a group column first."))
      sample_ids <- methyl_sheet_sample_ids(sheet, colnames(src$mat))
      common <- intersect(colnames(src$mat), sample_ids)
      validate(need(length(common) >= 10, "Fewer than 10 samples match between the matrix and the phenotype file."))
      mat0 <- src$mat[, common, drop = FALSE]
      ph0 <- as.data.frame(sheet)[match(common, sample_ids), , drop = FALSE]
      grp_raw <- as.character(ph0[[input$fs_group_col]])
      keep_s <- !is.na(grp_raw) & grp_raw %in% c(input$fs_ref_group, input$fs_comp_group)
      validate(need(sum(keep_s) >= 10, "Fewer than 10 samples in the chosen reference/comparison groups."))
      mat1 <- mat0[, keep_s, drop = FALSE]
      grp <- factor(grp_raw[keep_s], levels = c(input$fs_ref_group, input$fs_comp_group))
      validate(need(all(table(grp) >= 3), "Each group needs at least 3 samples."))

      beta_scale <- if (identical(src$scale, "m")) 2^mat1 / (1 + 2^mat1) else mat1

      f_samp <- methyl_fs_sample_missing_ok(beta_scale, (input$fs_max_na_sample %||% 20) / 100)
      beta_scale <- beta_scale[, f_samp$keep, drop = FALSE]
      grp <- grp[f_samp$keep]
      validate(need(ncol(beta_scale) >= 6, "Fewer than 6 samples remain after the sample-missingness filter."))

      filters <- list()
      filters[["Probe missingness"]] <- methyl_filter_missing(beta_scale, (input$fs_max_na_probe %||% 5) / 100)
      var_filter <- switch(input$fs_var_metric %||% "variance",
        sd = methyl_filter_sd(beta_scale, input$fs_var_min %||% 0),
        iqr = methyl_fs_filter_iqr(beta_scale, input$fs_var_min %||% 0),
        methyl_filter_variance(beta_scale, input$fs_var_min %||% 0))
      filters[["Variance/SD/IQR"]] <- var_filter
      ar <- anno_result()
      if (isTRUE(input$fs_filter_snp)) filters[["SNP-associated"]] <- methyl_filter_snp(beta_scale, ar)
      if (isTRUE(input$fs_filter_sexchr)) filters[["Sex-chromosome"]] <- methyl_filter_sex_chr(beta_scale, ar, "remove_xy")
      if (isTRUE(input$fs_filter_crossreactive)) filters[["Cross-reactive"]] <- methyl_filter_cross_reactive(beta_scale, fs_exclusion_ids())
      if (isTRUE(input$fs_filter_maf)) filters[["MAF"]] <- methyl_filter_maf(beta_scale, fs_maf_table(), input$fs_max_maf %||% 0.05)

      n_probes_start <- nrow(beta_scale)
      keep_probe <- rep(TRUE, n_probes_start)
      for (f in filters) keep_probe <- keep_probe & f$keep
      validate(need(sum(keep_probe) >= 10,
        sprintf("After filtering, only %d CpGs remain. Reduce the filtering thresholds or increase the number of retained features.", sum(keep_probe))))
      cascade <- methyl_probe_retention_cascade(n_probes_start, filters)
      filter_notes <- vapply(filters, function(f) f$note, character(1))
      beta_scale <- beta_scale[keep_probe, , drop = FALSE]

      ## Top-variance cap: a separate step from the probe filters above (it
      ## subsets rows rather than producing a same-length keep mask), so it's
      ## appended to the cascade table as one extra row instead of folded
      ## into the filters list that cascade/keep_probe were built from.
      cap <- if (identical(input$fs_max_probes %||% "5000", "custom")) (input$fs_max_probes_custom %||% 5000) else as.numeric(input$fs_max_probes %||% 5000)
      n_before_cap <- nrow(beta_scale)
      beta_scale <- methyl_fs_cap_top_variance(beta_scale, round(cap))
      cap_note <- sprintf("Capped to the top %d most-variable CpGs (of %d).", nrow(beta_scale), n_before_cap)
      cascade <- rbind(cascade, data.frame(step = "Top-variance cap", retained = nrow(beta_scale), removed = n_before_cap - nrow(beta_scale)))
      filter_notes <- c(filter_notes, "Top-variance cap" = cap_note)

      imputed_method <- NA_character_
      if (isTRUE(input$fs_impute) && anyNA(beta_scale)) {
        imp <- methyl_fs_impute(beta_scale, input$fs_impute_method %||% "median")
        beta_scale <- imp$mat
        imputed_method <- imp$method_used
      }

      m_scale <- methyl_beta_to_mvalue(beta_scale)

      list(beta = beta_scale, m = m_scale, grp = grp, ref_group = input$fs_ref_group, comp_group = input$fs_comp_group,
           group_col = input$fs_group_col, array_type = src$array_type, source_label = src$source_label,
           n_probes = nrow(beta_scale), n_samples = ncol(beta_scale), cascade = cascade,
           filter_notes = filter_notes, imputed_method = imputed_method,
           anno = anno_result(), run_at = Sys.time())
    })

    observeEvent(input$filters_run_btn, filters_has_run(TRUE), ignoreInit = TRUE)

    output$fs_filter_result_ui <- renderUI({
      if (!filters_has_run()) return(mod_methyl_fs_empty_note())
      r <- fs_filter_result()
      tagList(
        div(class = "card",
          div(class = "card-title", icon("circle-check"), "Filter result"),
          div(class = "methyl-stats-row", fluidRow(
            valueBox(format(r$n_probes, big.mark = ","), "CpGs retained", icon = icon("dna"), color = "purple", width = 4),
            valueBox(r$n_samples, sprintf("Samples (%s vs %s)", r$comp_group, r$ref_group), icon = icon("users"), color = "light-blue", width = 4),
            valueBox(if (is.na(r$imputed_method)) "None" else r$imputed_method, "Imputation", icon = icon("wand-magic-sparkles"), color = "black", width = 4)
          )),
          tags$ul(lapply(r$filter_notes, tags$li)),
          withSpinner(plotOutput(ns("fs_cascade_plot"), height = 260), color = "#2563EB", type = 6)
        )
      )
    })

    output$fs_cascade_plot <- renderPlot({ req(filters_has_run()); methyl_plot_cascade(fs_filter_result()$cascade) })

    ## ---------------------------------------------------------------------
    ## Stage 2: Univariate Selection
    ## ---------------------------------------------------------------------

    output$uni_ui <- renderUI({
      if (!filters_has_run()) return(mod_methyl_fs_need_filters_note())
      div(class = "card",
        div(class = "card-title", icon("chart-simple"), "Univariate Selection"),
        fluidRow(
          column(4,
            selectInput(ns("uni_method"), "Statistical method",
                        choices = c("Moderated t-test (limma)" = "moderated_t", "t-test" = "t_test", "ANOVA (>2 groups)" = "anova",
                                    "Wilcoxon rank-sum" = "wilcoxon", "Kruskal-Wallis" = "kruskal",
                                    "Linear regression (continuous phenotype)" = "linear_regression",
                                    "Pearson correlation (continuous)" = "pearson", "Spearman correlation (continuous)" = "spearman"),
                        selected = "moderated_t"),
            conditionalPanel(condition = sprintf("['linear_regression','pearson','spearman'].includes(input['%s'])", ns("uni_method")),
              p(class = "empty-note", icon("circle-info"),
                "Data & Filters currently defines a reference/comparison group, not an independent continuous phenotype column - this method runs against that group contrast numerically coded (0/1), equivalent to a two-group test.")),
            uiOutput(ns("uni_covariate_ui"))
          ),
          column(4,
            radioButtons(ns("uni_rule"), "Selection rule", choices = c("FDR threshold" = "fdr", "p-value threshold" = "pvalue", "Top N" = "top_n"), selected = "fdr"),
            numericInput(ns("uni_threshold"), "Threshold", value = 0.05, min = 0, max = 1, step = 0.01),
            numericInput(ns("uni_top_n"), "Top N", value = 100, min = 5, max = 5000, step = 5)
          ),
          column(4,
            numericInput(ns("uni_dbeta_min"), "Minimum |Delta-Beta| (categorical methods only)", value = 0, min = 0, max = 1, step = 0.01)
          )
        ),
        actionButton(ns("uni_run_btn"), "Run Univariate Selection", icon = icon("play"), class = "btn-primary"),
        uiOutput(ns("uni_result_ui"))
      )
    })

    output$uni_covariate_ui <- renderUI({
      r <- tryCatch(fs_filter_result(), error = function(e) NULL)
      if (is.null(r) || is.null(fs_active_source()$sheet)) return(NULL)
      cand <- setdiff(colnames(fs_active_source()$sheet), r$group_col)
      if (length(cand) == 0) return(NULL)
      checkboxGroupInput(ns("uni_covariates"), "Covariate adjustment (linear-model methods only)", choices = cand, selected = character(0))
    })

    uni_has_run <- reactiveVal(FALSE)
    observeEvent(fs_filter_result(), uni_has_run(FALSE))
    observeEvent(input$uni_run_btn, uni_has_run(TRUE), ignoreInit = TRUE)

    fs_uni_result <- eventReactive(input$uni_run_btn, {
      r <- fs_filter_result()
      is_linear <- input$uni_method %in% c("moderated_t", "t_test", "anova", "linear_regression", "pearson")
      cov_df <- NULL
      if (is_linear && length(input$uni_covariates %||% character(0)) > 0) {
        sheet <- fs_active_source()$sheet
        ids <- methyl_sheet_sample_ids(sheet, colnames(r$beta))
        cov_df <- as.data.frame(sheet)[match(colnames(r$beta), ids), input$uni_covariates, drop = FALSE]
      }
      ranked <- if (is_linear) {
        methyl_fs_univariate_linear(r$m, r$grp, covariates = cov_df, mode = input$uni_method)
      } else {
        methyl_fs_univariate_rank(r$beta, r$grp, mode = input$uni_method)
      }
      sel <- methyl_fs_univariate_select(ranked, r$beta, r$grp, rule = input$uni_rule, threshold = input$uni_threshold,
                                          dbeta_min = input$uni_dbeta_min, top_n = input$uni_top_n)
      list(ranked = sel$ranked, selected_ids = sel$selected_ids, method = input$uni_method, run_at = Sys.time())
    })

    output$uni_result_ui <- renderUI({
      if (!uni_has_run()) return(mod_methyl_fs_empty_note())
      r <- fs_uni_result()
      tagList(
        div(class = "card",
          div(class = "card-title", icon("circle-check"), "Univariate result"),
          p(strong(length(r$selected_ids)), sprintf(" of %d CpGs selected (%s).", nrow(r$ranked), r$method)),
          div(class = "table-toolbar", downloadButton(ns("uni_download"), "Ranked table (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns("uni_table"))
        )
      )
    })

    output$uni_table <- DT::renderDataTable({
      req(uni_has_run())
      DT::datatable(fs_uni_result()$ranked, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatSignif(columns = intersect(c("statistic", "p", "fdr", "dbeta", "logfc"), colnames(fs_uni_result()$ranked)), digits = 4)
    })
    outputOptions(output, "uni_table", suspendWhenHidden = FALSE)

    output$uni_download <- downloadHandler(
      filename = function() "univariate_selection.csv",
      content = function(file) utils::write.csv(fs_uni_result()$ranked, file, row.names = FALSE)
    )

    ## ---------------------------------------------------------------------
    ## Stage 3: Regularization (LASSO / Elastic Net)
    ## ---------------------------------------------------------------------

    output$reg_ui <- renderUI({
      if (!filters_has_run()) return(mod_methyl_fs_need_filters_note())
      div(class = "card",
        div(class = "card-title", icon("chart-line"), "Regularization - LASSO / Elastic Net"),
        fluidRow(
          column(4,
            sliderInput(ns("reg_alpha"), "Alpha (1 = LASSO, 0 = Ridge, in-between = Elastic Net)", min = 0, max = 1, value = 1, step = 0.05),
            radioButtons(ns("reg_lambda_choice"), "Lambda", choices = c("lambda.min" = "lambda.min", "lambda.1se (sparser)" = "lambda.1se"), selected = "lambda.min")
          ),
          column(4,
            numericInput(ns("reg_nfolds"), "Cross-validation folds", value = 10, min = 3, max = 10, step = 1),
            numericInput(ns("reg_nlambda"), "Number of lambda values (nlambda)", value = 100, min = 20, max = 300, step = 10),
            numericInput(ns("reg_seed"), "Random seed", value = 1234, min = 1, step = 1)
          ),
          column(4,
            checkboxInput(ns("reg_standardize"), "Standardize predictors", value = TRUE),
            radioButtons(ns("reg_class_weight"), "Class weighting", choices = c("Equal (default)" = "equal", "Balanced (inverse frequency)" = "balanced"), selected = "equal"),
            numericInput(ns("reg_max_selected"), "Max selected CpGs (optional)", value = 200, min = 5, max = 5000, step = 5),
            numericInput(ns("reg_coef_threshold"), "Coefficient magnitude threshold", value = 0, min = 0, step = 0.001)
          )
        ),
        actionButton(ns("reg_run_btn"), "Run LASSO / Elastic Net", icon = icon("play"), class = "btn-primary"),
        uiOutput(ns("reg_result_ui"))
      )
    })

    reg_has_run <- reactiveVal(FALSE)
    observeEvent(fs_filter_result(), reg_has_run(FALSE))
    observeEvent(input$reg_run_btn, reg_has_run(TRUE), ignoreInit = TRUE)

    fs_reg_result <- eventReactive(input$reg_run_btn, {
      r <- fs_filter_result()
      X <- t(r$m)
      w <- if (identical(input$reg_class_weight, "balanced")) {
        tb <- table(r$grp); ww <- max(tb) / tb; unname(ww[as.character(r$grp)])
      } else NULL
      fit <- methyl_fs_lasso_fit(X, r$grp, alpha = input$reg_alpha, nfolds = input$reg_nfolds, lambda_choice = input$reg_lambda_choice,
                                  nlambda = input$reg_nlambda, standardize = input$reg_standardize, weights = w,
                                  seed = input$reg_seed, max_selected = input$reg_max_selected, coef_threshold = input$reg_coef_threshold)
      fit$n_input <- ncol(X); fit$run_at <- Sys.time()
      fit
    })

    output$reg_result_ui <- renderUI({
      if (!reg_has_run()) return(mod_methyl_fs_empty_note())
      r <- fs_reg_result()
      tagList(
        div(class = "card",
          div(class = "card-title", icon("circle-check"), "Regularization result"),
          p(strong(length(r$selected_ids)), sprintf(" of %d candidate CpGs selected (lambda = %.5g).", r$n_input, r$lambda_used)),
          withSpinner(plotOutput(ns("reg_plot"), height = 300), color = "#2563EB", type = 6),
          div(class = "table-toolbar", downloadButton(ns("reg_download"), "Selected CpGs (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns("reg_table"))
        )
      )
    })

    output$reg_plot <- renderPlot({ req(reg_has_run()); plot(fs_reg_result()$cv) })

    output$reg_table <- DT::renderDataTable({
      req(reg_has_run())
      DT::datatable(fs_reg_result()$ranked, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatSignif(columns = c("coefficient", "abs_coefficient"), digits = 4)
    })
    outputOptions(output, "reg_table", suspendWhenHidden = FALSE)

    output$reg_download <- downloadHandler(
      filename = function() "regularization_selected_cpgs.csv",
      content = function(file) utils::write.csv(fs_reg_result()$ranked, file, row.names = FALSE)
    )

    ## ---------------------------------------------------------------------
    ## Stage 4: Tree-Based Selection (Random Forest)
    ## ---------------------------------------------------------------------

    output$rf_ui <- renderUI({
      if (!filters_has_run()) return(mod_methyl_fs_need_filters_note())
      div(class = "card",
        div(class = "card-title", icon("tree"), "Tree-Based Selection - Random Forest"),
        fluidRow(
          column(4,
            numericInput(ns("rf_ntree"), "Number of trees", value = 1000, min = 100, max = 5000, step = 100),
            radioButtons(ns("rf_mtry_mode"), "mtry", choices = c("Auto (sqrt(p))" = "auto", "Manual" = "manual"), selected = "auto"),
            conditionalPanel(condition = sprintf("input['%s'] == 'manual'", ns("rf_mtry_mode")), numericInput(ns("rf_mtry_manual"), "mtry value", value = 10, min = 1, step = 1))
          ),
          column(4,
            numericInput(ns("rf_nodesize"), "Minimum node size", value = 1, min = 1, max = 50, step = 1),
            checkboxInput(ns("rf_maxnodes_unlimited"), "Unlimited tree depth", value = TRUE),
            conditionalPanel(condition = sprintf("!input['%s']", ns("rf_maxnodes_unlimited")), numericInput(ns("rf_maxnodes"), "Max terminal nodes", value = 20, min = 2, step = 1)),
            checkboxInput(ns("rf_replace"), "Sample with replacement", value = TRUE)
          ),
          column(4,
            radioButtons(ns("rf_importance_type"), "Importance metric", choices = c("Mean Decrease Gini" = "gini", "Mean Decrease Accuracy" = "accuracy"), selected = "gini"),
            radioButtons(ns("rf_rule"), "Selection rule", choices = c("Top N" = "top_n", "Threshold" = "threshold", "Percentile" = "percentile"), selected = "top_n"),
            numericInput(ns("rf_top_n"), "N / threshold / percentile value", value = 50, min = 0, step = 1),
            numericInput(ns("rf_seed"), "Random seed", value = 1234, min = 1, step = 1)
          )
        ),
        actionButton(ns("rf_run_btn"), "Run Random Forest", icon = icon("play"), class = "btn-primary"),
        uiOutput(ns("rf_result_ui"))
      )
    })

    rf_has_run <- reactiveVal(FALSE)
    observeEvent(fs_filter_result(), rf_has_run(FALSE))
    observeEvent(input$rf_run_btn, rf_has_run(TRUE), ignoreInit = TRUE)

    fs_rf_result <- eventReactive(input$rf_run_btn, {
      r <- fs_filter_result()
      X <- t(r$m)
      maxnodes <- if (isTRUE(input$rf_maxnodes_unlimited)) NULL else input$rf_maxnodes
      fit <- methyl_fs_rf_fit(X, r$grp, ntree = input$rf_ntree, mtry_mode = input$rf_mtry_mode, mtry_manual = input$rf_mtry_manual,
                               nodesize = input$rf_nodesize, maxnodes = maxnodes, replace = input$rf_replace, seed = input$rf_seed,
                               importance_type = input$rf_importance_type, rule = input$rf_rule, top_n = input$rf_top_n, threshold = input$rf_top_n)
      fit$n_input <- ncol(X); fit$run_at <- Sys.time()
      fit
    })

    output$rf_result_ui <- renderUI({
      if (!rf_has_run()) return(mod_methyl_fs_empty_note())
      r <- fs_rf_result()
      tagList(
        div(class = "card",
          div(class = "card-title", icon("circle-check"), "Random Forest result"),
          p(strong(length(r$selected_ids)), sprintf(" of %d candidate CpGs selected (mtry = %d, %s).", r$n_input, r$mtry, r$importance_type)),
          withSpinner(plotOutput(ns("rf_plot"), height = 340), color = "#2563EB", type = 6),
          div(class = "table-toolbar", downloadButton(ns("rf_download"), "Ranked table (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns("rf_table"))
        )
      )
    })

    output$rf_plot <- renderPlot({
      req(rf_has_run())
      r <- fs_rf_result()
      df <- utils::head(r$ranked, 25)
      ggplot(df, aes(x = reorder(cpg, importance), y = importance)) +
        geom_col(fill = ARTHOMIX_COLORS$blue) + coord_flip() +
        labs(x = NULL, y = r$importance_type) + theme_arthomix()
    })

    output$rf_table <- DT::renderDataTable({
      req(rf_has_run())
      df <- fs_rf_result()$ranked
      df$selected <- df$cpg %in% fs_rf_result()$selected_ids
      DT::datatable(df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatSignif(columns = "importance", digits = 4)
    })
    outputOptions(output, "rf_table", suspendWhenHidden = FALSE)

    output$rf_download <- downloadHandler(
      filename = function() "random_forest_ranked_cpgs.csv",
      content = function(file) utils::write.csv(fs_rf_result()$ranked, file, row.names = FALSE)
    )

    ## ---------------------------------------------------------------------
    ## Stage 5: RFE / Wrapper Selection
    ## ---------------------------------------------------------------------

    output$rfe_ui <- renderUI({
      if (!filters_has_run()) return(mod_methyl_fs_need_filters_note())
      r <- fs_filter_result()
      div(class = "card",
        div(class = "card-title", icon("arrows-to-dot"), "RFE / Wrapper Selection"),
        fluidRow(
          column(4, selectInput(ns("rfe_flavor"), "Method", choices = c("Random Forest RFE" = "rf", "Logistic Regression RFE" = "logistic", "SVM-RFE" = "svm"), selected = "rf")),
          column(4, textInput(ns("rfe_sizes"), "Feature subset sizes (comma-separated)", value = paste(methyl_fs_rfe_sizes(NULL, r$n_probes), collapse = ", "))),
          column(4, numericInput(ns("rfe_folds"), "Cross-validation folds", value = 5, min = 3, max = 10, step = 1))
        ),
        conditionalPanel(condition = sprintf("input['%s'] == 'svm'", ns("rfe_flavor")),
          fluidRow(
            column(4, numericInput(ns("rfe_svm_cost"), "SVM cost (C)", value = 1, min = 0.001, step = 0.1)),
            column(4, numericInput(ns("rfe_svm_tolerance"), "Convergence tolerance", value = 0.001, min = 0.00001, step = 0.0001)),
            column(4, radioButtons(ns("rfe_svm_panel_mode"), "Panel size", choices = c("Auto: minimize CV error" = "auto", "Manual" = "manual"), selected = "auto"),
                   conditionalPanel(condition = sprintf("input['%s'] == 'manual'", ns("rfe_svm_panel_mode")), numericInput(ns("rfe_svm_manual_k"), "k", value = 10, min = 1, step = 1)))
          )
        ),
        numericInput(ns("rfe_seed"), "Random seed", value = 1234, min = 1, step = 1),
        actionButton(ns("rfe_run_btn"), "Run RFE", icon = icon("play"), class = "btn-primary"),
        uiOutput(ns("rfe_result_ui"))
      )
    })

    rfe_has_run <- reactiveVal(FALSE)
    observeEvent(fs_filter_result(), rfe_has_run(FALSE))
    observeEvent(input$rfe_run_btn, rfe_has_run(TRUE), ignoreInit = TRUE)

    fs_rfe_result <- eventReactive(input$rfe_run_btn, {
      r <- fs_filter_result()
      X <- t(r$m)
      sizes <- methyl_fs_rfe_sizes(input$rfe_sizes, ncol(X))
      fit <- if (identical(input$rfe_flavor, "svm")) {
        methyl_fs_svm_rfe_run(X, r$grp, cost = input$rfe_svm_cost, tolerance = input$rfe_svm_tolerance, k_folds = input$rfe_folds,
                               seed = input$rfe_seed, panel_mode = input$rfe_svm_panel_mode, manual_k = input$rfe_svm_manual_k)
      } else {
        methyl_fs_rfe_run(X, r$grp, sizes = sizes, k = input$rfe_folds, flavor = input$rfe_flavor, seed = input$rfe_seed)
      }
      fit$flavor <- input$rfe_flavor; fit$n_input <- ncol(X); fit$run_at <- Sys.time()
      fit
    })

    output$rfe_result_ui <- renderUI({
      if (!rfe_has_run()) return(mod_methyl_fs_empty_note())
      r <- fs_rfe_result()
      tagList(
        div(class = "card",
          div(class = "card-title", icon("circle-check"), "RFE result"),
          p(strong(length(r$selected_ids)), sprintf(" of %d candidate CpGs selected (optimal size = %s).", r$n_input, r$optimal_size)),
          withSpinner(plotOutput(ns("rfe_plot"), height = 300), color = "#2563EB", type = 6),
          div(class = "table-toolbar", downloadButton(ns("rfe_download"), "Ranked table (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns("rfe_table"))
        )
      )
    })

    output$rfe_plot <- renderPlot({
      req(rfe_has_run())
      r <- fs_rfe_result()
      if (identical(r$flavor, "svm")) {
        df <- data.frame(k = r$curve$k, error = r$curve$err)
        ggplot(df, aes(x = k, y = error)) + geom_line(color = ARTHOMIX_COLORS$blue) + geom_point(color = ARTHOMIX_COLORS$blue) +
          geom_vline(xintercept = r$optimal_size, linetype = "dashed", color = "grey50") +
          labs(x = "Top-k ranked CpGs", y = "CV classification error") + theme_arthomix()
      } else {
        ggplot(r$perf, aes(x = size, y = performance)) + geom_line(color = ARTHOMIX_COLORS$blue) + geom_point(color = ARTHOMIX_COLORS$blue) +
          geom_vline(xintercept = r$optimal_size, linetype = "dashed", color = "grey50") +
          labs(x = "Feature subset size", y = "CV performance") + theme_arthomix()
      }
    })

    output$rfe_table <- DT::renderDataTable({
      req(rfe_has_run())
      DT::datatable(fs_rfe_result()$ranked, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact")
    })
    outputOptions(output, "rfe_table", suspendWhenHidden = FALSE)

    output$rfe_download <- downloadHandler(
      filename = function() "rfe_ranked_cpgs.csv",
      content = function(file) utils::write.csv(fs_rfe_result()$ranked, file, row.names = FALSE)
    )

    ## ---------------------------------------------------------------------
    ## Stage 6: Stability Selection
    ## ---------------------------------------------------------------------

    output$stab_ui <- renderUI({
      if (!filters_has_run()) return(mod_methyl_fs_need_filters_note())
      div(class = "card",
        div(class = "card-title", icon("repeat"), "Stability Selection"),
        p(class = "submodule-desc", "Repeatedly resamples a LASSO base selector and tabulates how often each CpG is chosen - important for methylation data, where hundreds of thousands of CpGs can be highly correlated."),
        fluidRow(
          column(4, radioButtons(ns("stab_type"), "Resampling scheme", choices = c("Bootstrap" = "bootstrap", "Repeated k-fold" = "repeated_kfold", "Subsampling" = "subsampling"), selected = "bootstrap")),
          column(4,
            numericInput(ns("stab_n_resamples"), "Number of resamples", value = 50, min = 10, max = 500, step = 10),
            numericInput(ns("stab_fraction"), "Sample fraction (bootstrap/subsampling)", value = 0.8, min = 0.1, max = 1, step = 0.05)
          ),
          column(4,
            numericInput(ns("stab_k"), "k (repeated k-fold)", value = 5, min = 2, max = 10, step = 1),
            numericInput(ns("stab_repeats"), "Repeats (repeated k-fold)", value = 5, min = 1, max = 20, step = 1),
            numericInput(ns("stab_freq_threshold"), "Minimum selection frequency", value = 0.7, min = 0.05, max = 1, step = 0.05),
            numericInput(ns("stab_seed"), "Random seed", value = 1234, min = 1, step = 1)
          )
        ),
        actionButton(ns("stab_run_btn"), "Run Stability Selection", icon = icon("play"), class = "btn-primary"),
        uiOutput(ns("stab_result_ui"))
      )
    })

    stab_has_run <- reactiveVal(FALSE)
    observeEvent(fs_filter_result(), stab_has_run(FALSE))
    observeEvent(input$stab_run_btn, stab_has_run(TRUE), ignoreInit = TRUE)

    fs_stab_result <- eventReactive(input$stab_run_btn, {
      r <- fs_filter_result()
      X <- t(r$m)
      fit <- methyl_fs_stability_run(X, r$grp, type = input$stab_type, n_resamples = input$stab_n_resamples, fraction = input$stab_fraction,
                                      k = input$stab_k, repeats = input$stab_repeats, seed = input$stab_seed, freq_threshold = input$stab_freq_threshold)
      fit$n_input <- ncol(X); fit$run_at <- Sys.time()
      fit
    })

    output$stab_result_ui <- renderUI({
      if (!stab_has_run()) return(mod_methyl_fs_empty_note())
      r <- fs_stab_result()
      tagList(
        div(class = "card",
          div(class = "card-title", icon("circle-check"), "Stability Selection result"),
          p(strong(length(r$selected_ids)), sprintf(" of %d candidate CpGs selected (>= %.0f%% of %d resamples).", r$n_input, input$stab_freq_threshold * 100, r$n_resamples)),
          withSpinner(plotOutput(ns("stab_plot"), height = 300), color = "#2563EB", type = 6),
          div(class = "table-toolbar", downloadButton(ns("stab_download"), "Ranked table (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns("stab_table"))
        )
      )
    })

    output$stab_plot <- renderPlot({
      req(stab_has_run())
      df <- utils::head(fs_stab_result()$ranked, 30)
      ggplot(df, aes(x = reorder(cpg, selection_frequency), y = selection_frequency)) +
        geom_col(fill = ARTHOMIX_COLORS$aqua) + coord_flip() +
        geom_hline(yintercept = input$stab_freq_threshold, linetype = "dashed", color = "grey50") +
        labs(x = NULL, y = "Selection frequency") + theme_arthomix()
    })

    output$stab_table <- DT::renderDataTable({
      req(stab_has_run())
      DT::datatable(fs_stab_result()$ranked, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatSignif(columns = c("selection_frequency", "stability_score"), digits = 4)
    })
    outputOptions(output, "stab_table", suspendWhenHidden = FALSE)

    output$stab_download <- downloadHandler(
      filename = function() "stability_selection.csv",
      content = function(file) utils::write.csv(fs_stab_result()$ranked, file, row.names = FALSE)
    )

    ## ---------------------------------------------------------------------
    ## Stage 7: Consensus / Overlap (+ correlation reduction)
    ## ---------------------------------------------------------------------

    fs_available_methods <- reactive({
      avail <- c(Univariate = uni_has_run(), Regularization = reg_has_run(), TreeBased = rf_has_run(), RFE = rfe_has_run(), Stability = stab_has_run())
      names(avail)[avail]
    })

    fs_method_ids <- reactive({
      out <- list()
      if (uni_has_run()) out$Univariate <- fs_uni_result()$selected_ids
      if (reg_has_run()) out$Regularization <- fs_reg_result()$selected_ids
      if (rf_has_run()) out$TreeBased <- fs_rf_result()$selected_ids
      if (rfe_has_run()) out$RFE <- fs_rfe_result()$selected_ids
      if (stab_has_run()) out$Stability <- fs_stab_result()$selected_ids
      out
    })

    output$consensus_ui <- renderUI({
      if (!filters_has_run()) return(mod_methyl_fs_need_filters_note())
      avail <- fs_available_methods()
      if (length(avail) == 0) {
        return(mod_methyl_fs_empty_note("Run at least one method (Univariate/Regularization/Tree-Based/RFE/Stability) before computing consensus."))
      }
      div(class = "card",
        div(class = "card-title", icon("object-group"), "Consensus / Overlap"),
        fluidRow(
          column(4,
            checkboxGroupInput(ns("consensus_methods"), "Methods to include", choices = avail, selected = avail),
            numericInput(ns("consensus_min_methods"), "Minimum number of methods", value = min(2, length(avail)), min = 1, max = length(avail), step = 1)
          ),
          column(4, uiOutput(ns("consensus_weight_ui"))),
          column(4,
            checkboxInput(ns("consensus_use_weighted"), "Also require a minimum weighted score", value = FALSE),
            conditionalPanel(condition = sprintf("input['%s']", ns("consensus_use_weighted")), numericInput(ns("consensus_min_weighted"), "Minimum weighted score", value = 0.5, min = 0, max = 1, step = 0.05))
          )
        ),
        actionButton(ns("consensus_run_btn"), "Run Consensus Selection", icon = icon("play"), class = "btn-primary"),
        uiOutput(ns("consensus_result_ui"))
      )
    })

    output$consensus_weight_ui <- renderUI({
      avail <- fs_available_methods()
      req(length(avail) > 0)
      tagList(lapply(avail, function(m) numericInput(ns(paste0("consensus_w_", m)), sprintf("Weight: %s", m), value = 1, min = 0, max = 10, step = 0.5)))
    })

    consensus_has_run <- reactiveVal(FALSE)
    observeEvent(fs_filter_result(), consensus_has_run(FALSE))
    observeEvent(input$consensus_run_btn, consensus_has_run(TRUE), ignoreInit = TRUE)

    fs_consensus_result <- eventReactive(input$consensus_run_btn, {
      ids <- fs_method_ids()
      chosen <- intersect(input$consensus_methods %||% names(ids), names(ids))
      validate(need(length(chosen) >= 1, "Pick at least one method."))
      ids <- ids[chosen]
      weights <- stats::setNames(vapply(chosen, function(m) input[[paste0("consensus_w_", m)]] %||% 1, numeric(1)), chosen)
      ct <- methyl_fs_consensus_table(ids, weights = weights)
      sel <- methyl_fs_consensus_select(ct$table, min_methods = input$consensus_min_methods,
                                         min_weighted = if (isTRUE(input$consensus_use_weighted)) input$consensus_min_weighted else NULL)
      list(table = ct$table, methods = chosen, weights = weights, selected_ids = sel,
           min_methods = input$consensus_min_methods, run_at = Sys.time())
    })

    observeEvent(fs_consensus_result(), {
      r <- fs_consensus_result()
      if (!is.null(results)) {
        results$featureselection <- list(
          selected_cpgs = r$selected_ids, n_selected = length(r$selected_ids),
          method_counts = vapply(r$methods, function(m) sum(r$table[[m]] == 1), integer(1)),
          consensus_rule = sprintf(">= %d of %d methods", r$min_methods, length(r$methods))
        )
      }
    })

    output$consensus_result_ui <- renderUI({
      if (!consensus_has_run()) return(mod_methyl_fs_empty_note())
      r <- fs_consensus_result()
      tagList(
        div(class = "card",
          div(class = "card-title", icon("circle-check"), "Consensus result"),
          p(strong(length(r$selected_ids)), sprintf(" CpGs meet the consensus threshold (>= %d of %d methods), across %s.",
                                                       r$min_methods, length(r$methods), paste(r$methods, collapse = ", "))),
          fluidRow(
            column(6, tags$h5("Venn diagram"), withSpinner(plotOutput(ns("consensus_venn"), height = 320), color = "#2563EB", type = 6)),
            column(6, tags$h5("UpSet-style plot"), withSpinner(plotOutput(ns("consensus_upset"), height = 320), color = "#2563EB", type = 6))
          ),
          fluidRow(
            column(6, tags$h5("Consensus ranking"), withSpinner(plotOutput(ns("consensus_rank_plot"), height = 320), color = "#2563EB", type = 6)),
            column(6, tags$h5("Method comparison heatmap"), withSpinner(plotOutput(ns("consensus_heatmap"), height = 320), color = "#2563EB", type = 6))
          ),
          div(class = "table-toolbar", downloadButton(ns("consensus_download"), "Intersection table (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns("consensus_table_out"))
        ),
        div(class = "card",
          div(class = "card-title", icon("link-slash"), "Correlation-based CpG reduction (optional)"),
          p(class = "submodule-desc", "Reduces the consensus panel for collinearity - only applied on click, and only to the (small) consensus panel above."),
          fluidRow(
            column(4, radioButtons(ns("corr_method"), "Correlation method", choices = c("Pearson" = "pearson", "Spearman" = "spearman"), selected = "pearson")),
            column(4, numericInput(ns("corr_threshold"), "Correlation threshold", value = 0.8, min = 0.5, max = 0.99, step = 0.05)),
            column(4, checkboxInput(ns("use_corr_reduced"), "Use the correlation-reduced panel downstream", value = FALSE))
          ),
          actionButton(ns("corr_reduce_btn"), "Reduce for collinearity", icon = icon("compress"), class = "btn-default btn-sm"),
          shinyjs::hidden(div(id = ns("corr_wrap"), uiOutput(ns("corr_result_ui"))))
        )
      )
    })

    output$consensus_venn <- renderPlot({
      req(consensus_has_run())
      r <- fs_consensus_result()
      ids <- fs_method_ids()[r$methods]
      draw_overlap_venn(ids, title = sprintf("%d-CpG consensus", length(r$selected_ids)))
    })

    output$consensus_upset <- renderPlot({
      req(consensus_has_run())
      methyl_fs_upset_plot(fs_method_ids()[fs_consensus_result()$methods])
    })

    output$consensus_rank_plot <- renderPlot({ req(consensus_has_run()); methyl_fs_consensus_rank_plot(fs_consensus_result()$table) })
    output$consensus_heatmap <- renderPlot({ req(consensus_has_run()); methyl_fs_method_heatmap(fs_consensus_result()$table, fs_consensus_result()$methods) })

    output$consensus_table_out <- DT::renderDataTable({
      req(consensus_has_run())
      DT::datatable(fs_consensus_result()$table, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatSignif(columns = "weighted_score", digits = 4)
    })
    outputOptions(output, "consensus_table_out", suspendWhenHidden = FALSE)

    output$consensus_download <- downloadHandler(
      filename = function() "consensus_intersection_table.csv",
      content = function(file) utils::write.csv(fs_consensus_result()$table, file, row.names = FALSE)
    )

    corr_has_run <- reactiveVal(FALSE)
    fs_corr_result <- eventReactive(input$corr_reduce_btn, {
      r <- fs_filter_result(); cr <- fs_consensus_result()
      score <- stats::setNames(cr$table$weighted_score, cr$table$cpg)
      methyl_fs_correlation_reduce(r$beta, cr$selected_ids, score, method = input$corr_method, r_threshold = input$corr_threshold)
    })
    observeEvent(input$corr_reduce_btn, { corr_has_run(TRUE); shinyjs::show(id = "corr_wrap") }, ignoreInit = TRUE)

    output$corr_result_ui <- renderUI({
      req(corr_has_run())
      r <- fs_corr_result()
      cr <- fs_consensus_result()
      tagList(
        p(sprintf("%d of %d CpGs removed for collinearity (|r| > %.2f); %d remain.",
                   length(cr$selected_ids) - length(r$reduced_ids), length(cr$selected_ids), input$corr_threshold, length(r$reduced_ids))),
        DT::dataTableOutput(ns("corr_dropped_table"))
      )
    })

    output$corr_dropped_table <- DT::renderDataTable({
      req(corr_has_run())
      DT::datatable(fs_corr_result()$dropped, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 5), class = "stripe hover compact")
    })
    outputOptions(output, "corr_dropped_table", suspendWhenHidden = FALSE)

    ## ---------------------------------------------------------------------
    ## Stage 8: Selected Features
    ## ---------------------------------------------------------------------

    fs_final_panel <- reactive({
      req(consensus_has_run())
      cr <- fs_consensus_result()
      ids <- if (corr_has_run() && isTRUE(input$use_corr_reduced)) fs_corr_result()$reduced_ids else cr$selected_ids
      list(ids = ids, table = cr$table[cr$table$cpg %in% ids, , drop = FALSE])
    })

    output$selected_ui <- renderUI({
      if (!consensus_has_run()) return(mod_methyl_fs_empty_note("Run Consensus Selection first."))
      panel <- fs_final_panel()
      r <- fs_filter_result()
      div(class = "card",
        div(class = "card-title", icon("list-check"), "Selected Features"),
        p(strong(length(panel$ids)), " final CpGs."),
        div(class = "table-toolbar",
            downloadButton(ns("selected_download_csv"), "CSV", class = "btn-sm"),
            downloadButton(ns("selected_download_tsv"), "TSV", class = "btn-sm"),
            downloadButton(ns("selected_copy"), "Copy IDs (TXT)", class = "btn-sm")),
        DT::dataTableOutput(ns("selected_table")),
        shinyjs::hidden(div(id = ns("cpg_detail_wrap"), uiOutput(ns("cpg_detail_ui"))))
      )
    })

    fs_selected_table_full <- reactive({
      panel <- fs_final_panel()
      r <- fs_filter_result()
      anno_df <- methyl_fs_annotate_panel(panel$ids, r$anno)
      merge(panel$table, anno_df, by = "cpg", all.x = TRUE, sort = FALSE)[order(-panel$table$weighted_score[match(panel$table$cpg, panel$table$cpg)]), , drop = FALSE]
    })

    output$selected_table <- DT::renderDataTable({
      req(consensus_has_run())
      DT::datatable(fs_selected_table_full(), rownames = FALSE, filter = "top", selection = "single",
                    options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact")
    })
    outputOptions(output, "selected_table", suspendWhenHidden = FALSE)

    observeEvent(input$selected_table_rows_selected, {
      req(length(input$selected_table_rows_selected) > 0)
      shinyjs::show(id = "cpg_detail_wrap")
    })

    output$cpg_detail_ui <- renderUI({
      req(length(input$selected_table_rows_selected) > 0)
      df <- fs_selected_table_full()
      row <- df[input$selected_table_rows_selected, , drop = FALSE]
      cpg <- row$cpg[1]
      r <- fs_filter_result()
      tagList(
        div(class = "card-title", icon("magnifying-glass"), sprintf("Detail: %s", cpg)),
        p(strong("Methods selecting it: "), paste(names(fs_method_ids())[vapply(fs_method_ids(), function(s) cpg %in% s, logical(1))], collapse = ", ")),
        p(strong("Consensus score: "), sprintf("%.3f (%d methods)", row$weighted_score[1], row$n_methods[1])),
        withSpinner(plotOutput(ns("cpg_detail_plot"), height = 260), color = "#2563EB", type = 6),
        DT::dataTableOutput(ns("cpg_detail_table"))
      )
    })

    output$cpg_detail_plot <- renderPlot({
      req(length(input$selected_table_rows_selected) > 0)
      cpg <- fs_selected_table_full()$cpg[input$selected_table_rows_selected]
      r <- fs_filter_result()
      req(cpg %in% rownames(r$beta))
      df <- data.frame(sample = colnames(r$beta), beta = r$beta[cpg, ], group = as.character(r$grp))
      ggplot(df, aes(x = group, y = beta, fill = group)) +
        geom_violin(alpha = 0.5) + geom_boxplot(width = 0.15, outlier.size = 0.6) +
        scale_fill_manual(values = arthomix_pair(df$group)) +
        labs(x = NULL, y = "Beta value", fill = NULL) + theme_arthomix()
    })

    output$cpg_detail_table <- DT::renderDataTable({
      req(length(input$selected_table_rows_selected) > 0)
      cpg <- fs_selected_table_full()$cpg[input$selected_table_rows_selected]
      r <- fs_filter_result()
      req(cpg %in% rownames(r$beta))
      DT::datatable(data.frame(sample = colnames(r$beta), group = as.character(r$grp), beta = r$beta[cpg, ], row.names = NULL),
                    rownames = FALSE, options = list(scrollX = TRUE, pageLength = 5), class = "stripe hover compact")
    })
    outputOptions(output, "cpg_detail_table", suspendWhenHidden = FALSE)

    output$selected_download_csv <- downloadHandler(
      filename = function() "selected_features.csv",
      content = function(file) utils::write.csv(fs_selected_table_full(), file, row.names = FALSE)
    )
    output$selected_download_tsv <- downloadHandler(
      filename = function() "selected_features.tsv",
      content = function(file) utils::write.table(fs_selected_table_full(), file, sep = "\t", row.names = FALSE, quote = FALSE)
    )
    output$selected_copy <- downloadHandler(
      filename = function() "selected_cpg_ids.txt",
      content = function(file) writeLines(fs_final_panel()$ids, file)
    )

    ## ---------------------------------------------------------------------
    ## Stage 9: Model & Export (validation + save/load)
    ## ---------------------------------------------------------------------

    output$export_ui <- renderUI({
      if (!consensus_has_run()) return(mod_methyl_fs_empty_note("Run Consensus Selection first."))
      tagList(
        div(class = "card",
          div(class = "card-title", icon("vial"), "Validate selected panel"),
          fluidRow(
            column(4, radioButtons(ns("validate_mode"), "Validation mode",
                                    choices = c("Frozen panel (fast, optimistic)" = "frozen", "Leakage-safe (reselect per fold)" = "nested"), selected = "frozen")),
            column(4, radioButtons(ns("validate_classifier"), "Classifier", choices = c("Logistic Regression" = "glm", "Random Forest" = "rf", "SVM" = "svm"), selected = "glm")),
            column(4, numericInput(ns("validate_k"), "Folds (k)", value = 5, min = 3, max = 10, step = 1), numericInput(ns("validate_repeats"), "Repeats", value = 1, min = 1, max = 10, step = 1))
          ),
          conditionalPanel(condition = sprintf("input['%s'] == 'frozen'", ns("validate_mode")),
            p(class = "empty-note", icon("triangle-exclamation"), "This mode evaluates the panel that was already chosen using this same data, so the reported AUC is optimistic - use \"Leakage-safe\" for an honest estimate.")),
          conditionalPanel(condition = sprintf("input['%s'] == 'nested'", ns("validate_mode")),
            p(class = "empty-note", icon("circle-info"), "Leakage-safe mode reselects the panel using the Univariate + Regularization steps inside every outer fold; it does not rerun Tree-Based/RFE/Stability/Consensus per fold, since a full 5-method ensemble refit per fold is impractical in a live session.")),
          actionButton(ns("validate_run_btn"), "Run Validation", icon = icon("play"), class = "btn-primary"),
          uiOutput(ns("validate_result_ui"))
        ),
        div(class = "card",
          div(class = "card-title", icon("floppy-disk"), "Save Model as RDS"),
          p(class = "submodule-desc", "Bundles the full configuration - filters, every method's parameters/results, consensus config, validation metrics, fitted model, and session info - not just the CpG list."),
          downloadButton(ns("save_rds"), "Save Model as RDS", class = "btn-primary btn-sm"),
          div(class = "table-toolbar",
              downloadButton(ns("export_ranking"), "Feature-ranking table (CSV)", class = "btn-sm"),
              downloadButton(ns("export_consensus"), "Consensus table (CSV)", class = "btn-sm"),
              downloadButton(ns("export_summary"), "Analysis summary (TXT)", class = "btn-sm"))
        ),
        div(class = "card",
          div(class = "card-title", icon("upload"), "Load Previous RDS Model"),
          fileInput(ns("load_model_file"), "Load a previously saved model (.rds)"),
          shinyjs::hidden(div(id = ns("loaded_model_wrap"), uiOutput(ns("loaded_model_ui"))))
        )
      )
    })

    validate_has_run <- reactiveVal(FALSE)
    observeEvent(fs_final_panel(), validate_has_run(FALSE))
    observeEvent(input$validate_run_btn, validate_has_run(TRUE), ignoreInit = TRUE)

    fs_validate_result <- eventReactive(input$validate_run_btn, {
      r <- fs_filter_result(); panel <- fs_final_panel()
      validate(need(length(panel$ids) >= 2, "Need at least 2 selected CpGs to validate."))
      if (identical(input$validate_mode, "nested")) {
        res <- methyl_fs_validate_nested(r$beta, r$m, r$grp, uni_params = list(rule = "top_n", top_n = 100),
                                          lasso_params = list(alpha = 1), classifier = input$validate_classifier,
                                          outer_k = input$validate_k, repeats = input$validate_repeats)
        list(mode = "nested", mean_auc = res$mean_auc, sd_auc = res$sd_auc, per_fold = res$per_fold, run_at = Sys.time())
      } else {
        X <- t(r$m[panel$ids, , drop = FALSE])
        res <- methyl_fs_validate_frozen(X, r$grp, classifier = input$validate_classifier, k = input$validate_k, repeats = input$validate_repeats)
        list(mode = "frozen", mean_auc = res$mean_auc, per_fold = res$resample_results, run_at = Sys.time())
      }
    })

    output$validate_result_ui <- renderUI({
      if (!validate_has_run()) return(mod_methyl_fs_empty_note())
      r <- fs_validate_result()
      tagList(
        p(strong("Mean AUC: "), sprintf("%.3f", r$mean_auc), if (!is.null(r$sd_auc)) sprintf(" (SD %.3f)", r$sd_auc) else NULL,
          sprintf(" - %s mode.", if (identical(r$mode, "nested")) "leakage-safe" else "frozen-panel")),
        DT::dataTableOutput(ns("validate_table"))
      )
    })

    output$validate_table <- DT::renderDataTable({
      req(validate_has_run())
      DT::datatable(fs_validate_result()$per_fold, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact")
    })
    outputOptions(output, "validate_table", suspendWhenHidden = FALSE)

    fs_model_export <- reactive({
      r <- fs_filter_result(); panel <- fs_final_panel(); cr <- fs_consensus_result()
      list(
        schema_version = "1.0", module = "mod_methyl_featureselection", saved_at = Sys.time(),
        dataset = list(source_label = r$source_label, array_type = r$array_type, n_probes_input = r$n_probes,
                        n_samples = r$n_samples, group_col = r$group_col, reference_level = r$ref_group, comparison_level = r$comp_group),
        filters = list(notes = r$filter_notes, cascade = r$cascade, imputed_method = r$imputed_method),
        methods = list(
          univariate = if (uni_has_run()) list(ran = TRUE, selected_ids = fs_uni_result()$selected_ids, ranked = fs_uni_result()$ranked) else list(ran = FALSE),
          regularization = if (reg_has_run()) list(ran = TRUE, selected_ids = fs_reg_result()$selected_ids, coefficients = fs_reg_result()$coefficients) else list(ran = FALSE),
          tree_based = if (rf_has_run()) list(ran = TRUE, selected_ids = fs_rf_result()$selected_ids, ranked = fs_rf_result()$ranked) else list(ran = FALSE),
          rfe = if (rfe_has_run()) list(ran = TRUE, selected_ids = fs_rfe_result()$selected_ids, optimal_size = fs_rfe_result()$optimal_size) else list(ran = FALSE),
          stability = if (stab_has_run()) list(ran = TRUE, selected_ids = fs_stab_result()$selected_ids, ranked = fs_stab_result()$ranked) else list(ran = FALSE)
        ),
        consensus = list(methods_used = cr$methods, min_methods_required = cr$min_methods, weights = cr$weights, selected_ids = cr$selected_ids, table = cr$table),
        correlation_reduction = list(applied = isTRUE(input$use_corr_reduced) && corr_has_run(),
                                       r_threshold = if (corr_has_run()) input$corr_threshold else NA, reduced_ids = if (corr_has_run()) fs_corr_result()$reduced_ids else NULL),
        final_panel = list(cpg_ids = panel$ids, n = length(panel$ids), annotation = fs_selected_table_full()),
        validation = if (validate_has_run()) fs_validate_result() else list(mode = NA),
        session_info = list(r_version = getRversion(), package_versions = list(
          glmnet = as.character(utils::packageVersion("glmnet")), randomForest = as.character(utils::packageVersion("randomForest")),
          caret = as.character(utils::packageVersion("caret")), e1071 = as.character(utils::packageVersion("e1071"))))
      )
    })

    output$save_rds <- downloadHandler(
      filename = function() sprintf("methyl_featureselection_%s.rds", format(Sys.Date(), "%Y%m%d")),
      content = function(file) saveRDS(fs_model_export(), file)
    )
    output$export_ranking <- downloadHandler(
      filename = function() "feature_ranking.csv",
      content = function(file) utils::write.csv(fs_selected_table_full(), file, row.names = FALSE)
    )
    output$export_consensus <- downloadHandler(
      filename = function() "consensus_table.csv",
      content = function(file) utils::write.csv(fs_consensus_result()$table, file, row.names = FALSE)
    )
    output$export_summary <- downloadHandler(
      filename = function() "analysis_summary.txt",
      content = function(file) {
        r <- fs_filter_result(); panel <- fs_final_panel()
        lines <- c(
          "Methylomics ML Feature Selection - analysis summary",
          sprintf("Dataset: %s", r$source_label),
          sprintf("CpGs after filtering: %d, Samples: %d", r$n_probes, r$n_samples),
          sprintf("Comparison: %s vs %s (column: %s)", r$comp_group, r$ref_group, r$group_col),
          sprintf("Final selected CpGs: %d", length(panel$ids)),
          if (validate_has_run()) sprintf("Validation mean AUC: %.3f (%s mode)", fs_validate_result()$mean_auc, fs_validate_result()$mode) else "Validation: not run",
          sprintf("Generated: %s", format(Sys.time()))
        )
        writeLines(lines, file)
      }
    )

    fs_loaded_model <- eventReactive(input$load_model_file, {
      x <- tryCatch(readRDS(input$load_model_file$datapath), error = function(e) NULL)
      validate(need(!is.null(x) && identical(x$module, "mod_methyl_featureselection"),
                    "This doesn't look like a saved ML Feature Selection model from this app."))
      x
    }, ignoreInit = TRUE)

    observeEvent(input$load_model_file, shinyjs::show(id = "loaded_model_wrap"), ignoreInit = TRUE)

    output$loaded_model_ui <- renderUI({
      x <- tryCatch(fs_loaded_model(), error = function(e) NULL)
      req(x)
      tagList(
        p(class = "empty-note", icon("circle-info"), "Loaded for reference/comparison only - does not modify the current run."),
        p(strong("Saved at: "), format(x$saved_at)),
        p(strong("Source: "), x$dataset$source_label, sprintf(" (%d CpGs, %d samples)", x$dataset$n_probes_input, x$dataset$n_samples)),
        p(strong("Final panel: "), sprintf("%d CpGs", x$final_panel$n)),
        if (!identical(x$validation$mode, NA) && !is.null(x$validation$mode)) p(strong("Validation: "), sprintf("mean AUC %.3f (%s)", x$validation$mean_auc, x$validation$mode)) else NULL,
        DT::dataTableOutput(ns("loaded_model_table"))
      )
    })

    output$loaded_model_table <- DT::renderDataTable({
      x <- tryCatch(fs_loaded_model(), error = function(e) NULL)
      req(x)
      DT::datatable(x$final_panel$annotation, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact")
    })
    outputOptions(output, "loaded_model_table", suspendWhenHidden = FALSE)
  })
}
