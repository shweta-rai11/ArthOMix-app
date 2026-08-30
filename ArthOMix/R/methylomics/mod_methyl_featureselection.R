## R/methylomics/mod_methyl_featureselection.R
## CpG feature-selection workflow: filter probes, run Univariate/LASSO/RF/RFE/
## Stability selection independently, combine into a weighted Consensus panel
## (with optional correlation-based reduction), then validate and export.
##
## Reuses filter/parsing/annotation/plotting helpers from qc.R / parse_upload.R
## / annotation.R / global.R without modifying them.

mod_methyl_featureselection_config <- list(
  id = "featureselection", title = "ML Feature Selection", icon = "sliders", group = "Biomarker modeling",
  description = "CpG feature selection using univariate, regularization, tree-based, RFE, and stability selection methods"
)

## =============================================================================
## Helpers - data source / scale / grouping
## =============================================================================

## .rds companion to qc.R's methyl_parse_matrix() (CSV/TSV only); same list(ok, mat, error) shape.
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

## Beta-vs-M-value heuristic: checks whether the 0.1-99.9% quantile of a sample
## sits inside [-0.05, 1.05]. Only used for uploads; preloaded input_scale is authoritative.
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

## IQR companion to qc.R's methyl_filter_variance()/methyl_filter_sd(); same list(keep, note) shape.
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

## Per-probe median imputation by default; uses k-NN (impute::impute.knn) when that
## optional package is installed, falling back to median otherwise.
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

## Linear-model engine (t-test / moderated t-test / ANOVA / linear regression / Pearson):
## different design/contrast choices on the same limma::lmFit()/eBayes() fit on M-values,
## optionally with covariate adjustment.
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

## Rank-based engine (Wilcoxon / Kruskal-Wallis / Spearman); no covariate adjustment.
## Base R apply() is fine here since the panel is already capped by Data & Filters.
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
  two_group <- FALSE
  if (!is_continuous) {
    grp <- factor(as.character(y))
    lv <- levels(grp)
    two_group <- length(lv) == 2
    if (two_group) {
      mean1 <- rowMeans(beta_mat[df$cpg, grp == lv[1], drop = FALSE], na.rm = TRUE)
      mean2 <- rowMeans(beta_mat[df$cpg, grp == lv[2], drop = FALSE], na.rm = TRUE)
      df$mean_beta_ref <- mean1; df$mean_beta_comp <- mean2
      df$dbeta <- mean2 - mean1
      df$direction <- ifelse(df$dbeta > 0, "hyper", "hypo")
    } else {
      ## >2 groups: reference-vs-comparison delta-beta isn't well-defined, leave NA.
      df$mean_beta_ref <- NA_real_; df$mean_beta_comp <- NA_real_
      df$dbeta <- NA_real_; df$direction <- NA_character_
    }
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
  if (!is.na(dbeta_min) && dbeta_min > 0 && two_group) keep <- keep & !is.na(df$dbeta) & abs(df$dbeta) >= dbeta_min
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
  ## randomForest() errors if sampsize/maxnodes/classwt are passed as explicit NULL,
  ## so unset optional args are only added to the call when actually provided.
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

## Ported from the transcriptomics module's SVM-RFE. Backward elimination by
## smallest squared linear-SVM weight, then a k-fold CV error curve to pick panel size.
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

## Generates resample index sets (bootstrap / repeated k-fold / subsampling) for
## Meinshausen-Buhlmann stability selection against a fixed reference-lambda LASSO.
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
  ## Tracks which resamples produced a fit, so failed/degenerate ones are excluded
  ## from the frequency denominator instead of deflating selection_frequency.
  valid <- rep(FALSE, length(idx_sets))
  for (i in seq_along(idx_sets)) {
    idx <- idx_sets[[i]]
    if (length(unique(as.character(y[idx]))) < 2) next
    fit <- tryCatch(glmnet::glmnet(X[idx, , drop = FALSE], y[idx], family = "binomial", alpha = alpha, lambda = ref_lambda),
                     error = function(e) NULL)
    if (is.null(fit)) next
    co <- stats::coef(fit)[-1, 1]
    sel_mat[names(co)[co != 0], i] <- 1L
    valid[i] <- TRUE
  }
  validate(need(any(valid), "Every resample failed to fit (e.g. too few samples per class in a resample) - reduce the sample fraction or number of resamples, or check class balance."))
  sel_mat <- sel_mat[, valid, drop = FALSE]
  freq <- rowMeans(sel_mat)
  mean_rank <- apply(sel_mat, 1, function(r) if (sum(r) == 0) NA_real_ else mean(which(r == 1)))
  stability_score <- freq
  ranked <- data.frame(cpg = names(freq), selection_frequency = as.numeric(freq),
                        stability_score = as.numeric(stability_score), row.names = NULL)
  ranked <- ranked[order(-ranked$selection_frequency), , drop = FALSE]
  list(n_resamples = sum(valid), n_resamples_requested = length(idx_sets), ref_lambda = ref_lambda, ranked = ranked,
       selected_ids = ranked$cpg[ranked$selection_frequency >= freq_threshold])
}

## =============================================================================
## Helpers - Consensus / Overlap / Correlation reduction
## =============================================================================

methyl_fs_consensus_table <- function(id_lists, weights = NULL) {
  all_ids <- Reduce(union, id_lists)
  methods <- names(id_lists)
  w <- weights %||% stats::setNames(rep(1, length(methods)), methods)
  w <- w[methods]
  if (length(all_ids) == 0) {
    ## Every method selected zero CpGs (plausible on a small/noisy panel) -
    ## still return the full expected column set (one 0/1 column per method,
    ## n_methods, weighted_score), just with zero rows, so every downstream
    ## consumer (Venn/UpSet/rank plot/heatmap/table) can keep assuming those
    ## columns exist instead of special-casing an empty consensus.
    empty <- stats::setNames(data.frame(cpg = character(0), stringsAsFactors = FALSE),
                              "cpg")
    for (m in methods) empty[[m]] <- integer(0)
    empty$n_methods <- integer(0)
    empty$weighted_score <- numeric(0)
    return(list(table = empty, all_ids = character(0), methods = methods, weights = w))
  }
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

## Greedy pairwise-correlation pruning: repeatedly drops the lower-scoring member
## of the highest-|r| pair above r_threshold, until no pair exceeds it.
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

## Same chr/pos/gene join pattern as mod_methyl_dmp.R:648-655. CpG-island/genomic-region
## columns left NA - shared annotation.R only exposes chr/pos/Type/SNP-rs/gene.
methyl_fs_annotate_panel <- function(ids, anno_result) {
  ## data.frame() can't recycle a length-1 NA against a length-0 `ids` (a
  ## legitimate empty panel - e.g. every method selected zero CpGs) without
  ## an explicit rep(), which throws "differing number of rows: 0, 1" instead.
  n <- length(ids)
  df <- data.frame(cpg = ids, chr = rep(NA_character_, n), pos = rep(NA_real_, n), gene = rep(NA_character_, n),
                    cpg_island = rep(NA_character_, n), genomic_region = rep(NA_character_, n), stringsAsFactors = FALSE)
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

## Leakage-safe: outer k-fold loop reselecting the panel (Univariate + Regularization
## only) inside each training fold, mirroring the diagnostic module's CV shape.
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

## Sex-stratification: female/male fit separately, pooled always available -
## mirrors transcriptomics' three-build parallel-panel pattern
## (R/transcriptomics/mod_featureselection.R), adapted to this file's own
## sex-detection helpers (mod_methyl_dmp_sex_col/_choices) and its existing
## one-Run-button-per-stage pipeline shape.
FS_SEXES <- c("female", "male", "pooled")

## Arranges up to 3 per-sex result cards: Female + Male side-by-side when both
## are present, Pooled always full-width below - the methylomics-card-styled
## analog of transcriptomics' fs-pair-row layout.
mod_methyl_fs_sex_layout <- function(cards) {
  fm <- cards[intersect(c("female", "male"), names(cards))]
  pooled <- cards[["pooled"]]
  tagList(
    if (length(fm) == 2) fluidRow(column(6, fm[[1]]), column(6, fm[[2]]))
    else if (length(fm) == 1) fluidRow(column(12, fm[[1]])),
    if (!is.null(pooled)) fluidRow(column(12, pooled))
  )
}

## =============================================================================
## Server
## =============================================================================

mod_methyl_featureselection_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## ---------------------------------------------------------------------
    ## Stage 0: source resolution - always live, never Run-gated; heavy computation
    ## stays behind explicit Run buttons.
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
        validate(need(!is.null(dataset$beta), "No dataset is loaded yet - open the Methylomics Dataset tab and load one (preloaded, upload, or GEO fetch), or switch to \"Upload my own\" above."))
        v <- dataset$beta[is.finite(dataset$beta)]
        frac_out <- if (length(v) > 0) mean((if (length(v) > 20000) v[sample.int(length(v), 20000)] else v) < -0.05 | (if (length(v) > 20000) v[sample.int(length(v), 20000)] else v) > 1.05) else 0
        list(mat = dataset$beta, sheet = dataset$sample_sheet, scale = dataset$input_scale %||% "beta",
             frac_out_of_range = frac_out, array_type = dataset$array_type %||% "450K",
             source_label = dataset$source %||% "Dataset tab", preloaded = isTRUE(dataset$preloaded))
      }
    })

    anno_result <- reactive({ methyl_get_annotation(fs_active_source()$array_type) })

    fs_group_choices <- reactive({
      src <- fs_active_source()
      if (is.null(src$sheet)) return(character(0))
      colnames(src$sheet)
    })

    ## Sex detection reuses this app's own DMP/DMR/WGCNA idiom
    ## (mod_methyl_dmp_sex_col()/mod_methyl_dmp_sex_choices(), mod_methyl_dmp.R:102-123)
    ## rather than inventing a bespoke sex-guessing reactive like transcriptomics does.
    fs_sex_col <- reactive({ mod_methyl_dmp_sex_col(fs_active_source()$sheet) })
    fs_sex_choices <- reactive({ mod_methyl_dmp_sex_choices(fs_active_source()$sheet, fs_sex_col()) })
    fs_sex_available <- reactive({ length(fs_sex_choices()) > 1 })

    ## "female"/"male" -> the raw sex-column value to filter on; "pooled" -> NULL
    ## (skips the sex filter entirely, pooling every sample regardless of sex).
    fs_sex_value_for <- function(sex_label) {
      if (identical(sex_label, "pooled")) return(NULL)
      ch <- fs_sex_choices()
      rest <- ch[-1]
      if (identical(sex_label, "female") && length(rest) >= 1) return(unname(rest[1]))
      if (identical(sex_label, "male")   && length(rest) >= 2) return(unname(rest[2]))
      NA_character_
    }

    fs_sex_display_label <- function(sex_label) {
      if (identical(sex_label, "pooled")) return("Pooled (all samples)")
      ch <- fs_sex_choices()
      nm <- names(ch)[-1]
      if (identical(sex_label, "female") && length(nm) >= 1) return(nm[1])
      if (identical(sex_label, "male")   && length(nm) >= 2) return(nm[2])
      tools::toTitleCase(sex_label)
    }

    ## ---------------------------------------------------------------------
    ## Stage 1: Data & Filters
    ## ---------------------------------------------------------------------

    output$filters_ui <- renderUI({
      tagList(
        div(class = "card",
          div(class = "card-title", icon("filter"), "Data source"),
          ## Gated dynamically on dataset$beta rather than the static
          ## METH_RAW_DATA_AVAILABLE flag - this choice reads whatever is
          ## actually loaded on the Methylomics Dataset tab (preloaded,
          ## uploaded, or GEO-fetched alike), so it should only appear/default
          ## once something has actually been loaded there.
          radioButtons(ns("fs_source"), NULL,
                       choices = if (!is.null(dataset$beta)) c("Use dataset loaded on the Dataset tab" = "preloaded", "Upload my own" = "upload") else c("Upload my own" = "upload"),
                       selected = if (!is.null(dataset$beta)) "preloaded" else "upload", inline = TRUE),
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
        tags$h5(icon("venus-mars"), " Sex"),
        if (!fs_sex_available())
          p(class = "empty-note", icon("circle-info"), "No usable sex information was found for this dataset - only \"Run All (pooled)\" is available.")
        else
          p(class = "empty-note", icon("circle-info"),
            sprintf("%s and %s fit separately below; \"Run All (pooled)\" pools every sample regardless of sex - sex is never used as a covariate here.",
                    fs_sex_display_label("female"), fs_sex_display_label("male"))),
        div(style = "display:flex; gap:8px; flex-wrap: wrap;",
            if (fs_sex_available()) actionButton(ns("run_female_btn"), sprintf("Run %s", fs_sex_display_label("female")), icon = icon("play"), class = "btn-primary btn-sm"),
            if (fs_sex_available()) actionButton(ns("run_male_btn"), sprintf("Run %s", fs_sex_display_label("male")), icon = icon("play"), class = "btn-primary btn-sm"),
            actionButton(ns("run_pooled_btn"), "Run All (pooled)", icon = icon("play"), class = "btn-primary btn-sm")
        )
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

    ## Staleness: every sex's build goes stale when the source dataset changes,
    ## cleared the moment that sex's own button is clicked again - mirrors
    ## transcriptomics' fs_stale (R/transcriptomics/mod_featureselection.R:868-871).
    fs_stale <- reactiveValues(female = FALSE, male = FALSE, pooled = FALSE)
    observeEvent(dataset$beta, { fs_stale$female <- TRUE; fs_stale$male <- TRUE; fs_stale$pooled <- TRUE }, ignoreNULL = TRUE)
    observeEvent(input$fs_source, { fs_stale$female <- TRUE; fs_stale$male <- TRUE; fs_stale$pooled <- TRUE })
    observeEvent(input$run_female_btn, fs_stale$female <- FALSE, ignoreInit = TRUE)
    observeEvent(input$run_male_btn,   fs_stale$male   <- FALSE, ignoreInit = TRUE)
    observeEvent(input$run_pooled_btn, fs_stale$pooled <- FALSE, ignoreInit = TRUE)

    ## The one filter pipeline (missingness/variance/methylation-specific/
    ## cross-reactive/MAF/imputation) - unchanged from before sex stratification,
    ## just parameterized by which sex-value to restrict the sample set to first
    ## (sex_value = NULL, the "pooled" case, skips that restriction entirely).
    fs_build_filters <- function(sex_label, sex_value) {
      src <- fs_active_source()
      sheet <- src$sheet
      validate(need(!is.null(sheet) && !is.null(input$fs_group_col), "Load phenotype/sample metadata and pick a group column first."))
      ## Guards factor(..., levels = c(ref, comp)) below against duplicated levels;
      ## nothing in the UI prevents both selectInputs pointing at the same value.
      validate(need(!is.null(input$fs_ref_group) && !is.null(input$fs_comp_group) && !identical(input$fs_ref_group, input$fs_comp_group),
                    "Reference and comparison group must be two different levels."))
      sample_ids <- methyl_sheet_sample_ids(sheet, colnames(src$mat))
      common <- intersect(colnames(src$mat), sample_ids)
      validate(need(length(common) >= 10, "Fewer than 10 samples match between the matrix and the phenotype file."))
      mat0 <- src$mat[, common, drop = FALSE]
      ph0 <- as.data.frame(sheet)[match(common, sample_ids), , drop = FALSE]

      if (!is.null(sex_value)) {
        validate(need(!is.na(sex_value), sprintf("No usable %s stratum is available in this dataset.", sex_label)))
        sc <- fs_sex_col()
        validate(need(!is.null(sc), "No sex column available to subset on."))
        keep_sex <- !is.na(ph0[[sc]]) & as.character(ph0[[sc]]) == sex_value
        validate(need(sum(keep_sex) >= 10, sprintf("Fewer than 10 samples remain after restricting to %s.", sex_label)))
        mat0 <- mat0[, keep_sex, drop = FALSE]
        ph0 <- ph0[keep_sex, , drop = FALSE]
      }

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

      ## Top-variance cap subsets rows rather than producing a keep mask, so it's
      ## appended to the cascade table as an extra row instead of folded into filters.
      cap <- if (identical(input$fs_max_probes %||% "5000", "custom")) (input$fs_max_probes_custom %||% 5000) else as.numeric(input$fs_max_probes %||% 5000)
      n_before_cap <- nrow(beta_scale)
      beta_scale <- methyl_fs_cap_top_variance(beta_scale, round(cap))
      cap_note <- sprintf("Capped to the top %d most-variable CpGs (of %d).", nrow(beta_scale), n_before_cap)
      cascade <- rbind(cascade, data.frame(step = "Top-variance cap", retained = nrow(beta_scale), removed = n_before_cap - nrow(beta_scale)))
      filter_notes <- c(filter_notes, "Top-variance cap" = cap_note)

      imputed_method <- NA_character_
      impute_fallback <- FALSE
      if (isTRUE(input$fs_impute) && anyNA(beta_scale)) {
        imp <- methyl_fs_impute(beta_scale, input$fs_impute_method %||% "median")
        beta_scale <- imp$mat
        imputed_method <- imp$method_used
        ## methyl_fs_impute() silently degrades k-NN -> median when `impute` isn't
        ## installed - flagged here so the result card can say so explicitly.
        impute_fallback <- identical(input$fs_impute_method, "knn") && !grepl("k-NN", imputed_method, fixed = TRUE)
      }

      m_scale <- methyl_beta_to_mvalue(beta_scale)

      list(beta = beta_scale, m = m_scale, grp = grp, ref_group = input$fs_ref_group, comp_group = input$fs_comp_group,
           group_col = input$fs_group_col, array_type = src$array_type, source_label = src$source_label,
           n_probes = nrow(beta_scale), n_samples = ncol(beta_scale), cascade = cascade,
           filter_notes = filter_notes, imputed_method = imputed_method, impute_fallback = impute_fallback,
           anno = anno_result(), run_at = Sys.time(), sex_label = sex_label)
    }

    fs_filter_result_female <- eventReactive(input$run_female_btn, { fs_build_filters("female", fs_sex_value_for("female")) }, ignoreInit = TRUE)
    fs_filter_result_male   <- eventReactive(input$run_male_btn,   { fs_build_filters("male",   fs_sex_value_for("male")) },   ignoreInit = TRUE)
    fs_filter_result_pooled <- eventReactive(input$run_pooled_btn, { fs_build_filters("pooled", NULL) }, ignoreInit = TRUE)

    ## The shared "does this sex currently have a live, non-stale Stage-1
    ## build?" accessor every downstream stage reads instead of a separate
    ## has_run flag - mirrors transcriptomics' fs_result_or_null()
    ## (R/transcriptomics/mod_featureselection.R:884-899). Calling an
    ## eventReactive before its button has ever fired just raises a cheap
    ## req()-style silent stop, caught here as a plain NULL.
    fs_filter_result_for <- function(sex_label) {
      if (isTRUE(fs_stale[[sex_label]])) return(NULL)
      tryCatch(
        switch(sex_label, female = fs_filter_result_female(), male = fs_filter_result_male(), pooled = fs_filter_result_pooled()),
        error = function(e) NULL
      )
    }
    fs_built_sexes <- reactive({ Filter(function(sx) !is.null(fs_filter_result_for(sx)), FS_SEXES) })

    output$fs_filter_result_ui <- renderUI({
      built <- fs_built_sexes()
      if (length(built) == 0) return(mod_methyl_fs_empty_note())
      cards <- stats::setNames(lapply(built, function(sx) {
        r <- fs_filter_result_for(sx)
        div(class = "card",
          div(class = "card-title", icon("circle-check"), sprintf("Filter result - %s", fs_sex_display_label(sx))),
          div(class = "methyl-stats-row", fluidRow(
            valueBox(format(r$n_probes, big.mark = ","), "CpGs retained", icon = icon("dna"), color = "purple", width = 4),
            valueBox(r$n_samples, sprintf("Samples (%s vs %s)", r$comp_group, r$ref_group), icon = icon("users"), color = "light-blue", width = 4),
            valueBox(if (is.na(r$imputed_method)) "None" else r$imputed_method, "Imputation", icon = icon("wand-magic-sparkles"), color = "black", width = 4)
          )),
          if (isTRUE(r$impute_fallback))
            p(class = "empty-note", icon("triangle-exclamation"),
              "k-NN imputation was requested but the \"impute\" package is not installed in this deployment - fell back to per-probe median imputation.")
          else NULL,
          tags$ul(lapply(r$filter_notes, tags$li)),
          withSpinner(plotOutput(ns(paste0("fs_cascade_plot_", sx)), height = 260), color = "#2563EB", type = 6)
        )
      }), built)
      mod_methyl_fs_sex_layout(cards)
    })

    lapply(FS_SEXES, function(sx) {
      output[[paste0("fs_cascade_plot_", sx)]] <- renderPlot({
        r <- fs_filter_result_for(sx); req(r)
        methyl_plot_cascade(r$cascade)
      })
    })

    ## ---------------------------------------------------------------------
    ## Shared technique-stage factory (Stages 2-6) - the technique-layer half
    ## of transcriptomics' build-layer/technique-layer split
    ## (register_sex_technique_outputs(), R/transcriptomics/mod_featureselection.R:1170),
    ## adapted to this file's existing one-Run-button-per-stage convention:
    ## clicking a stage's own Run button fits `fit_fn` against every sex Stage 1
    ## currently has a live build for, instead of tripling every stage's button.
    ## ---------------------------------------------------------------------

    fs_register_technique_stage <- function(run_btn_name, fit_fn) {
      store <- reactiveValues(female = NULL, male = NULL, pooled = NULL)
      lapply(FS_SEXES, function(sx) {
        observeEvent(fs_filter_result_for(sx), store[[sx]] <- NULL, ignoreNULL = FALSE)
      })
      observeEvent(input[[run_btn_name]], {
        for (sx in fs_built_sexes()) {
          r <- fs_filter_result_for(sx)
          store[[sx]] <- tryCatch(fit_fn(r), error = function(e) e)
        }
      }, ignoreInit = TRUE)
      store
    }

    ## Renders a technique stage's up-to-3 per-sex cards (or the empty note).
    ## `card_fn(sx, fit)` builds one sex's `.card` - `fit` may be a caught
    ## `condition` (a real validate()/error failure) rather than a real result;
    ## card_fn is expected to check `inherits(fit, "condition")` and show
    ## conditionMessage(fit) in that case, so a real failure is never confused
    ## with "not run yet" (mirrors transcriptomics' fs_result_error_msg()).
    fs_render_technique_ui <- function(store, card_fn, need_run_note) {
      present <- Filter(function(sx) !is.null(store[[sx]]), FS_SEXES)
      if (length(present) == 0) return(mod_methyl_fs_empty_note(need_run_note))
      cards <- stats::setNames(lapply(present, function(sx) card_fn(sx, store[[sx]])), present)
      mod_methyl_fs_sex_layout(cards)
    }

    ## Shared "run this stage first" prompt shown when no sex has been built yet.
    fs_need_build_note <- "Build at least one sex (Female/Male/Run All pooled) in Data & Filters, then Run this stage."

    ## ---------------------------------------------------------------------
    ## Stage 2: Univariate Selection
    ## ---------------------------------------------------------------------

    output$uni_ui <- renderUI({
      if (length(fs_built_sexes()) == 0) return(mod_methyl_fs_need_filters_note())
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
                "Data & Filters currently defines a reference/comparison group, not an independent continuous phenotype column - this method runs against that group contrast numerically coded, equivalent to a two-group test.")),
            conditionalPanel(condition = sprintf("['anova','kruskal'].includes(input['%s'])", ns("uni_method")),
              p(class = "empty-note", icon("circle-info"),
                "If the Group/phenotype column picked in Data & Filters has more than the two Reference/Comparison levels among the retained samples, this runs a genuine multi-group test across every level of that column; otherwise it falls back to the two-group Reference-vs-Comparison contrast.")),
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
      ## Covariate candidates don't depend on sex - any built stratum's group_col works.
      built <- fs_built_sexes()
      req(length(built) > 0)
      r <- fs_filter_result_for(built[1])
      if (is.null(r) || is.null(fs_active_source()$sheet)) return(NULL)
      cand <- setdiff(colnames(fs_active_source()$sheet), r$group_col)
      if (length(cand) == 0) return(NULL)
      checkboxGroupInput(ns("uni_covariates"), "Covariate adjustment (linear-model methods only)", choices = cand, selected = character(0))
    })

    ## r$beta/r$m/r$grp are restricted to the two Reference/Comparison levels, but
    ## ANOVA/Kruskal-Wallis can use >2; rebuild against the full phenotype column for those.
    fs_uni_multigroup_data <- function(r, method) {
      fallback <- list(beta = r$beta, m = r$m, y = r$grp)
      if (!(method %in% c("anova", "kruskal"))) return(fallback)
      src <- tryCatch(fs_active_source(), error = function(e) NULL)
      if (is.null(src) || is.null(src$sheet) || is.null(r$group_col) || !(r$group_col %in% colnames(src$sheet))) return(fallback)
      sample_ids <- methyl_sheet_sample_ids(src$sheet, colnames(src$mat))
      common <- intersect(colnames(src$mat), sample_ids)
      if (length(common) < 6) return(fallback)
      ph <- as.data.frame(src$sheet)[match(common, sample_ids), , drop = FALSE]
      grp_raw <- as.character(ph[[r$group_col]])
      eligible <- !is.na(grp_raw) & nzchar(grp_raw)
      if (sum(eligible) < 6 || length(unique(grp_raw[eligible])) <= 2) return(fallback)
      probes <- intersect(rownames(r$beta), rownames(src$mat))
      mat_sub <- src$mat[probes, common[eligible], drop = FALSE]
      beta_sub <- if (identical(src$scale, "m")) 2^mat_sub / (1 + 2^mat_sub) else mat_sub
      samp_ok <- methyl_fs_sample_missing_ok(beta_sub, (input$fs_max_na_sample %||% 20) / 100)
      beta_sub <- beta_sub[, samp_ok$keep, drop = FALSE]
      y <- factor(grp_raw[eligible][samp_ok$keep])
      if (nlevels(y) <= 2 || ncol(beta_sub) < 6) return(fallback)
      if (anyNA(beta_sub)) beta_sub <- methyl_fs_impute(beta_sub, input$fs_impute_method %||% "median")$mat
      list(beta = beta_sub, m = methyl_beta_to_mvalue(beta_sub), y = y)
    }

    fs_uni_fit_one <- function(r) {
      is_linear <- input$uni_method %in% c("moderated_t", "t_test", "anova", "linear_regression", "pearson")
      ud <- fs_uni_multigroup_data(r, input$uni_method)
      y <- ud$y
      cov_df <- NULL
      if (is_linear && length(input$uni_covariates %||% character(0)) > 0) {
        sheet <- fs_active_source()$sheet
        ids <- methyl_sheet_sample_ids(sheet, colnames(ud$beta))
        cov_df <- as.data.frame(sheet)[match(colnames(ud$beta), ids), input$uni_covariates, drop = FALSE]
      }
      ranked <- if (is_linear) {
        methyl_fs_univariate_linear(ud$m, y, covariates = cov_df, mode = input$uni_method)
      } else {
        methyl_fs_univariate_rank(ud$beta, y, mode = input$uni_method)
      }
      sel <- methyl_fs_univariate_select(ranked, ud$beta, y, rule = input$uni_rule, threshold = input$uni_threshold,
                                          dbeta_min = input$uni_dbeta_min, top_n = input$uni_top_n)
      list(ranked = sel$ranked, selected_ids = sel$selected_ids, method = input$uni_method,
           n_groups = nlevels(factor(as.character(y))), run_at = Sys.time())
    }

    fs_uni_store <- fs_register_technique_stage("uni_run_btn", fs_uni_fit_one)

    lapply(FS_SEXES, function(sx) {
      output[[paste0("uni_table_", sx)]] <- DT::renderDataTable({
        r <- fs_uni_store[[sx]]; req(r); req(!inherits(r, "condition"))
        DT::datatable(r$ranked, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
          DT::formatSignif(columns = intersect(c("statistic", "p", "fdr", "dbeta", "logfc"), colnames(r$ranked)), digits = 4)
      })
      outputOptions(output, paste0("uni_table_", sx), suspendWhenHidden = FALSE)
      output[[paste0("uni_download_", sx)]] <- downloadHandler(
        filename = function() sprintf("univariate_selection_%s.csv", sx),
        content = function(file) utils::write.csv(fs_uni_store[[sx]]$ranked, file, row.names = FALSE)
      )
    })

    output$uni_result_ui <- renderUI({
      fs_render_technique_ui(fs_uni_store, function(sx, r) {
        if (inherits(r, "condition")) {
          return(div(class = "card", div(class = "card-title", icon("triangle-exclamation"), sprintf("Univariate result - %s", fs_sex_display_label(sx))),
                      p(class = "empty-note", conditionMessage(r))))
        }
        div(class = "card",
          div(class = "card-title", icon("circle-check"), sprintf("Univariate result - %s", fs_sex_display_label(sx))),
          p(strong(length(r$selected_ids)), sprintf(" of %d CpGs selected (%s%s).", nrow(r$ranked), r$method,
                                                       if (r$method %in% c("anova", "kruskal")) sprintf(", %d groups", r$n_groups) else "")),
          div(class = "table-toolbar", downloadButton(ns(paste0("uni_download_", sx)), "Ranked table (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns(paste0("uni_table_", sx)))
        )
      }, fs_need_build_note)
    })

    ## ---------------------------------------------------------------------
    ## Stage 3: Regularization (LASSO / Elastic Net)
    ## ---------------------------------------------------------------------

    output$reg_ui <- renderUI({
      if (length(fs_built_sexes()) == 0) return(mod_methyl_fs_need_filters_note())
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

    fs_reg_fit_one <- function(r) {
      X <- t(r$m)
      w <- if (identical(input$reg_class_weight, "balanced")) {
        tb <- table(r$grp); ww <- max(tb) / tb; unname(ww[as.character(r$grp)])
      } else NULL
      fit <- methyl_fs_lasso_fit(X, r$grp, alpha = input$reg_alpha, nfolds = input$reg_nfolds, lambda_choice = input$reg_lambda_choice,
                                  nlambda = input$reg_nlambda, standardize = input$reg_standardize, weights = w,
                                  seed = input$reg_seed, max_selected = input$reg_max_selected, coef_threshold = input$reg_coef_threshold)
      fit$n_input <- ncol(X); fit$run_at <- Sys.time()
      fit
    }

    fs_reg_store <- fs_register_technique_stage("reg_run_btn", fs_reg_fit_one)

    lapply(FS_SEXES, function(sx) {
      output[[paste0("reg_plot_", sx)]] <- renderPlot({
        r <- fs_reg_store[[sx]]; req(r); req(!inherits(r, "condition")); plot(r$cv)
      })
      output[[paste0("reg_table_", sx)]] <- DT::renderDataTable({
        r <- fs_reg_store[[sx]]; req(r); req(!inherits(r, "condition"))
        DT::datatable(r$ranked, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
          DT::formatSignif(columns = c("coefficient", "abs_coefficient"), digits = 4)
      })
      outputOptions(output, paste0("reg_table_", sx), suspendWhenHidden = FALSE)
      output[[paste0("reg_download_", sx)]] <- downloadHandler(
        filename = function() sprintf("regularization_selected_cpgs_%s.csv", sx),
        content = function(file) utils::write.csv(fs_reg_store[[sx]]$ranked, file, row.names = FALSE)
      )
    })

    output$reg_result_ui <- renderUI({
      fs_render_technique_ui(fs_reg_store, function(sx, r) {
        if (inherits(r, "condition")) {
          return(div(class = "card", div(class = "card-title", icon("triangle-exclamation"), sprintf("Regularization result - %s", fs_sex_display_label(sx))),
                      p(class = "empty-note", conditionMessage(r))))
        }
        div(class = "card",
          div(class = "card-title", icon("circle-check"), sprintf("Regularization result - %s", fs_sex_display_label(sx))),
          p(strong(length(r$selected_ids)), sprintf(" of %d candidate CpGs selected (lambda = %.5g).", r$n_input, r$lambda_used)),
          withSpinner(plotOutput(ns(paste0("reg_plot_", sx)), height = 300), color = "#2563EB", type = 6),
          div(class = "table-toolbar", downloadButton(ns(paste0("reg_download_", sx)), "Selected CpGs (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns(paste0("reg_table_", sx)))
        )
      }, fs_need_build_note)
    })

    ## ---------------------------------------------------------------------
    ## Stage 4: Tree-Based Selection (Random Forest)
    ## ---------------------------------------------------------------------

    output$rf_ui <- renderUI({
      if (length(fs_built_sexes()) == 0) return(mod_methyl_fs_need_filters_note())
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
            checkboxInput(ns("rf_replace"), "Sample with replacement", value = TRUE),
            radioButtons(ns("rf_class_weight"), "Class weighting", choices = c("Equal (default)" = "equal", "Balanced (inverse frequency)" = "balanced"), selected = "equal")
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

    fs_rf_fit_one <- function(r) {
      X <- t(r$m)
      maxnodes <- if (isTRUE(input$rf_maxnodes_unlimited)) NULL else input$rf_maxnodes
      ## classwt is one weight per class (named vector), unlike glmnet's per-observation weights.
      classwt <- if (identical(input$rf_class_weight, "balanced")) {
        tb <- table(r$grp); stats::setNames(as.numeric(max(tb) / tb), names(tb))
      } else NULL
      fit <- methyl_fs_rf_fit(X, r$grp, ntree = input$rf_ntree, mtry_mode = input$rf_mtry_mode, mtry_manual = input$rf_mtry_manual,
                               nodesize = input$rf_nodesize, maxnodes = maxnodes, replace = input$rf_replace, classwt = classwt, seed = input$rf_seed,
                               importance_type = input$rf_importance_type, rule = input$rf_rule, top_n = input$rf_top_n, threshold = input$rf_top_n)
      fit$n_input <- ncol(X); fit$run_at <- Sys.time()
      fit
    }

    fs_rf_store <- fs_register_technique_stage("rf_run_btn", fs_rf_fit_one)

    lapply(FS_SEXES, function(sx) {
      output[[paste0("rf_plot_", sx)]] <- renderPlot({
        r <- fs_rf_store[[sx]]; req(r); req(!inherits(r, "condition"))
        df <- utils::head(r$ranked, 25)
        ggplot(df, aes(x = reorder(cpg, importance), y = importance)) +
          geom_col(fill = ARTHOMIX_COLORS$blue) + coord_flip() +
          labs(x = NULL, y = r$importance_type) + theme_arthomix()
      })
      output[[paste0("rf_table_", sx)]] <- DT::renderDataTable({
        r <- fs_rf_store[[sx]]; req(r); req(!inherits(r, "condition"))
        df <- r$ranked
        df$selected <- df$cpg %in% r$selected_ids
        DT::datatable(df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
          DT::formatSignif(columns = "importance", digits = 4)
      })
      outputOptions(output, paste0("rf_table_", sx), suspendWhenHidden = FALSE)
      output[[paste0("rf_download_", sx)]] <- downloadHandler(
        filename = function() sprintf("random_forest_ranked_cpgs_%s.csv", sx),
        content = function(file) utils::write.csv(fs_rf_store[[sx]]$ranked, file, row.names = FALSE)
      )
    })

    output$rf_result_ui <- renderUI({
      fs_render_technique_ui(fs_rf_store, function(sx, r) {
        if (inherits(r, "condition")) {
          return(div(class = "card", div(class = "card-title", icon("triangle-exclamation"), sprintf("Random Forest result - %s", fs_sex_display_label(sx))),
                      p(class = "empty-note", conditionMessage(r))))
        }
        div(class = "card",
          div(class = "card-title", icon("circle-check"), sprintf("Random Forest result - %s", fs_sex_display_label(sx))),
          p(strong(length(r$selected_ids)), sprintf(" of %d candidate CpGs selected (mtry = %d, %s).", r$n_input, r$mtry, r$importance_type)),
          withSpinner(plotOutput(ns(paste0("rf_plot_", sx)), height = 340), color = "#2563EB", type = 6),
          div(class = "table-toolbar", downloadButton(ns(paste0("rf_download_", sx)), "Ranked table (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns(paste0("rf_table_", sx)))
        )
      }, fs_need_build_note)
    })

    ## ---------------------------------------------------------------------
    ## Stage 5: RFE / Wrapper Selection
    ## ---------------------------------------------------------------------

    output$rfe_ui <- renderUI({
      built <- fs_built_sexes()
      if (length(built) == 0) return(mod_methyl_fs_need_filters_note())
      r <- fs_filter_result_for(built[1])
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

    fs_rfe_fit_one <- function(r) {
      X <- t(r$m)
      sizes <- methyl_fs_rfe_sizes(input$rfe_sizes, ncol(X))
      ## methyl_fs_rfe_sizes() silently drops sizes <= 0 or > ncol(X); surfaced
      ## here so the result card can flag dropped sizes to the user.
      entered <- suppressWarnings(as.integer(trimws(strsplit(input$rfe_sizes %||% "", ",")[[1]])))
      dropped_sizes <- setdiff(entered[!is.na(entered)], sizes)
      fit <- if (identical(input$rfe_flavor, "svm")) {
        methyl_fs_svm_rfe_run(X, r$grp, cost = input$rfe_svm_cost, tolerance = input$rfe_svm_tolerance, k_folds = input$rfe_folds,
                               seed = input$rfe_seed, panel_mode = input$rfe_svm_panel_mode, manual_k = input$rfe_svm_manual_k)
      } else {
        methyl_fs_rfe_run(X, r$grp, sizes = sizes, k = input$rfe_folds, flavor = input$rfe_flavor, seed = input$rfe_seed)
      }
      fit$flavor <- input$rfe_flavor; fit$n_input <- ncol(X); fit$run_at <- Sys.time(); fit$dropped_sizes <- dropped_sizes
      fit
    }

    fs_rfe_store <- fs_register_technique_stage("rfe_run_btn", fs_rfe_fit_one)

    lapply(FS_SEXES, function(sx) {
      output[[paste0("rfe_plot_", sx)]] <- renderPlot({
        r <- fs_rfe_store[[sx]]; req(r); req(!inherits(r, "condition"))
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
      output[[paste0("rfe_table_", sx)]] <- DT::renderDataTable({
        r <- fs_rfe_store[[sx]]; req(r); req(!inherits(r, "condition"))
        DT::datatable(r$ranked, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact")
      })
      outputOptions(output, paste0("rfe_table_", sx), suspendWhenHidden = FALSE)
      output[[paste0("rfe_download_", sx)]] <- downloadHandler(
        filename = function() sprintf("rfe_ranked_cpgs_%s.csv", sx),
        content = function(file) utils::write.csv(fs_rfe_store[[sx]]$ranked, file, row.names = FALSE)
      )
    })

    output$rfe_result_ui <- renderUI({
      fs_render_technique_ui(fs_rfe_store, function(sx, r) {
        if (inherits(r, "condition")) {
          return(div(class = "card", div(class = "card-title", icon("triangle-exclamation"), sprintf("RFE result - %s", fs_sex_display_label(sx))),
                      p(class = "empty-note", conditionMessage(r))))
        }
        div(class = "card",
          div(class = "card-title", icon("circle-check"), sprintf("RFE result - %s", fs_sex_display_label(sx))),
          p(strong(length(r$selected_ids)), sprintf(" of %d candidate CpGs selected (optimal size = %s).", r$n_input, r$optimal_size)),
          if (length(r$dropped_sizes) > 0)
            p(class = "empty-note", icon("triangle-exclamation"),
              sprintf("%d requested size(s) (%s) exceeded the number of candidate CpGs and were dropped.",
                      length(r$dropped_sizes), paste(sort(r$dropped_sizes), collapse = ", ")))
          else NULL,
          withSpinner(plotOutput(ns(paste0("rfe_plot_", sx)), height = 300), color = "#2563EB", type = 6),
          div(class = "table-toolbar", downloadButton(ns(paste0("rfe_download_", sx)), "Ranked table (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns(paste0("rfe_table_", sx)))
        )
      }, fs_need_build_note)
    })

    ## ---------------------------------------------------------------------
    ## Stage 6: Stability Selection
    ## ---------------------------------------------------------------------

    output$stab_ui <- renderUI({
      if (length(fs_built_sexes()) == 0) return(mod_methyl_fs_need_filters_note())
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

    fs_stab_fit_one <- function(r) {
      X <- t(r$m)
      fit <- methyl_fs_stability_run(X, r$grp, type = input$stab_type, n_resamples = input$stab_n_resamples, fraction = input$stab_fraction,
                                      k = input$stab_k, repeats = input$stab_repeats, seed = input$stab_seed, freq_threshold = input$stab_freq_threshold)
      fit$n_input <- ncol(X); fit$run_at <- Sys.time()
      fit
    }

    fs_stab_store <- fs_register_technique_stage("stab_run_btn", fs_stab_fit_one)

    lapply(FS_SEXES, function(sx) {
      output[[paste0("stab_plot_", sx)]] <- renderPlot({
        r <- fs_stab_store[[sx]]; req(r); req(!inherits(r, "condition"))
        df <- utils::head(r$ranked, 30)
        ggplot(df, aes(x = reorder(cpg, selection_frequency), y = selection_frequency)) +
          geom_col(fill = ARTHOMIX_COLORS$aqua) + coord_flip() +
          geom_hline(yintercept = input$stab_freq_threshold, linetype = "dashed", color = "grey50") +
          labs(x = NULL, y = "Selection frequency") + theme_arthomix()
      })
      output[[paste0("stab_table_", sx)]] <- DT::renderDataTable({
        r <- fs_stab_store[[sx]]; req(r); req(!inherits(r, "condition"))
        DT::datatable(r$ranked, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
          DT::formatSignif(columns = c("selection_frequency", "stability_score"), digits = 4)
      })
      outputOptions(output, paste0("stab_table_", sx), suspendWhenHidden = FALSE)
      output[[paste0("stab_download_", sx)]] <- downloadHandler(
        filename = function() sprintf("stability_selection_%s.csv", sx),
        content = function(file) utils::write.csv(fs_stab_store[[sx]]$ranked, file, row.names = FALSE)
      )
    })

    output$stab_result_ui <- renderUI({
      fs_render_technique_ui(fs_stab_store, function(sx, r) {
        if (inherits(r, "condition")) {
          return(div(class = "card", div(class = "card-title", icon("triangle-exclamation"), sprintf("Stability Selection result - %s", fs_sex_display_label(sx))),
                      p(class = "empty-note", conditionMessage(r))))
        }
        div(class = "card",
          div(class = "card-title", icon("circle-check"), sprintf("Stability Selection result - %s", fs_sex_display_label(sx))),
          p(strong(length(r$selected_ids)), sprintf(" of %d candidate CpGs selected (>= %.0f%% of %d valid resamples).", r$n_input, input$stab_freq_threshold * 100, r$n_resamples)),
          if (r$n_resamples < r$n_resamples_requested)
            p(class = "empty-note", icon("triangle-exclamation"),
              sprintf("%d of %d requested resamples were degenerate (a single class) or failed to fit, and were excluded from the frequency calculation above.",
                      r$n_resamples_requested - r$n_resamples, r$n_resamples_requested))
          else NULL,
          withSpinner(plotOutput(ns(paste0("stab_plot_", sx)), height = 300), color = "#2563EB", type = 6),
          div(class = "table-toolbar", downloadButton(ns(paste0("stab_download_", sx)), "Ranked table (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns(paste0("stab_table_", sx)))
        )
      }, fs_need_build_note)
    })

    ## ---------------------------------------------------------------------
    ## Stage 7: Consensus / Overlap (+ correlation reduction)
    ## ---------------------------------------------------------------------

    ## Per-sex method availability/ids - the direct analog of transcriptomics'
    ## "Overlap" tab, computed independently for each sex using only that
    ## sex's own upstream method results.
    fs_available_methods_for <- function(sx) {
      avail <- c(Univariate = !is.null(fs_uni_store[[sx]]) && !inherits(fs_uni_store[[sx]], "condition"),
                 Regularization = !is.null(fs_reg_store[[sx]]) && !inherits(fs_reg_store[[sx]], "condition"),
                 TreeBased = !is.null(fs_rf_store[[sx]]) && !inherits(fs_rf_store[[sx]], "condition"),
                 RFE = !is.null(fs_rfe_store[[sx]]) && !inherits(fs_rfe_store[[sx]], "condition"),
                 Stability = !is.null(fs_stab_store[[sx]]) && !inherits(fs_stab_store[[sx]], "condition"))
      names(avail)[avail]
    }
    fs_method_ids_for <- function(sx) {
      out <- list()
      if ("Univariate" %in% fs_available_methods_for(sx)) out$Univariate <- fs_uni_store[[sx]]$selected_ids
      if ("Regularization" %in% fs_available_methods_for(sx)) out$Regularization <- fs_reg_store[[sx]]$selected_ids
      if ("TreeBased" %in% fs_available_methods_for(sx)) out$TreeBased <- fs_rf_store[[sx]]$selected_ids
      if ("RFE" %in% fs_available_methods_for(sx)) out$RFE <- fs_rfe_store[[sx]]$selected_ids
      if ("Stability" %in% fs_available_methods_for(sx)) out$Stability <- fs_stab_store[[sx]]$selected_ids
      out
    }
    ## Union across every currently-built sex, just to populate the shared
    ## method-picker/weights UI - the actual consensus computed per sex below
    ## only ever uses whichever of these that specific sex actually has.
    fs_consensus_avail_union <- reactive({
      out <- character(0)
      for (sx in fs_built_sexes()) out <- union(out, fs_available_methods_for(sx))
      out
    })

    output$consensus_ui <- renderUI({
      if (length(fs_built_sexes()) == 0) return(mod_methyl_fs_need_filters_note())
      avail <- fs_consensus_avail_union()
      if (length(avail) == 0) {
        return(mod_methyl_fs_empty_note("Run at least one method (Univariate/Regularization/Tree-Based/RFE/Stability) for at least one sex before computing consensus."))
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
        p(class = "empty-note", icon("circle-info"), "Computed independently per sex, using whichever of the checked methods that sex actually has results for."),
        actionButton(ns("consensus_run_btn"), "Run Consensus Selection", icon = icon("play"), class = "btn-primary"),
        uiOutput(ns("consensus_result_ui"))
      )
    })

    output$consensus_weight_ui <- renderUI({
      avail <- fs_consensus_avail_union()
      req(length(avail) > 0)
      tagList(lapply(avail, function(m) numericInput(ns(paste0("consensus_w_", m)), sprintf("Weight: %s", m), value = 1, min = 0, max = 10, step = 0.5)))
    })

    fs_consensus_fit_for <- function(sx) {
      ids_all <- fs_method_ids_for(sx)
      chosen <- intersect(input$consensus_methods %||% names(ids_all), names(ids_all))
      if (length(chosen) == 0) return(simpleError(sprintf("None of the checked methods have a result for %s.", fs_sex_display_label(sx))))
      ids <- ids_all[chosen]
      weights <- stats::setNames(vapply(chosen, function(m) input[[paste0("consensus_w_", m)]] %||% 1, numeric(1)), chosen)
      ct <- methyl_fs_consensus_table(ids, weights = weights)
      min_methods <- min(input$consensus_min_methods %||% 1, length(chosen))
      sel <- methyl_fs_consensus_select(ct$table, min_methods = min_methods,
                                         min_weighted = if (isTRUE(input$consensus_use_weighted)) input$consensus_min_weighted else NULL)
      list(table = ct$table, methods = chosen, weights = weights, selected_ids = sel, min_methods = min_methods, run_at = Sys.time())
    }

    fs_consensus_store <- reactiveValues(female = NULL, male = NULL, pooled = NULL)
    lapply(FS_SEXES, function(sx) {
      observeEvent(fs_filter_result_for(sx), fs_consensus_store[[sx]] <- NULL, ignoreNULL = FALSE)
    })
    observeEvent(input$consensus_run_btn, {
      for (sx in fs_built_sexes()) fs_consensus_store[[sx]] <- tryCatch(fs_consensus_fit_for(sx), error = function(e) e)
      ## Shares one result with other modules via the app-wide `results` cache -
      ## prefers pooled (a single-target contract other modules already expect),
      ## else the first sex that actually produced a result.
      if (!is.null(results)) {
        for (sx in c("pooled", fs_built_sexes())) {
          r <- fs_consensus_store[[sx]]
          if (!is.null(r) && !inherits(r, "condition")) {
            results$featureselection <- list(
              selected_cpgs = r$selected_ids, n_selected = length(r$selected_ids),
              method_counts = vapply(r$methods, function(m) sum(r$table[[m]] == 1), integer(1)),
              consensus_rule = sprintf(">= %d of %d methods", r$min_methods, length(r$methods)), sex_label = sx
            )
            break
          }
        }
      }
    }, ignoreInit = TRUE)

    lapply(FS_SEXES, function(sx) {
      output[[paste0("consensus_venn_", sx)]] <- renderPlot({
        r <- fs_consensus_store[[sx]]; req(r); req(!inherits(r, "condition"))
        draw_overlap_venn(fs_method_ids_for(sx)[r$methods], title = sprintf("%d-CpG consensus", length(r$selected_ids)))
      })
      output[[paste0("consensus_upset_", sx)]] <- renderPlot({
        r <- fs_consensus_store[[sx]]; req(r); req(!inherits(r, "condition"))
        methyl_fs_upset_plot(fs_method_ids_for(sx)[r$methods])
      })
      output[[paste0("consensus_rank_plot_", sx)]] <- renderPlot({
        r <- fs_consensus_store[[sx]]; req(r); req(!inherits(r, "condition"))
        validate(need(nrow(r$table) > 0, "No CpGs met the consensus threshold for this stratum."))
        methyl_fs_consensus_rank_plot(r$table)
      })
      output[[paste0("consensus_heatmap_", sx)]] <- renderPlot({
        r <- fs_consensus_store[[sx]]; req(r); req(!inherits(r, "condition"))
        validate(need(nrow(r$table) > 0, "No CpGs met the consensus threshold for this stratum."))
        methyl_fs_method_heatmap(r$table, r$methods)
      })
      output[[paste0("consensus_table_out_", sx)]] <- DT::renderDataTable({
        r <- fs_consensus_store[[sx]]; req(r); req(!inherits(r, "condition"))
        DT::datatable(r$table, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
          DT::formatSignif(columns = "weighted_score", digits = 4)
      })
      outputOptions(output, paste0("consensus_table_out_", sx), suspendWhenHidden = FALSE)
      output[[paste0("consensus_download_", sx)]] <- downloadHandler(
        filename = function() sprintf("consensus_intersection_table_%s.csv", sx),
        content = function(file) utils::write.csv(fs_consensus_store[[sx]]$table, file, row.names = FALSE)
      )
    })

    output$consensus_result_ui <- renderUI({
      fs_render_technique_ui(fs_consensus_store, function(sx, r) {
        if (inherits(r, "condition")) {
          return(div(class = "card", div(class = "card-title", icon("triangle-exclamation"), sprintf("Consensus - %s", fs_sex_display_label(sx))),
                      p(class = "empty-note", conditionMessage(r))))
        }
        div(class = "card",
          div(class = "card-title", icon("circle-check"), sprintf("Consensus result - %s", fs_sex_display_label(sx))),
          p(strong(length(r$selected_ids)), sprintf(" CpGs meet the consensus threshold (>= %d of %d methods), across %s.",
                                                       r$min_methods, length(r$methods), paste(r$methods, collapse = ", "))),
          fluidRow(
            column(6, tags$h5("Venn diagram"), withSpinner(plotOutput(ns(paste0("consensus_venn_", sx)), height = 320), color = "#2563EB", type = 6)),
            column(6, tags$h5("UpSet-style plot"), withSpinner(plotOutput(ns(paste0("consensus_upset_", sx)), height = 320), color = "#2563EB", type = 6))
          ),
          fluidRow(
            column(6, tags$h5("Consensus ranking"), withSpinner(plotOutput(ns(paste0("consensus_rank_plot_", sx)), height = 320), color = "#2563EB", type = 6)),
            column(6, tags$h5("Method comparison heatmap"), withSpinner(plotOutput(ns(paste0("consensus_heatmap_", sx)), height = 320), color = "#2563EB", type = 6))
          ),
          div(class = "table-toolbar", downloadButton(ns(paste0("consensus_download_", sx)), "Intersection table (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns(paste0("consensus_table_out_", sx)))
        )
      }, "Run at least one method for at least one sex, then Run Consensus Selection.")
    })

    ## ---------------------------------------------------------------------
    ## Stage 8: Selected Features
    ## ---------------------------------------------------------------------

    ## A single trained panel/model is inherently one target - unlike
    ## transcriptomics, which has no downstream "finalize one model" step, this
    ## Stratum selector is the methylomics-specific bridge from 3 parallel
    ## Consensus results down to the one stratum that gets finalized/exported.
    fs_consensus_built_sexes <- reactive({
      Filter(function(sx) !is.null(fs_consensus_store[[sx]]) && !inherits(fs_consensus_store[[sx]], "condition"), FS_SEXES)
    })

    output$fs_stratum_select_ui <- renderUI({
      built <- fs_consensus_built_sexes()
      req(length(built) > 0)
      selectInput(ns("fs_stratum"), "Stratum to finalize / export",
                  choices = stats::setNames(built, vapply(built, fs_sex_display_label, character(1))), selected = built[1])
    })

    ## req() here (not a silent NA fall-through) is load-bearing: this is read
    ## by observeEvent(fs_final_panel(), ...) below, which runs eagerly from
    ## session start regardless of whether this tab's UI has ever been opened.
    ## Before any stratum is built, `built` is character(0) and `built[1]` is
    ## NA - indexing fs_consensus_store[[NA_character_]] throws a real R error
    ## ("key must be not be \"\" or NA"), and an error thrown inside an
    ## observeEvent gets retried on every reactive flush instead of failing
    ## once, which looked like the whole app hanging/dimmed (empirically
    ## verified against a real Playwright-driven upload before this fix).
    fs_selected_stratum <- reactive({
      built <- fs_consensus_built_sexes()
      req(length(built) > 0)
      s <- input$fs_stratum %||% built[1]
      if (!(s %in% built)) built[1] else s
    })

    corr_has_run <- reactiveVal(FALSE)
    observeEvent(fs_selected_stratum(), corr_has_run(FALSE))
    fs_corr_result <- eventReactive(input$corr_reduce_btn, {
      sx <- fs_selected_stratum()
      r <- fs_filter_result_for(sx); cr <- fs_consensus_store[[sx]]
      validate(need(!is.null(r) && !is.null(cr) && !inherits(cr, "condition"), "Run Consensus Selection for this stratum first."))
      score <- stats::setNames(cr$table$weighted_score, cr$table$cpg)
      methyl_fs_correlation_reduce(r$beta, cr$selected_ids, score, method = input$corr_method, r_threshold = input$corr_threshold)
    })
    observeEvent(input$corr_reduce_btn, { corr_has_run(TRUE); shinyjs::show(id = "corr_wrap") }, ignoreInit = TRUE)

    fs_final_panel <- reactive({
      sx <- fs_selected_stratum()
      cr <- fs_consensus_store[[sx]]
      req(!is.null(cr), !inherits(cr, "condition"))
      ids <- if (corr_has_run() && isTRUE(input$use_corr_reduced)) fs_corr_result()$reduced_ids else cr$selected_ids
      list(ids = ids, table = cr$table[cr$table$cpg %in% ids, , drop = FALSE], sex_label = sx)
    })

    output$selected_ui <- renderUI({
      built <- fs_consensus_built_sexes()
      if (length(built) == 0) return(mod_methyl_fs_empty_note("Run Consensus Selection first."))
      panel <- fs_final_panel()
      tagList(
        div(class = "card",
          div(class = "card-title", icon("sliders"), "Stratum"),
          uiOutput(ns("fs_stratum_select_ui")),
          div(class = "card-title", icon("link-slash"), "Correlation-based CpG reduction (optional)"),
          p(class = "submodule-desc", "Reduces the selected stratum's Consensus panel for collinearity - only applied on click."),
          fluidRow(
            column(4, radioButtons(ns("corr_method"), "Correlation method", choices = c("Pearson" = "pearson", "Spearman" = "spearman"), selected = "pearson")),
            column(4, numericInput(ns("corr_threshold"), "Correlation threshold", value = 0.8, min = 0.5, max = 0.99, step = 0.05)),
            column(4, checkboxInput(ns("use_corr_reduced"), "Use the correlation-reduced panel downstream", value = FALSE))
          ),
          actionButton(ns("corr_reduce_btn"), "Reduce for collinearity", icon = icon("compress"), class = "btn-default btn-sm"),
          shinyjs::hidden(div(id = ns("corr_wrap"), uiOutput(ns("corr_result_ui"))))
        ),
        div(class = "card",
          div(class = "card-title", icon("list-check"), sprintf("Selected Features - %s", fs_sex_display_label(panel$sex_label))),
          p(strong(length(panel$ids)), " final CpGs."),
          div(class = "table-toolbar",
              downloadButton(ns("selected_download_csv"), "CSV", class = "btn-sm"),
              downloadButton(ns("selected_download_tsv"), "TSV", class = "btn-sm"),
              downloadButton(ns("selected_copy"), "Copy IDs (TXT)", class = "btn-sm")),
          DT::dataTableOutput(ns("selected_table")),
          shinyjs::hidden(div(id = ns("cpg_detail_wrap"), uiOutput(ns("cpg_detail_ui"))))
        )
      )
    })

    output$corr_result_ui <- renderUI({
      req(corr_has_run())
      r <- fs_corr_result()
      cr <- fs_consensus_store[[fs_selected_stratum()]]
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

    fs_selected_table_full <- reactive({
      panel <- fs_final_panel()
      r <- fs_filter_result_for(panel$sex_label)
      anno_df <- methyl_fs_annotate_panel(panel$ids, r$anno)
      merge(panel$table, anno_df, by = "cpg", all.x = TRUE, sort = FALSE)[order(-panel$table$weighted_score[match(panel$table$cpg, panel$table$cpg)]), , drop = FALSE]
    })

    output$selected_table <- DT::renderDataTable({
      req(length(fs_consensus_built_sexes()) > 0)
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
      sx <- fs_selected_stratum()
      tagList(
        div(class = "card-title", icon("magnifying-glass"), sprintf("Detail: %s", cpg)),
        p(strong("Methods selecting it: "), paste(names(fs_method_ids_for(sx))[vapply(fs_method_ids_for(sx), function(s) cpg %in% s, logical(1))], collapse = ", ")),
        p(strong("Consensus score: "), sprintf("%.3f (%d methods)", row$weighted_score[1], row$n_methods[1])),
        withSpinner(plotOutput(ns("cpg_detail_plot"), height = 260), color = "#2563EB", type = 6),
        DT::dataTableOutput(ns("cpg_detail_table"))
      )
    })

    output$cpg_detail_plot <- renderPlot({
      req(length(input$selected_table_rows_selected) > 0)
      cpg <- fs_selected_table_full()$cpg[input$selected_table_rows_selected]
      r <- fs_filter_result_for(fs_selected_stratum())
      req(r, cpg %in% rownames(r$beta))
      df <- data.frame(sample = colnames(r$beta), beta = r$beta[cpg, ], group = as.character(r$grp))
      ggplot(df, aes(x = group, y = beta, fill = group)) +
        geom_violin(alpha = 0.5) + geom_boxplot(width = 0.15, outlier.size = 0.6) +
        scale_fill_manual(values = arthomix_pair(df$group)) +
        labs(x = NULL, y = "Beta value", fill = NULL) + theme_arthomix()
    })

    output$cpg_detail_table <- DT::renderDataTable({
      req(length(input$selected_table_rows_selected) > 0)
      cpg <- fs_selected_table_full()$cpg[input$selected_table_rows_selected]
      r <- fs_filter_result_for(fs_selected_stratum())
      req(r, cpg %in% rownames(r$beta))
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
      if (length(fs_consensus_built_sexes()) == 0) return(mod_methyl_fs_empty_note("Run Consensus Selection first."))
      tagList(
        p(class = "empty-note", icon("circle-info"),
          sprintf("Validating and exporting the \"%s\" stratum selected in Selected Features.", fs_sex_display_label(fs_selected_stratum()))),
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
            p(class = "empty-note", icon("circle-info"), "Leakage-safe mode reselects the panel using the Univariate + Regularization steps inside every outer fold - not literally the same panel shown in Selected Features - and it does not rerun Tree-Based/RFE/Stability/Consensus per fold, since a full 5-method ensemble refit per fold is impractical in a live session."),
            p(class = "empty-note", icon("triangle-exclamation"), "Data & Filters' own missingness/variance filtering, the top-variance cap, and imputation are still computed once on the full dataset before this outer cross-validation runs, not redone per fold - so \"leakage-safe\" here means the supervised feature-selection step is redone per fold, not that the estimate is completely free of leakage.")),
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
      r <- fs_filter_result_for(fs_selected_stratum()); panel <- fs_final_panel()
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
      sx <- fs_selected_stratum()
      r <- fs_filter_result_for(sx); panel <- fs_final_panel(); cr <- fs_consensus_store[[sx]]
      list(
        schema_version = "1.0", module = "mod_methyl_featureselection", saved_at = Sys.time(), sex_label = sx,
        dataset = list(source_label = r$source_label, array_type = r$array_type, n_probes_input = r$n_probes,
                        n_samples = r$n_samples, group_col = r$group_col, reference_level = r$ref_group, comparison_level = r$comp_group),
        filters = list(notes = r$filter_notes, cascade = r$cascade, imputed_method = r$imputed_method),
        methods = list(
          univariate = if (!is.null(fs_uni_store[[sx]]) && !inherits(fs_uni_store[[sx]], "condition")) list(ran = TRUE, selected_ids = fs_uni_store[[sx]]$selected_ids, ranked = fs_uni_store[[sx]]$ranked) else list(ran = FALSE),
          regularization = if (!is.null(fs_reg_store[[sx]]) && !inherits(fs_reg_store[[sx]], "condition")) list(ran = TRUE, selected_ids = fs_reg_store[[sx]]$selected_ids, coefficients = fs_reg_store[[sx]]$coefficients) else list(ran = FALSE),
          tree_based = if (!is.null(fs_rf_store[[sx]]) && !inherits(fs_rf_store[[sx]], "condition")) list(ran = TRUE, selected_ids = fs_rf_store[[sx]]$selected_ids, ranked = fs_rf_store[[sx]]$ranked) else list(ran = FALSE),
          rfe = if (!is.null(fs_rfe_store[[sx]]) && !inherits(fs_rfe_store[[sx]], "condition")) list(ran = TRUE, selected_ids = fs_rfe_store[[sx]]$selected_ids, optimal_size = fs_rfe_store[[sx]]$optimal_size) else list(ran = FALSE),
          stability = if (!is.null(fs_stab_store[[sx]]) && !inherits(fs_stab_store[[sx]], "condition")) list(ran = TRUE, selected_ids = fs_stab_store[[sx]]$selected_ids, ranked = fs_stab_store[[sx]]$ranked) else list(ran = FALSE)
        ),
        consensus = list(methods_used = cr$methods, min_methods_required = cr$min_methods, weights = cr$weights, selected_ids = cr$selected_ids, table = cr$table),
        correlation_reduction = list(applied = isTRUE(input$use_corr_reduced) && corr_has_run(),
                                       r_threshold = if (corr_has_run()) input$corr_threshold else NA, reduced_ids = if (corr_has_run()) fs_corr_result()$reduced_ids else NULL),
        final_panel = list(cpg_ids = panel$ids, n = length(panel$ids), annotation = fs_selected_table_full()),
        validation = if (validate_has_run()) fs_validate_result() else list(mode = NA),
        session_info = list(r_version = getRversion(), package_versions = list(
          glmnet = as.character(utils::packageVersion("glmnet")), randomForest = as.character(utils::packageVersion("randomForest")),
          caret = as.character(utils::packageVersion("caret")), e1071 = as.character(utils::packageVersion("e1071")),
          limma = as.character(utils::packageVersion("limma")), pROC = as.character(utils::packageVersion("pROC"))))
      )
    })

    output$save_rds <- downloadHandler(
      filename = function() sprintf("methyl_featureselection_%s_%s.rds", fs_selected_stratum(), format(Sys.Date(), "%Y%m%d")),
      content = function(file) saveRDS(fs_model_export(), file)
    )
    output$export_ranking <- downloadHandler(
      filename = function() "feature_ranking.csv",
      content = function(file) utils::write.csv(fs_selected_table_full(), file, row.names = FALSE)
    )
    output$export_consensus <- downloadHandler(
      filename = function() "consensus_table.csv",
      content = function(file) utils::write.csv(fs_consensus_store[[fs_selected_stratum()]]$table, file, row.names = FALSE)
    )
    output$export_summary <- downloadHandler(
      filename = function() "analysis_summary.txt",
      content = function(file) {
        r <- fs_filter_result_for(fs_selected_stratum()); panel <- fs_final_panel()
        lines <- c(
          "Methylomics ML Feature Selection - analysis summary",
          sprintf("Stratum: %s", fs_sex_display_label(panel$sex_label)),
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
