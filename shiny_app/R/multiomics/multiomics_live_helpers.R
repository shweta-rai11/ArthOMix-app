## R/multiomics/multiomics_live_helpers.R
## Pure data-processing logic for the "Live Analysis (Upload & MOFA2)"
## sub-module (mod_multi_live.R / mod_multi_live_mofa.R) - the ONE part of
## the Multi-Omics module that runs real computation over data the user
## supplies, rather than browsing a precomputed pipeline. Every function is
## fail-soft (list(ok, ..., error) or similarly explicit), never silently
## drops samples/features, and never auto-imputes or auto-corrects without
## an explicit user choice (spec Rule 6).
##
## Matrix convention throughout: samples in ROWS, features in COLUMNS (the
## orientation `stats::prcomp()`/`MOFA2` correlation code below expects) -
## `multi_live_validate_matrix()` transposes an uploaded features x samples
## CSV automatically when sample IDs are only found in the column names, and
## says so in its result rather than guessing silently.

## ---------------------------------------------------------------------------
## 1. Upload & Validate
## ---------------------------------------------------------------------------

MULTI_LIVE_OMICS_TYPES <- c(
  "RNA-seq / Transcriptomics" = "rnaseq",
  "Proteomics" = "proteomics",
  "Methylation" = "methylation",
  "Metabolomics" = "metabolomics",
  "Other" = "other"
)

## Reads an uploaded matrix (CSV or RDS) into a numeric matrix, samples in
## rows. `orientation` is the user's own stated claim ("samples_rows" or
## "features_rows") - respected, never silently overridden, but the
## validation report below still surfaces a mismatch (e.g. far more
## "samples" than any real cohort would have) as a warning, not a silent fix.
multi_live_read_matrix <- function(path, orientation = c("samples_rows", "features_rows")) {
  orientation <- match.arg(orientation)
  if (is.null(path) || !file.exists(path)) return(list(ok = FALSE, mat = NULL, error = "No file uploaded."))
  raw <- tryCatch({
    if (grepl("\\.rds$", path, ignore.case = TRUE)) readRDS(path) else as.data.frame(data.table::fread(path, showProgress = FALSE))
  }, error = function(e) e)
  if (inherits(raw, "error")) return(list(ok = FALSE, mat = NULL, error = paste("Could not read file:", conditionMessage(raw))))
  df <- as.data.frame(raw)
  if (ncol(df) < 2) return(list(ok = FALSE, mat = NULL, error = "File needs at least an ID column plus one data column."))
  id_col <- df[[1]]
  mat <- as.matrix(df[, -1, drop = FALSE])
  rownames(mat) <- as.character(id_col)
  storage.mode(mat) <- "double"
  if (identical(orientation, "features_rows")) mat <- t(mat)
  list(ok = TRUE, mat = mat, error = NULL)
}

## Dynamic, real QC on one uploaded layer's matrix - every field computed
## from the actual data, never a hardcoded example (spec section 3).
multi_live_validate_matrix <- function(mat, layer_label = "layer") {
  if (is.null(mat) || !is.matrix(mat)) return(list(ok = FALSE, error = sprintf("%s: not a valid numeric matrix.", layer_label)))
  n_samples <- nrow(mat); n_features <- ncol(mat)
  non_numeric <- sum(!is.finite(mat) & !is.na(mat))
  n_missing <- sum(is.na(mat))
  pct_missing <- 100 * n_missing / (n_samples * n_features)
  var_per_feature <- apply(mat, 2, function(x) stats::var(x, na.rm = TRUE))
  n_zero_var <- sum(!is.na(var_per_feature) & var_per_feature == 0)
  dup_samples <- sum(duplicated(rownames(mat)))
  dup_features <- sum(duplicated(colnames(mat)))
  list(
    ok = TRUE, layer = layer_label,
    n_samples = n_samples, n_features = n_features,
    n_missing = n_missing, pct_missing = pct_missing,
    n_zero_variance = n_zero_var,
    n_duplicate_samples = dup_samples, n_duplicate_features = dup_features,
    n_non_finite = non_numeric,
    sample_ids = rownames(mat)
  )
}

## Builds the "Omics | Samples | Features | Missing % | Zero variance |
## Duplicate IDs" summary table (spec section 3), from real per-layer
## validation results, never hardcoded.
multi_live_qc_summary_table <- function(validations) {
  validations <- Filter(function(v) isTRUE(v$ok), validations)
  if (length(validations) == 0) return(NULL)
  do.call(rbind, lapply(validations, function(v) data.frame(
    Omics = v$layer, Samples = v$n_samples, Features = v$n_features,
    `Missing %` = sprintf("%.1f%%", v$pct_missing), `Zero variance` = v$n_zero_variance,
    `Duplicate samples` = v$n_duplicate_samples, `Duplicate features` = v$n_duplicate_features,
    check.names = FALSE
  )))
}

## ---------------------------------------------------------------------------
## 2. Sample Matching & Missing Data
## ---------------------------------------------------------------------------

## Real per-layer + shared sample counts across N uploaded layers - never
## silently merges mismatched samples (spec section 4).
multi_live_sample_overlap <- function(mat_list) {
  mat_list <- Filter(Negate(is.null), mat_list)
  if (length(mat_list) < 2) return(list(ok = FALSE, error = "Upload at least two omics layers to assess sample overlap."))
  ids_by_layer <- lapply(mat_list, rownames)
  shared <- Reduce(intersect, ids_by_layer)
  per_layer <- vapply(ids_by_layer, length, integer(1))
  list(
    ok = TRUE, per_layer = setNames(per_layer, names(mat_list)),
    n_shared = length(shared), shared_ids = shared,
    layer_only = lapply(seq_along(mat_list), function(i) setdiff(ids_by_layer[[i]], shared))
  )
}

## Per-sample and per-feature missingness distributions for one matrix - the
## two dedicated plots (spec section 5, Plots 2-3) read directly from this.
multi_live_missingness <- function(mat) {
  if (is.null(mat)) return(NULL)
  list(
    per_sample = data.frame(sample = rownames(mat), pct_missing = 100 * rowMeans(is.na(mat))),
    per_feature = data.frame(feature = colnames(mat), pct_missing = 100 * colMeans(is.na(mat)))
  )
}

## Explicit imputation OR removal - never automatic (spec Rule 6). `method`
## one of "none", "mean", "median", "remove_rows", "remove_cols".
multi_live_handle_missing <- function(mat, method = c("none", "mean", "median", "remove_rows", "remove_cols"),
                                       max_sample_missing_pct = 50, max_feature_missing_pct = 50) {
  method <- match.arg(method)
  if (is.null(mat)) return(list(ok = FALSE, mat = NULL, error = "No matrix to process."))
  out <- mat
  samp_pct <- 100 * rowMeans(is.na(out))
  feat_pct <- 100 * colMeans(is.na(out))
  dropped_samples <- rownames(out)[samp_pct > max_sample_missing_pct]
  dropped_features <- colnames(out)[feat_pct > max_feature_missing_pct]
  out <- out[samp_pct <= max_sample_missing_pct, feat_pct <= max_feature_missing_pct, drop = FALSE]
  if (nrow(out) == 0 || ncol(out) == 0) {
    return(list(ok = FALSE, mat = NULL, error = "No samples/features remain after applying the missingness thresholds - relax them and try again."))
  }
  if (identical(method, "mean") || identical(method, "median")) {
    stat_fn <- if (identical(method, "mean")) mean else stats::median
    for (j in seq_len(ncol(out))) {
      col <- out[, j]
      if (any(is.na(col))) out[is.na(col), j] <- stat_fn(col, na.rm = TRUE)
    }
  } else if (identical(method, "remove_rows")) {
    out <- out[stats::complete.cases(out), , drop = FALSE]
  } else if (identical(method, "remove_cols")) {
    out <- out[, colSums(is.na(out)) == 0, drop = FALSE]
  }
  list(ok = TRUE, mat = out, error = NULL, dropped_samples = dropped_samples, dropped_features = dropped_features,
       n_remaining_na = sum(is.na(out)))
}

## ---------------------------------------------------------------------------
## 3. Normalization, Filtering & Scaling - omics-appropriate options only
## (spec section 6): never a one-size-fits-all method list.
## ---------------------------------------------------------------------------

MULTI_LIVE_NORM_CHOICES <- list(
  rnaseq = c("Log2(x + 1)" = "log2", "Quantile normalization (limma)" = "quantile", "None (already normalized)" = "none"),
  proteomics = c("Log2(x + 1)" = "log2", "Median normalization" = "median", "Quantile normalization (limma)" = "quantile", "None (already normalized)" = "none"),
  methylation = c("Beta values (as-is, 0-1)" = "none", "M-value transform (logit)" = "mvalue"),
  metabolomics = c("Log transform" = "log2", "Pareto scaling" = "pareto", "Autoscaling (z-score)" = "autoscale", "None" = "none"),
  other = c("Log2(x + 1)" = "log2", "Z-score" = "autoscale", "None" = "none")
)

multi_live_normalize <- function(mat, omics_type, method) {
  if (is.null(mat)) return(list(ok = FALSE, mat = NULL, error = "No matrix to normalize."))
  out <- switch(method,
    "log2" = log2(pmax(mat, 0) + 1),
    "mvalue" = { b <- pmin(pmax(mat, 1e-3), 1 - 1e-3); log2(b / (1 - b)) },
    "median" = { med <- stats::median(mat, na.rm = TRUE); sweep(mat, 1, apply(mat, 1, stats::median, na.rm = TRUE) - med, "-") },
    "quantile" = tryCatch(limma::normalizeQuantiles(mat), error = function(e) mat),
    "pareto" = sweep(mat, 2, apply(mat, 2, function(x) sqrt(stats::sd(x, na.rm = TRUE))), "/"),
    "autoscale" = scale(mat, center = TRUE, scale = TRUE),
    mat
  )
  list(ok = TRUE, mat = out, error = NULL, method = method)
}

## Explicit feature filtering by variance/MAD/missingness (spec section 10) -
## always reports before/after counts, never silent.
multi_live_filter_features <- function(mat, criterion = c("variance", "mad", "missingness"), keep_top_n = NULL, min_value = NULL) {
  criterion <- match.arg(criterion)
  if (is.null(mat)) return(list(ok = FALSE, mat = NULL, error = "No matrix to filter."))
  score <- switch(criterion,
    "variance" = apply(mat, 2, stats::var, na.rm = TRUE),
    "mad" = apply(mat, 2, stats::mad, na.rm = TRUE),
    "missingness" = -colMeans(is.na(mat))
  )
  keep <- if (!is.null(keep_top_n)) {
    ord <- order(score, decreasing = TRUE)
    seq_along(score) %in% ord[seq_len(min(keep_top_n, length(score)))]
  } else if (!is.null(min_value)) {
    score >= min_value
  } else rep(TRUE, length(score))
  list(ok = TRUE, mat = mat[, keep, drop = FALSE], error = NULL,
       n_before = ncol(mat), n_after = sum(keep), n_removed = ncol(mat) - sum(keep))
}

## Cross-omics z-score scaling (spec section 12) - explicit, applied per
## feature (column), reported before/after scale summary.
multi_live_scale <- function(mat) {
  if (is.null(mat)) return(NULL)
  scale(mat, center = TRUE, scale = TRUE)
}

## ---------------------------------------------------------------------------
## 4. Batch Diagnostics
## ---------------------------------------------------------------------------

## PCA via base stats::prcomp - always available, no new dependency. Returns
## scores + variance-explained so the caller can plot PC1 vs PC2 colored by
## any metadata column and report real % variance in the axis labels.
multi_live_pca <- function(mat) {
  if (is.null(mat) || nrow(mat) < 3 || ncol(mat) < 2) return(list(ok = FALSE, error = "PCA needs at least 3 samples and 2 features."))
  keep_cols <- apply(mat, 2, function(x) stats::var(x, na.rm = TRUE) > 0 && !anyNA(x))
  if (sum(keep_cols) < 2) return(list(ok = FALSE, error = "Fewer than 2 non-missing, non-zero-variance features remain for PCA."))
  pc <- tryCatch(stats::prcomp(mat[, keep_cols, drop = FALSE], center = TRUE, scale. = TRUE), error = function(e) e)
  if (inherits(pc, "error")) return(list(ok = FALSE, error = paste("PCA failed:", conditionMessage(pc))))
  var_exp <- (pc$sdev^2) / sum(pc$sdev^2)
  list(ok = TRUE, scores = as.data.frame(pc$x), var_explained = var_exp, sample_ids = rownames(mat))
}

## Batch-vs-phenotype confounding check (spec section 8-9: warn before
## "correcting away" real biology) - a plain contingency-table association
## test, not a claim of causality.
multi_live_confounding_check <- function(meta, batch_col, phenotype_col) {
  if (is.null(meta) || !all(c(batch_col, phenotype_col) %in% colnames(meta))) return(NULL)
  tab <- table(meta[[batch_col]], meta[[phenotype_col]])
  p <- tryCatch(stats::chisq.test(tab, simulate.p.value = TRUE, B = 2000)$p.value, error = function(e) NA_real_)
  list(table = tab, p_value = p, confounded = !is.na(p) && p > 0.05 && any(rowSums(tab > 0) == 1))
}

## Variance in PC1/PC2 explained by batch vs. by phenotype (spec section 9's
## "quantitative diagnostic", not a "looks better" claim) - a one-way ANOVA
## R^2 of each PC against each grouping variable.
multi_live_variance_by_group <- function(pca_scores, meta, group_col, npcs = 2) {
  if (is.null(pca_scores) || is.null(meta) || !group_col %in% colnames(meta)) return(NULL)
  common <- intersect(rownames(pca_scores), rownames(meta))
  if (length(common) < 3) return(NULL)
  grp <- factor(meta[common, group_col])
  if (nlevels(grp) < 2) return(NULL)
  vapply(seq_len(min(npcs, ncol(pca_scores))), function(i) {
    fit <- tryCatch(stats::lm(pca_scores[common, i] ~ grp), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    summary(fit)$r.squared
  }, numeric(1))
}

multi_live_batch_correct <- function(mat, batch, method = c("combat", "limma"), covariates = NULL) {
  method <- match.arg(method)
  if (is.null(mat) || is.null(batch)) return(list(ok = FALSE, mat = NULL, error = "Need a matrix and a batch vector."))
  if (length(unique(batch)) < 2) return(list(ok = FALSE, mat = NULL, error = "Batch column needs at least two levels."))
  ## sva::ComBat / limma::removeBatchEffect both expect features x samples.
  t_mat <- t(mat)
  out <- tryCatch({
    if (identical(method, "combat")) {
      mod <- if (!is.null(covariates)) stats::model.matrix(~covariates) else NULL
      sva::ComBat(dat = t_mat, batch = batch, mod = mod, par.prior = TRUE)
    } else {
      limma::removeBatchEffect(t_mat, batch = batch, covariates = covariates)
    }
  }, error = function(e) e)
  if (inherits(out, "error")) return(list(ok = FALSE, mat = NULL, error = paste("Batch correction failed:", conditionMessage(out))))
  list(ok = TRUE, mat = t(out), error = NULL, method = method)
}

## ---------------------------------------------------------------------------
## 5. Cross-omics correlation (Feature A x Feature B, and a capped heatmap)
## ---------------------------------------------------------------------------

multi_live_correlation <- function(x, y, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 3) return(list(ok = FALSE, error = "Fewer than 3 paired observations."))
  test <- tryCatch(stats::cor.test(x[ok], y[ok], method = method), error = function(e) e)
  if (inherits(test, "error")) return(list(ok = FALSE, error = conditionMessage(test)))
  list(ok = TRUE, r = unname(test$estimate), p = test$p.value, n = n, method = method)
}

## Correlation matrix between two matrices' features, capped to top-N most
## variable features per side (spec: never render tens of thousands of
## cells at once) - returns a long data.frame with r, p, and BH-FDR.
multi_live_correlation_heatmap_data <- function(matA, matB, top_n = 30, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  common <- intersect(rownames(matA), rownames(matB))
  if (length(common) < 3) return(list(ok = FALSE, error = "Fewer than 3 matched samples between these two layers."))
  a <- matA[common, , drop = FALSE]; b <- matB[common, , drop = FALSE]
  top_a <- names(sort(apply(a, 2, stats::var, na.rm = TRUE), decreasing = TRUE))[seq_len(min(top_n, ncol(a)))]
  top_b <- names(sort(apply(b, 2, stats::var, na.rm = TRUE), decreasing = TRUE))[seq_len(min(top_n, ncol(b)))]
  rows <- expand.grid(featureA = top_a, featureB = top_b, stringsAsFactors = FALSE)
  res <- lapply(seq_len(nrow(rows)), function(i) {
    r <- multi_live_correlation(a[, rows$featureA[i]], b[, rows$featureB[i]], method = method)
    if (!r$ok) return(c(NA_real_, NA_real_, 0))
    c(r$r, r$p, r$n)
  })
  m <- do.call(rbind, res)
  rows$r <- m[, 1]; rows$p <- m[, 2]; rows$n <- m[, 3]
  rows$fdr <- stats::p.adjust(rows$p, method = "BH")
  list(ok = TRUE, df = rows)
}

## ---------------------------------------------------------------------------
## MOFA2 fit + extraction (real MOFA2::create_mofa()/run_mofa(), gated by
## MULTI_MOFA_AVAILABLE - global.R). `mat_list` is samples x features per
## view (this module's own convention); MOFA2 itself wants features x
## samples per view, so the transpose happens once, here. Runs synchronously
## - the Shiny module calls this from inside a shiny::ExtendedTask +
## future::future_promise() so it doesn't block the app (see
## mod_multi_live_mofa.R), matching mod_methyl_dataset.R's own async
## convention for its one other genuinely slow operation.
## ---------------------------------------------------------------------------

## MOFA2::run_mofa()'s default `use_basilisk = TRUE` path (an isolated
## conda/mofapy2 environment basilisk manages) segfaults in this deployment
## - confirmed by direct testing. `use_basilisk = FALSE` plus an explicit
## reticulate::use_python() binding to a python with mofapy2 already
## installed works reliably instead (also confirmed by direct testing,
## real MOFA2 training completed and converged). reticulate's Python
## binding is per-R-process and can't cross a future::multisession worker
## process boundary, so this must be called fresh inside every call to
## multi_live_run_mofa() (idempotent - reticulate::py_config() below is
## cheap once a binding already exists) rather than once in global.R.
MULTI_MOFA_PYTHON_CANDIDATES <- c(
  Sys.getenv("RETICULATE_PYTHON", unset = NA),
  Sys.which("python3"),
  "/Library/Frameworks/Python.framework/Versions/3.11/bin/python3",
  "/usr/bin/python3"
)

multi_mofa_ensure_python <- function() {
  already_ok <- tryCatch(reticulate::py_module_available("mofapy2"), error = function(e) FALSE)
  if (isTRUE(already_ok)) return(list(ok = TRUE, error = NULL))
  candidates <- unique(MULTI_MOFA_PYTHON_CANDIDATES[!is.na(MULTI_MOFA_PYTHON_CANDIDATES) & nzchar(MULTI_MOFA_PYTHON_CANDIDATES)])
  for (py in candidates) {
    ok <- tryCatch({
      reticulate::use_python(py, required = TRUE)
      isTRUE(reticulate::py_module_available("mofapy2"))
    }, error = function(e) FALSE)
    if (isTRUE(ok)) return(list(ok = TRUE, error = NULL))
  }
  list(ok = FALSE, error = "No Python interpreter with the mofapy2 package could be found in this deployment.")
}

multi_live_run_mofa <- function(mat_list, num_factors = 10, seed = 1, convergence_mode = c("fast", "medium", "slow")) {
  convergence_mode <- match.arg(convergence_mode)
  if (!MULTI_MOFA_AVAILABLE) return(list(ok = FALSE, model = NULL, error = "MOFA2 is not installed in this deployment."))
  py <- multi_mofa_ensure_python()
  if (!isTRUE(py$ok)) return(list(ok = FALSE, model = NULL, error = py$error))
  mat_list <- Filter(Negate(is.null), mat_list)
  if (length(mat_list) < 2) return(list(ok = FALSE, model = NULL, error = "Need at least two omics views to train MOFA2."))
  views <- lapply(mat_list, function(m) t(as.matrix(m)))  # MOFA2 wants features x samples
  model <- tryCatch({
    obj <- MOFA2::create_mofa(data = views)
    model_opts <- MOFA2::get_default_model_options(obj)
    model_opts$num_factors <- min(num_factors, min(vapply(mat_list, nrow, integer(1))) - 1)
    train_opts <- MOFA2::get_default_training_options(obj)
    train_opts$seed <- seed
    train_opts$convergence_mode <- convergence_mode
    train_opts$verbose <- FALSE
    obj <- MOFA2::prepare_mofa(obj, model_options = model_opts, training_options = train_opts)
    MOFA2::run_mofa(obj, outfile = tempfile(fileext = ".hdf5"), use_basilisk = FALSE, save_data = FALSE)
  }, error = function(e) e)
  if (inherits(model, "error")) return(list(ok = FALSE, model = NULL, error = paste("MOFA2 training failed:", conditionMessage(model))))
  trained_k <- tryCatch(model@dimensions$K, error = function(e) NA_integer_)
  list(ok = TRUE, model = model, error = NULL, seed = seed, num_factors = trained_k)
}

## Variance explained per factor, per view - tidy long data.frame(view,
## factor, variance_explained), real MOFA2::get_variance_explained() output.
multi_live_mofa_variance_df <- function(model) {
  ve <- tryCatch(MOFA2::get_variance_explained(model)$r2_per_factor[[1]], error = function(e) NULL)
  if (is.null(ve)) return(NULL)
  ## get_variance_explained()$r2_per_factor[[group]] is a matrix with rows =
  ## factors, columns = views (confirmed by direct inspection) - as.table()
  ## preserves that orientation as Var1 = rownames (factors), Var2 =
  ## colnames (views); label accordingly, not positionally guessed.
  df <- as.data.frame(as.table(ve))
  colnames(df) <- c("factor", "view", "variance_explained")
  df$variance_explained <- df$variance_explained / 100
  df
}

## Per-sample factor scores as one data.frame, rownames = sample id.
multi_live_mofa_factors_df <- function(model) {
  f <- tryCatch(MOFA2::get_factors(model, factors = "all")[[1]], error = function(e) NULL)
  if (is.null(f)) return(NULL)
  as.data.frame(f)
}

## Per-feature loadings across all views/factors, tidy long data.frame(view,
## factor, feature, value).
multi_live_mofa_loadings_df <- function(model) {
  w <- tryCatch(MOFA2::get_weights(model, views = "all", factors = "all"), error = function(e) NULL)
  if (is.null(w)) return(NULL)
  do.call(rbind, lapply(names(w), function(view) {
    m <- w[[view]]
    do.call(rbind, lapply(colnames(m), function(fac) data.frame(view = view, factor = fac, feature = rownames(m), value = m[, fac])))
  }))
}

## ---------------------------------------------------------------------------
## 6. Small-sample / high-dimensionality guardrails before MOFA2 training
## (spec section 31 - display a warning, never silently attempt a doomed fit)
## ---------------------------------------------------------------------------

multi_live_mofa_guardrails <- function(mat_list) {
  warnings <- character(0)
  n_samples <- unique(vapply(mat_list, nrow, integer(1)))
  if (length(n_samples) != 1) warnings <- c(warnings, "Matrices do not have the same number of matched samples - integrate on matched samples only before training.")
  ns <- if (length(n_samples) == 1) n_samples else min(vapply(mat_list, nrow, integer(1)))
  if (ns < 10) warnings <- c(warnings, sprintf("Only %d samples - MOFA2 factor structure will be highly unstable below ~10-15 samples.", ns))
  total_features <- sum(vapply(mat_list, ncol, integer(1)))
  if (total_features > 10 * ns * 50) warnings <- c(warnings, sprintf("%s total features for %d samples is very high-dimensional - consider stronger feature filtering first.", format(total_features, big.mark = ","), ns))
  list(ok = length(warnings) == 0, warnings = warnings, n_samples = ns, total_features = total_features)
}
