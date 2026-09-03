## R/multiomics/01_Data_Workspace/multiomics_dataset_helpers.R
## Pure data-processing logic for the "Live Analysis (Upload & MOFA2)"
## sub-module (mod_multi_mofa.R / mod_multi_mofa_engine.R) - the ONE part of

MULTI_LIVE_OMICS_TYPES <- c(
  "RNA-seq / Transcriptomics" = "rnaseq",
  "Proteomics" = "proteomics",
  "Methylation" = "methylation",
  "Metabolomics" = "metabolomics",
  "miRNA" = "mirna",
  "Genomics" = "genomics",
  "Microbiome" = "microbiome",
  "Other" = "other"
)

MULTI_LIVE_FEATURE_ID_LABELS <- c(
  rnaseq = "Gene ID", proteomics = "Protein ID", methylation = "Probe ID",
  metabolomics = "Compound ID", mirna = "miRNA ID", genomics = "Variant/Gene ID",
  microbiome = "Taxon ID", other = "Feature ID"
)

multi_live_read_matrix <- function(path, orientation = c("samples_rows", "features_rows"), filename = NULL) {
  orientation <- match.arg(orientation)
  if (is.null(path) || !file.exists(path)) return(list(ok = FALSE, mat = NULL, error = "No file uploaded."))
  ext <- tolower(tools::file_ext(filename %||% path))
  raw <- tryCatch({
    if (identical(ext, "rds")) {
      loaded <- safe_read_rds(path)
      if (!isTRUE(loaded$ok)) stop(loaded$error)
      loaded$value
    } else if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("openxlsx", quietly = TRUE)) {
        stop("The openxlsx package is not installed in this deployment - export as CSV/TSV instead.")
      }
      openxlsx::read.xlsx(path, sheet = 1)
    } else {
      as.data.frame(data.table::fread(path, showProgress = FALSE))
    }
  }, error = function(e) e)
  if (inherits(raw, "error")) return(list(ok = FALSE, mat = NULL, error = paste("Could not read file:", conditionMessage(raw))))
  df <- as.data.frame(raw)
  if (ncol(df) < 2) return(list(ok = FALSE, mat = NULL, error = "File needs at least an ID column plus one data column."))
  id_col <- df[[1]]
  data_cols <- df[, -1, drop = FALSE]
  was_blank_or_na <- vapply(data_cols, function(col) is.na(col) | !nzchar(trimws(as.character(col))), logical(nrow(data_cols)))
  mat <- as.matrix(data_cols)
  rownames(mat) <- as.character(id_col)
  storage.mode(mat) <- "double"
  n_coerced_na <- sum(is.na(mat) & !was_blank_or_na)
  if (identical(orientation, "features_rows")) mat <- t(mat)
  list(ok = TRUE, mat = mat, error = NULL, n_coerced_na = n_coerced_na)
}

multi_live_detect_table_shape <- function(df) {
  if (is.null(df) || ncol(df) < 3 || nrow(df) < 2) return(list(shape = "wide", confident = FALSE, reason = NULL))
  is_num <- vapply(df, function(col) is.numeric(col) || (is.character(col) && suppressWarnings(!any(is.na(as.numeric(col[!is.na(col) & nzchar(trimws(col))]))))), logical(1))
  n_numeric <- sum(is_num)
  n_id_like <- ncol(df) - n_numeric
  if (n_numeric > 1 && n_numeric >= ncol(df) - 1) return(list(shape = "wide", confident = FALSE, reason = NULL))
  if (n_id_like < 2 || n_numeric == 0) return(list(shape = "wide", confident = FALSE, reason = NULL))
  id_like_cols <- names(df)[!is_num]
  uniq_counts <- vapply(id_like_cols, function(nm) length(unique(df[[nm]])), integer(1))
  if (nrow(df) > 1.5 * max(uniq_counts) && n_numeric <= 2) {
    return(list(shape = "long", confident = TRUE,
                reason = sprintf("%s rows, %d numeric column(s), %d ID column(s) - looks like long/tidy format, not a samples x features matrix.",
                                  format(nrow(df), big.mark = ","), n_numeric, n_id_like)))
  }
  list(shape = "wide", confident = FALSE, reason = NULL)
}

MULTI_LIVE_LONG_FIELD_PATTERNS <- list(
  feature = c("^gene[_ .]?(id|symbol)?$", "^protein[_ .]?id$", "^probe[_ .]?id$", "^compound[_ .]?id$",
              "^taxon[_ .]?id$", "^feature[_ .]?id$", "^feature$", "^symbol$", "^variant[_ .]?id$"),
  sample = c("^sample[_ .]?id$", "^sample$", "^subject[_ .]?id$", "^patient[_ .]?id$", "^sample[_ .]?name$"),
  value = c("^tpm$", "^fpkm$", "^rpkm$", "^counts?$", "^expression$", "^value$", "^intensity$",
            "^abundance$", "^beta$", "^beta[_ .]?value$", "^signal$", "^level$", "^measurement$"),
  group = c("^group$", "^condition$", "^treatment$", "^status$", "^phenotype$", "^class$", "^arm$", "^response$")
)

multi_live_detect_long_columns <- function(df) {
  cols <- colnames(df)
  claimed <- character(0)
  find <- function(field) {
    hit <- NA_character_
    for (p in MULTI_LIVE_LONG_FIELD_PATTERNS[[field]]) {
      m <- grep(p, tolower(trimws(setdiff(cols, claimed))), perl = TRUE)
      if (length(m) >= 1) { hit <- setdiff(cols, claimed)[m[1]]; break }
    }
    if (!is.na(hit)) claimed <<- c(claimed, hit)
    hit
  }
  feature_col <- find("feature"); sample_col <- find("sample")
  value_col <- find("value"); group_col <- find("group")
  if (is.na(value_col)) {
    remaining <- setdiff(cols, claimed)
    is_num <- vapply(df[remaining], function(x) is.numeric(x) || (is.character(x) && suppressWarnings(!any(is.na(as.numeric(x[!is.na(x) & nzchar(trimws(x))]))))), logical(1))
    if (any(is_num)) value_col <- remaining[is_num][1]
  }
  warnings <- character(0)
  if (is.na(feature_col)) warnings <- c(warnings, "No feature/gene/protein ID column detected by name - select it manually below.")
  if (is.na(sample_col)) warnings <- c(warnings, "No sample ID column detected by name - select it manually below.")
  if (is.na(value_col)) warnings <- c(warnings, "No numeric measurement column detected - select it manually below.")
  list(feature_col = feature_col, sample_col = sample_col, value_col = value_col, group_col = group_col, warnings = warnings)
}

multi_live_pivot_long <- function(df, feature_col, sample_col, value_col, group_col = NULL) {
  if (is.null(df) || is.na(feature_col) || is.na(sample_col) || is.na(value_col)) {
    return(list(ok = FALSE, mat = NULL, group_df = NULL, error = "Feature, Sample, and Value columns must all be selected before this dataset can be used."))
  }
  if (!all(c(feature_col, sample_col, value_col) %in% colnames(df))) {
    return(list(ok = FALSE, mat = NULL, group_df = NULL, error = "One or more selected columns are not present in this file."))
  }
  feature <- as.character(df[[feature_col]]); sample <- as.character(df[[sample_col]])
  value <- suppressWarnings(as.numeric(df[[value_col]]))
  keep <- !is.na(feature) & nzchar(feature) & !is.na(sample) & nzchar(sample)
  if (sum(!is.na(value[keep])) == 0) return(list(ok = FALSE, mat = NULL, group_df = NULL, error = sprintf("The selected Value column (\"%s\") has no valid numeric values.", value_col)))
  feature <- feature[keep]; sample <- sample[keep]; value <- value[keep]

  key <- paste(sample, feature, sep = "\r")
  n_dup <- sum(duplicated(key))
  mat <- tapply(value, list(sample, feature), FUN = mean, na.rm = TRUE)
  mat <- as.matrix(mat)

  group_df <- NULL
  if (!is.null(group_col) && !identical(group_col, "(none)") && group_col %in% colnames(df)) {
    grp <- as.character(df[[group_col]])[keep]
    by_sample <- split(grp, sample)
    n_vals <- vapply(by_sample, function(g) length(unique(g[!is.na(g) & nzchar(g)])), integer(1))
    if (any(n_vals > 1)) {
      group_df <- list(df = NULL, warning = sprintf("The Group/Condition column (\"%s\") has more than one value for %d sample(s) - not used as sample metadata.", group_col, sum(n_vals > 1)))
    } else {
      group_val <- vapply(by_sample, function(g) { u <- unique(g[!is.na(g) & nzchar(g)]); if (length(u) == 0) NA_character_ else u[1] }, character(1))
      gdf <- data.frame(row.names = names(group_val), stringsAsFactors = FALSE)
      gdf[[group_col]] <- unname(group_val)
      group_df <- list(df = gdf, warning = NULL)
    }
  }

  list(ok = TRUE, mat = mat, group_df = group_df,
       error = NULL, n_duplicate_pairs = n_dup,
       note = sprintf("Pivoted %s measurement(s) into %d sample(s) x %d feature(s).%s",
                       format(length(value), big.mark = ","), nrow(mat), ncol(mat),
                       if (n_dup > 0) sprintf(" %d duplicate (sample, feature) pair(s) were averaged.", n_dup) else ""))
}

multi_live_detect_orientation <- function(df) {
  if (is.null(df) || ncol(df) < 2 || nrow(df) < 1) return(list(suggested = "samples_rows", confident = FALSE, reason = NULL))
  id_col <- as.character(df[[1]])
  header <- colnames(df)[-1]
  more_cols_than_rows <- length(header) > nrow(df) * 3
  header_id_like <- mean(grepl("^[A-Za-z0-9_.-]+$", header)) > 0.9 && length(unique(header)) == length(header)
  id_col_id_like <- mean(grepl("^[A-Za-z0-9_.-]+$", id_col)) > 0.9 && length(unique(id_col)) == length(id_col)
  if (more_cols_than_rows && header_id_like && !id_col_id_like) {
    return(list(suggested = "features_rows", confident = TRUE,
                reason = sprintf("%s columns vs. %s rows, and the column headers look like sample identifiers - this looks like Feature x Sample orientation.",
                                  format(length(header), big.mark = ","), format(nrow(df), big.mark = ","))))
  }
  list(suggested = "samples_rows", confident = FALSE, reason = NULL)
}

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

multi_live_detect_omics_type <- function(mat) {
  if (is.null(mat) || !is.matrix(mat) || ncol(mat) == 0) return(list(ok = FALSE, detected = "unclassifiable", id_type = NA_character_, value_type = NA_character_, reason = "No feature columns to inspect."))
  id_type <- mcc_detect_id_type(colnames(mat))
  value_type <- mcc_detect_methylation_value_type(mat)
  detected <- if (identical(id_type, "Illumina CpG probe ID")) "methylation"
    else if (id_type %in% c("Gene symbol", "Ensembl Gene ID", "Entrez ID")) "rnaseq"
    else "unclassifiable"
  corroborated <- !(
    (identical(detected, "methylation") && !value_type %in% c("beta", "M-value")) ||
    (identical(detected, "rnaseq") && identical(value_type, "beta"))
  )
  list(
    ok = !identical(detected, "unclassifiable"),
    detected = detected, id_type = id_type, value_type = value_type,
    corroborated = corroborated,
    reason = sprintf("Feature IDs look like %s; value range looks like %s.", id_type, value_type)
  )
}

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

multi_live_sample_overlap <- function(mat_list) {
  mat_list <- Filter(Negate(is.null), mat_list)
  if (length(mat_list) < 2) return(list(ok = FALSE, error = "Upload at least two omics layers to assess sample overlap."))

  ids_by_layer <- lapply(mat_list, rownames)
  per_layer <- vapply(ids_by_layer, length, integer(1))

  dup_exact_by_layer <- lapply(ids_by_layer, function(x) duplicated(x) | duplicated(x, fromLast = TRUE))
  unambig_ids_by_layer <- Map(function(ids, dup) unique(ids[!dup]), ids_by_layer, dup_exact_by_layer)
  exact_shared <- Reduce(intersect, unambig_ids_by_layer)
  exact_dup_ids <- unique(unlist(Map(function(ids, dup) unique(ids[dup]), ids_by_layer, dup_exact_by_layer)))

  norm_by_layer <- lapply(ids_by_layer, ch_normalize_id)
  dup_by_layer <- lapply(norm_by_layer, function(x) duplicated(x) | duplicated(x, fromLast = TRUE))
  unambig_norm_by_layer <- Map(function(nrm, dup) unique(nrm[!dup & nzchar(nrm)]), norm_by_layer, dup_by_layer)
  shared_norm <- Reduce(intersect, unambig_norm_by_layer)
  canonical <- stats::setNames(vapply(shared_norm, function(nrm) ids_by_layer[[1]][norm_by_layer[[1]] == nrm][1], character(1)), shared_norm)

  mats_out <- Map(function(m, ids, nrm) {
    rn <- ids
    norm_hit <- !(rn %in% exact_shared) & nrm %in% shared_norm
    rn[norm_hit] <- unname(canonical[nrm[norm_hit]])
    rownames(m) <- rn
    m
  }, mat_list, ids_by_layer, norm_by_layer)

  shared <- union(exact_shared, unname(canonical[shared_norm]))
  ambiguous <- unique(c(exact_dup_ids, unlist(Map(function(nrm, dup) unique(nrm[dup & nzchar(nrm)]), norm_by_layer, dup_by_layer))))

  list(
    ok = TRUE, mats = mats_out, per_layer = setNames(per_layer, names(mat_list)),
    n_shared = length(shared), shared_ids = shared,
    layer_only = lapply(seq_along(mat_list), function(i) {
      ids <- ids_by_layer[[i]]; nrm <- norm_by_layer[[i]]
      ids[!(ids %in% exact_shared) & !(nrm %in% shared_norm)]
    }),
    n_ambiguous = length(ambiguous), ambiguous_ids = ambiguous
  )
}

multi_live_missingness <- function(mat) {
  if (is.null(mat)) return(NULL)
  list(
    per_sample = data.frame(sample = rownames(mat), pct_missing = 100 * rowMeans(is.na(mat))),
    per_feature = data.frame(feature = colnames(mat), pct_missing = 100 * colMeans(is.na(mat)))
  )
}

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

MULTI_LIVE_NORM_CHOICES <- list(
  rnaseq = c("Log2(x + 1)" = "log2", "Quantile normalization (limma)" = "quantile", "None (already normalized)" = "none"),
  proteomics = c("Log2(x + 1)" = "log2", "Median normalization" = "median", "Quantile normalization (limma)" = "quantile", "None (already normalized)" = "none"),
  methylation = c("Beta values (as-is, 0-1)" = "none", "M-value transform (logit)" = "mvalue"),
  metabolomics = c("Log transform" = "log2", "Pareto scaling" = "pareto", "Autoscaling (z-score)" = "autoscale", "None" = "none"),
  mirna = c("Log2(x + 1)" = "log2", "Quantile normalization (limma)" = "quantile", "None (already normalized)" = "none"),
  genomics = c("None (already coded/normalized)" = "none", "Log2(x + 1)" = "log2"),
  microbiome = c("Log transform" = "log2", "Autoscaling (z-score)" = "autoscale", "None" = "none"),
  other = c("Log2(x + 1)" = "log2", "Z-score" = "autoscale", "None" = "none")
)

multi_live_normalize <- function(mat, omics_type, method) {
  if (is.null(mat)) return(list(ok = FALSE, mat = NULL, error = "No matrix to normalize."))
  out <- switch(method,
    "log2" = log2(pmax(mat, 0) + 1),
    "mvalue" = { b <- pmin(pmax(mat, 1e-3), 1 - 1e-3); log2(b / (1 - b)) },
    "median" = { med <- stats::median(mat, na.rm = TRUE); sweep(mat, 1, apply(mat, 1, stats::median, na.rm = TRUE) - med, "-") },
    "quantile" = tryCatch(t(limma::normalizeQuantiles(t(mat))), error = function(e) mat),
    "pareto" = sweep(mat, 2, apply(mat, 2, function(x) sqrt(stats::sd(x, na.rm = TRUE))), "/"),
    "autoscale" = scale(mat, center = TRUE, scale = TRUE),
    mat
  )
  list(ok = TRUE, mat = out, error = NULL, method = method)
}

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

multi_live_scale <- function(mat) {
  if (is.null(mat)) return(NULL)
  scale(mat, center = TRUE, scale = TRUE)
}

multi_live_pca <- function(mat) {
  if (is.null(mat) || nrow(mat) < 3 || ncol(mat) < 2) return(list(ok = FALSE, error = "PCA needs at least 3 samples and 2 features."))
  keep_cols <- apply(mat, 2, function(x) stats::var(x, na.rm = TRUE) > 0 && !anyNA(x))
  if (sum(keep_cols) < 2) return(list(ok = FALSE, error = "Fewer than 2 non-missing, non-zero-variance features remain for PCA."))
  pc <- tryCatch(stats::prcomp(mat[, keep_cols, drop = FALSE], center = TRUE, scale. = TRUE), error = function(e) e)
  if (inherits(pc, "error")) return(list(ok = FALSE, error = paste("PCA failed:", conditionMessage(pc))))
  var_exp <- (pc$sdev^2) / sum(pc$sdev^2)
  list(ok = TRUE, scores = as.data.frame(pc$x), var_explained = var_exp, sample_ids = rownames(mat))
}

multi_live_confounding_check <- function(meta, batch_col, phenotype_col) {
  if (is.null(meta) || !all(c(batch_col, phenotype_col) %in% colnames(meta))) return(NULL)
  tab <- table(meta[[batch_col]], meta[[phenotype_col]])
  p <- tryCatch(stats::chisq.test(tab, simulate.p.value = TRUE, B = 2000)$p.value, error = function(e) NA_real_)
  list(table = tab, p_value = p, confounded = any(rowSums(tab > 0) == 1))
}

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

multi_live_sample_correlation_data <- function(mat, method = c("pearson", "spearman"), max_samples = 60) {
  method <- match.arg(method)
  if (is.null(mat) || nrow(mat) < 3) return(list(ok = FALSE, error = "Fewer than 3 samples to correlate."))
  truncated <- nrow(mat) > max_samples
  if (truncated) mat <- mat[seq_len(max_samples), , drop = FALSE]
  cm <- tryCatch(stats::cor(t(mat), method = method, use = "pairwise.complete.obs"), error = function(e) NULL)
  if (is.null(cm)) return(list(ok = FALSE, error = "Sample correlation could not be computed."))
  df <- as.data.frame(as.table(cm))
  colnames(df) <- c("featureA", "featureB", "r")
  list(ok = TRUE, df = df, truncated = truncated, n_samples = nrow(mat))
}

multi_live_correlation <- function(x, y, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 3) return(list(ok = FALSE, error = "Fewer than 3 paired observations."))
  test <- tryCatch(stats::cor.test(x[ok], y[ok], method = method), error = function(e) e)
  if (inherits(test, "error")) return(list(ok = FALSE, error = conditionMessage(test)))
  list(ok = TRUE, r = unname(test$estimate), p = test$p.value, n = n, method = method)
}

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
  mofa_outfile <- tempfile(fileext = ".hdf5")
  on.exit(unlink(mofa_outfile), add = TRUE)
  model <- tryCatch({
    obj <- MOFA2::create_mofa(data = views)
    model_opts <- MOFA2::get_default_model_options(obj)
    model_opts$num_factors <- min(num_factors, min(vapply(mat_list, nrow, integer(1))) - 1)
    train_opts <- MOFA2::get_default_training_options(obj)
    train_opts$seed <- seed
    train_opts$convergence_mode <- convergence_mode
    train_opts$verbose <- FALSE
    obj <- MOFA2::prepare_mofa(obj, model_options = model_opts, training_options = train_opts)
    MOFA2::run_mofa(obj, outfile = mofa_outfile, use_basilisk = FALSE, save_data = FALSE)
  }, error = function(e) e)
  if (inherits(model, "error")) return(list(ok = FALSE, model = NULL, error = paste("MOFA2 training failed:", conditionMessage(model))))
  trained_k <- tryCatch(model@dimensions$K, error = function(e) NA_integer_)
  list(ok = TRUE, model = model, error = NULL, seed = seed, num_factors = trained_k)
}

multi_live_mofa_variance_df <- function(model) {
  ve <- tryCatch(MOFA2::get_variance_explained(model)$r2_per_factor[[1]], error = function(e) NULL)
  if (is.null(ve)) return(NULL)
  df <- as.data.frame(as.table(ve))
  colnames(df) <- c("factor", "view", "variance_explained")
  df$variance_explained <- df$variance_explained / 100
  df
}

multi_live_mofa_factors_df <- function(model) {
  f <- tryCatch(MOFA2::get_factors(model, factors = "all")[[1]], error = function(e) NULL)
  if (is.null(f)) return(NULL)
  as.data.frame(f)
}

multi_live_mofa_loadings_df <- function(model) {
  w <- tryCatch(MOFA2::get_weights(model, views = "all", factors = "all"), error = function(e) NULL)
  if (is.null(w)) return(NULL)
  do.call(rbind, lapply(names(w), function(view) {
    m <- w[[view]]
    do.call(rbind, lapply(colnames(m), function(fac) data.frame(view = view, factor = fac, feature = rownames(m), value = m[, fac])))
  }))
}

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

multi_dataset_status <- function(validation, n_shared = NULL, n_own = NULL) {
  if (is.null(validation) || !isTRUE(validation$ok)) {
    return(list(level = "not_compatible", label = "Not Compatible",
                reasons = "The selected file does not contain a usable numeric feature matrix."))
  }
  reasons <- character(0)
  if (validation$n_samples < 3) {
    return(list(level = "not_compatible", label = "Not Compatible",
                reasons = "Fewer than 3 samples detected - too few for matched-sample integration."))
  }
  if (validation$pct_missing > 20) reasons <- c(reasons, sprintf("%.0f%% missing values - review missing-value handling before analysis.", validation$pct_missing))
  if (validation$n_duplicate_samples > 0) reasons <- c(reasons, sprintf("%d duplicate sample ID(s).", validation$n_duplicate_samples))
  if (validation$n_duplicate_features > 0) reasons <- c(reasons, sprintf("%d duplicate feature ID(s).", validation$n_duplicate_features))
  if (validation$n_zero_variance > 0) reasons <- c(reasons, sprintf("%d zero-variance feature(s).", validation$n_zero_variance))
  if (!is.null(n_shared) && !is.null(n_own) && isTRUE(n_own > 0) && (n_shared / n_own) < 0.5) {
    reasons <- c(reasons, sprintf("Only %d of this dataset's %d samples match the other selected datasets - provide a sample mapping or review identifiers.", n_shared, n_own))
  }
  if (length(reasons) > 0) return(list(level = "review", label = "Review Required", reasons = reasons))
  list(level = "ready", label = "Ready", reasons = character(0))
}

multi_dataset_compatibility <- function(validations, overlap = NULL, has_metadata = FALSE) {
  validations <- Filter(Negate(is.null), validations)
  labels <- names(validations)
  per_layer <- lapply(labels, function(nm) {
    v <- validations[[nm]]
    n_shared <- if (!is.null(overlap) && isTRUE(overlap$ok)) overlap$n_shared else NULL
    n_own <- if (!is.null(overlap) && isTRUE(overlap$ok)) overlap$per_layer[[nm]] else NULL
    list(label = nm, status = multi_dataset_status(v, n_shared, n_own))
  })
  names(per_layer) <- labels
  sample_matching_ok <- !is.null(overlap) && isTRUE(overlap$ok) && overlap$n_shared >= 3
  levels <- vapply(per_layer, function(p) p$status$level, character(1))
  overall_label <- if (length(levels) == 0) "NO DATASETS SELECTED"
    else if (any(levels == "not_compatible")) "NOT READY - one or more datasets cannot be used"
    else if (!sample_matching_ok) "REVIEW REQUIRED - insufficient sample matching"
    else if (any(levels == "review")) "READY WITH REVIEW"
    else "READY"
  list(per_layer = per_layer, sample_matching_ok = sample_matching_ok,
       has_metadata = isTRUE(has_metadata), overall_label = overall_label)
}

multi_geo_layer_fetch <- function(accession) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    return(list(ok = FALSE, error = "The GEOquery package is not installed in this deployment. Install it with BiocManager::install(\"GEOquery\"), or upload this dataset as a file instead."))
  }
  acc <- toupper(trimws(accession %||% ""))
  if (!grepl("^GSE[0-9]+$", acc)) return(list(ok = FALSE, error = "Enter a valid GEO Series accession, e.g. GSE12345."))
  gset <- tryCatch(suppressMessages(GEOquery::getGEO(acc, GSEMatrix = TRUE)), error = function(e) e)
  if (inherits(gset, "error") || !is.list(gset) || length(gset) == 0) {
    return(list(ok = FALSE, error = sprintf("Could not fetch %s from GEO - check the accession is a Series (GSExxxxx), not a Sample (GSM) or Platform (GPL) ID.", acc)))
  }
  list(ok = TRUE, accession = acc, platforms = gset, error = NULL)
}

multi_geo_platform_matrix <- function(eset, collapse_genes = TRUE) {
  ex <- tryCatch(Biobase::exprs(eset), error = function(e) NULL)
  if (is.null(ex) || nrow(ex) == 0 || ncol(ex) == 0) {
    return(list(ok = FALSE, error = "No expression matrix here - download the raw file and upload it instead."))
  }
  used_collapse <- FALSE
  if (isTRUE(collapse_genes)) {
    collapsed <- tryCatch(collapse_probes_to_genes(eset), error = function(e) NULL)
    if (!is.null(collapsed) && nrow(collapsed) > 0 && isTRUE(attr(collapsed, "collapsed"))) { ex <- collapsed; used_collapse <- TRUE }
  }
  mat <- t(as.matrix(ex))
  storage.mode(mat) <- "double"
  meta <- tryCatch(as.data.frame(Biobase::pData(eset)), error = function(e) NULL)
  list(ok = TRUE, mat = mat, meta = meta,
       platform = tryCatch(Biobase::annotation(eset), error = function(e) NA_character_),
       collapsed = used_collapse)
}
