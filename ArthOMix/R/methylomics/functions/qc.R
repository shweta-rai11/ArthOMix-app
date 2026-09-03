## R/methylomics/functions/qc.R
## Probe- and sample-level QC for the Methylomics module (mod_methyl_qc.R).
## Each probe filter returns list(keep = <logical>, note = <summary>) instead of filtering in place,

methyl_row_vars <- function(m) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowVars(m, na.rm = TRUE)
  } else {
    apply(m, 1, function(x) stats::var(x, na.rm = TRUE))
  }
}

methyl_filter_missing <- function(mat, max_na_frac = 0) {
  na_frac <- rowMeans(is.na(mat))
  keep <- na_frac <= max_na_frac
  list(keep = keep, note = sprintf("%d probe(s) exceed %.0f%% missing values across samples.", sum(!keep), max_na_frac * 100))
}

methyl_filter_variance <- function(mat, min_variance = 0) {
  v <- methyl_row_vars(mat)
  v[!is.finite(v)] <- 0
  keep <- v >= min_variance
  list(keep = keep, note = sprintf("%d probe(s) below the variance threshold (%.4g).", sum(!keep), min_variance))
}

methyl_filter_sd <- function(mat, min_sd = 0) {
  sd <- sqrt(pmax(methyl_row_vars(mat), 0))
  sd[!is.finite(sd)] <- 0
  keep <- sd >= min_sd
  list(keep = keep, note = sprintf("%d probe(s) below the standard-deviation threshold (%.4g).", sum(!keep), min_sd))
}

methyl_filter_mean_range <- function(mat, lo, hi) {
  m <- rowMeans(mat, na.rm = TRUE)
  keep <- !is.na(m) & m >= lo & m <= hi
  list(keep = keep, note = sprintf("%d probe(s) with mean value outside [%.3g, %.3g].", sum(!keep), lo, hi))
}

methyl_filter_non_cpg <- function(mat) {
  keep <- methyl_probe_is_cpg(rownames(mat))
  list(keep = keep, note = sprintf("%d non-CpG (CpH / control) probe(s) removed by ID prefix.", sum(!keep)))
}

methyl_filter_snp <- function(mat, anno_result) {
  if (!isTRUE(anno_result$ok)) return(list(keep = rep(TRUE, nrow(mat)), note = anno_result$reason))
  a <- anno_result$anno
  ids <- rownames(mat)
  hit <- ids %in% rownames(a)
  has_snp <- rep(FALSE, length(ids))
  cols <- intersect(c("Probe_rs", "CpG_rs", "SBE_rs"), colnames(a))
  if (length(cols) > 0 && any(hit)) {
    sub <- a[ids[hit], cols, drop = FALSE]
    has_snp[hit] <- apply(sub, 1, function(r) any(!is.na(r) & nzchar(as.character(r))))
  }
  keep <- !has_snp
  list(keep = keep, note = sprintf("%d probe(s) overlap a known SNP (manifest Probe_rs/CpG_rs/SBE_rs).", sum(has_snp)))
}

methyl_filter_sex_chr <- function(mat, anno_result, mode = "remove_xy") {
  if (identical(mode, "keep")) {
    return(list(keep = rep(TRUE, nrow(mat)), note = "Sex-chromosome probes kept (no filtering)."))
  }
  if (!isTRUE(anno_result$ok)) return(list(keep = rep(TRUE, nrow(mat)), note = anno_result$reason))
  a <- anno_result$anno
  ids <- rownames(mat)
  hit <- ids %in% rownames(a)
  chr <- rep(NA_character_, length(ids))
  chr[hit] <- a[ids[hit], "chr"]
  target_chr <- if (identical(mode, "remove_y_only")) c("chrY", "Y") else c("chrX", "chrY", "X", "Y")
  is_sex <- !is.na(chr) & chr %in% target_chr
  keep <- !is_sex
  label <- if (identical(mode, "remove_y_only")) "chrY" else "chrX/chrY"
  list(keep = keep, note = sprintf("%d probe(s) on %s removed.", sum(is_sex), label))
}

methyl_filter_cross_reactive <- function(mat, exclusion_ids = NULL) {
  if (is.null(exclusion_ids) || length(exclusion_ids) == 0) {
    return(list(keep = rep(TRUE, nrow(mat)), note =
      "No cross-reactive probe list is bundled in this deployment - upload a probe-exclusion list (one probe ID per line, e.g. a published Chen et al. 2013 / Pidsley et al. 2016 / McCartney et al. 2016 list) to enable this filter."))
  }
  keep <- !(rownames(mat) %in% exclusion_ids)
  list(keep = keep, note = sprintf("%d probe(s) removed via the uploaded exclusion list.", sum(!keep)))
}

methyl_parse_maf_list <- function(datapath, filename) {
  df <- tryCatch(as.data.frame(data.table::fread(datapath, showProgress = FALSE)), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0 || ncol(df) < 2) {
    return(list(ok = FALSE, error = "Could not parse this file as a probe_id,maf table (CSV/TSV, probe ID in the first column, allele frequency in the second)."))
  }
  maf_col <- intersect(c("maf", "MAF", "af", "AF"), colnames(df))
  maf_vals <- suppressWarnings(as.numeric(df[[if (length(maf_col) > 0) maf_col[1] else 2]]))
  if (all(is.na(maf_vals))) {
    return(list(ok = FALSE, error = "No numeric MAF/allele-frequency column found (expected a column named maf/MAF/af/AF, or a numeric second column)."))
  }
  ids <- as.character(df[[1]])
  ok <- !is.na(maf_vals) & nzchar(ids)
  list(ok = TRUE, maf = stats::setNames(maf_vals[ok], ids[ok]))
}

methyl_filter_maf <- function(mat, maf_table = NULL, max_maf = 0.05) {
  if (is.null(maf_table) || length(maf_table) == 0) {
    return(list(keep = rep(TRUE, nrow(mat)), note =
      "No MAF table uploaded - upload a probe_id,maf list (e.g. a published population allele-frequency table for this array type) to enable this filter."))
  }
  ids <- rownames(mat)
  hit <- ids %in% names(maf_table)
  above <- rep(FALSE, length(ids))
  above[hit] <- maf_table[ids[hit]] > max_maf
  keep <- !above
  list(keep = keep, note = sprintf("%d probe(s) exceed MAF %.3g in the uploaded table (%d probe(s) had no match and were kept).", sum(above), max_maf, sum(!hit)))
}

methyl_filter_detection_p <- function(mat, detp, threshold = 0.01) {
  if (is.null(detp)) {
    return(list(keep = rep(TRUE, nrow(mat)), note =
      "Detection p-values require raw IDAT input - not available for an uploaded beta/M-value matrix."))
  }
  common <- intersect(rownames(mat), rownames(detp))
  frac_fail <- stats::setNames(rep(0, nrow(mat)), rownames(mat))
  frac_fail[common] <- rowMeans(detp[common, , drop = FALSE] > threshold, na.rm = TRUE)
  keep <- frac_fail == 0
  list(keep = keep, note = sprintf("%d probe(s) fail detection p < %.3g in at least one sample.", sum(!keep), threshold))
}

methyl_filter_beadcount <- function(mat, beadcount, threshold = 3) {
  if (is.null(beadcount)) {
    return(list(keep = rep(TRUE, nrow(mat)), note =
      "Bead counts require raw IDAT input - not available for an uploaded beta/M-value matrix."))
  }
  common <- intersect(rownames(mat), rownames(beadcount))
  low <- stats::setNames(rep(FALSE, nrow(mat)), rownames(mat))
  low[common] <- rowMeans(beadcount[common, , drop = FALSE] < threshold, na.rm = TRUE) > 0
  keep <- !low
  list(keep = keep, note = sprintf("%d probe(s) have bead count < %d in at least one sample.", sum(low), threshold))
}

methyl_sample_call_rate <- function(mat) {
  1 - colMeans(is.na(mat))
}

methyl_sample_failed_probe_pct <- function(mat, detp, threshold = 0.01) {
  if (is.null(detp)) {
    return(list(ok = FALSE, reason = "Failed-probe percentage requires raw IDAT input - not available for an uploaded beta/M-value matrix."))
  }
  common_probes <- intersect(rownames(mat), rownames(detp))
  common_samples <- intersect(colnames(mat), colnames(detp))
  if (length(common_probes) == 0 || length(common_samples) == 0) {
    return(list(ok = FALSE, reason = "No overlapping probes/samples between the detection p-value matrix and this run."))
  }
  pct <- stats::setNames(rep(NA_real_, ncol(mat)), colnames(mat))
  pct[common_samples] <- colMeans(detp[common_probes, common_samples, drop = FALSE] > threshold, na.rm = TRUE) * 100
  list(ok = TRUE, pct = pct)
}

methyl_sample_low_intensity <- function(median_int_result, min_intensity = 10) {
  if (!isTRUE(median_int_result$ok)) {
    return(list(ok = FALSE, reason = median_int_result$reason %||% "Median intensity requires raw IDAT input."))
  }
  d <- median_int_result$detail
  score <- (d$med_meth_log2 + d$med_unmeth_log2) / 2
  list(ok = TRUE, score = stats::setNames(score, d$sample), low = stats::setNames(score < min_intensity, d$sample))
}

methyl_sample_outliers_correlation <- function(mat, n_features = 5000, k = 3) {
  cr <- methyl_sample_correlation(mat, n_features = n_features)
  if (!isTRUE(cr$ok)) return(list(ok = FALSE, reason = cr$reason))
  diag(cr$cor) <- NA
  mean_cor <- rowMeans(cr$cor, na.rm = TRUE)
  med <- stats::median(mean_cor); mad_val <- stats::mad(mean_cor)
  lo <- med - k * mad_val
  list(ok = TRUE, mean_correlation = mean_cor, outlier = mean_cor < lo, threshold = lo)
}

methyl_sample_outliers_pca <- function(mat, n_features = 5000, sd_threshold = 3) {
  m <- stats::na.omit(mat)
  if (nrow(m) < 10 || ncol(m) < 4) {
    return(list(ok = FALSE, reason = "Not enough complete-case probes/samples for PCA-based outlier detection."))
  }
  v <- methyl_row_vars(m)
  top <- order(v, decreasing = TRUE)[seq_len(min(n_features, nrow(m)))]
  pc <- stats::prcomp(t(m[top, , drop = FALSE]), scale. = TRUE)
  pc1 <- pc$x[, 1]; pc2 <- pc$x[, 2]
  dist <- sqrt(((pc1 - mean(pc1)) / stats::sd(pc1))^2 + ((pc2 - mean(pc2)) / stats::sd(pc2))^2)
  list(ok = TRUE, scores = pc$x[, seq_len(min(4, ncol(pc$x))), drop = FALSE],
       outlier = dist > sd_threshold, distance = dist)
}

methyl_sample_outliers_hclust <- function(mat, n_features = 5000, height_frac = 0.5) {
  m <- stats::na.omit(mat)
  if (nrow(m) < 10 || ncol(m) < 4) {
    return(list(ok = FALSE, reason = "Not enough complete-case probes/samples for hierarchical-clustering outlier detection."))
  }
  v <- methyl_row_vars(m)
  top <- order(v, decreasing = TRUE)[seq_len(min(n_features, nrow(m)))]
  d <- stats::dist(t(m[top, , drop = FALSE]))
  hc <- stats::hclust(d, method = "average")
  cl <- stats::cutree(hc, h = max(hc$height) * height_frac)
  tbl <- table(cl)
  singleton <- names(tbl[tbl == 1])
  list(ok = TRUE, hc = hc, outlier = cl %in% as.integer(singleton))
}

methyl_sample_outliers_mahalanobis <- function(mat, n_features = 5000, n_pcs = 10, alpha = 0.01) {
  m <- stats::na.omit(mat)
  if (nrow(m) < 10 || ncol(m) < 4) {
    return(list(ok = FALSE, reason = "Not enough complete-case probes/samples for Mahalanobis-distance outlier detection."))
  }
  v <- methyl_row_vars(m)
  top <- order(v, decreasing = TRUE)[seq_len(min(n_features, nrow(m)))]
  pc <- stats::prcomp(t(m[top, , drop = FALSE]), scale. = TRUE)
  k <- min(n_pcs, ncol(pc$x), ncol(m) - 2)
  if (k < 2) {
    return(list(ok = FALSE, reason = "Not enough samples relative to components for Mahalanobis-distance outlier detection."))
  }
  scores <- pc$x[, seq_len(k), drop = FALSE]
  d2 <- tryCatch(stats::mahalanobis(scores, colMeans(scores), stats::cov(scores)), error = function(e) NULL)
  if (is.null(d2)) {
    return(list(ok = FALSE, reason = "The covariance matrix of the top components is singular - try fewer components or more samples."))
  }
  thresh <- stats::qchisq(1 - alpha, df = k)
  list(ok = TRUE, distance2 = d2, threshold = thresh, outlier = d2 > thresh)
}

methyl_cluster_sex <- function(y, reported = NULL) {
  cluster <- NULL
  if (length(y) >= 6 && length(unique(y)) >= 2) {
    km <- tryCatch(stats::kmeans(y, centers = 2, nstart = 25), error = function(e) NULL)
    if (!is.null(km)) cluster <- km$cluster
  }
  if (is.null(cluster)) cluster <- ifelse(y > stats::median(y, na.rm = TRUE), 2L, 1L)
  hi <- as.integer(names(which.max(tapply(y, cluster, mean))))

  hi_sex <- "M"; lo_sex <- "F"; direction_assumed <- TRUE
  if (!is.null(reported)) {
    agree_hi <- as.character(reported)[cluster == hi]
    agree_hi <- agree_hi[!is.na(agree_hi) & agree_hi %in% c("M", "F")]
    if (length(agree_hi) > 0) {
      hi_sex <- names(sort(table(agree_hi), decreasing = TRUE))[1]
      lo_sex <- setdiff(c("M", "F"), hi_sex)[1]
      direction_assumed <- FALSE
    }
  }
  list(sex = ifelse(cluster == hi, hi_sex, lo_sex), direction_assumed = direction_assumed)
}

methyl_sex_check_attach_mismatch <- function(detail, reported_sex, method) {
  n_mismatch <- NA_integer_
  if (!is.null(reported_sex)) {
    rs <- as.character(reported_sex[detail$sample])
    detail$reported_sex <- rs
    detail$sex_mismatch <- !is.na(rs) & rs != detail$predicted_sex
    n_mismatch <- sum(detail$sex_mismatch, na.rm = TRUE)
  }
  list(ok = TRUE, method = method, detail = detail, n_mismatch = n_mismatch)
}

methyl_sex_check <- function(mat, anno_result, rg_set = NULL, reported_sex = NULL) {
  if (!is.null(rg_set) && requireNamespace("minfi", quietly = TRUE)) {
    sex <- tryCatch({
      gmset <- minfi::mapToGenome(minfi::preprocessRaw(rg_set))
      minfi::getSex(gmset)
    }, error = function(e) NULL)
    if (!is.null(sex)) {
      df <- as.data.frame(sex)
      detail <- data.frame(sample = rownames(df), predicted_sex = df$predictedSex,
                            chrX_median = df$xMed, chrY_median = df$yMed, row.names = NULL)
      return(methyl_sex_check_attach_mismatch(detail, reported_sex,
             "minfi::getSex() - copy-number based, from raw IDAT intensities"))
    }
  }
  if (!isTRUE(anno_result$ok)) {
    return(list(ok = FALSE, reason = "Sex check requires either raw IDAT input, or manifest annotation for this array type (unavailable here)."))
  }
  a <- anno_result$anno
  ids <- rownames(mat)
  hit <- ids %in% rownames(a)
  chr <- rep(NA_character_, length(ids)); chr[hit] <- a[ids[hit], "chr"]
  x_probes <- ids[!is.na(chr) & chr %in% c("chrX", "X")]
  y_probes <- ids[!is.na(chr) & chr %in% c("chrY", "Y")]
  if (length(x_probes) < 10 || length(y_probes) < 5) {
    return(list(ok = FALSE, reason = "Too few annotated chrX/chrY probes in this matrix to estimate sex."))
  }
  mean_x <- colMeans(mat[x_probes, , drop = FALSE], na.rm = TRUE)
  mean_y <- colMeans(mat[y_probes, , drop = FALSE], na.rm = TRUE)

  sample_ids <- colnames(mat)
  clustered <- methyl_cluster_sex(mean_y, if (!is.null(reported_sex)) reported_sex[sample_ids] else NULL)
  method <- sprintf(
    "chrY-methylation two-cluster split (k-means) - weaker than the raw-intensity method; upload IDAT for minfi::getSex(). Cluster-to-sex direction %s.",
    if (clustered$direction_assumed) "assumed (higher chrY beta = male) - upload a sample sheet with a sex column to resolve this by concordance instead"
    else "resolved by majority concordance with the uploaded sample sheet's reported sex"
  )
  detail <- data.frame(sample = sample_ids, mean_chrX = mean_x, mean_chrY = mean_y,
                        predicted_sex = clustered$sex, row.names = NULL)
  methyl_sex_check_attach_mismatch(detail, reported_sex, method)
}

methyl_sheet_sample_ids <- function(sheet, all_ids) {
  id_col <- intersect(c("sample", "Sample", "sample_id", "Sample_ID"), colnames(sheet))
  if (length(id_col) > 0) return(as.character(sheet[[id_col[1]]]))
  if (nrow(sheet) == length(all_ids)) return(all_ids)
  rownames(sheet)
}

methyl_qc_subgroup_filter <- function(mat, sheet, group_col, level, min_n = 3) {
  all_ids <- colnames(mat)
  has_col <- !is.null(sheet) && !is.null(group_col) && nzchar(group_col) && group_col %in% colnames(sheet)
  group_labels <- stats::setNames(rep(NA_character_, length(all_ids)), all_ids)
  if (has_col) {
    sample_ids <- methyl_sheet_sample_ids(sheet, all_ids)
    group_labels <- stats::setNames(as.character(sheet[[group_col]]), sample_ids)[all_ids]
    names(group_labels) <- all_ids
  }
  if (!has_col || is.null(level) || identical(level, "__all__") || !nzchar(level)) {
    return(list(mat = mat, label = sprintf("All samples (n=%d)", length(all_ids)),
                included = all_ids, excluded = character(0),
                group_labels = group_labels, low_n = length(all_ids) < min_n))
  }
  keep <- !is.na(group_labels) & group_labels == level
  included <- all_ids[keep]
  list(mat = mat[, included, drop = FALSE],
       label = sprintf("%s (n=%d)", level, length(included)),
       included = included, excluded = setdiff(all_ids, included),
       group_labels = group_labels,
       low_n = length(included) < min_n)
}

methyl_apply_manual_exclude <- function(subgroup, excluded_ids) {
  if (is.null(excluded_ids) || length(excluded_ids) == 0) return(subgroup)
  removed <- intersect(subgroup$included, excluded_ids)
  keep <- setdiff(subgroup$included, excluded_ids)
  subgroup$mat <- subgroup$mat[, keep, drop = FALSE]
  subgroup$excluded <- union(subgroup$excluded, removed)
  subgroup$included <- keep
  if (length(removed) > 0) {
    subgroup$label <- sprintf("%s, %d manually excluded", subgroup$label, length(removed))
  }
  subgroup$low_n <- length(keep) < 3
  subgroup
}

methyl_outlier_score_table <- function(sample_qc) {
  flag_cols <- intersect(c("pca_outlier", "hclust_outlier", "mahalanobis_outlier", "correlation_outlier", "iqr_outlier"), colnames(sample_qc))
  if (length(flag_cols) == 0) {
    return(data.frame(sample = sample_qc$sample, outlier_score = NA_integer_, n_methods = 0L))
  }
  flags <- as.matrix(sample_qc[, flag_cols, drop = FALSE])
  score <- rowSums(flags == TRUE, na.rm = TRUE)
  out <- data.frame(sample = sample_qc$sample, sample_qc[, flag_cols, drop = FALSE],
                     outlier_score = score, n_methods = length(flag_cols), row.names = NULL)
  out[order(-out$outlier_score), , drop = FALSE]
}

methyl_probe_retention_cascade <- function(n_probes_start, filters) {
  if (length(filters) == 0) {
    return(data.frame(step = "No filters applied", retained = n_probes_start, removed = 0L))
  }
  keep <- rep(TRUE, n_probes_start)
  steps <- character(0); retained <- integer(0)
  for (nm in names(filters)) {
    keep <- keep & filters[[nm]]$keep
    steps <- c(steps, nm)
    retained <- c(retained, sum(keep))
  }
  data.frame(step = c("Start", steps), retained = c(n_probes_start, retained),
             removed = c(NA_integer_, -diff(c(n_probes_start, retained))))
}

methyl_beta_to_mvalue <- function(beta, eps = 1e-4) {
  b <- pmin(pmax(beta, eps), 1 - eps)
  log2(b / (1 - b))
}

methyl_mvalue_to_beta <- function(mvalue) {
  2^mvalue / (1 + 2^mvalue)
}

methyl_pca_scores <- function(mat, n_features = 5000, n_pcs = 10) {
  m <- stats::na.omit(mat)
  if (nrow(m) < 10 || ncol(m) < 4) {
    return(list(ok = FALSE, reason = "Not enough complete-case probes/samples for PCA."))
  }
  v <- methyl_row_vars(m)
  top <- order(v, decreasing = TRUE)[seq_len(min(n_features, nrow(m)))]
  pc <- stats::prcomp(t(m[top, , drop = FALSE]), scale. = TRUE)
  k <- min(n_pcs, ncol(pc$x))
  var_explained <- (pc$sdev^2 / sum(pc$sdev^2))[seq_len(k)]
  list(ok = TRUE, scores = pc$x[, seq_len(k), drop = FALSE], var_explained = var_explained)
}

methyl_sample_correlation <- function(mat, n_features = 5000) {
  m <- stats::na.omit(mat)
  if (nrow(m) < 10 || ncol(m) < 2) {
    return(list(ok = FALSE, reason = "Not enough complete-case probes/samples for a correlation heatmap."))
  }
  v <- methyl_row_vars(m)
  top <- order(v, decreasing = TRUE)[seq_len(min(n_features, nrow(m)))]
  list(ok = TRUE, cor = stats::cor(m[top, , drop = FALSE]))
}

methyl_beta_density_sample <- function(mat, n_probes = 5000, seed_probes = NULL) {
  ids <- if (!is.null(seed_probes)) seed_probes else {
    n <- min(n_probes, nrow(mat))
    rownames(mat)[order(stats::runif(nrow(mat)))[seq_len(n)]]
  }
  sub <- mat[intersect(ids, rownames(mat)), , drop = FALSE]
  df <- as.data.frame(as.table(sub))
  colnames(df) <- c("probe", "sample", "beta")
  df[!is.na(df$beta), , drop = FALSE]
}

methyl_mds_scores <- function(mat, n_features = 5000, k = 2) {
  m <- stats::na.omit(mat)
  if (nrow(m) < 10 || ncol(m) < 4) {
    return(list(ok = FALSE, reason = "Not enough complete-case probes/samples for MDS."))
  }
  v <- methyl_row_vars(m)
  top <- order(v, decreasing = TRUE)[seq_len(min(n_features, nrow(m)))]
  d <- stats::dist(t(m[top, , drop = FALSE]))
  k_used <- min(k, ncol(m) - 1)
  fit <- tryCatch(stats::cmdscale(d, k = k_used), error = function(e) NULL)
  if (is.null(fit)) return(list(ok = FALSE, reason = "MDS failed to converge on this distance matrix."))
  colnames(fit) <- paste0("Dim", seq_len(ncol(fit)))
  list(ok = TRUE, scores = fit)
}

methyl_mean_sd_table <- function(mat, n_probes = 20000) {
  ids <- rownames(mat)
  if (length(ids) > n_probes) ids <- ids[order(stats::runif(length(ids)))[seq_len(n_probes)]]
  m <- mat[ids, , drop = FALSE]
  means <- rowMeans(m, na.rm = TRUE)
  sds <- sqrt(pmax(methyl_row_vars(m), 0))
  data.frame(probe = ids, mean = means, sd = sds, row.names = NULL)
}

methyl_control_probe_matrix <- function(rg_set) {
  if (is.null(rg_set) || !requireNamespace("minfi", quietly = TRUE)) {
    return(list(ok = FALSE, reason = "Control-probe intensities require raw IDAT input."))
  }
  ctrl <- tryCatch(minfi::getProbeInfo(rg_set, type = "Control"), error = function(e) e)
  if (inherits(ctrl, "error") || is.null(ctrl) || nrow(ctrl) == 0) {
    return(list(ok = FALSE, reason = "Could not extract control-probe information from this IDAT upload."))
  }
  green <- tryCatch(minfi::getGreen(rg_set), error = function(e) NULL)
  red <- tryCatch(minfi::getRed(rg_set), error = function(e) NULL)
  if (is.null(green) || is.null(red)) {
    return(list(ok = FALSE, reason = "Could not read raw green/red channel intensities for control probes."))
  }
  addr <- intersect(ctrl$Address, intersect(rownames(green), rownames(red)))
  if (length(addr) == 0) return(list(ok = FALSE, reason = "No control-probe addresses matched this array's raw intensity data."))
  ctrl <- ctrl[match(addr, ctrl$Address), , drop = FALSE]
  intensity <- (green[addr, , drop = FALSE] + red[addr, , drop = FALSE]) / 2
  rownames(intensity) <- sprintf("%s: %s", ctrl$Type, ctrl$ExtendedType)
  list(ok = TRUE, mat = log2(pmax(intensity, 1)), types = ctrl$Type)
}

methyl_outlier_diagnostic_table <- function(sample_qc, pca_detail, mahal_detail = NULL) {
  df <- data.frame(sample = sample_qc$sample, stringsAsFactors = FALSE)
  if ("call_rate" %in% colnames(sample_qc)) df$call_rate <- sample_qc$call_rate
  if (isTRUE(pca_detail$ok)) {
    df$pca_distance <- pca_detail$distance[df$sample]
  }
  if (!is.null(mahal_detail) && isTRUE(mahal_detail$ok)) {
    df$mahalanobis_distance2 <- mahal_detail$distance2[df$sample]
  }
  df
}

METHYL_BATCH_COL_PATTERN <- "batch|chip|plate|slide|sentrix|array_id|scan_date|^run$"

METHYL_SEX_COL_CANDIDATES <- c("sex", "Sex", "SEX", "gender", "Gender")

methyl_batch_columns <- function(sheet) {
  if (is.null(sheet)) return(character(0))
  grep(METHYL_BATCH_COL_PATTERN, colnames(sheet), ignore.case = TRUE, value = TRUE)
}

methyl_guess_batch_column <- function(sheet) {
  hit <- methyl_batch_columns(sheet)
  if (length(hit) == 0) NULL else hit[1]
}

methyl_batch_correct_combat <- function(mat, batch, input_scale = "beta") {
  if (!requireNamespace("sva", quietly = TRUE)) {
    return(list(ok = FALSE, reason = "The sva package (ComBat) is not available in this deployment."))
  }
  batch <- as.character(batch)
  tbl <- table(batch)
  if (length(tbl) < 2) {
    return(list(ok = FALSE, reason = "Batch correction needs at least 2 distinct batch levels; only 1 was found."))
  }
  if (any(tbl < 2)) {
    return(list(ok = FALSE, reason = sprintf(
      "Batch level(s) with fewer than 2 samples cannot be adjusted by ComBat: %s.",
      paste(names(tbl[tbl < 2]), collapse = ", "))))
  }
  m <- if (identical(input_scale, "m")) mat else methyl_beta_to_mvalue(mat)
  keep <- stats::complete.cases(m)
  corrected_m <- tryCatch(
    sva::ComBat(dat = m[keep, , drop = FALSE], batch = batch, par.prior = TRUE, mean.only = FALSE),
    error = function(e) NULL
  )
  if (is.null(corrected_m)) {
    return(list(ok = FALSE, reason = "ComBat failed to converge on this matrix/batch split."))
  }
  corrected <- m
  corrected[keep, ] <- if (identical(input_scale, "m")) corrected_m else 2^corrected_m / (1 + 2^corrected_m)
  corrected[!keep, ] <- mat[!keep, , drop = FALSE]
  list(ok = TRUE, corrected = corrected, batch = batch, n_batches = length(tbl))
}

methyl_batch_correct_ruvm <- function(mat, rg_set, group, k = 1, input_scale = "beta") {
  if (is.null(rg_set)) {
    return(list(ok = FALSE, reason = "RUVm requires raw IDAT input (its internal negative-control probes aren't present in a beta/M-value matrix upload or the preloaded dataset's bundled matrix)."))
  }
  if (!requireNamespace("missMethyl", quietly = TRUE) || !requireNamespace("minfi", quietly = TRUE)) {
    return(list(ok = FALSE, reason = "The missMethyl package (RUVm) is not available in this deployment."))
  }
  group <- as.character(group)
  if (length(unique(stats::na.omit(group))) < 2) {
    return(list(ok = FALSE, reason = "RUVm needs at least 2 levels in the selected \"factor of interest\" column to know what biological signal to protect."))
  }
  incs <- tryCatch(missMethyl::getINCs(rg_set), error = function(e) e)
  if (inherits(incs, "error") || is.null(incs) || nrow(incs) == 0) {
    return(list(ok = FALSE, reason = "Could not extract internal negative-control probes from this IDAT upload (missMethyl::getINCs)."))
  }
  common_samples <- intersect(colnames(mat), intersect(colnames(incs), names(group)[!is.na(group)]))
  if (length(common_samples) < 4) {
    return(list(ok = FALSE, reason = "Fewer than 4 samples have both a matrix column, a control-probe column, and a non-missing factor-of-interest value."))
  }
  mat_sub <- mat[, common_samples, drop = FALSE]
  m <- if (identical(input_scale, "m")) mat_sub else methyl_beta_to_mvalue(mat_sub)
  keep <- stats::complete.cases(m)
  grp <- group[common_samples]
  design <- stats::model.matrix(~grp)
  mc <- rbind(m[keep, , drop = FALSE], as.matrix(incs[, common_samples, drop = FALSE]))
  ctl <- c(rep(FALSE, sum(keep)), rep(TRUE, nrow(incs)))
  fit <- tryCatch(missMethyl::RUVfit(Y = t(mc), X = design[, 2, drop = FALSE], ctl = ctl, k = k, method = "ruv4"), error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(ok = FALSE, reason = paste("RUVfit failed to converge on this matrix/design:", conditionMessage(fit))))
  }
  adj <- tryCatch(missMethyl::RUVadj(Y = t(mc), fit = fit), error = function(e) e)
  if (inherits(adj, "error")) {
    return(list(ok = FALSE, reason = paste("RUVadj failed:", conditionMessage(adj))))
  }
  corrected_m <- t(adj$newY)[seq_len(sum(keep)), , drop = FALSE]
  corrected <- m
  corrected[keep, ] <- if (identical(input_scale, "m")) corrected_m else 2^corrected_m / (1 + 2^corrected_m)
  corrected[!keep, ] <- mat_sub[!keep, , drop = FALSE]
  list(ok = TRUE, corrected = corrected, group = grp, k = k, n_controls = nrow(incs))
}

methyl_plot_cascade <- function(cascade_df) {
  df <- cascade_df
  df$step <- factor(df$step, levels = df$step)
  ggplot(df, aes(x = step, y = retained)) +
    geom_col(fill = ARTHOMIX_COLORS$blue, width = 0.6) +
    geom_text(aes(label = format(retained, big.mark = ",")), vjust = -0.4, size = 3.4, color = ARTHOMIX_COLORS$ink) +
    labs(x = NULL, y = "Probes remaining") + theme_arthomix() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

methyl_plot_detp_heatmap <- function(long_df) {
  ggplot(long_df, aes(x = sample, y = probe, fill = detection_p)) +
    geom_tile() +
    scale_fill_gradient(low = ARTHOMIX_COLORS$blue, high = ARTHOMIX_STATUS$critical) +
    labs(x = NULL, y = NULL, fill = "Detection p") + theme_arthomix() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6), axis.text.y = element_blank(), axis.ticks = element_blank())
}

methyl_plot_beadcount_dist <- function(vals, threshold = 3) {
  ggplot(data.frame(beadcount = vals[is.finite(vals)]), aes(x = beadcount)) +
    geom_histogram(bins = 40, fill = ARTHOMIX_COLORS$aqua, color = NA) +
    geom_vline(xintercept = threshold, color = ARTHOMIX_STATUS$critical, linetype = "dashed") +
    labs(x = "Bead count", y = "Probe x sample count") + theme_arthomix()
}

.methyl_stage_fill <- c("Before filtering" = ARTHOMIX_COLORS$blue, "After filtering" = ARTHOMIX_COLORS$orange,
                        "Before normalization" = ARTHOMIX_COLORS$blue, "After normalization" = ARTHOMIX_COLORS$orange)

methyl_plot_density <- function(density_df, x_label = "Beta value") {
  ggplot(density_df, aes(x = beta, color = stage)) +
    geom_density(linewidth = 0.6) +
    scale_color_manual(values = .methyl_stage_fill) +
    labs(x = x_label, y = "Density", color = NULL) + theme_arthomix()
}

methyl_plot_boxplot <- function(long_df) {
  ggplot(long_df, aes(x = sample, y = beta, fill = stage)) +
    geom_boxplot(outlier.size = 0.4, linewidth = 0.3) +
    scale_fill_manual(values = .methyl_stage_fill) +
    labs(x = NULL, y = "Beta value", fill = NULL) + theme_arthomix() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
}

methyl_plot_violin <- function(long_df) {
  ggplot(long_df, aes(x = sample, y = beta, fill = stage)) +
    geom_violin(scale = "width", linewidth = 0.2, alpha = 0.85) +
    scale_fill_manual(values = .methyl_stage_fill) +
    labs(x = NULL, y = "Beta value", fill = NULL) + theme_arthomix() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
}

methyl_plot_mean_sd <- function(mean_sd_df) {
  ggplot(mean_sd_df, aes(x = mean, y = sd)) +
    geom_point(color = ARTHOMIX_COLORS$blue, alpha = 0.25, size = 0.8) +
    geom_smooth(color = ARTHOMIX_STATUS$critical, se = FALSE, linewidth = 0.7, method = "loess", formula = y ~ x) +
    labs(x = "Mean", y = "Standard deviation") + theme_arthomix()
}

methyl_plot_scatter2d <- function(df, x_lab, y_lab, color_lab = NULL, palette = NULL) {
  gg <- ggplot(df, aes(x = x, y = y, color = color, text = text)) +
    geom_point(size = 2.5, alpha = 0.85) +
    labs(x = x_lab, y = y_lab, color = color_lab) + theme_arthomix()
  if (!is.null(palette)) gg <- gg + scale_color_manual(values = palette)
  gg
}

methyl_plot_corr_heatmap <- function(cor_matrix) {
  ord <- stats::hclust(stats::as.dist(1 - cor_matrix), method = "average")$order
  cm <- cor_matrix[ord, ord]
  df <- as.data.frame(as.table(cm)); colnames(df) <- c("sample_x", "sample_y", "correlation")
  df$sample_x <- factor(df$sample_x, levels = rownames(cm))
  df$sample_y <- factor(df$sample_y, levels = rev(rownames(cm)))
  rng <- range(df$correlation, na.rm = TRUE)
  ggplot(df, aes(x = sample_x, y = sample_y, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = ARTHOMIX_COLORS$red, mid = "white", high = ARTHOMIX_COLORS$blue, midpoint = mean(rng), limits = rng) +
    labs(x = NULL, y = NULL, fill = "r") + theme_arthomix() +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
}

methyl_plot_outlier_diagnostic <- function(diag_df, metric_col, y_lab) {
  df <- diag_df[order(-diag_df[[metric_col]]), , drop = FALSE]
  df$sample <- factor(df$sample, levels = df$sample)
  ggplot(df, aes(x = sample, y = .data[[metric_col]])) +
    geom_col(fill = ARTHOMIX_COLORS$blue, width = 0.6) +
    labs(x = NULL, y = y_lab) + theme_arthomix() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
}

methyl_check_beta_range <- function(mat, input_scale) {
  rng <- range(mat, na.rm = TRUE)
  out_of_range <- identical(input_scale, "beta") && is.finite(rng[1]) && is.finite(rng[2]) &&
    (rng[1] < -0.05 || rng[2] > 1.05)
  list(range = rng, out_of_range = out_of_range)
}

methyl_qc_status_badge <- function(overview) {
  reasons <- character(0)
  status <- "pass"
  cr <- overview$median_call_rate
  if (!is.na(cr) && cr < 0.90) {
    status <- "fail"
    reasons <- c(reasons, sprintf("Median sample call rate is low (%.1f%%).", cr * 100))
  }
  miss <- overview$overall_missing_pct
  if (!is.na(miss) && miss > 10 && status == "pass") {
    status <- "warning"
    reasons <- c(reasons, sprintf("Overall missingness across the loaded matrix is somewhat high (%.1f%%).", miss))
  }
  if (isTRUE(overview$range_check$out_of_range)) {
    if (status == "pass") status <- "warning"
    reasons <- c(reasons, sprintf(
      "Matrix is labeled as beta values (0-1) but observed values range [%.3g, %.3g] - check the input scale selected on the Dataset tab.",
      overview$range_check$range[1], overview$range_check$range[2]))
  }
  if (length(reasons) == 0) reasons <- "No issues detected in this basic pass - run the other QC tabs below for deeper checks."
  list(status = status, reasons = reasons)
}

methyl_qc_r_code <- function(settings) {
  s <- settings
  lines <- c(
    "## Equivalent Bioconductor QC pipeline for the settings used in this run",
    "library(minfi); library(ChAMP); library(wateRmelon)",
    "",
    "## Probe filtering",
    if (isTRUE(s$f_detp)) sprintf("detP <- detectionP(rgSet); failed <- detP > %s", format(s$detp_thresh %||% 0.01)),
    if (isTRUE(s$f_beadcount)) sprintf("beadcount <- beadcount(rgSet); low_beads <- beadcount < %s", format(s$beadcount_thresh %||% 3)),
    if (isTRUE(s$f_snp)) "mset <- champ.filter(beta = beta, pd = pd, filterSNPHit = TRUE)",
    if (!identical(s$sexchr_mode %||% "keep", "keep")) sprintf("mset <- champ.filter(mset, filterXY = TRUE) # mode: %s", s$sexchr_mode %||% "remove_xy"),
    if (isTRUE(s$f_noncpg)) "mset <- champ.filter(mset, filterNoCG = TRUE)",
    if (isTRUE(s$f_maf)) sprintf("beta <- beta[!(rownames(beta) %%in%% names(maf)[maf > %s]), ] # uploaded MAF table", format(s$maf_max %||% 0.05)),
    if (isTRUE(s$f_missing)) sprintf("beta <- beta[rowMeans(is.na(beta)) <= %s, ]", format(s$missing_max %||% 0)),
    if (isTRUE(s$f_variance)) sprintf("beta <- beta[apply(beta, 1, var, na.rm = TRUE) >= %s, ]", format(s$variance_min %||% 0)),
    if (isTRUE(s$f_sd)) sprintf("beta <- beta[apply(beta, 1, sd, na.rm = TRUE) >= %s, ]", format(s$sd_min %||% 0)),
    "",
    "## Sample QC, Sex QC, Outlier QC, Batch QC, and Normalization each have their",
    "## own independent settings/code, shown in their own tab once run - this panel",
    "## only reflects the Probe QC filters above."
  )
  paste(Filter(Negate(is.null), lines), collapse = "\n")
}

methyl_qc_summary_table <- function(overview = NULL, sample_qc = NULL, probe_qc = NULL,
                                     sex_qc = NULL, outlier_qc = NULL, batch_qc = NULL) {
  rows <- list()
  add <- function(metric, value) rows[[length(rows) + 1]] <<- data.frame(metric = metric, value = as.character(value), stringsAsFactors = FALSE)

  if (!is.null(overview)) {
    add("overview_samples", overview$n_samples)
    add("overview_cpgs", overview$n_probes)
    add("overview_median_call_rate", sprintf("%.4f", overview$median_call_rate))
    add("overview_missingness_pct", sprintf("%.2f", overview$overall_missing_pct))
  } else add("overview", "not run")

  if (!is.null(sample_qc)) {
    add("sample_qc_n_samples", ncol(sample_qc$mat))
    add("sample_qc_below_call_rate_min", sum(sample_qc$sample_qc$call_rate_flag, na.rm = TRUE))
  } else add("sample_qc", "not run")

  if (!is.null(probe_qc)) {
    add("probe_qc_cpgs_in", nrow(probe_qc$mat))
    add("probe_qc_cpgs_kept", nrow(probe_qc$filtered))
    add("probe_qc_cpgs_removed", nrow(probe_qc$mat) - nrow(probe_qc$filtered))
    add("probe_qc_removed_pct", sprintf("%.2f", 100 * (nrow(probe_qc$mat) - nrow(probe_qc$filtered)) / nrow(probe_qc$mat)))
  } else add("probe_qc", "not run")

  if (!is.null(sex_qc)) {
    if (isTRUE(sex_qc$sex$ok)) {
      tbl <- table(sex_qc$sex$detail$predicted_sex)
      add("sex_qc_predicted_distribution", paste(sprintf("%s:%d", names(tbl), as.integer(tbl)), collapse = ", "))
      add("sex_qc_mismatches", if (is.na(sex_qc$sex$n_mismatch)) "n/a" else sex_qc$sex$n_mismatch)
    } else add("sex_qc", sprintf("ran, no result: %s", sex_qc$sex$reason))
  } else add("sex_qc", "not run")

  if (!is.null(outlier_qc)) {
    add("outlier_qc_flagged_samples", sum(outlier_qc$outlier_scores$outlier_score > 0, na.rm = TRUE))
  } else add("outlier_qc", "not run")

  if (!is.null(batch_qc)) {
    if (isTRUE(batch_qc$out$ok)) {
      add("batch_qc_method", batch_qc$method)
      add("batch_qc_n_groups", length(unique(stats::na.omit(as.character(batch_qc$batch)))))
    } else add("batch_qc", sprintf("ran, no result: %s", batch_qc$out$reason))
  } else add("batch_qc", "not run")

  do.call(rbind, rows)
}

methyl_rmd_safe_text <- function(x) {
  x <- as.character(x %||% "")
  x[is.na(x)] <- ""
  x <- gsub("`", "'", x, fixed = TRUE)
  x <- gsub("[\r\n]+", " ", x)
  x
}

methyl_qc_report_plots <- function(methyl_dataset, probe_qc = NULL, outlier_qc = NULL) {
  plots <- list()
  skipped <- character(0)

  if (!is.null(probe_qc)) {
    plots[["probe_filtering_cascade"]] <- methyl_plot_cascade(probe_qc$cascade)

    x_lab <- if (identical(methyl_dataset$input_scale, "beta")) "Beta value" else "M-value"
    probe_ids <- rownames(probe_qc$mat)[order(stats::runif(nrow(probe_qc$mat)))[seq_len(min(5000, nrow(probe_qc$mat)))]]
    before <- methyl_beta_density_sample(probe_qc$mat, seed_probes = probe_ids); before$stage <- "Before filtering"
    after_ids <- intersect(probe_ids, rownames(probe_qc$filtered))
    density_df <- before
    if (length(after_ids) > 0) {
      after <- methyl_beta_density_sample(probe_qc$filtered, seed_probes = after_ids); after$stage <- "After filtering"
      density_df <- rbind(before, after)
    }
    plots[["beta_density"]] <- methyl_plot_density(density_df, x_lab)
    plots[["beta_boxplot"]] <- methyl_plot_boxplot(density_df)
    plots[["beta_violin"]] <- methyl_plot_violin(density_df)

    cr <- methyl_sample_correlation(probe_qc$filtered)
    if (isTRUE(cr$ok)) plots[["sample_correlation_heatmap"]] <- methyl_plot_corr_heatmap(cr$cor) else skipped <- c(skipped, sprintf("Sample correlation heatmap: %s", cr$reason))

    pca <- methyl_pca_scores(probe_qc$filtered)
    if (isTRUE(pca$ok)) {
      df <- data.frame(x = pca$scores[, 1], y = pca$scores[, 2], color = "Sample", text = rownames(pca$scores))
      plots[["pca_2d"]] <- methyl_plot_scatter2d(df, "PC1", "PC2", NULL, palette = c("Sample" = ARTHOMIX_COLORS$blue))
    } else skipped <- c(skipped, sprintf("PCA: %s", pca$reason))

    mds <- methyl_mds_scores(probe_qc$filtered)
    if (isTRUE(mds$ok)) {
      df <- data.frame(x = mds$scores[, 1], y = mds$scores[, 2], color = "Sample", text = rownames(mds$scores))
      plots[["mds"]] <- methyl_plot_scatter2d(df, "Dim1", "Dim2", NULL, palette = c("Sample" = ARTHOMIX_COLORS$blue))
    } else skipped <- c(skipped, sprintf("MDS: %s", mds$reason))

    plots[["mean_sd"]] <- methyl_plot_mean_sd(methyl_mean_sd_table(probe_qc$filtered))
  } else skipped <- c(skipped, "Probe-QC-derived figures (cascade, density, boxplot, violin, correlation heatmap, PCA, MDS, mean-SD) were skipped - run Probe QC first.")

  if (!is.null(outlier_qc) && isTRUE(outlier_qc$pca_detail$ok)) {
    diag <- methyl_outlier_diagnostic_table(outlier_qc$sample_qc, outlier_qc$pca_detail, outlier_qc$mahal_detail)
    if ("pca_distance" %in% colnames(diag)) plots[["outlier_diagnostic"]] <- methyl_plot_outlier_diagnostic(diag, "pca_distance", "PCA distance from centroid")
  } else skipped <- c(skipped, "Outlier diagnostic plot was skipped - run Outlier Detection first.")

  list(plots = plots, skipped = skipped)
}

methyl_qc_report_html <- function(methyl_dataset, summary_df, plots, skipped, subtitle = NULL) {
  if (!requireNamespace("htmltools", quietly = TRUE) || !requireNamespace("base64enc", quietly = TRUE)) {
    return(list(ok = FALSE, reason = "The htmltools/base64enc packages are required to build the HTML report."))
  }
  fig_tags <- lapply(names(plots), function(nm) {
    f <- tempfile(fileext = ".png")
    ggplot2::ggsave(f, plots[[nm]], width = 8, height = 5, dpi = 110)
    uri <- base64enc::dataURI(file = f, mime = "image/png")
    unlink(f)
    htmltools::tags$div(
      htmltools::tags$h3(gsub("_", " ", nm)),
      htmltools::tags$img(src = uri, style = "max-width:100%;")
    )
  })
  summary_rows <- lapply(seq_len(nrow(summary_df)), function(i)
    htmltools::tags$tr(htmltools::tags$td(summary_df$metric[i]), htmltools::tags$td(as.character(summary_df$value[i]))))
  doc <- htmltools::tags$html(
    htmltools::tags$head(htmltools::tags$title("Methylomics QC report"),
      htmltools::tags$style("body{font-family:sans-serif;margin:24px;} table{border-collapse:collapse;} td{border:1px solid #ddd;padding:4px 10px;} h3{margin-top:32px;}")),
    htmltools::tags$body(
      htmltools::tags$h1("Methylomics quality-control report"),
      htmltools::tags$p(subtitle %||% sprintf("Source: %s.", methyl_dataset$source %||% "n/a")),
      htmltools::tags$h2("Summary"),
      htmltools::tags$table(summary_rows),
      htmltools::tags$h2("Figures"),
      fig_tags,
      if (length(skipped) > 0) htmltools::tagList(
        htmltools::tags$h2("Skipped / not yet run"),
        htmltools::tags$ul(lapply(skipped, htmltools::tags$li))
      )
    )
  )
  out <- tempfile(fileext = ".html")
  htmltools::save_html(doc, file = out)
  list(ok = TRUE, path = out)
}

methyl_qc_report_zip <- function(plots) {
  if (length(plots) == 0) return(list(ok = FALSE, reason = "No figures were available to export for this run."))
  dir <- tempfile("methyl_qc_figs_"); dir.create(dir)
  for (nm in names(plots)) {
    ggplot2::ggsave(file.path(dir, paste0(nm, ".png")), plots[[nm]], width = 8, height = 5, dpi = 150)
  }
  zip_path <- tempfile(fileext = ".zip")
  old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
  setwd(dir)
  status <- tryCatch(utils::zip(zip_path, list.files(dir, pattern = "\\.png$"), flags = "-q"), error = function(e) -1L)
  setwd(old_wd)
  if (!file.exists(zip_path) || (is.numeric(status) && status != 0)) {
    return(list(ok = FALSE, reason = "Could not create a ZIP archive - the `zip` command-line tool may not be available in this deployment."))
  }
  list(ok = TRUE, path = zip_path)
}
