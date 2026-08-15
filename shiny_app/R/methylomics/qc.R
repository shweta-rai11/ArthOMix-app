## R/methylomics/qc.R
## Probe- and sample-level QC for the Methylomics module (shiny_app/R/mod_methyl_qc.R).
## Every probe filter below returns list(keep = <logical over rownames(mat)>,
## note = <one-line summary>) rather than filtering in place, so the caller
## can show how many probes each individual filter would remove before
## combining them into one final keep vector. Filters that need raw
## intensities (detection p-value, beadcount) or manifest annotation
## (SNP/sex-chromosome) return keep = all-TRUE with an explanatory note
## when that input isn't available, rather than silently doing nothing or
## guessing.

## ---- Probe filters --------------------------------------------------------

## Row-wise variance, vectorized via matrixStats (a compiled, single-pass
## computation) when available, falling back to base R's apply() otherwise.
## Every top-variance-probe selection in this file (probe filtering, PCA,
## hierarchical clustering, Mahalanobis distance, the correlation heatmap)
## goes through this - at 400k+ probes, apply(m, 1, var) alone takes long
## enough (several seconds, paid on every recompute) to make the live QC
## tool feel broken; matrixStats::rowVars() is the same computation roughly
## two orders of magnitude faster, with no change in result.
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

## Standard-deviation companion to methyl_filter_variance() above - same
## computation (methyl_row_vars(), sqrt'd), offered as its own independent
## optional filter since "minimum SD" is the more familiar unit for most
## methylation QC workflows even though it's mathematically the same cut.
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

## Only meaningful for an actual Illumina array upload - see
## METHYL_ARRAY_TYPES_ILLUMINA in annotation.R; the caller decides whether
## to offer this filter for the selected array type.
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

## `mode`: "remove_xy" (drop both chrX and chrY, the original behavior),
## "remove_y_only" (drop chrY only - keeps chrX, useful for studies that
## want to retain X-linked signal but still avoid Y's near-binary
## presence/absence-by-sex artifact), or "keep" (no-op, included so callers
## can route every sex-chromosome choice through this one function instead
## of branching around it).
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

## No published cross-reactive probe blacklist (Chen et al. 2013 for 450K;
## Pidsley et al. 2016 / McCartney et al. 2016 for EPIC) is bundled in this
## deployment - fabricating one would violate this project's own
## evidence-based-methods requirement. Instead this accepts a user-supplied
## exclusion list (methyl_parse_probe_list()), which can be exactly one of
## those published lists downloaded and uploaded by the user.
methyl_filter_cross_reactive <- function(mat, exclusion_ids = NULL) {
  if (is.null(exclusion_ids) || length(exclusion_ids) == 0) {
    return(list(keep = rep(TRUE, nrow(mat)), note =
      "No cross-reactive probe list is bundled in this deployment - upload a probe-exclusion list (one probe ID per line, e.g. a published Chen et al. 2013 / Pidsley et al. 2016 / McCartney et al. 2016 list) to enable this filter."))
  }
  keep <- !(rownames(mat) %in% exclusion_ids)
  list(keep = keep, note = sprintf("%d probe(s) removed via the uploaded exclusion list.", sum(!keep)))
}

## Optional probe_id,maf exclusion list - same shape/parsing convention as
## methyl_parse_probe_list() (shiny_app/R/methylomics/parse_upload.R), but
## also reads a numeric MAF column (by name if present - "maf"/"MAF"/"af"/
## "AF" - else the second field positionally) so methyl_filter_maf() below
## can threshold on it. No population-MAF table is bundled in this
## deployment (see methyl_filter_cross_reactive()'s comment above for why:
## this app's own manifest annotation has no allele-frequency column), so -
## exactly like the cross-reactive filter - this only activates against a
## user-uploaded list (e.g. a downloaded 1000 Genomes/gnomAD MAF-by-probe
## table for this array type).
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

## Removes probes whose uploaded MAF exceeds `max_maf` - i.e. keeps probes
## that are either below the threshold or absent from the uploaded table
## (absence isn't evidence of a common SNP, so it's not penalized).
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

## Removes a probe if its bead count falls below `threshold` in ANY sample -
## stricter than ChAMP's own default (champ.filter(beadCutoff=0.05) removes
## a probe only when it fails in >5% of samples), a deliberate choice here
## rather than an attempt to reproduce ChAMP's exact rule.
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

## ---- Sample-level QC -------------------------------------------------------

methyl_sample_call_rate <- function(mat) {
  1 - colMeans(is.na(mat))
}

## Per-sample failed-probe percentage: the fraction of probes with
## detection p-value > `threshold` in that sample - the sample-level
## transpose of methyl_filter_detection_p()'s per-probe "fails in at least
## one sample" logic above. Requires raw IDAT detection p-values, same
## degrade-with-a-reason pattern as bisulfite conversion/median intensity
## in idat_metrics.R (a matrix-only upload never has `detp`).
methyl_sample_failed_probe_pct <- function(mat, detp, threshold = 0.01) {
  if (is.null(detp)) {
    return(list(ok = FALSE, reason = "Failed-probe percentage requires raw IDAT input - not available for an uploaded beta/M-value matrix."))
  }
  common_probes <- intersect(rownames(mat), rownames(detp))
  common_samples <- intersect(colnames(mat), colnames(detp))
  if (length(common_probes) == 0 || length(common_samples) == 0) {
    return(list(ok = FALSE, reason = "No overlapping probes/samples between the detection p-value matrix and this run."))
  }
  ## NA (not all(colnames(mat)) - detp may not cover every stratum sample if
  ## the sample sheet's IDs don't perfectly match, same tolerance the
  ## per-probe detection-p filter already has via its own intersect().
  pct <- stats::setNames(rep(NA_real_, ncol(mat)), colnames(mat))
  pct[common_samples] <- colMeans(detp[common_probes, common_samples, drop = FALSE] > threshold, na.rm = TRUE) * 100
  list(ok = TRUE, pct = pct)
}

## Per-sample low-intensity flag, built on idat_metrics.R's
## methyl_median_intensity() (minfi::getQC()'s median log2 methylated/
## unmethylated intensity) - a sample is "low intensity" when the average
## of its two medians falls below `min_intensity`, the same signal minfi's
## own QC tutorial plots (mMed vs uMed) to catch failed arrays. IDAT-only,
## same degrade pattern as everything else that needs raw intensities.
methyl_sample_low_intensity <- function(median_int_result, min_intensity = 10) {
  if (!isTRUE(median_int_result$ok)) {
    return(list(ok = FALSE, reason = median_int_result$reason %||% "Median intensity requires raw IDAT input."))
  }
  d <- median_int_result$detail
  score <- (d$med_meth_log2 + d$med_unmeth_log2) / 2
  list(ok = TRUE, score = stats::setNames(score, d$sample), low = stats::setNames(score < min_intensity, d$sample))
}

methyl_sample_outliers_iqr <- function(x, k = 1.5) {
  q <- stats::quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  lo <- q[1] - k * iqr; hi <- q[2] + k * iqr
  x < lo | x > hi
}

## Correlation-based outlier flag: a sample whose mean pairwise correlation
## with every other sample is more than `k` MADs below the cohort median is
## flagged - a sample that simply doesn't look like the rest of the cohort,
## independent of any PCA/clustering projection. Same family of idea as
## WGCNA's own sample-network "connectivity Z-score" outlier check (Oldham/
## Horvath), but median/MAD here rather than WGCNA's mean/SD - a deliberate
## choice, not an attempt to reproduce WGCNA's exact formula: MAD is itself
## robust to the very low-correlation outliers this is trying to detect,
## where an SD computed across the same samples would already be inflated
## by them. Built on the same top-variance-probe correlation matrix
## methyl_sample_correlation() (below) already computes for the
## Visualizations tab's heatmap, so this is a re-derivation of that same
## signal into a per-sample flag rather than a second, differently-computed
## correlation.
methyl_sample_outliers_correlation <- function(mat, n_features = 5000, k = 3) {
  cr <- methyl_sample_correlation(mat, n_features = n_features)
  if (!isTRUE(cr$ok)) return(list(ok = FALSE, reason = cr$reason))
  diag(cr$cor) <- NA
  mean_cor <- rowMeans(cr$cor, na.rm = TRUE)
  med <- stats::median(mean_cor); mad_val <- stats::mad(mean_cor)
  lo <- med - k * mad_val
  list(ok = TRUE, mean_correlation = mean_cor, outlier = mean_cor < lo, threshold = lo)
}

## PCA-distance outlier flag: top-variance probes -> prcomp -> Euclidean
## distance from the centroid on the first two PCs, in units of that PC's
## own SD. A common, easy-to-explain ad hoc heuristic in array QC tooling
## (the same idea arrayQualityMetrics-style PCA plots visualize), but it is
## NOT a formal statistical test - it only uses 2 of the retained PCs and
## has no citable significance threshold behind the "3 SD" default the way
## a proper Hotelling's T^2 test would. Treat it as a quick visual-style
## flag, not a rigorous outlier test; methyl_sample_outliers_mahalanobis()
## below is the more principled multivariate alternative.
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

## Average-linkage hierarchical clustering on the same top-variance
## probes; any sample left in a singleton cluster at height_frac of the
## tree's max height is flagged.
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

## Mahalanobis distance computed on the top PCA components (not the raw
## probes - p >> n there, so the covariance matrix is singular), flagged
## against a chi-squared threshold at `alpha` - correct for a CLASSICAL
## Mahalanobis distance under multivariate normality (df = number of PCs
## used). Note this uses the classical (non-robust) mean/covariance
## (stats::cov()), which is itself sensitive to the very outliers it's
## trying to detect (masking/swamping) - a robust estimator (e.g. MCD, per
## Rousseeuw & Van Driessen 1999, as advocated for outlier detection by
## Filzmoser et al. 2005) would be more resistant, at the cost of needing
## more samples than this module can assume. Treat this as a reasonable,
## disclosed simplification, not a claim of being the most robust method
## available - methyl_sample_outliers_pca()/_correlation() above offer
## different, simpler heuristics to cross-check flagged samples against.
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

## Splits `y` (mean chrY probe beta per sample) into two groups and labels
## the higher-mean group "M" - unless `reported` (a same-length vector of
## reported sex, NA where unknown) lets each cluster instead be relabeled to
## whichever sex it agrees with most often. Uses k-means (nstart=25) rather
## than a single global-median threshold whenever there are enough samples
## for it to be meaningful: a median split misclassifies the *minority*
## cluster's mislabeled tail whenever the cohort is imbalanced enough that
## the global median falls inside the majority cluster instead of the gap
## between the two clusters - exactly what happens on this app's own
## reference cohort (492 female / 197 male: a median split alone
## misclassifies 147 of 492 females as male, despite chrY beta being
## cleanly bimodal with zero overlap between the true clusters). See
## script01_dataload_QC/01_loadandQC.R's own sex_km/sex_pred columns,
## which hit this exact failure mode and already fix it this same way -
## this mirrors that, so the live tool no longer disagrees with the
## reference pipeline's methodology on an imbalanced upload.
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

## Appends reported_sex/sex_mismatch columns to a sex-check detail table
## when a reported-sex vector (named by sample, from the sample sheet) is
## available, so the caller gets the same "concordance vs metadata" flag
## the preloaded default analysis already shows, instead of just a bare
## prediction with nothing to compare it against.
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

## Estimated sex from the methylation profile. When a raw RGChannelSet is
## available, uses minfi::getSex() (copy-number based on median chrX/chrY
## total-intensity log2 ratios, not beta - the standard method most
## published pipelines use). Otherwise falls back to a chrY-beta clustering
## heuristic via methyl_cluster_sex() - only mean chrY beta actually drives
## the cluster assignment (mean chrX beta is computed alongside it purely
## for the detail table/X-vs-Y scatter plot, not for the prediction itself)
## - clearly weaker than the raw-intensity method and labeled as such.
## `reported_sex`, if given
## (named by sample, e.g. from the sample sheet's sex column), resolves
## which methylation cluster is "male" by majority concordance instead of
## assuming higher chrY beta = male, and adds a sex_mismatch flag to the
## result - this is what mod_methyl_qc.R passes in for a live analysis.
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

## ---- Sample subgroup (sex/group) filtering --------------------------------

## Resolves sample-sheet rows to matrix column IDs: matches on a
## sample/Sample/sample_id/Sample_ID column when present, else assumes the
## sheet's row order matches `all_ids`'s order (a sheet exported alongside
## its matrix, or e.g. minfi's own Sample_Name-based convention, not in
## that column-name list). Deliberately NOT rownames(sheet) for that
## fallback: a freshly data.table::fread()-parsed data.frame gets sequential
## integer rownames ("1","2",...), which never equal a real sample ID, so
## using them silently turned every "row order" case into zero matches
## instead of the row-order match the caller actually wants. Only falls
## back to rownames(sheet) when the row count doesn't even match
## length(all_ids), since there's no safe way to align them positionally
## at that point (this preserves the old, already-broken behavior only for
## a case that was never going to align anyway).
methyl_sheet_sample_ids <- function(sheet, all_ids) {
  id_col <- intersect(c("sample", "Sample", "sample_id", "Sample_ID"), colnames(sheet))
  if (length(id_col) > 0) return(as.character(sheet[[id_col[1]]]))
  if (nrow(sheet) == length(all_ids)) return(all_ids)
  rownames(sheet)
}

## Resolves a "subgroup column" + "stratum" selection against
## methyl_dataset$sample_sheet into the matching subset of `mat`'s columns,
## so every downstream QC computation (call rate, outlier detection, PCA,
## sex check) runs on the chosen stratum instead of always the full cohort -
## this is what mod_methyl_qc.R's live filtering tool was missing, which is
## why picking "Female" vs "Male" there previously had no effect on any
## result. Returns the unfiltered matrix unchanged when no sheet, column,
## or the "All" stratum is selected, so this is a no-op unless a user has
## loaded a sheet and picked one specific stratum.
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

## Removes a manually-excluded sample set (e.g. from the Sample QC tab's
## inclusion/exclusion table, or a Sex QC "Exclude" action on a discordant
## sample) from an already-subgroup-filtered result, so one exclusion
## mechanism feeds every downstream computation (PCA, outlier detection, sex
## check, exports) instead of each tab tracking its own list.
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

## ---- Outlier scoring, probe-retention cascade, misc export helpers -------

## One row per sample, one column per outlier-detection method actually run
## (whichever of pca_outlier/hclust_outlier/mahalanobis_outlier/iqr_outlier
## columns are present in `sample_qc`), plus a summary `outlier_score` -
## how many of the selected methods flagged that sample. Used by the Outlier
## Detection tab so several independent methods collapse into one ranked
## table instead of one column each with no combined view.
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

## Sequential (cascade) probe retention: applies each named filter in
## `filters` (in the order given, same order the caller ran them in) as a
## running AND, and reports how many probes remain after each step - this is
## what the Probe QC tab's retention flowchart plots, distinct from the
## "probes removed by this filter alone" counts already shown in the filter
## summary table.
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

## Beta -> M-value (logit) transform for the "download M-value matrix" export
## - beta is clipped away from the exact 0/1 boundary first since log2(0/1)
## and log2(x/0) are undefined/infinite, matching how minfi::logit2()
## and lumi::beta2m() handle this same edge case.
methyl_beta_to_mvalue <- function(beta, eps = 1e-4) {
  b <- pmin(pmax(beta, eps), 1 - eps)
  log2(b / (1 - b))
}

## Top-variance-probe PCA scores, factored out of the outlier-detection
## functions above (which each independently redo this) for the new call
## sites that need PCA scores alone, without the outlier-flagging logic:
## the Visualizations tab's color-by-metadata PCA, and Batch Correction's
## before/after comparison.
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

## Sample-by-sample correlation matrix on the same top-variance-probe subset
## as the PCA/outlier functions above, for the Visualizations tab's
## correlation heatmap - computing this on all ~412k probes would be both
## slow and dominated by uninformative, near-constant probes.
methyl_sample_correlation <- function(mat, n_features = 5000) {
  m <- stats::na.omit(mat)
  if (nrow(m) < 10 || ncol(m) < 2) {
    return(list(ok = FALSE, reason = "Not enough complete-case probes/samples for a correlation heatmap."))
  }
  v <- methyl_row_vars(m)
  top <- order(v, decreasing = TRUE)[seq_len(min(n_features, nrow(m)))]
  list(ok = TRUE, cor = stats::cor(m[top, , drop = FALSE]))
}

## A random subset of probes' beta values in long format, one row per
## probe x sample, for the Visualizations tab's density plot - melting the
## full ~412k x n matrix would be both slow to plot and visually identical
## to a representative sample at this scale.
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

## Classical (metric) MDS via stats::cmdscale on the same top-variance-probe
## Euclidean-distance convention methyl_pca_scores() uses for PCA - a
## distinct sample-structure view from PCA (MDS optimizes pairwise
## distances directly rather than variance-maximizing axes), offered
## alongside it rather than instead of it since they can disagree on a
## noisy dataset.
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

## Per-probe mean vs. standard deviation, for the classic "mean-SD plot"
## (Huber et al. 2002) used to spot whether variance is confounded with
## mean intensity - matrixStats-accelerated the same way every other
## per-row summary in this file is (see methyl_row_vars()'s own comment).
methyl_mean_sd_table <- function(mat, n_probes = 20000) {
  ids <- rownames(mat)
  if (length(ids) > n_probes) ids <- ids[order(stats::runif(length(ids)))[seq_len(n_probes)]]
  m <- mat[ids, , drop = FALSE]
  means <- rowMeans(m, na.rm = TRUE)
  sds <- sqrt(pmax(methyl_row_vars(m), 0))
  data.frame(probe = ids, mean = means, sd = sds, row.names = NULL)
}

## Illumina control-probe (staining/extension/hybridization/target
## removal/specificity/non-polymorphic/negative, etc.) intensities for the
## control-probe heatmap - IDAT-only, since these are raw-intensity probes
## a beta/M-value matrix never carries.
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

## One row per sample combining every distance-based outlier signal already
## computed elsewhere in this file (PCA distance from centroid, Mahalanobis
## distance, call rate) into a single ranked diagnostic table/plot, instead
## of each metric only being visible in its own separate tab.
methyl_outlier_diagnostic_table <- function(sample_qc, pca_detail, mahal_detail = NULL) {
  df <- data.frame(sample = sample_qc$sample, stringsAsFactors = FALSE)
  ## `sample_qc` here is Outlier QC's own independent per-sample data frame
  ## (sample + flag columns only) - it doesn't carry call_rate the way
  ## Sample QC's own sample_qc data frame does, since the two tabs compute
  ## independently now. Included only when present so this function works
  ## for either caller's shape.
  if ("call_rate" %in% colnames(sample_qc)) df$call_rate <- sample_qc$call_rate
  if (isTRUE(pca_detail$ok)) {
    df$pca_distance <- pca_detail$distance[df$sample]
  }
  if (!is.null(mahal_detail) && isTRUE(mahal_detail$ok)) {
    df$mahalanobis_distance2 <- mahal_detail$distance2[df$sample]
  }
  df
}

## Best-guess batch/chip/plate column in a sample sheet, same convention as
## the sex-column guess already used for subgroup filtering
## (mod_methyl_qc.R's `default_col`). Deliberately returns NULL rather than
## guessing wrong when nothing matches - Batch Correction is disabled with an
## explanatory message in that case rather than silently offered against a
## column that isn't actually a batch identifier.
methyl_guess_batch_column <- function(sheet) {
  if (is.null(sheet)) return(NULL)
  hit <- grep("batch|chip|plate|slide|sentrix|array_id|scan_date|^run$",
              colnames(sheet), ignore.case = TRUE, value = TRUE)
  if (length(hit) == 0) NULL else hit[1]
}

## ComBat batch correction (sva::ComBat) on a beta-value matrix - logit'd to
## M-values first (ComBat assumes a roughly Gaussian, unbounded outcome;
## beta values are bounded to [0,1] and heavily bimodal near the ends, which
## is exactly what M-values correct for - Du et al. 2010 (BMC Bioinformatics,
## "Comparison of Beta-value and M-value methods") is the standard citation
## for M-values being the better-behaved, more homoscedastic scale for this
## kind of parametric/statistical modeling, which is why ComBat/SVA-style
## corrections on methylation arrays are conventionally run on M-values).
## Corrected M-values are converted back to beta for consistency with the
## rest of this module's beta-scale display. Requires >=2 samples in every
## batch level ComBat is asked to adjust, same requirement sva::ComBat()
## itself enforces.
methyl_batch_correct_combat <- function(mat, batch) {
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
  m <- methyl_beta_to_mvalue(mat)
  keep <- stats::complete.cases(m)
  corrected_m <- tryCatch(
    sva::ComBat(dat = m[keep, , drop = FALSE], batch = batch, par.prior = TRUE, mean.only = FALSE),
    error = function(e) NULL
  )
  if (is.null(corrected_m)) {
    return(list(ok = FALSE, reason = "ComBat failed to converge on this matrix/batch split."))
  }
  corrected_beta <- m
  corrected_beta[keep, ] <- 2^corrected_m / (1 + 2^corrected_m)
  corrected_beta[!keep, ] <- mat[!keep, , drop = FALSE]
  list(ok = TRUE, corrected = corrected_beta, batch = batch, n_batches = length(tbl))
}

## RUVm (Maksimovic et al. 2015) batch/unwanted-variation correction via
## missMethyl::RUVfit() - unlike ComBat, RUVm doesn't take an explicit batch
## label; it estimates unwanted variation empirically from the array's own
## internal negative-control probes (missMethyl::getINCs()), conditioned on
## a `group` factor of interest to protect (the biological signal RUVm must
## NOT remove). This is why it needs a raw RGChannelSet (control probes
## aren't in a beta/M-value matrix upload, and aren't in the preloaded
## dataset either - only the derived beta matrix is bundled there, see
## mod_methyl_dataset.R) rather than being usable everywhere ComBat is.
## `k` is the number of unwanted-variation factors to remove - missMethyl's
## own vignette uses k=1 as a starting point, adjustable here for the same
## reason every other method's own defaults are adjustable in this module.
methyl_batch_correct_ruvm <- function(mat, rg_set, group, k = 1) {
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
  m <- methyl_beta_to_mvalue(mat[, common_samples, drop = FALSE])
  keep <- stats::complete.cases(m)
  grp <- group[common_samples]
  design <- stats::model.matrix(~grp)
  mc <- rbind(m[keep, , drop = FALSE], as.matrix(incs[, common_samples, drop = FALSE]))
  ctl <- c(rep(FALSE, sum(keep)), rep(TRUE, nrow(incs)))
  ## method must be explicitly "ruv4" (or "ruv2") - RUVfit()'s own default
  ## method ("inv") ignores `k` entirely (it ignores the RUV-4/RUV-2 factor
  ## count and instead directly inverts the control-probe regression), which
  ## would make the "Unwanted-variation factors (k)" control above silently
  ## do nothing. ruv4 is the k-driven variant Maksimovic et al. 2015's own
  ## RUVm examples use.
  fit <- tryCatch(missMethyl::RUVfit(Y = t(mc), X = design[, 2, drop = FALSE], ctl = ctl, k = k, method = "ruv4"), error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(ok = FALSE, reason = paste("RUVfit failed to converge on this matrix/design:", conditionMessage(fit))))
  }
  adj <- tryCatch(missMethyl::RUVadj(Y = t(mc), fit = fit), error = function(e) e)
  if (inherits(adj, "error")) {
    return(list(ok = FALSE, reason = paste("RUVadj failed:", conditionMessage(adj))))
  }
  corrected_m <- t(adj$newY)[seq_len(sum(keep)), , drop = FALSE]
  corrected_beta <- m
  corrected_beta[keep, ] <- 2^corrected_m / (1 + 2^corrected_m)
  corrected_beta[!keep, ] <- mat[!keep, common_samples, drop = FALSE]
  list(ok = TRUE, corrected = corrected_beta, group = grp, k = k, n_controls = nrow(incs))
}

## ---- Shared ggplot builders -------------------------------------------
## One function per figure, called by both mod_methyl_qc.R's live
## renderPlot/renderPlotly outputs (wrapped in plotly::ggplotly() for the
## interactive ones) and methyl_qc_report_plots() below (for the HTML
## report + figures ZIP) - a figure's actual drawing logic lives in exactly
## one place instead of being re-implemented for the static/interactive/
## exported versions separately.

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

## Shared before/after color mapping - covers both the probe-filtering
## stage labels ("Before/After filtering") and the normalization-preview
## stage labels ("Before/After normalization"), so methyl_plot_density()/
## methyl_plot_boxplot()/methyl_plot_violin() all recognize whichever pair
## of labels the caller's long-format data frame actually uses.
.methyl_stage_fill <- c("Before filtering" = ARTHOMIX_COLORS$blue, "After filtering" = ARTHOMIX_COLORS$orange,
                        "Before normalization" = ARTHOMIX_COLORS$blue, "After normalization" = ARTHOMIX_COLORS$orange)

methyl_plot_density <- function(density_df, x_label = "Beta value") {
  ggplot(density_df, aes(x = beta, color = stage)) +
    geom_density(linewidth = 0.6) +
    scale_color_manual(values = .methyl_stage_fill) +
    labs(x = x_label, y = "Density", color = NULL) + theme_arthomix()
}

## Sample-wise beta distribution, one boxplot/violin per sample, colored by
## before/after stage - `long_df` has probe/sample/beta/stage columns, same
## shape methyl_beta_density_sample() already produces (stage added by the
## caller, same convention as methyl_plot_density() above).

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

## Generic 2D scatter for PCA/MDS - `df` needs x/y/color/text columns
## (caller names the axes via `x_lab`/`y_lab`); shared by PCA (2D), MDS, and
## the outlier/sex-check/batch-correction scatters that were already
## hand-rolling this same shape inline in mod_methyl_qc.R.
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

## Pass/Warning/Fail summary for the Overview tab's own self-contained
## "Run Overview QC" pass - deliberately built ONLY from `overview`'s own
## lightweight call-rate/missingness numbers, not from Sample/Probe/Sex/
## Outlier QC's results: those are independent, separately-run methods now
## (see mod_methyl_qc.R), and Overview must not silently depend on any of
## them having been run, or its own badge would change every time another
## tab is re-run without Overview itself being re-run - exactly the "one
## global result" anti-pattern this module was restructured to avoid.
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
  if (length(reasons) == 0) reasons <- "No issues detected in this basic pass - run the other QC tabs below for deeper checks."
  list(status = status, reasons = reasons)
}

## Plain-text minfi/ChAMP-equivalent R code reflecting the exact filter
## settings used in a Probe QC run (`settings` is probe_qc_result()$settings
## in mod_methyl_qc.R), for the Reports & Export tab's "Copy R code"
## reproducibility panel - a static template filled in from the actual
## settings, not a live-executed pipeline.
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

## ---- Live QC summary / report / export --------------------------------

## Single-row-per-metric summary assembled from whichever independent QC
## tabs (mod_methyl_qc.R) have actually been run this session - each of
## `overview`/`sample_qc`/`probe_qc`/`sex_qc`/`outlier_qc`/`batch_qc` is
## either that tab's own result object, or NULL if its Run button hasn't
## been clicked yet. A section that hasn't run shows as "not run" rather
## than being silently omitted, so the exported table always has the same
## shape and never implies a method ran when it didn't - the "Live QC
## Summary" export.
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

  ## sex_qc/batch_qc distinguish "not run" (NULL - tab's button never
  ## clicked) from "ran, but had nothing to report" (non-NULL with
  ## ok=FALSE - e.g. too few chrX/chrY probes, or no batch levels) - the
  ## two are different facts and collapsing them into one "not run" label
  ## would misreport a method that DID run as if it never had.
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

## Assembles every report/export figure as a plain ggplot object (never a
## plotly widget - htmltools::save_html()/ggplot2::ggsave() both need a
## static image), reusing the same methyl_plot_*() builders the live tabs
## use. Each figure group is gated on the independent tab that produces its
## underlying data (`probe_qc`/`outlier_qc`, both possibly NULL if that
## tab's Run button hasn't been clicked) - a figure that can't be built
## because its source analysis hasn't run yet is skipped with a message
## saying so, same as a figure skipped for a genuine data-availability
## reason (e.g. no detection p-values on a matrix upload).
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

## Self-contained HTML QC report - no rmarkdown/pandoc dependency (see
## qc.R's header note on this file's package-availability constraints):
## htmltools::save_html() with base64 data-URI <img> tags produces one HTML
## file that opens standalone in any browser, with the actual PNGs baked
## in rather than referenced as sibling files a user could lose track of.
## `summary_df` is the output of methyl_qc_summary_table(); `subtitle` a
## one-line string (dataset source + stratum) the caller builds since this
## function no longer has one unified result object to read either from.
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

## Every report figure as a standalone PNG, zipped - utils::zip() shells
## out to a `zip` binary on PATH (present by default on macOS/Linux, not
## guaranteed on Windows), so this degrades with a reason rather than
## erroring when it isn't found, same as every other optional-tool path in
## this file.
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
