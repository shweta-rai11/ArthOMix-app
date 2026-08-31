## R/methylomics/normalization.R
## Methylation normalization methods for the Normalization sub-module
## (mod_methyl_normalization.R). Defaults follow each method's own published docs:
##   - Noob:                 minfi::preprocessNoob(offset=15, dyeMethod="single")
##   - Funnorm:               minfi::preprocessFunnorm(nPCs=2)
##   - SWAN:                  minfi::preprocessSWAN()
##   - Stratified quantile:   minfi::preprocessQuantile() (Touleimat & Tost 2012)
##   - Dasen:                 wateRmelon::dasen() (wateRmelon's recommended default)
##   - BMIQ:                  wateRmelon::BMIQ(nfit=50000, nL=3, ...) (Teschendorff et al. 2013 defaults)
##   - PBC:                   ChAMP:::DoPBC() (Dedeurwaerder et al. 2011); called directly
##                             since champ.norm() has unwanted dir.create()/setwd() side effects
##   - Quantile (plain):      limma::normalizeQuantiles() - simple baseline, works on any input
##
## Noob/Funnorm/SWAN/Dasen/Stratified-quantile need raw IDAT (RGChannelSet/MethylSet);
## BMIQ/PBC need Type I/II annotation (450K/EPIC) but work from a beta matrix; plain
## quantile works everywhere. Each function degrades to list(ok = FALSE, reason = ...).

## Type I/II design vector (1 = Type I, 2 = Type II) aligned to `probe_ids`, used
## by BMIQ/PBC. Probes missing from the manifest are dropped, not guessed at.
methyl_design_vector <- function(probe_ids, anno_result) {
  if (!isTRUE(anno_result$ok)) {
    return(list(ok = FALSE, reason = anno_result$reason))
  }
  a <- anno_result$anno
  hit <- probe_ids %in% rownames(a)
  if (sum(hit) == 0) {
    return(list(ok = FALSE, reason = "None of these probe IDs were found in the manifest annotation."))
  }
  design_v <- rep(NA_real_, length(probe_ids))
  design_v[hit] <- ifelse(a[probe_ids[hit], "Type"] == "I", 1, 2)
  names(design_v) <- probe_ids
  list(ok = TRUE, design_v = design_v, n_dropped = sum(!hit))
}

## ---- IDAT-only methods (need an RGChannelSet) -----------------------------

methyl_norm_noob <- function(rg_set, offset = 15, dye_method = "single") {
  if (is.null(rg_set) || !requireNamespace("minfi", quietly = TRUE)) {
    return(list(ok = FALSE, reason = "Noob requires raw IDAT input."))
  }
  mset <- tryCatch(minfi::preprocessNoob(rg_set, offset = offset, dyeMethod = dye_method), error = function(e) e)
  if (inherits(mset, "error")) return(list(ok = FALSE, reason = paste("preprocessNoob failed:", conditionMessage(mset))))
  list(ok = TRUE, beta = minfi::getBeta(mset), note = sprintf("minfi::preprocessNoob(offset=%d, dyeMethod=\"%s\")", offset, dye_method))
}

methyl_norm_funnorm <- function(rg_set, n_pcs = 2) {
  if (is.null(rg_set) || !requireNamespace("minfi", quietly = TRUE)) {
    return(list(ok = FALSE, reason = "Functional normalization requires raw IDAT input."))
  }
  gmset <- tryCatch(minfi::preprocessFunnorm(rg_set, nPCs = n_pcs), error = function(e) e)
  if (inherits(gmset, "error")) return(list(ok = FALSE, reason = paste("preprocessFunnorm failed:", conditionMessage(gmset))))
  list(ok = TRUE, beta = minfi::getBeta(gmset), note = sprintf("minfi::preprocessFunnorm(nPCs=%d) - recommended when samples span distinct biological groups/tissues.", n_pcs))
}

methyl_norm_swan <- function(rg_set, mset) {
  if (is.null(rg_set) || is.null(mset) || !requireNamespace("minfi", quietly = TRUE)) {
    return(list(ok = FALSE, reason = "SWAN requires raw IDAT input."))
  }
  sset <- tryCatch(minfi::preprocessSWAN(rg_set, mset), error = function(e) e)
  if (inherits(sset, "error")) return(list(ok = FALSE, reason = paste("preprocessSWAN failed:", conditionMessage(sset))))
  list(ok = TRUE, beta = minfi::getBeta(sset), note = "minfi::preprocessSWAN()")
}

methyl_norm_stratified_quantile <- function(rg_set) {
  if (is.null(rg_set) || !requireNamespace("minfi", quietly = TRUE)) {
    return(list(ok = FALSE, reason = "Stratified quantile normalization requires raw IDAT input."))
  }
  gmset <- tryCatch(minfi::preprocessQuantile(rg_set), error = function(e) e)
  if (inherits(gmset, "error")) return(list(ok = FALSE, reason = paste("preprocessQuantile failed:", conditionMessage(gmset))))
  list(ok = TRUE, beta = minfi::getBeta(gmset), note = "minfi::preprocessQuantile() - minfi's re-implementation of Touleimat & Tost (2012) subset/stratified quantile normalization. Recommended for a single tissue/cell type without large global methylation differences between samples.")
}

methyl_norm_dasen <- function(mset) {
  if (is.null(mset) || !requireNamespace("wateRmelon", quietly = TRUE)) {
    return(list(ok = FALSE, reason = "Dasen requires raw IDAT input."))
  }
  ## wateRmelon::dasen() returns a MethylSet (not a beta matrix directly,
  ## despite some of its own documentation examples implying otherwise) -
  ## minfi::getBeta() is required to get an actual beta matrix out of it;
  ## as.matrix() on the MethylSet itself errors (no coercion method).
  mset_out <- tryCatch(wateRmelon::dasen(mset), error = function(e) e)
  if (inherits(mset_out, "error")) return(list(ok = FALSE, reason = paste("wateRmelon::dasen() failed:", conditionMessage(mset_out))))
  beta <- tryCatch(minfi::getBeta(mset_out), error = function(e) e)
  if (inherits(beta, "error")) return(list(ok = FALSE, reason = paste("Could not extract beta values from wateRmelon::dasen()'s result:", conditionMessage(beta))))
  list(ok = TRUE, beta = as.matrix(beta), note = "wateRmelon::dasen() - equalizes Type I/II background before dasen normalization; wateRmelon's own documentation recommends this as the default choice.")
}

## ---- Matrix-based methods (need Type I/II annotation, not raw IDAT) ------

## Per-sample BMIQ (Teschendorff et al. 2013), via wateRmelon::BMIQ() one sample
## at a time. A failed sample keeps its original beta values rather than NA.
methyl_norm_bmiq <- function(mat, anno_result, nfit = 50000, nL = 3, tol = 0.001) {
  if (!requireNamespace("wateRmelon", quietly = TRUE)) {
    return(list(ok = FALSE, reason = "wateRmelon is not installed."))
  }
  dv <- methyl_design_vector(rownames(mat), anno_result)
  if (!isTRUE(dv$ok)) return(list(ok = FALSE, reason = dv$reason))
  keep <- !is.na(dv$design_v)
  m <- mat[keep, , drop = FALSE]
  design_v <- dv$design_v[keep]
  m[m == 0] <- 1e-6; m[m == 1] <- 1 - 1e-6  ## BMIQ is undefined at the beta boundary

  ## BMIQ samples `nfit` probes per type without replacement; caps nfit to what's
  ## actually available to avoid an outright error on a small/pre-filtered upload.
  n_type1 <- sum(design_v == 1, na.rm = TRUE); n_type2 <- sum(design_v == 2, na.rm = TRUE)
  nfit_used <- max(1, min(nfit, n_type1, n_type2))

  out <- m
  failed <- character(0)
  for (s in colnames(m)) {
    res <- tryCatch(wateRmelon::BMIQ(m[, s], design_v, nfit = nfit_used, nL = nL, tol = tol, plots = FALSE, pri = FALSE), error = function(e) e)
    if (inherits(res, "error")) failed <- c(failed, s) else out[, s] <- res$nbeta
  }
  note <- sprintf("wateRmelon::BMIQ(nfit=%d%s) per sample, on %d Type I/II-annotated probe(s) (%d dropped, no manifest match).%s",
                  nfit_used, if (nfit_used < nfit) sprintf(", reduced from %d - fewer probes of one Infinium type than that", nfit) else "",
                  sum(keep), dv$n_dropped,
                  if (length(failed) > 0) sprintf(" %d sample(s) failed and were left unnormalized: %s.", length(failed), paste(failed, collapse = ", ")) else "")
  ## Column-level flag (named by sample) distinguishing actually-BMIQ-corrected
  ## columns from ones left as raw input after a failed fit - `note`'s prose
  ## already says this, but a caller/UI needs a structured flag to act on it
  ## (e.g. badge affected sample columns) rather than parsing free text.
  sample_ok <- stats::setNames(!(colnames(out) %in% failed), colnames(out))
  list(ok = TRUE, beta = out, note = note, failed_samples = failed, sample_ok = sample_ok)
}

## Peak-based correction (Dedeurwaerder et al. 2011) via ChAMP's unexported
## DoPBC() - see file header.
methyl_norm_pbc <- function(mat, anno_result) {
  if (!requireNamespace("ChAMP", quietly = TRUE)) {
    return(list(ok = FALSE, reason = "ChAMP is not installed."))
  }
  dv <- methyl_design_vector(rownames(mat), anno_result)
  if (!isTRUE(dv$ok)) return(list(ok = FALSE, reason = dv$reason))
  keep <- !is.na(dv$design_v)
  m <- mat[keep, , drop = FALSE]
  m[m == 0] <- 1e-6; m[m == 1] <- 1 - 1e-6

  out <- tryCatch(ChAMP:::DoPBC(m, dv$design_v[keep]), error = function(e) e)
  if (inherits(out, "error")) return(list(ok = FALSE, reason = paste("PBC failed:", conditionMessage(out))))
  rownames(out) <- rownames(m); colnames(out) <- colnames(m)
  list(ok = TRUE, beta = out, note = sprintf("Peak-based correction (Dedeurwaerder et al. 2011), on %d Type I/II-annotated probe(s) (%d dropped, no manifest match).", sum(keep), dv$n_dropped))
}

## ---- Universal fallback ----------------------------------------------------

methyl_norm_quantile <- function(mat) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    return(list(ok = FALSE, reason = "limma is not installed."))
  }
  out <- tryCatch(limma::normalizeQuantiles(mat), error = function(e) e)
  if (inherits(out, "error")) return(list(ok = FALSE, reason = paste("Quantile normalization failed:", conditionMessage(out))))
  list(ok = TRUE, beta = out, note = "limma::normalizeQuantiles() - a plain (non-stratified) quantile normalization; works on any input but doesn't correct for the Type I/II probe-design bias the other methods do.")
}

## ---- Sequential (two-step) workflows --------------------------------------
## Background/dye-bias correction (Noob) and probe-design normalization (BMIQ/SWAN)
## are distinct steps - these combos run one of each in sequence.

methyl_norm_noob_bmiq <- function(rg_set, anno_result, offset = 15, dye_method = "single", nfit = 50000, nL = 3, tol = 0.001) {
  step1 <- methyl_norm_noob(rg_set, offset, dye_method)
  if (!isTRUE(step1$ok)) return(step1)
  step2 <- methyl_norm_bmiq(step1$beta, anno_result, nfit = nfit, nL = nL, tol = tol)
  if (!isTRUE(step2$ok)) return(list(ok = FALSE, reason = sprintf("Noob succeeded but the BMIQ step failed: %s", step2$reason)))
  list(ok = TRUE, beta = step2$beta, note = sprintf("Sequential Noob -> BMIQ. Noob: %s BMIQ: %s", step1$note, step2$note),
       failed_samples = step2$failed_samples, sample_ok = step2$sample_ok)
}

methyl_norm_noob_swan <- function(rg_set, offset = 15, dye_method = "single") {
  if (is.null(rg_set) || !requireNamespace("minfi", quietly = TRUE)) {
    return(list(ok = FALSE, reason = "Noob + SWAN requires raw IDAT input."))
  }
  mset <- tryCatch(minfi::preprocessNoob(rg_set, offset = offset, dyeMethod = dye_method), error = function(e) e)
  if (inherits(mset, "error")) return(list(ok = FALSE, reason = paste("preprocessNoob failed:", conditionMessage(mset))))
  sset <- tryCatch(minfi::preprocessSWAN(rg_set, mset), error = function(e) e)
  if (inherits(sset, "error")) return(list(ok = FALSE, reason = paste("preprocessSWAN failed:", conditionMessage(sset))))
  list(ok = TRUE, beta = minfi::getBeta(sset),
       note = sprintf("Sequential Noob (offset=%d, dyeMethod=\"%s\") -> SWAN.", offset, dye_method))
}

## ---- Method explanations for the Normalization tab's per-method info panel ----
## `category` distinguishes background/technical correction from probe-design normalization.
METHYL_NORM_METHOD_INFO <- list(
  noob = list(category = "Background / technical correction", text =
    "Noob (normal-exponential out-of-band) performs background correction and dye-bias correction using out-of-band Infinium I probe intensities. Appropriate whenever raw methylation-array intensity data are available; on its own it does not address Type I/II probe-design distribution differences."),
  funnorm = list(category = "Background / technical correction", text =
    "Functional normalization uses the array's built-in control probes (summarized via PCA) to estimate and remove unwanted technical variation, on top of Noob-style background correction. Useful when samples span distinct biological groups or tissues, since - unlike quantile-based methods - it does not assume similar global methylation across samples."),
  swan = list(category = "Probe-design / distribution normalization", text =
    "SWAN (subset-quantile within-array normalization) matches the beta-value distributions of Type I and Type II probes with similar CpG density within each array. Addresses Infinium probe-design effects; it does not perform background/dye correction on its own."),
  stratified_quantile = list(category = "Probe-design / distribution normalization", text =
    "Stratified (subset) quantile normalization - minfi's re-implementation of Touleimat & Tost (2012) - quantile-normalizes probe subsets stratified by region/probe type. Recommended for a single tissue/cell type without large expected global methylation differences between samples; it can distort real, large biological differences if they are present."),
  dasen = list(category = "Probe-design / distribution normalization", text =
    "Dasen separately quantile-normalizes methylated/unmethylated intensities within each probe type, then recombines them - wateRmelon's own documentation describes it as their recommended default for Illumina methylation-array data."),
  bmiq = list(category = "Probe-design / distribution normalization", text =
    "BMIQ (beta-mixture quantile normalization, Teschendorff et al. 2013) fits a three-state beta-mixture model to Type II probes and transforms them onto the Type I distribution. Targets Type I/II probe-design bias specifically - it is not a general-purpose normalization method for every methylation dataset, and requires Type I/II probe-design annotation."),
  pbc = list(category = "Probe-design / distribution normalization", text =
    "Peak-based correction (Dedeurwaerder et al. 2011) aligns the density peaks of the Type I and Type II probe distributions. An alternative probe-design correction to BMIQ/SWAN; also requires Type I/II probe-design annotation."),
  quantile = list(category = "Universal baseline", text =
    "Plain (non-stratified) quantile normalization forces every sample's beta-value distribution to match the same reference distribution. Works on any input, but makes a strong assumption that samples don't differ substantially in overall methylation - use carefully when biologically meaningful global methylation differences are expected, since it can remove real signal along with technical noise."),
  noob_bmiq = list(category = "Sequential workflow", text =
    "Runs Noob (background/dye-bias correction) first, then BMIQ (Type I/II probe-design correction) on the result - a two-step workflow that separately addresses two different sources of technical variation rather than treating them as one operation."),
  noob_swan = list(category = "Sequential workflow", text =
    "Runs Noob (background/dye-bias correction) first, then SWAN (Type I/II probe-design correction) on the result - the same two-step logic as Noob + BMIQ, using SWAN instead of BMIQ for the probe-design correction step.")
)

## ---- Extended (Normalization-tab-only) annotation -------------------------
## methyl_get_annotation() only carries chr/pos/Type/SNP/gene. This tab also needs
## island-relation and gene-region columns (from a different data object), so this
## builds its own cached extension rather than changing the shared function.
.methyl_norm_anno_cache <- new.env(parent = emptyenv())

methyl_get_norm_annotation <- function(array_type) {
  base <- methyl_get_annotation(array_type)
  if (!isTRUE(base$ok)) return(base)
  pkg <- METHYL_ANNOTATION_PACKAGES[[array_type]]
  cached <- .methyl_norm_anno_cache[[pkg]]
  if (!is.null(cached)) return(list(ok = TRUE, anno = cached, reason = NULL))

  ext <- tryCatch({
    e <- new.env(parent = emptyenv())
    utils::data(list = c("Islands.UCSC", "Other"), package = pkg, envir = e)
    isl <- as.data.frame(e$Islands.UCSC)
    oth <- as.data.frame(e$Other)
    a <- base$anno
    ids <- rownames(a)
    a$island_relation <- isl[ids, "Relation_to_Island"]
    ## Same first-token simplification methyl_get_annotation() uses for `gene`.
    a$gene_region <- vapply(strsplit(oth[ids, "UCSC_RefGene_Group"], ";"), function(g) {
      if (length(g) == 0 || !nzchar(g[1])) NA_character_ else g[1]
    }, character(1))
    a
  }, error = function(e) e)

  ## Degrades to the base annotation (no island_relation/gene_region
  ## columns) rather than failing outright - callers check for the
  ## column's presence before offering that specific filter.
  if (inherits(ext, "error")) return(list(ok = TRUE, anno = base$anno, reason = NULL))
  .methyl_norm_anno_cache[[pkg]] <- ext
  list(ok = TRUE, anno = ext, reason = NULL)
}

## ---- Normalization-tab-only probe filters ----------------------------------
## Same "return keep + note" convention as qc.R's filters; kept here since these
## depend on methyl_get_norm_annotation()'s extra columns.

methyl_filter_island_relation <- function(mat, anno_result, keep_categories) {
  a <- anno_result$anno
  if (is.null(a) || !("island_relation" %in% colnames(a)) || length(keep_categories) == 0) {
    return(list(keep = rep(TRUE, nrow(mat)), note = "CpG island-relation annotation is unavailable for this array type, or no category was selected."))
  }
  ids <- rownames(mat)
  hit <- ids %in% rownames(a)
  rel <- rep(NA_character_, length(ids)); rel[hit] <- as.character(a[ids[hit], "island_relation"])
  removed <- hit & !(rel %in% keep_categories)
  list(keep = !removed, note = sprintf("%d probe(s) outside the selected island-relation categories removed (%d unannotated probe(s) kept).", sum(removed), sum(!hit)))
}

methyl_filter_gene_region <- function(mat, anno_result, keep_regions) {
  a <- anno_result$anno
  if (is.null(a) || !("gene_region" %in% colnames(a)) || length(keep_regions) == 0) {
    return(list(keep = rep(TRUE, nrow(mat)), note = "Gene-region annotation is unavailable for this array type, or no region was selected."))
  }
  ids <- rownames(mat)
  hit <- ids %in% rownames(a)
  reg <- rep(NA_character_, length(ids)); reg[hit] <- as.character(a[ids[hit], "gene_region"])
  removed <- hit & !(reg %in% keep_regions)
  list(keep = !removed, note = sprintf("%d probe(s) outside the selected gene-region categories removed (%d intergenic/unannotated probe(s) kept).", sum(removed), sum(!hit)))
}

methyl_filter_chromosome <- function(mat, anno_result, exclude_chr) {
  if (length(exclude_chr) == 0) return(list(keep = rep(TRUE, nrow(mat)), note = "No chromosome exclusion selected."))
  if (!isTRUE(anno_result$ok)) return(list(keep = rep(TRUE, nrow(mat)), note = anno_result$reason))
  a <- anno_result$anno
  ids <- rownames(mat)
  hit <- ids %in% rownames(a)
  chr <- rep(NA_character_, length(ids)); chr[hit] <- a[ids[hit], "chr"]
  removed <- hit & chr %in% exclude_chr
  list(keep = !removed, note = sprintf("%d probe(s) on excluded chromosome(s) (%s) removed.", sum(removed), paste(exclude_chr, collapse = ", ")))
}

## Sample-missingness filter for the "remove low-quality samples" option. Unlike
## QC's sample-QC tab (report-only), this one actually drops columns before normalizing.
methyl_filter_samples_missingness <- function(mat, max_na_frac) {
  if (is.null(max_na_frac)) return(list(keep = rep(TRUE, ncol(mat)), note = "No sample-missingness threshold applied."))
  na_frac <- colMeans(is.na(mat))
  keep <- na_frac <= max_na_frac
  list(keep = keep, note = sprintf("%d sample(s) exceed %.0f%% missing values across probes and were removed.", sum(!keep), max_na_frac * 100))
}

## ---- Automatic diagnostics --------------------------------------------
## Summary numbers the Normalization tab shows before any filtering/normalization runs.

methyl_norm_diagnostics <- function(mat, dataset, anno_result = NULL) {
  vals <- mat[is.finite(mat)]
  representation_short <- if (identical(dataset$input_scale, "beta")) "Beta values"
                           else if (identical(dataset$input_scale, "m")) "M-values"
                           else "Unknown"
  representation <- if (identical(dataset$input_scale, "beta")) "Beta values (0-1 methylation proportion)"
                     else if (identical(dataset$input_scale, "m")) "M-values (logit-transformed)"
                     else "Unknown representation"
  list(
    n_probes = nrow(mat), n_samples = ncol(mat),
    n_missing = sum(is.na(mat)), pct_missing = 100 * mean(is.na(mat)),
    n_infinite = sum(is.infinite(mat)),
    value_min = if (length(vals) > 0) min(vals) else NA_real_,
    value_max = if (length(vals) > 0) max(vals) else NA_real_,
    median = if (length(vals) > 0) stats::median(vals) else NA_real_,
    mean = if (length(vals) > 0) mean(vals) else NA_real_,
    sd = if (length(vals) > 0) stats::sd(vals) else NA_real_,
    representation = representation, representation_short = representation_short,
    platform = dataset$array_type %||% "Unknown",
    has_raw_intensity = !is.null(dataset$rg_set),
    has_beta = identical(dataset$input_scale, "beta"),
    has_mvalue = identical(dataset$input_scale, "m"),
    has_probe_type_info = isTRUE(anno_result$ok) && !is.null(anno_result$anno) && "Type" %in% colnames(anno_result$anno)
  )
}

## Kolmogorov-Smirnov distance between pooled Type I and Type II beta distributions -
## the signal BMIQ/SWAN correct for (Teschendorff et al. 2013). Used both for the
## "already normalized?" status and to quantify before/after improvement. Sampled
## (not exhaustive) for speed, consistent with qc.R's sampling elsewhere in this module.
methyl_type_bias_stat <- function(mat, anno_result) {
  if (!isTRUE(anno_result$ok) || is.null(anno_result$anno) || !("Type" %in% colnames(anno_result$anno))) {
    return(list(ok = FALSE, reason = "No Type I/II probe-design annotation available for this array type."))
  }
  a <- anno_result$anno
  ids <- intersect(rownames(mat), rownames(a))
  if (length(ids) < 100) return(list(ok = FALSE, reason = "Too few annotated probes to estimate Type I/II distribution bias."))
  type1_ids <- ids[a[ids, "Type"] == "I"]
  type2_ids <- ids[a[ids, "Type"] == "II"]
  if (length(type1_ids) < 50 || length(type2_ids) < 50) {
    return(list(ok = FALSE, reason = "Too few Type I or Type II annotated probes to estimate distribution bias."))
  }
  n_each <- min(5000, length(type1_ids), length(type2_ids))
  t1 <- type1_ids[order(stats::runif(length(type1_ids)))[seq_len(n_each)]]
  t2 <- type2_ids[order(stats::runif(length(type2_ids)))[seq_len(n_each)]]
  v1 <- mat[t1, , drop = FALSE]; v1 <- v1[is.finite(v1)]
  v2 <- mat[t2, , drop = FALSE]; v2 <- v2[is.finite(v2)]
  if (length(v1) < 100 || length(v2) < 100) return(list(ok = FALSE, reason = "Too few finite values to compare Type I/II distributions."))
  ks <- suppressWarnings(stats::ks.test(v1, v2))
  list(ok = TRUE, ks_stat = unname(ks$statistic), n_type1 = length(t1), n_type2 = length(t2),
       median_type1 = stats::median(v1), median_type2 = stats::median(v2))
}

## ---- Normalization-status detection ---------------------------------------
## Deliberately narrower than "has this been normalized" - that's not answerable
## from a beta matrix alone (e.g. this app's preloaded Liu et al. 2013 cohort is
## genuinely normalized upstream but still shows a large Type I/II gap here).
##   - "raw": a structural fact - rg_set present means beta came only from
##     minfi::preprocessRaw(), no processing at all.
##   - "bias_detected"/"no_bias_detected": describe only the Type I/II signal
##     from methyl_type_bias_stat(), never "requires/already normalized".
##   - "unknown": signal unavailable or inconclusive.
methyl_norm_status <- function(mat, dataset, anno_result) {
  if (!is.null(dataset$rg_set)) {
    return(list(status = "raw", message =
      "Raw, unnormalized data - normalization is recommended before downstream analysis.",
      bias = NULL))
  }
  if (!identical(dataset$input_scale, "beta")) {
    return(list(status = "unknown", message =
      "M-value scale - no probe-design-bias signature available; review diagnostics below.",
      bias = NULL))
  }
  bias <- methyl_type_bias_stat(mat, anno_result)
  if (!isTRUE(bias$ok)) {
    return(list(status = "unknown", message = sprintf(
      "Unable to determine automatically (%s).", bias$reason), bias = NULL))
  }
  if (bias$ks_stat < 0.03) {
    return(list(status = "no_bias_detected", message = sprintf(
      "No Type I/II probe-design bias detected (KS = %.3f).", bias$ks_stat),
      bias = bias))
  }
  if (bias$ks_stat > 0.06) {
    return(list(status = "bias_detected", message = sprintf(
      "Type I/II probe-design bias detected (KS = %.3f) - a BMIQ/SWAN-style correction would address it.", bias$ks_stat),
      bias = bias))
  }
  list(status = "unknown", message = sprintf(
    "Type I/II distribution difference is inconclusive (KS = %.3f).", bias$ks_stat),
    bias = bias)
}

## Advisory recommendation text - always non-binding, never restricts the method picker.
methyl_norm_recommendation <- function(dataset, status, available_methods) {
  if (!is.null(dataset$rg_set)) {
    return("Your dataset contains raw Illumina methylation-array intensity data. A Noob-based workflow (optionally paired with BMIQ or SWAN for probe-design correction) is available and is a reasonable default for correcting background and dye-bias effects before downstream analysis.")
  }
  if (identical(status$status, "no_bias_detected")) {
    return("Your dataset shows no evidence of uncorrected Type I/II probe-design bias. Re-normalizing is usually unnecessary and can occasionally distort an already-corrected distribution - consider keeping the current normalization unless you have a specific reason to reprocess it.")
  }
  if (identical(status$status, "bias_detected") && "bmiq" %in% available_methods) {
    return("Your dataset shows evidence of uncorrected Type I/II probe-design bias. Raw-intensity preprocessing methods such as Noob cannot be applied to a beta/M-value matrix directly; BMIQ or PBC (both beta-value-based, probe-design-aware methods) may be considered.")
  }
  if ("bmiq" %in% available_methods) {
    return("Your dataset contains beta values with Type I/II probe-design annotation available, but no raw intensity channels. Raw-intensity preprocessing methods such as Noob cannot be applied directly; BMIQ or PBC (both beta-value-based, probe-design-aware methods) may be considered.")
  }
  "Your dataset contains beta or M-values without raw intensity channels or Type I/II probe-design annotation. Raw-intensity and probe-design-aware methods are unavailable here; plain quantile normalization is the only method compatible with this input."
}

## ---- Post-run validation ----------------------------------------------
## Compares before/after technical-variation and structure metrics, not just
## whether the algorithm completed (section "Normalization validation").
methyl_norm_validation <- function(before, after, anno_result, group_labels = NULL) {
  out <- list(
    var_before = mean(methyl_row_vars(before), na.rm = TRUE),
    var_after = mean(methyl_row_vars(after), na.rm = TRUE),
    missing_before = 100 * mean(is.na(before)),
    missing_after = 100 * mean(is.na(after))
  )
  cb <- methyl_sample_correlation(before); ca <- methyl_sample_correlation(after)
  out$mean_cor_before <- if (isTRUE(cb$ok)) mean(cb$cor[upper.tri(cb$cor)]) else NA_real_
  out$mean_cor_after <- if (isTRUE(ca$ok)) mean(ca$cor[upper.tri(ca$cor)]) else NA_real_
  bb <- methyl_type_bias_stat(before, anno_result); ba <- methyl_type_bias_stat(after, anno_result)
  out$ks_before <- if (isTRUE(bb$ok)) bb$ks_stat else NA_real_
  out$ks_after <- if (isTRUE(ba$ok)) ba$ks_stat else NA_real_

  out$signal_check <- NULL
  if (!is.null(group_labels) && length(unique(stats::na.omit(group_labels))) >= 2) {
    r2_for <- function(mat_side) {
      pca <- methyl_pca_scores(mat_side)
      if (!isTRUE(pca$ok)) return(NA_real_)
      g <- group_labels[rownames(pca$scores)]
      df <- data.frame(pc1 = pca$scores[, 1], g = g)
      df <- df[!is.na(df$g), , drop = FALSE]
      if (length(unique(df$g)) < 2 || nrow(df) < 4) return(NA_real_)
      fit <- tryCatch(stats::lm(pc1 ~ g, data = df), error = function(e) NULL)
      if (is.null(fit)) NA_real_ else summary(fit)$r.squared
    }
    r2_before <- r2_for(before); r2_after <- r2_for(after)
    if (!is.na(r2_before) && !is.na(r2_after)) {
      drop_frac <- if (r2_before > 1e-6) (r2_before - r2_after) / r2_before else NA_real_
      flagged <- !is.na(drop_frac) && drop_frac > 0.5 && r2_before > 0.05
      out$signal_check <- list(r2_before = r2_before, r2_after = r2_after, flagged = flagged)
    }
  }
  out
}

## Plain-language interpretation of methyl_norm_validation()'s numbers - a
## biological-signal-preservation warning takes priority over a favorable
## technical-variation readout.
methyl_norm_interpretation <- function(v) {
  if (isTRUE(v$signal_check$flagged)) {
    return(list(status = "warning", text = sprintf(
      "Normalized, but PC1's association with the group column dropped from R² = %.3f to %.3f - review before selecting this result.",
      v$signal_check$r2_before, v$signal_check$r2_after)))
  }
  improved <- character(0)
  if (!is.na(v$ks_before) && !is.na(v$ks_after) && v$ks_after < v$ks_before) improved <- c(improved, "reduced Type I/II probe-design bias")
  if (!is.na(v$mean_cor_before) && !is.na(v$mean_cor_after) && v$mean_cor_after > v$mean_cor_before) improved <- c(improved, "increased sample-to-sample correlation")
  if (length(improved) > 0) {
    return(list(status = "pass", text = sprintf(
      "Normalization completed successfully. QC diagnostics indicate %s while preserving overall sample structure.", paste(improved, collapse = " and "))))
  }
  list(status = "neutral", text =
    "Normalization completed. QC diagnostics show limited change on the metrics tracked here - consider comparing an alternative compatible normalization method.")
}

## ---- Reproducibility record -------------------------------------------
methyl_norm_processing_record <- function(r) {
  lines <- c(
    sprintf("Normalization method: %s", r$method_label),
    sprintf("Input dataset: %s", r$dataset_source %||% "Unlabeled dataset"),
    sprintf("Date/time: %s", format(r$run_at, "%Y-%m-%d %H:%M:%S")),
    sprintf("Samples: %d -> %d", r$n_samples_before, r$n_samples_after),
    sprintf("Probes: %s -> %s", format(r$n_probes_before, big.mark = ","), format(r$n_probes_after, big.mark = ",")),
    "Filters applied:",
    if (length(r$filter_notes) > 0) paste0("  - ", r$filter_notes) else "  - none",
    sprintf("Method note: %s", r$note)
  )
  paste(lines, collapse = "\n")
}
