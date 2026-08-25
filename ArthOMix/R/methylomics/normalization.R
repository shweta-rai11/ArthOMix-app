## R/methylomics/normalization.R
## Methylation normalization methods for the Methylomics "Normalization"
## sub-module (ArthOMix/R/methylomics/mod_methyl_normalization.R). Method
## choice and default parameters are taken from each method's own published
## documentation, not guessed:
##   - Noob:                     minfi::preprocessNoob(offset=15, dyeMethod="single")
##   - Functional normalization: minfi::preprocessFunnorm(nPCs=2)
##   - SWAN:                     minfi::preprocessSWAN()
##   - Stratified quantile:      minfi::preprocessQuantile() - minfi's own
##                                re-implementation of Touleimat & Tost 2012's
##                                subset/stratified quantile method
##   - Dasen:                    wateRmelon::dasen() - wateRmelon's own docs
##                                describe this as "the recommended method"
##   - BMIQ:                     wateRmelon::BMIQ(nfit=50000, nL=3, ...) - the
##                                published Teschendorff et al. 2013 defaults
##   - PBC:                      the Dedeurwaerder et al. 2011 peak-based
##                                correction method as implemented in
##                                ChAMP:::DoPBC() - not exported by ChAMP, but
##                                this is the real, tested implementation
##                                champ.norm(method="PBC") itself calls; safer
##                                to call it directly than to re-derive PBC's
##                                density-peak logic from scratch, and
##                                champ.norm() itself has unwanted side
##                                effects here (dir.create()/setwd() into a
##                                results folder) that calling DoPBC directly
##                                avoids
##   - Quantile (plain):         limma::normalizeQuantiles() - a simple,
##                                non-stratified baseline that works on any
##                                input, unlike every method above
##
## Noob/Funnorm/SWAN/Dasen/Stratified-quantile all need raw IDAT intensities
## (an RGChannelSet/MethylSet) and are unavailable for a matrix-only upload;
## BMIQ/PBC need Type I/II probe-design annotation (450K/EPIC only) but work
## from a beta matrix alone; plain quantile normalization works everywhere.
## Every function below degrades to list(ok = FALSE, reason = ...) rather
## than guessing when its required input isn't available.

## Type I/II design vector (1 = Type I, 2 = Type II) aligned to `probe_ids`,
## used by BMIQ/PBC - both need to know each probe's Infinium chemistry to
## correct for the type-II compression artifact. Probes missing from the
## manifest (shouldn't happen for a genuine 450K/EPIC upload, but possible
## for a subsetted or custom probe list) are dropped, not guessed at.
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

## Per-sample BMIQ (Teschendorff et al. 2013 beta-mixture quantile
## normalization) - wateRmelon::BMIQ() operates one sample at a time.
## Documented to occasionally fail on a poor-quality sample; failures are
## reported per-sample and that sample's ORIGINAL beta values are kept
## rather than dropping it or filling in NA. `nL`/`tol` are exposed
## (wateRmelon::BMIQ()'s own "number of states" and convergence-tolerance
## parameters) since the Normalization tab's BMIQ panel offers them, not
## guessed defaults.
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

  ## BMIQ samples `nfit` probes PER TYPE without replacement to fit its
  ## beta-mixture model - on a small or pre-filtered upload (fewer probes
  ## of one Infinium type than the 50000 default), that sample() call
  ## errors outright rather than gracefully using what's available.
  ## wateRmelon's own docs note smaller nfit values are an expected,
  ## sanctioned tradeoff ("~10000 will make BMIQ run faster at the expense
  ## of a small loss in accuracy"), so capping to what's actually there is
  ## the correct fix, not a workaround.
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
  list(ok = TRUE, beta = out, note = note, failed_samples = failed)
}

## Peak-based correction (Dedeurwaerder et al. 2011), via ChAMP's own
## (unexported) reference implementation - see this file's header comment.
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
## The literature draws a real distinction between background/dye-bias
## correction (what Noob does) and probe-design/distribution normalization
## (what BMIQ/SWAN do) - these two combos run one of each in sequence,
## rather than treating the two kinds of correction as interchangeable.

methyl_norm_noob_bmiq <- function(rg_set, anno_result, offset = 15, dye_method = "single", nfit = 50000, nL = 3, tol = 0.001) {
  step1 <- methyl_norm_noob(rg_set, offset, dye_method)
  if (!isTRUE(step1$ok)) return(step1)
  step2 <- methyl_norm_bmiq(step1$beta, anno_result, nfit = nfit, nL = nL, tol = tol)
  if (!isTRUE(step2$ok)) return(list(ok = FALSE, reason = sprintf("Noob succeeded but the BMIQ step failed: %s", step2$reason)))
  list(ok = TRUE, beta = step2$beta, note = sprintf("Sequential Noob -> BMIQ. Noob: %s BMIQ: %s", step1$note, step2$note))
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

## ---- Method explanations, for the Normalization tab's per-method info ----
## panel (section "Method-specific explanations"). `category` distinguishes
## background/technical correction from probe-design/distribution
## normalization, per this file's own header comment on why those are not
## the same operation.
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
## methyl_get_annotation() (annotation.R) only carries chr/pos/Type/SNP/gene -
## enough for the Quality Control module's own filters. The Normalization
## tab's CpG-island-relation and gene-region filters need two more columns
## (Relation_to_Island, UCSC_RefGene_Group) that live in the SAME underlying
## annotation packages but a DIFFERENT data object (Islands.UCSC/Other) than
## methyl_get_annotation() loads - rather than editing that shared function
## (used by Quality Control) to carry columns only this tab needs, this
## builds its own extension on top of it, cached separately, so
## annotation.R's own behavior/return shape for every other caller is
## completely unchanged.
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
    ## Same "first token of the semicolon-separated list" simplification
    ## methyl_get_annotation() already uses for `gene` - one representative
    ## gene-region label per probe.
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
## Same "return keep + note, never filter in place" convention as qc.R's own
## filters (methyl_filter_snp() etc.) - kept here rather than in qc.R since
## they depend on methyl_get_norm_annotation()'s extra columns, which QC's
## own filters don't use.

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

## Sample-level missingness filter for the Normalization tab's own "remove
## low-quality samples" option - distinct from Quality Control's sample-QC
## tab (which reports what it would remove but never touches the shared
## dataset); this one genuinely drops columns before normalization runs,
## since a sample too sparse to normalize meaningfully shouldn't silently
## ride along.
methyl_filter_samples_missingness <- function(mat, max_na_frac) {
  if (is.null(max_na_frac)) return(list(keep = rep(TRUE, ncol(mat)), note = "No sample-missingness threshold applied."))
  na_frac <- colMeans(is.na(mat))
  keep <- na_frac <= max_na_frac
  list(keep = keep, note = sprintf("%d sample(s) exceed %.0f%% missing values across probes and were removed.", sum(!keep), max_na_frac * 100))
}

## ---- Automatic diagnostics --------------------------------------------
## The compact summary + "Show Details" numbers the Normalization tab shows
## for both the preloaded and uploaded data pathways before any filtering or
## normalization runs - purely descriptive of what's already loaded.

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

## Kolmogorov-Smirnov distance between the pooled Type I and Type II
## beta-value distributions - the actual signal both BMIQ and SWAN exist to
## correct for (Teschendorff et al. 2013's own motivating figure is exactly
## this Type I vs Type II distributional gap). Used both to give the
## "already normalized?" status a real, checkable number instead of a guess,
## and to quantify before/after improvement once a method has been run.
## Sampled (not exhaustive) for speed on a 450K/EPIC-sized matrix; every
## other top-variance-probe computation in this module (qc.R) samples the
## same way, so this is consistent with the rest of the app, not a new
## shortcut.
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
## Deliberately narrower than "has this data been normalized at all", because
## that broader question is NOT answerable from a beta/M-value matrix alone -
## a dataset can be fully corrected for background fluorescence, dye bias,
## and batch effects (a real normalization pipeline) and still show a large
## residual Type I/II gap, because that specific correction (what BMIQ/SWAN
## do) is a separate, narrower step many published normalization pipelines
## never applied - this is exactly the case for this app's own preloaded
## cohort (Liu et al. 2013's GEO-deposited, author-normalized beta values;
## see the Dataset tab), which reliably shows a large Kolmogorov-Smirnov gap
## here despite being genuinely normalized upstream. So:
##   - "raw" is reserved for the one thing this function can state as a
##     STRUCTURAL FACT, not a guess: mod_methyl_dataset.R's IDAT path only
##     ever derives beta via minfi::preprocessRaw(), never normalizes, so
##     rg_set being present means no processing has happened at all.
##   - "bias_detected" / "no_bias_detected" describe ONLY the Type I/II
##     probe-design signal (what methyl_type_bias_stat() measures) - never
##     phrased as "requires normalization" or "already normalized", since
##     this function has no way to know whether background/dye/batch
##     correction happened upstream of the beta matrix it was given.
##   - "unknown" covers every case where even that narrower signal isn't
##     available or is inconclusive.
methyl_norm_status <- function(mat, dataset, anno_result) {
  if (!is.null(dataset$rg_set)) {
    return(list(status = "raw", message =
      "This dataset was derived directly from raw IDAT intensities and has not been normalized. A normalization step is recommended before downstream analysis.",
      bias = NULL))
  }
  if (!identical(dataset$input_scale, "beta")) {
    return(list(status = "unknown", message =
      "This matrix is on the M-value scale, which does not carry a straightforward probe-design-bias signature. Review the diagnostics below before deciding whether to normalize.",
      bias = NULL))
  }
  bias <- methyl_type_bias_stat(mat, anno_result)
  if (!isTRUE(bias$ok)) {
    return(list(status = "unknown", message = sprintf(
      "Unable to determine automatically (%s). Please review the data diagnostics before proceeding.", bias$reason), bias = NULL))
  }
  if (bias$ks_stat < 0.03) {
    return(list(status = "no_bias_detected", message = sprintf(
      "Type I and Type II probe beta-value distributions are closely aligned (Kolmogorov-Smirnov statistic = %.3f) - no evidence of uncorrected Type I/II probe-design bias. This speaks only to that specific signal, not to whether background/dye-bias/batch correction was applied upstream (not determinable from a beta matrix alone).", bias$ks_stat),
      bias = bias))
  }
  if (bias$ks_stat > 0.06) {
    return(list(status = "bias_detected", message = sprintf(
      "Type I and Type II probe beta-value distributions differ substantially (Kolmogorov-Smirnov statistic = %.3f), consistent with uncorrected Type I/II probe-design bias - the kind BMIQ/SWAN specifically address. This does not necessarily mean no normalization was ever applied: a dataset can be fully corrected for background/dye/batch effects and still show this residual difference if a probe-design-aware step was never run on it.", bias$ks_stat),
      bias = bias))
  }
  list(status = "unknown", message = sprintf(
    "The Type I/Type II distribution difference (Kolmogorov-Smirnov statistic = %.3f) is inconclusive. Please review the data diagnostics before proceeding.", bias$ks_stat),
    bias = bias)
}

## Advisory recommendation text (section "Smart method recommendation") -
## always non-binding; the method picker below it never restricts the user
## to whatever this suggests.
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
## Computed once a normalization result exists, from the SAME before/after
## matrices the density plot already shows - not a claim that the algorithm
## merely completed, but an actual comparison of technical-variation and
## structure metrics (section "Normalization validation").
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

## Plain-language interpretation of methyl_norm_validation()'s numbers -
## never declares success merely because the algorithm completed (section
## "Normalization validation"); a biological-signal-preservation warning
## (section "Preserve biological signal") takes priority over an otherwise
## favorable technical-variation readout.
methyl_norm_interpretation <- function(v) {
  if (isTRUE(v$signal_check$flagged)) {
    return(list(status = "warning", text = sprintf(
      "Normalization completed, but a potential biological signal change was detected: PC1's association with the selected group column dropped from R² = %.3f to R² = %.3f. Review the before/after PCA and group-level distributions before selecting this result for downstream analysis.",
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
