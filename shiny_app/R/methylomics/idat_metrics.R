## R/methylomics/idat_metrics.R
## QC metrics computable only from raw IDAT intensities (an RGChannelSet),
## not from an uploaded beta/M-value matrix: detection p-values, bead
## counts, bisulfite conversion efficiency, and median methylated/
## unmethylated intensity. Every function here degrades to
## list(ok = FALSE, reason = ...) when its input or required package isn't
## available, rather than throwing into a raw Shiny error screen.

## Raw beta matrix (minfi::preprocessRaw(), no normalization - normalization
## methods live in a later Methylomics sub-module) plus the per-probe
## detection p-value and bead-count matrices used by the probe filters in
## R/methylomics/qc.R.
methyl_idat_derive <- function(rg_set) {
  if (is.null(rg_set) || !requireNamespace("minfi", quietly = TRUE)) {
    return(list(ok = FALSE, reason = "minfi is required to process raw IDAT input."))
  }
  mset <- tryCatch(minfi::preprocessRaw(rg_set), error = function(e) e)
  if (inherits(mset, "error")) {
    return(list(ok = FALSE, reason = paste("Could not derive methylated/unmethylated intensities:", conditionMessage(mset))))
  }
  beta <- tryCatch(minfi::getBeta(mset), error = function(e) NULL)
  if (is.null(beta)) return(list(ok = FALSE, reason = "Could not compute beta values from this IDAT upload."))

  detp <- tryCatch(minfi::detectionP(rg_set), error = function(e) NULL)
  beadcount <- tryCatch(minfi::beadcount(rg_set), error = function(e) NULL)

  list(ok = TRUE, mset = mset, beta = beta, detp = detp, beadcount = beadcount)
}

## Per-sample bisulfite conversion efficiency (%) from the array's built-in
## conversion control probes - wateRmelon::bscon(), the standard method for
## this metric on Illumina methylation arrays.
methyl_bisulfite_conversion <- function(rg_set) {
  if (is.null(rg_set) || !requireNamespace("wateRmelon", quietly = TRUE)) {
    return(list(ok = FALSE, reason = "Bisulfite conversion efficiency requires raw IDAT input and the wateRmelon package."))
  }
  pct <- tryCatch(wateRmelon::bscon(rg_set), error = function(e) NULL)
  if (is.null(pct)) return(list(ok = FALSE, reason = "Could not compute bisulfite conversion efficiency (wateRmelon::bscon) for this upload."))
  list(ok = TRUE, pct = pct)
}

## Per-sample median methylated / unmethylated intensity - minfi::getQC(),
## the standard low-intensity-sample screen (Illumina/minfi's own tutorial
## flags samples with median log2 intensity well below the cohort).
methyl_median_intensity <- function(mset) {
  if (is.null(mset) || !requireNamespace("minfi", quietly = TRUE)) {
    return(list(ok = FALSE, reason = "Median intensity requires raw IDAT input."))
  }
  qc <- tryCatch(minfi::getQC(mset), error = function(e) NULL)
  if (is.null(qc)) return(list(ok = FALSE, reason = "Could not compute median intensities for this IDAT upload."))
  qc <- as.data.frame(qc)
  list(ok = TRUE, detail = data.frame(sample = rownames(qc), med_meth_log2 = qc$mMed, med_unmeth_log2 = qc$uMed, row.names = NULL))
}
