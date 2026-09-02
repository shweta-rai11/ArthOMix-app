## R/transcriptomics/expression_type.R
## Shared expression-matrix scale/type heuristics and the declare-then-verify
## upload validator for the Transcriptomics module. Promoted out of
## mod_dge.R (which used to keep looks_like_raw_counts()/
## looks_like_normalized_totals() as local closures, and only checked them at
## analysis run-time) so mod_dataset.R's own upload path, mod_dge.R's
## decoupled "Upload your own data" path, and mod_deconvolution.R's run gate
## can all share one implementation instead of drifting duplicates -
## mod_deconvolution.R previously kept an independent near-copy of
## looks_like_raw_counts() under the name is_linear_counts(), now removed in
## favor of calling looks_like_raw_counts() directly.
##
## Mirrors R/methylomics/parse_upload.R's methyl_validate_matrix_upload():
## never throws, always returns list(ok=, mat=, error=, note=), hard-blocks
## on an unambiguous declared-vs-actual mismatch, and only warns (does not
## block) on an ambiguous one.

## Flags raw-count-like data: non-negative with a wide linear range (99th
## percentile > 100). Not an integer check, since count matrices like RSEM
## output are often fractional.
looks_like_raw_counts <- function(m) {
  vals <- as.numeric(m)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0 || any(vals < 0)) return(FALSE)
  q99 <- suppressWarnings(stats::quantile(vals[vals > 0], 0.99, na.rm = TRUE))
  isTRUE(!is.na(q99) && q99 > 100)
}

## Flags TPM/FPKM/CPM-normalized data: per-sample column sums cluster tightly
## around a fixed target (1e2/1e4/1e6), unlike raw counts' varying library sizes.
looks_like_normalized_totals <- function(m) {
  csums <- colSums(m, na.rm = TRUE)
  csums <- csums[is.finite(csums) & csums > 0]
  if (length(csums) < 2) return(FALSE)
  cv <- stats::sd(csums) / mean(csums)
  pinned <- any(vapply(c(1e2, 1e4, 1e6), function(target) abs(mean(csums) - target) / target < 0.05, logical(1)))
  isTRUE(cv < 0.05 && pinned)
}

## Light heuristic (not a full schema validator): does this "expression
## matrix" actually look like a differential-results table (genes keyed to
## logFC/p-value/padj/baseMean/AveExpr columns) rather than a sample-level
## expression matrix? Checked against the parsed matrix's column names -
## for a real results table uploaded via the same "feature ID in column 1,
## rest of the columns numeric" convention this app's upload parsers use,
## those trailing columns are exactly these stat names instead of sample IDs.
tx_looks_like_results_table <- function(mat) {
  cn <- colnames(mat)
  if (is.null(cn) || length(cn) == 0) return(FALSE)
  results_terms <- "log2?fc|log2foldchange|p\\.?value|padj|p\\.adj|adj\\.?p|basemean|ave ?expr|avgexpr"
  any(grepl(results_terms, cn, ignore.case = TRUE))
}

## Fail-soft declare-then-verify check for an uploaded expression matrix:
## does it actually look like what `declared_type` ("raw" / "normalized" /
## "logtransformed") says it is? Hard-blocks on an unambiguous mismatch
## (negative or TPM/FPKM/CPM-pinned values declared "raw"; raw-count-shaped
## values declared "normalized"/"logtransformed"; a results-table-shaped
## input regardless of declared type), and only warns (does not block) on an
## ambiguous case (e.g. declared "raw" but the value range doesn't clearly
## look count-like - could just be a small/filtered gene panel). Never
## throws; always returns list(ok, mat, error, note).
tx_validate_expr_upload <- function(mat, declared_type) {
  if (isTRUE(tx_looks_like_results_table(mat))) {
    return(list(ok = FALSE, error = "This looks like a differential-expression results table (column names like logFC/p-value/padj/baseMean), not a sample-level expression matrix - each row should be a gene and each column a sample. Upload the underlying expression matrix instead."))
  }

  notes <- character(0)
  if (!is.null(ncol(mat)) && ncol(mat) < 3) {
    notes <- c(notes, "This matrix has fewer than 3 sample columns - double check this is really a sample-level expression matrix, not a summary/results table.")
  }

  vals <- mat[is.finite(mat)]
  if (length(vals) == 0) {
    return(list(ok = FALSE, error = "No finite numeric values found in this matrix - check it isn't entirely NA/blank."))
  }

  has_negative <- any(vals < 0)
  is_raw <- looks_like_raw_counts(mat)
  is_norm_totals <- looks_like_normalized_totals(mat)

  if (identical(declared_type, "raw")) {
    if (has_negative) {
      return(list(ok = FALSE, error = "\"Raw counts\" is selected as the data type, but this matrix has negative values - raw sequencing counts can't be negative. This looks like normalized or log-transformed data instead; change \"Data type\" above, or upload the actual raw count matrix."))
    }
    if (is_norm_totals) {
      return(list(ok = FALSE, error = "\"Raw counts\" is selected as the data type, but this matrix's per-sample totals are tightly pinned near a fixed value (e.g. ~1e6) - the signature of TPM/FPKM/CPM-normalized expression, not raw sequencing counts. Change \"Data type\" above to \"Normalized\", or upload the actual raw count matrix."))
    }
    if (!is_raw) {
      notes <- c(notes, "Note: this data doesn't show the usual wide dynamic range of raw sequencing counts (99th percentile of values <= 100) - double-check \"Data type\" above if downstream results look off.")
    }
  } else if (identical(declared_type, "normalized")) {
    if (is_raw && !is_norm_totals) {
      return(list(ok = FALSE, error = "\"Normalized (TPM/FPKM/CPM)\" is selected as the data type, but this data looks like raw, un-normalized sequencing counts (wide, unpinned per-sample totals), not normalized expression. Change \"Data type\" above to \"Raw counts\", or upload the actual normalized matrix."))
    }
    if (has_negative) {
      notes <- c(notes, "Note: this matrix has negative values, which is unusual for TPM/FPKM/CPM-normalized data (though possible after further transformation) - double-check \"Data type\" above if this is actually log-transformed data.")
    }
  } else if (identical(declared_type, "logtransformed")) {
    if (is_raw && !is_norm_totals) {
      return(list(ok = FALSE, error = "\"Already log-transformed\" is selected as the data type, but this data looks like raw, un-normalized sequencing counts, not log-transformed values. Change \"Data type\" above to \"Raw counts\", or upload the actual log-transformed matrix."))
    }
  }
  ## No declared_type at all (NA/unset, e.g. an older upload path or a
  ## programmatic caller) - accept as-is, heuristic-only inference downstream.

  list(ok = TRUE, mat = mat, note = if (length(notes) > 0) paste(notes, collapse = " ") else NULL)
}
