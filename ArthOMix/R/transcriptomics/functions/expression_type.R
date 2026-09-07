## R/transcriptomics/functions/expression_type.R
## Shared expression-matrix scale/type heuristics and upload validator.

looks_like_raw_counts <- function(m) {
  vals <- as.numeric(m)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0 || any(vals < 0)) return(FALSE)
  q99 <- arthomix_quiet(stats::quantile(vals[vals > 0], 0.99, na.rm = TRUE))
  isTRUE(!is.na(q99) && q99 > 100)
}

looks_like_normalized_totals <- function(m) {
  csums <- colSums(m, na.rm = TRUE)
  csums <- csums[is.finite(csums) & csums > 0]
  if (length(csums) < 2) return(FALSE)
  cv <- stats::sd(csums) / mean(csums)
  pinned <- any(vapply(c(1e2, 1e4, 1e6), function(target) abs(mean(csums) - target) / target < 0.05, logical(1)))
  isTRUE(cv < 0.05 && pinned)
}

deconv_is_linear_scale <- function(declared_type, expr) {
  if (!is.null(declared_type) && !is.na(declared_type) && nzchar(declared_type)) {
    declared_type %in% c("raw", "normalized")
  } else {
    looks_like_raw_counts(expr) || looks_like_normalized_totals(expr)
  }
}

tx_looks_like_results_table <- function(mat) {
  cn <- colnames(mat)
  if (is.null(cn) || length(cn) == 0) return(FALSE)
  results_terms <- "log2?fc|log2foldchange|p\\.?value|padj|p\\.adj|adj\\.?p|basemean|ave ?expr|avgexpr"
  any(grepl(results_terms, cn, ignore.case = TRUE))
}

tx_validate_expr_upload <- function(mat, declared_type) {
  if (isTRUE(tx_looks_like_results_table(mat))) {
    return(list(ok = FALSE, error = "This looks like a differential-expression results table (columns like logFC/p-value/padj/baseMean), not a sample-level expression matrix. Each row should be a gene and each column a sample. Upload the expression matrix instead."))
  }

  notes <- character(0)
  if (!is.null(ncol(mat)) && ncol(mat) < 3) {
    notes <- c(notes, "This matrix has fewer than 3 sample columns. Double check this is a sample-level expression matrix, not a summary/results table.")
  }

  vals <- mat[is.finite(mat)]
  if (length(vals) == 0) {
    return(list(ok = FALSE, error = "No finite numeric values found in this matrix. Check it isn't entirely NA/blank."))
  }

  has_negative <- any(vals < 0)
  is_raw <- looks_like_raw_counts(mat)
  is_norm_totals <- looks_like_normalized_totals(mat)

  if (identical(declared_type, "raw")) {
    if (has_negative) {
      return(list(ok = FALSE, error = "\"Raw counts\" is selected, but this matrix has negative values. Raw counts can't be negative - this looks like normalized or log-transformed data. Change \"Data type\" above, or upload the actual raw count matrix."))
    }
    if (is_norm_totals) {
      return(list(ok = FALSE, error = "\"Raw counts\" is selected, but per-sample totals are tightly pinned near a fixed value (e.g. ~1e6). That's the signature of TPM/FPKM/CPM-normalized data, not raw counts. Change \"Data type\" to \"Normalized\", or upload the actual raw count matrix."))
    }
    if (!is_raw) {
      notes <- c(notes, "Note: this data lacks the usual wide dynamic range of raw sequencing counts (99th percentile <= 100). Double-check \"Data type\" above if downstream results look off.")
    }
  } else if (identical(declared_type, "normalized")) {
    if (is_raw && !is_norm_totals) {
      return(list(ok = FALSE, error = "\"Normalized (TPM/FPKM/CPM)\" is selected, but this looks like raw, un-normalized sequencing counts (wide, unpinned per-sample totals). Change \"Data type\" to \"Raw counts\", or upload the actual normalized matrix."))
    }
    if (has_negative) {
      notes <- c(notes, "Note: this matrix has negative values, unusual for TPM/FPKM/CPM data (though possible after further transformation). Double-check \"Data type\" above if this is actually log-transformed data.")
    }
  } else if (identical(declared_type, "logtransformed")) {
    if (is_raw && !is_norm_totals) {
      return(list(ok = FALSE, error = "\"Already log-transformed\" is selected, but this looks like raw, un-normalized sequencing counts, not log-transformed values. Change \"Data type\" to \"Raw counts\", or upload the actual log-transformed matrix."))
    }
  }

  list(ok = TRUE, mat = mat, note = if (length(notes) > 0) paste(notes, collapse = " ") else NULL)
}
