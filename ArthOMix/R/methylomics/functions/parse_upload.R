## R/methylomics/functions/parse_upload.R
## Upload parsing for the Methylomics Dataset tab. Never throws - each parser
## returns list(ok = FALSE, error = <message>) instead, same fail-soft pattern

methyl_parse_matrix <- function(datapath, filename) {
  df <- tryCatch(
    as.data.frame(data.table::fread(datapath, showProgress = FALSE,
                                     na.strings = c("NA", "", "NaN", "null", "NULL", "#N/A"))),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) == 0 || ncol(df) < 2) {
    return(list(ok = FALSE, error = "Could not parse this file as a probe-by-sample matrix (CSV/TSV, one probe-ID column plus at least one numeric sample column)."))
  }
  probe_ids <- as.character(df[[1]])
  if (any(duplicated(probe_ids))) {
    return(list(ok = FALSE, error = sprintf(
      "%d duplicated probe ID(s) in the first column - each row must be a unique probe.",
      sum(duplicated(probe_ids))
    )))
  }
  mat <- tryCatch({
    m <- as.matrix(df[, -1, drop = FALSE])
    storage.mode(m) <- "double"
    rownames(m) <- probe_ids
    m
  }, error = function(e) NULL)
  if (is.null(mat)) {
    return(list(ok = FALSE, error = "Every column after the first must be numeric."))
  }
  list(ok = TRUE, mat = mat)
}

methyl_parse_sample_sheet <- function(datapath, filename) {
  df <- tryCatch(as.data.frame(data.table::fread(datapath, showProgress = FALSE)), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) {
    return(list(ok = FALSE, error = "Could not parse this file as a sample sheet (CSV/TSV, one row per sample)."))
  }
  list(ok = TRUE, df = df)
}

methyl_detect_orientation <- function(mat) {
  probe_pattern <- "^(cg|ch\\.|rs)[0-9]"
  row_hits <- mean(grepl(probe_pattern, rownames(mat), ignore.case = TRUE))
  col_hits <- mean(grepl(probe_pattern, colnames(mat) %||% character(0), ignore.case = TRUE))
  list(row_hits = row_hits, col_hits = col_hits, transposed = col_hits > row_hits && col_hits > 0.5)
}

methyl_validate_matrix_upload <- function(mat, declared_scale) {
  orient <- methyl_detect_orientation(mat)
  notes <- character(0)
  if (isTRUE(orient$transposed)) {
    mat <- t(mat)
    notes <- c(notes, "Detected this matrix was oriented samples x probes (column headers looked like probe IDs, the first column like sample names) - transposed automatically to probes x samples.")
  }
  vals <- mat[is.finite(mat)]
  if (length(vals) == 0) {
    return(list(ok = FALSE, error = "No finite numeric values found in this matrix - check it isn't entirely NA/blank."))
  }
  frac_in_unit <- mean(vals >= -0.05 & vals <= 1.05)
  if (identical(declared_scale, "beta") && frac_in_unit < 0.95) {
    return(list(ok = FALSE, error = sprintf(
      "\"Beta values (0-1)\" is selected as the input scale, but %.0f%% of the values in this matrix fall outside 0-1 - this looks like M-values or a non-methylation dataset, not beta values. Switch \"Input scale\" to M-values if that's what this is, or double check this file is really methylation data.",
      100 * (1 - frac_in_unit)
    )))
  }
  if (identical(declared_scale, "m") && frac_in_unit > 0.99 && stats::sd(vals) < 0.5) {
    notes <- c(notes, "Note: these values look like they could already be beta values (all within 0-1, low spread) rather than M-values - double-check \"Input scale\" above if downstream results look off.")
  }
  list(ok = TRUE, mat = mat, note = if (length(notes) > 0) paste(notes, collapse = " ") else NULL)
}

methyl_parse_probe_list <- function(datapath, filename) {
  lines <- tryCatch(readLines(datapath, warn = FALSE), error = function(e) NULL)
  if (is.null(lines)) return(list(ok = FALSE, error = "Could not read this file.", ids = character(0)))
  first_field <- sub(",.*$|\\t.*$", "", trimws(lines))
  ids <- unique(first_field[nzchar(first_field)])
  if (length(ids) == 0) return(list(ok = FALSE, error = "No probe IDs found in this file.", ids = character(0)))
  list(ok = TRUE, ids = ids)
}

methyl_read_idat <- function(files) {
  if (!requireNamespace("minfi", quietly = TRUE)) {
    return(list(ok = FALSE, error = "The minfi package is not installed in this deployment - IDAT processing is unavailable; upload a beta or M-value matrix instead."))
  }
  is_idat <- grepl("\\.idat(\\.gz)?$", files$name, ignore.case = TRUE)
  if (!any(is_idat)) {
    return(list(ok = FALSE, error = "No .idat files found in this upload."))
  }
  files <- files[is_idat, , drop = FALSE]
  safe_names <- basename(files$name)
  has_grn <- any(grepl("_Grn\\.idat", safe_names, ignore.case = TRUE))
  has_red <- any(grepl("_Red\\.idat", safe_names, ignore.case = TRUE))
  if (!has_grn || !has_red) {
    return(list(ok = FALSE, error = "IDAT upload must include both a _Grn.idat and a _Red.idat file for every sample."))
  }

  tmp_dir <- tempfile("methyl_idat_")
  dir.create(tmp_dir, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  ok_copy <- file.copy(files$datapath, file.path(tmp_dir, safe_names), overwrite = TRUE)
  if (!all(ok_copy)) {
    return(list(ok = FALSE, error = "Could not stage one or more uploaded IDAT files for processing."))
  }

  rg <- tryCatch(minfi::read.metharray.exp(base = tmp_dir, force = TRUE), error = function(e) e)
  if (inherits(rg, "error")) {
    return(list(ok = FALSE, error = paste("Could not read the uploaded IDAT files:", conditionMessage(rg))))
  }
  list(ok = TRUE, rg = rg)
}
