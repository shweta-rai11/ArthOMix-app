## R/methylomics/parse_upload.R
## Upload parsing for the Methylomics Dataset tab (ArthOMix/R/mod_methyl_dataset.R).
## Matrix/sample-sheet parsing never throws - each returns
## list(ok = FALSE, error = <message>) instead, the same fail-soft sentinel
## pattern ArthOMix/R/mod_dataset.R and mod_preprocessing_explore.R already
## use elsewhere in this app, so a bad upload shows a clean message instead
## of a raw Shiny error screen.

## One probe-by-sample matrix (beta or M-values): CSV/TSV, first column the
## probe ID, everything else numeric.
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

## Sample sheet / phenotype metadata: CSV/TSV, one row per sample.
methyl_parse_sample_sheet <- function(datapath, filename) {
  df <- tryCatch(as.data.frame(data.table::fread(datapath, showProgress = FALSE)), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) {
    return(list(ok = FALSE, error = "Could not parse this file as a sample sheet (CSV/TSV, one row per sample)."))
  }
  list(ok = TRUE, df = df)
}

## Optional cross-reactive/custom probe-exclusion list: one probe ID per
## line (plain text, or the first column of a CSV/TSV).
methyl_parse_probe_list <- function(datapath, filename) {
  lines <- tryCatch(readLines(datapath, warn = FALSE), error = function(e) NULL)
  if (is.null(lines)) return(list(ok = FALSE, error = "Could not read this file.", ids = character(0)))
  first_field <- sub(",.*$|\\t.*$", "", trimws(lines))
  ids <- unique(first_field[nzchar(first_field)])
  if (length(ids) == 0) return(list(ok = FALSE, error = "No probe IDs found in this file.", ids = character(0)))
  list(ok = TRUE, ids = ids)
}

## `files` is a fileInput() data.frame (name/datapath/...) for a multi-file
## IDAT upload. Shiny stores each upload under a randomized tmp path, not
## its original name, but minfi pairs *_Grn.idat/*_Red.idat by matching
## basename prefix - so every file is first copied into one throwaway temp
## directory under its own (path-sanitized) original name before minfi ever
## sees it.
methyl_read_idat <- function(files) {
  if (!requireNamespace("minfi", quietly = TRUE)) {
    return(list(ok = FALSE, error = "The minfi package is not installed in this deployment - IDAT processing is unavailable; upload a beta or M-value matrix instead."))
  }
  is_idat <- grepl("\\.idat(\\.gz)?$", files$name, ignore.case = TRUE)
  if (!any(is_idat)) {
    return(list(ok = FALSE, error = "No .idat files found in this upload."))
  }
  files <- files[is_idat, , drop = FALSE]
  ## basename() strips any client-supplied path component (path-traversal
  ## protection - fileInput's reported `name` is client-controlled text).
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
