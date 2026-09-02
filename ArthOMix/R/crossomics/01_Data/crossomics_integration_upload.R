## R/crossomics/01_Data/crossomics_integration_upload.R
## Upload parsing for the "Expression and Methylation" Cross-Omics sub-module -
## CSV/TSV/TXT via data.table::fread, XLSX via openxlsx (both already
## installed in this deployment, called by namespace so no global.R library()
## edit is needed). Same fail-soft list(ok = FALSE, error = <message>)
## sentinel as R/methylomics/parse_upload.R, so a bad upload shows a clean
## message instead of a raw Shiny error screen.

CX_UPLOAD_NA_STRINGS <- c("NA", "", "NaN", "null", "NULL", "#N/A")

## Reads a gene/CpG-level table (not a matrix) from any of csv/tsv/txt/xlsx.
## Returns list(ok, df, error).
cx_read_table <- function(datapath, filename) {
  ext <- tolower(tools::file_ext(filename))
  if (ext %in% c("csv", "tsv", "txt")) {
    df <- tryCatch(
      as.data.frame(data.table::fread(datapath, showProgress = FALSE, na.strings = CX_UPLOAD_NA_STRINGS)),
      error = function(e) NULL
    )
    if (is.null(df) || nrow(df) == 0 || ncol(df) < 2) {
      return(list(ok = FALSE, df = NULL, error = "Could not parse this file as a table (CSV/TSV/TXT, at least one ID column plus one data column)."))
    }
    return(list(ok = TRUE, df = df, error = NULL))
  }
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      return(list(ok = FALSE, df = NULL, error = "The openxlsx package is not installed in this deployment - XLSX upload is unavailable; export your file as CSV/TSV instead."))
    }
    df <- tryCatch(openxlsx::read.xlsx(datapath, sheet = 1), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0 || ncol(df) < 2) {
      return(list(ok = FALSE, df = NULL, error = "Could not parse this XLSX file's first sheet as a table."))
    }
    return(list(ok = TRUE, df = df, error = NULL))
  }
  list(ok = FALSE, df = NULL, error = sprintf("Unsupported file type \".%s\" - upload CSV, TSV, TXT, or XLSX.", ext))
}

## Bundles the read + auto-detect steps together, since every call site needs
## both immediately after a fileInput fires.
cx_read_and_detect <- function(datapath, filename, kind = c("expression", "methylation")) {
  kind <- match.arg(kind)
  res <- cx_read_table(datapath, filename)
  if (!res$ok) return(res)
  mapping <- cx_detect_columns(res$df, kind = kind)
  list(ok = TRUE, df = res$df, mapping = mapping, error = NULL)
}
