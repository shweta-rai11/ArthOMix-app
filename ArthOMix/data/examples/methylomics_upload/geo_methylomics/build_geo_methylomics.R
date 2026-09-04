#!/usr/bin/env Rscript
## =============================================================================
## geo_to_arthomix_upload.R
##
## Convert an NCBI GEO series (GSExxxxx) into the CSV files that ArthOMix's
## "Upload your own data" tabs accept (Transcriptomics, Methylomics and the
## Multi-omics Data Workspace).
##
## Output per GEO series (written to --out, default data/uploads/<GSE>):
##
##   <GSE>_exp.csv     expression: one row per gene, one column per sample.
##                     Column 1 is "gene" (probes collapsed to gene symbols with
##                     the same MaxMean rule the app uses) or "probe" when the
##                     platform table carries no symbols.
##   <GSE>_meth.csv    methylation: one row per CpG ("cpg" in column 1), one
##                     column per sample; beta or M-values exactly as deposited,
##                     restricted to the --top-cpgs most variable probes.
##   <GSE>_sample.csv  sample metadata: sample (= the matrix column names),
##                     title, group, sex (F/M/NA), every GEO characteristic,
##                     platform, dataset.
##
## Two layers from the same patients (a GEO SuperSeries) can be linked with
## --link so both matrices share one sample key and the Multi-omics workspace
## can match them with "Exact Sample ID":
##
##   <SUPER>_exp.csv, <SUPER>_meth.csv, <SUPER>_sample.csv
##
## In the app, upload each matrix with orientation
## "Feature and Sample (first column = feature ID)" and the _sample.csv as the
## sample metadata (first column = sample ID). Map Sample ID = sample,
## Group = group, Sex = sex.
##
## Usage (run from the ArthOMix app folder, or use absolute paths):
##
##   # list the sample characteristics and titles of a series first
##   Rscript data/examples/geo_to_arthomix_upload.R --gse GSE89252 --show
##
##   # one expression series; group/sex columns guessed unless given
##   Rscript data/examples/geo_to_arthomix_upload.R --gse GSE89252 \
##       --group "clinical activity" --sex sex --out data/uploads/GSE89252
##
##   # a methylation series (detected from the platform; force with --layer)
##   Rscript data/examples/geo_to_arthomix_upload.R --gse GSE89251 --top-cpgs 20000
##
##   # a series whose series-matrix has no values (RNA-seq counts, EPIC arrays):
##   # download the processed supplementary file from the GEO page first
##   Rscript data/examples/geo_to_arthomix_upload.R --gse GSE201753 \
##       --supp downloads/GSE201753_CD14_ReadCount.xlsx --group classification
##
##   # link an expression and a methylation series from the same patients
##   Rscript data/examples/geo_to_arthomix_upload.R \
##       --link GSE89252,GSE89251 --super GSE89253 --key-exp title --key-meth title
##
##   # keys that need a regex (capture groups are joined with "_")
##   Rscript data/examples/geo_to_arthomix_upload.R \
##       --link GSE32863,GSE32861 --super GSE32867 \
##       --key-exp "([0-9A-Za-z]+_[NT])" --key-meth "([0-9A-Za-z]+_[NT])"
##
##   # keys that need code: prefix with "R:"; `meta` is that layer's sample table
##   Rscript data/examples/geo_to_arthomix_upload.R \
##       --link GSE50101,GSE50222 --super GSE50387 \
##       --key-exp  "R: paste0(sub('^CS_([0-9]+[HP])_.*$', '\\1', meta$title), '_', ifelse(grepl('_DS_', meta$title), 'during', 'outside'))" \
##       --key-meth "R: paste0(sub('^([0-9]+)_.*$', '\\1', meta$title), ifelse(grepl('healthy', meta$title), 'H', 'P'), '_', ifelse(grepl('during', meta$title), 'during', 'outside'))"
##
##   # add --verify to any run to replay the app's own upload checks on the
##   # output folder, or check an existing folder on its own:
##   Rscript data/examples/geo_to_arthomix_upload.R --check data/examples/multiomics_upload/geo_multiomics
##
## Other options: --platform GPLxxxx (series with several platforms; for --link
## use --platform-exp / --platform-meth), --supp-exp / --supp-meth, --group-meth /
## --sex-meth, --cache DIR (where GEO downloads are kept, default a temp dir),
## --app-dir DIR (ArthOMix app folder for --verify/--check; default: two levels
## above this script).
##
## From an R session:
##   source("data/examples/geo_to_arthomix_upload.R")
##   lay <- geo_layer("GSE89252", group_col = "clinical activity")
##   write_layer(lay, "data/uploads/GSE89252")
##
## Requirements: GEOquery, Biobase, data.table (all in renv.lock). WGCNA is used
## for probe collapsing when installed; an equivalent MaxMean fallback runs
## otherwise. readxl is needed only for .xlsx supplementary files.
## =============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a)) b else a
msg <- function(...) cat(sprintf(...), "\n", sep = "")

script_path <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) return(normalizePath(f[1]))
  for (i in rev(seq_len(sys.nframe()))) { of <- sys.frame(i)$ofile; if (!is.null(of)) return(normalizePath(of)) }
  NA_character_
}
SCRIPT_PATH <- script_path()

## ---- 1. fetch ---------------------------------------------------------------

#' Download one GEO series (values + sample table, no platform table) and return
#' its ExpressionSet. A series spanning several arrays returns several
#' platforms: pick one with `platform`, otherwise they are listed and we stop.
geo_fetch <- function(gse, platform = NULL, cache = tempdir()) {
  gse <- toupper(trimws(gse))
  stopifnot(grepl("^GSE[0-9]+$", gse))
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  msg("Fetching %s from GEO (cache: %s) ...", gse, cache)
  gset <- suppressMessages(getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE, destdir = cache))
  plats <- vapply(gset, annotation, character(1))
  if (length(gset) > 1) {
    if (is.null(platform) || !platform %in% plats) {
      msg("%s spans %d platforms: %s", gse, length(gset), paste(plats, collapse = ", "))
      msg("Re-run with --platform <GPL>. (A SuperSeries lists both its expression and its methylation platform here; convert each SubSeries separately, then --link them.)")
      stop("platform not chosen")
    }
    eset <- gset[[which(plats == platform)[1]]]
  } else {
    eset <- gset[[1]]
  }
  msg("  %s: %s features x %s samples on %s", gse, format(nrow(eset), big.mark = ","), ncol(eset), annotation(eset))
  eset
}

#' Platform annotation table (GEO's curated .annot file when it exists, else
#' the submitter's table) - the same source the app's GEO fetch uses.
geo_platform_table <- function(gpl, cache = tempdir()) {
  g <- tryCatch(suppressMessages(getGEO(gpl, destdir = cache, AnnotGPL = TRUE)), error = function(e) NULL)
  if (is.null(g)) g <- tryCatch(suppressMessages(getGEO(gpl, destdir = cache, AnnotGPL = FALSE)), error = function(e) NULL)
  if (is.null(g)) return(NULL)
  tb <- Table(g)
  if (!"ID" %in% colnames(tb)) return(NULL)
  tb
}

## ---- 2. sample metadata -----------------------------------------------------

norm_sex <- function(x) {
  s <- toupper(substr(gsub("[^A-Za-z]", "", as.character(x)), 1, 1))
  s[!s %in% c("F", "M")] <- NA_character_
  s
}
clean_name <- function(x) gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(sub(":ch1$", "", x))))
clean_title <- function(x) gsub("\\s+", "_", trimws(sub("\\s*[\\[(][^\\])]*[\\])]$", "", trimws(as.character(x)), perl = TRUE)))

#' Build the _sample.csv table from an ExpressionSet's phenoData.
#' `group_col` / `sex_col` are GEO characteristic names (as printed by --show);
#' when NULL they are guessed and the guess is printed.
geo_sample_meta <- function(eset, gse, group_col = NULL, sex_col = NULL) {
  pd <- pData(eset)
  ch <- grep(":ch1$", colnames(pd), value = TRUE)
  chars <- if (length(ch)) as.data.frame(lapply(pd[ch], function(v) trimws(as.character(v))), stringsAsFactors = FALSE, check.names = FALSE)
           else data.frame(row.names = seq_len(nrow(pd)))
  colnames(chars) <- clean_name(ch)
  find_col <- function(spec, rx, label) {
    if (!is.null(spec)) {
      hit <- colnames(chars)[colnames(chars) == clean_name(spec)]
      if (length(hit) == 0) stop(sprintf("%s column '%s' not found. Available: %s", label, spec, paste(colnames(chars), collapse = ", ")))
      return(hit[1])
    }
    hit <- grep(rx, colnames(chars), value = TRUE)
    if (length(hit) == 0) { msg("  !! no %s column guessed - pass --%s <characteristic> (see --show)", label, tolower(label)); return(NA_character_) }
    msg("  %s column guessed: '%s' (override with --%s)", label, hit[1], tolower(label))
    hit[1]
  }
  gcol <- find_col(group_col, "^(group|disease|status|condition|diagnosis|phenotype|treatment|tissue|clinical|subtype|classification)", "Group")
  scol <- find_col(sex_col, "^(sex|gender)$", "Sex")
  meta <- data.frame(sample = rownames(pd), title = as.character(pd$title),
                     group = if (is.na(gcol)) NA_character_ else chars[[gcol]],
                     sex = if (is.na(scol)) NA_character_ else norm_sex(chars[[scol]]),
                     stringsAsFactors = FALSE, check.names = FALSE)
  extra <- chars[, setdiff(colnames(chars), c("sample", "title", "group", "sex")), drop = FALSE]
  meta <- cbind(meta, extra)
  meta$platform <- annotation(eset); meta$dataset <- gse
  if (!is.na(gcol)) msg("  group: %s", paste(names(table(meta$group)), table(meta$group), sep = "=", collapse = ", "))
  msg("  sex: %s", paste(names(table(meta$sex, useNA = "ifany")), table(meta$sex, useNA = "ifany"), sep = "=", collapse = ", "))
  rownames(meta) <- NULL
  meta
}

## ---- 3. expression: probes -> gene symbols ----------------------------------

#' Find a gene-symbol vector in a platform table. Broader than the app's own
#' search (which only accepts "gene symbol"): Illumina uses "Symbol", Agilent
#' "GENE_SYMBOL", Affymetrix HTA/Clariom a "gene_assignment" string.
geo_probe_symbols <- function(fd) {
  col <- grep("^(gene[ ._]?symbol|symbol|ilmn_gene|gene|orf)$", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(col)) return(as.character(fd[[col]]))
  ga <- grep("^gene_assignment$", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(ga)) {
    first <- sub(" ///.*$", "", as.character(fd[[ga]]))
    parts <- strsplit(first, " // ", fixed = TRUE)
    return(vapply(parts, function(p) if (length(p) >= 2) trimws(p[2]) else NA_character_, character(1)))
  }
  NULL
}

#' Collapse a probe x sample matrix to gene x sample with the MaxMean rule (the
#' probe with the highest mean represents the gene) - identical to
#' collapse_probes_to_genes() in the app's global.R.
collapse_to_genes <- function(ex, sym) {
  keep <- !is.na(sym) & sym != "" & !grepl("///", sym)
  ex <- ex[keep, , drop = FALSE]; sym <- sym[keep]
  ok <- rowSums(is.na(ex)) < ncol(ex)
  ex <- ex[ok, , drop = FALSE]; sym <- sym[ok]
  if (nrow(ex) == 0) return(NULL)
  if (requireNamespace("WGCNA", quietly = TRUE)) {
    out <- WGCNA::collapseRows(ex, rowGroup = sym, rowID = rownames(ex), method = "MaxMean",
                               connectivityBasedCollapsing = FALSE, connectivityPower = 1,
                               selectFewestMissing = TRUE, thresholdCombine = NA)$datETcollapsed
  } else {
    mm <- rowMeans(ex, na.rm = TRUE); o <- order(sym, -mm)
    ex <- ex[o, , drop = FALSE]; sym <- sym[o]
    out <- ex[!duplicated(sym), , drop = FALSE]; rownames(out) <- sym[!duplicated(sym)]
  }
  out[order(rownames(out)), , drop = FALSE]
}

looks_like_probe_ids <- function(ids) mean(grepl("^[0-9]+_(at|st)$|^ILMN_|^A_[0-9]+_P|^TC[0-9]{2}|^[0-9]{6,}$", ids)) > 0.5

## ---- 4. supplementary processed matrices ------------------------------------

#' Read a processed matrix that GEO ships only as a supplementary file
#' (RNA-seq count tables, EPIC "processed_data" tables) and rename its columns
#' to GSM ids using the series' sample table: exact GSM, sample title, or the
#' Sentrix id in the IDAT file name. Column names are also tried with a
#' trailing "_token" removed and reduced to their last "_token".
geo_supp_matrix <- function(path, eset) {
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", path)))
  dt <- if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("readxl is needed for .xlsx supplementary files")
    as.data.table(readxl::read_excel(path, sheet = 1))
  } else {
    fread(path, showProgress = FALSE)
  }
  cn <- colnames(dt)
  keep <- seq_along(cn) > 1 & !grepl("p_?val|pval|detection", cn, ignore.case = TRUE)
  m <- as.matrix(dt[, keep, with = FALSE]); rownames(m) <- as.character(dt[[1]])
  suppressWarnings(storage.mode(m) <- "double")
  pd <- pData(eset); gsm <- rownames(pd); title <- trimws(as.character(pd$title))
  supp_cols <- grep("^supplementary_file", colnames(pd), value = TRUE)
  supp <- if (length(supp_cols)) do.call(paste, c(lapply(pd[supp_cols], as.character), sep = " ")) else rep("", length(gsm))
  sentrix <- rep(NA_character_, length(gsm))
  has <- grepl("[0-9]{10,}_R[0-9]{2}C[0-9]{2}", supp)
  sentrix[has] <- regmatches(supp, regexpr("[0-9]{10,}_R[0-9]{2}C[0-9]{2}", supp))
  keys <- list(gsm = gsm, title = title, title_last = sub("^.*_", "", title), sentrix = sentrix)
  cols <- colnames(m)
  variants <- list(as_is = cols, no_suffix = sub("_[^_]+$", "", cols), last_token = sub("^.*_", "", cols))
  for (kn in names(keys)) for (vn in names(variants)) {
    hit <- match(variants[[vn]], keys[[kn]])
    if (mean(!is.na(hit)) > 0.9 && !any(duplicated(hit[!is.na(hit)]))) {
      msg("  supplementary columns (%s) mapped to GSM via %s: %d of %d", vn, kn, sum(!is.na(hit)), length(cols))
      m <- m[, !is.na(hit), drop = FALSE]; colnames(m) <- gsm[hit[!is.na(hit)]]
      return(m)
    }
  }
  stop(sprintf("Could not map the supplementary file's columns to GSM ids. First columns: %s; first titles: %s. Rename the columns to GSM ids or sample titles and retry.",
               paste(head(cols, 4), collapse = ", "), paste(head(title, 4), collapse = ", ")))
}

## ---- 5. one layer -----------------------------------------------------------

METH_PLATFORMS <- c("GPL8490", "GPL13534", "GPL16304", "GPL21145", "GPL23976", "GPL33022")

#' Fetch one GEO series and turn it into an upload-ready layer:
#' list(kind = "expression"|"methylation", mat (feature x GSM), meta, gse, note).
geo_layer <- function(gse, layer = c("auto", "expression", "methylation"), platform = NULL,
                      group_col = NULL, sex_col = NULL, supp = NULL, top_cpgs = 20000, cache = tempdir()) {
  layer <- match.arg(layer)
  eset <- geo_fetch(gse, platform, cache)
  ex <- exprs(eset)
  if (!is.null(supp)) {
    msg("  reading supplementary matrix %s", supp)
    ex <- geo_supp_matrix(supp, eset)
    eset <- eset[, colnames(ex)]
  } else if (nrow(ex) == 0 || all(is.na(ex))) {
    stop(sprintf("%s has no values in its GEO series matrix. Download the processed supplementary file from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s, decompress nothing (gz is fine), and pass it with --supp.", gse, gse))
  }
  storage.mode(ex) <- "double"
  if (layer == "auto") {
    is_cpg <- mean(grepl("^cg[0-9]+$", rownames(ex))) > 0.5
    layer <- if (annotation(eset) %in% METH_PLATFORMS || is_cpg) "methylation" else "expression"
    msg("  layer detected: %s", layer)
  }
  meta <- geo_sample_meta(eset, gse, group_col, sex_col)
  note <- character(0)
  if (layer == "expression") {
    genes <- NULL
    if (is.null(supp) || looks_like_probe_ids(rownames(ex))) {
      tb <- geo_platform_table(annotation(eset), cache)
      if (!is.null(tb)) {
        sym <- geo_probe_symbols(tb)
        if (!is.null(sym)) genes <- collapse_to_genes(ex, sym[match(rownames(ex), as.character(tb$ID))])
      }
    }
    if (!is.null(genes)) {
      ex <- genes; msg("  probes collapsed to %s gene symbols (MaxMean)", format(nrow(ex), big.mark = ","))
    } else if (looks_like_probe_ids(rownames(ex))) {
      note <- c(note, "platform table has no gene symbols - file left at probe level (column 1 = probe); map probes to symbols with a Bioconductor annotation package if you need genes")
      msg("  !! %s", note)
    }
    ex <- ex[rowSums(is.na(ex)) < ncol(ex), , drop = FALSE]
    r <- range(ex, na.rm = TRUE)
    is_int <- all(ex[is.finite(ex)] == round(ex[is.finite(ex)]))
    scale <- if (r[1] >= 0 && r[2] > 1000 && is_int) "Raw counts"
      else if (r[2] <= 30) "Already log-transformed (or Normalized) - log2 scale" else "Normalized (linear scale, not log)"
    n0 <- sum(rowSums(ex != 0, na.rm = TRUE) == 0)
    msg("  value range %.2f .. %.2f -> declare as: %s%s", r[1], r[2], scale, if (n0) sprintf(" (%s all-zero genes kept; the app filters them)", format(n0, big.mark = ",")) else "")
    note <- c(note, sprintf("declare the data type as: %s", scale))
  } else {
    ex <- ex[grepl("^cg[0-9]+$", rownames(ex)), , drop = FALSE]
    na_frac <- rowMeans(is.na(ex))
    ex <- ex[na_frac <= 0.1, , drop = FALSE]
    r <- range(ex, na.rm = TRUE)
    vt <- if (r[1] >= -0.001 && r[2] <= 1.001) "beta" else if (r[1] < 0 && abs(r[2]) < 20) "M-value" else "unrecognised scale"
    msg("  %s CpGs (<=10%% missing), range %.3f .. %.3f -> %s values", format(nrow(ex), big.mark = ","), r[1], r[2], vt)
    if (top_cpgs > 0 && nrow(ex) > top_cpgs) {
      v <- apply(ex, 1, var, na.rm = TRUE); ex <- ex[order(-v)[seq_len(top_cpgs)], , drop = FALSE]
      msg("  kept the %s most variable CpGs (--top-cpgs 0 keeps all)", format(top_cpgs, big.mark = ","))
    }
    note <- c(note, if (vt == "beta") "beta values (0-1): Methylomics tab input scale = beta; Multi-omics normalization = 'Beta values (as-is)' or 'M-value transform'"
                    else "M-values (already logit-transformed): Methylomics tab input scale = M-value; Multi-omics normalization = 'Beta values (as-is)' - do NOT apply the M-value transform again")
  }
  stopifnot(identical(colnames(ex), meta$sample))
  list(kind = layer, mat = ex, meta = meta, gse = gse, note = note)
}

write_matrix_csv <- function(mat, id_col, path) {
  dt <- data.table(id = rownames(mat), as.data.table(round(mat, 5)))
  setnames(dt, "id", id_col)
  fwrite(dt, path)
}

#' Write <gse>_exp.csv / <gse>_meth.csv and <gse>_sample.csv.
write_layer <- function(lay, out, name = lay$gse) {
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  id_col <- if (lay$kind == "methylation") "cpg" else if (looks_like_probe_ids(rownames(lay$mat))) "probe" else "gene"
  suffix <- if (lay$kind == "methylation") "_meth.csv" else "_exp.csv"
  write_matrix_csv(lay$mat, id_col, file.path(out, paste0(name, suffix)))
  fwrite(lay$meta, file.path(out, paste0(name, "_sample.csv")))
  msg("  wrote %s (%s x %s) and %s to %s", paste0(name, suffix), format(nrow(lay$mat), big.mark = ","), ncol(lay$mat), paste0(name, "_sample.csv"), out)
  for (n in lay$note) msg("  note: %s", n)
  invisible(file.path(out, paste0(name, c(suffix, "_sample.csv"))))
}

## ---- 6. link two layers on a shared patient key -----------------------------

#' Turn a key spec into a function(meta) -> character keys.
#'   "title"           cleaned title ("13106 T0" -> "13106_T0"; a trailing
#'                     "[methylation]" / "(expression)" tag is dropped)
#'   "<column>"        a column of the sample table (e.g. "subject", "donor")
#'   "R: <expression>" R code evaluated with `meta` in scope
#'   anything else     a regex applied to title; capture groups joined with "_"
make_key_fun <- function(spec) {
  if (is.function(spec)) return(spec)
  if (identical(spec, "title")) return(function(meta) clean_title(meta$title))
  if (grepl("^R:", spec)) { expr <- parse(text = sub("^R:\\s*", "", spec)); return(function(meta) as.character(eval(expr, list(meta = meta)))) }
  function(meta) {
    if (spec %in% colnames(meta)) return(as.character(meta[[spec]]))
    m <- regmatches(meta$title, regexec(spec, meta$title))
    vapply(m, function(g) if (length(g) >= 2) paste(g[-1], collapse = "_") else if (length(g) == 1) g else NA_character_, character(1))
  }
}

#' Rename both layers' samples to a shared key and write <super>_exp.csv,
#' <super>_meth.csv, <super>_sample.csv. Samples present in only one layer are
#' kept (the app lists them as layer-only) and flagged in the sample table.
link_layers <- function(exp_lay, meth_lay, super, out, key_exp = "title", key_meth = "title") {
  kx <- make_key_fun(key_exp)(exp_lay$meta); km <- make_key_fun(key_meth)(meth_lay$meta)
  check_key <- function(k, lay, flag) {
    empty <- is.na(k) | k == ""
    if (any(empty)) stop(sprintf("%s: --key-%s gives no key for some samples, e.g. title '%s'", lay$gse, flag, lay$meta$title[empty][1]))
    if (any(duplicated(k))) stop(sprintf("%s: --key-%s is not unique (e.g. %s) - refine it, or include the time point / tissue in the key", lay$gse, flag, paste(head(unique(k[duplicated(k)]), 3), collapse = ", ")))
  }
  check_key(kx, exp_lay, "exp"); check_key(km, meth_lay, "meth")
  shared <- intersect(kx, km)
  msg("Linking %s (%d samples) with %s (%d samples): %d shared keys, %d expression-only, %d methylation-only",
      exp_lay$gse, length(kx), meth_lay$gse, length(km), length(shared), length(setdiff(kx, km)), length(setdiff(km, kx)))
  if (length(setdiff(kx, km))) msg("  expression-only: %s", paste(head(setdiff(kx, km), 8), collapse = ", "))
  if (length(setdiff(km, kx))) msg("  methylation-only: %s", paste(head(setdiff(km, kx), 8), collapse = ", "))
  if (length(shared) < 3) stop("fewer than 3 shared samples - the key specs disagree (use --show to print both layers' titles)")
  ex <- exp_lay$mat; colnames(ex) <- kx
  me <- meth_lay$mat; colnames(me) <- km
  mx <- exp_lay$meta; mx$gsm_exp <- mx$sample; mx$sample <- kx
  mm <- meth_lay$meta; mm$gsm_meth <- mm$sample; mm$sample <- km
  mm$group_meth <- mm$group; mm$sex_meth <- mm$sex
  mm <- mm[, c("sample", "gsm_meth", "group_meth", "sex_meth", setdiff(colnames(mm), c(colnames(mx), "gsm_meth", "group_meth", "sex_meth"))), drop = FALSE]
  meta <- merge(mx, mm, by = "sample", all = TRUE)
  for (f in c("group", "sex")) {   # fill from the methylation layer when the expression layer lacks it
    alt <- paste0(f, "_meth"); na <- is.na(meta[[f]]) & !is.na(meta[[alt]])
    meta[[f]][na] <- meta[[alt]][na]
    both <- !is.na(meta[[f]]) & !is.na(meta[[alt]])
    if (all(meta[[f]][both] == meta[[alt]][both])) meta[[alt]] <- NULL else msg("  !! %s differs between layers for %d samples - kept both columns", f, sum(meta[[f]][both] != meta[[alt]][both]))
  }
  meta$in_expression <- meta$sample %in% kx; meta$in_methylation <- meta$sample %in% km
  meta$superseries <- super
  first <- c("sample", "gsm_exp", "gsm_meth", "group", "sex")
  meta <- meta[order(!meta$in_expression, !meta$in_methylation, meta$sample), c(first, setdiff(colnames(meta), first)), drop = FALSE]
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  write_matrix_csv(ex, if (looks_like_probe_ids(rownames(ex))) "probe" else "gene", file.path(out, paste0(super, "_exp.csv")))
  write_matrix_csv(me, "cpg", file.path(out, paste0(super, "_meth.csv")))
  fwrite(meta, file.path(out, paste0(super, "_sample.csv")))
  msg("  wrote %s_exp.csv, %s_meth.csv, %s_sample.csv (%d metadata rows) to %s", super, super, super, nrow(meta), out)
  for (n in c(exp_lay$note, meth_lay$note)) msg("  note: %s", n)
  invisible(meta)
}

## ---- 7. replay the app's own upload checks on a folder ----------------------

#' Source the app's upload helpers and run the checks the Multi-omics Data
#' Workspace runs on upload: read (feature x sample), omics-type detection,
#' matrix validation/status, cross-layer sample overlap, metadata coverage.
arthomix_check_upload <- function(dir, app_dir = NULL) {
  app_dir <- app_dir %||% if (!is.na(SCRIPT_PATH)) normalizePath(file.path(dirname(SCRIPT_PATH), "..", "..")) else getwd()
  helpers <- file.path(app_dir, c("R/multiomics/01_Data_Workspace/multiomics_dataset_helpers.R",
                                  "R/multiomics/06_Gene_CpG_Mapping/multiomics_mapping_helpers.R",
                                  "R/multiomics/02_Cohort_Harmonization/cohort_harmonization_helpers.R"))
  if (!all(file.exists(helpers))) stop("app helpers not found under ", app_dir, " - pass --app-dir <path to the ArthOMix app folder>")
  env <- new.env(); env$`%||%` <- `%||%`
  for (h in helpers) sys.source(h, envir = env)
  files <- list.files(dir, pattern = "_(exp|meth)\\.csv$", full.names = TRUE)
  if (length(files) == 0) stop("no *_exp.csv / *_meth.csv in ", dir)
  ok_all <- TRUE
  for (prefix in unique(sub("_(exp|meth)\\.csv$", "", basename(files)))) {
    msg("\n== %s", prefix)
    mats <- list()
    for (kind in c("exp", "meth")) {
      f <- file.path(dir, paste0(prefix, "_", kind, ".csv")); if (!file.exists(f)) next
      r <- env$multi_live_read_matrix(f, orientation = "features_rows", filename = f)
      if (!isTRUE(r$ok)) { msg("  %s: READ FAILED - %s", basename(f), r$error); ok_all <- FALSE; next }
      det <- env$multi_live_detect_omics_type(r$mat)
      want <- if (kind == "exp") "rnaseq" else "methylation"
      v <- env$multi_live_validate_matrix(r$mat, basename(f)); st <- env$multi_dataset_status(v)
      msg("  %s: %s samples x %s features; detected %s (%s); status %s%s", basename(f), nrow(r$mat), format(ncol(r$mat), big.mark = ","),
          det$detected, det$reason, st$label, if (length(st$reasons)) paste0(" - ", paste(st$reasons, collapse = "; ")) else "")
      if (!identical(det$detected, want)) { msg("  !! the app would REJECT this file as %s", want); ok_all <- FALSE }
      if (r$n_coerced_na > 0) msg("  !! %d values could not be read as numbers", r$n_coerced_na)
      mats[[kind]] <- r$mat
    }
    sf <- file.path(dir, paste0(prefix, "_sample.csv"))
    if (file.exists(sf)) {
      meta <- as.data.frame(fread(sf)); ids <- as.character(meta[[1]])
      if (any(duplicated(ids))) { msg("  !! duplicated sample ids in %s", basename(sf)); ok_all <- FALSE }
      for (kind in names(mats)) {
        miss <- setdiff(rownames(mats[[kind]]), ids)
        if (length(miss)) { msg("  !! %d %s samples missing from %s: %s", length(miss), kind, basename(sf), paste(head(miss, 5), collapse = ", ")); ok_all <- FALSE }
      }
      msg("  %s: %d rows; columns: %s", basename(sf), nrow(meta), paste(head(colnames(meta), 8), collapse = ", "))
      if ("group" %in% colnames(meta)) msg("  group: %s", paste(names(table(meta$group, useNA = "ifany")), table(meta$group, useNA = "ifany"), sep = "=", collapse = ", "))
      if ("sex" %in% colnames(meta)) msg("  sex: %s", paste(names(table(meta$sex, useNA = "ifany")), table(meta$sex, useNA = "ifany"), sep = "=", collapse = ", "))
    } else msg("  (no %s_sample.csv)", prefix)
    if (length(mats) == 2) {
      ov <- env$multi_live_sample_overlap(mats)
      cmp <- env$multi_dataset_compatibility(lapply(mats, function(m) env$multi_live_validate_matrix(m)), ov, has_metadata = file.exists(sf))
      msg("  exact-ID overlap between layers: %d shared -> %s", ov$n_shared, cmp$overall_label)
      if (ov$n_shared < 3) ok_all <- FALSE
    }
  }
  msg("\n%s", if (ok_all) "ALL CHECKS PASSED" else "SOME CHECKS FAILED (see the !! lines)")
  invisible(ok_all)
}

## ---- 8. command line --------------------------------------------------------

FLAGS <- c("show", "verify")

parse_args <- function(args) {
  out <- list(); i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (!grepl("^--", a)) stop("unexpected argument: ", a)
    key <- sub("^--", "", a)
    if (key %in% FLAGS) { out[[key]] <- TRUE; i <- i + 1; next }
    if (i == length(args)) stop("missing value for --", key)
    out[[key]] <- args[i + 1]; i <- i + 2
  }
  out
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    cat(paste(grep("^## ", readLines(SCRIPT_PATH, n = 80), value = TRUE), collapse = "\n"), "\n")
    return(invisible())
  }
  o <- parse_args(args)
  cache <- o$cache %||% file.path(tempdir(), "geo_cache")
  top <- as.integer(o[["top-cpgs"]] %||% 20000)
  if (!is.null(o$check)) return(invisible(arthomix_check_upload(o$check, o[["app-dir"]])))
  if (!is.null(o$link)) {
    ids <- toupper(trimws(strsplit(o$link, ",")[[1]])); stopifnot(length(ids) == 2)
    super <- o$super %||% paste(ids, collapse = "_")
    out <- o$out %||% file.path("data", "uploads", super)
    ex <- geo_layer(ids[1], "expression", o[["platform-exp"]], o$group, o$sex, o[["supp-exp"]], top, cache)
    me <- geo_layer(ids[2], "methylation", o[["platform-meth"]], o[["group-meth"]], o[["sex-meth"]], o[["supp-meth"]], top, cache)
    if (isTRUE(o$show)) {
      msg("\nExpression titles e.g. %s", paste(head(ex$meta$title, 6), collapse = " | "))
      msg("Methylation titles e.g. %s", paste(head(me$meta$title, 6), collapse = " | "))
    }
    link_layers(ex, me, super, out, o[["key-exp"]] %||% "title", o[["key-meth"]] %||% "title")
  } else {
    if (is.null(o$gse)) stop("--gse GSExxxxx, --link GSEexp,GSEmeth or --check DIR is required (--help for usage)")
    out <- o$out %||% file.path("data", "uploads", toupper(o$gse))
    if (isTRUE(o$show)) {
      eset <- geo_fetch(o$gse, o$platform, cache); pd <- pData(eset)
      msg("\nSample characteristics available for --group / --sex:")
      for (ch in grep(":ch1$", colnames(pd), value = TRUE)) msg("  %-32s e.g. %s", clean_name(ch), paste(head(unique(trimws(as.character(pd[[ch]]))), 4), collapse = " | "))
      msg("Titles e.g. %s", paste(head(pd$title, 5), collapse = " | "))
      return(invisible())
    }
    lay <- geo_layer(o$gse, o$layer %||% "auto", o$platform, o$group, o$sex, o$supp, top, cache)
    write_layer(lay, out)
  }
  if (isTRUE(o$verify)) arthomix_check_upload(out, o[["app-dir"]])
  invisible()
}

if (sys.nframe() == 0L && !interactive()) main()
