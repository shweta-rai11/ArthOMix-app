#!/usr/bin/env Rscript
## Build ArthOMix transcriptomics upload example files from GEO, one dataset per call.
## Usage: Rscript build_geo_transcriptomics.R <GSE> <out_dir> <cache_dir>
## RNA-seq specs read the GEO supplementary count file from <cache_dir> (or the folder in env var GEO_SUPP_DIR).
suppressPackageStartupMessages({library(GEOquery); library(Biobase)})
a <- commandArgs(TRUE); gse <- a[1]; out <- a[2]; cache <- a[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE); dir.create(cache, recursive = TRUE, showWarnings = FALSE)
supp_dir <- Sys.getenv("GEO_SUPP_DIR", cache)

specs <- list(
  GSE110169 = list(platform = NULL, group = "disease", sex = "Sex"),
  GSE19314  = list(platform = NULL, group = "diagnosis", sex = "sex"),
  GSE17114  = list(platform = NULL, group = "affected status", sex = "gender",
                   recode = c("BD patient" = "BD", "control" = "Control")),
  GSE8650   = list(platform = "GPL96", group = "Illness", sex = "Gender", log2 = TRUE),
  GSE83456  = list(platform = NULL, group = "disease state", sex = "gender"),
  GSE39582  = list(platform = NULL, group = "chemotherapy.adjuvant", sex = "Sex",
                   recode = c(Y = "AdjuvantChemo", N = "NoChemo", "N/A" = NA)),
  GSE68465  = list(platform = NULL, group = "disease_state", sex = "Sex", log2 = TRUE,
                   recode = c("Lung Adenocarcinoma" = "Adenocarcinoma", Normal = "Normal")),
  GSE10846  = list(platform = NULL, group = "Chemotherapy", sex = "Gender",
                   recode = c("R-CHOP-Like Regimen" = "R-CHOP", "CHOP-Like Regimen" = "CHOP", "NA" = NA)),
  GSE112057 = list(platform = NULL, group = "disease state (diagnosis)", sex = "gender",
                   supp = "GSE112057_RawCounts_dataset.txt.gz",
                   recode = c("Crohn's Disease" = "CD", "Ulcerative Colitis" = "UC", "Polyarticular JIA" = "polyJIA",
                              "Oligoarticular JIA" = "oligoJIA", "Systemic JIA" = "sJIA", Control = "Control")),
  GSE178388 = list(platform = NULL, group = "diagnosis", sex = "gender",
                   supp = "GSE178388_PedBatch_all_counts_matrix.csv.gz",
                   recode = c(covid = "COVID", healthy = "Healthy", MISC = "MISC"))
)
sp <- specs[[gse]]; if (is.null(sp)) stop("no spec for ", gse)

## ---- fetch series matrix (metadata always; values for arrays) ----
gl <- suppressMessages(getGEO(gse, GSEMatrix = TRUE, destdir = cache, getGPL = is.null(sp$supp)))
eset <- if (!is.null(sp$platform)) gl[[grep(sp$platform, names(gl))[1]]] else gl[[1]]
pd <- pData(eset)

## ---- parse characteristics ourselves (robust to repeated keys like "Clinical info: X: Y") ----
chr_cols <- grep("^characteristics_ch1", colnames(pd), value = TRUE)
kv <- lapply(seq_len(nrow(pd)), function(i) {
  v <- unlist(pd[i, chr_cols]); v <- v[!is.na(v) & nzchar(trimws(v))]
  v <- sub("^Clinical info: ", "", v)
  key <- trimws(sub(":.*$", "", v)); val <- trimws(sub("^[^:]*:", "", v))
  setNames(val, key)
})
keys <- unique(unlist(lapply(kv, names)))
chars <- as.data.frame(lapply(keys, function(k) sapply(kv, function(x) unname(x[k])[1])), stringsAsFactors = FALSE)
colnames(chars) <- keys
clean <- function(x) { x <- tolower(gsub("[^A-Za-z0-9]+", "_", x)); gsub("^_|_$", "", x) }
group_raw <- chars[[sp$group]]; if (is.null(group_raw)) stop("group key not found: ", sp$group, " (have: ", paste(keys, collapse = " | "), ")")
group <- if (!is.null(sp$recode)) unname(sp$recode[group_raw]) else group_raw
sex <- toupper(substr(trimws(chars[[sp$sex]]), 1, 1)); sex[!sex %in% c("F", "M")] <- NA
extras <- chars[, setdiff(keys, c(sp$group, sp$sex)), drop = FALSE]; colnames(extras) <- clean(colnames(extras))
extras <- extras[, !duplicated(colnames(extras)), drop = FALSE]
colnames(extras)[grepl("^timepoint", colnames(extras))] <- "timepoint"
clash <- colnames(extras) %in% c("sample", "gsm", "title", "group", "sex", "group_geo", "platform", "dataset")
colnames(extras)[clash] <- paste0(colnames(extras)[clash], "_geo")

meta <- data.frame(sample = rownames(pd), gsm = rownames(pd), title = pd$title, group = group, sex = sex,
                   group_geo = group_raw, extras, platform = annotation(eset), dataset = gse,
                   check.names = FALSE, stringsAsFactors = FALSE)

## ---- expression matrix ----
maxmean_collapse <- function(mat, sym) {                              # same rule as global.R collapse_probes_to_genes()
  keep <- !is.na(sym) & sym != "" & !grepl("///", sym)
  mat <- mat[keep, , drop = FALSE]; sym <- sym[keep]
  o <- order(sym, -rowMeans(mat, na.rm = TRUE)); mat <- mat[o, , drop = FALSE]; sym <- sym[o]
  mat <- mat[!duplicated(sym), , drop = FALSE]; rownames(mat) <- sym[!duplicated(sym)]; mat
}
if (is.null(sp$supp)) {
  mat <- exprs(eset); fd <- fData(eset)
  if (nrow(mat) == 0) stop(gse, ": empty series matrix")
  sym_col <- grep("^(gene[ ._]?symbol|symbol)$", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(sym_col)) stop("no symbol column in GPL: ", paste(colnames(fd), collapse = " | "))
  n_probe <- nrow(mat)
  mat <- maxmean_collapse(mat, trimws(as.character(fd[[sym_col]])))
  if (isTRUE(sp$log2)) { mat <- log2(pmax(mat, 0) + 1); cat(gse, ": MAS5 unlogged values -> log2(x+1)\n") }
  cat(sprintf("%s: %d probes -> %d genes (symbol column '%s')\n", gse, n_probe, nrow(mat), sym_col))
} else {
  f <- file.path(supp_dir, sp$supp)
  if (gse == "GSE112057") {
    raw <- read.delim(gzfile(f), check.names = FALSE, stringsAsFactors = FALSE)
    mat <- as.matrix(raw[, -1]); n_in <- nrow(mat)
    mat <- maxmean_collapse(mat, trimws(as.character(raw[[1]])))      # deposited file repeats some symbols
    cat(sprintf("%s: %d rows -> %d unique symbols\n", gse, n_in, nrow(mat)))
    key <- sub("_.*$", "", meta$title)                                   # "Sample001_Crohn's Disease" -> Sample001
  } else if (gse == "GSE178388") {
    raw <- read.csv(gzfile(f), check.names = FALSE, stringsAsFactors = FALSE)
    mat <- as.matrix(raw[, -(1:7)]); ens <- sub("\\..*$", "", raw$Geneid)
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    sym <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, keys = ens, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
    n_ens <- nrow(mat); mat <- maxmean_collapse(mat, unname(sym))
    cat(sprintf("%s: %d Ensembl ids -> %d symbols\n", gse, n_ens, nrow(mat)))
    key <- sub("^[A-Z][0-9]+_", "", meta$title)                              # "C1_3748.1" -> 3748.1
  }
  idx <- match(key, colnames(mat))
  if (any(is.na(idx))) stop("unmatched samples: ", paste(meta$title[is.na(idx)], collapse = ", "))
  mat <- mat[, idx, drop = FALSE]; colnames(mat) <- meta$gsm
  mode(mat) <- "numeric"
}
stopifnot(identical(colnames(mat), meta$sample))

write.csv(setNames(data.frame(rownames(mat), round(mat, 5), check.names = FALSE), c("gene", colnames(mat))),
          file.path(out, paste0(gse, "_exp.csv")), row.names = FALSE)
write.csv(meta, file.path(out, paste0(gse, "_sample.csv")), row.names = FALSE, na = "NA")
cat(sprintf("%s: %d genes x %d samples; group: %s; sex: %s; range %.2f..%.2f\n", gse, nrow(mat), ncol(mat),
            paste(names(table(group, useNA = "ifany")), table(group, useNA = "ifany"), collapse = ", "),
            paste(names(table(sex, useNA = "ifany")), table(sex, useNA = "ifany"), collapse = ", "),
            min(mat, na.rm = TRUE), max(mat, na.rm = TRUE)))
