#!/usr/bin/env Rscript
## Convert any NCBI GEO series into ArthOMix upload CSVs (Feature x Sample matrix + sample metadata).
## Usage:  Rscript geo_to_arthomix_upload.R GSE89252 "clinical activity" sex [out_dir]
##         args: GSE id, group characteristic, sex characteristic (names as on the GEO sample page), output folder
## Writes: <GSE>_exp.csv (gene x sample) or <GSE>_meth.csv (cpg x sample), and <GSE>_sample.csv (sample, group, sex, ...).
## In the app: upload the matrix as "Feature and Sample (first column = feature ID)", metadata first column = sample ID.
suppressPackageStartupMessages({library(GEOquery); library(Biobase)})
a <- commandArgs(TRUE); gse <- a[1]; group_col <- a[2]; sex_col <- a[3]; out <- if (length(a) > 3) a[4] else file.path("data/uploads", gse)

eset <- getGEO(gse, GSEMatrix = TRUE, destdir = tempdir())[[1]]   # [[1]]: first platform; a SuperSeries lists several
mat <- exprs(eset); pd <- pData(eset); fd <- fData(eset)
if (nrow(mat) == 0) stop(gse, " has no values in GEO's series matrix - download its processed supplementary file instead")
is_meth <- mean(grepl("^cg[0-9]+$", rownames(mat))) > 0.5

if (is_meth) {                                                       # methylation: complete CpGs, top 20,000 by variance
  mat <- mat[grepl("^cg", rownames(mat)) & rowSums(is.na(mat)) == 0, ]
  mat <- mat[order(-apply(mat, 1, var))[seq_len(min(20000, nrow(mat)))], ]
} else {                                                             # expression: collapse probes to genes (MaxMean, as in the app)
  sym_col <- grep("^(gene[ ._]?symbol|symbol|gene_assignment)$", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(sym_col)) {
    sym <- sub("^.*? // (.*?) //.*$", "\\1", as.character(fd[[sym_col]]))
    keep <- !is.na(sym) & sym != "" & !grepl("///", sym)
    o <- order(sym[keep], -rowMeans(mat[keep, , drop = FALSE], na.rm = TRUE))
    mat <- mat[keep, ][o, ]; sym <- sym[keep][o]
    mat <- mat[!duplicated(sym), ]; rownames(mat) <- sym[!duplicated(sym)]
  }
}

ch <- function(nm) if (is.na(nm) || is.null(pd[[paste0(nm, ":ch1")]])) NA else trimws(as.character(pd[[paste0(nm, ":ch1")]]))
sex <- toupper(substr(ch(sex_col), 1, 1)); sex[!sex %in% c("F", "M")] <- NA
meta <- data.frame(sample = colnames(mat), title = pd$title, group = ch(group_col), sex = sex,
                   pd[grep(":ch1$", colnames(pd))], platform = annotation(eset), check.names = FALSE)

dir.create(out, recursive = TRUE, showWarnings = FALSE)
id <- if (is_meth) "cpg" else "gene"
write.csv(setNames(data.frame(rownames(mat), round(mat, 5), check.names = FALSE), c(id, colnames(mat))),
          file.path(out, paste0(gse, if (is_meth) "_meth.csv" else "_exp.csv")), row.names = FALSE)
write.csv(meta, file.path(out, paste0(gse, "_sample.csv")), row.names = FALSE)
cat(sprintf("%s: %d %ss x %d samples -> %s\n", gse, nrow(mat), id, ncol(mat), out))

## To use two layers from the same patients in the Multi-omics workspace, run this once per SubSeries and
## rename the sample columns of both matrices to a shared key (e.g. the GEO title), or upload as-is and
## choose "Patient ID (from metadata)" matching with `title` as the patient column.
