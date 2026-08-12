## merge_and_export.R - builds a single merged, gene-level test dataset for
## the ArthOMix shiny app's "Upload your own data" section, replicating
## Chen et al. 2021/2022 (Clin Rheumatol, PMID 34767108) as closely as
## possible: all three GEO series (GSE55235, GSE55457, GSE12021) combined
## exactly as the paper's stated N implies (35 RA + 29 HC = 64 arrays),
## including the 21 GSE12021 arrays that are in fact identical to 21 of the
## GSE55457 arrays (both series were MAS5-renormalized together by the
## depositing lab for cross-study comparability, and GSE55457 is a
## reprocessing that incorporates the older GSE12021 arrays) - the paper's
## reported sample count only adds up if these are kept, not deduplicated,
## so this script keeps them to match what the paper actually analyzed.
##
## Steps: read the three probe-level exports -> log2 transform (raw MAS5
## linear-scale signal) -> quantile-normalize across all 64 arrays together
## (removes array-to-array scale differences without removing the biological
## group signal the way a full ComBat batch correction on 3 studies of
## unequal group balance risks doing) -> collapse probes to gene symbols
## with WGCNA::collapseRows (method = "MaxMean", the package's own
## recommended approach: for each gene, keep the single probe with the
## highest mean expression across arrays) -> write one merged expression
## matrix + one merged metadata file, in exactly the format
## shiny_app/R/mod_dataset.R's upload form expects (genes x samples,
## sample/group/... metadata).
suppressMessages({
  library(WGCNA)
  library(preprocessCore)
})

IN_DIR  <- "upload_csv_probelevel"
OUT_DIR <- "upload_csv_merged"
dir.create(OUT_DIR, showWarnings = FALSE)

gse_ids <- c("GSE55235", "GSE55457", "GSE12021")

read_expr <- function(id) {
  df <- read.csv(file.path(IN_DIR, paste0(id, "_probe_expression_matrix.csv")),
                  row.names = 1, check.names = FALSE)
  as.matrix(df)
}
read_meta <- function(id) {
  read.csv(file.path(IN_DIR, paste0(id, "_sample_metadata.csv")), stringsAsFactors = FALSE)
}

exprs_list <- lapply(gse_ids, read_expr)
meta_list  <- lapply(gse_ids, read_meta)
names(exprs_list) <- gse_ids

## All three exports are GPL96 only (GSE12021's GPL97 arm was already
## excluded in export.R), so the probe sets should already match exactly -
## intersect just guards against any incidental row-order/subset mismatch.
common_probes <- Reduce(intersect, lapply(exprs_list, rownames))
cat(sprintf("Common probes across all 3 series: %d\n", length(common_probes)))

merged_expr <- do.call(cbind, lapply(exprs_list, function(m) m[common_probes, , drop = FALSE]))
merged_meta <- do.call(rbind, meta_list)
stopifnot(identical(colnames(merged_expr), merged_meta$sample))
cat(sprintf("Merged: %d probes x %d samples (%d RA, %d HC)\n",
            nrow(merged_expr), ncol(merged_expr),
            sum(merged_meta$group == "RA"), sum(merged_meta$group == "HC")))

## ---- log2 transform (raw values are linear-scale MAS5 signal) ----
merged_expr[merged_expr < 1] <- 1
log_expr <- log2(merged_expr)

## ---- quantile-normalize across all 64 arrays jointly ----
qn <- preprocessCore::normalize.quantiles(log_expr)
dimnames(qn) <- dimnames(log_expr)

## ---- collapse probes -> genes: keep the max-mean probe per gene ----
annot <- read.csv(file.path(IN_DIR, "GPL96_probe_to_gene_annotation.csv"), stringsAsFactors = FALSE)
annot <- annot[match(rownames(qn), annot$probe), ]
stopifnot(identical(annot$probe, rownames(qn)))
keep_probe <- !is.na(annot$gene_symbol) & nzchar(annot$gene_symbol)
cat(sprintf("Probes with a mapped gene symbol: %d of %d\n", sum(keep_probe), nrow(qn)))

cr <- WGCNA::collapseRows(
  qn[keep_probe, , drop = FALSE],
  rowGroup = annot$gene_symbol[keep_probe],
  rowID = rownames(qn)[keep_probe],
  method = "MaxMean"
)
gene_expr <- cr$datETcollapsed
cat(sprintf("After probe->gene collapse: %d genes x %d samples\n", nrow(gene_expr), ncol(gene_expr)))

## ---- write outputs in the app's expected upload format ----
out_expr <- data.frame(gene = rownames(gene_expr), gene_expr, check.names = FALSE)
write.csv(out_expr, file.path(OUT_DIR, "chen2021_merged_expression_matrix.csv"), row.names = FALSE)
write.csv(merged_meta, file.path(OUT_DIR, "chen2021_merged_sample_metadata.csv"), row.names = FALSE)

cat("\nDone. Upload these two files in the app's \"Upload your own data\" form:\n")
cat(" -", file.path(OUT_DIR, "chen2021_merged_expression_matrix.csv"), "\n")
cat(" -", file.path(OUT_DIR, "chen2021_merged_sample_metadata.csv"), "\n")
