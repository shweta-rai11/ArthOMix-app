## export.R - probe-level exports + GPL96 annotation for replicating
## Chen et al. 2021/2022 (Clin Rheumatol, PMID 34767108) in ArthOMix.
suppressMessages(library(GEOquery))
OUT_DIR <- "upload_csv_probelevel"
dir.create(OUT_DIR, showWarnings = FALSE)

## ---- annotation: first-listed symbol for multi-mapped probes (closest
## available match to the paper's stated post-merge gene count) ----
gpl <- getGEO("GPL96", destdir = "raw")
tbl <- GEOquery::Table(gpl)
sym <- as.character(tbl[["Gene Symbol"]])
sym_first <- vapply(strsplit(sym, "///"), function(x) trimws(x[1]), character(1))
sym_first[is.na(sym_first) | sym_first == ""] <- NA_character_
annot <- data.frame(probe = as.character(tbl$ID), gene_symbol = sym_first, stringsAsFactors = FALSE)
write.csv(annot, file.path(OUT_DIR, "GPL96_probe_to_gene_annotation.csv"), row.names = FALSE)
cat(sprintf("Annotation: %d probes, %d unique genes.\n", nrow(annot), length(unique(na.omit(annot$gene_symbol)))))

export_probe_expr <- function(gse_id, eset, group_fn, sex_field = NULL) {
  pd  <- Biobase::pData(eset)
  grp <- group_fn(pd)
  keep <- grp %in% c("HC", "RA")
  eset_f <- eset[, keep]
  ex <- Biobase::exprs(eset_f)
  sex <- if (!is.null(sex_field) && sex_field %in% colnames(pd)) {
    s <- toupper(substr(as.character(pd[[sex_field]][keep]), 1, 1))
    ifelse(s %in% c("F", "M"), s, NA_character_)
  } else NA_character_
  meta <- data.frame(sample = colnames(ex), group = grp[keep], sex = sex, dataset = gse_id, stringsAsFactors = FALSE)
  write.csv(data.frame(probe = rownames(ex), ex, check.names = FALSE),
            file.path(OUT_DIR, paste0(gse_id, "_probe_expression_matrix.csv")), row.names = FALSE)
  write.csv(meta, file.path(OUT_DIR, paste0(gse_id, "_sample_metadata.csv")), row.names = FALSE)
  cat(sprintf("%s: %d probes x %d samples (%d HC, %d RA)\n", gse_id, nrow(ex), ncol(ex), sum(grp[keep]=="HC"), sum(grp[keep]=="RA")))
}

g1 <- getGEO("GSE55235", destdir = "raw", GSEMatrix = TRUE, AnnotGPL = TRUE)
export_probe_expr("GSE55235", g1[[1]], function(pd) {
  ds <- as.character(pd[["disease state:ch1"]])
  ifelse(grepl("healthy control", ds, ignore.case = TRUE), "HC",
  ifelse(grepl("^rheumatoid arthritis$", ds, ignore.case = TRUE), "RA", "OA"))
})

g2 <- getGEO("GSE55457", destdir = "raw", GSEMatrix = TRUE, AnnotGPL = TRUE)
export_probe_expr("GSE55457", g2[[1]], function(pd) {
  cs <- as.character(pd[["clinical status:ch1"]])
  ifelse(grepl("normal control", cs, ignore.case = TRUE), "HC",
  ifelse(grepl("rheumatoid arthritis", cs, ignore.case = TRUE), "RA", "OA"))
}, sex_field = "gender:ch1")

g3 <- getGEO("GSE12021", destdir = "raw", GSEMatrix = TRUE, AnnotGPL = TRUE)
idx3 <- which(vapply(g3, function(e) Biobase::annotation(e), character(1)) == "GPL96")
export_probe_expr("GSE12021", g3[[idx3]], function(pd) {
  t <- as.character(pd$title)
  ifelse(grepl("^Normal", t, ignore.case = TRUE), "HC",
  ifelse(grepl("^Rheumatoid arthritis", t, ignore.case = TRUE), "RA", "OA"))
}, sex_field = "sex:ch1")

cat("\nDone. Files are in", normalizePath(OUT_DIR), "\n")
