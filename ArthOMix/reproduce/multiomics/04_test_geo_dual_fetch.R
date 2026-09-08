## Verifies the new multi_geo_dual_fetch() + geo_title_patient_num fix against
## the real-world case that motivated it: GSE138746 (RNA-seq) + GSE138653
## (methylation), Tao et al. 2021 - two independently-submitted, UNLINKED GEO
## series (multi_geo_autosplit_fetch() cannot find them; GEO's own
## series-relation metadata lists zero linked sub-series for either).
##
## Checks, in order:
##  1. multi_geo_autosplit_fetch() on either accession alone still correctly
##     fails with the documented "doesn't look like a SuperSeries" message
##     (confirms this really is the unlinked case, not a regression).
##  2. multi_geo_dual_fetch() succeeds fetching both, REGARDLESS of which
##     accession is passed first (order-independence).
##  3. geo_title_patient_num is correctly derived on both sides.
##  4. mo_apply_matching(..., "patient_id", patient_col="geo_title_patient_num")
##     - the exact function the Shiny "Sample Matching" tab calls - actually
##     matches samples between the two independently-fetched matrices.
##
## Run from the ArthOMix/ app directory:
##   Rscript reproduce/multiomics/04_test_geo_dual_fetch.R

suppressMessages(suppressWarnings(
  shiny::loadSupport(".", renv = globalenv(), globalrenv = globalenv())
))

cat("=== Step 1: confirm autosplit still correctly fails for the unlinked case ===\n")
r1 <- multi_geo_autosplit_fetch("GSE138746")
cat(sprintf("autosplit(GSE138746): ok=%s  %s\n", r1$ok, if (!isTRUE(r1$ok)) r1$error else ""))
stopifnot(!isTRUE(r1$ok))

cat("\n=== Step 2: dual_fetch succeeds, order-independent ===\n")
r2a <- multi_geo_dual_fetch("GSE138746", "GSE138653")
r2b <- multi_geo_dual_fetch("GSE138653", "GSE138746")
cat(sprintf("dual_fetch(746, 653): ok=%s  expr=%s (%d x %d)  meth=%s (%d x %d)\n",
            r2a$ok, r2a$expression$accession, nrow(r2a$expression$mat), ncol(r2a$expression$mat),
            r2a$methylation$accession, nrow(r2a$methylation$mat), ncol(r2a$methylation$mat)))
cat(sprintf("dual_fetch(653, 746): ok=%s  expr=%s (%d x %d)  meth=%s (%d x %d)\n",
            r2b$ok, r2b$expression$accession, nrow(r2b$expression$mat), ncol(r2b$expression$mat),
            r2b$methylation$accession, nrow(r2b$methylation$mat), ncol(r2b$methylation$mat)))
stopifnot(isTRUE(r2a$ok), isTRUE(r2b$ok))
stopifnot(identical(r2a$expression$accession, r2b$expression$accession))
stopifnot(identical(r2a$methylation$accession, r2b$methylation$accession))

cat("\n=== Step 3: geo_title_patient_num derived correctly ===\n")
n_expr_id <- sum(!is.na(r2a$expression$meta$geo_title_patient_num))
n_meth_id <- sum(!is.na(r2a$methylation$meta$geo_title_patient_num))
cat(sprintf("Expression: %d / %d samples got a derived patient number\n", n_expr_id, nrow(r2a$expression$mat)))
cat(sprintf("Methylation: %d / %d samples got a derived patient number\n", n_meth_id, nrow(r2a$methylation$mat)))
stopifnot(n_expr_id > 0, n_meth_id > 0)

cat("\n=== Step 4: mo_apply_matching() - the exact function the Shiny UI calls ===\n")
## Restrict expression to PBMC only (the same cell type as methylation) - a real
## user would do this via the app's own filtering; done directly here since
## this script isn't driving the UI.
pbmc_ids <- rownames(r2a$expression$meta)[r2a$expression$meta[["cell type:ch1"]] == "PBMC"]
expr_pbmc <- r2a$expression$mat[intersect(pbmc_ids, rownames(r2a$expression$mat)), , drop = FALSE]

## raw$meta in the live app is the row-stacked union of both series' metadata
## (mo_merge_sample_meta), which is exactly what mo_apply_matching() expects:
## each GSM row keeps its own series' column values, rownamed by GSM.
combined_meta <- rbind(
  cbind(r2a$expression$meta[rownames(expr_pbmc), c("title", "geo_title_patient_num"), drop = FALSE]),
  cbind(r2a$methylation$meta[, c("title", "geo_title_patient_num"), drop = FALSE])
)

mats <- list(Transcriptomics = expr_pbmc, Methylomics = r2a$methylation$mat)
matched <- mo_apply_matching(mats, "patient_id", meta = combined_meta, patient_col = "geo_title_patient_num")
n_expr_matched <- nrow(matched$mats$Transcriptomics)
n_meth_matched <- nrow(matched$mats$Methylomics)
n_common <- length(intersect(rownames(matched$mats$Transcriptomics), rownames(matched$mats$Methylomics)))
cat(sprintf("After matching: Transcriptomics %d samples, Methylomics %d samples, %d with a shared patient ID (usable pair)\n",
            n_expr_matched, n_meth_matched, n_common))
stopifnot(n_common >= 70)

cat("\n=== ALL CHECKS PASSED: the GEO auto-fetch fix works end-to-end for a real unlinked two-series dataset ===\n")
