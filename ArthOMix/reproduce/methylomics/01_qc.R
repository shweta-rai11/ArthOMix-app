## Regenerates the methylomics QC stage (PCA + chrY sex-check) from the
## bundled raw beta matrix, and compares it against the precomputed
## sample_metadata_qc.csv / qc_pca_sexcheck.csv tables. See
## reproduce/methylomics/README.md for the reproducibility framing.
##
## Run from the ArthOMix/ app directory:
##   Rscript reproduce/methylomics/01_qc.R

## Replicates exactly how shiny::runApp() bootstraps this app: global.R, then
## every file under R/ (loadSupport() is what Shiny calls internally before
## ui.R/server.R run) - so every module/helper function is available exactly
## as it is in the live app, with no manual re-sourcing/order guessing.
suppressMessages(suppressWarnings(
  shiny::loadSupport(".", renv = globalenv(), globalrenv = globalenv())
))

beta <- readRDS("data/preloaded/methylomics/matrix/beta_raw.rds")
pheno <- as.data.frame(readRDS("data/preloaded/methylomics/matrix/pheno.rds"))
rownames(pheno) <- pheno$gsm
cat(sprintf("beta_raw.rds: %d probes x %d samples\n", nrow(beta), ncol(beta)))

anno <- methyl_get_annotation("450K")
stopifnot(isTRUE(anno$ok))

## GENUINE GAP, not a bug: beta_raw.rds as bundled has already had every
## chrX/chrY probe excluded (confirmed: zero probes annotate to chrX/chrY out
## of 412,492) - the chrY-based sex-check that produced chrY_mean/sex_mismatch
## in the precomputed table must have run on an earlier, non-bundled version
## of the matrix (pre sex-chromosome filtering). It cannot be reproduced from
## what's actually in this repo, and is reported here as NOT-REPRODUCIBLE
## rather than silently skipped.
x_probes <- sum(!is.na(anno$anno[rownames(beta), "chr"]) & anno$anno[rownames(beta), "chr"] %in% c("chrX", "X"))
y_probes <- sum(!is.na(anno$anno[rownames(beta), "chr"]) & anno$anno[rownames(beta), "chr"] %in% c("chrY", "Y"))
cat(sprintf("chrX/chrY probes present in beta_raw.rds: %d / %d (both should be >0 for the sex-check to be reproducible)\n", x_probes, y_probes))

pca <- methyl_sample_outliers_pca(beta, n_features = 5000, sd_threshold = 3)
stopifnot(isTRUE(pca$ok))

regen <- data.frame(gsm = rownames(pca$scores), PC1_regen = pca$scores[, 1],
                     PC2_regen = pca$scores[, 2], outlier_regen = pca$outlier)

precomp <- read.csv("data/preloaded/methylomics/tables/script01_dataload_QC/tables/sample_metadata_qc.csv")
cmp <- merge(precomp, regen, by = "gsm")
stopifnot(nrow(cmp) == nrow(precomp))

r_pc1 <- cor(cmp$PC1, cmp$PC1_regen)
r_pc2 <- cor(cmp$PC2, cmp$PC2_regen)
## PC1/PC2 sign is arbitrary in PCA (prcomp can flip sign run to run) - compare
## the absolute correlation, and separately whichever orientation the
## precomputed table happened to use.
outlier_agree <- mean(cmp$outlier == cmp$outlier_regen)

cat(sprintf("\n=== QC regeneration vs precomputed sample_metadata_qc.csv (n=%d samples) ===\n", nrow(cmp)))
cat(sprintf("PC1 correlation (sign-sensitive): r = %.4f  (|r| = %.4f)\n", r_pc1, abs(r_pc1)))
cat(sprintf("PC2 correlation (sign-sensitive): r = %.4f  (|r| = %.4f)\n", r_pc2, abs(r_pc2)))
cat(sprintf("outlier flag agreement:            %.1f%% (%d/%d)\n", 100 * outlier_agree, sum(cmp$outlier == cmp$outlier_regen), nrow(cmp)))
cat("sex-check (chrY_mean/sex_mismatch):  NOT REPRODUCIBLE - beta_raw.rds has zero chrX/chrY probes (excluded upstream of what's bundled)\n")

write.csv(cmp, "reproduce/methylomics/output/01_qc_comparison.csv", row.names = FALSE)
saveRDS(list(r_pc1 = r_pc1, r_pc2 = r_pc2, outlier_agree = outlier_agree, n = nrow(cmp),
             sexcheck_reproducible = FALSE, sexcheck_reason = "beta_raw.rds has 0 chrX/chrY probes"),
        "reproduce/methylomics/output/01_qc_summary.rds")
cat("\nWrote reproduce/methylomics/output/01_qc_comparison.csv\n")
