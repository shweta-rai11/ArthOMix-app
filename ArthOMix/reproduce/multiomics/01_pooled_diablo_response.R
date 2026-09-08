## Independent sanity check of the app's headline Multi-Omics finding: does a
## live DIABLO nested-CV fit, run on freshly-fetched real GEO data from the
## Tao et al. 2021 cohort (GSE138746 RNA-seq + GSE138653 methylation, ground-
## truth-confirmed via NCBI E-utilities - see session notes), recover
## meaningful discrimination for anti-TNF response, closer to the paper's own
## published 0.79-0.89 accuracy range, or does it also land near-chance like
## the app's bundled sex-stratified numbers (female 0.411 / male 0.565,
## Table8_benchmark_vs_published.csv)?
##
## Uses the app's own leakage-safe engine (mss_nested_cv - feature selection
## refit inside every fold, out-of-fold AUROC) in its built-in POOLED mode (no
## sex split - sex labels for this cohort could not be reliably matched
## between the two independently-submitted GEO series; see session notes).
## Response binarized EULAR-standard: Good+Moderate = Responder, No = NonResponder.
##
## Data: data/examples/multiomics_upload/gse138746_pbmc_rnaseq_counts.csv (RNA-seq,
## 18333 genes x 79 patients) + gse138653_methylation_beta_top2000.csv (2000 CpGs x
## 79 patients) + gse138747_sample_metadata.csv (drug/response) - pre-existing in
## this repo, independently legitimacy-checked this session (realistic count/beta
## ranges, drug counts 37 ADA/42 ETN match the paper's 38/42 almost exactly).
##
## Run from the ArthOMix/ app directory:
##   Rscript reproduce/multiomics/01_pooled_diablo_response.R

suppressMessages(suppressWarnings(
  shiny::loadSupport(".", renv = globalenv(), globalrenv = globalenv())
))

rna <- as.matrix(read.csv("data/examples/multiomics_upload/gse138746_pbmc_rnaseq_counts.csv", check.names = FALSE, row.names = 1))
meth <- as.matrix(read.csv("data/examples/multiomics_upload/gse138653_methylation_beta_top2000.csv", check.names = FALSE, row.names = 1))
meta <- read.csv("data/examples/multiomics_upload/gse138747_sample_metadata.csv")
rownames(meta) <- meta$sample

common <- Reduce(intersect, list(colnames(rna), colnames(meth), meta$sample))
cat(sprintf("Patients with expression + methylation + metadata: %d\n", length(common)))

## log2(CPM+1) - RNA-seq raw counts aren't directly usable by limma's
## Gaussian moderated-t feature ranking or DIABLO's sPLS-DA.
libsize <- colSums(rna[, common])
expr <- t(log2(sweep(rna[, common], 2, libsize, "/") * 1e6 + 1))
meth_t <- t(meth[, common])

response_raw <- meta[common, "response"]
outcome <- factor(ifelse(response_raw == "NoResponse", "NonResponder", "Responder"),
                   levels = c("NonResponder", "Responder"))
names(outcome) <- common
cat(sprintf("Outcome: %s\n", paste(names(table(outcome)), table(outcome), sep = "=", collapse = ", ")))

rownames(expr) <- common; rownames(meth_t) <- common

cat("\nRunning mss_nested_cv (DIABLO engine, pooled, 5x5 CV, leakage-safe in-fold feature selection)...\n")
res <- mss_nested_cv(expr, meth_t, outcome, covariate = NULL, engine = "diablo")

if (!isTRUE(res$ok)) {
  cat("FAILED:", res$error, "\n")
} else {
  perf <- res$performance
  cat(sprintf("\n=== Pooled DIABLO nested-CV result (n=%d) ===\n", perf$n))
  cat(sprintf("AUROC: %.3f [%.3f, %.3f]  (%s)\n", perf$auroc, perf$ci_lo, perf$ci_hi, perf$auroc_call))
  cat(sprintf("\nFor comparison (from the app's own bundled Table8_benchmark_vs_published.csv):\n"))
  cat("  Published (Tao et al. 2021, per-drug single-omics accuracy): 0.79-0.89\n")
  cat("  App's own sex-stratified integrated-model accuracy: female=0.411, male=0.565\n")
  write.csv(res$scores, "reproduce/multiomics/output/01_pooled_diablo_scores.csv", row.names = FALSE)
  saveRDS(res, "reproduce/multiomics/output/01_pooled_diablo_result.rds")
}
