## Extends 01_pooled_diablo_response.R (which tested ONE arm: pooled DIABLO
## response prediction) to cover the full set of scenarios the app's bundled
## tables report: drug-stratified x sex-stratified DIABLO (Table29-equivalent),
## sex-stratified response DIABLO (Table34-equivalent), the same two with the
## RF engine (Table37/39-equivalent), and SNF joint clustering per drug
## (Table_SNFjoint-equivalent). All on the same freshly-fetched, independently
## legitimacy-checked GSE138746 (RNA-seq) + GSE138653 (methylation) data as
## 01_pooled_diablo_response.R.
##
## Sex labels (blocked in the first pass) are now resolved: GSE138746 and
## GSE138653 sample titles share an identical trailing patient number
## (e.g. "PBMC_E_n_01" / "DNA_E_n_01") that exactly matches the pre-existing
## PT## metadata's patient_id (100% drug/response agreement on merge, verified
## this session - see reproduce/multiomics/output/patient_sex_lookup.csv).
##
## Run from the ArthOMix/ app directory (after 01_pooled_diablo_response.R):
##   Rscript reproduce/multiomics/02_all_arms.R

suppressMessages(suppressWarnings(
  shiny::loadSupport(".", renv = globalenv(), globalrenv = globalenv())
))

rna <- as.matrix(read.csv("data/examples/multiomics_upload/gse138746_pbmc_rnaseq_counts.csv", check.names = FALSE, row.names = 1))
meth <- as.matrix(read.csv("data/examples/multiomics_upload/gse138653_methylation_beta_top2000.csv", check.names = FALSE, row.names = 1))
meta0 <- read.csv("data/examples/multiomics_upload/gse138747_sample_metadata.csv")
sexlk <- read.csv("reproduce/multiomics/output/patient_sex_lookup.csv")
meta <- merge(meta0, sexlk[, c("sample", "sex")], by = "sample")
rownames(meta) <- meta$sample

common <- Reduce(intersect, list(colnames(rna), colnames(meth), meta$sample))
libsize <- colSums(rna[, common])
expr <- t(log2(sweep(rna[, common], 2, libsize, "/") * 1e6 + 1))
meth_t <- t(meth[, common])
rownames(expr) <- common; rownames(meth_t) <- common

m <- meta[common, ]
m$response_bin <- ifelse(m$response == "NoResponse", "NonResponder", "Responder")
cat(sprintf("Full cohort: n=%d, sex: %s\n", nrow(m), paste(names(table(m$sex)), table(m$sex), sep = "=", collapse = ", ")))

run_arm <- function(label, ids, engine) {
  sub_meta <- m[ids, , drop = FALSE]
  res <- mss_run_stratified(expr[ids, , drop = FALSE], meth_t[ids, , drop = FALSE], sub_meta,
                             outcome_col = "response_bin", sex_mode = "both", engine = engine)
  if (!isTRUE(res$ok)) {
    cat(sprintf("[%s / %s] FAILED: %s\n", label, engine, res$error))
    return(NULL)
  }
  perf <- res$performance
  perf$arm <- label
  for (i in seq_len(nrow(perf))) {
    cat(sprintf("[%s / %s / %s] n=%d AUROC=%.3f [%.3f,%.3f] (%s)\n",
                label, engine, perf$stratum[i], perf$n[i], perf$auroc[i], perf$ci_lo[i], perf$ci_hi[i], perf$auroc_call[i]))
  }
  perf
}

all_perf <- list()

cat("\n=== Response prediction, sex-stratified, pooled across drugs (Table34/39-equivalent) ===\n")
all_perf$response_diablo <- run_arm("response_pooled_drug", rownames(m), "diablo")
all_perf$response_rf <- run_arm("response_pooled_drug", rownames(m), "rf")

for (drug in c("Adalimumab", "Etanercept")) {
  ids <- rownames(m)[m$drug == drug]
  cat(sprintf("\n=== Drug-sex arm: %s, n=%d (Table29/37-equivalent) ===\n", drug, length(ids)))
  all_perf[[paste0(drug, "_diablo")]] <- run_arm(paste0("drugsex_", drug), ids, "diablo")
  all_perf[[paste0(drug, "_rf")]] <- run_arm(paste0("drugsex_", drug), ids, "rf")
}

perf_all <- do.call(rbind, Filter(Negate(is.null), all_perf))
write.csv(perf_all, "reproduce/multiomics/output/02_all_arms_performance.csv", row.names = FALSE)

cat("\n=== SNF joint clustering, per drug (Table_SNFjoint-equivalent) ===\n")
snf_results <- list()
for (drug in c("Adalimumab", "Etanercept")) {
  ids <- rownames(m)[m$drug == drug]
  layers <- list(expression = expr[ids, , drop = FALSE], methylation = meth_t[ids, , drop = FALSE])
  snf <- sfc_snf_run(layers, params = list())
  if (!isTRUE(snf$ok)) { cat(sprintf("[SNF %s] FAILED: %s\n", drug, snf$error)); next }
  response_named <- stats::setNames(m$response_bin, rownames(m))
  assoc <- sfc_test_categorical(snf$clusters, response_named)
  if (!isTRUE(assoc$ok)) { cat(sprintf("[SNF %s] clustering ok, association test FAILED: %s\n", drug, assoc$error)); next }
  cat(sprintf("[SNF %s] n=%d, %d clusters, response association p=%.4f (Fisher), Cramer's V=%.3f\n",
              drug, length(snf$clusters), snf$params$n_clusters, assoc$p_value, assoc$effect))
  snf_results[[drug]] <- list(drug = drug, n = length(snf$clusters), n_clusters = snf$params$n_clusters,
                               fisher_p = assoc$p_value, cramers_v = assoc$effect)
}
saveRDS(snf_results, "reproduce/multiomics/output/02_snf_results.rds")

cat("\n=== FULL SUMMARY: all tested arms ===\n")
print(perf_all[, c("arm", "stratum", "n", "auroc", "ci_lo", "ci_hi", "auroc_call", "engine")])
n_chance <- sum(perf_all$auroc_call == "chance")
n_above <- sum(perf_all$auroc_call == "above_chance")
n_below <- sum(perf_all$auroc_call == "below_chance")
cat(sprintf("\n%d arms tested: %d chance, %d above chance, %d below chance\n", nrow(perf_all), n_chance, n_above, n_below))
