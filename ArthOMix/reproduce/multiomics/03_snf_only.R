## Runs just the SNF joint-clustering portion of 02_all_arms.R (fast, was
## silently broken there due to a names()-dropping bug and a $p vs $p_value
## field mismatch, both fixed here and back-ported to 02_all_arms.R).
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
response_named <- stats::setNames(m$response_bin, rownames(m))

snf_results <- list()
for (drug in c("Adalimumab", "Etanercept")) {
  ids <- rownames(m)[m$drug == drug]
  layers <- list(expression = expr[ids, , drop = FALSE], methylation = meth_t[ids, , drop = FALSE])
  snf <- sfc_snf_run(layers, params = list())
  if (!isTRUE(snf$ok)) { cat(sprintf("[SNF %s] FAILED: %s\n", drug, snf$error)); next }
  assoc <- sfc_test_categorical(snf$clusters, response_named)
  if (!isTRUE(assoc$ok)) { cat(sprintf("[SNF %s] clustering ok, association FAILED: %s\n", drug, assoc$error)); next }
  cat(sprintf("[SNF %s] n=%d, %d clusters (sizes: %s), response association p=%.4f (Fisher), Cramer's V=%.3f\n",
              drug, length(snf$clusters), snf$params$n_clusters, paste(table(snf$clusters), collapse = ","),
              assoc$p_value, assoc$effect))
  snf_results[[drug]] <- list(drug = drug, n = length(snf$clusters), n_clusters = snf$params$n_clusters,
                               fisher_p = assoc$p_value, cramers_v = assoc$effect)
}
saveRDS(snf_results, "reproduce/multiomics/output/02_snf_results.rds")
