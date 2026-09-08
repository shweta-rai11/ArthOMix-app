## Proves multi_geo_dual_fetch() succeeds end-to-end when GEO actually has
## real values on both sides, using a real, independently-verified SuperSeries
## (GSE117931 -> sub-series GSE117928 + GSE117929, from
## data/examples/multiomics_upload/geo_multiomics/, already verified in a
## prior session) - fed to the NEW dual-accession path (as if they were two
## unlinked series a user found separately) rather than the existing
## autosplit path, to test the new code specifically.
suppressMessages(suppressWarnings(
  shiny::loadSupport(".", renv = globalenv(), globalrenv = globalenv())
))

cat("=== Improved error message for the RNA-seq-supplementary-file case ===\n")
r0 <- multi_geo_dual_fetch("GSE138746", "GSE138653")
cat(sprintf("ok=%s\n%s\n", r0$ok, r0$error))
stopifnot(!isTRUE(r0$ok), grepl("supplementary file", r0$error))

cat("\n=== Real success case: GSE117928 + GSE117929 (sub-series of verified SuperSeries GSE117931) ===\n")
r1 <- multi_geo_dual_fetch("GSE117928", "GSE117929")
if (!isTRUE(r1$ok)) {
  cat("FAILED:", r1$error, "\n")
} else {
  cat(sprintf("ok=TRUE\nexpression: %s, %d samples x %d features\nmethylation: %s, %d samples x %d features\n",
              r1$expression$accession, nrow(r1$expression$mat), ncol(r1$expression$mat),
              r1$methylation$accession, nrow(r1$methylation$mat), ncol(r1$methylation$mat)))
  n_id <- sum(!is.na(r1$expression$meta$geo_title_patient_num))
  cat(sprintf("geo_title_patient_num populated for %d / %d expression samples\n", n_id, nrow(r1$expression$mat)))

  ## Order independence check
  r1b <- multi_geo_dual_fetch("GSE117929", "GSE117928")
  stopifnot(isTRUE(r1b$ok), identical(r1$expression$accession, r1b$expression$accession))
  cat("Order-independence confirmed (GSE117929, GSE117928 gives the same expression/methylation assignment)\n")

  cat("\n=== ALL CHECKS PASSED ===\n")
}
