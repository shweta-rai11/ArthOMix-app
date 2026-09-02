## Module 4 (Cross-omics) - Expression and Methylation Integration
## sub-module, via testServer(): fully synchronous (no ExtendedTask), so
## this drives a REAL full "Run Integration" click on standardized
## Transcriptomics/Methylomics data mirroring what the Dataset tab
## publishes into the shared cross_dataset store - real gene-ID
## harmonization, real Stouffer's-Z methylation aggregation, real
## classification, and (when sample columns are present) real per-gene
## correlation - verifying the published cross_results$integration summary
## matches genuine computed counts.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "01_Data", "crossomics_integration_upload.R"))
source_from_app_root(file.path("R", "crossomics", "functions", "integration", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "crossomics", "02_Expression_Methylation_Integration", "mod_cross_integration.R"))

cx_dataset_fixture <- function() {
  ## 4 genes: A (Hyper+Down candidate), B (Hypo+Up candidate), C (expression-
  ## only significant), D (nothing significant) - a small, hand-worked
  ## fixture whose real classification can be checked by hand.
  expr_df <- data.frame(
    gene = c("A", "B", "C", "D"), log2fc = c(-2, 2, 3, 0.1),
    pvalue = c(0.001, 0.001, 0.001, 0.5), fdr = c(0.001, 0.001, 0.001, 0.5),
    stringsAsFactors = FALSE
  )
  meth_df <- data.frame(
    gene = c("A", "B", "C", "D"), cpg = c("cg1", "cg2", "cg3", "cg4"),
    dbeta = c(0.3, -0.3, 0.01, 0.01), pvalue = c(0.001, 0.001, 0.9, 0.9), fdr = c(0.001, 0.001, 0.9, 0.9),
    chr = NA, pos = NA, region_raw = NA, region = NA, region_fine = NA, island_context = NA,
    stringsAsFactors = FALSE
  )
  list(expr_df = expr_df, meth_df = meth_df)
}

test_that("Run Integration (real, synchronous) correctly classifies a hand-worked gene set and publishes matching summary counts to cross_results", {
  fx <- cx_dataset_fixture()
  cross_dataset <- shiny::reactiveValues(user_expr_df = fx$expr_df, user_expr_source = "Example data (FEMALE, sex-stratified DEG)",
                                           user_meth_df = fx$meth_df, user_meth_source = "Example data (FEMALE, sex-stratified DMP, SVA/bacon-adjusted)")
  cross_results <- shiny::reactiveValues()
  shiny::testServer(mod_cross_integration_server, args = list(id = "ci", cross_dataset = cross_dataset, cross_results = cross_results), {
    session$setInputs(expr_thresh = 1, expr_fdr_thresh = 0.05, meth_thresh = 0.1, meth_fdr_thresh = 0.05,
                        agg_method = "mean", cor_method = "pearson", padj_method = "BH")
    session$setInputs(run_integration = 0)
    session$setInputs(run_integration = 1)

    df <- integ$df
    expect_false(is.null(df))
    expect_equal(as.character(df$category[df$gene == "A"]), "Hyper + Down")
    expect_equal(as.character(df$category[df$gene == "B"]), "Hypo + Up")
    expect_equal(as.character(df$category[df$gene == "C"]), "Not significant")  ## expr-only, not both layers
    expect_equal(as.character(df$category[df$gene == "D"]), "Not significant")

    expect_false(is.null(cross_results$integration))
    expect_equal(cross_results$integration$summary$n_genes, 4L)
    expect_equal(cross_results$integration$summary$n_integrated, 2L)  ## A and B
    expect_equal(integ$params$sex_stratum, "FEMALE")
  })
})

test_that("Run Integration refuses cleanly (no crash, integ$df stays NULL) when Methylomics data hasn't been loaded on the Dataset tab", {
  cross_dataset <- shiny::reactiveValues(user_expr_df = data.frame(gene = "A", log2fc = 1, pvalue = 0.01, fdr = 0.01), user_meth_df = NULL)
  cross_results <- shiny::reactiveValues()
  shiny::testServer(mod_cross_integration_server, args = list(id = "ci", cross_dataset = cross_dataset, cross_results = cross_results), {
    session$setInputs(expr_thresh = 1, expr_fdr_thresh = 0.05, meth_thresh = 0.1, meth_fdr_thresh = 0.05,
                        agg_method = "mean", cor_method = "pearson", padj_method = "BH")
    session$setInputs(run_integration = 0)
    session$setInputs(run_integration = 1)
    expect_null(integ$df)
    expect_null(cross_results$integration)
  })
})

test_that("Run Integration computes a REAL per-gene sample-level correlation when both sides expose matched sample columns", {
  set.seed(2100)
  ids <- paste0("S", 1:10)
  expr_wide <- data.frame(gene = c("A", "B"), matrix(rnorm(20), 2, 10, dimnames = list(NULL, ids)), stringsAsFactors = FALSE)
  meth_wide <- expr_wide
  meth_wide[, ids] <- -expr_wide[, ids] + matrix(rnorm(20, sd = 0.05), 2, 10)  ## strong real negative correlation
  colnames(meth_wide)[1] <- "gene"

  expr_df <- data.frame(gene = c("A", "B"), log2fc = c(-2, 2), pvalue = c(0.001, 0.001), fdr = c(0.001, 0.001), stringsAsFactors = FALSE)
  meth_df <- data.frame(gene = c("A", "B"), cpg = c("cg1", "cg2"), dbeta = c(0.3, -0.3), pvalue = c(0.001, 0.001), fdr = c(0.001, 0.001),
                          chr = NA, pos = NA, region_raw = NA, region = NA, region_fine = NA, island_context = NA, stringsAsFactors = FALSE)

  cross_dataset <- shiny::reactiveValues(
    user_expr_df = expr_df, user_expr_source = "Uploaded: x.csv", user_expr_wide = expr_wide,
    user_expr_mapping = c(gene = "gene"), user_expr_sample_cols = ids,
    user_meth_df = meth_df, user_meth_source = "Uploaded: y.csv", user_meth_wide = meth_wide,
    user_meth_mapping = c(gene = "gene"), user_meth_sample_cols = ids
  )
  cross_results <- shiny::reactiveValues()
  shiny::testServer(mod_cross_integration_server, args = list(id = "ci", cross_dataset = cross_dataset, cross_results = cross_results), {
    session$setInputs(expr_thresh = 1, expr_fdr_thresh = 0.05, meth_thresh = 0.1, meth_fdr_thresh = 0.05,
                        agg_method = "mean", cor_method = "pearson", padj_method = "BH")
    session$setInputs(run_integration = 0)
    session$setInputs(run_integration = 1)

    expect_true(integ$pairing$paired)
    a_row <- integ$df[integ$df$gene == "A", ]
    expect_true(a_row$correlation_r < -0.9)
    expect_equal(as.character(a_row$evidence_level), "Strong candidate")  ## both sig, inverse, real significant negative correlation
  })
})
