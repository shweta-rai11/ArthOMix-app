## Regression guard for a gap found in the transcriptomics audit (2026-08-26):
## neither CIBERSORT (LM22, via IOBR) nor MCP-counter validate that
## rownames(expr) are actually HUGO gene symbols - both silently intersect

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "functions", "expression_type.R"))
source_from_app_root(file.path("R", "transcriptomics", "15_Immune_Deconvolution", "mod_deconvolution.R"))

test_that("real gene symbols score a high overlap percentage", {
  d0 <- load_default_dataset()
  pct <- deconv_gene_id_overlap_pct(rownames(d0$expr)[1:400])
  expect_gt(pct, 90)
})

test_that("Ensembl-ID-shaped identifiers score a near-zero overlap percentage", {
  set.seed(5)
  fake_ensembl <- sprintf("ENSG%011d", sample.int(9e8, 300))
  pct <- deconv_gene_id_overlap_pct(fake_ensembl)
  expect_lt(pct, 5)
})

test_that("deconv_gene_id_overlap_pct() sits right at the 10% block threshold for a mixed real/fake gene ID matrix", {
  set.seed(6)
  real_genes <- rownames(load_default_dataset()$expr)[1:10]
  fake_genes <- sprintf("ENSG%011d", sample.int(9e8, 90))
  pct <- deconv_gene_id_overlap_pct(c(real_genes, fake_genes))
  expect_equal(pct, 10)
})

test_that("result() blocks a mismatched-gene-ID matrix before running CIBERSORT/MCP-counter", {
  set.seed(5)
  n_genes <- 300; n_samples <- 20
  fake_ensembl <- sprintf("ENSG%011d", sample.int(9e8, n_genes))
  expr <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1), n_genes, n_samples,
                  dimnames = list(fake_ensembl, paste0("S", 1:n_samples)))
  meta <- data.frame(sample = colnames(expr), group = rep(c("HC", "RA"), each = n_samples / 2),
                      stringsAsFactors = FALSE)
  dataset <- shiny::reactiveValues(expr = expr, meta = meta, source = "Uploaded dataset: ensembl_test.csv")

  shiny::testServer(mod_deconvolution_server, args = list(id = "dc", dataset = dataset, results = shiny::reactiveValues()), {
    session$setInputs(group_col = "group", cib_perm = 10)
    session$setInputs(run_btn = 1)
    expect_error(result(), class = "validation")
  })
})
