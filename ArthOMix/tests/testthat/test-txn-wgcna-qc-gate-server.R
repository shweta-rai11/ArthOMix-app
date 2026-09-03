## Regression guard for a gap found in the transcriptomics audit (2026-08-26):
## mod_wgcna.R never called WGCNA::goodSamplesGenes() (or an equivalent) before
## network construction, unlike the sibling methylomics WGCNA module. The

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "05_WGCNA", "mod_wgcna.R"))

test_that("a zero-variance gene is flagged and dropped before network construction", {
  set.seed(1)
  n_genes <- 60; n_samples <- 20
  expr <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1), n_genes, n_samples,
                  dimnames = list(paste0("gene", 1:n_genes), paste0("S", 1:n_samples)))
  expr["gene1", ] <- 5
  meta <- data.frame(sample = colnames(expr), group = rep(c("HC", "RA"), each = n_samples / 2),
                      stringsAsFactors = FALSE)
  dataset <- shiny::reactiveValues(expr = expr, meta = meta, source = "Example dataset: synthetic QC test")

  shiny::testServer(mod_wgcna_server, args = list(id = "wg", dataset = dataset, results = shiny::reactiveValues()), {
    session$setInputs(gene_filter_method = "all", exclude_pattern = "", remove_outliers = FALSE)
    fi <- final_input()
    expect_true("gene1" %in% fi$removed_genes)
    expect_false("gene1" %in% colnames(fi$texpr))
    expect_equal(ncol(fi$texpr), n_genes - 1L)
  })
})
