## Module 1 (Transcriptomics) - WGCNA's sample/gene-count guardrails,
## complementing test-wgcna-qc-gate.R's zero-variance-gene regression guard
## with the surrounding boundary checks: <15 samples up front, <20 genes or

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "05_WGCNA", "mod_wgcna.R"))

wgcna_fixture <- function(n_genes = 60, n_samples = 20, seed = 90) {
  set.seed(seed)
  expr <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1), n_genes, n_samples,
                  dimnames = list(paste0("gene", 1:n_genes), paste0("S", 1:n_samples)))
  meta <- data.frame(sample = colnames(expr), group = rep(c("HC", "RA"), length.out = n_samples),
                      stringsAsFactors = FALSE)
  shiny::reactiveValues(expr = expr, meta = meta, source = "wgcna guardrail test",
                          source_type = "uploaded", is_bundled_reference = FALSE)
}

test_that("fewer than 15 samples in the dataset is rejected up front (qc())", {
  dataset <- wgcna_fixture(n_samples = 10)
  shiny::testServer(mod_wgcna_server, args = list(id = "wg", dataset = dataset, results = shiny::reactiveValues()), {
    session$setInputs(gene_filter_method = "all", exclude_pattern = "", remove_outliers = FALSE)
    err <- tryCatch(qc(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Fewer than 15 samples", conditionMessage(err)))
  })
})

test_that("exactly 15 samples clears the sample-count gate", {
  dataset <- wgcna_fixture(n_samples = 15)
  shiny::testServer(mod_wgcna_server, args = list(id = "wg", dataset = dataset, results = shiny::reactiveValues()), {
    session$setInputs(gene_filter_method = "all", exclude_pattern = "", remove_outliers = FALSE)
    d <- qc()
    expect_equal(nrow(d$meta), 15L)
  })
})

test_that("fewer than 20 genes remaining after a strict variance-percentile filter is rejected", {
  dataset <- wgcna_fixture(n_genes = 60, n_samples = 20)
  shiny::testServer(mod_wgcna_server, args = list(id = "wg", dataset = dataset, results = shiny::reactiveValues()), {
    session$setInputs(gene_filter_method = "var_pct", var_pct = 99.9, exclude_pattern = "")
    err <- tryCatch(gene_selection(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Fewer than .* genes remain after this filter", conditionMessage(err)))
  })
})

test_that("an invalid exclude_pattern regex is rejected with a clear message rather than a raw grepl error", {
  dataset <- wgcna_fixture()
  shiny::testServer(mod_wgcna_server, args = list(id = "wg", dataset = dataset, results = shiny::reactiveValues()), {
    session$setInputs(gene_filter_method = "all", exclude_pattern = "(unclosed[")
    err <- tryCatch(gene_selection(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Invalid regex pattern", conditionMessage(err)))
  })
})

test_that("removing outlier samples down to fewer than 15 remaining is rejected by final_input()", {
  set.seed(91)
  n_genes <- 60
  expr <- matrix(rnorm(n_genes * 16, mean = 8, sd = 1), n_genes, 16,
                  dimnames = list(paste0("gene", 1:n_genes), paste0("S", 1:16)))
  expr[, 12:16] <- expr[, 12:16] + 40
  meta <- data.frame(sample = colnames(expr), group = rep(c("HC", "RA"), length.out = 16), stringsAsFactors = FALSE)
  dataset <- shiny::reactiveValues(expr = expr, meta = meta, source = "outlier test",
                                     source_type = "uploaded", is_bundled_reference = FALSE)

  shiny::testServer(mod_wgcna_server, args = list(id = "wg", dataset = dataset, results = shiny::reactiveValues()), {
    session$setInputs(gene_filter_method = "all", exclude_pattern = "")
    tree <- sample_tree()
    cut_h <- tree$height[length(tree$height)] * 0.5
    session$setInputs(remove_outliers = TRUE, outlier_height = cut_h)
    err <- tryCatch(final_input(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Fewer than 15 samples remain after outlier removal", conditionMessage(err)))
  })
})
