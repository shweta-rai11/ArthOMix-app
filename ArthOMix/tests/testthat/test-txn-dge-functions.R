## Module 1 (Transcriptomics) - Differential Gene Expression tab's top-level
## pure upload-cleanup helper, dge_clean_expr_matrix().

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "functions", "expression_type.R"))
source_from_app_root(file.path("R", "transcriptomics", "04_Differential_Expression", "mod_dge.R"))

test_that("dge_clean_expr_matrix() collapses duplicate feature IDs to their mean expression", {
  m <- matrix(c(1, 3, 5, 7), nrow = 2, byrow = TRUE, dimnames = list(c("G1", "G1"), c("S1", "S2")))
  out <- dge_clean_expr_matrix(m)
  expect_equal(nrow(out$mat), 1L)
  expect_equal(unname(out$mat["G1", ]), c(mean(c(1, 5)), mean(c(3, 7))))
  expect_true(any(grepl("collapsed to their mean", out$notes)))
})

test_that("dge_clean_expr_matrix() drops all-missing rows and all-missing columns", {
  m <- matrix(c(NA, NA, 1, 2, 3, 4), nrow = 2, byrow = TRUE, dimnames = list(c("ALLNA", "G2"), c("S1", "S2", "S3")))
  m["G2", "S1"] <- NA
  out <- dge_clean_expr_matrix(m)
  expect_false("ALLNA" %in% rownames(out$mat))
  expect_true(any(grepl("no data \\(all missing\\) removed", out$notes)))
})

test_that("dge_clean_expr_matrix() drops zero-variance features", {
  m <- matrix(c(5, 5, 5, 1, 2, 3), nrow = 2, byrow = TRUE, dimnames = list(c("CONST", "VAR"), c("S1", "S2", "S3")))
  out <- dge_clean_expr_matrix(m)
  expect_false("CONST" %in% rownames(out$mat))
  expect_true("VAR" %in% rownames(out$mat))
  expect_true(any(grepl("zero-variance", out$notes)))
})

test_that("dge_clean_expr_matrix() is a no-op (besides type coercion) on an already-clean matrix", {
  fm <- fx_expr_meta(n_genes = 10, n_samples = 6, seed = 60)
  out <- dge_clean_expr_matrix(fm$expr)
  expect_equal(dim(out$mat), dim(fm$expr))
  expect_length(out$notes, 0)
})
