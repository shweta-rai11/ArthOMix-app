## Regression guard for the deconvolution scale bug found in the 2026-09-07
## defense audit: mod_deconvolution.R only treated declared_data_type == "raw"
## as already-linear-scale and exponentiated (2^x) everything else, including
## "normalized" (TPM/FPKM/CPM) - which is the UI's own default declared type
## and is already linear. That inflated CIBERSORT's input by ~10^15x.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "functions", "expression_type.R"))

test_that("raw counts are treated as linear-scale (no exponentiation)", {
  expect_true(deconv_is_linear_scale("raw", matrix(1:10, nrow = 2)))
})

test_that("normalized (TPM/FPKM/CPM) is treated as linear-scale (no exponentiation)", {
  expect_true(deconv_is_linear_scale("normalized", matrix(1:10, nrow = 2)))
})

test_that("log-transformed data is NOT treated as linear-scale (gets exponentiated)", {
  expect_false(deconv_is_linear_scale("logtransformed", matrix(1:10, nrow = 2)))
})

test_that("with no declared type, wide-dynamic-range values fall back to linear-scale", {
  set.seed(1)
  raw_like <- matrix(rpois(200, lambda = 500), nrow = 20)
  expect_true(deconv_is_linear_scale(NA_character_, raw_like))
})

test_that("with no declared type, pinned-total values fall back to linear-scale", {
  set.seed(1)
  norm_like <- matrix(runif(200, 0, 50), nrow = 20)
  norm_like <- sweep(norm_like, 2, colSums(norm_like) / 1e6, "/")
  expect_true(deconv_is_linear_scale(NA_character_, norm_like))
})

test_that("with no declared type, log2-scale values (small dynamic range, unpinned totals) are not linear", {
  set.seed(1)
  log_like <- matrix(rnorm(200, mean = 5, sd = 2), nrow = 20)
  expect_false(deconv_is_linear_scale(NA_character_, log_like))
})
