## Module 1 (Transcriptomics) - shared expression-matrix scale/type helpers
## (R/transcriptomics/functions/expression_type.R): looks_like_raw_counts()/
## looks_like_normalized_totals() (promoted out of mod_dge.R, and de-

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "functions", "expression_type.R"))

raw_counts_matrix <- function(seed = 300, n_genes = 40, n_samples = 8) {
  set.seed(seed)
  lambda <- sample(c(50, 500, 5000, 20000), n_genes, replace = TRUE)
  m <- matrix(rpois(n_genes * n_samples, lambda = lambda), n_genes, n_samples,
              dimnames = list(paste0("GENE", 1:n_genes), paste0("S", 1:n_samples)))
  m
}

tpm_matrix <- function(seed = 301, n_genes = 40, n_samples = 8) {
  m <- raw_counts_matrix(seed = seed, n_genes = n_genes, n_samples = n_samples)
  sweep(m, 2, colSums(m), FUN = "/") * 1e6
}

logtransformed_matrix <- function(seed = 302, n_genes = 40, n_samples = 8) {
  set.seed(seed)
  matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1.5), n_genes, n_samples,
         dimnames = list(paste0("GENE", 1:n_genes), paste0("S", 1:n_samples)))
}

test_that("looks_like_raw_counts() flags a raw Poisson-count-like matrix and rejects log-scale/negative data", {
  expect_true(looks_like_raw_counts(raw_counts_matrix()))
  expect_false(looks_like_raw_counts(logtransformed_matrix()))
  neg <- logtransformed_matrix(); neg[1, 1] <- -5
  expect_false(looks_like_raw_counts(neg))
})

test_that("looks_like_normalized_totals() flags per-sample totals tightly pinned near 1e6 (TPM-like) and rejects raw counts", {
  expect_true(looks_like_normalized_totals(tpm_matrix()))
  expect_false(looks_like_normalized_totals(raw_counts_matrix()))
})

test_that("tx_looks_like_results_table() flags DE-results-shaped column names and passes ordinary sample matrices", {
  results_like <- matrix(1, 5, 4, dimnames = list(paste0("G", 1:5), c("logFC", "P.Value", "padj", "baseMean")))
  expect_true(tx_looks_like_results_table(results_like))
  expect_false(tx_looks_like_results_table(raw_counts_matrix()))
})

test_that("tx_validate_expr_upload() accepts clean raw counts declared 'raw'", {
  res <- tx_validate_expr_upload(raw_counts_matrix(), declared_type = "raw")
  expect_true(res$ok)
  expect_null(res$note)
})

test_that("tx_validate_expr_upload() blocks TPM-normalized data declared 'raw'", {
  res <- tx_validate_expr_upload(tpm_matrix(), declared_type = "raw")
  expect_false(res$ok)
  expect_true(grepl("TPM/FPKM/CPM-normalized", res$error))
})

test_that("tx_validate_expr_upload() blocks negative/non-count values declared 'raw'", {
  m <- logtransformed_matrix()
  m[1, ] <- -abs(m[1, ]) - 1
  res <- tx_validate_expr_upload(m, declared_type = "raw")
  expect_false(res$ok)
  expect_true(grepl("negative values", res$error))
})

test_that("tx_validate_expr_upload() blocks a differential-results-table-shaped input regardless of declared type", {
  results_like <- matrix(c(1.2, -0.8, 0.5, 2.1, 0.01, 0.2, 0.03, 0.5, 0.02, 0.3, 0.04, 0.6, 100, 200, 300, 400),
                          4, 4, dimnames = list(paste0("G", 1:4), c("logFC", "P.Value", "padj", "baseMean")))
  res_raw <- tx_validate_expr_upload(results_like, declared_type = "raw")
  res_norm <- tx_validate_expr_upload(results_like, declared_type = "normalized")
  expect_false(res_raw$ok)
  expect_false(res_norm$ok)
  expect_true(grepl("differential-expression results table", res_raw$error))
})

test_that("tx_validate_expr_upload() blocks raw-count-shaped data declared 'normalized' or 'logtransformed'", {
  res_norm <- tx_validate_expr_upload(raw_counts_matrix(), declared_type = "normalized")
  expect_false(res_norm$ok)
  expect_true(grepl("raw, un-normalized sequencing counts", res_norm$error))

  res_log <- tx_validate_expr_upload(raw_counts_matrix(), declared_type = "logtransformed")
  expect_false(res_log$ok)
  expect_true(grepl("raw, un-normalized sequencing counts", res_log$error))
})

test_that("tx_validate_expr_upload() accepts TPM data declared 'normalized' and log-scale data declared 'logtransformed'", {
  res_norm <- tx_validate_expr_upload(tpm_matrix(), declared_type = "normalized")
  expect_true(res_norm$ok)

  res_log <- tx_validate_expr_upload(logtransformed_matrix(), declared_type = "logtransformed")
  expect_true(res_log$ok)
})

test_that("tx_validate_expr_upload() warns but does not block an ambiguous declared-'raw' matrix with a narrow value range", {
  set.seed(303)
  m <- matrix(rpois(40 * 6, lambda = 20), 40, 6, dimnames = list(paste0("G", 1:40), paste0("S", 1:6)))
  res <- tx_validate_expr_upload(m, declared_type = "raw")
  expect_true(res$ok)
  expect_true(grepl("wide dynamic range", res$note))
})

test_that("tx_validate_expr_upload() rejects an all-non-finite matrix", {
  m <- matrix(NA_real_, 5, 4, dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))
  res <- tx_validate_expr_upload(m, declared_type = "raw")
  expect_false(res$ok)
  expect_true(grepl("No finite numeric values", res$error))
})

test_that("tx_validate_expr_upload() is a pass-through (fail-soft) when declared_type is NA/unset", {
  res <- tx_validate_expr_upload(tpm_matrix(), declared_type = NA_character_)
  expect_true(res$ok)
})
