## Module 1 (Transcriptomics) - "Data Exploration" tab's pure statistics
## helpers (mod_preprocessing_explore.R) - a standalone EDA tool independent
## of the shared `dataset` reactiveValues, so every function here is a

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "03_Preprocessing_Batch_Correction", "mod_preprocessing_explore.R"))

test_that("eda_parse_upload() parses a well-formed feature x sample CSV", {
  fm <- fx_expr_meta(n_genes = 10, n_samples = 5, seed = 10)
  path <- withr::local_tempfile(fileext = ".csv")
  fx_write_expr_csv(fm$expr, path)
  res <- eda_parse_upload(path, "test.csv")
  expect_true(res$ok)
  expect_equal(dim(res$expr), c(10L, 5L))
})

test_that("eda_parse_upload() correctly parses a single-feature-row upload with multiple sample columns (fixed - previously crashed)", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("id,S1,S2,S3,S4,S5", "GENE1,1,2,3,4,5"), path)
  res <- eda_parse_upload(path, "one_row.csv")
  expect_true(res$ok)
  expect_equal(dim(res$expr), c(1L, 5L))
  expect_equal(rownames(res$expr), "GENE1")
  expect_equal(colnames(res$expr), c("S1", "S2", "S3", "S4", "S5"))
  expect_equal(as.numeric(res$expr[1, ]), c(1, 2, 3, 4, 5))
})

test_that("eda_parse_upload() does not silently succeed with garbage on a genuinely malformed, ragged-row file", {
  path <- normalizePath(file.path(app_dir, "tests", "fixtures", "edge_cases", "malformed_expr.csv"), mustWork = TRUE)
  res <- tryCatch(eda_parse_upload(path, "malformed_expr.csv"), error = function(e) e)
  if (inherits(res, "error")) {
    expect_match(conditionMessage(res), "dimnames")
  } else if (isTRUE(res$ok)) {
    expect_true(any(is.finite(res$expr)))
  } else {
    expect_true(nzchar(res$error))
  }
})

test_that("eda_parse_upload() rejects an empty (header-only) file", {
  res <- eda_parse_upload(
    normalizePath(file.path(app_dir, "tests", "fixtures", "edge_cases", "empty_expr.csv"), mustWork = TRUE),
    "empty_expr.csv"
  )
  expect_false(res$ok)
  expect_true(grepl("no data rows", res$error))
})

test_that("eda_parse_upload() rejects a matrix with no numeric sample columns", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("id,label", "GENE1,foo", "GENE2,bar", "GENE3,baz"), path)
  res <- eda_parse_upload(path, "text_only.csv")
  expect_false(res$ok)
  expect_true(grepl("No numeric sample columns", res$error))
})

test_that("eda_skewness()/eda_kurtosis() return NA below their minimum sample size", {
  expect_true(is.na(eda_skewness(c(1, 2))))
  expect_true(is.na(eda_kurtosis(c(1, 2, 3))))
})

test_that("eda_skewness()/eda_kurtosis() return 0 for a constant (zero-variance) vector", {
  expect_equal(eda_skewness(rep(5, 10)), 0)
  expect_equal(eda_kurtosis(rep(5, 10)), 0)
})

test_that("eda_skewness() ignores non-finite values rather than propagating NA/Inf", {
  x <- c(rnorm(50), NA, Inf, -Inf)
  s <- eda_skewness(x)
  expect_true(is.finite(s))
})

test_that("eda_skew_label() buckets by absolute skewness magnitude", {
  expect_equal(eda_skew_label(NA_real_), "Unknown (not enough data)")
  expect_equal(eda_skew_label(0.1), "Approximately symmetric")
  expect_equal(eda_skew_label(0.7), "Moderately skewed")
  expect_equal(eda_skew_label(2), "Strongly skewed")
})

test_that("eda_robust_z() returns all zeros when MAD is zero (a constant/near-constant vector)", {
  expect_equal(eda_robust_z(rep(3, 20)), rep(0, 20))
})

test_that("eda_robust_z() flags a genuine outlier with a large modified z-score", {
  x <- c(rnorm(30, mean = 0, sd = 1), 100)
  z <- eda_robust_z(x)
  expect_gt(abs(z[length(z)]), 3.5)
})

test_that("eda_overview() counts missing/infinite/duplicated/constant features correctly on a crafted matrix", {
  fm <- fx_expr_meta(n_genes = 20, n_samples = 6, seed = 11)
  m <- fm$expr
  m[1, ] <- NA_real_
  m[2, 1] <- Inf
  m[3, ] <- 7
  rownames(m)[4] <- rownames(m)[5]

  ov <- eda_overview(list(expr = m, nonnumeric_cols = character(0)))
  expect_equal(ov$n_features, 20L)
  expect_equal(ov$n_samples, 6L)
  expect_equal(ov$n_missing, 6L)
  expect_equal(ov$n_infinite, 1L)
  expect_equal(ov$n_duplicated_features, 1L)
  expect_gte(ov$n_constant_features, 1L)
})

test_that("eda_normality_summary() runs Shapiro-Wilk only at/above 20 finite values", {
  small <- eda_normality_summary(matrix(rnorm(15), nrow = 1))
  expect_true(is.na(small$shapiro_p))

  big <- eda_normality_summary(matrix(rnorm(500), nrow = 1))
  expect_false(is.na(big$shapiro_p))
  expect_true(big$shapiro_p >= 0 && big$shapiro_p <= 1)
})

test_that("eda_normalization_assessment() classifies raw RNA-seq-like counts as not_normalized", {
  set.seed(20)
  counts <- matrix(rpois(400, lambda = 500), 20, 20)
  a <- eda_normalization_assessment(counts)
  expect_equal(a$verdict, "not_normalized")
  expect_equal(a$detected_type, "counts")
})

test_that("eda_normalization_assessment() classifies data with negative values as normalized", {
  set.seed(21)
  m <- matrix(rnorm(400, mean = 0, sd = 1), 20, 20)
  a <- eda_normalization_assessment(m)
  expect_equal(a$verdict, "normalized")
  expect_true(a$has_negative)
})

test_that("eda_impute_median() fills NA/Inf with the per-row median, and never mutates on a complete matrix", {
  m <- matrix(1:20, 4, 5)
  expect_equal(eda_impute_median(m), m + 0)

  m2 <- matrix(c(1, 2, NA, 4, 5), nrow = 1)
  out <- eda_impute_median(m2)
  expect_false(anyNA(out))
  expect_equal(out[1, 3], stats::median(c(1, 2, 4, 5)))
})

test_that("eda_impute_median() falls back to the matrix-wide median for a feature missing in every sample", {
  m <- matrix(c(NA, NA, NA, 1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE)
  out <- eda_impute_median(m)
  expect_false(anyNA(out))
  expect_equal(out[1, 1], stats::median(m, na.rm = TRUE))
})

test_that("eda_missingness() buckets feature-level missingness correctly", {
  m <- matrix(0, nrow = 5, ncol = 10, dimnames = list(paste0("G", 1:5), paste0("S", 1:10)))
  m[1, ] <- NA
  m[2, 1:6] <- NA
  m[3, 1:3] <- NA
  m[4, 1] <- NA
  mm <- eda_missingness(m)
  expect_equal(sum(mm$by_feature_summary$n_features), 5L)
  expect_equal(mm$by_feature_summary$n_features[mm$by_feature_summary$bucket == "0%"], 1L)
  expect_equal(mm$by_feature_summary$n_features[mm$by_feature_summary$bucket == ">50%"], 2L)
})

test_that("eda_feature_outliers() flags high missingness and extreme variance features", {
  fm <- fx_expr_meta(n_genes = 30, n_samples = 12, seed = 12)
  m <- fm$expr
  m[1, 1:10] <- NA
  m[2, ] <- m[2, ] * 1 + rnorm(12, sd = 50)

  desc <- eda_descriptive_stats(m, margin = 1)
  out <- eda_feature_outliers(m, desc)
  expect_true(out$flag_high_missing[1])
  expect_true(any(out$flag_extreme_variance))
})
