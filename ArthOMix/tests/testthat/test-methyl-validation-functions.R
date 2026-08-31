## Module 2 (Methylomics) - External Validation's pure functions: feature-
## set alignment, sample-independence checking, the compatibility audit
## table, bootstrap CI, and small formatting/lookup helpers.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "mod_methyl_validation.R"))

## ---- vld_feature_alignment() -------------------------------------------------

test_that("vld_feature_alignment() correctly reports matched/missing/extra CpGs and overlap percentage", {
  out <- vld_feature_alignment(required_ids = c("cg1", "cg2", "cg3", "cg4"), available_ids = c("cg1", "cg2", "cg5"))
  expect_equal(out$matched, 2L)
  expect_setequal(out$missing, c("cg3", "cg4"))
  expect_setequal(out$extra, "cg5")
  expect_equal(out$overlap_pct, 50)
  expect_false(out$ok)  ## missing required CpGs
})

test_that("vld_feature_alignment() is ok=TRUE only when every required CpG is present", {
  out <- vld_feature_alignment(c("cg1", "cg2"), c("cg1", "cg2", "cg3"))
  expect_true(out$ok)
  expect_equal(out$overlap_pct, 100)
})

## ---- vld_sample_overlap() -----------------------------------------------------

test_that("vld_sample_overlap() detects shared sample IDs between training and validation cohorts", {
  out <- vld_sample_overlap(c("S1", "S2", "S3"), c("S3", "S4"))
  expect_equal(out$status, "overlap")
  expect_equal(out$n_overlap, 1L)
})

test_that("vld_sample_overlap() reports 'independent' when there is no shared ID, and 'unknown' when IDs are unavailable", {
  ind <- vld_sample_overlap(c("S1", "S2"), c("S3", "S4"))
  expect_equal(ind$status, "independent")
  unk <- vld_sample_overlap(character(0), c("S3", "S4"))
  expect_equal(unk$status, "unknown")
})

## ---- vld_compat_table() -------------------------------------------------------

test_that("vld_compat_table() marks the feature/CpG-set row 'Fail' when alignment is incomplete, 'Pass' otherwise", {
  model <- list(feature_ids = c("cg1", "cg2"), threshold = 0.42)
  cohort <- list(scale_declared = "beta", n_missing = 0, normalization_note = "n/a",
                   normalization_status = "Pass", platform = "EPIC", platform_status = "Pass")
  align_ok <- vld_feature_alignment(c("cg1", "cg2"), c("cg1", "cg2"))
  align_bad <- vld_feature_alignment(c("cg1", "cg2"), c("cg1"))

  tbl_ok <- vld_compat_table(model, cohort, align_ok)
  tbl_bad <- vld_compat_table(model, cohort, align_bad)
  expect_equal(tbl_ok$Status[tbl_ok$Component == "Feature / CpG set"], "Pass")
  expect_equal(tbl_bad$Status[tbl_bad$Component == "Feature / CpG set"], "Fail")
})

test_that("vld_compat_table() flags missing values at required CpGs as a Warning row", {
  model <- list(feature_ids = c("cg1"), threshold = 0.5)
  cohort <- list(scale_declared = "beta", n_missing = 3, normalization_note = "n/a",
                   normalization_status = "Pass", platform = "EPIC", platform_status = "Pass")
  align <- vld_feature_alignment("cg1", "cg1")
  tbl <- vld_compat_table(model, cohort, align)
  expect_equal(tbl$Status[tbl$Component == "Missing values"], "Warning")
})

## ---- vld_bootstrap_ci() -------------------------------------------------------

test_that("vld_bootstrap_ci() produces a valid, deterministic (fixed-seed) CI around sensitivity", {
  set.seed(300)
  y <- factor(rep(c("Class0", "Class1"), each = 30), levels = c("Class0", "Class1"))
  prob <- c(runif(30, 0.1, 0.4), runif(30, 0.6, 0.9))
  sens_fn <- function(yy, pp, thr) {
    pred_pos <- pp >= thr; obs_pos <- yy == "Class1"
    sum(pred_pos & obs_pos) / sum(obs_pos)
  }
  ci1 <- vld_bootstrap_ci(y, prob, threshold = 0.5, stat_fn = sens_fn, n_boot = 200, seed = 42)
  ci2 <- vld_bootstrap_ci(y, prob, threshold = 0.5, stat_fn = sens_fn, n_boot = 200, seed = 42)
  expect_identical(ci1, ci2)  ## fixed seed -> reproducible
  expect_true(ci1[1] <= ci1[2])
  expect_true(all(ci1 >= 0 & ci1 <= 1))
})

## ---- Small formatting/lookup helpers -------------------------------------------

test_that("vld_ci_label() formats a valid CI and reports 'NA' for a missing/NULL one", {
  expect_equal(vld_ci_label(c(0.1234, 0.5678)), "0.123-0.568")
  expect_equal(vld_ci_label(c(NA_real_, NA_real_)), "NA")
  expect_equal(vld_ci_label(NULL), "NA")
})

test_that("vld_which_latest() picks the most recent timestamp via character comparison, not numeric coercion", {
  ts <- c("2026-08-01 10:00:00", "2026-08-30 09:00:00", "2026-08-15 12:00:00")
  expect_equal(vld_which_latest(ts), 2L)
})
