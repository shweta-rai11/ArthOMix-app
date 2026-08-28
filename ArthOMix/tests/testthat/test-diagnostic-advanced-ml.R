## Regression guards for two bugs found in the transcriptomics audit
## (2026-08-26) in mod_diagnostic.R's Advanced ML "Compare Models" engine:
##
## 1. Reproducibility: the "Random seed" UI control (cfg$seed) used to only
##    seed the initial train/test split - outer-CV fold assignment, filter
##    fitting, caret tuning, and class-imbalance resampling were all
##    unseeded, so re-running the same comparison did not reproduce the same
##    numbers despite the seed control implying it would.
## 2. Leakage: unsupervised feature filters (variance/missingness/
##    correlation) used to be fit once on the whole outer-CV Train pool in
##    diag_adv_compare_models(), before folds were drawn - leaking each
##    fold's held-out rows into the feature set used to score that same
##    fold. Supervised filters were already correctly refit per fold; this
##    made the unsupervised branch match.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "mod_diagnostic.R"))

diag_adv_test_data <- function(seed = 99, n = 60, p = 150) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p, dimnames = list(paste0("S", 1:n), paste0("gene", 1:p)))
  y <- factor(rep(c("HC", "RA"), each = n / 2), levels = c("HC", "RA"))
  X[y == "RA", 1:5] <- X[y == "RA", 1:5] + 1.5
  list(expr_sub = t(X), y = y)
}

diag_adv_test_cfg <- function(seed) {
  list(
    test_frac = 0.3, seed = seed,
    filter = list(method = "variance", cap_n = 30),
    tune = list(mode = "manual", inner_method = "cv", inner_k = 3, manual_grid = NULL),
    outer = list(k = 3),
    imbalance = list(method = "none"),
    class_weight = list(mode = "equal", ratio = 1),
    threshold = 0.5
  )
}

test_that("Compare Models is fully reproducible given the same seed, config, and data", {
  d <- diag_adv_test_data()
  cfg <- diag_adv_test_cfg(4242)

  r1 <- suppressWarnings(diag_adv_compare_models(d$expr_sub, d$y, c("logistic"), cfg))
  r2 <- suppressWarnings(diag_adv_compare_models(d$expr_sub, d$y, c("logistic"), cfg))

  expect_equal(as.numeric(r1$results$logistic$pred_test), as.numeric(r2$results$logistic$pred_test))
  expect_identical(sort(r1$results$logistic$features_final), sort(r2$results$logistic$features_final))

  cv_auc1 <- vapply(r1$results$logistic$outer, function(o) o$metrics$ROC_AUC, numeric(1))
  cv_auc2 <- vapply(r2$results$logistic$outer, function(o) o$metrics$ROC_AUC, numeric(1))
  expect_equal(cv_auc1, cv_auc2)
})

test_that("a different seed changes the Compare Models result", {
  d <- diag_adv_test_data()
  r1 <- suppressWarnings(diag_adv_compare_models(d$expr_sub, d$y, c("logistic"), diag_adv_test_cfg(4242)))
  r2 <- suppressWarnings(diag_adv_compare_models(d$expr_sub, d$y, c("logistic"), diag_adv_test_cfg(777)))

  cv_auc1 <- vapply(r1$results$logistic$outer, function(o) o$metrics$ROC_AUC, numeric(1))
  cv_auc2 <- vapply(r2$results$logistic$outer, function(o) o$metrics$ROC_AUC, numeric(1))
  expect_false(isTRUE(all.equal(cv_auc1, cv_auc2)))
})

test_that("unsupervised filters are refit per outer fold, not fixed once on the whole Train pool", {
  d <- diag_adv_test_data()
  r <- suppressWarnings(diag_adv_compare_models(d$expr_sub, d$y, c("logistic"), diag_adv_test_cfg(4242)))
  feats_by_fold <- lapply(r$results$logistic$outer, function(o) sort(o$features))

  ## Under the pre-fix behavior, every fold's `features` was the same
  ## globally-prefiltered set (identical across folds, and identical to the
  ## final whole-Train fit's feature set). A real per-fold refit on
  ## genuinely different Train rows should disagree at least sometimes.
  expect_false(identical(feats_by_fold[[1]], feats_by_fold[[2]]))
  expect_false(identical(sort(r$results$logistic$features_final), feats_by_fold[[1]]))
})
