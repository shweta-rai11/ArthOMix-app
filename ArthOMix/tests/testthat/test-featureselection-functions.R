## Module 1 (Transcriptomics) - Feature Selection's pure/near-pure
## computational core: class-weighting helpers and the LASSO+RF+SVM-RFE
## per-sex fitting pipeline (fs_fit_sex). Manual mtry/cost params are used
## throughout to skip fs_fit_sex's own internal CV grid search (caret::train/
## e1071::tune), keeping these tests fast without changing which code path
## is exercised - fs_fit_sex's own "manual" branches are real production
## code, not a test-only shortcut.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "mod_featureselection.R"))

## ---- fs_class_weight_levels() / fs_obs_weights() --------------------------

test_that("fs_class_weight_levels() computes inverse-frequency weights in 'balanced' mode", {
  y <- factor(c(rep("HC", 8), rep("RA", 2)), levels = c("HC", "RA"))
  w <- fs_class_weight_levels(y, "balanced", ratio = NULL)
  expect_equal(unname(w["HC"]), 1)
  expect_equal(unname(w["RA"]), 4)  ## max(n)/n = 8/2
})

test_that("fs_class_weight_levels() uses the manual ratio for the second level, and 1 for both in 'equal' mode", {
  y <- factor(c("HC", "RA"), levels = c("HC", "RA"))
  w <- fs_class_weight_levels(y, "manual", ratio = 3)
  expect_equal(unname(w), c(1, 3))

  w_eq <- fs_class_weight_levels(y, "equal", ratio = NULL)
  expect_equal(unname(w_eq), c(1, 1))
})

test_that("fs_obs_weights() maps per-observation weights from fs_class_weight_levels()", {
  y <- factor(c("HC", "RA", "HC", "RA"), levels = c("HC", "RA"))
  w <- fs_obs_weights(y, "manual", ratio = 5)
  expect_equal(w, c(1, 5, 1, 5))
})

## ---- fs_svm_rfe_rank() / fs_svm_rfe_curve() --------------------------------

fs_separable_fixture <- function(n_per_group = 10, n_genes = 8, seed = 110) {
  set.seed(seed)
  n <- n_per_group * 2
  y <- factor(rep(c("HC", "RA"), each = n_per_group), levels = c("HC", "RA"))
  X <- matrix(rnorm(n * n_genes), n, n_genes, dimnames = list(NULL, paste0("G", 1:n_genes)))
  ## G1/G2 carry a large, clean separation; the rest are pure noise.
  X[y == "RA", 1] <- X[y == "RA", 1] + 8
  X[y == "RA", 2] <- X[y == "RA", 2] + 8
  list(X = X, y = y)
}

test_that("fs_svm_rfe_rank() returns a full permutation of the feature names, most-to-least important", {
  fx <- fs_separable_fixture()
  rank <- fs_svm_rfe_rank(fx$X, fx$y, cost = 1)
  expect_setequal(rank, colnames(fx$X))
  expect_length(rank, ncol(fx$X))
  ## rank[1] is the last feature standing when the RFE elimination loop
  ## finally stops (most important); rank[length(rank)] was the very first
  ## one dropped (least important) - confirmed against fs_fit_sex()'s own
  ## use of this ranking (`svm_rank[seq_len(curve$best)]` takes a top-k
  ## panel from the FRONT). The two truly informative genes should rank
  ## among the most important - at least one of them, given linear-SVM
  ## weight instability on redundant noise dimensions can occasionally
  ## displace one of two equally-strong true signals.
  expect_true(any(c("G1", "G2") %in% utils::head(rank, 3)))
})

test_that("fs_svm_rfe_curve() reports a CV error curve over panel sizes 1..p with a valid 'best' index", {
  fx <- fs_separable_fixture()
  rank <- fs_svm_rfe_rank(fx$X, fx$y, cost = 1)
  curve <- fs_svm_rfe_curve(fx$X, fx$y, rank, cost = 1, folds = 5)
  expect_equal(curve$k, seq_along(rank))
  expect_true(all(curve$err >= 0 & curve$err <= 1))
  expect_true(curve$best >= 1 && curve$best <= length(rank))
  expect_equal(curve$besterr, min(curve$err))
})

## ---- fs_fit_sex(): full LASSO + RF + SVM-RFE pipeline ----------------------

test_that("fs_fit_sex() (manual mtry/cost, skipping internal CV grids) returns a well-formed result with real signal recovered", {
  fx <- fs_separable_fixture(n_per_group = 10, n_genes = 10)
  fit <- fs_fit_sex(fx$X, fx$y, params = list(
    rf_mtry_mode = "manual", rf_mtry_manual = 3, rf_ntree = 100,
    svm_cost_mode = "manual", svm_cost_manual = 1,
    lasso_cv_folds = 5, rf_cv_folds = 5, svm_cv_folds = 5
  ))
  expect_true(all(c("lasso_genes", "rf_genes", "svm_genes", "sets", "consensus", "consensus_methods") %in% names(fit)))
  expect_equal(fit$n_input, 10L)
  expect_equal(fit$n_samples, 20L)
  expect_setequal(names(fit$sets), c("LASSO", "RandomForest", "SVM_RFE"))
  ## Real, injected signal (G1/G2) should be recoverable by at least the RF importance ranking.
  expect_true(any(c("G1", "G2") %in% fit$rf_genes))
})

test_that("fs_fit_sex() consensus is the intersection of the selected methods (or all three by default)", {
  fx <- fs_separable_fixture(n_per_group = 10, n_genes = 10)
  fit <- fs_fit_sex(fx$X, fx$y, params = list(
    rf_mtry_mode = "manual", rf_mtry_manual = 3, rf_ntree = 100,
    svm_cost_mode = "manual", svm_cost_manual = 1,
    consensus_methods = c("LASSO", "RandomForest")
  ))
  expect_equal(fit$consensus_methods, c("LASSO", "RandomForest"))
  expect_setequal(fit$consensus, intersect(fit$sets$LASSO, fit$sets$RandomForest))
})

test_that("fs_fit_sex() sanitizes non-syntactic gene symbols (e.g. HLA-A) and translates them back on the way out", {
  fx <- fs_separable_fixture(n_per_group = 10, n_genes = 4)
  colnames(fx$X) <- c("HLA-A", "IL-6", "TNF-alpha", "G4")
  fit <- fs_fit_sex(fx$X, fx$y, params = list(
    rf_mtry_mode = "manual", rf_mtry_manual = 2, rf_ntree = 100,
    svm_cost_mode = "manual", svm_cost_manual = 1
  ))
  expect_true(all(fit$svm_rank %in% c("HLA-A", "IL-6", "TNF-alpha", "G4")))
  expect_true(all(names(fit$gini) %in% c("HLA-A", "IL-6", "TNF-alpha", "G4")))
})
