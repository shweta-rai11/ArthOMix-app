## Module 2 (Methylomics) - Diagnostic Classifier's pure/near-pure
## functions: sex-label helpers, the pre-flight data-quality checklist,
## feature-set intersection, ROC bundling, threshold-picking strategies,

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "13_Diagnostic_Classifier", "mod_methyl_diagnostic.R"))

test_that("dxm_sex_label() maps the three recognized selections and title-cases anything else", {
  expect_equal(dxm_sex_label("all"), "All samples")
  expect_equal(dxm_sex_label("female"), "Female")
  expect_equal(dxm_sex_label("male"), "Male")
  expect_equal(dxm_sex_label(NULL), "Female")
})

test_that("dxm_normalize_sex() maps free-text F/M spellings to a single-letter code, NA otherwise", {
  expect_equal(dxm_normalize_sex(c("Female", "f", "M", "male", "unknown")), c("F", "F", "M", "M", NA_character_))
})

test_that("dxm_intersect_features() reports shared and unmatched feature IDs correctly", {
  out <- dxm_intersect_features(c("cg1", "cg2", "cg3"), c("cg2", "cg3", "cg4"))
  expect_setequal(out$shared, c("cg2", "cg3"))
  expect_setequal(out$unmatched, "cg1")
})

test_that("dxm_validate_checklist() flags duplicate CpGs, missing values, and class imbalance", {
  X <- data.frame(cg1 = c(1, 2, 3, 4, 5, 6), cg1b = c(1, 1, 1, 1, 1, 1))
  colnames(X) <- c("cg1", "cg1")
  X[1, 1] <- NA
  y <- factor(c(rep("Class0", 5), "Class1"), levels = c("Class0", "Class1"))
  dxm <- list(train_X = X, train_y = y, test_internal_X = X[1:2, , drop = FALSE], test_internal_y = y[1:2],
               ref_level = "HC", comp_level = "RA")
  chk <- dxm_validate_checklist(dxm)
  expect_equal(chk$Status[chk$Check == "Duplicated CpGs"], "WARN")
  expect_equal(chk$Status[chk$Check == "Missing methylation values"], "WARN")
  expect_equal(chk$Status[chk$Check == "Class balance (training)"], "WARN")
})

dxm_prob_fixture <- function(n_per_group = 20, seed = 290) {
  set.seed(seed)
  y <- factor(rep(c("Class0", "Class1"), each = n_per_group), levels = c("Class0", "Class1"))
  prob <- c(runif(n_per_group, 0.1, 0.5), runif(n_per_group, 0.5, 0.9))
  list(y = y, prob = prob)
}

test_that("dxm_roc_bundle() computes a real AUC with a CI bracketing it", {
  fx <- dxm_prob_fixture()
  rb <- dxm_roc_bundle(fx$y, fx$prob)
  expect_true(rb$auc > 0.9)
  expect_true(rb$ci_lo <= rb$auc && rb$auc <= rb$ci_hi)
})

test_that("dxm_pick_threshold() 'youden' picks the threshold maximizing sensitivity+specificity-1", {
  fx <- dxm_prob_fixture()
  rb <- dxm_roc_bundle(fx$y, fx$prob)
  thr <- dxm_pick_threshold("youden", rb)
  expect_true(thr > 0 && thr < 1)
  conf <- dxm_confusion(fx$y, fx$prob, thr)
  expect_gt(conf$balanced_accuracy, 0.85)
})

test_that("dxm_pick_threshold() falls back to 0.5 when roc_bundle is NULL", {
  expect_equal(dxm_pick_threshold("youden", NULL), 0.5)
})

test_that("dxm_confusion() computes sensitivity/specificity/precision/F1/MCC correctly on a hand-worked example", {
  y <- factor(c("Class0", "Class0", "Class1", "Class1", "Class1"), levels = c("Class0", "Class1"))
  prob <- c(0.1, 0.6, 0.9, 0.4, 0.8)
  conf <- dxm_confusion(y, prob, threshold = 0.5)
  expect_equal(unname(conf$TP), 2)
  expect_equal(unname(conf$FP), 1)
  expect_equal(unname(conf$FN), 1)
  expect_equal(unname(conf$TN), 1)
  expect_equal(conf$sensitivity, 2 / 3)
  expect_equal(conf$specificity, 1 / 2)
})

test_that("dxm_metrics_bundle() combines confusion-matrix metrics with AUC/CI/Brier/PR-AUC", {
  fx <- dxm_prob_fixture()
  rb <- dxm_roc_bundle(fx$y, fx$prob)
  m <- dxm_metrics_bundle(fx$y, fx$prob, threshold = 0.5, roc_bundle = rb)
  expect_true(all(c("sensitivity", "specificity", "auc", "brier", "pr_auc") %in% names(m)))
  expect_true(m$brier >= 0 && m$brier <= 1)
  expect_equal(m$auc, rb$auc)
})

test_that("dxm_calibration() produces a monotonic-ish binned calibration table and a valid Brier score", {
  fx <- dxm_prob_fixture(n_per_group = 50)
  cal <- dxm_calibration(fx$y, fx$prob, bins = 5)
  expect_true(all(c("mean_pred", "mean_obs", "n") %in% colnames(cal$table)))
  expect_true(cal$brier >= 0 && cal$brier <= 1)
  expect_true(sum(cal$table$n) == length(fx$prob))
})

test_that("dxm_overfitting_note() flags a large train-vs-test AUC gap as possible overfitting", {
  note_overfit <- dxm_overfitting_note(train_auc = 0.99, cv_auc = 0.95, test_auc = 0.6)
  expect_true(grepl("overfitting", note_overfit))
  note_ok <- dxm_overfitting_note(train_auc = 0.85, cv_auc = 0.83, test_auc = 0.82)
  expect_true(grepl("no strong evidence of overfitting", note_ok))
  expect_equal(dxm_overfitting_note(NA_real_, 0.8, 0.8), "Not enough completed evaluations yet to assess overfitting.")
})
