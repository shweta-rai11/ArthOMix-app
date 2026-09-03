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

## dxm_validate_nested() / dxm_attach_headline(): the leakage-safe headline
## metric added for the Diagnostic Classifier's default (no genuine held-out
## split) path. X is samples (rows) x CpGs (columns), matching dxm$full_X's
## own orientation in this module (the opposite of transcriptomics' expr
## matrix, which is genes x samples).
dxm_noise_leakage_fixture <- function(fixture_seed = 2, n = 60, n_cpgs = 400) {
  set.seed(fixture_seed)
  y <- factor(rep(c(DXM_NEG, DXM_POS), each = n / 2), levels = c(DXM_NEG, DXM_POS))
  m <- matrix(rnorm(n_cpgs * n), n, n_cpgs,
              dimnames = list(paste0("S", seq_len(n)), paste0("cg", seq_len(n_cpgs))))
  list(X = as.data.frame(m), y = y)
}

dxm_leaky_frozen_panel_cv_auc <- function(X, y, uni_top_n = 100, seed = 1234) {
  m_mat <- t(as.matrix(X))
  design <- stats::model.matrix(~y)
  uni_fit <- limma::eBayes(limma::lmFit(m_mat, design))
  tt <- limma::topTable(uni_fit, coef = 2, number = Inf, sort.by = "P")
  uni_cpgs <- rownames(tt)[seq_len(min(uni_top_n, nrow(tt)))]
  Xall <- as.matrix(X[, uni_cpgs, drop = FALSE])
  cv <- glmnet::cv.glmnet(Xall, y, family = "binomial", alpha = 1, nfolds = 5)
  co <- stats::coef(cv, s = "lambda.min")[-1, 1, drop = TRUE]
  panel <- names(co)[co != 0]
  if (length(panel) < 1) panel <- utils::head(uni_cpgs, 5)

  set.seed(seed)
  folds <- caret::createFolds(y, k = 5)
  aucs <- vapply(folds, function(te) {
    tr <- setdiff(seq_along(y), te)
    fit_df <- data.frame(y = y[tr], X[tr, panel, drop = FALSE], check.names = FALSE)
    model <- suppressWarnings(stats::glm(y ~ ., data = fit_df, family = stats::binomial))
    pred <- as.numeric(stats::predict(model, newdata = data.frame(X[te, panel, drop = FALSE], check.names = FALSE), type = "response"))
    rb <- dxm_roc_bundle(y[te], pred)
    if (is.null(rb)) NA_real_ else rb$auc
  }, numeric(1))
  mean(aucs, na.rm = TRUE)
}

test_that("dxm_validate_nested() scores a pure-noise CpG panel close to chance, well below the leaky frozen-panel CV-AUC", {
  fx <- dxm_noise_leakage_fixture()
  frozen_auc <- dxm_leaky_frozen_panel_cv_auc(fx$X, fx$y)
  nested <- dxm_validate_nested(fx$X, fx$y, outer_k = 5)

  expect_true(isTRUE(nested$pooled$available))
  nested_auc <- nested$pooled$auc

  expect_gt(frozen_auc, 0.65)
  expect_lt(nested_auc, 0.75)
  expect_gt(frozen_auc - nested_auc, 0.15)
})

test_that("dxm_validate_nested() is deterministic given a fixed seed", {
  fx <- dxm_noise_leakage_fixture(fixture_seed = 100)
  r1 <- dxm_validate_nested(fx$X, fx$y, outer_k = 5, seed = 4242)
  r2 <- dxm_validate_nested(fx$X, fx$y, outer_k = 5, seed = 4242)
  expect_identical(r1$pooled$auc, r2$pooled$auc)
  expect_identical(r1$per_fold, r2$per_fold)
})

test_that("dxm_validate_nested() fold-specific feature selection only ever sees that fold's training samples", {
  fx <- dxm_noise_leakage_fixture(fixture_seed = 7, n = 40, n_cpgs = 60)
  seen_cols <- new.env()
  seen_cols$calls <- list()
  real_lmFit <- limma::lmFit
  testthat::local_mocked_bindings(
    lmFit = function(object, ...) {
      seen_cols$calls[[length(seen_cols$calls) + 1]] <- colnames(object)
      real_lmFit(object, ...)
    },
    .package = "limma"
  )
  invisible(dxm_validate_nested(fx$X, fx$y, outer_k = 4, seed = 55))
  expect_gt(length(seen_cols$calls), 0)
  full_pool <- rownames(fx$X)
  for (cols in seen_cols$calls) {
    expect_lt(length(cols), length(full_pool))
    expect_true(all(cols %in% full_pool))
  }
})

test_that("dxm_attach_headline() marks the automatic nested-CV AUC as primary only when not leakage-safe and nested-CV succeeded", {
  nested_ok <- list(pooled = list(available = TRUE, auc = 0.55, ci_lo = 0.4, ci_hi = 0.7, n = 40),
                     per_fold = data.frame(), outer_k = 5, n_folds_completed = 5)

  out1 <- dxm_attach_headline(FALSE, nested_ok)
  expect_equal(out1$headline_metric, "nested_cv")
  expect_identical(out1$nested_cv, nested_ok)

  nested_unavailable <- list(pooled = list(available = FALSE, reason = "not enough CpGs"), per_fold = data.frame(), outer_k = NA_integer_)
  out2 <- dxm_attach_headline(FALSE, nested_unavailable)
  expect_equal(out2$headline_metric, "test_split")
  expect_identical(out2$nested_cv, nested_unavailable)

  out3 <- dxm_attach_headline(FALSE, NULL)
  expect_equal(out3$headline_metric, "test_split")
  expect_null(out3$nested_cv)

  out4 <- dxm_attach_headline(TRUE, nested_ok)
  expect_equal(out4$headline_metric, "test_split")
  expect_null(out4$nested_cv)
})
