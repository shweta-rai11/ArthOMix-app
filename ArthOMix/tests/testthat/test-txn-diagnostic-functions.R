## Module 1 (Transcriptomics) - Diagnostic Model's pure/near-pure
## computational core: z-scoring, class weighting, CV-AUC harness, and the
## full LR+Elastic-Net+RF+SVM per-sex fitting pipeline (diag_fit_sex). RF
## mtry and SVM cost use manual params to skip internal CV grid search
## (caret::train/e1071::tune), keeping tests fast without changing which
## code path runs - these "manual" branches are real production code.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "mod_diagnostic.R"))

## ---- diag_zrows() -----------------------------------------------------

test_that("diag_zrows() z-scores each row independently and zeroes out a constant row", {
  m <- matrix(c(1, 2, 3, 4, 5, 5, 5, 5), nrow = 2, byrow = TRUE)
  out <- diag_zrows(m)
  expect_equal(unname(out[1, ]), as.numeric(scale(c(1, 2, 3, 4))))
  expect_equal(unname(out[2, ]), rep(0, 4))
})

## ---- diag_class_weight_levels() / diag_obs_weights() -------------------

test_that("diag_class_weight_levels() mirrors the balanced/manual/equal modes", {
  y <- factor(c(rep("HC", 9), rep("RA", 3)), levels = c("HC", "RA"))
  w <- diag_class_weight_levels(y, "balanced", ratio = NULL)
  expect_equal(unname(w), c(1, 3))
  w2 <- diag_class_weight_levels(y, "manual", ratio = 2.5)
  expect_equal(unname(w2), c(1, 2.5))
  w3 <- diag_class_weight_levels(y, "equal", ratio = NULL)
  expect_equal(unname(w3), c(1, 1))
})

test_that("diag_obs_weights() maps per-observation weights correctly", {
  y <- factor(c("HC", "RA", "HC"), levels = c("HC", "RA"))
  expect_equal(diag_obs_weights(y, "manual", ratio = 4), c(1, 4, 1))
})

## ---- diag_split_train_test() -------------------------------------------

test_that("diag_split_train_test() is deterministic and produces a stratified, non-overlapping partition", {
  y <- factor(rep(c("HC", "RA"), each = 20), levels = c("HC", "RA"))
  s1 <- diag_split_train_test(y, test_frac = 0.3)
  s2 <- diag_split_train_test(y, test_frac = 0.3)
  expect_identical(s1, s2)  ## fixed seed -> reproducible
  expect_length(intersect(s1$train, s1$test), 0)
  expect_equal(length(s1$train) + length(s1$test), length(y))
  ## createDataPartition stratifies by class, so both groups should be represented in test.
  expect_setequal(unique(as.character(y[s1$test])), c("HC", "RA"))
})

## ---- diag_cv_auc() -------------------------------------------------------

test_that("diag_cv_auc() returns one AUC per fold, each in [0,1] (or NA), for a simple linear-discriminant refit", {
  set.seed(120)
  n <- 40
  y <- factor(rep(c("HC", "RA"), each = n / 2), levels = c("HC", "RA"))
  X <- matrix(rnorm(n * 3), n, 3)
  X[y == "RA", 1] <- X[y == "RA", 1] + 3  ## real signal on column 1
  aucs <- diag_cv_auc(
    X, y, n_folds = 5,
    refit_fn = function(Xtr, ytr) suppressWarnings(stats::glm(ytr ~ ., data = data.frame(ytr, Xtr), family = stats::binomial)),
    predict_fn = function(m, Xte) as.numeric(predict(m, newdata = data.frame(Xte), type = "response"))
  )
  expect_length(aucs, 5)
  expect_true(all(is.na(aucs) | (aucs >= 0 & aucs <= 1)))
  expect_gt(mean(aucs, na.rm = TRUE), 0.7)  ## real signal should be clearly detectable
})

## ---- diag_auc_ci() --------------------------------------------------------

test_that("diag_auc_ci() returns a CI that brackets the point AUC estimate", {
  set.seed(121)
  n <- 60
  y <- factor(rep(c("HC", "RA"), each = n / 2), levels = c("HC", "RA"))
  score <- rnorm(n) + ifelse(y == "RA", 2, 0)
  roc_obj <- pROC::roc(y, score, quiet = TRUE, levels = levels(y), direction = "<")
  ci <- diag_auc_ci(roc_obj)
  expect_true(all(c("auc", "lo", "hi") %in% names(ci)))
  expect_true(ci["lo"] <= ci["auc"] && ci["auc"] <= ci["hi"])
})

## ---- diag_perf_at_cutoff() -------------------------------------------------

test_that("diag_perf_at_cutoff() computes sensitivity/specificity/accuracy correctly on a crafted example", {
  y <- factor(c("HC", "HC", "RA", "RA", "RA"), levels = c("HC", "RA"))
  prob <- c(0.1, 0.6, 0.9, 0.4, 0.8)  ## threshold 0.5: predicted pos = obs 2,3,5; RA(3,4,5) obs pos
  out <- diag_perf_at_cutoff(prob, y, threshold = 0.5, positive_level = "RA")
  ## obs_pos = RA at indices 3,4,5; pred_pos = prob>=0.5 at indices 2,3,5
  expect_equal(out$sensitivity, 2 / 3)   ## correctly predicted RA: indices 3,5 of 3 RA
  expect_equal(out$specificity, 1 / 2)   ## correctly predicted non-RA: index 1 of 2 HC (index 2 is a false positive)
  expect_equal(out$accuracy, mean(c(FALSE, FALSE, TRUE, TRUE, TRUE)))
})

## ---- diag_fit_sex(): full LR + Elastic Net + RF + SVM pipeline ------------

diag_separable_fixture <- function(n_per_group = 15, n_genes = 8, seed = 122) {
  set.seed(seed)
  n <- n_per_group * 2
  y <- factor(rep(c("HC", "RA"), each = n_per_group), levels = c("HC", "RA"))
  expr <- matrix(rnorm(n_genes * n, mean = 8, sd = 1), n_genes, n, dimnames = list(paste0("G", 1:n_genes), paste0("S", 1:n)))
  expr[1, y == "RA"] <- expr[1, y == "RA"] + 5
  expr[2, y == "RA"] <- expr[2, y == "RA"] + 4
  list(expr = expr, y = y)
}

test_that("diag_fit_sex() (manual mtry/cost) returns all four model fits with valid AUCs and a held-out test evaluation", {
  fx <- diag_separable_fixture()
  fit <- diag_fit_sex(fx$expr, fx$y, params = list(
    rf_mtry_mode = "manual", rf_mtry_manual = 3, rf_ntree = 100,
    svm_cost_mode = "manual", svm_cost_manual = 1, test_frac = 0.3
  ))
  for (model_key in c("lr", "enet", "rf", "svm")) {
    m <- fit[[model_key]]
    expect_true(is.finite(m$full_auc) && m$full_auc >= 0 && m$full_auc <= 1, info = model_key)
    expect_true(isTRUE(m$test$available), info = model_key)
    expect_true(m$test$auc >= 0 && m$test$auc <= 1, info = model_key)
  }
  expect_equal(fit$n_input, 8L)
  expect_equal(fit$n_samples + fit$n_test, 30L)
})

test_that("diag_fit_sex() rejects a sex pool too small for the requested train/test split", {
  fx <- diag_separable_fixture(n_per_group = 5)  ## 10 total; 0.3 test -> ~7 train/3 test, below the 10-train floor
  err <- tryCatch(
    diag_fit_sex(fx$expr, fx$y, params = list(rf_mtry_mode = "manual", rf_mtry_manual = 2, svm_cost_mode = "manual", svm_cost_manual = 1)),
    error = function(e) e
  )
  expect_s3_class(err, "validation")
  expect_true(grepl("Not enough samples for a train/test split", conditionMessage(err)))
})

## ---- diag_validate_nested(): leakage-safe nested-CV validation ------------
##
## The Diagnostic Model's gene panel is chosen upstream (Candidate Gene
## Identification + Feature Selection's LASSO/RF/SVM-RFE consensus) using
## label information from every sample in the pool, with no train/test split
## anywhere in that discovery chain (Ambroise & McLachlan 2002, PNAS: feature-
## selection-before-split leakage). diag_validate_nested() fixes this for the
## Univariate + LASSO slice of that chain by reselecting the panel inside each
## outer training fold instead of once on the full pool.
##
## A pure-noise fixture (no real signal at all, more candidate genes than
## samples) makes the leakage's effect large and reliable to detect: with
## enough noise genes, *some* gene looks spuriously associated with the
## outcome purely by chance when ranked across the FULL sample pool - the
## textbook Ambroise-McLachlan "selection bias" demonstration. A "frozen"
## fixed panel selected once on the full pool, then merely refit per CV fold
## via the pipeline's own diag_cv_auc() (exactly the shape of the pre-existing
## Test-split/CV-AUC computation this fix addresses), reproduces that
## inflation. diag_validate_nested(), which reselects the panel using only
## each fold's own training labels, should score close to chance instead.

noise_leakage_fixture <- function(fixture_seed = 2, n = 60, n_genes = 400) {
  set.seed(fixture_seed)
  y <- factor(rep(c("HC", "RA"), each = n / 2), levels = c("HC", "RA"))
  expr <- matrix(rnorm(n_genes * n), n_genes, n,
                 dimnames = list(paste0("G", seq_len(n_genes)), paste0("S", seq_len(n))))
  list(expr = expr, y = y)
}

## Reproduces the leaky architecture being fixed: panel selected ONCE using
## every sample's label (full-pool univariate ranking + LASSO, mirroring
## fs_fit_sex()'s own LASSO step), then that fixed panel is merely refit per
## fold by the pipeline's existing diag_cv_auc() harness - i.e. exactly what
## Model Testing/Training already do to a panel chosen upstream on the full
## pool.
leaky_frozen_panel_cv_auc <- function(expr, y, uni_top_n = 100, seed = 1234) {
  design <- stats::model.matrix(~y)
  uni_fit <- limma::eBayes(limma::lmFit(expr, design))
  tt <- limma::topTable(uni_fit, coef = 2, number = Inf, sort.by = "P")
  uni_genes <- rownames(tt)[seq_len(min(uni_top_n, nrow(tt)))]
  Xall <- t(expr[uni_genes, , drop = FALSE])
  cv <- glmnet::cv.glmnet(Xall, y, family = "binomial", alpha = 1, nfolds = 5)
  co <- coef(cv, s = "lambda.min")[-1, 1, drop = TRUE]
  panel <- names(co)[co != 0]
  if (length(panel) < 1) panel <- utils::head(uni_genes, 5)
  Xraw <- t(expr[panel, , drop = FALSE])
  mean(diag_cv_auc(
    Xraw, y, n_folds = 5,
    refit_fn = function(Xtr, ytr) suppressWarnings(stats::glm(ytr ~ ., data = data.frame(ytr, Xtr), family = stats::binomial)),
    predict_fn = function(m, Xte) as.numeric(predict(m, newdata = data.frame(Xte), type = "response")),
    seed = seed
  ), na.rm = TRUE)
}

test_that("diag_validate_nested() scores a pure-noise panel close to chance, well below the leaky frozen-panel CV-AUC", {
  fx <- noise_leakage_fixture()
  frozen_auc <- leaky_frozen_panel_cv_auc(fx$expr, fx$y)
  nested <- diag_validate_nested(fx$expr, fx$y, outer_k = 5)

  expect_true(isTRUE(nested$pooled$available))
  nested_auc <- nested$pooled$auc

  ## The leaky, full-pool-selected panel is clearly inflated on pure noise...
  expect_gt(frozen_auc, 0.65)
  ## ...while reselecting the panel inside each training fold only pulls the
  ## honest estimate detectably back down toward chance (0.5), with a wide
  ## safety margin under the fixture's actual empirical gap (~0.44).
  expect_lt(nested_auc, 0.75)
  expect_gt(frozen_auc - nested_auc, 0.15)
})

test_that("diag_validate_nested() is deterministic given a fixed seed", {
  fx <- noise_leakage_fixture(fixture_seed = 100)
  r1 <- diag_validate_nested(fx$expr, fx$y, outer_k = 5, seed = 4242)
  r2 <- diag_validate_nested(fx$expr, fx$y, outer_k = 5, seed = 4242)
  expect_identical(r1$pooled$auc, r2$pooled$auc)
  expect_identical(r1$per_fold, r2$per_fold)
})

test_that("diag_validate_nested() fold-specific feature selection only ever sees that fold's training labels", {
  ## Wraps limma::lmFit so every call it makes records exactly which sample
  ## columns (by name) it was given. If diag_validate_nested() ever leaked
  ## test-fold samples into the per-fold univariate ranking step, at least
  ## one call would include every sample - the FULL pool - not a strict
  ## subset of it.
  fx <- noise_leakage_fixture(fixture_seed = 7, n = 40, n_genes = 60)
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
  invisible(diag_validate_nested(fx$expr, fx$y, outer_k = 4, seed = 55))
  expect_gt(length(seen_cols$calls), 0)
  full_pool <- colnames(fx$expr)
  for (cols in seen_cols$calls) {
    expect_lt(length(cols), length(full_pool))
    expect_true(all(cols %in% full_pool))
  }
})
