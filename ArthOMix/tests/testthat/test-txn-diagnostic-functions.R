## Module 1 (Transcriptomics) - Diagnostic Model's pure/near-pure
## computational core: z-scoring, class weighting, CV-AUC harness, and the
## full LR+Elastic-Net+RF+SVM per-sex fitting pipeline (diag_fit_sex). RF

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "10_Diagnostic_Model", "mod_diagnostic.R"))

test_that("diag_zrows() z-scores each row independently and zeroes out a constant row", {
  m <- matrix(c(1, 2, 3, 4, 5, 5, 5, 5), nrow = 2, byrow = TRUE)
  out <- diag_zrows(m)
  expect_equal(unname(out[1, ]), as.numeric(scale(c(1, 2, 3, 4))))
  expect_equal(unname(out[2, ]), rep(0, 4))
})

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

test_that("diag_split_train_test() is deterministic and produces a stratified, non-overlapping partition", {
  y <- factor(rep(c("HC", "RA"), each = 20), levels = c("HC", "RA"))
  s1 <- diag_split_train_test(y, test_frac = 0.3)
  s2 <- diag_split_train_test(y, test_frac = 0.3)
  expect_identical(s1, s2)
  expect_length(intersect(s1$train, s1$test), 0)
  expect_equal(length(s1$train) + length(s1$test), length(y))
  expect_setequal(unique(as.character(y[s1$test])), c("HC", "RA"))
})

test_that("diag_cv_auc() returns one AUC per fold, each in [0,1] (or NA), for a simple linear-discriminant refit", {
  set.seed(120)
  n <- 40
  y <- factor(rep(c("HC", "RA"), each = n / 2), levels = c("HC", "RA"))
  X <- matrix(rnorm(n * 3), n, 3)
  X[y == "RA", 1] <- X[y == "RA", 1] + 3
  aucs <- diag_cv_auc(
    X, y, n_folds = 5,
    refit_fn = function(Xtr, ytr) suppressWarnings(stats::glm(ytr ~ ., data = data.frame(ytr, Xtr), family = stats::binomial)),
    predict_fn = function(m, Xte) as.numeric(predict(m, newdata = data.frame(Xte), type = "response"))
  )
  expect_length(aucs, 5)
  expect_true(all(is.na(aucs) | (aucs >= 0 & aucs <= 1)))
  expect_gt(mean(aucs, na.rm = TRUE), 0.7)
})

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

test_that("diag_perf_at_cutoff() computes sensitivity/specificity/accuracy correctly on a crafted example", {
  y <- factor(c("HC", "HC", "RA", "RA", "RA"), levels = c("HC", "RA"))
  prob <- c(0.1, 0.6, 0.9, 0.4, 0.8)
  out <- diag_perf_at_cutoff(prob, y, threshold = 0.5, positive_level = "RA")
  expect_equal(out$sensitivity, 2 / 3)
  expect_equal(out$specificity, 1 / 2)
  expect_equal(out$accuracy, mean(c(FALSE, FALSE, TRUE, TRUE, TRUE)))
})

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
  fx <- diag_separable_fixture(n_per_group = 5)
  err <- tryCatch(
    diag_fit_sex(fx$expr, fx$y, params = list(rf_mtry_mode = "manual", rf_mtry_manual = 2, svm_cost_mode = "manual", svm_cost_manual = 1)),
    error = function(e) e
  )
  expect_s3_class(err, "validation")
  expect_true(grepl("Not enough samples for a train/test split", conditionMessage(err)))
})

noise_leakage_fixture <- function(fixture_seed = 2, n = 60, n_genes = 400) {
  set.seed(fixture_seed)
  y <- factor(rep(c("HC", "RA"), each = n / 2), levels = c("HC", "RA"))
  expr <- matrix(rnorm(n_genes * n), n_genes, n,
                 dimnames = list(paste0("G", seq_len(n_genes)), paste0("S", seq_len(n))))
  list(expr = expr, y = y)
}

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

  expect_gt(frozen_auc, 0.65)
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

test_that("diag_attach_headline() marks the automatic nested-CV AUC as the primary/headline metric only when the run is not leakage-safe and nested-CV succeeded", {
  fit_leaky <- list(leakage_safe = FALSE)
  nested_ok <- list(pooled = list(available = TRUE, auc = 0.55, ci_lo = 0.4, ci_hi = 0.7, n = 40),
                     per_fold = data.frame(), outer_k = 5, n_folds_completed = 5)

  out1 <- diag_attach_headline(fit_leaky, nested_ok)
  expect_equal(out1$headline_metric, "nested_cv")
  expect_identical(out1$nested_cv, nested_ok)

  nested_unavailable <- list(pooled = list(available = FALSE, reason = "not enough genes"), per_fold = data.frame(), outer_k = NA_integer_)
  out2 <- diag_attach_headline(fit_leaky, nested_unavailable)
  expect_equal(out2$headline_metric, "test_split")
  expect_identical(out2$nested_cv, nested_unavailable)

  out3 <- diag_attach_headline(fit_leaky, NULL)
  expect_equal(out3$headline_metric, "test_split")
  expect_null(out3$nested_cv)

  fit_safe <- list(leakage_safe = TRUE)
  out4 <- diag_attach_headline(fit_safe, nested_ok)
  expect_equal(out4$headline_metric, "test_split")
  expect_null(out4$nested_cv)
})
