## Module 3 (Multiomics) - the sex-stratified DIABLO/Random-Forest engine
## (multiomics_sexstratified_engine.R), the single implementation behind
## Biomarker Discovery's and Integration's "Sex-Stratified" tabs and the
## natural home for this project's required 4-case sex-stratification
## matrix (Male+Female pooled / Male-only / Female-only / no-sex-metadata).
## A direct, parameterized port of the source pipeline's own leakage-safe
## nested-CV scripts - every claim here is checked against REAL
## limma/mixOmics/randomForest/pROC computation on small synthetic data
## with a known planted signal, never a fabricated AUROC or selection
## frequency. Folds/repeats/ntree/top-K are overridden small for test
## runtime; the pipeline-derived defaults (MSS_DEFAULTS) are left
## untouched by these overrides.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_sexstratified_engine.R"))

## ---- mss_cohort_for_stratum() ------------------------------------------------

test_that("mss_cohort_for_stratum() intersects expr/meth/requested-stratum sample IDs and enforces the minimum stratum size", {
  ids <- paste0("S", 1:20)
  expr <- matrix(rnorm(200), 20, 10, dimnames = list(ids, paste0("g", 1:10)))
  meth <- matrix(rnorm(200), 20, 10, dimnames = list(ids, paste0("cg", 1:10)))
  out <- mss_cohort_for_stratum(expr, meth, NULL, ids[1:15])
  expect_true(out$ok)
  expect_equal(length(out$ids), 15L)

  too_small <- mss_cohort_for_stratum(expr, meth, NULL, ids[1:4])
  expect_false(too_small$ok)
  expect_true(grepl("at least 6", too_small$error))
})

## ---- mss_impute_fold_matrix() (train-only mean, applied to both train/test) ----

test_that("mss_impute_fold_matrix() imputes with the TRAINING fold's own mean, never the test fold's values", {
  X <- matrix(c(1, 2, NA, 100, 5, NA, 7, 200), 4, 2, dimnames = list(paste0("S", 1:4), c("f1", "f2")))
  ## f1: train rows 1-2 = (1,2) -> train_mean=1.5; test rows 3-4 = (NA,100) imputed with 1.5 for the NA.
  ## f2: train rows 1-2 = (5,NA) -> train_mean=5 (na.rm); test rows 3-4 = (7,200), no test NAs to touch.
  out <- mss_impute_fold_matrix(X, train_idx = 1:2, test_idx = 3:4)
  expect_equal(out$train["S2", "f1"], 2)  ## untouched (not NA)
  expect_equal(unname(out$test["S3", "f1"]), 1.5)  ## imputed with TRAIN mean, not test mean
  expect_equal(unname(out$train["S1", "f2"]), 5)  ## train's own NA imputed with train mean
})

test_that("mss_impute_fold_matrix() leaves a column NA in place when its training values are entirely NA (matches the pipeline's own guard)", {
  X <- matrix(c(NA, NA, 3, 4), 2, 2, dimnames = list(c("S1", "S2"), c("f1", "f2")))
  out <- mss_impute_fold_matrix(X, train_idx = 1, test_idx = 2)
  expect_true(is.na(out$train["S1", "f1"]))  ## train_mean for f1 is NA -> "next", left as NA
})

## ---- mss_limma_design() -------------------------------------------------------

test_that("mss_limma_design() builds a 2-column design without a covariate, 3-column with one", {
  outcome <- factor(rep(c("A", "B"), 5))
  d1 <- mss_limma_design(outcome, NULL)
  expect_equal(ncol(d1), 2L)
  cov <- factor(rep(c("X", "Y"), 5))
  d2 <- mss_limma_design(outcome, cov)
  expect_equal(ncol(d2), 3L)
})

## ---- mss_select_features_diablo() / mss_rank_features_pvalue() (real limma) ----

test_that("mss_select_features_diablo()/mss_rank_features_pvalue() (real limma moderated-t) rank a planted differential feature first", {
  set.seed(800)
  n <- 20
  outcome <- factor(rep(c("A", "B"), each = n / 2))
  expr <- matrix(rnorm(n * 10), n, 10, dimnames = list(NULL, paste0("g", 1:10)))
  expr[outcome == "B", 1] <- expr[outcome == "B", 1] + 8  ## g1: huge, unmistakable group difference
  meth <- matrix(rnorm(n * 10), n, 10, dimnames = list(NULL, paste0("cg", 1:10)))
  meth[outcome == "B", 1] <- meth[outcome == "B", 1] + 8

  sel <- mss_select_features_diablo(expr, meth, outcome, NULL, top_expr = 3, top_meth = 3)
  expect_equal(sel$expr[1], "g1")
  expect_equal(sel$meth[1], "cg1")

  ranked <- mss_rank_features_pvalue(expr, outcome, NULL, top_k = 3)
  expect_equal(ranked[1], "g1")
})

## ---- mss_diablo_fold() / mss_rf_fold() (real single-fold fits) ---------------

test_that("mss_diablo_fold() (real mixOmics::block.splsda, one fold) returns a valid held-out score and the in-fold selected features", {
  skip_if_not_installed("mixOmics")
  set.seed(810)
  n <- 20
  ids <- paste0("S", seq_len(n))
  outcome <- factor(rep(c("A", "B"), each = n / 2)); names(outcome) <- ids
  expr <- matrix(rnorm(n * 15), n, 15, dimnames = list(ids, paste0("g", 1:15)))
  expr[outcome == "B", 1:5] <- expr[outcome == "B", 1:5] + 4
  meth <- matrix(rnorm(n * 15), n, 15, dimnames = list(ids, paste0("cg", 1:15)))
  meth[outcome == "B", 1:5] <- meth[outcome == "B", 1:5] - 4

  train_idx <- ids[1:16]; test_idx <- ids[17:20]
  params <- utils::modifyList(MSS_DEFAULTS, list(top_expr = 8, top_meth = 8, keepx_min = 3, keepx_max = 8, ncomp = 1, min_selected_diablo = 3))
  res <- mss_diablo_fold(expr, meth, outcome, NULL, train_idx, test_idx, params)
  expect_true(!is.null(res))
  expect_equal(length(res$score), length(test_idx))
  expect_true(all(res$score %in% c(0, 1)))
  expect_true(length(res$selected$expr) >= 3)
})

test_that("mss_rf_fold() (real randomForest, one fold) returns a valid held-out probability score, and appends the covariate column when supplied", {
  skip_if_not_installed("randomForest")
  set.seed(820)
  n <- 20
  ids <- paste0("S", seq_len(n))
  outcome <- factor(rep(c("A", "B"), each = n / 2)); names(outcome) <- ids
  covariate <- factor(rep(c("X", "Y"), n / 2)); names(covariate) <- ids
  expr <- matrix(rnorm(n * 15), n, 15, dimnames = list(ids, paste0("g", 1:15)))
  expr[outcome == "B", 1:5] <- expr[outcome == "B", 1:5] + 4
  meth <- matrix(rnorm(n * 15), n, 15, dimnames = list(ids, paste0("cg", 1:15)))

  train_idx <- ids[1:16]; test_idx <- ids[17:20]
  params <- utils::modifyList(MSS_DEFAULTS, list(top_expr_rf = 8, top_meth_rf = 8, ntree = 50, min_selected_rf = 3))
  res <- mss_rf_fold(expr, meth, outcome, covariate, train_idx, test_idx, params)
  expect_true(!is.null(res))
  expect_equal(length(res$score), length(test_idx))
  expect_true(all(res$score >= 0 & res$score <= 1))
})

## ---- mss_nested_cv() (real small nested CV, RF engine) -----------------------

mss_planted_dataset <- function(n = 40, seed = 830) {
  set.seed(seed)
  ids <- paste0("S", seq_len(n))
  outcome <- factor(rep(c("HC", "RA"), each = n / 2)); names(outcome) <- ids
  expr <- matrix(rnorm(n * 20), n, 20, dimnames = list(ids, paste0("g", 1:20)))
  expr[outcome == "RA", 1:6] <- expr[outcome == "RA", 1:6] + 3
  meth <- matrix(rnorm(n * 20), n, 20, dimnames = list(ids, paste0("cg", 1:20)))
  meth[outcome == "RA", 1:6] <- meth[outcome == "RA", 1:6] - 3
  list(ids = ids, outcome = outcome, expr = expr, meth = meth)
}

test_that("mss_nested_cv() (real nested CV, RF engine) recovers a better-than-chance AUROC on a planted signal", {
  skip_if_not_installed("randomForest")
  fx <- mss_planted_dataset()
  params <- list(folds = 3, repeats = 1, top_expr_rf = 8, top_meth_rf = 8, ntree = 50, min_selected_rf = 3, min_train_rows_rf = 4)
  res <- mss_nested_cv(fx$expr, fx$meth, fx$outcome, covariate = NULL, engine = "rf", params = params)
  expect_true(res$ok)
  expect_true(res$performance$auroc > 0.6)
  expect_equal(res$performance$engine, "rf")
  expect_equal(nrow(res$scores), sum(!is.na(rowMeans(matrix(res$scores$score, ncol = 1), na.rm = TRUE))))
})

test_that("mss_nested_cv() refuses a non-binary outcome and a too-small minority class", {
  fx <- mss_planted_dataset(n = 30)
  y3 <- factor(rep(c("A", "B", "C"), 10)); names(y3) <- fx$ids
  out_multi <- mss_nested_cv(fx$expr, fx$meth, y3, engine = "rf", params = list(folds = 2, repeats = 1))
  expect_false(out_multi$ok)
  expect_true(grepl("two-class outcome", out_multi$error))

  ## Exactly 3 in the minority class clears the ">= 3" floor (the check is
  ## a strict "< 3" refusal) - 2 is needed to actually trigger it.
  y_imb <- factor(c(rep("A", 28), rep("B", 2))); names(y_imb) <- fx$ids
  out_imb <- mss_nested_cv(fx$expr, fx$meth, y_imb, engine = "rf", params = list(folds = 2, repeats = 1))
  expect_false(out_imb$ok)
  expect_true(grepl("at least 3 are needed per class", out_imb$error))
})

## ---- mss_selection_frequency() -----------------------------------------------

test_that("mss_selection_frequency() tallies real per-fold selected features into a correct frequency table", {
  fold_selected <- list(
    list(expr = c("g1", "g2"), meth = c("cg1")),
    list(expr = c("g1", "g3"), meth = c("cg1", "cg2")),
    list(expr = c("g1"), meth = c("cg2"))
  )
  out <- mss_selection_frequency(fold_selected, block = "expr")
  expect_equal(out$feature[1], "g1")
  expect_equal(out$times_selected[out$feature == "g1"], 3L)
  expect_equal(out$selection_frequency[out$feature == "g1"], 1)
  expect_equal(out$times_selected[out$feature == "g2"], 1L)
})

test_that("mss_selection_frequency() returns NULL with zero folds", {
  expect_null(mss_selection_frequency(list(), block = "expr"))
})

## ---- mss_full_cohort_panel() (real full-cohort RF importance) ----------------

test_that("mss_full_cohort_panel() (real full-cohort Random Forest refit) ranks the planted differential features with the highest importance", {
  skip_if_not_installed("randomForest")
  fx <- mss_planted_dataset(n = 30, seed = 840)
  params <- list(top_expr_rf = 8, top_meth_rf = 8, ntree = 50)
  out <- mss_full_cohort_panel(fx$expr, fx$meth, fx$outcome, covariate = NULL, engine = "rf", params = params)
  expect_true(out$ok)
  expect_true(any(grepl("live extension", out$note)))
  top_feature <- out$panel$feature[which.max(out$panel$importance)]
  expect_true(top_feature %in% paste0("g", 1:6))  ## one of the truly planted-signal genes ranks highest
})

## ---- mss_run_stratified(): the required 4-case sex-stratification matrix ----

mss_sex_dataset <- function(n_per_sex = 20, seed = 850, include_sex = TRUE) {
  set.seed(seed)
  n <- n_per_sex * 2
  ids <- paste0("S", seq_len(n))
  sex <- rep(c("Female", "Male"), each = n_per_sex)
  outcome <- rep(c("HC", "RA"), n / 2)  ## balanced within each sex half
  expr <- matrix(rnorm(n * 20), n, 20, dimnames = list(ids, paste0("g", 1:20)))
  expr[outcome == "RA", 1:6] <- expr[outcome == "RA", 1:6] + 3
  meth <- matrix(rnorm(n * 20), n, 20, dimnames = list(ids, paste0("cg", 1:20)))
  meth[outcome == "RA", 1:6] <- meth[outcome == "RA", 1:6] - 3
  meta <- data.frame(outcome = outcome, stringsAsFactors = FALSE, row.names = ids)
  if (include_sex) meta$sex <- sex
  list(ids = ids, expr = expr, meth = meth, meta = meta)
}

MSS_TEST_PARAMS <- list(folds = 3, repeats = 1, top_expr_rf = 8, top_meth_rf = 8, ntree = 40, min_selected_rf = 3, min_train_rows_rf = 4)

test_that("mss_run_stratified() Case 1: 'pooled' (Male+Female) runs one stratum labeled 'pooled' over all matched samples", {
  skip_if_not_installed("randomForest")
  fx <- mss_sex_dataset()
  out <- mss_run_stratified(fx$expr, fx$meth, fx$meta, "outcome", sex_mode = "pooled", engine = "rf", params = MSS_TEST_PARAMS)
  expect_true(out$ok)
  expect_equal(out$performance$stratum, "pooled")
  expect_equal(out$strata$pooled$n, nrow(fx$meta))
})

test_that("mss_run_stratified() Case 2: 'female' runs on the Female-only subset, detecting the real sex column", {
  skip_if_not_installed("randomForest")
  fx <- mss_sex_dataset()
  out <- mss_run_stratified(fx$expr, fx$meth, fx$meta, "outcome", sex_mode = "female", engine = "rf", params = MSS_TEST_PARAMS)
  expect_true(out$ok)
  expect_equal(out$sex_col, "sex")
  expect_equal(names(out$strata), "female")
  expect_equal(out$strata$female$n, sum(fx$meta$sex == "Female"))
})

test_that("mss_run_stratified() Case 3: 'male' runs on the Male-only subset independently of the female run", {
  skip_if_not_installed("randomForest")
  fx <- mss_sex_dataset()
  out <- mss_run_stratified(fx$expr, fx$meth, fx$meta, "outcome", sex_mode = "male", engine = "rf", params = MSS_TEST_PARAMS)
  expect_true(out$ok)
  expect_equal(names(out$strata), "male")
  expect_equal(out$strata$male$n, sum(fx$meta$sex == "Male"))
})

test_that("mss_run_stratified() Case 4: no sex/gender column in sample_meta -> 'both' mode fails honestly, 'pooled' still works", {
  skip_if_not_installed("randomForest")
  fx <- mss_sex_dataset(include_sex = FALSE)
  out_both <- mss_run_stratified(fx$expr, fx$meth, fx$meta, "outcome", sex_mode = "both", engine = "rf", params = MSS_TEST_PARAMS)
  expect_false(out_both$ok)
  expect_true(grepl("No sex/gender column detected", out_both$error))

  ## Pooled mode never depends on a sex column - the app's own documented
  ## fallback for a dataset with no sex metadata at all.
  out_pooled <- mss_run_stratified(fx$expr, fx$meth, fx$meta, "outcome", sex_mode = "pooled", engine = "rf", params = MSS_TEST_PARAMS)
  expect_true(out_pooled$ok)
  expect_null(out_pooled$sex_col)
})

test_that("mss_run_stratified() 'both' mode (real sex detection) runs Female AND Male as two independent strata with their own performance rows", {
  skip_if_not_installed("randomForest")
  fx <- mss_sex_dataset()
  out <- mss_run_stratified(fx$expr, fx$meth, fx$meta, "outcome", sex_mode = "both", engine = "rf", params = MSS_TEST_PARAMS)
  expect_true(out$ok)
  expect_setequal(tolower(names(out$strata)), c("female", "male"))
  expect_setequal(tolower(out$performance$stratum), c("female", "male"))
})

test_that("mss_run_stratified() requires a valid outcome column", {
  fx <- mss_sex_dataset()
  out <- mss_run_stratified(fx$expr, fx$meth, fx$meta, "not_a_real_column", sex_mode = "pooled", engine = "rf")
  expect_false(out$ok)
})

## ---- mss_panel_wide_by_sex() ---------------------------------------------------

test_that("mss_panel_wide_by_sex() pivots a long per-stratum panel into one row per feature with a column per stratum", {
  panels <- data.frame(
    stratum = c("female", "female", "male"),
    view = c("expression", "methylation", "expression"),
    feature = c("g1", "cg1", "g1"),
    loading = c(0.5, -0.3, 0.7),
    stringsAsFactors = FALSE
  )
  out <- mss_panel_wide_by_sex(panels)
  expect_true(all(c("female", "male") %in% colnames(out)))
  g1_row <- out[out$feature == "g1", ]
  expect_equal(g1_row$female, 0.5)
  expect_equal(g1_row$male, 0.7)
})

test_that("mss_panel_wide_by_sex() returns NULL for an empty/NULL panels table", {
  expect_null(mss_panel_wide_by_sex(NULL))
  expect_null(mss_panel_wide_by_sex(data.frame()))
})
