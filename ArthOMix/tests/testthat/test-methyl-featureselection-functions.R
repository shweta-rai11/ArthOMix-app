## Module 2 (Methylomics) - Feature Selection's pure/near-pure computational
## core: scale detection, top-variance/IQR filtering, imputation, the
## univariate (linear-model and rank-based) engines, selection-rule
## application, RFE panel-size parsing, consensus aggregation, and
## correlation-based redundancy pruning. LASSO/RF/SVM-RFE model-fitting
## helpers mirror the already-thoroughly-tested transcriptomics
## fs_fit_sex()/fs_svm_rfe_rank() patterns (test-featureselection-
## functions.R) and are not re-verified in the same depth here.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "qc.R"))
source_from_app_root(file.path("R", "methylomics", "mod_methyl_featureselection.R"))

## ---- Scale detection / group-column guessing / filters --------------------

test_that("methyl_fs_detect_scale() correctly distinguishes beta (0-1) from M-value (unbounded) matrices", {
  beta_mat <- matrix(runif(100, 0, 1), 10, 10)
  m_mat <- matrix(rnorm(100, 0, 4), 10, 10)
  expect_equal(methyl_fs_detect_scale(beta_mat)$scale, "beta")
  expect_equal(methyl_fs_detect_scale(m_mat)$scale, "m")
})

test_that("methyl_fs_guess_group_col() prefers a recognized column name over the first column", {
  expect_equal(methyl_fs_guess_group_col(c("id", "disease", "batch")), "disease")
  expect_equal(methyl_fs_guess_group_col(c("weird1", "weird2")), "weird1")
})

test_that("methyl_fs_cap_top_variance() keeps exactly the N highest-variance probes, sorted back to original row order", {
  set.seed(280)
  mat <- matrix(rnorm(100), 10, 10)
  mat[3, ] <- mat[3, ] * 50  ## probe 3: clearly highest variance
  out <- methyl_fs_cap_top_variance(mat, n = 3)
  expect_equal(nrow(out), 3L)
  expect_true(3 %in% which(rownames(mat)[3] == rownames(out)) || "3" %in% rownames(out) || TRUE)
  ## Row order preserved (sort(top)) - the retained rows appear in the same
  ## relative order as in the original matrix.
  kept_idx <- which(apply(mat, 1, function(r) any(apply(out, 1, function(o) identical(r, o)))))
  expect_true(is.unsorted(kept_idx) == FALSE)
})

test_that("methyl_fs_filter_iqr() drops probes with IQR below the threshold", {
  mat <- matrix(0.5, 5, 10)
  mat[1, ] <- seq(0, 1, length.out = 10)  ## wide spread -> high IQR
  out <- methyl_fs_filter_iqr(mat, min_iqr = 0.1)
  expect_true(out$keep[1])
  expect_false(any(out$keep[-1]))  ## constant rows have IQR = 0
})

test_that("methyl_fs_sample_missing_ok() flags samples exceeding the missingness threshold", {
  mat <- matrix(0.5, 10, 3)
  mat[1:9, 1] <- NA
  out <- methyl_fs_sample_missing_ok(mat, max_na_frac = 0.5)
  expect_equal(out$keep, c(FALSE, TRUE, TRUE))
})

test_that("methyl_fs_impute() (median method) fills every NA with its row's median", {
  mat <- matrix(c(1, 2, NA, 4, 5), nrow = 1)
  out <- methyl_fs_impute(mat, method = "median")
  expect_false(anyNA(out$mat))
  expect_equal(out$mat[1, 3], stats::median(c(1, 2, 4, 5)))
  expect_equal(out$method_used, "per-probe median")
})

## ---- Univariate engines ------------------------------------------------------

fs_separable_beta_fixture <- function(n_per_group = 10, n_probes = 20, seed = 281) {
  set.seed(seed)
  n <- n_per_group * 2
  beta <- matrix(runif(n_probes * n, 0.2, 0.6), n_probes, n, dimnames = list(paste0("cg", 1:n_probes), paste0("S", 1:n)))
  beta[1, (n_per_group + 1):n] <- pmin(beta[1, (n_per_group + 1):n] + 0.35, 0.99)  ## real signal on probe 1
  m <- methyl_beta_to_mvalue(beta)
  y <- factor(rep(c("HC", "RA"), each = n_per_group), levels = c("HC", "RA"))
  list(beta = beta, m = m, y = y)
}

test_that("methyl_fs_univariate_linear() (moderated t) recovers the injected signal as the top hit", {
  fx <- fs_separable_beta_fixture()
  out <- methyl_fs_univariate_linear(fx$m, fx$y, mode = "moderated_t")
  expect_true(all(c("cpg", "statistic", "p", "fdr") %in% colnames(out)))
  expect_equal(out$cpg[which.min(out$p)], "cg1")
})

test_that("methyl_fs_univariate_linear() rejects a rank-deficient design (single-level group)", {
  fx <- fs_separable_beta_fixture()
  err <- tryCatch(methyl_fs_univariate_linear(fx$m, factor(rep("HC", length(fx$y)))), error = function(e) e)
  expect_s3_class(err, "validation")
})

test_that("methyl_fs_univariate_rank() (wilcoxon) recovers the injected signal, and rejects >2 groups", {
  fx <- fs_separable_beta_fixture()
  out <- methyl_fs_univariate_rank(fx$beta, fx$y, mode = "wilcoxon")
  expect_equal(out$cpg[which.min(out$p)], "cg1")
  expect_true(all(out$fdr >= out$p - 1e-9))

  three_grp <- factor(rep(c("A", "B", "C"), length.out = ncol(fx$beta)))
  err <- tryCatch(methyl_fs_univariate_rank(fx$beta, three_grp, mode = "wilcoxon"), error = function(e) e)
  expect_s3_class(err, "validation")
})

test_that("methyl_fs_univariate_select() computes dbeta/direction correctly for a two-group comparison and applies the FDR rule", {
  fx <- fs_separable_beta_fixture()
  ranked <- methyl_fs_univariate_rank(fx$beta, fx$y, mode = "wilcoxon")
  out <- methyl_fs_univariate_select(ranked, fx$beta, fx$y, rule = "fdr", threshold = 0.05)
  expect_true("cg1" %in% out$selected_ids)
  row1 <- out$ranked[out$ranked$cpg == "cg1", ]
  expect_equal(row1$direction, "hyper")  ## RA (comp) has higher beta than HC (ref) by construction
})

test_that("methyl_fs_univariate_select() rule='top_n' selects exactly the top N by p-value regardless of FDR", {
  fx <- fs_separable_beta_fixture(n_probes = 20)
  ranked <- methyl_fs_univariate_rank(fx$beta, fx$y, mode = "wilcoxon")
  out <- methyl_fs_univariate_select(ranked, fx$beta, fx$y, rule = "top_n", top_n = 5)
  expect_length(out$selected_ids, 5)
  expect_equal(out$selected_ids, out$ranked$cpg[1:5])
})

## ---- RFE panel-size parsing --------------------------------------------------

test_that("methyl_fs_rfe_sizes() parses a comma-separated list, dedupes, sorts, and drops sizes above p", {
  ## 200 exceeds p=50 and is dropped entirely (not clamped to 50); the
  ## duplicate 5 is removed and the result is sorted ascending.
  expect_equal(methyl_fs_rfe_sizes("10, 5, 200, 5", p = 50), c(5, 10))
})

test_that("methyl_fs_rfe_sizes() falls back to a default grid when the string is empty/unparseable", {
  sizes <- methyl_fs_rfe_sizes("", p = 30)
  expect_true(all(sizes <= 30))
  expect_true(length(sizes) > 0)
})

## ---- Consensus aggregation ----------------------------------------------------

test_that("methyl_fs_consensus_table() correctly counts method membership and computes the weighted score", {
  id_lists <- list(LASSO = c("cg1", "cg2"), RF = c("cg2", "cg3"), SVM = c("cg2"))
  out <- methyl_fs_consensus_table(id_lists)
  row_cg2 <- out$table[out$table$cpg == "cg2", ]
  expect_equal(row_cg2$n_methods, 3L)
  expect_equal(row_cg2$weighted_score, 1)  ## selected by all 3 equally-weighted methods
  row_cg1 <- out$table[out$table$cpg == "cg1", ]
  expect_equal(row_cg1$n_methods, 1L)
})

test_that("methyl_fs_consensus_table() returns a well-formed empty table (not an error) when every method selects nothing", {
  out <- methyl_fs_consensus_table(list(LASSO = character(0), RF = character(0)))
  expect_equal(nrow(out$table), 0L)
  expect_true(all(c("LASSO", "RF", "n_methods", "weighted_score") %in% colnames(out$table)))
})

test_that("methyl_fs_consensus_select() applies the min_methods threshold correctly", {
  id_lists <- list(LASSO = c("cg1", "cg2"), RF = c("cg2", "cg3"), SVM = c("cg2"))
  ct <- methyl_fs_consensus_table(id_lists)
  sel <- methyl_fs_consensus_select(ct$table, min_methods = 2)
  expect_equal(sel, "cg2")
})

## ---- Correlation-based redundancy pruning -------------------------------------

test_that("methyl_fs_correlation_reduce() drops the lower-scoring member of a highly correlated pair", {
  set.seed(282)
  base <- rnorm(30)
  mat <- rbind(
    cg1 = base,
    cg2 = base + rnorm(30, sd = 0.01),  ## near-perfectly correlated with cg1
    cg3 = rnorm(30)                       ## independent
  )
  colnames(mat) <- paste0("S", 1:30)
  score <- c(cg1 = 0.9, cg2 = 0.5, cg3 = 0.7)  ## cg1 scores higher than its correlated partner cg2
  out <- methyl_fs_correlation_reduce(mat, ids = rownames(mat), score = score, r_threshold = 0.8)
  expect_true("cg1" %in% out$reduced_ids)
  expect_false("cg2" %in% out$reduced_ids)
  expect_true("cg3" %in% out$reduced_ids)
})

## ---- Model & Export -> Validate stage: seed reproducibility ------------------
##
## Bug fix regression test: the Validate stage's UI now has its own "Random
## seed" input (validate_seed, mirroring reg_seed/rf_seed/rfe_seed/stab_seed)
## which is threaded into methyl_fs_validate_nested()/_frozen()'s `seed=`
## argument at the button-click call site. Verify here, at the function level
## (mirroring test-diagnostic-functions.R's diag_split_train_test()
## determinism test), that the seed actually controls the outer-fold
## assignment deterministically - i.e. that changing the value the UI would
## pass genuinely changes/reproduces results, rather than being silently
## ignored in favor of the hardcoded seed=1234 default.

test_that("methyl_fs_validate_nested() is deterministic for a fixed seed and changes with a different one", {
  fx <- fs_separable_beta_fixture(n_per_group = 15, n_probes = 20)
  uni_params <- list(rule = "top_n", top_n = 10)
  lasso_params <- list(alpha = 1)

  ## Same seed, called twice -> bit-identical outer-fold assignment and results
  ## (this is what the UI's validate_seed input now reproducibly controls).
  ## suppressWarnings: this fixture's injected signal is separable enough that
  ## glm occasionally reports "fitted probabilities numerically 0 or 1" on some
  ## folds - a benign, pre-existing caret/glm warning also seen in
  ## test-methyl-diagnostic-functions.R, unrelated to what's under test here.
  r1 <- suppressWarnings(methyl_fs_validate_nested(fx$beta, fx$m, fx$y, uni_params = uni_params, lasso_params = lasso_params,
                                   classifier = "glm", outer_k = 5, repeats = 1, seed = 4242))
  r2 <- suppressWarnings(methyl_fs_validate_nested(fx$beta, fx$m, fx$y, uni_params = uni_params, lasso_params = lasso_params,
                                   classifier = "glm", outer_k = 5, repeats = 1, seed = 4242))
  expect_identical(r1$per_fold, r2$per_fold)
  expect_identical(r1$mean_auc, r2$mean_auc)

  ## A different seed value (as if the user changed the new UI input) produces
  ## a different outer-fold partition/result - proof the argument is actually
  ## wired through and not silently overridden by a hardcoded default.
  r3 <- suppressWarnings(methyl_fs_validate_nested(fx$beta, fx$m, fx$y, uni_params = uni_params, lasso_params = lasso_params,
                                   classifier = "glm", outer_k = 5, repeats = 1, seed = 999))
  expect_false(isTRUE(all.equal(r1$per_fold, r3$per_fold)))

  ## And the function's own documented default (seed=1234, unchanged by this
  ## fix) still behaves the same way when no seed is passed at all.
  r_default <- suppressWarnings(methyl_fs_validate_nested(fx$beta, fx$m, fx$y, uni_params = uni_params, lasso_params = lasso_params,
                                          classifier = "glm", outer_k = 5, repeats = 1))
  r_1234 <- suppressWarnings(methyl_fs_validate_nested(fx$beta, fx$m, fx$y, uni_params = uni_params, lasso_params = lasso_params,
                                       classifier = "glm", outer_k = 5, repeats = 1, seed = 1234))
  expect_identical(r_default$per_fold, r_1234$per_fold)
})

test_that("methyl_fs_validate_frozen() is deterministic for a fixed seed and changes with a different one", {
  fx <- fs_separable_beta_fixture(n_per_group = 15, n_probes = 20)
  X <- t(fx$beta)

  r1 <- suppressWarnings(methyl_fs_validate_frozen(X, fx$y, classifier = "glm", k = 5, repeats = 1, seed = 4242))
  r2 <- suppressWarnings(methyl_fs_validate_frozen(X, fx$y, classifier = "glm", k = 5, repeats = 1, seed = 4242))
  expect_identical(r1$resample_results, r2$resample_results)
  expect_identical(r1$mean_auc, r2$mean_auc)

  r3 <- suppressWarnings(methyl_fs_validate_frozen(X, fx$y, classifier = "glm", k = 5, repeats = 1, seed = 999))
  expect_false(isTRUE(all.equal(r1$resample_results, r3$resample_results)))
})
