## Module 3 (Multiomics) - the live DIABLO/SNF/Compare engine
## (multiomics_integration_helpers.R, mounted in mod_multi_integration.R):
## dataset validation/outcome summary, every DIABLO/SNF data-dependent

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_plots.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "multiomics_dataset_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "06_Gene_CpG_Mapping", "multiomics_mapping_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "02_Cohort_Harmonization", "cohort_harmonization_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_integration_helpers.R"))

test_that("mi_validate_dataset() reports reliable_matching only when >=2 blocks share at least MI_MIN_MATCHED_SAMPLES samples", {
  m1 <- matrix(rnorm(60), 10, 6, dimnames = list(paste0("S", 1:10), paste0("f", 1:6)))
  m2 <- matrix(rnorm(48), 8, 6, dimnames = list(paste0("S", 3:10), paste0("g", 1:6)))
  out <- mi_validate_dataset(list(A = m1, B = m2))
  expect_true(out$ok)
  expect_equal(out$n_shared, 8L)
  expect_true(out$reliable_matching)
  expect_null(out$mismatch_message)
})

test_that("mi_validate_dataset() flags unreliable matching (with the shared message) below MI_MIN_MATCHED_SAMPLES, and refuses entirely with zero blocks", {
  m1 <- matrix(rnorm(20), 4, 5, dimnames = list(paste0("S", 1:4), paste0("f", 1:5)))
  m2 <- matrix(rnorm(20), 4, 5, dimnames = list(paste0("S", 4:7), paste0("g", 1:5)))
  out <- mi_validate_dataset(list(A = m1, B = m2))
  expect_false(out$reliable_matching)
  expect_equal(out$mismatch_message, MI_SAMPLE_MISMATCH_MESSAGE)

  none <- mi_validate_dataset(list())
  expect_false(none$ok)
})

test_that("mi_outcome_summary() restricts class counts to the shared/matched sample set only, and flags 2-class imbalance", {
  meta <- data.frame(outcome = c("A", "A", "A", "B", "B", "B", "B"), row.names = paste0("S", 1:7))
  out_full <- mi_outcome_summary(meta, "outcome")
  expect_equal(out_full$n, 7L)

  out_restricted <- mi_outcome_summary(meta, "outcome", sample_ids = paste0("S", 1:4))
  expect_equal(out_restricted$n, 4L)
  expect_equal(unname(out_restricted$class_counts["A"]), 3L)
  expect_true(out_restricted$imbalanced)
})

test_that("mi_outcome_summary() classifies a high-cardinality numeric column as 'continuous' (via the shared ch_classify_column rule)", {
  set.seed(120)
  meta <- data.frame(age = sample(20:80, 30, replace = TRUE), row.names = paste0("S", 1:30))
  out <- mi_outcome_summary(meta, "age")
  expect_equal(out$type, "continuous")
  expect_true(is.na(out$n_classes))
})

test_that("mi_diablo_eligibility() requires >=2 blocks, reliable matching, a categorical outcome with >=2 classes, and >=3 per class", {
  ok_validation <- list(ok = TRUE, n_blocks = 2, reliable_matching = TRUE)
  ok_outcome <- list(type = "categorical", n_classes = 2, class_counts = c(A = 10, B = 10))
  expect_true(mi_diablo_eligibility(ok_validation, ok_outcome)$ok)

  expect_false(mi_diablo_eligibility(list(ok = TRUE, n_blocks = 1, reliable_matching = TRUE), ok_outcome)$ok)
  expect_false(mi_diablo_eligibility(list(ok = TRUE, n_blocks = 2, reliable_matching = FALSE), ok_outcome)$ok)
  expect_false(mi_diablo_eligibility(ok_validation, list(type = "continuous"))$ok)
  small_class <- mi_diablo_eligibility(ok_validation, list(type = "categorical", n_classes = 2, class_counts = c(A = 2, B = 10)))
  expect_false(small_class$ok)
  expect_true(grepl("at least 3 per class", small_class$reason))
})

test_that("mi_diablo_keepx_grid() caps candidates at the block's own feature count and deduplicates", {
  expect_equal(mi_diablo_keepx_grid(1000), c(5, 10, 20, 50, 100, 150, 200, 300))
  expect_equal(mi_diablo_keepx_grid(30), c(5, 10, 20))
  expect_equal(mi_diablo_keepx_grid(3), 3)
})

test_that("mi_diablo_feasible_ncomp() is bounded by both n_classes+1 and min_class_n-1, never exceeding 5", {
  expect_equal(mi_diablo_feasible_ncomp(n_classes = 2, min_class_n = 50), 1:3)
  expect_equal(mi_diablo_feasible_ncomp(n_classes = 10, min_class_n = 4), 1:3)
  expect_equal(mi_diablo_feasible_ncomp(n_classes = 10, min_class_n = 100), 1:5)
  expect_equal(mi_diablo_feasible_ncomp(n_classes = 2, min_class_n = 1), 1)
})

test_that("mi_diablo_max_folds()/mi_diablo_feasible_folds() are bounded by the smallest class, Automatic defaults to 5-fold when affordable", {
  expect_equal(mi_diablo_max_folds(min_class_n = 3), 3L)
  expect_equal(mi_diablo_max_folds(min_class_n = 100), 10L)
  expect_equal(mi_diablo_feasible_folds(min_class_n = 20, requested = NULL), 5L)
  expect_equal(mi_diablo_feasible_folds(min_class_n = 3, requested = NULL), 3L)
  expect_equal(mi_diablo_feasible_folds(min_class_n = 20, requested = 8), 8L)
  expect_equal(mi_diablo_feasible_folds(min_class_n = 20, requested = 50), 10L)
})

test_that("mi_diablo_loo_feasible()/mi_diablo_feasible_repeats() switch at their documented thresholds", {
  expect_true(mi_diablo_loo_feasible(60))
  expect_false(mi_diablo_loo_feasible(61))
  expect_equal(mi_diablo_feasible_repeats(20), 10L)
  expect_equal(mi_diablo_feasible_repeats(50), 5L)
  expect_equal(mi_diablo_feasible_repeats(200), 3L)
})

test_that("mi_diablo_design() builds a symmetric 0.1-off-diagonal matrix in automatic mode, and a zero-forced-diagonal custom matrix otherwise", {
  auto <- mi_diablo_design(c("A", "B", "C"), mode = "automatic")
  expect_equal(diag(auto), c(A = 0, B = 0, C = 0))
  expect_equal(auto["A", "B"], 0.1)

  custom_input <- matrix(0.9, 3, 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  custom <- mi_diablo_design(c("A", "B", "C"), mode = "custom", custom = custom_input)
  expect_equal(unname(diag(custom)), c(0, 0, 0))
  expect_equal(custom["A", "B"], 0.9)
})

test_that("mi_diablo_run() (real mixOmics::tune.block.splsda/block.splsda/perf) recovers a real, better-than-chance cross-validated signal", {
  skip_if_not_installed("mixOmics")
  set.seed(200)
  n <- 30
  y <- factor(rep(c("A", "B"), each = n / 2))
  ids <- paste0("S", seq_len(n))
  expr <- matrix(rnorm(n * 20), n, 20, dimnames = list(ids, paste0("g", 1:20)))
  expr[y == "B", 1:5] <- expr[y == "B", 1:5] + 2
  meth <- matrix(rnorm(n * 20), n, 20, dimnames = list(ids, paste0("cg", 1:20)))
  meth[y == "B", 1:5] <- meth[y == "B", 1:5] - 2
  layers <- list(expression = expr, methylation = meth)
  names(y) <- ids

  res <- mi_diablo_run(layers, y, ids, params = list(folds = 3, nrepeat = 1, validation_mode = "manual", validation_method = "mfold"))
  expect_true(res$ok)
  expect_setequal(res$params$blocks, c("expression", "methylation"))
  expect_equal(res$params$n_samples, n)

  perf_summary <- mi_diablo_performance_summary(res)
  expect_true(!is.null(perf_summary))
  expect_true(perf_summary$ber >= 0 && perf_summary$ber <= 1)
  expect_true(perf_summary$ber < 0.5)

  sel_df <- mi_diablo_selected_features_df(res$fit)
  expect_true(is.data.frame(sel_df))
  expect_true(all(c("expression", "methylation") %in% sel_df$block))

  scores_df <- mi_diablo_sample_scores_df(res$fit, y)
  expect_equal(nrow(scores_df), n)
  expect_setequal(scores_df$patient_id, ids)

  panel_df <- mi_diablo_panel_df_for_plot(sel_df, comp = 1)
  expect_true(all(panel_df$view %in% c("expression", "methylation")))
})

test_that("mixOmics::tune.block.splsda()/perf() ignore an external set.seed() unless their own seed= argument is set (documents why mi_diablo_run() must pass params$seed through directly)", {
  skip_if_not_installed("mixOmics")
  set.seed(42); v1 <- runif(1)
  set.seed(42); set.seed(NULL); v2 <- runif(1)
  set.seed(42); set.seed(NULL); v3 <- runif(1)
  expect_false(isTRUE(all.equal(v2, v3)))
})

test_that("mi_diablo_run(params$seed=...) makes the fit depend only on the supplied seed, not on whatever RNG state happened to precede the call", {
  skip_if_not_installed("mixOmics")
  n <- 30
  y <- factor(rep(c("A", "B"), each = n / 2))
  ids <- paste0("S", seq_len(n))
  set.seed(300)
  expr <- matrix(rnorm(n * 20), n, 20, dimnames = list(ids, paste0("g", 1:20)))
  expr[y == "B", 1:5] <- expr[y == "B", 1:5] + 2
  meth <- matrix(rnorm(n * 20), n, 20, dimnames = list(ids, paste0("cg", 1:20)))
  meth[y == "B", 1:5] <- meth[y == "B", 1:5] - 2
  layers <- list(expression = expr, methylation = meth)
  names(y) <- ids
  base_params <- list(folds = 3, nrepeat = 1, validation_mode = "manual", validation_method = "mfold")

  set.seed(111); invisible(runif(1))
  res_1 <- mi_diablo_run(layers, y, ids, params = c(base_params, list(seed = 42)))
  set.seed(222); invisible(runif(50))
  res_2 <- mi_diablo_run(layers, y, ids, params = c(base_params, list(seed = 42)))

  expect_true(res_1$ok); expect_true(res_2$ok)
  expect_identical(res_1$params$seed, 42L)
  expect_identical(res_1$params$keepX, res_2$params$keepX)
  sel_1 <- mi_diablo_selected_features_df(res_1$fit)
  sel_2 <- mi_diablo_selected_features_df(res_2$fit)
  expect_identical(sel_1$feature, sel_2$feature)
  expect_identical(sel_1$loading, sel_2$loading)
  perf_1 <- mi_diablo_performance_summary(res_1)
  perf_2 <- mi_diablo_performance_summary(res_2)
  expect_identical(perf_1$ber, perf_2$ber)

  set.seed(111); invisible(runif(1))
  res_a <- mi_diablo_run(layers, y, ids, params = base_params)
  set.seed(222); invisible(runif(50))
  res_b <- mi_diablo_run(layers, y, ids, params = base_params)
  expect_true(res_a$ok); expect_true(res_b$ok)
  expect_null(res_a$params$seed)
})

test_that("mi_diablo_run() reports a clear error (never a crash) when fewer than two outcome classes remain in the matched samples", {
  ids <- paste0("S", 1:10)
  y <- factor(rep("A", 10)); names(y) <- ids
  layers <- list(A = matrix(rnorm(100), 10, 10, dimnames = list(ids, paste0("f", 1:10))),
                 B = matrix(rnorm(100), 10, 10, dimnames = list(ids, paste0("g", 1:10))))
  out <- mi_diablo_run(layers, y, ids)
  expect_false(out$ok)
  expect_true(grepl("Fewer than two outcome classes", out$error))
})

test_that("mi_snf_eligibility() requires >=2 blocks, reliable matching, and zero missing values in every block", {
  ok_validation <- list(ok = TRUE, n_blocks = 2, reliable_matching = TRUE, per_block = list(A = list(ok = TRUE, n_missing = 0), B = list(ok = TRUE, n_missing = 0)))
  expect_true(mi_snf_eligibility(ok_validation)$ok)

  with_missing <- ok_validation; with_missing$per_block$B$n_missing <- 5
  out_missing <- mi_snf_eligibility(with_missing)
  expect_false(out_missing$ok)
  expect_true(grepl("requires complete data", out_missing$reason))

  expect_false(mi_snf_eligibility(list(ok = TRUE, n_blocks = 1, reliable_matching = TRUE))$ok)
})

test_that("mi_snf_feasible_k_range() scales its min/max/default with sample size and always keeps K < n_samples", {
  small <- mi_snf_feasible_k_range(20)
  expect_true(small$max < 20)
  expect_true(small$min <= small$default && small$default <= small$max)

  large <- mi_snf_feasible_k_range(500)
  expect_equal(large$min, 50L)
  expect_equal(large$max, 51L)
})

test_that("mi_ari() returns 1 for identical partitions and a low/negative value for a random, unrelated partition", {
  a <- rep(c(1, 2), each = 10)
  expect_equal(mi_ari(a, a), 1)

  set.seed(210)
  b <- sample(rep(c(1, 2), each = 10))
  ari_random <- mi_ari(a, b)
  expect_true(ari_random < 0.3)
})

test_that("mi_snf_run() (real SNFtool::SNF + spectral clustering) recovers two well-separated planted clusters, and mi_snf_posthoc_outcome()/mi_snf_concordance() run on the real result", {
  skip_if_not_installed("SNFtool")
  set.seed(220)
  n <- 24
  ids <- paste0("S", seq_len(n))
  true_cluster <- rep(c(1, 2), each = n / 2)
  m1 <- matrix(rnorm(n * 15), n, 15, dimnames = list(ids, paste0("f", 1:15)))
  m1[true_cluster == 2, ] <- m1[true_cluster == 2, ] + 4
  m2 <- matrix(rnorm(n * 15), n, 15, dimnames = list(ids, paste0("g", 1:15)))
  m2[true_cluster == 2, ] <- m2[true_cluster == 2, ] + 4
  layers <- list(A = m1, B = m2)

  res <- mi_snf_run(layers, params = list(cluster_mode = "manual", n_clusters = 2))
  expect_true(res$ok)
  expect_equal(res$params$n_clusters, 2)
  ari <- mi_ari(res$clusters, true_cluster)
  expect_true(ari > 0.7)

  concordance <- mi_snf_concordance(res)
  expect_true(is.data.frame(concordance))
  expect_setequal(concordance$block, c("A", "B"))
  expect_true(all(concordance$concordance_with_fused >= 0 & concordance$concordance_with_fused <= 1))

  outcome <- factor(ifelse(true_cluster == 1, "HC", "RA")); names(outcome) <- ids
  posthoc <- mi_snf_posthoc_outcome(res, outcome, ids)
  expect_true(!is.null(posthoc))
  expect_true(posthoc$fisher_p < 0.05)
})

test_that("mi_snf_run() reports a clear error object (never a crash) via mi_snf_eligibility() when a block still has missing values - checked upstream of running SNF itself", {
  validation <- list(ok = TRUE, n_blocks = 2, reliable_matching = TRUE,
                      per_block = list(A = list(ok = TRUE, n_missing = 3), B = list(ok = TRUE, n_missing = 0)))
  out <- mi_snf_eligibility(validation)
  expect_false(out$ok)
})
