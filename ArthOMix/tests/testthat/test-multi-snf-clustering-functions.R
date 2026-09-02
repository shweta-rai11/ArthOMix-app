## Module 3 (Multiomics) - SNF Clustering / Patient Stratification's own
## pure functions (snf_clustering_helpers.R): single-block validation/
## eligibility (the "Single-Omics Clustering" fallback mi_validate_dataset()
## doesn't handle), data-type detection + explicit preprocessing chain,
## the single-omics SNF fallback (real SNFtool affinity/spectral clustering,
## no fusion step) alongside delegation to the already-tested mi_snf_run()
## for >=2 blocks, clinical-variable detection, real categorical/continuous/
## survival association tests, resampling-based stability (real repeated
## SNF reruns + ARI), parameter sensitivity, and per-feature cluster
## association ranking - never a fabricated stability/sensitivity verdict.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_plots.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "multiomics_dataset_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "06_Gene_CpG_Concordance", "multiomics_concordance_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "02_Cohort_Harmonization", "cohort_harmonization_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_integration_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "04_SNF_Clustering", "snf_clustering_helpers.R"))

## ---- sfc_validate_dataset() / sfc_eligibility() -----------------------------

test_that("sfc_validate_dataset() handles the single-block case directly (mi_validate_dataset() itself requires >=2)", {
  mat <- matrix(rnorm(100), 20, 5, dimnames = list(paste0("S", 1:20), paste0("f", 1:5)))
  out <- sfc_validate_dataset(list(A = mat))
  expect_true(out$ok)
  expect_equal(out$n_blocks, 1L)
  expect_true(out$reliable_matching)
  expect_equal(out$n_shared, 20L)
})

test_that("sfc_validate_dataset() flags too few samples in a single block, and delegates to mi_validate_dataset() for >=2 blocks", {
  small <- matrix(rnorm(10), 2, 5, dimnames = list(c("S1", "S2"), paste0("f", 1:5)))
  out_small <- sfc_validate_dataset(list(A = small))
  expect_false(out_small$reliable_matching)

  m1 <- matrix(rnorm(100), 20, 5, dimnames = list(paste0("S", 1:20), paste0("f", 1:5)))
  m2 <- matrix(rnorm(100), 20, 5, dimnames = list(paste0("S", 1:20), paste0("g", 1:5)))
  out_multi <- sfc_validate_dataset(list(A = m1, B = m2))
  expect_equal(out_multi$n_blocks, 2L)
})

test_that("sfc_eligibility() resolves 'single_omics' vs 'multi_omics_snf' mode, and refuses when missing values remain", {
  ok_single <- list(ok = TRUE, n_blocks = 1, reliable_matching = TRUE, per_block = list(A = list(ok = TRUE, n_missing = 0)))
  expect_equal(sfc_eligibility(ok_single)$mode, "single_omics")

  ok_multi <- list(ok = TRUE, n_blocks = 2, reliable_matching = TRUE, per_block = list(A = list(ok = TRUE, n_missing = 0), B = list(ok = TRUE, n_missing = 0)))
  expect_equal(sfc_eligibility(ok_multi)$mode, "multi_omics_snf")

  with_missing <- list(ok = TRUE, n_blocks = 1, reliable_matching = TRUE, per_block = list(A = list(ok = TRUE, n_missing = 3)))
  out <- sfc_eligibility(with_missing)
  expect_false(out$ok)
  expect_true(grepl("cannot proceed until they are handled", out$reason))
})

## ---- sfc_detect_data_type() / sfc_transform_choices() -----------------------

test_that("sfc_detect_data_type() distinguishes binary/proportion/count/continuous shapes", {
  expect_equal(sfc_detect_data_type(matrix(c(0, 1, 0, 1), 2, 2))$type, "binary")
  expect_equal(sfc_detect_data_type(matrix(runif(20, 0, 1), 4, 5))$type, "proportion_0_1")
  expect_equal(sfc_detect_data_type(matrix(rpois(20, 50), 4, 5))$type, "count")
  expect_equal(sfc_detect_data_type(matrix(rnorm(20, 5, 3), 4, 5))$type, "continuous")
})

test_that("sfc_transform_choices() uses the richer MULTI_LIVE_NORM_CHOICES when omics_type is declared, the value-shape heuristic otherwise", {
  mat <- matrix(runif(20, 0, 1), 4, 5)
  declared <- sfc_transform_choices(mat, omics_type = "methylation")
  expect_equal(declared$type, "methylation")
  expect_equal(declared$choices, MULTI_LIVE_NORM_CHOICES[["methylation"]])

  undeclared <- sfc_transform_choices(mat, omics_type = NULL)
  expect_equal(undeclared$type, "proportion_0_1")
  expect_true("mvalue" %in% names(undeclared$choices))
})

## ---- sfc_preprocess_block() (real chained missing/normalize/filter) --------

test_that("sfc_preprocess_block() hard-stops when missing values exist and missing_method='none' (SNF hard-requires complete data)", {
  mat <- matrix(c(1, NA, 3, 4), 2, 2)
  out <- sfc_preprocess_block(mat, missing_method = "none")
  expect_false(out$ok)
  expect_true(grepl("no missing-value handling was selected", out$error))
})

test_that("sfc_preprocess_block() chains missing-value handling, transform, and feature filtering in order, logging each step", {
  set.seed(600)
  mat <- matrix(rpois(20 * 10, 100), 20, 10, dimnames = list(paste0("S", 1:20), paste0("f", 1:10)))
  mat[1, 1] <- NA  ## one isolated missing cell, well under any drop threshold
  out <- sfc_preprocess_block(mat, transform = "log2", missing_method = "mean", filter_criterion = "variance", filter_top_n = 5)
  expect_true(out$ok)
  expect_equal(ncol(out$mat), 5L)
  expect_false(anyNA(out$mat))
  expect_true(any(grepl("Missing values: mean", out$log)))
  expect_true(any(grepl("Transformation applied: log2", out$log)))
  expect_true(any(grepl("Feature filter \\(variance\\): kept top 5 of 10", out$log)))
})

test_that("sfc_preprocess_block() refuses when too few samples/features remain after preprocessing", {
  mat <- matrix(rnorm(4), 2, 2, dimnames = list(c("S1", "S2"), c("f1", "f2")))
  out <- sfc_preprocess_block(mat, missing_method = "none")
  expect_false(out$ok)  ## fewer than MI_MIN_MATCHED_SAMPLES rows, no NAs so it reaches the size check
  expect_true(grepl("Too few samples or features", out$error))
})

## ---- sfc_snf_run() single-omics fallback (real SNFtool) + multi-block delegation ----

test_that("sfc_snf_run() single-block fallback (real SNFtool affinity+spectral clustering, no fusion) recovers a planted 2-cluster structure", {
  skip_if_not_installed("SNFtool")
  set.seed(610)
  n <- 24
  ids <- paste0("S", seq_len(n))
  true_cluster <- rep(c(1, 2), each = n / 2)
  m1 <- matrix(rnorm(n * 15), n, 15, dimnames = list(ids, paste0("f", 1:15)))
  m1[true_cluster == 2, ] <- m1[true_cluster == 2, ] + 4

  res <- sfc_snf_run(list(A = m1), params = list(cluster_mode = "manual", n_clusters = 2))
  expect_true(res$ok)
  expect_equal(res$params$mode, "single_omics")
  expect_true(is.na(res$params$t))  ## no fusion step for a single block
  ari <- mi_ari(res$clusters, true_cluster)
  expect_true(ari > 0.7)
})

test_that("sfc_snf_run() delegates to mi_snf_run() for >=2 blocks and tags params$mode = 'multi_omics_snf'", {
  set.seed(620)
  n <- 20
  ids <- paste0("S", seq_len(n))
  m1 <- matrix(rnorm(n * 10), n, 10, dimnames = list(ids, paste0("f", 1:10)))
  m2 <- matrix(rnorm(n * 10), n, 10, dimnames = list(ids, paste0("g", 1:10)))
  res <- sfc_snf_run(list(A = m1, B = m2), params = list(cluster_mode = "manual", n_clusters = 2))
  expect_true(res$ok)
  expect_equal(res$params$mode, "multi_omics_snf")
})

test_that("sfc_snf_run() refuses with zero omics blocks", {
  out <- sfc_snf_run(list())
  expect_false(out$ok)
})

## ---- sfc_detect_clinical() ---------------------------------------------------

test_that("sfc_detect_clinical() classifies categorical vs. continuous columns and detects a real survival pair by name+shape", {
  meta <- data.frame(
    sex = rep(c("F", "M"), 10), age = rnorm(20, 55, 10),
    os_time = rexp(20, 0.01), os_event = rep(c(0, 1), 10),
    row.names = paste0("S", 1:20)
  )
  out <- sfc_detect_clinical(meta)
  expect_true("sex" %in% out$categorical)
  expect_true("age" %in% out$continuous)
  expect_equal(out$survival$time_col, "os_time")
  expect_equal(out$survival$event_col, "os_event")
  ## Survival columns are excluded from the plain categorical/continuous lists once claimed.
  expect_false("os_time" %in% out$continuous)
  expect_false("os_event" %in% out$categorical)
})

test_that("sfc_detect_clinical() reports survival unavailable when only a time-like OR only an event-like column exists, never guessing the other half", {
  meta_time_only <- data.frame(duration = rexp(10, 0.01), row.names = paste0("S", 1:10))
  expect_null(sfc_detect_clinical(meta_time_only)$survival)

  meta_event_only <- data.frame(death = rep(c(0, 1), 5), row.names = paste0("S", 1:10))
  expect_null(sfc_detect_clinical(meta_event_only)$survival)
})

## ---- sfc_test_categorical() / sfc_test_continuous() / sfc_test_survival() (real tests) ----

test_that("sfc_test_categorical() runs a real Fisher's exact test and reports Cramer's V", {
  set.seed(630)
  ids <- paste0("S", 1:20)
  clusters <- setNames(rep(c(1, 2), each = 10), ids)
  x <- setNames(c(rep("A", 9), "B", rep("B", 9), "A"), ids)  ## strongly associated with cluster
  out <- sfc_test_categorical(clusters, x)
  expect_true(out$ok)
  expect_equal(out$test, "Fisher's exact test")
  expect_true(out$p_value < 0.01)
})

test_that("sfc_test_continuous() runs a real Kruskal-Wallis test and reports per-cluster medians", {
  set.seed(640)
  ids <- paste0("S", 1:30)
  clusters <- setNames(rep(c(1, 2, 3), each = 10), ids)
  x <- setNames(c(rnorm(10, 0), rnorm(10, 10), rnorm(10, 20)), ids)  ## clearly separated by cluster
  out <- sfc_test_continuous(clusters, x)
  expect_true(out$ok)
  expect_true(out$p_value < 0.001)
  expect_equal(nrow(out$summary), 3L)
})

test_that("sfc_test_survival() runs a real Kaplan-Meier/log-rank/Cox fit and detects a strong survival difference between 2 clusters", {
  skip_if_not_installed("survival")
  set.seed(650)
  ids <- paste0("S", 1:40)
  clusters <- setNames(rep(c(1, 2), each = 20), ids)
  time <- setNames(c(stats::rexp(20, rate = 0.1), stats::rexp(20, rate = 0.01)), ids)  ## cluster 1 much shorter survival
  event <- setNames(rep(1, 40), ids)
  out <- sfc_test_survival(clusters, time, event)
  expect_true(out$ok)
  expect_true(out$logrank_p < 0.05)
  expect_true(!is.null(out$hr))
  expect_true(out$hr$hr != 1)
})

test_that("sfc_test_categorical()/sfc_test_continuous() refuse with fewer than 6 matched, non-missing observations", {
  clusters <- setNames(c(1, 2), c("S1", "S2"))
  expect_false(sfc_test_categorical(clusters, setNames(c("A", "B"), c("S1", "S2")))$ok)
  expect_false(sfc_test_continuous(clusters, setNames(c(1, 2), c("S1", "S2")))$ok)
})

## ---- sfc_clinical_run() (BH-FDR across the user's own selected variables) ---

test_that("sfc_clinical_run() BH-corrects p-values across exactly the variables tested, in the same order", {
  set.seed(660)
  ids <- paste0("S", 1:20)
  meta <- data.frame(
    v1 = c(rep("A", 9), "B", rep("B", 9), "A"),   ## strongly associated
    v2 = sample(c("A", "B"), 20, replace = TRUE),  ## unrelated to cluster
    row.names = ids
  )
  clusters <- setNames(rep(c(1, 2), each = 10), ids)
  out <- sfc_clinical_run(clusters, meta, vars = c("v1", "v2"), kind = "categorical")
  raw_p <- vapply(out, `[[`, numeric(1), "p_value")
  expected_fdr <- stats::p.adjust(raw_p, method = "BH")
  actual_fdr <- vapply(out, `[[`, numeric(1), "p_fdr")
  expect_equal(unname(actual_fdr), unname(expected_fdr))
})

## ---- sfc_stability_run() (real repeated resampling + ARI) -------------------

test_that("sfc_stability_run() (real repeated SNF reruns + ARI) reports high stability for a strongly-separated planted structure", {
  skip_if_not_installed("SNFtool")
  set.seed(670)
  n <- 30
  ids <- paste0("S", seq_len(n))
  true_cluster <- rep(c(1, 2), each = n / 2)
  m1 <- matrix(rnorm(n * 15), n, 15, dimnames = list(ids, paste0("f", 1:15)))
  m1[true_cluster == 2, ] <- m1[true_cluster == 2, ] + 5
  ref <- setNames(true_cluster, ids)
  out <- sfc_stability_run(list(A = m1), ref, params = list(cluster_mode = "manual", n_clusters = 2), n_resamples = 8, subsample_frac = 0.8, seed = 1)
  expect_true(out$ok)
  expect_true(out$mean_ari > 0.5)
  expect_equal(out$verdict, sfc_stability_verdict(out$mean_ari))
})

test_that("sfc_stability_verdict() classifies at the documented 0.5/0.75 thresholds", {
  expect_equal(sfc_stability_verdict(0.8), "Stable")
  expect_equal(sfc_stability_verdict(0.6), "Moderately stable")
  expect_equal(sfc_stability_verdict(0.3), "Unstable")
  expect_equal(sfc_stability_verdict(NA_real_), "Not computable")
})

test_that("sfc_stability_run() refuses for a cohort too small to subsample meaningfully", {
  ## n=5 at subsample_frac=0.8 -> n_sub=round(4)=4, which is below the
  ## MI_MIN_MATCHED_SAMPLES(3)+2=5 floor - genuinely too small, unlike n=6
  ## (n_sub=5, which clears the floor and runs a real, if noisy, check).
  ids <- paste0("S", 1:5)
  ref <- setNames(c(1, 1, 2, 2, 2), ids)
  m1 <- matrix(rnorm(50), 5, 10, dimnames = list(ids, paste0("f", 1:10)))
  out <- sfc_stability_run(list(A = m1), ref, params = list(), n_resamples = 5, subsample_frac = 0.8)
  expect_false(out$ok)
})

## ---- sfc_sensitivity_run() (real re-clustering at parameter extremes) ------

test_that("sfc_sensitivity_run() reruns clustering at K/Alpha low/high (and T for >=2 blocks) and reports ARI vs. the reference", {
  skip_if_not_installed("SNFtool")
  set.seed(680)
  n <- 24
  ids <- paste0("S", seq_len(n))
  true_cluster <- rep(c(1, 2), each = n / 2)
  m1 <- matrix(rnorm(n * 12), n, 12, dimnames = list(ids, paste0("f", 1:12)))
  m1[true_cluster == 2, ] <- m1[true_cluster == 2, ] + 4
  m2 <- matrix(rnorm(n * 12), n, 12, dimnames = list(ids, paste0("g", 1:12)))
  m2[true_cluster == 2, ] <- m2[true_cluster == 2, ] + 4
  ref <- setNames(true_cluster, ids)

  base_params <- list(cluster_mode = "manual", n_clusters = 2)
  out <- sfc_sensitivity_run(list(A = m1, B = m2), base_params, ref, seed = 1)
  expect_true(is.data.frame(out$detail))
  expect_setequal(unique(out$detail$parameter), c("K", "Alpha", "T"))  ## T included since >=2 blocks
  expect_true(all(out$summary$sensitivity %in% c("Low sensitivity", "High sensitivity", "Not computable")))
})

test_that("sfc_sensitivity_run() skips T for a single omics block (no fusion step exists to vary)", {
  skip_if_not_installed("SNFtool")
  set.seed(690)
  n <- 20
  ids <- paste0("S", seq_len(n))
  m1 <- matrix(rnorm(n * 10), n, 10, dimnames = list(ids, paste0("f", 1:10)))
  ref <- setNames(rep(c(1, 2), each = n / 2), ids)
  out <- sfc_sensitivity_run(list(A = m1), list(cluster_mode = "manual", n_clusters = 2), ref, seed = 1)
  expect_setequal(unique(out$detail$parameter), c("K", "Alpha"))
})

## ---- sfc_feature_ranking() (real per-feature Kruskal-Wallis) -----------------

test_that("sfc_feature_ranking() ranks a real cluster-associated feature above an unrelated one, real BH-FDR", {
  set.seed(700)
  ids <- paste0("S", 1:20)
  clusters <- setNames(rep(c(1, 2), each = 10), ids)
  mat <- matrix(rnorm(20 * 5), 20, 5, dimnames = list(ids, paste0("f", 1:5)))
  mat[clusters == 2, 1] <- mat[clusters == 2, 1] + 6  ## f1 strongly associated with cluster
  out <- sfc_feature_ranking(mat, clusters, "TestBlock", top_n = 5)
  expect_true(out$ok)
  expect_equal(out$table$feature[1], "f1")  ## strongest association ranks first (lowest p)
  expect_equal(out$n_features_tested, 5L)
})

test_that("sfc_feature_ranking() refuses with too few matched samples or fewer than 2 clusters", {
  mat <- matrix(rnorm(10), 2, 5, dimnames = list(c("S1", "S2"), paste0("f", 1:5)))
  clusters <- setNames(c(1, 1), c("S1", "S2"))
  out <- sfc_feature_ranking(mat, clusters, "X")
  expect_false(out$ok)
})

## ---- sfc_summary_lines() -----------------------------------------------------

test_that("sfc_summary_lines() reports every real parameter used, including 'not applicable' for T on a single-omics run", {
  res <- list(
    params = list(blocks = "A", mode = "single_omics", n_samples = 20, n_clusters = 2, cluster_mode = "automatic",
                   k = 8, k_mode = "automatic", alpha = NA_real_, alpha_mode = "automatic", t = NA_integer_),
    clusters = setNames(rep(c(1, 2), each = 10), paste0("S", 1:20))
  )
  lines <- sfc_summary_lines(res)
  text <- paste(lines, collapse = "\n")
  expect_true(grepl("Single-Omics Clustering", text))
  expect_true(grepl("not applicable \\(single-omics\\)", text))
  expect_true(grepl("Cluster stability: Not computed", text))
})
