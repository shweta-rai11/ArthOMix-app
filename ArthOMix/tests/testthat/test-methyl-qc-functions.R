## Module 2 (Methylomics) - QC's pure probe/sample-level filter and scoring
## functions (qc.R). Plotting/report-generation helpers (methyl_plot_*,
## methyl_qc_report_*) are intentionally not unit-tested here (presentation,

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "functions", "qc.R"))

qc_fixture_mat <- function(n_probes = 20, n_samples = 8, seed = 210) {
  set.seed(seed)
  matrix(runif(n_probes * n_samples, 0, 1), n_probes, n_samples,
          dimnames = list(paste0("cg", 10000000 + 1:n_probes), paste0("S", 1:n_samples)))
}

test_that("methyl_filter_missing() keeps only probes at/under the missingness threshold", {
  mat <- qc_fixture_mat()
  mat[1, 1:6] <- NA
  out <- methyl_filter_missing(mat, max_na_frac = 0.5)
  expect_false(out$keep[1])
  expect_true(all(out$keep[-1]))
})

test_that("methyl_filter_variance()/methyl_filter_sd() agree on which probes are dropped (SD = sqrt(variance))", {
  mat <- qc_fixture_mat()
  mat[1, ] <- 0.5
  v_out <- methyl_filter_variance(mat, min_variance = 0.001)
  sd_out <- methyl_filter_sd(mat, min_sd = sqrt(0.001))
  expect_false(v_out$keep[1])
  expect_identical(v_out$keep, sd_out$keep)
})

test_that("methyl_filter_mean_range() drops probes whose mean beta falls outside the given range", {
  mat <- qc_fixture_mat()
  mat[1, ] <- 0.01
  out <- methyl_filter_mean_range(mat, lo = 0.1, hi = 0.9)
  expect_false(out$keep[1])
})

test_that("methyl_filter_detection_p()/methyl_filter_beadcount() are no-ops (keep=all-TRUE) when raw IDAT data isn't available", {
  mat <- qc_fixture_mat()
  dp_out <- methyl_filter_detection_p(mat, detp = NULL)
  bc_out <- methyl_filter_beadcount(mat, beadcount = NULL)
  expect_true(all(dp_out$keep))
  expect_true(all(bc_out$keep))
  expect_true(grepl("raw IDAT input", dp_out$note))
})

test_that("methyl_filter_detection_p() drops a probe failing detection p in any one sample", {
  mat <- qc_fixture_mat(n_probes = 5, n_samples = 4)
  detp <- matrix(0.001, 5, 4, dimnames = dimnames(mat))
  detp[2, 3] <- 0.5
  out <- methyl_filter_detection_p(mat, detp, threshold = 0.01)
  expect_equal(unname(out$keep), c(TRUE, FALSE, TRUE, TRUE, TRUE))
})

test_that("methyl_filter_beadcount() drops a probe with low bead count in any one sample", {
  mat <- qc_fixture_mat(n_probes = 5, n_samples = 4)
  beadcount <- matrix(10, 5, 4, dimnames = dimnames(mat))
  beadcount[3, 1] <- 1
  out <- methyl_filter_beadcount(mat, beadcount, threshold = 3)
  expect_equal(unname(out$keep), c(TRUE, TRUE, FALSE, TRUE, TRUE))
})

test_that("methyl_filter_cross_reactive()/methyl_filter_maf() are no-ops without an uploaded list/table, and filter correctly when provided", {
  mat <- qc_fixture_mat(n_probes = 5, n_samples = 4)
  expect_true(all(methyl_filter_cross_reactive(mat, NULL)$keep))
  out <- methyl_filter_cross_reactive(mat, exclusion_ids = rownames(mat)[2])
  expect_equal(which(!out$keep), 2)

  expect_true(all(methyl_filter_maf(mat, NULL)$keep))
  maf_table <- stats::setNames(c(0.1, 0.01), rownames(mat)[1:2])
  out2 <- methyl_filter_maf(mat, maf_table, max_maf = 0.05)
  expect_equal(which(!out2$keep), 1)
})

test_that("methyl_sample_call_rate() computes 1 - fraction missing per sample", {
  mat <- qc_fixture_mat(n_probes = 10, n_samples = 3)
  mat[1:5, 1] <- NA
  cr <- methyl_sample_call_rate(mat)
  expect_equal(unname(cr[1]), 0.5)
  expect_equal(unname(cr[2]), 1)
})

test_that("methyl_sample_failed_probe_pct() computes per-sample failed-probe percentage and flags no-overlap gracefully", {
  mat <- qc_fixture_mat(n_probes = 5, n_samples = 3)
  detp <- matrix(0.001, 5, 3, dimnames = dimnames(mat))
  detp[1:3, 2] <- 0.5
  out <- methyl_sample_failed_probe_pct(mat, detp)
  expect_true(out$ok)
  expect_equal(unname(out$pct["S2"]), 60)

  no_overlap <- methyl_sample_failed_probe_pct(mat, matrix(0.001, 5, 3, dimnames = list(paste0("zz", 1:5), paste0("Z", 1:3))))
  expect_false(no_overlap$ok)
})

test_that("methyl_cluster_sex() separates two well-separated groups and assigns higher-value cluster to M by default", {
  set.seed(211)
  y <- c(rnorm(10, mean = 0.2, sd = 0.02), rnorm(10, mean = 0.8, sd = 0.02))
  out <- methyl_cluster_sex(y)
  expect_true(out$direction_assumed)
  expect_true(all(out$sex[11:20] == "M"))
  expect_true(all(out$sex[1:10] == "F"))
})

test_that("methyl_cluster_sex() resolves cluster-to-sex direction by majority concordance when reported sex is given", {
  set.seed(212)
  y <- c(rnorm(10, mean = 0.2, sd = 0.02), rnorm(10, mean = 0.8, sd = 0.02))
  reported <- c(rep("M", 10), rep("F", 10))
  out <- methyl_cluster_sex(y, reported = reported)
  expect_false(out$direction_assumed)
  expect_true(all(out$sex[11:20] == "F"))
})

test_that("methyl_sex_check() (no rg_set/minfi) falls back to chrX/chrY beta clustering and flags reported-sex mismatches", {
  set.seed(213)
  n_probes <- 30
  probe_ids <- paste0("cg", 1:n_probes)
  chr <- c(rep("chrX", 15), rep("chrY", 10), rep("chr1", 5))
  anno <- data.frame(chr = chr, row.names = probe_ids)
  anno_result <- list(ok = TRUE, anno = anno)

  n_samples <- 10
  mat <- matrix(runif(n_probes * n_samples, 0.3, 0.5), n_probes, n_samples,
                 dimnames = list(probe_ids, paste0("S", 1:n_samples)))
  y_idx <- which(chr == "chrY")
  mat[y_idx, 1:5] <- runif(length(y_idx) * 5, 0.7, 0.9)

  reported_sex <- stats::setNames(c(rep("M", 5), rep("F", 4), "M"), colnames(mat))
  res <- methyl_sex_check(mat, anno_result, rg_set = NULL, reported_sex = reported_sex)
  expect_true(res$ok)
  expect_true(all(res$detail$predicted_sex[1:5] == "M"))
  expect_true(all(res$detail$predicted_sex[6:9] == "F"))
  expect_true(res$detail$sex_mismatch[res$detail$sample == "S10"])
})

test_that("methyl_sex_check() reports failure gracefully when too few chrX/chrY probes are annotated", {
  mat <- qc_fixture_mat(n_probes = 5, n_samples = 4)
  anno_result <- list(ok = TRUE, anno = data.frame(chr = rep("chr1", 5), row.names = rownames(mat)))
  res <- methyl_sex_check(mat, anno_result)
  expect_false(res$ok)
  expect_true(grepl("Too few annotated chrX/chrY probes", res$reason))
})

test_that("methyl_qc_subgroup_filter() restricts to one stratum and flags low_n below min_n", {
  mat <- qc_fixture_mat(n_probes = 5, n_samples = 6)
  sheet <- data.frame(sample = colnames(mat), group = c("A", "A", "B", "B", "B", "B"), stringsAsFactors = FALSE)
  out_a <- methyl_qc_subgroup_filter(mat, sheet, "group", "A", min_n = 3)
  expect_equal(ncol(out_a$mat), 2L)
  expect_true(out_a$low_n)

  out_b <- methyl_qc_subgroup_filter(mat, sheet, "group", "B", min_n = 3)
  expect_equal(ncol(out_b$mat), 4L)
  expect_false(out_b$low_n)

  out_all <- methyl_qc_subgroup_filter(mat, sheet, "group", "__all__", min_n = 3)
  expect_equal(ncol(out_all$mat), 6L)
})

test_that("methyl_apply_manual_exclude() removes the given samples and updates low_n/label", {
  mat <- qc_fixture_mat(n_probes = 5, n_samples = 5)
  sheet <- data.frame(sample = colnames(mat), group = "A", stringsAsFactors = FALSE)
  subgroup <- methyl_qc_subgroup_filter(mat, sheet, "group", "__all__", min_n = 3)
  out <- methyl_apply_manual_exclude(subgroup, excluded_ids = colnames(mat)[1:3])
  expect_equal(ncol(out$mat), 2L)
  expect_true(out$low_n)
  expect_true(grepl("manually excluded", out$label))
})

test_that("methyl_probe_retention_cascade() reports cumulative AND-combined retention across sequential filters", {
  n <- 100
  f1 <- list(keep = c(rep(TRUE, 80), rep(FALSE, 20)))
  f2 <- list(keep = c(rep(FALSE, 10), rep(TRUE, 90)))
  cascade <- methyl_probe_retention_cascade(n, list(missingness = f1, variance = f2))
  expect_equal(cascade$retained, c(100, 80, 70))
})

test_that("methyl_beta_to_mvalue() matches the standard logit2 transform and clips away from 0/1", {
  beta <- c(0.5, 0, 1, 0.25)
  m <- methyl_beta_to_mvalue(beta)
  expect_equal(m[1], 0)
  expect_true(is.finite(m[2]) && is.finite(m[3]))
  expect_equal(m[4], log2(0.25 / 0.75))
})
