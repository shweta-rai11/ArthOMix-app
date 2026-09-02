## Module 2 (Methylomics) - Normalization tab's beta-matrix-level pure
## functions (normalization.R): the plain quantile method (works on any
## input, no raw IDAT needed), the Normalization-tab-specific probe/sample
## filters, and the automatic diagnostics/status/recommendation logic.
## Raw-IDAT-only methods (Noob/Funnorm/SWAN/Dasen/stratified-quantile) and
## BMIQ/PBC's real numerical fits are not exercised here - no committed
## sample IDAT fixture exists in this project, and BMIQ's mixture-model fit
## needs realistic bimodal beta distributions to converge meaningfully
## rather than a synthetic uniform fixture (documented gap, not silently
## skipped).

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "functions", "annotation.R"))
source_from_app_root(file.path("R", "methylomics", "functions", "normalization.R"))

## ---- methyl_norm_quantile() (real limma::normalizeQuantiles) --------------

test_that("methyl_norm_quantile() makes every sample's quantiles identical afterward", {
  set.seed(220)
  mat <- matrix(c(rexp(20, rate = 2), rexp(20, rate = 5)), 20, 2)  ## deliberately different distributions per sample
  out <- methyl_norm_quantile(mat)
  expect_true(out$ok)
  q1 <- sort(out$beta[, 1]); q2 <- sort(out$beta[, 2])
  expect_equal(q1, q2, tolerance = 1e-8)
})

## ---- methyl_design_vector() (real 450K annotation) ------------------------

test_that("methyl_design_vector() correctly labels Type I/II probes from the real 450K manifest and drops unmatched IDs", {
  anno_result <- methyl_get_annotation("450K")
  real_probes <- rownames(anno_result$anno)[1:20]
  probe_ids <- c(real_probes, "cg_not_a_real_probe")
  out <- methyl_design_vector(probe_ids, anno_result)
  expect_true(out$ok)
  expect_equal(out$n_dropped, 1L)
  expect_true(all(out$design_v[real_probes] %in% c(1, 2)))
  expect_true(is.na(out$design_v["cg_not_a_real_probe"]))
})

test_that("methyl_design_vector() fails soft when the annotation itself is unavailable", {
  out <- methyl_design_vector(c("cg001"), list(ok = FALSE, reason = "no manifest"))
  expect_false(out$ok)
})

## ---- Normalization-tab probe/sample filters --------------------------------

test_that("methyl_filter_chromosome() removes probes on the excluded chromosome(s) using real annotation", {
  anno_result <- methyl_get_annotation("450K")
  a <- anno_result$anno
  chrY_probes <- rownames(a)[which(a$chr == "chrY")][1:3]
  chr1_probes <- rownames(a)[which(a$chr == "chr1")][1:3]
  mat <- matrix(0.5, 6, 2, dimnames = list(c(chrY_probes, chr1_probes), c("S1", "S2")))
  out <- methyl_filter_chromosome(mat, anno_result, exclude_chr = "chrY")
  expect_equal(out$keep, c(rep(FALSE, 3), rep(TRUE, 3)))
})

test_that("methyl_filter_chromosome() is a no-op when no exclusion is selected", {
  mat <- matrix(0.5, 4, 2)
  out <- methyl_filter_chromosome(mat, list(ok = TRUE, anno = data.frame()), exclude_chr = character(0))
  expect_true(all(out$keep))
})

test_that("methyl_filter_island_relation()/methyl_filter_gene_region() filter by the requested categories and are no-ops without them", {
  a <- data.frame(island_relation = c("Island", "Island", "OpenSea"), gene_region = c("TSS200", "Body", "Body"),
                    row.names = c("cg1", "cg2", "cg3"))
  mat <- matrix(0.5, 3, 2, dimnames = list(rownames(a), c("S1", "S2")))
  anno_result <- list(ok = TRUE, anno = a)

  out_isl <- methyl_filter_island_relation(mat, anno_result, keep_categories = "Island")
  expect_equal(out_isl$keep, c(TRUE, TRUE, FALSE))

  out_reg <- methyl_filter_gene_region(mat, anno_result, keep_regions = "Body")
  expect_equal(out_reg$keep, c(FALSE, TRUE, TRUE))

  out_noop <- methyl_filter_island_relation(mat, anno_result, keep_categories = character(0))
  expect_true(all(out_noop$keep))
})

test_that("methyl_filter_samples_missingness() drops samples over the missingness threshold and no-ops when threshold is NULL", {
  mat <- matrix(0.5, 10, 3)
  mat[1:9, 1] <- NA  ## sample 1: 90% missing
  out <- methyl_filter_samples_missingness(mat, max_na_frac = 0.5)
  expect_equal(out$keep, c(FALSE, TRUE, TRUE))
  expect_true(all(methyl_filter_samples_missingness(mat, NULL)$keep))
})

## ---- Diagnostics / status / recommendation ---------------------------------

test_that("methyl_norm_diagnostics() reports correct summary stats for a crafted matrix", {
  mat <- matrix(c(0.1, 0.5, NA, Inf, 0.9, 0.3), nrow = 2)
  dataset <- list(input_scale = "beta")
  d <- methyl_norm_diagnostics(mat, dataset)
  expect_equal(d$n_probes, 2L)
  expect_equal(d$n_samples, 3L)
  expect_equal(d$n_missing, 1L)
  expect_equal(d$n_infinite, 1L)
  expect_equal(d$value_min, 0.1)
})

test_that("methyl_norm_status() flags raw IDAT data as needing normalization regardless of anything else", {
  status <- methyl_norm_status(matrix(0.5, 2, 2), dataset = list(rg_set = "not really null"), anno_result = NULL)
  expect_equal(status$status, "raw")
})

test_that("methyl_norm_status() reports 'unknown' for M-value-scale data (no bias signature available)", {
  status <- methyl_norm_status(matrix(0.5, 2, 2), dataset = list(rg_set = NULL, input_scale = "m"), anno_result = NULL)
  expect_equal(status$status, "unknown")
})

test_that("methyl_norm_recommendation() gives non-binding, context-specific text for each dataset/status combination", {
  raw_rec <- methyl_norm_recommendation(list(rg_set = "x"), status = list(status = "unknown"), available_methods = character(0))
  expect_true(grepl("raw Illumina methylation-array", raw_rec))

  no_bias_rec <- methyl_norm_recommendation(list(rg_set = NULL), status = list(status = "no_bias_detected"), available_methods = character(0))
  expect_true(grepl("no evidence of uncorrected", no_bias_rec))

  bias_rec <- methyl_norm_recommendation(list(rg_set = NULL), status = list(status = "bias_detected"), available_methods = "bmiq")
  expect_true(grepl("evidence of uncorrected Type I/II probe-design bias", bias_rec))

  fallback_rec <- methyl_norm_recommendation(list(rg_set = NULL), status = list(status = "unknown"), available_methods = character(0))
  expect_true(grepl("plain quantile normalization is the only method", fallback_rec))
})
