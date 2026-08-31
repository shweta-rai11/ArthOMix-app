## Module 2 (Methylomics) - Cell-Type Deconvolution's pure functions (scale
## detection/conversion, marker ranking/selection, overlap QC,
## reconstruction validation) plus a real EpiDISH::epidish() deconvolution
## run on a synthetic bulk mixture with a KNOWN true composition - the
## scientific-contract check this submodule is built around.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "parse_upload.R"))
source_from_app_root(file.path("R", "methylomics", "annotation.R"))
source_from_app_root(file.path("R", "methylomics", "celltype.R"))

## ---- Scale detection / conversion --------------------------------------------

test_that("methyl_ct_detect_scale() distinguishes beta, percent, and M-value scales", {
  expect_equal(methyl_ct_detect_scale(matrix(runif(100, 0, 1), 10, 10))$scale, "beta")
  expect_equal(methyl_ct_detect_scale(matrix(runif(100, 0, 100), 10, 10))$scale, "percent")
  expect_equal(methyl_ct_detect_scale(matrix(rnorm(100, 0, 4), 10, 10))$scale, "m")
})

test_that("methyl_ct_pct_to_beta()/methyl_ct_m_to_beta() are correct, inverse-consistent transforms", {
  expect_equal(methyl_ct_pct_to_beta(c(0, 50, 100)), c(0, 0.5, 1))
  beta <- c(0.2, 0.5, 0.8)
  m <- log2(beta / (1 - beta))
  expect_equal(methyl_ct_m_to_beta(m), beta, tolerance = 1e-10)
})

## ---- Marker ranking / selection ------------------------------------------------

ct_ref_fixture <- function() {
  ## 5 CpGs x 3 cell types - cg1/cg2 are strong hyper-markers for CellA,
  ## cg3 a strong hypo-marker for CellB, cg4/cg5 non-specific (similar across types).
  matrix(c(
    0.9, 0.1, 0.1,   ## cg1: CellA high
    0.85, 0.15, 0.1,  ## cg2: CellA high
    0.1, 0.9, 0.85,   ## cg3: CellB/CellC high, CellA low -> hypo marker for CellA... let's design cleanly below
    0.5, 0.5, 0.5,    ## cg4: non-specific
    0.5, 0.52, 0.48   ## cg5: non-specific
  ), nrow = 5, byrow = TRUE, dimnames = list(paste0("cg", 1:5), c("CellA", "CellB", "CellC")))
}

test_that("methyl_ct_marker_rank() correctly identifies each CpG's strongest-marked cell type and direction", {
  ref <- ct_ref_fixture()
  rank_df <- methyl_ct_marker_rank(ref)
  expect_equal(rank_df$cell_type[rank_df$cpg == "cg1"], "CellA")
  expect_equal(rank_df$direction[rank_df$cpg == "cg1"], "hyper")
  ## cg3: CellA is far lower than CellB/CellC -> a hypo marker for CellA.
  expect_equal(rank_df$cell_type[rank_df$cpg == "cg3"], "CellA")
  expect_equal(rank_df$direction[rank_df$cpg == "cg3"], "hypo")
  ## Non-specific CpGs should have much smaller effect sizes than the true markers.
  expect_true(rank_df$effect[rank_df$cpg == "cg1"] > rank_df$effect[rank_df$cpg == "cg4"])
})

test_that("methyl_ct_marker_rank() errors clearly with fewer than 2 cell types", {
  ref <- matrix(0.5, 5, 1, dimnames = list(paste0("cg", 1:5), "OnlyOne"))
  expect_error(methyl_ct_marker_rank(ref), "at least 2 cell types")
})

test_that("methyl_ct_select_markers() filters by minimum effect size and direction", {
  ref <- ct_ref_fixture()
  rank_df <- methyl_ct_marker_rank(ref)
  strong_only <- methyl_ct_select_markers(rank_df, effect_min = 0.3)
  expect_true(all(strong_only$effect >= 0.3))
  hyper_only <- methyl_ct_select_markers(rank_df, direction = "hyper")
  expect_true(all(hyper_only$direction == "hyper"))
})

test_that("methyl_ct_top_n_balanced() caps the total and balances picks roughly evenly across cell types", {
  set.seed(310)
  df <- data.frame(cpg = paste0("cg", 1:30), cell_type = rep(c("A", "B", "C"), each = 10),
                     effect = runif(30, 0.1, 1))
  out <- methyl_ct_top_n_balanced(df, sort_col = "effect", top_n = 9)
  expect_equal(nrow(out), 9L)
  tab <- table(out$cell_type)
  expect_true(all(tab >= 2))  ## roughly balanced, not all from one type
})

## ---- Overlap QC ---------------------------------------------------------------

test_that("methyl_ct_overlap_qc() reports matched/missing marker CpGs and percentage", {
  out <- methyl_ct_overlap_qc(marker_ids = c("cg1", "cg2", "cg3"), working_ids = c("cg1", "cg2", "cg9"))
  expect_equal(out$n_matched, 2L)
  expect_equal(out$n_missing, 1L)
  expect_equal(out$pct_matched, 200 / 3, tolerance = 1e-6)
})

test_that("methyl_ct_overlap_by_type() reports per-cell-type match percentages", {
  marker_df <- data.frame(cpg = c("cg1", "cg2", "cg3", "cg4"), cell_type = c("A", "A", "B", "B"))
  out <- methyl_ct_overlap_by_type(marker_df, working_ids = c("cg1", "cg3"))
  expect_equal(out$pct_matched[out$cell_type == "A"], 50)
  expect_equal(out$pct_matched[out$cell_type == "B"], 50)
})

## ---- Real EpiDISH deconvolution on a synthetic known-composition mixture -----

test_that("methyl_ct_run_epidish() (real EpiDISH::epidish, RPC) approximately recovers a known true cell-type mixture", {
  set.seed(311)
  n_markers <- 200
  ## Two well-separated cell-type reference profiles, with realistic per-CpG
  ## noise - a perfectly piecewise-constant reference (no within-cluster
  ## variability at all) makes RPC's internal MASS::rlm() fit singular
  ## ("'x' is singular: singular fits are not implemented in 'rlm'"),
  ## confirmed directly against EpiDISH; this is a property of the fixture,
  ## not a real bug, so realistic noise is added instead of switching methods.
  ref_mat <- matrix(NA_real_, n_markers, 2, dimnames = list(paste0("cg", 1:n_markers), c("CellA", "CellB")))
  ref_mat[, "CellA"] <- pmin(pmax(rnorm(n_markers, mean = c(rep(0.9, n_markers / 2), rep(0.1, n_markers / 2)), sd = 0.02), 0), 1)
  ref_mat[, "CellB"] <- pmin(pmax(rnorm(n_markers, mean = c(rep(0.1, n_markers / 2), rep(0.9, n_markers / 2)), sd = 0.02), 0), 1)

  ## One synthetic bulk sample: a known 70% CellA / 30% CellB mixture, plus small noise.
  true_frac <- c(CellA = 0.7, CellB = 0.3)
  bulk <- ref_mat %*% matrix(true_frac, ncol = 1)
  bulk <- bulk + rnorm(n_markers, sd = 0.01)
  bulk <- pmin(pmax(bulk, 0), 1)
  colnames(bulk) <- "S1"

  res <- methyl_ct_run_epidish(bulk, ref_mat, method = "RPC")
  expect_true(res$ok)
  est <- res$fractions["S1", ]
  ## Absolute-difference bound (not expect_equal(..., tolerance=)): the intent
  ## here is "estimated fraction within 0.05 of the true fraction", but
  ## expect_equal()'s `tolerance` is interpreted as a *relative* difference
  ## under testthat's all.equal()-based edition 2 and differently again under
  ## edition 3's waldo::compare(); for a value as small as 0.3 that distinction
  ## flips a 0.025 absolute (and passing) deviation into a test failure under
  ## edition 3. expect_lt() on the raw absolute difference is edition-agnostic
  ## and states the actual intended bound directly (verified against real
  ## testthat 3.3.2: this file fails under Config/testthat/edition: 3 with the
  ## old expect_equal(tolerance=) form, and passes under both editions with
  ## this form).
  expect_lt(abs(unname(est["CellA"]) - 0.7), 0.05)
  expect_lt(abs(unname(est["CellB"]) - 0.3), 0.05)
  expect_equal(sum(est), 1, tolerance = 1e-6)  ## EpiDISH's own sum-to-one constraint
})

test_that("methyl_ct_run_epidish() fails soft with too few overlapping marker CpGs", {
  bulk <- matrix(0.5, 3, 1, dimnames = list(paste0("cg", 1:3), "S1"))
  ref <- matrix(0.5, 3, 2, dimnames = list(paste0("zz", 1:3), c("A", "B")))
  res <- methyl_ct_run_epidish(bulk, ref, method = "RPC")
  expect_false(res$ok)
  expect_true(grepl("too few for deconvolution", res$reason))
})

## ---- Reconstruction validation --------------------------------------------------

test_that("methyl_ct_reconstruct()/methyl_ct_validation_metrics() report near-perfect fit when fractions exactly match the true mixture", {
  ref_mat <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, dimnames = list(c("cg1", "cg2"), c("CellA", "CellB")))
  fractions <- matrix(c(0.7, 0.3), 1, 2, dimnames = list("S1", c("CellA", "CellB")))
  observed <- methyl_ct_reconstruct(ref_mat, fractions)  ## build "observed" as the exact reconstruction
  rec <- methyl_ct_reconstruct(ref_mat, fractions)
  colnames(observed) <- "S1"; colnames(rec) <- "S1"

  metrics <- methyl_ct_validation_metrics(observed, rec)
  expect_true(metrics$ok)
  expect_equal(metrics$overall$rmse, 0, tolerance = 1e-10)
  expect_equal(metrics$overall$r2, 1, tolerance = 1e-8)
})

test_that("methyl_ct_validation_metrics() fails soft with fewer than 2 overlapping CpGs", {
  observed <- matrix(0.5, 1, 1, dimnames = list("cg1", "S1"))
  reconstructed <- matrix(0.5, 1, 1, dimnames = list("cg1", "S1"))
  out <- methyl_ct_validation_metrics(observed, reconstructed)
  expect_false(out$ok)
})
