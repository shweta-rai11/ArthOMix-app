## Regression coverage for a gap found in the transcriptomics audit
## (2026-09-03): every existing Preprocessing/Batch-Correction test sets
## skip_combat = TRUE, so ComBat / limma::removeBatchEffect / SVA are never
## actually executed by any test - a real regression in the batch-correction
## engines themselves would go undetected. These tests drive the real
## run_btn -> result() path (skip_combat = FALSE) on a synthetic two-batch
## fixture with a KNOWN injected batch offset and a KNOWN
## injected, batch-orthogonal group signal, and check: (a) dimensions are
## preserved, (b) no NaN/Inf is introduced, (c) the injected batch effect is
## measurably reduced, and (d) the real biological (group) signal survives
## correction - distinguishing genuine correction from mere data destruction.
##
## Covers: ComBat (default) and limma::removeBatchEffect. SVA and
## ComBat-seq/TMM are not covered here (time-boxed) - see the module's own
## existing tests for TMM-specific validation gates, which are unaffected by
## this gap.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "functions", "expression_type.R"))
source_from_app_root(file.path("R", "transcriptomics", "01_Data", "mod_dataset.R"))
source_from_app_root(file.path("R", "transcriptomics", "03_Preprocessing_Batch_Correction", "mod_preprocessing_explore.R"))
source_from_app_root(file.path("R", "transcriptomics", "03_Preprocessing_Batch_Correction", "mod_preprocessing.R"))

## Two "batches" (a single synthetic dataset with its own "batch" column),
## each with an orthogonal 50/50 HC/RA group split. Every gene gets a fixed
## additive offset in batch B (the injected, known batch effect); a subset of
## "signal" genes additionally get a fixed additive offset in RA vs HC, in
## BOTH batches (the injected, known, batch-orthogonal biological signal).
##
## Built as one combined expr/meta pair (with meta$batch already set) rather
## than two separately-uploaded sources, so it can be fed in through the
## "Currently Loaded Dataset" preloaded-cohort path: merge_inputs() ->
## merged() passes a single preloaded source's meta straight through
## (mod_preprocessing.R's merged(), the length(lst) == 1 branch), preserving
## an already-populated "batch" column untouched.
pp_make_batch_signal_fixture <- function(seed_a = 101, seed_b = 102, batch_offset = 4,
                                          n_genes = 80, n_signal = 10, n_per_group = 8,
                                          group_effect = 3) {
  make_one <- function(sample_prefix, seed, offset) {
    set.seed(seed)
    genes <- paste0("GENE", seq_len(n_genes))
    n <- n_per_group * 2
    samples <- paste0(sample_prefix, seq_len(n))
    group <- rep(c("HC", "RA"), each = n_per_group)

    m <- matrix(rnorm(n_genes * n, mean = 8, sd = 0.4), n_genes, n, dimnames = list(genes, samples))
    m <- m + offset
    signal_genes <- genes[seq_len(n_signal)]
    m[signal_genes, group == "RA"] <- m[signal_genes, group == "RA"] + group_effect
    list(expr = m, samples = samples, group = group, signal_genes = signal_genes,
         non_signal_genes = setdiff(genes, signal_genes))
  }
  a <- make_one("A", seed_a, 0)
  b <- make_one("B", seed_b, batch_offset)
  expr <- cbind(a$expr, b$expr)
  meta <- data.frame(
    sample = c(a$samples, b$samples),
    group = c(a$group, b$group),
    batch = rep(c("A", "B"), c(length(a$samples), length(b$samples))),
    stringsAsFactors = FALSE
  )
  list(expr = expr, meta = meta, fx1 = a, fx2 = b)
}

## Loads the two-batch fixture as the "Currently Loaded Dataset" preloaded
## source and merges it through the module's real preloaded+merge UI path,
## then runs batch correction with skip_combat = FALSE and the given
## correction_method, returning result() plus the fixture's ground truth.
run_pp_batch_correction <- function(correction_method = "combat") {
  fx <- pp_make_batch_signal_fixture(batch_offset = 4)
  dataset <- shiny::reactiveValues(expr = fx$expr, meta = fx$meta,
                                    source = "Synthetic two-batch fixture", source_type = "preloaded")

  out <- NULL
  shiny::testServer(mod_preprocessing_server, args = list(id = "pp", dataset = dataset, results = shiny::reactiveValues()), {
    session$setInputs(preloaded_selected = "__current__", preloaded_log2 = "skip")
    session$setInputs(preloaded_run = 1)
    session$setInputs(merge_mode = "own")
    session$setInputs(merge_btn = 1)
    m <- merged()
    ## Sanity-check the fixture itself merged as expected: a single preloaded
    ## source whose meta already carried a real "batch" column, preserved as-is.
    testthat::expect_true("batch" %in% colnames(m$meta))
    testthat::expect_length(unique(m$meta$batch), 2L)

    session$setInputs(color_by = "group", batch_col = "batch", norm_method = "skip",
                       skip_combat = FALSE, correction_method = correction_method,
                       protect_cols = "group", mad_k = 3, min_pct = 0, variance_pct = 0)
    session$setInputs(run_btn = 1)
    out <<- result()
  })
  list(res = out, fx1 = fx$fx1, fx2 = fx$fx2)
}

## Mean expression, per fixture-defined gene set, restricted to the given
## batch/group subset of a corrected-or-uncorrected expression matrix.
bc_group_mean <- function(expr, fx1, fx2, genes, which_batch, which_group) {
  samples <- if (which_batch == "A") fx1$samples[fx1$group == which_group] else fx2$samples[fx2$group == which_group]
  mean(expr[genes, samples, drop = FALSE])
}

test_that("ComBat batch correction preserves matrix dimensions and introduces no NaN/Inf", {
  out <- run_pp_batch_correction("combat")
  res <- out$res
  expect_identical(dim(res$expr_combat), dim(res$expr_qnorm))
  expect_identical(sort(rownames(res$expr_combat)), sort(rownames(res$expr_qnorm)))
  expect_identical(sort(colnames(res$expr_combat)), sort(colnames(res$expr_qnorm)))
  expect_true(all(is.finite(as.matrix(res$expr_combat))))
  expect_null(res$combat_fallback_note)
  expect_false(isTRUE(res$skip_combat))
})

test_that("ComBat measurably reduces the injected between-batch offset on non-signal genes", {
  out <- run_pp_batch_correction("combat")
  res <- out$res; fx1 <- out$fx1; fx2 <- out$fx2
  non_signal <- fx1$non_signal_genes

  before_gap <- abs(mean(res$expr_qnorm[non_signal, fx2$samples]) - mean(res$expr_qnorm[non_signal, fx1$samples]))
  after_gap  <- abs(mean(res$expr_combat[non_signal, fx2$samples]) - mean(res$expr_combat[non_signal, fx1$samples]))

  ## A batch_offset of 4 was injected; ComBat should shrink the observed gap
  ## by at least an order of magnitude (not merely "some" reduction, and not
  ## a no-op that would leave the ~4-unit gap intact).
  expect_gt(before_gap, 3)
  expect_lt(after_gap, before_gap / 5)
  expect_lt(after_gap, 0.5)
})

test_that("ComBat preserves the real, batch-orthogonal biological (group) signal on signal genes", {
  out <- run_pp_batch_correction("combat")
  res <- out$res; fx1 <- out$fx1; fx2 <- out$fx2
  sig <- fx1$signal_genes

  ## Group effect (RA - HC) computed within each batch separately, before and
  ## after correction - a group_effect of 3 was injected identically in both
  ## batches (orthogonal to the batch offset), so correction should not
  ## remove it, whichever batch it is measured in.
  before_a <- bc_group_mean(res$expr_qnorm, fx1, fx2, sig, "A", "RA") - bc_group_mean(res$expr_qnorm, fx1, fx2, sig, "A", "HC")
  after_a  <- bc_group_mean(res$expr_combat, fx1, fx2, sig, "A", "RA") - bc_group_mean(res$expr_combat, fx1, fx2, sig, "A", "HC")
  before_b <- bc_group_mean(res$expr_qnorm, fx1, fx2, sig, "B", "RA") - bc_group_mean(res$expr_qnorm, fx1, fx2, sig, "B", "HC")
  after_b  <- bc_group_mean(res$expr_combat, fx1, fx2, sig, "B", "RA") - bc_group_mean(res$expr_combat, fx1, fx2, sig, "B", "HC")

  expect_gt(before_a, 2); expect_gt(before_b, 2)
  ## The corrected group effect should stay within the same ballpark as the
  ## injected 3-unit effect in both batches - correction, not destruction.
  expect_gt(after_a, 1.5); expect_lt(after_a, 4.5)
  expect_gt(after_b, 1.5); expect_lt(after_b, 4.5)
})

test_that("limma::removeBatchEffect also preserves dimensions, reduces batch offset, and preserves group signal", {
  out <- run_pp_batch_correction("limma")
  res <- out$res; fx1 <- out$fx1; fx2 <- out$fx2
  expect_identical(dim(res$expr_combat), dim(res$expr_qnorm))
  expect_true(all(is.finite(as.matrix(res$expr_combat))))

  non_signal <- fx1$non_signal_genes
  before_gap <- abs(mean(res$expr_qnorm[non_signal, fx2$samples]) - mean(res$expr_qnorm[non_signal, fx1$samples]))
  after_gap  <- abs(mean(res$expr_combat[non_signal, fx2$samples]) - mean(res$expr_combat[non_signal, fx1$samples]))
  expect_gt(before_gap, 3)
  expect_lt(after_gap, 0.5)

  sig <- fx1$signal_genes
  after_a <- bc_group_mean(res$expr_combat, fx1, fx2, sig, "A", "RA") - bc_group_mean(res$expr_combat, fx1, fx2, sig, "A", "HC")
  expect_gt(after_a, 1.5)
})

test_that("disabling batch correction (skip_combat = TRUE) leaves the injected batch offset fully intact, as a negative control", {
  fx <- pp_make_batch_signal_fixture(batch_offset = 4)
  fx1 <- fx$fx1; fx2 <- fx$fx2
  dataset <- shiny::reactiveValues(expr = fx$expr, meta = fx$meta,
                                    source = "Synthetic two-batch fixture", source_type = "preloaded")

  shiny::testServer(mod_preprocessing_server, args = list(id = "pp", dataset = dataset, results = shiny::reactiveValues()), {
    session$setInputs(preloaded_selected = "__current__", preloaded_log2 = "skip")
    session$setInputs(preloaded_run = 1)
    session$setInputs(merge_mode = "own")
    session$setInputs(merge_btn = 1)

    session$setInputs(color_by = "group", batch_col = "batch", norm_method = "skip",
                       skip_combat = TRUE, mad_k = 3, min_pct = 0, variance_pct = 0)
    session$setInputs(run_btn = 1)
    res <- result()
    expect_true(isTRUE(res$skip_combat))
    expect_identical(res$expr_combat, res$expr_qnorm)
    non_signal <- fx1$non_signal_genes
    gap <- abs(mean(res$expr_combat[non_signal, fx2$samples]) - mean(res$expr_combat[non_signal, fx1$samples]))
    expect_gt(gap, 3)
  })
})
