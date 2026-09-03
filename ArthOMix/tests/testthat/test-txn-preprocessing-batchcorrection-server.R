## Coverage gap found in the 2026-09 test audit: every existing Preprocessing
## test that reaches "Run normalisation and batch correction" sets
## skip_combat=TRUE / norm_method="skip", so ComBat, limma::removeBatchEffect
## and SVA (mod_preprocessing.R's run_combat()/run_limma()/run_sva(), inline
## closures inside the result <- eventReactive(input$run_btn, ...) body around
## line 1350-1600) are never actually executed by any test. These closures
## aren't factored out as standalone functions, so they're driven here via
## shiny::testServer() against mod_preprocessing_server(), following the
## fixture/session-mocking idioms in test-txn-preprocessing-multi-upload-server.R
## and test-txn-preprocessing-dedup-server.R (single already-loaded dataset,
## fed through the "__current__" preloaded-read path, then "own" merge).
##
## fx_batch_signal_data() (helper-fixtures.R) builds a balanced batch x group
## design with a KNOWN injected batch offset on every gene and a KNOWN
## injected group effect on a subset of "signal" genes, orthogonal to batch.
## That lets every test here assert real numerical claims instead of smoke
## checks: dimensions preserved, no NaN/Inf introduced, the batch effect
## shrinks by a large factor post-correction, and the biological signal
## survives (doesn't just get destroyed along with the batch effect).

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "functions", "expression_type.R"))
source_from_app_root(file.path("R", "transcriptomics", "01_Data", "mod_dataset.R"))
source_from_app_root(file.path("R", "transcriptomics", "03_Preprocessing_Batch_Correction", "mod_preprocessing_explore.R"))
source_from_app_root(file.path("R", "transcriptomics", "03_Preprocessing_Batch_Correction", "mod_preprocessing.R"))

## Runs the Preprocessing module up through "Run normalisation and batch
## correction" for a single already-loaded dataset (fx$expr/fx$meta), with
## the given correction method, and returns the run's result() list.
## `out <<-` inside shiny::testServer()'s expr assigns into this function's
## own frame - the documented way to pull a value back out of a testServer block.
pp_run_batch_correction <- function(fx, correction_method = "combat", extra_inputs = list()) {
  dataset <- shiny::reactiveValues(expr = fx$expr, meta = fx$meta,
                                    source = "Uploaded dataset: batch_fx.csv",
                                    source_type = "uploaded")
  out <- NULL
  shiny::testServer(mod_preprocessing_server, args = list(id = "pp", dataset = dataset, results = shiny::reactiveValues()), {
    session$setInputs(preloaded_selected = "__current__", preloaded_log2 = "skip")
    session$setInputs(preloaded_run = 1)
    session$setInputs(merge_mode = "own")
    session$setInputs(merge_btn = 1)

    inputs <- utils::modifyList(list(
      color_by = "group", batch_col = "batch", norm_method = "skip",
      skip_combat = FALSE, protect_cols = "group",
      correction_method = correction_method,
      mad_k = 3, min_pct = 0, variance_pct = 0
    ), extra_inputs)
    do.call(session$setInputs, inputs)
    session$setInputs(run_btn = 1)
    out <<- result()
  })
  out
}

## Mean(batch2) - mean(batch1) across every sample and gene: with the balanced
## 2x2 (batch x group) design in fx_batch_signal_data(), the group effect
## contributes equally to both batches and cancels out of this contrast, so it
## isolates the injected batch effect alone.
batch_contrast <- function(expr, meta) {
  b1 <- meta$sample[meta$batch == "batch1"]
  b2 <- meta$sample[meta$batch == "batch2"]
  mean(expr[, b2]) - mean(expr[, b1])
}

## Mean(RA) - mean(HC) restricted to the signal genes: batch is balanced
## within each group here too, so this isolates the injected group effect.
group_contrast <- function(expr, meta, signal_genes) {
  ra <- meta$sample[meta$group == "RA"]
  hc <- meta$sample[meta$group == "HC"]
  mean(expr[signal_genes, ra]) - mean(expr[signal_genes, hc])
}

test_that("ComBat removes the injected batch effect while preserving dimensions, finiteness and biological signal", {
  fx <- fx_batch_signal_data(seed = 21)
  res <- pp_run_batch_correction(fx, correction_method = "combat")

  expect_equal(dim(res$expr_combat), dim(fx$expr))
  expect_true(all(is.finite(res$expr_combat)))

  batch_before <- batch_contrast(res$expr_qnorm, res$meta)
  batch_after  <- batch_contrast(res$expr_combat, res$meta)
  expect_gt(abs(batch_before), 2)             # sanity: the injected offset (3) is actually there pre-correction
  expect_lt(abs(batch_after), abs(batch_before) * 0.15)  # ComBat should remove the vast majority of it

  group_before <- group_contrast(res$expr_qnorm, res$meta, fx$signal_genes)
  group_after  <- group_contrast(res$expr_combat, res$meta, fx$signal_genes)
  expect_gt(abs(group_before), 1)             # sanity: the injected group effect (2) is actually there
  ## Protecting "group" in the model matrix should leave the real signal close to intact,
  ## clearly distinguishing correction from mere data destruction.
  expect_gt(abs(group_after), abs(group_before) * 0.6)
  expect_equal(sign(group_after), sign(fx$group_effect))
})

test_that("limma::removeBatchEffect removes the injected batch effect while preserving dimensions and biological signal", {
  fx <- fx_batch_signal_data(seed = 22)
  res <- pp_run_batch_correction(fx, correction_method = "limma",
                                  extra_inputs = list(show_advanced = TRUE))

  expect_equal(dim(res$expr_combat), dim(fx$expr))
  expect_true(all(is.finite(res$expr_combat)))
  expect_equal(res$correction_method, "limma")

  batch_before <- batch_contrast(res$expr_qnorm, res$meta)
  batch_after  <- batch_contrast(res$expr_combat, res$meta)
  expect_gt(abs(batch_before), 2)
  ## removeBatchEffect() with an explicit design matrix is an exact linear
  ## projection, so the batch contrast should collapse close to zero.
  expect_lt(abs(batch_after), abs(batch_before) * 0.1)

  group_before <- group_contrast(res$expr_qnorm, res$meta, fx$signal_genes)
  group_after  <- group_contrast(res$expr_combat, res$meta, fx$signal_genes)
  expect_gt(abs(group_before), 1)
  expect_gt(abs(group_after), abs(group_before) * 0.6)
  expect_equal(sign(group_after), sign(fx$group_effect))
})

test_that("SVA (surrogate variable analysis) reduces the batch-correlated signal without destroying every gene", {
  ## SVA doesn't use the batch column directly - it estimates hidden sources of
  ## variation from the data itself. With a single, strong, cleanly two-level
  ## confound (batch) that is NOT part of the protected model (only "group"
  ## is protected), its top surrogate variable should track batch closely
  ## enough that regressing it out collapses the batch contrast substantially.
  fx <- fx_batch_signal_data(seed = 23, batch_effect = 4)
  res <- pp_run_batch_correction(fx, correction_method = "sva",
                                  extra_inputs = list(show_advanced = TRUE, sva_n_sv = 1))

  expect_equal(dim(res$expr_combat), dim(fx$expr))
  expect_true(all(is.finite(res$expr_combat)))
  expect_equal(res$correction_method, "sva")

  batch_before <- batch_contrast(res$expr_qnorm, res$meta)
  batch_after  <- batch_contrast(res$expr_combat, res$meta)
  expect_gt(abs(batch_before), 3)
  expect_lt(abs(batch_after), abs(batch_before) * 0.5)

  ## SVA is intentionally not asserted to preserve the group signal as tightly
  ## as ComBat/limma above - it estimates surrogate variables from the data
  ## rather than being told the batch labels, so the strength of protection
  ## for "group" is weaker by construction. It's still checked for gross
  ## data destruction (a real, nonzero signal must remain).
  group_after <- group_contrast(res$expr_combat, res$meta, fx$signal_genes)
  expect_gt(abs(group_after), 0.3)
})
