## Drives ComBat, removeBatchEffect and SVA (run_combat()/run_limma()/run_sva() in mod_preprocessing.R) via testServer().
## Uses fx_batch_signal_data() to assert dims preserved, no NaN/Inf, batch effect shrinks and group signal survives.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "functions", "expression_type.R"))
source_from_app_root(file.path("R", "transcriptomics", "01_Data", "mod_dataset.R"))
source_from_app_root(file.path("R", "transcriptomics", "03_Preprocessing_Batch_Correction", "mod_preprocessing_explore.R"))
source_from_app_root(file.path("R", "transcriptomics", "03_Preprocessing_Batch_Correction", "mod_preprocessing.R"))

## Runs Preprocessing through "Run normalisation and batch correction" for one dataset and returns result() (via `out <<-`).
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

## Mean(batch2) - mean(batch1) over all samples and genes; the balanced 2x2 design cancels the group effect.
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
  ## SVA finds batch without labels: with a strong batch outside the protected "group" model, its top SV should track it.
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

  ## SVA is only checked for gross data destruction: it is not told the batch labels, so "group" protection is weaker.
  group_after <- group_contrast(res$expr_combat, res$meta, fx$signal_genes)
  expect_gt(abs(group_after), 0.3)
})

## Guard: the non-TMM branch must not quantile-normalise raw RNA-seq counts ("auto" or "quantile"); it needs TMM+log2-CPM.
fx_raw_counts_batch_data <- function(n_genes = 60, n_per_cell = 5, seed = 1) {
  set.seed(seed)
  genes <- sprintf("GENE%03d", seq_len(n_genes))
  design <- expand.grid(batch = c("batch1", "batch2"), group = c("HC", "RA"), stringsAsFactors = FALSE)
  design <- design[rep(seq_len(nrow(design)), each = n_per_cell), ]
  samples <- sprintf("S%02d", seq_len(nrow(design)))
  design$sample <- samples
  ## Wide dynamic range, mostly-integer values - the signature of raw RNA-seq counts.
  expr <- matrix(stats::rnbinom(n_genes * nrow(design), mu = 500, size = 2), n_genes, nrow(design),
                 dimnames = list(genes, samples))
  meta <- data.frame(sample = samples, group = design$group, batch = design$batch, stringsAsFactors = FALSE)
  list(expr = expr, meta = meta)
}

test_that("normalisation left at 'auto' refuses to quantile-normalise raw-count-like data", {
  fx <- fx_raw_counts_batch_data(seed = 30)
  res <- tryCatch(
    pp_run_batch_correction(fx, correction_method = "combat",
                             extra_inputs = list(norm_method = "auto", skip_combat = TRUE)),
    error = function(e) e
  )
  expect_true(inherits(res, "shiny.silent.error") || inherits(res, "validation"))
  expect_true(grepl("raw, un-normalised values", conditionMessage(res)))
})

test_that("normalisation explicitly forced to 'quantile' still refuses raw-count-like data", {
  fx <- fx_raw_counts_batch_data(seed = 31)
  res <- tryCatch(
    pp_run_batch_correction(fx, correction_method = "combat",
                             extra_inputs = list(norm_method = "quantile", skip_combat = TRUE)),
    error = function(e) e
  )
  expect_true(inherits(res, "shiny.silent.error") || inherits(res, "validation"))
  expect_true(grepl("TMM", conditionMessage(res)))
})

test_that("declared_data_type == 'raw' blocks quantile normalisation even on a matrix that doesn't look count-like by itself", {
  fx <- fx_batch_signal_data(seed = 32)  # log-scale fixture, would NOT trigger the heuristic alone
  dataset <- shiny::reactiveValues(expr = fx$expr, meta = fx$meta,
                                    source = "Uploaded dataset: batch_fx.csv",
                                    source_type = "uploaded", declared_data_type = "raw")
  out <- tryCatch({
    res <- NULL
    shiny::testServer(mod_preprocessing_server, args = list(id = "pp", dataset = dataset, results = shiny::reactiveValues()), {
      session$setInputs(preloaded_selected = "__current__", preloaded_log2 = "skip")
      session$setInputs(preloaded_run = 1)
      session$setInputs(merge_mode = "own")
      session$setInputs(merge_btn = 1)
      session$setInputs(color_by = "group", batch_col = "batch", norm_method = "auto",
                         skip_combat = TRUE, protect_cols = "group", correction_method = "combat",
                         mad_k = 3, min_pct = 0, variance_pct = 0)
      session$setInputs(run_btn = 1)
      res <<- result()
    })
    res
  }, error = function(e) e)
  expect_true(inherits(out, "shiny.silent.error") || inherits(out, "validation"))
})
