## Module 3 (Multiomics) - Biomarker Discovery sub-module, via testServer():
## data-source dispatch, real block-role guessing (from layer_meta$omics_type),
## sample matching/outcome-class selection restricted to the matched set,
## the unsupervised variance prefilter, final eligibility, and the
## synchronous `validate(need(mb_elig()$ok, ...))` gate on "Run analysis".
##
## KNOWN LIMITATION (disclosed - see test-multi-integration-server.R's
## header): the actual DIABLO fit dispatches through a real
## `future::multisession` ExtendedTask, not awaited here. `mi_diablo_run()`/
## `mb_cv_roc()` already have real, end-to-end coverage in
## test-multi-integration-live-functions.R / test-multi-biomarker-helpers-
## functions.R.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_plots.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "multiomics_dataset_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "multiomics_dataset_plots.R"))
source_from_app_root(file.path("R", "multiomics", "02_Cohort_Harmonization", "cohort_harmonization_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_integration_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_integration_plots.R"))
source_from_app_root(file.path("R", "multiomics", "03_DIABLO_SNF_Integration", "mod_multi_integration.R"))  ## shared mi_warn()/mi_ok()/mi_stop()
source_from_app_root(file.path("R", "multiomics", "05_Biomarker_Discovery", "multiomics_biomarker_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "05_Biomarker_Discovery", "multiomics_biomarker_plots.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_sexstratified_engine.R"))
source_from_app_root(file.path("R", "multiomics", "05_Biomarker_Discovery", "mod_multi_biomarker.R"))

mb_live_fixture <- function(n = 20, seed = 1300) {
  set.seed(seed)
  ids <- paste0("S", seq_len(n))
  outcome <- rep(c("HC", "RA"), each = n / 2)
  expr <- matrix(rnorm(n * 15), n, 15, dimnames = list(ids, paste0("g", 1:15)))
  meth <- matrix(rnorm(n * 15), n, 15, dimnames = list(ids, paste0("cg", 1:15)))
  meta <- data.frame(outcome = outcome, row.names = ids, stringsAsFactors = FALSE)
  list(ids = ids, expr = expr, meth = meth, meta = meta)
}

test_that("block_role_ui's guessing logic correctly identifies Transcriptomics/Methylomics roles from real layer_meta$omics_type", {
  fx <- mb_live_fixture()
  multi_dataset <- shiny::reactiveValues(
    active = TRUE, source = "upload", layers = list(RNA_layer = fx$expr, Meth_layer = fx$meth),
    layer_meta = list(RNA_layer = list(omics_type = "rnaseq"), Meth_layer = list(omics_type = "methylation")),
    sample_meta = fx$meta
  )
  shiny::testServer(mod_multi_biomarker_server, args = list(id = "mb", multi_dataset = multi_dataset), {
    session$setInputs(data_source = "active")
    d <- mb_dataset()
    expect_true(d$ok)
    expect_equal(unname(d$omics_type["RNA_layer"]), "rnaseq")
    expect_equal(unname(d$omics_type["Meth_layer"]), "methylation")
  })
})

test_that("mb_val_raw()/mb_outcome_summary_raw() compute real sample-matching and outcome-class counts once both block roles are assigned", {
  fx <- mb_live_fixture()
  multi_dataset <- shiny::reactiveValues(
    active = TRUE, source = "upload", layers = list(Transcriptomics = fx$expr, Methylomics = fx$meth),
    layer_meta = list(Transcriptomics = list(omics_type = "rnaseq"), Methylomics = list(omics_type = "methylation")),
    sample_meta = fx$meta
  )
  shiny::testServer(mod_multi_biomarker_server, args = list(id = "mb", multi_dataset = multi_dataset), {
    session$setInputs(data_source = "active", transcript_block = "Transcriptomics", methyl_block = "Methylomics", outcome_col = "outcome")
    v <- mb_val_raw()
    expect_true(v$ok)
    expect_equal(v$n_shared, 20L)
    o <- mb_outcome_summary_raw()
    expect_equal(o$type, "categorical")
    expect_equal(o$n_classes, 2L)
  })
})

test_that("mb_eligible_ids() restricts to the classes_selected subset, and mb_final_layers() applies the real variance prefilter to exactly that subset", {
  fx <- mb_live_fixture()
  multi_dataset <- shiny::reactiveValues(
    active = TRUE, source = "upload", layers = list(Transcriptomics = fx$expr, Methylomics = fx$meth),
    layer_meta = list(Transcriptomics = list(omics_type = "rnaseq"), Methylomics = list(omics_type = "methylation")),
    sample_meta = fx$meta
  )
  shiny::testServer(mod_multi_biomarker_server, args = list(id = "mb", multi_dataset = multi_dataset), {
    session$setInputs(data_source = "active", transcript_block = "Transcriptomics", methyl_block = "Methylomics", outcome_col = "outcome")
    session$setInputs(classes_selected = "RA")  ## HC-only samples excluded
    ids <- mb_eligible_ids()
    expect_equal(length(ids), 10L)
    expect_true(all(fx$meta[ids, "outcome"] == "RA"))

    session$setInputs(max_features_t = 5, max_features_m = 5)
    layers <- mb_final_layers()
    expect_equal(ncol(layers$Transcriptomics), 5L)
    expect_equal(nrow(layers$Transcriptomics), 10L)
  })
})

test_that("mb_elig() correctly refuses when the smallest selected outcome class has fewer than 3 samples", {
  fx <- mb_live_fixture(n = 20)
  ## Recode so RA has only 2 members among the matched set (original: rows
  ## 1-10 HC, 11-20 RA - collapse all but 2 of the RA rows down to HC).
  fx$meta$outcome[13:20] <- "HC"
  multi_dataset <- shiny::reactiveValues(
    active = TRUE, source = "upload", layers = list(Transcriptomics = fx$expr, Methylomics = fx$meth),
    layer_meta = list(Transcriptomics = list(omics_type = "rnaseq"), Methylomics = list(omics_type = "methylation")),
    sample_meta = fx$meta
  )
  shiny::testServer(mod_multi_biomarker_server, args = list(id = "mb", multi_dataset = multi_dataset), {
    session$setInputs(data_source = "active", transcript_block = "Transcriptomics", methyl_block = "Methylomics", outcome_col = "outcome")
    session$setInputs(max_features_t = 15, max_features_m = 15)
    out <- mb_elig()
    expect_false(out$ok)
    expect_true(grepl("at least 3 per class", out$reason))
  })
})

test_that("clicking 'Run analysis' while ineligible is blocked by validate() before dispatch - mb_state$submitted stays FALSE", {
  fx <- mb_live_fixture(n = 20)
  fx$meta$outcome[3:20] <- "HC"  ## RA left with only 2 -> ineligible
  multi_dataset <- shiny::reactiveValues(
    active = TRUE, source = "upload", layers = list(Transcriptomics = fx$expr, Methylomics = fx$meth),
    layer_meta = list(Transcriptomics = list(omics_type = "rnaseq"), Methylomics = list(omics_type = "methylation")),
    sample_meta = fx$meta
  )
  shiny::testServer(mod_multi_biomarker_server, args = list(id = "mb", multi_dataset = multi_dataset), {
    session$setInputs(data_source = "active", transcript_block = "Transcriptomics", methyl_block = "Methylomics", outcome_col = "outcome")
    session$setInputs(max_features_t = 15, max_features_m = 15, ncomp = 1, folds = 3, nrepeat = 1)
    expect_false(mb_elig()$ok)
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)
    expect_false(isTRUE(mb_state$submitted))
    expect_null(mb_state$result)
  })
})

test_that("clicking 'Run analysis' while eligible synchronously flips mb_state$submitted to TRUE, snapshotting the real outcome/layers used", {
  fx <- mb_live_fixture()
  multi_dataset <- shiny::reactiveValues(
    active = TRUE, source = "upload", layers = list(Transcriptomics = fx$expr, Methylomics = fx$meth),
    layer_meta = list(Transcriptomics = list(omics_type = "rnaseq"), Methylomics = list(omics_type = "methylation")),
    sample_meta = fx$meta
  )
  shiny::testServer(mod_multi_biomarker_server, args = list(id = "mb", multi_dataset = multi_dataset), {
    session$setInputs(data_source = "active", transcript_block = "Transcriptomics", methyl_block = "Methylomics", outcome_col = "outcome")
    session$setInputs(max_features_t = 15, max_features_m = 15, ncomp = 1, folds = 3, nrepeat = 1, keepx_t = "8", keepx_m = "8")
    expect_true(mb_elig()$ok)
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)
    expect_true(isTRUE(mb_state$submitted))
    expect_equal(length(mb_state$outcome_used), 20L)
  })
})
