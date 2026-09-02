## Module 3 (Multiomics) - Multi-omics Integration (DIABLO & SNF) sub-module,
## via testServer(): data-source dispatch (preloaded vs. Active Multi-Omics
## Dataset), real validation/outcome-summary/eligibility reactive wiring,
## and the synchronous portion of the Run DIABLO/SNF click handlers (the
## `validate(need(eligibility$ok, ...))` gate, which fires and can block the
## run BEFORE any async dispatch happens).
##
## KNOWN LIMITATION (disclosed, not silently skipped): DIABLO/SNF/Sex-
## Stratified all dispatch their actual computation through
## `shiny::ExtendedTask` + `promises::future_promise()` against a REAL
## `future::plan(future::multisession, workers = 2)` (global.R) - genuine
## separate OS worker processes, deferred one flush cycle via
## `session$onFlushed(..., once = TRUE)`. There is no cache to pre-warm
## here (unlike mod_methyl_dataset.R's async preloaded-matrix load) since
## every run is freshly parameterized - driving a real cross-process
## future_promise to resolution through testServer() was judged
## impractical/flaky to do reliably, so the actual RESULT of a successful
## DIABLO/SNF run through the Shiny module is not verified here. The
## underlying computation itself (`mi_diablo_run()`/`mi_snf_run()`) already
## has thorough, real, end-to-end coverage in
## test-multi-integration-live-functions.R - this file verifies the
## reactive wiring AROUND that computation instead.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "multiomics", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_plots.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_dataset_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_dataset_plots.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_integration_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_integration_plots.R"))
source_from_app_root(file.path("R", "multiomics", "cohort_harmonization_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_sexstratified_engine.R"))
source_from_app_root(file.path("R", "multiomics", "mod_multi_integration.R"))

mi_live_fixture <- function(n = 20, seed = 1100) {
  set.seed(seed)
  ids <- paste0("S", seq_len(n))
  outcome <- factor(rep(c("HC", "RA"), each = n / 2))
  expr <- matrix(rnorm(n * 10), n, 10, dimnames = list(ids, paste0("g", 1:10)))
  meth <- matrix(rnorm(n * 10), n, 10, dimnames = list(ids, paste0("cg", 1:10)))
  meta <- data.frame(outcome = as.character(outcome), row.names = ids, stringsAsFactors = FALSE)
  list(ids = ids, expr = expr, meth = meth, meta = meta)
}

test_that("mi_dataset() dispatches to the Active Multi-Omics Dataset by default, and reports a clear error when none is active", {
  multi_dataset <- shiny::reactiveValues(active = FALSE)
  shiny::testServer(mod_multi_integration_server, args = list(id = "mi", multi_dataset = multi_dataset), {
    session$setInputs(data_source = "active")
    d <- mi_dataset()
    expect_false(d$ok)
    expect_true(grepl("No Active Multi-Omics Dataset", d$error))
  })
})

test_that("mi_dataset()/mi_val()/mi_outcome() compute real validation and outcome-class summaries on a live 2-block dataset", {
  fx <- mi_live_fixture()
  multi_dataset <- shiny::reactiveValues(active = TRUE, source = "upload", layers = list(Transcriptomics = fx$expr, Methylomics = fx$meth), sample_meta = fx$meta)
  shiny::testServer(mod_multi_integration_server, args = list(id = "mi", multi_dataset = multi_dataset), {
    session$setInputs(data_source = "active")
    d <- mi_dataset()
    expect_true(d$ok)
    v <- mi_val()
    expect_true(v$ok)
    expect_equal(v$n_shared, 20L)
    expect_true(v$reliable_matching)

    session$setInputs(outcome_col = "outcome")
    o <- mi_outcome()
    expect_equal(o$type, "categorical")
    expect_equal(o$n_classes, 2L)
    expect_equal(o$n, 20L)
  })
})

test_that("diablo_elig()/snf_elig() correctly refuse a single-block dataset (both require >=2 blocks)", {
  fx <- mi_live_fixture()
  multi_dataset <- shiny::reactiveValues(active = TRUE, source = "upload", layers = list(Transcriptomics = fx$expr), sample_meta = fx$meta)
  shiny::testServer(mod_multi_integration_server, args = list(id = "mi", multi_dataset = multi_dataset), {
    session$setInputs(data_source = "active", outcome_col = "outcome")
    expect_false(diablo_elig()$ok)
    expect_false(snf_elig()$ok)
  })
})

test_that("clicking 'Run DIABLO' on an ineligible dataset is blocked by validate() BEFORE any async dispatch - diablo_state$submitted stays FALSE", {
  fx <- mi_live_fixture()
  ## Single block -> ineligible for DIABLO (needs >=2).
  multi_dataset <- shiny::reactiveValues(active = TRUE, source = "upload", layers = list(Transcriptomics = fx$expr), sample_meta = fx$meta)
  shiny::testServer(mod_multi_integration_server, args = list(id = "mi", multi_dataset = multi_dataset), {
    session$setInputs(data_source = "active", outcome_col = "outcome")
    expect_false(diablo_elig()$ok)
    session$setInputs(d_run_btn = 0)
    session$setInputs(d_run_btn = 1)
    ## validate(need(...)) throws inside the observer before diablo_state$submitted
    ## is ever set - it must still read its initial FALSE.
    expect_false(isTRUE(diablo_state$submitted))
    expect_null(diablo_state$result)
  })
})

test_that("clicking 'Run SNF' on an ineligible dataset (missing values present) is blocked before dispatch - snf_state$submitted stays FALSE", {
  fx <- mi_live_fixture()
  meth_with_na <- fx$meth
  meth_with_na[1, 1] <- NA
  multi_dataset <- shiny::reactiveValues(active = TRUE, source = "upload", layers = list(Transcriptomics = fx$expr, Methylomics = meth_with_na), sample_meta = fx$meta)
  shiny::testServer(mod_multi_integration_server, args = list(id = "mi", multi_dataset = multi_dataset), {
    session$setInputs(data_source = "active", outcome_col = "outcome")
    expect_false(snf_elig()$ok)
    session$setInputs(s_blocks = c("Transcriptomics", "Methylomics"))
    session$setInputs(s_run_btn = 0)
    session$setInputs(s_run_btn = 1)
    expect_false(isTRUE(snf_state$submitted))
    expect_null(snf_state$result)
  })
})

test_that("clicking 'Run DIABLO' on an ELIGIBLE dataset synchronously flips diablo_state$submitted to TRUE before the (deferred, async) computation itself runs", {
  fx <- mi_live_fixture()
  multi_dataset <- shiny::reactiveValues(active = TRUE, source = "upload", layers = list(Transcriptomics = fx$expr, Methylomics = fx$meth), sample_meta = fx$meta)
  shiny::testServer(mod_multi_integration_server, args = list(id = "mi", multi_dataset = multi_dataset), {
    session$setInputs(data_source = "active", outcome_col = "outcome")
    expect_true(diablo_elig()$ok)
    session$setInputs(d_blocks = c("Transcriptomics", "Methylomics"), d_ncomp = 1, d_folds = 3, d_nrepeat = 1)
    session$setInputs(d_run_btn = 0)
    session$setInputs(d_run_btn = 1)
    ## The observer's own body sets submitted <- TRUE synchronously; the
    ## actual mixOmics fit is deferred via session$onFlushed() to a real
    ## multisession worker and is not awaited here (see file header).
    expect_true(isTRUE(diablo_state$submitted))
  })
})
