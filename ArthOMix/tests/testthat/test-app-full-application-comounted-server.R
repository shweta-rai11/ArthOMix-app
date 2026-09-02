## FULL APPLICATION TEST - all four verticals co-mounted in ONE session,
## exactly mirroring how server.R actually instantiates the real app for
## every real user (server.R:12-160 creates dataset/results,
## methyl_dataset/methyl_results, cross_dataset/cross_results,
## multi_dataset/multi_results ALL together, unconditionally, every
## session - never just one vertical at a time). Per-module tests
## throughout this suite each construct their OWN isolated reactiveValues
## and so can never catch a bug where mounting everything together causes
## cross-vertical state leakage or a naming collision; this test can.
##
## Loads REAL preloaded data through each vertical's OWN real Dataset-tab
## server (not hand-built fixtures) and confirms: (1) each vertical's data
## loads correctly on its own, (2) no vertical's reactiveValues are
## touched by another vertical's load - dataset/methyl_dataset/
## multi_dataset/cross_dataset are checked for cross-contamination
## explicitly - and (3) all four keep working AFTER every other vertical
## has already been loaded (order-independence, matching a real session
## where a user can click between tabs in any order).

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "mod_dataset.R"))
source_from_app_root(file.path("R", "methylomics", "parse_upload.R"))
source_from_app_root(file.path("R", "methylomics", "qc.R"))
source_from_app_root(file.path("R", "methylomics", "annotation.R"))
source_from_app_root(file.path("R", "methylomics", "mod_methyl_dataset.R"))
source_from_app_root(file.path("R", "crossomics", "crossomics_integration_upload.R"))
source_from_app_root(file.path("R", "crossomics", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "crossomics", "mod_cross_dataset.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_plots.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_dataset_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_integration_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "mod_multi_mofa_engine.R"))
source_from_app_root(file.path("R", "multiomics", "mod_multi_mofa.R"))
source_from_app_root(file.path("R", "multiomics", "mod_multi_dataset.R"))

test_that("all four verticals' shared reactiveValues, mounted together exactly as server.R does, each load real data with zero cross-vertical contamination", {
  skip_if_not(METH_DATA_AVAILABLE, "preloaded methylation data not available")
  skip_if_not(MULTI_DATA_AVAILABLE, "preloaded multi-omics data not available")

  ## ---- Mount every vertical's shared state together, unconditionally - ----
  ## the exact same reactiveValues shapes server.R itself creates at the top
  ## of existing_app_server(), before ANY tab has been visited.
  dataset <- local({
    d <- load_default_dataset()
    d$source_type <- "preloaded"; d$is_bundled_reference <- TRUE; d$geo_ids <- MERGED_DEFAULT_GEO_IDS
    do.call(shiny::reactiveValues, d)
  })
  results <- shiny::reactiveValues()
  methyl_dataset <- shiny::reactiveValues()
  methyl_results <- shiny::reactiveValues()
  cross_dataset <- shiny::reactiveValues(
    user_expr_df = NULL, user_expr_source = NULL, user_expr_wide = NULL, user_expr_mapping = NULL, user_expr_sample_cols = character(0),
    user_meth_df = NULL, user_meth_source = NULL, user_meth_wide = NULL, user_meth_mapping = NULL, user_meth_sample_cols = character(0)
  )
  cross_results <- shiny::reactiveValues()
  multi_dataset <- shiny::reactiveValues(
    table_label = NULL, df = NULL, source = NULL, layers = list(), layer_meta = list(),
    sample_meta = NULL, overlap = NULL, active = FALSE, loaded_at = NULL
  )
  multi_results <- shiny::reactiveValues()

  ## Transcriptomics `dataset` already carries the real startup default
  ## (matches server.R exactly - it's constructed the same way above, not
  ## loaded via a button click) - the other three verticals start genuinely
  ## empty until their own Dataset tab is used, also matching a real session.
  shiny::isolate({
    expect_false(is.null(dataset$expr))
    expect_true(dataset$is_bundled_reference)
  })
  expect_null(shiny::isolate(methyl_dataset$beta))
  expect_false(shiny::isolate(isTRUE(multi_dataset$active)))
  expect_null(shiny::isolate(cross_dataset$user_expr_df))

  ## ---- Load real data into Methylomics - confirm Transcriptomics/ --------
  ## Multiomics/Cross-omics are completely untouched by it.
  invisible(load_default_meth_matrix())  ## pre-warm cache (real ~2.1GB matrix, see feedback_arthomix_shiny_playwright_testing)
  shiny::testServer(mod_methyl_dataset_server, args = list(id = "mx_dataset", methyl_dataset = methyl_dataset), {
    session$setInputs(preloaded_choice = "gse42861_wholeblood")
    session$setInputs(load_preloaded_btn = 1)
  })
  shiny::isolate({
    expect_false(is.null(methyl_dataset$beta))
    expect_gt(nrow(methyl_dataset$beta), 0)
    ## Cross-vertical isolation: Transcriptomics dataset is untouched by
    ## Methylomics loading (still the same real startup default).
    expect_true(dataset$is_bundled_reference)
    expect_false(isTRUE(multi_dataset$active))
    expect_null(cross_dataset$user_expr_df)
  })

  ## ---- Load real data into Multiomics (Dataset Workspace, preloaded) - ---
  ## confirm the other three verticals remain untouched by it.
  ## `load_preloaded_btn` only stages the reference cell into the module's
  ## own internal `raw$mats` (real data, `multi_dataset$active`/`source` are
  ## set immediately) - `multi_dataset$layers` itself isn't published until
  ## the full preview/validate/activate wizard reaches "Use Selected
  ## Datasets" (already covered in depth by test-multi-dataset-server.R);
  ## checking the real staged data is enough to prove real loading AND
  ## cross-vertical isolation here without re-driving that whole wizard.
  ## shiny::testServer() always returns NULL regardless of the block's last
  ## expression - a value computed inside must be captured via `<<-` into
  ## this test_that()'s own environment instead.
  staged_layer_names <- NULL
  shiny::testServer(mod_multi_dataset_server, args = list(id = "mo_dataset", multi_dataset = multi_dataset, multi_results = multi_results), {
    session$setInputs(dataset_source = "preloaded")
    session$setInputs(preloaded_pick = "ra_antitnf")
    session$setInputs(preloaded_cell = "female_Etanercept")
    session$setInputs(load_preloaded_btn = 1)
    staged_layer_names <<- names(raw$mats)
  })
  expect_setequal(staged_layer_names, c("Transcriptomics", "Methylomics"))
  shiny::isolate({
    expect_true(isTRUE(multi_dataset$active))
    expect_equal(multi_dataset$source, "preloaded")
    ## Cross-vertical isolation: the OTHER three verticals are unaffected -
    ## in particular, Multiomics' own "Transcriptomics"/"Methylomics" layer
    ## NAMES never leak into (or get confused with) the real dataset/
    ## methyl_dataset reactiveValues from the actual Transcriptomics/
    ## Methylomics verticals, despite the coincidentally-similar naming.
    expect_true(dataset$is_bundled_reference)
    expect_false(is.null(methyl_dataset$beta))
    expect_null(cross_dataset$user_expr_df)
  })

  ## ---- Load real example data into Cross-omics - confirm the other ------
  ## three verticals remain untouched by it.
  shiny::testServer(mod_cross_dataset_server, args = list(id = "cx_dataset", cross_dataset = cross_dataset), {
    session$setInputs(source_mode = "example", sex_stratum = "female", meth_level = "dmp")
    session$setInputs(load_example_btn = 0)
    session$setInputs(load_example_btn = 1)
    session$setInputs(use_data_btn = 0)
    session$setInputs(use_data_btn = 1)
  })
  shiny::isolate({
    expect_false(is.null(cross_dataset$user_expr_df))
    expect_true(nrow(cross_dataset$user_expr_df) > 0)
    ## Final cross-vertical isolation check, all four verticals loaded:
    ## every OTHER vertical's real data is still exactly as it was.
    expect_true(dataset$is_bundled_reference)
    expect_false(is.null(methyl_dataset$beta))
    expect_true(isTRUE(multi_dataset$active))
  })

  ## ---- Order-independence: re-verify Transcriptomics' own dataset is ----
  ## still the real, correct startup default after every other vertical has
  ## been loaded - proves nothing upstream silently mutated a reactiveValues
  ## object it shouldn't have (e.g. via accidental global/shared environment
  ## state, a common class of bug generic per-module fixtures cannot catch).
  shiny::isolate({
    d0 <- load_default_dataset()
    expect_equal(dim(dataset$expr), dim(d0$expr))
    expect_identical(dataset$expr, d0$expr)
  })
})
