## Module 3 (Multiomics) - Dataset Workspace's reactive server logic, via
## testServer(): the preloaded branch (real mi_preloaded_cell_dataset()
## adapter, no separate special-cased path), the pipeline-switch reset
## (spec: switching sources must never leave the previous pipeline's data
## visible), and the precomputed-table browser (Load button).

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "multiomics", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_plots.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_live_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_integration_live_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "mod_multi_live_mofa.R"))
source_from_app_root(file.path("R", "multiomics", "mod_multi_live.R"))
source_from_app_root(file.path("R", "multiomics", "mod_multi_dataset.R"))

test_that("loading the preloaded reference dataset populates raw$mats with a real 2-layer (Transcriptomics + Methylomics) dataset and publishes multi_dataset$source/active", {
  skip_if_not(MULTI_DATA_AVAILABLE, "multi-omics preloaded data not available in this deployment")
  multi_dataset <- shiny::reactiveValues()
  multi_results <- shiny::reactiveValues()
  shiny::testServer(mod_multi_dataset_server, args = list(id = "ds", multi_dataset = multi_dataset, multi_results = multi_results), {
    session$setInputs(dataset_source = "preloaded")
    session$setInputs(preloaded_pick = "ra_antitnf")
    session$setInputs(preloaded_cell = "female_Etanercept")
    session$setInputs(load_preloaded_btn = 1)

    expect_equal(multi_dataset$source, "preloaded")
    expect_true(multi_dataset$active)
    expect_setequal(names(raw$mats), c("Transcriptomics", "Methylomics"))
    expect_gt(nrow(raw$mats$Transcriptomics), 0)
    expect_gt(nrow(raw$mats$Methylomics), 0)
  })
})

test_that("switching dataset_source resets the previously-published Active dataset and in-progress staging state", {
  skip_if_not(MULTI_DATA_AVAILABLE, "multi-omics preloaded data not available in this deployment")
  multi_dataset <- shiny::reactiveValues()
  multi_results <- shiny::reactiveValues()
  shiny::testServer(mod_multi_dataset_server, args = list(id = "ds", multi_dataset = multi_dataset, multi_results = multi_results), {
    session$setInputs(dataset_source = "preloaded")
    session$setInputs(preloaded_pick = "ra_antitnf")
    session$setInputs(preloaded_cell = "female_Etanercept")
    session$setInputs(load_preloaded_btn = 1)
    expect_true(isTRUE(multi_dataset$active))
    expect_gt(length(raw$mats), 0)

    ## Switching away from "preloaded" must clear both the published dataset
    ## and the staged raw$mats - never leave stale data visible under a
    ## different source selection.
    session$setInputs(dataset_source = "upload")
    expect_false(isTRUE(multi_dataset$active))
    expect_null(multi_dataset$source)
    expect_equal(length(raw$mats), 0)
  })
})

test_that("the precomputed-table browser (table_pick + Load) publishes a real registry table's label into multi_dataset$table_label", {
  skip_if_not(MULTI_DATA_AVAILABLE, "multi-omics preloaded data not available in this deployment")
  multi_dataset <- shiny::reactiveValues()
  multi_results <- shiny::reactiveValues()
  shiny::testServer(mod_multi_dataset_server, args = list(id = "ds", multi_dataset = multi_dataset, multi_results = multi_results), {
    session$setInputs(dataset_source = "preloaded")
    session$setInputs(table_pick = "RNA-seq QC summary")
    ## loaded_table is a bare-button eventReactive(ignoreInit=TRUE) - prime
    ## to 0 first (see feedback_shiny_testserver_ignoreinit_actionbutton_priming.md).
    session$setInputs(load_table_btn = 0)
    session$setInputs(load_table_btn = 1)

    expect_equal(multi_dataset$table_label, "RNA-seq QC summary")
    res <- loaded_table()
    expect_true(res$ok)
    expect_gt(nrow(res$df), 0)
  })
})

test_that("the precomputed-table browser does NOT publish table_label when the chosen table fails to load", {
  multi_dataset <- shiny::reactiveValues()
  multi_results <- shiny::reactiveValues()
  shiny::testServer(mod_multi_dataset_server, args = list(id = "ds", multi_dataset = multi_dataset, multi_results = multi_results), {
    session$setInputs(dataset_source = "preloaded")
    session$setInputs(table_pick = "Not a real registry label")
    session$setInputs(load_table_btn = 0)
    session$setInputs(load_table_btn = 1)

    expect_null(multi_dataset$table_label)
  })
})
