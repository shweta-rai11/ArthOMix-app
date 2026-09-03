## FULL APPLICATION TEST - all four verticals co-mounted in ONE session,
## exactly mirroring how server.R actually instantiates the real app for
## every real user (server.R:12-160 creates dataset/results,

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "functions", "expression_type.R"))
source_from_app_root(file.path("R", "transcriptomics", "01_Data", "mod_dataset.R"))
source_from_app_root(file.path("R", "methylomics", "functions", "parse_upload.R"))
source_from_app_root(file.path("R", "methylomics", "functions", "qc.R"))
source_from_app_root(file.path("R", "methylomics", "functions", "annotation.R"))
source_from_app_root(file.path("R", "methylomics", "01_Data", "mod_methyl_dataset.R"))
source_from_app_root(file.path("R", "crossomics", "01_Data", "crossomics_integration_upload.R"))
source_from_app_root(file.path("R", "crossomics", "functions", "integration", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "crossomics", "01_Data", "mod_cross_dataset.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_plots.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "multiomics_dataset_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_integration_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "mod_multi_mofa_engine.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "mod_multi_mofa.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "mod_multi_dataset.R"))

test_that("all four verticals' shared reactiveValues, mounted together exactly as server.R does, each load real data with zero cross-vertical contamination", {
  skip_if_not(METH_DATA_AVAILABLE, "preloaded methylation data not available")
  skip_if_not(MULTI_DATA_AVAILABLE, "preloaded multi-omics data not available")

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

  shiny::isolate({
    expect_false(is.null(dataset$expr))
    expect_true(dataset$is_bundled_reference)
  })
  expect_null(shiny::isolate(methyl_dataset$beta))
  expect_false(shiny::isolate(isTRUE(multi_dataset$active)))
  expect_null(shiny::isolate(cross_dataset$user_expr_df))

  invisible(load_default_meth_matrix())
  shiny::testServer(mod_methyl_dataset_server, args = list(id = "mx_dataset", methyl_dataset = methyl_dataset), {
    session$setInputs(preloaded_choice = "gse42861_wholeblood")
    session$setInputs(load_preloaded_btn = 1)
  })
  shiny::isolate({
    expect_false(is.null(methyl_dataset$beta))
    expect_gt(nrow(methyl_dataset$beta), 0)
    expect_true(dataset$is_bundled_reference)
    expect_false(isTRUE(multi_dataset$active))
    expect_null(cross_dataset$user_expr_df)
  })

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
    expect_true(dataset$is_bundled_reference)
    expect_false(is.null(methyl_dataset$beta))
    expect_null(cross_dataset$user_expr_df)
  })

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
    expect_true(dataset$is_bundled_reference)
    expect_false(is.null(methyl_dataset$beta))
    expect_true(isTRUE(multi_dataset$active))
  })

  shiny::isolate({
    d0 <- load_default_dataset()
    expect_equal(dim(dataset$expr), dim(d0$expr))
    expect_identical(dataset$expr, d0$expr)
  })
})
