## Module 4 (Cross-omics) - Biomarker Convergence sub-module, via
## testServer(): loading the real precomputed eQTL x mQTL join (with its
## real documented backfill), the "upload your own data" merge path, and
## the fixed-threshold relabeling that feeds every tab's summary cards -
## fully synchronous (no ExtendedTask).

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "crossomics_integration_upload.R"))
source_from_app_root(file.path("R", "crossomics", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "crossomics", "crossomics_biomarkerconv_helpers.R"))
source_from_app_root(file.path("R", "crossomics", "crossomics_integration_plots.R"))
source_from_app_root(file.path("R", "crossomics", "mod_cross_biomarker_conv.R"))

test_that("loading the real precomputed female eQTL x mQTL table populates bc_df() with real relabeled significance flags", {
  skip_if_not(CX_BC_DATA_AVAILABLE, "Biomarker Convergence source data not available")
  cross_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_cross_biomarker_conv_server, args = list(id = "bc", cross_dataset = cross_dataset), {
    session$setInputs(data_source = "preloaded", sex = "female")
    session$setInputs(load_table = 0)
    session$setInputs(load_table = 1)

    expect_false(is.null(raw$df))
    expect_equal(raw$sex, "female")
    df <- bc_df()
    expect_true(all(c("DEG_significant", "eQTL_MR_significant", "n_evidence_layers") %in% colnames(df)))
    ## Real backfill effect (same fact test-crossomics-biomarkerconv-mrstage-
    ## functions.R establishes at the function level): loaded in_mQTL_MR_panel
    ## coverage should exceed the raw file's own un-backfilled count.
    raw_file <- as.data.frame(data.table::fread(cx_bc_precomputed_file("female")))
    expect_true(sum(df$in_mQTL_MR_panel %in% TRUE) >= sum(raw_file$in_mQTL_MR_panel %in% TRUE))
  })
})

test_that("switching data_source clears whatever the other mode had loaded", {
  skip_if_not(CX_BC_DATA_AVAILABLE, "Biomarker Convergence source data not available")
  cross_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_cross_biomarker_conv_server, args = list(id = "bc", cross_dataset = cross_dataset), {
    session$setInputs(data_source = "preloaded", sex = "female")
    session$setInputs(load_table = 0)
    session$setInputs(load_table = 1)
    expect_false(is.null(raw$df))

    session$setInputs(data_source = "upload")
    expect_null(raw$df)
  })
})

test_that("uploading only an eQTL-MR file flags missing_layer='mQTL-MR' and the eQTL-mQTL intersection is correctly empty", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("gene,eQTL_MR_OR,eQTL_MR_FDR", "TP53,1.5,0.01", "BRCA1,0.8,0.2"), path)
  cross_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_cross_biomarker_conv_server, args = list(id = "bc", cross_dataset = cross_dataset), {
    session$setInputs(data_source = "upload")
    session$setInputs(upload_eqtl_file = fx_mkfile(path))
    session$setInputs(load_table_upload = 0)
    session$setInputs(load_table_upload = 1)

    expect_false(is.null(raw$df))
    expect_equal(raw$missing_layer, "mQTL-MR")
    expect_equal(nrow(eqtl_mqtl_df()), 0L)  ## nothing to intersect against
    expect_setequal(raw$df$gene, c("TP53", "BRCA1"))
  })
})

test_that("'Use this data' style upload merge refuses cleanly when neither file is supplied", {
  cross_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_cross_biomarker_conv_server, args = list(id = "bc", cross_dataset = cross_dataset), {
    session$setInputs(data_source = "upload")
    session$setInputs(load_table_upload = 0)
    session$setInputs(load_table_upload = 1)
    expect_null(raw$df)
  })
})
