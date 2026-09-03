## Module 4 (Cross-omics) - Dataset tab, via testServer(): loading real
## example DEG/DMP/DMR data, uploading + auto-standardizing real CSV files,
## the "Use this data" hand-off into the shared cross_dataset store, source-

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "01_Data", "crossomics_integration_upload.R"))
source_from_app_root(file.path("R", "crossomics", "functions", "integration", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "crossomics", "01_Data", "mod_cross_dataset.R"))

test_that("loading real example data (female DEG + DMP) standardizes both tables and previews real counts", {
  skip_if_not(METH_DATA_AVAILABLE, "preloaded methylation data not available")
  cross_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_cross_dataset_server, args = list(id = "ds", cross_dataset = cross_dataset), {
    session$setInputs(source_mode = "example", sex_stratum = "female", meth_level = "dmp")
    session$setInputs(load_example_btn = 0)
    session$setInputs(load_example_btn = 1)

    ed <- expr_data(); md <- meth_data()
    expect_false(is.null(ed))
    expect_true(nrow(ed$df) > 0)
    expect_true(all(c("gene", "log2fc", "pvalue", "fdr") %in% colnames(ed$df)))
    expect_false(is.null(md))
    expect_true(all(c("gene", "cpg", "dbeta") %in% colnames(md$df)))
    expect_true(grepl("FEMALE", ed$source))
  })
})

test_that("uploading a real CSV auto-detects and standardizes it identically to the example-data path", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("gene_symbol,log2FC,adj.P.Val", "TP53,2.5,0.001", "BRCA1,-1.2,0.02"), path)
  cross_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_cross_dataset_server, args = list(id = "ds", cross_dataset = cross_dataset), {
    session$setInputs(source_mode = "upload")
    session$setInputs(expr_file = fx_mkfile(path))
    ed <- expr_data()
    expect_false(is.null(ed))
    expect_setequal(ed$df$gene, c("TP53", "BRCA1"))
    expect_equal(ed$mapping[["gene"]], "gene_symbol")
  })
})

test_that("'Use this data' publishes the standardized tables into the shared cross_dataset store, computing real sample-column detection", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("gene_symbol,log2FC,S1,S2,S3", "TP53,2.5,1.1,1.2,1.3", "BRCA1,-1.2,2.1,2.2,2.3"), path)
  cross_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_cross_dataset_server, args = list(id = "ds", cross_dataset = cross_dataset), {
    session$setInputs(source_mode = "upload")
    session$setInputs(expr_file = fx_mkfile(path))
    session$setInputs(use_data_btn = 0)
    session$setInputs(use_data_btn = 1)

    expect_false(is.null(cross_dataset$user_expr_df))
    expect_setequal(cross_dataset$user_expr_df$gene, c("TP53", "BRCA1"))
    expect_setequal(cross_dataset$user_expr_sample_cols, c("S1", "S2", "S3"))
  })
})

test_that("switching source_mode clears whatever the other mode had staged, and 'Clear' resets the published cross_dataset store", {
  skip_if_not(METH_DATA_AVAILABLE, "preloaded methylation data not available")
  cross_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_cross_dataset_server, args = list(id = "ds", cross_dataset = cross_dataset), {
    session$setInputs(source_mode = "example", sex_stratum = "female", meth_level = "dmp")
    session$setInputs(load_example_btn = 0)
    session$setInputs(load_example_btn = 1)
    expect_false(is.null(expr_data()))

    session$setInputs(source_mode = "upload")
    expect_null(expr_data())

    session$setInputs(load_example_btn = 2)
    session$setInputs(use_data_btn = 0)
    session$setInputs(use_data_btn = 1)
    expect_false(is.null(cross_dataset$user_expr_df))

    session$setInputs(clear_btn = 0)
    session$setInputs(clear_btn = 1)
    expect_null(cross_dataset$user_expr_df)
    expect_null(expr_data())
  })
})

test_that("'Use this data' refuses (via validate()) when nothing has been loaded or uploaded yet", {
  cross_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_cross_dataset_server, args = list(id = "ds", cross_dataset = cross_dataset), {
    session$setInputs(source_mode = "upload")
    session$setInputs(use_data_btn = 0)
    session$setInputs(use_data_btn = 1)
    expect_null(cross_dataset$user_expr_df)
    expect_null(cross_dataset$user_meth_df)
  })
})
