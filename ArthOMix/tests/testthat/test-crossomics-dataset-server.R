## Module 4 (Cross-omics) - Dataset tab, via testServer(): loading live
## session DGE/DMP results, uploading + auto-standardizing real CSV files,
## the "Use this data" hand-off into the shared cross_dataset store, source-

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "01_Data", "crossomics_integration_upload.R"))
source_from_app_root(file.path("R", "crossomics", "functions", "integration", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "crossomics", "01_Data", "mod_cross_dataset.R"))

## Minimal live-session fixtures: a Transcriptomics DGE run (as save_dge_result() in
## mod_dge.R stores it) and a Methylomics DMP run (as mod_methyl_dmp.R publishes it).
fx_results_live <- function() shiny::reactiveValues(dge_runs = list(
  run1 = list(contrast = "RA vs HC (Female)", method = "limma", n_samples = 40,
              table = data.frame(gene = c("TP53", "BRCA1", "EGFR"), logFC = c(2.5, -1.2, 0.1),
                                 adj.P.Val = c(0.001, 0.02, 0.9), direction = c("Up", "Down", "NS"),
                                 stringsAsFactors = FALSE))
))
fx_methyl_results_live <- function(array_type = "WGBS") shiny::reactiveValues(
  dmp = list(comparison = "RA vs HC (Female)", n_probes = 4L, n_sig = 2L, array_type = array_type),
  dmp_table = data.frame(cpg = c("cg1", "cg2", "cg3", "cg4"), gene = c("TP53", "BRCA1", NA, "EGFR"),
                         dbeta = c(-0.3, 0.3, 0.2, 0.01), p_raw = c(0.001, 0.001, 0.001, 0.9),
                         fdr = c(0.001, 0.001, 0.001, 0.9), chr = "chr1", pos = 1:4,
                         direction = c("hypo", "hyper", "hyper", "hyper"), stringsAsFactors = FALSE)
)

test_that("'My analysis results' loads the chosen live DGE run + latest live DMP run, standardized, dropping un-annotated CpGs", {
  cross_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_cross_dataset_server,
                    args = list(id = "ds", cross_dataset = cross_dataset, results = fx_results_live(), methyl_results = fx_methyl_results_live()), {
    session$setInputs(source_mode = "live", live_dge_run = "run1")
    session$setInputs(load_live_btn = 0)
    session$setInputs(load_live_btn = 1)

    ed <- expr_data(); md <- meth_data()
    expect_false(is.null(ed))
    expect_setequal(ed$df$gene, c("TP53", "BRCA1", "EGFR"))
    expect_true(all(c("gene", "log2fc", "pvalue", "fdr") %in% colnames(ed$df)))
    expect_equal(ed$df$fdr[ed$df$gene == "TP53"], 0.001)
    expect_true(grepl("RA vs HC \\(Female\\)", ed$source))
    expect_false(is.null(md))
    expect_true(all(c("gene", "cpg", "dbeta", "region", "island_context") %in% colnames(md$df)))
    expect_setequal(md$df$cpg, c("cg1", "cg2", "cg4"))
    expect_equal(md$platform, "WGBS (no manifest annotation available)")
  })
})

test_that("'My analysis results' with only one of the two live runs available stages nothing (the load button is not offered)", {
  cross_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_cross_dataset_server,
                    args = list(id = "ds", cross_dataset = cross_dataset, results = fx_results_live(), methyl_results = shiny::reactiveValues()), {
    session$setInputs(source_mode = "live", live_dge_run = "run1")
    expect_null(live_dmp_run())
    expect_null(expr_data()); expect_null(meth_data())
  })
})

test_that("uploading a real CSV auto-detects and standardizes it identically to the live-results path", {
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

test_that("switching source_mode clears whatever the other mode had staged, and 'Clear' resets the published cross_dataset store (incl. the platform label)", {
  cross_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_cross_dataset_server,
                    args = list(id = "ds", cross_dataset = cross_dataset, results = fx_results_live(), methyl_results = fx_methyl_results_live()), {
    session$setInputs(source_mode = "live", live_dge_run = "run1")
    session$setInputs(load_live_btn = 0)
    session$setInputs(load_live_btn = 1)
    expect_false(is.null(expr_data()))

    session$setInputs(source_mode = "upload")
    expect_null(expr_data())

    session$setInputs(source_mode = "live")
    session$setInputs(load_live_btn = 2)
    session$setInputs(use_data_btn = 0)
    session$setInputs(use_data_btn = 1)
    expect_false(is.null(cross_dataset$user_expr_df))
    expect_equal(cross_dataset$user_meth_platform, "WGBS (no manifest annotation available)")

    session$setInputs(clear_btn = 0)
    session$setInputs(clear_btn = 1)
    expect_null(cross_dataset$user_expr_df)
    expect_null(cross_dataset$user_meth_platform)
    expect_null(expr_data())
  })
})

test_that("a methylation upload with beta values but no Δβ column is refused, not accepted with an all-NA Δβ", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("cpg,gene,beta,fdr", "cg1,TP53,0.8,0.01", "cg2,BRCA1,0.2,0.02"), path)
  cross_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_cross_dataset_server, args = list(id = "ds", cross_dataset = cross_dataset), {
    session$setInputs(source_mode = "upload")
    session$setInputs(meth_file = fx_mkfile(path))
    expect_null(meth_data())
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
