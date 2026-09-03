## Module 1 (Transcriptomics) - Differential Gene Expression tab's server
## logic, via testServer(): contrast validation, method-vs-data-scale
## gating (limma vs DESeq2), and results$dge/dge_runs accumulation.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "01_Data", "crossomics_integration_upload.R"))
source_from_app_root(file.path("R", "provenance.R"))
source_from_app_root(file.path("R", "transcriptomics", "functions", "expression_type.R"))
source_from_app_root(file.path("R", "transcriptomics", "04_Differential_Expression", "mod_dge.R"))

dge_fixture_dataset <- function(n_per_group = 6, seed = 70) {
  set.seed(seed)
  n <- n_per_group * 2
  genes <- paste0("GENE", 1:40)
  samples <- paste0("S", 1:n)
  grp <- rep(c("HC", "RA"), each = n_per_group)
  m <- matrix(rnorm(40 * n, mean = 8, sd = 1.2), 40, n, dimnames = list(genes, samples))
  m[1:5, grp == "RA"] <- m[1:5, grp == "RA"] + 2
  meta <- data.frame(sample = samples, group = grp, sex = rep(c("F", "M"), length.out = n), stringsAsFactors = FALSE)
  shiny::reactiveValues(expr = m, meta = meta, source = "dge test cohort", source_type = "uploaded",
                          is_bundled_reference = FALSE, geo_ids = character(0))
}

dge_counts_dataset <- function(n_per_group = 6, seed = 71) {
  set.seed(seed)
  n <- n_per_group * 2
  genes <- paste0("GENE", 1:40)
  samples <- paste0("S", 1:n)
  grp <- rep(c("HC", "RA"), each = n_per_group)
  lambda <- ifelse(rep(grp, each = 40) == "RA", 900, 500)
  m <- matrix(rpois(40 * n, lambda = lambda), 40, n, dimnames = list(genes, samples))
  meta <- data.frame(sample = samples, group = grp, stringsAsFactors = FALSE)
  shiny::reactiveValues(expr = m, meta = meta, source = "dge counts cohort", source_type = "uploaded",
                          is_bundled_reference = FALSE, geo_ids = character(0))
}

test_that("a limma fit on the pipeline dataset succeeds, produces a well-formed DEG table, and is saved to results$dge/dge_runs", {
  dataset <- dge_fixture_dataset()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    session$setInputs(data_source = "pipeline", method = "limma", contrast_col = "group",
                        ref_group = "HC", comp_group = "RA", padj_cut = 0.05, lfc_cut = 0.1)
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)

    tbl <- fit_result()$table
    expect_true(all(c("gene", "logFC", "P.Value", "adj.P.Val") %in% colnames(tbl)))
    expect_gt(nrow(tbl), 0)
    expect_true(all(tbl$P.Value >= 0 & tbl$P.Value <= 1, na.rm = TRUE))

    expect_false(is.null(results$dge))
    expect_equal(results$dge$method, "limma")
    expect_length(results$dge_runs, 1)
  })
})

test_that("a DESeq2 fit on raw-count-like data succeeds and produces the expected column names", {
  dataset <- dge_counts_dataset()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    session$setInputs(data_source = "pipeline", method = "deseq2", contrast_col = "group",
                        ref_group = "HC", comp_group = "RA", padj_cut = 0.05, lfc_cut = 0.1)
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)

    tbl <- fit_result()$table
    expect_true(all(c("gene", "logFC", "P.Value", "adj.P.Val") %in% colnames(tbl)))
    expect_gt(nrow(tbl), 0)
  })
})

test_that("DESeq2 is rejected on already-normalised (log-scale, negative-free but non-count) data", {
  dataset <- dge_fixture_dataset()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    session$setInputs(data_source = "pipeline", method = "deseq2", contrast_col = "group",
                        ref_group = "HC", comp_group = "RA", padj_cut = 0.05, lfc_cut = 0.1)
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)
    err <- tryCatch(fit_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("DESeq2 needs raw", conditionMessage(err)))
  })
})

test_that("limma is rejected on raw, un-normalised count data", {
  dataset <- dge_counts_dataset()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    session$setInputs(data_source = "pipeline", method = "limma", contrast_col = "group",
                        ref_group = "HC", comp_group = "RA", padj_cut = 0.05, lfc_cut = 0.1)
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)
    err <- tryCatch(fit_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("raw, non-negative sequencing counts", conditionMessage(err)))
  })
})

test_that("DESeq2 is rejected on CPM-like normalised-totals data even though it is non-negative and wide-range", {
  dataset <- dge_counts_dataset()
  m <- shiny::isolate(dataset$expr)
  cpm <- sweep(m, 2, colSums(m), FUN = "/") * 1e6
  shiny::isolate(dataset$expr <- cpm)
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    session$setInputs(data_source = "pipeline", method = "deseq2", contrast_col = "group",
                        ref_group = "HC", comp_group = "RA", padj_cut = 0.05, lfc_cut = 0.1)
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)
    err <- tryCatch(fit_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("TPM/FPKM/CPM-normalised", conditionMessage(err)))
  })
})

test_that("identical reference and comparison levels are rejected", {
  dataset <- dge_fixture_dataset()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    session$setInputs(data_source = "pipeline", method = "limma", contrast_col = "group",
                        ref_group = "HC", comp_group = "HC", padj_cut = 0.05, lfc_cut = 0.1)
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)
    err <- tryCatch(fit_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("must be different", conditionMessage(err)))
  })
})

test_that("fewer than 6 total matching samples is rejected", {
  dataset <- dge_fixture_dataset(n_per_group = 6)
  shiny::isolate(dataset$meta$group[5:12] <- NA)
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    session$setInputs(data_source = "pipeline", method = "limma", contrast_col = "group",
                        ref_group = "HC", comp_group = "RA", padj_cut = 0.05, lfc_cut = 0.1)
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)
    err <- tryCatch(fit_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Fewer than 6 samples", conditionMessage(err)))
  })
})

test_that("fewer than 3 samples in one contrast level is rejected even with >=6 total samples", {
  dataset <- dge_fixture_dataset(n_per_group = 6)
  shiny::isolate(dataset$meta$group <- c(rep("HC", 10), rep("RA", 2)))
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    session$setInputs(data_source = "pipeline", method = "limma", contrast_col = "group",
                        ref_group = "HC", comp_group = "RA", padj_cut = 0.05, lfc_cut = 0.1)
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)
    err <- tryCatch(fit_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("at least 3 samples", conditionMessage(err)))
  })
})

test_that("a declared_data_type = 'raw' on the shared dataset lets DESeq2 run even when the live heuristic would be ambiguous, and blocks limma", {
  set.seed(73)
  n <- 12
  genes <- paste0("GENE", 1:40)
  samples <- paste0("S", 1:n)
  grp <- rep(c("HC", "RA"), each = n / 2)
  m <- matrix(rpois(40 * n, lambda = 20), 40, n, dimnames = list(genes, samples))
  meta <- data.frame(sample = samples, group = grp, stringsAsFactors = FALSE)
  dataset <- shiny::reactiveValues(expr = m, meta = meta, source = "small-panel counts",
                                     source_type = "uploaded", is_bundled_reference = FALSE,
                                     geo_ids = character(0), declared_data_type = "raw")
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    session$setInputs(data_source = "pipeline", method = "deseq2", contrast_col = "group",
                        ref_group = "HC", comp_group = "RA", padj_cut = 0.05, lfc_cut = 0.1)
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)
    tbl <- fit_result()$table
    expect_true(all(c("gene", "logFC", "P.Value", "adj.P.Val") %in% colnames(tbl)))

    session$setInputs(method = "limma")
    session$setInputs(run_btn = 2)
    err <- tryCatch(fit_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("raw, non-negative sequencing counts", conditionMessage(err)))
  })
})

test_that("declaring the wrong data type on this module's own upload path blocks at upload time, before any fit is attempted", {
  fm <- fx_expr_meta(n_genes = 40, n_samples = 12, seed = 74)
  counts_like <- round(exp(fm$expr))
  tpm_like <- sweep(counts_like, 2, colSums(counts_like), FUN = "/") * 1e6
  dir <- withr::local_tempdir()
  expr_path <- file.path(dir, "expr.csv")
  meta_path <- file.path(dir, "meta.csv")
  fx_write_expr_csv(tpm_like, expr_path)
  fx_write_meta_csv(fm$meta, meta_path)

  dataset <- shiny::reactiveValues(expr = NULL, meta = NULL, source = NULL, source_type = "preloaded")
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    session$setInputs(data_source = "upload")
    session$setInputs(dge_expr_file = fx_mkfile(expr_path))
    session$setInputs(dge_meta_file = fx_mkfile(meta_path))
    session$setInputs(map_sample_id = "sample")
    session$setInputs(dge_declared_data_type = "raw")

    expect_true(grepl("TPM/FPKM/CPM-normalized", fx_html_text(output$upload_type_warning_ui)))

    src <- cur_source()
    expect_equal(src$mode, "pipeline")
  })
})

test_that("run_dge_now() runs a contrast directly (no run_btn click) and writes the same results$dge/dge_runs as the button path", {
  dataset <- dge_fixture_dataset()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    session$setInputs(data_source = "pipeline")
    out <- run_dge_now(contrast_col = "group", ref_group = "HC", comp_group = "RA", method = "limma",
                        padj_cut = 0.05, lfc_cut = 0.1)

    expect_null(session$input$run_btn)
    expect_false(is.null(results$dge))
    expect_equal(results$dge$method, "limma")
    expect_length(results$dge_runs, 1)
    expect_equal(out, results$dge)
  })
})

test_that("run_dge_now() propagates a validate() failure as a plain catchable error, same message as the button path", {
  dataset <- dge_fixture_dataset()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    session$setInputs(data_source = "pipeline")
    err <- tryCatch(
      run_dge_now(contrast_col = "group", ref_group = "HC", comp_group = "HC", method = "limma"),
      error = function(e) conditionMessage(e)
    )
    expect_true(is.character(err))
    expect_true(grepl("must be different", err))
    expect_null(results$dge)
  })
})

test_that("uploading its own expr+meta pair (bypassing the shared dataset) runs a fit against the uploaded data", {
  fm <- fx_expr_meta(n_genes = 40, n_samples = 12, seed = 72)
  dir <- withr::local_tempdir()
  expr_path <- file.path(dir, "expr.csv")
  meta_path <- file.path(dir, "meta.csv")
  fx_write_expr_csv(fm$expr, expr_path)
  fx_write_meta_csv(fm$meta, meta_path)

  dataset <- shiny::reactiveValues(expr = NULL, meta = NULL, source = NULL, source_type = "preloaded")
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    session$setInputs(data_source = "upload")
    session$setInputs(dge_expr_file = fx_mkfile(expr_path))
    session$setInputs(dge_meta_file = fx_mkfile(meta_path))
    session$setInputs(map_sample_id = "sample")
    session$setInputs(method = "limma", contrast_col = "group", ref_group = "HC", comp_group = "RA",
                        padj_cut = 0.05, lfc_cut = 0.1)
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)

    src <- cur_source()
    expect_equal(src$mode, "upload")
    tbl <- fit_result()$table
    expect_gt(nrow(tbl), 0)
  })
})
