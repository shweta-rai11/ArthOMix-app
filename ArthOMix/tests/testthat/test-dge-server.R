## Module 1 (Transcriptomics) - Differential Gene Expression tab's server
## logic, via testServer(): contrast validation, method-vs-data-scale
## gating (limma vs DESeq2), and results$dge/dge_runs accumulation.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "crossomics_integration_upload.R"))
source_from_app_root(file.path("R", "transcriptomics", "mod_dge.R"))

dge_fixture_dataset <- function(n_per_group = 6, seed = 70) {
  set.seed(seed)
  n <- n_per_group * 2
  genes <- paste0("GENE", 1:40)
  samples <- paste0("S", 1:n)
  grp <- rep(c("HC", "RA"), each = n_per_group)
  m <- matrix(rnorm(40 * n, mean = 8, sd = 1.2), 40, n, dimnames = list(genes, samples))
  ## Inject real signal into the first 5 genes so significance isn't purely noise-driven.
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
    session$setInputs(run_btn = 0)  ## prime: real actionButtons start at 0 in the browser DOM; testServer doesn't do this automatically, and fit_result()'s ignoreInit=TRUE eventReactive needs that starting value to correctly treat the next setInputs as a real click
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
    session$setInputs(run_btn = 0)  ## prime: real actionButtons start at 0 in the browser DOM; testServer doesn't do this automatically, and fit_result()'s ignoreInit=TRUE eventReactive needs that starting value to correctly treat the next setInputs as a real click
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
    session$setInputs(run_btn = 0)  ## prime: real actionButtons start at 0 in the browser DOM; testServer doesn't do this automatically, and fit_result()'s ignoreInit=TRUE eventReactive needs that starting value to correctly treat the next setInputs as a real click
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
    session$setInputs(run_btn = 0)  ## prime: real actionButtons start at 0 in the browser DOM; testServer doesn't do this automatically, and fit_result()'s ignoreInit=TRUE eventReactive needs that starting value to correctly treat the next setInputs as a real click
    session$setInputs(run_btn = 1)
    err <- tryCatch(fit_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("raw, non-negative sequencing counts", conditionMessage(err)))
  })
})

test_that("DESeq2 is rejected on CPM-like normalised-totals data even though it is non-negative and wide-range", {
  dataset <- dge_counts_dataset()
  ## Rescale each sample to sum to 1e6 (CPM), preserving relative gene signal.
  m <- shiny::isolate(dataset$expr)
  cpm <- sweep(m, 2, colSums(m), FUN = "/") * 1e6
  shiny::isolate(dataset$expr <- cpm)
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    session$setInputs(data_source = "pipeline", method = "deseq2", contrast_col = "group",
                        ref_group = "HC", comp_group = "RA", padj_cut = 0.05, lfc_cut = 0.1)
    session$setInputs(run_btn = 0)  ## prime: real actionButtons start at 0 in the browser DOM; testServer doesn't do this automatically, and fit_result()'s ignoreInit=TRUE eventReactive needs that starting value to correctly treat the next setInputs as a real click
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
    session$setInputs(run_btn = 0)  ## prime: real actionButtons start at 0 in the browser DOM; testServer doesn't do this automatically, and fit_result()'s ignoreInit=TRUE eventReactive needs that starting value to correctly treat the next setInputs as a real click
    session$setInputs(run_btn = 1)
    err <- tryCatch(fit_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("must be different", conditionMessage(err)))
  })
})

test_that("fewer than 6 total matching samples is rejected", {
  dataset <- dge_fixture_dataset(n_per_group = 6)
  ## Restrict to only 4 samples by pruning the metadata's group labels.
  shiny::isolate(dataset$meta$group[5:12] <- NA)
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    session$setInputs(data_source = "pipeline", method = "limma", contrast_col = "group",
                        ref_group = "HC", comp_group = "RA", padj_cut = 0.05, lfc_cut = 0.1)
    session$setInputs(run_btn = 0)  ## prime: real actionButtons start at 0 in the browser DOM; testServer doesn't do this automatically, and fit_result()'s ignoreInit=TRUE eventReactive needs that starting value to correctly treat the next setInputs as a real click
    session$setInputs(run_btn = 1)
    err <- tryCatch(fit_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Fewer than 6 samples", conditionMessage(err)))
  })
})

test_that("fewer than 3 samples in one contrast level is rejected even with >=6 total samples", {
  dataset <- dge_fixture_dataset(n_per_group = 6)
  ## 10 HC, 2 RA = 12 total (passes the >=6 check) but RA has < 3.
  shiny::isolate(dataset$meta$group <- c(rep("HC", 10), rep("RA", 2)))
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    session$setInputs(data_source = "pipeline", method = "limma", contrast_col = "group",
                        ref_group = "HC", comp_group = "RA", padj_cut = 0.05, lfc_cut = 0.1)
    session$setInputs(run_btn = 0)  ## prime: real actionButtons start at 0 in the browser DOM; testServer doesn't do this automatically, and fit_result()'s ignoreInit=TRUE eventReactive needs that starting value to correctly treat the next setInputs as a real click
    session$setInputs(run_btn = 1)
    err <- tryCatch(fit_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("at least 3 samples", conditionMessage(err)))
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
    session$setInputs(run_btn = 0)  ## prime: real actionButtons start at 0 in the browser DOM; testServer doesn't do this automatically, and fit_result()'s ignoreInit=TRUE eventReactive needs that starting value to correctly treat the next setInputs as a real click
    session$setInputs(run_btn = 1)

    src <- cur_source()
    expect_equal(src$mode, "upload")
    tbl <- fit_result()$table
    expect_gt(nrow(tbl), 0)
  })
})
