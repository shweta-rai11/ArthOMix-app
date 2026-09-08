## Module 1 (Transcriptomics) - Overview and Datasets tab: QC/outlier
## detection and normalisation-check reactives, via testServer().

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "functions", "expression_type.R"))
source_from_app_root(file.path("R", "transcriptomics", "02_Overview", "mod_overview.R"))

fixture_dataset <- function() {
  fm <- fx_expr_meta(n_genes = 200, n_samples = 20, seed = 50)
  shiny::reactiveValues(expr = fm$expr, meta = fm$meta, source = "test cohort",
                          source_type = "uploaded", is_bundled_reference = FALSE, geo_ids = character(0))
}

test_that("run_qc_btn computes sample_qc() against the active dataset and flags/unflags it as stale on dataset change", {
  dataset <- fixture_dataset()
  shiny::testServer(mod_overview_server, args = list(id = "ov", dataset = dataset, results = NULL), {
    session$setInputs(qc_source = "active", mad_k = 3)
    session$setInputs(run_qc_btn = 1)
    qc <- sample_qc()
    expect_equal(nrow(qc), 20L)
    expect_true(all(c("flag_signal", "flag_detected", "flag_cor") %in% colnames(qc)))

    dataset$expr <- fx_expr_meta(n_genes = 200, n_samples = 20, seed = 999)$expr
    session$flushReact()
    expect_true(grepl("Dataset changed", fx_html_text(output$qc_summary_ui)))
  })
})

test_that("running outlier detection writes back to shared results so ArthoChat's grounding sees it as run", {
  ## Regression guard for the 2026-09-07 defense audit finding: mod_overview_server
  ## received results but never wrote to it, so ArthoChat's context builder
  ## (which keys off results[["overview"]] being non-NULL) permanently reported
  ## this sub-module as "NOT YET RUN IN THIS SESSION" even after a real QC run.
  dataset <- fixture_dataset()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_overview_server, args = list(id = "ov", dataset = dataset, results = results), {
    session$setInputs(qc_source = "active", mad_k = 3)
    session$setInputs(run_qc_btn = 0)
    session$setInputs(run_qc_btn = 1)
    sample_qc()
    ov <- shiny::isolate(results$overview)
    expect_false(is.null(ov))
    expect_equal(ov$n_samples, 20L)
  })
})

fixture_raw_counts_dataset <- function() {
  set.seed(52)
  genes <- sprintf("GENE%03d", seq_len(200))
  samples <- sprintf("S%02d", seq_len(20))
  ## Wide dynamic range, mostly-integer values - the signature of raw RNA-seq counts,
  ## the exact input quantile normalisation must never be silently applied to (RED
  ## finding: mod_overview.R's "Apply quantile normalisation" button + styling offered
  ## it as the recommended fix for this exact shape of data).
  expr <- matrix(rnbinom(200 * 20, mu = 500, size = 2), 200, 20, dimnames = list(genes, samples))
  meta <- data.frame(sample = samples, group = rep(c("HC", "RA"), length.out = 20), stringsAsFactors = FALSE)
  shiny::reactiveValues(expr = expr, meta = meta, source = "raw counts test cohort",
                          source_type = "uploaded", is_bundled_reference = FALSE, geo_ids = character(0),
                          declared_data_type = NA_character_)
}

test_that("raw-count-like data is never offered quantile normalisation as the recommended fix", {
  dataset <- fixture_raw_counts_dataset()
  shiny::testServer(mod_overview_server, args = list(id = "ov", dataset = dataset, results = NULL), {
    session$setInputs(qc_source = "active", mad_k = 3)
    session$setInputs(run_norm_btn = 1)

    summary_html <- fx_html_text(output$norm_summary_ui)
    expect_true(grepl("raw, un-normalised sequencing counts", summary_html))
    expect_true(grepl("TMM", summary_html))

    apply_ui_html <- fx_html_text(output$norm_apply_ui)
    expect_true(grepl("disabled here to prevent misuse", apply_ui_html))
    expect_false(grepl("Apply quantile normalisation", apply_ui_html))
  })
})

test_that("norm_apply_result refuses to quantile-normalise raw-count-like data even if the button is driven directly", {
  dataset <- fixture_raw_counts_dataset()
  shiny::testServer(mod_overview_server, args = list(id = "ov", dataset = dataset, results = NULL), {
    session$setInputs(qc_source = "active", mad_k = 3)
    session$setInputs(run_norm_btn = 1)
    session$setInputs(apply_norm_btn = 1)

    res <- tryCatch(norm_apply_result(), error = function(e) e)
    expect_true(inherits(res, "shiny.silent.error") || inherits(res, "validation"))
  })
})

test_that("declared_data_type == 'raw' blocks quantile normalisation even when the matrix itself doesn't look count-like", {
  fm <- fx_expr_meta(n_genes = 200, n_samples = 20, seed = 53)
  dataset <- shiny::reactiveValues(expr = fm$expr, meta = fm$meta, source = "test cohort",
                                     source_type = "uploaded", is_bundled_reference = FALSE, geo_ids = character(0),
                                     declared_data_type = "raw")
  shiny::testServer(mod_overview_server, args = list(id = "ov", dataset = dataset, results = NULL), {
    session$setInputs(qc_source = "active", mad_k = 3)
    session$setInputs(run_norm_btn = 1)
    expect_true(norm_target_raw_like())
  })
})

test_that("an artificially injected outlier sample is actually flagged by sample_qc()", {
  fm <- fx_expr_meta(n_genes = 200, n_samples = 20, seed = 51)
  fm$expr[, 1] <- fm$expr[, 1] + 40
  dataset <- shiny::reactiveValues(expr = fm$expr, meta = fm$meta, source = "test cohort",
                                     source_type = "uploaded", is_bundled_reference = FALSE, geo_ids = character(0))
  shiny::testServer(mod_overview_server, args = list(id = "ov", dataset = dataset, results = NULL), {
    session$setInputs(qc_source = "active", mad_k = 3)
    session$setInputs(run_qc_btn = 1)
    qc <- sample_qc()
    expect_true(qc$flag_signal[1])
  })
})

test_that("run_norm_btn's diagnostic correctly distinguishes an unnormalised (raw counts) cohort from an already-normalised one", {
  set.seed(52)
  raw_counts <- matrix(rpois(200 * 20, lambda = 600), 200, 20,
                        dimnames = list(paste0("G", 1:200), paste0("S", 1:20)))
  meta <- data.frame(sample = colnames(raw_counts), group = rep(c("HC", "RA"), 10), stringsAsFactors = FALSE)
  dataset <- shiny::reactiveValues(expr = raw_counts, meta = meta, source = "raw counts test",
                                     source_type = "uploaded", is_bundled_reference = FALSE, geo_ids = character(0))
  shiny::testServer(mod_overview_server, args = list(id = "ov", dataset = dataset, results = NULL), {
    session$setInputs(qc_source = "active")
    session$setInputs(run_norm_btn = 1)
    d <- norm_check()$diag
    expect_true(needs_quantile_norm(d))
  })
})

test_that("adopting the normalised version writes back into the shared dataset$expr and relabels dataset$source", {
  dataset <- fixture_dataset()
  shiny::testServer(mod_overview_server, args = list(id = "ov", dataset = dataset, results = NULL), {
    session$setInputs(qc_source = "active", norm_color_by = "group")
    session$setInputs(run_norm_btn = 1)
    session$setInputs(apply_norm_btn = 1)
    session$setInputs(adopt_norm_btn = 1)

    expect_identical(dataset$expr, norm_apply_result()$expr_after)
    expect_true(grepl("quantile-normalised", dataset$source))
  })
})

test_that("adopting shows a persistent confirmation banner, and a repeat adopt doesn't double the source suffix", {
  dataset <- fixture_dataset()
  shiny::testServer(mod_overview_server, args = list(id = "ov", dataset = dataset, results = NULL), {
    session$setInputs(qc_source = "active", norm_color_by = "group")
    session$setInputs(run_norm_btn = 1)
    session$setInputs(apply_norm_btn = 1)
    session$setInputs(adopt_norm_btn = 1)

    banner <- paste(as.character(output$adopt_confirmation_ui), collapse = "")
    expect_true(grepl("Done - every sub-module", banner))

    session$setInputs(run_norm_btn = 2)
    session$setInputs(apply_norm_btn = 2)
    session$setInputs(adopt_norm_btn = 2)
    expect_equal(dataset$source, "test cohort (quantile-normalised)")

    session$setInputs(qc_source = "GSE15573")
    expect_error(as.character(output$adopt_confirmation_ui))
  })
})

test_that("adopt_norm_btn is a no-op when qc_source is not 'active' (a read-only reference source is selected)", {
  dataset <- fixture_dataset()
  original_expr <- shiny::isolate(dataset$expr)
  shiny::testServer(mod_overview_server, args = list(id = "ov", dataset = dataset, results = NULL), {
    session$setInputs(qc_source = "GSE15573", norm_color_by = "group")
    session$setInputs(run_norm_btn = 1)
    session$setInputs(apply_norm_btn = 1)
    session$setInputs(adopt_norm_btn = 1)
    expect_identical(dataset$expr, original_expr)
  })
})
