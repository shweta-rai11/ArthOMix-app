## Module 1 (Transcriptomics) - Overview and Datasets tab: QC/outlier
## detection and normalisation-check reactives, via testServer().

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "mod_overview.R"))

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

    ## Changing the active dataset invalidates qc_target() -> marks qc_stale().
    dataset$expr <- fx_expr_meta(n_genes = 200, n_samples = 20, seed = 999)$expr
    session$flushReact()
    ## qc_stale is internal, but its externally-observable effect is
    ## qc_summary_ui showing the "changed - re-run" note instead of a fresh
    ## valueBox; check via the rendered output rather than internal state.
    expect_true(grepl("Dataset changed", fx_html_text(output$qc_summary_ui)))
  })
})

test_that("an artificially injected outlier sample is actually flagged by sample_qc()", {
  fm <- fx_expr_meta(n_genes = 200, n_samples = 20, seed = 51)
  fm$expr[, 1] <- fm$expr[, 1] + 40   ## grossly inflate one sample's signal
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

test_that("adopt_norm_btn is a no-op when qc_source is not 'active' (a read-only reference source is selected)", {
  dataset <- fixture_dataset()
  original_expr <- shiny::isolate(dataset$expr)
  shiny::testServer(mod_overview_server, args = list(id = "ov", dataset = dataset, results = NULL), {
    session$setInputs(qc_source = "GSE15573", norm_color_by = "group")
    session$setInputs(run_norm_btn = 1)
    session$setInputs(apply_norm_btn = 1)
    session$setInputs(adopt_norm_btn = 1)
    ## req(identical(input$qc_source, "active")) in the observer should block
    ## this entirely - the shared dataset must be untouched.
    expect_identical(dataset$expr, original_expr)
  })
})
