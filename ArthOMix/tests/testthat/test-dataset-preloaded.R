## Module 1 (Transcriptomics) - Dataset tab, "Switch to preloaded data" path.
## Covers the top-level catalog helpers in mod_dataset.R (individual_dataset_entry,
## preloaded_choices, entry_geo_ids - all plain functions, not inside the
## moduleServer closure, so directly unit-testable) plus the load_preloaded_btn
## observer's effect on the shared `dataset` reactiveValues via testServer().

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "mod_dataset.R"))

## ---- Unit: catalog helpers -----------------------------------------------

test_that("individual_dataset_entry() uses the display label when known, else falls back to the raw GSE ID", {
  known <- individual_dataset_entry("GSE93272")
  expect_equal(known$id, "GSE93272")
  expect_equal(known$label, "Whole Blood Training Cohort A")

  unknown <- individual_dataset_entry("GSE99999999")
  expect_equal(unknown$label, "GSE99999999")
})

test_that("PRELOADED_DATASETS catalog contains the default merged entry plus one entry per GEO_SOURCES accession", {
  ids <- vapply(PRELOADED_DATASETS, `[[`, character(1), "id")
  expect_true("__default_merged__" %in% ids)
  for (gse in vapply(GEO_SOURCES, `[[`, character(1), "gse")) {
    expect_true(gse %in% ids, info = gse)
  }
  expect_length(ids, length(GEO_SOURCES) + 1L)
})

test_that("preloaded_choices() returns a named vector of ids keyed by their display label", {
  ch <- preloaded_choices()
  expect_true(is.character(ch))
  expect_true("Merged Data" %in% names(ch))
  expect_equal(unname(ch[["Merged Data"]]), "__default_merged__")
})

test_that("entry_geo_ids() maps the merged entry to both training GSEs, and an individual entry to just itself", {
  expect_equal(entry_geo_ids(default_dataset_entry), c("GSE93272", "GSE110169"))
  indiv <- individual_dataset_entry("GSE15573")
  expect_equal(entry_geo_ids(indiv), "GSE15573")
})

test_that("default_dataset_entry$load() returns the same shape load_default_dataset() does", {
  d <- default_dataset_entry$load()
  expect_true(all(c("expr", "meta", "source") %in% names(d)))
  expect_equal(dim(d$expr), dim(load_default_dataset()$expr))
})

## ---- testServer: load_preloaded_btn observer -----------------------------

test_that("loading the default merged entry activates dataset immediately and marks it as the bundled reference", {
  dataset <- shiny::reactiveValues()
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(preloaded_choice = "__default_merged__")
    session$setInputs(load_preloaded_btn = 1)

    expect_false(is.null(dataset$expr))
    expect_equal(dataset$source_type, "preloaded")
    expect_true(dataset$is_bundled_reference)
    expect_equal(dataset$geo_ids, c("GSE93272", "GSE110169"))
    ## staged_* mirrors the active dataset for this path (see mod_dataset.R's
    ## own comment: the preloaded/GEO/upload handlers activate immediately,
    ## not just stage).
    expect_identical(dataset$staged_expr, dataset$expr)
  })
})

test_that("loading an individual raw GEO entry is NOT marked as the bundled reference and uses its own GSE id", {
  dataset <- shiny::reactiveValues()
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(preloaded_choice = "GSE15573")
    session$setInputs(load_preloaded_btn = 1)

    expect_false(is.null(dataset$expr))
    expect_equal(dataset$source_type, "preloaded")
    expect_false(dataset$is_bundled_reference)
    expect_equal(dataset$geo_ids, "GSE15573")
  })
})

test_that("switching preloaded picks overwrites the previous dataset rather than merging with it", {
  dataset <- shiny::reactiveValues()
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset), {
    session$setInputs(preloaded_choice = "__default_merged__")
    session$setInputs(load_preloaded_btn = 1)
    n_merged <- ncol(dataset$expr)

    session$setInputs(preloaded_choice = "GSE15573")
    session$setInputs(load_preloaded_btn = 2)
    expect_false(dataset$is_bundled_reference)
    expect_equal(dataset$geo_ids, "GSE15573")
    ## Different cohort -> not the same sample count as the merged default
    ## (regression guard against a stale-state bug where switching picks
    ## silently kept the previous dataset$expr in place).
    expect_false(identical(ncol(dataset$expr), n_merged))
  })
})
