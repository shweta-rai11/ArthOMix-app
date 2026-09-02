## Regression guard for a mislabeling bug found in the transcriptomics audit
## (2026-08-26): mod_preprocessing.R's "Use this as the active dataset" button
## used to hardcode dataset$source to start with "Uploaded dataset", even when
## every input to the run was a bundled/preloaded cohort. mod_wgcna.R,
## mod_candidates.R, and mod_overview.R all key real behavior off a
## `grepl("^Uploaded dataset", dataset$source)` check, so that mislabeling
## silently flipped WGCNA's default methodology profile (and revealed an
## upload-only UI panel in Candidates) for data nobody actually uploaded.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "functions", "expression_type.R"))
source_from_app_root(file.path("R", "transcriptomics", "01_Data", "mod_dataset.R"))
source_from_app_root(file.path("R", "transcriptomics", "03_Preprocessing_Batch_Correction", "mod_preprocessing_explore.R"))
source_from_app_root(file.path("R", "transcriptomics", "03_Preprocessing_Batch_Correction", "mod_preprocessing.R"))

test_that("activating a preprocessing run built entirely from bundled cohorts is not labelled 'Uploaded dataset'", {
  d0 <- load_default_dataset()
  dataset <- shiny::reactiveValues(expr = d0$expr, meta = d0$meta, source = d0$source, source_type = "preloaded")

  shiny::testServer(mod_preprocessing_server, args = list(id = "pp", dataset = dataset, results = shiny::reactiveValues()), {
    session$setInputs(preloaded_selected = "__default_merged__", preloaded_log2 = "auto")
    session$setInputs(preloaded_run = 1)
    session$setInputs(merge_mode = "own")
    session$setInputs(merge_btn = 1)
    expect_false(grepl("Uploaded dataset:|NCBI GEO:", merged()$sources))

    session$setInputs(color_by = "group", batch_col = "group", norm_method = "skip",
                       skip_combat = TRUE, mad_k = 3, min_pct = 0, variance_pct = 0)
    session$setInputs(run_btn = 1)
    session$setInputs(activate_btn = 1)

    expect_false(grepl("^Uploaded dataset", dataset$source))
    expect_true(grepl("^Preloaded dataset", dataset$source))
  })
})

test_that("activating a run that includes a genuinely uploaded 'currently loaded dataset' IS labelled uploaded", {
  d0 <- load_default_dataset()
  dataset <- shiny::reactiveValues(
    expr = d0$expr, meta = d0$meta, source = d0$source, source_type = "preloaded",
    staged_expr = d0$expr, staged_meta = d0$meta,
    staged_source = "Uploaded dataset: my_expr.csv + my_meta.csv"
  )

  shiny::testServer(mod_preprocessing_server, args = list(id = "pp", dataset = dataset, results = shiny::reactiveValues()), {
    session$setInputs(preloaded_selected = "__current__", preloaded_log2 = "auto")
    session$setInputs(preloaded_run = 1)
    session$setInputs(merge_mode = "own")
    session$setInputs(merge_btn = 1)
    expect_true(grepl("Uploaded dataset:", merged()$sources))

    session$setInputs(color_by = "group", batch_col = "group", norm_method = "skip",
                       skip_combat = TRUE, mad_k = 3, min_pct = 0, variance_pct = 0)
    session$setInputs(run_btn = 1)
    session$setInputs(activate_btn = 1)

    expect_true(grepl("^Uploaded dataset", dataset$source))
  })
})
