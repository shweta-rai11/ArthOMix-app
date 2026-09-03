## Module 3 (Multiomics) - Results Summary & Reproducibility sub-module.
##
## FIXED (was a KNOWN BUG): mod_multi_summary_config/_ui/_server are fully

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "09_Results_Summary", "mod_multi_summary.R"))

test_that("mod_multi_summary is now referenced in MULTI_MODULES (R/submodules_registry.R) - the Results Summary tab is reachable (fixed - was unreachable)", {
  registry_src <- paste(readLines(file.path(app_dir, "R", "submodules_registry.R")), collapse = "\n")
  expect_true(grepl("mod_multi_summary", registry_src, fixed = TRUE))
  expect_true(grepl("mod_multi_biomarkercard_config", registry_src, fixed = TRUE))
})

test_that("dashboard_ui reactive rollup reports real counts from whatever multi_results actually holds, 'Not loaded' otherwise", {
  multi_results <- shiny::reactiveValues(
    overview = list(harmonization = list(n_matched = 42)),
    biomarker = list(df = data.frame(feature = c("g1", "g2", "g1"))),
    concordance = list(df = data.frame(gene_symbol = c("A", "B", "C")))
  )
  shiny::testServer(mod_multi_summary_server, args = list(id = "sm", multi_dataset = shiny::reactiveValues(), multi_results = multi_results), {
    html <- fx_html_text(output$dashboard_ui)
    expect_true(grepl("42", html))
    expect_true(grepl("1\\.2em[^>]*>2</div>", html))
    expect_true(grepl("None loaded", html))
  })
})
