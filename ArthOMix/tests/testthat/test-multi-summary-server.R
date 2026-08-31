## Module 3 (Multiomics) - Results Summary & Reproducibility sub-module.
##
## KNOWN BUG (reported, not fixed - see module report): mod_multi_summary_config/
## _ui/_server are fully implemented but NEVER referenced in
## R/submodules_registry.R's MULTI_MODULES list - the single source of
## truth every mounting mechanism uses (server.R:155's server instantiation,
## ui.R:1548's sub-module grid, the sidebar nav/search/jump-to helpers, all
## keyed off MULTI_MODULES). Confirmed by grepping the whole R/ tree for any
## reference to `mod_multi_summary` outside its own definition file, and for
## the literal id "summary" anywhere else in the registry - neither exists.
## This means the "Results Summary & Reproducibility" tab (session
## dashboard, software-versions table, downloadable session bundle) is
## completely unreachable in the running app - a real, user-facing feature
## gap, not just an edge case. The tests below verify the module's own
## logic works correctly in isolation (it does), which only makes the
## registration gap more clearly a bug rather than untested/broken code.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "multiomics", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "mod_multi_summary.R"))

test_that("KNOWN BUG: mod_multi_summary is never referenced in MULTI_MODULES (R/submodules_registry.R) - the Results Summary tab is unreachable", {
  ## Textual check only (not sourcing the registry file) - MULTI_MODULES'
  ## list literal references every OTHER mod_multi_*_config/_ui/_server by
  ## name, which would each need their own file sourced first; the bug is
  ## the literal absence of "mod_multi_summary" from that source file, which
  ## a plain text search over the real, on-disk registry proves directly.
  registry_src <- paste(readLines(file.path(app_dir, "R", "submodules_registry.R")), collapse = "\n")
  expect_false(grepl("mod_multi_summary", registry_src, fixed = TRUE))  ## documents the current (buggy) state
  expect_true(grepl("mod_multi_biomarkercard_config", registry_src, fixed = TRUE))  ## sanity: the file IS the real registry
})

test_that("dashboard_ui reactive rollup reports real counts from whatever multi_results actually holds, 'Not loaded' otherwise", {
  multi_results <- shiny::reactiveValues(
    overview = list(harmonization = list(n_matched = 42)),
    biomarker = list(df = data.frame(feature = c("g1", "g2", "g1"))),
    concordance = list(df = data.frame(gene_symbol = c("A", "B", "C")))
  )
  shiny::testServer(mod_multi_summary_server, args = list(id = "sm", multi_dataset = shiny::reactiveValues(), multi_results = multi_results), {
    html <- fx_html_text(output$dashboard_ui)
    expect_true(grepl("42", html))       ## matched samples
    ## 2 unique biomarker features (g1 deduplicated) - anchored to the
    ## specific card div so this can't accidentally match a stray "2"
    ## digit inside "42" (matched samples) elsewhere in the same page.
    expect_true(grepl("1\\.2em[^>]*>2</div>", html))
    expect_true(grepl("None loaded", html))  ## no integration cell loaded
  })
})
