## Module 1 (Transcriptomics) - UI/E2E: opens every one of the 16 TX_MODULES
## submodule cards (via their "Add" toggle - Sub-modules tab mechanism) and
## confirms each renders with no shiny-output-error, then verifies one real
## cross-module data-flow handoff (Dataset -> DGE's contrast-column picker)
## end to end through the actual browser, not just testServer().

skip_if_not_installed("shinytest2")
skip_if_not_installed("chromote")

test_that("every Transcriptomics sub-module tab opens and renders with no output error", {
  skip_if(!nzchar(Sys.getenv("ARTHOMIX_TEST_EMAIL")), "no test Supabase account configured (ARTHOMIX_TEST_EMAIL/ARTHOMIX_TEST_PASSWORD)")

  app <- new_app_driver(
    name = "arthomix-tx-submodules",
    height = 900, width = 1400,
    timeout = 90 * 1000,
    load_timeout = 90 * 1000
  )
  on.exit(app$stop(), add = TRUE)
  login_test_user(app)

  app$set_inputs(sidebar_tabs = "transcriptomics")
  app$wait_for_idle(timeout = 20 * 1000)

  ## Every TX_MODULES config id, in registry order (see R/submodules_registry.R).
  tx_ids <- c("overview", "preprocessing", "dge", "wgcna", "candidates", "mr", "coloc",
              "featureselection", "diagnostic", "interaction", "crosstissue", "crossancestry",
              "enrichment", "deconvolution", "nomogram", "biomarkercard")

  app$set_inputs(tx_menu = "Sub-modules")
  app$wait_for_idle(timeout = 20 * 1000)

  for (hid in tx_ids) {
    app$click(paste0("sm_toggle_", hid))
    app$wait_for_idle(timeout = 30 * 1000)
    html <- app$get_html("body")
    expect_false(grepl("shiny-output-error", html, fixed = TRUE), info = hid)
  }
})

test_that("Dataset -> DGE data flow: loading the default preloaded dataset populates DGE's contrast-column picker with real metadata columns", {
  skip_if(!nzchar(Sys.getenv("ARTHOMIX_TEST_EMAIL")), "no test Supabase account configured (ARTHOMIX_TEST_EMAIL/ARTHOMIX_TEST_PASSWORD)")

  app <- new_app_driver(
    name = "arthomix-tx-dge-flow",
    height = 900, width = 1400,
    timeout = 90 * 1000,
    load_timeout = 90 * 1000
  )
  on.exit(app$stop(), add = TRUE)
  login_test_user(app)

  app$set_inputs(sidebar_tabs = "transcriptomics")
  app$wait_for_idle(timeout = 20 * 1000)

  ## The default dataset is already active at startup (server.R initializes
  ## it from load_default_dataset()) - explicitly (re-)loading it via the
  ## Dataset tab's own "Load this dataset" button still proves the full
  ## click -> server -> shared reactiveValues round-trip, not just that the
  ## startup default happens to already be populated.
  app$set_inputs(tx_menu = "Dataset")
  app$wait_for_idle(timeout = 20 * 1000)
  app$set_inputs(`tx_dataset-preloaded_choice` = "__default_merged__")
  app$click("tx_dataset-load_preloaded_btn")
  app$wait_for_idle(timeout = 30 * 1000)
  ## wait_for_idle() alone can return before this renderUI actually reaches
  ## the DOM under load (confirmed live - see helper-setup.R's
  ## wait_for_html_containing() comment) - poll for the real content.
  load_msg <- wait_for_html_containing(app, "#tx_dataset-preloaded_load_message", "Loaded", timeout = 30)
  expect_true(grepl("Loaded", load_msg, fixed = TRUE))

  ## Open Differential Expression and confirm its contrast-column picker
  ## reflects the just-loaded dataset's real metadata (not empty/stale).
  app$set_inputs(tx_menu = "Sub-modules")
  app$wait_for_idle(timeout = 20 * 1000)
  app$click("sm_toggle_dge")
  app$wait_for_idle(timeout = 30 * 1000)

  html <- app$get_html("body")
  expect_false(grepl("shiny-output-error", html, fixed = TRUE))
  contrast_col <- wait_for_input_value(app, "tx_dge-contrast_col", timeout = 30)
  expect_true(!is.null(contrast_col) && nzchar(contrast_col))
})
