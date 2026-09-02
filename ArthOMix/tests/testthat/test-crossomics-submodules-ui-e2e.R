## Module 4 (Cross-omics) - UI/E2E: opens every one of the 3 CX_MODULES
## submodule cards (via their "Add" toggle - Sub-modules tab mechanism,
## "cx_" id_prefix) and confirms each renders with no shiny-output-error,
## then verifies one real cross-module data-flow handoff (Dataset tab ->
## Expression and Methylation Integration's real "Run Integration") end to
## end through the actual browser, not just testServer(). Mirrors
## test-transcriptomics-submodules-ui.R's / test-multiomics-submodules-ui.R's
## approach.
##
## Like every other AppDriver-based E2E test in this project, this requires
## a confirmed Supabase test account (ARTHOMIX_TEST_EMAIL/
## ARTHOMIX_TEST_PASSWORD) - skips cleanly without one rather than failing.

skip_if_not_installed("shinytest2")
skip_if_not_installed("chromote")

test_that("every Cross-omics sub-module tab opens and renders with no output error", {
  skip_if(!nzchar(Sys.getenv("ARTHOMIX_TEST_EMAIL")), "no test Supabase account configured (ARTHOMIX_TEST_EMAIL/ARTHOMIX_TEST_PASSWORD)")

  app <- new_app_driver(
    name = "arthomix-cx-submodules",
    height = 900, width = 1400,
    timeout = 90 * 1000,
    load_timeout = 90 * 1000
  )
  on.exit(app$stop(), add = TRUE)
  login_test_user(app)

  app$set_inputs(sidebar_tabs = "crossomics")
  app$wait_for_idle(timeout = 20 * 1000)

  ## Every CX_MODULES config id, in registry order (see R/submodules_registry.R).
  cx_ids <- c("integration", "biomarkerconv", "mrstage")

  app$set_inputs(cx_menu = "Sub-modules")
  app$wait_for_idle(timeout = 20 * 1000)

  for (hid in cx_ids) {
    app$click(paste0("cx_sm_toggle_", hid))
    app$wait_for_idle(timeout = 30 * 1000)
    html <- app$get_html("body")
    expect_false(grepl("shiny-output-error", html, fixed = TRUE), info = hid)
  }
})

test_that("Dataset -> Expression and Methylation Integration data flow: loading real example DEG/DMP data and running Integration produces real results with no output error", {
  skip_if(!nzchar(Sys.getenv("ARTHOMIX_TEST_EMAIL")), "no test Supabase account configured (ARTHOMIX_TEST_EMAIL/ARTHOMIX_TEST_PASSWORD)")

  app <- new_app_driver(
    name = "arthomix-cx-integration-flow",
    height = 900, width = 1400,
    timeout = 90 * 1000,
    load_timeout = 90 * 1000
  )
  on.exit(app$stop(), add = TRUE)
  login_test_user(app)

  app$set_inputs(sidebar_tabs = "crossomics")
  app$wait_for_idle(timeout = 20 * 1000)
  app$set_inputs(cx_menu = "Dataset")
  app$wait_for_idle(timeout = 20 * 1000)
  app$set_inputs(`cx_dataset-source_mode` = "example", `cx_dataset-sex_stratum` = "female", `cx_dataset-meth_level` = "dmp")
  app$click("cx_dataset-load_example_btn")
  app$wait_for_idle(timeout = 30 * 1000)
  app$click("cx_dataset-use_data_btn")
  app$wait_for_idle(timeout = 10 * 1000)

  app$set_inputs(cx_menu = "Sub-modules")
  app$wait_for_idle(timeout = 20 * 1000)
  app$click("cx_sm_toggle_integration")
  app$wait_for_idle(timeout = 20 * 1000)
  ## The Integration card's UI (including this button) is inserted via
  ## insertTab() - wait_for_idle() alone can return before the client
  ## actually finishes binding it (confirmed live for the equivalent DGE
  ## contrast-column case - see helper-setup.R's retry_click() comment),
  ## so app$click() can hit "Cannot find HTML element ... shiny-bound-input"
  ## on the very next line. Retry instead of a single click attempt.
  retry_click(app, "cx_integration-run_integration", timeout = 30)
  app$wait_for_idle(timeout = 60 * 1000)

  html <- app$get_html("body")
  expect_false(grepl("shiny-output-error", html, fixed = TRUE))
  expect_true(grepl("Integration", html, fixed = TRUE))
})
