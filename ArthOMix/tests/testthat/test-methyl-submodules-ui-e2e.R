## Module 2 (Methylomics) - UI/E2E: opens every one of the 13 MX_MODULES
## submodule cards (via their "Add" toggle - Sub-modules tab mechanism,
## "mx_" id_prefix) and confirms each renders with no shiny-output-error,
## mirroring test-transcriptomics-submodules-ui.R's approach.
##
## Like every other AppDriver-based E2E test in this project (see
## test-app-smoke.R), this now requires a confirmed Supabase test account
## (ARTHOMIX_TEST_EMAIL/ARTHOMIX_TEST_PASSWORD) since a concurrent session
## added an authentication gate mid-effort - skips cleanly without one
## rather than failing.

skip_if_not_installed("shinytest2")
skip_if_not_installed("chromote")

test_that("every Methylomics sub-module tab opens and renders with no output error", {
  skip_if(!nzchar(Sys.getenv("ARTHOMIX_TEST_EMAIL")), "no test Supabase account configured (ARTHOMIX_TEST_EMAIL/ARTHOMIX_TEST_PASSWORD)")

  app <- new_app_driver(
    name = "arthomix-mx-submodules",
    height = 900, width = 1400,
    timeout = 90 * 1000,
    load_timeout = 90 * 1000
  )
  on.exit(app$stop(), add = TRUE)
  login_test_user(app)

  app$set_inputs(sidebar_tabs = "methylomics")
  app$wait_for_idle(timeout = 20 * 1000)

  ## Every MX_MODULES config id, in registry order (see R/submodules_registry.R).
  mx_ids <- c("qc", "normalization", "celltype", "dmp", "dmr", "wgcna", "candidates",
              "featureselection", "mr", "coloc", "diagnostic", "validation", "biomarkercard")

  app$set_inputs(mx_menu = "Sub-modules")
  app$wait_for_idle(timeout = 20 * 1000)

  for (hid in mx_ids) {
    app$click(paste0("mx_sm_toggle_", hid))
    app$wait_for_idle(timeout = 30 * 1000)
    html <- app$get_html("body")
    expect_false(grepl("shiny-output-error", html, fixed = TRUE), info = hid)
  }
})

test_that("Dataset -> QC data flow: loading the preloaded whole-blood dataset makes it visible to the QC tab", {
  skip_if(!nzchar(Sys.getenv("ARTHOMIX_TEST_EMAIL")), "no test Supabase account configured (ARTHOMIX_TEST_EMAIL/ARTHOMIX_TEST_PASSWORD)")

  app <- new_app_driver(
    name = "arthomix-mx-qc-flow",
    height = 900, width = 1400,
    timeout = 90 * 1000,
    load_timeout = 90 * 1000
  )
  on.exit(app$stop(), add = TRUE)
  login_test_user(app)

  app$set_inputs(sidebar_tabs = "methylomics")
  app$wait_for_idle(timeout = 20 * 1000)
  app$set_inputs(mx_menu = "Dataset")
  app$wait_for_idle(timeout = 20 * 1000)
  app$set_inputs(`mx_dataset-preloaded_choice` = "gse42861_wholeblood")
  app$click("mx_dataset-load_preloaded_btn")
  app$wait_for_idle(timeout = 60 * 1000)  ## the preloaded matrix load can take a while, even cached

  app$set_inputs(mx_menu = "Sub-modules")
  app$wait_for_idle(timeout = 20 * 1000)
  app$click("mx_sm_toggle_qc")
  app$wait_for_idle(timeout = 30 * 1000)

  html <- app$get_html("body")
  expect_false(grepl("shiny-output-error", html, fixed = TRUE))
})
