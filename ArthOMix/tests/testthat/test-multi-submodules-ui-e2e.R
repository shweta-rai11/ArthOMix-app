## Module 3 (Multiomics) - UI/E2E: opens every one of the 7 MULTI_MODULES
## submodule cards (via their "Add" toggle - Sub-modules tab mechanism,
## "mo_" id_prefix) and confirms each renders with no shiny-output-error,

skip_if_not_installed("shinytest2")
skip_if_not_installed("chromote")

test_that("every Multiomics sub-module tab opens and renders with no output error", {
  skip_if(!nzchar(Sys.getenv("ARTHOMIX_TEST_EMAIL")), "no test Supabase account configured (ARTHOMIX_TEST_EMAIL/ARTHOMIX_TEST_PASSWORD)")

  app <- new_app_driver(
    name = "arthomix-mo-submodules",
    height = 900, width = 1400,
    timeout = 90 * 1000,
    load_timeout = 90 * 1000
  )
  on.exit(app$stop(), add = TRUE)
  login_test_user(app)

  app$set_inputs(sidebar_tabs = "multiomics")
  app$wait_for_idle(timeout = 20 * 1000)

  mo_ids <- c("overview", "integration", "stratification", "biomarker", "mapping", "pathway", "biomarkercard")

  app$set_inputs(mo_menu = "Sub-modules")
  app$wait_for_idle(timeout = 20 * 1000)

  for (hid in mo_ids) {
    app$click(paste0("mo_sm_toggle_", hid))
    app$wait_for_idle(timeout = 30 * 1000)
    html <- app$get_html("body")
    expect_false(grepl("shiny-output-error", html, fixed = TRUE), info = hid)
  }
})

test_that("Dataset Workspace -> Overview data flow: loading the preloaded RA anti-TNF cell makes it visible to Cohort Harmonization", {
  skip_if(!nzchar(Sys.getenv("ARTHOMIX_TEST_EMAIL")), "no test Supabase account configured (ARTHOMIX_TEST_EMAIL/ARTHOMIX_TEST_PASSWORD)")

  app <- new_app_driver(
    name = "arthomix-mo-overview-flow",
    height = 900, width = 1400,
    timeout = 90 * 1000,
    load_timeout = 90 * 1000
  )
  on.exit(app$stop(), add = TRUE)
  login_test_user(app)

  app$set_inputs(sidebar_tabs = "multiomics")
  app$wait_for_idle(timeout = 20 * 1000)
  app$set_inputs(mo_menu = "Dataset")
  app$wait_for_idle(timeout = 20 * 1000)
  app$set_inputs(`mo_dataset-preloaded_pick` = "ra_antitnf")
  app$set_inputs(`mo_dataset-preloaded_cell` = "female_Etanercept")
  app$click("mo_dataset-load_preloaded_btn")
  app$wait_for_idle(timeout = 60 * 1000)

  app$set_inputs(mo_menu = "Sub-modules")
  app$wait_for_idle(timeout = 20 * 1000)
  app$click("mo_sm_toggle_overview")
  app$wait_for_idle(timeout = 30 * 1000)

  html <- app$get_html("body")
  expect_false(grepl("shiny-output-error", html, fixed = TRUE))
  expect_true(grepl("Preloaded Dataset", html, fixed = TRUE))
})
