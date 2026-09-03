## End-to-end smoke test: launches the real app in a headless browser and
## visits each top-level omics module, asserting no Shiny output error is
## rendered. This is the app-level counterpart to test-data-loaders.R's

skip_if_not_installed("shinytest2")
skip_if_not_installed("chromote")

test_that("the app boots and every top-level module renders with no output error", {
  skip_if(!nzchar(Sys.getenv("ARTHOMIX_TEST_EMAIL")), "no test Supabase account configured (ARTHOMIX_TEST_EMAIL/ARTHOMIX_TEST_PASSWORD)")

  app <- new_app_driver(
    name = "arthomix-smoke",
    height = 900, width = 1400,
    timeout = 60 * 1000,
    load_timeout = 90 * 1000
  )
  on.exit(app$stop(), add = TRUE)
  login_test_user(app)

  expect_no_error_in_dom <- function(label) {
    app$wait_for_idle(timeout = 30 * 1000)
    html <- app$get_html("body")
    expect_false(grepl("shiny-output-error", html, fixed = TRUE), info = label)
  }

  expect_no_error_in_dom("home")

  for (tab in c("transcriptomics", "methylomics", "crossomics", "multiomics")) {
    app$set_inputs(sidebar_tabs = tab)
    expect_no_error_in_dom(tab)
  }
})
