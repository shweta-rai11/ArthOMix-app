## End-to-end: the "Start an analysis" guide drives the existing Transcriptomics pages in a real browser.
## It must use the pages' own inputs and run button, so the stored DE result matches a manual run.

skip_if_not_installed("shinytest2")
skip_if_not_installed("chromote")

wf_bar_text <- function(app) gsub("\\s+", " ", tryCatch(app$get_text("#wf_bar"), error = function(e) ""))
wf_wait_bar <- function(app, pattern, timeout = 300) {
  deadline <- Sys.time() + timeout
  repeat {
    txt <- wf_bar_text(app)
    if (grepl(pattern, txt, fixed = TRUE) || Sys.time() >= deadline) return(txt)
    Sys.sleep(1)
  }
}

test_that("the guide walks Transcriptomics from the home page to validation using the existing pages", {
  app <- new_app_driver(name = "workflow-guide", height = 1000, width = 1400,
                        timeout = 120 * 1000, load_timeout = 600 * 1000)
  on.exit(app$stop(), add = TRUE)
  try(app$wait_for_idle(timeout = 60 * 1000), silent = TRUE)

  home <- app$get_html(".home-hero")
  expect_match(home, "Start an analysis", fixed = TRUE)
  expect_match(home, "Explore Modules", fixed = TRUE)
  expect_equal(lengths(regmatches(home, gregexpr("home-hero-flow-step", home))), 7L)

  app$click("home_start_analysis")
  Sys.sleep(2)
  ## wait_ = FALSE: picking a radio in the dialog updates no output, so the default wait would just time out.
  app$set_inputs(wf_layer = "transcriptomics", wf_data_choice = "example", wait_ = FALSE)
  app$click("wf_start_confirm")

  txt <- wf_wait_bar(app, "Step 2 of 7")
  expect_match(txt, "Step 2 of 7", fixed = TRUE)
  expect_match(txt, "Groups (group):", fixed = TRUE)
  expect_equal(app$get_value(input = "tx_menu"), "Overview and Datasets")

  app$click("wf_next")
  txt <- wf_wait_bar(app, "Reference: ")
  expect_equal(app$get_value(input = "tx_menu"), "Differential Expression")
  expect_equal(app$get_value(input = "tx_dge-contrast_col"), "group")

  app$click("wf_next")
  wf_wait_bar(app, "Step 4 of 7")
  app$set_inputs(wf_sex_mode = "pooled")
  app$click("wf_next")
  wf_wait_bar(app, "Step 5 of 7")
  Sys.sleep(6)
  expect_equal(app$get_value(input = "tx_dge-covariate_col"), "sex")
  expect_equal(app$get_value(input = "tx_dge-covariate_mode"), "adjust")

  app$click("wf_run")
  txt <- wf_wait_bar(app, "Analysis finished", timeout = 900)
  expect_match(txt, "genes tested", fixed = TRUE)
  expect_match(txt, "adjusted for sex", fixed = TRUE)

  app$click("wf_next")
  wf_wait_bar(app, "Step 6 of 7")
  Sys.sleep(3)
  expect_equal(app$get_value(input = "tx_menu"), "Candidate Gene Identification")

  app$click("wf_next")
  wf_wait_bar(app, "Step 7 of 7")
  Sys.sleep(6)
  expect_equal(app$get_value(input = "tx_menu"), "Diagnostic Model")

  app$click("wf_next")
  expect_match(wf_wait_bar(app, "Guided analysis complete", timeout = 30), "Guided analysis complete", fixed = TRUE)
  expect_false(grepl("shiny-output-error-error", app$get_html("#wf_bar"), fixed = TRUE))
})
