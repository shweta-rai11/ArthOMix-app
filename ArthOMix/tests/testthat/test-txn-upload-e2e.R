## End-to-end test of the Transcriptomics upload pipeline (Dataset tab ->
## "Upload your own data"), using the chen2021 replication fixtures bundled
## at data/examples/transcriptomics_upload/merged/ (originally authored in

skip_if_not_installed("shinytest2")
skip_if_not_installed("chromote")

test_that("uploading the chen2021 merged fixture completes the full upload -> map -> load flow", {
  skip_if(!nzchar(Sys.getenv("ARTHOMIX_TEST_EMAIL")), "no test Supabase account configured (ARTHOMIX_TEST_EMAIL/ARTHOMIX_TEST_PASSWORD)")

  expr_path <- get_example_path("transcriptomics_upload", "merged", "chen2021_merged_expression_matrix.csv")
  meta_path <- get_example_path("transcriptomics_upload", "merged", "chen2021_merged_sample_metadata.csv")
  expect_true(file.exists(expr_path))
  expect_true(file.exists(meta_path))

  app <- new_app_driver(
    name = "arthomix-upload-tx",
    height = 900, width = 1400,
    timeout = 60 * 1000,
    load_timeout = 90 * 1000
  )
  on.exit(app$stop(), add = TRUE)
  login_test_user(app)

  app$set_inputs(sidebar_tabs = "transcriptomics")
  app$wait_for_idle(timeout = 20 * 1000)
  app$set_inputs(tx_menu = "Dataset")
  app$wait_for_idle(timeout = 20 * 1000)

  app$upload_file(`tx_dataset-expr_file` = expr_path)
  app$upload_file(`tx_dataset-meta_file` = meta_path)
  app$wait_for_idle(timeout = 30 * 1000)

  map_id <- wait_for_input_value(app, "tx_dataset-map_id", timeout = 30)
  map_group <- wait_for_input_value(app, "tx_dataset-map_group", timeout = 30)
  expect_equal(map_id, "sample")
  expect_equal(map_group, "group")

  retry_click(app, "tx_dataset-load_btn", timeout = 30)
  app$wait_for_idle(timeout = 30 * 1000)

  html <- app$get_html("body")
  expect_false(grepl("shiny-output-error", html, fixed = TRUE))
  load_message <- wait_for_html_containing(app, "#tx_dataset-load_message", "Loaded", timeout = 30)
  expect_true(grepl("Loaded", load_message, fixed = TRUE))
  expect_false(grepl("Could not load", load_message, fixed = TRUE))

  app$set_inputs(tx_menu = "Sub-modules")
  app$wait_for_idle(timeout = 20 * 1000)
  app$click("sm_toggle_overview")
  app$wait_for_idle(timeout = 30 * 1000)
  app$set_inputs(`tx_overview-tabs` = "Metadata")
  app$wait_for_idle(timeout = 20 * 1000)
  provenance_html <- wait_for_html_containing(app, "#tx_overview-understand_ui", "Uploaded dataset:", timeout = 30)
  expect_true(grepl("Uploaded dataset:", provenance_html, fixed = TRUE))
})
