## End-to-end test of the Transcriptomics upload pipeline (Dataset tab ->
## "Upload your own data"), using the chen2021 replication fixtures bundled
## at data/examples/transcriptomics_upload/merged/ (originally authored in
## replicate_chen2021/, which stays outside data/ as the fixtures' source -
## see R/transcriptomics/mod_dataset.R for the actual upload/mapping/load
## logic under test here).

skip_if_not_installed("shinytest2")
skip_if_not_installed("chromote")

test_that("uploading the chen2021 merged fixture completes the full upload -> map -> load flow", {
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

  ## Navigate: top nav -> Transcriptomics -> "Dataset" inner tab.
  app$set_inputs(sidebar_tabs = "transcriptomics")
  app$wait_for_idle(timeout = 20 * 1000)
  app$set_inputs(tx_menu = "Dataset")
  app$wait_for_idle(timeout = 20 * 1000)

  ## STEP 1: choose files (fileInput ids are namespaced under "tx_dataset").
  app$upload_file(`tx_dataset-expr_file` = expr_path)
  app$upload_file(`tx_dataset-meta_file` = meta_path)
  app$wait_for_idle(timeout = 30 * 1000)

  ## STEP 2: column mapping renders and auto-guesses "sample"/"group" from
  ## the fixture's own column names - just confirm the guess landed right
  ## rather than re-picking it, since that guessing logic isn't what this
  ## test is exercising.
  mapping <- app$get_values(input = c("tx_dataset-map_id", "tx_dataset-map_group"))$input
  expect_equal(mapping[["tx_dataset-map_id"]], "sample")
  expect_equal(mapping[["tx_dataset-map_group"]], "group")

  ## STEP 3: confirm.
  app$click("tx_dataset-load_btn")
  app$wait_for_idle(timeout = 30 * 1000)

  html <- app$get_html("body")
  expect_false(grepl("shiny-output-error", html, fixed = TRUE))
  expect_true(grepl("Loaded", app$get_html("#tx_dataset-load_message")))
})
