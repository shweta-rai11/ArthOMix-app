## End-to-end test of the Transcriptomics upload pipeline (Dataset tab ->
## "Upload your own data"), using the chen2021 replication fixtures bundled
## at data/examples/transcriptomics_upload/merged/ (originally authored in
## replicate_chen2021/, which stays outside data/ as the fixtures' source -
## see R/transcriptomics/01_Data/mod_dataset.R for the actual upload/mapping/load
## logic under test here).

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
  ## mod_dataset.R's load_btn observer renders "Loaded <n> genes across <n>
  ## samples..." on success (this pipeline is now immediately the active
  ## dataset app-wide, not just a preview - see its own comment), or "Could
  ## not load this dataset: ..." on failure.
  load_message <- app$get_html("#tx_dataset-load_message")
  expect_true(grepl("Loaded", load_message, fixed = TRUE))
  expect_false(grepl("Could not load", load_message, fixed = TRUE))

  ## Isolation: the active dataset (and its provenance) switched immediately -
  ## no detour through Preprocessing needed. Every downstream sub-module reads
  ## this same shared reactiveValues object, so checking it here is equivalent
  ## to checking any of them.
  header_html <- app$get_html("body")
  expect_true(grepl("Your own data", header_html, fixed = TRUE))
})
