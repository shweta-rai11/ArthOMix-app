## Module 4 (Cross-omics) - UI/E2E: opens every one of the 3 CX_MODULES
## submodule cards (via their "Add" toggle - Sub-modules tab mechanism,
## "cx_" id_prefix) and confirms each renders with no shiny-output-error,

skip_if_not_installed("shinytest2")
skip_if_not_installed("chromote")

test_that("every Cross-omics sub-module tab opens and renders with no output error", {

  app <- new_app_driver(
    name = "arthomix-cx-submodules",
    height = 900, width = 1400,
    timeout = 90 * 1000,
    load_timeout = 90 * 1000
  )
  on.exit(app$stop(), add = TRUE)
  app$wait_for_idle(timeout = 20 * 1000)

  app$set_inputs(sidebar_tabs = "crossomics")
  app$wait_for_idle(timeout = 20 * 1000)

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

test_that("Dataset -> Expression and Methylation Integration data flow: uploading DEG/DMP tables and running Integration produces real results with no output error", {
  expr_path <- tempfile(fileext = ".csv"); meth_path <- tempfile(fileext = ".csv")
  writeLines(c("gene_symbol,log2FC,adj.P.Val", "TP53,2.5,0.001", "BRCA1,-1.2,0.02", "EGFR,0.1,0.9", "MYC,1.4,0.01"), expr_path)
  writeLines(c("cpg,gene,delta_beta,fdr", "cg1,TP53,-0.3,0.001", "cg2,TP53,-0.2,0.01", "cg3,BRCA1,0.3,0.001", "cg4,EGFR,0.2,0.001", "cg5,MYC,0.01,0.9"), meth_path)

  app <- new_app_driver(
    name = "arthomix-cx-integration-flow",
    height = 900, width = 1400,
    timeout = 90 * 1000,
    load_timeout = 90 * 1000
  )
  on.exit(app$stop(), add = TRUE)
  app$wait_for_idle(timeout = 20 * 1000)

  app$set_inputs(sidebar_tabs = "crossomics")
  app$wait_for_idle(timeout = 20 * 1000)
  app$set_inputs(cx_menu = "Dataset")
  app$wait_for_idle(timeout = 20 * 1000)
  app$set_inputs(`cx_dataset-source_mode` = "upload")
  app$wait_for_idle(timeout = 10 * 1000)
  app$upload_file(`cx_dataset-expr_file` = expr_path)
  app$wait_for_idle(timeout = 20 * 1000)
  app$upload_file(`cx_dataset-meth_file` = meth_path)
  app$wait_for_idle(timeout = 20 * 1000)
  app$click("cx_dataset-use_data_btn")
  app$wait_for_idle(timeout = 10 * 1000)

  app$set_inputs(cx_menu = "Sub-modules")
  app$wait_for_idle(timeout = 20 * 1000)
  app$click("cx_sm_toggle_integration")
  app$wait_for_idle(timeout = 20 * 1000)
  retry_click(app, "cx_integration-run_integration", timeout = 30)
  app$wait_for_idle(timeout = 60 * 1000)

  html <- app$get_html("body")
  expect_false(grepl("shiny-output-error", html, fixed = TRUE))
  expect_true(grepl("Integration", html, fixed = TRUE))
})
