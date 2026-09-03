## Module 1 (Transcriptomics) - Cross-Tissue Validation's metadata-aware
## tissue-type validation, exercised through the REAL mod_crosstissue_server
## reactive graph via shiny::testServer() (not just the pure functions in

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "10_Diagnostic_Model", "mod_diagnostic.R"))
source_from_app_root(file.path("R", "transcriptomics", "12_Cross_Tissue_Validation", "mod_crosstissue.R"))

ct_val_upload_fixture <- function(tissue_value, seed = 321, n_per_cell = 4) {
  set.seed(seed)
  n <- n_per_cell * 4
  genes <- paste0("G", 1:120)
  samples <- paste0("VS", 1:n)
  grp <- rep(rep(c("HC", "RA"), each = n_per_cell), 2)
  sex <- rep(c("F", "M"), each = n_per_cell * 2)
  counts <- matrix(rpois(120 * n, lambda = 200), 120, n, dimnames = list(genes, samples))
  expr_df <- data.frame(gene = genes, counts, check.names = FALSE)
  meta_df <- data.frame(sample = samples, sex = sex, group = grp, tissue = tissue_value, stringsAsFactors = FALSE)

  expr_path <- tempfile(fileext = ".csv"); meta_path <- tempfile(fileext = ".csv")
  write.csv(expr_df, expr_path, row.names = FALSE)
  write.csv(meta_df, meta_path, row.names = FALSE)
  list(expr_path = expr_path, meta_path = meta_path, genes = genes)
}

blood_training_dataset <- function() {
  set.seed(11)
  n <- 20
  genes <- paste0("G", 1:120)
  samples <- paste0("BS", 1:n)
  expr <- matrix(rnorm(120 * n, mean = 8, sd = 1.5), 120, n, dimnames = list(genes, samples))
  meta <- data.frame(sample = samples, group = rep(c("HC", "RA"), length.out = n),
                      sex = rep(c("F", "M"), length.out = n),
                      tissue = "Whole Blood", stringsAsFactors = FALSE)
  shiny::reactiveValues(expr = expr, meta = meta, source = "test blood training cohort",
                        source_type = "uploaded", is_bundled_reference = FALSE, geo_ids = character(0))
}

test_that("bundled synovium validation dataset passes the tissue gate against a blood training dataset", {
  dataset <- blood_training_dataset()
  results <- shiny::reactiveValues()

  shiny::testServer(mod_crosstissue_server, args = list(id = "ct", dataset = dataset, results = results), {
    tt <- shiny::isolate(ct_training_tissue())
    vt <- shiny::isolate(ct_validation_tissue())
    expect_equal(tt$classification, "blood")
    expect_equal(vt$classification, "non-blood")
    gate <- shiny::isolate(ct_cross_tissue_gate())
    expect_true(gate$valid)

    html <- paste(as.character(shiny::isolate(output$tissue_validation_status)), collapse = " ")
    expect_match(html, "Valid cross-tissue dataset")
    expect_match(html, "Blood-derived")
    expect_match(html, "Non-blood")
  })
})

test_that("an uploaded blood-derived validation dataset is rejected before any DE fitting runs, and Run produces no result", {
  dataset <- blood_training_dataset()
  results <- shiny::reactiveValues()
  fx <- ct_val_upload_fixture("Peripheral Blood")

  shiny::testServer(mod_crosstissue_server, args = list(id = "ct", dataset = dataset, results = results), {
    session$setInputs(val_source = "upload")
    session$setInputs(val_expr_file = fx_mkfile(fx$expr_path))
    session$setInputs(val_meta_file = fx_mkfile(fx$meta_path))
    session$setInputs(val_map_id = "sample", val_map_sex = "sex", val_map_group = "group")
    session$setInputs(val_ref_group = "HC", val_comp_group = "RA")

    vt <- shiny::isolate(ct_validation_tissue())
    expect_equal(vt$classification, "blood")
    gate <- shiny::isolate(ct_cross_tissue_gate())
    expect_false(gate$valid)
    expect_match(gate$status, "both datasets are blood-derived")

    err <- tryCatch(val_active(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_match(conditionMessage(err), "blood-derived")

    session$setInputs(run_pooled_btn_disc = 0)
    session$flushReact()
    session$setInputs(run_pooled_btn_disc = 1)
    run_err <- tryCatch(ct_result_pooled(), error = function(e) e)
    expect_s3_class(run_err, c("validation", "shiny.silent.error"))
    expect_null(shiny::isolate(results$crosstissue$pooled))

    html <- paste(as.character(shiny::isolate(output$tissue_validation_status)), collapse = " ")
    expect_match(html, "Rejected")
  })
})

test_that("an uploaded non-blood-tissue validation dataset passes the gate and Run produces a real fitted result", {
  dataset <- blood_training_dataset()
  results <- shiny::reactiveValues()
  fx <- ct_val_upload_fixture("Synovial tissue")
  panel <- paste(fx$genes[1:5], collapse = "\n")

  shiny::testServer(mod_crosstissue_server, args = list(id = "ct", dataset = dataset, results = results), {
    session$setInputs(val_source = "upload")
    session$setInputs(val_expr_file = fx_mkfile(fx$expr_path))
    session$setInputs(val_meta_file = fx_mkfile(fx$meta_path))
    session$setInputs(val_map_id = "sample", val_map_sex = "sex", val_map_group = "group")
    session$setInputs(val_ref_group = "HC", val_comp_group = "RA")
    session$setInputs(panel_source = "own", gene_list = panel)

    vt <- shiny::isolate(ct_validation_tissue())
    expect_equal(vt$classification, "non-blood")
    gate <- shiny::isolate(ct_cross_tissue_gate())
    expect_true(gate$valid)

    va <- val_active()
    expect_true(all(c("logcpm", "grp", "sex", "tt") %in% names(va)))

    session$setInputs(run_pooled_btn_disc = 0)
    session$flushReact()
    session$setInputs(run_pooled_btn_disc = 1)
    r <- tryCatch(ct_result_pooled(), error = function(e) e)
    if (inherits(r, "shiny.silent.error")) {
      session$setInputs(run_pooled_btn_disc = 2)
      r <- tryCatch(ct_result_pooled(), error = function(e) e)
    }
    expect_false(inherits(r, "error"))
    expect_true(all(c("lr", "enet", "rf", "svm", "discovery") %in% names(r)))
    expect_true(is.numeric(r$lr$full_auc))

    html <- paste(as.character(shiny::isolate(output$tissue_validation_status)), collapse = " ")
    expect_match(html, "Valid cross-tissue dataset")
  })
})
