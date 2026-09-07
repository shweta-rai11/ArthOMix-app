## Module 1 (Transcriptomics) - Diagnostic Model's External Validation tab, via testServer():
## the bundled GSE15573 blood cohort is scored with the frozen models trained on the bundled
## reference cohort (no refitting), the result is persisted to results$diagnostic[[sex]]$external
## for the Biomarker Card's external-validation tier, and a provenance record is pushed.
## Added 2026-09-05 in response to the transcriptomics audit (findings 2 and 3).

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "provenance.R"))
source_from_app_root(file.path("R", "transcriptomics", "10_Diagnostic_Model", "mod_diagnostic.R"))

bundled_dataset_rv <- function() {
  d <- load_default_dataset()
  shiny::reactiveValues(expr = d$expr, meta = d$meta, source = d$source %||% "Example dataset: sex-stratified RA blood cohort",
                          source_type = "preloaded", is_bundled_reference = TRUE, geo_ids = c("GSE93272", "GSE110169"))
}

test_that("External Validation scores the frozen female models on bundled GSE15573, persists the result for the Biomarker Card and pushes a provenance record", {
  skip_if_not(file.exists(file.path(RAW_DIR, "GSE15573_raw.rds")), "bundled GSE15573 raw data not present")
  dataset <- bundled_dataset_rv()
  results <- shiny::reactiveValues()
  ext_cohort <- diag_bundled_external_cohort()
  panel <- head(intersect(intersect(c("AIF1", "VPS52", "GNL1", "AP4B1", "CSNK2B", "LSM2", "BRD2", "WDR46", "HLA-DMA", "HLA-DPA1"),
                                    rownames(shiny::isolate(dataset$expr))), rownames(ext_cohort$expr)), 6)
  expect_gte(length(panel), 4)

  shiny::testServer(mod_diagnostic_server, args = list(id = "diag", dataset = dataset, results = results), {
    arthomix_provenance_clear(session)
    session$setInputs(
      panel_source = "own", gene_list = paste(panel, collapse = "\n"),
      ref_group = "HC", comp_group = "RA", test_frac_pct = 30, class_weight_mode = "equal",
      lr_cv_folds = 5, enet_cv_folds = 5, enet_alpha_grid = "0.5", enet_lambda_choice = "lambda.min", enet_nlambda = 30, enet_type_measure = "auc",
      rf_cv_folds = 3, rf_ntree = 100, rf_mtry_mode = "manual", rf_mtry_manual = 2, rf_nodesize = 1, rf_maxnodes_unlimited = TRUE,
      svm_cv_folds = 3, svm_kernel = "linear", svm_cost_mode = "manual", svm_cost_manual = 1, svm_gamma_mode = "auto", svm_tolerance = 0.001,
      run_female_btn = 0
    )
    session$setInputs(run_female_btn = 1)
    fit <- diag_result_female()
    expect_setequal(fit$genes, panel)
    expect_false(is.null(results$diagnostic$female))
    expect_null(results$diagnostic$female$external)

    ## External Validation on the bundled cohort, female panel restricted to female external samples
    session$setInputs(ext_source = "bundled", ext_panel_choice = "female", ext_match_sex = TRUE,
                      ext_hub_auc_thr = 0.85, ext_hub_p_thr = 0.05, run_ext_btn = 0)
    groups <- sort(unique(ext_data()$meta$group))
    expect_setequal(groups, c("HC", "RA"))
    session$setInputs(ext_ref_group = "HC", ext_comp_group = "RA")
    session$setInputs(run_ext_btn = 1)
    r <- ext_result()
    expect_equal(r$cohort_source, "bundled")
    expect_match(r$cohort_label, "GSE15573")
    expect_true(isTRUE(r$sex_restricted))
    expect_equal(r$n_ref, 10L); expect_equal(r$n_comp, 14L)   # GSE15573 females: 10 HC / 14 RA
    expect_setequal(r$genes, panel)
    expect_false(is.null(r$ext_models))
    expect_length(r$ext_models$models, 4)
    for (mm in r$ext_models$models) {
      expect_true(isTRUE(mm$available))
      expect_true(is.finite(mm$auc) && mm$auc >= 0 && mm$auc <= 1)
    }

    ## persisted for the Biomarker Card
    ex <- results$diagnostic$female$external
    expect_false(is.null(ex))
    expect_true(isTRUE(ex$models_scored))
    expect_equal(ex$cohort_source, "bundled")
    expect_length(ex$models, 4)
    expect_true(all(vapply(ex$models, function(m) is.finite(m$auc), logical(1))))
    expect_true(all(c("lr", "enet", "rf", "svm") %in% vapply(ex$models, `[[`, character(1), "key")))

    ## provenance: one training record + one external-validation record
    recs <- arthomix_provenance_records(session)
    mods <- vapply(recs, `[[`, character(1), "module")
    expect_true("mod_diagnostic" %in% mods)
    expect_true("mod_diagnostic_external_validation" %in% mods)
    ext_rec <- recs[[which(mods == "mod_diagnostic_external_validation")[1]]]
    expect_equal(ext_rec$params$stratum, "female")
    expect_equal(ext_rec$params$external_source, "bundled")
    expect_true(isTRUE(ext_rec$params$sex_restricted))
    expect_length(ext_rec$params$model_auc, 4)

    ## retraining invalidates the stale external result
    session$setInputs(run_female_btn = 2)
    invisible(diag_result_female())
    expect_null(results$diagnostic$female$external)
  })
})

test_that("External Validation persists nothing when no trained model exists for the chosen panel (per-gene only is not model validation)", {
  skip_if_not(file.exists(file.path(RAW_DIR, "GSE15573_raw.rds")), "bundled GSE15573 raw data not present")
  dataset <- bundled_dataset_rv()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_diagnostic_server, args = list(id = "diag", dataset = dataset, results = results), {
    session$setInputs(panel_source = "own", gene_list = "AIF1\nVPS52\nGNL1\nBRD2",
                      ext_source = "bundled", ext_panel_choice = "pooled", ext_match_sex = TRUE,
                      ext_hub_auc_thr = 0.85, ext_hub_p_thr = 0.05, run_ext_btn = 0)
    session$setInputs(ext_ref_group = "HC", ext_comp_group = "RA")
    session$setInputs(run_ext_btn = 1)
    r <- ext_result()
    expect_null(r$ext_models)
    expect_false(isTRUE(r$sex_restricted))
    expect_equal(r$n_ref + r$n_comp, 33L)
    expect_true(all(is.finite(r$gr$auc)))
    expect_null(results$diagnostic)
  })
})

test_that("External Validation scores the frozen models on the sealed validation samples reserved on the Dataset tab (pipeline-wide leakage-safe)", {
  source_from_app_root(file.path("R", "transcriptomics", "functions", "expression_type.R"))
  source_from_app_root(file.path("R", "transcriptomics", "01_Data", "mod_dataset.R"))
  d <- load_default_dataset()
  ids <- tx_reserve_validation_ids(d$meta, 0.3, seed = 1234, stratify_sex = TRUE)
  keep <- setdiff(colnames(d$expr), ids)
  dataset <- shiny::reactiveValues(
    expr = d$expr[, keep, drop = FALSE], meta = d$meta[match(keep, d$meta$sample), , drop = FALSE],
    source = paste0(d$source %||% "Example dataset", sprintf(" [%d validation samples reserved]", length(ids))),
    source_type = "preloaded", is_bundled_reference = FALSE, geo_ids = c("GSE93272", "GSE110169"),
    reserved_ids = ids, reserved_expr = d$expr[, ids, drop = FALSE], reserved_meta = d$meta[match(ids, d$meta$sample), , drop = FALSE],
    reserved_info = list(frac = 0.3, seed = 1234, stratify_sex = TRUE, was_bundled_reference = TRUE, resealed = 0L)
  )
  results <- shiny::reactiveValues()
  panel <- c("AIF1", "VPS52", "GNL1", "AP4B1", "CSNK2B", "LSM2")
  expect_true(all(panel %in% rownames(d$expr)))

  shiny::testServer(mod_diagnostic_server, args = list(id = "diag", dataset = dataset, results = results), {
    session$setInputs(
      panel_source = "own", gene_list = paste(panel, collapse = "\n"),
      ref_group = "HC", comp_group = "RA", test_frac_pct = 30, class_weight_mode = "equal",
      lr_cv_folds = 5, enet_cv_folds = 5, enet_alpha_grid = "0.5", enet_lambda_choice = "lambda.min", enet_nlambda = 30, enet_type_measure = "auc",
      rf_cv_folds = 3, rf_ntree = 100, rf_mtry_mode = "manual", rf_mtry_manual = 2, rf_nodesize = 1, rf_maxnodes_unlimited = TRUE,
      svm_cv_folds = 3, svm_kernel = "linear", svm_cost_mode = "manual", svm_cost_manual = 1, svm_gamma_mode = "auto", svm_tolerance = 0.001,
      run_pooled_btn = 0
    )
    session$setInputs(run_pooled_btn = 1)
    fit <- diag_result_pooled()
    ## the model never saw a reserved sample
    expect_length(intersect(ids, colnames(dataset$expr)), 0)
    expect_equal(fit$n_samples + fit$n_test, length(keep))

    session$setInputs(ext_source = "reserved", ext_panel_choice = "pooled", ext_match_sex = TRUE,
                      ext_hub_auc_thr = 0.85, ext_hub_p_thr = 0.05, run_ext_btn = 0)
    expect_setequal(sort(unique(ext_data()$meta$group)), c("HC", "RA"))
    session$setInputs(ext_ref_group = "HC", ext_comp_group = "RA")
    session$setInputs(run_ext_btn = 1)
    r <- ext_result()
    expect_equal(r$cohort_source, "reserved")
    expect_match(r$cohort_label, "sealed validation samples")
    expect_equal(r$n_ref + r$n_comp, length(ids))
    expect_false(isTRUE(r$sex_restricted))
    expect_false(is.null(r$ext_models))
    for (mm in r$ext_models$models) {
      expect_true(isTRUE(mm$available))
      expect_true(is.finite(mm$auc))
    }
    ex <- results$diagnostic$pooled$external
    expect_equal(ex$cohort_source, "reserved")
    expect_true(isTRUE(ex$models_scored))
    recs <- arthomix_provenance_records(session)
    ext_rec <- Filter(function(x) identical(x$module, "mod_diagnostic_external_validation"), recs)
    expect_length(ext_rec, 1)
    expect_equal(ext_rec[[1]]$params$external_source, "reserved")
  })
})
