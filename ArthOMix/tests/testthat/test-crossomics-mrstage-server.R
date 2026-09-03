## Module 4 (Cross-omics) - Cross-Omics MR sub-module, via testServer():
## loading the real precomputed MR-stage instrument results, the upload
## paths (MR results + standalone evidence file), the fast-path reuse of

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "01_Data", "crossomics_integration_upload.R"))
source_from_app_root(file.path("R", "crossomics", "functions", "integration", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "crossomics", "functions", "biomarker_convergence", "crossomics_biomarkerconv_helpers.R"))
source_from_app_root(file.path("R", "crossomics", "04_Cross_Omics_MR", "crossomics_mrstage_helpers.R"))
source_from_app_root(file.path("R", "crossomics", "functions", "integration", "crossomics_integration_plots.R"))
source_from_app_root(file.path("R", "crossomics", "04_Cross_Omics_MR", "mod_cross_mr_stage.R"))

test_that("loading the real precomputed MR-stage results populates mrs$df with real instrument data", {
  skip_if_not(exists("CX_MR_DATA_AVAILABLE") && isTRUE(CX_MR_DATA_AVAILABLE), "MR-stage source data not available")
  cross_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_cross_mr_stage_server, args = list(id = "mr", cross_dataset = cross_dataset), {
    session$setInputs(mr_source = "preloaded")
    session$setInputs(load_mr = 0)
    session$setInputs(load_mr = 1)
    expect_false(is.null(mrs$df))
    expect_true(nrow(mrs$df) > 0)
    expect_true(is.logical(mrs$df$steiger_dir))
  })
})

test_that("join_df() falls back to a fresh cx_bc_load_precomputed() when no matching cross_results$biomarkerconv is published", {
  skip_if_not(CX_BC_DATA_AVAILABLE, "Biomarker Convergence source data not available")
  skip_if_not(exists("CX_MR_DATA_AVAILABLE") && isTRUE(CX_MR_DATA_AVAILABLE), "MR-stage source data not available")
  cross_dataset <- shiny::reactiveValues()
  cross_results <- shiny::reactiveValues()
  shiny::testServer(mod_cross_mr_stage_server, args = list(id = "mr", cross_dataset = cross_dataset, cross_results = cross_results), {
    session$setInputs(mr_source = "preloaded", sex = "female", evidence_source = "preloaded")
    session$setInputs(load_mr = 0)
    session$setInputs(load_mr = 1)
    jd <- join_df()
    expect_false(is.null(jd))
    expect_true(all(c("DEG_significant", "eQTL_MR_significant") %in% colnames(jd)))
  })
})

test_that("join_df() takes the fast path and reuses Biomarker Convergence's already-published table when the sex matches", {
  marker_df <- data.frame(gene = "MARKER_GENE_XYZ", DEG_significant = TRUE, DMP_genomewide_significant = FALSE,
                            DMR_significant = FALSE, mQTL_MR_significant = FALSE, eQTL_MR_significant = FALSE, stringsAsFactors = FALSE)
  cross_dataset <- shiny::reactiveValues()
  cross_results <- shiny::reactiveValues(biomarkerconv = list(df = marker_df, sex = "female", run_at = "x"))
  shiny::testServer(mod_cross_mr_stage_server, args = list(id = "mr", cross_dataset = cross_dataset, cross_results = cross_results), {
    session$setInputs(mr_source = "preloaded", sex = "female", evidence_source = "preloaded")
    jd <- join_df()
    expect_identical(jd, marker_df)
  })
})

test_that("join_df() does NOT reuse the published table when the requested sex differs", {
  marker_df <- data.frame(gene = "MARKER_GENE_XYZ", DEG_significant = TRUE)
  cross_dataset <- shiny::reactiveValues()
  cross_results <- shiny::reactiveValues(biomarkerconv = list(df = marker_df, sex = "male", run_at = "x"))
  shiny::testServer(mod_cross_mr_stage_server, args = list(id = "mr", cross_dataset = cross_dataset, cross_results = cross_results), {
    session$setInputs(mr_source = "preloaded", sex = "female", evidence_source = "preloaded")
    jd <- join_df()
    expect_false(identical(jd, marker_df))
  })
})

test_that("categories() (real 5-category classification) runs end-to-end on real MR + real join data", {
  skip_if_not(CX_BC_DATA_AVAILABLE, "Biomarker Convergence source data not available")
  skip_if_not(exists("CX_MR_DATA_AVAILABLE") && isTRUE(CX_MR_DATA_AVAILABLE), "MR-stage source data not available")
  cross_dataset <- shiny::reactiveValues()
  cross_results <- shiny::reactiveValues()
  shiny::testServer(mod_cross_mr_stage_server, args = list(id = "mr", cross_dataset = cross_dataset, cross_results = cross_results), {
    session$setInputs(mr_source = "preloaded", sex = "female", evidence_source = "preloaded")
    session$setInputs(load_mr = 0)
    session$setInputs(load_mr = 1)
    cats <- categories()
    expect_equal(length(cats), 5L)
    expect_true(all(vapply(cats, is.data.frame, logical(1))))
  })
})

test_that("uploading an MR results file and switching mr_source clears any prior load", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("gene,pval", "TP53,0.001", "BRCA1,0.5"), path)
  cross_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_cross_mr_stage_server, args = list(id = "mr", cross_dataset = cross_dataset), {
    session$setInputs(mr_source = "upload")
    session$setInputs(upload_mr_file = fx_mkfile(path))
    session$setInputs(load_mr_upload = 0)
    session$setInputs(load_mr_upload = 1)
    expect_false(is.null(mrs$df))
    expect_setequal(mrs$df$gene, c("TP53", "BRCA1"))

    session$setInputs(mr_source = "preloaded")
    expect_null(mrs$df)
  })
})

test_that("uploading a standalone evidence file computes a real relabeled table via the shared cx_bc_relabel() thresholds", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("gene,DEG_adjP,eQTL_MR_FDR", "TP53,0.001,0.01", "BRCA1,0.9,0.9"), path)
  cross_dataset <- shiny::reactiveValues()
  shiny::testServer(mod_cross_mr_stage_server, args = list(id = "mr", cross_dataset = cross_dataset), {
    session$setInputs(evidence_source = "upload")
    session$setInputs(upload_evidence_file = fx_mkfile(path))
    session$setInputs(load_evidence_upload = 0)
    session$setInputs(load_evidence_upload = 1)

    expect_false(is.null(uploaded_evidence$df))
    tp53 <- uploaded_evidence$df[uploaded_evidence$df$gene == "TP53", ]
    expect_true(tp53$DEG_significant)
    expect_true(tp53$eQTL_MR_significant)
  })
})
