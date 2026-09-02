## Module 3 (Multiomics) - Overview / Cohort Harmonization sub-module,
## via testServer(): the "no active dataset" gate, and the real
## "Analyze Cohort" pipeline (real modality overlap/matching computed on a
## live 2-layer dataset) publishing correct n_total/n_matched into the
## shared multi_results$overview every other sub-module's QC scorecard/
## summary table reads.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_plots.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "multiomics_dataset_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "multiomics_dataset_plots.R"))
source_from_app_root(file.path("R", "multiomics", "06_Gene_CpG_Concordance", "multiomics_concordance_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "02_Cohort_Harmonization", "cohort_harmonization_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "02_Cohort_Harmonization", "cohort_harmonization_plots.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_integration_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "02_Cohort_Harmonization", "mod_multi_overview.R"))

test_that("body_ui shows the 'no active dataset' empty state when multi_dataset isn't active yet", {
  multi_dataset <- shiny::reactiveValues(active = FALSE)
  multi_results <- shiny::reactiveValues()
  shiny::testServer(mod_multi_overview_server, args = list(id = "ov", multi_dataset = multi_dataset, multi_results = multi_results), {
    html <- fx_html_text(output$body_ui)
    expect_true(grepl("No Active Multi-Omics Dataset", html))
  })
})

test_that("clicking 'Analyze Cohort' on a real 2-layer live dataset computes real matched-sample counts and publishes them to multi_results$overview", {
  set.seed(1000)
  ids_a <- paste0("S", 1:20)
  ids_b <- paste0("S", 11:25)  ## 10 shared (S11-S20) out of 20/15
  expr <- matrix(rnorm(200), 20, 10, dimnames = list(ids_a, paste0("g", 1:10)))
  meth <- matrix(rnorm(150), 15, 10, dimnames = list(ids_b, paste0("cg", 1:10)))

  multi_dataset <- shiny::reactiveValues(
    active = TRUE, source = "upload",
    layers = list(Transcriptomics = expr, Methylomics = meth),
    layer_meta = list(Transcriptomics = list(omics_type = "rnaseq"), Methylomics = list(omics_type = "methylation")),
    sample_meta = NULL
  )
  multi_results <- shiny::reactiveValues()

  shiny::testServer(mod_multi_overview_server, args = list(id = "ov", multi_dataset = multi_dataset, multi_results = multi_results), {
    session$setInputs(sel_modalities = c("Transcriptomics", "Methylomics"), min_overlap = 3)
    session$setInputs(analyze_btn = 0)
    session$setInputs(analyze_btn = 1)

    h <- harmonization_result()
    expect_true(h$ok)
    expect_equal(h$n_matched, 10L)   ## real intersection: S11-S20
    expect_equal(h$n_total, 25L)     ## real union: S1-S25
    expect_equal(h$matched_summary$status, "Partially matched")

    expect_false(is.null(multi_results$overview))
    expect_true(multi_results$overview$harmonization$ok)
    expect_equal(multi_results$overview$harmonization$n_matched, 10L)
    expect_equal(multi_results$overview$harmonization$n_total, 25L)
    expect_null(multi_results$overview$summary36)  ## only populated for the preloaded source
  })
})

test_that("harmonization_result() restricts to only the modalities selected in sel_modalities", {
  ids <- paste0("S", 1:10)
  expr <- matrix(rnorm(100), 10, 10, dimnames = list(ids, paste0("g", 1:10)))
  meth <- matrix(rnorm(50), 5, 10, dimnames = list(ids[1:5], paste0("cg", 1:10)))

  multi_dataset <- shiny::reactiveValues(
    active = TRUE, source = "upload",
    layers = list(Transcriptomics = expr, Methylomics = meth),
    layer_meta = list(Transcriptomics = list(omics_type = "rnaseq"), Methylomics = list(omics_type = "methylation")),
    sample_meta = NULL
  )
  multi_results <- shiny::reactiveValues()

  shiny::testServer(mod_multi_overview_server, args = list(id = "ov", multi_dataset = multi_dataset, multi_results = multi_results), {
    ## Only Transcriptomics selected -> n_matched should equal ITS OWN full
    ## sample count (10), not the 5-sample cross-modality intersection.
    session$setInputs(sel_modalities = "Transcriptomics", min_overlap = 3)
    session$setInputs(analyze_btn = 0)
    session$setInputs(analyze_btn = 1)

    h <- harmonization_result()
    expect_equal(names(h$descriptors), "Transcriptomics")
    expect_equal(h$n_matched, 10L)
  })
})
