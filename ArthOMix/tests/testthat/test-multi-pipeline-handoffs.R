## Module 3 (Multiomics) - pipeline hand-offs between sub-modules (added 2026-09-06):
## Cohort Harmonisation -> downstream (readiness tier, phenotype default), Integration -> SNF
## Clustering (fused network carried over unchanged), and provenance records from every stage.

suppressWarnings(suppressMessages(source_from_app_root("global.R")))
for (f in list.files(file.path(app_dir, "R", "multiomics"), pattern = "[.]R$", recursive = TRUE, full.names = TRUE)) suppressWarnings(suppressMessages(source(f)))
## the modules read this flag from the global environment; force the synchronous run path for testServer
.ph_async_old <- get0("ARTHOMIX_ASYNC_AVAILABLE", envir = globalenv(), ifnotfound = NULL)
assign("ARTHOMIX_ASYNC_AVAILABLE", FALSE, envir = globalenv())
withr::defer(if (!is.null(.ph_async_old)) assign("ARTHOMIX_ASYNC_AVAILABLE", .ph_async_old, envir = globalenv()), teardown_env())

ph_fixture <- function() {
  skip_if_not(MULTI_DATA_AVAILABLE, "multi-omics preloaded data not available in this deployment")
  ds <- mi_preloaded_cell_dataset("female_response")
  meta <- ds$sample_meta; meta$batch <- rep(c("b1", "b2"), length.out = nrow(meta))
  list(
    multi_dataset = shiny::reactiveValues(source = "upload", active = TRUE, layers = ds$layers, sample_meta = meta,
      layer_meta = list(Transcriptomics = list(omics_type = "rnaseq", status = list(level = "ready")), Methylomics = list(omics_type = "methylation", status = list(level = "ready")))),
    multi_results = shiny::reactiveValues()
  )
}
strip_html <- function(h) gsub("\\s+", " ", gsub("<[^>]+>", " ", as.character(h)))

test_that("Cohort Harmonisation publishes a readiness tier, matched IDs and phenotype/batch candidates, and writes a provenance record", {
  fx <- ph_fixture()
  shiny::testServer(mod_multi_overview_server, args = list(id = "ov", multi_dataset = fx$multi_dataset, multi_results = fx$multi_results), {
    session$setInputs(sel_modalities = c("Transcriptomics", "Methylomics"), min_overlap = 3)
    session$setInputs(analyze_btn = 1)
    h <- multi_results$overview$harmonization
    expect_true(h$ok)
    expect_equal(h$n_matched, 56)
    expect_equal(h$readiness$label, "Ready")
    expect_setequal(h$matched_ids, rownames(shiny::isolate(fx$multi_dataset$layers$Transcriptomics)))
    expect_true("outcome" %in% h$phenotype_candidates)
    expect_true("batch" %in% h$batch_candidates)
    recs <- arthomix_provenance_records(session)
    expect_equal(recs[[length(recs)]]$module, "mod_multi_overview")
  })
})

test_that("Integration reads the harmonisation hand-off (note + outcome default) and publishes its fused network with a provenance record", {
  fx <- ph_fixture()
  fx$multi_results$overview <- list(harmonization = list(ok = TRUE, n_total = 56, n_matched = 56, modalities = c("Transcriptomics", "Methylomics"),
    readiness = list(level = "ready", label = "Ready", reason = "56 matched samples."), n_unmatched_ids = 0L, phenotype_candidates = "outcome", batch_candidates = "batch"))
  shiny::testServer(mod_multi_integration_server, args = list(id = "it", multi_dataset = fx$multi_dataset, multi_results = fx$multi_results), {
    session$setInputs(data_source = "active")
    expect_match(strip_html(output$source_note$html), "Cohort Harmonisation: Ready")
    expect_equal(multi_harmonisation_outcome_default(multi_results, c("batch", "outcome")), "outcome")
    session$setInputs(outcome_col = "outcome", s_blocks = c("Transcriptomics", "Methylomics"), s_standardize = TRUE, s_k_auto = FALSE, s_k = 8, s_alpha_auto = FALSE, s_alpha = 0.5,
                      s_t_auto = FALSE, s_t = 20, s_cluster_auto = FALSE, s_n_clusters = 2, s_cluster_method = "spectral", s_seed = 1)
    session$setInputs(s_run_btn = 0); session$setInputs(s_run_btn = 1)
    expect_true(isTRUE(snf_state$result$ok))
    pub <- multi_results$integration$snf
    expect_true(isTRUE(pub$result$ok))
    expect_setequal(names(pub$layers), c("Transcriptomics", "Methylomics"))
    expect_equal(nrow(pub$result$W), 56)
    expect_true(any(vapply(arthomix_provenance_records(session), function(r) identical(r$module, "mod_multi_integration_snf"), logical(1))))
  })
})

test_that("SNF Clustering can carry the Integration tab's fused network unchanged, partition it with its own settings, and stress-test it with the same SNF parameters", {
  fx <- ph_fixture()
  layers <- shiny::isolate(fx$multi_dataset$layers)
  snf <- mi_snf_run(layers, list(standardize = TRUE, k_mode = "manual", k = 8, alpha_mode = "manual", alpha = 0.5, t_mode = "manual", t = 20, cluster_mode = "manual", n_clusters = 2, cluster_method = "spectral", seed = 1))
  fx$multi_results$integration <- list(snf = list(result = snf, layers = layers, sample_meta = shiny::isolate(fx$multi_dataset$sample_meta), dataset_label = "fixture", run_at = Sys.time()))
  shiny::testServer(mod_multi_stratification_server, args = list(id = "sc", multi_dataset = fx$multi_dataset, multi_results = fx$multi_results), {
    session$setInputs(data_source = "integration")
    d <- sc_dataset(); expect_true(d$ok); expect_false(is.null(d$carried))
    expect_match(strip_html(output$preproc_ui$html), "Not applicable")
    session$setInputs(blocks = c("Transcriptomics", "Methylomics"), n_clusters = 3, cluster_auto = FALSE, cluster_method = "spectral", seed = 1)
    session$setInputs(run_btn = 0); session$setInputs(run_btn = 1)
    r <- state$result
    expect_true(isTRUE(r$ok))
    expect_equal(r$params$mode, "carried_fused_network")
    expect_equal(r$W, snf$W)
    expect_equal(length(unique(r$clusters)), 3)
    expect_equal(c(r$params$k, r$params$alpha, r$params$t), c(8, 0.5, 20))
    st <- state$stability
    expect_true(isTRUE(st$ok)); expect_gt(st$n_distinct_subsamples, 1)
    rec <- Filter(function(x) identical(x$module, "mod_multi_stratification"), arthomix_provenance_records(session))
    expect_true(length(rec) >= 1)
    expect_true(isTRUE(rec[[length(rec)]]$params$carried_from_integration))
  })
  ## without a published network the option refuses cleanly
  fx2 <- ph_fixture()
  shiny::testServer(mod_multi_stratification_server, args = list(id = "sc", multi_dataset = fx2$multi_dataset, multi_results = fx2$multi_results), {
    session$setInputs(data_source = "integration")
    d <- sc_dataset(); expect_false(d$ok); expect_match(d$error, "Integration")
  })
})

test_that("multi_pipeline_progress() reflects which stages have published results", {
  md <- list(active = TRUE, layers = list(a = 1, b = 2))
  mr <- list(overview = list(harmonization = list(ok = TRUE)), integration = list(snf = list(result = list(ok = TRUE))), biomarker = list(df = data.frame(x = 1)))
  pr <- multi_pipeline_progress(md, mr)
  expect_true(pr[["Dataset Workspace"]]); expect_true(pr[["Cohort Harmonisation"]]); expect_true(pr[["Integration (DIABLO / SNF)"]])
  expect_false(pr[["SNF Clustering"]]); expect_true(pr[["Biomarker Discovery"]]); expect_false(pr[["Gene-CpG Mapping"]])
})

test_that("activating a dataset in the Dataset Workspace writes a provenance record with the matching and preprocessing choices", {
  skip_if_not(MULTI_DATA_AVAILABLE, "multi-omics preloaded data not available in this deployment")
  multi_dataset <- shiny::reactiveValues(); multi_results <- shiny::reactiveValues()
  shiny::testServer(mod_multi_dataset_server, args = list(id = "ds", multi_dataset = multi_dataset, multi_results = multi_results), {
    session$setInputs(dataset_source = "preloaded")
    session$setInputs(preloaded_pick = "ra_antitnf", preloaded_cell = "female_Etanercept")
    session$setInputs(load_preloaded_btn = 1)
    session$setInputs(impute_method = "none", max_sample_missing = 100, max_feature_missing = 100, norm_Transcriptomics = "none", norm_Methylomics = "none")
    session$setInputs(preprocess_btn = 1)
    session$setInputs(active_layers = c("Transcriptomics", "Methylomics"), activate_btn = 1)
    recs <- Filter(function(x) identical(x$module, "mod_multi_dataset"), arthomix_provenance_records(session))
    expect_equal(length(recs), 1)
    expect_equal(recs[[1]]$params$n_matched, 29)
    expect_true(isTRUE(recs[[1]]$params$same_patient_check))
    expect_equal(recs[[1]]$params$normalisation$Transcriptomics, "none")
  })
})
