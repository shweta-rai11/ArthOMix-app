## Module 3 (Multiomics) - SNF Clustering / Patient Stratification
## sub-module, via testServer(): data-source dispatch, raw validation,
## the real per-block preprocessing chain (dynamic block-scoped input IDs),

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_plots.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "multiomics_dataset_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "multiomics_dataset_plots.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_integration_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_integration_plots.R"))
source_from_app_root(file.path("R", "multiomics", "02_Cohort_Harmonization", "cohort_harmonization_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "03_DIABLO_SNF_Integration", "mod_multi_integration.R"))
source_from_app_root(file.path("R", "multiomics", "04_SNF_Clustering", "snf_clustering_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "04_SNF_Clustering", "snf_clustering_plots.R"))
source_from_app_root(file.path("R", "multiomics", "04_SNF_Clustering", "mod_multi_stratification.R"))

sc_live_fixture <- function(n = 20, seed = 1200) {
  set.seed(seed)
  ids <- paste0("S", seq_len(n))
  expr <- matrix(rnorm(n * 10), n, 10, dimnames = list(ids, paste0("g", 1:10)))
  meth <- matrix(rnorm(n * 10), n, 10, dimnames = list(ids, paste0("cg", 1:10)))
  list(ids = ids, expr = expr, meth = meth)
}

test_that("sc_dataset()/sc_val_raw() compute real validation on the Active Multi-Omics Dataset", {
  skip_if_not(MULTI_SNF_LIVE_AVAILABLE, "SNFtool not installed")
  fx <- sc_live_fixture()
  multi_dataset <- shiny::reactiveValues(active = TRUE, source = "upload", layers = list(Transcriptomics = fx$expr, Methylomics = fx$meth), sample_meta = NULL)
  shiny::testServer(mod_multi_stratification_server, args = list(id = "sc", multi_dataset = multi_dataset), {
    session$setInputs(data_source = "active", blocks = c("Transcriptomics", "Methylomics"))
    d <- sc_dataset()
    expect_true(d$ok)
    v <- sc_val_raw()
    expect_true(v$ok)
    expect_equal(v$n_shared, 20L)
  })
})

test_that("sc_ready() applies the real per-block preprocessing chain (missing-value handling + log2 transform) via dynamic block-scoped input IDs", {
  skip_if_not(MULTI_SNF_LIVE_AVAILABLE, "SNFtool not installed")
  fx <- sc_live_fixture()
  expr_with_na <- fx$expr
  expr_with_na[1, 1] <- NA
  multi_dataset <- shiny::reactiveValues(active = TRUE, source = "upload", layers = list(Transcriptomics = expr_with_na, Methylomics = fx$meth), sample_meta = NULL)
  shiny::testServer(mod_multi_stratification_server, args = list(id = "sc", multi_dataset = multi_dataset), {
    session$setInputs(data_source = "active", blocks = c("Transcriptomics", "Methylomics"))
    session$setInputs(
      transform_Transcriptomics = "none", missing_Transcriptomics = "mean", misspct_Transcriptomics = 50,
      filtcrit_Transcriptomics = "none",
      transform_Methylomics = "none", missing_Methylomics = "none", misspct_Methylomics = 50,
      filtcrit_Methylomics = "none"
    )
    r <- sc_ready()
    expect_equal(length(r$errors), 0L)
    expect_false(anyNA(r$layers$Transcriptomics))
    expect_equal(dim(r$layers$Methylomics), dim(fx$meth))
  })
})

test_that("sc_elig() refuses when a selected block's own missing values are left unhandled ('none' with real NAs)", {
  skip_if_not(MULTI_SNF_LIVE_AVAILABLE, "SNFtool not installed")
  fx <- sc_live_fixture()
  expr_with_na <- fx$expr
  expr_with_na[1, 1] <- NA
  multi_dataset <- shiny::reactiveValues(active = TRUE, source = "upload", layers = list(Transcriptomics = expr_with_na, Methylomics = fx$meth), sample_meta = NULL)
  shiny::testServer(mod_multi_stratification_server, args = list(id = "sc", multi_dataset = multi_dataset), {
    session$setInputs(data_source = "active", blocks = c("Transcriptomics", "Methylomics"))
    session$setInputs(
      transform_Transcriptomics = "none", missing_Transcriptomics = "none", misspct_Transcriptomics = 50, filtcrit_Transcriptomics = "none",
      transform_Methylomics = "none", missing_Methylomics = "none", misspct_Methylomics = 50, filtcrit_Methylomics = "none"
    )
    expect_false(sc_elig()$ok)
  })
})

test_that("clicking 'Run SNF Clustering' while ineligible is blocked by validate() before any async dispatch - state$submitted stays FALSE", {
  skip_if_not(MULTI_SNF_LIVE_AVAILABLE, "SNFtool not installed")
  fx <- sc_live_fixture()
  expr_with_na <- fx$expr
  expr_with_na[1, 1] <- NA
  multi_dataset <- shiny::reactiveValues(active = TRUE, source = "upload", layers = list(Transcriptomics = expr_with_na, Methylomics = fx$meth), sample_meta = NULL)
  shiny::testServer(mod_multi_stratification_server, args = list(id = "sc", multi_dataset = multi_dataset), {
    session$setInputs(data_source = "active", blocks = c("Transcriptomics", "Methylomics"))
    session$setInputs(
      transform_Transcriptomics = "none", missing_Transcriptomics = "none", misspct_Transcriptomics = 50, filtcrit_Transcriptomics = "none",
      transform_Methylomics = "none", missing_Methylomics = "none", misspct_Methylomics = 50, filtcrit_Methylomics = "none"
    )
    expect_false(sc_elig()$ok)
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)
    expect_false(isTRUE(state$submitted))
    expect_null(state$result)
  })
})

test_that("clicking 'Run SNF Clustering' while eligible synchronously flips state$submitted to TRUE, snapshotting the real preprocessed layers", {
  skip_if_not(MULTI_SNF_LIVE_AVAILABLE, "SNFtool not installed")
  fx <- sc_live_fixture()
  multi_dataset <- shiny::reactiveValues(active = TRUE, source = "upload", layers = list(Transcriptomics = fx$expr, Methylomics = fx$meth), sample_meta = NULL)
  shiny::testServer(mod_multi_stratification_server, args = list(id = "sc", multi_dataset = multi_dataset), {
    session$setInputs(data_source = "active", blocks = c("Transcriptomics", "Methylomics"))
    session$setInputs(
      transform_Transcriptomics = "none", missing_Transcriptomics = "none", misspct_Transcriptomics = 50, filtcrit_Transcriptomics = "none",
      transform_Methylomics = "none", missing_Methylomics = "none", misspct_Methylomics = 50, filtcrit_Methylomics = "none",
      k_auto = TRUE, alpha_auto = TRUE, t_auto = TRUE, cluster_auto = TRUE
    )
    expect_true(sc_elig()$ok)
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)
    expect_true(isTRUE(state$submitted))
    expect_equal(dim(state$layers_used$Transcriptomics), dim(fx$expr))
  })
})
