## Module 2 (Methylomics) - QC tab's Sample QC / Probe QC run reactives via
## testServer(), wiring qc.R's already-unit-tested filter/scoring functions
## against a real methyl_dataset-shaped fixture.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "functions", "qc.R"))
source_from_app_root(file.path("R", "methylomics", "functions", "annotation.R"))
source_from_app_root(file.path("R", "methylomics", "functions", "idat_metrics.R"))
source_from_app_root(file.path("R", "methylomics", "02_Quality_Control", "mod_methyl_qc.R"))

methyl_qc_fixture <- function(n_probes = 30, n_samples = 12, seed = 240) {
  set.seed(seed)
  mat <- matrix(runif(n_probes * n_samples, 0, 1), n_probes, n_samples,
                 dimnames = list(paste0("cg", 10000000 + 1:n_probes), paste0("S", 1:n_samples)))
  sheet <- data.frame(sample = colnames(mat), group = rep(c("HC", "RA"), n_samples / 2),
                        sex = rep(c("F", "M"), each = n_samples / 2), stringsAsFactors = FALSE)
  shiny::reactiveValues(beta = mat, sample_sheet = sheet, input_scale = "beta", array_type = "EPIC",
                          rg_set = NULL, mset = NULL, detp = NULL, beadcount = NULL,
                          preloaded = FALSE, source_type = "uploaded")
}

test_that("Sample QC computes a call-rate table for every sample in the active stratum, flagging low call rate correctly", {
  methyl_dataset <- methyl_qc_fixture()
  shiny::isolate(methyl_dataset$beta[1:20, 1] <- NA)
  methyl_results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_qc_server, args = list(id = "qc", methyl_dataset = methyl_dataset, methyl_results = methyl_results), {
    session$setInputs(call_rate_min = 0.95)
    session$setInputs(run_sample_qc_btn = 1)

    r <- sample_qc_result()
    expect_equal(nrow(r$sample_qc), 12L)
    expect_true(r$sample_qc$call_rate_flag[r$sample_qc$sample == "S1"])
    expect_false(any(r$sample_qc$call_rate_flag[r$sample_qc$sample != "S1"]))
  })
})

test_that("current_subgroup() restricts Sample QC to one sex stratum when selected", {
  methyl_dataset <- methyl_qc_fixture()
  methyl_results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_qc_server, args = list(id = "qc", methyl_dataset = methyl_dataset, methyl_results = methyl_results), {
    session$setInputs(live_group_col = "sex", live_stratum = "F")
    sg <- current_subgroup()
    expect_equal(ncol(sg$mat), 6L)
    expect_true(all(sg$included %in% colnames(shiny::isolate(methyl_dataset$beta))[1:6]))
  })
})

test_that("Probe QC's missing-value filter removes probes exceeding the missingness threshold from the active stratum", {
  methyl_dataset <- methyl_qc_fixture()
  shiny::isolate(methyl_dataset$beta[1, ] <- NA)
  methyl_results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_qc_server, args = list(id = "qc", methyl_dataset = methyl_dataset, methyl_results = methyl_results), {
    session$setInputs(f_missing = TRUE, missing_max = 0, f_detp = FALSE, f_beadcount = FALSE,
                        f_snp = FALSE, f_noncpg = FALSE, sexchr_mode = "keep",
                        f_crossreactive = FALSE, f_maf = FALSE, f_variance = FALSE, f_sd = FALSE, f_meanrange = FALSE)
    session$setInputs(run_probe_qc_btn = 1)

    r <- probe_qc_result()
    expect_equal(nrow(r$mat), 30L)
    expect_false("cg10000001" %in% rownames(r$filtered))
    expect_equal(nrow(r$filtered), 29L)
  })
})
