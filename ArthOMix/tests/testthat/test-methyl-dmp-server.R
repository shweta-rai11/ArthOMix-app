## Module 2 (Methylomics) - DMP tab's live limma differential-methylation
## fit, via testServer(), verifying the scientific/data-contract (cpg/t/
## p_raw/fdr/dbeta columns, valid ranges) and key validation gates.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "functions", "qc.R"))
source_from_app_root(file.path("R", "methylomics", "functions", "annotation.R"))
source_from_app_root(file.path("R", "methylomics", "functions", "normalization.R"))
source_from_app_root(file.path("R", "provenance.R"))
source_from_app_root(file.path("R", "methylomics", "05_Differential_Methylation_Position", "mod_methyl_dmp.R"))

dmp_fixture_dataset <- function(n_per_group = 10, seed = 260) {
  set.seed(seed)
  n <- n_per_group * 2
  m <- matrix(runif(60 * n, 0.2, 0.8), 60, n, dimnames = list(paste0("cg", 10000000 + 1:60), paste0("S", 1:n)))
  m[1:5, (n_per_group + 1):n] <- pmin(m[1:5, (n_per_group + 1):n] + 0.3, 0.99)
  sheet <- data.frame(sample = colnames(m), group = rep(c("HC", "RA"), each = n_per_group),
                        sex = rep(c("F", "M"), length.out = n), stringsAsFactors = FALSE)
  shiny::reactiveValues(beta = m, sample_sheet = sheet, input_scale = "beta", array_type = "EPIC",
                          rg_set = NULL, mset = NULL, detp = NULL, beadcount = NULL,
                          preloaded = FALSE, source_type = "uploaded", source = "dmp test")
}

test_that("the live DMP fit recovers real signal and produces a well-formed, scientifically valid result table", {
  methyl_dataset <- dmp_fixture_dataset()
  methyl_results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_dmp_server, args = list(id = "dmp", methyl_dataset = methyl_dataset, methyl_results = methyl_results), {
    session$setInputs(live_group_col = "group", live_ref = "HC", live_comp = "RA", live_sex = "__all__",
                        live_min_valid_pct = 80, live_min_variance = 0, live_snp_filter = FALSE, live_covariates = character(0))
    session$setInputs(live_run_btn = 1)

    r <- live_result()
    expect_true(all(c("cpg", "t", "p_raw", "fdr", "dbeta", "direction") %in% colnames(r$df)))
    expect_true(all(r$df$p_raw >= 0 & r$df$p_raw <= 1, na.rm = TRUE))
    expect_true(all(r$df$fdr >= r$df$p_raw - 1e-9, na.rm = TRUE))
    expect_equal(r$n_ref, 10L)
    expect_equal(r$n_comp, 10L)
    top5 <- r$df$cpg[order(r$df$p_raw)][1:5]
    expect_true(sum(top5 %in% paste0("cg", 10000001:10000005)) >= 3)

    expect_false(is.null(methyl_results$dmp))
    expect_true(grepl("RA vs HC", methyl_results$dmp$comparison))
  })
})

test_that("identical reference and comparison groups are rejected", {
  methyl_dataset <- dmp_fixture_dataset()
  methyl_results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_dmp_server, args = list(id = "dmp", methyl_dataset = methyl_dataset, methyl_results = methyl_results), {
    session$setInputs(live_group_col = "group", live_ref = "HC", live_comp = "HC", live_sex = "__all__")
    session$setInputs(live_run_btn = 1)
    err <- tryCatch(live_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("must be different", conditionMessage(err)))
  })
})

test_that("fewer than 3 samples in one group is rejected even with enough total samples", {
  methyl_dataset <- dmp_fixture_dataset(n_per_group = 10)
  shiny::isolate(methyl_dataset$sample_sheet$group <- c(rep("HC", 18), rep("RA", 2)))
  methyl_results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_dmp_server, args = list(id = "dmp", methyl_dataset = methyl_dataset, methyl_results = methyl_results), {
    session$setInputs(live_group_col = "group", live_ref = "HC", live_comp = "RA", live_sex = "__all__")
    session$setInputs(live_run_btn = 1)
    err <- tryCatch(live_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Each group needs at least 3 samples", conditionMessage(err)))
  })
})

test_that("restricting to a single sex with too few remaining samples is rejected", {
  methyl_dataset <- dmp_fixture_dataset(n_per_group = 10)
  shiny::isolate(methyl_dataset$sample_sheet$sex <- c(rep("M", 8), rep("F", 2), rep("M", 8), rep("F", 2)))
  methyl_results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_dmp_server, args = list(id = "dmp", methyl_dataset = methyl_dataset, methyl_results = methyl_results), {
    session$setInputs(live_group_col = "group", live_ref = "HC", live_comp = "RA", live_sex = "F")
    session$setInputs(live_run_btn = 1)
    err <- tryCatch(live_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Fewer than 6 samples remain after restricting to sex", conditionMessage(err)))
  })
})
