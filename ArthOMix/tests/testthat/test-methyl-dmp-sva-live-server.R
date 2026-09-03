## Module 2 (Methylomics) - SVA tab's live engine (sva::sva() + bacon::bacon()
## on top of the shared live-model plumbing), via testServer(). Verifies this
## now works for Upload and GEO-fetched data (methyl_dataset$preloaded ==

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "functions", "qc.R"))
source_from_app_root(file.path("R", "methylomics", "functions", "annotation.R"))
source_from_app_root(file.path("R", "methylomics", "functions", "normalization.R"))
source_from_app_root(file.path("R", "provenance.R"))
source_from_app_root(file.path("R", "methylomics", "05_Differential_Methylation_Position", "mod_methyl_dmp.R"))

sva_fixture_dataset <- function(n_per_group = 20, n_probes = 300, seed = 411, source_type = "uploaded") {
  set.seed(seed)
  n <- n_per_group * 2
  m <- matrix(runif(n_probes * n, 0.2, 0.8), n_probes, n,
              dimnames = list(paste0("cg", 20000000 + seq_len(n_probes)), paste0("S", seq_len(n))))
  m[1:5, (n_per_group + 1):n] <- pmin(m[1:5, (n_per_group + 1):n] + 0.35, 0.99)
  batch <- rep(c("plate1", "plate2"), length.out = n)[sample.int(n)]
  batch_shift <- ifelse(batch == "plate2", 0.15, 0)
  m[6:150, ] <- pmin(pmax(sweep(m[6:150, ], 2, batch_shift, `+`), 0.01), 0.99)

  sheet <- data.frame(sample = colnames(m), group = rep(c("HC", "RA"), each = n_per_group),
                        sex = rep(c("F", "M"), length.out = n), batch = batch, stringsAsFactors = FALSE)
  shiny::reactiveValues(beta = m, sample_sheet = sheet, input_scale = "beta", array_type = "EPIC",
                          rg_set = NULL, mset = NULL, detp = NULL, beadcount = NULL,
                          preloaded = FALSE, source_type = source_type, source = "sva-live test")
}

test_that("the live SVA-adjusted engine runs on uploaded data and produces a well-formed, bacon-corrected result table", {
  methyl_dataset <- sva_fixture_dataset(source_type = "uploaded")
  methyl_results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_dmp_server, args = list(id = "dmp", methyl_dataset = methyl_dataset, methyl_results = methyl_results), {
    session$setInputs(svalive_group_col = "group", svalive_ref = "HC", svalive_comp = "RA", svalive_sex = "__all__",
                        svalive_min_valid_pct = 80, svalive_min_variance = 0, svalive_snp_filter = FALSE,
                        svalive_covariates = character(0))
    session$setInputs(svalive_run_btn = 1)

    r <- svalive_result()
    expect_true(all(c("cpg", "t", "p_raw", "p_bacon", "fdr", "dbeta", "direction") %in% colnames(r$df)))
    expect_true(all(r$df$p_raw >= 0 & r$df$p_raw <= 1, na.rm = TRUE))
    expect_true(all(r$df$p_bacon >= 0 & r$df$p_bacon <= 1, na.rm = TRUE))
    expect_true(all(r$df$fdr >= r$df$p_bacon - 1e-9, na.rm = TRUE))
    expect_equal(r$n_ref, 20L)
    expect_equal(r$n_comp, 20L)
    expect_true(is.numeric(r$n_sv) && r$n_sv >= 0)
    expect_true(is.finite(r$lambda_bacon) && r$lambda_bacon > 0)
    expect_true(is.finite(r$bias_bacon))

    top10 <- r$df$cpg[order(r$df$fdr)][1:10]
    expect_true(sum(top10 %in% paste0("cg", 20000001:20000005)) >= 2)
  })
})

test_that("the live SVA-adjusted engine also runs for a live GEO-fetched dataset (not just Upload)", {
  methyl_dataset <- sva_fixture_dataset(source_type = "geo", seed = 99)
  methyl_results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_dmp_server, args = list(id = "dmp", methyl_dataset = methyl_dataset, methyl_results = methyl_results), {
    session$setInputs(svalive_group_col = "group", svalive_ref = "HC", svalive_comp = "RA", svalive_sex = "__all__",
                        svalive_min_valid_pct = 80, svalive_min_variance = 0, svalive_snp_filter = FALSE,
                        svalive_covariates = character(0))
    session$setInputs(svalive_run_btn = 1)

    r <- svalive_result()
    expect_true(all(c("cpg", "t", "p_raw", "p_bacon", "fdr", "dbeta") %in% colnames(r$df)))
    expect_equal(r$n_ref, 20L)
    expect_equal(r$n_comp, 20L)
  })
})

test_that("identical reference and comparison groups are rejected by the live SVA engine", {
  methyl_dataset <- sva_fixture_dataset()
  methyl_results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_dmp_server, args = list(id = "dmp", methyl_dataset = methyl_dataset, methyl_results = methyl_results), {
    session$setInputs(svalive_group_col = "group", svalive_ref = "HC", svalive_comp = "HC", svalive_sex = "__all__")
    session$setInputs(svalive_run_btn = 1)
    err <- tryCatch(svalive_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("must be different", conditionMessage(err)))
  })
})

test_that("fewer than 3 samples in one group is rejected by the live SVA engine", {
  methyl_dataset <- sva_fixture_dataset(n_per_group = 20)
  shiny::isolate(methyl_dataset$sample_sheet$group <- c(rep("HC", 38), rep("RA", 2)))
  methyl_results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_dmp_server, args = list(id = "dmp", methyl_dataset = methyl_dataset, methyl_results = methyl_results), {
    session$setInputs(svalive_group_col = "group", svalive_ref = "HC", svalive_comp = "RA", svalive_sex = "__all__")
    session$setInputs(svalive_run_btn = 1)
    err <- tryCatch(svalive_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Each group needs at least 3 samples", conditionMessage(err)))
  })
})

test_that("mod_methyl_sva_fit() estimates surrogate variables and returns a full-rank design", {
  set.seed(7)
  n <- 40
  grp <- factor(rep(c("HC", "RA"), each = 20), levels = c("HC", "RA"))
  m <- matrix(rnorm(500 * n), 500, n)
  batch <- rep(c(0, 1), length.out = n)[sample.int(n)]
  m[1:100, ] <- sweep(m[1:100, ], 2, batch * 2, `+`)

  out <- mod_methyl_sva_fit(m, grp, cov_df = NULL, sv_topn = 500)
  expect_true(is.matrix(out$design))
  expect_equal(nrow(out$design), n)
  expect_equal(qr(out$design)$rank, ncol(out$design))
  expect_true(out$n_sv >= 0)
  expect_equal(ncol(out$design), 2 + out$n_sv)
})
