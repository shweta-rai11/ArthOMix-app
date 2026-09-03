## Module 2 (Methylomics) - Sex Interaction Analysis: a live limma
## group*sex interaction model fit on M-values, via testServer(). Structural
## clone of tests/testthat/test-txn-interaction-server.R's 8 cases, adapted

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "functions", "qc.R"))
source_from_app_root(file.path("R", "methylomics", "functions", "annotation.R"))
source_from_app_root(file.path("R", "methylomics", "functions", "normalization.R"))
source_from_app_root(file.path("R", "methylomics", "05_Differential_Methylation_Position", "mod_methyl_dmp.R"))
source_from_app_root(file.path("R", "methylomics", "07_Sex_Interaction_Analysis", "mod_methyl_interaction.R"))

methyl_interaction_fixture <- function(n_per_cell = 4, sex_col = TRUE, sexes = c("F", "M"), seed = 380) {
  set.seed(seed)
  n_groups <- 2; n_sexes <- length(sexes)
  n <- n_per_cell * n_groups * n_sexes
  cpgs <- paste0("cg", 10000000 + 1:30)
  samples <- paste0("S", 1:n)
  grp <- rep(rep(c("HC", "RA"), each = n_per_cell), times = n_sexes)
  sx  <- rep(sexes, each = n_per_cell * n_groups)
  m <- matrix(stats::runif(30 * n, 0.2, 0.8), 30, n, dimnames = list(cpgs, samples))
  sheet <- data.frame(sample = samples, group = grp, stringsAsFactors = FALSE)
  if (sex_col) sheet$sex <- sx
  shiny::reactiveValues(beta = m, sample_sheet = sheet, input_scale = "beta", array_type = "EPIC",
                          rg_set = NULL, mset = NULL, detp = NULL, beadcount = NULL,
                          preloaded = FALSE, source_type = "uploaded", source = "methyl interaction test")
}

run_click <- function(session) {
  session$setInputs(run_btn = 0)
  session$setInputs(run_btn = 1)
}

test_that("Case 1 (M+F present): the interaction model fits on M-values and reports a coefficient name, significance counts, and a beta-scale delta", {
  dataset <- methyl_interaction_fixture()
  methyl_results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_interaction_server, args = list(id = "int", methyl_dataset = dataset, methyl_results = methyl_results), {
    session$setInputs(group_col = "group", ref_group = "HC", comp_group = "RA", ref_sex = "F", comp_sex = "M", padj_cut = 0.05)
    run_click(session)
    res <- fit_result()
    expect_true(grepl(":", res$coef_name))
    df <- sig_table()
    expect_true(all(c("cpg", "logFC", "adj.P.Val", "significant") %in% colnames(df)))
    expect_equal(nrow(df), 30L)

    expect_false(isTRUE(res$is_m_scale))
    expect_true("dbeta_interaction" %in% colnames(df))
    expect_true(all(abs(df$dbeta_interaction) <= 1, na.rm = TRUE))

    expect_false(is.null(methyl_results$interaction))
    expect_true(grepl(":", methyl_results$interaction$coef_name))
    expect_equal(methyl_results$interaction$n_tested, 30L)
  })
})

test_that("Case 2 (male-only data): the controls UI blocks with 'Needs a sex column with at least two values', not a false interaction result", {
  dataset <- methyl_interaction_fixture(sexes = "M")
  shiny::testServer(mod_methyl_interaction_server, args = list(id = "int", methyl_dataset = dataset, methyl_results = shiny::reactiveValues()), {
    err <- tryCatch(output$controls, error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Needs a sex column with at least two values", conditionMessage(err)))
  })
})

test_that("Case 3 (female-only data): same guard applies", {
  dataset <- methyl_interaction_fixture(sexes = "F")
  shiny::testServer(mod_methyl_interaction_server, args = list(id = "int", methyl_dataset = dataset, methyl_results = shiny::reactiveValues()), {
    err <- tryCatch(output$controls, error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Needs a sex column with at least two values", conditionMessage(err)))
  })
})

test_that("Case 4 (no sex metadata column at all): the pooled/fallback behavior is the same explicit block, not a silent pooled analysis", {
  dataset <- methyl_interaction_fixture(sex_col = FALSE)
  shiny::testServer(mod_methyl_interaction_server, args = list(id = "int", methyl_dataset = dataset, methyl_results = shiny::reactiveValues()), {
    err <- tryCatch(output$controls, error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Needs a sex column with at least two values", conditionMessage(err)))
  })
})

test_that("fewer than two group values blocks the group-level UI", {
  dataset <- methyl_interaction_fixture()
  shiny::isolate(dataset$sample_sheet$group <- "HC")
  shiny::testServer(mod_methyl_interaction_server, args = list(id = "int", methyl_dataset = dataset, methyl_results = shiny::reactiveValues()), {
    session$setInputs(group_col = "group")
    err <- tryCatch(output$group_level_ui, error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Needs at least two group values", conditionMessage(err)))
  })
})

test_that("identical reference/comparison group or sex is rejected", {
  dataset <- methyl_interaction_fixture()
  shiny::testServer(mod_methyl_interaction_server, args = list(id = "int", methyl_dataset = dataset, methyl_results = shiny::reactiveValues()), {
    session$setInputs(group_col = "group", ref_group = "HC", comp_group = "HC", ref_sex = "F", comp_sex = "M", padj_cut = 0.05)
    run_click(session)
    err <- tryCatch(fit_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Reference and comparison group must be different", conditionMessage(err)))
  })
})

test_that("fewer than 12 total matching samples across both sexes/groups is rejected", {
  dataset <- methyl_interaction_fixture(n_per_cell = 2)
  shiny::testServer(mod_methyl_interaction_server, args = list(id = "int", methyl_dataset = dataset, methyl_results = shiny::reactiveValues()), {
    session$setInputs(group_col = "group", ref_group = "HC", comp_group = "RA", ref_sex = "F", comp_sex = "M", padj_cut = 0.05)
    run_click(session)
    err <- tryCatch(fit_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Fewer than 12 samples", conditionMessage(err)))
  })
})

test_that("a group-by-sex cell with fewer than 2 samples is rejected even when the total sample count is high enough", {
  dataset <- methyl_interaction_fixture(n_per_cell = 4)
  shiny::isolate({
    idx <- which(dataset$sample_sheet$group == "HC" & dataset$sample_sheet$sex == "F")
    dataset$sample_sheet$group[idx[-1]] <- "RA"
  })
  shiny::testServer(mod_methyl_interaction_server, args = list(id = "int", methyl_dataset = dataset, methyl_results = shiny::reactiveValues()), {
    session$setInputs(group_col = "group", ref_group = "HC", comp_group = "RA", ref_sex = "F", comp_sex = "M", padj_cut = 0.05)
    run_click(session)
    err <- tryCatch(fit_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Each group-by-sex cell needs at least 2 samples", conditionMessage(err)))
  })
})
