## Module 1 (Transcriptomics) - Sex Interaction Analysis: a live limma
## group*sex interaction model, via testServer(). This submodule is
## inherently a "both sexes present" analysis (Case 1 of the sex-
## stratification test matrix); Cases 2-4 (male-only/female-only/no-sex-
## metadata) are covered here as the "controls UI correctly blocks running
## at all" behavior this module is designed to have in those situations.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "11_Sex_Interaction_Analysis", "mod_interaction.R"))

interaction_fixture <- function(n_per_cell = 4, sex_col = TRUE, sexes = c("F", "M"), seed = 80) {
  set.seed(seed)
  n_groups <- 2; n_sexes <- length(sexes)
  n <- n_per_cell * n_groups * n_sexes
  genes <- paste0("GENE", 1:25)
  samples <- paste0("S", 1:n)
  grp <- rep(rep(c("HC", "RA"), each = n_per_cell), times = n_sexes)
  sx  <- rep(sexes, each = n_per_cell * n_groups)
  m <- matrix(rnorm(25 * n, mean = 8, sd = 1.2), 25, n, dimnames = list(genes, samples))
  meta <- data.frame(sample = samples, group = grp, stringsAsFactors = FALSE)
  if (sex_col) meta$sex <- sx
  shiny::reactiveValues(expr = m, meta = meta, source = "interaction test", source_type = "uploaded",
                          is_bundled_reference = FALSE, geo_ids = character(0))
}

## Prime helper for the ignoreInit=TRUE bare-actionButton gotcha (see
## feedback_shiny_testserver_ignoreinit_actionbutton_priming memory).
run_click <- function(session) {
  session$setInputs(run_btn = 0)
  session$setInputs(run_btn = 1)
}

test_that("Case 1 (M+F present): the interaction model fits and reports a coefficient name and significance counts", {
  dataset <- interaction_fixture()
  shiny::testServer(mod_interaction_server, args = list(id = "int", dataset = dataset), {
    session$setInputs(ref_group = "HC", comp_group = "RA", ref_sex = "F", comp_sex = "M", padj_cut = 0.05)
    run_click(session)
    res <- fit_result()
    expect_true(grepl(":", res$coef_name))  ## interaction term name contains the ':' between grp/sx
    df <- sig_table()
    expect_true(all(c("gene", "logFC", "adj.P.Val", "significant") %in% colnames(df)))
    expect_equal(nrow(df), 25L)
  })
})

test_that("Case 2 (male-only data): the controls UI blocks with 'needs a sex column with at least two values', not a false interaction result", {
  dataset <- interaction_fixture(sexes = "M")
  shiny::testServer(mod_interaction_server, args = list(id = "int", dataset = dataset), {
    err <- tryCatch(output$controls, error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Needs a sex column with at least two values", conditionMessage(err)))
  })
})

test_that("Case 3 (female-only data): same guard applies", {
  dataset <- interaction_fixture(sexes = "F")
  shiny::testServer(mod_interaction_server, args = list(id = "int", dataset = dataset), {
    err <- tryCatch(output$controls, error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Needs a sex column with at least two values", conditionMessage(err)))
  })
})

test_that("Case 4 (no sex metadata column at all): the pooled/fallback behavior is the same explicit block, not a silent pooled analysis", {
  dataset <- interaction_fixture(sex_col = FALSE)
  shiny::testServer(mod_interaction_server, args = list(id = "int", dataset = dataset), {
    err <- tryCatch(output$controls, error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Needs a sex column with at least two values", conditionMessage(err)))
  })
})

test_that("fewer than two group values blocks the controls UI", {
  dataset <- interaction_fixture()
  shiny::isolate(dataset$meta$group <- "HC")
  shiny::testServer(mod_interaction_server, args = list(id = "int", dataset = dataset), {
    err <- tryCatch(output$controls, error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Needs at least two group values", conditionMessage(err)))
  })
})

test_that("identical reference/comparison group or sex is rejected", {
  dataset <- interaction_fixture()
  shiny::testServer(mod_interaction_server, args = list(id = "int", dataset = dataset), {
    session$setInputs(ref_group = "HC", comp_group = "HC", ref_sex = "F", comp_sex = "M", padj_cut = 0.05)
    run_click(session)
    err <- tryCatch(fit_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Reference and comparison group must be different", conditionMessage(err)))
  })
})

test_that("fewer than 12 total matching samples across both sexes/groups is rejected", {
  dataset <- interaction_fixture(n_per_cell = 2)  ## 2*2*2 = 8 total, below 12
  shiny::testServer(mod_interaction_server, args = list(id = "int", dataset = dataset), {
    session$setInputs(ref_group = "HC", comp_group = "RA", ref_sex = "F", comp_sex = "M", padj_cut = 0.05)
    run_click(session)
    err <- tryCatch(fit_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Fewer than 12 samples", conditionMessage(err)))
  })
})

test_that("a group-by-sex cell with fewer than 2 samples is rejected even when the total sample count is high enough", {
  dataset <- interaction_fixture(n_per_cell = 4)
  ## Collapse one cell (HC/F) down to a single sample by reassigning the rest to RA/F.
  shiny::isolate({
    idx <- which(dataset$meta$group == "HC" & dataset$meta$sex == "F")
    dataset$meta$group[idx[-1]] <- "RA"
  })
  shiny::testServer(mod_interaction_server, args = list(id = "int", dataset = dataset), {
    session$setInputs(ref_group = "HC", comp_group = "RA", ref_sex = "F", comp_sex = "M", padj_cut = 0.05)
    run_click(session)
    err <- tryCatch(fit_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Each group-by-sex cell needs at least 2 samples", conditionMessage(err)))
  })
})
