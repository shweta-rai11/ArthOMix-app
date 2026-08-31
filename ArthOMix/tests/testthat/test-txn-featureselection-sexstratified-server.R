## Module 1 (Transcriptomics) - Feature Selection's per-sex server flow via
## testServer(), using data_source="expr" (its own self-contained upload,
## bypassing the shared `dataset`/WGCNA/DGE dependencies entirely) and
## gene_source="custom" (a pasted gene list) - the fastest, most isolated
## path through fs_build_sex() to results$featureselection. Manual RF/SVM
## params keep each run's internal model fitting fast.
##
## Note on button clicks: unlike mod_dge.R's fit_result (which needs
## session$setInputs(run_btn = 0) primed before the real "= 1" click - see
## feedback_shiny_testserver_ignoreinit_actionbutton_priming), a single
## plain session$setInputs(run_female_btn = 1) here already fires
## fs_result_female() exactly once (confirmed empirically) - this module's
## results$featureselement_runs is appended to by a SEPARATE
## observeEvent(fs_result_female(), ...) (ignoreInit = FALSE, keyed on the
## eventReactive's own value rather than the raw button), which changes the
## priming arithmetic; adding the DGE-style prime step here double-counts.
## Always verify empirically per module rather than assuming one pattern.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "mod_featureselection.R"))

## sex_n: either a single count (applied to both sexes) or a named
## c(F = , M = ) vector for asymmetric per-sex sample counts. Each sex gets
## 2x its count (split evenly HC/RA).
fs_upload_fixture <- function(sex_n = 6, sex_col = TRUE, sexes = c("F", "M"), seed = 111) {
  set.seed(seed)
  if (length(sex_n) == 1) sex_n <- setNames(rep(sex_n, length(sexes)), sexes)
  genes <- paste0("GENE", 1:15)
  grp <- unlist(lapply(sexes, function(s) rep(c("HC", "RA"), each = sex_n[[s]])))
  sx <- unlist(lapply(sexes, function(s) rep(s, sex_n[[s]] * 2)))
  n <- length(grp)
  samples <- paste0("S", 1:n)
  expr <- matrix(rnorm(15 * n, mean = 8, sd = 1.2), 15, n, dimnames = list(genes, samples))
  ## Real separating signal on GENE1/GENE2 for the RA group, present in both sexes.
  expr[1:2, grp == "RA"] <- expr[1:2, grp == "RA"] + 3
  meta <- data.frame(sample = samples, group = grp, stringsAsFactors = FALSE)
  if (sex_col) meta$sex <- sx else meta$no_sex_signal <- "unknown"  ## a genuinely single-valued stand-in column

  ## Not withr::local_tempdir(): called from this helper (not directly inside
  ## a test_that() block), its cleanup would defer to THIS function's own
  ## frame and delete the directory the moment fs_upload_fixture() returns,
  ## before the caller ever reads the files back.
  dir <- tempfile("fs_upload_fixture_")
  dir.create(dir)
  expr_path <- file.path(dir, "expr.csv")
  meta_path <- file.path(dir, "meta.csv")
  write.csv(data.frame(gene = rownames(expr), expr, check.names = FALSE), expr_path, row.names = FALSE)
  write.csv(meta, meta_path, row.names = FALSE)
  list(expr_path = expr_path, meta_path = meta_path, genes = genes)
}

fs_common_inputs <- function(session, fx, map_sex = "sex") {
  session$setInputs(data_source = "expr")
  session$setInputs(expr_file = fx_mkfile(fx$expr_path))
  session$setInputs(meta_file = fx_mkfile(fx$meta_path))
  session$setInputs(map_id = "sample", map_group = "group", map_sex = map_sex)
  session$setInputs(gene_source = "custom", gene_list = paste(fx$genes, collapse = ","))
  session$setInputs(ref_group = "HC", comp_group = "RA")
  session$setInputs(class_weight_mode = "equal",
                      rf_mtry_mode = "manual", rf_mtry_manual = 3, rf_ntree = 100,
                      svm_cost_mode = "manual", svm_cost_manual = 1)
}

test_that("Case 1 (M+F present): running Female then Male populates results$featureselection for both sexes with recoverable signal", {
  fx <- fs_upload_fixture()
  dataset <- shiny::reactiveValues(expr = NULL, meta = NULL, source_type = "preloaded")
  results <- shiny::reactiveValues()
  shiny::testServer(mod_featureselection_server, args = list(id = "fs", dataset = dataset, results = results), {
    fs_common_inputs(session, fx)
    session$setInputs(run_female_btn = 1)
    session$setInputs(run_male_btn = 1)

    r_f <- fs_result_female()
    r_m <- fs_result_male()
    expect_true(any(c("GENE1", "GENE2") %in% r_f$rf_genes))
    expect_true(any(c("GENE1", "GENE2") %in% r_m$rf_genes))

    expect_false(is.null(results$featureselection$female))
    expect_false(is.null(results$featureselection$male))
    expect_length(results$featureselection_runs, 2)
  })
})

test_that("Case 4 (no sex metadata column at all): sex_levels() blocks with an explicit message, but the sex-less Pooled run still works", {
  fx <- fs_upload_fixture(sex_col = FALSE)
  dataset <- shiny::reactiveValues(expr = NULL, meta = NULL, source_type = "preloaded")
  results <- shiny::reactiveValues()
  shiny::testServer(mod_featureselection_server, args = list(id = "fs", dataset = dataset, results = results), {
    session$setInputs(data_source = "expr")
    session$setInputs(expr_file = fx_mkfile(fx$expr_path))
    session$setInputs(meta_file = fx_mkfile(fx$meta_path))
    ## "no_sex_signal" is a genuinely single-valued stand-in for "no usable
    ## sex column" - mapping to "group" instead would trivially pass
    ## sex_levels()'s >=2-distinct-values check (HC/RA), which isn't a real
    ## sex signal but numerically looks like one.
    session$setInputs(map_id = "sample", map_group = "group", map_sex = "no_sex_signal")
    err <- tryCatch(sex_levels(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("sex column needs at least two distinct values", conditionMessage(err)))

    session$setInputs(gene_source = "custom", gene_list = paste(fx$genes, collapse = ","))
    session$setInputs(ref_group = "HC", comp_group = "RA")
    session$setInputs(class_weight_mode = "equal",
                        rf_mtry_mode = "manual", rf_mtry_manual = 3, rf_ntree = 100,
                        svm_cost_mode = "manual", svm_cost_manual = 1)
    session$setInputs(run_pooled_btn = 1)
    r_p <- fs_result_pooled()
    expect_true(any(c("GENE1", "GENE2") %in% r_p$rf_genes))
    expect_false(is.null(results$featureselection$pooled))
  })
})

test_that("fewer than 10 samples for a given sex is rejected before any model is fit, even with >=20 total samples", {
  ## Female stratum: 4 samples (below the 10-sample floor). Male stratum: 20
  ## samples, so the TOTAL (24) clears source_expr_meta()'s own >=20-overlap
  ## gate - isolating the per-sex floor from that earlier, unrelated gate.
  fx <- fs_upload_fixture(sex_n = c(F = 2, M = 10))
  dataset <- shiny::reactiveValues(expr = NULL, meta = NULL, source_type = "preloaded")
  results <- shiny::reactiveValues()
  shiny::testServer(mod_featureselection_server, args = list(id = "fs", dataset = dataset, results = results), {
    fs_common_inputs(session, fx)
    session$setInputs(run_female_btn = 1)
    err <- tryCatch(fs_result_female(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Fewer than 10 female samples", conditionMessage(err)))
  })
})

test_that("fewer than 3 candidate genes present in the expression matrix is rejected", {
  fx <- fs_upload_fixture()
  dataset <- shiny::reactiveValues(expr = NULL, meta = NULL, source_type = "preloaded")
  results <- shiny::reactiveValues()
  shiny::testServer(mod_featureselection_server, args = list(id = "fs", dataset = dataset, results = results), {
    fs_common_inputs(session, fx)
    session$setInputs(gene_list = "NOTAGENE1,NOTAGENE2")  ## none present in the matrix
    session$setInputs(run_female_btn = 1)
    err <- tryCatch(fs_result_female(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("Fewer than 3 female candidate genes", conditionMessage(err)))
  })
})
