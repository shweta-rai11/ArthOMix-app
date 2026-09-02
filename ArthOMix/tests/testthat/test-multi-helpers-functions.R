## Module 3 (Multiomics) - multiomics_helpers.R's pure functions: the shared
## table/fit loaders, the six analysis-cell lookup, the sex-normalization/
## filtering primitives every live sub-module (Concordance's
## mcc_sex_candidates/mcc_sex_groups, Biomarker/Integration's ss_sex_col,
## multiomics_sexstratified_engine.R's mss_run_stratified()) delegates to
## rather than reimplementing, plus the Overview tab's QC scorecard/summary
## table and the Results Summary report/package-versions helpers.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_helpers.R"))

## ---- multi_read_table() / multi_read_registry_table() ------------------------

test_that("multi_read_table() reads a real precomputed table off disk", {
  skip_if_not(MULTI_DATA_AVAILABLE, "multi-omics preloaded data not available in this deployment")
  out <- multi_read_table(MULTI_TABLE_REGISTRY[["RNA-seq QC summary"]])
  expect_true(out$ok)
  expect_true(is.data.frame(out$df))
  expect_gt(nrow(out$df), 0)
  expect_null(out$error)
})

test_that("multi_read_table() fails soft (never errors) for a missing or NULL path", {
  out_missing <- multi_read_table(file.path(tempdir(), "does-not-exist-xyz.csv"))
  expect_false(out_missing$ok)
  expect_null(out_missing$df)
  expect_true(nzchar(out_missing$error))

  out_null <- multi_read_table(NULL)
  expect_false(out_null$ok)
  expect_true(grepl("no path", out_null$error))
})

test_that("multi_read_table() fails soft on an empty (header-only) CSV", {
  path <- tempfile(fileext = ".csv")
  writeLines("a,b,c", path)
  ## multi_read_table() only checks MULTI_DATA_AVAILABLE, not whether path is
  ## under MULTI_DATA_ROOT, so an arbitrary existing empty file exercises the
  ## nrow==0 branch directly regardless of deployment data availability.
  skip_if_not(MULTI_DATA_AVAILABLE, "multi-omics preloaded data not available in this deployment")
  out <- multi_read_table(path)
  expect_false(out$ok)
  expect_true(grepl("no rows", out$error))
})

test_that("multi_read_registry_table() resolves a MULTI_TABLE_REGISTRY label before reading", {
  skip_if_not(MULTI_DATA_AVAILABLE, "multi-omics preloaded data not available in this deployment")
  out <- multi_read_registry_table("RNA-seq QC summary")
  expect_true(out$ok)
  out_bad <- multi_read_registry_table("Not a real registry label")
  expect_false(out_bad$ok)
})

## ---- MULTI_CELLS / multi_cell_by_key() ----------------------------------------

test_that("MULTI_CELLS defines exactly the six documented sex x drug/response cohorts", {
  keys <- vapply(MULTI_CELLS, `[[`, character(1), "key")
  expect_setequal(keys, c("female_Adalimumab", "male_Adalimumab", "female_Etanercept",
                            "male_Etanercept", "female_response", "male_response"))
  ## has_snf is FALSE for both Adalimumab cells (SNF was never run for
  ## Adalimumab - Table22 contains only Etanercept rows) and for the two
  ## drug-pooled response cells (SNF was not run drug-pooled), TRUE only for
  ## the two Etanercept drug x sex cells.
  by_key <- setNames(MULTI_CELLS, keys)
  expect_false(by_key$female_Adalimumab$has_snf)
  expect_false(by_key$male_Adalimumab$has_snf)
  expect_true(by_key$female_Etanercept$has_snf)
  expect_true(by_key$male_Etanercept$has_snf)
  expect_false(by_key$female_response$has_snf)
  expect_false(by_key$male_response$has_snf)
  ## Drug-pooled cells carry no drug value.
  expect_true(is.na(by_key$female_response$drug))
  expect_true(is.na(by_key$male_response$drug))
})

test_that("multi_cell_by_key() returns the matching cell, or NULL for an unknown key", {
  cell <- multi_cell_by_key("female_Etanercept")
  expect_equal(cell$sex, "female")
  expect_equal(cell$drug, "Etanercept")
  expect_null(multi_cell_by_key("not_a_real_key"))
})

## ---- multi_norm_sex() ----------------------------------------------------------

test_that("multi_norm_sex() normalizes every recognized female/male spelling, case-insensitively", {
  expect_equal(multi_norm_sex(c("Female", "f", "FEMALE", "female")), rep("female", 4))
  expect_equal(multi_norm_sex(c("Male", "m", "MALE", "male")), rep("male", 4))
})

test_that("multi_norm_sex() passes through unrecognized values unchanged (not fabricated as female/male)", {
  expect_equal(multi_norm_sex(c("unknown", "", NA_character_)), c("unknown", "", NA_character_))
})

## ---- multi_filter_cell() -------------------------------------------------------

test_that("multi_filter_cell() subsets by both sex and drug when both columns are present", {
  df <- data.frame(id = 1:4, sex = c("Female", "female", "Male", "male"),
                     drug = c("Adalimumab", "Etanercept", "Adalimumab", "Etanercept"))
  out <- multi_filter_cell(df, sex = "female", drug = "Etanercept")
  expect_equal(out$id, 2L)
})

test_that("multi_filter_cell() subsets by sex alone when drug is NA (drug-pooled cells)", {
  df <- data.frame(id = 1:4, sex = c("female", "female", "male", "male"))
  out <- multi_filter_cell(df, sex = "female", drug = NA_character_)
  expect_setequal(out$id, c(1L, 2L))
})

test_that("multi_filter_cell() does not drop rows when the sex or drug column is absent from df", {
  df <- data.frame(id = 1:3, value = c(10, 20, 30))  ## no sex/drug columns at all
  out <- multi_filter_cell(df, sex = "female", drug = "Etanercept")
  expect_equal(nrow(out), 3L)  ## nothing to filter on -> passthrough, not an empty result
})

## ---- multi_sex_candidates() / multi_sex_groups() ------------------------------

test_that("multi_sex_candidates() finds sex/gender-named columns via the shared regex, case-insensitively", {
  meta <- data.frame(Sample_ID = 1:3, Sex = c("F", "M", "F"), Gender = c("Female", "Male", "Female"),
                       age = c(40, 50, 60))
  cands <- multi_sex_candidates(meta)
  expect_setequal(cands, c("Sex", "Gender"))
})

test_that("multi_sex_candidates() returns character(0) when no sex/gender column exists", {
  meta <- data.frame(id = 1:3, age = c(40, 50, 60), batch = c("A", "B", "A"))
  expect_equal(multi_sex_candidates(meta), character(0))
})

test_that("multi_sex_groups() splits sample IDs (keyed by rownames) into named groups by the chosen sex column's values", {
  meta <- data.frame(sex = c("F", "F", "M", "M", "F", "M"), row.names = paste0("S", 1:6))
  groups <- multi_sex_groups(meta, sex_col = "sex", sample_ids = rownames(meta))
  expect_setequal(names(groups), c("F", "M"))
  expect_setequal(groups[["F"]], c("S1", "S2", "S5"))
  expect_setequal(groups[["M"]], c("S3", "S4", "S6"))
})

test_that("multi_sex_groups() restricts to only the requested sample_ids, ignoring unrelated metadata rows", {
  meta <- data.frame(sex = c("F", "F", "M", "M", "F", "M"), row.names = paste0("S", 1:6))
  groups <- multi_sex_groups(meta, sex_col = "sex", sample_ids = c("S1", "S3", "S4"))
  expect_setequal(unlist(groups), c("S1", "S3", "S4"))
})

test_that("multi_sex_groups() returns NULL when sex_col is absent or sample_meta is NULL", {
  meta <- data.frame(sex = c("F", "M"), row.names = c("S1", "S2"))
  expect_null(multi_sex_groups(meta, sex_col = "not_a_col", sample_ids = c("S1", "S2")))
  expect_null(multi_sex_groups(NULL, sex_col = "sex", sample_ids = c("S1", "S2")))
})

## ---- multi_diablo_fit() / multi_diablo_variance_df() ---------------------------

test_that("multi_diablo_fit() loads a real saved mixOmics block.splsda fit for a known cell", {
  skip_if_not(MULTI_DATA_AVAILABLE, "multi-omics preloaded data not available in this deployment")
  cell <- multi_cell_by_key("female_Etanercept")
  skip_if_not(file.exists(MULTI_DIABLO_FIT_REGISTRY[[cell$key]]), "no saved DIABLO fit file for this cell in this deployment")
  out <- multi_diablo_fit(cell)
  expect_true(out$ok)
  expect_true(inherits(out$fit, "block.splsda"))
  expect_null(out$error)
})

test_that("multi_diablo_fit() fails soft with a clear error when no fit is registered for a cell's key", {
  fake_cell <- list(key = "not_a_real_cell_key")
  out <- multi_diablo_fit(fake_cell)
  expect_false(out$ok)
  expect_null(out$fit)
  expect_true(nzchar(out$error))
})

test_that("multi_diablo_variance_df() tidies prop_expl_var into one data.frame per real omics block, excluding the internal Y block", {
  skip_if_not(MULTI_DATA_AVAILABLE, "multi-omics preloaded data not available in this deployment")
  cell <- multi_cell_by_key("female_Etanercept")
  fit_out <- multi_diablo_fit(cell)
  skip_if_not(fit_out$ok, "no saved DIABLO fit available to tidy")
  df <- multi_diablo_variance_df(fit_out$fit)
  expect_true(is.data.frame(df))
  expect_false("Y" %in% df$block)
  expect_true(all(c("expression", "methylation") %in% df$block))
  expect_true(all(df$variance_explained >= 0 & df$variance_explained <= 1))
})

test_that("multi_diablo_variance_df() returns NULL when prop_expl_var is missing or has no real blocks", {
  expect_null(multi_diablo_variance_df(list()))
  expect_null(multi_diablo_variance_df(list(prop_expl_var = list(Y = c(comp1 = 0.5)))))
})

## ---- multi_qc_scorecard() -------------------------------------------------------

test_that("multi_qc_scorecard() reports 'Data availability' honestly from MULTI_DATA_AVAILABLE, and 'warn' status for every unloaded sub-module", {
  scorecard <- multi_qc_scorecard(list())  ## nothing loaded this session
  by_label <- setNames(scorecard, vapply(scorecard, `[[`, character(1), "label"))
  expect_equal(by_label[["Data availability"]]$status, if (MULTI_DATA_AVAILABLE) "pass" else "fail")
  expect_equal(by_label[["Integration cell loaded"]]$status, "warn")
  expect_equal(by_label[["Patient stratification loaded"]]$status, "warn")
  expect_equal(by_label[["Biomarker Discovery signature loaded"]]$status, "warn")
  expect_equal(by_label[["Pathway enrichment loaded"]]$status, "warn")
})

test_that("multi_qc_scorecard() flips loaded-submodule items to 'pass' once multi_results carries their state", {
  r <- list(integration = list(cell = list(label = "Female - Etanercept (response)")),
            stratification = list(drug = "Etanercept"),
            biomarker = list(df = data.frame(feature = c("g1", "g2", "g1"))),
            pathway = list(df = data.frame(term = c("t1", "t2"))))
  scorecard <- multi_qc_scorecard(r)
  by_label <- setNames(scorecard, vapply(scorecard, `[[`, character(1), "label"))
  expect_equal(by_label[["Integration cell loaded"]]$status, "pass")
  expect_true(grepl("Female - Etanercept", by_label[["Integration cell loaded"]]$detail))
  expect_equal(by_label[["Patient stratification loaded"]]$status, "pass")
  expect_equal(by_label[["Biomarker Discovery signature loaded"]]$status, "pass")
  expect_true(grepl("^2 selected features", by_label[["Biomarker Discovery signature loaded"]]$detail))  ## unique() of g1,g2,g1 -> 2
  expect_equal(by_label[["Pathway enrichment loaded"]]$status, "pass")
})

test_that("multi_qc_scorecard() 'Model performance (honesty check)' is 'warn' when zero results exclude chance, 'pass' when at least one does", {
  r_none <- list(overview = list(summary36 = data.frame(excludes_chance = c(FALSE, FALSE, NA))))
  r_some <- list(overview = list(summary36 = data.frame(excludes_chance = c(TRUE, FALSE))))
  sc_none <- multi_qc_scorecard(r_none)
  sc_some <- multi_qc_scorecard(r_some)
  by_label_none <- setNames(sc_none, vapply(sc_none, `[[`, character(1), "label"))
  by_label_some <- setNames(sc_some, vapply(sc_some, `[[`, character(1), "label"))
  expect_equal(by_label_none[["Model performance (honesty check)"]]$status, "warn")
  expect_equal(by_label_some[["Model performance (honesty check)"]]$status, "pass")
})

## ---- multi_analysis_summary_table() --------------------------------------------

test_that("multi_analysis_summary_table() shows 'Not loaded'/'None loaded' placeholders, never a fabricated value, when nothing is loaded", {
  tbl <- multi_analysis_summary_table(multi_dataset = list(), multi_results = list())
  by_param <- setNames(tbl$Result, tbl$Parameter)
  expect_equal(by_param[["Samples analyzed (matched)"]], "Not loaded")
  expect_equal(by_param[["Active dataset table"]], "None loaded")
  expect_equal(by_param[["Active Multi-Omics Dataset source"]], "None selected yet")
  expect_equal(by_param[["Integration cell"]], "Not loaded")
  expect_equal(by_param[["SNF cohort loaded"]], "Not loaded")
})

test_that("multi_analysis_summary_table() reflects real loaded state once multi_dataset/multi_results carry it", {
  md <- list(active = TRUE, source = "geo", table_label = "GSE12345 (uploaded)")
  r <- list(overview = list(harmonization = list(n_matched = 42)),
            integration = list(cell = list(label = "Male - Adalimumab (response)"), snf_perf = NULL),
            stratification = list(drug = "Etanercept"))
  tbl <- multi_analysis_summary_table(md, r)
  by_param <- setNames(tbl$Result, tbl$Parameter)
  expect_equal(by_param[["Samples analyzed (matched)"]], "42")
  expect_equal(by_param[["Active dataset table"]], "GSE12345 (uploaded)")
  expect_equal(by_param[["Active Multi-Omics Dataset source"]], "NCBI GEO")
  expect_equal(by_param[["Integration cell"]], "Male - Adalimumab (response)")
  expect_equal(by_param[["Integration method(s)"]], "DIABLO")  ## no snf_perf -> DIABLO only
  expect_equal(by_param[["SNF cohort loaded"]], "Etanercept")
})

test_that("multi_analysis_summary_table() reports 'DIABLO + SNF' once snf_perf is present", {
  r <- list(integration = list(cell = list(label = "x"), snf_perf = data.frame(auc = 0.8)))
  tbl <- multi_analysis_summary_table(list(), r)
  by_param <- setNames(tbl$Result, tbl$Parameter)
  expect_equal(by_param[["Integration method(s)"]], "DIABLO + SNF")
})

## ---- multi_concordance_add_fdr() -----------------------------------------------

test_that("multi_concordance_add_fdr() adds BH-adjusted FDR columns for whichever raw p-value columns are present", {
  df <- data.frame(gene = c("g1", "g2", "g3"), expr_p = c(0.001, 0.04, 0.5), meth_p = c(0.002, 0.03, 0.6))
  out <- multi_concordance_add_fdr(df)
  expect_equal(out$expr_fdr, stats::p.adjust(df$expr_p, method = "BH"))
  expect_equal(out$meth_fdr, stats::p.adjust(df$meth_p, method = "BH"))
})

test_that("multi_concordance_add_fdr() leaves df untouched (no fabricated columns) when the raw p-value columns are absent, and passes NULL through", {
  df <- data.frame(gene = c("g1", "g2"), score = c(1, 2))
  out <- multi_concordance_add_fdr(df)
  expect_false("expr_fdr" %in% colnames(out))
  expect_false("meth_fdr" %in% colnames(out))
  expect_null(multi_concordance_add_fdr(NULL))
})

## ---- multi_active_dataset_banner() ---------------------------------------------

test_that("multi_active_dataset_banner() shows the 'no active dataset' note when nothing is active", {
  html <- htmltools::doRenderTags(multi_active_dataset_banner(list()))
  expect_true(grepl("No active Multi-Omics dataset", html))

  html_inactive <- htmltools::doRenderTags(multi_active_dataset_banner(list(source = "upload", active = FALSE)))
  expect_true(grepl("No active Multi-Omics dataset", html_inactive))
})

test_that("multi_active_dataset_banner() distinguishes preloaded (results shown) from upload/GEO (no stored results yet)", {
  html_preloaded <- htmltools::doRenderTags(multi_active_dataset_banner(list(source = "preloaded", active = TRUE)))
  expect_true(grepl("Preloaded Dataset", html_preloaded))
  expect_true(grepl("Existing results are available", html_preloaded))

  html_upload <- htmltools::doRenderTags(multi_active_dataset_banner(list(source = "upload", active = TRUE)))
  expect_true(grepl("User Upload", html_upload))
  expect_true(grepl("No stored results", html_upload))

  html_geo <- htmltools::doRenderTags(multi_active_dataset_banner(list(source = "geo", active = TRUE)))
  expect_true(grepl("NCBI GEO", html_geo))
})

## ---- multi_package_versions() ---------------------------------------------------

test_that("multi_package_versions() reports real installed package versions, 'not installed' for anything absent", {
  df <- multi_package_versions()
  expect_true(all(c("mixOmics", "SNFtool", "limma") %in% df$Package))
  ## Every reported version is either a real version string or the literal not-installed marker - never blank/NA.
  expect_true(all(nzchar(df$Version)))
  limma_row <- df[df$Package == "limma", ]
  expect_equal(limma_row$Version, as.character(utils::packageVersion("limma")))
})

## ---- multi_build_report() -------------------------------------------------------

test_that("multi_build_report() marks every sub-module 'not loaded this session' when multi_results is empty", {
  lines <- multi_build_report(list())
  text <- paste(lines, collapse = "\n")
  expect_true(grepl("## overview", text, fixed = TRUE))
  expect_equal(sum(grepl("(not loaded this session)", lines, fixed = TRUE)), 9L)  ## all 9 tracked sub-modules
})

test_that("multi_build_report() marks a loaded sub-module distinctly from the unloaded ones, and always lists limitations/reproducibility scripts", {
  lines <- multi_build_report(list(integration = list(cell = list(label = "x"))))
  text <- paste(lines, collapse = "\n")
  expect_true(grepl("Loaded - see the accompanying CSV", text, fixed = TRUE))
  expect_true(grepl("## Known limitations", text, fixed = TRUE))
  expect_true(grepl("## Reproducibility - source scripts", text, fixed = TRUE))
  ## overview still shows not-loaded since only "integration" was supplied.
  overview_idx <- which(lines == "## overview")
  expect_equal(lines[overview_idx + 1], "(not loaded this session)")
})
