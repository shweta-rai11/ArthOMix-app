## Guided "Start an analysis" workflow: the pure helpers that decide what each step may do.

suppressPackageStartupMessages(library(shiny))
source_from_app_root(file.path("R", "methylomics", "05_Differential_Methylation_Position", "mod_methyl_dmp.R"))
source_from_app_root(file.path("R", "workflow_guide.R"))

wf_meta <- function(n_hc_f = 10, n_hc_m = 6, n_ra_f = 12, n_ra_m = 5) {
  data.frame(
    sample = paste0("s", seq_len(n_hc_f + n_hc_m + n_ra_f + n_ra_m)),
    group = c(rep("HC", n_hc_f + n_hc_m), rep("RA", n_ra_f + n_ra_m)),
    sex = c(rep("F", n_hc_f), rep("M", n_hc_m), rep("F", n_ra_f), rep("M", n_ra_m)),
    stringsAsFactors = FALSE
  )
}

test_that("the tracker has the seven steps the reviewer asked for, in order", {
  expect_equal(vapply(WF_STEPS, `[[`, character(1), "label"),
               c("Choose data", "Check samples", "Define comparison", "Choose sex analysis",
                 "Run analysis", "Explore biomarkers", "Validate"))
})

test_that("every step destination is an existing sub-module id", {
  tx_ids <- c("overview", "dge", "interaction", "candidates", "biomarkercard", "diagnostic", "crosstissue")
  mx_ids <- c("qc", "dmp", "interaction", "candidates", "biomarkercard", "validation")
  tx <- WF_LAYERS$transcriptomics; mx <- WF_LAYERS$methylomics
  expect_true(all(c(tx$qc, tx$compare, tx$interaction, tx$candidates, tx$card, vapply(tx$validate, `[[`, "", "id")) %in% tx_ids))
  expect_true(all(c(mx$qc, mx$compare, mx$interaction, mx$candidates, mx$card, vapply(mx$validate, `[[`, "", "id")) %in% mx_ids))
})

test_that("wf_check_samples reports sample, group and sex counts for a complete dataset", {
  m <- wf_meta()
  chk <- wf_check_samples(m, wf_group_col(m), wf_sex_col(m))
  expect_true(chk$ok)
  expect_true(chk$sex_ok)
  expect_equal(chk$n, 33L)
  expect_equal(as.integer(chk$group_counts[c("HC", "RA")]), c(16L, 17L))
  expect_equal(as.integer(chk$sex_counts[c("F", "M")]), c(22L, 11L))
})

test_that("wf_check_samples blocks a dataset with no group column and explains why", {
  m <- wf_meta(); m$group <- NULL
  chk <- wf_check_samples(m, wf_group_col(m), wf_sex_col(m))
  expect_false(chk$ok)
  expect_null(wf_group_col(m))
  expect_match(chk$problems[1], "no group column")
})

test_that("wf_check_samples blocks a single-group dataset", {
  m <- wf_meta(); m$group <- "RA"
  chk <- wf_check_samples(m, "group", "sex")
  expect_false(chk$ok)
  expect_match(chk$problems[1], "At least two comparison groups are required")
})

test_that("wf_check_samples allows a dataset without sex but says only pooled is possible", {
  m <- wf_meta(); m$sex <- NULL
  chk <- wf_check_samples(m, "group", wf_sex_col(m))
  expect_true(chk$ok)
  expect_false(chk$sex_ok)
  expect_match(paste(chk$notes, collapse = " "), "Only the pooled analysis")
  m2 <- wf_meta(); m2$sex <- NA_character_
  expect_false(wf_check_samples(m2, "group", "sex")$sex_ok)
})

test_that("wf_check_samples rejects an empty or missing dataset", {
  expect_false(wf_check_samples(NULL, NULL, NULL)$ok)
  expect_match(wf_check_samples(NULL, NULL, NULL)$problems, "Please select a valid dataset")
})

test_that("wf_default_contrast puts a control-like group first and otherwise keeps the page default", {
  expect_equal(wf_default_contrast(c("RA", "HC", "RA")), list(ref = "HC", comp = "RA"))
  expect_equal(wf_default_contrast(c("RA", "Control")), list(ref = "Control", comp = "RA"))
  expect_equal(wf_default_contrast(c("Zeta", "Healthy")), list(ref = "Healthy", comp = "Zeta"))
  expect_equal(wf_default_contrast(c("B", "A")), list(ref = "A", comp = "B"))
  expect_null(wf_default_contrast(c("RA", NA)))
})

test_that("wf_check_contrast enables all three sex designs when every stratum is large enough", {
  cc <- wf_check_contrast(wf_meta(), "group", "HC", "RA", "sex", "transcriptomics")
  expect_true(cc$ok)
  expect_equal(cc$n_ref, 16L); expect_equal(cc$n_comp, 17L)
  expect_true(all(cc$modes))
  expect_equal(unname(cc$sex_levels), c("F", "M"))
  expect_equal(names(cc$sex_levels), c("Female", "Male"))
})

test_that("wf_check_contrast turns off stratified when one sex has too few samples in a group", {
  cc <- wf_check_contrast(wf_meta(n_ra_m = 2), "group", "HC", "RA", "sex", "methylomics")
  expect_true(cc$ok)
  expect_false(cc$modes[["stratified"]])
  expect_true(cc$modes[["sex_specific"]])
  expect_match(cc$mode_notes[["stratified"]], "Each sex needs")
})

test_that("wf_check_contrast turns off sex-specific when a sex-by-group cell is empty", {
  cc <- wf_check_contrast(wf_meta(n_hc_m = 0), "group", "HC", "RA", "sex", "methylomics")
  expect_false(cc$modes[["sex_specific"]])
  expect_false(cc$modes[["stratified"]])
})

test_that("wf_check_contrast allows only pooled without sex metadata", {
  m <- wf_meta(); m$sex <- NULL
  cc <- wf_check_contrast(m, "group", "HC", "RA", NULL, "transcriptomics")
  expect_true(cc$ok)
  expect_equal(cc$modes, c(pooled = TRUE, stratified = FALSE, sex_specific = FALSE))
})

test_that("wf_check_contrast: transcriptomics sex-specific needs the standard group column", {
  m <- wf_meta(); m$diagnosis <- m$group
  cc <- wf_check_contrast(m, "diagnosis", "HC", "RA", "sex", "transcriptomics")
  expect_false(cc$modes[["sex_specific"]])
  expect_true(cc$modes[["stratified"]])
})

test_that("wf_check_contrast rejects identical groups and too-small groups", {
  expect_match(wf_check_contrast(wf_meta(), "group", "RA", "RA", "sex", "transcriptomics")$problems, "must be different")
  small <- wf_meta(n_hc_f = 1, n_hc_m = 1)
  expect_match(wf_check_contrast(small, "group", "HC", "RA", "sex", "transcriptomics")$problems, "At least two comparison groups")
  expect_false(wf_check_contrast(wf_meta(), NULL, NULL, NULL, "sex", "transcriptomics")$ok)
})

test_that("wf_run_plan maps pooled, stratified and sex-specific onto the existing analyses", {
  lv <- c(Female = "F", Male = "M")
  p <- wf_run_plan("transcriptomics", "pooled", "sex", lv)
  expect_length(p, 1); expect_equal(p[[1]]$kind, "dge"); expect_true(p[[1]]$adjust_sex); expect_null(p[[1]]$sex_level)
  s <- wf_run_plan("methylomics", "stratified", "sex", lv)
  expect_length(s, 2)
  expect_equal(vapply(s, `[[`, "", "kind"), c("dmp", "dmp"))
  expect_equal(vapply(s, `[[`, "", "sex_level"), c("F", "M"))
  expect_equal(vapply(s, `[[`, "", "label"), c("Female samples only", "Male samples only"))
  x <- wf_run_plan("transcriptomics", "sex_specific", "sex", lv)
  expect_equal(x[[1]]$kind, "interaction")
  np <- wf_run_plan("transcriptomics", "pooled", NULL, character(0))
  expect_false(np[[1]]$adjust_sex)
})

test_that("wf_summary_line reports only stored values and never invents missing ones", {
  dge <- list(contrast = "RA vs HC (group) adjusted for sex", n_samples = 183, n_tested = 15763,
              n_significant = 1200, n_up = 700, n_down = 500)
  expect_equal(wf_summary_line("dge", dge, padj = 0.05, lfc = 0.5),
               "RA vs HC (group) adjusted for sex: 183 samples, 15,763 genes tested, 1,200 significant (700 up, 500 down) at adjusted p < 0.05 and |log2FC| > 0.5.")
  dmp <- list(comparison = "RA vs Control (F)", n_probes = 400000, n_sig = 12)
  expect_equal(wf_summary_line("dmp", dmp), "RA vs Control (F): 400,000 CpGs tested, 12 significant at FDR < 0.05.")
  expect_null(wf_summary_line("dge", NULL))
})

test_that("tracker states: done, current and upcoming", {
  expect_equal(wf_step_state(1, 3, c(1L, 2L), TRUE), "done")
  expect_equal(wf_step_state(3, 3, c(1L, 2L), TRUE), "current")
  expect_equal(wf_step_state(5, 3, c(1L, 2L), TRUE), "future")
  expect_equal(wf_step_state(1, 1, integer(0), FALSE), "future")
  html <- as.character(wf_home_tracker_ui(3L, c(1L, 2L), TRUE))
  expect_equal(lengths(regmatches(html, gregexpr("home-hero-flow-step", html))), 7L)
  expect_equal(lengths(regmatches(html, gregexpr("wf-state-done", html))), 2L)
  expect_equal(lengths(regmatches(html, gregexpr("wf-state-current", html))), 1L)
})

test_that("wf_check_contrast on an unreadable comparison has no per-mode notes (the bar must not index them)", {
  cc <- wf_check_contrast(wf_meta(), NULL, NULL, NULL, "sex", "transcriptomics")
  expect_false(cc$ok)
  expect_length(cc$mode_notes, 0)
  expect_false("stratified" %in% names(cc$mode_notes))
})
