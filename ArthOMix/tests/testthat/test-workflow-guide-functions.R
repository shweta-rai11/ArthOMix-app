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

test_that("the tracker has the eight steps the reviewer asked for, in order", {
  ## "Define comparison" is folded into "Choose analysis" - the comparison defaults from the
  ## dataset and the sex design needs it resolved anyway, so locking both together in one step
  ## removes a step that was otherwise just duplicating what Choose analysis already needed.
  ## Check samples comes before Choose analysis - see the cohort composition before picking a design.
  expect_equal(vapply(WF_STEPS, `[[`, character(1), "label"),
               c("Choose data", "Check samples", "Choose analysis",
                 "Run analysis", "Explore biomarkers", "Biomarker modeling", "Validate",
                 "Interpret biomarkers"))
  expect_equal(WF_N_STEPS, 8L)
  expect_equal(vapply(WF_STEPS, `[[`, integer(1), "n"), 1:8)
  ## Step 5 no longer sends people to the Biomarker Card; that is step 8's job.
  expect_false(grepl("Biomarker Card", WF_STEPS[[5]]$blurb, fixed = TRUE))
  expect_match(WF_STEPS[[8]]$blurb, "Biomarker Card", fixed = TRUE)
  expect_true(all(nzchar(vapply(WF_STEPS, `[[`, character(1), "icon"))))
})

test_that("every step destination is an existing sub-module id", {
  tx_ids <- c("overview", "preprocessing", "dge", "wgcna", "candidates", "mr", "coloc",
              "featureselection", "diagnostic", "interaction", "crosstissue", "crossancestry",
              "enrichment", "deconvolution", "biomarkercard")
  mx_ids <- c("qc", "dmp", "interaction", "candidates", "featureselection", "diagnostic",
              "biomarkercard", "validation")
  ids <- function(lay) c(lay$qc, lay$compare, lay$interaction, lay$candidates, lay$card,
                         vapply(lay$explore, `[[`, "", "id"), vapply(lay$model, `[[`, "", "id"),
                         vapply(lay$validate, `[[`, "", "id"))
  expect_true(all(ids(WF_LAYERS$transcriptomics) %in% tx_ids))
  expect_true(all(ids(WF_LAYERS$methylomics) %in% mx_ids))
  ## Transcriptomics: feature selection and the classifier in step 7, the three validation routes
  ## in step 8, the card in step 9.
  tx <- WF_LAYERS$transcriptomics
  expect_equal(vapply(tx$model, `[[`, "", "id"), c("featureselection", "diagnostic"))
  expect_equal(vapply(tx$validate, `[[`, "", "id"), c("diagnostic", "crosstissue", "crossancestry"))
  expect_equal(tx$validate[[1]]$tabset, "main_tabs")
  expect_equal(tx$validate[[1]]$tab, "External Validation")
  expect_equal(tx$card, "biomarkercard")
  ## Methylomics keeps every id it used before.
  expect_equal(WF_LAYERS$methylomics$validate[[1]]$id, "validation")
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
  expect_equal(lengths(regmatches(html, gregexpr("home-hero-flow-step", html))), 8L)
  expect_equal(lengths(regmatches(html, gregexpr("wf-state-done", html))), 2L)
  expect_equal(lengths(regmatches(html, gregexpr("wf-state-current", html))), 1L)
})

test_that("wf_check_contrast on an unreadable comparison has no per-mode notes (the bar must not index them)", {
  cc <- wf_check_contrast(wf_meta(), NULL, NULL, NULL, "sex", "transcriptomics")
  expect_false(cc$ok)
  expect_length(cc$mode_notes, 0)
  expect_false("stratified" %in% names(cc$mode_notes))
})

## ---------------------------------------------------------------------------------------------------
## TX_STEP_MAP and the pure helpers the step-driven Transcriptomics shell is built on.
## No Shiny: dataset, results and settings are plain lists here, exactly as the server passes them.
## ---------------------------------------------------------------------------------------------------

## Titles are normally read from each module's own config; the registry is not loaded in this file, so
## the tests pass the same titles in explicitly rather than asserting on bare ids.
wf_test_titles <- function() {
  wf_step_titles(registry = list(
    overview         = list(config = list(title = "Overview and Datasets")),
    preprocessing    = list(config = list(title = "Preprocessing and Batch Correction")),
    dge              = list(config = list(title = "Differential Expression")),
    wgcna            = list(config = list(title = "WGCNA Co-expression Network")),
    mr               = list(config = list(title = "Mendelian Randomization")),
    coloc            = list(config = list(title = "Colocalization")),
    interaction      = list(config = list(title = "Sex Interaction Analysis")),
    deconvolution    = list(config = list(title = "Immune Deconvolution")),
    candidates       = list(config = list(title = "Candidate Gene Identification")),
    enrichment       = list(config = list(title = "Functional Enrichment")),
    featureselection = list(config = list(title = "Feature Selection")),
    diagnostic       = list(config = list(title = "Diagnostic Model")),
    crosstissue      = list(config = list(title = "Cross-Tissue Replication")),
    crossancestry    = list(config = list(title = "Cross-Ancestry MR Replication")),
    biomarkercard    = list(config = list(title = "Biomarker Card"))
  ))
}

## The seven dataset/result states the shell has to get right, as plain lists.
wf_ds <- function(ok = TRUE, sex_ok = TRUE, load_id = 1L) list(ok = ok, sex_ok = sex_ok, load_id = load_id)

wf_state <- function(results = list(), contrast = list(col = "group", ref = "HC", comp = "RA"),
                     sex_mode = "pooled", ds = wf_ds(), stamp_keys = character(0)) {
  stamps <- stats::setNames(as.list(seq_along(stamp_keys)), stamp_keys)
  set <- list(contrast = contrast, sex_mode = sex_mode, stamps = stamps)
  set$fingerprints <- stats::setNames(
    lapply(stamp_keys, function(k) wf_fingerprint(k, ds, set)), stamp_keys)
  list(ds = ds, results = results, settings = set)
}

wf_st <- function(st, key) {
  wf_module_status(key, st$ds, st$results, st$settings, titles = wf_test_titles())
}
wf_all <- function(st) {
  vapply(names(TX_STEP_MAP), function(k) wf_st(st, k)$status, character(1))
}

test_that("TX_STEP_MAP covers the Dataset tab's 15 sub-modules plus the second Diagnostic entry", {
  expect_equal(length(TX_STEP_MAP), 16L)
  expect_equal(length(unique(vapply(TX_STEP_MAP, `[[`, character(1), "id"))), 15L)
  expect_equal(names(TX_STEP_MAP), TX_STEP_MAP_ORDER)
  ## Steps are non-decreasing along TX_STEP_MAP_ORDER, so that order is already pipeline order.
  steps <- vapply(TX_STEP_MAP, `[[`, integer(1), "step")
  expect_false(is.unsorted(steps))
  expect_equal(sort(unique(unname(steps))), c(2L, 4L, 5L, 6L, 7L, 8L))
  expect_true(all(steps >= 2L & steps <= 8L))
})

test_that("the map encodes the dependencies the modules' own code enforces", {
  expect_equal(TX_STEP_MAP$candidates$hard, list("dge", "wgcna"))
  expect_equal(TX_STEP_MAP$candidates$soft, c("mr", "coloc"))
  expect_equal(TX_STEP_MAP$featureselection$hard, list("candidates"))
  expect_equal(TX_STEP_MAP$diagnostic$hard, list(c("featureselection", "candidates")))
  expect_equal(TX_STEP_MAP$diagnostic_external$hard, list("diagnostic"))
  expect_equal(TX_STEP_MAP$crosstissue$hard, list("featureselection"))
  expect_equal(TX_STEP_MAP$interaction$hard, list("dataset", "sex"))
  expect_equal(TX_STEP_MAP$dge$hard, list("dataset", "contrast"))
  ## mr, coloc, crossancestry and enrichment run off bundled data or a pasted list: never locked.
  for (k in c("mr", "coloc", "crossancestry", "enrichment")) expect_length(TX_STEP_MAP[[k]]$hard, 0)
  ## Preprocessing stores into results$study_meta, and External Validation into diagnostic[[sex]]$external.
  expect_equal(TX_STEP_MAP$preprocessing$slot, "study_meta")
  expect_equal(TX_STEP_MAP$dge$slot, "dge_runs")
  expect_equal(TX_STEP_MAP$diagnostic_external$slot, "diagnostic")
  expect_equal(TX_STEP_MAP$diagnostic_external$subfield, "external")
})

test_that("diagnostic appears in step 6 and step 7 as one tab with two entries", {
  expect_equal(TX_STEP_MAP$diagnostic$step, 6L)
  expect_equal(TX_STEP_MAP$diagnostic_external$step, 7L)
  expect_equal(TX_STEP_MAP$diagnostic$id, TX_STEP_MAP$diagnostic_external$id)
  expect_false(TX_STEP_MAP$diagnostic$alias)
  expect_true(TX_STEP_MAP$diagnostic_external$alias)
  expect_equal(TX_STEP_MAP$diagnostic_external$tabset, "main_tabs")
  expect_equal(TX_STEP_MAP$diagnostic_external$tab, "External Validation")
  ## Exactly one alias, so the catalogue still shows one card per module.
  expect_equal(sum(vapply(TX_STEP_MAP, `[[`, logical(1), "alias")), 1L)
})

test_that("empty session: only the dataset-free and dataset-only modules are ready", {
  st <- wf_state(contrast = NULL)
  s <- wf_all(st)
  expect_equal(unname(s[c("overview", "preprocessing", "wgcna", "deconvolution", "mr", "coloc",
                          "crossancestry", "enrichment", "biomarkercard", "interaction")]),
               rep("ready", 10))
  expect_equal(unname(s[c("dge", "candidates", "featureselection", "diagnostic",
                          "diagnostic_external", "crosstissue")]), rep("locked", 6))
  expect_equal(wf_st(st, "dge")$reason, "Needs: a comparison (Step 3).")
  expect_equal(wf_st(st, "candidates")$reason, "Needs: Differential Expression and WGCNA Co-expression Network.")
  expect_equal(wf_st(st, "diagnostic")$reason, "Needs: Feature Selection or Candidate Gene Identification.")
  expect_equal(wf_st(st, "diagnostic_external")$reason, "Needs: Diagnostic Model.")
})

test_that("no dataset locks everything that needs one, and no sex column locks Sex Interaction", {
  none <- wf_state(contrast = NULL, ds = wf_ds(ok = FALSE, sex_ok = FALSE))
  expect_equal(wf_st(none, "overview")$status, "locked")
  expect_equal(wf_st(none, "overview")$reason, "Needs: a loaded dataset.")
  expect_equal(wf_st(none, "biomarkercard")$status, "locked")
  ## Bundled-data modules are still runnable with nothing loaded at all.
  expect_equal(wf_st(none, "mr")$status, "ready")
  expect_equal(wf_st(none, "crossancestry")$status, "ready")
  nosex <- wf_state(ds = wf_ds(sex_ok = FALSE))
  expect_equal(wf_st(nosex, "interaction")$status, "locked")
  expect_match(wf_st(nosex, "interaction")$reason, "sex column with at least two values")
  expect_equal(wf_st(nosex, "dge")$status, "ready")
})

test_that("DGE only: candidates still waits for WGCNA", {
  st <- wf_state(results = list(dge_runs = list(r1 = list())), stamp_keys = "dge")
  expect_equal(wf_st(st, "dge")$status, "done")
  expect_equal(wf_st(st, "wgcna")$status, "ready")
  expect_equal(wf_st(st, "candidates")$status, "locked")
  expect_equal(wf_st(st, "candidates")$reason, "Needs: WGCNA Co-expression Network.")
})

test_that("DGE + WGCNA: candidates unlocks, and its downstream is still locked", {
  st <- wf_state(results = list(dge_runs = list(r1 = list()), wgcna = list(module_genes = list(blue = "A"))),
                 stamp_keys = c("dge", "wgcna"))
  expect_equal(wf_st(st, "candidates")$status, "ready")
  expect_equal(wf_st(st, "featureselection")$status, "locked")
  expect_equal(wf_st(st, "featureselection")$reason, "Needs: Candidate Gene Identification.")
  expect_equal(wf_st(st, "diagnostic")$status, "locked")
})

test_that("panel built: feature selection done unlocks the classifier and cross-tissue", {
  st <- wf_state(
    results = list(dge_runs = list(r1 = list()), wgcna = list(module_genes = list(blue = "A")),
                   candidates = list(final = list(genes = "A")),
                   featureselection = list(pooled = list(consensus_genes = c("A", "B")))),
    stamp_keys = c("dge", "wgcna", "candidates", "featureselection"))
  expect_equal(wf_st(st, "featureselection")$status, "done")
  expect_equal(wf_st(st, "diagnostic")$status, "ready")
  expect_equal(wf_st(st, "crosstissue")$status, "ready")
  expect_equal(wf_st(st, "diagnostic_external")$status, "locked")
})

test_that("model trained: External Validation unlocks, and is done only once it has scored a cohort", {
  trained <- list(dge_runs = list(r1 = list()), wgcna = list(module_genes = list(blue = "A")),
                  candidates = list(final = list(genes = "A")),
                  featureselection = list(pooled = list(consensus_genes = c("A", "B"))),
                  diagnostic = list(pooled = list(auc = 0.8)))
  keys <- c("dge", "wgcna", "candidates", "featureselection", "diagnostic")
  st <- wf_state(results = trained, stamp_keys = keys)
  expect_equal(wf_st(st, "diagnostic")$status, "done")
  ## A trained model with no external run yet: unlocked, but not done.
  expect_equal(wf_st(st, "diagnostic_external")$status, "ready")
  scored <- trained
  scored$diagnostic$pooled$external <- list(cohort = "GSE15573", models_scored = TRUE)
  st2 <- wf_state(results = scored, stamp_keys = c(keys, "diagnostic_external"))
  expect_equal(wf_st(st2, "diagnostic_external")$status, "done")
  ## The two entries read the same slot but report separately.
  expect_equal(wf_st(st2, "diagnostic")$status, "done")
})

test_that("changing the contrast makes DGE stale and everything downstream with it", {
  res <- list(dge_runs = list(r1 = list()), wgcna = list(module_genes = list(blue = "A")),
              candidates = list(final = list(genes = "A")),
              featureselection = list(pooled = list(consensus_genes = c("A", "B"))),
              diagnostic = list(pooled = list(auc = 0.8, external = list(cohort = "GSE15573"))),
              crosstissue = list(pooled = list()), biomarkercard = list(gene = "A"))
  keys <- c("dge", "wgcna", "candidates", "featureselection", "diagnostic",
            "diagnostic_external", "crosstissue", "biomarkercard")
  st <- wf_state(results = res, stamp_keys = keys)
  expect_equal(unname(wf_all(st)[keys]), rep("done", length(keys)))

  st2 <- st
  st2$settings$contrast <- list(col = "group", ref = "HC", comp = "OA")
  s <- wf_all(st2)
  expect_equal(unname(s[c("dge", "candidates", "featureselection", "diagnostic",
                          "diagnostic_external", "biomarkercard")]), rep("stale", 6))
  ## Cross-tissue does not fingerprint the contrast itself, but its upstream panel is stale.
  expect_equal(unname(s[["crosstissue"]]), "stale")
  ## WGCNA is unsupervised: a different contrast does not invalidate it.
  expect_equal(unname(s[["wgcna"]]), "done")
  expect_equal(wf_st(st2, "dge")$reason, "Stale - re-run. the comparison changed since this ran.")
  expect_match(wf_st(st2, "candidates")$reason, "^Stale - re-run\\.")
})

test_that("changing the sex design makes the sex-aware results stale, WGCNA included only via upstream", {
  res <- list(dge_runs = list(r1 = list()), wgcna = list(module_genes = list(blue = "A")),
              candidates = list(final = list(genes = "A")),
              deconvolution = list(fractions = 1))
  st <- wf_state(results = res, stamp_keys = c("dge", "wgcna", "candidates", "deconvolution"))
  st2 <- st; st2$settings$sex_mode <- "stratified"
  s <- wf_all(st2)
  expect_equal(unname(s[["dge"]]), "stale")
  expect_equal(unname(s[["candidates"]]), "stale")
  expect_equal(unname(s[["wgcna"]]), "done")
  expect_equal(unname(s[["deconvolution"]]), "done")
  expect_equal(wf_st(st2, "dge")$reason, "Stale - re-run. the sex design changed since this ran.")
})

test_that("loading a different dataset makes every dataset-backed result stale", {
  res <- list(overview = list(n = 1), dge_runs = list(r1 = list()), wgcna = list(module_genes = list(blue = "A")),
              mr = list(hits = 1))
  st <- wf_state(results = res, stamp_keys = c("overview", "dge", "wgcna", "mr"))
  st2 <- st; st2$ds <- wf_ds(load_id = 2L)
  s <- wf_all(st2)
  expect_equal(unname(s[c("overview", "dge", "wgcna")]), rep("stale", 3))
  expect_equal(wf_st(st2, "wgcna")$reason, "Stale - re-run. the loaded dataset changed since this ran.")
  ## Bundled-data MR has no dataset in its fingerprint.
  expect_equal(unname(s[["mr"]]), "done")
})

test_that("re-running feature selection makes the classifier and everything after it stale", {
  res <- list(dge_runs = list(r1 = list()), wgcna = list(module_genes = list(blue = "A")),
              candidates = list(final = list(genes = "A")),
              featureselection = list(pooled = list(consensus_genes = c("A", "B"))),
              diagnostic = list(pooled = list(auc = 0.8, external = list(cohort = "GSE15573"))),
              crosstissue = list(pooled = list()), crossancestry = list(female = list()),
              biomarkercard = list(gene = "A"))
  keys <- c("dge", "wgcna", "candidates", "featureselection", "diagnostic",
            "diagnostic_external", "crosstissue", "crossancestry", "biomarkercard")
  st <- wf_state(results = res, stamp_keys = keys)
  expect_equal(unname(wf_all(st)[keys]), rep("done", length(keys)))
  ## A new Feature Selection run gets a new stamp; nothing else re-runs.
  st2 <- st
  st2$settings$stamps$featureselection <- 99L
  st2$settings$fingerprints$featureselection <- wf_fingerprint("featureselection", st2$ds, st2$settings)
  s <- wf_all(st2)
  expect_equal(unname(s[["featureselection"]]), "done")
  expect_equal(unname(s[c("diagnostic", "diagnostic_external", "crosstissue", "biomarkercard")]),
               rep("stale", 4))
  expect_equal(unname(s[c("dge", "wgcna", "candidates")]), rep("done", 3))
  ## Cross-Ancestry MR reads no upstream result, so it is untouched.
  expect_equal(unname(s[["crossancestry"]]), "done")
  expect_match(wf_st(st2, "diagnostic")$reason, "Feature Selection changed since this ran", fixed = TRUE)
})

test_that("wf_required_closure adds what a selection cannot run without, and says why", {
  titles <- wf_test_titles()
  cl <- wf_required_closure("candidates", titles = titles)
  expect_equal(cl$keys, c("dge", "wgcna", "candidates"))
  expect_equal(cl$added$dge, "Added automatically: Candidate Gene Identification needs Differential Expression.")
  expect_equal(cl$added$wgcna, "Added automatically: Candidate Gene Identification needs WGCNA Co-expression Network.")
  ## Any-of group: Feature Selection is the first alternative the Diagnostic Model can use.
  expect_true("featureselection" %in% wf_required_closure("diagnostic", titles = titles)$keys)
  ## An already-selected alternative is not duplicated.
  expect_false("featureselection" %in% names(wf_required_closure(c("diagnostic", "candidates"), titles = titles)$added))
  ## Selecting only External Validation pulls the whole chain behind it.
  ev <- wf_required_closure("diagnostic_external", titles = titles)
  expect_equal(ev$keys, c("dge", "wgcna", "candidates", "featureselection", "diagnostic", "diagnostic_external"))
  ## Modules with no hard requirements add nothing.
  expect_length(wf_required_closure(c("mr", "coloc", "crossancestry", "enrichment"))$added, 0)
  expect_equal(wf_required_closure(character(0))$keys, character(0))
})

test_that("wf_run_order puts mr/coloc before candidates and feature selection before the classifier", {
  ## Deliberately reversed input: the order must come from the dependencies, not the input order.
  o <- wf_run_order(c("biomarkercard", "diagnostic", "featureselection", "candidates",
                      "coloc", "mr", "wgcna", "dge"))
  expect_lt(match("mr", o), match("candidates", o))
  expect_lt(match("coloc", o), match("candidates", o))
  expect_lt(match("dge", o), match("candidates", o))
  expect_lt(match("wgcna", o), match("candidates", o))
  expect_lt(match("candidates", o), match("featureselection", o))
  expect_lt(match("featureselection", o), match("diagnostic", o))
  expect_equal(o[length(o)], "biomarkercard")
  expect_setequal(o, c("biomarkercard", "diagnostic", "featureselection", "candidates",
                       "coloc", "mr", "wgcna", "dge"))
  expect_equal(wf_run_order(character(0)), character(0))
  expect_equal(wf_run_order("dge"), "dge")
  ## Unknown keys are dropped rather than breaking the order.
  expect_equal(wf_run_order(c("dge", "not_a_module")), "dge")
})

test_that("each preset resolves to a closed, dependency-ordered plan", {
  titles <- wf_test_titles()
  std <- wf_preset_plan("standard", titles = titles)
  expect_equal(std$keys, c("dge", "wgcna", "candidates", "featureselection", "diagnostic",
                           "diagnostic_external", "biomarkercard"))
  expect_equal(std$sex_mode, "pooled")

  sx <- wf_preset_plan("sex_differences", titles = titles)
  expect_equal(sx$sex_mode, "sex_specific")
  ## Sex-specific adds the interaction model and still runs the main DGE alongside it.
  expect_true(all(c("interaction", "dge") %in% sx$keys))
  expect_true(all(c("candidates", "featureselection", "crosstissue") %in% sx$keys))

  full <- wf_preset_plan("full", titles = titles)
  expect_setequal(full$keys, TX_STEP_MAP_ORDER)
  expect_equal(full$keys, wf_run_order(TX_STEP_MAP_ORDER))

  cust <- wf_preset_plan("custom", titles = titles)
  expect_equal(cust$keys, character(0))
  ## Custom with one late pick pulls its whole chain in, in order, and explains each addition.
  picked <- wf_preset_plan("custom", extra = "crosstissue", titles = titles)
  expect_equal(picked$keys, c("dge", "wgcna", "candidates", "featureselection", "crosstissue"))
  expect_true(all(grepl("^Added automatically", unlist(picked$added))))
  ## Sex-specific on a custom pick still brings the interaction model in.
  expect_true("interaction" %in% wf_preset_plan("custom", sex_mode = "sex_specific", extra = "dge")$keys)
})

test_that("every preset plan is runnable in the order it gives", {
  for (nm in names(TX_PRESETS)) {
    plan <- wf_preset_plan(nm, titles = wf_test_titles())
    for (i in seq_along(plan$keys)) {
      hard <- intersect(unlist(TX_STEP_MAP[[plan$keys[i]]]$hard), TX_STEP_MAP_ORDER)
      ## Any-of groups are satisfied by at least one earlier entry; hard singletons by that one.
      for (grp in TX_STEP_MAP[[plan$keys[i]]]$hard) {
        grp <- intersect(grp, TX_STEP_MAP_ORDER)
        if (!length(grp)) next
        expect_true(any(grp %in% utils::head(plan$keys, i - 1)),
                    info = sprintf("%s: %s not scheduled before %s", nm, paste(grp, collapse = "/"), plan$keys[i]))
      }
    }
  }
})

test_that("a step's rolled-up status is the worst thing in it", {
  expect_equal(wf_step_status(c("done", "done")), "done")
  expect_equal(wf_step_status(c("done", "stale")), "stale")
  expect_equal(wf_step_status(c("locked", "locked")), "locked")
  expect_equal(wf_step_status(c("locked", "ready")), "ready")
  expect_equal(wf_step_status(c("done", "ready")), "ready")
  expect_equal(wf_step_status(character(0)), "ready")
  expect_equal(sort(names(WF_STATUS_LABEL)), c("done", "locked", "ready", "stale"))
})

test_that("a result with no recorded fingerprint is reported done, not falsely stale", {
  st <- wf_state(results = list(wgcna = list(module_genes = list(blue = "A"))))
  st$settings$fingerprints <- list()
  expect_equal(wf_st(st, "wgcna")$status, "done")
})
