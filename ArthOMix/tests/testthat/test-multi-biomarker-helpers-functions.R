## Module 3 (Multiomics) - Biomarker Discovery's own pure functions
## (multiomics_biomarker_helpers.R): block role selection, class-label
## friendliness (display-only, never changing the filtered value), the

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_plots.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "multiomics_dataset_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "06_Gene_CpG_Mapping", "multiomics_mapping_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "02_Cohort_Harmonization", "cohort_harmonization_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_integration_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "05_Biomarker_Discovery", "multiomics_biomarker_helpers.R"))

test_that("mb_select_blocks() re-labels the user-chosen roles to fixed Transcriptomics/Methylomics names", {
  layers <- list(RNA = matrix(1, 2, 2), Meth = matrix(2, 2, 2), Extra = matrix(3, 2, 2))
  out <- mb_select_blocks(layers, transcript_block = "RNA", methyl_block = "Meth")
  expect_equal(names(out), c("Transcriptomics", "Methylomics"))
  expect_identical(out$Transcriptomics, layers$RNA)
})

test_that("mb_select_blocks() returns NULL when the same block is assigned to both roles, or a role doesn't exist", {
  layers <- list(RNA = matrix(1, 2, 2), Meth = matrix(2, 2, 2))
  expect_null(mb_select_blocks(layers, "RNA", "RNA"))
  expect_null(mb_select_blocks(layers, "RNA", "NotThere"))
  expect_null(mb_select_blocks(NULL, "RNA", "Meth"))
})

test_that("mb_friendly_class_label() maps known responder/non-responder abbreviations, case/whitespace-insensitively", {
  expect_equal(unname(mb_friendly_class_label(c("resp", " NON ", "Response"))), c("Responder", "Non-responder", "Responder"))
})

test_that("mb_friendly_class_label() passes through any value not on the synonym list unchanged (never invents a label)", {
  expect_equal(unname(mb_friendly_class_label(c("HC", "RA", "Class1"))), c("HC", "RA", "Class1"))
})

test_that("mb_variance_prefilter() keeps exactly the top-N highest-variance columns without reading any outcome", {
  set.seed(300)
  mat <- matrix(rnorm(60), 6, 10, dimnames = list(paste0("S", 1:6), paste0("f", 1:10)))
  vars <- apply(mat, 2, var)
  top4 <- names(sort(vars, decreasing = TRUE))[1:4]
  out <- mb_variance_prefilter(mat, max_features = 4)
  expect_setequal(colnames(out), top4)
})

test_that("mb_variance_prefilter() passes the matrix through unchanged when it's already within the cap", {
  mat <- matrix(rnorm(20), 5, 4)
  expect_identical(mb_variance_prefilter(mat, max_features = 10), mat)
  expect_identical(mb_variance_prefilter(mat, max_features = NULL), mat)
})

test_that("mb_default_max_features() caps at 5000 but never exceeds the actual feature count", {
  expect_equal(mb_default_max_features(100), 100L)
  expect_equal(mb_default_max_features(20000), 5000L)
})

test_that("mb_data_check_table() reports 'Not available' rows when validation itself failed", {
  tbl <- mb_data_check_table(list(ok = FALSE), NULL, list(ok = FALSE, reason = "no data"))
  expect_equal(tbl$Status[tbl$Check == "Transcriptomics"], "Not available")
  expect_equal(tbl$Status[tbl$Check == "DIABLO eligibility"], "no data")
})

test_that("mb_data_check_table() reports real per-block sample/feature counts, missingness severity, and outcome class summary", {
  validation <- list(
    ok = TRUE, reliable_matching = TRUE, n_shared = 25,
    per_block = list(
      Transcriptomics = list(ok = TRUE, n_samples = 25, n_features = 500, pct_missing = 1.2),
      Methylomics = list(ok = TRUE, n_samples = 25, n_features = 800, pct_missing = 8.5)
    )
  )
  outcome_summary <- list(type = "categorical", n_classes = 2, n = 25, imbalanced = FALSE)
  eligibility <- list(ok = TRUE)
  tbl <- mb_data_check_table(validation, outcome_summary, eligibility)
  expect_true(grepl("25.*500", tbl$Status[tbl$Check == "Transcriptomics"]))
  expect_true(grepl("High \\(8.5%", tbl$Status[tbl$Check == "Missing values"]))
  expect_true(grepl("2 classes", tbl$Status[tbl$Check == "Outcome"]))
  expect_equal(tbl$Status[tbl$Check == "DIABLO eligibility"], "Ready")
})

test_that("mb_data_check_table() reports a continuous outcome as not usable for DIABLO as-is", {
  validation <- list(ok = TRUE, reliable_matching = TRUE, n_shared = 25,
                       per_block = list(Transcriptomics = list(ok = TRUE, n_samples = 25, n_features = 500, pct_missing = 0),
                                          Methylomics = list(ok = TRUE, n_samples = 25, n_features = 500, pct_missing = 0)))
  outcome_summary <- list(type = "continuous", n = 25)
  tbl <- mb_data_check_table(validation, outcome_summary, list(ok = FALSE, reason = "x"))
  expect_true(grepl("not usable for DIABLO", tbl$Status[tbl$Check == "Outcome"]))
})

test_that("mb_stability_category() cuts selection frequency at the documented 0.5/0.8 thresholds", {
  out <- mb_stability_category(c(0.1, 0.49, 0.5, 0.79, 0.8, 1.0))
  expect_equal(as.character(out), c("Low stability", "Low stability", "Moderately stable", "Moderately stable", "Stable", "Stable"))
})

test_that("mb_component_correlation() computes the real correlation between two blocks' first-component scores", {
  ids <- paste0("S", 1:10)
  v1 <- rnorm(10); v2 <- v1 * 0.8 + rnorm(10, sd = 0.1)
  fit <- list(variates = list(
    Transcriptomics = matrix(v1, 10, 1, dimnames = list(ids, "comp1")),
    Methylomics = matrix(v2, 10, 1, dimnames = list(ids, "comp1")),
    Y = matrix(0, 10, 1, dimnames = list(ids, "comp1"))
  ))
  out <- mb_component_correlation(fit)
  expect_equal(out$r, stats::cor(v1, v2))
  expect_equal(out$n, 10L)
})

test_that("mb_component_correlation() returns NULL when there aren't exactly 2 real (non-Y) blocks", {
  fit_one <- list(variates = list(Transcriptomics = matrix(1, 5, 1), Y = matrix(0, 5, 1)))
  expect_null(mb_component_correlation(fit_one))
})

test_that("mb_reproducibility_table() reports every real parameter used, never a placeholder", {
  diablo_res <- list(params = list(
    n_samples = 30, classes = c("HC", "RA"), ncomp = 2,
    keepX = list(Transcriptomics = c(10, 20), Methylomics = c(15, 25)),
    design = matrix(c(0, 0.1, 0.1, 0), 2, 2, dimnames = list(c("Transcriptomics", "Methylomics"), c("Transcriptomics", "Methylomics"))),
    distance = "max.dist", validation_method = "Mfold", folds = 5, nrepeat = 3
  ))
  tbl <- mb_reproducibility_table(diablo_res, dataset_label = "My Dataset", preprocessing_note = NULL, seed = 42)
  by_param <- setNames(tbl$Value, tbl$Parameter)
  expect_equal(by_param[["Data source"]], "My Dataset")
  expect_equal(by_param[["Samples analyzed (matched)"]], "30")
  expect_equal(by_param[["keepX (Transcriptomics)"]], "10,20")
  expect_equal(by_param[["Design (Transcriptomics <-> Methylomics)"]], "0.10")
  expect_equal(by_param[["Random seed"]], "42")
})

test_that("mb_matched_sample_table() lists every matched sample id with its real outcome value, NA when no outcome is available", {
  ids <- c("S1", "S2")
  outcome <- c(S1 = "RA", S2 = "HC")
  tbl <- mb_matched_sample_table(list(ok = TRUE), outcome, ids)
  expect_equal(tbl$outcome, c("RA", "HC"))
  tbl_no_outcome <- mb_matched_sample_table(list(ok = TRUE), NULL, ids)
  expect_true(all(is.na(tbl_no_outcome$outcome)))
})

test_that("mb_signature_table()/mb_cv_roc() run on a real small mixOmics fit: real stability categories and a real, better-than-chance pooled OOF AUC", {
  skip_if_not_installed("mixOmics")
  set.seed(400)
  n <- 30
  y <- factor(rep(c("HC", "RA"), each = n / 2))
  ids <- paste0("S", seq_len(n))
  names(y) <- ids
  expr <- matrix(rnorm(n * 15), n, 15, dimnames = list(ids, paste0("g", 1:15)))
  expr[y == "RA", 1:4] <- expr[y == "RA", 1:4] + 2.5
  meth <- matrix(rnorm(n * 15), n, 15, dimnames = list(ids, paste0("cg", 1:15)))
  meth[y == "RA", 1:4] <- meth[y == "RA", 1:4] - 2.5
  layers <- list(Transcriptomics = expr, Methylomics = meth)

  res <- mi_diablo_run(layers, y, ids, params = list(folds = 3, nrepeat = 2, validation_mode = "manual", validation_method = "mfold"))
  expect_true(res$ok)

  sig <- mb_signature_table(res)
  expect_true(is.data.frame(sig))
  expect_true(all(c("omics", "component", "feature", "stability_category") %in% colnames(sig)))
  expect_true(all(sig$stability_category %in% c("Low stability", "Moderately stable", "Stable", "Not available (needs >1 CV repeat)")))

  roc_res <- mb_cv_roc(res$fit$X, y, res$params, seed = 1)
  expect_true(!is.null(roc_res))
  expect_true(roc_res$auc >= 0 && roc_res$auc <= 1)
  expect_true(roc_res$auc > 0.6)
})

test_that("mb_cv_roc() returns NULL (never a crash) for a non-binary outcome", {
  y3 <- factor(rep(c("A", "B", "C"), each = 4)); names(y3) <- paste0("S", 1:12)
  out <- mb_cv_roc(list(), y3, list(folds = 3))
  expect_null(out)
})

## ---- Sealed hold-out and external panel transfer (added 2026-09-06) ----

mb_toy_cohort <- function(n = 40, seed = 7, gene_ids = paste0("G", 1:30), cpg_ids = paste0("cg", sprintf("%08d", 1:30))) {
  set.seed(seed)
  ids <- sprintf("S%02d", seq_len(n))
  y <- factor(rep(c("A", "B"), length.out = n), levels = c("A", "B"))
  e <- matrix(rnorm(n * length(gene_ids)), n, length(gene_ids), dimnames = list(ids, gene_ids))
  m <- matrix(rnorm(n * length(cpg_ids)), n, length(cpg_ids), dimnames = list(ids, cpg_ids))
  e[y == "B", 1:6] <- e[y == "B", 1:6] + 2.5
  m[y == "B", 1:6] <- m[y == "B", 1:6] - 2.5
  list(layers = list(Transcriptomics = e, Methylomics = m), outcome = stats::setNames(y, ids), ids = ids)
}

test_that("mb_holdout_split() seals a stratified share of the samples and refuses splits that leave too few per class", {
  cohort <- mb_toy_cohort()
  sp <- mb_holdout_split(cohort$outcome, cohort$ids, 0.25, seed = 1)
  expect_true(sp$ok)
  expect_equal(length(intersect(sp$train_ids, sp$test_ids)), 0)
  expect_setequal(c(sp$train_ids, sp$test_ids), cohort$ids)
  expect_equal(length(sp$test_ids), 10)
  expect_true(all(table(cohort$outcome[sp$test_ids]) >= 1))
  none <- mb_holdout_split(cohort$outcome, cohort$ids, 0, seed = 1)
  expect_true(none$ok); expect_equal(length(none$test_ids), 0)
  small <- mb_toy_cohort(n = 8)
  expect_false(mb_holdout_split(small$outcome, small$ids, 0.5, seed = 1)$ok)
})

test_that("mb_holdout_evaluate() scores sealed samples with the trained model and reports an AUROC with a DeLong CI", {
  cohort <- mb_toy_cohort()
  sp <- mb_holdout_split(cohort$outcome, cohort$ids, 0.25, seed = 1)
  params <- list(ncomp_mode = "manual", ncomp = 1, keepx_mode = "manual", keepx_manual = list(Transcriptomics = 8, Methylomics = 8),
                 validation_mode = "manual", validation_method = "mfold", folds = 3, nrepeat = 1, distance = "max.dist", scale = TRUE, seed = 1)
  dr <- mi_diablo_run(cohort$layers, cohort$outcome, sp$train_ids, params)
  expect_true(dr$ok)
  ev <- mb_holdout_evaluate(dr, cohort$layers, cohort$outcome, sp$test_ids)
  expect_true(ev$ok)
  expect_equal(ev$n, length(sp$test_ids))
  expect_true(ev$auc >= 0 && ev$auc <= 1)
  expect_true(ev$ci_lo <= ev$auc && ev$auc <= ev$ci_hi)
  expect_gt(ev$auc, 0.8)
  expect_setequal(ev$scores$sample_id, sp$test_ids)
})

test_that("mi_diablo_align_keepx() recycles a short keepX to ncomp so no component is silently non-sparse", {
  X <- list(A = matrix(0, 5, 40), B = matrix(0, 5, 12))
  k <- mi_diablo_align_keepx(list(A = c(20, 10), B = c(20, 10)), X, ncomp = 3)
  expect_equal(k$A, c(20, 10, 10))
  expect_equal(k$B, c(12, 10, 10))
  expect_equal(mi_diablo_align_keepx(list(A = 5, B = 5), X, ncomp = 2)$A, c(5, 5))
})

test_that("mb_panel_transfer() refits the selected panel on the discovery samples and scores an external cohort, refusing low feature coverage", {
  disc <- mb_toy_cohort(n = 40, seed = 7)
  ext <- mb_toy_cohort(n = 30, seed = 99)
  rownames(ext$layers$Transcriptomics) <- rownames(ext$layers$Methylomics) <- names(ext$outcome) <- sprintf("E%02d", 1:30)
  params <- list(ncomp_mode = "manual", ncomp = 1, keepx_mode = "manual", keepx_manual = list(Transcriptomics = 8, Methylomics = 8),
                 validation_mode = "manual", validation_method = "mfold", folds = 3, nrepeat = 1, distance = "max.dist", scale = TRUE, seed = 1)
  dr <- mi_diablo_run(disc$layers, disc$outcome, disc$ids, params)
  expect_true(dr$ok)
  tr <- mb_panel_transfer(dr, disc$layers, disc$ids, disc$outcome, ext$layers, ext$outcome)
  expect_true(tr$ok)
  expect_equal(tr$n, 30)
  expect_equal(tr$overall_coverage, 1)
  expect_gt(tr$auc, 0.8)
  ## drop most of the panel features from the external cohort -> refused under the default coverage guard
  ext_sparse <- ext$layers
  ext_sparse$Transcriptomics <- ext_sparse$Transcriptomics[, 25:30, drop = FALSE]
  ext_sparse$Methylomics <- ext_sparse$Methylomics[, 25:30, drop = FALSE]
  bad <- mb_panel_transfer(dr, disc$layers, disc$ids, disc$outcome, ext_sparse, ext$outcome)
  expect_false(bad$ok)
  expect_match(bad$error, "panel features")
})

test_that("mb_map_external_outcome() maps chosen external levels onto the discovery classes and drops the rest", {
  meta <- data.frame(row.names = c("a", "b", "c", "d"), group = c("flare", "ID", "NO ID", "flare"))
  out <- mb_map_external_outcome(meta, "group", pos_levels = "ID", neg_levels = "flare", pos_class = "resp", neg_class = "non resp")
  expect_equal(unname(out), c("non resp", "resp", NA, "non resp"))
  expect_equal(names(out), c("a", "b", "c", "d"))
})
