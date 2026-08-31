## Module 3 (Multiomics) - Cohort Harmonization sub-module's pure functions
## (cohort_harmonization_helpers.R, used by mod_multi_overview.R): the
## per-sample master table, modality descriptors (live matrix-backed vs.
## preloaded availability-table-only), pairwise overlap, ID harmonization
## status classification, candidate-column detection, analysis-cell
## enumeration/readiness, and the leakage-safe nested-CV binary-outcome
## evaluator (real glmnet/caret/pROC computation - never a fabricated AUC).

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "multiomics", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "cohort_harmonization_helpers.R"))

## ---- ch_sample_master_table() ---------------------------------------------

test_that("ch_sample_master_table() marks Present/Missing per modality and counts modalities present", {
  id_sets <- list(Transcriptomics = c("S1", "S2", "S3"), Methylomics = c("S2", "S3", "S4"))
  out <- ch_sample_master_table(id_sets)
  row_s1 <- out[out$`Sample ID` == "S1", ]
  row_s2 <- out[out$`Sample ID` == "S2", ]
  expect_equal(row_s1$Transcriptomics, "Present")
  expect_equal(row_s1$Methylomics, "Missing")
  expect_equal(row_s1$`Modalities present`, 1L)
  expect_equal(row_s2$`Modalities present`, 2L)
  expect_equal(nrow(out), 4L)  ## union of S1-S4
})

test_that("ch_sample_master_table() attaches metadata columns by matching rownames, NA for samples absent from meta", {
  id_sets <- list(Transcriptomics = c("S1", "S2"))
  meta <- data.frame(sex = c("F", "M"), row.names = c("S1", "S3"))  ## S2 has no meta row, S3 isn't in id_sets
  out <- ch_sample_master_table(id_sets, meta)
  expect_equal(out$sex[out$`Sample ID` == "S1"], "F")
  expect_true(is.na(out$sex[out$`Sample ID` == "S2"]))
})

test_that("ch_sample_master_table() returns NULL when every id set is empty", {
  expect_null(ch_sample_master_table(list(Transcriptomics = character(0))))
})

## ---- ch_value_scale() -------------------------------------------------------

test_that("ch_value_scale() identifies beta-value (0-1) matrices", {
  mat <- matrix(runif(100, 0, 1), 10, 10)
  expect_equal(ch_value_scale(mat), "Beta-values (0-1 scale, likely)")
})

test_that("ch_value_scale() identifies standardized/z-scored matrices", {
  set.seed(1)
  mat <- matrix(rnorm(200, 0, 1), 20, 10)
  expect_equal(ch_value_scale(mat), "Standardized/z-scored (likely)")
})

test_that("ch_value_scale() identifies raw-count-like matrices (non-negative integers, large range)", {
  set.seed(2)
  mat <- matrix(rpois(200, lambda = 500), 20, 10)
  expect_equal(ch_value_scale(mat), "Raw counts (likely)")
})

test_that("ch_value_scale() returns 'Unknown' for a non-matrix, NULL, or too-small input, never a guess", {
  expect_equal(ch_value_scale(NULL), "Unknown")
  expect_equal(ch_value_scale(data.frame(a = 1:5)), "Unknown")
  expect_equal(ch_value_scale(matrix(runif(5), 1, 5)), "Unknown")  ## < 10 finite values
})

## ---- ch_modality_descriptors_live() / ch_modality_descriptors() -------------

test_that("ch_modality_descriptors_live() builds one descriptor per layer with real matrix-derived n_samples/n_features/value_scale", {
  expr <- matrix(rnorm(60), nrow = 6, ncol = 10, dimnames = list(paste0("S", 1:6), paste0("g", 1:10)))
  md <- list(layers = list(expression = expr), layer_meta = list(expression = list(omics_type = "rnaseq", processing = "Normalized")))
  out <- ch_modality_descriptors_live(md)
  expect_equal(names(out), "expression")
  expect_equal(out$expression$n_samples, 6L)
  expect_equal(out$expression$n_features, 10L)
  expect_true(out$expression$has_raw_matrix)
  expect_equal(out$expression$omics_type, "rnaseq")
  expect_equal(out$expression$processing, "Normalized")
})

test_that("ch_modality_descriptors_live() returns an empty list when no layers exist", {
  expect_equal(ch_modality_descriptors_live(list()), list())
  expect_equal(ch_modality_descriptors_live(list(layers = list())), list())
})

test_that("ch_modality_descriptors() dispatches to preloaded vs. live based on multi_dataset$source, and returns empty when inactive", {
  expect_equal(ch_modality_descriptors(NULL), list())
  expect_equal(ch_modality_descriptors(list(active = FALSE, source = "upload")), list())

  skip_if_not(MULTI_DATA_AVAILABLE, "multi-omics preloaded data not available in this deployment")
  out_preloaded <- ch_modality_descriptors(list(active = TRUE, source = "preloaded"))
  expect_true(all(c("Transcriptomics", "Methylomics") %in% names(out_preloaded)))
  expect_false(out_preloaded$Transcriptomics$has_raw_matrix)  ## no raw matrix bundled for preloaded
})

## ---- ch_modality_descriptors_preloaded() (real registry tables) ------------

test_that("ch_modality_descriptors_preloaded() reports real sample counts from the actual patient matching table", {
  skip_if_not(MULTI_DATA_AVAILABLE, "multi-omics preloaded data not available in this deployment")
  matching <- multi_read_registry_table("Patient sample matching (all 80 patients)")
  skip_if_not(matching$ok, "patient matching table not available")
  expected_rna_n <- sum(matching$df$RNA_available_PBMC %in% c(TRUE, "TRUE", "Yes", "yes", 1))
  expected_meth_n <- sum(matching$df$methylation_available %in% c(TRUE, "TRUE", "Yes", "yes", 1))

  out <- ch_modality_descriptors_preloaded()
  expect_equal(out$Transcriptomics$n_samples, expected_rna_n)
  expect_equal(out$Methylomics$n_samples, expected_meth_n)
  expect_true(grepl("no raw matrix bundled", out$Transcriptomics$value_scale))
})

## ---- ch_pairwise_overlap_matrix() -------------------------------------------

test_that("ch_pairwise_overlap_matrix() computes a symmetric NxN intersection-count matrix with correct diagonal", {
  sets <- list(A = c("S1", "S2", "S3"), B = c("S2", "S3", "S4"), C = c("S3"))
  m <- ch_pairwise_overlap_matrix(sets)
  expect_equal(m["A", "A"], 3L)  ## self-overlap = own size
  expect_equal(m["A", "B"], 2L)
  expect_equal(m["B", "A"], 2L)  ## symmetric
  expect_equal(m["A", "C"], 1L)
  expect_equal(m["C", "C"], 1L)
})

test_that("ch_pairwise_overlap_matrix() returns NULL when every set is NULL/absent", {
  expect_null(ch_pairwise_overlap_matrix(list(A = NULL, B = NULL)))
})

## ---- ch_id_harmonization_table() --------------------------------------------

test_that("ch_id_harmonization_table() marks 'Exact match' for identical IDs shared across >=2 modalities", {
  sets <- list(RNA = c("S1", "S2"), Meth = c("S1", "S3"))
  out <- ch_id_harmonization_table(sets)
  s1_rows <- out[out$Original == "S1", ]
  expect_true(all(s1_rows$Status == "Exact match"))
})

test_that("ch_id_harmonization_table() marks 'Normalized match' when only case/whitespace differ across modalities", {
  sets <- list(RNA = "  Sample1", Meth = "sample1")
  out <- ch_id_harmonization_table(sets)
  expect_true(all(out$Status == "Normalized match"))
})

test_that("ch_id_harmonization_table() marks 'Unmatched' for an ID present in only one modality", {
  sets <- list(RNA = c("S1", "S2"), Meth = c("S1"))
  out <- ch_id_harmonization_table(sets)
  expect_equal(out$Status[out$Modality == "RNA" & out$Original == "S2"], "Unmatched")
})

test_that("ch_id_harmonization_table() marks 'Duplicate' (repeated within one modality only) vs. 'Ambiguous' (also appears elsewhere)", {
  dup_only <- ch_id_harmonization_table(list(RNA = c("S1", "S1", "S2")))
  expect_equal(unique(dup_only$Status[dup_only$Original == "S1"]), "Duplicate")

  ambiguous <- ch_id_harmonization_table(list(RNA = c("S1", "S1"), Meth = c("S1")))
  expect_equal(unique(ambiguous$Status[ambiguous$Modality == "RNA"]), "Ambiguous")
})

## FIXED (was a KNOWN BUG): blank/empty identifiers were silently left with
## an empty Status instead of "Invalid", because `by_norm[[nrm]]` looked
## the group up via `[[""]]` - base R's `[[` always returns NULL for a
## zero-length-string name lookup, even when a list element is genuinely
## named "". ch_id_harmonization_table() now iterates by_norm by position
## instead of by name, sidestepping that lookup entirely.
test_that("ch_id_harmonization_table() marks blank/empty identifiers 'Invalid' (fixed - was silently left blank)", {
  out <- ch_id_harmonization_table(list(RNA = c("", "  ", "S2")))
  blank_rows <- out[out$Original %in% c("", "  "), ]
  expect_true(all(blank_rows$Status == "Invalid"))
  expect_true(all(blank_rows$Reason == "Empty or missing identifier."))
})

test_that("ch_id_harmonization_table() returns NULL when every modality's id set is NULL", {
  expect_null(ch_id_harmonization_table(list(RNA = NULL)))
})

## ---- ch_classify_column() / ch_classify_metadata_columns() ------------------

test_that("ch_classify_column() flags near-all-unique columns as 'identifier'", {
  expect_equal(ch_classify_column(paste0("ID_", 1:50)), "identifier")
})

test_that("ch_classify_column() flags high-cardinality numeric columns as 'continuous', but only below the 'near-all-unique' identifier threshold", {
  ## A column where every value is unique (e.g. rnorm() with no repeats) is
  ## classified "identifier" instead - the near-all-unique check (>90% of
  ## values unique) is checked first and fires regardless of numeric type.
  ## Realistic bounded-range continuous data (e.g. integer ages) has real
  ## repeats, so it correctly falls through to "continuous".
  set.seed(5)
  v <- sample(20:80, 50, replace = TRUE)
  expect_equal(ch_classify_column(v), "continuous")
})

test_that("ch_classify_column() classifies a numeric column with no repeated values as 'identifier', not 'continuous' (near-all-unique wins first)", {
  set.seed(3)
  expect_equal(ch_classify_column(rnorm(50)), "identifier")
})

test_that("ch_classify_column() defaults to 'categorical' for low-cardinality columns", {
  expect_equal(ch_classify_column(rep(c("A", "B"), 25)), "categorical")
  expect_equal(ch_classify_column(character(0)), "categorical")  ## all-empty -> categorical, not an error
})

test_that("ch_classify_metadata_columns() suggests a keyword-matched column over a plain balanced categorical one", {
  meta <- data.frame(
    patient_id = paste0("P", 1:20),
    batch = rep(c("A", "B"), 10),
    treatment_response = rep(c("responder", "non-responder"), 10),
    stringsAsFactors = FALSE
  )
  out <- ch_classify_metadata_columns(meta)
  expect_equal(out$table$type[out$table$column == "patient_id"], "identifier")
  expect_equal(out$suggested_default, "treatment_response")  ## keyword hit beats plain "batch"
})

test_that("ch_classify_metadata_columns() falls back to the first balanced categorical column when no keyword column exists", {
  meta <- data.frame(patient_id = paste0("P", 1:20), site = rep(c("X", "Y"), 10), stringsAsFactors = FALSE)
  out <- ch_classify_metadata_columns(meta)
  expect_equal(out$suggested_default, "site")
})

test_that("ch_classify_metadata_columns() still falls back to the only categorical column even when it's heavily imbalanced (>90% one level) - a starting pick, not an eligibility gate", {
  meta <- data.frame(patient_id = paste0("P", 1:20), rare_flag = c(rep("no", 19), "yes"), stringsAsFactors = FALSE)
  out <- ch_classify_metadata_columns(meta)
  expect_equal(out$suggested_default, "rare_flag")
})

test_that("ch_classify_metadata_columns() returns NULL suggested_default only when there is no categorical column at all", {
  meta <- data.frame(patient_id = paste0("P", 1:20), age = rnorm(20), stringsAsFactors = FALSE)
  out <- ch_classify_metadata_columns(meta)
  expect_null(out$suggested_default)
})

test_that("ch_classify_metadata_columns() returns NULL table/suggestion for empty metadata", {
  out <- ch_classify_metadata_columns(NULL)
  expect_null(out$table)
  expect_null(out$suggested_default)
})

## ---- ch_detect_candidate_columns() -------------------------------------------

test_that("ch_detect_candidate_columns('batch') matches batch/cohort/study/platform/site-named columns and excludes ID-like ones", {
  meta <- data.frame(sample_id = 1:5, batch_id = 1:5, study_site = 1:5, age = 1:5)
  out <- ch_detect_candidate_columns(meta, "batch")
  expect_setequal(out, c("batch_id", "study_site"))
})

test_that("ch_detect_candidate_columns('phenotype') returns only non-ID categorical columns", {
  meta <- data.frame(patient_id = paste0("P", 1:10), sex = rep(c("F", "M"), 5), age = rnorm(10))
  out <- ch_detect_candidate_columns(meta, "phenotype")
  expect_equal(out, "sex")  ## age is continuous, patient_id is id-like
})

## ---- ch_matched_sample_summary() ---------------------------------------------

test_that("ch_matched_sample_summary() reports 'Matched' when the intersection equals the smaller modality's full size", {
  out <- ch_matched_sample_summary(list(RNA = c("S1", "S2", "S3"), Meth = c("S1", "S2")))
  expect_equal(out$n_matched, 2L)
  expect_equal(out$status, "Matched")
  expect_true(grepl("2 matched", out$sentence))
})

test_that("ch_matched_sample_summary() reports 'Partially matched' when overlap is a strict subset of every modality", {
  out <- ch_matched_sample_summary(list(RNA = c("S1", "S2", "S3"), Meth = c("S2", "S3", "S4")))
  expect_equal(out$n_matched, 2L)
  expect_equal(out$status, "Partially matched")
})

test_that("ch_matched_sample_summary() reports 'Unmatched' with zero overlap, and 'Single modality' with only one", {
  unmatched <- ch_matched_sample_summary(list(RNA = c("S1"), Meth = c("S2")))
  expect_equal(unmatched$status, "Unmatched")
  single <- ch_matched_sample_summary(list(RNA = c("S1", "S2")))
  expect_equal(single$status, "Single modality")
})

test_that("ch_matched_sample_summary() handles zero modalities without fabricating a match", {
  out <- ch_matched_sample_summary(list())
  expect_equal(out$status, "Unmatched")
  expect_equal(out$n_matched, 0L)
})

## ---- ch_analysis_cells() -------------------------------------------------------

test_that("ch_analysis_cells() enumerates every subset (singles + pairs + full) for <= max_full_subsets modalities", {
  sets <- list(RNA = paste0("S", 1:10), Meth = paste0("S", 6:15))
  out <- ch_analysis_cells(sets, pheno_available = TRUE, min_integration = 3, min_prediction = 6)
  expect_null(out$omitted_note)
  labels <- vapply(out$cells, `[[`, character(1), "label")
  expect_setequal(labels, c("RNA", "Meth", "RNA + Meth"))
  fused <- out$cells[[which(labels == "RNA + Meth")]]
  expect_equal(fused$n_matched, 5L)  ## S6-S10
  expect_true("Unsupervised integration" %in% fused$methods)  ## 5 >= min_integration(3)
  expect_false("Supervised prediction" %in% fused$methods)  ## 5 < min_prediction(6)
})

test_that("ch_analysis_cells() reports 'None feasible' when a cell has too few matched samples for any method", {
  sets <- list(RNA = c("S1", "S2"), Meth = c("S3", "S4"))  ## zero overlap
  out <- ch_analysis_cells(sets, pheno_available = TRUE)
  fused <- out$cells[[which(vapply(out$cells, `[[`, character(1), "label") == "RNA + Meth")]]
  expect_equal(fused$methods, "None feasible")
})

test_that("ch_analysis_cells() caps combinations and states an explicit omitted_note beyond max_full_subsets modalities", {
  sets <- setNames(lapply(1:5, function(i) paste0("S", 1:10)), paste0("M", 1:5))
  out <- ch_analysis_cells(sets, max_full_subsets = 4)
  expect_true(!is.null(out$omitted_note))
  expect_true(grepl("omitted for brevity", out$omitted_note))
  ## Singles (5) + pairs (choose(5,2)=10) + the one full combination = 16, not the full 2^5-1=31.
  expect_equal(length(out$cells), 16L)
})

test_that("ch_analysis_cells() returns an empty cell list for zero modalities", {
  out <- ch_analysis_cells(list())
  expect_equal(out$cells, list())
  expect_null(out$omitted_note)
})

## ---- ch_integration_readiness() -------------------------------------------------

test_that("ch_integration_readiness() classifies single-modality cells as 'single' regardless of sample size", {
  cell <- list(modalities = "RNA", n_matched = 100, label = "RNA")
  out <- ch_integration_readiness(cell)
  expect_equal(out$level, "single")
})

test_that("ch_integration_readiness() classifies Ready/Limited/Not suitable by matched-sample count thresholds", {
  base_cell <- list(modalities = c("RNA", "Meth"), label = "RNA + Meth")
  ready <- ch_integration_readiness(c(base_cell, list(n_matched = 15)), min_ready = 10, min_limited = 3)
  limited <- ch_integration_readiness(c(base_cell, list(n_matched = 5)), min_ready = 10, min_limited = 3)
  not_suitable <- ch_integration_readiness(c(base_cell, list(n_matched = 2)), min_ready = 10, min_limited = 3)
  expect_equal(ready$level, "ready")
  expect_equal(limited$level, "limited")
  expect_equal(not_suitable$level, "not_suitable")
})

## ---- ch_fold_predict_view() / ch_evaluate_binary_outcome() (real glmnet/caret/pROC) ----

test_that("ch_fold_predict_view() falls back to 0.5 (no fabricated confidence) when the training fold has a single outcome class", {
  set.seed(10)
  X <- matrix(rnorm(100), 10, 10)
  y <- factor(rep("A", 10), levels = c("A", "B"))  ## single class in the whole vector
  out <- ch_fold_predict_view(X, y, train_idx = 1:8, test_idx = 9:10, max_features = 5)
  expect_equal(out, rep(0.5, 2))
})

test_that("ch_evaluate_binary_outcome() (real glmnet/caret nested CV) recovers strong AUC on a clearly separable synthetic signal, beating the majority baseline", {
  set.seed(42)
  n <- 60
  y <- factor(rep(c("A", "B"), each = n / 2))
  ## A real, moderately separable signal - not perfectly separable, so glmnet's
  ## regularization has genuine work to do rather than trivially memorizing.
  X <- matrix(rnorm(n * 30), n, 30, dimnames = list(paste0("S", 1:n), paste0("f", 1:30)))
  X[y == "B", 1:5] <- X[y == "B", 1:5] + 1.5
  rownames(X) <- paste0("S", 1:n)
  names(y) <- paste0("S", 1:n)

  out <- ch_evaluate_binary_outcome(list(expression = X), y, k_folds = 5, seed = 1)
  expect_true(out$ok)
  expect_equal(out$n, n)
  expect_true(out$per_view_auc$expression > 0.6)  ## real signal should beat chance clearly
  expect_equal(out$majority_baseline, 0.5)
})

test_that("ch_evaluate_binary_outcome() refuses (never silently downgrades) with fewer than 10 matched samples", {
  y <- factor(rep(c("A", "B"), each = 3)); names(y) <- paste0("S", 1:6)
  X <- matrix(rnorm(60), 6, 10, dimnames = list(paste0("S", 1:6), paste0("f", 1:10)))
  out <- ch_evaluate_binary_outcome(list(expression = X), y)
  expect_false(out$ok)
  expect_true(grepl("at least 10", out$error))
})

test_that("ch_evaluate_binary_outcome() refuses a non-binary outcome and a too-small minority class", {
  y3 <- factor(rep(c("A", "B", "C"), each = 5)); names(y3) <- paste0("S", 1:15)
  X15 <- matrix(rnorm(150), 15, 10, dimnames = list(paste0("S", 1:15), paste0("f", 1:10)))
  out_multi <- ch_evaluate_binary_outcome(list(expression = X15), y3)
  expect_false(out_multi$ok)
  expect_true(grepl("only binary outcomes", out_multi$error))

  y_imb <- factor(c(rep("A", 18), rep("B", 2))); names(y_imb) <- paste0("S", 1:20)
  X20 <- matrix(rnorm(200), 20, 10, dimnames = list(paste0("S", 1:20), paste0("f", 1:10)))
  out_imb <- ch_evaluate_binary_outcome(list(expression = X20), y_imb)
  expect_false(out_imb$ok)
  expect_true(grepl("at least 3 per class", out_imb$error))
})

## FIXED (was a KNOWN BUG): the function's own documented guard ("Duplicate
## sample IDs detected across the matched samples - resolve before
## evaluation.") was unreachable dead code, since `intersect()` always
## de-duplicates its own output before the old `any(duplicated(common))`
## check could ever see a duplicate. The check now looks for duplicates in
## each input's own raw rownames/names (restricted to the matched set)
## directly, before any intersect() call can hide them.
test_that("ch_evaluate_binary_outcome() refuses when sample IDs are duplicated (fixed - was silently proceeding)", {
  ids <- c("S1", rep("S2", 5), paste0("S", 3:15))  ## S2 duplicated 5x, 15 rows total
  y <- factor(rep(c("A", "B"), length.out = length(ids))); names(y) <- ids
  X <- matrix(rnorm(length(ids) * 10), length(ids), 10, dimnames = list(ids, paste0("f", 1:10)))
  out <- ch_evaluate_binary_outcome(list(expression = X), y)
  expect_false(out$ok)
  expect_true(grepl("Duplicate sample IDs", out$error))
})

test_that("ch_evaluate_binary_outcome() still runs normally on genuinely unique sample IDs after the duplicate-check fix", {
  set.seed(2200)
  n <- 20
  ids <- paste0("S", seq_len(n))
  y <- factor(rep(c("A", "B"), each = n / 2)); names(y) <- ids
  X <- matrix(rnorm(n * 10), n, 10, dimnames = list(ids, paste0("f", 1:10)))
  X[y == "B", 1:3] <- X[y == "B", 1:3] + 3
  out <- ch_evaluate_binary_outcome(list(expression = X), y)
  expect_true(out$ok)
  expect_equal(out$n, n)
})
