## Module 3 (Multiomics) - the "Live Analysis (Upload & MOFA2)" engine's
## core pure functions (multiomics_dataset_helpers.R): upload parsing (wide +
## long/tidy pivot), orientation/omics-type detection, matrix validation,
## sample overlap/missingness, explicit (never-automatic) missing-data
## handling/normalization/filtering/scaling, PCA/confounding/batch-
## correction diagnostics (real prcomp/sva::ComBat/limma::removeBatchEffect),
## cross-omics correlation, MOFA2 guardrails, and dataset-compatibility
## verdicts. This is the ONE part of Multiomics that computes on
## user-supplied data rather than browsing a precomputed cohort, so every
## claim here is checked against real computation, never invented.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "multiomics", "06_Gene_CpG_Concordance", "multiomics_concordance_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "multiomics_dataset_helpers.R"))

## ---- multi_live_read_matrix() ---------------------------------------------

test_that("multi_live_read_matrix() reads a samples-in-rows CSV correctly and counts coerced-NA cells", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("sample,f1,f2,f3", "S1,1.5,2.5,x", "S2,3.5,4.5,5.5"), path)
  out <- multi_live_read_matrix(path, orientation = "samples_rows")
  expect_true(out$ok)
  expect_equal(dim(out$mat), c(2L, 3L))
  expect_equal(rownames(out$mat), c("S1", "S2"))
  expect_equal(out$n_coerced_na, 1L)  ## "x" -> NA
})

test_that("multi_live_read_matrix() transposes when orientation='features_rows'", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("feature,S1,S2,S3", "f1,1,2,3", "f2,4,5,6"), path)
  out <- multi_live_read_matrix(path, orientation = "features_rows")
  expect_equal(dim(out$mat), c(3L, 2L))  ## transposed to samples x features
  expect_setequal(rownames(out$mat), c("S1", "S2", "S3"))
})

test_that("multi_live_read_matrix() fails soft for a missing file or a file with fewer than 2 columns", {
  out_missing <- multi_live_read_matrix(file.path(tempdir(), "does-not-exist.csv"))
  expect_false(out_missing$ok)
  path <- tempfile(fileext = ".csv")
  writeLines(c("sample", "S1", "S2"), path)
  out_onecol <- multi_live_read_matrix(path)
  expect_false(out_onecol$ok)
  expect_true(grepl("at least an ID column", out_onecol$error))
})

## ---- multi_live_detect_table_shape() / detect_long_columns() / pivot_long() ----

test_that("multi_live_detect_table_shape() confidently detects a long/tidy table (many rows, few numeric columns, repeated IDs)", {
  df <- data.frame(sample = rep(paste0("S", 1:5), each = 20), gene = rep(paste0("g", 1:20), 5), value = rnorm(100))
  out <- multi_live_detect_table_shape(df)
  expect_equal(out$shape, "long")
  expect_true(out$confident)
})

test_that("multi_live_detect_table_shape() calls a small/near-square, mostly-numeric table 'wide' (not confidently, since that's the default)", {
  df <- data.frame(id = paste0("g", 1:5), S1 = rnorm(5), S2 = rnorm(5), S3 = rnorm(5))
  out <- multi_live_detect_table_shape(df)
  expect_equal(out$shape, "wide")
})

test_that("multi_live_detect_long_columns() finds feature/sample/value columns by name, warns when a role is undetected", {
  df <- data.frame(gene_symbol = c("A", "B"), sample_id = c("S1", "S2"), tpm = c(1.1, 2.2), condition = c("HC", "RA"))
  out <- multi_live_detect_long_columns(df)
  expect_equal(out$feature_col, "gene_symbol")
  expect_equal(out$sample_col, "sample_id")
  expect_equal(out$value_col, "tpm")
  expect_equal(out$group_col, "condition")
  expect_equal(length(out$warnings), 0L)
})

test_that("multi_live_detect_long_columns() falls back to the first remaining numeric column for value, and warns when feature/sample can't be named", {
  df <- data.frame(colA = c("x1", "x2"), colB = c("y1", "y2"), measurement_score = c(1.1, 2.2))
  out <- multi_live_detect_long_columns(df)
  expect_equal(out$value_col, "measurement_score")
  expect_true(is.na(out$feature_col))
  expect_true(is.na(out$sample_col))
  expect_equal(length(out$warnings), 2L)
})

test_that("multi_live_pivot_long() pivots into the correct wide matrix, averaging duplicate (sample,feature) pairs and reporting the count", {
  df <- data.frame(
    feature = c("g1", "g1", "g2", "g1", "g2"),
    sample = c("S1", "S1", "S1", "S2", "S2"),
    value = c(1, 3, 5, 7, 9)  ## (g1,S1) duplicated: 1 and 3 -> averaged to 2
  )
  out <- multi_live_pivot_long(df, "feature", "sample", "value")
  expect_true(out$ok)
  expect_equal(out$mat["S1", "g1"], 2)
  expect_equal(out$mat["S1", "g2"], 5)
  expect_equal(out$mat["S2", "g1"], 7)
  expect_equal(out$n_duplicate_pairs, 1L)
})

test_that("multi_live_pivot_long() builds group_df from a per-sample-consistent group column, and warns (without df) when it varies within a sample", {
  df_ok <- data.frame(feature = c("g1", "g2"), sample = c("S1", "S1"), value = c(1, 2), grp = c("RA", "RA"))
  out_ok <- multi_live_pivot_long(df_ok, "feature", "sample", "value", group_col = "grp")
  expect_equal(out_ok$group_df$df["S1", "grp"], "RA")

  df_bad <- data.frame(feature = c("g1", "g2"), sample = c("S1", "S1"), value = c(1, 2), grp = c("RA", "HC"))
  out_bad <- multi_live_pivot_long(df_bad, "feature", "sample", "value", group_col = "grp")
  expect_null(out_bad$group_df$df)
  expect_true(nzchar(out_bad$group_df$warning))
})

test_that("multi_live_pivot_long() fails soft when required columns are unselected or the value column has no numeric data", {
  df <- data.frame(feature = "g1", sample = "S1", value = "not_a_number")
  out_missing_cols <- multi_live_pivot_long(df, NA, "sample", "value")
  expect_false(out_missing_cols$ok)
  out_no_numeric <- multi_live_pivot_long(df, "feature", "sample", "value")
  expect_false(out_no_numeric$ok)
  expect_true(grepl("no valid numeric values", out_no_numeric$error))
})

## ---- multi_live_detect_orientation() ---------------------------------------

test_that("multi_live_detect_orientation() confidently suggests 'features_rows' when columns vastly outnumber rows, headers look like sample IDs, and the first column does NOT (e.g. descriptive gene names)", {
  n_samples <- 40; n_features <- 5
  df <- as.data.frame(matrix(rnorm(n_features * n_samples), n_features, n_samples))
  colnames(df) <- paste0("Sample_", 1:n_samples)
  ## A first column of descriptive names (spaces) fails the same
  ## alnum-only regex the header IDs pass - the function distinguishes
  ## "features_rows" from "samples_rows" precisely on this asymmetry, not
  ## on row/column count alone.
  df <- cbind(feature = paste("Gene description", 1:n_features), df)
  out <- multi_live_detect_orientation(df)
  expect_equal(out$suggested, "features_rows")
  expect_true(out$confident)
})

test_that("multi_live_detect_orientation() defaults to (unconfident) 'samples_rows' for an ordinary wide table", {
  df <- data.frame(sample = paste0("S", 1:10), f1 = rnorm(10), f2 = rnorm(10))
  out <- multi_live_detect_orientation(df)
  expect_equal(out$suggested, "samples_rows")
  expect_false(out$confident)
})

## ---- multi_live_validate_matrix() ------------------------------------------

test_that("multi_live_validate_matrix() reports real missingness/zero-variance/duplicate counts from the actual matrix", {
  mat <- matrix(c(1, 2, NA, 4, 5, 5, 5, 5, 7, 8, 9, 10), nrow = 3, dimnames = list(c("S1", "S2", "S1"), c("f1", "f2", "f3", "f3")))
  out <- multi_live_validate_matrix(mat, "TestLayer")
  expect_true(out$ok)
  expect_equal(out$n_samples, 3L)
  expect_equal(out$n_features, 4L)
  expect_equal(out$n_missing, 1L)
  expect_equal(out$n_duplicate_samples, 1L)  ## "S1" appears twice
  expect_equal(out$n_duplicate_features, 1L)  ## "f3" appears twice
})

test_that("multi_live_validate_matrix() fails soft for a NULL or non-matrix input", {
  expect_false(multi_live_validate_matrix(NULL)$ok)
  expect_false(multi_live_validate_matrix(data.frame(a = 1))$ok)
})

## ---- multi_live_detect_omics_type() (real mcc_detect_id_type/value_type) ---

test_that("multi_live_detect_omics_type() identifies methylation from Illumina CpG probe IDs with beta-scale values", {
  mat <- matrix(runif(30, 0, 1), 5, 6, dimnames = list(paste0("S", 1:5), paste0("cg", 10000000 + 1:6)))
  out <- multi_live_detect_omics_type(mat)
  expect_equal(out$detected, "methylation")
  expect_true(out$corroborated)
})

test_that("multi_live_detect_omics_type() flags a mismatch (not corroborated) when methylation-looking IDs carry non-beta-scale values", {
  ## mean 50 / sd 3 keeps every value far above 0 (>1.001, so never "beta")
  ## and never negative (so never "M-value" either) with no realistic
  ## chance of flakiness - unlike rnorm(30, 5, 2), which is only ~2.5 SD
  ## from 0 and occasionally dips negative, accidentally classifying as a
  ## valid M-value and making this test pass for the wrong reason.
  mat <- matrix(rnorm(30, 50, 3), 5, 6, dimnames = list(paste0("S", 1:5), paste0("cg", 10000000 + 1:6)))
  out <- multi_live_detect_omics_type(mat)
  expect_equal(out$detected, "methylation")
  expect_false(out$corroborated)
})

test_that("multi_live_detect_omics_type() identifies rnaseq from gene-symbol-like feature IDs", {
  mat <- matrix(rpois(30, 100), 5, 6, dimnames = list(paste0("S", 1:5), c("TP53", "BRCA1", "EGFR", "MYC", "KRAS", "PTEN")))
  out <- multi_live_detect_omics_type(mat)
  expect_equal(out$detected, "rnaseq")
})

## ---- multi_live_qc_summary_table() / multi_live_sample_overlap() ----------

test_that("multi_live_qc_summary_table() tabulates only the ok=TRUE validations, real per-layer counts", {
  validations <- list(
    expr = multi_live_validate_matrix(matrix(rnorm(20), 4, 5, dimnames = list(paste0("S", 1:4), paste0("f", 1:5))), "expr"),
    bad = list(ok = FALSE)
  )
  tbl <- multi_live_qc_summary_table(validations)
  expect_equal(nrow(tbl), 1L)
  expect_equal(tbl$Omics, "expr")
})

test_that("multi_live_sample_overlap() computes real shared/layer-only sample sets across >=2 layers", {
  m1 <- matrix(1, 4, 2, dimnames = list(c("S1", "S2", "S3", "S4"), c("f1", "f2")))
  m2 <- matrix(1, 3, 2, dimnames = list(c("S2", "S3", "S5"), c("g1", "g2")))
  out <- multi_live_sample_overlap(list(A = m1, B = m2))
  expect_true(out$ok)
  expect_setequal(out$shared_ids, c("S2", "S3"))
  expect_equal(out$n_shared, 2L)
  expect_setequal(out$layer_only[[1]], c("S1", "S4"))
})

test_that("multi_live_sample_overlap() refuses with fewer than 2 layers", {
  out <- multi_live_sample_overlap(list(A = matrix(1, 2, 2)))
  expect_false(out$ok)
})

## ---- multi_live_missingness() / multi_live_handle_missing() ---------------

test_that("multi_live_missingness() computes real per-sample and per-feature missingness percentages", {
  ## Column-major fill: f1=(1,NA), f2=(3,NA), f3=(NA,6).
  mat <- matrix(c(1, NA, 3, NA, NA, 6), 2, 3, dimnames = list(c("S1", "S2"), c("f1", "f2", "f3")))
  out <- multi_live_missingness(mat)
  expect_equal(out$per_sample$pct_missing[out$per_sample$sample == "S1"], 100 / 3, tolerance = 1e-6)
  expect_equal(out$per_feature$pct_missing[out$per_feature$feature == "f2"], 50)
})

test_that("multi_live_handle_missing() drops samples/features over the missingness threshold, then imputes/removes as explicitly requested", {
  ## f2 is 60% missing (3/5 rows) -> dropped by the default 50% feature
  ## threshold. f3's single NA (row S4) sits in a DIFFERENT row than f2's
  ## NAs, so every row's own missingness stays at/under 33% (1 of 3
  ## columns) - no row gets dropped, and f3's NA survives the threshold
  ## step for "none" to leave in place and "mean" to actually impute.
  mat <- rbind(
    S1 = c(f1 = 1,  f2 = NA, f3 = 9),
    S2 = c(f1 = 2,  f2 = NA, f3 = 10),
    S3 = c(f1 = 3,  f2 = NA, f3 = 11),
    S4 = c(f1 = 4,  f2 = 8,  f3 = NA),
    S5 = c(f1 = 5,  f2 = 9,  f3 = 13)
  )
  out_mean <- multi_live_handle_missing(mat, method = "mean")
  expect_true(out_mean$ok)
  expect_false("f2" %in% colnames(out_mean$mat))
  expect_equal(nrow(out_mean$mat), 5L)  ## no row exceeds the 50% row threshold
  expect_equal(out_mean$n_remaining_na, 0L)  ## mean-imputed

  out_none <- multi_live_handle_missing(mat, method = "none")
  expect_equal(out_none$n_remaining_na, 1L)  ## "none" leaves f3's one remaining NA as-is
})

test_that("multi_live_handle_missing() refuses when every sample/feature is dropped by the thresholds", {
  mat <- matrix(NA_real_, 3, 3, dimnames = list(paste0("S", 1:3), paste0("f", 1:3)))
  out <- multi_live_handle_missing(mat, method = "none")
  expect_false(out$ok)
})

## ---- multi_live_normalize() / multi_live_filter_features() / multi_live_scale() ----

test_that("multi_live_normalize() log2 transform matches a hand-computed value, and 'none'/unknown falls through to the input matrix unchanged", {
  mat <- matrix(c(0, 1, 3, 7), 2, 2)
  out <- multi_live_normalize(mat, "rnaseq", "log2")
  expect_equal(out$mat[1, 1], log2(0 + 1))
  expect_equal(out$mat[2, 2], log2(7 + 1))

  out_none <- multi_live_normalize(mat, "rnaseq", "none")
  expect_identical(out_none$mat, mat)
})

test_that("multi_live_normalize() mvalue transform is the correct logit of a beta matrix", {
  mat <- matrix(c(0.1, 0.5, 0.9, 0.99999), 2, 2)  ## last value clamped to 1-1e-3 before logit
  out <- multi_live_normalize(mat, "methylation", "mvalue")
  expect_equal(out$mat[1, 1], log2(0.1 / 0.9), tolerance = 1e-8)
  b_clamped <- 1 - 1e-3
  expect_equal(out$mat[2, 2], log2(b_clamped / (1 - b_clamped)), tolerance = 1e-8)
})

test_that("multi_live_filter_features() keeps exactly the top-N highest-variance features and reports before/after counts", {
  set.seed(50)
  mat <- matrix(rnorm(50), 5, 10)
  colnames(mat) <- paste0("f", 1:10)
  vars <- apply(mat, 2, var)
  top3 <- names(sort(vars, decreasing = TRUE))[1:3]
  out <- multi_live_filter_features(mat, criterion = "variance", keep_top_n = 3)
  expect_setequal(colnames(out$mat), top3)
  expect_equal(out$n_before, 10L)
  expect_equal(out$n_after, 3L)
  expect_equal(out$n_removed, 7L)
})

test_that("multi_live_scale() z-scores every column to mean 0, sd 1", {
  mat <- matrix(rnorm(40, mean = 10, sd = 5), 10, 4)
  out <- multi_live_scale(mat)
  expect_equal(unname(colMeans(out)), rep(0, 4), tolerance = 1e-10)
  expect_equal(unname(apply(out, 2, sd)), rep(1, 4), tolerance = 1e-10)
})

## ---- multi_live_pca() (real stats::prcomp) ---------------------------------

test_that("multi_live_pca() runs real PCA and reports variance-explained fractions summing to 1", {
  set.seed(60)
  mat <- matrix(rnorm(100), 10, 10, dimnames = list(paste0("S", 1:10), paste0("f", 1:10)))
  out <- multi_live_pca(mat)
  expect_true(out$ok)
  expect_equal(nrow(out$scores), 10L)
  expect_equal(sum(out$var_explained), 1, tolerance = 1e-8)
})

test_that("multi_live_pca() refuses with fewer than 3 samples or 2 features, or all-zero-variance features", {
  expect_false(multi_live_pca(matrix(1, 2, 5))$ok)
  expect_false(multi_live_pca(matrix(1, 5, 1))$ok)
  zero_var <- matrix(5, 5, 5)  ## every column constant -> zero variance
  expect_false(multi_live_pca(zero_var)$ok)
})

## ---- multi_live_confounding_check() / multi_live_variance_by_group() ------

test_that("multi_live_confounding_check() flags complete confounding when a batch level maps to exactly one phenotype level", {
  meta <- data.frame(batch = c("B1", "B1", "B1", "B2", "B2", "B2"), phenotype = c("HC", "HC", "HC", "RA", "RA", "RA"))
  out <- multi_live_confounding_check(meta, "batch", "phenotype")
  expect_true(out$confounded)
})

test_that("multi_live_confounding_check() does not flag confounding when every batch contains both phenotype levels", {
  ## Each batch level must contain >=2 distinct phenotype values to avoid
  ## the "one batch level maps to exactly one phenotype level" flag -
  ## rep(c("B1","B2"),6) paired positionally with rep(c("HC","RA"),6) is
  ## itself perfectly confounded (every B1 row is HC, every B2 row is RA);
  ## an explicit within-batch phenotype mix is needed instead.
  meta <- data.frame(batch = rep(c("B1", "B2"), each = 6), phenotype = rep(c("HC", "RA"), 6))
  out <- multi_live_confounding_check(meta, "batch", "phenotype")
  expect_false(out$confounded)
})

test_that("multi_live_variance_by_group() computes a real one-way ANOVA R^2 of PCs against a grouping variable", {
  set.seed(70)
  grp <- rep(c("A", "B"), each = 10)
  pc1 <- c(rnorm(10, -3), rnorm(10, 3))  ## strongly separated by group
  scores <- data.frame(PC1 = pc1, PC2 = rnorm(20), row.names = paste0("S", 1:20))
  meta <- data.frame(group = grp, row.names = paste0("S", 1:20))
  out <- multi_live_variance_by_group(scores, meta, "group", npcs = 2)
  expect_true(out[1] > 0.5)  ## PC1 strongly explained by group
})

## ---- multi_live_batch_correct() (real sva::ComBat / limma::removeBatchEffect) ----

test_that("multi_live_batch_correct() (real ComBat) removes a known additive batch shift while preserving matrix dimensions", {
  set.seed(80)
  n <- 20
  mat <- matrix(rnorm(n * 30), n, 30, dimnames = list(paste0("S", 1:n), paste0("f", 1:30)))
  batch <- rep(c("B1", "B2"), each = n / 2)
  mat[batch == "B2", ] <- mat[batch == "B2", ] + 5  ## known additive shift
  out <- multi_live_batch_correct(mat, batch, method = "combat")
  expect_true(out$ok)
  expect_equal(dim(out$mat), dim(mat))
  ## After correction, the two batches' column means should be far closer than before.
  before_gap <- mean(abs(colMeans(mat[batch == "B1", ]) - colMeans(mat[batch == "B2", ])))
  after_gap <- mean(abs(colMeans(out$mat[batch == "B1", ]) - colMeans(out$mat[batch == "B2", ])))
  expect_true(after_gap < before_gap * 0.5)
})

test_that("multi_live_batch_correct() refuses with a single-level batch vector", {
  mat <- matrix(rnorm(20), 4, 5)
  out <- multi_live_batch_correct(mat, rep("B1", 4), method = "combat")
  expect_false(out$ok)
})

## ---- multi_live_sample_correlation_data() / multi_live_correlation() / heatmap_data() ----

test_that("multi_live_correlation() computes a real Pearson correlation and p-value matching cor.test()", {
  set.seed(90)
  x <- rnorm(20); y <- x * 2 + rnorm(20, sd = 0.1)
  out <- multi_live_correlation(x, y, method = "pearson")
  ref <- stats::cor.test(x, y, method = "pearson")
  expect_equal(out$r, unname(ref$estimate))
  expect_equal(out$p, ref$p.value)
  expect_true(out$r > 0.9)
})

test_that("multi_live_correlation() refuses with fewer than 3 finite paired observations", {
  out <- multi_live_correlation(c(1, NA, NA), c(1, 2, NA), method = "pearson")
  expect_false(out$ok)
})

test_that("multi_live_correlation_heatmap_data() caps to top-N variable features per side and BH-adjusts p-values", {
  set.seed(100)
  n <- 15
  matA <- matrix(rnorm(n * 10), n, 10, dimnames = list(paste0("S", 1:n), paste0("a", 1:10)))
  matB <- matrix(rnorm(n * 8), n, 8, dimnames = list(paste0("S", 1:n), paste0("b", 1:8)))
  out <- multi_live_correlation_heatmap_data(matA, matB, top_n = 5)
  expect_true(out$ok)
  expect_equal(nrow(out$df), 25L)  ## 5 x 5
  expect_equal(out$df$fdr, stats::p.adjust(out$df$p, method = "BH"))
})

test_that("multi_live_correlation_heatmap_data() refuses with fewer than 3 matched samples", {
  matA <- matrix(rnorm(20), 2, 10, dimnames = list(c("S1", "S2"), paste0("a", 1:10)))
  matB <- matrix(rnorm(20), 2, 10, dimnames = list(c("S1", "S3"), paste0("b", 1:10)))
  out <- multi_live_correlation_heatmap_data(matA, matB)
  expect_false(out$ok)
})

## ---- multi_live_mofa_guardrails() -------------------------------------------

test_that("multi_live_mofa_guardrails() warns on mismatched sample counts and small-sample instability", {
  mat_list <- list(A = matrix(1, 5, 20), B = matrix(1, 5, 20))
  out_small <- multi_live_mofa_guardrails(mat_list)
  expect_false(out_small$ok)
  expect_true(any(grepl("highly unstable", out_small$warnings)))

  mismatched <- list(A = matrix(1, 5, 20), B = matrix(1, 8, 20))
  out_mismatch <- multi_live_mofa_guardrails(mismatched)
  expect_true(any(grepl("do not have the same number", out_mismatch$warnings)))
})

test_that("multi_live_mofa_guardrails() passes clean with adequate sample size and reasonable feature count", {
  mat_list <- list(A = matrix(1, 20, 50), B = matrix(1, 20, 50))
  out <- multi_live_mofa_guardrails(mat_list)
  expect_true(out$ok)
  expect_equal(length(out$warnings), 0L)
})

## ---- multi_dataset_status() / multi_dataset_compatibility() ---------------

test_that("multi_dataset_status() reports 'Not Compatible' for an invalid matrix or too few samples, distinct reasons for each", {
  expect_equal(multi_dataset_status(list(ok = FALSE))$level, "not_compatible")
  few_samples <- multi_dataset_status(list(ok = TRUE, n_samples = 2, pct_missing = 0, n_duplicate_samples = 0, n_duplicate_features = 0, n_zero_variance = 0))
  expect_equal(few_samples$level, "not_compatible")
  expect_true(grepl("Fewer than 3 samples", few_samples$reasons))
})

test_that("multi_dataset_status() reports 'Review Required' with specific reasons for high missingness/duplicates/zero-variance, 'Ready' otherwise", {
  review <- multi_dataset_status(list(ok = TRUE, n_samples = 20, pct_missing = 30, n_duplicate_samples = 1, n_duplicate_features = 0, n_zero_variance = 2))
  expect_equal(review$level, "review")
  expect_equal(length(review$reasons), 3L)  ## missingness + dup samples + zero-variance

  ready <- multi_dataset_status(list(ok = TRUE, n_samples = 20, pct_missing = 0, n_duplicate_samples = 0, n_duplicate_features = 0, n_zero_variance = 0))
  expect_equal(ready$level, "ready")
})

test_that("multi_dataset_status() flags poor cross-layer sample matching (less than 50% of this layer's samples are shared)", {
  out <- multi_dataset_status(list(ok = TRUE, n_samples = 20, pct_missing = 0, n_duplicate_samples = 0, n_duplicate_features = 0, n_zero_variance = 0), n_shared = 5, n_own = 20)
  expect_equal(out$level, "review")
  expect_true(grepl("Only 5 of this dataset's 20 samples", out$reasons))
})

test_that("multi_dataset_compatibility() rolls up per-layer verdicts into the correct overall_label at each severity level", {
  v_ready <- list(ok = TRUE, n_samples = 20, pct_missing = 0, n_duplicate_samples = 0, n_duplicate_features = 0, n_zero_variance = 0)
  v_bad <- list(ok = FALSE)

  none <- multi_dataset_compatibility(list())
  expect_equal(none$overall_label, "NO DATASETS SELECTED")

  not_ready <- multi_dataset_compatibility(list(A = v_ready, B = v_bad))
  expect_equal(not_ready$overall_label, "NOT READY - one or more datasets cannot be used")

  ## Ready per-layer, but no overlap object supplied -> sample_matching_ok is FALSE.
  review_matching <- multi_dataset_compatibility(list(A = v_ready, B = v_ready))
  expect_equal(review_matching$overall_label, "REVIEW REQUIRED - insufficient sample matching")

  good_overlap <- list(ok = TRUE, n_shared = 15, per_layer = list(A = 20, B = 20))
  ready <- multi_dataset_compatibility(list(A = v_ready, B = v_ready), overlap = good_overlap)
  expect_equal(ready$overall_label, "READY")
})
