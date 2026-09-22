## Guard: the single-CpG internal-validation AUC must not re-estimate direction on pooled out-of-fold labels.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "functions", "annotation.R"))
source_from_app_root(file.path("R", "crossomics", "functions", "integration", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "methylomics", "15_Biomarker_Analysis", "mod_methyl_biomarkercard.R"))

test_that("bc_single_cpg_cv does not turn pure noise into a >= 0.5 AUC", {
  set.seed(12)
  y <- rep(c("Case", "Control"), each = 20)
  ## Draw every CpG up front: bc_single_cpg_cv() calls set.seed() internally, so drawing inside the loop
  ## would score the same vector 60 times.
  X <- matrix(runif(40 * 60), nrow = 40)
  aucs <- vapply(seq_len(60), function(i) bc_single_cpg_cv(X[, i], y, "Case", "Control", k = 5)$auc, numeric(1))
  expect_true(all(is.finite(aucs)))
  ## Pooled-label "auto" gave only 6 of 60 below 0.5 and a mean of 0.58; honest scoring gives ~26 and ~0.51.
  expect_gte(sum(aucs < 0.5), 20)
  expect_lt(mean(aucs), 0.53)
})

test_that("bc_single_cpg_cv still reports a clearly separable, hypo-methylated CpG as separable", {
  set.seed(3)
  x <- c(rnorm(20, mean = 0.2, sd = 0.03), rnorm(20, mean = 0.6, sd = 0.03))
  y <- c(rep("Case", 20), rep("Control", 20))
  res <- bc_single_cpg_cv(x, y, "Case", "Control", k = 5)
  expect_true(res$ok)
  expect_gt(res$auc, 0.9)
})
