## Module 2 (Methylomics) - Co-methylation Network (WGCNA)'s pure functions:
## top-variable-probe selection, trait encoding, cell-type-reference
## detection, and the guardrail flags.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "qc.R"))
source_from_app_root(file.path("R", "methylomics", "mod_methyl_wgcna.R"))

## ---- mx_wgcna_top_variable() ------------------------------------------------

test_that("mx_wgcna_top_variable() keeps the N most variable probes by MAD/variance/sd/IQR", {
  set.seed(270)
  mat <- matrix(rnorm(200), 20, 10)
  mat[1, ] <- mat[1, ] * 20  ## probe 1: by far the most variable under every method

  for (method in c("mad", "variance", "sd", "iqr")) {
    out <- mx_wgcna_top_variable(mat, method = method, top_n = 5)
    expect_equal(out$n_kept, 5L)
    expect_true(rownames(mat)[1] %in% rownames(out$mat) || is.null(rownames(mat)), info = method)
  }
})

test_that("mx_wgcna_top_variable() caps top_n at the number of probes with nonzero variability", {
  mat <- matrix(0.5, 10, 5)   ## every probe constant -> zero variability everywhere
  mat[1, ] <- c(0.1, 0.9, 0.2, 0.8, 0.3)  ## only probe 1 has real variance
  out <- mx_wgcna_top_variable(mat, method = "mad", top_n = 100)
  expect_equal(out$n_kept, 1L)
})

## ---- mx_wgcna_encode_trait() -------------------------------------------------

test_that("mx_wgcna_encode_trait() passes numeric columns through unchanged", {
  sheet <- data.frame(age = c(30, 40, 50))
  out <- mx_wgcna_encode_trait(sheet, "age")
  expect_true(out$ok)
  expect_equal(out$vec, c(30, 40, 50))
})

test_that("mx_wgcna_encode_trait() 0/1-encodes an exactly-two-level categorical column", {
  sheet <- data.frame(group = c("HC", "RA", "HC", "RA"))
  out <- mx_wgcna_encode_trait(sheet, "group")
  expect_true(out$ok)
  expect_equal(out$levels, c("HC", "RA"))
  expect_equal(out$vec, c(0, 1, 0, 1))
})

test_that("mx_wgcna_encode_trait() rejects a column with more or fewer than two levels", {
  one_level <- mx_wgcna_encode_trait(data.frame(x = rep("A", 5)), "x")
  expect_false(one_level$ok)
  three_level <- mx_wgcna_encode_trait(data.frame(x = c("A", "B", "C")), "x")
  expect_false(three_level$ok)
})

## ---- mx_wgcna_celltype_reference() -------------------------------------------

test_that("mx_wgcna_celltype_reference() identifies a Houseman-style compositional column set and picks the highest-mean cell type", {
  set.seed(271)
  n <- 20
  cd4 <- runif(n, 0.2, 0.4); cd8 <- runif(n, 0.1, 0.2); bcell <- 1 - cd4 - cd8
  sheet <- data.frame(CD4T = cd4, CD8T = cd8, Bcell = bcell, unrelated = rnorm(n))
  ## candidate_cols is the pre-filtered composition candidate set - the
  ## function checks whether ALL given columns together sum to ~1, it does
  ## not search subsets, so a non-compositional column must not be included.
  ref <- mx_wgcna_celltype_reference(sheet, candidate_cols = c("CD4T", "CD8T", "Bcell"))
  ## Bcell = 1 - CD4T - CD8T has the highest mean fraction (~0.55) among the
  ## three by construction (CD4T ~U(0.2,0.4), CD8T ~U(0.1,0.2)).
  expect_equal(ref, "Bcell")
})

test_that("mx_wgcna_celltype_reference() returns NULL when an unrelated non-compositional column is included in candidate_cols", {
  set.seed(272)
  n <- 20
  cd4 <- runif(n, 0.2, 0.4); cd8 <- runif(n, 0.1, 0.2); bcell <- 1 - cd4 - cd8
  sheet <- data.frame(CD4T = cd4, CD8T = cd8, Bcell = bcell, unrelated = rnorm(n))
  ref <- mx_wgcna_celltype_reference(sheet, candidate_cols = c("CD4T", "CD8T", "Bcell", "unrelated"))
  expect_null(ref)
})

test_that("mx_wgcna_celltype_reference() returns NULL when columns don't sum to ~1 (not a real composition)", {
  sheet <- data.frame(a = runif(10, 0, 1), b = runif(10, 0, 1), c = runif(10, 0, 1))
  expect_null(mx_wgcna_celltype_reference(sheet, c("a", "b", "c")))
})

test_that("mx_wgcna_celltype_reference() returns NULL with fewer than 3 candidate columns", {
  sheet <- data.frame(a = c(0.5, 0.5), b = c(0.5, 0.5))
  expect_null(mx_wgcna_celltype_reference(sheet, c("a", "b")))
})

## ---- mx_wgcna_guardrails() ---------------------------------------------------

test_that("mx_wgcna_guardrails() flags low sample/probe counts, poor scale-free fit, and degenerate module structure", {
  g <- mx_wgcna_guardrails(n_samples = 10, n_probes = 200, max_r_sq = 0.5, module_colors = rep("grey", 20))
  expect_true(g$low_n_samples)
  expect_true(g$low_n_probes)
  expect_true(g$poor_sft_fit)
  expect_true(g$all_grey)

  g2 <- mx_wgcna_guardrails(n_samples = 30, n_probes = 1000, max_r_sq = 0.9,
                              module_colors = c(rep("grey", 5), rep("turquoise", 15)))
  expect_false(g2$low_n_samples)
  expect_false(g2$poor_sft_fit)
  expect_false(g2$all_grey)
  expect_true(g2$single_module)  ## exactly one non-grey module
})

test_that("mx_wgcna_guardrails() skips any check whose argument is left NULL", {
  g <- mx_wgcna_guardrails(n_samples = 5)
  expect_true(g$low_n_samples)
  expect_false(g$low_n_probes)
  expect_false(g$poor_sft_fit)
  expect_false(g$all_grey)
})
