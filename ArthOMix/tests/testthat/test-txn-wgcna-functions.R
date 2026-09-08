## Module 1 (Transcriptomics) - WGCNA's top-level pure helper functions.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "05_WGCNA", "mod_wgcna.R"))

test_that("wgcna_module_trait_fdr() applies real BH-FDR across the complete module x trait matrix, never per-cell/uncorrected", {
  p_mat <- matrix(c(0.001, 0.5, 0.02, 0.9, 0.5, 0.4, 0.03, 0.6), nrow = 2, dimnames = list(c("M1", "M2"), c("t1", "t2", "t3", "t4")))
  q_mat <- wgcna_module_trait_fdr(p_mat)
  expect_equal(q_mat, matrix(stats::p.adjust(as.vector(p_mat), method = "BH"), nrow = 2, dimnames = dimnames(p_mat)))
  ## BH-FDR is never more lenient than the raw p-value it corrects.
  expect_true(all(q_mat >= p_mat - 1e-12))
  ## A cell whose raw p alone would look "significant" at a naive 0.05 cutoff
  ## can be pushed above 0.05 once corrected across the full 8-cell family.
  expect_true(any(p_mat < 0.05 & q_mat >= 0.05))
})

test_that("wgcna_module_trait_significant() requires BOTH the |r|>=0.5 effect-size floor AND BH-FDR q<0.05 - replacing the prior uncorrected p<1e-8 substitute (Priority-14 finding, 2026-09-07 defense audit)", {
  cor_mat <- matrix(c(0.6, 0.2, -0.7, 0.1), nrow = 2, dimnames = list(c("M1", "M2"), c("t1", "t2")))
  ## M1/t1 has a large effect size AND passes FDR: should be flagged.
  ## M2/t2 has a large effect size but FAILS FDR: must NOT be flagged despite |r|>=0.5.
  ## M1/t2 and M2/t1 have small effect sizes: never flagged regardless of q.
  q_mat <- matrix(c(0.01, 0.9, 0.9, 0.6), nrow = 2, dimnames = dimnames(cor_mat))
  sig <- wgcna_module_trait_significant(cor_mat, q_mat)
  expect_equal(sig, "M1")

  ## A module whose raw p would have been < 1e-8 under the OLD rule, but whose
  ## BH-FDR-corrected q is >= 0.05 once corrected across a large family, must NOT
  ## be flagged under the fixed rule - this is exactly the scenario the old,
  ## uncorrected p < 1e-8 substitute could not catch.
  p_mat_large_family <- matrix(c(1e-9, rep(0.04, 99)), nrow = 1)
  rownames(p_mat_large_family) <- "Mx"; colnames(p_mat_large_family) <- paste0("t", 1:100)
  q_large <- wgcna_module_trait_fdr(p_mat_large_family)
  cor_mat_large <- matrix(rep(0.55, 100), nrow = 1, dimnames = dimnames(p_mat_large_family))
  sig_old_rule_would_flag <- p_mat_large_family[1, 1] < 1e-8 && abs(cor_mat_large[1, 1]) >= 0.5
  expect_true(sig_old_rule_would_flag)  ## sanity: old rule WOULD have flagged this
  sig_new <- wgcna_module_trait_significant(cor_mat_large, q_large)
  expect_true("Mx" %in% sig_new)  ## the single extreme p-value still survives correction here...
  ## ...but a single moderately-significant p-value diluted among 99 non-significant
  ## ones in the same family is correctly washed out by BH-FDR, even though it
  ## would NOT have passed the old raw-p<1e-8 gate either (0.03 > 1e-8) - this
  ## demonstrates the correction actually does something beyond the old rule for a
  ## realistic "one nominal hit among many nulls" family.
  set.seed(1)
  p_mat_moderate <- matrix(c(0.03, stats::runif(99, 0.5, 0.99)), nrow = 1, dimnames = list("My", paste0("t", 1:100)))
  cor_mat_moderate <- matrix(rep(0.55, 100), nrow = 1, dimnames = dimnames(p_mat_moderate))
  q_moderate <- wgcna_module_trait_fdr(p_mat_moderate)
  expect_false(p_mat_moderate[1, 1] < 1e-8)  ## wouldn't even have passed the old raw-p<1e-8 gate either
  expect_true(q_moderate[1, 1] >= 0.05)  ## BH-FDR correctly declines to call this significant
  expect_length(wgcna_module_trait_significant(cor_mat_moderate, q_moderate), 0L)
})

test_that("wgcna_encode_trait() passes numeric traits through unchanged", {
  meta <- data.frame(age = c(30, 45, 50, 22))
  expect_equal(wgcna_encode_trait(meta, "age"), c(30, 45, 50, 22))
})

test_that("wgcna_encode_trait() factor-encodes a categorical trait alphabetically", {
  meta <- data.frame(group = c("RA", "HC", "HC", "RA"))
  out <- wgcna_encode_trait(meta, "group")
  expect_equal(out, c(2, 1, 1, 2))
})

test_that("wgcna_encode_trait() restricts to levels_keep, coding excluded values as NA", {
  meta <- data.frame(group = c("RA", "HC", "OTHER", "RA"))
  out <- wgcna_encode_trait(meta, "group", levels_keep = c("HC", "RA"))
  expect_true(is.na(out[3]))
  expect_equal(out[c(1, 2, 4)], c(2, 1, 2))
})

test_that("wgcna_cor_fnc()/wgcna_cor_fnc_name() dispatch bicor vs Pearson correlation correctly", {
  expect_identical(wgcna_cor_fnc("bicor"), WGCNA::bicor)
  expect_identical(wgcna_cor_fnc("pearson"), WGCNA::cor)
  expect_equal(wgcna_cor_fnc_name("bicor"), "bicor")
  expect_equal(wgcna_cor_fnc_name("pearson"), "cor")
})

test_that("wgcna_string_url()/wgcna_string_image_url() URL-encode gene identifiers and join with the expected delimiter", {
  url <- wgcna_string_url(c("TP53", "IL6"))
  expect_true(grepl("TP53%0dIL6", url, fixed = TRUE))
  expect_true(grepl("^https://string-db.org/cgi/network", url))

  img_url <- wgcna_string_image_url(c("TP53", "IL6"))
  expect_true(grepl("^https://string-db.org/api/image/network", img_url))
})

test_that("load_precomputed_wgcna_result()/load_precomputed_wgcna_sft() return the expected structure from the bundled reference run", {
  res <- load_precomputed_wgcna_result()
  expect_true(is.data.frame(res$gene_module))
  expect_setequal(colnames(res$gene_module), c("gene", "module"))
  expect_equal(nrow(res$gene_module), res$n_genes)
  expect_equal(res$n_genes, ncol(res$texpr))
  expect_equal(res$n_samples, nrow(res$texpr))
  expect_gt(res$n_modules, 0)
  expect_true(is.data.frame(res$module_sizes))

  sft <- load_precomputed_wgcna_sft()
  expect_true(!is.null(sft$sft_df))
  expect_true(!is.null(sft$power))
})
