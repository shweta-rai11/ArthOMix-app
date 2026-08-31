## Module 1 (Transcriptomics) - WGCNA's top-level pure helper functions.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "mod_wgcna.R"))

test_that("wgcna_encode_trait() passes numeric traits through unchanged", {
  meta <- data.frame(age = c(30, 45, 50, 22))
  expect_equal(wgcna_encode_trait(meta, "age"), c(30, 45, 50, 22))
})

test_that("wgcna_encode_trait() factor-encodes a categorical trait alphabetically", {
  meta <- data.frame(group = c("RA", "HC", "HC", "RA"))
  out <- wgcna_encode_trait(meta, "group")
  expect_equal(out, c(2, 1, 1, 2))  ## HC=1, RA=2 alphabetically
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
