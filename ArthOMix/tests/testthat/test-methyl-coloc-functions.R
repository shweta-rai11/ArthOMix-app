## Module 2 (Methylomics) - Colocalization's pure functions: the LD-matrix
## alignment helper and the posterior-probability interpretation/verdict
## logic (the single source of truth every hypothesis-classification view

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "12_Colocalization", "mod_methyl_coloc.R"))

test_that(".mcol_prep_ld() extracts the requested SNP submatrix from a raw LD table", {
  raw_df <- data.frame(snp = c("rs1", "rs2", "rs3"), rs1 = c(1, 0.8, 0.1), rs2 = c(0.8, 1, 0.2), rs3 = c(0.1, 0.2, 1))
  out <- .mcol_prep_ld(raw_df, snp_ids = c("rs1", "rs2"))
  expect_equal(dim(out), c(2L, 2L))
  expect_equal(unname(out["rs1", "rs2"]), 0.8)
})

test_that(".mcol_prep_ld() returns NULL for malformed input or too few matched SNPs", {
  expect_null(.mcol_prep_ld(NULL, c("rs1")))
  expect_null(.mcol_prep_ld(data.frame(a = 1), c("rs1")))
  raw_df <- data.frame(snp = c("rs1", "rs2"), rs1 = c(1, 0.5), rs2 = c(0.5, 1))
  expect_null(.mcol_prep_ld(raw_df, snp_ids = "rs9"))
})

test_that(".mcol_verdict() classifies coloc-supported / coloc-refuted / inconclusive by the PP.H4/H3 threshold", {
  expect_equal(.mcol_verdict(h3 = 0.05, h4 = 0.85), "coloc-supported (shared causal variant)")
  expect_equal(.mcol_verdict(h3 = 0.85, h4 = 0.05), "coloc-refuted (distinct causal variants)")
  expect_true(grepl("inconclusive", .mcol_verdict(h3 = 0.3, h4 = 0.3)))
})

test_that(".mcol_verdict() respects a custom threshold", {
  expect_equal(.mcol_verdict(h3 = 0.1, h4 = 0.55, threshold = 0.5), "coloc-supported (shared causal variant)")
  expect_true(grepl("inconclusive", .mcol_verdict(h3 = 0.1, h4 = 0.55, threshold = 0.8)))
})

test_that(".mcol_interpret() identifies the strongest hypothesis and reports its posterior probability", {
  out <- .mcol_interpret(h0 = 0.01, h1 = 0.02, h2 = 0.02, h3 = 0.05, h4 = 0.9)
  txt <- paste(sapply(out, function(x) if (inherits(x, "shiny.tag")) htmltools::doRenderTags(x) else ""), collapse = " ")
  expect_true(grepl("Strongest-supported hypothesis: H4", txt))
  expect_true(grepl("90.0%", txt))
  expect_true(grepl("Strong evidence", txt))
})

test_that(".mcol_interpret() gives the 'moderate evidence' framing for H4 just below the project's 0.8 threshold", {
  out <- .mcol_interpret(h0 = 0.05, h1 = 0.05, h2 = 0.05, h3 = 0.1, h4 = 0.75)
  txt <- paste(sapply(out, function(x) if (inherits(x, "shiny.tag")) htmltools::doRenderTags(x) else ""), collapse = " ")
  expect_true(grepl("Moderate evidence", txt))
})
