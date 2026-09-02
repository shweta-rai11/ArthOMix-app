## Module 2 (Methylomics) - Candidate CpGs' pure functions: chromosome-name
## normalization, DMR filtering, and the WGCNA-module x DMR genomic overlap
## join (real GenomicRanges::findOverlaps() computation, the scientific
## heart of this submodule).

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "09_Candidate_CpGs", "mod_methyl_candidates.R"))

## ---- mcd_norm_chr() ----------------------------------------------------------

test_that("mcd_norm_chr() normalizes bare chromosome numbers/letters and leaves already-prefixed names alone", {
  expect_equal(mcd_norm_chr(c("1", "X", "chr2", "CHR3", "y")), c("chr1", "chrX", "chr2", "chr3", "chrY"))
})

## ---- mcd_filter_dmrs() -------------------------------------------------------

test_that("mcd_filter_dmrs() applies FDR, delta-beta, CpG-count, direction, and chromosome filters together", {
  dmr <- data.frame(
    dmr_id = paste0("d", 1:5),
    dmr_fdr = c(0.001, 0.001, 0.001, 0.001, 0.2),
    delta_beta = c(0.2, 0.01, 0.2, -0.2, 0.2),
    n_cpgs = c(5, 5, 1, 5, 5),
    direction = c("hyper", "hyper", "hyper", "hypo", "hyper"),
    chr = c("chr1", "chr1", "chr1", "chr2", "chr1"),
    stringsAsFactors = FALSE
  )
  out <- mcd_filter_dmrs(dmr, fdr_max = 0.05, dbeta_min = 0.1, mincpgs_min = 3, direction = "hyper", chr_restrict = "chr1")
  expect_equal(out$dmr_id, "d1")  ## d2 fails dbeta, d3 fails cpg count, d4 fails direction+chr, d5 fails fdr
})

test_that("mcd_filter_dmrs() with no thresholds set is a no-op", {
  dmr <- data.frame(dmr_id = c("d1", "d2"), dmr_fdr = c(0.001, 0.9))
  out <- mcd_filter_dmrs(dmr)
  expect_equal(nrow(out), 2L)
})

## ---- mcd_compute_overlap() (real GenomicRanges join) --------------------------

test_that("mcd_compute_overlap() correctly joins WGCNA module CpGs to overlapping DMRs by genomic position", {
  module_assign <- data.frame(cpg = c("cg1", "cg2", "cg3", "cg4"), module = c("turquoise", "turquoise", "blue", "blue"),
                                stringsAsFactors = FALSE)
  annotation <- data.frame(cpg = c("cg1", "cg2", "cg3", "cg4"), chr = c("chr1", "chr1", "chr1", "chr2"),
                             pos = c(1010, 1050, 9999, 1010), stringsAsFactors = FALSE)
  ## Only one DMR, on chr1:1000-1100 - overlaps cg1/cg2 but not cg3 (outside)
  ## or cg4 (wrong chromosome).
  dmr_filtered <- data.frame(dmr_id = "d1", chr = "chr1", start = 1000, end = 1100,
                               delta_beta = 0.3, n_cpgs = 5, gene = "TP53", direction = "hyper",
                               stringsAsFactors = FALSE)

  out <- mcd_compute_overlap(module_assign, annotation, dmr_filtered)
  expect_setequal(out$joined$cpg, c("cg1", "cg2"))
  expect_true(all(out$joined$dmr_gene == "TP53" | out$joined$gene == "TP53"))
  expect_equal(nrow(out$cpg_universe), 4L)  ## all 4 CpGs have resolvable coordinates
})

test_that("mcd_compute_overlap() respects the flank argument, rescuing a CpG just outside the raw DMR boundary", {
  module_assign <- data.frame(cpg = "cg1", module = "turquoise", stringsAsFactors = FALSE)
  annotation <- data.frame(cpg = "cg1", chr = "chr1", pos = 1150, stringsAsFactors = FALSE)  ## 50bp past the DMR end
  dmr_filtered <- data.frame(dmr_id = "d1", chr = "chr1", start = 1000, end = 1100,
                               delta_beta = 0.2, n_cpgs = 3, gene = NA_character_, direction = "hyper",
                               stringsAsFactors = FALSE)

  no_flank <- mcd_compute_overlap(module_assign, annotation, dmr_filtered, flank = 0)
  expect_equal(nrow(no_flank$joined), 0L)

  with_flank <- mcd_compute_overlap(module_assign, annotation, dmr_filtered, flank = 100)
  expect_equal(nrow(with_flank$joined), 1L)
})

test_that("mcd_compute_overlap() returns an empty joined table (not an error) when no DMRs or no resolvable CpGs are given", {
  module_assign <- data.frame(cpg = "cg1", module = "turquoise", stringsAsFactors = FALSE)
  annotation <- data.frame(cpg = "cg1", chr = "chr1", pos = 1010, stringsAsFactors = FALSE)
  empty_dmr <- data.frame(dmr_id = character(0), chr = character(0), start = numeric(0), end = numeric(0),
                            delta_beta = numeric(0), n_cpgs = numeric(0), gene = character(0), direction = character(0))
  out <- mcd_compute_overlap(module_assign, annotation, empty_dmr)
  expect_equal(nrow(out$joined), 0L)
})
