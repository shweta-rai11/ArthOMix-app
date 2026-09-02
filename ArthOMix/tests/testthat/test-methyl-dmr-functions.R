## Module 2 (Methylomics) - DMR tab's pure functions: the region-level
## filter (extends DMP's filter with region-width/CpG-count thresholds),
## sex-label lookup, and the region x sample heatmap matrix builder (real
## GenomicRanges::findOverlaps() computation). The live DMRcate-calling
## server reactive itself is not exercised here - it needs a realistic
## genome-wide correlated-CpG-cluster structure to find any region at all,
## which a small synthetic fixture cannot meaningfully provide; documented
## as a coverage gap rather than tested against a fixture that would prove
## nothing.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "05_Differential_Methylation_Position", "mod_methyl_dmp.R"))
source_from_app_root(file.path("R", "methylomics", "06_Differential_Methylation_Region", "mod_methyl_dmr.R"))

## ---- mod_methyl_dmr_filter() ------------------------------------------------

test_that("mod_methyl_dmr_filter() applies the DMP-style FDR/effect filter plus region CpG-count and width thresholds", {
  df <- data.frame(
    region = c("r1", "r2", "r3", "r4"),
    fdr = c(0.001, 0.001, 0.001, 0.001),
    dbeta = c(0.2, 0.2, 0.2, 0.2),
    no.cpgs = c(5, 1, 5, 5),
    width = c(500, 500, 50, 5000),
    stringsAsFactors = FALSE
  )
  out <- mod_methyl_dmr_filter(df, "fdr", "dbeta", fdr_max = 0.05, effect_min = 0.05,
                                 direction = "any", min_cpgs = 2, min_width = 100, max_width = 2000)
  expect_equal(out$region, "r1")  ## r2 fails cpg count, r3 fails min width, r4 fails max width
})

## ---- mod_methyl_dmr_sex_label() ---------------------------------------------

test_that("mod_methyl_dmr_sex_label() matches DMP's own sex-choice labels for consistent banners", {
  sheet <- data.frame(sex = c("F", "M"))
  expect_equal(mod_methyl_dmr_sex_label(sheet, "sex", "F"), "Female")
  expect_equal(mod_methyl_dmr_sex_label(sheet, "sex", "M"), "Male")
  expect_equal(mod_methyl_dmr_sex_label(sheet, "sex", "__all__"), "All samples")
})

test_that("mod_methyl_dmr_sex_label() falls back to the raw value when no matching label is found", {
  expect_equal(mod_methyl_dmr_sex_label(data.frame(sex = character(0)), "sex", "weird_value"), "weird_value")
})

## ---- mod_methyl_dmr_topplot() -----------------------------------------------

test_that("mod_methyl_dmr_topplot() ranks by region-level FDR and labels regions by coordinates (+gene when annotated)", {
  dt <- data.frame(
    seqnames = c("chr1", "chr2"), start = c(1000, 2000), end = c(1100, 2100),
    dmr_fdr = c(0.001, 0.5), meandiff = c(0.1, -0.1),
    overlapping.genes = c("TP53", NA_character_), stringsAsFactors = FALSE
  )
  out <- mod_methyl_dmr_topplot(dt, n = 5)
  expect_true(any(grepl("TP53", levels(out$data$label))))
})

test_that("mod_methyl_dmr_topplot() errors clearly when no region has both fdr and meandiff", {
  dt <- data.frame(seqnames = "chr1", start = 1, end = 2, dmr_fdr = NA_real_, meandiff = NA_real_,
                     overlapping.genes = NA_character_, stringsAsFactors = FALSE)
  expect_error(mod_methyl_dmr_topplot(dt), class = "validation")
})

## ---- mod_methyl_dmr_heatmap() (real GenomicRanges overlap computation) -----

test_that("mod_methyl_dmr_heatmap() correctly averages constituent-probe beta values per region via real genomic overlap", {
  ## 2 regions on chr1: region A (1000-1100) contains probes at 1010/1050;
  ## region B (5000-5100) contains probe at 5050. A 4th probe (9999) is
  ## outside both regions and must be excluded from both region means.
  sig_dt <- data.frame(seqnames = c("chr1", "chr1"), start = c(1000, 5000), end = c(1100, 5100),
                         dmr_fdr = c(0.001, 0.002), dmr_id = c("DMR1", "DMR2"), stringsAsFactors = FALSE)
  probe_chr <- c("chr1", "chr1", "chr1", "chr1")
  probe_pos <- c(1010, 1050, 5050, 9999)
  beta_scale <- matrix(c(0.2, 0.4, 0.9, 0.1,   ## S1: probe1=0.2, probe2=0.4, probe3=0.9, probe4=0.1
                          0.6, 0.8, 0.3, 0.5),  ## S2
                         nrow = 4, ncol = 2, dimnames = list(NULL, c("S1", "S2")))
  grp <- factor(c("HC", "RA"), levels = c("HC", "RA"))

  gg <- mod_methyl_dmr_heatmap(sig_dt, beta_scale, probe_chr, probe_pos, grp)
  plot_data <- gg$data
  ## Region DMR1, sample S1 = mean(probe1, probe2) for S1 = mean(0.2, 0.4) = 0.3
  expect_equal(plot_data$beta[plot_data$region == "DMR1" & plot_data$sample == "S1"], 0.3)
  ## Region DMR2, sample S1 = probe3 alone = 0.9
  expect_equal(plot_data$beta[plot_data$region == "DMR2" & plot_data$sample == "S1"], 0.9)
})

test_that("mod_methyl_dmr_heatmap() errors clearly when there are no significant regions to show", {
  empty_dt <- data.frame(seqnames = character(0), start = numeric(0), end = numeric(0), dmr_fdr = numeric(0))
  expect_error(mod_methyl_dmr_heatmap(empty_dt, matrix(numeric(0), 0, 0), character(0), numeric(0), factor()), class = "validation")
})

## ---- methyl_dmr_engine_pkgs_ok() ---------------------------------------------

test_that("methyl_dmr_engine_pkgs_ok() reports the real installed-package status truthfully", {
  ok <- methyl_dmr_engine_pkgs_ok()
  expect_type(ok, "logical")
  ## DMRcate/bumphunter are in this project's dependency surface - expect them installed here.
  expect_true(ok)
})
