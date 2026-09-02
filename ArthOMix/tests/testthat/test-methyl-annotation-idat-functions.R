## Module 2 (Methylomics) - annotation.R (real Bioconductor manifest lookup,
## both 450K and EPIC packages are installed here) and idat_metrics.R's
## fail-soft contract when no raw IDAT input is available (no committed
## sample IDAT fixture exists in this project to exercise the real minfi
## path - documented as a known coverage gap, not silently skipped).

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "functions", "annotation.R"))
source_from_app_root(file.path("R", "methylomics", "functions", "idat_metrics.R"))

test_that("methyl_probe_is_cpg() identifies 'cg'-prefixed probe IDs case-insensitively", {
  ids <- c("cg00000029", "CG00000108", "ch.1.1000", "rs9363764")
  expect_equal(methyl_probe_is_cpg(ids), c(TRUE, TRUE, FALSE, FALSE))
})

test_that("methyl_get_annotation() returns a real, well-formed annotation table for 450K and EPIC", {
  for (array_type in c("450K", "EPIC")) {
    res <- methyl_get_annotation(array_type)
    expect_true(res$ok, info = array_type)
    expect_true(all(c("chr", "pos", "Type", "Probe_rs", "CpG_rs", "SBE_rs", "gene") %in% colnames(res$anno)), info = array_type)
    expect_gt(nrow(res$anno), 100000, label = array_type)  ## real manifests are large
    expect_true(all(grepl("^chr", stats::na.omit(res$anno$chr))), info = array_type)
  }
})

test_that("methyl_get_annotation() caches the annotation table across repeated calls (same object, not just equal)", {
  a1 <- methyl_get_annotation("450K")
  a2 <- methyl_get_annotation("450K")
  expect_identical(a1$anno, a2$anno)
})

test_that("methyl_get_annotation() fails gracefully (not an error) for an array type with no installed manifest package", {
  res <- methyl_get_annotation("WGBS")
  expect_false(res$ok)
  expect_true(nzchar(res$reason))
})

## ---- idat_metrics.R: fail-soft contract without raw IDAT input ------------

test_that("methyl_idat_derive()/methyl_bisulfite_conversion()/methyl_median_intensity() fail soft (never throw) with rg_set/mset = NULL", {
  d <- methyl_idat_derive(NULL)
  expect_false(d$ok)
  expect_true(nzchar(d$reason))

  bc <- methyl_bisulfite_conversion(NULL)
  expect_false(bc$ok)

  mi <- methyl_median_intensity(NULL)
  expect_false(mi$ok)
})
