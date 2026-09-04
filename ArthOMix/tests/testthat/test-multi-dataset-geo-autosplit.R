suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "multiomics", "06_Gene_CpG_Mapping", "multiomics_mapping_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "multiomics_dataset_helpers.R"))

fx_rnaseq_eset <- function(n_genes = 20, n_samples = 8, seed = 1) {
  set.seed(seed)
  genes <- sprintf("GENE%d", seq_len(n_genes))
  samples <- sprintf("GSM%04d", 2000 + seq_len(n_samples))
  mat <- matrix(runif(n_genes * n_samples, 3, 12), n_genes, n_samples, dimnames = list(genes, samples))
  pdat <- data.frame(title = sprintf("patient_%d", seq_len(n_samples)), geo_accession = samples,
                      row.names = samples, stringsAsFactors = FALSE)
  eset <- Biobase::ExpressionSet(assayData = mat, phenoData = Biobase::AnnotatedDataFrame(pdat))
  Biobase::annotation(eset) <- "GPL_RNA"
  eset
}

fx_methylation_eset <- function(n_probes = 20, n_samples = 8, seed = 2) {
  set.seed(seed)
  probes <- sprintf("cg%08d", seq_len(n_probes))
  samples <- sprintf("GSM%04d", 3000 + seq_len(n_samples))
  mat <- matrix(runif(n_probes * n_samples, 0, 1), n_probes, n_samples, dimnames = list(probes, samples))
  pdat <- data.frame(title = sprintf("patient_%d", seq_len(n_samples)), geo_accession = samples,
                      row.names = samples, stringsAsFactors = FALSE)
  eset <- Biobase::ExpressionSet(assayData = mat, phenoData = Biobase::AnnotatedDataFrame(pdat))
  Biobase::annotation(eset) <- "GPL_METH"
  eset
}

fx_empty_eset <- function(n_samples = 8) {
  samples <- sprintf("GSM%04d", 4000 + seq_len(n_samples))
  mat <- matrix(numeric(0), 0, n_samples, dimnames = list(NULL, samples))
  pdat <- data.frame(title = sprintf("patient_%d", seq_len(n_samples)), row.names = samples, stringsAsFactors = FALSE)
  eset <- Biobase::ExpressionSet(assayData = mat, phenoData = Biobase::AnnotatedDataFrame(pdat))
  Biobase::annotation(eset) <- "GPL_EMPTY"
  eset
}

geo_quick_lines <- function(acc, relations = character(0)) {
  c(sprintf("^SERIES = %s", acc), sprintf("!Series_geo_accession = %s", acc),
    sprintf("!Series_relation = %s", relations))
}

test_that("multi_geo_series_relation() rejects a malformed accession without any network call", {
  called <- FALSE
  testthat::local_mocked_bindings(readLines = function(...) { called <<- TRUE; character(0) }, .package = "base")
  res <- multi_geo_series_relation("not-an-accession")
  expect_false(res$ok)
  expect_false(called)
})

test_that("multi_geo_series_relation() parses 'SuperSeries of:' relation lines into sub-series accessions", {
  testthat::local_mocked_bindings(
    readLines = function(...) geo_quick_lines("GSE201754", c("SuperSeries of: GSE201752", "SuperSeries of: GSE201753", "BioProject: https://example.com/1")),
    .package = "base"
  )
  res <- multi_geo_series_relation("gse201754")
  expect_true(res$ok)
  expect_setequal(res$subseries, c("GSE201752", "GSE201753"))
})

test_that("multi_geo_series_relation() reports a clear error for an accession GEO doesn't recognize as a Series", {
  testthat::local_mocked_bindings(readLines = function(...) c("<!DOCTYPE HTML PUBLIC>", "<HTML>"), .package = "base")
  res <- multi_geo_series_relation("GSE99999999")
  expect_false(res$ok)
  expect_match(res$error, "Series")
})

test_that("multi_geo_series_relation() reports a clear error when the network call itself fails", {
  testthat::local_mocked_bindings(readLines = function(...) stop("connection timed out"), .package = "base")
  res <- multi_geo_series_relation("GSE12345")
  expect_false(res$ok)
  expect_match(res$error, "NCBI GEO")
})

test_that("multi_geo_autosplit_fetch() splits a single accession that itself carries both platforms directly, without needing relation lookup", {
  relation_called <- FALSE
  testthat::local_mocked_bindings(readLines = function(...) { relation_called <<- TRUE; character(0) }, .package = "base")
  testthat::local_mocked_bindings(
    getGEO = function(...) list(GPL_RNA = fx_rnaseq_eset(), GPL_METH = fx_methylation_eset()),
    .package = "GEOquery"
  )
  res <- multi_geo_autosplit_fetch("GSE11111")
  expect_true(res$ok)
  expect_equal(res$expression$accession, "GSE11111")
  expect_equal(res$methylation$accession, "GSE11111")
  expect_equal(dim(res$expression$mat), c(8L, 20L))
  expect_equal(dim(res$methylation$mat), c(8L, 20L))
  expect_false(relation_called)
})

test_that("multi_geo_autosplit_fetch() falls back to SuperSeries relation discovery when the accession alone is single-platform", {
  testthat::local_mocked_bindings(
    readLines = function(...) geo_quick_lines("GSE201754", c("SuperSeries of: GSE201752", "SuperSeries of: GSE201753")),
    .package = "base"
  )
  testthat::local_mocked_bindings(
    getGEO = function(acc, ...) {
      acc <- toupper(acc)
      if (identical(acc, "GSE201754")) list(GPL_ONE = fx_rnaseq_eset())     # the SuperSeries itself: single platform, not a pair
      else if (identical(acc, "GSE201752")) list(GPL_METH = fx_methylation_eset())
      else if (identical(acc, "GSE201753")) list(GPL_RNA = fx_rnaseq_eset())
      else stop("unexpected accession in test: ", acc)
    },
    .package = "GEOquery"
  )
  res <- multi_geo_autosplit_fetch("GSE201754")
  expect_true(res$ok)
  expect_equal(res$expression$accession, "GSE201753")
  expect_equal(res$methylation$accession, "GSE201752")
})

test_that("multi_geo_autosplit_fetch() errors out (does not fabricate a split) when GEO lists fewer than two linked sub-series", {
  testthat::local_mocked_bindings(readLines = function(...) geo_quick_lines("GSE1", character(0)), .package = "base")
  testthat::local_mocked_bindings(getGEO = function(...) list(GPL_ONE = fx_rnaseq_eset()), .package = "GEOquery")
  res <- multi_geo_autosplit_fetch("GSE1")
  expect_false(res$ok)
  expect_match(res$error, "SuperSeries")
})

test_that("multi_geo_autosplit_fetch() errors out when a linked sub-series' own series matrix is empty (data only in a supplementary file)", {
  testthat::local_mocked_bindings(
    readLines = function(...) geo_quick_lines("GSE201754", c("SuperSeries of: GSE201752", "SuperSeries of: GSE201753")),
    .package = "base"
  )
  testthat::local_mocked_bindings(
    getGEO = function(acc, ...) {
      acc <- toupper(acc)
      if (identical(acc, "GSE201754")) list(GPL_ONE = fx_rnaseq_eset())
      else if (identical(acc, "GSE201752")) list(GPL_EMPTY = fx_empty_eset())
      else if (identical(acc, "GSE201753")) list(GPL_RNA = fx_rnaseq_eset())
      else stop("unexpected accession in test: ", acc)
    },
    .package = "GEOquery"
  )
  res <- multi_geo_autosplit_fetch("GSE201754")
  expect_false(res$ok)
  expect_match(res$error, "Upload Dataset")
})

test_that("multi_geo_autosplit_fetch() errors out when both linked sub-series classify as the same omics type (ambiguous, not guessed)", {
  testthat::local_mocked_bindings(
    readLines = function(...) geo_quick_lines("GSE2", c("SuperSeries of: GSE3", "SuperSeries of: GSE4")),
    .package = "base"
  )
  testthat::local_mocked_bindings(
    getGEO = function(acc, ...) {
      acc <- toupper(acc)
      if (identical(acc, "GSE2")) list(GPL_ONE = fx_rnaseq_eset())
      else list(GPL_RNA = fx_rnaseq_eset())   # both GSE3 and GSE4 look like RNA-seq
    },
    .package = "GEOquery"
  )
  res <- multi_geo_autosplit_fetch("GSE2")
  expect_false(res$ok)
  expect_match(res$error, "resolved to exactly one")
})
