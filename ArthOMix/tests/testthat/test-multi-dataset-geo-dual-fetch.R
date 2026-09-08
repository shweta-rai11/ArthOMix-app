suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "multiomics", "06_Gene_CpG_Mapping", "multiomics_mapping_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "multiomics_dataset_helpers.R"))

fx_rnaseq_eset <- function(n_genes = 20, n_samples = 8, seed = 1, titles = NULL) {
  set.seed(seed)
  genes <- sprintf("GENE%d", seq_len(n_genes))
  samples <- sprintf("GSM%04d", 2000 + seq_len(n_samples))
  mat <- matrix(runif(n_genes * n_samples, 3, 12), n_genes, n_samples, dimnames = list(genes, samples))
  pdat <- data.frame(title = titles %||% sprintf("patient_%d", seq_len(n_samples)), geo_accession = samples,
                      row.names = samples, stringsAsFactors = FALSE)
  eset <- Biobase::ExpressionSet(assayData = mat, phenoData = Biobase::AnnotatedDataFrame(pdat))
  Biobase::annotation(eset) <- "GPL_RNA"
  eset
}

fx_methylation_eset <- function(n_probes = 20, n_samples = 8, seed = 2, titles = NULL) {
  set.seed(seed)
  probes <- sprintf("cg%08d", seq_len(n_probes))
  samples <- sprintf("GSM%04d", 3000 + seq_len(n_samples))
  mat <- matrix(runif(n_probes * n_samples, 0, 1), n_probes, n_samples, dimnames = list(probes, samples))
  pdat <- data.frame(title = titles %||% sprintf("patient_%d", seq_len(n_samples)), geo_accession = samples,
                      row.names = samples, stringsAsFactors = FALSE)
  eset <- Biobase::ExpressionSet(assayData = mat, phenoData = Biobase::AnnotatedDataFrame(pdat))
  Biobase::annotation(eset) <- "GPL_METH"
  eset
}

fx_empty_eset <- function(n_samples = 8) {
  samples <- sprintf("GSM%04d", 5000 + seq_len(n_samples))
  mat <- matrix(numeric(0), 0, n_samples, dimnames = list(NULL, samples))
  pdat <- data.frame(title = sprintf("patient_%d", seq_len(n_samples)), row.names = samples, stringsAsFactors = FALSE)
  eset <- Biobase::ExpressionSet(assayData = mat, phenoData = Biobase::AnnotatedDataFrame(pdat))
  Biobase::annotation(eset) <- "GPL_EMPTY"
  eset
}

## ---- multi_geo_derive_title_patient_num() ----

test_that("multi_geo_derive_title_patient_num() extracts a zero-padded trailing number from title", {
  meta <- data.frame(title = c("PBMC_E_n_01", "PBMC_A_m_02", "DNA_E_n_9"), stringsAsFactors = FALSE)
  out <- multi_geo_derive_title_patient_num(meta)
  expect_equal(out, c("01", "02", "09"))
})

test_that("multi_geo_derive_title_patient_num() returns NA for titles with no trailing digits", {
  meta <- data.frame(title = c("control_sample", "PBMC_07", "no_number_here"), stringsAsFactors = FALSE)
  out <- multi_geo_derive_title_patient_num(meta)
  expect_equal(out, c(NA_character_, "07", NA_character_))
})

test_that("multi_geo_derive_title_patient_num() returns all-NA when there is no title column, without erroring", {
  meta <- data.frame(other_col = c("a", "b"), stringsAsFactors = FALSE)
  out <- multi_geo_derive_title_patient_num(meta)
  expect_equal(out, c(NA_character_, NA_character_))
})

test_that("multi_geo_derive_title_patient_num() handles NULL metadata without erroring", {
  expect_equal(multi_geo_derive_title_patient_num(NULL), character(0))
})

## ---- multi_geo_dual_fetch() ----

test_that("multi_geo_dual_fetch() rejects malformed or identical accessions without any network call", {
  called <- FALSE
  testthat::local_mocked_bindings(getGEO = function(...) { called <<- TRUE; NULL }, .package = "GEOquery")
  res1 <- multi_geo_dual_fetch("not-an-accession", "GSE12345")
  expect_false(res1$ok)
  res2 <- multi_geo_dual_fetch("GSE12345", "GSE12345")
  expect_false(res2$ok)
  expect_match(res2$error, "different accessions")
  expect_false(called)
})

test_that("multi_geo_dual_fetch() succeeds and classifies correctly regardless of accession order", {
  testthat::local_mocked_bindings(
    getGEO = function(acc, ...) {
      acc <- toupper(acc)
      if (identical(acc, "GSE100")) list(GPL_RNA = fx_rnaseq_eset(n_samples = 6, titles = sprintf("PBMC_E_n_%02d", 1:6)))
      else if (identical(acc, "GSE200")) list(GPL_METH = fx_methylation_eset(n_samples = 6, titles = sprintf("DNA_E_n_%02d", 1:6)))
      else stop("unexpected accession in test: ", acc)
    },
    .package = "GEOquery"
  )
  res_ab <- multi_geo_dual_fetch("GSE100", "GSE200")
  res_ba <- multi_geo_dual_fetch("GSE200", "GSE100")

  expect_true(res_ab$ok); expect_true(res_ba$ok)
  expect_equal(res_ab$expression$accession, "GSE100")
  expect_equal(res_ab$methylation$accession, "GSE200")
  expect_equal(res_ba$expression$accession, "GSE100")
  expect_equal(res_ba$methylation$accession, "GSE200")
  expect_equal(dim(res_ab$expression$mat), c(6L, 20L))
  expect_equal(dim(res_ab$methylation$mat), c(6L, 20L))
})

test_that("multi_geo_dual_fetch() attaches geo_title_patient_num to both fetched layers, matched by shared trailing number", {
  testthat::local_mocked_bindings(
    getGEO = function(acc, ...) {
      acc <- toupper(acc)
      if (identical(acc, "GSE100")) list(GPL_RNA = fx_rnaseq_eset(n_samples = 4, titles = sprintf("PBMC_E_n_%02d", 1:4)))
      else list(GPL_METH = fx_methylation_eset(n_samples = 4, titles = sprintf("DNA_E_n_%02d", 1:4)))
    },
    .package = "GEOquery"
  )
  res <- multi_geo_dual_fetch("GSE100", "GSE200")
  expect_true(res$ok)
  expect_equal(res$expression$meta$geo_title_patient_num, c("01", "02", "03", "04"))
  expect_equal(res$methylation$meta$geo_title_patient_num, c("01", "02", "03", "04"))
})

test_that("multi_geo_dual_fetch() gives a specific, actionable error when a series has no values embedded (RNA-seq supplementary-file case)", {
  testthat::local_mocked_bindings(
    getGEO = function(acc, ...) {
      acc <- toupper(acc)
      if (identical(acc, "GSE100")) list(GPL_EMPTY = fx_empty_eset())
      else list(GPL_METH = fx_methylation_eset())
    },
    .package = "GEOquery"
  )
  res <- multi_geo_dual_fetch("GSE100", "GSE200")
  expect_false(res$ok)
  expect_match(res$error, "supplementary file")
  expect_match(res$error, "Upload Dataset")
})

test_that("multi_geo_dual_fetch() errors out when both accessions classify as the same omics type", {
  testthat::local_mocked_bindings(
    getGEO = function(acc, ...) list(GPL_RNA = fx_rnaseq_eset()),
    .package = "GEOquery"
  )
  res <- multi_geo_dual_fetch("GSE100", "GSE200")
  expect_false(res$ok)
  expect_match(res$error, "could not resolve", ignore.case = TRUE)
})

test_that("multi_geo_dual_fetch() surfaces which accession failed when only one of the two can be fetched", {
  testthat::local_mocked_bindings(
    getGEO = function(acc, ...) {
      acc <- toupper(acc)
      if (identical(acc, "GSE100")) stop("network error")
      else list(GPL_METH = fx_methylation_eset())
    },
    .package = "GEOquery"
  )
  res <- multi_geo_dual_fetch("GSE100", "GSE200")
  expect_false(res$ok)
  expect_match(res$error, "GSE100")
})
