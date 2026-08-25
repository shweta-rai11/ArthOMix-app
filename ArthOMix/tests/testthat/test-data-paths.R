## Fast, no-package-load tests: every path constant data_paths.R defines
## must resolve to something that actually exists under ArthOMix/data/.
## This is the fastest possible check that the app is fully self-contained -
## no external "../<folder>" dependency should ever appear here again.

test_that("DATA_DIR resolves under the app root, not an external folder", {
  expect_true(dir.exists(DATA_DIR))
  expect_true(startsWith(DATA_DIR, app_dir))
  expect_false(grepl("Research_Q[234]|Research_05_multiomics|(^|/)methylomics(/|$)", DATA_DIR))
})

test_that("transcriptomics preloaded roots exist", {
  expect_true(dir.exists(DATA_ROOT))
  expect_true(dir.exists(TABLES_DIR))
  expect_true(dir.exists(PROCESSED_DIR))
  expect_true(dir.exists(PROCESSED_NEW_DIR))
  expect_true(dir.exists(FIGURES_DIR)) # empty placeholder, but addResourcePath() requires it to exist
  expect_true(dir.exists(GENE_PANELS_DIR))
})

test_that("transcriptomics default dataset files exist", {
  expect_true(file.exists(DEFAULT_EXPR_RDS))
  expect_true(file.exists(DEFAULT_META_CSV))
  expect_true(file.exists(MR_PRIMARY_OBJECTS_RDS))
  expect_true(file.exists(COLOC_REGIONS_RDS))
  expect_true(file.exists(VAL_SYNOVIUM_RDS))
  expect_true(file.exists(DGE_RESULTS_RDS))
  expect_true(file.exists(PROJECT_CHAPTER_MD))
})

test_that("methylomics preloaded tables root and every script0N subfolder exist", {
  expect_true(METH_DATA_AVAILABLE)
  expect_true(dir.exists(METH_DATA_ROOT))
  for (d in c("script01_dataload_QC", "script03_dmp_sexstratified", "script03_dmp_sva_sexstratified",
              "script04_dmr_sexstratified", "script05_wgcna_sexstratified",
              "script07_ml_feature_selection", "script08_mendelian_randomization",
              "script09_diagnostic_classifier")) {
    expect_true(dir.exists(file.path(METH_DATA_ROOT, d)), info = d)
  }
  expect_true(file.exists(METH_QC_PHENO_CSV))
  expect_true(file.exists(METH_QC_PCA_SEXCHECK_CSV))
  expect_true(dir.exists(METH_DMP_PLAIN_DIR))
  expect_true(dir.exists(METH_DMP_SVA_DIR))
  expect_true(dir.exists(METH_DMR_DIR))
  expect_true(dir.exists(METH_WGCNA_DIR))
  expect_true(dir.exists(METH_DIAGNOSTIC_VOTES_DIR))
  expect_true(dir.exists(METH_MR_DIR))
  expect_true(dir.exists(METH_DIAGNOSTIC_DIR))
})

test_that("methylomics live beta matrix + diagnostic panels exist", {
  expect_true(METH_RAW_DATA_AVAILABLE)
  expect_true(file.exists(METH_BETA_RAW_RDS))
  expect_true(file.exists(METH_PHENO_RDS))
  expect_true(METH_DIAG_DATA_AVAILABLE)
  expect_true(file.exists(METH_DIAG_INTERNAL_RDS))
  expect_true(file.exists(METH_DIAG_EXTERNAL_RDS))
})

test_that("cross-omics table registry paths all exist", {
  expect_true(CX_DATA_AVAILABLE)
  for (label in names(CX_TABLE_REGISTRY)) {
    expect_true(file.exists(CX_TABLE_REGISTRY[[label]]), info = label)
  }
})

test_that("multi-omics table + fit registry paths all exist", {
  expect_true(MULTI_DATA_AVAILABLE)
  for (label in names(MULTI_TABLE_REGISTRY)) {
    expect_true(file.exists(MULTI_TABLE_REGISTRY[[label]]), info = label)
  }
  for (key in names(MULTI_DIABLO_FIT_REGISTRY)) {
    expect_true(file.exists(MULTI_DIABLO_FIT_REGISTRY[[key]]), info = key)
  }
})

test_that("reference and annotation files exist", {
  expect_true(file.exists(get_reference_path("cytoBandIdeo_hg19.txt.gz")))
  expect_true(length(list.files(GENE_PANELS_DIR, pattern = "\\.txt$")) >= 1)
})

test_that("regenerable cache dirs were created under data/.cache, not data/preloaded", {
  expect_true(dir.exists(COLLAPSED_CACHE_DIR))
  expect_true(dir.exists(WGCNA_CACHE_DIR))
  expect_true(dir.exists(METH_WGCNA_CACHE_DIR))
  expect_true(grepl("/data/\\.cache/", COLLAPSED_CACHE_DIR))
  expect_true(grepl("/data/\\.cache/", WGCNA_CACHE_DIR))
  expect_true(grepl("/data/\\.cache/", METH_WGCNA_CACHE_DIR))
})

test_that("upload test fixtures (data/examples/) exist", {
  merged_dir <- get_example_path("transcriptomics_upload", "merged")
  probelevel_dir <- get_example_path("transcriptomics_upload", "probelevel")
  expect_true(file.exists(file.path(merged_dir, "chen2021_merged_expression_matrix.csv")))
  expect_true(file.exists(file.path(merged_dir, "chen2021_merged_sample_metadata.csv")))
  expect_true(dir.exists(probelevel_dir))
  expect_gte(length(list.files(probelevel_dir, pattern = "\\.csv$")), 7)
})
