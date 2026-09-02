## Sources the real global.R once (all ~25 analysis packages + every
## load_default_*() function), then exercises every data-loading path each
## of the four omics modules actually uses at runtime - this is the "is
## every dataset each module needs actually present in data/" check, made
## repeatable instead of a one-off manual pass. Known-good shapes below were
## captured by hand against the migrated data and are pinned as regression
## guards, not just NULL-checks.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))

## ---- Transcriptomics ----------------------------------------------------

test_that("default transcriptomics dataset loads with the expected shape", {
  d <- load_default_dataset()
  expect_equal(dim(d$expr), c(15763L, 183L))
  expect_equal(nrow(d$meta), 183L)
})

test_that("transcriptomics precomputed WGCNA result loads (regression guard for the DATA_ROOT/data/processed bug)", {
  expect_true(file.exists(file.path(PROCESSED_DIR, "wgcna_results.rds")))
  w <- readRDS(file.path(PROCESSED_DIR, "wgcna_results.rds"))
  expect_true(is.list(w))
  net_files <- list.files(PROCESSED_DIR, pattern = "^wgcna_net_.*\\.rds$")
  expect_gte(length(net_files), 1)
})

test_that("transcriptomics genetics-submodule default objects load", {
  mr_primary <- readRDS(MR_PRIMARY_OBJECTS_RDS)
  expect_true(is.list(mr_primary))
  expect_true(all(c("primary", "res_or", "het", "pleio", "inst", "dat") %in% names(mr_primary)))
  expect_gt(nrow(mr_primary$primary), 0)
  expect_true(!is.null(readRDS(COLOC_REGIONS_RDS)))
  expect_true(!is.null(readRDS(VAL_SYNOVIUM_RDS)))
  expect_true(!is.null(readRDS(DGE_RESULTS_RDS)))
})

test_that("transcriptomics results/tables files used across submodules are readable", {
  for (f in c("MR35_crossancestry_female.csv", "MR35_crossancestry_male.csv",
              "WGCNA_05_gene_module_assignment.csv", "WGCNA_07_hub_genes_only.csv",
              "FS_venn_membership.csv",
              "FS_input_female.csv", "FS_input_female_noMHC.csv",
              "FS_input_male.csv", "FS_input_male_noMHC.csv",
              "MR_MHC_sensitivity_female.csv", "MR_MHC_sensitivity_male.csv",
              "candidates_female_disease.csv", "candidates_male_disease.csv",
              "DEG_female_full.csv", "DEG_male_full.csv", "DEG_all_full.csv")) {
    d <- read_table_safe(f)
    expect_true(!is.null(d), info = f)
    expect_gt(nrow(d), 0, label = f)
  }
  ## Excel exports are passed through via file.copy(), not fread() - just
  ## confirm the source files the download handler copies actually exist.
  expect_true(file.exists(file.path(TABLES_DIR, "MR_female_all_tables.xlsx")))
  expect_true(file.exists(file.path(TABLES_DIR, "MR_male_all_tables.xlsx")))
})

test_that("gene panels are discoverable and loadable", {
  panels <- list_gene_panels()
  expect_gte(length(panels), 1)
  genes <- load_gene_panel(panels[[1]])
  expect_gt(length(genes), 0)
})

test_that("project methodology lookup (ArthOChat tool) finds real content", {
  txt <- project_methods("wgcna")
  expect_true(grepl("WGCNA", txt))
  expect_false(grepl("No module matched", txt))
})

## ---- Methylomics ----------------------------------------------------

test_that("methylomics DMP (plain + SVA) and DMR load with the expected row counts", {
  dmp_plain <- load_default_dmp("plain", "female")
  dmp_sva   <- load_default_dmp("sva", "female")
  dmr       <- load_default_dmr("female")
  expect_equal(nrow(dmp_plain), 412492L)
  expect_equal(nrow(dmp_sva), 412492L)
  expect_gt(nrow(dmr), 0)
})

test_that("methylomics QC/pheno tables load with the expected sample count", {
  pheno <- load_default_meth_pheno()
  sexcheck <- load_default_meth_qc_sexcheck()
  expect_equal(nrow(pheno), 689L)
  expect_gt(nrow(sexcheck), 0)
})

test_that("methylomics WGCNA reference tables load for both sexes, merged and unmerged", {
  for (sex in c("female", "male")) {
    for (merged in c(FALSE, TRUE)) {
      mt <- load_default_wgcna_module_trait(sex, merged)
      ma <- load_default_wgcna_module_assignment(sex, merged)
      expect_true(!is.null(mt), info = paste(sex, merged))
      expect_true(!is.null(ma), info = paste(sex, merged))
    }
  }
  expect_true(!is.null(load_default_dmr_biomarker_panel("female")))
  expect_true(!is.null(load_default_dmr_biomarker_panel("male")))
})

test_that("methylomics MR + coloc reference tables load", {
  expect_true(!is.null(load_default_mr_estimates("female")))
  expect_true(!is.null(load_default_mr_estimates("male")))
  expect_true(!is.null(load_default_mr_harmonised()))
  expect_true(!is.null(load_default_mr_steiger("female")))
  expect_true(!is.null(load_default_mr_instrument_counts()))
  expect_true(!is.null(load_default_meth_coloc_results()))
})

test_that("methylomics diagnostic classifier reference + live panels load", {
  expect_true(!is.null(load_default_diagnostic_ensemble_votes("female")))
  expect_true(!is.null(load_default_diagnostic_ensemble_votes("male")))
  expect_true(!is.null(load_default_diagnostic_panel_auc("female")))
  expect_true(!is.null(load_default_diagnostic_perprobe_auc("female")))

  tt <- load_default_diagnostic_train_test()
  expect_setequal(names(tt), c("internal", "external"))
})

test_that("methylomics live beta matrix loads with the expected shape", {
  mm <- load_default_meth_matrix()
  expect_equal(dim(mm$beta), c(412492L, 689L))
  expect_equal(nrow(mm$pheno), 689L)
})

test_that("methylomics biomarker card's cytoband file loads", {
  source(file.path(app_dir, "R", "methylomics", "15_Biomarker_Analysis", "mod_methyl_biomarkercard.R"), local = TRUE)
  cb <- bc_cytoband_hg19()
  expect_gt(nrow(cb), 0)
})

test_that("methylomics methodology lookup finds real content", {
  txt <- project_methods_methylomics("dmp")
  expect_false(grepl("No Methylomics module matched", txt))
})

## ---- Cross-Omics ----------------------------------------------------

test_that("every cross-omics registry table loads", {
  for (label in names(CX_TABLE_REGISTRY)) {
    d <- load_default_cx_table(label)
    expect_true(!is.null(d), info = label)
    expect_gt(nrow(d), 0, label = label)
  }
})

## ---- Multi-Omics ----------------------------------------------------
## multi_read_registry_table() lives in R/multiomics/ (sourced by Shiny's
## loadSupport() at app boot, not by global.R alone) - read the registry
## paths directly here since what's under test is data completeness, not
## that wrapper's own logic (already covered by the app-smoke/upload tests
## and the original module audit).

test_that("every multi-omics table registry entry is a readable, non-empty CSV", {
  for (label in names(MULTI_TABLE_REGISTRY)) {
    d <- as.data.frame(data.table::fread(MULTI_TABLE_REGISTRY[[label]], showProgress = FALSE))
    expect_gt(nrow(d), 0, label = label)
  }
})

test_that("every multi-omics DIABLO fit object is a readable mixOmics fit", {
  for (key in names(MULTI_DIABLO_FIT_REGISTRY)) {
    fit <- readRDS(MULTI_DIABLO_FIT_REGISTRY[[key]])
    expect_true(inherits(fit, "sgccda"), info = key)
  }
})
