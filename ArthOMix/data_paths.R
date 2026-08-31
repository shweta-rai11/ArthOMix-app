## data_paths.R
## Centralized data-path configuration for ArthOMix.
##
## Every dataset the app reads at runtime lives under ArthOMix/data/ - nothing
## outside this app directory is required. Layout:
##   data/preloaded/{transcriptomics,methylomics,cross_omics,multiomics}/
##   data/reference/     - small standalone reference files (cytoband, ...)
##   data/annotations/   - curated gene panels etc.
##   data/examples/      - reserved for a bundled example-upload dataset
##                          (none currently ships - see dataset_manifest.csv)
##   data/uploads/       - reserved; NOT written to at runtime. User uploads
##                          go through Shiny's own per-session tempfile
##                          mechanism (input$file$datapath), which is already
##                          correctly isolated per session - redirecting them
##                          into a shared folder here would only add a
##                          cross-session leak risk with no benefit.
##   data/.cache/         - regenerable WGCNA / probe-collapse caches, kept
##                          out of preloaded/ so they're clearly distinguished
##                          from bundled, versioned data.
##
## This file is sourced explicitly, by path, from the very first line of
## global.R. It deliberately lives at the app ROOT (sibling of global.R/ui.R/
## server.R), NOT inside R/ - Shiny's loadSupport() auto-sources every
## top-level R/*.R file a second time, in its own child environment, and a
## file in R/ would end up sourced twice.
##
## Every constant/registry below keeps its ORIGINAL name (previously defined
## inline in global.R, pointing at sibling folders one level above ArthOMix/)
## so every load_default_*() function and every direct-by-name reference
## elsewhere in this app keeps working unchanged - only the paths they
## resolve to have moved.

DATA_DIR <- normalizePath(file.path(getwd(), "data"), mustWork = FALSE)

get_data_path       <- function(...) file.path(DATA_DIR, ...)
get_preloaded_path  <- function(...) file.path(DATA_DIR, "preloaded", ...)
get_reference_path  <- function(...) file.path(DATA_DIR, "reference", ...)
get_annotation_path <- function(...) file.path(DATA_DIR, "annotations", ...)
get_example_path    <- function(...) file.path(DATA_DIR, "examples", ...)
get_upload_path      <- function(...) file.path(DATA_DIR, "uploads", ...)

## ---------------------------------------------------------------------------
## Transcriptomics preloaded data
## ---------------------------------------------------------------------------
DATA_ROOT <- get_preloaded_path("transcriptomics")

if (!dir.exists(DATA_ROOT)) {
  stop(
    "Cannot find ArthOMix/data/preloaded/transcriptomics/. This app expects ",
    "all preloaded data to be bundled inside ArthOMix/data/ - see ",
    "data/dataset_manifest.csv, or re-run the data migration."
  )
}

TABLES_DIR  <- file.path(DATA_ROOT, "results", "tables")
FIGURES_DIR <- file.path(DATA_ROOT, "results", "figures", "by_section")
PROCESSED_DIR <- file.path(DATA_ROOT, "processed")
PROCESSED_NEW_DIR <- file.path(PROCESSED_DIR, "new")
RAW_DIR <- file.path(DATA_ROOT, "raw")

GENE_PANELS_DIR <- get_annotation_path("gene_panels")

DEFAULT_EXPR_RDS <- file.path(PROCESSED_DIR, "combined_expr_batchcorrected.rds")
DEFAULT_META_CSV <- file.path(PROCESSED_DIR, "combined_meta.csv")
MR_INSTRUMENTS_RDS <- file.path(PROCESSED_NEW_DIR, "MR_instruments.rds")
MR_PRIMARY_OBJECTS_RDS <- file.path(PROCESSED_NEW_DIR, "MR_primary_objects.rds")
COLOC_REGIONS_RDS <- file.path(PROCESSED_NEW_DIR, "coloc_regions.rds")
VAL_SYNOVIUM_RDS <- file.path(PROCESSED_NEW_DIR, "val_synovium.rds")
DGE_RESULTS_RDS <- file.path(PROCESSED_DIR, "dge_results.rds")
MR35_CROSSANCESTRY_FEMALE_CSV <- "MR35_crossancestry_female.csv"
MR35_CROSSANCESTRY_MALE_CSV   <- "MR35_crossancestry_male.csv"

PROJECT_CHAPTER_MD <- file.path(DATA_ROOT, "Chapter_2_subchapter2_sexstratified.md")

## Regenerable caches - deliberately outside data/preloaded/. These folders
## hold derivatives of USER-UPLOADED data (get_collapsed_genes(),
## get_or_compute_wgcna_blocks()/get_or_compute_meth_wgcna_blocks() in
## global.R all saveRDS() their result here, keyed by a content digest, with
## no expiry), so two things matter for a shared/multi-user deployment:
##   1. Created with owner-only permissions (mode 0700 - ignored on Windows,
##      a no-op there rather than an error) rather than the process umask
##      default, so other local accounts on a shared machine can't read
##      cached derivatives of another user's uploaded expression/methylation
##      data.
##   2. Swept for stale files at every app startup (arthomix_cleanup_stale_
##      cache() below) so this doesn't grow unbounded forever - there is no
##      other cleanup path anywhere in the app.
## This is deliberately just a TTL sweep, not a change to the cache
## key/lookup logic in global.R.

## Best-effort removal of cache files older than `max_age_days`. A single
## unreadable/unremovable file (permissions, concurrent access from another
## session) is skipped rather than aborting the whole sweep - this is
## opportunistic housekeeping, not a correctness-critical path.
arthomix_cleanup_stale_cache <- function(dir, max_age_days = 14) {
  if (!dir.exists(dir)) return(invisible(NULL))
  files <- list.files(dir, full.names = TRUE, recursive = TRUE, no.. = TRUE)
  if (length(files) == 0) return(invisible(NULL))
  cutoff <- Sys.time() - as.numeric(max_age_days) * 86400
  for (f in files) {
    info <- tryCatch(file.info(f), error = function(e) NULL)
    if (is.null(info) || is.na(info$mtime) || isTRUE(info$isdir)) next
    if (info$mtime < cutoff) {
      tryCatch(unlink(f), error = function(e) NULL)
    }
  }
  invisible(NULL)
}

COLLAPSED_CACHE_DIR <- get_data_path(".cache", "transcriptomics")
dir.create(COLLAPSED_CACHE_DIR, showWarnings = FALSE, recursive = TRUE, mode = "0700")
arthomix_cleanup_stale_cache(COLLAPSED_CACHE_DIR)
WGCNA_CACHE_DIR <- get_data_path(".cache", "transcriptomics", "wgcna")
dir.create(WGCNA_CACHE_DIR, showWarnings = FALSE, recursive = TRUE, mode = "0700")
arthomix_cleanup_stale_cache(WGCNA_CACHE_DIR)

## ---------------------------------------------------------------------------
## Methylomics: preloaded pipeline result tables
## ---------------------------------------------------------------------------
METH_DATA_ROOT <- get_preloaded_path("methylomics", "tables")
METH_DATA_AVAILABLE <- dir.exists(METH_DATA_ROOT)

METH_QC_PHENO_CSV      <- file.path(METH_DATA_ROOT, "script01_dataload_QC", "tables", "sample_metadata_qc.csv")
METH_QC_PCA_SEXCHECK_CSV <- file.path(METH_DATA_ROOT, "script01_dataload_QC", "tables", "qc_pca_sexcheck.csv")
METH_DMP_PLAIN_DIR  <- file.path(METH_DATA_ROOT, "script03_dmp_sexstratified", "tables")
METH_DMP_SVA_DIR    <- file.path(METH_DATA_ROOT, "script03_dmp_sva_sexstratified", "tables")
METH_DMR_DIR        <- file.path(METH_DATA_ROOT, "script04_dmr_sexstratified", "tables")
METH_WGCNA_DIR       <- file.path(METH_DATA_ROOT, "script05_wgcna_sexstratified", "tables")
METH_DIAGNOSTIC_VOTES_DIR <- file.path(METH_DATA_ROOT, "script07_ml_feature_selection", "tables")
METH_MR_DIR          <- file.path(METH_DATA_ROOT, "script08_mendelian_randomization", "tables")
METH_DIAGNOSTIC_DIR  <- file.path(METH_DATA_ROOT, "script09_diagnostic_classifier", "tables")

## ---------------------------------------------------------------------------
## Methylomics: the actual preloaded beta matrix (live QC, not just
## reproduced tables) + diagnostic-classifier train/test panels
## ---------------------------------------------------------------------------
METH_RAW_DATA_ROOT <- normalizePath(
  Sys.getenv("METH_RAW_DATA_ROOT", get_preloaded_path("methylomics", "matrix")),
  mustWork = FALSE
)
METH_BETA_RAW_RDS  <- file.path(METH_RAW_DATA_ROOT, "beta_raw.rds")
METH_PHENO_RDS     <- file.path(METH_RAW_DATA_ROOT, "pheno.rds")
METH_RAW_DATA_AVAILABLE <- file.exists(METH_BETA_RAW_RDS) && file.exists(METH_PHENO_RDS)

METH_DIAG_INTERNAL_RDS <- file.path(METH_RAW_DATA_ROOT, "gse42861_internal_panel_celltype_adjusted.rds")
METH_DIAG_EXTERNAL_RDS <- file.path(METH_RAW_DATA_ROOT, "gse111942_external_panel.rds")
METH_DIAG_DATA_AVAILABLE <- file.exists(METH_DIAG_INTERNAL_RDS) && file.exists(METH_DIAG_EXTERNAL_RDS)

METH_WGCNA_CACHE_DIR <- get_data_path(".cache", "methylomics", "wgcna")
dir.create(METH_WGCNA_CACHE_DIR, showWarnings = FALSE, recursive = TRUE, mode = "0700")
arthomix_cleanup_stale_cache(METH_WGCNA_CACHE_DIR)

## ---------------------------------------------------------------------------
## Cross-Omics: precomputed biomarker-convergence / cross-omics-MR tables
## ---------------------------------------------------------------------------
CX_DATA_ROOT <- get_preloaded_path("cross_omics")
CX_DATA_AVAILABLE <- dir.exists(CX_DATA_ROOT)
CX_RESULTS_DIR <- file.path(CX_DATA_ROOT, "tables")

CX_TABLE_REGISTRY <- list(
  "Combined convergence table (all layers, both sexes)" = file.path(CX_RESULTS_DIR, "MASTER_cross_omics_all_layers.csv"),
  "Tier 1 candidate biomarkers" = file.path(CX_RESULTS_DIR, "TIER1_MASTER_verified.csv"),
  "eQTL x mQTL combined" = file.path(CX_RESULTS_DIR, "cross_omics_eQTL_mQTL_combined.csv"),
  "eQTL x mQTL - female" = file.path(CX_RESULTS_DIR, "cross_omics_eQTL_mQTL_female.csv"),
  "eQTL x mQTL - male" = file.path(CX_RESULTS_DIR, "cross_omics_eQTL_mQTL_male.csv"),
  "Cross-omics MR - per-gene summary" = file.path(CX_RESULTS_DIR, "mr_stage_cross_omics_per_gene_summary.csv"),
  "Cross-omics MR - credible genes (full)" = file.path(CX_RESULTS_DIR, "mr_stage_cross_omics_credible_genes_full.csv"),
  "Cross-omics MR - eQTL-significant genes' mQTL instruments" = file.path(CX_RESULTS_DIR, "mr_stage_eqtl_significant_genes_mqtl_mr.csv")
)

## ---------------------------------------------------------------------------
## Multi-Omics: DIABLO/SNF pipeline tables + fits
## ---------------------------------------------------------------------------
MULTI_DATA_ROOT <- get_preloaded_path("multiomics")
MULTI_DATA_AVAILABLE <- dir.exists(MULTI_DATA_ROOT)
MULTI_RESULTS_ROOT <- file.path(MULTI_DATA_ROOT, "fits")
MULTI_RESULTS_DIR <- file.path(MULTI_DATA_ROOT, "tables")
MULTI_METADATA_DIR <- file.path(MULTI_DATA_ROOT, "tables")
MULTI_SUMMARY_DIR <- file.path(MULTI_DATA_ROOT, "tables", "summary")
MULTI_ADA_DIR <- file.path(MULTI_DATA_ROOT, "tables", "adalimumab")
MULTI_ETN_DIR <- file.path(MULTI_DATA_ROOT, "tables", "etanercept")

MULTI_TABLE_REGISTRY <- list(
  "Patient sample matching (all 80 patients)" = file.path(MULTI_METADATA_DIR, "patient_sample_matching_table.csv"),
  "RNA-seq QC summary" = file.path(MULTI_RESULTS_DIR, "RNAseq_QC_summary.csv"),
  "Methylation QC summary" = file.path(MULTI_RESULTS_DIR, "Methylation_QC_summary.csv"),
  "Benchmark vs. published (leakage-safe nested CV, Models A/B/C)" = file.path(MULTI_RESULTS_DIR, "Table8_benchmark_vs_published.csv"),
  "DIABLO performance - drug x sex (response)" = file.path(MULTI_RESULTS_DIR, "Table29_diablo_drugsex_performance.csv"),
  "DIABLO scores - drug x sex (response)" = file.path(MULTI_RESULTS_DIR, "Table29b_diablo_drugsex_scores.csv"),
  "DIABLO panel - drug x sex (response)" = file.path(MULTI_RESULTS_DIR, "Table30_diablo_drugsex_panel.csv"),
  "DIABLO performance - response (drug-pooled)" = file.path(MULTI_RESULTS_DIR, "Table34_diablo_response_sexstratified_performance.csv"),
  "DIABLO scores - response (drug-pooled)" = file.path(MULTI_RESULTS_DIR, "Table34b_diablo_response_sexstratified_scores.csv"),
  "DIABLO panel - response (drug-pooled)" = file.path(MULTI_RESULTS_DIR, "Table35_diablo_response_sexstratified_panel.csv"),
  "SNF performance - drug x sex" = file.path(MULTI_RESULTS_DIR, "Table22_snf_integration_performance.csv"),
  "Random Forest performance - drug x sex" = file.path(MULTI_RESULTS_DIR, "Table37_rf_drugsex_performance.csv"),
  ## Table38-41 - the pipeline's own RF/drug-type-DIABLO tables that never
  ## made it into earlier deployments (raw-CSV browse-only, via the Dataset
  ## Workspace's registry table picker - the live equivalents are on
  ## Integration/Biomarker Discovery's "Sex-Stratified" tab). Table40 here
  ## is deliberately distinctly labeled from the "Candidate multi-omics
  ## biomarkers" Table40 below - same pipeline table NUMBER, different file,
  ## different analysis (drug-type-outcome DIABLO vs. candidate-biomarker
  ## convergence).
  "Random Forest performance - single-omics, sex-pooled" = file.path(MULTI_RESULTS_DIR, "Table38_rf_single_omics_sexpooled_performance.csv"),
  "Random Forest performance - response (drug-pooled, sex-stratified)" = file.path(MULTI_RESULTS_DIR, "Table39_rf_response_sexstratified_performance.csv"),
  "DIABLO performance - drug type (sex-stratified)" = file.path(MULTI_RESULTS_DIR, "Table40_diablo_drugtype_sexstratified_performance.csv"),
  "DIABLO panel - drug type (sex-stratified)" = file.path(MULTI_RESULTS_DIR, "Table41_diablo_drugtype_sexstratified_panel.csv"),
  "SNF patient clusters - Adalimumab" = file.path(MULTI_ADA_DIR, "Table_SNFjoint_cluster_assignments_adalimumab.csv"),
  "SNF cluster-response association - Adalimumab" = file.path(MULTI_ADA_DIR, "Table_SNFjoint_cluster_response_association.csv"),
  "SNF concordance NMI - Adalimumab" = file.path(MULTI_ADA_DIR, "Table_SNFjoint_concordanceNMI_adalimumab.csv"),
  "SNF patient clusters - Etanercept" = file.path(MULTI_ETN_DIR, "Table_SNFjoint_cluster_assignments_etanercept.csv"),
  "SNF cluster-response association - Etanercept" = file.path(MULTI_ETN_DIR, "Table_SNFjoint_cluster_response_association.csv"),
  "SNF concordance NMI - Etanercept" = file.path(MULTI_ETN_DIR, "Table_SNFjoint_concordanceNMI_etanercept.csv"),
  "Master six-part summary (integrated vs single-omics)" = file.path(MULTI_SUMMARY_DIR, "Table36_master_six_part_summary.csv"),
  "Candidate multi-omics biomarkers - drug x sex (Etanercept panel)" = file.path(MULTI_SUMMARY_DIR, "Table40_candidate_multiomics_biomarkers_male_female.csv"),
  "Candidate multi-omics biomarkers - response (drug-pooled)" = file.path(MULTI_SUMMARY_DIR, "Table44b_candidate_multiomics_biomarkers_response_male_female.csv"),
  "Gene <-> CpG concordance - drug x sex (Etanercept panel)" = file.path(MULTI_SUMMARY_DIR, "Table42_gene_cpg_concordance_male_female_ETN.csv"),
  "Gene <-> CpG concordance - response (drug-pooled)" = file.path(MULTI_SUMMARY_DIR, "Table45_gene_cpg_concordance_male_female_response.csv"),
  "Pathway enrichment - drug x sex (Etanercept panel)" = file.path(MULTI_SUMMARY_DIR, "Table43_pathway_enrichment_male_female_ETN_panels.csv"),
  ## Genome-wide response-driven DEG/DMP (Table3/4, script 05) restricted to
  ## the gene/CpG IDs that actually appear in the candidate-biomarker/
  ## concordance tables above - the pipeline's own real per-feature nominal
  ## p_response/fdr_response (and logFC_response/delta_M_response), not a
  ## small-panel-recomputed statistic. Table4 full (~1.7M CpG rows, 230MB) is
  ## never shipped or loaded whole; this is a pre-filtered lookup built once
  ## against the ~1,600 candidate CpGs/genes actually used by this app.
  "Genome-wide DEG lookup - response-driven, by sex (candidates only)" = file.path(MULTI_SUMMARY_DIR, "Table3_response_driven_DEG_candidates.csv"),
  "Genome-wide DMP lookup - response-driven, by sex (candidates only)" = file.path(MULTI_SUMMARY_DIR, "Table4_response_driven_DMP_candidates.csv")
)

MULTI_DIABLO_FIT_REGISTRY <- list(
  female_Adalimumab = file.path(MULTI_RESULTS_ROOT, "diablo_drugsex_FEMALE_Adalimumab_fit.rds"),
  male_Adalimumab = file.path(MULTI_RESULTS_ROOT, "diablo_drugsex_MALE_Adalimumab_fit.rds"),
  female_Etanercept = file.path(MULTI_RESULTS_ROOT, "diablo_drugsex_FEMALE_Etanercept_fit.rds"),
  male_Etanercept = file.path(MULTI_RESULTS_ROOT, "diablo_drugsex_MALE_Etanercept_fit.rds"),
  female_response = file.path(MULTI_RESULTS_ROOT, "diablo_response_female_fit.rds"),
  male_response = file.path(MULTI_RESULTS_ROOT, "diablo_response_male_fit.rds")
)
