## global.R
## ArthOMix Explorer
## Shared setup: packages, data paths, and the module/submodule registry.

## Shiny's own default (5MB) is too small for a real expression matrix
## upload (e.g. a genome-wide RNA-seq counts CSV) - raise it for every
## fileInput in the app (Dataset tab, Preprocessing, Feature Selection's
## own upload path, etc.), not just one.
options(shiny.maxRequestSize = 200 * 1024^2)

## Background execution for genuinely slow, blocking operations - currently
## only the Methylomics Dataset tab's preloaded ~2.1GB live beta matrix read
## (mod_methyl_dataset.R). Shiny is single-threaded per R process by
## default, so without this, that one readRDS() call froze the *entire*
## app - every tab, every open browser session - for however long the read
## took (often 30s-2min+). shiny::ExtendedTask (>=1.8) is the supported way
## to run such a call in the background while keeping even the session that
## triggered it fully responsive; it needs a `future` plan set once,
## process-wide. Guarded by requireNamespace() since `future`/`promises`
## aren't otherwise required by this app - if either is missing, the
## affected loads simply fall back to blocking (see mod_methyl_dataset.R).
ARTHOMIX_ASYNC_AVAILABLE <- requireNamespace("future", quietly = TRUE) && requireNamespace("promises", quietly = TRUE)
if (ARTHOMIX_ASYNC_AVAILABLE) {
  future::plan(future::multisession, workers = 2)
}

library(shiny)
library(shinydashboard)
library(shinythemes)
library(shinyWidgets)
library(shinyjs)
library(shinycssloaders)
library(DT)
library(ggplot2)
library(ggrepel)
library(plotly)
library(dplyr)
library(tidyr)
library(data.table)

## Analysis packages, loaded once here so every submodule can use them
## without repeating library() calls.
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(WGCNA)
  library(glmnet)
  library(randomForest)
  library(e1071)
  library(caret)
  library(pROC)
  library(sva)
  library(VennDiagram)
  library(ggVennDiagram)
  library(MCPcounter)
  ## IOBR: needed for CIBERSORT (LM22) deconvolution in mod_deconvolution.R -
  ## must be library()-attached, not just called via IOBR:: - deconvo_cibersort()
  ## looks up the bundled lm22 reference by bare name, which only resolves once
  ## the package's lazy data is attached, not merely namespace-loaded.
  library(IOBR)
  library(rms)
  library(MendelianRandomization)
  library(coloc)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  ## igraph/ggraph: render the WGCNA hub-gene network as an in-app figure
  ## (mod_wgcna.R's Step 5). igraph::union()/groups()/simplify()/etc. mask
  ## same-named base/dplyr functions once attached, but checked directly:
  ## igraph's union() falls back to base::union() behavior on plain vectors
  ## (confirmed against mod_featureselection.R's own union(lasso_genes,
  ## rf_genes) call), and none of the others are called unqualified anywhere
  ## in this app.
  library(igraph)
  library(ggraph)
})

## WGCNA's correlation/TOM computation (mod_wgcna.R's blockwiseModules and
## pickSoftThreshold calls) is single-threaded by default - on a full,
## unfiltered gene set (this project's own default, all 15,763 network
## genes rather than a variance-filtered subset) that's a large, slow
## computation on one core even though the machine typically has several
## idle. Enabled once per R session, not per-call, since that's what
## WGCNA::enableWGCNAThreads() is designed for. Leaves 2 cores free for the
## Shiny process itself and the OS rather than claiming every core.
suppressMessages(
  WGCNA::enableWGCNAThreads(nThreads = max(1, parallel::detectCores() - 2))
)

## AI research assistant: chat client + Shiny chat UI. Kept separate from
## the analysis packages above since they're not statistical tooling.
library(ellmer)
library(shinychat)

## Two of the above packages define functions that mask ones other
## submodules depend on, once every submodule file is sourced into the
## global environment (an unqualified call in that code follows the normal
## search() path, not a package namespace):
##  - rms defines its own validate(), masking shiny::validate(), which every
##    submodule relies on for validate(need(...)) input checks.
##  - org.Hs.eg.db/clusterProfiler pull in Bioconductor's IRanges/S4Vectors,
##    which define their own cor(), masking WGCNA::cor() (needed internally
##    by WGCNA::blockwiseModules) since they attach after WGCNA.
## Restore both as the versions found by unqualified calls.
validate <- shiny::validate
cor <- WGCNA::cor

## ---------------------------------------------------------------------------
## Data location
## ---------------------------------------------------------------------------
## The app reads precomputed pipeline outputs directly from the research
## project folder instead of duplicating them. If that folder is renamed or
## moved, update PROJECT_DIR_NAME below.

PROJECT_DIR_NAME <- "Research_Q2_TRANSCRIPTOMICS_sexstratified_COPY"
DATA_ROOT <- normalizePath(file.path("..", PROJECT_DIR_NAME), mustWork = FALSE)

if (!dir.exists(DATA_ROOT)) {
  stop(
    "Cannot find the data folder '", PROJECT_DIR_NAME, "' next to shiny_app/. ",
    "Keep shiny_app/ and ", PROJECT_DIR_NAME, "/ inside the same parent folder, ",
    "or edit PROJECT_DIR_NAME in global.R."
  )
}

TABLES_DIR  <- file.path(DATA_ROOT, "results", "tables")
FIGURES_DIR <- file.path(DATA_ROOT, "results", "figures", "by_section")
PROCESSED_DIR <- file.path(DATA_ROOT, "data", "processed")
PROCESSED_NEW_DIR <- file.path(PROCESSED_DIR, "new")

## Bundled reference gene panels (one gene symbol per line, .txt) that any
## submodule can offer as an optional third set to intersect against - e.g.
## Candidate Gene Identification narrowing a WGCNA-module/DEG overlap down
## to a disease-process-specific panel, the way a published study's own
## curated gene list (ferroptosis, autophagy, whatever) would. Not
## hardcoded to any one panel: drop any "<label>.txt" file into this folder
## and it's picked up automatically, one gene symbol per line.
GENE_PANELS_DIR <- file.path(DATA_ROOT, "data", "gene_panels")

list_gene_panels <- function() {
  files <- list.files(GENE_PANELS_DIR, pattern = "\\.txt$", full.names = TRUE)
  if (length(files) == 0) return(character(0))
  labels <- tools::toTitleCase(gsub("_", " ", tools::file_path_sans_ext(basename(files))))
  stats::setNames(files, labels)
}

load_gene_panel <- function(path) {
  genes <- trimws(readLines(path, warn = FALSE))
  unique(genes[nzchar(genes)])
}

## Serve the precomputed figures as static files at /figures/<name>.png
addResourcePath("figures", FIGURES_DIR)

## ---------------------------------------------------------------------------
## Default dataset and saved objects
## ---------------------------------------------------------------------------
## The app opens with the user's own cohort loaded. Any submodule that needs
## it reads from the shared `dataset` reactiveValues (set up in server.R),
## not from these paths directly, so a user upload on the Dataset tab
## transparently replaces it everywhere.

DEFAULT_EXPR_RDS <- file.path(PROCESSED_DIR, "combined_expr_batchcorrected.rds")
DEFAULT_META_CSV <- file.path(PROCESSED_DIR, "combined_meta.csv")

## Genetics submodules (MR, colocalisation) need exposure/outcome summary
## statistics that an expression matrix cannot supply. These default objects
## let those modules run live out of the box; a matching upload replaces them.
MR_INSTRUMENTS_RDS   <- file.path(PROCESSED_NEW_DIR, "MR_instruments.rds")
## The project's own cached end-to-end MR run (scripts/00_shared/10_MR.R):
## already cis-restricted (+/-1Mb of the gene body, GRCh37), MHC-flagged, and
## harmonised (action=2) against the Okada 2014 RA GWAS - covering 1,701
## genes with a usable (mr_keep) instrument, not just the 33 mod_mr.R could previously reach through
## coloc_regions.rds. mod_mr.R's live advanced filters (p-value, F-statistic,
## MHC handling) subset this cached object rather than re-deriving cis
## restriction from scratch, since that requires the EnsDb.Hsapiens.v75
## GRCh37 annotation this app does not otherwise load.
MR_PRIMARY_OBJECTS_RDS <- file.path(PROCESSED_NEW_DIR, "MR_primary_objects.rds")
COLOC_REGIONS_RDS    <- file.path(PROCESSED_NEW_DIR, "coloc_regions.rds")
VAL_SYNOVIUM_RDS     <- file.path(PROCESSED_NEW_DIR, "val_synovium.rds")
## Sex-stratified blood differential-expression results (limma, RA vs HC) -
## this project's own dge_results.rds (scripts/00_shared, consumed by
## 20_testing_synovium_external.R as `D`). mod_crosstissue.R reads this as
## the BUNDLED fallback for a panel gene's blood direction of effect (sign of
## logFC), used to assess direction concordance with synovium - the live
## alternative is whichever of this session's own DGE runs
## (results$dge_runs, from mod_dge.R) matches that sex.
DGE_RESULTS_RDS      <- file.path(PROCESSED_DIR, "dge_results.rds")
MR35_CROSSANCESTRY_FEMALE_CSV <- "MR35_crossancestry_female.csv"
MR35_CROSSANCESTRY_MALE_CSV   <- "MR35_crossancestry_male.csv"

load_default_dataset <- function() {
  expr <- readRDS(DEFAULT_EXPR_RDS)
  meta <- as.data.frame(data.table::fread(DEFAULT_META_CSV, showProgress = FALSE))
  common <- intersect(colnames(expr), meta$sample)
  expr <- expr[, common, drop = FALSE]
  meta <- meta[match(common, meta$sample), , drop = FALSE]
  list(expr = expr, meta = meta, source = "Example dataset: sex-stratified RA blood cohort (GSE93272 + GSE110169)")
}

## ---------------------------------------------------------------------------
## Raw GEO sources
## ---------------------------------------------------------------------------
## The example cohort's two training datasets, and two more held out for
## validation, all as raw ExpressionSets - read once per session and
## cached, since these matrices are large and read-only. The Preprocessing
## submodule builds live, interactive views from these (and from your own
## uploads), rather than reading any precomputed pipeline output.

RAW_DIR <- file.path(DATA_ROOT, "data", "raw")

.arthomix_cache <- new.env(parent = emptyenv())

get_raw_eset <- function(gse_id) {
  key <- paste0("eset_", gse_id)
  if (is.null(.arthomix_cache[[key]])) {
    path <- file.path(RAW_DIR, paste0(gse_id, "_raw.rds"))
    if (!file.exists(path)) return(NULL)
    e <- readRDS(path)
    if (is(e, "list")) e <- e[[1]]
    .arthomix_cache[[key]] <- e
  }
  .arthomix_cache[[key]]
}

## Two-tier cache, same key (GSE ID) either way: an in-memory tier
## (.arthomix_cache, same as get_raw_eset() above) for repeat calls within
## one running process, and an on-disk RDS tier under data/cache/ so the
## expensive step - WGCNA::collapseRows over the full probe-level matrix,
## tens of thousands of rows - is only ever paid ONCE per dataset, ever,
## rather than once per fresh R session/app restart. mod_preprocessing.R
## calls this for the same two training GSEs from three different places
## (per-source preprocessing, the merge overlap Venn, the live example
## merge); without caching, a single click through that tab could pay this
## cost 3-4 times in a row, and since Shiny runs single-threaded, the whole
## app is unresponsive to every user for the entire time.
COLLAPSED_CACHE_DIR <- file.path(DATA_ROOT, "data", "cache")
dir.create(COLLAPSED_CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

## Same two-tier cache, generalised for mod_wgcna.R's blockwiseModules
## result (the single slowest step in this app, especially at the
## project's own default of "all genes", not a variance-filtered subset -
## see mod_wgcna.R Step 1/2/3). Keyed on a digest of the exact expression
## matrix plus every setting that changes the result, so ANY combination a
## user picks - a quick exploratory top-2000-gene run just as much as the
## full 15,763-gene default - is computed at most once per machine ever,
## no matter how many times "Detect modules" is clicked or the app is
## restarted. Nothing about which filters/settings are choosable changes;
## this only short-circuits recomputing a combination already paid for.
WGCNA_CACHE_DIR <- file.path(DATA_ROOT, "data", "cache", "wgcna")
dir.create(WGCNA_CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

get_or_compute_wgcna_blocks <- function(key_parts, compute_fn) {
  cache_key <- paste0("wgcna_", digest::digest(key_parts, algo = "xxhash64"))
  if (!is.null(.arthomix_cache[[cache_key]])) return(.arthomix_cache[[cache_key]])
  disk_path <- file.path(WGCNA_CACHE_DIR, paste0(cache_key, ".rds"))
  result <- if (file.exists(disk_path)) {
    readRDS(disk_path)
  } else {
    result <- compute_fn()
    ## Best-effort: a read-only data folder shouldn't break the
    ## already-computed in-memory result for this session.
    tryCatch(saveRDS(result, disk_path), error = function(e) NULL)
    result
  }
  .arthomix_cache[[cache_key]] <- result
  result
}

get_collapsed_genes <- function(gse_id) {
  key <- paste0("collapsed_", gse_id)
  if (is.null(.arthomix_cache[[key]])) {
    disk_path <- file.path(COLLAPSED_CACHE_DIR, paste0(gse_id, "_collapsed.rds"))
    if (file.exists(disk_path)) {
      .arthomix_cache[[key]] <- readRDS(disk_path)
    } else {
      eset <- get_raw_eset(gse_id)
      if (is.null(eset)) return(NULL)
      collapsed <- collapse_probes_to_genes(eset)
      .arthomix_cache[[key]] <- collapsed
      ## Best-effort: a read-only data folder shouldn't break the
      ## already-computed in-memory result for this session.
      tryCatch(saveRDS(collapsed, disk_path), error = function(e) NULL)
    }
  }
  .arthomix_cache[[key]]
}

## ---------------------------------------------------------------------------
## Methylomics: preloaded default analysis (GSE42861, Liu et al. 2013)
## ---------------------------------------------------------------------------
## The methylomics counterpart to DATA_ROOT/load_default_dataset() above -
## Research_Q3_METHYLOMICS_sexstratified_COPY/methylomics/, a nested project
## (own .gitignore entry, not tracked by this repo) with one script0N_*/
## stage per analysis step, each with its own tables/figures and a
## METHODS_*.md documenting exactly what was run. Unlike the transcriptomics
## DATA_ROOT, this one is optional - the app must still run (Methylomics
## upload-only) if a user hasn't copied this folder in, so nothing here
## calls stop() the way DATA_ROOT's missing-folder check does.
##
## Only the already-computed result tables are read here, not the ~2.1GB
## QC'd beta matrix those tables were derived from (data/processed/, on the
## source machine only, excluded when this folder was copied in - see its
## own .gitignore comment) - every table below already carries genome-wide
## results for all 412,492 retained probes, so interactive FDR/delta-beta
## filtering works directly off them with no need to hold the full matrix
## in memory.
METH_DATA_ROOT <- normalizePath(file.path("..", "Research_Q3_METHYLOMICS_sexstratified_COPY", "methylomics"), mustWork = FALSE)
METH_DATA_AVAILABLE <- dir.exists(METH_DATA_ROOT)

METH_QC_PHENO_CSV      <- file.path(METH_DATA_ROOT, "script01_dataload_QC", "tables", "sample_metadata_qc.csv")
METH_QC_PCA_SEXCHECK_CSV <- file.path(METH_DATA_ROOT, "script01_dataload_QC", "tables", "qc_pca_sexcheck.csv")
METH_DMP_PLAIN_DIR  <- file.path(METH_DATA_ROOT, "script03_dmp_sexstratified", "tables")
METH_DMP_SVA_DIR    <- file.path(METH_DATA_ROOT, "script03_dmp_sva_sexstratified", "tables")
METH_DMR_DIR        <- file.path(METH_DATA_ROOT, "script04_dmr_sexstratified", "tables")

## Sequential probe-filtering cascade (METHODS_load_qc.md Table 2.X.3) -
## the actual logged counts from the completed run, not recomputed here
## (recomputing it needs the ~2.1GB QC'd beta matrix, excluded when this
## folder was copied in - see METH_DATA_ROOT's own comment above).
METH_QC_PROBE_CASCADE <- data.frame(
  step = c("Raw (post-parse)", "cg-prefix only", "MASK_general removal", "Multi-hit removal", "Sex-chromosome removal", "Missingness > 5% removal"),
  retained = c(485577L, 482421L, 422520L, 422520L, 412492L, 412492L),
  removed  = c(NA_integer_, 3156L, 59901L, 0L, 10028L, 0L)
)

## Reads one sex-stratum's DMP results for one pipeline stage ("plain" - the
## unadjusted model, which found zero genome-wide-significant CpGs due to
## residual inflation - or "sva", the surrogate-variable-adjusted follow-up
## that resolves it; see METHODS_dmp_sexstratified.md /
## METHODS_dmp_sva_sexstratified.md). Returns NULL rather than erroring if
## METH_DATA_AVAILABLE is FALSE, so mod_methyl_dmp.R can show an empty
## state instead of crashing when this folder hasn't been copied in.
load_default_dmp <- function(stage = c("plain", "sva"), sex = c("female", "male")) {
  stage <- match.arg(stage); sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  dir <- if (stage == "plain") METH_DMP_PLAIN_DIR else METH_DMP_SVA_DIR
  path <- file.path(dir, sprintf("dmp_%s_full.csv", sex))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

## Reads one sex-stratum's precomputed DMR (region-level) results
## (script04_dmr_sexstratified/tables/dmr_{sex}_full.csv) - DMRcate region
## calls (lambda=1000, C=2) on the SVA-adjusted, bacon-corrected per-CpG
## statistics from load_default_dmp("sva", sex), with an explicit
## Benjamini-Hochberg region-level FDR (dmr_fdr, on the Stouffer statistic)
## already computed; see METHODS_dmr_sexstratified.md. Returns NULL rather
## than erroring if METH_DATA_AVAILABLE is FALSE, mirroring load_default_dmp().
load_default_dmr <- function(sex = c("female", "male")) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_DMR_DIR, sprintf("dmr_%s_full.csv", sex))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

## The QC-derived sample metadata (690 rows incl. header - 689 GSE42861
## samples) - group counts, sex, PCA outlier flags, chrY sex-check - for
## the dataset-context panel above the DMP results, not full QC replication.
load_default_meth_pheno <- function() {
  if (!METH_DATA_AVAILABLE || !file.exists(METH_QC_PHENO_CSV)) return(NULL)
  as.data.frame(data.table::fread(METH_QC_PHENO_CSV, showProgress = FALSE))
}

## Per-sample PCA coordinates, MAD-outlier flag, and chrY-based sex-check
## result (script01_dataload_QC/tables/qc_pca_sexcheck.csv) - the QC
## module's default-analysis view reads this directly rather than the
## summary already merged into load_default_meth_pheno(), since this one
## additionally retains the k-means cluster assignment and predicted sex.
load_default_meth_qc_sexcheck <- function() {
  if (!METH_DATA_AVAILABLE || !file.exists(METH_QC_PCA_SEXCHECK_CSV)) return(NULL)
  as.data.frame(data.table::fread(METH_QC_PCA_SEXCHECK_CSV, showProgress = FALSE))
}

## ---------------------------------------------------------------------------
## Methylomics: the actual preloaded beta matrix (live QC, not just reproduced
## tables)
## ---------------------------------------------------------------------------
## METH_DATA_ROOT above only carries the completed pipeline's OUTPUT tables -
## the ~2.1GB QC'd beta matrix that produced them was deliberately excluded
## when that folder was copied in. It still exists at its original location
## on the source machine (outside this repo's directory tree entirely, not
## just excluded via .gitignore within it - unlike METH_DATA_ROOT/DATA_ROOT,
## which are both a sibling of shiny_app/ and so reachable by a single "..",
## this one has no such relative relationship to shiny_app/'s location).
## Read directly from where it already lives instead of duplicating a
## multi-gigabyte file into this repo, same idea as RAW_DIR above. Optional,
## like METH_DATA_ROOT - the app must still run (upload-only Methylomics QC)
## on a machine that only has the METH_DATA_ROOT results copy; set the
## METH_RAW_DATA_ROOT environment variable to point at this folder on any
## other machine.
METH_RAW_DATA_ROOT <- normalizePath(
  Sys.getenv("METH_RAW_DATA_ROOT", "/Users/swetarai/ArthOMix/methylomics/data/processed"),
  mustWork = FALSE
)
METH_BETA_RAW_RDS  <- file.path(METH_RAW_DATA_ROOT, "beta_raw.rds")
METH_PHENO_RDS     <- file.path(METH_RAW_DATA_ROOT, "pheno.rds")
METH_RAW_DATA_AVAILABLE <- file.exists(METH_BETA_RAW_RDS) && file.exists(METH_PHENO_RDS)

## Reads the real preloaded beta matrix + its phenotype table, once per
## running process (cached in .arthomix_cache, same two-tier idea as
## get_raw_eset() above - this file is large enough that re-reading it on
## every "Load preloaded dataset" click would make the Dataset tab feel
## broken). Adds a `sample` column (= pheno$gsm) to the phenotype table so
## methyl_sheet_sample_ids() in qc.R resolves samples by ID rather than by
## row-order coincidence, and orders the matrix's columns to match the
## phenotype table exactly (both are keyed on GSM accession).
load_default_meth_matrix <- function() {
  key <- "meth_default_matrix"
  if (!is.null(.arthomix_cache[[key]])) return(.arthomix_cache[[key]])
  if (!METH_RAW_DATA_AVAILABLE) return(NULL)

  beta  <- readRDS(METH_BETA_RAW_RDS)
  pheno <- readRDS(METH_PHENO_RDS)
  pheno$sample <- as.character(pheno$gsm)

  common <- intersect(colnames(beta), pheno$sample)
  beta  <- beta[, common, drop = FALSE]
  pheno <- pheno[match(common, pheno$sample), , drop = FALSE]

  result <- list(beta = beta, pheno = pheno)
  .arthomix_cache[[key]] <- result
  result
}

## ---------------------------------------------------------------------------
## Cross-Omics: precomputed biomarker-convergence tables
## ---------------------------------------------------------------------------
## Research_Q4_cross_Omics_sexstratified_COPY/cross_Omics_Sexstratified_COPY/,
## a nested project (own .gitignore entry, not tracked by this repo) that
## joins the completed Transcriptomics eQTL-MR/DEG results with the
## completed Methylomics mQTL-MR/DMP results (scripts/01_cross_omics_eqtl_
## mqtl_biomarkers.R) and runs a follow-up cross-omics MR stage (scripts/
## 02_mr_stage_cross_omics.R). Only the already-computed result tables are
## read here - same graceful-degradation contract as METH_DATA_ROOT above
## (returns NULL rather than erroring when CX_DATA_AVAILABLE is FALSE),
## since Cross-Omics is meant to work even on a deployment that hasn't
## copied this folder in.
CX_DATA_ROOT <- normalizePath(file.path("..", "Research_Q4_cross_Omics_sexstratified_COPY", "cross_Omics_Sexstratified_COPY"), mustWork = FALSE)
CX_DATA_AVAILABLE <- dir.exists(CX_DATA_ROOT)
CX_RESULTS_DIR <- file.path(CX_DATA_ROOT, "results")

## Named registry of every table the Dataset tab lets a user browse - label
## -> path, in display order. Kept as one list rather than one
## load_default_cx_<x>() function per table (unlike load_default_dmp()/
## load_default_dmr() above), since every one of these is a flat CSV read
## with no stage/sex argument branching - a single generic loader below
## covers all of them.
CX_TABLE_REGISTRY <- list(
  ## Stage 1 (01_cross_omics_eqtl_mqtl_biomarkers.R) - joins the eQTL-MR
  ## gene panel, mQTL-MR CpG panel, and raw DEG/DMP tables into one
  ## biomarker/annotation table per sex, plus a combined and a Tier-1
  ## (highest-confidence) subset.
  "Master convergence table (all layers, both sexes)" = file.path(CX_RESULTS_DIR, "MASTER_cross_omics_all_layers.csv"),
  "Tier 1 verified biomarkers" = file.path(CX_RESULTS_DIR, "TIER1_MASTER_verified.csv"),
  "eQTL x mQTL combined" = file.path(CX_RESULTS_DIR, "cross_omics_eQTL_mQTL_combined.csv"),
  "eQTL x mQTL - female" = file.path(CX_RESULTS_DIR, "cross_omics_eQTL_mQTL_female.csv"),
  "eQTL x mQTL - male" = file.path(CX_RESULTS_DIR, "cross_omics_eQTL_mQTL_male.csv"),
  ## Stage 2 (02_mr_stage_cross_omics.R) - a follow-up single-instrument
  ## mQTL-MR analysis for the eQTL-MR-significant genes, since the eQTL-MR
  ## gene panel and mQTL-MR CpG panel share zero genes by original design.
  "Cross-omics MR - per-gene summary" = file.path(CX_RESULTS_DIR, "mr_stage_cross_omics_per_gene_summary.csv"),
  "Cross-omics MR - credible genes (full)" = file.path(CX_RESULTS_DIR, "mr_stage_cross_omics_credible_genes_full.csv"),
  "Cross-omics MR - eQTL-significant genes' mQTL instruments" = file.path(CX_RESULTS_DIR, "mr_stage_eqtl_significant_genes_mqtl_mr.csv")
)

## Reads one of the tables above by its display label. Returns NULL rather
## than erroring if CX_DATA_AVAILABLE is FALSE or the specific file isn't
## present, mirroring load_default_dmp()/load_default_dmr() above.
load_default_cx_table <- function(label) {
  path <- CX_TABLE_REGISTRY[[label]]
  if (is.null(path) || !CX_DATA_AVAILABLE || !file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

## ---------------------------------------------------------------------------
## Multi-Omics: paths into Research_05_multiomics_sexstratified - a fully
## executed, independently audited (see its own AUDIT.md) DIABLO/SNF
## multi-omics pipeline over an independent RA anti-TNF cohort (Tao et al.
## 2021), producing real per-cell integration performance, patient
## clustering, gene<->CpG concordance, and pathway enrichment tables. Same
## fail-soft, "returns NULL/nothing to browse rather than erroring" contract
## as CX_DATA_ROOT/CX_DATA_AVAILABLE above, since this module must work even
## on a deployment that hasn't copied this folder in.
## ---------------------------------------------------------------------------
MULTI_DATA_ROOT <- normalizePath(file.path("..", "Research_05_multiomics_sexstratified"), mustWork = FALSE)
MULTI_DATA_AVAILABLE <- dir.exists(MULTI_DATA_ROOT)
MULTI_RESULTS_ROOT <- file.path(MULTI_DATA_ROOT, "results")
MULTI_RESULTS_DIR <- file.path(MULTI_RESULTS_ROOT, "tables")
MULTI_METADATA_DIR <- file.path(MULTI_DATA_ROOT, "metadata")
MULTI_SUMMARY_DIR <- file.path(MULTI_DATA_ROOT, "analyses", "07_cross_analysis_summary", "results", "tables")
MULTI_ADA_DIR <- file.path(MULTI_DATA_ROOT, "analyses", "01_female_male_adalimumab", "results", "tables")
MULTI_ETN_DIR <- file.path(MULTI_DATA_ROOT, "analyses", "02_female_male_etanercept", "results", "tables")

## Named registry of every table the Multi-Omics Dataset tab lets a user
## browse, and every sub-module reads a specific entry from by label (see
## multiomics_helpers.R::multi_read_registry_table()) - label -> path, one
## flat list rather than per-table loader functions, same convention
## CX_TABLE_REGISTRY (above) uses, since every one of these is a flat CSV
## read with no stage/sex argument branching needed at this level (per-cell
## filtering happens downstream, in multiomics_helpers.R).
MULTI_TABLE_REGISTRY <- list(
  "Patient sample matching (all 80 patients)" = file.path(MULTI_METADATA_DIR, "patient_sample_matching_table.csv"),
  "RNA-seq QC summary" = file.path(MULTI_RESULTS_DIR, "RNAseq_QC_summary.csv"),
  "Methylation QC summary" = file.path(MULTI_RESULTS_DIR, "Methylation_QC_summary.csv"),
  "Benchmark vs. published (leakage-safe nested CV, Models A/B/C)" = file.path(MULTI_RESULTS_DIR, "Table8_benchmark_vs_published.csv"),
  "DIABLO performance — drug x sex (response)" = file.path(MULTI_RESULTS_DIR, "Table29_diablo_drugsex_performance.csv"),
  "DIABLO scores — drug x sex (response)" = file.path(MULTI_RESULTS_DIR, "Table29b_diablo_drugsex_scores.csv"),
  "DIABLO panel — drug x sex (response)" = file.path(MULTI_RESULTS_DIR, "Table30_diablo_drugsex_panel.csv"),
  "DIABLO performance — response (drug-pooled)" = file.path(MULTI_RESULTS_DIR, "Table34_diablo_response_sexstratified_performance.csv"),
  "DIABLO scores — response (drug-pooled)" = file.path(MULTI_RESULTS_DIR, "Table34b_diablo_response_sexstratified_scores.csv"),
  "DIABLO panel — response (drug-pooled)" = file.path(MULTI_RESULTS_DIR, "Table35_diablo_response_sexstratified_panel.csv"),
  "SNF performance — drug x sex" = file.path(MULTI_RESULTS_DIR, "Table22_snf_integration_performance.csv"),
  "Random Forest performance — drug x sex" = file.path(MULTI_RESULTS_DIR, "Table37_rf_drugsex_performance.csv"),
  "SNF patient clusters — Adalimumab" = file.path(MULTI_ADA_DIR, "Table_SNFjoint_cluster_assignments_adalimumab.csv"),
  "SNF cluster-response association — Adalimumab" = file.path(MULTI_ADA_DIR, "Table_SNFjoint_cluster_response_association.csv"),
  "SNF concordance NMI — Adalimumab" = file.path(MULTI_ADA_DIR, "Table_SNFjoint_concordanceNMI_adalimumab.csv"),
  "SNF patient clusters — Etanercept" = file.path(MULTI_ETN_DIR, "Table_SNFjoint_cluster_assignments_etanercept.csv"),
  "SNF cluster-response association — Etanercept" = file.path(MULTI_ETN_DIR, "Table_SNFjoint_cluster_response_association.csv"),
  "SNF concordance NMI — Etanercept" = file.path(MULTI_ETN_DIR, "Table_SNFjoint_concordanceNMI_etanercept.csv"),
  "Master six-part summary (integrated vs single-omics)" = file.path(MULTI_SUMMARY_DIR, "Table36_master_six_part_summary.csv"),
  "Candidate multi-omics biomarkers — drug x sex (Etanercept panel)" = file.path(MULTI_SUMMARY_DIR, "Table40_candidate_multiomics_biomarkers_male_female.csv"),
  "Candidate multi-omics biomarkers — response (drug-pooled)" = file.path(MULTI_SUMMARY_DIR, "Table44b_candidate_multiomics_biomarkers_response_male_female.csv"),
  "Gene <-> CpG concordance — drug x sex (Etanercept panel)" = file.path(MULTI_SUMMARY_DIR, "Table42_gene_cpg_concordance_male_female_ETN.csv"),
  "Gene <-> CpG concordance — response (drug-pooled)" = file.path(MULTI_SUMMARY_DIR, "Table45_gene_cpg_concordance_male_female_response.csv"),
  "Pathway enrichment — drug x sex (Etanercept panel)" = file.path(MULTI_SUMMARY_DIR, "Table43_pathway_enrichment_male_female_ETN_panels.csv")
)

## Saved block.splsda fit objects (mixOmics), keyed by MULTI_CELLS' own
## `key` - read via plain readRDS() and list-indexed (e.g. $prop_expl_var)
## rather than requiring library(mixOmics)/mixOmics:: S3 method dispatch,
## since only the fit's own already-computed numeric slots are read, never
## re-fit or re-predicted. Used by multi_diablo_fit() (multiomics_helpers.R)
## to surface the model's real per-component, per-block variance explained -
## an existing mixOmics output this module didn't surface until now.
MULTI_DIABLO_FIT_REGISTRY <- list(
  female_Adalimumab = file.path(MULTI_RESULTS_ROOT, "diablo_drugsex_FEMALE_Adalimumab_fit.rds"),
  male_Adalimumab = file.path(MULTI_RESULTS_ROOT, "diablo_drugsex_MALE_Adalimumab_fit.rds"),
  female_Etanercept = file.path(MULTI_RESULTS_ROOT, "diablo_drugsex_FEMALE_Etanercept_fit.rds"),
  male_Etanercept = file.path(MULTI_RESULTS_ROOT, "diablo_drugsex_MALE_Etanercept_fit.rds"),
  female_response = file.path(MULTI_RESULTS_ROOT, "diablo_response_female_fit.rds"),
  male_response = file.path(MULTI_RESULTS_ROOT, "diablo_response_male_fit.rds")
)

## Whether MOFA2 (and its Python backend, mofapy2, via reticulate/basilisk)
## can actually run in this deployment - gates the Multi-Omics "Live
## Analysis" sub-module's MOFA2 integration step, same
## requireNamespace()-gated-optional-dependency convention
## ARTHOMIX_ASYNC_AVAILABLE/openxlsx already use elsewhere in this app.
## Checked once at app start rather than per-session, since package
## availability doesn't change while the app process is running.
MULTI_MOFA_AVAILABLE <- requireNamespace("MOFA2", quietly = TRUE)

## ---------------------------------------------------------------------------
## Methylomics WGCNA: dedicated cache + optional published-reference loaders
## ---------------------------------------------------------------------------
## get_or_compute_wgcna_blocks() above hardcodes WGCNA_CACHE_DIR, a
## transcriptomics-scoped path (DATA_ROOT/data/cache/wgcna) resolved as a
## free variable, with no cache_dir parameter - reusing it as-is for
## methylomics would write co-methylation network results into the
## transcriptomics cache folder. This is a parallel, methylomics-scoped
## copy (not a modification of the function above), rooted under
## METH_RAW_DATA_ROOT when that folder exists, else falling back next to
## shiny_app/ itself so caching still works on an upload-only deployment.
METH_WGCNA_CACHE_DIR <- file.path(
  if (dir.exists(METH_RAW_DATA_ROOT)) METH_RAW_DATA_ROOT else normalizePath(".", mustWork = FALSE),
  "cache", "wgcna_methyl"
)
dir.create(METH_WGCNA_CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

get_or_compute_meth_wgcna_blocks <- function(key_parts, compute_fn) {
  cache_key <- paste0("methwgcna_", digest::digest(key_parts, algo = "xxhash64"))
  if (!is.null(.arthomix_cache[[cache_key]])) return(.arthomix_cache[[cache_key]])
  disk_path <- file.path(METH_WGCNA_CACHE_DIR, paste0(cache_key, ".rds"))
  result <- if (file.exists(disk_path)) {
    readRDS(disk_path)
  } else {
    result <- compute_fn()
    ## Best-effort: a read-only data folder shouldn't break the
    ## already-computed in-memory result for this session.
    tryCatch(saveRDS(result, disk_path), error = function(e) NULL)
    result
  }
  .arthomix_cache[[cache_key]] <- result
  result
}

## script05_wgcna_sexstratified/ (Research_Q3_METHYLOMICS_sexstratified_COPY)
## - the published, sex-stratified co-methylation network analysis this
## module's Methylomics WGCNA "Compare with published results" panel reads
## from. Table-only, same graceful-degradation contract as load_default_dmp()/
## load_default_dmr() above (returns NULL rather than erroring when
## METH_DATA_AVAILABLE is FALSE or a specific file isn't present) - this is
## reference display only, never the live compute path.
METH_WGCNA_DIR <- file.path(METH_DATA_ROOT, "script05_wgcna_sexstratified", "tables")

## Per-sex module-eigengene-vs-RA-status correlation/p/FDR table
## (module_trait_{sex}[_merged10].csv) from the published run.
load_default_wgcna_module_trait <- function(sex = c("female", "male"), merged = FALSE) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_WGCNA_DIR, sprintf("module_trait_%s%s.csv", sex, if (merged) "_merged10" else ""))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

## Per-sex CpG-to-module assignment table (module_assignment_{sex}[_merged10].csv)
## from the published run.
load_default_wgcna_module_assignment <- function(sex = c("female", "male"), merged = FALSE) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_WGCNA_DIR, sprintf("module_assignment_%s%s.csv", sex, if (merged) "_merged10" else ""))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

## Per-sex DMP/DMR biomarker panel (script04_dmr_sexstratified/tables/
## biomarker_panel_{sex}.csv) - reuses the already-defined METH_DMR_DIR
## above. This is what the published WGCNA script's own enrichment step
## (Fisher's exact test of each significant module's CpGs against this
## panel) tests against; the Functional Enrichment stage below reuses it
## rather than inventing a separate GO/KEGG path.
load_default_dmr_biomarker_panel <- function(sex = c("female", "male")) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_DMR_DIR, sprintf("biomarker_panel_%s.csv", sex))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

## ---------------------------------------------------------------------------
## Methylomics Mendelian Randomization (script08_mendelian_randomization) -
## same table-only, graceful-degradation contract as load_default_dmp()/
## load_default_dmr() above. mod_methyl_mr.R's "Preloaded Data" route reads
## these directly as its computational source of truth (cis-mQTL/GoDMC
## instruments, already clumped r2<0.001/10000kb, already harmonised
## action=2, against the Ishigaki et al. 2022 RA GWAS) rather than
## re-deriving instrument selection/clumping/harmonisation live - see
## METHODS_mendelian_randomization.md for the full pipeline this reproduces.
## script08d_mr_coloc.R (colocalization) is a separate analysis belonging to
## mod_methyl_coloc.R, not read here.
## ---------------------------------------------------------------------------
METH_MR_DIR <- file.path(METH_DATA_ROOT, "script08_mendelian_randomization", "tables")

## Per-sex MR estimates (one row per CpG x method actually run under
## script08b/c's tiered hierarchy: n=1 Wald ratio, n=2 IVW-only, n>=3 the
## full IVW+Egger+weighted-median+mode+simple-mode set).
load_default_mr_estimates <- function(sex = c("female", "male")) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_MR_DIR, sprintf("mr_estimates_%s.csv", sex))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

## The single shared harmonised exposure/outcome dataset (one row per
## harmonised SNP x CpG x outcome triple, both sexes' candidate CpGs
## together) - already clumped and already harmonise_data(action=2)'d, so
## sensitivity analyses / single-SNP estimates / plots can be recomputed
## live from this directly with TwoSampleMR's own functions without
## re-running GoDMC extraction, clumping, or harmonisation.
load_default_mr_harmonised <- function() {
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_MR_DIR, "mr_harmonised_all_cpgs.csv")
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

## Per-sex Steiger directionality flags (SNP, exposure=CpG, steiger_dir,
## steiger_pval), already computed once upstream in 08_mr_instrument_prep.R.
load_default_mr_steiger <- function(sex = c("female", "male")) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_MR_DIR, sprintf("mr_steiger_%s.csv", sex))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

## Per-CpG instrument count after LD clumping (cpg, n_snp) - both sexes'
## candidate CpGs together; the "before clumping -> after clumping" summary
## mod_methyl_mr.R's LD Clumping tab shows for the Preloaded route reads
## this recorded count rather than re-calling ieugwasr::ld_clump() live.
load_default_mr_instrument_counts <- function() {
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_MR_DIR, "instrument_counts.csv")
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

## Precomputed coloc.abf() summary (script08d_mr_coloc.R) - PP.H0-H4 +
## verdict per CpG carried into MR that had >=10 GoDMC candidate SNPs in
## its cis window. mod_methyl_coloc.R's Preloaded route reads this
## directly - no per-SNP GoDMC/RA-GWAS region data is bundled with this
## deployment (see that script's own header comment), so coloc.abf()
## itself cannot be re-run live for these CpGs, only looked up.
load_default_meth_coloc_results <- function() {
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_MR_DIR, "coloc_results.csv")
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

## ---------------------------------------------------------------------------
## Methylomics Diagnostic Classifier (script09_diagnostic_classifier) -
## table-only reference loaders + the two small panel-CpG-only RDS objects
## the module trains/evaluates on live.
## ---------------------------------------------------------------------------
## script09's own panel CpGs come from script07's majority-vote ensemble
## (n_votes >= 2 of LASSO/Boruta/SVM-RFE) - reused here as-is rather than
## re-deriving it, same table-only, graceful-degradation contract as
## load_default_dmp()/load_default_dmr() above.
METH_DIAGNOSTIC_VOTES_DIR <- file.path(METH_DATA_ROOT, "script07_ml_feature_selection", "tables")
METH_DIAGNOSTIC_DIR       <- file.path(METH_DATA_ROOT, "script09_diagnostic_classifier", "tables")

load_default_diagnostic_ensemble_votes <- function(sex = c("female", "male")) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_DIAGNOSTIC_VOTES_DIR, sprintf("ensemble_votes_%s.csv", sex))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

## The published run's own combined-panel and per-probe AUC tables
## (train/internal-test/external-test, with DeLong CIs) - for an optional
## "compare with published results" panel, same idea as mod_methyl_wgcna.R's.
## Never the live compute path.
load_default_diagnostic_panel_auc <- function(sex = c("female", "male")) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_DIAGNOSTIC_DIR, sprintf("diagnostic_panel_auc_%s.csv", sex))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

load_default_diagnostic_perprobe_auc <- function(sex = c("female", "male")) {
  sex <- match.arg(sex)
  if (!METH_DATA_AVAILABLE) return(NULL)
  path <- file.path(METH_DIAGNOSTIC_DIR, sprintf("diagnostic_perprobe_auc_%s.csv", sex))
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

## ---------------------------------------------------------------------------
## Methylomics Diagnostic Classifier: the actual train/test panel data (live
## fit, not just reproduced tables)
## ---------------------------------------------------------------------------
## Unlike METH_BETA_RAW_RDS/METH_PHENO_RDS above (the main pipeline's ~2.1GB
## QC'd matrix, at every retained probe), these two objects are the small
## panel-CpG-only (21 CpGs total) inputs script09 itself trains/evaluates on:
##  - gse42861_internal_panel_celltype_adjusted.rds: GSE42861 reprocessed
##    from raw IDATs through minfi::preprocessNoob() (to be comparable with
##    the external cohort below, since the main pipeline's own normalization
##    is NOT what script09 uses) and granulocyte-adjusted
##    (lm(M ~ Neutro + Eosino) per CpG). All 689 samples, both sexes; script09
##    itself takes a stratified 75/25 train/internal-test split of this per
##    sex-stratum (set.seed(42), caret::createDataPartition(y, p = 0.75)) -
##    reproduced identically in mod_methyl_diagnostic.R's dxm_get_active_data().
##  - gse111942_external_panel.rds: an independent external cohort, Noob-
##    processed the same way, subset to the panel CpGs. All 43 samples are
##    female, so external validation exists only for the female stratum.
## Both live at METH_RAW_DATA_ROOT alongside METH_BETA_RAW_RDS/METH_PHENO_RDS,
## same optional/graceful contract as METH_RAW_DATA_AVAILABLE above - a
## deployment without this folder still runs (Upload-only Diagnostic
## Classifier).
METH_DIAG_INTERNAL_RDS <- file.path(METH_RAW_DATA_ROOT, "gse42861_internal_panel_celltype_adjusted.rds")
METH_DIAG_EXTERNAL_RDS <- file.path(METH_RAW_DATA_ROOT, "gse111942_external_panel.rds")
METH_DIAG_DATA_AVAILABLE <- file.exists(METH_DIAG_INTERNAL_RDS) && file.exists(METH_DIAG_EXTERNAL_RDS)

## Reads both objects once per running process (same two-tier idea as
## load_default_meth_matrix() above, in-memory cache only since these files
## are already small - 131KB/7.6KB - so there's no separate disk-cache tier
## to add). Returns NULL rather than erroring when METH_DIAG_DATA_AVAILABLE
## is FALSE, mirroring every other load_default_*() in this file.
load_default_diagnostic_train_test <- function() {
  key <- "diag_train_test"
  if (!is.null(.arthomix_cache[[key]])) return(.arthomix_cache[[key]])
  if (!METH_DIAG_DATA_AVAILABLE) return(NULL)

  internal <- readRDS(METH_DIAG_INTERNAL_RDS)
  external <- readRDS(METH_DIAG_EXTERNAL_RDS)
  internal$pheno <- as.data.frame(internal$pheno)
  external$pheno <- as.data.frame(external$pheno)

  result <- list(internal = internal, external = external)
  .arthomix_cache[[key]] <- result
  result
}

GEO_SOURCES <- list(
  list(gse = "GSE93272",  role = "Training (whole blood)",   used_in = "Merged into the example cohort"),
  list(gse = "GSE110169", role = "Training (whole blood)",   used_in = "Merged into the example cohort"),
  list(gse = "GSE15573",  role = "Validation (blood, PBMC)", used_in = "Used later, for cross-ancestry validation"),
  list(gse = "GSE89408",  role = "Validation (synovium)",    used_in = "Used later, for cross-tissue validation")
)
geo_link <- function(gse) paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gse)

## Collapses a raw ExpressionSet's probe-level matrix to one row per gene
## symbol (WGCNA::collapseRows, MaxMean rule - the probe with the highest
## mean expression represents its gene) - matching
## scripts/00_shared/03_normalize_batch.R's collapse_to_genes(). Two
## datasets profiled on different microarray platforms almost never share
## probe IDs, so merging them by raw probe ID yields close to nothing in
## common; merging by gene symbol instead is what actually lets the example
## pipeline's own GSE93272 (GPL570) and GSE110169 (GPL13667) combine.
## Reference: Miller JA et al. BMC Bioinformatics 2011;12:322.
collapse_probes_to_genes <- function(eset) {
  fd <- Biobase::fData(eset)
  col <- grep("^gene[ ._]?symbol$", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
  ex <- Biobase::exprs(eset)
  if (is.na(col)) return(ex)
  sym <- as.character(fd[[col]])
  keep <- !is.na(sym) & sym != "" & !grepl("///", sym)
  ex <- ex[keep, , drop = FALSE]; sym <- sym[keep]
  ok <- rowSums(is.na(ex)) < ncol(ex)
  ex <- ex[ok, , drop = FALSE]; sym <- sym[ok]
  if (nrow(ex) == 0) return(ex)
  WGCNA::collapseRows(
    ex, rowGroup = sym, rowID = rownames(ex), method = "MaxMean",
    connectivityBasedCollapsing = FALSE, connectivityPower = 1,
    selectFewestMissing = TRUE, thresholdCombine = NA
  )$datETcollapsed
}

## Disease-group + sex harmonisation from a raw GEO ExpressionSet's
## phenotype data, matching scripts/00_shared/eda.R's harmonize().
eset_harmonize_meta <- function(eset, name) {
  pd <- Biobase::pData(eset)
  clin <- grep(":ch1$", colnames(pd), value = TRUE)
  discol <- grep("disease|status", clin, value = TRUE, ignore.case = TRUE)[1]
  sexcol <- grep("gender|sex", clin, value = TRUE, ignore.case = TRUE)[1]
  dis <- if (!is.na(discol)) as.character(pd[[discol]]) else rep(NA_character_, nrow(pd))
  grp <- ifelse(grepl("sle|lupus", dis, ignore.case = TRUE), "SLE",
          ifelse(grepl("normal|healthy|control", dis, ignore.case = TRUE), "HC",
          ifelse(grepl("rheumatoid|(^|[^a-z])ra([^a-z]|$)", dis, ignore.case = TRUE), "RA", "other")))
  sx <- if (!is.na(sexcol)) toupper(substr(gsub("[^A-Za-z]", "", as.character(pd[[sexcol]])), 1, 1)) else rep(NA_character_, nrow(pd))
  sx[!sx %in% c("F", "M")] <- NA
  data.frame(sample = rownames(pd), dataset = name, group = grp, sex = sx, stringsAsFactors = FALSE)
}

## Loads one individual raw GEO dataset (not merged, not batch-corrected) as
## a plain expr/meta pair, for QC that looks at each source dataset on its
## own instead of the already-merged/batch-corrected working dataset.
## Cached after first read, same as get_raw_eset(), since these matrices are
## large and read-only.
load_individual_dataset <- function(gse_id) {
  key <- paste0("indiv_", gse_id)
  if (!is.null(.arthomix_cache[[key]])) return(.arthomix_cache[[key]])

  result <- if (identical(gse_id, "GSE89408")) {
    ## GSE89408's raw ExpressionSet has no exprs() - counts live in a
    ## separate file; group/sex come from the eset's pData, matching
    ## scripts/goal2_sex_stratified/20_testing_synovium_external.R.
    path <- file.path(RAW_DIR, "GSE89408_counts.txt.gz")
    eset <- get_raw_eset("GSE89408")
    if (!file.exists(path) || is.null(eset)) {
      NULL
    } else {
      cnt <- as.matrix(utils::read.delim(gzfile(path), row.names = 1, check.names = FALSE))
      pref <- gsub("_[0-9]+$", "", colnames(cnt))
      dis_map <- c(normal_tissue = "HC", RA_tissue = "RA", OA_tissue = "other",
                    AG_tissue = "other", undiff_tissue = "other")
      grp <- unname(dis_map[pref])
      sx <- Biobase::pData(eset)[["Sex:ch1"]]
      meta <- data.frame(sample = colnames(cnt), dataset = "GSE89408",
                          group = grp, sex = sx, stringsAsFactors = FALSE)
      keep <- !is.na(meta$group)
      list(expr = cnt[, keep, drop = FALSE], meta = meta[keep, , drop = FALSE],
           label = "GSE89408 (synovium, RNA-seq raw counts)")
    }
  } else {
    eset <- get_raw_eset(gse_id)
    if (is.null(eset) || nrow(Biobase::exprs(eset)) == 0) {
      NULL
    } else {
      expr <- Biobase::exprs(eset)
      meta <- eset_harmonize_meta(eset, gse_id)
      keep <- !is.na(meta$group)
      list(expr = expr[, keep, drop = FALSE], meta = meta[keep, , drop = FALSE],
           label = paste0(gse_id, " (", Biobase::annotation(eset), ", microarray, raw probes)"))
    }
  }
  .arthomix_cache[[key]] <- result
  result
}

## ---------------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

## AI research assistant runs against a local Ollama server instead of the
## Anthropic API - no API key, no per-token billing. Set OLLAMA_BASE_URL to
## point at a non-default Ollama host; otherwise localhost:11434 is assumed.
## Not qwen2.5-coder: verified directly against Ollama's /api/chat that
## qwen2.5-coder returns tool calls as free-text JSON it invents (never
## populates the structured `tool_calls` field), so PubMed lookups would
## silently never fire.
## qwen2.5:1.5b (current default, chosen for speed/size over qwen3:8b) DOES
## populate real structured tool_calls - verified directly - but is
## noticeably less consistent about choosing to call a tool at all: a
## casually-phrased question ("search PubMed for X") got a plain chat reply
## asking for clarification instead of a tool call, while a more directive
## prompt matching this app's actual system prompt ("search before
## answering") correctly triggered pubmed_search. Worth re-testing this way
## against /api/chat before changing the model again - small models vary a
## lot here.
ARTHOMIX_OLLAMA_MODEL <- "qwen2.5:1.5b"
ollama_base_url <- function() Sys.getenv("OLLAMA_BASE_URL", "http://localhost:11434")

## Reachability + model-pulled check in one request: /api/tags lists every
## model Ollama has locally, so a substring match confirms both that the
## server is up and that ARTHOMIX_OLLAMA_MODEL has actually been pulled.
ollama_available <- function() {
  tryCatch({
    con <- url(paste0(ollama_base_url(), "/api/tags"), open = "rb")
    on.exit(close(con), add = TRUE)
    txt <- paste(readLines(con, warn = FALSE), collapse = "")
    grepl(ARTHOMIX_OLLAMA_MODEL, txt, fixed = TRUE)
  }, error = function(e) FALSE)
}

## PubMed lookup tool for the assistant (see mod_assistant.R), so it can back
## scientific claims with real citations instead of inventing them. NCBI's
## E-utilities are free, keyless, and rate-limit-friendly at this volume.
pubmed_search <- function(query, max_results = 5) {
  max_results <- min(max(as.integer(max_results %||% 5), 1L), 10L)
  esearch <- httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi") %>%
    httr2::req_url_query(db = "pubmed", term = query, retmode = "json", retmax = max_results, sort = "relevance") %>%
    httr2::req_url_query(tool = "arthomix-explorer") %>%
    httr2::req_perform() %>%
    httr2::resp_body_json()
  ids <- unlist(esearch$esearchresult$idlist)
  if (length(ids) == 0) return("No PubMed results found for this query.")

  esummary <- httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi") %>%
    httr2::req_url_query(db = "pubmed", id = paste(ids, collapse = ","), retmode = "json") %>%
    httr2::req_perform() %>%
    httr2::resp_body_json()
  recs <- esummary$result

  refs <- lapply(ids, function(pmid) {
    r <- recs[[pmid]]
    if (is.null(r)) return(NULL)
    authors <- vapply(r$authors, function(a) a$name %||% "?", character(1))
    author_str <- if (length(authors) > 1) paste0(authors[1], " et al.") else if (length(authors) == 1) authors[1] else "Unknown author"
    year <- substr(r$pubdate %||% "", 1, 4)
    sprintf("%s (%s). %s. %s. PMID: %s (https://pubmed.ncbi.nlm.nih.gov/%s)",
            author_str, year, r$title, r$fulljournalname %||% r$source %||% "", pmid, pmid)
  })
  refs <- Filter(Negate(is.null), refs)
  if (length(refs) == 0) return("No PubMed results found for this query.")
  paste(refs, collapse = "\n")
}

## ---------------------------------------------------------------------------
## OpenGWAS dataset catalogue search (for ArthOChat + the MR module's
## upload-your-own-GWAS hint)
## ---------------------------------------------------------------------------
## "Which GWAS/eQTL dataset should I use for trait X" is a different question
## from pubmed_search's "what does the literature say", and needs a
## different, keyed API: OpenGWAS (api.opengwas.io, the same catalogue
## scripts/00_shared/10_MR.R itself drew eQTLGen and Okada 2014 from) requires
## a free personal access token, unlike PubMed's keyless E-utilities.
## Degrades to a clear instructional message instead of a raw HTTP error when
## no token is configured, exactly like ollama_available()'s own reachability
## check - so this is always safe to register as a tool even when unusable.

opengwas_token_configured <- function() nzchar(Sys.getenv("OPENGWAS_JWT", ""))

gwas_catalog_search <- function(query, max_results = 10) {
  if (!opengwas_token_configured()) {
    return(paste(
      "OpenGWAS dataset search isn't configured in this deployment - it needs a free",
      "personal access token. Sign in at https://api.opengwas.io/profile/, copy the JWT",
      "shown there, and set it as the OPENGWAS_JWT environment variable before launching",
      "this app (Sys.setenv(OPENGWAS_JWT = \"...\") or in .Renviron), then restart the app.",
      "In the meantime, the Mendelian Randomization tab's \"Upload your own GWAS summary",
      "statistics\" option lets you run MR directly on your own files instead of an",
      "OpenGWAS ID."
    ))
  }
  max_results <- min(max(as.integer(max_results %||% 10), 1L), 25L)

  ## The full catalogue (~50k datasets) is fetched at most once per session
  ## and cached, same two-tier pattern as get_raw_eset() above - repeat
  ## searches (a user or ArthOChat trying several trait terms) filter the
  ## cached table locally instead of re-querying the API every time.
  cat_df <- .arthomix_cache[["opengwas_catalog"]]
  if (is.null(cat_df)) {
    cat_df <- tryCatch(as.data.frame(ieugwasr::gwasinfo()), error = function(e) NULL)
    if (is.null(cat_df) || !nrow(cat_df)) {
      return("Could not reach the OpenGWAS catalogue (network issue, or an invalid/expired token). Try again shortly, or use the Mendelian Randomization tab's upload option instead.")
    }
    .arthomix_cache[["opengwas_catalog"]] <- cat_df
  }

  hit <- rep(FALSE, nrow(cat_df))
  if ("trait" %in% names(cat_df)) hit <- hit | grepl(query, cat_df$trait, ignore.case = TRUE)
  if ("consortium" %in% names(cat_df)) hit <- hit | grepl(query, cat_df$consortium, ignore.case = TRUE)
  if ("author" %in% names(cat_df)) hit <- hit | grepl(query, cat_df$author, ignore.case = TRUE)
  hits <- cat_df[hit, , drop = FALSE]
  if (!nrow(hits)) return(sprintf("No OpenGWAS datasets found matching \"%s\". Try a broader term (e.g. a disease name without qualifiers).", query))
  if ("sample_size" %in% names(hits)) hits <- hits[order(-hits$sample_size), , drop = FALSE]
  hits <- utils::head(hits, max_results)

  col_or_na <- function(df, col) if (col %in% names(df)) as.character(df[[col]]) else rep(NA_character_, nrow(df))
  lines <- sprintf(
    "%s | %s | population: %s | N: %s | year: %s",
    col_or_na(hits, "id"), col_or_na(hits, "trait"), col_or_na(hits, "population"),
    col_or_na(hits, "sample_size"), col_or_na(hits, "year")
  )
  paste(c(sprintf("Top %d OpenGWAS datasets matching \"%s\" (id | trait | population | N | year):", length(lines), query), lines), collapse = "\n")
}

## ---------------------------------------------------------------------------
## Project methodology lookup tool (for ArthOChat)
## ---------------------------------------------------------------------------
## The thesis this app is built on wrote up its own methodology per analysis
## section - a narrative in Chapter_2_subchapter2_sexstratified.md (every
## section, "**2.X Title**" markers) and, for the more involved sections, an
## expanded satellite doc in results/ ending in a curated "References cited"
## bibliography. This is a more authoritative source for "how does THIS
## project do X" than a live PubMed search - it's what the pipeline actually
## did, in the author's own words, plus the exact papers that justify it.
## project_methods() surfaces that directly instead of leaving the assistant
## to reconstruct it from memory.

PROJECT_CHAPTER_MD <- file.path(DATA_ROOT, "Chapter_2_subchapter2_sexstratified.md")

## One row per sub-module: id/title match the TX_MODULES config for that
## module (kept as a static table, not read from TX_MODULES, so this helper
## has no load-order dependency on submodules_registry.R), section is the
## thesis's own numbering, and satellite is the expanded write-up's filename
## in results/ when one exists (NULL for sections only covered in the
## chapter narrative) - both per METHODS_00_INDEX.md, the canonical map.
## `aliases` covers spelling variants (this thesis writes British "-isation";
## sub-module titles elsewhere in the app use American "-ization"), common
## abbreviations, and technique names a user is more likely to type than the
## formal section title.
ARTHOMIX_METHODS_TOPICS <- list(
  list(id = "overview",         title = "Overview and Datasets",       section = "2.1",  satellite = NULL,
       aliases = c("dataset", "datasets", "geo", "cohort")),
  list(id = "preprocessing",    title = "Preprocessing and Batch Correction", section = "2.2", satellite = NULL,
       aliases = c("preprocessing", "pre-processing", "normalisation", "normalization", "batch correction", "batch effect", "combat", "quantile normalisation")),
  list(id = "dge",              title = "Differential Expression",     section = "2.3",  satellite = NULL,
       aliases = c("dge", "differential gene expression", "limma", "deg")),
  list(id = "wgcna",            title = "WGCNA Co-expression Network", section = "2.4",  satellite = "METHODS_2.4_WGCNA_expanded.md",
       aliases = c("wgcna", "co-expression", "coexpression", "gene network", "module detection")),
  list(id = "candidates",       title = "Candidate Gene Identification", section = "2.5", satellite = NULL,
       aliases = c("candidate gene", "candidate genes")),
  list(id = "mr",               title = "Mendelian Randomization",     section = "2.6",  satellite = "METHODS_2.6_mendelian_randomisation.md",
       aliases = c("mr", "mendelian randomisation", "mendelian randomization", "instrumental variable", "causal inference", "eqtl")),
  list(id = "coloc",            title = "Colocalization",              section = "2.7",  satellite = "METHODS_2.7_colocalisation.md",
       aliases = c("coloc", "colocalisation", "colocalization", "coloc.abf", "coloc.susie")),
  list(id = "featureselection", title = "Feature Selection",           section = "2.8",  satellite = "METHODS_2.8_feature_selection.md",
       aliases = c("feature selection", "lasso", "random forest", "svm-rfe", "svm rfe", "biomarker panel")),
  list(id = "diagnostic",       title = "Diagnostic Model",            section = "2.9",  satellite = "METHODS_2.9_diagnostic_model.md",
       aliases = c("diagnostic model", "classifier", "auc", "roc", "cross-validation")),
  list(id = "interaction",      title = "Sex Interaction Analysis",    section = "2.10", satellite = NULL,
       aliases = c("interaction", "sex interaction", "sex-differential", "sex differential")),
  list(id = "crosstissue",      title = "Cross-Tissue Validation",     section = "2.11", satellite = "METHODS_2.11_crosstissue.md",
       aliases = c("cross-tissue", "cross tissue", "synovium", "synovial")),
  list(id = "crossancestry",    title = "Cross-Ancestry Validation",   section = "2.12", satellite = "METHODS_2.12_crossancestry.md",
       aliases = c("cross-ancestry", "cross ancestry", "ancestry")),
  list(id = "enrichment",       title = "Functional Enrichment",       section = "2.13", satellite = "METHODS_2.13_functional_enrichment.md",
       aliases = c("enrichment", "go term", "kegg", "pathway analysis", "gene ontology")),
  list(id = "deconvolution",    title = "Immune Deconvolution",        section = "2.14", satellite = "METHODS_2.14_deconvolution.md",
       aliases = c("deconvolution", "immune cell", "cell composition", "mcp-counter", "mcpcounter")),
  list(id = "nomogram",         title = "Clinical Utility Nomogram",   section = "2.15", satellite = "METHODS_2.15_clinical_utility_nomogram.md",
       aliases = c("nomogram", "clinical utility", "decision curve"))
)

## Text between this section's "**2.X Title**" marker and the next one (or
## EOF). The trailing space in the pattern is load-bearing: matching "2.1"
## must not also match "2.10" - "**2.1 " only matches when the character
## after "2.1" is a space, which is true for section 2.1 but not 2.10.
extract_chapter_section <- function(section) {
  if (!file.exists(PROJECT_CHAPTER_MD)) return(NULL)
  lines <- readLines(PROJECT_CHAPTER_MD, warn = FALSE)
  marker <- paste0("^\\*\\*", gsub("\\.", "\\\\.", section), " ")
  start <- grep(marker, lines)
  if (length(start) == 0) return(NULL)
  start <- start[1]
  all_markers <- grep("^\\*\\*2\\.", lines)
  later <- all_markers[all_markers > start]
  end <- if (length(later) > 0) min(later) - 1 else length(lines)
  paste(lines[start:end], collapse = "\n")
}

## Everything from "## References cited..." to EOF in a satellite doc - by
## convention (confirmed across every satellite doc) the last section, so
## nothing follows it to accidentally include.
extract_references_section <- function(path) {
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path, warn = FALSE)
  start <- grep("^## References cited", lines, ignore.case = TRUE)
  if (length(start) == 0) return(NULL)
  paste(lines[start[1]:length(lines)], collapse = "\n")
}

.truncate_with_note <- function(txt, max_chars) {
  if (nchar(txt) <= max_chars) return(txt)
  paste0(substr(txt, 1, max_chars), "\n[...truncated - ask a more specific follow-up for more detail...]")
}

project_methods <- function(module) {
  q <- tolower(trimws(module %||% ""))
  hit <- Find(function(t) {
    q == tolower(t$id) || q == t$section ||
      grepl(q, tolower(t$title), fixed = TRUE) || grepl(tolower(t$title), q, fixed = TRUE) ||
      any(vapply(t$aliases, function(a) grepl(a, q, fixed = TRUE) || grepl(q, a, fixed = TRUE), logical(1)))
  }, ARTHOMIX_METHODS_TOPICS)

  if (is.null(hit)) {
    return(paste0(
      "No module matched \"", module, "\". Available topics: ",
      paste(vapply(ARTHOMIX_METHODS_TOPICS, function(t) sprintf("%s (Section %s)", t$title, t$section), character(1)), collapse = "; ")
    ))
  }

  narrative <- extract_chapter_section(hit$section)
  refs <- if (!is.null(hit$satellite)) extract_references_section(file.path(DATA_ROOT, "results", hit$satellite)) else NULL

  parts <- c(
    sprintf("## %s (Section %s) - this project's own methodology", hit$title, hit$section),
    if (!is.null(narrative)) .truncate_with_note(narrative, 3000) else "(No chapter write-up found for this section.)",
    "",
    if (!is.null(refs)) .truncate_with_note(refs, 3000) else "(No dedicated reference list for this section - use pubmed_search for supporting literature.)"
  )
  paste(parts, collapse = "\n")
}

## ---------------------------------------------------------------------------
## Methylomics methodology lookup tool (registered on the shared ArthOChat
## drawer - see R/transcriptomics/mod_arthochat.R - so methylomics
## methodology questions stay grounded even though the Methylomics sidebar
## only links to that one shared drawer, rather than running its own
## separate scoped chat session)
## ---------------------------------------------------------------------------
## Mirrors project_methods() above, but against the methylomics pipeline's
## own per-stage write-ups (METH_DATA_ROOT/script0N_*/METHODS_*.md - see
## METH_DATA_ROOT's own comment) instead of the transcriptomics thesis
## chapter. `files` is relative to METH_DATA_ROOT and may be more than one
## path (e.g. DMPs has a plain and an SVA-adjusted write-up); modules with no
## dedicated script (e.g. colocalisation) are simply left out, so
## project_methods_methylomics() falls through to its "no module matched"
## branch and the assistant relies on pubmed_search instead.
METH_METHODS_TOPICS <- list(
  list(id = "qc",             title = "Quality Control",
       files = "script01_dataload_QC/METHODS_load_qc.md",
       aliases = c("qc", "quality control", "detection p-value", "bead count", "sex check", "probe filtering")),
  list(id = "normalization",  title = "Normalization",
       files = "script01_dataload_QC/METHODS_load_qc.md",
       aliases = c("normalization", "normalisation", "dasen", "bmiq", "noob", "funnorm", "swan", "pbc", "quantile normalisation")),
  list(id = "celltype",       title = "Cell-Type Deconvolution",
       files = "script02_celltype/METHODS_celltype.md",
       aliases = c("celltype", "cell type", "cell composition", "houseman", "deconvolution")),
  list(id = "dmp",            title = "Differentially Methylated Positions (DMPs)",
       files = c("script03_dmp_sexstratified/METHODS_dmp_sexstratified.md", "script03_dmp_sva_sexstratified/METHODS_dmp_sva_sexstratified.md"),
       aliases = c("dmp", "differentially methylated position", "differential methylation", "sva", "delta beta")),
  list(id = "dmr",            title = "Differentially Methylated Regions (DMRs)",
       files = "script04_dmr_sexstratified/METHODS_dmr_sexstratified.md",
       aliases = c("dmr", "differentially methylated region", "dmrcate", "comb-p", "bumphunter")),
  list(id = "wgcna",          title = "WGCNA (Co-Methylation Network)",
       files = "script05_wgcna_sexstratified/METHODS_wgcna_sexstratified.md",
       aliases = c("wgcna", "co-methylation", "comethylation", "module eigengene", "soft threshold", "network")),
  list(id = "candidates",     title = "Candidate CpGs (Module-DMR Overlap)",
       files = "script06_module_dmr_venn/METHODS_module_dmr_venn.md",
       aliases = c("candidate cpg", "candidate cpgs", "module-dmr", "venn", "overlap")),
  list(id = "featureselection", title = "ML Feature Selection",
       files = "script07_ml_feature_selection/METHODS_ml_feature_selection.md",
       aliases = c("feature selection", "lasso", "random forest", "svm-rfe", "svm rfe", "elastic net", "boruta")),
  list(id = "mr",             title = "Mendelian Randomization",
       files = "script08_mendelian_randomization/METHODS_mendelian_randomization.md",
       aliases = c("mr", "mendelian randomisation", "mendelian randomization", "mqtl", "instrumental variable", "causal inference")),
  list(id = "diagnostic",     title = "Diagnostic Classifier",
       files = "script09_diagnostic_classifier/METHODS_diagnostic_classifier.md",
       aliases = c("diagnostic model", "diagnostic classifier", "classifier", "auc", "roc", "panel"))
)

project_methods_methylomics <- function(module) {
  q <- tolower(trimws(module %||% ""))
  hit <- Find(function(t) {
    q == tolower(t$id) ||
      grepl(q, tolower(t$title), fixed = TRUE) || grepl(tolower(t$title), q, fixed = TRUE) ||
      any(vapply(t$aliases, function(a) grepl(a, q, fixed = TRUE) || grepl(q, a, fixed = TRUE), logical(1)))
  }, METH_METHODS_TOPICS)

  if (is.null(hit)) {
    return(paste0(
      "No Methylomics module matched \"", module, "\". Available topics: ",
      paste(vapply(METH_METHODS_TOPICS, function(t) t$title, character(1)), collapse = "; ")
    ))
  }

  paths <- file.path(METH_DATA_ROOT, hit$files)
  bodies <- Filter(Negate(is.null), lapply(paths, function(p) {
    if (!file.exists(p)) return(NULL)
    paste(readLines(p, warn = FALSE), collapse = "\n")
  }))

  header <- sprintf("## %s - this project's own methodology", hit$title)
  if (length(bodies) == 0) {
    return(paste(header, "(No write-up available for this stage in this deployment.)", sep = "\n"))
  }
  paste(header, "", .truncate_with_note(paste(bodies, collapse = "\n\n---\n\n"), 4000), sep = "\n")
}

read_table_safe <- function(filename, dir = TABLES_DIR) {
  path <- file.path(dir, filename)
  if (!file.exists(path)) return(NULL)
  as.data.frame(data.table::fread(path, showProgress = FALSE))
}

figure_exists <- function(filename) {
  file.exists(file.path(FIGURES_DIR, filename))
}

## ---------------------------------------------------------------------------
## GWAS/eQTL summary-statistics upload helpers - shared by mod_mr.R's
## upload-your-own-exposure/outcome mode, mod_coloc.R's upload-your-own-GWAS
## mode, and mod_crossancestry.R's upload-a-replacement-replication/transfer-
## GWAS mode - one column-mapping UI and one column-guessing table instead of
## three independent copies that could drift apart on which header spellings
## each recognises.
## ---------------------------------------------------------------------------

GWAS_COL_PATTERNS <- list(
  snp  = c("^snp$", "^rsid$", "^rs_?id$", "variant"),
  ea   = c("^effect_allele$", "^ea$", "^a1$", "^alt$"),
  oa   = c("^other_allele$", "^oa$", "^a2$", "^ref$", "^nea$"),
  beta = c("^beta$", "^b$", "^effect$"),
  se   = c("^se$", "^standard_error$", "^stderr$"),
  pval = c("^p$", "^pval$", "^p_value$", "^pvalue$"),
  eaf  = c("^eaf$", "^freq$", "^maf$", "^effect_allele_freq"),
  n    = c("^n$", "^samplesize$", "^sample_size$", "^total_n$")
)

guess_gwas_col <- function(cols, patterns) {
  for (p in patterns) { hit <- grep(p, cols, ignore.case = TRUE, value = TRUE); if (length(hit)) return(hit[1]) }
  NA_character_
}

read_uploaded_table <- function(path) tryCatch(as.data.frame(data.table::fread(path, showProgress = FALSE)), error = function(e) NULL)

## Each uploaded file should be parsed from disk exactly ONCE per caller
## (wrap in a `reactive()` keyed on the fileInput, as mod_mr.R does with
## exp_df_r/out_df_r) - this function itself just renders the mapping UI from
## whatever data.frame it's handed.
##
## `ns` is the caller's own `session$ns` (the generated input ids must live in
## whichever module renders this). `extra_fields` appends additional required
## column pickers after the core seven every caller needs - e.g. mod_coloc.R's
## coloc.abf() also needs a sample-size column mod_mr.R's MR estimator never
## did (`extra_fields = "n"`).
gwas_col_map_ui <- function(ns, file_input, df_reactive, prefix, label, extra_fields = character(0)) {
  extra_labels <- c(n = "Sample size (N)")
  renderUI({
    fi <- file_input()
    req(fi)
    df <- df_reactive()
    validate(need(!is.null(df) && ncol(df) >= 4, sprintf("Could not read the %s file as a delimited table with at least 4 columns.", label)))
    cols <- colnames(df)
    guess <- function(key) { g <- guess_gwas_col(cols, GWAS_COL_PATTERNS[[key]]); if (is.na(g)) cols[1] else g }
    tagList(
      div(class = "empty-note", style = "font-size:12.5px; margin: 4px 0; padding: 8px 12px;",
        icon("circle-check"), strong(sprintf(" Uploaded: %s", fi$name)),
        sprintf(" - %s rows, %d columns. Columns auto-matched below - check them.", format(nrow(df), big.mark = ","), ncol(df))),
      fluidRow(
        column(6, selectInput(ns(paste0(prefix, "_snp")), "SNP ID", cols, selected = guess("snp"))),
        column(6, selectInput(ns(paste0(prefix, "_beta")), "Beta (effect size)", cols, selected = guess("beta")))
      ),
      fluidRow(
        column(6, selectInput(ns(paste0(prefix, "_se")), "Standard error", cols, selected = guess("se"))),
        column(6, selectInput(ns(paste0(prefix, "_pval")), "p-value", cols, selected = guess("pval")))
      ),
      fluidRow(
        column(6, selectInput(ns(paste0(prefix, "_ea")), "Effect allele", cols, selected = guess("ea"))),
        column(6, selectInput(ns(paste0(prefix, "_oa")), "Other allele", cols, selected = guess("oa")))
      ),
      selectInput(ns(paste0(prefix, "_eaf")), "Effect allele frequency (optional - improves palindromic SNP resolution)",
                  c("(none)" = "", cols), selected = { g <- guess_gwas_col(cols, GWAS_COL_PATTERNS$eaf); if (is.na(g)) "" else g }),
      lapply(extra_fields, function(key) selectInput(ns(paste0(prefix, "_", key)), extra_labels[[key]] %||% key, cols, selected = guess(key)))
    )
  })
}

## ---------------------------------------------------------------------------
## Shared MR estimator hierarchy - moved here (unchanged) from mod_mr.R so
## mod_crossancestry.R's live-upload arm can fit the SAME estimator hierarchy
## (Section 2.6.5: >=3 SNPs -> IVW + MR-Egger + weighted median, 2 SNPs -> IVW
## alone, 1 SNP -> Wald ratio) instead of a second, drift-prone copy. Pure
## function of its arguments - no reactive/session dependency - so this is a
## straight relocation, not a rewrite.
## ---------------------------------------------------------------------------
estimate_mr_set <- function(d, include_mode = FALSE, full = TRUE, ci_level = 0.95, run_presso = FALSE) {
  n_snp <- nrow(d)
  methods <- list()
  heterogeneity <- NULL
  pleiotropy <- NULL
  presso <- NULL
  alpha <- 1 - ci_level
  z <- stats::qnorm(1 - alpha / 2)

  if (n_snp == 1) {
    b <- d$beta.outcome[1] / d$beta.exposure[1]
    se <- abs(d$se.outcome[1] / d$beta.exposure[1])
    methods[["Wald ratio"]] <- c(estimate = b, se = se, ci_low = b - z * se, ci_high = b + z * se,
                                  p = 2 * stats::pnorm(-abs(b / se)))
    primary_method <- "Wald ratio"
  } else {
    mrobj <- MendelianRandomization::mr_input(
      bx = d$beta.exposure, bxse = d$se.exposure,
      by = d$beta.outcome, byse = d$se.outcome,
      snps = d$SNP, exposure = unique(d$gene)[1], outcome = "RA risk (Okada 2014 GWAS)"
    )
    ## model = "random": MendelianRandomization::mr_ivw()'s own default
    ## ("fixed") understates the SE relative to the bundled dataset's own
    ## cached values. The upstream pipeline's IVW is a random-effects
    ## (multiplicative dispersion, phi = max(1, Q/(n-1))) estimator by
    ## default - confirmed by reproducing the bundled cached b/se/p
    ## bit-for-bit across multiple genes (e.g. LPCAT2: se 0.0248 here vs
    ## 0.0220 under "fixed"; b is identical either way, only the SE/p
    ## differ). Without this, every p-value and CI shown here would be
    ## too narrow relative to the bundled reference numbers.
    ivw <- MendelianRandomization::mr_ivw(mrobj, model = "random", alpha = alpha)
    methods[["IVW"]] <- c(estimate = ivw@Estimate, se = ivw@StdError, ci_low = ivw@CILower, ci_high = ivw@CIUpper, p = ivw@Pvalue)
    primary_method <- "IVW"
    if (n_snp >= 3 && full) {
      ## Bootstrap SE (10,000 draws, package default seed) - the bundled
      ## dataset's own weighted-median SE is also a bootstrap draw, made
      ## inside a sequential full-gene-set run, so it cannot be
      ## reproduced bit-for-bit gene-by-gene here. Never the primary
      ## estimate in the default hierarchy, so this is a disclosed
      ## limitation, not a correctness issue for the number that
      ## actually drives conclusions.
      med <- MendelianRandomization::mr_median(mrobj, alpha = alpha)
      methods[["Weighted median"]] <- c(estimate = med@Estimate, se = med@StdError, ci_low = med@CILower, ci_high = med@CIUpper, p = med@Pvalue)
      egg <- MendelianRandomization::mr_egger(mrobj, alpha = alpha)
      ## Egger's Estimate/StdError.Est/Intercept/StdError.Int all match
      ## the bundled cached numbers exactly, but its own Pvalue.Est/
      ## Pvalue.Int don't: TwoSampleMR::mr_egger_regression gets
      ## p-values from base R's summary(lm(...)) - a t-distribution on
      ## n_snp-2 residual df - and neither of this package's own
      ## distribution options ("normal", "t-dist") reproduces that df.
      ## Recomputed directly here instead of trusting either.
      egger_df <- n_snp - 2
      egger_slope_p <- 2 * stats::pt(-abs(egg@Estimate / egg@StdError.Est), df = egger_df)
      egger_int_p <- 2 * stats::pt(-abs(egg@Intercept / egg@StdError.Int), df = egger_df)
      methods[["MR-Egger"]] <- c(estimate = egg@Estimate, se = egg@StdError.Est, ci_low = egg@CILower.Est, ci_high = egg@CIUpper.Est, p = egger_slope_p)
      heterogeneity <- list(Q = unname(ivw@Heter.Stat[1]), Q_df = n_snp - 1, Q_pval = unname(ivw@Heter.Stat[2]))
      pleiotropy <- list(intercept = egg@Intercept, se = egg@StdError.Int, p = egger_int_p, I2 = egg@I.sq)
      if (isTRUE(include_mode)) {
        mode_est <- tryCatch(MendelianRandomization::mr_mbe(mrobj, alpha = alpha), error = function(e) NULL)
        if (!is.null(mode_est)) {
          methods[["Weighted mode"]] <- c(estimate = mode_est@Estimate, se = mode_est@StdError,
                                           ci_low = mode_est@CILower, ci_high = mode_est@CIUpper, p = mode_est@Pvalue)
        }
      }
    }
    ## MR-PRESSO (Verbanck et al. 2018): simulation-based global test for
    ## horizontal pleiotropy plus, when the global test is significant,
    ## an outlier-corrected re-estimate excluding the flagged SNP(s).
    ## Needs >=4 instruments by design (one more than MR-Egger, since it
    ## must leave >=3 after excluding a candidate outlier); optional and
    ## off by default since NbDistribution=1000 simulations add a
    ## noticeable delay per run.
    if (n_snp >= 4 && full && isTRUE(run_presso)) {
      presso <- tryCatch({
        pr <- MRPRESSO::mr_presso(
          BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome", SdExposure = "se.exposure", data = d,
          OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
          NbDistribution = 1000, SignifThreshold = alpha, seed = 2024
        )
        mrr <- pr$`Main MR results`
        raw <- mrr[mrr$`MR Analysis` == "Raw", ]
        corr <- mrr[mrr$`MR Analysis` == "Outlier-corrected", ]
        list(
          global_p = pr$`MR-PRESSO results`$`Global Test`$Pvalue,
          raw_estimate = raw$`Causal Estimate`[1], raw_p = raw$`P-value`[1],
          corrected_estimate = if (!is.na(corr$`Causal Estimate`[1])) corr$`Causal Estimate`[1] else NA_real_,
          corrected_p = if (!is.na(corr$`P-value`[1])) corr$`P-value`[1] else NA_real_,
          n_outliers = length(pr$`MR-PRESSO results`$`Outlier Test`$Pvalue[
            !is.na(suppressWarnings(as.numeric(pr$`MR-PRESSO results`$`Outlier Test`$Pvalue))) &
            suppressWarnings(as.numeric(pr$`MR-PRESSO results`$`Outlier Test`$Pvalue)) < alpha
          ])
        )
      }, error = function(e) list(error = conditionMessage(e)))
    }
  }

  res_table <- as.data.frame(do.call(rbind, methods))
  res_table <- data.frame(method = names(methods), res_table, row.names = NULL, check.names = FALSE)
  res_table$primary <- res_table$method == primary_method

  list(res_table = res_table, n_snp = n_snp, heterogeneity = heterogeneity, pleiotropy = pleiotropy,
       primary_method = primary_method, presso = presso, ci_level = ci_level)
}

## ---------------------------------------------------------------------------
## Shared cis-eQTL instrument table (Okada 2014 RA GWAS harmonisation) -
## MR_primary_objects.rds, loaded once per caller. Both mod_mr.R (its
## bundled-dataset gene picker and batch screen) and mod_crossancestry.R (its
## live-upload replication/transfer arm, which re-harmonises these SAME
## exposure-side rows against a *different* outcome GWAS) need this table, so
## it's loaded and relabel-fixed here once instead of twice.
##
## `$dat$gene` is stale for 21 SNPs / 37 genes that share an eQTL with a
## neighbouring gene - the pair-matched relabelling fix was written into the
## upstream pipeline but this cached object predates it. The correct label is
## recomputed here from `$inst` (SNP + exposure-dataset ID pair) rather than
## trusted from the cached column, and exactly what changed is returned in
## `relabel_check` for callers to disclose rather than fix silently.
## ---------------------------------------------------------------------------
load_mr_instrument_table <- function() {
  mr_obj <- readRDS(MR_PRIMARY_OBJECTS_RDS)

  dat <- as.data.frame(mr_obj$dat)
  dat$gene <- mr_obj$inst$gene[match(
    paste(dat$SNP, dat$id.exposure),
    paste(mr_obj$inst$SNP, mr_obj$inst$id.exposure)
  )]
  dat <- dat[dat$mr_keep, , drop = FALSE]

  inst_dt <- as.data.table(mr_obj$inst)
  dat_dt  <- as.data.table(mr_obj$dat)
  chk <- merge(dat_dt[, .(SNP, id.exposure, gene_cached = gene)],
               inst_dt[, .(SNP, id.exposure, gene_correct = gene)],
               by = c("SNP", "id.exposure"))
  mism <- chk[gene_cached != gene_correct]
  relabel_check <- list(n_snp = nrow(mism), genes = unique(c(mism$gene_cached, mism$gene_correct)))

  list(dat = dat, primary = as.data.frame(mr_obj$primary), relabel_check = relabel_check)
}

## ---------------------------------------------------------------------------
## Shared chart styling
## ---------------------------------------------------------------------------
## One fixed, colorblind-safe-ordered palette used everywhere instead of
## ggplot's default hue wheel, so live figures across submodules read as one
## system rather than each picking its own colors.

ARTHOMIX_COLORS <- list(
  blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", yellow = "#eda100",
  magenta = "#e87ba4", violet = "#4a3aa7", red = "#e34948",
  ink = "#0b0b0b", ink_secondary = "#52514e", ink_muted = "#898781",
  grid = "#e1e0d9", axis = "#c3c2b7"
)

## Reserved for state, never reused for an ordinary series.
ARTHOMIX_STATUS <- list(good = "#0ca30c", warning = "#fab219", critical = "#d03b3b")

## Fixed hue order for categorical series (group, sex, ...); a factor with
## more levels than this just runs out of distinct colors rather than
## cycling ggplot's rainbow.
arthomix_pair <- function(levels) {
  pal <- c(ARTHOMIX_COLORS$blue, ARTHOMIX_COLORS$red, ARTHOMIX_COLORS$aqua,
            ARTHOMIX_COLORS$orange, ARTHOMIX_COLORS$violet, ARTHOMIX_COLORS$magenta,
            ARTHOMIX_COLORS$yellow)
  levels <- as.character(sort(unique(levels)))
  setNames(pal[seq_along(levels)], levels)
}

theme_arthomix <- function(base_size = 12) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = ARTHOMIX_COLORS$grid, linewidth = 0.3),
      axis.line = element_line(color = ARTHOMIX_COLORS$axis, linewidth = 0.3),
      axis.ticks = element_line(color = ARTHOMIX_COLORS$axis, linewidth = 0.3),
      axis.text = element_text(color = ARTHOMIX_COLORS$ink_secondary),
      axis.title = element_text(color = ARTHOMIX_COLORS$ink_secondary),
      plot.title = element_text(color = ARTHOMIX_COLORS$ink, size = base_size, hjust = 0),
      plot.subtitle = element_text(color = ARTHOMIX_COLORS$ink_muted, size = base_size - 2, hjust = 0),
      legend.position = "bottom",
      legend.title = element_text(color = ARTHOMIX_COLORS$ink_secondary, size = base_size - 2),
      legend.text = element_text(color = ARTHOMIX_COLORS$ink_secondary, size = base_size - 2),
      strip.background = element_blank(),
      strip.text = element_text(color = ARTHOMIX_COLORS$ink, size = base_size - 1)
    )
}

## ---------------------------------------------------------------------------
## Sample-level QC, shared by the Dataset tab and the Preprocessing submodule
## ---------------------------------------------------------------------------
## Three standard array/RNA-seq QC signals per sample: total signal, number
## of detected features, and mean correlation to every other sample (a
## sample that doesn't correlate with the rest of the cohort is the
## clearest sign of a technical outlier). Flags are relative to the
## dataset's own median +/- mad_k median absolute deviations, not a fixed
## cutoff, so the same function works whether `expr` holds raw counts,
## log-CPM or log2 microarray intensities.

compute_sample_qc <- function(expr, mad_k = 3, top_n_cor = 2000) {
  detect_cutoff <- stats::quantile(expr, 0.25, na.rm = TRUE)
  signal   <- colSums(expr, na.rm = TRUE)
  detected <- colSums(expr > detect_cutoff, na.rm = TRUE)

  gene_var <- apply(expr, 1, stats::var, na.rm = TRUE)
  top_idx  <- order(gene_var, decreasing = TRUE)[seq_len(min(top_n_cor, nrow(expr)))]
  sample_cor <- cor(expr[top_idx, , drop = FALSE])
  mean_cor <- (colSums(sample_cor) - 1) / (ncol(sample_cor) - 1)

  is_outlier <- function(x, low_only = FALSE) {
    m <- stats::median(x); s <- stats::mad(x)
    if (s == 0) return(rep(FALSE, length(x)))
    if (low_only) (m - x) > mad_k * s else abs(x - m) > mad_k * s
  }

  data.frame(
    sample = colnames(expr),
    signal = signal,
    detected = detected,
    mean_cor = mean_cor,
    flag_signal = is_outlier(signal),
    flag_detected = is_outlier(detected),
    flag_cor = is_outlier(mean_cor, low_only = TRUE),
    row.names = NULL
  )
}

## ---------------------------------------------------------------------------
## THE NORMALISATION IMPLEMENTATION - audited 2026-08-12 (see Data
## Exploration & Normalisation tab in mod_preprocessing.R for the user-
## facing side of this audit). Summary, so this stays easy to find later:
##
##   Where raw data is loaded:    mod_dataset.R (Dataset tab: preloaded /
##                                 GEO fetch / upload) and mod_preprocessing.R
##                                 ("Preprocessing" tab, per-source). Either
##                                 way, the result is a genes(or probes) x
##                                 samples numeric matrix (`expr`) kept
##                                 strictly separate from the sample-level
##                                 metadata data.frame (`meta`) - phenotype/
##                                 clinical columns are never columns of
##                                 `expr`, so they can't accidentally be
##                                 normalised as if they were expression
##                                 values.
##   Where preprocessing starts:  mod_preprocessing.R's per-source filters
##                                 (missing-value tolerance, expression
##                                 pattern exclusion, log2 decision) run
##                                 BEFORE datasets are merged; a second,
##                                 merged-matrix filter (expression/variance
##                                 percentile) runs in the "Batch correction"
##                                 tab, immediately before normalisation.
##                                 Order is explicit and enforced by data
##                                 flow, not just convention: Filter -> (log2
##                                 Transform, applied per-source pre-merge) ->
##                                 Normalise -> Batch-correct.
##   Normalisation function:      this file, below - needs_quantile_norm()
##                                 decides IF; the actual transform is either
##                                 `limma::normalizeBetweenArrays(..., method
##                                 = "quantile")` (microarray/log-scale data)
##                                 or `edgeR::calcNormFactors(..., method =
##                                 "TMM")` + `edgeR::cpm(..., log = TRUE)`
##                                 (raw RNA-seq counts) - see the `norm_method`
##                                 branch in mod_preprocessing.R's `result`
##                                 reactive (Batch correction tab).
##   Data it receives/returns:    a numeric expression matrix in, the same
##                                 shape out (never touches `meta`).
##   Correctness verdict:         both methods are appropriate for their
##                                 respective data type and were NOT
##                                 rewritten - quantile normalisation for
##                                 microarray/already-log-scale data, TMM for
##                                 raw counts, in the right order relative to
##                                 filtering. The auto-detect heuristic
##                                 (needs_quantile_norm(), below) already
##                                 guards against silently re-normalising
##                                 data that looks already normalised.
##   Gap found and addressed:     none of this was surfaced to the user
##                                 BEFORE they picked a method - no raw-data
##                                 preview, no visual "does this look already
##                                 normalised" check, no missing/zero/
##                                 infinite/duplicate-feature audit, and no
##                                 record of exactly what was applied for
##                                 reproducibility. That gap is what the new
##                                 Data Exploration & Normalisation tab (and
##                                 the small helpers just below) fills - it
##                                 calls this same summarize_norm_diagnostics()/
##                                 needs_quantile_norm() and the same
##                                 limma::normalizeBetweenArrays()/edgeR TMM
##                                 calls Batch Correction already uses, rather
##                                 than reimplementing normalisation a second
##                                 time.
## ---------------------------------------------------------------------------

## Per-sample scale diagnostics, before/after normalisation. Same metrics
## and same "needs quantile normalisation" decision rule the thesis
## pipeline uses in scripts/00_shared/03_normalize_batch.R: normalise if
## the data are still on a linear scale, or per-sample medians/IQRs still
## disagree by more than 0.5 on the log scale.
summarize_norm_diagnostics <- function(m) {
  sample_medians <- apply(m, 2, stats::median, na.rm = TRUE)
  sample_iqr <- apply(m, 2, stats::IQR, na.rm = TRUE)
  data.frame(
    n_samples = ncol(m), n_genes = nrow(m),
    max_value = max(m, na.rm = TRUE), min_value = min(m, na.rm = TRUE),
    median_sd = stats::sd(sample_medians, na.rm = TRUE), iqr_sd = stats::sd(sample_iqr, na.rm = TRUE),
    median_range = diff(range(sample_medians, na.rm = TRUE)),
    iqr_range = diff(range(sample_iqr, na.rm = TRUE))
  )
}

needs_quantile_norm <- function(diag) {
  diag$max_value > 100 || diag$median_sd > 0.5 || diag$iqr_sd > 0.5
}

## ---------------------------------------------------------------------------
## Data Exploration & Normalisation tab helpers (mod_preprocessing.R) - raw-
## data audit and data-type detection only. None of these apply a
## normalisation themselves beyond what's delegated straight to the same
## limma/edgeR calls Batch Correction already uses (see apply_chosen_norm()
## below) - this section exists to decide WHETHER and WHICH, and to surface
## raw-data quality issues, not to reimplement normalisation math.
## ---------------------------------------------------------------------------

## Raw-data health counts a scientist would want before deciding anything:
## missingness, zeros (routine and expected in RNA-seq counts, worth a
## second look in microarray/log data), non-finite values (never expected -
## a sign of a broken upload, e.g. a ratio/fold-change column mistaken for
## an expression matrix) and duplicated feature IDs (silently ambiguous for
## any downstream row-name-keyed merge/lookup). None of this was checked
## anywhere in the app before this tab.
expr_raw_health <- function(expr) {
  m <- as.matrix(expr)
  vals <- as.numeric(m)
  list(
    n_features = nrow(m), n_samples = ncol(m),
    n_missing = sum(is.na(vals)),
    n_zero = sum(vals == 0, na.rm = TRUE),
    n_infinite = sum(is.infinite(vals)),
    n_duplicated_features = sum(duplicated(rownames(m)))
  )
}

## Best-effort guess at what kind of expression data this is, from the
## values alone (no filename/platform metadata to go on) - drives which
## normalisation methods the Data Exploration & Normalisation tab offers.
## Mirrors the distinction Batch Correction's own TMM-vs-quantile radio
## already draws by hand, as a reusable, named function:
##   "counts"              - non-negative, (near-)integer-valued, large
##                            range - raw RNA-seq counts. TMM + log2-CPM is
##                            appropriate; quantile-normalising counts
##                            directly is not (it assumes a roughly
##                            continuous, already-comparable distribution
##                            shape, which raw count data doesn't have).
##   "already_normalised"  - negative values present (only possible after a
##                            log/z-score-like transform), or per-sample
##                            medians/IQRs already agree
##                            (needs_quantile_norm() is FALSE) - most likely
##                            already log-transformed and/or normalised
##                            upstream, e.g. a public series' already-
##                            processed matrix.
##   "expression"           - anything else: linear- or log-scale
##                            expression that hasn't been shown to already
##                            agree across samples - the normal case
##                            quantile normalisation targets.
detect_expr_data_type <- function(expr) {
  m <- as.matrix(expr)
  finite_vals <- m[is.finite(m)]
  if (length(finite_vals) == 0) return("expression")
  has_negative <- any(finite_vals < 0)
  ## Threshold is 0.6, not ~1: real "raw counts" files are routinely
  ## expected-count estimates from an EM-based quantifier (RSEM, Salmon,
  ## kallisto), not literal integer read tallies - GSE89408's own bundled
  ## counts file (scripts/00_shared/eda.R's source) is ~89% near-integer,
  ## not 100%. Verified against a genuinely continuous synthetic matrix
  ## (~0% near-integer) to confirm 0.6 keeps real microarray/log-scale data
  ## out of this bucket while still catching expected-count RNA-seq data.
  looks_like_counts <- !has_negative &&
    mean(abs(finite_vals - round(finite_vals)) < 1e-6) > 0.6 &&
    max(finite_vals) > 100
  if (looks_like_counts) return("counts")
  diag <- summarize_norm_diagnostics(m)
  if (has_negative || !needs_quantile_norm(diag)) return("already_normalised")
  "expression"
}

## Filter -> Transform, in that explicit order, kept as one small function
## so the Data Exploration & Normalisation tab can't accidentally apply them
## the other way round. Filtering first means a normalisation method that
## estimates parameters from the data (quantile ranks, TMM size factors)
## never sees near-empty features that would only distort those estimates;
## log2 is applied last, evaluating the filter thresholds on the scale a
## user actually set them on (e.g. "raw count >= 10"), not on already-
## logged values. Missing values remaining after the row-level filter are
## median-imputed, the same simple technique mod_pp_source_server() (the
## per-dataset Preprocessing step) already uses - not a new one invented
## here.
filter_and_transform_expr <- function(expr, min_expr = 0, min_sample_frac = 0,
                                        max_na_pct = 100, log2_transform = FALSE) {
  m <- as.matrix(expr)
  ## Infinite values are never a valid expression measurement; treat them as
  ## missing so they fall through the same missing-value handling below
  ## instead of silently winning every ">= min_expr" comparison, or
  ## corrupting max()/quantile()-based diagnostics downstream. Nothing
  ## upstream of this tab currently guards against this (see
  ## expr_raw_health()'s n_infinite count, surfaced as its own warning) -
  ## this is the one place that actually neutralises it before it reaches
  ## filtering, transformation or normalisation.
  m[is.infinite(m)] <- NA
  keep_na <- (rowMeans(is.na(m)) * 100) <= max_na_pct
  above <- m >= min_expr
  above[is.na(above)] <- FALSE
  keep_expr <- (rowMeans(above) * 100) >= min_sample_frac
  m <- m[keep_na & keep_expr, , drop = FALSE]
  if (anyNA(m)) {
    row_med <- apply(m, 1, stats::median, na.rm = TRUE)
    na_idx <- which(is.na(m), arr.ind = TRUE)
    m[na_idx] <- row_med[na_idx[, 1]]
  }
  if (isTRUE(log2_transform)) {
    m[m <= 0] <- NA
    m <- log2(m)
    m <- m[stats::complete.cases(m), , drop = FALSE]
  }
  m
}

## Normalise (or deliberately don't) - the SAME limma::normalizeBetweenArrays
## / edgeR TMM calls Batch Correction's `result` reactive already uses, not
## a reimplementation. method: "none" (used as filtered/transformed, no
## further scaling - what "Do not normalise" in the tab maps to), "quantile"
## (microarray/log-scale expression) or "tmm" (raw RNA-seq counts).
apply_chosen_norm <- function(expr, method = c("none", "quantile", "tmm")) {
  method <- match.arg(method)
  m <- as.matrix(expr)
  if (identical(method, "none")) {
    return(list(expr = m, label = "None - used as filtered/transformed, no normalisation applied"))
  }
  if (identical(method, "tmm")) {
    validate(need(all(m >= 0, na.rm = TRUE),
                  "TMM normalisation expects raw, non-negative counts, but this data has negative values."))
    counts <- round(m)
    storage.mode(counts) <- "integer"
    dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts), method = "TMM")
    return(list(expr = edgeR::cpm(dge, log = TRUE, prior.count = 1),
                label = "TMM (edgeR::calcNormFactors) plus log2-CPM"))
  }
  mtx <- limma::normalizeBetweenArrays(m, method = "quantile")
  rownames(mtx) <- rownames(m); colnames(mtx) <- colnames(m)
  list(expr = mtx, label = "Quantile normalisation (limma::normalizeBetweenArrays)")
}

## Plain-language read of a before/after normalisation comparison - QC
## evidence, not proof (a histogram looking better is not the same as
## "correct"). Five fixed messages, matching what's actually measurable
## from summarize_norm_diagnostics()/needs_quantile_norm() rather than
## claiming more than the data shows.
normalisation_assessment <- function(diag_before, diag_after, method) {
  if (identical(method, "none")) return("No substantial transformation was applied.")
  if (is.null(diag_after) || inherits(diag_after, "error")) {
    return("Normalisation could not be performed because of invalid input values.")
  }
  before_ok <- !needs_quantile_norm(diag_before)
  after_ok  <- !needs_quantile_norm(diag_after)
  if (before_ok) {
    "Data appear already normalised - this step had little additional effect."
  } else if (after_ok) {
    "Distributions are more comparable across samples."
  } else {
    "Large differences between sample distributions remain."
  }
}

## Headline data-processing status for the Data Exploration & Normalisation
## tab's status panel - a coarser, user-facing classification built ONLY on
## top of detect_expr_data_type()/needs_quantile_norm() (never re-derives
## its own detection rule, so this can't disagree with the value shown
## elsewhere in that tab). Exactly three states, matched to what's actually
## knowable from the matrix alone: "normalized" (already normalised-
## looking - detect_expr_data_type()'s "already_normalised" bucket),
## "not_normalized" (raw counts or expression that hasn't been shown to
## already agree across samples), or "unknown" (fewer than 2 samples, or no
## finite values at all - nothing to compare distributions across, so the
## question doesn't apply and guessing would be dishonest).
explore_data_status <- function(expr, dt) {
  m <- as.matrix(expr)
  finite_n <- sum(is.finite(m))
  if (ncol(m) < 2 || finite_n == 0) {
    return(list(
      state = "unknown", label = "Unknown",
      detail = "Fewer than two samples (or no finite expression values) are available, so between-sample normalization cannot be assessed."
    ))
  }
  if (identical(dt, "already_normalised")) {
    return(list(
      state = "normalized", label = "Already normalized",
      detail = "The active expression matrix has already undergone transformation/normalization. Additional normalization is not required unless explicitly selected below."
    ))
  }
  list(
    state = "not_normalized", label = "Not normalized",
    detail = "The active expression matrix has not been normalized. A normalization method can be selected below before downstream analysis."
  )
}

## PCA on a genes-in-rows expression matrix, dropping zero-variance genes
## first (prcomp's scale.=TRUE would otherwise error on them). Shared by
## Preprocessing (before/after ComBat) and Overview's QC (a single
## snapshot of whatever is currently loaded).
pca_of <- function(m, n_pc = 5) {
  m <- m[apply(m, 1, stats::sd) > 0, , drop = FALSE]
  p <- prcomp(t(m), scale. = TRUE)
  var_exp <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
  n_pc <- min(n_pc, ncol(p$x))
  df <- as.data.frame(p$x[, seq_len(n_pc), drop = FALSE])
  colnames(df) <- paste0("PC", seq_len(n_pc))
  df$sample <- rownames(df)
  list(df = df, var_exp = var_exp, n_pc = n_pc)
}

## Scree plot: % variance explained per PC, matching xOmicsShiny's
## fviz_eig()-based "Eigenvalues" view in its QC Plots module.
scree_plot <- function(var_exp, max_pc = 10) {
  n <- min(max_pc, length(var_exp))
  df <- data.frame(pc = factor(paste0("PC", seq_len(n)), levels = paste0("PC", seq_len(n))),
                    pct = var_exp[seq_len(n)])
  ggplot(df, aes(x = pc, y = pct)) +
    geom_col(fill = ARTHOMIX_COLORS$blue, width = 0.6) +
    geom_text(aes(label = paste0(pct, "%")), vjust = -0.5, size = 3, color = ARTHOMIX_COLORS$ink_secondary) +
    labs(x = NULL, y = "% variance explained") +
    theme_arthomix() +
    theme(panel.grid.major.x = element_blank())
}

## PCA scatter with a selectable PC pair, group confidence ellipses and
## optional sample labels - the same controls xOmicsShiny's PCA Plot tab
## offers (color-by, PC choice, ellipsoid, sample labels), minus
## 3D/loadings/convex-hull, which are out of scope here.
plot_pca_advanced <- function(pca_obj, meta, color_by, pc_x = 1, pc_y = 2,
                                title_suffix = NULL, show_ellipse = TRUE, show_labels = FALSE) {
  xcol <- paste0("PC", pc_x); ycol <- paste0("PC", pc_y)
  df <- merge(pca_obj$df[, c("sample", xcol, ycol)], meta[, c("sample", color_by)], by = "sample")
  xlab <- paste0(xcol, " (", pca_obj$var_exp[pc_x], "%)")
  ylab <- paste0(ycol, " (", pca_obj$var_exp[pc_y], "%)")
  is_num <- is.numeric(df[[color_by]])

  p <- ggplot(df, aes(x = .data[[xcol]], y = .data[[ycol]], color = .data[[color_by]])) +
    geom_point(size = 2.2, alpha = 0.85) +
    labs(x = xlab, y = ylab, color = color_by, title = title_suffix) +
    theme_arthomix()

  if (!is_num && show_ellipse) {
    ok_groups <- names(which(table(df[[color_by]]) >= 4))
    if (length(ok_groups) > 0) {
      p <- p + stat_ellipse(data = df[df[[color_by]] %in% ok_groups, , drop = FALSE],
                              type = "norm", level = 0.68, linewidth = 0.4, show.legend = FALSE)
    }
  }
  if (show_labels) {
    p <- p + ggrepel::geom_text_repel(aes(label = sample), size = 2.6, show.legend = FALSE, max.overlaps = 15)
  }

  if (is_num) p + scale_color_gradient(low = "#cde2fb", high = ARTHOMIX_COLORS$blue)
  else p + scale_color_manual(values = arthomix_pair(df[[color_by]]))
}

## Feature (gene/probe) overlap across 2+ datasets - the same idea as
## scripts/00_shared/03_dataset_gene_overlap_venn.R's ggVennDiagram call,
## generalised to any number of input datasets instead of just the two
## training platforms. `sets` is a named list of character vectors (feature
## IDs measured on each dataset). Area-proportional circles/ellipses for 2-4
## sets (ggVennDiagram); beyond that, circles stop being readable, so it
## falls back to an UpSet-style bar chart of exact-combination sizes, capped
## at the `max_bars` largest combinations so 5-6 sets don't produce 60+ bars.
## Every non-empty "belongs to exactly these sets and no others" region,
## with its size - the same exact-combination arithmetic a proportional
## Venn encodes visually, but as a plain table any downstream code (a
## downloadable region table, a >7-set fallback) can use directly.
overlap_region_sizes <- function(sets, max_regions = NULL) {
  nms <- names(sets)
  combos <- unlist(lapply(seq_along(nms), function(k) utils::combn(nms, k, simplify = FALSE)), recursive = FALSE)
  sizes <- vapply(combos, function(cmb) {
    inc <- Reduce(intersect, sets[cmb])
    excl_from <- setdiff(nms, cmb)
    if (length(excl_from) > 0) inc <- setdiff(inc, Reduce(union, sets[excl_from]))
    length(inc)
  }, numeric(1))
  df <- data.frame(combination = vapply(combos, paste, character(1), collapse = " & "), n_features = sizes)
  df <- df[df$n_features > 0, , drop = FALSE]
  df <- df[order(-df$n_features, df$combination), , drop = FALSE]
  if (!is.null(max_regions) && nrow(df) > max_regions) df <- df[seq_len(max_regions), ]
  rownames(df) <- NULL
  df
}

## Feature (gene/probe) overlap across 2-7 datasets - ggVennDiagram handles
## every set count actually offered in this app (MAX_PP_SOURCES = 6) as a
## true proportional Venn, not just 2-4, so there's no separate bar-chart
## fallback to keep in sync; `label = "both"` shows count and % of the total
## per region, and `set_color` gives each dataset its own outline/label
## colour from the app's own palette instead of ggVennDiagram's default
## black, so which circle is which is legible without cross-referencing a
## legend. Above 7 sets (not reachable from this app's own UI) falls back to
## overlap_region_sizes() as a bar chart, rather than erroring.
## fill_low/fill_high override the region-size gradient's two end colors
## (default: pale blue to ARTHOMIX_COLORS$blue) - e.g. mod_candidates.R
## passes a green gradient for its Female panel and a brown one for Male,
## matching this project's own reference figures
## (fig_venn_female_disease_candidates.png / fig_venn_male_disease_candidates.png
## in Research_Q2_TRANSCRIPTOMICS_sexstratified_COPY/results/figures/), which
## use the same single-hue-by-region-size scheme, just one hue per sex so
## the two panels read as visibly distinct at a glance. Every other caller
## (Preprocessing's dataset-gene-overlap Venn) omits both and keeps the
## original blue.
draw_overlap_venn <- function(sets, title = NULL, max_bars = 30, fill_low = "#EAF3FB", fill_high = ARTHOMIX_COLORS$blue) {
  sets <- sets[lengths(sets) > 0]
  validate(need(length(sets) >= 2, "Need at least two non-empty datasets to compare."))

  common_n <- length(Reduce(intersect, sets))
  auto_title <- title %||% sprintf("Common to all %d datasets: %s features", length(sets), format(common_n, big.mark = ","))
  pal <- c(ARTHOMIX_COLORS$blue, ARTHOMIX_COLORS$orange, ARTHOMIX_COLORS$aqua,
            ARTHOMIX_COLORS$violet, ARTHOMIX_COLORS$magenta, ARTHOMIX_COLORS$yellow, ARTHOMIX_COLORS$red)

  if (length(sets) <= 7) {
    p <- ggVennDiagram::ggVennDiagram(
      sets, label = "both", label_alpha = 0, label_size = 3.3, label_percent_digit = 1,
      edge_size = 0.9, set_size = 4.2, set_color = pal[seq_along(sets)]
    )
    ## ggVennDiagram centers each set-name label (hjust = 0.5) on an anchor
    ## point just outside its own circle - fine when circles are roughly
    ## same-sized and spread apart, but when one set is mostly contained in
    ## the other (a highly asymmetric 2-set overlap, e.g. 65%/14%/21%), that
    ## anchor sits close enough to the OTHER circle's edge that half the
    ## centered label text swings back over it. Right-aligning (hjust = 1)
    ## keeps the anchor in place but grows the label outward, away from both
    ## circles, instead - identified by its `label` aesthetic mapping to
    ## `.data$name` (the category-name layer ggVennDiagram itself adds),
    ## not the region count/percent labels. Scoped to exactly 2 sets: with
    ## 3+ circles the default centered anchors are already clear of each
    ## other, and forcing hjust = 1 there instead dragged adjacent set
    ## names (e.g. "LASSO" / "RandomForest") into one another at the top
    ## of the diagram.
    two_set_asymmetric <- length(sets) == 2
    if (two_set_asymmetric) {
      for (i in seq_along(p$layers)) {
        l <- p$layers[[i]]
        if (inherits(l$geom, "GeomText") && !is.null(l$mapping$label) &&
            grepl("name", rlang::as_label(l$mapping$label), fixed = TRUE)) {
          p$layers[[i]]$aes_params$hjust <- 1
        }
      }
    }
    ## Every set count uses the same modest scale expansion now - it only
    ## needs to clear the circles themselves, not reserve room for growing
    ## labels (that used to be a 110%/35% expansion for the 2-set case,
    ## which shrank the circles down to a sliver of the panel and, applied
    ## to every set count, is what crowded 3+-set region-count labels into
    ## illegible overlapping text).
    p <- p +
      scale_fill_gradient(low = fill_low, high = fill_high, name = "Features") +
      scale_x_continuous(expand = expansion(mult = 0.12)) +
      scale_y_continuous(expand = expansion(mult = 0.12)) +
      ## coord_equal(), not coord_cartesian(): ggVennDiagram's circles are
      ## only actually round if the x/y scales stay 1:1 - coord_cartesian()
      ## would silently replace ggVennDiagram's own coord_equal() (ggplot2
      ## keeps only one coordinate system per plot), stretching the circles
      ## to whatever aspect ratio the plot panel happens to render at.
      coord_equal(clip = "off") +
      labs(title = auto_title) +
      theme(legend.position = "right", legend.text = element_text(size = 9),
            plot.title = element_text(face = "bold", size = 13, hjust = 0))
    ## The 2-set hjust = 1 fix above grows both category-name labels
    ## leftward, outside the panel, with no upper bound - scale expansion
    ## (a FRACTION of the data range) can't reserve enough absolute room
    ## for that reliably, which is what was clipping long names like
    ## "WGCNA module background" down to "ule background" even after the
    ## expansion was made generous. A fixed left plot margin (an absolute
    ## amount, in points, independent of the data range) reserves real
    ## room for it instead, with clip = "off" above letting the label draw
    ## into that margin - without shrinking the circles the way growing
    ## the scale expansion further did.
    if (two_set_asymmetric) {
      p <- p + theme(plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = 5.5, l = 130, unit = "pt"))
    }
    p
  } else {
    df <- overlap_region_sizes(sets, max_regions = max_bars)
    ggplot(df, aes(x = reorder(combination, n_features), y = n_features)) +
      geom_col(fill = ARTHOMIX_COLORS$blue, width = 0.7) +
      geom_text(aes(label = format(n_features, big.mark = ",")), hjust = -0.15, size = 3, color = ARTHOMIX_COLORS$ink_secondary) +
      coord_flip(clip = "off") +
      labs(x = NULL, y = "Features in exactly this combination of datasets", title = auto_title) +
      theme_arthomix() +
      theme(panel.grid.major.y = element_blank())
  }
}

## A thin, rank-ordered bar chart shared by every QC metric: neutral blue
## for samples within range, status red for flagged ones, no per-sample
## x-axis labels (illegible past ~30 samples, and the point is the shape
## of the distribution plus which bars are red).
qc_bar_plot <- function(df, y_col, flag_col, y_label, subtitle = NULL) {
  df$.flag <- ifelse(df[[flag_col]], "Flagged", "Within range")
  ggplot(df, aes(x = reorder(sample, .data[[y_col]]), y = .data[[y_col]], fill = .flag)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = c("Within range" = ARTHOMIX_COLORS$blue, "Flagged" = ARTHOMIX_STATUS$critical)) +
    labs(x = NULL, y = y_label, fill = NULL, subtitle = subtitle) +
    theme_arthomix() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank())
}

## ---------------------------------------------------------------------------
## Top level module registry (Home > Modules > one card per omics layer)
## ---------------------------------------------------------------------------

MODULE_REGISTRY <- list(
  list(
    id = "transcriptomics", tab = "transcriptomics",
    title = "Transcriptomics",
    tagline = "The core pipeline: sex-stratified differential expression, WGCNA, Mendelian randomisation and colocalisation, biomarker panel selection, and cross-tissue / cross-ancestry validation.",
    icon = "dna", status = "available", kind = "Single-omics"
  ),
  list(
    id = "methylomics", tab = "methylomics",
    title = "Methylomics",
    tagline = "Analyze DNA methylation data to identify differentially methylated sites and regions, discover methylation patterns, perform feature selection and network analysis, and generate biologically relevant insights.",
    icon = "circle-nodes", status = "available", kind = "Single-omics"
  ),
  list(
    id = "crossomics", tab = "crossomics",
    title = "Cross-Omics",
    tagline = "Browse the biomarker-convergence tables joining Transcriptomics' eQTL-MR genes and Methylomics' mQTL-MR CpGs; live convergence and cross-omics MR sub-modules are coming next.",
    icon = "arrows-left-right", status = "available", kind = "Multi-omics"
  ),
  list(
    id = "multiomics", tab = "multiomics",
    title = "Multi-Omics",
    tagline = "Browse the pipeline's own DIABLO and SNF integration, joint biomarker discovery, gene<->CpG concordance, and pathway enrichment across transcriptomics + methylomics.",
    icon = "layer-group", status = "available", kind = "Multi-omics"
  ),
  list(
    id = "arthochat", tab = "arthochat",
    title = "ArthOChat",
    tagline = "Ask anything about your dataset, your results so far, or the underlying methodology, grounded in a live PubMed search with citations. Sees every omics module, not just one.",
    icon = "comments", status = "available", kind = "Assistant"
  )
)

## A single, app-wide ArthOChat session lives on its own top-level tab (see
## ui.R's navbarPage and server.R) rather than a separate chat instance
## nested inside every sub-module - one shared history that already sees
## whatever dataset/results are current, regardless of which module surfaces
## this link. Drops straight into the top-level "arthochat" tab via a client
## side click on that navbar link (no server-side session plumbing needed to
## reach across module boundaries), matching the tab's `value` set in ui.R.
## ArthOChat lives in a persistent slide-out drawer (see ui.R), not a
## navbarPage tab, so opening it from anywhere - a submodule shortcut, the
## header button, the Modules landing card, a server-side search match -
## never navigates away from whatever the user is currently looking at.
## Every "open ArthOChat" trigger in the app shares this exact JS so
## there's one behavior to keep in sync, not several copies that could
## drift: client-side triggers (onclick=) use ARTHOCHAT_DRAWER_OPEN_JS
## directly; server-side ones (e.g. server.R's header-search handler) run
## ARTHOCHAT_DRAWER_OPEN_JS_STATEMENT via shinyjs::runjs().
ARTHOCHAT_DRAWER_OPEN_JS_STATEMENT <- "document.getElementById('arthochat_drawer').classList.add('open')"
ARTHOCHAT_DRAWER_OPEN_JS <- paste0(ARTHOCHAT_DRAWER_OPEN_JS_STATEMENT, "; return false;")

arthochat_shortcut_ui <- function(hint = NULL, compact = FALSE) {
  jump_link <- tags$a(
    "Open ArthOChat", href = "#", class = paste("btn btn-primary btn-sm", if (compact) "btn-block"),
    onclick = ARTHOCHAT_DRAWER_OPEN_JS
  )
  if (compact) {
    return(div(
      class = "arthochat-shortcut arthochat-shortcut-compact",
      div(icon("comments"), strong(" Ask ArthOChat")),
      p(hint %||% "Have a methodology question? Ask ArthOChat."),
      jump_link
    ))
  }
  div(
    class = "arthochat-shortcut",
    icon("comments"),
    span(
      strong("Have a methodology question? "),
      hint %||% "Ask ArthOChat. It can see this dataset and cites PubMed for any scientific claim.",
      " "
    ),
    jump_link
  )
}

## ---------------------------------------------------------------------------
## Transcriptomics submodules
## Each submodule's config, UI wrapper and server wrapper live in their own
## file: R/transcriptomics/mod_<id>.R (Methylomics' own submodules live the
## same way under R/methylomics/ - see R/0_load_omics_modules.R for how
## both subfolders get sourced, since shiny's own R/ autoloader doesn't
## recurse into them). They are assembled into TX_MODULES by
## R/submodules_registry.R, which is sourced after all of them (see the
## sourcing note at the top of that file). Section numbers refer to
## results/METHODS_00_INDEX.md in the data folder.
## ---------------------------------------------------------------------------
