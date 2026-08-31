## R/crossomics/crossomics_mrstage_helpers.R
## Pure data-processing logic for the "Cross-Omics MR" Cross-Omics
## sub-module (mod_cross_mr_stage.R).
##
## The MR analysis itself (single-instrument Wald ratio, GoDMC cis-mQTL
## exposure -> Ishigaki et al. 2022 RA GWAS outcome) was already run once by
## the pipeline's own cross_Omics_Sexstratified_COPY/scripts/02_mr_stage_
## cross_omics.R, and its output already sits in
## cross_Omics_Sexstratified_COPY/results/mr_stage_eqtl_significant_genes_
## mqtl_mr.csv. This module never re-runs that computation - it loads that
## file (instant) and INTEGRATES it with the DEG/DMP/DMR evidence from
## crossomics_biomarkerconv_helpers.R's precomputed join, letting the
## significance/MHC/Steiger filters that decide "credible" evidence and Tier
## be reconfigured live (cx_mr_classify_tier()) by relabeling already-
## computed values - never a new statistical test, and never a file read
## from outside cross_Omics_Sexstratified_COPY (the raw GoDMC/Ishigaki
## inputs that produced this file live elsewhere and are not touched here).
##
## MR estimates are genuine causal estimates under the standard IV
## assumptions (relevance, independence, exclusion restriction) - worded as
## such throughout, never claimed as proven causation.
##
## CAVEAT (not previously documented here): this file is one row per
## CpG-instrument, so an individual gene with only one CpG-instrument indeed
## cannot be heterogeneity-tested. But several genes in the precomputed file
## have many independent CpG-instruments each (the MHC genes especially -
## HLA-DRB1 alone has dozens) - for those genes, standard per-gene
## heterogeneity/pleiotropy testing (Cochran's Q, MR-Egger intercept across
## that gene's own instruments) is statistically possible and is exactly
## what the field-standard STROBE-MR checklist expects before treating a
## multi-instrument gene's estimate as robust. This module does not
## currently aggregate per gene or run that test - each CpG-instrument is
## still shown/filtered as an independent row - so "impossible to test" was
## an overstatement for the genes where it matters most. No instrument-
## strength (F-statistic) figure is available in the precomputed file
## either. Both are real, currently-unaddressed gaps; the MHC-region
## caveat above is the substitute safeguard this module ships with instead.

CX_MR_PRECOMPUTED_FILE <- file.path(CX_RESULTS_DIR, "mr_stage_eqtl_significant_genes_mqtl_mr.csv")
CX_MR_DATA_AVAILABLE <- CX_DATA_AVAILABLE && file.exists(CX_MR_PRECOMPUTED_FILE)
CX_MR_OUTCOME_ID <- "GCST90132223"
CX_MR_OUTCOME_NAME <- "Rheumatoid arthritis (Ishigaki et al. 2022, GCST90132223)"

## Volcano-plot direction levels (mod_cross_mr_stage.R's build_volcano_plot()) -
## b is the MR log-odds estimate, so b > 0 means the exposure increases RA
## risk ("Up") and b < 0 means it's protective ("Down"), the same up/down
## framing as a differential-expression volcano. Direction is only
## meaningful among significant points; every non-significant point is
## grouped as "Not significant" regardless of its sign.
CX_MR_VOLCANO_LEVELS <- c("Up (risk, b > 0)", "Down (protective, b < 0)", "Not significant")

## ---------------------------------------------------------------------------
## Precomputed MR results (the module's only data source) - the pipeline's
## own already-run Wald-ratio output, one row per CpG-instrument: cpg,
## gene, SNP, b, se, pval, OR, OR_lo, OR_hi, FDR, mr_significant,
## steiger_dir, steiger_pval.
## ---------------------------------------------------------------------------

cx_mr_load_precomputed <- function() {
  if (!CX_MR_DATA_AVAILABLE) {
    return(list(ok = FALSE, df = NULL, error = "The pipeline's precomputed Cross-Omics MR results are not available in this deployment."))
  }
  df <- tryCatch(as.data.frame(data.table::fread(CX_MR_PRECOMPUTED_FILE, showProgress = FALSE)), error = function(e) e)
  if (inherits(df, "error")) return(list(ok = FALSE, df = NULL, error = paste("Could not read the precomputed MR results:", conditionMessage(df))))
  if (nrow(df) == 0) return(list(ok = FALSE, df = NULL, error = "The precomputed MR results file has no rows."))
  ## steiger_dir is stored as a string ("TRUE"/"FALSE") in the CSV; every
  ## downstream comparison expects a logical.
  if (!is.logical(df$steiger_dir)) df$steiger_dir <- as.logical(df$steiger_dir)
  list(ok = TRUE, df = df, error = NULL)
}

## ---------------------------------------------------------------------------
## "Upload your own data" - MR instrument results, as an alternative to the
## precomputed file above. Same list(ok, df, error) contract, same
## cx_read_table() (crossomics_integration_upload.R) parser every other
## upload path in this app uses. Only "gene" and "pval" are required - every
## other column (cpg, SNP, nsnp, b, se, OR, OR_lo, OR_hi, FDR, steiger_dir,
## steiger_pval) is optional and simply reads as unavailable if omitted,
## never fabricated. FDR is recomputed via BH from pval when not supplied -
## a standard multiple-testing correction on the caller's own p-values, not
## a new statistical test.
## ---------------------------------------------------------------------------

CX_MR_REQUIRED_UPLOAD_COLS <- c("gene", "pval")

cx_mr_load_upload <- function(datapath, filename) {
  res <- cx_read_table(datapath, filename)
  if (!res$ok) return(list(ok = FALSE, df = NULL, error = res$error))
  df <- res$df
  missing <- setdiff(CX_MR_REQUIRED_UPLOAD_COLS, colnames(df))
  if (length(missing) > 0) {
    return(list(ok = FALSE, df = NULL, error = sprintf(
      "This MR results file is missing required column(s): %s. One row per CpG-instrument; \"gene\" and \"pval\" are required.",
      paste(missing, collapse = ", "))))
  }
  if (nrow(df) == 0) return(list(ok = FALSE, df = NULL, error = "The uploaded MR results file has no rows."))
  df$gene <- as.character(df$gene)
  df <- df[!is.na(df$gene) & nzchar(df$gene), , drop = FALSE]
  if (nrow(df) == 0) return(list(ok = FALSE, df = NULL, error = "No row in the uploaded MR results file has a non-empty gene value."))
  df$pval <- suppressWarnings(as.numeric(df$pval))
  for (cl in intersect(c("b", "se", "OR", "OR_lo", "OR_hi", "FDR", "nsnp", "steiger_pval"), colnames(df))) {
    df[[cl]] <- suppressWarnings(as.numeric(df[[cl]]))
  }
  if (!"FDR" %in% colnames(df) || all(is.na(df$FDR))) {
    ok <- !is.na(df$pval)
    df$FDR <- NA_real_
    df$FDR[ok] <- stats::p.adjust(df$pval[ok], method = "BH")
  }
  if (!"steiger_dir" %in% colnames(df)) df$steiger_dir <- NA
  else if (!is.logical(df$steiger_dir)) df$steiger_dir <- as.logical(df$steiger_dir)
  if (!"cpg" %in% colnames(df)) df$cpg <- NA_character_
  if (!"SNP" %in% colnames(df)) df$SNP <- NA_character_
  list(ok = TRUE, df = df, error = NULL)
}

## ---------------------------------------------------------------------------
## Evidence-combination categories - replaces this module's original Tier
## 1/2/3 priority ranking (which followed cross_Omics_Sexstratified_COPY/
## results/CROSS_OMICS_REPORT.md sections 4.2.6-4.2.7) with 5 independently-
## checked evidence combinations, one per tab, requested explicitly rather
## than derived from the pipeline report. A gene can match more than one
## category (e.g. both DEG-eQTL and DMP-mQTL at once) - these are NOT
## mutually exclusive tiers, unlike the scheme they replace.
##
## Every significance flag used below (DEG_significant, DMP_genomewide_
## significant, DMR_significant, mQTL_MR_significant, eQTL_MR_significant)
## is join_df's own column, already computed by cx_bc_relabel()
## (crossomics_biomarkerconv_helpers.R) at Biomarker Convergence's own
## default thresholds (CX_BC_DEFAULT_PARAMS) - the exact same flags that
## module's own eQTL-MR/mQTL-MR/eQTL-mQTL tabs use, including the mQTL-MR
## panel-membership backfill documented there. There is no separate
## "credible" mQTL-MR definition here anymore (that lived behind this
## module's old FDR/MHC/Steiger filter inputs, removed along with Tier
## 1/2/3) - every category below is pure relabeling of join_df, so it's
## just as instant as the scheme it replaces, but reconfigured from the
## Biomarker Convergence tab's thresholds, not from filters on this page.
## ---------------------------------------------------------------------------

## eQTL_MHC_region is added to every category below that includes eQTL-MR evidence.
## The MHC region (chr6, ~25-34Mb) is the single most notorious horizontal-pleiotropy
## hotspot in autoimmune-disease genetics - extreme LD across dozens of genes makes any
## one gene's instrument a poor proxy for that gene specifically, so an MR estimate
## landing here is far more likely to reflect regional confounding than the gene's own
## causal effect. The flag was already computed upstream (Biomarker Convergence's own
## precomputed pipeline column) but, before this fix, was carried through to the UI as a
## silent display column with no warning and no way to see at a glance how much of a
## category's "convergent evidence" sits in this one hotspot - see the MHC callout in
## mod_cross_mr_stage.R's per-category UI.
CX_MR_CATEGORIES <- list(
  list(id = "deg_dmp_qtl", tab = "DEG-DMP-QTL",
       rule_text = "DEG significant AND DMP genome-wide-significant AND mQTL-MR significant AND eQTL-MR significant.",
       cols = c("gene", "DEG_logFC", "DEG_direction", "DEG_adjP", "DMP_top_cpg", "DMP_dbeta", "DMP_direction", "DMP_fdr_bacon",
                "mQTL_candidate_cpg", "mQTL_MR_beta", "mQTL_MR_pval", "eQTL_MR_OR", "eQTL_MR_direction", "eQTL_MR_FDR", "eQTL_MHC_region")),
  list(id = "deg_dmr_qtl", tab = "DEG-DMR-QTL",
       rule_text = "DEG significant AND DMR significant AND mQTL-MR significant AND eQTL-MR significant.",
       cols = c("gene", "DEG_logFC", "DEG_direction", "DEG_adjP", "DMR_id", "DMR_meandiff", "DMR_direction", "DMR_fdr",
                "mQTL_candidate_cpg", "mQTL_MR_beta", "mQTL_MR_pval", "eQTL_MR_OR", "eQTL_MR_direction", "eQTL_MR_FDR", "eQTL_MHC_region")),
  list(id = "deg_eqtl", tab = "DEG-eQTL",
       rule_text = "DEG significant AND eQTL-MR significant.",
       cols = c("gene", "DEG_logFC", "DEG_direction", "DEG_adjP", "eQTL_MR_OR", "eQTL_MR_direction", "eQTL_MR_FDR", "eQTL_MHC_region")),
  list(id = "dmp_mqtl", tab = "DMP-mQTL",
       rule_text = "DMP genome-wide-significant AND mQTL-MR significant.",
       cols = c("gene", "DMP_top_cpg", "DMP_dbeta", "DMP_direction", "DMP_fdr_bacon", "mQTL_candidate_cpg", "mQTL_MR_beta", "mQTL_MR_pval")),
  list(id = "dmr_mqtl", tab = "DMR-mQTL",
       rule_text = "DMR significant AND mQTL-MR significant.",
       cols = c("gene", "DMR_id", "DMR_meandiff", "DMR_direction", "DMR_fdr", "mQTL_candidate_cpg", "mQTL_MR_beta", "mQTL_MR_pval"))
)

## ---------------------------------------------------------------------------
## "Upload your own data" - a gene-level evidence table, as an alternative
## to Biomarker Convergence's own precomputed table (either instead of, or
## even when Biomarker Convergence itself has no upload option for DEG/DMP/
## DMR). Only "gene" is required; every DEG/DMP/DMR/mQTL-MR/eQTL-MR column
## below is optional and defaults to NA (never fabricated) if omitted -
## missing columns simply mean that layer reads as "not significant/not
## evaluated" once run through cx_bc_relabel() (crossomics_biomarkerconv_
## helpers.R, reused unchanged here - same column names, same default
## thresholds, so an uploaded table plugs into cx_mr_classify_categories()
## exactly like the preloaded one does).
## ---------------------------------------------------------------------------

CX_MR_EVIDENCE_UPLOAD_NUMERIC_COLS <- c("DEG_logFC", "DEG_adjP", "DMP_dbeta", "DMP_fdr_bacon", "DMR_meandiff", "DMR_fdr",
                                         "mQTL_MR_beta", "mQTL_MR_pval", "eQTL_MR_OR", "eQTL_MR_FDR")
CX_MR_EVIDENCE_UPLOAD_TEXT_COLS <- c("DEG_direction", "DMP_top_cpg", "DMP_direction", "DMR_id", "DMR_direction",
                                      "mQTL_candidate_cpg", "eQTL_MR_direction")

cx_mr_load_evidence_upload <- function(datapath, filename) {
  res <- cx_read_table(datapath, filename)
  if (!res$ok) return(list(ok = FALSE, df = NULL, error = res$error))
  df <- res$df
  if (!"gene" %in% colnames(df)) {
    return(list(ok = FALSE, df = NULL, error = "This evidence file is missing the required \"gene\" column. One row per gene."))
  }
  if (nrow(df) == 0) return(list(ok = FALSE, df = NULL, error = "The uploaded evidence file has no rows."))
  df$gene <- as.character(df$gene)
  df <- df[!is.na(df$gene) & nzchar(df$gene), , drop = FALSE]
  if (nrow(df) == 0) return(list(ok = FALSE, df = NULL, error = "No row in the uploaded evidence file has a non-empty gene value."))
  for (cl in intersect(CX_MR_EVIDENCE_UPLOAD_NUMERIC_COLS, colnames(df))) df[[cl]] <- suppressWarnings(as.numeric(df[[cl]]))
  for (cl in CX_MR_EVIDENCE_UPLOAD_NUMERIC_COLS) if (!cl %in% colnames(df)) df[[cl]] <- NA_real_
  for (cl in CX_MR_EVIDENCE_UPLOAD_TEXT_COLS) if (!cl %in% colnames(df)) df[[cl]] <- NA_character_
  ## No panel-membership column to read (unlike Biomarker Convergence's own
  ## eQTL upload path, and unlike the preloaded pipeline table where
  ## in_eQTL_MR_panel was already FDR<0.05-filtered upstream before this
  ## app ever saw it - see cx_bc_relabel()'s own comment on
  ## eQTL_MR_significant). cx_bc_relabel() reads in_eQTL_MR_panel directly
  ## as "eQTL-MR significant" with no threshold of its own, so it has to be
  ## computed AT that same 0.05 threshold here - merely having a non-NA
  ## eQTL_MR_FDR value is not enough (a gene with eQTL_MR_FDR = 0.9 must
  ## NOT count as "in panel").
  df$in_eQTL_MR_panel <- !is.na(df$eQTL_MR_FDR) & df$eQTL_MR_FDR < 0.05
  list(ok = TRUE, df = df, error = NULL)
}

## `join_df` is cx_bc_relabel()'d Biomarker Convergence data for the
## selected sex. Returns NULL if join_df is NULL (evidence not available -
## matches every other cx_*() loader's fail-soft contract). Otherwise
## returns a named list (by CX_MR_CATEGORIES id) of data.frames, each
## pre-filtered to the genes matching that category's rule and restricted
## to that category's own relevant columns (only columns actually present
## in join_df are kept - never fabricates a column the data doesn't have).
cx_mr_classify_categories <- function(join_df) {
  if (is.null(join_df)) return(NULL)
  df <- join_df
  flags <- list(
    deg  = df$DEG_significant %in% TRUE,
    dmp  = df$DMP_genomewide_significant %in% TRUE,
    dmr  = df$DMR_significant %in% TRUE,
    mqtl = df$mQTL_MR_significant %in% TRUE,
    eqtl = df$eQTL_MR_significant %in% TRUE
  )
  matches <- list(
    deg_dmp_qtl = flags$deg & flags$dmp & flags$mqtl & flags$eqtl,
    deg_dmr_qtl = flags$deg & flags$dmr & flags$mqtl & flags$eqtl,
    deg_eqtl    = flags$deg & flags$eqtl,
    dmp_mqtl    = flags$dmp & flags$mqtl,
    dmr_mqtl    = flags$dmr & flags$mqtl
  )
  setNames(
    lapply(CX_MR_CATEGORIES, function(cat) {
      cols <- intersect(cat$cols, colnames(df))
      df[matches[[cat$id]] %in% TRUE, cols, drop = FALSE]
    }),
    vapply(CX_MR_CATEGORIES, function(cat) cat$id, character(1))
  )
}

