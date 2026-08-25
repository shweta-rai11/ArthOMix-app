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
## such throughout, never claimed as proven causation, and explicitly
## flagged as impossible to validate via heterogeneity with only one
## instrument per exposure.

CX_MR_PRECOMPUTED_FILE <- file.path(CX_RESULTS_DIR, "mr_stage_eqtl_significant_genes_mqtl_mr.csv")
CX_MR_DATA_AVAILABLE <- CX_DATA_AVAILABLE && file.exists(CX_MR_PRECOMPUTED_FILE)
CX_MR_OUTCOME_ID <- "GCST90132223"
CX_MR_OUTCOME_NAME <- "Rheumatoid arthritis (Ishigaki et al. 2022, GCST90132223)"

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
  ## downstream comparison (cx_mr_classify_tier()) expects a logical.
  if (!is.logical(df$steiger_dir)) df$steiger_dir <- as.logical(df$steiger_dir)
  list(ok = TRUE, df = df, error = NULL)
}

## ---------------------------------------------------------------------------
## Tier classification - NOT computed by either pipeline script; implemented
## here from the exact rule stated in
## cross_Omics_Sexstratified_COPY/results/CROSS_OMICS_REPORT.md sections
## 4.2.6-4.2.7:
##
##   "credible" mQTL-MR = BH-FDR significant AND gene not MHC-flagged AND
##   Steiger filtering supports the intended (methylation -> disease)
##   direction.
##
##   Tier 1 = DEG significant AND (DMP or DMR significant) AND credible
##   mQTL-MR.
##   Tier 2 = DEG significant AND exactly one of {(DMP or DMR significant),
##   credible mQTL-MR}.
##   Tier 3 = everything else.
##
## This is quoted verbatim in the module's own UI (never a black box).
## ---------------------------------------------------------------------------

CX_MR_TIER_LEVELS <- c("Tier 1", "Tier 2", "Tier 3")

CX_MR_TIER_DEFAULT_FILTERS <- list(fdr_cutoff = 0.05, require_steiger = TRUE, exclude_mhc = TRUE)

## Rendered from the actual filter values in use, not a fixed string - so
## the UI always shows the rule that was really applied, not just the
## report's own original defaults.
cx_mr_tier_rule_text <- function(filters = CX_MR_TIER_DEFAULT_FILTERS) {
  f <- utils::modifyList(CX_MR_TIER_DEFAULT_FILTERS, filters %||% list())
  c(
    sprintf("\"Credible\" mQTL-MR evidence = FDR < %s (Benjamini-Hochberg, as computed by the pipeline's own MR run)%s%s.",
            f$fdr_cutoff,
            if (isTRUE(f$exclude_mhc)) " AND the gene is not MHC-flagged (per the eQTL-MR panel's own MHC_gene flag)" else " (MHC-flagged genes are NOT excluded at this setting)",
            if (isTRUE(f$require_steiger)) " AND Steiger filtering supports the methylation -> disease direction" else " (Steiger direction is NOT required at this setting)"),
    "Tier 1: DEG significant AND (DMP genome-wide-significant OR DMR significant) AND credible mQTL-MR evidence.",
    "Tier 2: DEG significant AND exactly one of - (DMP or DMR significant), credible mQTL-MR evidence.",
    "Tier 3: every other gene in the join (including genes with DEG evidence only, or with methylation evidence but no qualifying mQTL-MR support).",
    "Source: cross_Omics_Sexstratified_COPY/results/CROSS_OMICS_REPORT.md, sections 4.2.6-4.2.7. Implemented here in code for the first time - neither pipeline script computes Tier; this module reproduces the documented rule, on whichever filter values are currently set (defaults match the report's own)."
  )
}

## `mr_df` is cx_mr_load_precomputed()$df (needs gene, FDR, steiger_dir -
## `mr_significant` is recomputed here at `filters$fdr_cutoff` rather than
## trusted as-is, so relabeling significance is instant). `join_df` is
## cx_bc_load_precomputed()$df (optionally cx_bc_relabel()'d) for the SAME
## sex (needs gene, DEG_significant, methylation_significant,
## eQTL_MHC_region) - if NULL, DEG/methylation evidence is reported "Not
## available" rather than guessed, and every gene falls back to Tier 3 (the
## only tier reachable without that evidence). `filters` lets the "credible
## mQTL-MR" definition itself be reconfigured (FDR cutoff, whether
## MHC-flagged genes are excluded, whether Steiger direction is required) -
## all cheap relabeling of already-computed values, so this whole function
## is fast enough to call on every filter change, no "Run" button needed.
cx_mr_classify_tier <- function(mr_df, join_df = NULL, filters = CX_MR_TIER_DEFAULT_FILTERS) {
  f <- utils::modifyList(CX_MR_TIER_DEFAULT_FILTERS, filters %||% list())
  mr_best <- cx_bc_dedup_min(mr_df, "gene", "pval")  ## best (lowest-p) instrument per gene, matching the pipeline's own per-gene summary
  mr_best$mr_significant <- !is.na(mr_best$FDR) & mr_best$FDR < f$fdr_cutoff
  mr_best$credible_mQTL_MR <- mr_best$mr_significant
  if (isTRUE(f$require_steiger)) mr_best$credible_mQTL_MR <- mr_best$credible_mQTL_MR & !(mr_best$steiger_dir %in% FALSE)

  if (is.null(join_df)) {
    out <- data.frame(gene = mr_best$gene, DEG_significant = NA, methylation_significant = NA,
                       credible_mQTL_MR = mr_best$credible_mQTL_MR, tier = "Tier 3",
                       tier_evidence_available = FALSE, stringsAsFactors = FALSE)
    return(out)
  }
  if (isTRUE(f$exclude_mhc)) {
    eqtl_mhc <- unique(join_df[, c("gene", "eQTL_MHC_region"), drop = FALSE])
    mr_best <- merge(mr_best, eqtl_mhc, by = "gene", all.x = TRUE)
    mr_best$credible_mQTL_MR <- mr_best$credible_mQTL_MR & !(mr_best$eQTL_MHC_region %in% TRUE)
  }

  ev <- unique(join_df[, c("gene", "DEG_significant", "methylation_significant"), drop = FALSE])
  out <- merge(ev, mr_best[, c("gene", "credible_mQTL_MR"), drop = FALSE], by = "gene", all.x = TRUE)
  out$credible_mQTL_MR[is.na(out$credible_mQTL_MR)] <- FALSE

  n_secondary <- rowSums(cbind(out$methylation_significant %in% TRUE, out$credible_mQTL_MR %in% TRUE))
  out$tier <- "Tier 3"
  out$tier[out$DEG_significant %in% TRUE & n_secondary == 1] <- "Tier 2"
  out$tier[out$DEG_significant %in% TRUE & n_secondary >= 2] <- "Tier 1"
  out$tier <- factor(out$tier, levels = CX_MR_TIER_LEVELS)
  out$tier_evidence_available <- TRUE
  out
}

cx_mr_build_provenance <- function(filters, n_instruments, join_sex, run_at) {
  f <- utils::modifyList(CX_MR_TIER_DEFAULT_FILTERS, filters %||% list())
  c(
    sprintf("Source: cross_Omics_Sexstratified_COPY/results/%s - the pipeline's own already-run MR results (no re-computation performed here).", basename(CX_MR_PRECOMPUTED_FILE)),
    "Exposure: GoDMC (Min et al. 2021) cis-mQTL summary statistics, one lead cis-SNP per CpG.",
    sprintf("Outcome: %s.", CX_MR_OUTCOME_NAME),
    "Method: single-instrument Wald ratio (TwoSampleMR::mr_wald_ratio) - one lead cis-SNP per CpG; not IVW/MR-Egger, which need >=2 instruments per exposure.",
    sprintf("Instruments in the precomputed table: %s", format(n_instruments, big.mark = ",")),
    sprintf("\"Credible\" evidence FDR threshold: %s", f$fdr_cutoff),
    sprintf("MHC-flagged genes excluded from \"credible\": %s", if (isTRUE(f$exclude_mhc)) "yes" else "no"),
    sprintf("Steiger direction required for \"credible\": %s", if (isTRUE(f$require_steiger)) "yes" else "no"),
    sprintf("DEG/DMP/DMR evidence joined from: Biomarker Convergence's %s precomputed table.", toupper(join_sex %||% "(not available)")),
    "These are Mendelian randomization estimates, valid under the standard instrumental-variable assumptions (relevance, independence, exclusion restriction). A single instrument per exposure cannot be tested for validity via heterogeneity - interpret individual estimates cautiously, as an association consistent with a causal effect, not as proof of one.",
    sprintf("Run at: %s", run_at %||% "(not run yet)")
  )
}
