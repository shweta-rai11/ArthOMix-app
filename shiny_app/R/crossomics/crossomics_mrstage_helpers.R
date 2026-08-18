## R/crossomics/crossomics_mrstage_helpers.R
## Pure data-processing logic for the "Cross-Omics MR" Cross-Omics
## sub-module (mod_cross_mr_stage.R) - a live, reconfigurable reimplementation
## of cross_Omics_Sexstratified_COPY/scripts/02_mr_stage_cross_omics.R:
## single-instrument (Wald-ratio) Mendelian Randomization for the eQTL-MR-
## significant genes, GoDMC cis-mQTL exposure -> Ishigaki et al. 2022 RA GWAS
## outcome, via TwoSampleMR.
##
## Scope boundary, stated here and in the UI: the instrument panel is the
## pipeline's own PRE-EXTRACTED lead-SNP set (one lead cis-SNP per CpG,
## already filtered from the raw ~6.3GB GoDMC association file upstream) -
## this module never reads that raw file live. "Reconfigurable" means the
## F-statistic/harmonisation/Steiger-handling parameters change on this
## fixed instrument panel, not that arbitrary genes can be added.
##
## MR estimates are genuine causal estimates under the standard IV
## assumptions (relevance, independence, exclusion restriction) - worded as
## such throughout, never claimed as proven causation, and explicitly
## flagged as impossible to validate via heterogeneity with only one
## instrument per exposure.

CX_MR_INSTRUMENT_FILE <- file.path(CX_DATA_ROOT, "data", "eqtl_sig_genes_lead_snps.csv")
CX_MR_RSID_LOOKUP_FILE <- file.path(CX_DATA_ROOT, "data", "lead_snp_rsid_lookup.csv")
CX_MR_OUTCOME_FILE <- file.path(CX_BC_Q3_DIR, "data", "raw", "ishigaki_ra_gwas", "GCST90132223_buildGRCh37.tsv.gz")
CX_MR_OUTCOME_ID <- "GCST90132223"
CX_MR_OUTCOME_NAME <- "Rheumatoid arthritis (Ishigaki et al. 2022, GCST90132223)"
CX_MR_RA_NCASE <- 22350L
CX_MR_RA_NCONTROL <- 74823L

CX_MR_DATA_AVAILABLE <- file.exists(CX_MR_INSTRUMENT_FILE) && file.exists(CX_MR_RSID_LOOKUP_FILE) && file.exists(CX_MR_OUTCOME_FILE) &&
  requireNamespace("TwoSampleMR", quietly = TRUE)

## ---------------------------------------------------------------------------
## Instrument panel (spec: GoDMC lead cis-SNP per eQTL-MR-significant gene's
## candidate CpG, pre-extracted upstream - see header)
## ---------------------------------------------------------------------------

## Returns list(ok, df, error, n_before_fstat, n_after_fstat). `df` has one
## row per CpG-instrument with a resolvable rsID and F-stat >= min_f_stat.
cx_mr_load_instruments <- function(min_f_stat = 10) {
  if (!CX_MR_DATA_AVAILABLE) {
    return(list(ok = FALSE, df = NULL, error = "Cross-Omics MR source data (GoDMC lead-SNP instrument panel and/or the Ishigaki RA GWAS outcome file) is not available in this deployment.", n_before_fstat = 0, n_after_fstat = 0))
  }
  lead <- tryCatch(as.data.frame(data.table::fread(CX_MR_INSTRUMENT_FILE, showProgress = FALSE)), error = function(e) e)
  if (inherits(lead, "error")) return(list(ok = FALSE, df = NULL, error = paste("Could not read the instrument panel:", conditionMessage(lead)), n_before_fstat = 0, n_after_fstat = 0))
  rsid <- tryCatch(as.data.frame(data.table::fread(CX_MR_RSID_LOOKUP_FILE, showProgress = FALSE)), error = function(e) e)
  if (inherits(rsid, "error")) return(list(ok = FALSE, df = NULL, error = paste("Could not read the rsID lookup:", conditionMessage(rsid)), n_before_fstat = 0, n_after_fstat = 0))

  m <- merge(lead, rsid[, c("name", "rsid"), drop = FALSE], by.x = "snp", by.y = "name")
  ## Same rule as the upstream pipeline's own instrument-prep step: an
  ## unresolvable (non-"rs"-prefixed) rsID can't be matched to the outcome
  ## GWAS file, so those instruments are dropped rather than guessed at.
  m <- m[grepl("^rs", m$rsid), , drop = FALSE]
  n_before <- nrow(m)
  m$F_stat <- (m$beta_a1 / m$se)^2
  m <- m[!is.na(m$F_stat) & m$F_stat >= min_f_stat, , drop = FALSE]
  n_after <- nrow(m)
  if (n_after == 0) return(list(ok = FALSE, df = NULL, error = sprintf("No instruments pass the F-statistic threshold (F >= %s).", min_f_stat), n_before_fstat = n_before, n_after_fstat = 0))
  list(ok = TRUE, df = m, error = NULL, n_before_fstat = n_before, n_after_fstat = n_after)
}

## ---------------------------------------------------------------------------
## Wald-ratio MR (spec: single instrument per CpG, TwoSampleMR::mr_wald_ratio,
## exactly mirroring script 02's harmonisation/Steiger/FDR choices)
## ---------------------------------------------------------------------------

CX_MR_DEFAULT_PARAMS <- list(min_f_stat = 10, harmonise_action = 2, steiger_mode = "flag")

## `instruments_df` is cx_mr_load_instruments()$df. Returns
## list(ok, df, error, n_tested, n_harmonised). `df` has one row per
## CpG-instrument: cpg, gene, SNP, b, se, pval, OR, OR_lo, OR_hi, FDR,
## mr_significant, steiger_dir, steiger_pval.
cx_mr_run_wald_ratio <- function(instruments_df, params = CX_MR_DEFAULT_PARAMS) {
  p <- utils::modifyList(CX_MR_DEFAULT_PARAMS, params %||% list())
  steiger_mode <- if (identical(p$steiger_mode, "drop")) "drop" else "flag"

  ## One exposure per CpG (phenotype_col/id_col = "cpg"), NOT per gene -
  ## confirmed by hand that keying on gene instead reproduces a real
  ## harmonise_data() crash from duplicate phenotype names (a gene can have
  ## more than one candidate CpG instrument); this matches script 02's own
  ## design, which assumes and asserts exactly one instrument per exposure.
  exp_dat <- tryCatch(TwoSampleMR::format_data(
    instruments_df, type = "exposure", phenotype_col = "cpg", snp_col = "rsid",
    beta_col = "beta_a1", se_col = "se", effect_allele_col = "allele1", other_allele_col = "allele2",
    eaf_col = "freq_a1", pval_col = "pval", samplesize_col = "samplesize", id_col = "cpg"
  ), error = function(e) e)
  if (inherits(exp_dat, "error")) return(list(ok = FALSE, df = NULL, error = paste("Could not format exposure (GoDMC instrument) data:", conditionMessage(exp_dat)), n_tested = 0, n_harmonised = 0))
  n_tested <- nrow(exp_dat)

  ## Reads the ~324MB outcome file FILTERED to just these instrument SNPs
  ## (TwoSampleMR::read_outcome_data()'s own snps= argument) - never the
  ## whole file into memory unfiltered. This step alone takes ~60-90s
  ## (confirmed by hand) - the caller should show progress/expect a wait.
  out_dat <- tryCatch(TwoSampleMR::read_outcome_data(
    snps = exp_dat$SNP, filename = CX_MR_OUTCOME_FILE, sep = "\t",
    snp_col = "variant_id", chr_col = "chromosome", pos_col = "base_pair_location",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    beta_col = "beta", se_col = "standard_error", pval_col = "p_value"
  ), error = function(e) e)
  if (inherits(out_dat, "error")) return(list(ok = FALSE, df = NULL, error = paste("Could not read the outcome (Ishigaki RA GWAS) data:", conditionMessage(out_dat)), n_tested = n_tested, n_harmonised = 0))
  if (nrow(out_dat) == 0) return(list(ok = FALSE, df = NULL, error = "None of the instrument SNPs were found in the outcome GWAS file.", n_tested = n_tested, n_harmonised = 0))
  out_dat$outcome <- CX_MR_OUTCOME_NAME
  out_dat$id.outcome <- CX_MR_OUTCOME_ID
  out_dat$samplesize.outcome <- CX_MR_RA_NCASE + CX_MR_RA_NCONTROL

  ## action = 2 (default, matches script 02): palindromic SNPs are dropped
  ## rather than allele-frequency-inferred, because the outcome file has no
  ## per-SNP EAF to infer strand from.
  harm <- tryCatch(TwoSampleMR::harmonise_data(exp_dat, out_dat, action = p$harmonise_action), error = function(e) e)
  if (inherits(harm, "error")) return(list(ok = FALSE, df = NULL, error = paste("Harmonisation failed:", conditionMessage(harm)), n_tested = n_tested, n_harmonised = 0))
  harm <- harm[harm$mr_keep %in% TRUE, , drop = FALSE]
  n_harmonised <- nrow(harm)
  if (n_harmonised == 0) return(list(ok = FALSE, df = NULL, error = "No instrument-outcome pairs survived harmonisation (all failed allele/strand matching).", n_tested = n_tested, n_harmonised = 0))

  ## Steiger filtering FLAGS (never silently drops) instruments where the
  ## outcome explains more variance than the exposure - i.e. potentially
  ## outcome-driven rather than exposure-driven - matching script 02's own
  ## "retained, flagged not dropped" design. Dropping is offered as an
  ## explicit, separate opt-in (steiger_mode = "drop").
  harm_sg <- tryCatch(TwoSampleMR::steiger_filtering(harm), error = function(e) e)
  steiger_ok <- !inherits(harm_sg, "error")
  if (steiger_ok) harm <- harm_sg

  res <- tryCatch(TwoSampleMR::mr(harm, method_list = "mr_wald_ratio"), error = function(e) e)
  if (inherits(res, "error")) return(list(ok = FALSE, df = NULL, error = paste("MR computation failed:", conditionMessage(res)), n_tested = n_tested, n_harmonised = n_harmonised))
  if (nrow(res) == 0) return(list(ok = FALSE, df = NULL, error = "No Wald-ratio estimates could be computed from the harmonised data.", n_tested = n_tested, n_harmonised = n_harmonised))

  res$OR <- exp(res$b)
  res$OR_lo <- exp(res$b - 1.96 * res$se)
  res$OR_hi <- exp(res$b + 1.96 * res$se)
  res$FDR <- stats::p.adjust(res$pval, method = "BH")
  res$mr_significant <- !is.na(res$FDR) & res$FDR < 0.05
  ## NOT res$id.exposure - confirmed by hand that TwoSampleMR::format_data()
  ## ignores id_col and auto-generates a random id.exposure regardless; the
  ## `exposure` column (from phenotype_col) is what actually carries the CpG
  ## text through format_data() -> harmonise_data() -> mr() unchanged.
  res$cpg <- res$exposure

  gene_lookup <- unique(instruments_df[, c("cpg", "gene"), drop = FALSE])
  res <- merge(res, gene_lookup, by = "cpg", all.x = TRUE)
  if (steiger_ok && all(c("exposure", "steiger_dir", "steiger_pval") %in% colnames(harm))) {
    sd_lookup <- unique(harm[, c("exposure", "steiger_dir", "steiger_pval"), drop = FALSE])
    res <- merge(res, sd_lookup, by.x = "cpg", by.y = "exposure", all.x = TRUE)
  } else {
    res$steiger_dir <- NA; res$steiger_pval <- NA_real_
  }
  if (identical(steiger_mode, "drop") && steiger_ok) {
    res <- res[res$steiger_dir %in% TRUE, , drop = FALSE]
    if (nrow(res) == 0) return(list(ok = FALSE, df = NULL, error = "No instruments remain after dropping Steiger-failing instruments.", n_tested = n_tested, n_harmonised = n_harmonised))
  }

  list(ok = TRUE, df = res, error = NULL, n_tested = n_tested, n_harmonised = n_harmonised, steiger_ok = steiger_ok)
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

CX_MR_TIER_RULE_TEXT <- c(
  "\"Credible\" mQTL-MR evidence = BH-FDR significant (FDR < 0.05, computed over this run's tested instruments) AND the gene is not MHC-flagged (per the eQTL-MR panel's own MHC_gene flag) AND Steiger filtering supports the methylation -> disease direction.",
  "Tier 1: DEG significant AND (DMP genome-wide-significant OR DMR significant) AND credible mQTL-MR evidence.",
  "Tier 2: DEG significant AND exactly one of - (DMP or DMR significant), credible mQTL-MR evidence.",
  "Tier 3: every other gene in the join (including genes with DEG evidence only, or with methylation evidence but no qualifying mQTL-MR support).",
  "Source: cross_Omics_Sexstratified_COPY/results/CROSS_OMICS_REPORT.md, sections 4.2.6-4.2.7. Implemented here in code for the first time - neither pipeline script computes Tier; this module reproduces the documented rule exactly, on the current run's own thresholds."
)

## `mr_df` is cx_mr_run_wald_ratio()$df (needs gene, FDR, mr_significant,
## steiger_dir). `join_df` is cx_bc_build_join()$df for the SAME sex (needs
## gene, DEG_significant, methylation_significant, eQTL_MHC_region) - if
## NULL, DEG/methylation evidence is reported "Not available" rather than
## guessed, and every gene falls back to Tier 3 (the only tier reachable
## without that evidence).
cx_mr_classify_tier <- function(mr_df, join_df = NULL) {
  mr_best <- cx_bc_dedup_min(mr_df, "gene", "pval")  ## best (lowest-p) instrument per gene, matching script 02's own per-gene summary
  mr_best$credible_mQTL_MR <- mr_best$mr_significant %in% TRUE & !(mr_best$steiger_dir %in% FALSE)

  if (is.null(join_df)) {
    out <- data.frame(gene = mr_best$gene, DEG_significant = NA, methylation_significant = NA,
                       credible_mQTL_MR = mr_best$credible_mQTL_MR, tier = "Tier 3",
                       tier_evidence_available = FALSE, stringsAsFactors = FALSE)
    return(out)
  }
  eqtl_mhc <- unique(join_df[, c("gene", "eQTL_MHC_region"), drop = FALSE])
  mr_best <- merge(mr_best, eqtl_mhc, by = "gene", all.x = TRUE)
  mr_best$credible_mQTL_MR <- mr_best$credible_mQTL_MR & !(mr_best$eQTL_MHC_region %in% TRUE)

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

cx_mr_build_provenance <- function(params, n_tested, n_harmonised, run_at) {
  c(
    "Exposure: GoDMC (Min et al. 2021) cis-mQTL summary statistics, pipeline's own pre-extracted lead-SNP-per-CpG panel (data/eqtl_sig_genes_lead_snps.csv) - this module never reads the raw genome-wide GoDMC association file live.",
    sprintf("Outcome: Rheumatoid arthritis GWAS, Ishigaki et al. 2022 (%s) - N cases = %s, N controls = %s.", CX_MR_OUTCOME_ID, format(CX_MR_RA_NCASE, big.mark = ","), format(CX_MR_RA_NCONTROL, big.mark = ",")),
    "Method: single-instrument Wald ratio (TwoSampleMR::mr_wald_ratio) - one lead cis-SNP per CpG; not IVW/MR-Egger, which need >=2 instruments per exposure.",
    sprintf("Minimum instrument F-statistic: %s", params$min_f_stat),
    sprintf("Harmonisation action: %s (2 = drop palindromic SNPs rather than infer strand from allele frequency, since the outcome file has no per-SNP EAF)", params$harmonise_action),
    sprintf("Steiger filtering: %s", if (identical(params$steiger_mode, "drop")) "instruments failing the methylation->disease direction test are dropped" else "instruments are flagged, not dropped (matches the original script)"),
    sprintf("Multiple-testing adjustment: Benjamini-Hochberg FDR, computed over the %s instruments actually tested in this run", format(n_harmonised, big.mark = ",")),
    sprintf("Instruments tested: %s -> harmonised/kept: %s", format(n_tested, big.mark = ","), format(n_harmonised, big.mark = ",")),
    "These are Mendelian randomization estimates, valid under the standard instrumental-variable assumptions (relevance, independence, exclusion restriction). A single instrument per exposure cannot be tested for validity via heterogeneity - interpret individual estimates cautiously, as an association consistent with a causal effect, not as proof of one.",
    sprintf("Run at: %s", run_at %||% "(not run yet)")
  )
}
