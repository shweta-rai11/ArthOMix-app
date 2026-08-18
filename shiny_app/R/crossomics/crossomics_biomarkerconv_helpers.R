## R/crossomics/crossomics_biomarkerconv_helpers.R
## Pure data-processing logic for the "Biomarker Convergence" Cross-Omics
## sub-module (mod_cross_biomarker_conv.R).
##
## The join itself (eQTL-MR gene panel x mQTL-MR CpG panel x raw DEG/DMP/DMR
## tables, per sex) was already run once by the pipeline's own
## cross_Omics_Sexstratified_COPY/scripts/01_cross_omics_eqtl_mqtl_biomarkers.R,
## and its output already sits in
## cross_Omics_Sexstratified_COPY/results/cross_omics_eQTL_mQTL_{female,male}.csv
## - every field this module needs (DEG_logFC/adjP, DMP_dbeta/fdr_bacon,
## DMR_meandiff/fdr, mQTL_MR_beta/pval, eQTL_MR_OR/FDR/MHC flag) is already
## in that one file. This module's job is to load that file (instant) and
## let the significance thresholds that decide "*_significant" be
## reconfigured live by RELABELING the retained raw p/FDR/effect-size
## columns - never re-joining, never re-deriving evidence, and never
## reaching into any file outside cross_Omics_Sexstratified_COPY.
##
## Every loader is fail-soft (list(ok, df, error)) - never errors, never
## fabricates.

## ---------------------------------------------------------------------------
## Precomputed join loader (CX_DATA_ROOT/CX_RESULTS_DIR/CX_DATA_AVAILABLE
## already defined in global.R:439-441, resolving to
## cross_Omics_Sexstratified_COPY itself)
## ---------------------------------------------------------------------------

CX_BC_DATA_AVAILABLE <- CX_DATA_AVAILABLE

cx_bc_precomputed_file <- function(sex) file.path(CX_RESULTS_DIR, sprintf("cross_omics_eQTL_mQTL_%s.csv", sex))

## Returns list(ok, df, error). `df` is the pipeline's own already-joined,
## one-row-per-gene table for `sex`, read as-is (no join/dedup logic here -
## that already happened upstream and is exactly what's on disk).
cx_bc_load_precomputed <- function(sex = c("female", "male")) {
  sex <- match.arg(sex)
  if (!CX_BC_DATA_AVAILABLE) return(list(ok = FALSE, df = NULL, error = "Biomarker Convergence source data is not available in this deployment."))
  path <- cx_bc_precomputed_file(sex)
  if (!file.exists(path)) return(list(ok = FALSE, df = NULL, error = sprintf("Precomputed eQTL x mQTL table for %s is not available (%s).", sex, basename(path))))
  df <- tryCatch(as.data.frame(data.table::fread(path, showProgress = FALSE)), error = function(e) e)
  if (inherits(df, "error")) return(list(ok = FALSE, df = NULL, error = paste("Could not read the precomputed join:", conditionMessage(df))))
  if (nrow(df) == 0) return(list(ok = FALSE, df = NULL, error = "The precomputed join has no rows."))
  list(ok = TRUE, df = df, error = NULL)
}

## Generic "keep one row per key, the one with the smallest order_col" dedup
## - still needed by crossomics_mrstage_helpers.R's Tier classification
## (best CpG-instrument per gene).
cx_bc_dedup_min <- function(df, key_col, order_col) {
  if (!any(duplicated(df[[key_col]]))) return(df)
  by_key <- split(seq_len(nrow(df)), df[[key_col]])
  idx <- vapply(by_key, function(ix) {
    ov <- df[[order_col]][ix]
    ix[if (all(is.na(ov))) 1L else which.min(ifelse(is.na(ov), Inf, ov))]
  }, integer(1))
  df[idx, , drop = FALSE]
}

## ---------------------------------------------------------------------------
## Relabeling (spec: every threshold the original script hardcoded is a
## live parameter here, applied to the ALREADY-JOINED table's own retained
## raw values - relabeling is cheap, so this runs on every filter change,
## no "Run" button needed)
## ---------------------------------------------------------------------------

## DMP_suggestive_only is intentionally not reconfigurable here: the
## precomputed table retains DMP_fdr_bacon (adjusted) but not the
## underlying raw p-value, so its pipeline-computed value is displayed
## as-is (see cx_bc_relabel()) rather than faking a threshold this data
## can't actually support - never invents a value the source doesn't have.
CX_BC_DEFAULT_PARAMS <- list(
  deg_fdr = 0.05, dmp_genomewide_fdr = 0.05,
  mqtl_sig_basis = "nominal_p", mqtl_sig_cutoff = 0.05, dmr_fdr = 0.05
)

## `df` is cx_bc_load_precomputed()$df. Returns the same table with
## *_significant/n_evidence_layers columns recomputed at `params` - no join,
## no new statistics, only relabeling from already-retained raw values.
cx_bc_relabel <- function(df, params = CX_BC_DEFAULT_PARAMS) {
  p <- utils::modifyList(CX_BC_DEFAULT_PARAMS, params %||% list())
  out <- df

  out$DEG_significant <- !is.na(out$DEG_adjP) & out$DEG_adjP < p$deg_fdr
  out$DMP_genomewide_significant <- !is.na(out$DMP_fdr_bacon) & out$DMP_fdr_bacon < p$dmp_genomewide_fdr
  out$DMR_significant <- !is.na(out$DMR_fdr) & out$DMR_fdr < p$dmr_fdr

  mqtl_stat <- if (identical(p$mqtl_sig_basis, "fdr")) {
    ok <- !is.na(out$mQTL_MR_pval)
    fdr <- rep(NA_real_, nrow(out))
    fdr[ok] <- stats::p.adjust(out$mQTL_MR_pval[ok], method = "BH")
    fdr
  } else out$mQTL_MR_pval
  out$mQTL_MR_significant <- !is.na(mqtl_stat) & mqtl_stat < p$mqtl_sig_cutoff

  ## eQTL-MR panel membership is already FDR<0.05-filtered upstream (that
  ## filtering happened before this table was even built) - "significant"
  ## here just means "is in the panel at all", not a separate threshold.
  out$eQTL_MR_significant <- out$in_eQTL_MR_panel %in% TRUE

  out$methylation_significant <- out$DMP_genomewide_significant %in% TRUE | out$DMR_significant %in% TRUE
  out$n_evidence_layers <- rowSums(cbind(out$eQTL_MR_significant %in% TRUE, out$DEG_significant %in% TRUE,
                                          out$methylation_significant %in% TRUE, out$mQTL_MR_significant %in% TRUE))
  out
}

cx_bc_build_provenance <- function(sex, params, run_at) {
  p <- utils::modifyList(CX_BC_DEFAULT_PARAMS, params %||% list())
  c(
    sprintf("Sex: %s", toupper(sex)),
    sprintf("Source: cross_Omics_Sexstratified_COPY/results/%s - the pipeline's own already-joined eQTL-MR x mQTL-MR x DEG x DMP x DMR table (no re-join performed here).", basename(cx_bc_precomputed_file(sex))),
    sprintf("DEG significance threshold: FDR < %s (relabeled from the table's own DEG_adjP)", p$deg_fdr),
    sprintf("DMP genome-wide-significant threshold: FDR < %s (relabeled from DMP_fdr_bacon)", p$dmp_genomewide_fdr),
    "DMP suggestive flag: shown as computed by the pipeline (raw P < 1e-4) - not reconfigurable here, since the underlying raw p-value isn't retained in this table.",
    sprintf("DMR significant threshold: FDR < %s (relabeled from DMR_fdr)", p$dmr_fdr),
    sprintf("mQTL-MR significance basis: %s < %s (relabeled from mQTL_MR_pval)", if (identical(p$mqtl_sig_basis, "fdr")) "BH-FDR (recomputed over this table's tested CpGs)" else "nominal P", p$mqtl_sig_cutoff),
    "No new statistical test is performed at this stage - only relabeling of already-computed per-layer results at the thresholds above.",
    sprintf("Run at: %s", run_at %||% "(not run yet)")
  )
}
