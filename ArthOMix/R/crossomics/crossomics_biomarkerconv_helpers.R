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
## columns - never re-joining, never re-deriving evidence.
##
## One known defect in that per-sex file itself: in_eQTL_MR_panel and
## in_mQTL_MR_panel are mutually exclusive for every gene (0 overlap),
## because it was built by stacking a separately-exported eQTL-MR block and
## an mQTL-MR block without merging genes present in both - confirmed
## against the sibling MASTER_cross_omics_all_layers.csv (same folder),
## which retains both eQTL_MR_OR and mQTL_MR_cpg per gene and shows most
## "eQTL-only" genes (30/32 female, 23/25 male) genuinely also have mQTL-MR
## results. cx_bc_load_precomputed() backfills exactly that gap - see
## cx_bc_backfill_mqtl_from_mrstage() below - from
## mr_stage_eqtl_significant_genes_mqtl_mr.csv (crossomics_mrstage_helpers.R's
## own precomputed source, already restricted to eQTL-MR-significant genes),
## which is why this file is the one exception to "never reaching into any
## file outside this table": it's the pipeline's own sibling precomputed
## result, read here strictly read-only (no new statistics computed), and it
## only fills gaps the source table left NA/FALSE - it never overwrites a
## gene the table already correctly flagged.
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
## that already happened upstream and is exactly what's on disk), then
## passed through cx_bc_backfill_mqtl_from_mrstage() to correct the
## in_eQTL_MR_panel/in_mQTL_MR_panel data defect documented at the top of
## this file - a gap-fill, not a re-join.
cx_bc_load_precomputed <- function(sex = c("female", "male", "combined")) {
  sex <- match.arg(sex)
  if (!CX_BC_DATA_AVAILABLE) return(list(ok = FALSE, df = NULL, error = "Biomarker Convergence source data is not available in this deployment."))
  path <- cx_bc_precomputed_file(sex)
  if (!file.exists(path)) return(list(ok = FALSE, df = NULL, error = sprintf("Precomputed eQTL x mQTL table for %s is not available (%s).", sex, basename(path))))
  df <- tryCatch(as.data.frame(data.table::fread(path, showProgress = FALSE)), error = function(e) e)
  if (inherits(df, "error")) return(list(ok = FALSE, df = NULL, error = paste("Could not read the precomputed join:", conditionMessage(df))))
  if (nrow(df) == 0) return(list(ok = FALSE, df = NULL, error = "The precomputed join has no rows."))
  df <- cx_bc_backfill_mqtl_from_mrstage(df)
  list(ok = TRUE, df = df, error = NULL)
}

## Fixes the in_eQTL_MR_panel/in_mQTL_MR_panel data-prep defect documented at
## the top of this file. `df` is cross_omics_eQTL_mQTL_{sex}.csv, read as-is.
## For every gene in `df` that is NOT already flagged in_mQTL_MR_panel = TRUE
## but has a real mQTL-MR instrument result in
## mr_stage_eqtl_significant_genes_mqtl_mr.csv (CX_MR_PRECOMPUTED_FILE,
## crossomics_mrstage_helpers.R's own precomputed source - one row per
## CpG-instrument, already restricted to eQTL-MR-significant genes, which is
## exactly the population this join dropped mQTL evidence for), this fills
## in_mQTL_MR_panel/mQTL_candidate_cpg/mQTL_MR_beta/mQTL_MR_pval/
## mQTL_instrument_available from that gene's best (lowest-p) CpG instrument
## (cx_bc_dedup_min() - the same per-gene reduction crossomics_mrstage_
## helpers.R's own Tier classification already uses). mQTL_cpg_chr/
## mQTL_cpg_pos_hg19 are left NA for backfilled rows - that file doesn't
## carry CpG coordinates, and this never fabricates a value the source
## doesn't have. mQTL_MR_significant is left as-is here; cx_bc_relabel()
## recomputes it for every row unconditionally from mQTL_MR_pval, so the
## newly-filled value is picked up automatically on the next relabel.
## Fail-soft: if the MR source isn't available, `df` is returned unchanged.
cx_bc_backfill_mqtl_from_mrstage <- function(df) {
  if (!exists("CX_MR_DATA_AVAILABLE", inherits = TRUE) || !isTRUE(CX_MR_DATA_AVAILABLE)) return(df)
  mr <- tryCatch(as.data.frame(data.table::fread(CX_MR_PRECOMPUTED_FILE, showProgress = FALSE)), error = function(e) NULL)
  if (is.null(mr) || nrow(mr) == 0) return(df)
  mr_best <- cx_bc_dedup_min(mr, "gene", "pval")

  needs_backfill <- !(df$in_mQTL_MR_panel %in% TRUE) & df$gene %in% mr_best$gene
  if (!any(needs_backfill)) return(df)

  m <- mr_best[match(df$gene[needs_backfill], mr_best$gene), , drop = FALSE]
  df$in_mQTL_MR_panel[needs_backfill] <- TRUE
  df$mQTL_candidate_cpg[needs_backfill] <- m$cpg
  df$mQTL_MR_beta[needs_backfill] <- m$b
  df$mQTL_MR_pval[needs_backfill] <- m$pval
  df$mQTL_instrument_available[needs_backfill] <- TRUE
  df
}

## ---------------------------------------------------------------------------
## "Upload your own data" - eQTL-MR and mQTL-MR results as two SEPARATE
## files (each optional independently), merged live by gene into the same
## one-row-per-gene shape cx_bc_relabel() expects - a plain outer join, not
## a re-run of either MR analysis. Neither file needs DEG/DMP/DMR columns
## (this module never re-derives those; see the Dataset tab / Expression x
## Methylation module for that) - cx_bc_merge_eqtl_mqtl() below adds them as
## all-NA so DEG/methylation significance simply reads as "not evaluated"
## for uploaded data, never fabricated.
## ---------------------------------------------------------------------------

CX_BC_REQUIRED_EQTL_COLS <- "gene"
CX_BC_REQUIRED_MQTL_COLS <- c("gene", "mQTL_MR_pval")

## Same list(ok, df, error) contract as cx_bc_load_precomputed(). Reuses
## cx_read_table() (crossomics_integration_upload.R) for CSV/TSV/TXT/XLSX
## parsing - the same reader "Upload your own data" uses on the Dataset tab.
## Every gene present in this file is treated as in_eQTL_MR_panel = TRUE by
## cx_bc_merge_eqtl_mqtl() - only `gene` is required; OR/p/FDR/MHC columns
## are carried through if present, shown as NA if not.
cx_bc_load_eqtl_upload <- function(datapath, filename) {
  res <- cx_read_table(datapath, filename)
  if (!res$ok) return(list(ok = FALSE, df = NULL, error = res$error))
  df <- res$df
  missing <- setdiff(CX_BC_REQUIRED_EQTL_COLS, colnames(df))
  if (length(missing) > 0) {
    return(list(ok = FALSE, df = NULL, error = sprintf(
      "This eQTL-MR file is missing required column(s): %s. One row per gene; only a \"gene\" column is required.",
      paste(missing, collapse = ", "))))
  }
  if (nrow(df) == 0) return(list(ok = FALSE, df = NULL, error = "The uploaded eQTL-MR file has no rows."))
  df$gene <- as.character(df$gene)
  df <- df[!is.na(df$gene) & nzchar(df$gene), , drop = FALSE]
  if (nrow(df) == 0) return(list(ok = FALSE, df = NULL, error = "No row in the uploaded eQTL-MR file has a non-empty gene value."))
  for (cl in intersect(c("eQTL_MR_OR", "eQTL_MR_pval", "eQTL_MR_FDR"), colnames(df))) {
    df[[cl]] <- suppressWarnings(as.numeric(df[[cl]]))
  }
  list(ok = TRUE, df = df, error = NULL)
}

## Same contract, for the mQTL-MR file - here `mQTL_MR_pval` is also
## required (it's what cx_bc_relabel() needs to compute significance for
## this layer at all; every other mQTL_* column is optional).
cx_bc_load_mqtl_upload <- function(datapath, filename) {
  res <- cx_read_table(datapath, filename)
  if (!res$ok) return(list(ok = FALSE, df = NULL, error = res$error))
  df <- res$df
  missing <- setdiff(CX_BC_REQUIRED_MQTL_COLS, colnames(df))
  if (length(missing) > 0) {
    return(list(ok = FALSE, df = NULL, error = sprintf(
      "This mQTL-MR file is missing required column(s): %s. One row per gene; \"gene\" and \"mQTL_MR_pval\" are required.",
      paste(missing, collapse = ", "))))
  }
  if (nrow(df) == 0) return(list(ok = FALSE, df = NULL, error = "The uploaded mQTL-MR file has no rows."))
  df$gene <- as.character(df$gene)
  df <- df[!is.na(df$gene) & nzchar(df$gene), , drop = FALSE]
  if (nrow(df) == 0) return(list(ok = FALSE, df = NULL, error = "No row in the uploaded mQTL-MR file has a non-empty gene value."))
  df$mQTL_MR_pval <- suppressWarnings(as.numeric(df$mQTL_MR_pval))
  for (cl in intersect(c("mQTL_MR_beta", "mQTL_cpg_pos_hg19"), colnames(df))) {
    df[[cl]] <- suppressWarnings(as.numeric(df[[cl]]))
  }
  list(ok = TRUE, df = df, error = NULL)
}

## Merges the (at most) two separately-uploaded layers by gene - a full
## outer join, purely mechanical (never re-runs either MR analysis). At
## least one of eqtl_df/mqtl_df must be supplied. in_eQTL_MR_panel/
## in_mQTL_MR_panel are set from which file(s) a gene actually appears in;
## DEG_adjP/DMP_fdr_bacon/DMR_fdr and (if the mQTL-MR file was omitted)
## mQTL_MR_pval are added as all-NA so cx_bc_relabel() still has every
## column it reads, computing "not significant/not evaluated" for evidence
## this upload path genuinely doesn't have - never invented.
cx_bc_merge_eqtl_mqtl <- function(eqtl_df, mqtl_df) {
  if (is.null(eqtl_df) && is.null(mqtl_df)) {
    return(list(ok = FALSE, df = NULL, error = "Upload at least one file (eQTL-MR and/or mQTL-MR)."))
  }
  genes <- unique(c(if (!is.null(eqtl_df)) eqtl_df$gene, if (!is.null(mqtl_df)) mqtl_df$gene))
  out <- data.frame(gene = genes, stringsAsFactors = FALSE)
  out$in_eQTL_MR_panel <- if (!is.null(eqtl_df)) out$gene %in% eqtl_df$gene else FALSE
  out$in_mQTL_MR_panel <- if (!is.null(mqtl_df)) out$gene %in% mqtl_df$gene else FALSE
  if (!is.null(eqtl_df)) out <- merge(out, eqtl_df, by = "gene", all.x = TRUE)
  if (!is.null(mqtl_df)) out <- merge(out, mqtl_df, by = "gene", all.x = TRUE)
  for (cl in c("DEG_adjP", "DMP_fdr_bacon", "DMR_fdr")) if (!cl %in% colnames(out)) out[[cl]] <- NA_real_
  if (!"mQTL_MR_pval" %in% colnames(out)) out$mQTL_MR_pval <- NA_real_
  list(ok = TRUE, df = out, error = NULL)
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

