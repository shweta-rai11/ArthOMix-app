## R/crossomics/functions/biomarker_convergence/crossomics_biomarkerconv_helpers.R
## Data-processing logic for the "Biomarker Convergence" sub-module.

CX_BC_DATA_AVAILABLE <- CX_DATA_AVAILABLE

cx_bc_precomputed_file <- function(sex) file.path(CX_RESULTS_DIR, sprintf("cross_omics_eQTL_mQTL_%s.csv", sex))

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

## A gene can have up to several dozen CpG-level mQTL-MR tests in the MR-stage
## source file (one per nearby CpG instrument). Naively adopting the minimum
## p-value across those tests as "the" gene-level p-value is an uncorrected
## multiple-comparisons/winner's-curse selection: it systematically understates
## the true p-value and can flip a gene from "no instrument evidence" to
## "significant" purely because many CpGs, not one, were tested.
##
## Correction (two-stage, deterministic, no permutation required):
##   1. Within-gene: Bonferroni-correct the selected minimum p-value by the
##      number of CpG-instrument tests actually run for that gene
##      (adjusted_p = min(1, p_min * n_tests)). Bonferroni is used rather than a
##      permutation/minP approach because it is valid regardless of dependence
##      between nearby CpGs (which are frequently correlated), deterministic,
##      and directly answers "what is the chance of seeing a p-value this small
##      by chance, having tested n_tests instruments for this gene".
##   2. Across-gene: this corrected value is then eligible for BH-FDR across the
##      full displayed gene family via cx_bc_relabel(mqtl_sig_basis = "fdr")
##      (now the default - see CX_BC_DEFAULT_PARAMS), giving a standard
##      hierarchical (within-gene family-wise, then across-gene FDR) correction.
## The raw, uncorrected minimum p-value and the instrument count are retained in
## mQTL_MR_pval_raw_min / mQTL_instruments_tested for transparency - they are
## informational only and are never used for significance calls.
cx_bc_backfill_mqtl_from_mrstage <- function(df) {
  if (!"mQTL_MR_pval_raw_min" %in% colnames(df)) df$mQTL_MR_pval_raw_min <- NA_real_
  if (!"mQTL_instruments_tested" %in% colnames(df)) df$mQTL_instruments_tested <- NA_integer_

  if (!exists("CX_MR_DATA_AVAILABLE", inherits = TRUE) || !isTRUE(CX_MR_DATA_AVAILABLE)) return(df)
  mr <- tryCatch(as.data.frame(data.table::fread(CX_MR_PRECOMPUTED_FILE, showProgress = FALSE)), error = function(e) NULL)
  if (is.null(mr) || nrow(mr) == 0) return(df)

  n_tested <- table(mr$gene)
  mr_best <- cx_bc_dedup_min(mr, "gene", "pval")
  mr_best$n_instruments_tested <- as.integer(n_tested[mr_best$gene])
  mr_best$mQTL_MR_pval_gene_corrected <- pmin(1, mr_best$pval * mr_best$n_instruments_tested)

  needs_backfill <- !(df$in_mQTL_MR_panel %in% TRUE) & df$gene %in% mr_best$gene
  if (!any(needs_backfill)) return(df)

  m <- mr_best[match(df$gene[needs_backfill], mr_best$gene), , drop = FALSE]
  df$in_mQTL_MR_panel[needs_backfill] <- TRUE
  df$mQTL_candidate_cpg[needs_backfill] <- m$cpg
  df$mQTL_MR_beta[needs_backfill] <- m$b
  df$mQTL_MR_pval_raw_min[needs_backfill] <- m$pval
  df$mQTL_instruments_tested[needs_backfill] <- m$n_instruments_tested
  df$mQTL_MR_pval[needs_backfill] <- m$mQTL_MR_pval_gene_corrected
  df$mQTL_instrument_available[needs_backfill] <- TRUE
  df
}

CX_BC_REQUIRED_EQTL_COLS <- "gene"
CX_BC_REQUIRED_MQTL_COLS <- c("gene", "mQTL_MR_pval")

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

cx_bc_dedup_min <- function(df, key_col, order_col) {
  if (!any(duplicated(df[[key_col]]))) return(df)
  by_key <- split(seq_len(nrow(df)), df[[key_col]])
  idx <- vapply(by_key, function(ix) {
    ov <- df[[order_col]][ix]
    ix[if (all(is.na(ov))) 1L else which.min(ifelse(is.na(ov), Inf, ov))]
  }, integer(1))
  df[idx, , drop = FALSE]
}

## Fixed thresholds (no UI): the precomputed join was built at these cut-offs.
## mqtl_sig_basis defaults to "fdr" (not the raw nominal p) because mQTL_MR_pval
## for backfilled genes is itself already a within-gene Bonferroni-corrected
## value (see cx_bc_backfill_mqtl_from_mrstage()); applying BH-FDR across the
## full displayed gene family on top of that gives the required second stage of
## correction (within-gene family-wise, then across-gene FDR).
CX_BC_DEFAULT_PARAMS <- list(
  deg_fdr = 0.05, dmp_genomewide_fdr = 0.05,
  mqtl_sig_basis = "fdr", mqtl_sig_cutoff = 0.05, dmr_fdr = 0.05,
  eqtl_fdr = 0.05
)

CX_BC_THRESHOLD_TEXT <- paste0(
  "Fixed thresholds: DEG FDR < 0.05; DMP (bacon) FDR < 0.05; DMR FDR < 0.05. ",
  "mQTL-MR: genes with more than one CpG instrument are first Bonferroni-corrected within-gene ",
  "(p_gene = min(1, min(p) x n_instruments_tested), to correct the minimum-p-value selection across up to dozens ",
  "of CpGs per gene); the resulting p-values are then BH-FDR-corrected across the full displayed gene set, FDR < 0.05. ",
  "eQTL-MR: uses the table's own significance flag if present, else FDR < 0.05 if supplied, else panel membership."
)

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

  ## eQTL-MR: use the pipeline's own flag if present, else FDR threshold, else panel membership.
  out$eQTL_MR_significant <- if ("eQTL_MR_significant" %in% colnames(df)) {
    out$in_eQTL_MR_panel %in% TRUE & as.logical(df$eQTL_MR_significant) %in% TRUE
  } else if ("eQTL_MR_FDR" %in% colnames(out) && any(!is.na(out$eQTL_MR_FDR))) {
    out$in_eQTL_MR_panel %in% TRUE & !is.na(out$eQTL_MR_FDR) & out$eQTL_MR_FDR < p$eqtl_fdr
  } else {
    out$in_eQTL_MR_panel %in% TRUE
  }

  out$methylation_significant <- out$DMP_genomewide_significant %in% TRUE | out$DMR_significant %in% TRUE
  out$n_evidence_layers <- rowSums(cbind(out$eQTL_MR_significant %in% TRUE, out$DEG_significant %in% TRUE,
                                          out$methylation_significant %in% TRUE, out$mQTL_MR_significant %in% TRUE))
  out
}

