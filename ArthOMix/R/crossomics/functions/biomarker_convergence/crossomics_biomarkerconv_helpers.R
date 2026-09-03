## R/crossomics/functions/biomarker_convergence/crossomics_biomarkerconv_helpers.R
## Pure data-processing logic for the "Biomarker Convergence" Cross-Omics
## sub-module (mod_cross_biomarker_conv.R).

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

CX_BC_DEFAULT_PARAMS <- list(
  deg_fdr = 0.05, dmp_genomewide_fdr = 0.05,
  mqtl_sig_basis = "nominal_p", mqtl_sig_cutoff = 0.05, dmr_fdr = 0.05
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

  out$eQTL_MR_significant <- out$in_eQTL_MR_panel %in% TRUE

  out$methylation_significant <- out$DMP_genomewide_significant %in% TRUE | out$DMR_significant %in% TRUE
  out$n_evidence_layers <- rowSums(cbind(out$eQTL_MR_significant %in% TRUE, out$DEG_significant %in% TRUE,
                                          out$methylation_significant %in% TRUE, out$mQTL_MR_significant %in% TRUE))
  out
}

