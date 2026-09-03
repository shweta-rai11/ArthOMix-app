## R/crossomics/04_Cross_Omics_MR/crossomics_mrstage_helpers.R
## Pure data-processing logic for the "Cross-Omics MR" Cross-Omics
## sub-module (mod_cross_mr_stage.R).

CX_MR_PRECOMPUTED_FILE <- file.path(CX_RESULTS_DIR, "mr_stage_eqtl_significant_genes_mqtl_mr.csv")
CX_MR_DATA_AVAILABLE <- CX_DATA_AVAILABLE && file.exists(CX_MR_PRECOMPUTED_FILE)
CX_MR_OUTCOME_ID <- "GCST90132223"
CX_MR_OUTCOME_NAME <- "Rheumatoid arthritis (Ishigaki et al. 2022, GCST90132223)"

CX_MR_VOLCANO_LEVELS <- c("Up (risk, b > 0)", "Down (protective, b < 0)", "Not significant")

cx_mr_load_precomputed <- function() {
  if (!CX_MR_DATA_AVAILABLE) {
    return(list(ok = FALSE, df = NULL, error = "The pipeline's precomputed Cross-Omics MR results are not available in this deployment."))
  }
  df <- tryCatch(as.data.frame(data.table::fread(CX_MR_PRECOMPUTED_FILE, showProgress = FALSE)), error = function(e) e)
  if (inherits(df, "error")) return(list(ok = FALSE, df = NULL, error = paste("Could not read the precomputed MR results:", conditionMessage(df))))
  if (nrow(df) == 0) return(list(ok = FALSE, df = NULL, error = "The precomputed MR results file has no rows."))
  if (!is.logical(df$steiger_dir)) df$steiger_dir <- as.logical(df$steiger_dir)
  list(ok = TRUE, df = df, error = NULL)
}

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
  df$in_eQTL_MR_panel <- !is.na(df$eQTL_MR_FDR) & df$eQTL_MR_FDR < 0.05
  list(ok = TRUE, df = df, error = NULL)
}

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

