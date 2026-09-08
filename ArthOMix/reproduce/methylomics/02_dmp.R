## Regenerates the methylomics DMP (SVA-adjusted) stage from beta_raw.rds +
## pheno.rds, using the app's own live functions (mod_methyl_dmp_prepare_subset,
## mod_methyl_sva_fit, methyl_chunked_lmfit, limma, bacon), and compares
## against the precomputed dmp_{female,male}_full.csv tables.
##
## Model, per the app's own UI description (mod_methyl_dmp.R:450): sex-stratified
## limma on M-values (group + age + smoking + cell-type), bacon-corrected,
## BH-FDR. Cell-type fractions are estimated live via EpiDISH CP on the
## blood7 (Reinius/Houseman) reference - the module's own UI defaults
## (mod_methyl_celltype.R:636,652). The live UI's true default covariate
## selection is EMPTY (checkboxGroupInput selected=character(0)) - the
## "group + age + smoking + cell-type" set is this script's explicit choice,
## matching the documented default-table methodology, not a UI default.
##
## Run from the ArthOMix/ app directory:
##   Rscript reproduce/methylomics/02_dmp.R

suppressMessages(suppressWarnings(
  shiny::loadSupport(".", renv = globalenv(), globalrenv = globalenv())
))

beta <- readRDS("data/preloaded/methylomics/matrix/beta_raw.rds")
pheno <- as.data.frame(readRDS("data/preloaded/methylomics/matrix/pheno.rds"))
rownames(pheno) <- pheno$gsm
colnames(beta) <- pheno$gsm[match(colnames(beta), pheno$gsm)]
cat(sprintf("beta_raw.rds: %d probes x %d samples\n", nrow(beta), ncol(beta)))
cat("group levels:", paste(table(pheno$group), collapse = " / "), "\n")

## --- Cell-type deconvolution (live EpiDISH CP, blood7 reference) ---
ref_mat <- methyl_ct_get_reference("blood7")
stopifnot(!is.null(ref_mat))
ct <- methyl_ct_run_epidish(beta, ref_mat, method = "CP")
stopifnot(isTRUE(ct$ok))
fractions <- ct$fractions
cat(sprintf("Cell-type fractions: %d samples x %d types (%s)\n", nrow(fractions), ncol(fractions), paste(colnames(fractions), collapse = ",")))

run_dmp_sex <- function(sex_choice) {
  methyl_dataset <- list(beta = beta, sample_sheet = pheno, input_scale = "beta")
  ct_choices <- mod_methyl_dmp_celltype_choices(fractions)
  covariate_cols <- c("age", "smoking", unname(ct_choices))

  sub <- mod_methyl_dmp_prepare_subset(
    methyl_dataset, sex_choice = sex_choice, sex_col_name = "sex",
    group_col = "group", ref = "Control", comp = "RA",
    covariate_cols = covariate_cols,
    min_valid_pct = 80, min_variance = 0, snp_filter = FALSE,
    anno = methyl_get_annotation("450K"), celltype_fractions = fractions
  )
  m <- sub$m; beta_scale <- sub$beta_scale; grp <- sub$grp; cov_df <- sub$cov_df
  cat(sprintf("  [%s] n_ref(Control)=%d n_comp(RA)=%d, %d probes after filters\n", sex_choice, sub$n_ref, sub$n_comp, nrow(m)))

  sv_fit <- mod_methyl_sva_fit(m, grp, cov_df)
  design <- sv_fit$design
  stopifnot(qr(design)$rank == ncol(design))
  cat(sprintf("  [%s] %d surrogate variables retained\n", sex_choice, sv_fit$n_sv))

  fit <- methyl_chunked_lmfit(m, design)
  cm <- limma::makeContrasts(contrasts = "RA-Control", levels = design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm))
  tt <- limma::topTable(fit2, number = Inf, sort.by = "none")

  bc <- bacon::bacon(teststatistics = tt$t, verbose = FALSE)
  p_bacon <- as.numeric(bacon::pval(bc))
  fdr_bacon <- stats::p.adjust(p_bacon, method = "BH")

  beta_ref <- rowMeans(beta_scale[, grp == "Control", drop = FALSE], na.rm = TRUE)
  beta_comp <- rowMeans(beta_scale[, grp == "RA", drop = FALSE], na.rm = TRUE)
  dbeta <- (beta_comp - beta_ref)[rownames(tt)]

  data.frame(cpg = rownames(tt), logFC_M_regen = tt$logFC, t_regen = tt$t,
             p_raw_regen = tt$P.Value, p_bacon_regen = p_bacon,
             fdr_bacon_regen = fdr_bacon, dbeta_regen = dbeta,
             row.names = NULL, stringsAsFactors = FALSE)
}

compare_sex <- function(sex_choice, precomp_path) {
  cat(sprintf("\n--- %s ---\n", sex_choice))
  regen <- run_dmp_sex(sex_choice)
  precomp <- read.csv(precomp_path)
  cmp <- merge(precomp, regen, by = "cpg")
  cat(sprintf("Matched %d / %d precomputed CpGs\n", nrow(cmp), nrow(precomp)))

  r_t <- cor(cmp$t, cmp$t_regen, use = "complete.obs")
  r_dbeta <- cor(cmp$dbeta, cmp$dbeta_regen, use = "complete.obs")
  n_sig_precomp <- sum(cmp$fdr_bacon < 0.05, na.rm = TRUE)
  n_sig_regen <- sum(cmp$fdr_bacon_regen < 0.05, na.rm = TRUE)
  sig_overlap <- sum(cmp$fdr_bacon < 0.05 & cmp$fdr_bacon_regen < 0.05, na.rm = TRUE)

  cat(sprintf("t-statistic correlation:   r = %.4f\n", r_t))
  cat(sprintf("dbeta correlation:         r = %.4f\n", r_dbeta))
  cat(sprintf("Significant CpGs (FDR<0.05): precomputed=%d, regenerated=%d, overlap=%d\n",
              n_sig_precomp, n_sig_regen, sig_overlap))

  write.csv(cmp, sprintf("reproduce/methylomics/output/02_dmp_%s_comparison.csv", sex_choice), row.names = FALSE)
  list(sex = sex_choice, r_t = r_t, r_dbeta = r_dbeta, n_sig_precomp = n_sig_precomp,
       n_sig_regen = n_sig_regen, sig_overlap = sig_overlap, n_matched = nrow(cmp))
}

res_f <- compare_sex("F", "data/preloaded/methylomics/tables/script03_dmp_sva_sexstratified/tables/dmp_female_full.csv")
res_m <- compare_sex("M", "data/preloaded/methylomics/tables/script03_dmp_sva_sexstratified/tables/dmp_male_full.csv")

saveRDS(list(female = res_f, male = res_m), "reproduce/methylomics/output/02_dmp_summary.rds")
cat("\n=== DMP regeneration summary ===\n")
print(rbind(as.data.frame(res_f), as.data.frame(res_m)))
