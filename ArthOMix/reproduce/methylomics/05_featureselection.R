## Regenerates the methylomics ML Feature Selection stage (LASSO stability
## selection + Random Forest importance) using the module's own function-level
## defaults - see reproduce/methylomics/README.md. LOWER CONFIDENCE than
## 02_dmp.R/03_dmr.R/04_wgcna.R: this module has NO preloaded-specific
## hyperparameter defaults anywhere in mod_methyl_featureselection.R (checked;
## zero `if (isTRUE(...preloaded...))` branches for any ML hyperparameter), so
## there is no code evidence these generic defaults are what the original
## pipeline actually used. Reported honestly as an approximation, not a proof.
##
## Recipe replicated from mod_methyl_featureselection.R:786-865 (fs_build_filters):
## 70/30 holdout (seed 1234, BEFORE any feature filtering - leakage-safe, same
## as the live app), then on the training partition only: 5%/20% probe/sample
## missingness filter, top-5000-by-variance cap (beta scale). LASSO stability
## (methyl_fs_stability_run defaults: bootstrap, 50 resamples, alpha=1, seed
## 1234) and RF (methyl_fs_rf_fit defaults: ntree=1000, mtry=sqrt(p), nodesize=1,
## replace=TRUE, seed=1234, gini importance).
##
## Run from the ArthOMix/ app directory:
##   Rscript reproduce/methylomics/05_featureselection.R

suppressMessages(suppressWarnings(
  shiny::loadSupport(".", renv = globalenv(), globalrenv = globalenv())
))

beta <- readRDS("data/preloaded/methylomics/matrix/beta_raw.rds")
pheno <- as.data.frame(readRDS("data/preloaded/methylomics/matrix/pheno.rds"))
rownames(pheno) <- pheno$gsm

run_fs_sex <- function(sex_letter, sex_label, lasso_freq_path, rf_imp_path) {
  cat(sprintf("\n--- %s ---\n", sex_label))
  sub <- methyl_qc_subgroup_filter(beta, pheno, "sex", sex_letter, min_n = 10)
  mat0 <- sub$mat
  ph0 <- pheno[colnames(mat0), ]
  grp_raw <- ph0$group
  keep_s <- grp_raw %in% c("Control", "RA")
  mat1 <- mat0[, keep_s, drop = FALSE]
  grp_full <- factor(grp_raw[keep_s], levels = c("Control", "RA"))

  set.seed(1234)
  train_idx <- as.integer(caret::createDataPartition(grp_full, p = 0.7, list = FALSE)[, 1])
  mat1 <- mat1[, train_idx, drop = FALSE]
  grp <- grp_full[train_idx]
  cat(sprintf("Training partition (70%%, held-out excluded): %d samples (%s)\n", ncol(mat1), paste(table(grp), collapse = "/")))

  beta_scale <- mat1
  f_probe <- methyl_filter_missing(beta_scale, 0.05)
  beta_scale <- beta_scale[f_probe$keep, , drop = FALSE]
  beta_scale <- methyl_fs_cap_top_variance(beta_scale, 5000)
  cat(sprintf("Candidate CpGs after top-5000-variance cap: %d\n", nrow(beta_scale)))

  X <- t(beta_scale)

  cat("Running LASSO stability selection (50 bootstrap resamples)...\n")
  stab <- methyl_fs_stability_run(X, grp, type = "bootstrap", n_resamples = 50, seed = 1234, alpha = 1)
  freq_regen <- data.frame(cpg = stab$ranked$cpg, freq_regen = stab$ranked$selection_frequency, row.names = NULL)

  lasso_precomp <- read.csv(lasso_freq_path)
  cmp_lasso <- merge(lasso_precomp, freq_regen, by = "cpg", all.x = TRUE)
  cmp_lasso$freq_regen[is.na(cmp_lasso$freq_regen)] <- 0
  r_lasso <- cor(cmp_lasso$freq, cmp_lasso$freq_regen)
  cat(sprintf("LASSO selection-frequency correlation (precomputed vs regenerated, n=%d): r = %.4f\n", nrow(cmp_lasso), r_lasso))

  cat("Running Random Forest (ntree=1000)...\n")
  rf <- methyl_fs_rf_fit(X, grp, ntree = 1000, seed = 1234)
  imp_regen <- data.frame(cpg = names(rf$importance), importance_regen = as.numeric(rf$importance), row.names = NULL)

  rf_precomp <- read.csv(rf_imp_path)
  imp_col <- intersect(c("MeanDecreaseGini", "MeanDecreaseAccuracy"), colnames(rf_precomp))[1]
  cmp_rf <- merge(rf_precomp, imp_regen, by = "cpg")
  ## rf_precomp has far fewer rows than this script's 5000-CpG top-variance
  ## candidate pool (e.g. 12 for the female stratum, close to the DMP
  ## significant-CpG count) - the real Feature Selection candidate universe was
  ## almost certainly the DMP-significant panel, not a fresh top-variance
  ## selection. Report the mismatch honestly rather than force a number.
  if (nrow(cmp_rf) < 3) {
    cat(sprintf("RF importance: NOT COMPARABLE - only %d/%d precomputed CpGs fall inside this script's 5000-CpG candidate pool (precomputed table looks DMP-panel-derived, not a fresh top-variance selection - see header note).\n", nrow(cmp_rf), nrow(rf_precomp)))
    r_rf <- NA_real_; top50_overlap <- NA_integer_
  } else {
    r_rf <- cor(cmp_rf[[imp_col]], cmp_rf$importance_regen)
    cat(sprintf("RF importance correlation (precomputed vs regenerated, matched=%d/%d): r = %.4f\n", nrow(cmp_rf), nrow(rf_precomp), r_rf))
    top_n <- min(50, nrow(rf_precomp))
    top50_precomp <- head(rf_precomp$cpg[order(-rf_precomp[[imp_col]])], top_n)
    top50_regen <- head(names(sort(rf$importance, decreasing = TRUE)), top_n)
    top50_overlap <- length(intersect(top50_precomp, top50_regen))
    cat(sprintf("RF top-%d CpG overlap: %d / %d\n", top_n, top50_overlap, top_n))
  }

  write.csv(cmp_lasso, sprintf("reproduce/methylomics/output/05_fs_%s_lasso_comparison.csv", sex_letter), row.names = FALSE)
  write.csv(cmp_rf, sprintf("reproduce/methylomics/output/05_fs_%s_rf_comparison.csv", sex_letter), row.names = FALSE)
  list(sex = sex_label, n_train = ncol(mat1), n_candidates = nrow(beta_scale),
       r_lasso_freq = r_lasso, r_rf_importance = r_rf, rf_top50_overlap = top50_overlap)
}

res_f <- run_fs_sex("F", "Female",
  "data/preloaded/methylomics/tables/script07_ml_feature_selection/tables/lasso_selection_frequency_female.csv",
  "data/preloaded/methylomics/tables/script07_ml_feature_selection/tables/rf_importance_female.csv")
res_m <- run_fs_sex("M", "Male",
  "data/preloaded/methylomics/tables/script07_ml_feature_selection/tables/lasso_selection_frequency_male.csv",
  "data/preloaded/methylomics/tables/script07_ml_feature_selection/tables/rf_importance_male.csv")

saveRDS(list(female = res_f, male = res_m), "reproduce/methylomics/output/05_fs_summary.rds")
cat("\n=== Feature Selection regeneration summary (LOWER CONFIDENCE - see header) ===\n")
print(rbind(as.data.frame(res_f), as.data.frame(res_m)))
