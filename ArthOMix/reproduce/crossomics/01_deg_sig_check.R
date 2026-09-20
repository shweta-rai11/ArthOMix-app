## Checks the precomputed DEG_sig (MASTER_cross_omics_all_layers.csv) against significance recomputed from 01_deg.R.
## Run from the app directory after 01_deg.R: Rscript reproduce/crossomics/01_deg_sig_check.R

deg_f <- read.csv("data/preloaded/transcriptomics/results/tables/DEG_female_full.csv")
deg_m <- read.csv("data/preloaded/transcriptomics/results/tables/DEG_male_full.csv")
cx <- read.csv("data/preloaded/cross_omics/tables/MASTER_cross_omics_all_layers.csv")

compute_sig <- function(deg_df) {
  ## DEG_sig uses DEG_full.csv's own `sig` column (|logFC|>=0.1), not the UI's default |logFC|>=0.5.
  data.frame(gene = deg_df$gene, sig_regen = as.logical(deg_df$sig))
}

check_sex <- function(sex_label, deg_df) {
  cat(sprintf("\n--- %s ---\n", sex_label))
  sig_regen <- compute_sig(deg_df)
  cx_sex <- cx[tolower(cx$sex) == tolower(sex_label), c("gene", "DEG_sig")]
  cx_sex$DEG_sig_bool <- as.logical(cx_sex$DEG_sig)
  cmp <- merge(cx_sex, sig_regen, by = "gene")
  cat(sprintf("Matched %d / %d precomputed Cross-Omics rows\n", nrow(cmp), nrow(cx_sex)))
  agree <- mean(cmp$DEG_sig_bool == cmp$sig_regen, na.rm = TRUE)
  n_disagree <- sum(cmp$DEG_sig_bool != cmp$sig_regen, na.rm = TRUE)
  cat(sprintf("DEG_sig agreement: %.2f%% (%d/%d), %d disagreements\n",
              100 * agree, sum(cmp$DEG_sig_bool == cmp$sig_regen, na.rm = TRUE), nrow(cmp), n_disagree))
  if (n_disagree > 0 && n_disagree <= 10) {
    print(cmp[cmp$DEG_sig_bool != cmp$sig_regen, c("gene", "DEG_sig_bool", "sig_regen")])
  }
  list(sex = sex_label, n_matched = nrow(cmp), agree_pct = 100 * agree, n_disagree = n_disagree)
}

res_f <- check_sex("female", deg_f)
res_m <- check_sex("male", deg_m)
cat("\n=== Cross-Omics DEG_sig reproduction summary ===\n")
print(rbind(as.data.frame(res_f), as.data.frame(res_m)))
