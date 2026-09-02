# 07_Mendelian_Randomization

`TX_MODULES` stage `id = "mr"`. Counterpart to `R/methylomics/11_Mendelian_Randomization/mod_methyl_mr.R` (eQTL MR here vs. mQTL MR there).

- **File**: `mod_mr.R`.
- **Input**: candidate genes + GWAS/eQTL summary statistics (upload or default instrument table).
- **Main analysis**: two-sample Mendelian Randomization.
- **Output**: `results$mr`.
- **UI**: Transcriptomics → Sub-modules → Mendelian Randomization.
