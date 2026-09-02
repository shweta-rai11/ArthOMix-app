# 08_Colocalization

`TX_MODULES` stage `id = "coloc"`. Counterpart to `R/methylomics/12_Colocalization/mod_methyl_coloc.R` (eQTL colocalization here vs. mQTL there).

- **File**: `mod_coloc.R`.
- **Input**: candidate genes + eQTL/GWAS summary statistics.
- **Main analysis**: Bayesian colocalization at candidate loci.
- **Output**: `results$coloc`.
- **UI**: Transcriptomics → Sub-modules → Colocalization.
