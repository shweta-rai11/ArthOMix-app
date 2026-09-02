# 12_Colocalization

`MX_MODULES` stage `id = "coloc"`.

- **File**: `mod_methyl_coloc.R` (self-contained).
- **Input**: candidate CpGs + mQTL/GWAS summary statistics.
- **Main analysis**: Bayesian colocalization (mQTL vs GWAS) at candidate loci.
- **Output**: `results$coloc`.
- **UI**: Methylomics → Sub-modules → Colocalization.
