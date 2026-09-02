# 08_WGCNA_Co_Methylation_Network

`MX_MODULES` stage `id = "wgcna"`.

- **File**: `mod_methyl_wgcna.R`. Calls `mod_methyl_dmp_sex_col/_choices/_covariate_cols` from `../05_Differential_Methylation_Position/mod_methyl_dmp.R`.
- **Input**: normalized `methyl_dataset`.
- **Main analysis**: WGCNA co-methylation network (soft-threshold, modules, module-trait correlation).
- **Output**: `results$wgcna`.
- **UI**: Methylomics → Sub-modules → WGCNA Co-Methylation Network.
- **Dependencies**: `../functions/qc.R`, `../functions/normalization.R` (`methyl_get_norm_annotation`).
