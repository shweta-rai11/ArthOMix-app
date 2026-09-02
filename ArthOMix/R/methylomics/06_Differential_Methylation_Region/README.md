# 06_Differential_Methylation_Region

`MX_MODULES` stage `id = "dmr"`.

- **File**: `mod_methyl_dmr.R`. Calls shared helpers physically defined in `../05_Differential_Methylation_Position/mod_methyl_dmp.R` (`mod_methyl_dmp_sex_col/_choices/_covariate_cols/_filter/_volcano/_betadist`, `mod_methyl_lambda_gc`, `mod_methyl_qq_plot`) — safe across folders since every file under `R/methylomics/` is still sourced into one shared environment (`R/0_load_omics_modules.R`).
- **Input**: normalized `methyl_dataset` (+ optionally `results$dmp`).
- **Main analysis**: region-level differential methylation (comb-p/DMRcate-style region calling over per-CpG statistics).
- **Output**: `results$dmr` (per-region effect size, p-value, FDR).
- **UI**: Methylomics → Sub-modules → Differential Methylation Region.
- **Dependencies**: `../functions/qc.R`, `../functions/annotation.R`, `../functions/normalization.R`.
