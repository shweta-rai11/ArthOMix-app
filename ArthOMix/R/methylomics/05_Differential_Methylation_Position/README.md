# 05_Differential_Methylation_Position

`MX_MODULES` stage `id = "dmp"`.

- **File**: `mod_methyl_dmp.R`. Also defines helpers reused directly by `../06_Differential_Methylation_Region/`, `../07_Sex_Interaction_Analysis/`, `../08_WGCNA_Co_Methylation_Network/`, and `../10_ML_Feature_Selection/` — see `../README.md`'s "Cross-module helper coupling" section before splitting or moving this file further.
- **Input**: normalized `methyl_dataset`.
- **Main analysis**: per-CpG differential methylation (limma-style linear model), with SVA/genomic-control correction (`mod_methyl_sva_fit`, `mod_methyl_lambda_gc`).
- **Output**: `results$dmp` (per-CpG effect size, p-value, FDR).
- **UI**: Methylomics → Sub-modules → Differential Methylation Position.
- **Dependencies**: `../functions/qc.R`, `../functions/annotation.R`, `../functions/normalization.R` (`methyl_norm_status`).
