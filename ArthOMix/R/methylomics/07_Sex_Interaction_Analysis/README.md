# 07_Sex_Interaction_Analysis

`MX_MODULES` stage `id = "interaction"`. Mirrors `R/transcriptomics/mod_interaction.R`'s group×sex interaction design, applied to methylation data.

- **File**: `mod_methyl_interaction.R`.
- **Input**: normalized `methyl_dataset`.
- **Main analysis**: live limma `group*sex` interaction model per CpG, using `methyl_chunked_lmfit()` and `mod_methyl_dmp_sex_col()` from `../05_Differential_Methylation_Position/mod_methyl_dmp.R`.
- **Output**: `results$interaction`.
- **UI**: Methylomics → Sub-modules → Sex Interaction Analysis.
- **Dependencies**: `../functions/qc.R` (`methyl_sheet_sample_ids`); cross-module helpers from `mod_methyl_dmp.R`.
