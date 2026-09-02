# 11_Sex_Interaction_Analysis

`TX_MODULES` stage `id = "interaction"`. Mirrored by `R/methylomics/07_Sex_Interaction_Analysis/mod_methyl_interaction.R`.

- **File**: `mod_interaction.R`.
- **Input**: processed `dataset`.
- **Main analysis**: live limma `group*sex` interaction model.
- **Output**: `results$interaction` (mirrors `mod_dge.R`'s `results$dge` contract).
- **UI**: Transcriptomics → Sub-modules → Sex Interaction Analysis.
