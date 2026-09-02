# 03_Biomarker_Convergence

`CX_MODULES` stage `id = "biomarkerconv"`, titled **"Biomarker Convergence"** (group "Data") in the Cross-Omics Sub-modules grid.

- **File**: `mod_cross_biomarker_conv.R` (Shiny module trio).
- **Input**: `cross_results$integration` from `../02_Expression_Methylation_Integration/`.
- **Main analysis**: ranks/tiers candidate genes by strength of cross-omics evidence (expression + methylation direction agreement) using `cx_bc_*` helpers in `../functions/biomarker_convergence/crossomics_biomarkerconv_helpers.R`.
- **Output**: `cross_results$biomarkerconv` — ranked/tiered candidate table, reused as an input by `../04_Cross_Omics_MR/` (`cx_bc_load_precomputed`/`cx_bc_relabel`) and by the Transcriptomics/Multiomics Biomarker Card modules (`cx_classify_evidence()`, see `../README.md`).
- **UI**: Cross-Omics → Sub-modules → "Biomarker Convergence" card.
- **Dependencies**: `../functions/biomarker_convergence/crossomics_biomarkerconv_helpers.R`; `cx_empty_state` from `../functions/integration/crossomics_integration_plots.R`.
