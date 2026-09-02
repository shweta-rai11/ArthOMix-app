# 05_Biomarker_Discovery

`MULTI_MODULES` stage `id = "biomarker"`, title **"Biomarker Discovery"**.

- **Files**: `mod_multi_biomarker.R`, `multiomics_biomarker_helpers.R` (`mb_*`), `multiomics_biomarker_plots.R`.
- **Input**: `multi_dataset`, optionally `multi_results$stratification`.
- **Main analysis**: joint cross-layer biomarker discovery, sex-stratified (via the shared nested-CV engine).
- **Output**: `multi_results$biomarker`.
- **UI**: Multiomics → Sub-modules → Biomarker Discovery.
- **Dependencies**: `../functions/multiomics_sexstratified_engine.R` (`mss_*`).
