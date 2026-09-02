# 04_SNF_Clustering

`MULTI_MODULES` stage `id = "stratification"`, title **"SNF Clustering"** (the task template guessed "Patient Stratification" — this is the real title from `R/submodules_registry.R`).

- **Files**: `mod_multi_stratification.R`, `snf_clustering_helpers.R` (`sfc_*`), `snf_clustering_plots.R`.
- **Input**: `multi_results$integration`'s fused SNF network from `../03_DIABLO_SNF_Integration/`.
- **Main analysis**: consensus clustering on the fused network to identify patient subgroups.
- **Output**: `multi_results$stratification` (cluster assignments, stability metrics).
- **UI**: Multiomics → Sub-modules → SNF Clustering.
- **Dependencies**: `snf_clustering_helpers.R`/`_plots.R` (same folder), which delegate to `../functions/multiomics_integration_helpers.R`/`_plots.R` (`mi_snf_*`, `mi_ari`) for every ≥2-block case — only the single-omics (1-block) fallback is unique to this folder's helpers. See `../README.md` and `../../REFACTORING_NOTES.md`.
