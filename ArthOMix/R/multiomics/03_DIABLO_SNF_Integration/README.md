# 03_DIABLO_SNF_Integration

`MULTI_MODULES` stage `id = "integration"`, title **"Multi-omics Integration (DIABLO & SNF)"**.

- **File**: `mod_multi_integration.R`.
- **Input**: `multi_dataset` (matched samples across layers).
- **Main analysis**: DIABLO discriminant analysis across omics layers, and SNF (Similarity Network Fusion) to build a fused patient similarity network.
- **Output**: `multi_results$integration` (DIABLO loadings/performance, fused SNF network) — the fused network feeds `../04_SNF_Clustering/`.
- **UI**: Multiomics → Sub-modules → Multi-omics Integration (DIABLO & SNF).
- **Dependencies**: `../functions/multiomics_integration_helpers.R` (`mi_diablo_*`, `mi_snf_*`), `../functions/multiomics_integration_plots.R`, `../functions/multiomics_plots.R`, `../functions/multiomics_sexstratified_engine.R` (`mss_run_stratified`).
