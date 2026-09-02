# 09_Results_Summary

`MULTI_MODULES` stage `id = "summary"`, title **"Results Summary & Reproducibility"**.

- **File**: `mod_multi_summary.R`.
- **Input**: `multi_dataset`, `multi_results` (all upstream stages).
- **Main analysis**: none — a reproducibility/summary roll-up, not a statistical stage.
- **Output**: package version report (`multi_package_versions`), analysis summary table (`multi_analysis_summary_table`) across every stage run this session.
- **UI**: Multiomics → Sub-modules → Results Summary & Reproducibility.
- **Dependencies**: `../functions/multiomics_helpers.R`.
