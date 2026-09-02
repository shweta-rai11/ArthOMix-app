# functions

Shared helpers used across multiple Multiomics stages.

- **`multiomics_helpers.R`** — cross-cutting app-wide utilities (`multi_read_table`, `multi_active_dataset_banner`, `multi_qc_scorecard`, `multi_package_versions`, `multi_build_report`). Used by nearly every stage.
- **`multiomics_plots.R`** — shared plotting utilities (`multi_empty_state`, `multi_plot_or_empty`, `multi_diablo_*_plot`). Used across nearly every stage.
- **`multiomics_sexstratified_engine.R`** (`mss_*`) — shared sex-stratified nested-CV engine, used by `../03_DIABLO_SNF_Integration/` and `../05_Biomarker_Discovery/`.
- **`multiomics_integration_helpers.R`** / **`multiomics_integration_plots.R`** (`mi_*`) — DIABLO/SNF implementation, primarily for `../03_DIABLO_SNF_Integration/` but also called from `../04_SNF_Clustering/` (duplicate SNF logic — see `../README.md`).

See `../README.md` for the full per-stage dependency table.
