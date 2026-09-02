# 01_Data

The Cross-Omics "Dataset" tab — the data-loading entry point for the whole vertical. Not one of the `CX_MODULES` analytical stages; every stage below reads from the shared reactiveValues this tab publishes.

- **Files**: `mod_cross_dataset.R` (Shiny module trio, `id = "dataset"`), `crossomics_integration_upload.R` (`cx_read_table()`/`cx_read_and_detect()` — table parsing/type detection for custom uploads).
- **Input**: preloaded transcriptomics DEG + methylomics DMP/DMR result tables, or user-uploaded equivalents; sample-level methylation matrix for sample matching.
- **Main operation**: loads/validates the DEG and DMR/DMP tables, detects sample/feature ID columns (`cx_detect_sample_columns`), standardizes gene/CpG identifiers (`cx_standardize_*`, defined in `../functions/integration/crossomics_integration_helpers.R`).
- **Output**: the shared `cross_dataset` reactiveValues (standardized DEG table, DMR/DMP table, sample metadata) consumed by every stage in this vertical.
- **UI**: `ui.R`'s `crossomicsUI()` → "Dataset" tab. Mounted directly in `server.R` (`mod_cross_dataset_server("cx_dataset", cross_dataset)`), outside the `CX_MODULES` loop.
- **Dependencies**: `../functions/integration/crossomics_integration_helpers.R` (`cx_load_default_deg/dmr/methylation`, `cx_standardize_*`), `crossomics_integration_upload.R` (same folder).
