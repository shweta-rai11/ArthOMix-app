# 01_Data

The Methylomics "Dataset" tab — data-loading entry point, not an `MX_MODULES` stage.

- **File**: `mod_methyl_dataset.R`.
- **Input**: preloaded GSE42861 catalog selection, GEO fetch by GPL platform, or user upload (beta/M-value matrix + sample sheet, or raw IDAT pairs).
- **Main operation**: loads/validates the methylation matrix and sample metadata, detects array platform, harmonizes sample IDs.
- **Output**: the shared `methyl_dataset` reactiveValues every stage below reads from.
- **UI**: Methylomics → "Dataset" tab. Mounted directly in `server.R` (`mod_methyl_dataset_server("mx_dataset", methyl_dataset)`), outside the `MX_MODULES` loop.
- **Dependencies**: `../functions/idat_metrics.R`, `../functions/parse_upload.R`.
