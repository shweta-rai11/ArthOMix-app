# 01_Data

The Transcriptomics "Dataset" tab — data-loading entry point, not a `TX_MODULES` stage.

- **File**: `mod_dataset.R`.
- **Input**: upload (expression matrix + sample metadata), GEO fetch, or preloaded dataset.
- **Main operation**: loads/validates the expression matrix and sample metadata, matches sample IDs, rejects duplicate/blank IDs.
- **Output**: the shared `dataset` reactiveValues every stage below reads from.
- **UI**: Transcriptomics → "Dataset" tab. Mounted directly in `server.R` (`mod_dataset_server("tx_dataset", dataset)`), outside the `TX_MODULES` loop.
- **Dependencies**: `../functions/expression_type.R` (raw-vs-normalized detection).
