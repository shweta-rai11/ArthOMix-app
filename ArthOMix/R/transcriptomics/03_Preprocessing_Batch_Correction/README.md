# 03_Preprocessing_Batch_Correction

`TX_MODULES` stage `id = "preprocessing"`.

- **Files**: `mod_preprocessing.R` (also defines the reusable `mod_pp_source_ui/_server` per-dataset source picker used elsewhere in this module), `mod_preprocessing_explore.R` (the "Data Exploration" sub-tab's EDA statistics/plot functions — `mod_data_exploration_ui/_server`, mounted inside `mod_preprocessing.R` but not itself a `TX_MODULES` entry).
- **Input**: `dataset` (raw, possibly multi-source).
- **Main analysis**: multi-dataset merge, batch correction, normalization, exploratory data analysis (skewness/kurtosis/PCA/correlation/dendrogram/QQ).
- **Output**: `results$preprocessing`; updates `dataset` in place with the processed expression matrix.
- **UI**: Transcriptomics → Sub-modules → Preprocessing and Batch Correction.
