# 03_Normalisation

`MX_MODULES` stage `id = "normalization"`.

- **File**: `mod_methyl_normalization.R`. Runs normalization methods in a `callr::r_bg()` background worker (`methyl_norm_bg_worker()`/`methyl_norm_compare_bg_worker()`) — these workers `source()` `../functions/annotation.R`, `../functions/qc.R`, `../functions/normalization.R`, and this file itself by explicit path (a fresh `callr` process has none of the app's already-loaded files), so any future move of those four files must update the `for (f in c(...))` loop near the top of `methyl_norm_bg_worker()`.
- **Input**: `methyl_dataset` (post-QC).
- **Main analysis**: noob / funnorm / SWAN / dasen / BMIQ / PBC / quantile normalization, with diagnostics and status/recommendation.
- **Output**: `results$normalization`; normalized beta/M-value matrix written back to `methyl_dataset`.
- **UI**: Methylomics → Sub-modules → Normalisation.
- **Dependencies**: `../functions/normalization.R` (method implementations), `../functions/qc.R`, `../functions/annotation.R`, `../functions/parse_upload.R`.
