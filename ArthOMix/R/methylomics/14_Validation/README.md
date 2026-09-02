# 14_Validation

`MX_MODULES` stage `id = "validation"`.

- **File**: `mod_methyl_validation.R`.
- **Input**: `results$diagnostic` (trained model), independent validation cohort (upload or preloaded).
- **Main analysis**: validates the Diagnostic Classifier's trained models on an independent cohort.
- **Output**: `results$validation` (external validation performance metrics).
- **UI**: Methylomics → Sub-modules → Validation.
- **Dependencies**: `../functions/qc.R`, `../functions/parse_upload.R`.
