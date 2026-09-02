# 12_Cross_Tissue_Validation

`TX_MODULES` stage `id = "crosstissue"`.

- **File**: `mod_crosstissue.R`.
- **Input**: `results$diagnostic` (trained model artifacts).
- **Main analysis**: cross-tissue validation of the diagnostic model, four-classifier panel, by sex.
- **Output**: `results$crosstissue`.
- **UI**: Transcriptomics → Sub-modules → Cross-Tissue Validation.
