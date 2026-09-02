# 10_Diagnostic_Model

`TX_MODULES` stage `id = "diagnostic"`.

- **File**: `mod_diagnostic.R`.
- **Input**: `results$featureselection` (selected gene panel).
- **Main analysis**: logistic / elastic-net / RF / SVM diagnostic classifier panel, by sex.
- **Output**: `results$diagnostic` (trained model artifacts, performance metrics).
- **UI**: Transcriptomics → Sub-modules → Diagnostic Model.
