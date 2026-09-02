# 13_Diagnostic_Classifier

`MX_MODULES` stage `id = "diagnostic"`.

- **File**: `mod_methyl_diagnostic.R`.
- **Input**: `results$featureselection` (selected CpG panel).
- **Main analysis**: logistic / elastic-net / RF / SVM diagnostic classifier panel, by sex.
- **Output**: `results$diagnostic` (trained model artifacts, performance metrics).
- **UI**: Methylomics → Sub-modules → Diagnostic Classifier.
- **Dependencies**: `../functions/qc.R`, `../functions/parse_upload.R`.
