# 09_Feature_Selection

`TX_MODULES` stage `id = "featureselection"`. Counterpart to `R/methylomics/10_ML_Feature_Selection/mod_methyl_featureselection.R`.

- **File**: `mod_featureselection.R`.
- **Input**: candidate genes.
- **Main analysis**: LASSO / random forest / SVM-RFE feature selection, by sex.
- **Output**: `results$featureselection` (selected gene panel).
- **UI**: Transcriptomics → Sub-modules → Feature Selection.
