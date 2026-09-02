# 10_ML_Feature_Selection

`MX_MODULES` stage `id = "featureselection"`.

- **File**: `mod_methyl_featureselection.R`. Calls `mod_methyl_dmp_sex_col/_choices` from `../05_Differential_Methylation_Position/mod_methyl_dmp.R`.
- **Input**: `results$candidates` (or DMP/DMR results directly).
- **Main analysis**: LASSO / random forest / SVM-RFE feature selection over candidate CpGs, by sex.
- **Output**: `results$featureselection` (selected CpG panel).
- **UI**: Methylomics → Sub-modules → ML Feature Selection.
- **Dependencies**: `../functions/qc.R`, `../functions/annotation.R`, `../functions/parse_upload.R`.
