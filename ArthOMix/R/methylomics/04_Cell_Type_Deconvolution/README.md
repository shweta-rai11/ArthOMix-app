# 04_Cell_Type_Deconvolution

`MX_MODULES` stage `id = "celltype"`.

- **File**: `mod_methyl_celltype.R` (includes the cell-type engine functions, `methyl_ct_*`, merged directly into this module).
- **Input**: normalized `methyl_dataset` beta matrix.
- **Main analysis**: estimates blood/tissue cell-type composition from methylation signal.
- **Output**: `results$celltype` (per-sample estimated cell-type proportions).
- **UI**: Methylomics → Sub-modules → Cell Type Deconvolution.
- **Dependencies**: `../functions/qc.R`, `../functions/annotation.R`, `../functions/parse_upload.R`.
