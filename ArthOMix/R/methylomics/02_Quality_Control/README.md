# 02_Quality_Control

`MX_MODULES` stage `id = "qc"`.

- **File**: `mod_methyl_qc.R`.
- **Input**: the shared `methyl_dataset` reactiveValues from `../01_Data/`.
- **Main analysis**: sample/probe-level QC — detection p-values, intensity outliers, sex-check, batch diagnostics.
- **Output**: `results$qc` (QC-passed sample/probe sets, flags).
- **UI**: Methylomics → Sub-modules → Quality Control.
- **Dependencies**: `../functions/qc.R`, `../functions/idat_metrics.R`, `../functions/annotation.R`, `../functions/parse_upload.R`.
