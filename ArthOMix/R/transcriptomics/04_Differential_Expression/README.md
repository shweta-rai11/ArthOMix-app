# 04_Differential_Expression

`TX_MODULES` stage `id = "dge"`.

- **File**: `mod_dge.R`. Sex-stratified analysis is handled inline (group/sex selectors), not a separate file.
- **Input**: processed `dataset` from `../03_Preprocessing_Batch_Correction/`.
- **Main analysis**: limma/DESeq2 differential expression.
- **Output**: `results$dge` (per-gene log2FC, p-value, FDR).
- **UI**: Transcriptomics → Sub-modules → Differential Expression.
- **Dependencies**: `../functions/expression_type.R` (decoupled "Upload your own data" path).
