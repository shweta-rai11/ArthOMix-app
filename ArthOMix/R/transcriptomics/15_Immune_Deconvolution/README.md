# 15_Immune_Deconvolution

`TX_MODULES` stage `id = "deconvolution"`.

- **File**: `mod_deconvolution.R`.
- **Input**: `dataset$expr` (must pass a raw-counts run gate).
- **Main analysis**: immune cell composition estimation (CIBERSORT/LM22, MCP-counter).
- **Output**: `results$deconvolution`.
- **UI**: Transcriptomics → Sub-modules → Immune Deconvolution.
- **Dependencies**: `../functions/expression_type.R` (`looks_like_raw_counts()` run gate).
