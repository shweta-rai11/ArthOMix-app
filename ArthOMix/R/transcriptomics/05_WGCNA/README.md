# 05_WGCNA

`TX_MODULES` stage `id = "wgcna"`.

- **File**: `mod_wgcna.R`.
- **Input**: processed `dataset`.
- **Main analysis**: WGCNA co-expression network (soft-threshold selection, module detection, module-trait correlation, hub genes, Cytoscape export).
- **Output**: `results$wgcna`.
- **UI**: Transcriptomics → Sub-modules → WGCNA.
