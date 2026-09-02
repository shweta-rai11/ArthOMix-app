# 14_Functional_Enrichment

`TX_MODULES` stage `id = "enrichment"`.

- **File**: `mod_enrichment.R`.
- **Input**: a gene list (from upload or upstream `results$*`).
- **Main analysis**: live GO/KEGG/Reactome over-representation analysis.
- **Output**: `results$enrichment`.
- **UI**: Transcriptomics → Sub-modules → Functional Enrichment.
