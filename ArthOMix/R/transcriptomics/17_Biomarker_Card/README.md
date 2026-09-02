# 17_Biomarker_Card

`TX_MODULES` stage `id = "biomarkercard"`. Largest file in the vertical.

- **File**: `mod_biomarkercard.R`.
- **Input**: a selected gene (from any upstream stage's results, or manual search).
- **Main analysis**: single-gene/panel biomarker lookup — live GO/KEGG/Reactome/WikiPathways/Open Targets/HPA/STRING/DGIdb/PubMed queries.
- **Output**: downloadable biomarker report; no `results$*` write-back.
- **UI**: Transcriptomics → Sub-modules → Biomarker Card.
- **Dependencies**: `R/crossomics/functions/integration/crossomics_integration_helpers.R` (evidence classification, shared with the Methylomics and Multiomics Biomarker Cards).
