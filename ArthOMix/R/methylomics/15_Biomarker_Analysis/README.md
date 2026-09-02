# 15_Biomarker_Analysis

`MX_MODULES` stage `id = "biomarkercard"`. Largest file in the vertical.

- **File**: `mod_methyl_biomarkercard.R`.
- **Input**: a selected CpG/gene (from any upstream stage's results, or manual search).
- **Main analysis**: single-CpG/panel biomarker lookup — live annotation, evidence classification, cross-omics context.
- **Output**: downloadable biomarker report; no `results$*` write-back (read-only interpretation layer).
- **UI**: Methylomics → Sub-modules → Biomarker Analysis.
- **Dependencies**: `../functions/annotation.R`, `../functions/parse_upload.R`; cross-vertical `cx_harmonize_gene_ids()` (`R/crossomics/functions/integration/crossomics_integration_helpers.R`) and `mp_get_wikipathways_termgene()` (`R/multiomics/multiomics_pathway_helpers.R`).
