# 07_Pathways

`MULTI_MODULES` stage `id = "pathway"`, title **"Pathways"**.

- **Files**: `mod_multi_pathway.R`, `multiomics_pathway_helpers.R` (`mp_*`, including `mp_get_wikipathways_termgene()` — also called cross-vertical from `R/methylomics/15_Biomarker_Analysis/mod_methyl_biomarkercard.R` and `R/transcriptomics/mod_biomarkercard.R`), `multiomics_pathway_plots.R`.
- **Input**: `multi_results$mapping` or `multi_results$biomarker`.
- **Main analysis**: pathway-level integration/enrichment over discovered cross-omics biomarkers.
- **Output**: `multi_results$pathway`.
- **UI**: Multiomics → Sub-modules → Pathways.
