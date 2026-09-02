# 06_Gene_CpG_Concordance

`MULTI_MODULES` stage `id = "concordance"`, title **"Gene–CpG Concordance"**.

- **Files**: `mod_multi_concordance.R`, `multiomics_concordance_helpers.R` (`mcc_*`), `multiomics_concordance_plots.R`.
- **Input**: `multi_dataset`'s matched expression + methylation layers, optionally `multi_results$biomarker`.
- **Main analysis**: per gene-CpG pair expression/methylation correlation and evidence classification.
- **Output**: `multi_results$concordance` — the one table joining a candidate gene's expression evidence to its CpG's methylation evidence, consumed by `../08_Biomarker_Card/`.
- **UI**: Multiomics → Sub-modules → Gene–CpG Concordance.
- **Note**: `multiomics_concordance_plots.R`'s network plot reuses the same technique as `R/crossomics/functions/integration/crossomics_integration_plots.R::cx_gene_cpg_network_plot()`.
