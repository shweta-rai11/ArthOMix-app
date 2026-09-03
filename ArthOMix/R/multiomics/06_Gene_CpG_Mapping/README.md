# 06_Gene_CpG_Mapping

`MULTI_MODULES` stage `id = "mapping"`, title **"Gene–CpG Mapping"**.

- **Files**: `mod_multi_mapping.R`, `multiomics_mapping_helpers.R` (`mcc_*`), `multiomics_mapping_plots.R`.
- **Input**: `multi_dataset`'s matched expression + methylation layers, optionally `multi_results$biomarker`.
- **Main analysis**: per gene-CpG pair expression/methylation correlation and evidence classification.
- **Output**: `multi_results$mapping` — the one table joining a candidate gene's expression evidence to its CpG's methylation evidence, consumed by `../08_Biomarker_Card/`.
- **UI**: Multiomics → Sub-modules → Gene–CpG Mapping.
- **Note**: `multiomics_mapping_plots.R`'s network plot reuses the same technique as `R/crossomics/functions/integration/crossomics_integration_plots.R::cx_gene_cpg_network_plot()`.
