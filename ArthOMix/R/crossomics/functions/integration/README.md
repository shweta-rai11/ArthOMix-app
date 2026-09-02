# functions/integration

Shared helpers/plots used by the Cross-Omics Dataset tab (`../../01_Data/`) and the Expression-Methylation Integration stage (`../../02_Expression_Methylation_Integration/`) — and, for `cx_empty_state` specifically, by Biomarker Convergence and Cross-Omics MR as well.

- **`crossomics_integration_helpers.R`** — the largest file in the vertical. Gene/CpG ID standardization (`cx_standardize_*`), evidence classification (`cx_classify*`, `cx_classify_evidence`), gene ID harmonization (`cx_harmonize_gene_ids`, also called from `R/transcriptomics/mod_biomarkercard.R` and `R/multiomics/mod_multi_biomarkercard.R`), default-data loaders (`cx_load_default_deg/dmr/methylation`).
- **`crossomics_integration_plots.R`** — `cx_empty_state` (shared empty-state UI across all 3 analytical stages) plus integration-specific plots including `cx_gene_cpg_network_plot()` (also reused by `R/multiomics/multiomics_concordance_plots.R`).

Not stage-specific by name (no `01_`/`02_` prefix) because these functions are genuinely called from more than one numbered stage folder — see `../../README.md`.
