# 01_Data_Workspace

The Multiomics "Dataset Workspace" tab — data-loading entry point, not a `MULTI_MODULES` stage.

- **Files**: `mod_multi_dataset.R` (largest file in the vertical — source-selection UI, upload parsing, GEO fetch, sample harmonization/matching, provenance), `multiomics_dataset_helpers.R`, `multiomics_dataset_plots.R`, plus `mod_multi_mofa.R`/`mod_multi_mofa_engine.R` (a real, live MOFA2 factor analysis run on the Active Multi-Omics Dataset, mounted directly inside this tab as "Integrated Analysis (MOFA2)" rather than as its own `MULTI_MODULES` entry).
- **Input**: upload, GEO, or preloaded multi-omics fit.
- **Main operation**: builds the "Active Multi-Omics Dataset" — matched samples across layers, validated and harmonized.
- **Output**: the shared `multi_dataset`/`multi_results` reactiveValues every stage below reads from.
- **UI**: Multiomics → "Dataset Workspace" tab. Mounted directly in `server.R` (`mod_multi_dataset_server(...)`), outside the `MULTI_MODULES` loop.
