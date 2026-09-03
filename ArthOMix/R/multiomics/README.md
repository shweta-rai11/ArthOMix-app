# Multiomics

Code backing the app's **Multiomics** tab. Distinct from **Cross-Omics** (`R/crossomics/`): Multiomics jointly models raw multi-layer data (DIABLO/SNF/MOFA2 on matched expression+methylation+other layers), while Cross-Omics correlates two already-computed single-omics result sets. Every stage's config/UI/server trio is registered in `R/modules_index.R`'s `MULTI_MODULES` list (8 entries), consumed generically by `ui.R`'s `multiomicsUI()`/`build_submodule_grid()` and `server.R`'s `lapply(MULTI_MODULES, ...)` loop.

## Pipeline

```
01_Data_Workspace (upload / GEO / preloaded fit; sample harmonization; optional live MOFA2)
   ↓
02_Cohort_Harmonization  (cohort/sample-level harmonization summary)
   ↓
03_DIABLO_SNF_Integration  (DIABLO discriminant analysis + SNF network fusion)
   ↓
04_SNF_Clustering  (consensus clustering on the fused network → patient subgroups)
   ↓
05_Biomarker_Discovery  (joint cross-layer biomarker discovery, sex-stratified)
   ↓
06_Gene_CpG_Mapping  (per gene-CpG pair expression/methylation correlation)
   ↓
07_Pathways  (pathway-level integration over discovered biomarkers)
   ↓
08_Biomarker_Card  (read-only integrated per-biomarker interpretation)
   ↓
09_Results_Summary  (reproducibility: package versions, run summary)
```

## Stage table

| Folder | `MULTI_MODULES` id / title | Main file(s) | Depends on (`functions/`) |
|---|---|---|---|
| `01_Data_Workspace/` | Dataset tab (not a `MULTI_MODULES` entry) | `mod_multi_dataset.R` (largest file in the vertical), `multiomics_dataset_helpers.R`, `multiomics_dataset_plots.R`, `mod_multi_mofa.R` + `mod_multi_mofa_engine.R` (live MOFA2, mounted *inside* this tab, not a separate `MULTI_MODULES` stage) | — |
| `02_Cohort_Harmonization/` | `overview` / "Cohort Harmonization" | `mod_multi_overview.R`, `cohort_harmonization_helpers.R`, `cohort_harmonization_plots.R` | — |
| `03_DIABLO_SNF_Integration/` | `integration` / "Multi-omics Integration (DIABLO & SNF)" | `mod_multi_integration.R` | `multiomics_integration_helpers.R`, `multiomics_integration_plots.R`, `multiomics_plots.R`, `multiomics_sexstratified_engine.R` |
| `04_SNF_Clustering/` | `stratification` / "SNF Clustering" | `mod_multi_stratification.R`, `snf_clustering_helpers.R`, `snf_clustering_plots.R` | `multiomics_integration_helpers.R`/`_plots.R` (see duplicate-SNF-logic note below) |
| `05_Biomarker_Discovery/` | `biomarker` / "Biomarker Discovery" | `mod_multi_biomarker.R`, `multiomics_biomarker_helpers.R`, `multiomics_biomarker_plots.R` | `multiomics_sexstratified_engine.R` |
| `06_Gene_CpG_Mapping/` | `mapping` / "Gene–CpG Mapping" | `mod_multi_mapping.R`, `multiomics_mapping_helpers.R`, `multiomics_mapping_plots.R` | — (reuses a plotting technique from `R/crossomics/functions/integration/crossomics_integration_plots.R::cx_gene_cpg_network_plot()`) |
| `07_Pathways/` | `pathway` / "Pathways" | `mod_multi_pathway.R`, `multiomics_pathway_helpers.R`, `multiomics_pathway_plots.R` | — |
| `08_Biomarker_Card/` | `biomarkercard` / "Biomarker Card" | `mod_multi_biomarkercard.R` | reads `multi_results$mapping$df`; calls `cx_classify_evidence()` (Cross-Omics) |
| `09_Results_Summary/` | `summary` / "Results Summary & Reproducibility" | `mod_multi_summary.R` | `multiomics_helpers.R` (`multi_package_versions`, `multi_analysis_summary_table`) |

## `functions/` (shared across multiple stages)

- **`multiomics_helpers.R`** — cross-cutting app-wide multiomics utilities (`multi_read_table`, `multi_active_dataset_banner`, `multi_qc_scorecard`, `multi_package_versions`, `multi_build_report`). Used by nearly every `mod_multi_*.R`.
- **`multiomics_plots.R`** — shared plotting utilities (`multi_empty_state`, `multi_plot_or_empty`, `multi_diablo_*_plot`). Used across nearly every stage.
- **`multiomics_sexstratified_engine.R`** — shared sex-stratified nested-CV engine (`mss_*`). Used by `03_DIABLO_SNF_Integration/` and `05_Biomarker_Discovery/`.
- **`multiomics_integration_helpers.R`** / **`multiomics_integration_plots.R`** — DIABLO/SNF implementation (`mi_diablo_*`, `mi_snf_*`). Primarily `03_DIABLO_SNF_Integration/`'s, but also called from `04_SNF_Clustering/` (see below) — placed in `functions/` rather than folder-scoped to Integration alone.

## SNF Clustering vs. DIABLO & SNF Integration — not duplicate logic

`04_SNF_Clustering/mod_multi_stratification.R` calls both its own `sfc_snf_*` helpers (`snf_clustering_helpers.R`, same folder) and `mi_snf_*`/`mi_ari` from `functions/multiomics_integration_helpers.R`/`_plots.R`. This looks like two parallel SNF code paths but isn't one: `sfc_snf_run()` delegates directly to `mi_snf_run()` for every ≥2-block case; the only logic unique to `snf_clustering_helpers.R` is a single-omics (1-block) fallback that `mi_snf_run()` structurally cannot serve. See `REFACTORING_NOTES.md`.

## Naming note

Folder names/order follow the **actual** `MULTI_MODULES` registry titles rather than the task template's guessed names — e.g. the `stratification` stage's real title is "SNF Clustering" (not "Patient Stratification"), and there is no separate "Data Quality Control"/"Data Preprocessing" stage as distinct code in this vertical (those concerns live inside the Dataset Workspace tab and Cohort Harmonization).

See `../../CODE_MAP.md` and `../../PUBLICATION_PIPELINE.md`.
