# Cross-Omics

Code backing the app's **Cross-Omics** tab (`ui.R`'s `crossomicsUI()`, mounted at `ui.R:1620`), which integrates transcriptomics (DEG) and methylomics (DMP/DMR) results computed elsewhere in the app. Cross-Omics is distinct from **Multiomics** (`R/multiomics/`): it correlates two *already-computed, single-omics* result sets (gene-level expression + CpG-level methylation), rather than jointly modeling raw multi-layer data (DIABLO/SNF/MOFA2).

## Pipeline

```
01_Data  (upload / preloaded DEG+DMR / sample & feature matching)
   ↓
02_Expression_Methylation_Integration  (gene↔CpG classification: Hyper+Down, Hypo+Up, ...)
   ↓
03_Biomarker_Convergence  (evidence-tier ranking across expression+methylation)
   ↓
04_Cross_Omics_MR  (Mendelian Randomization on convergent hits)
```

This mirrors the 3 real analytical stages registered in `R/modules_index.R`'s `CX_MODULES` (`integration`, `biomarkerconv`, `mrstage`), plus the Dataset tab (`01_Data/`, not itself a `CX_MODULES` entry — see `R/multiomics/README.md` for why dataset tabs are structured this way in every vertical).

## Folders

- **`01_Data/`** — the "Dataset" tab: upload/preloaded DEG+DMR/sample-methylation loading, sample/feature matching. Publishes the shared `cross_dataset` reactiveValues every stage below reads from.
- **`02_Expression_Methylation_Integration/`** — "Expression and Methylation" stage: classifies each gene/CpG pair by concordance direction.
- **`03_Biomarker_Convergence/`** — "Biomarker Convergence" stage: ranks candidates by cross-omics evidence strength.
- **`04_Cross_Omics_MR/`** — "Cross-Omics MR" stage: two-sample Mendelian Randomization on convergent candidates.
- **`functions/integration/`** — helpers/plots shared by the Dataset tab and the Integration stage (and `cx_empty_state`, reused by Biomarker Convergence and Cross-Omics MR too).
- **`functions/biomarker_convergence/`** — helpers shared by the Biomarker Convergence and Cross-Omics MR stages.

## Dependencies outside this folder

- `R/modules_index.R` assembles `CX_MODULES`/`CX_MODULES_BY_ID` from each stage's `mod_cross_*_config/_ui/_server` trio, and `build_cx_context()` for ArthOChat.
- `server.R` mounts `mod_cross_dataset_server` directly, then loops `CX_MODULES` (the `integration`/`mrstage` entries take extra arguments — see `server.R`'s special-cased calls).
- Reads Transcriptomics (`results$dge`) and Methylomics (`methyl_results$dmp`/`dmr`) results directly, passed in from `server.R`.

## Known cross-vertical reuse (not moved here, referenced from elsewhere)

- `R/transcriptomics/mod_biomarkercard.R` and `R/multiomics/mod_multi_biomarkercard.R` call `cx_classify_evidence()`/`cx_harmonize_gene_ids()` from `functions/integration/crossomics_integration_helpers.R`.
- `R/multiomics/multiomics_mapping_plots.R` reuses `cx_gene_cpg_network_plot()` from `functions/integration/crossomics_integration_plots.R`.

See `../../CODE_MAP.md` for the full app-wide code map and `../../PUBLICATION_PIPELINE.md` for how this stage sequence maps to the broader scientific workflow.
