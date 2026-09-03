# Transcriptomics

Code backing the app's **Transcriptomics** tab. Every stage's config/UI/server trio is registered in `R/modules_index.R`'s `TX_MODULES` list (16 entries, in pipeline order), consumed generically by `ui.R`'s `build_submodule_grid()` and `server.R`'s `lapply(TX_MODULES, ...)` loop.

Unlike Cross-Omics/Methylomics/Multiomics, this vertical is the oldest in the app and was never split into per-stage `<stage>_helpers.R`/`<stage>_plots.R` files — each `mod_<id>.R` bundles config+UI+server+statistics+plots in one file (several 1500-2500+ lines). Per this reorganization's own rule ("organization, not rewriting"), **files were moved as-is, not split** — splitting a working file into multiple scripts is a functional-risk change with no organizational necessity once each stage already has its own numbered folder.

## Pipeline

```
01_Data (upload / GEO / preloaded)
   ↓
02_Overview (QC dashboard + pipeline status roll-up)
   ↓
03_Preprocessing_Batch_Correction
   ↓
04_Differential_Expression → 05_WGCNA → 06_Candidate_Gene_Identification
   ↓
07_Mendelian_Randomization → 08_Colocalization
   ↓
09_Feature_Selection → 10_Diagnostic_Model → 11_Sex_Interaction_Analysis
   ↓
12_Cross_Tissue_Validation → 13_Cross_Ancestral_Validation
   ↓
14_Functional_Enrichment → 15_Immune_Deconvolution
   ↓
16_Nomogram
   ↓
17_Biomarker_Card
```

## Stage table

| Folder | `TX_MODULES` id / title | Main file |
|---|---|---|
| `01_Data/` | Dataset tab (not a `TX_MODULES` entry) | `mod_dataset.R` |
| `02_Overview/` | `overview` | `mod_overview.R` |
| `03_Preprocessing_Batch_Correction/` | `preprocessing` | `mod_preprocessing.R` (+ `mod_preprocessing_explore.R`, its EDA sub-tab) |
| `04_Differential_Expression/` | `dge` | `mod_dge.R` |
| `05_WGCNA/` | `wgcna` | `mod_wgcna.R` |
| `06_Candidate_Gene_Identification/` | `candidates` | `mod_candidates.R` |
| `07_Mendelian_Randomization/` | `mr` | `mod_mr.R` |
| `08_Colocalization/` | `coloc` | `mod_coloc.R` |
| `09_Feature_Selection/` | `featureselection` | `mod_featureselection.R` |
| `10_Diagnostic_Model/` | `diagnostic` | `mod_diagnostic.R` |
| `11_Sex_Interaction_Analysis/` | `interaction` | `mod_interaction.R` |
| `12_Cross_Tissue_Validation/` | `crosstissue` | `mod_crosstissue.R` |
| `13_Cross_Ancestral_Validation/` | `crossancestry` | `mod_crossancestry.R` |
| `14_Functional_Enrichment/` | `enrichment` | `mod_enrichment.R` |
| `15_Immune_Deconvolution/` | `deconvolution` | `mod_deconvolution.R` |
| `16_Nomogram/` | `nomogram` | `mod_nomogram.R` |
| `17_Biomarker_Card/` | `biomarkercard` | `mod_biomarkercard.R` |
| `functions/` | shared (`expression_type.R`) | — |

## `functions/`

- **`expression_type.R`** — expression-matrix scale/type heuristics (raw-counts vs. normalized detection) and the declare-then-verify upload validator. Genuinely shared across 3 stages: `01_Data/mod_dataset.R`'s upload path, `04_Differential_Expression/mod_dge.R`'s decoupled upload path, and `15_Immune_Deconvolution/mod_deconvolution.R`'s run gate. Mirrors `R/methylomics/functions/parse_upload.R`'s validation pattern.

## Note

`R/shared/mod_arthochat.R` (the AI assistant) used to live in this folder but was relocated during the Cross-Omics reorganization pass — it's app-wide (spans all 4 verticals), not transcriptomics-specific. See `../shared/mod_arthochat.R` and `../../CODE_MAP.md`.

See `../../CODE_MAP.md` and `../../PUBLICATION_PIPELINE.md`.
