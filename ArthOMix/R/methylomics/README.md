# Methylomics

Code backing the app's **Methylomics** tab. Every stage's config/UI/server trio is registered in `R/submodules_registry.R`'s `MX_MODULES` list (14 entries, in pipeline order), which `ui.R`'s `methylomicsUI()`/`build_submodule_grid()` and `server.R`'s generic `lapply(MX_MODULES, ...)` loop consume without hardcoding any one stage by name.

## Pipeline

```
01_Data (preloaded GSE42861 / GEO fetch by platform / upload of matrix+sample sheet+IDAT)
   ↓
02_Quality_Control → 03_Normalisation → 04_Cell_Type_Deconvolution
   ↓
05_Differential_Methylation_Position → 06_Differential_Methylation_Region → 07_Sex_Interaction_Analysis
   ↓
08_WGCNA_Co_Methylation_Network → 09_Candidate_CpGs
   ↓
10_ML_Feature_Selection → 11_Mendelian_Randomization → 12_Colocalization
   ↓
13_Diagnostic_Classifier → 14_Validation
   ↓
15_Biomarker_Analysis
```

## Stage table

| Folder | `MX_MODULES` id / title | Main file | Depends on (`functions/`) |
|---|---|---|---|
| `01_Data/` | Dataset tab (not an `MX_MODULES` entry) | `mod_methyl_dataset.R` | `idat_metrics.R`, `parse_upload.R` |
| `02_Quality_Control/` | `qc` | `mod_methyl_qc.R` | `qc.R`, `idat_metrics.R`, `annotation.R`, `parse_upload.R` |
| `03_Normalisation/` | `normalization` | `mod_methyl_normalization.R` | `normalization.R`, `qc.R`, `annotation.R`, `parse_upload.R` |
| `04_Cell_Type_Deconvolution/` | `celltype` | `mod_methyl_celltype.R` | `qc.R`, `annotation.R`, `parse_upload.R` |
| `05_Differential_Methylation_Position/` | `dmp` | `mod_methyl_dmp.R` | `qc.R`, `annotation.R`, `normalization.R` |
| `06_Differential_Methylation_Region/` | `dmr` | `mod_methyl_dmr.R` | `qc.R`, `annotation.R`, `normalization.R`; + `mod_methyl_dmp.R`'s shared helpers (see note below) |
| `07_Sex_Interaction_Analysis/` | `interaction` | `mod_methyl_interaction.R` | `qc.R`; + `mod_methyl_dmp.R`'s shared helpers |
| `08_WGCNA_Co_Methylation_Network/` | `wgcna` | `mod_methyl_wgcna.R` | `qc.R`, `normalization.R`; + `mod_methyl_dmp.R`'s shared helpers |
| `09_Candidate_CpGs/` | `candidates` | `mod_methyl_candidates.R` | self-contained (reads `results$dmr`/`results$wgcna`) |
| `10_ML_Feature_Selection/` | `featureselection` | `mod_methyl_featureselection.R` | `qc.R`, `annotation.R`, `parse_upload.R`; + `mod_methyl_dmp.R`'s shared helpers |
| `11_Mendelian_Randomization/` | `mr` | `mod_methyl_mr.R` | `annotation.R` only |
| `12_Colocalization/` | `coloc` | `mod_methyl_coloc.R` | self-contained |
| `13_Diagnostic_Classifier/` | `diagnostic` | `mod_methyl_diagnostic.R` | `qc.R`, `parse_upload.R` |
| `14_Validation/` | `validation` | `mod_methyl_validation.R` | `qc.R`, `parse_upload.R` |
| `15_Biomarker_Analysis/` | `biomarkercard` | `mod_methyl_biomarkercard.R` | `annotation.R`, `parse_upload.R`; + Cross-Omics/Multiomics (see below) |

## `functions/` (shared across multiple stages)

- **`qc.R`** — despite the name, the general-purpose statistics/plotting library for the whole vertical: filtering, batch correction (ComBat/RUVm), PCA/MDS, beta↔M conversion, sex-check, outlier detection, plots, QC report generation. Used by 9 of the 14 stages.
- **`annotation.R`** — array annotation lookup (`methyl_get_annotation`, `methyl_probe_is_cpg`). Used by 8 stages.
- **`parse_upload.R`** — matrix/sample-sheet/probe-list/IDAT upload parsing and validation. Used by almost every stage that accepts custom uploads.
- **`normalization.R`** — normalization method engine (noob/funnorm/SWAN/dasen/BMIQ/PBC/quantile), diagnostics, status detection. Used by Normalisation directly and by DMP/DMR/WGCNA for status checks.
- **`idat_metrics.R`** — IDAT-derived QC metrics. Used by the Dataset tab and QC.

## Cross-module helper coupling (not extracted, documented per `REFACTORING_NOTES.md`)

`mod_methyl_dmp.R` (`05_Differential_Methylation_Position/`) defines several helpers used by other stages directly rather than through `functions/`: `mod_methyl_dmp_sex_col`, `_sex_choices`, `_covariate_cols`, `_filter`, `_volcano`, `_betadist`, plus `mod_methyl_lambda_gc`, `mod_methyl_qq_plot`, `methyl_chunked_lmfit`, `mod_methyl_sva_fit`. `06_Differential_Methylation_Region/`, `07_Sex_Interaction_Analysis/`, `08_WGCNA_Co_Methylation_Network/`, and `10_ML_Feature_Selection/` all call into these. Moving `mod_methyl_dmp.R` doesn't break this — every file under `R/methylomics/` is still sourced into the same shared environment regardless of subfolder (see `R/0_load_omics_modules.R`) — but it means DMP isn't purely a leaf stage; keep this in mind before ever splitting `mod_methyl_dmp.R` itself.

## Cross-vertical coupling

`mod_methyl_biomarkercard.R` (`15_Biomarker_Analysis/`) calls `cx_harmonize_gene_ids()` (`R/crossomics/functions/integration/crossomics_integration_helpers.R`) and `mp_get_wikipathways_termgene()` (`R/multiomics/multiomics_pathway_helpers.R`) directly — Methylomics isn't fully self-contained.

See `../../CODE_MAP.md` for the app-wide code map and `../../PUBLICATION_PIPELINE.md` for how this stage sequence maps to the scientific workflow.
