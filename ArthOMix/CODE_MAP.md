# ArthOMix Code Map

Maps every analysis in the app to its code location, UI entry point, and dependencies. This is a living document — being filled in vertical-by-vertical as the publication-ready reorganization proceeds (see `REFACTORING_NOTES.md` for what's not yet done).

```
ArthOMix
├── Transcriptomics   (R/transcriptomics/)      — not yet reorganized, see below
├── Methylomics       (R/methylomics/)          — reorganized
├── Cross-Omics       (R/crossomics/)           — reorganized
└── Multiomics        (R/multiomics/)           — reorganized
```

## Shared app infrastructure (not part of any one vertical)

| Location | Purpose |
|---|---|
| `ui.R`, `server.R`, `global.R`, `data_paths.R` | App root, required by Shiny's own directory convention (`shiny::runApp()` auto-sources `ui.R`/`server.R`/`global.R`; `data_paths.R` is explicitly `source()`d once from `global.R` before `R/` loads, to avoid double-sourcing). |
| `R/0_load_omics_modules.R` | Recursively sources every `.R` file under each of the 4 vertical folders (works around `shiny:::loadSupport()` only scanning `R/*.R` non-recursively). |
| `R/0a_load_auth_modules.R` | Sources `R/auth/*.R`. |
| `R/0b_load_shared_modules.R` | Sources `R/shared/*.R` (cross-vertical modules — currently just ArthOChat). |
| `R/submodules_registry.R` | Assembles `TX_MODULES`/`MX_MODULES`/`CX_MODULES`/`MULTI_MODULES` (and `*_MODULES_BY_ID`) from every `mod_*_config/_ui/_server` trio; defines the `build_tx_context()`/`build_mx_context()`/`build_cx_context()`/`build_mo_context()` functions ArthOChat uses to summarize each vertical's live results. |
| `R/ui_shell.R` | Shared SaaS-dashboard shell: `app_header()`, `omics_sidebar()`, `pipeline_summary_ui()` — presentational only, reused by all 4 verticals. |
| `R/auth/` | Supabase-backed sign-up/login/password-reset (`mod_auth_ui.R`, `mod_auth_server.R`, `auth_api.R`). |
| `R/shared/mod_arthochat.R` | The AI assistant ("ArthOChat"), app-wide (spans all 4 verticals' context builders) — relocated here from `R/transcriptomics/` since it isn't transcriptomics-specific. |
| `R/provenance.R` | Shared provenance-manifest helpers (`arthomix_provenance_record()`, `arthomix_provenance_download_handler()`) used across modules. |

## Cross-Omics

See `R/crossomics/README.md` for the full narrative. Summary:

| Folder | Stage (`CX_MODULES` id / title) | Main file(s) | UI |
|---|---|---|---|
| `R/crossomics/01_Data/` | Dataset tab (not a `CX_MODULES` entry) | `mod_cross_dataset.R`, `crossomics_integration_upload.R` | Cross-Omics → Dataset |
| `R/crossomics/02_Expression_Methylation_Integration/` | `integration` / "Expression and Methylation" | `mod_cross_integration.R` | Cross-Omics → Sub-modules → Expression and Methylation |
| `R/crossomics/03_Biomarker_Convergence/` | `biomarkerconv` / "Biomarker Convergence" | `mod_cross_biomarker_conv.R` | Cross-Omics → Sub-modules → Biomarker Convergence |
| `R/crossomics/04_Cross_Omics_MR/` | `mrstage` / "Cross-Omics MR" | `mod_cross_mr_stage.R`, `crossomics_mrstage_helpers.R` | Cross-Omics → Sub-modules → Cross-Omics MR |
| `R/crossomics/functions/integration/` | shared (Dataset + Integration; `cx_empty_state` used everywhere) | `crossomics_integration_helpers.R`, `crossomics_integration_plots.R` | — |
| `R/crossomics/functions/biomarker_convergence/` | shared (Biomarker Convergence + Cross-Omics MR) | `crossomics_biomarkerconv_helpers.R` | — |

**Publication relevance**: this pipeline is what supports any figure/table showing gene-methylation concordance classification, cross-omics biomarker ranking, or MR estimates on cross-omics-convergent candidates.

## Methylomics

See `R/methylomics/README.md` for the full narrative and dependency table. Summary:

| Folder | Stage (`MX_MODULES` id) | Main file |
|---|---|---|
| `R/methylomics/01_Data/` | Dataset tab (not an `MX_MODULES` entry) | `mod_methyl_dataset.R` |
| `R/methylomics/02_Quality_Control/` | `qc` | `mod_methyl_qc.R` |
| `R/methylomics/03_Normalisation/` | `normalization` | `mod_methyl_normalization.R` |
| `R/methylomics/04_Cell_Type_Deconvolution/` | `celltype` | `mod_methyl_celltype.R` |
| `R/methylomics/05_Differential_Methylation_Position/` | `dmp` | `mod_methyl_dmp.R` |
| `R/methylomics/06_Differential_Methylation_Region/` | `dmr` | `mod_methyl_dmr.R` |
| `R/methylomics/07_Sex_Interaction_Analysis/` | `interaction` | `mod_methyl_interaction.R` |
| `R/methylomics/08_WGCNA_Co_Methylation_Network/` | `wgcna` | `mod_methyl_wgcna.R` |
| `R/methylomics/09_Candidate_CpGs/` | `candidates` | `mod_methyl_candidates.R` |
| `R/methylomics/10_ML_Feature_Selection/` | `featureselection` | `mod_methyl_featureselection.R` |
| `R/methylomics/11_Mendelian_Randomization/` | `mr` | `mod_methyl_mr.R` |
| `R/methylomics/12_Colocalization/` | `coloc` | `mod_methyl_coloc.R` |
| `R/methylomics/13_Diagnostic_Classifier/` | `diagnostic` | `mod_methyl_diagnostic.R` |
| `R/methylomics/14_Validation/` | `validation` | `mod_methyl_validation.R` |
| `R/methylomics/15_Biomarker_Analysis/` | `biomarkercard` | `mod_methyl_biomarkercard.R` |
| `R/methylomics/functions/` | shared (`qc.R`, `annotation.R`, `parse_upload.R`, `normalization.R`, `idat_metrics.R`) | — |

## Multiomics

See `R/multiomics/README.md` for the full narrative and dependency table. Summary:

| Folder | Stage (`MULTI_MODULES` id) | Main file(s) |
|---|---|---|
| `R/multiomics/01_Data_Workspace/` | Dataset Workspace tab (not a `MULTI_MODULES` entry; hosts live MOFA2) | `mod_multi_dataset.R`, `mod_multi_mofa.R`, `mod_multi_mofa_engine.R` |
| `R/multiomics/02_Cohort_Harmonization/` | `overview` | `mod_multi_overview.R` |
| `R/multiomics/03_DIABLO_SNF_Integration/` | `integration` | `mod_multi_integration.R` |
| `R/multiomics/04_SNF_Clustering/` | `stratification` | `mod_multi_stratification.R` |
| `R/multiomics/05_Biomarker_Discovery/` | `biomarker` | `mod_multi_biomarker.R` |
| `R/multiomics/06_Gene_CpG_Concordance/` | `concordance` | `mod_multi_concordance.R` |
| `R/multiomics/07_Pathways/` | `pathway` | `mod_multi_pathway.R` |
| `R/multiomics/08_Biomarker_Card/` | `biomarkercard` | `mod_multi_biomarkercard.R` |
| `R/multiomics/09_Results_Summary/` | `summary` | `mod_multi_summary.R` |
| `R/multiomics/functions/` | shared (`multiomics_helpers.R`, `multiomics_plots.R`, `multiomics_sexstratified_engine.R`, `multiomics_integration_helpers.R`, `multiomics_integration_plots.R`) | — |

## Transcriptomics

Not yet reorganized into numbered stage folders — still a flat directory (`R/transcriptomics/`), with a `TX_MODULES` registry in `R/submodules_registry.R` giving the authoritative stage list, titles, and order. Each stage's file can be found by grepping `R/submodules_registry.R` for the stage's `config` object name (e.g. `mod_dge_config`) and locating the matching `mod_*.R` file in `R/transcriptomics/`.
