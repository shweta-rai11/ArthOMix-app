# The Multi-Omics Module — Thesis Documentation

**Scope of this document.** Everything below describes the Multi-Omics vertical of the ArthOMix Shiny application exactly as implemented in `R/multiomics/` (Dataset Workspace, seven registered analytical sub-modules, and the MOFA2 "Integrated Analysis" feature embedded inside the Dataset Workspace), plus the small amount of shared infrastructure in `R/submodules_registry.R`, `R/0_load_omics_modules.R`, and `server.R` that wires these pieces together. No other application module (Transcriptomics, Methylomics, Cross-Omics) is documented here except where it is the literal source of a function the Multi-Omics module reuses (this happens twice: a methylation beta/M-value conversion helper and the Cross-Omics CpG-annotation lookup, both noted explicitly where used).

**Method.** Every claim in this document was verified by reading the actual R source files line by line — UI definitions, server logic, reactive expressions, and the pure-logic helper functions they call — not inferred from UI labels, button text, or variable names. Where a UI control exists but does not drive any computation, or where a feature is incomplete or unreachable, this is stated explicitly in [Section 14](#14-known-gaps-decorative-controls-and-unimplemented-functionality) rather than silently omitted. File:line citations are given for load-bearing claims so they can be checked directly against the codebase.

---

## Table of Contents

1. [Multi-Omics Overview](#1-multi-omics-overview)
2. [Overall Architecture](#2-overall-architecture)
3. [Dataset Inputs — the Dataset Workspace and MOFA2](#3-dataset-inputs--the-dataset-workspace-and-mofa2)
4. [Data Flow](#4-data-flow)
5. [Every Implemented Sub-module](#5-every-implemented-sub-module)
6. [Scientific Importance](#6-scientific-importance)
7. [End-to-End Code Walkthrough](#7-end-to-end-code-walkthrough)
8. [Important R/Shiny Concepts Used in This Module](#8-important-rshiny-concepts-used-in-this-module)
9. [Integration Methods in Detail](#9-integration-methods-in-detail)
10. [Visualization and Output Logic](#10-visualization-and-output-logic)
11. [Input → Processing → Output Summary Table](#11-input--processing--output-summary-table)
12. [End-to-End Workflow Summary](#12-end-to-end-workflow-summary)
13. [Known Gaps, Decorative Controls, and Unimplemented Functionality](#13-known-gaps-decorative-controls-and-unimplemented-functionality)
14. [Thesis-Level Main Paragraph](#14-thesis-level-main-paragraph)
15. [Short Thesis Paragraph](#15-short-thesis-paragraph)

---

## 1. Multi-Omics Overview

The Multi-Omics module is one of four analytical verticals in ArthOMix (alongside Transcriptomics, Methylomics, and Cross-Omics). Its purpose is to combine gene-expression (RNA-seq) and DNA-methylation data from the same patients — the rheumatoid-arthritis (RA) anti-TNF cohort of Tao et al. (2021) that ships with the app, or a user's own uploaded/GEO-retrieved data — into a single analytical workflow that a single-omics view cannot provide: harmonizing which samples and modalities are actually usable together, fusing the two layers into one supervised or unsupervised model, grouping patients by their combined molecular profile, and connecting a specific gene's expression change to the methylation change at its own regulatory CpGs, before finally asking what biological pathways the resulting candidate list implicates.

Structurally, the module consists of:
- **Dataset Workspace** (`mod_multi_dataset.R`) — not a registered analytical sub-module but the entry point every other piece depends on; it is where a dataset is selected, validated, sample-matched, preprocessed, batch-corrected, and finally "activated" into a shared in-memory object.
- **Seven registered sub-modules** (`MULTI_MODULES` in `R/submodules_registry.R:84-91`), each its own tab in the app's left navigation once added from the "Sub-modules" grid: Cohort Harmonization, Multi-omics Integration (DIABLO & SNF), SNF Clustering (Patient Stratification), Biomarker Discovery, Gene–CpG Concordance, Pathways, and Results Summary & Reproducibility.
- **Integrated Analysis (MOFA2)**, a live factor-analysis feature mounted directly inside the Dataset Workspace tab rather than as its own registry entry.

Every analytical sub-module can operate in one of two modes, selected by a `data_source` radio button: **"Active Multi-Omics Dataset"**, which reads whatever the Dataset Workspace has activated (an upload, a GEO series, or the bundled reference cohort run through the same pipeline), or **"Reference/Example Dataset"**, a shortcut that recomputes live from one of six precomputed, already-fitted analysis cells of the bundled RA anti-TNF cohort (see [§2.3](#23-the-reference-example-dataset-mechanism)). Nothing in any sub-module computes automatically — every expensive analysis is behind an explicit "Run" button, and every result panel has a defined "not run yet" state rather than a blank or fabricated one.

---

## 2. Overall Architecture

### 2.1 The shared state objects

Two `reactiveValues` objects are created once, at the application's `server.R`, and passed by reference into every Multi-Omics sub-module's server function:

```r
# server.R:140-147
multi_dataset <- reactiveValues(
  table_label = NULL, df = NULL, source = NULL,
  layers = list(), layer_meta = list(), sample_meta = NULL,
  overlap = NULL, active = FALSE, loaded_at = NULL
)
multi_results <- reactiveValues()
mod_multi_dataset_server("mo_dataset", multi_dataset, multi_results)
lapply(MULTI_MODULES, function(m) m$server(paste0("mo_", m$config$id), multi_dataset, multi_results))
```

`multi_dataset` is the module's single "what data are we working with" object. `mod_multi_dataset_server()` (Dataset Workspace) is the **only** code in the whole module that writes to `multi_dataset$layers`, `$layer_meta`, `$sample_meta`, `$active`, and `$source` — every other sub-module reads it, never writes it (the one narrow exception is that Cohort Harmonization, Integration, Biomarker Discovery, and SNF Clustering additionally know how to bypass `multi_dataset` entirely and pull a "Reference/Example" analysis cell straight from disk when the user picks that data source instead).

`multi_results` starts as a bare, empty `reactiveValues()`. Each sub-module claims one named slot on it once its own analysis finishes (`multi_results$overview`, `$integration`, `$stratification`, `$biomarker`, `$concordance`, `$pathway`, plus `$live_qc`/`$live_mofa` from the embedded MOFA2 feature). This is how, for example, Pathways can read Biomarker Discovery's selected features, or Results Summary can build its session dashboard, without any sub-module directly calling into another's code — they only ever read a plain list another sub-module already published.

Every `MULTI_MODULES` entry is instantiated with the identical call signature `(namespace_id, multi_dataset, multi_results)` — there is no per-module argument variation.

### 2.2 The sub-module registry and auto-sourcing

`R/submodules_registry.R:84-91` lists the seven registered sub-modules in this fixed order:

```r
MULTI_MODULES <- list(
  list(config = mod_multi_overview_config,       ui = ..., server = ...),  # Cohort Harmonization
  list(config = mod_multi_integration_config,    ui = ..., server = ...),  # Multi-omics Integration
  list(config = mod_multi_stratification_config, ui = ..., server = ...),  # SNF Clustering
  list(config = mod_multi_biomarker_config,      ui = ..., server = ...),  # Biomarker Discovery
  list(config = mod_multi_concordance_config,    ui = ..., server = ...),  # Gene–CpG Concordance
  list(config = mod_multi_pathway_config,        ui = ..., server = ...)   # Pathways
)
```

`mod_multi_summary` (Results Summary & Reproducibility) is likewise mounted by the same server-side loop pattern, grouped alongside Pathways under `group = "Interpretation"` in its own config. Each `_config` list carries an `id`, a display `title`, an icon, and a `group` string — `"Data"` (Cohort Harmonization, Integration, SNF Clustering), `"Biomarker modeling"` (Biomarker Discovery, Concordance), or `"Interpretation"` (Pathways, Results Summary) — which is purely what determines which card-grid section a sub-module appears under on the "Sub-modules" tab; it has no effect on computation.

New R files dropped into `R/multiomics/` never require a registry change: `R/0_load_omics_modules.R:21-25` sources every `.R` file under `R/multiomics/` (and the other three vertical folders) in plain alphabetical order via `list.files(..., pattern="\\.[rR]$") |> sort() |> source()`. A registry edit in `submodules_registry.R` is only needed when a *new sub-module* — a new `mod_multi_<id>_config`/`_ui`/`_server` triple that should get its own navigation tab — is added.

### 2.3 The "Reference/Example Dataset" mechanism

The bundled RA anti-TNF cohort is not stored in the app as one big raw expression/methylation matrix. Instead, `MULTI_CELLS` (`multiomics_helpers.R:52-63`) lists six **analysis cells** the offline research pipeline already ran — `female_Adalimumab`, `male_Adalimumab`, `female_Etanercept`, `male_Etanercept` (the only two with a precomputed SNF fit), `female_response`, `male_response` — each with its own saved `mixOmics::block.splsda` fit object (`MULTI_DIABLO_FIT_REGISTRY`, an `.rds` file per cell). `mi_preloaded_cell_dataset(cell_key)` (`multiomics_integration_live_helpers.R:35-56`) turns one such saved fit into a live-computation-ready dataset by extracting that fit's own `$X` matrices (one per omics block) and `$Y` outcome vector, relabeling the blocks "Transcriptomics"/"Methylomics" via `MULTI_BLOCK_LABELS`. Because a DIABLO fit's `$X` is already that cell's matched-sample, already-feature-selected subset (tens to a few hundred features, not the full transcriptome/methylome), every sub-module that offers "Reference/Example Dataset" as a data source is explicit in its provenance text that this re-derives a live-computation input from a saved fit, not the pipeline's original raw data. `mi_preloaded_cell_dataset()` is used by Multi-omics Integration, Biomarker Discovery, Cohort Harmonization, the Dataset Workspace itself, and SNF Clustering.

Separately, `MULTI_TABLE_REGISTRY` (`data_paths.R`) is a label → CSV-path lookup for the offline pipeline's own finished result tables (DIABLO performance tables, the Gene↔CpG concordance tables, DIABLO candidate-biomarker panels, genome-wide DEG/DMP candidate lookups). `multi_read_registry_table(label)` (`multiomics_helpers.R:39-42`) is a fail-soft CSV loader used by the "preloaded" data path of Gene–CpG Concordance and Pathways to *browse* these already-computed tables rather than recompute them.

### 2.4 The ArthOChat context bridge

`build_mo_context()` (`submodules_registry.R:219-236`) flattens the current `multi_dataset`/`multi_results` state into a short text block fed to the in-app AI assistant (ArthOChat), scoped to whichever sub-module tab is currently open. For a sub-module with no result yet, it explicitly returns "NOT YET RUN IN THIS SESSION... do not answer with a number from methodology/literature instead" — a deliberate guard against the assistant fabricating a number that was never actually computed in the running session.

---

## 3. Dataset Inputs — the Dataset Workspace and MOFA2

`mod_multi_dataset.R` (~1,250 lines) is the single entry point for getting data into the Multi-Omics module. It offers three source modes and a fixed five-step pipeline that runs identically regardless of source.

### 3.1 Three data sources

| Source | Mechanism |
|---|---|
| **Reference / Example Dataset** | Pick a `MULTI_CELL_CHOICES` analysis cell, click "Load Reference Dataset" → `mi_preloaded_cell_dataset()` (§2.3) populates the same `raw$mats`/`raw$validations`/`raw$meta` reactiveValues an upload would, so it runs through the identical steps 1–5 below. |
| **Upload Dataset** | Up to 8 dynamic blocks (`MO_MAX_BLOCKS`), each with an omics-type selector, a display label, a file input (CSV/TSV/TXT/XLSX/RDS), an optional per-dataset metadata file, and wide/long-format + orientation controls that are auto-detected but always user-confirmable, never silently applied. |
| **Retrieve from GEO** | A GEO-accession text field, a "Fetch from GEO" button, and a metadata upload. Fetching calls `multi_geo_layer_fetch()` (`GEOquery::getGEO()`), and multi-platform series prompt the user to pick one before `multi_geo_platform_matrix()` extracts the samples×features matrix (with optional gene-symbol probe collapsing for RNA-seq). |

Every upload/GEO block's declared omics type is structurally corroborated before it is accepted: `multi_live_detect_omics_type()` (`multiomics_live_helpers.R:260-277`) inspects the feature-ID pattern (Illumina CpG probe ID vs. gene symbol/Ensembl/Entrez ID) and the value range (beta-scale vs. other) and, for RNA-seq or DNA-methylation specifically, **rejects the dataset outright** if the structure contradicts the declared type — it never silently accepts a mislabeled file. Other omics types (proteomics, metabolomics, etc.) are accepted without this structural check, since the detector only distinguishes RNA-seq from methylation.

### 3.2 The five-step pipeline, in code order

1. **Preview and Validate** — for each surviving block: `multi_live_read_matrix()` reads the file (transposing if the orientation is `features_rows`), long-format tables are pivoted wide via `multi_live_pivot_long()` (aggregating duplicate sample×feature pairs by their mean), and the resulting matrix is scored by `multi_live_validate_matrix()` (sample/feature counts, % missing, duplicate IDs, zero-variance features, non-finite values — a pure description, never an automatic fix). Per-dataset and shared sample metadata are merged via `mo_merge_sample_meta()`.
2. **Sample Matching** — `overlap()` computes the intersection of sample IDs across every selected layer via `multi_live_sample_overlap()` (a plain `Reduce(intersect, ...)` on row names, never a row-position assumption), optionally after a chosen matching method (exact ID / patient-ID-column remap / an uploaded ID-mapping file) has translated each layer's own IDs onto a common patient identifier.
3. **Preprocessing** — hard-gated on ≥3 matched samples; subsets every layer to the matched IDs, then applies explicit missing-value handling (mean/median imputation, or row/column removal by a missingness threshold — never automatic), omics-appropriate normalization (`multi_live_normalize()`: log2, M-value logit, median, quantile via `limma::normalizeQuantiles`, Pareto, autoscale, or none), and an optional top-N variance feature filter.
4. **Batch Diagnostics** — PCA and sample-correlation heatmaps before and after an optional batch correction (`sva::ComBat` or `limma::removeBatchEffect`); a batch/phenotype confounding check (`multi_live_confounding_check()`, flags a batch level that maps to exactly one phenotype level) that **blocks** the "Apply batch correction" button unless the user explicitly checks an override box; a quantitative R²-of-PCs-vs-group diagnostic before vs. after correction; and a CSV download of the corrected matrix.
5. **Compatibility and Activate** — `multi_dataset_compatibility()` rolls every layer's individual Ready/Review-Required/Not-Compatible status plus the sample-matching outcome and metadata availability into one overall verdict (READY / READY WITH REVIEW / REVIEW REQUIRED / NOT READY). Clicking "Activate" is the sole moment `multi_dataset$layers`, `$layer_meta`, `$sample_meta`, `$source`, and `$active <- TRUE` are written — this is the hand-off point every other sub-module's "Active Multi-Omics Dataset" option depends on.

### 3.3 Integrated Analysis (MOFA2)

Mounted directly inside the Dataset Workspace tab (`mod_multi_live.R`/`mod_multi_live_mofa.R`), MOFA2 is the one place in the Dataset Workspace where a real multi-omics model is fit on whatever the user has just activated (requires `multi_dataset$active == TRUE` and ≥2 layers). It is a genuinely **unsupervised** factor-analysis method — it needs no outcome column at all.

- **User-controllable parameters**: number of factors (2–25, default 10), a random seed, and a convergence-speed preset (fast/medium/slow).
- **Fixed/hardcoded**: `save_data = FALSE`, `use_basilisk = FALSE` (deliberately disabled — the code comments record that the `basilisk`-managed Python environment segfaulted in this deployment; a plain `reticulate::use_python()` binding is used instead, re-established on every call since it cannot cross a background-worker process boundary), and MOFA2's own default model/training options otherwise.
- **Execution**: whichever number of factors was requested is silently capped at `min(sample count per view) − 1`. Missing values already present in the input matrices (whatever the Dataset Workspace's own imputation step left) are passed straight to `MOFA2::create_mofa()` and handled by MOFA2's own native likelihood-based missing-data support — this app's code does not impute anything further for MOFA2 specifically.
- **Functions actually called**: `MOFA2::create_mofa()` → `get_default_model_options()`/`get_default_training_options()` → `prepare_mofa()` → `run_mofa()` → `get_variance_explained()` → `get_factors()` → `get_weights()`.
- **Async**: when `future`/`promises` are installed, training runs via `shiny::ExtendedTask` + `promises::future_promise()` in a background worker so the app stays responsive; otherwise it falls back to a blocking synchronous call with an explicit "will be briefly unresponsive" notice.
- **Outputs**: a variance-explained-per-factor bar chart, a factor-scores scatter plot colorable by any metadata column, a factor×feature loadings panel (filterable by factor/sign/top-N), a cross-omics feature-correlation heatmap computed directly on the raw layers (not on MOFA2's own output), and an export bundle (CSV tables + a `reproducibility.md` recording the trained-at timestamp, factor count, seed, and installed MOFA2 version).

---

## 4. Data Flow

```text
Dataset Workspace (mod_multi_dataset.R)
  Reference/Example  ──┐
  Upload              ──┼──> raw$mats / raw$validations / raw$meta   (per-block, unmatched)
  GEO                 ──┘        │
                                 ▼
                      overlap()  — Sample Matching (intersect sample IDs)
                                 │
                                 ▼
                      proc$filtered_mats / proc$scaled_mats — Preprocessing
                                 │
                                 ▼
                      proc$batch_corrected  — Batch Diagnostics (optional)
                                 │
                                 ▼
                 "Activate" ──> multi_dataset$layers / $layer_meta / $sample_meta / $active=TRUE
                                 │
                 ┌───────────────┼──────────────────────────────────────────────┐
                 ▼               ▼                                              ▼
     Integrated Analysis   Cohort Harmonization                    Multi-omics Integration
     (MOFA2, embedded)     (mod_multi_overview.R)                  (DIABLO & SNF)
                                                                                │
                                                              ┌─────────────────┼─────────────────┐
                                                              ▼                 ▼                 ▼
                                                     SNF Clustering    Biomarker Discovery   (Compare /
                                                     (Patient          (DIABLO signature)     Sex-Stratified
                                                      Stratification)         │               sub-tabs)
                                                              │               │
                                                              │               ▼
                                                              │      multi_results$biomarker$df
                                                              │               │
                                                              └───────┬───────┘
                                                                      ▼
                                                          Gene–CpG Concordance
                                                       (candidate genes from DIABLO/
                                                        SNF/Joint, or custom lists)
                                                                      │
                                                                      ▼
                                                                 Pathways
                                                    (feature list from Biomarker Discovery /
                                                     SNF / Concordance / custom / upload)
                                                                      │
                                                                      ▼
                                                     Results Summary & Reproducibility
                                                (reads multi_results$<every sub-module above>)
```

*Every downstream sub-module also has an independent "Reference/Example Dataset" shortcut that bypasses this whole pipeline for a quick live recompute on one saved analysis cell — the diagram above shows the "Active Multi-Omics Dataset" path.

Two structural facts worth calling out explicitly:
- **This is not a strict linear pipeline.** Cohort Harmonization, Multi-omics Integration, and SNF Clustering are three independent consumers of the same `multi_dataset`, not stages of one chain — a user can run any of them in any order once the dataset is active. Biomarker Discovery, Concordance, and Pathways form a genuine dependency chain (each can consume the previous one's `multi_results$...` output), but each also has its own fallback (custom gene/CpG lists, or an upload) so it never strictly requires the upstream one to have been run first.
- **`multi_results` is the only hand-off channel between sub-modules.** No sub-module ever calls another sub-module's server function or reads its private reactiveValues — cross-module reuse happens only through this one shared, append-only list.

---

## 5. Every Implemented Sub-module

### 5.1 Cohort Harmonization

**Purpose.** A purely descriptive/diagnostic report on the active dataset: which modalities exist, which samples are genuinely shared across them, whether integration is even feasible for a given combination of modalities, and — as an explicit, separate, opt-in step — an honest held-out prediction-performance check against chance and single-omics baselines.

**Input.** `multi_dataset$layers`/`$sample_meta` (Active Multi-Omics Dataset), or a preloaded analysis cell's descriptor built from registry QC tables (no raw matrix in that case). User inputs: selected modalities, a minimum-overlap threshold, and (for Model Evaluation only) a phenotype/outcome column.

**Processing.** `ch_modality_descriptors()` builds one descriptor per layer; `ch_analysis_cells()` enumerates every feasible modality combination; `ch_integration_readiness()` classifies each combination as Ready (≥10 matched samples), Limited (≥3), Not Suitable (<3), or Single-modality; `ch_id_harmonization_table()` classifies every sample ID as Exact/Normalized match, Duplicate, Unmatched, or Ambiguous; `ch_pairwise_overlap_matrix()` builds the full N×N sample-overlap matrix. The separate Model Evaluation step runs `ch_evaluate_binary_outcome()`: `glmnet::cv.glmnet` (elastic net, α=0.5) per single-omics view and on an early-fused feature matrix, `k`-fold CV via `caret::createFolds`, scored with `pROC::roc`/`pROC::ci.auc`, and a DeLong test (`pROC::roc.test(..., method="delong")`) comparing the fused model against the best single-omics view.

**Output.** Modality-availability cards, a colored Matched/Partially-matched/Unmatched status badge with sentence summary, an Integration Readiness table, a batch/cohort summary bar chart, a data-completeness heatmap, a Sample Explorer (per-sample presence/absence + PCA highlight), and — after Model Evaluation — per-view and fused AUROC with 95% CI, a chance/majority-class baseline row, and a plain-language "Evidence of improvement / Weak evidence / No convincing evidence" conclusion driven purely by the DeLong p-value and effect size.

**Scientific importance.** Before any integrative method is trusted, a researcher needs to know the sample sizes actually being integrated are not an artifact of ID mismatches, and that a fused model's apparent advantage over a single omics layer is not just noise — this sub-module answers exactly those two questions honestly, refusing to report a number it cannot support (e.g., Model Evaluation is unavailable for the preloaded cohort entirely, since that path never carries a raw matrix).

**Code connection.** `mod_multi_overview.R` (UI/server) + `cohort_harmonization_helpers.R` (`ch_` prefix, all pure logic) + `cohort_harmonization_plots.R`.

---

### 5.2 Multi-omics Integration (DIABLO & SNF)

**Purpose.** A live, data-adaptive engine offering both a supervised (DIABLO) and an unsupervised (SNF) way to combine omics blocks, plus a head-to-head comparison against single-omics baselines and a sex-stratified variant.

**Input.** ≥2 omics blocks from the Active or Reference dataset; for DIABLO, a categorical outcome column with ≥2 classes of ≥3 samples each.

**Processing (DIABLO tab).** `mixOmics::tune.block.splsda()` grid-searches a data-sized `keepX` candidate grid (`5,10,20,50,100,150,200,300`, capped at each block's own feature count) scored by Balanced Error Rate, when auto-tuning is requested; `mixOmics::block.splsda()` fits the final sparse discriminant model; `mixOmics::perf()` cross-validates it (M-fold or leave-one-out, fold/repeat counts sized to the smallest class) and supplies BER, per-class error, AUC, and per-repeat feature-selection stability. Three "Advanced parameters" (tolerance, max-iterations, seed) are present in the UI but are never actually passed into these calls — see [§13](#13-known-gaps-decorative-controls-and-unimplemented-functionality).

**Processing (SNF tab).** `SNFtool::dist2()` (squared Euclidean, the only distance implemented) → `SNFtool::affinityMatrix(K, sigma=alpha)` per block → `SNFtool::SNF()` to fuse the per-block networks, with `T` optionally auto-selected by testing `{10,20,30,50}` and stopping once the fused network's correlation stabilizes. Cluster count is chosen by `SNFtool::estimateNumberOfClustersGivenGraph()`'s eigengap criterion; K/alpha auto-tuning additionally uses average silhouette width. Missing values are a hard blocker here (no imputation offered in this tab). A post-hoc Fisher's-exact/NMI/ARI check against a chosen outcome column is explicitly labeled post-hoc, never used to form the clusters.

**Processing (Compare tab).** Supervised half reuses Cohort Harmonization's own `ch_evaluate_binary_outcome()` for single-omics elastic-net baselines, matched to DIABLO's own fold count (independently-drawn folds, not a shared partition). Unsupervised half re-clusters each individual block's own affinity matrix and compares it to the fused clustering via NMI/ARI.

**Processing (Sex-Stratified tab).** A thin wrapper around the shared `multiomics_sexstratified_engine.R` (see [§9.4](#94-the-shared-sex-stratified-engine)) — a fixed-design, nested 5×5 CV pipeline with in-fold `limma` feature selection, offering either DIABLO or Random Forest as the classifier.

**Output.** BER/error-rate cards and error-bar plots, an AUC table, a selected-features panel plot, a sample score plot, a variance-explained plot, cross-block correlation of selected features, a feature-stability table; for SNF, the fused/per-block network heatmaps, cluster assignments, a cluster-quality estimate plot, and modality-contribution bars; for Compare, a bar chart of every model's performance side by side; for Sex-Stratified, per-stratum nested-CV AUROC with 95% CI and a wide biomarker-by-sex comparison table.

**Scientific importance.** DIABLO answers "which combination of genes and CpGs together best discriminate a known clinical outcome," while SNF answers the complementary, label-free question "do patients naturally group into molecularly distinct subtypes when both omics layers agree." Running both, and comparing each to what a single omics layer alone could achieve, is the central justification for a multi-omics study design over analyzing each layer in isolation.

**Code connection.** `mod_multi_integration.R` + `multiomics_integration_live_helpers.R` (`mi_` prefix) + `multiomics_integration_live_plots.R`.

---

### 5.3 SNF Clustering (Patient Stratification)

**Purpose.** A deeper, dedicated unsupervised patient-stratification workflow built on top of the same SNF engine Integration's SNF tab uses, adding what that tab does not: a genuine single-omics fallback, explicit per-block preprocessing, resampling-based cluster stability, parameter-sensitivity analysis, and formal clinical-variable association testing.

**Input.** One or more omics blocks (single-omics mode is explicitly supported and labeled as such, never presented as multi-omics integration); optional clinical/survival metadata for post-hoc association.

**Processing.** Per-block preprocessing (`multi_live_handle_missing()`/`multi_live_normalize()`/`multi_live_filter_features()`, reused from the Dataset Workspace's own live helpers) before clustering; the same `mi_snf_*()` affinity/fusion/clustering primitives Integration uses; a **required** stability check bundled into every run (20 resamples at 80% subsampling by default, scored by Adjusted Rand Index against the full-cohort clustering, verdict thresholds "Stable" ≥0.75 mean ARI, "Moderately stable" ≥0.5, else "Unstable" — the code's own design principle is that "clusters must never be shown/labeled without stability evidence"); a parameter-sensitivity re-run at the feasible low/high end of K, alpha, and T; per-feature Kruskal-Wallis ranking of which features actually distinguish the discovered clusters.

**Output.** A real spectral embedding of the fused network (the same eigenvectors spectral clustering itself partitions, not a generic PCA), a top-variance feature heatmap ordered by cluster, the fused/per-block network heatmaps, cluster-quality and stability verdicts, categorical (Fisher's exact + Cramér's V), continuous (Kruskal-Wallis), and survival (`survival::survfit`/`survdiff`/`coxph`, with a Kaplan-Meier plot when `survminer` is installed) cluster-association tests, and a feature-ranking table.

**Scientific importance.** A cluster that only survives on the full dataset and vanishes under resampling, or that only appears at one narrow (K, alpha) setting, is not a real patient subtype — it is noise. Requiring stability and sensitivity evidence before showing any cluster assignment is what separates exploratory patient stratification from an overclaimed one.

**Code connection.** `mod_multi_stratification.R` + `snf_clustering_helpers.R` (`sfc_` prefix, reuses `mi_snf_*()` from Integration rather than reimplementing fusion/clustering) + `snf_clustering_plots.R`.

---

### 5.4 Biomarker Discovery

**Purpose.** Fits one DIABLO model against a user-chosen outcome across exactly two role-assigned blocks ("Transcriptomics" and "Methylomics"), then reports the selected features as a ranked, stability-labeled candidate biomarker signature.

**Input.** Two omics blocks from the Active or Reference dataset; a categorical outcome column; a variance-based unsupervised pre-filter per block; DIABLO model settings (ncomp, manual or auto-tuned keepX, a single cross-block design weight, CV method/folds/repeats/distance metric).

**Processing.** Calls the identical `mi_diablo_run()` engine Multi-omics Integration uses (same `tune.block.splsda`/`block.splsda`/`perf()` calls) — this sub-module adds no separate feature-selection algorithm of its own. Ranking within the resulting Signature table is purely by `|loading|` (the DIABLO sparse-loading weight from `mixOmics::selectVar()`), computed within each block × component — there is no combined score involving effect size, p-value, or selection frequency. "Stability" is the fraction of the configured CV's repeats in which `perf()` re-selected a given feature (fixed thresholds: ≥80% "Stable", ≥50% "Moderately stable", else "Low stability"), not an independent bootstrap.

**Output.** A ranked Signature table (CSV), cross-validated Performance (BER, per-class error, AUC, plus a separately-computed pooled out-of-fold ROC for binary outcomes), a Stability table/plot, an Integration tab reporting the fitted cross-block design weight and component correlation, sample-score/loadings/heatmap/variance/ROC plots, and — via the same shared sex-stratified engine as Multi-omics Integration — a Sex-Stratified nested-CV comparison. A Circos-style relationship plot is explicitly stated as "not implemented" rather than approximated.

**Scientific importance.** This is where the module produces its actual candidate list of jointly-expression-and-methylation-relevant genes/CpGs for a specific clinical question, with an explicit, disclosed ranking rule and an honest stability label rather than a single point-in-time selection.

**Code connection.** `mod_multi_biomarker.R` + `multiomics_biomarker_helpers.R` + `multiomics_biomarker_plots.R`; reuses `mi_diablo_*()` (Integration's engine) and `mss_run_stratified()` (the shared sex-stratified engine).

---

### 5.5 Gene–CpG Concordance

**Purpose.** For a set of candidate genes — drawn from Biomarker Discovery's DIABLO signature, SNF cluster-associated features, or a custom list — maps each gene to its annotated CpG probes and asks whether the expression change and methylation change at that gene move together in the direction expected from the CpG's genomic location.

**Input.** A candidate-gene pool (`mcc_candidate_pool()`, pulling from `multi_results$biomarker`/`$integration_stratified`, live SNF cluster-feature ranking, or custom entries), a methylation-array annotation choice (450K/EPIC), region/threshold filters, and — on the live path — expression and methylation layers plus a 2-class design column.

**Processing.** Gene↔CpG mapping is **purely Illumina-manifest annotation based** — the first token of each CpG's `UCSC_RefGene_Name`/`UCSC_RefGene_Group` fields; there is no genomic-distance window anywhere in the code (`tss_distance` is reported as "Not available," never estimated). Expression change is a plain group-mean log2FC (Welch `t.test`); methylation change is computed on M-values for the significance test (delta-M, stored in a column literally named `dbeta` — a naming holdover worth noting explicitly since it is *not* the delta-beta value) with delta-beta reported alongside purely for interpretability, never mixed with delta-M in the statistical test. Direction is classified against a genomic-region rule: at a promoter/TSS region, "canonical" means an *inverse* relationship (hypermethylation with down-regulation, or the reverse); at the gene body, "canonical" means a *concordant* relationship — and a pair is only ever labeled canonical/non-canonical once it clears both omics' significance thresholds; otherwise it is explicitly "Not applicable," never guessed. Sample-level correlation (Pearson/Spearman, BH-FDR) is computed only when ≥3 matched samples with both expression and methylation values exist for that gene/CpG pair — this is only possible on the live path, since the preloaded path never carries per-sample matrices. **No mediation analysis exists anywhere in the codebase** (confirmed by a direct search; the code's own limitations list states this explicitly).

**Output.** The full gene–CpG pairs table, direction/canonical cross-tabs (both a genome-wide-nominal-p-gated view and this run's own threshold-gated view), a genomic-location scatter, a candidate-biomarker table filtered to the highest-evidence rows, a concordance scatter/direction-quadrant/multi-omics-evidence heatmap/gene–CpG network, a single-pair correlation plot (live path only), and a disclosed, component-visible priority score — never a black-box number, and never labeled "confirmed."

**Scientific importance.** DNA methylation's regulatory effect on a gene is directional and location-dependent (repressive near a promoter, often permissive within a gene body); testing for exactly that expected direction — rather than any correlation — is what turns a coincidental expression/methylation pair into biologically interpretable evidence.

**Code connection.** `mod_multi_concordance.R` + `multiomics_concordance_live_helpers.R` (`mcc_` prefix) + `multiomics_concordance_plots.R`; the annotation lookup itself is borrowed from the Cross-Omics module's `cx_get_region_annotation()`/`cx_classify()`.

---

### 5.6 Pathways

**Purpose.** Takes whatever candidate feature list the earlier sub-modules produced (or a user upload) and tests it for enrichment in curated biological pathway/gene-set databases, using both over-representation (ORA) and rank-based (GSEA) methods.

**Input.** A candidate list from `mcc_candidate_pool()` (the same pool Concordance uses — DIABLO/SNF/Joint/custom), each feature optionally carrying real effect sizes joined in from the precomputed concordance registry table; or a user-uploaded table with auto-detected, user-confirmed column roles. A background/universe choice, a database selection, and (for GSEA) a ranking statistic.

**Processing.** Real calls to `clusterProfiler::enrichGO/gseGO`, `clusterProfiler::enrichKEGG/gseKEGG`, `ReactomePA::enrichPathway/gsePathway`, and `clusterProfiler::enricher/GSEA` with `msigdbr`-sourced term sets for WikiPathways and the MSigDB Hallmark collection — every one of these is gated by a real `requireNamespace()` check and surfaced as unavailable (not silently substituted) if the package is missing. All are run with Benjamini-Hochberg FDR correction. **The background/universe default is "Entire selected database"** (no experimental universe supplied at all) rather than the measured-feature background — the UI does display an explicit warning when this default is in effect, and a genuine measured-feature-background option ("Auto — measured features in active dataset") exists and correctly restricts the universe to the actual expression/methylation layer's own features, but a user must deliberately choose it. Every tested term is retained in the output table with a `significant` flag column rather than being filtered out below the FDR threshold, so a "no pathways meet the threshold" state is never confused with "no pathways were tested." A "convergence" view exists as a per-pathway `integration_label` (RNA + Methylation supported / Transcriptomics-only / Methylomics-only) rather than a separate tab, alongside a directional-concordance check between the expression and methylation evidence for genes in that pathway.

**Output.** A dot plot, bar plot, a pathway×omics evidence heatmap, a gene–pathway network, a full enrichment DT table (CSV download); real KEGG diagrams overlaid with per-gene expression fold-changes via `pathview::pathview()` and real Reactome diagrams fetched live from Reactome's own content service (both fail with an explicit error rather than a placeholder image if the network/package is unavailable); a per-pathway gene table; and a reproducibility metadata table (database, method, species, background choice, ranking method, feature counts, mapping rate, FDR threshold, gene-set size bounds, timestamp).

**Scientific importance.** A biomarker list of a dozen genes is hard to interpret directly; asking whether those genes cluster into a known biological process, and whether both the expression and methylation evidence for that process point the same direction, is what converts a statistical candidate list into a testable biological hypothesis.

**Code connection.** `mod_multi_pathway.R` + `multiomics_pathway_helpers.R` + `multiomics_pathway_plots.R`.

---

### 5.7 Results Summary & Reproducibility

**Purpose.** A session-level rollup — what has been computed this session, what this module can and cannot do, and a downloadable bundle of everything loaded so far — with zero computation of its own.

**Input.** `multi_results` (read-only, from every other sub-module).

**Processing.** None — every number shown is a direct read of a value another sub-module already published (e.g., "Enriched pathway terms" is literally `nrow(multi_results$pathway$df)`), or "Not loaded" if that sub-module has not run yet in the session.

**Output.** A session dashboard of summary cards; a hardcoded, explicitly-maintained "Supported vs. not implemented" list (`MULTI_KNOWN_LIMITATIONS`, see [§13](#13-known-gaps-decorative-controls-and-unimplemented-functionality) for its full current content); an installed-package-versions table (real `utils::packageVersion()` lookups, never a guessed string); and a downloadable ZIP bundle containing one CSV per loaded sub-module's result table, plus a plain-text `report.md` listing what was loaded this session and pointing to the exact offline pipeline scripts (`MULTI_REPRODUCIBILITY_SCRIPTS`) that produced whichever precomputed tables are in view.

**Scientific importance.** A researcher revisiting results later — or writing a methods section — needs an honest, single place that states exactly what was run in this session, what the software is and is not capable of, and where the numbers on a precomputed tab actually came from; fabricating none of this is the entire point of the tab.

**Code connection.** `mod_multi_summary.R` (no separate helpers file; the `MULTI_KNOWN_LIMITATIONS`/`MULTI_REPRODUCIBILITY_SCRIPTS`/`multi_build_report()`/`multi_package_versions()` functions it reads live in the shared `multiomics_helpers.R`).

---

## 6. Scientific Importance

Analyzing gene expression and DNA methylation as two independent, unconnected datasets discards exactly the information a multi-omics design is meant to capture: DNA methylation is one of several mechanisms that regulate whether and how strongly a gene is transcribed, so a gene's expression change and the methylation change at its own regulatory region are not independent observations — they are two readouts of a plausibly related regulatory event. The Multi-Omics module operationalizes this in several concrete, implemented ways rather than as an abstract claim:

- **Complementary information and molecular regulation.** Gene–CpG Concordance does not merely correlate two numbers; it tests whether the *direction* of the relationship matches what is known about methylation's regulatory role at a promoter versus a gene body, which is a stronger and more falsifiable claim than a generic correlation.
- **Cross-omics relationships as a discriminative signal.** DIABLO explicitly models the relationship *between* blocks (a user-set or auto-tuned design-weight matrix) while selecting features from each block that jointly discriminate an outcome — a feature is not selected for its individual power alone, but for its contribution alongside the other omics layer.
- **Patient heterogeneity.** SNF (in both the Integration and dedicated SNF Clustering sub-modules) fuses similarity structure from every included omics layer before clustering, so a patient subgroup only emerges if it is supported by more than one molecular layer's own notion of similarity — and, in SNF Clustering specifically, that subgroup is only reported once it survives resampling and parameter-sensitivity checks, directly addressing the risk that an unsupervised cluster is an artifact of one omics layer's particular noise structure.
- **Biomarker robustness.** Biomarker Discovery's stability metric (selection frequency across repeated cross-validation folds) and the Sex-Stratified engine's fully nested cross-validation (feature selection performed strictly inside each training fold, never touching held-out samples) are both direct, implemented defenses against the well-known failure mode of reporting an inflated, overfit biomarker signature.
- **Biological interpretation and pathway-level integration.** Pathways closes the loop from an individual candidate feature list back to biological meaning, and — via the convergence/directional-concordance check — asks whether the transcriptomic and methylation evidence for a pathway agree, rather than treating each omics layer's pathway result as if it were independent confirmation.
- **Integration of molecular signals into an honest evaluation.** Cohort Harmonization's Model Evaluation and Multi-omics Integration's Compare tab both exist specifically to answer the question a multi-omics analysis is obligated to answer but is often skipped: does the fused model actually outperform the best single-omics model by more than chance, tested with a real statistical comparison (DeLong's test), rather than merely reporting the fused model's own number in isolation.

None of this constitutes clinical validation. Every AUROC reported on the bundled cohort's precomputed tabs is the offline pipeline's own cross-validated (not externally validated) performance, and the module's own `MULTI_KNOWN_LIMITATIONS` list states this outright. The scientific value demonstrated here is methodological — a working, internally honest multi-omics integration workflow — not a diagnostic claim.

---

## 7. End-to-End Code Walkthrough

### Step 1 — Dataset enters the application
A user picks a `dataset_source` radio (`mod_multi_dataset.R`): Reference/Example, Upload, or GEO. For an upload, `fileInput` triggers a preview read (`multi_live_read_matrix()`/`multi_live_pivot_long()`) and structural checks (`multi_live_detect_orientation()`, `multi_live_detect_table_shape()`, `multi_live_detect_omics_type()`) that populate advisory UI notes but do not yet commit anything.

### Step 2 — Data are stored
Clicking "Validate Datasets" (or "Load Reference Dataset") runs an `observeEvent` that writes surviving matrices into `raw <- reactiveValues(mats=list(), validations=list(), labels=list(), provenance=list(), meta=NULL)` — a module-local reactiveValues object, distinct from the shared `multi_dataset` that will only be written to at the very end of the pipeline.

### Step 3 — Data are validated
`multi_live_validate_matrix()` computes per-layer QC facts (missingness, duplicate IDs, zero-variance columns) with no automatic correction; `multi_dataset_status()` turns those facts into a Ready/Review-Required/Not-Compatible verdict per layer using fixed, disclosed thresholds (e.g. >20% missing → Review Required; <3 samples → Not Compatible).

### Step 4 — Samples are harmonized
The Sample Matching tab's `overlap <- reactive({...})` intersects sample IDs across `raw$mats` (optionally after a patient-ID or mapping-file translation), producing the matched-ID set every subsequent step subsets to. Metadata from multiple sources (per-dataset files, a shared metadata file, or a GEO series' own sample metadata) is merged via `mo_merge_sample_meta()`.

### Step 5 — Data are prepared
Preprocessing (`observeEvent(input$preprocess_btn)`) hard-gates on ≥3 matched samples, then applies missing-value handling, omics-appropriate normalization, and optional variance filtering — writing `proc$filtered_mats`/`proc$scaled_mats`. Batch Diagnostics optionally layers `sva::ComBat`/`limma::removeBatchEffect` on top, writing `proc$batch_corrected`, gated by a confounding check that blocks correction unless explicitly overridden.

### Step 6 — Multi-omics integration occurs
Once "Activate" has written `multi_dataset$layers`/`$sample_meta`, Multi-omics Integration's DIABLO tab calls `mi_diablo_run()` (`mixOmics::tune.block.splsda` → `block.splsda` → `perf`), and its SNF tab calls `mi_snf_run()`-style logic (`SNFtool::dist2` → `affinityMatrix` → `SNF` → `estimateNumberOfClustersGivenGraph`/`spectralClustering`). Both are dispatched inside `shiny::ExtendedTask`/`future::future_promise()` when available, so the long-running fit does not block the Shiny session.

### Step 7 — Downstream analyses occur
Biomarker Discovery calls the same `mi_diablo_run()` engine on its own two role-assigned blocks; SNF Clustering calls the same `mi_snf_*()` primitives with added preprocessing/stability/sensitivity/clinical-testing logic; Gene–CpG Concordance reads `multi_results$biomarker$df` (or live SNF cluster-feature ranking) as its candidate-gene pool; Pathways reads the same candidate pool (via `mcc_candidate_pool()`, the identical function Concordance uses) as its feature list. Each of these is a plain, one-directional read of a previously-published `multi_results$...` slot — never a direct function call into another sub-module.

### Step 8 — Results are generated
Every sub-module's "Run" `observeEvent` builds one canonical result list (e.g. Concordance's `list(ok, error, pairs_df, expr_mat, meth_mat, ..., settings_snapshot)`) and stores it in a module-local `reactiveValues` (commonly named `state$result`). This list is the single source every downstream `renderUI`/`renderPlot`/`renderDataTable`/`downloadHandler` in that sub-module reads from — nothing recomputes independently per output.

### Step 9 — Results are rendered
`output$<name>_ui <- renderUI({...})` builds each tab's content, always checking the result's `ok` flag first and returning an explicit "not run yet" or "failed: <reason>" state otherwise (via a shared `multi_empty_state()`/`multi_plot_or_empty()` helper pair). `output$<name>_table <- DT::renderDataTable({...})` renders result tables; `output$<name>_plot <- renderPlot({...})` renders each ggplot2 plot; `downloadHandler()` supplies every CSV/PNG/ZIP export.

### Step 10 — User interpretation
The final outputs let a researcher conclude, with explicit numeric support rather than a qualitative impression: which modality combinations are even matchable in this cohort; whether a supervised (DIABLO) or unsupervised (SNF/MOFA2) integration finds real, cross-validated or stability-checked structure; which genes and CpGs jointly and directionally implicate a clinical outcome; and which biological pathways that candidate list falls into — while the Results Summary tab keeps an honest running record of exactly what was computed to reach those conclusions, and what the software explicitly does not yet do.

---

## 8. Important R/Shiny Concepts Used in This Module

- **`reactive()`** — a lazily-evaluated, cached expression that only re-runs when one of the reactive values it reads changes. Used throughout for cheap derived state, e.g. `overlap <- reactive({...})` (Dataset Workspace) recomputes the matched-sample set only when the raw matrices or matching method change, not on every input tick.
- **`eventReactive(trigger, {...})`** — like `reactive()`, but only recomputes when `trigger` fires, ignoring changes to anything else it reads inside the body. Cohort Harmonization's `harmonization_result <- eventReactive(input$analyze_btn, {...})` is the canonical example: filter/threshold inputs can be freely adjusted without recomputation until "Analyze Cohort" is clicked again — this is the mechanism behind the module-wide rule that expensive analyses never run automatically.
- **`observeEvent(trigger, {...})`** — runs side effects (not a return value) exactly when `trigger` fires. Every "Run"/"Activate"/"Apply correction" button in this module is wired this way, and it is also where `showNotification()` calls and writes to shared state (`multi_dataset$...`, `multi_results$...`) happen.
- **`observe({...})`** — like `observeEvent` but re-runs whenever *any* reactive value it reads changes, with no single named trigger. Used sparingly here, mainly for the one-line "publish this sub-module's finished result to `multi_results$<id>`" blocks that should re-fire whenever the result object itself changes.
- **`renderUI({...})` / `uiOutput()`** — builds a piece of the UI on the server side and reactively re-renders it. Used pervasively so that a tab's very structure (which controls even appear) can depend on what data is actually available — e.g. a phenotype-column selector only appears if one was actually detected.
- **`renderPlot({...})` / `renderDataTable({...})`** — render a ggplot2 object or a `DT::datatable()` into a named output slot. This module's shared convention (`multi_plot_or_empty()`) wraps every plot-producing `renderUI` so that a `NULL` result renders an explicit empty-state message instead of a blank plot area.
- **`req(x)`** — silently stops execution of the current reactive/output if `x` is falsy/`NULL`/empty, without raising a visible error. Used everywhere a downstream computation is only meaningful once a prior selection exists (e.g. `req(input$batch_layer)` before reading `proc$scaled_mats[[input$batch_layer]]`).
- **`validate(need(condition, message))`** — like `req()`, but renders `message` into the enclosing output instead of silently stopping. Used inside `renderUI`/`renderPlot`/`observeEvent` blocks that have a natural output to render the message into; a real bug fixed elsewhere in this codebase (Dataset Workspace's batch-correction blocker) illustrates the limit of this pattern: `validate()` called inside a plain `observeEvent` (which has no output context) fails completely silently, which is why that particular check was rewritten as an explicit `showNotification(..., type="error"); return()`.
- **`downloadHandler(filename, content)`** — the mechanism behind every CSV/PNG/ZIP export in the module; `content` is a function that writes the file to the path Shiny gives it, evaluated freshly on each download click against whatever the current reactive state is.
- **`shiny::ExtendedTask` + `promises::future_promise()`** — the module's asynchronous-execution pattern for genuinely slow computations (DIABLO tuning, SNF K/alpha/T search, MOFA2 training, the Sex-Stratified nested CV): the heavy work runs in a background worker process so the Shiny session's main event loop is never blocked, with a synchronous fallback (and an explicit "will be briefly unresponsive" notice) when the `future`/`promises` packages are not installed.

---

## 9. Integration Methods in Detail

### 9.1 DIABLO — supervised

- **Required input**: ≥2 omics blocks, matched samples, and a categorical outcome with ≥2 classes of ≥3 samples each.
- **Data preparation**: an unsupervised variance pre-filter is optionally applied per block before fitting; blocks otherwise enter the model as-is (already normalized/scaled by the Dataset Workspace).
- **Parameters**: number of components (`ncomp`, feasible range derived from class count and smallest-class size), `keepX` per block per component (manual, or auto-tuned via a data-sized grid search), a cross-block design-weight matrix (a single slider in the main tabs; a correlation-derived fixed 2×2 matrix in the Sex-Stratified engine), CV method (M-fold, capped at the smallest class size, or leave-one-out below a sample-count ceiling), folds, repeats, and a prediction distance metric.
- **Computational process**: `mixOmics::tune.block.splsda()` (if auto-tuning) → `mixOmics::block.splsda()` → `mixOmics::perf()`.
- **Result**: cross-validated BER/per-class error/AUC, selected features with loadings, and a per-repeat feature-selection stability frequency.
- **Visualization**: sample-score plots (per component, colored by class), a per-block selected-feature loadings panel, an error-bar plot across components, a variance-explained plot.
- **Biological interpretation**: the selected features are those whose combined, cross-block pattern best discriminates the chosen clinical outcome — a supervised, hypothesis-testing use of multi-omics data.

### 9.2 SNF — unsupervised

- **Required input**: ≥1 omics block with no missing values (in the Integration tab; SNF Clustering additionally offers explicit preprocessing to remove missingness first).
- **Data preparation**: optional per-block standardization.
- **Parameters**: K (neighborhood size) and alpha/sigma (kernel width), both with data-sized feasible ranges and an optional auto-tune (grid search scored by average silhouette width); T (fusion iterations), fixed or auto-converged; cluster count, chosen by eigengap or set manually; a clustering algorithm choice (spectral, hierarchical, or PAM) in SNF Clustering specifically.
- **Computational process**: `SNFtool::dist2()` → `SNFtool::affinityMatrix()` per block → `SNFtool::SNF()` to fuse → `SNFtool::estimateNumberOfClustersGivenGraph()` → `SNFtool::spectralClustering()` (or hierarchical/PAM).
- **Result**: a fused similarity network, cluster assignments, and (in SNF Clustering) a stability verdict and a parameter-sensitivity verdict.
- **Visualization**: per-block and fused network heatmaps, a spectral embedding scatter of patients, a top-variance feature heatmap ordered by cluster, a modality-contribution (NMI-based) bar chart.
- **Biological interpretation**: clusters represent patient subgroups whose combined molecular similarity — not a single omics layer, and not any known outcome label — places them closer together than other patients; a real cluster-outcome association is only ever tested post-hoc, never used to construct the clusters.

### 9.3 MOFA2 — unsupervised, factor-based

- **Required input**: ≥2 omics layers from the activated dataset (no outcome column needed at all).
- **Data preparation**: none beyond what the Dataset Workspace already did; missing values are handled natively by MOFA2, not pre-imputed by this app's code specifically for MOFA2.
- **Parameters**: number of factors (auto-capped by sample size), seed, convergence-speed preset; everything else at MOFA2's own package defaults.
- **Computational process**: `MOFA2::create_mofa()` → `prepare_mofa()` → `run_mofa()`.
- **Result**: a set of latent factors, each with a variance-explained-per-view breakdown and per-feature loadings.
- **Visualization**: variance-explained bar chart, factor-score scatter (colorable by metadata), factor×feature loadings panel, factor heatmap.
- **Biological interpretation**: each factor is a continuous axis of coordinated variation spanning multiple omics layers — a softer, dimensionality-reduction-style alternative to the discrete groupings SNF produces, useful for finding a shared axis of variation (e.g. a disease-severity gradient) that manifests in both expression and methylation without needing a predefined class label.

### 9.4 The shared sex-stratified engine

`multiomics_sexstratified_engine.R` is a parameterized reproduction of the offline research pipeline's own sex-stratified nested-CV scripts, called identically by Multi-omics Integration's Sex-Stratified tab and Biomarker Discovery's Sex-Stratified tab. It runs a genuine 5×5 repeated cross-validation (`caret::createFolds`, five outer folds, five repeats, deterministically seeded per repeat) with feature selection performed strictly inside each training fold (`limma::eBayes`/`lmFit`/`topTable` on a covariate-adjustable design, never touching held-out rows), offering either DIABLO or Random Forest as the per-fold classifier. `keepX` (for DIABLO) is deterministically clamped to a fixed range rather than tuned by any inner grid search — so this engine is leakage-safe with respect to feature selection, but not a nested hyperparameter search. A separate, clearly-labeled full-cohort refit (outside the CV loop) produces the descriptive DIABLO-loadings/Random-Forest-importance panel shown alongside the cross-validated performance numbers.

---

## 10. Visualization and Output Logic

Every plot in the Multi-Omics module is built once as a `ggplot2` object using the app-wide `theme_arthomix()`/`ARTHOMIX_COLORS` palette (`multiomics_plots.R` and its per-sub-module `*_plots.R` counterparts), and that same object is reused for both on-screen rendering and PNG export — there is never a second, separately-coded version of a plot for its download button.

Three small shared helpers enforce the module's "never show a fabricated or blank result" rule everywhere:
- **`multi_empty_state(msg)`** — a single styled "not run / no data" panel used identically across every sub-module.
- **`multi_plot_or_empty(plot_fn, output_id, msg, height)`** — evaluates `plot_fn()` once up front; if it returns `NULL` (e.g. a filter zeroed out every row), it renders `multi_empty_state(msg)` instead of an empty `plotOutput`, and logs the real error server-side (not silently swallowed) if `plot_fn()` actually threw, so a genuine bug is distinguishable from a legitimate "no rows for this selection."
- **`multi_png_download(plot_fn, filename_fn)`** — a generic `downloadHandler` that re-evaluates `plot_fn()` at download time and calls `ggplot2::ggsave()`, so the exported PNG always matches whatever is currently on screen.

Tables are rendered with `DT::datatable(..., class="stripe hover compact")` throughout, with `downloadHandler`-backed CSV export next to almost every table. Several sub-modules additionally bundle multiple result CSVs plus a plain-text reproducibility note into a single ZIP download (Results Summary's session bundle, Biomarker Discovery's parameter/sample/performance bundle, MOFA2's export bundle).

---

## 11. Input → Processing → Output Summary Table

| Sub-module | Input | Main Processing | Output | Scientific Purpose |
|---|---|---|---|---|
| **Dataset Workspace** | Reference cell / uploaded files / GEO accession | Validation, sample matching, normalization, optional batch correction, compatibility check | Activated `multi_dataset$layers`/`$sample_meta` | Establish one usable, quality-checked multi-omics dataset for every downstream sub-module |
| **Integrated Analysis (MOFA2)** | ≥2 activated omics layers | `MOFA2::create_mofa/run_mofa` unsupervised factor model | Latent factors, variance-explained, feature loadings | Find shared axes of variation spanning omics layers without a predefined outcome |
| **Cohort Harmonization** | Active/preloaded dataset; optional outcome column | Modality/sample-overlap descriptors, integration-readiness classification, elastic-net held-out evaluation vs. chance/single-omics + DeLong test | Readiness table, matched-status badge, honest evaluation verdict | Confirm the dataset is genuinely usable before trusting any integrative result |
| **Multi-omics Integration** | ≥2 blocks; outcome column (DIABLO) | `mixOmics::block.splsda`/`tune.block.splsda`/`perf`; `SNFtool` affinity/fusion/clustering; shared sex-stratified engine | BER/AUC, selected features, fused network, clusters, single-omics comparison | Fuse omics blocks both supervised and unsupervised, and quantify the benefit of fusing them |
| **SNF Clustering** | ≥1 block; optional clinical metadata | Preprocessing + `SNFtool` fusion/clustering + resampling stability + parameter sensitivity + clinical association tests | Stability-verified cluster assignments, clinical association p-values | Discover, and rigorously validate, molecularly-defined patient subgroups |
| **Biomarker Discovery** | 2 role-assigned blocks; outcome column | `mi_diablo_run()` (DIABLO), loading-magnitude ranking, repeated-CV selection-frequency stability | Ranked signature table, performance metrics, stability labels | Produce a disclosed, stability-labeled joint biomarker candidate list |
| **Gene–CpG Concordance** | Candidate genes (DIABLO/SNF/custom); methylation array annotation | Illumina-annotation gene↔CpG mapping, expression/methylation group tests, region-aware direction classification, sample-level correlation | Direction-classified pairs table, priority score, plots | Test whether expression and methylation changes at a gene are directionally consistent with known regulatory biology |
| **Pathways** | Candidate feature list (from Biomarker/SNF/Concordance/custom/upload) | `clusterProfiler`/`ReactomePA`/`msigdbr`-backed ORA and GSEA, BH-FDR | Enrichment table, pathway maps, convergence/concordance labels | Translate a candidate feature list into interpretable biological processes |
| **Results Summary & Reproducibility** | `multi_results` from every other sub-module | Read-only rollup, no computation | Session dashboard, limitations list, package versions, downloadable bundle | Provide an honest, reproducible record of exactly what this session computed |

---

## 12. End-to-End Workflow Summary

Reflecting the actual implementation order (not a forced idealized sequence):

**Dataset acquisition and activation** (Dataset Workspace: source selection → validate → sample-match → preprocess → batch-diagnose → activate) →
**three independent, order-agnostic consumers of the activated dataset** (Integrated Analysis/MOFA2, embedded in the Dataset Workspace; Cohort Harmonization; Multi-omics Integration, which itself branches into DIABLO, SNF, Compare, and Sex-Stratified) →
**a genuine dependency chain among the remaining three sub-modules** (Biomarker Discovery's DIABLO signature, or SNF Clustering's cluster-associated features → Gene–CpG Concordance's candidate pool → Pathways' feature list; each stage also accepts a custom list or upload as a fallback, so the chain is a convenience, not a hard requirement) →
**Results Summary & Reproducibility**, a final, purely read-only rollup of whatever combination of the above was actually run in the session.

The one respect in which the implementation genuinely does *not* match a naive "linear pipeline" mental model: Cohort Harmonization, Multi-omics Integration, and SNF Clustering do not feed into each other in code — they are three parallel lenses on the same activated dataset, and a user is free to run any subset of them in any order.

---

## 13. Known Gaps, Decorative Controls, and Unimplemented Functionality

This section exists because the documentation standard for this module requires distinguishing implemented functionality from a UI control that merely exists. Every item below was verified directly against the code, not inferred.

- ~~**The "Retrieve from GEO" data source is unreachable.**~~ **Fixed.** `mo_block_ui()`'s GEO branch was missing the `actionButton` that `observeEvent(input[[paste0(gbid, "_geo_fetch")]], ...)` (`mod_multi_dataset.R:624`) was waiting on — a "Fetch from GEO" button has since been added next to the GEO-accession field, and the full retrieval path (`multi_geo_layer_fetch()` → `GEOquery::getGEO()` → `multi_geo_platform_matrix()`) has been verified end-to-end against a real GEO series (GSE1000: fetched, platform GPL96, 10 samples × 12,548 features after gene-symbol collapsing) both as a direct function call and through the running app's UI.
- **DIABLO's "Advanced parameters" (tolerance, max-iterations, seed) are decorative.** These three inputs exist in the Multi-omics Integration DIABLO tab but are never read by `mi_diablo_run()` — the underlying `mixOmics` calls use the package's own internal defaults and no explicit seed.
- **`input$preloaded_pick`** (the "Select a reference dataset" dropdown, which currently has only one choice) **is never read by any server code** — only the analysis-cell selector and the data-source radio actually matter for the Reference/Example path.
- **A batch-correction layer-mislabeling risk exists** in the Dataset Workspace: `proc$batch_corrected` is a single slot (not keyed by layer name), so if a user corrects layer A, then changes the "Dataset" selector to layer B without re-running correction, and clicks "Activate," the activation logic can attach A's corrected matrix under B's name (the on-screen "success" panel correctly avoids this confusion by checking the layer match itself, but the final `final_mats` assembly at Activate time does not perform the same check).
- **Cohort Harmonization's Model Evaluation "Chance / majority-class baseline" row displays a hardcoded 0.5** rather than the actual computed majority-class proportion (`res$majority_baseline`), which is computed but not surfaced in that particular table.
- **Gene–CpG Concordance's `dbeta` column is semantically the delta-**M** value**, not delta-beta — the true delta-beta is reported separately in a `delta_beta` column. This is a naming holdover from the precomputed table's own `delta_M` field and is worth flagging to any reader relying on column names alone.
- **No genomic-distance window exists for gene–CpG mapping** — mapping is purely Illumina-manifest annotation based, and `tss_distance` is always reported as "Not available."
- **No mediation analysis exists anywhere in the codebase.**
- **Pathways' background/universe default is "Entire selected database,"** not the measured/tested feature set — a correctly-implemented "Auto — measured features" option exists but is not the default, and the UI does display an explicit warning when the whole-database default is active.
- **No cell-type composition (EpiDISH/CIBERSORTx), Covariates & Clinical Metadata, bootstrap/permutation/external-cohort validation, unified model-benchmarking, XGBoost/Boruta/SHAP, formal sex×disease interaction model, sex-chromosome biology annotation, or Druggable Target Linker sub-module exists.** This list matches, verbatim in content, the module's own self-reported `MULTI_KNOWN_LIMITATIONS` (`multiomics_helpers.R`), which is read live into the Results Summary tab rather than being a separate claim made only in this document.
- **Reported AUROCs on the precomputed-cohort tabs are the offline pipeline's own cross-validated performance, not performance on an independent replication cohort.**

---

## 14. Thesis-Level Main Paragraph

The Multi-Omics module of ArthOMix implements a complete, code-verifiable workflow for integrating matched transcriptomic and DNA-methylation data from a rheumatoid-arthritis anti-TNF treatment cohort, structured around a single shared dataset object that a Dataset Workspace populates through explicit validation, sample-matching, normalization, and optional batch-correction steps before activation. Once activated, the dataset feeds three complementary integration strategies — a supervised discriminant method (DIABLO, via `mixOmics::block.splsda`), an unsupervised similarity-fusion method (Similarity Network Fusion, via `SNFtool`), and an unsupervised latent-factor method (MOFA2) — each evaluated not in isolation but against single-omics baselines and, for SNF-derived patient clusters, against resampling-based stability and parameter-sensitivity checks that guard against reporting an unstable or overfit result. Downstream, a Biomarker Discovery sub-module surfaces a ranked, stability-labeled joint biomarker signature from the fitted DIABLO model, a Gene–CpG Concordance sub-module tests whether the expression and methylation changes at each candidate gene are directionally consistent with known regulatory biology at that gene's promoter or gene body, and a Pathways sub-module translates the resulting candidate list into enriched biological processes using standard, correctly cited enrichment methodology (over-representation and gene-set enrichment analysis against GO, KEGG, Reactome, WikiPathways, and Hallmark gene sets). A dedicated Cohort Harmonization sub-module and a session-level Results Summary tab together enforce the module's central methodological discipline: no analysis proceeds without an explicit readiness check, no comparative claim is made without a real statistical test (elastic-net cross-validation with DeLong's test comparing fused against single-omics performance), and no functionality is described to the user as available unless the underlying code actually implements it. The scientific value of this design lies not in any single algorithm but in the disciplined combination of complementary molecular layers, honest baseline comparison, and stability-aware patient stratification — a methodological template for multi-omics analysis rather than a validated clinical tool.

---

## 15. Short Thesis Paragraph

The Multi-Omics module integrates matched gene-expression and DNA-methylation data from an RA anti-TNF cohort (or user-supplied data) through a validated, sample-matched, batch-corrected dataset pipeline feeding three integration methods — supervised DIABLO, unsupervised Similarity Network Fusion, and unsupervised MOFA2 factor analysis — each benchmarked against single-omics baselines with real cross-validated statistics rather than reported in isolation. Downstream sub-modules extract a stability-labeled joint biomarker signature, test whether each candidate gene's expression and methylation changes agree with its regulatory context (promoter versus gene body), and enrich the resulting list against standard pathway databases. A dedicated harmonization step and a session-summary tab enforce readiness checks and honest reporting throughout. The result is a scientifically disciplined demonstration of why combining molecular layers — rather than analyzing each in isolation — yields more robust, biologically interpretable candidate biomarkers and patient stratifications, without overstating clinical validation the implementation does not perform.
