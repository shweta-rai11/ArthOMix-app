# `mod_preprocessing.R` + `mod_preprocessing_explore.R` — Full Teaching, Audit, and Thesis-Documentation Notes

Files:
- `ArthOMix/R/transcriptomics/mod_preprocessing.R` (2,109 lines) — Preprocessing, Merge Datasets, Batch Correction tabs.
- `ArthOMix/R/transcriptomics/mod_preprocessing_explore.R` (1,060 lines) — the nested "Data Exploration" tab (its own Shiny module, `id = "eda"`).
- Supporting, pre-audited helpers this module calls into but does not define, all in `ArthOMix/global.R`: `load_default_dataset()`, `get_raw_eset()`, `get_collapsed_genes()`, `collapse_probes_to_genes()` (the `ExpressionSet` version), `eset_harmonize_meta()`, `load_individual_dataset()`, `merged_training_subset()`, `tx_parse_expr_matrix_rds()`, `compute_sample_qc()`, `summarize_norm_diagnostics()`, `needs_quantile_norm()`, `expr_raw_health()`, `detect_expr_data_type()`, `pca_of()`, `scree_plot()`, `plot_pca_advanced()`, `overlap_region_sizes()`, `draw_overlap_venn()`, `ARTHOMIX_COLORS`, `arthomix_pair()`, `theme_arthomix()`, `%||%`. And from `mod_dataset.R`: `default_dataset_entry`, `preloaded_choices()`.
Prepared: 2026-08-25.

**A note on method, stated up front rather than left implicit:** the brief asks for a separate `Code / What / Purpose / Inputs / Output / Data transformation / Why scientifically / Shiny behavior / Beginner / Advanced / Validation` block for *every single line* of ~3,170 lines of code. That is not a viable document — it would run to tens of thousands of entries, most of them repeating the same explanation for lines that are one logical statement split across several lines of R source (a very common style in this codebase; e.g. a single `tagList(...)` call spans 40 lines). This file instead teaches **every logically distinct statement, exactly once, in the order it executes**, grouping only lines that are one indivisible expression (an `if/else` chain, one `tagList()`, one pipe of `%>%` calls). Every code line in both files is covered by some block below — nothing is skipped — but adjacent lines that do one thing together are taught together, the same convention already used successfully in `mod_dataset_teaching_notes.md` and `mod_overview_teaching_notes.md`. Where the code repeats an already-taught pattern verbatim (e.g. the "empty-note / status" `renderUI` idiom appears well over a dozen times), later occurrences are taught by cross-reference, not re-explained from scratch.

---

## PART 1 — MODULE STRUCTURE

### 1.1 What this module is

`mod_preprocessing.R` is the **Transcriptomics → Preprocessing and Batch Correction** submodule (`mod_preprocessing_config$id = "preprocessing"`, section "2.2" in `submodules_registry.R`). Its own header comment (L1–20) is the most reliable one-paragraph description in the codebase, and is quoted here because everything below expands on it, not contradicts it:

> "Everything on this page runs live, on whatever data is loaded: each source dataset is cleaned individually first (log2 scale check, sample filters), then datasets are merged onto their shared genes/probes (Venn/overlap diagram), then the merged cohort is normalised and batch-corrected (ComBat, limma::removeBatchEffect, or ComBat-seq)."

This is the single most architecturally important fact about the module: **nothing here reads a precomputed result file**. Every number, plot, and table on this page is computed in the current R session, from whatever data the user picked, the moment a button is clicked. This is a deliberate contrast with modules like DGE or MR, which read audited precomputed pipeline output for the bundled cohort — Preprocessing exists specifically so a user can *rebuild* that pipeline (or run their own) interactively and see every intermediate step.

### 1.2 The four subtabs

`mod_preprocessing_ui()` (L621–659) builds one `tabsetPanel(id = ns("tabs"))` with four `tabPanel()`s, in this fixed order:

| # | `tabPanel` value | Title (icon + label) | Purpose | Backing `renderUI` output |
|---|---|---|---|---|
| 1 | `"Preprocessing"` | ⚙ broom / "Preprocessing" | Load and individually clean one or more raw source datasets (sample filters, feature filters, log2 decision) | `output$preprocessing_tab_ui` |
| 2 | `"Merge datasets"` | code-merge / "Merge Datasets" | Intersect features across every preprocessed source, visualize the overlap (Venn/UpSet), merge into one matrix | `output$merge_tab_ui` |
| 3 | `"Batch correction"` | wand-magic-sparkles / "Batch Correction" | Normalize the merged matrix and remove batch effects (ComBat / limma / ComBat-seq / SVA), with full QC before/after | `output$batch_tab_ui` |
| 4 | `"Explore"` | magnifying-glass-chart / "Data Exploration" | A **standalone**, independent EDA tool — its own upload, unrelated to the pipeline in tabs 1–3 | `mod_data_exploration_ui(ns("eda"))` — a fully separate nested Shiny module |

Tabs 1–3 are a strict linear pipeline that all read and write **one shared chain of reactive values inside `mod_preprocessing_server()`'s closure**: tab 1's per-source `result` reactives feed tab 2's `merge_inputs()`, tab 2's `merged()` feeds tab 3's `result` (batch correction) reactive, and tab 3's "Use this as the active dataset" button is the *only* place this module writes to the app-wide `dataset` reactiveValues. Tab 4 is architecturally severed from this chain on purpose (see §1.3).

### 1.3 How the subtabs are connected — the actual data-flow

```
                          ┌─────────────────────────────────────────────┐
                          │   dataset (app-wide reactiveValues,         │
                          │   defined in server.R, shared everywhere)   │
                          │   $expr / $meta / $source  (ACTIVE)         │
                          │   $staged_expr/_meta/_source (PREVIEW, from │
                          │   mod_dataset.R's Dataset tab)               │
                          └───────────────┬───────────────┬─────────────┘
                                          │ read (staged ⊕ active)      │ write (ONLY from
                                          │                             │ tab 3's activate_btn)
                                          ▼                             │
   TAB 1 — Preprocessing                                                │
   ┌──────────────────────────────────────────────────┐                │
   │ mod_pp_source_server() × up to MAX_PP_SOURCES(6)  │                │
   │   source: upload | preloaded GEO | "current"      │                │
   │   → raw_pair() → sample filters → feature filters │                │
   │   → log2 decision → result() [eventReactive: run] │                │
   │                                                    │                │
   │ pp_preloaded_read() (checkbox multi-select path)  │                │
   │   → preloaded_results_val [observeEvent:           │                │
   │     preloaded_run]                                │                │
   └───────────────────────┬────────────────────────────┘                │
                            │ merge_inputs() reads preloaded_results()   │
                            ▼                                            │
   TAB 2 — Merge datasets                                                │
   ┌──────────────────────────────────────────────────┐                 │
   │ selected_lst() → optional pp_collapse_probes_to_  │                 │
   │  genes() → overlap_sets()/draw_overlap_venn()     │                 │
   │  → merged() [eventReactive: merge_btn /            │                 │
   │    merge_use_example_btn]                          │                 │
   │  (alternate: example_live_merge(), the "Merge the  │                 │
   │  example pipeline's training datasets" radio path) │                 │
   └───────────────────────┬────────────────────────────┘                │
                            │ merged()                                    │
                            ▼                                             │
   TAB 3 — Batch correction                                               │
   ┌──────────────────────────────────────────────────┐                  │
   │ settings_ui() reads merged()$meta's columns        │                  │
   │  → filter/normalize (TMM or quantile) → outlier    │                  │
   │  exclusion → ComBat / limma / SVA / ComBat-seq      │                  │
   │  → result() [eventReactive: run_btn]                │                  │
   │  → activate_btn ─────────────────────────────────────────────────────┘
   └────────────────────────────────────────────────────┘

   TAB 4 — Data Exploration (SEVERED from the chain above)
   ┌──────────────────────────────────────────────────┐
   │ mod_data_exploration_server("eda") — its own       │
   │  fileInput, its own eda_result() eventReactive,    │
   │  never reads `dataset`, never reads tabs 1-3's      │
   │  state, never writes anything back anywhere.        │
   └──────────────────────────────────────────────────┘
```

Three points worth stating explicitly, because they are easy to get wrong just from using the app and are the kind of implementation detail a thesis Methods section should get right:

1. **"Currently loaded dataset" always means the Dataset tab's *staged* preview if one exists, falling back to the *active* dataset only if nothing has been staged.** This rule (`dataset$staged_expr %||% dataset$expr`) is duplicated in three independent places in this file (`pp_preloaded_read()` L145–147, `mod_pp_source_server`'s `current_source()` L303–306, and nowhere else needs it) rather than centralized in one function — a real, disclosable coupling risk (see §7's audit notes).
2. **Nothing advances automatically.** Every step of the pipeline requires its own explicit button click — `input$run` per source, `input$preloaded_run` for the bundled-cohort path, `input$merge_btn`/`input$merge_use_example_btn` to merge, `input$run_btn` to normalize/batch-correct, `input$activate_btn` to make the result the app-wide active dataset. `pp_progress()` (L673–681) exists purely to track this state for the right-column stepper UI; it has no effect on the pipeline itself.
3. **The Merge tab has two entirely different code paths** depending on `input$merge_mode`: `"own"` (build from whatever was preprocessed in tab 1) or `"example"` (rebuild the bundled RA training cohort live, from its two raw GEO sources, bypassing tab 1 entirely). These are documented separately in §3.

### 1.4 UI functions vs. server functions

| Function | Role |
|---|---|
| `mod_pp_field_hint(text)` | UI-only helper (a hover tooltip fragment) |
| `mod_pp_source_ui(id, default_gse, n_sources_id)` | UI builder for one source's upload/filter block |
| `mod_pp_source_server(id, default_label, default_gse, dataset)` | Server logic for one source block; returns a `result` reactive |
| `pp_tab_title(ic, label)` | UI-only helper (icon + label tab title) |
| `mod_preprocessing_ui(id)` | Top-level UI builder — the 4-tab `tabsetPanel` |
| `mod_preprocessing_server(id, dataset, results)` | Top-level server — owns tabs 1–3's logic and mounts tab 4 |
| `mod_data_exploration_ui(id)` | UI builder for the Data Exploration tab |
| `mod_data_exploration_server(id)` | Server logic for the Data Exploration tab |

Everything else defined at file scope in `mod_preprocessing.R` (`pp_guess_col`, `pp_collapse_probes_to_genes`, `pp_cohort_label`, `pp_cohort_choices`, `pp_preloaded_read`, and the constants `MAX_PP_SOURCES`/`PP_COHORT_LABELS`/`PP_MERGED_COHORT_LABEL`) is a **pure helper or data-processing function** — no `input`/`output`/`session`, callable and testable outside of any Shiny session. `mod_preprocessing_explore.R` is almost entirely pure helper/statistical/plotting functions (see §1.6) plus exactly two Shiny-aware functions (`mod_data_exploration_ui`, `mod_data_exploration_server`).

### 1.5 Function inventory — `mod_preprocessing.R`

**Top-level (file-scope) functions: 13.** Inside the two `moduleServer()` bodies (`mod_pp_source_server`, `mod_preprocessing_server`), there are approximately **70 additional named reactive constructs** — `reactive()`/`eventReactive()` expressions, `observe()`/`observeEvent()` blocks, and `output$... <- render*(...)` bindings — each of which is itself a function in R's implementation, but conventionally discussed as "a reactive," "an observer," or "an output," not as a free-standing function. The table below inventories every one of them.

| Function / reactive construct | Type | Purpose | Input | Output | Used by |
|---|---|---|---|---|---|
| `mod_pp_field_hint()` | UI helper | Renders a hover-tooltip `<span>` | `text` string | Shiny tag | `mod_pp_source_ui` |
| `pp_guess_col()` | Utility | Case-insensitive exact→substring column-name guesser | `cols, exact, contains, fallback` | one column name | `pp_collapse_probes_to_genes`, batch-upload UI |
| `pp_collapse_probes_to_genes()` | Data-processing | Collapse a probe-level matrix to one row per gene via an uploaded annotation table | `expr` matrix, `annot` data.frame, `method` | collapsed matrix | Merge tab's optional probe-collapsing step |
| `pp_cohort_label()` | Utility | GEO accession → display label | `id` string | string | `pp_cohort_choices`, several UI blocks |
| `pp_cohort_choices()` | Utility | Builds the `selectInput`-ready named vector of bundled cohorts | none | named character vector | `mod_pp_source_ui`, `preprocessing_tab_ui` |
| `pp_preloaded_read()` | Data-processing + validation | Reads + log2-normalizes one bundled/"current" dataset for the checkbox multi-select path | `choice_id, log2_choice, dataset` | `list(label, expr, meta, n_samples_before/after, n_genes_before/after, log2_applied)` | `observeEvent(input$preloaded_run)` |
| `mod_pp_source_ui()` | UI | Full per-source upload/filter panel | `id, default_gse, n_sources_id` | Shiny tags | `mod_preprocessing_ui` (indirectly, via `pp_sources`) |
| `mod_pp_source_server()` | Server (moduleServer) | Everything for one source: load → filter → preprocess | `id, default_label, default_gse, dataset` | a `result` reactive | `pp_sources` (instantiated 6×) |
| `pp_tab_title()` | UI helper | icon+label tab title | `ic, label` | Shiny tags | `mod_preprocessing_ui` |
| `mod_preprocessing_ui()` | UI | The whole 4-tab page | `id` | Shiny tags | called once from `ui_shell.R`/the transcriptomics tab registry |
| `mod_preprocessing_server()` | Server (moduleServer) | Everything for tabs 1–3, plus mounting tab 4 | `id, dataset, results` | (side effects only — no return value used) | called once from `server.R` |

**Reactive expressions / observers / outputs inside `mod_pp_source_server()` (≈20 per instance, instantiated 6×):** `output$source_type_ui`, `source_type`, `use_preloaded`, `use_upload`, `use_current`, `current_source`, `output$current_note`, `meta_raw`, `expr_raw_preview`, `output$upload_preview_ui`, `guess_col` (a locally-scoped duplicate of `pp_guess_col`), `output$colmap`, `raw_pair`, `output$group_filter_ui`, `output$filter_val_ui`, `output$numeric_filter_ui`, `output$num_filter_range_ui`, `output$dedup_col_ui`, `output$log2_ui`, `result` (`eventReactive(input$run, ...)`), `output$status_ui`.

**Reactive expressions / observers / outputs inside `mod_preprocessing_server()` (≈50):** `pp_progress`, `output$pipeline_summary`, `available_example_groups`, `pp_source_default_gse`, `pp_sources` (a `lapply` instantiating 6 `mod_pp_source_server` calls), `preloaded_results_val`/`observeEvent(input$preloaded_run)`/`preloaded_results`, `output$preloaded_status_ui`, `output$preprocessing_tab_ui`, `merge_inputs`, `output$merge_tab_ui`, `output$merge_example_ui`, `output$merge_venn_example_ui`, `output$merge_example_composition_table`, `example_overlap_sets`, `example_merge_from_raw`, `example_live_merge`, `venn_plot_example_obj`/`output$venn_plot_example`/`output$download_venn_example_png`, `output$venn_table_example`, `venn_regions_example`/`output$venn_region_table_example`/`output$download_venn_example`, `output$merge_select_ui`, `collapse_annot`, `selected_lst`, `output$merge_venn_ui`, `overlap_sets`, `venn_plot_custom_obj`/`output$venn_plot_custom`/`output$download_venn_custom_png`, `output$venn_table_custom`, `venn_regions_custom`/`output$venn_region_table_custom`/`output$download_venn_custom`, `merged` (`eventReactive`), `output$merge_summary_ui`, `output$download_merged_expr`/`_meta`/`_rds`, `output$merge_composition_table`, `active_meta_df`, `output$settings_ui`, `output$ref_batch_ui`, `observeEvent(list(input$batch_col, input$batch_col2))`, `result` (`eventReactive(input$run_btn, ...)` — the batch-correction pipeline itself), `output$vb_samples`/`vb_genes_kept`/`vb_genes_dropped`/`vb_flagged`, `output$decisions_ui`, `output$signal_plot`/`detected_plot`/`cor_plot`, `dist_summary`, `output$dist_plot`, `output$pca_before`/`pca_after`/`scree_plot`, `assoc_pvalue`, `output$summary_ui`, `output$norm_table`/`download_norm`, `qc_table_display`/`output$qc_table`/`download_qc`, `pca_table`/`output$pca_table`/`download_pca`, `output$activate_ui`/`observeEvent(input$activate_btn)`, `bc_section` (a plain non-reactive helper defined inline), `output$results_top_ui`, `output$results_rest_ui`, `batch_content` (a plain `tagList`, not a reactive), `output$batch_tab_ui`, the `outputOptions(..., suspendWhenHidden = FALSE)` loop, and the final `mod_data_exploration_server("eda")` call that mounts tab 4.

### 1.6 Function inventory — `mod_preprocessing_explore.R`

**39 top-level pure functions** (no Shiny reactivity — safe to call and unit-test outside any session), organized by role:

| Category | Functions |
|---|---|
| Upload parsing | `eda_parse_upload()` |
| Pure statistics | `eda_skewness()`, `eda_kurtosis()`, `eda_robust_z()`, `eda_skew_label()` |
| Dataset/feature/sample summaries | `eda_overview()`, `eda_descriptive_stats()`, `eda_normality_summary()`, `eda_normalization_assessment()` |
| Data prep for structure analysis | `eda_impute_median()`, `eda_prep_for_structure()`, `eda_pca()`, `eda_sample_correlation()`, `eda_hclust()` |
| Outlier detection | `eda_sample_outliers()`, `eda_feature_outliers()` |
| Missingness / variance / transform diagnostics | `eda_missingness()`, `eda_mean_variance_df()`, `eda_transform_diagnostic()` |
| Final synthesis | `eda_final_summary()` |
| Plot builders (pure: data in, ggplot/plotly object out) | `eda_value_axis_label()`, `eda_hist_plot()`, `eda_density_plot()`, `eda_box_plot()`, `eda_box_plot_interactive()`, `eda_violin_plot()`, `eda_sample_density_plot()`, `eda_qq_plot()`, `eda_meanvar_plot()`, `eda_missing_bar_plot()`, `eda_transform_diag_plot()`, `eda_corr_heatmap_plot()`, `eda_dendro_plot()`, `eda_pca_plot()`, `eda_scree_plot()` |
| UI-composition helpers | `eda_section_card()`, `eda_status_panel_ui()`, `eda_summary_card_ui()`, `eda_upload_info_ui()` |
| Shiny module functions | `mod_data_exploration_ui()`, `mod_data_exploration_server()` |

**Inside `mod_data_exploration_server()`: ≈35 reactive/observer/output constructs** — `raw_data` (`reactiveVal`), `observeEvent(input$raw_file)`, `output$body_ui`, `output$head_preview_table`, `eda_result` (`eventReactive(input$run_btn, ...)` — the single pipeline that computes everything), and one `output$..._ui`/`output$..._plot`/`output$..._table` pair per lettered section (A through M: overview, descriptive stats, distribution, normality, normalization status, outliers, PCA, correlation, missingness, mean-variance, transform diagnostic, summary, plus `download_summary` and the final `results_ui` assembler).

### 1.7 Grand total

- **Top-level named functions across both files: 13 + 39 + 2 (the module UI/server pair already counted in the 13) = 52 distinct top-level functions**, of which 8 are UI/Shiny-module functions and 44 are plain R functions (utility, data-processing, validation, statistical, or plotting).
- **Reactive expressions, observers, and render-outputs across both files: ≈70 (mod_preprocessing.R) + ≈35 (mod_preprocessing_explore.R) ≈ 105.**
- **Grand total: ≈157 named callable units** implementing this module. This is the number quoted for "how many functions are in this module" — the brief's category list (UI / server / helper / data-processing / validation / statistical / plotting / reactive / utility) is applied per-function throughout Parts 2–5 below rather than repeated as one flat list here, since several functions genuinely belong to more than one category (e.g. `pp_preloaded_read()` is simultaneously data-processing *and* validation).

---

## PART 2 — SUBTAB: "Preprocessing"

### 2.A Purpose

**Simple:** Before you can compare or combine two datasets, each one needs to be cleaned up on its own — remove samples that don't belong (wrong diagnosis, missing values in the field you care about, obvious duplicates), remove genes that are mostly missing, and make sure the numbers are on the right scale (log2 or not). This tab does that, one dataset at a time.

**Intermediate:** Every gene-expression dataset arrives with its own quirks: some samples might be technical replicates, some might have failed QC upstream and need excluding by a clinical field, some platforms report raw linear intensities while others report already-logged values. Doing this cleaning *before* merging is important because merging first and filtering second would let one bad dataset's problems (e.g. an unlogged linear-scale platform sitting next to a log2 one) contaminate the combined matrix in a way that's much harder to diagnose after the fact. This tab runs entirely on **raw, single-platform, not-yet-merged** data — its own UI note says so explicitly (L231–232) — so every filter here acts on one homogeneous dataset.

**Advanced:** This is the "Filter → (per-source log2 Transform)" stage of the pipeline global.R's own audit comment (L1497–1507) describes as enforced by data flow, not convention: sample/feature filters and the log2 decision happen here, strictly before the Merge tab's feature intersection and before Batch Correction's normalization. This ordering matters statistically — deciding "does this look log2 already?" (`q99 > 100` heuristic, §2.E) on a per-source basis, before any cross-platform merge, means the decision is made on a genuinely single-scale, single-platform value distribution rather than on a matrix that already conflates two platforms' dynamic ranges. Skipping or reordering this step — e.g. merging two datasets first and only then checking "does this look logged" on the combined matrix — risks a decision that's correct for neither platform alone (e.g. a q99 that averages a log2 platform's ~15 and a linear platform's ~4,000 into something that crosses the threshold in either direction depending on relative sample counts).

**What happens if this step is skipped or done incorrectly:** If a linear-scale microarray dataset is never log2'd, its expression values (often in the thousands) dominate any downstream variance-based step (PCA, quantile normalization, ComBat's empirical Bayes shrinkage) purely because of scale, not real biological signal — every batch-correction and clustering result downstream becomes an artifact of that one dataset's units, not biology. If a sample filter is set to the wrong diagnosis group, an entirely different comparison than intended gets merged and batch-corrected without any downstream step noticing (there is no re-validation later that the merged cohort's group composition matches what was intended).

### 2.B UI inventory

Tab 1's UI is built from **two independent boxes**, both produced by `output$preprocessing_tab_ui` (L791–809) and `mod_pp_source_ui()` (L205–269, instantiated once per "Dataset N" block, though the checkbox path below is what's actually shown on the tab as currently wired — see §2.D for why `mod_pp_source_ui`/`mod_pp_source_server` exist as a second, parallel per-source path used by the Merge tab's "own data" option, not directly rendered on this tab's visible box).

**Box 1 — "Preloaded Data" (`output$preprocessing_tab_ui`, what the user actually sees on this tab):**

| Element | Input ID (post-`ns()`) | What the user sees | What it controls | Server code that reacts |
|---|---|---|---|---|
| `checkboxGroupInput(ns("preloaded_selected"), ...)` | `...-preloaded_selected` | A row of checkboxes: the 4 bundled cohorts (professional labels via `pp_cohort_choices()`) + "Currently Loaded Dataset" | Which cohorts get loaded+preprocessed on the next click | `observeEvent(input$preloaded_run)` |
| `radioButtons(ns("preloaded_log2"), ...)` | `...-preloaded_log2` | 3 choices: Auto-detect / Force log2 / Skip | The log2 decision applied identically to every checked cohort | `pp_preloaded_read()`'s `log2_choice` argument |
| `actionButton(ns("preloaded_run"), ...)` | `...-preloaded_run` | "Load and Preprocess Selected Cohorts" | Triggers the whole read+filter+log2 pipeline for every checked cohort | `observeEvent(input$preloaded_run)` |
| `uiOutput(ns("preloaded_status_ui"))` | — (output) | Per-cohort green "kept N/M samples..." or red error line | Read-only | `output$preloaded_status_ui` |

**Box 2 — per-source panel, `mod_pp_source_ui()` (used by the Merge tab's "own data" path, and structurally identical to what an individual `pp_sources[[i]]` instance renders):**

| Element | Input ID | What the user sees | Controls | Reacts via |
|---|---|---|---|---|
| `radioButtons(ns("source_type"), ...)` (built in `output$source_type_ui`) | `...-source_type` | 3 choices: "Upload files" / "A bundled cohort" / "Currently loaded dataset" | Which of 3 code branches `raw_pair()` takes | `source_type`, `use_preloaded`/`use_upload`/`use_current` |
| `uiOutput(ns("current_note"))` | — | "Using X: N genes × M samples" or a warning | Read-only | `output$current_note` |
| `selectInput(ns("preloaded_choice"), ...)` | `...-preloaded_choice` | Cohort dropdown (hidden via `shinyjs::hidden()` when `default_gse` is fixed, e.g. Dataset 1/2's training-cohort slots) | Which bundled GSE to load | `raw_pair()` |
| `textInput(ns("label"), ...)` | `...-label` | Free-text dataset label | Display name used everywhere downstream | `raw_pair()`'s `label` field |
| `fileInput(ns("expr_file"), ...)` | `...-expr_file` | Expression matrix upload (.csv/.rds) | Raw matrix source | `expr_raw_preview()` |
| `fileInput(ns("meta_file"), ...)` | `...-meta_file` | Metadata upload (.csv/.rds) | Raw metadata source | `meta_raw()` |
| `uiOutput(ns("upload_preview_ui"))` | — | "Read N features × M samples..." | Read-only | `output$upload_preview_ui` |
| `uiOutput(ns("colmap"))` (4 nested `selectInput`s: `map_id`/`map_group`/`map_sex`/`map_batch`) | — | Column-mapping dropdowns, auto-guessed | Which metadata columns become `sample`/`group`/`sex`/`batch` | `raw_pair()` |
| `checkboxInput(ns("dedup"), ...)` + `uiOutput(ns("dedup_col_ui"))` | `...-dedup`, `...-dedup_col` | "Deduplicate samples by an ID column" toggle + column picker | Whether/how `result()` drops duplicate rows | `result()` |
| `uiOutput(ns("group_filter_ui"))` (→ `filter_col`, `filter_val_ui` → `filter_vals`) | `...-filter_col`, `...-filter_vals` | "Keep only samples where [column] equals one of [values]" | Categorical sample filter | `result()` |
| `uiOutput(ns("numeric_filter_ui"))` (→ `num_filter_col`, `num_filter_range_ui` → `num_filter_range`) | `...-num_filter_col`, `...-num_filter_range` | "...and within this numeric range" slider | Numeric sample filter (e.g. age, RIN) | `result()` |
| `uiOutput(ns("log2_ui"))` | `...-log2` | 3-way log2 radio, default differs by source type (see §2.E) | Log2 decision override | `result()` |
| `sliderInput(ns("max_na_pct"), ...)` | `...-max_na_pct` | 0–80%, step 5, default 0 | Feature missing-data tolerance | `result()` |
| `textInput(ns("exclude_pattern"), ...)` | `...-exclude_pattern` | Free-text regex | Feature-ID exclusion pattern (e.g. `^AFFX`) | `result()` |
| `actionButton(ns("run"), ...)` | `...-run` | "Preprocess this dataset" | Triggers `result` (`eventReactive`) | `result` |
| `uiOutput(ns("status_ui"))` | — | Green/red summary line | Read-only | `output$status_ui` |

### 2.C Function inventory — Preprocessing tab

| Function | Type | Purpose | Input | Output | Used by |
|---|---|---|---|---|---|
| `mod_pp_field_hint()` | UI helper | Hover tooltip markup | `text` | Shiny tag | `mod_pp_source_ui` |
| `pp_guess_col()` | Utility | Column-name guesser (exact→substring→fallback) | `cols, exact, contains, fallback` | 1 column name | `pp_collapse_probes_to_genes`; also duplicated locally as `guess_col()` inside `mod_pp_source_server` |
| `PP_COHORT_LABELS`, `PP_MERGED_COHORT_LABEL` | Constants | GEO accession → professional display label map | — | named vector / string | `pp_cohort_label()` |
| `pp_cohort_label()` | Utility | id → display label, with fallback | `id` | string | UI labels throughout tabs 1–2 |
| `pp_cohort_choices()` | Utility | Builds the bundled-cohort dropdown/checklist choices | none | named vector | `preprocessing_tab_ui`, `mod_pp_source_ui` |
| `pp_preloaded_read()` | Data-processing + validation | Read one bundled/"current" cohort, harmonize `batch` column, auto-log2, median-impute | `choice_id, log2_choice, dataset` | `list(label, expr, meta, n_*_before/after, log2_applied)` | `observeEvent(input$preloaded_run)` |
| `mod_pp_source_ui()` | UI | One source's full panel | `id, default_gse, n_sources_id` | Shiny tags | `pp_sources` instantiation (indirectly; see note below) |
| `mod_pp_source_server()` | Server | One source's load→filter→preprocess pipeline | `id, default_label, default_gse, dataset` | `result` reactive | `pp_sources <- lapply(1:6, ...)` |
| `output$source_type_ui` | Reactive output | Radio choices for data source, remembers current selection across re-renders | `input$source_type` (isolated) | `radioButtons` | UI |
| `source_type`/`use_preloaded`/`use_upload`/`use_current` | Reactive expressions | Thin wrappers turning the radio value into 3 booleans | `input$source_type` | reactive value | gate every downstream block |
| `current_source` | Reactive expression | Reads staged-or-active `dataset` | `dataset` (shared) | `list(expr, meta, label)` | `current_note`, `raw_pair` |
| `output$current_note` | Reactive output | "Using X: N × M" note | `current_source()` | Shiny tag | UI |
| `meta_raw` | Reactive expression | Parses uploaded metadata file (CSV/RDS) | `input$meta_file` | data.frame | `upload_preview_ui`, `colmap`, `raw_pair` |
| `expr_raw_preview` | Reactive expression | Parses uploaded expression file, preview only | `input$expr_file` | matrix | `upload_preview_ui` |
| `output$upload_preview_ui` | Reactive output | "Read N features × M samples..." or parse error | `expr_raw_preview()`, `meta_raw()` | Shiny tag | UI |
| `guess_col()` (local) | Utility | Identical logic to `pp_guess_col()`, scoped inside the server closure | same signature | 1 column name | `output$colmap` |
| `output$colmap` | Reactive output | 4 auto-guessed mapping dropdowns | `colnames(meta_raw())` | Shiny tags | UI, feeding `input$map_*` |
| `raw_pair` | Reactive expression | The single point where all 3 source types converge into one `(expr, meta, label)` shape, pre-filter | `source_type()`, and whichever inputs that branch needs | `list(expr, meta, label)` | every UI-population `output$*_filter_ui`, `result()` |
| `output$group_filter_ui` | Reactive output | Categorical filter column picker | `raw_pair()$meta` columns | Shiny tags | UI, feeds `input$filter_col` |
| `output$filter_val_ui` | Reactive output | Multi-select of the chosen column's distinct values, all pre-checked | `input$filter_col`, `raw_pair()` | `checkboxGroupInput` | UI, feeds `input$filter_vals` |
| `output$numeric_filter_ui` | Reactive output | Numeric-column picker (only shown if ≥1 numeric metadata column exists) | `raw_pair()$meta` | Shiny tags | UI |
| `output$num_filter_range_ui` | Reactive output | Range slider bounded by the chosen column's actual min/max | `input$num_filter_col`, `raw_pair()` | `sliderInput` | UI, feeds `input$num_filter_range` |
| `output$dedup_col_ui` | Reactive output | Dedup-by column picker | `raw_pair()$meta` | Shiny tags | UI |
| `output$log2_ui` | Reactive output | 3-way log2 radio; **default differs by path** (`"skip"` for upload, `"auto"` otherwise) | `use_upload()` | `radioButtons` | UI, feeds `input$log2` |
| `result` (`eventReactive(input$run, ...)`) | Data-processing + statistical | The actual filter/transform pipeline for one source (§2.D, §2.E) | `raw_pair()`, every filter input | `list(label, expr, meta, n_*_before/after, log2_applied)` | `output$status_ui`; externally, `pp_sources[[i]]` |
| `output$status_ui` | Reactive output | "N of M samples kept, N of M features kept..." or an error | `result()` | Shiny tag | UI |
| `pp_source_default_gse()` | Utility | Maps source index 1/2 → the two training GSEs, else `NULL` | `i` | string or `NULL` | `pp_sources` construction |
| `pp_sources` | `lapply` of 6 `mod_pp_source_server()` calls | Instantiates all 6 per-source server instances up front | `dataset` | list of 6 `result` reactives | never directly read elsewhere in this file (kept alive for state persistence; see note below) |
| `preloaded_results_val`/`observeEvent(input$preloaded_run)` | ReactiveVal + observer | Runs `pp_preloaded_read()` over every checked cohort, catching per-cohort errors independently | `input$preloaded_selected`, `input$preloaded_log2` | side effect: sets `preloaded_results_val` | `preloaded_results()`, `merge_inputs()` |
| `preloaded_results` | Plain accessor function | Thin wrapper around the reactiveVal | none | list of per-cohort results | `output$preloaded_status_ui`, `merge_inputs` |
| `output$preloaded_status_ui` | Reactive output | Per-cohort success/error line | `preloaded_results()` | Shiny tags | UI |
| `output$preprocessing_tab_ui` | Reactive output | Assembles the visible "Preloaded Data" box | (static markup + `pp_cohort_choices()`) | Shiny tags | tab 1's `tabPanel` |

**Note on `pp_sources`/`mod_pp_source_ui`/`mod_pp_source_server` vs. the visible tab:** reading `output$preprocessing_tab_ui` (L791–809) shows that what's actually rendered on the "Preprocessing" tab today is **only** the "Preloaded Data" checkbox box driven by `pp_preloaded_read()` — `mod_pp_source_ui()`'s richer per-source panel (with its own upload/filters/log2/run button) is never inserted into `mod_preprocessing_ui()`'s tab 1 markup at all. `pp_sources` is still instantiated (`lapply(seq_len(MAX_PP_SOURCES), ...)`, L735–738) and its 6 `mod_pp_source_server()` calls do register the reactive machinery described above (including `pp_progress()` reading `input$preloaded_selected`, not `pp_sources`), but nothing in `mod_preprocessing_ui()` ever calls `mod_pp_source_ui()`. This is exactly the kind of "server logic with no corresponding UI control" the audit brief (§12) asks to be flagged explicitly — see §7's **DESIGN ISSUE** entry for the full discussion; the practical consequence is that a user cannot currently reach the "own data" per-source upload+filter workflow this rich panel implements from tab 1 at all, only indirectly through the Merge tab, which is discussed next.

### 2.D Line-by-line teaching

#### Block 1 — `mod_pp_field_hint()` (L41–45)

```r
mod_pp_field_hint <- function(text) {
  tags$span(class = "field-hint", tabindex = "0",
            icon("circle-info"),
            tags$span(class = "field-hint-box", text))
}
```
**What/Purpose:** A pure UI-generating function — no reactivity. Returns a `<span>` wrapping a small info icon plus a second, CSS-hidden-until-hover `<span>` carrying the actual tooltip text. `tabindex = "0"` makes the icon keyboard-focusable, so the tooltip is also reachable without a mouse (an accessibility detail, not a Shiny-specific one). **Data flow:** `text` string in → nested Shiny tags out. **Shiny behavior:** none — this runs once, at UI-definition time, producing static HTML; there is no server-side counterpart. **Beginner:** this is exactly like writing a small reusable HTML snippet as an R function, the same way you'd write a Python function that returns an HTML string. **Advanced:** the actual hover behavior is pure CSS (`.field-hint*` rules in `www/custom.css`, referenced in the file's own comment) — deliberately avoiding a JS tooltip library dependency for something this simple. **Validation:** visually inspect the rendered page; there's no reactive value to unit-test.

#### Block 2 — `pp_guess_col()` (L50–56)

```r
pp_guess_col <- function(cols, exact, contains = exact, fallback = cols[1]) {
  hit <- cols[tolower(cols) %in% tolower(exact)]
  if (length(hit) > 0) return(hit[1])
  hit <- cols[grepl(paste(contains, collapse = "|"), cols, ignore.case = TRUE)]
  if (length(hit) > 0) return(hit[1])
  fallback
}
```
This is **line-for-line identical** to `mod_dataset.R`'s own `guess_col()` (already fully taught in `mod_dataset_teaching_notes.md` §4 Block I) and to the locally-scoped copy inside `mod_pp_source_server` (below). All three implement the same two-tier strategy: case-insensitive **exact** match first (`tolower(cols) %in% tolower(exact)`), then case-insensitive **substring** match (`grepl(paste(contains, collapse="|"), ...)` — building one alternation regex from the candidate list), then a caller-supplied `fallback`. **Why duplicated three times instead of shared:** `mod_dataset.R`'s copy is private to that file's `moduleServer()` closure; this file's module-scope copy (`pp_guess_col`) exists so `pp_collapse_probes_to_genes()` (file-scope, needs it before any `moduleServer()` runs) can use it; and `mod_pp_source_server()`'s own local `guess_col()` (L363–369) is a third, byte-identical copy — the file's own comment at L47–49 explains this was a deliberate pull-out ("Same name-based guess mod_pp_source_server's own colmap uses internally... pulled out to module scope so the batch-upload path below can reuse it"), but the original in-module copy was left in place rather than replaced with a call to the new module-scope one. **DESIGN ISSUE** (flagged per the audit brief's requirement, not silently fixed): three independent copies of the same ~6-line function is a real drift risk — a future bug fix to one copy (e.g. handling of `NA`-valued column names) would need to be applied three times, and nothing enforces that.

#### Block 3 — `pp_collapse_probes_to_genes()` (L58–90)

```r
pp_collapse_probes_to_genes <- function(expr, annot, method = c("median", "maxmean", "mean")) {
  method <- match.arg(method)
  cols <- colnames(annot)
  probe_col <- pp_guess_col(cols, c("probe", "probe_id", "probeid", "probeset", "probesetid", "id"))
  gene_col  <- pp_guess_col(cols, c("gene_symbol", "genesymbol", "symbol", "gene"),
                             fallback = if (length(cols) >= 2) cols[2] else cols[1])
```
**What:** `match.arg(method)` is base R's idiom for validating that a string argument is one of a fixed set of choices, defaulting to the first (`"median"`) if the caller passes nothing — it will `stop()` with a clear error if given anything else, rather than silently misbehaving. `pp_guess_col()` (Block 2) is then used twice: once to find which column of the user-uploaded annotation file holds probe IDs, once for gene symbols (falling back to the annotation file's 2nd column, on the assumption a 2-column file is "probe, gene" in that order if name-guessing fails outright). **Why this matters scientifically:** this is the mapping table any probe-to-gene collapse depends on entirely — a wrong `gene_col` guess (e.g. picking a "platform" column instead of "gene_symbol") would silently produce a nonsense collapsed matrix with technically-valid-looking but biologically meaningless row names, and nothing downstream would catch it. **Validation:** upload a small (10-row) annotation file with known column names, run this function directly (`pp_collapse_probes_to_genes(test_expr, test_annot)`), and check `probe_col`/`gene_col` by inspecting which columns actually get used (add a `browser()` or check the returned matrix's structure against hand-computed expectations).

```r
  map <- stats::setNames(as.character(annot[[gene_col]]), as.character(annot[[probe_col]]))
  sym <- unname(map[rownames(expr)])
  keep <- !is.na(sym) & sym != "" & !grepl("///", sym, fixed = TRUE)
  validate(need(any(keep), "None of this expression matrix's row IDs matched the annotation file's probe-ID column, or every match was to more than one gene. Check the annotation file's columns."))
  ex <- expr[keep, , drop = FALSE]; sym <- sym[keep]
```
**What:** `setNames(values, names)` builds a **named character vector acting as a hash map**: `map["1007_s_at"]` → the gene symbol for that probe. `map[rownames(expr)]` then does a **vectorized lookup** — for every row of the expression matrix, look up its gene symbol by probe ID; any probe not present in the annotation file returns `NA` automatically (R's named-vector indexing behavior, not a special case coded here). `unname()` strips the resulting vector's names (which would otherwise just be the same probe IDs again, since that's what was used to index). **The `keep` filter, in words:** drop any probe with (a) no matching annotation row (`is.na(sym)`), (b) an empty-string symbol, or (c) a symbol containing `"///"` — the `///`-delimiter convention some Affymetrix annotation files use to record "this probe maps ambiguously to multiple gene symbols" (e.g. `"MIR4640///DDR1"`). **Why never split or duplicated:** the file's own header comment (L67–68) states this explicitly — an ambiguous probe is dropped entirely rather than assigned to either gene, which is the statistically conservative choice (assigning it to one gene arbitrarily would inject noise into that gene's aggregate signal; duplicating it into both genes would inject spurious correlation between two genes that shouldn't be correlated by measurement). **`validate(need(any(keep), ...))`:** the Shiny-only "stop this reactive and show this exact message" idiom (fully taught in `mod_dataset_teaching_notes.md` Block C) — guards against a total column-mapping failure (annotation file used the wrong ID system entirely) producing a silently-empty matrix instead of a clear, actionable error. **Data transformation:** `expr` (P probes × N samples) → `ex` (P' ≤ P probes × N samples, unchanged values, just row-subsetted) plus a parallel `sym` vector of gene symbols, one per remaining row.

```r
  if (method == "maxmean") {
    row_mean <- rowMeans(ex, na.rm = TRUE)
    keep_idx <- tapply(seq_along(sym), sym, function(idx) idx[which.max(row_mean[idx])])
    out <- ex[unlist(keep_idx), , drop = FALSE]
    rownames(out) <- names(keep_idx)
    return(out[order(rownames(out)), , drop = FALSE])
  }
  agg_fn <- if (method == "median") stats::median else base::mean
  apply(ex, 2, function(col_vals) tapply(col_vals, sym, agg_fn, na.rm = TRUE))
```
**What the three collapsing methods do, precisely:**
- **`"maxmean"`** (the branch above): for each gene symbol, among every probe mapped to it, keep only the **one probe with the highest across-sample mean expression** as that gene's representative row — every other probe mapping to the same gene is discarded outright. `tapply(seq_along(sym), sym, function(idx) idx[which.max(row_mean[idx])])` is the mechanism: group row-indices by gene symbol, and within each group pick the index of the row with the largest `row_mean`. This is **ArthOMix's own convention** for its bundled preloaded datasets (same logic as `global.R`'s `collapse_probes_to_genes()`, used on `ExpressionSet` objects — see §2.G) — it is a real, if debatable, choice: it discards information from every non-maximal probe rather than combining it.
- **`"median"`** (the default; falls into the `agg_fn <- stats::median` branch): for each gene, at each sample, take the **median value across every probe mapping to that gene**. `tapply(col_vals, sym, agg_fn, na.rm = TRUE)` does this per column (via the outer `apply(ex, 2, ...)`) — group the one sample's values by gene symbol, aggregate each group. The file's own comment attributes this to "most published methods, e.g. Zhu et al. 2021."
- **`"mean"`**: identical mechanism, `base::mean` instead of `stats::median`.

**Why the method choice is scientifically consequential, not cosmetic:** median/mean genuinely combine every probe's signal (more robust to a single noisy/mismeasured probe, but can dilute a strong signal that only one specific probe captures — e.g. if probes target different transcript isoforms with different expression levels); maxmean instead trusts one "best" probe entirely and discards the rest (more sensitive to that one probe's idiosyncratic noise, but preserves the sharpest available signal rather than averaging it away). Which is "correct" depends on why the probes differ in the first place (isoform-specific measurement vs. redundant technical replication of the same target) — information this function has no way to know and does not attempt to infer. **Beginner:** `tapply(values, groups, function)` is base R's "group by, then apply" primitive — the same idea as SQL's `GROUP BY` + aggregate, or pandas' `.groupby().agg()`. **Validation:** on a small synthetic case (3 probes, 2 mapping to gene A with known values, 1 to gene B), compute all three methods by hand and compare to the function's output.

#### Block 4 — Bundled-cohort labels (L109–128)

```r
PP_COHORT_LABELS <- c(
  "GSE93272"  = "Whole Blood Training Cohort A",
  "GSE110169" = "Whole Blood Training Cohort B",
  "GSE15573"  = "PBMC Validation Cohort",
  "GSE89408"  = "Synovial Tissue Validation Cohort"
)
PP_MERGED_COHORT_LABEL <- "Whole Blood Training Cohort (Merged)"

pp_cohort_label <- function(id) {
  if (identical(id, default_dataset_entry$id)) return(PP_MERGED_COHORT_LABEL)
  if (id %in% names(PP_COHORT_LABELS)) unname(PP_COHORT_LABELS[[id]]) else id
}

pp_cohort_choices <- function() {
  ids <- unname(preloaded_choices())
  stats::setNames(ids, vapply(ids, pp_cohort_label, character(1)))
}
```
**What/Why:** `PP_COHORT_LABELS` is a second, independent labeling table from `mod_dataset.R`'s own `INDIVIDUAL_DATASET_LABELS` (already taught in `mod_dataset_teaching_notes.md` Block C) — same GEO accessions, deliberately *different* label text ("Whole Blood Training Cohort A" here vs. "Whole Blood Training Cohort A" — actually identical wording in this case, but the file's own comment at L101–108 explains the intent: this tab is "a working cohort picker for someone running the pipeline, so it reads as tissue + role rather than an accession lookup," while `mod_dataset.R`/Overview/Datasets deliberately keep the raw GEO accession visible for traceability). `pp_cohort_label()` checks the merged-cohort sentinel ID first (`identical(id, default_dataset_entry$id)` — the same `"__default_merged__"` sentinel from `mod_dataset.R`, imported by reference since both files are sourced into the same global environment), then falls back to the `PP_COHORT_LABELS` lookup, then to the raw `id` itself if neither matches — the comment at L106–108 explains this fallback is deliberately a `%in%` check rather than `%||%`, because a plain-vector `[[` on a missing name **errors** in R rather than returning `NULL` the way a list's `[[` would, so `%||%` (which only catches `NULL`) wouldn't actually protect against an unmapped GSE ID here. `pp_cohort_choices()` reuses `mod_dataset.R`'s `preloaded_choices()` (its 5-entry id↔label vector) but discards its labels entirely (`unname(preloaded_choices())` keeps only the ids) and rebuilds the *names* via `pp_cohort_label()` — i.e., same set of loadable datasets, this tab's own label text. **Why this matters for reproducibility/thesis writing:** these two independent labeling tables (`INDIVIDUAL_DATASET_LABELS` in `mod_dataset.R` and `PP_COHORT_LABELS` here) must be kept in sync by hand if a new bundled GEO source is ever added — the file's own comment explicitly flags this ("Add an entry here if a new bundled source is ever added to GEO_SOURCES... anything missing just falls back to its raw ID"), which is graceful degradation, not silent breakage, but is still a coupling worth disclosing (§7).

#### Block 5 — `pp_preloaded_read()` (L139–203)

This is the single most scientifically important function on this tab — it is the entire "read a raw cohort, log2 it, impute it" pipeline for the checkbox-driven UI path that's actually shown on screen.

```r
pp_preloaded_read <- function(choice_id, log2_choice, dataset = NULL) {
  if (identical(choice_id, "__current__")) {
    use_expr <- dataset$staged_expr %||% dataset$expr
    use_meta <- dataset$staged_meta %||% dataset$meta
    use_label <- dataset$staged_source %||% dataset$source
    validate(need(!is.null(dataset) && !is.null(use_expr),
                  "No dataset is currently loaded. Preview one on the Dataset tab first."))
    expr <- use_expr
    meta <- use_meta
    label <- use_label %||% "Currently Loaded Dataset"
```
**What:** The "Currently Loaded Dataset" branch. `%||%` (R's null-coalescing operator, defined once in `global.R`, used pervasively) returns its left side unless that's `NULL`, in which case it returns the right side — so `dataset$staged_expr %||% dataset$expr` reads as "prefer whatever's staged on the Dataset tab; otherwise fall back to whatever's already active app-wide." This exact 3-line pattern is duplicated (not shared via a helper) in `mod_pp_source_server`'s `current_source()` — see the DESIGN ISSUE noted in §1.3. **Validation guard:** `validate(need(!is.null(dataset) && !is.null(use_expr), ...))` — a plain-language message rather than Shiny's default error screen, guiding the user back to the Dataset tab.

```r
  } else {
    gse <- choice_id
    if (identical(gse, default_dataset_entry$id)) {
      d <- load_default_dataset()
      expr <- d$expr; meta <- d$meta
    } else if (identical(gse, "GSE89408")) {
      d <- load_individual_dataset(gse)
      validate(need(!is.null(d), paste("Raw data for", gse, "was not found on disk.")))
      expr <- d$expr; meta <- d$meta
    } else {
      eset <- get_raw_eset(gse)
      validate(need(!is.null(eset), paste("Raw file for", gse, "not found on disk.")))
      expr <- get_collapsed_genes(gse)
      meta <- eset_harmonize_meta(eset, gse)
      keep <- !is.na(meta$group)
      meta <- meta[keep, , drop = FALSE]
      expr <- expr[, meta$sample, drop = FALSE]
    }
    label <- pp_cohort_label(gse)
  }
```
**What — three genuinely different read paths, by platform type:**
1. **Merged/default cohort** (`gse == default_dataset_entry$id`, the `"__default_merged__"` sentinel): reads the bundled, already-merged-and-batch-corrected cohort directly via `load_default_dataset()` (`global.R` L180–187) — this is *not* a raw single-platform dataset at all, despite living in a "raw" cohort picker; picking "Merged Data" here effectively short-circuits the whole point of this tab (see §7's SCIENTIFIC RISK note).
2. **GSE89408** (RNA-seq, synovial tissue): `load_individual_dataset()` (`global.R` L772+) reads raw counts directly — already gene-level (no probe-to-gene collapse needed, since RNA-seq quantifies transcripts, not array probes), so `expr`/`meta` come back ready to use.
3. **Every other bundled GSE** (microarray): `get_raw_eset(gse)` reads the cached raw `ExpressionSet`; `get_collapsed_genes(gse)` (a two-tier in-memory + on-disk cache around `global.R`'s `collapse_probes_to_genes()`, the `ExpressionSet`-flavored twin of this file's own `pp_collapse_probes_to_genes()`, using the fixed `"MaxMean"` rule via `WGCNA::collapseRows()`) does the probe→gene collapse; `eset_harmonize_meta()` (`global.R` L731–743) derives `sample`/`dataset`/`group`/`sex` from the raw `ExpressionSet`'s phenotype columns via regex matching on column names containing `disease|status` and `gender|sex`; then `keep <- !is.na(meta$group)` drops any sample whose diagnosis couldn't be classified into HC/RA/SLE/"other" at all (**this is itself an unlabeled sample filter** — see §2.E filter #5).

**Why the branching matters scientifically:** the collapse method for the two training microarray cohorts is **hardcoded to `"MaxMean"`** here (via `get_collapsed_genes()`/`collapse_probes_to_genes()` in `global.R`), whereas the Merge tab's *optional* probe-collapsing step (for user-uploaded probe-level data) lets the user pick median/maxmean/mean — i.e., the bundled cohorts' probe-to-gene decision is fixed and cannot be changed from this UI, while a user's own uploaded probe-level data gets a real choice. This is a legitimate, disclosable asymmetry between "reproduce the example pipeline exactly as built" and "run your own analysis your way."

```r
  if (!"batch" %in% colnames(meta)) meta$batch <- NA_character_

  n_samples_before <- ncol(expr); n_genes_before <- nrow(expr)
  q99 <- suppressWarnings(stats::quantile(as.numeric(expr[expr > 0]), 0.99, na.rm = TRUE))
  needs_log <- if (identical(log2_choice, "force")) TRUE
               else if (identical(log2_choice, "skip")) FALSE
               else isTRUE(!is.na(q99) && q99 > 100)
  if (needs_log) {
    expr[expr <= 0] <- NA
    expr <- log2(expr)
    expr <- expr[stats::complete.cases(expr), , drop = FALSE]
  }
```
**What:** First, guarantee a `batch` column exists (as `NA` if nothing else set it) — every downstream consumer of this dataset's `meta` can now safely reference `meta$batch` without a `NULL`-column error, even if batch correction is never actually applied to it. Then the **auto-detect log2 heuristic**, this codebase's single recurring pattern (identical logic also appears in `mod_pp_source_server`'s `result()` and in `example_live_merge()`): compute the **99th percentile of strictly-positive values** (`expr[expr > 0]` — zeros and negatives excluded from the quantile computation, since a log2-scale dataset can legitimately contain negative values, e.g. after mean-centering, but the "is this linear-scale" signal specifically comes from *how large the positive tail is*), and if that q99 exceeds 100, conclude "this looks like linear-scale data, needs log2." `identical(log2_choice, "force")`/`"skip"` let the user override this heuristic explicitly; `"auto"` (or anything else) falls through to the heuristic. **Why q99, not max:** a single extreme outlier value (a saturated probe, a sequencing artifact) would make `max()` an unreliable signal; the 99th percentile is far more robust to that one bad data point while still capturing "is the bulk of the data's upper range linear-scale-sized." **Why 100 specifically:** log2-scale microarray/RNA-seq expression values conventionally fall in roughly the 0–20 range (2^20 ≈ 1,000,000, already an enormous raw expression value), so a q99 above 100 is a strong signal the data hasn't been logged yet — this is a heuristic threshold, not a statistically derived one, and the file's own later comment (L504–510, discussing the upload path's *different* default) explicitly documents its known failure mode: it "only looks at value magnitude... can't tell a large-valued RAW RNA-seq count matrix (correct DESeq2 input, should stay unlogged) apart from large-valued already-normalised data that genuinely needs logging."

**The transform itself, once triggered:** `expr[expr <= 0] <- NA` — log2 is undefined for zero and negative values, so these are converted to missing rather than producing `-Inf`/`NaN`; `expr <- log2(expr)`; `expr <- expr[stats::complete.cases(expr), , drop = FALSE]` — **any gene with even one remaining `NA` (from a zero/negative value in even one sample) is dropped entirely**, a hard `complete.cases()` filter, not an imputation, at this specific point in the pipeline.

```r
  if (anyNA(expr)) {
    row_med <- apply(expr, 1, stats::median, na.rm = TRUE)
    na_idx <- which(is.na(expr), arr.ind = TRUE)
    expr[na_idx] <- row_med[na_idx[, 1]]
  }

  list(label = label, expr = as.matrix(expr), meta = meta,
       n_samples_before = n_samples_before, n_samples_after = ncol(expr),
       n_genes_before = n_genes_before, n_genes_after = nrow(expr),
       log2_applied = needs_log)
}
```
**What:** A second, distinct missing-data handling step — this one **per-gene median imputation** for any residual `NA`s (the comment at L188–192 explains this specifically covers the case where the log2 branch was *skipped* because the data was auto-detected as already-logged, yet still contained some missingness from upstream — e.g. a GEO series with genuinely absent measurements for some probe/sample combinations). `which(is.na(expr), arr.ind = TRUE)` returns a 2-column matrix of (row, column) indices for every `NA` cell; `row_med[na_idx[, 1]]` looks up each `NA` cell's own row's median (computed once, ignoring `NA`s, via `apply(expr, 1, median, na.rm=TRUE)`) and `expr[na_idx] <- ...` writes all of them back in one vectorized assignment. **Why median, not mean, not zero:** median imputation is robust to the same skew/outlier concerns as `q99` above, and is the same technique `global.R`'s own audited `filter_and_transform_expr()` uses — the file's own comment states this explicitly, i.e. this is a deliberately consistent choice across the codebase, not a one-off. **The return value:** a `list` with the cleaned `expr`/`meta`, the source `label`, before/after sample and gene counts, and whether log2 was actually applied — this exact shape is the contract every downstream consumer (`observeEvent(input$preloaded_run)`, `merge_inputs()`, `example_live_merge()`) expects.

#### Block 6 — `mod_pp_source_ui()`/`mod_pp_source_server()` (L205–603)

Structurally, this pair implements the same three-way "upload / preloaded / current" source pattern as `pp_preloaded_read()` above, but wraps it in a full Shiny module with its own per-source sample/feature filter UI (§2.B already tabulates every input; §2.C every reactive). The parts worth teaching that are genuinely new relative to Block 5 (not simply "the same pattern again"):

```r
    output$source_type_ui <- renderUI({
      type_choices <- c("Upload files" = "upload", "A bundled cohort" = "preloaded",
                         "Currently loaded dataset" = "current")
      current <- isolate(input$source_type)
      selected <- if (!is.null(current) && current %in% type_choices) {
        current
      } else if (!is.null(default_gse)) {
        "preloaded"
      } else {
        "upload"
      }
      radioButtons(ns("source_type"), "Data source", choices = type_choices, selected = selected)
    })
```
**What's new here vs. a plain static `radioButtons()`:** this radio group is itself rebuilt via `renderUI()` rather than declared once in `mod_pp_source_ui()`, specifically so it can **remember the user's last selection** (`isolate(input$source_type)`) across whatever re-triggers this output, instead of always resetting to a hardcoded default. `isolate()` is the Shiny primitive that reads a reactive value *without* creating a dependency on it — critical here, because if this block read `input$source_type` normally, it would create a circular dependency (the output that defines `input$source_type`'s own UI depending on `input$source_type`'s current value) that would either error or infinite-loop. **Shiny behavior:** this is a defensive/UX pattern for a `renderUI` block whose own re-render could otherwise clobber user state — worth citing if the thesis discusses Shiny implementation patterns specifically.

```r
    result <- eventReactive(input$run, {
      pair <- raw_pair()
      expr <- pair$expr; meta <- pair$meta; label <- pair$label
      n_samples_before <- ncol(expr); n_genes_before <- nrow(expr)
      ...
    })
```
The full filter cascade inside this `eventReactive` — categorical filter, numeric filter, dedup, exclusion-pattern feature filter, missing-data tolerance, log2 — is taught filter-by-filter in **§2.E**, since it is functionally near-identical to (and, for the log2 step, literally the same heuristic as) `pp_preloaded_read()`'s pipeline, just parameterized by this module's own UI inputs instead of the checkbox-panel's shared controls.

### 2.E Filters in depth — Preprocessing tab

Every filter that actually exists in `mod_pp_source_server`'s `result` `eventReactive` (L519–582), in the order it is applied:

**Filter 1 — Categorical sample filter (`input$filter_col`/`input$filter_vals`)**
```r
if (!is.null(input$filter_col) && !identical(input$filter_col, "(no filter)") && length(input$filter_vals) > 0) {
  keep <- as.character(meta[[input$filter_col]]) %in% input$filter_vals
  meta <- meta[keep, , drop = FALSE]
  expr <- expr[, meta$sample, drop = FALSE]
}
```
- **What is filtered:** samples (columns of `expr`, rows of `meta`).
- **Criterion:** the chosen metadata column's value must be one of the checked values in `input$filter_vals` (a `checkboxGroupInput`, pre-populated with *every* distinct value checked by default, so the initial state is a no-op filter until the user unchecks something).
- **Threshold:** none numeric — a set-membership test.
- **Why it exists:** the most common real use is restricting to a diagnosis subset (e.g. "HC and RA only," excluding an "other"/"SLE" group present in a bundled GEO series but out of scope for a given comparison).
- **Bias risk:** entirely user-directed; the risk is a user narrowing the cohort in a way that changes what the resulting comparison actually tests, without necessarily updating every downstream label to say so — this filter changes *sample composition* invisibly to anything reading `dataset$source`'s text label unless the user also edits the per-source `label` field.
- **Validation:** compare `nrow(meta)` before/after in the browser console or by checking `result()$n_samples_before` vs. `n_samples_after` in the rendered status line.

**Filter 2 — Numeric sample filter (`input$num_filter_col`/`input$num_filter_range`)**
```r
if (!is.null(input$num_filter_col) && !identical(input$num_filter_col, "(no filter)") && !is.null(input$num_filter_range)) {
  v <- suppressWarnings(as.numeric(meta[[input$num_filter_col]]))
  keep <- !is.na(v) & v >= input$num_filter_range[1] & v <= input$num_filter_range[2]
  meta <- meta[keep, , drop = FALSE]
  expr <- expr[, meta$sample, drop = FALSE]
}
```
- **What/criterion:** samples whose value in the chosen numeric metadata column (e.g. age, RIN score) falls within an inclusive `[min, max]` range set via a slider whose own bounds are `floor(min)`/`ceiling(max)` of that column's actual observed range (`output$num_filter_range_ui`, L487–494).
- **Note:** any sample whose value is genuinely non-numeric/unparseable (`as.numeric()` → `NA`, via `suppressWarnings`) is **dropped**, not kept — the `!is.na(v)` term of `keep` excludes it. This is a real, silent-if-unnoticed data-quality gate: a metadata column with a handful of `"N/A"` text entries mixed into otherwise-numeric values would drop exactly those rows from this filter, with no separate count/warning surfaced anywhere for *this specific reason* (contrast with `expr_raw_health()`'s explicit missingness counts, which this filter doesn't call).
- **Threshold:** entirely user-set via the slider; there is no scientifically-motivated default (the slider always starts at the column's full observed range, i.e. no filtering).
- **Validation:** same before/after count comparison as Filter 1; additionally spot-check that the slider's displayed min/max genuinely match `range(meta[[col]], na.rm=TRUE)`.

**Filter 3 — Sample deduplication (`input$dedup`/`input$dedup_col`)**
```r
if (isTRUE(input$dedup) && !is.null(input$dedup_col) && input$dedup_col %in% colnames(meta)) {
  keep <- !duplicated(meta[[input$dedup_col]])
  meta <- meta[keep, , drop = FALSE]
  expr <- expr[, meta$sample, drop = FALSE]
}
```
- **What is filtered:** samples whose chosen ID column's value has already appeared earlier in the metadata table — `duplicated()` marks every occurrence **after the first** as `TRUE`, so `!duplicated(...)` keeps only first occurrences (not, e.g., an average of duplicates).
- **Why:** guards against a real GEO/upload data-quality issue — a sample accidentally listed twice under two `GSM` IDs, or a technical replicate a user doesn't want double-counted as independent biological signal.
- **Bias risk:** "keep first occurrence" is an arbitrary tie-break — if two rows genuinely represent different technical runs of the same biological sample with different quality, this filter has no way to prefer the better one; it purely depends on row order in the source file.
- **Validation:** check `sum(duplicated(meta[[dedup_col]]))` before running, and confirm it equals `n_samples_before - n_samples_after` (assuming no other filter ran in the same pass).

**Filter 4 — Feature exclusion pattern (`input$exclude_pattern`)**
```r
pattern <- trimws(input$exclude_pattern %||% "")
if (nzchar(pattern)) {
  matched <- tryCatch(grepl(pattern, rownames(expr), perl = TRUE), error = function(e) NULL)
  validate(need(!is.null(matched), paste("Invalid regex pattern:", pattern)))
  expr <- expr[!matched, , drop = FALSE]
  validate(need(nrow(expr) > 0, "The exclusion pattern matched every feature. Check the regular expression and try again."))
}
```
- **What is filtered:** features (rows of `expr`) whose row name (gene symbol or probe ID) matches a user-supplied Perl-compatible regular expression (`perl = TRUE` enables PCRE syntax, a superset of POSIX ERE — e.g. lookahead assertions, if a user needed them).
- **Criterion/threshold:** pattern match, `TRUE`/`FALSE` per feature — no numeric threshold. `nzchar(pattern)` (does the trimmed string have nonzero length) gates the whole block, so an empty pattern is a genuine no-op, not "match everything."
- **Scientific rationale:** the UI's own example, `^AFFX` (Affymetrix's convention for internal control/QC probes, not real gene measurements), is the canonical use case — control probes should never enter a differential-expression or clustering analysis as if they were biological features.
- **Robustness:** wrapped in `tryCatch` specifically because a malformed regex from free-text input would otherwise throw an uncaught R error inside a reactive, crashing that computation with an opaque message; `validate(need(...))` turns it into the same plain-language error style used everywhere else in this app. A second `validate()` guards the degenerate case of a pattern matching literally every feature (e.g. `.` un-escaped, matching any single character, thus every non-empty row name) — a real failure mode for a novice user unfamiliar with regex metacharacters.
- **Validation:** test with a known pattern against a small synthetic row-name vector and confirm `grepl()`'s output by hand.

**Filter 5 — Missing-data tolerance (`input$max_na_pct`, default 0)**
```r
na_pct <- rowMeans(is.na(expr)) * 100
expr <- expr[na_pct <= (input$max_na_pct %||% 0), , drop = FALSE]
validate(need(nrow(expr) > 0, "No features remain within the missing-data tolerance. Raise the missing-data slider and try again."))
if (anyNA(expr)) {
  row_med <- apply(expr, 1, stats::median, na.rm = TRUE)
  na_idx <- which(is.na(expr), arr.ind = TRUE)
  expr[na_idx] <- row_med[na_idx[, 1]]
}
```
- **What is filtered:** features (rows) whose fraction of samples with a missing value exceeds the slider's threshold (0–80%, step 5, default **0** — meaning, at the default, this behaves as a hard `complete.cases()` filter, dropping any feature with even a single missing sample).
- **Why a slider at all, not always hard-drop:** the file's own comment (L552–556) states this directly — this is "strictly more permissive than the old hard `complete.cases()` when the slider is above 0, identical to it at the default of 0," i.e. this is documented as a deliberate relaxation of a previously stricter rule, giving the user explicit control over the missingness/feature-retention tradeoff instead of a fixed, silent cutoff.
- **What happens after filtering:** any missingness that survives *within* the tolerance (e.g. a feature missing in 10% of samples, under a 20% tolerance) is **median-imputed**, per-gene, the same technique used throughout this module (Block 5, Filter above).
- **Could this remove biologically meaningful information?** Yes, directly and by design — a feature that is genuinely undetectable in a biologically meaningful subgroup (e.g. a transcript only expressed in one disease subtype) would show up as high missingness specifically *because* it's biologically informative, and a strict tolerance would discard exactly that feature. This is a real, disclosable tension between data-completeness convenience and true-positive preservation, worth naming explicitly in a thesis Limitations section rather than treating the default (0%) as self-evidently correct.
- **Inappropriate threshold example:** setting this to 80% (the slider's max) on a dataset with substantial platform-level missingness would let through features that are missing in 4 out of 5 samples, then median-impute the majority of that feature's own values from the minority that remain — statistically dubious (the "median" of one real value and four imputed copies of itself is just that one value, restated four times, which then enters every downstream variance-based computation as if it were four independent observations).
- **Validation:** compare `result()$n_genes_before` vs. `n_genes_after` in the status line at a few different slider settings on a test dataset with known missingness, and confirm the counts move monotonically (higher tolerance → equal or more genes retained, never fewer).

**Filter 6 — Log2 transform decision (`input$log2`)**
Same heuristic as `pp_preloaded_read()` (Block 5): `q99 > 100` on positive values, forceable/skippable by radio choice — see that block's full teaching for the mechanism. The one genuinely new detail on this path is the **default differing by source type** (L501–517):
```r
output$log2_ui <- renderUI({
  default <- if (use_upload()) "skip" else "auto"
  ...
})
```
- **Why:** the comment (L502–510) documents a real, previously-encountered bug this default change fixes — "auto" only looks at value magnitude, so a genuinely raw RNA-seq count matrix (correctly meant to stay unlogged for downstream DESeq2/TMM use) with large integer counts was being auto-log2'd, then rejected by DESeq2 for having negative/non-integer values after the fact. For the "preloaded"/"current" paths, `"auto"` remains the default (unchanged, since the bundled cohorts' scale is already known/audited); for a fresh upload, the safer default is `"skip"`, with the user free to override either way.
- **This is a real example of a threshold/heuristic being data-type-dependent** rather than one-size-fits-all, and is worth citing directly if the thesis discusses the log2 decision's implementation, since it demonstrates the app doesn't apply one blind rule everywhere.

**Distinguishing the filter categories requested by the brief, for this tab specifically:**
- **UI filters:** all six above are UI-driven (slider/checkbox/text-input controlled); none are hardcoded.
- **Data-quality filters:** Filter 3 (dedup), Filter 5 (missingness).
- **Statistical filters:** none on this tab are variance/significance-based (that's Batch Correction's job, §4.E) — this tab's filters are compositional (which samples/features to include at all) and scale-related (log2), not inferential.
- **Feature filters:** Filter 4, Filter 5, Filter 6 (log2 affects features symmetrically but isn't a subsetting filter per se).
- **Sample filters:** Filter 1, Filter 2, Filter 3.
- **Missing-data filters:** Filter 5.
- **Batch-related filters:** none on this tab — batch is only assigned a placeholder column (`meta$batch <- NA_character_` if absent) here, never filtered on.
- **User-selected filters:** all six (no automatic, non-overridable filter exists on this tab).

### 2.F Results — what the user should see

- **`output$status_ui`** (per-source path) / **`output$preloaded_status_ui`** (checkbox path): a green `✓` line reading `"{label}: {n_samples_after} of {n_samples_before} samples kept, {n_genes_after} of {n_genes_before} features kept{, log2-transformed}"`, or a red `⚠` line with the specific error message from whichever `validate()` fired first.
  - **What it means:** a direct, human-readable audit of exactly what this one preprocessing run did.
  - **Correct behavior:** the "after" counts should never exceed the "before" counts (filters only remove, never add, rows/columns); the log2 phrase should appear if and only if the auto-detect/force logic actually ran `log2()`.
  - **Bug indicator:** an "after" count *greater* than "before" would indicate a real logic error (e.g. a `merge()` accidentally duplicating rows) — should never happen given the code as read; if observed in practice, it is a genuine regression worth investigating against this exact reference (this file never introduces new rows/columns, only subsets or transforms existing ones).

### 2.G Validation checklist — Preprocessing tab

**Input validation:**
- Upload a correctly-shaped expression matrix + metadata pair → confirm `expr_raw_preview()`/`meta_raw()` parse without error and the preview text shows plausible feature/sample counts.
- Upload a metadata CSV with a missing/misnamed sample-ID column → confirm the `map_id` dropdown falls back sensibly (`guess_col()`'s `fallback = cols[1]`) rather than crashing, and that a downstream sample-ID mismatch is caught by `raw_pair()`'s own `validate(need(length(common) >= 3, ...))` guard (for the upload path).
- Upload an expression matrix in the wrong orientation (samples in rows) → expect either an outright parse failure or, more insidiously, a matrix that "works" but is scientifically nonsensical (features and samples swapped) — this is not detected anywhere in the code; a matrix-orientation sanity check (e.g. "does `colnames(expr)` look like sample IDs or like gene symbols?") is not implemented (see §7).
- Set the missing-data slider to 0% on a dataset with no missingness at all → confirm feature counts before/after are identical.

**Functional validation:**
- Click "Preprocess this dataset" (or "Load and Preprocess Selected Cohorts") and confirm `result()`/`preloaded_results()` actually populates (not `NULL`, not an error object).
- Change one filter (e.g. narrow the categorical group filter) and re-run → confirm `n_samples_after` changes accordingly and by the expected amount (cross-check against `table(meta[[filter_col]])` computed independently, e.g. in an R console on the same raw file).
- Toggle the log2 radio between "Auto," "Force," and "Skip" on the same dataset → confirm the resulting `expr` matrix's value range changes exactly as expected (e.g. max value drops by orders of magnitude under "Force" on genuinely linear-scale data).

**Scientific validation:**
- After preprocessing a known bundled cohort (e.g. GSE93272), spot-check a handful of genes' post-log2 values against an independent read of the same raw file outside the app (e.g. in a plain R script calling `get_raw_eset()`/`collapse_probes_to_genes()` directly) to confirm the app's pipeline reproduces the same numbers.
- Confirm the sample count after the "keep only HC/RA" implicit filter (Block 5's `keep <- !is.na(meta$group)`) matches the expected cohort composition documented elsewhere in the project (e.g. `Chapter_2_subchapter2_sexstratified.md`'s stated 183-sample RA-vs-HC training cohort, referenced directly in `mod_preprocessing_server`'s own comment at L713–719).

**Reproducibility validation:**
- Re-run the identical filter/log2 configuration twice in the same session → confirm byte-identical output (there is no randomness anywhere in this tab's pipeline — median imputation, filtering, and log2 are all deterministic given the same inputs, so this should always hold; if it doesn't, that is a bug).
- Confirm package versions are not silently load-bearing here — this tab uses only base R, `stats`, and `data.table::fread()`, none of which have version-sensitive default behavior likely to change a numeric result across reasonable version ranges (contrast with Batch Correction's `sva::ComBat`/`sva::sva`, which are more version-sensitive — see §4.G).

### 2.H Thesis-ready interpretation — Preprocessing tab

**Methodological description:** Each raw source dataset is individually quality-controlled before any cross-dataset merge: samples are optionally restricted by categorical diagnosis/covariate criteria and/or a numeric range, deduplicated by a chosen identifier, and features are optionally excluded by identifier pattern (e.g. array control probes) and filtered by a user-set missing-data tolerance (default: complete cases only), with any residual missingness within tolerance resolved by per-gene median imputation. A log2 transformation is applied automatically when the 99th percentile of positive expression values exceeds 100 (linear-scale heuristic), or can be forced/skipped explicitly by the user; the automatic default differs between bundled cohorts ("auto") and fresh uploads ("skip"), reflecting the greater risk of misclassifying raw RNA-seq counts as needing log2 on unfamiliar data.

**Computational implementation:** Implemented as a live Shiny reactive pipeline (`mod_preprocessing.R`, `pp_preloaded_read()` and `mod_pp_source_server`'s `result` `eventReactive`) rather than a precomputed script — every filter/transform recomputes from the raw source file(s) on each user-triggered "Preprocess"/"Load and Preprocess" click, with no caching of intermediate filtered states.

**Parameters:** Missing-data tolerance 0–80% (default 0%, step 5%); log2 mode auto/force/skip (heuristic threshold: q99 of positive values > 100); feature-exclusion pattern (free-text regex, default none); categorical/numeric sample filters (no defaults — user-set per run).

**Validation:** See §2.G's checklist; no automated unit tests exist in the codebase for this tab's functions as of this writing — validation is presently a manual, ad hoc process a user or developer performs by reading the status line and cross-checking counts.

**Expected results:** A cleaned `(expr, meta, label)` triple per source, with disclosed before/after sample and gene counts and whether log2 was applied, ready to be selected in the Merge tab.

**Limitations:** The log2 auto-detect heuristic is a single-threshold rule on the 99th percentile, not a formal statistical test, and is documented in the code itself as unable to distinguish "raw RNA-seq counts (correctly left unlogged)" from "already-normalised, genuinely-needs-logging data" by magnitude alone. Sample deduplication breaks ties arbitrarily by file row order. The missing-data tolerance filter can, at higher settings, retain features whose surviving values are majority-imputed rather than majority-measured. Three independent, byte-identical copies of the column-name-guessing helper (`pp_guess_col`, the local `guess_col` inside `mod_pp_source_server`, and `mod_dataset.R`'s own `guess_col`) exist in the codebase rather than one shared function, a maintainability risk rather than a correctness one today.

---

## PART 3 — SUBTAB: "Merge Datasets"

### 3.A Purpose

**Simple:** Once each dataset is individually cleaned, this tab combines them into one matrix so they can be analyzed together. But two datasets can only be combined on the genes they both actually measured — this tab shows you exactly which genes that is, and how much overlap there is, before doing the merge.

**Intermediate:** Different expression platforms measure different, only partially-overlapping sets of genes (a microarray only reports what's printed on its chip; an RNA-seq experiment reports whatever transcripts were detected). Merging by row name means the merged matrix can only contain genes present as a row in *every* dataset being merged — anything unique to one platform is necessarily dropped. This tab visualizes exactly how much is shared vs. lost (a Venn diagram / region-size table) before committing to the merge, so the tradeoff is visible rather than silent.

**Advanced:** This is a **feature-intersection merge**, the simplest and most conservative of several possible cross-platform integration strategies (contrast with, e.g., a union merge with `NA` for unmeasured genes, or an imputation-based cross-platform harmonization). The conservative choice trades away potentially-informative platform-unique genes for the guarantee that every retained gene has a real, directly-comparable measurement in every source dataset — no cross-platform imputation uncertainty is introduced. The tab enforces a **minimum shared-feature floor of 20** (`validate(need(length(common) >= 20, ...))`, both in the "own data" merge and the "example" rebuild path) as a sanity check against a merge that would otherwise proceed on essentially no shared measurement space — this is a hard-coded, non-adjustable threshold (see §7).

**What happens if this step is skipped or done incorrectly:** Skipping the merge and running Batch Correction directly on one un-merged dataset is fully supported (the code explicitly handles the `length(lst) == 1` case in `merged`, L1259–1264) and does nothing wrong per se, but the batch-correction step then has nothing to correct for (a single dataset has no "batch" of the relevant kind unless it has internal technical batches, e.g. separate scan dates). Doing the merge on a *wrong* feature-ID basis — e.g. merging two datasets where one uses gene symbols and the other still uses raw probe IDs, without collapsing first — would silently intersect to almost nothing (probe IDs from different platforms essentially never coincide), triggering the 20-feature floor's `validate()` error, which is itself the safety net for this exact mistake.

### 3.B UI inventory

`output$merge_tab_ui` (L829–885) is built around one top-level `radioButtons(ns("merge_mode"), ...)` choosing between two entirely different code paths, each with its own `conditionalPanel`:

| Element | Input ID | What the user sees | Controls | Reacts via |
|---|---|---|---|---|
| `radioButtons(ns("merge_mode"), ...)` | `...-merge_mode` | Two options: "Merge the example pipeline's training datasets" / "Merge your own data" | Which of two independent code branches runs | `output$merge_example_ui` vs. everything else below |
| *(mode = "example")* `checkboxGroupInput(ns("example_groups"), ...)` | `...-example_groups` | Diagnosis groups to include, defaulting to HC+RA only | `example_live_merge()`'s group filter | `example_live_merge()` |
| *(mode = "example")* `actionButton(ns("merge_use_example_btn"), ...)` | `...-merge_use_example_btn` | "Merge these datasets" | Triggers `merged` and `merge_venn_example_ui` | `merged`, `output$merge_venn_example_ui` |
| *(mode = "own")* `checkboxInput(ns("collapse_probes"), ...)` | `...-collapse_probes` | "My selected data is at probe level - collapse to one row per gene before merging" | Whether the optional probe-collapse step runs | `selected_lst()` |
| *(mode = "own", collapse on)* `fileInput(ns("collapse_annot_file"), ...)` | `...-collapse_annot_file` | Annotation-file upload | Probe→gene mapping source | `collapse_annot()` |
| *(mode = "own", collapse on)* `radioButtons(ns("collapse_method"), ...)` | `...-collapse_method` | median / maxmean / mean | Aggregation method for `pp_collapse_probes_to_genes()` | `selected_lst()` |
| *(mode = "own")* `uiOutput(ns("merge_select_ui"))` (→ `checkboxGroupInput` `merge_selected`) | `...-merge_selected` | Which preprocessed datasets to include (all checked by default) | `selected_lst()` | `overlap_sets()`, `merged` |
| *(mode = "own")* `uiOutput(ns("merge_venn_ui"))` | — | Feature-overlap Venn/region table (or a guidance message) | Read-only | `output$venn_plot_custom`, `output$venn_table_custom`, `output$venn_region_table_custom` |
| *(mode = "own")* `downloadButton(ns("download_venn_custom_png"))` / `downloadButton(ns("download_venn_custom"))` | — | PNG / CSV download of the overlap diagram/table | — | `downloadHandler` |
| *(mode = "own")* `actionButton(ns("merge_btn"), ...)` | `...-merge_btn` | "Merge datasets" | Triggers `merged` | `merged` |
| `uiOutput(ns("merge_summary_ui"))` | — | Post-merge summary: sample/feature counts, duplicate-feature warning, composition table, downloads | Read-only | `output$merge_summary_ui` |

### 3.C Function inventory — Merge Datasets tab

| Function | Type | Purpose | Input | Output | Used by |
|---|---|---|---|---|---|
| `merge_inputs` | Reactive expression | Validates every checkbox-loaded preloaded cohort succeeded, unwraps to a plain list | `preloaded_results()` | list of `(expr, meta, label)` | `output$merge_select_ui`, `selected_lst` |
| `output$merge_tab_ui` | Reactive output | Assembles the whole tab (mode radio + both conditional panels) | static markup | Shiny tags | tab 2's `tabPanel` |
| `output$merge_example_ui` | Reactive output | The "example" path's group-selection + merge button | `available_example_groups()` | Shiny tags | UI |
| `output$merge_venn_example_ui` | Reactive output | Post-click: either the fast (already-merged) summary, or a full Venn/region-table for the from-raw rebuild | `input$merge_use_example_btn`, `example_merge_from_raw()`, `example_live_merge()` | Shiny tags | UI |
| `output$merge_example_composition_table` | Reactive output | Dataset × Group sample-count table (fast path only) | `example_live_merge()` | `DT::datatable` | UI |
| `example_overlap_sets` | Reactive expression | Per-training-GSE gene-symbol sets, for the Venn | `get_collapsed_genes()` × 2 | named list of character vectors | `venn_plot_example_obj`, `venn_table_example`, `venn_regions_example` |
| `example_merge_from_raw` | Reactive expression (constant `FALSE`) | Switch: rebuild from raw probes vs. reuse the bundled merged cohort | none | boolean | `example_live_merge`, `merge_venn_example_ui` |
| `example_live_merge` | Data-processing + statistical | The actual "example" merge — either a direct read of the bundled cohort, or a from-scratch rebuild | `input$example_groups`, `example_merge_from_raw()` | `list(expr, meta, sources, n_dup_features)` | `merged`, `merge_venn_example_ui` |
| `venn_plot_example_obj`/`output$venn_plot_example`/`output$download_venn_example_png` | Plotting | Renders/downloads the training-cohort feature-overlap Venn | `example_overlap_sets()` | ggplot object / PNG file | UI |
| `output$venn_table_example` | Plotting/output | Per-dataset + common feature counts | `example_overlap_sets()` | `DT::datatable` | UI |
| `venn_regions_example`/`output$venn_region_table_example`/`output$download_venn_example` | Statistical + output | Every exact set-combination's feature count | `example_overlap_sets()` | data.frame / CSV | UI |
| `output$merge_select_ui` | Reactive output | Which-datasets-to-include checklist | `merge_inputs()` | Shiny tags | UI |
| `collapse_annot` | Reactive expression + validation | Parses the uploaded probe→gene annotation file | `input$collapse_annot_file` | data.frame | `selected_lst` |
| `selected_lst` | Reactive expression + data-processing | Applies dataset selection + optional probe collapse | `merge_inputs()`, `input$merge_selected`, `input$collapse_probes` | list of `(expr, meta, label)` | `overlap_sets`, `merged` |
| `output$merge_venn_ui` | Reactive output | Feature-overlap section (own-data path) | `selected_lst()` | Shiny tags | UI |
| `overlap_sets` | Reactive expression | Per-selected-dataset feature-name sets | `selected_lst()` | named list | `venn_plot_custom_obj`, `venn_table_custom`, `venn_regions_custom` |
| `venn_plot_custom_obj`/`output$venn_plot_custom`/`output$download_venn_custom_png` | Plotting | Own-data feature-overlap Venn | `overlap_sets()` | ggplot / PNG | UI |
| `output$venn_table_custom` | Plotting/output | Own-data per-dataset + common counts | `overlap_sets()` | `DT::datatable` | UI |
| `venn_regions_custom`/`output$venn_region_table_custom`/`output$download_venn_custom` | Statistical + output | Own-data exact-combination counts | `overlap_sets()` | data.frame / CSV | UI |
| `merged` (`eventReactive`) | Data-processing + validation | **The merge itself** — either delegates to `example_live_merge()`, or performs the intersect+cbind+rbind for the own-data path | `input$merge_mode`, `selected_lst()` or `example_live_merge()` | `list(expr, meta, sources, n_dup_features)` | Batch Correction tab's `active_meta_df`, `settings_ui`, `result` |
| `output$merge_summary_ui` | Reactive output | Post-merge summary card | `merged()` | Shiny tags | UI |
| `output$download_merged_expr`/`_meta`/`_rds` | Output/download | Export the merged matrix/metadata | `merged()` | CSV/RDS file | UI |
| `output$merge_composition_table` | Plotting/output | Dataset × group sample counts (own-data path) | `merged()` | `DT::datatable` | UI |
| `available_example_groups` | Reactive expression | Every non-missing diagnosis group actually present in the two training GEO sources | `get_raw_eset()`/`eset_harmonize_meta()` × 2 | character vector | `merge_example_ui` |

### 3.D Line-by-line teaching

#### Block 1 — `merge_inputs` (L818–827)

```r
merge_inputs <- reactive({
  res <- preloaded_results()
  validate(need(length(res) > 0, "Load at least one dataset in the Preprocessing tab first."))
  failed <- Filter(Negate(function(r) isTRUE(r$ok)), res)
  validate(need(length(failed) == 0,
                sprintf("%s failed to load: %s. Fix and re-run before merging.",
                        paste(vapply(failed, `[[`, character(1), "label"), collapse = ", "),
                        paste(vapply(failed, `[[`, character(1), "error"), collapse = "; "))))
  lapply(res, `[[`, "value")
})
```
**What:** The bridge from tab 1's checkbox-driven results into tab 2. `Filter(Negate(f), x)` is base R's "keep elements where `f` returns `FALSE`" idiom (`Negate()` wraps a predicate function so it returns the opposite boolean); here it isolates every per-cohort result that did **not** succeed (`!isTRUE(r$ok)`). The first `validate()` guards against merging with nothing loaded at all; the second refuses to proceed if *any* selected cohort failed — an **all-or-nothing gate**, not a partial-merge-with-what-succeeded fallback — with an error message listing every failed cohort's label and specific error, concatenated via `sprintf`+`paste(..., collapse=...)`. `lapply(res, `[[`, "value")` then strips away the `ok`/`error` bookkeeping and returns just the successfully-loaded `(label, expr, meta, ...)` payloads. **Why all-or-nothing:** a partial merge would silently proceed on a subset the user selected but didn't get, which is a worse failure mode than forcing the user to notice and fix the broken cohort first. **Validation:** intentionally break one cohort (e.g. temporarily rename its raw file) and confirm the merge is blocked with a message naming exactly that cohort, not a generic failure.

#### Block 2 — the "example" path: `available_example_groups`, `example_merge_from_raw`, `example_live_merge` (L720–727, 985–1096, 1012)

```r
available_example_groups <- reactive({
  grps <- unlist(lapply(PP_TRAINING_GEO_IDS, function(gse) {
    eset <- get_raw_eset(gse)
    if (is.null(eset)) return(character(0))
    unique(stats::na.omit(eset_harmonize_meta(eset, gse)$group))
  }))
  sort(unique(grps))
})
```
**What:** A cheap, metadata-only scan (no expression matrix read at all — only `get_raw_eset()`'s cached `ExpressionSet` object and its phenotype table) across both training GEO sources, collecting every diagnosis group `eset_harmonize_meta()` actually assigns (recall from §2.D Block 5 that this can be `"HC"`, `"RA"`, `"SLE"`, or `"other"`). **Why this matters:** GSE110169 carries SLE samples alongside RA/HC that GSE93272 doesn't have — this function's output (`c("HC", "RA", "SLE")`, typically) is what lets `merge_example_ui` default the group checklist to exactly `intersect(c("HC","RA"), groups)`, i.e. HC+RA only, matching the actual published training cohort composition, while still surfacing SLE as selectable for a different comparison. **Advanced/scientific point, stated directly in the code's own comment (L710–719):** this default is deliberately **named generically** (`setdiff(groups, default_groups)`, not a hardcoded `"SLE"` string) — "the exclusion itself is by group identity (anything that isn't HC/RA), not by that specific label, so this stays accurate if the training sources or their diagnosis categories ever change." This is a defensible, forward-compatible design choice worth citing if the thesis discusses maintainability.

```r
example_merge_from_raw <- reactive(FALSE)
```
**What:** A one-line reactive returning a hardcoded constant. This is the single most consequential *hidden* switch on this tab, and deserves to be flagged directly (this is a **DESIGN ISSUE**, see §7): despite its UI presenting "Merge the example pipeline's training datasets" as if it always rebuilds the merge live from raw data, this switch — currently permanently `FALSE` — makes that path instead **reuse the bundled, already-merged-and-batch-corrected cohort directly**, skipping the from-raw rebuild entirely. The code's own extensive comment (L997–1033) explains why: rebuilding from raw pays a genuinely slow one-time cost (`WGCNA::collapseRows()` over ~54,000 and ~49,000 probes) that the already-computed, already-validated merged cohort makes redundant for *this specific pair of training datasets* — but the practical consequence is that **Batch Correction, when fed via this path, has almost no real batch effect left to remove**, since the input is already corrected. The comment states this outcome explicitly and even gives the one-line code change needed to flip it back to a genuine from-raw rebuild (`!is.null(get_raw_eset(...)) && !is.null(get_raw_eset(...))`), but as shipped, this is not what runs by default.

```r
example_live_merge <- reactive({
  if (!example_merge_from_raw()) {
    d <- load_default_dataset()
    meta <- d$meta
    if (!"batch" %in% colnames(meta) || all(is.na(meta$batch))) meta$batch <- meta$dataset
    return(list(expr = d$expr, meta = meta, sources = PP_TRAINING_COHORT_LABEL, n_dup_features = 0L))
  }
  ...
```
**What (the FALSE/fast branch, which is what actually runs):** reads the bundled merged cohort via `load_default_dataset()`, then ensures a usable `batch` column exists — falling back to the per-sample `dataset` column (i.e., "which of the two source GEO series this sample came from") if no real `batch` column is present or it's entirely `NA`. Returns immediately with `n_dup_features = 0L` (no duplicate-feature accounting needed, since nothing was actually merged in this session).

```r
  parts <- lapply(PP_TRAINING_GEO_IDS, function(gse) {
    eset <- get_raw_eset(gse)
    validate(need(!is.null(eset), paste("Raw file for", gse, "not found on disk.")))
    expr <- get_collapsed_genes(gse)
    meta <- eset_harmonize_meta(eset, gse)
    wanted_groups <- input$example_groups %||% available_example_groups()
    keep <- !is.na(meta$group) & meta$group %in% wanted_groups
    meta <- meta[keep, , drop = FALSE]
    expr <- expr[, meta$sample, drop = FALSE]
    q99 <- suppressWarnings(stats::quantile(as.numeric(expr[expr > 0]), 0.99, na.rm = TRUE))
    if (isTRUE(!is.na(q99) && q99 > 100)) { expr[expr <= 0] <- NA; expr <- log2(expr); expr <- expr[stats::complete.cases(expr), , drop = FALSE] }
    if (anyNA(expr)) { row_med <- apply(expr, 1, stats::median, na.rm = TRUE); na_idx <- which(is.na(expr), arr.ind = TRUE); expr[na_idx] <- row_med[na_idx[, 1]] }
    if (!"batch" %in% colnames(meta)) meta$batch <- NA_character_
    list(expr = expr, meta = meta, label = pp_cohort_label(gse))
  })
  common <- Reduce(intersect, lapply(parts, function(x) rownames(x$expr)))
  validate(need(length(common) >= 20, "Fewer than 20 common genes between the two training datasets."))
  n_dup_features <- sum(vapply(parts, function(x) expr_raw_health(x$expr)$n_duplicated_features, integer(1)))
  merged_expr <- do.call(cbind, lapply(parts, function(x) x$expr[common, , drop = FALSE]))
  metas <- lapply(parts, function(x) { m <- x$meta; m$dataset <- x$label; m })
  all_cols <- unique(unlist(lapply(metas, colnames)))
  metas <- lapply(metas, function(m) { missing_cols <- setdiff(all_cols, colnames(m)); for (cl in missing_cols) m[[cl]] <- NA; m[, all_cols, drop = FALSE] })
  merged_meta <- do.call(rbind, metas)
  if (!"batch" %in% colnames(merged_meta) || all(is.na(merged_meta$batch))) merged_meta$batch <- merged_meta$dataset
  stopifnot(identical(colnames(merged_expr), merged_meta$sample))
  list(expr = merged_expr, meta = merged_meta, sources = PP_TRAINING_COHORT_LABEL, n_dup_features = n_dup_features)
})
```
**What (the TRUE/from-raw branch — currently unreachable via the UI, but the genuine "rebuild the pipeline live" implementation):** for each training GSE, read raw → collapse to gene symbol → filter to the wanted diagnosis groups → auto-log2 → median-impute → tag a placeholder `batch` — i.e., **exactly the same per-source pipeline as `pp_preloaded_read()`**, just inlined rather than called (a second, drift-prone copy of that logic — worth noting as a **DESIGN ISSUE**, consistent with §2.D Block 2's guess-column triplication). Then the actual merge: `Reduce(intersect, ...)` across both datasets' row names gives the common gene set (`Reduce` with a 2-argument function like `intersect` applied pairwise across a list — the general "fold" pattern, here just intersecting exactly two sets since there are always exactly two training GSEs); the 20-feature floor `validate()`; **`do.call(cbind, ...)`** column-binds the two datasets' expression matrices side-by-side, each first subset to exactly the common gene rows in the same order (`x$expr[common, , drop = FALSE]` — indexing by name guarantees row alignment regardless of each dataset's original internal row order); metadata is column-harmonized (any column present in one dataset's metadata but not the other gets filled with `NA` in the other) via a `setdiff`+fill loop, then `do.call(rbind, ...)` stacks the two metadata data.frames; a `batch` column is set to the per-sample source dataset if none exists; and a final `stopifnot(identical(colnames(merged_expr), merged_meta$sample))` is a **hard, non-recoverable assertion** (not a `validate()` — this would crash the whole reactive context with an ugly error if it ever failed) that the expression matrix's column order and the metadata's `sample` column are in lockstep, a critical invariant every downstream PCA/batch-correction/plotting call implicitly relies on.

**Why `stopifnot` here and `validate()` everywhere else:** this specific check is a genuine internal-consistency invariant that should be mathematically guaranteed by the code immediately above it (cbind preserves column order; the metadata construction preserves row order in lockstep) — if it ever fails, that means the code itself has a bug, not that the user did something wrong, so a hard crash (surfacing to a developer, not phrased as user guidance) is the more honest signal. This is a subtle but real distinction in how this codebase uses `stopifnot()` vs. `validate(need(...))`, worth noting if the thesis discusses error-handling conventions.

#### Block 3 — the "own data" path: `collapse_annot`, `selected_lst`, `overlap_sets`, `merged` (L1149–1176, 1211–1215, 1253–1287)

```r
collapse_annot <- reactive({
  validate(need(!is.null(input$collapse_annot_file), "..."))
  path <- input$collapse_annot_file$datapath
  if (grepl("\\.rds$", input$collapse_annot_file$name, ignore.case = TRUE)) {
    d <- readRDS(path); validate(need(is.data.frame(d), "...")); as.data.frame(d)
  } else { as.data.frame(data.table::fread(path, showProgress = FALSE)) }
})
```
Identical CSV/RDS-branch parsing pattern already fully taught for `meta_raw()` in §2.D Block 6 and in `mod_dataset_teaching_notes.md` — not re-taught here (see the file-level methodology note at the top of this document).

```r
selected_lst <- reactive({
  lst <- merge_inputs()
  labels <- vapply(lst, `[[`, character(1), "label")
  sel <- if (length(lst) < 2) labels else (input$merge_selected %||% labels)
  validate(need(length(sel) >= 1, "Select at least one dataset to merge."))
  lst <- lst[labels %in% sel]
  if (isTRUE(input$collapse_probes)) {
    annot <- collapse_annot()
    lst <- lapply(lst, function(x) { x$expr <- pp_collapse_probes_to_genes(x$expr, annot, input$collapse_method %||% "median"); x })
  }
  lst
})
```
**What:** Filters `merge_inputs()`'s full list down to just the checked datasets (`labels %in% sel`), then, if the collapse checkbox is on, applies `pp_collapse_probes_to_genes()` (§2.D Block 3) to **every** selected dataset **using the same single uploaded annotation file for all of them** — the UI's own note (L858) explicitly warns this is only correct "for merging same-platform sources," directing the user to collapse different-platform datasets separately before uploading if they differ. **Scientific risk if ignored:** applying one platform's probe-to-gene annotation file to a dataset actually measured on a *different* platform would silently produce a wrong (or empty, triggering the downstream 20-feature floor) mapping — nothing in this function detects a platform mismatch.

```r
overlap_sets <- reactive({
  lst <- selected_lst()
  validate(need(length(lst) >= 2, "Fewer than two datasets are selected, so there is nothing to compare."))
  setNames(lapply(lst, function(x) rownames(x$expr)), vapply(lst, `[[`, character(1), "label"))
})
```
**What:** Builds the named-list-of-character-vectors input `draw_overlap_venn()`/`overlap_region_sizes()` (both defined in `global.R`, L1829 and L1865, taught in §3.E below) expect — one entry per selected dataset, its value being that dataset's current row names (gene symbols, or probe IDs if collapsing wasn't applied). **Data flow:** `selected_lst()` (post-collapse, if applicable) → `rownames()` per dataset → named list → feeds both the visual Venn and the exact-region-count table.

```r
merged <- eventReactive(list(input$merge_btn, input$merge_use_example_btn), {
  if (identical(input$merge_mode, "example")) return(example_live_merge())
  lst <- selected_lst()
  if (length(lst) == 1) {
    x <- lst[[1]]; meta <- x$meta
    if (!"dataset" %in% colnames(meta)) meta$dataset <- x$label
    return(list(expr = x$expr, meta = meta, sources = x$label, n_dup_features = expr_raw_health(x$expr)$n_duplicated_features))
  }
  sets <- lapply(lst, function(x) rownames(x$expr))
  common <- Reduce(intersect, sets)
  validate(need(length(common) >= 20, "..."))
  n_dup_features <- sum(vapply(lst, function(x) expr_raw_health(x$expr)$n_duplicated_features, integer(1)))
  merged_expr <- do.call(cbind, lapply(lst, function(x) x$expr[common, , drop = FALSE]))
  metas <- lapply(lst, function(x) { m <- x$meta; m$dataset <- x$label; m })
  all_cols <- unique(unlist(lapply(metas, colnames)))
  metas <- lapply(metas, function(m) { missing_cols <- setdiff(all_cols, colnames(m)); for (cl in missing_cols) m[[cl]] <- NA; m[, all_cols, drop = FALSE] })
  merged_meta <- tryCatch(do.call(rbind, metas), error = function(e) { validate("Could not combine metadata across datasets. A column with the same name has a different type in different datasets, for example numeric in one and text in another. Rename or fix that column, then preprocess again.") })
  if (!"batch" %in% colnames(merged_meta) || all(is.na(merged_meta$batch))) merged_meta$batch <- merged_meta$dataset
  stopifnot(identical(colnames(merged_expr), merged_meta$sample))
  list(expr = merged_expr, meta = merged_meta, sources = paste(vapply(lst, `[[`, character(1), "label"), collapse = " + "), n_dup_features = n_dup_features)
}, ignoreInit = TRUE)
```
**What's new here vs. `example_live_merge()`'s own-data twin (Block 2):** this is a genuinely general N-dataset merge (not hardcoded to exactly 2), triggered by **either** of two different buttons (`eventReactive(list(a, b), ...)` fires whenever *either* changes — this is Shiny's way of making one reactive expression respond to multiple independent event sources) and dispatching on `input$merge_mode` right at the top to decide which of the two totally different pipelines (`example_live_merge()` vs. this own-data logic) actually runs. Two details genuinely new relative to the example path:
1. **The single-dataset special case** (`length(lst) == 1`): if only one dataset is selected, "merging" is a pass-through — no intersection, no cbind, just tagging a `dataset` column onto the metadata if one doesn't already exist, so downstream code (which always expects a `dataset` column) doesn't break. `n_dup_features` here still runs `expr_raw_health()` (a real duplicate-row-name check *within* that single dataset), even though no cross-dataset merging happened.
2. **`tryCatch` around `do.call(rbind, metas)`:** unlike the example path's unguarded `rbind`, this one specifically catches a real, disclosed failure mode — two datasets sharing a metadata column name but with incompatible types (numeric in one, character in the other), which `rbind.data.frame` would otherwise error on with a raw, unhelpful message; wrapped here into the same plain-language `validate()` style as everywhere else. This is a genuine robustness improvement the "example" path (dealing with only two known, pre-audited datasets) doesn't need but the general user-upload path does.

### 3.E Filters and thresholds in depth — Merge Datasets tab

**Filter 1 — Which datasets to include (`input$merge_selected` / `input$example_groups`)**
A pure user-selection filter, not a statistical one — every successfully-preprocessed dataset is checked by default; unchecking one removes it from both the overlap diagram and the merge itself (own-data path), or, for the example path, `input$example_groups` filters by *diagnosis group* rather than by dataset. This is the only "sample-selection-adjacent" filter on this tab (contrast with tab 1's per-sample filters, which act on individual samples within one dataset).

**Filter 2 — Feature intersection (the merge itself)**
- **What is filtered:** features (genes/probes) — kept only if present as a row name in *every* selected dataset.
- **Threshold:** none numeric per se, but a **floor of ≥20 resulting common features** is enforced (`validate(need(length(common) >= 20, ...))`), applied identically in `example_live_merge()` and `merged`'s own-data branch.
- **Why 20:** the code offers no derivation for this number — it reads as a sanity-check floor ("clearly too few to be a meaningful merge") rather than a principled statistical minimum. **This is a hard-coded, non-adjustable threshold** — flagged explicitly here as required by the audit brief (§7 **DESIGN ISSUE**): a real analysis with, say, 25 shared genes would pass this gate but might still be scientifically under-powered for many downstream uses (WGCNA module detection, robust PCA), and the app gives no further guidance once past this bar.
- **Scientific rationale:** conservative, no-imputation cross-platform merging — every retained feature has a real, independently-measured value in every dataset, at the cost of dropping every platform-unique feature.
- **Bias risk:** platforms with more comprehensive gene coverage (e.g. a newer, larger microarray, or RNA-seq's full-transcriptome view) get pulled down to the coverage of the most restrictive platform in the merge — the resulting analysis is capped at the *intersection*, not the *union*, of biological knowledge available.
- **Validation:** the Venn diagram (`draw_overlap_venn()`) and region table (`overlap_region_sizes()`) directly visualize and enumerate this filter's effect — cross-check the diagram's "common to all" region size against `length(Reduce(intersect, sets))` computed independently.
- **Reporting in a thesis:** report both the per-platform feature count *and* the final merged/common count, with the percentage retained from each source, exactly as the Venn diagram already presents it — this is directly citable, code-verified information.

**Filter 3 — Optional probe-to-gene collapse (own-data path only)**
Already taught in full in §2.D Block 3 (`pp_collapse_probes_to_genes()`); applied here, per-dataset, before the intersection filter runs, so the intersection operates on gene symbols rather than platform-specific probe IDs when this option is used — this ordering (collapse, *then* intersect) is essential, since intersecting first on raw probe IDs across different platforms would yield close to nothing in common.

### 3.F Batch correction context set up here (not yet applied)

The Merge tab does **not** itself run any batch correction — but it is the tab that determines the `batch` column Batch Correction will default to. Both merge paths set `meta$batch <- meta$dataset` **only if no real `batch` column exists or it's entirely `NA`** — i.e., "which source dataset a sample came from" becomes the batch proxy by default whenever nothing more specific (e.g. a scan-date or processing-lot column from the original metadata) was available. This is a reasonable, common convention (dataset-of-origin is often the dominant source of batch effect in a cross-cohort merge), but it is worth being explicit in a thesis Methods section that "batch," as corrected for downstream, defaults to "which dataset this sample came from" unless the uploaded/bundled metadata provided something more granular.

### 3.G Results — what the user should see

- **The Venn diagram / region table:** area-proportional circles (2–7 sets, via `ggVennDiagram`) or, beyond 7, a bar chart (this app's own configured maximum of `MAX_PP_SOURCES = 6` never reaches that fallback). **Colors** are one fixed hue per dataset from `ARTHOMIX_COLORS`; **labels** show both count and percentage per region (`label = "both"`); the **region table** is the same information as a sortable, filterable, downloadable data table — every exact "belongs to these datasets and no others" combination, largest first.
  - **Correct behavior:** the "common to all" region's size should exactly equal `length(Reduce(intersect, sets))`; per-dataset totals in the table should equal `lengths(sets)`.
  - **Bug indicator:** a region-table row with a combination naming datasets that weren't actually selected, or counts that don't sum consistently across regions, would indicate a real defect in `overlap_region_sizes()`'s combinatorial logic.
- **`output$merge_summary_ui`:** post-merge, a green line with final sample × feature counts and source list, a yellow duplicate-feature-count warning if any row-name collisions were silently resolved by "keep first occurrence" during the `x$expr[common, , drop = FALSE]` subsetting, download buttons, and a dataset × group composition table.
  - **What the duplicate-feature warning means:** row-name-keyed matrix subsetting in R (`expr[common, ]`) silently returns **only the first matching row** when `rownames()` contains duplicates — this count (via `expr_raw_health()$n_duplicated_features`, computed *before* the merge, on each input dataset) discloses how many rows were affected, since the merge itself doesn't detect or report this at the point it happens.

### 3.H Validation checklist — Merge Datasets tab

**Input validation:** merging with only one dataset preprocessed (confirm the single-dataset pass-through path in `merged`); merging with datasets on different feature-ID systems without collapsing (confirm the <20-feature floor catches it); merging with a metadata column of conflicting type across datasets in the own-data path (confirm the `tryCatch`-wrapped `rbind` error message appears rather than a raw R error).

**Functional validation:** toggling `input$merge_selected` and confirming `overlap_sets()`/the Venn update accordingly; confirming `merged()$expr`'s column count equals the sum of each included dataset's sample count, and its row count equals `length(common)` exactly.

**Scientific validation:** for the bundled example cohort, confirm the merged sample/gene counts match the project's own documented training-cohort composition (referenced in `mod_preprocessing_server`'s own comments, §2.D Block 2); confirm `example_merge_from_raw()`'s current `FALSE` setting is understood before citing "live rebuild from raw data" in a thesis methodology section — as shipped, the default example-merge path reuses the precomputed cohort, not a fresh recomputation (§3.D Block 2).

**Reproducibility validation:** re-running the identical dataset selection and collapse settings should produce byte-identical merged output (no randomness in this tab's logic).

### 3.I Thesis-ready interpretation — Merge Datasets tab

**Methodological description:** Datasets are combined via feature-ID intersection — only genes/probes present in every selected dataset are retained — with an enforced minimum of 20 common features. Sample metadata is column-harmonized across datasets (missing columns filled with `NA`) and row-bound; a `batch` column defaults to dataset-of-origin when no more specific batch variable was supplied upstream. An optional, user-invoked probe-to-gene collapse (median/maxmean/mean across probes per gene, via an uploaded annotation file) precedes the intersection when the input data is still probe-level.

**Computational implementation:** Live, on-demand recomputation (`mod_preprocessing.R`'s `merged`/`example_live_merge` `eventReactive`s) — no caching of a previous merge across sessions; every click of "Merge datasets" re-executes the full intersect/cbind/rbind pipeline.

**Parameters:** Minimum common-feature floor: 20 (hard-coded, non-adjustable). Probe-collapse method (when used): median (default) / maxmean / mean.

**Validation:** See §3.H; the Venn diagram and region table are themselves a built-in, always-visible validation artifact for the feature-intersection step specifically.

**Expected results:** A single merged `(expr, meta, sources, n_dup_features)` object whose sample count is the sum of all included datasets' samples and whose gene count is the size of their feature-name intersection (≥20 by construction).

**Limitations:** The "Merge the example pipeline's training datasets" UI option, as currently configured (`example_merge_from_raw()` hardcoded to `FALSE`), does not perform a genuine from-scratch rebuild by default — it reuses the bundled, already-batch-corrected cohort, so any Batch Correction run downstream of this specific path finds little or no real batch effect left to demonstrate correcting. The 20-feature merge floor has no stated statistical derivation. The feature-intersection strategy is conservative by construction and necessarily discards every platform-unique gene, which should be disclosed as a real information loss, not just a technical detail, whenever reporting merged-cohort gene counts.

---

## PART 4 — SUBTAB: "Batch Correction"

### 4.A Purpose

**Simple:** When you combine data from two different experiments (different labs, different days, different equipment), some of the differences you see between samples aren't biology — they're just an artifact of which batch a sample happened to be processed in. This tab tries to mathematically remove that artifact while keeping the real biological differences (e.g. disease vs. healthy) intact, and shows you, with pictures, whether it worked.

**Intermediate:** A **batch effect** is systematic, non-biological variation introduced by technical factors correlated with when/where/how a sample was processed — a different microarray lot, a different sequencing run, a different scanner calibration, a different day. Because samples from the same GEO series are usually processed together, "which dataset a sample came from" is very often *confounded* with "which technical batch it was in," which is why this tab defaults the batch variable to dataset-of-origin (§3.F). The correction methods here (ComBat, limma, SVA, ComBat-seq) work by modeling and removing batch-associated shifts in each gene's mean (and, for ComBat, variance) **while holding a set of user-chosen biological covariates fixed** — i.e., explicitly protecting real signal (disease group, sex) from being treated as if it were technical noise.

**Advanced:** This is the classic **empirical Bayes batch-correction problem** as formalized by Johnson, Li & Rabinovic (2007) for ComBat, extended by Leek et al. (2012) for SVA (unknown/latent sources of unwanted variation instead of a labeled batch), and by Zhang, Parmigiani & Johnson (2020) for ComBat-seq (a negative-binomial model appropriate for raw RNA-seq counts rather than continuous log-expression). All four methods here share the same conceptual model — partition observed variance into "of interest" (`mod`, the protected biological covariates' design matrix) and "batch" (the column(s) chosen for correction) — but differ in how they estimate and remove the batch component, and in what data scale they expect (see §4.F for the full per-method breakdown). Getting this step wrong in either direction is a real, common failure mode in multi-cohort genomics: **under-correction** leaves residual technical signal that masquerades as biological signal in any downstream test; **over-correction** (failing to protect a genuinely confounded biological covariate) can remove real disease signal along with the batch effect, especially dangerous when disease status happens to correlate strongly with batch (e.g. if one entire GEO series is disease-only and another is control-only).

**What happens if this step is skipped or done incorrectly:** Skipping it (`input$skip_combat = TRUE`) is fully supported and simply passes normalized data through unchanged — appropriate for a genuinely single-batch dataset, inappropriate (and silently misleading, since the "after" plots would then just repeat the "before" plots) for a genuinely multi-batch merged cohort. Doing it incorrectly — e.g. correcting for a batch column that is itself the variable of biological interest, or failing to protect a covariate confounded with batch — can either strip real biological signal or leave it partially intact in a way that's hard to detect after the fact; this is exactly why the tab shows PC1-vs-batch association p-values and PCA plots both before and after (§4.H), rather than asking the user to trust the correction blindly.

### 4.B UI inventory

`output$settings_ui` (L1349–1437) is the single largest UI-generating block in the module. Every input listed below only appears once `merged()` has succeeded (an early `is.null(m)` guard shows a guidance message otherwise).

| Element | Input ID | What the user sees | Controls | Reacts via |
|---|---|---|---|---|
| `selectInput(ns("color_by"), ...)` | `...-color_by` | Any merged-metadata column | Which variable colors every PCA plot | `plot_pca_advanced()` |
| `selectInput(ns("pc_x"))`/`selectInput(ns("pc_y"))` | `...-pc_x`/`...-pc_y` | PC1–PC5 dropdowns, default PC1/PC2 | Which principal-component pair is plotted | `plot_pca_advanced()` |
| `checkboxInput(ns("show_ellipse"))` | `...-show_ellipse` | Toggle, default on | Draw 68% confidence ellipses per group | `plot_pca_advanced()` |
| `checkboxInput(ns("show_labels"))` | `...-show_labels` | Toggle, default off | Label each point with its sample ID | `plot_pca_advanced()` |
| `selectInput(ns("batch_col"), ...)` | `...-batch_col` | Any merged-metadata column, defaulting to `batch`/`batch_full`/`dataset` if present | **The batch variable being corrected for** | the entire `result` pipeline |
| `radioButtons(ns("norm_method"), ...)` | `...-norm_method` | Auto-detect / Skip / Quantile / TMM+log2-CPM | Which normalization runs | `result` |
| `sliderInput(ns("min_pct"), ...)` (quantile/skip/auto branches only) | `...-min_pct` | 0–90%, step 5, default 0 | Expression-percentile gene filter | `result` |
| `radioButtons(ns("tmm_correction_stage"), ...)` (TMM branch only) | `...-tmm_correction_stage` | "After TMM" (post) / "Before TMM" (pre, ComBat-seq) | Whether batch correction runs on raw counts (ComBat-seq) or normalized log2-CPM | `result` |
| `selectInput(ns("protect_cols"), ...)`, multiple | `...-protect_cols` | Any merged-metadata column(s), defaulting to `group`/`sex` if present | Biological covariates excluded from correction | `result` (via `mod`) |
| `checkboxInput(ns("skip_combat"))` | `...-skip_combat` | "Skip batch correction (normalise only)" | Bypasses every correction method entirely | `result` |
| `checkboxInput(ns("show_advanced"))` | `...-show_advanced` | Reveals the block below | UI visibility only | — |
| `selectInput(ns("correction_method"), ...)` (advanced) | `...-correction_method` | ComBat / limma::removeBatchEffect / SVA | Which correction algorithm runs | `result` |
| `radioButtons(ns("combat_prior"), ...)` (ComBat only) | `...-combat_prior` | Parametric / Non-parametric | ComBat's empirical Bayes prior form | `run_combat()` |
| `checkboxInput(ns("combat_mean_only"), ...)` (ComBat only) | `...-combat_mean_only` | "Adjust batch mean only" | Whether ComBat also adjusts batch-wise variance | `run_combat()` |
| `uiOutput(ns("ref_batch_ui"))` (ComBat only, → `selectInput` `ref_batch`) | `...-ref_batch` | A specific batch level, or "(none)" | Reference-batch mode vs. pooled-average mode | `run_combat()` |
| `numericInput(ns("sva_n_sv"), ...)` (SVA only) | `...-sva_n_sv` | 0 (auto) or a fixed integer, 0–20 | Number of surrogate variables | `run_sva()` |
| `selectInput(ns("batch_col2"), ...)` | `...-batch_col2` | "(none)" or any column | Optional second batch column, combined into an interaction | `result` |
| `sliderInput(ns("variance_pct"), ...)` | `...-variance_pct` | 0–90%, step 5, default 0 | Additional variance-percentile gene filter | `result` |
| `checkboxInput(ns("exclude_outliers"), ...)` | `...-exclude_outliers` | "Exclude samples flagged as QC outliers before correcting" | Pre-correction outlier removal | `result` |
| `sliderInput(ns("mad_k"), ...)` | `...-mad_k` | 2–6, step 0.5, default 3 | MAD-based outlier sensitivity (used both for exclusion and post-hoc flagging) | `compute_sample_qc()` |
| `actionButton(ns("run_btn"), ...)` | `...-run_btn` | "Run normalisation and batch correction" | Triggers `result` | `result` |

Downstream, results-side outputs: `valueBoxOutput`s (`vb_samples`/`vb_genes_kept`/`vb_genes_dropped`/`vb_flagged`), `uiOutput(ns("decisions_ui"))` (plain-language pipeline summary), six `plotOutput`s (`signal_plot`/`detected_plot`/`cor_plot`/`dist_plot`/`scree_plot`/`pca_before`/`pca_after`), `uiOutput(ns("summary_ui"))` (PC1-vs-batch p-values), three `DT::dataTableOutput`s with matching `downloadButton`s (`norm_table`/`qc_table`/`pca_table`), and finally `uiOutput(ns("activate_ui"))` → `actionButton(ns("activate_btn"))` — "Use this as the active dataset app-wide."

### 4.C Function inventory — Batch Correction tab

| Function | Type | Purpose | Input | Output | Used by |
|---|---|---|---|---|---|
| `active_meta_df` | Reactive expression | Just `merged()$meta`, isolated for `settings_ui`'s column choices | `merged()` | data.frame | `settings_ui`, `ref_batch_ui`, `observeEvent` |
| `output$settings_ui` | Reactive output | The whole settings panel (§4.B) | `merged()` | Shiny tags | UI |
| `output$ref_batch_ui` | Reactive output | Reference-batch dropdown, options = observed levels of `batch_col` | `input$batch_col`, `active_meta_df()` | Shiny tags | UI |
| `observeEvent(list(input$batch_col, input$batch_col2))` | Observer | Keeps `protect_cols`' choices from ever including the currently-chosen batch column(s) | `input$batch_col`, `input$batch_col2` | side effect: `updateSelectInput` | UI consistency |
| `result` (`eventReactive(input$run_btn, ...)`) | Statistical + data-processing | **The full normalize→filter→outlier-exclude→correct pipeline** | every Settings input | large `list(...)` (§4.D) | every results-panel output |
| `output$vb_samples`/`vb_genes_kept`/`vb_genes_dropped`/`vb_flagged` | Plotting/output | Headline value boxes | `result()` | `valueBox` | UI |
| `output$decisions_ui` | Reactive output | Plain-language recap of exactly what ran | `result()` | Shiny tags | UI |
| `output$signal_plot`/`detected_plot`/`cor_plot` | Plotting | Per-sample QC bar charts (signal, detected features, cohort correlation), post-correction | `result()$qc`, `qc_bar_plot()` (defined in `global.R`) | ggplot | UI |
| `dist_summary()` | Statistical helper | Per-sample expression quantile summary (for the box plot) | `expr, meta` | data.frame | `output$dist_plot` |
| `output$dist_plot` | Plotting | Before/after normalization box plots, faceted, colored by group | `result()$expr_prenorm`/`expr_qnorm`, `dist_summary()` | ggplot | UI |
| `output$pca_before`/`pca_after`/`scree_plot` | Plotting | PCA scatter (pre/post correction) + scree plot | `result()$before`/`after`, `plot_pca_advanced()`, `scree_plot()` (global.R) | ggplot | UI |
| `assoc_pvalue()` | Statistical | One-way ANOVA p-value, PC1 ~ batch | `pc1` vector, `batch` vector | numeric p-value | `output$summary_ui` |
| `output$summary_ui` | Reactive output | Before/after PC1-vs-batch association text | `result()`, `assoc_pvalue()` | Shiny tags | UI |
| `output$norm_table`/`download_norm` | Plotting/output | Before/after normalization diagnostics table | `result()$norm_diag` | `DT::datatable`/CSV | UI |
| `qc_table_display`/`output$qc_table`/`download_qc` | Statistical + output | Flagged-samples table with reasons | `result()$qc` | data.frame/`DT::datatable`/CSV | UI |
| `pca_table`/`output$pca_table`/`download_pca` | Output | Merged before+after PCA coordinate table | `result()$before`/`after` | data.frame/`DT::datatable`/CSV | UI |
| `output$activate_ui`/`observeEvent(input$activate_btn)` | UI + observer (the module's only write to app-wide `dataset`) | "Use this as the active dataset" | `result()` | side effect: `dataset$expr/meta/source` set | app-wide |
| `bc_section()` | UI helper | Consistent section-header markup, replacing `box()` for this tab specifically | `icon_name, title, ..., desc` | Shiny tags | every results section |
| `output$results_top_ui`/`results_rest_ui` | Reactive output | Assembles the results column | `result()` | Shiny tags | UI |
| `output$batch_tab_ui` | Reactive output | Assembles the whole tab (settings + results columns) | static markup | Shiny tags | tab 3's `tabPanel` |

### 4.D Line-by-line teaching — the `result` pipeline (L1463–1715)

This `eventReactive(input$run_btn, ...)` is the scientific core of the entire module. It is taught here in the exact order it executes.

#### Step 0 — dispatch and setup (L1463–1475)

```r
result <- eventReactive(input$run_btn, {
  req(input$batch_col, input$color_by)
  {
    m <- merged()
    expr <- m$expr; meta <- m$meta; sources <- m$sources
    norm_method <- input$norm_method %||% "auto"
    skip_combat <- isTRUE(input$skip_combat)
    already_corrected <- FALSE
```
**What:** `req()` blocks execution until both a batch column and a color-by column are actually chosen (both have UI defaults, so this is mostly a startup-ordering guard, not something a user typically triggers). The rest reads `merged()` fresh at the top of every run (so a fresh merge upstream is always reflected) and establishes two flags used throughout the rest of the function: `skip_combat` (does *any* correction run at all) and `already_corrected` (set `TRUE` only by the ComBat-seq branch, since that branch produces its final corrected output directly, bypassing the standard post-normalization correction block entirely — see Step 2).

#### Step 1 — normalization branch: TMM vs. quantile/skip (L1476–1577)

This is a genuine **fork in the pipeline** depending on `input$norm_method == "tmm"` or not — the two branches operate on fundamentally different data types (raw counts vs. continuous expression) and are taught separately.

**Branch A — TMM + log2-CPM (raw RNA-seq counts):**
```r
if (identical(norm_method, "tmm")) {
  validate(need(all(expr >= 0, na.rm = TRUE), "TMM normalisation expects raw, non-negative counts, but this data has negative values, which suggests it is already log-transformed. Preprocess this dataset again with log2 set to \"Skip\"."))
  validate(need("group" %in% colnames(meta), "TMM normalisation needs a group column to filter low-count genes by."))
  grp <- factor(meta$group)
  validate(need(length(unique(na.omit(grp))) >= 2, "TMM normalisation needs at least two group levels."))
  counts <- round(as.matrix(expr)); storage.mode(counts) <- "integer"
  dge0 <- edgeR::DGEList(counts = counts)
  keepg <- edgeR::filterByExpr(dge0, group = grp)
  n_before <- nrow(expr)
  validate(need(sum(keepg) >= 50, "Fewer than 50 genes pass edgeR's expression filter for this group split."))
  counts_f <- counts[keepg, , drop = FALSE]
```
**What:** Three input-validity guards specific to count data (non-negativity, a usable `group` column, ≥2 group levels — TMM's filtering step below needs a real grouping variable), then `round()` + `storage.mode(counts) <- "integer"` **coerces the expression matrix to literal integer counts** (any residual floating-point noise from upstream processing is truncated). `edgeR::DGEList(counts = counts)` wraps this as edgeR's standard container object. **`edgeR::filterByExpr(dge0, group = grp)`** is the **low-expression gene filter** — edgeR's own recommended, group-aware heuristic (not a flat "keep genes with mean count above X" rule): it estimates, per gene, whether counts-per-million are large enough in a sufficient number of samples *within at least one group* to be considered reliably detected, adapting its effective threshold to each sample's library size. A ≥50-gene floor is then enforced post-filter — an explicit, hard-coded minimum, analogous in spirit to the Merge tab's 20-feature floor.
```r
  tmm_stage <- input$tmm_correction_stage %||% "post"
  if (!skip_combat && identical(tmm_stage, "pre")) {
    ... # ComBat-seq branch — taught in full in §4.F
  } else {
    dge <- edgeR::calcNormFactors(dge0[keepg, , keep.lib.sizes = FALSE], method = "TMM")
    expr_prenorm <- counts_f
    expr_qnorm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
    norm_label <- "TMM (edgeR::calcNormFactors) plus log2-CPM"
  }
  needs_log <- FALSE; q99 <- NA_real_; apply_qnorm <- TRUE
```
**What (the "post" sub-branch, i.e. batch-correct *after* normalizing, the more common choice):** `edgeR::calcNormFactors(..., method = "TMM")` computes **Trimmed Mean of M-values** normalization factors — a robust, per-sample scaling factor that corrects for differences in RNA composition between samples (not just total library size), the standard normalization for RNA-seq count data (Robinson & Oshlack 2010). `edgeR::cpm(dge, log = TRUE, prior.count = 1)` then converts the TMM-normalized counts to **log2 counts-per-million**, with `prior.count = 1` adding a small pseudocount before logging so that zero counts don't produce `-Inf`. **Why TMM+log2-CPM specifically, not quantile normalization, for count data:** quantile normalization assumes the underlying value distributions are already roughly continuous and comparable in shape across samples — a raw count matrix's distribution (dominated by many zero/near-zero values, discrete, highly skewed) doesn't satisfy that assumption, whereas TMM's composition-bias correction is specifically designed for count data's actual statistical properties (this exact reasoning is stated directly in the UI's own field descriptions).

**Branch B — quantile/skip/auto (continuous expression, e.g. microarray or already-log2 RNA-seq):**
```r
} else {
  n_before <- nrow(expr)
  gene_mean <- rowMeans(expr, na.rm = TRUE)
  gene_var  <- apply(expr, 1, stats::var, na.rm = TRUE)
  mean_cutoff <- stats::quantile(gene_mean, input$min_pct / 100, na.rm = TRUE)
  var_cutoff  <- stats::quantile(gene_var, (input$variance_pct %||% 0) / 100, na.rm = TRUE)
  keep <- !is.na(gene_mean) & !is.na(gene_var) & gene_mean >= mean_cutoff & gene_var > 0 & gene_var >= var_cutoff
  expr <- expr[keep, , drop = FALSE]
  validate(need(nrow(expr) >= 50, "Fewer than 50 genes remain after filtering. Lower the expression/variance percentile cutoffs."))
  expr_prenorm <- expr
  needs_log <- FALSE; q99 <- NA_real_
```
**What:** This is the **low-expression + low-variance gene filter** for continuous data — two independent percentile-based cutoffs, both user-set via sliders (`min_pct`, `variance_pct`, both default 0, meaning no filtering by default beyond the hard `gene_var > 0` zero-variance exclusion, which always applies regardless of slider position). `gene_mean >= mean_cutoff` keeps genes whose across-sample mean expression is at or above the `min_pct`-th percentile of all genes' means; `gene_var >= var_cutoff` does the same for variance. **Zero-variance genes are always excluded** (`gene_var > 0`), independent of the slider — a gene with literally no variation across the entire merged cohort carries no information for any comparative analysis and would additionally break `prcomp(scale.=TRUE)`'s per-gene scaling later in the pipeline (division by zero). Same ≥50-gene floor as the TMM branch.
```r
  diag_before <- summarize_norm_diagnostics(expr_prenorm)
  apply_qnorm <- switch(norm_method, skip = FALSE, quantile = TRUE, needs_quantile_norm(diag_before))
  expr_qnorm <- if (apply_qnorm) {
    mtx <- limma::normalizeBetweenArrays(as.matrix(expr_prenorm), method = "quantile")
    rownames(mtx) <- rownames(expr_prenorm); colnames(mtx) <- colnames(expr_prenorm)
    mtx
  } else { as.matrix(expr_prenorm) }
  norm_label <- switch(norm_method, skip = "None, used as loaded", quantile = "Quantile normalisation (forced)",
    if (apply_qnorm) "Quantile normalisation (auto-detected as needed)" else "None, auto-detected as already normalised")
```
**What:** `summarize_norm_diagnostics()` and `needs_quantile_norm()` (both `global.R`, fully described in §1's helper list and taught in detail in §4.F) compute the same before/after-comparability diagnostic used throughout the app. `switch(norm_method, skip=FALSE, quantile=TRUE, needs_quantile_norm(diag_before))` is R's `switch()` with a **fall-through default**: if `norm_method` matches neither `"skip"` nor `"quantile"` literally (i.e., it's `"auto"`), the unnamed final argument is evaluated and returned — the auto-detect heuristic. **`limma::normalizeBetweenArrays(..., method = "quantile")`** is the actual normalization when triggered: it forces every sample's expression-value distribution to have an identical set of quantiles (the classic microarray cross-sample normalization technique, Bolstad et al. 2003) — appropriate here specifically because, unlike raw counts, continuous log-scale expression values are expected to be roughly comparably distributed across samples once technical scaling differences are removed, which is exactly what forcing identical quantiles enforces.

#### Step 2 — outlier exclusion and batch correction proper (L1579–1691)

This entire block is skipped (`if (!already_corrected)`) when the TMM-branch's ComBat-seq path already ran, since that path produces its final `expr_combat` directly on raw counts, before any of this post-normalization logic applies.

```r
if (isTRUE(input$exclude_outliers)) {
  qc_pre <- compute_sample_qc(expr_qnorm, mad_k = input$mad_k)
  flagged <- qc_pre$sample[qc_pre$flag_signal | qc_pre$flag_detected | qc_pre$flag_cor]
  if (length(flagged) > 0) {
    validate(need(ncol(expr_qnorm) - length(flagged) >= 6, "..."))
    keep_samples <- setdiff(colnames(expr_qnorm), flagged)
    expr_prenorm <- expr_prenorm[, keep_samples, drop = FALSE]
    expr_qnorm   <- expr_qnorm[, keep_samples, drop = FALSE]
    meta <- meta[match(keep_samples, meta$sample), , drop = FALSE]
    n_excluded_outliers <- length(flagged)
  }
}
```
**What — the pre-correction outlier filter (optional, off by default):** `compute_sample_qc()` (`global.R` L1453–1479, taught fully in §4.F) computes three robust, MAD-based per-sample QC metrics on the *already-normalized* matrix, and any sample flagged on **any** of the three (`flag_signal | flag_detected | flag_cor`, a logical OR — one flag is enough) is dropped **before** correction runs, rather than only being flagged afterward for inspection. A ≥6-remaining-sample floor guards against excluding so many samples that too little data remains for a meaningful correction. **Why doing this before rather than after correction matters:** an extreme outlier sample can distort ComBat's per-batch mean/variance estimates for every other sample in its batch, so removing it before fitting the correction (rather than after) prevents that distortion from ever entering the model — this is a real, defensible design choice, though it does mean the "before" QC plots and PCA the user sees are computed on the *already-excluded* sample set if this option is on, not the full original cohort (worth noting explicitly in a thesis Methods section if this option was used).

```r
batch_primary <- as.character(meta[[input$batch_col]])
use_batch2 <- !identical(input$batch_col2 %||% "(none)", "(none)") && (input$batch_col2 %in% colnames(meta))
batch <- if (use_batch2) paste(batch_primary, as.character(meta[[input$batch_col2]]), sep = "_") else batch_primary
if (!skip_combat) {
  validate(need(length(unique(na.omit(batch))) >= 2, "..."))
  validate(need(all(table(batch) >= 2), "..."))
}
```
**What — building the batch label vector:** if a second batch column is chosen, the two are **concatenated into a single interaction label** (`paste(a, b, sep="_")`) — e.g. `"GSE93272_lot1"`, `"GSE93272_lot2"`, etc. — rather than modeled as two separate factors; ComBat/limma/SVA below then see this as one categorical batch variable with more, finer-grained levels. Two validity guards (skipped when correction itself is skipped): batch needs ≥2 distinct levels to correct for at all, and **every level needs ≥2 samples** — ComBat's per-batch mean/variance estimation is undefined (or degenerate) for a singleton batch.

```r
protect <- intersect(input$protect_cols %||% character(0), colnames(meta))
protect <- protect[vapply(protect, function(cl) length(unique(na.omit(meta[[cl]]))) >= 2, logical(1))]
batch_cols_used <- c(input$batch_col, if (use_batch2) input$batch_col2 else NULL)
protect_dropped_for_batch <- intersect(protect, batch_cols_used)
protect <- setdiff(protect, batch_cols_used)
mod <- if (length(protect) > 0) {
  meta_mod <- meta
  for (cl in protect) meta_mod[[cl]] <- ifelse(is.na(meta_mod[[cl]]), "Unknown", meta_mod[[cl]])
  stats::model.matrix(stats::as.formula(paste("~", paste(protect, collapse = " + "))), data = meta_mod)
} else { NULL }
```
**What — building the protected-covariate design matrix `mod`, the mechanism that keeps biological signal from being erased:** `protect` starts as the user's chosen covariates, first restricted to columns that actually have ≥2 distinct non-missing values (a constant column contributes nothing to a design matrix and would make `model.matrix()` either drop it silently or, worse, produce a rank-deficient design). Then — **this is the collinearity guard the file's own comment (L1613–1621) explains directly** — any covariate that is also (part of) the batch column being corrected for is identified (`protect_dropped_for_batch`) and **removed from the protected set** rather than left in: "Protecting the exact column being corrected for is a degenerate design - batch and mod become collinear, and depending on which fallback ComBat/limma lands on, that can either silently strip the 'protected' signal (mod gets absorbed into batch) or refuse to remove any batch effect at all." This is a genuinely important, non-obvious statistical safeguard — without it, a user could accidentally select the same column as both "batch" and "protect," producing an ill-defined correction. Missing values in a protected covariate are recoded to the literal string `"Unknown"` (rather than dropped or left as `NA`, which `model.matrix()` would otherwise drop rows for, silently shrinking the sample set) — so a sample with unknown sex, say, is protected as its own "Unknown" category rather than removed from the analysis. `model.matrix(~ protect1 + protect2, data = meta_mod)` builds the standard R design matrix (dummy/indicator columns for each factor level, an intercept column) that ComBat/limma's `mod=`/`design=` arguments expect — this is the exact mechanism, textbook linear-model design-matrix construction, by which "protect this covariate" becomes a concrete numeric object the correction algorithm can use.

### 4.E Filters in depth — Batch Correction tab

**Filter 1 — Low-expression gene filter (percentile-based, quantile/skip/auto branch): `input$min_pct`**
- **What/criterion:** genes whose across-sample mean expression is below the `min_pct`-th percentile of all genes' means.
- **Threshold:** user-set slider, 0–90%, default 0 (no-op).
- **Before/after:** disclosed via `vb_genes_dropped` (a `valueBox`) and `norm_diag`'s before/after `n_genes` rows.
- **Bias risk:** could remove low-but-genuinely-expressed genes that are nonetheless biologically important (e.g. a transcription factor with naturally low steady-state expression) — a percentile cutoff has no way to distinguish "low because unimportant" from "low because biologically rare but real."
- **Inappropriate threshold:** setting this near 90% would discard the vast majority of the transcriptome, likely well past the point of scientific usefulness for most downstream analyses (WGCNA, DGE) that expect a broad gene universe.

**Filter 2 — Low-variance gene filter (percentile-based): `input$variance_pct`**
- Same percentile mechanism as Filter 1, applied to `apply(expr, 1, var, na.rm=TRUE)` instead of `rowMeans`. Default 0 (no-op beyond the always-on zero-variance exclusion below).
- **Scientific rationale:** a gene with little to no variance across the cohort cannot contribute to any comparative statistic (differential expression, clustering, PCA loading) — filtering it out reduces the multiple-testing burden and computational cost downstream without discarding informative signal, *provided* the cutoff is conservative.
- **Bias risk:** the same as Filter 1 — a percentile-based variance filter can remove a gene that varies little overall but sharply within one biologically meaningful subgroup.

**Filter 3 — Zero-variance gene exclusion (always on, no slider): `gene_var > 0`**
- **Not user-adjustable** — applied unconditionally in the quantile/skip/auto branch. Necessary for `prcomp(scale.=TRUE)` (used in `pca_of()`) to work at all, and scientifically uncontroversial: a gene with literally zero variance carries zero statistical information in this cohort.

**Filter 4 — edgeR `filterByExpr()` (TMM/count branch only)**
- **What/criterion:** a data-adaptive, per-gene, group-aware minimum-expression rule (not a single flat threshold) — a gene is kept if its counts-per-million are large enough, in enough samples within at least one group, given each sample's library size.
- **Threshold:** edgeR's own internal default parameters (not exposed in this UI at all) — this is a real "hidden threshold" (flagged per the audit brief's list of things to look for) in the sense that a user cannot see or adjust edgeR's internal CPM/sample-count cutoffs from this tab, only observe the resulting gene count via `vb_genes_dropped`.
- **Floor:** ≥50 genes must survive, hard-coded.
- **Why group-aware matters:** a flat "average CPM above X" rule would systematically under-detect genes expressed specifically in a minority group (e.g. a small disease subgroup) if that group's average gets diluted by a larger, non-expressing group — edgeR's group-aware filter avoids this by checking "enough samples *within at least one group*," not "enough samples overall."

**Filter 5 — Pre-correction outlier sample exclusion: `input$exclude_outliers` + `input$mad_k`**
- **What is filtered:** samples (columns), flagged by `compute_sample_qc()`'s three robust MAD-based metrics (§4.F), before correction runs.
- **Threshold:** `mad_k` — how many median absolute deviations from the cohort median counts as "outlier," 2–6, step 0.5, default 3 (a conventional robust-statistics choice; 3 MADs is a common convention, though not universal — the file offers no citation for this specific default).
- **Floor:** ≥6 samples must remain after exclusion.
- **Off by default:** this filter does not run unless explicitly enabled — the default pipeline flags outliers post-hoc (in the QC table) but does not remove them.
- **Bias risk:** a real biological outlier (e.g. a genuinely unusual but correctly-measured patient) is indistinguishable, by this purely statistical criterion, from a technical artifact — removing it before correction assumes it's the latter.

**Filter 6 — protected-covariate/batch-column exclusivity (`protect_dropped_for_batch`)**
Not user-adjustable directly — an automatic consequence of the batch-column and protect-column choices, described in full in §4.D above. Included here because it *is* a real filter on the *set of protected covariates*, silently narrowing `protect` whenever it overlaps `batch_cols_used`, and disclosed to the user via `decisions_ui`'s warning line rather than left invisible.

### 4.F Batch correction methods — exact implementation, one by one

**What is a batch, concretely, in this codebase:** whatever column the user assigns to `input$batch_col` (optionally combined with a second column into an interaction label via `input$batch_col2`) — defaulting to `batch`/`batch_full`/`dataset` if one of those exists in the merged metadata, and to dataset-of-origin specifically if the Merge tab set `meta$batch <- meta$dataset` (§3.F). **Which variables represent biological conditions:** whatever the user selects in `input$protect_cols`, defaulting to `group`/`sex` if present. **How biological variables are protected:** via the `mod` design matrix (§4.D), passed to every correction method's `mod=`/`design=`/`covar_mod=` argument so the algorithm can distinguish "this shift is associated with the protected covariate, leave it" from "this shift is associated with batch, remove it."

**Method 1 — ComBat (default, `sva::ComBat`)**
```r
run_combat <- function(b, use_mod = TRUE, use_ref = TRUE) {
  sva::ComBat(dat = expr_qnorm, batch = b, mod = if (use_mod) mod else NULL,
              par.prior = identical(combat_prior, "param"), mean.only = combat_mean_only,
              ref.batch = if (use_ref) ref_batch else NULL)
}
...
expr_combat <- ... tryCatch(run_combat(batch, use_mod = TRUE, use_ref = TRUE),
  error = function(e) tryCatch(run_combat(batch_primary, use_mod = TRUE, use_ref = FALSE),
    error = function(e2) run_combat(batch_primary, use_mod = FALSE, use_ref = FALSE)))
```
- **Package/function:** `sva::ComBat` (Johnson, Li & Rabinovic, *Biostatistics* 2007) — an **empirical Bayes** method for removing batch effects from microarray/continuous expression data.
- **Input matrix:** `expr_qnorm` — the already-filtered, already-normalized matrix (genes × samples), on the log/continuous scale ComBat's Gaussian model assumes.
- **Batch variable:** `batch` (or `batch_primary` on fallback — see below).
- **Covariates:** `mod`, the protected-covariate design matrix (or `NULL` on fallback).
- **Output:** a batch-adjusted matrix, same dimensions.
- **Statistical assumptions:** per-gene, per-batch location (and, unless `mean.only`, scale) shifts are modeled with an empirical Bayes prior that borrows strength across genes — this is specifically why ComBat performs well even with small per-batch sample sizes, at the cost of assuming batch effects are broadly similar in *form* across genes (the prior's whole point is shrinkage toward a common estimate).
- **Two tunable knobs surfaced in this UI:**
  - `par.prior` (`combat_prior`, "Parametric" default vs. "Non-parametric"): the parametric prior is faster and is ComBat's default; the non-parametric prior is described in the UI as "slower, more robust for small or uneven batches" — this matches the original ComBat paper's own guidance.
  - `mean.only` (`combat_mean_only`, off by default): when on, ComBat adjusts only each batch's mean shift, leaving per-batch variance untouched — a real, disclosed choice for when a user believes batch effects are location-only, not scale.
- **Reference batch (`ref.batch`, optional):** when set, every other batch is shifted **to match this one specific batch's own distribution**, rather than to a pooled cross-batch average — useful when one batch is considered the "gold standard" (e.g. the platform an established reference cohort was measured on).
- **The fallback ladder (`tryCatch` nested three deep):** this is a real, disclosed robustness mechanism, not silent degradation — if the full call (with `mod` and `ref.batch`) fails (a common ComBat failure mode: a covariate combination that creates a singular/non-invertible per-batch design), the code retries first without `ref.batch`, then without `mod` entirely. **This exactly mirrors** the project's own upstream pipeline script (`scripts/00_shared/03_normalize_batch.R`'s `run_combat()`, per the code's own comment at L1682–1684) — i.e., this is not an ad hoc addition but a reproduction of an already-used pattern. **Scientific caveat worth stating directly:** if the fallback ladder is actually triggered, the correction that ran is *not* the one the user configured (e.g. it ran without protecting the chosen covariates) — the UI's `decisions_ui` panel reports what was ultimately used (`res$protect`, `res$combat_prior`, etc.), so this is disclosed, but a careless reading of only the "Correction method: ComBat" headline without checking the "Protected: ..." line could miss that a fallback occurred.

**Method 2 — limma::removeBatchEffect**
```r
run_limma <- function(b) {
  design <- if (!is.null(mod)) mod else matrix(1, ncol(expr_qnorm), 1)
  limma::removeBatchEffect(expr_qnorm, batch = b, design = design)
}
```
- **Package/function:** `limma::removeBatchEffect` (Ritchie et al. 2015, the `limma` package).
- **Statistical model:** a **simple linear-model adjustment** — fits and subtracts a batch-associated linear effect per gene, given a `design` matrix specifying what to protect (falling back to a plain intercept-only matrix, i.e. no explicit protection, if no covariates are selected — this is a real difference from ComBat's own behavior: `removeBatchEffect` has no fallback ladder or empirical Bayes shrinkage; it is the most direct, least assumption-laden of the four options).
- **When preferred over ComBat:** the UI's own label calls this "simple linear adjustment" — no borrowing-of-strength across genes, no distributional prior, so it is a more transparent but potentially less stable choice for datasets with small per-batch sample counts (where ComBat's empirical Bayes shrinkage genuinely helps).

**Method 3 — SVA (Surrogate Variable Analysis, `sva::sva`/`sva::num.sv`)**
```r
run_sva <- function() {
  mod_full <- if (!is.null(mod)) mod else matrix(1, ncol(expr_qnorm), 1)
  mod0 <- matrix(1, ncol(expr_qnorm), 1)
  n_genes <- nrow(expr_qnorm); vfilt <- if (n_genes > 2000) 2000L else NULL
  n_sv <- as.integer(input$sva_n_sv %||% 0)
  if (n_sv <= 0) n_sv <- tryCatch(sva::num.sv(as.matrix(expr_qnorm), mod_full, method = "be", vfilter = vfilt), error = function(e) NA_integer_)
  n_sv <- if (is.na(n_sv)) 1L else max(1L, min(n_sv, ncol(expr_qnorm) - ncol(mod_full) - 1L, 20L))
  sv_obj <- sva::sva(as.matrix(expr_qnorm), mod_full, mod0, n.sv = n_sv, vfilter = vfilt)
  validate(need(sv_obj$n.sv >= 1, "..."))
  limma::removeBatchEffect(expr_qnorm, covariates = sv_obj$sv, design = mod_full)
}
```
- **Package/functions:** `sva::num.sv` (estimate how many surrogate variables to extract), `sva::sva` (extract them), `limma::removeBatchEffect` (regress them out) — Leek, Johnson, Parker, Jaffe & Storey, *Bioinformatics* 2012.
- **What makes SVA fundamentally different from the other three methods:** SVA **does not take a batch label at all**. It estimates latent sources of unwanted variation directly from the residual structure of the data itself, after protecting the covariates in `mod_full` — i.e. it is the tool of choice specifically when "the real source of batch effects is unknown or only partly captured by a column you have" (the UI's own description, L1423). This is the method to reach for when a labeled `batch`/`dataset` column doesn't fully capture the technical variation actually present (e.g. unrecorded processing-day effects within one nominal "batch").
- **`vfilter` (feature-count cap for SV estimation, 2000 if more genes are present):** both `num.sv()`'s permutation-based `"be"` method and `sva()`'s own eigendecomposition are computationally expensive per iteration and scale with gene count; the code's own comment states this cap follows "the sva package vignette for large expression matrices" — a documented, not invented, performance practice.
- **`n_sv` bounds:** `max(1L, min(n_sv, ncol(expr_qnorm) - ncol(mod_full) - 1L, 20L))` — at least 1, at most 20, and never more than the degrees of freedom actually available (`n_samples - n_protected_covariate_columns - 1`), a real statistical necessity (SVA cannot meaningfully extract more surrogate variables than there is residual degrees of freedom to support).
- **Final step:** `limma::removeBatchEffect(expr_qnorm, covariates = sv_obj$sv, design = mod_full)` — this is the sva package's own documented recipe (the code's comment cites this directly) for producing SVA-adjusted data for downstream visualization/use: once the surrogate variables are estimated, they're regressed out exactly the way a *known* batch label would be, via the same `removeBatchEffect` function Method 2 uses directly.

**Method 4 — ComBat-seq (`sva::ComBat_seq`, raw-count branch, "pre-TMM" stage only)**
```r
counts_adj <- tryCatch(
  sva::ComBat_seq(counts = counts_f, batch = cs_batch, group = as.character(grp), covar_mod = cs_covar_mod),
  error = function(e) sva::ComBat_seq(counts = counts_f, batch = cs_batch_primary, group = as.character(grp))
)
dge_before <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts_f), method = "TMM")
dge_after  <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts_adj), method = "TMM")
expr_prenorm <- counts_f
expr_qnorm  <- edgeR::cpm(dge_before, log = TRUE, prior.count = 1)
expr_combat <- edgeR::cpm(dge_after,  log = TRUE, prior.count = 1)
```
- **Package/function:** `sva::ComBat_seq` (Zhang, Parmigiani & Johnson, *Briefings in Bioinformatics* 2020, cited directly in the code).
- **What's fundamentally different:** this is the **only** method here that operates on **raw counts**, using a **negative-binomial** model appropriate for count data, rather than on continuous/log-scale expression — chosen specifically for "batch effects that are themselves count-scale (library-size/depth-driven) rather than effects that only show up after log-CPM" (the code's own comment). It runs **before** TMM normalization, not after — the UI's radio choice explicitly frames this as "Before TMM: ComBat-seq on raw counts... better for large, count-driven batch effects" vs. the standard "After TMM: ComBat or limma on log2-CPM."
- **Batch/group/covariates:** `batch = cs_batch` (built the same interaction-label way as the standard path); `group = as.character(grp)` — **ComBat-seq always protects the group column directly** as a required argument (the UI's own note states this: "It always protects the group column directly"), distinct from the standard path's optional, user-chosen `protect_cols`; `covar_mod` — a `model.matrix()` built from whatever *other* protected covariates remain after removing `group` and the batch columns from the protect set.
- **Fallback:** a single-level `tryCatch` (not the 3-deep ladder ComBat has) — if the full call with `covar_mod` fails, retry with just `batch`/`group`.
- **After correction:** the corrected counts are **then** TMM-normalized and log2-CPM'd (`dge_after`), and — importantly — the "before" comparison (`expr_qnorm`, used for every "before batch correction" plot/statistic elsewhere in the pipeline) is **TMM-normalized *uncorrected* counts**, not raw counts, so the before/after comparison is apples-to-apples on the same final scale (log2-CPM), differing only in whether ComBat-seq ran beforehand.
- **What ComBat-seq ignores, per the UI's own explicit disclosure:** "ComBat-seq (before TMM) ignores the correction method, prior, reference batch and exclude-outliers options below" — i.e. selecting this path silently overrides several other Settings inputs, a real UI/logic coupling worth being aware of before citing "the batch correction settings used" in a thesis without checking which path actually ran.

**Summary table — method selection logic:**

| Data type | Correction timing | Method | Trigger |
|---|---|---|---|
| Continuous/log-scale (microarray, normalized RNA-seq) | After normalization | ComBat (default) | `norm_method != "tmm"`, `correction_method == "combat"` (or unset) |
| Continuous/log-scale | After normalization | limma::removeBatchEffect | `correction_method == "limma"` |
| Continuous/log-scale | After normalization | SVA (unknown/latent sources) | `correction_method == "sva"` |
| Raw RNA-seq counts | After TMM (on log2-CPM) | ComBat / limma / SVA, same as above | `norm_method == "tmm"`, `tmm_correction_stage == "post"` |
| Raw RNA-seq counts | Before TMM (on raw counts) | ComBat-seq | `norm_method == "tmm"`, `tmm_correction_stage == "pre"` |

**Normalization diagnostics — `summarize_norm_diagnostics()`/`needs_quantile_norm()` (global.R L1550–1564), used throughout §4.D/4.F:**
```r
summarize_norm_diagnostics <- function(m) {
  sample_medians <- apply(m, 2, stats::median, na.rm = TRUE)
  sample_iqr <- apply(m, 2, stats::IQR, na.rm = TRUE)
  data.frame(n_samples=ncol(m), n_genes=nrow(m), max_value=max(m, na.rm=TRUE), min_value=min(m, na.rm=TRUE),
             median_sd=stats::sd(sample_medians, na.rm=TRUE), iqr_sd=stats::sd(sample_iqr, na.rm=TRUE),
             median_range=diff(range(sample_medians, na.rm=TRUE)), iqr_range=diff(range(sample_iqr, na.rm=TRUE)))
}
needs_quantile_norm <- function(diag) diag$max_value > 100 || diag$median_sd > 0.5 || diag$iqr_sd > 0.5
```
**What/why:** computes, per sample, the median and IQR of its own expression values, then summarizes the **spread of those per-sample summaries across the whole cohort** (`sd`/`range` of the per-sample medians, and separately of the per-sample IQRs) — a direct, quantitative measure of "how differently-distributed are these samples from each other." `needs_quantile_norm()` is the auto-detect decision rule used in Step 1: normalize if the raw max value is implausibly large for log-scale data (`> 100`, the same threshold family as the log2 heuristic elsewhere in this module) **or** if per-sample medians/IQRs disagree by more than 0.5 (on whatever scale the data is currently on). **This exact rule is stated in the code's own comment (L1545–1549) to match `scripts/00_shared/03_normalize_batch.R`'s own thesis-pipeline logic** — i.e., this heuristic is not invented for the live Shiny tool; it reproduces the upstream, already-used pipeline decision rule, which is directly relevant if the thesis needs to argue the live and precomputed tracks use consistent logic.

**Sample-level QC — `compute_sample_qc()` (global.R L1453–1479), used in Step 2's outlier exclusion and throughout the QC plots/table:**
```r
compute_sample_qc <- function(expr, mad_k = 3, top_n_cor = 2000) {
  detect_cutoff <- stats::quantile(expr, 0.25, na.rm = TRUE)
  signal   <- colSums(expr, na.rm = TRUE)
  detected <- colSums(expr > detect_cutoff, na.rm = TRUE)
  gene_var <- apply(expr, 1, stats::var, na.rm = TRUE)
  top_idx  <- order(gene_var, decreasing = TRUE)[seq_len(min(top_n_cor, nrow(expr)))]
  sample_cor <- cor(expr[top_idx, , drop = FALSE])
  mean_cor <- (colSums(sample_cor) - 1) / (ncol(sample_cor) - 1)
  is_outlier <- function(x, low_only = FALSE) {
    m <- stats::median(x); s <- stats::mad(x)
    if (s == 0) return(rep(FALSE, length(x)))
    if (low_only) (m - x) > mad_k * s else abs(x - m) > mad_k * s
  }
  data.frame(sample=colnames(expr), signal=signal, detected=detected, mean_cor=mean_cor,
             flag_signal=is_outlier(signal), flag_detected=is_outlier(detected), flag_cor=is_outlier(mean_cor, low_only = TRUE))
}
```
**Three independent QC metrics, each a genuinely different signal:**
1. **`signal`** — total summed expression per sample (`colSums`). An unusually high or low total (flagged both directions, `abs(x - m) > mad_k * s`) can indicate a scanning/library-prep problem.
2. **`detected`** — count of features above the cohort's 25th-percentile value (`detect_cutoff`), a proxy for "how many genes registered a real signal in this sample" (flagged both directions).
3. **`mean_cor`** — each sample's average pairwise correlation to every other sample, computed **only on the top 2,000 most-variable genes** (`top_n_cor`, restricting to genes most likely to carry real biological signal rather than noise) — `(colSums(sample_cor) - 1) / (ncol(sample_cor) - 1)` excludes each sample's self-correlation (always 1) from its own average. Flagged **only when low** (`low_only = TRUE`, i.e. `(m - x) > mad_k * s`, not `abs(...)`) — a sample unusually *dissimilar* to the rest of the cohort is a real outlier signal; a sample unusually *similar* to everyone else is not a problem worth flagging.
- **`is_outlier()`'s robust-statistics mechanism:** median + MAD (median absolute deviation) rather than mean + SD — the standard robust alternative, far less sensitive to the very outliers being detected (an extreme value inflates a mean/SD estimate far more than a median/MAD one, which is precisely the failure mode a naive mean-based outlier detector would have). `if (s == 0) return(rep(FALSE, ...))` guards the degenerate case of zero MAD (e.g. every sample has an identical `signal` value) rather than dividing by zero.
- **Where this exact function is reused:** it is not local to Batch Correction — it is also the mechanism behind `mod_dataset.R`'s Dataset-tab QC (per that file's own teaching notes) and this file's own EDA tab (`eda_sample_outliers()`, wrapping this same function — §5.D), i.e. one shared, audited implementation used consistently across the app rather than three independent reimplementations, in contrast to the column-guessing helper's actual triplication (§2.D Block 2).

### 4.G Validation checklist — Batch Correction tab

**Input validation:** run with a batch column that has only 1 level (confirm the `≥2 levels` `validate()` fires); run with a batch level containing only 1 sample (confirm the `all(table(batch) >= 2)` guard fires); select the same column for both `batch_col` and `protect_cols` (confirm the `observeEvent` keeps it out of `protect_cols`'s choices, and that the defensive `setdiff()` inside `result()` still holds even if that observer somehow hasn't fired yet); run TMM normalization on already-log2 (negative-containing) data (confirm the non-negativity `validate()` fires with the exact guidance to re-preprocess with log2 set to "Skip").

**Functional validation:** confirm `result()$expr_prenorm`, `expr_qnorm`, and `expr_combat` all have the same row count (post-filter gene count) and, unless outlier exclusion ran, the same column count as `merged()$expr`; confirm `vb_genes_dropped` equals `n_before - n_after` computed independently; confirm toggling `skip_combat` on makes `pca_before`/`pca_after` visually and numerically identical (`before$df` should equal `after$df` exactly, since `expr_combat <- expr_qnorm` when skipped).

**Scientific validation — does correction actually reduce the batch effect while preserving biology:**
- **The tab's own built-in check:** `assoc_pvalue()`'s one-way ANOVA of PC1 against the batch column, before vs. after — a **larger** p-value after correction (less of PC1's variance explained by batch) is the expected signature of successful correction. Cross-check this by computing the same ANOVA independently outside the app on the returned `pca_table` CSV download.
- **Confirm biology is preserved, not just batch removed:** run the same before/after PCA colored by the protected covariate (`color_by` set to `group`, say) rather than by batch — group separation should be visually comparable (not diminished) before and after correction, if the protected-covariate mechanism worked as intended. The app does not compute a second ANOVA (PC1-vs-group) automatically; this is a manual cross-check worth performing and reporting.
- **Are important features being removed unexpectedly:** check `vb_genes_dropped` against the specific filter percentiles chosen; a much larger drop than the chosen percentile sliders would suggest implies the `filterByExpr()`/zero-variance filters are doing more work than the visible sliders alone (expected — they are additional, not user-adjustable, filters, per §4.E).

**Reproducibility validation:** ComBat's parametric-prior path (`par.prior = TRUE`) is deterministic given identical input; SVA's `sva::sva()` call, however, **is not guaranteed bit-for-bit deterministic across package versions or even across runs**, since its iterative surrogate-variable estimation can depend on algorithmic details not fully pinned by a fixed seed in this code (no `set.seed()` call appears anywhere in this file) — **this is a real reproducibility caveat worth disclosing explicitly** if SVA is used and cited in a thesis: re-running the same SVA correction on the same data may not produce bit-identical output, only closely comparable output, unless the exact package version and any implicit seed state are also controlled for. Package versions are potentially load-bearing for `sva::ComBat`/`sva::ComBat_seq`/`sva::sva` specifically (these functions have had documented behavioral changes across Bioconductor releases); this is not something the code itself pins or checks.

### 4.H Results — what the user should see

- **`vb_samples`/`vb_genes_kept`/`vb_genes_dropped`/`vb_flagged`:** four headline numbers. Correct behavior: `genes_kept + genes_dropped == n_before` (the pre-filter gene count); `vb_flagged` should be 0 or small for a well-behaved cohort, and should be 0 always for TMM's branch if `n_excluded_outliers > 0` already removed the flagged samples pre-correction (they'd no longer be present to flag again in the post-correction QC).
- **`decisions_ui`:** a plain-language recap — normalization method, correction method (with its specific sub-configuration: prior, mean-only, batch2 interaction, reference batch), protected covariates (or an explicit statement that none were protected/none qualified), any covariate dropped for batch-column overlap, and how many outliers were excluded pre-correction. **This is the single most important text block for thesis reporting** — it states in one place exactly what ran, generated fresh from the actual `result()` object, not a static description.
- **Per-sample bar plots (signal/detected/correlation):** x-axis = sample (unlabeled, ordered), y-axis = the metric; flagged samples should visually stand out; **expected pattern:** most samples clustered near a common level, with genuine outliers as the visually distinct minority — a pattern where *most* samples are "flagged" would indicate the MAD-based threshold itself, or `mad_k`, needs reconsidering, not that most of the cohort is actually anomalous.
- **`dist_plot` (before/after normalization boxplots):** x-axis = sample (unlabeled), y-axis = expression value, faceted "Before normalisation" / "After normalisation", fill = group. **Expected pattern after normalization:** box heights (medians/IQRs) should visibly align across samples; **unexpected/bug-indicating pattern:** boxes remaining wildly different heights after normalization was actually applied (not skipped) would suggest either normalization didn't run as configured, or the underlying data has some more fundamental issue (e.g. one dataset with a corrupted subset of values).
- **`pca_before`/`pca_after`:** x/y = the chosen PC pair, color = `color_by`, optional 68%-confidence ellipses per group (only drawn for groups with ≥4 samples — `stat_ellipse()` requires enough points to estimate a covariance ellipse meaningfully) and optional sample-ID labels. **Expected pattern when coloring by batch:** batch-driven clustering visible "before," reduced or absent "after." **Expected pattern when coloring by a protected biological covariate:** any real separation present "before" should remain visible "after" — its *disappearance* would be the visual signature of over-correction.
- **`scree_plot`:** bar chart, one bar per PC (up to 10), % variance explained, computed **after** batch correction. A sharp drop after PC2–PC3 (the plot's own description) indicates most systematic structure is captured by the low-order components actually being plotted.
- **`summary_ui`:** the PC1-vs-batch ANOVA p-value pair (§4.G) — the tab's single quantitative correction-effectiveness statistic.
- **`norm_table`:** two rows ("before normalisation"/"after normalisation"), columns = `n_samples, n_genes, max_value, min_value, median_sd, iqr_sd, median_range, iqr_range` — the exact `summarize_norm_diagnostics()` output, before/after. Large "after" `median_sd`/`iqr_sd` values indicate normalization did not fully align the samples (the same interpretation guidance the table's own description gives).
- **`qc_table`:** only flagged samples shown by default (an empty table, `df[0,]`, if none are flagged) — sample, group, signal, detected, correlation, and a `reason` column built by `apply()`ing over the three flag columns and pasting together whichever human-readable reasons apply.
- **`pca_table`:** every sample's PC1–PC5 coordinates, before and after, suffixed accordingly — the rawest, most directly reusable artifact for independent statistical cross-checking outside the app.

### 4.I Thesis-ready interpretation — Batch Correction tab

**Methodological description:** The merged expression matrix is filtered (percentile-based expression/variance cutoffs for continuous data, or edgeR's group-aware `filterByExpr()` for raw counts, always excluding zero-variance genes), then normalized (quantile normalization via `limma::normalizeBetweenArrays` for continuous/log-scale data, or TMM via `edgeR::calcNormFactors` + log2-CPM for raw counts, auto-detected via a per-sample median/IQR agreement heuristic unless forced). Optionally, samples flagged by robust (median/MAD-based) signal, detection-rate, or cohort-correlation outlier criteria are excluded before correction. Batch effects are then removed using one of four methods selected by the user — ComBat (empirical Bayes, parametric or non-parametric prior, optional mean-only mode, optional reference-batch anchoring), limma's `removeBatchEffect` (direct linear adjustment), SVA (latent/unknown source estimation, no batch label required), or ComBat-seq (negative-binomial correction on raw counts, applied before TMM normalization for count-driven batch effects) — each explicitly protecting a user-chosen set of biological covariates via a `model.matrix()`-derived design matrix, with an automatic safeguard preventing a covariate from being both the batch variable and a protected covariate simultaneously.

**Computational implementation:** `mod_preprocessing.R`'s `result` `eventReactive`, wrapping `sva::ComBat`/`sva::ComBat_seq`/`sva::sva`/`limma::removeBatchEffect`/`limma::normalizeBetweenArrays`/`edgeR::calcNormFactors`/`edgeR::filterByExpr`/`edgeR::cpm` — no reimplementation of any correction mathematics; every call is a direct, documented invocation of the corresponding Bioconductor package function.

**Parameters:** Expression/variance percentile filters (0–90%, default 0); MAD outlier sensitivity (2–6, default 3); ComBat prior (parametric default); mean-only toggle (off by default); optional reference batch (none by default); optional second (interaction) batch column (none by default); SVA surrogate-variable count (0 = auto-estimated via `sva::num.sv(method="be")`, capped at 20).

**Validation:** See §4.G; the tab's own before/after PC1-vs-batch ANOVA is the built-in quantitative correction-effectiveness check, and should be supplemented (manually, outside the app) with a group-vs-PC1/PC2 check to confirm biological signal was preserved, not just that batch signal was removed.

**Expected results:** A normalized, batch-corrected expression matrix; before/after PCA and normalization-diagnostic tables; a per-sample QC flag table; and, on confirmation, promotion to the app-wide active dataset for every other transcriptomics submodule.

**Limitations:** The low-expression/low-variance gene filters are percentile-based, not informed by any formal statistical test of expression reliability, and can remove biologically real but rare/subtle signal. `edgeR::filterByExpr()`'s internal CPM/sample-count thresholds are not exposed or adjustable in this UI. The ComBat correction path includes an automatic fallback ladder that can silently run without the user's chosen `mod`/`ref.batch` if the fully-specified call fails — disclosed via `decisions_ui`, but easy to miss if that panel isn't read carefully. SVA's estimation is not guaranteed bit-for-bit reproducible across runs or package versions, and no explicit random seed is set anywhere in this pipeline. ComBat-seq silently overrides several other Settings inputs (correction method, prior, reference batch, outlier exclusion) when selected, a real UI/logic coupling that should be disclosed whenever this specific path is used and reported.

---

## PART 5 — SUBTAB: "Data Exploration" (`mod_preprocessing_explore.R`)

### 5.A Purpose

**Simple:** Before you trust any dataset enough to preprocess or analyze it, it helps to just *look* at it first — how many genes and samples, how much is missing, does it look already log-transformed, are there any obviously weird samples. This tab is a self-contained "look at my raw file first" tool, completely separate from everything else on this page.

**Intermediate:** This is a classic **exploratory data analysis (EDA)** workflow, applied specifically to omics-scale feature-by-sample matrices: descriptive statistics, distributional diagnostics, an automated normalization-status assessment, unsupervised structure (PCA, sample correlation/clustering), outlier flagging, missing-data and low-variance summaries, and a mean-variance relationship check — each with its own visualization and a written interpretation, ending in one synthesized plain-language verdict and a recommended next step. Unlike tabs 1–3, this tab **never writes anything back** — its own UI text states this explicitly, twice (module-level intro text and the "Diagnostic only" note on the transformation-comparison section).

**Advanced:** Architecturally, this is the one tab in the module that is **not** part of the Filter→Merge→Normalize→Batch-correct pipeline and does not share the `dataset` reactiveValues at all — it has its own `fileInput`, its own `eventReactive` pipeline (`eda_result`), and touches no shared app state. This is a deliberate isolation: a user can safely explore a completely unrelated file (a candidate new dataset, a sanity check on someone else's export) without any risk of it interacting with whatever is currently staged/active in the rest of the app. Statistically, this tab performs genuine inferential diagnostics (Shapiro-Wilk normality testing, skewness/kurtosis, robust Iglewicz-Hoaglin outlier z-scores) that no other tab in the module runs — it is deliberately more statistically rigorous about *characterizing* the data than the pipeline tabs are, precisely because its job is diagnosis, not transformation.

**What happens if this step is skipped:** Nothing about the pipeline breaks if this tab is never used — it's optional, diagnostic-only. What is lost is an independent, pre-pipeline sanity check: e.g. discovering that an uploaded file is already normalized (making a second normalization pass in Batch Correction redundant or even harmful) or that several samples show multiple simultaneous outlier signals (worth investigating before they silently distort a merge/correction elsewhere).

### 5.B UI inventory

`mod_data_exploration_ui()` (L712–719) is a thin shell — one static info `div` plus `withSpinner(uiOutput(ns("body_ui")))`. Nearly everything else is generated dynamically by `output$body_ui` and, after a run, `output$results_ui`.

| Element | Input/Output ID | What the user sees | Controls | Reacts via |
|---|---|---|---|---|
| `fileInput(ns("raw_file"), ...)` | `...-raw_file` | Upload box, accepts .csv/.tsv/.txt only (no RDS, no Excel — stated explicitly) | The single dataset this whole tab operates on | `observeEvent(input$raw_file)` |
| *(after a valid parse)* `eda_upload_info_ui()`'s value boxes + `DT::dataTableOutput(ns("head_preview_table"))` | — | Immediate (no button) structural summary: row/column counts, detected ID column, non-numeric columns set aside, first-rows preview | Read-only | `output$head_preview_table` |
| `actionButton(ns("run_btn"), ...)` | `...-run_btn` | "Run Exploratory Data Analysis" | Triggers the entire `eda_result` pipeline | `eda_result` |
| `uiOutput(ns("results_ui"))` | — | 13 lettered result boxes (A–M), each with its own tabbed sub-outputs | Read-only | assembled once `eda_result()` succeeds |
| `downloadButton(ns("download_summary"), ...)` (inside box M) | — | CSV export of the final plain-language summary | — | `downloadHandler` |

### 5.C Function inventory — Data Exploration tab

*(Every function is listed and categorized in §1.6; this table adds the specific per-function role within the EDA pipeline's actual execution order.)*

| Function | Type | Purpose | Input | Output |
|---|---|---|---|---|
| `eda_parse_upload()` | Data-processing + validation | Parses any delimited upload into a feature×sample numeric matrix, auto-detecting the ID column and non-numeric columns; never throws | `datapath, filename` | `list(ok, expr, id_col_name, nonnumeric_cols, ...)` or `list(ok=FALSE, error=...)` |
| `eda_skewness()`, `eda_kurtosis()` | Statistical | Third/fourth standardized moments | numeric vector | numeric |
| `eda_robust_z()` | Statistical | Iglewicz-Hoaglin modified (median/MAD-based) z-score | numeric vector | numeric vector |
| `eda_skew_label()` | Utility | Numeric skewness → plain-language label | skewness value | string |
| `eda_overview()` | Statistical | Dataset-level summary (dimensions, missingness, duplicates, constant/near-zero-variance features, pooled descriptive stats) | parsed upload | `list(...)` |
| `eda_descriptive_stats()` | Statistical | Per-feature or per-sample descriptive statistics (n, mean, median, sd, var, quantiles, skewness, kurtosis, CV) | matrix, `margin` | data.frame |
| `eda_normality_summary()` | Statistical | Pooled skewness/kurtosis + capped Shapiro-Wilk test | matrix | `list(...)` |
| `eda_normalization_assessment()` | Statistical + validation | Evidence-based verdict: normalized / not normalized / inconclusive, reusing `detect_expr_data_type()`/`summarize_norm_diagnostics()`/`needs_quantile_norm()` | matrix | `list(verdict, label, evidence, ...)` |
| `eda_impute_median()` | Data-processing | Diagnostic-only median imputation, matrix-wide fallback for all-missing rows | matrix | matrix (never written back) |
| `eda_prep_for_structure()` | Data-processing | Impute + restrict to top-variance features, for PCA/correlation | matrix | `list(sub, n_features_used, n_features_total)` or `NULL` |
| `eda_pca()` | Statistical | PCA on the prepared subset | matrix | `list(scores, var_exp, n_features_used/total)` or `NULL` |
| `eda_sample_correlation()` | Statistical | Pairwise sample correlation on the prepared subset | matrix | `list(cor, n_features_used/total)` or `NULL` |
| `eda_hclust()` | Statistical | Average-linkage hierarchical clustering on `1 - correlation` distance | correlation matrix | `hclust` object |
| `eda_sample_outliers()` | Statistical | Combines `compute_sample_qc()` (global.R) with a PCA-distance-from-centroid robust-z flag | matrix, `pca` | data.frame with `n_flags` |
| `eda_feature_outliers()` | Statistical | Extreme-variance / extreme-skew / high-missingness feature flags | matrix, descriptive-stats df | data.frame |
| `eda_missingness()` | Statistical | Per-sample missing counts + per-feature missingness bucket histogram | matrix | `list(by_sample, by_feature_summary)` |
| `eda_mean_variance_df()` | Statistical | Per-feature mean vs. variance | matrix | data.frame |
| `eda_transform_diagnostic()` | Statistical | Raw-vs-log2 skewness comparison, diagnostic only | matrix | `list(can_log, raw_sample, log_sample, skew_raw, skew_log)` |
| `eda_final_summary()` | Statistical + validation | Synthesizes every section above into one quality/next-steps verdict | overview, norm_assess, normality, samp_outliers, feat_outliers | `list(quality, normalization, distribution, outliers, missing, variance, next_steps)` |
| 14 plot builders (`eda_hist_plot()` … `eda_scree_plot()`) | Plotting | Pure `data → ggplot/plotly` functions, one per visualization | varies | ggplot/plotly object |
| `eda_section_card()`, `eda_status_panel_ui()`, `eda_summary_card_ui()`, `eda_upload_info_ui()` | UI composition | Consistent card/panel markup | varies | Shiny tags |
| `mod_data_exploration_ui()`, `mod_data_exploration_server()` | UI / Server | The Shiny module itself | `id` | tags / side effects |

### 5.D Line-by-line teaching — the most scientifically load-bearing functions

#### `eda_parse_upload()` (L39–94) — auto-detecting structure from an arbitrary file

```r
df <- tryCatch(as.data.frame(data.table::fread(datapath, showProgress = FALSE,
                 na.strings = c("NA", "", "NaN", "null", "NULL", "#N/A"))), error = function(e) NULL)
if (is.null(df)) return(list(ok = FALSE, error = "..."))
if (nrow(df) == 0 || ncol(df) < 2) return(list(ok = FALSE, error = "..."))
```
**What/why:** unlike every other upload path in this module (which all assume a fixed "first column = ID, rest = numeric" shape), this parser must handle a genuinely unknown file, so `na.strings` is set defensively wide — recognizing not just R's own `"NA"` but also empty strings, `"NaN"`, and the literal text `"null"`/`"NULL"`/`"#N/A"` (an Excel-export convention) as missing-value markers a real-world CSV export might use. Two early `validate`-free (plain `if`/`return`) guards catch a totally unparseable file or a degenerate (too few rows/columns) one — this function **never throws**, always returning a `list(ok=, ...)` sentinel object instead, the same pattern taught in `mod_dataset_teaching_notes.md` for `tryCatch`-based error handling, but pushed one level further: even a parse failure is captured and returned as data, not allowed to propagate as an R condition at all.

```r
first_col <- df[[1]]
first_num <- suppressWarnings(as.numeric(as.character(first_col)))
id_is_char <- mean(is.na(first_num)) > 0.5
if (id_is_char) {
  ids <- as.character(first_col); rest <- df[, -1, drop = FALSE]; id_col_name <- colnames(df)[1]
} else {
  ids <- paste0("row_", seq_len(nrow(df))); rest <- df; id_col_name <- NULL
}
```
**What:** the **ID-column detection heuristic** — coerce the first column to numeric and check what fraction fails (`is.na`). If **more than half** fail to parse as numeric, the first column is treated as a genuine identifier column (gene/probe/CpG names) and split off; otherwise, every column (including the first) is treated as numeric data, and synthetic `row_1, row_2, ...` IDs are generated instead. **Why >50%, not "any failure":** a real ID column will be *entirely* non-numeric in the overwhelming majority of cases (gene symbols, probe IDs), so a simple any-failure rule would be equivalent in practice — the `>0.5` threshold is a deliberately generous majority-vote rule that tolerates a handful of numeric-looking IDs (e.g. Entrez gene IDs, which are literally numbers) still being correctly classified as an ID column, as long as most of the column isn't numeric.

```r
is_num_col <- vapply(rest, function(col) {
  if (is.numeric(col)) return(TRUE)
  v <- suppressWarnings(as.numeric(as.character(col)))
  blank <- is.na(col) | (is.character(col) & trimws(as.character(col)) == "")
  mean(is.na(v) & !blank) < 0.2
}, logical(1))
```
**What:** for every remaining column, decide whether it counts as a numeric "sample" column. A column that's already R-typed numeric passes trivially; otherwise, attempt coercion and check the failure rate **excluding legitimate blanks** (`!blank` — a truly empty cell shouldn't count against a column's "is this numeric" case, since blanks are handled separately as missing values, not as evidence of non-numeric-ness). **Why `< 0.2`, a 20% failure tolerance, rather than requiring 100% clean parsing:** real-world files often have a stray non-numeric artifact (a footnote marker, a formatting error) in an otherwise-numeric column; this threshold tolerates that without discarding the whole column, while still excluding genuinely categorical/text columns (which would fail far more than 20% of the time). This is a materially more permissive numeric-detection rule than anywhere else in the module — worth citing directly if the thesis discusses this tab's specific engineering tradeoffs versus the stricter, single-purpose parsers elsewhere (e.g. `tx_parse_expr_matrix_rds()`, which requires *all* columns after the first to be numeric with no tolerance).

#### `eda_normalization_assessment()` (L232–273) — turning raw diagnostics into a plain-language, evidenced verdict

```r
dt <- detect_expr_data_type(m)
diag <- summarize_norm_diagnostics(m)
differs <- needs_quantile_norm(diag)
```
**What:** this function is explicitly built **on top of** the same three already-audited `global.R` primitives Batch Correction itself uses (§4.F) — the code's own comment (L225–231) states this directly: "Built on top of the app's own already-audited detection primitives... rather than re-deriving the same heuristics a second time." This is a genuinely good design decision worth citing: the EDA tab's normalization verdict cannot disagree with what Batch Correction's own auto-detect logic would decide, because it's calling the identical functions, not a parallel reimplementation.

```r
if (identical(dt, "counts")) {
  verdict <- "not_normalized"; label <- "Likely not normalized (raw counts)"
  evidence <- c(evidence, sprintf("%.0f%% of finite values are at or near integers and non-negative, with a maximum value of %s...", frac_integer * 100, ...), "Raw counts of this kind have not yet been adjusted for library size or composition...")
} else if (has_negative) {
  verdict <- "normalized"; label <- "Likely normalized / transformed"
  evidence <- c(evidence, "Negative values are present, which is only possible after a log-ratio, z-score, or similarly centered transform...")
  ...
} else if (!differs) {
  verdict <- "normalized"; label <- "Likely normalized"
  ...
} else {
  verdict <- "inconclusive"; label <- "Inconclusive"
  ...
}
```
**What — the four-branch verdict logic, each branch producing both a verdict *and* a written evidence trail:** `detect_expr_data_type()` (`global.R` L1616–1635, already documented in Batch Correction's own filter table, §4.E) does the heavy lifting of the "counts vs. already-normalized vs. plain expression" three-way classification; this function adds two more evidence checks specifically for the "normalized" verdict (negative values present — only possible post-transform; or `!differs`, i.e. per-sample medians/IQRs already agree) and a genuinely honest fourth outcome, **"inconclusive,"** for when the evidence is mixed (continuous-looking values, but per-sample distributions still disagree) — this is the branch worth highlighting: the function is explicitly designed to say "I don't know" rather than force a binary answer when the signal doesn't support one. **Why this matters for scientific honesty in a thesis:** an automated verdict tool that always returns a confident yes/no despite genuinely ambiguous evidence would be a worse tool than one that discloses its own uncertainty — this is a direct, code-verifiable example of that principle being implemented, not just claimed.

#### `eda_sample_outliers()`/`eda_feature_outliers()` (L326–363) — combining multiple independent QC signals

```r
qc <- tryCatch(compute_sample_qc(eda_impute_median(m)), error = function(e) NULL)
...
qc$pca_distance <- NA_real_; qc$flag_pca <- FALSE
if (!is.null(pca) && ncol(pca$scores) >= 2) {
  scores <- pca$scores[, 1:2, drop = FALSE]
  center <- colMeans(scores)
  dist <- sqrt(rowSums((scores - matrix(center, nrow(scores), ncol(scores), byrow = TRUE))^2))
  z <- eda_robust_z(dist)
  idx <- match(qc$sample, rownames(scores))
  qc$pca_distance <- dist[idx]
  qc$flag_pca <- ifelse(is.na(idx), FALSE, abs(z[idx]) > 3.5)
}
qc$n_flags <- rowSums(qc[, c("flag_signal", "flag_detected", "flag_cor", "flag_pca")], na.rm = TRUE)
```
**What — a fourth, genuinely new outlier signal added on top of `compute_sample_qc()`'s existing three (§4.F):** Euclidean distance from each sample's PC1/PC2 coordinates to the cohort centroid (`colMeans`), converted to a **robust z-score** (`eda_robust_z()` — the same Iglewicz-Hoaglin median/MAD formula, `0.6745 * (x - med) / mad`, the standard scaling constant that makes this comparable to a normal-theory z-score under normality) and flagged beyond `|z| > 3.5` — a slightly more conservative threshold than Batch Correction's user-adjustable `mad_k` default of 3, and, unlike that setting, **not exposed to the user at all here** (a hard-coded constant). **Why combine four independent signals into one `n_flags` count rather than picking the single "best" one:** each metric can catch a different failure mode a sample might exhibit (unusual total signal, unusual detection rate, unusual dissimilarity to the cohort, unusual position in reduced-dimension space) — a sample flagged by multiple independent signals simultaneously is a much stronger candidate for genuine investigation than one flagged by only one, which is exactly how `eda_final_summary()` (below) uses this count (`n_flags >= 2` as its own "worth investigating" bar).

`eda_feature_outliers()` mirrors this at the feature level with three signals — extreme (robust-z) variance, extreme skewness (`|skewness| > 2`, a fixed, uncited threshold), and high missingness (`> 20%`, likewise fixed and uncited) — each independently flaggable, none combined into a single feature-level score the way samples are (features are ranked in the results table by raw `n_flags` sum, computed inline in `output$feat_outlier_table` rather than stored on the data.frame itself).

### 5.E Filters — Data Exploration tab

This tab **flags, but never filters or removes**, anything — this is a deliberate, stated design choice (the module intro text: "nothing you upload here is written back anywhere"; the outlier section's own description: "Nothing is removed automatically - flagged rows/columns are for investigation, not automatic exclusion"). The "filters" present are therefore all **diagnostic thresholds that gate a flag, not a removal**:

| Threshold | Value | Adjustable? | What it flags |
|---|---|---|---|
| Sample-level `mad_k` (via `compute_sample_qc()`) | 3 (function default, not exposed here) | No — hard-coded default, not surfaced in this tab's UI (unlike Batch Correction's own `mad_k` slider) | `flag_signal`, `flag_detected`, `flag_cor` |
| PCA-distance robust z | 3.5 | No | `flag_pca` |
| Feature variance robust z | 3.5 | No | `flag_extreme_variance` |
| Feature skewness | `\|skew\| > 2` | No | `flag_extreme_skew` |
| Feature missingness | `> 20%` | No | `flag_high_missing` |
| "Investigate" bar (`eda_final_summary`) | `n_flags >= 2` (samples) | No | the plain-language "outliers" summary line |
| Data-quality verdict thresholds (`eda_final_summary`) | `pct_missing > 20`/`5`; duplicated features `> 5%` of total | No | "Needs attention" / "Moderate" / "Good" quality label |

None of these are re-derivable from the UI — they are only visible by reading the source, which is exactly why they are enumerated here for thesis-methods citation purposes; if this tab's outlier/quality verdicts are quoted in a thesis, these are the exact numeric thresholds behind them.

### 5.F Results — what the user should see

Thirteen lettered result boxes (A–M), each following the brief's own "result → visualization → interpretation" shape (`eda_section_card()`'s standard structure):

- **A. Dataset overview:** 8 value boxes (samples, features, % missing, infinite values, duplicated feature IDs, duplicated sample columns, constant features, near-zero-variance features) + a 9-statistic pooled table (mean/median/sd/var/min/max/Q1/Q3/IQR). Correct behavior: `n_missing`'s displayed percentage should match `pct_missing` in Filter 5's `eda_final_summary()` verdict, since both come from the same `eda_overview()` call.
- **B. Descriptive statistics:** two tabs, feature-level and sample-level, each a full `eda_descriptive_stats()` table (14 columns including skewness/kurtosis/CV).
- **C. Distribution analysis:** 5 tabs — histogram, density, per-sample density overlay, boxplot, violin — all on the pooled or per-sample value distribution, sampled down via the `EDA_MAX_*` constants (L23–28) when the matrix is large, to keep rendering responsive.
- **D. Normality/distribution assessment:** Q-Q plot + skewness/kurtosis/Shapiro-Wilk numbers, with an explicit caveat that Shapiro-Wilk "becomes extremely sensitive at large n... a significant p-value here is common for real data and is not, by itself, evidence of a meaningful problem" — a genuinely important, code-stated methodological caution (high-dimensional molecular data is essentially never exactly Gaussian at the whole-matrix level, and this box says so directly rather than implying otherwise).
- **E. Normalization status assessment:** the four-verdict panel from §5.D, with its full evidence list and four fact chips (detected type, % near-integer, negative values present, between-sample agreement).
- **F. Outlier detection:** an interactive (plotly) boxplot colored by outlier status, plus full sample- and feature-level flag tables, sorted worst-first.
- **G. PCA/sample structure:** interactive PCA scatter (colored by outlier status) + scree plot, computed on the top-variance-feature subset (count disclosed in the panel text).
- **H. Sample correlation/distance:** a hierarchically-ordered correlation heatmap + dendrogram, same top-variance subset.
- **I. Missing data analysis:** a per-sample missingness bar chart + a per-feature missingness bucket table (0% / >0–5% / >5–20% / >20–50% / >50%).
- **J. Low-variance features:** a plain-text count of constant + near-zero-variance features (no plot — a deliberately minimal section).
- **K. Mean-variance relationship:** a scatter of per-feature mean vs. variance with a LOESS smooth — the panel's own description states the expected pattern for count-like data (variance increasing with mean) as the signal that a variance-stabilizing transform may be worth considering.
- **L. Before/after transformation diagnostic:** side-by-side raw vs. log2 density plots (diagnostic only, explicitly labeled as not modifying the uploaded data), with a skewness-comparison sentence.
- **M. EDA summary:** the final synthesized verdict card (`eda_final_summary()`'s output) plus a CSV download of the same.

**Bug/regression indicators specific to this tab:** any of the 13 boxes rendering with no content at all (rather than a `validate()`-styled guidance message) after a successful `run_btn` click would indicate a broken reactive dependency (`eda_result()` failing silently); a normalization verdict that contradicts Batch Correction's own auto-detect decision on the *same* uploaded file would indicate the shared-primitive design (§5.D) has been broken by a future edit to one path but not the other.

### 5.G Validation checklist — Data Exploration tab

**Input validation:** upload a non-delimited or binary file → confirm `eda_parse_upload()`'s `tryCatch`-wrapped `fread()` returns the clean parse-failure message, not a raw R error; upload a file with zero data rows → confirm the explicit row/column-count guard fires; upload a file with a numeric-looking first column (e.g. numeric Entrez IDs) → confirm the `>50%`-non-numeric heuristic still classifies it correctly (or, if it doesn't, confirm the fallback `row_N` synthetic-ID behavior is at least not silently wrong — it would misassign the true ID column as a data column instead).

**Functional validation:** confirm `eda_result()` populates all of `overview`, `feat_stats`, `samp_stats`, `normality`, `norm_assess`, `pca`, `corr`, `feat_outliers`, `samp_outliers`, `missingness`, `meanvar_df`, `transform_diag`, `summary` after a successful run; confirm `pca`/`corr` gracefully return `NULL` (triggering the "not enough data" guidance message, not a crash) on a deliberately tiny (< 3 sample) test file.

**Scientific validation:** upload a dataset with a known, deliberately-introduced batch effect (e.g. the two raw training GEO series concatenated *without* the app's own gene-symbol collapse, so probe-level features barely overlap) and confirm the PCA/correlation sections visually surface the expected structure; upload a matrix already known to be quantile-normalized and confirm the normalization-assessment verdict correctly reports "Likely normalized," cross-checked against the same file's verdict if run instead through Batch Correction's own auto-detect logic.

**Reproducibility validation:** re-running `eda_result()` on the identical uploaded file should be deterministic **except** for the capped-sampling steps (`EDA_MAX_POOLED_VALUES`, `EDA_MAX_SHAPIRO_N`, `EDA_MAX_VIOLIN_FEATURES`, `EDA_MAX_DENSITY_SAMPLES`, all of which call `sample()` without a fixed seed) — **this is a real, disclosable non-determinism**: on a large matrix, the histogram/density/Shapiro-Wilk/violin/mean-variance-scatter visuals and statistics are computed on a *different random subsample* each run, so two runs on the same file can show slightly different Shapiro-Wilk p-values or plot shapes purely from sampling variability, not from any change in the underlying data. This should be disclosed in a thesis Methods/Limitations section wherever this tab's diagnostics are cited for a large (> `EDA_MAX_POOLED_VALUES` = 200,000-value) matrix.

### 5.H Thesis-ready interpretation — Data Exploration tab

**Methodological description:** An independent, non-destructive exploratory data analysis workflow applied to any user-uploaded feature-by-sample matrix: automatic identifier-column detection, descriptive statistics at both feature and sample level, distributional and normality diagnostics (skewness, kurtosis, capped Shapiro-Wilk), an evidenced normalization-status verdict reusing the same detection logic as the app's Batch Correction pipeline, unsupervised structure analysis (PCA and hierarchical sample clustering on the top-variance feature subset), multi-signal sample- and feature-level outlier flagging, missing-data and low-variance characterization, a mean-variance relationship check, and a synthesized plain-language quality summary with a recommended next step.

**Computational implementation:** A fully separate nested Shiny module (`mod_data_exploration_ui`/`mod_data_exploration_server`, `mod_preprocessing_explore.R`) with its own upload and its own single `eventReactive` pipeline (`eda_result`); shares no reactive state with the rest of the Preprocessing/Merge/Batch-correction pipeline and never mutates the uploaded data.

**Parameters:** All outlier/quality thresholds (§5.E) are fixed constants, not user-adjustable from this tab's UI; performance caps on sampled computations (`EDA_MAX_*` constants, L23–28) bound how much of a very large matrix is actually used for the pooled/violin/mean-variance/Shapiro-Wilk diagnostics.

**Validation:** See §5.G; the tab's normalization verdict is directly cross-checkable against Batch Correction's own auto-detect decision on the same file, since both call the identical `global.R` primitives.

**Expected results:** A 13-section diagnostic report culminating in one synthesized quality/next-step verdict, exportable as a summary CSV.

**Limitations:** Every outlier/quality threshold in this tab is a fixed, uncited constant, not user-adjustable or independently derived for the specific dataset being explored. Several diagnostics (Shapiro-Wilk, histogram/density/violin sampling) are computed on a random, unseeded subsample when the matrix exceeds fixed size caps, making exact reported values non-reproducible run-to-run on large datasets, though the qualitative conclusions should remain stable. The ID-column and numeric-column detection heuristics are majority-vote rules (>50% / <20% failure tolerance respectively) that can, in principle, misclassify an unusual but valid file layout.

---

## PART 6 — CODE AUDIT: BUGS, RISKS, AND DESIGN ISSUES

*(This section consolidates every finding cross-referenced above as "§7" — that citation style is kept in Parts 2–5 for internal consistency, but resolves to this Part 6. Each entry is labeled per the brief's own required taxonomy: **BUG** / **POTENTIAL BUG** / **SCIENTIFIC RISK** / **DESIGN ISSUE** / **VALIDATED CORRECT**.)*

| # | Location | Label | Finding |
|---|---|---|---|
| 1 | `mod_pp_source_ui()`/`mod_pp_source_server()`, whole functions (L205–603) vs. `output$preprocessing_tab_ui` (L791–809) | **DESIGN ISSUE** | `mod_pp_source_ui()` — the rich per-source panel with its own upload/sample-filter/feature-filter/log2/run UI — is never inserted into `mod_preprocessing_ui()`'s tab 1 markup. Six `mod_pp_source_server()` instances are still instantiated (`pp_sources`, L735–738) and their reactive machinery genuinely runs, but with no UI for a user to reach it, this is server logic with no corresponding UI control, per the audit brief's own checklist item. Practical consequence: the module's only reachable per-source upload+filter workflow is the simpler checkbox-driven `pp_preloaded_read()` path, which has no equivalent of the rich panel's categorical/numeric sample filters, dedup, or per-source log2 override for an *uploaded* dataset specifically (that richer control set only exists in the unreachable panel). |
| 2 | `dataset$staged_expr %||% dataset$expr` fallback rule, duplicated in `pp_preloaded_read()` (L145–147) and `mod_pp_source_server`'s `current_source()` (L303–306) | **DESIGN ISSUE** | The same 3-line "prefer staged, fall back to active" rule is written out twice rather than shared via one helper. A future change to this rule (e.g. adding a third fallback tier) would need to be applied in both places by hand, with nothing in the code enforcing that. |
| 3 | `pp_guess_col()` (module scope), the local `guess_col()` inside `mod_pp_source_server` (L363–369), and `mod_dataset.R`'s own `guess_col()` | **DESIGN ISSUE** | Three byte-identical ~6-line implementations of the same column-name-guessing heuristic exist in the codebase. The comment at L47–49 documents that the module-scope copy was deliberately pulled out for the checkbox-upload path to reuse, but the original in-module copy was left in place rather than replaced with a call to the new shared one — the drift risk is real (a fix to one copy's edge-case handling would not automatically propagate). |
| 4 | `PP_COHORT_LABELS` (this file, L109–114) vs. `INDIVIDUAL_DATASET_LABELS` (`mod_dataset.R`) | **DESIGN ISSUE** | Two independent GEO-accession → display-label tables must be kept in sync by hand whenever a new bundled source is added to `GEO_SOURCES`. The fallback (raw accession ID) is graceful, not a crash, but a silently stale label in one of the two tables is a real, easy-to-miss possibility. |
| 5 | `pp_preloaded_read()`'s "Merged Data" branch (L155–160) | **SCIENTIFIC RISK** | Selecting "Merged Data" from what is presented as a "raw, single-platform data - not yet merged or normalised" cohort picker (the UI's own text at L231–232) actually reads the bundled, already-merged-and-batch-corrected cohort directly via `load_default_dataset()`. A user relying on the panel's own stated description without checking which specific option they picked could unknowingly carry an already-corrected dataset through this tab's per-source filters and into a "fresh" merge/batch-correction run, believing it to be raw. |
| 6 | `example_merge_from_raw` (L1012), hardcoded `reactive(FALSE)` | **DESIGN ISSUE** | "Merge the example pipeline's training datasets" presents itself in the UI as reconstructing the training cohort live from its two raw GEO sources, but as shipped (this switch fixed at `FALSE`) it instead reuses the bundled, already-batch-corrected cohort directly, skipping the from-raw rebuild. This is disclosed in the code's own extensive comment and in `merge_venn_example_ui`'s rendered text when this path is taken, but a user who doesn't read that disclosure and later runs Batch Correction on this merge would see little or no batch effect to remove and could mistakenly conclude the correction step is broken or unnecessary, rather than that the input was already corrected. |
| 7 | `example_live_merge()`'s TRUE-branch (L1041–1078) vs. `pp_preloaded_read()` (L139–203) | **DESIGN ISSUE** | The from-raw rebuild path (currently unreachable via the UI per finding #6, but present and maintained in the code) re-implements the same read→collapse→filter→log2→impute sequence as `pp_preloaded_read()`, inline, rather than calling it. A second, independently-maintained copy of the same pipeline logic. |
| 8 | 20-common-feature floor, `merged`/`example_live_merge` (both branches) | **DESIGN ISSUE** | A hard-coded, non-adjustable minimum with no stated statistical derivation in the code or its comments. A merge passing this floor with, e.g., 22 shared genes is not flagged as potentially under-powered for any specific downstream use. |
| 9 | `edgeR::filterByExpr()`'s internal CPM/sample-count thresholds (used inside the `result` pipeline's TMM branch) | **DESIGN ISSUE (hidden threshold)** | This is a real, data-adaptive filter, but its internal parameters are not exposed anywhere in this tab's UI — a user can observe the resulting gene count (`vb_genes_dropped`) but cannot see or adjust what drove it, unlike every other filter on this tab, which has a visible slider. |
| 10 | ComBat's 3-deep `tryCatch` fallback ladder (`run_combat`, L1635–1687) | **DESIGN ISSUE** (disclosed, not silent) | If the fully-specified ComBat call (with `mod` and `ref.batch`) fails, the code silently retries with progressively less specification, ultimately without `mod` at all — meaning the correction that actually ran may not have protected the user's chosen covariates. This is disclosed via `decisions_ui`'s "Protected: ..." line reflecting what actually ran, not what was requested, but nothing forces a user to notice a fallback occurred versus reading only the correction-method headline. |
| 11 | ComBat-seq path overriding other Settings inputs (L1398, code at L1494–1543) | **DESIGN ISSUE** | Selecting "Before TMM: ComBat-seq" silently ignores the correction-method dropdown, ComBat prior, reference batch, and exclude-outliers settings — stated explicitly in the UI's own info text, but a real coupling between one radio choice and several otherwise-independent-looking controls. |
| 12 | `sva::sva()`/`sva::num.sv()` calls (`run_sva()`, L1654–1673) | **SCIENTIFIC RISK / reproducibility** | No `set.seed()` call anywhere in this file. SVA's surrogate-variable estimation is not guaranteed bit-for-bit deterministic run-to-run; re-running the identical SVA configuration on identical data may produce closely comparable but not identical corrected values. |
| 13 | `EDA_MAX_POOLED_VALUES`/`EDA_MAX_SHAPIRO_N`/`EDA_MAX_VIOLIN_FEATURES`/`EDA_MAX_DENSITY_SAMPLES` sampling caps (`mod_preprocessing_explore.R`, throughout) | **SCIENTIFIC RISK / reproducibility** | Every capped diagnostic (`sample()` calls with no fixed seed) is computed on a different random subsample each run on a large matrix — reported Shapiro-Wilk statistics and several plots are not bit-for-bit reproducible run-to-run on data exceeding these caps, though the qualitative interpretation should typically remain stable. |
| 14 | Matrix-orientation validation, everywhere in the module | **POTENTIAL BUG (absence of a check)** | No function in this module checks whether an uploaded expression matrix is genuinely oriented features-in-rows/samples-in-columns versus the transpose. A transposed upload would likely still "work" mechanically (produce a matrix, pass through filters) while being scientifically nonsensical (treating samples as if they were genes and vice versa), and nothing in the code would flag this. |
| 15 | Duplicate-feature handling on merge (`x$expr[common, , drop = FALSE]`, both merge paths) | **DESIGN ISSUE (disclosed)** | R's row-name-keyed subsetting silently keeps only the first matching row for a duplicated feature ID. This is *counted* and disclosed to the user (`n_dup_features`, surfaced as a warning in `merge_summary_ui`), which is the important distinguishing fact — this is not a silent failure, but the specific row kept ("first" by original file order) is still an arbitrary tie-break with no scientific justification for why that particular duplicate, rather than another, survives. |
| 16 | Batch-vs-protect-column collinearity guard (`protect_dropped_for_batch`, `result` pipeline, L1613–1621) | **VALIDATED CORRECT** | A genuinely important, correctly-implemented safeguard: a column cannot simultaneously be the batch variable being corrected for and a covariate protected from correction, preventing a degenerate/collinear design matrix. Disclosed to the user via `decisions_ui`'s own warning line when it triggers. |
| 17 | `outputOptions(output, ..., suspendWhenHidden = FALSE)` loop for Batch Correction's outputs (L2093–2095) | **VALIDATED CORRECT** | A deliberate, documented fix for a real Shiny footgun: a `uiOutput` nested inside a not-yet-visited `tabPanel` defaults to a suspended state and can miss its first render trigger, requiring the user to leave and revisit the tab to see content that should already be there. Forcing `suspendWhenHidden = FALSE` on this tab's own outputs (and mirroring the same fix already applied in `mod_wgcna.R`, per the comment) is the correct, minimal countermeasure. |
| 18 | Data Exploration tab's isolation from the shared `dataset` reactiveValues (`mod_data_exploration_server(id)` takes only `id`, no `dataset` argument) | **VALIDATED CORRECT** | Confirmed by reading the function signature and every line inside it: no read or write of `dataset` occurs anywhere in `mod_preprocessing_explore.R`. The module's own claim of full independence from the rest of the app's state is accurate as implemented. |
| 19 | `pp_preloaded_read()`/`example_live_merge()`'s per-source auto-log2 heuristic (`q99 > 100` on positive values) | **SCIENTIFIC RISK (disclosed in-code)** | The code's own comment (L504–510) already discloses this heuristic's known failure mode (cannot distinguish raw RNA-seq counts from already-normalised large-valued data by magnitude alone) — included here for completeness since it is a genuine, code-acknowledged limitation, not a silent one, and is exactly the kind of threshold a thesis reviewer would reasonably probe. |

---

## PART 7 — FINAL LEARNING SUMMARY

### 7.A Complete function map

```
mod_preprocessing_config (registry entry)
│
├── mod_preprocessing_ui(id) ── UI ── 4× tabPanel (Preprocessing / Merge datasets / Batch correction / Explore)
│
└── mod_preprocessing_server(id, dataset, results)
    │
    ├── [Tab 1: Preprocessing]
    │     pp_cohort_choices() ← pp_cohort_label() ← PP_COHORT_LABELS, default_dataset_entry (mod_dataset.R)
    │     pp_preloaded_read() ← load_default_dataset() / load_individual_dataset() / get_raw_eset()
    │                            + get_collapsed_genes() + eset_harmonize_meta()      (global.R)
    │     observeEvent(preloaded_run) → preloaded_results_val → preloaded_results()
    │     [unreachable via UI, but instantiated:] pp_sources[1..6] = mod_pp_source_server()
    │          each: raw_pair() → group/numeric/dedup filters → log2_ui → result (eventReactive)
    │
    ├── [Tab 2: Merge datasets]
    │     merge_inputs() ← preloaded_results()
    │     mode "example": available_example_groups(), example_merge_from_raw(), example_live_merge()
    │          ← get_raw_eset()/get_collapsed_genes()/eset_harmonize_meta() (2 training GSEs)
    │     mode "own": collapse_annot(), selected_lst() ← pp_collapse_probes_to_genes()
    │          overlap_sets() → draw_overlap_venn()/overlap_region_sizes()   (global.R)
    │     merged (eventReactive) → dispatches to example_live_merge() or own-data intersect+cbind+rbind
    │
    ├── [Tab 3: Batch correction]
    │     active_meta_df() ← merged()
    │     result (eventReactive) ─┬─ TMM branch: edgeR::filterByExpr/calcNormFactors/cpm
    │                             │     → optional sva::ComBat_seq (pre-TMM)
    │                             ├─ quantile/skip/auto branch: mean/var percentile filter
    │                             │     → summarize_norm_diagnostics()/needs_quantile_norm() (global.R)
    │                             │     → limma::normalizeBetweenArrays (if triggered)
    │                             ├─ optional pre-correction outlier exclusion: compute_sample_qc() (global.R)
    │                             ├─ mod ← model.matrix() from protect_cols, batch-collinearity guard
    │                             └─ correction: sva::ComBat | limma::removeBatchEffect | sva::sva+num.sv
    │     pca_of()/plot_pca_advanced()/scree_plot() (global.R) → before/after PCA + scree
    │     assoc_pvalue() → PC1-vs-batch ANOVA, before/after
    │     activate_btn → dataset$expr/meta/source (the ONLY write to shared state in this module)
    │
    └── [Tab 4: Explore] — mod_data_exploration_server("eda")  (fully isolated, own file)
          eda_parse_upload() → eda_result (eventReactive):
             eda_overview, eda_descriptive_stats, eda_normality_summary,
             eda_normalization_assessment ← detect_expr_data_type()/summarize_norm_diagnostics()/
                                              needs_quantile_norm() (global.R, SAME calls as Tab 3)
             eda_pca/eda_sample_correlation/eda_hclust,
             eda_sample_outliers ← compute_sample_qc() (global.R, SAME call as Tab 3's outlier exclusion),
             eda_feature_outliers, eda_missingness, eda_mean_variance_df, eda_transform_diagnostic
             → eda_final_summary() → 13-section results_ui + downloadable CSV
```

**The two threads that tie the whole module together, worth stating explicitly:** (1) the auto-detect log2/normalization heuristic (`q99 > 100` / `summarize_norm_diagnostics()`/`needs_quantile_norm()`) recurs in `pp_preloaded_read()`, `mod_pp_source_server`'s `result`, `example_live_merge()`, Batch Correction's `result`, and (via the shared `global.R` functions) the EDA tab's `eda_normalization_assessment()` — one heuristic family, applied consistently everywhere in the module rather than reinvented per tab; (2) `compute_sample_qc()` (`global.R`) is the single shared implementation of robust sample-level QC used identically by Batch Correction's optional pre-correction outlier exclusion and the EDA tab's outlier flagging — a genuine example of this codebase reusing one audited implementation rather than duplicating it (in contrast to the column-guessing helper's actual triplication, Part 6 finding #3).

### 7.B Complete data-flow map

**Input → Validation → Filtering → Transformation → Normalization → Batch Correction → Output**, using only steps genuinely present in this module:

```
INPUT            Bundled GEO source (raw ExpressionSet / raw counts) via get_raw_eset()/
                 load_individual_dataset(), OR user upload (CSV/RDS) via mod_pp_source_server's
                 expr_raw_preview()/meta_raw(), OR "currently loaded" (Dataset tab's staged/active
                 dataset)
                        │
VALIDATION       File-shape checks (ncol>=2, parseable), column-mapping validate()s (own-data path:
                 >=3 matching sample IDs), per-source result()'s ncol(expr)>=3 sample floor
                        │
FILTERING        [Tab 1, own-data path only, currently unreachable via UI per audit finding #1]
                 categorical/numeric sample filters → dedup → feature exclusion pattern →
                 missing-data-tolerance feature filter (median-impute residual gaps)
                 [Tab 1, checkbox path, what's actually used]: implicit !is.na(group) sample filter
                 (microarray branch only)
                        │
TRANSFORMATION   Auto-detect/force/skip log2 (q99>100 heuristic on positive values), per source,
                 pre-merge
                        │
                 [Tab 2] Probe-to-gene collapse (optional, own-data path: median/maxmean/mean) →
                 feature-ID intersection across selected sources (>=20 floor) → cbind/rbind merge
                        │
NORMALIZATION    [Tab 3] Expression/variance percentile gene filter (continuous data) OR
                 edgeR::filterByExpr (count data) → limma::normalizeBetweenArrays quantile
                 normalization (continuous, auto-detected via summarize_norm_diagnostics()/
                 needs_quantile_norm()) OR edgeR TMM + log2-CPM (count data)
                        │
BATCH CORRECTION [Tab 3] optional pre-correction outlier exclusion (compute_sample_qc(), MAD-based)
                 → protected-covariate design matrix (model.matrix(), batch-collinearity guard) →
                 sva::ComBat (default) | limma::removeBatchEffect | sva::sva (unknown sources) |
                 sva::ComBat_seq (raw counts, pre-TMM)
                        │
OUTPUT           Corrected expression matrix + metadata + before/after PCA/QC diagnostics →
                 (on explicit "Use this as the active dataset" click) dataset$expr/meta/source,
                 read by every other transcriptomics submodule
```

*(Tab 4/Data Exploration is not part of this pipeline at all — it is a parallel, dead-end diagnostic branch: `Input → [EDA diagnostics] → (nothing written anywhere)`.)*

### 7.C What I learned

**R programming concepts:** closures and lazy-evaluation factory functions (the `load()` field pattern shared with `mod_dataset.R`); `tapply()`/`vapply()`/`Reduce()` as R's group-by/type-safe-map/fold primitives; the `drop = FALSE` gotcha when subsetting matrices/data.frames down to one row or column; `%||%` as a hand-rolled null-coalescing operator; `match.arg()` for validating a fixed-choice string argument; `stopifnot()` vs. `validate(need(...))` as two genuinely different error-handling registers (internal-invariant crash vs. user-facing guidance); `setNames()`/named-vector indexing as R's lightweight hash-map idiom; `switch()` with an unnamed fall-through default.

**Shiny programming concepts:** the module (`moduleServer`/`NS`) namespacing pattern at scale (6 simultaneous instances of one module in `pp_sources`); `renderUI`/`uiOutput` for UI that depends on data not known until runtime (dynamically-generated dropdowns whose *choices* are a just-uploaded file's own column names); `eventReactive` as the mechanism for "only recompute on an explicit button click," used consistently at every pipeline stage (never auto-run); `eventReactive(list(a, b), ...)` responding to either of two independent triggers; `isolate()` to break a circular self-referential UI-rebuild dependency; `outputOptions(suspendWhenHidden = FALSE)` as the fix for a `uiOutput` inside a not-yet-visited tab missing its first render.

**Statistical/bioinformatics concepts:** the log2/linear-scale detection heuristic and why it can't distinguish raw counts from unlogged continuous data by magnitude alone; quantile normalization's assumption of comparably-shaped distributions vs. TMM's composition-bias correction for count data, and why the two are not interchangeable; edgeR's group-aware `filterByExpr()` vs. a flat expression-percentile cutoff; robust (median/MAD) vs. classical (mean/SD) outlier statistics, and the Iglewicz-Hoaglin modified z-score specifically; empirical-Bayes shrinkage (ComBat) vs. direct linear adjustment (limma) vs. latent-factor estimation (SVA) vs. negative-binomial count-scale correction (ComBat-seq) as four genuinely different statistical strategies for the same underlying batch-effect problem; why protecting a biological covariate from batch correction requires a design matrix, and why that covariate cannot also be the batch variable itself; feature-intersection cross-platform merging as the conservative alternative to imputation-based integration; PCA-based visual and ANOVA-based quantitative batch-effect assessment as complementary, not redundant, evidence.

**Reproducibility/software-engineering concepts:** the value (and cost) of duplicated helper logic across a codebase (three copies of the same column-guessing function vs. one shared, audited QC function reused verbatim across two tabs); the difference between a disclosed limitation (documented in a code comment and surfaced to the user) and a silent one; why unseeded random sampling inside a diagnostic tool is a genuine, citable reproducibility caveat even when it doesn't affect correctness.

### 7.D Most important functions to understand, ranked

**Beginner (start here — pure, small, illustrate one R/Shiny idea each):**
1. `pp_guess_col()` — case-insensitive exact→substring→fallback string matching.
2. `mod_pp_field_hint()` — a trivial UI-generating function.
3. `pp_cohort_label()`/`pp_cohort_choices()` — named-vector lookup with a fallback.
4. `eda_skewness()`/`eda_kurtosis()`/`eda_robust_z()` — small, self-contained statistical formulas.

**Intermediate (the module's own data-processing core):**
5. `pp_collapse_probes_to_genes()` — three genuinely different aggregation strategies for the same problem.
6. `pp_preloaded_read()` — the full per-source read→log2→impute pipeline, in one function.
7. `merge_inputs()`/`selected_lst()`/`overlap_sets()` — the Merge tab's data-flow spine.
8. `eda_parse_upload()`/`eda_normalization_assessment()` — heuristic-based structure/verdict inference from raw data alone.
9. `compute_sample_qc()` (global.R) — the shared robust QC implementation used by two different tabs.

**Advanced (the statistically and architecturally densest logic in the module):**
10. `merged`/`example_live_merge()` — the full intersect+cbind+rbind merge with its collinearity/consistency guards.
11. Batch Correction's `result` `eventReactive` (L1463–1715) — the single largest, most consequential function in the module: normalization branch selection, outlier exclusion, protected-covariate design matrix construction, and the four-way correction-method dispatch (ComBat/limma/SVA/ComBat-seq), including ComBat's fallback ladder.
12. `run_sva()` — surrogate-variable estimation and its degrees-of-freedom bounding logic.
13. `eda_sample_outliers()`/`eda_feature_outliers()` — combining multiple independent statistical signals into one flag count.

### 7.E Most important filters, at a glance

| Filter | Tab | Axis | Adjustable | Default |
|---|---|---|---|---|
| Categorical sample filter | Preprocessing (per-source panel) | Sample | Yes | none |
| Numeric range sample filter | Preprocessing (per-source panel) | Sample | Yes | none |
| Sample deduplication | Preprocessing (per-source panel) | Sample | Yes | off |
| Feature exclusion pattern (regex) | Preprocessing (per-source panel) | Feature | Yes | none |
| Missing-data tolerance | Preprocessing (per-source panel) | Feature | Yes, 0–80% | 0% (= complete cases) |
| Log2 auto/force/skip | Preprocessing (both paths) | Value scale | Yes | auto (preloaded/current) / skip (upload) |
| Implicit `!is.na(group)` | Preprocessing (checkbox/microarray path) | Sample | No | always on |
| Feature-ID intersection | Merge Datasets | Feature | No (which datasets, yes; the intersection itself, no) | — |
| Optional probe→gene collapse | Merge Datasets | Feature (aggregation) | Yes | off |
| Expression-percentile gene filter | Batch Correction | Feature | Yes, 0–90% | 0% |
| Variance-percentile gene filter | Batch Correction | Feature | Yes, 0–90% | 0% |
| Zero-variance exclusion | Batch Correction | Feature | No | always on |
| `edgeR::filterByExpr()` | Batch Correction (TMM branch) | Feature | No (internal thresholds hidden) | always on, ≥50-gene floor |
| Pre-correction outlier exclusion | Batch Correction | Sample | Yes (on/off + `mad_k`) | off |
| Batch/protect-column exclusivity | Batch Correction | Covariate set | No (automatic) | always on |
| (Diagnostic-only flags, never removals) | Data Exploration | Sample & feature | No (thresholds fixed) | always on |

### 7.F Practical validation checklist (for validating against real datasets)

1. **Round-trip a known dataset:** load a bundled cohort (e.g. GSE93272) through Preprocessing → Merge → Batch Correction, and independently recompute at least the gene/sample counts and one spot-checked expression value outside the app (a plain R script calling the same `global.R` functions directly) to confirm the app reproduces them.
2. **Confirm every filter's before/after counts are internally consistent:** `n_after <= n_before` always; `genes_kept + genes_dropped == n_before` in Batch Correction's value boxes.
3. **Confirm the log2 heuristic's decision is sane on data of a known scale:** feed genuinely raw RNA-seq counts through with "Auto" and confirm log2 is *not* applied (the documented failure mode, §2.E Filter 6); feed genuinely linear-scale microarray intensities and confirm it *is*.
4. **Confirm the merge floor and Venn/region table agree:** the diagram's "common to all" region size should exactly equal the merged matrix's row count.
5. **Confirm batch correction actually reduces the batch signal without erasing the biological signal:** PC1-vs-batch ANOVA p-value should increase after correction; a parallel PC-vs-group check (performed manually, since the app doesn't automate this) should show the group separation is not diminished.
6. **Deliberately trigger every `validate()` guard at least once** (batch column with 1 level, <20 common features, wrong file format, mismatched sample IDs) and confirm each produces its specific, correct plain-language message rather than a raw R error or a silent wrong result.
7. **Re-run an identical configuration twice** and confirm deterministic reproduction for every method except SVA and the EDA tab's capped-sampling diagnostics (Part 6, findings #12–13), where run-to-run variation is expected and should be disclosed, not treated as a bug.
8. **Cross-check the EDA tab's normalization verdict against Batch Correction's own auto-detect decision** on the identical uploaded file — they call the same underlying functions and should never disagree.

### 7.G Thesis documentation checklist

For every use of this module's output cited in a thesis, document:
- **Which specific path was taken** at each of the three ambiguous forks: (a) Preprocessing's "Merged Data" option vs. a genuine raw single-platform load (Part 6 finding #5); (b) Merge's "example" mode's `example_merge_from_raw()` state — bundled-cohort reuse (as shipped) vs. a genuine from-raw rebuild (Part 6 finding #6); (c) Batch Correction's normalization method and whether the auto-detect heuristic or an explicit override was used.
- **The exact filter thresholds used** at every stage (missing-data tolerance %, expression/variance percentile cutoffs, `mad_k`, minimum-common-feature floor — noting the last is fixed at 20, not user-set), quoted from the `decisions_ui`/status-line text generated at run time, not from memory or from the UI's default values.
- **The exact batch-correction method and its full configuration** — method (ComBat/limma/SVA/ComBat-seq), ComBat's prior type and mean-only setting if applicable, reference batch if used, which covariates were protected (and any that were automatically dropped for batch-column overlap), and whether the ComBat fallback ladder was triggered (checkable by comparing the requested vs. reported protected-covariate list).
- **Both the before/after PC1-vs-batch p-value and an independent group-preservation check**, not just the correction method's name, as evidence the correction achieved its intended effect.
- **The final merged/corrected sample and gene counts**, alongside the pre-filter counts, so the cumulative effect of every filtering stage is auditable end-to-end.
- **Package versions** for `sva`, `limma`, `edgeR` at the time of the reported analysis, given SVA's documented run-to-run and cross-version reproducibility caveats (Part 6 finding #12).
- **A disclosed acknowledgment of the conservative, intersection-based merge strategy** and what it necessarily discards (every platform-unique feature) whenever reporting a merged cohort's final gene count.
- **If the Data Exploration tab's diagnostics are cited**, note whether the analyzed matrix exceeded the tab's internal sampling caps (`EDA_MAX_POOLED_VALUES` = 200,000 values, `EDA_MAX_SHAPIRO_N` = 5,000), since exact statistics reported for a larger matrix are drawn from an unseeded random subsample.





