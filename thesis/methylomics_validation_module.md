# Methylomics — Validation Module

**Scope of this document.** This document describes exactly one submodule of the Methylomics section of the ArthOMix Shiny application: the **Validation** module ("flask-vial" icon, group "Biomarker modeling"). Every statement below is derived from reading the actual source code. It does **not** describe, and must not be confused with, Transcriptomics Validation, Cross-Tissue Validation, Cross-Ancestry Validation, the Transcriptomics Diagnostic Model, Multi-Omics, Cross-Omics, or any other Methylomics submodule — none of those files were read for this document.

**A note on tab count.** This module implements **11** sub-tabs, not 12. The count is verified directly against `mod_methyl_validation_ui()`. The module's own header comment states it "mirrors `mod_methyl_diagnostic.R`'s per-model tab layout" — and the sibling Diagnostic Classifier module does have 12 tabs (Datasets, Feature Source, Filters & Parameters, 6 model tabs, Model Comparison, Test Internal Data, Export). Validation collapses Diagnostic Classifier's two setup tabs ("Feature Source", "Filters & Parameters" — both meaningless here, since Validation never selects features or trains anything) into a single **Compatibility** tab, and renames "Test Internal Data" to "Test External Data". That yields 2 setup tabs + 6 per-algorithm tabs + 3 trailing tabs = **11**, confirmed by direct line-by-line reading of the `tabsetPanel()` call. Per the documentation brief's own Rule 5 ("code is the authority") and Rule 3 ("no invented tabs"), this document reports 11 tabs rather than fabricating a 12th.

**Primary source file:** [`ArthOMix/R/methylomics/mod_methyl_validation.R`](../ArthOMix/R/methylomics/mod_methyl_validation.R) (777 lines) — contains the module's config, UI, server, and all module-local (`vld_*`) helper functions in one file. There is no separate helpers file.

**Directly reused code from a sibling module** (called by name, not by data file — the app loads every `R/**/*.R` file into one global environment):
- [`ArthOMix/R/methylomics/mod_methyl_diagnostic.R`](../ArthOMix/R/methylomics/mod_methyl_diagnostic.R) — `dxm_*` helper functions and the `DXM_MODEL_SPECS` model registry. The Validation module's own header states: "Reuses `mod_methyl_diagnostic.R`'s `dxm_*` helpers and `DXM_MODEL_SPECS` registry." This is the one code-confirmed dependency on another Methylomics submodule (Diagnostic Classifier); no transcriptomics code is referenced anywhere in this file.
- [`ArthOMix/global.R`](../ArthOMix/global.R) — `load_default_diagnostic_train_test()` (line 679).
- [`ArthOMix/data_paths.R`](../ArthOMix/data_paths.R) — `METH_DIAG_INTERNAL_RDS`, `METH_DIAG_EXTERNAL_RDS`, `METH_DIAG_DATA_AVAILABLE`.
- [`ArthOMix/R/methylomics/parse_upload.R`](../ArthOMix/R/methylomics/parse_upload.R) — `methyl_parse_matrix()`, `methyl_parse_sample_sheet()` (upload-mode parsing).
- [`ArthOMix/R/methylomics/qc.R`](../ArthOMix/R/methylomics/qc.R) — `methyl_sheet_sample_ids()`.
- [`ArthOMix/R/submodules_registry.R`](../ArthOMix/R/submodules_registry.R) — registers the module into the Methylomics tab grid.

---

## 1. Scope and module identity

**Registry entry** (`ArthOMix/R/submodules_registry.R:39-53`, `MX_MODULES` list, entry 12 of 13, immediately after Diagnostic Classifier):
```r
list(config = mod_methyl_validation_config, ui = mod_methyl_validation_ui, server = mod_methyl_validation_server)
```

**Module config** (`mod_methyl_validation.R:8-11`):
```r
mod_methyl_validation_config <- list(
  id = "validation", title = "Validation", icon = "flask-vial", group = "Biomarker modeling",
  description = "External validation of the Diagnostic Classifier's trained models on an independent cohort"
)
```
- **Exact module title:** "Validation"
- **Module id:** `validation`
- **Group:** "Biomarker modeling" (same group as Diagnostic Classifier and Biomarker Card; there is no separate "Validation" group in Methylomics, unlike Transcriptomics)
- **Shiny namespace:** `"mx_validation"` (`ArthOMix/server.R:95`: `lapply(MX_MODULES, function(m) m$server(paste0("mx_", m$config$id), methyl_dataset, methyl_results))`)

**Module purpose (code fact):** the module applies models that were already trained and internally tested in the Diagnostic Classifier module (`results$diagnostic_models`) to an independent external cohort, without any refitting, retuning, or threshold re-optimization. This is stated explicitly in the file's own header comment (lines 1-6):

> "Applies the Diagnostic Classifier's already-trained models (`results$diagnostic_models`) to an independent EXTERNAL cohort - never retrains... every 'test' evaluation here is scored on the external validation cohort loaded on the Datasets tab, not the internal train/test split."

The module performs **no model training of its own** — no `caret::train()`, `xgboost::xgb.train()`, `glm()`, or any other fitting call appears anywhere in `mod_methyl_validation.R`. Its only computational operations are: loading/filtering an external cohort, checking CpG-feature compatibility, calling `predict()` on an already-fitted model object, and computing performance metrics (ROC/AUC, confusion-matrix statistics, bootstrap confidence intervals, calibration) on the resulting predictions.

---

## 2. Web-application implementation

`mod_methyl_validation_ui(id)` (`mod_methyl_validation.R:303-367`) builds one `tabsetPanel(id = ns("main_tabs"), type = "tabs")`. `mod_methyl_validation_server(id, dataset, results)` (`mod_methyl_validation.R:371-777`) is invoked as `mod_methyl_validation_server("mx_validation", methyl_dataset, methyl_results)`, receiving the shared `methyl_results` reactiveValues (its `diagnostic_models` field is written exclusively by the Diagnostic Classifier module) as its `results` argument.

Reactive data flow: `avail_models <- reactive({ results$diagnostic_models %||% list() })` (line 388) is the single source of every trained model the module can evaluate. A `cohort` `reactiveValues` object (lines 375-381) holds the currently loaded external validation matrix/labels once the user loads one on the Datasets tab; every other tab reads from `cohort` and `avail_models()`. Each of the six per-algorithm tabs owns its own `vms` (`reactiveValues`, lines 642-648) holding that algorithm's external-run state, and a shared `vruns` (`reactiveValues`, line 383) accumulates one completed run per algorithm for the Model Comparison/Export tabs. User controls are `radioButtons`, `fileInput`, `selectInput`/`selectizeInput`, and `actionButton`/`downloadButton`; validation/error messages are surfaced via `showNotification()` (through the module-local `vld_notify_fail()` wrapper, lines 107-111) inside `observeEvent()` handlers, and via a single `validate(need(...))` call (line 700) inside a render context. Tables render via `DT::renderDataTable`/`DT::datatable`; figures via `renderPlot`; dynamic tab bodies via `renderUI`/`uiOutput`, each wrapped in `shinycssloaders::withSpinner()`.

The six per-algorithm tabs are generated programmatically — `vld_render_model_panel()` (UI, lines 165-230) and `vld_register_model_server()` (server, lines 233-299) are each called once per entry of `DXM_MODEL_SPECS` (`mod_methyl_diagnostic.R:437-530`), so their tab labels are exactly that registry's `label` fields, with input/output IDs of the form `<mid>_...` (e.g. `lr_run_btn`, `svm_ext_metrics_table`).

---

## 3. Data loading and data sources

**Preloaded data.** `METH_DIAG_DATA_AVAILABLE` (`data_paths.R:111`) gates the preloaded path:
```r
METH_RAW_DATA_ROOT <- normalizePath(Sys.getenv("METH_RAW_DATA_ROOT", get_preloaded_path("methylomics", "matrix")), ...)
METH_DIAG_INTERNAL_RDS <- file.path(METH_RAW_DATA_ROOT, "gse42861_internal_panel_celltype_adjusted.rds")
METH_DIAG_EXTERNAL_RDS <- file.path(METH_RAW_DATA_ROOT, "gse111942_external_panel.rds")
METH_DIAG_DATA_AVAILABLE <- file.exists(METH_DIAG_INTERNAL_RDS) && file.exists(METH_DIAG_EXTERNAL_RDS)
```
Both objects are read once per process by `load_default_diagnostic_train_test()` (`global.R:679-692`), which caches the result in `.arthomix_cache`, and returns `list(internal = ..., external = ...)`. The Validation module's "Load Preloaded External Cohort" handler uses only `dd$external` (`mod_methyl_validation.R:453-455`). Provenance, from `global.R:649-670`:
- `gse42861_internal_panel_celltype_adjusted.rds` — GSE42861 reprocessed via `minfi::preprocessNoob()`, granulocyte-adjusted (`lm(M ~ Neutro + Eosino)` per CpG), 689 samples, both sexes; this is the internal cohort the Diagnostic Classifier trains and internally tests on (75/25 split). **Never opened by the Validation module directly.**
- `gse111942_external_panel.rds` — an independent cohort, Noob-processed the same way, subset to the panel CpGs, 43 samples, **all female**. This is what the Validation module's preloaded path loads. The UI states this explicitly (`mod_methyl_validation.R:324`): *"43 samples, all female, Noob-renormalized and granulocyte-adjusted."*

If `METH_DIAG_DATA_AVAILABLE` is `FALSE`, the preloaded path shows an error notification ("The preloaded external cohort isn't available in this deployment.", line 451) and the deployment falls back to Upload-only.

**Uploaded data.** Two `fileInput`s: `val_upload_matrix` ("Methylation matrix (CSV/TSV, probes x samples)") and `val_upload_sheet` ("Sample sheet (CSV/TSV)"). Parsing:
- `methyl_parse_matrix()` (`parse_upload.R:8-34`) — `data.table::fread()`; first column is the probe ID (must be unique, checked at line 522: `sum(duplicated(rownames(mat))) == 0`); remaining columns coerced to a numeric matrix.
- `methyl_parse_sample_sheet()` (`parse_upload.R:37-43`) — `data.table::fread()` into a data frame, one row per sample.
- `methyl_sheet_sample_ids()` (`qc.R:355-360`) — resolves sample IDs from a `sample`/`Sample`/`sample_id`/`Sample_ID` column if present, else assumes row order matches matrix column order.

**In-memory object consumed from shared state.** `results$diagnostic_models` — a list keyed by `paste(mid, m$analysis_type, paste(m$feature_ids, collapse="|"))`, populated exclusively by the Diagnostic Classifier module. This is the *only* source of trained models the Validation module can evaluate; it never trains its own.

**How data are stored/passed through the module:** the loaded external cohort lives in the module-local `cohort` reactiveValues (`m` = feature × sample M-value matrix, `y` = a two-level factor, plus label/provenance metadata); a completed per-algorithm external run's outputs live in that algorithm's `vms` reactiveValues; completed runs across all algorithms accumulate in the shared `vruns` reactiveValues for cross-model comparison and export.

---

## 4. Input data

**Methylation representation.** The module works exclusively in **M-values** (logit of beta: `log2(b / (1-b))`, clipped to `[1e-6, 1-1e-6]` via `dxm_beta_to_m()`, `mod_methyl_diagnostic.R:26-29`). The preloaded cohort is stored as beta values and is always converted; an uploaded matrix is converted unless the user declares it is already M-values (`val_upload_scale` radio button, "Beta-value"/"M-value", default "Beta-value").

**Features.** Rows of the validation matrix are **CpG probe IDs** (rownames). The module never selects, ranks, or filters CpGs itself — the required CpG set for any given model is fixed by whatever was trained in Diagnostic Classifier (`m$feature_ids`), and the Validation module only checks whether the loaded cohort's row names contain that exact set.

**Samples / phenotype.** Columns are samples. The outcome label is a two-level factor built by comparing a phenotype/class column to the reference model's `ref_level`/`comp_level` (internally recoded to `Class0`/`Class1` — `DXM_POS`/`DXM_NEG`, `mod_methyl_diagnostic.R:67-68`). For the preloaded cohort, the phenotype column is `ext$pheno$group`; for an uploaded cohort, it is a user-selected sample-sheet column (`val_upload_pheno_col`, auto-guessed by regex `"class|phenotype|group|status|disease"`).

**Sex.** An optional stratification field, not a modeling covariate: `radioButtons` (`cohort_sex_stratum` / `val_upload_sex_stratum`, choices "All samples"/"Female"/"Male", default "All samples") filter which samples are loaded into the cohort. For uploads, sex values are read from a user-selected column (`val_upload_sex_col`, auto-guessed by regex `"sex|gender"`) and normalized to "F"/"M" by `dxm_normalize_sex()`.

**User-defined cutoff.** None. See §7.

**Train/test ratio.** Not applicable — the module performs no splitting (see §6).

**Classifier parameters.** None are exposed here; all six classifiers' hyperparameters are set in the Diagnostic Classifier module and are already baked into the fitted model objects this module consumes.

Assumptions explicitly **not** made because the code does not perform them: no normalization, imputation, or annotation-file lookup happens inside this module (missing-value samples are dropped, never imputed — see §6); CpGs are not mapped to genes anywhere in this file.

---

## 5. All 11 tabs

### Tab 1 — "Datasets"

**Purpose:** lists every trained model published by Diagnostic Classifier and loads the external validation cohort (preloaded or uploaded).

**Input data:** `results$diagnostic_models`; either the preloaded `gse111942_external_panel.rds` cohort or an uploaded matrix + sample sheet.

**User controls:**
- `cohort_source` — `radioButtons`, "Preloaded external validation cohort" / "Upload an independent cohort", default `"preloaded"`
- `cohort_sex_stratum` — `radioButtons`, "All samples"/"Female"/"Male", default `"all"` (preloaded path)
- `load_preloaded_btn` — `actionButton`, "Load Preloaded External Cohort"
- `val_upload_matrix` — `fileInput`, methylation matrix (CSV/TSV)
- `val_upload_scale` — `radioButtons`, "Beta-value"/"M-value", default `"beta"`
- `val_upload_sheet` — `fileInput`, sample sheet (CSV/TSV)
- `val_upload_pheno_col` — dynamic `selectInput`, auto-guessed phenotype column
- `val_upload_sex_stratum` — `radioButtons`, default `"all"`
- `val_upload_sex_col` — dynamic `selectInput`, auto-guessed sex column
- `load_upload_btn` — `actionButton`, "Load Uploaded Cohort"

**Processing:** on `load_preloaded_btn`: requires ≥1 trained model and `METH_DIAG_DATA_AVAILABLE`; filters the preloaded external panel by sex stratum and by the most-recently-tested model's class labels; requires >5 matching samples; converts beta→M via `dxm_beta_to_m()`; builds the outcome factor. On `load_upload_btn`: parses the uploaded matrix/sheet; requires ≥3 common sample IDs between matrix and sheet; requires both the reference model's class labels present in the chosen phenotype column; optional sex filter (requires ≥1 match); requires ≥3 samples after class/sex filtering; rejects duplicated CpG rownames; converts to M-values unless the user declared M-value scale.

**Output data:** trained-model summary table; a training-vs-validation-cohort comparison table; a sample-ID-overlap check between training and validation cohorts (only when real sample IDs are available); a per-session load-history table.

**Output format:** tables (DT), text/alert UI.

**Execution condition:** model summary is automatic once models exist; cohort loading requires clicking `load_preloaded_btn` or `load_upload_btn`; downstream comparison/overlap tables appear once a cohort is loaded.

---

### Tab 2 — "Compatibility"

**Purpose:** audits CpG-feature overlap between every trained model and the loaded validation cohort.

**Input data:** `avail_models()`, `cohort$m` (rownames = CpGs available in the validation cohort).

**User controls:** `compat_model_select` — `selectInput`, populated dynamically with every trained model.

**Processing:** for every model, `vld_feature_alignment(m$feature_ids, rownames(cohort$m))` computes intersect/setdiff of required vs. available CpGs and an overlap percentage; `ok` requires 100% of required CpGs present.

**Output data:** an overview table (Model, Analysis type, Required CpGs, Matched, Missing, Overlap %, Status) across all models; a 3-state banner (all-blocked / some-blocked / all-pass); a per-model KPI row (Required CpGs, Available, Matched, Missing, Overlap %); a 7-row compatibility audit table (`vld_compat_table()`) comparing "Training Pipeline" to "Validation Requirement/Applied" for: Data representation, Feature/CpG set, Scaling, Normalization/batch correction, Missing values, Classification threshold, Platform/array type — each row carries a Pass/Fail/Warning status.

**Output format:** tables (DT), alert banner, valueBoxes.

**Execution condition:** automatic once a cohort is loaded; the detailed per-model audit table updates on selecting a model in `compat_model_select`.

---

### Tabs 3–8 — Per-algorithm tabs: "Logistic Regression", "Elastic Net", "Support Vector Machine", "Random Forest", "Gradient Boosting / XGBoost", "k-Nearest Neighbors"

These six tabs share identical structure (generated by `vld_render_model_panel()`/`vld_register_model_server()`, parameterized by `DXM_MODEL_SPECS` entries `lr`, `enet`, `svm`, `rf`, `gbm`, `knn`).

**Purpose:** score that algorithm's already-trained model against the loaded external validation cohort and report performance.

**Input data:** the trained model object for that algorithm (`m$fit$model`, `m$feature_ids`, `m$threshold`, `m$train_metrics`, `m$cv_roc`) from `results$diagnostic_models`; the loaded `cohort` (`m`/`y`).

**User controls:**
- `<mid>_run_select` — `selectInput`, "Trained run" (shown only if more than one trained run exists for that algorithm)
- `<mid>_run_btn` — `actionButton`, "Run External Validation"
- `<mid>_roc_btn` — `actionButton`, "Generate ROC/AUC"
- `<mid>_calib_btn` — `actionButton`, "Generate Calibration"
- `<mid>_roc_png` — `downloadButton`, "Download ROC plot (PNG)"

**Processing:** on `<mid>_run_btn`, `vld_do_run_external()` (`mod_methyl_validation.R:116-160`) runs the module's core scoring routine (full logic in §6): guardrail checks, `predict()` via `dxm_predict_prob(m$fit$model, Xv)` (no refitting), ROC/AUC via `dxm_roc_bundle()` (`pROC::roc`), a metrics bundle via `dxm_metrics_bundle()` (includes `PRROC::pr.curve` PR-AUC), a confusion matrix at the model's fixed training-derived threshold via `dxm_confusion()`, and three bootstrap confidence intervals (`vld_bootstrap_ci()`, 1000 resamples each, seeds `m$seed`, `m$seed+1`, `m$seed+2`) for sensitivity/specificity/accuracy. `<mid>_roc_btn` reveals the ROC plots (no recomputation). `<mid>_calib_btn` computes calibration via `dxm_calibration()` (`glm(y ~ prob, family = binomial())` for slope/intercept, plus a 10-bin reliability table and Brier score).

**Output data:** training metrics table and confusion matrix (as already published by Diagnostic Classifier); external-test metrics table (Accuracy, Balanced accuracy, Sensitivity, Specificity, Precision, NPV, F1, AUC with 95% CI, Brier, PR-AUC) and confusion matrix; 4 valueBoxes for training (Training AUC, Mean CV AUC ± SD, Classification threshold, Feature count) and 4 for external test (External AUC with 95% DeLong CI, Sensitivity with 95% bootstrap CI, Specificity with 95% bootstrap CI, Validated N split by class); a training-vs-external overfitting note (`dxm_overfitting_note()`); Training ROC plot, Cross-Validated ROC plot, External Test ROC plot; a calibration plot; a downloadable ROC PNG.

**Output format:** valueBoxes, DT tables, ggplot2 plots (`renderPlot`), a downloadable PNG file.

**Execution condition:** Setup box and Results: Training box appear automatically once a trained model exists; Results: External Test / ROC-AUC / Diagnostics boxes appear only after `<mid>_run_btn` is clicked and succeeds; the ROC plots require an additional click of `<mid>_roc_btn`; the calibration plot requires an additional click of `<mid>_calib_btn`.

---

### Tab 9 — "Model Comparison"

**Purpose:** compares external-validation performance across every algorithm that has completed a run.

**Input data:** `vruns` (one entry per completed external-validation run, across algorithms).

**User controls:**
- `compare_select` — `selectizeInput`, multi-select of completed runs
- `compare_curve` — `selectInput`, "External Test"/"Training"/"Cross-Validated", default "External Test"
- `compare_roc_btn` — `actionButton`, "Generate ROC Comparison"
- `compare_download` — `downloadButton`, "Download comparison (CSV)"

**Processing:** builds a comparison table across all completed runs (Model, Feature set, Train AUC, CV AUC, External AUC, N validated, Threshold, Ran at); on `compare_roc_btn`, `validate(need(length(sel) > 0, "Select at least one run to compare."))` gates the action, then builds ROC bundles for the selected runs/curve type and plots them via `dxm_plot_roc_compare()`.

**Output data:** comparison table; overlaid ROC comparison plot; a CSV export.

**Output format:** DT table, ggplot2 plot, downloadable CSV (`methylomics_validation_model_comparison.csv`).

**Execution condition:** table is automatic once ≥1 run has completed; the comparison ROC plot requires selecting ≥1 run and clicking `compare_roc_btn`.

---

### Tab 10 — "Test External Data"

**Purpose:** a read-only snapshot of the currently loaded external cohort and its CpG-feature availability per trained model.

**Input data:** `cohort` (loaded matrix/labels), `avail_models()`.

**User controls:** none beyond what is set on the Datasets tab.

**Processing:** reports sample counts by class, CpG rows available, missing-value count in the cohort matrix; for each trained model, re-runs `vld_feature_alignment()` to report Training features / Available in external cohort / Shared features / Unmatched (dropped) features.

**Output data:** summary text/list; a per-model feature-availability table.

**Output format:** text summary, DT table.

**Execution condition:** automatic once a cohort is loaded; otherwise shows a muted prompt to load one on the Datasets tab.

---

### Tab 11 — "Export"

**Purpose:** downloads of accumulated external-validation results.

**Input data:** `vruns` (completed runs), `avail_models()`, `cohort`.

**User controls:**
- `export_metrics_csv` — `downloadButton`, "Download all metrics (CSV)"
- `export_featavail_csv` — `downloadButton`, "Download feature availability per model (CSV)"

**Processing:** on download, assembles the same per-run metrics columns as the Model Comparison CSV, and a per-model feature-availability table (features required/available/matched/missing).

**Output data:** two CSV files.

**Output format:** downloadable files — `methylomics_validation_metrics.csv` (model, analysis_type, n_features, features, threshold, train_auc, cv_auc, external_auc, n_validated, ran_at); `methylomics_validation_feature_availability.csv` (model, analysis_type, features_required, features_available, features_matched, features_missing).

**Execution condition:** available once ≥1 external-validation run has completed; otherwise shows a muted prompt.

---

## 6. Validation workflow

The module's actual pipeline, as implemented, is:

**Data loading** (Datasets tab: preloaded `gse111942_external_panel.rds` or uploaded matrix + sheet) → **input validation** (sample-count minimums, class-label match to the reference model, duplicate-CpG check, sex-stratum availability) → **user selection/parameterization** (cohort source, sex stratum, phenotype/sex column choice for uploads, which trained run to evaluate) → **preprocessing** (beta→M-value conversion via `dxm_beta_to_m()`; complete-case filtering — samples missing any required CpG are dropped, never imputed) → **methylomics feature preparation** (subset the validation matrix to exactly the trained model's `feature_ids`, via `vld_feature_alignment()`, requiring 100% overlap) → **model application** (`dxm_predict_prob(m$fit$model, Xv)` — predict-only, no refitting) → **statistical evaluation** (`pROC::roc()`/`pROC::ci.auc()` for AUC and its DeLong CI; `dxm_confusion()` for sensitivity/specificity/precision/NPV/F1/accuracy/balanced-accuracy/MCC at the model's fixed threshold; `PRROC::pr.curve()` for PR-AUC; `vld_bootstrap_ci()` — 1000-resample percentile bootstrap — for sensitivity/specificity/accuracy CIs; `dxm_calibration()` — `glm(y ~ prob, family = binomial())` — for calibration slope/intercept and a reliability table) → **visualization/table generation** (ROC plots, calibration plot, metrics/confusion tables, cross-model comparison table/plot) → **output** (on-screen tables/plots plus CSV/PNG downloads).

**No train/test split is performed by this module.** It is architecturally forbidden from splitting or refitting by design (header comment: "never retrains"). The only train/test split in the methylomics diagnostic pipeline lives entirely inside the Diagnostic Classifier module (`mod_methyl_diagnostic.R:899-965`): `train_frac` defaults to `0.75` (range 0.5-0.9), a fixed random seed (`dxm_seed`, default `42`), and a stratified split via `caret::createDataPartition(y_all, p = train_frac, list = FALSE)`. That split, and the model fitted on its training partition, are already embedded in the `m$fit` object the Validation module consumes — Validation applies that fitted pipeline to a wholly separate cohort file (`gse111942_external_panel.rds`), never to a re-split of the training data.

---

## 7. Diagnostic classifier distinction (models consumed by this module)

The Validation module does not train classifiers; it evaluates the six classifiers already fitted in the Methylomics Diagnostic Classifier module:

- **Features entering the classifier:** the fixed CpG panel(s) chosen upstream in Diagnostic Classifier (single CpG or a multi-CpG panel from WGCNA/Feature Selection/manual pick) — never re-selected by Validation.
- **Outcome predicted:** a binary class label (`m$ref_level` vs. `m$comp_level`, internally `Class0`/`Class1`).
- **Classifiers implemented (6):** Logistic Regression, Elastic Net, Support Vector Machine, Random Forest, Gradient Boosting/XGBoost, k-Nearest Neighbors — the exact set in `DXM_MODEL_SPECS`.
- **Training vs. validation separation:** strictly separate. Training (and the 75/25 internal test) happen only in Diagnostic Classifier, on `gse42861_internal_panel_celltype_adjusted.rds`. Validation evaluates the resulting fitted objects only against `gse111942_external_panel.rds` (or an uploaded cohort) — an entirely different, independent set of samples.
- **Train/test ratio and whether it's user-adjustable:** 75% training / 25% internal test, set by `train_frac` (`numericInput`, default `0.75`, range 0.5-0.9, step 0.05) in Diagnostic Classifier — adjustable there, not in Validation, and not re-run by Validation.
- **Cutoff / classification threshold:** there is **no user-adjustable threshold anywhere in the Validation module**. The threshold used is always the fixed `m$threshold` computed once in Diagnostic Classifier via `dxm_pick_threshold()` from a `threshold_strategy` choice ("Default (0.50)", "Youden's J", "Sensitivity-focused ≥0.90 sens", "Specificity-focused ≥0.90 spec") — computed strictly from the **training ROC only**, never recomputed on internal-test or external-validation data. The Compatibility tab's audit table states this design guarantee verbatim: *"Reused unchanged - never re-optimized on validation labels."*
- **Cutoff source:** training data only, never the validation/test set — confirmed both by the threshold-strategy code (`dxm_pick_threshold`, operates on the training ROC bundle) and by the Compatibility audit table's explicit claim.
- **Prediction output:** class probabilities from `dxm_predict_prob()`, thresholded into a binary prediction for the confusion matrix.
- **Confusion matrix:** yes — computed for both the training set (as published by Diagnostic Classifier) and the external test set, via `dxm_confusion()`.
- **Sensitivity, Specificity, Accuracy, Precision (as NPV/precision), Balanced accuracy, F1, MCC:** all present, computed by `dxm_confusion()`.
- **Recall:** identical to sensitivity in this codebase (`dxm_metrics_display` labels it "Sensitivity (recall)").
- **ROC / AUC:** present — `pROC::roc()`, plotted for Training, Cross-Validated, and External Test curves separately.
- **AUC 95% confidence interval:** present — `pROC::ci.auc(r, method = "delong")`, falling back to `method = "bootstrap", boot.n = 1000` if DeLong fails.
- **Confidence intervals for sensitivity/specificity/accuracy:** present — ordinary (with-replacement) 1000-resample percentile bootstrap, `vld_bootstrap_ci()`, requiring ≥50 usable resamples or returning `NA`.
- **Calibration:** present — `dxm_calibration()`: a 10-bin reliability table, Brier score, and a logistic-calibration slope/intercept from `glm(y ~ prob, family = binomial())`.
- **PR-AUC:** present — `PRROC::pr.curve()`, included in the metrics bundle.
- **Any other metric:** none found beyond the above.

**Not implemented:** `cutpointr` is not used anywhere in this module or in Diagnostic Classifier; `caret::confusionMatrix()` is never called (confusion-matrix statistics are computed by hand in `dxm_confusion()`); no `glm()`/logistic-regression fitting occurs inside Validation itself (the only `glm()` call in this module's dependency chain is the calibration-slope fit, which does not classify samples).

External validation, as implemented, is **internal validation of Diagnostic Classifier's own training/internal-test split, followed by application to one specific held-out external file** — the module never claims "external validation" against multiple independent cohorts; it applies to exactly one bundled external file (or one user upload) at a time.

---

## 8. Validation strategy

The Validation module implements exactly one validation procedure: apply an already-fitted, already-thresholded model to a cohort that was never used in its training. Distinguishing the roles precisely, as implemented:

- **Training data:** `gse42861_internal_panel_celltype_adjusted.rds`, split 75%/25% inside Diagnostic Classifier (`train_frac`, default 0.75, `caret::createDataPartition`, `set.seed(42)` by default). This is an **internal** train/internal-test procedure, entirely inside Diagnostic Classifier — it is not the Validation module's own internal validation, but Validation's fitted-model input depends on it.
- **Independent validation data (this module's own scope):** `gse111942_external_panel.rds` (or an uploaded cohort) — a file distinct from the training file, drawn from a different GEO accession, never used for training or tuning. The cohort-comparison table's own label for this is literal: *"Independent validation - never used for training or tuning"* (`mod_methyl_validation.R:550`).
- **Fitted model:** `m$fit$model` — the exact `caret::train`/`xgb.Booster` object published by Diagnostic Classifier; reused verbatim via `predict()`, never refit.
- **Prediction:** class probabilities from `dxm_predict_prob()`.
- **Cutoff:** `m$threshold` — fixed at training time, never recomputed on validation data (§7).
- **Final classification:** probability ≥ `m$threshold` → positive class, via `dxm_confusion()`.
- **Validation metrics:** AUC (with DeLong or bootstrap 95% CI), sensitivity/specificity/accuracy (each with a 1000-resample bootstrap 95% CI), precision, NPV, F1, balanced accuracy, MCC, PR-AUC, Brier score, calibration slope/intercept.

Because the "training" cohort and the "validation" cohort here are two different, independently collected datasets (different GEO accessions) rather than two random splits of one dataset, this procedure is more accurately described as **evaluation on an independent external cohort using a fixed, previously-trained model** — not a k-fold or repeated train/test resampling scheme. No cross-tissue, cross-ancestry, or transcriptomic validation is implemented anywhere in this file; the only other dataset ever referenced is the training cohort embedded in the consumed model object.

---

## 9. Output data

Consolidated across all tabs (see §5 for per-tab detail):

- Trained-model summary table (Algorithm, Analysis type, Features, Contrast, Stratum, Train AUC, CV AUC, Internal-test AUC, Tested at)
- Training-vs-validation-cohort comparison table
- Sample-ID overlap / independence check
- Cohort load-history table
- Compatibility overview table (per model: Required CpGs, Matched, Missing, Overlap %, Status) and a 7-row detailed compatibility audit table
- Per-algorithm (×6): training metrics table, external-test metrics table (Accuracy, Balanced accuracy, Sensitivity, Specificity, Precision, NPV, F1, AUC + 95% CI, Brier, PR-AUC), training and external confusion matrices, Training/Cross-Validated/External ROC plots, calibration plot, downloadable ROC PNG
- Model Comparison table (Model, Feature set, Train AUC, CV AUC, External AUC, N validated, Threshold, Ran at) and overlaid ROC-comparison plot
- Test External Data snapshot (sample counts, CpG row count, missing-value count) and per-model feature-availability table
- Two Export CSVs (all metrics; feature availability per model)

All outputs are either DT tables, ggplot2 plots rendered via `renderPlot`, dynamic UI text/alerts via `renderUI`, or downloadable CSV/PNG files via `downloadHandler`. No model-object or raw-prediction downloads are implemented — only summary metrics and plots are exportable.

---

## 10. XomicShiny-style short implementation description

The Methylomics Validation module is a Shiny submodule (`mod_methyl_validation.R`, id `validation`, namespace `mx_validation`) that evaluates classifiers already trained by the Diagnostic Classifier submodule against an independent methylation cohort, without refitting. Methylation data enter the module either as a bundled preloaded external cohort (`gse111942_external_panel.rds`; 43 samples, all female, Noob-normalized, granulocyte-adjusted, subset to trained-panel CpGs) or as a user-uploaded probe × sample beta-/M-value matrix (CSV/TSV) with an accompanying sample sheet (CSV/TSV), both parsed with `data.table::fread()`. The interface exposes 11 tabs: Datasets (model discovery and cohort loading), Compatibility (CpG-feature overlap audit against each trained model), one tab per algorithm (Logistic Regression, Elastic Net, Support Vector Machine, Random Forest, Gradient Boosting/XGBoost, k-Nearest Neighbors), Model Comparison, Test External Data, and Export. For each algorithm, clicking "Run External Validation" calls `predict()` on the stored fitted model at its fixed, training-derived classification threshold, then computes ROC/AUC (`pROC`, DeLong 95% CI), a confusion matrix and its derived statistics (sensitivity, specificity, precision, NPV, F1, balanced accuracy, MCC), PR-AUC (`PRROC`), bootstrap 95% confidence intervals (1000 resamples) for sensitivity/specificity/accuracy, and, on request, a calibration curve (10-bin reliability table, Brier score, logistic calibration slope/intercept). Results render as DT tables and ggplot2 plots and can be exported as CSV (per-model and cross-model comparison) or PNG (ROC curve); no functionality beyond these code-implemented computations is exposed to the user.

---

## 11. Thesis paragraph

The Methylomics Validation module was implemented as a Shiny-based interactive workflow for applying the Diagnostic Classifier submodule's already-trained CpG-based classifiers to an independent, previously unseen methylation cohort, without any refitting or threshold re-optimization. The module accepts either a bundled external validation cohort (GSE111942, 43 samples, all female, Noob-normalized and granulocyte-adjusted, subset to the trained CpG panel) or a user-uploaded methylation matrix and sample sheet, and provides controls for cohort source, sex stratum, and phenotype/sex column selection. Across its 11 code-defined tabs, the implementation supports discovery of every trained model, an automated audit of CpG-feature compatibility between each model and the loaded cohort, per-algorithm external scoring for six classifiers (logistic regression, elastic net, support vector machine, random forest, gradient boosting/XGBoost, and k-nearest neighbors), cross-model performance comparison, and a read-only snapshot of the loaded validation data. The processing workflow consists of beta-to-M-value conversion, complete-case filtering of samples with missing required CpGs, and exact CpG-feature alignment to each trained model's fixed panel, followed by prediction with the stored fitted model at its training-derived classification threshold and evaluation via ROC/AUC (with DeLong confidence intervals), a confusion matrix and its derived statistics, PR-AUC, bootstrap confidence intervals for sensitivity/specificity/accuracy, and calibration analysis. The resulting outputs include per-model metrics and confusion-matrix tables, ROC and calibration plots, a cross-model comparison table and overlaid ROC plot, and downloadable CSV/PNG exports, allowing users to assess how each trained methylation classifier generalizes to an independent cohort under a fixed, previously-established decision threshold.

---

## 12. Code-to-thesis mapping

| Component | Code location | Input | Processing | Output |
|---|---|---|---|---|
| Tab 1 — Datasets | `mod_methyl_validation.R:312-346` (UI), `:407-583` (server) | `results$diagnostic_models`; preloaded `gse111942_external_panel.rds` or uploaded matrix+sheet | Sex/class filtering, beta→M conversion, sample matching, `vld_sample_overlap()` | Model list table, cohort comparison table, load-history table |
| Tab 2 — Compatibility | `mod_methyl_validation.R:348-358` (UI), `:587-637` (server) | `avail_models()`, `cohort$m` rownames | `vld_feature_alignment()`, `vld_compat_table()` | Overview + detailed audit tables, status banner, KPI valueBoxes |
| Tab 3 — Logistic Regression | `mod_methyl_validation.R:165-299` (generic), spec `mod_methyl_diagnostic.R:438-453` | `m$fit$model` (glmnet/glm), `cohort` | `dxm_predict_prob()`, `dxm_roc_bundle()`, `dxm_confusion()`, `vld_bootstrap_ci()`, `dxm_calibration()` | Metrics/confusion tables, ROC/calibration plots, ROC PNG |
| Tab 4 — Elastic Net | same functions, spec `mod_methyl_diagnostic.R:454-476` | `m$fit$model` (glmnet/glm) | same | same |
| Tab 5 — Support Vector Machine | same functions, spec `mod_methyl_diagnostic.R:477-493` | `m$fit$model` (svmLinear/Radial/Poly) | same | same |
| Tab 6 — Random Forest | same functions, spec `mod_methyl_diagnostic.R:494-508` | `m$fit$model` (rf) | same | same |
| Tab 7 — Gradient Boosting / XGBoost | same functions, spec `mod_methyl_diagnostic.R:509-520` | `m$fit$model` (xgb.Booster) | `dxm_predict_prob()` xgb.DMatrix path | same |
| Tab 8 — k-Nearest Neighbors | same functions, spec `mod_methyl_diagnostic.R:521-529` | `m$fit$model` (knn) | same | same |
| Tab 9 — Model Comparison | `mod_methyl_validation.R:362` (UI), `:653-711` (server) | `vruns` (all completed runs) | Aggregation, `dxm_plot_roc_compare()` | Comparison table, comparison ROC plot, CSV export |
| Tab 10 — Test External Data | `mod_methyl_validation.R:363` (UI), `:715-739` (server) | `cohort`, `avail_models()` | `vld_feature_alignment()` per model | Cohort snapshot text, feature-availability table |
| Tab 11 — Export | `mod_methyl_validation.R:364` (UI), `:743-773` (server) | `vruns`, `avail_models()`, `cohort` | Table assembly | Two CSV downloads |

---

## 13. Documentation limitations

- Exact byte sizes of the preloaded RDS files were not independently re-verified for this document (a sibling document, `methylomics_Diagnostic_Classifier.md`, cites 131KB/7.6KB for the two files; this was not re-checked here since it is outside the Validation module's own code).
- The internal contents of `gse111942_external_panel.rds` (e.g. exact CpG identifiers, exact sample IDs) were not inspected — only the code paths that read and transform it were read.
- Whether any deployment currently has `METH_DIAG_DATA_AVAILABLE == TRUE` (i.e., whether the preloaded files actually exist on disk in this environment) was not checked; this document describes the code path, not a specific runtime's file presence.
