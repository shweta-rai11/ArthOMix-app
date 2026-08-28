# Methylomics — Diagnostic Classifier Module

**Scope of this document.** This document describes exactly one submodule of the Methylomics section of the ArthOMix Shiny application: the **Diagnostic Classifier** ("stethoscope" icon, group "Biomarker modeling"). All statements are derived from reading the actual source code, the actual bundled data files, and the actual methods write-up shipped alongside the underlying reference pipeline. It does **not** describe, and must not be confused with, the unrelated Transcriptomics `mod_diagnostic.R` tab — that file was not read for this document and none of its content, classifiers, gene panels, or results appear here.

**Primary source file:** [`ArthOMix/R/methylomics/mod_methyl_diagnostic.R`](../ArthOMix/R/methylomics/mod_methyl_diagnostic.R) (1,297 lines) — this is the entire module; there is no second file defining Diagnostic-Classifier-specific UI or server logic. The file's own header comment states this explicitly:

> "Diagnostic Classifier submodule (script09_diagnostic_classifier). Methylomics only - transcriptomics' `mod_diagnostic.R` is a separate, unrelated tab."

**Directly reused helper functions** (called by name, not sourced separately — the app loads every `R/**/*.R` file into one global environment):
- [`ArthOMix/global.R`](../ArthOMix/global.R) — `load_default_diagnostic_train_test()` (line 679), `load_default_wgcna_module_assignment()` (line 515), `load_default_diagnostic_ensemble_votes()` (line 621).
- [`ArthOMix/data_paths.R`](../ArthOMix/data_paths.R) — `METH_DIAG_INTERNAL_RDS`, `METH_DIAG_EXTERNAL_RDS`, `METH_DIAG_DATA_AVAILABLE`, `METH_DIAGNOSTIC_VOTES_DIR`, `METH_DIAGNOSTIC_DIR`.
- [`ArthOMix/R/methylomics/parse_upload.R`](../ArthOMix/R/methylomics/parse_upload.R) — `methyl_parse_matrix()`, `methyl_parse_sample_sheet()`, `methyl_parse_probe_list()` (upload-mode parsing only).
- [`ArthOMix/R/methylomics/qc.R`](../ArthOMix/R/methylomics/qc.R) — `methyl_sheet_sample_ids()` (sample-ID resolution for uploaded sample sheets).
- [`ArthOMix/R/submodules_registry.R`](../ArthOMix/R/submodules_registry.R) — registers the module into the Methylomics tab grid (line 50).

**Bundled data referenced in this document:**
- `data/preloaded/methylomics/raw/…/gse42861_internal_panel_celltype_adjusted.rds` and `gse111942_external_panel.rds` (paths resolved via `METH_RAW_DATA_ROOT`) — the actual matrices this module trains and internally tests on.
- `data/preloaded/methylomics/tables/script07_ml_feature_selection/tables/ensemble_votes_{female,male}.csv` — consumed live by this module's "Feature Selection" feature source.
- `data/preloaded/methylomics/tables/script05_wgcna_sexstratified/tables/module_assignment_{female,male}[_merged10].csv` — consumed live by this module's "WGCNA" feature source.
- `data/preloaded/methylomics/tables/script09_diagnostic_classifier/tables/{diagnostic_panel_auc,diagnostic_perprobe_auc,diagnostic_perprobe_all_algorithms,diagnostic_panel_importance}_{female,male}.csv` and `METHODS_diagnostic_classifier.md` — bundled reference outputs from the underlying research pipeline. **These CSVs are not read by any R code in this deployment** (verified by repository-wide grep; see §12). They document the analysis the live module reproduces the *engine* of, not files the module opens.

---

## 1. Scope and module identity

**Registry entry** ([`submodules_registry.R:50`](../ArthOMix/R/submodules_registry.R#L50)):
```r
list(config = mod_methyl_diagnostic_config, ui = mod_methyl_diagnostic_ui, server = mod_methyl_diagnostic_server)
```

**Module config** ([`mod_methyl_diagnostic.R:15-18`](../ArthOMix/R/methylomics/mod_methyl_diagnostic.R#L15-L18)):
```r
mod_methyl_diagnostic_config <- list(
  id = "diagnostic", title = "Diagnostic Classifier", icon = "stethoscope", group = "Biomarker modeling",
  description = "Diagnostic classifiers (logistic regression, elastic net, SVM, random forest, XGBoost, kNN) on single CpGs or a panel, with cross-validated tuning."
)
```
- **Exact module title:** "Diagnostic Classifier"
- **Module id:** `diagnostic`
- **Group:** "Biomarker modeling"

**Module purpose (code fact):** the module lets a user, entirely inside the running app, train and internally test one of six binary classifiers on either a single CpG or a combined CpG panel, using either a bundled sex-stratified whole-blood RA-vs-control panel or their own uploaded methylation matrix + sample sheet. It is a **live-fitting** module (`caret::train()` / `xgboost::xgb.train()` run inside the Shiny session on button click), not a static reproduction of precomputed results — this is the opposite design from the module's own bundled `METHODS_diagnostic_classifier.md`, which documents a fixed six-algorithm sweep already run offline (script09). The module's file-header comment states this design choice explicitly: "Preloaded Data reproduces script09's train/test workflow, not the main pipeline's `beta_raw.rds`/`pheno.rds`."

---

## 2. Web-application implementation

### 2.1 UI structure (exact sub-tab names)

`mod_methyl_diagnostic_ui()` ([`mod_methyl_diagnostic.R:806-823`](../ArthOMix/R/methylomics/mod_methyl_diagnostic.R#L806-L823)) builds one `tabsetPanel` with these tabs, in this exact order:

1. **Datasets**
2. **Feature Source**
3. **Filters & Parameters**
4. **Logistic Regression** *(one tab per classifier — see §2.2)*
5. **Elastic Net**
6. **Support Vector Machine**
7. **Random Forest**
8. **Gradient Boosting / XGBoost**
9. **k-Nearest Neighbors**
10. **Model Comparison**
11. **Test Internal Data**
12. **Export**

The six model tabs (4–9) are generated programmatically from `DXM_MODEL_SPECS` ([`mod_methyl_diagnostic.R:437-530`](../ArthOMix/R/methylomics/mod_methyl_diagnostic.R#L437-L530)), one `tabPanel(spec$label, …)` per entry, so their titles are exactly the `label` fields of that registry (no invented names).

### 2.2 Model registry

`DXM_MODEL_SPECS` (one entry per model tab, `mod_methyl_diagnostic.R:437-530`):

| id | Tab label | `kind` | Underlying engine |
|---|---|---|---|
| `lr` | Logistic Regression | `caret` | `caret::train(method = "glmnet")` with `alpha` fixed (L1 or L2), tuning only `lambda`; or plain `glm` in single-CpG mode |
| `enet` | Elastic Net | `caret` | `caret::train(method = "glmnet")`, tuning both `alpha` and `lambda` (auto search or user grid); or plain `glm` in single-CpG mode |
| `svm` | Support Vector Machine | `caret` | `caret::train(method = "svmLinear"/"svmRadial"/"svmPoly")`, `preProcess = c("center","scale")` |
| `rf` | Random Forest | `caret` | `caret::train(method = "rf")` |
| `gbm` | Gradient Boosting / XGBoost | `xgb` (native) | `xgboost::xgb.cv()` + `xgboost::xgb.train()` directly (not via caret) |
| `knn` | k-Nearest Neighbors | `caret` | `caret::train(method = "knn")`, `preProcess = c("center","scale")` |

These six are the only classifiers this module implements; no others (e.g. neural network) exist in the live app, even though the bundled `METHODS_diagnostic_classifier.md` for the offline script09 pipeline additionally reports an ANN — that is a fact about the offline reference pipeline, not about this Shiny module (see §12).

---

## 3. Features

- Live, in-session model fitting (not lookups of precomputed tables) for all six classifiers.
- Two data-source modes: **Preloaded** sex-stratified whole-blood panel, or **Upload your own dataset**.
- Two feature-panel modes: **Single CpG** or **Combined CpG Panel**.
- Four feature-panel *sources*: Methylomics WGCNA, Methylomics Feature Selection, an uploaded/preloaded panel file, or manual CpG search-and-select.
- Class-imbalance handling: none, class weighting (fold-safe up-sampling), or SMOTE (fold-safe, train-fold-only).
- Configurable cross-validation (folds, repeats, grid vs. randomized search) and four classification-threshold strategies.
- Per-model results: training/CV/internal-test AUC with DeLong or bootstrap CIs, confusion matrices, ROC plots, calibration plots, learning curves.
- Cross-model comparison (ROC overlay, tabular comparison, single-CpG-vs-panel comparison for one model).
- CSV export of metrics and the active feature panel.
- Publishes each tested model's fitted object into a shared `results$diagnostic_models` list for the separate Validation submodule to apply, unrefit, to an external cohort (§9.3).

---

## 4. Inputs

### 4.1 Datasets tab

| Input ID | Control | Choices | Default |
|---|---|---|---|
| `data_mode` | `radioButtons` | "Preloaded whole-blood dataset (sex-stratified RA)" / "Upload your own dataset" | `preloaded` if `METH_DIAG_DATA_AVAILABLE` else `upload` |
| `sex_stratum` (preloaded only) | `radioButtons` | "All samples" / "Female" / "Male" | `female` |
| `upload_matrix` (upload only) | `fileInput` | CSV/TSV, CpG rows × sample columns | — |
| `upload_scale` (upload only) | `radioButtons` | "Beta-value" / "M-value" | `beta` |
| `upload_sheet` (upload only) | `fileInput` | CSV/TSV, one row per sample | — |
| `upload_pheno_col` (upload only) | `selectInput`, auto-populated from sheet columns | any column; auto-guessed via regex `class|phenotype|group|status|disease` | first matching column, else first column |
| `upload_sex_stratum` (upload only) | `radioButtons` | "All samples" / "Female" / "Male" | `all` |
| `upload_sex_col` (upload only, shown only if a sex stratum ≠ all) | `selectInput`, auto-populated; auto-guessed via `sex|gender` | first matching column, else first column |
| `ref_level` | `textInput` | free text | `"Control"` |
| `comp_level` | `textInput` | free text | `"RA"` |
| `train_frac` | `numericInput` | 0.5–0.9, step 0.05 | `0.75` |
| `dxm_seed` | `numericInput` | ≥1 | `42` |
| `validate_btn` | `actionButton` ("Validate Data") | — | — |

### 4.2 Feature Source tab (`req(dxm$validated)`)

| Input ID | Control | Notes |
|---|---|---|
| `feature_source` | `radioButtons` | "Methylomics WGCNA" (`wgcna`) / "Methylomics Feature Selection" (`featureselection`) / "Uploaded/Preloaded Feature Set" (`uploaded`) / "Manual CpG Selection" (`manual`); default `featureselection` |
| `wgcna_module` (preloaded mode) | `selectInput` | "All modules" + distinct module/color values from the loaded assignment table |
| `wgcna_top_n` (preloaded mode) | `numericInput` | 0 = all; default `0` |
| `wgcna_upload` (upload mode) | `fileInput` | WGCNA module-assignment CSV |
| `fs_min_votes` (preloaded mode) | `numericInput` | 1–3, default `2` |
| `fs_top_n` (preloaded mode) | `numericInput` | 0 = all; default `0` |
| `fs_upload` (upload mode) | `fileInput` | Feature Selection module's own RDS export |
| `panel_upload` | `fileInput` | CSV/TXT list of CpG IDs |
| `manual_cpg_select` | `selectizeInput`, multiple | choices = `dxm$all_cpgs` |
| `{wgcna,fs,panel,manual}_load_btn` | `actionButton` | one per source |

### 4.3 Filters & Parameters tab (`req(dxm$validated)`)

| Input ID | Control | Choices | Default |
|---|---|---|---|
| `analysis_type` | `radioButtons` | "Single CpG" (`single`) / "Combined CpG Panel" (`combined`) | `combined` |
| `single_cpg` (single mode) | `selectInput` | `feat$selected` if loaded, else `dxm$all_cpgs` | first choice |
| `imbalance_mode` | `radioButtons` | "None" / "Class weighting (fold-safe up-sampling)" / "SMOTE (fold-safe)" | `none` |
| `cv_folds` | `numericInput` | 3–20 | `10` |
| `cv_repeats` | `numericInput` | 1–10 | `1` |
| `search_method` | `radioButtons` | "Grid Search" / "Randomized Search" | `grid` |
| `threshold_strategy` | `selectInput` | "Default (0.50)" / "Youden's J" / "Sensitivity-focused (≥0.90 sens)" / "Specificity-focused (≥0.90 spec)" | `default` |
| `run_model_choice` | `selectInput` | any of the six model labels | first model (`lr`) |
| `run_selected_model_btn` | `actionButton` ("Run Model") | — | — |

### 4.4 Per-model tabs — parameter inputs (`params_ui`, `mod_methyl_diagnostic.R:437-530`)

| Model | Inputs | Defaults |
|---|---|---|
| Logistic Regression | Penalty (L2/L1 radio); `C` (comma list); max iterations; tolerance | `0` (L2); `"0.01, 0.1, 1, 10"`; `1e5`; `1e-7` |
| Elastic Net | Automatic hyperparameter search (checkbox); if off: `alpha` and `lambda` comma lists; search size (if auto); max iterations; tolerance | `TRUE`; `"0, 0.2, 0.4, 0.6, 0.8, 1"` / `"0.0001, 0.001, 0.01, 0.1, 1"`; `10`; `1e5`; `1e-7` |
| Support Vector Machine | Kernel radio (Linear/RBF/Polynomial); then `C`; sigma (RBF); degree + scale (Polynomial) | `radial`; `"0.25, 0.5, 1, 2, 4"`; `"0.01, 0.05, 0.1"`; `"2, 3"` / `"0.01, 0.1"` |
| Random Forest | Number of trees; `mtry` comma list (blank = 1..p); min terminal node size (0 = default); max terminal nodes (0 = unlimited) | `500`; `""`; `0`; `0` |
| Gradient Boosting / XGBoost | Max boosting rounds; early-stopping rounds; `eta`, `max_depth`, `min_child_weight`, `subsample`, `colsample_bytree`, `gamma` (comma lists) | `200`; `20`; `"0.05, 0.1, 0.3"`; `"2, 3, 4"`; `"1, 3"`; `"0.8, 1"`; `"0.8, 1"`; `"0"` |
| k-Nearest Neighbors | Number of neighbors `k` (comma list) | `"3, 5, 7, 9, 11, 15, 21"` |

Every model tab additionally has a **"Run Model"** / **"Run Single-CpG Diagnostic Analysis"** button, and, once fitted, a **"Run Test Evaluation"**, **"Generate ROC/AUC"**, **"Download ROC plot (PNG)"**, **"Generate Calibration"**, and **"Generate Learning Curve"** button/download.

---

## 5. Processing / data flow

### 5.1 Full pipeline (exact implementation, `dxm_get`-style flow through the module)

```
Preloaded RDS (internal panel, GSE42861 re-Noob'd + granulocyte-adjusted, 21 panel CpGs)
  or Uploaded CSV matrix + sample sheet
        │  filter by sex_stratum + {ref_level, comp_level}
        ▼
  beta-value matrix, samples × CpGs (after transpose)
        │  dxm_beta_to_m(): clip beta to [1e-6, 1-1e-6], M = log2(b/(1-b))
        ▼
  M-value matrix (dxm$all_cpgs = full candidate CpG set)
        │  set.seed(dxm_seed); caret::createDataPartition(y, p = train_frac)
        ▼
  ┌─────────────────────┐        ┌──────────────────────────┐
  │ dxm$train_X/train_y │        │ dxm$test_internal_X/_y   │  (held out, touched only at "Run Test Evaluation")
  └─────────────────────┘        └──────────────────────────┘
        │  Feature Source tab selects CpG id(s) → feat$selected
        │  (WGCNA / Feature-Selection / uploaded panel / manual)
        ▼
  Xtr = dxm$train_X[, ids]      (single CpG or full panel)
        │  dxm_cv_control(): caret::trainControl(method="repeatedcv", ..., sampling = up/SMOTE/NULL)
        ▼
  spec$fit(Xtr, ytr, input, mid, ctrl, seed)   -- one of 6 model specs (§2.2)
        │
        ├─► fitted model (ms$fit)
        ├─► training-set predicted probabilities → ms$train_roc, ms$train_metrics
        ├─► out-of-fold CV predictions (from caret's own `$pred`, or a dedicated
        │      k-fold loop for the native XGBoost path) → ms$cv_roc
        └─► ms$threshold = dxm_pick_threshold(strategy, ms$train_roc)   [fit on train/CV ROC only]
        │
        │  user clicks "Run Test Evaluation"
        ▼
  Xte = dxm$test_internal_X[, ids]  → predict with the already-fitted model, apply ms$threshold
        ▼
  ms$test_internal_roc / ms$test_internal_metrics / ms$confusion_test
        │
        ├─► "ROC/AUC" tab: train, CV, and test ROC curves
        ├─► "Diagnostics": confusion matrices, calibration (from pooled CV predictions), learning curve
        ├─► `runs[[key]]` entry → Model Comparison tab
        └─► results$diagnostic / results$diagnostic_models[[key]] (fitted model + train data + threshold,
              consumed only by the separate Validation submodule, never refit there)
```

### 5.2 Preloaded data path in detail (`observeEvent(input$validate_btn, …)`, preloaded branch, `mod_methyl_diagnostic.R:906-925`)

1. `load_default_diagnostic_train_test()` reads two RDS objects once per R process (in-memory cached): `METH_DIAG_INTERNAL_RDS` (`gse42861_internal_panel_celltype_adjusted.rds`) and `METH_DIAG_EXTERNAL_RDS` (`gse111942_external_panel.rds`). Only the `internal` object is used by this module; `external` is read by the separate Validation submodule.
2. Sex filter: `sex_code <- switch(sex_sel, male="M", female="F", NA)`; `NA` (i.e. "all") keeps every row.
3. Class filter: `internal$pheno$group %in% c(ref_lab, comp_lab)`.
4. Guard: `validate(need(sum(keep) > 20, …))` — fewer than 20 matching samples blocks validation.
5. `y_all <- factor(ifelse(pheno_sub$group == comp_lab, "Class1", "Class0"), levels=c("Class0","Class1"))`.
6. `set.seed(seed); train_idx <- caret::createDataPartition(y_all, p = train_frac, list = FALSE)[, 1]` — a **stratified** split (this is what `createDataPartition` on a factor `y` does).
7. `Xm <- as.data.frame(t(dxm_beta_to_m(beta_sub)))`, `rownames(Xm) <- pheno_sub$gsm`.
8. `dxm$train_X`/`train_y` = rows `train_idx`; `dxm$test_internal_X`/`test_internal_y` = the complement.

### 5.3 Upload data path in detail (same event, upload branch, `mod_methyl_diagnostic.R:926-965`)

Parses the matrix (`methyl_parse_matrix()` — CSV/TSV, first column = probe ID, remaining columns numeric) and the sample sheet (`methyl_parse_sample_sheet()` — any CSV/TSV, one row per sample), matches sample IDs (`methyl_sheet_sample_ids()` — matches a `sample`/`Sample`/`sample_id`/`Sample_ID` column if present, else assumes row order matches, else falls back to row names), requires ≥10 matched samples, requires both `ref_level` and `comp_level` to literally appear in the chosen phenotype column, optionally filters by a user-declared sex column (normalized via `dxm_normalize_sex()`, first letter upper-cased to `F`/`M`), requires ≥6 samples in the smaller class, converts beta→M via `dxm_beta_to_m()` only if `upload_scale == "beta"`, then applies the identical `set.seed(seed); caret::createDataPartition(y_all, p = train_frac)` split as the preloaded path.

### 5.4 Methylation/CpG representation used for fitting

**M-values** (`dxm_beta_to_m()`, `mod_methyl_diagnostic.R:26-29`): `b <- pmin(pmax(beta, 1e-6), 1-1e-6); log2(b/(1-b))`. The module's own header comment is explicit that this is **not** the covariate-adjusted (age/smoking/cell-type-residualized) representation used upstream for CpG selection in other modules — it is raw M-values, matching what the bundled `METHODS_diagnostic_classifier.md` (§2.GG.1) states for the offline pipeline this reproduces, and matching a deployable classifier's realistic access to only raw measured methylation.

### 5.5 Preprocessing performed inside the module

- Beta→M transform with numerical clipping (above).
- `dxm_validate_checklist()` (`mod_methyl_diagnostic.R:74-104`) — a diagnostic report only, computed after the split, never itself altering the data: sample counts, duplicated sample/CpG IDs, missing values, non-numeric columns, constant features, near-zero-variance features (`caret::nearZeroVar`), class-balance ratio, missing phenotype labels, train/test feature overlap, minimum per-class training size.
- Model-specific `preProcess = c("center","scale")` inside `caret::train()` for SVM and kNN only (fit on the training fold(s) by caret's own machinery, never on test data).
- Fold-safe class-imbalance handling: `"weighted"` → caret `sampling = "up"`; `"smote"` → a custom `dxm_smote_fold()` sampling function passed to `caret::trainControl(sampling = list(func = dxm_smote_fold, first = TRUE))`, so resampling happens **inside** each CV fold, never on the held-out data. SMOTE is explicitly disabled for the native-XGBoost path (`dxm_fit_xgb_native()`) with a `validate(need(...))` guard, because that path bypasses caret's `trainControl(sampling=)` hook entirely.

### 5.6 Feature selection performed inside the module

The module performs **no CpG feature selection of its own**. It only lets the user pick which already-computed CpG set to load as the modeling panel, from one of four sources (§4.2):
- **WGCNA** — reads that module's own published per-sex `module_assignment_{sex}[_merged10].csv` table (or an uploaded equivalent), optionally filtered to one module and/or top-N.
- **Feature Selection** — reads that module's own published per-sex `ensemble_votes_{sex}.csv` table (or an uploaded RDS export from that module's "Save Model as RDS"), filtered by minimum vote count (of 3 selection methods) and/or top-N.
- **Uploaded/Preloaded Feature Set** — a plain user-supplied CpG-ID list.
- **Manual CpG Selection** — free-text search/multi-select over `dxm$all_cpgs`.

In every case the selected IDs are intersected with the columns actually present in the active train/test matrices before fitting (`dxm_active_ids()`, `mod_methyl_diagnostic.R:655-660`; `dxm_do_run_model()`, line 623).

---

## 6. Models / classifiers

See §2.2 for the full registry. All six are genuinely implemented and independently runnable; there is no ensemble/voting step across them inside this module (cross-model comparison is visual/tabular only, on the Model Comparison tab — no combined prediction is produced). Hyperparameter tuning for the five caret-based models uses `metric = "ROC"` inside `caret::train()`; the native XGBoost path selects hyperparameters by mean cross-validated AUC read from `xgboost::xgb.cv()$evaluation_log$test_auc_mean` at the early-stopped best iteration (`dxm_fit_xgb_native()`, `mod_methyl_diagnostic.R:165-204`). The code comment explains this departure from `caret::train(method="xgbTree")`: that wrapper "returns NA ROC under xgboost >=3.2's objective-in-params API break."

Single-CpG mode substitutes plain unregularized `glm` for both `lr` and `enet` (`ncol(X) < 2` branch in each spec's `fit`), since `glmnet` requires ≥2 predictor columns.

---

## 7. Validation strategy

### 7.1 Train/internal-test split

- **Mechanism:** `caret::createDataPartition(y_all, p = train_frac, list = FALSE)[, 1]` — stratified by class label (this is `createDataPartition`'s documented behavior on a two-level factor).
- **Fraction:** user-configurable `train_frac` (`numericInput`, range 0.5–0.9, step 0.05), **default 0.75**.
- **Fixed vs. regenerated:** regenerated every time "Validate Data" is clicked, but deterministically — `set.seed(input$dxm_seed %||% 42)` is called immediately before the partition call, so the same seed + same fraction + same filtered sample set always reproduces the same split.
- **Seed:** yes, `dxm_seed` (default `42`), user-editable.

### 7.2 Cross-validation

- Performed via `caret::trainControl(method = "repeatedcv", number = cv_folds, repeats = cv_repeats, classProbs = TRUE, summaryFunction = caret::twoClassSummary, savePredictions = "final", search = grid/random, sampling = ...)` (`dxm_cv_control()`, `mod_methyl_diagnostic.R:128-139`), defaults 10 folds × 1 repeat.
- For the native XGBoost path, CV is instead `xgboost::xgb.cv(..., nfold = cv_folds, stratified = TRUE, early_stopping_rounds = ...)` during hyperparameter search, and a separate manual `caret::createFolds()` loop (`dxm_xgb_cv_roc()`) to reconstruct out-of-fold ROC curves for display, since `xgb.cv` itself does not return per-fold prediction probabilities in the form this module needs.
- **CV occurs only on the training split** — `dxm_do_run_model()` (`mod_methyl_diagnostic.R:622-653`) is only ever called with `Xtr <- dxm$train_X[, ids]`, `ytr <- dxm$train_y`; `dxm$test_internal_X`/`_y` are not referenced anywhere inside model fitting or CV.

### 7.3 Held-out internal test

- `dxm$test_internal_X`/`_y` are produced once at the "Validate Data" step and are **never used for model selection, tuning, or feature-panel choice** — the only place they are read is the per-model "Run Test Evaluation" button handler (`mod_methyl_diagnostic.R:683-730`), which applies the already-fitted model and the already-chosen threshold (fit on train/CV data — §7.4) to them.
- This matches the requirement that "the validation/test set is untouched during model selection": model fitting, hyperparameter tuning, and threshold selection all read only `dxm$train_X`/`train_y` and caret's own CV predictions.

### 7.4 Cutoff / threshold handling

- **User-selectable strategy** (`threshold_strategy` on the Filters & Parameters tab): Default (0.50), Youden's J, Sensitivity-focused (≥0.90 sensitivity), Specificity-focused (≥0.90 specificity).
- **Computed from:** the **training-set** ROC bundle (`ms$train_roc`, built from training-set predicted probabilities) — `dxm_pick_threshold(input$threshold_strategy, ms$train_roc)` (`mod_methyl_diagnostic.R:643`). The code comment for `dxm_pick_threshold()` states explicitly: "Threshold strategies - always computed from train/CV, never from test."
- **Applied to:** both the training confusion matrix and the internal-test confusion matrix (same fixed `ms$threshold` value, never re-picked on test data).
- Sensitivity/specificity strategies fall back to 0.5 if no coordinate on the ROC curve reaches the 0.90 target (avoids `which.max()` silently returning the first coordinate of an all-`FALSE` condition).

### 7.5 Which AUC is reported where

Three distinct AUC values are always shown side-by-side once available, never conflated:
1. **Training AUC** (`ms$train_roc$auc`) — in-sample; labelled "Training AUC" on its own value box.
2. **Mean CV AUC** (`ms$cv_roc$mean_auc` ± SD, from out-of-fold predictions at the best tuned hyperparameters) — labelled "Mean CV AUC (± SD, N folds)".
3. **Internal-test AUC** (`ms$test_internal_metrics$auc`, with a DeLong/bootstrap 95% CI via `pROC::ci.auc`) — labelled "Test AUC (95% CI …)", shown only after "Run Test Evaluation" and computed on the untouched held-out split.

No AUC in this module is computed against an external cohort — the module's own comment on the fitted-model export states: "External-cohort evaluation lives only in that submodule [Validation] - this file only ever produces internal-test results."

---

## 8. Cutoff/threshold implementation

(Cross-reference — see §7.4 for the full mechanism.) The four selectable strategies and their exact formulas, from `dxm_pick_threshold()` (`mod_methyl_diagnostic.R:270-287`):

| Strategy | Formula |
|---|---|
| Default | fixed `0.5` |
| Youden's J | `threshold` at `which.max(sensitivity + specificity - 1)` over `pROC::coords(r, "all", …)` |
| Sensitivity-focused | first threshold (by `which.max` over an ordered condition) where `sensitivity >= 0.9`; else `0.5` |
| Specificity-focused | first threshold where `specificity >= 0.9`; else `0.5` |

---

## 9. Outputs

### 9.1 Plots (`ggplot2`, `dxm_plot_*` functions, `mod_methyl_diagnostic.R:358-427`)

| Plot | Function | Data | Where shown |
|---|---|---|---|
| Training ROC | `dxm_plot_roc()` | `ms$train_roc` | ROC/AUC tab, left panel |
| Cross-validated ROC (per-fold + mean) | `dxm_plot_cv_roc()` | `ms$cv_roc` | ROC/AUC tab, right panel |
| Test ROC | `dxm_plot_roc()` | `ms$test_internal_roc` | ROC/AUC tab, below |
| Calibration curve | `dxm_plot_calibration()` | pooled out-of-fold CV predictions, binned into 10 bins | Diagnostics tab (on demand) |
| Learning curve (train vs. CV AUC across sample fractions 0.4/0.6/0.8/1.0) | `dxm_plot_learning_curve()` | re-fits at each fraction via `dxm_learning_curve()` | Diagnostics tab (on demand) |
| ROC comparison across selected runs | `dxm_plot_roc_compare()` | selected entries of `runs` | Model Comparison tab (on demand) |

### 9.2 Tables

| Table | Content |
|---|---|
| Validation summary | `dxm_validate_checklist()` output — Check / Status (OK/WARN/FAIL) / Detail |
| Selected features table | the loaded CpG panel, with source-specific annotation columns (e.g. WGCNA module, vote counts) |
| Class balance table | Class / N / Percent for the training split |
| Training metrics table | `dxm_metrics_display()` — Accuracy, Balanced accuracy, Sensitivity, Specificity, Precision, NPV, F1, MCC, ROC-AUC, AUC 95% CI, PR-AUC, Brier score, N |
| Test metrics table | same metric set, on internal-test data |
| Training / Test confusion matrices | 2×2 predicted×actual counts |
| Model Comparison table | one row per tested run: Model, Feature set, Train/CV/Test AUC, Threshold, Ran-at timestamp |
| Single CpG vs. Combined Panel table | filtered Model Comparison rows for one chosen model |

### 9.3 Downloadable outputs

| Download | Handler | Content |
|---|---|---|
| ROC plot (PNG), per model | `{mid}_roc_png` | last-rendered ROC ggplot, 7×6in @150dpi |
| Comparison (CSV) | `compare_download` | model/analysis_type/n_features/features/threshold/train_auc/cv_auc/cv_auc_sd/test_auc/ran_at for every run |
| All metrics (CSV) | `export_metrics_csv` | identical column set to the comparison CSV |
| Selected feature panel (CSV) | `export_panel_csv` | the active `feat$table` (or a bare `cpg` column) |

### 9.4 In-memory handoff to the Validation submodule (not a file)

On "Run Test Evaluation", if a shared `results` object was passed to the module, two things are written (`mod_methyl_diagnostic.R:697-728`):
- `results$diagnostic` — a short last-run summary (model label, analysis type, feature count, sex stratum, mode, train/CV/test AUC).
- `results$diagnostic_models[[key]]` — the **full fitted-model artifact**: the `caret`/`xgboost` fit object itself, exact feature ID order, the training-derived threshold, training data/labels, class table, and all training/CV/test metrics — keyed by `paste(model_id, analysis_type, feature_ids)` so every distinct model/panel combination tested in a session is retained, not only the most recent.

This object is read only by `mod_methyl_validation.R` (`avail_models()` / `ref_model()`, lines 388–396) to apply — never refit, never re-tuned, never re-thresholded — against an externally loaded cohort. That downstream application is out of scope for this document.

---

## 10. Sub-tab documentation

### 10.1 Datasets
- **Purpose:** choose and validate the source dataset (preloaded panel or upload), class labels, train/test split fraction, and seed.
- **Inputs:** `data_mode`, `sex_stratum`/`upload_*`, `ref_level`, `comp_level`, `train_frac`, `dxm_seed`.
- **Computation:** loads/filters the data, performs the beta→M transform, the stratified train/test split, and the validation checklist (§5.2–5.5, §7.1).
- **Outputs:** Validation summary box (table + one-line dataset description).
- **Reactive trigger:** `observeEvent(input$validate_btn, …)`.
- **Relevant functions:** `dxm_beta_to_m`, `dxm_validate_checklist`, `methyl_parse_matrix`, `methyl_parse_sample_sheet`, `methyl_sheet_sample_ids`, `dxm_normalize_sex`, `load_default_diagnostic_train_test`.

### 10.2 Feature Source
- **Purpose:** load the CpG set to model, from one of four sources.
- **Inputs:** `feature_source` plus source-specific controls (§4.2).
- **Computation:** reads a published table/RDS/upload, filters/intersects with `dxm$all_cpgs`.
- **Outputs:** Selected features table.
- **Reactive trigger:** one `observeEvent` per source's load button (`wgcna_load_btn`, `fs_load_btn`, `panel_load_btn`, `manual_load_btn`).
- **Relevant functions:** `dxm_load_wgcna_for_sex`, `dxm_load_fs_votes_for_sex`, `load_default_wgcna_module_assignment`, `load_default_diagnostic_ensemble_votes`.

### 10.3 Filters & Parameters
- **Purpose:** choose single-CpG vs. combined-panel analysis, imbalance handling, CV design, threshold strategy, and dispatch a run.
- **Inputs:** §4.3.
- **Computation:** none itself beyond building `dxm_cv_control()`; the "Run Model" button here calls the same `dxm_do_run_model()` as each model tab's own button.
- **Outputs:** Class balance table.
- **Reactive trigger:** `observeEvent(input$run_selected_model_btn, …)`.
- **Relevant functions:** `dxm_cv_control`, `dxm_do_run_model`.

### 10.4 Per-model tabs (Logistic Regression / Elastic Net / Support Vector Machine / Random Forest / Gradient Boosting-XGBoost / k-Nearest Neighbors)
- **Purpose:** configure, run, and inspect one specific classifier.
- **Inputs:** model-specific hyperparameter controls (§4.4) plus the shared "Run Model"/"Run Test Evaluation"/"Generate ROC/AUC"/"Generate Calibration"/"Generate Learning Curve" buttons.
- **Computation:** `spec$fit()` (model training), `dxm_predict_prob`, `dxm_roc_bundle`, `dxm_cv_roc_from_fit`/`dxm_xgb_cv_roc`, `dxm_pick_threshold`, `dxm_metrics_bundle`, `dxm_confusion`, `dxm_calibration`, `dxm_learning_curve`.
- **Outputs:** value boxes (Training AUC, Mean CV AUC, threshold, feature count; then Test AUC/sensitivity/specificity/N), training/test metrics tables, confusion matrices, ROC/calibration/learning-curve plots, ROC PNG download.
- **Reactive trigger:** `observeEvent(input[[paste0(mid,"_run_btn")]], …)`, `..._test_btn`, `..._roc_btn`, `..._calib_btn`, `..._lc_btn`.
- **Relevant functions:** `dxm_register_model_server`, `dxm_render_model_panel`, `dxm_do_run_model`.

### 10.5 Model Comparison
- **Purpose:** compare all tested model/panel combinations within the session.
- **Inputs:** `compare_select` (runs to compare), `compare_curve` (Test/Training/Cross-Validated), `svp_model` (for single-CpG-vs-panel).
- **Computation:** subsets `runs` (a `reactiveValues` accumulator keyed by model+analysis-type+features) and re-plots the chosen ROC bundle.
- **Outputs:** Comparison table (+ CSV download), ROC comparison plot, Single-CpG-vs-Panel table.
- **Reactive trigger:** `observeEvent(input$compare_roc_btn, …)`; the comparison table itself is reactive to `runs` directly.
- **Relevant functions:** `dxm_plot_roc_compare`.
- **Note (code fact, shown to the user in-app):** the module itself warns that Test AUC for WGCNA/Feature-Selection-sourced panels "is drawn from the same cohort the … CpG panels were originally selected on, so it can be optimistically biased for those two feature sources."

### 10.6 Test Internal Data
- **Purpose:** show train/test sample counts and feature overlap diagnostics.
- **Inputs:** none (read-only).
- **Computation:** `dxm_intersect_features(colnames(train_X), colnames(test_internal_X))`.
- **Outputs:** counts of training features, test features, shared features, unmatched (dropped) features.
- **Reactive trigger:** `renderUI` reactive to `dxm$validated`.
- **Relevant functions:** `dxm_intersect_features`.

### 10.7 Export
- **Purpose:** download session results.
- **Inputs:** none beyond the two download buttons.
- **Computation:** serializes `runs`/`feat$table` to CSV.
- **Outputs:** "Download all metrics (CSV)", "Download selected feature panel (CSV)".
- **Reactive trigger:** `downloadHandler` for each button; the tab itself only renders once at least one run exists.
- **Relevant functions:** none beyond base `utils::write.csv`.

---

## 11. Function-by-function code mapping

| Function | Defined | Input | Does | Why used here | Returns | UI/output depending on it |
|---|---|---|---|---|---|---|
| `dxm_beta_to_m` | `mod_methyl_diagnostic.R:26` | beta matrix | clips to [1e-6,1-1e-6], `log2(b/(1-b))` | convert beta→M for fitting, matching script09 | M-value matrix | all model fitting |
| `dxm_parse_num_list` | `:31` | comma text, default | parses numeric list or falls back | lets users type hyperparameter grids | numeric vector | all `params_ui`/`fit` closures |
| `dxm_sex_label` | `:39` | sex code | maps `all/female/male` → display label | consistent labeling | string | Datasets/Feature Source captions |
| `dxm_normalize_sex` | `:43` | free-text sex column | first-letter upper-case → `F`/`M`/`NA` | tolerate arbitrary uploaded sex labels | character vector | upload-mode sex filter |
| `dxm_load_wgcna_for_sex` | `:50` | sex selector | loads/unions published WGCNA table(s) | feature source "wgcna" (preloaded) | data.frame or NULL | Feature Source (wgcna) |
| `dxm_load_fs_votes_for_sex` | `:56` | sex selector | loads/unions/dedups published ensemble-vote table(s) | feature source "featureselection" (preloaded) | data.frame or NULL | Feature Source (featureselection) |
| `dxm_validate_checklist` | `:74` | `dxm` reactive values | 11 data-quality checks | surfaced-to-user QC gate before modeling | data.frame (Check/Status/Detail) | Datasets validation table |
| `dxm_intersect_features` | `:106` | train/test column names | set intersection/difference | shows feature drift between splits | list(train,test,shared,unmatched) | Test Internal Data tab |
| `dxm_smote_fold` | `:116` | fold's X,y | `smotefamily::SMOTE()` inside one fold | fold-safe oversampling | list(x,y) | caret `trainControl(sampling=)` |
| `dxm_cv_control` | `:128` | `input` | builds `caret::trainControl` | shared CV/tuning/imbalance config | `trainControl` object | every model's `fit` |
| `dxm_fit_caret` | `:147` | X,y,grid/length,ctrl,seed,preProcess | wraps `caret::train()` | shared caret entry point | list(model, kind="caret") | LR/ENet/SVM/RF/kNN |
| `dxm_xgb_grid` | `:154` | `input`, model id | builds XGBoost hyperparameter grid | grid search space for native XGBoost | data.frame | `dxm_fit_xgb_native` |
| `dxm_fit_xgb_native` | `:165` | X,y,input,mid,ctrl,seed | manual grid search via `xgb.cv`/`xgb.train` | works around caret+xgboost≥3.2 NA-ROC bug | list(model,kind="xgb",params,nrounds,best) | Gradient Boosting/XGBoost tab |
| `dxm_predict_prob` | `:206` | fit, X | `predict(..., type="prob")` or `xgb` predict | unify prediction interface across engines | numeric vector (P[Class1]) | all metric/ROC computation |
| `dxm_roc_bundle` | `:219` | y, prob | `pROC::roc()` + `ci.auc()` (DeLong, bootstrap fallback) + coords | one shared ROC/AUC/CI object | list(roc,auc,ci_lo,ci_hi,n,coords) | ROC plots, metrics tables, threshold picking |
| `dxm_cv_roc_from_fit` | `:231` | caret `train` fit | reuses caret's own out-of-fold `$pred` at `bestTune` | avoids re-running a second CV loop | list(folds,overall,pooled,mean_auc,sd_auc,n_folds) | CV ROC plot, calibration |
| `dxm_xgb_cv_roc` | `:245` | X,y,params,nrounds,folds,seed | manual `createFolds` + refit per fold | XGBoost has no caret `$pred` equivalent here | same shape as above | CV ROC plot for XGBoost |
| `dxm_pick_threshold` | `:270` | strategy, ROC bundle | one of 4 threshold rules, from train/CV only | user-configurable cutoff, leakage-safe | numeric threshold | training/test confusion matrices |
| `dxm_confusion` | `:289` | y, prob, threshold | 2×2 table + sens/spec/prec/NPV/F1/acc/bal.acc/MCC | standard classification metrics | list incl. `table` | metrics/confusion tables |
| `dxm_metrics_bundle` | `:301` | y, prob, threshold, ROC bundle | adds Brier + `PRROC::pr.curve` PR-AUC to confusion output | full metrics panel | list | metrics tables |
| `dxm_calibration` | `:311` | y, prob, bins | 10-bin mean-predicted vs. mean-observed + `glm` slope/intercept | calibration assessment | list(table,brier,slope,intercept) | Calibration plot |
| `dxm_overfitting_note` | `:328` | train/CV/test AUC | templated overfitting sentence | plain-language interpretation aid | character string | "Results: Test Internal Data" box |
| `dxm_learning_curve` | `:340` | X,y,fit_one,fracs,seed | re-fits at increasing sample fractions | diagnoses sample-size sensitivity | data.frame(frac,n,train_auc,cv_auc) | Learning Curve plot |
| `dxm_plot_roc`, `dxm_plot_cv_roc`, `dxm_plot_calibration`, `dxm_plot_learning_curve`, `dxm_plot_roc_compare` | `:360-427` | respective bundles | `ggplot2` renderers | visualization | `ggplot` object | ROC/AUC, Diagnostics, Model Comparison tabs |
| `dxm_metrics_display` | `:532` | metrics list | formats into a display table | consistent Metric/Value table | data.frame | training/test metrics tables |
| `dxm_render_model_panel` | `:549` | model id/spec/reactive state | assembles the entire per-model tab UI | shared UI builder across 6 tabs | `tagList` | each model tab |
| `dxm_do_run_model` | `:622` | mid,spec,input,dxm,feat,ms | orchestrates one full fit (fit→train ROC→CV ROC→threshold→train metrics) | single entry point used by both the model tab's own button and the Filters & Parameters dispatcher | invisible(TRUE/FALSE) | any "Run Model" button |
| `dxm_active_ids` | `:655` | mid,input,feat | resolves which CpG IDs are "active" (single vs. panel) | shared feature-selection logic | character vector | `dxm_do_run_model` |
| `dxm_register_model_server` | `:662` | full model context | wires all per-model `observeEvent`s/outputs | shared server logic across 6 tabs | (side effect) | every model tab's reactivity |

---

## 12. Leakage / QC audit

All statements below are demonstrated directly by the code cited (§5–§7); nothing here is inferred beyond what is shown.

1. **Split generation:** `caret::createDataPartition(y_all, p = train_frac, list = FALSE)[, 1]` — **stratified** by class label (this is `createDataPartition`'s documented contract for a two-level factor input).
2. **Test fraction:** exactly `train_frac`, a user-editable `numericInput` (range 0.5–0.9, step 0.05), default **0.75**.
3. **Fixed vs. regenerated:** the split is regenerated on every "Validate Data" click, but is fully deterministic given `(seed, train_frac, filtered sample set)` because `set.seed(seed)` precedes the partition call each time.
4. **Seed used:** yes — `dxm_seed`, default `42`, user-editable.
5. **Transformations fitted on training data only:** the only per-CpG transform ("center","scale") is applied inside `caret::train(preProcess=...)` for SVM and kNN, which caret fits on the training folds and applies to held-out folds/test internally — not computed on the pooled train+test matrix beforehand. The beta→M transform (`dxm_beta_to_m`) is a fixed, parameter-free, per-sample function (no fitted statistic), so there is no leakage risk from it regardless of when it is applied.
6. **Scaling parameters frozen:** yes, by construction of `caret::train(preProcess=...)` — center/scale statistics come from the training fold only.
7. **Feature selection relative to the split:** the CpG panel (WGCNA / Feature Selection / uploaded / manual) is chosen **before** the module's own train/test split is used for fitting, but the panel itself was derived **upstream, in separate modules, from the same GSE42861 cohort** the training split is drawn from. The module surfaces this directly to the user on the Model Comparison tab: "Test AUC … is drawn from the same cohort the WGCNA- and Feature-Selection-derived CpG panels were originally selected on, so it can be optimistically biased for those two feature sources — it is internal-test performance only." This is a genuine, code-acknowledged optimism risk for those two feature sources specifically (not for "Uploaded/Preloaded Feature Set" or "Manual CpG Selection", which do not depend on this cohort's own labels for CpG choice within this module). No CpG feature selection is performed by this module itself using any part of the internal-test split.
8. **Cross-validation scope:** `dxm_do_run_model()` builds and tunes exclusively from `dxm$train_X`/`dxm$train_y`; `dxm$test_internal_X`/`_y` are referenced nowhere in fitting/CV/threshold code — confirmed by direct reading of `mod_methyl_diagnostic.R:622-653` and the `_test_btn` handler at `:683-730`, the only place `test_internal_X`/`_y` appear.
9. **Held-out set untouched during model selection:** confirmed — hyperparameter tuning (`metric="ROC"` inside caret, or CV-AUC-based grid search for XGBoost) and threshold selection (`dxm_pick_threshold`, explicitly "never from test") both use only training/CV data; the internal-test set is read only after the user explicitly clicks "Run Test Evaluation" on an already-finalized fit.
10. **User-selectable cutoff:** yes (§7.4/§8), computed once from train/CV data and then applied unchanged to both training and test predictions.
11. **Reported AUC provenance:** the app always labels and displays three separate AUCs — Training (in-sample), Mean CV (out-of-fold on training data), and Test (held-out internal split, with a 95% CI) — never presenting one as a stand-in for another (§7.5).
12. **Bundled `script09_diagnostic_classifier` tables are not read by this module.** Repository-wide search confirms `diagnostic_panel_auc_{sex}.csv` and `diagnostic_perprobe_auc_{sex}.csv` have loader functions defined in `global.R` (`load_default_diagnostic_panel_auc`, `load_default_diagnostic_perprobe_auc`) that are **never called anywhere** in the app's R code, and `diagnostic_perprobe_all_algorithms_{sex}.csv` / `diagnostic_panel_importance_{sex}.csv` have **no loader function at all**. These four files, and `METHODS_diagnostic_classifier.md`, are reference artifacts of the offline research pipeline (script09) documenting a fixed six-algorithm sweep with an external validation cohort (GSE111942) — a design this live Shiny module's *engine* reproduces (same M-value representation, same panel CpGs, same 75/25 stratified split, same seed default), but the live module itself never opens these specific files and has no external-cohort AUC of its own (external evaluation is delegated entirely to the separate Validation submodule). No leakage is demonstrated in the live module by this fact; it only means the "external test AUC" figures in `METHODS_diagnostic_classifier.md` (§`2.GG.3` — e.g. Random Forest external AUC 0.791 for the female panel) describe the offline pipeline, not a number this module computes or displays.
13. **No inflation from the CV/threshold path onto the reported test metrics:** `ms$threshold` is computed once (from training/CV ROC) at "Run Model" time and stored in `ms$fit`'s companion reactive `ms$threshold`; "Run Test Evaluation" reads this same value (`dxm_confusion(yte, ms$test_internal_prob, ms$threshold)`) rather than re-optimizing it against the test labels.

**Overall assessment:** the live module's train/CV/test separation is implemented correctly and leak-safe for model fitting, tuning, and thresholding. The one caveat the code itself flags is feature-panel provenance for the WGCNA and Feature-Selection sources, which draws on the same cohort used for this module's own split — an optimism risk for those two panel sources' internal-test AUC, explicitly disclosed to the user in-app rather than hidden.

---

## 13. Code-to-thesis mapping

| Thesis document element | Source in code |
|---|---|
| Module identity, title, group | `mod_methyl_diagnostic_config` (`mod_methyl_diagnostic.R:15-18`) |
| Registry wiring | `submodules_registry.R:50` |
| Preloaded internal panel provenance (GSE42861, re-Noob'd, granulocyte-adjusted, 21 CpGs, both sexes, 689 samples) | code comments `global.R:653-670`; RDS paths `data_paths.R:109-111` |
| 6-classifier registry | `DXM_MODEL_SPECS`, `mod_methyl_diagnostic.R:437-530` |
| Stratified 75/25 split, seed 42 default | `mod_methyl_diagnostic.R:868, 904, 920` (UI defaults + `createDataPartition` call) |
| Threshold strategies | `dxm_pick_threshold`, `mod_methyl_diagnostic.R:270-287` |
| Fold-safe imbalance handling | `dxm_cv_control`, `dxm_smote_fold`, `mod_methyl_diagnostic.R:116-139` |
| Feature-source cohort-overlap caveat | in-app text, `mod_methyl_diagnostic.R:1185-1186` |
| Handoff to Validation submodule | `mod_methyl_diagnostic.R:703-728`; consumed at `mod_methyl_validation.R:388-396` |
| Unused bundled reference tables | absence of any caller for `load_default_diagnostic_panel_auc`/`load_default_diagnostic_perprobe_auc`, and absence of any loader for `diagnostic_perprobe_all_algorithms_*`/`diagnostic_panel_importance_*` (verified by repository grep, §12.12) |

---

## 14. Short XomicShiny-style implementation description

**Input.** The Diagnostic Classifier submodule accepts either a bundled, sex-stratified whole-blood methylation panel (GSE42861, re-processed via `minfi::preprocessNoob()` and granulocyte-fraction-adjusted upstream, restricted to 21 candidate CpGs) or a user-uploaded methylation matrix (beta- or M-value) with a matching sample sheet, together with user-declared reference/comparison class labels, a training-fraction and random-seed setting, and a chosen CpG feature panel — drawn from the app's own WGCNA module-assignment table, its own Feature Selection ensemble-vote table, an uploaded CpG list, or manual search — for either single-CpG or combined-panel analysis.

**Processing / model.** Beta values are converted to clipped M-values and split into training and internal-test sets via a seeded, class-stratified `caret::createDataPartition()` at a user-set fraction (default 75/25). One of six binary classifiers — logistic regression, elastic-net logistic regression, a support-vector machine, a random forest, gradient-boosted trees (native XGBoost), or k-nearest neighbors — is fit on the training split, with grid- or random-search hyperparameter tuning inside repeated k-fold cross-validation (`caret::trainControl`), optional fold-safe class-imbalance handling (up-sampling or SMOTE, applied only within each training fold), and a user-selectable classification threshold (default 0.50, Youden's J, or a sensitivity-/specificity-targeted cutoff) computed exclusively from training/cross-validated predictions.

**Validation.** Model selection, tuning, and threshold choice never touch the held-out internal-test split; that split is scored only once the user explicitly requests it, using the already-fixed model and threshold. Discrimination is reported separately as in-sample training AUC, mean cross-validated AUC (± SD across folds), and internal-test AUC with a DeLong or bootstrap 95% confidence interval, alongside calibration curves, learning curves, and full confusion-matrix-derived metrics (sensitivity, specificity, precision, NPV, F1, MCC, Brier score, PR-AUC).

**Output.** Results are rendered as per-model ROC/calibration/learning-curve plots and metrics tables, an across-model comparison table and ROC overlay, and CSV exports of metrics and the active feature panel; each tested model, together with its exact feature set and training-derived threshold, is additionally handed off in-memory to this application's separate Validation submodule for evaluation against an independently loaded cohort, without being refit or re-thresholded there.
