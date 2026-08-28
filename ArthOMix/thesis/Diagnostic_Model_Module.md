# Diagnostic Model Module — `mod_diagnostic.R`

**Source file:** `ArthOMix/R/transcriptomics/mod_diagnostic.R` (2,985 lines)
**Registration:** `mod_diagnostic_config` (id = `"diagnostic"`, group = "Biomarker modeling", title = "Diagnostic Model", icon = `"stethoscope"`), wired into the app via `submodules_registry.R`.
Prepared: 2026-08-25

This document is derived **exclusively** from the code in `mod_diagnostic.R`. Anything not present in this code is not stated. Where the module's own in-code comments describe intent or scope decisions, they are quoted or paraphrased and attributed as such, not presented as independent claims.

---

## 1. Module Purpose

The module's own header comment (lines 1–66) states the scope directly:

> "'Your analysis' fits FOUR live classifiers — plain logistic regression, elastic-net logistic regression, random forest and SVM — on a user-chosen gene panel and two-group contrast, sex-stratified exactly like Feature Selection (female and male are ALWAYS fit completely separately, never combined with sex as a covariate)."

Two tabs correspond to two evaluation tiers, per the header comment:

- **"Model Training"** — each sex's own available samples are split **once** (stratified, seed 1234) into Train/Test at a user-chosen ratio (default 70:30). Fitting and tuning happen on Train only; the tab reports Train's full-fit ROC/AUC and a k-fold cross-validated AUC by fold. Test is never touched here.
- **"Model Testing (Internal)"** — the same Train-fit model, applied once to the Test split held out above. "Internal" means held out from the same loaded cohort, not a separate dataset. This tab does not score any separately bundled external validation or cross-tissue holdout used by the project's own offline scripts — that is explicitly out of scope for "a live, point-and-click training tool that only ever has the one dataset currently loaded on the Dataset tab to work with" (header comment).

The header also documents which of the four classifiers are the project's own methodology versus offered as additional choices:

- **Elastic net** and **plain logistic regression** are described as the project's own methodology (referencing the project's offline scripts `15_model_training_elasticnet.R`, `16_model_training_final_panel.R`, `14_model_training_nested_cv.R`, `16b_model_training_final_panel_noMHC.R`, `16d_nested_cv_reconciliation.R`), with the header stating both should be reported, "neither is 'the' winner."
- **Random forest** and **SVM** are stated as not part of the project's own diagnostic-model methodology (only ever used as feature *selectors* elsewhere) — offered here as additional classifier choices on the same gene panel.

The header documents two further scope simplifications relative to the project's own offline scripts: (1) the module does not reproduce the project's own **nested** (feature-selection-redone-per-fold) cross-validation; (2) hyperparameters (alpha/mtry/cost) are tuned once per sex on Train and then reused unchanged inside every k-fold CV refit, rather than re-tuned per fold — stated as "a deliberate simplification to keep a live 'Run Female' / 'Run Male' click fast," making the reported CV number "an upper bound."

A second, independent modeling environment — **"Advanced ML Modeling"** — is layered on top of the four-model engine, gated behind a checkbox unchecked by default ("Additive-only: unchecked by default, so nothing about the four models/panels above changes unless a user explicitly opts in," line 1465). It offers 15 classifiers, configurable feature filtering, hyperparameter tuning, and a nested-CV validation design, and is implemented with entirely separate `diag_adv_*`/`DIAG_ADV_*`-prefixed code that does not read or modify the four-model engine's own reactives (line 630–634, 2530–2536).

---

## 2. Module Structure

### Top-level UI (`mod_diagnostic_ui`, lines 1416–1535)

A left sidebar (3 columns) and a main tab area (9 columns):

**Sidebar — "Gene panel & samples" box:**
- Gene panel source (`radioButtons`, id `panel_source`): "Follow this project's pipeline (recommended)" (`project`), "Paste my own gene list" (`own`), "A WGCNA module from this session" (`wgcna`).
- Contrast controls (`uiOutput("contrast_controls")`, rendered server-side — see §4).
- Saved-run status list (`uiOutput("saved_runs_ui")`).
- "Enable Advanced ML Modeling" checkbox (`adv_ml_enable`, default `FALSE`).

**Main area — `tabsetPanel(id = "main_tabs")` with three tabs:**

| Tab | Contents |
|---|---|
| **Model Training** | (conditionally) the Advanced ML Modeling panel at the top; a shared model-parameters box (`model_params_ui`); a `tabsetPanel(id = "train_sex_tabs")` with **Female / Male / Pooled (all)** sub-tabs, each built by `mod_diagnostic_training_sex_panel()`. |
| **Model Testing (Internal)** | a `tabsetPanel(id = "test_sex_tabs")` with **Female / Male / Pooled (all)** sub-tabs, each built by `mod_diagnostic_testing_sex_panel()`. |
| **External Validation** | `mod_diagnostic_external_panel()` — a single panel, not sex-split. |

A "References" box (`uiOutput("references_box_ui")`) sits below the tab area and appears once any sex has been run.

"Pooled (all)" is documented as pooling every sample regardless of sex "for replicating a published method that never stratified by sex, same convention as Feature Selection's own 'Run All (pooled)'" (line 1499–1502).

### UI-builder function inventory

| Function | Role |
|---|---|
| `mod_diagnostic_training_panel(ns, prefix, title)` | One model's full Training-tab box (KPI tiles, 3 plots, performance table + downloads) — instantiated once per model × sex |
| `mod_diagnostic_testing_panel(ns, prefix, title)` | One model's full Testing-tab box (KPI tiles, ROC plot, performance table + download) |
| `mod_diagnostic_params_box(ns, prefix, method_label, defaults_desc, ...)` | Generic wrapper box for one model's parameter controls |
| `mod_diagnostic_generoc_box_sex(ns, sex_label, mode, title)` | Per-gene ROC/AUC box (plot, hub-gene threshold inputs, table, downloads) |
| `mod_diagnostic_training_sex_panel(ns, sex_label)` | One sex's entire Training sub-tab: Run button, 4 model pills (`tabsetPanel(type="pills")`), the per-gene ROC box, a training comparison table, a result-line summary — hidden via `conditionalPanel` until that sex's Run button (Training or Testing) has been clicked |
| `mod_diagnostic_testing_sex_panel(ns, sex_label)` | One sex's entire Testing sub-tab, same reveal pattern |
| `mod_diagnostic_ui(id)` | Top-level UI, described above |
| `mod_diagnostic_external_panel(ns)` | External Validation tab UI: file uploads, column mapping, panel/group pickers, hub-gene thresholds, Run button |
| `diag_adv_manual_ui_one(ns, model_key)` | One Advanced-ML model's manual-hyperparameter input row, driven by `DIAG_ADV_MANUAL_PARAM_SPECS` |

### Server function (`mod_diagnostic_server(id, dataset, results)`, lines 1571–2985)

Reads the app-wide `dataset` (expression matrix + sample metadata) and the app-wide `results` list (for `results$featureselection` and `results$wgcna`), and writes `results$diagnostic` once a sex's four-model run completes.

---

## 3. Statistical / Computational Engine — Four-Model (Core) Path

### 3.1 Preprocessing shared by every run (`diag_fit_sex`, lines 361–616)

For one sex's full available sample pool:

1. **Split:** `diag_split_train_test()` — `caret::createDataPartition()`, stratified, seed 1234, at the user-chosen train:test ratio (`test_frac`, default 0.3). A single fixed split, not resampled.
2. **Validation guards:** at least 10 train + 4 test samples overall; at least 3 samples per group in Train; both groups present in Test — enforced via `shiny::validate(need(...))`, surfaced as inline error text rather than a crash.
3. **Class weighting** (`diag_class_weight_levels()` / `diag_obs_weights()`): shared by all four models — `"equal"` (default, all weights 1, i.e. the project's own unweighted convention), `"balanced"` (inverse class-frequency), or `"manual"` (a fixed comparison:reference ratio).
4. **Scaling:** Train is per-gene (row) z-scored via `diag_zrows()` (subtract mean, divide by SD, both computed on Train only). Test is standardized using **Train's own** per-gene mean/SD (documented as leakage-free supervised transfer for an in-dataset split). The k-fold CV loop instead re-computes fold-train-only mean/SD at each fold (`diag_cv_auc()`).

### 3.2 The four models, each fit on Train only

- **Elastic net** (`glmnet::cv.glmnet`, family binomial): alpha swept over a user-editable grid (default `0.1, 0.3, 0.5, 0.7, 0.9, 1.0`), each alpha's own lambda chosen by glmnet's internal CV; the alpha minimizing (or, for the AUC metric, maximizing) the chosen CV metric (deviance / AUC / misclassification error) is selected. Final lambda is `lambda.min` (default) or `lambda.1se`.
- **Random forest** (`randomForest::randomForest`): `mtry` auto-tuned via `caret::train(method = "rf")` over a grid derived from the panel size (`1, 2, sqrt(p), p/3, p/2, p`), or set manually; `ntree`, `nodesize`, `maxnodes`, and `classwt` (class weighting) are additional user controls.
- **SVM** (`e1071::svm`): cost auto-tuned via `e1071::tune()` over a user-editable grid (default `0.01…16`), or set manually; kernel choice linear (default, matching the project's own SVM-RFE convention)/radial/polynomial; `scale = FALSE` throughout because the input is already gene-wise z-scored; gamma/degree/tolerance/class-weights all exposed.
- **Logistic regression** (`stats::glm(y ~ ., family = binomial)`): unpenalized, on every gene in the panel, no hyperparameters — described in the header as exactly the classifier the project's own nested-CV and downstream validation scripts score.

For every model: a full-Train-fit ROC/AUC and Youden-optimal cutoff (`pROC::coords(..., best.method = "youden")`); a k-fold CV AUC by fold (`diag_cv_auc()`, default 5 folds per model, each user-adjustable, re-scaling per fold on that fold's own train rows); a one-time score of the Test split at the **locked** Train Youden cutoff (sensitivity/specificity/accuracy — never re-optimized on Test, per the module's own leakage-avoidance comment at lines 286–297); and, where applicable, the actual hyperparameter-search grid and its scores (`tuning_search`), captured for the "explore hyperparameter tuning" plot.

### 3.3 AUC confidence intervals (`diag_auc_ci`)

DeLong's method when n ≥ 20; stratified bootstrap (seed 1234, 2000 reps) when n < 20 — matching, per the code comment, the project's own `auc_ci()` convention (Carpenter & Bithell 2000). A "[separation, n=…]" note (`diag_separation_note`) is appended whenever the CI's lower bound collapses to ≥0.999, flagging a small-n perfect-separation artefact rather than genuine discriminative strength.

### 3.4 Per-gene univariate diagnostics (`diag_gene_roc`, lines 324–342)

For every gene in the panel, independently of the multivariate models: a univariate ROC/AUC on that gene's own raw expression (`pROC::roc(..., direction = "auto")`, so a down-regulated marker's AUC is reported in the direction that is ≥ 0.5) and a two-sided Wilcoxon rank-sum test p-value between the two groups. Computed separately on Train and Test. A user-adjustable "hub gene" rule (default AUC ≥ 0.85 and P < 0.05, matching the module's stated Chen et al. 2021/2022 convention for this panel type) is applied reactively to flag rows without re-running any model.

### 3.5 Overfitting flag (four-model path)

In the Training-tab KPI block (lines 2253–2261): if Train AUC ≥ 0.95 **and** the gap between Train AUC and mean CV AUC is ≥ 0.25, an inline warning is shown ("Likely overfitting… Trust the CV-AUC"), noting logistic regression (no regularization path) as the most exposed of the four techniques. Purely informational; changes no setting.

---

## 4. Server-Side Data Flow (Four-Model Path)

1. `sex_levels()` detects which metadata value in `dataset$meta$sex` means female/male (prefix-matching `"f"`/`"m"`, falling back to sorted distinct values).
2. `contrast_controls` (`renderUI`) exposes reference/comparison group pickers (from `dataset$meta$group`), the Train:Test ratio slider (`test_frac_pct`, 10–50% test, default 30%), and the class-weighting radio/ratio.
3. Gene panel resolution, identical logic for Female/Male/Pooled and reused by External Validation:
   - `project_panel_genes(sex_label)` — prefers a live `results$featureselection[[sex_label]]$consensus_genes`; falls back to a bundled `FS_input_<sex>.csv`.
   - `own_panel_genes(sex_label)` — parses the pasted gene-list textarea, same list for every sex.
   - `wgcna_panel_genes(sex_label)` — reads `results$wgcna$module_genes[[selected module]]`, same module for every sex; requires WGCNA Step 3 (Modules) to have been run first.
4. `diag_build_sex(sex_label, sex_value)` filters `dataset$meta` to the chosen sex (or all samples, if `sex_value = NULL`, for "pooled") and the two chosen groups, intersects candidate genes with the loaded expression matrix's rows, then calls `diag_fit_sex()`.
5. Each sex has its own `eventReactive` (`diag_result_female`/`_male`/`_pooled`), triggered by **two** button IDs each (one on the Training tab, one on the Testing tab) via a shared `reactiveVal` trigger, so a run can be started from either tab.
6. On completion, `save_result()` writes a compact summary (gene count, sample count, each model's full-fit and mean CV AUC, gene list) into `results$diagnostic[[sex_label]]` and fires a `showNotification()`.
7. `active_model_pill` tracks whichever model "pill" (Logistic Regression / Elastic Net / Random Forest / SVM) was last clicked in *either* sex's Training tab, driving a single shared `model_params_ui` parameter box (one box for both sexes, not duplicated).
8. `register_sex_model_outputs(sex_label, res)` is called once per sex ("female", "male", "pooled") and, inside it, once per model (`DIAG_TECHNIQUES`), wiring every KPI tile, plot, table, and download for that sex × model combination, plus the per-sex training/testing comparison tables and the per-gene ROC box (registered once per sex, Train mode only).

---

## 5. External Validation Tab

Purpose, per the module's own comment (lines 1516–1524, 1687–1695): validate the **same gene panel** configured on the left sidebar against a genuinely separate, user-uploaded cohort — explicitly scoped as the same kind of check as "the paper's own Fig. 3d/4a validation step (expression boxplots + per-gene ROC/AUC on the external cohort)," i.e. per-gene, not a refit multivariate model.

**Inputs:** expression matrix upload (CSV/RDS) and sample metadata upload (CSV/RDS); sample-ID and group/diagnosis column mapping (auto-populated dropdowns); reference/comparison group pickers (from the uploaded metadata's mapped group column); a panel-sex selector (Pooled/Female/Male, reusing whichever `panel_source` is set on the left); hub-gene AUC/P thresholds (defaults 0.85 / 0.05).

**Processing:** on "Run external validation" (`eventReactive(input$run_ext_btn)`), samples are restricted to the two chosen groups (each requiring ≥3 samples), the panel's genes are intersected with the uploaded expression matrix's rows (requiring ≥1 present gene), and `diag_gene_roc()` — the identical per-gene AUC/Wilcoxon-P function used for the Train/Test per-gene boxes — is run on the external cohort.

**Outputs:** a status line (genes present, sample counts, panel-source note); a per-gene AUC/P table (CSV-downloadable), with hub-gene rows visually flagged; an expression-by-group boxplot faceted per gene (capped to the top 24 genes by AUC, full table always available regardless of the plot cap).

---

## 6. Advanced ML Modeling (Optional, Additive Environment)

Enabled only via the sidebar checkbox; renders as a box at the top of the Model Training tab (`adv_ml_panel_ui`). Described in-app as: "A second, independent modeling environment: pick any combination of models, a feature-filtering strategy, a hyperparameter-tuning method and a validation strategy, then compare them side by side."

### 6.1 Configuration surface

| Section | Controls |
|---|---|
| Scope | Sample scope (Female/Male/Pooled); feature universe (the same gene panel as the left sidebar, or all genes in the loaded dataset); random seed |
| 1. Feature filtering | One of: no filtering, variance filter, missingness filter, correlation/redundancy filter, univariate statistical filter (raw P or BH-FDR), fold-change/effect-size filter, feature-importance filter, LASSO-based selection, recursive feature elimination — each with its own parameter controls (e.g. top-N cap, correlation cutoff, P/FDR threshold, target feature count) |
| 2. Models | A `checkboxGroupInput` over all 15 registered models (see §6.2); default selection Logistic Regression, Elastic Net, Random Forest, SVM (Linear) |
| 3. Hyperparameter tuning | Manual / Grid search / Random search / Automatic (default); grid size or random-search-iteration count where applicable; a manual-parameter UI generated per selected model from `DIAG_ADV_MANUAL_PARAM_SPECS` |
| 4. Validation strategy | Outer CV folds (performance estimate, default 5); inner CV folds (hyperparameter tuning, default 5); Test-set size (10–50%, default 30%); classification threshold (0.05–0.95, default 0.5) |
| 5. Class imbalance | Resampling on training folds only (none/up/down/SMOTE); class weighting (equal/balanced/manual ratio) |

A pre-flight compute-cost estimate (`diag_adv_estimate_fits()` — `models × (outer_k+1) × (inner_k × tune_length + 1)`) and the current class distribution are both shown before "Compare Models" is clickable. The UI explicitly states LightGBM is unavailable ("no CRAN package available in this R environment") and that Bayesian hyperparameter search is unavailable ("no supported package installed"), with random search offered as the closest alternative.

### 6.2 The 15-model registry (`DIAG_ADV_MODEL_REGISTRY`)

| Group | Models |
|---|---|
| Linear | Logistic Regression (`glm`), Ridge Logistic Regression (`glmnet`, α pinned to 0), LASSO Logistic Regression (`glmnet`, α pinned to 1), Elastic Net Logistic Regression (`glmnet`) |
| SVM | Linear, RBF, Polynomial (caret's `svmLinear`/`svmRadial`/`svmPoly`) |
| Trees | Random Forest (`rf`), Extra Trees (`ranger`, split rule pinned to `extratrees`), Gradient Boosting (`gbm`), XGBoost (custom `caret` modelInfo, see §6.3), AdaBoost (`AdaBoost.M1`), Decision Tree (`rpart`) |
| Other | Naive Bayes, k-Nearest Neighbors |

Each registry entry records how class weighting reaches that model's underlying fitter (`weight_mode`: `native` weights=, `extra_arg` such as `class.weights`/`classwt`, or `unsupported` — e.g. AdaBoost, Naive Bayes, kNN do not support class weighting here) and whether it exposes coefficients vs. an importance measure, driving the drill-down's "Coefficients" vs. "Feature importance" panel choice.

### 6.3 Custom XGBoost integration (`DIAG_ADV_XGB_MODELINFO`)

A hand-written `caret` `modelInfo` list is substituted for caret's own built-in `"xgbTree"`, because — per the code comment — caret's built-in method "errors against xgboost's rewritten ≥2.1 R API (confirmed live in this environment)." It defines its own `grid()` (grid or random search over `nrounds`, `max_depth`, `eta`, `min_child_weight`, `subsample`, `colsample_bytree`, `gamma`, `reg_lambda`, `reg_alpha`), `fit()` (wraps `xgboost::xgboost()`, returning the booster inside a plain list rather than mutating the ALTREP-backed booster object directly, again noted as a live-environment workaround), `predict()`/`prob()`, and `varImp()` via `xgboost::xgb.importance()`.

### 6.4 Feature filtering (`diag_adv_apply_filter`)

Nine methods, split by whether they use the outcome label:

- **Unsupervised** (computed once on the whole Train set, never see `y`): none, variance (rank by variance, optional top-N cap), missingness (drop features above a max-NA-fraction, rank by missingness), correlation (`caret::findCorrelation`, drop one of any pair above a cutoff, default 0.9).
- **Supervised** (re-fit inside every outer CV fold, fold-train only): univariate (Wilcoxon rank-sum P, raw or BH-FDR thresholded), fold-change (mean-difference effect size on z-scored data), importance (random forest `MeanDecreaseGini` ranking), LASSO (`glmnet` α=1, keep nonzero coefficients at `lambda.1se`), RFE (iterative random-forest importance elimination down to a target count, starting-pool capped at 150 features for tractability inside outer CV).

### 6.5 Nested cross-validation pipeline (`diag_adv_run_model`, `diag_adv_compare_models`)

Per the module's own design note (lines 636–648):

- **One fixed Train/Test split** (`diag_split_train_test()`, the same helper as the four-model engine), shared across every selected model.
- **Unsupervised filters** computed once on Train, before the per-model loop (cannot leak class information).
- **Supervised filters, hyperparameter tuning, scaling, and imbalance correction** are all fit inside each outer fold's own training rows only (`caret::createFolds`), with `caret::trainControl(sampling=, preProcess=)` handling the composition of resampling + scaling per inner-CV resample so the inner loop stays leakage-free without duplicating that logic manually. Test is scored **exactly once**, at the end, by one final pipeline fit on the whole of Train — never used to pick a model, filter, or hyperparameter.
- Grid construction (`diag_adv_build_grid`) delegates to each model's own `caret::getModelInfo()$grid()` for "automatic"/"grid"/"random" modes, with small-sample safety patches (e.g. capping `gbm`'s `n.minobsinnode`, `knn`'s `k`, and `AdaBoost`'s `mfinal`/grid size) documented as fixes for errors "confirmed live" in this environment; "manual" mode uses one user-specified hyperparameter row and skips the inner tuning loop entirely.
- Class-imbalance correction (`diag_adv_balance`): up-sampling/down-sampling (`caret::upSample`/`downSample`) or SMOTE (`smotefamily::SMOTE`), applied **after** scaling ("scale-before-balance," because SMOTE's nearest-neighbour interpolation is scale-dependent).
- Per-fold and final-fit metrics (`diag_adv_metrics`): Accuracy, Balanced Accuracy, Sensitivity, Specificity, Precision, F1, MCC, Kappa, ROC-AUC (`pROC`), PR-AUC (`PRROC::pr.curve`), LogLoss (`MLmetrics::LogLoss`), Brier score.

### 6.6 Overfitting/robustness flags (`diag_adv_overfitting_flags`)

Five independent, always-computed, plain-language rule checks: (1) Train-to-Test AUC gap > 0.15; (2) CV AUC vs. Test AUC disagreement > 0.15; (3) fold-AUC SD > 0.12 (unstable CV); (4) retained-feature-to-training-sample ratio > 0.2; (5) Train AUC > 0.995 ("suspiciously close to perfect — check for leakage"). Never hidden; the comparison table's "Overfitting flags" column is simply the count of these that fired for that model.

### 6.7 Results / drill-down

- **Comparison table** (`adv_comparison_df`): one row per selected model — Features retained, mean CV AUC ± SD, Test/Train AUC, Accuracy, Balanced Accuracy, Sensitivity, Specificity, Precision, F1, MCC, Kappa, PR-AUC, LogLoss, Brier, Overfitting-flag count. CSV-downloadable, as is the full analysis configuration (RDS: config, per-model best hyperparameters, per-model retained features).
- **Per-model detail** (`adv_drill_ui`, selected via a model dropdown): KPI tiles (CV AUC, Test AUC, Train AUC, features retained); the overfitting-flag list or a "no warnings" note; ROC (Train/Test/CV overlay), Precision-Recall (Test), and calibration (reliability curve binned into 10 bins, Test predictions vs. observed frequency, with Brier score) plots; a confusion-matrix heatmap and a predicted-probability histogram by group (both Test); a live **threshold explorer** — a slider that recomputes sensitivity/specificity/PPV/NPV/F1/predicted-class counts instantly from the already-stored Test predictions, no refit; a coefficients table (linear models, via `stats::coef`/`glmnet::coef`) or an importance table (tree/SVM models, via `caret::varImp`); the retained-feature list; and downloads for coefficients/importance (CSV), retained features (CSV), Train+Test predictions with probabilities (CSV), and the trained final-fit model object (`.rds`).

**Implementation note (code-fidelity, not user-facing):** the calibration plot's reliability curve is computed directly from the final model's Test-set predictions (`diag_adv_calibration_curve(res$pred_test, res$obs_test, ...)`); the module also defines a Platt-scaling calibrator (`diag_adv_platt_calibrate()`, fit on out-of-fold outer-CV Train predictions) and a standalone threshold-sweep table builder (`diag_adv_threshold_table()`), but neither function is invoked from any registered output — the wired-up "Threshold explorer" instead recomputes metrics directly via `diag_adv_metrics()` at the slider's current value. Both functions are present in the source but currently dead code with respect to the rendered UI.

---

## 7. Output Inventory (Complete)

### Model Training tab, per sex × per model (Logistic Regression / Elastic Net / Random Forest / SVM)

- KPI tiles: Train AUC, k-fold CV AUC (± SD), selected hyperparameter (`diag_hyperparam_value`) — plus an inline overfitting warning when triggered (§3.5).
- ROC plot: Train vs Test overlay with a neutral CV-AUC (mean ± SD) annotation (`diag_roc_plot_traintest`).
- CV-AUC-by-fold bar plot (dashed line at the mean).
- Hyperparameter-tuning-grid plot: alpha vs. CV metric (elastic net), mtry vs. CV ROC (random forest), cost vs. CV error on a log scale (SVM); a plain message for logistic regression (no tuning) or when a manual value bypassed the search.
- Performance table (Train full-fit AUC/threshold/sensitivity/specificity/accuracy; CV mean/SD AUC) + CSV download.
- Trained-model `.rds` download (`build_model_bundle`), bundling the model object, hyperparameters, gene list, group labels, and a written-out scoring recipe (how to z-score new data and call `predict()` for each model type).

Per sex additionally: a per-gene ROC/AUC box (plot faceted per gene, Train+Test curves overlaid, capped at 24 genes by Train AUC; a full table with the hub-gene rule and CSV downloads for all genes and hub-genes-only), a training model-comparison table (all four models' Train/CV AUC side by side), and a one-line run summary.

### Model Testing (Internal) tab, per sex × per model

- KPI tiles: Test-split AUC (n) or "N/A" with a reason if unavailable; number of Test samples held out.
- Test ROC plot (publication-style single curve with AUC/CI).
- Performance table (Test AUC + CI, n/n-positive, sensitivity/specificity/accuracy at the locked training cutoff) + CSV download.

Per sex additionally: a testing model-comparison table (all four models' Test AUC side by side, with the separation-artefact note where applicable).

### External Validation tab

Status line; per-gene AUC/P table (CSV download); expression-by-group boxplot (top 24 genes by AUC).

### Advanced ML Modeling (when enabled)

Comparison table (CSV) + configuration (RDS); per-model drill-down with ROC/PR/calibration/confusion-matrix/probability-distribution plots, a threshold-explorer table, a coefficients-or-importance table (CSV), a retained-features table (CSV), a predictions-and-probabilities CSV, and a trained-model `.rds`.

### References box

A static list of method citations shown once any sex has been run: Hosmer/Lemeshow/Sturdivant (logistic regression), Friedman/Hastie/Tibshirani 2010 and Zou/Hastie 2005 (elastic net/glmnet), Breiman 2001 (random forests), Cortes/Vapnik 1995 (SVM), DeLong 1988 / Carpenter & Bithell 2000 (ROC/AUC CIs), Kuhn 2008 (`caret`).

---

## 8. Workflow Summary

1. **Configure** (sidebar): choose a gene-panel source (project pipeline / pasted list / WGCNA module) and a reference/comparison group contrast; set the Train:Test ratio and class-weighting mode.
2. **Run per sex**: click "Run Female"/"Run Male"/"Run Pooled" from either the Training or Testing tab. `diag_build_sex()` assembles that sex's samples and gene panel; `diag_fit_sex()` splits Train/Test once, fits all four models on Train, computes each model's CV-AUC and per-gene ROC, and scores Test once at the locked Train cutoff.
3. **Inspect**: the Training tab shows Train-fit performance, CV stability, and hyperparameter search per model; the Testing tab shows the once-only held-out Test performance; both expose per-gene univariate diagnostics and a hub-gene rule.
4. **Persist**: each completed sex's summary is written into `results$diagnostic[[sex_label]]` for reuse elsewhere in the app; trained models can be exported as `.rds`.
5. **Optionally validate externally**: upload a separate cohort, map its columns, and check whether the same panel's per-gene AUC/P holds up outside the training data.
6. **Optionally go deeper**: enable Advanced ML Modeling, configure feature filtering / model selection / tuning / validation / imbalance handling, click "Compare Models" to run the nested-CV pipeline across up to 15 classifiers, and drill into any one model's ROC/PR/calibration/confusion-matrix/threshold behavior and interpretability outputs.

---

## 9. Concise Thesis Subsection (as delivered)

See the chat response for the publication-style subsection intended for direct inclusion in the thesis. It is reproduced verbatim below for archival purposes.

> **Diagnostic Model**
>
> *Purpose.* The Diagnostic Model submodule fits and evaluates supervised classifiers that discriminate two user-selected sample groups on a user-chosen gene panel, run separately for female, male, and pooled samples. It reports both a within-sample fit and a genuinely held-out test performance, and offers an independent per-gene univariate diagnostic and an optional expanded modeling environment.
>
> *Web-app implementation.* Implemented as `mod_diagnostic.R` ("Diagnostic Model," Biomarker modeling group). A left sidebar sets the gene panel, group contrast, train:test ratio, and class-weighting; the main area has three tabs — Model Training, Model Testing (Internal), and External Validation — the first two each split into Female/Male/Pooled sub-tabs with their own Run buttons; results are hidden until a run completes. An optional "Advanced ML Modeling" panel (off by default) adds a second, independent 15-model comparison environment inside the Training tab.
>
> *Inputs.* Gene panel source (this project's live/bundled Feature Selection panel, a pasted gene list, or a WGCNA module), reference/comparison group, train:test split ratio (default 70:30), class-weighting mode (equal/balanced/manual ratio), and per-model parameters (cross-validation folds; elastic net alpha grid, lambda choice, and CV metric; random forest ntree/mtry/nodesize/maxnodes; SVM kernel, cost, gamma, degree, tolerance). The External Validation tab additionally takes an uploaded expression matrix and metadata file, column mapping, its own group contrast, and hub-gene AUC/P thresholds.
>
> *Processing.* For each sex, samples are split once (stratified, fixed seed) into Train/Test. On Train, logistic regression (unpenalized `glm`), elastic net (`glmnet`, CV-tuned alpha and lambda), random forest (CV-tuned mtry), and SVM (CV-tuned cost) are each fit on gene-wise z-scored expression, with each model's k-fold cross-validated AUC computed using fold-local rescaling. Test is scored exactly once, standardized with Train's own gene-wise mean/SD, at the Youden cutoff already fixed on Train — never re-optimized on Test. A per-gene univariate ROC/AUC and Wilcoxon rank-sum P are computed independently of the multivariate models, on both Train and Test, with a user-adjustable "hub gene" AUC/P rule. External Validation applies this same per-gene ROC/Wilcoxon procedure to an uploaded, independent cohort rather than refitting a multivariate model. The optional Advanced ML Modeling environment runs a nested cross-validation pipeline — outer folds for performance estimation, inner folds for hyperparameter tuning, both re-fit per fold — over up to 15 classifiers (linear, SVM, tree-based, and other families), with selectable feature filtering, class-imbalance correction, and manual/grid/random/automatic hyperparameter search, reporting a full metric panel (accuracy, sensitivity, specificity, precision, F1, MCC, Kappa, ROC-AUC, PR-AUC, log-loss, Brier score) plus rule-based overfitting flags.
>
> *Outputs.* KPI tiles for Train AUC, cross-validated AUC, and chosen hyperparameter; Train-vs-Test ROC overlays and per-gene ROC panels; a cross-validated-AUC-by-fold plot and a hyperparameter-search plot; downloadable performance tables and trained-model files; per-sex model-comparison tables for both Training and Testing; per-gene AUC/P tables with hub-gene flagging and CSV export; External Validation per-gene tables and expression boxplots. The Advanced ML environment adds a sortable model-comparison table and a per-model drill-down (ROC, precision-recall, calibration, confusion-matrix, and probability-distribution plots; a live threshold explorer; coefficient/importance tables; and full prediction, feature, and model exports).
>
> *Sub-tabs.*
> - **Model Training** (Female/Male/Pooled): Input — the configured panel/contrast/split; Processing — fits all four models on Train and computes their CV-AUC and per-gene diagnostics; Output — KPI tiles, ROC/CV/tuning plots, performance tables, model downloads, and a per-model comparison table.
> - **Model Testing (Internal)** (Female/Male/Pooled): Input — the same Train-fit models; Processing — one-time scoring of the held-out Test split at the fixed training cutoff; Output — Test-AUC KPI tiles, Test ROC plots, performance tables, and a comparison table.
> - **External Validation**: Input — an uploaded, independent expression/metadata cohort and the chosen panel; Processing — per-gene ROC/AUC and Wilcoxon-P on the new cohort; Output — a per-gene results table and expression-by-group boxplots.
>
> *Workflow.* The user selects a gene panel and group contrast, then runs a sex from either the Training or Testing tab; the module builds that sex's Train/Test split, fits and cross-validates the four models on Train, and scores Test once, populating both tabs and saving a summary for reuse elsewhere in the application. The user may separately validate the same panel against an uploaded external cohort, or opt into the Advanced ML environment to compare a wider model set under a configurable nested-cross-validation protocol.
