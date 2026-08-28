# Cross-Tissue Validation Module — `mod_crosstissue.R`

**Source file:** `ArthOMix/R/transcriptomics/mod_crosstissue.R` (1,480 lines)
**Registration:** `mod_crosstissue_config` (id = `"crosstissue"`, group = "Validation", title = "Cross-Tissue Validation", icon = `"shuffle"`), wired into the app via `submodules_registry.R`.
**Section:** 2.11 (`global.R`: `id = "crosstissue", title = "Cross-Tissue Validation", section = "2.11", satellite = "METHODS_2.11_crosstissue.md"`), positioned after Sex Interaction Analysis (2.10) and before Cross-Ancestry Validation (2.12).
Prepared: 2026-08-25

This document is derived from the code in `mod_crosstissue.R`, cross-checked against `METHODS_2.11_crosstissue.md` (the project's own methods write-up for this section) and the bundled `val_synovium.rds` object it reads. Anything not present in this code is not stated as implemented. Sections marked **Code-confirmed implementation** describe only what the code does; sections marked **Interpretation** provide framing/context and are labeled as such.

---

## 1. Module Purpose (Code-confirmed implementation)

The module header (lines 1–76) states its scope: "'Your analysis' evaluates a user-chosen gene panel in the independent RA synovium dataset (GSE89408, `val_synovium.rds`): sex-stratified discovery (synovium log2FC/significance, direction concordance with blood) plus a full four-classifier panel model (logistic regression, elastic net, random forest, SVM)."

The header is explicit that this is **not** a train/test split like Diagnostic Model (`mod_diagnostic.R`): the synovial cohort has no natural held-out partition of its own — it *is* the held-out compartment relative to blood. Each classifier is therefore fit **once** on the full synovium sex-subset (an apparent/resubstitution AUC, reported as an optimistic upper bound only — Harrell, Lee and Mark, 1996) and separately scored by pooling out-of-fold predictions across an outer cross-validation into one ROC/AUC (the project's own headline synovium estimate).

**What transfers from blood, and what does not** (header, lines 27–33): only the *identity* of the panel genes transfers from blood — never a blood-fitted model's coefficients. Every classifier is refit from scratch within the validation tissue data. This is the same distinction drawn in `METHODS_2.11_crosstissue.md` §2.11.1: "The locked blood model, with its blood-derived coefficients, was not transported into synovium; both the panel-level and the gene-level synovial analyses re-estimate their parameters within the synovial data themselves."

---

## 2. Terminology (Code-confirmed implementation)

| Term | Value |
|---|---|
| Module title | `Cross-Tissue Validation` |
| Analysis description | Sex-stratified cross-tissue validation |
| Validation tissue | Synovium |
| Validation dataset (bundled) | GSE89408 |
| Bundled validation object | `val_synovium.rds` (`VAL_SYNOVIUM_RDS`, `data_paths.R:68`) |
| Key per-gene outputs | Synovium log2FC, statistical significance (adjusted P), direction concordance with blood, sex-stratified results |

This module is distinct from, and must not be confused with:
- **Cross-Ancestry Validation** (Section 2.12) — a separate module.
- **Mendelian Randomization** — a separate module (`mod_coloc.R`/`mod_candidates.R` family).
- **Diagnostic Model** (`mod_diagnostic.R`) — internal train/test evaluation on the *loaded blood dataset*; see §10 below.

---

## 3. Module Structure

### 3.1 UI (`mod_crosstissue_ui`, lines 550–692)

A left sidebar (3 columns) and a main tab area (9 columns).

**Sidebar:**
- **"Validation dataset" box** (new, Section 5/6 of this rewrite): `radioButtons(ns("val_source"), ...)` — "Use preloaded validation dataset (Synovium, GSE89408)" (`preloaded`, default) vs. "Upload my own validation dataset" (`upload`). When `upload` is selected, a `conditionalPanel` reveals `fileInput(ns("val_expr_file"))`, `fileInput(ns("val_meta_file"))`, and `uiOutput(ns("val_column_mapping"))`.
- **"Gene panel & synovium contrast" box**: gene-panel source `radioButtons(ns("panel_source"))` — "Follow this project's pipeline (recommended)" (`project`) vs. "Paste my own gene list" (`own`); a blood-direction status line (`uiOutput("blood_direction_ui")`); saved-run status list (`uiOutput("saved_runs_ui")`).
- **"Advanced filters" box**: significance threshold (`numericInput("sig_cutoff")`, default 0.05); gene-AUC orientation (`radioButtons("orient_view")` — "Best-direction (discovery)" vs. "Train-fixed (cross-dataset)"); outer cross-validation fold count (`numericInput("cv_folds")`, default 10) and fold-stratification choice (`radioButtons("stratified_folds")` — stratified by disease status, default, vs. simple random, matching the project's own offline script).

**Main area — `tabsetPanel(id = "main_tabs")`, four tabs, in this exact order:**

| Tab (exact label from code) | Contents |
|---|---|
| **`Synovium Discovery & Concordance`** | Per-gene synovium discovery — KPI tiles, direction-concordance scatter, gene-AUC lollipop, a per-gene table — inside a `tabsetPanel(id = "disc_sex_tabs")` with **Female**/**Male** sub-tabs, each with its own Run button. |
| **`Panel Classifier - Full Fit`** | The four-classifier apparent (resubstitution) fit — `tabsetPanel(id = "full_sex_tabs")` (Female/Male), each with a `tabsetPanel(..., type = "pills")` over Logistic Regression / Elastic Net / Random Forest / SVM, plus a full-fit comparison table and a right-hand model-parameters box. |
| **`Panel Classifier - Cross-Validated`** | The same four classifiers, scored out-of-fold — `tabsetPanel(id = "cv_sex_tabs")` (Female/Male), each with the same four model pills, plus a cross-validated comparison table. |
| **`Cross-Dataset Comparison`** | Read-only: this module's synovium AUCs lined up against Diagnostic Model's saved blood AUCs (`results$diagnostic`) for the same sex, when available. |

A "References" box (`uiOutput("references_box_ui")`) appears below the tab area once any sex has been run.

These four tab labels are unchanged from the pre-existing code and are **not** renamed by this update — they already describe their implemented functionality precisely (synovium discovery + concordance; apparent panel-classifier fit; cross-validated panel-classifier fit; cross-dataset AUC comparison).

### 3.2 UI-builder function inventory

| Function | Role |
|---|---|
| `mod_crosstissue_params_box(ns, prefix, method_label, defaults_desc, body)` | Generic wrapper box for one classifier's parameter controls |
| `mod_crosstissue_fullfit_panel(ns, prefix, roc_height)` | One model's Full-Fit box (KPI tiles, ROC/CV-fold/tuning plots, performance table + downloads) |
| `mod_crosstissue_cv_panel(ns, prefix, roc_height)` | One model's Cross-Validated box (KPI tiles, pooled ROC plot, performance table + download) |
| `mod_crosstissue_discovery_sex_panel(ns, sex_label)` | One sex's Discovery sub-tab: Run button, KPI tiles, concordance/gene-AUC plots, per-gene table |
| `mod_crosstissue_fullfit_sex_panel(ns, sex_label)` | One sex's Full-Fit sub-tab: Run button, four model pills, comparison table, result line |
| `mod_crosstissue_cv_sex_panel(ns, sex_label)` | One sex's Cross-Validated sub-tab: Run button, four model pills, comparison table |
| `mod_crosstissue_crossdata_sex_panel(ns, sex_label)` | One sex's Cross-Dataset Comparison block: note, bar plot, table |
| `mod_crosstissue_ui(id)` | Top-level UI, described above |

### 3.3 Server function (`mod_crosstissue_server(id, dataset, results)`, lines 694–1480)

Reads the app-wide `results` list (`results$featureselection`, `results$dge_runs`, `results$diagnostic`) and writes `results$crosstissue` once a sex's run completes. Does **not** read the app-wide `dataset` reactiveValues for the validation-tissue expression data — by design, since the whole point of this module is a compartment `dataset` never holds (blood is loaded there; the validation tissue is either the bundled `val_synovium.rds` or a separate upload).

---

## 4. Data Sources (Code-confirmed implementation)

### 4.1 Option A — Preloaded validation dataset

`val_bundled <- { v <- readRDS(VAL_SYNOVIUM_RDS); v$tt <- as.data.frame(v$tt); v }` (once per session). `VAL_SYNOVIUM_RDS` resolves to `data/preloaded/transcriptomics/processed/new/val_synovium.rds` (`data_paths.R:68`).

Verified object shape (`val_synovium.rds`): `logcpm` (19,433 genes × 180 samples, log2-CPM, TMM-normalised), `grp` (factor, levels `Normal` (reference) / `RA` (comparison), 28/152), `sex` (character `"F"`/`"M"`, 120/60), `tt` (per-gene limma-voom result: `logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `B`, `gene` — one table, shared across both sexes, from a sex-adjusted fit), `fsig`/`msig` (bundled female/male consensus gene panels, 6 genes each — the project's own fallback candidate panel).

This bundled workflow is unchanged by this update and continues to work exactly as before.

### 4.2 Option B — User-uploaded validation dataset (added by this update)

Selected via `radioButtons(ns("val_source"))` = `"upload"`. Inputs, matching the smallest structure the existing computation actually needs:

- **Expression matrix** — `fileInput(ns("val_expr_file"))`, CSV or RDS, raw RNA-seq counts, genes in rows / samples in columns (CSV: first column = gene ID). Parsed by `val_expr_raw()` (lines 727–738), reusing the app-wide `tx_parse_expr_matrix_rds()` helper (`global.R:1196`) for RDS and the same `data.table::fread` + first-column-as-rownames pattern Diagnostic Model's External Validation tab uses for CSV.
- **Sample metadata** — `fileInput(ns("val_meta_file"))`, CSV or RDS, parsed by `val_meta_raw()` — same inline RDS/CSV pattern as Diagnostic Model's `ext_meta_raw()`.
- **Column mapping** — `output$val_column_mapping` (auto-populated dropdowns once both files are present): Sample ID column, Sex column, Group column (`val_map_id`/`val_map_sex`/`val_map_group`).
- **Reference/comparison group pick** — `output$val_group_pick_ui`: once a group column is chosen, its distinct values populate `val_ref_group` (e.g. healthy/control) and `val_comp_group` (e.g. disease) selectors — the same "map your own columns" pattern as Diagnostic Model's External Validation tab (`ext_ref_group`/`ext_comp_group`).

No gene-identifier/annotation upload is required beyond the expression matrix's own rownames, and no separate "reference information for comparison with blood" file is requested — blood direction of effect is supplied by `ct_blood_direction()` (§6 below), unchanged, from the app's own live/bundled blood DE, not from the validation-tissue upload.

`val_uploaded()` (reactive) calls `ct_build_uploaded_val(expr, meta, id_col, sex_col, group_col, ref_group, comp_group)` (pure function, lines 407–428), which:
1. Matches expression-matrix sample columns to the metadata's ID column (`validate`: ≥12 matched samples).
2. Derives `sex` as `"F"`/`"M"` from the sex column by prefix match (`^f`/`^m`, case-insensitive) — `validate`s every matched sample resolves.
3. Restricts to the two chosen groups and builds `grp <- factor(..., levels = c(ref_group, comp_group))` (`validate`: ≥12 samples after filtering; ≥4 samples per sex).
4. Calls `ct_voom_de_table(expr, grp, sex)` (pure function, lines 382–405) to compute the per-gene table and log-CPM matrix.

`ct_voom_de_table()` runs the **identical** pipeline `METHODS_2.11_crosstissue.md` §2.11.3 documents for the bundled dataset — `edgeR::DGEList` → `edgeR::filterByExpr(group = grp)` (`validate`: ≥50 genes retained) → `edgeR::calcNormFactors(method = "TMM")` → `stats::model.matrix(~ grp + sex)` (sex-adjusted design, matching "The design matrix contained sex and disease group, so that the disease coefficient is estimated adjusted for sex") → `limma::voom()` → `limma::lmFit()` → `limma::eBayes()` → `limma::topTable(coef = "grp<comparison level>", sort.by = "none")`. `logcpm` is `voom()`'s own `$E` (already log2-CPM on the TMM-normalised library). `fsig`/`msig` are returned as empty (`character(0)`) — an uploaded cohort has no bundled consensus panel of its own; see §5.

This is not a new statistical method: it is the same filterByExpr/TMM/voom/limma/eBayes pipeline already named in this module's own References box (line 1028, "Synovium DE (filterByExpr, TMM, voom, limma, eBayes)") and in `METHODS_2.11_crosstissue.md` §2.11.3, applied live to a user-supplied cohort instead of being read from a pre-computed bundled file.

### 4.3 `val_active()` — single point of dispatch

```r
val_active <- reactive({
  if (identical(input$val_source %||% "preloaded", "upload")) val_uploaded() else val_bundled
})
```

Every downstream function (`ct_project_panel_genes`, `ct_build_sex`, and — through it — `ct_discovery_table`/`ct_fit_sex`) reads `val_active()` rather than a hard-coded bundled object. Because `ct_build_uploaded_val()` returns the exact same field set (`logcpm`, `grp`, `sex`, `tt`, `fsig`, `msig`) as the bundled `val_synovium.rds`, **no downstream analysis function was changed** to add upload support — the uploaded data enters the identical `ct_discovery_table()`/`ct_fit_sex()` code path as the bundled cohort.

---

## 5. Gene-Panel Workflow (Code-confirmed implementation, unchanged)

`ct_project_panel_genes(sex_label)`:
1. Prefers a live `results$featureselection[[sex_label]]$consensus_genes` (this session's own Feature Selection run), if ≥2 genes.
2. Otherwise falls back to `val_active()$fsig` (female) / `$msig` (male) — the bundled project consensus panel. For an uploaded validation cohort this is empty by construction (§4.2), so this fallback only ever fires for the preloaded dataset; the status message explains this explicitly when it fails.

`ct_own_panel_genes(sex_label)`: parses `textAreaInput(ns("gene_list"))` (comma/whitespace/newline-separated), same list for both sexes.

No automatic gene-discovery or new panel-generation algorithm was added. A user validating their own upload without a live Feature Selection panel must paste a gene list — the same mechanism already available for the bundled dataset.

---

## 6. Analysis Workflow (Code-confirmed implementation)

Sequence, per `ct_build_sex(sex_label)` (lines 893–920):

1. **Selected gene panel** → `cand <- ct_project_panel_genes(sex_label)` or `ct_own_panel_genes(sex_label)`.
2. **Validation tissue expression data** → `va <- val_active()`; samples restricted to the chosen sex (`va$sex == sex_code`); `validate`: both groups present with ≥4 samples each.
3. **Sex-stratified gene-level evaluation** → `ct_discovery_table(genes_req, sex_code, va, bd)` (lines 146–167): for each requested gene, looks up `va$tt` (shared, sex-adjusted DE fit) for `syn_log2FC`/`syn_adjP`; looks up blood's `logFC` for the same gene from `ct_blood_direction(sex_label)`; flags `concordant <- sign(syn_log2FC) == sign(blood_log2FC)`; computes `ct_gene_auc()` on this sex's own `logcpm` subset under both the best-direction and blood-train-fixed orientation conventions (§2.11.6 of the methods doc).
4. **Tissue log2FC / statistical evidence** → `syn_log2FC`/`syn_adjP` columns above.
5. **Comparison with blood direction** → `concordant` column above; `ct_biomarker_flag()` (lines 124–128) additionally requires `syn_adjP < sig_cutoff` and `auc_bestdir >= CT_BIOMARKER_AUC_MIN` (0.70) for a gene to count as a "validated cross-tissue biomarker" — the single definition every KPI tile, plot and table column reads from.
6. **Panel-classifier fit** → `ct_fit_sex(expr_sub, y_full, params)` (lines 228–363) fits logistic regression, elastic net, random forest and SVM once on the full sex-subset (apparent AUC) and via `ct_cv_eval()` (lines 180–221) across outer cross-validation folds (pooled out-of-fold AUC + per-fold AUC vector).
7. **Validation results** → returned `fit` list carries `discovery`, per-model apparent/CV results, `dataset_label` (`"synovium, GSE89408"` or `"user-uploaded validation cohort"`), and `grp_levels` (the active reference/comparison labels) for provenance in downstream tables/downloads.

`ct_blood_direction(sex_label)` (unchanged) supplies the blood side of the comparison independently of the validation-tissue data source: live `results$dge_runs` for that sex if this session has already run Differential Expression, else the bundled `dge_results.rds`.

---

## 7. Outputs (Code-confirmed implementation)

### 7.1 Tissue-level effect / statistical evidence / direction comparison / sex-stratified results (Discovery & Concordance tab, per sex)

- KPI tiles (`output$<sex>_disc_stats`): validated cross-tissue biomarker count; panel genes present in the validation dataset (n/N); direction-concordant-with-blood count (n/N); median AUC among biomarkers.
- **Direction concordance plot** (`output$<sex>_concordance_plot`): scatter of blood log2FC (x) vs. synovium log2FC (y), validated biomarkers highlighted/labelled.
- **Gene-AUC plot** (`output$<sex>_geneauc_plot`): ranked lollipop of each present gene's synovium AUC (best-direction or train-fixed, per `input$orient_view`), reference lines at 0.5 (chance) and 0.70 (biomarker cutoff).
- **Per-gene table** (`output$<sex>_disc_table`, `DT::renderDataTable`): gene, cross-tissue-biomarker flag, present, synovium log2FC, synovium adjusted P, blood log2FC, concordant, AUC (best-direction), AUC (train-fixed). CSV download (`output$<sex>_disc_download`).

### 7.2 Panel Classifier — Full Fit (apparent), per sex × model

KPI tiles (apparent AUC, k-fold CV AUC ± SD, selected hyperparameter); apparent ROC plot; CV-AUC-by-fold bar plot; hyperparameter-tuning-grid plot (or a note when not applicable); performance table + CSV download; trained-model `.rds` download (`build_model_bundle()` — bundles the model object, hyperparameters, gene list, and a written-out scoring recipe; `contrast` field now reads `"<comparison> vs <reference> (<dataset_label>)"`, e.g. `"RA vs Normal (synovium, GSE89408)"` or `"RA vs Normal (user-uploaded validation cohort)"`). Per-sex full-fit comparison table across all four models, and a one-line result summary.

### 7.3 Panel Classifier — Cross-Validated (pooled out-of-fold), per sex × model

KPI tiles (pooled out-of-fold AUC ± 95% CI, or "N/A" with a reason; number of samples with an out-of-fold prediction); pooled ROC plot; performance table + CSV download. Per-sex cross-validated comparison table across all four models.

### 7.4 Cross-Dataset Comparison (read-only), per sex

A note (whether `results$diagnostic[[sex]]` — Diagnostic Model's own saved blood AUCs — is available this session, and the panel-gene overlap count); a grouped bar chart (Blood train full-fit / Blood train k-fold CV / Synovium apparent / Synovium pooled CV, per model); a matching table.

### 7.5 References box

Static citation list (filterByExpr/TMM/voom/limma/eBayes provenance; multiple-testing correction; ROC/AUC CIs; elastic net/RF/SVM/caret; the selection-bias caveat on a fixed, blood-derived panel), shown once any sex has been run.

---

## 8. Distinction from Diagnostic Model (Interpretation, code-confirmed)

**Diagnostic Model** (`mod_diagnostic.R`): fits classifiers on the currently loaded blood `dataset`, with a single stratified Train/Test split per sex — an internal, within-cohort held-out evaluation.

**Cross-Tissue Validation** (`mod_crosstissue.R`): evaluates the *same panel identity* (not the blood model's coefficients) in a genuinely separate tissue cohort — bundled synovium (GSE89408) or an uploaded validation-tissue dataset. There is no Train/Test split here because, per the module's own header, this cohort *is* the held-out compartment; every classifier is refit from scratch within it.

These are not interchangeable and this update does not merge them: no Diagnostic Model functionality (Train/Test splitting, `results$diagnostic`'s own write path, Advanced ML Modeling) was added to or duplicated inside `mod_crosstissue.R`.

---

## 9. Distinction About Transfer Performance (Interpretation, code-confirmed)

Because every classifier in `ct_fit_sex()` is refit on the validation-tissue data (§6, step 6) rather than scored with blood-fitted coefficients, the reported AUCs — apparent and pooled cross-validated alike — are **within-dataset re-estimated estimates of the panel's discrimination in synovium**, not a locked-coefficient, out-of-sample transfer score of the blood classifier. `METHODS_2.11_crosstissue.md` §2.11.1 states this directly: "An estimate obtained in this manner cannot be quoted as out-of-sample performance of the blood panel, and is not so quoted here." The Cross-Dataset Comparison tab's own description text makes the same point in-app: "Synovium AUC alongside Diagnostic Model's saved blood AUC — not a transfer of the blood model."

---

## 10. Code-to-Thesis Mapping

| Code element | Thesis description |
|---|---|
| `mod_crosstissue.R` | Cross-Tissue Validation web module |
| `mod_crosstissue_ui` | User interface for tissue validation |
| `mod_crosstissue_server` | Reactive validation workflow |
| `val_bundled` (`VAL_SYNOVIUM_RDS`) | Preloaded validation input (Synovium, GSE89408) |
| `val_uploaded()` / `ct_build_uploaded_val()` | User-provided validation dataset |
| `val_active()` | Single dispatch point between the two data sources, feeding identical downstream code |
| `ct_voom_de_table()` | Sex-adjusted limma-voom DE, computed live for an uploaded cohort using the same pipeline documented for the bundled dataset |
| `panel_source` radio / `ct_project_panel_genes()` / `ct_own_panel_genes()` | Selected genes evaluated in the validation tissue |
| Female/Male Run buttons, `ct_build_sex(sex_label)` | Separate evaluation by sex |
| `va$tt$logFC` (via `ct_discovery_table`) | Tissue-specific effect estimate (synovium log2FC) |
| `va$tt$adj.P.Val` (via `ct_discovery_table`) | Statistical evidence (BH-adjusted significance) |
| `concordant` (via `ct_discovery_table`, vs. `ct_blood_direction()`) | Direction concordance with blood |
| `ct_fit_sex()` apparent fit | Full-fit (apparent, resubstitution) panel-classifier AUC |
| `ct_cv_eval()` pooled result | Cross-validated (out-of-fold) panel-classifier AUC |
| Concordance scatter, gene-AUC lollipop, ROC/CV-fold/tuning plots, Cross-Dataset bar chart | Visual validation results |
| Discovery table, Full-Fit/CV performance tables, comparison tables, Cross-Dataset table | Tabular validation results |

---

## 11. Concise Thesis Subsection (XomicShiny-style, implementation-focused)

> **Cross-Tissue Validation**
>
> *Purpose.* The Cross-Tissue Validation submodule evaluates a user-selected gene panel's identity — never a blood-fitted model's coefficients — in an independent validation-tissue cohort, sex-stratified. It reports each panel gene's tissue-specific effect, its statistical significance, and whether its direction of association agrees with blood, alongside a four-classifier panel model refit entirely within the validation cohort.
>
> *Web-app implementation.* Implemented as `mod_crosstissue.R` ("Cross-Tissue Validation," Validation group, Section 2.11). A left sidebar chooses the validation dataset (bundled Synovium/GSE89408, or an uploaded cohort), the gene panel (this project's live/bundled consensus panel, or a pasted list), and shared filtering/cross-validation settings. The main area has four tabs — Synovium Discovery & Concordance, Panel Classifier - Full Fit, Panel Classifier - Cross-Validated, and Cross-Dataset Comparison — the first three each split into Female/Male sub-tabs with their own Run button.
>
> *Inputs.* Validation-dataset source (preloaded GSE89408, or an uploaded raw-count expression matrix plus sample metadata with sample-ID/sex/group column mapping and reference/comparison group selection); gene-panel source; a Benjamini-Hochberg significance threshold; a gene-AUC orientation convention (best-direction vs. blood-train-fixed); and outer cross-validation fold count/stratification, shared across all four classifiers.
>
> *Processing.* For an uploaded cohort, a sex-adjusted limma-voom differential-expression pipeline (`edgeR::filterByExpr`, TMM normalisation, `voom`, `limma`, `eBayes`) — the same pipeline underlying the bundled dataset — is run once to produce the per-gene effect table and log-CPM matrix that every downstream calculation reads, so uploaded and bundled data enter identical analysis code. Per sex, each panel gene's tissue log2FC, adjusted P-value, and concordance with the corresponding blood log2FC are read from this table; logistic regression, elastic net, random forest and SVM are each fit once on the full validation-tissue sex-subset (apparent AUC) and separately evaluated by outer cross-validation, pooling out-of-fold predictions into a single ROC/AUC per model.
>
> *Outputs.* Per-gene KPI tiles, a direction-concordance scatter, and a ranked gene-AUC plot with a per-gene table (Discovery & Concordance); per-model apparent-fit KPI tiles, ROC/CV-fold/tuning plots, performance tables and downloadable trained-model objects (Full Fit); per-model pooled out-of-fold KPI tiles, ROC plot and performance table (Cross-Validated); and a read-only bar chart/table lining up this module's synovium AUCs against Diagnostic Model's saved blood AUCs (Cross-Dataset Comparison).
>
> *Sub-tabs.*
> - **Synovium Discovery & Concordance** (Female/Male): Input — the selected panel and validation-tissue data; Processing — per-gene tissue log2FC/significance and blood-direction concordance; Output — KPI tiles, concordance/gene-AUC plots, per-gene table.
> - **Panel Classifier - Full Fit** (Female/Male): Input — the same panel/data, all samples; Processing — apparent (resubstitution) fit of four classifiers; Output — KPI tiles, ROC/CV-fold/tuning plots, performance tables, model downloads.
> - **Panel Classifier - Cross-Validated** (Female/Male): Input — the same four classifiers; Processing — outer cross-validation, pooled out-of-fold ROC/AUC; Output — KPI tiles, pooled ROC plot, performance table.
> - **Cross-Dataset Comparison**: Input — this module's synovium AUCs and Diagnostic Model's saved blood AUCs; Processing — none (presentation only); Output — grouped bar chart and table.
>
> *Workflow.* The user chooses a validation dataset (bundled or uploaded) and a gene panel, then runs a sex from any of the first three tabs; the module builds that sex's validation-tissue subset, computes per-gene discovery/concordance, fits and cross-validates the four classifiers, and saves a summary (including which dataset source was used) for reuse elsewhere in the application.

---

## 12. Before → After Summary of Textual UI Renames

No existing UI label was renamed by this update — the four main tab labels, all sub-tab labels, and the sidebar box titles already accurately described their implemented functionality. The only textual changes are **additions** (new controls) and **generalisations** of two static description strings that would otherwise have become misleading once a second data source existed:

| Element | Before | After | Source |
|---|---|---|---|
| Sidebar box (new) | *(did not exist)* | `"Validation dataset"` box with the preloaded-vs-upload radio and upload controls | `mod_crosstissue_ui`, new box before "Gene panel & synovium contrast" |
| Gene-panel box description | `"RA vs Normal synovium (GSE89408), sex-stratified."` | `"RA vs Normal synovium (GSE89408, or your uploaded validation cohort), sex-stratified."` | `mod_crosstissue_ui` |
| Module config `description` | `"...evaluating a chosen gene panel in the independent RA synovial tissue cohort."` | `"...evaluating a chosen gene panel in the independent RA synovial tissue cohort - bundled (GSE89408) or your own uploaded validation cohort."` | `mod_crosstissue_config` |
| `project_source_ui` fallback note | *(no upload-specific case existed)* | Added a case: `"No live Feature Selection panel, and an uploaded validation cohort has no bundled consensus panel - run Feature Selection first, or paste a gene list instead."` | `mod_crosstissue_server`, `output$project_source_ui` |
| Downloaded model bundle `contrast` field | Hard-coded `"RA vs Normal (synovium, GSE89408)"` | Dynamic: `"<comparison> vs <reference> (<dataset_label>)"`, e.g. `"RA vs Normal (user-uploaded validation cohort)"` | `build_model_bundle()` |
| Validation-subset error messages | `"...synovium subset needs at least 4 samples in each group (RA / Normal)."` / `"...present in the synovium dataset."` | Generalised to the active dataset's own group labels / `"validation dataset"` so the message stays accurate for an uploaded cohort with different group names | `ct_build_sex()` |

All four main tabs (`Synovium Discovery & Concordance`, `Panel Classifier - Full Fit`, `Panel Classifier - Cross-Validated`, `Cross-Dataset Comparison`), their Female/Male sub-tab structure, the four-classifier pill tabs, every plot, table, KPI tile, and download button are unchanged.
