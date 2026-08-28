# Methylomics Quality Control Module: `mod_methyl_qc.R`

**Source files:**
- `ArthOMix/R/methylomics/mod_methyl_qc.R` (1,619 lines) — UI + server for the "Quality Control" sub-module.
- `ArthOMix/R/methylomics/qc.R` (1,157 lines) — probe/sample QC computation functions, plot builders, batch-correction, and report/export helpers called by `mod_methyl_qc.R`.
- Supporting files read for this document: `ArthOMix/R/methylomics/idat_metrics.R` (54 lines), `ArthOMix/R/methylomics/annotation.R` (102 lines), `ArthOMix/R/methylomics/parse_upload.R` (97 lines), `ArthOMix/R/methylomics/mod_methyl_dataset.R` (298 lines), `ArthOMix/R/methylomics/normalization.R` (boundary check only), `ArthOMix/R/submodules_registry.R`, `ArthOMix/global.R`.

**Registration:** `mod_methyl_qc_config` — `id = "qc"`, `title = "Quality Control"`, `icon = "magnifying-glass-chart"`, `group = "Data"` (`mod_methyl_qc.R:40-43`). Registered first in `MX_MODULES` (`submodules_registry.R:39-40`), the Methylomics sub-module list, invoked as `mod_methyl_qc_server("mx_qc", methyl_dataset, methyl_results)` via `lapply(MX_MODULES, ...)` at `server.R:95`.

Prepared: 2026-08-26.

This document is derived **exclusively** from the code cited above. Every non-trivial technical claim carries a `file:line` citation. Where a claim could not be verified by reading the code, that is stated explicitly rather than inferred. Two label conventions are used throughout to keep general methylomics-QC knowledge visually distinct from what this codebase actually implements:

- **Scientific background:** a statement about methylation-array QC in general (textbook/literature knowledge), not a claim about this code.
- **Code evidence:** a statement about what `mod_methyl_qc.R`/`qc.R` actually do, always with a citation.

---

## 1. Module overview

**Scientific background:** Quality control in DNA methylation array studies (Illumina 450K/EPIC or bisulfite sequencing) typically comprises two layers: *sample-level* QC (call rate, detection p-values, bisulfite-conversion efficiency, intensity, predicted-vs-reported sex, outlier detection via PCA/clustering) and *probe-level* QC (removal of probes with poor detection, low bead count, SNP overlap, cross-reactivity, non-CpG design, sex-chromosome location, or excessive missingness). A QC pass typically precedes normalization and is a prerequisite for any downstream differential-methylation or region-based analysis; well-established Bioconductor tools for this are `minfi`, `ChAMP`, and `wateRmelon`.

**Code evidence:** `mod_methyl_qc.R`'s own header comment states its scope directly: it "reads the shared `methyl_dataset` reactiveValues populated by the Dataset tab ... and runs every probe- and sample-level QC check from `R/methylomics/qc.R` and `R/methylomics/idat_metrics.R` against it. Nothing here removes probes or samples from `methyl_dataset` itself; it reports what each filter WOULD remove and exposes the resulting filtered matrix for download" (`mod_methyl_qc.R:1-10`). This claim was verified independently (see §12 and §14, Finding C-1) by grepping the file for `methyl_dataset$... <-` assignments: none exist anywhere in `mod_methyl_qc.R`, and `methyl_results` (the shared results store for the whole Methylomics group, `server.R:93`) is accepted as a server function parameter (`mod_methyl_qc.R:71`) but never referenced again anywhere in the 1,619-line file — confirmed by an exhaustive grep for `methyl_results` in that file, which returns only the signature line.

**Accepted data formats and structure.** The module reads five main fields off the shared `methyl_dataset` reactiveValues object (initialized at `server.R:82-92`, populated by `mod_methyl_dataset.R`):
- `methyl_dataset$beta` — a numeric probe-by-sample matrix (rows = probes/CpGs, columns = samples), on either the beta (0–1) or M-value scale depending on `methyl_dataset$input_scale` (`"beta"` or `"m"`, set at `mod_methyl_dataset.R:92,208,278`). **Code evidence: row/column orientation is probes-in-rows, samples-in-columns** — confirmed by every filter function in `qc.R` operating with `rowMeans()`/`rowVars()` for probe-level statistics and `colMeans()` for sample-level statistics (e.g. `methyl_filter_missing()` at `qc.R:31` uses `rowMeans(is.na(mat))`; `methyl_sample_call_rate()` at `qc.R:192-194` uses `1 - colMeans(is.na(mat))`), and by the upload parser `methyl_parse_matrix()` which reads "first column the probe ID... everything else numeric" and sets `rownames(m) <- probe_ids` (`parse_upload.R:9-37`).
- `methyl_dataset$sample_sheet` — an optional data.frame of per-sample metadata (any columns; sex/group/batch columns are auto-detected by name pattern).
- `methyl_dataset$rg_set` — an optional `minfi::RGChannelSet` (raw IDAT intensities), present only for an IDAT upload (`mod_methyl_dataset.R:281`) or the preloaded dataset's live matrix path (`NULL` even then — see below).
- `methyl_dataset$mset`, `methyl_dataset$detp`, `methyl_dataset$beadcount` — derived-from-IDAT objects (`minfi::MethylSet`, detection p-value matrix, bead-count matrix), all `NULL` unless an IDAT upload was performed (`mod_methyl_dataset.R:277-284`; `methyl_idat_derive()`, `idat_metrics.R:13-28`).
- `methyl_dataset$array_type`, `methyl_dataset$preloaded`, `methyl_dataset$source` — dataset-provenance metadata used to decide which filters/plots to offer and which explanatory text to show.

Three data-source paths exist on the Dataset tab (`mod_methyl_dataset.R:17-20`): a preloaded whole-blood dataset, an uploaded beta/M-value matrix (CSV/TSV), and uploaded raw IDAT files. **Code evidence:** even in the preloaded path, `methyl_dataset$rg_set`/`mset`/`detp`/`beadcount` are explicitly set to `NULL` (`mod_methyl_dataset.R:95-98`), because "only the derived beta matrix is bundled, not the raw intensities" (`mod_methyl_dataset.R:103`) — so every raw-IDAT-only QC signal (detection p-value, bead count, bisulfite conversion, `minfi::getSex()`, control-probe heatmap) is unavailable for the preloaded dataset even when its live beta matrix is loaded, and only available for a user's own IDAT upload.

**What is configurable.** Every filter/detector exposes its own threshold(s) as a Shiny input with an explicit default (catalogued in full in §9's Parameter inventory table). No threshold is silently fixed except where noted in §14 (Code-level audit).

**What is produced.** Per tab: filtered/flagged tables, diagnostic plots (mostly `plotly`-wrapped `ggplot2`), a probe-retention cascade, a downloadable filtered beta/M-value matrix, a QC summary CSV, and a self-contained HTML/PDF report with baked-in figures (`mod_methyl_qc.R:1543-1595`; `methyl_qc_report_html()`, `qc.R:1099-1134`).

**How QC connects downstream.** **Code evidence:** QC's filtered matrix and flags are *not* automatically consumed by any other module — confirmed by grepping the whole `R/methylomics/` tree and `server.R` for reads of `probe_qc_result`/`sample_qc_result`/etc. or writes into `methyl_results` from this file (none found; see §12). The only connection to the rest of the app is (a) the CSV/HTML/PDF/ZIP downloads a user can manually re-upload elsewhere, and (b) the separate `mod_methyl_normalization` sub-module, which independently reads `methyl_dataset$input_scale` to decide whether it is looking at beta or M-values (`normalization.R:345-349,430-432`) — this is a boundary the QC module does not cross: QC never normalizes, and normalization is documented in this module's own sub-module, not here.

---

## 2. Tab count

**Number of live QC tabs: 8** — Overview, Sample QC, Probe QC, Sex QC, Batch QC, Outlier QC, Visualizations, Reports & Export, enumerated in one `tabsetPanel` at `mod_methyl_qc.R:416-434`.

In addition, a **separate, non-tab "historical pipeline reference" section** sits below the live tabset (`output$default_ui_wrap`/`output$default_ui`, `mod_methyl_qc.R:182-253`), collapsed by default behind an `actionLink` toggle (`mod_methyl_qc.R:185-195`). It reproduces a completed, offline QC run on the preloaded GSE42861 whole-blood dataset by reading three precomputed objects — `load_default_meth_pheno()`, `load_default_meth_qc_sexcheck()` (both `global.R:340-353`, reading CSVs under `METH_DATA_ROOT`), and `METH_QC_PROBE_CASCADE` (a hard-coded `data.frame` literal at `global.R:300-304`) — and does **not** recompute anything from a live matrix. It is documented separately in §3.0 below, not counted among the 8 tabs, per the module's own UI structure and its own header comment: "Kept as-is above the live tool: the historical, exact pipeline numbers plus a genuinely live, recomputable tool underneath it, not one replacing the other" (`mod_methyl_qc.R:62-69`).

---

## 3. Per-tab documentation

### 3.0 Historical pipeline reference (non-tab section)

**Purpose.** Shows the user the *exact*, already-completed QC numbers from the offline pipeline run that produced the preloaded whole-blood dataset — a fixed historical record, not a live tool.

**Visibility gate.** Only rendered when `methyl_dataset$preloaded` is `TRUE` (`mod_methyl_qc.R:183,198`); otherwise shows "No default analysis loaded" (`mod_methyl_qc.R:199-202`).

**Data sources (all precomputed, none recomputed):**
1. **Cohort composition** (`mod_methyl_qc.R:210-215`): a `radioButtons("qc_sex", ...)` stratum picker (Female/Male/Both, `mod_methyl_qc.R:212-213`) filters `load_default_meth_pheno()`'s output into a `group × sex` cross-tab rendered by `output$qc_cohort_table` (`mod_methyl_qc.R:234-240`), via `table(df$group, df$sex)`.
2. **PCA outliers & chrY sex-check** (`mod_methyl_qc.R:216-225`): reads `load_default_meth_qc_sexcheck()` and displays three `valueBox`es — count of `sexcheck$outlier` (">5 MAD, within-sex", `mod_methyl_qc.R:219`), count of `sexcheck$sex_mismatch` (`mod_methyl_qc.R:220`), and total samples checked (`mod_methyl_qc.R:221`) — plus a filtered table (`output$qc_outlier_table`, `mod_methyl_qc.R:242-248`) of only the flagged rows (`df$outlier | df$sex_mismatch`). An explanatory note states outlier detection was run "separately within each sex" because "whole-blood methylation separates strongly by sex on early PCs, so pooling would misflag one sex as 'outliers' relative to the other" (`mod_methyl_qc.R:223`).
3. **Probe-filtering cascade** (`mod_methyl_qc.R:226-230`): renders the hard-coded `METH_QC_PROBE_CASCADE` data frame (`global.R:300-304`) verbatim via `output$qc_cascade_table` (`mod_methyl_qc.R:250-253`), described in its own UI text as "ChAMP-style: cg-prefix restriction, Zhou et al. 2017 MASK_general, multi-hit removal, sex-chromosome removal, >5% missingness" applied as five sequential steps on all 689 samples (`mod_methyl_qc.R:228`). The exact retained/removed counts (`global.R:301-303`): Raw 485,577 → cg-prefix-only 482,421 (−3,156) → MASK_general removal 422,520 (−59,901) → multi-hit removal 422,520 (−0) → sex-chromosome removal 412,492 (−10,028) → missingness>5% removal 412,492 (−0).

**Code evidence this is reproduction, not recomputation:** `global.R`'s own comment on `METH_QC_PROBE_CASCADE` states "the actual logged counts from the completed run, not recomputed here (recomputing it needs the ~2.1GB QC'd beta matrix, excluded when this folder was copied in)" (`global.R:296-299`). `load_default_meth_pheno()`/`load_default_meth_qc_sexcheck()` both simply `data.table::fread()` a CSV under `METH_DATA_ROOT` and return `NULL` if `METH_DATA_AVAILABLE` is `FALSE` or the file is missing (`global.R:340-353`) — neither function performs any statistical computation.

**Interaction with the live tool:** the Sex QC tab's own live check explicitly tells the user, when the sex-chromosome probes have already been stripped by this historical cascade, that "the preloaded whole-blood dataset's bundled matrix has already had sex-chromosome probes removed by the original pipeline's own QC cascade... so there are no chrX/chrY probes left here to predict sex from" (`mod_methyl_qc.R:884-885`) — a direct, code-confirmed dependency of the live Sex QC tab's *degraded* behavior on this historical section's own filtering, even though no data or object is actually shared between them.

---

### 3.1 Tab 1: Overview

**Purpose.** A lightweight, always-visible dataset summary plus a self-contained, independent "basic QC pass" — deliberately *not* an aggregate of the other seven tabs' results. Code comment: "Overview must not silently change just because some other tab was re-run" (`mod_methyl_qc.R:440-443`).

**Input data.** `methyl_dataset$beta` (required, gates the whole tab via `req(methyl_dataset$beta)`, `mod_methyl_qc.R:446`), `methyl_dataset$sample_sheet` (optional — used only to populate the subgroup-column selector and to detect batch/sex columns by regex), `current_subgroup()` (the shared stratum-filtered, exclusion-applied sample set).

**Required/optional columns.** No sample-sheet column is required. If present, `sheet` columns are scanned with `grep("batch|chip|plate|slide|sentrix|array_id|scan_date|^run$", ..., ignore.case = TRUE)` for a batch/chip hint (`mod_methyl_qc.R:448-449`) and intersected against `c("sex","Sex","SEX","gender","Gender")` for a sex-column hint (`mod_methyl_qc.R:450`) — both purely informational text, not filters, on this tab.

**Parameters and defaults.** `live_group_col` (`selectInput`, default = auto-detected sex/gender column if present, else blank = "All samples", `mod_methyl_qc.R:462-463`); `live_stratum` (`radioButtons`, dynamically populated per-level counts, default `"__all__"`, `mod_methyl_qc.R:319-335`). These two inputs are physically defined inside Overview's UI (`output$live_stratum_ui`) but are read by every other tab via `current_subgroup()`, since a `tabsetPanel` renders all tab bodies up front regardless of which tab is selected (`mod_methyl_qc.R:398-406`).

**Reactive dependencies.** `output$overview_ui` is a `renderUI` gated only on `methyl_dataset$beta` (`mod_methyl_qc.R:445-446`) — it re-renders whenever the dataset or `current_subgroup()`'s inputs change, since it inline-calls `sg <- current_subgroup()` (`mod_methyl_qc.R:451`). The actual "Basic QC pass" computation is a separate `eventReactive(input$run_overview_btn, ...)` (`overview_result`, `mod_methyl_qc.R:493-502`), gated behind its own button, independent of the UI's live stratum re-render.

**Functions called.** `methyl_sample_call_rate(mat)` (`qc.R:192-194`, called at `mod_methyl_qc.R:497`); `methyl_qc_status_badge(o)` (`qc.R:929-944`, called at `mod_methyl_qc.R:508`).

**Statistical operations.** Per-sample call rate = `1 - colMeans(is.na(mat))` (`qc.R:193`); overall missingness = `100 * mean(is.na(mat))` (`mod_methyl_qc.R:499`); median call rate across the current stratum (`mod_methyl_qc.R:500`).

**Pass/fail decision logic — `methyl_qc_status_badge()` (`qc.R:929-944`).** This is the *only* place in the entire QC module that computes an explicit pass/warning/fail verdict. Its logic, read in full:
1. Start `status = "pass"`, `reasons = character(0)` (`qc.R:930-931`).
2. If `median_call_rate < 0.90` (not `NA`): `status <- "fail"`, add a reason (`qc.R:932-935`).
3. Else if `overall_missing_pct > 10` (not `NA`) and status is still `"pass"`: `status <- "warning"`, add a reason (`qc.R:937-940`).
4. If no reasons were added, add a generic "No issues detected in this basic pass — run the other QC tabs below for deeper checks." (`qc.R:942`).

Both thresholds (0.90 call rate, 10% missingness) are **hard-coded literals inside `methyl_qc_status_badge()` — not exposed as any Shiny input** (see §14, Finding C-2). Note the `status == "pass"` guard in step 3 means the function can never *downgrade* a fail back to a warning, but a warning can still be issued even when call rate is fine, purely from missingness — the two checks are independent OR-conditions, not a combined score.

**Plots/tables.** No plots. Three `valueBox`es (status, median call rate, overall missingness — `mod_methyl_qc.R:515-519`) and a `tags$ul` of the badge's `reasons` (`mod_methyl_qc.R:520`).

**Output objects.** `overview_result()` — a list with `n_samples`, `n_probes`, `overall_missing_pct`, `median_call_rate`, `subgroup`, `run_at` (`mod_methyl_qc.R:498-501`). Feeds only `current_qc_pieces()$overview` for the Reports & Export summary export (`mod_methyl_qc.R:1442`) — never read by any other of the 8 tabs' own computations.

**Interpretation.** A "Fail" status flags a stratum-wide low-call-rate problem; "Warning" flags moderate global missingness; "Pass" means neither trip-wire fired, explicitly not a guarantee of overall data quality (per the badge's own fallback text).

**Connections to other tabs / downstream.** None computationally — `overview_has_run` only resets other tabs' `*_has_run` flags when a *new dataset* loads (`mod_methyl_qc.R:288-292`), which is UI hygiene, not a data dependency. `current_subgroup()` (defined once, read by every tab) is shared, but that is data scoping, not an Overview-specific output.

---

### 3.2 Tab 2: Sample QC

**Purpose.** Per-sample call-rate, failed-probe-percentage, low-intensity, bisulfite-conversion, and median-intensity diagnostics, plus the manual sample inclusion/exclusion mechanism shared by every other tab.

**Input data.** `current_subgroup()$mat` (the stratum-filtered beta/M-value matrix); `current_rg_subset()` / `current_mset_subset()` (raw-IDAT-derived subsets, `NULL` unless IDAT was uploaded, `mod_methyl_qc.R:347-356`); `methyl_dataset$detp`.

**Required/optional columns.** No sample-sheet column required for call-rate filtering (works "from the beta/M-value matrix alone", `mod_methyl_qc.R:557`). `has_idat <- !is.null(methyl_dataset$rg_set) && !is.null(methyl_dataset$detp)` (`mod_methyl_qc.R:529`) gates whether failed-probe-% and min-intensity filters, bisulfite conversion, and median intensity are offered at all — otherwise an explanatory note is shown (`mod_methyl_qc.R:556-557`).

**Parameters and defaults** (`mod_methyl_qc.R:541-554`):
| Input ID | Type | Default | Notes |
|---|---|---|---|
| `call_rate_min` | numericInput | 0.95 | range 0–1, step 0.01 |
| `sample_detp_thresh` | numericInput | 0.01 | detection-p threshold used for failed-probe % |
| `f_failed_probe_pct` | checkboxInput | FALSE | gates `failed_probe_pct_max` |
| `failed_probe_pct_max` | numericInput | 5 | percent, range 0–100 |
| `f_min_intensity` | checkboxInput | FALSE | gates `min_intensity_thresh` |
| `min_intensity_thresh` | numericInput | 10 | median log2 intensity |

**Reactive dependencies.** `sample_qc_result` — `eventReactive(input$run_sample_qc_btn, ...)` (`mod_methyl_qc.R:574-608`) — reads `current_subgroup()`, `current_rg_subset()`, `current_mset_subset()`, `methyl_dataset$detp`, and the six inputs above, all captured at button-click time only.

**Functions called (with citations).** `methyl_sample_call_rate(mat)` (`qc.R:192-194`, called `mod_methyl_qc.R:580`); `methyl_bisulfite_conversion(rg)` (`idat_metrics.R:33-40`, wraps `wateRmelon::bscon()`, called `mod_methyl_qc.R:586`); `methyl_median_intensity(mset)` (`idat_metrics.R:45-53`, wraps `minfi::getQC()`, called `mod_methyl_qc.R:587`); `methyl_sample_failed_probe_pct(mat, detp, threshold)` (`qc.R:202-217`, called `mod_methyl_qc.R:590`); `methyl_sample_low_intensity(median_int, min_intensity)` (`qc.R:225-232`, called `mod_methyl_qc.R:595`).

**Statistical operations.** Call rate: `1 - colMeans(is.na(mat))`. Failed-probe %: `colMeans(detp[common_probes, common_samples] > threshold, na.rm = TRUE) * 100` restricted to probes/samples common to `mat` and `detp` (`qc.R:214-215`). Low-intensity score: `(med_meth_log2 + med_unmeth_log2) / 2`, flagged when `< min_intensity` (`qc.R:230-231`) — mirrors `minfi`'s own `mMed`/`uMed` QC-tutorial convention per the function's comment (`qc.R:219-224`). Bisulfite conversion: `wateRmelon::bscon(rg_set)` percentage, no further transformation.

**Plots/tables.** `sample_qc_table` (DT, call rate + optional flags, `mod_methyl_qc.R:640-645`); `bisulfite_table` (DT, IDAT-only, `mod_methyl_qc.R:647-653`); `median_int_table` (DT, IDAT-only, `mod_methyl_qc.R:655-660`). No plots on this tab.

**Output objects.** `sample_qc_result()` — list with `mat`, `subgroup`, `sample_qc` (data.frame), `bisulfite`, `median_int`, `settings`, `run_at` (`mod_methyl_qc.R:602-607`).

**Interpretation.** A sample below `call_rate_min` or above `failed_probe_pct_max`, or below `min_intensity_thresh`, is flagged in its own boolean column but **never removed automatically** — flags are display-only on this tab (manual exclusion is a separate, explicit user action below).

**Manual inclusion/exclusion mechanism** (physically part of this tab's UI, `mod_methyl_qc.R:560-569`, but consumed globally): `stratum_all_samples()` (`mod_methyl_qc.R:363-367`) lists every sample in the current stratum *before* exclusion, so an excluded sample stays selectable for re-inclusion. `manual_table` (DT, multi-row-select) plus `manual_apply_btn`/`manual_clear_btn` write to the shared `manual_exclude` `reactiveVal` (`mod_methyl_qc.R:280,379-385`), which `current_subgroup()` (via `methyl_apply_manual_exclude()`, `qc.R:507-519`) applies for *every* subsequent method's *next* run.

**Connections to other tabs.** The manual-exclusion table/buttons are shared UI wired into `current_subgroup()`, so excluding a sample here changes the sample set every other tab reads the next time *that tab's own* Run button is clicked — not retroactively for already-computed results. Sex QC's discordant-sample "Exclude" action (`mod_methyl_qc.R:929-935`) and Outlier QC's "Apply Sample Exclusions" (`mod_methyl_qc.R:1155-1160`) write to the identical `manual_exclude` reactiveVal, so all three exclusion mechanisms converge on one shared list, confirmed by all three calling `manual_exclude(union(manual_exclude(), sel))`.

---

### 3.3 Tab 3: Probe QC

**Purpose.** The module's central probe-filtering tool: up to eleven independent, composable filters, each returning `list(keep=<logical>, note=<string>)`, combined into one AND'd keep vector and a sequential retention cascade.

**Input data.** `current_subgroup()$mat`; `anno_result()` (a `reactive` wrapping `methyl_get_annotation(methyl_dataset$array_type)`, `mod_methyl_qc.R:260-263`); `methyl_dataset$detp`, `methyl_dataset$beadcount`; optional user-uploaded cross-reactive-probe list and MAF table.

**Required/optional columns.** No sample-sheet columns required; probe-level filters operate on `rownames(mat)` alone or against Bioconductor manifest annotation.

**Parameters and defaults** (`mod_methyl_qc.R:669-714`) — full detail in §9's parameter table; summarized here:
| Filter | Input(s) | Default enabled | Default threshold |
|---|---|---|---|
| Detection p-value | `f_detp`/`detp_thresh` | TRUE iff IDAT present | 0.01 |
| Bead count | `f_beadcount`/`beadcount_thresh` | TRUE iff IDAT + beadcount present | 3 |
| Missing values | `f_missing`/`missing_max` | TRUE | 0 (i.e. zero tolerance) |
| SNP-overlap | `f_snp` | TRUE (Illumina arrays only) | n/a |
| Non-CpG (CpH) | `f_noncpg` | TRUE (Illumina arrays only) | n/a |
| Sex-chromosome | `sexchr_mode` | "Keep all" | keep / remove_xy / remove_y_only |
| Cross-reactive (uploaded) | `f_crossreactive`/file | FALSE | n/a, needs upload |
| MAF (uploaded) | `f_maf`/file/`maf_max` | FALSE | 0.05, needs upload |
| Variance | `f_variance`/`variance_min` | FALSE | 0 |
| Standard deviation | `f_sd`/`sd_min` | FALSE | 0 |
| Mean-range | `f_meanrange`/`mean_lo`,`mean_hi` | FALSE | 0.01–0.99 (beta) or −6–6 (M) |

**Reactive dependencies.** `probe_qc_result` — `eventReactive(input$run_probe_qc_btn, ...)` (`mod_methyl_qc.R:720-777`) — requires `length(sg$included) >= 3` (`mod_methyl_qc.R:722-723`), the only tab with an explicit minimum-sample validation gate on Probe QC itself.

**Functions called (with citations), applied conditionally in this literal order** (`mod_methyl_qc.R:741-752`): `methyl_filter_detection_p()` (`qc.R:162-172`) → `methyl_filter_beadcount()` (`qc.R:178-188`) → `methyl_filter_snp()` (`qc.R:68-81`) → `methyl_filter_non_cpg()` (`qc.R:63-66`) → `methyl_filter_sex_chr()` (`qc.R:89-104`) → `methyl_filter_cross_reactive()` (`qc.R:112-119`) → `methyl_filter_maf()` (`qc.R:149-160`) → `methyl_filter_missing()` (`qc.R:30-34`) → `methyl_filter_variance()` (`qc.R:36-41`) → `methyl_filter_sd()` (`qc.R:47-52`) → `methyl_filter_mean_range()` (`qc.R:54-58`). Combined: `keep <- rep(TRUE, nrow(mat)); for (f in filters) keep <- keep & f$keep` (`mod_methyl_qc.R:754-755`) — this is an unordered AND across all *enabled* filters (not order-sensitive for the final `filtered` matrix), separately from the order-sensitive cascade below. Also: `methyl_probe_retention_cascade(nrow(mat), filters)` (`qc.R:547-560`, called `mod_methyl_qc.R:759`) and, for upload parsing, `methyl_parse_probe_list()` (`parse_upload.R:50-57`) and `methyl_parse_maf_list()` (`qc.R:131-144`).

**Statistical operations.** Detection-p filter: probe removed if it fails `p > threshold` in **any** sample (`qc.R:169-171`, `frac_fail == 0` required to keep). Bead-count filter: probe removed if bead count `< threshold` in **any** sample (`qc.R:185-187`) — see §14 Finding I-1 for the deliberate ChAMP-default divergence. SNP filter: probe removed if any of `Probe_rs`/`CpG_rs`/`SBE_rs` manifest columns is non-NA/non-empty (`qc.R:74-78`). Non-CpG filter: `grepl("^cg", probe_ids, ignore.case=TRUE)` via `methyl_probe_is_cpg()` (`annotation.R:99-101`). Sex-chromosome filter: manifest `chr` column matched against `chrX`/`chrY`/`X`/`Y` depending on mode (`qc.R:97-103`). Variance/SD/mean-range: `methyl_row_vars()` (matrixStats-accelerated, `qc.R:22-28`) thresholded directly.

**Plots/tables.** `filter_table` (DT: filter name, probes removed, note — `mod_methyl_qc.R:807-816`); `cascade_plot` (plotly-wrapped `methyl_plot_cascade()`, a bar chart of probes remaining after each sequential filter, `qc.R:818-826`, called `mod_methyl_qc.R:818-820`).

**Output objects.** `probe_qc_result()` — list with `mat` (pre-filter), `filters` (named list of filter results), `keep` (combined logical), `filtered` (post-filter matrix), `cascade` (data.frame), `subgroup`, `anno`, `settings`, `run_at` (`mod_methyl_qc.R:775-776`). This is the *only* tab whose result (`filtered`) other tabs (Visualizations) read.

**Interpretation.** The valueBox row (Probes in / Probes kept / Probes removed / Samples, `mod_methyl_qc.R:788-791`) and the cascade plot together show how much of the panel survives; the filter table breaks that down per individual filter.

**Connections to other tabs.** **Code evidence — the one genuine cross-tab data dependency in the whole module:** the Visualizations tab's eight probe-QC-derived plots (PCA 2D/3D, MDS, density, boxplot, violin, correlation heatmap, mean-SD) all read `probe_qc_result()$filtered` directly (e.g. `mod_methyl_qc.R:1305,1330,1341-1348,1373,1382`), gated on `probe_qc_has_run()` (`mod_methyl_qc.R:1248`) — this is the single exception to the module's otherwise fully independent per-tab design (see §11).

**Connection downstream.** `dl_filtered_beta`/`dl_filtered_mvalue` (`mod_methyl_qc.R:1492-1507`) and `dl_filter_summary` (`mod_methyl_qc.R:1508-1518`) let a user download the filtered matrix for re-use outside the app; nothing inside ArthOMix consumes it automatically.

---

### 3.4 Tab 4: Sex QC

**Purpose.** Predicts each sample's genetic/methylation-inferred sex and compares it against the sample sheet's reported sex, flagging discordant samples for optional manual exclusion.

**Input data.** `current_subgroup()$mat`; `anno_result()`; `current_rg_subset()` (raw IDAT, optional); `methyl_dataset$sample_sheet`'s sex/gender column, if present.

**Required/optional columns.** Sex/gender column optional — intersected against `c("sex","Sex","SEX","gender","Gender")` (`mod_methyl_qc.R:857`); if absent, the tab still predicts sex but shows "nothing to compare it against" (`mod_methyl_qc.R:879-880`).

**Parameters.** None exposed as tunable inputs — a single `run_sex_qc_btn` (`mod_methyl_qc.R:841`).

**Reactive dependencies.** `sex_qc_result` — `eventReactive(input$run_sex_qc_btn, ...)` (`mod_methyl_qc.R:847-863`) — requires `length(sg$included) >= 1` (`mod_methyl_qc.R:849`).

**Functions called.** `methyl_sex_check(mat, anno, rg, reported_sex)` (`qc.R:407-446`, called `mod_methyl_qc.R:861`), which internally dispatches to `methyl_cluster_sex()` (`qc.R:356-376`) and `methyl_sex_check_attach_mismatch()` (`qc.R:383-392`); `methyl_sheet_sample_ids()` (`qc.R:463-468`, resolves sheet rows to matrix column IDs).

**Statistical operations — two mutually exclusive code paths inside `methyl_sex_check()`:**
1. **Raw-IDAT path** (preferred, used when `rg_set` is non-`NULL` and `minfi` is installed): `minfi::mapToGenome(minfi::preprocessRaw(rg_set))` then `minfi::getSex(gmset)` (`qc.R:409-411`) — a copy-number-based method using median chrX/chrY total-intensity log2 ratios, described in the function's own comment as "the standard method most published pipelines use" (`qc.R:394-396`).
2. **Beta-heuristic fallback** (used otherwise, e.g. for the preloaded dataset or a matrix-only upload with manifest annotation available): mean chrX and chrY beta per sample computed from manifest-annotated probes (`qc.R:433-434`), then `methyl_cluster_sex(mean_y, reported_sex)` clusters `mean_y` into two groups via **k-means (`stats::kmeans(y, centers=2, nstart=25)`)** when ≥6 samples and ≥2 unique values exist, else falls back to a global median split (`qc.R:358-362`). The higher-mean cluster is labeled "M" *unless* a reported-sex vector is supplied, in which case each cluster is relabeled by **majority concordance** with reported sex instead of assuming direction (`qc.R:365-374`).

**Code-evidence detail on the k-means-over-median-split design choice:** the function's own comment (`qc.R:340-355`) states a median split "misclassifies the *minority* cluster's mislabeled tail whenever the cohort is imbalanced enough that the global median falls inside the majority cluster instead of the gap between the two clusters — exactly what happens on this app's own reference cohort (492 female / 197 male: a median split alone misclassifies 147 of 492 females as male, despite chrY beta being cleanly bimodal with zero overlap between the true clusters)." This is presented as a fix that "mirrors" the same correction already applied in the reference pipeline's own `sex_km`/`sex_pred` columns (`qc.R:353-355`) — i.e. this is a documented, code-cited real-data failure mode this implementation deliberately avoids, not a generic claim.

**Plots/tables.** `sex_table` (DT of predicted sex, `mod_methyl_qc.R:901-904`); `sex_scatter` (plotly, X-methylation vs. Y-methylation scatter colored by predicted sex, shaped by concordance, `mod_methyl_qc.R:906-919`); `discordant_table` (DT, only mismatched rows, selectable, `mod_methyl_qc.R:922-927`).

**Output objects.** `sex_qc_result()` — list `sex` (from `methyl_sex_check()`: `ok`, `method`, `detail` data.frame, `n_mismatch`), `subgroup`, `run_at` (`mod_methyl_qc.R:862`).

**Interpretation.** `n_mismatch > 0` flags possible sample mislabeling. The tab distinguishes "ran but not comparable" (no sex column, `n_mismatch = NA`) from "ran and concordant" (`n_mismatch = 0`) from "ran and discordant" (`n_mismatch > 0`) (`mod_methyl_qc.R:875-880`).

**Connections to other tabs.** The "Exclude selected discordant samples" button (`mod_methyl_qc.R:896,929-935`) writes into the same shared `manual_exclude` reactiveVal as Sample QC's manual table and Outlier QC's exclusion button — the only cross-tab interaction, and it is data-scoping (sample exclusion), not a shared analysis result.

---

### 3.5 Tab 5: Batch QC

**Purpose.** PCA-visualized batch/unwanted-variation correction via two alternative methods, ComBat (needs a labeled batch column) or RUVm (needs raw IDAT control probes, no batch label required).

**Input data.** `current_subgroup()$mat`; `methyl_dataset$sample_sheet` (batch/chip column, ComBat); `current_rg_subset()` (control probes, RUVm).

**Required/optional columns.** ComBat requires a detected batch column (regex `batch|chip|plate|slide|sentrix|array_id|scan_date|^run$` against sheet column names, `mod_methyl_qc.R:944-945`); if none exists and RUVm is also unavailable, the whole tab shows only an explanatory card (`mod_methyl_qc.R:947-953`). RUVm requires a user-chosen "factor of interest" column (any sheet column, `mod_methyl_qc.R:964`) and `ruvm_available()` (`rg_set` non-`NULL` and `missMethyl` installed, `mod_methyl_qc.R:394-396`).

**Parameters and defaults.** `batch_method` (`selectInput`, "ComBat" default if a batch column exists else "RUVm", `mod_methyl_qc.R:958-960`); `batch_col` (ComBat, no explicit default beyond first detected column); `ruvm_group_col` (RUVm, first sheet column); `ruvm_k` (numericInput, default 1, min 1, max 10, `mod_methyl_qc.R:966`).

**Reactive dependencies.** `batch_qc_result` — `eventReactive(input$run_batch_btn, ...)` (`mod_methyl_qc.R:974-994`) — validates `input$batch_method %in% c("combat","ruvm")` and `length(sg$included) >= 1` (`mod_methyl_qc.R:975,977`).

**Functions called.** `methyl_pca_scores(mat)` (`qc.R:576-587`, both before and after correction, `mod_methyl_qc.R:980,992`); `methyl_sheet_sample_ids()` (`qc.R:463-468`); `methyl_batch_correct_combat(mat, batch)` (`qc.R:723-750`) or `methyl_batch_correct_ruvm(mat, rg, group, k)` (`qc.R:764-808`), dispatched on `input$batch_method` (`mod_methyl_qc.R:986,990`).

**Statistical operations — ComBat path (`qc.R:723-750`):** requires ≥2 distinct batch levels each with ≥2 samples (`qc.R:729-736`, matching `sva::ComBat()`'s own requirement). Beta values are logit-transformed to M-values first via `methyl_beta_to_mvalue()` (`qc.R:566-569`, clipped at `[1e-4, 1-1e-4]`), because — per the function's comment citing Du et al. 2010 — "ComBat assumes a roughly Gaussian, unbounded outcome; beta values are bounded... M-values are the better-behaved, more homoscedastic scale" (`qc.R:713-719`). `sva::ComBat(dat=m[keep,], batch=batch, par.prior=TRUE, mean.only=FALSE)` is applied only to complete-case rows (`qc.R:738-742`); incomplete rows are passed through unchanged (`qc.R:748`). Result converted back to beta (`qc.R:747`).

**Statistical operations — RUVm path (`qc.R:764-808`):** requires ≥2 levels in the factor-of-interest column (`qc.R:772-773`) and `missMethyl::getINCs(rg_set)` (internal negative-control probes, `qc.R:775`). Fits `missMethyl::RUVfit(Y=t(mc), X=design[,2], ctl=ctl, k=k, method="ruv4")` (`qc.R:795`) then `missMethyl::RUVadj()` (`qc.R:799`). **Code evidence on the explicit `method="ruv4"` choice:** the function's comment states `RUVfit()`'s own default method `"inv"` "ignores `k` entirely... which would make the 'Unwanted-variation factors (k)' control above silently do nothing" — `ruv4` is chosen specifically because it is the k-driven variant (`qc.R:789-794`), a deliberate, disclosed API-behavior fix rather than an oversight.

**Plots/tables.** `batch_pca_before`/`batch_pca_after` (side-by-side plotly PCA scatter colored by batch/factor, via shared helper `.batch_pca_plot()`, `mod_methyl_qc.R:1020-1039`); `batch_variance_table` (DT, per-PC variance explained before vs. after, `mod_methyl_qc.R:1041-1049`).

**Output objects.** `batch_qc_result()` — list `out` (correction result: `ok`, `corrected`, etc.), `before`/`after` (PCA score objects), `batch`, `method`, `subgroup`, `run_at` (`mod_methyl_qc.R:993`).

**Interpretation.** A successful correction should visually tighten the batch/factor separation on PC1–PC2 in the "after" panel relative to "before"; the variance-explained table gives a numeric companion.

**Connections to other tabs.** None — this tab reads `current_subgroup()$mat` directly, independent of Probe QC. Code comment confirms: "Runs on the current stratum's raw matrix, independently of Probe QC — this tab never depends on whether/how Probe QC was run" (`mod_methyl_qc.R:938-939`). The corrected matrix (`out$corrected`) is held only inside `batch_qc_result()` and is **not** offered as a download or fed into any other tab — see §14, Finding M-1.

---

### 3.6 Tab 6: Outlier QC

**Purpose.** Four independent sample-outlier detection methods (PCA-distance, hierarchical-clustering singleton, correlation-based, Mahalanobis-distance), collapsed into one ranked "outlier score" table, with explicit, separate flag-vs-exclude steps.

**Input data.** `current_subgroup()$mat`.

**Required/optional columns.** None.

**Parameters and defaults** (`mod_methyl_qc.R:1063-1074`):
| Input ID | Type | Default | Notes |
|---|---|---|---|
| `outlier_methods` | checkboxGroupInput | `c("pca","hclust")` | of pca/hclust/correlation/mahalanobis |
| `pca_sd` | numericInput | 3 | SD distance threshold, shown only if `pca` selected |
| `corr_k` | numericInput | 3 | MAD multiplier, shown only if `correlation` selected |
| `mahal_alpha` | numericInput | 0.01 | chi-squared alpha, shown only if `mahalanobis` selected |

Note `hclust` has no tunable threshold exposed in the UI even though `methyl_sample_outliers_hclust()` accepts a `height_frac` parameter defaulted to 0.5 in the function signature (`qc.R:291`) — see §14, Finding L-1.

**Reactive dependencies.** `outlier_qc_result` — `eventReactive(input$run_outlier_btn, ...)` (`mod_methyl_qc.R:1083-1112`) — requires `length(sg$included) >= 4` (`mod_methyl_qc.R:1085-1086`) and at least one method selected (`mod_methyl_qc.R:1089`).

**Functions called.** `methyl_sample_outliers_pca(mat, sd_threshold)` (`qc.R:274-286`, called unconditionally at `mod_methyl_qc.R:1093` — see note below); `methyl_sample_outliers_hclust(mat)` (`qc.R:291-304`, called unconditionally at `mod_methyl_qc.R:1097`); `methyl_sample_outliers_correlation(mat, k)` (`qc.R:255-263`, called unconditionally at `mod_methyl_qc.R:1100`); `methyl_sample_outliers_mahalanobis(mat, alpha)` (`qc.R:319-338`, called unconditionally at `mod_methyl_qc.R:1103`); `methyl_outlier_score_table(sample_qc)` (`qc.R:529-539`, called `mod_methyl_qc.R:1106`); `methyl_outlier_diagnostic_table()` (`qc.R:681-696`, used by the diagnostic plot).

**Code evidence — all four detectors always run, only their columns are exposed:** every detector function call at `mod_methyl_qc.R:1093,1097,1100,1103` executes unconditionally regardless of `methods_sel`; the `if ("pca" %in% methods_sel) ...` / etc. guards (`mod_methyl_qc.R:1094,1098,1101,1104`) only control whether that method's flag column is *added to the displayed table* — see §14, Finding L-2 for the performance implication.

**Statistical operations, per method:**
- **PCA-distance** (`qc.R:274-286`): top 5,000 variance-ranked, complete-case probes → `prcomp(..., scale.=TRUE)` → Euclidean distance from centroid on PC1–PC2 in units of each PC's own SD, flagged when `> sd_threshold`. Function's own comment: "a common, easy-to-explain ad hoc heuristic... but it is NOT a formal statistical test — it only uses 2 of the retained PCs and has no citable significance threshold" (`qc.R:265-273`).
- **Hierarchical clustering** (`qc.R:291-304`): average-linkage `hclust()` on Euclidean distance over the same top-5,000-variance probe subset; `cutree(hc, h = max(height)*0.5)`; any sample in a singleton cluster is flagged.
- **Correlation-based** (`qc.R:255-263`): mean pairwise correlation with every other sample (top-5,000-variance probes); flagged when `< median(mean_cor) - k*mad(mean_cor)`. Comment notes this uses median/MAD (robust) rather than WGCNA's mean/SD convention, "a deliberate choice, not an attempt to reproduce WGCNA's exact formula" (`qc.R:241-254`).
- **Mahalanobis distance** (`qc.R:319-338`): top-10 PCA components (not raw probes, since p≫n there) → `stats::mahalanobis(scores, colMeans(scores), stats::cov(scores))` → flagged against `qchisq(1-alpha, df=k)`. **Code-disclosed limitation, verbatim from the function's own comment:** "this uses the classical (non-robust) mean/covariance (`stats::cov()`), which is itself sensitive to the very outliers it's trying to detect (masking/swamping) — a robust estimator (e.g. MCD...) would be more resistant, at the cost of needing more samples than this module can assume. Treat this as a reasonable, disclosed simplification, not a claim of being the most robust method available" (`qc.R:306-318`).

**Plots/tables.** `outlier_score_table` (DT, ranked by `outlier_score` = count of selected methods flagging each sample, `mod_methyl_qc.R:1149-1153`); `pca_outlier_plot` (lazy-loaded plotly PC1/PC2 scatter, flagged vs. within-range, `mod_methyl_qc.R:1162-1173`); `outlier_diagnostic_plot` (lazy-loaded, PCA-distance-from-centroid bar chart via `methyl_plot_outlier_diagnostic()`, `qc.R:912-919`, `mod_methyl_qc.R:1176-1184`); `dendro_plot` (lazy-loaded base-graphics dendrogram via `ggdendro::dendro_data()`, `mod_methyl_qc.R:1186-1201`, only shown if `hclust` was among the selected methods).

**Output objects.** `outlier_qc_result()` — list `mat`, `subgroup`, `sample_qc` (per-sample flags), `outlier_scores`, `pca_detail`, `mahal_detail`, `corr_detail`, `hc_res`, `methods`, `run_at` (`mod_methyl_qc.R:1109-1111`).

**Flag-vs-remove distinction — explicit code confirmation:** the tab's own header comment states "Detection and removal are two separate, explicit steps — this tab only ever FLAGS samples; nothing is excluded until the user selects rows and clicks 'Apply Sample Exclusions' below" (`mod_methyl_qc.R:1052-1054`). Verified: `observeEvent(input$apply_outlier_exclusions_btn, ...)` (`mod_methyl_qc.R:1155-1160`) is the only code path that writes to `manual_exclude`; no automatic exclusion occurs from computing `outlier_qc_result()` alone.

**Connections to other tabs.** Reads `current_subgroup()$mat` directly (independent of Probe QC / Sample QC results); writes to the same shared `manual_exclude` reactiveVal as Sample QC and Sex QC.

---

### 3.7 Tab 7: Visualizations

**Purpose.** A pure re-plotting tab: every card here either re-visualizes data another tab (almost exclusively Probe QC) has already produced, or visualizes raw-IDAT diagnostics directly from `methyl_dataset`. Its own header comment: "nothing here computes a NEW analysis of its own; it only re-plots what Probe QC / raw IDAT diagnostics already have available" (`mod_methyl_qc.R:1204-1207`).

**Input data.** `probe_qc_result()$filtered`/`$mat` (for the eight Probe-QC-gated plots); `methyl_dataset$detp`/`$beadcount`; `current_rg_subset()` (control-probe heatmap); `current_subgroup()$mat` (detection-p / bead-count plots, restricted to the current stratum).

**Required/optional columns.** None beyond what Probe QC/IDAT already require. `viz_color_by` (`selectInput`, populated from all `sample_sheet` columns plus a "no coloring" option, `mod_methyl_qc.R:1256-1261`) is fully optional and purely cosmetic.

**Parameters.** `viz_color_by` only (default `"__none__"`). Individual plots have no other tunable parameters; each is a "Generate" button (`lazy_plot_ui()`, see §10) rather than an auto-recompute.

**Reactive dependencies — three independent gating layers, verified in code:**
1. **Dataset-level shell** (`output$viz_ui`, `mod_methyl_qc.R:1216-1240`): depends only on `methyl_dataset$beta` and dataset-level facts (`has_ctrl`/`has_detp`/`has_beadcount`), explicitly *not* on `probe_qc_has_run()` or any `plot_shown` flag, per its own comment (`mod_methyl_qc.R:1209-1215`) — this shell renders once per dataset load and is never rebuilt by a Run/Generate click.
2. **Probe-QC-gated sub-shell** (`output$viz_probe_qc_gate`, `mod_methyl_qc.R:1247-1283`): re-renders *only* when `probe_qc_has_run()` flips, offering a "Go to Probe QC" button (`mod_methyl_qc.R:1253`, wired via `updateTabsetPanel`, `mod_methyl_qc.R:1285-1287`) if it hasn't run yet, else the eight PCA/MDS/density/etc. Generate buttons.
3. **Per-plot lazy rendering** (`lazy_plot_ui()`/`plot_shown` reactiveValues, `mod_methyl_qc.R:131-168`): each plot's own `render*` output is gated via `req(isTRUE(plot_shown[[pid]]))` — clicking one plot's Generate button never re-renders any other plot's already-drawn output (this was the specific bug the module's header comment describes fixing, `mod_methyl_qc.R:93-122`).

**Functions called.** `methyl_pca_scores()` (`qc.R:576-587`, 2D & 3D); `methyl_mds_scores()` (`qc.R:624-637`); `methyl_beta_density_sample()` (`qc.R:607-616`, density/boxplot/violin share one `viz_density_df` reactive, `mod_methyl_qc.R:1338-1349`); `methyl_sample_correlation()` (`qc.R:593-601`); `methyl_mean_sd_table()` (`qc.R:643-650`); `methyl_control_probe_matrix()` (`qc.R:656-675`); plotting builders `methyl_plot_scatter2d()`, `methyl_plot_density()`, `methyl_plot_boxplot()`, `methyl_plot_violin()`, `methyl_plot_corr_heatmap()`, `methyl_plot_mean_sd()`, `methyl_plot_detp_heatmap()`, `methyl_plot_beadcount_dist()` (all `qc.R:818-919`, cited individually in §5).

**Plots (11 total, plus a WebGL 3D scatter not built through `ggplot2`).** PCA 2D (`viz_pca_plot`), PCA 3D via `plotly::plot_ly(..., type="scatter3d")` — "SVG export isn't supported for 3D plots (PNG only, via the camera icon)" (`mod_methyl_qc.R:1265`), MDS (`mds_plot`), beta-density before/after (`density_plot`), per-sample boxplot (`boxplot_plot`), per-sample violin (`violin_plot`), sample correlation heatmap (`corr_heatmap`, hierarchically ordered via `hclust()` inside `methyl_plot_corr_heatmap()`, `qc.R:899`), mean-SD plot (`mean_sd_plot`, with a loess trend line), detection-p heatmap (`detp_heatmap`, IDAT-only, 150 randomly-sampled probes, `mod_methyl_qc.R:1392`), bead-count histogram (`beadcount_dist`, IDAT-only), control-probe heatmap (`control_heatmap`, IDAT-only, top 60 highest-variance control probes, `mod_methyl_qc.R:1415`).

**Output objects.** No persisted result object — every plot is computed fresh inside its own `render*` block from `probe_qc_result()`/`current_subgroup()`/`methyl_dataset` each time its Generate button is clicked (not cached across clicks).

**Interpretation.** Standard sample-structure/QC-diagnostic reading: PCA/MDS clustering by metadata (`viz_color_by`) suggests batch or biological structure; density/boxplot/violin before-vs-after show the effect of Probe QC's filtering on the beta-value distribution; the correlation heatmap and mean-SD plot are classic array-QC diagnostics (Huber et al. 2002 convention for mean-SD, cited in `qc.R:639-641`).

**Connections to other tabs.** The Probe-QC dependency (§3.3 above) is the module's one real cross-tab data flow. All other plots are dataset-level (IDAT diagnostics) or stratum-level (`current_subgroup()`), independent of every other tab's *results*.

---

### 3.8 Tab 8: Reports & Export

**Purpose.** Aggregates whatever subset of the other seven tabs has actually been run this session into downloadable artifacts: a CSV summary, filtered matrices, a self-contained HTML report, an optional PDF report, a figures ZIP, a reproducibility settings table, and an equivalent-Bioconductor-code snippet.

**Input data.** `current_qc_pieces()` — a `reactive` resolving each of `overview`/`sample_qc`/`probe_qc`/`sex_qc`/`outlier_qc`/`batch_qc` to that tab's result object if `*_has_run()` is `TRUE`, else `NULL` (`mod_methyl_qc.R:1441-1448`) — "so 'which tabs have run' is resolved exactly once per call" (comment, `mod_methyl_qc.R:1439-1440`).

**Required/optional columns.** None directly; entirely derived from other tabs' already-validated inputs.

**Parameters.** None — every export is a `downloadButton`/`downloadHandler` pair; no numeric/threshold inputs live on this tab itself.

**Reactive dependencies.** `pdf_report_available()` (`reactive`, checks `rmarkdown::pandoc_available()` and a LaTeX toolchain via `Sys.which()`/`tinytex::is_tinytex()`, `mod_methyl_qc.R:1431-1437`); `report_plots()` (`reactive`, builds figures once per download rather than per handler, `mod_methyl_qc.R:1538-1541`).

**Functions called (with citations).** `methyl_qc_summary_table()` (`qc.R:985-1034`, one row per metric, explicitly distinguishing "not run" from "ran, no result" for sex/batch QC, `qc.R:1009-1013`); `methyl_qc_report_plots()` (`qc.R:1045-1089`, rebuilds every figure as a plain `ggplot` object — "never a plotly widget" since `htmltools::save_html()`/`ggplot2::ggsave()` both need a static image, `qc.R:1036-1038` — gated per-figure-group on `probe_qc`/`outlier_qc` being non-`NULL`, with each unavailable group added to a `skipped` character vector with a reason); `methyl_qc_report_html()` (`qc.R:1099-1134`, base64-embeds every PNG via `base64enc::dataURI()` into one self-contained HTML file via `htmltools::save_html()` — explicitly built to avoid an `rmarkdown`/`pandoc` dependency, `qc.R:1091-1095`); `methyl_qc_report_zip()` (`qc.R:1141-1157`, shells out to a `zip` binary via `utils::zip()`, degrading with a reason if unavailable); `methyl_qc_r_code()` (`qc.R:951-973`, a static text template reflecting Probe QC's actual settings, explicitly "not a live-executed pipeline", `qc.R:947-950`); `methyl_beta_to_mvalue()` (`qc.R:566-569`, for the M-value download).

**Downloads offered, each independently gated:**
| Download | Handler | Gate | Citation |
|---|---|---|---|
| Filtered beta matrix (CSV) | `dl_filtered_beta` | `probe_qc_has_run()` | `mod_methyl_qc.R:1492-1499` |
| Filtered M-value matrix (CSV) | `dl_filtered_mvalue` | `probe_qc_has_run()` | `mod_methyl_qc.R:1500-1507` |
| Probe filter summary (CSV) | `dl_filter_summary` | `probe_qc_has_run()` | `mod_methyl_qc.R:1508-1518` |
| Sample QC table (CSV) | `dl_sample_qc` | `sample_qc_has_run()` | `mod_methyl_qc.R:1519-1525` |
| QC summary (CSV) | `dl_qc_summary` | always available (rows report "not run" per-piece) | `mod_methyl_qc.R:1526-1534` |
| QC report (HTML) | `dl_report_html` | always available | `mod_methyl_qc.R:1543-1555` |
| QC report (PDF) | `dl_report_pdf` | `pdf_report_available()` | `mod_methyl_qc.R:1557-1585` |
| All figures (ZIP) | `dl_figures_zip` | always available (degrades per `methyl_qc_report_zip()`) | `mod_methyl_qc.R:1587-1595` |
| Probe QC reproducibility table | n/a (DT) | `probe_qc_has_run()` | `mod_methyl_qc.R:1597-1604` |
| Copy R code (Probe QC) | `copy_code_btn` (clipboard JS) | `probe_qc_has_run()` | `mod_methyl_qc.R:1606-1616` |

**Statistical operations.** None new — this tab performs no computation of its own beyond re-deriving M-values (`methyl_beta_to_mvalue()`) and assembling already-computed pieces into export formats.

**Plots/tables.** The HTML/PDF reports embed every figure `methyl_qc_report_plots()` was able to build (probe-filtering cascade, beta density/boxplot/violin, correlation heatmap, PCA 2D, MDS, mean-SD, outlier PCA-distance diagnostic — up to 9 figures, fewer if Probe QC/Outlier QC haven't run) plus a metric/value summary table and a "Skipped / not yet run" list (`qc.R:1121-1128`).

**Output objects.** Files only (no in-app result object persists past the download).

**Interpretation.** The reproducibility table and R-code panel exist specifically so a user can hand-verify or manually reproduce this session's Probe QC filtering outside the app in `minfi`/`ChAMP`/`wateRmelon`.

**Connections to other tabs.** This tab is a pure sink — it reads every other tab's result object (via `current_qc_pieces()`) but writes nothing back to any of them; it is the one place in the module where all eight tabs' state is read simultaneously, though only for export, never to alter any tab's own behavior.

---

## 4. End-to-end pipeline diagram

**Code evidence — the architecture is NOT a strict sequential pipeline.** This was verified directly by reading every `eventReactive` in `mod_methyl_qc.R` for reads of another tab's `*_result()`: the only such read anywhere in the file is Visualizations reading `probe_qc_result()` (§3.3/§3.7, `mod_methyl_qc.R:1305,1330,1341,1373,1382`). No other tab reads another tab's result object. `current_subgroup()` (`mod_methyl_qc.R:342-346`) and `manual_exclude()` (`mod_methyl_qc.R:280`) are the only reactives read by more than two tabs, and both are explicitly "data scoping ('which samples are we looking at'), not a QC method's output" per the file's own header comment (`mod_methyl_qc.R:33-39`). This confirms the header comment's claim rather than merely repeating it.

Rendered as a diagram (arrows are genuine code-level reads; dotted arrows are the shared, always-live scoping state, not a QC result):

```
methyl_dataset (Dataset tab, mod_methyl_dataset.R)
   │  $beta  $sample_sheet  $rg_set  $mset  $detp  $beadcount  $array_type  $input_scale  $preloaded  $source
   ▼
current_subgroup()  ◄────────────────────────────────┐  (stratum filter + manual_exclude(), always live)
   │                                                  │
   ├──► Overview        (own eventReactive, own button)   │
   ├──► Sample QC        (own eventReactive, own button) ──┤ writes manual_exclude()
   ├──► Probe QC         (own eventReactive, own button) ──┼───────────────┐
   ├──► Sex QC           (own eventReactive, own button) ──┤ writes manual_exclude()
   ├──► Batch QC         (own eventReactive, own button)   │
   ├──► Outlier QC        (own eventReactive, own button) ──┘ writes manual_exclude()
   │                                                        │
   ▼                                                        ▼
Visualizations  ◄── reads probe_qc_result()$filtered / $mat (the ONE genuine cross-tab dependency)
   │
   ▼
Reports & Export  ◄── reads ALL SIX of {overview, sample_qc, probe_qc, sex_qc, outlier_qc, batch_qc}_result()
                       via current_qc_pieces(), read-only, export-only — no tab is altered by this read
```

Each of the six left-column tabs is triggered **only** by its own `actionButton` (verified per tab in §3), reads `current_subgroup()` fresh at the moment its own button is clicked, and does not react to any other tab's button, result, or run state — except that every `*_has_run` flag is reset to `FALSE` (not recomputed) when `methyl_dataset$beta` changes identity, i.e. a new dataset is loaded (`mod_methyl_qc.R:288-292`). Batch QC and Outlier QC/Sex QC/Sample QC/Overview do not write to or read from each other under any circumstance found in the code.

---

## 5. Function-by-function documentation

Every `methyl_*` function in `qc.R`, `idat_metrics.R`, `annotation.R`, and `parse_upload.R` that the QC module uses (or defines but never calls, flagged as such), plus the key base-R/package functions each wraps.

### 5.1 `qc.R` — probe filters

**`methyl_row_vars(m)`** — `qc.R:22-28`. Purpose: row-wise variance, `matrixStats::rowVars(m, na.rm=TRUE)` when available, else `apply(m, 1, var, na.rm=TRUE)`. Called by: every top-variance-probe selection in the file (`methyl_filter_variance`, `methyl_filter_sd`, `methyl_sample_outliers_pca/hclust/mahalanobis`, `methyl_pca_scores`, `methyl_sample_correlation`, `methyl_mds_scores`, `methyl_mean_sd_table`, the control-probe heatmap's top-60 selection). QC significance: the single performance-critical primitive in the file — its comment states `apply()` alone "takes long enough (several seconds, paid on every recompute) to make the live QC tool feel broken" at 400k+ probes (`qc.R:14-21`). Audit: correct, appropriately optimized, degrades gracefully without `matrixStats`.

**`methyl_filter_missing(mat, max_na_frac=0)`** — `qc.R:30-34`. Removes probes whose `rowMeans(is.na(mat))` exceeds `max_na_frac`. Default UI threshold 0 (zero tolerance, `mod_methyl_qc.R:682`). Audit: correct, straightforward.

**`methyl_filter_variance(mat, min_variance=0)`** — `qc.R:36-41`. Keeps probes with `methyl_row_vars(mat) >= min_variance`; non-finite variances coerced to 0 before comparison (`qc.R:38`), so an all-NA probe is treated as zero-variance and removed by any positive threshold. Audit: correct; the non-finite coercion is a sensible, disclosed edge-case handling (not silent — the probe would already be caught by `methyl_filter_missing` in most configurations).

**`methyl_filter_sd(mat, min_sd=0)`** — `qc.R:47-52`. `sqrt(pmax(methyl_row_vars(mat), 0))`, same edge-case handling as above. Explicitly documented as "the same computation... offered as its own independent optional filter since 'minimum SD' is the more familiar unit" (`qc.R:43-46`) — i.e. `methyl_filter_variance` and `methyl_filter_sd` are intentionally redundant/duplicate in substance (mathematically the same cut, different units), not an accidental duplication (see §14, Finding L-3 for the UI implication of offering both).

**`methyl_filter_mean_range(mat, lo, hi)`** — `qc.R:54-58`. Keeps probes with `rowMeans(mat, na.rm=TRUE)` inside `[lo, hi]`; `NA` means are removed (`!is.na(m) & ...`, `qc.R:56`). Audit: correct.

**`methyl_filter_non_cpg(mat)`** — `qc.R:63-66`. Wraps `methyl_probe_is_cpg()` (`annotation.R:99-101`, `grepl("^cg", ids, ignore.case=TRUE)`). Illumina-array-only by design (caller gates on `is_illumina_array()`, `mod_methyl_qc.R:685-692,744`). Audit: correct for its stated scope; purely ID-prefix-based, no manifest lookup, so it works even without annotation-package availability.

**`methyl_filter_snp(mat, anno_result)`** — `qc.R:68-81`. Degrades to all-`TRUE` with the annotation's own failure reason if `anno_result$ok` is `FALSE` (`qc.R:69`). Otherwise flags a probe if any of the manifest's `Probe_rs`/`CpG_rs`/`SBE_rs` columns is non-`NA` and non-empty (`qc.R:74-78`). Audit: correct implementation of "does this probe overlap a known dbSNP variant per the array's own manifest", but see §13 — this is a manifest-based SNP check, not a population-MAF-based common-SNP check (that is `methyl_filter_maf`, a separate, user-upload-only filter).

**`methyl_filter_sex_chr(mat, anno_result, mode="remove_xy")`** — `qc.R:89-104`. Three modes: `"keep"` (no-op, all-`TRUE`, `qc.R:90-92`), `"remove_xy"` (drops chrX+chrY), `"remove_y_only"` (drops chrY only, "useful for studies that want to retain X-linked signal", `qc.R:85-86`). Degrades gracefully without annotation. Audit: correct; the three-mode design (rather than a boolean) is a genuine, disclosed usability choice (`qc.R:83-88`).

**`methyl_filter_cross_reactive(mat, exclusion_ids=NULL)`** — `qc.R:112-119`. No bundled blacklist by design — comment explicitly states "fabricating one would violate this project's own evidence-based-methods requirement" (`qc.R:106-111`), citing Chen et al. 2013 (450K) and Pidsley et al. 2016 / McCartney et al. 2016 (EPIC) as the published lists a user could upload instead. Requires a user-uploaded exclusion list parsed by `methyl_parse_probe_list()`. No-op (all-`TRUE`) with an explanatory note otherwise. Audit: a defensible, disclosed scientific-integrity choice rather than an implementation gap.

**`methyl_parse_maf_list(datapath, filename)`** — `qc.R:131-144`. Parses a `probe_id,maf` CSV/TSV via `data.table::fread()`; MAF column found by name (`maf`/`MAF`/`af`/`AF`) or positionally (second column). Returns `ok=FALSE` with a specific error message if unparseable or no numeric MAF column found. Audit: correct, fails soft.

**`methyl_filter_maf(mat, maf_table=NULL, max_maf=0.05)`** — `qc.R:149-160`. Removes probes whose uploaded MAF exceeds `max_maf`; probes absent from the table are kept ("absence isn't evidence of a common SNP, so it's not penalized", `qc.R:146-148`). No bundled population-MAF table, same design rationale as `methyl_filter_cross_reactive`. Audit: correct, conservative default behavior for unmatched probes.

**`methyl_filter_detection_p(mat, detp, threshold=0.01)`** — `qc.R:162-172`. IDAT-only; degrades to all-`TRUE` otherwise (`qc.R:163-166`). A probe is removed if it fails `p > threshold` in **any** sample among `intersect(rownames(mat), rownames(detp))` (`qc.R:167-171`). Audit: standard, matches `minfi`/`ChAMP`-style detection-p filtering conventions (though see §14 for the "any sample" strictness discussed alongside beadcount below).

**`methyl_filter_beadcount(mat, beadcount, threshold=3)`** — `qc.R:178-188`. IDAT-only. A probe is removed if its bead count is `< threshold` in **any** sample (`qc.R:184-186`). **Code-disclosed, deliberate divergence from ChAMP:** the function's own comment states this is "stricter than ChAMP's own default (`champ.filter(beadCutoff=0.05)` removes a probe only when it fails in >5% of samples), a deliberate choice here rather than an attempt to reproduce ChAMP's exact rule" (`qc.R:174-177`). Flagged as an Informational audit finding (Finding I-1, §14) purely because it may surprise a user expecting ChAMP-equivalent numeric output, not because it is scientifically wrong — "any sample fails" is itself a legitimate, more conservative QC convention.

### 5.2 `qc.R` — sample-level QC

**`methyl_sample_call_rate(mat)`** — `qc.R:192-194`. `1 - colMeans(is.na(mat))`. Called by Overview and Sample QC. Audit: correct, trivial.

**`methyl_sample_failed_probe_pct(mat, detp, threshold=0.01)`** — `qc.R:202-217`. IDAT-only. Per-sample % of probes with `detection_p > threshold`, restricted to probes/samples common between `mat` and `detp`; samples not covered by `detp` are left `NA` rather than silently zero (`qc.R:214-215`). Called by Sample QC. Audit: correct, deliberately tolerant of ID mismatches, consistent with the per-probe filter's own `intersect()` convention.

**`methyl_sample_low_intensity(median_int_result, min_intensity=10)`** — `qc.R:225-232`. Score = mean of median methylated/unmethylated log2 intensity; flagged when `< min_intensity`. Built on `methyl_median_intensity()` (`idat_metrics.R`), mirroring `minfi`'s own `mMed`/`uMed` QC-tutorial convention per its comment (`qc.R:219-224`). Audit: correct, standard.

**`methyl_sample_outliers_iqr(x, k=1.5)`** — `qc.R:234-239`. Standard Tukey IQR-fence outlier flag on a numeric vector `x`. **Audit finding: defined but never called anywhere in the codebase** (confirmed by exhaustive grep across `R/methylomics/*.R`) — dead code (see §14, Finding L-4).

**`methyl_sample_outliers_correlation(mat, n_features=5000, k=3)`** — `qc.R:255-263`. Mean pairwise correlation per sample (top-variance probes via `methyl_sample_correlation()`), flagged when `< median - k*mad`. Called by Outlier QC. Audit: correct, robust-statistic choice explicitly justified in its own comment (`qc.R:241-254`).

**`methyl_sample_outliers_pca(mat, n_features=5000, sd_threshold=3)`** — `qc.R:274-286`. `prcomp(t(top-variance probes), scale.=TRUE)`; Euclidean centroid distance on PC1–PC2 in SD units. Called by Outlier QC. Requires ≥10 complete-case probes and ≥4 samples (`qc.R:276-277`), else degrades with a reason. Audit: correct but explicitly non-rigorous per its own comment (`qc.R:265-273`, quoted in full in §3.6).

**`methyl_sample_outliers_hclust(mat, n_features=5000, height_frac=0.5)`** — `qc.R:291-304`. Average-linkage `hclust()`, `cutree()` at `height_frac` of max height, singleton clusters flagged. Called by Outlier QC (with the default `height_frac`, unexposed in UI). Audit: correct implementation of a standard heuristic; `height_frac` is a reasonable but arbitrary constant not surfaced for tuning (Finding L-1, §14).

**`methyl_sample_outliers_mahalanobis(mat, n_features=5000, n_pcs=10, alpha=0.01)`** — `qc.R:319-338`. Classical Mahalanobis distance on top-10 PCA components, `qchisq(1-alpha, df=k)` threshold. Called by Outlier QC. Requires `k >= 2` components relative to sample count (`qc.R:327-330`); degrades if the covariance matrix is singular (`qc.R:332-335`). Audit: statistically correct for a *classical* (non-robust) Mahalanobis distance; the non-robustness is explicitly disclosed in the function's own comment (quoted in §3.6) as a known masking/swamping risk, citing Rousseeuw & Van Driessen 1999 and Filzmoser et al. 2005 as the more-robust alternative this deployment does not implement.

**`methyl_cluster_sex(y, reported=NULL)`** — `qc.R:356-376`. `stats::kmeans(y, centers=2, nstart=25)` if `length(y)>=6` and `length(unique(y))>=2`, else a global median split. Higher-mean cluster assumed male unless `reported` resolves direction by majority concordance. Called by `methyl_sex_check()`'s fallback path. Audit: correct, and the k-means-over-median-split design choice is backed by a specific, code-cited real-cohort failure mode (492F/197M, 147 misclassified females under a naive median split) — see §3.4 for the full quotation.

**`methyl_sex_check_attach_mismatch(detail, reported_sex, method)`** — `qc.R:383-392`. Appends `reported_sex`/`sex_mismatch` columns when a reported-sex vector is supplied; `n_mismatch` stays `NA` otherwise (rather than 0, correctly distinguishing "not applicable" from "zero mismatches"). Called by `methyl_sex_check()`. Audit: correct.

**`methyl_sex_check(mat, anno_result, rg_set=NULL, reported_sex=NULL)`** — `qc.R:407-446`. Dispatches to `minfi::getSex()` (raw IDAT path) or the chrX/chrY beta-clustering fallback, full detail in §3.4. Called by the Sex QC tab. Audit: correct two-tier design with an explicitly weaker fallback correctly labeled as such in its own output text (`qc.R:438-442`).

**`methyl_sheet_sample_ids(sheet, all_ids)`** — `qc.R:463-468`. Resolves sample-sheet rows to matrix column IDs: matches a `sample`/`Sample`/`sample_id`/`Sample_ID` column if present; else, **only if `nrow(sheet) == length(all_ids)`**, assumes row order matches; else falls back to `rownames(sheet)`. Called throughout (subgroup filtering, sex check, batch correction, color-by metadata). **Code-disclosed prior-bug fix:** the comment explains this deliberately avoids `rownames(sheet)` as the row-order fallback because "a freshly `data.table::fread()`-parsed data.frame gets sequential integer rownames... which never equal a real sample ID, so using them silently turned every 'row order' case into zero matches" (`qc.R:454-462`) — i.e. this function's current form is itself a documented bug fix over a prior, more naive implementation. Audit: correct as written; the final `rownames(sheet)` fallback (`qc.R:467`) is acknowledged in the same comment as "preserving the old, already-broken behavior only for a case that was never going to align anyway" — a knowingly inert fallback, not a hidden risk.

**`methyl_qc_subgroup_filter(mat, sheet, group_col, level, min_n=3)`** — `qc.R:479-500`. Resolves a subgroup column + stratum level into the matching column subset of `mat`. No-op (returns everything, labeled "All samples (n=...)") when no sheet/column/level is selected. Sets `low_n` when the resulting sample count is below `min_n` (default 3). Called once, by `current_subgroup()` (`mod_methyl_qc.R:344`) and `stratum_all_samples()` (`mod_methyl_qc.R:365`). Audit: correct; the docstring notes this function's introduction fixed a prior bug where "picking 'Female' vs 'Male' there previously had no effect on any result" (`qc.R:473-476`) — again, a disclosed historical fix, not a currently open issue.

**`methyl_apply_manual_exclude(subgroup, excluded_ids)`** — `qc.R:507-519`. Subtracts `excluded_ids` from an already-subgroup-filtered result; updates `mat`/`included`/`excluded`/`label`/`low_n` (recomputed against a hard-coded threshold of 3, `qc.R:517` — same constant as `methyl_qc_subgroup_filter`'s own `min_n` default, but not parameterized through from it). Called once, by `current_subgroup()` (`mod_methyl_qc.R:345`). Audit: correct.

### 5.3 `qc.R` — outlier scoring, retention cascade, transforms

**`methyl_outlier_score_table(sample_qc)`** — `qc.R:529-539`. One row per sample, `outlier_score` = count of TRUE flags across whichever of `pca_outlier`/`hclust_outlier`/`mahalanobis_outlier`/`correlation_outlier`/`iqr_outlier` columns are present, sorted descending. Called by Outlier QC (`mod_methyl_qc.R:1106`). Note the `iqr_outlier` column name is anticipated in `flag_cols` (`qc.R:530`) but no caller ever populates it, since `methyl_sample_outliers_iqr()` is never called (see Finding L-4) — this is dead provision for a column that can never appear.

**`methyl_probe_retention_cascade(n_probes_start, filters)`** — `qc.R:547-560`. Applies each named filter in `filters` **in the order given** as a running AND, recording probes retained after each cumulative step — distinct from the "probes removed by this filter alone" counts in the filter-summary table (`mod_methyl_qc.R:807-816`), which are computed independently per filter without the cumulative AND. Called by Probe QC (`mod_methyl_qc.R:759`) and reused inside `methyl_qc_report_plots()` for the exported cascade figure. Audit: correct, order in `filters` (an R named list, insertion-ordered) is exactly the literal `if` sequence at `mod_methyl_qc.R:741-752`.

**`methyl_beta_to_mvalue(beta, eps=1e-4)`** — `qc.R:566-569`. Logit transform, `beta` clipped to `[eps, 1-eps]` first to avoid `log2(0)`/`log2(Inf)`, "matching how `minfi::logit2()` and `lumi::beta2m()` handle this same edge case" (`qc.R:563-565`). Called by the M-value download handler, ComBat, and RUVm. Audit: correct, standard.

### 5.4 `qc.R` — structure/visualization helpers

**`methyl_pca_scores(mat, n_features=5000, n_pcs=10)`** — `qc.R:576-587`. Top-variance-probe `prcomp(..., scale.=TRUE)`; returns scores for up to 10 PCs and `var_explained`. Requires ≥10 complete-case probes, ≥4 samples. Called by Batch QC (before/after) and Visualizations (PCA 2D/3D). Audit: correct; factored out specifically to avoid re-deriving PCA independently in multiple call sites (`qc.R:571-575`).

**`methyl_sample_correlation(mat, n_features=5000)`** — `qc.R:593-601`. Top-variance-probe `stats::cor()`. Requires ≥10 complete-case probes, ≥2 samples. Called by the correlation-outlier detector and the Visualizations correlation heatmap. Audit: correct.

**`methyl_beta_density_sample(mat, n_probes=5000, seed_probes=NULL)`** — `qc.R:607-616`. Random subsample of probes (or a caller-supplied fixed set via `seed_probes`, used to keep before/after comparisons on the *same* probe IDs) melted to long format. Called by Visualizations' density/boxplot/violin and the exported report figures. Audit: correct; `seed_probes` reuse across before/after is a sound comparability choice.

**`methyl_mds_scores(mat, n_features=5000, k=2)`** — `qc.R:624-637`. `stats::cmdscale()` (classical/metric MDS) on Euclidean distance over top-variance probes. Called by Visualizations' MDS plot. Audit: correct; explicitly offered as "a distinct sample-structure view from PCA... since they can disagree on a noisy dataset" (`qc.R:618-623`).

**`methyl_mean_sd_table(mat, n_probes=20000)`** — `qc.R:643-650`. Per-probe mean vs. SD on a random subsample of up to 20,000 probes. Called by Visualizations' mean-SD plot and the exported report. Audit: correct, references Huber et al. 2002's mean-SD-plot convention.

**`methyl_control_probe_matrix(rg_set)`** — `qc.R:656-675`. `minfi::getProbeInfo(rg_set, type="Control")` + raw green/red channel averaging, log2-transformed. IDAT-only. Called by Visualizations' control-probe heatmap. Audit: correct.

**`methyl_outlier_diagnostic_table(sample_qc, pca_detail, mahal_detail=NULL)`** — `qc.R:681-696`. Combines call rate (if present), PCA distance, Mahalanobis distance² into one per-sample table. Called by Outlier QC's diagnostic plot and the exported report figure. Audit: correct; explicitly documented as tolerant of either caller's differently-shaped `sample_qc` input (`qc.R:683-687`).

**`methyl_guess_batch_column(sheet)`** — `qc.R:704-709`. Same batch/chip regex as `mod_methyl_qc.R:944-945` uses inline. **Audit finding: defined but never called** — Batch QC's UI re-implements the identical `grep()` pattern directly rather than calling this function (compare `qc.R:706-707` to `mod_methyl_qc.R:944-945`, character-for-character identical regex) — genuine duplication with a provided-but-unused helper (see §14, Finding L-5).

**`methyl_batch_correct_combat(mat, batch)`** — `qc.R:723-750`, **`methyl_batch_correct_ruvm(mat, rg_set, group, k=1)`** — `qc.R:764-808`. Both documented in full in §3.5. Called by Batch QC only. Audit: both correct implementations of their respective published methods, each with a specific, code-disclosed design rationale (M-value scale for ComBat; `method="ruv4"` for RUVm) rather than default-API usage.

### 5.5 `qc.R` — plot builders

All in the "Shared ggplot builders" section (`qc.R:810-919`), explicitly built so "a figure's actual drawing logic lives in exactly one place instead of being re-implemented for the static/interactive/exported versions separately" (`qc.R:811-816`) — verified true: every plot builder below is called from both a live `render*` output in `mod_methyl_qc.R` (wrapped in `plotly::ggplotly()` where interactive) and from `methyl_qc_report_plots()` (`qc.R:1045-1089`) for the static export, with no separate plotting code path for either.

| Function | Location | Purpose | Called by (live) | Called by (export) |
|---|---|---|---|---|
| `methyl_plot_cascade` | `qc.R:818-826` | bar chart, probes retained per cascade step | `mod_methyl_qc.R:820` | `qc.R:1050` |
| `methyl_plot_detp_heatmap` | `qc.R:828-834` | probe×sample detection-p tile heatmap | `mod_methyl_qc.R:1395` | not exported |
| `methyl_plot_beadcount_dist` | `qc.R:836-841` | bead-count histogram + threshold line | `mod_methyl_qc.R:1406` | not exported |
| `methyl_plot_density` | `qc.R:851-856` | beta-value density, before/after stage | `mod_methyl_qc.R:1354` | `qc.R:1061` |
| `methyl_plot_boxplot` | `qc.R:863-869` | per-sample beta boxplot, before/after stage | `mod_methyl_qc.R:1360` | `qc.R:1062` |
| `methyl_plot_violin` | `qc.R:871-877` | per-sample beta violin, before/after stage | `mod_methyl_qc.R:1366` | `qc.R:1063` |
| `methyl_plot_mean_sd` | `qc.R:879-884` | mean vs. SD scatter + loess line | `mod_methyl_qc.R:1382` | `qc.R:1080` |
| `methyl_plot_scatter2d` | `qc.R:890-896` | generic 2D scatter (PCA/MDS/sex/batch/outlier) | multiple sites | `qc.R:1071,1077` |
| `methyl_plot_corr_heatmap` | `qc.R:898-910` | hierarchically-ordered sample correlation heatmap | `mod_methyl_qc.R:1375` | `qc.R:1066` |
| `methyl_plot_outlier_diagnostic` | `qc.R:912-919` | ranked bar chart of a distance metric | `mod_methyl_qc.R:1182` | `qc.R:1085` |

`.methyl_stage_fill` (`qc.R:848-849`) is a shared, dot-prefixed (internal, non-exported-convention) named-color vector recognized by `methyl_plot_density`/`methyl_plot_boxplot`/`methyl_plot_violin` for both the "Before/After filtering" and "Before/After normalization" label pairs, "so ... all recognize whichever pair of labels the caller's long-format data frame actually uses" (`qc.R:843-847`) — a cross-module-aware design choice, since "Before/After normalization" is used by the separate normalization sub-module, not QC itself.

### 5.6 `qc.R` — status, code-gen, summary, and report/export helpers

**`methyl_qc_status_badge(overview)`** — `qc.R:929-944`. The module's only pass/warning/fail decision function, full logic in §3.1. Called only by Overview (`mod_methyl_qc.R:508`). Audit: correct as designed; both thresholds (0.90, 10%) are hard-coded (Finding C-2, §14).

**`methyl_qc_r_code(settings)`** — `qc.R:951-973`. Builds a static, template-filled `minfi`/`ChAMP`/`wateRmelon` code snippet from Probe QC's actual settings object; explicitly "a static template filled in from the actual settings, not a live-executed pipeline" (`qc.R:947-950`). Called by the Reports & Export tab (`mod_methyl_qc.R:1608`). Audit: correct as a documentation aid; not itself validated to run (it is not executed by the app), so its Bioconductor-equivalence is asserted, not tested — reasonable given its stated purpose.

**`methyl_qc_summary_table(overview=NULL, sample_qc=NULL, probe_qc=NULL, sex_qc=NULL, outlier_qc=NULL, batch_qc=NULL)`** — `qc.R:985-1034`. One row per metric across all six independently-run tabs; explicitly distinguishes "not run" (`NULL`) from "ran, no result" (non-`NULL` with an internal `ok=FALSE`) for sex/batch QC specifically (`qc.R:1009-1013`), so the exported CSV never misrepresents a method that ran but produced nothing as one that was never attempted. Called by `dl_qc_summary` and both report handlers. Audit: correct, careful design around a real reporting-accuracy risk.

**`methyl_qc_report_plots(methyl_dataset, probe_qc=NULL, outlier_qc=NULL)`** — `qc.R:1045-1089`. Rebuilds every static export figure from `probe_qc`/`outlier_qc` (both possibly `NULL`), populating a `skipped` vector with a specific reason for each figure group that cannot be built — either because its source tab hasn't run, or (independently) because the underlying computation itself returned `ok=FALSE` (e.g. `cr <- methyl_sample_correlation(...)`, then `if (isTRUE(cr$ok)) ... else skipped <- c(skipped, sprintf("Sample correlation heatmap: %s", cr$reason))`, `qc.R:1065-1066`). Called by `report_plots()` (`mod_methyl_qc.R:1538-1541`), shared across all three export handlers so figures are built once per download. Audit: correct, and its double-source-of-skip logic (tab-not-run vs. computation-failed) is more granular than a single boolean would allow.

**`methyl_qc_report_html(methyl_dataset, summary_df, plots, skipped, subtitle=NULL)`** — `qc.R:1099-1134`. Base64-embeds every figure into one self-contained HTML file via `htmltools`/`base64enc`; degrades with a reason if either package is missing. Called by `dl_report_html`. Audit: correct; the "self-contained, no external file references" design is explicitly chosen over `rmarkdown` to avoid a pandoc dependency for the HTML path specifically (the PDF path, by contrast, does require `rmarkdown`+LaTeX, `mod_methyl_qc.R:1557-1585`).

**`methyl_qc_report_zip(plots)`** — `qc.R:1141-1157`. Shells out to a `zip` binary on `PATH` via `utils::zip()`; degrades with a reason ("not guaranteed on Windows", `qc.R:1136-1140`) if unavailable or if `status != 0`. Called by `dl_figures_zip`. Audit: correct, appropriately defensive around an external-tool dependency.

### 5.7 `idat_metrics.R`

**`methyl_idat_derive(rg_set)`** — `idat_metrics.R:13-28`. `minfi::preprocessRaw()` → `minfi::getBeta()`, plus `minfi::detectionP()` and `minfi::beadcount()`. Called once, from the Dataset tab's IDAT-upload handler (`mod_methyl_dataset.R:253`), not from within `mod_methyl_qc.R` itself — QC only ever *reads* the resulting `methyl_dataset$mset`/`$detp`/`$beadcount`, it never derives them. Audit: correct; explicitly "no normalization — normalization methods live in a later Methylomics sub-module" (`idat_metrics.R:9-10`), the clearest single-line statement of the QC/normalization boundary in the codebase.

**`methyl_bisulfite_conversion(rg_set)`** — `idat_metrics.R:33-40`. Wraps `wateRmelon::bscon()`. Called by Sample QC (`mod_methyl_qc.R:586`). Audit: correct, standard citation for this metric on Illumina arrays.

**`methyl_median_intensity(mset)`** — `idat_metrics.R:45-53`. Wraps `minfi::getQC()`, reshapes to a `sample`/`med_meth_log2`/`med_unmeth_log2` data.frame. Called by Sample QC (`mod_methyl_qc.R:587`) and, via `methyl_sample_low_intensity()`, feeds the low-intensity flag. Audit: correct.

### 5.8 `annotation.R`

**`methyl_get_annotation(array_type)`** — `annotation.R:48-93`. Reads `Locations`/`Manifest`/`Other`/`SNPs.147CommonSingle` data objects directly out of the relevant Bioconductor annotation package's namespace (450K or EPIC only — `METHYL_ANNOTATION_PACKAGES`, `annotation.R:18-21`), rather than via `minfi::getAnnotation()`. **Code-disclosed rationale for this unusual approach:** `getAnnotation()`'s own `IlluminaMethylationAnnotation` wrapper "triggers an internal `Biobase::updateObject()` step that fails unless the package is `library()`-attached... and attaching it pulls in `Biostrings`, which masks `base::strsplit()` for the entire R session — a real risk to every OTHER module in this app that calls `strsplit()` unqualified" (`annotation.R:38-47`) — a specific, app-wide side-effect avoidance, not a generic implementation choice. Cached per array type in `.methyl_anno_cache` for the process lifetime (`annotation.R:23-26,56-57,91`). Called by `anno_result()` (`mod_methyl_qc.R:260-263`), which every annotation-dependent filter and the beta-heuristic sex check read. Audit: correct and well-justified; EPICv2/WGBS/RRBS/Custom array have no annotation package configured and degrade with an explicit reason (`annotation.R:50-55`).

**`methyl_probe_is_cpg(probe_ids)`** — `annotation.R:99-101`. `grepl("^cg", probe_ids, ignore.case=TRUE)`. Called by `methyl_filter_non_cpg()`. Audit: correct, ID-prefix-only (no manifest dependency), consistent with standard Illumina probe-naming convention (`cg`=CpG, `ch`=non-CpG/CpH, `rs`/`nv`=SNP control probes, per `annotation.R:95-98`).

### 5.9 `parse_upload.R`

**`methyl_parse_matrix(datapath, filename)`** — `parse_upload.R:11-37`. `data.table::fread()`; first column = probe ID (checked for duplicates, `parse_upload.R:21-26`), remaining columns coerced to numeric via `storage.mode(m) <- "double"` (`parse_upload.R:28-29`). Fails soft with `ok=FALSE` on parse failure, duplicate IDs, or non-numeric data columns. Called by the Dataset tab's matrix-upload path (`mod_methyl_dataset.R:174`), not by QC directly. Audit: correct; documented in §6/§7 below for how it shapes QC's input assumptions.

**`methyl_parse_sample_sheet(datapath, filename)`** — `parse_upload.R:40-46`. Minimal `data.table::fread()` wrapper, no column validation beyond "at least one row". Called by the Dataset tab (`mod_methyl_dataset.R:178,233`). Audit: correct but minimal — see §7 (no required-column check at all).

**`methyl_parse_probe_list(datapath, filename)`** — `parse_upload.R:50-57`. Plain-text or first-CSV/TSV-column probe-ID list, deduplicated. Called by Probe QC's cross-reactive-probe upload handler (`mod_methyl_qc.R:729`). Audit: correct.

**`methyl_read_idat(files)`** — `parse_upload.R:65-96`. Validates both `_Grn.idat` and `_Red.idat` files are present (`parse_upload.R:77-81`), stages uploads under sanitized `basename()`s in a throwaway temp directory (explicit path-traversal protection against Shiny's client-controlled `name` field, `parse_upload.R:74-76`) before calling `minfi::read.metharray.exp()`. Called by the Dataset tab's IDAT-upload path (`mod_methyl_dataset.R:229,246`), not QC directly. Audit: correct, defensively written.

### 5.10 Key base-R / package functions relied upon

| Function | Package | Used for | QC significance |
|---|---|---|---|
| `stats::prcomp` | stats | PCA (scores, PCA-outlier, Mahalanobis basis) | Dimensionality reduction underlying most sample-structure diagnostics |
| `stats::hclust`/`cutree` | stats | hierarchical-clustering outlier detection, dendrogram, correlation-heatmap ordering | Average-linkage clustering throughout |
| `stats::mahalanobis`/`cov` | stats | Mahalanobis distance-based outlier flag | Classical (non-robust) multivariate distance |
| `stats::kmeans` | stats | sex clustering (`methyl_cluster_sex`) | 2-cluster split, `nstart=25` |
| `stats::cor` | stats | sample correlation matrix/heatmap, correlation-outlier detection | |
| `stats::cmdscale` | stats | MDS | Classical/metric MDS |
| `stats::qchisq` | stats | Mahalanobis threshold | chi-squared quantile, `df=k` PCs |
| `stats::mad`/`median` | stats | correlation-outlier threshold | Robust dispersion/center |
| `sva::ComBat` | sva | batch correction | Empirical-Bayes batch adjustment on M-values |
| `missMethyl::RUVfit`/`RUVadj`/`getINCs` | missMethyl | RUVm batch correction | Control-probe-based unwanted-variation removal |
| `minfi::getSex` | minfi | raw-IDAT sex check | Copy-number-based, chrX/chrY intensity |
| `minfi::getQC` | minfi | median intensity | mMed/uMed |
| `minfi::detectionP`/`beadcount` | minfi | detection p-value / bead count matrices | Derived at Dataset-tab load time |
| `wateRmelon::bscon` | wateRmelon | bisulfite conversion efficiency | Control-probe-based |
| `matrixStats::rowVars` | matrixStats | fast row variance | Performance-critical, ~2 orders of magnitude faster than `apply()` |
| `plotly::ggplotly`/`plot_ly` | plotly | interactive wrapping of every `ggplot2` figure; native 3D scatter for PCA(3D) | |

---

## 6. Input data documentation

**Formats accepted.** Three, chosen on the Dataset tab, not inside QC itself: (a) preloaded whole-blood dataset (`mod_methyl_dataset.R:18`), (b) CSV/TSV beta or M-value matrix upload (`mod_methyl_dataset.R:19,46-53`), (c) raw `.idat`/`.idat.gz` file upload, requiring both `_Grn` and `_Red` channels per sample (`mod_methyl_dataset.R:20,56-64`; `parse_upload.R:77-81`).

**Beta vs. M-value handling.** `methyl_dataset$input_scale` (`"beta"` or `"m"`) is set once, at load time, and never changed by QC. **Code evidence: QC does not itself branch its filtering logic on scale for most filters** — `methyl_filter_variance`/`_sd`/`_missing` operate identically regardless of scale (they are scale-agnostic statistics). The one place scale is respected is the mean-range filter's *default* bounds, which switch between `0.01`–`0.99` (beta) and `-6`–`6` (M-value) based on `identical(methyl_dataset$input_scale, "beta")` (`mod_methyl_qc.R:709-710`) — a UI-default convenience, not a computational branch (a user can still type any bounds regardless of scale). Batch correction (ComBat/RUVm) *does* branch on scale internally: both always logit-transform to M-values first via `methyl_beta_to_mvalue()` before correcting, then convert back to beta for display (`qc.R:737,747,783,805`) — meaning if the input was already M-values, `methyl_beta_to_mvalue()` (which assumes a `[0,1]` input) would silently misbehave. **This is a genuine, code-confirmed gap: `methyl_batch_correct_combat`/`methyl_batch_correct_ruvm` are never guarded by an `input_scale` check** — see §14, Finding H-1.

**Orientation.** Rows = probes/CpGs, columns = samples — confirmed in §1 and consistently throughout every function in `qc.R` (`rowMeans`/`rowVars` for probe stats, `colMeans` for sample stats). The upload parser enforces this shape by construction: `methyl_parse_matrix()` always treats the first column as the probe ID and every other column as a sample (`parse_upload.R:20,28`) — there is no option to upload a transposed (sample-by-probe) matrix, and **no runtime check exists to detect and reject a transposed upload** (see §7).

**Required vs. optional.** Required: `methyl_dataset$beta` (every tab's `renderUI`/`eventReactive` gates on `req(methyl_dataset$beta)`, e.g. `mod_methyl_qc.R:446,528,665,836,942,1057,1217,1451`). Optional: `sample_sheet` (its absence degrades subgroup filtering, sex comparison, and batch/RUVm to "unavailable" states with explanatory text, never a crash — confirmed at, e.g., `mod_methyl_qc.R:464,556-557,947-953`); `rg_set`/`mset`/`detp`/`beadcount` (all IDAT-derived, `NULL` for a matrix upload or the preloaded dataset).

**Missing values.** **Code evidence: `NA`s are permitted throughout and are handled explicitly, not rejected at upload.** `methyl_parse_matrix()` recognizes `"NA"`, `""`, `"NaN"`, `"null"`, `"NULL"`, `"#N/A"` as missing on read (`parse_upload.R:14`). Downstream, every statistical function that could be `NA`-sensitive explicitly passes `na.rm=TRUE` (row means/vars) or uses `stats::na.omit(mat)` (PCA/clustering/MDS/outlier functions, e.g. `qc.R:275,292,320,577,594,625`) to drop incomplete cases before matrix algebra that cannot tolerate `NA` (`prcomp`, `dist`, `cov`). There is no upstream requirement that a matrix be complete.

**Beta-value range validation.** **Code evidence: NOT validated anywhere in the pipeline.** An uploaded matrix declared as `"beta"` scale is never checked to actually fall within `[0,1]` — `methyl_parse_matrix()` only checks that values are numeric (`parse_upload.R:27-35`), and no QC filter or the Dataset tab itself calls `range()`/`min()`/`max()` on the uploaded matrix to sanity-check it against its declared scale. This was confirmed by grepping `mod_methyl_dataset.R` and `parse_upload.R` for any bounds check — none exists. See §7 and §14, Finding H-2.

---

## 7. Input validation audit

| Check | Status | Evidence |
|---|---|---|
| File type (CSV/TSV vs. other) | Partially implemented | `methyl_parse_matrix()`/`methyl_parse_sample_sheet()` both attempt `data.table::fread()` regardless of extension (`parse_upload.R:13,41`) and fail soft on unparseable content — there is no extension allow-list check beyond the `fileInput(accept=...)` browser hint (`mod_methyl_dataset.R:49,61`), which is a client-side suggestion only, not a server-side enforcement. |
| Empty files | Implemented | `methyl_parse_matrix()` checks `nrow(df)==0 \|\| ncol(df)<2` (`parse_upload.R:17`); `methyl_parse_sample_sheet()` checks `nrow(df)==0` (`parse_upload.R:42`); both return a clean `ok=FALSE` error rather than crashing. |
| Missing/required columns in sample sheet | Not implemented | `methyl_parse_sample_sheet()` performs **no** column-name validation at all (`parse_upload.R:40-46`) — any sheet with ≥1 row is accepted; a sheet with no sample-ID-like column silently falls through to `methyl_sheet_sample_ids()`'s row-order fallback (`qc.R:463-468`), which only works if row counts happen to match. |
| Non-numeric matrix values | Implemented | `methyl_parse_matrix()` explicitly attempts `storage.mode(m) <- "double"` in a `tryCatch` and returns `ok=FALSE, error="Every column after the first must be numeric."` on failure (`parse_upload.R:27-35`). |
| NAs in the matrix | Implemented (permitted, not rejected) | Recognized on read (`parse_upload.R:14`) and handled with `na.rm=TRUE`/`na.omit()` throughout `qc.R` — see §6. |
| Duplicate probe IDs | Implemented | `methyl_parse_matrix()` explicitly checks `any(duplicated(probe_ids))` and rejects with a count-specific error (`parse_upload.R:21-26`). |
| Duplicate sample IDs | Not implemented | No check anywhere for duplicate column names in the uploaded matrix, nor duplicate sample IDs in the sample sheet — confirmed by the absence of any `duplicated(colnames(...))`-style check in `parse_upload.R` or `mod_methyl_dataset.R`. |
| Invalid beta-value range (outside [0,1]) | Not implemented | See §6 — no bounds check exists anywhere in the upload or QC path. |
| Wrong matrix orientation (samples in rows) | Not implemented | No shape-sanity check (e.g. "probe IDs should look like `cgXXXXXXXX`") exists; a transposed upload would silently be treated as probes-in-rows and every downstream statistic would be computed on the wrong axis without error. |
| Insufficient samples | Partially implemented, per-tab, not at upload | No check at Dataset-tab load time. Individual QC tabs validate their own minimums at run time: Probe QC requires ≥3 samples in the current stratum (`mod_methyl_qc.R:722-723`), Outlier QC requires ≥4 (`mod_methyl_qc.R:1085-1086`), several structure functions (`methyl_pca_scores`, `methyl_sample_outliers_pca/hclust/mahalanobis`, `methyl_mds_scores`) independently require ≥10 complete-case probes and ≥4 samples (`qc.R:276-277,293-294,321-322,578-579,595,626-627`) and degrade with a reason rather than erroring. `methyl_qc_subgroup_filter()`/`methyl_apply_manual_exclude()` also set a `low_n` flag at `<3` samples, surfaced as a UI warning on the Sample QC tab (`mod_methyl_qc.R:535-536`) but **not** blocking any computation. |
| Insufficient CpGs | Not explicitly implemented | No minimum-probe-count check exists anywhere; an extremely small matrix would simply propagate through until a downstream function (e.g. `prcomp`) errors or degrades via its own `nrow(m) < 10` guard. |
| Missing metadata / group column | Implemented (graceful degradation) | Every metadata-dependent feature (subgroup stratification, sex comparison, batch column detection, RUVm factor-of-interest) explicitly checks for `is.null(sheet)`/column presence and shows an explanatory "not available" message rather than failing (e.g. `mod_methyl_qc.R:464,556-557,947-953,964-965`). |
| Array-type / annotation availability | Implemented | `methyl_get_annotation()` returns a specific `ok=FALSE` reason for unsupported array types (`annotation.R:50-55`), and every annotation-dependent UI element is conditionally hidden or shows an explanatory note for non-Illumina array types (`mod_methyl_qc.R:685-692`) or missing annotation packages. |
| IDAT pairing (Grn/Red) | Implemented | `methyl_read_idat()` explicitly checks both `_Grn.idat` and `_Red.idat` are present before calling `minfi::read.metharray.exp()` (`parse_upload.R:77-81`). |
| Path traversal in uploaded filenames | Implemented | `basename()` applied to every uploaded IDAT filename before it is used to construct a file path (`parse_upload.R:74-76`), explicitly called out as "path-traversal protection" in the code's own comment. |

---

## 8. Per-QC-analysis scientific documentation

For each analysis: what is measured, why it matters, how it is calculated in this codebase, whether the threshold is configurable, and — critically — whether the code makes an actual pass/fail decision or only visualizes/flags.

**Overview / basic QC pass.** *Measures:* median sample call rate, overall matrix missingness. *Why it matters (scientific background):* call rate and missingness are the most basic sanity checks on a methylation matrix before any further processing. *How calculated:* §3.1/§5.6. *Configurable:* No — both thresholds (0.90, 10%) are hard-coded inside `methyl_qc_status_badge()` (`qc.R:933,938`), not exposed as inputs. *Decision:* **Yes — this is the one analysis in the whole module that computes an explicit Pass/Warning/Fail verdict**, via `methyl_qc_status_badge()` (`qc.R:929-944`).

**Sample call rate (Sample QC).** *Measures:* per-sample fraction of non-missing probes. *How calculated:* `1 - colMeans(is.na(mat))` (`qc.R:192-194`). *Configurable:* `call_rate_min` (default 0.95, `mod_methyl_qc.R:541`). *Decision:* Flag-only — `call_rate_flag` is a boolean column in the displayed table (`mod_methyl_qc.R:582`); no sample is ever removed from `methyl_dataset` by this flag.

**Failed-probe percentage (Sample QC).** *Measures:* per-sample % of probes with `detection_p > threshold`. *Requires:* raw IDAT (`methyl_dataset$detp`). *Configurable:* `sample_detp_thresh` (default 0.01) and `failed_probe_pct_max` (default 5%, only applied if `f_failed_probe_pct` is checked, default off, `mod_methyl_qc.R:546-548`). *Decision:* Flag-only.

**Minimum signal intensity (Sample QC).** *Measures:* mean of median log2 methylated/unmethylated intensity per sample. *Requires:* raw IDAT. *Configurable:* `min_intensity_thresh` (default 10, only applied if `f_min_intensity` checked, default off, `mod_methyl_qc.R:551-553`). *Decision:* Flag-only.

**Bisulfite conversion efficiency (Sample QC).** *Measures:* % successful bisulfite conversion via `wateRmelon::bscon()`. *Requires:* raw IDAT. *Configurable:* No threshold exposed at all — the table is display-only, with no boolean flag column and no numeric cutoff input anywhere in the UI (confirmed: `bisulfite_table`, `mod_methyl_qc.R:647-653`, has no accompanying flag logic). *Decision:* **Visualization only — not even a flag**, unlike every other Sample QC metric.

**Probe filters (Probe QC — 11 filters).** *Measures:* per-probe detection failure, low bead count, missingness, SNP overlap, non-CpG design, sex-chromosome location, cross-reactivity (upload), population MAF (upload), variance, SD, mean-range. *How calculated:* §5.1. *Configurable:* every threshold has a numeric input with an explicit default (full table in §9); every filter can be individually enabled/disabled. *Decision:* **Yes, but the removal only applies to the exported/downloaded matrix (`probe_qc_result()$filtered`), never to `methyl_dataset$beta` itself** — the underlying shared dataset is never mutated (confirmed in §1/§12), so from the whole-app's perspective this is "compute and show what would be removed," even though within the tab's own scope the `filtered` matrix genuinely has fewer rows than `mat`.

**Good vs. bad Probe QC outcome (scientific background, contextualized to this code):** a "good" run typically retains the large majority of probes after standard filters (illustrated by the historical cascade's 412,492/485,577 ≈ 85% final retention, `global.R:300-304`) with each individual filter removing a modest fraction; a "bad" outcome — e.g. a filter removing a very large fraction of probes — is only visible by inspection of the filter-summary table/cascade plot; the code itself applies **no** its-own-outlier-check on the retention numbers (no "this filter removed suspiciously many probes" warning exists).

**Sex check.** *Measures:* predicted-vs-reported sex concordance. *How calculated:* §3.4/§5.2. *Configurable:* No threshold — deterministic method dispatch (raw-IDAT vs. beta-heuristic), no tunable input. *Decision:* Flag-only (`sex_mismatch` boolean per sample); exclusion requires an explicit separate button click (`mod_methyl_qc.R:896`).

**Batch correction (Batch QC — ComBat/RUVm).** *Measures:* whether batch/factor-driven structure is visible on PC1–PC2 before vs. after correction. *Configurable:* `batch_col` (ComBat, a categorical choice, not a threshold); `ruvm_group_col`, `ruvm_k` (default 1). *Decision:* **No pass/fail — visualization + a variance-explained table only**; the corrected matrix is neither auto-applied nor separately downloadable (Finding M-1, §14).

**Outlier detection (Outlier QC — 4 methods).** *Measures:* PCA-centroid distance, hierarchical-clustering singleton status, mean-correlation deviation, Mahalanobis distance². *Configurable:* `pca_sd` (default 3), `corr_k` (default 3), `mahal_alpha` (default 0.01); `hclust`'s `height_frac` is not exposed (Finding L-1). *Decision:* **Flag-only, by explicit design** ("Detection and removal are two separate, explicit steps", `mod_methyl_qc.R:1052-1054`) — samples are only removed via the separate "Apply Sample Exclusions" button, which writes to the shared `manual_exclude` list, not automatically from any detector's output.

**Cell-type-composition-related QC:** **Not applicable to this sub-module** — cell-type deconvolution (e.g. Houseman-style reference-based estimation) is not implemented anywhere in `mod_methyl_qc.R`/`qc.R`; it is out of scope for this document (a separate, registry-listed `mod_methyl_celltype` sub-module exists per `submodules_registry.R:42` but was not audited here, per the task's file scope).

---

## 9. Parameter inventory table

Every `numericInput`/`checkboxInput`/`selectInput`/`radioButtons`/`checkboxGroupInput` in `mod_methyl_qc.R`, tab by tab.

### Overview
| Parameter | Default | Allowed values | Function affected | Effect | Output affected |
|---|---|---|---|---|---|
| `live_group_col` (`mod_methyl_qc.R:462-463`) | auto-detected sex column, else `""` | any `sample_sheet` column, or "All samples" | `methyl_qc_subgroup_filter()` | selects the stratification column read by every tab | `current_subgroup()`, hence every tab |
| `live_stratum` (`mod_methyl_qc.R:319-335`) | `"__all__"` | levels of `live_group_col`, dynamically populated | `methyl_qc_subgroup_filter()` | selects the stratum level | `current_subgroup()`, hence every tab |

### Sample QC
| Parameter | Default | Allowed values | Function affected | Effect | Output affected |
|---|---|---|---|---|---|
| `call_rate_min` | 0.95 | 0–1, step 0.01 | flag logic inline (`mod_methyl_qc.R:582`) | sets `call_rate_flag` threshold | `sample_qc_table` |
| `sample_detp_thresh` | 0.01 | 0–1, step 0.01 | `methyl_sample_failed_probe_pct()` | detection-p cutoff for failed-probe % | `sample_qc_table` |
| `f_failed_probe_pct` | FALSE | boolean | flag logic inline | enables `failed_probe_pct_flag` column | `sample_qc_table` |
| `failed_probe_pct_max` | 5 | 0–100, step 1 | flag logic inline | sets the failed-probe % threshold | `sample_qc_table` |
| `f_min_intensity` | FALSE | boolean | flag logic inline | enables `low_intensity_flag` column | `sample_qc_table` |
| `min_intensity_thresh` | 10 | numeric, step 0.5 | `methyl_sample_low_intensity()` | sets the median-log2-intensity threshold | `sample_qc_table` |

### Probe QC
| Parameter | Default | Allowed values | Function affected | Effect | Output affected |
|---|---|---|---|---|---|
| `f_detp` | TRUE iff IDAT present | boolean | `methyl_filter_detection_p()` | enable/disable filter | `filtered`, `cascade`, `filter_table` |
| `detp_thresh` | 0.01 | 0–1, step 0.01 | `methyl_filter_detection_p()` | detection-p cutoff | as above |
| `f_beadcount` | TRUE iff IDAT+beadcount present | boolean | `methyl_filter_beadcount()` | enable/disable filter | as above |
| `beadcount_thresh` | 3 | ≥1, step 1 | `methyl_filter_beadcount()` | min bead count | as above |
| `f_missing` | TRUE | boolean | `methyl_filter_missing()` | enable/disable filter | as above |
| `missing_max` | 0 | 0–1, step 0.05 | `methyl_filter_missing()` | max allowed missing fraction | as above |
| `f_snp` | TRUE (Illumina only) | boolean | `methyl_filter_snp()` | enable/disable filter | as above |
| `f_noncpg` | TRUE (Illumina only) | boolean | `methyl_filter_non_cpg()` | enable/disable filter | as above |
| `sexchr_mode` | `"keep"` | keep / remove_xy / remove_y_only | `methyl_filter_sex_chr()` | sex-chromosome handling | as above |
| `f_crossreactive` | FALSE | boolean | `methyl_filter_cross_reactive()` | enable/disable filter | as above |
| `crossreactive_file` | none | uploaded CSV/TXT/TSV | `methyl_parse_probe_list()` | supplies exclusion IDs | as above |
| `f_maf` | FALSE | boolean | `methyl_filter_maf()` | enable/disable filter | as above |
| `maf_file` | none | uploaded CSV/TXT/TSV | `methyl_parse_maf_list()` | supplies MAF table | as above |
| `maf_max` | 0.05 | 0–0.5, step 0.01 | `methyl_filter_maf()` | max allowed MAF | as above |
| `f_variance` | FALSE | boolean | `methyl_filter_variance()` | enable/disable filter | as above |
| `variance_min` | 0 | ≥0, step 0.001 | `methyl_filter_variance()` | min variance | as above |
| `f_sd` | FALSE | boolean | `methyl_filter_sd()` | enable/disable filter | as above |
| `sd_min` | 0 | ≥0, step 0.01 | `methyl_filter_sd()` | min SD | as above |
| `f_meanrange` | FALSE | boolean | `methyl_filter_mean_range()` | enable/disable filter | as above |
| `mean_lo`/`mean_hi` | 0.01/0.99 (beta) or −6/6 (M) | numeric | `methyl_filter_mean_range()` | mean-value bounds | as above |

### Sex QC
No tunable parameters — one Run button only (`mod_methyl_qc.R:841`).

### Batch QC
| Parameter | Default | Allowed values | Function affected | Effect | Output affected |
|---|---|---|---|---|---|
| `batch_method` | "combat" if batch col found, else "ruvm" | combat / ruvm | dispatch (`mod_methyl_qc.R:986,990`) | selects correction method | `batch_qc_result()` |
| `batch_col` | first detected batch column | any `sample_sheet` column | `methyl_batch_correct_combat()` | supplies batch labels | as above |
| `ruvm_group_col` | first `sample_sheet` column | any `sample_sheet` column | `methyl_batch_correct_ruvm()` | factor of interest to protect | as above |
| `ruvm_k` | 1 | 1–10, step 1 | `methyl_batch_correct_ruvm()` (`RUVfit(k=...)`) | # unwanted-variation factors | as above |

### Outlier QC
| Parameter | Default | Allowed values | Function affected | Effect | Output affected |
|---|---|---|---|---|---|
| `outlier_methods` | `c("pca","hclust")` | subset of pca/hclust/correlation/mahalanobis | dispatch (`mod_methyl_qc.R:1094,1098,1101,1104`) | which flag columns are shown | `outlier_score_table` |
| `pca_sd` | 3 | ≥1, step 0.5 | `methyl_sample_outliers_pca()` | SD-distance threshold | PCA flag, plots |
| `corr_k` | 3 | ≥1, step 0.5 | `methyl_sample_outliers_correlation()` | MAD multiplier | correlation flag |
| `mahal_alpha` | 0.01 | 0.001–0.2, step 0.001 | `methyl_sample_outliers_mahalanobis()` | chi-squared alpha | Mahalanobis flag |

### Visualizations
| Parameter | Default | Allowed values | Function affected | Effect | Output affected |
|---|---|---|---|---|---|
| `viz_color_by` | `"__none__"` | any `sample_sheet` column | `.viz_color_vec()` (`mod_methyl_qc.R:1289-1299`) | point coloring | PCA 2D/3D, MDS plots only |

### Reports & Export
No tunable parameters — export buttons only.

---

## 10. Reactive-logic documentation

**The has-run-gate pattern.** Each of the six computable tabs (Overview, Sample QC, Probe QC, Sex QC, Batch QC, Outlier QC) follows an identical three-piece structure:
1. A `reactiveVal(FALSE)` flag, e.g. `overview_has_run` (`mod_methyl_qc.R:281-286`).
2. An `eventReactive(input$run_*_btn, {...})` that performs the actual computation, e.g. `overview_result` (`mod_methyl_qc.R:493-502`) — this is what "button-driven" means concretely: the reactive body only re-executes when its bound `actionButton`'s value changes, ignoring every other reactive dependency read inside it (Shiny's `eventReactive` semantics).
3. An `observeEvent(result(), has_run_flag(TRUE))` (e.g. `mod_methyl_qc.R:503`) that flips the flag once the computation completes.

**Why a separate gate output exists.** `register_has_run_gate(gate_id, has_run_flag_fn, result_output_id, not_run_message)` (`qc.R`... — actually `mod_methyl_qc.R:124-129`) creates a small, independent `renderUI` per tab that reads *only* the `has_run` flag and decides between a placeholder message and `uiOutput(result_output_id)`. This exists specifically to fix a documented prior performance bug: previously each tab's *entire* `renderUI` (controls, inputs, and the just-clicked button itself) read the `has_run` flag inline, so clicking Run tore down and rebuilt the whole tab, "which is slow, and can eat the very click that triggered it if the browser replaces the button's DOM node while the click is still being processed" (`mod_methyl_qc.R:96-109`). After the fix, a tab's controls-and-button `renderUI` (e.g. `output$overview_ui`) never reads its own `has_run` flag at all — only the tiny gate output does.

**Automatic vs. button-driven vs. upload-triggered — a full accounting:**
- **Never automatic:** none of the six tabs' `eventReactive`s fire on dataset load, parameter change, or stratum change — confirmed per-tab in §3, and matching the file's own header claim verbatim (`mod_methyl_qc.R:16-31`).
- **Button-driven:** all six tabs' main computations, gated on their own `actionButton` only.
- **Upload-triggered (but not QC-computation-triggered):** loading a dataset (any of the three paths on the Dataset tab) only ever assigns `methyl_dataset$...` fields (`mod_methyl_dataset.R:91-100,207-218,277-286`); the one QC-side reaction to this is the `observeEvent(methyl_dataset$beta, {...})` at `mod_methyl_qc.R:288-292`, which **resets** (does not compute) `manual_exclude` and all six `*_has_run` flags back to their empty/`FALSE` state, "UI hygiene, not an auto-run: nothing is recomputed by that reset, results just stop being shown until re-run" (`mod_methyl_qc.R:29-31`).
- **Always-live, ungated reactives (data scoping, not analysis):** `current_subgroup()` (`mod_methyl_qc.R:342-346`), `current_rg_subset()`/`current_mset_subset()` (`mod_methyl_qc.R:347-356`), `stratum_all_samples()` (`mod_methyl_qc.R:363-367`), `anno_result()`/`is_illumina_array()` (`mod_methyl_qc.R:260-268`) — all plain `reactive()`s, recomputed on every relevant input change, deliberately not gated behind any button because they are "cheap (index/character-vector work, no real computation)" (`mod_methyl_qc.R:339-341`).

**Lazy plot rendering (Visualizations, and the Outlier QC plots).** A second, independent gating mechanism, `lazy_plot_ui()`/`plot_shown` (`mod_methyl_qc.R:131-168`), controls 14 individual plots (`lazy_plot_ids`, `mod_methyl_qc.R:132-134`) across the Outlier QC and Visualizations tabs. Each plot's "Generate" button and its (initially `shinyjs::hidden()`) output container are built together, once, unconditionally (`mod_methyl_qc.R:156-168`); clicking Generate only flips `plot_shown[[pid]] <- TRUE` and toggles CSS visibility via `shinyjs::show()`/`hide()` (`mod_methyl_qc.R:138-141`) — the plot's own `render*` block is separately gated via `req(isTRUE(plot_shown[[pid]]))` (e.g. `mod_methyl_qc.R:1163,1302`). A `shinyjs::delay(200, ...Plotly.Plots.resize...)` nudge (`mod_methyl_qc.R:146-149`) works around `plotly` rendering blank inside a container that was `display:none` at creation time. This is UI/rendering plumbing, not a QC method, per the task's framing — documented here for completeness, not analyzed as a statistical procedure.

---

## 11. Tab-to-tab relationship section

Restating and formalizing the finding from §4: **eight tabs, effectively one genuine cross-tab data dependency.**

```
             ┌─────────────┐
             │  Overview   │  (self-contained; reads current_subgroup() only)
             └─────────────┘
             ┌─────────────┐        writes ──┐
             │ Sample QC   │  (self-contained)│
             └─────────────┘                  │
             ┌─────────────┐        ══════▶ manual_exclude()  ◀── read by
             │ Probe QC    │  (self-contained)│  every tab via current_subgroup()
             └──────┬──────┘                  │
                    │ $filtered               │
                    ▼                         │
             ┌─────────────┐                  │
             │Visualizations│                 │
             └─────────────┘                  │
             ┌─────────────┐        writes ──┤
             │  Sex QC     │  (self-contained, exclude button)
             └─────────────┘                  │
             ┌─────────────┐                  │
             │ Batch QC    │  (self-contained)│
             └─────────────┘                  │
             ┌─────────────┐        writes ──┘
             │ Outlier QC  │  (self-contained, exclude button)
             └─────────────┘
             ┌─────────────┐
             │Reports&Export│ (reads ALL SIX *_result() objects, read-only sink)
             └─────────────┘
```

**Verified empirically, not assumed:**
- `current_subgroup()` (`mod_methyl_qc.R:342-346`) is read by every tab's `eventReactive`, but it is itself built purely from `methyl_dataset$beta`/`$sample_sheet` and the shared UI inputs (`live_group_col`/`live_stratum`) plus `manual_exclude()` — never from another tab's *result*.
- `manual_exclude()` (`reactiveVal`, `mod_methyl_qc.R:280`) is written by three independent UI actions — Sample QC's "Apply exclusions" (`mod_methyl_qc.R:379-382`), Sex QC's "Exclude selected discordant samples" (`mod_methyl_qc.R:929-935`), Outlier QC's "Apply Sample Exclusions" (`mod_methyl_qc.R:1155-1160`) — each via `manual_exclude(union(manual_exclude(), sel))`, i.e. exclusions accumulate across all three mechanisms into one shared set, never overwritten by one tab's action alone.
- No `eventReactive` in the file reads another tab's `*_result()` **except** Visualizations reading `probe_qc_result()` (§3.3/§3.7) and Reports & Export reading all six via `current_qc_pieces()` (§3.8) — both confirmed by direct reading of every `eventReactive`/`reactive`/`render*` block in the 1,619-line file, not by trusting the header comment.
- The header comment's own framing — "The one deliberate exception is that every method reads the SAME `manual_exclude()` sample-exclusion set, since that is data scoping... not a QC method's output" (`mod_methyl_qc.R:36-39`) — is accurate as written; this document independently confirmed it rather than merely restating it.

---

## 12. Output documentation

For each output type produced by the module, whether it is visualization-only or data-transforming, and whether anything outside `mod_methyl_qc.R` consumes it.

| Output | Type/dimensions | Source | Visualization-only or data-transforming | Consumed downstream (outside this file)? |
|---|---|---|---|---|
| `overview_result()` | list (scalars + subgroup) | Overview eventReactive | Data-transforming (computes summary stats) but result lives only in-tab | No — only feeds `current_qc_pieces()$overview` for export |
| `sample_qc_result()` | list incl. data.frame `sample_qc` | Sample QC eventReactive | Data-transforming | No |
| `probe_qc_result()` | list incl. matrices `mat`/`filtered` | Probe QC eventReactive | Data-transforming — `filtered` is a genuinely reduced-row matrix | Yes, but only *within this file*: Visualizations tab (§3.3/3.7) and Reports & Export |
| `sex_qc_result()` | list incl. data.frame `sex$detail` | Sex QC eventReactive | Data-transforming (predictions + concordance flag) | No |
| `batch_qc_result()` | list incl. `out$corrected` matrix | Batch QC eventReactive | Data-transforming (produces a corrected matrix) | **No — `out$corrected` is never downloaded or read by any other tab** (Finding M-1) |
| `outlier_qc_result()` | list incl. data.frame `outlier_scores` | Outlier QC eventReactive | Data-transforming (flags) | No |
| Filtered beta/M-value CSV downloads | file | `dl_filtered_beta`/`dl_filtered_mvalue` | Data export | User-driven only, outside the app |
| QC report HTML/PDF | file | `dl_report_html`/`dl_report_pdf` | Static export (baked-in PNGs, not interactive) | User-driven only |
| Figures ZIP | file | `dl_figures_zip` | Static export | User-driven only |
| `methyl_results` (shared reactiveValues) | n/a | — | — | **Never written** — confirmed by exhaustive grep; the parameter is accepted (`mod_methyl_qc.R:71`) but unused throughout the file |
| `methyl_dataset` (shared reactiveValues) | n/a | — | — | **Never mutated by this module** — confirmed by exhaustive grep for `methyl_dataset$... <-` assignments in `mod_methyl_qc.R`: none exist |

**Downstream-module check.** Grepped `R/methylomics/*.R` and `server.R` for any reference to `probe_qc_result`, `sample_qc_result`, `sex_qc_result`, `batch_qc_result`, `outlier_qc_result`, or `overview_result` outside `mod_methyl_qc.R` itself: none found. The header comment's claim — QC "never removes probes or samples from `methyl_dataset` itself; it reports what each filter WOULD remove and exposes the resulting filtered matrix for download" (`mod_methyl_qc.R:6-9`) — is fully confirmed: every `eventReactive` in the file was individually inspected in §3, and none contains a `methyl_dataset$... <-` assignment.

---

## 13. Scientific and computational audit

For each area: Implemented / Not implemented / Not applicable / Potential enhancement, strictly evidence-based.

| Area | Status | Evidence |
|---|---|---|
| Sample call-rate QC | Implemented | `methyl_sample_call_rate()`, `qc.R:192-194` |
| Detection p-value probe filtering | Implemented (IDAT-only) | `methyl_filter_detection_p()`, `qc.R:162-172` |
| Detection p-value sample-level failed-probe % | Implemented (IDAT-only) | `methyl_sample_failed_probe_pct()`, `qc.R:202-217` |
| Bead-count probe filtering | Implemented (IDAT-only, stricter than ChAMP default) | `methyl_filter_beadcount()`, `qc.R:178-188` |
| Bisulfite conversion efficiency | Implemented (IDAT-only, visualization only, no flag) | `methyl_bisulfite_conversion()`, `idat_metrics.R:33-40` |
| Missingness filtering (probe-level) | Implemented | `methyl_filter_missing()`, `qc.R:30-34` |
| Beta/M-value distribution QC (density/boxplot/violin) | Implemented | `methyl_plot_density/boxplot/violin()`, `qc.R:851-877` |
| Outlier detection — PCA | Implemented (explicitly non-rigorous, disclosed) | `methyl_sample_outliers_pca()`, `qc.R:274-286` |
| Outlier detection — hierarchical clustering | Implemented | `methyl_sample_outliers_hclust()`, `qc.R:291-304` |
| Outlier detection — correlation-based | Implemented | `methyl_sample_outliers_correlation()`, `qc.R:255-263` |
| Outlier detection — Mahalanobis (robust/MCD variant) | Not implemented (classical variant only, disclosed) | `methyl_sample_outliers_mahalanobis()` uses `stats::cov()`, `qc.R:332`; no MCD estimator anywhere in the file |
| Outlier detection — IQR-based univariate | Not implemented (dead code) | `methyl_sample_outliers_iqr()` defined (`qc.R:234-239`) but never called |
| PCA sample-structure visualization | Implemented | `methyl_pca_scores()`, `qc.R:576-587`; PCA 2D/3D plots |
| MDS sample-structure visualization | Implemented | `methyl_mds_scores()`, `qc.R:624-637` |
| Sex check — raw-intensity (copy-number) method | Implemented (IDAT-only) | `minfi::getSex()` via `methyl_sex_check()`, `qc.R:408-419` |
| Sex check — beta-based fallback | Implemented | `methyl_cluster_sex()`, `qc.R:356-376` |
| SNP-overlapping probe removal (manifest-based) | Implemented (450K/EPIC only) | `methyl_filter_snp()`, `qc.R:68-81` |
| Cross-reactive probe removal (published blacklist) | Not implemented (deliberately, upload-only) | `methyl_filter_cross_reactive()`, `qc.R:112-119`, no bundled list |
| Population-MAF-based probe filtering | Not implemented (deliberately, upload-only) | `methyl_filter_maf()`, `qc.R:149-160`, no bundled table |
| Non-CpG (CpH) probe removal | Implemented (Illumina arrays only) | `methyl_filter_non_cpg()`, `qc.R:63-66` |
| Sex-chromosome probe removal | Implemented, three modes | `methyl_filter_sex_chr()`, `qc.R:89-104` |
| Batch correction — ComBat | Implemented | `methyl_batch_correct_combat()`, `qc.R:723-750` |
| Batch correction — RUVm | Implemented (IDAT-only) | `methyl_batch_correct_ruvm()`, `qc.R:764-808` |
| Normalization | Not applicable to this module | `idat_metrics.R:9-10` explicitly states "no normalization — normalization methods live in a later Methylomics sub-module"; confirmed no normalization function exists in `qc.R` |
| Cell-type deconvolution / composition QC | Not applicable to this module | Not present in `mod_methyl_qc.R`/`qc.R`; a separate registry entry (`mod_methyl_celltype`) exists but is out of this document's scope |
| Beta-value range validation | Not implemented | See §6/§7 |
| Explicit pass/fail decision-making | Implemented, narrowly (Overview only) | `methyl_qc_status_badge()`, `qc.R:929-944` |
| Automated probe/sample removal from the shared dataset | Not implemented (by design) | Confirmed in §12 — the module only ever reports/exports, never mutates `methyl_dataset` |
| Robust (MCD-based) multivariate outlier detection | Potential enhancement | Explicitly flagged as a known limitation in the code's own comment, `qc.R:306-318` |
| Bundled cross-reactive/MAF probe lists | Potential enhancement (deliberately deferred) | `qc.R:106-111,126-130` — the code explains why this was not done rather than treating it as an oversight |

---

## 14. Code-level audit

Findings only where genuinely evidence-based; no finding is manufactured. Severities: Critical / High / Moderate / Low / Informational.

**Finding C-1 (Informational) — `methyl_results` parameter accepted but never used.**
*Evidence:* `mod_methyl_qc_server <- function(id, methyl_dataset, methyl_results)` (`mod_methyl_qc.R:71`); exhaustive grep for `methyl_results` in the file returns only this one line.
*Why it matters:* Identical pattern to the documented `mod_interaction.R` `results` parameter (house-style reference, §2 item 28) — a shared results store other Methylomics sub-modules could theoretically read from is passed in but this module writes nothing to it, so no cross-sub-module handoff of QC results exists via this mechanism.
*Consequence:* None functionally (the parameter being unused causes no bug); purely a documentation/API-clarity note.
*Recommended action:* None required; flagged for completeness only.

**Finding C-2 (Low) — Overview's pass/fail thresholds are hard-coded, not user-configurable.**
*Evidence:* `if (!is.na(cr) && cr < 0.90)` (`qc.R:933`); `if (!is.na(miss) && miss > 10 ...)` (`qc.R:938`) — neither `0.90` nor `10` is passed in from a Shiny input; both are literals inside `methyl_qc_status_badge()`.
*Why it matters:* A user cannot adjust what counts as a "Fail" vs. "Warning" on the Overview tab, unlike every threshold on Sample QC/Probe QC/Outlier QC, which are all exposed as `numericInput`s.
*Consequence:* Minor inconsistency in the module's otherwise pervasive "everything is a tunable input" design; not a correctness bug.
*Recommended action:* None mandated by this audit — noted as a design-consistency observation.

**Finding H-1 (High) — Batch correction functions are never guarded by an `input_scale` check.**
*Evidence:* `methyl_batch_correct_combat()` (`qc.R:723-750`) and `methyl_batch_correct_ruvm()` (`qc.R:764-808`) both unconditionally call `methyl_beta_to_mvalue(mat)` (`qc.R:737,783`), which clips its input to `[1e-4, 1-1e-4]` and applies `log2(b/(1-b))` (`qc.R:566-569`) — an operation that assumes `mat` is already on the `[0,1]` beta scale. Neither `batch_qc_result` (`mod_methyl_qc.R:974-994`) nor either correction function checks `methyl_dataset$input_scale` before calling this.
*Why it matters:* If a user uploaded an M-value matrix (`input_scale == "m"`), M-values (which are unbounded, typically roughly in `[-6,6]`) would be silently clamped into `[1e-4, 1-1e-4]` by `pmin(pmax(beta, eps), 1-eps)` (`qc.R:567`) — a nonsensical transform that would corrupt the corrected output without any warning or error.
*Consequence:* Silently wrong batch-corrected values for any M-value-scale dataset run through Batch QC; the resulting `out$corrected` matrix (even though not currently downloadable, per Finding M-1) would also corrupt the "after" PCA plot's apparent batch-correction quality.
*Recommended action:* Guard `methyl_batch_correct_combat`/`methyl_batch_correct_ruvm` (or their caller at `mod_methyl_qc.R:982-991`) with an explicit `input_scale` check, converting from M-value to beta first if needed, or refusing with a clear message.

**Finding H-2 (High) — No beta-value range validation anywhere in the pipeline.**
*Evidence:* See §6/§7 — `methyl_parse_matrix()` (`parse_upload.R:11-37`) validates only numeric-ness and duplicate IDs, never a `[0,1]` bounds check; no other function in `qc.R`/`mod_methyl_qc.R` performs one either.
*Why it matters:* A user who selects "Beta values (0-1)" (`mod_methyl_dataset.R:48`) but actually uploads M-values (or an already-standardized/z-scored matrix) would have every beta-scale-assuming computation in the module (mean-range filter defaults, `methyl_beta_to_mvalue()`'s clipping, the density-plot x-axis label) silently produce misleading output with no warning.
*Consequence:* Same class of risk as Finding H-1, but triggerable purely by a mislabeled upload rather than a genuine M-value dataset routed into Batch QC specifically.
*Recommended action:* Add a lightweight range/sanity check at matrix-load time (e.g. warn if `input_scale=="beta"` but `range(mat, na.rm=TRUE)` falls well outside `[0,1]`), surfaced as a non-blocking warning consistent with the module's fail-soft conventions elsewhere.

**Finding I-1 (Informational) — Bead-count filter is deliberately stricter than ChAMP's default.**
*Evidence:* `qc.R:174-188`, quoted in full in §5.1 — removes a probe if it fails in **any** sample, vs. ChAMP's `champ.filter(beadCutoff=0.05)` (>5% of samples).
*Why it matters:* A user familiar with ChAMP's numeric output may be surprised this tool removes more probes at the "same" bead-count threshold.
*Consequence:* None — this is a disclosed, intentional design choice documented in the code's own comment, not an unintentional bug.
*Recommended action:* None required by this audit; the in-app note already explains the ChAMP divergence to the user only indirectly (via the generated R-code snippet, `qc.R:959`), not directly in the Probe QC tab's own UI text — a minor UX-clarity opportunity, not a correctness issue.

**Finding L-1 (Low) — `hclust` outlier detection's `height_frac` parameter is not exposed as a tunable input.**
*Evidence:* `methyl_sample_outliers_hclust(mat, n_features=5000, height_frac=0.5)` (`qc.R:291`) is called at `mod_methyl_qc.R:1097` with no `height_frac` argument, so the default `0.5` is always used; no corresponding `numericInput` exists in the Outlier QC UI (`mod_methyl_qc.R:1056-1081`), unlike `pca_sd`/`corr_k`/`mahal_alpha`, which are all exposed.
*Why it matters:* Inconsistent with the module's otherwise pervasive "every threshold is a numericInput" pattern for the other three outlier methods.
*Consequence:* A user cannot tune how aggressively the dendrogram-singleton rule flags samples.
*Recommended action:* None mandated; noted as a minor UI-completeness gap.

**Finding L-2 (Low/Performance) — All four outlier detectors always run, regardless of which methods are selected.**
*Evidence:* `mod_methyl_qc.R:1093,1097,1100,1103` call `methyl_sample_outliers_pca/hclust/correlation/mahalanobis` unconditionally inside `outlier_qc_result`; the `methods_sel` checkbox group (`mod_methyl_qc.R:1094,1098,1101,1104`) only gates whether each result's flag column is attached to the displayed table, not whether the underlying computation executes.
*Why it matters:* Selecting only "PCA-based detection" still pays the full computational cost of hierarchical clustering, correlation-matrix computation, and Mahalanobis-distance PCA — each of which independently repeats a top-variance-probe selection and (for PCA/Mahalanobis) a separate `prcomp()` call, rather than sharing one PCA computation across methods that both need it.
*Consequence:* Unnecessary compute time on every Outlier QC run, worse as sample/probe count grows; no incorrect output results from this, purely a performance inefficiency.
*Recommended action:* Gate each detector call on `method %in% methods_sel` before computing, and/or factor out one shared `prcomp()` call for PCA-outlier and Mahalanobis (both currently call `methyl_pca_scores()`-equivalent logic independently, `qc.R:280-281` and `qc.R:325-326`).

**Finding L-3 (Low) — `methyl_filter_variance` and `methyl_filter_sd` are mathematically redundant.**
*Evidence:* `qc.R:36-52`, both documented in §5.1; `sd = sqrt(variance)`, so `sd >= min_sd` and `variance >= min_sd^2` are the same cut expressed in different units. Both are independently offered as separate checkboxes in the Probe QC UI (`mod_methyl_qc.R:700-705`).
*Why it matters:* A user could enable both filters simultaneously with mismatched thresholds under the mistaken impression they measure different things, without realizing they're the same statistic in different units — the code's own comment acknowledges this ("same computation... offered as its own independent optional filter", `qc.R:43-46`) as a deliberate usability choice, not an accident.
*Consequence:* No incorrect output (both filters compute correctly), only a possible source of user confusion.
*Recommended action:* None mandated — already a disclosed, intentional design choice.

**Finding L-4 (Low) — `methyl_sample_outliers_iqr()` is dead code.**
*Evidence:* Defined at `qc.R:234-239`; exhaustive grep across `R/methylomics/*.R` finds no caller.
*Consequence:* No functional impact; `methyl_outlier_score_table()`'s `flag_cols` (`qc.R:530`) anticipates an `iqr_outlier` column that can never be populated.
*Recommended action:* Either wire it into the Outlier QC checkbox group as a fifth method, or remove it as unused.

**Finding L-5 (Low) — `methyl_guess_batch_column()` is defined but its logic is duplicated inline instead of called.**
*Evidence:* `qc.R:704-709` defines the function; `mod_methyl_qc.R:944-945` re-implements the identical regex (`"batch|chip|plate|slide|sentrix|array_id|scan_date|^run$"`, `ignore.case=TRUE`) inline via `grep()` rather than calling `methyl_guess_batch_column(sheet)`.
*Consequence:* Pure code duplication; both copies are currently in sync, but a future edit to one regex without the other would silently desynchronize Overview's batch-column detection note (`mod_methyl_qc.R:448-449`, which also duplicates the same pattern a third time) from Batch QC's actual gating logic.
*Recommended action:* Replace both inline `grep()` calls with `methyl_guess_batch_column(sheet)`.

**Finding M-1 (Moderate) — Batch-corrected matrix is never downloadable or reusable.**
*Evidence:* `batch_qc_result()$out$corrected` (`qc.R:749,807`) is read only by `.batch_pca_plot()`/`methyl_pca_scores()` for the "after" PCA plot (`mod_methyl_qc.R:992,1034-1036`) and by `batch_variance_table` (`mod_methyl_qc.R:1041-1049`); no `downloadHandler` anywhere in the Reports & Export tab (§3.8's download table) offers the corrected matrix, unlike Probe QC's `filtered` matrix, which has two dedicated CSV downloads.
*Why it matters:* A user who runs Batch QC and is satisfied with the correction has no way to export the corrected values from the app — only Probe QC's filtered matrix is exportable.
*Consequence:* Batch QC is effectively a diagnostic-only tab in practice (visualize whether correction would help) despite computing a full corrected matrix internally.
*Recommended action:* Add a `dl_batch_corrected` download handler analogous to `dl_filtered_beta`, gated on `batch_qc_has_run()`.

**No findings of Critical severity** were identified — no evidence of data corruption affecting the *displayed* results under normal (correctly-labeled beta-scale) usage, no crash-inducing code path found in the sections read, and the module's core "never mutates the shared dataset" claim was independently verified true.

**Severity breakdown:** Critical: 0. High: 2 (H-1, H-2). Moderate: 1 (M-1). Low: 4 (L-1 through L-5, five items — see note). Informational: 2 (C-1, I-1). Consistency note (C-2): 1 (Low). Total distinct findings: 10.

---

## 15. End-to-end user workflow narrative

Implemented steps only, in the order a user would actually encounter them:

1. **Dataset selection** (Dataset tab, `mod_methyl_dataset.R`, outside this module): choose preloaded / matrix upload / IDAT upload; click the corresponding load button. `methyl_dataset` fields populate; no QC computation is triggered by this step (`mod_methyl_qc.R:288-292` only resets flags).
2. **Navigate to Quality Control** (`mx_qc`, registered `mod_methyl_qc_config`, `mod_methyl_qc.R:40-43`). `output$body_ui` renders the 8-tab `tabsetPanel` once `methyl_dataset$beta` is non-`NULL` (`mod_methyl_qc.R:408-436`), else a prompt to load data first (`mod_methyl_qc.R:409-413`).
3. **Overview tab (opens first):** review dataset facts (probe/sample counts, scale); optionally pick a `live_group_col`/`live_stratum` to restrict every subsequent tab to a subgroup (`mod_methyl_qc.R:456-467`); click "Run Overview QC" for the basic call-rate/missingness pass/fail badge.
4. **Sample QC:** set `call_rate_min` and optional failed-probe-%/min-intensity filters; click "Run Sample QC"; review flags; optionally select rows in the manual-inclusion/exclusion table and click "Apply exclusions" — this narrows the sample scope every tab reads on its *next* run.
5. **Probe QC:** enable/disable and configure up to 11 filters; click "Run Probe QC"; review the filter-summary table and retention cascade; optionally download the filtered beta/M-value matrix.
6. **Sex QC:** click "Run Sex QC"; review predicted-vs-reported sex concordance and the X-vs-Y scatter; optionally select discordant samples and click "Exclude selected discordant samples".
7. **Batch QC:** choose ComBat or RUVm and its inputs; click "Run Batch QC"; compare before/after PCA plots and the variance-explained table (no export of the corrected matrix is currently possible, per Finding M-1).
8. **Outlier QC:** select one or more detection methods and their thresholds; click "Run Outlier Detection"; review the ranked outlier-score table; optionally select rows and click "Apply Sample Exclusions".
9. **Visualizations:** if Probe QC has been run, click individual "Generate" buttons for PCA/MDS/density/boxplot/violin/correlation-heatmap/mean-SD plots (optionally colored by a metadata column); if IDAT was uploaded, also generate detection-p, bead-count, and control-probe plots.
10. **Reports & Export:** download whichever CSVs correspond to tabs already run; download the self-contained HTML report (always available) or PDF report (only if a LaTeX toolchain is present); download the all-figures ZIP; review the Probe QC reproducibility table and copy the equivalent Bioconductor R-code snippet.

At every step, re-running an earlier tab after excluding samples on a later tab **does not automatically happen** — a user must revisit and re-click each tab's own Run button to see updated results against the current sample scope, per the module's explicit no-auto-recompute design (§10).

---

## 16. Function inventory table

| Function | File | Type | Used by | Purpose | Input | Output |
|---|---|---|---|---|---|---|
| `methyl_row_vars` | qc.R:22-28 | probe filter primitive | 9+ call sites | fast row variance | matrix | numeric vector |
| `methyl_filter_missing` | qc.R:30-34 | probe filter | Probe QC | drop high-missingness probes | matrix, `max_na_frac` | keep/note list |
| `methyl_filter_variance` | qc.R:36-41 | probe filter | Probe QC | drop low-variance probes | matrix, `min_variance` | keep/note list |
| `methyl_filter_sd` | qc.R:47-52 | probe filter | Probe QC | drop low-SD probes | matrix, `min_sd` | keep/note list |
| `methyl_filter_mean_range` | qc.R:54-58 | probe filter | Probe QC | drop out-of-range-mean probes | matrix, lo, hi | keep/note list |
| `methyl_filter_non_cpg` | qc.R:63-66 | probe filter | Probe QC | drop non-CpG (CpH) probes | matrix | keep/note list |
| `methyl_filter_snp` | qc.R:68-81 | probe filter | Probe QC | drop SNP-overlapping probes | matrix, annotation | keep/note list |
| `methyl_filter_sex_chr` | qc.R:89-104 | probe filter | Probe QC | drop chrX/chrY probes | matrix, annotation, mode | keep/note list |
| `methyl_filter_cross_reactive` | qc.R:112-119 | probe filter | Probe QC | drop uploaded exclusion-list probes | matrix, IDs | keep/note list |
| `methyl_parse_maf_list` | qc.R:131-144 | upload parser | Probe QC | parse probe_id,maf table | file path | ok/maf/error |
| `methyl_filter_maf` | qc.R:149-160 | probe filter | Probe QC | drop high-MAF probes | matrix, maf table, max | keep/note list |
| `methyl_filter_detection_p` | qc.R:162-172 | probe filter | Probe QC | drop failed-detection probes | matrix, detp, threshold | keep/note list |
| `methyl_filter_beadcount` | qc.R:178-188 | probe filter | Probe QC | drop low-bead-count probes | matrix, beadcount, threshold | keep/note list |
| `methyl_sample_call_rate` | qc.R:192-194 | sample metric | Overview, Sample QC | per-sample call rate | matrix | numeric vector |
| `methyl_sample_failed_probe_pct` | qc.R:202-217 | sample metric | Sample QC | per-sample failed-probe % | matrix, detp, threshold | ok/pct list |
| `methyl_sample_low_intensity` | qc.R:225-232 | sample metric | Sample QC | flag low-intensity samples | median-intensity result, threshold | ok/score/low list |
| `methyl_sample_outliers_iqr` | qc.R:234-239 | outlier detector | **none (dead code)** | IQR-fence outlier flag | numeric vector, k | logical vector |
| `methyl_sample_outliers_correlation` | qc.R:255-263 | outlier detector | Outlier QC | correlation-based outlier flag | matrix, n_features, k | ok/outlier list |
| `methyl_sample_outliers_pca` | qc.R:274-286 | outlier detector | Outlier QC | PCA-distance outlier flag | matrix, n_features, sd_threshold | ok/scores/outlier list |
| `methyl_sample_outliers_hclust` | qc.R:291-304 | outlier detector | Outlier QC | clustering-singleton outlier flag | matrix, n_features, height_frac | ok/hc/outlier list |
| `methyl_sample_outliers_mahalanobis` | qc.R:319-338 | outlier detector | Outlier QC | Mahalanobis-distance outlier flag | matrix, n_features, n_pcs, alpha | ok/distance2/outlier list |
| `methyl_cluster_sex` | qc.R:356-376 | sex prediction | `methyl_sex_check` | k-means/median-split sex clustering | chrY beta, reported sex | sex/direction list |
| `methyl_sex_check_attach_mismatch` | qc.R:383-392 | sex QC helper | `methyl_sex_check` | attach concordance flag | detail df, reported sex | detail list |
| `methyl_sex_check` | qc.R:407-446 | sex QC dispatcher | Sex QC | predict + compare sex | matrix, annotation, rg_set, reported | ok/method/detail list |
| `methyl_sheet_sample_ids` | qc.R:463-468 | ID resolution | subgroup/sex/batch/viz | resolve sheet rows to matrix IDs | sheet, all_ids | character vector |
| `methyl_qc_subgroup_filter` | qc.R:479-500 | data scoping | `current_subgroup` | stratum filtering | matrix, sheet, col, level | subgroup list |
| `methyl_apply_manual_exclude` | qc.R:507-519 | data scoping | `current_subgroup` | manual exclusion | subgroup, excluded IDs | subgroup list |
| `methyl_outlier_score_table` | qc.R:529-539 | aggregation | Outlier QC | combine outlier flags into a score | sample_qc df | ranked df |
| `methyl_probe_retention_cascade` | qc.R:547-560 | aggregation | Probe QC, reports | sequential retention counts | n_start, filters | df |
| `methyl_beta_to_mvalue` | qc.R:566-569 | transform | export, batch correction | beta→M-value logit | beta matrix, eps | M-value matrix |
| `methyl_pca_scores` | qc.R:576-587 | structure analysis | Batch QC, Visualizations | top-variance PCA | matrix, n_features, n_pcs | ok/scores/var list |
| `methyl_sample_correlation` | qc.R:593-601 | structure analysis | Outlier QC, Visualizations | sample correlation matrix | matrix, n_features | ok/cor list |
| `methyl_beta_density_sample` | qc.R:607-616 | data reshape | Visualizations, reports | long-format probe subsample | matrix, n_probes, seed_probes | df |
| `methyl_mds_scores` | qc.R:624-637 | structure analysis | Visualizations | classical MDS | matrix, n_features, k | ok/scores list |
| `methyl_mean_sd_table` | qc.R:643-650 | diagnostic | Visualizations, reports | per-probe mean vs SD | matrix, n_probes | df |
| `methyl_control_probe_matrix` | qc.R:656-675 | IDAT diagnostic | Visualizations | control-probe intensities | rg_set | ok/mat/types list |
| `methyl_outlier_diagnostic_table` | qc.R:681-696 | aggregation | Outlier QC, reports | combine distance metrics | sample_qc, pca_detail, mahal_detail | df |
| `methyl_guess_batch_column` | qc.R:704-709 | metadata heuristic | **none (dead code — logic duplicated inline)** | guess batch column | sheet | character or NULL |
| `methyl_batch_correct_combat` | qc.R:723-750 | batch correction | Batch QC | ComBat on M-values | matrix, batch | ok/corrected list |
| `methyl_batch_correct_ruvm` | qc.R:764-808 | batch correction | Batch QC | RUVm via control probes | matrix, rg_set, group, k | ok/corrected list |
| `methyl_plot_cascade` | qc.R:818-826 | plot builder | Probe QC, reports | retention-cascade bar chart | cascade df | ggplot |
| `methyl_plot_detp_heatmap` | qc.R:828-834 | plot builder | Visualizations | detection-p heatmap | long df | ggplot |
| `methyl_plot_beadcount_dist` | qc.R:836-841 | plot builder | Visualizations | bead-count histogram | values, threshold | ggplot |
| `methyl_plot_density` | qc.R:851-856 | plot builder | Visualizations, reports | beta density | density df | ggplot |
| `methyl_plot_boxplot` | qc.R:863-869 | plot builder | Visualizations, reports | per-sample boxplot | long df | ggplot |
| `methyl_plot_violin` | qc.R:871-877 | plot builder | Visualizations, reports | per-sample violin | long df | ggplot |
| `methyl_plot_mean_sd` | qc.R:879-884 | plot builder | Visualizations, reports | mean-SD scatter | mean_sd df | ggplot |
| `methyl_plot_scatter2d` | qc.R:890-896 | plot builder | multiple | generic 2D scatter | x/y/color/text df | ggplot |
| `methyl_plot_corr_heatmap` | qc.R:898-910 | plot builder | Visualizations, reports | correlation heatmap | cor matrix | ggplot |
| `methyl_plot_outlier_diagnostic` | qc.R:912-919 | plot builder | Outlier QC, reports | ranked distance bar chart | diag df, metric | ggplot |
| `methyl_qc_status_badge` | qc.R:929-944 | decision function | Overview | pass/warning/fail verdict | overview result | status/reasons list |
| `methyl_qc_r_code` | qc.R:951-973 | codegen | Reports & Export | Bioconductor-equivalent snippet | probe QC settings | text |
| `methyl_qc_summary_table` | qc.R:985-1034 | aggregation | Reports & Export | cross-tab metric summary | 6 result objects | df |
| `methyl_qc_report_plots` | qc.R:1045-1089 | aggregation | Reports & Export | rebuild export figures | dataset, probe_qc, outlier_qc | plots/skipped list |
| `methyl_qc_report_html` | qc.R:1099-1134 | export builder | Reports & Export | self-contained HTML report | dataset, summary, plots | ok/path list |
| `methyl_qc_report_zip` | qc.R:1141-1157 | export builder | Reports & Export | figures ZIP | plots | ok/path list |
| `methyl_idat_derive` | idat_metrics.R:13-28 | IDAT processing | Dataset tab (not QC directly) | derive beta/detp/beadcount | rg_set | ok/beta/detp/beadcount list |
| `methyl_bisulfite_conversion` | idat_metrics.R:33-40 | sample metric | Sample QC | bisulfite conversion % | rg_set | ok/pct list |
| `methyl_median_intensity` | idat_metrics.R:45-53 | sample metric | Sample QC | median intensity | mset | ok/detail list |
| `methyl_get_annotation` | annotation.R:48-93 | annotation loader | Probe QC, Sex QC | manifest annotation | array_type | ok/anno list |
| `methyl_probe_is_cpg` | annotation.R:99-101 | probe classifier | `methyl_filter_non_cpg` | ID-prefix CpG check | probe IDs | logical vector |
| `methyl_parse_matrix` | parse_upload.R:11-37 | upload parser | Dataset tab (not QC directly) | parse beta/M matrix | file path | ok/mat/error |
| `methyl_parse_sample_sheet` | parse_upload.R:40-46 | upload parser | Dataset tab (not QC directly) | parse sample sheet | file path | ok/df/error |
| `methyl_parse_probe_list` | parse_upload.R:50-57 | upload parser | Probe QC | parse probe-exclusion list | file path | ok/ids/error |
| `methyl_read_idat` | parse_upload.R:65-96 | upload parser | Dataset tab (not QC directly) | read raw IDAT | fileInput df | ok/rg/error |

**Total functions inventoried: 62** (58 called by QC directly or indirectly; 2 confirmed dead code — `methyl_sample_outliers_iqr`, `methyl_guess_batch_column`; 2 called only from the Dataset tab, not QC itself, but documented since QC consumes their output — `methyl_idat_derive`, and the three `parse_upload.R` matrix/sheet/IDAT parsers are grouped as related-but-external).

---

## 17. Tab inventory table

| Tab | Purpose | Input | Main functions | Processing | Output |
|---|---|---|---|---|---|
| Overview | Dataset summary + basic pass/fail QC | `methyl_dataset$beta`, `sample_sheet`, `current_subgroup()` | `methyl_sample_call_rate`, `methyl_qc_status_badge` | call rate + missingness → hard-coded threshold verdict | status badge, 3 valueBoxes, reasons list |
| Sample QC | Per-sample call rate, intensity, conversion diagnostics; manual exclusion UI | `current_subgroup()$mat`, `rg_set`/`mset`/`detp` | `methyl_sample_call_rate`, `methyl_sample_failed_probe_pct`, `methyl_sample_low_intensity`, `methyl_bisulfite_conversion`, `methyl_median_intensity` | flag computation only, no removal | 3 DT tables, manual-exclusion table |
| Probe QC | 11-filter probe-level filtering | `current_subgroup()$mat`, annotation, `detp`/`beadcount`, uploads | 11 `methyl_filter_*` functions, `methyl_probe_retention_cascade` | AND-combine enabled filters; separate cumulative cascade | filter table, cascade plot, filtered matrix + downloads |
| Sex QC | Predicted-vs-reported sex concordance | `current_subgroup()$mat`, annotation, `rg_set`, sheet sex column | `methyl_sex_check`, `methyl_cluster_sex` | `minfi::getSex()` or chrY k-means clustering | DT table, X-vs-Y scatter, discordant table + exclude button |
| Batch QC | ComBat/RUVm batch correction visualization | `current_subgroup()$mat`, sheet batch column, `rg_set` | `methyl_batch_correct_combat`, `methyl_batch_correct_ruvm`, `methyl_pca_scores` | logit to M-value, correct, PCA before/after | before/after PCA plots, variance table (no export) |
| Outlier QC | 4-method sample-outlier detection | `current_subgroup()$mat` | `methyl_sample_outliers_pca/hclust/correlation/mahalanobis`, `methyl_outlier_score_table` | all 4 detectors always computed; displayed subset selectable | ranked score table, PCA/diagnostic/dendrogram plots, exclude button |
| Visualizations | Re-plots Probe-QC-filtered data + raw IDAT diagnostics | `probe_qc_result()`, `methyl_dataset$detp/beadcount`, `current_rg_subset()` | `methyl_pca_scores`, `methyl_mds_scores`, `methyl_sample_correlation`, `methyl_mean_sd_table`, `methyl_control_probe_matrix`, 8 plot builders | lazy (button-per-plot) rendering only | 11 distinct plots |
| Reports & Export | Aggregate export of whichever tabs have run | all 6 other tabs' result objects via `current_qc_pieces()` | `methyl_qc_summary_table`, `methyl_qc_report_plots`, `methyl_qc_report_html`, `methyl_qc_report_zip`, `methyl_qc_r_code` | assembles already-computed pieces into files | 10 distinct downloads/panels (§3.8 table) |

---

## 18. Data-flow summary diagram

This reflects the actual architecture verified in §4/§11 — a shared, stratum-filtered dataset with independent per-tab analyses, not a strict pipeline:

```
┌────────────────────────────────────────────────────────────────────────┐
│  methyl_dataset (reactiveValues, populated ONLY by the Dataset tab)    │
│  beta · sample_sheet · rg_set · mset · detp · beadcount ·              │
│  array_type · input_scale · preloaded · source                        │
└───────────────────────────────┬──────────────────────────────────────┘
                                 │  read-only by every tab below
                                 ▼
                  ┌───────────────────────────┐
                  │   current_subgroup()      │◄── manual_exclude() (shared,
                  │  (stratum + exclusions)   │    written by 3 tabs' buttons)
                  └──────────────┬────────────┘
     ┌──────────┬────────────────┼────────────────┬──────────┬───────────┐
     ▼          ▼                ▼                 ▼          ▼           ▼
 Overview   Sample QC        Probe QC           Sex QC    Batch QC   Outlier QC
 (button)   (button)         (button)          (button)   (button)   (button)
     │          │                │ .filtered       │          │           │
     │          │                ▼                 │          │           │
     │          │         Visualizations            │          │           │
     │          │        (button-per-plot)          │          │           │
     └──────────┴────────────────┴──────────────────┴──────────┴───────────┘
                                 │
                                 ▼  current_qc_pieces() — read-only aggregation
                       Reports & Export
                    (CSV / HTML / PDF / ZIP downloads)
                                 │
                                 ▼
                 outside the app only (manual re-upload/re-use);
              methyl_dataset is NEVER mutated; methyl_results is NEVER written
```

---

## 19. Closing paragraph (XomicsShiny-referenced style)

> The Methylomics Quality Control sub-module (`mod_methyl_qc.R`/`qc.R`, registered as `"qc"` in the "Data" group) provides **eight independent, button-driven live tabs** — Overview, Sample QC, Probe QC, Sex QC, Batch QC, Outlier QC, Visualizations, and Reports & Export — plus a separate, collapsed-by-default historical reference section that reproduces (never recomputes) the completed offline QC run behind the app's preloaded whole-blood dataset. Its input is the shared `methyl_dataset` object populated by the Dataset tab: a probe-by-sample beta- or M-value matrix, optional sample metadata, and, for a raw-IDAT upload only, the underlying intensity data (`RGChannelSet`, `MethylSet`, detection p-values, bead counts). Each tab independently reads a common, always-live stratum-and-exclusion filter (`current_subgroup()`) rather than each other's results — the one exception being that Visualizations re-plots Probe QC's filtered matrix once Probe QC has been run. Major operations include eleven composable probe-level filters (detection p-value, bead count, missingness, SNP overlap, non-CpG design, sex-chromosome location, variance/SD/mean-range, and optional user-uploaded cross-reactive-probe and MAF exclusion lists), four sample-outlier detectors (PCA-distance, hierarchical-clustering, correlation-based, and classical Mahalanobis distance), a raw-intensity or beta-clustering sex check, and two batch-correction methods (ComBat and RUVm) built on `minfi`, `sva`, `missMethyl`, and `wateRmelon`. Quality control matters here because every filter and detector is explicit about what it *would* remove or flag without ever altering the shared dataset itself — samples are excluded only through a single, explicit, user-triggered mechanism shared across three tabs, and probes are excluded only in the exported, downloadable matrix, never retroactively. Outputs are per-tab diagnostic tables and plots, a probe-retention cascade, downloadable filtered beta/M-value matrices, a cross-tab QC summary CSV, and a self-contained HTML report (PDF optional, environment-dependent) with every figure baked in. Together, the eight tabs give a user a complete, independently re-runnable methylation QC audit — sample-level, probe-level, sex-concordance, batch-structure, and outlier — organized so that no tab's output silently depends on another's having been run, while still sharing one consistent, exclusion-aware view of "which samples are currently in scope."

---

*End of document. Prepared by direct code inspection of `mod_methyl_qc.R` (1,619 lines, read in full), `qc.R` (1,157 lines, read in full), `idat_metrics.R`, `annotation.R`, `parse_upload.R`, `mod_methyl_dataset.R`, relevant sections of `normalization.R`, `submodules_registry.R`, `server.R`, and `global.R`. No application code was modified in the course of this audit.*
