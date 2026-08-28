# Methylomics Normalization Module: `mod_methyl_normalization.R`

**Source files:**
- `ArthOMix/R/methylomics/mod_methyl_normalization.R` (855 lines) — UI + server for the "Normalization" sub-module.
- `ArthOMix/R/methylomics/normalization.R` (551 lines) — normalization-method implementations, diagnostics, status-detection, validation, and Normalization-tab-only annotation/filter helpers called by `mod_methyl_normalization.R`.
- Supporting files read for this document: `ArthOMix/R/methylomics/qc.R` (probe/sample filters and plot builders shared with Quality Control), `ArthOMix/R/methylomics/annotation.R` (102 lines), `ArthOMix/R/methylomics/parse_upload.R` (97 lines), `ArthOMix/R/methylomics/mod_methyl_dataset.R` (298 lines), `ArthOMix/R/methylomics/idat_metrics.R` (raw-IDAT beta derivation only), `ArthOMix/R/submodules_registry.R`, `ArthOMix/global.R` (preloaded-data-loading and package-attach sections), and — for cross-checking the QC↔Normalization boundary — `ArthOMix/R/methylomics/mod_methyl_qc.R`.

**Registration:** `mod_methyl_normalization_config` — `id = "normalization"`, `title = "Normalization"`, `icon = "wave-square"`, `group = "Data"` (`mod_methyl_normalization.R:36-39`). Registered second in `MX_MODULES`, immediately after Quality Control (`submodules_registry.R:40-41`), invoked as `mod_methyl_normalization_server("mx_normalization", methyl_dataset, methyl_results)` via the same `lapply(MX_MODULES, ...)` loop QC is invoked from (`server.R:95`).

Prepared: 2026-08-26.

This document is derived **exclusively** from the code cited above and was scoped strictly to Methylomics → Normalization; no other module was modified or redesigned in the course of this audit. Every non-trivial technical claim carries a `file:line` citation. Where a claim could not be verified from the inspected code, that is stated explicitly as **Not determinable from the inspected implementation** rather than inferred. Two label conventions, matching the companion `methylomics_quality_control.md` document's own convention, are used throughout:

- **Scientific background:** a statement about methylation-array normalization in general (textbook/literature knowledge), not a claim about this code.
- **Code evidence:** a statement about what `mod_methyl_normalization.R`/`normalization.R` actually do, always with a citation.

---

## 1. Module overview

**Scientific background.** Illumina methylation-array measurements (450K/EPIC) carry several well-characterized sources of technical variation on top of true biological methylation signal: background fluorescence and dye-bias differences between the Cy3/Cy5 channels, and — specific to the Infinium chemistry — a systematic distributional difference between Type I and Type II probe designs, because the two chemistries measure methylation with different dynamic range and noise characteristics at the same true methylation level. Normalization methods correct one or both of these effects. This is scientifically distinct from *quality control* (removing individual probes/samples that fail a quality threshold) and from *batch correction* (removing site/plate/scan-date effects, typically requiring a known batch label) — a dataset can be well-normalized and still contain batch effects, or vice versa, since these are three independent technical-variation sources with different signatures and different corrective methods.

**Code evidence — the module's own header states this exact distinction.** `mod_methyl_normalization.R`'s header comment explains the two data pathways it serves and states the deliberate choice behind each (`mod_methyl_normalization.R:1-34`), and `normalization.R`'s header separately draws the line between "background/dye-bias correction" (Noob) and "probe-design/distribution normalization" (BMIQ/SWAN/PBC/Dasen/stratified quantile), stating explicitly that the two "sequential (two-step)" combination methods (Noob+BMIQ, Noob+SWAN) exist specifically because "these two kinds of correction" are not interchangeable (`normalization.R:185-190`). Batch correction (ComBat/RUVm) is **not** implemented in this module at all — it lives entirely in the Quality Control module's own Batch QC tab (`qc.R:723-808`, documented in `methylomics_quality_control.md` §3.5); Normalization contains no batch-label input and no ComBat/RUVm call anywhere in `mod_methyl_normalization.R` or `normalization.R` (verified by exhaustive grep for `ComBat`/`RUVm`/`batch` in both files — zero hits).

**Two data pathways, one shared diagnostics engine (code-confirmed, from the module's own header).** Both the preloaded whole-blood dataset and an uploaded dataset run through the identical `methyl_norm_diagnostics()` and `methyl_norm_status()` functions (`mod_methyl_normalization.R:7-10,105-112`). What differs is what happens next:
- The **preloaded** dataset was analyzed from its original author-normalized data (Liu et al. 2013, GEO-deposited) rather than reprocessed in this app — "a deliberate choice that preserves direct comparability with the originally published beta values rather than risking systematic differences from an independently re-implemented normalization" (`mod_methyl_normalization.R:11-17`). This tab therefore never offers a live method picker, filters, or a Run button for that dataset.
- An **uploaded** dataset gets the full live workflow: filters, a method picker whose choices depend on what was actually uploaded, a Compare Methods panel, and progressively-revealed before/after results (`mod_methyl_normalization.R:18-27`).

**Input data (code-confirmed).** The module reads six fields off the shared `methyl_dataset` reactiveValues object (initialized `server.R:82-92`, populated only by `mod_methyl_dataset.R`):
- `methyl_dataset$beta` — the probe-by-sample matrix, rows = probes/CpGs, columns = samples (verified §9 below), on either the beta (0–1) or M-value scale per `methyl_dataset$input_scale` (`"beta"` / `"m"`, set at `mod_methyl_dataset.R:92,208,278`).
- `methyl_dataset$array_type` — `"450K"`/`"EPIC"`/`"EPICv2"`/`"WGBS"`/`"RRBS"`/`"Custom array"` (`annotation.R:10`), used to decide which methods/filters are offered.
- `methyl_dataset$rg_set` — a `minfi::RGChannelSet`, present only for an IDAT upload (`mod_methyl_dataset.R:281`); always `NULL` for the preloaded dataset and for a matrix upload (`mod_methyl_dataset.R:95,211`).
- `methyl_dataset$mset` — a `minfi::MethylSet` derived from raw IDAT (`minfi::preprocessRaw()`, `idat_metrics.R:17`), needed by SWAN/Dasen.
- `methyl_dataset$detp` — a raw-IDAT detection-p-value matrix (`minfi::detectionP()`, `idat_metrics.R:24`), used only by the Filters tab's optional sample-detection-rate filter.
- `methyl_dataset$sample_sheet`, `methyl_dataset$preloaded`, `methyl_dataset$source` — provenance/metadata, used to gate the biological-signal-preservation check and to label the promoted dataset.

No file is read directly by this module — every input arrives pre-parsed from the Dataset tab (`mod_methyl_dataset.R`, out of scope for this audit); Normalization itself performs no file I/O beyond an optional cross-reactive-probe-list upload on its own Filters tab (`mod_methyl_normalization.R:359`, `methyl_parse_probe_list()`, `parse_upload.R:50-57`).

**Example structure (matches the actual implementation, verified §9):**
```
Rows    → methylation probes/CpGs (Illumina cg-prefixed IDs, or genomic-coordinate IDs for WGBS/RRBS)
Columns → samples
Cells   → beta values (0-1 methylation proportion) or M-values (logit-transformed), per methyl_dataset$input_scale
```

**What is produced.** A candidate normalized matrix, before/after diagnostic plots and statistics, a filtered-probe/removed-sample list, a plain-text processing record, and — only on explicit user action — promotion of that candidate to the shared `methyl_dataset$beta`.

**How this connects to the rest of the app.** The module writes to two shared objects: `methyl_results$normalization` (a small summary list, written unconditionally once a run completes, `mod_methyl_normalization.R:575-579`) and, only on the explicit "Use this as the active Methylomics dataset" button, `methyl_dataset$beta`/`$input_scale`/`$source` (`mod_methyl_normalization.R:581-594`). **Code evidence:** `methyl_results$normalization` is written but never read anywhere else in the codebase (`grep -rn "methyl_results\$normalization" R/` returns only the one write site, `mod_methyl_normalization.R:577`) — see §14, Finding N-1.

---

## 2. Tab count

**Number of Normalization tabs: this depends on which data pathway is active — the tab count itself is not fixed.**

- **Preloaded dataset:** **zero** sub-tabs. `preloaded_ui()` renders two stacked cards (diagnostics, status) directly, with no `tabsetPanel` at all (`mod_methyl_normalization.R:184-218`). "'Re-normalize' and 'Compare normalization methods' are intentionally not offered for this dataset" (`mod_methyl_normalization.R:211`).
- **Uploaded dataset, and `status$status != "no_bias_detected"` (or the user has clicked "Re-normalize"/"Compare normalization methods" after an "already normalized" verdict):** **3 live sub-tabs**, in one `tabsetPanel(id = ns("config_tabs"), type = "tabs")` (`mod_methyl_normalization.R:302-312`), in this exact order and with these exact displayed labels:

| # | Tab (exact label from code) | `tabPanel` value |
|---|---|---|
| 1 | **Filters** | `"Filters"` |
| 2 | **Method & Run** | `"Method & Run"` |
| 3 | **Compare Methods** | `"Compare Methods"` |

- **Uploaded dataset, `status$status == "no_bias_detected"`, and the user has not yet chosen an action:** **zero** sub-tabs — only the diagnostics card, the status card, and three choice buttons ("Keep current normalization" / "Re-normalize" / "Compare normalization methods") are shown; the 3-tab `tabsetPanel` above is conditionally omitted entirely (`show_live_workflow()`, `mod_methyl_normalization.R:258-262,288-315`).

A `tabsetPanel` renders every `tabPanel`'s body up front (tab switching is CSS visibility, not conditional server-side rendering) — the module's own comment states this explicitly and explains why it matters: "filters set on the Filters tab stay live and are read by both Method & Run and Compare Methods regardless of which tab is currently open" (`mod_methyl_normalization.R:293-301`). This was verified: `filters_ui()`'s inputs (`input$f_probe_missing`, etc., `mod_methyl_normalization.R:319-385`) are read by `build_probe_filters()` (`mod_methyl_normalization.R:396-406`), which is called from `run_full()` (`mod_methyl_normalization.R:522-558`) — the single shared code path both "Run Normalization" (Method & Run tab) and "Compare Methods" (Compare Methods tab) invoke, so the same live filter state feeds both regardless of which tab is currently visible.

No nested (grandchild) tabs exist anywhere in this module — confirmed by grepping `mod_methyl_normalization.R` for `tabsetPanel`: exactly one call site (`mod_methyl_normalization.R:303`).

---

## 3. Per-pathway / per-tab documentation

### 3.0 Preloaded-dataset branch (not a tab — `preloaded_ui()`)

**Purpose.** Show the preloaded whole-blood cohort's automatic diagnostics and a documented, non-recomputed "already normalized (by its original authors)" status — never offer live reprocessing.

**Input data.** `methyl_dataset$beta` (the preloaded, author-normalized beta matrix, gated by `isTRUE(methyl_dataset$preloaded)` at `mod_methyl_normalization.R:171`); if the live matrix isn't bundled in a given deployment (`methyl_dataset$beta` is `NULL` even though `preloaded == TRUE`), a plain "No separate normalization step for this dataset" card is shown instead, citing the Dataset tab as the source of truth (`mod_methyl_normalization.R:185-191`).

**Functions called.** `methyl_norm_diagnostics(mat, dataset, anno_result)` (`normalization.R:343-367`) → `diagnostics_card()` (`mod_methyl_normalization.R:114-128`); `methyl_norm_status(mat, dataset, anno_result)` (`normalization.R:424-453`) → read directly (`st$bias`), **not** through the shared `status_card()` helper.

**Code evidence — the status badge is deliberately overridden for this dataset.** The preloaded branch is the one place in the module that does *not* call `status_card()`; its own comment explains why: `methyl_norm_status()`'s badge reflects only the narrow Type I/II probe-design-bias signal, which is "a DIFFERENT question from 'was this dataset normalized'" — for this specific cohort the ground-truth answer to the latter is already documented (author-normalized before GEO deposition), so "a heuristic should [not] second-guess or appear to contradict" it (`mod_methyl_normalization.R:193-204`). Instead a fixed green "Using existing (author) normalization" badge is shown (`mod_methyl_normalization.R:207-211`), and the heuristic's own KS-statistic reading is demoted to a secondary "Technical note" (`mod_methyl_normalization.R:212-215`).

**Statistical operations.** Identical to §5's `methyl_norm_diagnostics()`/`methyl_type_bias_stat()` — no different computation for this pathway, only different framing of the result.

**Plots/tables.** None — value boxes and prose only (`mod_methyl_normalization.R:206-217`).

**Output objects.** None persisted; this branch performs no run, no promotion, and offers no downloads.

**Connection to other tabs.** None — this branch is a dead end by design; the module's header states the preloaded cohort's dataset-tab-documented provenance is the authority, not this tab (`mod_methyl_normalization.R:11-17`).

**Audit findings for this branch.** None — the deliberate suppression of a live re-normalization option here is well-reasoned and disclosed (see §19 for the corresponding QC-relationship note). One minor observation: the "Technical note" text states GEO-deposited beta values from this era "typically received background/dye/batch correction but not necessarily a probe-design-aware step" (`mod_methyl_normalization.R:213`) — this is presented as general historical context about the field, correctly distinguished from a claim about this specific cohort's own processing history, which the text does not otherwise assert.

---

### 3.1 Tab 1: Filters

**Purpose.** Collect optional, non-destructive probe- and sample-level filter configuration, applied only when "Run Normalization" or "Compare Methods" is clicked — "nothing here changes the loaded dataset by itself" (`mod_methyl_normalization.R:331`).

**Input data.** `methyl_dataset$sample_sheet` (optional, populates the biological-signal-preservation group-column selector); `norm_anno_result()` (island/gene-region columns, `normalization.R:252-282`); `anno_result()` (chromosome column, `annotation.R:48-93`).

**User inputs, by section (`mod_methyl_normalization.R:329-384`):**

*General:*
| Input ID | Type | Default | Purpose |
|---|---|---|---|
| `f_probe_missing` / `probe_missing_max` | checkboxInput / numericInput | off / 0.05 | max missing fraction per probe |
| `f_sample_missing` / `sample_missing_max` | checkboxInput / numericInput | off / 0.1 | max missing fraction per sample |
| `f_sample_detrate` / `detp_thresh`, `min_det_rate` | checkboxInput / numericInput ×2 | off / 0.01, 95 | detection-p threshold and minimum sample detection rate; IDAT-only, disabled with an explanatory note otherwise (`mod_methyl_normalization.R:349`) |

*Methylation-specific (Illumina array types only — hidden entirely otherwise, `mod_methyl_normalization.R:352-378`):*
| Input ID | Type | Default | Purpose |
|---|---|---|---|
| `f_snp` | checkboxInput | off | remove SNP-overlapping probes |
| `f_crossreactive` / `crossreactive_file` | checkboxInput / fileInput | off / none | remove probes in a user-uploaded exclusion list |
| `f_sexchr` | checkboxInput | off | remove chrX/chrY probes |
| `f_chr` / `exclude_chr` | checkboxInput / selectizeInput | off / none | exclude specific chromosome(s) |
| `f_island` / `island_categories` | checkboxInput / checkboxGroupInput | off / all categories pre-selected | keep only selected CpG-island-relation categories (hidden if annotation unavailable) |
| `f_generegion` / `gene_regions` | checkboxInput / checkboxGroupInput | off / all regions pre-selected | keep only selected gene-region categories (hidden if annotation unavailable) |

*Biological-signal-preservation check (shown only if a sample sheet with ≥1 column is loaded):* `group_col_check` (`selectInput`, default `""` = "(none)", `mod_methyl_normalization.R:381-383`) — see §7 Finding N-3 for why this input's default matters.

**Reactive dependencies.** No `eventReactive`/button on this tab itself — every input here is read live, at click time, by `build_probe_filters()` (probe-level, `mod_methyl_normalization.R:396-406`) and `sample_scope()` (sample-level, `mod_methyl_normalization.R:489-509`), both called from `run_full()` (`mod_methyl_normalization.R:522-558`).

**Functions called (with citations).** `methyl_filter_missing()` (`qc.R:30-34`), `methyl_filter_snp()` (`qc.R:68-81`), `methyl_filter_sex_chr(mode="remove_xy")` (`qc.R:89-104`), `methyl_filter_chromosome()` (`normalization.R:314-323`, Normalization-tab-only), `methyl_filter_cross_reactive()` (`qc.R:112-119`), `methyl_filter_island_relation()` / `methyl_filter_gene_region()` (`normalization.R:290-300`, `302-312`, both Normalization-tab-only), `methyl_filter_samples_missingness()` (`normalization.R:331-336`, Normalization-tab-only), `methyl_sample_failed_probe_pct()` (`qc.R:202-217`), `methyl_parse_probe_list()` (`parse_upload.R:50-57`).

**Statistical operations.** Identical primitives to Quality Control's own probe filters where shared (`methyl_filter_missing`/`methyl_filter_snp`/`methyl_filter_sex_chr`/`methyl_filter_cross_reactive` are the literal same functions, imported from `qc.R`, not re-implemented) — see §19 for the significance of this sharing. Three filters are unique to this tab: chromosome exclusion by name (`methyl_filter_chromosome()`), CpG-island-relation and gene-region filters (using a second, wider annotation object `methyl_get_norm_annotation()` that adds `Relation_to_Island`/`UCSC_RefGene_Group` columns Quality Control's own `methyl_get_annotation()` does not carry, `normalization.R:239-249`), and sample-level missingness (`methyl_filter_samples_missingness()`) — this one, unlike everything on Quality Control's own Sample QC tab, **genuinely drops columns before normalization runs** rather than only flagging them, per its own comment: "a sample too sparse to normalize meaningfully shouldn't silently ride along" (`normalization.R:325-330`).

**Plots/tables.** None on this tab — it is configuration-only.

**Output objects.** None directly; every input is read by `build_probe_filters()`/`sample_scope()` at Run-time.

**Interpretation.** Not applicable — no result to interpret on this tab.

**Connection to other tabs.** Feeds both Method & Run and Compare Methods identically (§2).

**Audit findings for this tab.** See §7 for the method-ordering rationale (probe filters applied before vs. after the normalization step, by design) and §14 Finding L-6 (sex-chromosome/SNP/chromosome-exclusion checkboxes remain visibly checkable for EPICv2/Custom-array datasets even though no manifest annotation exists for those array types, silently no-opping rather than disabling).

---

### 3.2 Tab 2: Method & Run

**Purpose.** Pick one normalization method (from only the methods actually compatible with the loaded data), configure its parameters, and run it — producing a promotable candidate result plus progressively-revealed before/after diagnostics.

**Input data.** `methyl_dataset$beta` (all method types, for the "before" comparison and for matrix-based methods); `methyl_dataset$rg_set` (raw-intensity methods); `methyl_dataset$mset` (SWAN/Dasen); `anno_result()` (BMIQ/PBC/Noob+BMIQ Type I/II design vector).

**User inputs — method picker and per-method parameters (`method_ui()`, `mod_methyl_normalization.R:412-444`):**
| Input ID | Type | Applies to | Default |
|---|---|---|---|
| `method` | radioButtons | all | `default_method()` (see below) |
| `noob_offset` | numericInput | Noob, Noob+SWAN, Noob+BMIQ | 15 |
| `noob_dye` | radioButtons | Noob, Noob+SWAN, Noob+BMIQ | `"single"` (ssNoob) |
| `funnorm_npcs` | numericInput | Functional normalization | 2 |
| `bmiq_nfit` | numericInput | BMIQ, Noob+BMIQ | 50000 |
| `bmiq_nl` | numericInput | BMIQ, Noob+BMIQ | 3 |
| `bmiq_tol` | numericInput | BMIQ, Noob+BMIQ | 0.001 |
| `run_btn` | actionButton | all | — |

SWAN, Dasen, Stratified quantile, PBC, plain quantile, and "no normalization" expose no tunable parameters beyond the method choice itself, and the UI states this explicitly (`mod_methyl_normalization.R:439-440`).

**Method availability logic — `available_methods()` (`mod_methyl_normalization.R:222-230`), a pure function of what data is actually loaded:**
```
has_idat()               → Noob, Functional normalization, SWAN, Dasen, Stratified quantile, Noob+SWAN
has_idat() & manifest     → + Noob+BMIQ
manifest & is_beta_scale  → + BMIQ, PBC
always                    → No normalization / keep current data, Quantile normalization (plain)
```
Every branch's un-met precondition is shown to the user as an explanatory note (`mod_methyl_normalization.R:417-422`) rather than silently omitting the option with no explanation.

**Default-method logic — `default_method()` (`mod_methyl_normalization.R:232-238`):** if the automatic status check reads "no bias detected" and "none" is available, default to "none"; else prefer "noob" if available, else "bmiq" if available, else the last-listed available method (which — given `available_methods()`'s construction order — is always at least "quantile").

**Reactive dependencies.** `norm_result <- eventReactive(input$run_btn, ...)` (`mod_methyl_normalization.R:565-573`) — `req(input$method)`; `validate(need(!is.null(methyl_dataset$beta), ...))`; wraps `run_full(input$method)` in `withProgress()`. `norm_has_run` (`reactiveVal`, `mod_methyl_normalization.R:562-563`) resets to `FALSE` whenever `methyl_dataset$beta` changes (new dataset loaded) and flips `TRUE` in an `observeEvent(norm_result(), ...)` (`mod_methyl_normalization.R:575-579`) that also writes `methyl_results$normalization`.

**Functions called (with citations) — dispatched by `run_one_method()` (`mod_methyl_normalization.R:470-486`):** `methyl_norm_noob()` (`normalization.R:60-67`), `methyl_norm_funnorm()` (`normalization.R:69-76`), `methyl_norm_swan()` (`normalization.R:78-85`), `methyl_norm_dasen()` (`normalization.R:96-109`), `methyl_norm_stratified_quantile()` (`normalization.R:87-94`), `methyl_norm_noob_swan()` (`normalization.R:199-209`), `methyl_norm_noob_bmiq()` (`normalization.R:191-197`), `methyl_norm_bmiq()` (`normalization.R:121-154`), `methyl_norm_pbc()` (`normalization.R:158-172`), `methyl_norm_quantile()` (`normalization.R:176-183`). The orchestration function `run_full()` (`mod_methyl_normalization.R:522-558`) applies `sample_scope()` first, then dispatches to `run_one_method()`, then `build_probe_filters()` (order depends on method type — see §7/§9), then `methyl_norm_validation()` (`normalization.R:479-513`).

**Statistical operations.** Documented per-method in §7.

**Plots/tables (progressively revealed, each behind its own `.methyl_norm_toggle_ui()` "Show ..." link, `mod_methyl_normalization.R:69-78,658-663`):**
- **Show before vs. after:** `density_plot` (`methyl_plot_density()`, `qc.R:863-868`), `boxplot_plot` (`methyl_plot_boxplot()`, `qc.R:875-881`) — both auto-rendered; plus two button-gated additions: `bvsa_pca_plot` (before/after PCA via `methyl_pca_scores()`, `qc.R:569-580`, faceted by stage) and `bvsa_corr_before_plot`/`bvsa_corr_after_plot` (side-by-side correlation heatmaps via `methyl_sample_correlation()`+`methyl_plot_corr_heatmap()`, `qc.R:586-594,910-919`).
- **Show normalization statistics:** `stats_body` — a DT table of mean row variance, mean sample-sample correlation, Type I/II KS statistic, and missingness, before vs. after (`mod_methyl_normalization.R:741-756`), plus a PC1~group R² row if the biological-signal check ran.
- **Show filtered probes / removed samples:** two DT tables (`mod_methyl_normalization.R:758-775`).
- **Show processing details:** a plain-text `methyl_norm_processing_record()` dump (`normalization.R:538-550`, rendered `mod_methyl_normalization.R:777-788`).

**Output objects.** `norm_result()` — list with `ok`, `method`, `method_label`, `before`, `after`, `note`, `removed_probes`, `removed_samples`, `filter_notes`, `n_probes_before`, `n_samples_before`, `validation`, `run_at` (`mod_methyl_normalization.R:553-557`).

**Downloads (all gated on `norm_result()` existing, `mod_methyl_normalization.R:596-622`):** normalized matrix (CSV), filtered probe list (CSV), removed sample list (CSV), processing record (TXT). None require a separate confirmation step beyond having run the method once.

**Interpretation.** `methyl_norm_interpretation()` (`normalization.R:520-535`) — see §7's dedicated discussion; never declares success merely because the algorithm completed.

**Connection to other tabs.** Reads Filters-tab inputs live (§2); "Use this as the active Methylomics dataset" (`promote_btn`, `mod_methyl_normalization.R:581-594`) is the only action anywhere in this module that mutates `methyl_dataset$beta`/`$input_scale`/`$source`, following the same explicit-promotion pattern as `mod_preprocessing.R`'s own "Use this as the active dataset" button (module header comment, `mod_methyl_normalization.R:29-34`).

**Audit findings for this tab.** See §7's per-method findings (BMIQ per-sample-failure surfacing, §7 Finding), §14 Finding H-3 (opt-in biological-signal check), and §14 Finding N-1 (`methyl_results$normalization` write is never read).

---

### 3.3 Tab 3: Compare Methods

**Purpose.** Run two or more compatible methods side by side, under the identical filter configuration, and compare their summary statistics and beta-value density distributions.

**Input data.** Same as Method & Run, plus `input$compare_methods_sel` (`checkboxGroupInput`, choices = `available_methods()`, `mod_methyl_normalization.R:461`).

**User inputs.** `compare_methods_sel` (no default selection); `compare_btn` (`actionButton`).

**Reactive dependencies.** `compare_result <- eventReactive(input$compare_btn, ...)` (`mod_methyl_normalization.R:792-805`) — `validate(need(length(sel) >= 2, ...))`; iterates `run_full(sel[i])` for each selected method inside one `withProgress()`.

**Functions called.** The identical `run_full()` orchestration Method & Run uses (`mod_methyl_normalization.R:798`) — Compare Methods is not a separate implementation, it is `n` independent calls to the exact same code path documented in §3.2, so every method-specific function cited there applies here unchanged.

**Statistical operations.** Identical per-method statistics as §3.2/§7, computed once per selected method; no new statistic is introduced by this tab.

**Plots/tables.** `compare_table` (DT: method, samples, probes, mean_variance, mean_correlation, type_ks_stat — `mod_methyl_normalization.R:823-832`); `compare_density_plot` (overlaid density curves — one "Before normalization" reference curve plus one "after" curve per successfully-run method, `mod_methyl_normalization.R:834-852`). Methods that failed to run (e.g., BMIQ selected on a dataset with too few Type I/II-annotated probes) are reported by name and reason in a separate warning line rather than silently omitted (`mod_methyl_normalization.R:813-815`).

**Output objects.** `compare_result()` — a named list (by method key) of `run_full()`'s own return list, each with a `key` field added (`mod_methyl_normalization.R:792-805`). Not persisted beyond the current reactive graph; not downloadable (no `downloadHandler` reads `compare_result()` anywhere in the file — a genuinely absent feature, not a bug, since the task's own Method & Run downloads already cover the single-method case).

**Interpretation.** No `methyl_norm_interpretation()`-style status text is generated for Compare Methods — the table's raw numbers (and the density-plot shapes) are left for the user to compare themselves, unlike Method & Run's own pass/warning/neutral banner.

**Connection to other tabs.** Reads Filters-tab state identically to Method & Run (§2); the two tabs never interact with each other's results (`norm_result()` and `compare_result()` are two independent `eventReactive`s, each gated on its own button).

**Audit findings for this tab.** §14 Finding M-2 (the computed biological-signal-preservation check, `validation$signal_check`, is silently dropped from `compare_table` even when the user configured `group_col_check` on the Filters tab) and §14 Finding L-7 (`compare_density_plot`'s single shared "before" baseline, `res[[1]]$before`, does not account for each compared method's own, potentially slightly different, filtered probe set).

---

## 4. End-to-end pipeline diagram

This reflects the actual implementation verified in §3 above — **not** the illustrative example in the task brief, which does not match this code (there is no separate "input matrix is prepared" step distinct from filtering, and filtering/normalization order is method-type-dependent, not fixed):

```
Dataset tab (mod_methyl_dataset.R, out of scope)
   │  populates methyl_dataset (beta · rg_set · mset · detp · array_type · input_scale · preloaded · source)
   ▼
Normalization tab loads ──► preloaded == TRUE? ──yes──► preloaded_ui(): diagnostics + fixed
   │                                                     "author-normalized" status card. DEAD END —
   │no                                                   no live workflow offered (§3.0).
   ▼
methyl_norm_diagnostics() + methyl_norm_status()  (identical for both pathways, §5)
   │
   ▼
status == "no_bias_detected"?
   │yes, no choice made yet        │yes, "Keep" chosen      │yes, "Re-normalize"/"Compare" chosen, or status != "no_bias_detected"
   ▼                               ▼                         ▼
show 3 choice buttons only    "No changes made" note,   3-tab live workflow (Filters / Method & Run / Compare Methods)
                               dead end                       │
                                                               ▼
                                          Filters tab: probe/sample filter CONFIGURATION only (no computation)
                                                               │
                              ┌────────────────────────────────┴────────────────────────────────┐
                              ▼ (Method & Run: "Run Normalization")                               ▼ (Compare Methods: "Compare Methods")
                    run_full(one method)                                            run_full(method) × N selected, independently
                              │
              ┌───────────────┴────────────────┐
   matrix method (BMIQ/PBC/quantile/none):      raw-intensity method (Noob/Funnorm/SWAN/Dasen/
   sample_scope() → probe filters → method      Stratified quantile/Noob+SWAN/Noob+BMIQ):
              │                                 sample_scope() → method → probe filters (on the
              ▼                                 method's OWN output beta, not the raw input)
   common_probes = intersect(before, after) ◄──────────────────┘
              │
              ▼
   methyl_norm_validation(before, after, anno, group_labels)   ◄── group_labels only if group_col_check set (opt-in)
              │
              ▼
   norm_result() / compare_result()  →  progressive plots/tables, downloads
              │
              ▼ (explicit user action only — "Use this as the active Methylomics dataset")
   methyl_dataset$beta <- after   (methyl_dataset$input_scale relabeled "beta" unless method ∈ {quantile, none})
              │
              ▼
   every Methylomics sub-module below (Celltype, DMP, DMR, WGCNA, Candidates, MR, Coloc, Diagnostic, Biomarker Card)
   now reads the promoted matrix, since they all share the same methyl_dataset reactiveValues object
```

---

## 5. Function-by-function documentation

### 5A. Shiny/UI functions (as used in this module specifically)

- **`tabsetPanel(id = ns("config_tabs"), type = "tabs")`** (`mod_methyl_normalization.R:303-311`) — hosts the 3 live sub-tabs; renders all three bodies up front (§2), which is why Filters state is live for the other two tabs without any explicit cross-tab reactive wiring.
- **`radioButtons(ns("method"), ...)`** (`mod_methyl_normalization.R:423`) — single-select method picker; choices/selected are both reactive (`available_methods()`/`default_method()`), so the option list itself changes when a different dataset is loaded.
- **`conditionalPanel(condition = cond_any(...))`** (`mod_methyl_normalization.R:410,425,431,433,439`) — client-side JS condition strings (built by the small `cond_any()` helper) that show/hide each method's parameter block without a server round-trip.
- **`eventReactive(input$run_btn, ...)` / `eventReactive(input$compare_btn, ...)`** (`mod_methyl_normalization.R:565-573,792-805`) — the module's only two compute-triggering reactives; every other `reactive()` in the file (`anno_result`, `diag_result`, `status_result`, `available_methods`, `default_method`, `sample_scope`, etc.) recomputes automatically whenever its own dependencies change, with no button gate.
- **`shinyjs::toggle()`** (`mod_methyl_normalization.R:668`) — powers the "Show ..." progressive-reveal links; a single `lapply` wires all five toggle sections (`diag_details`, `bvsa`, `stats`, `filtered`, `details`) through one shared observer pattern rather than five separate ones.
- **`downloadHandler()`** ×4 (`mod_methyl_normalization.R:596-622`) — each reads `norm_result()` directly inside its `content` function, so a download always reflects the most recently completed run, not a snapshot taken at click time.
- **`withProgress()` / `incProgress()`** (`mod_methyl_normalization.R:568,795-797`) — user-visible progress messages during `run_full()`/the Compare Methods loop; purely UX, no effect on the computation.

### 5B. Data-processing functions

- **`methyl_design_vector(probe_ids, anno_result)`** (`normalization.R:43-56`) — resolves each probe ID to its Infinium Type I/II design code (`1`/`2`) from manifest annotation; probes absent from the manifest are dropped (`NA`), not guessed. Used by BMIQ and PBC.
- **`build_probe_filters(m)`** (`mod_methyl_normalization.R:396-406`) — reads every enabled Filters-tab checkbox and returns a named list of `keep`/`note` results; combined via a running logical AND in `run_full()` (`mod_methyl_normalization.R:529-530,536-538`).
- **`sample_scope()`** (`mod_methyl_normalization.R:489-509`) — applies the two sample-level filters (missingness, detection rate) and subsets `beta`/`rg_set`/`mset` consistently by column, so every downstream computation for a given run sees the same sample set.
- **`methyl_get_norm_annotation(array_type)`** (`normalization.R:252-282`) — extends the base manifest annotation with `island_relation`/`gene_region` columns from a *different* underlying data object (`Islands.UCSC`/`Other`) than `methyl_get_annotation()` reads, cached separately (`.methyl_norm_anno_cache`, `normalization.R:250`) so the shared QC annotation cache/behavior is untouched.

### 5C. Statistical/normalization functions

Documented per-method in §7 (the task's own required method-by-method audit); shared statistical primitives:
- **`methyl_type_bias_stat(mat, anno_result)`** (`normalization.R:379-400`) — Kolmogorov-Smirnov distance between pooled Type I and Type II beta distributions, on up to 5,000 randomly-sampled probes per type (`normalization.R:391-393`). This is the one quantitative signal both the "already normalized?" heuristic and the before/after validation table are built on.
- **`methyl_norm_status(mat, dataset, anno_result)`** (`normalization.R:424-453`) — see §5D below; a rule-based classifier over `methyl_type_bias_stat()`'s KS statistic (thresholds `0.03`/`0.06`, hard-coded, `normalization.R:440,445`) plus two structural checks (`rg_set` present ⇒ `"raw"`; not on beta scale ⇒ `"unknown"`).
- **`methyl_norm_validation(before, after, anno_result, group_labels)`** (`normalization.R:479-513`) — recomputes mean row variance (`methyl_row_vars()`, `qc.R:22-28`), mean pairwise sample correlation (`methyl_sample_correlation()`, `qc.R:586-594`), the Type I/II KS statistic, and missingness, each before *and* after; optionally an `lm(pc1 ~ group)` R² comparison (§7 discusses the opt-in nature of this last piece).
- **`methyl_norm_interpretation(v)`** (`normalization.R:520-535`) — turns the validation numbers into one of three status strings (`pass`/`warning`/`neutral`) plus a plain-language sentence; the `warning` branch is the only one that can override an otherwise-favorable technical readout, but only fires if `v$signal_check$flagged` is non-`NULL` (§7, §14 Finding H-3).

### 5D. Visualization functions

All visualization functions this module calls are defined in `qc.R` and reused verbatim (not re-implemented) — the same `.methyl_stage_fill` before/after color convention (`qc.R:855-861`) already used by Quality Control's own probe-filtering before/after plots is recognized here for "Before normalization"/"After normalization" stage labels without any change to that shared code:
- `methyl_plot_density()` (`qc.R:863-868`), `methyl_plot_boxplot()` (`qc.R:875-881`) — before/after distribution shape.
- `methyl_plot_scatter2d()` (`qc.R:902-908`) — reused for the before/after PCA facet plot.
- `methyl_plot_corr_heatmap()` (`qc.R:910-919`) — reused for the separate before/after correlation heatmaps.
- `methyl_pca_scores()` (`qc.R:569-580`) and `methyl_sample_correlation()` (`qc.R:586-594`) — the underlying computations behind the two plots above.

---

## 6. What "normalization" means scientifically, tied to this implementation

**Why methylation data require normalization (scientific background).** Beta values are computed as `M / (M + U + 100)` from methylated (`M`) and unmethylated (`U`) probe intensities. Type I and Type II Infinium probe chemistries have different intensity distributions at the same underlying methylation level, background fluorescence and dye bias differ between arrays/scans, and Illumina array processing batches (chip, position, reagent lot, scan date) introduce additional systematic offsets. None of these are biological methylation differences, yet all of them can shift beta values enough to be mistaken for real signal in a downstream comparison.

**Biological vs. technical variation, and why normalization must not remove real signal.** A well-designed normalization method targets a specific, characterized technical artifact (Type I/II design bias, dye bias, background) and leaves the biological methylation signal — differences between disease groups, tissues, cell types, or individuals — intact. **Code evidence** this principle is taken seriously in this implementation: `methyl_norm_validation()`'s optional `signal_check` (`normalization.R:493-511`) exists specifically to catch the failure mode where a normalization method (most plausibly plain quantile normalization, which "makes a strong assumption that samples don't differ substantially in overall methylation" per its own method-info text, `normalization.R:232`) removes real between-group separation along with technical noise, by comparing a chosen group column's PC1 association (`lm(pc1 ~ group)` R²) before vs. after. `methyl_norm_interpretation()` gives this check veto power over an otherwise favorable readout (`normalization.R:521-525`) — but see §14 Finding H-3 for why this safeguard is opt-in rather than automatic.

**Normalization vs. quality control vs. batch correction vs. transformation — how this codebase actually keeps them separate:**
- **Normalization vs. QC:** documented in full in §19; this module and Quality Control's Probe/Sample QC tabs share several literal filter functions (`methyl_filter_missing`/`_snp`/`_sex_chr`/`_cross_reactive`) but compute and apply them completely independently — QC's filtering is report-only, this module's is enforced at run time.
- **Normalization vs. batch correction:** entirely separate modules (§1) — ComBat/RUVm exist only in Quality Control's Batch QC tab.
- **Normalization vs. transformation:** the beta↔M-value logit transform (`methyl_beta_to_mvalue()`, `qc.R:559-562`) is a Quality Control tab feature (its download option), not something this module performs on its own initiative — Normalization methods here either operate on beta values directly (BMIQ/PBC/quantile/matrix path) or derive fresh beta values from raw intensities via `minfi::getBeta()` (every IDAT-based method); no method in this file converts scale as a side effect except quantile/none, which explicitly preserve whatever scale (beta or M) was already loaded (`mod_methyl_normalization.R:584-591`).

**Beta values vs. M-values, as this code treats them.** Beta values (0–1, biologically interpretable as % methylation) are used directly by BMIQ/PBC (both require values strictly inside `(0,1)`, clipped away from the exact boundary at `normalization.R:130,166`, since the beta-mixture/peak-based models are undefined at 0/1). `methyl_norm_diagnostics()` and every method-availability check consult `methyl_dataset$input_scale` to decide whether BMIQ/PBC are even offered (`is_beta_scale()`, `mod_methyl_normalization.R:100`) — they are hidden, with an explanatory note, for an M-value-scale upload (`mod_methyl_normalization.R:421-422`).

---

## 7. Audit of the actual normalization methods implemented

Nine methods plus two "no-op/universal" options are implemented, matching exactly the choices exposed by `METHYL_NORM_METHODS_IDAT`/`_MATRIX`/`_UNIVERSAL`/`_COMBO_SWAN`/`_COMBO_BMIQ` (`mod_methyl_normalization.R:41-48`). No method beyond these eleven exists in the code; none of the UI labels claim a method the underlying function does not perform (spot-checked against each function's own docstring-style comment, below).

| Method (UI label) | Function | Package call | Input | Output |
|---|---|---|---|---|
| Noob | `methyl_norm_noob()` | `minfi::preprocessNoob(offset, dyeMethod)` | `RGChannelSet` | beta via `minfi::getBeta()` |
| Functional normalization | `methyl_norm_funnorm()` | `minfi::preprocessFunnorm(nPCs)` | `RGChannelSet` | beta |
| SWAN | `methyl_norm_swan()` | `minfi::preprocessSWAN(rg, mset)` | `RGChannelSet` + `MethylSet` | beta |
| Dasen | `methyl_norm_dasen()` | `wateRmelon::dasen(mset)` → `minfi::getBeta()` | `MethylSet` | beta |
| Stratified quantile | `methyl_norm_stratified_quantile()` | `minfi::preprocessQuantile()` | `RGChannelSet` | beta |
| BMIQ | `methyl_norm_bmiq()` | `wateRmelon::BMIQ()`, per sample | beta matrix + design vector | beta |
| PBC | `methyl_norm_pbc()` | `ChAMP:::DoPBC()` | beta matrix + design vector | beta |
| Quantile (plain) | `methyl_norm_quantile()` | `limma::normalizeQuantiles()` | beta or M-value matrix | same scale |
| No normalization | inline `list(ok=TRUE, beta=mat_in, ...)` | — | matrix | unchanged matrix |
| Noob + BMIQ | `methyl_norm_noob_bmiq()` | Noob then BMIQ, sequential | `RGChannelSet` | beta |
| Noob + SWAN | `methyl_norm_noob_swan()` | Noob then SWAN, sequential | `RGChannelSet` | beta |

**1. Noob** (`normalization.R:60-67`). Background/dye-bias correction via `minfi::preprocessNoob(rg_set, offset, dyeMethod)`. Statistical principle: normal-exponential deconvolution of the out-of-band signal to estimate and subtract background fluorescence, plus a dye-bias equalization step. Appropriate here because it is the standard first-pass correction whenever raw intensities are available, and — per the method-info text — "on its own it does not address Type I/II probe-design distribution differences" (`normalization.R:217-218`), correctly not oversold as a complete normalization. Parameters `offset` (15) and `dyeMethod` ("single"/ssNoob) are both `minfi::preprocessNoob()`'s own documented arguments, passed through unmodified (`mod_methyl_normalization.R:472`). **Audit: correct, no findings.**

**2. Functional normalization** (`normalization.R:69-76`). `minfi::preprocessFunnorm(rg_set, nPCs)`. Statistical principle: uses the array's built-in control probes, summarized via PCA, as a covariate to regress out unwanted technical variation on top of Noob-style background correction — unlike quantile-based methods, it does not assume similar global methylation across samples, so it is specifically recommended (by minfi's own documentation, cited in the method-info text) "when samples span distinct biological groups or tissues" (`normalization.R:219-220`). **Audit: correct, no findings.**

**3. SWAN** (`normalization.R:78-85`). `minfi::preprocessSWAN(rg_set, mset)`. Statistical principle: subset-quantile within-array normalization — matches beta-value distributions between Type I and Type II probes with similar CpG density, addressing probe-design bias specifically without doing background correction itself (correctly disclosed, `normalization.R:221-222`). Requires both `rg_set` and a pre-existing `mset` — the module supplies `sc$mset`, itself derived once at IDAT-load time via `minfi::preprocessRaw()` (`idat_metrics.R:17`), not re-derived per run. **Audit: correct, no findings.**

**4. Stratified (subset) quantile normalization** (`normalization.R:87-94`). `minfi::preprocessQuantile(rg_set)` — minfi's own re-implementation of Touleimat & Tost (2012). Statistical principle: quantile-normalizes probe subsets stratified by region/probe type. The method-info text correctly flags its own assumption: "Recommended for a single tissue/cell type without large expected global methylation differences between samples; it can distort real, large biological differences if they are present" (`normalization.R:223-224`) — this is the same caveat plain quantile normalization gets, appropriately applied here too since `preprocessQuantile()` is fundamentally a quantile method, just a stratified one. **Audit: correct, no findings.**

**5. Dasen** (`normalization.R:96-109`). `wateRmelon::dasen(mset)`. Statistical principle: separately quantile-normalizes methylated/unmethylated intensities within each probe type, then recombines. **Code evidence — a real API-behavior fix, disclosed in the code's own comment:** `wateRmelon::dasen()` returns a `MethylSet`, not a beta matrix, "despite some of its own documentation examples implying otherwise" — `minfi::getBeta()` is required as an explicit second step, and `as.matrix()` on the raw result errors with no coercion method (`normalization.R:100-103`). This is exactly the kind of implementation detail the task brief asks to be verified rather than assumed from a comment; it was independently checked here by reading `methyl_norm_dasen()`'s own two-step `tryCatch()` chain (`normalization.R:104-108`), which confirms the claim: `wateRmelon::dasen(mset)`'s result is passed through a *second*, separately-error-handled `minfi::getBeta()` call, not treated as already being a matrix. **Audit: correct, no findings.**

**6. BMIQ** (`normalization.R:121-154`). `wateRmelon::BMIQ()`, called **once per sample** in a loop (`normalization.R:145-148`). Statistical principle: fits a three-state beta-mixture model to Type II probes and transforms them onto the Type I distribution (Teschendorff et al. 2013). Requires Type I/II design annotation (`methyl_design_vector()`) — probes without a manifest match are dropped, not guessed (`normalization.R:125-129`). Boundary values are clipped (`normalization.R:130`) since BMIQ's mixture model is undefined at exactly 0 or 1.
  - **`nfit` capping is a correct, disclosed fix, not a bug.** BMIQ samples `nfit` probes *per type* to fit its model; on a small or heavily pre-filtered dataset, the naive default (50,000) can exceed the number of probes of one Infinium type actually present, which would otherwise error. The code caps `nfit_used <- max(1, min(nfit, n_type1, n_type2))` (`normalization.R:140-141`) and cites wateRmelon's own documentation sanctioning a smaller `nfit` as an accuracy/speed tradeoff (`normalization.R:136-139`). **Audit: correct, no findings on this point.**
  - **AUDIT FINDING — HIGH — per-sample BMIQ failures are silently mixed into the promoted matrix without prominent surfacing.** *Location:* `normalization.R:143-153`, consumed at `mod_methyl_normalization.R:581-594` (the promote action). *What the code does:* if `wateRmelon::BMIQ()` errors for a given sample, that sample's *original, unnormalized* beta values are copied through unchanged (`out <- m` initialized before the loop; only successfully-normalized columns are overwritten, `normalization.R:143,147`), and the failure is recorded only in the free-text `note` field ("N sample(s) failed and were left unnormalized: ...", `normalization.R:150-152`). *What it should do, or at minimum surface:* the module's headline "Normalization summary" card (`mod_methyl_normalization.R:632-656`) and its pass/warning/neutral interpretation banner (`methyl_norm_interpretation()`, `normalization.R:520-535`) never inspect `res$note`/`failed_samples` — a per-sample BMIQ failure is visible **only** inside the collapsed-by-default "Show processing details" toggle (`mod_methyl_normalization.R:777-788`) or the downloadable processing-record TXT file (`mod_methyl_normalization.R:611-622`), neither of which is shown by default. *Why it matters:* a user could promote a "BMIQ-normalized" matrix to `methyl_dataset$beta` (feeding every downstream Methylomics sub-module) that actually contains a silent mixture of normalized and unnormalized samples, without ever seeing that fact unless they specifically expand a collapsed section. *Potential scientific consequence:* mixed normalization states within one matrix can introduce a spurious systematic difference between the (accidentally) unnormalized sample(s) and the rest of the cohort — exactly the kind of technical artifact normalization exists to remove — that could be misread as a real biological effect in downstream differential-methylation analysis. *Potential user-facing consequence:* undermines trust in the "Use this as the active dataset" promotion once discovered, and is easy to miss for a first-time user who doesn't expand every collapsed section. *Recommended correction (not applied — documentation only, per this audit's scope):* surface `failed_samples`/`note` prominently in the results-summary card itself (e.g., a dedicated warning banner, not just embedded prose text), and/or block promotion (or require explicit acknowledgment) when `length(failed_samples) > 0`.

**7. PBC (peak-based correction)** (`normalization.R:158-172`). `ChAMP:::DoPBC()` — Dedeurwaerder et al. 2011. Statistical principle: aligns the density peaks of Type I and Type II probe distributions, an alternative to BMIQ/SWAN for the same probe-design-bias problem. Same boundary clipping and design-vector dependency as BMIQ.
  - **AUDIT FINDING — LOW/INFORMATIONAL, disclosed — reliance on an unexported ChAMP internal function.** *Location:* `normalization.R:168` (`ChAMP:::DoPBC(m, ...)`, triple-colon access). *Code says vs. why:* the file's own header explains this choice at length: `champ.norm(method="PBC")` (ChAMP's exported, documented entry point) has "unwanted side effects here (`dir.create()`/`setwd()` into a results folder) that calling `DoPBC` directly avoids," and `DoPBC()` is "the real, tested implementation `champ.norm()` itself calls" (`normalization.R:16-26`) — i.e., this is a disclosed, reasoned tradeoff (avoiding filesystem side effects and a results-directory dependency inappropriate for a multi-user Shiny server) rather than an oversight. *Why it still matters as an audit item:* an unexported function (`:::`) carries no API stability guarantee across ChAMP versions/CRAN or Bioconductor updates — a future ChAMP release could rename or change `DoPBC()`'s signature with no deprecation warning, silently breaking PBC in this app. *Consequence:* none under the currently pinned/installed ChAMP version; a latent upgrade-fragility risk only. *Recommended action:* none mandated by this audit given the disclosed, reasoned tradeoff; worth a version pin or an integration test if ChAMP is ever upgraded.

**8. Quantile normalization (plain)** (`normalization.R:176-183`). `limma::normalizeQuantiles()`. Statistical principle: forces every sample's value distribution to match one reference distribution — the only method here with no Type I/II or raw-intensity requirement, and the only universal fallback besides "none." Its own method-info text is explicit about the tradeoff: "makes a strong assumption that samples don't differ substantially in overall methylation... can remove real signal along with technical noise" (`normalization.R:231-232`) — an honest, non-oversold description. **Audit: correct, appropriately caveated in-app.**

**9. Noob + BMIQ / Noob + SWAN (sequential combos)** (`normalization.R:191-209`). Two-step pipelines: Noob (background/dye correction) then BMIQ or SWAN (probe-design correction) on the result. The file's header explains the scientific rationale for treating these as two separable steps rather than one operation (`normalization.R:185-190`), and each combo's failure handling correctly attributes which step failed ("Noob succeeded but the BMIQ step failed: ...", `normalization.R:195`). **Audit: correct, no findings** beyond BMIQ's own per-sample-failure surfacing gap (Finding above), which propagates unchanged into Noob+BMIQ since it calls the identical `methyl_norm_bmiq()`.

**"No normalization / keep current data"** (`mod_methyl_normalization.R:483`) — a documented pass-through, not disguised as a method; its note text is explicit ("matrix kept exactly as loaded"). Verified: `before`/`after` are numerically identical for this choice once both are restricted to the shared filtered probe set (traced through `run_full()`'s `common_probes` logic, `mod_methyl_normalization.R:544-546`) — the before/after diagnostics for "none" genuinely isolate the effect of the Filters tab alone, with no normalization-attributable change mixed in, which is the scientifically correct behavior for this option.

**Method-parameter passing — verified correct for every exposed parameter.** `noob_offset`/`noob_dye` → `preprocessNoob(offset=, dyeMethod=)` (`normalization.R:64`); `funnorm_npcs` → `preprocessFunnorm(nPCs=)` (`normalization.R:73`); `bmiq_nfit`/`bmiq_nl`/`bmiq_tol` → `BMIQ(nfit=, nL=, tol=)` (`normalization.R:146`) — every UI-exposed numeric input traces to the correct, matching argument name in the underlying package call, with no silent mismatch found.

**AUDIT FINDING — HIGH — the biological-signal-preservation check is opt-in and off by default, yet the "pass" status text implies structural preservation regardless.** *Location:* `group_col_check` selectInput, default `""` (`mod_methyl_normalization.R:381-383`); `methyl_norm_validation()`'s `signal_check` stays `NULL` unless `group_labels` is supplied (`normalization.R:493-494`); `methyl_norm_interpretation()`'s "pass" branch text (`normalization.R:529-531`). *What the code does:* unless a user actively selects a sample-sheet column on the Filters tab, `methyl_norm_interpretation()` can only ever return `"pass"` or `"neutral"` — never `"warning"` — because the only branch that can produce `"warning"` is gated on `isTRUE(v$signal_check$flagged)`, which requires `signal_check` to be non-`NULL` in the first place (`normalization.R:521-525`). When "pass" fires from technical metrics alone (reduced KS statistic and/or increased mean correlation), its text reads: "Normalization completed successfully. QC diagnostics indicate [...] **while preserving overall sample structure**" (`normalization.R:530-531`, emphasis added) — language a user could reasonably read as "biological signal was checked and preserved," when in fact no biological grouping was ever compared. *Why it matters:* the one safeguard this module has against a normalization method silently erasing real disease/control or sex/tissue separation (exactly the risk plain quantile normalization's own method-info text warns about) requires an extra, non-default, easy-to-skip step, and the default-path success message does not indicate that step was skipped. *Potential scientific consequence:* a user who runs Method & Run with default settings and no sample sheet loaded (or a sheet loaded but `group_col_check` left at "(none)") could promote a normalization result that has, in fact, compressed real biological variation, while seeing a green "pass" banner. *Recommended action:* either make the "pass" text conditional on `signal_check` having actually run ("...while reducing technical variation; biological-signal preservation was not checked — set a group column on the Filters tab to verify"), or default `group_col_check` to the first available sample-sheet column when one exists.

---

## 8. Matrix orientation — verified at every stage

**Convention (confirmed, consistent throughout):** `probes × samples` — rows are probes/CpGs, columns are samples, for every matrix this module touches. Verified as follows:

| Function | Input dims | Output dims | Rows | Columns |
|---|---|---|---|---|
| `methyl_parse_matrix()` (`parse_upload.R:11-37`) | file | probes × samples | probe ID (col 1) | sample columns |
| `methyl_idat_derive()` → `minfi::getBeta()` (`idat_metrics.R:13-28`) | RGChannelSet | probes × samples | minfi's own convention (Illumina probe IDs) | array/sample IDs |
| `methyl_filter_missing/_snp/_sex_chr` etc. (`qc.R:30-119`) | matrix | same rows, fewer allowed | `rowMeans()`/probe-indexed logic throughout | untouched |
| `methyl_filter_samples_missingness()` (`normalization.R:331-336`) | matrix | fewer columns | untouched | `colMeans()`-indexed logic |
| Every `methyl_norm_*()` normalization function | probes × samples (or RGChannelSet) | probes × samples | `minfi::getBeta()`'s own row convention preserved | preserved |
| `run_full()`'s `before`/`after` (`mod_methyl_normalization.R:544-546`) | probes × samples | probes × samples, row-subset to `common_probes` | probe IDs via `rownames()` | sample IDs via `colnames()`, unchanged from `sc$mat` |

**No orientation error found.** Specifically checked and cleared:
- **Normalizing samples instead of probes:** every method operates on the array's own probe axis internally (via `minfi`/`wateRmelon`/`ChAMP`/`limma`'s own established conventions for these object types) — no manual transpose exists anywhere in `normalization.R`, confirmed by grepping for `t(` in the file (zero hits).
- **Losing row/column names:** `methyl_norm_pbc()` explicitly re-attaches `rownames`/`colnames` after the `ChAMP:::DoPBC()` call (`normalization.R:170`), the only method whose underlying call does not already preserve them itself — every other method inherits names automatically from `minfi::getBeta()` or, for quantile/BMIQ, from operating in-place on a named input matrix.
- **Silently coercing values / converting numeric to character:** `methyl_parse_matrix()` explicitly sets `storage.mode(m) <- "double"` and validates numeric-ness (`parse_upload.R:28-34`) before this module ever sees the matrix.
- **Accidentally normalizing metadata:** `methyl_dataset$sample_sheet` is read only for column names/values (group labels, filter categories) — never passed into any `methyl_norm_*()` function as data.
- **Mismatching sample order and metadata:** `group_labels_for_check()` (`mod_methyl_normalization.R:511-516`) explicitly re-indexes sheet values by resolved sample ID (`methyl_sheet_sample_ids()`, `qc.R:456-461`) and subsets `[sample_ids]` at the end — not relying on row-order coincidence.

---

## 9. Reactivity and data-dependency audit

**Reactive graph (§4's diagram is the data-flow view; this is the Shiny-reactivity view):**
- `reactive()`s with no button gate, recomputed automatically on any dependency change: `anno_result`, `norm_anno_result`, `manifest_available`, `has_idat`, `is_beta_scale`, `is_illumina` (`mod_methyl_normalization.R:90-103`), `diag_result`, `status_result` (`mod_methyl_normalization.R:105-112`), `available_methods`, `default_method`, `recommendation_text` (`mod_methyl_normalization.R:222-242`), `show_live_workflow` (`mod_methyl_normalization.R:258-262`), `cross_reactive_ids` (`mod_methyl_normalization.R:387-391`), `sample_scope` (`mod_methyl_normalization.R:489-509`).
- `reactiveVal`s: `norm_choice` (`mod_methyl_normalization.R:251`, reset to `NULL` whenever `methyl_dataset$beta` changes, set by four `observeEvent`s on the choice buttons, `mod_methyl_normalization.R:253-256`); `norm_has_run` (`mod_methyl_normalization.R:562`, reset to `FALSE` on a new dataset).
- `eventReactive()`s (the only two genuinely button-gated computations): `norm_result` (`input$run_btn`), `compare_result` (`input$compare_btn`).
- `observeEvent()`s with side effects on shared state: writing `methyl_results$normalization` (`mod_methyl_normalization.R:575-579`); writing `methyl_dataset$beta`/`$input_scale`/`$source` (`mod_methyl_normalization.R:581-594`, the promotion action).

**Does a change in one input automatically recalculate results, or does the user have to click Run?** Both patterns coexist, by design: the diagnostics/status cards and the method picker's *available choices* update live as soon as the dataset or its annotation resolves (no button); the actual normalization computation (and Compare Methods) strictly requires its own explicit button click, and — critically — changing a Filters-tab input *after* a run does **not** retroactively update `norm_result()`; the user must click "Run Normalization" again for the new filter configuration to take effect. This mirrors the QC module's own no-auto-recompute convention (`methylomics_quality_control.md` §10) and was verified the same way: `norm_result`/`compare_result` are `eventReactive`s keyed only on their respective buttons, not on any Filters-tab input.

**Shared reactive objects across the module.** `sample_scope()` and `build_probe_filters()` are the two functions genuinely shared between the Method & Run and Compare Methods code paths (both go through `run_full()`); no other cross-tab reactive dependency exists within this module. `methyl_dataset` and `methyl_results` are the only reactives shared with the rest of the app.

**No caching across Compare Methods' N method runs.** Each `run_full(sel[i])` call independently re-derives `sample_scope()` and re-applies `build_probe_filters()` from scratch (`mod_methyl_normalization.R:796-801`), even though every compared method shares the identical filter configuration and (for matrix methods) an identical pre-normalization input. This is a performance inefficiency, not a correctness bug — flagged as §14 Finding L-8.

---

## 10. Input validation audit

| Mechanism | Where used | Assessment |
|---|---|---|
| `req()` | `anno_result`/`norm_anno_result` (`mod_methyl_normalization.R:91,95`); `diag_result`/`status_result` (`mod_methyl_normalization.R:106,110`); `norm_result` (`method`, `mod_methyl_normalization.R:566`) | Standard Shiny idiom: silently blocks (not errors) until the dependency exists — appropriate for "nothing loaded yet" states. |
| `validate(need(...))` | Dataset presence (`mod_methyl_normalization.R:567`); Compare Methods selection count (`mod_methyl_normalization.R:794`); every `run_full()` failure path (insufficient samples/probes, `mod_methyl_normalization.R:524,542`); every before/after PCA/correlation plot's own minimum-sample-size guard (`mod_methyl_normalization.R:723,731,737`) | Every genuinely-reachable failure mode inside `run_full()` and its downstream plots surfaces a specific, human-readable message rather than a raw R error — good coverage. |
| `tryCatch()` | Every `methyl_norm_*()` method function wraps its own package call (`normalization.R:64,73,82,91,104,106,146,168,180`); `sample_scope()`'s `rg_set`/`mset` subsetting (`mod_methyl_normalization.R:506-507`) | Consistent — no method function can throw an uncaught error into the Shiny session; every failure degrades to `list(ok=FALSE, reason=...)`. |
| Explicit numeric-range validation on user-entered parameters (`noob_offset`, `bmiq_nfit`, `bmiq_nl`, `bmiq_tol`, `funnorm_npcs`) | `numericInput(..., min=, max=, step=)` UI-level bounds only (`mod_methyl_normalization.R:427-437`) | **No server-side re-validation** — a `numericInput`'s `min`/`max` are soft client-side hints in Shiny (a user can still type an out-of-range value, or clear the field to `NA`), and no `req()`/`validate()` in `run_one_method()`/`run_full()` checks these values before passing them straight into the underlying package call (`mod_methyl_normalization.R:472-480`). See §14 Finding L-9. |
| Impossible parameter combinations | `available_methods()`'s construction (§7) already prevents selecting an incompatible method (e.g., BMIQ without manifest annotation) at the UI level — the `method` `radioButtons` choices are always a subset of what's actually runnable. | Sufficient — this is the module's real safeguard against invalid method/data combinations, and it is enforced structurally (choices never include an incompatible option) rather than only by post-hoc validation. |
| Missing-value handling | Every filter function treats `NA`/unannotated probes conservatively (kept, not guessed as failing) — e.g., `methyl_filter_snp()`'s `hit`-gated logic (`qc.R:72-77`), `methyl_design_vector()`'s explicit `NA` for unmatched probes (`normalization.R:52-53`). BMIQ/PBC's Type-vector `NA`s are filtered out entirely before the method runs (`normalization.R:127-129,164-165`). | Consistent, conservative, and disclosed via each filter's own `note` text — sufficient for this module's scope. |

**Are these safeguards scientifically sufficient?** Largely yes for structural/crash-prevention purposes (no code path was found that could feed genuinely incompatible data into a method function). The one gap worth naming is the biological-signal-preservation check being opt-in (§7's High finding) — that is a *scientific* validation gap, not merely a defensive-programming one, since it means the module's one real check on preserved biological signal is not exercised unless a user specifically asks for it.

---

## 11. Output data documentation

| Output | Type | Source | Processing applied | Meaning | Downstream use |
|---|---|---|---|---|---|
| `download_normalized` (CSV) | probe × sample beta/M-value matrix | `norm_result()$after` | Sample-filtered, method-normalized, probe-filtered | The candidate normalized matrix, after this run | Manual re-use outside the app only |
| `download_removed_probes` (CSV) | one-column probe-ID list | `norm_result()$removed_probes` | `setdiff(rownames(sc$mat), common_probes)` | Probes removed by the selected filters for this run | Manual audit/re-use |
| `download_removed_samples` (CSV) | one-column sample-ID list | `norm_result()$removed_samples` | From `sample_scope()`'s filtering | Samples dropped before normalization even ran | Manual audit/re-use |
| `download_processing_record` (TXT) | plain text | `methyl_norm_processing_record()` | Method label, dataset source, timestamp, sample/probe counts, filter notes, method note | Human-readable reproducibility record | Manual re-derivation outside the app (e.g., in a raw `minfi`/`wateRmelon` script) |
| Promoted `methyl_dataset$beta` | probe × sample matrix | `norm_result()$after` (Method & Run only — Compare Methods never promotes) | Identical to `download_normalized`'s content | Becomes the shared "active" dataset | **Every** downstream Methylomics sub-module (Celltype, DMP, DMR, WGCNA, Candidates, MR, Coloc, Diagnostic, Biomarker Card), since they all read the same `methyl_dataset` reactiveValues object |
| `methyl_results$normalization` | small list (`method`, `note`, `n_probes`, `n_samples`) | Written unconditionally once any run completes | None | A summary record | **Not read anywhere in the codebase** — see §14 Finding N-1 |
| `compare_table`/`compare_density_plot` | DT / ggplot | `compare_result()` | Per-method summary statistics | Side-by-side method comparison | Display only — not downloadable, not promotable |

All outputs are **before-or-after clearly labeled** in the UI (valueBoxes explicitly show "N → M" pairs, `mod_methyl_normalization.R:638-640`) — no output was found that presents post-normalization data without indicating that a transformation occurred. Metadata (`sample_sheet`) is never written to by this module and is not included in any of the above outputs.

---

## 12. Reproducibility audit

**Deterministic components (verified):** Every quantile-based, mixture-alignment-free method (Noob, Funnorm, SWAN, Dasen, Stratified quantile, plain quantile) is a deterministic function of its input — no `sample()`/`runif()`/`kmeans()` call exists in any of `methyl_norm_noob/_funnorm/_swan/_dasen/_stratified_quantile/_quantile()` (`normalization.R:60-183`, verified by inspection).

**AUDIT FINDING — MODERATE — no `set.seed()` anywhere in this module, while BMIQ's own fitting procedure is stochastic.** *Location:* `normalization.R:121-154` (`methyl_norm_bmiq`), and by inheritance `normalization.R:158-172` (PBC calls the same design-vector machinery, though `ChAMP:::DoPBC()`'s own internal determinism was **not** independently verified — not determinable from the inspected implementation, since `DoPBC()`'s source is outside this codebase). *Evidence:* `wateRmelon::BMIQ()`'s own documented behavior (cited in this file's header comment, `normalization.R:14-15,132-139`) samples `nfit` probes per Infinium type to fit its three-state beta-mixture model — a random subsample whose selection this app does not seed. Grepped for `set.seed` across every `.R` file in `R/methylomics/`: it is called explicitly in `mod_methyl_celltype.R`, `mod_methyl_diagnostic.R`, and `mod_methyl_featureselection.R` (each for their own stochastic model-fitting), but **not once** in `normalization.R` or `mod_methyl_normalization.R`. *Why it matters:* re-running BMIQ (or Noob+BMIQ) on the identical input, filters, and parameters can produce numerically different normalized beta values from run to run, purely from the different `nfit`-probe subsample each run draws — a genuine reproducibility gap for the one method family in this module whose fitting procedure is not analytically deterministic. *Consequence:* a user who downloads the "Processing details" record and later re-runs the identical configuration for an audit trail or publication supplementary material cannot expect bit-identical output for BMIQ/Noob+BMIQ, even though the record's parameter list would suggest otherwise. *Recommended action:* call `set.seed()` with a fixed or user-configurable seed immediately before the BMIQ per-sample loop, and record the seed value in `methyl_norm_processing_record()`'s output.

**Plotting-only, non-data-affecting stochasticity (Low/Informational, does not affect the promoted matrix or any statistic).** `beta_long_before_after()` (`mod_methyl_normalization.R:696-708`) and `compare_density_plot` (`mod_methyl_normalization.R:834-852`) both subsample up to 20,000/15,000 probes via unseeded `sample()` calls purely for plot rendering — the *displayed* density/boxplot curves' exact probe subset changes on every re-render of the same `norm_result()`, though the underlying values, statistics, and downloads are unaffected (no `sample()` call exists anywhere in the actual normalization or validation computation path).

**Package versions.** Not pinned anywhere in the inspected files (no `renv.lock`/`DESCRIPTION` version constraint was part of this audit's scope) — **not determinable from the inspected implementation** whether results are version-locked across deployments; this is consistent with the equivalent gap already noted for Quality Control.

**Deterministic parameter defaults.** All confirmed literal, hard-coded defaults, not derived from data (`normalization.R`'s function signatures and `mod_methyl_normalization.R:427-437`'s `numericInput` defaults) — no default is itself randomly initialized.

**Session state / saved outputs.** Neither `norm_result()` nor `compare_result()` persists across a browser session reload (standard Shiny reactive-value lifetime); the only durable artifacts are the four manual downloads and, if promoted, the mutated `methyl_dataset$beta` for the remainder of that session.

---

## 13. Normalisation vs. Quality Control — the module boundary, verified in code

**What QC produces that Normalization consumes.** **Nothing, at the data level.** Quality Control's filtered matrix (`probe_qc_result()$filtered`, documented in `methylomics_quality_control.md` §3.3) is never read by this module — grepping `mod_methyl_normalization.R` and `normalization.R` for `probe_qc_result`/`sample_qc_result`/any QC-tab-specific reactive name returns zero hits. Both modules read only the same upstream `methyl_dataset$beta`, independently.

**Is QC filtering actually enforced before normalization runs? No — code-verified.** Quality Control's own filtering is confirmed (in `methylomics_quality_control.md` §1/§14) to never mutate `methyl_dataset$beta` — it is a report/export-only tool. This was independently re-confirmed for this audit by grepping `mod_methyl_qc.R` for `methyl_dataset$beta\s*<-` / `promote` / "active dataset": **zero matches**. Since the *only* code path anywhere in the Methylomics group that writes `methyl_dataset$beta` is this Normalization module's own `promote_btn` handler (`mod_methyl_normalization.R:583`), a user can navigate directly to Normalization and run any compatible method on the completely unfiltered, as-loaded dataset — with no probe-level or sample-level QC ever having been applied, and no gate in this module's own code that checks whether Quality Control has been run.

**Does Normalization independently perform QC-shaped filtering? Yes — and this is the actual (non-obvious) mechanism by which any filtering happens before a promoted normalization.** The Filters tab (§3.1) reimplements — by calling the *same shared functions* (`methyl_filter_missing`, `methyl_filter_snp`, `methyl_filter_sex_chr`, `methyl_filter_cross_reactive`, all imported unchanged from `qc.R`) plus three Normalization-tab-only additions (chromosome exclusion, island-relation, gene-region) — a second, independent probe/sample filtering pass, entirely separate from whatever (if anything) was configured on the Quality Control tab. Because both tabs read from `methyl_dataset$beta` independently and neither writes back to it except via Normalization's own promotion, **a user's Quality Control filter selections have zero effect on what Normalization actually filters** unless the user manually re-selects equivalent filters on the Normalization tab's own Filters sub-tab.

**Can normalization be run on unfiltered data? Yes, structurally always possible**, and it is in fact the default path: if a user opens Normalization, picks a method, and clicks Run without ever touching the Filters tab, every `f_*` checkbox defaults to `FALSE` (`mod_methyl_normalization.R:335,338,343,356,357,360,361,367,372` — every filter checkbox's `value = FALSE`), so `build_probe_filters()`/`sample_scope()` return effectively no-op filter lists and `run_full()` operates on the complete, as-loaded probe/sample set.

**Is this a scientific risk? A genuine one, moderate severity, disclosed here rather than in either module's own UI.** Running BMIQ/SWAN/etc. on data still containing failed probes (low detection rate, SNP-overlapping, cross-reactive) does not make the normalization *computation* itself incorrect — every method function operates correctly on whatever probe/sample set it receives — but it means the normalized result can still carry forward quality problems a QC pass would have caught, and any Compare Methods statistics (mean variance, mean correlation, KS statistic) are then computed over a probe/sample set that has not been quality-screened, regardless of which normalization method looks "best" by those metrics.

**Is there duplicated processing?** Yes, by design rather than accident — `methyl_filter_missing`/`_snp`/`_sex_chr`/`_cross_reactive` are literally the same function objects in both modules (imported from the same `qc.R` definitions, not two independently-written copies with potential drift), so there is no risk of the *filter logic itself* silently diverging between QC and Normalization the way, for example, `methyl_guess_batch_column()` was found to be duplicated-and-drift-risked within QC alone (`methylomics_quality_control.md` §14 Finding L-5). The duplication here is in *invocation* (each tab independently decides whether/when to apply these filters, with separate UI state), not in *implementation*.

**AUDIT FINDING — MODERATE — no cross-module signal indicating whether Quality Control has been run, in either direction.** *Location:* absence of any read of a QC-tab reactive/flag anywhere in `mod_methyl_normalization.R` (verified by exhaustive grep, above), and — symmetrically — no equivalent read of Normalization's own `norm_has_run`/status in `mod_methyl_qc.R`. *Why it matters:* neither module can tell a user "you have not yet run Quality Control" or "you are about to normalize unfiltered data" — the only signal Normalization gives is its own generic Type I/II bias status card, which (per `methyl_norm_status()`'s own deliberately-narrow scope, `normalization.R:402-423`) says nothing about probe/sample quality. *Consequence:* a user following a linear "QC tab, then Normalization tab" workflow (the order the tabs are listed in, `submodules_registry.R:40-41`) could reasonably assume their QC choices carry forward, when they do not. *Recommended action:* a lightweight, informational note on the Normalization Filters tab surfacing whether Quality Control's Probe QC has been run this session (e.g., reading a shared flag, without necessarily forcing any filter to be re-applied automatically) — no code change is made here, this is a documentation-and-audit finding only, per this task's scope restriction.

---

## 14. Data-leakage audit

**Not a machine-learning workflow — the classical train/test data-leakage concept does not directly apply**, and this audit does not manufacture a leakage finding where none exists. Normalization parameters here are the fitted state internal to each `minfi`/`wateRmelon` function call (e.g., BMIQ's per-sample beta-mixture fit, quantile normalization's reference distribution) — every method is fit and applied to the **same** set of samples in one pass; there is no train/test split, no held-out validation set, and no parameter reuse across separate sample subsets anywhere in this module.

**The closest analogous concern — and it is answered cleanly.** Because BMIQ fits its model *per sample independently* (the loop at `normalization.R:145-148` calls `wateRmelon::BMIQ()` once per column, with no cross-sample parameter sharing), there is no possibility of one sample's normalization "leaking" information from another sample's data — each sample's fit uses only its own beta values. Noob/Funnorm/SWAN/Dasen/quantile methods, by contrast, do pool information *across* samples in the same run (quantile normalization by construction, Funnorm via its cross-sample control-probe PCA) — this is expected, standard normalization-method behavior (these methods are explicitly designed to use the whole batch), not a leakage bug, since there is no train/test boundary being crossed.

---

## 15. Common methylomics-error checklist — findings only

Checked against the task's own list; items not flagged were verified correct and are not repeated as separate findings elsewhere:

- **Incorrect beta-value handling:** not found — boundary clipping (BMIQ/PBC) is correct and standard.
- **Incorrect M-value conversion:** not applicable — this module never converts scale itself (§6); no finding.
- **Inappropriate log transformation:** not applicable — none performed by this module.
- **Incorrect matrix orientation:** not found (§9).
- **Normalizing metadata:** not found (§9).
- **Loss of CpG/sample identifiers:** not found, one method (PBC) requires and receives an explicit fix (§9).
- **Incorrect NA handling:** not found — conservative throughout (§10).
- **Zero-variance / duplicate probes:** duplicate probe IDs are already rejected at upload time (`methyl_parse_matrix()`, `parse_upload.R:21-26`, out of this module's own scope but upstream-verified); zero-variance probes are not specifically filtered by this module (no `methyl_filter_variance`/`_sd` equivalent exists on the Normalization Filters tab, unlike Probe QC, which does offer both) — **not a bug**, since normalization methods here do not require non-zero variance to run correctly, but worth noting this specific QC filter type is not duplicated onto Normalization's own Filters tab the way missingness/SNP/sex-chr/cross-reactive are.
- **Inappropriate normalization for array type:** not found — method availability is correctly gated by `array_type`/manifest-annotation availability throughout (§7).
- **Batch effects confused with normalization:** not found — cleanly separated modules (§1, §6).
- **Biological groups unintentionally removed:** the real, flagged risk — §7's High finding (opt-in signal-preservation check).
- **Parameters in UI but unused in analysis, or vice versa:** not found — every exposed numeric input traces to a matching underlying-package argument (§7); no undisclosed hidden parameter was found in any `methyl_norm_*()` signature that lacks a corresponding UI control, except the deliberately-not-exposed defaults on SWAN/Dasen/Stratified-quantile/PBC/quantile (which genuinely have no tunable parameters in their respective package functions).
- **Labels not matching algorithms:** not found — every method-info text (`METHYL_NORM_METHOD_INFO`, `normalization.R:216-237`) was checked against its function's actual package call and found accurate, including the deliberately-hedged claims (e.g., Noob "does not address Type I/II... on its own").
- **Inconsistent preprocessing between preloaded and uploaded datasets:** by design, and disclosed (§3.0) — the preloaded dataset is never reprocessed through this module at all.

---

## 16. Function inventory table

| Function | File | Type | Used by | Purpose | Input | Output |
|---|---|---|---|---|---|---|
| `methyl_design_vector` | normalization.R:43-56 | data-processing | BMIQ, PBC | resolve Type I/II design per probe | probe IDs, `anno_result` | keep-vector + design codes |
| `methyl_norm_noob` | normalization.R:60-67 | normalization method | Method & Run, Compare Methods | background/dye correction | `RGChannelSet` | beta + note |
| `methyl_norm_funnorm` | normalization.R:69-76 | normalization method | same | control-probe-based normalization | `RGChannelSet` | beta + note |
| `methyl_norm_swan` | normalization.R:78-85 | normalization method | same | Type I/II distribution matching | `RGChannelSet` + `MethylSet` | beta + note |
| `methyl_norm_stratified_quantile` | normalization.R:87-94 | normalization method | same | stratified quantile norm | `RGChannelSet` | beta + note |
| `methyl_norm_dasen` | normalization.R:96-109 | normalization method | same | wateRmelon default method | `MethylSet` | beta + note |
| `methyl_norm_bmiq` | normalization.R:121-154 | normalization method | same | per-sample beta-mixture correction | beta matrix, design vector | beta + note + `failed_samples` |
| `methyl_norm_pbc` | normalization.R:158-172 | normalization method | same | peak-based correction | beta matrix, design vector | beta + note |
| `methyl_norm_quantile` | normalization.R:176-183 | normalization method | same | universal quantile baseline | any-scale matrix | same-scale matrix + note |
| `methyl_norm_noob_bmiq` | normalization.R:191-197 | sequential combo | same | Noob then BMIQ | `RGChannelSet` | beta + note |
| `methyl_norm_noob_swan` | normalization.R:199-209 | sequential combo | same | Noob then SWAN | `RGChannelSet` | beta + note |
| `methyl_get_norm_annotation` | normalization.R:252-282 | data-processing | Filters tab | extend base annotation with island/region columns | array type | annotation data.frame |
| `methyl_filter_island_relation` | normalization.R:290-300 | probe filter | Filters tab | keep selected island-relation categories | matrix, annotation, categories | keep/note |
| `methyl_filter_gene_region` | normalization.R:302-312 | probe filter | Filters tab | keep selected gene-region categories | matrix, annotation, regions | keep/note |
| `methyl_filter_chromosome` | normalization.R:314-323 | probe filter | Filters tab | exclude named chromosome(s) | matrix, annotation, chr list | keep/note |
| `methyl_filter_samples_missingness` | normalization.R:331-336 | sample filter | Filters tab (via `sample_scope`) | drop sparse samples pre-normalization | matrix, threshold | keep/note |
| `methyl_norm_diagnostics` | normalization.R:343-367 | statistical | both pathways | dataset summary | matrix, dataset, annotation | diagnostics list |
| `methyl_type_bias_stat` | normalization.R:379-400 | statistical | status, validation | Type I/II KS distance | matrix, annotation | KS stat + medians |
| `methyl_norm_status` | normalization.R:424-453 | statistical/classifier | both pathways | raw/bias-detected/no-bias/unknown | matrix, dataset, annotation | status + message |
| `methyl_norm_recommendation` | normalization.R:458-472 | advisory text | uploaded pathway | non-binding method suggestion | dataset, status, methods | text |
| `methyl_norm_validation` | normalization.R:479-513 | statistical | Method & Run, Compare Methods | before/after metrics + optional signal check | before, after, annotation, group labels | validation list |
| `methyl_norm_interpretation` | normalization.R:520-535 | statistical/classifier | Method & Run | pass/warning/neutral verdict | validation list | status + text |
| `methyl_norm_processing_record` | normalization.R:538-550 | reporting | Method & Run downloads | reproducibility text record | run metadata | formatted string |
| `available_methods` | mod_methyl_normalization.R:222-230 | data-processing | Method & Run, Compare Methods | compatible-method list | dataset state | named vector |
| `default_method` | mod_methyl_normalization.R:232-238 | data-processing | Method & Run | pre-selected method | available methods, status | method key |
| `build_probe_filters` | mod_methyl_normalization.R:396-406 | data-processing | `run_full` | assemble enabled probe filters | matrix | named filter-result list |
| `sample_scope` | mod_methyl_normalization.R:489-509 | data-processing | `run_full` | apply sample-level filters | dataset | filtered mat/rg/mset |
| `run_one_method` | mod_methyl_normalization.R:470-486 | dispatch | `run_full` | method-key → function call | method key, inputs | method result |
| `run_full` | mod_methyl_normalization.R:522-558 | orchestration | Run button, Compare Methods | one complete filter→method→validate run | method key | full result list |
| `group_labels_for_check` | mod_methyl_normalization.R:511-516 | data-processing | `run_full` | resolve group labels for signal check | sample IDs | named character vector |

*(Shared, imported-unchanged functions from `qc.R`/`annotation.R`/`parse_upload.R` — `methyl_filter_missing`, `methyl_filter_snp`, `methyl_filter_sex_chr`, `methyl_filter_cross_reactive`, `methyl_sample_failed_probe_pct`, `methyl_sheet_sample_ids`, `methyl_row_vars`, `methyl_pca_scores`, `methyl_sample_correlation`, `methyl_plot_density`, `methyl_plot_boxplot`, `methyl_plot_scatter2d`, `methyl_plot_corr_heatmap`, `methyl_get_annotation`, `methyl_parse_probe_list` — are documented in full in `methylomics_quality_control.md` §5 and are not re-documented here beyond their call sites already cited above.)*

---

## 17. Tab inventory table

| Tab | Purpose | Inputs | Main functions | Processing | Outputs | Depends on |
|---|---|---|---|---|---|---|
| *(preloaded, not a tab)* | Show fixed author-normalized status | none | `methyl_norm_diagnostics`, `methyl_norm_status` | descriptive only, no computation offered | value boxes, technical note | Dataset tab's preloaded load |
| Filters | Configure optional probe/sample filters | 11 filter checkboxes + thresholds, 1 group-check selector | `build_probe_filters` inputs only (no computation on this tab) | none — configuration only | none directly | none |
| Method & Run | Pick, configure, and run one normalization method | method radio + per-method numeric/radio params, Run button | `run_full`, `run_one_method`, every `methyl_norm_*` | filter→method→filter (order by method type)→validate | `norm_result()`, 4 downloads, promotion | Filters tab (live) |
| Compare Methods | Run N≥2 methods under identical filters | checkbox-group method selection, Compare button | `run_full` × N | identical to Method & Run, repeated per method | comparison table + density plot | Filters tab (live) |

---

## 18. Package and dependency audit

| Package | Function(s) used | Explicit namespace? | Actually loaded/attached elsewhere? | Scientific role |
|---|---|---|---|---|
| `minfi` | `preprocessNoob`, `preprocessFunnorm`, `preprocessSWAN`, `preprocessQuantile`, `getBeta` | Yes, every call (`minfi::`) | Not `library()`-attached globally (deliberate, per `annotation.R:39-47`'s own explanation of why — attaching would mask `base::strsplit()` app-wide via `Biostrings`); accessed only via explicit namespace throughout | Raw-IDAT-based normalization methods |
| `wateRmelon` | `dasen`, `BMIQ` | Yes, every call | Not attached globally | Dasen and BMIQ implementations |
| `ChAMP` | `:::DoPBC` (unexported) | Yes (`ChAMP:::`) | Not attached globally | PBC implementation (§7 Finding, disclosed) |
| `limma` | `normalizeQuantiles` | Yes (`limma::`), though redundant | **Yes** — `library(limma)` at `global.R:74` | Plain quantile normalization baseline |
| `matrixStats` | `rowVars` (via shared `methyl_row_vars`, `qc.R:22-28`) | Yes, guarded by `requireNamespace` with a base-R fallback | Not attached globally | Fast row-variance computation for validation metrics |
| `tidyr` | `pivot_longer` (`mod_methyl_normalization.R:699,842`) | Implicit (no `tidyr::` prefix — relies on the package being attached elsewhere in the app) | **Not determinable from the inspected implementation** whether `library(tidyr)` is called in `global.R` — not searched as part of this scoped audit; flagged only as an observation, not a finding, since it did not error in any code path reviewed | Reshaping before/after beta values into long format for `ggplot2` |
| `DT` | `datatable`, `formatRound` (`mod_methyl_normalization.R:755,769,773,830-831`) | Yes (`DT::`) | — | Interactive result tables |
| `ggplot2` | `aes`, `geom_density`, `scale_color_manual`, `labs`, `facet_wrap` (`mod_methyl_normalization.R:726,848-851`) | Implicit (unqualified `ggplot`/`aes`/etc. calls, and via the shared `qc.R` plot builders which are also unqualified) | Consistent with the rest of the app's plotting code (out of this module's own scope to verify globally) | Compare Methods density overlay, before/after PCA facet |
| `shinyjs` | `useShinyjs`, `hidden`, `toggle`, `show` (`mod_methyl_normalization.R:60,74,680,684,693,694,668`) | Yes (`shinyjs::`) | — | Progressive-reveal toggle sections |
| `shinycssloaders` | `withSpinner` (`mod_methyl_normalization.R:61,308,309,310,675,676,680,684,818`) | Implicit (unqualified `withSpinner`) | — | Loading spinners on async-feeling sections |

**No namespace conflicts found.** Because `minfi`/`wateRmelon`/`ChAMP` are deliberately never `library()`-attached (only accessed via `::`/`:::`), there is no risk of these three heavy Bioconductor packages' own function-masking side effects reaching the rest of the app — the same design rationale `annotation.R`'s header already documents for the shared annotation-loading code (`annotation.R:39-47`) is consistently honored in every normalization-method function in this file.

**No deprecated-function usage was identified** in the functions actually called (all are current, documented entry points of their respective packages as of the versions implied by the code's own citations).

**Package versions are not pinned** in any file this audit inspected — see §12.

---

## 19. Thesis Implementation Paragraph

> The Methylomics Normalization sub-module (`mod_methyl_normalization.R`/`normalization.R`, registered as `"normalization"` in the "Data" group, immediately after Quality Control) presents either a fixed, non-recomputed status card for the app's preloaded, author-normalized whole-blood cohort, or, for any uploaded dataset, three live tabs — Filters, Method & Run, and Compare Methods. It receives the shared `methyl_dataset` object populated by the Dataset tab: a probe-by-sample beta- or M-value matrix, and, where a raw-IDAT upload was used, the underlying `RGChannelSet`/`MethylSet`/detection-p-value data. Eleven method choices are implemented and offered only when their data requirements are actually met: `minfi::preprocessNoob`, `preprocessFunnorm`, `preprocessSWAN`, and `preprocessQuantile`, and `wateRmelon::dasen` for raw-intensity input; `wateRmelon::BMIQ` and ChAMP's peak-based correction for beta-value input with Type I/II probe-design annotation; `limma::normalizeQuantiles` as a universal, scale-preserving baseline; and two sequential Noob-then-BMIQ/SWAN combination workflows. These methods matter because they target scientifically distinct sources of technical variation in Illumina methylation-array data — background fluorescence and dye bias (Noob), and residual Type I/II probe-design distributional bias (BMIQ, SWAN, PBC, Dasen, stratified quantile) — rather than treating "normalization" as one interchangeable operation, and each choice is quantified before and after with row variance, sample-sample correlation, and the Kolmogorov-Smirnov Type I/II distributional statistic that also drives the module's own automatic "already normalized?" status detection. The resulting candidate matrix is promoted to the shared active dataset only by explicit user action, never automatically, supporting every downstream Methylomics sub-module — cell-type deconvolution, differential-methylation, region-based, and network analyses — that reads from the same shared dataset object.

### 19.1 Thesis Methods-section register

> Methylation normalization was implemented as a dedicated sub-module of the Methylomics pipeline, positioned immediately downstream of quality control. For the study's preloaded whole-blood cohort, the module deliberately performed no re-normalization: because this dataset had already been normalized by its original authors prior to deposition, reprocessing it independently was judged more likely to introduce a systematic discrepancy from the published values than to improve upon them, and the module instead reported descriptive diagnostics and a fixed provenance statement. For any user-supplied dataset, normalization proceeded through three stages — optional probe- and sample-level filtering, selection and execution of a normalization method, and, where desired, direct comparison of multiple methods under identical filtering conditions. Eleven normalization strategies were made available, restricted at each stage to those methods whose data requirements were actually satisfied by the uploaded material: background and dye-bias correction (Noob), control-probe-based functional normalization, and two probe-design correction methods requiring raw signal intensities (SWAN, Dasen) or Type I/II manifest annotation alone (BMIQ, peak-based correction), together with a stratified quantile method, a scale-preserving quantile baseline applicable to any input, and two sequential Noob-initiated workflows. This selection reflected the recognized distinction in the methylation-array literature between background/dye-bias correction and correction of the systematic distributional difference between Infinium Type I and Type II probe designs, rather than treating normalization as a single undifferentiated procedure. The effect of each applied method was quantified by comparing probe-wise variance, inter-sample correlation, and the Kolmogorov–Smirnov distance between Type I and Type II beta-value distributions before and against after normalization, with an optional check of whether a user-defined biological grouping's association with the leading principal component was preserved. A normalized matrix was adopted as the dataset used by all subsequent stages of the Methylomics pipeline — cell-type deconvolution, differential methylation, regional analysis, and network-based modeling — only upon explicit confirmation by the user, ensuring that no normalization result altered downstream analyses without deliberate, auditable approval.

### 19.2 Thesis Methods-section register, with input/output stated explicitly

> Methylation normalization was implemented as a dedicated sub-module of the Methylomics pipeline, positioned immediately downstream of quality control. **Input** to the module was the probe-by-sample methylation matrix held in the shared dataset object — beta values (proportion methylated, bounded 0–1) or M-values (logit-transformed), together with, for a raw-array upload, the underlying green/red channel signal intensities, methylated/unmethylated intensity sets, and detection p-values; optional sample-level metadata supplied a grouping variable against which biological-signal preservation could be checked. For the study's preloaded whole-blood cohort, this input was already normalized by its original authors prior to deposition, and the module deliberately performed no re-normalization, reporting descriptive diagnostics and a fixed provenance statement instead. For any user-supplied dataset, the module offered optional probe- and sample-level filtering followed by one of eleven normalization strategies, applied only where the input actually met that method's requirements: background and dye-bias correction, control-probe-based functional normalization, and probe-design correction methods operating on either raw signal intensities or on beta values with Type I/II manifest annotation, together with a scale-preserving quantile baseline and two sequential combination workflows. This selection reflected the distinction, established in the methylation-array literature, between background/dye-bias correction and correction of the systematic distributional difference between Infinium Type I and Type II probe designs. **Output** was a normalized probe-by-sample beta matrix of identical structure to the input, accompanied by a quantitative before/after comparison — probe-wise variance, inter-sample correlation, the Kolmogorov–Smirnov distance between Type I and Type II beta-value distributions, and, where a grouping variable was supplied, its association with the leading principal component — together with the list of probes and samples removed by filtering and a plain-text processing record. This normalized matrix was adopted as the active dataset for all subsequent stages of the Methylomics pipeline — cell-type deconvolution, differential methylation, regional analysis, and network-based modeling — only upon explicit confirmation by the user, ensuring that no normalization result altered downstream analyses without deliberate, auditable approval.

---

## 20. XomicShiny-style implementation paragraph

> Normalization: raw or as-uploaded methylation data were optionally normalized using one of eleven implemented methods (Noob, Functional normalization, SWAN, Dasen, stratified quantile, BMIQ, peak-based correction, plain quantile, Noob+SWAN, Noob+BMIQ, or no normalization), selected according to data availability (raw intensity vs. beta-value-only, with or without Type I/II probe-design annotation). Background and dye-bias correction (Noob) and probe-design bias correction (BMIQ/SWAN/PBC/Dasen/stratified quantile) were treated as distinct correction targets, consistent with their respective methodological literature. Normalization performance was assessed by row-wise variance, mean sample-sample correlation, and the Kolmogorov-Smirnov statistic between Type I and Type II probe beta-value distributions, computed before and after each run; an optional check compared a user-selected sample-grouping variable's association with the first principal component before and after normalization, to flag potential loss of biological signal. Normalized beta values were only adopted as the active dataset upon explicit user confirmation.

---

## 21. Short tab-based thesis paragraph

> The Normalization submodule contains three live tabs for an uploaded dataset (none for the preloaded, author-normalized cohort, which is shown a fixed status instead). The Filters tab accepts optional probe- and sample-level filter configuration and applies nothing by itself. The Method & Run tab accepts a chosen normalization method and its parameters, applies `run_full()` (filter → method → filter, order dependent on method type), and produces a normalized beta matrix plus before/after diagnostics, which the user may promote to the shared active dataset. The Compare Methods tab applies the identical `run_full()` procedure to two or more selected methods under the same filter configuration and produces a side-by-side comparison table and density plot. Together, these steps generate the promoted normalized matrix used by every downstream Methylomics sub-module.

---

## 22. Audit Summary

### Correctly implemented

- Matrix orientation (probes × rows, samples × columns) is consistent and correct throughout every function this module calls (§9).
- Every normalization method's exposed parameter maps correctly to its underlying package argument, with no mislabeled UI text found (§7).
- Method availability is correctly and defensively gated by actual data/annotation presence — no path was found where an incompatible method could be selected (§7, §10).
- The "no normalization" pass-through option genuinely leaves values unchanged, verified via the shared `common_probes` restriction logic (§7).
- The before/after diagnostic comparison genuinely isolates the normalization method's own effect from the Filters tab's effect, for every method type, including the matrix-vs-raw-intensity filter-ordering difference (§7, verified in detail rather than assumed).
- Boundary handling for BMIQ/PBC (clipping away from exact 0/1) and the `nfit` auto-capping for small datasets are both correct, disclosed, and appropriately reasoned (§7).
- Dasen's `MethylSet`-not-beta-matrix API quirk is correctly handled with an explicit second extraction step (§7).
- Failure handling is comprehensive and consistent — every method function degrades to `list(ok=FALSE, reason=...)` rather than throwing (§10).
- `minfi`/`wateRmelon`/`ChAMP` are deliberately never globally attached, avoiding namespace-masking risk to the rest of the app (§18).
- The promotion mechanism (explicit button, never automatic) correctly follows the same reversible-preview pattern used elsewhere in the app (§3.2).

### Potential improvements

- `bmiq_nfit`/`bmiq_nl`/`bmiq_tol`/`funnorm_npcs`/`noob_offset` have no server-side range re-validation beyond the `numericInput` client-side hint (§10, Finding L-9).
- Compare Methods independently re-derives `sample_scope()`/`build_probe_filters()` per method with no caching (§9, Finding L-8).
- `compare_density_plot` uses one shared "before" baseline (`res[[1]]$before`) rather than each compared method's own filtered "before" (§3.3, Finding L-7).
- `methyl_results$normalization` is written but never read anywhere in the codebase (§1, §11, Finding N-1).
- Filter checkboxes for SNP/sex-chromosome/chromosome-exclusion remain checkable (and silently no-op) for array types with no manifest annotation, such as EPICv2/Custom array, rather than being disabled (§3.1, Finding L-6).

### Scientific concerns

- **The biological-signal-preservation check is opt-in and off by default, yet the module's own "pass" status text reads as if structural preservation were verified regardless (§7, High).** This is the most consequential finding in this audit: it affects the trustworthiness of the module's own success messaging, not merely a UI nicety.
- **Quality Control's filtering is never enforced before Normalization, and Normalization independently reimplements an overlapping-but-separately-configured filter set (§13, Moderate).** Not a bug in either module individually, but a real, code-verified gap in the intended pipeline order the tab registration order (`submodules_registry.R:40-41`) otherwise suggests.
- **No `set.seed()` anywhere in this module despite BMIQ's stochastic per-sample model fit (§12, Moderate).** A genuine reproducibility gap for one specific method family.

### High-severity findings

```
Finding: Per-sample BMIQ failures are silently mixed into the promoted matrix without prominent surfacing.
Location: normalization.R:143-153 (methyl_norm_bmiq); consumed at mod_methyl_normalization.R:581-594.
Evidence: A failed sample's ORIGINAL, unnormalized beta values are copied through unchanged; the failure
is recorded only in a free-text note visible solely inside a collapsed "Show processing details" toggle
or a downloadable TXT record — never in the headline summary card or the pass/warning/neutral verdict.
Why it matters: A user can promote a "BMIQ-normalized" matrix that actually mixes normalized and
unnormalized samples, feeding every downstream Methylomics sub-module, without seeing this fact by default.
Potential consequence: A spurious systematic difference between the unnormalized sample(s) and the rest of
the cohort, mistakable for real biological signal in downstream differential-methylation analysis.
Recommended action: Surface failed_samples/note prominently in the results-summary card; consider requiring
explicit acknowledgment (or blocking promotion) when failures occurred. (Documentation-only per this audit's
scope — no code was changed.)
```

```
Finding: The biological-signal-preservation check is opt-in (off by default), but the "pass" status text
implies overall structural preservation was verified regardless of whether the check ran.
Location: group_col_check default "" at mod_methyl_normalization.R:381-383; the gating logic at
normalization.R:493-494,521-525; the "pass" text at normalization.R:529-531.
Evidence: methyl_norm_interpretation() can only return "warning" if v$signal_check is non-NULL, which
requires the user to have selected a sample-sheet group column on the Filters tab; the "pass" branch text
("...while preserving overall sample structure") fires from technical metrics alone and does not
distinguish "checked and preserved" from "never checked."
Why it matters: This is the module's one real safeguard against normalization erasing genuine biological
signal (a risk its own quantile-normalization method-info text explicitly names), and it can be silently
skipped by simply not touching one selectInput.
Potential consequence: A promoted normalization result that has, in fact, compressed real biological
variation, accompanied by a green "pass" banner that does not disclose the check was skipped.
Recommended action: Make the "pass" text conditional on signal_check having actually run, or default
group_col_check to an available sample-sheet column. (Documentation-only per this audit's scope.)
```

---

*End of document. Prepared by direct code inspection of `mod_methyl_normalization.R` (855 lines, read in full) and `normalization.R` (551 lines, read in full), plus targeted reads of `qc.R`, `annotation.R`, `parse_upload.R`, `mod_methyl_dataset.R`, `idat_metrics.R`, `submodules_registry.R`, `global.R`, and `mod_methyl_qc.R` (for the QC↔Normalization boundary in §13). Scope was strictly limited to Methylomics → Normalization per the task's own instructions; no application code was modified, and no UI, styling, navigation, or other module was altered in the course of this audit.*
