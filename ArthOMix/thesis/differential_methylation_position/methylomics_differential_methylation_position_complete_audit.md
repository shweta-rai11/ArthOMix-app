# Methylomics Differential Methylation Position (DMP) — Complete Audit and Teaching Documentation

**Source file (the entire submodule lives in one file):**
- `ArthOMix/R/methylomics/mod_methyl_dmp.R` (875 lines) — UI + server + every helper function for the "Differential Methylation (DMPs)" sub-module.

**Supporting files read to trace the full execution path (not modified, not documented in their own right except where they feed directly into DMP):**
- `ArthOMix/global.R` — `load_default_dmp()` (:313-320), `load_default_meth_pheno()` (:340-343), package `library()` calls (:57-122), `ARTHOMIX_COLORS`/`ARTHOMIX_STATUS`/`arthomix_pair()`/`theme_arthomix()` (:1417-1450ish), `METH_DATA_AVAILABLE`/`METH_DMP_PLAIN_DIR`/`METH_DMP_SVA_DIR` (defined in `data_paths.R`, sourced from `global.R`).
- `ArthOMix/R/methylomics/annotation.R` (102 lines) — `methyl_get_annotation()` (:48-93), `METHYL_ANNOTATION_PACKAGES` (:18-21).
- `ArthOMix/R/methylomics/qc.R` — `methyl_filter_missing()` (:30-34), `methyl_filter_variance()` (:36-41), `methyl_filter_snp()` (:68-81), `methyl_row_vars()` (:22-28), `methyl_sheet_sample_ids()` (:456-465).
- `ArthOMix/R/methylomics/normalization.R` — `methyl_norm_status()` (:424-455).
- `ArthOMix/R/submodules_registry.R` (:42-45) — registration of `mod_methyl_dmp_config`/`_ui`/`_server` into `MX_MODULES`.
- `ArthOMix/server.R` (:93-95) — module invocation: `mod_methyl_dmp_server("mx_dmp", methyl_dataset, methyl_results)`.
- `ArthOMix/data/preloaded/methylomics/tables/script03_dmp_sexstratified/METHODS_dmp_sexstratified.md` and `.../script03_dmp_sva_sexstratified/METHODS_dmp_sva_sexstratified.md` — the offline research pipeline's own methods write-up for the CSV tables the "SVA" sub-tab reads and re-displays. **These describe an offline analysis script, not code that runs inside the Shiny app** — see §3 and §10 for the distinction.
- `ArthOMix/DESCRIPTION` — declared package dependency surface, used in §24 (Reproducibility Audit).

Prepared: 2026-08-26.

This document is derived **exclusively** from the code cited above and is scoped strictly to Methylomics → Differential Methylation (DMPs). No other module, tab, file, or line of application code was modified in the course of this audit. Every non-trivial technical claim carries a `file:line` citation. Where a claim could not be verified from the inspected code, that is stated explicitly as **Not implemented in the current code** or **Not determinable from the inspected implementation**, never inferred. Two label conventions are used throughout, matching the house style of the companion `methylomics_quality_control.md` / `methylomics_normalization.md` documents in this same `thesis/` folder:

- **Scientific background** — a general statistical-genomics statement (textbook/literature knowledge), not a claim about this code.
- **Code evidence** — a statement about what `mod_methyl_dmp.R` (or a function it calls) actually does, always with a citation.

---

## 1. Module Overview

**Registration (code evidence).** `mod_methyl_dmp_config` declares `id = "dmp"`, `title = "Differential Methylation (DMPs)"`, `icon = "chart-scatter"`, `group = "Data"` (`mod_methyl_dmp.R:42-45`). It is the fourth entry in `MX_MODULES`, after Quality Control, Normalization, and Cell-Type Deconvolution (`submodules_registry.R:42-45`), and is invoked from `server.R` as `mod_methyl_dmp_server("mx_dmp", methyl_dataset, methyl_results)` inside the same `lapply(MX_MODULES, ...)` loop every other Methylomics sub-module is invoked from (`server.R:93-95`).

**Two independent analysis paths inside one module (code evidence, module's own header comment, `mod_methyl_dmp.R:1-40`).** The module implements two structurally different things, not two views of the same computation:

1. **"SVA" sub-tab — a reproduction/browser, not a live computation.** It reads four already-computed CSV result tables (`load_default_dmp("plain"|"sva", "female"|"male")`, `global.R:313-320`) that were produced by an **offline R script outside the Shiny app** (documented in `METHODS_dmp_sexstratified.md` / `METHODS_dmp_sva_sexstratified.md`), and lets the user filter/plot/download that already-computed table. It never calls `limma::lmFit()`, `sva::sva()`, or `bacon::bacon()` itself — those ran once, outside this app, to produce the CSV files this tab reads. This is a critical distinction for the audit and is treated as such throughout this document (see §10).
2. **"DMP" sub-tab — a genuine, configurable, live statistical engine.** It fits `limma::lmFit()` → `limma::eBayes()` → `limma::topTable()` on whatever beta/M-value matrix is currently loaded in `methyl_dataset$beta` (an upload, or — only in a deployment where the ~2.1GB raw matrix was copied in alongside the app — the preloaded dataset's own live matrix), against a user-chosen group/reference/comparison/covariate/sex configuration. **This path applies no SVA or bacon correction** — that is an explicit, code-documented limitation of this specific tab, not an oversight this audit is the first to notice (`mod_methyl_dmp.R:757-759` shows the app's own UI warning about this).

**What a "DMP" is (scientific background).** A Differentially Methylated Position is a single CpG dinucleotide whose methylation level (conventionally the *beta value*, the proportion of methylated alleles at that site, bounded 0–1) differs, on average, between two groups of samples by more than would be expected under the null hypothesis of no true difference, after accounting for multiple-testing burden across every CpG tested simultaneously. This module tests one CpG at a time (a "probe-wise" or "single-site" EWAS design), as distinct from a Differentially Methylated *Region* (DMR) analysis, which aggregates statistics across several spatially adjacent CpGs — DMR calling is implemented in the separate `mod_methyl_dmr.R` file and is explicitly out of scope for this document.

---

## 2. Scope of the Analysis

This document covers only:
- `mod_methyl_dmp_config`, `mod_methyl_dmp_ui`, `mod_methyl_dmp_server` and every function defined in `mod_methyl_dmp.R` (`mod_methyl_dmp_filter`, `mod_methyl_dmp_volcano`, `mod_methyl_lambda_gc`, `mod_methyl_qq_plot`, `methyl_chunked_lmfit`, `mod_methyl_dmp_sex_col`, `mod_methyl_dmp_sex_choices`, `mod_methyl_dmp_covariate_cols`, `mod_methyl_dmp_manhattan`, `mod_methyl_dmp_topplot`, `mod_methyl_dmp_betadist`).
- The external functions the above call, only as far as needed to explain what enters/leaves the DMP module (`load_default_dmp`, `load_default_meth_pheno`, `methyl_get_annotation`, `methyl_filter_missing`, `methyl_filter_variance`, `methyl_filter_snp`, `methyl_sheet_sample_ids`, `methyl_norm_status`).

Explicitly **not** covered, per the task's scope instructions: Methylomics Dataset (upload/parsing), Quality Control, Normalization, DMR, or any Transcriptomics/Proteomics/Cross-Omics/Multi-Omics module, even where the DMP module reads their shared reactive state (`methyl_dataset`) or is itself read by another module (`methyl_results$dmp`, consumed by `mod_cross_integration.R:244-250` for a one-line status display — noted in §17, not documented further, per scope).

---

## 3. Source Code Files

| File | Lines | Role in DMP |
|---|---|---|
| `ArthOMix/R/methylomics/mod_methyl_dmp.R` | 875 | Entire submodule: config, UI, all 11 helper functions, and the `moduleServer()` body. |
| `ArthOMix/global.R` | (partial) | Defines `load_default_dmp()`/`load_default_meth_pheno()` (precomputed-CSV readers for the SVA tab), package attachment (`library(limma)`, `library(sva)`), shared plot theme/palette. |
| `ArthOMix/R/methylomics/annotation.R` | 102 | `methyl_get_annotation()` — Illumina manifest lookup (chromosome/position/gene/SNP-overlap) used by the live DMP engine's annotation and SNP-filter steps. |
| `ArthOMix/R/methylomics/qc.R` | (partial) | `methyl_filter_missing()`, `methyl_filter_variance()`, `methyl_filter_snp()`, `methyl_sheet_sample_ids()` — probe/sample-matching utilities shared with the QC module, reused unmodified by the live DMP engine. |
| `ArthOMix/R/methylomics/normalization.R` | (partial) | `methyl_norm_status()` — read-only diagnostic message shown in the live DMP tab's header card; DMP does not call any normalization method itself. |
| `ArthOMix/R/submodules_registry.R` | :42-45 | Registers the module into `MX_MODULES`. |
| `ArthOMix/server.R` | :93-95 | Invokes `mod_methyl_dmp_server()`. |

No other file defines, overrides, or duplicates any function used by this module (verified by `grep -rn` for every function name above across `ArthOMix/R/`).

---

## 4. Package Dependencies

**Code evidence — what is actually attached at app start (`global.R:57-122`), relevant to DMP:**

| Package | `library()` call | Used by DMP for |
|---|---|---|
| `limma` | `global.R:74` | `stats::model.matrix()` (base R, not limma), `limma::lmFit()` (wrapped by `methyl_chunked_lmfit()`), `limma::makeContrasts()`, `limma::contrasts.fit()`, `limma::eBayes()`, `limma::topTable()` — the entire live-engine statistical core. |
| `sva` | `global.R:82` | **Not called anywhere in `mod_methyl_dmp.R`** (verified: zero occurrences of `sva::` or a bare `sva(` call in the file). It is attached at app start for other modules; the "SVA" sub-tab's name refers to the *offline* pipeline that produced the CSVs it reads (§10), not to a live call from this file. |
| `data.table` | `global.R:69` | Used indirectly — `load_default_dmp()`/`load_default_meth_pheno()` call `data.table::fread()` (`global.R:319,342`). |
| `ggplot2`, `ggrepel` (unused here) | `global.R:64-65` | All six plotting helpers (`mod_methyl_dmp_volcano`, `mod_methyl_qq_plot`, `mod_methyl_dmp_manhattan`, `mod_methyl_dmp_topplot`, `mod_methyl_dmp_betadist`) build `ggplot2` objects. `ggrepel` is attached app-wide but **not used** by any DMP plot (no `geom_text_repel`/`geom_label_repel` call in `mod_methyl_dmp.R`). |
| `DT` | `global.R:63` | Both results tables (`default_table`, `live_table`) use `DT::renderDataTable()`/`DT::datatable()`/`DT::formatSignif()`. |
| `shiny`, `shinycssloaders` | `global.R:57,62` | UI framework; `withSpinner()` wraps every `uiOutput`/`plotOutput` that can take noticeable time. |
| `matrixStats` (optional, via `requireNamespace`) | not `library()`-attached; checked at call time | `methyl_row_vars()` (`qc.R:22-28`, called by `methyl_filter_variance()`) uses `matrixStats::rowVars()` when installed, falling back to base `apply()` otherwise. |

**`bacon` is not a declared dependency anywhere (code evidence, not an oversight of this audit).** `bacon::bacon()` is named in `METHODS_dmp_sva_sexstratified.md` (the offline pipeline's own methods text) but does **not** appear in `ArthOMix/DESCRIPTION`'s `Imports:` list, is not `library()`-attached in `global.R`, and is never called from `mod_methyl_dmp.R`. This is consistent with — not contradicting — the file's own architecture: bacon-correction happened once, offline, to produce the `fdr_bacon`/`p_bacon` columns already present in the CSVs the "SVA" tab reads (`dmp_female_full.csv` header: `cpg,logFC_M,dbeta,t,p_raw,p_bacon,fdr_bacon`); the app itself never needs the `bacon` package installed to display those columns. The live "DMP" engine, which genuinely does run inside the app, applies no bias/inflation correction at all (§10, §13).

**Full app-wide `Imports:` surface** is documented in `ArthOMix/DESCRIPTION`; only the packages named above are exercised by code paths this document traces.

---

## 5. DMP Submodule Tab Structure

**Total number of tabs/sub-tabs: 2** — a `tabsetPanel` with exactly two `tabPanel`s, `"SVA"` and `"DMP"` (`mod_methyl_dmp.R:57-67`).

```r
tabsetPanel(
  id = ns("dmp_subtabs"), type = "tabs",
  tabPanel("SVA", br(), withSpinner(uiOutput(ns("default_ui")), color = "#2563EB", type = 6)),
  tabPanel("DMP", br(), withSpinner(uiOutput(ns("live_ui")), color = "#2563EB", type = 6))
)
```

**Code evidence — this is purely a UI reorganization, per the module's own comment (`mod_methyl_dmp.R:47-56`):** the two tabs wrap the exact same two `uiOutput`s (`default_ui`, `live_ui`) the module already had stacked vertically before this tab split; no output ID, reactive, or server logic differs between "having two tabs" and "having one long page." Because Shiny's `tabsetPanel` only suspends a hidden tab's *rendering*, not its underlying `eventReactive` results, switching between "SVA" and "DMP" never forces either tab's analysis to recompute (`mod_methyl_dmp.R:53-56`).

| # | Tab (exact label) | Gated on | Renders |
|---|---|---|---|
| 1 | **SVA** | `methyl_dataset$preloaded` is `TRUE` (`req()` at top of `default_ui`, `mod_methyl_dmp.R:283`) | Reproduction/browser of four precomputed CSVs (plain + SVA-adjusted, female + male). |
| 2 | **DMP** | `methyl_dataset$beta` is not `NULL` (`mod_methyl_dmp.R:439`) | Configurable live `limma` engine on whatever matrix is currently loaded. |

Both can be visible together (a preloaded dataset whose own live matrix is also available in this deployment); either can be visible alone; **if neither condition holds** (no dataset loaded at all), "SVA" renders nothing (not even a placeholder — `req()` silently blocks, `mod_methyl_dmp.R:283`) and "DMP" renders an explanatory upload prompt (`mod_methyl_dmp.R:439-448`).

Below, §6.1 documents the **SVA** tab and §6.2 documents the **DMP** tab, each against the full 20-point checklist requested.

---

## 6. Tab-by-Tab Analysis

### 6.1 Tab 1 — SVA

#### Purpose

Reproduces, filters, plots, and lets the user download the app's own bundled sex-stratified DMP analysis of the preloaded GSE42861 whole-blood RA-vs-Control dataset — specifically the surrogate-variable-adjusted, bacon-corrected pipeline stage (`METHODS_dmp_sva_sexstratified.md`), which the module's own explanatory card states is "the panel actually used downstream" (`mod_methyl_dmp.R:315-316`).

#### User Input

| Control | Type | Default | Required/Optional |
|---|---|---|---|
| `sva_sex` (Stratum) | `radioButtons`, inline, `Female`/`Male` | `"female"` | Required |
| `sva_fdr` (FDR threshold) | `numericInput`, 0–1, step 0.01 | `0.05` | Required |
| `sva_dbeta` (Min \|Δβ\|) | `numericInput`, 0–1, step 0.01 | `0` | Required |
| `sva_direction` | `radioButtons`, inline, `Any`/`Hypermethylated`/`Hypomethylated` | `"any"` | Required |
| `sva_run_btn` ("View results") | `actionButton` | — | Must be clicked before any result renders |

Source: `mod_methyl_dmp.R:298-305`. There is **no** group/reference/comparison selector, no covariate selector, and no dataset selector on this tab — the underlying analysis (which comparison, which covariates, which model) is entirely fixed by the offline pipeline; only *post-hoc filtering* of its already-computed output is configurable here.

#### Input Data Structure

`default_data()` (`mod_methyl_dmp.R:269-276`) assembles a list of five already-computed, flat CSV-derived `data.frame`s, none of which are matrices of probes × samples — the DMP tables are **one row per CpG, already collapsed across samples**:

- `pheno` — `load_default_meth_pheno()`, one row per sample (689 GSE42861 samples), columns include `group` (`"RA"`/`"Control"`) and `sex` (`"F"`/`"M"`) (`mod_methyl_dmp.R:272`, used only to compute the descriptive counts shown in the header, `mod_methyl_dmp.R:286-287`).
- `plain_f`, `plain_m`, `sva_f`, `sva_m` — one `data.frame` per {pipeline stage} × {sex}, each one row per CpG. Verified column header of `dmp_female_full.csv` (SVA stage): `cpg,logFC_M,dbeta,t,p_raw,p_bacon,fdr_bacon`.

#### Processing

1. `default_data()` loads all five tables unconditionally whenever `METH_DATA_AVAILABLE` is `TRUE` (`mod_methyl_dmp.R:270-276`); this is a `reactive()`, so it is cached and only re-evaluated if its (non-existent) reactive dependencies change — in practice it runs once per session.
2. **`plain_f`/`plain_m` are loaded but only ever used to compute two summary counts** (`n_sig_plain_f`, `n_sig_plain_m`, `mod_methyl_dmp.R:290-291`) shown in static prose (`mod_methyl_dmp.R:312-317`). The plain-stage *table itself* is never filterable, plottable, or downloadable from this UI — only the SVA-stage table is (see §22, Finding INFO-1).
3. On clicking "View results", `sva_run` (an `eventReactive`, `mod_methyl_dmp.R:328-337`) selects `d$sva_f` or `d$sva_m` per the chosen stratum and packages the chosen FDR/Δβ/direction filter values.
4. `default_filtered` (`mod_methyl_dmp.R:342-345`) applies `mod_methyl_dmp_filter()` (§9) to that table.
5. No sample-matching, no group assignment, no design matrix, no model fitting happens on this tab — all of that already happened in the offline pipeline that produced the CSV.

#### Functions Used

See §9 for full documentation of `load_default_dmp`, `mod_methyl_dmp_filter`, `mod_methyl_dmp_volcano`.

#### Statistical Analysis

**Not run live.** The statistics displayed (`t`, `p_raw`, `p_bacon`, `fdr_bacon`, `dbeta`, `logFC_M`) are read verbatim from CSV. Per `METHODS_dmp_sva_sexstratified.md`, the offline pipeline that produced them: fit `limma::lmFit()` on M-values with design `~ group + age + smoking + <6 cell-type fractions> + <SVA surrogate variables>` per sex stratum, moderated with `limma::eBayes()`, passed the `group` coefficient's moderated t-statistic to `bacon::bacon()` for bias/inflation-corrected p-values, then applied Benjamini-Hochberg FDR within each stratum independently. **This document does not re-verify that offline pipeline's own correctness** — it is out of scope (it is not part of the Shiny application's own code); it is described here only so the columns this tab displays are traceable to a documented method rather than invented.

#### Output

Filtered `data.frame` with columns `cpg, logFC_M, dbeta, t, p_raw, p_bacon, fdr_bacon` (`mod_methyl_dmp.R:381-384`), no `chr`/`pos`/`gene` annotation columns (this stage's precomputed table does not carry them, unlike the live DMP tab's output — see §14).

#### Figures

One figure: **volcano plot** (`output$default_volcano`, `mod_methyl_dmp.R:354-357`), via the shared `mod_methyl_dmp_volcano()` helper (§9), x = `dbeta`, y = `-log10(fdr_bacon)`, colored by whether a point passes the currently-selected FDR/Δβ thresholds. Rendered only after "View results" is clicked (`mod_methyl_dmp.R:347-352`) — **generated after statistical filtering is defined but drawn over the full (unfiltered) table**, with significance encoded by point color rather than by removing non-significant points.

#### Downloads

One `downloadHandler`, `download_default` (`mod_methyl_dmp.R:399-402`): writes `default_filtered()` (the FDR/Δβ/direction-filtered subset only, never the full table) to CSV, filename `gse42861_dmp_sva_<sex>_fdr<threshold>.csv`.

---

### 6.2 Tab 2 — DMP

#### Purpose

A fully configurable, live differential methylation model fit against whatever beta/M-value matrix + sample sheet is currently loaded on the Dataset tab — an upload, or (only when this deployment has the raw matrix copied in) the preloaded dataset's own live matrix. Unlike the SVA tab, nothing about this tab's group/reference/comparison/covariate configuration is hardcoded to GSE42861 (`mod_methyl_dmp.R:29-36`).

#### User Input

| Control | Type | Default | Required/Optional |
|---|---|---|---|
| `live_sex` | `radioButtons`, dynamic choices from the sample sheet's sex column, plus a fixed `"All samples"` option | `"__all__"` | Required (always has a value) |
| `live_group_col` | `selectInput`, choices = every sample-sheet column | first of `group`/`Group`/`disease`/`Disease` found, else first column | Required |
| `live_ref` / `live_comp` (Reference/Comparison group) | `selectInput` × 2, choices = distinct values of the chosen group column | 1st / 2nd distinct level | Required, and must differ |
| `live_fdr` (FDR threshold) | `numericInput`, 0–1, step 0.01 | `0.05` | Required |
| `live_dbeta` (Absolute Δβ threshold) | `numericInput`, 0–1, step 0.01 | `0` | Required |
| `live_direction` | `radioButtons`, `All DMPs`/`Hypermethylated`/`Hypomethylated` | `"any"` | Required |
| `live_min_valid_pct` (Min valid sample %) | `numericInput`, 0–100, step 5 | `80` | Required |
| `live_min_variance` (Min methylation variance) | `numericInput`, ≥0 | `0` | Optional (0 = no-op filter) |
| `live_snp_filter` | `checkboxInput` — only shown when array annotation is available | `FALSE` | Optional |
| `live_covariates` | `checkboxGroupInput`, dynamic choices from remaining sample-sheet columns | none selected | Optional |
| `live_run_btn` ("Run DMP Analysis") | `actionButton` | — | Must be clicked; nothing computes automatically |

Source: `mod_methyl_dmp.R:468-498`. Sex is excluded from the covariate list whenever a single-sex subset is already selected, since it would then be constant within the model (`mod_methyl_dmp.R:519`, comment `:516-518`).

#### Input Data Structure

- `methyl_dataset$beta` — probe × sample matrix. **Rows = CpG probes** (rownames are Illumina probe IDs, e.g. `cg#######`), **columns = samples** (colnames are sample IDs) — confirmed by every downstream operation indexing rows with probe filters (`keep_probe`, `mod_methyl_dmp.R:620-635`) and columns with sample subsetting (`beta0[, keep_sex, ...]`, `mod_methyl_dmp.R:562`).
- `methyl_dataset$sample_sheet` — one row per sample, arbitrary phenotype columns; matched to matrix columns via `methyl_sheet_sample_ids()` (§9), not assumed to be pre-aligned or same-order (`mod_methyl_dmp.R:541-543`).
- `methyl_dataset$input_scale` — `"beta"` or `"m"`, determines whether the loaded matrix is on the 0–1 beta scale or the logit M-value scale (`mod_methyl_dmp.R:615`).
- `methyl_dataset$array_type` — feeds `methyl_get_annotation()` for chromosome/position/gene/SNP lookup (`mod_methyl_dmp.R:421-424`).

#### Processing

The full 22-step pipeline inside `live_result` (`mod_methyl_dmp.R:533-705`), triggered only by `input$live_run_btn`:

1. **Validation gate** (`mod_methyl_dmp.R:534-539`): matrix loaded, sheet loaded, group column chosen and valid, reference/comparison chosen and different — any failure shows a `validate(need(...))` message in place of results, no partial computation happens.
2. **Sample matching** (`mod_methyl_dmp.R:541-543`): `methyl_sheet_sample_ids()` resolves sheet rows to matrix column IDs; `intersect()` with the matrix's own column names; requires ≥6 common samples.
3. **Sheet coercion** (`mod_methyl_dmp.R:544-552`): the sheet is coerced to a plain `data.frame` — the preloaded dataset's own sheet is a `data.table`, whose single-bracket column selection by a character vector means something different than for a plain `data.frame`; coercing once avoids that ambiguity for every downstream `ph0`/`ph1` access.
4. **Sex subset** (`mod_methyl_dmp.R:555-566`), only if `live_sex != "__all__"`: requires a detected sex column; requires ≥6 samples remain after restricting.
5. **Group/reference/comparison subset** (`mod_methyl_dmp.R:569-573`): keeps only rows whose group value is exactly the reference or comparison level; the surviving `grp` is a two-level `factor(levels = c(ref, comp))` — reference level ordered first.
6. **Memory cleanup** (`mod_methyl_dmp.R:574-580`): `rm(beta0)` — explicitly documented as necessary because a full-size intermediate copy of a genome-wide (~400k-probe) matrix is otherwise held alive for no reason, one of several redundant copies the file's own comments say previously exceeded R's vector memory limit on the real preloaded matrix.
7. **Optional covariates, complete-cases only** (`mod_methyl_dmp.R:582-594`): if any covariate columns are selected, rows with any missing value in *any* selected covariate are dropped entirely (both from the matrix and the phenotype); requires ≥6 complete-case samples remain; every non-numeric covariate column is coerced to a `factor`.
8. **Minimum group size check** (`mod_methyl_dmp.R:596-600`): both reference and comparison groups need ≥3 samples after all the above subsetting.
9. **Probe filters, on the beta scale, before the M-value transform** (`mod_methyl_dmp.R:602-637`): if the loaded matrix is on the M-value scale, a full-matrix-sized beta-scale copy is computed once (`beta_scale_full`) purely to run the filters, then discarded; if the matrix is already beta-scale, no extra copy is made. Three filters, each independently reported:
   - Missingness (`methyl_filter_missing()`, threshold = `1 - live_min_valid_pct/100`).
   - Variance (`methyl_filter_variance()`, threshold = `live_min_variance`, `0` by default = no probes removed).
   - SNP overlap (`methyl_filter_snp()`), only if the checkbox is ticked **and** manifest annotation is available for the array type.
   Requires ≥10 probes survive combined. The code's own comment (`mod_methyl_dmp.R:602-614`) states this ordering (filter-before-transform) was deliberately chosen over the reverse to avoid holding 3+ full-size matrix copies simultaneously, which previously exceeded R's 16GB vector memory limit on the real 412,492-probe × 689-sample preloaded matrix.
10. **M-value / beta-scale derivation on the filtered subset only** (`mod_methyl_dmp.R:639-641`): if the input was beta-scale, M = `log2(beta/(1-beta))` with beta clipped to `[1e-6, 1-1e-6]` first (undefined at exact 0/1); if the input was already M-value scale, an inverse-logit beta-scale matrix is derived instead, for the descriptive Δβ/plot use below — the model itself always fits on `m` (whichever scale that variable holds after this step, i.e. always the M-value scale, since if the input was already M it is used directly and if it was beta it was just converted).
11. **Design matrix** (`mod_methyl_dmp.R:644-658`): `design_grp <- model.matrix(~0 + grp)` — a no-intercept, means-parameterized two-column group design; if covariates are selected, `model.matrix()` is built from a backtick-quoted formula over the covariate columns, its own intercept column dropped, and `cbind()`-appended to `design_grp`. The combined design's rank is checked (`qr(design)$rank == ncol(design)`); a rank-deficient design (e.g. a covariate collinear with group) blocks the run with an explanatory message rather than silently fitting a degenerate model.
12. **Chunked model fit** (`mod_methyl_dmp.R:660`): `methyl_chunked_lmfit(m, design)` (§9) — row-chunked `limma::lmFit()`, verified (per the function's own header comment) to be bit-for-bit identical to a whole-matrix fit.
13. **Contrast + moderation** (`mod_methyl_dmp.R:661-664`): `limma::makeContrasts(comp - ref, levels = design)`, `limma::contrasts.fit()`, `limma::eBayes()`. Both steps are wrapped in `tryCatch()` with a `validate(need(...))` fallback message rather than an uncaught R error.
14. **`topTable()`** (`mod_methyl_dmp.R:665`): `limma::topTable(fit2, number = Inf, sort.by = "P")` — every tested probe returned, sorted by raw p-value.
15. **Descriptive Δβ** (`mod_methyl_dmp.R:667-669`): `dbeta` = row mean beta (comparison group) minus row mean beta (reference group), computed independently of the M-value model coefficient, on `beta_scale` (always the 0–1 scale regardless of which scale the model itself fit on).
16. **Result assembly** (`mod_methyl_dmp.R:671-675`): one `data.frame`, columns `cpg, t, p_raw, fdr, dbeta, ref_mean_beta, comp_mean_beta`.
17. **Annotation** (`mod_methyl_dmp.R:676-684`): if `methyl_get_annotation()` succeeded for the loaded array type, `chr`/`pos`/`gene` columns are filled in by row-matching against the manifest; otherwise they remain `NA` for every row (never fabricated).
18. **Direction label** (`mod_methyl_dmp.R:685`): `"hyper"` if `dbeta > 0`, else `"hypo"` — purely descriptive, computed for every tested probe regardless of significance.
19. **Design-formula string** (`mod_methyl_dmp.R:687-688`): a human-readable `"Methylation ~ <group_col> [+ covariates]"` string built for display only, not used to refit anything.
20. **Genomic inflation diagnostic** (`mod_methyl_dmp.R:702`): `mod_methyl_lambda_gc(df$p_raw)` computed on the full, unfiltered raw p-value vector.
21. Result `list()` returned (`mod_methyl_dmp.R:690-704`) carries the table plus every configuration/sample-size/filter value needed by the results UI — nothing is recomputed later from `input$...` directly except the *post-hoc* FDR/Δβ/direction display filter (§6.2 "Statistical Analysis" below).
22. `observeEvent(live_result(), ...)` (`mod_methyl_dmp.R:707-713`) publishes a three-field summary (`comparison`, `n_probes`, `n_sig`) to the shared `methyl_results$dmp` object and flips `live_has_run` to `TRUE`.

#### Functions Used

`methyl_sheet_sample_ids`, `methyl_filter_missing`, `methyl_filter_variance`, `methyl_filter_snp`, `methyl_chunked_lmfit` (wrapping `limma::lmFit`), `stats::model.matrix`, `qr`, `limma::makeContrasts`, `limma::contrasts.fit`, `limma::eBayes`, `limma::topTable`, `methyl_get_annotation`, `mod_methyl_lambda_gc` — all fully documented in §9.

#### Statistical Analysis

**Design matrix:** means-parameterized, `~0 + group [+ covariates]` — two named coefficient columns for the group levels (no shared intercept absorbed into a baseline), plus one column per covariate level/unit. **Contrast:** `comparison_group − reference_group`, i.e. a simple two-group difference-of-means test on the M-value scale, moderated. **Test:** `limma`'s moderated t-statistic — an ordinary per-probe t-test whose standard error is empirical-Bayes-shrunk toward a common prior across all tested probes (`eBayes()`), which increases power and stability especially with small per-group sample sizes, at the cost of assuming that error variances are exchangeable across probes (a standard, well-established assumption for array-based EWAS, not unique to this app). **Multiple testing:** `topTable()`'s `adj.P.Val` column, Benjamini-Hochberg FDR by default (limma's own default `topTable(..., adjust.method = "BH")` is not overridden anywhere in this call — verified: no `adjust.method` argument is passed at `mod_methyl_dmp.R:665`, so limma's own default `"BH"` applies). **No bias/inflation correction** (no bacon, no genomic-control rescaling of p-values) is applied — `mod_methyl_lambda_gc()` is computed and displayed as a *diagnostic*, but the p-values/FDR values shown are the raw `eBayes()`/`topTable()` output, unmodified by that diagnostic (§13).

**Post-run, non-refitting filters:** the FDR threshold, Δβ threshold, and direction filter remain live-adjustable *after* a run completes without re-fitting the model (`live_filtered`, `mod_methyl_dmp.R:718-722`) — only changing sex/group/covariate/probe-filter controls and re-clicking "Run DMP Analysis" changes what was actually fit.

#### Output

`live_result()$df` columns: `cpg, t, p_raw, fdr, dbeta, ref_mean_beta, comp_mean_beta, chr, pos, gene, direction` (`mod_methyl_dmp.R:671-685`). The results table adds a computed `significant` (`"Yes"`/`"No"`) column at render time only (`mod_methyl_dmp.R:842`), not stored in the underlying data.

#### Figures

Six, all only after "Run DMP Analysis" is clicked and only for the currently-loaded, currently-fit result:

| Figure | Function | x | y | Notes |
|---|---|---|---|---|
| QQ plot | `mod_methyl_qq_plot` | expected `-log10(p)` | observed `-log10(p)` | Raw p-values, unfiltered — genomic-inflation diagnostic. |
| Volcano | `mod_methyl_dmp_volcano` | Δβ | `-log10(fdr)` | Colored by pass/fail of current FDR+Δβ thresholds; drawn over the full table, not the filtered subset. |
| Manhattan | `mod_methyl_dmp_manhattan` | genome position, ordered by chr | `-log10(fdr)` | Only shown if annotation succeeded; drops any CpG missing chr/pos. |
| Top DMPs bar chart | `mod_methyl_dmp_topplot` | Δβ | CpG (± gene) label | Ranked by FDR, \|Δβ\|, or FDR/\|Δβ\| combined; top N ∈ {10,20,50,100}, user-selectable. |
| β-value distribution | `mod_methyl_dmp_betadist` | CpG | β value | Boxplots, one box per group, restricted to the CpGs shown in the top-DMP plot. |

#### Downloads

Two `downloadHandler`s (`mod_methyl_dmp.R:849-873`):
1. `download_live` — `live_filtered()` (FDR/Δβ/direction-filtered subset) to CSV, `dmp_<comp>_vs_<ref>_<sex>.csv`.
2. `download_live_config` — a 19-row parameter/value CSV documenting the run's exact configuration (dataset, sex, groups, sample sizes, covariates, method string, design formula, all four filter values, probe counts, significant count, timestamp) — `dmp_analysis_configuration.csv`. This is the module's reproducibility export.

---

## 7. Input Data

**Scientific background.** A methylation-array or bisulfite-sequencing experiment produces, per CpG site per sample, a *beta value* — the proportion of methylated to total (methylated + unmethylated) signal at that site, mathematically bounded to [0, 1] (0 = fully unmethylated, 1 = fully methylated). This is the field's standard primary measurement unit.

**Code evidence.** The DMP module never reads a raw file itself; it consumes `methyl_dataset$beta` and `methyl_dataset$sample_sheet`, populated entirely by the (out-of-scope) Dataset tab. `methyl_dataset$input_scale` (`"beta"` or `"m"`) tells the DMP engine which scale that matrix is already on (`mod_methyl_dmp.R:615`) — the module does not independently detect the scale from the data's numeric range.

## 8. Data Structures

- **`methyl_dataset$beta`**: numeric matrix, rows = CpG probes (Illumina `cg#######`-style IDs), columns = samples. Confirmed by row-indexed probe filters (`keep_probe`) and column-indexed sample subsetting throughout `live_result` (`mod_methyl_dmp.R:562,571,589,635`).
- **`methyl_dataset$sample_sheet`**: one row per sample, arbitrary columns; matched to matrix columns by `methyl_sheet_sample_ids()`, not by row order alone.
- **SVA-tab CSVs**: one row per CpG, no sample dimension at all (already-summarized statistics only) — columns `cpg, logFC_M, dbeta, t, p_raw, p_bacon, fdr_bacon`.
- **Live-tab result `df`**: one row per tested CpG, columns `cpg, t, p_raw, fdr, dbeta, ref_mean_beta, comp_mean_beta, chr, pos, gene, direction`.

---

## 9. Function-by-Function Documentation

Functions defined inside `mod_methyl_dmp.R` receive the full requested template. External functions this module calls are documented in a compact reference table (§9.13) — per the task's own instruction not to blindly document every package/base-R function, only those genuinely on this module's execution path.

### `mod_methyl_dmp_filter()`

**Package:** none (project function). **Location:** `mod_methyl_dmp.R:71-76`.
**Purpose:** Applies the FDR / minimum-\|effect\| / direction filter shared by both sub-tabs, parameterized by which column names to test.
**Input:** `df` (a DMP results table), `fdr_col`/`effect_col` (column names to test), `fdr_max`, `effect_min`, `direction` (`NULL`/`"hyper"`/`"hypo"`).
**Operation:** Builds a logical `keep` vector: FDR column non-`NA` and ≤ `fdr_max`, AND `abs(effect_col) >= effect_min`; if `direction` is `"hyper"`, additionally requires `effect_col > 0`; if `"hypo"`, requires `effect_col < 0`.
**Output:** `df[keep, ]` — same columns, fewer rows.
**Role in DMP Pipeline:** The single filtering implementation reused for the SVA tab's table/volcano and the DMP tab's table/volcano/value-boxes — one definition, so both tabs' filter semantics cannot silently diverge.
**Audit:** No implementation issue identified from the inspected code. Note: rows where `fdr_col` is `NA` are always excluded, never treated as "unknown/pass" — an appropriate default for a significance filter.

### `mod_methyl_dmp_volcano()`

**Package:** `ggplot2` (project wrapper). **Location:** `mod_methyl_dmp.R:78-87`.
**Purpose:** Builds the shared volcano plot (effect size vs. significance) used by both sub-tabs.
**Input:** `df`, `effect_col`, `p_col` (the column plotted on -log10, either `fdr_bacon` on the SVA tab or `fdr` on the DMP tab — see Audit below), `effect_label`, `fdr_max`, `effect_min`.
**Operation:** Computes `neglog10p = -log10(pmax(p_col, .Machine$double.xmin))` (floors at machine epsilon to avoid `-log10(0) = Inf`); computes a `sig` logical using the *same* pass/fail rule as `mod_methyl_dmp_filter()` but re-implemented inline rather than calling it; plots `geom_point()` colored by `sig`.
**Output:** a `ggplot` object.
**Role in DMP Pipeline:** Primary at-a-glance significance/effect-size visualization for both tabs.
**Audit:** **Naming note, not a defect.** The function's `p_col` parameter is always fed an *adjusted* p-value (FDR) at both call sites (`mod_methyl_dmp.R:356,818`), never a raw p-value, so the y-axis is genuinely `-log10(FDR)` in both uses despite the generic parameter name — verified correct at both call sites, not misleading in practice. The significance rule is duplicated (here and in `mod_methyl_dmp_filter()`) rather than the plot calling the filter function directly; both implementations were compared line-by-line and are logically identical, so this is a minor maintainability observation (duplicated logic, two places to keep in sync), not a correctness issue.

### `mod_methyl_lambda_gc()`

**Package:** `stats` (project wrapper). **Location:** `mod_methyl_dmp.R:97-101`.
**Purpose:** Computes the genomic inflation factor λ, the standard median-chi-square diagnostic for whether a genome-wide test's p-value distribution is systematically inflated relative to the null.
**Input:** `p` — a numeric vector of raw p-values, required by the function's own header comment (`mod_methyl_dmp.R:89-96`) to be the **full, unfiltered** per-CpG vector from the just-fitted model, not a post-filtering subset (a subset would bias the estimate).
**Operation:** `stats::median(stats::qchisq(1 - p, df = 1)) / stats::qchisq(0.5, df = 1)` — converts each p-value to its equivalent 1-df chi-square statistic, takes the median across all tested probes, and divides by the chi-square distribution's own theoretical median. A value of 1.0 indicates a perfectly calibrated null; values above ~1.1 are conventionally read as meaningful inflation.
**Output:** a single numeric λ, or `NA_real_` if no finite p ∈ (0,1] is present.
**Role in DMP Pipeline:** The live DMP engine's only calibration diagnostic, since (unlike the offline SVA pipeline) it applies no bias/inflation correction of its own — this is how a user of the live engine is told the p-values may be optimistic.
**Audit:** No implementation issue identified from the inspected code; the formula matches the standard genomic-control definition. Correctly documented (in the function's own header) as needing the *pre-filter* p-vector, and correctly called that way at its one call site (`mod_methyl_dmp.R:702`, on `df$p_raw` before any FDR/Δβ filtering).

### `mod_methyl_qq_plot()`

**Package:** `ggplot2` (project wrapper). **Location:** `mod_methyl_dmp.R:106-116`.
**Purpose:** The standard observed-vs-expected QQ companion plot to λ.
**Input:** `p` — raw p-value vector (same full, unfiltered vector as `mod_methyl_lambda_gc`).
**Operation:** Sorts finite p ∈ (0,1], builds `observed = -log10(p)` against `expected = -log10(stats::ppoints(n))` (the theoretical uniform-order-statistic expectation), plots points plus a dashed y = x reference line.
**Output:** a `ggplot` object, or `NULL` if no valid p-values.
**Role in DMP Pipeline:** Visual companion to λ; systematic departure above the diagonal indicates inflation beyond what λ alone summarizes as one number.
**Audit:** No implementation issue identified from the inspected code.

### `methyl_chunked_lmfit()`

**Package:** `limma` (project wrapper). **Location:** `mod_methyl_dmp.R:143-157` (rationale comment `:120-142`).
**Purpose:** A memory-safe, row-chunked replacement for a direct `limma::lmFit(m, design)` call on a full genome-wide matrix.
**Input:** `m` (probe × sample matrix — M-values in this module's only call site), `design` (the model matrix), `chunk_size` (default `20000` rows).
**Operation:** If `nrow(m) <= chunk_size`, delegates straight to `limma::lmFit()`. Otherwise splits rows into chunks, fits each chunk independently via `limma::lmFit()`, then reassembles one combined `MArrayLM` object by concatenating each fit's per-row 1-D fields (`df.residual`, `sigma`, `Amean`) with `c()` and per-row 2-D fields (`coefficients`, `stdev.unscaled`) with `rbind()`.
**Output:** an `MArrayLM` object, structurally identical to what `limma::lmFit()` would have returned on the whole matrix at once.
**Role in DMP Pipeline:** Makes the live engine usable on a real genome-wide array (400k+ probes); the function's own comment states a direct whole-matrix `lmFit()` reproduced a real "vector memory limit of 16.0 Gb reached" crash against the preloaded dataset's actual 412,492-probe × 689-sample matrix on a 16GB machine.
**Audit:** No implementation issue identified from the inspected code. The correctness argument (chunked-vs-whole-matrix `lmFit()` fits are combinable because each row's fit is independent given a shared design, and `eBayes()`'s moderation runs once, after combining, over the real per-row summary statistics either way) is sound and is stated in the code's own comment as having been verified bit-for-bit identical on synthetic data before adoption (`mod_methyl_dmp.R:140-142`) — this document did not independently re-run that verification (out of scope for a read-only audit) but the reasoning is standard and the claim is specific and falsifiable, not vague.

### `mod_methyl_dmp_sex_col()`

**Package:** none. **Location:** `mod_methyl_dmp.R:163-168`.
**Purpose:** Detects which sample-sheet column (if any) encodes sex/gender, by name only.
**Input:** `sheet` (sample sheet data.frame).
**Operation:** `intersect(c("sex","Sex","SEX","gender","Gender"), colnames(sheet))`, returns the first match or `NULL`.
**Output:** a single column name, or `NULL`.
**Role in DMP Pipeline:** Drives whether the Sex radio control offers Female/Male options at all, and whether sex can be offered as a covariate.
**Audit:** No implementation issue identified from the inspected code. The candidate list is intentionally the same one already used by `mod_methyl_qc.R`'s sex-check panel (per the function's own comment, `mod_methyl_dmp.R:159-162`) so "which column is sex" is one consistent guess across the Methylomics module — a reasonable design choice, not verified by this document beyond confirming the candidate list is in fact identical (out of scope: `mod_methyl_qc.R` itself).

### `mod_methyl_dmp_sex_choices()`

**Package:** none. **Location:** `mod_methyl_dmp.R:176-187`.
**Purpose:** Builds the "All samples / Female / Male" (or fallback) radio-button choice list from whatever the detected sex column's *actual values* are, never from a hardcoded assumption.
**Input:** `sheet`, `sex_col` (from the function above).
**Operation:** Always includes `"All samples" = "__all__"`. If no sex column, returns only that. Otherwise takes the distinct non-missing values; if exactly one level case-insensitively matches `^f(emale)?$` and exactly one matches `^m(ale)?$`, labels them `"Female"`/`"Male"`; otherwise falls back to using the raw values themselves as both label and value (never mislabeling an unfamiliar encoding).
**Output:** a named character vector suitable for `radioButtons(choices = ...)`.
**Role in DMP Pipeline:** Ensures the sex control works for any dataset's own sex-coding convention, not just GSE42861's `"F"`/`"M"`.
**Audit:** No implementation issue identified from the inspected code. Correctly generalizes beyond a hardcoded two-level assumption (e.g. a sheet with `"Female"`/`"Male"`/`"Unknown"` — three levels — falls into the raw-values fallback rather than mis-mapping "Unknown" onto Male or Female).

### `mod_methyl_dmp_covariate_cols()`

**Package:** none. **Location:** `mod_methyl_dmp.R:193-203`.
**Purpose:** Determines which sample-sheet columns are eligible to offer as covariates.
**Input:** `sheet`, `exclude` (ID columns, the chosen group column, and — conditionally — the sex column).
**Operation:** For every non-excluded column: excluded if it has fewer than 2 distinct non-missing values (constant, uninformative); excluded if non-numeric **and** every value is unique (free-text/ID-like); a numeric column is always kept even if every value happens to be unique (continuous covariates like age or a cell-type fraction are expected to vary continuously).
**Output:** a character vector of eligible column names.
**Role in DMP Pipeline:** Populates the "Covariates (optional)" checkbox group.
**Audit:** No implementation issue identified from the inspected code.

### `mod_methyl_dmp_manhattan()`

**Package:** `ggplot2` (project wrapper). **Location:** `mod_methyl_dmp.R:205-223`.
**Purpose:** Genome-wide Manhattan-style significance plot.
**Input:** `df` (must carry `chr`, `pos`, `fdr` columns), `fdr_max`.
**Operation:** Drops any row missing `chr`/`pos`/`fdr`; coerces `chr` to an integer by stripping a leading `"chr"` prefix (drops any row where that doesn't parse to an integer — i.e. this plot silently omits chrX/chrY/chrM/scaffold-named probes, since `"X"`, `"Y"`, `"M"` do not parse as integers via `as.integer()`); `validate(need(...))` blocks with an explanatory message if zero rows survive; sorts by chromosome then position, assigns a sequential plot index, colors alternating chromosomes, overplots FDR-significant points in a highlight color, adds a dashed significance threshold line.
**Output:** a `ggplot` object (or a `validate()` block).
**Role in DMP Pipeline:** The DMP tab's genome-wide overview figure; only shown when annotation succeeded.
**Audit:** **Informational — silent scope narrowing, not a bug.** Because `chr_num <- as.integer(gsub("^chr", ...))` and `is.na(d$chr_num)` rows are dropped without a user-facing count of how many were dropped for that specific reason (as opposed to missing chr/pos/fdr entirely, which *does* share the same `validate()` gate but not a distinct count), a CpG on chrX/chrY (both of which annotate with `chr = "chrX"`/`"chrY"` per `annotation.R:80-86`, i.e. never a bare integer) is always excluded from this specific plot even when it has a valid chr/pos/fdr, with no separate note distinguishing "no annotation at all" from "annotated but on a sex chromosome, which this plot's x-axis (autosome index) cannot place." The results *table* and the *volcano* plot are unaffected — this is scoped to the Manhattan plot only. See §21, Finding LOW-1.

### `mod_methyl_dmp_topplot()`

**Package:** `ggplot2` (project wrapper). **Location:** `mod_methyl_dmp.R:225-243`.
**Purpose:** Horizontal bar chart of the top-N ranked CpGs by the user's chosen ranking rule.
**Input:** `df`, `rank_by` (`"fdr"`/`"dbeta"`/`"combined"`), `n`.
**Operation:** Drops rows missing `fdr` or `dbeta`; `validate(need(...))` if none remain; orders by FDR ascending (ties broken by \|Δβ\| descending), or \|Δβ\| descending (ties by FDR), or `fdr / max(|dbeta|, 1e-6)` ascending (a combined score that rewards low FDR *and* large effect, `"combined"`); takes the top `n`; builds a `gene (cpg)` label when a gene is annotated, else the bare CpG ID; bars colored by direction of Δβ.
**Output:** `list(plot = <ggplot>, cpgs = <character vector of the plotted CpG IDs>)` — note this returns *two* things, not just a plot.
**Role in DMP Pipeline:** Both a visualization and the CpG-selection mechanism for the next plot (β-value distribution), which is fed the `cpgs` element.
**Audit:** No implementation issue identified from the inspected code. The `"combined"` ranking's `pmax(abs(dbeta), 1e-6)` denominator guard correctly avoids division-by-zero for a CpG with exactly zero Δβ.

### `mod_methyl_dmp_betadist()`

**Package:** `ggplot2` (project wrapper). **Location:** `mod_methyl_dmp.R:245-261`.
**Purpose:** Per-group boxplots of raw β values for the top-ranked CpGs (from the function above).
**Input:** `beta_mat` (the analysis's own beta-scale matrix), `cpgs` (top-ranked CpG IDs), `grp` (the two-level group factor).
**Operation:** Intersects `cpgs` with `rownames(beta_mat)` (defensive — a top-ranked CpG could in principle be absent, e.g. if the beta-scale matrix and the ranked table ever diverged, though in this module's own call site they do not); `validate(need(...))` if the intersection is empty; reshapes to long format; boxplots per CpG, dodged/filled by group.
**Output:** a `ggplot` object.
**Role in DMP Pipeline:** Lets the user visually confirm that a statistically significant Δβ corresponds to a real, visible separation in raw methylation values, not just a numerically small but "significant" shift.
**Audit:** No implementation issue identified from the inspected code.

### `mod_methyl_dmp_server()` (orchestrator)

**Package:** `shiny` (project function). **Location:** `mod_methyl_dmp.R:263-875`.
**Purpose:** The `moduleServer()` body wiring every reactive above to the UI; not a single-purpose function, so documented here as an orchestrator rather than against the single-function template.
**Key reactive graph (SVA tab):** `default_data()` (reactive) → `sva_run()` (eventReactive on the Run button) → `default_filtered()` (reactive) → `default_volcano`/`sva_valueboxes_ui`/`default_table`/`download_default` (outputs).
**Key reactive graph (DMP tab):** `sex_col()`/`id_cols()`/`anno_result()`/`dataset_norm_status()` (cheap reactives, recomputed whenever their inputs change) → `live_result()` (eventReactive on the Run button; the only expensive step) → `live_filtered()`/`live_sig_count()`/`live_top()` (cheap reactives over the already-fit result) → seven outputs (`live_results_ui`, `live_qq`, `live_volcano`, `live_manhattan`, `live_topplot`, `live_betadist`, `live_table`) plus two `downloadHandler`s.
**Audit:** See §20–§22 for the scientific/code audit of this orchestration; two specific defensive-programming patterns are worth noting here as **positive** findings, not problems: (1) `live_has_run` is explicitly reset to `FALSE` whenever `methyl_dataset$beta` changes identity (`mod_methyl_dmp.R:530-531`), preventing a stale result from a previous dataset from being left on screen after a new dataset loads; (2) both `default_table` and `live_table` explicitly set `outputOptions(..., suspendWhenHidden = FALSE)` (`mod_methyl_dmp.R:397,847`) to work around a reproduced, documented Shiny visibility-detection bug specific to a table nested two `renderUI()` levels deep — without this, the table's data would compute correctly server-side but never actually reach the browser, a real (already-diagnosed) bug this code fixes rather than a hypothetical one.

### 9.13 External functions (compact reference)

| Function | Location | Role in DMP |
|---|---|---|
| `load_default_dmp(stage, sex)` | `global.R:313-320` | Reads one of the four precomputed SVA-tab CSVs; returns `NULL` gracefully if the preloaded results folder isn't present in this deployment. |
| `load_default_meth_pheno()` | `global.R:340-343` | Reads the preloaded cohort's sample metadata for the SVA tab's header counts. |
| `methyl_get_annotation(array_type)` | `annotation.R:48-93` | Illumina manifest lookup (chr/pos/gene/SNP columns), cached per array type for the process lifetime; returns `list(ok=FALSE, reason=...)` gracefully for array types with no installed annotation package (EPICv2/WGBS/RRBS/Custom). |
| `methyl_filter_missing(mat, max_na_frac)` | `qc.R:30-34` | Row-wise missingness filter, reused unmodified from the QC module. |
| `methyl_filter_variance(mat, min_variance)` | `qc.R:36-41` | Row-wise variance filter; delegates to `methyl_row_vars()` (`qc.R:22-28`, `matrixStats::rowVars()` when available). |
| `methyl_filter_snp(mat, anno_result)` | `qc.R:68-81` | Flags probes overlapping a manifest-annotated SNP (`Probe_rs`/`CpG_rs`/`SBE_rs`); a no-op with an explanatory note if annotation is unavailable. |
| `methyl_sheet_sample_ids(sheet, all_ids)` | `qc.R:456-465` | Resolves sample-sheet rows to matrix column IDs by column-name match, falling back to positional row order only when row counts match. |
| `methyl_norm_status(mat, dataset, anno_result)` | `normalization.R:424-455` | Read-only diagnostic message shown in the live DMP tab's header card (`mod_methyl_dmp.R:426-429,467`); DMP does not act on this status (no gating, no forced normalization) — purely informational. |
| `stats::model.matrix`, `qr` | base R | Design matrix construction and rank check. |
| `limma::lmFit`, `limma::makeContrasts`, `limma::contrasts.fit`, `limma::eBayes`, `limma::topTable` | `limma` | The statistical core — see §12. |

---

## 10. Statistical Methodology

**General concept vs. this implementation, made explicit per tab:**

| Concept | SVA tab | DMP tab (live) |
|---|---|---|
| Where the model is fit | Offline, outside the app (documented in `METHODS_dmp_sva_sexstratified.md`) | Live, inside this Shiny session, on Run-button click |
| Response scale | M-values (offline pipeline) | M-values (`mod_methyl_dmp.R:660`, fit on `m`) |
| Design | `~ group + age + smoking + 6 cell-type fractions + SVA surrogate variables` (fixed, per the offline methods doc) | `~0 + group [+ user-chosen covariates]` (user-configurable) |
| Bias/inflation correction | `bacon::bacon()`, offline | **None** — raw `eBayes()`/`topTable()` output only |
| Multiple testing | Benjamini-Hochberg, on bacon-corrected p-values, offline | Benjamini-Hochberg (limma's `topTable()` default), on raw moderated p-values |
| Effect size | Δβ (beta scale), offline | Δβ (beta scale), computed live at `mod_methyl_dmp.R:667-669` |

**Scientific background — why bacon/SVA matter.** Genome-wide array-based EWAS are prone to test-statistic inflation from unmodeled technical or population-structure confounding; a genomic inflation factor (λ) noticeably above 1 signals that the raw p-value distribution is not calibrated to the null, which — left uncorrected — produces an excess of false-positive "significant" CpGs. SVA (Leek & Storey, 2007) and bacon (van Iterson et al., 2017) are two established, complementary responses to this: SVA estimates and adjusts for latent confounding structure *before* model fitting; bacon rescales the resulting test statistics' empirical null distribution *after* fitting.

**Code evidence — the live DMP engine's own UI states this limitation to the user, not just to this audit.** `mod_methyl_dmp.R:757-759`: when `lambda_gc > 1.1`, the app displays *"This live engine does not apply SVA/bacon correction... Treat the p-values/FDR above as potentially optimistic and cross-check against the default/SVA-corrected panel before drawing conclusions."* This is treated in this document as a **documented, intentional design boundary**, not an undisclosed defect (see §13 for the audit's own independent assessment of the consequence).

---

## 11. Data Processing

Covered in full, step-by-step, in §6.2 ("Processing", 22 steps). Summary of the processing *order*, which the code's own comments state was deliberately chosen for memory safety on a genome-wide array:

```
Load beta matrix + sheet
      → match samples (intersect, ≥6)
      → subset by sex (optional, ≥6 remain)
      → subset by group ∈ {ref, comp}
      → free the pre-subset matrix (rm)
      → subset by covariate complete-cases (optional, ≥6 remain)
      → check ≥3 samples per group
      → filter probes (missingness / variance / SNP) — BEFORE the M-value transform, on the smaller already-subsetted matrix
      → derive M-values (and a beta-scale copy) from the filtered subset only
```

---

## 12. Differential Methylation Calculation

Design: `model.matrix(~0 + grp [+ covariates])`. Fit: `methyl_chunked_lmfit()` (chunked `limma::lmFit()`). Contrast: `comparison − reference`. Moderation: `limma::eBayes()`. Extraction: `limma::topTable(..., number = Inf, sort.by = "P")`. Effect size: independently computed Δβ on the beta scale. All at `mod_methyl_dmp.R:644-669`; fully detailed in §6.2 steps 11–15.

## 13. Multiple-Testing Correction

**SVA tab:** Benjamini-Hochberg FDR, applied offline (per `METHODS_dmp_sva_sexstratified.md`), to bacon-corrected p-values, independently per sex stratum. Displayed as `fdr_bacon`.

**DMP tab (live):** Benjamini-Hochberg FDR via `limma::topTable()`'s own default `adjust.method` (not overridden — verified by inspecting the exact call at `mod_methyl_dmp.R:665`, which passes only `number` and `sort.by`), applied to the raw moderated p-values from `eBayes()`. Displayed as `fdr`. **No genomic-control or bacon rescaling precedes this correction** — the FDR is computed directly on whatever p-value distribution `eBayes()` produced, which the app's own UI (§10) tells the user may be inflated.

**Threshold used by the implementation:** user-adjustable `live_fdr`/`sva_fdr`, defaulting to the field's conventional `0.05` in both tabs (`mod_methyl_dmp.R:301,480`) — not hardcoded, not silently different from what the UI displays.

## 14. Annotation

**DMP tab only** (the SVA tab's precomputed CSVs carry no chr/pos/gene columns). `methyl_get_annotation()` (`annotation.R:48-93`) reads the Bioconductor annotation package's own `Locations`/`Manifest`/`Other`/`SNPs.147CommonSingle` data objects directly (deliberately bypassing `minfi::getAnnotation()`'s wrapper object for a documented namespace-pollution reason, `annotation.R:38-47`), producing chromosome, position, probe type, three SNP-overlap columns, and a single representative gene symbol per probe (first token of the semicolon-separated UCSC RefGene list). Available only for `450K` and `EPIC` (v1) array types in this deployment (`annotation.R:18-21`); `EPICv2`/`WGBS`/`RRBS`/`Custom array` degrade gracefully to `ok = FALSE` with an explanatory reason string, never a guess.

## 15. Visualization

Fully covered per-tab in §6.1/§6.2 ("Figures"). Six distinct figure types across both tabs (volcano ×2 tabs, QQ, Manhattan, top-DMP bar chart, β-distribution boxplot).

## 16. Results and Downloads

Fully covered per-tab in §6.1/§6.2 ("Output"/"Downloads"). Three `downloadHandler`s total: `download_default`, `download_live`, `download_live_config`.

---

## 17. Tab-to-Tab Connections

**Are the two sub-tabs sharing computation?** No. "SVA" and "DMP" are structurally independent — different data sources (precomputed CSV vs. a live matrix), different reactive chains (`sva_run`/`default_*` vs. `live_result`/`live_*`), different UI outputs, no shared reactive object between them beyond the fact that both live inside the same `moduleServer()` closure (`mod_methyl_dmp.R:263-875`). Switching the `tabsetPanel` between them never triggers either side's computation (`mod_methyl_dmp.R:53-56`).

**Does one tab need the other to have run first?** No. Each tab renders (or shows its own empty-state message) purely from `methyl_dataset` — "SVA" from `methyl_dataset$preloaded`, "DMP" from `methyl_dataset$beta` — with no dependency on the other tab's button having been clicked.

**Is anything reactive vs. click-gated?** Both tabs' *result-producing* computation is strictly click-gated via `eventReactive` (`sva_run` on `sva_run_btn`; `live_result` on `live_run_btn`) — changing a filter/sex/group/covariate control alone never silently recomputes or reruns the model. Within the DMP tab specifically, the *display* filters (FDR/Δβ/direction: `live_filtered`, `live_sig_count`) **are** ordinary `reactive()`s, so they *do* update live as those three controls change, without needing another click — but only over the already-fit `live_result()`, never by refitting.

**What happens if an upstream result is missing?** "SVA": `req(isTRUE(methyl_dataset$preloaded))` then `req(d)` (`mod_methyl_dmp.R:283,285`) — renders nothing at all if either fails, no placeholder card. "DMP": explicit conditional cards for "no matrix loaded" and "no sample sheet" (`mod_methyl_dmp.R:439-454`) — always shows *something* explaining what's missing, never a blank panel.

**Cross-module connection (noted, not documented further — out of scope per §2).** `methyl_results$dmp` is written unconditionally after every successful "DMP" tab run (`mod_methyl_dmp.R:707-713`) and is read by `mod_cross_integration.R:244-250` (Cross-Omics module, out of scope) purely to display a one-line "Methylomics DMP: `<comparison>`, `<n_probes>` probes, `<n_sig>` significant" status string. This is a one-directional, read-only status publication — the Cross-Omics module does not feed anything back into, or otherwise alter, the DMP submodule.

## 18. End-to-End Data Flow

**SVA tab:**
```
Precomputed CSV (offline limma+SVA+bacon pipeline, outside this app)
      ↓  load_default_dmp(stage, sex)
Stratum + FDR + Δβ + Direction selection (UI controls, not yet applied)
      ↓  click "View results"  (eventReactive gate)
sva_run()  — selects the SVA-stage table for the chosen stratum
      ↓
default_filtered()  — mod_methyl_dmp_filter()
      ↓
Volcano plot  +  Value boxes  +  DT table  +  CSV download
```

**DMP tab:**
```
methyl_dataset$beta + $sample_sheet  (upload, or preloaded live matrix if present)
      ↓
Sex / Group column / Reference / Comparison / Filters / Covariates  (UI controls, not yet applied)
      ↓  click "Run DMP Analysis"  (eventReactive gate)
Sample matching → sex subset → group/ref/comp subset → covariate complete-case subset
      ↓
Probe filters (missingness / variance / optional SNP)  — on beta scale, before M-transform
      ↓
M-value (and beta-scale) derivation on the filtered subset
      ↓
Design matrix (~0 + group [+ covariates])  →  rank check
      ↓
methyl_chunked_lmfit()  →  makeContrasts()  →  contrasts.fit()  →  eBayes()
      ↓
topTable()  (raw p, BH-adjusted FDR)  +  independently computed Δβ
      ↓
Annotation (chr/pos/gene, if available)  +  direction label  +  λ diagnostic
      ↓
live_result()  — published to methyl_results$dmp
      ↓
FDR/Δβ/Direction display filter (live, no refit)
      ↓
QQ + Volcano + Manhattan + Top-DMP + β-distribution plots  +  DT table  +  2 CSV downloads
```

No normalization, batch correction, or SVA step exists anywhere in the DMP-tab flow above — their absence is deliberate and stated in-app (§10), not an omission this audit is disclosing for the first time.

## 19. End-to-End Pipeline (User Perspective)

1. User opens Methylomics → Differential Methylation (DMPs).
2. If a preloaded dataset is active: "SVA" tab is available immediately, showing dataset-level counts.
3. User picks Stratum/FDR/Δβ/Direction and clicks "View results" → filtered table, volcano, value boxes, download appear.
4. Separately/additionally: if any beta/M matrix (upload or live preloaded matrix) is loaded, "DMP" tab becomes available.
5. User picks Sex, Group column, Reference/Comparison groups, filters, optional covariates.
6. User clicks "Run DMP Analysis" → sample matching, subsetting, probe filtering, M-value transform, design matrix, `limma` fit, contrast, moderation, `topTable()`, annotation, λ diagnostic all execute in one `eventReactive`.
7. Results UI appears: configuration summary, inflation diagnostic + QQ plot, summary value boxes, volcano, Manhattan (if annotated), top-DMP bar chart, β-distribution boxplot, results table.
8. User adjusts FDR/Δβ/Direction post-run (no refit) to explore the already-fit result.
9. User downloads the filtered CSV and/or the analysis-configuration CSV.

---

## 20. Scientific Code Audit

**Data handling.** Matrix orientation (probes × samples) is correct and consistently assumed throughout (§8). Sample matching is name-based, not positional-only (`methyl_sheet_sample_ids()`), with an explicit, narrow positional fallback only when row counts already match. Duplicate sample IDs are **not explicitly checked** — `intersect()`/`match()` against a sheet with a duplicated sample ID would silently succeed with an arbitrary/duplicated match rather than erroring (see §21, Finding MED-1). Missing values are handled via an explicit, user-controlled missingness filter (not median/mean imputation) and via covariate complete-case restriction, both count-reported to the user. Zero-variance probes are removable via the variance filter but the filter defaults to `0` (no-op) — a genuinely zero-variance probe would still enter the model fit if the user leaves this at default (not typically fatal to `limma::lmFit()`, but see §21, Finding LOW-2).

**Statistical design.** Reference/comparison group choice is fully user-controlled and explicit (not inferred). Contrast direction (`comp - ref`) is unambiguous and displayed. `eBayes()`/`topTable()` moderated-t is an appropriate, standard choice for array-based single-CpG EWAS (Maksimovic, Phipson & Oshlack, 2016 — cited by the offline pipeline's own methods doc for the same reason). Covariate handling is correct (complete-case, formula-based, rank-checked). Small-group behavior is guarded (`≥3` per group, `≥6` overall) but not "small" in an absolute statistical-power sense — a 3-vs-3 comparison is permitted and will fit, but is likely underpowered; the UI does warn "results may be underpowered" only below `n < 10` per group (`mod_methyl_dmp.R:747-748`), not below the hard `n ≥ 3` gate itself.

**Multiple testing.** Raw and adjusted p-values are both computed and both retained in the output table (`p_raw`, `fdr`) — never only one. The FDR method (BH) is standard and appropriate for genome-wide EWAS. Filtering (FDR + Δβ + direction) is applied identically via one shared function on both tabs, avoiding drift between them.

**Biological interpretation.** CpG-to-gene annotation is applied only when a manifest is actually available for the array type; never guessed. Direction of methylation change (`hyper`/`hypo`) is computed directly from the sign of a real, independently-computed Δβ, not inferred from the model's M-value coefficient sign (a good practice, since a probe's Δβ and ΔM signs can occasionally disagree near the scale's extremes due to the nonlinearity of the beta↔M transform — computing Δβ directly, as this code does, avoids that edge case entirely rather than assuming M-scale direction always matches beta-scale direction).

**Visualization.** Volcano/Manhattan/top-DMP plots all correctly encode effect size and significance from the actual result table (not recomputed or approximated for plotting). The Manhattan plot's silent chrX/chrY exclusion (§9, `mod_methyl_dmp_manhattan` Audit) is the one visualization-accuracy caveat identified.

## 21. Code Audit Findings

### HIGH

None identified. No finding in this audit rises to "can materially change scientific results or cause the analysis to fail" without already being an explicitly surfaced, in-app-documented limitation (the absence of bias/inflation correction in the live engine, §10/§13, is exactly this kind of consequential-but-disclosed limitation, categorized as MEDIUM below because the app itself already tells the user, in the results UI, to treat the numbers as potentially optimistic).

### MEDIUM

**MED-1 — No duplicate-sample-ID check before matching.**
- File/Function: `mod_methyl_dmp.R:541-543`, `methyl_sheet_sample_ids()` (`qc.R:456-465`).
- What the code does: `intersect(colnames(beta), methyl_sheet_sample_ids(sheet, colnames(beta)))`, then `match(common, sample_ids)` to align `ph0` to the matrix's column order.
- Why it matters: `match()` returns the *first* matching index for a duplicated key. If the sample sheet contains two rows with the same sample ID (e.g. a re-scanned or accidentally duplicated GEO/upload row), `ph0` would silently take the first row's phenotype for that sample without any error or warning, potentially assigning the wrong group/covariate values to that sample.
- Scientific consequence: A silently mis-assigned sample can shift group means and inflate or deflate the moderated t-statistic for every probe, in a way that would not be visible from the output alone.
- Recommended correction: Add an explicit `validate(need(!anyDuplicated(sample_ids_matched), "Duplicate sample IDs detected in the sample sheet — resolve before running a DMP analysis."))` before the `match()` step.
- Should this be made now: This is an audit/documentation task; no correction was made. Flagging for the developer/thesis author to decide.

**MED-2 — Live engine's disclosed lack of bias/inflation correction has no enforced cross-check.**
- File/Function: `mod_methyl_dmp.R:660-665` (model fit), `:757-759` (the app's own warning text).
- What the code does: Fits and reports FDR-significant CpGs with no bacon/SVA/genomic-control adjustment; only *displays* a λ-based warning when λ > 1.1, with no gate that prevents viewing/downloading results regardless of λ.
- Why it matters: A user who does not read the warning text can download and cite a DMP list from an inflated model with no software-enforced friction.
- Scientific consequence: Potential inclusion of false-positive CpGs in a downstream candidate/biomarker list if the live engine is used on a dataset/comparison with real unmodeled confounding, and the warning is missed.
- Recommended correction: None recommended beyond what already exists — a hard block on high-λ results would over-correct (λ > 1.1 is a soft heuristic, not proof of a false result, and some genuinely large, diffuse biological effects can also elevate λ). The existing soft warning is a reasonable, standard field practice; documented here as a **known, disclosed limitation** rather than a defect requiring a code change.
- Should this be made now: No — this is the module's own intentional, disclosed design choice, not an audit-discovered bug.

### LOW

**LOW-1 — Manhattan plot silently drops sex-chromosome CpGs with no distinguishing count shown.**
- File/Function: `mod_methyl_dmp_manhattan()`, `mod_methyl_dmp.R:207-209`.
- What the code does: `as.integer(gsub("^chr", ...))` converts `"chrX"`/`"chrY"` to `NA`, which is then dropped by the same `!is.na(d$chr_num)` filter used for genuinely unannotated CpGs, with one combined `validate()` message and no separate count for "annotated but on a sex chromosome" vs. "not annotated at all."
- Why it matters: A user could misread "no chrX/Y points on the Manhattan plot" as "no chrX/Y association was tested," when in fact chrX/Y CpGs were tested and are visible in the results table and volcano plot — only this one figure omits them.
- Scientific consequence: Interpretation risk only (a plot-reading pitfall), not a computation error — the underlying statistics are correct and available elsewhere in the same tab.
- Recommended correction: A one-line caption/footnote on the Manhattan plot noting that sex-chromosome CpGs are excluded from this specific figure.
- Should this be made now: No — documentation/audit task; flagged for the developer's discretion.

**LOW-2 — Default variance filter threshold is a no-op.**
- File/Function: `mod_methyl_dmp.R:485`, `methyl_filter_variance()` (`qc.R:36-41`).
- What the code does: `live_min_variance` defaults to `0`, and `methyl_filter_variance(mat, 0)` keeps every probe with variance ≥ 0, i.e. every probe including any that happen to be exactly zero-variance within the current sample subset.
- Why it matters: A zero-variance probe (identical value in every sample of the current subset) contributes no information and, in rare edge cases, can affect the numerical stability of a linear-model fit.
- Scientific consequence: Minor — `limma::lmFit()` is generally robust to this, and such probes would simply return an uninformative (near-infinite standard error, non-significant) result rather than crashing the fit.
- Recommended correction: None necessary; documented as a design choice (defaulting to "no filter" rather than a nonzero default) worth being aware of, not a defect.
- Should this be made now: No.

### INFORMATIONAL

**INFO-1 — "Plain" (unadjusted) precomputed table is loaded but never independently viewable.**
- File/Function: `mod_methyl_dmp.R:273,290-291,315-316`.
- Observation: `default_data()` loads `plain_f`/`plain_m` for both sexes, but the only use of those two tables anywhere in the file is to compute `n_sig_plain_f`/`n_sig_plain_m` for one sentence of static prose. There is no UI path to view, filter, plot, or download the plain-stage table itself from this module.
- Not necessarily a problem: this is explicitly consistent with the module's own stated design ("the SVA-adjusted model above resolves it and is the panel actually used downstream," `mod_methyl_dmp.R:315-316`) — the plain stage is intentionally background context, not a parallel interactive panel. Documented here purely as a **Declared vs. Actually executed** distinction per the task's own instruction, not as a defect.

**INFO-2 — `sva` and `bacon` are named by the tab/pipeline but not called live.**
- Already covered in full in §4 and §10; repeated here as a formal finding for completeness: the "SVA" tab name refers to an *offline* pipeline stage, not a live `sva::sva()`/`bacon::bacon()` call inside `mod_methyl_dmp.R`. No mislabeling was found — the module's own comments are explicit and accurate about this distinction (`mod_methyl_dmp.R:10-23`); this finding exists to make sure a reader of *this document alone*, without having read the source comments, does not draw the wrong conclusion from the tab's name.

## 22. Unused / Dead / Misleading Code

- **`ggrepel`** is attached app-wide (`global.R:65`) but genuinely unused by any function in `mod_methyl_dmp.R` (verified: zero `geom_text_repel`/`geom_label_repel`/`ggrepel::` occurrences). Not dead code *within this file* (it's a shared app-wide library load, not something this module itself declared and failed to use) — noted for completeness only.
- **Every UI control has corresponding server logic and vice versa** — verified by cross-checking every `input$live_*`/`input$sva_*` reference against its `ns(...)` declaration in the UI functions; no orphaned input or output was found.
- **No download button is missing an implementation** — all three `downloadButton`/`downloadHandler` pairs are present and wired (`default_table`↔`download_default`; `live_table`↔`download_live`/`download_live_config`).
- **No analysis option silently executes a different method than labeled** — the "Model" text shown to the user (`mod_methyl_dmp.R:741`, `tags$code(r$design_formula)`, "limma moderated t-test, eBayes") accurately reflects what was actually fit; the design formula string is built from the same `input$live_group_col`/`covariate_cols` values actually used in `model.matrix()`, not a separately-maintained (and potentially stale) copy.
- Both `default_table`/`live_table`'s `outputOptions(..., suspendWhenHidden = FALSE)` calls (§9, `mod_methyl_dmp_server` Audit) are themselves evidence of a *previously* real dead-output bug (data computed server-side, never sent to the browser) that has already been fixed in the inspected code — not a currently-live issue, documented here because the code's own comment (`mod_methyl_dmp.R:385-396`) explicitly discusses it and a thesis reader deserves to know it was found and fixed, not merely absent.

## 23. Error Handling and Validation

Every user-facing failure mode in the live engine uses `validate(need(...))`, which renders a plain-text message in place of the output rather than an uncaught R error/stack trace — verified present for: no matrix loaded, no sample sheet, no group column chosen, no ref/comp chosen, ref == comp, <6 matched samples, <6 samples after sex restriction, <6 complete-case samples after covariate restriction, <3 samples in either group, <10 probes surviving filters, rank-deficient design, `makeContrasts()` failure (e.g. invalid group-name characters), `eBayes()`/`contrasts.fit()` failure. This is a comprehensive, layered validation chain — no step downstream of a failed check is reachable (`mod_methyl_dmp.R:534-658`). The SVA tab has a much shorter validation surface, appropriate to its narrower scope (post-hoc filtering of an already-valid precomputed table): `req(d)` for data availability (`mod_methyl_dmp.R:285`).

## 24. Reproducibility Audit

- **Package dependencies:** documented in §4; `ArthOMix/DESCRIPTION` pins minimum versions app-wide (e.g. `limma (>= 3.62.2)`), though this DESCRIPTION file is explicitly stated (its own header, `ArthOMix/DESCRIPTION:6-9`) to document the dependency surface only, not to drive an actual package build — exact reproducibility depends on the accompanying `renv.lock`, not inspected as part of this file-scoped audit.
- **Random seeds:** **Not implemented in the current code.** No `set.seed()` call exists anywhere in `mod_methyl_dmp.R`. This is scientifically appropriate here: nothing in the live DMP engine's own computation (`lmFit`/`eBayes`/`topTable`) is stochastic — no permutation, bootstrap, or randomized algorithm is used, so no seed is needed for exact reproducibility of a given run's numeric output.
- **Hard-coded thresholds:** all statistical thresholds (FDR, Δβ, min valid %, min variance) are user-adjustable `numericInput`s with documented defaults (§6.1/§6.2 tables), not hardcoded into the computation.
- **Hard-coded group names:** none — group/reference/comparison are fully derived from whatever the loaded sample sheet contains (§6.2, code evidence `mod_methyl_dmp.R:29-36`).
- **Hard-coded paths:** none within `mod_methyl_dmp.R` itself; `METH_DMP_PLAIN_DIR`/`METH_DMP_SVA_DIR` are defined in `data_paths.R` (out of scope), read only via `load_default_dmp()`.
- **The `download_live_config` export is this module's own reproducibility mechanism** (§6.2 Downloads item 2) — it records every parameter needed to describe a live run (dataset, sex, groups, sample sizes, covariates, method, design formula, all filter thresholds, probe/significance counts, timestamp) in one CSV, specifically so a given run's configuration is not lost once the session ends.
- **Can the analysis be reproduced from the code alone?** The **live DMP tab**: yes, given the same input matrix/sheet and the same configuration (recorded in the config-export CSV), the computation is deterministic. The **SVA tab**: **no** — its numbers originate from an offline script not present in `mod_methyl_dmp.R` or anywhere else in the Shiny application's own R/ directory; only its documented methods (`METHODS_dmp_sva_sexstratified.md`) are available to a reader of this codebase, not the script itself.

## 25. Educational Plain-English Walkthrough

The SVA tab is a **viewer**, not a calculator: someone already ran a careful statistical analysis — comparing RA patients to healthy controls, separately for men and women, adjusting for age, smoking, blood-cell-type mixture, and hidden technical noise (via SVA), then double-checking the result wasn't inflated (via bacon) — and saved the results to a spreadsheet. This tab reads that spreadsheet, lets you pick a sex, set how strict you want to be about significance (FDR) and how big a methylation change you care about (Δβ), and shows you a plot and table of only the CpGs that pass those cutoffs. You can download exactly what you're looking at.

The DMP tab is a **calculator**: you provide your own methylation matrix and phenotype sheet (or use the app's preloaded one, if the deployment has it), and you tell the app which two groups to compare, whether to look at everyone or just one sex, which extra variables (like age) to adjust for, and how strict your quality filters should be. When you click Run, the app matches your samples to your matrix, throws out probes with too much missing data or too little variation, converts the 0–1 methylation values to a statistically better-behaved scale (M-values), fits a linear model per CpG comparing your two groups, uses `limma`'s empirical-Bayes trick to stabilize the per-CpG variance estimates (this is what makes it more reliable than a plain t-test with few samples), and reports a raw p-value and an FDR-adjusted p-value for every CpG tested. It also tells you, honestly, whether the overall result pattern looks inflated (the λ diagnostic) — because, unlike the SVA tab, this live calculator does not apply the extra correction steps that would normally be used to fix that. You then get six different plots and a searchable, downloadable results table to explore what it found.

## 26. Input → Function → Output Summary

| Tab | Input | Input Type | Function(s) | Processing | Output | Output Type |
|---|---|---|---|---|---|---|
| SVA | Stratum, FDR, Δβ, Direction (+ "View results" click) | UI controls | `load_default_dmp`, `mod_methyl_dmp_filter`, `mod_methyl_dmp_volcano` | Select precomputed table for stratum → filter | Filtered per-CpG table, volcano plot, value boxes, CSV | `data.frame` (7 cols), `ggplot`, CSV |
| DMP | Beta/M matrix, sample sheet, Sex, Group column, Ref/Comp, FDR, Δβ, Direction, Min valid %, Min variance, SNP filter, Covariates (+ "Run" click) | Matrix + UI controls | `methyl_sheet_sample_ids`, `methyl_filter_missing`, `methyl_filter_variance`, `methyl_filter_snp`, `methyl_chunked_lmfit`, `limma::makeContrasts/contrasts.fit/eBayes/topTable`, `methyl_get_annotation`, `mod_methyl_lambda_gc` | Match → subset → filter probes → M-transform → design → fit → contrast → moderate → extract → annotate | Per-CpG results table (`cpg,t,p_raw,fdr,dbeta,ref_mean_beta,comp_mean_beta,chr,pos,gene,direction`), λ + QQ, volcano, Manhattan, top-DMP chart, β-distribution chart, 2 CSVs | `data.frame` (11 cols), 5×`ggplot`, 2×CSV |

---

## 27. Thesis Implementation Paragraph

The Differential Methylation (DMPs) sub-module of the Methylomics component provides two complementary analysis pathways, organized as two tabs ("SVA" and "DMP"). The first reproduces and interactively filters a precomputed, sex-stratified, surrogate-variable-adjusted, bacon-corrected `limma` differential methylation analysis of the bundled whole-blood cohort. The second is a live, fully configurable engine that fits `limma::lmFit()` with empirical-Bayes moderation (`eBayes()`) to a user-selected or preloaded beta/M-value matrix, using a user-defined reference/comparison contrast, optional covariates, and configurable missingness/variance/SNP probe filters, before reporting Benjamini-Hochberg FDR-adjusted significance and beta-scale effect size per CpG. Differential methylation is identified as a moderated t-test on M-values per probe; significance is FDR-controlled; a genomic-inflation diagnostic (λ) is reported because the live engine, unlike the offline stage it complements, applies no bias-correction step. Outputs include annotated per-CpG result tables, volcano, QQ, Manhattan, and effect-size visualizations, and full reproducibility exports. This dual design lets a user compare a validated reference analysis against an independently configurable model on the same underlying biology.

## 28. Very Short Thesis Version

This sub-module identifies differentially methylated CpG positions between two sample groups using `limma`'s moderated t-test framework (`lmFit`/`eBayes`/`topTable`) with Benjamini-Hochberg FDR correction, applied either to a precomputed, surrogate-variable-adjusted reference analysis of the bundled cohort or live to a user-loaded beta/M-value matrix with user-defined groups and covariates. It reports per-CpG t-statistics, raw and adjusted p-values, and beta-scale effect sizes, annotated to chromosome, position, and gene where array-manifest data permit, supporting hypothesis generation for disease-associated epigenetic variation.

---

## 29. Overall Assessment

**What the module does well.** A clear, code-documented separation between a validated reference reproduction and a genuinely general-purpose live engine; extensive, layered input validation with human-readable failure messages rather than crashes; a memory-conscious processing order (filter-before-transform, chunked model fitting, explicit intermediate-object cleanup) specifically engineered against a real, previously-reproduced out-of-memory failure on the genome-wide preloaded matrix; an honest, in-UI disclosure of the live engine's own calibration limitation rather than a silent omission; a full reproducibility export of every run's configuration.

**What the module actually implements.** Single-CpG (not region-level) differential methylation testing via `limma`'s standard moderated-t framework; BH FDR correction; Δβ and direction as descriptive effect sizes; optional covariate adjustment and sex stratification; missingness/variance/SNP probe filtering; manifest-based annotation for 450K/EPIC arrays only.

**Important assumptions.** That the loaded matrix's rows are genuinely probes/CpGs and columns are genuinely samples (never independently verified against, e.g., a transposed upload). That `methyl_dataset$input_scale` correctly reflects the matrix's true scale (not independently re-derived from the data's numeric range). That sample-sheet IDs are unique (§21, MED-1).

**Important limitations.** No SVA/bacon/genomic-control correction in the live engine (disclosed in-app). No batch-effect handling anywhere in this file. Manhattan plot silently excludes sex-chromosome CpGs. The "plain" (unadjusted) precomputed stage is loaded but not independently browsable.

**What should be verified before thesis submission.** (1) Confirm with the offline pipeline's own author/logs that the `METHODS_dmp_sva_sexstratified.md` methods text this document relies on for §6.1/§10 is still an accurate description of the currently-bundled CSV files (this audit did not re-run that offline pipeline). (2) Decide, and state explicitly in the thesis text, whether the live engine's disclosed lack of bias correction is an acceptable limitation for whatever specific comparisons the thesis reports from it, or whether those specific results should instead be cross-validated against a bacon/SVA-corrected re-analysis. (3) If citing sample sizes or significant-CpG counts from the "SVA" tab in the thesis body, cite them from `METHODS_dmp_sva_sexstratified.md`'s own Results section (§2.AA.8, reproduced in part in §6.1 above) rather than re-deriving them, since that document is the authoritative source this tab merely displays.
