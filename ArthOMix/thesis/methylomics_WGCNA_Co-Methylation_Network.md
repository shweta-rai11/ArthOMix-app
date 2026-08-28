# Methylomics → WGCNA (Co-Methylation Network): Complete Technical Audit and Teaching Document

**Scope of this document.** This is a read-only code audit and teaching reference for exactly one submodule: **Methylomics → WGCNA (Co-Methylation Network)**. No application code was modified to produce this document. Every claim below is traceable to a specific file and line number; where general WGCNA theory is described for teaching purposes without a corresponding line of code, it is explicitly labeled as background knowledge, not an implementation claim.

**Primary source file (read in full, 1,097 lines):** [`ArthOMix/R/methylomics/mod_methyl_wgcna.R`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R)

**Dependency files read in full or in the relevant sections:**
- [`ArthOMix/R/methylomics/qc.R`](../ArthOMix/R/methylomics/qc.R) — `methyl_row_vars()`, `methyl_filter_missing()`, `methyl_sheet_sample_ids()`, `methyl_qc_subgroup_filter()`
- [`ArthOMix/R/methylomics/mod_methyl_dmp.R`](../ArthOMix/R/methylomics/mod_methyl_dmp.R) — `methyl_chunked_lmfit()`, `mod_methyl_dmp_sex_col()`, `mod_methyl_dmp_sex_choices()`, `mod_methyl_dmp_covariate_cols()` (reused by WGCNA, not redefined)
- [`ArthOMix/R/methylomics/normalization.R`](../ArthOMix/R/methylomics/normalization.R) — `methyl_get_norm_annotation()`
- [`ArthOMix/R/methylomics/mod_methyl_dataset.R`](../ArthOMix/R/methylomics/mod_methyl_dataset.R) — the shared `methyl_dataset` reactiveValues structure this module reads from
- [`ArthOMix/global.R`](../ArthOMix/global.R) — `get_or_compute_meth_wgcna_blocks()`, `load_default_wgcna_module_trait()`, `load_default_wgcna_module_assignment()`, `load_default_dmr_biomarker_panel()`, `ARTHOMIX_COLORS`, `ARTHOMIX_STATUS`, `theme_arthomix()`
- [`ArthOMix/data_paths.R`](../ArthOMix/data_paths.R) — `METH_DATA_AVAILABLE`, `METH_RAW_DATA_AVAILABLE`, `METH_WGCNA_DIR`, `METH_WGCNA_CACHE_DIR`
- [`ArthOMix/R/submodules_registry.R`](../ArthOMix/R/submodules_registry.R) — module registration
- [`ArthOMix/server.R`](../ArthOMix/server.R) — invocation of `MX_MODULES`
- `data/preloaded/methylomics/tables/script05_wgcna_sexstratified/METHODS_wgcna_sexstratified.md` (103 lines) — the thesis chapter documenting the **published, offline** sex-stratified WGCNA analysis this module's live tool is designed to reproduce, and whose static output the "Compare with published results" panel displays
- The published run's own output tables (`module_trait_female.csv`, `module_assignment_female.csv`, `biomarker_panel_female.csv`, etc.) were inspected directly to verify column names referenced in the code

**Registration.** `mod_methyl_wgcna_config` — `id = "wgcna"`, `title = "WGCNA (Co-Methylation Network)"`, `icon = "circle-nodes"`, `group = "Network"` ([`mod_methyl_wgcna.R:42-45`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L42-L45)). Registered 6th of 12 in `MX_MODULES`, between DMR (Differentially Methylated Regions) and Candidates ([`submodules_registry.R:45`](../ArthOMix/R/submodules_registry.R#L45)). Invoked generically, like every other Methylomics submodule, via `lapply(MX_MODULES, function(m) m$server(paste0("mx_", m$config$id), methyl_dataset, methyl_results))` in [`server.R:95`](../ArthOMix/server.R#L95) — i.e. `mod_methyl_wgcna_server("mx_wgcna", methyl_dataset, methyl_results)`.

Note on scope: unlike the registry's own header comment, which describes several Methylomics submodules past DMR as "registry-only placeholder scaffolds... queued up to be built out one at a time," `mod_methyl_wgcna.R` is a complete, 1,097-line, fully live implementation — not a placeholder. This audit treats it as such.

---

## 1. Module Overview

### 1.1 What WGCNA is (general background)

**Weighted Gene/Genome Co-expression/Co-methylation Network Analysis (WGCNA)** (Zhang & Horvath, 2005; Langfelder & Horvath, 2008) is a systems-biology method that, given a matrix of some molecular feature measured across many samples, groups the features into **modules** — sets of features whose values rise and fall together across samples — on the basis of their pairwise correlation structure alone, without reference to any outcome label. Originally developed for gene-expression microarrays, the same algorithm and R/Bioconductor package have been applied directly to DNA methylation array data under the same name, with CpG probes substituted for genes as the network's nodes (e.g. Men et al., 2017, cited in the published methods document this module reproduces). This is exactly what the term **co-methylation network** means in this context, and exactly what this submodule implements: **CpGs are the network's nodes; samples are the observations used to estimate their pairwise similarity.**

The core idea in four steps, stated generally:
1. Compute a correlation matrix between every pair of features across all samples.
2. Raise that correlation to a power (the **soft threshold**) to produce a weighted **adjacency matrix**, exaggerating strong correlations relative to weak ones so the resulting network has an approximately scale-free (few highly-connected hub nodes, many weakly-connected nodes) topology, a documented property of many real biological networks.
3. Convert adjacency into **topological overlap** (TOM), a similarity measure that also accounts for shared neighbors, then hierarchically cluster `1 - TOM` to detect **modules** — branches of the resulting dendrogram, cut by a tree-cutting algorithm.
4. Summarize each module by its **eigengene** (the first principal component of its member features' values across samples) and correlate that one summary value per sample against a phenotype/trait of interest — the **module-trait relationship**.

### 1.2 Why analyze CpGs this way

A single-CpG differential-methylation test (as run in this application's DMP submodule) or a spatially-local, adjacent-CpG-cluster test (DMR submodule) each ask a different, narrower question than WGCNA does. WGCNA instead asks: are there whole *sets* of CpGs, possibly scattered across the genome, that vary together as a coordinated unit across samples — and does that coordinated unit, as a whole, track a phenotype? This can surface a distributed signal that no single CpG carries strongly on its own, and it collapses the multiple-testing burden from hundreds of thousands of individual probes down to the handful-to-low-hundreds of modules actually detected. This rationale is stated explicitly in the published methods document this module's defaults are designed to reproduce (`METHODS_wgcna_sexstratified.md`, Section 2.CC.1) and is reflected directly in the implementation: module detection groups CpGs (`net$colors`, [`mod_methyl_wgcna.R:601`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L601)) purely from their correlation structure, with no trait column involved until the separate Module-Trait Analysis tab.

### 1.3 What a module eigengene represents — implemented

**Module eigengenes ARE implemented.** `net$MEs` (from `WGCNA::blockwiseModules()`, [`mod_methyl_wgcna.R:591-599`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L591-L599)) is a samples × modules matrix, one column per detected module (named `ME<color>`), each column being the first principal component of that module's constituent CpGs' methylation values across samples — a single per-sample number summarizing the whole module's coordinated methylation level. This is the object the Module-Trait Analysis tab correlates against phenotype ([`mod_methyl_wgcna.R:696-702`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L696-L702)) and the object Hub CpGs correlates individual CpGs against to compute module membership (kME) ([`mod_methyl_wgcna.R:831`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L831)).

### 1.4 Module-trait relationships — implemented

**Implemented**, on the Module-Trait Analysis tab: each module eigengene is correlated (Pearson or Spearman, user's choice) against a single sample-sheet column encoded as numeric-or-binary, with Student's-t-derived p-values (`WGCNA::corPvalueStudent()`) and BH/Bonferroni correction restricted to real (non-grey) modules. See Section 5.5 and Section 13 below for full detail.

### 1.5 Biological questions this submodule can answer

Based strictly on what is implemented: (a) do groups of CpGs across the genome show coordinated methylation patterns in this cohort, independent of any phenotype; (b) does any such group's aggregate methylation level associate with a chosen binary or continuous trait (e.g. disease status); (c) within an associated module, which individual CpGs are most centrally connected (**hub CpGs**, by module membership kME and intramodular connectivity kWithin); (d) do a module's CpGs overlap, more than chance would predict, with an independently-derived single-CpG/region biomarker panel (the Functional Enrichment tab's Fisher's-exact-test check). It does **not** by itself establish causation, and the code's own p-value/FDR machinery only ever supports an *association* claim (Section 17).

### 1.6 Position within the Methylomics workflow

WGCNA sits after Quality Control, Normalization, Cell-Type Deconvolution, DMP (single-CpG), and DMR (region-based) in the `MX_MODULES` registry order, and before Candidates, Feature Selection, MR, Coloc, Diagnostic, and Biomarker Card ([`submodules_registry.R:38-49`](../ArthOMix/R/submodules_registry.R#L38-L49)). It reads the same shared `methyl_dataset` reactiveValues object (`$beta`, `$sample_sheet`, `$input_scale`, `$array_type`, `$preloaded`, `$source`) that the Dataset tab populates ([`mod_methyl_dataset.R:91-100, 208-219`](../ArthOMix/R/methylomics/mod_methyl_dataset.R#L91-L100)) — it does not read anything from the DMP or DMR submodules' own reactive results, and it writes nothing back to `methyl_results` (a `grep` of the file confirms no `results$` assignment anywhere in `mod_methyl_wgcna.R` — see Section 19, Reproducibility Audit). Its only outbound coupling to another submodule's output is the Functional Enrichment tab's read of the DMP/DMR biomarker panel via `load_default_dmr_biomarker_panel()`, which is a static reference table on disk (from the published pipeline), not a live read of the DMR submodule's current reactive state.

### 1.7 Implemented vs. not implemented — explicit summary

| Implemented in this submodule | General WGCNA concept **not** implemented here |
|---|---|
| Missingness + variability (MAD/variance/SD/IQR) CpG filtering | Batch-effect correction (ComBat, SVA) — deliberately excluded, see Section 8 |
| Optional `limma`-based covariate residualization | True reference-free deconvolution of confounders |
| `WGCNA::goodSamplesGenes()` gate | Outlier-sample removal by a WGCNA-specific method (e.g. `WGCNA::adjacency`-based Z-score cutting); only a descriptive dendrogram/PCA is shown (Section 5.2) |
| `WGCNA::pickSoftThreshold()` with a "first-power-reaching-cutoff, else best-observed" rule | Automatic re-widening of the tested power range when no power reaches the fit target — the app reports the shortfall but does not search further on its own |
| `WGCNA::blockwiseModules()` (signed/unsigned/hybrid, TOM types, PAM stage, deep split, merge cut height) | `WGCNA::TOMplot()` / whole-network TOM heatmap visualization — `saveTOMs = FALSE` ([`mod_methyl_wgcna.R:597`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L597)), so no TOM object is retained to plot |
| `WGCNA::moduleEigengenes()` (internally, via `blockwiseModules()`'s `$MEs`) and `WGCNA::orderMEs()` | A standalone `moduleEigengenes()` call exposed as its own tab/output |
| Module–trait correlation + `corPvalueStudent()` + BH/Bonferroni | Multi-level (>2-category, non-numeric) trait support — explicitly rejected with a validation message (Section 5.5) |
| Hub-CpG ranking by kME (module membership) and kWithin (intramodular connectivity via `WGCNA::intramodularConnectivity.fromExpr()`) | Network export (e.g. Cytoscape edge/node lists via `WGCNA::exportNetworkToCytoscape()`) |
| Fisher's-exact-test enrichment of significant modules against a static DMP/DMR biomarker panel | GO/KEGG-style functional enrichment from a CpG→gene mapping |
| A "Compare with published results" static reference panel (preloaded dataset only) | Post-hoc module consolidation (`WGCNA::mergeCloseModules()` search over a `cutHeight` grid, described in the published methods Section 2.CC.6-2.CC.7) — this exists in the **published, offline** pipeline's methodology document but has **no corresponding live control** anywhere in `mod_methyl_wgcna.R` |

---

## 2. Number of WGCNA Sub-Tabs

The UI is assembled as a single `tabsetPanel(id = ns("subtabs"), ...)` inside `main_ui()` ([`mod_methyl_wgcna.R:1056-1065`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L1056-L1065)), containing exactly seven `tabPanel()` entries.

> **Number of WGCNA sub-tabs: 7**

| # | Tab value | Title (with icon) | Line |
|---|---|---|---|
| 1 | `data` | Data & Filtering | [1058](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L1058) |
| 2 | `sampleqc` | Sample QC | [1059](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L1059) |
| 3 | `power` | Soft Threshold | [1060](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L1060) |
| 4 | `modules` | Network & Modules | [1061](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L1061) |
| 5 | `traits` | Module-Trait Analysis | [1062](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L1062) |
| 6 | `hubs` | Hub CpGs | [1063](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L1063) |
| 7 | `export` | Results & Export | [1064](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L1064) |

Below the tab set, one additional, **always-visible, cross-cutting control** exists that is *not* itself a tab: the **Sex Stratum** `radioButtons` card ([`mod_methyl_wgcna.R:1067-1074`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L1067-L1074)), offering "All samples" plus whatever Female/Male (or raw) levels are detected in the sample sheet's sex column. This selection is read at the very start of Data & Filtering's own computation (Section 8.1), not as a separate pipeline stage, because the published analysis this tool reproduces stratifies by sex *before* residualizing and ranking CpGs by variability, not after ([`mod_methyl_wgcna.R:9-12, 251-256`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L9-L12)).

Also always visible with no button: a **status strip** (`status_ui()`, [`mod_methyl_wgcna.R:187-197`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L187-L197)) reporting the loaded dataset's dimensions, or a notice that only metadata (no live matrix) is available.

---

## 3. Two Data Pathways, One Shared Pipeline

The module's own header comment ([`mod_methyl_wgcna.R:8-26`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L8-L26)) states this design contract explicitly, and the code matches it: there is exactly **one** live computational pipeline (Data & Filtering → Sample QC → Soft Threshold → Network & Modules → Module-Trait → Hub CpGs → Results & Export). Two data sources feed it, differing only in (a) the starting matrix and (b) each control's *default* value — every control stays user-editable regardless of source.

| | Preloaded pathway | Uploaded pathway |
|---|---|---|
| Source | `methyl_dataset$beta` from the whole-blood GSE42861 cohort (Liu et al. 2013), loaded via the Methylomics Dataset tab's "Preloaded whole-blood dataset" option | `methyl_dataset$beta` from a user-uploaded beta/M-value CSV/TSV or derived from uploaded IDAT files |
| Default top-N variable CpGs | 20,000 ([`mod_methyl_wgcna.R:239`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L239)) | 5,000 |
| Default residualization covariates | Age + smoking + cell-type columns minus the auto-detected reference cell type ([`mod_methyl_wgcna.R:173-183`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L173-L183)) | none pre-ticked |
| Default correlation method | Pearson ([`mod_methyl_wgcna.R:435`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L435)) | bicor (biweight midcorrelation) |
| Default custom power vector | ON, `1,2,3,4,5,6,7,8,9,10,12,14,16,18,20` — the published script's own exact tested powers ([`mod_methyl_wgcna.R:439, 446`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L439-L446)) | OFF (auto range 1–20 step 1) |
| Default minimum module size | 20 ([`mod_methyl_wgcna.R:557`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L557)) | 30 (general WGCNA best practice) |
| Default sex stratum | Female, if a Female level is detected ([`mod_methyl_wgcna.R:144-148`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L144-L148)) | "All samples" |

This means: the preloaded pathway's *defaults* are tuned to reproduce the specific published analysis in `METHODS_wgcna_sexstratified.md` as closely as a live, re-run pipeline can (Section 4 covers the live/offline distinction), while the uploaded pathway's defaults follow generic WGCNA-on-methylation practice. A user can freely cross these — e.g. run the preloaded dataset with bicor and no residualization — because nothing in the server code hard-switches behavior by `preloaded`; only the *initial* UI-control values differ (verified: every `value = if (isTRUE(methyl_dataset$preloaded)) ... else ...` construct in the file sets a UI default, never branches the actual computation logic itself).

### 3.1 Live vs. published/offline results — a critical distinction

Two entirely separate result sources exist and must not be conflated:

1. **The live pipeline** (this module's own 7 tabs) recomputes everything from `methyl_dataset$beta` on demand, whenever the deployment has a live matrix available (`METH_RAW_DATA_AVAILABLE`, verified true in this deployment — both `data/preloaded/methylomics/matrix/beta_raw.rds` and `pheno.rds` exist on disk).
2. **The published, offline analysis** described in `METHODS_wgcna_sexstratified.md` was run once, outside this Shiny app, by a separate R script (`script05_wgcna_sexstratified/05_wgcna_sexstratified.R`, not present in this repository — only its methods narrative and output CSVs are bundled). Its results are shown, read-only, in the "Compare with published results" panel on the Results & Export tab (Section 5.7), explicitly labeled "not this run's live output" ([`mod_methyl_wgcna.R:948`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L948)).

Because power selection depends on `pickSoftThreshold()`'s random-order-independent but numerically-derived fit statistics, and because `blockwiseModules()`'s block preclustering uses `randomSeed = 1234` ([`mod_methyl_wgcna.R:597`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L597)) rather than whatever seed the original offline script used, a live re-run with identical UI settings to the published defaults is **expected to closely approximate, but is not guaranteed to numerically reproduce exactly**, the static reference tables. The UI never claims otherwise; it only claims to reproduce the *methodology* (see the header comment, [`mod_methyl_wgcna.R:14-23`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L14-L23)).

---

## 4. Tab-by-Tab Documentation

### 4.1 Tab 1 — Data & Filtering

**Purpose.** Turns the raw `methyl_dataset$beta` matrix into the single filtered, optionally residualized, top-variable-CpG matrix (`f$mat`) that every downstream tab consumes. Nothing is computed until the user clicks "Build filtered matrix."

**Input data.** `methyl_dataset$beta` (CpG-rows × sample-columns numeric matrix — beta values `[0,1]` or, if `methyl_dataset$input_scale == "m"`, unbounded M-values), `methyl_dataset$sample_sheet` (optional phenotype data.frame), and the current **Sex Stratum** selection (`input$sex_stratum`, from the cross-cutting control described in Section 2).

**User inputs** ([`mod_methyl_wgcna.R:230-245`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L230-L245)):

| Input ID | What it is | Default | Range | Used by |
|---|---|---|---|---|
| `max_probe_missing` | Max % missing values per CpG | 5 | 0–100, step 1 | `methyl_filter_missing()` |
| `force_transpose` | Checkbox: treat the loaded matrix as sample-rows × CpG-columns and transpose it | `FALSE` | — | `t(mat_full)` before any other step |
| `var_method` | Variability statistic for ranking CpGs | `"mad"` | mad / variance / sd / iqr | `mx_wgcna_top_variable()` |
| `top_n` | Number of most-variable CpGs to keep | 20,000 (preloaded) / 5,000 (uploaded) | ≥100, step 500 | `mx_wgcna_top_variable()` |
| `resid_covariates` | Checkbox group of sample-sheet columns to regress out | age+smoking+cell-type minus reference (preloaded) / none (uploaded) | any column `mod_methyl_dmp_covariate_cols()` returns as a candidate | `stats::model.matrix()` design, `methyl_chunked_lmfit()` |
| `filter_btn` | Action button that triggers the whole stage | — | — | `eventReactive` trigger for `mx_wgcna_filtered()` |

Note: the trait/disease-status column is deliberately never offered as a candidate covariate ([`mod_methyl_wgcna.R:165, 243`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L165)) — `covariate_choices()` excludes `trait_col_default()` and the active sex column from the candidate list — so any group-associated signal is preserved in the residuals for the later Module-Trait Analysis tab to detect, rather than being regressed away.

**Functions used.**

- **`methyl_qc_subgroup_filter()`** (`qc.R:472`) — Subsets `mat_full`'s *columns* (samples) to the chosen sex stratum by matching a sample-sheet ID column (or row order) against the matrix's column names; returns the unfiltered matrix unchanged when "All samples" is chosen. Why needed: makes stratified sex analysis possible without a separate upstream step, and reproduces the published pipeline's per-sex-first design.
- **`methyl_filter_missing()`** (`qc.R:30`) — Computes `rowMeans(is.na(mat))` and keeps CpGs at or below the missingness threshold. Why needed: standard CpG-level QC before any downstream statistic depends on complete data.
- **`methyl_chunked_lmfit()`** (`mod_methyl_dmp.R:143`) — A drop-in, bit-for-bit-verified chunked wrapper around `limma::lmFit()` that fits the model in row (CpG) chunks of 20,000 to avoid holding more than one full-size matrix copy in memory at once. Input: the M-value matrix and a covariate design matrix. Output: an `lmFit`-shaped object (`$coefficients`, etc.) assembled by row-binding/column-binding per-chunk fits.
- **`stats::model.matrix()`** — Builds the covariate design matrix from the selected `resid_covariates` columns of the sample sheet. A rank check (`qr(design)$rank == ncol(design)`) rejects a rank-deficient design (e.g. cell-type fractions summing to a constant) before fitting.
- **`mx_wgcna_top_variable()`** (own file, [`mod_methyl_wgcna.R:65-77`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L65-L77)) — Ranks CpGs (rows) by the chosen variability statistic (via `mx_wgcna_row_mads()`, `methyl_row_vars()`, or `stats::IQR()`) and keeps the top N. Why needed: WGCNA needs a "broad, correlation-structure-rich" but computationally tractable probe set (published methods, Section 2.CC.2); this is the top-N-by-rank companion to `qc.R`'s fixed-threshold filters.
- **`WGCNA::goodSamplesGenes()`** ([`mod_methyl_wgcna.R:343`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L343)) — Called on `t(m)` (samples × CpGs, WGCNA's expected orientation — see Section 9). Flags and the code then removes any probe/sample combination with near-zero variance or excess missingness that would otherwise break network construction. Verified present in this module, unlike a documented gap in the *transcriptomics* WGCNA module caught by `tests/testthat/test-wgcna-qc-gate.R` (a regression test for the sibling `mod_wgcna.R`, not this file — see Section 15).
- **`mx_wgcna_celltype_reference()`** (own file, [`mod_methyl_wgcna.R:103-113`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L103-L113)) — Detects cell-type-composition covariate columns generically (≥3 numeric columns whose row sums are ≈1) and picks the one with the highest mean fraction as the implicit reference to leave unticked, avoiding a rank-deficient design. Never hardcodes literal column names, so an uploaded dataset's own differently-named cell-type estimates are handled the same way.

**Processing sequence** (inside `mx_wgcna_filtered <- eventReactive(input$filter_btn, {...})`, [`mod_methyl_wgcna.R:266-358`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L266-L358)):
1. Optionally transpose the raw matrix if `force_transpose` is checked; validate rows ≥ columns afterward.
2. Subset to the chosen sex stratum (`methyl_qc_subgroup_filter()`); require ≥6 samples remaining.
3. Scale check: if not already flagged as M-value, verify the value range is within `[-0.1, 1.1]`; reject with a specific message if it looks like a 0–100 percentage scale or an out-of-range/M-value-like scale. Clip beta values into `(1e-6, 1-1e-6)` to avoid `log2(0/1)` blow-up in the next step.
4. Missingness filter (`methyl_filter_missing()`); require ≥50 CpGs remaining.
5. Logit-transform beta → M-value (`log2(mat/(1-mat))`) unless the input was already M-value scale.
6. Optional covariate residualization: build the design matrix, check rank, fit via `methyl_chunked_lmfit()`, subtract fitted covariate effects from `m` in 20,000-row chunks.
7. Rank by the chosen variability statistic and keep the top N (`mx_wgcna_top_variable()`).
8. `WGCNA::goodSamplesGenes()` gate; require ≥6 samples and ≥20 CpGs remaining after any removal.

**Output data.** A list: `mat` (the final CpG × sample filtered matrix), `stratum_label`, `n_probes_in`/`n_probes_kept`/`n_samples_kept`, `missing_note`, `resid_note`, `resid_covariates`, `var_method`, `top_n`, and `gsg_note`. This is the sole input every later tab depends on, directly or transitively.

---

### 4.2 Tab 2 — Sample QC

**Purpose.** A read-only diagnostic view of the filtered matrix's sample-level structure — no button, no parameters, purely descriptive.

**Input data.** `mx_wgcna_filtered()`'s `$mat` — the same filtered matrix Data & Filtering produced. If that has not been run yet, the tab shows "Build the filtered matrix on 'Data & Filtering' first" ([`mod_methyl_wgcna.R:385`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L385)).

**User inputs.** None. This tab is purely a `reactive`-driven view (technically wrapped in `tryCatch(mx_wgcna_filtered(), error = ...)` rather than a fresh `reactive()`, since it reuses the upstream `eventReactive`'s cached value without adding its own invalidation trigger).

**Functions used.**
- **`stats::hclust(stats::dist(t(f$mat)), method = "average")`** — Average-linkage hierarchical clustering of samples by Euclidean distance, on `t(f$mat)` (samples in rows, so `dist()` computes pairwise sample distances). Output: a dendrogram plotted with base `graphics::plot()`.
- **`stats::prcomp(t(f$mat), scale. = TRUE)`** — Standard PCA with per-CpG scaling, samples as observations (rows). Output: `PC1`/`PC2` scatter via `ggplot2`.
- **`colMeans(is.na(f$mat))`** — Per-sample missingness percentage, plotted as a bar chart.

**Processing.** No transformation beyond what Data & Filtering already produced; this tab only computes distance/PCA/missingness summaries for display.

**Output data.** Three plots: a sample dendrogram, a PC1-vs-PC2 scatter, and a per-sample missingness bar chart. None are stored for downstream reuse — they are diagnostic only.

**Guardrail shown.** If `f$n_samples_kept < 15`, a warning notes that "module detection is not reliable below ~15 samples" ([`mod_methyl_wgcna.R:390`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L390), via `mx_wgcna_guardrails()$low_n_samples`). Note: this same call also computes `low_n_probes` (CpGs < 500) but that flag is **never read or displayed anywhere** — see Section 16 (Audit Table) and Section 18 (Findings).

---

### 4.3 Tab 3 — Soft Threshold

**Purpose.** Selects the soft-thresholding power that will be used to build the weighted network, by evaluating scale-free-topology fit across a range of candidate powers.

**Input data.** `mx_wgcna_filtered()$mat`, transposed to `texpr` (samples × CpGs — see Section 9 for why).

**User inputs** ([`mod_methyl_wgcna.R:430-453`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L430-L453)):

| Input ID | What it is | Default | Range/choices | Used by |
|---|---|---|---|---|
| `network_type` | Signed / signed hybrid / unsigned | Signed | 3 choices | `WGCNA::pickSoftThreshold(networkType=...)` |
| `cor_method` | bicor vs Pearson | Pearson (preloaded) / bicor (uploaded) | 2 choices | `corFnc` argument |
| `power_mode` | Automatic vs manual power selection | Automatic | 2 choices | selects `auto_power` vs `manual_power` |
| `use_custom_powers` | Checkbox: supply an explicit power list instead of a range | ON (preloaded) / OFF (uploaded) | — | switches between `custom_powers` and `power_min/max/step` |
| `r_sq_cutoff` | Target scale-free R² | 0.85 | 0–1, step 0.01 | `RsquaredCut` argument; also the "reached cutoff" comparison |
| `manual_power` | Power to use if `power_mode == "manual"` | 6 | 1–30 | overrides `auto_power` |
| `custom_powers` | Comma-separated explicit power list | `1,2,3,4,5,6,7,8,9,10,12,14,16,18,20` (preloaded) / empty (uploaded) | free text, parsed and validated ≥2 valid values | `powerVector` argument |
| `power_min`/`power_max`/`power_step` | Range-based power list (when custom powers is off) | 1 / 20 / 1 | 1–30 / 2–30 / 1–5 | `seq()` → `powerVector` |
| `power_btn` | Action button | — | — | `eventReactive` trigger for `mx_wgcna_sft()` |

**Functions used.**
- **`WGCNA::pickSoftThreshold()`** ([`mod_methyl_wgcna.R:473-474`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L473-L474)). Input: `texpr` (samples × CpGs), the power vector, network type, correlation function name (`"bicor"` or `"cor"`), and the R² cutoff. What it does: for each candidate power, fits a linear model of log(connectivity) vs. log(connectivity rank) and reports the fit R² (`SFT.R.sq`) and mean/median/max connectivity — the standard scale-free-topology diagnostic. Output: `$fitIndices`, a data.frame with one row per tested power. Why used: identifies the power at which the resulting network best approximates a scale-free topology, the theoretical justification for weighting correlations nonlinearly in WGCNA.
- **`mx_wgcna_top_variable()`/preceding filtering** — already covered; not re-invoked here.
- **Power-selection rule** (own logic, [`mod_methyl_wgcna.R:476-482`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L476-L482)): `auto_power <- if (reached) min(fi$Power[fi$SFT.R.sq >= r_sq_cutoff]) else fi$Power[which.max(fi$SFT.R.sq)]`. This is **not** generic WGCNA-package behavior — it is the exact rule documented in `METHODS_wgcna_sexstratified.md` Section 2.CC.3 ("the power actually used ... is computed directly from each stratum's own `pickSoftThreshold()` fit-index table as the single power ... with the highest observed SFT.R²"), reimplemented faithfully as a first-power-reaching-cutoff-else-best-observed-fit rule, never a hardcoded literal power.

**Processing.** Build the power vector from either the custom list or the min/max/step range → run `pickSoftThreshold()` once → apply the auto/manual selection rule → package the fit table, chosen power, and whether the cutoff was actually reached.

**Output data.**
- `fit_indices`: the full per-power fit table (Power, SFT.R.sq, slope, truncated.R.sq, mean/median/max connectivity) — displayed as a `DT::dataTableOutput` and downloadable as CSV.
- Two `ggplot2` plots: R² vs. power (with the cutoff line and chosen-power marker) and mean connectivity vs. power.
- `power`, `auto_power`, `reached_cutoff`, `max_r_sq`, and the chosen `network_type`/`cor_method`/`r_sq_cutoff` — all carried forward into the Network & Modules tab as the pre-filled `net_power` default and the fixed `network_type`/`cor_method` used in the actual network build.

**Guardrail.** If `max_r_sq < 0.7`, a warning explains that methylation networks may not reach conventional transcriptomic-style fit thresholds and suggests widening the power range or evaluating on connectivity/interpretability instead ([`mod_methyl_wgcna.R:502-503`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L502-L503)) — a direct, honest reflection of the published pipeline's own finding (max R² of 0.756/0.673, both below 0.85, `METHODS_wgcna_sexstratified.md` Section 2.CC.5).

---

### 4.4 Tab 4 — Network & Modules

**Purpose.** Builds the actual weighted co-methylation network and detects modules — the computational core of the submodule.

**Input data.** `mx_wgcna_filtered()$mat` (transposed again to `texpr`) and `mx_wgcna_sft()`'s chosen power/network type/correlation method.

**User inputs** ([`mod_methyl_wgcna.R:550-567`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L550-L567)):

| Input ID | What it is | Default | Range | Maps to `blockwiseModules()` argument |
|---|---|---|---|---|
| `net_power` | Soft-thresholding power (pre-filled from Soft Threshold, editable) | `sft$power` | 1–30 | `power` |
| `tom_type` | TOM variant | Signed | 6 choices (signed/unsigned/Nowick variants ×2) | `TOMType` |
| `max_block_size` | Max probes per computational block | 5,000 | ≥500, step 500 | `maxBlockSize` |
| `min_module_size` | Minimum CpGs per module | 20 (preloaded) / 30 (uploaded) | ≥2 | `minModuleSize` |
| `deep_split` | Dendrogram cut sensitivity | 2 | 0–4 | `deepSplit` |
| `merge_cut_height` | Eigengene-correlation merge threshold | 0.25 | 0–1, step 0.01 | `mergeCutHeight` |
| `pam_stage` | Use PAM refinement | ON | checkbox | `pamStage` |
| `pam_respects_dendro` | PAM respects dendrogram structure | ON | checkbox | `pamRespectsDendro` |
| `reassign_threshold` | Module reassignment p-value threshold | 1e-6 | 0–1 | `reassignThreshold` |
| `min_kme_to_stay` | Minimum kME for a CpG to remain in its module | 0.3 | 0–1, step 0.01 | `minKMEtoStay` |
| `min_core_kme` | Minimum "core" kME | 0.5 | 0–1, step 0.01 | `minCoreKME` |
| `modules_btn` | Action button | — | — | `eventReactive` trigger for `mx_wgcna_net()` |

Every one of the eleven numeric/select defaults above (`min_module_size=20`, `merge_cut_height=0.25`, `deep_split=2`, `max_block_size=5000`, `pam_stage=TRUE`, `pam_respects_dendro=TRUE`, `reassign_threshold=1e-6`, `min_kme_to_stay=0.3`, `min_core_kme=0.5`) matches, verbatim, either the published script's own stated parameters or `WGCNA::blockwiseModules()`'s own package defaults — `METHODS_wgcna_sexstratified.md` Section 2.CC.3 explicitly states "all other arguments at package defaults."

**Functions used.**
- **`WGCNA::blockwiseModules()`** ([`mod_methyl_wgcna.R:591-599`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L591-L599)). Input: `texpr` (samples × CpGs), power, network type, TOM type, correlation type, deep split, min module size, merge cut height, max block size, PAM options, reassignment/kME thresholds, `numericLabels = FALSE` (module colors, not integers), `saveTOMs = FALSE`, `randomSeed = 1234`. What it does: pre-clusters CpGs into blocks (projective k-means) if the probe count exceeds `maxBlockSize`, computes the weighted adjacency and topological overlap matrix within each block, hierarchically clusters `1-TOM`, dynamically cuts the tree into modules, optionally applies PAM-based reassignment of "grey"-adjacent CpGs, then merges modules with highly correlated eigengenes across blocks. Output: `$colors` (named module-color vector), `$MEs` (module eigengenes), `$dendrograms`, `$blockGenes`. Why used: this single function performs network construction, TOM computation, and module detection together — the standard, memory-efficient WGCNA entry point for datasets too large to build one whole-matrix TOM.
- **`get_or_compute_meth_wgcna_blocks()`** (`global.R:477-492`) — Wraps the `blockwiseModules()` call in a content-addressed cache: a key is computed via `digest::digest(key_parts, algo = "xxhash64")` over every input that affects the result (the actual `texpr` matrix plus every network parameter), checked first against an in-session memory cache (`.arthomix_cache`), then against a disk cache under `METH_WGCNA_CACHE_DIR`. Why used: `blockwiseModules()` is one of the most expensive calls in the app; identical inputs+parameters are never recomputed, whether across tab re-visits, button re-clicks with unchanged settings, or (via the disk cache) across separate app sessions.

**Processing.** Transpose the filtered matrix → resolve the correlation type string (`"bicor"` or `"pearson"`) from the Soft Threshold tab's choice → assemble a `key_parts` list of every parameter → call `get_or_compute_meth_wgcna_blocks()`, which either returns a cached result or runs `blockwiseModules()` → attach CpG names to `net$colors` → tabulate module sizes.

**Output data.**
- `module_colors`: a named character vector, one color label per CpG (including "grey" for unassigned CpGs).
- `MEs`: samples × modules eigengene matrix.
- `module_sizes`: a small data.frame of module → CpG count, shown as both a horizontal bar chart and a `DT` table, and downloadable as "Module assignment (CSV)".
- A dendrogram plot (`WGCNA::plotDendroAndColors()`) showing the first block's clustering colored by final module assignment.
- Guardrail messages if all CpGs fell into "grey" (no real modules) or only one real module was found ([`mod_methyl_wgcna.R:624-625`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L624-L625), via `mx_wgcna_guardrails()$all_grey`/`$single_module`).

This is the single most important handoff point in the submodule: every later tab (Module-Trait, Hub CpGs, Results & Export's network summary, Functional Enrichment) depends directly on `mx_wgcna_net()`'s returned list.

---

### 4.5 Tab 5 — Module-Trait Analysis

**Purpose.** Tests whether each module's eigengene associates with a chosen phenotype/trait column.

**Input data.** `mx_wgcna_net()`'s `$MEs` and `$texpr` (for sample IDs), plus `methyl_dataset$sample_sheet` for the trait column itself.

**User inputs** ([`mod_methyl_wgcna.R:674-680`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L674-L680)):

| Input ID | What it is | Default | Range | Used by |
|---|---|---|---|---|
| `trait_col` | Sample-sheet column to correlate against | first of `group`/`Group`/`disease`/`Disease` if present, else first column | any column | `mx_wgcna_encode_trait()` |
| `mt_cor_method` | Pearson vs Spearman | Pearson | 2 choices | correlation function selection |
| `mt_correction` | BH vs Bonferroni | BH | 2 choices | `stats::p.adjust()` method |
| `mt_sig_thr` | FDR significance threshold (for display/counting, not for filtering the table) | 0.05 | 0–1, step 0.01 | significant-module count, heatmap/enrichment gating |
| `traits_btn` | Action button | — | — | `eventReactive` trigger for `mx_wgcna_module_trait()` |

**Functions used.**
- **`WGCNA::orderMEs()`** ([`mod_methyl_wgcna.R:696`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L696)) — Reorders module-eigengene columns into a canonical, related-module-adjacent order (standard WGCNA display convention) after subsetting rows to samples common between the network and the sample sheet.
- **`mx_wgcna_encode_trait()`** (own file, [`mod_methyl_wgcna.R:84-93`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L84-L93)) — Numeric columns pass through unchanged; a column with exactly two non-missing levels is coded 0/1; anything else (e.g. a 3+ category factor) returns `ok = FALSE` with an explicit message that module-trait correlation "needs a numeric column or exactly two levels." This is a real, deliberate **restriction** (Section 1.7), not a silent failure.
- **`WGCNA::cor`** (Pearson path) or **`stats::cor(..., method = "spearman")`** (Spearman path) — Computes the module-eigengene-vs-trait correlation for every module column at once (`corfn(as.matrix(MEs_all), trait_vec, use = "p")`).
- **`WGCNA::corPvalueStudent()`** ([`mod_methyl_wgcna.R:704`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L704)) — Converts a correlation coefficient plus sample size into a two-sided p-value under the standard Student's-t approximation for Pearson correlation significance. Applied uniformly with the single `n_used` value, which is valid here because every module's correlation was computed against the *same* trait vector with the *same* missingness pattern (Section 8.2 discusses this assumption).
- **`stats::p.adjust()`** (BH or Bonferroni) — Applied **only to real (non-grey) modules** ([`mod_methyl_wgcna.R:706-709`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L706-L709)), matching the corrected methodology in `METHODS_wgcna_sexstratified.md` Section 2.CC.3 (which documents and fixes an earlier bug where grey was mistakenly included in the correction).

**Processing.** Match samples between the network and the sample sheet (`methyl_sheet_sample_ids()`) → reorder eigengenes → encode the chosen trait → correlate every module against it in one vectorized call → compute p-values → BH/Bonferroni-correct across real modules only → sort by raw p-value.

**Output data.**
- A results table: `module`, `n_cpgs`, `cor`, `p_value`, `fdr` — downloadable as CSV.
- A heatmap (`mt_heatmap_df()`/`renderPlot`) — for a two-level trait, one column per level with the "other" level's column being the exact negative of the correlation (mathematically exact, since a two-level indicator and its complement are perfectly anti-correlated: `cor(x, 1-v) = -cor(x, v)`); for a numeric trait, a single column labeled by the trait's own name.
- A summary line reporting how many of the real modules reached the significance threshold.

---

### 4.6 Tab 6 — Hub CpGs

**Purpose.** Within one selected module, ranks its member CpGs by how centrally they belong to that module (module membership, kME) and how connected they are within the whole network (intramodular connectivity, kWithin), for follow-up (e.g. candidate biomarker or functional-annotation work).

**Input data.** `mx_wgcna_net()`'s `$texpr`, `$MEs`, `$module_colors`, `$cor_type`, `$network_type`, `$power`; and `methyl_dataset$array_type` for annotation.

**User inputs** ([`mod_methyl_wgcna.R:793-798`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L793-L798)):

| Input ID | What it is | Default | Range | Used by |
|---|---|---|---|---|
| `hub_module` | Which detected (non-grey) module to inspect | (first alphabetically) | any real module color | filters `net$module_colors` |
| `hub_kme_thr` | Minimum absolute kME to be listed as a hub | 0.7 | 0–1, step 0.05 | final filter |
| `hub_top_n` | How many top hub CpGs to return | 20 | 10/20/50/100/Custom | `utils::head()` |
| `hub_top_n_custom` | Free-entry count when "Custom" is chosen | 20 | ≥1 | `utils::head()` |
| `hubs_btn` | Action button | — | — | `eventReactive` trigger for `mx_wgcna_hubs()` |

**Functions used.**
- **`WGCNA::bicor()` or `WGCNA::cor()`** ([`mod_methyl_wgcna.R:830-831`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L830-L831)) — Computes kME (module membership): the correlation of every CpG's methylation profile against the selected module's own eigengene, using whichever correlation function (`bicor`/Pearson) was chosen back on the Soft Threshold tab (`net$cor_type`), for consistency with how the network itself was built.
- **`WGCNA::intramodularConnectivity.fromExpr()`** ([`mod_methyl_wgcna.R:816-822`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L816-L822)) — Computes each CpG's total, within-module, and out-of-module connectivity directly from the expression-style matrix (`net$texpr`), the module assignment, the correlation function, network type, and power — i.e., the same adjacency definition used during network construction. This is wrapped in a **plain `reactive()`, not part of the `hubs_btn` eventReactive**, specifically so it is computed once per network build and reused across every module the user switches to, rather than recomputed on every "Compute hub CpGs" click — the code comment explicitly notes this mirrors an identical optimization already used in the sibling transcriptomics WGCNA module ([`mod_methyl_wgcna.R:806-812`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L806-L812)).
- **`methyl_get_norm_annotation()`** (`normalization.R:252`) — Best-effort CpG→genomic-location/gene annotation join (chr, pos, gene, island_relation, gene_region), available only for array types with a manifest package (450K/EPIC); degrades gracefully (no extra columns, `annotation_available = FALSE`) otherwise.

**Processing.** Compute kME for every CpG in the selected module against that module's eigengene → look up each CpG's pre-computed kWithin → join genomic annotation when available → filter to `|kME| >= hub_kme_thr` → sort descending by `|kME|` → keep the top N.

**Output data.** A table (`cpg`, `kME`, `kWithin`, `abs_kME`, plus annotation columns when available), downloadable as "Hub CpGs (CSV)" with a filename that includes the module name.

---

### 4.7 Tab 7 — Results & Export

**Purpose.** Consolidates a run's parameters into one summary table, provides bulk CSV downloads, and hosts two secondary panels: comparison against the published static reference, and optional Functional Enrichment.

**Input data.** `mx_wgcna_filtered()` and, if computed, `mx_wgcna_net()` and `mx_wgcna_module_trait()`.

**User inputs.**
- `ref_sex` (radio: Female/Male) — which published-reference sex table to display in "Compare with published results," pre-selected to match the current live sex-stratum selection when it maps cleanly to Male ([`mod_methyl_wgcna.R:949-950`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L949-L950)).
- `enrich_btn` — action button triggering `mx_wgcna_enrichment()`.
- Three `downloadButton`s: filtered CpG list, module assignment, analysis parameters.

**Functions used.**
- **Analysis Summary table** ([`mod_methyl_wgcna.R:881-908`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L881-L908)) — Plain R list/data.frame assembly of every parameter used so far (dataset source, stratum, filtering stats, and, once computed, every network parameter and the final module count). No package function beyond `DT::datatable()`.
- **`load_default_wgcna_module_trait()`** (`global.R:505-511`) — Loads the published, offline run's per-sex module–trait CSV directly from disk (`METH_WGCNA_DIR/module_trait_{sex}.csv`), returning `NULL` gracefully if the reference data folder isn't bundled in this deployment.
- **`load_default_dmr_biomarker_panel()`** (`global.R:529-535`) — Loads the published DMP/DMR biomarker panel CSV (`cpg`, `gene`, `dmr_fdr`, `tier` columns, verified against the actual bundled file) for the Functional Enrichment test.
- **`stats::fisher.test()`** ([`mod_methyl_wgcna.R:991-997`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L991-L997)) — Two-sided Fisher's exact test on a 2×2 contingency table (in-module/in-panel × in-module/out-of-panel × out-of-module/in-panel × out-of-module/out-of-panel), one test per significant module. Input: counts of CpGs in each of the four categories, restricted to the CpGs present in the *current live network's own background* (`names(net$module_colors)`), matching the published methodology's stated background ("relative to the full probe background tested in this stage," `METHODS_wgcna_sexstratified.md` Section 2.CC.4). Output: odds ratio and p-value per module.
- **`stats::p.adjust(method = "BH")`** — BH-corrects the enrichment p-values across the tested (already-significant) modules.

**Processing (Functional Enrichment sub-flow, gated behind Module-Trait Analysis having been run):** resolve which published biomarker panel sex to test against from the current sex-stratum selection (see Section 18 for a flagged asymmetry in this mapping) → identify modules with `fdr < mt_sig_thr` and `module != "grey"` → for each, build the 2×2 table against the panel → Fisher's exact test → BH-correct → sort by p-value.

**Output data.**
- Analysis Summary table (parameters), Export section (3 CSV downloads), the published-reference comparison table (static), and — conditionally — the Functional Enrichment results table (module, n_cpgs, n_panel_in_module, n_panel_total, odds_ratio, p_value, fdr), downloadable as CSV.

---

## 5. Complete Function Inventory

### 5.1 Data handling / matrix manipulation

| Function | Defined in | Purpose |
|---|---|---|
| `methyl_qc_subgroup_filter()` | `qc.R:472` | Subsets matrix columns to a sex/group stratum |
| `methyl_sheet_sample_ids()` | `qc.R:456` | Resolves sample-sheet rows to matrix column IDs (by ID column or row order) |
| `methyl_filter_missing()` | `qc.R:30` | Row-wise (CpG) missingness filter |
| `mx_wgcna_top_variable()` | `mod_methyl_wgcna.R:65` | Top-N-by-variability CpG selection |
| `mx_wgcna_row_mads()` | `mod_methyl_wgcna.R:53` | Row-wise MAD, `matrixStats`-accelerated with a base-R fallback |
| `methyl_row_vars()` | `qc.R:22` | Row-wise variance, `matrixStats`-accelerated with a base-R fallback |
| `t()` (base R) | — | Matrix transpose, used at 5 distinct points to switch between CpG-rows and WGCNA's expected sample-rows orientation (Section 9) |

### 5.2 Normalization / preprocessing

| Function | Defined in | Purpose |
|---|---|---|
| `log2(mat/(1-mat))` (inline, no dedicated function) | `mod_methyl_wgcna.R:303` | Beta → M-value logit transform |
| `stats::model.matrix()` | base R | Builds the covariate design matrix for residualization |
| `methyl_chunked_lmfit()` | `mod_methyl_dmp.R:143` | Memory-chunked `limma::lmFit()` wrapper |
| in-place residual subtraction loop | `mod_methyl_wgcna.R:325-332` | Subtracts fitted covariate effects in row chunks |
| `mx_wgcna_celltype_reference()` | `mod_methyl_wgcna.R:103` | Auto-detects the implicit cell-type reference column |

**No quantile normalization, functional normalization, or BMIQ is performed anywhere in this file** — confirmed by inspection; those are Normalization-tab (`mod_methyl_normalization.R`) responsibilities upstream of the shared `methyl_dataset$beta`, not something WGCNA re-applies. See Section 8 for the full normalization audit.

### 5.3 WGCNA package functions

| Function | Purpose | Called at |
|---|---|---|
| `WGCNA::goodSamplesGenes()` | Flags near-zero-variance/excess-missingness probes and samples | [`mod_methyl_wgcna.R:343`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L343) |
| `WGCNA::pickSoftThreshold()` | Scale-free-topology fit across candidate powers | [`mod_methyl_wgcna.R:473`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L473) |
| `WGCNA::blockwiseModules()` | Network construction, TOM, module detection, in one call | [`mod_methyl_wgcna.R:591`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L591) |
| `WGCNA::plotDendroAndColors()` | Dendrogram + module-color bar plot | [`mod_methyl_wgcna.R:638`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L638) |
| `WGCNA::orderMEs()` | Canonical module-eigengene column ordering | [`mod_methyl_wgcna.R:696`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L696) |
| `WGCNA::cor` | Pearson correlation (WGCNA's own NA-robust implementation) | [`mod_methyl_wgcna.R:701, 830`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L701) |
| `WGCNA::bicor` | Biweight midcorrelation (robust to outliers) | [`mod_methyl_wgcna.R:830`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L830) |
| `WGCNA::corPvalueStudent()` | Student's-t p-value for a correlation coefficient | [`mod_methyl_wgcna.R:704`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L704) |
| `WGCNA::intramodularConnectivity.fromExpr()` | Total/within/out-of-module connectivity | [`mod_methyl_wgcna.R:817`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L817) |

Functions the task prompt names as commonly relevant but **not present** in this file: `WGCNA::TOMsimilarity()`, `WGCNA::TOMplot()`, `WGCNA::moduleEigengenes()` (called only implicitly inside `blockwiseModules()`, never as a standalone tab-level call), `WGCNA::mergeCloseModules()` (used only in the *published offline* script's post-hoc consolidation step, Section 1.7 — no live control for it exists here).

### 5.4 Statistical functions

| Function | Purpose |
|---|---|
| `stats::cor(method = "spearman")` | Spearman module-trait correlation alternative |
| `stats::p.adjust(method = "BH"/"bonferroni")` | Multiple-testing correction, restricted to real modules |
| `stats::fisher.test()` | Two-sided exact test for module/biomarker-panel enrichment |
| `stats::prcomp()` | Sample PCA (Sample QC tab) |
| `stats::hclust()` / `stats::dist()` | Sample dendrogram (Sample QC tab) |
| `stats::mad()` (fallback path only) | Row-wise MAD when `matrixStats` unavailable |
| `stats::IQR()` | IQR variability-ranking option |
| `stats::complete.cases()` | Drops samples with any missing covariate before residualization |
| `qr(design)$rank` | Rank-deficiency check on the residualization design matrix |

### 5.5 Visualization functions

| Function | Renders |
|---|---|
| `graphics::plot()` on an `hclust` object | Sample dendrogram (Sample QC) |
| `ggplot2::ggplot() + geom_point()` | PCA scatter (Sample QC), R²-vs-power and connectivity-vs-power (Soft Threshold) |
| `ggplot2::ggplot() + geom_col()` | Missingness bar chart (Sample QC), module-size bar chart (Network & Modules) |
| `WGCNA::plotDendroAndColors()` | Module dendrogram with color bar (Network & Modules) |
| `ggplot2::ggplot() + geom_tile() + geom_text()` | Module-trait correlation heatmap (Module-Trait Analysis) |
| `DT::datatable()` (with `DT::formatRound()`) | Every tabular output across all 7 tabs |

Every plot uses the app-wide `theme_arthomix()` and `ARTHOMIX_COLORS`/`ARTHOMIX_STATUS` palette (`global.R:1417-1438`), consistent with the rest of the application — this module introduces no new visual styling.

### 5.6 Shiny/reactive control functions

| Construct | Role in this module |
|---|---|
| `eventReactive()` (×6: filter, sft, net, module_trait, hubs, enrichment) | Gates each stage strictly behind its own action button; nothing computes speculatively |
| `reactive()` (×1: `mx_wgcna_intramod_conn`) | Deliberately *not* an `eventReactive` — computed once per network, cached and reused across module switches |
| `renderUI()` / `uiOutput()` (per-tab, `tabout_*`) | Each tab body is its own output pair rather than an inline function call — so triggering one stage's computation only invalidates that one Shiny output, not the whole `tabsetPanel` |
| `tryCatch(<eventReactive>(), error = function(e) NULL)` | The idiom used everywhere a downstream tab reads an upstream stage — converts "not run yet"/`validate()` errors into a graceful `NULL` rather than a crash, driving the "not run yet" empty-note UI |
| `validate(need(...))` | Used extensively inside `mx_wgcna_filtered()` and `mx_wgcna_module_trait()` for user-facing, specific rejection messages (Section 15) |
| `req()` | Used in `renderPlot`/`eventReactive` bodies to silently halt when a prerequisite input/reactive is absent |
| `downloadHandler()` (×8) | Soft-threshold CSV, module assignment (×2, one per tab), hub CpGs, filtered CpG list, analysis parameters, module-trait table, enrichment table |
| `withProgress()`/`incProgress()` | Wraps the two most expensive stages (filtering, network construction) with a progress bar |
| `shinyjs::useShinyjs()` | Declared in the UI (`mod_methyl_wgcna_ui`) but not used for any conditional show/hide logic within this file — see Section 18 |

---

## 6. Function-Level Audit Table

| Function | Package/Source | Input | What it does | Output | Why used | Audit observation |
|---|---|---|---|---|---|---|
| `methyl_qc_subgroup_filter()` | App (`qc.R:472`) | Matrix, sample sheet, group column, level | Subsets samples to a sex/group stratum | Filtered matrix + label | Enables sex-stratified analysis matching the published pipeline | No implementation issue identified from the inspected code. |
| `methyl_filter_missing()` | App (`qc.R:30`) | Matrix, max NA fraction | Row-wise missingness filter | Keep vector + note | Standard CpG QC | No implementation issue identified from the inspected code. |
| `mx_wgcna_top_variable()` | App (`mod_methyl_wgcna.R:65`) | Matrix, method, top_n | Ranks and keeps top-N variable CpGs | Filtered matrix + stats | Reduces the probe set to a computationally tractable, information-rich subset | `stat[!is.finite(stat)] <- 0` silently zeroes non-finite variability values rather than excluding them outright; a CpG with, e.g., an `NA` MAD is pushed to the bottom of the ranking (safe) rather than removed, so it can never be spuriously *selected*, but it is retained in the pool being ranked rather than flagged. Minor, and safe as implemented for its actual purpose (ranking, not filtering). |
| `methyl_chunked_lmfit()` | App (`mod_methyl_dmp.R:143`) | Matrix, design matrix, chunk size | Chunked `limma::lmFit()`, bit-for-bit verified against a whole-matrix fit per its own code comment | `lmFit`-shaped object | Avoids >1 full-size matrix copy in memory at 400k+ probe scale | No implementation issue identified from the inspected code (verification claim is stated in the source comment, not independently re-derived in this audit). |
| `WGCNA::goodSamplesGenes()` | `WGCNA` package | `t(m)` (samples × CpGs) | Flags near-zero-variance/excess-missing probes+samples | `$goodGenes`, `$goodSamples` logical vectors | Required pre-network-construction QC gate; a documented gap in the sibling transcriptomics module (caught by `test-wgcna-qc-gate.R`) makes its presence here specifically load-bearing | Correctly implemented: called before `blockwiseModules()`, its flagged rows/columns are actually applied (`m <- m[gsg$goodGenes, gsg$goodSamples]`), and a user-facing note reports how many were removed. |
| `WGCNA::pickSoftThreshold()` | `WGCNA` package | `texpr`, power vector, network type, corFnc, R² cutoff | Fits scale-free-topology statistics per candidate power | Fit-index table | Determines the network's soft-thresholding power | No implementation issue identified; the auto-selection rule built around this call's output (first-reaching-cutoff, else best-observed) is a faithful, non-hardcoded reimplementation of the published methodology. |
| `WGCNA::blockwiseModules()` | `WGCNA` package | `texpr`, 11 network/module parameters | Builds the network, TOM, and modules in one memory-managed call | `$colors`, `$MEs`, `$dendrograms` | Core module-detection engine | `randomSeed = 1234` is hardcoded, not user-exposed — see Findings (Reproducibility, Minor). `saveTOMs = FALSE` means no TOM object persists for a whole-network heatmap view — a deliberate memory-vs-feature trade-off, not a bug. |
| `get_or_compute_meth_wgcna_blocks()` | App (`global.R:477`) | `key_parts` list, compute function | Content-addressed memory+disk cache keyed by `digest::digest(..., "xxhash64")` over every parameter (including the actual matrix) | Cached or freshly computed `blockwiseModules()` result | Avoids re-running the single most expensive computation in the app when inputs are unchanged | No implementation issue identified — the cache key includes the matrix itself, so a different filtered input can never collide with a stale cache entry from a different filter/stratum combination. |
| `mx_wgcna_encode_trait()` | App (`mod_methyl_wgcna.R:84`) | Sample sheet, column name | Encodes a trait as numeric or 0/1 for correlation | Encoded vector or an `ok=FALSE` rejection | Module-trait correlation is a linear-correlation step, not a multi-group test | Correctly implemented and honestly restrictive: explicitly rejects >2-level categorical traits rather than silently mis-encoding them (e.g. via an arbitrary integer factor code, which would produce a meaningless "correlation"). |
| `WGCNA::corPvalueStudent()` | `WGCNA` package | Correlation vector, n | Student's-t p-value for each correlation | p-value vector | Standard WGCNA module-trait significance test | Depends on the assumption that every module's correlation was computed from the same `n_used` samples — true here because all modules are correlated against the *same* trait vector, verified by reading the code; would be invalid if the code ever computed `n_used` per-module while trait missingness varied by module (it does not). |
| `stats::fisher.test()` | base R | 2×2 contingency counts | Exact test of module/biomarker-panel co-occurrence | Odds ratio, p-value | Reproduces the published pipeline's convergent-evidence check between WGCNA modules and the independent DMP/DMR panel | No implementation issue identified; correctly restricted to modules already passing the Module-Trait FDR threshold, avoiding a second full multiple-testing burden over all modules. |
| `mx_wgcna_guardrails()` | App (`mod_methyl_wgcna.R:117`) | n_samples, n_probes, max_r_sq, module_colors | Computes 5 boolean scientific-guardrail flags | Named logical list | "Warn, don't silently force a meaningless analysis through" | `low_n_probes` is computed at one call site (`mod_methyl_wgcna.R:386`) but its value is **never read or displayed anywhere in the UI** — a genuine dead-validation gap (Findings, Minor). |
| `methyl_get_norm_annotation()` | App (`normalization.R:252`) | Array type | Extends base annotation with island/gene-region columns, cached | Annotation data.frame or `ok=FALSE` | Enables genomic interpretation of hub CpGs | Degrades gracefully (returns base annotation without the extra columns) rather than erroring when the extension data can't be loaded — safe as implemented. |
| `WGCNA::intramodularConnectivity.fromExpr()` | `WGCNA` package | `texpr`, module colors, corFnc, network type, power | Per-CpG total/within/out-of-module connectivity | Data.frame incl. `kWithin` | Ranks hub CpGs by network centrality, not just module-eigengene correlation | Wrapped in a plain `reactive()` for reuse across module switches — a deliberate performance choice, correctly implemented (verified: it is not re-triggered by `hubs_btn`, only by a change to `mx_wgcna_net()` itself). |

---

## 7. Normalization and Preprocessing Audit

Answering the 20 audit questions directly against the code:

1. **What data enter WGCNA?** `methyl_dataset$beta` — whatever the Dataset tab (preloaded matrix or user upload) populated. WGCNA itself never touches raw IDAT intensities.
2. **Beta or M-values?** The matrix arrives as beta values (default, `[0,1]`) *or* M-values, per `methyl_dataset$input_scale`. Inside the WGCNA pipeline it is **always converted to M-value** before ranking/residualizing/networking (`m <- if (is_m_scale) mat else log2(mat/(1-mat))`, [`mod_methyl_wgcna.R:303`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L303)) — the network itself is built on M-values, not raw beta values, regardless of the input scale.
3. **Is normalization performed before WGCNA?** **No new normalization is performed inside this file.** WGCNA consumes whatever normalization state `methyl_dataset$beta` is already in when it arrives from the Dataset/Normalization tabs upstream; this module does not call any of `mod_methyl_normalization.R`'s quantile/functional/BMIQ routines.
4. **If yes, which method?** N/A — see #3. The only transformation genuinely performed *within* this file is the beta→M-value logit transform, which is a scale conversion, not a normalization method.
5. **Which function performs it?** The inline expression at [`mod_methyl_wgcna.R:303`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L303); no dedicated helper function.
6. **Is the matrix oriented correctly?** Yes throughout — see Section 9's full trace.
7. **Are samples rows or columns at each stage?** CpG-rows/sample-columns in `f$mat` (the "genomics convention"); transposed to sample-rows/CpG-columns (`texpr`, the "WGCNA convention") at every point that calls a WGCNA network function. Full trace in Section 9.
8. **Are CpGs rows or columns at each stage?** The exact inverse of #7 at each point — see Section 9.
9. **Is the matrix transposed before WGCNA?** Yes, explicitly and repeatedly (`t(m)` at `goodSamplesGenes()`, `t(f$mat)` at Soft Threshold, Network & Modules, and Sample QC's PCA/dendrogram) — see Section 9 for why each transpose is necessary and correctly placed.
10. **Are missing values handled?** Yes: `methyl_filter_missing()` (row-wise, threshold-based) before the M-value transform, and `WGCNA::goodSamplesGenes()` (which also accounts for missingness, not just variance) after residualization/ranking.
11. **Are low-variance CpGs removed?** Yes, but as a **top-N selection**, not a fixed-threshold minimum-variance cut — `mx_wgcna_top_variable()` always keeps exactly `top_n` CpGs (or fewer, if fewer have positive variability), ranked by the user's chosen statistic, rather than excluding CpGs below an absolute variance floor. `qc.R`'s dedicated `methyl_filter_variance()`/`methyl_filter_sd()` (fixed-threshold filters) exist in the codebase but are **not called anywhere in this file** — this module uses only the top-N ranking approach.
12. **Are problematic samples detected?** Only descriptively (Sample QC's dendrogram/PCA/missingness plots) and via `WGCNA::goodSamplesGenes()`'s automatic removal — there is no interactive "exclude this sample" control on this tab (contrast with the QC submodule's own manual-exclusion mechanism, `methyl_apply_manual_exclude()`, which this file does not call).
13. **Is sample clustering performed?** Yes, for diagnostic display only (`stats::hclust()` in Sample QC) — not used to drive any automated outlier removal.
14. **Are outlier samples removed?** Only via `WGCNA::goodSamplesGenes()`'s variance/missingness-based criterion; no distance- or Z-score-based outlier-sample removal step exists.
15. **Is batch correction performed?** **No.** No `sva`, `ComBat`, or similar batch-correction call appears anywhere in this file.
16. **Are biological covariates protected?** Yes, explicitly: the trait/disease-status column is excluded from the residualization covariate candidate list ([`mod_methyl_wgcna.R:165`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L165)), so it can never be accidentally regressed out before the Module-Trait Analysis tab needs it.
17. **Is filtering performed before network construction?** Yes — missingness filter, then (optional) residualization, then top-N variability filter, then `goodSamplesGenes()` — all inside Data & Filtering, strictly upstream of Soft Threshold/Network & Modules.
18. **Could preprocessing unintentionally remove biological signal?** Potentially, in one specific, code-documented way: residualizing against covariates that happen to correlate with the trait of interest (e.g. if a cell-type fraction genuinely differs by disease status) would remove some real trait-associated signal along with the intended confound — a standard, general bias-variance trade-off in confound adjustment, not a coding defect. The published methods document (Section 2.CC.2) discusses an analogous, related trade-off explicitly (favoring known-covariate adjustment over latent-factor/SVA adjustment specifically *because* the latter removes more real co-methylation signal) — this is cited as the scientific rationale for the module's design, not evidence of a bug.
19. **Are preprocessing choices exposed to the user?** Yes, essentially all of them: missingness threshold, transpose override, variability method and top-N count, and the full covariate checkbox list are all live UI controls, not hardcoded.
20. **Are defaults scientifically reasonable?** For the preloaded pathway, the defaults are drawn directly from a documented, published, peer-review-style methods chapter (`METHODS_wgcna_sexstratified.md`) with its own stated rationale for each choice; for the uploaded pathway, the defaults (bicor, signed, min module size 30, no residualization) reflect general WGCNA-on-methylation good practice, though the code offers no dataset-specific tuning guidance beyond the guardrail warnings (Section 5's `mx_wgcna_guardrails()`).

**Summary distinction:** what this code does — a specific residualize-then-MAD-rank-then-`goodSamplesGenes()` pipeline, all user-adjustable, with no batch correction and no outlier-sample removal beyond WGCNA's own variance gate — versus what a "standard" WGCNA-on-methylation workflow might additionally include in some published pipelines (e.g. a Z-score-based sample outlier cut via `WGCNA::adjacency()`+hierarchical clustering height threshold, or ComBat batch correction) is documented above as background context, not implemented here, and the code does not claim otherwise anywhere in its comments or UI text.

---

## 8. Matrix Orientation Audit

```text
methyl_dataset$beta  (as loaded: CpG-rows × sample-columns; parse_upload.R sets
                       rownames(m) <- probe_ids, every other column a sample)
        ↓  [optional t() if force_transpose is checked]
mat_full              CpG-rows × sample-columns  (validated: nrow >= ncol)
        ↓  methyl_qc_subgroup_filter()  -- subsets COLUMNS (samples) to sex stratum
mat (subset)          CpG-rows × sample-columns
        ↓  scale check + clip, then log2(mat/(1-mat)) if beta-scale
m                     CpG-rows × sample-columns   (now M-values)
        ↓  optional: design <- model.matrix(sample-sheet covariates)
        ↓  methyl_chunked_lmfit(m, design); residual subtraction m[rows,] - coefs %*% t(design)
m (residualized)      CpG-rows × sample-columns   (limma convention: features=rows, samples=cols -- matches)
        ↓  mx_wgcna_top_variable() -- ranks/keeps rows (CpGs)
m (top-N)             CpG-rows × sample-columns
        ↓  WGCNA::goodSamplesGenes(t(m))  <-- TRANSPOSE #1: sample-rows × CpG-columns
        ↓  (WGCNA's own convention: rows=samples/observations, columns=genes/probes)
        ↓  m <- m[gsg$goodGenes, gsg$goodSamples]  -- goodGenes indexes m's ROWS (CpGs, = t(m)'s columns);
        ↓                                             goodSamples indexes m's COLUMNS (samples, = t(m)'s rows)
f$mat  (FINAL filtered matrix)     CpG-rows × sample-columns
        ↓  [Sample QC tab]  t(f$mat) -> sample-rows × CpG-columns, for stats::dist()/prcomp() (rows=observations)
        ↓  [Soft Threshold tab]  texpr <- t(f$mat)  <-- TRANSPOSE #2: sample-rows × CpG-columns
        ↓  WGCNA::pickSoftThreshold(texpr, ...)   (WGCNA convention: rows=samples, cols=genes -- matches)
        ↓  [Network & Modules tab]  texpr <- t(f$mat)  <-- TRANSPOSE #3 (recomputed independently)
        ↓  WGCNA::blockwiseModules(texpr, ...)   (WGCNA convention -- matches)
net$colors            named by colnames(texpr) = CpG IDs      -- one label per CpG, correct
net$MEs               samples(rows) × modules(columns)         -- module eigengenes, correct
        ↓  [Module-Trait tab]  MEs_all <- orderMEs(net$MEs[common samples, ])
        ↓  corfn(as.matrix(MEs_all), trait_vec, use="p")  -- correlates each MODULE COLUMN vs. one per-sample trait vector
module_cor            length = number of modules               -- correct: one correlation per module
        ↓  [Hub CpGs tab]  kme <- cor_fnc_r(net$texpr, net$MEs[[me_col]], use="p")
        ↓                  -- correlates EVERY CpG COLUMN of texpr against one per-sample module-eigengene vector
kme                    length = number of CpGs in texpr         -- correct: one kME per CpG
```

**Explicit statements:**
- **Rows = CpGs, columns = samples** at every point *before* a WGCNA network function is called (matches `limma`/general genomics convention, and is what `methyl_chunked_lmfit()`/`limma::lmFit()` requires).
- **Rows = samples, columns = CpGs** at every point *at or after* a WGCNA network function is called (`pickSoftThreshold`, `blockwiseModules`, `goodSamplesGenes`, `intramodularConnectivity.fromExpr`), which is WGCNA's own documented required input convention.
- **Correlations are calculated between CpGs** (kME: one CpG-column of `texpr` at a time against one module-eigengene vector) **and between modules** (module-trait: one module-column of `MEs_all` at a time against one trait vector) — **never between samples**. Verified by reading every `cor()`/`bicor()`/`WGCNA::cor` call site in the file; none correlate two sample-vectors against each other.
- **No orientation bug was identified.** Every transpose is at a WGCNA-API boundary, is locally documented in the surrounding code (explicitly for `goodSamplesGenes()`, implicitly-but-consistently for the network functions), and the post-transpose indexing (`gsg$goodGenes`/`gsg$goodSamples` applied back onto the pre-transpose matrix's rows/columns respectively) is internally consistent.

---

## 9. Soft-Threshold Selection — Explained Against the Code

**What soft-thresholding power means (background).** Raising a correlation matrix to a power `β` before using it as an adjacency matrix exaggerates the difference between strong and weak correlations. A well-chosen `β` produces a network whose node-degree distribution approximates a power law (scale-free) — a property WGCNA's authors argue is both biologically plausible and mathematically useful (a few "hub" nodes, many peripheral ones).

**Which function selects it.** `WGCNA::pickSoftThreshold()`, called once per "Run Soft Threshold Analysis" click ([`mod_methyl_wgcna.R:473-474`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L473-L474)).

**Which candidate powers are tested.** User-controlled: either an explicit comma-separated list (default for the preloaded pathway: `1,2,3,4,5,6,7,8,9,10,12,14,16,18,20`, exactly the published pipeline's own tested range) or a `seq(min, max, step)` range (default for the uploaded pathway: `1:20`).

**What scale-free topology means (background).** A log-log plot of node degree (connectivity) vs. its rank should be approximately linear if the network is scale-free; `SFT.R.sq` is the R² of that linear fit — the code's proxy for "how scale-free is the network at this power."

**How mean connectivity is evaluated.** `pickSoftThreshold()`'s `mean.k.` column is plotted directly (`sft_k_plot`) alongside the R² plot, so a user can see the connectivity-vs-power trade-off (higher power → sparser, lower-connectivity network) alongside the fit statistic.

**Which power is ultimately selected — automatic or user-defined?** Both, by choice: `power_mode` defaults to "Automatic," which applies the rule described above (Section 4.3); switching to "Manual" lets the user override with any integer 1–30 via `manual_power`.

**What happens if no appropriate power is found (i.e. the R² cutoff is never reached)?** The code does **not** error or block progress. `auto_power` falls back to `fi$Power[which.max(fi$SFT.R.sq)]` — the single best-observed-fit power among those actually tested — and the UI displays a specific, non-alarmist explanation: "No strong scale-free topology fit was detected... Methylation networks may not reach conventional transcriptomic-style fit thresholds; consider widening the power range... or evaluating the network on connectivity/interpretability rather than forcing a threshold" ([`mod_methyl_wgcna.R:502-503`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L502-L503)). This exact scenario is not hypothetical: the published reference run for both sexes never reached the 0.85 target at any tested power (max 0.756 female / 0.673 male, both at power 20), so this fallback path is the one the reference analysis itself actually took (`METHODS_wgcna_sexstratified.md` Section 2.CC.5).

---

## 10. Network Construction — Parameter-by-Parameter

All parameters below feed a single `WGCNA::blockwiseModules()` call ([`mod_methyl_wgcna.R:591-599`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L591-L599)).

| Parameter | Actual code value (default) | Meaning | Effect on the analysis |
|---|---|---|---|
| `power` | From Soft Threshold, editable (`net_power`) | Soft-thresholding exponent | Controls how sharply strong vs. weak correlations are separated in the adjacency matrix |
| `networkType` | `"signed"` | Signed network — positive and negative correlations are NOT collapsed together | Prevents hyper- and hypomethylated CpGs moving in opposite directions from being merged into the same module (the WGCNA authors' own documented recommendation, reiterated in `METHODS_wgcna_sexstratified.md` Section 2.CC.3) |
| `corType` | `"bicor"` or `"pearson"`, from the Soft Threshold tab | Correlation measure | bicor down-weights outlier samples (robust); Pearson is the classical measure and matches the published preloaded-pathway default |
| `TOMType` | `"signed"` (6 choices offered) | Topological overlap type | Determines whether shared-neighbor overlap also respects correlation sign |
| `deepSplit` | `2` (0–4) | Dynamic-tree-cut sensitivity | Higher values split the dendrogram into more, smaller modules |
| `minModuleSize` | `20` (preloaded) / `30` (uploaded) | Smallest allowed module | Modules below this size are merged into "grey" (unassigned) |
| `mergeCutHeight` | `0.25` | Eigengene-correlation distance threshold for merging modules | Modules whose eigengenes correlate above `1 - 0.25 = 0.75` are merged into one |
| `maxBlockSize` | `5000` | Max CpGs processed together in one TOM computation | Controls peak memory use; the published pipeline explicitly reduced this from an attempted 20,000 after that size exceeded available RAM (`METHODS_wgcna_sexstratified.md` Section 2.CC.3), and the app's own in-UI note repeats this exact warning ([`mod_methyl_wgcna.R:568`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L568)) |
| `pamStage` | `TRUE` | Enables PAM-based cluster refinement | Reassigns some borderline CpGs from "grey" into a real module based on distance to module centroids |
| `pamRespectsDendro` | `TRUE` | PAM refinement restricted to the same dendrogram branch | Prevents PAM from moving CpGs across unrelated branches |
| `reassignThreshold` | `1e-6` | p-value threshold for reassigning a CpG to a better-fitting module | WGCNA package default; very permissive numerically but standard |
| `minKMEtoStay` | `0.3` | Minimum module membership for a CpG to remain assigned | CpGs weakly correlated with their assigned module's eigengene are dropped to "grey" |
| `minCoreKME` | `0.5` | Minimum kME for a "core" (seed) CpG within a proposed module | Governs whether a candidate module is confirmed at all |
| `numericLabels` | `FALSE` | Module labels reported as colors | Matches the color-coded conventions used throughout the tabs (dendrogram, heatmap, tables) |
| `saveTOMs` | `FALSE` | TOM matrices not written to disk/kept | Reduces memory/disk footprint; trade-off is no whole-network TOM heatmap feature exists |
| `randomSeed` | `1234` (hardcoded) | Seed for the projective k-means preclustering step | Makes repeated runs with identical inputs deterministic within this app, but is not user-configurable and is not necessarily the same seed the original published offline script used |
| `verbose` | `0` | Suppresses WGCNA's own console logging | Cosmetic only |

---

## 11. Module Detection — Step by Step

1. **CpG similarity.** Pairwise correlation (bicor or Pearson) between every pair of CpGs in the filtered set.
2. **Adjacency.** Correlation raised to the chosen `power`, signed per `networkType`.
3. **Topological overlap.** TOM computed from the adjacency matrix (per-block, since the probe count exceeds `maxBlockSize`), which additionally rewards CpGs sharing similar network neighbors, not just direct correlation.
4. **Clustering.** Average-linkage hierarchical clustering of `1 - TOM` within each block.
5. **Dynamic tree cutting.** `deepSplit = 2` controls how aggressively the dendrogram is cut into candidate modules.
6. **Module assignment.** Each CpG receives a color label (`net$colors`), with "grey" reserved for unassigned CpGs.
7. **Module merging.** Candidate modules across blocks whose eigengenes correlate above `1 - mergeCutHeight = 0.75` are merged into one final module — this is the mechanism by which per-block preclustering (necessitated by `maxBlockSize`) is reconciled back into genome-wide modules, with the documented limitation that two CpGs placed in different blocks by the preclustering step can only end up in the same final module if their block-level modules are subsequently merged (`METHODS_wgcna_sexstratified.md` Section 2.CC.3) — a limitation of the underlying WGCNA algorithm, inherited here, not something this app's own code could avoid while still using `blockwiseModules()` at this scale.
8. **Final module count.** Reported directly as `length(setdiff(unique(net$module_colors), "grey"))` ([`mod_methyl_wgcna.R:617`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L617)).

**Biological meaning of the final assignment.** Each non-grey module represents a set of CpGs whose methylation levels rise and fall together across the analyzed samples, independent of any phenotype label — a candidate "co-methylation unit." "Grey" CpGs are explicitly **not** a module; they are the residual bucket of CpGs that did not fit any detected pattern strongly enough, and the code consistently excludes "grey" from every "module count," correlation, correction, and enrichment step (verified across the Module-Trait, Hub CpGs, and Functional Enrichment tabs).

---

## 12. Module–Trait Relationships — Full Detail

- **Which traits are used:** exactly one, user-selected sample-sheet column (`input$trait_col`) per run — the module does not test multiple traits simultaneously.
- **Where the traits come from:** `methyl_dataset$sample_sheet`, the same phenotype table the Dataset tab loaded.
- **How samples are matched:** `methyl_sheet_sample_ids()` resolves sheet rows to the network's sample IDs (by ID column when present, else row order), then `intersect()`s with `rownames(net$texpr)`; requires ≥3 matching samples.
- **Which correlation is calculated:** Pearson (via `WGCNA::cor`) or Spearman (via `stats::cor(method="spearman")`), user's choice — no other correlation type is offered here (contrast with bicor being available for network construction but not for this step).
- **How p-values are calculated:** `WGCNA::corPvalueStudent()`, the standard Student's-t approximation for testing whether an observed correlation differs from zero, applied uniformly with the single matched sample count.
- **Multiple-testing correction:** BH (default) or Bonferroni, applied strictly to real (non-grey) modules only.
- **Heatmap generation:** `mt_heatmap_df()`/`renderPlot` — for a two-level trait, both the trait level and its exact algebraic negative are shown as two heatmap columns; module rows are ordered by correlation; label text size and total plot height both scale with the number of real modules so labels never overlap regardless of module count (an explicit fix noted in the published pipeline's own history, `METHODS_wgcna_sexstratified.md`'s provenance note, and replicated as a first-class feature here rather than merely referenced).
- **Significance thresholds:** the `mt_sig_thr` input (default 0.05) governs the "N of M modules significant" summary count, the heatmap's implicit framing, and — critically — which modules are eligible for the Functional Enrichment tab (only `fdr < mt_sig_thr` modules are tested there).
- **Interpretation of positive/negative associations:** a positive correlation means higher module-eigengene values associate with higher-coded trait values (e.g., disease=1 vs. control=0); a negative correlation means the opposite. The code does not itself state a biological direction (e.g. "hypermethylated in disease") — that interpretation requires knowing which allele of the trait was coded 1 (visible in the table's own `trait_levels`, surfaced in the heatmap's column labels) and is left to the user/reader, appropriately, since the app cannot know the disease-relevant direction for an arbitrary uploaded trait.
- **Module membership (kME) vs. module–trait correlation — the code's own distinction:** kME (Hub CpGs tab) measures how strongly one *CpG* correlates with its own module's eigengene (a within-module concept); module–trait correlation (this tab) measures how strongly one *module's* eigengene correlates with the phenotype (a between-module-and-outcome concept). The two are computed by structurally identical code (`corfn(matrix, vector, use="p")`) but against different pairs of objects — verified by comparing [`mod_methyl_wgcna.R:702`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L702) (module columns × trait vector) against [`mod_methyl_wgcna.R:831`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L831) (CpG columns × one module-eigengene vector).
- **Gene/CpG significance (as a distinct WGCNA concept from kME):** **not implemented as its own named quantity.** WGCNA literature sometimes defines "CpG/gene significance" as the direct correlation of an individual CpG (not the module eigengene) against the trait. This module never computes that value directly; the closest proxy a user could construct is combining the Hub CpGs table's kME (CpG-vs-module) with the Module-Trait table's cor (module-vs-trait) for the same module, but this is an inference the user would have to make manually — the code does not compute or display a direct CpG-vs-trait correlation anywhere.

---

## 13. How the WGCNA Tabs Are Connected

```text
methyl_dataset$beta / $sample_sheet / $input_scale / $array_type / $preloaded
        │  (Dataset tab; shared reactiveValues, read-only from here on)
        ▼
[Tab 1: Data & Filtering]  --(filter_btn)-->  mx_wgcna_filtered()
        │  produces: $mat (CpG × sample), $stratum_label, filtering/residualization notes
        ▼
[Tab 2: Sample QC]  (reads mx_wgcna_filtered() directly; no button, no new reactive stored)
        │
        ▼
[Tab 3: Soft Threshold]  --(power_btn)-->  mx_wgcna_sft()
        │  requires: mx_wgcna_filtered() != NULL, else validate() blocks with a specific message
        │  produces: $power, $network_type, $cor_method, $fit_indices
        ▼
[Tab 4: Network & Modules]  --(modules_btn)-->  mx_wgcna_net()
        │  requires: mx_wgcna_filtered() AND mx_wgcna_sft(), else validate() blocks
        │  produces: $module_colors, $MEs, $texpr, $cor_type, $network_type, $power
        ├──────────────────────────────┬───────────────────────────────┐
        ▼                              ▼                               ▼
[Tab 5: Module-Trait]           [Tab 6: Hub CpGs]              [Tab 7: Results & Export]
--(traits_btn)-->                --(hubs_btn)-->                (Analysis Summary always
mx_wgcna_module_trait()          mx_wgcna_hubs()                 available once mx_wgcna_filtered()
  requires: mx_wgcna_net()         requires: mx_wgcna_net()       exists; Export CSVs require
  + methyl_dataset$sample_sheet    (+ mx_wgcna_intramod_conn(),   mx_wgcna_net(); "Compare with
  produces: $table (module,        a plain reactive() computed    published" needs $preloaded;
  cor, p_value, fdr), $MEs         once per net, reused across    Functional Enrichment needs
        │                          module switches)               mx_wgcna_module_trait())
        ▼                          produces: $table (cpg, kME,
[Tab 7: Functional Enrichment]     kWithin, annotation)
--(enrich_btn)-->
mx_wgcna_enrichment()
  requires: mx_wgcna_module_trait() AND mx_wgcna_net()
  produces: $table (module, odds_ratio, p_value, fdr)
```

**What happens if the upstream analysis has not been run — for every transition:** every downstream tab's `renderUI` wraps its corresponding `eventReactive()` call in `tryCatch(..., error = function(e) NULL)`; the `validate(need(...))` calls inside each `eventReactive` throw a Shiny "validation error" (not a crash) when a prerequisite is missing, which `tryCatch` converts to `NULL`, and the tab then renders a specific, named empty-note (e.g. "Run WGCNA (Network & Modules) before continuing.", [`mod_methyl_wgcna.R:671`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L671)) rather than an error page or a blank screen. This behavior was verified at every one of the 6 `eventReactive` chains listed above.

---

## 14. Reactive Dependency Architecture

- **Which reactive expressions depend on which inputs:** each of the six `eventReactive()`s is triggered *only* by its own action button (`filter_btn`, `power_btn`, `modules_btn`, `traits_btn`, `hubs_btn`, `enrich_btn`). Every other input read inside an `eventReactive`'s body (e.g. `input$max_probe_missing`, `input$var_method`, `input$resid_covariates`) is read via Shiny's standard `isolate()`-under-the-hood `eventReactive` semantics — meaning **changing those inputs alone does not re-trigger computation**; only the next click of the associated button does, at which point the *current* values of every isolated input are used.
- **Whether calculations occur automatically or only after a button click:** exclusively after a button click, for all six computational stages. Sample QC (Tab 2) is the sole exception — it is not gated behind its own button because it performs no new computation, only reads the already-computed `mx_wgcna_filtered()`.
- **Whether outputs depend on previous tabs:** yes, strictly hierarchically, as diagrammed in Section 13 — there is no tab whose computation can run without its documented upstream dependency.
- **Whether cached results are reused:** yes, at two levels — Shiny's own `eventReactive` memoization (the same click's result is reused across every `output$...` that reads it within one reactive flush) and, uniquely for the Network & Modules stage, the explicit content-addressed disk+memory cache in `get_or_compute_meth_wgcna_blocks()` (Section 6).
- **Whether changing an input invalidates downstream results:** **no, not automatically, and this is a deliberate, documented design choice** ([`mod_methyl_wgcna.R:36-40`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L36-L40): "Nothing computes or renders speculatively"). Changing, e.g., `top_n` after already running Network & Modules does **not** clear or recompute `mx_wgcna_net()`; the user must re-click "Build filtered matrix," then "Run Soft Threshold Analysis," then "Run WGCNA," in that order, for the change to propagate.
- **Whether stale results can remain visible:** **yes, by design.** If a user changes an upstream input (e.g. `resid_covariates`) after already viewing Network & Modules results, those results remain on screen, computed from the *old* filtered matrix, until the user manually re-runs every affected downstream stage. No on-screen indicator flags that displayed results are now based on superseded inputs. This is the same idiom the sibling transcriptomics WGCNA module (`mod_wgcna.R`) uses, per this file's own header comment, so it is a consistent, intentional app-wide convention rather than an isolated oversight — but it remains a genuine usability risk (Section 18, Minor Issues).
- **Whether `req()`/validation prevents invalid execution:** yes — every `eventReactive` begins with one or more `validate(need(...))` calls checking its specific prerequisites (upstream stage completed, minimum sample/CpG counts, valid trait encoding, etc.) before doing any real computation.

---

## 15. Input Validation and Error Handling Audit

| Scenario | Status | Where |
|---|---|---|
| No live beta matrix available | **Implemented** | `validate(need(!is.null(methyl_dataset$beta), ...))`, [`mod_methyl_wgcna.R:267`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L267); also gated at the whole-module level (`body_ui`, [`mod_methyl_wgcna.R:1087-1092`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L1087-L1092)) |
| Matrix orientation looks wrong (more columns than rows) | **Implemented** (advisory) | `orientation_check_ui`, [`mod_methyl_wgcna.R:201-209`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L201-L209) — a warning, not a hard block; also a hard block *after* an attempted transpose, [`mod_methyl_wgcna.R:272`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L272) |
| Too few samples in the chosen stratum | **Implemented** | `ncol(mat) >= 6` after stratum filtering, [`mod_methyl_wgcna.R:279`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L279); a softer ≥15 guardrail warning elsewhere (not a hard block) |
| Non-beta-scale values (e.g. 0–100 percentages, out-of-range) | **Implemented**, with distinct messages per likely cause | [`mod_methyl_wgcna.R:284-293`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L284-L293) |
| Too much missingness leaving <50 CpGs | **Implemented** | [`mod_methyl_wgcna.R:301`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L301) |
| Too few samples matching the sample sheet for residualization | **Implemented** (<10 samples) | [`mod_methyl_wgcna.R:312`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L312) |
| Missing covariate values | **Implemented** (<10 complete cases) | [`mod_methyl_wgcna.R:317`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L317) |
| Rank-deficient residualization design | **Implemented** | `qr(design)$rank == ncol(design)` check, [`mod_methyl_wgcna.R:323`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L323) |
| Constant/near-zero-variance CpGs or samples | **Implemented** | `WGCNA::goodSamplesGenes()`, [`mod_methyl_wgcna.R:343-349`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L343-L349) |
| Too few samples/CpGs after all filtering | **Implemented** | `ncol(m) >= 6`, `nrow(m) >= 20`, [`mod_methyl_wgcna.R:350-351`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L350-L351) |
| Duplicate CpGs / duplicate samples | **Missing validation** — no explicit `duplicated()` check anywhere in this file. Relies implicitly on upstream `parse_upload.R`/Dataset-tab behavior, not re-verified here. | — |
| NA/Inf values reaching the correlation/network step | **Partially implemented** — `use = "p"` (pairwise-complete) is passed to every correlation call, and `goodSamplesGenes()` catches excess missingness, but no explicit `is.finite()` sweep exists immediately before `blockwiseModules()` itself | — |
| Missing/invalid trait values for Module-Trait | **Implemented** | `mx_wgcna_encode_trait()` rejection message; `sum(!is.na(trait_vec)) >= 3` check, [`mod_methyl_wgcna.R:700`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L700) |
| Incompatible sample-sheet/matrix dimensions | **Implemented** | `length(common) >= 3` (Module-Trait) / `>= 10` (residualization) checks after `intersect()` |
| Invalid WGCNA parameters (e.g. `min_module_size` set to 0) | **Missing validation** — `numericInput` widgets have a declared `min=` attribute (e.g. `min_module_size` has `min=2`), which constrains the browser UI control but is **not independently re-validated server-side**; a manipulated/edge-case input value could theoretically reach `blockwiseModules()` unchecked | — |
| No modules detected (all grey) | **Implemented** (informational, not blocking) | `mx_wgcna_guardrails()$all_grey`, [`mod_methyl_wgcna.R:624`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L624) |
| Only one module detected | **Implemented** (informational) | `mx_wgcna_guardrails()$single_module`, [`mod_methyl_wgcna.R:625`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L625) |
| Low sample count (<15) | **Implemented** (informational, two separate call sites) | Sample QC and Data & Filtering summary |
| Low CpG count (<500) | **Computed but not surfaced** — see Section 6/18 | `mx_wgcna_guardrails()`'s `low_n_probes` field |
| Poor scale-free fit (<0.7 R²) | **Implemented** (informational) | [`mod_methyl_wgcna.R:502-503`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L502-L503) |

---

## 16. Scientific and Statistical Audit

### Data assumptions
- **Methylation representation.** M-values (logit-transformed beta) are used for every quantitative step (residualization, variability ranking, correlation) — a well-established choice for linear-model-style analysis of methylation data, since M-values are approximately homoscedastic across their range whereas beta values compress toward 0 and 1. This is a reasonable, standard assumption, matching the same rationale documented in this codebase's DMP submodule (`methyl_chunked_lmfit()`'s own use in `mod_methyl_dmp.R`).
- **Sample size.** The hard minimum enforced is 6 samples (nearly unusable for a real network); the *warned* minimum is 15 (`mx_wgcna_guardrails()$low_n_samples`). WGCNA methodological guidance typically recommends considerably more (often ≥20-30) for stable module detection — the app's 15-sample threshold is a pragmatic UI guardrail, not a claim that 15 samples is statistically ideal, and the warning text itself says only that results become "not reliable," not that they are invalid.
- **Feature count.** A 500-CpG warned minimum is computed but never shown (Section 18) — a scientific assumption (WGCNA needs "a broad, correlation-structure-rich" input, `METHODS_wgcna_sexstratified.md` Section 2.CC.2) that currently has no working UI enforcement.
- **Missingness.** Handled at both the CpG level (`methyl_filter_missing()`) and, redundantly but safely, again by `goodSamplesGenes()`.
- **Variance.** CpGs are always selected by top-N variability, never filtered by an absolute variance floor within this file — meaning a dataset where even the "most variable" CpGs are nearly constant (e.g. all near 0 or 1) would still produce a full-size `top_n` matrix, just one WGCNA might then struggle to find real structure in; `goodSamplesGenes()` is the only downstream safety net for that scenario.

### Network assumptions
- **Correlation choice** (bicor vs Pearson) is user-exposed and defaults sensibly per pathway (Pearson to match the published methodology when reproducing it; bicor, a more outlier-robust default, for a fresh/unknown uploaded dataset).
- **Network type** defaults to signed, matching the WGCNA authors' own documented recommendation for exactly the reason stated in this file's UI text and the published methods document (avoiding conflation of hyper-/hypomethylation directions).
- **Scale-free topology.** The code neither assumes nor asserts that a good fit will be found — it explicitly handles and explains the case where it is not (Section 9), which is scientifically honest given that DNA methylation networks are not guaranteed, and in this deployment's own published reference run did not, reach the classical 0.85 R² target.
- **Soft threshold and module size** are both fully user-adjustable with scientifically grounded defaults (Section 10).

### Statistical assumptions
- **Correlation assumptions.** Pearson (and `corPvalueStudent()`'s significance test) assumes an approximately linear, bivariate-normal relationship; Spearman is offered as a rank-based alternative for the module-trait step specifically, but **not** for network construction itself (bicor is the offered robust alternative there instead — a different robustness strategy for a different purpose, and a defensible design, not an inconsistency).
- **Trait encoding.** Restricted to numeric or exactly-two-level categorical traits (Section 4.5/12) — a genuine, stated methodological limitation of a linear-correlation-based module-trait test, not an oversight; it is the correct restriction given the method (a multi-category ANOVA-style test would be a different, unimplemented feature).
- **Multiple testing.** BH/Bonferroni, correctly restricted to real (non-grey) modules for the module-trait step, and again correctly restricted to only the *already-significant* subset of modules for the downstream enrichment step (avoiding compounding the correction unnecessarily across the full module set a second time).
- **P-value interpretation.** `corPvalueStudent()`'s p-value tests whether a correlation differs from zero; it says nothing about effect size or biological importance on its own, and the app displays the correlation coefficient itself alongside every p-value/FDR, allowing a reader to judge both together.

### Biological interpretation
- **Are modules appropriately interpreted as co-methylated CpGs?** Yes — nowhere in the UI or code does the module claim a module represents anything beyond a correlation-derived grouping; "co-methylation module" is the term used consistently in the tab's own title and description.
- **Do module–trait associations establish association rather than causation?** The code computes and displays only correlation coefficients, p-values, and FDR — no causal-inference method (e.g. Mendelian randomization, which *is* implemented as its own, separate submodule in this application) is invoked here. The UI text itself never uses causal language ("causes," "leads to") about module-trait results; it uses "correlated," "associated." This is scientifically appropriate framing for what is implemented.
- **Does the implementation support the biological conclusions claimed by the UI?** Based on the text actually present in the module (status messages, tab descriptions, guardrail warnings), yes — every claim inspected is either a direct restatement of a computed statistic or an honest caveat about a limitation (poor fit, low sample size, restricted trait types). No instance of unsupported or overstated claims (e.g. "this proves," "this guarantees") was found in the UI strings read during this audit.

---

## 17. Findings

### Critical Issues

None identified. No finding in this audit was assessed as substantially invalidating results or breaking execution under normal use.

### Major Issues

**Finding: `low_n_probes` guardrail is computed but never displayed.**
- **Location:** [`mod_methyl_wgcna.R:120`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L120) (definition), [`mod_methyl_wgcna.R:386`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L386) (only call site passing `n_probes`)
- **Relevant function:** `mx_wgcna_guardrails()`
- **What the code does:** Computes `low_n_probes = n_probes < 500` inside the Sample QC tab's guardrail call, but the returned list's `$low_n_probes` field is never read by any `if`, `p()`, or conditional rendering anywhere in the file (verified by an exhaustive grep for `low_n_probes` and `gr$` across the whole file — the only match for the former is its own definition line).
- **Why it matters:** A user running WGCNA on a very small, aggressively-filtered CpG set (e.g. `top_n` reduced far below 500, or a naturally CpG-poor uploaded dataset) receives no warning that the network is being built on a probe count the code's own author judged too small to be a "broad, correlation-structure-rich" input (the phrase used in the published methodology this module reproduces).
- **Potential consequence:** A user could interpret a small-probe-count, potentially unstable or uninformative module structure as equally trustworthy as a full-scale run, with no in-app signal to prompt caution.
- **Recommended correction (not applied — audit only):** Surface `gr$low_n_probes` as a warning message, analogous to the existing `low_n_samples` warning at the same call site.

*Classified as Major rather than Critical because it is a missing warning, not an incorrect computation — the module still runs and produces internally consistent output regardless of CpG count; only user awareness of a real scientific caveat is affected.*

### Minor Issues

**Finding: Stale downstream results can remain visible after an upstream input change.**
- **Location:** Architectural — applies to the whole `eventReactive` chain, e.g. [`mod_methyl_wgcna.R:266, 460, 576, 687, 824, 979`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L266)
- **Relevant function:** every `eventReactive()` in the file
- **What the code does:** Each stage recomputes only when its own action button is clicked; changing an upstream input alone does not invalidate or flag already-displayed downstream results.
- **Why it matters:** This is an explicit, documented design choice (header comment, [`mod_methyl_wgcna.R:36-40`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L36-L40)) shared with the sibling transcriptomics WGCNA module, made specifically to avoid speculative computation of an expensive pipeline. It is not a defect in the ordinary sense.
- **Potential consequence:** A user who changes, e.g., `top_n` or `resid_covariates`, then navigates directly to Network & Modules without re-clicking every intervening button, sees results computed from the *previous* filtered matrix, with no on-screen indication that inputs have since changed.
- **Recommended correction (not applied — audit only):** An optional, purely cosmetic "inputs have changed since this was last run" banner, without altering the deliberate no-auto-recompute behavior.

**Finding: `randomSeed = 1234` is hardcoded and not user-exposed.**
- **Location:** [`mod_methyl_wgcna.R:597`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L597)
- **Relevant function:** `WGCNA::blockwiseModules()` call
- **What the code does:** Fixes the random seed governing `blockwiseModules()`'s internal projective k-means preclustering step to the literal value 1234, for every run, regardless of dataset or parameters.
- **Why it matters:** Makes repeated runs of *this app* with identical inputs reproducible (a positive), but the value is not the same seed the original published offline script used (that seed is not stated in `METHODS_wgcna_sexstratified.md` and is not recoverable from this repository), so exact numeric reproduction of the published reference tables via the live tool is not guaranteed even with matching UI settings — this is already correctly caveated in Section 3.1 of this document and is not claimed otherwise anywhere in the app's own UI text.
- **Potential consequence:** Minor — affects exact-reproduction expectations only, not correctness of the live analysis on its own terms.
- **Recommended correction (not applied — audit only):** Expose the seed as an advanced, optional input if exact multi-run comparability across different sessions/users becomes a requirement.

**Finding: No explicit duplicate-CpG/duplicate-sample check within this file.**
- **Location:** Whole-file — no `duplicated()` or `anyDuplicated()` call appears anywhere in `mod_methyl_wgcna.R`.
- **Relevant function:** N/A (absence of a function)
- **What the code does:** Relies entirely on upstream data-loading code (`parse_upload.R`, the Dataset tab) to have already produced a matrix with unique row/column names; performs no defensive re-check itself.
- **Why it matters:** If a duplicate CpG ID or duplicate sample ID were somehow present in `methyl_dataset$beta` when it reaches this module, several R indexing operations (`match()`, `[` by name) used throughout this file would silently resolve to the *first* matching row/column rather than erroring, which could misassign a sample's phenotype or a CpG's annotation.
- **Potential consequence:** Silent misassignment in an edge case that this audit did not find any current code path capable of introducing (upstream parsing was not the subject of this audit and is documented separately for the Dataset tab).
- **Recommended correction (not applied — audit only):** A defensive `stopifnot`/`validate` check if this scenario is judged plausible given how the Dataset tab actually parses uploads.

**Finding: The "All samples combined" sex stratum defaults Functional Enrichment to the female biomarker panel.**
- **Location:** [`mod_methyl_wgcna.R:983-984`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L983-L984)
- **Relevant function:** `mx_wgcna_enrichment()`
- **What the code does:** `panel_sex <- if (identical(sex_choice, unname(sex_choices_r()["Male"]))) "male" else "female"` — the *only* two branches are "male" and "female"; if the current `sex_stratum` is `"__all__"` (All samples combined), the `else` branch fires and the female panel is loaded, with no combined/sex-agnostic biomarker panel available to load instead (none exists in the bundled reference data).
- **Why it matters:** A user who runs the "All samples combined" pathway (explicitly noted elsewhere in this same file as "not part of the published methodology," [`mod_methyl_wgcna.R:1071`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L1071)) and then runs Functional Enrichment gets no explicit statement that their combined-sex modules are being tested against a female-only reference panel specifically (as opposed to, say, an average of both, or an error stating no appropriate panel exists).
- **Potential consequence:** A user could misinterpret a combined-sex enrichment result as sex-agnostic when it is actually benchmarked against a sex-specific reference panel chosen by a silent default rather than an explicit choice.
- **Recommended correction (not applied — audit only):** Either surface the resolved `panel_sex` more prominently in the enrichment result header (it is already shown — "Tested against the %s biomarker panel," [`mod_methyl_wgcna.R:1008`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L1008) — so the information IS present, just not flagged as a consequential default choice), or disable/caveat the enrichment button specifically for the "All samples" stratum.

*Downgraded to Minor on reflection: the resolved panel sex IS shown to the user in the result panel itself, just not called out as a non-obvious default at selection time — this is a UX clarity gap, not a hidden or silent miscalculation.*

**Finding: No unit/regression test exists specifically for `mod_methyl_wgcna.R`.**
- **Location:** N/A — absence confirmed via `find tests -iname "*wgcna*"`, which returns only `tests/testthat/test-wgcna-qc-gate.R`, a regression test for the **sibling transcriptomics** module `R/transcriptomics/mod_wgcna.R`, not this file.
- **Relevant function:** N/A
- **What the code does:** N/A
- **Why it matters:** The transcriptomics WGCNA module's own `goodSamplesGenes()` gap (that test's stated motivation) was caught by comparison against this methylomics module's already-correct behavior — but this methylomics module itself has no equivalent automated guard against a future regression of its own `goodSamplesGenes()` call, or any of its other validated behaviors (rank-deficiency check, trait-encoding rejection, guardrail thresholds).
- **Potential consequence:** A future refactor of this file could silently remove or weaken one of the validated behaviors documented in Section 15 without any automated test failing.
- **Recommended correction (not applied — audit only):** A dedicated `test-methyl-wgcna-*.R` suite mirroring the transcriptomics module's own regression-test pattern.

### No Issue / Correctly Implemented

- Matrix orientation handling across all 5 transpose points (Section 8) — no bug identified.
- The soft-threshold auto-selection rule (first-power-reaching-cutoff, else best-observed) — a faithful, non-hardcoded reimplementation of the published methodology (Section 9).
- Module-trait BH/Bonferroni correction restricted to real modules — correctly reproduces a fix the published pipeline's own methods document says was needed and applied (Section 12, `METHODS_wgcna_sexstratified.md` Section 2.CC.3's footnote).
- The content-addressed `get_or_compute_meth_wgcna_blocks()` cache — key correctly includes the actual input matrix, not just parameters, so no stale-cache collision risk across different filtered inputs (Section 6).
- `WGCNA::goodSamplesGenes()` is present and its output is actually applied to the matrix — this exact gap was found *missing* in the sibling transcriptomics module by a dedicated regression test, making its correct presence here specifically verified as load-bearing, not incidental (Section 4.1, Section 15).
- Trait encoding's explicit rejection of >2-level categorical columns, rather than a silent, statistically meaningless numeric coercion (Section 12).
- The Functional Enrichment step's background universe (the live network's own CpG set) and correction scope (already-significant modules only) both correctly match the published methodology's stated approach (Section 4.7, `METHODS_wgcna_sexstratified.md` Section 2.CC.4).
- No causal-language overreach was found anywhere in the module's UI text (Section 16).

---

## 18. Reproducibility Audit

**R packages used (directly, by this file):** `WGCNA`, `stats`, `limma` (via `methyl_chunked_lmfit()`), `matrixStats` (optional, with a base-R fallback), `digest` (via `get_or_compute_meth_wgcna_blocks()`), `data.table` (via the reference-table loaders), `DT`, `ggplot2`, `shiny`, `shinyjs` (declared but not functionally used within this file — no `shinyjs::show`/`hide`/`toggle` call appears anywhere in `mod_methyl_wgcna.R`, so `useShinyjs()` at [`mod_methyl_wgcna.R:130`](../ArthOMix/R/methylomics/mod_methyl_wgcna.R#L130) is currently a no-op dependency declaration within this specific file — worth noting as a minor code-cleanliness observation, not a functional issue since the app-wide `shinyjs` initialization elsewhere is unaffected).

**Package functions:** enumerated fully in Sections 5 and 6.

**Parameters and defaults:** enumerated fully in Sections 4 (per tab) and 10 (network construction specifically).

**Random seeds:** exactly one, `randomSeed = 1234`, hardcoded (Section 18 Findings, Minor).

**User-configurable parameters:** the overwhelming majority of scientifically consequential choices are exposed as UI controls — missingness threshold, variability method/count, residualization covariates, network type, correlation method, power selection mode and candidates, R² cutoff, TOM type, block size, module size, deep split, merge cut height, PAM options, reassignment/kME thresholds, trait column, correlation method for module-trait, multiple-testing method, significance threshold, hub module/threshold/count.

**Hard-coded thresholds not exposed to the UI:** the `randomSeed`; the 6-sample and 20-CpG hard minimums inside `mx_wgcna_filtered()`; the 10-complete-case minimum for residualization; the 3-sample minimum for module-trait matching; the `reassignThreshold`/`minKMEtoStay`/`minCoreKME` *are* exposed (Section 10), contrary to what might be assumed — they are not hidden defaults.

**Data dependencies:** the live pipeline needs only `methyl_dataset$beta` (and, optionally, `$sample_sheet` for residualization/module-trait/enrichment); the "Compare with published results" and Functional Enrichment panels additionally need the bundled `data/preloaded/methylomics/tables/script05_wgcna_sexstratified/` and `script04_dmr_sexstratified/` reference folders, gracefully degrading (via `load_default_*()`'s `NULL`-return contract) when absent.

**Required metadata:** none is strictly required for Data & Filtering through Network & Modules to run (a sample sheet is optional); a sample sheet becomes required only for residualization, Module-Trait Analysis, and (indirectly, for the published-panel comparison to be meaningful) Functional Enrichment.

**Downloadable results:** 8 distinct CSV exports (Section 5.6), collectively sufficient to reconstruct every table shown in the UI outside of the app itself.

**Hidden assumptions:** the single most consequential one is the assumption, embedded in the "first-power-reaching-cutoff, else best-observed" rule, that the *tested* power range is wide enough to contain a meaningful answer — the code does not itself widen the range automatically if the best-observed fit is still poor (Section 9), and a user must notice the guardrail warning and manually adjust the range.

**Is the analysis reproducible from the UI and code?** For the **live pathway**: yes, deterministically, given the same `methyl_dataset$beta`/`$sample_sheet` and the same UI parameter values (the hardcoded `randomSeed` makes `blockwiseModules()`'s own preclustering step deterministic within this app). For **exact numeric reproduction of the published, offline reference tables**: not guaranteed (Section 3.1, Section 18 Findings) — the app is explicit that the "Compare with published results" panel shows a separate, static, previously-computed result, not a live recomputation target.

---

## 19. End-to-End WGCNA Pipeline

| Stage | Input | Function | Processing | Output | Used by |
|---|---|---|---|---|---|
| Data source | Dataset tab selection | `mod_methyl_dataset_server()` (upstream module) | Loads preloaded or uploaded matrix + sheet | `methyl_dataset$beta`, `$sample_sheet` | Every stage below |
| Sex stratification | `methyl_dataset$beta`, `$sample_sheet`, `sex_stratum` input | `methyl_qc_subgroup_filter()` | Subsets samples (columns) to chosen stratum | Stratum-subset matrix + label | Missingness filter |
| Missingness filter | Stratum-subset matrix | `methyl_filter_missing()` | Row-wise NA-fraction threshold | Reduced CpG set | M-value transform |
| Scale conversion | Reduced CpG set | inline `log2(mat/(1-mat))` | Beta → M-value (unless already M-scale) | M-value matrix | Residualization / variability ranking |
| Residualization (optional) | M-value matrix, `resid_covariates` | `methyl_chunked_lmfit()`, in-place subtraction | Regresses out selected covariates | Residualized M-value matrix | Variability ranking |
| Variability ranking | (Residualized) M-value matrix | `mx_wgcna_top_variable()` | Ranks CpGs by MAD/variance/SD/IQR, keeps top N | Top-N CpG matrix | `goodSamplesGenes()` gate |
| Sample/probe QC gate | Top-N CpG matrix (transposed) | `WGCNA::goodSamplesGenes()` | Flags/removes near-zero-variance or excess-missing rows/columns | **Final filtered matrix `f$mat`** | Every subsequent tab |
| Soft-threshold scan | `t(f$mat)` | `WGCNA::pickSoftThreshold()` | Fits scale-free-topology R²/connectivity per candidate power | Fit table, chosen power | Network construction |
| Network + module detection | `t(f$mat)`, chosen power/parameters | `WGCNA::blockwiseModules()` (cached via `get_or_compute_meth_wgcna_blocks()`) | Adjacency → TOM → clustering → dynamic cut → merge | `$colors`, `$MEs`, `$dendrograms` | Module-Trait, Hub CpGs, Results & Export |
| Module–trait correlation | `$MEs`, sample sheet trait column | `WGCNA::orderMEs()`, `WGCNA::cor`/`stats::cor`, `WGCNA::corPvalueStudent()`, `stats::p.adjust()` | Correlates each module eigengene against the encoded trait; BH/Bonferroni-corrects | Module-trait table, heatmap | Functional Enrichment; Results & Export |
| Hub CpG ranking | `$texpr`, `$MEs`, chosen module | `WGCNA::bicor`/`WGCNA::cor`, `WGCNA::intramodularConnectivity.fromExpr()`, `methyl_get_norm_annotation()` | Computes kME + kWithin, joins annotation, filters/ranks | Hub CpG table | Standalone (no further in-app consumer) |
| Functional enrichment | Significant modules, DMP/DMR biomarker panel | `stats::fisher.test()`, `stats::p.adjust()` | 2×2 exact test per significant module vs. panel | Enrichment table | Standalone (Results & Export display) |
| Export/summary | All of the above | `DT::datatable()`, `downloadHandler()` | Aggregates parameters and provides CSV downloads | Summary table, 8 CSV files | End user (no in-app consumer) |

---

## 20. Educational Explanation of Major Concepts

### Soft-thresholding power
- **What is it?** An exponent applied to correlation values before treating them as network-edge weights.
- **Why is it needed?** Raw correlations produce a network where almost every pair of features has some nonzero weight; raising to a power suppresses weak correlations relative to strong ones, targeting an approximately scale-free network topology.
- **What does this application do?** Tests a user-chosen set of candidate powers with `WGCNA::pickSoftThreshold()` and selects one automatically (first to reach the target fit, else the single best-observed) or lets the user override manually.
- **What goes into it?** The filtered, sample-rows × CpG-columns matrix, and a list of candidate integer powers.
- **What comes out?** A fit-index table (one row per power) and one chosen power.
- **How should the result be interpreted?** A high R² at the chosen power supports treating the resulting network as approximately scale-free; a low R² (as in this deployment's own published reference run) means that assumption is only weakly supported, and results should be read with correspondingly more caution about network-level claims, though module detection can still proceed.
- **What should the user be careful about?** Not over-trusting a chosen power just because it is "the recommended one" when the underlying fit was poor — the app itself warns about this, but the user must read the warning.

### Topological overlap and module detection
- **What is it?** A measure of similarity between two nodes based not just on their direct correlation but on how similar their sets of network neighbors are; used to build a dissimilarity measure (`1 - TOM`) for hierarchical clustering.
- **Why is it needed?** Direct correlation alone can be noisy for individual pairs; TOM is more robust because it aggregates information across many neighboring relationships.
- **What does this application do?** Computes TOM internally as part of the single `WGCNA::blockwiseModules()` call, per computational block, then hierarchically clusters and dynamically cuts the resulting tree into modules.
- **What goes into it?** The adjacency matrix (correlation raised to the chosen power).
- **What comes out?** A module-color label per CpG, plus per-block dendrograms.
- **How should the result be interpreted?** Each non-grey color is one candidate co-methylation module; "grey" specifically means "did not fit any detected module," not a module of its own.
- **What should the user be careful about?** Very large or very small modules deserve extra scrutiny (the app already warns if everything collapses to one module or all-grey); module composition near block boundaries is influenced by the preclustering step's own partition, a known limitation inherited from `blockwiseModules()`'s block-based algorithm at this probe-count scale.

### Module eigengenes and module-trait correlation
- **What is it?** A module eigengene is one number per sample summarizing that module's overall methylation level (its first principal component); module-trait correlation tests whether that one-number-per-sample summary associates with a phenotype.
- **Why is it needed?** Testing one summary number per module, instead of hundreds/thousands of individual CpGs, dramatically reduces the multiple-testing burden while still capturing coordinated, module-level signal.
- **What does this application do?** Computes eigengenes automatically inside `blockwiseModules()`, then lets the user pick exactly one trait column and correlation method to test every module's eigengene against it, with BH/Bonferroni correction restricted to real modules.
- **What goes into it?** The eigengene matrix and one encoded trait vector.
- **What comes out?** One correlation, p-value, and FDR per real module.
- **How should the result be interpreted?** A significant, strongly positive or negative correlation suggests that module's coordinated methylation pattern tracks the chosen trait in this cohort; it is an association, not evidence of causation, and the app's own text never claims otherwise.
- **What should the user be careful about?** The trait-encoding restriction (numeric or exactly two levels only); the direction of a correlation depends on how the trait's two levels happened to be alphabetically/factor-ordered, which the heatmap's own column labels make explicit but a reader must still check.

### Hub CpGs
- **What is it?** CpGs that are both strongly correlated with their own module's eigengene (high kME, i.e. "module membership") and highly connected within the network overall (high kWithin, "intramodular connectivity").
- **Why is it needed?** Not every CpG in a module contributes equally; hub CpGs are the most representative/central members, useful candidates for follow-up validation or mechanistic interpretation.
- **What does this application do?** Computes kME by correlating every CpG against the selected module's eigengene, and kWithin via `WGCNA::intramodularConnectivity.fromExpr()`, then filters and ranks by a user-chosen |kME| threshold and top-N count.
- **What goes into it?** The full network's expression-style matrix, the module-eigengene matrix, and the chosen module.
- **What comes out?** A ranked table of CpGs with kME, kWithin, and (when available) genomic annotation.
- **How should the result be interpreted?** A high-kME CpG is one whose own methylation pattern closely tracks its module's overall pattern — a good representative, not necessarily one individually significant for any trait (that is a separate, unimplemented "CpG significance" concept, Section 12).
- **What should the user be careful about?** Genomic annotation (gene, chromosome, position) is only available for array types with a bundled manifest (450K/EPIC); it is silently absent (with a stated reason) for other array types or unannotated uploads.

---

## 21. Code-to-Concept Mapping

| Code | Scientific concept | Practical role in this module |
|---|---|---|
| `methyl_qc_subgroup_filter()` | Stratified analysis design | Splits samples by sex before any downstream computation, matching the published pipeline's per-sex-first methodology |
| `mx_wgcna_top_variable()` | Feature/probe selection for network construction | Reduces the genome-wide CpG set to a tractable, information-rich subset |
| `methyl_chunked_lmfit()` | Covariate adjustment (confound removal) | Removes known-covariate effects (age, smoking, cell type) before ranking CpGs by residual variability, per the published pipeline's explicit design rationale (Section 2.CC.2 of `METHODS_wgcna_sexstratified.md`) |
| `WGCNA::goodSamplesGenes()` | Network-readiness QC | Guarantees no zero-variance/excess-missing probe or sample reaches network construction |
| `WGCNA::pickSoftThreshold()` | Soft-thresholding power selection | Chooses the exponent producing the best available approximation to scale-free topology |
| `WGCNA::blockwiseModules()` | Network construction + module detection | Builds the weighted co-methylation network and identifies co-methylation modules, in one memory-managed call |
| `net$MEs` | Module eigengenes | One per-sample summary value per module, the unit correlated against phenotype |
| `WGCNA::orderMEs()` / `WGCNA::cor` (module-trait step) | Module-trait association | Links module-level co-methylation patterns to a chosen clinical/phenotypic trait |
| `WGCNA::corPvalueStudent()` + `stats::p.adjust()` | Statistical significance of module-trait association | Determines which modules' trait association exceeds chance, after multiple-testing correction |
| `WGCNA::bicor`/`WGCNA::cor` (hub CpG step) | Module membership (kME) | Quantifies how representative a CpG is of its assigned module |
| `WGCNA::intramodularConnectivity.fromExpr()` | Intramodular connectivity (kWithin) | Quantifies a CpG's network centrality within its module |
| `stats::fisher.test()` | Convergent-evidence enrichment testing | Tests whether a significant module's CpGs overlap, more than chance, with an independently derived biomarker panel |
| `get_or_compute_meth_wgcna_blocks()` | Reproducible, memoized computation | Ensures identical inputs+parameters never re-run the expensive network-construction step |

---

## 22. Thesis-Oriented Tab Summary

| Tab | Input | Main Function | Processing | Output | Scientific Purpose |
|---|---|---|---|---|---|
| Data & Filtering | `methyl_dataset$beta`/`$sample_sheet`, sex stratum, filter/residualization settings | `methyl_qc_subgroup_filter()`, `methyl_filter_missing()`, `methyl_chunked_lmfit()`, `mx_wgcna_top_variable()`, `WGCNA::goodSamplesGenes()` | Stratify → scale-check → missingness-filter → M-value transform → optional residualization → top-N variability ranking → QC gate | Final filtered CpG × sample matrix | Prepares a tractable, confound-aware, network-ready CpG set |
| Sample QC | Filtered matrix | `stats::hclust()`, `stats::prcomp()` | Sample dendrogram, PCA, per-sample missingness | Diagnostic plots | Visual check of sample-level structure before network construction |
| Soft Threshold | Filtered matrix (transposed) | `WGCNA::pickSoftThreshold()` | Fits scale-free-topology statistics across candidate powers | Fit table, chosen power | Selects the network's weighting exponent |
| Network & Modules | Filtered matrix, chosen power | `WGCNA::blockwiseModules()` | Adjacency → TOM → clustering → dynamic cut → merge | Module assignments, eigengenes, dendrogram | Constructs the co-methylation network and detects modules |
| Module-Trait Analysis | Module eigengenes, sample-sheet trait column | `WGCNA::corPvalueStudent()`, `stats::p.adjust()` | Correlate each module eigengene against the trait; correct for multiple testing | Module-trait table, heatmap | Tests which co-methylation modules associate with the phenotype |
| Hub CpGs | Network matrix, module eigengenes, chosen module | `WGCNA::intramodularConnectivity.fromExpr()`, `WGCNA::bicor`/`cor` | Compute kME and kWithin; filter/rank | Ranked hub-CpG table with annotation | Identifies the most representative, central CpGs within a module |
| Results & Export | All prior stages | `stats::fisher.test()`, `DT`/`downloadHandler()` | Summarize parameters; test enrichment against an external panel; export | Summary table, enrichment table, 8 CSV downloads | Consolidates, contextualizes (against published results and an independent biomarker panel), and exports the full run |

---

## 23. Thesis Implementation Paragraph

The Methylomics WGCNA (Co-Methylation Network) submodule is implemented as a seven-tab pipeline (Data & Filtering, Sample QC, Soft Threshold, Network & Modules, Module-Trait Analysis, Hub CpGs, Results & Export) that treats individual CpG probes as network nodes and samples as observations. Starting from a beta- or M-value methylation matrix — either a preloaded whole-blood cohort or a user-uploaded dataset — the pipeline applies missingness filtering, optional `limma`-based covariate residualization, and top-N variability ranking before a `WGCNA::goodSamplesGenes()` quality gate. Co-methylation modules are constructed with `WGCNA::pickSoftThreshold()` for power selection and `WGCNA::blockwiseModules()` for network construction, topological-overlap computation, and dynamic module detection in a single memory-managed call, with results cached by a content-addressed key over both the input matrix and every network parameter. Module-trait relationships are evaluated by correlating each module's eigengene against a user-selected phenotype column, with Student's-t-derived p-values and Benjamini-Hochberg or Bonferroni correction restricted to real (non-grey) modules; hub CpGs within a module are ranked by module membership (kME) and intramodular connectivity (kWithin). Outputs include per-stage diagnostic plots, module and module-trait tables, a hub-CpG table with optional genomic annotation, an optional Fisher's-exact-test enrichment against an independent DMP/DMR biomarker panel, and eight downloadable CSV exports. This analysis is important for methylomics because it recovers coordinated, distributed methylation signal across many weakly-associated CpGs that neither single-CpG nor region-based testing is designed to detect, while collapsing the genome-wide multiple-testing burden to the small number of modules actually found.

## Tab-Based Thesis Description

Analysis begins on the Data & Filtering tab, where the shared methylation matrix and, when available, an optional sample sheet are stratified by sex, checked for scale and orientation, filtered for missingness, optionally residualized against known covariates via a chunked `limma` fit, and reduced to the top most-variable CpGs, with a final `WGCNA::goodSamplesGenes()` pass removing any residual near-zero-variance probe or sample; this produces the single filtered matrix every later tab depends on. Network parameters are then determined on the Soft Threshold tab, where `WGCNA::pickSoftThreshold()` evaluates scale-free-topology fit across a user-defined range of candidate powers and either the first power reaching a target fit, or the single best-observed power if none does, is carried forward. Co-methylation modules are generated on the Network & Modules tab by `WGCNA::blockwiseModules()`, which builds the weighted adjacency network, computes topological overlap, clusters CpGs, and merges closely related modules in one call, yielding a module-color assignment for every CpG and a module-eigengene matrix summarizing each module's coordinated methylation pattern per sample. Downstream tabs consume this same module structure in two complementary ways: Module-Trait Analysis correlates each module's eigengene against a chosen phenotype column, applying Student's-t p-values and multiple-testing correction restricted to real modules to determine which modules associate with the trait; Hub CpGs instead ranks the individual CpGs within one selected module by their module membership and network connectivity, surfacing the most representative members for follow-up. A final Results & Export tab consolidates every parameter used, offers an optional Fisher's-exact-test comparison of significant modules against an independently derived biomarker panel as convergent evidence, and provides CSV exports of every intermediate and final result. Taken together, this workflow is useful for identifying coordinated methylation patterns distributed across many CpGs that single-probe or single-region testing would not detect as a unit, while directly linking any such pattern back to a phenotype of interest through an explicit, correctable statistical test.
