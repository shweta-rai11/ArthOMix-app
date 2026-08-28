# Methylomics — ML Feature Selection

**Source file:** `ArthOMix/R/methylomics/mod_methyl_featureselection.R` (1,721 lines — the entirety of this submodule's UI, server, and helper logic lives in this one file)
**Registered as:** `mod_methyl_featureselection_config` / `_ui` / `_server`, id `"featureselection"`, title **"ML Feature Selection"**, group **"Biomarker modeling"** (`R/submodules_registry.R:47`, inside the `MX_MODULES` list)
**Parent module:** Methylomics (`MX_MODULES`, wired generically in `server.R:82-95`)

---

## Post-Audit Fixes Applied (2026-08-26)

Following the audit below, the underlying code was revised to correct every High-severity finding and several Moderate/Minor ones. This document otherwise still describes the original implementation as read; each corrected finding is marked **[FIXED]** inline with a short note on the change, so the audit trail (what was found → what was done) stays intact rather than being silently edited away. Summary of what changed in `mod_methyl_featureselection.R`:

1. **Crash on identical Reference/Comparison group [High → FIXED].** Added a `validate(need(...))` guard so selecting the same level for both now shows a clean message instead of an uncaught `factor()` "duplicated levels" error.
2. **"ANOVA (>2 groups)"/Kruskal-Wallis was structurally unreachable [High → FIXED].** Added `fs_uni_multigroup_data()`, which rebuilds the already-chosen probe panel against every sample with a non-missing value in the chosen Group/phenotype column (not just the two Reference/Comparison levels), re-applying the same sample-missingness/imputation logic Data & Filters already uses. ANOVA and Kruskal-Wallis now run a genuine multi-group test whenever the phenotype column actually has more than two levels among the eligible samples (verified with synthetic data: a true 4-group ANOVA now returns an F-like statistic, while the existing 2-group behavior is unchanged). `methyl_fs_univariate_select()` was also updated so its Δβ/direction effect-size columns are correctly left `NA` for a genuine multi-group run instead of silently comparing only the first/last level.
3. **Preloaded-cohort option gated on the wrong availability flag [High-relevance for the "preloaded data" ask → FIXED].** The "Preloaded whole-blood cohort" radio choice was gated on `METH_DATA_AVAILABLE` (only the lightweight precomputed-tables flag) instead of `METH_RAW_DATA_AVAILABLE` (the flag that actually determines whether the raw beta matrix — and hence `dataset$beta` — can ever be loaded). A deployment with the tables folder but not the ~2.1GB raw matrix would have offered a "Preloaded" option that could never actually work. Now gated on the correct flag, and the not-yet-loaded message was corrected to distinguish "not loaded yet" from "not available in this deployment."
4. **Random Forest class-imbalance weighting unreachable from the UI [Moderate → FIXED].** Added a "Class weighting" control (Equal/Balanced) to the Tree-Based Selection tab, wired to `randomForest`'s existing `classwt` parameter (which previously had no UI control setting it).
5. **Stability Selection denominator bug [Moderate → FIXED].** Degenerate (single-class) or failed-fit resamples previously still counted toward the denominator when computing each CpG's selection frequency, silently deflating every score. Now excluded from the denominator, with the result card explicitly reporting how many of the requested resamples were valid (confirmed with a stress test: an imbalanced cohort correctly excluded 5 of 60 resamples and adjusted the denominator to 55).
6. **RFE silent truncation [Minor → FIXED].** Requested feature-subset sizes exceeding the candidate panel size are now reported to the user instead of being silently dropped.
7. **k-NN imputation silent fallback [Minor → FIXED].** The result card now explicitly warns when k-NN was requested but the `impute` package isn't installed, rather than only being visible in the small "Imputation" value box.
8. **Incomplete reproducibility metadata [Minor → FIXED].** The `.rds` export's package-version record now includes `limma` and `pROC` (previously omitted despite being used in every run).
9. **"Leakage-safe" validation disclosure [Moderate → partially addressed].** The UI copy now explicitly discloses that Data & Filters' own missingness/variance filtering, top-variance cap, and imputation are computed once on the full dataset before the nested cross-validation runs, not redone per fold — and that the nested mode validates a re-derived Univariate+LASSO panel, not literally the panel in Selected Features. The underlying residual leakage itself was judged too architecturally invasive to eliminate without a live test environment to verify against (it would require re-plumbing Stage 1's filtering to run per-fold); this is left as documented follow-up work rather than attempted blind.

All changes were verified by parsing the modified file with R and running targeted functional smoke tests against synthetic data (2-group and genuine 4-group univariate runs, a stability-selection stress test with an imbalanced cohort, and a Random-Forest run with class weights) before being considered complete — not by static reading alone.

---

## 1. Module Overview

| | |
|---|---|
| Module | Methylomics |
| Submodule | ML Feature Selection |
| Purpose | Interactive, leakage-aware selection of a small panel of CpG probes from a much larger methylation feature space, using up to five independent statistical/ML selection methods combined into an overlap-based consensus |
| Scientific objective | Reduce hundreds of thousands of CpG sites to a biologically interpretable, statistically defensible panel associated with a two-group phenotype contrast (e.g. disease vs. control) |
| Computational objective | Provide a live, parameter-configurable, re-runnable Shiny workflow (not a static/precomputed report) that a user can drive end-to-end: filter → select (×5 methods) → combine → validate → export |
| Expected input | A CpG-by-sample beta or M-value matrix, either the app's own preloaded whole-blood cohort or a user-uploaded matrix, plus (usually) a phenotype/sample-sheet file defining a two-level group column |
| Expected output | A ranked/selected CpG list per method, a consensus CpG panel, per-CpG plots and tables, a cross-validated AUC estimate, and a single reproducible `.rds` bundle |
| Main approaches implemented | Univariate testing (6 sub-methods), LASSO/Elastic-Net regularization, Random Forest importance, three flavors of Recursive Feature Elimination (RF-RFE, Logistic-RFE, SVM-RFE), Meinshausen–Bühlmann-style stability selection, overlap/weighted consensus, optional pairwise-correlation redundancy reduction, and two validation strategies (frozen-panel and a partially leakage-safe nested design) |
| Number of sub-tabs | **9** (`mod_methyl_featureselection.R:594-604`) |
| Overall workflow | Data & Filters → five independent, explicitly Run-gated selection methods → Consensus/Overlap (+ optional correlation reduction) → Selected Features → Model & Export (validation + save/load) |

### Why methylomics needs feature selection (plain academic English)

A methylation array (450K/EPIC) or bisulfite-sequencing experiment measures methylation at hundreds of thousands of CpG sites per sample, while a typical cohort has, at best, a few hundred samples. This "wide, short" data shape — far more features (`p`) than samples (`n`) — causes several concrete problems that ML feature selection is designed to address:

- **Curse of dimensionality**: with `p >> n`, almost any classifier can achieve perfect separation on the training data by chance alone, so a model trained on all CpGs will overfit and generalize poorly.
- **Redundancy and correlation**: nearby CpGs (same CpG island, same gene promoter) are frequently co-methylated, so many features carry near-duplicate information; including all of them adds noise and instability without adding real signal.
- **Noise**: methylation beta values are bounded, heteroscedastic (variance depends on the mean), and technically noisy (batch, array, bisulfite-conversion effects), so many probes vary for purely technical reasons unrelated to phenotype.
- **Interpretability**: a clinically or biologically useful biomarker panel needs to be small — tens of CpGs, not hundreds of thousands — for a downstream classifier, qPCR/pyrosequencing validation assay, or mechanistic follow-up to be feasible.

Feature selection is therefore not merely a performance optimization; it is what turns an intractable, uninterpretable `p >> n` problem into a small, testable, biologically motivated CpG panel. This module implements that reduction as a multi-method, consensus-driven, user-auditable workflow rather than a single black-box selector — a design choice that is itself scientifically defensible: no single feature-selection method is universally "correct" for high-dimensional, correlated methylation data, so requiring several independent methods to agree is a hedge against any one method's idiosyncrasies (e.g. LASSO's tendency to pick one representative CpG out of a correlated cluster somewhat arbitrarily, or Random Forest's known bias toward high-cardinality/high-variance predictors).

---

## 2. Scientific Purpose

The submodule exists to answer: *"Which CpG sites best distinguish two phenotype groups (e.g. rheumatoid arthritis vs. control) in a way that is robust to the choice of statistical method, robust to resampling of the same cohort, and — optionally — robust to the optimism of evaluating a panel on the same data used to select it?"* Every one of the five methods, the consensus step, and the validation step is a different angle on this same question: univariate testing asks "does this CpG differ marginally by group?"; regularization and tree-based/RFE ask "does this CpG carry independent, non-redundant, multivariate predictive information?"; stability selection asks "is this CpG chosen consistently across resamples, or only by chance in one fit?"; and the consensus/correlation-reduction/validation stages ask "how many of these independent answers agree, and how well does the resulting panel generalize?"

## 3. Computational Purpose

Each of the five selection methods needs the same starting matrix (filtered/imputed beta and M-value matrices restricted to a two-group contrast) but a different representation and fitting routine (a design matrix for `limma`, a transposed `sample × probe` matrix for `glmnet`/`randomForest`/`caret`/`e1071`). The module's computational contribution is to (a) build that shared, filtered matrix exactly once (Stage 1: Data & Filters), (b) let every method reuse it independently and be individually re-run without invalidating the others, (c) reconcile the resulting CpG-ID sets into one consensus table, and (d) package the whole configuration — not just the final CpG list — into one reproducible artifact.

## 4. Number and Names of Sub-Tabs

Exactly **9** `tabPanel`s inside one `tabsetPanel(id = ns("fs_subtabs"))` (`mod_methyl_featureselection.R:594-604`):

1. Data & Filters
2. Univariate Selection
3. LASSO
4. Tree-Based Selection
5. RFE / Wrapper Selection
6. Stability Selection
7. Consensus / Overlap
8. Selected Features
9. Model & Export

## 5. Overall Architecture

- **Single-file module.** UI (`mod_methyl_featureselection_ui`, lines 590-607) and server (`mod_methyl_featureselection_server`, lines 621-1721) live in the same file; every tab's UI is rendered lazily via `uiOutput()` + `renderUI()` rather than being built eagerly, so tabs the user hasn't opened don't force their inputs to exist.
- **Explicit Run-button gating.** No selection method fires automatically on input change. Every stage has its own `actionButton` (`filters_run_btn`, `uni_run_btn`, `reg_run_btn`, `rf_run_btn`, `rfe_run_btn`, `stab_run_btn`, `consensus_run_btn`, `corr_reduce_btn`, `validate_run_btn`) driving an `eventReactive()`, and a paired `*_has_run` `reactiveVal(FALSE)` that gates whether that stage's result UI renders at all. This is a deliberate, comment-documented design choice (line 6-8: "each behind its own explicit 'Run' button") so a slow method (Random Forest, RFE, Stability Selection) never recomputes just because the user tweaked an unrelated slider.
- **Cascading invalidation.** Every downstream stage's `*_has_run` flag is reset to `FALSE` whenever the upstream `fs_filter_result()` changes (i.e. whenever "Run Filters" is re-clicked) — see `observeEvent(fs_filter_result(), X_has_run(FALSE))` at lines 939, 1017, 1094, 1174, 1259, 1356, 1604. Re-running an individual method's own parameters without re-clicking its own Run button does **not** invalidate it — the previous result stays displayed until the user re-runs that specific method (this is intentional, not a staleness bug: see §15).
- **Two mutually exclusive data sources**, resolved once per session into one shape by `fs_active_source()` (lines 659-673): the app's preloaded whole-blood cohort (`dataset$beta`, a shared `reactiveValues` populated by the Methylomics Dataset tab), or a user upload parsed by this module's own parsers.
- **Shared results channel.** On every successful Consensus run, the module writes a small summary (`selected_cpgs`, `n_selected`, per-method counts, the consensus rule text) into the shared `results$featureselection` slot (lines 1372-1381), the same `methyl_results` `reactiveValues` object every other Methylomics submodule reads from (`server.R:93`). As of this reading, no other Methylomics submodule actually reads `results$featureselection` back (a repo-wide `grep` found no consumer) — the module year currently writes to a shared channel nothing downstream yet reads.

## 6. Input Data

Two mutually exclusive routes, chosen via `radioButtons(ns("fs_source"), ...)` (lines 691-693):

**A. Preloaded whole-blood cohort** (`fs_active_source()`, lines 665-672): `dataset$beta` — populated by the Methylomics Dataset tab (`mod_methyl_dataset.R`) from `load_default_meth_matrix()`/`load_default_meth_pheno()` (`global.R:340-385`+). This is the GSE42861 (Liu et al. 2013) whole-blood rheumatoid-arthritis cohort: 689 samples, 412,492 QC'd CpGs on the Illumina 450K array, beta scale (`global.R:275-343`; the exact same QC'd cascade documented for the Methylomics QC submodule). Only available when `METH_DATA_AVAILABLE` (and the raw matrix specifically, `METH_RAW_DATA_AVAILABLE`) is `TRUE` for the deployment; the UI's `radioButtons` choice list itself drops the "Preloaded" option entirely when unavailable (line 692).

**B. Upload my own** (lines 694-711): a beta/M-value matrix (CSV/TSV/TXT/RDS, probe IDs in the first column or rownames, samples as columns — `methyl_fs_parse_matrix_any()`, lines 52-55, dispatching to `methyl_parse_matrix()` in `parse_upload.R:11-37` for text formats or the module's own `methyl_fs_parse_matrix_rds()` for `.rds`), plus three optional files: a phenotype/sample-metadata sheet (`methyl_parse_sample_sheet()`, `parse_upload.R:40-46`), a cross-reactive probe exclusion list (`methyl_parse_probe_list()`, `parse_upload.R:50-57`), and a probe-MAF table (`methyl_parse_maf_list()`, `qc.R:131-144`).

For an upload, the array type (for annotation-dependent filters) and value scale (beta vs. M-value) are user-declared, with an "Auto-detect" option (`methyl_fs_detect_scale()`, lines 60-68) that samples up to 20,000 finite values and checks whether the 0.1–99.9% quantile range sits inside `[-0.05, 1.05]`. Both routes converge on the same shape: `list(mat, sheet, scale, frac_out_of_range, array_type, source_label, preloaded)`.

**Dimensions.** Preloaded: 412,492 CpGs × 689 samples (before this module's own filters run). Upload: fully user-dependent — the module makes no assumption about the number of probes or samples beyond the hard minimums enforced in Stage 1 (§7).

---

## 7. Tab 1 — Data & Filters

### Purpose
Establishes the one shared, filtered, imputed feature matrix (beta + derived M-value) and the two-group phenotype contrast that every other tab in this submodule reuses without recomputation.

### User Inputs
| Input ID | Control | Label | Default | Notes |
|---|---|---|---|---|
| `fs_source` | radioButtons | Data source | `preloaded` (or `upload` if preloaded unavailable) | Required |
| `fs_matrix_file` | fileInput | Beta/M-value matrix | — | Required only in upload mode |
| `fs_upload_array_type` | selectInput | Array type | `450K` | From `METHYL_ARRAY_TYPES` |
| `fs_scale_mode` | radioButtons | Value scale | `auto` | auto / force-beta / force-M |
| `fs_sheet_file` | fileInput | Phenotype file | — | Optional but required for any group-based method |
| `fs_exclusion_file` | fileInput | Cross-reactive exclusion list | — | Optional |
| `fs_maf_file` | fileInput | MAF table | — | Optional |
| `fs_group_col` | selectInput | Group/phenotype column | best-guess via `methyl_fs_guess_group_col()` | Required once a sheet is loaded |
| `fs_ref_group` / `fs_comp_group` | selectInput | Reference / comparison level | first two distinct levels | Required; **nothing prevents selecting the same level for both — see §14 edge case** |
| `fs_max_na_probe` | numericInput | Max missingness per CpG (%) | 5 | 0–100 |
| `fs_max_na_sample` | numericInput | Max missingness per sample (%) | 20 | 0–100 |
| `fs_impute` | checkboxInput | Impute remaining missing values | TRUE | |
| `fs_impute_method` | radioButtons | Imputation method | `median` | median / k-NN (if `impute` package installed) |
| `fs_var_metric` | radioButtons | Variance metric | `variance` | variance / SD / IQR |
| `fs_var_min` | numericInput | Minimum threshold | 0 | |
| `fs_max_probes` | selectInput | Top-variance cap | `5000` | 500/1,000/5,000/10,000/Custom |
| `fs_filter_snp` | checkboxInput | Remove SNP-associated probes | FALSE | Requires array annotation |
| `fs_filter_sexchr` | checkboxInput | Remove sex-chromosome probes | FALSE | Requires array annotation |
| `fs_filter_crossreactive` | checkboxInput | Apply uploaded exclusion list | FALSE | No-op without an uploaded list |
| `fs_filter_maf` | checkboxInput | Apply uploaded MAF filter | FALSE | No-op without an uploaded table |
| `fs_max_maf` | numericInput | Max MAF | 0.05 | 0–0.5 |
| `filters_run_btn` | actionButton | Run Filters | — | Fires `fs_filter_result()` |

### Input Data
Comes from `fs_active_source()` (§6). Dimensions: `nrow(mat) × ncol(mat)` probes × samples, unfiltered.

### Processing (execution order, `fs_filter_result`, lines 802-872)
1. Resolve sample IDs between the matrix and the phenotype sheet (`methyl_sheet_sample_ids()`, `qc.R:456-461`); require ≥10 overlapping samples.
2. Subset to samples whose group-column value is exactly the chosen reference or comparison level; require ≥10 such samples and ≥3 per group.
3. Convert to the beta scale if the detected/declared scale is M-value (`2^m / (1 + 2^m)`).
4. **Sample-axis** missingness filter (`methyl_fs_sample_missing_ok()`, lines 94-98); require ≥6 samples remaining.
5. **Probe-axis** filters, each independently computed and then ANDed together: missingness (`methyl_filter_missing()`), variance/SD/IQR (`methyl_filter_variance()`/`methyl_filter_sd()`/`methyl_fs_filter_iqr()`), and — if checked — SNP-associated (`methyl_filter_snp()`), sex-chromosome (`methyl_filter_sex_chr()`), cross-reactive (`methyl_filter_cross_reactive()`), MAF (`methyl_filter_maf()`). Require ≥10 probes surviving.
6. Build a retention cascade table (`methyl_probe_retention_cascade()`).
7. **Top-variance cap**: subset to the top-N most-variable surviving probes (`methyl_fs_cap_top_variance()`), appended as an extra cascade row.
8. **Imputation** (only if any `NA` remain): per-probe median, or k-NN if the optional `impute` package is installed and selected (`methyl_fs_impute()`, lines 104-114 — silently falls back to median if `impute` isn't installed, with no user-facing warning that the requested k-NN was not actually used).
9. Derive the M-value matrix from the now-complete beta matrix (`methyl_beta_to_mvalue()`).

### Functions Used
See the [Function-by-Function Code Audit](#function-by-function-code-audit) below for every function in this list, individually documented: `methyl_fs_parse_matrix_any`, `methyl_fs_detect_scale`, `methyl_fs_guess_group_col`, `methyl_sheet_sample_ids`, `methyl_fs_sample_missing_ok`, `methyl_filter_missing`, `methyl_filter_variance`/`methyl_filter_sd`/`methyl_fs_filter_iqr`, `methyl_filter_snp`, `methyl_filter_sex_chr`, `methyl_filter_cross_reactive`, `methyl_filter_maf`, `methyl_probe_retention_cascade`, `methyl_fs_cap_top_variance`, `methyl_fs_impute`, `methyl_beta_to_mvalue`, `methyl_get_annotation`.

### Output
Three value boxes (CpGs retained, samples/contrast, imputation method), a bulleted list of per-filter notes, and a retention-cascade bar chart (`methyl_plot_cascade`, `output$fs_cascade_plot`). The underlying `fs_filter_result()` object (`beta`, `m`, `grp`, cascade, filter notes) is what every subsequent tab consumes — nothing in this tab is itself downloadable.

### Scientific Interpretation
This tab operationalizes standard methylation QC practice — remove probes/samples with excessive missingness, remove near-constant (low-variance) probes that cannot carry group-discriminating signal, optionally remove technically unreliable probe classes (SNP-overlapping, cross-reactive, sex-linked, common-variant) — before any supervised selection method sees the data. All of these filters are **unsupervised** (they never look at the group label to decide which probes/samples to drop), which is the methodologically correct order relative to supervised selection.

### Code Audit
- **Correct / appropriate:** filter order (missingness → variance → optional QC filters → top-variance cap → impute → M-value), all filters returning an explicit `note` string, minimum-sample-size guards at every subsetting step, imputation deliberately happening after — not before — the statistics that determine which probes survive.
- **Minor: k-NN imputation silently degrades to median with no user-facing notice when the `impute` package isn't installed (line 105-108) — the returned `method_used` string does correctly report which method actually ran, but only if the user opens the result card to check it. [FIXED — see "Post-Audit Fixes Applied" above.]**
- **Moderate — no cross-reactive/MAF reference list is bundled.** `methyl_filter_cross_reactive()`/`methyl_filter_maf()` are pure pass-through no-ops unless the user supplies their own exclusion/MAF list (`qc.R:106-119`, `146-160`); this is explicitly, deliberately disclosed in the source comments ("fabricating one would violate this project's own evidence-based-methods requirement") rather than silently omitted, so it is a documented scope limitation, not a bug.
- **High — unvalidated identical reference/comparison group selection can crash the reactive. [FIXED — see "Post-Audit Fixes Applied" above.]** Nothing checks `input$fs_ref_group != input$fs_comp_group` before line 815's `factor(grp_raw[keep_s], levels = c(input$fs_ref_group, input$fs_comp_group))`. If a user selects the same level for both (nothing in the UI prevents this — `fs_ref_group`/`fs_comp_group` are two independent `selectInput`s over the same level list), R's `factor()` throws `"duplicated levels in factor are not allowed"` — an uncaught error rather than a clean `validate()`/`need()` message, surfacing as Shiny's generic red error box instead of this module's usual inline guidance. See §14.
- **Minor — top-variance cap is architecturally separate from the other filters** (documented in-code at lines 847-850) — correct as implemented, just worth noting for a reader tracing the cascade table, since it is appended as an extra row rather than folded into the same `filters` list the `keep_probe` mask was built from.

---

## 8. Tab 2 — Univariate Selection

### Purpose
Tests each candidate CpG's marginal (single-probe-at-a-time) association with the phenotype contrast, independent of every other probe — the fastest and most classically interpretable of the five methods, and the only one offering explicit covariate adjustment.

### User Inputs
| Input ID | Control | Label | Default |
|---|---|---|---|
| `uni_method` | selectInput | Statistical method | `moderated_t` |
| `uni_covariates` | checkboxGroupInput | Covariate adjustment (linear-model methods only) | none selected |
| `uni_rule` | radioButtons | Selection rule (FDR / p-value / Top N) | `fdr` |
| `uni_threshold` | numericInput | Threshold | 0.05 |
| `uni_top_n` | numericInput | Top N | 100 |
| `uni_dbeta_min` | numericInput | Minimum \|Δβ\| (categorical methods only) | 0 |
| `uni_run_btn` | actionButton | Run Univariate Selection | — |

Eight `uni_method` choices are offered: Moderated t-test (limma), t-test, ANOVA (>2 groups), Wilcoxon rank-sum, Kruskal-Wallis, Linear regression, Pearson correlation, Spearman correlation. A conditional note (lines 911-913) already discloses that the three "continuous phenotype" options actually run against the same reference/comparison contrast, numerically coded — see the audit finding below for the more serious consequence of that same structural fact.

### Input Data
`fs_filter_result()$m` (linear-model methods) or `$beta` (rank-based methods), and `$grp` — always exactly a **two-level** factor, because Data & Filters enforces a binary reference/comparison contrast (§7). Optional covariates come from the active source's phenotype sheet, matched by sample ID (lines 946-950).

### Processing
`is_linear <- uni_method %in% c("moderated_t","t_test","anova","linear_regression","pearson")` selects the engine (line 944): linear-model methods go through `methyl_fs_univariate_linear()`; rank-based methods (`wilcoxon`, `kruskal`, `spearman`) go through `methyl_fs_univariate_rank()`. The ranked table is then thresholded by `methyl_fs_univariate_select()`.

### Functions Used

#### `methyl_fs_univariate_linear()` (lines 124-151)
**Source/package:** custom, on top of `limma::lmFit()`/`limma::eBayes()`/`limma::topTable()`.
**Purpose:** one shared engine for t-test / moderated t-test / ANOVA / linear regression / Pearson correlation — all are just different design-matrix/coefficient choices on the same `limma` fit.
**Arguments:** `m_mat` (M-values), `y` (group factor or numeric), `covariates` (optional data.frame), `mode`.
**Processing:** builds a design matrix (`~grp + covariates` for categorical modes, `~y + covariates` for continuous modes), checks the design is full rank (`qr(design)$rank == ncol(design)`), fits `limma::lmFit()` then `limma::eBayes(trend = (mode == "moderated_t"))`, and reads off `topTable()`.
**Output:** a data.frame of `cpg, statistic, p, fdr, logfc`.
**Downstream dependency:** `methyl_fs_univariate_select()`, then Consensus.
**Audit:** correctly guards against rank-deficient designs with a clean `validate()` message (line 140) rather than letting `lmFit()` throw; correctly chooses `trend=TRUE` only for the moderated-t option (an array-appropriate empirical-Bayes trend, matching `limma`'s own guidance for methylation-scale data — see §12).

#### `methyl_fs_univariate_rank()` (lines 156-185)
**Source/package:** custom, on top of base R `stats::wilcox.test()`/`stats::kruskal.test()`/`stats::cor.test()`.
**Purpose:** shared engine for the three non-parametric/rank-based options.
**Processing:** iterates `apply(mat, 1, ...)` per probe (not vectorized — acceptable at the post-filter candidate-panel scale, per the file's own comment at line 154-155, since this only ever runs on the already-capped ≤10,000-probe matrix, not a genome-wide scan).
**Output:** `cpg, statistic, p, fdr` (BH-adjusted).
**Audit:** correct and defensively wrapped in `tryCatch()` per probe so one degenerate probe (e.g. all-identical values, breaking `wilcox.test`) doesn't abort the whole scan — that probe simply gets `NA` statistics instead.

#### `methyl_fs_univariate_select()` (lines 190-215)
**Purpose:** applies the chosen selection rule (FDR/p-value/Top-N) and attaches Δβ/direction effect-size columns computed directly from the beta-scale matrix (not from the M-value fit statistics), so effect sizes stay in the interpretable [0,1] units even when the test itself ran on M-values.
**Audit:** correctly recognizes the "is this actually a continuous phenotype" case via `is_continuous <- !anyNA(yf) && length(unique(as.character(y))) > 6` — but because `y` here is always the same 2-level `grp` factor (never a genuine continuous phenotype column; see below), `is_continuous` is always `FALSE` in practice, and `df$dbeta`/`df$direction` are always computed from the categorical branch (mean-difference between the two groups) regardless of which `uni_method` was selected — which is actually the right outcome, since there is no genuine continuous phenotype anywhere in this workflow to begin with.

### Output
A ranked/annotated table (`DT::dataTableOutput(ns("uni_table"))`), a "N of M CpGs selected" summary, and a CSV download (`uni_download`).

### Scientific Interpretation
A small, FDR- or Δβ-significant p-value means this CpG's methylation level differs between the two groups more than expected by chance, marginally (ignoring every other CpG). It does **not** mean the CpG is causal, mechanistically involved, or even non-redundant with other selected CpGs — two perfectly co-methylated CpGs will both be "significant" here even though they carry identical information.

### Code Audit
- **Correct / appropriate:** shared `limma` engine for all parametric options is efficient and statistically standard for methylation-array differential analysis; BH FDR correction applied consistently; per-probe `tryCatch` in the rank engine prevents one bad probe from crashing the whole run.
- **High — the "ANOVA (>2 groups)" option is structurally unreachable and silently degrades to a two-group test. [FIXED — see "Post-Audit Fixes Applied" above; ANOVA/Kruskal-Wallis now genuinely test all levels of the phenotype column when more than two are present.]** Data & Filters (§7) always constructs `grp` as a factor with **exactly two** levels (`levels = c(ref_group, comp_group)`, line 815), regardless of how many distinct values the original phenotype column had. Inside `methyl_fs_univariate_linear()`'s `anova` branch, `coef_name <- grep("^grp", colnames(design), value = TRUE)` (line 138) therefore always has length 1, so the caller's branch condition `identical(mode, "anova") && length(coef_name) > 1` (line 143) is never `TRUE`, and the code always falls through to the same single-coefficient `topTable(..., sort.by = "P")` path used by `t_test`/`moderated_t`. **In the current implementation, choosing "ANOVA (>2 groups)" in the dropdown produces results structurally identical in kind to a t-test**, not a genuine multi-group F-test — the label promises functionality the surrounding data model (which enforces a binary contrast one tab upstream) never lets this code path exercise. This should be verified against product intent: either the label should be corrected/removed until Data & Filters supports selecting more than two groups, or Data & Filters should be extended to allow a genuine multi-level contrast.
- **Minor — "Kruskal-Wallis (>2 groups)" has the same structural limitation** but is not misleading in the same way: `kruskal.test()` run on exactly two groups is a valid (if underpowered relative to Wilcoxon) statistical test, so no code path goes silently unreachable, it is simply redundant with the Wilcoxon option under the current binary-only data model.
- **Minor — the "continuous phenotype" methods do not use a genuine continuous phenotype.** Linear regression / Pearson / Spearman all operate on `as.numeric(grp)` — the same binary contrast numerically coded as factor levels (1/2, not literally 0/1 as the UI note states) — which the UI already partially discloses (lines 911-913) but understates: these are not "linear regression against a continuous phenotype", they are a 2-group test wearing continuous-method clothing, and are mathematically equivalent in direction/significance to the corresponding categorical test (an affine 1/2-vs-0/1 recoding does not change a Pearson correlation's sign or p-value). Not a computational bug, but a real usability/interpretation risk: a reader of the Selected Features / exported CSV could reasonably believe a genuine continuous predictor was modeled.

---

## 9. Tab 3 — LASSO (Regularization)

### Purpose
Fits an L1/elastic-net-penalized logistic regression across all candidate CpGs simultaneously, so the surviving non-zero coefficients represent probes carrying independent (not merely marginal) discriminative signal after accounting for correlation with other probes.

### User Inputs
| Input ID | Control | Label | Default |
|---|---|---|---|
| `reg_alpha` | sliderInput | Alpha (1=LASSO, 0=Ridge) | 1 |
| `reg_lambda_choice` | radioButtons | Lambda (`lambda.min`/`lambda.1se`) | `lambda.min` |
| `reg_nfolds` | numericInput | CV folds | 10 |
| `reg_nlambda` | numericInput | Number of lambda values | 100 |
| `reg_seed` | numericInput | Random seed | 1234 |
| `reg_standardize` | checkboxInput | Standardize predictors | TRUE |
| `reg_class_weight` | radioButtons | Class weighting (equal/balanced) | equal |
| `reg_max_selected` | numericInput | Max selected CpGs | 200 |
| `reg_coef_threshold` | numericInput | Coefficient magnitude threshold | 0 |
| `reg_run_btn` | actionButton | Run LASSO / Elastic Net | — |

### Input Data
`X <- t(fs_filter_result()$m)` — a `sample × CpG` matrix of M-values (transposed from the probe × sample storage orientation used everywhere else), and `r$grp` as the binomial outcome.

### Processing / Functions Used

#### `methyl_fs_lasso_fit()` (lines 221-241)
**Package:** `glmnet` (`glmnet::cv.glmnet`, `stats::coef`).
**Arguments:** `alpha` (elastic-net mixing), `nfolds` (clamped: `max(3, min(nfolds, min(table(y))))` — defensively bounded by the smaller class's size so `cv.glmnet` never receives more folds than the rarer class can support), `lambda_choice`, `weights` (from "balanced" class-weighting, computed as `max(table)/table` per class), `seed`, `max_selected`, `coef_threshold`.
**Processing:** `set.seed(seed)`; `glmnet::cv.glmnet(X, y, family="binomial", ...)`; reads coefficients at the chosen lambda, drops zero/below-threshold coefficients, optionally caps to the top-`max_selected` by `|coefficient|`.
**Output:** `list(cv, lambda_used, coefficients, selected_ids, ranked)`.
**Downstream dependency:** Consensus, Selected Features, `fs_model_export`.
**Audit:** correct, appropriately defensive `nfolds` clamp, `tryCatch` around the `cv.glmnet` call surfaced as a clean `validate()` message rather than a raw error.

### Output
A `glmnet` CV-error plot (`plot(fs_reg_result()$cv)` — the package's own base-R plot method, not a `ggplot`, so it does not follow this app's `theme_arthomix()` styling, unlike every other plot in the module), a ranked coefficient table, and a CSV download.

### Scientific Interpretation
A non-zero LASSO coefficient at the CV-selected lambda means this CpG contributes independent predictive signal to the group classification *given the other CpGs also in the model* — a materially different claim from univariate significance. Because L1 penalization tends to pick one representative from a cluster of correlated CpGs somewhat arbitrarily, a CpG's *absence* from the LASSO panel does not mean it is uninformative, only that a correlated neighbor was chosen instead.

### Code Audit
- **Correct / appropriate:** `family="binomial"` is the right choice given the data model's binary contrast; fold-count and class-weight handling are both sound; seed is user-visible and passed through to `fs_model_export` for reproducibility.
- **Minor:** the raw `glmnet::plot.cv.glmnet()` base-R plot breaks visual consistency with the rest of the module's `ggplot2`/`theme_arthomix()` styling — cosmetic only.

---

## 10. Tab 4 — Tree-Based Selection (Random Forest)

### Purpose
Ranks CpGs by their contribution to a Random Forest classifier's split quality — a non-linear, interaction-aware alternative to the two linear methods above.

### User Inputs
`rf_ntree` (1000), `rf_mtry_mode`/`rf_mtry_manual` (auto=`sqrt(p)`), `rf_nodesize` (1), `rf_maxnodes_unlimited`/`rf_maxnodes`, `rf_replace` (TRUE), `rf_importance_type` (Gini/Accuracy), `rf_rule` (Top N / Threshold / Percentile), `rf_top_n`, `rf_seed` (1234), `rf_run_btn`.

### Input Data
`X <- t(fs_filter_result()$m)`, `r$grp`.

### Functions Used

#### `methyl_fs_rf_fit()` (lines 247-282)
**Package:** `randomForest` (`randomForest::randomForest`).
**Processing:** resolves `mtry` (auto `sqrt(p)` or manual, clamped to `[1, p]`); builds an argument list conditionally — `maxnodes`/`classwt`/`sampsize` are only added to the call when actually supplied, because `randomForest()`'s own `missing(arg)`-style internal defaults break if an explicit `NULL` is passed instead of the argument being entirely absent (documented in-code at lines 256-260, a correct and non-obvious defensive pattern); `set.seed(seed)`; fits via `do.call()`; ranks by the chosen importance column (`MeanDecreaseGini`/`MeanDecreaseAccuracy`); applies the selection rule.
**Output:** `list(rf, mtry, importance_type, importance, selected_ids, ranked)`.
**Audit:** the `do.call()`-with-conditional-args pattern is correct and necessary; the function signature accepts `classwt` and `sample_fraction` parameters for class-imbalance handling, but **neither is exposed anywhere in `rf_ui`** (lines 1065-1091 offer no corresponding `numericInput`) — see §15, Minor.

### Output
A horizontal bar chart of the top 25 CpGs by importance (`ggplot`, `theme_arthomix()`), a full ranked table with a `selected` flag column, and a CSV download.

### Scientific Interpretation
Higher importance means this CpG contributes more to reducing classification impurity (Gini) or accuracy loss (permutation) across the forest's trees — a measure that, unlike the LASSO coefficient, can capture non-linear and interaction effects, but is also known to be biased toward high-variance/high-cardinality predictors, which is one motivation for cross-checking against the other four methods via Consensus.

### Code Audit
- **Correct / appropriate:** conditional-argument `do.call()` pattern, sensible clamping of `mtry`, both importance metrics correctly require `importance=TRUE` (always set).
- **Minor — dead/unreachable parameters. [PARTIALLY FIXED — see "Post-Audit Fixes Applied" above.]** `methyl_fs_rf_fit()`'s `classwt` and `sample_fraction` arguments exist specifically to correct for class imbalance, but no UI control ever sets them (always `NULL`), so this capability is present in the helper function but not reachable by a user of this tab — a real gap if the underlying cohort has an imbalanced group split. A "Class weighting" control was added for `classwt`; `sample_fraction` remains unexposed (a rarer, more advanced knob than the class-weighting fix this audit prioritized).

---

## 11. Tab 5 — RFE / Wrapper Selection

### Purpose
Wrapper-style selection: repeatedly refits a model on shrinking feature subsets and picks the subset size that empirically minimizes cross-validated error, rather than relying on a single importance score's implicit threshold.

### User Inputs
`rfe_flavor` (Random Forest RFE / Logistic Regression RFE / SVM-RFE), `rfe_sizes` (comma-separated, pre-filled by `methyl_fs_rfe_sizes()`), `rfe_folds` (5), and — only for the SVM flavor — `rfe_svm_cost`, `rfe_svm_tolerance`, `rfe_svm_panel_mode` (auto-minimize-CV-error / manual), `rfe_svm_manual_k`; `rfe_seed` (1234); `rfe_run_btn`.

### Input Data
`X <- t(fs_filter_result()$m)`, `r$grp`.

### Functions Used

#### `methyl_fs_rfe_sizes()` (lines 288-293)
Parses/validates the user's comma-separated size list against `p`; falls back to a fixed default ladder (`10,25,50,100,250,500,1000,2000,5000`, capped at `p`) if parsing yields nothing usable. **Sizes exceeding `p` are silently dropped with no user-facing note** (unlike every filter in Stage 1, which always emits an explicit `note` string) — a minor inconsistency (§15).

#### `methyl_fs_rfe_run()` (lines 295-306)
**Package:** `caret` (`caret::rfe`, `caret::rfeControl`, `caret::rfFuncs`/`caret::lrFuncs`, `caret::predictors`).
**Processing:** `caret::rfeControl(functions = <rfFuncs|lrFuncs>, method="cv", number=max(3,k))`; `set.seed(seed)`; `caret::rfe(x, y, sizes, rfeControl=ctrl)`.
**Output:** `list(fit, sizes, optimal_size, selected_ids, perf, ranked)`.

#### `methyl_fs_svm_rfe_rank()` / `methyl_fs_svm_rfe_curve()` / `methyl_fs_svm_rfe_run()` (lines 313-346)
**Package:** `e1071` (`e1071::svm`, linear kernel).
**Provenance:** explicitly noted as "ported near-verbatim from the transcriptomics module's own SVM-RFE" (lines 308-312) — a deliberate, disclosed code-reuse decision, not independent duplication.
**Processing:** classic backward elimination — repeatedly fits a linear SVM on the remaining features, computes each feature's squared weight (`(t(coefs) %*% SV)^2`), drops the smallest, and records the elimination order as a full ranking (most-to-least important); then sweeps `k = 1..p` along that ranking, refitting a `k`-feature linear SVM with `cross = folds`-fold CV at each step, and picks the `k` minimizing CV classification error (or a user-specified manual `k`).
**Audit:** computationally the most expensive of the three RFE flavors (`p` full-SVM refits for ranking, then another `p` CV-fit sweep) — appropriate only because it runs on the already-capped candidate panel (≤10,000 probes after Stage 1), not the full CpG space; correctly reuses a proven pattern from elsewhere in the codebase rather than reimplementing SVM-RFE independently.

### Output
Either a CV-error-vs-`k` curve (SVM flavor) or a performance-vs-subset-size curve (RF/Logistic flavor), both with the chosen optimal size marked by a dashed vertical line; a ranked table; a CSV download.

### Scientific Interpretation
The RFE-selected subset size is the one that empirically minimized cross-validated error *for the specific wrapped model* (RF, logistic, or linear SVM) — a data-driven alternative to an arbitrary fixed threshold, but one that inherits whatever biases the wrapped model itself has (e.g. a linear-SVM RFE cannot detect non-linear discriminative CpGs that a Random-Forest-wrapped RFE might).

### Code Audit
- **Correct / appropriate:** all three flavors correctly seed before their respective CV/fitting; `caret::rfe`'s own internal CV avoids hand-rolled fold logic for the RF/Logistic flavors.
- **Minor: silent truncation of out-of-range `rfe_sizes` entries with no user-facing note (unlike Stage 1's filters). [FIXED — see "Post-Audit Fixes Applied" above.]**

---

## 12. Tab 6 — Stability Selection

### Purpose
Addresses a well-known weakness of any single LASSO fit on correlated, high-dimensional data — the selected set can be unstable under small perturbations of the data — by repeatedly resampling and tabulating how often each CpG is selected by a fixed-lambda LASSO base selector (Meinshausen–Bühlmann-style stability selection).

### User Inputs
`stab_type` (Bootstrap / Repeated k-fold / Subsampling), `stab_n_resamples` (50), `stab_fraction` (0.8, bootstrap/subsampling), `stab_k` (5) / `stab_repeats` (5, repeated k-fold only), `stab_freq_threshold` (0.7), `stab_seed` (1234), `stab_run_btn`.

### Input Data
`X <- t(fs_filter_result()$m)`, `r$grp`.

### Functions Used

#### `methyl_fs_stability_resample_indices()` (lines 356-369)
Generates the resample index sets for whichever of the three schemes is chosen — `sample.int(n, n, replace=TRUE)` for bootstrap, a fixed-fraction `sample.int(n, sz, replace=FALSE)` for subsampling, or `caret::createMultiFolds()` for repeated k-fold.

#### `methyl_fs_stability_run()` (lines 371-397)
**Package:** `glmnet`.
**Processing:** fits one *reference* `cv.glmnet()` on the full data to fix a single reference lambda (`ref_lambda`); for each resample, refits `glmnet::glmnet()` at that fixed lambda on the resampled subset and records which CpGs have non-zero coefficients; `freq <- rowMeans(sel_mat)` is each CpG's selection frequency across all `length(idx_sets)` resamples; the panel is every CpG with `freq >= freq_threshold`.
**Output:** `list(n_resamples, ref_lambda, ranked, selected_ids)`.
**Downstream dependency:** Consensus, Selected Features, `fs_model_export`.
**Audit:** this is a documented, deliberate **simplification** of full Meinshausen–Bühlmann stability selection — the original formulation sweeps a lambda path and aggregates selection probability across it; this implementation fixes one reference lambda from a single initial CV fit and resamples only at that lambda. That is a reasonable, disclosed engineering trade-off for interactive runtime, but is a genuine methodological simplification worth stating explicitly in any thesis description of "stability selection" here, since the term has a more specific meaning in the original literature than what is implemented.

### Output
A bar chart of the top 30 CpGs by selection frequency with the threshold marked, a ranked table, a CSV download.

### Scientific Interpretation
A high selection frequency means this CpG's LASSO selection is robust to resampling noise, not an artifact of the particular sample composition in one fit — a stability, not an effect-size or significance, claim.

### Code Audit
- **Correct / appropriate:** the general resampling-loop abstraction correctly unifies three resampling schemes behind one interface; `set.seed()` applied once before generating index sets, ensuring reproducibility of the resample sets themselves.
- **Moderate — degenerate resamples silently count as "zero features selected" rather than being excluded from the denominator. [FIXED — see "Post-Audit Fixes Applied" above; verified with a stress test excluding 5 of 60 resamples correctly.]** When a resample's `y[idx]` contains only one class (line 382, `next`) or the per-resample `glmnet::glmnet()` fit fails (line 385, `next`), that resample's column in `sel_mat` is left as all-zeros (the pre-allocated default), but `n_resamples <- length(idx_sets)` (line 395) — the *requested* resample count — is what both the returned metadata and the result UI's "≥X% of N resamples" text (line 1277) report to the user, and `freq <- rowMeans(sel_mat)` (line 389) divides by that same total, not by the count of resamples that actually produced a valid fit. Every CpG's `selection_frequency`/`stability_score` is therefore silently deflated by however many resamples degenerated to a single class or failed to fit, and the user has no visibility into how many of the "N resamples" reported were actually valid. This is most likely to matter with a small, imbalanced cohort and a high `stab_fraction`/small `stab_k`, where single-class subsamples become more probable.

---

## 13. Tab 7 — Consensus / Overlap

### Purpose
Combines the CpG-ID lists from however many of the five methods the user has actually run into one ranked, weighted overlap table, and offers an optional secondary pass to prune the resulting panel for pairwise collinearity.

### User Inputs
`consensus_methods` (checkboxGroupInput, defaults to every method run so far), `consensus_min_methods` (default `min(2, n_available)`), one `consensus_w_<method>` numericInput per available method (weight, default 1), `consensus_use_weighted`/`consensus_min_weighted` (optional additional weighted-score gate), `consensus_run_btn`; then, in a second card, `corr_method` (Pearson/Spearman), `corr_threshold` (0.8), `use_corr_reduced` (checkbox), `corr_reduce_btn`.

### Input Data
`fs_method_ids()` (lines 1315-1323) — a named list of `selected_ids` character vectors, one per method that has actually been run (`uni_has_run()`, `reg_has_run()`, `rf_has_run()`, `rfe_has_run()`, `stab_has_run()`), so the Consensus tab dynamically reflects however many (1–5) methods have actually completed.

### Functions Used

#### `methyl_fs_consensus_table()` (lines 403-416)
Builds a `cpg × method` binary membership table over the union of all selected IDs, then `n_methods` (row-sum) and `weighted_score` (weighted row-sum ÷ total weight).

#### `methyl_fs_consensus_select()` (lines 418-422)
Keeps CpGs with `n_methods >= min_methods`, OR (if enabled) `weighted_score >= min_weighted`.

#### `methyl_fs_upset_plot()` (lines 427-458)
A hand-rolled UpSet-style plot (intersection-size bar chart + dot-matrix membership panel, combined via `patchwork::wrap_plots`) — explicitly built without an `UpSetR`/`ComplexUpset` dependency since at most 5 sets are ever involved (comment at line 424-426).

#### `methyl_fs_consensus_rank_plot()` / `methyl_fs_method_heatmap()` (lines 460-479)
Bar chart of the top-30-by-weighted-score CpGs; a method-membership heatmap for the top 60 CpGs (via `stats::reshape()` to long format).

#### `draw_overlap_venn()` (`global.R:1884`+, shared helper)
**Package:** `ggVennDiagram`.
Used for the Venn diagram (2–7 sets); a shared plotting helper reused across several modules in the app (candidates, preprocessing), not written specifically for this one.

#### `methyl_fs_correlation_reduce()` (lines 484-504)
**Purpose:** greedy pairwise-correlation pruning of the (already small) consensus panel — repeatedly finds the highest-|r| pair above `r_threshold` and drops whichever member has the lower consensus `weighted_score`, until no remaining pair exceeds the threshold.
**Processing:** `stats::cor(t(mat[ids,]), method, use="pairwise.complete.obs")`.
**Output:** `list(reduced_ids, dropped)` — a log of every kept/dropped pair and their correlation.
**Audit:** correct greedy-elimination logic with a proper termination guard (`length(active) < 2`); explicitly scoped to only the small consensus panel (not run genome-wide), which is the only regime where an `O(n²)`-per-iteration greedy loop is tractable — correctly documented as such in-code.

### Output
Venn diagram, UpSet plot, consensus-rank bar chart, method-comparison heatmap, the full intersection table (CSV download), and — if run — the correlation-reduction summary and dropped-pairs table.

### Scientific Interpretation
`n_methods`/`weighted_score` quantify cross-method agreement, not statistical significance or effect size — a CpG selected by all five methods is more *robustly* associated with the phenotype under this workflow's specific method choices, not necessarily more *strongly* associated in a magnitude sense. Correlation reduction addresses redundancy (two CpGs telling the same story) but does not address confounding or causality.

### Code Audit
- **Correct / appropriate:** consensus logic is simple, transparent, and fully traceable (a plain weighted union-count, no hidden heuristics); correlation reduction is correctly scoped to the small post-consensus panel only.
- **Minor — "consensus" from a single method is technically reachable.** `consensus_min_methods` defaults to `min(2, length(avail))`, so if only one method has been run, the default requires only `n_methods >= 1`, meaning the "consensus panel" is then just that one method's own selected list. This is a defensible design choice (nothing stops a user from wanting to preview a single method's panel through the Consensus/Selected-Features UI) but should be interpreted with that caveat — "Consensus" in the strict sense implies ≥2 independent methods agreeing.

---

## 14. Tab 8 — Selected Features

### Purpose
Presents the final CpG panel — post-consensus, and post-correlation-reduction if the user opted into it — as one interactive, annotated, exportable table, with a per-CpG drill-down.

### User Inputs
None beyond selecting a table row (`input$selected_table_rows_selected`, `DT` single-row selection) and the three download buttons (CSV/TSV/copy-as-TXT).

### Input Data
`fs_final_panel()` (lines 1472-1477): the correlation-reduced ID list if `corr_has_run() && input$use_corr_reduced`, else the raw consensus-selected IDs; joined against `methyl_fs_annotate_panel()`'s chr/pos/gene columns (from `methyl_get_annotation()`, `annotation.R:48`+ — only populated for array types with a bundled Bioconductor manifest, i.e. 450K/EPIC; every other array type gets `NA` annotation columns, correctly disclosed by the function's own construction rather than guessed).

### Processing
`fs_selected_table_full()` (lines 1495-1500) merges the consensus table (per-method flags, `n_methods`, `weighted_score`) with the annotation data.frame by `cpg`.

### Output
The annotated, filterable/sortable `DT` table; on row selection, a per-CpG detail panel showing which methods selected it, its consensus score, a violin+boxplot of its beta value by group (`ggplot`, `theme_arthomix()`), and a small per-sample table; CSV/TSV/plain-text-ID downloads.

### Scientific Interpretation
This is the panel the rest of the workflow (validation, export) treats as "the answer" — the practical, reduced, biologically inspectable output of the whole submodule. Note explicitly what it does **not** mean (see §17 below): panel membership is agreement-based evidence of statistical/ML association, not proof of biological mechanism, causality, or clinical utility.

### Code Audit
- **Correct / appropriate:** row-selection drill-down is defensively guarded (`req(cpg %in% rownames(r$beta))`, lines 1533, 1545) against the case where a user changes Data & Filters after having already selected a row, avoiding a stale-selection crash.
- No issues identified beyond those already covered in Tab 7's audit (this tab is a thin, correctly-guarded presentation layer over `fs_final_panel()`).

---

## 15. Tab 9 — Model & Export

### Purpose
Estimates how well the final panel actually generalizes (two validation strategies of differing rigor), and packages the entire run — not just the CpG list — into one reproducible artifact.

### User Inputs
`validate_mode` (Frozen panel / Leakage-safe), `validate_classifier` (Logistic/RF/SVM), `validate_k` (5), `validate_repeats` (1), `validate_run_btn`; then a Save-as-RDS button plus three CSV/TXT export buttons; then a Load-Previous-RDS `fileInput`.

### Input Data
`fs_filter_result()` and `fs_final_panel()`.

### Functions Used

#### `methyl_fs_validate_frozen()` (lines 531-543)
**Package:** `caret` (`caret::trainControl`, `caret::train`), `pROC`.
**Processing:** standard k-fold (or repeated-k-fold) CV of a classifier (`glm`/`rf`/`svmLinear`) trained on the panel's M-values, `metric="ROC"`, `savePredictions="final"`; AUC computed via `pROC::roc()`/`pROC::auc()` on the pooled out-of-fold predictions.
**Audit:** the UI itself explicitly and correctly discloses this mode's optimism (line 1580: *"the panel that was already chosen using this same data, so the reported AUC is optimistic"*) — a rare and valuable example of the code self-disclosing a leakage risk in its own interface text rather than leaving it implicit.

#### `methyl_fs_validate_nested()` (lines 549-584)
**Purpose:** a partially leakage-safe alternative — reselects a panel *inside every outer CV fold*, using only that fold's training data, rather than validating the one fixed panel chosen from the full dataset.
**Processing:** for each outer fold, refits Univariate selection (fixed to `moderated_t`, `top_n=100`) and LASSO (`alpha=1`) on the training split only, then trains the chosen classifier with `trainControl(method="none")` on that fold-specific panel and evaluates on the held-out fold.
**Audit — High, code-evidenced partial leakage in the "leakage-safe" mode. [DISCLOSURE IMPROVED — see "Post-Audit Fixes Applied" above; the underlying residual leakage was judged too invasive to eliminate without a live test environment, so this remains a documented follow-up rather than a full fix.]** The function receives `beta_mat`/`m_mat` as **already-filtered, already-imputed, already top-variance-capped** matrices — `fs_filter_result()` (Stage 1) is computed exactly once, outside and before the outer-fold loop (`r <- fs_filter_result()` at line 1608, then `methyl_fs_validate_nested(r$beta, r$m, r$grp, ...)` at line 1611). This means every outer fold's "training" data was already influenced, before the split, by summary statistics (variance, missingness) computed across **all** samples, including whichever samples land in that fold's held-out test set, and by a single global per-probe median-imputation pass likewise computed across all samples. Only the *supervised* re-selection stage (Univariate + LASSO) is genuinely redone per fold — the *unsupervised* preprocessing stage is not. The UI's own disclosure of this mode's limits (line 1582: *"it does not rerun Tree-Based/RFE/Stability/Consensus per fold"*) is accurate as far as it goes but does not mention this additional, narrower leakage source from the shared Stage-1 matrix. In practice this is a much smaller leakage risk than the Frozen-panel mode's (unsupervised variance/missingness filtering is a milder form of information leak than reusing the exact same supervised-selected panel), but it means "Leakage-safe" is a relative, not an absolute, description of this validation mode, and that distinction is not made explicit to the user.
- **Related — the nested mode validates a different, re-derived panel, not literally the Selected-Features panel.** `methyl_fs_validate_nested()` never receives `fs_final_panel()$ids` as an argument; it always re-derives its own panel per fold via a fixed Univariate(top-100)+LASSO recipe, regardless of which methods the user actually ran or how many CpGs ended up in the real consensus panel. A user could reasonably read "Leakage-safe validation: mean AUC 0.81" as an honest estimate of *the exact panel shown in the Selected Features tab*, when it is in fact an estimate of a different, always-univariate-plus-LASSO-derived proxy panel of a different size. This should be clearly communicated in any UI copy or thesis description that references this feature.

### Output
Mean (± SD, nested mode) AUC, a per-fold results table, and — separately — a `.rds` bundle (`fs_model_export()`, lines 1638-1661: dataset info, filter cascade/notes, every method's own parameters and results, consensus config, correlation-reduction config, the final annotated panel, validation results, and session/package-version info for `glmnet`/`randomForest`/`caret`/`e1071`), plus standalone CSV/TXT exports (ranking table, consensus table, plain-text analysis summary) and a "Load Previous RDS Model" viewer (explicitly disclosed as reference/comparison only, never mutating the active session — verified correct: `fs_loaded_model()` writes to no shared reactive state).

### Scientific Interpretation
Mean AUC quantifies how well a classifier built on this panel discriminates the two groups on held-out data — the single most important number for judging whether the whole feature-selection exercise produced something useful, but its trustworthiness depends entirely on which validation mode was used and on the caveats above.

### Code Audit — additional findings for this tab
- **Minor — incomplete reproducibility metadata. [FIXED — see "Post-Audit Fixes Applied" above.]** `session_info$package_versions` in `fs_model_export()` (lines 1657-1659) recorded only `glmnet`/`randomForest`/`caret`/`e1071`, omitting `limma` and `pROC`, both of which are used in every single run (Univariate Selection and any validation, respectively) — a small gap in an otherwise thorough reproducibility bundle. Both are now included.
- **Correct / appropriate:** the export bundle is genuinely comprehensive relative to its own stated goal ("not just the CpG list" — card description, line 1588) and does include full per-method parameter/result payloads, matching that claim; the "Load Previous RDS" feature correctly validates the file's provenance (`identical(x$module, "mod_methyl_featureselection")`) before displaying anything, with a clean `validate()` message otherwise.

---

## Function-by-Function Code Audit

*(Every function already documented in its owning tab's section above is cross-referenced here rather than repeated verbatim; this section adds the remaining shared/helper functions not already covered tab-by-tab.)*

### `methyl_fs_parse_matrix_rds()` (lines 32-50)
**Source:** `mod_methyl_featureselection.R`. **Package:** base R (`readRDS`).
**Purpose:** `.rds` companion to `qc.R`'s CSV/TSV-only `methyl_parse_matrix()`, so upload parsing doesn't need to branch on file type at call sites.
**Processing:** accepts either a bare numeric matrix, or a data.frame with probe IDs in the first column (coerced, with a duplicate-ID check).
**Audit:** correct, matches the exact `list(ok, mat, error)` contract of its CSV/TSV sibling; the transcriptomics module's own `.rds` upload parser (`global.R:1196`, `tx_parse_expr_matrix_rds`) is explicitly documented as mirroring this same function — a deliberate, disclosed cross-module pattern reuse, not independent duplication drift.

### `methyl_fs_parse_matrix_any()` (lines 52-55)
Trivial dispatcher on file extension. Correct.

### `methyl_fs_detect_scale()` (lines 60-68)
**Purpose:** beta-vs-M-value heuristic for uploads only (never used for the preloaded dataset, whose scale is authoritative from `dataset$input_scale`).
**Processing:** samples ≤20,000 finite values, checks whether the 0.1–99.9% quantile sits inside `[-0.05, 1.05]`.
**Audit:** a reasonable, appropriately-labeled heuristic (never claims certainty — the result UI additionally warns if >0.1% of sampled values fall outside the expected beta range, line 726-729).

### `methyl_fs_guess_group_col()` (lines 70-73)
Best-effort column-name guess (`group`/`disease`/`status`/`phenotype`, case variants), falling back to the first column. Correct, low-risk (only sets a default `selectInput` value the user can immediately override).

### `methyl_fs_cap_top_variance()` (lines 75-80)
Thin wrapper around `methyl_row_vars()` (`qc.R:22-28`) selecting the top-N by variance. Correct; the resulting rows are re-sorted back to original probe order (`sort(top)`) so downstream displays aren't scrambled.

### `methyl_fs_filter_iqr()` (lines 85-91)
**Package:** `matrixStats::rowQuantiles`.
IQR companion to `qc.R`'s SD/variance filters, deliberately kept local to this module rather than added to the shared `qc.R` filter set, since IQR thresholding is specific to this workflow (documented in-code, lines 82-84). Correct.

### `methyl_fs_sample_missing_ok()` (lines 94-98)
Sample-axis missingness — correctly noted in-code as having no equivalent in `qc.R` (which is probe-axis only). Correct.

### `methyl_fs_impute()` (lines 104-114)
Already covered in §7's audit (silent k-NN→median degradation).

### `methyl_fs_annotate_panel()` (lines 514-525)
Already covered in §14. Correctly leaves CpG-island/genomic-region columns as `NA` rather than fabricating them, since the shared `annotation.R` only exposes chr/pos/Type/SNP-rs/gene (documented in-code, lines 510-513).

### Shared cross-file dependencies (not modified by this module, reused as-is)
| Function | File | Role here |
|---|---|---|
| `methyl_row_vars` | `qc.R:22-28` | variance computation (matrixStats-accelerated) underlying `methyl_filter_variance`/`methyl_filter_sd`/`methyl_fs_cap_top_variance` |
| `methyl_filter_missing`, `methyl_filter_variance`, `methyl_filter_sd` | `qc.R:30-52` | Stage 1 probe filters |
| `methyl_filter_snp`, `methyl_filter_sex_chr` | `qc.R:68-104` | optional annotation-dependent probe filters |
| `methyl_filter_cross_reactive`, `methyl_filter_maf`, `methyl_parse_maf_list` | `qc.R:112-160` | optional upload-dependent probe filters |
| `methyl_sheet_sample_ids` | `qc.R:456-461` | sample/sheet ID resolution |
| `methyl_probe_retention_cascade` | `qc.R:540-553` | cascade table for the filter-result plot |
| `methyl_beta_to_mvalue` | `qc.R:559-562` | beta→M-value logit transform |
| `methyl_plot_cascade` | `qc.R:830-838` | shared cascade bar-chart builder |
| `methyl_parse_matrix`, `methyl_parse_sample_sheet`, `methyl_parse_probe_list` | `parse_upload.R:11-57` | CSV/TSV upload parsing |
| `methyl_get_annotation` | `annotation.R:48`+ | manifest annotation (chr/pos/gene/SNP) for 450K/EPIC |
| `ARTHOMIX_COLORS`, `arthomix_pair`, `theme_arthomix`, `draw_overlap_venn` | `global.R:1417-1455`, `1884`+ | shared app-wide plot styling |
| `METH_DATA_AVAILABLE`, `METH_RAW_DATA_AVAILABLE`, `load_default_meth_matrix`/`_pheno` | `data_paths.R`, `global.R:340-390`+ | preloaded-dataset availability/loading |

None of these shared functions are modified by this module — confirming the file header's own claim (lines 15-19): "every filter/parsing/annotation/plotting helper below is reused, never modified."

---

## Shiny Reactive Logic

| Pattern | Where used | Input → Reactive calculation → Processing → Output |
|---|---|---|
| `eventReactive()` | `fs_filter_result`, `fs_uni_result`, `fs_reg_result`, `fs_rf_result`, `fs_rfe_result`, `fs_stab_result`, `fs_consensus_result`, `fs_corr_result`, `fs_validate_result`, `fs_loaded_model` | Each fires only on its own `actionButton` click, never on parameter change alone — the core mechanism implementing the "explicit Run per stage" design (§5) |
| `reactive()` | `fs_own_matrix`, `fs_own_sheet`, `fs_exclusion_ids`, `fs_maf_table`, `fs_active_source`, `anno_result`, `fs_group_choices`, `fs_available_methods`, `fs_method_ids`, `fs_final_panel`, `fs_selected_table_full`, `fs_model_export` | Always-live derived values (source resolution, dynamic option lists, panel assembly) — deliberately *not* Run-gated, since these are cheap/pure derivations rather than model fits (matches the file's own stated design principle, lines 626-629) |
| `reactiveVal()` | `filters_has_run`, `uni_has_run`, `reg_has_run`, `rf_has_run`, `rfe_has_run`, `stab_has_run`, `consensus_has_run`, `corr_has_run`, `validate_has_run` | One boolean gate per stage, toggled `TRUE` on that stage's Run button and reset `FALSE` when its upstream dependency changes — this is what makes every result card show `mod_methyl_fs_empty_note()`/`mod_methyl_fs_need_filters_note()` instead of stale content |
| `observeEvent()` | resets each `*_has_run` on `fs_filter_result()` change; sets each `*_has_run(TRUE)` on its Run button; shows/hides `shinyjs` panels (correlation-reduction result, CpG detail, loaded-model viewer) | Standard state-machine wiring; consistently uses `ignoreInit = TRUE` on the "set TRUE" observers so a fresh module load doesn't spuriously mark a stage as already run |
| `req()` | throughout (e.g. `fs_own_matrix`, drill-down plot/table guards) | Silently halts a reactive rather than erroring when a required input isn't yet available |
| `validate(need(...))` | throughout every fit/parsing function | Converts expected failure modes (too few samples, rank-deficient design, `glmnet`/`randomForest`/`caret::rfe` fit failures) into clean, in-app messages instead of raw R errors — see §14 for the one case (`factor()` with duplicate levels) this pattern does **not** cover |
| `renderUI()` | one per tab (`filters_ui`, `uni_ui`, `reg_ui`, `rf_ui`, `rfe_ui`, `stab_ui`, `consensus_ui`, `selected_ui`, `export_ui`) plus per-result sub-panels | Lazy tab construction — a tab's inputs don't exist until its `tabPanel` is actually opened |
| `downloadHandler()` | every CSV/TSV/TXT/RDS export across all 9 tabs | Straightforward `utils::write.csv`/`write.table`/`writeLines`/`saveRDS` content functions, no dynamic filename collisions observed |
| `outputOptions(..., suspendWhenHidden = FALSE)` | every `DT::renderDataTable` output | Ensures tables continue to render/update even while their tab isn't the active one — necessary because `DT` outputs otherwise suspend when hidden, which would break the "switch tabs, come back, data's still there" expectation |
| Module namespace | `NS(id)` / `ns()` throughout; server invoked as `paste0("mx_", m$config$id)` (`server.R:95`) | Standard Shiny module namespacing — no cross-module ID collisions identified |

---

## Normalisation and Preprocessing Audit

This module performs **no independent normalization of its own** — normalization (dasen/BMIQ/noob/funnorm/etc.) is the responsibility of the upstream Methylomics Normalization submodule (`mod_methyl_normalization.R`), and this module's preloaded-cohort route consumes an already-QC'd, already-normalized beta matrix (`dataset$beta`, sourced from the pipeline documented in `METHODS_load_qc.md`). What this module *does* perform, in execution order (§7):

| Operation | What / Why | Before or after selection? | Samples or features? | Leakage risk? | Appropriate for methylomics? | Matches UI label? | Reproducible? | Train/eval-consistent? |
|---|---|---|---|---|---|---|---|---|
| Sample-missingness filter | Drop samples with too many missing probes | Before | Samples | None (unsupervised, and applied identically regardless of eventual CV) | Yes | Yes | Yes (deterministic) | Yes — same matrix feeds every method |
| Probe missingness filter | Drop probes missing in too many samples | Before | Features | **See below** | Yes | Yes | Yes | Yes |
| Variance/SD/IQR filter | Drop near-constant probes | Before | Features | **See below** | Yes, standard practice | Yes | Yes | Yes |
| SNP / sex-chr / cross-reactive / MAF filters | Remove technically unreliable probe classes | Before | Features | None (unsupervised, annotation- or list-based) | Yes | Yes | Yes (deterministic given inputs) | Yes |
| Top-variance cap | Bound candidate panel size for runtime | Before | Features | **See below** | Reasonable engineering trade-off | Yes | Yes | Yes |
| Imputation (median or k-NN) | Fill remaining missing beta values | Before | Features (per-probe stat) | **See below** | Standard practice for these methods | Partially — silently degrades k-NN→median with no notice | Yes given fixed method | Yes |
| Beta → M-value transform | Logit transform for statistically better-behaved modeling scale | Before every method needing M-values | Values (elementwise) | None | Yes — matches standard practice (`minfi::logit2`/`lumi::beta2m` convention, per in-code comment) | Yes | Yes | Yes |

**On data leakage specifically:** every filtering/imputation operation above is *unsupervised* — none of them look at the group label `y` to decide which probes or samples to keep, so there is no leakage of the phenotype itself into feature selection. However, as detailed in §15's audit of `methyl_fs_validate_nested()`, **all of the above operations are computed once, globally, before the outer cross-validation loop even in the "Leakage-safe" validation mode** — meaning summary statistics (missingness fractions, variances, per-probe medians) are computed across samples that will later be split into training and held-out test folds. This is a real, code-verified partial leakage source, distinct from and milder than the (correctly self-disclosed) full leakage of the Frozen-panel validation mode, but neither the code nor its UI copy currently distinguishes "no supervised leakage" from "no leakage at all" — the latter is not achieved by either validation path.

---

## Machine-Learning Methodology Audit

| Method | Type | Supervised? | Response | Predictors | Hyperparameters (defaults) | Training/validation | Feature-importance mechanism | Threshold | Seed | Output |
|---|---|---|---|---|---|---|---|---|---|---|
| Univariate (6 sub-methods) | Statistical test, per-feature | Yes | 2-level group (or numerically-coded proxy for the 3 "continuous" options) | 1 CpG at a time | test choice, optional covariates | none (marginal test, no CV) | test statistic / p-value | FDR / p / Top-N | N/A (deterministic) | ranked table + selected IDs |
| LASSO / Elastic Net | Regularized logistic regression | Yes | 2-level group | all candidate CpGs jointly | alpha (1), lambda selection, nfolds (10), class weighting | internal `cv.glmnet` k-fold | non-zero coefficient magnitude | coefficient ≠ 0 (+ optional magnitude floor, + optional cap) | 1234 | selected IDs + coefficients |
| Random Forest | Ensemble tree classifier | Yes | 2-level group | all candidate CpGs jointly | ntree (1000), mtry (√p), nodesize (1), maxnodes (unlimited) | out-of-bag internally (no explicit CV split shown) | Gini or permutation-accuracy importance | Top-N / threshold / percentile | 1234 | ranked table + selected IDs |
| RFE (RF/Logistic/SVM) | Wrapper, iterative elimination | Yes | 2-level group | all candidate CpGs, shrinking | subset sizes, CV folds (5), (SVM: cost, tolerance) | `caret::rfe` internal CV (RF/Logistic) or manual k-fold CV sweep (SVM) | model-specific ranking + CV performance curve | empirical CV-error minimum (or manual) | 1234 | selected IDs + optimal size + curve |
| Stability Selection | Resampled LASSO ensemble | Yes | 2-level group | all candidate CpGs jointly | resampling scheme, n_resamples (50), fraction/k/repeats | resampling, no held-out test per se | selection frequency across resamples | frequency ≥ threshold (0.7) | 1234 | ranked table + selected IDs |
| Consensus | Rule-based aggregation | N/A | — | method membership | min_methods, weights, optional weighted floor | none | count/weighted overlap | ≥N methods (and/or weighted score) | N/A | consensus table + panel |
| Correlation reduction | Greedy redundancy pruning | No (unsupervised, post-hoc) | — | consensus panel only | correlation method, threshold (0.8) | none | pairwise |r| vs. consensus score | |r| > threshold | N/A | reduced panel |
| Validation (Frozen) | Classifier CV | Yes | 2-level group | final panel | classifier, k (5), repeats (1) | k-fold/repeated-k-fold CV via `caret::train` | N/A | N/A | inherits stage seeds; validation itself unseeded per-fold beyond `caret`'s own | mean AUC + per-fold table |
| Validation (Nested) | Classifier CV with per-fold reselection | Yes | 2-level group | re-derived (Uni+LASSO) panel per fold | outer_k (5), repeats (1), classifier | outer k-fold, inner reselection | N/A | N/A | fixed seed reused across folds (documented as intentional for a fixed-hyperparameter design) | mean ± SD AUC + per-fold table |

**Explicit checks against common pitfalls, evaluated against the code (not assumed):**
- **Leakage:** partially present, as detailed above — unsupervised Stage-1 preprocessing is never redone inside either validation mode's fold structure.
- **Overfitting:** mitigated by CV-based lambda selection (LASSO), CV-based RFE subset sizing, and out-of-fold AUC estimation (both validation modes) — standard, appropriate safeguards.
- **Class imbalance:** partially addressed — LASSO offers a reachable "balanced" class-weighting option; Random Forest's `classwt` exists in the helper function but is not exposed in the UI (§10, Minor finding).
- **Insufficient sample size:** explicitly guarded at multiple points in Stage 1 (≥10 total, ≥3 per group, ≥6 post-sample-filter) with clean `validate()` messages.
- **Inappropriate cross-validation:** fold counts are defensively clamped to the smaller class's size everywhere folds are used (LASSO, stability, validation) — appropriate.
- **Feature-selection leakage:** addressed by the nested validation mode for the *supervised* stage, not for the *unsupervised* Stage-1 filtering stage (see above).
- **Training/evaluation contamination:** the Frozen mode is contamination by design (self-disclosed); the Nested mode has the narrower, code-verified contamination described above.
- **Unstable rankings:** directly the subject of the Stability Selection method, which exists specifically to quantify this.
- **Arbitrary thresholds:** most thresholds are user-configurable rather than hard-coded, with sensible defaults (FDR 0.05, correlation 0.8, stability frequency 0.7); RFE's subset-size choice is explicitly data-driven (CV-error minimization) rather than arbitrary.
- **Missing seed:** every stochastic method exposes and uses a seed input; the purely deterministic methods (Univariate) correctly have no seed input at all.
- **Incorrect response encoding:** the 2-level `grp` factor is consistently constructed with `levels = c(ref_group, comp_group)` and reused identically across every method — no inconsistent encoding was found between methods.
- **Incorrect feature orientation:** every method-facing helper transposes to the `sample × feature` orientation `glmnet`/`randomForest`/`caret`/`e1071` expect (`X <- t(r$m)`), consistently across Regularization/Tree-Based/RFE/Stability — no orientation bugs identified.

---

## Tab-to-Tab Data Flow

```text
Data & Filters (Tab 1)
   |  writes: fs_filter_result()  {beta, m, grp, cascade, ...}
   |  (Run-gated: filters_run_btn; invalidates every downstream *_has_run)
   v
   +-----------------------------------------------------------------+
   |         (each of the next five tabs consumes fs_filter_result() |
   |          independently and in any order; none of them depend    |
   |          on each other)                                         |
   v
Univariate Selection (Tab 2)  --\
LASSO (Tab 3)                   \
Tree-Based Selection (Tab 4)     |--> fs_method_ids()  (named list of
RFE / Wrapper Selection (Tab 5) /       selected_ids, only for methods
Stability Selection (Tab 6)  --/        that have actually been Run)
   |
   v
Consensus / Overlap (Tab 7)
   |  writes: fs_consensus_result()  {table, methods, weights, selected_ids}
   |  side effect: results$featureselection  <- summary  (shared, currently unread)
   |  optional: Correlation reduction  ->  fs_corr_result()
   v
Selected Features (Tab 8)
   |  reads: fs_final_panel()  = corr-reduced IDs (if opted in) else consensus IDs
   v
Model & Export (Tab 9)
   |  Validation reads: fs_final_panel() (Frozen mode) OR re-derives its own
   |                     panel independently from fs_filter_result() (Nested mode)
   |  Export bundles: fs_filter_result() + every method's own result objects +
   |                   fs_consensus_result() + fs_corr_result() + fs_final_panel() +
   |                   fs_validate_result()
   v
  .rds / CSV / TXT downloads
```

Every tab operates on **exactly one** shared `fs_filter_result()` snapshot at a time (Tabs 2–6 are otherwise fully independent of each other — none reads another's output), and every downstream `*_has_run` flag resets whenever that upstream snapshot changes, so stale cross-tab combinations cannot silently persist after a Data & Filters re-run. The one asymmetry in the diagram is Tab 9's Nested validation mode, which bypasses the Tabs 2–8 pipeline entirely for its own internal per-fold panel selection (§15).

---

## End-to-End ML Feature Selection Pipeline

```text
Preloaded cohort (GSE42861, 689 samples, 412,492 CpGs, beta, 450K)
        or
Uploaded beta/M-value matrix (+ optional phenotype/exclusion/MAF files)
        ↓  (fs_active_source, lines 659-673)
Scale detection / declaration (beta vs. M-value)
        ↓
Sample-ID resolution + reference/comparison group subsetting  (fs_filter_result, 802-816)
        ↓
Sample-missingness filter  (821-823)
        ↓
Probe filters: missingness, variance/SD/IQR, [SNP, sex-chr, cross-reactive, MAF]  (825-843)
        ↓
Top-variance cap  (851-856)
        ↓
Imputation (median or k-NN)  (858-863)
        ↓
Beta → M-value derivation  (865)
        ↓
   +---------------------------------------------------------------+
   |  Univariate  |  LASSO  |  Random Forest  |  RFE  |  Stability  |
   |  (independent, each explicitly Run, each producing selected_ids)
   +---------------------------------------------------------------+
        ↓
Consensus / overlap table (weighted union-count across whichever methods ran)
        ↓
[optional] Greedy pairwise-correlation redundancy reduction
        ↓
Final CpG panel (Selected Features tab)
        ↓
   +-------------------------------------------+
   |  Frozen-panel CV (optimistic, disclosed)   |
   |    or                                      |
   |  Nested CV (re-selects Uni+LASSO per fold, |
   |    partial leakage from Stage-1 preprocessing) |
   +-------------------------------------------+
        ↓
Reproducible export: .rds bundle (full config + results) + CSV/TSV/TXT exports
```

---

## UI-to-Code Mapping

| UI Component | Input/Output ID | Code Location | Function/Reactive | Purpose |
|---|---|---|---|---|
| Data source radio | `fs_source` | 691-693 | `fs_active_source()` | Choose preloaded vs. upload |
| Matrix upload | `fs_matrix_file` | 698 | `fs_own_matrix()` | Load a beta/M-value matrix |
| Phenotype upload | `fs_sheet_file` | 704 | `fs_own_sheet()` | Load sample metadata |
| Exclusion-list upload | `fs_exclusion_file` | 706 | `fs_exclusion_ids()` | Cross-reactive filter input |
| MAF-table upload | `fs_maf_file` | 708 | `fs_maf_table()` | MAF filter input |
| "Run Filters" button | `filters_run_btn` | 782 | `fs_filter_result()` | Build the shared filtered matrix |
| Filter result cascade plot | `fs_cascade_plot` | 893 | `methyl_plot_cascade()` | Visualize probe retention |
| Univariate method dropdown | `uni_method` | 905-910 | `fs_uni_result()` | Choose test |
| "Run Univariate Selection" | `uni_run_btn` | 925 | `fs_uni_result()` | Execute test |
| Univariate table + download | `uni_table` / `uni_download` | 974-984 | `DT::renderDataTable` / `downloadHandler` | Display/export ranked CpGs |
| LASSO alpha slider | `reg_alpha` | 996 | `fs_reg_result()` | LASSO↔Ridge mixing |
| "Run LASSO / Elastic Net" | `reg_run_btn` | 1011 | `fs_reg_result()` | Execute `cv.glmnet` |
| RF importance-type radio | `rf_importance_type` | 1082 | `fs_rf_result()` | Gini vs. Accuracy |
| "Run Random Forest" | `rf_run_btn` | 1088 | `fs_rf_result()` | Execute `randomForest` |
| RFE flavor selector | `rfe_flavor` | 1155 | `fs_rfe_result()` | RF/Logistic/SVM RFE |
| "Run RFE" | `rfe_run_btn` | 1168 | `fs_rfe_result()` | Execute chosen RFE flavor |
| Stability resampling-scheme radio | `stab_type` | 1241 | `fs_stab_result()` | Bootstrap/k-fold/subsampling |
| "Run Stability Selection" | `stab_run_btn` | 1253 | `fs_stab_result()` | Execute resampling loop |
| Consensus method checkboxes | `consensus_methods` | 1335 | `fs_consensus_result()` | Which methods count |
| "Run Consensus Selection" | `consensus_run_btn` | 1344 | `fs_consensus_result()` | Build consensus table; also sets `results$featureselection` |
| Venn / UpSet / rank / heatmap plots | `consensus_venn`/`consensus_upset`/`consensus_rank_plot`/`consensus_heatmap` | 1416-1429 | `draw_overlap_venn()` / `methyl_fs_upset_plot()` / `methyl_fs_consensus_rank_plot()` / `methyl_fs_method_heatmap()` | Visualize overlap |
| "Reduce for collinearity" | `corr_reduce_btn` | 1410 | `fs_corr_result()` | Greedy correlation pruning |
| Selected-features table | `selected_table` | 1502-1507 | `fs_selected_table_full()` | Final panel display |
| Row-click CpG detail | `selected_table_rows_selected` | 1509-1549 | `cpg_detail_plot`/`cpg_detail_table` | Per-CpG beta distribution |
| CSV/TSV/TXT downloads (final panel) | `selected_download_csv`/`_tsv`/`selected_copy` | 1551-1562 | `downloadHandler` | Export the panel |
| Validation-mode radio | `validate_mode` | 1574 | `fs_validate_result()` | Frozen vs. Nested |
| "Run Validation" | `validate_run_btn` | 1583 | `fs_validate_result()` | Execute chosen validation |
| "Save Model as RDS" | `save_rds` | 1589 | `fs_model_export()` | Full-configuration export |
| "Load Previous RDS Model" | `load_model_file` | 1597 | `fs_loaded_model()` | Reference-only comparison view |

---

## Input/Output Data Analysis

| Tab | Input Data | Input Dimensions | Main Function/Algorithm | Processing | Output Data | Biological/Statistical Meaning |
|---|---|---|---|---|---|---|
| Data & Filters | raw beta/M matrix + phenotype sheet | `p_raw × n_raw` probes/samples (dataset-dependent) | filter cascade + imputation | missingness/variance/QC filtering, top-variance cap, imputation, M-value derivation | filtered `beta`/`m` matrices, `grp` factor | Establishes the analysis-ready feature space and phenotype contrast |
| Univariate Selection | filtered `m`/`beta`, `grp` | `p_filtered × n_filtered` | `limma` fit or rank test | per-probe marginal test | ranked table + selected CpG IDs | Marginal group association per CpG |
| LASSO | filtered `m` (transposed), `grp` | `n_filtered × p_filtered` | `glmnet::cv.glmnet` | penalized joint logistic fit | selected CpG IDs + coefficients | Non-redundant, jointly predictive CpGs |
| Tree-Based Selection | filtered `m` (transposed), `grp` | `n_filtered × p_filtered` | `randomForest` | ensemble tree fit + importance | ranked table + selected CpG IDs | Non-linear/interaction-aware importance |
| RFE / Wrapper Selection | filtered `m` (transposed), `grp` | `n_filtered × p_filtered` (shrinking) | `caret::rfe` or custom SVM-RFE | iterative elimination + CV curve | selected CpG IDs + optimal size | CV-error-minimizing panel size for the wrapped model |
| Stability Selection | filtered `m` (transposed), `grp` | `n_filtered × p_filtered`, resampled | resampled fixed-lambda `glmnet` | repeated resample + tabulate | ranked table (selection frequency) + selected CpG IDs | Robustness of LASSO selection to resampling |
| Consensus / Overlap | up to 5 `selected_ids` lists | union of selected CpGs (≤ candidate panel size) | weighted union-count | overlap table + optional correlation pruning | consensus table + panel | Cross-method agreement |
| Selected Features | consensus/correlation-reduced panel + annotation | final panel size (typically tens of CpGs) | table join | annotate with chr/pos/gene | annotated final panel | The biomarker candidate panel |
| Model & Export | filtered matrices + final panel | panel size × `n_filtered` | `caret::train` (+ nested reselection) | CV classifier fit/eval | mean AUC + per-fold table; `.rds`/CSV/TXT exports | Generalization estimate of the panel; full reproducible record |

**Sample/feature vocabulary used consistently throughout:** samples = array/sequencing samples (columns); CpG sites/probes = features (rows); methylation values = beta ([0,1]) or M-value (logit-transformed, unbounded); phenotype/class labels = the two-level `grp` factor derived from the user-chosen phenotype column and reference/comparison levels; metadata = the full phenotype sheet, of which the group column is one field; selected features = the per-method `selected_ids`; model-derived statistics = p/FDR/statistic (univariate), coefficient (LASSO), importance (RF), rank/optimal_size (RFE), selection_frequency (Stability), n_methods/weighted_score (Consensus), AUC (Validation).

---

## Error Handling and Edge Cases

| Scenario | Code behavior |
|---|---|
| No data loaded | `fs_active_source()` calls `validate(need(!is.null(dataset$beta), ...))` for the preloaded route; upload route's `req(input$fs_matrix_file)` silently halts until a file is provided — clean either way |
| Missing values in matrix | Handled by the missingness filters + optional imputation (§7); if imputation is disabled and NAs remain, downstream `limma`/`glmnet`/`randomForest`/`caret` calls would receive `NA`s — **not explicitly tested/guarded beyond the filter thresholds themselves**, since M-values are always derived from the (by-default-imputed) beta matrix |
| Too few samples | Explicit `validate(need(...))` guards at ≥10 total, ≥3/group, ≥6 post-sample-filter (§7) |
| Too few features | `validate(need(sum(keep_probe) >= 10, ...))` with an actionable message (line 841-842) |
| Only one class exists | Multiple guards: `nlevels(yf) >= 2` (univariate), and Stage 1's ≥3-per-group check structurally prevents a single-class contrast from ever reaching the selection methods — **except** the identical-ref/comp-group edge case below, which bypasses this guard entirely |
| Selected feature count exceeds available features | LASSO's `max_selected` and RF's `top_n` are both clamped (`min(p, ...)`) — cannot request more than exist |
| Invalid thresholds entered | Numeric inputs generally have `min`/`max`/`step` bounds in the UI; no server-side re-validation beyond what the fitting functions themselves reject via `validate()` |
| An ML model fails to fit | Every `glmnet`/`randomForest`/`caret::rfe`/`caret::train` call is wrapped in `tryCatch(..., error = function(e) validate(need(FALSE, paste(..., conditionMessage(e)))))` — converts to a clean in-app message, not a crash |
| A required column is absent | `fs_group_levels_ui` guards with `req(input$fs_group_col %in% colnames(src$sheet))`; covariate selection (`uni_covariate_ui`) only offers columns that actually exist |
| Phenotype variable missing entirely | `fs_filter_result()`'s first `validate(need(!is.null(sheet) && !is.null(input$fs_group_col), ...))` catches this cleanly |
| User changes inputs after results generated | Each stage's own parameters can be changed freely without invalidating its last result until its own Run button is clicked again (intentional, §5); changing **upstream** Data & Filters parameters does correctly invalidate every downstream stage |
| Empty result produced | Every selection function's `selected_ids` can legitimately be an empty character vector (e.g. an overly strict FDR threshold) — downstream tabs (`fs_available_methods`, Consensus) correctly treat a method with zero selections as still "available" (it ran), and `methyl_fs_consensus_table()` handles an all-empty `id_lists` case explicitly (line 406) |
| **Identical reference and comparison group selected** | **Not guarded.** `factor(grp_raw[keep_s], levels = c(ref, comp))` with `ref == comp` throws R's built-in `"duplicated levels in factor are not allowed"` error, which is **not** wrapped in `tryCatch`/`validate()` at that call site (line 815) — surfaces as Shiny's generic uncaught-error screen rather than this module's usual clean inline message. See §7's audit for full detail. |
| Consensus with only one method run | Reachable by design (`consensus_min_methods` defaults to `min(2, n_available)` = 1 when only one method ran) — produces a "consensus" that is really just that one method's own list; not blocked, and not mislabeled loudly enough for a casual user to notice (§13) |

---

## Quality-Control Audit

### A. Correct / Appropriate
- Explicit Run-button gating per stage with correctly cascading invalidation, avoiding both unwanted recomputation and stale cross-tab results.
- Every probe/sample filter returns an explicit `keep` mask *and* a human-readable `note`, making the retention cascade fully auditable.
- Consistent, correct `sample × feature` transposition for every `glmnet`/`randomForest`/`caret`/`e1071` call.
- Fold counts defensively clamped to the smaller class's size everywhere CV is used.
- Two independently reasoned validation modes, with the Frozen mode's optimism *explicitly self-disclosed in the UI text* — a genuinely good practice.
- Comprehensive, well-structured `.rds` export bundle that lives up to its "not just the CpG list" claim.
- `tryCatch`-wrapped model fits converted to clean `validate()` messages throughout, rather than raw crashes, with the one documented exception below.
- No modification of any shared/reused helper function from `qc.R`/`parse_upload.R`/`annotation.R`/`global.R` — the module's own claim of isolation (file header, lines 15-19) is verified true by this reading.

### B. Minor Issues
- k-NN imputation silently degrades to median with no user-facing notice when the `impute` package isn't installed.
- Silent truncation of out-of-range RFE subset sizes with no explanatory note (inconsistent with Stage 1's always-noted filters).
- `glmnet`'s raw base-R CV plot breaks visual consistency with the rest of the module's `ggplot2`/`theme_arthomix()` styling.
- `session_info$package_versions` in the export bundle omits `limma` and `pROC`, both used in every run.
- "ANOVA (>2 groups)" and "Kruskal-Wallis (>2 groups)" labels describe functionality that the binary-only Data & Filters contrast can never actually exercise beyond 2 groups.

### C. Moderate Issues
- Random Forest's `classwt`/`sample_fraction` class-imbalance parameters exist in the helper function but are not reachable from the UI.
- Stability Selection's `selection_frequency`/`stability_score` are silently deflated by degenerate (single-class or failed-fit) resamples, whose all-zero columns still count toward the reported denominator.
- "Stability Selection" as implemented is a documented simplification (single fixed reference lambda) of the term's stricter meaning in the Meinshausen–Bühlmann literature — worth flagging explicitly for a thesis reader who may expect the full lambda-path formulation.
- The nested "Leakage-safe" validation mode validates a re-derived Univariate+LASSO panel, not literally the panel shown in Selected Features — a real risk of over-interpreting the reported AUC as describing the actual final panel.

### D. High-Severity Issues
- **[FIXED] Selecting the same reference and comparison group crashes the reactive with an uncaught R error** (`factor()` duplicate-levels error), rather than the clean `validate()` message this module uses everywhere else. Concrete failure scenario: a user with a phenotype column containing levels `c("Control","RA")` opens the Reference/Comparison dropdowns and, whether by misclick or by a column with only one distinct value surviving upstream filtering, selects `"RA"` for both — the app then throws a raw error the first time "Run Filters" is clicked. A `validate(need(...))` guard now catches this with a clean message.
- **[FIXED] The "ANOVA (>2 groups)" option is functionally a two-group test**, because the enclosing workflow structurally enforces a binary contrast one tab upstream; the F-test code path it is meant to exercise is provably unreachable given the current Data & Filters design. This is a genuine correctness-of-labeling issue: the exported CSV/RDS would record `"anova"` as the method used, while the actual statistic computed is the same t-based statistic as `"t_test"`/`"moderated_t"`. ANOVA/Kruskal-Wallis now genuinely test every level of the chosen phenotype column when more than two are present among the eligible samples, verified with synthetic multi-group data.
- **[DISCLOSURE IMPROVED, not fully fixed] Partial data leakage persists in the "Leakage-safe" nested validation mode**, sourced from Stage-1's global (pre-fold-split) unsupervised filtering and imputation. This does not invalidate the mode's value relative to the Frozen mode (which has a much larger, self-disclosed leakage source), but the "Leakage-safe" name is not a fully accurate description of what the code does. The UI copy now explicitly says so; eliminating the residual leakage itself would require re-plumbing Stage 1's filtering to run per outer fold, judged too invasive to attempt without a live test environment and left as documented follow-up work.

Also fixed in this pass, beyond the items already flagged inline above: the "Preloaded whole-blood cohort" data-source option was gated on the wrong availability flag (`METH_DATA_AVAILABLE`, the lightweight-tables flag, instead of `METH_RAW_DATA_AVAILABLE`, the flag that actually governs whether the raw matrix can load) — a deployment with the pipeline tables but not the raw matrix would have offered a "Preloaded" choice that could never work. See "Post-Audit Fixes Applied" at the top of this document for the complete list and how each was verified.

---

## Teaching Notes: How to Understand This Module

**What is a feature?** In this module, a feature is one CpG site — one location in the genome where DNA methylation is measured, typically reported as a beta value between 0 (fully unmethylated) and 1 (fully methylated).

**What is a CpG?** A cytosine-guanine dinucleotide where DNA methylation predominantly occurs in mammalian genomes. Methylation arrays (like the Illumina 450K/EPIC platforms this app is built around) measure the methylation level at hundreds of thousands of specific CpG sites simultaneously.

**What is a predictor?** In the ML sense used throughout this module, a predictor is a feature (a CpG's methylation value) used as an input to a statistical test or machine-learning model that tries to explain or predict an outcome.

**What is the outcome?** Here, the outcome is always the two-level group label built from the user's chosen phenotype column and the two selected levels (e.g. "RA" vs. "Control") — every method in this module is a *supervised* method with respect to this same binary outcome (Consensus and correlation reduction are the two purely unsupervised, post-hoc exceptions).

**What is feature selection?** The process of choosing a subset of features (here, CpGs) that are most relevant to the outcome, discarding the rest — as distinct from *feature extraction* (which would combine features into new derived variables, e.g. principal components) or *dimensionality reduction* more generally. This module performs pure feature *selection*: every CpG in the final panel is one of the original measured CpGs, never a synthetic combination.

**Why is feature selection important for methylomics specifically?** Because methylation arrays measure far more CpGs than most cohorts have samples, and many CpGs are highly correlated with their neighbors — see §1's fuller explanation.

**Difference between feature selection and prediction.** Feature selection asks "which variables matter?"; prediction asks "how accurately can I forecast the outcome for a new sample?" These are related but distinct goals — a panel can be excellent for prediction (small, high-AUC) while individual CpGs within it are not necessarily the most *biologically* important ones, and vice versa. This module's five selection methods target the *selection* question directly; its Model & Export tab's validation step is the module's one nod to the *prediction* question, used here purely to sanity-check the panel's usefulness, not as an end in itself.

**Difference between feature importance and statistical significance.** A univariate p-value (Tab 2) answers "is this CpG's group difference bigger than expected by chance alone?" — a significance question. A LASSO coefficient or Random Forest importance score (Tabs 3–4) answers "how much does this CpG contribute to a joint/ensemble model's ability to separate the groups?" — an importance question. A CpG can be statistically significant but redundant (importance near zero once correlated CpGs are already in the model), or important to a model without being significant on its own (e.g. only informative in combination with another CpG). This is precisely why the module runs both kinds of methods and reconciles them via Consensus rather than trusting either alone.

**Why preprocessing matters.** Raw methylation data contains missing values, near-constant probes, and (optionally) technically unreliable probe classes (SNP-overlapping, cross-reactive, sex-linked). Feeding all of that directly into a statistical or ML method wastes computation on uninformative features and can introduce spurious associations from technical artifacts rather than biology — hence Stage 1's filtering cascade before any selection method runs.

**Why validation matters.** A feature-selection method chosen and evaluated on the exact same data it selected from will always look better than it truly is — the model has, in effect, already "seen the answer." The Model & Export tab's two validation modes exist to quantify (Frozen) or partially correct for (Nested) this optimism, so that the reported AUC means something closer to "how well would this panel work on a *new* sample" rather than "how well does this panel fit the data it was chosen from."

**Why leakage is dangerous.** Leakage — letting information from what should be held-out test data influence model/feature choices made on the training data — silently inflates every performance metric computed afterward, in a way that is invisible from the metric alone (an 0.95 AUC computed with leakage looks identical to a genuine 0.95 AUC until tested on truly independent data). This is exactly why this module's own UI text explicitly warns about the Frozen-panel mode's optimism, and why this audit's most important finding (§15) is that even the Nested mode is not fully leakage-free.

---

## Scientific Interpretation of Results

### Computational meaning
Each method's output is a numeric score (p-value, FDR, coefficient, importance, selection frequency) computed by a specific, traceable algorithm on the module's shared filtered matrix, and a resulting `selected_ids` list built by thresholding that score.

### Statistical meaning
- Univariate p/FDR: probability of observing this CpG's group difference (or more extreme) under the null hypothesis of no association, corrected for multiple testing where FDR is used.
- LASSO coefficient: the CpG's independent contribution to a penalized joint logistic model at the CV-selected regularization strength.
- RF importance: contribution to classification split quality / accuracy across the forest.
- RFE optimal size / SVM-RFE curve minimum: the empirically CV-error-minimizing feature-subset size for the specific wrapped model.
- Stability selection frequency: fraction of resamples in which this CpG survived a fixed-lambda LASSO fit.
- Consensus weighted_score / n_methods: how many (weighted) of the run methods independently selected this CpG.
- Validation AUC: probability the classifier ranks a randomly chosen positive sample above a randomly chosen negative sample, on out-of-fold data.

### Biological meaning
A CpG appearing in the final panel is a candidate site whose methylation level is statistically/computationally associated with the phenotype contrast under this specific cohort, filtering choices, and method ensemble — a starting point for further biological investigation (e.g. checking the annotated gene, CpG-island context if available, prior literature).

### What it does NOT mean
- **Not causality.** Association (even multi-method-consensus association) does not establish that methylation at this CpG causes the phenotype, or vice versa (reverse causation and confounding are both still possible) — this module performs no Mendelian-randomization-style causal inference (that is the separate, dedicated Mendelian Randomization submodule elsewhere in this app).
- **Not biological validation.** No wet-lab or orthogonal-platform confirmation (e.g. pyrosequencing) is performed by this module; a computational panel is a hypothesis-generating output, not a validated biomarker.
- **Not clinical utility.** An AUC computed on this cohort's held-out folds does not establish clinical-grade sensitivity/specificity, calibration, or utility in a different population — external validation on an independent cohort would be required for that claim, which this module cannot itself provide.
- **Not mechanistic evidence.** Neither high importance nor high consensus implies any specific molecular mechanism connecting methylation at that CpG to the phenotype.

---

## Limitations

- No genuine continuous-phenotype support anywhere in the workflow — every non-multi-group method, including the three nominally "continuous" univariate options, ultimately operates on the same binary reference/comparison contrast (§8). Unchanged by the post-audit fixes; this remains a real scope limitation.
- **[FIXED]** ~~No true multi-group (>2-level) contrast is reachable, despite one univariate option's label implying it (§8).~~ ANOVA and Kruskal-Wallis now genuinely test every level of the chosen phenotype column when more than two are present among the eligible samples (verified with synthetic 4-group data) — see "Post-Audit Fixes Applied." Every other method (LASSO/RF/RFE/Stability/Validation) remains binary-only by design, since `glmnet`'s `family="binomial"` and `pROC`'s AUC both require exactly two classes.
- No published cross-reactive-probe or population-MAF reference list is bundled — both filters are inert unless the user supplies their own list, a deliberate and disclosed scope choice rather than an oversight (§7).
- Stability Selection implements a simplified, single-reference-lambda version of the eponymous method (§12). Unchanged — a documented design trade-off, not a bug.
- **[DISCLOSURE FIXED, leakage itself not eliminated]** Neither validation mode is fully leakage-free; only the Frozen mode's leakage was originally self-disclosed in the UI text. The Nested mode's UI copy now discloses its residual Stage-1 leakage explicitly too, but the leakage itself remains (would require re-plumbing Stage 1's filtering to run per outer fold — left as follow-up work).
- **[PARTIALLY FIXED]** Class-imbalance correction was only reachable for LASSO, not Random Forest, despite the latter's helper function supporting it (§10) — a "Class weighting" control was added for Random Forest's `classwt`.
- No automated test coverage for this submodule was found in the repository at the time of this reading (`tests/testthat/` contains no `featureselection`-specific test file for the methylomics module), unlike several other recently-touched submodules in this codebase. This code-fix pass verified its changes with standalone functional smoke tests against synthetic data rather than adding to the repository's test suite; adding a proper `testthat` file for this submodule remains outstanding.

---

## Thesis Implementation Paragraph

The Methylomics ML Feature Selection submodule is organised into nine tabs that together carry out a complete feature selection workflow. The first tab, Data and Filters, prepares the input data by applying quality control filters commonly used in methylation studies, including filters for missing values, low variance probes, SNP associated probes, sex chromosome probes, cross reactive probes, and minor allele frequency, along with missing value imputation. The next five tabs each apply one feature selection method to this filtered data. Univariate Selection uses statistical tests such as the t-test, moderated t-test, ANOVA, Wilcoxon, Kruskal-Wallis, and correlation based methods. LASSO uses regularized logistic regression to select CpGs. Tree-Based Selection uses a Random Forest model to rank CpGs by importance. RFE / Wrapper Selection uses recursive feature elimination. Stability Selection uses repeated resampling to find CpGs that are selected consistently. The Consensus / Overlap tab then combines the CpGs chosen by these methods into one agreed panel, based on how many methods selected each CpG, with an option to remove highly correlated CpGs from this panel. The Selected Features tab displays this final panel along with its gene and chromosome annotation and allows individual CpGs to be inspected. The last tab, Model and Export, tests how well the selected panel performs using cross validation and allows the complete analysis, including every filter, method, and result, to be saved as one reproducible file. The module accepts input data in the form of a CpG by sample methylation matrix, given as beta values or M-values, together with a phenotype file describing the sample groups. This input can either be the app's own preloaded whole blood cohort or a matrix and phenotype file uploaded by the user, and both are processed in the same way. The output of the module is a ranked list of candidate CpGs from each method, a final consensus panel of selected CpGs with biological annotation, a cross validated estimate of how well this panel classifies the two sample groups, and a downloadable file containing the complete analysis for reproducibility.

## Very Brief Thesis Statement

Genome-wide methylation data (up to 412,492 CpGs) was reduced to a small, interpretable candidate panel using five independent machine-learning and statistical feature-selection methods — univariate testing, LASSO/elastic-net regularization, Random Forest importance, recursive feature elimination, and resampling-based stability selection — combined via an overlap-based consensus with optional correlation-based redundancy pruning. The resulting panel was evaluated with cross-validated classification AUC under two validation strategies of differing methodological rigor, producing a prioritized, biologically annotatable set of candidate CpGs together with a fully reproducible record of every filtering, selection, and validation parameter used to derive it.

---

## Final Audit Summary

**What the module does:** Provides a fully interactive, nine-tab CpG feature-selection workflow for a two-group methylation phenotype contrast, from raw/preloaded matrix through five independent selection methods to a validated, exportable consensus panel.

**How many tabs it contains:** 9 — Data & Filters, Univariate Selection, LASSO, Tree-Based Selection, RFE / Wrapper Selection, Stability Selection, Consensus / Overlap, Selected Features, Model & Export.

**Main input:** A CpG × sample beta or M-value matrix (preloaded GSE42861 whole-blood cohort or user upload) plus an optional phenotype sheet defining a two-level group contrast.

**Main processing steps:** Sample/probe QC filtering → imputation → beta/M-value derivation → five parallel, independently-run selection methods → weighted consensus → optional correlation-based redundancy reduction → cross-validated evaluation → reproducible export.

**Main ML methods:** Univariate statistical testing (`limma` + base-R rank tests), LASSO/Elastic Net (`glmnet`), Random Forest (`randomForest`), three RFE flavors (`caret` + custom SVM-RFE via `e1071`), resampled stability selection (`glmnet`), plus `caret`-based cross-validated classification for panel evaluation.

**Main outputs:** Per-method ranked/selected CpG tables and plots, a consensus intersection table with Venn/UpSet/rank/heatmap visualizations, a final annotated CpG panel, cross-validated AUC estimates, and a comprehensive `.rds` reproducibility bundle.

**Normalisation/preprocessing approach:** No independent normalization (delegated to the upstream Methylomics Normalization submodule / the preloaded cohort's own QC pipeline); this module's own preprocessing is unsupervised filtering (missingness, variance/SD/IQR, optional annotation/list-based filters), a top-variance candidate cap, per-probe imputation, and a beta→M-value logit transform — all correctly sequenced and consistently applied across every downstream method.

**Strong aspects:** Genuinely multi-method, transparent, fully re-runnable design; consistently clean error handling via `validate()`/`tryCatch` almost everywhere; an unusually honest, self-disclosed acknowledgment of the Frozen validation mode's optimism directly in the UI copy; a comprehensive, genuinely-reproducible export bundle; correct, verified reuse (never modification) of shared QC/annotation/plotting infrastructure.

**Scientific concerns (as originally audited; see "Post-Audit Fixes Applied" for what changed):** A structurally unreachable "ANOVA (>2 groups)" option that silently computed a two-group statistic instead — **fixed**, ANOVA/Kruskal-Wallis now genuinely test all levels of the phenotype column when more than two are present; a simplified (single-lambda) implementation of stability selection relative to the term's stricter published meaning — unchanged, a documented design trade-off rather than a bug; the "Leakage-safe" validation mode's name overstating what it actually guarantees — disclosure corrected, residual leakage itself intentionally left as follow-up work.

**Software concerns (as originally audited):** Two class-imbalance-handling parameters (Random Forest's `classwt`/`sample_fraction`) present in helper code but unreachable from the UI — `classwt` **fixed** (a Class Weighting control was added), `sample_fraction` still unexposed; silent truncation/degradation in two places (RFE size list, k-NN imputation fallback) without user-facing notice — **both fixed**; incomplete package-version metadata in the reproducibility export — **fixed**. Also newly found and fixed in this pass: the "Preloaded whole-blood cohort" data-source option was gated on the wrong availability flag, which could have offered a preloaded option that could never actually load in a deployment shipping only the lightweight pipeline tables.

**Reproducibility concerns:** Every stochastic method is correctly seeded and the `.rds` export captures the full configuration; the package-version record now includes `limma`/`pROC` (previously missing); no automated test suite currently covers this submodule — this remains a genuine gap, and none was added as part of this code-fix pass (which relied on standalone functional smoke tests against synthetic data instead — see "Post-Audit Fixes Applied").

**Most important issue to verify:** Whether the "Leakage-safe" nested validation mode's residual Stage-1 leakage (unsupervised filtering/imputation computed globally, before the outer fold split) is acceptable for this project's intended use of the reported AUC — the UI now discloses this explicitly, but eliminating it would require re-plumbing Stage 1's filtering to run per outer fold, which was judged too invasive to attempt without a live test environment to verify against.

**Overall assessment:** A well-engineered, unusually transparent, and methodologically thoughtful feature-selection workflow — its explicit multi-method consensus design and self-disclosed validation caveats are genuine scientific-software strengths. Of the three concrete, code-verified High-severity gaps originally found (a mislabeled/unreachable method option, an unguarded crash edge case, and a partially-leaky "leakage-safe" validation path), the first two are now fixed and verified with functional smoke tests, and the third has its disclosure corrected with the residual leakage itself left as documented follow-up work requiring a live test environment to safely re-architect. Both the preloaded-cohort and user-upload data paths were confirmed to share the same downstream code with no source-specific special-casing, and the one preloaded-specific gating bug found (wrong availability flag) has also been fixed.
