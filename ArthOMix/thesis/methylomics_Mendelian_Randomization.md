# Methylomics — Mendelian Randomization Module: Code Audit and Thesis Documentation

**Source file:** `ArthOMix/R/methylomics/mod_methyl_mr.R` (1,210 lines) — UI + server for the "Mendelian Randomization" sub-module.

**Helper functions/objects this module calls, defined elsewhere (traced and read in full for this audit):**
- `ArthOMix/global.R` — `load_default_mr_estimates()` (:554), `load_default_mr_harmonised()` (:568), `load_default_mr_instrument_counts()` (:589), `read_uploaded_table()` (:1186), `guess_gwas_col()` (:1181), `GWAS_COL_PATTERNS` (:1170), `gwas_col_map_ui()` (:1229), `ARTHOMIX_COLORS` (:1417), `theme_arthomix()` (:1438).
- `ArthOMix/R/methylomics/annotation.R` — `methyl_get_annotation()` (:48).
- `ArthOMix/data_paths.R` — `METH_DATA_AVAILABLE` (:85), `METH_MR_DIR`.
- `ArthOMix/R/submodules_registry.R` — module registration (:48).
- `ArthOMix/server.R` — module invocation (:95).
- Not called by this module, despite existing as a sibling shared helper: `estimate_mr_set()` (`global.R:1269`) and `load_mr_instrument_table()` (`global.R:1389`), which belong to `R/transcriptomics/mod_mr.R` and `mod_crossancestry.R`'s eQTL-MR engine — see §6.1 for why this matters.
- Precomputed data read by the "Use Preloaded Data" route: `ArthOMix/data/preloaded/methylomics/tables/script08_mendelian_randomization/tables/{mr_estimates_female.csv, mr_estimates_male.csv, mr_harmonised_all_cpgs.csv, instrument_counts.csv, mr_steiger_female.csv, mr_steiger_male.csv}`, and the pipeline write-up `METHODS_mendelian_randomization.md` in the same folder, all inspected directly for this audit.

Prepared: 2026-08-26.

This document is derived **exclusively** from the code and data files cited above. Every non-trivial technical claim carries a `file:line` citation or a direct data check. Two label conventions, matching this project's other methylomics thesis documents (e.g. `methylomics_quality_control.md`), are used throughout:

- **Scientific background:** general Mendelian randomization / methylation-QTL knowledge (textbook/literature), not a claim about this code.
- **Code evidence:** a claim about what `mod_methyl_mr.R` (and the helpers it calls) actually does, always with a citation.

Per the audit brief, this document does **not** propose redesigning the MR workflow, and no code outside `mod_methyl_mr.R` and the helper functions it calls was modified to produce it. Findings that would require touching another module (e.g. `mod_methyl_coloc.R`, `R/transcriptomics/mod_mr.R`) are noted for context only, never acted on.

---

## 1. Module purpose

**Scientific background.** Observational analyses — differential methylation testing, region-level calling, co-methylation network membership, ensemble machine-learning feature selection — can only show that methylation at a CpG site is *associated* with disease status. They cannot distinguish whether methylation differences cause disease, are caused by disease (reverse causation, e.g. inflammation altering blood cell composition), or merely share an unmeasured confounder (e.g. smoking) with disease risk. Mendelian randomization (MR) addresses this by using germline genetic variants that influence methylation at a CpG — cis-acting methylation quantitative trait loci (cis-mQTLs) — as **instrumental variables**. Because genotype is fixed at conception, it cannot be caused by adult disease status and is far less susceptible to the confounding that affects a plain association test. If a genetic variant that raises methylation at a CpG also, and only through that CpG, raises disease risk, this is evidence *consistent with* a causal effect of methylation on disease — not proof of one, and the accepted academic phrasing throughout this document follows that convention.

**Code evidence — declared scope.** The module's own header comment states it is "a complete two-sample MR workflow for METHYLATION exposures: cis/trans mQTL instruments for a CpG -> a GWAS outcome. Not the transcriptomics eQTL/gene-expression MR tool (`R/transcriptomics/mod_mr.R`, untouched)" (`mod_methyl_mr.R:1-7`). The registration object states the same scope for the module menu: "Two-sample Mendelian randomisation of mQTL instruments against a GWAS outcome, per CpG or across a panel. Uses the bundled GoDMC/RA-GWAS data by default, or your own uploaded summary statistics." (`mod_methyl_mr.R:63-66`).

**Two data routes**, chosen first via a `radioButtons("data_source", ...)` on Tab 1 (`mod_methyl_mr.R:98-105`), and nothing downstream renders until one is chosen and validated:

1. **"Use Preloaded Data"** — reproduces the completed `script08_mendelian_randomization` run: already-selected, already-LD-clumped (r²<0.001, 10,000 kb), already `harmonise_data(action=2)`'d GoDMC cis-mQTL instruments for the script07 majority-vote CpG panel, harmonised against the Ishigaki et al. (2022) rheumatoid-arthritis (RA) GWAS (`mod_methyl_mr.R:9-24`). Instrument selection and LD clumping are **not** re-run live in this route — the raw GoDMC/RA-GWAS source files are not bundled with the deployment. MR estimation, sensitivity analyses, single-SNP estimates, and every plot **are** computed live from the cached, already-harmonised table (`mr_harmonised_all_cpgs.csv`), using the same `TwoSampleMR` functions the original pipeline script called.
2. **"Upload Dataset"** — a fully live pipeline on the user's own exposure(mQTL)/outcome(GWAS) summary-statistic files: instrument selection, LD clumping (`ieugwasr::ld_clump()` against the OpenGWAS API — no local PLINK reference is bundled), harmonisation, MR estimation, sensitivity analysis, and plots all run live with user-adjustable parameters (`mod_methyl_mr.R:25-30`).

**Code evidence — stage gating.** "Nothing computes or renders ahead of an explicit click, at every stage: Data -> Filters & Instruments -> LD Clumping -> Harmonisation -> MR Analysis -> Sensitivity -> Results -> Plots. Changing a stage's own defining inputs invalidates every stage after it" (`mod_methyl_mr.R:32-36`) — implemented as a `reactiveValues` has-run flag set (`stage_flags`, `mod_methyl_mr.R:305`) plus a generic `invalidate_from()` helper (`mod_methyl_mr.R:306-309`) that clears every downstream stage flag. This pattern is verified accurate throughout §5 below.

**Exposure.** DNA methylation (beta value) at a single CpG probe, instrumented by that CpG's cis-mQTL SNPs.

**Outcome.** Rheumatoid arthritis (RA) case/control status on the Preloaded route (Ishigaki et al. 2022 GWAS); any user-uploaded binary or continuous GWAS trait on the Upload route.

---

## 2. Tab count

> The Mendelian Randomization module contains **8 sub-tabs**, as defined by the module UI code.

Defined in one `tabsetPanel(id = ns("mr_tabs"), type = "tabs", ...)` at `mod_methyl_mr.R:76-86`:

| # | UI label (exact) | Underlying UI function | Underlying `outputId`(s) driving its body |
|---|---|---|---|
| 1 | `"1. Data"` | `mmr_data_ui()` (`:92-146`) | `data_validation_ui` |
| 2 | `"2. Filters & Instruments"` | `mmr_filters_ui()` (`:150-152`) | `filters_tab_body`, `instruments_summary_ui` |
| 3 | `"3. LD Clumping"` | `mmr_clump_ui()` (`:185-187`) | `clump_tab_body`, `clump_summary_ui` |
| 4 | `"4. Harmonisation"` | `mmr_harmonise_ui()` (`:209-211`) | `harmonise_tab_body`, `harmonise_summary_ui` |
| 5 | `"5. MR Analysis"` | `mmr_analysis_ui()` (`:233-235`) | `analysis_tab_body`, `mr_summary_ui` |
| 6 | `"6. Sensitivity"` | `mmr_sensitivity_ui()` (`:278-280`) | `sensitivity_tab_body`, `sensitivity_results_ui` |
| 7 | `"7. Results"` | `mmr_results_ui()` (`:284-286`) | `results_tab_body` |
| 8 | `"8. Plots"` | `mmr_plots_ui()` (`:290-292`) | `plots_tab_body`, plus a nested `tabsetPanel` of 6 plot types (`:1123-1137`: Scatter, Forest, Funnel, Leave-one-out, Single-SNP forest, Manhattan overview) |

No tab's UI label diverges from its underlying function name in a misleading way — the numeric prefixes ("1.", "2.", ...) directly communicate the intended sequential dependency, which §5 confirms is real (each tab's body is gated behind the previous stage's `stage_flags` entry).

---

## 3. End-to-End Mendelian Randomization Pipeline

**Conceptual MR pipeline, annotated against what this code actually implements:**

```text
User input (data source + parameters)
   |
Input validation                         <- IMPLEMENTED (Tab 1: build_data_state(), :378-399)
   v
Data loading                             <- IMPLEMENTED (2 routes: preloaded CSV reads, or fread() on upload)
   v
Exposure/outcome preparation             <- IMPLEMENTED (TwoSampleMR::format_data() for Upload; already-done for Preloaded)
   v
Instrument selection (p-value, cis/trans,
MAF, min-F, min/max count)               <- IMPLEMENTED, Tab 2 (build_instruments_state(), :479-565)
   v
LD clumping                              <- IMPLEMENTED for Upload (ieugwasr::ld_clump(), live);
                                             LOOKUP ONLY for Preloaded (reads instrument_counts.csv,
                                             does not re-run clumping) - Tab 3 (:634-676)
   v
Allele harmonisation                     <- IMPLEMENTED, Tab 4 (TwoSampleMR::harmonise_data() for Upload;
                                             reads pre-harmonised mr_keep flags for Preloaded, :707-753)
   v
MR estimation (tiered by instrument
count: 1 SNP -> Wald ratio, 2 -> IVW,
>=3 -> user-selected method set)         <- IMPLEMENTED, Tab 5 (TwoSampleMR::mr(), :797-849)
   v
Primary causal estimate                  <- PARTIALLY IMPLEMENTED - see Finding H-1 (§10): no
                                             explicit "primary method" flag is carried on the results
                                             table; the multiple-testing step instead picks whichever
                                             method has the smallest p-value per CpG
   v
Sensitivity analysis:
  - Heterogeneity (Cochran's Q)          <- IMPLEMENTED, >=3 instruments only (TwoSampleMR::mr_heterogeneity())
  - Pleiotropy (MR-Egger intercept)      <- IMPLEMENTED, >=3 instruments only (TwoSampleMR::mr_pleiotropy_test())
  - Leave-one-out                        <- IMPLEMENTED, >=3 instruments only (TwoSampleMR::mr_leaveoneout())
  - Single-SNP estimates                 <- IMPLEMENTED, always (TwoSampleMR::mr_singlesnp())
  - Steiger / directionality             <- IMPLEMENTED but only as a CpG-aggregate test
                                             (TwoSampleMR::directionality_test()) - see Finding M-3 (§10);
                                             the ALREADY-COMPUTED per-instrument Steiger flags shipped in
                                             the cached data are never surfaced
  - MR-PRESSO (outlier/global pleiotropy)<- COMPUTED but never displayed - see Finding M-1 (§10)
   v
Multiple-testing correction              <- IMPLEMENTED across CpGs (BH or Bonferroni), but on a
                                             min-p-selected row per CpG - see Finding H-1 (§10)
   v
Result filtering / annotation            <- IMPLEMENTED (gene/chr/pos joined from the 450K manifest,
                                             cpg_annotation(), :977-986)
   v
Tables / plots                           <- IMPLEMENTED (2 result tables, 5 QC/sensitivity tables,
                                             6 plot types, all downloadable)
   v
Final interpretation                     <- NOT IMPLEMENTED IN CODE (no automated causal-claim
                                             verdict is generated; interpretation is left to the reader,
                                             consistent with the guidance in §9 not to overstate MR)
```

Two steps a standard MR write-up would list are **conceptually relevant but not separately implemented as code steps** in this live module:

- **Colocalization** (LD/haplotype-pleiotropy check) — implemented for the Preloaded route's data, but as a *lookup* into `coloc_results.csv` inside the sibling **`mod_methyl_coloc.R`** module, not inside `mod_methyl_mr.R` at all. This document does not audit `mod_methyl_coloc.R` (out of scope per the audit brief) beyond noting the boundary in §11.
- **F-statistic-based instrument-strength reporting** — implemented, but as a threshold *filter* (Tab 2) plus a diagnostics *table* (Tab 6), not as a single aggregate "instrument strength" summary statistic per CpG the way some MR papers report a single mean F.

---

## 4. Relationship Between MR Tabs

**Code evidence — dependency mechanism.** Every tab body is a `renderUI()` gated by `req(stage_flags$<previous stage>)` (e.g. `mod_methyl_mr.R:602` `req(stage_flags$data)` gates Tab 2's body; `:679` `req(stage_flags$instruments)` gates Tab 3; `:756` `req(stage_flags$clump)` gates Tab 4; `:853` `req(stage_flags$harmonise)` gates Tab 5; `:903` `req(stage_flags$mr)` gates Tab 6's controls; `:992` and `:1116` `req(stage_flags$mr)` gate Tabs 7 and 8). A tab whose upstream stage has not run yet shows nothing (an empty/blank `uiOutput`), not an error message — a user who skips straight to Tab 7 before running Tab 5 simply sees a blank tab with no explanation of what to do next.

```text
1. Data  ---data_state()--->  2. Filters & Instruments  ---instruments_state()--->  3. LD Clumping
                                                                                            |
                                                                                     clump_state()
                                                                                            v
7. Results  <---mr_state()---  5. MR Analysis  <---harmonise_state()---  4. Harmonisation
     |                              |
     v                              v (mr_state()$dat, $cpgs)
8. Plots                      6. Sensitivity  ---sensitivity_state()---> (Tab 7's QC summary only)
```

- **Shared reactive state:** a single linear chain of `reactiveVal()` objects — `data_state()` (`:376`) -> `instruments_state()` (`:477`) -> `clump_state()` (`:634`) -> `harmonise_state()` (`:707`) -> `mr_state()` (`:787`) -> `sensitivity_state()` (`:869`). Each is written only by its own stage's action-button `observeEvent()` and read by every later stage and by that stage's own summary UI.
- **Independent tab-level calculations:** none — every stage strictly consumes the previous stage's `reactiveVal()`, confirmed by reading each `build_*_state()`/`observeEvent()` body; there is no tab that recomputes anything from raw `input$` values without going through the chain.
- **Sequential dependency, enforced explicitly:** confirmed via `invalidate_from()` calls attached to every stage-defining input (e.g. `:337-339` any Data-tab input invalidates from "instruments" onward; `:597-599` any Filters input invalidates from "clump" onward; `:676` any Clumping input invalidates from "harmonise" onward; `:753` the harmonisation-strategy radio invalidates "mr"; `:850` MR method/CI/PRESSO inputs invalidate "sensitivity"). Changing an upstream parameter after downstream stages have already run correctly hides (does not merely grey out) those downstream results until their action button is clicked again — verified by `stage_flags$<x> <- FALSE` inside `invalidate_from()` (`:306-309`) combined with every downstream `renderUI`'s `req(stage_flags$<x>)` gate.
- **A button must be clicked before downstream results appear:** yes, at every stage (`validate_btn`, `instruments_btn`, `clump_btn`, `harmonise_btn`, `mr_run_btn`, `sensitivity_run_btn`) — confirmed, no `reactive()` in the chain recomputes automatically on input change; all writes to the `reactiveVal()` chain happen inside `observeEvent(input$<stage>_btn, ...)` handlers.
- **What happens if an upstream step has no valid results:** the relevant `build_*_state()` function calls `req()`/`validate(need(...))` (e.g. `:394` `validate(need(!is.null(exp_raw) && !is.null(out_raw), ...))`; `:527` `validate(need(nrow(exp_fmt) > 0, ...))`; `:735-736` a harmonisation-specific message about SNP-ID/allele-column mismatches) which shows a Shiny `validate()` message in place of the tab body and returns from the `observeEvent()` without setting the stage flag — the chain simply does not advance, and no partial/incorrect downstream state is produced. This was traced through all six `observeEvent(input$..._btn, ...)` bodies (`:401-407`, `:567-596`, `:636-675`, `:709-752`, `:797-849`, `:871-900`) and confirmed consistent throughout.
- **Duplicated calculations:** one confirmed instance — the Preloaded route recomputes `F_stat_recomputed <- (beta.exposure/se.exposure)^2` (`:489`) even though an equivalent `F_stat` column is already present in the cached `mr_harmonised_all_cpgs.csv` (verified bit-for-bit: cached `F_stat = 402.120396339169` for SNP `rs106111`/CpG `cg25598086` equals `(-0.2122724/0.0105856)^2` recomputed independently). Not a correctness issue (same formula), but redundant I/O-vs-compute; see §10, Finding L-2.
- **Hidden dependency — this module vs. the rest of the app:** `mod_methyl_mr_server(id, methyl_dataset, methyl_results = NULL)` (`:298`) is called as `mod_methyl_mr_server("mx_mr", methyl_dataset, methyl_results)` from the shared Methylomics module loop (`server.R:95`), receiving the same `methyl_dataset`/`methyl_results` `reactiveValues` objects every other Methylomics sub-module receives. **Neither parameter is referenced anywhere else in the 1,210-line file** (confirmed by an exhaustive `grep -n "methyl_dataset\|methyl_results"` returning only the signature line at `:298`). Practically: the MR module's own "Data" tab is **completely decoupled** from whatever dataset a user loaded on the Methylomics **Dataset** tab (QC-filtered matrix, normalized matrix, uploaded matrix — none of it is visible here); and no MR result is written back into `methyl_results` for a sibling module to read reactively. This mirrors the same architectural pattern documented for `mod_methyl_qc.R` in `methylomics_quality_control.md` §1 — a project-wide convention, not unique to this module — and it is consistent with `mod_methyl_coloc.R` independently re-reading `load_default_mr_harmonised()` from disk (`global.R:443`, inside `mod_methyl_coloc.R`) rather than sharing any reactive with this module.
- **Outputs not actually connected to later analysis:** the MR-PRESSO results computed in Tab 5 (§10, Finding M-1) and the cached per-instrument Steiger flags (§10, Finding M-3) are the two confirmed cases of computed-or-available-but-never-connected data.

---

## 5. Reactive programming audit

**Reactive primitives used and their role** (all in `mod_methyl_mr.R` unless noted):

| Construct | Where | Role |
|---|---|---|
| `reactiveValues()` | `stage_flags` (`:305`), `show_flags` (`:1142`) | Has-run flags per pipeline stage / per plot type |
| `reactiveVal(NULL)` | `data_state`, `instruments_state`, `clump_state`, `harmonise_state`, `mr_state`, `sensitivity_state` | One-shot snapshot of each stage's computed result, replaced wholesale on each run |
| `reactive()` | `pre_panel_cpgs` (`:314`), `exp_df_r`/`out_df_r` (`:345-346`), `cpg_sub`/`cpg_res` (`:1150-1151`) | Lazily (re)computed derived values with automatic dependency tracking — file parsing and plot-input slicing |
| `eventReactive()` | none | Not used in this module |
| `observeEvent()` | 6 stage-advance handlers (`validate_btn`, `instruments_btn`, `clump_btn`, `harmonise_btn`, `mr_run_btn`, `sensitivity_run_btn`) + 6 plot `show_*` handlers (`:1143-1148`) + several parameter-invalidation observers | Explicit, click-triggered state transitions, and downstream invalidation on defining-input change |
| `observe()` | none | Not used |
| `req()` | throughout `build_*_state()` functions and every gated `renderUI()` | Silently halts a reactive/output when a prerequisite is missing (e.g. no file uploaded yet) |
| `validate()`/`need()` | `build_data_state()`, `build_instruments_state()` (implicitly via `req`), harmonisation (`:731,735,746`), plot builders (`:1153-1162`) | User-facing inline error messages replacing the tab body, instead of a stack trace |
| `isolate()` | none | Not used — every input read that should NOT immediately trigger recomputation is instead handled via the has-run-flag/invalidate pattern rather than `isolate()` |
| `outputOptions(..., suspendWhenHidden = FALSE)` | attached to every `renderUI` gating a downstream tab (e.g. `:472, 607, 629, 688, 702, 763, 782, 856, 864, 911, 936, ...`) | Forces Shiny to keep computing a tab's body even while a *different* tab is the one currently visible — necessary here because the has-run flags and `reactiveVal` chain must update regardless of which of the 8 tabs is on screen, or a user could click "Run MR" on Tab 5, switch to Tab 7 immediately, and see stale/blank output until Tab 5 is revisited |

**What triggers each calculation:** exclusively the six named action buttons (§4). No calculation reruns purely because an input changed — input changes only ever flip a `stage_flags$<x> <- FALSE` (via `invalidate_from()`), never recompute anything themselves. This is a deliberate, consistently-applied idiom shared with `mod_methyl_dmp.R`/`mod_methyl_normalization.R` per the module's own header comment (`:36`).

**Whether results are recalculated unnecessarily:** no evidence of this — each `reactiveVal()` write happens once per button click, and `reactive()` values (`exp_df_r`/`out_df_r`) are Shiny-cached automatically and invalidated only when their own `fileInput` changes.

**Whether stale results can remain visible:** no — because every downstream `renderUI` is gated on the corresponding `stage_flags` entry, and every stage-defining input change flips that flag to `FALSE` before the user can see a stale table; the tab reverts to its empty/blank state rather than showing an outdated result next to a changed input. One exception worth flagging: the **F-statistic diagnostics table** on the Sensitivity tab (`output$f_stat_table`, `:938-945`) recomputes its own `weak` column directly from the *live current value* of `input$f_min_f` (`:942`) rather than the F-stat threshold that was actually in effect when Instrument Selection (Tab 2) last ran — so if a user changes the "Minimum F-statistic" numeric input on Tab 2 *without* re-clicking "Select Instruments", the Sensitivity tab's weak/not-weak flag for each SNP updates to reflect the new, not-yet-applied threshold, while the instruments that were actually retained into MR still reflect the *old* threshold. This is a narrow, low-impact staleness case — see §10, Finding L-3.

---

## 6. Complete Function Inventory

### 6.1 MR-specific statistical functions (detailed)

For each: statistical concept, code call site, inputs/outputs, and an audit verdict.

#### `TwoSampleMR::mr()`
**Package:** `TwoSampleMR` (v0.7.8 confirmed installed in this deployment).
**Location:** `mod_methyl_mr.R:811`, inside the per-CpG `lapply()` in the `mr_run_btn` handler.
**Purpose:** Runs one or more MR estimators (whichever the caller lists in `method_list`) on a single exposure/outcome harmonised data.frame and returns one result row per method.
**Input:** `sub` — the harmonised rows for one CpG (columns `beta.exposure`, `se.exposure`, `beta.outcome`, `se.outcome`, allele columns, etc., the standard `TwoSampleMR` harmonised-data schema); `method_list` — a character vector of method function names selected by `tier_methods()` (`:789-795`).
**Input data type:** data.frame, one row per instrument SNP.
**What it does:** Dispatches to the named estimator function(s) internally (`TwoSampleMR:::mr_ivw`, `:::mr_egger_regression`, etc.) and aggregates results.
**Output:** data.frame with columns `id.exposure, id.outcome, outcome, exposure, method, nsnp, b, se, pval` — one row per method actually run.
**Role in MR:** This is the module's single MR-computation entry point for every CpG on every run (both routes) — the only estimator engine this module uses, distinct from the `MendelianRandomization`-package-based `estimate_mr_set()` helper used by the sibling transcriptomics MR module (`global.R:1269-1372`). **Audit observation:** using two different MR engines (`TwoSampleMR::mr()` here vs. `MendelianRandomization::mr_ivw()`/`mr_egger()`/`mr_median()` in the transcriptomics module) across the same application is not a bug — each engine is legitimate and each module's engine matches its own pipeline's original engine (this one matches `script08`'s own `TwoSampleMR`-based methods per `METHODS_mendelian_randomization.md`) — but it does mean the two MR modules in this app are not numerically directly comparable estimator-for-estimator (e.g. default IVW model: `TwoSampleMR::mr_ivw`'s default is a fixed-effects model, whereas `estimate_mr_set()` explicitly overrides `MendelianRandomization::mr_ivw(model = "random")` — `global.R:1299` — for a documented, verified reason). This is worth being explicit about in any thesis Methods section that describes both MR modules together.

#### `tier_methods()`
**Package:** local (`mod_methyl_mr.R:789-795`).
**Purpose:** Implements the project's own tiered estimator-selection rule.
**Input:** `n_snp` (integer instrument count for one CpG), `selected` (the user's checked methods from Tab 5).
**What it does:** `n_snp == 1` -> `"mr_wald_ratio"` only; `n_snp == 2` -> `"mr_ivw"` only; `n_snp >= 3` -> intersects the user's selection with the full estimable set (`mr_ivw, mr_egger_regression, mr_weighted_median, mr_weighted_mode, mr_simple_mode`, plus `mr_raps`/`mr_penalised_weighted_median` if available).
**Output:** character vector of method function names, passed straight to `TwoSampleMR::mr()`.
**Role in MR:** directly implements the documented `METHODS_mendelian_randomization.md` §2.FF.2 tiered-hierarchy table (1 SNP -> Wald ratio; 2 -> IVW only; ≥3 -> full set) — **verified to match the published pipeline's own rule exactly.**
**Audit observation:** correct and faithfully implemented. Its one limitation: it does not encode which of the ≥3-instrument tier's several methods is the *primary* one — see Finding H-1 (§10) for the downstream consequence.

#### `mr_wald_ratio` (via `TwoSampleMR::mr()`)
**Scientific background:** with exactly one instrument, the causal estimate is simply `beta_outcome / beta_exposure` (the "Wald ratio"), with an approximate SE via the delta method. It is the simplest, and statistically weakest, MR estimator — a single instrument provides no internal check on pleiotropy.
**Code evidence:** invoked automatically for any CpG with `n_snp == 1` (`tier_methods()`, `:790`). In the current bundled Preloaded dataset, 6 of the 8 CpGs carrying any instrument fall into exactly this tier (verified directly against `mr_harmonised_all_cpgs.csv`: 6 CpGs have exactly one instrument row, 2 have exactly two).
**Audit observation:** correctly triggered only at `n_snp == 1`, matching the documented pipeline rule.

#### `mr_ivw` (Inverse-Variance Weighted)
**Scientific background:** the standard primary MR estimator with ≥2 instruments — a precision-weighted (inverse-variance-weighted) meta-analysis of each instrument's own Wald ratio, equivalent to a weighted linear regression of outcome SNP effects on exposure SNP effects through the origin.
**Code evidence:** run for every CpG with `n_snp >= 2` (`:791-794`), and it is the *only* method run when `n_snp == 2` (`:791`). Uses `TwoSampleMR`'s own default model (fixed-effects) — **note this differs from the transcriptomics module's explicit random-effects override** (§6.1 above).
**Audit observation:** correctly the sole/first method at every instrument-count tier, matching accepted practice and the documented pipeline.

#### `mr_egger_regression` (MR-Egger)
**Scientific background:** a weighted regression of outcome effects on exposure effects that, unlike IVW, does not force the regression line through the origin — its intercept estimates the average pleiotropic effect across instruments (the "InSIDE" assumption), providing a test of directional horizontal pleiotropy at the cost of much lower statistical power and requiring several instruments to be even approximately estimable.
**Code evidence:** included in the ≥3-instrument tier's estimable set (`:792`) and offered as a checked-by-default method on Tab 5 (`mmr_analysis_controls()`, `:254,261`).
**Audit observation:** correctly gated to ≥3 instruments only (never attempted at 1 or 2, where it is not meaningfully estimable) — matches accepted guidance and the documented pipeline. In the bundled Preloaded dataset, **no CpG reaches the ≥3-instrument tier at all** (verified: max instrument count is 2), so MR-Egger, weighted median/mode, simple mode, Cochran's Q, the MR-Egger intercept test, and leave-one-out are all mechanically inert for the Preloaded route's default panel — this is documented candidly in the pipeline's own write-up (`METHODS_mendelian_randomization.md` §2.FF.3: "none reached the three-instrument threshold... so that sensitivity suite is not reported for any CpG in this screen") and the module's own Sensitivity tab explains this live (`:906,921,924,927`).

#### `mr_weighted_median`
**Scientific background:** a consistent causal-effect estimator even if up to 50% of the instruments' total weight comes from invalid (pleiotropic) instruments, provided the majority (by weight) are valid — more robust than IVW to a minority of bad instruments, at a cost of somewhat lower power.
**Code evidence:** in the ≥3-instrument estimable set (`:792`); selected by default (`:261`).
**Audit observation:** correctly gated; standard implementation via `TwoSampleMR::mr_weighted_median`, no local reimplementation.

#### `mr_weighted_mode` / `mr_simple_mode`
**Scientific background:** "zero modal pleiotropy assumption" estimators (Hartwig, Davey Smith & Bowden 2017) — consistent if the most common causal-effect estimate among the instruments arises from valid instruments, even if that is not a majority. Weighted mode additionally weights by instrument precision; simple mode does not.
**Code evidence:** both in the ≥3-instrument estimable set (`:792`), both selected by default (`:261`).
**Audit observation:** correctly gated; no local recomputation of the mode-finding kernel — delegated entirely to `TwoSampleMR`.

#### `mr_raps` (MR-RAPS, robust adjusted profile score) / `mr_penalised_weighted_median`
**Scientific background:** additional robust estimators — MR-RAPS explicitly models measurement error in the exposure effect and downweights outliers via a robust loss function; penalised weighted median further penalises instruments whose individual Wald ratio departs strongly from the weighted-median estimate.
**Code evidence:** offered **only if** `exists("mr_raps"/"mr_penalised_weighted_median", where = asNamespace("TwoSampleMR"), inherits = FALSE)` returns `TRUE` (`mmr_robust_method_choices()`, `:245-250`) — i.e. the checkbox for each is added to the UI only when the installed package version genuinely provides it, rather than being hardcoded and silently failing at run time on an older install. **Verified for this audit:** both functions exist in the installed `TwoSampleMR` 0.7.8 (`exists(..., where = asNamespace("TwoSampleMR"))` returns `TRUE` for both), so both checkboxes are live and functional in this deployment, not dead UI.
**Audit observation:** correctly implemented defensive gating — a genuinely good practice this module follows consistently (see also the explicit non-inclusion of contamination mixture and Radial MR, `:264-265`, disclosed in-UI as unavailable/architecturally incompatible rather than silently omitted).

#### `TwoSampleMR::mr_heterogeneity()` (Cochran's Q)
**Scientific background:** tests whether the individual instruments' Wald-ratio estimates are more dispersed than expected from their standard errors alone — excess heterogeneity is a signal of horizontal pleiotropy (instruments not all estimating the same causal effect) or genuine effect-modifier variation across instruments.
**Code evidence:** `mod_methyl_mr.R:880`, run only for CpGs with `n_snp >= 3` (`cpgs_ge3`, `:875`), across `method_list = c("mr_ivw", "mr_egger_regression")` — i.e. Q is computed separately for the IVW and the MR-Egger fit, standard practice.
**Output:** displayed verbatim via `DT::datatable()` in the "Heterogeneity (Cochran's Q)" box (`:920-922, 947-950`); no local recomputation or reformatting of the numbers.
**Audit observation:** correctly gated to ≥3 instruments; correctly explained in-UI when unavailable (`:921`, "No CpG in this run has >=3 instruments - heterogeneity is not computable").

#### `TwoSampleMR::mr_pleiotropy_test()` (MR-Egger intercept)
**Scientific background:** formal significance test of the MR-Egger intercept from the same regression fit as `mr_egger_regression` — a non-zero intercept indicates average directional (unbalanced) horizontal pleiotropy across the instrument set, one operationalisation of the MR "exclusion restriction" assumption (see §7).
**Code evidence:** `mod_methyl_mr.R:881`, same ≥3-instrument gate as Cochran's Q. Displayed verbatim (`:923-925, 952-955`).
**Audit observation:** correctly gated and correctly labelled "Horizontal pleiotropy (MR-Egger intercept)" — the UI label matches exactly what the function tests, no overstatement.

#### `TwoSampleMR::mr_leaveoneout()`
**Scientific background:** re-runs IVW after excluding each instrument in turn, one at a time — a single influential/outlying instrument that swings the pooled estimate is visible as a leave-one-out estimate far from the rest.
**Code evidence:** `:882`, same ≥3-instrument gate. Plot version: `build_loo()` (`:1156`), gated additionally at plot-build time by `validate(need(nrow(cpg_sub()) >= 3, "Needs >=3 instruments."))`.
**Audit observation:** correctly gated; table and plot both consistent with the same ≥3 threshold.

#### `TwoSampleMR::mr_singlesnp()`
**Scientific background:** reports each individual instrument's own Wald-ratio estimate alongside the pooled method estimate(s) — the standard input to both a forest plot and a funnel plot of MR results.
**Code evidence:** `:884` (Sensitivity tab, run on the *full* multi-CpG harmonised set `dat`, unconditionally — no instrument-count gate, since a single-SNP estimate is always computable even for a 1-instrument CpG); also called independently, per-CpG, inside `build_forest()` and `build_funnel()` (`:1154-1155`) on `cpg_sub()`.
**Audit observation:** correctly ungated (always available, matching the tab's own claim "Single-SNP estimates and F-statistics are always shown", `:906`). **However**, see Finding M-2 (§10): the "Forest" and "Single-SNP forest" plot tabs both call the identical `build_forest()` function and therefore render byte-identical plots under two different tab labels.

#### `TwoSampleMR::directionality_test()` (Steiger / directionality)
**Scientific background:** a CpG-level test of whether the exposure (methylation) or the outcome (RA) explains more of the variance at each instrument, using each side's approximate sample size — a check against reverse causation (the outcome, not the exposure, is upstream of the instrument-phenotype relationship).
**Code evidence:** `:889`, gated on `has_n` — both `samplesize.exposure` and `samplesize.outcome` columns present and not all-NA (`:887-888`). Rendered as "Steiger directionality" (`:931-933`, `:967-971`).
**Input data type:** requires per-SNP sample-size columns, present in the Preloaded route's cached data (verified: `samplesize.exposure`/`samplesize.outcome` columns exist and are populated in `mr_harmonised_all_cpgs.csv`) but only present in an Upload-route run if the user explicitly mapped an `n` column for both exposure and outcome (`gwas_col_map_ui(..., extra_fields = "n")`, `:348-349`).
**Audit observation — confirmed limitation (see Finding M-3, §10):** this computes an aggregate, CpG-level directionality statistic. It is **not** the same computation as the per-instrument `steiger_dir`/`steiger_pval` columns already present in the cached `mr_harmonised_all_cpgs.csv` and in the dedicated `mr_steiger_{sex}.csv` files (both read directly from disk and inspected for this audit — confirmed present, e.g. `steiger_dir = TRUE, steiger_pval = 6.15e-70` for SNP `rs106111`), which were computed once upstream as part of the original `script08` pipeline (documented in `METHODS_mendelian_randomization.md` §2.FF.2, "Steiger filtering"). **Neither of those two pre-existing columns/files is read anywhere in `mod_methyl_mr.R`** (confirmed by grep for `"steiger"` in the module — the only hits are the label string and the live `directionality_test()` call). The Sensitivity tab therefore shows a real, correctly-implemented, but *different* directionality statistic under a heading a reader could easily mistake for the pipeline's own already-computed per-instrument Steiger filter.

#### `MRPRESSO::mr_presso()`
**Scientific background:** a simulation-based test (Verbanck et al. 2018) for horizontal pleiotropy via outlier detection — flags individual instruments as statistical outliers relative to the rest, tests a "global" null of no horizontal pleiotropy across the instrument set, and optionally re-estimates the causal effect after excluding flagged outliers. Requires ≥4 instruments (one more than MR-Egger, since it must be able to drop a candidate outlier and still have ≥3 remaining).
**Code evidence:** `mod_methyl_mr.R:829-843`, run only when the user checks "Also run MR-PRESSO outlier test" (default **off**, `:267`) **and** `requireNamespace("MRPRESSO", quietly = TRUE)` succeeds (verified installed in this deployment) **and** the CpG has `nrow(sub) >= 4` instruments (`:832`).
**Output:** a named list keyed by CpG, each holding only `global_p` (the MR-PRESSO global test p-value) — stored in `mr_state()$presso` (`:846`).
**Audit observation — CONFIRMED issue, see Finding M-1 (§10):** `mr_state()$presso` is never read again anywhere in the file after being set (confirmed by exhaustive grep for `"presso"` — every remaining hit is the UI checkbox/tooltip or the computation itself). No table, value box, or download surfaces it. A user who opts into this specifically-labelled, specifically-costly ("adds runtime per CpG", per its own tooltip, `:268`) computation receives no visible result for it whatsoever.

### 6.2 Data loading, validation, and formatting functions

| Function | Package | Location | Purpose | Audit status |
|---|---|---|---|---|
| `load_default_mr_estimates(sex)` | local (`global.R:554`) | called `mod_methyl_mr.R:316-320` | Reads `mr_estimates_{sex}.csv` to build the CpG picker list for the Preloaded route | Implemented; see Finding M-4 (§10) re: the picker caption |
| `load_default_mr_harmonised()` | local (`global.R:568`) | called `:383` | Reads the single shared already-harmonised instrument table (`mr_harmonised_all_cpgs.csv`) — the actual computational source of truth for the Preloaded route | Implemented, verified against the on-disk CSV |
| `load_default_mr_instrument_counts()` | local (`global.R:589`) | called `:642` | Reads `instrument_counts.csv` (per-CpG post-clumping SNP counts from the original run) for the Preloaded route's Clumping-tab summary | Implemented; see Finding M-5 (§10) re: potential drift from live filtering |
| `METH_DATA_AVAILABLE` | local constant (`data_paths.R:85`) | gates `:111,326,380` etc. | `TRUE` iff `METH_DATA_ROOT` exists on disk — the Preloaded route degrades to an explicit "not available" note (`:117-118`) rather than crashing, if this deployment hasn't bundled the results folder | Implemented, consistent with the rest of the Methylomics module family |
| `read_uploaded_table(path)` | local (`global.R:1186`) | called `:345-346` | `data.table::fread()` wrapped in `tryCatch`, returns `NULL` on any parse failure instead of throwing | Implemented; downstream `validate(need(!is.null(...)))` (`:394`) turns a `NULL` into a clean user-facing message |
| `guess_gwas_col(cols, patterns)` / `GWAS_COL_PATTERNS` | local (`global.R:1181,1170`) | used inside `gwas_col_map_ui()` and directly at `:358-369` for the CpG/chr/pos/gene extra columns | Regex-based best-guess column auto-detection (e.g. `^snp$|^rsid$` for the SNP-ID column) | Implemented; always a suggestion the user can override via the rendered `selectInput`s, never silently trusted |
| `gwas_col_map_ui(ns, file_input, df_reactive, prefix, label, extra_fields)` | local (`global.R:1229`) | `output$exp_map_ui`/`output$out_map_ui` (`:348-349`) | Shared column-mapping UI (SNP/beta/SE/p/EA/OA/EAF, plus optional extra fields) reused across `mod_mr.R`, `mod_coloc.R`, and this module | Implemented; the CpG-ID and chr/pos/gene extra columns specific to methylation are layered on top locally (`:353-371`) rather than pushed into the shared generic helper, a reasonable scope boundary |
| `TwoSampleMR::format_data()` | `TwoSampleMR` | `:526,730` | Converts an arbitrary uploaded data.frame plus explicit column names into the package's standard exposure/outcome schema (adds `beta.exposure`, `se.exposure`, etc.) | Implemented for the Upload route only; not called for Preloaded (already-formatted cached data is used as-is) |
| `TwoSampleMR::harmonise_data()` | `TwoSampleMR` | `:734` | Aligns exposure and outcome effect alleles, flags/resolves palindromic SNPs, sets the `mr_keep` QC flag | Implemented for the Upload route with a user-selectable `action` (1/2/3, `:213-223`); for the Preloaded route the cached data's own pre-existing `mr_keep`/`palindromic`/`ambiguous`/`remove` flags are read directly instead of re-running harmonisation (`:711-720`) |
| `ieugwasr::ld_clump()` | `ieugwasr` | `:649-656` | Live LD-clumping against the OpenGWAS API reference panel, population-selectable | Implemented for Upload route only; Preloaded route explicitly does not re-run this (`:641-645`) — see §9 |
| `cpg_annotation(cpgs)` | local (`:977-986`) | called by `results_tab_body`, `build_primary_results_df`, `build_manhattan` | Looks up chr/pos/gene for a vector of CpG IDs from the 450K Illumina manifest via `methyl_get_annotation("450K")` | Implemented; degrades to an all-`NA` data.frame if the annotation package/data isn't available (`:983-985`), rather than erroring |
| `methyl_get_annotation(array_type)` | local (`annotation.R:48`) | called by `cpg_annotation()` | Reads the Bioconductor 450K manifest's `Locations`/`Manifest`/`Other`/SNP-overlap objects directly (bypassing `minfi::getAnnotation()` for a documented namespace-conflict reason, `annotation.R:38-47`) | Implemented; `"450K"` is hardcoded in this module's call (`:978`) — a CpG genuinely only present on the EPIC array (not 450K) would silently fail to annotate (returns `NA` gene/chr/pos for that row, via the same `[cpgs, ...]` indexing returning `NA` for unmatched rownames) rather than erroring; see Finding L-1 (§10) |

### 6.3 Local UI/plotting/output helper functions

| Function | Location | Purpose | Audit status |
|---|---|---|---|
| `.mmr_tip(text)` | `:43` | Renders a small hover-tooltip icon next to a control | Cosmetic; correct |
| `.mmr_badge(label, color)` | `:45-48` | Builds an inline HTML colour badge (defined but not invoked anywhere else in the file — confirmed by grep for `.mmr_badge(` returning only the definition) | **Dead code** — defined, never called; see Finding L-4 (§10) |
| `.mmr_stage_order` | `:50` | The 6-element stage-name vector `invalidate_from()` walks | Correct, matches the 6 `stage_flags` entries exactly |
| `invalidate_from(stage)` | `:306-309` | Clears every `stage_flags` entry from `stage` to the end of `.mmr_stage_order` | Correct; audited in §4/§5 |
| `mmr_robust_method_choices()` | `:245-250` | Namespace-existence-gated optional method list | Correct; audited in §6.1 |
| `build_data_state()` / `build_instruments_state()` | `:378-399` / `:479-565` | Pure(ish) builder functions called from inside the corresponding `observeEvent`, wrapped in `tryCatch` at the call site (`:402,568`) so a `validate()`/`stop()` inside them is swallowed rather than crashing the button handler | Correct pattern, consistently applied |
| `build_primary_results_df()` | `:1018-1026` | Joins `cpg_annotation()` onto `mr_state()$results` | Correct; used by both the on-screen table and the CSV download, so what's displayed and what's downloaded are always identical (verified: `output$dl_mr_results` at `:1071-1074` calls the same function) |
| `adjusted_results_df()` | `:1041-1048` | Multiple-testing correction | Correct mechanically (`stats::p.adjust`), but see Finding H-1 (§10) for the min-p selection issue upstream of it |
| `build_scatter/forest/funnel/loo/manhattan()` | `:1153-1169` | One plot-builder function per plot type, each wrapped in `validate(need(...))` for its own minimum-instrument-count precondition | Correct, with the one exception noted in Finding M-2 |
| `render_plot_slot(show_flag, build_fn, dl_name)` | `:1171-1179` | Generic "click Show plot to render" wrapper, avoiding rendering all 6 plot types on tab load | Correct; note the `build_fn` parameter is accepted but **unused inside the function body** — the actual plot is rendered by a separately-bound `output$plot_<name> <- renderPlot(build_fn())` call outside this helper (`:1187-1192`), so `render_plot_slot()`'s own `build_fn` argument does nothing (the function only uses `show_flag`/`dl_name`); this is a harmless but slightly misleading parameter, not a functional bug, since every call site (`:1180-1185`) still passes the correct builder as its second positional argument out of habit/readability rather than necessity |
| `make_plot_dl(build_fn, base_name)` | `:1194-1199` | Generic PNG/PDF/SVG plot-download handler factory via `ggplot2::ggsave()` | Correct; format taken live from `input$plot_format` at download time |

### 6.4 Generic R / Shiny / package functions used (not individually detailed — listed for completeness per the audit brief)

**Shiny core:** `moduleServer`, `NS`, `reactiveValues`, `reactiveVal`, `reactive`, `observeEvent`, `req`, `validate`, `need`, `renderUI`, `uiOutput`, `renderPlot`, `plotOutput`, `downloadHandler`, `downloadButton`, `outputOptions`, `withProgress`/`incProgress`, `fileInput`, `radioButtons`, `checkboxInput`, `checkboxGroupInput`, `numericInput`, `selectInput`, `textInput`, `actionButton`, `tabsetPanel`, `tabPanel`, `conditionalPanel`, `fluidRow`/`column`, `tags$*`, `icon`, `br`, `div`, `p`, `strong`.
**DT:** `DT::renderDataTable`, `DT::dataTableOutput`, `DT::datatable`.
**dplyr:** `dplyr::tibble` (only, to build the `ld_clump()` input, `:650`).
**ggplot2/ggrepel:** `ggplot(), aes(), geom_point(), geom_hline(), labs(), ggrepel::geom_text_repel(), ggplot2::ggsave()` — all confined to `build_manhattan()`/`make_plot_dl()`.
**data.table/utils:** `data.table::fread` (inside `read_uploaded_table`), `utils::capture.output`, `write.csv`, `writeLines`.
**Base R/stats:** `if/else`, `switch`, `tryCatch`, `lapply`/`do.call`/`rbind`/`split`/`which.min`, `stats::qnorm`, `stats::pnorm`, `stats::pt`, `stats::p.adjust`, `stats::complete.cases`, `round`, `signif`, `grepl`, `intersect`, `duplicated`, `table`, `sprintf`, `%||%` (project-wide null-coalescing helper). These are used exactly as their documentation describes; none are individually audited as they carry no MR-specific behaviour.

---

## 7. Mendelian Randomization Assumption Audit

### 7.1 Relevance (instrument strongly associated with the exposure)

**Code evidence:** addressed on two fronts. (1) A genome-wide-significance p-value filter on the instrument-exposure association, default `5e-8` (`MMR_DEFAULT_CIS_PVAL`, `:58`, applied at `:491,530`), user-editable. (2) An explicit F-statistic filter, default `F >= 10` (`MMR_DEFAULT_MIN_F`, `:59`, the conventional weak-instrument threshold — computed as `(beta.exposure/se.exposure)^2` at `:489,529` and again per-instrument in the Sensitivity tab's diagnostics table, `:893-894`), with a checkbox (`f_exclude_weak`, default **checked**, `:177`) controlling whether weak instruments are actually dropped from downstream MR or merely flagged. **Weak instruments are flagged, not silently removed, when the checkbox is off** — verified: `d$weak_instrument <- weak` is always set (`:496,533`), and it is only subtracted from `retained` when `isTRUE(input$f_exclude_weak)` (`:508,560`), so a user can inspect exactly which instruments were considered weak even if they chose to keep them.
**Audit verdict: relevance is addressed**, with a conventional, disclosed, and user-adjustable threshold.

### 7.2 Independence (instrument independent of confounders)

**Scientific background:** true independence from unmeasured confounders cannot be verified from summary statistics alone — it is an assumption, supported (not proven) by using a randomly-inherited genetic variant, and partially defended in practice by LD-clumping (removing statistically non-independent, correlated variants within one exposure's own instrument set) and by using an ancestry-matched reference population for that clumping.
**Code evidence:** LD clumping is implemented live for the Upload route (`ieugwasr::ld_clump()`, `:649-656`) with a user-selectable population (`clump_pop`, EUR/AFR/AMR/EAS/SAS, `:194`) and configurable r²/window/index-SNP-p thresholds (defaults `r²<0.001`, `10,000 kb`, `p=1`, `:198-200`, matching `MMR_DEFAULT_CLUMP_R2`/`MMR_DEFAULT_CLUMP_KB`). The Preloaded route does not re-run clumping but reflects the original pipeline's own choice of the same r²/kb thresholds (`METHODS_mendelian_randomization.md` §2.FF.2). **This clumping only removes LD-correlated SNPs among an individual CpG's own candidate instrument set** — it does **not**, and cannot, verify that a retained instrument is free of confounding via population stratification, assortative mating, or dynastic effects; the module makes no claim that it does (no UI text overstates LD clumping as "independence guaranteed").
**Audit verdict: independence is partially addressed** (LD-correlation removal, ancestry-matched clumping population selectable) but, correctly, **not claimed to be fully guaranteed** anywhere in the UI or code comments — consistent with accepted MR practice, where independence remains an assumption rather than something summary-statistic-level analysis can fully verify.

### 7.3 Exclusion restriction (instrument affects the outcome only through the exposure)

**Code evidence:** addressed, but **only for CpGs reaching the ≥3-instrument tier**, via the MR-Egger intercept test (`TwoSampleMR::mr_pleiotropy_test()`, §6.1) and, optionally, MR-PRESSO's global pleiotropy/outlier test (§6.1) — though the latter's result is never surfaced to the user (Finding M-1). For CpGs with 1 or 2 instruments — the majority of the bundled Preloaded panel (8/8 of the CpGs that carry any instrument at all fall below the ≥3 tier, per §6.1) — **no pleiotropy/exclusion-restriction check runs at all**, a limitation the Sensitivity tab discloses explicitly in its own header text (`:906`) rather than silently omitting.
**Code evidence — what this module does *not* do:** biological/functional filtering of instruments (e.g. excluding a SNP known to lie in a gene coding region unrelated to methylation), cis-restriction as a pleiotropy-defence (implemented, but as an *instrument-selection* filter — see §8 — not framed anywhere in the code as an exclusion-restriction defence specifically), and colocalization (implemented, but in the separate `mod_methyl_coloc.R` module, out of this document's scope — see §11).
**Audit verdict: exclusion restriction is partially addressed**, gated strictly to sufficiently-instrumented CpGs, and the module is honest in-UI about when the check is unavailable. It does not, and does not claim to, guarantee the assumption holds for sparse-instrument CpGs.

---

## 8. Instrument Selection Audit

| Parameter | Present? | Default | User-configurable? | Code location |
|---|---|---|---|---|
| p-value threshold | Yes | `5e-8` | Yes (`f_pval`) | `:157-158, 491, 530` |
| Genome-wide significance | Yes (same as above) | `5e-8` | Yes | — |
| Cis/trans restriction | Yes | `"cis"` | Yes (`f_region_mode`); **informational-only for the Preloaded route** (see below) | `:159-160, 497-502, 539-554` |
| Cis window | Yes | `±1,000 kb` (`MMR_DEFAULT_CIS_WINDOW_BP`) | Yes (`f_cis_window`, in kb) | `:57, 164-165, 485` |
| LD clumping (r²) | Yes | `0.001` | Yes, Upload route only (`clump_r2`) | `:60, 198, 651` |
| LD clumping window (kb) | Yes | `10,000` | Yes, Upload route only (`clump_kb`) | `:61, 199, 652` |
| MAF/EAF filter | Yes | `0` (i.e. **no filter by default**) | Yes (`f_maf`) | `:168, 494-495, 531-532` |
| Palindromic-SNP handling | Yes, via `TwoSampleMR::harmonise_data()`'s own logic (see §9) | strategy "2" (infer strand from EAF) | Yes, Upload route only (`harmonise_action`) | `:213-223, 734` |
| Missing-SNP handling | Rows with missing key fields simply fail `format_data()`/`harmonise_data()`'s own row-level checks | n/a | n/a | `:526, 734` |
| Duplicate-SNP handling | Not explicitly deduplicated by this module's own code before formatting; duplicate counts are surfaced diagnostically on the Upload route's Data Validation panel (`"Duplicate SNP rows (exposure)"`, `:461`) but not automatically removed | n/a | No explicit dedup step | `:440, 461` |
| Weak-instrument filtering | Yes | `F >= 10` | Yes (`f_min_f`, `f_exclude_weak`) | `:59, 175-177, 496, 508, 533, 560` |
| F-statistics | Yes, computed and shown | n/a | n/a (always computed) | `:489, 529, 893-894` |
| Exposure sample size | Optional column mapping, Upload route only (`exp_n`) | none | Yes | `:523, 349` |
| Outcome sample size | Optional column mapping, Upload route only (`out_n`) | none | Yes | `:729, 349` |
| Minimum instruments per CpG | Yes | `1` | Yes (`f_min_instruments`) | `:172, 588-591` |
| Maximum instruments per CpG | Yes (keeps strongest-by-p-value, excess flagged not deleted) | `NA` (no cap) | Yes (`f_max_instruments`) | `:173, 576-587` |

**If a parameter is hard-coded:** none of the above are hard-coded without a UI control — every threshold in this table has a corresponding editable input. The one genuinely fixed value is the annotation array type inside `cpg_annotation()` (`"450K"`, `:978`), which is not an MR-selection parameter but affects gene/chr/pos labelling only (Finding L-1).

**If a parameter is absent:** "No explicit filtering step was identified" applies to two things: (1) automatic SNP deduplication before formatting (diagnosed, not auto-fixed); (2) any biological/annotation-based instrument exclusion (e.g. excluding SNPs in the MHC region, or known pleiotropic hotspots) — not implemented anywhere in this module, and not claimed to be.

**Preloaded-route caveat, stated plainly in code and honoured throughout:** cis/trans restriction, LD clumping, and harmonisation strategy are all **fixed** for the Preloaded route to whatever the original `script08` run used (cis ±1 Mb, r²<0.001/10,000 kb, harmonisation action 2) — the corresponding Tab 2/3/4 controls remain visible and interactive, but for that route their practical effect is either informational-only (cis/trans toggle, `region_mode_note` disclosed at `:624`) or fully inert (LD-clumping numeric inputs have no effect on the Preloaded route's `n_after`, since clumping is a lookup not a live recompute, `:641-645`) or explicitly disabled with a note (harmonisation strategy radio remains selectable but is accompanied by an explicit "was harmonised with strategy 2... switch to Upload Dataset to harmonise your own data with a different strategy" notice when `editable = FALSE`, `:213,224-225,757`).

---

## 9. Exposure and Outcome Data Audit

### Exposure (both routes)
DNA methylation at a CpG probe, instrumented by cis-mQTL SNPs. Required fields after formatting (`TwoSampleMR` standard schema): `SNP`, `beta.exposure`, `se.exposure`, `pval.exposure`, `effect_allele.exposure`, `other_allele.exposure`, `exposure` (= CpG ID). Optional: `eaf.exposure`, `samplesize.exposure`, `chr.exposure`/`pos.exposure` (needed only for the cis/trans filter on the Upload route). Upload-route column mapping is user-driven via `gwas_col_map_ui()` (core 7 fields) plus a methylation-specific extra block (`:353-371`: CpG-ID column, optional SNP chr/pos, optional CpG chr/pos, optional gene). The CpG-ID column doubles as both `phenotype_col` and `id_col` when calling `TwoSampleMR::format_data()` (`:520`).

### Outcome
- **Preloaded route:** rheumatoid arthritis case/control status, Ishigaki et al. (2022) European-ancestry GWAS (22,350 cases / 74,823 controls), read once as part of the cached `mr_harmonised_all_cpgs.csv` (outcome label taken from `harm$outcome[1]`, `:388`). Binary-outcome checkbox defaults **checked** (`pre_binary_outcome`, `:116`), correctly reflecting that this outcome is in fact case/control.
- **Upload route:** any user-supplied GWAS, with a free-text label (`outcome_label`, defaulting `"Uploaded outcome"`, `:132`) and a binary-outcome checkbox defaulting **unchecked** (`up_binary_outcome`, `:133`) — the module makes no assumption about outcome trait type until the user states one.

### Required identifiers (per `TwoSampleMR::format_data()`, verified against the actual call sites)
Effect allele, other allele, beta, standard error, p-value, and SNP ID are always required (both exposure and outcome, `:390-392` gates `req()` on all of them before the Data tab will validate). Sample size and EAF are optional but recommended — EAF is specifically needed to resolve palindromic SNPs during harmonisation (disclosed in the UI tooltip, `:1254`, "improves palindromic SNP resolution"). Chromosome/position are optional and needed only if the cis/trans filter is to be applied on the Upload route.

### What is *not* required
Ancestry/population labels beyond the clumping-population selector; a genome-build declaration (the Data Validation panel attempts a best-effort build guess from position ranges, `:444-447`, explicitly labelled as not reliably distinguishable between GRCh37/GRCh38 from position alone — an honest disclosed limitation, not a silent assumption).

---

## 10. Scientific and Statistical Audit Findings

Findings are ranked by severity: **HIGH** (could change the scientific interpretation of a result), **MEDIUM** (affects reproducibility/robustness/interpretation), **LOW** (usability/minor implementation), **INFORMATIONAL** (not an error, important for understanding the implementation).

### HIGH

**H-1. Multiple-testing correction is applied to a min-p-selected row per CpG, not a pre-specified primary method.**
- **File/location:** `mod_methyl_mr.R:1041-1048` (`adjusted_results_df()`).
- **What the code does:** `primary_per_cpg <- do.call(rbind, lapply(split(r, r$exposure), function(x) x[which.min(x$pval), , drop = FALSE]))` — for each CpG, picks whichever *method's* result row has the smallest raw p-value (out of up to 7 methods run at the ≥3-instrument tier: IVW, MR-Egger, weighted median, weighted mode, simple mode, and up to 2 robust methods), then applies BH/Bonferroni correction (`stats::p.adjust`) across CpGs using only that minimum p-value per CpG.
- **Why it matters:** selecting the smallest p-value among several correlated estimators computed on the *same* instrument set, before correcting for multiplicity, is a well-known source of anti-conservative bias (the effective number of independent tests per CpG is understated) — the resulting "Adjusted p" / "Significant" verdict in the "Adjusted results" table (`:1000-1006`) can overstate significance for CpGs with several estimable methods relative to CpGs (like every CpG in the current bundled Preloaded panel) that only ever produce one method's result.
- **Expected behaviour:** correction across CpGs using a single, pre-specified primary method's p-value per CpG (e.g. IVW when ≥2 instruments, Wald ratio at 1 — mirroring the tiered hierarchy's own documented "primary" concept, and the sibling `estimate_mr_set()` helper's explicit `res_table$primary <- res_table$method == primary_method` flag, `global.R:1368`, which this module does not replicate).
- **Actual behaviour:** as described above.
- **Consequence:** for any CpG panel where multiple CpGs reach the ≥3-instrument tier, the multiple-testing-corrected significance calls in the Results tab may not be reliable. (Note: this cannot currently be observed in the bundled Preloaded dataset, since every CpG there has ≤2 instruments and therefore only ever produces one method's result row — but it will manifest as soon as a CpG panel or an Upload-route dataset produces a ≥3-instrument CpG.)
- **Recommended correction (not applied, per audit scope):** designate and carry a `primary` flag on each CpG's result rows (the tiered method already unambiguously determines which method is "the" pre-specified one at 1 and 2 instruments; at ≥3, IVW is the conventional default primary per this project's own pipeline write-up) and correct across that one row per CpG.

### MEDIUM

**M-1. MR-PRESSO results are computed but never displayed, tabulated, or downloaded.**
- **File/location:** `mod_methyl_mr.R:828-846` (computation); confirmed absent from `:852-1210` (every subsequent output).
- **Code does:** stores `presso` (a per-CpG list of `global_p` values) into `mr_state()$presso` when the user checks "Also run MR-PRESSO outlier test."
- **Consequence:** the computation (explicitly flagged in its own tooltip as adding runtime, `:268`) runs for no visible benefit; a user cannot see the one number (`global_p`) it produces anywhere in the app.

**M-2. "Forest" and "Single-SNP forest" plot tabs are implemented identically.**
- **File/location:** `mod_methyl_mr.R:1127-1128,1133-1134` (tab definitions); `:1154,1184,1191,1204` (both bound to `build_forest`).
- **Code does:** both tabs call the same `build_forest()` function, both plot outputs (`plot_forest`/`plot_snp_forest`) render `build_forest()`, both downloads (`dl_forest`/`dl_snp_forest`) save `build_forest()`'s output under different filenames.
- **Consequence:** two apparently distinct plot views render byte-identical figures; whatever visual distinction ("single-SNP" implying, e.g., no pooled-method rows, vs. the full forest including pooled estimates) a reader would reasonably infer from the two different tab labels does not exist in the actual output.

**M-3. Cached per-instrument Steiger flags are never surfaced; the live "Steiger directionality" output is a different, aggregate statistic.**
- **File/location:** cached data — `mr_harmonised_all_cpgs.csv` (`steiger_dir`/`steiger_pval` columns) and `mr_steiger_{sex}.csv`, both verified present on disk; module code — `mod_methyl_mr.R:889,931-933,967-971` (live `TwoSampleMR::directionality_test()` only).
- **Consequence:** a reader who sees "Steiger directionality" in the Sensitivity tab and is aware (from `METHODS_mendelian_randomization.md` §2.FF.2) that the original pipeline performed instrument-level Steiger filtering could reasonably but incorrectly assume they are looking at that same per-instrument result; they are seeing a freshly-computed, CpG-aggregate statistic instead, and the original per-instrument flags are not accessible anywhere in the module's UI.

**M-4. The Preloaded-route CpG picker's caption overstates what it lists.**
- **File/location:** `mod_methyl_mr.R:325-334` (`pre_panel_cpgs()`, `output$pre_cpg_ui`).
- **Code does:** builds the picker's choices from `load_default_mr_estimates(sex)$exposure` — verified against `mr_estimates_female.csv`/`mr_estimates_male.csv` on disk to contain **6** and **2** unique CpGs respectively — while the picker's own caption text reads "the script07 majority-vote panel (n_votes >= 2) for this stratum" (`:332`), and `METHODS_mendelian_randomization.md` §2.FF.2 documents the actual majority-vote panels as **12** (female) and **9** (male) CpGs.
- **Consequence:** the 13 majority-vote CpGs (9 with zero GoDMC candidate association at all, 4 more failing the cis-significance/clumping/F-stat filters — both counts independently confirmed in `METHODS_mendelian_randomization.md` §2.FF.3) are silently absent from the picker with no on-screen explanation that the list has been narrowed to instrument-bearing CpGs only, not the full majority-vote panel the caption names.

**M-5. Preloaded-route clumping "after" count can diverge from the instrument count actually carried into MR.**
- **File/location:** `mod_methyl_mr.R:641-645`.
- **Code does:** `n_after = sum(counts$n_snp[counts$cpg %in% unique(d$exposure)])`, i.e. re-sums the *original* per-CpG instrument counts recorded in `instrument_counts.csv` for whichever CpGs remain, rather than `nrow(d)` (the row count of instruments actually retained by the live Filters & Instruments tab).
- **Consequence:** if a user raises the p-value/F-stat/MAF threshold above the pipeline's own defaults on the Preloaded route in a way that removes *some but not all* of a CpG's original instruments, the "After clumping" value box will still show that CpG's full original instrument count, not the smaller number of instruments that actually reach Harmonisation/MR.

### LOW

**L-1. `cpg_annotation()` hardcodes the 450K array for gene/chr/pos lookup.**
`mod_methyl_mr.R:978`. A CpG unique to the EPIC array (not present on 450K) would resolve to `NA` gene/chr/pos in the Results table and be silently dropped from the Manhattan-overview plot (which requires non-`NA` position, `:1162`) rather than erroring — a narrow edge case given the bundled panel is itself 450K-array-derived.

**L-2. Preloaded-route F-statistic is recomputed rather than read from the cached column.**
`mod_methyl_mr.R:489`; verified numerically identical to the pre-existing `F_stat` column in `mr_harmonised_all_cpgs.csv`. No correctness impact, minor redundancy only.

**L-3. F-stat "weak" flag on the Sensitivity diagnostics table can go stale relative to the instruments actually selected.**
`mod_methyl_mr.R:942` reads `input$f_min_f` live at render time rather than the threshold value that was in effect when "Select Instruments" (Tab 2) was last clicked — if the user changes that numeric input afterward without re-running Tab 2, the "Weak (< threshold)" column can disagree with which instruments the current MR run actually excluded.

**L-4. `.mmr_badge()` is defined but never called.**
`mod_methyl_mr.R:45-48`; confirmed via grep — dead code with no functional impact.

### INFORMATIONAL

**I-1. This module accepts, but never reads or writes, `methyl_dataset`/`methyl_results`.**
`mod_methyl_mr.R:298`; confirmed by exhaustive grep. Not a bug — a consistent, project-wide Methylomics architecture convention (also documented for `mod_methyl_qc.R`) — but important for understanding data flow: this module's own "Data" tab is fully decoupled from the shared Dataset tab, and no other module reactively consumes this module's results.

**I-2. Two MR computational engines coexist in this codebase.** `TwoSampleMR::mr()` here vs. `MendelianRandomization`-package-based `estimate_mr_set()` in the transcriptomics/cross-ancestry modules (`global.R:1269-1372`) — each matches its own module's original pipeline, but the two are not drop-in numerically equivalent (e.g. differing IVW default model), worth stating explicitly in any combined thesis Methods text.

**I-3. Robust methods (MR-RAPS, penalised weighted median) are conditionally offered and confirmed live in this deployment**, not dead UI — `TwoSampleMR` 0.7.8 provides both (verified via `exists()` against the installed namespace for this audit).

**I-4. Upload-route LD clumping depends on the OpenGWAS public API,** which as of the currently installed `ieugwasr` requires authentication; a failed/unauthenticated clumping call is caught and degrades gracefully to an "api_error" status with unclumped instruments carried forward and disclosed to the user (`:658-660`), rather than crashing — correct defensive handling of a dependency genuinely outside this module's control.

**I-5. `docs/user_guide/28-methylomics-mr.Rmd`, the project's own end-user documentation chapter for this module, is an unfilled placeholder** (bracketed template text only, `docs/user_guide/28-methylomics-mr.Rmd:1-38`) — there is no independently-authored end-user description of this module's behaviour to cross-check the implementation against; this thesis document is, at the time of writing, the first substantive documentation of the module's actual behaviour.

**I-6. Colocalization, the MR pipeline's documented complement for sparse-instrument CpGs (`METHODS_mendelian_randomization.md` §2.FF.1/§2.FF.3), lives in the separate `mod_methyl_coloc.R` module,** not in `mod_methyl_mr.R` — out of scope for this audit, noted only so a reader does not expect to find it here.

---

## 11. Common MR implementation problems — systematic checklist (per audit brief §18)

| Problem | Present in this module? |
|---|---|
| Incorrect exposure/outcome orientation | Not observed — exposure (`type="exposure"`, CpG) and outcome (`type="outcome"`, GWAS) are consistently labelled through `format_data()`/cached-data column suffixes (`.exposure`/`.outcome`) throughout |
| Incorrect beta direction | Not observed — betas pass through `TwoSampleMR`'s own harmonisation unmodified; no local sign-flipping code exists in this module |
| Missing allele harmonisation | Not present — harmonisation is always run (live, Upload route) or already-done-and-read (Preloaded route) before MR |
| Incorrect palindromic-SNP handling | Not observed — delegated entirely to `TwoSampleMR::harmonise_data()`'s own `action` parameter, user-selectable on the Upload route (§9) |
| Weak instruments | Addressed via F-stat filter (§7.1); not eliminated by design (flag vs. drop is user-controlled) |
| Inadequate instrument count | Disclosed, not hidden — every ≥3-instrument-gated output states plainly when it is unavailable |
| LD-correlated instruments | Addressed via clumping (§7.2), live for Upload, lookup for Preloaded |
| Duplicated SNPs | Diagnosed (Upload-route validation panel) but not auto-removed — see §8 |
| Incorrect SE handling | Not observed in this module's own code (no local SE recomputation outside the documented, verified F-stat/CI formulas at `:820-823`) |
| Incorrect OR conversion | `results$or <- exp(results$b)` etc. (`:825`) — standard, correct log-odds-to-OR transform, applied only when `binary_outcome` is set |
| Incorrect CI calculation | `ci_low/high <- b -/+ z * se` with `z <- qnorm(1 - alpha/2)` (`:820-823`) — standard normal-approximation CI, correctly derived from the user-selected `mr_ci_level` |
| Incorrect p-value calculation | Delegated to `TwoSampleMR::mr()` internals — not recomputed locally |
| Inappropriate estimator selection | Addressed by `tier_methods()` (§6.1) — correctly matches the documented pipeline hierarchy |
| Inappropriate MR-Egger use with too few instruments | Not present — gated to ≥3 (§6.1) |
| Multiple-testing problems | **Present — see Finding H-1** |
| Selection bias | The Preloaded-route candidate panel itself is a machine-learning-selected (majority-vote) CpG set from an upstream module, not derived within this one; this module makes no additional selective-reporting choice beyond H-1 |
| Horizontal pleiotropy | Addressed for ≥3-instrument CpGs only (§7.3); undetectable for sparser CpGs, disclosed |
| Sample overlap | Not checked in code (exposure = GoDMC blood-mQTL cohort, outcome = Ishigaki et al. RA GWAS — independent-cohort by construction of the bundled data; not independently re-verified as part of this code audit since it is a property of the source data, not the code) |
| Population/ancestry mismatch | Partially addressed — clumping population is user-selectable (Upload route); the bundled Preloaded outcome/exposure pairing's ancestry match is a data-provenance property discussed in `METHODS_mendelian_randomization.md` §2.FF.1, not something this module's code verifies at run time |
| Inconsistent genome builds | Not verified programmatically for the Upload route — the Data Validation panel's "Detected genome build" line (`:444-447`) is explicitly disclosed as a rough guess, not a hard check; Preloaded route is fixed to GRCh37 throughout (stated, `:420`) |
| Incorrect SNP identifiers | Not independently re-validated by this module beyond requiring a `SNP` column to exist |
| API failures | Handled gracefully for `ieugwasr::ld_clump()` (Finding I-4) |
| Missing GWAS data | Handled via `METH_DATA_AVAILABLE` gating and per-field `req()`/`validate()` |
| Silently dropped variants | The one confirmed near-miss is Finding M-4 (CpG-level, not variant-level) — no confirmed instance of variant-level silent dropping without an accompanying counted/disclosed exclusion flag |
| Failure to report excluded variants | Not present — every exclusion category (`excluded_pval`, `excluded_region`, `excluded_maf`, `weak_instrument`, `excluded_maxcap`) is retained as a column and surfaced in the Instrument Selection summary (`:609-628`) |
| Results shown before successful analysis | Not present — the stage-flag/`req()` gating (§4/§5) consistently prevents this |
| Stale reactive results | Addressed generally (§5); one narrow exception, Finding L-3 |
| UI claims that do not match the statistical implementation | **Present — Findings M-3 and M-4** |

---

## 12. Data flow with concrete objects (per audit brief §20)

**Preloaded route:**
```text
input$data_source = "preloaded", input$pre_sex, input$pre_cpgs
    v
load_default_mr_harmonised()                 [global.R:568]
    v
data_state()  (list: mode, cpgs, harmonised, binary_outcome, outcome_label)
    v
instruments_state()  (adds excluded_pval/region/maf, weak_instrument, retained columns)
    v
clump_state()  (d = filtered rows; n_after via instrument_counts.csv lookup)
    v
harmonise_state()  (harmonised = cs$d[mr_keep], summary counts)
    v
mr_state()  (results = TwoSampleMR::mr() output per CpG; dat = harmonise_state()$harmonised)
    v
sensitivity_state()  (het, pleio, loo, single_snp, directionality, f_tab)
    v
output$primary_results_table / output$adjusted_results_table / output$*_table (Sensitivity)
    v
output$dl_mr_results / dl_harmonised / dl_sensitivity / dl_qc_report / dl_full_report
```

**Upload route (exposure side shown; outcome side mirrors it via `out_*` inputs):**
```text
input$exp_file
    v
exp_df_r()  [read_uploaded_table(), global.R:1186]
    v
output$exp_map_ui  [gwas_col_map_ui(), global.R:1229]  ->  input$exp_snp/exp_beta/exp_se/exp_pval/exp_ea/exp_oa/exp_eaf/exp_n
output$exp_extra_map_ui  [mod_methyl_mr.R:353-371]      ->  input$exp_cpg/exp_snp_chr/exp_snp_pos/exp_cpg_chr/exp_cpg_pos/exp_gene
    v
data_state()$exp_raw
    v
TwoSampleMR::format_data(...)  ->  exp_fmt   [instruments_state(), mod_methyl_mr.R:526]
    v
instruments_state()$d  (filtered, retained-flagged)
    v
ieugwasr::ld_clump()  ->  clump_state()$d    [mod_methyl_mr.R:649]
    v
TwoSampleMR::harmonise_data(cs$d, out_fmt)  ->  harmonise_state()$harmonised
    v
mr_state() -> sensitivity_state() -> Results/Plots tabs (same as Preloaded route from here)
```

---

## 13. Tab-by-tab academic table (per audit brief §21)

| Tab | Purpose | Input Data | Key Functions | Statistical Method | Output Data | Connection |
|---|---|---|---|---|---|---|
| 1. Data | Choose and validate the exposure/outcome data source | Preloaded cached CSVs, or user-uploaded exposure (mQTL) + outcome (GWAS) files | `build_data_state()`, `read_uploaded_table()`, `gwas_col_map_ui()` | None (data validation only) | `data_state()` | Feeds Tab 2; any change here invalidates every later tab |
| 2. Filters & Instruments | Select genome-wide-significant, sufficiently strong, optionally cis-restricted instruments | `data_state()` | `build_instruments_state()` | p-value threshold, F-statistic (instrument-strength) filter, MAF filter | `instruments_state()` | Feeds Tab 3; invalidates Tabs 3–8 on change |
| 3. LD Clumping | Remove LD-correlated instruments per CpG | `instruments_state()$d` (retained rows) | `ieugwasr::ld_clump()` (Upload) / `load_default_mr_instrument_counts()` lookup (Preloaded) | LD clumping (r², window, index-SNP p) | `clump_state()` | Feeds Tab 4; invalidates Tabs 4–8 on change |
| 4. Harmonisation | Align exposure/outcome effect alleles; resolve/exclude palindromic SNPs | `clump_state()$d` + outcome data | `TwoSampleMR::harmonise_data()` (Upload) / cached `mr_keep` flags (Preloaded) | Allele harmonisation (strategy 1/2/3) | `harmonise_state()` | Feeds Tab 5; invalidates Tabs 5–8 on change |
| 5. MR Analysis | Estimate the causal effect per CpG, tiered by instrument count | `harmonise_state()$harmonised` | `tier_methods()`, `TwoSampleMR::mr()`, optional `MRPRESSO::mr_presso()` | Wald ratio / IVW / MR-Egger / weighted median / weighted mode / simple mode / MR-RAPS / penalised weighted median | `mr_state()` | Feeds Tabs 6, 7, 8; invalidates Tabs 6–8 on change |
| 6. Sensitivity | Assess heterogeneity, pleiotropy, robustness, and directionality | `mr_state()$dat` | `mr_heterogeneity()`, `mr_pleiotropy_test()`, `mr_leaveoneout()`, `mr_singlesnp()`, `directionality_test()` | Cochran's Q, MR-Egger intercept, leave-one-out, single-SNP, Steiger/directionality (CpG-aggregate) | `sensitivity_state()` | Feeds Tab 7's QC summary only |
| 7. Results | Present, annotate, correct, and export MR estimates | `mr_state()`, `sensitivity_state()` (QC summary only) | `cpg_annotation()`, `build_primary_results_df()`, `adjusted_results_df()` | Benjamini–Hochberg / Bonferroni correction | Primary + adjusted results tables, QC summary, 5 downloads | Terminal display tab |
| 8. Plots | Visualise MR estimates and diagnostics per CpG or across the panel | `mr_state()`, `cpg_sub()`/`cpg_res()` | `mr_scatter_plot()`, `mr_forest_plot()`, `mr_funnel_plot()`, `mr_leaveoneout_plot()`, local `build_manhattan()` | Visual representations of the same estimates computed in Tabs 5–6 | 6 plot types, each downloadable (PNG/PDF/SVG) | Terminal display tab |

---

## 14. Function inventory table (per audit brief §22)

| Function | Package | Purpose | Input | Output | Used in | Scientific role | Audit status |
|---|---|---|---|---|---|---|---|
| `TwoSampleMR::mr()` | TwoSampleMR | Multi-method MR estimation | Harmonised data.frame, method list | Per-method result rows | Tab 5 | Core causal-effect estimator | Implemented, correct |
| `tier_methods()` | local | Instrument-count-tiered method selection | `n_snp`, `selected` | Method name vector | Tab 5 | Encodes the project's tiered estimator rule | Implemented, correct |
| `TwoSampleMR::mr_heterogeneity()` | TwoSampleMR | Cochran's Q | Harmonised data (≥3 SNP CpGs) | Q, df, p per method | Tab 6 | Heterogeneity/pleiotropy signal | Implemented, correctly gated |
| `TwoSampleMR::mr_pleiotropy_test()` | TwoSampleMR | MR-Egger intercept test | Harmonised data (≥3 SNP CpGs) | Intercept, SE, p | Tab 6 | Directional pleiotropy test | Implemented, correctly gated |
| `TwoSampleMR::mr_leaveoneout()` | TwoSampleMR | Leave-one-out re-estimation | Harmonised data (≥3 SNP CpGs) | Per-SNP-excluded estimate | Tab 6, Tab 8 | Influential-instrument check | Implemented, correctly gated |
| `TwoSampleMR::mr_singlesnp()` | TwoSampleMR | Single-instrument estimates | Harmonised data | Per-SNP + pooled estimates | Tab 6, Tab 8 | Underlies forest/funnel plots | Implemented; see Finding M-2 |
| `TwoSampleMR::directionality_test()` | TwoSampleMR | Aggregate Steiger test | Harmonised data w/ sample sizes | Steiger p, correct-direction flag | Tab 6 | Reverse-causation check | Implemented; see Finding M-3 |
| `MRPRESSO::mr_presso()` | MRPRESSO | Outlier/global pleiotropy test | Harmonised data (≥4 SNP), optional | Global p, outlier list, corrected estimate | Tab 5 (computed only) | Simulation-based pleiotropy test | Computed but not surfaced — Finding M-1 |
| `TwoSampleMR::format_data()` | TwoSampleMR | Standardise uploaded summary stats | Raw data.frame + column names | Standard exposure/outcome schema | Tabs 2, 4 (Upload route) | Data preparation | Implemented |
| `TwoSampleMR::harmonise_data()` | TwoSampleMR | Allele harmonisation | Exposure + outcome formatted data | Harmonised data.frame with `mr_keep` | Tab 4 (Upload route) | Ensures effect alleles align | Implemented |
| `ieugwasr::ld_clump()` | ieugwasr | LD clumping via OpenGWAS API | rsid/pval/id tibble | Clumped SNP list | Tab 3 (Upload route) | Instrument independence | Implemented; API-dependent (Finding I-4) |
| `mr_scatter_plot()` | TwoSampleMR | Exposure-vs-outcome scatter with method fit lines | MR results + harmonised data | ggplot object | Tab 8 (Scatter) | Visual causal-effect summary | Implemented |
| `mr_forest_plot()` | TwoSampleMR | Forest plot of single-SNP + pooled estimates | `mr_singlesnp()` output | ggplot object | Tab 8 (Forest, Single-SNP forest) | Visual instrument-level summary | Implemented; duplicated across 2 tabs (Finding M-2) |
| `mr_funnel_plot()` | TwoSampleMR | Funnel plot (precision vs. estimate) | `mr_singlesnp()` output | ggplot object | Tab 8 (Funnel) | Visual asymmetry/pleiotropy check | Implemented |
| `mr_leaveoneout_plot()` | TwoSampleMR | Leave-one-out forest plot | `mr_leaveoneout()` output | ggplot object | Tab 8 (Leave-one-out) | Visual influential-instrument check | Implemented |
| `build_manhattan()` | local | Cross-CpG -log10(p) overview | `mr_state()$results`, `cpg_annotation()` | ggplot object | Tab 8 (Manhattan overview) | Panel-wide result overview | Implemented |
| `cpg_annotation()` | local | CpG -> gene/chr/pos lookup | CpG ID vector | data.frame keyed by CpG | Tabs 7, 8 | Biological labelling | Implemented; 450K-only (Finding L-1) |
| `adjusted_results_df()` | local | Multiple-testing correction | `mr_state()$results` | Corrected per-CpG table | Tab 7 | FDR/Bonferroni control | Implemented; min-p issue (Finding H-1) |
| `stats::p.adjust()` | stats | BH/Bonferroni correction | p-value vector | Adjusted p-value vector | Tab 7 | Multiple-testing control | Implemented correctly (mechanically) |
| `req()`, `validate()`, `need()` | shiny | Guard/error-display | Reactive prerequisites | Halted reactive / inline message | Throughout | Reactive-flow safety | Implemented, consistently |
| `DT::datatable()` | DT | Interactive result tables | data.frame | Rendered table widget | Tabs 6, 7 | Result presentation | Implemented |
| `downloadHandler()` | shiny | CSV/TXT/plot export | Reactive result object | File on disk, streamed to browser | Tabs 7, 8 | Reproducible export | Implemented; downloads match displayed data (verified for MR results CSV) |

(Base-language constructs — `if`/`else`, `{}`, vector indexing, `%||%` — are explained narratively above per the audit brief's own instruction not to tabulate them individually.)

---

## 15. Code-to-Thesis Mapping (per audit brief §23)

| Code implementation | Computational purpose | Statistical meaning | Scientific interpretation | Thesis wording |
|---|---|---|---|---|
| `tier_methods()` (`:789-795`) | Choose which estimator(s) to run per CpG | Matches estimator choice to statistical power available from the instrument count | Avoids reporting an MR-Egger/heterogeneity result that would not be meaningfully estimable | "MR estimators were applied per CpG according to a pre-specified hierarchy tiered by clumped instrument count (1 SNP: Wald ratio; 2: IVW; ≥3: the full multi-method set), following [the pipeline's own documented rule]." |
| `TwoSampleMR::mr(method_list = ...)` (`:811`) | Compute the causal point estimate(s) | Weighted regression / ratio estimator(s) of the exposure→outcome effect | Evidence consistent with (not proof of) a causal effect of CpG methylation on the outcome | "Two-sample MR was performed using the TwoSampleMR R package (Hemani et al., 2018)." |
| `f_min_f`/F-stat filter (`:489,529,894`) | Exclude/flag instruments with weak exposure association | Guards against weak-instrument bias | Ensures the instrument genuinely proxies the exposure before being trusted for causal inference | "Instruments were retained only where F ≥ 10, the conventional threshold for adequate instrument strength (Burgess et al., 2013)." |
| `ieugwasr::ld_clump()` / cached clumping (`:641-656`) | Remove statistically non-independent instruments | Enforces approximate instrument independence within the exposure's own candidate set | Supports (does not prove) the MR independence assumption | "Candidate instruments were LD-clumped (r² < 0.001, 10,000 kb) using the ieugwasr package (Hemani et al., 2018)." |
| `TwoSampleMR::harmonise_data()` (`:734`) | Align exposure/outcome effect alleles; resolve palindromic SNPs | Prevents a spurious sign-flip in the causal estimate from a mismatched allele coding | A prerequisite for any valid causal-effect estimate from summary statistics | "Exposure and outcome summary statistics were harmonised using the TwoSampleMR package's harmonisation procedure." |
| `mr_pleiotropy_test()`/`mr_heterogeneity()` (`:880-881`) | Test the exclusion-restriction assumption, where estimable | Non-zero Egger intercept / excess Q signals possible pleiotropy | A caveat that qualifies, rather than validates, an otherwise significant causal estimate | "For CpGs with ≥3 instruments, the MR-Egger intercept test and Cochran's Q statistic were used to assess directional and overall pleiotropy respectively." |
| `adjusted_results_df()` (`:1041-1048`) | Control the panel-wide false-discovery rate across CpGs | Benjamini–Hochberg / Bonferroni correction | Distinguishes a genuinely panel-wide-significant CpG from one significant only by chance across many CpGs tested | *Caveat required per Finding H-1*: correction is currently applied to each CpG's minimum p-value across its estimable methods, not a single pre-specified method — this should be corrected before being described in a thesis as a "primary-method-corrected" result, or the limitation stated explicitly if reported as-is. |

---

## 16. Reproducibility Information (per audit brief §24)

- **R packages actually invoked by this module (confirmed by call sites above):** `TwoSampleMR` (v0.7.8, confirmed installed for this audit), `ieugwasr`, `MRPRESSO`, `DT`, `ggplot2`, `ggrepel`, `data.table`, `shiny`/`shinydashboard`-family UI functions, `stats`, `utils`.
- **Not installed-version-pinned anywhere in this module's own code** — no `packageVersion()` check or `renv`/lockfile reference inside `mod_methyl_mr.R` itself; version reproducibility for a live Upload-route run therefore depends on whatever package versions are installed in the deployment environment at run time. (This is a deployment-level, not a module-level, reproducibility property — out of this file's scope to change.)
- **Random seeds:** the only stochastic step reachable from this module is `MRPRESSO::mr_presso()`'s outlier-test simulation, called with an explicit `seed = 2024` (`:837`) and `NbDistribution = 1000` — so MR-PRESSO's own result (when the checkbox is enabled) is deterministic across runs, even though it is never displayed (Finding M-1). No other step in the module is stochastic.
- **Hard-coded parameters:** `MMR_DEFAULT_CIS_WINDOW_BP = 1e6`, `MMR_DEFAULT_CIS_PVAL = 5e-8`, `MMR_DEFAULT_MIN_F = 10`, `MMR_DEFAULT_CLUMP_R2 = 0.001`, `MMR_DEFAULT_CLUMP_KB = 10000` (`:57-61`) — all are UI-editable *defaults*, not immovable constants; `MRPRESSO::mr_presso()`'s `NbDistribution = 1000`, `SignifThreshold = 0.05` (`:837`) and `seed = 2024` are not exposed to the UI at all.
- **Configurable parameters:** every value in the §8 instrument-selection table, plus `mr_ci_level` (CI width, default 0.95, `:270`), `padj_method` (BH/Bonferroni, default BH, `:1004`).
- **External APIs:** OpenGWAS (via `ieugwasr::ld_clump()`), Upload route only — a network dependency outside this module's or this deployment's control (Finding I-4).
- **Local data files (Preloaded route):** enumerated in the header of this document; all read-only, table-only lookups per the project-wide `METH_DATA_AVAILABLE` graceful-degradation contract.
- **Not available in the code and not invented here:** an explicit genome build declaration for Upload-route data (best-effort guess only, `:444-447`); a package-version manifest specific to this module.

---

## 17. Biological and Statistical Concepts, Plain-Language (per audit brief §19)

- **CpG:** a cytosine nucleotide immediately followed by a guanine in the DNA sequence — the site where DNA methylation (addition of a methyl group to the cytosine) is measured on the arrays this project uses. Each CpG is identified by a probe ID such as `cg25598086`.
- **DNA methylation / beta value:** the fraction of DNA molecules at a given CpG that carry a methyl mark in a sample, ranging 0 (fully unmethylated) to 1 (fully methylated); a quantitative, array-measurable trait like height or blood pressure, which is why it can itself be the "exposure" in an MR analysis.
- **Methylation quantitative trait locus (mQTL):** a genomic location where common genetic variation is statistically associated with methylation level at a nearby (cis) or distant (trans) CpG. A *cis*-mQTL lies close to (here, within ±1 Mb of) the CpG it influences.
- **Genetic instrument:** a genetic variant (typically a single-nucleotide polymorphism, SNP) used as a proxy for the exposure in an MR analysis, because it is (i) robustly associated with the exposure, (ii) not associated with known confounders, and (iii) affects the outcome only through the exposure. This module operationalises (i) via the p-value/F-stat filters (§7.1/§8), addresses (ii) partially via LD clumping (§7.2), and addresses (iii) only where instrument count permits (§7.3).
- **Summary statistics:** the aggregated result of a genome-wide association study (per-SNP effect size, standard error, p-value, and allele information) rather than individual-level genotype/phenotype data — this is what both the exposure (mQTL) and outcome (GWAS) files in this module consist of.
- **Why genetic instruments reduce confounding:** genotype is fixed at conception (Mendel's laws of independent assortment/segregation, hence "Mendelian" randomization), so it cannot be caused by, or systematically confounded with, adult lifestyle/environmental factors the way an observed methylation-disease association can be.
- **Why MR differs from an ordinary association test:** an ordinary DMP/DMR association test asks "is methylation at this CpG correlated with disease status?" — vulnerable to reverse causation and confounding. MR asks "does a genetic variant that only changes methylation at this CpG also change disease risk?" — a design intended to be robust to both.
- **Why pleiotropy is a concern:** if an instrument SNP affects the outcome through a pathway other than the exposure (e.g. it also happens to lie near a different, disease-relevant gene), the exclusion-restriction assumption is violated and the MR estimate is biased; MR-Egger and Cochran's Q (§7.3) are two of the checks this module implements specifically for this concern, when instrument count allows.
- **Why instrument strength (F-statistic) matters:** a weak instrument (low F) contributes little true signal relative to noise, biasing two-sample MR estimates toward the null and inflating variance — hence the F ≥ 10 default filter (§7.1, §8).
- **Why sensitivity analyses matter:** any single MR estimator makes its own specific assumption set (e.g. IVW assumes no pleiotropy at all; MR-Egger relaxes that but assumes pleiotropy is uncorrelated with instrument strength; weighted median/mode assume a majority/plurality of valid instruments) — agreement across several estimators strengthens confidence in a result; disagreement is itself informative, flagging that the estimate is sensitive to which assumption is chosen. This module implements this convergence check only where instrument count allows (§7.3), and is explicit in-UI about when it cannot.
- **Appropriate causal language:** consistent with the guidance in the audit brief, this document (and the module's own UI text, e.g. `:997`, "A single nominal p < 0.05 is not, on its own, evidence of a causal effect") uses phrasing such as "evidence consistent with a potential causal effect," never "proves causality."

---

## Thesis Implementation Paragraph

The Methylomics Mendelian Randomization sub-module (`mod_methyl_mr.R`) implements a two-sample MR workflow that tests whether DNA methylation at a candidate CpG is causally upstream of rheumatoid arthritis risk, using cis-acting methylation quantitative trait loci as genetic instruments. It is organised into **8 sequential sub-tabs** (Data, Filters & Instruments, LD Clumping, Harmonisation, MR Analysis, Sensitivity, Results, Plots), each gated behind an explicit user action so that no result is ever computed or displayed ahead of its prerequisite stage. The module offers two data routes: a "Use Preloaded Data" route that reproduces a completed pipeline run (GoDMC cis-mQTL instruments harmonised against the Ishigaki et al. 2022 rheumatoid arthritis GWAS) by recomputing MR estimation, sensitivity analyses, and plots live from the already-clumped, already-harmonised cached instrument table, and an "Upload Dataset" route that runs instrument selection, LD clumping (via the OpenGWAS API), harmonisation, and MR estimation entirely live on user-supplied summary statistics. MR estimation itself uses the `TwoSampleMR` R package's estimator suite (Wald ratio, inverse-variance-weighted, MR-Egger, weighted median, weighted mode, simple mode, and, where the installed package version provides them, MR-RAPS and penalised weighted median), applied under an instrument-count-tiered hierarchy consistent with the sparse-instrument reality of CpG-level mQTL data. Sensitivity analyses (Cochran's Q heterogeneity, MR-Egger intercept pleiotropy, leave-one-out, single-SNP estimation, and an aggregate directionality/Steiger test) are computed where the instrument count permits and explicitly disclosed as unavailable where it does not. The module is scientifically appropriate and largely faithful to its own documented pipeline, with a small number of confirmed implementation issues identified by this audit — most notably that its multiple-testing correction selects each CpG's minimum p-value across estimable methods rather than a single pre-specified primary method (§10, Finding H-1), and that an optional MR-PRESSO computation and two informational UI elements do not connect to any visible output (§10, Findings M-1–M-5) — none of which alters the module's basic MR machinery, but all of which should be weighed before quoting its multiple-testing-adjusted significance calls directly in a thesis results section.

## Thesis Paragraph — Tab-Based Implementation

Across its 8 tabs, the module proceeds as follows. The **Data** tab lets the user choose between the bundled, script08-reproducing preloaded dataset (with a CpG-panel picker restricted to the majority-vote CpGs that actually carried a usable genetic instrument in the original run) and an uploaded exposure(mQTL)/outcome(GWAS) file pair, validating both column mapping and basic data quality before proceeding. **Filters & Instruments** applies a genome-wide-significance p-value threshold, an optional cis/trans restriction, a minor-allele-frequency filter, and an instrument-strength (F-statistic ≥ 10 by default) filter, retaining every exclusion reason as an inspectable flag rather than silently discarding rows. **LD Clumping** removes statistically non-independent instruments — live, via the `ieugwasr` package's OpenGWAS-backed clumping, for uploaded data, or by reference to the original pipeline's own recorded clumping counts for the preloaded data. **Harmonisation** aligns exposure and outcome effect alleles and resolves or excludes palindromic variants, again either live (uploaded data, with a user-selectable harmonisation strategy) or by reading the pipeline's already-computed harmonisation flags (preloaded data). **MR Analysis** estimates the causal effect per CpG under a tiered estimator hierarchy — a single instrument yields a Wald ratio, two instruments yield inverse-variance-weighted estimation alone, and three or more instruments additionally support MR-Egger, weighted median, weighted mode, simple mode, and, where available, MR-RAPS and penalised weighted median. **Sensitivity** computes heterogeneity (Cochran's Q), pleiotropy (MR-Egger intercept), leave-one-out re-estimation, single-SNP estimates, and an aggregate directionality/Steiger test, each explicitly disclosed as unavailable for any CpG below its required instrument count. **Results** presents a per-CpG-per-method estimate table, a Benjamini–Hochberg/Bonferroni-corrected panel-wide table, a QC summary, and five downloadable exports; **Plots** renders scatter, forest, funnel, leave-one-out, single-SNP-forest, and cross-CpG Manhattan-style overview figures, each downloadable in PNG, PDF, or SVG. Each tab strictly consumes the previous tab's stored result via a chain of reactive values, and changing any stage-defining input correctly invalidates every downstream stage until its action button is clicked again — a consistently and correctly implemented sequential architecture, with the specific caveats on multiple-testing correction, one duplicated plot tab, one silently narrowed CpG picker, one dormant computation (MR-PRESSO), and one dormant cached column (per-instrument Steiger flags) documented in full in §10 of this audit.

---

## Final Audit Summary

### What is correctly implemented
- The 8-tab sequential pipeline architecture and its reactive-invalidation contract (§4, §5) — verified consistent throughout.
- The instrument-count-tiered MR estimator hierarchy, matching the documented pipeline exactly (§6.1, §8).
- F-statistic-based weak-instrument handling, with a genuine flag-vs-drop distinction (§7.1).
- LD clumping (live for uploads; faithful lookup for preloaded data) (§7.2, §9).
- Allele harmonisation, live or read-through, with disclosed QC counts (§4, §9).
- Sensitivity analyses (Q, Egger intercept, leave-one-out, single-SNP), correctly gated to instrument count and honestly disclosed when unavailable (§6.1, §7.3).
- Defensive, namespace-existence-gated optional methods (MR-RAPS, penalised weighted median) and graceful degradation on missing data/failed API calls (§6.1, Finding I-4).
- Confidence-interval and odds-ratio conversion arithmetic (§11 checklist).
- Download-matches-display consistency for the primary results export (§6.3).

### What is partially implemented
- Horizontal-pleiotropy / exclusion-restriction checking — only for CpGs reaching the ≥3-instrument tier (§7.3); genuinely disclosed, not hidden.
- Genome-build verification for uploaded data — a disclosed best-effort guess, not a hard check (§9).
- MR-PRESSO — fully computed, but with no user-facing output (§10, Finding M-1).

### Potential scientific limitations (inherited from the underlying study design, not code defects)
- The Preloaded route's outcome GWAS is not sex-stratified, though the candidate CpG panel was derived from sex-stratified evidence — a scope limitation documented in the pipeline's own write-up (`METHODS_mendelian_randomization.md` §2.FF.4) and unchanged by this module.
- Most CpGs yield only 1–2 independent instruments, limiting statistical power and the applicability of pleiotropy sensitivity checks — an inherent property of mQTL data density, not a code defect.

### Confirmed implementation issues
- **H-1** (min-p multiple-testing selection, §10) — the most consequential finding in this audit.
- **M-1** through **M-5** (§10) — MR-PRESSO not surfaced; duplicated Forest/Single-SNP-forest tabs; cached per-instrument Steiger flags unused in favour of a different aggregate statistic under a similar label; CpG-picker caption overstates the panel it lists; Preloaded-route clumping count display can diverge from the instruments actually used.
- **L-1** through **L-4** (§10) — minor, narrow-impact implementation notes.

### Recommended improvements (documentation only — no code changed as part of this audit)
1. Carry an explicit `primary` method flag through `mr_state()$results` and correct on that flag alone (resolves H-1).
2. Either surface MR-PRESSO's `global_p` (and ideally outlier/corrected-estimate fields) in a table, or remove the checkbox until it is wired to an output (resolves M-1).
3. Differentiate the "Single-SNP forest" tab from "Forest" (e.g. per-SNP rows only, no pooled-method rows) or merge the two tabs (resolves M-2).
4. Either surface the cached per-instrument `steiger_dir`/`steiger_pval` columns explicitly (e.g. as an additional table) or rename the live output to make clear it is a different, aggregate-level statistic (resolves M-3).
5. Correct the Preloaded CpG-picker's caption to state that it lists only the majority-vote CpGs with a usable instrument, and disclose the count of majority-vote CpGs excluded for lacking one (resolves M-4).
6. Base the Preloaded-route "after clumping" count on `nrow(d)` rather than the original `instrument_counts.csv` lookup, so it always reflects instruments actually reaching Harmonisation (resolves M-5).

### Overall assessment
**Mostly coherent, with minor limitations.** The module's core MR machinery — instrument selection, clumping, harmonisation, tiered estimation, and instrument-count-gated sensitivity analysis — is scientifically appropriate, closely and verifiably faithful to its own documented source pipeline, and consistently honest in-UI about where sparse instrumentation limits what can be checked. The one methodologically consequential issue (H-1, multiple-testing correction via min-p selection) does not currently manifest in the bundled Preloaded dataset (no CpG there reaches the ≥3-instrument tier) but would affect any richer instrument panel or Upload-route dataset, and should be corrected or explicitly caveated before this module's adjusted-significance output is quoted in a thesis results section. The remaining confirmed issues are narrow in scope (one dormant computation, one duplicated plot tab, two UI-text/display discrepancies) and do not undermine the module's basic scientific validity.
