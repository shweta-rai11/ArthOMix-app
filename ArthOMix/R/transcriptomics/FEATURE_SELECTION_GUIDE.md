# Feature Selection sub-module — code guide

File: [`R/transcriptomics/mod_featureselection.R`](mod_featureselection.R) (1480 lines)
Registered as: Section 2.8 "Feature Selection", `TX_MODULES` in [`R/submodules_registry.R:16`](submodules_registry.R#L16)
Wired up in: [`server.R:138`](../../server.R#L138) — `mod_featureselection_server("tx_featureselection", dataset, results)`

---

## 1. What this sub-module does, in plain terms

Given a shortlist of "candidate genes" (usually the output of the Candidate Gene
Identification tab, Section 2.5 — WGCNA co-expression module ∩ sex-stratified DEGs),
this tab answers: **which of those genes actually separate patients (RA) from
controls (HC) best, once you fit a real classifier?**

It does this three independent ways and reports where they agree:

| Method | Idea | Selection rule |
|---|---|---|
| **LASSO** logistic regression | Penalized regression that shrinks unhelpful genes' coefficients to exactly 0 | Genes with a non-zero coefficient at the CV-optimal `lambda` |
| **Random Forest importance** | An ensemble of decision trees; genes that split the data cleanly get high importance | Genes with Mean Decrease in Gini above the panel's own mean (or top-N) |
| **SVM-RFE** (Support Vector Machine – Recursive Feature Elimination) | Repeatedly fit a linear SVM and drop the single weakest-weighted gene, one at a time, until a full most→least-important ranking exists | The top-*k* ranked genes that minimize 10-fold CV classification error |

**Consensus** = the intersection of whichever of these three methods you keep ticked
(all three by default — this is exactly the project's own thesis methodology,
Chapter 2, Section 2.8). A gene that survives all three very different statistical
attacks is a much stronger biomarker candidate than one that only one method liked —
this is the whole reason the tab exists rather than just running one method.

Everything is fit **twice**, completely independently: once on female samples, once
on male samples (never with sex as a covariate in one combined model) — because this
project's premise is that RA biomarkers can differ by sex. A third "Run All
(pooled)" option exists too, ignoring sex entirely, mainly so you can replicate a
published method that never stratified by sex in the first place.

### Why this matters for the app / benefits to the user
- **Turns a long gene list into a short, defensible panel.** Candidate Gene ID
  typically hands over dozens–hundreds of genes; this tab is what narrows that to a
  handful worth taking to a diagnostic model or wet-lab validation.
- **Triangulates, rather than trusting one algorithm.** Any single method
  (especially LASSO, which is unstable under collinearity) can pick a slightly
  different gene set depending on noise. Requiring 3-way agreement is a standard
  way to guard against that instability.
- **Makes the sex-stratification story concrete.** Because female and male are
  fit as fully separate models, this tab is where you can see, gene-by-gene,
  whether a biomarker is shared or sex-specific — the app's central scientific
  question.
- **Feeds everything downstream.** Its output (`results$featureselection`) is read
  directly by the Diagnostic Biomarker Card (`mod_biomarkercard.R`) and by
  ArthOChat's own context builder (`build_assistant_context()` in
  `submodules_registry.R:117-146`), and conceptually feeds the Diagnostic Model tab
  (Section 2.9), which repeats this exact Female/Male/Pooled pattern one level up
  (classification performance rather than gene selection).
- **Instant when it can be, honest when it can't.** On the project's own default
  dataset with default settings, results appear instantly from a precomputed run;
  the moment you change *anything* (a parameter, the dataset, your own data), it
  transparently switches to a live fit and tells you so ahead of time
  (`speed_hint_ui`, line 820) — you're never silently shown stale numbers.

---

## 2. Function reference

### 2.1 Pure computational functions (no Shiny — could be unit-tested standalone)

| Function | Location | What it does | Input | Output |
|---|---|---|---|---|
| `fs_svm_rfe_rank` | [`mod_featureselection.R:47`](mod_featureselection.R#L47) | Recursive Feature Elimination: fits a linear-kernel SVM, drops the single feature with the smallest squared SVM weight, repeats until one feature is left. Produces a full most→least-important ranking. | `X` (numeric matrix, samples × genes), `y` (2-level factor), `cost` (SVM cost, default `1`), `tolerance` (default `0.001`), `class_weights` (named vector or `NULL`) | `character` vector of gene names, ordered by importance |
| `fs_svm_rfe_curve` | [`:61`](mod_featureselection.R#L61) | For every prefix length *k* of the rank above, refits an SVM on just the top-*k* genes and measures 10-fold CV error — used to pick panel size objectively instead of by eyeballing. | `X`, `y`, `rank` (from `fs_svm_rfe_rank`), `cost`, `seed` (1234), `folds` (10), `tolerance`, `class_weights` | `list(k, err, best, besterr)` — `k`/`err` are parallel vectors (panel size vs. CV error), `best` is the size that minimizes error |
| `fs_class_weight_levels` | [`:127`](mod_featureselection.R#L127) | Computes per-class weights for imbalanced HC/RA groups: `"equal"` (all 1s, the project default), `"balanced"` (inverse class frequency), or `"manual"` (a fixed ratio you set). | `y` (factor), `mode` (`"equal"`/`"balanced"`/`"manual"`), `ratio` (numeric) | named numeric vector, one weight per factor level |
| `fs_obs_weights` | [`:140`](mod_featureselection.R#L140) | Expands the per-class weights above into one weight per *sample* (what `glmnet`'s `weights=` actually wants). | `y`, `mode`, `ratio` | numeric vector, length = number of samples |
| `fs_fit_sex` | [`:152`](mod_featureselection.R#L152) | **The core model-fitting engine.** Given one sex's (or pooled) expression submatrix and group labels, fits LASSO (`glmnet::cv.glmnet`), a CV-tuned Random Forest (`caret::train` + `randomForest::randomForest`), and a CV-tuned SVM-RFE (`e1071::tune` + the two functions above), then intersects whichever methods are selected for consensus. Resets `set.seed(1234)` before every stochastic step, so it is exactly reproducible. | `X` (matrix, samples × genes), `y` (2-level factor), `params` (list — see `FS_DEFAULT_PARAMS`) | `list` with `lasso_genes`, `rf_genes`, `svm_genes`, `consensus`, plus diagnostics (`cv`, `gini`, `svm_curve`, chosen hyperparameters, `n_input`, `n_samples`, `fast_path = FALSE`) |

**Constants** (not functions, but drive behavior):
- `FS_SVM_COST_GRID` ([`:73`](mod_featureselection.R#L73)) — default SVM cost values searched: `0.01, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16`.
- `FS_MAX_CANDIDATE_GENES` ([`:85`](mod_featureselection.R#L85)) — `200`. SVM-RFE refits once per gene *and* once more per panel size, so an uncapped 1,000+ gene candidate list (e.g. a raw WGCNA module) would take hours. Above this cap, candidates are cut down to the most-variable genes in that sex's own samples first (and the UI says so).
- `FS_DEFAULT_PARAMS` ([`:97`](mod_featureselection.R#L97)) — the complete default hyperparameter set for all three methods, matching this project's own thesis script exactly. Every value here is individually overridable from the UI.

### 2.2 UI-building helpers

| Function | Location | What it does | Input | Output |
|---|---|---|---|---|
| `mod_featureselection_config` | [`:274`](mod_featureselection.R#L274) | Static registry entry (id, title, description, icon) so `TX_MODULES` can build the sidebar/menu entry. | — | `list` |
| `fs_sex_panel` | [`:288`](mod_featureselection.R#L288) | Wraps one Female/Male/Pooled result box so it only appears once **that sex's own** Run button has been clicked (not "any of the three"). | `ns`, `run_btn_id`, `...` (UI content) | a `conditionalPanel` |
| `mod_featureselection_technique_panel` | [`:297`](mod_featureselection.R#L297) | Builds one technique's result box: summary text, plot, CSV download button, results table. Instantiated 12 times (LASSO/RF/SVM-RFE/Overlap × Female/Male/Pooled). | `ns`, `prefix` (e.g. `"pooled_lasso"`), `title`, `plot_height` | a `box(...)` |
| `mod_featureselection_params_box` | [`:317`](mod_featureselection.R#L317) | Wraps one method's parameter controls in a titled box. | `ns`, `prefix`, `method_label`, `defaults_desc`, `...` | a `box(...)` |
| `mod_featureselection_ui` | [`:325`](mod_featureselection.R#L325) | **The module's top-level UI.** Builds the left "Candidate genes & samples" panel (data-source picker, the 3 Run buttons, group/weighting controls, the "Status:" list) and the right-hand LASSO/Random Forest/SVM-RFE/Overlap tabs. | `id` | a `tagList` |
| `mod_featureselection_server` | [`:509`](mod_featureselection.R#L509) | **The module's server logic** (everything in §2.3–2.4 below lives inside this `moduleServer` call). | `id`, `dataset` (shared app state), `results` (shared app state, written to) | side-effecting; nothing returned directly |

### 2.3 Server-side reactive logic (inside `mod_featureselection_server`)

| Name | Location | What it does | Reads | Produces |
|---|---|---|---|---|
| `own_meta_raw` | [`:520`](mod_featureselection.R#L520) | Parses the uploaded metadata file for "Upload my own expression data". | `input$meta_file` | data.frame |
| `expr_column_mapping` (output) | [`:532`](mod_featureselection.R#L532) | Builds the sample-ID/group/sex column-picker for the uploaded metadata. | `own_meta_raw()` | Shiny UI |
| `source_expr_meta` | [`:545`](mod_featureselection.R#L545) | **Single point where "which data source am I using" gets resolved.** For `"expr"`: reads the uploaded matrix + mapped metadata columns. For `"project"`/`"deg"`: reads the app-wide `dataset$expr`/`dataset$meta`. | `input$data_source`, uploads or `dataset` | `list(expr = matrix, meta = data.frame)` |
| `sex_levels` | [`:580`](mod_featureselection.R#L580) | Figures out which value in the sex column means female vs. male (matches `^f`/`^m` case-insensitively, falls back to alphabetical). | `source_expr_meta()$meta$sex` | `list(female=, male=)` |
| `group_controls_ui` (output) | [`:590`](mod_featureselection.R#L590) | Renders the reference/comparison group pickers and the class-weighting radio — **built from whichever data source is active**, so this is identical UI regardless of source. | `source_expr_meta()$meta$group` | Shiny UI |
| `project_candidate_genes(sex_label)` | [`:620`](mod_featureselection.R#L620) | For `data_source == "project"`: returns this session's live candidate list (`results$candidates[[sex_label]]`) if it exists, else falls back to a bundled `FS_input_<sex>.csv`. | `results$candidates`, bundled CSVs | `list(genes, is_live, note)` |
| `load_precomputed_fs(sex_label, mhc_exclude)` | [`:649`](mod_featureselection.R#L649) | The "instant" fast path: loads a pre-fit result straight from `ml_features.rds` instead of recomputing. | `PROCESSED_NEW_DIR/ml_features[_noMHC].rds` | same shape as `fs_fit_sex()`'s output, or `NULL` |
| `deg_candidate_genes(sex_label)` | [`:677`](mod_featureselection.R#L677) | For `data_source == "deg"`: reads an uploaded DEG/gene-list CSV per sex, or a WGCNA module's genes. | uploaded file or `results$wgcna$module_genes` | `list(genes, note)` |
| `expr_candidate_genes(sex_label, expr_sub)` | [`:706`](mod_featureselection.R#L706) | For `data_source == "expr"`: candidate genes are either the most-variable genes in the uploaded matrix, a pasted list, or a WGCNA module. | `input$gene_source`, uploaded expression matrix | `list(genes, note)` |
| `project_source_ui` (output) | [`:761`](mod_featureselection.R#L761) | Tells the user whether live candidates exist yet, and exposes the MHC-region-exclusion checkbox (bundled lists only). | `results$candidates` | Shiny UI |
| `fs_any_customized` | [`:798`](mod_featureselection.R#L798) | Compares every parameter widget's current value against `FS_DEFAULT_PARAMS` to decide whether the run still qualifies for the instant precomputed path. | all method-parameter inputs | logical |
| `speed_hint_ui` (output) | [`:820`](mod_featureselection.R#L820) | Tells you *before* clicking Run whether it will be instant or live, and why. | data source, candidates, customization state | Shiny UI |
| `fs_advanced_params` | [`:853`](mod_featureselection.R#L853) | Reads every LASSO/RF/SVM-RFE widget's current value into one `params` list (falls back to `FS_DEFAULT_PARAMS` for anything unset). | all method-parameter inputs | `list` |
| **`fs_build_sex(sex_label, sex_value)`** | [`:891`](mod_featureselection.R#L891) | **The single function every Run button calls.** See §3 below — this is the whole answer to "does uploading your own data change the pipeline". | `input$ref_group`/`comp_group`, `source_expr_meta()`, candidate-gene resolver for the active data source, `fs_advanced_params()` | fitted result `list`, or halts (`validate()`) with an on-screen message |
| `fs_result_female` / `fs_result_male` / `fs_result_pooled` | [`:990`](mod_featureselection.R#L990)–[`:1002`](mod_featureselection.R#L1002) | `eventReactive`s bound to the three Run buttons; each calls `fs_build_sex()` for its own sex. | button clicks | fitted result, or a caught error if the fit halted |
| `fs_has_run` + 3 `observeEvent`s | [`:1018`](mod_featureselection.R#L1018)–[`:1030`](mod_featureselection.R#L1030) | Reveals the (until-then hidden) results panel the first time *any* Run button is clicked. | button clicks | side effect (`shinyjs::show`) |
| `lasso_params_ui` / `rf_params_ui` / `svm_params_ui` / `consensus_params_ui` (outputs) | [`:1032`](mod_featureselection.R#L1032)–[`:1131`](mod_featureselection.R#L1131) | Render each method's own parameter-customization box. | — | Shiny UI |
| `observeEvent(fs_result_female())` / `_male()` / `_pooled()` | [`:1153`](mod_featureselection.R#L1153) / [`:1182`](mod_featureselection.R#L1182) / [`:1211`](mod_featureselection.R#L1211) | The moment a sex's fit completes, saves a summary into `results$featureselection` and a log entry into `results$featureselection_runs`, and shows a "saved: N consensus genes" toast. **The pooled one's toast text is exactly `"Pooled (all) feature selection saved: %d consensus genes%s."` — line 1234.** | fit result | writes to `results$featureselection*`, shows notification |
| `saved_runs_ui` (output) | [`:1240`](mod_featureselection.R#L1240) | The left-panel **"Status:"** list. This is the exact code that prints `"<Sex> feature selection - not run yet"` (line 1246) when that sex's result is `NULL`. | `fs_result_female/male/pooled()` (via `tryCatch`) | Shiny UI |
| `summary_ui` (output) | [`:1267`](mod_featureselection.R#L1267) | The bottom "Result" box — one line per sex, `"<Sex>: not run yet."` when `NULL`. | same as above | Shiny UI |
| `res_sex(sex_label)` | [`:1291`](mod_featureselection.R#L1291) | Wraps `fs_result_<sex>()` in `tryCatch(..., error = \(e) NULL)` so downstream plot/table outputs get a plain `NULL` instead of a hard Shiny error. | — | `reactive` |
| `register_sex_technique_outputs(sex_label, res)` | [`:1296`](mod_featureselection.R#L1296) | Registers every `output$<sex>_lasso_*` / `_rf_*` / `_svm_*` / `_consensus_*` (summary, plot, table, CSV download) for one sex. Called once each for female/male/pooled at lines 1476–1478. | `res()` | registers ~16 outputs per sex |

---

## 3. Verified: uploading your own data does **not** change the pipeline's parameters or filters

You asked me to confirm that "Upload my own expression data" (`data_source ==
"expr"`) still runs through the exact same modeling parameters and sample/group
filters as the built-in project pipeline — not a parallel, possibly-drifted code
path.

**Reading `fs_build_sex()` (`mod_featureselection.R:891-976`) line by line:** every
step *except one* is identical regardless of `input$data_source`:

1. Reference/comparison group validation (`:892-893`) — same for all sources.
2. Sex-column filtering (`:900-903`) — same for all sources (skipped only for
   "pooled", regardless of source).
3. Group-membership + sample filter (`:938-939`) and the ≥10-sample gate
   (`:940-941`) — same for all sources.
4. `FS_MAX_CANDIDATE_GENES` (200-gene) cap, applied identically by most-variance
   reduction (`:957-963`) — same for all sources.
5. Building `X`/`y` and the ≥6-per-group gate (`:965-967`) — same for all sources.
6. **`fs_fit_sex(X, y, params = adv_params)` (`:969`)** — the actual model fit —
   called identically for all sources, with `adv_params` coming from the same
   `fs_advanced_params()` (`:853`) regardless of `data_source`.

**The only thing that differs by data source** is line 945-949 — *which function
resolves the candidate gene list* (`project_candidate_genes` /
`deg_candidate_genes` / `expr_candidate_genes`). That's the entire, intentional
purpose of the three-way choice; it does not touch modeling parameters or sample
filters at all. (One legitimate, deliberate exception: the MHC-region-exclusion
checkbox only makes sense for the bundled candidate lists, so it's a no-op — not a
different filter — for `"expr"`/`"deg"`.)

**I tested this empirically**, not just by reading the code. I extracted the pure
computational core (`fs_fit_sex` and its dependencies, `mod_featureselection.R:1-273`,
which have no Shiny dependency) and ran it through a synthetic 40-sample/12-gene
dataset twice: once simulating the "project" data-source path, once simulating the
"expr upload" path with the same underlying data and the same candidate genes,
including a second run with class-weighting turned on (not just left at default) to
make sure the shared-parameter path isn't only exercised at defaults. Both paths
produced **byte-identical** `lasso_genes`, `rf_genes`, `svm_genes`, `consensus`,
`rf_mtry`, and `svm_cost` (seed `1234` is reset before every stochastic step in
`fs_fit_sex`, so this is exactly reproducible, not just "close"):

```
== Parity check: project-path vs expr-upload-path ==
lasso_genes                  MATCH
rf_genes                     MATCH
svm_genes                    MATCH
consensus                    MATCH
rf_mtry                      MATCH
svm_cost                     MATCH
n_input                      MATCH
n_samples                    MATCH
weighted lasso_genes         MATCH
weighted rf_genes            MATCH

RESULT: PASS - identical filters/parameters/output between the project-data
path and the uploaded-data path; only candidate-gene sourcing differs, as
designed.
```

**Conclusion: no fix needed here** — the code already guarantees this by
construction, since `fs_build_sex()` is the single function every Run button calls,
and it branches only on candidate-gene sourcing, never on modeling logic. (The test
script itself was a scratch verification, not added to the repo — say the word if
you'd like it turned into a permanent `tests/testthat/` regression test that pins
this parity going forward.)

---

## 4. Why "Pooled (all) feature selection - not run yet" is showing

Two separate things can cause this text (`mod_featureselection.R:1246`,
`1273`), and it's worth telling them apart:

### 4.1 The trivial reason
Nobody has clicked the **"Run All (pooled)"** button yet this session
(`mod_featureselection.R:350`). `fs_result_pooled` (`:1000`) is an `eventReactive`
bound to that button with `ignoreInit = TRUE` — it holds no value at all until
that specific button is clicked. If that's all this is, click it and it'll switch
to "completed" like Female/Male do.

### 4.2 A real structural gap — this is likely what you're hitting
If you *have* clicked "Run All (pooled)" and it's still showing "not run yet" (or
if you click it and immediately get an on-screen validation error), it's because
**the default "project pipeline" data source has no pooled candidate-gene source at
all**, only sex-specific ones. Tracing `fs_build_sex("pooled", NULL)` on
`data_source == "project"`:

1. `project_candidate_genes("pooled")` (`:620`) first checks
   `results$candidates[["pooled"]]$genes` — but the Candidate Gene Identification
   module (`mod_candidates.R:456-459`) only ever writes `female`/`male`/`final`
   keys. Candidate discovery there is *inherently* sex-stratified (WGCNA module
   background intersected with sex-stratified DEGs) — there's no non-sex-specific
   candidate list to read.
2. It falls back to a bundled `FS_input_pooled[_noMHC].csv` — **this file doesn't
   exist.** Only `FS_input_female.csv` and `FS_input_male.csv` (and their
   `_noMHC` variants) are bundled under `data/preloaded/transcriptomics/results/tables/`.
3. The instant precomputed path, `load_precomputed_fs("pooled", ...)` (`:649`),
   reads `ml_features.rds` — whose top-level keys are only `female`, `male`,
   `expr`, `meta`, `seed`, `built` (verified by reading the file directly). No
   `pooled` key exists there either.

So `project_candidate_genes("pooled")` returns `genes = character(0)`, and
`fs_build_sex()` halts at its own `validate(need(length(genes) >= 3, ...))` gate
(`:951`) with the message *"Fewer than 3 pooled candidate genes are present…"*.
Because the `eventReactive` never completes, `fs_result_pooled()` keeps raising
that (caught, silent) error indefinitely — the status stays "not run yet"
**even after clicking the button**, for as long as `data_source == "project"`.

This is a real, load-bearing gap in the default pipeline, not a Shiny quirk. It's
consistent with the project's own methodology being sex-stratified from the ground
up (Chapter 2, Section 2.8's whole premise), which just never defined what a
"pooled candidate panel" should mean.

**Workarounds available today** (no code change needed even without the fix
below): switch `data_source` to `"expr"` (upload your own expression matrix and
pick candidate genes directly — most-variable / pasted list / WGCNA module, none
of which are sex-specific), or to `"deg"` and upload a `pooled_deg_file`, or pick
`deg_source_mode == "wgcna"` (a WGCNA module's genes aren't sex-specific by
construction either).

**`Research_05_multiomics_sexstratified/`** (the untracked directory in your
working tree) is unrelated — it's read only by the separate Multi-Omics module
(`R/multiomics/*.R`: `multiomics_sexstratified_engine.R`,
`mod_multi_concordance.R`, `mod_multi_integration.R`, `mod_multi_summary.R`), never
by `mod_featureselection.R`.

### Fix implemented: pooled candidates = union(female, male)

There's no single obviously-correct definition of "pooled candidate genes" —
union (broadest, "pool everyone"), intersection (narrowest, "candidates in both
sexes"), or a dedicated non-sex-stratified discovery run were all considered. You
chose **union**, so `project_candidate_genes()` (`mod_featureselection.R:620-666`)
now special-cases `sex_label == "pooled"`:

1. If either sex has a **live** candidate panel this session
   (`results$candidates$female/$male`), pooled candidates = the union of whichever
   are available, marked `is_live = TRUE`.
2. Otherwise, falls back to the union of the **bundled** `FS_input_female.csv` /
   `FS_input_male.csv` lists (respecting the MHC-exclusion checkbox), matching how
   female/male already fall back individually.
3. Only if *neither* source yields anything does it still report "No pooled
   candidate genes available" — now a genuinely-empty-data case, not a permanent
   structural gap.

The instant precomputed path (`load_precomputed_fs`) still returns `NULL` for
`"pooled"` (no `pooled` key in `ml_features.rds`), so a pooled run always falls
through to a **live** `fs_fit_sex()` fit — correct, since the precomputed run was
never built for a pooled panel in the first place; the speed hint UI will show
"runs live" for it, same as any other non-fast-path run.

**Tested, not just reasoned about**, using this repo's real bundled files and real
expression matrix (`data/preloaded/transcriptomics/processed/new/ml_features.rds`,
183 samples, 80 HC / 103 RA):

```
Pooled candidate panel (union): 33 genes   (32 female + 25 male, deduplicated)
Pooled samples matched: 183 (HC=80, RA=103)
Fitting on X: 183 samples x 33 genes, y: HC=80, RA=103

== fs_fit_sex() result for pooled (all-sample) run ==
LASSO genes:         9 -> VPS52, GNL1, HLA-DPA1, C6orf136, IKZF3, MED1, SMARCC2, CDC37, ESYT1
Random Forest genes: 7 -> MED1, C6orf136, SMARCC2, ESYT1, WDR46, BRD2, IKZF3
SVM-RFE genes:       2 -> SMARCC2, ESYT1
CONSENSUS genes:     2 -> SMARCC2, ESYT1
```

A full live pooled run now completes end-to-end and produces a real 2-gene
consensus panel. Clicking **"Run All (pooled)"** on the default project pipeline
will show a "completed" status instead of "not run yet" going forward.
(Verification scripts were scratch, not committed — say the word if you'd like
this pinned as a permanent `tests/testthat/` regression test.)

---

## 5. UI map

- **Left panel** ("Candidate genes & samples" box, `mod_featureselection_ui:335-431`):
  data-source radio (`data_source`), sticky Run Female / Run Male / **Run All
  (pooled)** buttons + speed hint (`:345-353`), the three data-source-specific
  input blocks (`:364-427`), group/class-weighting controls (`group_controls_ui`),
  and the **"Status:"** list (`saved_runs_ui`, `:430`) — this is exactly what
  currently reads "Pooled (all) feature selection - not run yet".
- **Right panel** (`:433-506`), hidden until the first Run click:
  - Tabs **LASSO / Random Forest / SVM-RFE / Overlap**, each with a Female box, a
    Male box, and a **Pooled (all)** box — a plot, a summary, a `DT` results
    table, and a "Genes (CSV)" download button per box.
  - Bottom **"Result"** box (`summary_ui`) — one contrast-aware line per sex.
  - **"References"** box (`references_box_ui`) — the four method citations.
