# 03. Functions and Code Audit — Methylomics Cell-Type Deconvolution

Audit assessment classes used below: **Correct** (does what it claims, no caveats found), **Reasonable** (sound design choice, minor caveats noted), **Potential concern** (works but has an edge case worth watching), **Bug** (produces a wrong result under conditions that can occur in normal use — cross-referenced with `07_code_audit_findings.md`), **Unused** (declared but never called/read), **Redundant** (duplicates logic that exists elsewhere in the app, noted as a deliberate isolation choice), **Fragile** (works today but relies on an assumption that could silently break), **Missing validation**, **Scientifically questionable**.

---

## Part A — `celltype.R` (36 functions)

### A.1 Reference library registry

| Function | Input | Purpose / What it does | Output | Downstream use | Audit |
|---|---|---|---|---|---|
| `methyl_ct_reference_registry()` (`celltype.R:27-54`) | none | Lists the 7 EpiDISH built-in reference panels; for each, resolves the actual matrix via `getExportedValue("EpiDISH", s$object)` and reads `nrow`/`colnames` off the real object rather than hard-coding them | `list` of 7 specs, each with `id`/`label`/`object`/`tissue`/`available`/`ncpg`/`celltypes` | Populates the `ct_ref_id` dropdown and blood-only `ct_hepidish_ref2` dropdown | **Correct.** The choice of `getExportedValue()` over `get(..., envir = asNamespace(...))` is verified sound: `EpiDISH`'s centroid matrices are lazy-loaded data objects, and `getExportedValue()` is what the `::` operator itself calls to resolve lazy data without requiring `library(EpiDISH)`. Independently confirmed against the installed package in this environment: all 7 objects (`centDHSbloodDMC.m`, `centBloodSub.m`, `cent12CT.m`, `cent12CT450k.m`, `centEpiFibIC.m`, `centEpiFibFatIC.m`, `centCAB100i.m`) resolve successfully with `nrow`/`ncol`/`colnames` matching what this function would report — 7, 7, 12, 12, 3, 4, and 19 cell types respectively (a fixed property of the installed package version, not user data). |
| `methyl_ct_get_reference(ref_id)` (`celltype.R:56-61`) | `ref_id` string | Looks up one spec by id and re-resolves its matrix, returning `NULL` if unavailable | matrix or `NULL` | `active_reference_full()`, hepidish's second-stage reference lookup | **Correct.** Re-derives the registry on every call rather than caching it — harmless given the registry itself is cheap (7 lazy-data lookups), not a performance concern at this scale. |
| `methyl_ct_unavailable_methods()` (`celltype.R:68-77`) | none | Static list of MethylResolver / IDOL / reference-free, each with a stated installation reason | `list` of 3 | Rendered as a disabled-methods note on Tab 3 | **Correct.** Matches the "explain why, don't fake it" convention this app uses elsewhere (cited against `qc.R`'s `methyl_filter_cross_reactive()`/`methyl_filter_maf()`, verified to follow the identical pattern). |

### A.2 Scale detection, transforms, dataset summary

| Function | Input | Purpose | Output | Downstream use | Audit |
|---|---|---|---|---|---|
| `methyl_ct_detect_scale(mat)` (`celltype.R:90-104`) | numeric matrix | Classifies beta (0–1) vs. percent (0–100) vs. M-value (unbounded) from a sampled 0.1/50/99.9 percentile range (capped at 20,000 sampled values) | `list(scale, note)` | `ct_scale_detect()`, own-upload load handler | **Reasonable.** A genuine 3-way extension of `mod_methyl_featureselection.R`'s 2-way `methyl_fs_detect_scale()` (verified: that function's own quantile logic only distinguishes beta vs. M, `mod_methyl_featureselection.R:60-68`). Boundary tolerance (`-0.05`/`1.05`) matches the convention used elsewhere in the same file. Never transforms silently — only ever returns a classification, consumed by an explicit-click transform handler. |
| `methyl_ct_pct_to_beta(mat)` (`celltype.R:106`) | matrix | `mat / 100` | matrix | Transform handler | **Correct**, trivial. |
| `methyl_ct_m_to_beta(m)` (`celltype.R:109`) | matrix | `2^m / (1 + 2^m)`, the inverse logit | matrix | Transform handler; also applied to `dataset$beta` when `dataset$input_scale == "m"` (`ct_source()`) | **Correct** — exact algebraic inverse of `qc.R`'s `methyl_beta_to_mvalue()` (`log2(b/(1-b))`, clipped at `eps=1e-4`). No clipping needed on this direction since the sigmoid is bounded by construction. |
| `methyl_ct_working_summary(mat, array_type)` (`celltype.R:115-128`) | matrix, optional array type | Compact summary: n CpG, n sample, % missing, beta range, duplicated-ID count, chromosome count | `list` | `ct_source_summary_ui` | **Correct.** Chromosome count degrades to `NA` (not an error) when `methyl_get_annotation()` reports no manifest for the array type — matches the app-wide convention. |
| `methyl_ct_parse_custom_reference(datapath, filename)` (`celltype.R:135-143`) | file path/name | Parses via `methyl_parse_matrix()` (same CpG-rows/sample-columns parser, cell types instead of samples), then validates finite, `[-0.05, 1.05]`, ≥2 columns | `list(ok, mat)` or `list(ok=FALSE, error)` | `ct_custom_ref_raw()` | **Correct**, sound reuse of the existing matrix parser rather than a bespoke one; validation bounds are intentionally slightly wider than `[0,1]` to tolerate floating-point noise, consistent with how the rest of the app treats beta ranges. |

### A.3 Marker-CpG ranking (no fabricated p-values)

| Function | Input | Purpose | Output | Downstream use | Audit |
|---|---|---|---|---|---|
| `methyl_ct_marker_rank(ref_mat)` (`celltype.R:160-193`) | reference matrix (CpG × cell type) | For every CpG: the cell type whose centroid is most extreme vs. the other types' max/min (`eff_hyper`/`eff_hypo`), an `effect` size (the larger of the two directions), a `direction` label, a `specificity` score (effect ÷ SD of the other types' centroids, floor `1e-6` to avoid divide-by-zero), and `btw_type_var` (across-cell-type variance via `stats::var`) | `data.frame`, one row per CpG | `fs_result()`, and as the pre-feature-selection fallback marker set | **Reasonable / scientifically honest by design.** No p-value or FDR column is produced, and the code comment explains why: reference centroids are single mean-beta values per cell type with **no per-sample replicates**, so no real significance test is available — computing one would require fabricating a sample-size assumption. This is verified consistent with how the UI actually renders the result (`ct_fs_table` hard-codes `"n/a"`/`"n/a (no replicates)"` rather than computing anything, `mod_methyl_celltype.R:626`). See `05_statistical_methodology.md` for the underlying math. One subtlety: `specificity` uses the SD of only the *other* cell types' centroids (not the target's own value), so a reference with exactly 2 cell types has an "others" set of size 1, whose `sd()` is `NA` — this is guarded (`osd[!is.finite(osd) | osd == 0] <- 1e-6`, `celltype.R:173`), so specificity degrades to `effect / 1e-6` (a very large number) rather than erroring, for a 2-cell-type reference. Not incorrect, but the resulting specificity scale is not comparable across references with different cell-type counts — worth being aware of if specificity values are ever compared across runs with different references. |
| `methyl_ct_select_markers(rank_df, dbeta_min, effect_min, direction, specificity_mode, chr_allowed_ids)` (`celltype.R:203-217`) | ranked data.frame + filter params | Applies effect-size, direction, specificity (median-split), and chromosome-scope filters as a combined logical mask | filtered `data.frame` | `fs_result()` | **Correct.** Note `dbeta_min` and `effect_min` are both compared against the same `df$effect` column (`celltype.R:207-208`) — in practice the UI only ever sets one of them (`ct_fs_dbeta` maps to `dbeta_min`; `effect_min` is never passed a non-zero value anywhere in `mod_methyl_celltype.R`), so this redundancy is inert today but is a slightly confusing function signature (two parameters that currently always filter on the identical column). |
| `methyl_ct_top_n_balanced(df, sort_col, top_n)` (`celltype.R:225-239`) | ranked/filtered data.frame, sort column, cap | Caps to `top_n` while balancing roughly evenly across cell types (`ceiling(top_n / n_types)` per type, then trims global overshoot from the weakest picks) | `data.frame` | `fs_result()` | **Correct and well-designed.** Directly verified against its stated purpose: without this balancing, a top-N-by-effect-size cut would be dominated by whichever cell type has the largest overall separation from the others, silently starving markers for the remaining cell types — the per-type pre-allocation plus global trim avoids that. |
| `methyl_ct_chr_allowed_ids(cpg_ids, array_type, scope)` (`celltype.R:246-260`) | CpG IDs, array type, scope enum | Restricts to autosomes / autosomes+X / all, via `methyl_get_annotation()`; **unresolved CpGs (no manifest match) are kept, not dropped** | `list(ids, note)` | `fs_result()`'s `ct_fs_chr_scope` filter | **Correct**, and its "unresolved = kept" convention is confirmed to match `qc.R`'s `methyl_filter_maf()`'s identical treatment of unannotated probes (absence of a manifest hit is not evidence the probe sits on a sex chromosome). |

### A.4 Reference / working-matrix overlap QC

| Function | Input | Purpose | Output | Downstream use | Audit |
|---|---|---|---|---|---|
| `methyl_ct_overlap_qc(marker_ids, working_ids)` (`celltype.R:266-273`) | marker CpG IDs, working-matrix row names | Set-intersection overlap stats | `list(n_ref, n_matched, n_missing, pct_matched, matched, missing)` | `ct_overlap()`, `ct_overlap_ok()` (the Run-button gate) | **Correct**, simple set logic. Important scope caveat (see `07_code_audit_findings.md`): this checks CpG-**ID presence only**, not per-sample completeness — a marker CpG counted as "matched" here can still be dropped later inside `methyl_ct_run_epidish()`'s `stats::complete.cases()` step if it has even one `NA` among the retained samples. |
| `methyl_ct_overlap_by_type(marker_df, working_ids)` (`celltype.R:275-285`) | ranked/selected marker data.frame, working IDs | Same overlap stat broken out per cell type | `data.frame` | `ct_refqc_by_type_table` | **Correct.** |

### A.5 EpiDISH / hepidish wrappers

| Function | Input | Purpose | Output | Downstream use | Audit |
|---|---|---|---|---|---|
| `methyl_ct_run_epidish(beta_mat, ref_mat, method, maxit, nu.v, constraint)` (`celltype.R:299-319`) | filtered beta matrix, reference matrix, method/solver params | Thin validated wrapper around `EpiDISH::epidish()`: intersects CpG IDs (requires ≥10), drops any row with an `NA` in **any** sample via `stats::complete.cases()` (requires ≥10 remaining), then calls `epidish()` inside `tryCatch()` | `list(ok, fractions, method, n_markers_used, ref_used)` or `list(ok=FALSE, reason)` | `decon_result()`, `methyl_ct_compare_methods()` | **Correct, with one downstream-facing caveat.** Verified directly against the installed `EpiDISH::epidish()` source: it dispatches to `DoCP`/`DoRPC`/`DoCBS`, all deterministic (see `05_statistical_methodology.md`); the `constraint` argument is validated (`match.arg`) regardless of method but is **only actually consulted inside the `CP` branch** of `epidish()` itself — for `RPC`/`CBS` it is silently ignored by `EpiDISH`'s own dispatch logic (not a bug in this wrapper, but worth knowing when reading the Advanced Parameters panel; see `07`). The comment's claim that estimates are "already non-negative/sum-to-one for all three methods" is accurate for the *output* (verified: `DoRPC`/`DoCBS` explicitly clip negative coefficients and renormalize to sum 1 in their own code; `DoCP` under the `equality` constraint sums to exactly 1 by construction, and under `inequality` guarantees non-negativity and sum ≤ 1, not necessarily sum = 1). The claim of ">0.99 correlation on a synthetic known-mixture test" in the comment is an internal verification note this documentation set did not independently reproduce. |
| `methyl_ct_run_hepidish(beta_mat, ref1_mat, ref2_mat, ic_column, method, maxit, nu.v, constraint)` (`celltype.R:329-350`) | as above, plus a two-stage reference pair | Wrapper around `EpiDISH::hepidish()`: resolves `ic_column` to a column index, checks ≥10 CpG overlap for **both** stages independently, calls `hepidish()` in `tryCatch()` | `list(ok, fractions, method, n_markers_used)` (no `ref_used` field, unlike the single-stage wrapper — see `07`) | `decon_result()` (hepidish branch, always called with `method = "RPC"` hard-coded at the call site in `mod_methyl_celltype.R:718-720`, not by this function itself) | **Correct**, and matches `EpiDISH::hepidish()`'s verified source: it runs the chosen method independently on `ref1`/`ref2`, then replaces the `ic_column` fraction with `frac1[, ic_column] * frac2` (the top-level IC fraction multiplied into the second-stage sub-fractions) — exactly the "top-level split then further resolve one component" behavior the code comment describes. |

### A.6 Reconstruction validation

| Function | Input | Purpose | Output | Downstream use | Audit |
|---|---|---|---|---|---|
| `methyl_ct_reconstruct(ref_mat, fractions)` (`celltype.R:356-359`) | reference matrix, fractions matrix | `ref_mat[, ct] %*% t(fractions[, ct])`, restricted to the intersection of cell-type columns present in both | reconstructed beta matrix (CpG × sample) | `val_result()` | **Correct.** This is exactly the deconvolution mixture model run forward: reconstructed bulk methylation = weighted sum of reference centroids by estimated fraction. |
| `methyl_ct_validation_metrics(observed, reconstructed)` (`celltype.R:361-386`) | observed working matrix, reconstructed matrix | Aligns on common CpGs/samples, drops incomplete rows, computes overall Pearson correlation / RMSE / MAE / R² and the same 3 metrics per sample | `list(ok, overall, per_sample, observed, reconstructed, n_cpg, n_sample)` | `ct_val_result_ui`, `ct_val_plot`, `ct_val_table`, the export report | **Correct.** R² is computed as `1 - SS_res/SS_tot` against the observed matrix's own mean, the standard definition; correctly guards `<2` overlapping CpGs or `<1` sample. |

### A.7 Cross-method comparison

| Function | Input | Purpose | Output | Downstream use | Audit |
|---|---|---|---|---|---|
| `methyl_ct_compare_methods(beta_mat, ref_mat, methods, ...)` (`celltype.R:392-402`) | matrix, reference, method vector | Reruns `methyl_ct_run_epidish()` once per method on the identical inputs; requires ≥2 successful methods | `list(ok, fractions_by_method, failures)` | `cmpm_result()` | **Correct.** Since CP/RPC/CBS are all deterministic (verified in A.5), rerunning here reproduces the exact fractions Tab 4 would have gotten had the user picked that method originally — this is a meaningful, reproducible comparison, not an approximation. |
| `methyl_ct_method_correlation(fractions_by_method)` (`celltype.R:407-419`) | named list of fraction matrices | All-pairs Pearson correlation, flattened across samples × cell types | symmetric matrix | `ct_cmpmethods_corr_plot` | **Correct.** |
| `methyl_ct_method_agreement_summary(fractions_by_method)` (`celltype.R:422-433`) | same | Mean/max absolute difference per method pair | `data.frame` | `ct_cmpmethods_summary_table` | **Correct.** |

### A.8 Group / phenotype comparison

| Function | Input | Purpose | Output | Downstream use | Audit |
|---|---|---|---|---|---|
| `methyl_ct_group_stats(fractions, group)` (`celltype.R:445-472`) | fraction matrix, group label vector | Per cell type: `stats::wilcox.test()` (2 groups, rank-biserial effect size) or `stats::kruskal.test()` (>2 groups, epsilon-squared), then `stats::p.adjust(method="BH")` across cell types | `list(ok, table, test_used, levels)` | `cmp_result()` | **Correct, and deliberately not shared code.** The header comment states this reimplements the same auto-selected-by-group-count logic as `mod_deconvolution.R`'s `compute_group_stats()` locally rather than importing it, per this module's stated isolation policy (that transcriptomics module is explicitly off-limits to touch or depend on). **Maintainability note (not a functional bug):** this is a legitimate duplication risk — a future statistical-methodology fix (e.g. a different effect-size formula) applied to one copy would not propagate to the other. Given the explicit, documented isolation requirement, this is judged an acceptable, deliberate trade-off rather than an oversight. |

### A.9 Cell-composition ordination

| Function | Input | Purpose | Output | Downstream use | Audit |
|---|---|---|---|---|---|
| `methyl_ct_composition_pca(fractions, n_pcs)` (`celltype.R:484-493`) | fraction matrix (samples × cell types) | `stats::prcomp()` directly on the fraction matrix (scaled), after dropping zero-variance columns; requires ≥3 samples, ≥2 varying cell types | `list(ok, scores, var_explained)` | `build_ord_plot()` | **Correct, deliberately not `qc.R`'s `methyl_pca_scores()`.** Verified directly: `qc.R:569-580`'s `methyl_pca_scores()` enforces `nrow(m) < 10` (rows = probes in that CpG-scale context, after `na.omit`) before transposing to samples × top-variance-probes. Applied naively to a fraction matrix (rows = samples, not probes), that check would wrongly require ≥10 *samples* — a real risk for small cohorts or references with only a handful of cell types where the meaningful check should instead be "≥2 varying cell-type columns," which is exactly what this function checks instead. The reasoning is sound and independently verified against both implementations. |
| `methyl_ct_composition_mds(fractions, k)` (`celltype.R:495-504`) | fraction matrix, dimensions | `stats::dist()` + `stats::cmdscale()` directly on the fraction matrix; requires ≥3 samples | `list(ok, scores)` | `build_ord_plot()` | **Correct**, same reasoning as above vs. `qc.R`'s `methyl_mds_scores()` (which enforces `nrow(m) < 10 || ncol(m) < 4` on CpG-scale data, `qc.R:617-630`). |

### A.10 ggplot builders (12 functions, `celltype.R:510-679`)

All reuse `theme_arthomix()`/`ARTHOMIX_COLORS`/`ARTHOMIX_STATUS`/`arthomix_pair()` from `global.R`, verified present and behaving as described (`global.R:1417-1455`).

| Function | Chart | Audit |
|---|---|---|
| `methyl_ct_plot_marker_bar()` | CpGs retained per cell type, bar | Correct. |
| `methyl_ct_plot_marker_heatmap()` | marker × cell-type centroid heatmap, capped at 200 rows | Correct; the row cap is a rendering-performance choice, explicitly bounded rather than silently truncating without limit. |
| `methyl_ct_plot_marker_scatter()` | effect vs. specificity, colored by cell type via `arthomix_pair()` | Correct, **but see `arthomix_pair()`'s 7-color cap** below — a reference with more than 7 cell types (the registry's `blood12`/`blood12_450k`/`blood19` entries) will have some categories mapped to `NA` colors. |
| `methyl_ct_plot_stacked_bar()` | per-sample stacked composition | Correct; same `arthomix_pair()` cap caveat. |
| `methyl_ct_plot_heatmap()` | cell type × sample heatmap, optional row-normalize + `stats::hclust()` clustering | Correct; clustering only attempted when `nrow(m) > 2`/`ncol(m) > 2`, guarding `hclust()`'s minimum input size. |
| `methyl_ct_plot_box()` | box/violin of fractions per cell type, optional group split | Correct; same `arthomix_pair()` cap caveat when grouped. |
| `methyl_ct_plot_scores()` | PCA/MDS scatter, optional color-by and `ggrepel` labels | Correct. |
| `methyl_ct_plot_corr()` | cell-type correlation heatmap with printed `r` values, diverging scale | Correct. |
| `methyl_ct_plot_reconstruction()` | observed-vs-reconstructed scatter, capped at 20,000 sampled points | Correct; explicit, bounded downsampling for render performance. |
| `methyl_ct_plot_method_scatter()` | one method pair, one cell type, scatter | Correct; **declared but never called anywhere in `mod_methyl_celltype.R`** — see A.11 below. |
| `methyl_ct_plot_bland_altman()` | Bland-Altman across all cell types for a method pair, ±1.96 SD lines | Correct; standard Bland-Altman construction verified (`mean_diff`, `mean_diff ± 1.96*sd_diff`). Same `arthomix_pair()` cap caveat. |
| `methyl_ct_plot_group_diff()` | boxplot with BH-FDR significance stars per cell type | Correct; star thresholds (`***`/`**`/`*`/`ns` at 0.001/0.01/0.05) are the conventional cutoffs. Same `arthomix_pair()` cap caveat. |

### A.11 Cross-cutting note: `arthomix_pair()`'s 7-hue palette vs. this module's own up-to-19-cell-type references

`global.R:1430-1436` defines `arthomix_pair()` with exactly 7 base hues and states in its own comment that "a factor with more levels than this just runs out of distinct colors rather than cycling ggplot's rainbow" — confirmed by reading the implementation (`pal[seq_along(levels)]` returns `NA` for any index beyond 7). This module's own reference registry includes 3 entries with more than 7 cell types (`blood12`/`blood12_450k` = 12, `blood19` = 19, verified against the installed `EpiDISH` objects in A.1). Every composition figure that colors by cell type (marker scatter, stacked bar, box/violin-by-group, Bland-Altman, group-diff) will therefore show undefined (`NA`) fill/color for cell types beyond the 7th when one of these larger references is selected. This is a shared, pre-existing `global.R` limitation rather than a bug introduced by this module, but it is exercised more directly here than by most other modules because this module is one of the few that lets a user select a >7-category reference. See `07_code_audit_findings.md`.

### A.12 `methyl_ct_plot_method_scatter()` — declared but unused

Confirmed via full-file grep of `mod_methyl_celltype.R`: `methyl_ct_plot_method_scatter` (`celltype.R:630-637`, a per-cell-type single-method-pair scatter) is never referenced anywhere in the server code — `ct_cmpmethods_ba_plot` uses `methyl_ct_plot_bland_altman()` instead, which already covers all cell types for a chosen method pair in one figure. **Audit: Unused.** Not a bug (the Bland-Altman figure is arguably the more informative choice for this comparison), but a genuine dead function worth knowing about if this file is ever trimmed.

---

## Part B — `mod_methyl_celltype.R` server constructs, by tab

This module defines 3 local helper functions (`register_has_run_gate_local`, `make_plot_dl`, `plotly_safe`, `mod_methyl_celltype.R:292-316`) reused across every tab, plus roughly 70 reactive expressions/observers/render functions/download handlers. They are inventoried by tab below rather than as one flat list (see `02_tabs_and_ui.md` for the UI-input side of the same tabs); this section focuses on what each reactive computes and whether it does so correctly.

### Tab 1 (Data & QC)
- `own_raw`/`own_ready` (`reactiveVal`s) hold the module-scoped upload's raw and scale-resolved state. **Audit: Correct**, and specifically verified *not* to reproduce the two historical bugs the file's own comments describe as already fixed (a `renderUI` that both reads and writes the same `reactiveVal`, and a `paste0()` zero-length-vector recycling bug) — see `07_code_audit_findings.md` for confirmation these are resolved-historical, not live.
- `ct_source()` (`mod_methyl_celltype.R:400-412`) branches shared-vs-own and applies the M→beta conversion for the shared path. **Correct.**
- `ct_filtered()` (`mod_methyl_celltype.R:439-476`) chains `methyl_filter_missing()` → `methyl_filter_sex_chr()` → optional `methyl_filter_snp()` → optional `methyl_filter_cross_reactive()` → optional `methyl_filter_detection_p()` → `methyl_fs_sample_missing_ok()`, each filter applied against the progressively shrinking matrix, then calls `methyl_probe_retention_cascade(nrow(src$mat), cascade)`. **Audit: Bug** in the cascade-table construction specifically (not in the returned working matrix `mat`, which is built correctly step-by-step) — `methyl_probe_retention_cascade()` expects every `filters[[nm]]$keep` vector to be the same length as `n_probes_start`, but here each successive filter's `keep` vector has already shrunk to the size of the matrix at that step. See `07_code_audit_findings.md` for a reproduced, quantified example of the resulting wrong retention counts.
- `output$ct_qc_cascade_table` has `outputOptions(..., suspendWhenHidden = FALSE)` — **Correct/Reasonable**: keeps the table computed even while its tab isn't the active one, consistent with `tabsetPanel` rendering every panel body up front.

### Tab 3 (Reference & Method)
- `active_reference_full()`/`active_reference()` (`mod_methyl_celltype.R:506-525`) resolve the chosen reference and apply the cell-type-inclusion checkbox subset (registry path only), with a `validate(need(length(keep) >= 2, ...))` guard. **Correct.**
- `observeEvent(input$ct_method, ...)` (`mod_methyl_celltype.R:544-551`) repopulates the hepidish IC-column dropdown whenever the method changes to `hepidish`. **Correct**, but note it does *not* also fire when the *reference* changes while `ct_method` is already `hepidish` — if a user is on hepidish, then changes `ct_ref_id`/`ct_custom_ref_file` to a reference with different column names, `ct_hepidish_ic_col`'s choices are not automatically refreshed (the `observeEvent` only depends on `input$ct_method`). **Audit: Potential concern** — see `07_code_audit_findings.md`.

### Tab 2 (CpG Feature Selection)
- `fs_result()` (`mod_methyl_celltype.R:576-606`) implements the method dispatch (reference / variance / custom / disabled-dmc) described in `02_tabs_and_ui.md`. **Correct.**
- `ct_active_markers()`/`ct_active_marker_df()` (`mod_methyl_celltype.R:637-650`) fall back to the full active reference's CpGs when feature selection hasn't been run. **Correct**, and this fallback is exactly what makes Tab 4 usable without ever visiting Tab 2.
- `ct_overlap()`/`ct_overlap_ok()` (`mod_methyl_celltype.R:652-665`) compute the overlap gate. **Correct** as far as it goes, but see A.4 above: this is an ID-presence check, not a per-sample-completeness check.
- `observe({ ok <- ...; shinyjs::enable/disable("ct_run_decon_btn") })` (`mod_methyl_celltype.R:689-692`) is the client-side gate. **Audit: Correctly backed up server-side** — `decon_result()` independently re-derives `ref_use` from `active_reference()`/`ct_active_markers()` and re-validates `nrow(ref_use) >= 10` (`mod_methyl_celltype.R:702-709`), and `methyl_ct_run_epidish()`/`methyl_ct_run_hepidish()` each *also* independently re-check the CpG overlap (`length(common) < 10`) before running. This is a triple-redundant guard: even if the client-side `shinyjs::disable()` were bypassed (e.g. via browser devtools), clicking the button would still hit `eventReactive`'s own `validate()` and the wrapper functions' own checks, all of which fail closed with an explanatory message rather than running on too little data. **This is a positive finding**, not a gap.

### Tab 4 (Deconvolution)
- `decon_result()` (`mod_methyl_celltype.R:702-728`) is the estimation call itself, `set.seed()`'d immediately before it (see A.11/`07` for why this seed is currently inert for every available method). **Correct**, modulo the seed-inertness note.
- `observeEvent(decon_result(), ...)` (`mod_methyl_celltype.R:730-739`) writes `results$celltype <- list(method, cell_types, n_samples, n_markers_used, mean_fraction)`. **Correct** as a summary; see `04_data_flow_and_pipeline.md` for its one confirmed downstream reader (ArthOChat's context builder).

### Tab 5 (Cell Composition)
- `ct_group_vec(col)` (`mod_methyl_celltype.R:892-897`) resolves a sample-sheet column to a named vector aligned to `rownames(r$fractions)` via `methyl_sheet_sample_ids()`. **Correct.**
- `build_ord_plot()` (`mod_methyl_celltype.R:908-921`) is a plain function (not itself a reactive) called from both `renderPlot` and the PNG download handler, so the PCA/MDS figure and its download are guaranteed pixel-identical. **Correct/Reasonable pattern.**
- `ct_sample_order_vec()` (`mod_methyl_celltype.R:860-867`) supports "As in data"/"Alphabetical"/"By dominant cell type" ordering for the stacked bar. **Correct.**

### Tab 6 (Validation)
- `val_result()` (`mod_methyl_celltype.R:974-981`) and `cmpm_result()` (`mod_methyl_celltype.R:1013-1023`) both operate on `decon_result()$working_mat`/`ref_used` — i.e. the exact matrix/reference frozen at the moment Tab 4 was last run, not a live re-derivation from Tabs 1–3's current inputs. **Reasonable and correct as implemented**, but worth being explicit about for users: changing a QC filter after running Deconvolution does not retroactively change what Validation/Compare Methods show until Deconvolution is re-run (documented in `02_tabs_and_ui.md`'s "Connection to other tabs").

### Tab 7 (Export)
- All 7 `downloadHandler`s and the `observe({ shinyjs::toggleState(...) })` sync block (`mod_methyl_celltype.R:1062-1125`) are straightforward serializations of already-validated reactive results. **Correct.** The summary-report handler (`ct_export_report`, `mod_methyl_celltype.R:1097-1117`) builds its text conditionally per has-run flag, so a partial session (e.g. QC done, deconvolution not run) still produces a coherent, non-crashing report rather than an error.

---

## Part C — External package functions actually called

| Function | Package | Called from | Purpose | Verified behavior in this environment |
|---|---|---|---|---|
| `EpiDISH::epidish()` | EpiDISH (installed) | `methyl_ct_run_epidish()` | Single-stage reference-based deconvolution (CP/RPC/CBS) | Dispatches to internal `DoCP`/`DoRPC`/`DoCBS`; all three verified deterministic (no `runif`/`rnorm`/`sample` calls in their source) |
| `EpiDISH::hepidish()` | EpiDISH (installed) | `methyl_ct_run_hepidish()` | Two-stage hierarchical deconvolution | Verified: runs the chosen method on `ref1`/`ref2` independently, then substitutes `frac1[,ic] * frac2` for the IC column |
| `stats::wilcox.test()` | base R `stats` | `methyl_ct_group_stats()` | 2-group Wilcoxon rank-sum test | Standard library function |
| `stats::kruskal.test()` | base R `stats` | `methyl_ct_group_stats()` | >2-group Kruskal-Wallis test | Standard library function |
| `stats::p.adjust(method="BH")` | base R `stats` | `methyl_ct_group_stats()` | Benjamini-Hochberg FDR across cell types | Standard library function |
| `stats::prcomp()` | base R `stats` | `methyl_ct_composition_pca()` | PCA on the fraction matrix | Standard library function |
| `stats::cmdscale()` | base R `stats` | `methyl_ct_composition_mds()` | Classical MDS on a distance matrix | Standard library function |
| `stats::hclust()` | base R `stats` | `methyl_ct_plot_heatmap()` | Hierarchical clustering for heatmap row/column ordering | Standard library function |
| `stats::cor()` | base R `stats` | `methyl_ct_method_correlation()`, `ct_corr_plot`, `methyl_ct_validation_metrics()` | Pearson correlation | Standard library function |
| `stats::complete.cases()` | base R `stats` | `methyl_ct_run_epidish()` | Listwise-deletion mask for rows with any `NA` | Standard library function; see `07` for the downstream implication of this being applied *after* the ID-overlap QC check |
| `stats::sd()`, `stats::var()`, `stats::median()`, `stats::dist()`, `stats::quantile()` | base R `stats` | multiple | Standard summary/distance statistics | Standard library functions |
| `MASS::rlm()` (transitively, via `EpiDISH::DoRPC`) | MASS | not called directly by this codebase | Robust linear regression underlying RPC | Verified present in the installed `EpiDISH::DoRPC` source |
| `e1071::svm()` (transitively, via `EpiDISH::DoCBS`) | e1071 | not called directly by this codebase | nu-SVR underlying CBS | Verified present in the installed `EpiDISH::DoCBS` source |
| `quadprog::solve.QP()` (transitively, via `EpiDISH::DoCP`) | quadprog | not called directly by this codebase | Quadratic-programming solver underlying CP | Verified present in the installed `EpiDISH::DoCP` source |
