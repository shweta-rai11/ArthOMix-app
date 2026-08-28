# 07. Code Audit Findings — Methylomics Cell-Type Deconvolution

Scope: `mod_methyl_celltype.R` and `celltype.R`, plus their direct verified interactions with `qc.R`, `global.R`, `submodules_registry.R`, `server.R`, and the installed `EpiDISH` package. Two historical bugs are documented in the source's own comments as already fixed (a `renderUI` that both read and wrote the same `reactiveVal`, around `ct_scale_ui`/`own_ready`; and a `paste0()` zero-length-vector recycling bug, around `ct_ord_color`'s choices). Both were re-checked against the current code and confirmed resolved as described — they are **not** re-reported below as live findings, and this file's own review found no other unguarded instance of either pattern anywhere else in `mod_methyl_celltype.R`.

---

## Finding 1 — Probe-retention cascade table shows incorrect counts after the first QC filter

**Severity: Moderate**

**Evidence.** `ct_filtered()` (`mod_methyl_celltype.R:439-476`) applies QC filters sequentially, each against the matrix already narrowed by the previous filter:

```r
f1 <- methyl_filter_missing(mat, max_na_frac = miss_cpg); cascade[["Missing-CpG filter"]] <- f1
mat <- mat[f1$keep, , drop = FALSE]

f2 <- methyl_filter_sex_chr(mat, anno, mode = input$ct_qc_sexchr %||% "remove_xy"); cascade[["Chromosome scope"]] <- f2
mat <- mat[f2$keep, , drop = FALSE]
...
cascade_df <- methyl_probe_retention_cascade(nrow(src$mat), cascade)
```

`f2$keep` (and every subsequent filter's `$keep`) is therefore a logical vector whose length equals `nrow(mat)` *after* the previous filter, not `nrow(src$mat)` (the original count passed to `methyl_probe_retention_cascade()`). That function's contract (verified against its correct usage in `mod_methyl_qc.R:607-628`, where every filter is instead computed against the *same*, unfiltered matrix before being combined) requires every `filters[[nm]]$keep` to have length `n_probes_start`:

```r
methyl_probe_retention_cascade <- function(n_probes_start, filters) {
  keep <- rep(TRUE, n_probes_start)
  for (nm in names(filters)) {
    keep <- keep & filters[[nm]]$keep   # length mismatch from the 2nd filter onward
    ...
  }
}
```

Reproduced directly in this R environment with a minimal example (1,000 starting probes, a first filter keeping 890, a second filter's `keep` vector correctly sized to 890 as `ct_filtered()` would produce it): R's `&` recycles the shorter vector against the longer one with a "longer object length is not a multiple of shorter object length" warning, and the cascade table reports **798** probes retained after the second filter when the true count (verified independently) is **799** — an off-by-a-scrambled-amount error, not a consistent or predictable offset, because the misalignment is a recycling artifact, not a simple truncation.

**Impact.** The step-by-step "retained"/"removed" counts shown in `ct_qc_cascade_table` (Tab 1) are unreliable for every step after the first filter. This does **not** affect the actual working matrix used for deconvolution (`ct_filtered()$mat` is built correctly, filter-by-filter, independent of the cascade table), nor the top-line "Working matrix after filters: X CpGs × Y samples" message (computed directly from `nrow(f$mat)`, not from the cascade table). The practical effect is limited to a misleading intermediate QC-transparency table and R warnings appearing in server logs on every recompute.

**Recommended action (describe only, not implemented here).** Either compute every filter's `keep` mask against the same unfiltered matrix and combine them at the end (matching `mod_methyl_qc.R`'s pattern), or change the cascade-building call to pass each filter's own current matrix size instead of the original `nrow(src$mat)`.

---

## Finding 2 — The "Constraint" advanced parameter silently has no effect for RPC/CBS

**Severity: Moderate**

**Evidence.** `ct_adv_constraint` (Tab 3, Advanced Parameters) is passed as `constraint` to `EpiDISH::epidish()`/`hepidish()` regardless of which `ct_method` is selected (`mod_methyl_celltype.R:719-724`). Reading the installed `EpiDISH::epidish()` source directly:

```r
function (beta.m, ref.m, method = c("RPC", "CBS", "CP"), maxit = 50, nu.v = c(0.25, 0.5, 0.75), constraint = c("inequality", "equality")) {
    method <- match.arg(method); constraint <- match.arg(constraint)
    if (method == "RPC") { out.o <- DoRPC(beta.m, ref.m, maxit) }
    else if (method == "CBS") { out.o <- DoCBS(beta.m, ref.m, nu.v) }
    else if (method == "CP") { ... out.o <- DoCP(beta.m, ref.m, constraint) }
    ...
}
```

`constraint` is validated at the top (so an invalid string would still error) but is only ever forwarded into the actual computation inside the `method == "CP"` branch. `DoRPC()`/`DoCBS()` never receive it — both always clip negative coefficients to zero and renormalize to sum to 1 internally, independent of the `constraint` selector.

**Impact.** A user who selects RPC or CBS (2 of the 4 available methods, and 2 of the 3 methods hepidish can internally use — though hepidish is itself hard-coded to RPC, see Finding 4) and then changes "Constraint" from Inequality to Equality on the Advanced Parameters panel will see **no change in the result**, with no indication in the UI that this parameter is CP-only. This is a UI-transparency gap rather than a computational error (the underlying package behavior is exactly as documented for CP; this codebase's wrapper does not misrepresent it in code, only in the UI's undifferentiated presentation of the control).

**Recommended action.** Either hide/gray out the Constraint selector when `ct_method` is not `CP`, or add an inline note clarifying it only applies to CP.

---

## Finding 3 — "Random seed" advanced parameter is currently inert for every available method

**Severity: Low / Informational**

**Evidence.** `decon_result()` calls `set.seed(input$ct_adv_seed %||% 1234)` immediately before every estimation call (`mod_methyl_celltype.R:711-712`). Directly reading the installed `EpiDISH` package's internal `DoCP()` (`quadprog::solve.QP()`), `DoRPC()` (`MASS::rlm()`), and `DoCBS()` (`e1071::svm()`, linear-kernel nu-regression) confirms none of the three call any R random-number function — all are deterministic optimizations for a fixed input matrix.

**Impact.** Changing the seed value has no observable effect on results for CP, RPC, CBS, or hepidish (which always uses RPC internally — Finding 4). This is not a correctness bug (the seed doesn't need to do anything for a deterministic method, and setting it anyway is harmless), but the UI presents it as a meaningful "Advanced Parameter" without noting that none of the currently implemented methods are stochastic, which could give a false impression that re-running with a different seed is a way to check estimate stability.

**Recommended action.** Add a note that current methods are deterministic and the seed is reserved for forward compatibility (e.g. if a future reference-free or bootstrap-based method were added).

---

## Finding 4 — hepidish's internal EpiDISH method is hard-coded to RPC, not exposed or documented in the UI

**Severity: Low / Informational**

**Evidence.** `decon_result()`'s hepidish branch:

```r
res <- methyl_ct_run_hepidish(f$mat, ref_use, ref2, ic_column = input$ct_hepidish_ic_col,
                               method = "RPC", maxit = input$ct_adv_maxit %||% 50, nu.v = nu,
                               constraint = input$ct_adv_constraint %||% "inequality")
```

`method = "RPC"` is a literal, not read from any input. `methyl_ct_run_hepidish()`'s own signature does default to `method = c("RPC", "CBS", "CP")` (RPC first), so this call site's explicit `"RPC"` matches the function's own default, but there is no UI affordance anywhere on Tab 3 to choose CP or CBS specifically for the two-stage path, and no text stating that hepidish always uses RPC internally.

**Impact.** A user selecting "Two-stage (hepidish) — advanced" gets an RPC-based two-stage estimate with no way to request CP or CBS instead, and no way to discover this without reading the source. Not a bug — RPC is a reasonable default — but a missing piece of UI transparency.

**Recommended action.** Either surface a method selector specific to hepidish, or add explanatory text stating the internal method used.

---

## Finding 5 — Overlap-QC gate checks CpG-ID presence, not per-sample completeness; the final marker count used can be silently lower than what the gate reported

**Severity: Moderate**

**Evidence.** `methyl_ct_overlap_qc()` (`celltype.R:266-273`) is a pure set-intersection between the active marker IDs and `rownames(ct_filtered()$mat)` — it does not check whether each matched CpG has a non-missing value in every sample. `methyl_ct_run_epidish()` (`celltype.R:299-319`), however, applies `stats::complete.cases()` to the marker-restricted matrix *after* the ID intersection:

```r
common <- intersect(rownames(beta_mat), rownames(ref_mat))
...
bm <- beta_mat[common, , drop = FALSE]
...
complete <- stats::complete.cases(bm)
bm <- bm[complete, , drop = FALSE]
```

Any marker CpG with even one remaining `NA` in any one sample is dropped from **every** sample's estimation at this step, even though it was already counted as "matched" by the overlap-QC gate the Run button's enabled state depends on. Because Tab 1's "Missing-CpG threshold" filter (default 5%) deliberately tolerates some residual missingness per CpG rather than requiring zero, it is expected that some retained CpGs still carry scattered `NA`s, which this later listwise-deletion step then removes.

**Impact.** The percentage/count shown on Tab 3 ("X% of reference marker CpGs matched") can overstate the number of CpGs actually contributing to the final estimate. The *true* final count (`n_markers_used`) is correctly computed and surfaced downstream (Tab 4's result, the Export report, `results$celltype`), so this is not a silently-lost number overall — but there is no cross-reference back to the overlap-QC percentage explaining the discrepancy if one occurs, and no imputation option is offered before deconvolution (unlike `mod_methyl_featureselection.R`'s `methyl_fs_impute()`, which exists elsewhere in this app's Methylomics vertical but is not reused here).

**Recommended action.** Either have `ct_overlap()` also report a completeness-aware "usable" count alongside the ID-presence count, or add a one-line note near the overlap-QC message clarifying that per-sample missingness among matched CpGs is resolved by listwise deletion at run time and may reduce the effective marker count further.

---

## Finding 6 — `arthomix_pair()`'s 7-hue palette is exceeded by 3 of this module's own reference panels

**Severity: Low**

**Evidence.** `global.R:1430-1436`'s `arthomix_pair()` provides exactly 7 base hues and its own comment states levels beyond 7 simply run out of color rather than cycling. This module's reference registry (verified against the installed `EpiDISH` objects) includes `blood12`/`blood12_450k` (12 cell types) and `blood19` (19 cell types) — all three selectable from `ct_ref_id`. Every cell-type-colored figure in this module (`methyl_ct_plot_marker_scatter()`, `_stacked_bar()`, `_box()` when grouped, `_bland_altman()`, `_group_diff()`) calls `arthomix_pair()` on the active cell-type set.

**Impact.** Selecting a >7-cell-type reference will produce `NA`-valued color/fill entries for the 8th and later cell types in these figures — those categories render with ggplot's default missing-value handling (typically grey/blank) rather than a distinct color, degrading (but not crashing) the visualization. This is a shared, pre-existing limitation of `global.R`'s palette rather than a defect introduced by this module, but this module is one of only a few places in the app that lets a user select a >7-category grouping.

**Recommended action.** Extend `arthomix_pair()`'s base palette (a `global.R` change, out of scope for this module's own files) or have this module's larger reference panels use a different, higher-cardinality palette locally.

---

## Finding 7 — `ct_qc_beta_range` checkbox is declared but never read

**Severity: Low**

**Evidence.** `checkboxInput(ns("ct_qc_beta_range"), "Validate beta-value range (0-1)", value = TRUE)` (`mod_methyl_celltype.R:76`) is the only occurrence of `ct_qc_beta_range` anywhere in the file — confirmed by an exhaustive grep. No `input$ct_qc_beta_range` read exists in `ct_filtered()` or anywhere else.

**Impact.** The checkbox is fully inert — toggling it changes nothing about the QC pipeline's behavior. A user reading the UI would reasonably expect unchecking it to disable some beta-range validation step, but no such step is gated on this input anywhere.

**Recommended action.** Either wire it to an actual range-validation filter (e.g. flagging/excluding out-of-range values before the missingness/sex-chromosome filters), or remove the control if no such validation is intended.

---

## Finding 8 — hepidish's IC-column dropdown does not refresh when the reference changes while hepidish is already selected

**Severity: Low**

**Evidence.** `observeEvent(input$ct_method, { req(identical(input$ct_method, "hepidish")); ref <- tryCatch(active_reference_full(), error = function(e) NULL); if (!is.null(ref)) updateSelectInput(session, "ct_hepidish_ic_col", choices = colnames(ref), ...) })` (`mod_methyl_celltype.R:544-551`) is triggered only by `input$ct_method` changing — `active_reference_full()` is read inside the handler body but, since `observeEvent`'s handler expression executes isolated from reactivity by design, this does not add a reactive dependency on the reference itself.

**Impact.** If a user is already on `ct_method == "hepidish"` and then changes the reference selection (`ct_ref_id`, cell-type checkboxes, or uploads a different custom reference), `ct_hepidish_ic_col`'s dropdown keeps showing the *previous* reference's column names. If the previously chosen column name happens not to exist in the new reference, `methyl_ct_run_hepidish()`'s own guard (`match(ic_column, colnames(ref1_mat))` returning `NA`) will correctly fail closed with an explanatory message rather than silently misinterpreting the wrong column — so this is a UI staleness issue, not a silent-data-corruption issue, but it is a genuine, reproducible confusion point.

**Recommended action.** Also trigger the IC-column dropdown refresh reactively on `active_reference_full()` (e.g. via a plain `observe()` guarded on `ct_method == "hepidish"` rather than an `eventReactive`/`observeEvent` keyed only on the method input).

---

## Finding 9 — Server-side re-validation of the deconvolution gate is robust (positive finding)

**Severity: Informational (not a defect)**

**Evidence.** The client-side Run-button gate (`observe({ ok <- tryCatch(isTRUE(ct_overlap_ok()), error=function(e) FALSE); if (ok) shinyjs::enable(...) else shinyjs::disable(...) })`, `mod_methyl_celltype.R:689-692`) only controls whether the button is clickable in the browser. Independently, `decon_result()` re-derives `ref_use` from `active_reference()`/`ct_active_markers()` from scratch and re-validates `nrow(ref_use) >= 10` (`mod_methyl_celltype.R:702-709`), and `methyl_ct_run_epidish()`/`methyl_ct_run_hepidish()` each independently re-check the actual CpG overlap against the working matrix (`length(common) < 10`) before calling EpiDISH. This is three independent layers, two of them entirely server-side and unrelated to the client-side `shinyjs` state.

**Impact.** Even if the button's disabled attribute were bypassed client-side (e.g. via browser devtools) and the underlying Shiny input event fired anyway, the server-side `eventReactive` and wrapper functions would still block execution on too little data with a clear `validate()`/`reason` message rather than running EpiDISH on a degenerate input. This was explicitly checked per the audit's scope and found to be correctly implemented.

---

## Categories checked with no significant finding

- **Incorrect matrix orientation:** every filter/statistic function consistently treats rows as CpGs and columns as samples for the working/reference matrices, and rows as samples/columns as cell types for the fraction matrix; verified consistent throughout `celltype.R` and the reconstruct/validate math. No orientation bug found.
- **Wrong sample/feature matching:** `methyl_sheet_sample_ids()` is consistently used to align sample-sheet rows to fraction-matrix row names wherever a phenotype column is read (`ct_group_vec()`, the phenotype-linked export). No mismatch found.
- **Small-sample-size edge cases:** `methyl_ct_composition_pca()`/`_mds()` correctly guard `<3` samples / `<2` varying cell types; `methyl_ct_run_epidish()`/`_hepidish()` correctly guard `<10` overlapping/complete-case marker CpGs; `methyl_ct_validation_metrics()` guards `<2` CpGs/`<1` sample. All guards fail closed with an explanatory `validate()`/`reason` message rather than crashing or silently returning nonsense.
- **Mismatched-reference edge cases:** selecting a reference with too few retained cell-type columns (`<2`) is explicitly guarded via `validate(need(length(keep) >= 2, ...))` in `active_reference()`; an `ic_column` not present in the second-stage reference is explicitly guarded in `methyl_ct_run_hepidish()`.
- **Error-handling weaknesses:** every user-facing computation path (`own_matrix_parsed()`, `ct_custom_ref_raw()`, `decon_result()`, `val_result()`, `cmpm_result()`, `cmp_result()`) is wrapped in `tryCatch()`/`validate()` patterns that degrade to an explanatory message rather than an uncaught error. No unguarded computation reachable from a normal UI interaction was found to throw an unhandled error.
