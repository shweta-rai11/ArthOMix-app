# 04. Data Flow and Pipeline — Methylomics Cell-Type Deconvolution

All dimensions below are runtime-determined from the loaded dataset/reference; no sample count, CpG count, or cell-type count is fixed by this module's code except where a specific registry panel's dimensions are stated (verified against the installed `EpiDISH` package, itself a fixed property of the deployment, not of any particular user's uploaded data).

## 4.1 End-to-end object flow

```
                              ┌───────────────────────────────────────────┐
                              │  methyl_dataset$beta  (shared, Dataset tab) │
                              │  OR module-scoped upload (own_raw)          │
                              │  CpG rows × sample columns                  │
                              └───────────────────┬─────────────────────────┘
                                                   │  ct_source()  (mod_methyl_celltype.R:400-412)
                                                   │  - shared path: convert M→beta if input_scale=="m"
                                                   │  - own path: gated on own_ready() (post-transform)
                                                   ▼
                              ┌───────────────────────────────────────────┐
                              │  src$mat  — beta-scale working matrix       │
                              └───────────────────┬─────────────────────────┘
                                                   │  ct_filtered()  (mod_methyl_celltype.R:439-476)
                                                   │  methyl_filter_missing() → methyl_filter_sex_chr()
                                                   │  → [methyl_filter_snp()] → [methyl_filter_cross_reactive()]
                                                   │  → [methyl_filter_detection_p()] → methyl_fs_sample_missing_ok()
                                                   ▼
                              ┌───────────────────────────────────────────┐
                              │  ct_filtered()$mat — QC-filtered matrix     │
                              │  (fewer CpG rows, possibly fewer samples)   │
                              └───────┬───────────────────────────┬─────────┘
                                      │                            │
        Reference path               │                            │  consumed directly by decon_result()
        (independent of the above)   │                            │  and by fs_result() (for chr-scope lookup)
                                      │                            │
   ┌──────────────────────────┐      │                            │
   │ registry (7 EpiDISH       │      │                            │
   │ built-in panels) OR       │      │                            │
   │ custom uploaded reference │      │                            │
   │ CpG rows × cell-type cols │      │                            │
   └─────────────┬──────────────┘    │                            │
                  │ active_reference_full() → active_reference()   │
                  │ (mod_methyl_celltype.R:506-525;                │
                  │  optional cell-type-column subset)              │
                  ▼                                                 │
   ┌──────────────────────────┐                                     │
   │ ref  — active reference    │                                    │
   └─────────────┬──────────────┘                                    │
                  │ methyl_ct_marker_rank(ref)  (celltype.R:160-193) │
                  ▼                                                  │
   ┌──────────────────────────┐                                      │
   │ rank_df — per-CpG effect/  │                                     │
   │ specificity/direction table │                                    │
   └─────────────┬──────────────┘                                     │
                  │ methyl_ct_select_markers() + methyl_ct_top_n_balanced()
                  │ (fs_result(), mod_methyl_celltype.R:576-606; OPTIONAL —
                  │  skipped entirely if Tab 2 is never run)
                  ▼
   ┌──────────────────────────┐
   │ ct_active_markers()        │  = fs_result()$selected$cpg  if Tab 2 run,
   │ (CpG ID vector)             │  else rownames(active_reference())
   └─────────────┬──────────────┘
                  │
                  ├────────────────────────────────────────────────────────┐
                  ▼                                                        ▼
   ┌──────────────────────────┐                          ┌──────────────────────────────────┐
   │ ct_overlap()/ct_overlap_ok()│  ID-overlap QC gate     │  decon_result()  (mod_methyl_celltype.R:702-728) │
   │ vs. ct_filtered()$mat rows  │  (enables/disables the  │  ref_use = ref[intersect(rownames(ref), markers),] │
   │ (celltype.R:266-285)        │  Run Deconvolution btn) │  → methyl_ct_run_epidish() / methyl_ct_run_hepidish()│
   └──────────────────────────┘                          └──────────────────┬────────────────┘
                                                                             │
                                                     inside the wrapper: intersect(rownames(beta_mat), rownames(ref_mat))
                                                     → stats::complete.cases() listwise deletion (celltype.R:309-311)
                                                     → EpiDISH::epidish()/hepidish()
                                                                             ▼
                                              ┌──────────────────────────────────────────┐
                                              │  fractions  — samples × cell types matrix   │
                                              │  (non-negative, ~sum-to-one per row)         │
                                              └──────┬───────────┬───────────┬──────────────┘
                                                     │           │           │
                    results$celltype  ◄──────────────┘           │           │
                    (mod_methyl_celltype.R:730-739;               │           │
                     summary only: method, cell_types,            │           │
                     n_samples, n_markers_used, mean_fraction)     │           │
                                                                   │           │
      ┌────────────────────────────────────────────────────────────┘           │
      ▼                                                                        │
┌─────────────────────────────┐   ┌────────────────────────────┐              │
│ Cell Composition (Tab 5):     │   │ Validation (Tab 6):          │              │
│ stacked bar / heatmap / box /  │   │ methyl_ct_reconstruct(         │◄─────────────┘
│ PCA-MDS / correlation matrix /  │  │   ref_used, fractions)          │  (uses decon_result()$ref_used,
│ methyl_ct_group_stats() vs.     │   │ methyl_ct_validation_metrics(   │   $working_mat — frozen at the
│ uploaded phenotype column        │  │   working_mat, reconstructed)   │   moment Tab 4 last ran)
└─────────────────────────────┘   │ methyl_ct_compare_methods()      │
                                   │ (reruns epidish per method on    │
                                   │  the SAME working_mat/ref_used)  │
                                   └────────────────────────────┘
                                                    │
                                                    ▼
                                   ┌────────────────────────────┐
                                   │  Export (Tab 7): 6 CSV +       │
                                   │  1 text summary report          │
                                   └────────────────────────────┘
```

## 4.2 Stage-by-stage object identity and dimensions

| Stage | Object | Reactive/variable | Dimensions | Notes |
|---|---|---|---|---|
| Raw input | `src$mat` | `ct_source()` | rows = CpGs (data-dependent count), cols = samples (data-dependent count) | Beta-scale guaranteed by this point (own path gated on `own_ready()`; shared path converts M→beta if needed) |
| QC-filtered | `ct_filtered()$mat` | `ct_filtered()` | rows ≤ `nrow(src$mat)`, cols ≤ `ncol(src$mat)` | Exact counts data-dependent on which QC filters are enabled and their thresholds |
| Reference | `ref` / `active_reference()` | `active_reference_full()`/`active_reference()` | rows = reference's own CpG count (e.g. 333/188/600/600/716/491/1906 for the 7 registry panels respectively, verified against the installed `EpiDISH` objects — a fixed deployment property, or data-dependent for a custom upload), cols = 2..N cell types (2..7/7/12/12/3/4/19 for the registry panels, or user-chosen count for a custom upload) | Registry references are fixed by the installed EpiDISH version; a custom upload's dimensions are entirely data-dependent |
| Marker ranking | `rank_df` | inside `fs_result()`, via `methyl_ct_marker_rank(ref)` | rows = `nrow(ref)` (one row per reference CpG), cols = 8 (`cpg`,`cell_type`,`effect`,`direction`,`specificity`,`other_mean`,`btw_type_var`,`centroid_beta`) | Computed off reference centroids only — never touches the working matrix |
| Marker selection | `fs_result()$selected` | `fs_result()` (optional; skipped if Tab 2 unused) | rows ≤ `nrow(rank_df)`, capped at the user's `top_n` | Feeds `ct_active_markers()` |
| Active marker set | `ct_active_markers()` | `ct_active_markers()` | vector, length = `nrow(fs_result()$selected)` if Tab 2 run, else `nrow(active_reference())` | CpG ID vector only, not a matrix |
| Overlap QC | `ct_overlap()` | `ct_overlap()` | scalar/`data.frame` summary, not a matrix | Compares `ct_active_markers()` against `rownames(ct_filtered()$mat)` |
| Deconvolution input | `ref_use` (inside `decon_result()`) | `decon_result()` | rows = `length(intersect(rownames(ref), markers))`, cols = active reference's cell-type count | Re-derived independently of `ct_overlap()`, not read from it |
| Deconvolution input | effective marker set actually used by `epidish()` | inside `methyl_ct_run_epidish()` | rows = `length(intersect(rownames(ref_use), rownames(f$mat)))` further reduced by `stats::complete.cases()` | This is the **true** final marker count, reported back as `n_markers_used`; can be lower than `ct_overlap()`'s "matched" count because of the additional per-sample completeness requirement (see `07_code_audit_findings.md`) |
| Estimated fractions | `decon_result()$fractions` | `decon_result()` | rows = samples in `ct_filtered()$mat` (or fewer, if `complete.cases` also removed a whole sample's worth of data — it does not, since `complete.cases()` operates row-wise on CpGs, not columns), cols = active reference's cell-type count | Row sums are ~1 for all three EpiDISH methods (verified — see `03_functions_and_code_audit.md` A.5) |
| Group summary | `results$celltype` | `observeEvent(decon_result(), ...)` | small named list, not a matrix | See §4.3 |
| Reconstruction | `methyl_ct_reconstruct()`'s output | `val_result()` | rows = CpGs common to `ref_used`/`fractions`' cell types, cols = samples | Forward simulation of the mixture model |
| Validation metrics | `val_result()$overall`/`$per_sample` | `val_result()` | 4 scalars (overall) + one row per sample (per-sample) | — |
| Cross-method fractions | `cmpm_result()$fractions_by_method` | `cmpm_result()` | named list, one fraction matrix per method chosen (2 or 3) | All computed on the identical `decon_result()$working_mat`/`ref_used` |

## 4.3 The one confirmed downstream consumer of `results$celltype`

`observeEvent(decon_result(), ...)` (`mod_methyl_celltype.R:730-739`) writes:

```r
results$celltype <- list(
  method = r$method, cell_types = colnames(r$fractions), n_samples = nrow(r$fractions),
  n_markers_used = r$n_markers_used, mean_fraction = round(colMeans(r$fractions), 4)
)
```

This is written into `methyl_results` (`server.R:93`, the shared per-Methylomics-vertical results store), *not* into `methyl_dataset`. A whole-repository grep for `results$celltype`, `methyl_results$celltype`, and `\$celltype\b` (excluding the `mod_methyl_wgcna.R` helper function `mx_wgcna_celltype_reference()`, an unrelated same-word coincidence) confirms:

- **No other analysis submodule reads it.** `mod_methyl_biomarkercard.R` — the module most likely to want a cell-composition confounder — does not reference `celltype` anywhere.
- **It is read generically by ArthOChat's context builder.** `build_mx_context()` (`submodules_registry.R:163-181`) iterates every `MX_MODULES` id, including `"celltype"`, and calls `.format_results_block(MX_MODULES_BY_ID[["celltype"]]$config$title, methyl_results[["celltype"]])` (`submodules_registry.R:179`). `.format_results_block()` (`submodules_registry.R:125-139`) renders each named field as a `- name: value` bullet (values truncated to the first 20 `as.character()` elements) if `results$celltype` is non-`NULL`, or a "NOT YET RUN IN THIS SESSION" placeholder otherwise. `build_mx_context()` is invoked from the module that assembles ArthOChat's system-prompt context (`server.R:184-189`; the switch statement is defined in `submodules_registry.R` and called with `focus_id = submodule_id`). **Net effect:** the AI chat assistant can see and describe this module's estimated cell-type composition once it has been run, but no other quantitative analysis module in the app consumes it programmatically.

## 4.4 What is *not* connected

- Changing `ct_filtered()`'s QC filters or `active_reference()`'s reference/method selection **after** `decon_result()` has already run does not automatically invalidate Tabs 5/6 — they read `decon_result()$fractions`/`$working_mat`/`$ref_used`, which only update when the "Run Cell-Type Deconvolution" button is clicked again.
- No other Methylomics submodule (QC, Normalization, DMP, DMR, WGCNA, Candidates, Feature Selection, MR, Coloc, Diagnostic, Biomarker Card) reads any output of this module — confirmed by grepping each file in `R/methylomics/` for `celltype`; only `submodules_registry.R` (registry wiring) and this module's own two files reference it.
- The module does not write anything back to `methyl_dataset` itself — its outputs are entirely self-contained within `methyl_results$celltype` and the Export tab's downloads.
