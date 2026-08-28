# 02. Tabs and UI — Methylomics Cell-Type Deconvolution

**7 sub-tabs**, one `tabsetPanel(id = ns("ct_subtabs"), type = "tabs", ...)` (`mod_methyl_celltype.R:269-278`). Since a `tabsetPanel` renders every panel body up front in Shiny, all 7 tabs' UI elements exist in the DOM simultaneously; visibility is CSS-controlled by tab selection only, not by conditional creation. Reactive values (`decon_result()`, `fs_result()`, etc.) are the actual gates on whether each tab's content is meaningful, via the `register_has_run_gate_local()` pattern (`mod_methyl_celltype.R:292-297`).

---

## Tab 1 — "1. Data & QC" (`value = "dataqc"`, `mod_methyl_celltype_dataqc_ui()`, `mod_methyl_celltype.R:33-88`)

**Purpose:** select/upload the working methylation matrix, resolve its scale to beta, and apply the standard QC probe/sample filters this app reuses from `qc.R`, producing the working matrix every later tab operates on.

### Box "1. Data source"
| Input | Type | Choices | Default | Mandatory | Conditional | Server use |
|---|---|---|---|---|---|---|
| `ct_data_source` | `radioButtons` | "Use the dataset loaded on the Dataset tab" (`shared`) / "Upload a different matrix here" (`own`) | `shared` | Yes (always one selected) | — | Branches `ct_source()` (`mod_methyl_celltype.R:400-412`) |
| `ct_own_array_type` | `selectInput` | `METHYL_ARRAY_TYPES` (450K/EPIC/EPICv2/WGBS/RRBS/Custom array) | `EPIC` | Only relevant if `own` | shown when `ct_data_source == "own"` | Stored as `array_type` on the uploaded-raw list; used later for annotation lookups (`methyl_get_annotation()`) |
| `ct_own_matrix_file` | `fileInput` | .csv/.tsv/.txt | none | Yes, for `own` path | shown when `own` | Parsed by `methyl_parse_matrix()` in `own_matrix_parsed()` |
| `ct_own_sheet_file` | `fileInput` | .csv/.tsv/.txt | none (optional) | No | shown when `own` | Parsed by `methyl_parse_sample_sheet()` in `own_sheet_parsed()` |
| (preview) `ct_own_preview_ui` | `uiOutput` | — | — | — | shown when a file is selected | Shows row/column counts or a parse error |
| `ct_own_load_btn` | `actionButton` | — | — | click required to commit | shown when `own` | Fires `observeEvent`; sets `own_raw()` and (if already beta-scale) `own_ready()` |

### Box "2. Dataset preview" + scale card
- `ct_source_summary_ui` (`uiOutput`) — always visible; shows source label, CpG × sample counts, % missing, beta range, duplicated-probe count, and chromosome count (via `methyl_ct_working_summary()`, degrading to "unknown" when no manifest annotation resolves).
- `ct_scale_ui` (`uiOutput`, own-upload path only) — shows either a confirmation the matrix is on the beta scale, or a warning box with an `ct_apply_transform_btn` `actionButton` labeled dynamically `"Apply Transformation (% -> Beta)"` or `"Apply Transformation (M-value -> Beta)"`. Clicking it is the only way the working matrix's scale changes; nothing transforms automatically.

### Box "3. QC filters" (`mod_methyl_celltype.R:61-87`)
| Input | Type | Choices | Default | Server use |
|---|---|---|---|---|
| `ct_qc_missing_cpg` (+ `_custom`) | `selectInput` (+ conditional `numericInput`, min 0/max 1/step 0.01) | 0% / 1% / 5% / 10% / Custom | 5% | `methyl_filter_missing()` max-NA-per-probe threshold |
| `ct_qc_missing_sample` (+ `_custom`) | `selectInput` (+ conditional `numericInput`) | 0% / 5% / 10% / 20% / Custom | 10% | `methyl_fs_sample_missing_ok()` max-NA-per-sample threshold |
| `ct_qc_sexchr` | `selectInput` | Autosomes only (remove chrX/chrY) / Autosomes + chrX (remove chrY only) / Keep all | "Autosomes only" (`remove_xy`) | `methyl_filter_sex_chr()` mode |
| `ct_qc_beta_range` | `checkboxInput` | — | `TRUE` | **Declared but not read anywhere in the server code** — see `03_functions_and_code_audit.md` / `07_code_audit_findings.md` |
| `ct_qc_snp` | `checkboxInput` | — | `FALSE` | Gates a call to `methyl_filter_snp()` |
| `ct_qc_crossreactive` | `checkboxInput` | — | `FALSE` | Gates a call to `methyl_filter_cross_reactive()` |
| `ct_qc_crossreactive_file` | `fileInput` | .csv/.tsv/.txt | none | shown only when `ct_qc_crossreactive` is checked; parsed via `methyl_parse_probe_list()` |
| `ct_qc_detp_ui` → `ct_qc_detp` | `uiOutput` → conditional `checkboxInput` | — | `FALSE` | Only rendered as an actual checkbox when `ct_source()$detp` is non-NULL (i.e. an IDAT-derived dataset); otherwise shows an explanatory "not available" note. Gates `methyl_filter_detection_p()` |

**Output:** `ct_qc_cascade_ui` (probe-count summary text + `ct_qc_cascade_table`, a `DT::dataTableOutput` of the step-by-step retention cascade built by `methyl_probe_retention_cascade()`). **Note:** this cascade table's step-by-step "retained" counts are computed incorrectly for every step after the first — see Finding in `07_code_audit_findings.md`; the final CpG × sample summary line above the table is computed independently and correctly.

**Connection to other tabs:** every other tab's computation is downstream of `ct_filtered()$mat`, this tab's sole output object.

---

## Tab 2 — "2. CpG Feature Selection" (`value = "featsel"`, `mod_methyl_celltype_featsel_ui()`, `mod_methyl_celltype.R:90-130`)

**Purpose:** rank and select marker CpGs off the chosen reference's centroids (or a custom/uploaded list), restricting the CpG set used for deconvolution. Optional — deconvolution can run without ever visiting this tab (falls back to the full reference CpG set).

| Input | Type | Choices | Default | Conditional | Server use |
|---|---|---|---|---|---|
| `ct_fs_method` | `selectInput` | Reference-library markers (`reference`) / Variance-based CpGs (`variance`) / Custom CpG list (`custom`) / "Differential methylation markers / DMCs / DMRs (unavailable)" (`dmc_unavailable`) | `reference` | — | Branches `fs_result()`'s logic entirely |
| `ct_fs_dbeta_ui` → `ct_fs_dbeta` | `uiOutput` → `selectInput` | 0.05/0.10/0.15/0.20/Custom | 0.10 | shown only for `reference` method | Minimum `effect` threshold passed to `methyl_ct_select_markers()` |
| `ct_fs_topn` (+ `_custom`) | `selectInput` (+ conditional `numericInput`, min 1) | 50/100/200/333/500/1000/2000/Custom | 200 | — | `top_n` for `methyl_ct_top_n_balanced()` |
| `ct_fs_direction` | `selectInput` | Both / Hyper-methylated / Hypo-methylated | Both | — | `direction` filter in `methyl_ct_select_markers()` |
| `ct_fs_specificity` | `selectInput` | All available / Cell-type-specific / Shared | All | — | `specificity_mode` (median-split) in `methyl_ct_select_markers()` |
| `ct_fs_chr_scope` | `selectInput` | Autosomes only / Autosomes + X / All chromosomes | All | — | `methyl_ct_chr_allowed_ids()` scope |
| `ct_fs_custom_file` | `fileInput` | .csv/.tsv/.txt | none | shown only for `custom` method | Parsed via `methyl_parse_probe_list()`, intersected with the ranked CpG set |
| (static note) | `p` | — | — | — | States the FDR column is not applicable — no per-sample replicates |
| `ct_fs_run_btn` | `actionButton` | — | — | — | Fires `fs_result <- eventReactive(...)` |

**Outputs (all gated behind "has it run yet" via `fs_has_run`):** `ct_fs_bar_plot` (CpGs retained per cell type), `ct_fs_heatmap_plot` (marker × cell-type centroid heatmap, ≤200 rows shown), `ct_fs_scatter_plot` (effect vs. specificity, plotly), `ct_fs_table` (`DT`, with `P-value`/`FDR` columns hard-coded to `"n/a"`/`"n/a (no replicates)"` rather than a fabricated statistic), `ct_fs_download` (CSV of the selected marker table).

**Connection to other tabs:** if this tab has been run (`fs_has_run() == TRUE`), Tab 3's overlap-QC and Tab 4's deconvolution use `fs_result()$selected$cpg` as the active marker set (`ct_active_markers()`); otherwise they fall back to every CpG in the currently selected reference (`mod_methyl_celltype.R:637-650`).

---

## Tab 3 — "3. Reference & Method" (`value = "refmethod"`, `mod_methyl_celltype_refmethod_ui()`, `mod_methyl_celltype.R:132-209`)

**Purpose:** choose the reference matrix, the estimation method, advanced solver parameters, and see the reference/working-matrix CpG-overlap QC gate that ultimately enables or disables the "Run Deconvolution" button.

### Box "Reference library"
| Input | Type | Choices | Default | Conditional | Server use |
|---|---|---|---|---|---|
| `ct_ref_source` | `radioButtons` | Built-in reference library (`registry`) / Custom reference matrix (upload) (`custom`) | `registry` | — | Branches `active_reference_full()` |
| `ct_ref_id` | `selectInput` | 7 registry entries (labels from `methyl_ct_reference_registry()`) | `blood7` | shown for `registry` | `methyl_ct_get_reference(ct_ref_id)` |
| `ct_ref_celltypes_ui` → `ct_ref_celltypes_keep` | `uiOutput` → `checkboxGroupInput` | all cell-type columns of the selected reference | all selected | shown for `registry` (rendered dynamically once a reference resolves) | Subsets `active_reference()`'s columns; a `validate(need(length(keep) >= 2, ...))` blocks fewer than 2 |
| `ct_custom_ref_file` | `fileInput` | .csv/.tsv/.txt | none | shown for `custom` | `methyl_ct_parse_custom_reference()` |
| `ct_custom_ref_preview_ui` | `uiOutput` | — | — | shown for `custom` | Row/column-count preview or a validation error |

### Box "Deconvolution method"
| Input | Type | Choices | Default | Conditional | Server use |
|---|---|---|---|---|---|
| `ct_method` | `selectInput` | CP – Houseman constrained projection / EpiDISH RPC / EpiDISH CBS / Two-stage (hepidish) – advanced | first choice, `CP` (no `selected=` given, so Shiny defaults to the first item) | — | Branches `decon_result()`'s call between `methyl_ct_run_epidish()` and `methyl_ct_run_hepidish()` |
| `ct_hepidish_ic_col` | `selectInput` | dynamically populated with the active reference's column names | `"IC"` if present, else the last column | shown for `hepidish` | Passed as `ic_column`/`h.CT.idx` to hepidish |
| `ct_hepidish_ref2` | `selectInput` | blood-tissue registry entries only | `blood7` | shown for `hepidish` | Second-stage reference; note the hepidish call always fixes its **internal** EpiDISH sub-method to `"RPC"` (`mod_methyl_celltype.R:718-720`) — there is no UI control for CP/CBS inside hepidish; see `07_code_audit_findings.md` |
| (static list) | `div` | — | — | — | `methyl_ct_unavailable_methods()` rendered as a disabled-methods explanation list |

### Box "Reference-library QC"
| Input | Type | Choices | Default | Server use |
|---|---|---|---|---|
| `ct_overlap_threshold` | `selectInput` | 50% / 60% / 70% / 80% / 90% | 50% | Minimum fraction of active markers that must be present in the filtered working matrix (`ct_overlap_ok()`) to enable the Run button |

**Output:** `ct_refqc_ui` — matched-marker percentage message plus, when available, `ct_refqc_by_type_table` (per-cell-type overlap breakdown, `DT`).

### `<details>` "Advanced Parameters" (collapsed by default)
| Input | Type | Default | Server use |
|---|---|---|---|
| `ct_adv_maxit` | `numericInput` (min 10, max 500) | 50 | `maxit` for RPC's robust regression / hepidish's stage-1 RPC |
| `ct_adv_constraint` | `selectInput` | Inequality (non-negative) | `inequality` | Passed as `constraint` to `epidish()`/`hepidish()` — **only actually consumed by EpiDISH's own dispatch when `method == "CP"`**; see `07_code_audit_findings.md` |
| `ct_adv_seed` | `numericInput` | 1234 | `set.seed()` before the estimation call — **CP/RPC/CBS/hepidish are all deterministic (no RNG use found in `EpiDISH`'s `DoCP`/`DoRPC`/`DoCBS`)**, so this control currently has no observable effect on the result; see `07_code_audit_findings.md` |
| `ct_adv_nu1`/`nu2`/`nu3` | `numericInput` (min 0, max 1, step 0.05) each | 0.25/0.5/0.75 | `nu.v` for CBS's SVR tuning grid |

**Connection to other tabs:** feeds `active_reference()` (Tab 2's marker ranking and Tab 4's deconvolution both consume it) and gates whether Tab 4's Run button is enabled at all (via `ct_overlap_ok()`).

---

## Tab 4 — "4. Deconvolution" (`value = "deconv"`, `mod_methyl_celltype_deconv_ui()`, `mod_methyl_celltype.R:211-220`)

**Purpose:** the actual estimation run.

| Input | Type | Server use |
|---|---|---|
| `ct_run_decon_btn` | `actionButton`, disabled/enabled live by `shinyjs::disable`/`enable` based on `ct_overlap_ok()` (`mod_methyl_celltype.R:689-692`) | Fires `decon_result <- eventReactive(...)`, which independently re-derives `ref_use`/`markers` and re-validates (`nrow(ref_use) >= 10`) rather than trusting the button's disabled state alone |

**Outputs (behind `decon_has_run`):** `ct_decon_table` (long-format sample × cell-type × fraction, `DT`), `ct_decon_summary_table` (min/max/mean/median/sd per cell type), `ct_decon_download` (CSV).

**Side effect:** `observeEvent(decon_result(), ...)` writes a compact summary into the shared `results$celltype` object when a `results` reactiveValues was passed in (`mod_methyl_celltype.R:730-739`) — see `04_data_flow_and_pipeline.md` for its one confirmed downstream reader.

**Connection to other tabs:** every one of Tabs 5, 6, 7 is gated on `decon_has_run()`; none of them recompute deconvolution themselves (Compare Methods reruns EpiDISH per method but on the same `decon_result()$working_mat`/`ref_used`, not a fresh filter/reference selection).

---

## Tab 5 — "5. Cell Composition" (`value = "composition"`, `mod_methyl_celltype_composition_ui()`, `mod_methyl_celltype.R:222-224`)

**Purpose:** visualize and statistically compare the estimated fraction matrix. Entirely gated behind `decon_has_run()` (`ct_composition_gate`); its real UI (`ct_composition_ui`) is only rendered afterward.

### A. Stacked bar — cell composition per sample
| Input | Type | Choices | Default |
|---|---|---|---|
| `ct_bar_order` | `selectInput` | As in data / Alphabetical / By dominant cell type | As in data |
| `ct_bar_hide` | `checkboxGroupInput` | every cell type in the fraction matrix | all selected |
Output: `ct_bar_plot` (plotly-wrapped stacked bar, via `plotly_safe()`), `ct_bar_download` (PNG).

### B. Cell-type heatmap
| Input | Type | Default |
|---|---|---|
| `ct_heat_cluster_rows` | `checkboxInput` | `TRUE` |
| `ct_heat_cluster_cols` | `checkboxInput` | `TRUE` |
| `ct_heat_normalize` | `checkboxInput` (row-normalize) | `FALSE` |
Output: `ct_heat_plot`, `ct_heat_download` (PNG).

### C. Cell-type distribution
| Input | Type | Choices | Default |
|---|---|---|---|
| `ct_box_kind` | `radioButtons` | Boxplot / Violin | Boxplot |
| `ct_box_group` | `selectInput` | None + sample-sheet columns | None |
Output: `ct_box_plot`, `ct_box_download` (PNG).

### D. PCA / MDS of cell composition
| Input | Type | Choices | Default |
|---|---|---|---|
| `ct_ord_method` | `radioButtons` | PCA / MDS | PCA |
| `ct_ord_color` | `selectInput` | None + one entry per cell type + one entry per sample-sheet column (built with a length-guard around `sheet_cols`, see `01`/`07`) | None |
| `ct_ord_labels` | `checkboxInput` | `FALSE` |
Output: `ct_ord_plot` (built by `methyl_ct_composition_pca()`/`methyl_ct_composition_mds()`, deliberately not `qc.R`'s CpG-scale PCA/MDS — see `03_functions_and_code_audit.md`), `ct_ord_download` (PNG).

### E. Correlation matrix between cell types
No inputs beyond the fraction matrix itself. Output: `ct_corr_plot` (Pearson correlation across cell types via `stats::cor()`), `ct_corr_download` (PNG).

### Group Comparison
| Input | Type | Choices | Default |
|---|---|---|---|
| `ct_cmp_group_col` | `selectInput` | sample-sheet columns, or a disabled "No sample sheet loaded" placeholder | — |
| `ct_run_cmp_btn` | `actionButton` | — | — |
Output (behind `cmp_has_run`): test-used message (Wilcoxon or Kruskal-Wallis, auto-selected by group count), `ct_cmp_plot` (boxplot with significance stars), `ct_cmp_table` (`DT`), `ct_cmp_download` (CSV).

**Connection to other tabs:** entirely downstream of Tab 4's `decon_result()`; contributes the "Phenotype-linked fractions" and "Group comparison results" exports on Tab 7.

---

## Tab 6 — "6. Validation" (`value = "validation"`, `mod_methyl_celltype_validation_ui()`, `mod_methyl_celltype.R:226-242`)

**Purpose:** two independent sanity checks on the deconvolution result — reconstruction accuracy, and cross-method agreement.

### Reconstruction validation
| Input | Type |
|---|---|
| `ct_run_val_btn` | `actionButton` |
Output (behind `val_has_run`): 4 summary cards (Correlation/RMSE/MAE/R², from `methyl_ct_validation_metrics()`), `ct_val_plot` (observed-vs-reconstructed scatter), `ct_val_download` (PNG), `ct_val_table` (per-sample metrics, `DT`).

### Compare Methods
| Input | Type | Choices | Default |
|---|---|---|---|
| `ct_cmpmethods_pick` | `checkboxGroupInput` | CP / RPC / CBS | all 3 selected |
| `ct_run_cmpmethods_btn` | `actionButton` | — | — |
Output (behind `cmpm_has_run`): `ct_cmpmethods_corr_plot` (method-correlation heatmap), a pair-selection (`ct_ba_a`/`ct_ba_b`, `selectInput`s) driving `ct_cmpmethods_ba_plot` (Bland-Altman), `ct_cmpmethods_summary_table` (mean/max absolute-difference per method pair).

**Connection to other tabs:** reruns EpiDISH on `decon_result()$working_mat`/`ref_used` (i.e. the exact matrix/reference from Tab 4's run, not re-derived from Tabs 1–3's current input state) — so if a user changes a filter on Tab 1 *after* running Tab 4, Tabs 5/6 still reflect the prior run until Tab 4 is re-run.

---

## Tab 7 — "7. Export" (`value = "export"`, `mod_methyl_celltype_export_ui()`, `mod_methyl_celltype.R:244-263`)

**Purpose:** download every intermediate/final object produced by this module.

| Download | Gate | Content |
|---|---|---|
| `ct_export_beta` | always enabled (`req(f)` inside the handler) | `ct_filtered()$mat`, the QC-filtered working matrix |
| `ct_export_markers` | `fs_has_run()` | `fs_result()$selected`, the selected marker table |
| `ct_export_ref` | `decon_has_run()` | `decon_result()$ref_used`, the reference actually used |
| `ct_export_fractions` | `decon_has_run()` | `decon_result()$fractions` |
| `ct_export_pheno_fractions` | `decon_has_run()` | fractions joined to sample-sheet columns via `methyl_sheet_sample_ids()` |
| `ct_export_comparison` | `cmp_has_run()` | `cmp_result()$stats$table` |
| `ct_export_report` | always enabled | a plain-text summary assembled from whichever steps have run |

Button enabled/disabled state for the has-run-gated downloads is kept in sync via one `observe({ shinyjs::toggleState(...) })` block (`mod_methyl_celltype.R:1119-1125`).

---

## Module-level entry points not tied to one tab

- Sidebar navigation: `server.R:541-543` defines `observeEvent(input$sidebar_nav_methylomics_celltype, { jump_to_mx_submodule("celltype", ...) })`, with a comment claiming this exists so "Quality Control's Cell Composition tab" can link to this module via a plain `actionLink`. **No such `actionLink`/input with id `sidebar_nav_methylomics_celltype` was found anywhere in `mod_methyl_qc.R` or any other file in the repository** (verified by an exhaustive grep) — this shortcut handler is effectively dead code; the module remains reachable only through the normal Methylomics Sub-modules grid. This is outside the two core files audited here (it lives in `server.R`), so it is reported as an informational finding rather than treated as part of this module's own contract.
