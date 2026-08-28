# Cross-Omics — End-to-End Function Documentation

Scope: the entire Cross-Omics module as mounted by `ui.R`'s `crossomicsUI()` — the Dataset tab
(`mod_cross_dataset.R`) and the three registered sub-modules in `CX_MODULES`
(`submodules_registry.R:65-69`): Expression and Methylation (`mod_cross_integration.R`), Biomarker
Convergence (`mod_cross_biomarker_conv.R`), and Cross-Omics MR (`mod_cross_mr_stage.R`), plus every
helper file they call (`crossomics_integration_upload.R`, `crossomics_integration_helpers.R`,
`crossomics_integration_plots.R`, `crossomics_biomarkerconv_helpers.R`,
`crossomics_mrstage_helpers.R`). No filenames or functions were assumed; every entry below was
located by reading the named file.

> **Number of Cross-Omics tabs/sub-tabs identified from the inspected UI code: 21** total
> `tabPanel()` UI elements across the module, organized in three levels (see Section 1).

---

## 1. Tab inventory (exact UI order, from the actual `tabPanel()`/`tabsetPanel()` calls)

### Level 1 — `crossomicsUI()` (`ui.R:1481-1521`), `tabsetPanel(id = "cx_menu")`

1. **Dataset** (`ui.R:1503`) — `mod_cross_dataset_ui("cx_dataset")`. No further tabs inside it.
2. **Sub-modules** (`ui.R:1504-1516`) — hosts `build_submodule_grid(CX_MODULES, ...)`, a card grid
   (not itself a `tabPanel` set) grouped by `CX_SUBMODULE_GROUP_ORDER = c("Data", "Genetics")`
   (`ui.R:1470`). Clicking a card opens that sub-module, each of which has its own internal
   `tabsetPanel` (Level 2 below).

### Level 2 — the three `CX_MODULES` sub-modules, in registry order (`submodules_registry.R:66-68`)

**A. Expression and Methylation** (`mod_cross_integration_config`, group "Data") — 7 tabs
(`mod_cross_integration.R:52-106`, `tabsetPanel(id = ns("result_tabs"))`):
1. Expression data
2. Methylation data
3. Integration (setup form + Advanced Filters, no internal sub-tabs of its own)
4. Quadrant plot
5. Heatmap
6. Network analysis
7. Export

**B. Biomarker Convergence** (`mod_cross_biomarker_conv_config`, group "Data") — 4 tabs
(`mod_cross_biomarker_conv.R:58-64`, `tabsetPanel(id = ns("result_tabs"))`):
1. eQTL-MR
2. mQTL-MR
3. eQTL-mQTL
4. Downloads

**C. Cross-Omics MR** (`mod_cross_mr_stage_config`, group "Genetics") — 8 tabs
(`mod_cross_mr_stage.R:24-89`, `tabsetPanel(id = ns("result_tabs"))`): 5 tabs generated
dynamically, one per entry of `CX_MR_CATEGORIES` (`crossomics_mrstage_helpers.R:125-143`), plus 3
static tabs:
1. DEG-DMP-QTL
2. DEG-DMR-QTL
3. DEG-eQTL
4. DMP-mQTL
5. DMR-mQTL
6. Results Table
7. Volcanoplot
8. Downloads

**Total**: 2 (Level 1) + 7 + 4 + 8 (Level 2) = **21 `tabPanel()` elements**, confirmed by a direct
`grep -n "tabPanel("` across `ui.R` and every file in `R/crossomics/`.

---

## 2. Tab connection map (actual objects, not generic placeholders)

```text
Dataset tab (mod_cross_dataset.R)
   |  writes: cross_dataset$user_expr_df / user_expr_source / user_expr_wide /
   |          user_expr_mapping / user_expr_sample_cols
   |          user_meth_df / user_meth_source / user_meth_wide /
   |          user_meth_mapping / user_meth_sample_cols
   v
cross_dataset (shared reactiveValues store, created in server.R:104-107)
   |
   v
Expression and Methylation (mod_cross_integration.R) — the ONLY sub-module that reads cross_dataset
   |  mirrors cross_dataset -> raw$expr_df / raw$meth_df / raw$expr_sample_cols / raw$meth_sample_cols
   |  "Run Integration" click:
   |    cx_harmonize_gene_ids() -> id-harmonized gene symbols on both sides
   |    cx_aggregate_methylation() -> one row per gene (meth_h -> agg$df, agg$cpg_level)
   |    merge(expr_j, meth_j, by="gene") -> joined
   |    cx_classify() -> classified (adds sig_expression/sig_methylation/category/...)
   |    [if raw$expr_sample_cols & raw$meth_sample_cols overlap >= 3]
   |        cx_gene_correlation() -> correlation_r/correlation_p/correlation_fdr merged in
   |    cx_classify_evidence() -> evidence_level
   |  writes: integ$df (the per-run integrated table), integ$cpg_level, integ$pairing,
   |          integ$provenance, integ_by_sex$<all|female|male>
   |  publishes: cross_results$integration = list(df=classified, summary=..., provenance=..., params=...)
   v
Quadrant plot / Heatmap / Network analysis / Export tabs
   (all read integ$df / filtered_df() directly — no further shared-store hand-off)

Biomarker Convergence (mod_cross_biomarker_conv.R) — INDEPENDENT of cross_dataset/Dataset tab
   |  "Load Table": cx_bc_load_precomputed(sex) -> cross_omics_eQTL_mQTL_{sex}.csv (CX_RESULTS_DIR)
   |                 + cx_bc_backfill_mqtl_from_mrstage() (reads mr_stage_eqtl_significant_genes_mqtl_mr.csv)
   |  OR "Merge & Load": cx_bc_load_eqtl_upload()/cx_bc_load_mqtl_upload() -> cx_bc_merge_eqtl_mqtl()
   |  bc_df() = cx_bc_relabel(raw$df, CX_BC_DEFAULT_PARAMS) -> *_significant / n_evidence_layers columns
   |  publishes: cross_results$biomarkerconv = list(df=bc_df(), sex=..., run_at=...)
   v
eQTL-MR / mQTL-MR / eQTL-mQTL / Downloads tabs (subset bc_df() by in_eQTL_MR_panel / in_mQTL_MR_panel)

Cross-Omics MR (mod_cross_mr_stage.R) — INDEPENDENT of cross_dataset/Dataset tab
   |  "Load MR Results": cx_mr_load_precomputed() -> mr_stage_eqtl_significant_genes_mqtl_mr.csv
   |  join_df(): reuses cross_results$biomarkerconv$df if already loaded for the same sex,
   |             else calls cx_bc_load_precomputed()+cx_bc_relabel() itself (no dependency on
   |             the Biomarker Convergence tab having been opened first)
   |  categories() = cx_mr_classify_categories(join_df()) -> one data.frame per CX_MR_CATEGORIES entry
   |  publishes: cross_results$mrstage = list(df=mrs$df, categories=categories(), run_at=...)
   v
DEG-DMP-QTL / DEG-DMR-QTL / DEG-eQTL / DMP-mQTL / DMR-mQTL / Results Table / Volcanoplot / Downloads
```

**Key structural fact, verified against `mod_cross_biomarker_conv.R` and `mod_cross_mr_stage.R` in
full**: only **Expression and Methylation** actually reads `cross_dataset`. Biomarker Convergence
and Cross-Omics MR each accept `cross_dataset` as a function parameter (matching the shared
`(id, cross_dataset, cross_results, ...)` call signature `server.R:122-130` uses for every
`CX_MODULES` entry) but never reference it in their server bodies — both load their own
independent precomputed-file or upload data instead. The Dataset tab therefore is the entry point
for exactly **one** of the three sub-modules, not all three.

---

## 3. Function documentation

Template per Section 19 of the audit brief. Functions are grouped by file, in the order they
appear in that file.

### 3.1 `mod_cross_dataset.R` — the Dataset tab

```text
Function: mod_cross_dataset_ui(id)
File: R/crossomics/mod_cross_dataset.R:48-96
Package: shiny (tagList/fluidRow/box/radioButtons/fileInput/actionButton/uiOutput)
Called from: ui.R:1503 (crossomicsUI())
Purpose: Builds the Dataset tab's UI — the "1. Data Source" panel and the "Preview" panel.
Input: a Shiny module id ("cx_dataset").
Input structure: character(1).
Processing: pure UI construction; no data processing.
Statistical/computational principle: none (UI layer only).
Output: a Shiny tagList of UI elements.
Downstream dependency: none (rendered once at app start).
Why it matters scientifically: defines the only two ways real data can enter the Cross-Omics
  module's "Expression and Methylation" analysis.
QC: none applicable (UI definition).
Audit status: PASS — matches the server logic it renders inputs for.
```

```text
Function: mod_cross_dataset_server(id, cross_dataset)
File: R/crossomics/mod_cross_dataset.R:98-263
Package: shiny (moduleServer/reactiveVal/observeEvent/renderUI/DT::renderDataTable/renderPlot)
Called from: server.R:109
Purpose: Loads example or uploaded Transcriptomics/Methylomics data, standardizes it, previews it,
  and on "Use this data" publishes it into the shared cross_dataset store.
Input: `id` (module id); `cross_dataset`, a `reactiveValues` object created once in server.R and
  shared with every CX_MODULES sub-module.
Input structure: `cross_dataset` fields are described in Section 1 above.
Processing: see Cross-Omics_Dataset_Tab_Code_Audit.md Sections 7-13 for the full traced logic.
Statistical/computational principle: none itself — delegates all statistics-adjacent work
  (standardization, region bucketing, sample-column detection) to crossomics_integration_helpers.R.
Output: side effects on `cross_dataset` (via observeEvent), plus 5 Shiny outputs (preview_ui,
  expr_table, expr_plot, meth_table, meth_plot).
Downstream dependency: mod_cross_integration.R's `observe()` block (line 142) reads cross_dataset.
Why it matters scientifically: is the sole staging point that determines what "gene", "log2fc",
  "dbeta", "pvalue", "fdr" mean for every downstream Expression-and-Methylation computation.
QC: none of its own (see Section 9 of the audit file) — relies entirely on upstream QC for
  Example data, and on cx_standardize_*()'s structural checks for uploads.
Audit status: PASS for correctness of the standardization hand-off; WARNING for the absence of a
  beta/Δβ plausibility check on upload (see audit findings).
```

### 3.2 `crossomics_integration_upload.R`

```text
Function: cx_read_table(datapath, filename)
File: R/crossomics/crossomics_integration_upload.R:13-36
Package: data.table (fread), openxlsx (read.xlsx), tools (file_ext)
Called from: cx_read_and_detect() (this file); cx_bc_load_eqtl_upload()/cx_bc_load_mqtl_upload()
  (crossomics_biomarkerconv_helpers.R); cx_mr_load_upload()/cx_mr_load_evidence_upload()
  (crossomics_mrstage_helpers.R) — i.e. every upload path in the Cross-Omics module reuses this
  one reader.
Purpose: Parses a CSV/TSV/TXT/XLSX file into a plain data.frame.
Input: `datapath` (temp file path from a Shiny fileInput), `filename` (original name, for
  extension detection).
Input structure: character(1) each.
Processing: dispatches on file extension; fread() with a fixed NA-string vocabulary for
  delimited text, openxlsx::read.xlsx(sheet=1) for Excel; rejects results with <2 columns or 0 rows.
Statistical/computational principle: none — pure I/O.
Output: list(ok, df, error) — fail-soft, never throws.
Downstream dependency: every standardize/detect function that follows in the calling module.
Why it matters scientifically: the single point of truth for "what counts as a parseable input
  file" across every Cross-Omics upload path.
QC: file-existence (implicit via tryCatch), row/column-count floor, unsupported-extension message.
  Does not check file size or encoding.
Audit status: PASS.
```

```text
Function: cx_read_and_detect(datapath, filename, kind)
File: R/crossomics/crossomics_integration_upload.R:40-46
Package: none beyond cx_read_table()/cx_detect_columns()
Called from: mod_cross_dataset.R (expr_file/meth_file observers)
Purpose: Bundles file reading with column auto-detection for one upload event.
Input: datapath, filename, kind ("expression"/"methylation").
Input structure: character(1) x3.
Processing: cx_read_table() then cx_detect_columns(df, kind).
Statistical/computational principle: none.
Output: list(ok, df, mapping, error).
Downstream dependency: cx_standardize_expression()/cx_standardize_methylation().
Why it matters scientifically: first point at which a raw upload's columns are interpreted as
  "gene"/"log2fc"/"dbeta"/etc.
QC: inherits cx_read_table()'s checks; adds none of its own.
Audit status: PASS.
```

### 3.3 `crossomics_integration_helpers.R`

```text
Function: cx_match_column(cols, patterns, exclude)
File: crossomics_integration_helpers.R:39-48
Purpose: Finds the first column name matching any of a set of case-insensitive regex patterns.
Input: cols (character vector of column names), patterns (character vector of regexes),
  exclude (already-claimed column names).
Processing: lower-cases/trims, iterates patterns in priority order, returns the first hit.
Output: one column name or NA_character_.
Called from: cx_detect_columns().
Audit status: PASS — deterministic, order-dependent by design (first pattern in CX_FIELD_PATTERNS
  wins), documented as such.
```

```text
Function: cx_detect_columns(df, kind)
File: crossomics_integration_helpers.R:55-67
Purpose: Builds a full field->column mapping (gene, log2fc/dbeta/beta, pvalue, fdr, chr, pos,
  region, island, cpg, sample_id) for an expression or methylation table.
Input: df (data.frame), kind ("expression"/"methylation").
Processing: greedy, order-independent per field — a column claimed by an earlier field is excluded
  from later fields, so "fdr" and "pvalue" cannot both grab the same ambiguous column.
Output: named character vector, one entry per canonical field (value or NA).
Called from: cx_read_and_detect().
Statistical/computational principle: pattern-matching heuristic, not a statistical procedure.
Audit status: WARNING (informational) — greedy first-match-wins detection can mis-assign a
  column when a file's headers are ambiguous or non-standard; the only downstream safeguard is
  cx_standardize_*()'s hard requirement for gene+log2fc/dbeta, so a wrong PVALUE/FDR assignment
  would not be caught, only a wrong GENE/LOG2FC/DBETA assignment (which fails standardization) would.
```

```text
Function: cx_detect_sample_columns(df, mapping)
File: crossomics_integration_helpers.R:79-85
Purpose: Identifies which columns of a wide upload look like per-sample numeric measurements.
Input: df (the original wide upload), mapping (the detected field mapping).
Processing: candidates = all columns not already claimed by `mapping` and not in the hardcoded
  CX_NON_SAMPLE_NUMERIC_NAMES exclusion list (aveexpr, t, b, z, score, n, meandiff, ...); keeps
  those that are numeric (or numeric-coercible character).
Output: character vector of candidate sample-column names.
Called from: mod_cross_dataset.R (use_data_btn handler).
Why it matters scientifically: the sole basis for whether sample-level correlation can later be
  computed (cx_detect_sample_pairing()).
Audit status: PASS with a documented limitation — a genuinely non-sample numeric column not in the
  hardcoded exclusion list (e.g. a custom summary statistic) would be misread as a sample.
```

```text
Function: cx_as_numeric_safe(x)
File: crossomics_integration_helpers.R:91
Purpose: as.numeric() with warnings suppressed, so a non-numeric cell becomes NA quietly instead
  of emitting a console warning per call.
Audit status: PASS (thin wrapper).
```

```text
Function: cx_standardize_expression(df, mapping)
File: crossomics_integration_helpers.R:96-110
Purpose: Builds the standardized gene, log2fc, pvalue, fdr table from a raw table + column mapping.
Input: df (raw data.frame), mapping (named character vector from cx_detect_columns() or a manual
  override).
Processing: requires gene+log2fc mappings present; numeric-coerces log2fc/pvalue/fdr; drops
  blank/NA-gene rows; requires >=1 remaining row; deduplicates by gene via cx_dedup_by_gene().
Statistical/computational principle: none beyond deterministic dedup.
Output: list(ok, df) or list(ok=FALSE, error).
Called from: mod_cross_dataset.R (both example and upload paths).
Downstream dependency: the resulting df becomes cross_dataset$user_expr_df.
QC: presence checks, numeric coercion, non-blank-gene filter, non-empty-result requirement.
Audit status: PASS.
```

```text
Function: cx_dedup_by_gene(df)
File: crossomics_integration_helpers.R:116-126
Purpose: Collapses duplicate gene symbols on the expression side to one row per gene.
Input: a standardized-shape expression data.frame with possible duplicate `gene` values.
Processing: for each duplicated gene, keeps the row with the smallest FDR; if all FDR are NA,
  falls back to smallest p-value; if both are all-NA, keeps the first row.
Statistical/computational principle: deterministic "most significant row wins" reduction —
  explicitly the same idea as the methylation side's own "min FDR" aggregation option, applied
  here unconditionally (the expression side has no user-selectable aggregation method).
Output: one row per gene.
Called from: cx_standardize_expression(); also re-applied in mod_cross_integration.R:197 after
  gene-ID harmonization (harmonizing IDs to a canonical symbol can newly create duplicates that
  did not exist before harmonization).
Audit status: PASS — deterministic, documented, consistent with the methylation-side convention.
```

```text
Function: cx_region_bucket(region_raw)
File: crossomics_integration_helpers.R:132-140
Purpose: Buckets a raw genomic-context string into "Promoter" / "Gene body" / "Other".
Input: character vector (Illumina UCSC_RefGene_Group vocabulary or free-text upload values).
Processing: dplyr::case_when() on lower-cased text: tss/promoter/5'utr/1stExon/upstream ->
  Promoter; body/exon/intron/3'utr/exonbnd -> Gene body; else Other; NA/blank -> NA.
Statistical/computational principle: none (vocabulary mapping).
Output: character vector, same length as input.
Called from: cx_standardize_methylation().
Audit status: PASS.
```

```text
Function: cx_region_fine(region_raw)
File: crossomics_integration_helpers.R:151-161
Purpose: Extracts the first, fine-grained UCSC_RefGene_Group token (TSS200/TSS1500/5'UTR/1stExon/
  Body/3'UTR/ExonBnd) from a semicolon-joined annotation string.
Processing: vectorized first-token extraction, trimws() applied once over the whole result vector
  (explicitly optimized — the code comments this used to cost ~16s per-element on a 319k-row table).
Output: character vector.
Called from: cx_standardize_methylation().
Audit status: PASS.
```

```text
Function: cx_standardize_methylation(df, mapping)
File: crossomics_integration_helpers.R:165-191
Purpose: Builds the standardized cpg, gene, dbeta, pvalue, fdr, chr, pos, region_raw, region,
  region_fine, island_context table.
Input: df, mapping.
Processing: requires a gene mapping; requires dbeta OR beta mapping; if no cpg mapping, synthesizes
  "row1","row2",... ids; numeric-coerces dbeta/pvalue/fdr/pos; derives region/region_fine via the
  two functions above; drops blank/NA-gene rows; requires >=1 remaining row.
Output: list(ok, df) or list(ok=FALSE, error).
Called from: mod_cross_dataset.R (both paths); cx_load_default_methylation();
  cx_load_default_dmr().
QC: presence checks, numeric coercion, non-blank-gene filter, non-empty-result requirement. No
  beta/Δβ plausible-range check (see Dataset Tab Audit Section 14).
Audit status: PASS structurally; WARNING for the missing value-range check.
```

```text
Function: CX_AGGREGATION_METHODS / cx_cpg_level_table(meth_std, meth_thresh, meth_fdr_thresh)
File: crossomics_integration_helpers.R:197-223
Purpose: Computes the un-collapsed, per-CpG "Level 2" table — direction (Hyper/Hypo) and a
  per-CpG significance flag — BEFORE any gene-level aggregation.
Input: the standardized methylation df; the Integration tab's Δβ and FDR thresholds.
Processing: sig_cpg = |dbeta|>=thresh & fdr<fdr_thresh (both non-NA); selects a fixed column subset.
Output: a data.frame, one row per CpG, never discarded (only ever summarized).
Called from: cx_aggregate_methylation().
Statistical/computational principle: simple two-condition significance flag, no correction here
  (correction happens on the aggregated, gene-level p-values downstream).
Audit status: PASS — this is the mechanism the module uses to satisfy "do not collapse CpGs for a
  gene without retaining their individual identities."
```

```text
Function: cx_cpg_counts_per_gene(cpg_level_df)
File: crossomics_integration_helpers.R:233-267
Purpose: Per-gene breakdown of how many CpGs map to it, by significance/direction/region/island.
Input: the CpG-level table above.
Processing: vectorized table()/tapply()/max.col() aggregation (explicitly chosen over a per-gene
  loop for performance — commented as ~45s -> <1s on ~19k genes/~320k CpGs).
Output: one row per gene: n_cpg_total, n_cpg_unique, n_cpg_significant, n_cpg_hyper/hypo_sig,
  n_cpg_hyper/hypo_all, per-region counts (TSS200/TSS1500/5'UTR/1stExon/Body/3'UTR), n_promoter/
  gene_body_region, per-island-context counts, primary_region (most frequent).
Called from: cx_aggregate_methylation() (merged onto the aggregated result).
Audit status: PASS — always computed from the full pre-aggregation CpG set, independent of which
  single aggregation method the user picked, so it is not itself distorted by that choice.
```

```text
Function: cx_aggregate_methylation(meth_std, method, meth_thresh, meth_fdr_thresh)
File: crossomics_integration_helpers.R:277-364
Purpose: Collapses multiple CpGs per gene into one gene-level dbeta/pvalue/fdr row, by one of 7
  user-selectable methods (CX_AGGREGATION_METHODS: mean, median, min_fdr, max_abs_dbeta,
  promoter_only, gene_body_only, promoter_and_body).
Input: standardized methylation df, aggregation method, both significance thresholds (used only to
  build the CpG-level significance flags, not the aggregation math itself).
Processing (traced in full):
  - For promoter/body-restricted methods, filters CpGs by `region` first.
  - For "min_fdr"/"max_abs_dbeta": picks ONE CpG's entire row (dbeta, pvalue, fdr, chr, pos,
    region all from that same CpG) — dbeta and its significance statistic are always co-derived
    from the same observation for these two methods.
  - For "mean"/"median" (and the three region-restricted variants, which reduce to mean once
    filtered): dbeta_v = mean()/median() of ALL CpGs' dbeta for that gene. Critically, pvalue_v is
    NOT taken as min(fdr)/min(pvalue) of a single CpG — the code explicitly computes a Stouffer's
    Z combination: each CpG's p-value is signed by its own dbeta's direction, converted to a
    z-score (qnorm), z-scores are averaged across all CpGs mapped to the gene, and the combined z
    is converted back to a two-sided p-value. For a single-CpG gene this exactly reduces to that
    CpG's own p-value. fdr is then recomputed via Benjamini-Hochberg ACROSS ALL GENES' combined
    p-values (the gene, not the CpG, is the multiple-testing unit for these methods).
  - Non-finite dbeta/fdr/pvalue/dbeta_mean/dbeta_median are coerced to NA.
  - cx_cpg_counts_per_gene() is merged on.
Statistical/computational principle: for mean/median, Stouffer's combined-probability method
  (directional, weighted equally per CpG) followed by gene-level Benjamini-Hochberg FDR; for
  min_fdr/max_abs_dbeta, single-observation selection (no combination).
Output: list(ok, df, note, cpg_level).
Called from: mod_cross_integration.R:201 (Run Integration).
Downstream dependency: integ$cpg_level (CpG-level detail modal); the joined/classified table.
Why it matters scientifically: this is the exact mechanism the audit brief (Section 27) asks to be
  checked for "reported effect and FDR may correspond to different CpGs." See Section 4 below for
  the explicit finding.
QC: none beyond the significance-threshold flags already described; no minimum-CpG-count floor.
Audit status: PASS for mean/median (statistically coherent pairing of effect and significance, see
  Section 4); PASS for min_fdr/max_abs_dbeta (single-CpG selection is internally consistent by
  construction, though it discards the other CpGs' effect sizes from the gene-level dbeta).
```

```text
Function: cx_classify(joined, expr_thresh, expr_fdr_thresh, meth_thresh, meth_fdr_thresh)
File: crossomics_integration_helpers.R:410-427
Purpose: Flags per-gene expression/methylation significance and assigns one of 5
  Hyper/Hypo x Up/Down categories (+ "Not significant").
Input: the expression-methylation joined table; both pairs of thresholds.
Processing: sig_expression/sig_methylation from |value|>=thresh & fdr<fdr_thresh (both non-NA);
  category assigned only when BOTH are significant, from the sign of dbeta/log2fc; direction
  labels derived purely from the sign of the already-supplied log2fc/dbeta (never reversed).
Statistical/computational principle: simple two-threshold, four-quadrant classification.
Output: joined df + sig_expression, sig_methylation, expression_direction, methylation_direction,
  category (factor), category_label (always phrased "potential"/"association", never causal).
Called from: mod_cross_integration.R:219.
Audit status: PASS — framing explicitly avoids causal language (CX_CATEGORY_LABELS,
  crossomics_integration_helpers.R:370-376).
```

```text
Function: cx_classify_evidence(classified_df, has_correlation)
File: crossomics_integration_helpers.R:450-466
Purpose: Assigns an Evidence Level tier (Strong/Moderate candidate, Expression-only,
  Methylation-only, Discordant, Insufficient evidence).
Input: the cx_classify()'d table; whether real sample-level correlation was computed.
Processing: "Strong candidate" additionally requires inverse direction AND a significant negative
  sample-level correlation (correlation_r<0, correlation_fdr<0.05) — and ONLY when
  has_correlation=TRUE; never inferred when correlation was not computed.
Output: an evidence_level factor column.
Called from: mod_cross_integration.R:232.
Audit status: PASS — explicitly guards against inferring "Strong candidate" without real
  sample-level evidence.
```

```text
Function: cx_filter_by_category / cx_filter_by_region / cx_filter_by_island /
          cx_filter_by_evidence / cx_filter_by_min_cpg / cx_filter_by_correlation_direction
File: crossomics_integration_helpers.R:473-522
Purpose: Six independent, composable, no-op-safe filters applied to the integrated table for the
  Expression data/Methylation data tabs' "Advanced Filters" and category-card clicks.
Input: the integrated df, plus one filter-specific selector value.
Processing: each is a pure in-memory subset; each returns df unchanged if its selector is
  empty/NULL or its required column is absent (never errors).
Output: filtered data.frame.
Called from: mod_cross_integration.R's filtered_df() reactive.
Audit status: PASS.
```

```text
Function: cx_build_gene_sample_matrix(df, id_col, sample_cols)
File: crossomics_integration_helpers.R:528-537
Purpose: Collapses a probe/CpG-or-gene-level wide matrix to one row per unique id (averaging
  duplicates via rowsum()/count), for genes x samples correlation input.
Called from: mod_cross_integration.R:223-224 (only reached when pairing$paired is TRUE).
Audit status: PASS.
```

```text
Function: cx_detect_sample_pairing(expr_samples, meth_samples, min_overlap=3)
File: crossomics_integration_helpers.R:550-557
Purpose: Determines whether real, paired sample-level data exists on both sides.
Processing: paired = TRUE only if both sides have >0 sample columns AND >=3 shared column names
  (intersect on raw column-name text — no normalization/fuzzy matching).
Output: list(paired, common_samples, n_expr, n_meth, n_common).
Called from: mod_cross_integration.R:221.
Why it matters scientifically: the single gate that prevents correlation/Evidence-Level-tiering
  from running on gene-level summary data with no real sample pairing.
Audit status: PASS — conservative by design (reports FALSE, never inferred TRUE, whenever either
  side lacks per-sample columns, which is always true for this app's own "Example data").
```

```text
Function: cx_gene_correlation(expr_mat, meth_mat, common_samples, method, min_n=3)
File: crossomics_integration_helpers.R:564-581
Purpose: Per-gene Pearson/Spearman correlation between matched expression and methylation samples.
Processing: stats::cor.test() per gene present in both matrices, restricted to common_samples,
  requiring >=min_n non-missing paired observations.
Output: data.frame(gene, r, p, n); FDR added by the caller.
Called from: mod_cross_integration.R:225.
Statistical/computational principle: standard parametric/rank correlation test, per gene, no
  correction for testing many genes until cx_adjust_p() is applied by the caller.
Audit status: PASS.
```

```text
Function: cx_adjust_p(p, method)
File: crossomics_integration_helpers.R:583-586
Purpose: Wraps stats::p.adjust() (BH or Bonferroni).
Audit status: PASS (thin wrapper).
```

```text
Function: cx_get_region_annotation(array_type="450K")
File: crossomics_integration_helpers.R:604-633
Purpose: Reads chr/pos/gene/UCSC_RefGene_Group/Relation_to_Island directly from the
  IlluminaHumanMethylation450kanno.ilmn12.hg19 (or EPIC) Bioconductor package, process-cached.
Processing: reads Locations+Other+Islands.UCSC via utils::data(); takes the FIRST semicolon-token
  of UCSC_RefGene_Name as the representative gene, and of UCSC_RefGene_Group as the representative
  region annotation.
Output: list(ok, anno, reason) — anno keyed by probe ID.
Called from: cx_load_default_methylation().
Why it matters scientifically: the CpG->gene annotation source for the Dataset tab's "Example
  data (CpG-level DMP)" path; a probe annotated to multiple genes loses every gene but the first.
QC: package-availability check; cache hit/miss; no annotation-completeness check beyond that.
Audit status: PASS with a disclosed limitation (multi-gene probes truncated to first gene — see
  Dataset Tab Audit Section 12), which the module itself surfaces to the user in its provenance text.
```

```text
Function: cx_load_default_deg(sex)
File: crossomics_integration_helpers.R:646-651
Purpose: Reads DATA_ROOT/results/tables/DEG_{sex}_full.csv.
Output: a raw data.frame, or NULL if the file is absent (never errors).
Called from: mod_cross_dataset.R (load_example_btn handler).
Audit status: PASS.
```

```text
Function: cx_load_default_methylation(sex, array_type="450K")
File: crossomics_integration_helpers.R:659-699
Purpose: Reads the preloaded SVA-adjusted, bacon-corrected DMP table for `sex` and annotates it to
  genes via cx_get_region_annotation().
Processing: guards on METH_DATA_AVAILABLE; calls load_default_dmp(stage="sva", sex); a single
  match() against the annotation table's rownames (explicitly optimized vs. per-column rowname
  indexing, commented as ~15-20s -> a few seconds); drops rows with no resolved gene; standardizes
  via cx_standardize_methylation().
Output: list(ok, df, error).
Called from: mod_cross_dataset.R (load_example_btn handler, meth_level="dmp").
Downstream dependency: becomes cross_dataset$user_meth_df for "Example data".
Why it matters scientifically: the sole path by which this app's own bundled, real (GSE42861-
  derived) sex-stratified DMP results reach the Cross-Omics module.
QC: file-existence, non-empty-result, gene-resolution requirements; no re-QC of the DMP statistics
  themselves (see Dataset Tab Audit Section 9).
Audit status: PASS.
```

```text
Function: cx_load_default_dmr(sex)
File: crossomics_integration_helpers.R:713-738
Purpose: Reads the preloaded DMRcate region table (METH_DMR_DIR/dmr_{sex}_full.csv) and
  standardizes it, using the table's own `overlapping.genes` column (no separate annotation call).
Output: list(ok, df, error).
Called from: mod_cross_dataset.R (load_example_btn handler, meth_level="dmr").
Audit status: PASS.
```

```text
Function: cx_build_provenance(params)
File: crossomics_integration_helpers.R:746-766
Purpose: Renders a human-readable list of every parameter/setting used in an Integration run
  (sex stratum, input mode, sources, thresholds, aggregation method, correlation method, sample
  matching status, annotation source(s), and the "first-gene-only" caveat when relevant).
Called from: mod_cross_integration.R:264.
Audit status: PASS — the mechanism by which Section 8's "clearly show how the aggregation was
  performed" and Section 14's QC-evidence traceability requirement are actually satisfied in the UI.
```

```text
Function: cx_build_id_lookup()
File: crossomics_integration_helpers.R:788-811
Purpose: One-time, process-cached SYMBOL->{ENTREZID,ENSEMBL} and ALIAS->SYMBOL tables from
  org.Hs.eg.db, via AnnotationDbi::keys()/select().
Audit status: PASS.
```

```text
Function: cx_harmonize_gene_ids(genes)
File: crossomics_integration_helpers.R:821-896
Purpose: Resolves a mixed vector of HGNC symbols / Entrez IDs / Ensembl Gene IDs to one canonical
  HGNC symbol per input, via EXACT matching only (no fuzzy string matching anywhere in this
  function).
Processing order: (1) exact symbol match, case-insensitive text only; (2) exact Entrez match;
  (3) exact Ensembl match (version suffix stripped, e.g. ".5"); (4) alias resolution — exactly one
  candidate symbol -> "alias_resolved", more than one -> "ambiguous" (every candidate listed, never
  guessed), zero -> stays "unmatched".
Output: list(ok, df(input_id, canonical_symbol, entrez_id, ensembl_id, match_type), summary).
Called from: mod_cross_integration.R:191 (Run Integration).
Statistical/computational principle: none (deterministic identifier resolution).
Audit status: PASS — the explicit design choice to prefer "ambiguous"/"unmatched" over a guess is
  directly verifiable in the code and matches the "no fuzzy matching" requirement.
```

```text
Function: cx_apply_harmonization(genes, harm_df)
File: crossomics_integration_helpers.R:903-908
Purpose: Rewrites gene IDs to their harmonized canonical symbol wherever an unambiguous match
  exists; leaves ambiguous/unmatched entries as their original text.
Called from: mod_cross_integration.R:196,198.
Audit status: PASS.
```

```text
Function: cx_compare_sexes(runs)
File: crossomics_integration_helpers.R:924-963
Purpose: Compares independently-run ALL/FEMALE/MALE Integration results gene-by-gene.
Processing: explicitly never computes or implies a formal sex x molecular-effect interaction test
  (no matched sample-level dataset with a sex covariate exists for that); differences are reported
  as "sex-stratified difference in significance", never "sex-specific biological effect".
Output: list(ok, table, summary, error); requires >=2 strata already run.
Called from: NOWHERE — confirmed by a repository-wide search (`grep -rl cx_compare_sexes
  R/`): the only files that reference this name are the file that defines it
  (crossomics_integration_helpers.R) and its own UI renderer, cx_sex_comparison_summary_ui()
  (crossomics_integration_plots.R). Neither is called from mod_cross_integration.R or from any
  other module server in the repository. There is no "Sex Comparison" tabPanel in the current
  Expression and Methylation UI (Section 1 above lists all 7 of its actual tabs).
Audit status: WARNING (dead code / stale intent) — implemented and internally correct in
  isolation, but not reachable from any current UI element. Do not describe a live "Sex
  Comparison" feature in thesis text without independently re-confirming it has since been wired in.
```

```text
Function: cx_validate_dataset(expr_df, meth_df, id_harmonization)
File: crossomics_integration_helpers.R:970-995
Purpose: Builds the plain-language "is this dataset ready" checklist (gene ID detected, log2FC
  detected, FDR detected / CpG ID detected, Δβ detected, gene annotation detected / gene overlap
  count / harmonization summary).
Called from: mod_cross_integration.R:236 — i.e. AFTER Run Integration, not on the Dataset tab.
Audit status: PASS in isolation; see Dataset Tab Audit Section 8/14 for the finding that it is not
  surfaced at the point where the data was actually loaded.
```

### 3.4 `crossomics_integration_plots.R`

```text
Function: cx_empty_state(message)
File: crossomics_integration_plots.R:14-16
Purpose: Standard "nothing to show yet" placeholder UI, reused by every result tab before a run.
Audit status: PASS.
```

```text
Function: cx_fmt_num(x, digits, sci) / cx_gene_detail_modal(row, pairing, cpg_rows)
File: crossomics_integration_plots.R:18-82
Purpose: Number formatting helper, and the per-gene detail modal (shows a clicked quadrant point's
  full row plus its individual, un-collapsed CpG rows via cpg_rows — i.e. the Level-2 detail the
  aggregation step does not discard).
Called from: mod_cross_integration.R's quadrant-click observer.
Audit status: PASS.
```

```text
Function: cx_quadrant_ggplot / cx_quadrant_plotly
File: crossomics_integration_plots.R:84-125
Purpose: Static (ggplot2, for download) and interactive (plotly, for on-screen) versions of the
  log2FC-vs-Δβ quadrant scatter, with optional gene highlighting/labeling and quadrant boundary
  lines at the current thresholds.
Called from: mod_cross_integration.R (quadrant_plot output and the 3 quadrant download handlers).
Statistical/computational principle: none beyond the already-computed classification; purely a
  visualization of cx_classify()'s output.
Audit status: PASS — the on-screen caption explicitly states "This shows a statistical association
  ... it does not establish that one causes the other" (mod_cross_integration.R:428).
```

```text
Function: cx_volcano_expr_plot / cx_volcano_meth_plot / cx_gene_correlation_plot /
          cx_genomic_view_plot
File: crossomics_integration_plots.R:127-251
Purpose: Additional plotting helpers (expression volcano, methylation volcano, one-gene
  expression-vs-methylation sample scatter, a genomic-position view).
Called from: NOWHERE — confirmed by a repository-wide search; none of these four names appears
  outside the file that defines them. mod_cross_integration.R's actual result tabs (Section 1) do
  not include a standalone Correlation or Genomic View tab, and its Expression/Methylation-data
  tabs render DT tables, not these plots.
Audit status: WARNING (dead code) — implemented but unreachable from any current UI element in
  this repository. Do not present these as active Cross-Omics outputs.
```

```text
Function: cx_gene_cpg_network_plot(sub_df)
File: crossomics_integration_plots.R:198-229 (approx.)
Purpose: Draws the static gene-CpG association network for the Network analysis tab.
Called from: mod_cross_integration.R's network_plot output.
Audit status: PASS (described in the module's own UI text as static/best-effort, since no
  interactive network-graphing package is installed in this deployment).
```

```text
Function: cx_run_pathway_enrichment(genes, universe, ontology)
File: crossomics_integration_plots.R:252-272
Purpose: GO (clusterProfiler::enrichGO) or KEGG (bitr + enrichKEGG) over-representation analysis.
Input: a gene set, a background universe, and an ontology choice.
Processing: KEGG path requires internet access to the KEGG REST API (explicitly caught and
  surfaced as an error if it fails); GO path uses org.Hs.eg.db SYMBOL keys directly.
Statistical/computational principle: hypergeometric/Fisher over-representation test with
  Benjamini-Hochberg correction (clusterProfiler defaults, pAdjustMethod="BH" explicitly set).
Output: list(ok, table, enrich_obj) or list(ok=FALSE, error).
Called from: NOWHERE — confirmed by a repository-wide search; no `tabPanel` named "Pathway"/
  "Enrichment" exists in mod_cross_integration.R's current tabsetPanel (Section 1 above lists the
  7 actually-mounted tabs, none of which is a Pathways tab).
Audit status: WARNING (dead code) — this function exists and is well-formed, but the Cross-Omics
  UI code as currently written does not expose a Pathways/Enrichment tab to call it from within
  Expression and Methylation. Any thesis claim of a live "Pathways" Cross-Omics tab must be
  treated as false unless the wiring is added.
```

```text
Function: cx_build_report(df, provenance)
File: crossomics_integration_plots.R:278-307
Purpose: Builds the Markdown "Analysis report" download (parameters, summary counts, top 20
  integrated genes by |log2FC|, and a fixed interpretation paragraph: "These results represent
  statistical associations ... They do not, on their own, establish a causal regulatory
  relationship").
Called from: mod_cross_integration.R's dl_report/dl_all download handlers.
Audit status: PASS.
```

```text
Function: cx_validation_checklist_ui(validation) / cx_sex_comparison_summary_ui(cmp) /
          cx_methodology_references_ui()
File: crossomics_integration_plots.R:313-379
Purpose: UI renderers for cx_validate_dataset()'s checklist, cx_compare_sexes()'s summary, and a
  static Methodology & References block.
Audit status: `cx_validation_checklist_ui` — same repository-wide search shows it is likewise not
  called from mod_cross_integration.R or any other module (WARNING, dead code, consistent with
  cx_validate_dataset() being computed at line 236 but never rendered anywhere with this
  renderer). `cx_sex_comparison_summary_ui` / `cx_methodology_references_ui` — WARNING, dead code,
  same basis as cx_compare_sexes() above.
```

### 3.5 `mod_cross_biomarker_conv.R` and `crossomics_biomarkerconv_helpers.R`

```text
Function: cx_bc_load_precomputed(sex)
File: crossomics_biomarkerconv_helpers.R:53-63
Purpose: Reads the pipeline's own already-joined per-sex eQTL-MR x mQTL-MR x DEG x DMP x DMR table
  (CX_RESULTS_DIR/cross_omics_eQTL_mQTL_{sex}.csv), then backfills a known data gap.
Called from: mod_cross_biomarker_conv.R (load_table handler); mod_cross_mr_stage.R (join_df()
  fallback path).
Audit status: PASS.
```

```text
Function: cx_bc_backfill_mqtl_from_mrstage(df)
File: crossomics_biomarkerconv_helpers.R:83-99
Purpose: Documents and fixes a specific, disclosed data-preparation defect: in the source CSV,
  in_eQTL_MR_panel and in_mQTL_MR_panel are mutually exclusive for every gene (0 overlap) because
  the file was built by stacking two separately-exported blocks rather than merging genes present
  in both. For genes not already flagged in_mQTL_MR_panel=TRUE that DO have a real mQTL instrument
  in the sibling mr_stage_eqtl_significant_genes_mqtl_mr.csv file, this fills in
  in_mQTL_MR_panel/mQTL_candidate_cpg/mQTL_MR_beta/mQTL_MR_pval/mQTL_instrument_available from that
  gene's best (lowest-p) CpG instrument. mQTL_cpg_chr/mQTL_cpg_pos_hg19 stay NA for backfilled rows
  (never fabricated, since the source that gap is filled from does not carry coordinates).
Statistical/computational principle: gap-filling from a sibling precomputed file, not a new
  statistical test; explicitly documented as the one exception to "never read outside this table."
Audit status: PASS — the defect and its fix are both disclosed in code comments and independently
  verifiable against the described sibling file; this is a documented data-preparation correction,
  not a silent alteration.
```

```text
Function: cx_bc_load_eqtl_upload(datapath, filename) / cx_bc_load_mqtl_upload(datapath, filename)
File: crossomics_biomarkerconv_helpers.R:121-163
Purpose: Parses an uploaded eQTL-MR or mQTL-MR results file (gene required for eQTL; gene +
  mQTL_MR_pval required for mQTL), numeric-coercing known columns.
Audit status: PASS.
```

```text
Function: cx_bc_merge_eqtl_mqtl(eqtl_df, mqtl_df)
File: crossomics_biomarkerconv_helpers.R:173-186
Purpose: Plain outer join of the (at most two) separately-uploaded layers by gene; adds
  DEG_adjP/DMP_fdr_bacon/DMR_fdr/mQTL_MR_pval as all-NA when the corresponding file/column is
  absent, so those layers read as "not evaluated," never fabricated.
Audit status: PASS.
```

```text
Function: cx_bc_dedup_min(df, key_col, order_col)
File: crossomics_biomarkerconv_helpers.R:191-199
Purpose: Generic "one row per key, smallest order_col wins" reducer, reused by the backfill
  function above and by Cross-Omics MR's own best-instrument-per-gene selection.
Audit status: PASS.
```

```text
Function: cx_bc_relabel(df, params=CX_BC_DEFAULT_PARAMS)
File: crossomics_biomarkerconv_helpers.R:221-246
Purpose: Recomputes DEG_significant/DMP_genomewide_significant/DMR_significant/
  mQTL_MR_significant/eQTL_MR_significant/methylation_significant/n_evidence_layers from the
  table's own already-retained raw p/FDR/effect-size columns, at whatever thresholds are passed
  (defaults: FDR<0.05 for DEG/DMP/DMR, nominal p<0.05 for mQTL-MR; eQTL-MR "significant" = "is in
  the panel at all," since that panel was already FDR<0.05-filtered upstream before this table was
  built). No join or new statistic is computed — pure relabeling of existing values.
Called from: mod_cross_biomarker_conv.R's bc_df(); mod_cross_mr_stage.R (uploaded-evidence path).
Audit status: PASS — genuinely a relabel, verified against the fact that no file I/O or merge
  occurs inside this function.
```

### 3.6 `mod_cross_mr_stage.R` and `crossomics_mrstage_helpers.R`

```text
Function: cx_mr_load_precomputed()
File: crossomics_mrstage_helpers.R:45-56
Purpose: Reads the pipeline's own precomputed single-instrument Wald-ratio MR results
  (CX_RESULTS_DIR/mr_stage_eqtl_significant_genes_mqtl_mr.csv; GoDMC cis-mQTL exposure ->
  Ishigaki et al. 2022 rheumatoid arthritis GWAS outcome, GCST90132223).
Audit status: PASS.
```

```text
Function: cx_mr_load_upload(datapath, filename) / cx_mr_load_evidence_upload(datapath, filename)
File: crossomics_mrstage_helpers.R:72-100, 163-188
Purpose: Parses an uploaded MR-instrument file (gene+pval required; FDR recomputed via BH from
  pval if not supplied) or an uploaded gene-level evidence file (only gene required; every
  DEG/DMP/DMR/mQTL/eQTL column optional, defaulting to NA — never fabricated — if omitted).
Audit status: PASS.
```

```text
Function: cx_mr_classify_categories(join_df)
File: crossomics_mrstage_helpers.R:197-221
Purpose: Checks every gene against 5 independent (not mutually exclusive) evidence combinations —
  DEG-DMP-QTL, DEG-DMR-QTL, DEG-eQTL, DMP-mQTL, DMR-mQTL — each a boolean AND of already-relabeled
  significance flags from join_df (Biomarker Convergence's own cx_bc_relabel() output).
Statistical/computational principle: pure boolean relabeling/filtering, no new statistical test;
  explicitly replaces an earlier Tier 1/2/3 mutually-exclusive priority ranking with these 5
  independent combinations (a gene can match more than one).
Called from: mod_cross_mr_stage.R's categories() reactive.
Audit status: PASS.
```

---

## 4. Special audit: gene-level methylation aggregation and the effect/significance-mismatch concern

Directly checking the scenario the audit brief describes (a gene-level effect size paired with a
significance statistic computed from a *different* CpG than the ones producing that effect):

- **`min_fdr` and `max_abs_dbeta` methods**: select one CpG's entire row — dbeta, pvalue, and fdr
  all come from that same CpG. No mismatch is possible for these two methods.
- **`mean`, `median`, and the three region-restricted variants (which reduce to `mean` once
  filtered)**: `dbeta_v` is the mean/median of every CpG's dbeta for that gene. The code's own
  comment (`crossomics_integration_helpers.R:311-322`) explicitly names the exact risk the audit
  brief describes ("its significance must be derived from that same set — not min(fdr)/
  min(pvalue), which can (and generically will) be driven by a single different CpG") and avoids
  it by computing a **directional Stouffer's Z combination across the same set of CpGs that
  produced dbeta_v**, then recomputing FDR via Benjamini-Hochberg across genes. The reported
  effect and the reported significance for a gene therefore come from the same underlying CpG set
  for these methods.

**Conclusion**: the specific failure mode described in the audit brief (mean effect reported
alongside a min-FDR-of-a-different-CpG significance value) is **not present** in the inspected
code for any of the 7 aggregation methods offered. This is recorded as a **PASS**, with the
qualification that the `min_fdr`/`max_abs_dbeta` methods still discard the other CpGs' effect
sizes from the gene-level `dbeta` (by design — they report one representative CpG, not a
combination), which is a different and disclosed trade-off, not a mismatch.

---

## 5. Overall audit result

### Functionality — PASS
Every function traced actually does what its name and surrounding comments say it does; no
invented or missing logic was found in the Dataset tab or its three downstream sub-modules.

### Data integrity — WARNING
No plausible-value-range check exists for an uploaded methylation column (beta/Δβ), and no
duplicate-ID check exists at standardization time; both are silent gaps rather than failures, but
both are real (see Dataset Tab Audit Section 14).

### Statistical integrity — PASS
The one concern the audit brief specifically asked to check (mean effect vs. different-CpG
significance) does not occur; Stouffer's Z combination correctly ties the significance statistic
to the same CpG set as the reported effect for the mean/median aggregation methods.

### Methylomics preprocessing compatibility — WARNING
The Dataset tab performs no preprocessing of its own and does not verify that an uploaded
methylation column is on a comparable scale to the "Example data" beta-scale Δβ; upstream QC for
the bundled data is real and documented (`METHODS_load_qc.md` etc.) but its executable source is
outside this repository, so it cannot itself be re-verified from the inspected code.

### Transcriptomics–methylomics matching — WARNING
Gene-level matching is exact-text plus an explicit, disclosed, non-fuzzy ID-harmonization pass
(PASS in isolation); sample-level matching is conservative and correctly gated
(`cx_detect_sample_pairing()`, PASS), but is structurally unreachable for this app's own bundled
"Example data" (both `expr_wide`/`meth_wide` are always `NULL` for it) — a real limitation worth
stating plainly in any thesis text that discusses sample-level correlation.

### Reproducibility — PASS
Every Integration run's exact parameters (thresholds, aggregation method, correlation method,
sample-matching status, annotation source, first-gene-only caveat, timestamp) are captured by
`cx_build_provenance()` and included in every downloadable output.

### Thesis documentation readiness — PASS
Every scientific statement in this pair of documents is traced to a specific file and line range;
no step is described as implemented unless it was directly located in the source.
