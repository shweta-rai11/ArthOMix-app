# Cross-Omics Dataset Tab — Code Audit

## 1. Scope

This document audits **only** the path: **Methylomics → Cross-Omics → Dataset tab**, i.e. the file
`ArthOMix/R/crossomics/mod_cross_dataset.R` and the pure helper functions it directly calls
(`crossomics_integration_upload.R`, and the loading/standardization functions in
`crossomics_integration_helpers.R`). Downstream sub-modules (Expression and Methylation,
Biomarker Convergence, Cross-Omics MR) are described only to the extent needed to show what the
Dataset tab hands off to them — their own internal analyses are documented in the companion file
`Cross-Omics_End_to_End_Function_Documentation.md`, not audited here.

No UI, styling, navigation, other Methylomics sub-modules, Transcriptomics, Multi-Omics, or any
other Cross-Omics tab was modified to produce this document. This is a read-only inspection.

All statements below are traced to source lines in the files listed in Section 2. Where the
inspected code does not perform an operation a methylomics pipeline might normally be expected to
perform, this document says so explicitly rather than assuming it happened.

---

## 2. Exact files inspected

| File | Role |
|---|---|
| `ArthOMix/R/crossomics/mod_cross_dataset.R` | The Dataset tab itself — UI (`mod_cross_dataset_ui`) and server (`mod_cross_dataset_server`) |
| `ArthOMix/R/crossomics/crossomics_integration_upload.R` | Upload file reading (`cx_read_table`) and auto-detection dispatch (`cx_read_and_detect`) for "Upload your own data" |
| `ArthOMix/R/crossomics/crossomics_integration_helpers.R` | Pure, Shiny-free logic: column auto-detection, standardization, gene-level methylation aggregation, sample-pairing, CpG→gene annotation, preloaded-data loaders, gene-ID harmonization, dataset validation |
| `ArthOMix/ui.R` (lines ~1460–1521, `crossomicsUI()`) | Where the Dataset tab is mounted as the first tab of the Cross-Omics module |
| `ArthOMix/server.R` (lines ~97–130) | Where `cross_dataset` (the shared store the Dataset tab writes to) is created and wired to `mod_cross_dataset_server` and to every `CX_MODULES` sub-module server |
| `ArthOMix/R/submodules_registry.R` (lines 57–71) | Confirms `CX_MODULES` (Expression and Methylation, Biomarker Convergence, Cross-Omics MR) is a separate registry from `mod_cross_dataset.R`, and does not include the Dataset tab itself |
| `ArthOMix/data_paths.R` (lines 84–91, 119–131) | `METH_DATA_ROOT`, `METH_DATA_AVAILABLE`, `METH_DMR_DIR`, `CX_DATA_ROOT`, `CX_DATA_AVAILABLE`, `CX_RESULTS_DIR`, `DATA_ROOT` — the path constants the Dataset tab's loaders resolve against |
| `ArthOMix/global.R` (lines 296–332) | `METH_QC_PROBE_CASCADE`, `load_default_dmp()`, `load_default_dmr()` — the functions `cx_load_default_methylation()`/`cx_load_default_dmr()` call to reach the actual preloaded methylation files |
| `ArthOMix/R/crossomics/mod_cross_integration.R` (lines 1–41, 110–155) | The one consumer of the Dataset tab's output (`cross_dataset$user_expr_*` / `user_meth_*`) |
| `ArthOMix/data/preloaded/methylomics/tables/script01_dataload_QC/METHODS_load_qc.md` | Documented provenance of the upstream methylation QC/probe-filtering that produced the beta-value matrix underlying the preloaded DMP/DMR tables |
| `ArthOMix/data/preloaded/methylomics/tables/script03_dmp_sva_sexstratified/METHODS_dmp_sva_sexstratified.md` | Documented provenance of the SVA-adjusted, bacon-corrected DMP statistics the Dataset tab's "Example data" methylation source reads |
| `ArthOMix/data/preloaded/methylomics/tables/script04_dmr_sexstratified/METHODS_dmr_sexstratified.md` | Documented provenance of the DMRcate region calls the Dataset tab's DMR option reads |
| `ArthOMix/tests/testthat/test-data-loaders.R` | Confirms what test coverage exists (or does not exist) for this code path |

No filenames were assumed. Every function named in this document was located by direct text search
(`grep`) or by reading the file in full.

---

## 3. Dataset tab architecture (UI → server → reactive state)

**UI (`mod_cross_dataset_ui`, `mod_cross_dataset.R:48-96`).** A single-page layout — no internal
tabs — with two panels:
- **"1. Data Source"** (left, `box`): a `radioButtons` (`source_mode`) choosing "Example data" or
  "Upload your own data". Conditional on the choice, either (a) two more `radioButtons`
  (`sex_stratum`: ALL/FEMALE/MALE; `meth_level`: CpG-level DMP or Region-level DMR) plus a
  "Load example data" `actionButton`, or (b) two `fileInput`s (Transcriptomics file, Methylomics
  file, both optional, accepting `.csv/.tsv/.txt/.xlsx`). Below both branches: "Use this data" and
  "Clear" `actionButton`s.
- **"Preview"** (right, `box`): a `uiOutput("preview_ui")` that renders whatever has been loaded —
  identical rendering code path for Example data and for Upload, so the preview genuinely doubles
  as a worked example of the shape an upload should take.

**Server (`mod_cross_dataset_server`, `mod_cross_dataset.R:98-263`).** Two module-local
`reactiveVal`s, `expr_data` and `meth_data`, each holding `NULL` or
`list(df, source, raw, mapping)` — one shared shape regardless of which source mode produced it.
Four `observeEvent` blocks populate them (`load_example_btn`, `expr_file`, `meth_file`, plus a
`source_mode`-change observer that clears both). One `observeEvent(input$use_data_btn)` copies
`expr_data()`/`meth_data()` into the module-external `cross_dataset` `reactiveValues` object passed
in by `server.R`. One `observeEvent(input$clear_btn)` clears both the local reactiveVals and every
field in `cross_dataset`. The rest of the server (`preview_ui`, `expr_table`, `expr_plot`,
`meth_table`, `meth_plot`) is read-only rendering of `expr_data()`/`meth_data()`.

**Reactive state.** `cross_dataset` is created once in `server.R:104-107` as
`reactiveValues(user_expr_df, user_expr_source, user_expr_wide, user_expr_mapping,
user_expr_sample_cols, user_meth_df, user_meth_source, user_meth_wide, user_meth_mapping,
user_meth_sample_cols)`, all `NULL`/`character(0)` initially. It is passed by reference into
`mod_cross_dataset_server("cx_dataset", cross_dataset)` and into every `CX_MODULES` sub-module's
server call (`server.R:122-130`). Only the Dataset tab ever writes to it; every sub-module only
reads it. `mod_cross_integration.R:142-154` mirrors it into its own `raw` `reactiveValues` inside a
plain `observe()` block with no debounce/throttle — any change on the Dataset tab (a new
"Load example data" click, a new file upload, "Use this data", or "Clear") propagates to the
Expression and Methylation tab on the very next reactive flush, with no re-upload step required
there.

Confirmed by `submodules_registry.R:65-69`: `CX_MODULES` lists only `integration`,
`biomarkerconv`, and `mrstage` — **not** `dataset`. The Dataset tab is mounted directly in
`ui.R`'s `crossomicsUI()` (`ui.R:1503`) as the first `tabPanel` of the top-level
`cx_menu` `tabsetPanel`, structurally separate from the `CX_MODULES` sub-module grid that lives
under the "Sub-modules" tab (`ui.R:1515`).

---

## 4. Every input

| Input ID | Type | Choices / accepted formats | Effect |
|---|---|---|---|
| `source_mode` | `radioButtons` | `example` / `upload` | Switches which panel is shown; changing it clears both `expr_data`/`meth_data` (`mod_cross_dataset.R:173`) |
| `sex_stratum` | `radioButtons` (example mode only) | `all` / `female` / `male` | Passed to `cx_load_default_deg(sex=)` and `cx_load_default_methylation(sex=)`/`cx_load_default_dmr(sex=)` |
| `meth_level` | `radioButtons` (example mode only) | `dmp` / `dmr` | Chooses `cx_load_default_methylation()` (CpG-level) vs. `cx_load_default_dmr()` (region-level) |
| `load_example_btn` | `actionButton` | — | Triggers the example-data load `observeEvent` |
| `expr_file` | `fileInput` (upload mode only) | `.csv/.tsv/.txt/.xlsx` | Triggers `cx_read_and_detect(kind="expression")` → `cx_standardize_expression()` on file selection |
| `meth_file` | `fileInput` (upload mode only) | `.csv/.tsv/.txt/.xlsx` | Triggers `cx_read_and_detect(kind="methylation")` → `cx_standardize_methylation()` on file selection |
| `use_data_btn` | `actionButton` | — | Publishes `expr_data()`/`meth_data()` into `cross_dataset` |
| `clear_btn` | `actionButton` | — | Clears local state and `cross_dataset` |

There is **no manual column-mapping UI** on this tab. If auto-detection fails, the tab shows an
error/warning notification and explicitly directs the user to "Expression and Methylation's own
Upload option, which supports manual column mapping" (`mod_cross_dataset.R:146,161`) — it never
guesses at an ambiguous column itself.

---

## 5. Every output

| Output ID | Renderer | Content |
|---|---|---|
| `preview_ui` | `renderUI` | Header text + table + plot for whichever of `expr_data()`/`meth_data()` is non-`NULL`; an `empty-note` prompt otherwise |
| `expr_table` | `DT::renderDataTable` | `expr_data()$df` (gene, log2fc, pvalue, fdr) |
| `expr_plot` | `renderPlot` | Histogram of `log2fc` (ggplot2, 40 bins) |
| `meth_table` | `DT::renderDataTable` | `meth_data()$df` (cpg, gene, dbeta, pvalue, fdr, chr, pos, region, region_fine, island_context) |
| `meth_plot` | `renderPlot` | Histogram of `dbeta` (ggplot2, 40 bins) |

There is no download handler on the Dataset tab itself — data leaves this tab only via
`cross_dataset` for the downstream sub-modules (which have their own CSV/TSV/XLSX/report
downloads, documented in the End-to-End file).

---

## 6. Every function (called directly by the Dataset tab)

| Function | File | Called by | Purpose |
|---|---|---|---|
| `cx_load_default_deg()` | `crossomics_integration_helpers.R:646-651` | `load_example_btn` handler | Reads `DEG_{sex}_full.csv` |
| `cx_load_default_methylation()` | `crossomics_integration_helpers.R:659-699` | `load_example_btn` handler (when `meth_level="dmp"`) | Reads the preloaded SVA/bacon-corrected DMP table and annotates it to genes |
| `cx_load_default_dmr()` | `crossomics_integration_helpers.R:713-738` | `load_example_btn` handler (when `meth_level="dmr"`) | Reads the preloaded DMRcate region table |
| `cx_standardize_expression()` | `crossomics_integration_helpers.R:96-110` | example-data handler and `expr_file` handler | Builds the standardized `gene, log2fc, pvalue, fdr` table |
| `cx_standardize_methylation()` | `crossomics_integration_helpers.R:165-191` | example-data handler and `meth_file` handler | Builds the standardized `cpg, gene, dbeta, pvalue, fdr, chr, pos, region*, island_context` table |
| `cx_read_and_detect()` | `crossomics_integration_upload.R:40-46` | `expr_file`/`meth_file` handlers | Reads an uploaded file and auto-detects its column mapping |
| `cx_read_table()` | `crossomics_integration_upload.R:13-36` | `cx_read_and_detect()` | Parses CSV/TSV/TXT (`data.table::fread`) or XLSX (`openxlsx::read.xlsx`) into a data.frame |
| `cx_detect_columns()` | `crossomics_integration_helpers.R:55-67` | `cx_read_and_detect()` | Regex-matches column names to canonical fields (`CX_FIELD_PATTERNS`) |
| `cx_detect_sample_columns()` | `crossomics_integration_helpers.R:79-85` | `use_data_btn` handler | Flags numeric columns not already claimed by a mapped field as candidate per-sample columns |

Full plain-English descriptions of each are in Section 7 below and, for the complete
input/structure/output/QC/audit template, in `Cross-Omics_End_to_End_Function_Documentation.md`.

---

## 7. Dataset loading — exact logic

### "Example data" path

1. `cx_load_default_deg(sex)` (`crossomics_integration_helpers.R:646-651`): builds
   `file.path(CX_DEG_TABLE_DIR, sprintf("DEG_%s_full.csv", sex))` where
   `CX_DEG_TABLE_DIR = file.path(DATA_ROOT, "results", "tables")`. Returns `NULL` (not an error) if
   the file does not exist — "no file" and "wrong file" are not distinguished beyond that.
   `DATA_ROOT` is defined in `data_paths.R:45` as the Transcriptomics preloaded path, so this is
   reading a Transcriptomics DEG table, not a Methylomics file — the Cross-Omics module's
   Transcriptomics side of "Example data" always comes from the Transcriptomics pipeline's own
   output.
2. For methylation, if `meth_level == "dmr"`: `cx_load_default_dmr(sex)`
   (`crossomics_integration_helpers.R:713-738`) reads `METH_DMR_DIR/dmr_{sex}_full.csv`
   directly via `data.table::fread`, builds a synthetic `cpg` id as
   `"{seqnames}:{start}-{end}"`, and maps `dbeta = meandiff`, `pvalue = Stouffer`,
   `fdr = dmr_fdr` — this DMR table already carries its own `overlapping.genes` column, so no
   separate annotation lookup is performed for the DMR path.
   Otherwise (`meth_level == "dmp"`, the default): `cx_load_default_methylation(sex)`
   (`crossomics_integration_helpers.R:659-699`):
   - Guards on `METH_DATA_AVAILABLE` (a directory-existence flag from `data_paths.R:85`).
   - Calls `load_default_dmp(stage = "sva", sex = sex)` (`global.R:313-320`), which reads
     `METH_DMP_SVA_DIR/dmp_{sex}_full.csv`. **Only the `"sva"` (surrogate-variable-adjusted,
     bacon-corrected) stage is ever read here — the `"plain"` (unadjusted) stage exists in the
     data directory but this tab never reads it.**
   - This DMP table has no gene column (only `cpg`, `dbeta`, and the SVA/bacon-corrected
     `p_bacon`/`fdr_bacon` statistics) — `cx_get_region_annotation("450K")` is called to attach
     `gene`/`chr`/`pos`/`region_raw`/`island_context` per CpG (Section 9 below).
   - Rows whose annotation lookup did not resolve to a non-empty gene symbol are dropped
     (`crossomics_integration_helpers.R:691`).
3. Both results are passed to `cx_standardize_expression()`/`cx_standardize_methylation()`, exactly
   the same standardization functions the upload path uses (Section 8 below) — "Example data" is
   not a structurally different object, it is the same standardized shape produced by a different
   loader.
4. On success, `expr_data()`/`meth_data()` are set with `raw = NULL, mapping = NULL` — deliberately,
   since the preloaded tables are already gene/CpG-level with no per-sample columns to detect
   (`mod_cross_dataset.R:102-106`). This means **sample-level correlation is structurally
   unavailable for "Example data"** — confirmed in Section 11.

### "Upload your own data" path

1. `cx_read_and_detect(datapath, filename, kind)` → `cx_read_table()` parses the file
   (`fread`/`openxlsx::read.xlsx`); fails soft (`list(ok=FALSE, error=...)`) on an unparsable file,
   fewer than 2 columns, zero rows, or an unsupported extension.
2. `cx_detect_columns(df, kind)` regex-matches each canonical field (`gene`, `log2fc`/`dbeta`/
   `beta`, `pvalue`, `fdr`, `chr`, `pos`, `region`, `island`, `cpg` for methylation) against the
   file's column names, greedily claiming the first match per field and excluding already-claimed
   columns from later fields (Section 9's audit covers the false-positive risk this has).
3. `cx_standardize_expression()`/`cx_standardize_methylation()` build the same standardized shape as
   the example path. On failure (no gene column, no log2FC/Δβ column, or zero valid rows after
   dropping blank/`NA` genes) a `showNotification` explicitly tells the user to use Expression and
   Methylation's own manual-mapping Upload option instead — it never falls back to a guess.
4. `expr_data()$raw`/`meth_data()$raw` retain the **original wide upload**, and
   `cx_detect_sample_columns()` is run **only at "Use this data" time**
   (`mod_cross_dataset.R:184,191`), not at upload time — every numeric column not already claimed
   by a mapped metadata field, and not in the hardcoded `CX_NON_SAMPLE_NUMERIC_NAMES` exclusion
   list (`aveexpr`, `t`, `n`, `meandiff`, etc.), is treated as a per-sample column. This is the
   only path by which `user_expr_sample_cols`/`user_meth_sample_cols` become non-empty, which is
   in turn the only way sample-level correlation becomes available downstream (Section 11).

---

## 8. Data validation

The Dataset tab's own validation is narrow and happens inside the standardization functions, not as
a separate step:

- `cx_standardize_expression()`: requires a detected `gene` mapping and a detected `log2fc`
  mapping; rows with blank/`NA` gene are dropped; requires ≥1 remaining row.
- `cx_standardize_methylation()`: requires a detected `gene` mapping (methylation is always joined
  to expression by gene, so this is described in-code as "required to join with expression") and
  requires **either** a `dbeta` **or** a `beta` mapping; same blank/`NA`-gene drop and
  ≥1-row requirement.
- `cx_dedup_by_gene()` (called inside `cx_standardize_expression()`): if the same gene symbol
  appears more than once (e.g. multiple probes/transcripts), keeps only the row with the smallest
  FDR (falling back to smallest p-value, then the first row) — deterministic, but this is the
  **expression**-side dedup; the methylation side is deliberately **not** deduplicated at
  standardization time (multiple CpGs per gene are preserved and only collapsed later, at
  Integration-run time, by `cx_aggregate_methylation()` — see Section 9).

The **richer, plain-language checklist** the task's spec anticipates
(`cx_validate_dataset()`, `crossomics_integration_helpers.R:970-995`) exists in the codebase but is
**not called anywhere in `mod_cross_dataset.R`** — it is only invoked from
`mod_cross_integration.R:236`, after "Run Integration" on the Expression and Methylation tab, where
its result is stored as `integ$validation`. A repository-wide search additionally shows that
`integ$validation` is never read anywhere else in `mod_cross_integration.R`, and that its own UI
renderer, `cx_validation_checklist_ui()` (`crossomics_integration_plots.R:313`), is not called from
any file in the repository. So on the Dataset tab itself, "ready" is communicated only by whether
the Preview renders and whether a `showNotification` fired; the gene-identifier-detected /
log2FC-detected / FDR-detected / overlap checklist is **computed but never displayed anywhere in
the current UI** — a stronger finding than "one tab downstream," recorded in Section 14.

---

## 9. QC

**No probe-level or sample-level quality control is performed by the Dataset tab or by any function
it calls**, for either source mode:

- "Upload your own data": `cx_read_table()`/`cx_standardize_methylation()` perform only structural
  checks (parseable file, required columns present, numeric coercion, non-blank gene). There is no
  check for: detection p-value, bead count, SNP-associated probes, cross-reactive probes, sex
  chromosomes, missingness threshold, or plausible beta/Δβ range. **Not implemented in the
  inspected Dataset-tab code.**
- "Example data": the QC that produced the underlying `dmp_{sex}_full.csv`/`dmr_{sex}_full.csv`
  files happened **upstream, outside this application**, in an external pipeline whose *methods
  documentation* (not its executable source) is bundled in the repository at
  `data/preloaded/methylomics/tables/script01_dataload_QC/METHODS_load_qc.md`,
  `script03_dmp_sva_sexstratified/METHODS_dmp_sva_sexstratified.md`, and
  `script04_dmr_sexstratified/METHODS_dmr_sexstratified.md`. Per those documents (traced, not
  invented): a sequential 5-step probe-filtering cascade (`cg`-prefix restriction → Zhou et al.
  2017 `MASK_general` removal → multi-hit removal → sex-chromosome removal → >5% missingness
  removal, retaining 412,492 of 485,577 probes — the exact counts are also duplicated in this
  app's own `global.R:300-304` as `METH_QC_PROBE_CASCADE`), PCA/MAD-based sample-outlier detection
  within sex strata, and a chromosome-Y-methylation sex-mismatch check. **The executable R scripts
  that performed this QC are not present anywhere in this repository** — only their methods
  write-up and their final output tables are. The Dataset tab (and the rest of the Cross-Omics
  module) never re-runs, re-derives, or re-checks any of it; it treats the preloaded DMP/DMR files
  as already-QC'd, finished results.
- The live, interactive methylation QC tooling that **does** exist in this codebase
  (`R/methylomics/qc.R` — `methyl_filter_detection_p()`, `methyl_filter_beadcount()`,
  `methyl_filter_snp()`, `methyl_filter_sex_chr()`, `methyl_sex_check()`, etc., used by
  `mod_methyl_qc.R`) is **never called anywhere in `R/crossomics/`**. It belongs to the
  Methylomics module's own interactive QC tab, which operates on data a user uploads there, and is
  architecturally unconnected to the Cross-Omics Dataset tab.

---

## 10. Normalization

**No normalization or scale transformation is performed by the Dataset tab.**

- Expression: `log2fc` is passed through from whatever column was auto-detected/mapped as the
  fold-change column, with only `as.numeric()` coercion (`cx_as_numeric_safe()`,
  `crossomics_integration_helpers.R:91`). No re-normalization, re-centering, or re-scaling.
- Methylation: `dbeta`/`beta` likewise pass through `as.numeric()` only. Whether the uploaded
  values are on a beta scale, an M-value scale, or something else is **not checked or enforced by
  the code** — the field is labelled "Δβ (methylation change)" in the UI and downstream plots, but
  nothing validates that an uploaded numeric column is actually bounded in a beta-like range. **Not
  implemented in the inspected code**; this is flagged as an audit finding in Section 14.
- "Example data" methylation values are exactly the `dbeta` (beta-scale mean difference, RA vs.
  Control) and `p_bacon`/`fdr_bacon` columns already computed upstream by the external SVA/bacon
  pipeline (per `METHODS_dmp_sva_sexstratified.md` Section 2.AA.4, that pipeline itself modelled
  **M-values** with `limma::lmFit()`/`eBayes()` for the test statistic, while `dbeta` is a
  separately retained beta-scale descriptive effect size) — the Dataset tab reads this finished
  `dbeta` value as-is; it does not itself convert between beta and M-value at any point.
- No double-normalization is possible in this specific tab because no normalization step exists
  here at all to duplicate.

---

## 11. Sample matching

The Dataset tab **detects** which columns look like per-sample columns
(`cx_detect_sample_columns()`) but **does not perform any transcriptomics↔methylomics sample
matching itself** — that computation (`cx_detect_sample_pairing()`,
`crossomics_integration_helpers.R:550-557`) lives in, and is only invoked from,
`mod_cross_integration.R:221` at "Run Integration" time, one tab downstream.

What the Dataset tab does establish, and hands off via `cross_dataset`:

- `user_expr_sample_cols` / `user_meth_sample_cols`: character vectors of candidate per-sample
  column names, computed **only for the Upload path** (`raw` is `NULL` for Example data, so
  `cx_detect_sample_columns()` is never called for it — `mod_cross_dataset.R:184,191`). This means,
  as a direct consequence of the code: **"Example data" can never produce a sample-level-paired
  Integration run** — its `expr_wide`/`meth_wide` are always `NULL`, so
  `cx_detect_sample_pairing()` downstream always sees `n_expr = n_meth = 0` and reports
  `paired = FALSE` for it, regardless of sex stratum or thresholds.
- The identifier used for matching (once it does happen downstream) is the **column name itself**
  (`intersect(expr_samples, meth_samples)`, `crossomics_integration_helpers.R:551`) — no prefix/
  suffix stripping, no case-folding, no fuzzy matching is applied anywhere in this path.
- No metadata (sex, disease status, age) is joined at sample level anywhere in this module; the
  Dataset tab and Expression and Methylation module operate on gene/CpG-level summary statistics
  (log2FC, Δβ, p/FDR), not on subject-level clinical covariates.

---

## 12. Feature matching (gene ↔ CpG)

This is annotated but **not aggregated** on the Dataset tab itself — aggregation is a Run
Integration-time step (Section 9 of the End-to-End document covers `cx_aggregate_methylation()` in
full, including the audit of Section 27's mean-effect/min-FDR mismatch concern). What the Dataset
tab does is:

- **Upload path**: whatever `gene`/`cpg` columns the uploader supplied are used verbatim — no
  genomic coordinate lookup is performed. If the file has no `cpg` column, a synthetic row index
  id (`"row1"`, `"row2"`, …) is substituted (`crossomics_integration_helpers.R:175`) purely so the
  `cpg` field is never `NA`; this is not a real probe identifier.
- **Example data (DMP)**: `cx_get_region_annotation("450K")` (`crossomics_integration_helpers.R:
  604-633`) reads the `Locations`, `Other`, and `Islands.UCSC` objects directly from the
  `IlluminaHumanMethylation450kanno.ilmn12.hg19` Bioconductor package (not via
  `R/methylomics/annotation.R::methyl_get_annotation()`, which is a **separate, independent**
  reader of the same underlying package — deliberately, per that file's own comment, because
  `methyl_get_annotation()` does not expose `UCSC_RefGene_Group`, which this module needs for
  promoter/gene-body classification). For a CpG whose `UCSC_RefGene_Name` lists more than one
  gene (semicolon-separated, a real and common occurrence on the 450K array), **only the
  first-listed gene is used** (`strsplit(...)[[1]]`, `crossomics_integration_helpers.R:624-625`)
  — every other co-annotated gene for that probe is silently not represented in this table. This
  exact behavior is also disclosed to the user at Integration-run time via
  `cx_build_provenance()`'s "Note: CpG probes annotated to more than one gene ... are assigned to
  their first-listed gene only" (`crossomics_integration_helpers.R:761-762`).
- **Example data (DMR)**: uses the DMR table's own pre-computed `overlapping.genes` column
  directly (`crossomics_integration_helpers.R:724`) — no separate annotation lookup, and no
  "first gene only" truncation is applied here (a region can genuinely list multiple overlapping
  genes, retained as supplied).

No multi-CpG-per-gene **aggregation** happens on the Dataset tab — every CpG row from the
annotation lookup remains a separate row in `meth_data()$df`. Aggregation to one row per gene is
strictly a Run Integration-time operation performed downstream by
`cx_aggregate_methylation()`.

---

## 13. Data passed downstream

`cross_dataset$user_expr_df` / `user_meth_df` (the standardized data.frames), `user_expr_source` /
`user_meth_source` (a human-readable provenance string, e.g. `"Example data (FEMALE,
sex-stratified DEG)"` or `"Uploaded: mydata.csv"`), `user_expr_wide` / `user_meth_wide` (the
original wide upload, or `NULL` for Example data), `user_expr_mapping` / `user_meth_mapping` (the
detected column mapping, or `NULL` for Example data), and `user_expr_sample_cols` /
`user_meth_sample_cols` (candidate sample-column names, always `character(0)` for Example data).

The **only** downstream consumer is `mod_cross_integration.R` ("Expression and Methylation"),
mirrored live via a plain `observe()` block (`mod_cross_integration.R:142-154`). Biomarker
Convergence and Cross-Omics MR read `cross_dataset` in their function signatures but, per direct
inspection, **never reference `cross_dataset` in their server bodies** — both load their own,
independent precomputed/uploaded data (Section 14 finding).

---

## 14. Audit findings

**Confirmed correct:**
- The "Example data" and "Upload your own data" paths genuinely converge on one standardized shape
  before reaching `cross_dataset` — verified by tracing both branches to the same
  `cx_standardize_expression()`/`cx_standardize_methylation()` calls.
- Auto-detection failures are surfaced to the user with an explicit, specific message and a
  pointer to the manual-mapping alternative, rather than silently guessing or silently producing
  an empty/wrong result.
- `source_mode` switching clears stale state from the other mode (`mod_cross_dataset.R:173`),
  preventing an old upload from being silently carried into a "Use this data" click after switching
  to Example data (or vice versa).
- The vestigial `CX_TABLE_REGISTRY`/`load_default_cx_table()` biomarker-convergence table browser
  that an earlier version of this tab exposed has been fully removed from `mod_cross_dataset.R`
  (confirmed absent from the current file); the header comment says so and the only remaining
  consumers are `tests/test-data-loaders.R` and the unrelated Multi-Omics module — consistent with
  what the code shows.

**Potentially problematic (documented, not silently fixed):**
- **No beta/Δβ range check on upload.** An uploaded "methylation" file's numeric column is accepted
  as `dbeta`/`beta` with no check that its values are plausible for a methylation measurement
  (e.g. within [-1, 1] for Δβ, or [0, 1] for beta). A user could upload an M-value column, a raw
  intensity column, or an unrelated numeric column under a header the regex happens to match, and
  the tab would standardize and preview it without complaint.
- **Sample-level correlation is structurally unreachable for "Example data".** This is
  intentional and disclosed in code comments, not a bug, but it means any thesis claim of
  "sample-level methylation-expression correlation" can only ever apply to the Upload path with
  genuinely wide, per-sample files on both sides — never to this application's own bundled example
  data.
- **`cx_validate_dataset()`'s readiness checklist is computed but never displayed anywhere in the
  current UI.** It is not surfaced on the Dataset tab itself, and — per the repository-wide search
  above — its result (`integ$validation`) is also never read by any output on the Expression and
  Methylation tab that computes it, and its dedicated renderer (`cx_validation_checklist_ui()`) is
  called from nowhere. A user has no way to see this checklist in the running application today.

**Unsupported assumptions a reader should not make:**
- That the Dataset tab performs any of the QC/normalization steps a methylomics pipeline typically
  performs. It does not; those steps, for the bundled data, happened in an external pipeline
  documented only in the `METHODS_*.md` files, whose executable source is not part of this
  repository.
- That "Region-level (DMR)" and "CpG-level (DMP)" data are interchangeable inputs to the same
  downstream statistics. `cx_load_default_methylation()`'s `pvalue`/`fdr` are per-CpG bacon-
  corrected values; `cx_load_default_dmr()`'s `pvalue`/`fdr` are the DMRcate Stouffer statistic and
  its own region-level BH correction (Section 6 above and `METHODS_dmr_sexstratified.md` Section
  2.BB.4) — different statistical objects, both surfaced through the same `pvalue`/`fdr` column
  names.

**Missing safeguards:**
- No file-size or row-count ceiling is enforced on upload before parsing (`cx_read_table()` reads
  the entire file into memory via `fread`/`openxlsx::read.xlsx` unconditionally).
- No duplicate-CpG-ID check is performed at standardization time on the methylation upload path
  (duplicates are tolerated and only implicitly resolved much later, at aggregation time, by
  whichever `agg_method` is chosen on the Integration tab).

---

## 15. Scientific interpretation

The Dataset tab establishes, and only establishes, that a Transcriptomics differential-expression
table and/or a Methylomics differential-methylation table are present, are parseable, and have been
mapped onto a common `gene`/`log2fc`/`pvalue`/`fdr` (expression) or
`cpg`/`gene`/`dbeta`/`pvalue`/`fdr` (methylation) schema. It does **not** establish, and makes no
claim to establish: statistical significance (that is computed downstream in Expression and
Methylation, at user-set thresholds), sample-level pairing (available only for genuinely wide
uploads, never for the bundled example data), or any biological/regulatory relationship between a
gene's expression and its methylation (that is explicitly framed downstream as "potential",
"association", never causal — see `crossomics_integration_helpers.R:370-390` and
`crossomics_integration_plots.R`'s report text). The Dataset tab is correctly scoped as an
entry-point/staging step, not an analysis step.
