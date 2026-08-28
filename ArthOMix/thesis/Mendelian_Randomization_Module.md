# Mendelian Randomization Module — `mod_mr.R`

**Source file:** `ArthOMix/R/transcriptomics/mod_mr.R` (1,055 lines)
**Shared helper functions used by this module (defined in `global.R`):** `load_mr_instrument_table()`, `estimate_mr_set()`, `gwas_col_map_ui()`, `guess_gwas_col()`, `GWAS_COL_PATTERNS`, `read_uploaded_table()`, `read_table_safe()`, `ARTHOMIX_COLORS`, `ARTHOMIX_STATUS`, `theme_arthomix()`.
Prepared: 2026-08-25

This document is derived **exclusively** from the code in `mod_mr.R` and the `global.R` functions it explicitly calls. Anything not present in this code is marked *"Not implemented in the current module code."*

---

## 1. Module Purpose

According to the module's own header comment (`mod_mr.R:1-39`), this is:

> "A general-purpose two-sample MR tool. 'Your analysis' runs a live MR test for any gene against a bundled RA dataset, or against your own uploaded exposure/outcome GWAS summary statistics for any trait."

The module registration object (`mod_mr_config`, lines 46–51) states the same scope:

```r
mod_mr_config <- list(
  id = "mr", group = "Genetics",
  title = "Mendelian Randomization",
  description = "A general two-sample Mendelian randomisation tool: test any gene against a
  bundled RA dataset, or upload your own exposure/outcome GWAS summary statistics for any trait.
  Adjustable instrument filters, LD clumping, MR-PRESSO outlier testing, an MHC sensitivity
  analysis, and a bundled RA reference result included as a worked example.",
  icon = "route"
)
```

### Scientific/analytical task performed

The module performs **two-sample Mendelian randomization** (MR), estimating the causal effect of a gene's expression (the exposure, via cis-eQTL instruments) on rheumatoid arthritis risk (the outcome, via a GWAS), or, in upload mode, between any user-supplied exposure and outcome GWAS summary-statistics pair.

Two distinct data layers exist side by side, as stated in the header comment (lines 20–31):

1. **"Your analysis"** (top of page) — runs the estimator hierarchy **live**, gene-by-gene or as a batch screen, against either:
   - a bundled, pre-cached cis-eQTL instrument set (harmonised to the Okada 2014 RA GWAS), or
   - the user's own uploaded exposure/outcome files.
2. **The bundled RA reference results** ("MR female" / "MR male" tabs) — read directly from pre-computed CSVs (`MR_MHC_sensitivity_{sex}.csv`), already FDR-corrected and MHC-sensitivity-tested. **Nothing in this layer is recomputed by the module.**

### Bundled dataset provenance (as documented in code, lines 10–18)

- File: `MR_primary_objects.rds` (path resolved via `MR_PRIMARY_OBJECTS_RDS`, defined in `data_paths.R`).
- Cis-restricted (±1 Mb of the gene body, GRCh37), MHC-flagged, and harmonised (`action = 2`) cis-eQTL instruments.
- Outcome: Okada 2014 RA GWAS (OpenGWAS `ieu-a-832`).
- Exposure: eQTLGen whole-blood cis-eQTL (Vösa et al. 2021) — stated in the UI text at line 99.
- Covers 1,701 genes with at least one usable (`mr_keep`) harmonised instrument.
- The advanced filters (p-value, F-statistic, MHC handling) subset this **already cis-filtered, cached object live**; the code explicitly notes it does *not* re-derive cis restriction from scratch, since that would require a GRCh37 annotation the app does not otherwise load.

### Data-integrity correction (`relabel_check`, lines 32–39)

A documented data-quality check found that `$dat$gene` in the cached RDS is stale for 21 SNPs across 37 genes that share an eQTL with a neighbouring gene. The module recomputes the correct gene label from `$inst` (SNP + exposure-dataset-ID pair) at load time (`load_mr_instrument_table()`, see §7) rather than trusting the cached column, and surfaces exactly what was corrected to the user via a UI note (`data_quality_note_{sex}`, §6).

---

## 2. Module Structure

### Main tab

`mod_mr_ui(id)` (lines 53–63) builds one `tabsetPanel` (`id = ns("mr_top_tabs")`, `type = "tabs"`) with three sub-tabs:

```r
tabsetPanel(
  id = ns("mr_top_tabs"), type = "tabs",
  tabPanel("MR overall", br(), mod_mr_overall_ui(ns)),
  tabPanel("MR female",  br(), mod_mr_sex_tab_ui(ns, "female")),
  tabPanel("MR male",    br(), mod_mr_sex_tab_ui(ns, "male"))
)
```

### Every sub-tab implemented in the code

| Sub-tab (exact label) | UI builder | Role |
|---|---|---|
| **"MR overall"** | `mod_mr_overall_ui(ns)` | The interactive, live-analysis tool: run MR for one gene (bundled or uploaded data), plus an optional batch screen across all bundled genes. |
| **"MR female"** | `mod_mr_sex_tab_ui(ns, "female")` → `uiOutput(ns("female_tab_body"))` | Displays the bundled dataset's own pre-computed female-stratum MR reference results (read-only, from saved CSVs). |
| **"MR male"** | `mod_mr_sex_tab_ui(ns, "male")` → `uiOutput(ns("male_tab_body"))` | Same as above, for the male stratum. |

There is no fourth sub-tab, no nested tabset inside "MR overall", and no additional top-level module tab beyond these three. The female/male tabs are built from one shared factory function (`render_sex_tab_body(sx)`, lines 473–494) applied with `sx = "female"` and `sx = "male"`, "so the female/male panels can never drift apart" (comment, line 356–361).

### Server-side associations per sub-tab

**MR overall** (server logic in `mod_mr_server`, lines 217–1055):
- Data-source switch: `input$data_source` ("project" vs "upload")
- Gene picker (bundled mode): `output$gene_ui`
- Upload parsing: `exp_df_r`, `out_df_r` (reactives), `output$exp_map_ui`, `output$out_map_ui`, `output$upload_data_summary`
- Filters: `input$pval_cut`, `input$fstat_cut`, `input$max_snps`, `input$ci_level`, `input$alpha_cut`, `input$mhc_mode`, `input$include_mode`, `input$run_presso`
- Reset: `observeEvent(input$reset_btn, ...)`
- Run: `mr_result <- eventReactive(input$run_btn, ...)`, dispatching to `mr_result_project()` or `mr_result_uploaded()`
- Outputs: `output$result_box_ui`, `output$scatter_diag_ui`, `output$instrument_mr_tables_ui`, `output$summary_ui`, `output$instrument_table`, `output$mr_table`, `output$scatter_plot`, `output$diagnostics_ui`, `output$funnel_plot`, `output$loo_plot`, `output$download_mr`
- Batch screen: `batch_result <- eventReactive(input$run_batch_btn, ...)`, `output$batch_screen_section`, `output$batch_summary_ui`, `output$batch_table`, `output$download_batch`

**MR female / MR male** (shared factories, parameterised on `sx`):
- Reveal gate: `sex_shown$female` / `sex_shown$male` (`reactiveValues`), set by `observeEvent(input$show_female_btn, ...)` / `observeEvent(input$show_male_btn, ...)`
- Data read: `project_survivors[[sx]]` (built once at server-init from `read_table_safe()`)
- Stat boxes: `render_proj_stats(sx)` → `output$proj_stats_{sx}`
- Forest plot: `build_proj_forest(sx)` / `render_proj_forest(sx)` → `output$proj_forest_{sx}`; download via `make_dl_forest_png(sx)` → `output$dl_proj_forest_{sx}`
- Results table: `render_proj_table(sx)` → `output$proj_table_{sx}`
- Downloads: `make_dl_survivors(sx)` → `output$dl_proj_survivors_{sx}`; `make_dl_xlsx(sx)` → `output$dl_proj_xlsx_{sx}`
- Data-quality note: `render_data_quality_note()` → `output$data_quality_note_{sx}`
- Tab body assembly: `render_sex_tab_body(sx)` → `output${sx}_tab_body`

---

## 3. User Interface / Features (by sub-tab)

### 3.1 "MR overall"

Built by `mod_mr_overall_ui(ns)` (lines 74–215), laid out as a `fluidRow()` with a 4-column left panel (inputs/filters) and an 8-column right panel (results), followed by two full-width sections below the row.

#### Left column — "Exposure & outcome data" box (lines 83–126)

- **ArthOChat shortcut** (`arthochat_shortcut_ui(...)`, compact) — "New to Mendelian randomisation? Ask ArthOChat."
- **`radioButtons(ns("data_source"), ...)`** — two choices:
  - "Bundled RA dataset (default)" → value `"project"` (default selected)
  - "Upload your own GWAS summary statistics" → value `"upload"`
- **Conditional panel, `data_source == "project"`:**
  - `uiOutput(ns("gene_ui"))` — gene selector (rendered server-side, see §4)
  - An info note naming the GWAS sources: outcome = Okada et al. 2014 RA GWAS (OpenGWAS `ieu-a-832`, linked), exposure = eQTLGen whole-blood cis-eQTL (Vösa et al. 2021, linked)
- **Conditional panel, `data_source == "upload"`:**
  - Descriptive text: two delimited files (CSV/TSV), one row per SNP; suggests asking ArthOChat to find an outcome dataset in OpenGWAS
  - `textInput(ns("upload_label"), "Exposure name (for labelling only)", value = "Uploaded exposure")`
  - `fileInput(ns("exp_file"), "Exposure file", accept = c(".csv", ".tsv", ".txt"))`
  - `uiOutput(ns("exp_map_ui"))` — column-mapping controls, rendered once a file is uploaded
  - `fileInput(ns("out_file"), "Outcome file", accept = c(".csv", ".tsv", ".txt"))`
  - `uiOutput(ns("out_map_ui"))` — column-mapping controls for the outcome file
  - `uiOutput(ns("upload_data_summary"))` — live "Your data:" summary note
  - `checkboxInput(ns("do_clump"), "LD-clump exposure SNPs before harmonising (needs network access)", value = FALSE)` + info tooltip
  - Conditional panel, `do_clump == TRUE`:
    - `selectInput(ns("clump_r2"), "Clump r² threshold", c("0.1","0.01","0.001 (default)"), selected = "0.001")`
    - `selectInput(ns("clump_kb"), "Clump window (kb)", c("250","1,000","10,000 (default)"), selected = "10000")`
    - `selectInput(ns("clump_pop"), "LD reference population", c(EUR, EAS, SAS, AFR, AMR), selected = "EUR")`

#### Left column — "Instrument & analysis filters" box (lines 127–194)

- Info note: filters are pre-selected to recommended defaults; running with no changes should match the bundled reference estimate for the chosen gene.
- `selectInput(ns("pval_cut"), "eQTL p-value threshold", ...)` — choices `5e-8` (default), `5e-9`, `1e-10`, `1e-12`, with tooltip noting every cached instrument already clears the default (min F = 29.7)
- `selectInput(ns("fstat_cut"), "Minimum F-statistic", ...)` — choices `10` (default), `0`, `30`, `100`, with tooltip
- `selectInput(ns("max_snps"), "Cap instruments used per gene (strongest by p-value)", ...)` — choices `Inf` (default, "No limit"), `20`, `10`, `5`, `3`, with tooltip
- `selectInput(ns("ci_level"), "Confidence interval level", ...)` — choices `0.90`, `0.95` (default), `0.99`, with tooltip clarifying it only widens/narrows CIs, not the point estimate or p-value
- `numericInput(ns("alpha_cut"), "Significance threshold (α) for flagging results", value = 0.05, min = 0.0001, max = 0.5, step = 0.005)`, with tooltip
- Conditional panel, `data_source == "project"`:
  - `radioButtons(ns("mhc_mode"), "Extended MHC region (chr6:25-34Mb)", ...)` — "Include (default)" vs "Exclude (sensitivity)", with tooltip on HLA-DRB1 confounding
- `checkboxInput(ns("include_mode"), "Also compute weighted-mode (exploratory, ≥3 instruments)", value = FALSE)`
- `checkboxInput(ns("run_presso"), "Also run MR-PRESSO outlier test (≥4 instruments)", value = FALSE)`, with tooltip
- Buttons: `actionButton(ns("run_btn"), "Run MR", icon("play"), class = "btn-primary btn-sm")`, `actionButton(ns("reset_btn"), "Reset to defaults", icon("rotate-left"), class = "btn-default btn-sm")`

#### Right column (lines 196–208)

- `uiOutput(ns("scatter_diag_ui"))` — pre-click: an empty-state card ("Nothing run yet"); post-click: "SNP effect scatter" box + "Sensitivity diagnostics" box (funnel + leave-one-out plots)
- `uiOutput(ns("result_box_ui"))` — pre-click: `NULL`; post-click: "Result" box with the textual summary

#### Below the row (lines 210–213)

- `uiOutput(ns("instrument_mr_tables_ui"))` — pre-click: `NULL`; post-click: "Instruments used in this run" table box + "MR estimates" table box with a CSV download button
- `uiOutput(ns("batch_screen_section"))` — bundled-data-only; shows a "Run batch screen" button pre-click, and a full results box (summary text, download button, table) post-click

### 3.2 "MR female" / "MR male"

Both built by `render_sex_tab_body(sx)` (lines 473–494), identical structure:

- **Pre-click state:** descriptive paragraph ("Pre-computed MR results for {sex} differential-expression candidate genes against the bundled RA GWAS (Okada et al. 2014, OpenGWAS ieu-a-832)") + `actionButton(ns("show_{sex}_btn"), "Show {sex} MR results", icon("chart-bar"), class = "btn-primary btn-sm")`
- **Post-click state:**
  - `uiOutput(ns("data_quality_note_{sex}"))` — relabelling disclosure (only rendered if `relabel_check$n_snp > 0`)
  - `uiOutput(ns("proj_stats_{sex}"))` — four `valueBox()` tiles: "Candidate genes (DE screen)", "MR-tested (cis instrument found)", "Survive FDR < 0.05", "ROBUST after MHC exclusion"
  - `downloadButton(ns("dl_proj_forest_{sex}"), "Forest plot (PNG)")`
  - `plotOutput(ns("proj_forest_{sex}"), height = 560)` wrapped in `withSpinner(...)`
  - `downloadButton(ns("dl_proj_survivors_{sex}"), "Survivors (CSV)")`
  - `downloadButton(ns("dl_proj_xlsx_{sex}"), "Full TABLE1-4 (XLSX)")`
  - `DT::dataTableOutput(ns("proj_table_{sex}"))`

### Messages / warnings / errors present in the code

- `validate(need(...))` calls that surface as in-place messages when instrument sets are empty (e.g. `mr_result_for_gene`, `mr_result_uploaded`, `build_proj_forest`, `batch_result`)
- Upload-mode informational notes in `output$summary_ui`: LD-clumping outcome (`applied`, `applied_no_change`, `api_error`, `no_match`), MHC-instrument counts, relabelling disclosure, low-instrument-count notices (`< 3 SNP`, PRESSO needing ≥4)
- `data_quality_note_{sex}` — relabelling disclosure with a hoverable info tooltip
- `mr_info_tip()` — a reusable small hover tooltip helper (native browser `title` attribute via `icon("circle-info")`) attached to several filter controls

---

## 4. Inputs

### 4.1 Bundled-dataset mode (`data_source == "project"`)

| Input | Widget | Source of choices | Default |
|---|---|---|---|
| `gene` | `selectInput(..., selectize = TRUE)` | `available_genes <- sort(unique(mr_dat_all$gene))` (~1,701 genes) | `"LPCAT2"` if present in the list, else the first gene alphabetically (`gene_ui`, lines 338–349) |

No file upload is required in this mode; the exposure/outcome data comes entirely from the bundled, pre-cached `MR_primary_objects.rds`.

### 4.2 Upload mode (`data_source == "upload"`)

| Input | Widget | Format | Required? |
|---|---|---|---|
| `upload_label` | `textInput` | free text, used only as a display label | Optional (defaults to `"Uploaded exposure"` if blank) |
| `exp_file` | `fileInput` | `.csv`, `.tsv`, `.txt` (delimited; parsed by `data.table::fread` via `read_uploaded_table()`) | **Required** |
| `out_file` | `fileInput` | `.csv`, `.tsv`, `.txt` | **Required** |
| Column-mapping inputs (`exp_snp`, `exp_beta`, `exp_se`, `exp_pval`, `exp_ea`, `exp_oa`, `exp_eaf`, and the `out_*` equivalents) | `selectInput`, generated by `gwas_col_map_ui()` | Column names from the uploaded file | **Required except `*_eaf`** (effect allele frequency is optional) |
| `do_clump` | `checkboxInput` | — | Default `FALSE` |
| `clump_r2` | `selectInput` | `"0.1"`, `"0.01"`, `"0.001"` | Default `"0.001"` |
| `clump_kb` | `selectInput` | `"250"`, `"1000"`, `"10000"` | Default `"10000"` |
| `clump_pop` | `selectInput` | `"EUR"`, `"EAS"`, `"SAS"`, `"AFR"`, `"AMR"` | Default `"EUR"` |

Each uploaded file's `fileInput` requires at least 4 columns to render the mapping UI at all (`gwas_col_map_ui()`, `global.R:1235`: `validate(need(!is.null(df) && ncol(df) >= 4, ...))`).

**Column auto-guessing** (`guess_gwas_col()`, `GWAS_COL_PATTERNS`, `global.R:1170–1184`): regex patterns (case-insensitive) are matched against uploaded column names to pre-select the most likely mapping, per field:

```r
GWAS_COL_PATTERNS <- list(
  snp  = c("^snp$", "^rsid$", "^rs_?id$", "variant"),
  ea   = c("^effect_allele$", "^ea$", "^a1$", "^alt$"),
  oa   = c("^other_allele$", "^oa$", "^a2$", "^ref$", "^nea$"),
  beta = c("^beta$", "^b$", "^effect$"),
  se   = c("^se$", "^standard_error$", "^stderr$"),
  pval = c("^p$", "^pval$", "^p_value$", "^pvalue$"),
  eaf  = c("^eaf$", "^freq$", "^maf$", "^effect_allele_freq"),
  n    = c("^n$", "^samplesize$", "^sample_size$", "^total_n$")
)
```

If no pattern matches, `guess_gwas_col()` returns `NA`, and `gwas_col_map_ui()` falls back to the first column (`cols[1]`) for required fields, or `""` ("(none)") for the optional EAF field.

Note: the `n` (sample size) pattern exists in `GWAS_COL_PATTERNS` but is **not** requested by `mod_mr.R`'s calls to `gwas_col_map_ui()` (no `extra_fields` argument is passed) — it is used by a different module (`mod_coloc.R`, per the comment at `global.R:1226–1228`), not by MR.

### 4.3 Shared filter inputs (both modes)

| Input | Type | Choices / range | Default | Used by |
|---|---|---|---|---|
| `pval_cut` | `selectInput` | `5e-8`, `5e-9`, `1e-10`, `1e-12` | `5e-8` | Both `mr_result_for_gene()` and `mr_result_uploaded()`, and `batch_result()` |
| `fstat_cut` | `selectInput` | `10`, `0`, `30`, `100` | `10` | Same as above |
| `max_snps` | `selectInput` | `Inf`, `20`, `10`, `5`, `3` | `Inf` | `apply_max_snps_cap()` |
| `ci_level` | `selectInput` | `0.90`, `0.95`, `0.99` | `0.95` | `current_ci_level()` → `estimate_mr_set()` |
| `alpha_cut` | `numericInput` | `0.0001`–`0.5`, step `0.005` | `0.05` | Significance flag in `summary_ui`, batch-screen FDR-significance count |
| `mhc_mode` | `radioButtons` | `"include"`, `"exclude"` | `"include"` | `mr_result_for_gene()`, `batch_result()` (project mode only) |
| `include_mode` | `checkboxInput` | — | `FALSE` | `estimate_mr_set(..., include_mode = ...)` |
| `run_presso` | `checkboxInput` | — | `FALSE` | `estimate_mr_set(..., run_presso = ...)` |

`observeEvent(input$reset_btn, ...)` (lines 498–507) resets exactly these seven filters (`pval_cut`, `fstat_cut`, `max_snps`, `ci_level`, `alpha_cut`, `mhc_mode`, `include_mode`, `run_presso`) to the defaults listed above.

### 4.4 Validation

- `req(gene)` / `req(input$exp_file, input$out_file)` and `req()` on every required column-mapping input (`mr_result_uploaded()`, lines 533–535) — the run silently waits until these are non-empty.
- `validate(need(...))` checks (surfaced as in-place error text, not a hard crash):
  - No cached instrument rows for the selected gene (`mr_result_for_gene`, line 634)
  - No instruments remain after p-value/F-statistic/MHC filtering (`mr_result_for_gene`, line 641; `mr_result_uploaded`, line 568)
  - Could not read one of the uploaded files (line 539)
  - Harmonisation found no overlapping SNPs between exposure and outcome files (line 556)
  - No SNPs survived harmonisation, e.g. all dropped as unresolved palindromic variants (line 562)
  - No reference results available for a given sex (`build_proj_forest`, line 401)
  - Reference results file not found (`render_proj_table`, line 434)
  - Batch screen invoked while `data_source != "project"` (line 982)
  - No genes have any instrument left after batch-screen filters (line 991)

### 4.5 How inputs move from UI into analysis functions

```
input$data_source ── switches between two branches inside mr_result <- eventReactive(input$run_btn, ...)

Bundled path:
  input$gene, input$pval_cut, input$fstat_cut, input$mhc_mode, input$max_snps
    → mr_result_for_gene(gene)
      → subsets mr_dat_all (loaded once via load_mr_instrument_table())
      → apply_max_snps_cap(d)
      → estimate_mr_set(d, include_mode = input$include_mode, ci_level = current_ci_level(),
                         run_presso = input$run_presso)

Upload path:
  input$exp_file/out_file, input$exp_*/out_* column maps, input$upload_label,
  input$do_clump, input$clump_r2/kb/pop, input$pval_cut, input$fstat_cut, input$max_snps
    → mr_result_uploaded()
      → exp_df_r()/out_df_r() (cached reactive reads via read_uploaded_table())
      → TwoSampleMR::format_data() ×2 → TwoSampleMR::harmonise_data(action = 2)
      → filter by pval.exposure/Fstat
      → optional ieugwasr::ld_clump()
      → apply_max_snps_cap(d)
      → estimate_mr_set(d, include_mode = ..., ci_level = ..., run_presso = ...)
```

---

## 5. Processing / Functionality

### 5.1 Overall reactive flow (bundled-dataset gene run)

```
UI: input$gene, filters set  →  actionButton "Run MR" clicked (input$run_btn)
        ↓
eventReactive: mr_result() dispatches to mr_result_project() = mr_result_for_gene(input$gene)
        ↓
req(gene) — guard
        ↓
d0 <- mr_dat_all[mr_dat_all$gene == gene, ]           # subset the pre-loaded, cached, harmonised table
validate(need(nrow(d0) > 0, ...))                       # no instrument for this gene → stop with message
        ↓
d <- d0[pval.exposure <= pval_cut & Fstat >= fstat_cut, ]   # instrument filters
if (mhc_mode == "exclude") d <- d[!d$MHC, ]                # optional MHC exclusion
validate(need(nrow(d) >= 1, ...))                           # nothing left → stop with message
        ↓
d <- apply_max_snps_cap(d)                              # optional cap: strongest-by-p-value SNPs only
        ↓
est <- estimate_mr_set(d, include_mode, ci_level, run_presso)   # THE statistical core (global.R)
        ↓
list(gene, d, n_before, n_after, n_mhc_dropped/retained, mhc_mode, pval_cut, fstat_cut, est,
     project_reference, uploaded = FALSE)
        ↓
mr_has_run(TRUE) is set by a separate observeEvent(input$run_btn, ...)
        ↓
UI renders: result_box_ui / scatter_diag_ui / instrument_mr_tables_ui, each reading mr_result()
```

### 5.2 Overall reactive flow (upload mode)

```
UI: exp_file + out_file uploaded, columns mapped, filters set  →  "Run MR" clicked
        ↓
mr_result_uploaded():
  req(files, all column-mapping inputs)
  exp_raw <- exp_df_r(); out_raw <- out_df_r()             # cached fread() reads
  validate(need(both non-null, "Could not read one of the uploaded files."))
        ↓
  exp_fmt <- TwoSampleMR::format_data(exp_raw, type = "exposure",
                snp_col=, beta_col=, se_col=, pval_col=, effect_allele_col=,
                other_allele_col=, eaf_col = mapped-or-"eaf")
  out_fmt <- TwoSampleMR::format_data(out_raw, type = "outcome", ... same pattern)
  exp_fmt$exposure <- label ; out_fmt$outcome <- "Uploaded outcome"
        ↓
  dat_up <- TwoSampleMR::harmonise_data(exp_fmt, out_fmt, action = 2)
  validate(need(nrow > 0, "no overlapping SNPs..."))
        ↓
  dat_up$gene <- label
  dat_up$Fstat <- (beta.exposure / se.exposure)^2           # F-statistic computed directly, not from GWAS metadata
  dat_up$MHC <- FALSE                                         # no MHC concept for arbitrary uploaded traits
  dat_up <- dat_up[dat_up$mr_keep, ]                          # TwoSampleMR's own harmonisation-quality flag
  validate(need(nrow > 0, "No SNPs survived harmonisation..."))
        ↓
  d <- dat_up[pval.exposure <= pval_cut & Fstat >= fstat_cut, ]
  validate(need(nrow(d) >= 1, "No SNPs remain after filtering..."))
        ↓
  [optional] if input$do_clump && nrow(d) > 1:
     ieugwasr::ld_clump(data.frame(rsid=SNP, pval=pval.exposure, id=label),
                         clump_kb=, clump_r2=, clump_p = 1, pop=)
     → three explicit outcomes tracked: "api_error" (call errored), "no_match" (0 rows returned),
       "applied"/"applied_no_change" (succeeded); on any non-"applied" outcome the UNCLUMPED set is kept
        ↓
  d <- apply_max_snps_cap(d)
        ↓
  est <- estimate_mr_set(d, include_mode, ci_level, run_presso)
        ↓
  list(gene = label, d, n_before, n_after, n_mhc_dropped=0, n_mhc_retained=0, mhc_mode="n/a",
       pval_cut, fstat_cut, est, clump_status, n_before_clump, n_after_clump,
       project_reference = mr_primary_all[0,] (empty), uploaded = TRUE)
```

### 5.3 The estimator hierarchy — `estimate_mr_set()` (`global.R:1269–1372`)

This is the statistical core, shared verbatim with `mod_crossancestry.R`'s live-upload arm (per the comment at `global.R:1261–1268`).

**Signature:** `estimate_mr_set(d, include_mode = FALSE, full = TRUE, ci_level = 0.95, run_presso = FALSE)`

**Logic by instrument count (`n_snp <- nrow(d)`):**

- **`n_snp == 1`:** Wald ratio computed directly:
  ```r
  b  <- beta.outcome / beta.exposure
  se <- abs(se.outcome / beta.exposure)
  p  <- 2 * pnorm(-abs(b / se))
  ```
  `primary_method <- "Wald ratio"`.

- **`n_snp >= 2`:** builds an `MendelianRandomization::mr_input()` object from `beta.exposure`, `se.exposure`, `beta.outcome`, `se.outcome`, `SNP`.
  - **IVW** (`MendelianRandomization::mr_ivw(mrobj, model = "random", alpha = alpha)`) is always computed; `primary_method <- "IVW"`. The code comment (lines 1290–1298) documents that `model = "random"` (not the package default `"fixed"`) was chosen because it reproduces the bundled cached SE/p bit-for-bit (verified example: LPCAT2, SE 0.0248 vs 0.0220 under "fixed"; point estimate identical either way).
  - **If `n_snp >= 3` and `full = TRUE`:**
    - **Weighted median** (`MendelianRandomization::mr_median(mrobj, alpha = alpha)`) — a bootstrap SE (10,000 draws, package default seed); the code notes this cannot be reproduced bit-for-bit against the bundled reference (never the primary estimate, so this is a disclosed limitation).
    - **MR-Egger** (`MendelianRandomization::mr_egger(mrobj, alpha = alpha)`) — estimate/SE/intercept match the bundled cached numbers exactly, but the p-values are **recomputed manually** using a t-distribution on `n_snp - 2` degrees of freedom (`stats::pt(-abs(estimate/se), df = n_snp - 2)`), because neither of the package's own p-value options reproduces the reference df (documented in the code comment, lines 1313–1322).
    - **Heterogeneity**: `Q`, `Q_df = n_snp - 1`, `Q_pval` taken from `ivw@Heter.Stat`.
    - **Pleiotropy**: MR-Egger intercept, its SE, its recomputed p-value, and `I²` (from `egg@I.sq`).
    - **If `include_mode = TRUE`:** weighted mode via `MendelianRandomization::mr_mbe(mrobj, alpha = alpha)`, wrapped in `tryCatch` (silently omitted on error).
  - **If `n_snp >= 4` and `full = TRUE` and `run_presso = TRUE`:** `MRPRESSO::mr_presso()` — simulation-based test (`NbDistribution = 1000`, `OUTLIERtest = TRUE`, `DISTORTIONtest = TRUE`, `SignifThreshold = alpha`, `seed = 2024`). Returns `global_p`, `raw_estimate`/`raw_p`, `corrected_estimate`/`corrected_p` (when outliers are found), and `n_outliers` (count of per-instrument outlier p-values below `alpha`). Wrapped in `tryCatch`; on error, `presso$error` carries the message.

**Return value:** `list(res_table, n_snp, heterogeneity, pleiotropy, primary_method, presso, ci_level)`, where `res_table` is a data frame with one row per computed method (`method`, `estimate`, `se`, `ci_low`, `ci_high`, `p`, `primary` [logical, `TRUE` for the pre-specified primary method]).

**Primary-method hierarchy** (as stated in the UI text, line 718): *IVW > Wald ratio > weighted median > MR-Egger* — i.e., IVW is primary whenever ≥2 SNPs are available, Wald ratio when exactly 1 SNP is available; weighted median and MR-Egger are never chosen as primary, "fixed in advance, never chosen by which p-value looks best."

### 5.4 Instrument-cap logic — `apply_max_snps_cap()` (lines 520–524)

```r
apply_max_snps_cap <- function(d) {
  cap <- suppressWarnings(as.numeric(input$max_snps %||% "Inf"))
  if (is.na(cap) || !is.finite(cap) || nrow(d) <= cap) return(d)
  d[order(d$pval.exposure), , drop = FALSE][seq_len(cap), , drop = FALSE]
}
```
If a finite cap is set and the filtered instrument set exceeds it, only the `cap` strongest instruments (lowest `pval.exposure`) are retained.

### 5.5 Loading the bundled instrument table — `load_mr_instrument_table()` (`global.R:1389–1408`)

```r
load_mr_instrument_table <- function() {
  mr_obj <- readRDS(MR_PRIMARY_OBJECTS_RDS)

  dat <- as.data.frame(mr_obj$dat)
  dat$gene <- mr_obj$inst$gene[match(
    paste(dat$SNP, dat$id.exposure),
    paste(mr_obj$inst$SNP, mr_obj$inst$id.exposure)
  )]
  dat <- dat[dat$mr_keep, , drop = FALSE]

  # relabel_check: compares cached vs. correct gene label per (SNP, id.exposure) pair
  ...
  list(dat = dat, primary = as.data.frame(mr_obj$primary), relabel_check = relabel_check)
}
```
Called once at module-server startup (`mod_mr.R:230`): `mr_loaded <- load_mr_instrument_table()`. Produces `mr_dat_all` (per-SNP harmonised instrument rows, gene labels corrected, `mr_keep`-filtered), `mr_primary_all` (the bundled pipeline's own pre-computed primary estimate per gene, used for the "bundled reference pipeline" comparison text in `summary_ui`), and `relabel_check` (disclosure data).

### 5.6 Reference-results loading (MR female / MR male tabs)

At server startup (lines 269–287), independent of any button click:

```r
project_survivors <- lapply(list(female="female", male="male"), function(sx) {
  d <- read_table_safe(sprintf("MR_MHC_sensitivity_%s.csv", sx))
  req_cols <- c("gene","MHC_gene","nSNP_primary","method_primary","OR_primary",
                "OR_lo_primary","OR_hi_primary","p_primary","FDR_primary","verdict")
  if (is.null(d) || !all(req_cols %in% colnames(d))) return(NULL)
  d <- d[d$FDR_primary < 0.05, req_cols, drop = FALSE]
  d$direction <- ifelse(d$OR_primary > 1, "risk (OR>1)", "protective (OR<1)")
  d$verdict_short <- verdict_short(d$verdict)
  d[order(d$p_primary), , drop = FALSE]
})
```

`read_table_safe(filename, dir = TABLES_DIR)` (`global.R:1151–1155`) reads a CSV with `data.table::fread`, returning `NULL` if the file does not exist. No recomputation occurs — this is a direct read of the pipeline's saved output (per the code comment, lines 241–250, which also states this reproduces `FS_input_{sex}.csv` / `MR_causal_FDR_{sex}.csv` exactly: 32 female / 25 male rows, verified).

`project_n_tested` and `project_n_candidates` (lines 280–287) read the total row count of the same sensitivity file and the unique gene count of `candidates_{sex}_disease.csv` respectively, for the value-box tiles.

**Verdict classification** (`verdict_short()`, lines 251–259): maps the raw `verdict` column (string prefixes) to one of five short labels:
```r
dplyr::case_when(
  startsWith(v, "ROBUST")        ~ "ROBUST",
  startsWith(v, "UNTESTABLE")    ~ "UNTESTABLE w/o MHC",
  startsWith(v, "MHC-DEPENDENT") ~ "MHC-DEPENDENT",
  startsWith(v, "FDR-RANK")      ~ "FDR-rank only",
  TRUE                           ~ "ns in both"
)
```
Colour mapping (`verdict_color`) uses `ARTHOMIX_STATUS$good/warning/critical` for ROBUST/UNTESTABLE/MHC-DEPENDENT, and a neutral grey (`"#8A929C"`) for the other two.

### 5.7 Batch screen — `batch_result()` (lines 981–1016)

```
validate(need(data_source == "project", ...))          # bundled-data only
d_all <- mr_dat_all[pval.exposure <= pval_cut & Fstat >= fstat_cut, ]
if (mhc_mode == "exclude") d_all <- d_all[!d_all$MHC, ]
genes <- sort(unique(d_all$gene))
validate(need(length(genes) > 0, ...))
split_d <- split(d_all, d_all$gene)

for each gene:
  dg <- apply_max_snps_cap(split_d[[gene]])
  est <- tryCatch(estimate_mr_set(dg, include_mode = FALSE, full = FALSE, ci_level = current_ci_level()),
                   error = function(e) NULL)
  if not NULL: record gene, primary method, nSNP, estimate, se, p, OR = exp(estimate),
               OR_lo/OR_hi = exp(ci_low/ci_high), MHC_gene = any(dg$MHC), sensitivity_testable = n_snp >= 3

out$FDR <- stats::p.adjust(out$p, method = "BH")          # BH-FDR across genes tested THIS session, under current filters
out <- out[order(out$p), ]
```

Progress is surfaced via `withProgress()`/`setProgress()`, updated every 25 genes or on the last gene. The code comment (lines 973–980) explicitly distinguishes this session-local BH-FDR screen from the pre-computed, sex-stratified FDR reference in the "MR female"/"MR male" tabs — "kept explicitly distinct... in every label this produces." Note `full = FALSE` is passed to `estimate_mr_set()` here, so only the IVW/Wald-ratio primary estimate is computed per gene (no weighted median, MR-Egger, heterogeneity, pleiotropy, or PRESSO) — a batch-screen-specific performance choice visible directly in the call.

### 5.8 LD clumping (upload mode only)

```r
ieugwasr::ld_clump(
  data.frame(rsid = d$SNP, pval = d$pval.exposure, id = label),
  clump_kb = as.numeric(input$clump_kb), clump_r2 = as.numeric(input$clump_r2),
  clump_p = 1, pop = input$clump_pop
)
```
`clump_p = 1` is hardcoded (i.e., clumping is applied without an additional p-value gate beyond what has already been filtered). Three explicit outcomes are tracked and disclosed to the user (`clump_status`): `"api_error"` (the call itself errored — network/token issue), `"no_match"` (the call succeeded but returned zero rows — none of the uploaded SNP IDs matched the chosen population's LD reference panel), `"applied"` / `"applied_no_change"` (succeeded, with or without any SNPs actually removed). On `"api_error"` or `"no_match"`, the **unclumped** instrument set is kept and the run proceeds — never blocked.

---

## 6. Outputs (by sub-tab)

### 6.1 "MR overall" — per-run outputs

**Result box** (`output$result_box_ui` → `output$summary_ui`, lines 673–800), rendered only after `mr_has_run()` is `TRUE`:

- Upload-mode disclosure note (no MHC flag/bundled reference applies to custom data)
- LD-clumping status note (one of four messages depending on `clump_status`)
- Instrument count sentence: gene name, number of instruments after filtering, number before filtering ("cached" for bundled mode, "harmonised" for upload mode)
- MHC-exclusion note (count of MHC instruments excluded) or MHC-retention warning (count retained, "interpret this estimate as associated, not necessarily causal...")
- **Primary estimate line:** method name, `b` (point estimate, 4 dp), CI at the chosen level (4 dp), p-value (3 significant figures, scientific notation), and a coloured "(significant / not significant at α = ...)" tag (green if `p < alpha_cut`, grey otherwise) — colour uses `ARTHOMIX_STATUS$good`
- **Cochran's Q** line (only if `heterogeneity` is non-null): `Q` (2 dp), df, p-value
- **MR-Egger intercept** line (only if `pleiotropy` is non-null): intercept (4 dp), p-value, I²
- **MR-PRESSO** lines (only if `presso` is non-null and has no error): global test p-value, outlier count, and (if any outliers) the outlier-corrected estimate/p-value; if `presso$error` is set, a warning note with the error message instead
- Note if `n_snp < 3` (only IVW/Wald ratio shown; no Egger/median/Q/intercept)
- Note if `3 <= n_snp < 4` and PRESSO was requested (not run — needs ≥4)
- Relabelling caveat if the current gene is in `relabel_check$genes`
- **Bundled-reference comparison** (only if `project_reference` has rows, i.e. bundled mode and the gene has a cached primary estimate): "The bundled reference pipeline... estimated this gene at b = ..., p = ..., using N instrument(s) [(MHC-flagged)]"

**SNP effect scatter** (`output$scatter_plot`, via `scatter_plot_obj()`, lines 840–863):
- `ggplot`: x = `beta.exposure`, y = `beta.outcome`, coloured by MHC status (`"MHC region"` red / `"Non-MHC"` blue, via `ARTHOMIX_COLORS`)
- Horizontal/vertical error bars at ±1.96×SE for each axis
- Points, plus a dashed red line through the origin with slope = the primary estimate (`geom_abline(intercept = 0, slope = primary_row$estimate, ...)`)
- Title: gene name, SNP count, primary method name
- Axis labels: "SNP effect on expression (cis-eQTL beta)" / "SNP effect on RA risk (GWAS beta)"

**Sensitivity diagnostics** (`output$diagnostics_ui`), shown only when `n_snp >= 3`; otherwise a note explaining the instrument count is insufficient:
- **Funnel plot** (`output$funnel_plot`): per-SNP Wald ratio (`beta.outcome/beta.exposure`) on x, precision (`1/|se.outcome/beta.exposure|`) on y; vertical dashed red line at the primary estimate
- **Leave-one-out plot** (`output$loo_plot`, via `loo_result()`): re-fits `MendelianRandomization::mr_ivw(model="random")` with each SNP excluded in turn, plus the all-SNP estimate as a reference row labelled `"All SNPs (IVW)"`; horizontal error bars (95% CI), points coloured red for the all-SNP row, blue for each leave-one-out row

**Instrument table** (`output$instrument_table`, lines 807–826): columns present depend on what exists in `res$d` — `SNP`, `Chr` (`chr.exposure`), `Position (GRCh37)` (`pos.exposure`), `β exposure`, `SE exposure`, `eQTL/exposure p-value`, `F-statistic`, `In MHC`. Upload-mode data lacks `chr.exposure`/`pos.exposure` (not requested from `format_data()`) and always shows `MHC = No`, so those columns/rows differ from bundled-mode output accordingly — built by intersecting the candidate column list with `colnames(res$d)`. Rounding: p-value to 3 significant figures, β/SE to 4 dp, F-statistic to 1 dp; `In MHC` rendered as "Yes"/"No".

**MR estimates table** (`output$mr_table`, lines 828–838): `res_table` from `estimate_mr_set()`, renamed to `Method`, `Estimate (b)`, `SE`, `95% CI low`, `95% CI high`, `p-value`, `Project primary` (a "✓" mark for the pre-specified primary row). Rounding: estimate/SE/CI bounds to 4 dp, p-value to 3 significant figures.

**Downloads:**
- `output$download_mr` — CSV of `mr_result()$est$res_table`, filename `mr_<gene>.csv`
- No PNG download is offered for the scatter/funnel/LOO plots in "MR overall" — **not implemented in the current module code** (only the sex-stratified forest plot has a PNG download; see §6.2).

**Batch screen** (`output$batch_screen_section`, gated behind `batch_has_run()`):
- Summary sentence (`output$batch_summary_ui`): number of genes tested, number significant at the current FDR threshold
- Table (`output$batch_table`): columns renamed to `Gene`, `Method`, `nSNP`, `b`, `SE`, `p`, `OR`, `OR low`, `OR high`, `In MHC`, `≥3 SNP (pleiotropy testable)`, `FDR`; `DT` with `filter = "top"`
- Download: `output$download_batch` — CSV of the full batch table, filename `mr_batch_screen.csv`

### 6.2 "MR female" / "MR male" — outputs

Shown only after "Show {sex} MR results" is clicked:

- **Data-quality note** (only if the cached relabelling issue affected any SNPs): a hover tooltip with the count of relabelled SNPs/genes and, if applicable, which survivor genes were affected; plus an inline text summary
- **Value boxes** (four): candidate gene count (DE screen), MR-tested gene count, count surviving `FDR_primary < 0.05`, count of those with a `"ROBUST"` MHC-sensitivity verdict
- **Forest plot** (`build_proj_forest(sx)`): x = `OR_primary` (log scale), y = gene (ordered by OR), coloured by `verdict_short`; vertical dashed reference line at OR = 1; title states sex, gene count, and "FDR < 0.05 vs Okada 2014 RA GWAS"; subtitle: "Colour = verdict under MHC-instrument-exclusion sensitivity re-run"; uses `theme_arthomix()`
- **Results table** (`render_proj_table`): columns `Gene`, `Direction` ("risk (OR>1)"/"protective (OR<1)"), `OR` (3 dp), `95% CI` (formatted string), `p` (3 sig figs), `FDR` (3 sig figs), `nSNP`, `Method`, `In MHC` ("Yes"/"No"), `Verdict` (HTML-rendered coloured badge via `verdict_badge()`); `DT` with `filter = "top"`, `escape = FALSE`
- **Downloads:**
  - Forest plot PNG (`make_dl_forest_png`): `ggsave()` at 9×max(5, n_genes×0.16) inches, 300 dpi, white background, filename `MR_{sex}_forest_plot.png`
  - Survivors CSV (`make_dl_survivors`): `write.csv(project_survivors[[sx]], ...)`, filename `MR_{sex}_RA_reference_survivors.csv`
  - Full tables XLSX (`make_dl_xlsx`): `file.copy()` of a pre-existing file `MR_{sex}_all_tables.xlsx` from `TABLES_DIR`, filename `MR_{sex}_all_tables.xlsx` — this is a direct file copy, not a generated export

### 6.3 Cross-module output (`results$mr`)

Two `observeEvent`s write into the shared `results` object (passed into `mod_mr_server(id, dataset, results)`), used elsewhere in the app (consumer not documented here per the task's module-scope restriction, since it is outside `mod_mr.R`):

- On every completed single-gene run (`observeEvent(mr_result(), ...)`, lines 725–737): appends an entry to `results$mr$genes_tested[[gene]]` — `gene`, `method`, `n_snp`, `estimate` (4 dp), `p` (3 sig figs), `filters` (a formatted string of the active p-value/F-stat/MHC filters)
- On every completed batch run (`observeEvent(batch_result(), ...)`, lines 1018–1029): sets `results$mr$batch_screen` — `n_genes_tested`, `n_fdr_significant` (count with `FDR < alpha_cut`), `filters` (formatted string including `alpha`)

---

## 7. Code-Level Documentation

### `mr_info_tip(text)` — `mod_mr.R:44`
- **Purpose:** small hover tooltip for a filter's rationale
- **Inputs:** `text` (string)
- **Processing:** wraps a Font Awesome `circle-info` icon in a `tags$span` with a native `title` attribute (no extra JS/CSS)
- **Output:** a `shiny.tag` for inline placement next to a filter control
- **Where called:** next to `pval_cut`, `fstat_cut`, `max_snps`, `ci_level`, `alpha_cut`, `mhc_mode`, `do_clump`, `run_presso` (all in `mod_mr_overall_ui`)
- **UI element controlled:** none directly — purely decorative/informational

### `mod_mr_config` — `mod_mr.R:46–51`
- **Purpose:** module registration metadata (id, group, title, description, icon) consumed by the app's module registry (outside this file's scope)
- **Not itself a function; a static list**

### `mod_mr_ui(id)` — `mod_mr.R:53–63`
- **Purpose:** top-level UI entry point for the whole MR module
- **Inputs:** `id` (Shiny module id)
- **Processing:** creates the namespace (`ns <- NS(id)`), builds the three-tab `tabsetPanel`
- **Output:** a `tagList`
- **Called by:** the app's module-loading mechanism (outside this file)

### `mod_mr_sex_tab_ui(ns, sx)` — `mod_mr.R:70–72`
- **Purpose:** thin UI wrapper for a sex-stratified tab
- **Inputs:** `ns` (namespace function), `sx` (`"female"` or `"male"`)
- **Processing:** returns a single `uiOutput` bound to `{sx}_tab_body`
- **Output:** `shiny.tag`
- **Called by:** `mod_mr_ui()`

### `mod_mr_overall_ui(ns)` — `mod_mr.R:74–215`
- **Purpose:** builds the entire "MR overall" tab UI (documented in full in §3.1)
- **Inputs:** `ns`
- **Output:** `tagList`

### `mod_mr_server(id, dataset, results)` — `mod_mr.R:217–1055`
- **Purpose:** the module's server function; all reactive/statistical logic lives here
- **Inputs:** `id`, `dataset` (shared reactiveValues, **not read by this module** — no reference to `dataset$expr`/`dataset$meta` appears anywhere in the file, since MR uses its own bundled/uploaded GWAS data, not the app's expression matrix), `results` (shared reactiveValues, written to as described in §6.3)
- **Processing:** everything documented in §2–§6
- **Output:** none (side effects via `output$*` and `results$mr`)

### `estimate_mr_set(d, include_mode, full, ci_level, run_presso)` — `global.R:1269–1372`
Fully documented in §5.3. Called from `mr_result_for_gene()`, `mr_result_uploaded()`, and `batch_result()` (with `full = FALSE`).

### `load_mr_instrument_table()` — `global.R:1389–1408`
Fully documented in §5.5. Called once at server startup (`mod_mr.R:230`).

### `read_table_safe(filename, dir = TABLES_DIR)` — `global.R:1151–1155`
- **Purpose:** safe CSV reader; returns `NULL` instead of erroring if the file is absent
- **Inputs:** `filename`, `dir` (defaults to `TABLES_DIR`)
- **Processing:** `file.exists()` check, then `data.table::fread()`
- **Output:** `data.frame` or `NULL`
- **Called by:** `project_survivors` construction, `project_n_tested`, `project_n_candidates` (all at server startup)

### `guess_gwas_col(cols, patterns)` / `GWAS_COL_PATTERNS` — `global.R:1170–1184`
Fully documented in §4.2.

### `read_uploaded_table(path)` — `global.R:1186`
- **Purpose:** parse an uploaded delimited file
- **Processing:** `as.data.frame(data.table::fread(path, showProgress = FALSE))`, wrapped in `tryCatch` returning `NULL` on error
- **Called by:** `exp_df_r`, `out_df_r` reactives (`mod_mr.R:312–313`)

### `gwas_col_map_ui(ns, file_input, df_reactive, prefix, label, extra_fields)` — `global.R:1229–1259`
Fully documented in §3.1/§4.2. Called twice in this module: `output$exp_map_ui <- gwas_col_map_ui(ns, reactive(input$exp_file), exp_df_r, "exp", "Exposure file")` and the `out_*` equivalent — **no `extra_fields` passed**, so only the seven core fields (snp/beta/se/pval/ea/oa/eaf) are requested, not `n`.

### `verdict_short(v)` / `verdict_color` / `verdict_badge(v)` — `mod_mr.R:251–267`
Fully documented in §5.6.

### `project_survivors`, `project_n_tested`, `project_n_candidates` — `mod_mr.R:269–287`
Server-startup constants (not reactive; computed once, since they read static CSVs via `read_table_safe`).

### `render_data_quality_note()`, `render_proj_stats(sx)`, `build_proj_forest(sx)`, `render_proj_forest(sx)`, `make_dl_forest_png(sx)`, `render_proj_table(sx)`, `make_dl_survivors(sx)`, `make_dl_xlsx(sx)`, `render_sex_tab_body(sx)` — `mod_mr.R:362–496`
Factory functions applied to `sx = "female"` and `sx = "male"` to build each output pair without duplicating logic. Fully documented in §5.6/§6.2.

### `current_ci_level()` / `apply_max_snps_cap(d)` — `mod_mr.R:519–524`
Fully documented in §5.4.

### `mr_result_uploaded()` — `mod_mr.R:532–621`
Fully documented in §5.2.

### `mr_result_for_gene(gene)` / `mr_result_project()` — `mod_mr.R:627–659`
Fully documented in §5.1. `mr_result_project()` is a zero-argument wrapper calling `mr_result_for_gene(input$gene)`.

### `mr_result` (`eventReactive`) — `mod_mr.R:661–663`
```r
mr_result <- eventReactive(input$run_btn, {
  if (identical(input$data_source, "upload")) mr_result_uploaded() else mr_result_project()
})
```
The single dispatch point for both analysis branches; fires only on `input$run_btn`.

### `mr_has_run` (`reactiveVal`) — `mod_mr.R:670–671`
Gates whether any result UI renders at all, independent of `mr_result()`'s own value — set `TRUE` by the same `observeEvent(input$run_btn, ...)`.

### `scatter_plot_obj()`, `loo_result()`, `funnel_data()` — `mod_mr.R:840–901`
Fully documented in §6.1 (visualization/diagnostic reactives).

### `batch_result` (`eventReactive`) — `mod_mr.R:981–1016`
Fully documented in §5.7.

---

## 8. Line-by-Line Understanding — selected code blocks

### 8.1 UI code — data-source switch and conditional panels (`mod_mr.R:85–125`)

```r
radioButtons(
  ns("data_source"), NULL,
  choiceNames = list(
    tagList(icon("database"), " Bundled RA dataset (default)"),
    tagList(icon("upload"), " Upload your own GWAS summary statistics")
  ),
  choiceValues = list("project", "upload"), selected = "project"
),
conditionalPanel(
  condition = sprintf("input['%s'] == 'project'", ns("data_source")),
  uiOutput(ns("gene_ui")), ...
),
conditionalPanel(
  condition = sprintf("input['%s'] == 'upload'", ns("data_source")),
  ...
)
```
This is pure **UI code**. `radioButtons` defines the switch that determines which of the two mutually exclusive input panels is visible; `choiceNames`/`choiceValues` decouple the displayed label (with icon) from the underlying value used in server logic (`"project"`/`"upload"`). The two `conditionalPanel()` calls use a JS condition string referencing the namespaced input id (`input['mr-data_source']`), evaluated client-side — both panels exist in the DOM at all times, but only one is visible. This is why server code must still check `input$data_source` explicitly (client-side visibility does not imply the hidden inputs are unset or invalid).

### 8.2 Server/reactive code — cached, single-parse uploaded-file reads (`mod_mr.R:304–313`)

```r
exp_df_r <- reactive({ req(input$exp_file); read_uploaded_table(input$exp_file$datapath) })
out_df_r <- reactive({ req(input$out_file); read_uploaded_table(input$out_file$datapath) })
```
This is **server/reactive code**. Wrapping the file parse in `reactive()` means Shiny's reactive graph caches the result and only re-runs `read_uploaded_table()` when `input$exp_file` (or `input$out_file`) itself changes — not on every downstream read. The code comment (lines 304–311) states this was a deliberate fix: three consumers (column-mapping UI, "Your data" summary, and the MR run itself) previously each called the read function independently, tripling the parse cost on a large file.

### 8.3 Statistical/analysis code — IVW model choice (`global.R:1290–1299`)

```r
ivw <- MendelianRandomization::mr_ivw(mrobj, model = "random", alpha = alpha)
methods[["IVW"]] <- c(estimate = ivw@Estimate, se = ivw@StdError, ci_low = ivw@CILower,
                       ci_high = ivw@CIUpper, p = ivw@Pvalue)
primary_method <- "IVW"
```
This is **statistical/analysis code**. `model = "random"` selects a random-effects (multiplicative-dispersion) IVW estimator, deviating from the package's own default (`"fixed"`). The surrounding comment documents this was chosen empirically — reproducing the bundled reference pipeline's cached SE/p values bit-for-bit — rather than being an assumption. `alpha` here is `1 - ci_level` (computed at the top of `estimate_mr_set()`), so `mr_ivw()`'s own CI construction is driven directly by the UI's `ci_level` filter.

### 8.4 Statistical/analysis code — manual MR-Egger p-value recomputation (`global.R:1312–1323`)

```r
egg <- MendelianRandomization::mr_egger(mrobj, alpha = alpha)
egger_df <- n_snp - 2
egger_slope_p <- 2 * stats::pt(-abs(egg@Estimate / egg@StdError.Est), df = egger_df)
egger_int_p   <- 2 * stats::pt(-abs(egg@Intercept / egg@StdError.Int), df = egger_df)
methods[["MR-Egger"]] <- c(estimate = egg@Estimate, se = egg@StdError.Est,
                            ci_low = egg@CILower.Est, ci_high = egg@CIUpper.Est, p = egger_slope_p)
```
This overrides the package's own p-values with a manually computed two-sided t-test p-value on `n_snp - 2` residual degrees of freedom — matching what `TwoSampleMR::mr_egger_regression()` (via `summary(lm(...))`) would produce, per the code comment. This is a deliberate, documented correction, not an unexplained deviation.

### 8.5 Visualization code — SNP effect scatter with primary-estimate slope line (`mod_mr.R:845–856`)

```r
p <- ggplot(d, aes(x = beta.exposure, y = beta.outcome, color = mhc_label)) +
  geom_hline(yintercept = 0, color = "#d5d9de") +
  geom_vline(xintercept = 0, color = "#d5d9de") +
  geom_errorbar(aes(ymin = beta.outcome - 1.96 * se.outcome, ymax = beta.outcome + 1.96 * se.outcome), width = 0, alpha = 0.6) +
  geom_errorbarh(aes(xmin = beta.exposure - 1.96 * se.exposure, xmax = beta.exposure + 1.96 * se.exposure), height = 0, alpha = 0.6) +
  geom_point(size = 2.6) +
  geom_abline(intercept = 0, slope = primary_row$estimate, color = "#c0392b", linetype = "dashed") +
  ...
```
This is **visualization code**. `1.96 * se` is a fixed 95% error-bar width for the per-SNP effect estimates in the scatter — note this is independent of the user's `ci_level` filter (which only affects the estimate/CI reported in the "Result" panel and the leave-one-out plot's IVW refits, not this scatter's error bars). The `geom_abline` slope is drawn through the origin at the primary method's point estimate, giving a visual read of the causal-effect line the primary estimator fit through the instrument cloud.

### 8.6 Validation/error-handling code — three-way LD-clumping outcome (`mod_mr.R:588–609`)

```r
clump_status <- "off"
n_before_clump <- nrow(d)
if (isTRUE(input$do_clump) && nrow(d) > 1) {
  clumped <- tryCatch(
    ieugwasr::ld_clump(data.frame(rsid = d$SNP, pval = d$pval.exposure, id = label),
                        clump_kb = ..., clump_r2 = ..., clump_p = 1, pop = ...),
    error = function(e) NULL
  )
  if (is.null(clumped)) {
    clump_status <- "api_error"
  } else if (nrow(clumped) == 0) {
    clump_status <- "no_match"
  } else {
    d <- d[d$SNP %in% clumped$rsid, , drop = FALSE]
    clump_status <- if (nrow(d) < n_before_clump) "applied" else "applied_no_change"
  }
}
```
This is **validation/error-handling code**. It distinguishes a hard API/network failure (`tryCatch` caught an error → `NULL`) from a successful-but-empty response (`nrow(clumped) == 0`, e.g. no rsIDs matched the chosen population's LD panel) — two different failure modes that the code comment (lines 573–587) says were confirmed against the app's own test fixtures (synthetic rsIDs clump to zero, real cached rsIDs clump successfully). In every non-"applied" branch, the instrument set `d` is left unmodified (unclumped), so the run always proceeds rather than being blocked by a clumping failure.

### 8.7 Download/export code — forest-plot PNG sized to gene count (`mod_mr.R:422–428`)

```r
make_dl_forest_png <- function(sx) downloadHandler(
  filename = function() sprintf("MR_%s_forest_plot.png", sx),
  content = function(file) {
    n_genes <- max(1, nrow(project_survivors[[sx]]))
    ggsave(file, plot = build_proj_forest(sx), width = 9, height = max(5, n_genes * 0.16), dpi = 300, bg = "white", limitsize = FALSE)
  }
)
```
This is **download/export code**. The plot height scales with the number of surviving genes (`max(5, n_genes * 0.16)` inches), so a sex stratum with more FDR-significant genes gets a taller export rather than a fixed-size plot with overlapping gene labels; `limitsize = FALSE` explicitly permits `ggsave()` to exceed its normal maximum-dimension safety check for large gene counts. `build_proj_forest(sx)` is the exact same function used for the on-screen `renderPlot`, so the downloaded PNG and the on-screen plot can never diverge (stated design intent in the comment above this block, lines 393–398).

---

## 9. Explicitly Not Implemented in the Current Module Code

The following are noted here because they are plausible MR-related features that a reader might expect but which **do not appear anywhere in `mod_mr.R`**:

- No PNG/image download for the "MR overall" tab's SNP effect scatter, funnel plot, or leave-one-out plot (only the sex-stratified forest plot has a `downloadButton` for its plot).
- No re-derivation of cis-restriction from raw GWAS/eQTL summary statistics for the bundled dataset — the module explicitly subsets an already cis-filtered cached object rather than recomputing cis boundaries (stated directly in the header comment, `mod_mr.R:14–18`).
- No sample-size (`N`) column mapping is requested in this module's calls to `gwas_col_map_ui()` (the `extra_fields` argument is omitted), even though `GWAS_COL_PATTERNS` defines an `n` pattern (used by a different module).
- No steiger-filtering, no MR-RAPS, no contamination-mixture method, no MR-cML, no other `TwoSampleMR`/`MendelianRandomization` estimator beyond Wald ratio, IVW, weighted median, MR-Egger, and (optionally) weighted mode — these are the only methods `estimate_mr_set()` computes.
- No trans-eQTL instruments — only cis-restricted instruments are used (per the header comment).
- No multivariable MR (MVMR) — this module performs univariable two-sample MR only.
- No user control over the MR-PRESSO simulation count (`NbDistribution = 1000`) or seed (`2024`) — both are hardcoded.
- No re-fitting of the batch screen's secondary methods (median/Egger/heterogeneity/pleiotropy/PRESSO) — `batch_result()` calls `estimate_mr_set(..., full = FALSE)`, so only the primary point estimate is computed per gene across the whole screen.
