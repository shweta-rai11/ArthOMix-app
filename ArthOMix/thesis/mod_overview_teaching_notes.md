# `mod_overview.R` — Full Teaching, Audit, and Thesis-Documentation Notes

File: `ArthOMix/R/transcriptomics/mod_overview.R` (713 lines)
Prepared: 2026-08-25
Companion file: [`mod_dataset_teaching_notes.md`](mod_dataset_teaching_notes.md) — read that one first if you haven't; this file assumes you already know the `dataset` reactiveValues contract (`expr`/`meta`/`source`, `staged_expr`/`staged_meta`/`staged_source`) it established.

---

## 1. PURPOSE

### Why this file exists

`mod_overview.R` is the **"Overview and Datasets"** sub-module — the app's read-mostly window onto (a) the catalog of raw GEO series the whole example pipeline is built from, and (b) quality-control diagnostics for **any single dataset the user picks**, independent of whatever is currently the app-wide active dataset. It answers three different questions a scientist (or a thesis committee) would ask before trusting any downstream result: *where did this data come from*, *what does the raw sample metadata/expression matrix actually look like*, and *is this data clean enough to analyze*.

### The problem it solves

A pipeline built around one merged, batch-corrected training cohort still needs a place to (1) show provenance for the raw GEO series that cohort was built from, (2) let a user inspect *any one* of those sources — or their own uploaded data — completely unfiltered, before any merge/normalization decision has been made, and (3) run the same missing-value / outlier / normalization-need checks the pipeline's own preprocessing script (`scripts/00_shared/eda.R`, per the file's header comment) already runs internally, but interactively and on-demand, on whichever dataset the user is curious about. This file is that inspection layer — deliberately **read-only** with respect to the app's shared active dataset (only one narrow exception, discussed in §2.7).

### Role in the app / module hierarchy

- App: ArthOMix Explorer (`server.R`, `global.R`, `ui.R`)
  - Top-level omics area: **Transcriptomics** (`R/transcriptomics/`)
    - Sub-module: **Overview and Datasets** (`mod_overview.R`) ← this file — registered in `TX_MODULES` (`R/submodules_registry.R:9`), **first** in that list
    - Sibling, but architecturally distinct: **Dataset** (`mod_dataset.R`) — the tab where a user actually *loads/stages* data (see the companion notes file). **Important, code-verified naming trap:** the sidebar's "Datasets" nav item (`ui.R:1335`) and this module's own first sub-tab (also literally titled "Datasets", `mod_overview.R:33`) are **not** the same UI element as the always-visible "Dataset" tab built by `mod_dataset_ui()` (`ui.R:1360`). All three names are almost identical and refer to two different files with two different purposes — see §7, Issue O-1.

### Inputs it receives

- `id` — `"overview"` (from `mod_overview_config$id`), namespaced to `"tx_overview"` when instantiated (`server.R:141`: `m$server(paste0("tx_", m$config$id), dataset, results)`).
- `dataset` — the same shared `reactiveValues` object from `server.R:7-10` that `mod_dataset.R` writes to. This file **reads** `dataset$expr`/`dataset$meta`/`dataset$source` (to detect and offer an uploaded dataset) and, in exactly one place, **writes** `dataset$expr`/`dataset$source` (§2.7 — live quantile normalization, gated to the uploaded-data case only).
- `results` — passed in but **never referenced anywhere in this file's body** (verified: no `results$` or bare `results` token appears in `mod_overview_server`'s code). This is a real, verifiable gap — see §7, Issue O-9.
- Global data/functions from `global.R`: `GEO_SOURCES`, `get_raw_eset()`, `geo_link()`, `load_individual_dataset()`, `compute_sample_qc()`, `summarize_norm_diagnostics()`, `needs_quantile_norm()`, `pca_of()`, `scree_plot()`, `plot_pca_advanced()`, `qc_bar_plot()`, `theme_arthomix()`, `arthomix_pair()`, `ARTHOMIX_COLORS`, `ARTHOMIX_STATUS`.

### Outputs it produces

- No writes to any shared object except the one described in §2.7 (and even that only fires if the user both (a) has an uploaded dataset staged/active and (b) explicitly clicks "Use this normalised version for every sub-module").
- Everything else is pure, session-local UI: tables, plots, value boxes, downloadable CSVs — nothing here is read by any *other* sub-module. This module's outputs are terminal: a human reads them, no other R code consumes them (reinforced by the `results` gap in Issue O-9 — even the AI assistant can't see what this tab found).

### Downstream consumers

- **None**, functionally — this is this module's single most important architectural property, and the opposite of `mod_dataset.R`'s role. No other `mod_*.R` file reads anything this file computes.
- The one exception: if the "Normalise this dataset" → "Use this normalised version for every sub-module" button is clicked (only reachable when `dataset$source` already starts with `"Uploaded dataset:"`), `dataset$expr` is overwritten in place — at that point, **every** other Transcriptomics sub-module that reads `dataset$expr` is affected, immediately, without going through Preprocessing's "Use this as the active dataset" gate at all. This is a real exception to the "staged vs. active, Preprocessing is the only gate" architecture documented in the Dataset module notes — see §7, Issue O-2, it deserves explicit treatment in your Software Architecture chapter as a **second, narrower activation path**, not an oversight to hide.

### What would happen if this file stopped working

- The app would still fully function for every analysis: DGE, WGCNA, MR, etc. all read `dataset$expr/meta/source` directly and don't need this file. The user would simply lose the ability to browse the raw GEO catalog, look at unfiltered sample metadata/expression matrices, or run the interactive QC checks (missing values, outliers, normalization diagnostics, and the live "normalize now" feature for their own upload).
- Since `mod_overview_server()` is instantiated eagerly at session start (`server.R:141`, same `lapply(TX_MODULES, ...)` as every other sub-module), an error inside its `moduleServer()` setup body would very likely crash the whole Shiny session at startup, not just this tab — same risk profile as every sibling sub-module.

### Conceptual workflow

```
User adds "Overview and Datasets" from the Sub-modules grid (or is routed there by
a sidebar nav click that also pre-fills the grid's search box - see §12)
        |
User picks one of 4 top-level tabs: Datasets / Metadata / Expression data / QC
        |
(Metadata / Expression data / QC each independently ask: "Dataset to inspect?" -
 one of the 4 fixed GEO sources, or "Your uploaded data" if something was staged
 on the Dataset tab)
        |
Server resolves that choice to a plain {expr, meta, label} triple
 (resolve_qc_source(), reading either dataset$expr/meta directly, or
  load_individual_dataset() from global.R - never dataset$staged_*)
        |
Tab-specific view renders: full metadata table / full expression matrix /
 missing-value audit (always-on) / outlier detection (on demand, button) /
 normalization check (on demand, button) / metadata-driven sample filter (on demand)
        |
(QC tab, Normalised data sub-tab only): user can also click "Apply quantile
 normalisation" to see a live before/after comparison, and, only if the
 inspected dataset IS the user's own upload, "adopt" it - which overwrites
 dataset$expr directly, live, for every other sub-module
```

---

## 2. UI — organized by the four actual sub-tabs

`mod_overview_ui(id)` (L28-139) builds one `tabsetPanel(id = ns("tabs"), ...)` with **four `tabPanel`s**, matching this file's own header comment (L2-19) exactly:

1. **"Datasets"** (L32-36)
2. **"Metadata"** (L37-54)
3. **"Expression data"** (L55-70)
4. **"QC"** (L71-137) — itself a **nested** `tabsetPanel(id = ns("qc_tabs"), ...)` with four further tabs: **"Missing values"**, **"Outliers"**, **"Normalised data"**, **"Group"**.

This is the subtab structure the rest of this document follows, per your instruction to organize by the code's actual tabs rather than an invented grouping.

### 2.0 A mechanism shared by three of the four top-level tabs: the independent "Dataset to inspect" picker

Before walking each tab, one design decision recurs three times (Metadata, Expression data, QC) and is worth understanding once. Each of those three tabs has its **own**, independent `selectInput` (`qc_source_meta`, `qc_source_expr`, `qc_source`) offering the same choice set: the four fixed `GEO_SOURCES` entries, plus — only if `dataset$source` currently starts with `"Uploaded dataset:"` — a fifth "Your uploaded data" option. The file's own comment (L205-209) states the reason directly: a single shared picker, shown/hidden per tab via `conditionalPanel`, doesn't reliably bind once this module's UI is inserted dynamically via `insertTab()` (the Sub-modules grid's Add mechanism) — so three self-contained pickers were used instead of one shared one.

**Practical, testable consequence:** picking "GSE15573" on the Metadata tab does **not** carry over to the Expression data or QC tabs — each defaults independently to `GEO_SOURCES[[1]]$gse` (i.e., GSE93272) every time, until the user picks again on that specific tab. This is a genuine UX quirk to verify and document (§7, Issue O-3), not a bug in the sense of producing wrong data — it's an explicit, commented trade-off.

---

### 2.1 Sub-tab: "Datasets" (L32-36)

**UI elements:**

| Element | Type | Static or dynamic |
|---|---|---|
| Intro paragraph ("The example cohort is built from two GEO series...") | `p()`, static text | Static |
| `uiOutput(ns("sources_ui"))`, wrapped in `withSpinner()` | Placeholder for a grid of cards, one per `GEO_SOURCES` entry | **Dynamic** — populated by `output$sources_ui <- renderUI({...})` |

No inputs at all on this tab — it is purely informational, read-only, and requires no user action beyond opening the tab. `withSpinner(..., color = "#2c6fbb", type = 6)` (`shinycssloaders`) shows a themed loading spinner while the underlying `renderUI` computes (relevant here because `get_raw_eset()` reads a `.rds` file from disk per card, which is not instantaneous the first time).

### 2.2 Sub-tab: "Metadata" (L37-54)

| Element | Type | Static/Dynamic |
|---|---|---|
| `box()` containing `qc_source_ui_meta` + `qc_source_info_ui_meta` in a `fluidRow` | Panel | Static box, dynamic contents |
| `uiOutput(ns("qc_source_ui_meta"))` | The "Dataset to inspect" `selectInput` (§2.0) | Dynamic (rendered, not static — see §7 rationale) |
| `uiOutput(ns("qc_source_info_ui_meta"))` | One-line summary ("Label — N samples, M features") | Dynamic |
| Descriptive `p()` | "What's currently selected above, as a whole, before any filtering." | Static |
| `uiOutput(ns("understand_ui"))`, spinner-wrapped | 4 `valueBox`es (Samples / Groups / Features / Sex categories) + a "next steps" hint | Dynamic |
| `box(title = "Sample metadata")` containing a download button + `DT::dataTableOutput(ns("meta_table_full"))` | The full, **unfiltered** metadata table | Dynamic table; static "Download CSV" `downloadButton` |

### 2.3 Sub-tab: "Expression data" (L55-70)

| Element | Type | Static/Dynamic |
|---|---|---|
| `box()` with `qc_source_ui_expr` + `qc_source_info_ui_expr` | Same "Dataset to inspect" pattern, independent instance | Dynamic |
| Descriptive `p()` | "The actual expression matrix behind every check below..." | Static (slightly misleading wording — this text is scoped to the Expression data tab itself, not "every check below," which actually lives on the QC tab; see §7 Issue O-1 for the related mislabeled server-side comment) |
| `box(title = "Expression matrix")` with download button + `DT::dataTableOutput(ns("expr_table"))` | The full expression matrix, features (rows) × samples (columns) | Dynamic table (server-side DT processing — see §4) |

### 2.4 Sub-tab: "QC" — top-level wrapper (L71-79)

Same "Dataset to inspect" picker pattern a third time (`qc_source_ui`/`qc_source_info_ui`), followed by the nested `qc_tabs` tabsetPanel. This top-level picker is what every one of the four QC sub-tabs below actually reads from (`qc_target()`), **not** three separate pickers within QC itself — only Metadata/Expression data/QC have independent pickers; QC's own four inner tabs all share this one.

### 2.5 QC → "Missing values" (L82-93)

| Element | Type |
|---|---|
| Descriptive `p()` | Static |
| `box(title = "Percent missing by field")` with `plotOutput(ns("missing_plot"), height = 320)`, spinner-wrapped | Dynamic bar chart |
| `box(title = "By field")` with `DT::dataTableOutput(ns("missing_table"))` | Dynamic table |

No inputs — this check is **always on**, computed the moment the tab is opened (it depends only on `qc_target()`, not on any button).

### 2.6 QC → "Outliers" (L94-106)

| Element | Type |
|---|---|
| Descriptive `p()` | Static |
| `sliderInput(ns("mad_k"), "Outlier sensitivity...", min=2, max=6, value=3, step=0.5)` | Static input, dynamic value |
| `actionButton(ns("run_qc_btn"), "Run outlier detection", icon("play"))` | Static input |
| `uiOutput(ns("qc_summary_ui"))`, spinner-wrapped | Dynamic — a single `valueBox` (flagged-sample count) once run |
| `uiOutput(ns("qc_plots_ui"))`, spinner-wrapped | Dynamic — 3 small plots (signal/detected/correlation) + a flagged-samples table with its own download button, all appearing only after the button is clicked |

This is the first of two **on-demand, button-gated** computations in the module (the other is Normalised data). Nothing runs until `run_qc_btn` is clicked at least once.

### 2.7 QC → "Normalised data" (L107-126)

| Element | Type |
|---|---|
| Descriptive `p()` | Static |
| `uiOutput(ns("norm_color_by_ui"))` | Dynamic `selectInput` — "Color by", choices = every metadata column except `sample` |
| Two `selectInput`s side by side, `norm_pc_x`/`norm_pc_y` | Static inputs (choices fixed at PC1–PC5) |
| `checkboxInput(ns("norm_show_ellipse"), ..., value = TRUE)` | Static input |
| `checkboxInput(ns("norm_show_labels"), ..., value = FALSE)` | Static input |
| `actionButton(ns("run_norm_btn"), "Run normalisation check")` | Static input |
| `uiOutput(ns("norm_summary_ui"))`, spinner-wrapped | Dynamic — plain-language verdict text |
| `uiOutput(ns("norm_views_ui"))`, spinner-wrapped | Dynamic — per-sample boxplot, scree plot, PCA scatter + per-sample summary table |
| `uiOutput(ns("norm_apply_ui"))` | Dynamic — the "Normalise this dataset" sub-panel (its own action button, `apply_norm_btn`, plus, conditionally, `adopt_norm_btn`) |

This sub-tab has the module's only genuinely *mutating* control path (§1, "Outputs it produces"): clicking "Apply quantile normalisation" runs a live `limma::normalizeBetweenArrays()` call and shows a before/after comparison; if (and only if) the dataset currently under inspection is the user's own upload (`identical(input$qc_source, "uploaded")`), a further button, "Use this normalised version for every sub-module," appears and — if clicked — overwrites `dataset$expr` directly.

### 2.8 QC → "Group" (L127-135)

| Element | Type |
|---|---|
| Descriptive `p()` | Static |
| `uiOutput(ns("filters"))` | Dynamic — one input per filterable metadata column (a `pickerInput` multi-select for categorical columns, a `sliderInput` range for numeric columns), plus "Apply filters"/"Reset" buttons |
| `uiOutput(ns("filter_summary_ui"))` | Dynamic — value boxes (Samples selected / Groups / Datasets / Sex categories) |
| `uiOutput(ns("filtered_views_ui"))` | Dynamic — composition bar chart(s) + a filtered, still fully sortable/filterable metadata table with its own download button |

Also on-demand/button-gated (`apply_btn`), but with a `reset_btn` companion the other two on-demand tabs don't have.

---

## 3. FUNCTIONALITY (by logical section, in file order)

1. **Datasets tab server** (L147-175) — one `renderUI`, no reactivity beyond re-running if `GEO_SOURCES` itself changed (it never does at runtime — this block effectively runs once per session and is static in practice).
2. **The three independent "Dataset to inspect" picker/resolver trios** (L177-244) — `qc_source_choices()` (shared reactive), `resolve_qc_source()` (plain function, not reactive itself — called from inside three separate `reactive()`s), and three parallel `{ui, target, info}` triples, one per tab (`_meta`, `_expr`, un-suffixed for QC).
3. **Metadata tab server** (L248-272) — `understand_ui`, `meta_table_full`, `download_meta_full`; all depend on `qc_target_meta()`.
4. **QC → Missing values** (L276-302) — `missing_audit` (a `reactive()`, always-on), `missing_plot`, `missing_table`.
5. **QC → Outliers** (L306-365) — `sample_qc` (an `eventReactive` gated on `run_qc_btn`), plus its rendered summary/plots/table/download.
6. **QC → Normalised data — the check** (L369-462) — `norm_color_by_ui`, `norm_check` (an `eventReactive` gated on `run_norm_btn`), plus rendered summary/plots/table.
7. **QC → Normalised data — apply/adopt** (L464-546) — `norm_apply_result` (`eventReactive` on `apply_norm_btn`), rendered before/after comparison, and the conditional `adopt_norm_btn` `observeEvent` that is this file's one mutation of shared state.
8. **Expression data tab server** (L552-569) — `expr_table_data`, `expr_table` (server-side DT), `download_expr`. *(Physically located after the Normalise-and-Apply block in the file, under a comment that says "QC tab: browse the expression matrix" — see §7, Issue O-1: this actually feeds the "Expression data" tab, not "QC".)*
9. **QC → Group** (L573-711) — `filter_spec` (a `reactive()` that inspects `qc_target()$meta`'s columns generically), `filters` renderUI, `reset_btn` observer, `filtered_meta` (`eventReactive` on `apply_btn`), and its rendered summary/plots/table/download.

---

## 4. LINE-BY-LINE TEACHING

### Block A — Header comment (L1-19)

Documents the module's scope precisely: four tabs, each with its own dataset picker, explicitly **not** the merged/batch-corrected working dataset (that's Preprocessing's own view) except for the one uploaded-data exception. Worth quoting in your Implementation chapter almost verbatim as the module's own stated contract — it is, on inspection, accurate to what the code actually does.

### Block B — `mod_overview_config` (L21-26)

```r
mod_overview_config <- list(
  id = "overview", group = "Data",
  title = "Overview and Datasets",
  description = "...",
  icon = "table-cells"
)
```
Same shape as `mod_dataset_config` (see companion notes, Block B) — one of 16 identically-shaped config lists that `TX_MODULES` (`submodules_registry.R:9`) assembles into a single registry list, consumed by the Sub-modules grid (card title/description/icon) and by `jump_to_submodule()` (`server.R:335-346`, reads `cfg$title` to `updateTabsetPanel(session, "tx_menu", selected = cfg$title)`). The extra `group = "Data"` field (absent from `mod_dataset_config`) is **not read anywhere in this file** — grep the rest of the app if you need to confirm what, if anything, groups sub-modules by this field; **cannot confirm from this file alone**.

### Block C — `mod_overview_ui()`, "Datasets" tab (L28-36)

```r
mod_overview_ui <- function(id) {
  ns <- NS(id)
  tabsetPanel(
    id = ns("tabs"), type = "tabs",
    tabPanel(
      "Datasets", br(),
      p(class = "submodule-desc", "..."),
      withSpinner(uiOutput(ns("sources_ui")), color = "#2c6fbb", type = 6)
    ),
```
`ns <- NS(id)` — identical namespacing mechanism to `mod_dataset.R` (see that file's notes for the full explanation of `NS()`/`ns()`). `tabsetPanel(id = ns("tabs"), type = "tabs", ...)` — the outer container; giving it an `id` (rather than leaving it anonymous) is what lets server-side code later call `updateTabsetPanel(session, "tabs", selected = ...)` if needed (not actually done anywhere in this file, but is done externally: `jump_to_submodule()`'s `inner_tab` parameter targets exactly this kind of namespaced tabset ID pattern for *other* sub-modules — worth noting this module's own inner tabs are addressable the same way even though nothing currently addresses them). `withSpinner(uiOutput(...), color=, type=)` — `shinycssloaders`'s wrapper: shows a themed spinner (`type = 6` selects one of several built-in animation styles) over the placeholder until the wrapped output's first render completes; purely UX, no data implication.

### Block D — "Metadata" tab (L37-54)

```r
    tabPanel(
      "Metadata", br(),
      box(
        width = 12, status = "primary", solidHeader = FALSE,
        fluidRow(
          column(5, uiOutput(ns("qc_source_ui_meta"))),
          column(7, uiOutput(ns("qc_source_info_ui_meta")))
        )
      ),
```
`box(width = 12, status = "primary", solidHeader = FALSE, ...)` — a `shinydashboard`-style panel spanning the full 12-column grid width; `solidHeader = FALSE` means the title bar (there isn't one here — no `title=` argument) would be outline-styled rather than filled, moot since no title is given for this particular box. Two `uiOutput`s side by side in a 5/7 grid split: the picker itself gets less width than its own info line, a minor layout choice. Both are placeholders resolved by server code discussed in Block I.

```r
      p(class = "submodule-desc", "What's currently selected above, as a whole, before any filtering."),
      withSpinner(uiOutput(ns("understand_ui")), color = "#2c6fbb", type = 6),
      box(
        width = 12, title = "Sample metadata", status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "Every column is sortable and filterable (click a header, or use the box underneath it)."),
        div(class = "table-toolbar", downloadButton(ns("download_meta_full"), "Download CSV", class = "btn-sm")),
        withSpinner(DT::dataTableOutput(ns("meta_table_full")), color = "#2c6fbb", type = 6)
      )
    ),
```
`downloadButton(ns("download_meta_full"), "Download CSV", class = "btn-sm")` — a special Shiny input that, when clicked, triggers a browser file download by hitting a server-registered download handler (paired with `output$download_meta_full <- downloadHandler(...)` in the server, Block J) rather than firing a normal reactive event — no `input$download_meta_full` value is ever meaningfully read in server code; the button's *click* directly initiates an HTTP request Shiny intercepts and serves the `content` function's output as a file. `DT::dataTableOutput(ns("meta_table_full"))` — the placeholder for a `DT`-rendered (DataTables.js-backed) interactive table: sortable/searchable columns, pagination, per-column filter boxes (enabled server-side via `filter = "top"`, Block J) — a materially richer widget than a plain `tableOutput()`.

### Block E — "Expression data" tab (L55-70)

Structurally identical pattern to Metadata (own picker pair, descriptive text, a box with a download button and a `DT` table) — see the annotated comment on wording accuracy in §2.3/§7 Issue O-1.

### Block F — "QC" tab: outer picker + nested tabsetPanel (L71-138)

```r
    tabPanel(
      "QC", br(),
      box(
        width = 12, status = "primary", solidHeader = FALSE,
        fluidRow(
          column(5, uiOutput(ns("qc_source_ui"))),
          column(7, uiOutput(ns("qc_source_info_ui")))
        )
      ),
      tabsetPanel(
        id = ns("qc_tabs"), type = "tabs",
```
A **second, nested** `tabsetPanel`, given its own `id = ns("qc_tabs")` — this is what lets `jump_to_submodule("overview", inner_tab = "Datasets", ...)` (`server.R:340`, `paste0("tx_", mod_id, "-tabs")`) address a *specific inner tab*, though note the actual sidebar wiring (`server.R:352`) targets the **outer** `tabs` tabset's "Datasets" panel (the top-level sub-tab), not `qc_tabs` — `qc_tabs`'s four children (Missing values/Outliers/Normalised data/Group) have no corresponding sidebar shortcut in this codebase as inspected; **cannot confirm from this file alone** whether anything else targets `qc_tabs` specifically.

**Missing values** (L82-93), **Outliers** (L94-106), **Normalised data** (L107-126), **Group** (L127-135) — each already itemized in §2.5–§2.8; their UI-construction R syntax (`sliderInput`, `actionButton`, `checkboxInput`, `fluidRow`/`column` layout) is standard and mirrors patterns already explained for `mod_dataset.R`; the sections below focus on what's new — the *server* logic behind each.

---

### Block G — `mod_overview_server()` setup (L141-144)

```r
mod_overview_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
```
Same `moduleServer()`/`session$ns` pattern as `mod_dataset_server()`. Note the third parameter, `results = NULL`, given a default — meaning this function *can* be called with only two arguments — but `server.R:141`'s registry loop always passes three (`m$server(paste0("tx_", m$config$id), dataset, results)`) for every `TX_MODULES` entry uniformly, so the default is never actually exercised in this app; it exists purely so this function's *signature* matches every other analysis sub-module's `function(id, dataset, results)` shape, for the registry's generic `lapply()` call to work without special-casing. **`results` is dead weight in this specific file's body** — see Issue O-9.

### Block H — Datasets tab: `sources_ui` (L147-175)

```r
    output$sources_ui <- renderUI({
      cards <- lapply(GEO_SOURCES, function(s) {
        eset <- get_raw_eset(s$gse)
        div(
          class = "info-card",
          div(
            class = "module-card-title-row",
            h4(s$gse),
            tags$a(href = geo_link(s$gse), target = "_blank", rel = "noopener",
                    icon("up-right-from-square"), " NCBI GEO")
          ),
```
`lapply(GEO_SOURCES, function(s) {...})` — iterates the four-element `GEO_SOURCES` list from `global.R:694-699`, each `s` a list like `list(gse="GSE93272", role="Training (whole blood)", used_in="Merged into the example cohort")`. `get_raw_eset(s$gse)` — `global.R:201-211`: reads (and in-memory-caches, `.arthomix_cache`) `data/preloaded/transcriptomics/raw/<gse>_raw.rds` if present, returning `NULL` if the file doesn't exist on this deployment (this is the graceful-degradation path exercised for GSE93272/GSE110169 per the Dataset-module notes — their raw per-source `.rds` files were not part of the self-contained data migration, per that file's own comment chain). `tags$a(href = geo_link(s$gse), target = "_blank", rel = "noopener", ...)` — a genuine external hyperlink (not a Shiny input) to the live NCBI GEO record; `geo_link()` (`global.R:700`) is a one-line string-builder: `paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gse)`. `rel = "noopener"` is a standard security/perf hardening attribute for `target="_blank"` links (prevents the opened page from gaining a JS handle back to the opener window) — a small but genuinely correct detail.

```r
          if (!is.null(eset)) {
            tagList(
              p(class = "module-card-tagline",
                tryCatch(Biobase::experimentData(eset)@title, error = function(e) NULL)),
              p(strong("Role: "), s$role, br(), strong("Used for: "), s$used_in),
              p(strong("Platform: "), Biobase::annotation(eset), br(),
                strong("Samples: "), ncol(eset), ", ", strong("Probes: "), format(nrow(eset), big.mark = ","))
            )
          } else {
            tagList(
              p(strong("Role: "), s$role, br(), strong("Used for: "), s$used_in),
              div(class = "empty-note", icon("triangle-exclamation"), "Raw file not found on disk.")
            )
          }
        )
      })
      div(class = "module-grid", cards)
    })
```
Two-branch card body: if the raw `ExpressionSet` was found, show its title (`Biobase::experimentData(eset)@title` — an S4 slot access via `@`, reading the GEO series' own submitted title string, wrapped in a `tryCatch` because some `ExpressionSet`s may have empty/absent `experimentData`, in which case `NULL` is returned and `p()` simply renders an empty tag rather than erroring — a small, correct defensive touch), role/used-in text (from the static `GEO_SOURCES` list, not the eset), platform (`Biobase::annotation(eset)`, e.g. `"GPL570"`), and dimensions (`ncol(eset)` = samples, `nrow(eset)` = probes, since an `ExpressionSet`'s dimensions follow the same features-as-rows/samples-as-columns convention as a plain matrix). If not found, a shorter card with just role/used-in and a visible warning note — **this is the one place in the app that visibly and honestly tells the user "this raw file isn't available in this deployment,"** rather than silently omitting the card or crashing; a genuinely good transparency/reproducibility design choice worth citing.

**Reactivity note:** this `renderUI` reads `GEO_SOURCES` (a static, file-scope constant, never reassigned at runtime) and calls `get_raw_eset()` (whose internal cache means repeat calls are cheap) — there is no reactive input anywhere in this block, so despite being wrapped in `renderUI`, it effectively renders exactly once per session and never re-runs. This is a legitimate use of `renderUI` even for static-per-session content, since the underlying values (whether each `.rds` exists, its title/platform/dims) are only knowable at runtime, not at UI-definition time.

### Block I — The shared "Dataset to inspect" mechanism (L177-244)

```r
    qc_source_choices <- reactive({
      choices <- setNames(vapply(GEO_SOURCES, `[[`, character(1), "gse"),
                            vapply(GEO_SOURCES, function(s) sprintf("%s (%s, raw)", s$gse, s$role), character(1)))
      if (grepl("^Uploaded dataset:", dataset$source %||% "")) {
        choices <- c(setNames("uploaded", paste0("Your uploaded data (", dataset$source, ")")), choices)
      }
      choices
    })
```
A `reactive()` (not per-tab — genuinely shared by the three picker-rendering blocks below), producing a named vector: values are GSE accessions (or the literal string `"uploaded"`), names are display labels. `vapply(GEO_SOURCES, `[[`, character(1), "gse")` — same `` `[[` ``-as-a-function idiom explained in the Dataset-module notes (Block C's `vapply` explanation) — extracts all four `gse` fields at once. `grepl("^Uploaded dataset:", dataset$source %||% "")` — checks whether the shared `dataset$source` string currently starts (`^`) with the literal prefix `"Uploaded dataset:"`; this exact prefix is set in exactly one place in the whole app (`mod_dataset.R:286`, the upload-commit handler, once promoted to active by Preprocessing) — a **string-literal contract between two files with no shared constant**, worth flagging as a maintainability point (§7, Issue O-4) even though it's currently consistent. `dataset$source %||% ""` — defensive: before any dataset is ever loaded this could theoretically be `NULL` (it isn't in practice, since `server.R:7-10` initializes `dataset$source` via `load_default_dataset()` at app startup — so this is a defensive no-op in the current app, exactly analogous to the `%||%` use already seen in `mod_dataset.R`'s GEO-accession handling).

**Important, code-verified scope limitation:** this check only recognizes a dataset staged via **upload** (`mod_dataset.R`'s own file-upload path) as "uploaded" — a dataset currently active via the **preloaded-dataset** or **live-GEO-fetch** paths of `mod_dataset.R` does *not* match `"^Uploaded dataset:"` (those set `dataset$source` to a `load_default_dataset()`/`load_individual_dataset()` string or `"NCBI GEO: ..."` respectively — see the Dataset-module notes' Path 1/3 traces) and so will **not** appear as a 5th "Your uploaded data" choice here, even though the app's own active dataset genuinely differs from all four fixed `GEO_SOURCES`. This is a real, testable gap — see §7, Issue O-5.

```r
    resolve_qc_source <- function(source_id) {
      if (identical(source_id, "uploaded")) {
        req(dataset$expr, dataset$meta)
        return(list(expr = dataset$expr, meta = dataset$meta, label = paste0("Your uploaded data: ", dataset$source)))
      }
      d <- load_individual_dataset(source_id)
      validate(need(!is.null(d), paste("Raw data for", source_id, "was not found on disk.")))
      d
    }
```
A **plain function** (not itself a `reactive()`), called from inside three separate `reactive()` bodies below — each call re-evaluates it fresh; there is no shared memoization across the three tabs' resolvers (a mild inefficiency if the same GSE is picked on two tabs simultaneously — `load_individual_dataset()` itself is cached at the `global.R` level via `.arthomix_cache`, so the actual file I/O is not repeated, just this thin wrapper). `identical(source_id, "uploaded")` — again `identical()` over `==`, correct idiom (see Dataset-module notes' rationale). Uploaded branch: `req(dataset$expr, dataset$meta)` blocks if either is `NULL` (shouldn't happen given `qc_source_choices()` only offers `"uploaded"` when `dataset$source` already matches the upload prefix — but note `dataset$expr`/`dataset$meta` are the **active**, not staged, fields — i.e. this reads whatever `mod_preprocessing.R` most recently promoted, not whatever's merely staged/previewed on the Dataset tab — an important, correct distinction, consistent with this module's stated "never reads `dataset$staged_*`" design). Otherwise: `load_individual_dataset(source_id)` (`global.R:772-811` — reads the raw `ExpressionSet` or, for GSE89408, a separate counts file, or falls back to `merged_training_subset()` for GSE93272/GSE110169 when the raw file is absent), wrapped in `validate(need(!is.null(d), ...))` — the same Shiny "stop gracefully with a message" idiom seen throughout `mod_dataset.R`.

```r
    output$qc_source_ui_meta <- renderUI({
      selectInput(ns("qc_source_meta"), "Dataset to inspect", choices = qc_source_choices(),
                  selected = GEO_SOURCES[[1]]$gse, width = "100%")
    })
    qc_target_meta <- reactive({ req(input$qc_source_meta); resolve_qc_source(input$qc_source_meta) })
    output$qc_source_info_ui_meta <- renderUI({
      t <- qc_target_meta()
      div(class = "empty-note", icon("circle-info"), strong(t$label), " - ",
          format(nrow(t$meta), big.mark = ","), " samples, ",
          format(nrow(t$expr), big.mark = ","), " features.")
    })
```
This exact three-output pattern (`qc_source_ui_*` → `qc_target_*` → `qc_source_info_ui_*`) is written out **three times**, verbatim except for the `_meta`/`_expr`/(none) suffix — for Metadata, Expression data, and QC respectively. `selectInput(...)` is itself wrapped in `renderUI()` (rather than being static UI built once in `mod_overview_ui()`) — the file's own comment explains why (§2.0): eager static UI generation would need `choices` to be fixed at UI-definition time, but the "uploaded" 5th choice can only be known at runtime, and — per the comment — a `selectInput(choices=NULL)` + later `updateSelectInput()` pattern is unreliable here specifically because this module's server runs at session start, *before* the user has added its tab via the Sub-modules grid, so a fire-once `update*` call sent to a `<select>` that doesn't exist in the DOM yet is silently lost, permanently leaving the dropdown empty. Rendering the `selectInput` itself sidesteps this because `renderUI` only actually executes once its output is bound (i.e., once the tab is visible), guaranteeing the choices are populated whenever a user can actually see the control. **This is a subtle, correct, and non-obvious Shiny lesson** worth explaining fully in your Implementation chapter — it explains why this module (and, by the same architectural cause, several controls elsewhere in this app) prefers "always re-render the whole input" over "render once, then patch."

`selected = GEO_SOURCES[[1]]$gse` — always defaults to GSE93272 (index 1) specifically, regardless of what's currently active app-wide or what was picked on another tab; this is the concrete mechanism behind the §2.0 "doesn't remember across tabs" behavior. `qc_target_meta <- reactive({ req(input$qc_source_meta); resolve_qc_source(input$qc_source_meta) })` — the resolved `{expr, meta, label}` triple, recomputed any time the picker's value changes; this is what every other Metadata-tab output reads. The info-line `renderUI` formats sample/feature counts with `format(..., big.mark = ",")` for readability (e.g., `"15,763"`).

---

### Block J — Metadata tab server: `understand_ui`, `meta_table_full`, download (L246-272)

```r
    output$understand_ui <- renderUI({
      t <- qc_target_meta()
      meta <- t$meta
      tagList(
        p(strong("Source: "), t$label),
        fluidRow(
          valueBox(nrow(meta), "Samples", icon = icon("users"), color = "light-blue", width = 3),
          valueBox(length(unique(na.omit(meta$group))), "Groups", icon = icon("layer-group"), color = "purple", width = 3),
          valueBox(format(nrow(t$expr), big.mark = ","), "Features in matrix", icon = icon("dna"), color = "green", width = 3),
          valueBox(length(unique(na.omit(meta$sex))), "Sex categories", icon = icon("venus-mars"), color = "red", width = 3)
        ),
        p(class = "submodule-desc", "Next: browse the full metadata table below or the Expression data tab...")
      )
    })
```
`nrow(meta)` — sample count, directly from the metadata data frame's row count (**not** cross-checked against `ncol(t$expr)** — see §7, Issue O-6: for a dataset where `meta` and `expr` could in principle disagree in sample count, this value box would silently report the metadata's count only). `length(unique(na.omit(meta$group)))` — counts **distinct, non-missing** group labels: `na.omit()` strips `NA`s first, `unique()` then dedupes, `length()` counts — a correct, standard R idiom for "how many real categories are present," immune to `NA` inflating the count as a spurious extra "category." Same pattern for `meta$sex`. `valueBox(...)` — a `shinydashboard` big-number KPI tile; `color=`/`icon=`/`width=` are purely presentational (the `width` values sum to 12, filling the row exactly, per Bootstrap's grid).

`assumes `meta$group` and `meta$sex` columns exist** — true for every dataset this app can route through `resolve_qc_source()` (the four `load_individual_dataset()` branches all guarantee `group`/`sex` columns per `eset_harmonize_meta()`/`GSE89408`'s bespoke construction in `global.R`, and the uploaded-data path guarantees them via `mod_dataset.R`'s column-mapping harmonization) — so this is a safe assumption given the module's actual reachable inputs, not an unchecked one; worth stating this explicitly as a **verified invariant**, not "unchecked and therefore risky."

```r
    output$meta_table_full <- DT::renderDataTable({
      DT::datatable(qc_target_meta()$meta, rownames = FALSE, filter = "top",
                     options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_meta_full <- downloadHandler(
      filename = function() paste0(qc_target_meta()$meta$dataset[1] %||% "metadata", "_metadata.csv"),
      content = function(file) write.csv(qc_target_meta()$meta, file, row.names = FALSE)
    )
```
`DT::datatable(..., rownames = FALSE, filter = "top", options = list(pageLength = 10, scrollX = TRUE))` — `filter = "top"` adds a per-column filter input row directly under the header (this is what "or use the box underneath it," from the UI's own descriptive text, refers to); `rownames = FALSE` hides R's own row index column since it carries no information here (the `sample` column already identifies rows); `scrollX = TRUE` allows horizontal scrolling for a metadata frame with many columns rather than squeezing/wrapping. `class = "stripe hover compact"` are `DataTables`/Bootstrap CSS classes for visual styling only.

**Correction to an earlier draft of this document:** this call does **not** pass `server=` at all — but `DT::renderDataTable()`'s own default for that argument is `server = TRUE` (verified directly: `args(DT::renderDataTable)` shows `server = TRUE` as the default), not `FALSE`. So this table is, in fact, **already server-side processed**, exactly like the Expression data tab's `expr_table` (Block O) — there is **no actual client-side table anywhere in this module**. The only difference between this call and `expr_table`'s is that `expr_table` states `server = TRUE` explicitly (redundant with the default, but self-documenting given how large that specific matrix can get), while every other `DT::renderDataTable()` call in this file — this one, `missing_table`, `qc_table`, `norm_table`, and the filtered `meta_table` (Group tab) — relies on the identical default silently. This is a genuinely easy mistake to make reading this code casually (it certainly was, in an earlier pass of this very document), and is worth stating explicitly as a validation note: **do not infer client- vs. server-side DT behavior from the mere presence or absence of `server = TRUE` in the call** — check the package's actual default first.

`downloadHandler(filename = function() ..., content = function(file) ...)` — `filename` is evaluated once, when the download begins, to name the file the browser saves; `content` receives a temp file path (`file`) and must write the actual bytes to it (here, `write.csv(..., row.names = FALSE)`) — Shiny then streams that temp file to the browser. `qc_target_meta()$meta$dataset[1] %||% "metadata"` — attempts to name the download after the value in the metadata's own `dataset` column (first row), falling back to the literal string `"metadata"` if that column doesn't exist or is `NA` (`%||%` only substitutes on `NULL`, not `NA` — see §7, Issue O-7 for why this fallback may not always fire when you'd expect).

---

### Block K — QC → "Missing values" (L276-302)

```r
    missing_audit <- reactive({
      df <- qc_target()$meta
      miss <- data.frame(
        field = colnames(df),
        n_missing = vapply(df, function(x) sum(is.na(x) | x %in% c("", "NA", "N/A", "unknown")), integer(1)),
        stringsAsFactors = FALSE
      )
      miss$pct_missing <- round(100 * miss$n_missing / nrow(df), 1)
      miss$status <- ifelse(miss$n_missing == 0, "Complete", "Has missing")
      miss[order(miss$pct_missing), ]
    })
```
An **always-on** `reactive()` — no button, no `eventReactive`; it recomputes any time `qc_target()` changes (i.e., whenever the QC tab's own "Dataset to inspect" picker changes). `vapply(df, function(x) sum(is.na(x) | x %in% c("", "NA", "N/A", "unknown")), integer(1))` — iterates every **column** of the metadata data frame (`vapply` over a data frame iterates its columns, since a data frame is internally a list of equal-length vectors); for each column `x`, counts entries that are either a true `NA` **or** one of four literal sentinel strings often used to encode missingness in raw/exported data (`""`, `"NA"` as text, `"N/A"`, `"unknown"`) — a genuinely good practice, since a naive `sum(is.na(x))` alone would miss all four of those common "missing, but not R's `NA`" encodings, which are exactly the kind of thing that slips into GEO-derived phenotype columns. **Caveat worth stating explicitly:** `x %in% c(...)` performs an implicit `as.character()` coercion on `x` if it's not already character (`%in%` compares element-wise after coercing both sides to a common type) — for a genuinely numeric column, this comparison is harmless (no numeric value literally equals the string `"unknown"` after coercion) but slightly wasteful; not a correctness bug.

`round(100 * miss$n_missing / nrow(df), 1)` — percent missing per field, one decimal place. `ifelse(miss$n_missing == 0, "Complete", "Has missing")` — `ifelse()` is R's **vectorized** conditional (distinct from scalar `if`/`else` used elsewhere in this file) — appropriate here since it must classify every row (field) of `miss` at once, not a single value. `miss[order(miss$pct_missing), ]` — sorts ascending by percent missing (cleanest fields first) — feeds directly into the bar chart below, which further reorders by `reorder(field, pct_missing)` (so the sort here is largely redundant with that `reorder()` call for the plot specifically, but does matter for the *table*, which has no such reordering of its own beyond this).

```r
    output$missing_plot <- renderPlot({
      miss <- missing_audit()
      ggplot(miss, aes(x = reorder(field, pct_missing), y = pct_missing, fill = status)) +
        geom_col(width = 0.7) +
        scale_fill_manual(values = c("Complete" = ARTHOMIX_STATUS$good, "Has missing" = ARTHOMIX_STATUS$critical)) +
        coord_flip() +
        labs(x = NULL, y = "% missing", fill = NULL) +
        theme_arthomix(base_size = 13) +
        theme(panel.grid.major.y = element_blank())
    })
```
Standard `ggplot2` horizontal bar chart: `reorder(field, pct_missing)` sorts the categorical `field` axis by the numeric `pct_missing` value (ascending) rather than alphabetically; `coord_flip()` turns the resulting vertical bar chart horizontal (so long field names read left-to-right rather than being rotated); `scale_fill_manual(values = c("Complete" = ARTHOMIX_STATUS$good, "Has missing" = ARTHOMIX_STATUS$critical))` — maps the two `status` categories to this app's fixed semantic colors (`global.R:1410`: green `#0ca30c` for good, red `#d03b3b` for critical) — **not** the general-purpose categorical palette (`arthomix_pair()`), correctly, since "good vs. bad" is a status distinction, not an arbitrary category, and the file's own comment on `ARTHOMIX_STATUS` (`global.R:1409`) states it's "reserved for state, never reused for an ordinary series" — this usage respects that convention exactly. `theme_arthomix(base_size = 13)` — this app's shared `ggplot2` theme (`global.R:1423-1440`): minimal grid, muted axis text, bottom legend, consistent typography across every plot in the app; `theme(panel.grid.major.y = element_blank())` removes the (now-redundant, since bars are the y-categories) horizontal gridlines specifically for this flipped-coordinate plot.

`output$missing_table` — a second, smaller `DT::datatable` over the same `missing_audit()` data, restricted to 4 display columns, `dom = "tp"` (DataTables' layout-control string: only show the **t**able and **p**agination controls, omitting the search box/length selector/info line that the default `dom` string would include — appropriate for a small, fixed-size summary table where those extra controls add clutter, not value).

---

### Block L — QC → "Outliers" (L306-365)

```r
    sample_qc <- eventReactive(input$run_qc_btn, {
      t <- qc_target()
      qc <- compute_sample_qc(t$expr, mad_k = input$mad_k)
      merge(qc, t$meta[, intersect(c("sample", "group"), colnames(t$meta))], by = "sample", all.x = TRUE)
    })
```
`eventReactive(input$run_qc_btn, {...})` — the click-gated value-returning construct (same as `mod_dataset.R`'s `geo_fetch_result`), correctly chosen here: reading `input$mad_k` (the sensitivity slider) alone must **not** trigger recomputation on every drag — only an actual click of "Run outlier detection" should. `compute_sample_qc(t$expr, mad_k = input$mad_k)` — `global.R:1453-1479`: computes three per-sample QC metrics —

1. `signal <- colSums(expr, na.rm = TRUE)` — total expression signal per sample (sum over all features).
2. `detected <- colSums(expr > detect_cutoff, na.rm = TRUE)` — count of features above the dataset's own 25th-percentile value (`detect_cutoff <- stats::quantile(expr, 0.25, na.rm = TRUE)`) — a relative, data-driven "detected" threshold rather than a fixed absolute cutoff, which is what lets this one function work whether `expr` holds raw counts, log-CPM, or log2 microarray intensities (per its own comment).
3. `mean_cor` — mean pairwise correlation of each sample against every other sample, computed only over the top-`top_n_cor` (default 2000) most-variable features (`gene_var <- apply(expr, 1, stats::var, ...)`, `top_idx <- order(gene_var, decreasing = TRUE)[seq_len(min(top_n_cor, nrow(expr)))]`) — restricting to high-variance features is both a performance optimization (correlation over the full feature set on a large matrix is expensive) and, more importantly, a *signal* choice: low-variance (near-flat, uninformative) features would dilute the correlation structure that actually distinguishes technical outliers from the cohort.

Each metric is then flagged via `is_outlier(x, low_only=FALSE)`: `m <- median(x); s <- mad(x); abs(x - m) > mad_k * s` — a **robust** (median/MAD-based, not mean/SD-based) outlier rule, deliberately chosen because MAD is far less sensitive to the very outliers you're trying to detect than SD would be (an SD-based rule computed *including* the outliers can itself be inflated by them, weakening the test — a well-known robust-statistics rationale worth citing explicitly in your Methods). `flag_cor` uses `low_only = TRUE` (`(m - x) > mad_k * s`) — correctly one-sided: a sample correlating *unusually highly* with the rest of the cohort is not a QC concern the way unusually *low* correlation is (a sign of a technical outlier not resembling anything else in the batch), so only the low tail is flagged for this one metric, while `flag_signal`/`flag_detected` are flagged both-sided (`abs(...)`) — both unusually low **and** unusually high total signal/detected-feature counts are potential technical problems. `if (s == 0) return(rep(FALSE, length(x)))` — guards against a degenerate all-identical-values case where `mad()` is exactly zero, which would otherwise make **every** non-identical value "infinitely" far from the median in MAD units, flagging the entire cohort — correctly short-circuits to "nothing flagged" instead.

`merge(qc, t$meta[, intersect(c("sample","group"), colnames(t$meta))], by="sample", all.x=TRUE)` — joins the per-sample QC metrics back to the metadata's `sample`/`group` columns (only `group`, not the full metadata, and only if `group` actually exists — `intersect()` guards against a metadata frame missing it, though per §7 Issue O-6/J's own note, every reachable metadata source does have it). `all.x = TRUE` is a **left join**: every row of `qc` (one per expression-matrix sample) is kept even if no metadata match exists — for a sample present in `expr` but absent from `meta` (shouldn't happen given how `resolve_qc_source()` builds its inputs, but a `merge()` without `all.x=TRUE` would silently *drop* such a sample from the result instead of keeping it with `NA` group, a materially different and worse failure mode).

```r
    output$qc_summary_ui <- renderUI({
      if (!isTruthy(input$run_qc_btn) || input$run_qc_btn == 0) {
        return(div(class = "empty-note", icon("circle-info"), "Not run yet - click Run outlier detection to check for technical outliers."))
      }
      qc <- sample_qc()
      n_flagged <- sum(qc$flag_signal | qc$flag_detected | qc$flag_cor)
      valueBox(n_flagged, "Samples flagged", icon = icon("triangle-exclamation"),
                color = if (n_flagged > 0) "red" else "green", width = 12)
    })
```
`!isTruthy(input$run_qc_btn) || input$run_qc_btn == 0` — an explicit pre-click guard, checked *before* calling `sample_qc()` at all: `isTruthy()` is Shiny's own truthiness test (the same one `req()` uses internally) — an `actionButton`'s value is an integer counter starting at `0` before any click, and `isTruthy(0)` is `FALSE`, so in practice the `|| input$run_qc_btn == 0` clause is redundant with the first (both branches would catch the pre-click state identically) — a small, harmless bit of defensive redundancy, not a bug. `n_flagged <- sum(qc$flag_signal | qc$flag_detected | qc$flag_cor)` — `|` here is the **vectorized** logical OR (correct, since these are equal-length logical vectors, one entry per sample) — counts samples flagged by **any** of the three criteria (a sample flagged on two or three criteria simultaneously is still counted once, not double-counted — correct for "how many distinct samples were flagged," though note this single number can't distinguish "3 samples each flagged once" from "1 sample flagged on all three criteria," a nuance the flagged-samples table below (which shows a `reason` column) does resolve).

`output$qc_plots_ui` (L322-339) — a `renderUI` that only actually appears once `sample_qc()` has a value (guarded by `req(sample_qc())`, which itself won't resolve until the button's been clicked, since `sample_qc` is an `eventReactive`) — three side-by-side small boxed plots plus a flagged-samples table with its own download button. `output$signal_plot <- renderPlot(qc_bar_plot(sample_qc(), "signal", "flag_signal", "Total signal"))` and its two siblings all delegate to the single shared `qc_bar_plot()` helper (`global.R:1952-1960`) — a thin, reusable bar chart: neutral blue for in-range samples, status red for flagged ones, x-axis labels suppressed (illegible past ~30 samples per its own comment) — the same "shape of the distribution plus which bars are red" idea used throughout this app's QC visualizations.

```r
    qc_table_display <- reactive({
      df <- sample_qc()
      df$reason <- apply(df[, c("flag_signal", "flag_detected", "flag_cor")], 1, function(r) {
        reasons <- c("low/high signal", "low/high detected features", "low cohort correlation")[r]
        if (length(reasons) == 0) "" else paste(reasons, collapse = "; ")
      })
      df$signal <- round(df$signal, 1); df$mean_cor <- round(df$mean_cor, 3)
      df[, c("sample", intersect("group", colnames(df)), "signal", "detected", "mean_cor", "reason")]
    })
```
`apply(df[, c("flag_signal","flag_detected","flag_cor")], 1, function(r) {...})` — iterates **row-wise** (`MARGIN = 1`) over the three logical flag columns; for each sample's row `r` (a 3-element logical vector), `c("low/high signal", ...)[r]` uses **logical indexing** — a classic, idiomatic R trick: indexing a character vector by a same-length logical vector returns exactly the elements where the logical is `TRUE`, in order — so a sample with `r = c(TRUE, FALSE, TRUE)` yields `c("low/high signal", "low cohort correlation")`, then `paste(..., collapse="; ")` joins them into one human-readable string; a sample with no flags (`r = c(FALSE,FALSE,FALSE)`) yields a zero-length character vector, and the explicit `if (length(reasons)==0) ""` guard prevents `paste(character(0), collapse="; ")` (which would actually already correctly return `""`, making this guard technically redundant — but harmless and arguably clearer to read).

```r
    output$qc_table <- DT::renderDataTable({
      df <- qc_table_display()
      flagged_only <- df[df$reason != "", , drop = FALSE]
      DT::datatable(if (nrow(flagged_only) > 0) flagged_only else df[0, ], rownames = FALSE, filter = "top",
                     options = list(pageLength = 8, scrollX = TRUE), class = "stripe hover compact")
    })
```
**A deliberate, important display choice:** the on-screen table shows **only flagged samples** (`df[df$reason != "", ]`) — not the full per-sample QC table, even though the *download* (below) exports everything. `if (nrow(flagged_only) > 0) flagged_only else df[0, ]` — if nothing was flagged, shows an explicitly **empty** table (`df[0, ]` — zero rows, but the correct columns/types preserved) rather than falling back to showing all samples, which would otherwise be visually confusing (a table suddenly full of "clean" rows right after a "0 samples flagged" success message would look like a display bug). `output$download_qc <- downloadHandler(..., content = function(file) write.csv(qc_table_display(), file, row.names = FALSE))` — correctly exports the **full** `qc_table_display()` (every sample, flagged or not, with its `reason` column, blank for clean ones) — a materially different, and more complete, dataset than what's shown on screen; worth stating explicitly in your Methods that the visible table and the downloadable CSV are not the same view, by design.

---

### Block M — QC → "Normalised data": the check (L369-462)

```r
    output$norm_color_by_ui <- renderUI({
      cols <- setdiff(colnames(qc_target()$meta), "sample")
      selectInput(ns("norm_color_by"), "Color by", choices = cols,
                  selected = if ("group" %in% cols) "group" else cols[1])
    })
```
Same "render the input itself, don't pre-declare + `update*()`" rationale as Block I's three pickers (the file's own comment at L369-376 restates the identical reasoning). `setdiff(colnames(meta), "sample")` — every metadata column except the sample-ID column itself (coloring a PCA plot by sample ID would be meaningless — one color per point). `selected = if ("group" %in% cols) "group" else cols[1]` — sensible default: prefer `group` if present (it almost always is, per earlier notes), otherwise just the first available column, using scalar `if`/`else` as an expression (same idiom seen throughout this codebase).

```r
    norm_check <- eventReactive(input$run_norm_btn, {
      expr <- qc_target()$expr
      diag <- summarize_norm_diagnostics(expr)
      qs <- apply(expr, 2, stats::quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
      box_df <- data.frame(sample = colnames(expr), ymin = qs[1, ], lower = qs[2, ], middle = qs[3, ],
                            upper = qs[4, ], ymax = qs[5, ])
      persample <- data.frame(sample = colnames(expr), median = round(qs[3, ], 2),
                                q25 = round(qs[2, ], 2), q75 = round(qs[4, ], 2))
      list(diag = diag, box_df = box_df, persample = persample, pca = pca_of(expr))
    })
```
Click-gated (`eventReactive` on `run_norm_btn`), so changing `norm_color_by`/`norm_pc_x`/`norm_pc_y`/the checkboxes does **not** re-trigger this expensive block (PCA on a full gene-by-sample matrix is not free) — only affects how the *already-computed* result (`norm_check()`) is subsequently *displayed*, in the render functions below. `summarize_norm_diagnostics(expr)` — `global.R:1550-1560`: computes, per-sample, the median and IQR, then the **standard deviation of those per-sample medians/IQRs across samples** (`median_sd`, `iqr_sd`) plus the matrix's overall max/min — the actual diagnostic signal is "how much do samples disagree with each other," not any single sample's own spread. `apply(expr, 2, stats::quantile, probs = c(.05,.25,.5,.75,.95), na.rm=TRUE)` — column-wise (per-sample, `MARGIN=2`) five-number-style quantile summary, used to build both a box-and-whisker-plot-ready data frame (`box_df`, with `ymin`/`lower`/`middle`/`upper`/`ymax` — exactly the five aesthetics `ggplot2::geom_boxplot(stat="identity")` expects when you supply pre-computed quantiles rather than raw values) and a simpler per-sample summary table (`persample`). `pca_of(expr)` — `global.R:1760-1769`: drops zero-variance genes first (`m[apply(m, 1, sd) > 0, , drop=FALSE]` — necessary because `prcomp(scale.=TRUE)` would otherwise divide by a zero standard deviation and error), then runs `prcomp(t(m), scale.=TRUE)` (samples as rows, as `prcomp` expects observations-in-rows), computing `var_exp` (percent variance per component, rounded to 1 decimal) and returning the first `n_pc` (default 5) components as a data frame plus a `sample` column.

```r
    output$norm_summary_ui <- renderUI({
      if (!isTruthy(input$run_norm_btn) || input$run_norm_btn == 0) {
        return(div(class = "empty-note", icon("circle-info"), "Not run yet - ..."))
      }
      d <- norm_check()$diag
      needs <- needs_quantile_norm(d)
      tagList(
        p(icon("circle-info"), " Spread of per-sample medians (SD): ", strong(sprintf("%.3f", d$median_sd)),
          "; spread of per-sample IQRs (SD): ", strong(sprintf("%.3f", d$iqr_sd)),
          "; max value in matrix: ", strong(format(round(d$max_value, 1), big.mark = ",")), "."),
        p(class = "empty-note", icon(if (!needs) "check" else "triangle-exclamation"),
          if (!needs) "Samples look well aligned and log-scaled - quantile normalisation likely isn't needed."
          else "This looks like it needs quantile normalisation: ... See Normalise this dataset below.")
      )
    })
```
`needs_quantile_norm(diag)` — `global.R:1562-1564`, a **one-line, three-condition** decision rule: `diag$max_value > 100 || diag$median_sd > 0.5 || diag$iqr_sd > 0.5` — i.e., "the data still looks linear-scale (a log2 microarray/RNA-seq value essentially never exceeds ~20-30 in real biological data, so a max above 100 is a strong linear-scale signal), **or** per-sample medians disagree by more than 0.5 on whatever scale the data is on, **or** per-sample IQRs disagree by more than 0.5." This is explicitly documented (`global.R:1481-1543`, the "THE NORMALISATION IMPLEMENTATION" audit block) as **the same rule** the underlying research pipeline's own preprocessing script (`scripts/00_shared/03_normalize_batch.R`) uses — i.e., this interactive check is not a separate, potentially-inconsistent reimplementation, it calls the identical function the batch pipeline itself calls. This is a strong, citable point for a Reproducibility section: the live, on-demand QC check and the pipeline's own automated decision are provably the same code path, not two independently-maintained rules that could silently drift apart.

`output$norm_views_ui` (L410-429) — renders three sub-panels: a per-sample box plot (`norm_dist_plot`), a scree plot (`norm_scree_plot`), and a PCA scatter + summary table side by side (`norm_pca_plot`/`norm_table`).

```r
    output$norm_dist_plot <- renderPlot({
      req(norm_check(), input$norm_color_by)
      box_df <- norm_check()$box_df
      color_col <- input$norm_color_by
      meta <- qc_target()$meta
      join_cols <- intersect(c("sample", color_col), colnames(meta))
      box_df <- merge(box_df, meta[, join_cols, drop = FALSE], by = "sample", all.x = TRUE)
      if (!color_col %in% colnames(box_df)) { box_df[[color_col]] <- "all samples"; }
      ggplot(box_df, aes(x = reorder(sample, middle), ymin = ymin, lower = lower, middle = middle,
                           upper = upper, ymax = ymax, fill = .data[[color_col]])) +
        geom_boxplot(stat = "identity", width = 0.7, linewidth = 0.15) +
        scale_fill_manual(values = arthomix_pair(box_df[[color_col]])) +
        ...
```
`join_cols <- intersect(c("sample", color_col), colnames(meta))` — defensively only requests columns that actually exist in `meta` (guards against `color_col` somehow not existing, though `norm_color_by_ui`'s own choices are built directly from `colnames(meta)`, so this should always succeed in practice — belt-and-suspenders again). `if (!color_col %in% colnames(box_df)) { box_df[[color_col]] <- "all samples"; }` — a genuine fallback: if the merge somehow didn't bring in the color column (e.g., `join_cols` ended up missing it), fill it with a constant placeholder string so the subsequent `fill = .data[[color_col]]` aesthetic mapping and `scale_fill_manual()` call don't error on a nonexistent column — a small but real defensive-programming detail. `.data[[color_col]]` — `ggplot2`'s **tidy-eval** pronoun: `.data` lets you reference a column by a *string held in a variable* (`color_col`) inside an `aes()` call, which ordinarily expects bare, unquoted column names — this is the standard modern `ggplot2` idiom for "the column to plot is chosen dynamically at runtime," used correctly throughout this file (also in `plot_pca_advanced()`, `global.R:1797`). `geom_boxplot(stat = "identity", ...)` — `stat="identity"` tells `geom_boxplot()` to use the `ymin`/`lower`/`middle`/`upper`/`ymax` values **directly as given**, rather than its default behavior of computing quantiles itself from raw per-group values — necessary here since the quantiles were already computed by hand (`norm_check()`'s `apply(..., stats::quantile, ...)` call) rather than being derivable from a single "value" column `ggplot2` could summarize on its own. `arthomix_pair(box_df[[color_col]])` — this app's fixed categorical palette (`global.R:1415-1421`): sorts the unique levels alphabetically and assigns them, in order, from a fixed 7-color hue sequence — meaning **the same category always gets the same color across every plot in the app** (e.g., "RA" is always the same shade wherever it appears), a genuinely good visual-consistency practice worth naming explicitly, with the caveat that a factor with more than 7 levels "runs out" of distinct colors (per that function's own comment) rather than cycling — worth testing if any reachable metadata column here could exceed 7 levels (e.g., a `dataset` column with many source studies, or a continuous-looking categorical column) — see §7, Issue O-8.

`output$norm_scree_plot <- renderPlot({ scree_plot(norm_check()$pca$var_exp) })` and `output$norm_pca_plot <- renderPlot({ ... plot_pca_advanced(norm_check()$pca, qc_target()$meta, input$norm_color_by, pc_x=..., pc_y=..., show_ellipse=..., show_labels=...) })` — both delegate to shared `global.R` helpers already explained in full above (`scree_plot()`, `plot_pca_advanced()`); nothing new to add for their internals here beyond noting that `plot_pca_advanced()`'s own ellipse logic (`ok_groups <- names(which(table(df[[color_by]]) >= 4))`) silently **skips** drawing a confidence ellipse for any group with fewer than 4 samples (a `stat_ellipse()` correctness requirement, not a bug — an ellipse fit to 1-3 points is not statistically meaningful) — worth stating explicitly if you show this plot in your thesis with a small group present, so a missing ellipse for one group isn't misread as a rendering failure.

---

### Block N — QC → "Normalised data": apply / adopt live normalization (L464-546)

```r
    norm_apply_result <- eventReactive(input$apply_norm_btn, {
      expr <- as.matrix(qc_target()$expr)
      normalized <- limma::normalizeBetweenArrays(expr, method = "quantile")
      list(
        diag_before = summarize_norm_diagnostics(expr), diag_after = summarize_norm_diagnostics(normalized),
        box_before = norm_check()$box_df,
        box_after = { ... },
        expr_after = normalized
      )
    })
```
`limma::normalizeBetweenArrays(expr, method = "quantile")` — the exact same function call `mod_preprocessing.R`'s Batch Correction tab uses for microarray/log-scale data (per `global.R`'s own audit-comment block, L1508-1516) — i.e., this "try it live" button is not a toy/approximation, it runs the real normalization method the pipeline itself would apply, on whatever dataset is currently under inspection. `box_after` recomputes the same five-number-summary-per-sample logic inline (duplicating the pattern from `norm_check()` rather than factoring it into a shared helper — a minor code-duplication/maintainability note, not a correctness issue). Note this whole block silently **depends on `norm_check()` already having been run** (`box_before = norm_check()$box_df`) — if a user somehow reached `apply_norm_btn` without first clicking `run_norm_btn`, this would error inside the `eventReactive`, since `norm_check()` (also an `eventReactive`) has no cached value yet to return. In practice this can't happen because `norm_apply_ui` (which contains `apply_norm_btn`) is itself only rendered via `req(norm_check())` (Block below) — i.e., the "Apply" button doesn't even exist in the DOM until "Run normalisation check" has completed at least once — a correct, UI-enforced sequencing rather than a server-side guard, worth naming explicitly as the actual mechanism preventing this failure mode.

```r
    output$norm_apply_ui <- renderUI({
      req(norm_check())
      needs <- needs_quantile_norm(norm_check()$diag)
      tagList(
        box(
          width = 12, title = "Normalise this dataset", status = if (needs) "warning" else "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Runs the same quantile normalisation (limma::normalizeBetweenArrays) Preprocessing applies..."),
          actionButton(ns("apply_norm_btn"), "Apply quantile normalisation", ...),
          uiOutput(ns("norm_apply_result_ui"))
        )
      )
    })
```
`status = if (needs) "warning" else "primary"` — the box's header/border color changes to amber ("warning") specifically when `needs_quantile_norm()` says this dataset likely needs it — a nice, correct use of the same decision function to drive a *visual* cue consistently with the *textual* verdict shown just above in `norm_summary_ui` (both ultimately trace back to the identical `needs_quantile_norm(norm_check()$diag)` call, so they can never disagree with each other — a good internal-consistency property, worth verifying/citing as such).

```r
    output$norm_apply_result_ui <- renderUI({
      req(norm_apply_result())
      r <- norm_apply_result()
      tagList(
        br(),
        fluidRow(
          valueBox(sprintf("%.3f -> %.3f", r$diag_before$median_sd, r$diag_after$median_sd), "Spread of medians (SD)",
                    icon = icon("chart-line"), color = if (r$diag_after$median_sd < 0.5) "green" else "red", width = 6),
          valueBox(sprintf("%.3f -> %.3f", r$diag_before$iqr_sd, r$diag_after$iqr_sd), "Spread of IQRs (SD)",
                    icon = icon("chart-line"), color = if (r$diag_after$iqr_sd < 0.5) "green" else "red", width = 6)
        ),
        ...
        if (identical(input$qc_source, "uploaded")) {
          div(
            actionButton(ns("adopt_norm_btn"), "Use this normalised version for every sub-module", icon = icon("check"), class = "btn-success btn-sm"),
            uiOutput(ns("adopt_norm_msg"))
          )
        } else {
          p(class = "empty-note", icon("circle-info"),
            "This is one of the app's fixed reference datasets, so it stays read-only here - ... Upload your own data on the Dataset tab to normalise it and use the result app-wide.")
        }
      )
    })
```
Two before→after `valueBox`es, each independently colored green/red based on whether that specific metric now falls under the 0.5 threshold post-normalization — note these are **not** literally re-calling `needs_quantile_norm()` (which combines all three of its own conditions with `||`), but manually re-implementing the "< 0.5" half of that same logic per-metric — a minor, low-risk duplication (the two thresholds are numerically identical to the ones inside `needs_quantile_norm()`, so no actual disagreement is currently possible, but a future edit to one wouldn't automatically update the other — worth flagging as a maintainability note, §7 Issue O-4 extended).

**The activation gate:** `if (identical(input$qc_source, "uploaded")) { ... } else { ... }` — this single `if` is the entire mechanism deciding whether "adopt" is even offered. It checks the **QC tab's own picker value directly** (`input$qc_source`), not `norm_apply_result()` or any stored flag — meaning the check is re-evaluated fresh every time this `renderUI` runs, always reflecting whatever's *currently* selected in the picker, which (since `norm_apply_result()` itself was computed from `qc_target()`, itself derived from the *same* `input$qc_source` at the time "Apply" was clicked) will normally agree — but see §7, Issue O-2 for the scenario where a user changes the picker *after* running Apply but *before* this block re-renders, and what that implies.

```r
    observeEvent(input$adopt_norm_btn, {
      req(norm_apply_result(), identical(input$qc_source, "uploaded"))
      dataset$expr <- norm_apply_result()$expr_after
      dataset$source <- paste0(dataset$source, " (quantile-normalised)")
      output$adopt_norm_msg <- renderUI(
        div(class = "empty-note", icon("check"), "Done - every sub-module now reads the quantile-normalised version of your data.")
      )
    })
```
**The one and only place in this entire module that mutates shared, cross-module state.** `req(norm_apply_result(), identical(input$qc_source, "uploaded"))` — a second, server-side re-check of the exact same "uploaded" condition already gating the button's visibility (defense in depth, consistent with the pattern seen throughout `mod_dataset.R`). `dataset$expr <- norm_apply_result()$expr_after` — directly overwrites the shared, **active** expression matrix (not `staged_expr`) — this is precisely the "second, narrower activation path" flagged in §1: unlike every other write to `dataset$expr` in this app (which only ever happens via Preprocessing's own merge/batch-correct/"Use this as the active dataset" button, per the Dataset-module notes' architecture diagram), this one button, on this one tab, bypasses Preprocessing entirely for the single case of "you uploaded your own already-active data and want to normalize it in place." `dataset$source <- paste0(dataset$source, " (quantile-normalised)")` — appends a provenance suffix to the existing source string (rather than replacing it) — so the string keeps growing if a user somehow did this more than once (not tested here, but a plausible edge case — see §7).

---

### Block O — Expression data tab server (L552-569)

```r
    expr_table_data <- reactive({
      m <- qc_target_expr()$expr
      data.frame(feature = rownames(m), round(m, 3), check.names = FALSE, stringsAsFactors = FALSE)
    })

    output$expr_table <- DT::renderDataTable({
      df <- expr_table_data()
      DT::datatable(
        df, rownames = FALSE, extensions = "FixedColumns",
        options = list(pageLength = 15, dom = "lfrtip", scrollX = TRUE, fixedColumns = list(leftColumns = 1)),
        class = "stripe hover compact"
      )
    }, server = TRUE)

    output$download_expr <- downloadHandler(
      filename = function() "expression_matrix.csv",
      content = function(file) data.table::fwrite(expr_table_data(), file)
    )
```
`data.frame(feature = rownames(m), round(m, 3), check.names = FALSE, ...)` — converts the numeric matrix into a data frame with an explicit `feature` column (the gene/probe ID, otherwise only implicit as row names, which `DT` cannot display/search directly) prepended to every (rounded-to-3-decimals, for display readability) sample column; `check.names = FALSE` prevents R's default behavior of "sanitizing" column names (e.g., replacing a sample name like `"GSM-123"`'s hyphen with a dot) — important here since sample IDs are meaningful identifiers that must survive unmodified.

`DT::renderDataTable({...}, server = TRUE)` — as corrected just above (§4, Block J), `server = TRUE` is `DT::renderDataTable()`'s own default, so writing it explicitly here changes nothing behaviorally versus the other five `DT::renderDataTable()` calls in this file that omit it — **every one of them is already server-side**. What this explicit argument *does* signal is authorial intent: the file's own comment (L547-550, physically preceding this block) explains why server-side processing specifically matters here — "these matrices run to tens of thousands of rows, and rendering that client-side would ship the whole thing to the browser at once" — and even parenthetically notes "(the DT default)," confirming the code's own author was aware `server=TRUE` didn't need to be stated for it to apply. Writing it anyway, only on the one output where getting it wrong would actually matter at this app's data scale, is a reasonable defensive-documentation choice — worth citing as such, rather than (as an earlier draft of this document mistakenly did) as evidence of an actual client/server behavioral split across the module's tables. With `server = TRUE` (in force everywhere), `DT` requests only the currently-visible page of rows from the R server on each interaction (sort/filter/page), which is what makes a genome-wide feature-by-sample table tractable in the browser at all. `extensions = "FixedColumns"`, `fixedColumns = list(leftColumns = 1)` — a `DataTables` extension that pins the first column (`feature`) in place while horizontally scrolling through potentially dozens/hundreds of sample columns — without it, scrolling right would scroll the gene-ID column out of view too, making every visible value unattributable to a gene. `dom = "lfrtip"` — the default full `DataTables` layout string (length-selector, filter box, table, info, pagination — each letter is one control) spelled out explicitly here (as opposed to the trimmed `"tp"` used for the small Missing-values-by-field table) since a full-size, searchable expression matrix benefits from every control.

**File-organization note (already flagged in §2.3/§3):** despite living under a server-code comment that reads "QC tab: browse the expression matrix itself" (L547), this block's outputs (`expr_table`, `download_expr`) are consumed by the **"Expression data"** `tabPanel` (UI L64-69), not by any tab inside `qc_tabs`. This is a harmless but genuine internal-documentation inconsistency — see §7, Issue O-1.

---

### Block P — QC → "Group" (filtering) (L573-711)

```r
    filter_spec <- reactive({
      meta <- qc_target()$meta
      cols <- setdiff(colnames(meta), "sample")
      specs <- list()
      for (cl in cols) {
        x <- meta[[cl]]
        if (is.numeric(x)) {
          rng <- suppressWarnings(range(x, na.rm = TRUE))
          if (is.finite(rng[1]) && is.finite(rng[2]) && rng[1] < rng[2]) {
            specs[[cl]] <- list(type = "numeric", min = floor(rng[1]), max = ceiling(rng[2]))
          }
        } else {
          u <- sort(unique(na.omit(as.character(x))))
          if (length(u) >= 2 && length(u) <= 30) {
            specs[[cl]] <- list(type = "categorical", choices = u)
          }
        }
      }
      specs
    })
```
This is the module's most **generic, data-driven** piece of logic — it inspects the metadata's own columns at runtime and decides what kind of filter control (if any) makes sense for each, rather than hardcoding `group`/`sex` as the only filterable fields (the file's own comment, L571-577, states this explicitly). For each column: if numeric, `range(x, na.rm=TRUE)` (wrapped in `suppressWarnings()` since `range()` on an all-`NA` column emits a warning before returning `c(Inf, -Inf)`) — only becomes a filterable numeric spec if both bounds are finite **and** genuinely different (`rng[1] < rng[2]`) — a column that's numeric but constant (or entirely `NA`) is correctly excluded, since a slider with `min==max` (or undefined bounds) isn't a meaningful filter. `floor()`/`ceiling()` round the slider's bounds outward to whole numbers (a minor UX simplification — avoids a slider labeled e.g. "23.7 to 61.2"). If not numeric: `sort(unique(na.omit(as.character(x))))` — the distinct non-missing values as strings, sorted; becomes a categorical spec **only if between 2 and 30 distinct values** — `length(u) < 2` (a constant column) is excluded as not filterable, and `length(u) > 30` (e.g., a `sample` ID-like column, or a free-text field) is excluded per the file's own comment: "picking one of hundreds of individual sample IDs isn't a useful filter." Note `sample` itself is already excluded up front via `setdiff(colnames(meta), "sample")`, so the `> 30` cap is really a second line of defense against any *other* high-cardinality column (e.g., a numeric-looking ID that happened to be stored as character).

```r
    output$filters <- renderUI({
      specs <- filter_spec()
      validate(need(length(specs) > 0, "No filterable columns in the current metadata."))
      tagList(
        lapply(names(specs), function(cl) {
          s <- specs[[cl]]
          fid <- paste0("f_", cl)
          if (s$type == "categorical") {
            pickerInput(ns(fid), cl, choices = s$choices, selected = character(0),
                        multiple = TRUE, options = list(`actions-box` = TRUE, title = paste0("Choose ", cl, "... (optional)")))
          } else {
            sliderInput(ns(fid), cl, min = s$min, max = s$max, value = c(s$min, s$max))
          }
        }),
        fluidRow(
          column(6, actionButton(ns("apply_btn"), "Apply filters", ...)),
          column(6, actionButton(ns("reset_btn"), "Reset", ...))
        )
      )
    })
```
`fid <- paste0("f_", cl)` — dynamically constructs a **new** input ID per metadata column (e.g., `"f_group"`, `"f_age"`) — this is genuinely dynamic UI generation: the *set* of inputs that exist in this session is determined entirely by whatever columns the currently-inspected dataset's metadata happens to have, not a fixed list known at `mod_overview_ui()` definition time. `pickerInput(ns(fid), cl, choices = s$choices, selected = character(0), multiple = TRUE, options = list(`actions-box`=TRUE, ...))` — `shinyWidgets::pickerInput`, a richer multi-select widget than plain `selectInput(multiple=TRUE)`; `` `actions-box` = TRUE `` adds "select all / deselect all" convenience buttons; `selected = character(0)` means **no filter is pre-applied** — every categorical filter starts in an explicitly "show everything" state, consistent with the tab's own descriptive text ("Everything shows by default"). `sliderInput(ns(fid), cl, min=s$min, max=s$max, value=c(s$min, s$max))` — a **range** slider (two-handled, since `value` is a length-2 vector) defaulting to the *full* observed range — same "nothing filtered yet" default philosophy as the categorical case.

```r
    observeEvent(input$reset_btn, {
      specs <- filter_spec()
      for (cl in names(specs)) {
        s <- specs[[cl]]; fid <- paste0("f_", cl)
        if (s$type == "categorical") updatePickerInput(session, fid, selected = character(0))
        else updateSliderInput(session, fid, value = c(s$min, s$max))
      }
    })
```
`updatePickerInput`/`updateSliderInput` — this is the "update an existing, already-rendered input in place" pattern the file's own comments elsewhere (Block I) explicitly say is unsafe to use for the *initial* render of a dynamically-inserted module's controls — but it's **correct and safe here**, because by the time a user can click `reset_btn`, the filter inputs (`f_*`) must already exist and be bound in the DOM (the button that triggers this reset lives in the very same `renderUI` block that created those inputs, so they're guaranteed to already be live in the session) — a good contrast case to understand exactly when `update*()` functions are and aren't reliable in this app's specific module-loading architecture (eager server instantiation + deferred UI insertion via `insertTab`).

```r
    filtered_meta <- eventReactive(input$apply_btn, {
      meta <- qc_target()$meta
      specs <- filter_spec()
      for (cl in names(specs)) {
        s <- specs[[cl]]; val <- input[[paste0("f_", cl)]]
        if (s$type == "categorical") {
          if (length(val) > 0) meta <- meta[is.na(meta[[cl]]) | meta[[cl]] %in% val, , drop = FALSE]
        } else if (!is.null(val)) {
          meta <- meta[is.na(meta[[cl]]) | (meta[[cl]] >= val[1] & meta[[cl]] <= val[2]), , drop = FALSE]
        }
      }
      validate(need(nrow(meta) > 0, "No samples match the current filter combination - try Reset."))
      meta
    })
```
`input[[paste0("f_", cl)]]` — **dynamic input access by constructed name**: since the filter inputs' IDs were themselves built at runtime (`fid <- paste0("f_", cl)` above), reading them back requires the same string-building, via `input[[...]]` (double-bracket, computed-name access) rather than `input$literal_name` (which only works for a name known at write-time). **Every filter clause explicitly allows `NA` through** — `is.na(meta[[cl]]) | meta[[cl]] %in% val` (categorical) and `is.na(meta[[cl]]) | (meta[[cl]] >= val[1] & meta[[cl]] <= val[2])` (numeric) — meaning a sample with a **missing** value for a filtered column is never excluded by that filter, only samples with a present-but-non-matching value are. This is a genuine, non-obvious design decision worth stating explicitly and considering carefully in your Validation chapter: it means "filter to Sex = Female" will still include samples where sex is simply unrecorded, which may or may not be the behavior a user expects from a filter UI — see §7, Issue O-10 for the concrete testable consequence.

For categorical: `if (length(val) > 0)` — an **empty selection is treated as "no filter"** (matches the `pickerInput`'s "optional" framing — leaving a filter untouched, i.e. `character(0)` selected, does not narrow anything), whereas a non-empty selection restricts to exactly those levels (plus `NA`s, per above). For numeric: `!is.null(val)` — a slider's value is never actually `NULL` once rendered (it always has a two-element numeric vector), so this guard is more about defensive robustness than a reachable "unset" state in practice. `validate(need(nrow(meta) > 0, "No samples match the current filter combination - try Reset."))` — a genuinely useful, specific message pointing the user at the exact remedy (Reset), rather than a generic "no data" error.

`output$filter_summary_ui`, `output$filtered_views_ui`, `output$group_plot`, `output$dataset_plot`, `output$meta_table`, `output$download_meta` (L645-711) — render value boxes, a `group`×`sex`-faceted composition bar chart (`df %>% count(group, sex)`, `dplyr`'s pipe-and-count idiom, then `ggplot(..., fill=group) + facet_wrap(~sex)` — only facets if at least one non-missing `sex` value exists), an optional per-`dataset`-column bar chart (only rendered `if ("dataset" %in% colnames(filtered_meta()))` — correctly conditional, since not every reachable metadata source necessarily has a `dataset` column, e.g. a single individual GEO source's own harmonized metadata may or may not include one depending on `eset_harmonize_meta()`'s output — worth double-checking against `global.R:731-743`, which does in fact always add a `dataset` column for the eset-derived path, so in practice this condition is close to always true for this app's reachable inputs, but the code doesn't assume it), and a final filtered, still independently sortable/filterable `DT` table (server-side, per the correction in §4 Block J) plus its own CSV download — all structurally identical in construction to patterns already fully explained above (`DT::datatable(..., filter="top")`, `arthomix_pair()`-colored `ggplot2` bars, `downloadHandler`).

---

## 5. DATA FLOW (per sub-tab)

### Datasets tab
```
GEO_SOURCES (global.R, static list of 4) --> lapply --> per-source get_raw_eset(gse)
  (reads/caches data/preloaded/transcriptomics/raw/<gse>_raw.rds, or NULL if absent)
  --> one info-card per source: title/platform/dims (from the eset) + role/used_in
      (from GEO_SOURCES itself) + a live NCBI GEO hyperlink (geo_link())
  --> output$sources_ui (renderUI) --> module-grid of cards
```

### Metadata / Expression data / QC tabs (shared resolver)
```
User: picks a value in this tab's own "Dataset to inspect" selectInput
  --> input$qc_source_{meta|expr|<none>} (string: a GSE id, or "uploaded")
  --> qc_target_{meta|expr|<none>}() [reactive] --> resolve_qc_source(source_id)
        "uploaded"  --> dataset$expr / dataset$meta directly (the ACTIVE dataset, never staged_*)
        a GSE id    --> load_individual_dataset(gse_id) [global.R] --> raw ExpressionSet's
                        exprs()/harmonized pData(), a bespoke counts-file read (GSE89408), or a
                        merged_training_subset() fallback (when the raw file is absent)
  --> {expr: matrix, meta: data.frame, label: string} --> every tab-specific output below
```

### Metadata tab, specifically
```
qc_target_meta() --> output$understand_ui (value boxes: nrow(meta), unique(group), nrow(expr), unique(sex))
                  --> output$meta_table_full (server-side DT, DT's own default - see corrected §4 note; full unfiltered meta)
                  --> output$download_meta_full (CSV)
```

### Expression data tab, specifically
```
qc_target_expr()$expr --> expr_table_data() [rounds to 3dp, adds a `feature` column from rownames]
                        --> output$expr_table (server-side DT, FixedColumns)
                        --> output$download_expr (CSV, via data.table::fwrite)
```

### QC -> Missing values
```
qc_target()$meta --> missing_audit() [reactive, always-on] --> output$missing_plot (ggplot bar)
                                                              --> output$missing_table (DT)
```

### QC -> Outliers
```
User clicks "Run outlier detection" (input$run_qc_btn) with input$mad_k set
  --> sample_qc() [eventReactive] --> compute_sample_qc(expr, mad_k) [global.R] --> merge with meta$group
  --> output$qc_summary_ui (flagged-sample count valueBox)
  --> output${signal,detected,cor}_plot (qc_bar_plot(), global.R)
  --> qc_table_display() --> output$qc_table (flagged samples only, DT)
                          --> output$download_qc (CSV, all samples)
```

### QC -> Normalised data
```
User clicks "Run normalisation check" (input$run_norm_btn)
  --> norm_check() [eventReactive] --> summarize_norm_diagnostics(), per-sample quantiles, pca_of() [global.R]
  --> output$norm_summary_ui (needs_quantile_norm() verdict text)
  --> output$norm_dist_plot / norm_scree_plot / norm_pca_plot / norm_table
User clicks "Apply quantile normalisation" (input$apply_norm_btn)
  --> norm_apply_result() [eventReactive] --> limma::normalizeBetweenArrays(expr, "quantile")
  --> output$norm_apply_result_ui (before/after valueBoxes + plots)
  --> IF input$qc_source == "uploaded": "Use this normalised version..." button appears
       --> observeEvent(adopt_norm_btn) --> dataset$expr <- normalized; dataset$source <- ...(annotated)
           [THE ONLY WRITE TO SHARED STATE IN THIS ENTIRE FILE]
```

### QC -> Group
```
qc_target()$meta --> filter_spec() [reactive, inspects every column's type/cardinality]
                  --> output$filters (dynamic per-column pickerInput/sliderInput, IDs "f_<col>")
User sets some filters, clicks "Apply filters" (input$apply_btn)
  --> filtered_meta() [eventReactive] --> per-column filter loop (NA always passes) --> validate(>=1 row)
  --> output$filter_summary_ui (value boxes) / group_plot / dataset_plot / meta_table / download_meta
User clicks "Reset" (input$reset_btn) --> observeEvent --> update*Input() every f_<col> back to "show all"
```

---

## 6. RESULTS

| Output | What it shows | Source computation | Controlling input(s) | Updates when | If input missing/invalid | If dataset is degenerate |
|---|---|---|---|---|---|---|
| `sources_ui` | 4 provenance cards | `get_raw_eset()` per `GEO_SOURCES` entry | none | Effectively once per session | N/A | Missing raw file → explicit "Raw file not found on disk" card, not a crash |
| `understand_ui` | 4 KPI boxes | `nrow()`/`unique()`/`na.omit()` on `qc_target_meta()` | `qc_source_meta` | Picker changes | Blocked by `req()` inside `qc_target_meta()` | An all-`NA` `group`/`sex` column would report `0` groups/categories, not an error |
| `meta_table_full` | Full metadata table | `qc_target_meta()$meta`, unmodified | `qc_source_meta` | Picker changes | — | Empty metadata (0 rows) would render an empty table, no explicit guard |
| `expr_table` | Full expression matrix | `qc_target_expr()$expr`, rounded | `qc_source_expr` | Picker changes | — | 0-row matrix renders an empty table |
| `missing_plot`/`missing_table` | % missing per metadata field | `missing_audit()` | `qc_source` | Picker changes | — | Always computable, even on a tiny/degenerate dataset |
| `qc_summary_ui`/`qc_plots_ui`/`qc_table` | Flagged-sample count, 3 diagnostic plots, flagged-only table | `sample_qc()` via `compute_sample_qc()` | `run_qc_btn` (click), `mad_k` | Button clicked | "Not run yet" placeholder before first click | A dataset with `mad(x)==0` for a metric → 0 flags for that metric (guarded, not an error) |
| `norm_summary_ui`/`norm_views_ui` | Normalization-need verdict + 4 diagnostic views | `norm_check()` via `summarize_norm_diagnostics()`/`pca_of()` | `run_norm_btn` (click); display-only: `norm_color_by`, `norm_pc_x/y`, ellipse/label checkboxes | Button clicked (compute); display inputs change (re-render only) | "Not run yet" placeholder | A group with <4 samples silently gets no ellipse (`plot_pca_advanced()`'s own guard) |
| `norm_apply_result_ui` | Before/after normalization comparison | `norm_apply_result()` via `limma::normalizeBetweenArrays()` | `apply_norm_btn` (click) | Button clicked | UI only exists once `norm_check()` has run | — |
| `adopt_norm_msg` | "Done" confirmation | `observeEvent(adopt_norm_btn)` | `adopt_norm_btn` (click), gated on `qc_source == "uploaded"` | Button clicked | Button itself isn't rendered unless the gate is met | — |
| `filter_summary_ui`/`filtered_views_ui`/`meta_table`(filtered)/`group_plot`/`dataset_plot` | Filtered cohort composition + table | `filtered_meta()` | `apply_btn` (click); filter inputs `f_*` | Button clicked | "Set any filters..." placeholder before first click | `validate(nrow>0)` blocks an empty result with an actionable message |

### User perspective
Every sub-tab follows one of two patterns: **always-on** (Datasets, Metadata, Expression data, Missing values — computed the moment the tab/picker is set) or **click-to-run** (Outliers, Normalised data, Group filtering — nothing happens, and a neutral placeholder is shown, until an explicit button click). This distinction tracks computational cost directly: the always-on checks are cheap (table subsetting, simple percentages); the click-gated ones involve real numerical work (robust statistics over the full matrix, PCA, `limma` normalization) that shouldn't silently re-run on every minor input tweak.

### Scientific/computational perspective
Every diagnostic shown here (missing-value audit, MAD-based outlier flags, quantile-normalization need) is computed by the **same functions** the underlying pipeline's own automated preprocessing script calls — this module is best understood as an interactive, on-demand re-exposure of that pipeline's own internal QC logic, not a separately-invented "simplified" version of it. This is a genuinely strong point for your Methods/Reproducibility chapter: it directly means the numbers a user sees in this tab match, code-for-code, what actually decided the pipeline's own preprocessing choices — see the "THE NORMALISATION IMPLEMENTATION" audit block in `global.R:1481-1543`, which documents this equivalence explicitly and was evidently written to make exactly this claim auditable.

---

## 7. CODE VALIDATION

### Correctness — generally sound, with the same "silent-failure-mode" character as `mod_dataset.R`

- **Robust statistics choice (median/MAD)** for outlier detection — correct and well-justified; the one-sided `flag_cor` and the `mad(x)==0` degenerate-case guard are both genuinely careful details.
- **The normalization-need heuristic is provably shared**, not duplicated, between this tab and the pipeline's own batch-correction logic (`needs_quantile_norm()`/`summarize_norm_diagnostics()`, both called from `global.R`, both used identically here and in `mod_preprocessing.R`) — a real reproducibility strength worth citing directly.
- **Sample-alignment via `merge(..., all.x=TRUE)`** in the Outliers block is a correct left join, robust to a metadata/expression sample mismatch (keeps every expression-matrix sample, doesn't silently drop unmatched ones).
- **Every `DT` table in this module is server-side processed**, since `DT::renderDataTable()`'s own default is `server = TRUE` — `expr_table` (the one call that states it explicitly) is not actually a different processing mode from the other five tables in this file, just the one place the author chose to document the default in code, matched to where getting it wrong would matter most (the potentially-huge expression matrix). An earlier draft of this document mistakenly described this as a deliberate client/server split; it isn't — see the correction in §4, Block J.

### Issues found (documented, per your instruction, as observations — **not applied as fixes**)

**Issue O-1 — Naming/labeling drift between UI tab titles, server-code comments, and file structure.**
- Evidence: (a) the sidebar's "Datasets" nav item (`ui.R:1335`) and this module's first sub-tab, also titled "Datasets" (`mod_overview.R:33`), are easily confused with the separate, always-visible "Dataset" tab (`mod_dataset_ui`, `ui.R:1360`) — three near-identical names for two different files/purposes. (b) the server code implementing the **"Expression data"** tab's table is introduced by a comment reading **"QC tab: browse the expression matrix itself"** (`mod_overview.R:547-550`) — the comment names the wrong tab. (c) the "Expression data" tab's own descriptive text ("The actual expression matrix behind every check below") reads as if it's describing the QC tab's checks, when it's physically on a separate tab from all of them.
- Consequence: none for correctness (the code runs exactly as intended regardless of what its comments say) — but a real risk for **you**, writing a thesis from this code: if you describe "the QC tab's expression matrix browser" you'd be describing a UI element that doesn't exist under that tab. Verify tab attribution directly against the UI function (as this document does), not against server-code comments alone.
- Fix: not applied (comment-only issue, informational).

**Issue O-2 — The live "adopt normalised data" path is a second, narrower activation mechanism outside Preprocessing's documented single-gate design.**
- Evidence: `observeEvent(input$adopt_norm_btn, ...)` (L538-545) writes `dataset$expr`/`dataset$source` directly; every other write to those two fields in the app happens only via `mod_preprocessing.R`'s "Use this as the active dataset" button (per the companion Dataset-module notes' architecture diagram, §12).
- Consequence for your thesis: state this explicitly as an intentional **second activation path**, scoped narrowly (only reachable when the inspected dataset is already the active, user-uploaded one) — not an architectural inconsistency to hide, but a design fact your Software Architecture chapter should name and justify (the file's own comment, L464-470, gives the justification: the four fixed GEO sources are a read-only reference and already have Preprocessing as their real pipeline, so only the user's own upload — which has no other normalization path yet — gets this shortcut).
- A related, untested edge case: if a user clicks "Apply" and then changes the QC-tab picker to a *different* dataset **before** the "adopt" button's own `renderUI` (`norm_apply_result_ui`) has re-rendered, could a stale `norm_apply_result()` (computed from the *previous* selection) be adopted while the visible dataset info suggests a different one? Trace: `norm_apply_result()` is an `eventReactive` — its cached value only updates on a fresh `apply_norm_btn` click, so switching the picker without re-clicking "Apply" would leave `norm_apply_result()` referring to the *old* dataset while `input$qc_source` now reflects the *new* one; the `renderUI` at L499-526 re-evaluates on either change and would still show the old comparison (labeled generically, with no dataset-identifying text) with the *new* picker's "uploaded" gate check — if the new selection also happens to be "uploaded" (only possible if there's ever more than one thing recognized as "uploaded," which isn't currently possible per Issue O-5's finding that only one upload state exists at a time), this could adopt stale data. **Currently not exploitable** given only one dataset can ever be `"uploaded"` at a time, but worth stating as a latent coupling to watch if the "uploaded" concept is ever extended.
- Fix: not applied.

**Issue O-3 — The three independent "Dataset to inspect" pickers (Metadata/Expression data/QC) never share a selection.**
- Evidence: `selected = GEO_SOURCES[[1]]$gse` hardcoded identically in all three `renderUI` blocks (L211-212, L223-224, L235-236); no shared reactive value links them.
- Consequence: a user who picks "GSE15573" on Metadata, then switches to Expression data, sees GSE93272 again by default — an easy source of "why did the data change" confusion during real use or a live thesis-defense demo. This is an explicit, commented trade-off (§2.0), not an oversight, but worth rehearsing/anticipating if you demo this module live.
- Fix: not applied.

**Issue O-4 — Several string/threshold constants are duplicated across files or within this file with no single source of truth.**
- The `"^Uploaded dataset:"` prefix check (`qc_source_choices()`, L189) must stay in lockstep with the exact string `mod_dataset.R:286` writes; no shared constant links them.
- The `< 0.5` threshold appears three times with the same numeric meaning: once inside `needs_quantile_norm()` (twice, for `median_sd`/`iqr_sd`), and again, independently, in the two `norm_apply_result_ui` `valueBox` `color=` expressions (L506, L508) — a future change to the pipeline's threshold would need to be mirrored by hand in this display logic to stay consistent.
- Fix: not applied; flagged as a maintainability point for a software-engineering-quality discussion, not a present correctness bug (the values agree today).

**Issue O-5 — `qc_source_choices()` only recognizes an *uploaded* active dataset as offering a 5th "uploaded" option — a preloaded-individual-dataset or live-GEO-fetched active dataset is not recognized, even though it is genuinely different from all four fixed `GEO_SOURCES` entries.**
- Evidence: `grepl("^Uploaded dataset:", dataset$source %||% "")` (L189) — the only pattern checked; `mod_dataset.R`'s other two active-dataset-producing paths set `dataset$source` to strings starting `"Individual dataset: "` / `"NCBI GEO: "` / (the default) `"Example dataset: "`, none of which match.
- Consequence: if a user's actual active dataset (via Preprocessing) is, say, a live GEO fetch of some new series, this module's Metadata/Expression data/QC tabs offer no way to inspect *that exact* dataset — only the four fixed references. This is a real, testable functional gap, not merely a naming quirk.
- Test: activate a dataset via a live GEO fetch (Dataset tab → Fetch → Load → Preprocessing → "Use this as the active dataset"), then open Overview → Metadata; confirm the picker's choices are only the four fixed GSE ids (no 5th option), even though `dataset$source` clearly differs from all four.
- Fix: not applied.

**Issue O-6 — No cross-check that `meta`'s row count and `expr`'s column count agree, in any of the value-box/summary displays.**
- Evidence: `understand_ui` reports `nrow(meta)` (Samples) and, separately, `nrow(t$expr)` (Features) — but never checks/asserts `ncol(t$expr) == nrow(meta)`.
- Consequence: for a hypothetically malformed `resolve_qc_source()` result where the two disagree (not currently reachable through any of this app's own loaders, all of which already align expr/meta upstream — see the Dataset-module notes' correctness discussion of `intersect()`/`match()`), this display would silently report two inconsistent sample counts side by side rather than flagging the mismatch.
- Fix: not applied — noting this as a defensive gap worth testing against a deliberately malformed input if you want to characterize worst-case behavior for your thesis's Validation section, even though every *actual* code path in this app currently prevents the mismatch from arising.

**Issue O-7 — The Metadata tab's CSV download filename fallback (`%||% "metadata"`) does not fire on an `NA` `dataset` column value, only on a genuinely `NULL`/absent column.**
- Evidence: `qc_target_meta()$meta$dataset[1] %||% "metadata"` (L270) — `%||%` (`global.R:817`) only substitutes when its left side is `NULL`; if the `dataset` column exists but its first value happens to be `NA` (a data quality issue itself worth flagging, but distinct), `paste0(NA, "_metadata.csv")` produces the literal filename `"NA_metadata.csv"`, not the friendlier fallback.
- Fix: not applied.

**Issue O-8 — `arthomix_pair()`'s fixed 7-color palette silently runs out of distinct colors for any metadata column exceeding 7 unique levels used as a plot fill/color aesthetic (e.g., `norm_dist_plot`'s `norm_color_by`, or `group_plot`'s `group`).**
- Evidence: `global.R:1415-1421` — `pal[seq_along(levels)]` where `pal` has exactly 7 entries; a level beyond the 7th receives `NA` from `pal[8]`, which `ggplot2::scale_fill_manual()` would render as a missing/grey fill (or, depending on `ggplot2` version/settings, drop that level's data with a warning) rather than a distinguishable color.
- Consequence: if a user's metadata (uploaded data, or a `dataset`/`batch` column with many source studies) has more than 7 categories and is chosen as the "Color by" field, the resulting plot could be visually misleading (indistinguishable or missing colors for some categories) without any on-screen warning.
- Test: upload metadata with an 8+-level categorical column, select it as `norm_color_by`, and inspect `norm_dist_plot`'s legend/fill rendering.
- Fix: not applied.

**Issue O-9 — The `results` parameter is accepted by `mod_overview_server()` but never read or written anywhere in this file.**
- Evidence: `function(id, dataset, results = NULL)` (L141) — no further occurrence of the identifier `results` anywhere in the ~570 remaining lines of server code.
- Consequence: per the chapter outline's own stated architectural convention (populating a shared `results` object "right after its main reactive fires," so the ArthOChat assistant can see a module's findings), **this module currently contributes nothing to that shared context** — the assistant cannot answer a question like "were any samples flagged as outliers in the currently inspected dataset?" grounded in this tab's actual computed state, even after a user has run the check, because nothing here ever writes to `results`. This is a legitimate, code-verified **gap** (not a bug) worth naming explicitly in your Discussion/Limitations section, consistent with the "self-aware disclosure" pattern your project's `AUDIT.md`/thesis outline already models for other modules.
- Fix: not applied (would require deciding what, if anything, this read-mostly module should expose — a design decision, not a one-line patch).

**Issue O-10 — Every "Group" tab filter treats a missing (`NA`) value as automatically passing the filter, for both categorical and numeric columns.**
- Evidence: `is.na(meta[[cl]]) | meta[[cl]] %in% val` / `is.na(meta[[cl]]) | (meta[[cl]] >= val[1] & meta[[cl]] <= val[2])` (L636, L638).
- Consequence: filtering to, say, `sex = "F"` will still include samples whose `sex` is unrecorded — a specific, testable, and possibly surprising behavior for anyone reading a filtered sample count and assuming it means "exactly these criteria, and nothing else." Worth stating explicitly in your Methods/Validation write-up as the actual, verified filter semantics (inclusive-of-missing), rather than the perhaps more intuitive but incorrect assumption (exclusive-of-missing).
- Test: filter a dataset with at least one sample having `sex = NA` down to a single sex category; confirm that `NA`-sex sample still appears in `filtered_meta()`'s result and is counted in `filter_summary_ui`.
- Fix: not applied.

### None of the above are "the code doesn't run" bugs. As with `mod_dataset.R`, the issues found here are entirely in the category of **silent, well-defined-but-non-obvious behavior on edge cases** — exactly the material a thesis validation chapter should surface and characterize, not necessarily patch.

---

## 8. VALIDATION TESTS

| # | Test | Expected behavior | How to verify |
|---|---|---|---|
| 1 | Open the Datasets tab with all 4 raw `.rds` files present | 4 cards, each with title/platform/sample+probe counts and a working NCBI GEO link | `file.exists()` each `RAW_DIR/<gse>_raw.rds`; open tab and compare displayed counts to `dim(get_raw_eset(gse))` |
| 2 | Open the Datasets tab where GSE93272's raw file is absent (this deployment) | Card shows role/used-in text plus an explicit "Raw file not found on disk" warning, not a crash | Confirmed structurally by code reading; re-verify live in the running app |
| 3 | Metadata tab: switch the picker across all 5 possible choices (4 GSE + uploaded, if available) | `understand_ui`/`meta_table_full` update correctly each time; sample/feature counts match `dim()` of the resolved dataset | Compare displayed valueBoxes to `nrow()`/`ncol()` computed directly via `load_individual_dataset()`/`dataset$expr` in an R console |
| 4 | Confirm picker independence (Issue O-3) | Selecting a non-default source on Metadata does not change Expression data's or QC's own picker | Switch Metadata to GSE15573, then open Expression data and confirm it still shows GSE93272 by default |
| 5 | QC -> Outliers: run with `mad_k = 3` (default) vs. `mad_k = 6` | Fewer samples flagged at `mad_k = 6` (looser threshold) than at `mad_k = 3`, on the same dataset | Run both, compare `n_flagged` in `qc_summary_ui` |
| 6 | QC -> Outliers: a dataset with a genuinely aberrant sample (very low correlation to cohort) | That sample appears in the flagged table with `reason` including "low cohort correlation" | Requires a dataset with a known technical outlier, or a synthetic one with an intentionally shuffled/noised column |
| 7 | QC -> Normalised data: run the check on an already-normalized (e.g. merged/batch-corrected) dataset | `needs_quantile_norm()` returns FALSE; green "Samples look well aligned..." message; box does not turn amber | Run against the "Merged Data" preloaded default, which is already normalized/batch-corrected per its own pipeline |
| 8 | QC -> Normalised data: run the check on a raw, linear-scale (un-logged) individual dataset, if one is reachable | `needs_quantile_norm()` returns TRUE (`max_value > 100` likely trips first); amber box; "This looks like it needs..." message | Requires a genuinely linear-scale reachable input — verify against GSE89408 (raw RNA-seq counts) specifically, since it is the one individual source most likely to be linear-scale |
| 9 | QC -> Normalised data: click Apply on an uploaded dataset, confirm before/after diagnostics improve | `diag_after$median_sd`/`iqr_sd` <= `diag_before`'s; valueBoxes turn green if now < 0.5 | Upload a synthetic un-normalized dataset, run both steps, compare `r$diag_before`/`r$diag_after` numerically |
| 10 | Adopt-normalization gate (Issue O-2): confirm "Use this normalised version" only appears for an uploaded, currently-inspected dataset | Button absent for any of the 4 fixed GEO sources; present only when `qc_source == "uploaded"` | Toggle the QC picker between a GSE id and "uploaded" (with something staged) and observe the panel change |
| 11 | Adopt-normalization effect: click "Use this normalised version," then check another sub-module (e.g. DGE) | `dataset$expr` reflects the normalized matrix immediately, without visiting Preprocessing at all | After clicking adopt, open any other Transcriptomics sub-module and confirm its input matrix matches `norm_apply_result()$expr_after` |
| 12 | Group tab: filter semantics with missing data (Issue O-10) | A sample with `NA` in the filtered column remains in the result even when its filter value would otherwise exclude it | Prepare/identify a dataset with at least one `NA` in a filterable column; filter to a specific level; confirm the `NA` row survives |
| 13 | Group tab: Reset after Apply | All filter inputs return to "show everything" state; a subsequent Apply reproduces the full, unfiltered result | Apply a narrowing filter, click Reset, click Apply again, compare row count to the dataset's full `nrow(meta)` |
| 14 | Group tab: `filter_spec()` cardinality rule | A metadata column with exactly 30 unique non-missing values is offered as a categorical filter; one with 31 is not; one with 1 unique value is not | Construct/identify metadata columns at these exact boundary cardinalities and confirm which appear in `output$filters` |
| 15 | High-cardinality color-by column (Issue O-8) | A `norm_color_by` selection with >7 unique levels shows fewer than that many distinct colors (some missing/grey) | Upload metadata with an 8+-level categorical column, select it, inspect the rendered legend |
| 16 | Preloaded-individual or GEO-fetched active dataset is invisible to this module's pickers (Issue O-5) | The "uploaded" 5th choice never appears for a non-upload-sourced active dataset, even though it differs from all 4 fixed sources | Activate a dataset via the preloaded-individual or GEO-fetch path (through Preprocessing), then open Overview's pickers and confirm only 4 choices exist |

---

## 9. THESIS DOCUMENTATION

### Methodological concept
An interactive, on-demand quality-control and data-exploration layer for a transcriptomics analysis pipeline, re-exposing the same missing-value audit, robust (median/MAD) sample-outlier detection, and normalization-need diagnostic the underlying automated pipeline computes internally — applied, at the user's choice, to any one of four fixed reference GEO cohorts or their own uploaded/activated dataset, entirely independent of the app's shared "active dataset" state except for one explicit, narrowly-scoped live-normalization adoption path.

### Implementation
A Shiny module (`mod_overview_ui`/`mod_overview_server`) built from four top-level tabs (one nested further into four QC sub-tabs), using three architecturally-identical, independently-scoped "Dataset to inspect" resolvers (a deliberate trade-off against a single shared picker, made because of how this module's UI is dynamically mounted via `insertTab()`); always-on `reactive()`s for cheap/always-relevant checks (missing values) and click-gated `eventReactive()`s for computationally heavier, opt-in checks (outlier detection, normalization diagnostics); one runtime-generic, data-driven filter-builder (`filter_spec()`) that adapts its own UI to whatever columns a given dataset's metadata actually contains, rather than hardcoding `group`/`sex`.

### Inputs
Any of: the four fixed `GEO_SOURCES` reference datasets (via `load_individual_dataset()`), or the app's currently **active** dataset (`dataset$expr`/`meta`), when it was produced by the upload path specifically. Never reads the **staged** (`dataset$staged_*`) fields the Dataset tab writes — this module only ever looks at what's already fully activated elsewhere.

### Processing
Missing-value counting (including common non-`NA` missingness sentinels: `""`, `"NA"`, `"N/A"`, `"unknown"`); robust (median/MAD-based) per-sample outlier flagging across three signals (total signal, detected-feature count, mean correlation to cohort, the last one-sided); `limma::normalizeBetweenArrays(method="quantile")` diagnostics and, on request, live application; `prcomp`-based PCA with a zero-variance-gene filter; a fully generic, cardinality- and type-aware metadata filter-spec builder.

### Outputs
Read-only tables, plots, and summary statistics for human inspection; one narrowly-scoped exception where a live-normalized version of the user's own uploaded, currently-active expression matrix can be adopted app-wide with a single click, bypassing the Preprocessing module's otherwise-exclusive activation gate.

### Reproducibility — what to document
- That `needs_quantile_norm()`/`summarize_norm_diagnostics()` are the **identical functions** (not reimplementations) the underlying pipeline's own preprocessing script uses (`global.R:1481-1543`'s own audit comment is a directly citable source for this claim).
- The exact MAD-based outlier rule and its `mad_k` parameterization (default 3, user-adjustable 2–6) — record whatever `mad_k` value was actually used for any outlier-flagging result you report.
- The exact `limma::normalizeBetweenArrays(method="quantile")` call and parameters, matching Batch Correction's own.
- The one activation-bypass path (live-adopt for uploaded data) as a documented, deliberate exception to the otherwise-single-gated (Preprocessing-only) activation architecture.

### Limitations to acknowledge
- The three tab-specific dataset pickers do not share a selection (Issue O-3) — a UX, not correctness, limitation.
- A preloaded-individual or live-GEO-fetched active dataset cannot be selected for inspection here as itself — only as one of the four fixed references, or (if it happens to also be the upload case) the "uploaded" option (Issue O-5).
- Metadata filters (Group tab) always let missing values through regardless of the filter applied (Issue O-10) — a defined, but easily misread, semantic.
- This module currently contributes nothing to the shared `results` object the AI assistant reads from other sub-modules (Issue O-9) — a disclosed architectural gap, not a hidden one, once you've read this file.
- The fixed 7-color categorical palette (`arthomix_pair()`) can silently under-represent metadata columns with more than 7 unique levels used as a plot color/fill (Issue O-8).

### Validation
Demonstrate via the tests in §8 — particularly: (a) that the same normalization-need verdict this tab shows for any given dataset matches what the pipeline's own batch-correction logic would independently compute for that same matrix (both call `needs_quantile_norm()`), and (b) the documented, testable edge-case behaviors (Issues O-3, O-5, O-8, O-10) characterized against controlled inputs, exactly as recommended for the Dataset module's own validation chapter.

---

## 10. CODE → THESIS MAPPING

| Code section | Functionality | Scientific purpose | Thesis section | What to document |
|---|---|---|---|---|
| L147-175 (`sources_ui`) | Renders the 4-source GEO catalog, live from disk | Establishes cohort provenance/traceability | Data / Methods — Study Cohorts | Which sources have raw files bundled in this deployment vs. not; platform/sample/probe counts per source |
| L177-244 (picker/resolver trio) | Independent per-tab dataset selection + resolution | Lets a user inspect any reference/active dataset in isolation, without merge/batch-correction applied | Implementation — Interactive QC Architecture | The `insertTab`-driven rationale for 3 independent pickers instead of 1 shared one; the "uploaded-only" recognition gap (Issue O-5) |
| L276-302 (`missing_audit`) | Always-on missing-value audit, 4 sentinel patterns | Establishes data completeness before any statistical test | Methods — Data Quality Control | The 4 missingness sentinels checked; per-field percent-missing methodology |
| L306-365 (`sample_qc`/`compute_sample_qc`) | MAD-based outlier detection across 3 metrics | Flags technical outliers using the same rule the pipeline itself uses | Methods — Sample-Level QC | The 3 metrics, the one-sided correlation rule, the `mad_k` sensitivity parameter and its default/range |
| L383-462 (`norm_check`) | Normalization-need diagnostic + PCA/scree/boxplot views | Decides, transparently, whether quantile normalization is warranted | Methods — Normalization Assessment | The exact 3-condition `needs_quantile_norm()` rule; its provable identity with the pipeline's own batch-correction decision logic |
| L471-546 (`norm_apply_result`/`adopt_norm_btn`) | Live quantile normalization + narrow app-wide adoption path | Lets a user actually fix an un-normalized upload without a separate Preprocessing detour | Implementation — Activation Architecture / Reproducibility | This is a second, deliberate activation path outside Preprocessing's usual single gate — document explicitly, with its scope restriction (uploaded data only) |
| L573-643 (`filter_spec`/`filtered_meta`) | Generic, data-driven metadata filtering | Lets a user explore any cohort subset without hardcoded assumptions about which columns matter | Implementation — Generic UI Generation | The type/cardinality rules deciding numeric-slider vs. categorical-picker vs. not-filterable; the NA-always-passes filter semantics (Issue O-10) |

---

## 11. PACKAGES AND DEPENDENCIES

| Package | Function(s) used | Why needed here | Core/optional |
|---|---|---|---|
| `shiny` | `NS`, `moduleServer`, `reactive`, `eventReactive`, `observeEvent`, `renderUI`, `uiOutput`, `req`, `validate`, `need`, `isTruthy`, `selectInput`, `sliderInput`, `checkboxInput`, `actionButton`, `downloadButton`/`downloadHandler`, `tagList`, `fluidRow`/`column`, `icon` | Core reactive framework, same role as in `mod_dataset.R` | Core |
| `shinydashboard` | `box`, `valueBox`, `tabPanel`/`tabsetPanel` styling conventions | Panel/KPI-tile widgets used throughout | Core |
| `shinycssloaders` | `withSpinner` | Loading-state UX for every heavier `uiOutput` | Optional (cosmetic) |
| `shinyWidgets` | `pickerInput`, `updatePickerInput` | Richer multi-select filter widget (Group tab) | Core to that one feature |
| `DT` | `DT::datatable`, `DT::renderDataTable`, `DT::dataTableOutput` | Every interactive table in this module — all server-side processed, since `server = TRUE` is `DT::renderDataTable()`'s own default (verified: `args(DT::renderDataTable)`), whether or not the call states it explicitly | Core |
| `ggplot2` (+ `theme_arthomix`/`arthomix_pair`/`ARTHOMIX_COLORS`/`ARTHOMIX_STATUS` from `global.R`) | `ggplot`, `geom_col`, `geom_boxplot`, `scale_fill_manual`, `coord_flip`, `stat_ellipse` (via `plot_pca_advanced`) | Every plot in this module | Core |
| `dplyr` | `%>%`, `count()` | `group_plot`/`dataset_plot`'s tally-then-plot idiom | Core to those two plots |
| `limma` | `limma::normalizeBetweenArrays` | Live quantile normalization (check + apply) | Core to the Normalised data sub-tab |
| `Biobase` | `Biobase::experimentData`, `Biobase::annotation` | Extracting title/platform metadata for the Datasets tab's cards | Core to that one sub-tab |
| `data.table` | `data.table::fwrite` | Fast CSV export of the (potentially large) expression matrix | Core (performance) |
| `WGCNA` (indirectly, via `global.R`'s `compute_sample_qc()`'s `cor()`) | correlation computation | Restored/overridden `cor <- WGCNA::cor` at `global.R:135` to avoid a masking conflict from Bioconductor packages loaded later | Core, indirectly |

### Dependencies on other files
- `global.R`: `GEO_SOURCES`, `get_raw_eset()`, `geo_link()`, `load_individual_dataset()`, `compute_sample_qc()`, `summarize_norm_diagnostics()`, `needs_quantile_norm()`, `pca_of()`, `scree_plot()`, `plot_pca_advanced()`, `qc_bar_plot()`, `theme_arthomix()`, `arthomix_pair()`, `ARTHOMIX_COLORS`, `ARTHOMIX_STATUS`, `%||%`.
- `server.R`: constructs the shared `dataset` object this file both reads and (narrowly) writes; registers this module via `TX_MODULES`'s generic `lapply(..., dataset, results)` call.
- `mod_dataset.R`: the file whose exact `"Uploaded dataset:"` source-string prefix this file's `qc_source_choices()` pattern-matches against (a cross-file string contract with no shared constant — Issue O-4).
- `mod_preprocessing.R`: the module whose "Use this as the active dataset" button is this app's *usual* sole gate for `dataset$expr`/`meta`/`source` writes — this file's `adopt_norm_btn` is the one documented exception.

---

## 12. ARCHITECTURE

```
global.R (GEO_SOURCES, get_raw_eset(), load_individual_dataset(), every shared
          QC/normalization/plotting helper this file calls)
        |
        v
server.R
  dataset <- reactiveValues(expr=, meta=, source=)  [from load_default_dataset() at startup]
  results <- reactiveValues()                        [shared cross-module AI-context store]
        |
        +--> lapply(TX_MODULES, function(m) m$server(paste0("tx_", m$config$id), dataset, results))
        |      includes: mod_overview_server("tx_overview", dataset, results)  <- THIS FILE
        |      (results is accepted but never used inside this file - Issue O-9)
        |
        +--> "Sub-modules" grid: this module's tab is NOT mounted by default -
             a user must explicitly click "Add" on its card (or trigger the
             sidebar-nav "Overview"/"Datasets" shortcuts, which pre-fill the
             grid's search box and switch to it if not yet added -
             server.R:335-353's jump_to_submodule()) before its tabPanel is
             inserted into "tx_menu" via insertTab() - same lifecycle as
             every other TX_MODULES entry, DESPITE this module's own name
             ("Overview and Datasets") suggesting it might be a fixed,
             always-visible landing page. It is not - only mod_dataset.R's
             "Dataset" tab is hardcoded directly into tx_menu (ui.R:1360).
        |
        v
  THIS FILE reads dataset$expr/meta/source (active dataset, uploaded-only
  recognition) and GEO_SOURCES/load_individual_dataset() (the 4 fixed refs);
  writes dataset$expr/source in exactly one place (adopt_norm_btn), bypassing
  mod_preprocessing.R's otherwise-sole activation gate for that one case.
```

- **Mounting behavior, code-verified (not previously stated in the companion Dataset-module notes):** `added$ids <- reactiveValues(ids = character(0))` (`server.R:144`) starts empty every session — this module's tab genuinely is not visible until added, exactly like DGE/WGCNA/every other analysis sub-module, despite its name and its "first in `TX_MODULES`" registry position perhaps suggesting otherwise. This is a concrete, useful fact for your Implementation/Architecture chapter: **"Overview and Datasets" is a toggleable sub-module, not a fixed landing page** — only "Dataset" (`mod_dataset.R`) is always present.
- The same three-part `config`/`ui`/`server` registration pattern, and the same generic `insertTab()`/`removeTab()` toggle mechanism, documented in the companion notes' §12 Architecture diagram applies identically here — this file adds no new registration mechanism of its own.

### Files to inspect next, and why
1. **`mod_preprocessing.R`** — to confirm precisely how the Batch Correction tab's own `needs_quantile_norm()`/`limma::normalizeBetweenArrays()` calls compare line-for-line to this file's (already strongly suggested to be identical by `global.R`'s own audit comment, but worth a direct side-by-side read for full confidence before asserting equivalence in your thesis).
2. **`ui.R`** (already partially inspected in this session, specifically `TRANSCRIPTOMICS_SIDEBAR_NAV` and `transcriptomicsUI()`) — confirms the "Dataset" vs. "Overview and Datasets" tab distinction stated in Issue O-1; read in full if you want to trace every sidebar-nav-to-tab mapping precisely.
3. **`R/submodules_registry.R`** (already read in part) — confirms `mod_overview`'s position (first) in `TX_MODULES` and the shared shape every sub-module's config/ui/server trio follows.
4. **`AUDIT.md`** / project thesis outline (`thesis/chapter4_outline.md`, already in this repository) — for how (if at all) this module's absence from the `results` contract (Issue O-9) should be framed alongside that document's broader "which submodules are fully built vs. scaffolded" disclosure discipline — Overview is fully built, just architecturally silent toward the assistant, a distinct category worth naming precisely.

---

## 13. FINAL LEARNING SUMMARY

### A. What this file does (plain English)
It's the app's dataset-inspection and QC dashboard: browse where the four reference GEO cohorts come from, look at any one of them (or your own active data) completely unfiltered — full metadata table, full expression matrix — and run the same missing-value, outlier-detection, and normalization-need checks the underlying research pipeline itself runs, entirely on demand, with results shown only to the user (except for one narrow, explicit path that lets a live-normalized version of your own uploaded data become the app's active dataset directly).

### B. What you learned (R/Shiny concepts), beyond what the Dataset module already taught
- Why `renderUI`-wrapped `selectInput`s (re-rendered fully on every relevant change) are sometimes the *correct* choice over `update*Input()`, specifically because of *when* a dynamically-`insertTab()`-mounted module's server code actually starts running relative to the DOM elements it targets — and, in the Group tab's `reset_btn`, a contrasting case where `update*Input()` genuinely is safe, because the inputs it targets are guaranteed already-live.
- Fully data-driven/generic UI generation: building a *variable number* of inputs, with *dynamically constructed* IDs (`paste0("f_", cl)`), whose *type* (slider vs. multi-select) is chosen at runtime by inspecting the data itself (`is.numeric()`, cardinality) — and reading them back via computed-name (`input[[...]]`) access.
- `ggplot2`'s tidy-eval `.data[[var]]` pronoun for column-name-by-variable plotting.
- The distinction between an always-on `reactive()` (cheap, always relevant) and a click-gated `eventReactive()` (expensive or intentionally opt-in) as a *cost-driven* design choice, not an arbitrary one.
- Robust statistics (median/MAD) as a deliberate, justified alternative to mean/SD for outlier detection, and why the choice matters specifically when detecting the very outliers that would otherwise inflate a mean/SD-based rule.
- Why you cannot infer a package function's runtime behavior (here, `DT::renderDataTable()`'s server-side-vs-client-side processing) from whether a call states an argument explicitly — always check the documented default (`args()`/`?function`) before drawing a conclusion from its absence.

### C. What you should be able to explain at your defense
- Exactly which checks in this module are "the same code, reused" versus "a separate, potentially divergent reimplementation" relative to the underlying pipeline's own automated preprocessing (answer: `needs_quantile_norm()`/`summarize_norm_diagnostics()` are literally the same functions, confirmed by direct code reading, not by trusting a comment).
- Why this module deliberately uses three independent dataset pickers instead of one shared control, and what the concrete, user-visible consequence of that choice is (Issue O-3).
- The precise, narrow circumstances under which this module is allowed to change what every other sub-module analyzes (Issue O-2), and why that's scoped the way it is.
- The exact semantics of the Group tab's filters with respect to missing data (Issue O-10) — and why that's a defensible design choice, not an obviously wrong one, even though it may not be every user's first assumption.
- That this module is a **toggleable sub-module**, not a fixed landing page, despite its name — and what that implies about when a user actually sees it during a session.

### D. Validation status
- Correct: MAD-based outlier rule and its degenerate-case guard; the provable code-identity between this tab's normalization-need check and the pipeline's own; the `merge(..., all.x=TRUE)` sample-alignment join; the `renderUI`-for-dynamic-inputs pattern and its documented rationale.
- Corrected in this document (originally misstated, fixed after checking `args(DT::renderDataTable)` directly): every `DT` table in this module is server-side processed by default — there is no actual client/server split between `expr_table` and the other five tables, only a difference in whether the (identical) default is stated explicitly.
- Design trade-off, not a bug: three independent, non-synchronized dataset pickers (Issue O-3); the narrow live-adopt activation bypass (Issue O-2); NA-always-passes filter semantics (Issue O-10).
- Gap worth disclosing: this module contributes nothing to the shared `results`/AI-assistant context (Issue O-9); a preloaded-individual or GEO-fetched active dataset can't be selected for inspection as itself (Issue O-5); the fixed 7-color palette under-represents high-cardinality categorical columns (Issue O-8).
- Minor/cosmetic: a stale server-code comment mislabels which tab it implements (Issue O-1); a download-filename fallback doesn't cover the `NA`-value case (Issue O-7); a few numeric thresholds/string prefixes are duplicated rather than centralized (Issue O-4).
- Untested/needs live verification: the exact live-app behavior for Issues O-5, O-8, O-10, and O-2's stale-cache edge case, all of which are reasoned from code reading and should be exercised in the running app before being asserted as findings in your thesis proper.

### E. Questions to investigate further
1. Does `mod_preprocessing.R`'s own normalization-check UI call the exact same `norm_check`-equivalent logic, confirming full parity, or does it differ in any parameter?
2. Is `mod_overview_config`'s `group = "Data"` field read anywhere in the app (e.g., a grouped Sub-modules grid layout) or is it currently inert metadata? (**Cannot confirm from this file alone.**)
3. Does any other sub-module's `qc_tabs`-style nested tabset get addressed via `jump_to_submodule()`'s `inner_tab` parameter — i.e., is there a sidebar shortcut planned/possible for, say, jumping straight to QC → Outliers, that simply doesn't exist yet?
4. For Issue O-5: was the decision to only recognize the "uploaded" active-dataset case deliberate (e.g., because a preloaded-individual or GEO-fetched active dataset is *already* one of the choices this module offers directly, just not automatically pre-selected) — worth re-examining with fresh eyes, since a GEO-fetched *new* series (not one of the 4 fixed sources) genuinely has no dedicated entry point here at all.
5. Confirm live, in the running app, whether a >7-level categorical column (Issue O-8) produces a `ggplot2` warning, a silently missing color, or an error — the exact failure mode matters for how you'd phrase this in a Limitations section.

### F. Next file to study
Given this file and its companion (`mod_dataset_teaching_notes.md`) together cover the full "get data in, then look at it" lifecycle, the natural next step is **`mod_preprocessing.R`** — the file both this module and `mod_dataset.R` explicitly point the user toward ("go to Preprocessing and pick..."), and the only place `dataset$expr/meta/source` is normally allowed to change. It is also where you can directly verify the "same code, not a reimplementation" claim this document makes about `needs_quantile_norm()`/`summarize_norm_diagnostics()`/`limma::normalizeBetweenArrays()`.

---

## 14. COMBINED "OVERVIEW AND DATASETS" AREA — MODULE-LEVEL SUMMARY (spanning both `mod_dataset.R` and `mod_overview.R`)

Per the sidebar's own grouping comment (`ui.R:1327-1328`: *"Overview and Datasets is one outer tx_menu tab"*), and per this session's direct code reading, these two files jointly implement what a user experiences as **one conceptual area** even though they are two separate `mod_*.R` files with two separate mounting lifecycles. This section answers the brief's request for a module-level synthesis across both.

### Overall purpose of the Dataset module (`mod_dataset.R`)
The single ingestion point: converts a bundled reference cohort, an arbitrary user upload, or a live NCBI GEO fetch into one harmonized `{expr, meta}` pair, staged (never directly active) for Preprocessing to pick up.

### Overall purpose of the Overview module (`mod_overview.R`)
The inspection/QC layer: lets a user browse cohort provenance and run the pipeline's own missing-value/outlier/normalization diagnostics, interactively, against any fixed reference dataset or the currently-active dataset — read-only except for one narrowly-scoped live-normalization adoption shortcut.

### How their sub-tabs connect
- `mod_dataset.R` has no internal sub-tabs (three side-by-side boxes: preloaded / GEO fetch / upload) feeding one shared `staged_*` output.
- `mod_overview.R`'s four sub-tabs (Datasets / Metadata / Expression data / QC) are mutually independent **views**, not a pipeline — a user can open any one directly, in any order, and each resolves its own "Dataset to inspect" selection separately (§2.0).
- The genuine cross-file connections are: (1) `mod_overview.R`'s "Datasets" tab and `mod_dataset.R`'s preloaded catalog both ultimately read `GEO_SOURCES`/`get_raw_eset()`/`load_individual_dataset()` from `global.R` — the same single source of truth for "what reference cohorts exist"; (2) `mod_overview.R`'s uploaded-data recognition (`qc_source_choices()`) depends on the exact source-string convention `mod_dataset.R` (via Preprocessing's activation) establishes; (3) `mod_overview.R`'s one write to `dataset$expr` (adopt-normalize) is a second, narrower version of the same "make this the active dataset" action Preprocessing otherwise exclusively performs.

### Complete user workflow, end to end
```
1. (mod_dataset.R) Pick a source: preloaded / upload / GEO fetch -> preview -> "Load" -> dataset$staged_*
2. (mod_preprocessing.R, outside this document's scope) "Currently loaded dataset" reads staged_* ->
   merge / normalize / batch-correct -> "Use this as the active dataset" -> dataset$expr/meta/source
3. (mod_overview.R) Independently, at any point: browse GEO provenance (Datasets tab); inspect the
   active dataset OR any of the 4 fixed raw sources in full (Metadata/Expression data tabs); run
   missing-value/outlier/normalization QC on-demand (QC tab); optionally filter by any metadata
   column to explore a subset (QC -> Group)
4. (mod_overview.R, narrow exception) If the active dataset is the user's own upload and QC ->
   Normalised data flags it as needing normalization, the user may fix it in place, right here,
   without returning to Preprocessing at all
```

### Complete data flow, end to end
```
Bundled .rds/.csv (global.R paths) __or__ user-uploaded files __or__ GEOquery::getGEO()
        --> mod_dataset.R: format-detect, parse, harmonize columns, intersect+align samples
        --> dataset$staged_expr/staged_meta/staged_source
        --> mod_preprocessing.R (outside this document): merge/normalize/batch-correct
        --> dataset$expr/meta/source  [THE ACTIVE DATASET]
        --> every other Transcriptomics sub-module reads this
        --> IN PARALLEL, independent of the above: mod_overview.R reads dataset$expr/meta
            (if uploaded) OR global.R's GEO_SOURCES/load_individual_dataset() (the 4 fixed refs)
            for its own read-only inspection views
        --> ONE EXCEPTION: mod_overview.R's adopt_norm_btn writes dataset$expr/source directly,
            re-entering the same "active dataset" pool from a second path
```

### Key functions and dependencies shared by both files
`global.R`: `GEO_SOURCES`, `get_raw_eset()`, `load_individual_dataset()`, `load_default_dataset()`, `collapse_probes_to_genes()`, `%||%`, `theme_arthomix()`/`arthomix_pair()`/`ARTHOMIX_COLORS` (Overview only, for plots), `compute_sample_qc()`/`summarize_norm_diagnostics()`/`needs_quantile_norm()`/`pca_of()` (Overview only). Both files depend on the same shared `dataset` `reactiveValues` object constructed once in `server.R`.

### Important validation findings across both files (aggregate)
- **Correctness is generally sound in both files.** No finding in either file's validation section is a "produces wrong numbers" bug — every finding is a **silent, well-defined edge-case behavior**: heuristic column-guessing with no confidence signal (Dataset, Issue 4), no numeric-type validation on uploads (Dataset, Issues 1-2), non-synchronized dataset pickers (Overview, Issue O-3), NA-inclusive filter semantics (Overview, Issue O-10), and two real, disclosed **functional gaps**: this app's `results`/AI-assistant contract is not populated by Overview (Issue O-9), and a preloaded-individual or live-GEO-fetched active dataset has no direct inspection entry point in Overview (Issue O-5).
- **The strongest reproducibility claim available from direct code reading**, spanning both files: the normalization-need decision a user sees interactively in Overview's QC tab is the *same function call*, not a parallel reimplementation, as the one the underlying research pipeline's own batch-correction step uses — this is directly citable and independently verifiable by anyone re-reading `global.R:1481-1543`.
- **The most important data-provenance caveat**, established in the Dataset module and directly relevant whenever Overview's "Datasets"/QC tabs are used against GSE93272/GSE110169: those two "raw" sources, in this deployment, are not actually raw — they are `merged_training_subset()` fallbacks derived from the already-merged/batch-corrected cohort. Any QC finding you report against "GSE93272 individually" must carry this caveat.

### What should be documented in the thesis
1. The staged-vs-active dataset lifecycle (Dataset module) and the fact that Overview is a *separate*, largely read-only lifecycle layered on top of it, with one explicit, scoped exception.
2. The provable code-identity between Overview's interactive QC and the pipeline's own automated preprocessing decisions.
3. The GSE93272/GSE110169 raw-vs-fallback data-provenance caveat, wherever either file's handling of those two sources is discussed.
4. Every disclosed limitation above, framed — consistent with this project's own `AUDIT.md` precedent — as an honest, examined finding rather than an omission.

### Which parts of the implementation require additional testing
- Every item in each file's own §8 Validation Tests table, especially: the order-alignment invariant (Dataset), the MAD-outlier boundary behavior and the >7-category color palette (Overview), and both files' documented-but-not-live-verified edge cases (upload with non-numeric data; a GEO-fetched active dataset's invisibility to Overview's pickers).
- A live walkthrough of the exact user journey in §"Complete user workflow" above, on this specific deployment's actual bundled data, to confirm every count/label this document predicts from code reading (sample counts, feature counts, which raw files exist) matches what the running app actually shows — code reading establishes *what the code will do*, not a substitute for *having run it*.
