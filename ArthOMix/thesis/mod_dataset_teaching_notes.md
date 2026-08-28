# `mod_dataset.R` — Full Teaching, Audit, and Thesis-Documentation Notes

File: `ArthOMix/R/transcriptomics/mod_dataset.R` (452 lines)
Prepared: 2026-08-25

---

## 1. PURPOSE

### Why this file exists

`mod_dataset.R` is the **Transcriptomics → Dataset tab**. It is the single entry point through which any expression matrix + sample metadata pair — preloaded, uploaded, or fetched live from NCBI GEO — gets into the app in a shape every other transcriptomics sub-module can use (a numeric `expr` matrix, genes × samples, and a `meta` data frame with at least `sample`, `group`, `sex`, `batch` columns).

### The problem it solves

A multi-omics Shiny app built around one training cohort still needs to let a user (a) come back to that default cohort, (b) look at any of the four raw GEO sources on their own, (c) bring their own expression/metadata files, or (d) pull a new GEO series directly by accession — all through one consistent interface, without those four paths turning into four different downstream data shapes. This file is the harmonization layer that makes all four paths converge on the same `expr`/`meta`/`source` triple.

### Role in the app / module hierarchy

- App: ArthOMix Explorer (`server.R`, `global.R`, `ui.R`)
  - Top-level omics area: **Transcriptomics** (`R/transcriptomics/`)
    - Sub-module: **Dataset** (`mod_dataset.R`) ← this file
    - Sibling sub-modules that *read* what this one stages: Preprocessing (`mod_preprocessing.R`), and, once activated, DGE/WGCNA/MR/etc.

### Inputs it receives

- `id` — the Shiny module namespace string (`"tx_dataset"`, set in `server.R:76`).
- `dataset` — a `reactiveValues` object created once in `server.R:7-10`, shared by reference across the whole Transcriptomics area. This file both reads (`dataset` is passed in) and writes to it (`dataset$staged_expr` etc.).
- User inputs: uploaded files, a preloaded-dataset dropdown selection, a GEO accession string, and several column-mapping dropdowns.
- Package-level/global data: `GEO_SOURCES`, `load_default_dataset()`, `load_individual_dataset()`, `collapse_probes_to_genes()` — all defined in `global.R` and available because `global.R` is sourced before any module file runs.

### Outputs it produces

- Three fields on the shared `dataset` object: `dataset$staged_expr`, `dataset$staged_meta`, `dataset$staged_source`. **Never** `dataset$expr`/`meta`/`source` directly — that distinction is the whole design of this file (see the file's own header comment, lines 1-14).
- UI feedback: preview text (rows/columns read), a green success message, or a red error message, all rendered inline rather than as Shiny's default red error screen.

### Dependencies (other files/functions this module needs)

- `global.R`: `load_default_dataset()`, `load_individual_dataset()`, `GEO_SOURCES`, `collapse_probes_to_genes()`, the `%||%` null-coalescing helper.
- Packages: `shiny` (core), `shinyjs` (enable/disable buttons), `data.table` (fast CSV reads), `GEOquery` + `Biobase` (GEO fetch, optional at runtime), `WGCNA` (indirectly, inside `collapse_probes_to_genes()`).
- `server.R`: constructs `dataset` and calls `mod_dataset_server("tx_dataset", dataset)`.

### Downstream consumers

- `mod_preprocessing.R:145-152` and `:295-305`: reads `dataset$staged_expr %||% dataset$expr` (and the `_meta`/`_source` equivalents) — i.e. "prefer whatever was just previewed on the Dataset tab; otherwise fall back to whatever is currently active app-wide."
- `mod_preprocessing.R:1928-1938`: the **only** place in the app that ever writes to `dataset$expr`/`meta`/`source` after the app's initial load — the "Use this as the active dataset" button, gated behind Preprocessing's merge/batch-correction pipeline.
- `server.R:19-21, 28-41`: `dataset$source` drives the header badge and clears stale cross-module `results` whenever the *active* dataset changes (not the staged one).

### What would happen if this file stopped working

- No way to preview an individual raw GEO source, upload your own data, or fetch a new GEO series. The app would still boot fine on the default merged cohort (that load happens in `server.R`/`global.R`, independent of this file), but Preprocessing's "Currently loaded dataset" option would have nothing new to fall back on beyond the already-active dataset, since `dataset$staged_*` would simply never be set.
- If `mod_dataset_server()` errored during `moduleServer()` setup, the whole Shiny session would likely fail to start (module server functions run at session-init time), taking the entire app down, not just this tab.

### Conceptual workflow

```
User picks a source (preloaded dropdown / file upload / GEO accession)
        ↓
UI captures it as a reactive input (input$preloaded_choice / input$expr_file+meta_file / input$geo_accession)
        ↓
Server reads/parses the underlying data (reactive: expr_raw()/meta_raw(), or entry$load(), or GEOquery::getGEO())
        ↓
Preview renders immediately (upload_preview_ui / geo_fetch_status) — no button needed yet
        ↓
User maps metadata columns (sample ID / group / sex / batch) via guess_col()-defaulted dropdowns
        ↓
User clicks a "Load"-style button (load_preloaded_btn / load_btn / geo_load_btn)
        ↓
observeEvent validates, harmonizes column names, intersects sample IDs, writes dataset$staged_expr/meta/source
        ↓
Success/error message renders (output$load_message or output$preloaded_load_message)
        ↓
User is told to go to Preprocessing and pick "Currently loaded dataset" to actually use it
```

---

## 2. UI

`mod_dataset_ui(id)` (lines 79-130) builds three `box()` panels inside a two-column `fluidRow()`, plus one shared status area at the bottom.

| UI element | `inputId` (post-`ns()`) | What it is | What the user does | Server code that reacts |
|---|---|---|---|---|
| `selectInput(ns("preloaded_choice"), ...)` | `tx_dataset-preloaded_choice` | Dropdown of preloaded datasets | Pick "Merged Data" or one of 4 GSE cohorts | `output$preloaded_note` (renderUI), `observeEvent(input$load_preloaded_btn)` |
| `uiOutput(ns("preloaded_note"))` | — (output, not input) | Contextual info/warning box | Read-only, updates automatically | populated by `output$preloaded_note <- renderUI(...)` (L136-148) |
| `actionButton(ns("load_preloaded_btn"), ...)` | `tx_dataset-load_preloaded_btn` | "Load this dataset" button | Click to stage the selected preloaded dataset | `observeEvent(input$load_preloaded_btn)` (L242-255) |
| `uiOutput(ns("preloaded_load_message"))` | — | Inline success text next to that button | Read-only | set inside the same `observeEvent` above |
| `textInput(ns("geo_accession"), ...)` | `tx_dataset-geo_accession` | Free-text GEO Series accession box | Type e.g. `GSE12345` | read inside `geo_fetch_result` (L307-322) |
| `actionButton(ns("geo_fetch_btn"), ...)` | `tx_dataset-geo_fetch_btn` | "Fetch" button | Click to hit NCBI GEO | triggers `eventReactive(input$geo_fetch_btn, ...)` (L307) |
| `uiOutput(ns("geo_fetch_status"))` | — | Fetch result / error / preview | Read-only | `output$geo_fetch_status` (L364-380) |
| `uiOutput(ns("geo_platform_ui"))` | — (contains a dynamically-created input `geo_platform_choice`) | Platform picker, only shown if the series spans >1 platform | Pick a platform | `geo_eset` reactive (L332-341) |
| `uiOutput(ns("geo_column_mapping"))` | — (contains `geo_map_group`/`geo_map_sex`/`geo_map_batch`) | Column-mapping dropdowns for the fetched metadata | Confirm/adjust auto-guessed column mapping | `observeEvent(input$geo_load_btn)` (L409-448) |
| `actionButton(ns("geo_load_btn"), ...)` | `tx_dataset-geo_load_btn` | "Load this dataset" (GEO) | Click to stage the fetched dataset | `observeEvent(input$geo_load_btn)` (L409-448) |
| `fileInput(ns("expr_file"), ...)` | `tx_dataset-expr_file` | Expression matrix upload (.csv/.rds) | Choose a file from disk | `expr_raw` reactive (L165-177) |
| `fileInput(ns("meta_file"), ...)` | `tx_dataset-meta_file` | Sample metadata upload (.csv/.rds) | Choose a file from disk | `meta_raw` reactive (L150-160) |
| `uiOutput(ns("upload_preview_ui"))` | — | "Read N genes x M samples..." preview | Read-only, fires as soon as both files parse | L182-196 |
| `uiOutput(ns("column_mapping"))` | — (contains `map_id`/`map_group`/`map_sex`/`map_batch`) | Column-mapping dropdowns for the uploaded metadata | Confirm/adjust mapping | `observeEvent(input$load_btn)` (L257-298) |
| `actionButton(ns("load_btn"), ...)` | `tx_dataset-load_btn` | "Upload Data" — disabled until ready | Click once mapping is complete | `observeEvent(input$load_btn)`; enabled/disabled by `observe()` (L236-240) |
| `uiOutput(ns("load_message"))` | — | Shared success/error banner for both the upload path and the GEO path | Read-only | set in L279-297 and L434-447 |

### Dynamic UI mechanics

Almost every informational and mapping element here is a `uiOutput()`/`renderUI()` pair rather than static UI, because what should appear depends on data that only exists after a reactive step (a file was read, GEO was queried) — you cannot know the metadata's column names, or whether a GEO series has multiple platforms, at UI-definition time (`mod_dataset_ui()` runs once, before any session exists).

- `uiOutput(ns("x"))` in the UI function reserves a `<div>` placeholder tagged with the namespaced ID `tx_dataset-x`.
- `output$x <- renderUI({ ... })` in the server function is the reactive expression whose return value (a `tagList`/HTML tag) is pushed into that placeholder over websocket, any time a reactive value it reads (e.g. `input$meta_file`, or a reactive like `meta_raw()`) changes.
- `NS(id)` (line 80, `ns <- NS(id)`) returns a function that prefixes every ID with the module's `id` plus a dash, e.g. `ns("expr_file")` → `"tx_dataset-expr_file"`. This is Shiny's module namespacing: it lets `mod_dataset_ui`/`mod_dataset_server` be instantiated more than once in the same app (or reused in another app) without input/output ID collisions, because every ID is scoped under the caller-supplied `id`.
- Inside the server function, `ns <- session$ns` (line 134) is the *same kind* of function but derived from the running session rather than passed the literal `id` string — both produce identical namespacing; using `session$ns` inside the server is the idiomatic pattern because the server function doesn't have direct access to the original `id` argument's closure in all call patterns (here it technically does, since `id` is a formal argument of `mod_dataset_server`, but `session$ns` is the standard/portable convention taught by the Shiny docs and used consistently in this codebase).
- `input$whatever` inside the server function automatically resolves to the namespaced input, e.g. writing `input$expr_file` inside `mod_dataset_server` reads the value the browser sent for `tx_dataset-expr_file` — the module server never has to call `ns()` on `input$`/`output$` accesses, only when it constructs new IDs to hand back to the UI (as in `selectInput(ns("map_id"), ...)` inside a `renderUI()`).

Two dropdowns are themselves generated *inside* a `renderUI()` — `column_mapping` (L213-230) and `geo_column_mapping` (L386-401) — because their `choices` are literally the column names of a data frame that doesn't exist until a file is uploaded or GEO responds. This is a second, nested level of dynamism: `uiOutput` → `renderUI` produces more `selectInput`s, which themselves then generate new `input$...` reactive values (`input$map_id`, etc.) that did not exist in the session before that render happened.

---

## 3. FUNCTIONALITY (by logical section)

### 3.1 Module-level constants and helper builders (L16-77)

- **What:** A config list, a named character vector, two "entry" builder functions/objects, a combined list, and a formatter function — all defined **outside** `mod_dataset_server()`, at the top level of the file (i.e., evaluated once, when the app sources this file, not once per user session).
- **Why:** These are pure, session-independent data/functions: the catalog of what preloaded datasets exist doesn't change per user, so there's no reason to rebuild it inside `moduleServer()` on every session start. Defining them at file scope also lets `mod_preprocessing.R` reuse `preloaded_choices()` and `default_dataset_entry$id` (confirmed at `mod_preprocessing.R:126-128, 155`) — a private-inside-`moduleServer()` definition wouldn't be visible outside this file.
- **What it receives:** `GEO_SOURCES` (global.R), nothing else at define-time.
- **What it does:** Builds a list of five "loadable dataset" descriptors, each with `id`, `label`, and a `load()` closure that returns `list(expr=, meta=, source=)`.
- **What it returns:** `PRELOADED_DATASETS` (a list of 5 entries) and `preloaded_choices()` (a named vector suitable for `selectInput(choices=)`).
- **What depends on it:** `mod_dataset_ui()`'s `selectInput(ns("preloaded_choice"), choices = preloaded_choices())`; the `observeEvent(input$load_preloaded_btn)` handler, which does `Find(function(d) d$id == input$preloaded_choice, PRELOADED_DATASETS)`; and externally, `mod_preprocessing.R`'s `pp_cohort_choices()`.

### 3.2 UI builder — `mod_dataset_ui()` (L79-130)

- **What:** A plain R function returning a `tagList()` of HTML-generating Shiny tag functions. Not reactive — it runs once per session at UI-render time (technically, in this app's `ui.R`, `mod_dataset_ui("tx_dataset")` is likely called once when the whole app's UI is assembled — confirm by checking `ui.R`, see §12).
- **Why:** Shiny modules split UI generation (`*_ui`) from server logic (`*_server`) so both can be namespaced consistently and reused/tested independently.
- **What it receives:** `id`, the string that becomes this instance's namespace.
- **What it does:** Declares three input boxes (preloaded picker, GEO fetch, file upload) and one shared message area.
- **What it returns:** A nested tag tree (ultimately HTML) that Shiny inserts into the page.
- **What depends on it:** Whatever calls `mod_dataset_ui("tx_dataset")` in the app's top-level UI definition.

### 3.3 Server logic — `mod_dataset_server()` (L132-451)

This is one `moduleServer()` call containing:

1. **Reactive readers** (`meta_raw`, `expr_raw`, `geo_fetch_result`, `geo_eset`, `geo_expr_meta`) — lazily re-evaluated caches around parsing/fetching, so the same expensive operation (parsing a large CSV, hitting NCBI) isn't repeated by every UI element that needs the result.
2. **`renderUI` blocks** — pure presentation, driven by the reactives above and by raw `input$...` values.
3. **One `observe()` per form** (upload form, GEO form) — side-effect-only blocks that toggle button enabled/disabled state.
4. **Three `observeEvent()` "commit" handlers** (`load_preloaded_btn`, `load_btn`, `geo_load_btn`) — the only places that actually mutate the shared `dataset` object. Everything before this point is preview-only and has zero effect on any other module.
5. **One pure helper**, `guess_col()` — ordinary R function, not reactive, called from inside three different `renderUI` blocks to pick a sensible default dropdown selection.

---

## 4. LINE-BY-LINE TEACHING

### Block A — Header comment (L1-14)

```r
## R/mod_dataset.R
## The Dataset tab: lets the user PREVIEW a preloaded dataset, an uploaded
## file pair, or a live GEO fetch - a candidate, not yet the app's active
## analysis dataset. ...
```
**What this does:** Nothing executable — `##` starts a comment, R ignores the rest of the line. **Why it exists:** Documents the single most important architectural fact about this file (staged vs. active dataset) at the point someone is most likely to read it before editing. **R concept:** Comments (`#`); this codebase's convention is `##` for block/explanatory comments. **Thesis relevance:** High — this comment *is* your methodology's data-governance rule in prose; you can cite it near-verbatim when describing why "preview" and "activate" are deliberately separate operations (traceability / no-silent-mutation design).

### Block B — Module registration list (L16-19)

```r
mod_dataset_config <- list(
  id = "dataset", title = "Dataset", icon = "database",
  description = "Pick a preloaded dataset or upload your own - either way it's what every sub-module below reads from."
)
```
**What it does:** Creates an ordinary named `list`. **Why:** This is almost certainly consumed by a module registry (`TX_MODULES`, referenced in your thesis outline) that iterates over all transcriptomics sub-modules to build the sidebar/tab list, using `id`/`title`/`icon`/`description` uniformly across every sub-module file. **R concept:** `list()` with named elements — R's general-purpose key-value container; elements accessed via `$` (e.g. `mod_dataset_config$title`). **Data flow:** Static metadata → consumed at app-startup by whatever builds the navigation UI (not visible in this file — "Cannot determine from this file alone" exactly which function reads it; check `submodules_registry.R`, per your own thesis outline's citation of that file). **Thesis relevance:** Medium — useful for describing the app's declarative module-registry architecture, but not itself a computational step.

### Block C — Preloaded dataset labels and catalog (L21-77)

```r
INDIVIDUAL_DATASET_LABELS <- c(
  "GSE93272"  = "Whole Blood Training Cohort A",
  "GSE110169" = "Whole Blood Training Cohort B",
  "GSE15573"  = "PBMC Validation Cohort",
  "GSE89408"  = "Synovial Tissue Validation Cohort"
)
```
**What this does:** Builds a *named character vector* — the names (`"GSE93272"`, ...) are GEO accessions, the values are human-readable labels. **Why:** UI-facing labels are decoupled from GEO accessions so the dropdown reads as biology, not database IDs. **R concept:** `c(name = value, ...)` — the standard way to build a named vector; lookup by name works like a hash map: `INDIVIDUAL_DATASET_LABELS[["GSE93272"]]` → `"Whole Blood Training Cohort A"`. **Data flow:** Static constant → read inside `individual_dataset_entry()`. **Thesis relevance:** Low computationally, but directly relevant to your Methods section's cohort table (accession ↔ role ↔ tissue), since it's a second, independent place in the code that must stay consistent with `GEO_SOURCES` in `global.R` — worth explicitly cross-checking (see §7 Validation).

```r
individual_dataset_entry <- function(gse_id) {
  list(
    id = gse_id,
    label = if (gse_id %in% names(INDIVIDUAL_DATASET_LABELS)) INDIVIDUAL_DATASET_LABELS[[gse_id]] else gse_id,
    load = function() {
      d <- load_individual_dataset(gse_id)
      validate(need(!is.null(d), paste("Raw data for", gse_id, "was not found on disk.")))
      list(expr = d$expr, meta = d$meta, source = paste0("Individual dataset: ", d$label))
    }
  )
}
```
**What this is:** A *factory function* — given one GSE ID, it returns a 3-field descriptor list, where the third field (`load`) is itself a function (a **closure**) that "remembers" `gse_id` from its enclosing environment even after `individual_dataset_entry()` has returned.
**Why:** Deferred/lazy loading — building `PRELOADED_DATASETS` (below) must not eagerly read every GEO source's raw file from disk at app-startup; `load()` should only run when the user actually clicks "Load this dataset."
**R concept — closures:** In R, a function created inside another function captures its parent's environment. Here, every time `individual_dataset_entry("GSE93272")` is called, a *new* `load` function is created whose `gse_id` is permanently bound to `"GSE93272"`, distinct from any other `load` function built from a different call. This is what lets `PRELOADED_DATASETS` hold five different closures that each read a different file when eventually invoked.
**`validate(need(...))`:** A Shiny-specific idiom. `need(cond, message)` returns `NULL` if `cond` is `TRUE`, or a `shiny.silent.error`-carrying object with `message` if `cond` is `FALSE`/`NA`. `validate()` then either does nothing (cond met) or stops execution of the current reactive context and shows `message` in place of output, *without* triggering Shiny's generic red error screen and without a `try`/`tryCatch` block being visible at the call site. It only works inside a reactive context (an `observe`, `render*`, `reactive`, `eventReactive`) — calling it in plain code would just return an object nobody catches.
**What it returns:** `list(id=, label=, load=)`.
**Data flow:** `gse_id` string → looked up in `INDIVIDUAL_DATASET_LABELS` for a label → wrapped in a closure that, when called later, invokes `load_individual_dataset()` (global.R) and reshapes its result.
**Example:** `entry <- individual_dataset_entry("GSE15573"); entry$label` → `"PBMC Validation Cohort"`; `entry$load()` only executes `load_individual_dataset("GSE15573")` at the moment you call `entry$load()`, not at the moment `entry` was created.
**Thesis relevance:** High for a Software Architecture / Implementation subsection — this is a clean example of lazy evaluation via closures used specifically to avoid unnecessary disk I/O, worth a sentence of justification.

```r
default_dataset_entry <- list(
  id = "__default_merged__",
  label = "Merged Data",
  load = function() {
    d <- load_default_dataset()
    list(expr = d$expr, meta = d$meta, source = d$source)
  }
)
```
**What this does:** Same shape as an `individual_dataset_entry()` output, but hand-built (not via the factory) because it wraps `load_default_dataset()` instead. **Why a distinct sentinel `id`:** `"__default_merged__"` is deliberately not a real GEO accession, so it can never collide with a `GSE...` entry and is unambiguous to test against with `identical()` elsewhere (`mod_dataset.R:138`, `mod_preprocessing.R:155`). **R concept:** double-underscore-wrapped string as a private sentinel/constant — a plain convention, not a language feature; R has no true "private" identifiers here. **Thesis relevance:** Medium — this `id` string is effectively a magic constant duplicated in two files (`mod_dataset.R` and `mod_preprocessing.R:155`); worth a footnote acknowledging the coupling (see §7).

```r
PRELOADED_DATASETS <- c(
  list(default_dataset_entry),
  lapply(vapply(GEO_SOURCES, `[[`, character(1), "gse"), individual_dataset_entry)
)
```
**What this does:** Builds the final list of 5 entries: the merged default first, then one entry per GEO source in `GEO_SOURCES`.
**R concept — `vapply()`:** `vapply(X, FUN, FUN.VALUE, ...)` applies `FUN` to every element of `X` and returns a vector, *enforcing* that each result matches the type/length given by `FUN.VALUE` (here `character(1)`, i.e. exactly one string per call) — safer than `sapply()`, which can silently return a list if results aren't uniform. Here `` `[[` `` is used *as* `FUN`: for each `x` in `GEO_SOURCES` (each `x` is itself a list like `list(gse="GSE93272", role=..., used_in=...)`), it computes `` `[[`(x, "gse") ``, i.e. `x[["gse"]]`, extracting just the `"gse"` field. So this line extracts `c("GSE93272","GSE110169","GSE15573","GSE89408")` from `GEO_SOURCES`.
**R concept — `lapply()`:** Applies `individual_dataset_entry` to each of those four accession strings, returning a list of four entry-lists (built via closures, as explained above).
**R concept — `c()` on lists:** Concatenating two lists with `c()` produces one flat list containing all elements of both (here: 1 default entry + 4 individual entries = 5 total), as opposed to nesting one list inside the other.
**Data flow:** `GEO_SOURCES` (global.R, static) → 4 accession strings → 4 closures → combined with the 1 hand-built default entry → `PRELOADED_DATASETS`, a 5-element list.
**Example:** `PRELOADED_DATASETS[[1]]$id` → `"__default_merged__"`; `PRELOADED_DATASETS[[2]]$id` → `"GSE93272"`.
**Thesis relevance:** High — this single line is the mechanism that guarantees the dropdown's dataset catalog always mirrors `GEO_SOURCES`; if you add a fifth GEO source to `GEO_SOURCES`, this line alone (no other edit) surfaces it here. Worth citing as evidence of a single-source-of-truth design.

```r
preloaded_choices <- function() {
  setNames(vapply(PRELOADED_DATASETS, `[[`, character(1), "id"),
           vapply(PRELOADED_DATASETS, `[[`, character(1), "label"))
}
```
**What it does:** Builds a named vector where *values* are `id`s (what gets sent to the server as `input$preloaded_choice`) and *names* are `label`s (what the user sees in the dropdown) — this is exactly the `c(name = value)` shape `selectInput(choices=)` expects.
**R concept — `setNames(object, nm)`:** Returns `object` with its `names` attribute set to `nm`. Equivalent to `x <- vapply(...ids...); names(x) <- vapply(...labels...); x`, but as one expression.
**Data flow:** `PRELOADED_DATASETS` → parallel vectors of ids and labels → one named vector.
**What depends on it:** `mod_dataset_ui()`'s `selectInput(ns("preloaded_choice"), choices = preloaded_choices())` (L88), and externally `mod_preprocessing.R:126`.
**Thesis relevance:** Low on its own; it's UI plumbing, but note it in an implementation appendix as the pattern reused for Preprocessing's own dataset picker.

### Block D — `mod_dataset_ui()` (L79-130)

```r
mod_dataset_ui <- function(id) {
  ns <- NS(id)
  tagList(
```
**What:** Defines the UI-building function; `ns <- NS(id)` creates the namespacing function described in §2. `tagList(...)` is Shiny's way of grouping multiple top-level tags into one return value (a plain `list()` of tags wouldn't render correctly without it — `tagList` marks the object as HTML-renderable). **Data flow:** `id` string in → namespaced-ID-generating function `ns` out, closed over for the rest of the function body.

```r
    fluidRow(
      column(6, box(width = NULL, title = "Switch to preloaded data", ...
```
**What:** `fluidRow()`/`column(6, ...)` is Shiny's 12-unit Bootstrap grid — two `column(6, ...)` calls side by side make an even left/right split. `box(...)` is a `shinydashboard`/`bs4Dash`-style panel widget (title bar + body) — confirms this app uses a dashboard UI framework (consistent with `status = "primary"` styling elsewhere in the app). **Thesis relevance:** Purely presentational — mention only if your thesis documents UI/UX design choices, not computational methodology.

```r
selectInput(ns("preloaded_choice"), "Individual dataset",
            choices = preloaded_choices(), selected = character(0), width = "100%"),
```
**What this line does:** Renders a dropdown; `choices` comes from the function explained above; `selected = character(0)` means **no option is pre-selected** — the box starts blank/placeholder rather than defaulting to the first entry.
**Why `selected = character(0)` matters:** It forces an explicit user choice before `input$preloaded_choice` is truthy, so `req(input$preloaded_choice)` downstream (L137, L243) genuinely gates on "user has chosen something," not "user accepted whatever happened to be first in the list."
**R concept:** `character(0)` is a zero-length character vector — a valid, distinct value from `NULL` or `""`, and the idiomatic way to tell `selectInput` "nothing selected yet."
**Thesis relevance:** Worth one sentence — this is a deliberate UX/correctness choice (forces intent) rather than an oversight.

```r
uiOutput(ns("preloaded_note")),
div(style = "display: flex; align-items: center; gap: 10px; flex-wrap: wrap;",
    actionButton(ns("load_preloaded_btn"), "Load this dataset", icon = icon("rotate-left"), class = "btn-primary btn-sm"),
    uiOutput(ns("preloaded_load_message"), inline = TRUE))
```
**What:** A placeholder for the contextual note (see §2), a flex-styled `div` (plain CSS layout, no Shiny logic) holding the "Load this dataset" button and an inline success-message placeholder next to it. **`actionButton`:** Shiny's clickable-button input; its *value* (`input$load_preloaded_btn`) is an integer counter that increments by 1 every click — that counter, not any "boolean," is what `observeEvent()` watches (see Block G). **Thesis relevance:** Low.

```r
        box(
          width = NULL, title = "Fetch from NCBI GEO", status = "primary", solidHeader = FALSE,
          p(strong("GEO Series accession"), " (e.g. ", code("GSE12345"), "). Fetches the expression data and metadata from NCBI. For your own data, use \"Upload your own data\" instead."),
          fluidRow(
            column(4, textInput(ns("geo_accession"), NULL, placeholder = "GSE12345", width = "100%")),
            column(3, actionButton(ns("geo_fetch_btn"), "Fetch", icon = icon("download"), class = "btn-primary btn-sm", width = "100%"))
          ),
          uiOutput(ns("geo_fetch_status")),
          uiOutput(ns("geo_platform_ui")),
          uiOutput(ns("geo_column_mapping")),
          actionButton(ns("geo_load_btn"), "Load this dataset", icon = icon("upload"), class = "btn-primary btn-sm")
        )
```
**What this whole block does:** Builds the second box: a text box for a GEO accession, a Fetch button, and three placeholders that only populate after a successful fetch (status text, an optional platform picker, and the column-mapping dropdowns), plus a final commit button.
**Why `textInput(ns("geo_accession"), NULL, ...)`:** The second argument to `textInput()` is its visible `label`; `NULL` suppresses the label because the preceding `p()` already describes the field — avoids a redundant "GEO Series accession" label appearing twice.
**Why three separate `uiOutput`s rather than one combined one:** Each is independently reactive to a different upstream reactive (`geo_fetch_result()`, `geo_eset()`/`input$geo_platform_choice`, `geo_expr_meta()`) — combining them into one `renderUI` would force Shiny to re-render the whole block (including, e.g., re-drawing the platform selectInput and losing its current selection) any time *any* one of those changed, instead of only the piece that actually needs to change. This is a Shiny performance/UX pattern: keep independently-changing UI in independent outputs.
**Thesis relevance:** The design rationale for splitting `uiOutput`s (avoiding unnecessary re-render / preserving widget state) is a legitimate implementation-detail point if your thesis discusses UI responsiveness.

```r
      column(
        6,
        box(
          width = NULL, title = "Upload your own data", status = "primary", solidHeader = FALSE,
          div(class = "upload-step-label", "STEP 1 · Choose your files"),
          p(strong("Expression matrix"), " — CSV or RDS. Genes in rows, samples in columns; for CSV, the first column is the gene ID."),
          fileInput(ns("expr_file"), "Expression matrix", accept = c(".csv", ".rds", ".Rds")),
          p(strong("Sample metadata"), " — CSV or RDS data frame, one row per sample."),
          fileInput(ns("meta_file"), "Sample metadata", accept = c(".csv", ".rds", ".Rds")),
          uiOutput(ns("upload_preview_ui")),
          tags$hr(),
          div(class = "upload-step-label", "STEP 2 · Map the columns"),
          uiOutput(ns("column_mapping")),
          tags$hr(),
          div(class = "upload-step-label", "STEP 3 · Confirm"),
          actionButton(ns("load_btn"), "Upload Data", icon = icon("upload"), class = "btn-primary btn-sm"),
          div(class = "empty-note", style = "margin-top:6px;", icon("circle-info"),
              "Loads exactly what you provide, as-is - no merging, normalising, or batch correction. For that, use Preprocessing and Batch Correction instead.")
        )
      )
```
**What this does:** The third box, structured as an explicit 3-step wizard (visually, via `div(class="upload-step-label", ...)` labels — pure CSS/markup, no server logic tied to "steps" as a state machine; the "steps" are a UI convention, not a Shiny construct).
**`fileInput(ns("expr_file"), "Expression matrix", accept = c(".csv", ".rds", ".Rds"))`:** Renders a native browser file picker. `accept` is an HTML hint (restricts what the OS file dialog *offers* to select) — **it is not a server-side validation guarantee**; a user can still rename any file to `.csv` and upload it, so the actual format detection happens later in R via `grepl("\\.rds$", ...)` (see Block E) — worth flagging explicitly in your Validation section (§7).
**Data flow (`input$expr_file`):** When a file is selected, Shiny uploads it to a temp path server-side and sets `input$expr_file` to a one-row data frame with columns `name`, `size`, `type`, `datapath` — not the file's contents. Reading the contents is a separate step (`expr_raw()` reactive).
**Thesis relevance:** Medium — the expected file contract ("genes in rows, samples in columns, first CSV column = gene ID") is exactly the kind of assumption you should state explicitly and justify/cite in your Methods (this is the same convention as `DEFAULT_EXPR_RDS` in `global.R`).

```r
    ),
    uiOutput(ns("load_message"))
  )
}
```
**What:** Closes the two-column `fluidRow`, then adds one more placeholder *outside* both columns, spanning the full row width below both boxes — the shared final status message for both the upload and GEO-load paths (but *not* the preloaded path, which has its own inline message next to its own button — an inconsistency worth noting, see §7).

### Block E — `mod_dataset_server()` setup (L132-135)

```r
mod_dataset_server <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
```
**What `moduleServer(id, module)` does:** A Shiny function that (1) creates a namespaced `session` proxy scoped to `id`, and (2) calls your inner `function(input, output, session)` **once per client session**, with `input`/`output` already transparently namespaced — i.e., inside this function, `input$expr_file` really means the browser input tagged `tx_dataset-expr_file`, without you writing `ns()` anywhere on the read side.
**Why this pattern over the older `callModule()`:** `moduleServer()` is the modern (Shiny ≥1.5) idiom; it guarantees the UI and server share exactly one namespace definition path and is what all sibling `mod_*.R` files in this app use.
**`ns <- session$ns`:** Captured for later use inside `renderUI` blocks that construct *new* input IDs (e.g. `selectInput(ns("map_id"), ...)`), since those must be explicitly namespaced — only reads of existing `input$...ids` are auto-namespaced by the surrounding `moduleServer()`, not IDs you're actively creating.
**Data flow:** `id` ("tx_dataset") and `dataset` (the shared reactiveValues) enter as arguments; everything below runs inside this one function's closure, so every reactive/observer defined below shares access to both.

### Block F — `preloaded_note` (L136-148)

```r
    output$preloaded_note <- renderUI({
      req(input$preloaded_choice)
      if (identical(input$preloaded_choice, default_dataset_entry$id)) {
        ...
      } else if (input$preloaded_choice %in% c("GSE93272", "GSE110169")) {
        ...
      } else {
        ...
      }
    })
```
**What construct this is:** `output$name <- renderUI({...})` — a **reactive output binding**. The `{}` block is a reactive expression: Shiny tracks every reactive value read inside it (here, just `input$preloaded_choice`) and automatically re-runs the whole block whenever that value changes, pushing the new returned tag(s) to the browser.
**`req(input$preloaded_choice)`:** Short for "require" — if the argument is falsy/missing (here, before the user picks anything, `input$preloaded_choice` is `character(0)` because of `selected = character(0)` above), `req()` raises a special "silent" condition that stops the reactive *without* an error message shown to the user — the output is simply left blank/unrendered. This is the standard Shiny guard against rendering UI for inputs that don't exist yet.
**`identical(a, b)`:** Strict equality — checks type and value both match exactly (unlike `==`, which vectorizes and coerces types, and can behave surprisingly with `NA`/zero-length vectors). Using `identical()` against `default_dataset_entry$id` (the string sentinel from Block C) rather than `==` is the correct, idiomatic choice here.
**Three-branch `if/else if/else`:** Chooses one of three static informational messages based on which dataset is selected — a plain R conditional, not Shiny-specific, but its *branches* are Shiny tag objects (`p(class=..., icon(...), "...")`).
**What it does scientifically:** Communicates data-provenance caveats to the user — e.g., that GSE93272/GSE110169 in this dropdown are *not* their standalone raw files but samples filtered out of the already-merged/batch-corrected cohort (a fact only knowable from the code — see `merged_training_subset()` in `global.R`), versus the other two sources which really are raw, single-platform, probe-level.
**Thesis relevance:** High — this is a direct, code-verifiable statement of data provenance you should mirror in your Methods/Data section: **GSE93272 and GSE110169, when picked individually from this dropdown, are NOT independent raw reads — they are the subset of the merged, batch-corrected cohort belonging to that GEO series** (confirmed by tracing `load_individual_dataset()` → `merged_training_subset()` in `global.R:745-763`, given the raw per-source `.rds` files for those two are absent from this deployment). This is scientifically important: any per-source QC/analysis run against "GSE93272" from this dropdown is running on already-batch-corrected data, not the original platform-level signal.

### Block G — `meta_raw` and `expr_raw` reactives (L150-177)

```r
    meta_raw <- reactive({
      req(input$meta_file)
      path <- input$meta_file$datapath
      if (grepl("\\.rds$", input$meta_file$name, ignore.case = TRUE)) {
        d <- readRDS(path)
        validate(need(is.data.frame(d), "The uploaded metadata RDS file must contain a data frame."))
        as.data.frame(d)
      } else {
        as.data.frame(data.table::fread(path, showProgress = FALSE))
      }
    })
```
**What construct:** `reactive({...})` creates a **reactive expression** — unlike `renderUI`, it doesn't produce output for the browser; it produces a *cached, lazily-evaluated value* other reactive code can call as a function (`meta_raw()`). Its result is memoized: if called twice without any of its dependencies (`input$meta_file`) changing, the second call returns the cached value instantly without re-running the block.
**`req(input$meta_file)`:** Blocks until a metadata file has actually been chosen.
**`input$meta_file$datapath`:** The server-local temp file path Shiny copied the upload to (not the original filename/path on the user's machine, for security/portability).
**`grepl("\\.rds$", input$meta_file$name, ignore.case = TRUE)`:** Regex match on the *original filename* (not `datapath`, which is a random temp name without a `.rds` extension) to decide the parser branch. `\\.rds$` = literal dot, literal "rds", end of string; `ignore.case = TRUE` also matches `.Rds`/`.RDS`.
**RDS branch:** `readRDS(path)` deserializes an R object exactly as saved; `validate(need(is.data.frame(d), ...))` guards against an RDS file containing something other than a data frame (e.g. a list, a matrix, a vector) — without this check, later code like `meta[[input$map_id]]` would either error opaquely or silently misbehave.
**CSV branch:** `data.table::fread(path, showProgress=FALSE)` — a fast, C-backed CSV reader (chosen over base `read.csv` for speed on potentially large sample sheets); wrapped in `as.data.frame()` because `fread()` natively returns a `data.table`, and downstream code (`colnames()`, `meta[[col]]`, `meta$sample <- ...`) is written assuming plain `data.frame` semantics (a `data.table` has different `[`/`$<-` behavior in some edge cases).
**What it returns:** A `data.frame`, one row per sample, whatever columns the uploaded file had.
**What depends on it:** `upload_preview_ui`, `column_mapping`, and the `load_btn` handler.

```r
    expr_raw <- reactive({
      req(input$expr_file)
      path <- input$expr_file$datapath
      if (grepl("\\.rds$", input$expr_file$name, ignore.case = TRUE)) {
        readRDS(path)
      } else {
        m <- as.data.frame(data.table::fread(path, showProgress = FALSE))
        rn <- as.character(m[[1]])
        m <- as.matrix(m[, -1, drop = FALSE])
        rownames(m) <- rn
        m
      }
    })
```
**What's different from `meta_raw`:** The CSV branch here does real reshaping, not just a type cast, because an expression matrix must end up as a numeric **matrix with gene-ID rownames**, not a data frame with a gene-ID *column* (matrix algebra downstream — PCA, batch correction, WGCNA — needs `matrix`, and rownames double as the gene identifier index every merge/lookup elsewhere in the app keys on).
**`rn <- as.character(m[[1]])`:** Takes the first column (assumed to be the gene ID, per the UI's own instructions: "for CSV, the first column is the gene ID") and keeps it aside as a character vector.
**`m <- as.matrix(m[, -1, drop = FALSE])`:** `m[, -1]` drops column 1 by negative indexing (every column *except* the first); `drop = FALSE` prevents R's default behavior of silently converting a single-remaining-column data frame down to a bare vector (which would lose the sample-name column headers) — this is a classic, important R gotcha: `df[, -1]` on a 2-column data frame returns a vector, not a 1-column data frame, unless `drop = FALSE` is specified. `as.matrix()` then coerces everything to a numeric matrix (assuming all remaining columns are numeric — see §7 for what happens if they're not).
**`rownames(m) <- rn`:** Reattaches the earlier-saved gene IDs as row names.
**RDS branch:** Assumes the saved object is *already* a properly-shaped matrix — no validation at all here (contrast with `meta_raw`'s `validate(need(is.data.frame(d), ...))`) — an asymmetry worth flagging (§7).
**Thesis relevance:** High — this is the exact transformation you should describe formally in Methods when explaining "how is a user-uploaded CSV expression matrix parsed into the internal representation": *first column → row names (gene identifiers); remaining columns → numeric expression matrix, column order preserved as sample order.*

### Block H — `upload_preview_ui` (L179-196)

```r
    output$upload_preview_ui <- renderUI({
      req(input$expr_file, input$meta_file)
      preview <- tryCatch(
        list(expr = expr_raw(), meta = meta_raw()),
        error = function(e) e
      )
      if (inherits(preview, "error")) {
        return(div(class = "empty-note", icon("triangle-exclamation"),
                    paste("Could not read the uploaded file(s):", conditionMessage(preview))))
      }
      div(class = "empty-note", icon("circle-info"),
          sprintf("Read %s: %s features x %s samples. Read %s: %s rows. Map the columns below, then click \"Load dataset\".",
                  input$expr_file$name, format(nrow(preview$expr), big.mark = ","), ncol(preview$expr),
                  input$meta_file$name, nrow(preview$meta)))
    })
```
**What construct:** Another `renderUI`, depending on both `expr_raw()` and `meta_raw()` (and the raw `input$expr_file`/`input$meta_file` for filenames).
**`req(input$expr_file, input$meta_file)`:** `req()` accepts multiple arguments and requires *all* to be truthy — gate: don't attempt a preview until both files are present.
**`tryCatch(expr, error = function(e) e)`:** Base-R structured error handling. Unlike `validate(need(...))` (which is Shiny-only and designed to *replace* output with a message), `tryCatch` here is used to convert any thrown error (e.g. a malformed CSV, an RDS that isn't even readable) into an ordinary R object (the condition `e` itself) that the code can inspect afterward with `inherits(preview, "error")`, rather than letting it propagate and blow up the reactive graph. This is the **sentinel-object error pattern** the file's own header comment (L300-306) names explicitly, and it's reused four more times in this file (`geo_fetch_result`, `load_btn`'s handler, `geo_expr_meta`, `geo_load_btn`'s handler) — a deliberate, consistent idiom, not an accident.
**`inherits(preview, "error")`:** Checks whether the object's class vector includes `"error"` — true for anything `tryCatch`'s `error=` handler returns unmodified (all R conditions raised by `stop()` carry class `c("simpleError","error","condition")`).
**`conditionMessage(preview)`:** Extracts just the human-readable message string from an error/condition object.
**Success branch — `sprintf(...)`:** String templating; `%s` substitutes each following argument as a string, in order. `format(nrow(preview$expr), big.mark = ",")` formats the gene count with thousands separators (e.g. `15,763`) for readability.
**Why this exists as a *separate* preview step, before any "Load" button:** Per the file's own comment (L179-181), so a large upload doesn't look like nothing happened while it's actually just parsing — user feedback for an operation that can take real wall-clock time on a genome-wide matrix.
**Thesis relevance:** Medium — worth one sentence acknowledging this as a UX/robustness feature (early feedback + graceful degradation on malformed files), not a scientific computation itself.

### Block I — `guess_col()` (L198-211)

```r
    guess_col <- function(cols, exact, contains = exact, fallback = cols[1]) {
      hit <- cols[tolower(cols) %in% tolower(exact)]
      if (length(hit) > 0) return(hit[1])
      hit <- cols[grepl(paste(contains, collapse = "|"), cols, ignore.case = TRUE)]
      if (length(hit) > 0) return(hit[1])
      fallback
    }
```
**What construct:** An ordinary (non-reactive) function *defined inside* `moduleServer()`'s body — meaning a fresh copy is created per user session (harmless here since it's pure/stateless, but worth knowing: unlike the top-level constants in Block C, this is **not** shared across sessions/instantiations).
**Parameters:** `cols` — the actual column names present in the uploaded/fetched metadata; `exact` — a vector of candidate names to match verbatim (case-insensitively); `contains` — defaults to the same vector as `exact`, used for substring matching if no exact hit; `fallback` — defaults to the first column in `cols` if nothing else matches.
**Line 1 (`hit <- cols[tolower(cols) %in% tolower(exact)]`):** Case-insensitive **exact** match — lowercases both sides, then subsets `cols` to just the entries that appear in the lowercased `exact` list.
**`if (length(hit) > 0) return(hit[1])`:** If one or more exact matches exist, take the first (in case of duplicate-content column names, an edge case) and return immediately — an **early return**, a common R pattern for guard-clause-style functions.
**Line 2 (`grepl(paste(contains, collapse="|"), cols, ignore.case=TRUE)`):** Falls back to **substring** matching: `paste(contains, collapse="|")` turns e.g. `c("sample","sample_id","id")` into the regex `"sample|sample_id|id"`; `grepl(..., cols)` then tests each column name for whether it *contains* any of those substrings anywhere.
**Final line (`fallback`):** If neither exact nor substring matching found anything, return the fallback value (by default, just the first column — i.e., "guess something rather than nothing," but see §7 for why this default is scientifically risky).
**Why this exists at all (the comment at L198-204 explains it directly):** Without it, a plain `selectInput()` would default to whichever column happens to come first in the uploaded file — often a literal sample-ID column — silently producing a "group" variable with one unique value per sample (every downstream group-comparison would then be meaningless) instead of erroring loudly.
**Data flow:** `colnames(meta_raw())` (or `colnames(em$meta)` for GEO) → `guess_col()` → a single string → used as `selected=` in a `selectInput`.
**Example:** `guess_col(c("Sample_ID","Diagnosis","Sex"), c("sample","sample_id","id","geo_accession","accession"))` → exact match check: none of the candidate names equal any column name lowercased... actually `"sample_id"` lowercased is `"sample_id"`, and `tolower(cols)` includes `"sample_id"` → matches → returns `"Sample_ID"`.
**Thesis relevance:** High — this is exactly the kind of implementation choice a thesis defense committee might probe ("how do you guarantee correct column mapping on arbitrary user data?") — the honest, code-grounded answer is: *heuristic name-matching with a human-in-the-loop override* (the dropdown remains editable — this only sets the default), not fully automatic/guaranteed-correct inference. State this explicitly as a limitation (§7/§9).

### Block J — `column_mapping` renderUI (L213-230)

```r
    output$column_mapping <- renderUI({
      req(input$meta_file)
      cols <- colnames(meta_raw())
      tagList(
        selectInput(ns("map_id"), "Sample ID column", choices = cols,
                    selected = guess_col(cols, c("sample", "sample_id", "id", "geo_accession", "accession")),
                    selectize = FALSE),
        selectInput(ns("map_group"), "Group / diagnosis column", choices = cols,
                    selected = guess_col(cols, c("group", "diagnosis", "disease", "condition", "status", "phenotype")),
                    selectize = FALSE),
        selectInput(ns("map_sex"), "Sex column (optional)", choices = c("(none)", cols),
                    selected = guess_col(cols, c("sex", "gender"), fallback = "(none)"),
                    selectize = FALSE),
        selectInput(ns("map_batch"), "Batch column (optional)", choices = c("(none)", cols),
                    selected = guess_col(cols, c("batch", "cohort", "platform", "dataset"), fallback = "(none)"),
                    selectize = FALSE)
      )
    })
```
**What it does:** Dynamically builds four dropdowns whose `choices` are the *actual* uploaded metadata's column names, each pre-selected via `guess_col()`.
**Design detail:** `map_id`/`map_group` have `fallback = cols[1]` (the default) — i.e., these are *required* fields and will always default to *some* column, forcing the user to notice and correct a wrong guess rather than leaving the field empty. `map_sex`/`map_batch` instead pass `fallback = "(none)"` explicitly — these are *optional* fields, so guessing nothing is the safer wrong answer (an absent sex/batch column is a normal, valid state; an absent sample-ID or group column is not).
**`selectize = FALSE`:** Uses the plain HTML `<select>` widget instead of Shiny's default `selectize.js`-enhanced searchable dropdown — likely a deliberate simplicity/performance choice for a short, always-visible column list (searchability matters more for long lists).
**Data flow:** `meta_raw()` (a `data.frame`) → `colnames()` → `guess_col()` (four times, once per field) → four new reactive inputs (`input$map_id`, `input$map_group`, `input$map_sex`, `input$map_batch`) come into existence in the session the moment this UI renders.
**What depends on it:** The `observe()` enabling `load_btn` (Block K) and the `load_btn` `observeEvent` handler (Block L) both read `input$map_id`/`input$map_group`/etc. — but note these inputs **do not exist** until this `renderUI` has actually run at least once (i.e., until a metadata file has been uploaded) — this is why `req()`/null-checks are needed downstream rather than assuming these inputs are always present.

### Block K — enable/disable `load_btn` (L232-240)

```r
    observe({
      ready <- !is.null(input$expr_file) && !is.null(input$meta_file) &&
        !is.null(input$map_id) && !is.null(input$map_group)
      if (isTRUE(ready)) shinyjs::enable("load_btn") else shinyjs::disable("load_btn")
    })
```
**What construct:** `observe({...})` — like `reactive()`, it re-runs whenever any reactive value it reads changes, but unlike `reactive()`, it has **no return value that anything else can consume** — it exists purely for its *side effect* (here, enabling/disabling a button via `shinyjs`). This is the standard Shiny distinction: use `reactive()`/`renderUI()` when you need a *value*; use `observe()`/`observeEvent()` when you need an *action*.
**`!is.null(x)`:** Before any file is chosen, `input$expr_file` is `NULL` (not `character(0)` — different from the `preloaded_choice` dropdown's empty state); this checks presence.
**`&&`:** Scalar-only logical AND (as opposed to `&`, vectorized) — appropriate here since every operand is a single `TRUE`/`FALSE`; using `&&` also short-circuits (stops evaluating as soon as one operand is `FALSE`), a minor efficiency/robustness detail.
**`isTRUE(ready)`:** Defensive — guarantees the branch only fires on a literal scalar `TRUE`, not on `NA`, a zero-length logical, or anything else `&&` chain edge cases could theoretically produce.
**`shinyjs::enable("load_btn")` / `disable("load_btn")`:** Functions from the `shinyjs` package that toggle a DOM element's `disabled` attribute via injected JavaScript — note the ID passed is the **unnamespaced** string `"load_btn"`; `shinyjs` functions called from inside a `moduleServer()` are automatically namespace-aware in recent `shinyjs` versions (they resolve `"load_btn"` to `"tx_dataset-load_btn"` using the enclosing module's namespace), so this is correct as written, not a bug — but it's worth verifying against the installed `shinyjs` version (see §7/§11).
**Why this exists:** Pure UX guardrail — makes it visually obvious, without clicking, whether the form is ready, and prevents a click from firing `observeEvent(input$load_btn)` with missing required inputs (though that handler's own `req()` at L258 would also block it — this is a belt-and-suspenders pair, not the only line of defense).
**Thesis relevance:** Low scientifically; a defensible UI/robustness design point.

### Block L — `observeEvent(input$load_preloaded_btn, ...)` (L242-255)

```r
    observeEvent(input$load_preloaded_btn, {
      req(input$preloaded_choice)
      entry <- Find(function(d) d$id == input$preloaded_choice, PRELOADED_DATASETS)
      req(entry)
      d <- entry$load()
      dataset$staged_expr <- d$expr
      dataset$staged_meta <- d$meta
      dataset$staged_source <- d$source
      output$preloaded_load_message <- renderUI(
        span(style = "color: var(--color-success); font-size: 13px; font-weight: 600;", icon("check"), " ",
             sprintf("Previewed %s: %s genes across %s samples. This doesn't change what any sub-module runs on - go to Preprocessing and pick \"Currently loaded dataset\" to analyze it.",
                      entry$label, format(nrow(d$expr), big.mark = ","), ncol(d$expr)))
      )
    })
```
**What construct:** `observeEvent(eventExpr, handlerExpr)` — unlike plain `observe()`, this only re-runs when `eventExpr` (here, `input$load_preloaded_btn`, the click counter) *changes*, ignoring changes to any other reactive value read inside `handlerExpr` (e.g. `input$preloaded_choice` can change freely without re-triggering this block — only a fresh click does). This is precisely why a "Load"/"commit" action is written as `observeEvent` on the button, not as `observe()`/`reactive()` on the underlying data — those would fire the moment the user *touches* the dropdown, before they've decided to commit.
**`Find(f, list)`:** Base-R function; returns the *first* element of `list` for which `f(element)` is `TRUE`, or `NULL` if none match. Here it linearly searches `PRELOADED_DATASETS` for the entry whose `$id` matches the dropdown's current value.
**`req(entry)`:** Defensive guard — if `Find()` somehow returned `NULL` (shouldn't normally happen, since `preloaded_choice`'s choices are generated *from* `PRELOADED_DATASETS` itself, so every valid selection has a matching entry — but this guards against any future desync), stop here rather than erroring on `entry$load()`.
**`d <- entry$load()`:** Invokes the closure built back in Block C — this is the moment (and the *only* moment) the actual file read happens for the preloaded path, i.e., lazy evaluation realized.
**The three assignment lines — the core mutation:**
```r
dataset$staged_expr <- d$expr
dataset$staged_meta <- d$meta
dataset$staged_source <- d$source
```
**What happens when R executes this:** `dataset` is a `reactiveValues` object (created in `server.R:7-10`). Assigning `dataset$staged_expr <- d$expr` does two things: (1) stores `d$expr` under the key `"staged_expr"` inside `dataset` (note: `reactiveValues` objects can have new fields created dynamically by assignment — they don't need to be pre-declared, unlike a strict S4/R6 object); (2) **invalidates** every reactive consumer that has ever read `dataset$staged_expr` (right now, none has yet in this session — but `mod_preprocessing.R`'s `pp_preloaded_read()` will, the next time it runs), scheduling them to re-run.
**Why this matters architecturally:** Because `dataset` was created *once* in `server.R` and passed *by reference* into `mod_dataset_server(id, dataset)`, this write is visible to every other module that was also handed the same `dataset` object — cross-module communication via a shared reactive value, Shiny's standard pattern for module-to-module state sharing (as opposed to, say, global variables, which wouldn't be reactive/session-safe).
**`output$preloaded_load_message <- renderUI(...)`:** Note this is **not** wrapped in `{}` here — `renderUI(expr)` accepts either a `{}` block or a bare expression; since the whole thing is one `sprintf`/`span` call, no braces are needed. Sets the inline success message next to the "Load this dataset" button.
**Data flow summary for this block:** button click (`input$load_preloaded_btn` increments) → `input$preloaded_choice` (string) → `PRELOADED_DATASETS` lookup → `entry$load()` → `list(expr=, meta=, source=)` → three fields on shared `dataset` reactiveValues → downstream: `mod_preprocessing.R` reactively re-reads them next time it runs; also a rendered success message.
**Thesis relevance:** Very high — this is the concrete implementation of the "staging, not activation" design principle stated in the file header; you should walk through exactly this block in your Implementation chapter as the canonical example of that pattern, and note explicitly that **no validation of the loaded data's biological sanity happens here** (no check for, e.g., all-NA expression, degenerate group column, etc. — see §7).

### Block M — `observeEvent(input$load_btn, ...)` — the upload commit handler (L257-298)

```r
    observeEvent(input$load_btn, {
      req(input$expr_file, input$meta_file, input$map_id, input$map_group)

      result <- tryCatch({
        expr <- expr_raw()
        meta <- meta_raw()
        meta$sample <- as.character(meta[[input$map_id]])
        meta$group  <- as.character(meta[[input$map_group]])
        meta$sex    <- if (!identical(input$map_sex, "(none)")) as.character(meta[[input$map_sex]]) else NA_character_
        meta$batch  <- if (!identical(input$map_batch, "(none)")) as.character(meta[[input$map_batch]]) else NA_character_

        common <- intersect(colnames(expr), meta$sample)
        validate(need(
          length(common) >= 4,
          "Fewer than 4 sample IDs in the expression matrix match the metadata sample-ID column. Check the column mapping."
        ))

        expr <- expr[, common, drop = FALSE]
        meta <- meta[match(common, meta$sample), , drop = FALSE]
        list(expr = expr, meta = meta)
      }, error = function(e) e)
```
**`req(...)` with 4 arguments:** Requires all four to be present before proceeding — a second gate beyond the button-enabled state (defense in depth against, e.g., programmatic input manipulation or a race where the button briefly re-enables).
**`meta$sample <- as.character(meta[[input$map_id]])`:** Creates/overwrites a column literally named `sample` in the metadata data frame, populated from whichever column the user mapped as the ID column. `meta[[input$map_id]]` — double-bracket extraction by a *variable* column name (the string stored in `input$map_id`, e.g. `"Sample_ID"`) — returns that column as a vector; `as.character()` coerces it to string (defensive: guards against the ID column being read as a factor or numeric type, which would otherwise break later string-based operations like `intersect()`).
**Same pattern for `group`/`sex`/`batch`:** Each is harmonized into a fixed, standardized column name (`group`, `sex`, `batch`) regardless of what the source file called it — this is the **column-name harmonization** step that lets every downstream sub-module assume fixed column names rather than re-discovering them per dataset.
**`if (!identical(input$map_sex, "(none)")) ... else NA_character_`:** For the two *optional* fields, if the user left the mapping at the placeholder `"(none)"`, the harmonized column is filled with `NA_character_` (a typed, length-1 `NA` for a character vector — recycled to the full column length by R's standard vector-recycling rule when assigned into a data frame column) rather than being omitted entirely — meaning `meta$sex`/`meta$batch` **always exist** as columns after this block, just possibly all-`NA`. This matters: downstream code can safely assume `"sex" %in% colnames(meta)` is always `TRUE` for a dataset loaded through this path.
**`common <- intersect(colnames(expr), meta$sample)`:** Base-R set intersection — the sample IDs that appear in *both* the expression matrix's column names and the metadata's harmonized `sample` column. This is the critical **sample-matching** step: it does not assume the two files list samples in the same order, or even the same *set* — it computes the overlap.
**`validate(need(length(common) >= 4, "..."))`:** A scientifically meaningful minimum-sample-size gate — refuses to proceed with fewer than 4 matched samples. (Note: 4 is a hardcoded constant here and again at L270/L424 for the GEO path — worth checking why 4 specifically was chosen; likely a pragmatic "too few samples for any downstream statistic to be meaningful" floor rather than a value derived from a specific statistical power calculation — **cannot determine the exact rationale from this file alone**; check `mod_preprocessing.R`/other modules or ask the code's original author/your own notes if this constant needs formal justification for the thesis.)
**`expr <- expr[, common, drop = FALSE]`:** Subsets the expression matrix to only the matched samples, in the order `common` lists them (which, per `intersect()`'s documented behavior, follows the order of its *first* argument, i.e., `colnames(expr)`'s original order, filtered).
**`meta <- meta[match(common, meta$sample), , drop = FALSE]`:** **This is the line that actually aligns row order to `expr`'s column order.** `match(common, meta$sample)` returns, for each entry in `common`, the row-index in `meta$sample` where it's found — so `meta[that index vector, ]` reorders (and subsets) `meta`'s rows to be in exactly the same sample order as `expr`'s columns now are. Without this line, `expr` and `meta` could both contain the right *set* of samples but in different, mismatched *order* — a silent, dangerous correctness bug if omitted (row *i* of `meta` would describe a different sample than column *i* of `expr`). This is one of the most important lines in the whole file from a correctness standpoint.
**`list(expr = expr, meta = meta)`:** The success value, if nothing threw.
**`error = function(e) e`:** Same sentinel-object pattern as before — turns any error (e.g. `meta[[input$map_id]]` being `NULL` if the mapping references a nonexistent column, or `validate()`'s own thrown condition) into a plain returned object.

```r
      if (inherits(result, "error")) {
        output$load_message <- renderUI(
          div(class = "empty-note", icon("triangle-exclamation"), paste("Could not load this dataset:", conditionMessage(result)))
        )
      } else {
        dataset$staged_expr <- result$expr
        dataset$staged_meta <- result$meta
        dataset$staged_source <- paste0("Uploaded dataset: ", input$expr_file$name, " + ", input$meta_file$name)
        n_dup <- sum(duplicated(rownames(result$expr)))
        output$load_message <- renderUI(
          tagList(
            div(class = "empty-note", icon("check"),
                sprintf("Previewed %s genes across %s samples. This doesn't change what any sub-module runs on - go to Preprocessing and pick \"Currently loaded dataset\" to analyze it.",
                        format(nrow(result$expr), big.mark = ","), ncol(result$expr))),
            if (n_dup > 0) div(class = "empty-note", icon("triangle-exclamation"),
                sprintf("%d duplicated feature identifier(s) were detected in this dataset. All rows are kept here, but downstream row-name-keyed steps (e.g. the Preprocessing merge tab) will keep only the first occurrence of each - rename duplicates in your source file if this is unintended.", n_dup))
          )
        )
      }
    })
```
**Branch logic:** Mirrors the pattern used throughout the file: `tryCatch` result inspected via `inherits(x, "error")`, error branch renders a red warning, success branch commits the staged fields and renders a green confirmation.
**`n_dup <- sum(duplicated(rownames(result$expr)))`:** `duplicated(x)` returns a logical vector, `TRUE` at every position that repeats an earlier value; `sum()` on a logical vector counts the `TRUE`s (R coerces `TRUE`/`FALSE` to `1`/`0` under arithmetic). This *counts* duplicate gene/probe identifiers in the (now sample-filtered) expression matrix — a genuinely important QC signal, since duplicate row names later break any row-name-keyed join/merge (explicitly flagged in the message itself, referencing Preprocessing's merge tab). Note: this check runs, and the warning displays, but **the data is still staged as-is** — duplicates are not removed, deduplicated, or blocked here; this module only *warns*.
**`if (n_dup > 0) div(...)`:** An `if` without an `else`, used as an *expression* inside `tagList()` — when the condition is `FALSE`, `if` with no `else` returns `NULL` **invisibly**, and `tagList()` (like most Shiny tag-list functions) silently drops `NULL` elements when rendering — a common, idiomatic R/Shiny pattern for "conditionally include this tag."
**Data flow summary:** raw uploaded files → parsed (`expr_raw`/`meta_raw`) → harmonized metadata columns → sample-matched/reordered `expr`+`meta` → staged onto shared `dataset` → user-facing confirmation (with a duplicate-ID advisory if relevant).
**Thesis relevance:** Very high — this is the block to cite for "how are user-uploaded gene identifiers validated," and you should explicitly document, as a **known limitation**, that duplicate feature IDs are detected but not resolved at this stage (resolution is deferred to Preprocessing's merge step, which the message itself states keeps only the first occurrence — worth verifying that claim against `mod_preprocessing.R`'s actual merge code for full accuracy, rather than trusting the UI string alone).

### Block N — GEO fetch (`geo_fetch_result`) (L300-322)

```r
    geo_fetch_result <- eventReactive(input$geo_fetch_btn, {
      if (!requireNamespace("GEOquery", quietly = TRUE)) {
        return(simpleError("The GEOquery package is not installed in this deployment. Install it with BiocManager::install(\"GEOquery\") to enable fetching by GEO accession, or use \"Upload your own data\" instead."))
      }
      acc <- toupper(trimws(input$geo_accession %||% ""))
      if (!grepl("^GSE[0-9]+$", acc)) {
        return(simpleError("Enter a valid GEO Series accession, e.g. GSE12345."))
      }
      tryCatch({
        gset <- suppressMessages(GEOquery::getGEO(acc, GSEMatrix = TRUE))
        if (!is.list(gset) || length(gset) == 0) {
          stop(paste(acc, "returned no series matrix from GEO - check the accession is a Series (GSExxxxx), not a Sample (GSM) or Platform (GPL) ID."))
        }
        list(acc = acc, platforms = gset)
      }, error = function(e) e)
    })
```
**What construct:** `eventReactive(eventExpr, valueExpr)` — the "value-returning" counterpart to `observeEvent`: like a normal `reactive()`, its result (`geo_fetch_result()`) can be called by other reactive code, **but** it only *re-evaluates* `valueExpr` when `eventExpr` (`input$geo_fetch_btn`, the click counter) changes — reading `input$geo_accession` inside the body does *not* by itself trigger a re-fetch, exactly analogous to why `observeEvent` was used for the "commit" handlers above. This is the correct construct for "fetch on click, but expose the result as a reusable reactive value" (as opposed to `observeEvent`, which has no return value, or plain `reactive()`, which would re-fetch on every keystroke in the accession box).
**`requireNamespace("GEOquery", quietly = TRUE)`:** Checks whether the `GEOquery` package is installed **without** loading/attaching it (unlike `library()`) and without printing anything (`quietly=TRUE`) — the idiomatic way to make a package an *optional* runtime dependency: the app can be deployed without `GEOquery` installed, and this feature degrades gracefully to an explanatory error message instead of the whole app failing to start (which is what would happen if `library(GEOquery)` were called unconditionally at the top of a sourced file and the package were missing).
**`simpleError(msg)`:** Base R's constructor for a plain condition object of class `c("simpleError","error","condition")` — used here to construct an "error" object *directly*, without actually calling `stop()`/`tryCatch()`, since there's nothing to catch yet (this is a validation failure, not a caught runtime exception) — but it's returned in the same sentinel-object shape as everything else in this file, so downstream code can treat it identically via `inherits(x, "error")`.
**`acc <- toupper(trimws(input$geo_accession %||% ""))`:** `trimws()` strips leading/trailing whitespace; `toupper()` normalizes case; `input$geo_accession %||% ""` — the `%||%` operator (defined in `global.R:817`: `` `%||%` <- function(a, b) if (is.null(a)) b else a ``) substitutes `""` if the input is `NULL` (e.g., before the text box has ever been touched, though in practice `textInput` initializes to `""`, not `NULL`, so this is a defensive no-op in the common case, not dead code — it protects against any edge case where the input genuinely is `NULL`).
**`grepl("^GSE[0-9]+$", acc)`:** Validates the accession format strictly: must start with literal `GSE`, followed by one or more digits, and nothing else (`^`/`$` anchor to the whole string) — rejects `GSM...`/`GPL...` IDs and malformed input *before* ever hitting the network, which is both a UX improvement (instant feedback) and good API-citizenship (don't send obviously-invalid requests to NCBI).
**`GEOquery::getGEO(acc, GSEMatrix = TRUE)`:** The actual network call — downloads and parses the GEO Series Matrix file(s) for the accession, returning (per GEOquery's documented behavior) a **list of `ExpressionSet` objects**, one per platform the series used (most series use exactly one platform, hence a length-1 list; multi-platform series, e.g. a mixed Illumina+Affymetrix study, return more than one).
**`suppressMessages(...)`:** `GEOquery::getGEO()` is known to print verbose progress/status messages; this suppresses them from cluttering the R console/logs without suppressing warnings or errors.
**`if (!is.list(gset) || length(gset) == 0) stop(...)`:** Defensive check — even a syntactically valid `GSExxxxx` accession can return nothing useful (e.g., a superseries with no matrix file, or a withdrawn record); this converts that silent-empty-result case into an explicit, informative `stop()`, which the enclosing `tryCatch` then captures as the returned error object.
**`list(acc = acc, platforms = gset)`:** Success value — the normalized accession string plus the raw list of `ExpressionSet`s (still possibly length > 1 at this point).
**Thesis relevance:** High for a Data Acquisition/Reproducibility subsection — document precisely which GEOquery function and parameter (`getGEO(acc, GSEMatrix = TRUE)`) is used, since that determines exactly what data GEO gives you (series-matrix-derived processed expression values, **not** raw platform files/FASTQs) — an important scientific caveat: **live-fetched GEO data via this path is whatever normalized/processed matrix the original submitter deposited in the Series Matrix file**, which may differ in normalization method from dataset to dataset — worth an explicit limitation statement.

### Block O — platform selection (L324-341)

```r
    output$geo_platform_ui <- renderUI({
      res <- geo_fetch_result()
      req(res); req(!inherits(res, "error"))
      if (length(res$platforms) <= 1) return(NULL)
      selectInput(ns("geo_platform_choice"), "This series spans multiple platforms - pick one",
                  choices = names(res$platforms), width = "100%")
    })

    geo_eset <- reactive({
      res <- geo_fetch_result()
      req(res); req(!inherits(res, "error"))
      if (length(res$platforms) > 1) {
        req(input$geo_platform_choice)
        res$platforms[[input$geo_platform_choice]]
      } else {
        res$platforms[[1]]
      }
    })
```
**`req(res); req(!inherits(res, "error"))`:** Two separate `req()` calls on separate lines — first ensures a fetch has actually happened (`res` isn't `NULL`), second ensures it wasn't an error result; written as two statements rather than one `req(res, !inherits(res,"error"))` — functionally near-equivalent (both must be truthy to proceed) but slightly clearer to read stepwise. This pattern (`req(res); req(!inherits(res,"error"))`) recurs multiple times in this file — worth naming as a **recognized idiom**: "block until we have a valid (non-error) reactive value."
**`if (length(res$platforms) <= 1) return(NULL)`:** If the fetched series is single-platform, this `renderUI` returns `NULL` — Shiny interprets a `NULL` `renderUI` result as "render nothing," so the placeholder `<div>` simply stays empty; no platform picker appears.
**`selectInput(ns("geo_platform_choice"), ..., choices = names(res$platforms))`:** `GEOquery::getGEO()`'s returned list is named by platform accession (e.g. `"GSE12345_GPL570"`) — `names(res$platforms)` surfaces those as the picker's choices.
**`geo_eset` reactive:** Resolves to a single `ExpressionSet` — either the one the user picked (if multi-platform, gated by `req(input$geo_platform_choice)` so this doesn't proceed until a choice is made) or the sole platform (if single-platform, no user action needed). This is the reactive every downstream GEO-path step (`geo_expr_meta`) builds on.
**Data flow:** `geo_fetch_result()` (list of ExpressionSets) → (optional user platform choice) → `geo_eset()` (single ExpressionSet).

### Block P — `geo_expr_meta` (L343-362)

```r
    geo_expr_meta <- reactive({
      tryCatch({
        eset <- geo_eset()
        ex <- Biobase::exprs(eset)
        if (nrow(ex) == 0 || ncol(ex) == 0) {
          stop("This GEO series has no expression matrix in its series matrix file - common for RNA-seq series that only deposit raw counts as supplementary files. Download that file from the GEO page and use \"Upload your own data\" instead.")
        }
        collapsed <- tryCatch(collapse_probes_to_genes(eset), error = function(e) NULL)
        used_collapse <- !is.null(collapsed) && nrow(collapsed) > 0 && nrow(collapsed) < nrow(ex)
        expr <- if (used_collapse) collapsed else ex
        list(expr = expr, meta = as.data.frame(Biobase::pData(eset)),
             collapsed = used_collapse, platform = Biobase::annotation(eset))
      }, error = function(e) e)
    })
```
**`Biobase::exprs(eset)`:** Bioconductor's `Biobase` package accessor — extracts the numeric expression matrix (features × samples) from an `ExpressionSet` S4 object.
**Empty-matrix check:** Some GEO Series (notably RNA-seq studies that deposit only raw FASTQ/count files as *supplementary* files rather than in the series matrix itself) legitimately have a zero-row or zero-column `exprs()` matrix — this converts that case into an explicit, actionable error message rather than letting downstream code fail confusingly on an empty matrix.
**`collapsed <- tryCatch(collapse_probes_to_genes(eset), error = function(e) NULL)`:** A *nested* `tryCatch`, inside the outer one — calls `global.R`'s `collapse_probes_to_genes()` (which uses `WGCNA::collapseRows()` with the `"MaxMean"` rule to reduce a probe-level matrix to one row per gene symbol, matching the app's own bundled-source preprocessing script per its own comment in `global.R:702-710`), but if that itself errors (e.g., an unrecognized platform annotation schema), the error is swallowed to `NULL` rather than failing the whole reactive — because the collapse step is an *enhancement*, not a hard requirement (probe-level data is still usable).
**`used_collapse <- !is.null(collapsed) && nrow(collapsed) > 0 && nrow(collapsed) < nrow(ex)`:** Three conditions must all hold to actually use the collapsed result: it succeeded (`!is.null`), it's non-empty (`nrow > 0`), and it actually *reduced* row count relative to the raw matrix (`nrow(collapsed) < nrow(ex)`) — this last check is a sanity/no-op guard: `collapse_probes_to_genes()` itself falls back to returning the raw `exprs(eset)` matrix unchanged when it can't find a gene-symbol annotation column (see `global.R:715`, `if (is.na(col)) return(ex)`), so without this third condition, `used_collapse` could be `TRUE` even when no actual collapsing happened, misleadingly labeling probe-level data as "collapsed."
**`Biobase::pData(eset)`:** Extracts the sample-level phenotype/metadata data frame from the `ExpressionSet`.
**`Biobase::annotation(eset)`:** Extracts the platform accession string (e.g. `"GPL570"`).
**Return value:** `list(expr=, meta=, collapsed=<logical>, platform=<string>)` — note this is a *richer* shape than the other loaders' `list(expr=, meta=, source=)` — `geo_load_btn`'s handler (Block R) is responsible for reshaping this into the standard staged shape.
**Thesis relevance:** High — the collapse-or-not logic, and specifically the `used_collapse` triple condition, is worth documenting precisely, since it directly determines whether a live-fetched dataset ends up gene-level (comparable to the app's other sources) or probe-level (not directly comparable without further processing) — and the *user is told which happened* (`geo_fetch_status`, Block Q) rather than it happening silently.

### Block Q — `geo_fetch_status` (L364-380)

```r
    output$geo_fetch_status <- renderUI({
      req(input$geo_fetch_btn)
      res <- geo_fetch_result()
      if (inherits(res, "error")) {
        return(div(class = "empty-note", icon("triangle-exclamation"), paste("Could not fetch from GEO:", conditionMessage(res))))
      }
      em <- geo_expr_meta()
      if (inherits(em, "error")) {
        return(div(class = "empty-note", icon("triangle-exclamation"), conditionMessage(em)))
      }
      div(class = "empty-note", icon("circle-check"),
          sprintf("Fetched %s (platform %s): %s %s across %s samples.%s",
                  res$acc, em$platform, format(nrow(em$expr), big.mark = ","),
                  if (em$collapsed) "genes (collapsed from probes)" else "probes/features",
                  ncol(em$expr),
                  if (!em$collapsed) " No gene-symbol annotation found for this platform - left at probe/feature-ID level. You can still load it as-is, or use Preprocessing's own probe-collapse step afterward." else ""))
    })
```
**`req(input$geo_fetch_btn)`:** Gates on the button having been clicked at least once (its counter is `>= 1`) — before any click, this counter is `0`, which `req()` treats as falsy... actually note: `req()`'s default behavior treats `0` as **falsy** only for logical/`NULL`/empty values by its `cancelOutput`/cast rules — for an integer `0`, `req()` does treat it as falsy (its internal check is essentially `isTruthy()`, and `isTruthy(0)` is `FALSE` for numeric `0`). So this correctly blocks before the first click.
**Nested error handling — two sequential checks:** First checks whether the *fetch itself* failed (`geo_fetch_result()`), then, only if that succeeded, checks whether *extraction* (`geo_expr_meta()`) failed — two independently-failable steps get two independently-worded error messages, which is more informative than a single generic catch-all.
**Success message:** Uses two inline conditional (`if`/`else` as expression) substitutions inside one `sprintf()` call — `if (em$collapsed) "genes (collapsed from probes)" else "probes/features"` supplies one `%s` slot; the final `%s` slot is either an empty string or a full caveat sentence about missing gene-symbol annotation — demonstrating `if`/`else` used purely as a value-producing expression (not a control-flow statement) inside a function call argument list, a very idiomatic R pattern.
**Thesis relevance:** Medium — good evidence of transparent user communication of a scientifically consequential fact (probe- vs gene-level), but not itself a computational step.

### Block R — `geo_column_mapping`, enable/disable, and `geo_load_btn` (L382-448)

```r
    output$geo_column_mapping <- renderUI({
      em <- geo_expr_meta()
      req(em); req(!inherits(em, "error"))
      cols <- colnames(em$meta)
      tagList(
        selectInput(ns("geo_map_group"), "Group / diagnosis column", choices = cols,
                    selected = guess_col(cols, c("group", "diagnosis", "disease", "condition", "status", "phenotype", "characteristics_ch1")),
                    selectize = FALSE),
        ...
      )
    })
```
**Difference from the upload path's `column_mapping`:** No `map_id` dropdown here at all — the file's own comment (L382-385) explains why: `pData()`'s row names *are* the GSM sample accessions, and `exprs()` is already column-indexed by the same GSM IDs, so sample-ID mapping is structurally guaranteed correct, with nothing for the user to specify. Also note the group-column guess list includes `"characteristics_ch1"` — a GEO-specific column-naming convention (`Biobase::pData()` on a GEO-derived `ExpressionSet` often has columns literally named `characteristics_ch1`, `characteristics_ch1.1`, etc., holding free-text sample characteristics) not present in the generic upload-path guess list — a dataset-source-specific tuning of the same generic `guess_col()` helper.

```r
    observe({
      em <- tryCatch(geo_expr_meta(), error = function(e) e)
      ready <- !is.null(em) && !inherits(em, "error") && !is.null(input$geo_map_group)
      if (isTRUE(ready)) shinyjs::enable("geo_load_btn") else shinyjs::disable("geo_load_btn")
    })
```
**Difference from Block K's version:** Here `geo_expr_meta()` is called inside a `tryCatch` *within the `observe()` itself* — `geo_expr_meta()` is itself already wrapped in a `tryCatch` internally (Block P) and never throws, always returning either a valid list or an error-condition object — so this outer `tryCatch` is arguably redundant/defensive-in-excess (belt-and-suspenders again, or possibly a leftover from before `geo_expr_meta()` had its own internal `tryCatch` — **cannot determine intent from this file alone**; worth checking git history if you want the definitive reason for the thesis, otherwise just note it as harmless redundancy in your validation write-up).

```r
    observeEvent(input$geo_load_btn, {
      req(input$geo_map_group)

      result <- tryCatch({
        em <- geo_expr_meta()
        validate(need(!inherits(em, "error"), "No GEO data fetched yet."))
        expr <- em$expr
        meta <- em$meta
        meta$sample <- rownames(meta)
        meta$group  <- as.character(meta[[input$geo_map_group]])
        meta$sex    <- if (!identical(input$geo_map_sex, "(none)")) as.character(meta[[input$geo_map_sex]]) else NA_character_
        meta$batch  <- if (!identical(input$geo_map_batch, "(none)")) as.character(meta[[input$geo_map_batch]]) else NA_character_

        common <- intersect(colnames(expr), meta$sample)
        validate(need(
          length(common) >= 4,
          "Fewer than 4 sample IDs matched between the fetched expression matrix and metadata."
        ))
        expr <- expr[, common, drop = FALSE]
        meta <- meta[match(common, meta$sample), , drop = FALSE]
        label <- sprintf("%s (%s, %s)", geo_fetch_result()$acc, em$platform,
                          if (em$collapsed) "collapsed to genes" else "probe-level, raw")
        list(expr = expr, meta = meta, label = label)
      }, error = function(e) e)

      if (inherits(result, "error")) {
        output$load_message <- renderUI(
          div(class = "empty-note", icon("triangle-exclamation"), paste("Could not load this GEO dataset:", conditionMessage(result)))
        )
      } else {
        dataset$staged_expr <- result$expr
        dataset$staged_meta <- result$meta
        dataset$staged_source <- paste0("NCBI GEO: ", result$label)
        output$load_message <- renderUI(
          div(class = "empty-note", icon("check"),
              sprintf("Previewed %s genes across %s samples. This doesn't change what any sub-module runs on - go to Preprocessing and pick \"Currently loaded dataset\" to analyze it.",
                      format(nrow(result$expr), big.mark = ","), ncol(result$expr)))
        )
      }
    })
```
**Structurally identical logic to Block M (upload commit)**, with three GEO-specific differences:
1. `meta$sample <- rownames(meta)` instead of reading a mapped ID column (per the rationale in Block R's opening).
2. `validate(need(!inherits(em, "error"), "No GEO data fetched yet."))` — an extra guard specific to this path, re-checking that `geo_expr_meta()` didn't itself fail, even though the enabling `observe()` (just above) should already prevent the button from being clickable in that state — again, defense in depth.
3. `label` is built from `geo_fetch_result()$acc` + `em$platform` + collapse status, then used to build `dataset$staged_source <- paste0("NCBI GEO: ", result$label)` — a more information-dense provenance string than the upload path's (which just concatenates the two uploaded filenames) or the preloaded path's (a fixed descriptive string per entry).
**Thesis relevance:** High — same correctness points as Block M apply (sample intersection + `match()`-based reordering is the scientifically load-bearing step), plus this is a second, independent implementation of nearly the same logic — worth flagging in §7 as a candidate for refactoring/deduplication (not a correctness bug, but a maintainability observation you could legitimately note in a software-engineering-quality discussion).

---

## 5. DATA FLOW (full trace, all three paths)

### Path 1 — Preloaded dataset
```
User: selects a dataset in the dropdown (character, e.g. "GSE93272")
  → input$preloaded_choice (Shiny reactive input, string)
  → output$preloaded_note (renderUI) shows a provenance caveat (side-channel, informational only)
User: clicks "Load this dataset"
  → input$load_preloaded_btn (integer click counter) increments
  → observeEvent fires
  → Find() locates the matching PRELOADED_DATASETS entry (list)
  → entry$load() executes:
       - default entry → load_default_dataset() [global.R] → readRDS(DEFAULT_EXPR_RDS) (matrix) +
         data.table::fread(DEFAULT_META_CSV) (data.frame) → sample-intersected/reordered pair
       - individual entry → load_individual_dataset(gse_id) [global.R] → either a raw ExpressionSet's
         exprs()/harmonized pData(), a raw counts file (GSE89408), or a merged_training_subset()
         fallback (matrix + data.frame)
  → list(expr = <matrix>, meta = <data.frame>, source = <string>)
  → dataset$staged_expr / staged_meta / staged_source (reactiveValues fields; shared, session-scoped)
  → output$preloaded_load_message (renderUI) — user-visible confirmation
```

### Path 2 — Upload
```
User: chooses expr_file, meta_file (fileInput)
  → input$expr_file / input$meta_file (data.frame: name, size, type, datapath)
  → expr_raw()/meta_raw() (reactive) parse from datapath → matrix (expr) / data.frame (meta)
  → output$upload_preview_ui (renderUI) — row/col counts shown immediately
  → output$column_mapping (renderUI) — 4 selectInputs, defaults via guess_col()
User: (optionally) adjusts map_id/map_group/map_sex/map_batch
  → input$map_id / map_group / map_sex / map_batch (strings, or "(none)")
  → observe() → shinyjs::enable/disable("load_btn")
User: clicks "Upload Data"
  → input$load_btn increments → observeEvent fires
  → meta$sample/group/sex/batch harmonized (character vectors, NA_character_ if "(none)")
  → common <- intersect(colnames(expr), meta$sample) → character vector of matched sample IDs
  → validate(need(length(common) >= 4, ...)) — hard gate
  → expr <- expr[, common]; meta <- meta[match(common, meta$sample), ] — aligned pair
  → dataset$staged_expr / staged_meta / staged_source
  → output$load_message (renderUI) — success (+ duplicate-ID advisory) or error
```

### Path 3 — Live GEO fetch
```
User: types an accession, clicks "Fetch"
  → input$geo_accession (string) / input$geo_fetch_btn (counter)
  → geo_fetch_result() (eventReactive) → format-validates → GEOquery::getGEO() (network I/O)
    → list(acc = <string>, platforms = <list of ExpressionSet>) or an error-condition object
  → output$geo_platform_ui (renderUI) — shown only if multi-platform
User: (optionally) picks a platform
  → input$geo_platform_choice (string)
  → geo_eset() (reactive) → single ExpressionSet
  → geo_expr_meta() (reactive) → Biobase::exprs()/pData()/annotation() + optional
    collapse_probes_to_genes() [global.R, WGCNA::collapseRows] → list(expr=<matrix>, meta=<data.frame>,
    collapsed=<logical>, platform=<string>) or error object
  → output$geo_fetch_status (renderUI) — fetched N genes/probes x M samples, collapse status
  → output$geo_column_mapping (renderUI) — 3 selectInputs (no ID mapping needed)
User: clicks "Load this dataset" (GEO)
  → input$geo_load_btn increments → observeEvent fires
  → meta$sample <- rownames(meta); group/sex/batch harmonized
  → common <- intersect(...); validate(need(length(common) >= 4, ...))
  → expr/meta aligned via match()
  → dataset$staged_expr / staged_meta / staged_source (source string includes accession+platform+collapse status)
  → output$load_message (renderUI)
```

### Terminal stage (outside this file, for completeness)
```
dataset$staged_expr/meta/source (any of the 3 paths above)
  → mod_preprocessing.R: pp_preloaded_read("__current__", ...) reads
    dataset$staged_expr %||% dataset$expr (etc.) — staged takes priority if present
  → (user runs Preprocessing's merge/normalize/batch-correct pipeline)
  → user clicks "Use this as the active dataset" (mod_preprocessing.R:1928)
  → dataset$expr / meta / source overwritten — NOW every other sub-module (DGE, WGCNA, MR, ...)
    reads this as its input
```

---

## 6. RESULTS

| Output | What it shows | Source computation | Controlling input(s) | Reactive? | Updates when | If input missing | If invalid | If dataset empty | If filter/input changes |
|---|---|---|---|---|---|---|---|---|---|
| `preloaded_note` | Provenance caveat text | Static string keyed by `input$preloaded_choice` | `preloaded_choice` | Yes | dropdown selection changes | Renders nothing (`req()` blocks) | N/A (no invalid state possible — dropdown is closed-choice) | N/A | Text updates immediately, no data reload |
| `preloaded_load_message` | "Previewed X genes across Y samples..." | `entry$load()` result dimensions | `load_preloaded_btn` click + `preloaded_choice` | Yes (event-triggered) | Button clicked | N/A (`req()` on choice) | Would show whatever error `entry$load()`'s internal `validate()` raises (e.g. "Raw data ... not found on disk") | Message would report `0` genes/samples — not specially handled (see §7) | Only updates on click, not on dropdown change alone |
| `upload_preview_ui` | "Read N features x M samples..." | `nrow()`/`ncol()` of parsed `expr_raw()`/`meta_raw()` | `expr_file`, `meta_file` | Yes | either file changes | Blank (`req()`) | Shows "Could not read the uploaded file(s): <message>" | Would show `0` features/samples, no explicit block | N/A |
| `column_mapping` | 4 mapping dropdowns | `colnames(meta_raw())` + `guess_col()` | `meta_file` | Yes | metadata file changes | Not rendered | N/A | Dropdowns would have zero `choices` if `meta_raw()` had zero columns — edge case, see §7 | Regenerates fully (loses any manual override) any time `meta_file` changes |
| `load_message` (upload) | Success + gene/sample counts (+ duplicate-ID note) or error | `nrow()`/`ncol()` of final aligned `expr`; `duplicated(rownames(...))` | `load_btn` click | Yes (event) | Button clicked | N/A (`req()`) | "Could not load this dataset: <message>" — e.g., mapped column doesn't exist, or <4 matched samples | Blocked earlier by the `>= 4 common samples` gate before this renders success | Only on click |
| `geo_fetch_status` | Fetched N genes/probes x M samples, collapse status | `geo_fetch_result()` + `geo_expr_meta()` dimensions | `geo_fetch_btn` click, `geo_accession` | Yes (event via `req(input$geo_fetch_btn)`) | Fetch button clicked | Blank | "Could not fetch from GEO: <message>" or extraction error | "This GEO series has no expression matrix..." | Re-fetch requires a fresh click |
| `geo_load_btn`-triggered `load_message` | Same shape as upload success/error | Final aligned `expr`/`meta` dims | `geo_load_btn` click, `geo_map_group` | Yes (event) | Button clicked | N/A (`req()`) | "Could not load this GEO dataset: <message>" | Blocked by the `>= 4` gate | Only on click |

### User perspective
You pick one of three "sources," see an immediate, no-commitment preview of what was read, adjust column mapping if needed, then explicitly click a "Load" button to stage it. A clear, color-coded message always tells you whether it worked, and — deliberately — reminds you every single time that staging here does **not** change what any analysis module actually uses; you must go to Preprocessing and actively pick it up.

### Scientific/computational perspective
Every path converges on the same contract: a genes(rows) × samples(columns) numeric matrix, and a metadata data frame with harmonized `sample`/`group`/`sex`/`batch` columns, row-order-aligned to the matrix's column order via `match()`, gated on at least 4 overlapping samples. The three paths differ only in *where the matrix and metadata come from* (bundled RDS/CSV; user-uploaded files; a live NCBI GEO series-matrix fetch, optionally probe-to-gene collapsed via `WGCNA::collapseRows`).

---

## 7. CODE VALIDATION

### Correctness — generally sound, with specific caveats

- **Sample alignment (`intersect()` + `match()`)** — correct and is the single most scientifically important piece of logic in this file; verified independently in both the upload path (L268-275) and the GEO path (L422-428).
- **Lazy loading via closures (`individual_dataset_entry`)** — correct use of R closures; no shared-state bug (each closure captures its own `gse_id` by value at call time, not by reference to a loop variable — this file doesn't use a loop to build these, it uses `lapply()`, which is itself immune to the classic "closures over a loop variable" bug that a `for` loop would risk).
- **Reactive dependency structure** — appropriate use of `reactive()` (cached values), `eventReactive()`/`observeEvent()` (click-gated), and `observe()` (side-effect only) throughout; no obvious case of a `reactive()` that should have been an `eventReactive()` or vice versa.

### Issues found (documented per the Problem → Why → Evidence → Consequence → Test → Fix framework — **not applied**, per your instruction)

**Issue 1 — `expr_raw()`'s CSV branch has no numeric/type validation.**
- Problem: `as.matrix(m[, -1, drop = FALSE])` on a CSV with any non-numeric column (stray text, a mis-parsed header, an extra ID column) silently produces a **character matrix**, not a numeric one.
- Why it's a problem: A character matrix passes every check in this file (rownames get set, `common <- intersect(...)` still works on column names, `validate(need(length(common)>=4,...))` still passes) and gets staged successfully — the failure only surfaces later, deep in Preprocessing or an analysis module, as a confusing type-coercion error or silently-wrong arithmetic (e.g. string comparison instead of numeric).
- Evidence: `mod_dataset.R:171-175` — no `is.numeric()`/`mode()` check anywhere in `expr_raw()`.
- Consequence: A confusing downstream error far from its root cause; poor diagnosability.
- How to test: Upload a CSV expression matrix where one non-first column contains a text value (e.g., accidentally leaving a `"NA"` string or a stray annotation column); observe that `upload_preview_ui` still reports normal-looking dimensions with no warning.
- Possible fix (not applied): after `as.matrix()`, check `if (!is.numeric(m)) validate(need(FALSE, "Expression matrix must be entirely numeric (excluding the first ID column)."))`.

**Issue 2 — `expr_raw()`'s RDS branch has zero validation, unlike `meta_raw()`'s RDS branch.**
- Problem: `meta_raw()` validates `is.data.frame(d)` on an uploaded metadata RDS; `expr_raw()` has no equivalent check that the deserialized object is a matrix (or even numeric) at all.
- Why: Asymmetric defensiveness between two structurally parallel functions.
- Evidence: `mod_dataset.R:168-169` (`readRDS(path)` — no wrapping validation) vs. `L154-156` (`validate(need(is.data.frame(d), ...))`).
- Consequence: An RDS containing, say, a `data.frame` or a `list` instead of a `matrix` would fail later, opaquely, rather than at upload time with a clear message.
- Test: Upload an `.rds` file containing `saveRDS(data.frame(x=1:5), "bad.rds")` as the expression matrix; observe no error until something downstream calls e.g. `ncol()`/matrix-specific behavior that a data.frame *would* actually still satisfy in some cases — a genuinely more subtle failure than a hard crash, i.e. it might work by accident with wrong semantics.
- Possible fix (not applied): mirror `meta_raw()`'s pattern — `validate(need(is.matrix(d) || is.data.frame(d), "..."))`, coercing consistently.

**Issue 3 — The hardcoded minimum-sample threshold (`>= 4`) is duplicated in two places with no shared constant, and its scientific justification is not documented anywhere in this file.**
- Evidence: `mod_dataset.R:270` and `:424`, both literal `4`.
- Consequence for your thesis: You should be able to defend *why* 4, specifically, was chosen (a pragmatic floor "too few samples for any downstream statistic," or a specific test's minimum-cell-size requirement?) — **cannot determine the rationale from this file alone**; check whether it's documented elsewhere (e.g., a design doc, or inferred from the smallest group-comparison test used downstream) before stating a justification in your methodology chapter.
- Suggested validation test: attempt to load a dataset with exactly 3, exactly 4, and exactly 5 matched samples; confirm the boundary behaves as coded (`>= 4` means 4 passes, 3 fails).

**Issue 4 — `guess_col()`'s fallback for required fields (`map_id`, `map_group`) silently defaults to `cols[1]`, the first column, with no warning shown to the user if no match at all was found.**
- Problem: If a metadata file's group column is named something `guess_col()`'s candidate list doesn't cover (e.g. `"cohort_status_2024"`), the group dropdown silently defaults to whatever column happens to be first — potentially still the ID or an unrelated column — and nothing in the UI distinguishes "we found a confident match" from "we gave up and picked the first column."
- Why it's a problem: This is exactly the failure mode the file's own comment (L198-204) says `guess_col()` was built to *prevent* for the "no guess at all → position-based default" case — but the *fallback itself*, when triggered, reproduces the same risk one level down, just less often. A user who doesn't carefully check the dropdown before clicking "Upload Data" could stage a dataset with a meaningless `group` column, and the app would report success with a plausible-looking gene/sample count — no error anywhere.
- Evidence: `mod_dataset.R:205-211` (`fallback = cols[1]` default in the function signature is used, unmodified, by the `map_id`/`map_group` calls at L217-222 — no low-confidence indicator returned or displayed).
- Consequence: Silent, scientifically serious mis-mapping is possible and undetectable from the UI alone.
- Test: Upload metadata whose diagnosis/group column is named e.g. `"dx"` (not covered by the candidate list `c("group","diagnosis","disease","condition","status","phenotype")`, and not matched by substring either) with an ID-like column first; confirm the Group dropdown defaults to the wrong column and no warning appears.
- Possible fix (not applied): have `guess_col()` optionally return whether a real match was found (not just the fallback), and render a visible warning badge on the dropdown when a required field couldn't be confidently guessed.

**Issue 5 — `n_dup` (duplicate feature ID) check exists on the upload path but has no equivalent on the preloaded or GEO paths.**
- Evidence: `mod_dataset.R:287, 293-294` (upload only) — no analogous check in Block L (preloaded) or Block R (GEO).
- Consequence: A duplicate-row-name warning that's arguably just as relevant for a live GEO fetch (probe-to-gene collapse could in principle still leave duplicates depending on the annotation) is inconsistent across the three loading paths.
- Not necessarily a bug (the preloaded/collapsed paths may be much less likely to actually have duplicates), but worth verifying and noting as an inconsistency rather than assuming it's deliberate.

**Issue 6 — The `"__default_merged__"` sentinel ID and the `id`/`label` shape of `PRELOADED_DATASETS`/`preloaded_choices()` are duplicated conceptually in `mod_preprocessing.R` (`pp_cohort_choices()`, `pp_cohort_label()`, and a direct `identical(gse, default_dataset_entry$id)` check).**
- Not a correctness bug today (both files agree), but a maintainability/coupling point: if this file's sentinel or catalog shape changes, `mod_preprocessing.R` must be updated in lockstep, with no compiler/type system to catch a mismatch — purely a code-quality observation for a software-architecture discussion, not something to "fix" without being asked.

**Issue 7 — No explicit handling for a completely empty (0-row or 0-column) but successfully-parsed expression matrix in any of the three paths**, beyond the incidental `>= 4 common samples` gate (which would catch 0 samples, but a 0-*gene* matrix with valid sample columns would pass every check in this file and stage successfully with `nrow = 0`).
- Test: Craft/upload an expression CSV with a header row and gene-ID column but zero data rows below it; observe whether `upload_preview_ui` and the final success message report "Read 0 features x N samples" without any error.

**Issue 8 — Race/rapid-input consideration:** `observeEvent(input$load_btn, ...)` recomputes `expr_raw()`/`meta_raw()` (cheap, since they're memoized reactives, so no double file-parsing) but if a user changes `map_id`/`map_group` *while* a previous `load_btn` click's handler is still mid-execution (unlikely to matter in practice since Shiny is single-threaded per session — handlers run to completion before the next event is processed), there's no genuine race condition here. Included for completeness since your framework asks about rapid input changes — this is a case where the answer is "not a real risk, and here's why," which is itself worth stating explicitly rather than leaving unaddressed.

### None of the above are "the code doesn't run" bugs — the file is functionally coherent and internally consistent in its main path. The issues are primarily about **silent failure modes on malformed/edge-case input**, which is exactly the category your thesis's validation chapter should be interrogating.

---

## 8. VALIDATION TESTS

| # | Test | Expected behavior | How to verify | Pass/Fail (fill in) |
|---|---|---|---|---|
| 1 | Load "Merged Data" from the preloaded dropdown | `dataset$staged_expr`/`meta`/`source` populated; message reports the same gene/sample count as `load_default_dataset()` returns directly | In an R console with the app's environment loaded: `d <- load_default_dataset(); dim(d$expr)` vs. the UI's reported numbers | |
| 2 | Load "GSE93272" individually | Note text (Block F) correctly identifies it as a merged-cohort subset, not raw data (per `merged_training_subset()` fallback in `global.R`) — verify by checking whether `RAW_DIR/GSE93272_raw.rds` actually exists on disk | `file.exists(file.path(RAW_DIR, "GSE93272_raw.rds"))` — if `FALSE`, the fallback path is confirmed active | |
| 3 | Upload a valid expr/meta CSV pair with clean, well-named columns | `guess_col()` picks correct columns automatically; success message shows correct dimensions | Prepare a small synthetic CSV pair (5 genes x 6 samples, columns named `sample_id`, `diagnosis`, `sex`) and upload | |
| 4 | Upload expr/meta with only 3 overlapping sample IDs | Blocked with "Fewer than 4 sample IDs..." message; `dataset$staged_expr` unchanged from before the attempt | Deliberately mismatch 3+ sample names between the two files | |
| 5 | Upload expr matrix with duplicated gene IDs (rownames) | Loads successfully but shows the duplicate-ID advisory with the correct count | Craft a CSV with e.g. `"TP53"` appearing twice as the first-column value | |
| 6 | Upload metadata as an RDS containing a non-data-frame object | Rejected with "The uploaded metadata RDS file must contain a data frame." | `saveRDS(list(a=1), "bad_meta.rds")`, upload it | |
| 7 | Upload expr matrix RDS containing a non-numeric object (validates Issue 2 above) | Currently: no rejection at upload time — likely fails later or loads with wrong semantics | `saveRDS(data.frame(x=letters[1:5]), "bad_expr.rds")`, upload it, and trace what actually happens downstream | |
| 8 | Upload expr CSV with a non-numeric data column (validates Issue 1) | Currently: no rejection; matrix silently becomes character type | Add a stray text value in a numeric-looking data column | |
| 9 | Metadata with a group-like column named outside `guess_col()`'s candidate list (validates Issue 4) | Currently: silently defaults to `cols[1]`, no warning | Name the diagnosis column `"dx_2024"` and put a sample-ID-like column first | |
| 10 | Fetch an invalid GEO accession (e.g. `"GSM12345"` or `"not-an-id"`) | Blocked immediately with "Enter a valid GEO Series accession..." — no network call made | Type it into the accession box and click Fetch; confirm (e.g. via logging or a network monitor) that `GEOquery::getGEO()` was never actually invoked | |
| 11 | Fetch a real multi-platform GEO series | Platform picker appears; `geo_load_btn` stays disabled until a platform (and group column) is chosen | Requires a live network call to NCBI GEO and a known multi-platform accession | |
| 12 | Fetch a GEO series with no expression matrix in its series-matrix file (RNA-seq, counts-only) | "This GEO series has no expression matrix..." error shown | Requires a live network call to a known counts-only RNA-seq GEO series | |
| 13 | Zero-row (0 genes) but valid-column-count expression matrix (validates Issue 7) | Currently: no explicit block; loads "successfully" with 0 genes | Craft a header-only CSV with a valid gene-ID column name but no data rows | |
| 14 | Rapid repeated clicking of "Load this dataset" (preloaded) | No duplicate/garbled state; last click's result is what's staged | Click twice quickly in the running app and confirm the message/dimensions reflect a single coherent load | |
| 15 | GEOquery not installed in the deployment | "The GEOquery package is not installed..." shown instead of a crash | Test in an environment where `GEOquery` is intentionally absent, or simulate by checking the `requireNamespace()` branch logic directly | |

Example test-code skeleton for a local/offline test (tests 3-9, run outside Shiny against the reactive logic's underlying pure operations):
```r
# Minimal reproduction of the upload path's core alignment logic (L260-277),
# extracted for standalone testing outside the reactive graph.
test_align <- function(expr, meta, map_id, map_group) {
  meta$sample <- as.character(meta[[map_id]])
  meta$group  <- as.character(meta[[map_group]])
  common <- intersect(colnames(expr), meta$sample)
  stopifnot(length(common) >= 4)
  expr <- expr[, common, drop = FALSE]
  meta <- meta[match(common, meta$sample), , drop = FALSE]
  list(expr = expr, meta = meta)
}

expr <- matrix(rnorm(30), nrow = 5, dimnames = list(paste0("g", 1:5), paste0("S", 1:6)))
meta <- data.frame(sample_id = paste0("S", 1:6), dx = rep(c("RA", "HC"), 3))
out <- test_align(expr, meta, "sample_id", "dx")
stopifnot(identical(colnames(out$expr), out$meta$sample))  # order-alignment check
```

---

## 9. THESIS DOCUMENTATION

### Methodological concept
Multi-source dataset ingestion and harmonization for a transcriptomics analysis pipeline: reconciling three heterogeneous data-acquisition modes (bundled reference cohorts, arbitrary user uploads, live public-repository queries) into one internally consistent representation (gene × sample numeric matrix + harmonized sample metadata) before any downstream statistical analysis.

### Implementation
Implemented as a Shiny module (`mod_dataset_ui`/`mod_dataset_server`) using reactive programming: lazy, memoized reactive expressions for parsing (`expr_raw`, `meta_raw`, `geo_expr_meta`), event-gated reactive/observer pairs for user-triggered actions (`eventReactive`/`observeEvent`), and a shared, session-scoped `reactiveValues` object (`dataset`) as the cross-module state-sharing mechanism, deliberately split into "staged" vs. "active" fields to make dataset activation an explicit, auditable, single-button action elsewhere in the app (Preprocessing).

### Inputs
An uploaded CSV/RDS expression matrix (genes × samples) and CSV/RDS sample metadata (one row per sample); OR a GEO Series accession (`GSExxxxx`) resolved via `GEOquery::getGEO()`; OR one of five bundled reference datasets (a merged/batch-corrected training cohort, or four individual GEO sources) read from local `.rds`/`.csv` files established elsewhere in the codebase (`global.R`).

### Processing
Format detection by file extension; type coercion (data.frame → matrix with gene-symbol rownames for expression data); metadata column harmonization to a fixed schema (`sample`, `group`, `sex`, `batch`) via heuristic name-matching (`guess_col()`) with user-editable defaults; sample-set intersection (`intersect()`) and order-alignment (`match()`) between expression columns and metadata rows; a minimum-sample-count gate (≥4 matched samples); for GEO data, optional probe-to-gene collapsing via `WGCNA::collapseRows` (MaxMean rule) when platform gene-symbol annotation is available.

### Outputs
Three fields (`dataset$staged_expr`, `staged_meta`, `staged_source`) on a shared reactive object, consumed by the Preprocessing sub-module as a candidate dataset; plus user-facing preview/status/error messaging at every stage.

### Reproducibility — what to document
- Exact package versions: `shiny`, `shinyjs`, `data.table`, `GEOquery`, `Biobase`, `WGCNA` (record via `sessionInfo()`/`renv.lock`/`DESCRIPTION`, whichever this project uses — **cannot determine dependency-pinning mechanism from this file alone**).
- `GEOquery::getGEO(acc, GSEMatrix = TRUE)` — the exact function call and parameter, since it determines that only series-matrix-deposited (typically array-normalized) data is retrieved, never raw platform files.
- The `WGCNA::collapseRows(..., method = "MaxMean", connectivityBasedCollapsing = FALSE, ...)` exact parameterization (`global.R:711-727`), cited per Miller et al. 2011 (already referenced in that function's own comment) — reuse that citation.
- The hardcoded `>= 4` minimum-sample threshold, and your best-available justification for it (see Issue 3, §7 — verify before asserting a specific rationale).
- The fixed metadata schema (`sample`/`group`/`sex`/`batch`) as the harmonization target.

### Limitations to acknowledge
- Column-mapping defaults are heuristic (substring/exact name matching), not guaranteed-correct, and a wrong guess is not distinguished from a confident guess in the UI (Issue 4).
- No numeric-type validation on uploaded expression matrices (Issues 1-2) — malformed uploads can be staged without error and fail later, non-obviously.
- Individual "raw" GEO source selections for GSE93272/GSE110169 in this deployment are not actually raw, single-platform data, but a post-hoc subset of the already-merged/batch-corrected cohort — a data-provenance caveat that must be stated wherever you describe using "the individual training sources."
- GEO-fetched data reflects whatever normalization the original depositor applied to their series matrix — not independently re-normalized by this app.
- Duplicate feature-identifier detection exists only on the manual-upload path, not the other two.

### Validation
Demonstrate correctness via the tests in §8 — particularly the order-alignment property (`identical(colnames(expr), meta$sample)` after loading) as the single most important correctness invariant to show holds for all three ingestion paths, plus the documented failure-mode tests (Issues 1, 2, 4, 7) run against controlled malformed inputs to characterize — not necessarily "fix" — current behavior, which is itself valid thesis content (documenting known limitations you observed empirically).

---

## 10. CODE → THESIS MAPPING

| Code section | Functionality | Scientific purpose | Thesis section | What to document |
|---|---|---|---|---|
| L44-77 (`individual_dataset_entry`, `PRELOADED_DATASETS`) | Catalog + lazy loader for bundled reference cohorts | Defines the app's fixed cohort universe (1 merged + 4 individual GEO sources) | Data / Methods — Study Cohorts | The 4 GEO accessions, their roles (training/validation), and the GSE93272/GSE110169 raw-vs-fallback caveat |
| L150-177 (`meta_raw`, `expr_raw`) | File parsing + type coercion | Establishes the internal data representation contract | Implementation — Data Ingestion | Expected file format (genes-in-rows/first-column-ID convention), RDS vs CSV branches, validation gaps (Issues 1-2) |
| L198-211 (`guess_col`) | Heuristic column-name inference | Reduces user burden while avoiding silent positional mis-mapping | Implementation — Metadata Harmonization | Exact-then-substring matching algorithm; required vs. optional field fallback behavior; Issue 4 limitation |
| L257-298 (`load_btn` handler) | Sample-set intersection + order alignment + staging | Guarantees expr/meta correspondence before any statistic can be computed on them | Implementation — Data Integrity / Methods — QC | `intersect()`+`match()` alignment logic; the `>=4` sample gate; duplicate-ID detection |
| L300-322 (`geo_fetch_result`) | Live NCBI GEO query via GEOquery | Enables ad hoc validation on newly published cohorts, not just bundled ones | Methods — External Data Acquisition | `getGEO(acc, GSEMatrix=TRUE)` semantics; graceful degradation if GEOquery is absent |
| L343-362 (`geo_expr_meta`) | Probe-to-gene collapsing decision logic | Normalizes fetched data into the same gene-level shape used elsewhere in the pipeline | Methods — Probe/Gene-Level Harmonization | `collapse_probes_to_genes()`/`WGCNA::collapseRows` MaxMean rule; the `used_collapse` triple condition; Miller et al. 2011 citation |
| L1-14 (header comment) + `dataset$staged_*` vs `dataset$expr` split throughout | Staged-preview vs. active-dataset separation | Prevents any preview action from silently altering results already computed elsewhere in the app | Software Architecture — State Management / Reproducibility | The reactiveValues sharing pattern; why activation is a separate, explicit step gated in Preprocessing |

---

## 11. PACKAGES AND DEPENDENCIES

| Package | Function(s) used | What it does | Why needed here | Without it | Core / optional |
|---|---|---|---|---|---|
| `shiny` | `NS`, `moduleServer`, `reactive`, `reactiveValues` (via `dataset`, constructed in `server.R`), `observe`, `observeEvent`, `eventReactive`, `renderUI`, `uiOutput`, `req`, `validate`, `need`, `selectInput`, `textInput`, `fileInput`, `actionButton`, `tagList`, `div`, `p`, `span`, `icon` (from `shiny`'s bundled Font Awesome helper) | Core reactive-web-app framework | The entire file's structure | The module cannot exist at all | Core |
| `shinyjs` | `shinyjs::enable`, `shinyjs::disable` | JS-backed DOM manipulation from R (toggling button `disabled` state) | UX guardrail on `load_btn`/`geo_load_btn` | Buttons would always be clickable; would rely solely on internal `req()` guards (still functionally safe, just less obvious to the user) | Optional-ish (cosmetic/UX, not correctness-critical) |
| `data.table` | `data.table::fread` | Fast, C-backed delimited-file reader | Parses uploaded/bundled CSVs efficiently, including potentially large expression matrices | Would need base `read.csv`, likely slower on large files | Core (performance-relevant, not just a convenience) |
| `GEOquery` (Bioconductor) | `GEOquery::getGEO` | Downloads/parses NCBI GEO Series Matrix files into Bioconductor `ExpressionSet` objects | Powers the entire "Fetch from NCBI GEO" feature | That entire UI panel degrades to an explanatory error message (`requireNamespace()` check, L308-310) — rest of the app still works | **Optional** — explicitly checked and gracefully degraded, unlike the others |
| `Biobase` (Bioconductor) | `Biobase::exprs`, `Biobase::pData`, `Biobase::annotation`, `Biobase::fData` (the last one used inside `collapse_probes_to_genes()`, `global.R`) | Accessor functions for Bioconductor's `ExpressionSet` S4 class | Extracts expression matrix / sample metadata / platform / feature metadata from a fetched GEO object | GEO fetch path couldn't extract data even if `GEOquery` itself succeeded | Core to the GEO path (implicitly required alongside `GEOquery`, though not itself `requireNamespace()`-checked in this file — **cannot confirm from this file alone whether `Biobase` is separately guarded elsewhere**, e.g. at app startup) |
| `WGCNA` (via `global.R`'s `collapse_probes_to_genes()`, not called directly in this file) | `WGCNA::collapseRows` | Collapses a probe-level matrix to one row per gene (MaxMean rule) | Used indirectly through `geo_expr_meta()`'s call to `collapse_probes_to_genes()` | GEO-fetched microarray data would stay probe-level with no gene-symbol collapse option from this file's own logic | Core to that specific sub-feature, but only indirectly a dependency of this file |

### Dependencies on other files in this application
- `global.R`: `load_default_dataset()`, `load_individual_dataset()`, `GEO_SOURCES`, `collapse_probes_to_genes()`, `%||%`. All confirmed present and matching this file's usage by direct inspection.
- `server.R`: constructs the shared `dataset <- reactiveValues(...)` object and calls `mod_dataset_server("tx_dataset", dataset)` — this file cannot function standalone; it requires being instantiated with a pre-existing `reactiveValues` object shaped with (at least) `expr`/`meta`/`source`.
- `mod_preprocessing.R`: the actual consumer of `dataset$staged_*`; not a dependency *of* this file, but the essential next file to read to understand what happens to what this file produces (see §12).
- CSS classes referenced (`upload-step-label`, `empty-note`, `header-dataset-badge`) are defined elsewhere (likely a global stylesheet) — **cannot determine from this file alone** where.

---

## 12. ARCHITECTURE

```
global.R (sourced first: constants, GEO_SOURCES, load_default_dataset(),
          load_individual_dataset(), collapse_probes_to_genes(), %||%)
        │
        ▼
server.R
  dataset <- reactiveValues(expr=, meta=, source=)   [from load_default_dataset() at app startup]
        │
        ├──► mod_dataset_server("tx_dataset", dataset)   ← THIS FILE (writes dataset$staged_*)
        │
        └──► (other Transcriptomics sub-modules, e.g. mod_preprocessing.R, mod_dge.R, mod_wgcna.R, ...
              all also handed the same `dataset` object, per server.R's registration pattern)
                    │
                    ▼
              mod_preprocessing.R
                pp_preloaded_read("__current__", ...) reads dataset$staged_* %||% dataset$expr/meta/source
                → merge / normalize / batch-correct pipeline
                → "Use this as the active dataset" button writes dataset$expr/meta/source
                        │
                        ▼
              every other Transcriptomics sub-module (DGE, WGCNA, MR, feature selection, ...)
              reads dataset$expr/meta/source as its analysis input
```

- `ui.R`: calls `mod_dataset_ui("tx_dataset")` somewhere in the overall app UI assembly — **not inspected in this session**; read it next if you need to confirm exactly how/where this tab is mounted in the sidebar (see `mod_dataset_config`, Block B, which is very likely read by a registry-driven navigation builder rather than `ui.R` hardcoding this tab directly — **cannot confirm the exact mechanism without reading `ui.R`/`submodules_registry.R`**).
- The parallel structure in Methylomics (`methyl_dataset`, `mod_methyl_dataset_server`) and Cross-Omics (`cross_dataset`) confirms this "one shared reactiveValues per omics area, with a dedicated Dataset sub-module writing to it" is a repeated architectural pattern across the whole app, not unique to Transcriptomics — worth a sentence in your Software Architecture chapter describing it as the app's general convention.

### Files to inspect next, and why
1. **`mod_preprocessing.R`** (already partially inspected in this session) — the direct, sole consumer of `dataset$staged_*`, and the only place `dataset$expr/meta/source` is ever overwritten after startup. Essential to fully understand the "staging → activation" lifecycle this file only half-implements.
2. **`server.R`** (already partially inspected) — confirms `dataset`'s exact initial shape and how `mod_dataset_server` is wired in; also shows the header badge (`active_dataset_badge`) and `results` cache invalidation that react to `dataset$source` changing.
3. **`ui.R`** — to confirm exactly how `mod_dataset_ui("tx_dataset")` is mounted (directly, or via the `MODULE_REGISTRY`/`mod_dataset_config` pattern referenced in this file).
4. **`global.R`** in full (only the relevant excerpts were read here) — particularly `data_paths.R` (referenced but not opened this session) for where `DEFAULT_EXPR_RDS`/`DEFAULT_META_CSV`/`RAW_DIR` actually point, if you want to document the physical data provenance completely.

---

## 15. FINAL LEARNING SUMMARY

### A. What this file does (plain English)
It's the single "front door" for getting an expression matrix + sample metadata into the app, from three different sources (bundled cohorts, user uploads, live GEO fetches), reshaping all three into one consistent, aligned, harmonized representation, and staging that as a *candidate* dataset that another module (Preprocessing) must explicitly promote to "active" before it affects any real analysis anywhere else in the app.

### B. What you learned (R/Shiny concepts)
- Shiny modules: `NS()`, `moduleServer()`, `session$ns`, input/output namespacing.
- Reactive programming primitives and when each is correct: `reactive()` (cached value), `renderUI()`/`uiOutput()` (dynamic UI), `observe()` (side effect, any-dependency-change), `observeEvent()`/`eventReactive()` (click-gated action/value).
- `req()` and `validate(need(...))` as Shiny's two idioms for "stop gracefully" vs. "stop with a user-visible reason."
- The sentinel-object error pattern (`tryCatch(..., error = function(e) e)` + `inherits(x, "error")`) as a way to route errors through the reactive graph as ordinary data rather than exceptions.
- R closures (functions capturing their defining environment) as a lazy-evaluation mechanism (`individual_dataset_entry`).
- Core R data-shape operations: `intersect()`, `match()` for aligning two datasets by key; `vapply()`/`lapply()` for type-safe iteration; `drop = FALSE` as a correctness-critical guard against R's automatic dimension-dropping.
- The distinction between a shared, cross-module `reactiveValues` object (App/session state) versus a local `reactive()` (module-private cached computation).

### C. What you should be able to explain at your defense
- Why "staged" and "active" datasets are kept as separate fields, and what problem that separation prevents.
- Exactly how sample alignment between an expression matrix and its metadata is guaranteed (or *not* guaranteed, in the failure modes you identified) — the `intersect()`/`match()` mechanism, and the `>=4` sample floor.
- What "raw" actually means for each of the four bundled GEO sources in this dropdown — specifically, that two of them are not truly raw in this deployment.
- What a live GEO fetch via `GEOquery::getGEO()` actually retrieves (submitter-processed series-matrix values, not raw platform files), and when/why probe-to-gene collapsing is or isn't applied.
- The known, empirically-testable limitations you found (heuristic column mapping with no confidence signal; no numeric-type validation on uploads) and why they represent reasonable engineering trade-offs (usability vs. exhaustive validation) rather than oversights you're unaware of.

### D. Validation status
- ✅ Sample-alignment logic (`intersect()` + `match()`), both paths — correct and consistent.
- ✅ Lazy-loading closures — correct, no shared-state bug.
- ✅ Optional-dependency handling for `GEOquery` (`requireNamespace()`) — correct graceful degradation.
- ⚠️ Needs verification: exact rationale for the `>= 4` sample minimum; whether `Biobase` has its own `requireNamespace()` guard elsewhere; exact CSS class definitions' location; where `ui.R` mounts this module.
- ❌ Potential bug/gap: no numeric-type validation on uploaded expression matrices (CSV or RDS) — Issues 1-2.
- ❌ Potential bug/gap: silent low-confidence fallback in `guess_col()` for required fields, with no visible confidence signal — Issue 4.
- ❌ Potential bug/gap: no explicit handling of a zero-row (0-gene) but otherwise valid dataset — Issue 7.
- 🔬 Scientific decision requiring justification: the `>= 4` minimum matched-sample threshold.
- 🔬 Scientific decision requiring justification: MaxMean-rule probe collapsing (`WGCNA::collapseRows`) as the probe-to-gene reduction method for live GEO fetches — already partially justified by the Miller et al. 2011 citation in `global.R`, but worth restating explicitly for the GEO-fetch-specific use case.
- 📌 Documentation requirement: state clearly, wherever GSE93272/GSE110169 individual loading is discussed, that in this deployment they are derived from the merged/batch-corrected cohort, not independently raw.
- 📌 Documentation requirement: state the fixed metadata schema (`sample`/`group`/`sex`/`batch`) as this file's harmonization contract.

### E. Questions to investigate further
1. What is the actual git-history rationale for the `>= 4` sample minimum — was it chosen for a specific downstream statistical test's minimum group size?
2. Does any other module in this app perform numeric-type validation on `dataset$expr`/`staged_expr` before running arithmetic on it, effectively catching Issue 1/2 downstream (even if later than ideal)?
3. Is `Biobase` itself guarded by a `requireNamespace()` check anywhere (e.g., at app startup in `global.R`), or does the app assume it's always installed alongside `GEOquery`?
4. Where exactly does `ui.R` mount `mod_dataset_ui("tx_dataset")`, and is `mod_dataset_config` (Block B) actually consumed by a registry-driven UI builder, or currently unused/dead metadata?
5. Does `mod_preprocessing.R`'s merge step really "keep only the first occurrence" of duplicated feature IDs, exactly as this file's own advisory message (L294) claims — worth verifying directly in that file's merge logic rather than trusting the message string.

### F. Next file to study
**`mod_preprocessing.R`** — it is the direct, necessary continuation of this file's data flow: it reads `dataset$staged_*`, runs the merge/normalization/batch-correction pipeline, and contains the *only* code path that ever promotes a staged dataset to "active" (`dataset$expr/meta/source`). You cannot fully understand or defend `mod_dataset.R` in isolation without tracing what happens to its output immediately afterward — and several of this file's own UI messages (e.g. "go to Preprocessing and pick 'Currently loaded dataset'") are direct, literal pointers to exactly that file.
