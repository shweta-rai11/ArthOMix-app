# Sex Interaction Analysis Module: `mod_interaction.R`

**Source file:** `ArthOMix/R/transcriptomics/mod_interaction.R` (121 lines)
**Registration:** `mod_interaction_config` (id = `"interaction"`, group = "Biomarker modeling", title = "Sex Interaction Analysis", icon = `"venus-mars"`).
Prepared: 2026-08-25

This document is derived **exclusively** from the code in `mod_interaction.R`. Nothing here is inferred from the thesis methods, other transcriptomics modules, external literature, or general `limma` capability beyond what is actually called. Where the module's own header comment is quoted, it is attributed as such.

---

## 1. Internal audit mapping

`UI element → input ID → server observer/reactive → analysis function → output ID → displayed result`

| UI element | Input ID | Server reactive/observer | Analysis function called | Output ID | Displayed result |
|---|---|---|---|---|---|
| "Reference group" `selectInput` | `ref_group` | `fit_result` (`eventReactive`) | N/A (subsetting) | `summary_ui`, `int_table` | feeds the model design |
| "Comparison group" `selectInput` | `comp_group` | `fit_result` | N/A | `summary_ui`, `int_table` | feeds the model design |
| "Reference sex" `selectInput` | `ref_sex` | `fit_result` | N/A | `summary_ui`, `int_table` | feeds the model design |
| "Comparison sex" `selectInput` | `comp_sex` | `fit_result` | N/A | `summary_ui`, `int_table` | feeds the model design |
| (all four above rendered by) | `controls` | `output$controls` (`renderUI`) | reads `dataset$meta$group`, `dataset$meta$sex` | `ns("controls")` | the four `selectInput`s |
| "Adjusted p-value cutoff" `numericInput` | `padj_cut` | `sig_table` (`reactive`) | `mutate(significant = adj.P.Val < padj_cut)` | `int_table`, `summary_ui` | `significant` flag; significant-count text |
| "Run interaction model" `actionButton` | `run_btn` | `int_has_run` (`reactiveVal`, `observeEvent`); triggers `fit_result` (`eventReactive`) | `model.matrix`, `limma::lmFit`, `limma::eBayes`, `limma::topTable` | `summary_ui`, `int_table` | reveals result summary and table |
| "Result" box body | N/A | `output$summary_ui` (`renderUI`) | N/A | `ns("summary_ui")` | pre-run notice, or interaction term + gene counts |
| "Interaction result table" box | N/A | `output$int_table` (`DT::renderDataTable`) | N/A | `ns("int_table")` | `DT::datatable` of `sig_table()` |
| "Download CSV" `downloadButton` | `download_int` | `output$download_int` (`downloadHandler`) | `write.csv(sig_table(), ...)` | `ns("download_int")` | `sex_interaction.csv` file |

**No sub-tabs exist in this module.** The UI (`mod_interaction_ui`, lines 13–40) contains no `tabsetPanel`/`tabPanel`/`navset*` call; it is a single `fluidRow` of two boxes ("Model", "Result") followed by one full-width box ("Interaction result table").

---

## 2. Code-validation checklist

| # | Documented feature | Location in `mod_interaction.R` | Status |
|---|---|---|---|
| 1 | Single unified view, no sub-tabs | lines 15–39 (no `tabsetPanel`/`tabPanel` anywhere in file) | CONFIRMED |
| 2 | Reference/Comparison group `selectInput`s | lines 52–53 | CONFIRMED |
| 3 | Reference/Comparison sex `selectInput`s | lines 54–55 | CONFIRMED |
| 4 | Groups/sexes populated from `dataset$meta$group` / `dataset$meta$sex` | lines 47–48 | CONFIRMED |
| 5 | Validation: ≥2 distinct group values, ≥2 distinct sex values before controls render | lines 49–50 | CONFIRMED |
| 6 | `padj_cut` numeric input, default 0.05, range 0–1, step 0.01 | line 22 | CONFIRMED |
| 7 | "Run interaction model" action button | line 23 | CONFIRMED |
| 8 | `int_has_run` flag set on button click | lines 59–60 | CONFIRMED |
| 9 | Model fit gated behind button click (`eventReactive`) | line 62 | CONFIRMED |
| 10 | `req()` guard on the four selectors | line 63 | CONFIRMED |
| 11 | Validation: reference ≠ comparison, for group and for sex | lines 64–65 | CONFIRMED |
| 12 | Metadata filtered to the two chosen groups and two chosen sexes, missing sex excluded | lines 67–69 | CONFIRMED |
| 13 | Sample intersection between metadata and expression matrix columns | line 70 | CONFIRMED |
| 14 | Validation: ≥12 matching samples | line 71 | CONFIRMED |
| 15 | Metadata reordered / expression subset to matched samples | lines 72–73 | CONFIRMED |
| 16 | `grp`/`sx` factors built with reference level first | lines 75–76 | CONFIRMED |
| 17 | Validation: every group×sex cell has ≥2 samples | line 77 | CONFIRMED |
| 18 | `design <- model.matrix(~ grp * sx)` | line 79 | CONFIRMED |
| 19 | `limma::eBayes(limma::lmFit(expr, design))` | line 80 | CONFIRMED |
| 20 | Interaction coefficient = last column of design matrix | line 81 | CONFIRMED |
| 21 | `limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P")` | line 82 | CONFIRMED |
| 22 | `gene` column added from expression matrix row names, moved to first column | lines 83–85 | CONFIRMED |
| 23 | `significant` column = `adj.P.Val < padj_cut`, recomputed reactively without refit | lines 89–93 | CONFIRMED |
| 24 | Pre-run notice text in Result box | lines 96–99 | CONFIRMED |
| 25 | Post-run summary: interaction term name, genes tested, genes significant | lines 100–106 | CONFIRMED |
| 26 | Results `DT::datatable`: `rownames=FALSE`, `filter="top"`, `pageLength=15`, `scrollX=TRUE`, class `"stripe hover compact"`, gated on `int_has_run()` | lines 109–113 | CONFIRMED |
| 27 | CSV download of `sig_table()`, filename `sex_interaction.csv` | lines 115–118 | CONFIRMED |
| 28 | `results` argument accepted by `mod_interaction_server(id, dataset, results = NULL)` but never referenced in the function body | line 42 (signature); no other occurrence of `results` in file | CONFIRMED (unused parameter) |
| 29 | Registration metadata: id `"interaction"`, group "Biomarker modeling", title "Sex Interaction Analysis", icon `"venus-mars"` | lines 6–11 | CONFIRMED (registration only, not UI/server logic) |
| 30 | Header comment ("fits a live diagnosis-by-sex interaction model...group*sex interaction") matches the implemented `model.matrix(~ grp * sx)` / `limma` pipeline | lines 1–4, 9 vs. lines 79–82 | CONFIRMED: comment and code agree; no unimplemented comment content found |
| 31 | No plots (`ggplot`, `plotly`, base graphics) are generated anywhere in the module | full-file scan | CONFIRMED (absent) |
| 32 | No statistical test beyond `limma::lmFit`/`eBayes`/`topTable` (e.g., no separate per-sex contrasts, no Wilcoxon test, no volcano plot) | full-file scan | CONFIRMED (absent) |
| 33 | Only one downloadable output (`download_int`); no `.rds` or other export | full-file scan | CONFIRMED (absent) |

---

## 3. Thesis subsection (as delivered)

### 2.10 Sex Interaction Analysis

**Purpose and implementation.**
Implemented in `mod_interaction.R` as a single-view submodule (registered as id `"interaction"`, group "Biomarker modeling", title "Sex Interaction Analysis", icon `venus-mars`) that fits one live `limma` linear model with a group-by-sex interaction term on the currently loaded dataset and reports which genes show a statistically different group effect between two selected sex levels. The module implements no sub-tabs: the interface is a single "Model" configuration box, a "Result" summary box, and one full-width results-table box.

**Inputs.**
- Loaded dataset: `dataset$meta` (columns `group`, `sex`, `sample`) and `dataset$expr` (gene-by-sample expression matrix).
- Reference group (`ref_group`) and comparison group (`comp_group`) select inputs, populated from the distinct non-missing values of `dataset$meta$group`.
- Reference sex (`ref_sex`) and comparison sex (`comp_sex`) select inputs, populated from the distinct non-missing values of `dataset$meta$sex`.
- Adjusted p-value cutoff (`padj_cut`, numeric, default 0.05, range 0–1, step 0.01).
- "Run interaction model" action button (`run_btn`).
- The `results` argument accepted by `mod_interaction_server(id, dataset, results = NULL)` is not used anywhere in the module.

**Analysis workflow.**
On clicking `run_btn`, an `eventReactive` (`fit_result`) executes the following, exactly as coded:
1. Requires the four selectors to be set; validates that the reference and comparison differ for both group and sex.
2. Subsets `dataset$meta` to samples whose `group` is one of the two chosen groups and whose `sex` is non-missing and one of the two chosen sexes.
3. Intersects sample IDs between the filtered metadata and the columns of `dataset$expr`; validates at least 12 matching samples.
4. Reorders the metadata to match the retained samples and subsets the expression matrix to the same samples.
5. Builds `grp` and `sx` factors (reference level first) and validates that every group-by-sex combination contains at least 2 samples.
6. Fits `design <- model.matrix(~ grp * sx)`, then `limma::eBayes(limma::lmFit(expr, design))`.
7. Takes the interaction coefficient as the last column of the design matrix and extracts results with `limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P")`.
8. Adds a `gene` column from the expression matrix's row names.

A separate reactive (`sig_table`) appends a `significant` column (`adj.P.Val < padj_cut`) to the fitted table. Because this step is a plain reactive over `fit_result()` and `input$padj_cut`, changing the cutoff after a run re-flags significance without refitting the model.

**Sub-tabs and interface.**
The module has no sub-tabs. Its interface consists of three boxes:
- **Model**: the four group/sex select inputs, the `padj_cut` numeric input, and the `run_btn` button.
- **Result**: a spinner-wrapped text panel (`summary_ui`). Before a run it shows the notice "Not run yet. Set the groups and sexes on the left, then click 'Run interaction model'." After a run it reports the interaction coefficient name and the counts of genes tested and genes significant at the current cutoff.
- **Interaction result table**: a "Download CSV" button and a `DT` table, both gated on a run having completed (`int_has_run()`).

**Outputs.**
- A results table (`int_table`) rendered with `DT::datatable()` (top-row column filters, 15 rows per page, horizontal scrolling, "stripe hover compact" styling), listing `gene`, the full set of `limma::topTable` statistics, and the derived `significant` flag.
- A CSV download (`download_int`, filename `sex_interaction.csv`) of the same table.
- A text summary (`summary_ui`) giving the interaction term name, number of genes tested, and number of genes significant at the current adjusted-p cutoff.
- No plots are generated by this module.

**Implementation summary.**
`mod_interaction.R` registers as one submodule in the "Biomarker modeling" group and exposes a single, non-tabbed interface (a model-configuration box, a text-summary box, and a downloadable results table), driven end-to-end by one button-triggered `limma` group-by-sex interaction fit on the dataset currently loaded in the application.

---

## 4. XomicShiny-style paragraph (as delivered)

### 4.1 Short version

> The Sex Interaction Analysis module tests, on the dataset currently loaded in the app, whether a two-group comparison behaves differently between sexes, by fitting a group-by-sex interaction model with `limma`. Its input data is the loaded dataset's expression matrix and sample metadata (sample IDs, group labels, and sex labels), together with the user's chosen reference and comparison group, reference and comparison sex, and adjusted p-value threshold (default 0.05), selected via the "Model" panel and submitted by clicking "Run interaction model." Its output data is a per-gene results table listing each gene's interaction statistics and a significance flag at the chosen threshold, a text summary reporting the interaction term and the counts of genes tested and significant, and a downloadable CSV file of the results table; the module does not include any plots. This complements the app's other sex-agnostic group comparisons by testing the interaction term directly, so a gene can be flagged here as sex-differential even when a plain group comparison would not distinguish it, the question the module is built to answer per its own registration description ("which genes respond to the group difference differently in each sex").

### 4.2 Labeled version

> **Sex Interaction Analysis**
>
> *Purpose.* The Sex Interaction Analysis submodule fits a live linear model with a group-by-sex interaction term on the currently loaded dataset, identifying genes whose response to a two-group contrast differs between two selected sex levels.
>
> *Web-app implementation.* Implemented as `mod_interaction.R` ("Sex Interaction Analysis," Biomarker modeling group). The interface is a single, non-tabbed view: a "Model" box for configuration, a "Result" box for a text summary, and a full-width "Interaction result table" box.
>
> *Inputs.* Reference and comparison group, and reference and comparison sex, all populated from the loaded dataset's sample metadata; an adjusted p-value cutoff (default 0.05); and a "Run interaction model" button.
>
> *Processing.* On running, the module restricts samples to the two chosen groups and two chosen sexes with non-missing sex values, requires at least 12 matching samples and at least 2 samples in every group-by-sex combination, and fits `design <- model.matrix(~ group * sex)` via `limma::lmFit` and `limma::eBayes`. Results are extracted with `limma::topTable()` on the interaction coefficient (the last column of the design matrix), and a significance flag is applied reactively at the user-chosen adjusted p-value cutoff without refitting the model.
>
> *Outputs.* A text summary reporting the interaction term, the number of genes tested, and the number of genes significant at the current cutoff; a filterable, sortable table of all tested genes with their interaction statistics and significance flag; and a CSV download of that table. No plots are produced.
>
> *Workflow.* The user selects a reference/comparison group pair, a reference/comparison sex pair, and a p-value cutoff, then clicks "Run interaction model"; the module builds the matched sample subset, fits the interaction model once, and displays the summary and results table, both of which update instantly if the p-value cutoff is subsequently changed.
