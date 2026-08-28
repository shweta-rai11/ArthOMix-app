# Methylomics — Colocalization

**Source file:** `ArthOMix/R/methylomics/mod_methyl_coloc.R` (1,221 lines) — UI + server for the "Colocalisation" sub-module (registry id `coloc`, group "Genetics").

**Helper functions/objects this module calls, defined elsewhere (traced and read in full for this audit):**
- `ArthOMix/global.R` — `load_default_meth_coloc_results()` (:602), `load_default_mr_harmonised()` (:568), `read_uploaded_table()` (:1186), `guess_gwas_col()` (:1181), `GWAS_COL_PATTERNS` (:1170), `gwas_col_map_ui()` (:1229), `ARTHOMIX_COLORS` (:1417), `theme_arthomix()` (:1438), `` `%||%` `` (:817).
- `ArthOMix/data_paths.R` — `METH_DATA_AVAILABLE` (:85), `METH_MR_DIR` (:94).
- `ArthOMix/R/submodules_registry.R` — module registration in `MX_MODULES` (:49).
- `ArthOMix/server.R` — generic `lapply(MX_MODULES, ...)` dispatch (:95, :248, :278) — there is no colocalization-specific line in `server.R`; the module is invoked the same way every `MX_MODULES` entry is.
- Precomputed data read by the "Preloaded Data" route: `ArthOMix/data/preloaded/methylomics/tables/script08_mendelian_randomization/tables/coloc_results.csv` (7 rows) and `mr_harmonised_all_cpgs.csv`, plus the pipeline write-up `METHODS_mendelian_randomization.md` in the same folder — all inspected directly for this audit.
- External packages called directly: `coloc` (`coloc.abf()`, `coloc.susie()`, `sensitivity()`), `TwoSampleMR` (`format_data()`, `harmonise_data()`), `susieR` (indirectly, via `coloc::coloc.susie()`, only if installed), `DT`, `ggplot2`, `data.table` (via `read_uploaded_table()`), base `stats`/`utils`.

Prepared: 2026-08-26.

This document is derived **exclusively** from the code and data files cited above. Every non-trivial technical claim carries a `file:line` citation or a direct data check. Two label conventions, matching this project's other methylomics thesis documents (e.g. `methylomics_Mendelian_Randomization.md`, `methylomics_quality_control.md`), are used throughout:

- **Scientific background** — general colocalization / statistical-genetics knowledge (textbook/literature), not a claim about this code.
- **Code evidence** — a claim about what `mod_methyl_coloc.R` (and the helpers it calls) actually does, always with a citation.

Per the audit brief, this document does **not** propose redesigning the colocalization workflow, and no code outside `mod_methyl_coloc.R` and the helper functions it calls was modified to produce it. `R/transcriptomics/mod_coloc.R` (the eQTL/eGene colocalization module, already documented in `thesis/mod_coloc_teaching_notes.md`) is a structurally similar but entirely separate module — its own comment header states the two "deliberately" do not share compute helpers (`mod_methyl_coloc.R:41-43`) — and is out of scope here except where cited for contrast.

---

## 1. Overview

**Code evidence.** The module's own header comment states its scope precisely: "Bayesian colocalisation (coloc.abf, optionally coloc.susie for multiple causal signals) between a methylation-associated genetic signal (mQTL/CpG) and a GWAS trait's signal at one genomic region: do the two association patterns share a single causal variant (PP.H4) or reflect two distinct, LD-linked variants (PP.H3)? This is the methylomics counterpart to R/transcriptomics/mod_coloc.R (eQTL colocalisation, untouched by this file) - it is NOT that module and does not implement eQTL/eGene colocalisation." (`mod_methyl_coloc.R:4-11`). The registry description echoes this: "Tests whether an mQTL and a GWAS signal at a CpG's region share a causal variant, using coloc.abf/coloc.susie. Uses the bundled GoDMC/RA-GWAS data by default, or your own uploaded summary statistics." (`mod_methyl_coloc.R:110-113`).

**Two data routes**, chosen first via a `radioButtons("data_source", ...)` on Tab 1 (`mod_methyl_coloc.R:144-148`); nothing downstream renders until one is chosen and validated:

1. **"Preloaded Data"** — reproduces `script08_mendelian_randomization`'s own completed `coloc.abf()` run (referred to in comments as `08d_mr_coloc.R`): the GoDMC cis-mQTL signal vs. the Ishigaki et al. (2022) rheumatoid-arthritis (RA) GWAS signal, at each CpG carried into MR that had ≥10 GoDMC candidate SNPs in its ±1 Mb cis window (`mod_methyl_coloc.R:14-25`). **Only the per-CpG PP.H0–H4 summary table is bundled with this deployment** — the underlying per-SNP GoDMC/RA-GWAS region data is not — so `coloc.abf()` itself cannot be re-run live for these CpGs; results are looked up, not recomputed, and SNP-level output, regional plots, and prior-sensitivity analysis are unavailable for this route.
2. **"Upload Data"** — a fully live pipeline on the user's own methylation/mQTL and GWAS summary-statistic files: column mapping, validation/harmonisation (`TwoSampleMR::format_data()`/`harmonise_data()`, the same functions `mod_methyl_mr.R`'s own upload mode uses), genomic-region/association/variant filters, live `coloc.abf()` (and, only when the user additionally supplies an LD matrix for both datasets, `coloc.susie()`), SNP-level results, regional/comparison/posterior plots, and prior/parameter sensitivity all run live (`mod_methyl_coloc.R:26-33`).

**Code evidence — stage gating.** "Nothing computes or renders ahead of an explicit click: Validate Data -> Run Colocalisation -> Generate Regional Plot / Run Sensitivity Analysis, each gated behind a has-run flag exactly like `mod_methyl_mr.R`. Changing a stage's own defining inputs invalidates every stage after it." (`mod_methyl_coloc.R:35-38`), implemented as a `reactiveValues` flag set `stage_flags` (`:303`) and a generic `invalidate_from()` helper (`:304-307`) that clears every downstream stage flag. This is verified accurate throughout §19 below.

**Where this analysis sits in the wider Methylomics workflow.** Colocalization is registered as the tenth of twelve Methylomics sub-modules, immediately after Mendelian Randomization and before Diagnostic Classifier (`R/submodules_registry.R:49`, `MX_MODULES` list, entries 9–10). Per `METHODS_mendelian_randomization.md` (§2.FF.1, §2.FF.2), colocalization was run as "the intended complement specifically for CpGs with too few instruments to support a formal pleiotropy test" in the MR analysis — i.e., in the underlying research pipeline this reproduces, colocalization is not a free-standing analysis but a haplotype-pleiotropy check that sits downstream of, and depends on, the same candidate-CpG panel, GoDMC instrument extraction, and RA-GWAS outcome data used for MR.

**Distinguishing association → overlap → colocalization → causal interpretation.** The module's own interpretation text draws this line explicitly: "Colocalisation identifies statistical compatibility with a shared causal signal - it is not, by itself, proof of biological causality. A high PP.H4 does not mean a specific SNP causes the disease; it means the methylation- and disease-associated signals in this region are consistent with arising from the same underlying genetic variant." (`.mcol_interpret()`, `mod_methyl_coloc.R:99-100`). This document adopts the same distinction throughout — see §3 and §26.

---

## 2. Scientific Purpose

**Scientific background.** Two separate genome-wide association analyses can each report a significant signal at the same genomic locus — for example, a methylation quantitative trait locus (mQTL) analysis showing that a SNP is associated with methylation at a nearby CpG, and a GWAS showing that a (possibly different) SNP in the same region is associated with a disease. Because nearby variants are correlated with one another through linkage disequilibrium (LD), two apparently overlapping signals can be produced by **two distinct causal variants** that merely sit close together on the chromosome, rather than by **one shared causal variant** that affects both traits. Colocalization analysis (Giambartolomei et al., 2014) is a Bayesian statistical framework that uses only the two traits' own regional summary statistics (effect sizes, standard errors, allele frequencies, sample sizes) — not individual-level genotypes — to estimate the posterior probability of five mutually exclusive configurations at a locus:

| Hypothesis | Meaning |
|---|---|
| H0 | No causal variant for either trait in the region |
| H1 | A causal variant for trait 1 (methylation) only |
| H2 | A causal variant for trait 2 (the GWAS trait) only |
| H3 | Two distinct causal variants, one per trait, in LD with each other |
| H4 | One shared causal variant driving both traits |

**Biological question this submodule answers.** For a given CpG whose methylation level is influenced by a nearby genetic variant (an mQTL), and a GWAS trait (rheumatoid arthritis in the Preloaded route, or any user-supplied trait in the Upload route) with an association signal in the same genomic window, does the same underlying genetic variant plausibly drive both the methylation difference and the disease-association signal (PP.H4), or are they two separate, LD-linked genetic effects that only appear to overlap (PP.H3)?

**Molecular/genetic signals compared.** Dataset 1 is always a methylation-associated genetic signal (an mQTL: SNP–CpG methylation association summary statistics, treated by the module as a quantitative trait; `mod_methyl_coloc.R:167`, `:632-636`). Dataset 2 is always a GWAS trait signal (binary case-control, e.g. RA, or a user-specified quantitative trait; `mod_methyl_coloc.R:175-176`, `:637-640`). The module never compares two GWAS traits or two mQTL/eQTL signals — the analysis is architecturally fixed to "mQTL vs. GWAS," consistent with the header comment's statement that this is not eQTL/eGene colocalization (`mod_methyl_coloc.R:9-11`).

**What the analysis can and cannot conclude.** It can conclude, with a stated posterior probability, whether the two association *patterns observed in the region* are more consistent with a single shared variant or two distinct ones. It **cannot** identify which specific SNP is causal with certainty (only assign per-SNP posterior weight — see §14), cannot establish the direction or mechanism of any causal effect, cannot correct for LD misestimation because the module's `coloc.abf()` route supplies no LD matrix at all (single-causal-variant assumption; see §11), and — per the module's own text — is "not, by itself, proof of biological causality" (`mod_methyl_coloc.R:100`).

---

## 3. Colocalization in Simple Terms

Imagine two search teams, each independently scanning the same city block for "the loudest apartment." Team A (the mQTL study) is listening for methylation-changing noise; Team B (the GWAS) is listening for disease-risk noise. Each team reports a rough probability distribution over which apartment (variant) is the source. Colocalization asks: given both teams' probability maps over the *same* small block (the genomic region), is it more likely that **one** apartment is loud for both reasons at once (H4 — shared cause), or that **two different, neighbouring** apartments are each loud for their own separate reason, and because they're neighbours their maps happen to overlap (H3 — distinct causes, correlated by proximity/LD)? Because the "signal" going into this analysis is aggregate statistics (effect size, standard error) rather than a map of who lives where, colocalization can be run without ever genotyping individual people directly — it only needs each study's own already-computed regional association numbers.

The distinction from a plain **overlap** check is the key idea to teach: naively noticing that "there's a hit near this CpG in both the mQTL study and the GWAS" only shows *association* and *positional overlap*. Colocalization goes one step further by explicitly modelling the *shared-vs-distinct-variant* alternative using the two studies' full pattern of effect sizes across every SNP in the window, not just their single best hits — which is why coloc.abf needs SNP-level summary statistics across a region, not a single index SNP p-value.

---

## 4. Source-Code Architecture

**Code evidence.** The file is organised into four blocks, in this order:

1. **Local helpers (lines 40–108)** — deliberately not shared with `R/transcriptomics/mod_coloc.R` per the code's own comment (`:41-43`): `.mcol_tip()`, `.mcol_stage_order`, the `MCOL_DEFAULT_*` constants, `.mcol_prep_ld()`, `.mcol_interpret()`, `.mcol_verdict()`.
2. **`mod_methyl_coloc_config`** (`:110-113`) — the registry entry (id, title, icon, group, description).
3. **UI functions** (`:119-293`) — `mod_methyl_coloc_ui()` and one UI-builder function per tab (`mcol_data_ui()`, `mcol_filters_controls_preloaded()`/`_upload()`, and four `uiOutput()`-only stub functions for Results/Visualisation/Sensitivity/Export, whose actual content is built server-side).
4. **`mod_methyl_coloc_server()`** (`:299-1221`) — one `moduleServer()` closure containing every reactive expression, observer, and output for all six tabs.

This all-in-one-file, all-in-one-`moduleServer()` architecture matches every other `mod_methyl_*.R` submodule in this codebase (e.g. `mod_methyl_mr.R`) and is a project-wide convention, not something specific to Colocalization.

---

## 5. Number and Names of Colocalization Tabs

**Code evidence.** `mod_methyl_coloc_ui()` defines exactly **six** tabs in one `tabsetPanel(id = ns("coloc_tabs"), type = "tabs", ...)` (`mod_methyl_coloc.R:123-131`):

| # | UI label (exact) | UI-builder function |
|---|---|---|
| 1 | `"1. Data & Setup"` | `mcol_data_ui()` (`:137-194`) |
| 2 | `"2. Filters & Parameters"` | `mcol_filters_ui()` → `uiOutput("filters_tab_body")` (`:198`) |
| 3 | `"3. Results"` | `mcol_results_ui()` → `uiOutput("results_tab_body")` (`:281`) |
| 4 | `"4. Visualisation"` | `mcol_plots_ui()` → `uiOutput("plots_tab_body")` (`:285`) |
| 5 | `"5. Sensitivity Analysis"` | `mcol_sensitivity_ui()` → `uiOutput("sensitivity_tab_body")` (`:289`) |
| 6 | `"6. Export"` | `mcol_export_ui()` → `uiOutput("export_tab_body")` (`:293`) |

No nested `tabsetPanel` exists inside any of these six tabs (unlike, e.g., `mod_methyl_mr.R`'s Plots tab, which nests a further six plot-type tabs). Tabs 2–6 are entirely server-rendered (`uiOutput`) and dynamically switch content based on `input$data_source` and `stage_flags`; Tab 1 mixes static UI with two `conditionalPanel()`s gated on `input$data_source` (`:150,163`).

---

## 6. Tab-by-Tab Analysis

### 6.1 "1. Data & Setup"

**Purpose.** Choose a data route (Preloaded vs. Upload) and, for Upload, supply and column-map the two summary-statistic files.

**Scientific question addressed.** Which methylation/mQTL signal and which GWAS trait signal are being compared, and are their raw files well-formed enough to proceed?

**Preloaded-route UI (`:150-162`).** If `METH_DATA_AVAILABLE` is `TRUE`, shows a CpG multi-select (`output$pre_cpg_ui`, `:319-328`, populated from `unique(pre_coloc_tbl()$cpg)`) and a note that only the completed run's PP.H0–H4 summary is bundled. If `METH_DATA_AVAILABLE` is `FALSE`, shows a warning box and no CpG picker (`:159-160`).

**Upload-route UI (`:163-184`).** Two file boxes:
- *Dataset 1 (methylation/mQTL)*: `fileInput("meth_file")` (CSV/TSV/TXT), then a shared column-mapper (`gwas_col_map_ui()`, with `extra_fields = "n"` for sample size) rendered into `meth_map_ui`, plus a module-specific `meth_extra_map_ui` for CpG ID, optional SNP chromosome/position, optional gene annotation, and (if the file contains >1 CpG) a target-CpG selector.
- *Dataset 2 (GWAS)*: a free-text trait label, a binary/quantitative radio (`gwas_type`), a case-fraction slider shown only for binary traits, `fileInput("gwas_file")`, and the same shared column-mapper.

**Required columns (Upload route).** SNP ID, beta, SE, p-value, effect allele, other allele are required for both files (`build_validate_state_upload()`, `:459-461`, `req()`); EAF is optional but functionally required later for a quantitative methylation trait (see §8). CpG ID is required for the methylation file (`:466-467`).

**Validation requirements.** See §8.

**Default values.** `gwas_type = "cc"` (`:176`), `case_frac = 0.33` (`:179`).

**Functions called.** `read_uploaded_table()`, `gwas_col_map_ui()`, `guess_gwas_col()`, `req()`, `renderUI()`.

**Output.** A live preview `DT::dataTableOutput` of the first 20 rows of each uploaded file (Upload route) or of the selected CpGs' `(cpg, nsnps)` pair (Preloaded route) — explicitly "descriptive only" (`:405`) and never the PP.H0–H4 results themselves, plus a "Validate Data" action button (`:191`).

**Relationship with other tabs.** Defines everything downstream; any change to a defining input here invalidates the `"validate"` stage and everything after it (`:376-378`).

### 6.2 "2. Filters & Parameters"

**Purpose.** Set the region/association/variant filters and coloc priors (Upload route), or a minimum-instrument-count filter (Preloaded route, since priors/window are fixed upstream).

**Scientific question addressed.** Which subset of variants and which prior assumptions should feed into `coloc.abf()`/`coloc.susie()`?

**Body depends on route** (`output$filters_tab_body`, `:574-579`, `req(stage_flags$validate)` — the tab renders nothing until Tab 1's Validate Data has succeeded):

- **Preloaded** (`mcol_filters_controls_preloaded()`, `:200-209`): one `numericInput("pre_min_nsnps", value = 0)` and the "Run Colocalisation" button. A note states priors (p1=1e-4, p2=1e-4, p12=1e-5) and the ±1 Mb window were fixed upstream and are not re-run live.
- **Upload** (`mcol_filters_controls_upload()`, `:211-277`): see the full input table in §6.6 below — genomic-window restriction, minimum shared SNPs, mQTL/GWAS p-value thresholds, MAF bounds, minimum sample size, duplicate/ambiguous-SNP handling, p1/p2/p12 priors, an optional `coloc.susie()` toggle with LD file uploads, credible-set coverage/max-iterations, and (in a collapsible "Advanced parameters" `<details>`) the posterior-probability decision threshold.

**Validation requirements.** None directly on this tab; all numeric bounds are enforced downstream when `build_run_state_upload()`/`build_run_state_preloaded()` run (§14).

**Functions called.** `uiOutput()`, `numericInput()`, `checkboxInput()`, `radioButtons()`, `fileInput()`, `conditionalPanel()`.

**Relationship with other tabs.** Any change to a filter/prior/method input invalidates the `"run"` stage (`:581-586`), forcing Tab 3 back to its "not run yet" placeholder.

### 6.3 "3. Results"

**Purpose.** Display the completed colocalization result: PP.H0–H4, verdict, interpretation, and (Upload route only) SNP-level and (optionally) coloc.susie results.

**Scientific question addressed.** What did the analysis conclude, and how confident is that conclusion?

**Gating.** `output$results_tab_body` (`:749-826`) returns an "Not run yet…" placeholder unless `stage_flags$run` is `TRUE`.

**Preloaded content (`:756-780`).** A summary box (GWAS trait, CpGs tested, total SNPs, method = "coloc.abf (preloaded)"), a genome-build/prior note (GRCh37; fixed priors), a "Focus CpG" selector feeding a per-CpG interpretation panel (`output$pre_focus_ui`, `:828-847`), and a per-CpG results table with a `verdict` column computed client-side by `mapply(.mcol_verdict, tbl$PP.H3, tbl$PP.H4)`.

**Upload content (`:781-824`).** A summary box (GWAS trait, CpG, genomic region if available, variants analysed, shared SNPs used, method — "coloc.abf" or "coloc.abf + coloc.susie", sample sizes, priors), the lead variant (highest `SNP.PP.H4`), the five PP.H0–H4 values, a verdict line, an "Interpretation" box (`.mcol_interpret()`), a conditional coloc.susie credible-set table if that method ran, and a full SNP-level results table with a CSV download.

**Functions called.** `mapply()`, `.mcol_verdict()`, `.mcol_interpret()`, `DT::renderDataTable()`, `DT::formatRound()`/`formatSignif()`.

**Relationship with other tabs.** Feeds Tab 4 (plots use the same `run_state()`), Tab 5 (sensitivity re-uses `run_state()$h_full`/priors), and Tab 6 (every download reads `run_state()`).

### 6.4 "4. Visualisation"

**Purpose.** Generate plots from the completed run.

**Scientific question addressed.** How does the shared/distinct-signal evidence look spatially and per-SNP, not just as five summary numbers?

**Gating.** `output$plots_tab_body` (`:932-965`) requires `stage_flags$run`; nothing plots until the user additionally clicks "Generate Regional Plot" (`input$plot_btn`), which only sets `stage_flags$plot <- TRUE` (`:967`) — a client-side `conditionalPanel(condition = "input[...] > 0")` gate, the same idiom the code comment says mirrors `R/transcriptomics/mod_coloc.R` (`:875-876`).

**Preloaded route.** Only a bar plot of PP.H0–H4 for one selected CpG (`pre_pp_plot`), because "no per-SNP GoDMC/RA-GWAS region data is bundled to build a regional or posterior plot from" (`:939`).

**Upload route.** Three plots (see §17): Regional association (only if a SNP position column was mapped), Comparison of signals, and Posterior support by variant.

**Relationship with other tabs.** Purely a rendering layer over `run_state()`; produces no state anything else reads, except that its download handlers (Tab 6) call the same `build_*_plot()` functions again at export time.

### 6.5 "5. Sensitivity Analysis"

**Purpose.** Test how robust the H3/H4 conclusion is to the chosen p12 prior, and how much the shared-SNP count and PP.H3/H4 shift under alternative filter settings.

**Scientific question addressed.** Is the coloc verdict a stable feature of the data, or an artefact of one specific prior/filter choice?

**Gating.** `output$sensitivity_tab_body` (`:1007-1042`) requires `stage_flags$run`, and explicitly returns an early placeholder for the Preloaded route: "Sensitivity re-analysis needs per-SNP Bayes factors, which aren't bundled with the preloaded results… Use Upload Data…" (`:1010-1012`) — this tab performs **no compute at all** on the Preloaded route.

**Upload route — two independent sub-analyses, each behind its own button:**
1. *Prior sensitivity* (`sens_prior_btn`) — calls `coloc::sensitivity()` on the already-computed `abf_res` object across a grid of p12 values, plotted as PP.H3/PP.H4 vs. p12 (log scale).
2. *Parameter sensitivity* (`sens_param_btn`) — re-runs `coloc.abf()` from scratch (via the local `.mcol_rerun_with()` helper) at a fixed baseline plus perturbations of window size (×0.5/2/4), minimum MAF (0/0.01/0.05), and mQTL p-value threshold (1/0.05/5e-8), reporting `nsnps`/PP.H3/PP.H4 for each.

**Relationship with other tabs.** Reads `run_state()` (for `h_full`, priors, gwas_type) and `validate_state()` (for `case_frac`); writes only to its own local `sensitivity_state()`, read nowhere else.

### 6.6 Filters & Parameters — full input inventory (Upload route)

| Input ID | Type | Purpose | Default | Allowed values | Used by |
|---|---|---|---|---|---|
| `use_window` | checkbox | Restrict analysis to a genomic window | `FALSE` | on/off | `build_run_state_upload()` region filter |
| `window_center` | radio | Window centred on lead SNP or a manual position | `"lead"` | `lead`/`manual` | window filter |
| `window_pos` | numeric | Manual centre position (bp) | `NA` | any bp | window filter (only if `window_center="manual"`) |
| `window_kb` | numeric | Window half-width (kb) | `1000` (`MCOL_DEFAULT_WINDOW_KB`) | ≥1 | window filter |
| `f_min_shared` | numeric | Minimum shared SNPs required | `10` (`MCOL_DEFAULT_MIN_SHARED_SNPS`) | ≥3 | hard gate before `coloc.abf()` |
| `f_pval_meth` | numeric | mQTL p-value threshold | `1` (no filter) | [0,1] | association filter |
| `f_pval_gwas` | numeric | GWAS p-value threshold | `1` (no filter) | [0,1] | association filter |
| `f_maf_min` / `f_maf_max` | numeric | MAF bounds | `0` / `0.5` | [0,0.5] | variant filter |
| `f_min_n` | numeric | Minimum sample size (either dataset) | `0` | ≥0 | variant filter |
| `f_dedup` | checkbox | Remove duplicated SNPs | `TRUE` | on/off | dedup filter |
| `remove_ambiguous` | checkbox | Drop palindromic SNPs instead of EAF-resolving them | `FALSE` | on/off | selects `harmonise_data(action=3)` vs. default `action=2` |
| `p1` | numeric | Prior P(SNP associated with trait 1 only) | `1e-4` | [0,1] | `coloc.abf()`/`coloc.susie()` |
| `p2` | numeric | Prior P(SNP associated with trait 2 only) | `1e-4` | [0,1] | `coloc.abf()`/`coloc.susie()` |
| `p12` | numeric | Prior P(SNP associated with both traits) | `1e-5` | [0,1] | `coloc.abf()` only — **not** passed to `coloc.susie()` (see §26, finding A-1) |
| `use_susie` | checkbox | Also run multi-signal colocalization | `FALSE` | on/off | gates the `coloc.susie()` branch |
| `ld1_file` / `ld2_file` | file | LD (SNP×SNP correlation) matrices for dataset 1/2 | none | CSV/TSV | required for `coloc.susie()`; no LD is ever inferred |
| `susie_coverage` | numeric | SuSiE credible-set coverage | `0.95` | [0.5,0.99] | `coloc.susie()` |
| `susie_maxit` | numeric | SuSiE max iterations | `100` | ≥10 | `coloc.susie()` |
| `pp_threshold` | numeric | Posterior threshold for supported/refuted verdict | `0.8` (`MCOL_DEFAULT_PP_THRESHOLD`) | [0.5,0.99] | `.mcol_verdict()` |

---

## 7. Input Data

**File types and delimiters.** CSV/TSV/TXT, read via `data.table::fread()` inside `read_uploaded_table()` (`global.R:1186`), which auto-detects the delimiter (`fread`'s own default behaviour) and wraps the call in `tryCatch()`, returning `NULL` on failure rather than raising.

**Dataset 1 — Methylation/mQTL file.** UI text: "One row per SNP, optionally x CpG (long format). CSV/TSV. Methylation is treated as a quantitative trait - map the effect-allele-frequency column below (coloc.abf needs it to analyse a quantitative trait without a directly-known phenotype SD)." (`mod_methyl_coloc.R:167`). Required columns (mapped via `gwas_col_map_ui(extra_fields="n")` plus module-specific selectors): SNP ID, beta, SE, p-value, effect allele, other allele, sample size (N), CpG ID; optional: EAF, SNP chromosome, SNP position, gene annotation. Multiple CpGs in one file are supported — the user picks one target CpG per run (`meth_target_cpg_ui`, `:358-373`); "colocalisation is run for one CpG's region at a time" (`:371`).

**Dataset 2 — GWAS file.** Required columns: SNP ID, beta, SE, p-value, effect allele, other allele, sample size (N); optional: EAF. Trait type (binary/quantitative) and, for binary, the case fraction, are supplied by the user directly (not read from the file) because the module never asks for per-arm case/control counts, only a proportion (`:174-180`).

**Genome build.** The Preloaded route states GRCh37 explicitly (`mod_methyl_coloc.R:769`, matching `METHODS_mendelian_randomization.md`'s statement that CpG positions and both mQTL/outcome GWAS builds are hg19/GRCh37). The Upload route's genome build is **whatever the user's own files use** — the module performs no build check or liftover; see §23 (finding B-4).

**What enters → biological meaning → transformation → what leaves, per format:**

| Stage | What enters | Biological meaning | Transformation | What leaves |
|---|---|---|---|---|
| Raw meth file | delimited table, 1 row/SNP(×CpG) | candidate cis-mQTL association test results | `read_uploaded_table()`, column mapping | `data.frame`, unmapped |
| Raw GWAS file | delimited table, 1 row/SNP | GWAS summary statistics for a trait | `read_uploaded_table()`, column mapping | `data.frame`, unmapped |
| `TwoSampleMR::format_data()` (×2) | mapped data.frames + column-name args | — | standardises to `TwoSampleMR`'s canonical `beta.exposure`/`se.exposure`/... schema; assigns `exposure`/`outcome`/`id` labels | two `format_data()`-schema data.frames |
| `TwoSampleMR::harmonise_data()` | the two formatted data.frames | aligning effect alleles between two independently-generated association files | allele flipping, ambiguous/palindromic resolution or removal, `mr_keep` flagging | one harmonised data.frame, filtered to `mr_keep==TRUE` |
| Region/association/variant filters | harmonised data.frame | restrict analysis to the region and quality bar the user set | window/p-value/MAF/N/dedup subsetting (plain `data.frame` row filters, `:609-623`) | filtered data.frame `h` |
| `coloc::coloc.abf()` input lists (`d1`,`d2`) | filtered `h` | per-SNP effect estimates for two traits at one locus | reshaped into `coloc`'s named-list schema (`beta`,`varbeta`,`N`,`type`,`snp`,`MAF`/`s`) | `abf_res` (posterior probabilities + per-SNP Bayes factors) |

**Preloaded-route data provenance (no live file upload).** `load_default_meth_coloc_results()` reads `coloc_results.csv` — one row per CpG, columns `cpg, nsnps, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4, verdict` (verified directly: 7 data rows, header confirmed by `cat`). This is the **already-computed output** of an upstream script, not raw summary statistics — the module reads a finished analysis table, it does not run `coloc.abf()` for this route at all (`build_run_state_preloaded()`, `mod_methyl_coloc.R:593-602`, contains no call to `coloc::coloc.abf`).

---

## 8. Input Validation

Per `build_validate_state_upload()` (`:458-524`) and `build_validate_state_preloaded()` (`:436-456`):

| Checked | How | On failure |
|---|---|---|
| Files present, core columns mapped | `req(input$meth_file, input$gwas_file, input$meth_cpg, input$meth_snp, ...)` (`:459-461`) | Validate button silently produces nothing (Shiny `req()` semantics) until all mapped |
| File readable as a table | `read_uploaded_table()` returns non-`NULL` | `validate(need(...))` — "Could not read the … file." |
| Rows exist for the selected target CpG | `nrow(meth_sub) > 0` | "No rows for the selected CpG…" |
| SNP-ID overlap between files | `intersect(meth_sub[[snp]], gwas_raw[[snp]])` non-empty, **before** formatting | "No overlapping SNP IDs… check that both use the same SNP identifier convention (e.g. rsIDs)." |
| `format_data()` succeeds and yields rows | `tryCatch(...)`, `nrow(...) > 0` | "No usable rows after formatting… check column mapping." (per file) |
| `harmonise_data()` succeeds and yields rows | `tryCatch(...)`, `nrow(...) > 0` | "Harmonisation found no overlapping, alignable SNPs… check that both use the same SNP identifiers (rsIDs) and that allele columns are mapped correctly." |
| At least one variant survives `mr_keep==TRUE` | post-filter `nrow(harmonised) > 0` | "No variants survived harmonisation (all ambiguous/unresolvable, e.g. unresolved palindromic SNPs)." |
| ≥1 CpG selected (Preloaded) | `length(cpgs) > 0` | "Pick at least one CpG." |
| A coloc result exists for the selection (Preloaded) | `nrow(tbl) > 0` after filtering `pre_coloc_tbl()` | "No preloaded coloc result for the selected CpG(s)." |

At the **Run Colocalisation** stage (`build_run_state_upload()`), additional checks:

| Checked | How | On failure |
|---|---|---|
| Minimum shared SNPs after filtering | `nrow(h) >= min_shared` | "Fewer than N shared SNPs remain after filtering… colocalisation needs a minimally informative set of SNPs across the region. Loosen the filters or widen the genomic window." |
| Non-zero, non-negative SEs | `all(h$se.exposure > 0) && all(h$se.outcome > 0)` | "Some retained variants have a zero or negative standard error - cannot compute a Bayes factor for them." |
| MAF available for the (quantitative) mQTL dataset | `!is.null(d1$MAF)` | "coloc.abf needs a minor allele frequency to analyse a quantitative trait (methylation) without a directly-supplied phenotype SD - map an effect-allele-frequency column… and re-validate." |
| `coloc.abf()` itself does not error | `tryCatch()` around the call | "coloc.abf() failed: …" (raw package error message) |
| LD matrices align to the SNP set (coloc.susie only) | `.mcol_prep_ld()`, `length(common_ld) >= 3` | descriptive `susie_note`, method skipped (no crash) |
| `coloc.susie()` does not error / does return a credible set | `tryCatch()`; `nrow(summary)==0` check | descriptive `susie_note`, method skipped |

**What is explicitly *not* validated** (identified as audit findings, not invented deficiencies — see §23/§26):
- **Genome build consistency** between the uploaded mQTL file and the uploaded GWAS file, or against the ChAMPdata/hg19 build implicitly assumed by the CpG-position columns, is never checked.
- **Invalid p-values** (outside [0,1]) are *counted* (`invalid_p1`/`invalid_p2`, `:481-482`) and shown in the validation summary table but are **not filtered out or blocked** — rows with an out-of-range p-value still proceed into `format_data()`/`harmonise_data()`/`coloc.abf()`.
- **Duplicate SNPs** are counted at validate time (`dup1`, `dup2`, `:475-476`) but only removed later, at run time, and only if the user leaves `f_dedup = TRUE` (its default) — a user who unchecks it keeps duplicates in the coloc input.
- **Sample-size sanity** (e.g., a mapped N column containing zeros or absurd values) is never range-checked; only `NA` values are treated specially in the run-stage minimum-N filter (`:622-623`).
- **Strand orientation** beyond what `TwoSampleMR::harmonise_data()` itself performs is not independently checked by this module.

---

## 9. Complete Function Inventory

### 9.1 Application/custom functions (defined in `mod_methyl_coloc.R` unless noted)

| Function | Defined at | Role |
|---|---|---|
| `.mcol_tip()` | `:45` | Small info-icon tooltip UI helper |
| `.mcol_prep_ld()` | `:67-76` | Aligns an uploaded LD matrix to a target SNP set |
| `.mcol_interpret()` | `:81-102` | Builds the cautious-language interpretation panel from five PPs |
| `.mcol_verdict()` | `:104-108` | Classifies a result as supported/refuted/inconclusive from PP.H3/H4 + threshold |
| `mod_methyl_coloc_ui()` | `:119-133` | Top-level module UI (6-tab `tabsetPanel`) |
| `mcol_data_ui()` … `mcol_export_ui()` | `:137-293` | Per-tab UI builders |
| `mcol_filters_controls_preloaded()` / `_upload()` | `:200-277` | Route-specific Filters & Parameters bodies |
| `mod_methyl_coloc_server()` | `:299-1221` | Entire server logic |
| `build_validate_state_preloaded()` / `_upload()` | `:436-524` | Validate-stage business logic |
| `build_run_state_preloaded()` / `_upload()` | `:593-708` | Run-stage business logic (calls `coloc.abf`/`coloc.susie`) |
| `.mcol_rerun_with()` | `:984-1005` | Re-runs `coloc.abf()` under one perturbed filter, for parameter sensitivity |
| `build_pp_bar_plot()` | `:878-886` | PP.H0–H4 bar chart |
| `build_region_plot()` | `:888-902` | Regional -log10(p) plot for both traits |
| `build_comparison_plot()` | `:904-911` | mQTL vs. GWAS -log10(p) scatter, coloured by SNP.PP.H4 |
| `build_posterior_plot()` | `:913-930` | Per-SNP posterior-support plot |
| `make_plot_dl()` | `:1213-1216` | Generic `downloadHandler` factory for the three plots |
| `invalidate_from()` | `:304-307` | Clears `stage_flags` from a given stage onward |
| `load_default_meth_coloc_results()` | `global.R:602-607` | Reads `coloc_results.csv` |
| `load_default_mr_harmonised()` | `global.R:568-573` | Reads `mr_harmonised_all_cpgs.csv` (used for validation-summary context only) |
| `gwas_col_map_ui()` | `global.R:1229-1259` | Shared column-mapping UI (also used by `mod_methyl_mr.R`) |
| `guess_gwas_col()` | `global.R:1181-1184` | Regex-based column-name auto-detection |
| `read_uploaded_table()` | `global.R:1186` | Safe `fread()` wrapper |
| `theme_arthomix()` | `global.R:1438-1454` | Shared ggplot theme |

### 9.2 R / package functions actually called

| Function | Package | Used for |
|---|---|---|
| `coloc::coloc.abf()` | coloc 5.2.3 | Bayesian colocalization, single-causal-variant model |
| `coloc::coloc.susie()` | coloc 5.2.3 (requires `susieR`) | Multi-signal colocalization via SuSiE fine-mapping |
| `coloc::sensitivity()` | coloc 5.2.3 | Prior-sensitivity re-weighting of an existing `coloc.abf()` result |
| `TwoSampleMR::format_data()` | TwoSampleMR 0.7.8 | Standardise raw columns to exposure/outcome schema |
| `TwoSampleMR::harmonise_data()` | TwoSampleMR 0.7.8 | Allele harmonisation between exposure and outcome |
| `data.table::fread()` | data.table 1.18.2.1 | Fast delimited-file reading (via `read_uploaded_table()`) |
| `reactive()`, `reactiveVal()`, `reactiveValues()` | shiny | Reactive state |
| `observeEvent()` | shiny | Button/input-change handlers |
| `req()` | shiny | Silent gating on missing inputs |
| `validate()`, `need()` | shiny | User-facing error messages |
| `renderUI()`, `uiOutput()` | shiny | Dynamic UI |
| `DT::renderDataTable()`, `DT::dataTableOutput()`, `DT::datatable()`, `DT::formatRound()`, `DT::formatSignif()` | DT | Interactive tables |
| `renderPlot()`, `plotOutput()` | shiny | Static plots |
| `downloadHandler()`, `downloadButton()` | shiny | CSV/plot export |
| `updateTabsetPanel()` | shiny | Auto-switch to Results after a run |
| `showNotification()` | shiny | Toast messages on success/failure |
| `withSpinner()` | shinycssloaders | Loading spinner around plot outputs |
| `ggplot()`, `geom_col()`, `geom_point()`, `geom_vline()`, `geom_line()`, `facet_wrap()`, `scale_fill_manual()`, `scale_color_manual()`, `scale_color_gradient()`, `scale_x_log10()`, `coord_flip()`, `labs()` | ggplot2 | Plot construction |
| `ggplot2::ggsave()` | ggplot2 | Plot-file export |
| `write.csv()` | base utils | CSV export |
| `tryCatch()` | base | Error containment around every fallible call |
| `do.call()` | base | Dynamic-argument dispatch for `format_data()` and `.mcol_rerun_with()` |
| `mapply()` | base | Vectorised `.mcol_verdict()` over a table of CpGs |
| `stats::complete.cases()`, `stats::median()`, `stats::setNames()` | stats | Missing-value counting, sample-size aggregation, named-vector construction |
| `utils::head()`, `utils::modifyList()` | utils | Preview truncation; sensitivity-grid argument overrides |
| `duplicated()`, `intersect()`, `match()`, `order()`, `round()`, `sprintf()`, `pmin()`, `pmax()` | base | Row/column bookkeeping and formatting throughout |
| `requireNamespace()` | base | Optional-dependency check for `susieR` |

No `merge()`, `left_join()`/`inner_join()` (dplyr), or `subset()`/`filter()`/`select()`/`mutate()` (dplyr) calls exist anywhere in this file — all data-frame joins/filters are done with base-R bracket indexing (`h[condition, , drop = FALSE]`) and `match()`/`intersect()`. This is a deliberate, consistent choice across the file, not an omission.

---

## 10. Function-by-Function Explanation

### `coloc::coloc.abf()`

**Source:** package function (coloc 5.2.3).
**Where used:** `build_run_state_upload()` (`:642`) and `.mcol_rerun_with()` (`:1002`, the parameter-sensitivity re-run).
**Purpose:** Computes the five posterior probabilities (PP.H0–H4) and per-SNP log approximate Bayes factors for a pair of single-trait association datasets at one locus, under the single-causal-variant assumption.
**Input:** Two named lists (`dataset1`, `dataset2`), each with `beta`, `varbeta` (SE²), `N`, `type` (`"quant"`/`"cc"`), `snp`, plus `MAF` (quantitative) or `s` (case fraction, for `"cc"`); and scalar priors `p1`, `p2`, `p12`.
**Processing:** Approximates a Bayes factor per SNP per trait from `beta`/`varbeta`/`MAF` (or `s`), combines across SNPs under the specified priors to compute the five joint posterior probabilities, and records each SNP's own contribution (`internal.sum.lABF`, `SNP.PP.H4`).
**Output:** A list with `$summary` (named vector `nsnps, PP.H0.abf…PP.H4.abf`) and `$results` (per-SNP data.frame).
**Role in Colocalisation:** This is the module's core statistical engine for the Upload route.
**Scientific meaning:** Directly answers "shared vs. distinct causal variant" for this SNP set.
**Audit assessment:** Correctly implemented — called through a `tryCatch()`/`validate(need(...))` pair (`:642-643`) so a package-level failure surfaces as a readable message rather than crashing the app.
**Potential issue:** None found beyond the input-preparation issues noted in §23 (MAF/sdY handling asymmetry between the two datasets).

### `coloc::coloc.susie()`

**Source:** package function (coloc 5.2.3; requires `susieR`, checked via `requireNamespace()`, `:663`).
**Where used:** `build_run_state_upload()` (`:682-686`), only when `input$use_susie` is checked and both LD files are supplied and align.
**Purpose:** Relaxes the single-causal-variant assumption by first running SuSiE fine-mapping (using the supplied LD matrices) on each dataset independently to identify credible sets, then testing colocalization between every pair of credible sets across the two datasets.
**Input:** Two named lists like `coloc.abf()`'s but each additionally carrying an `LD` matrix; `susie.args = list(maxit, coverage)`; `p1`, `p2`, and **`p12 = MCOL_DEFAULT_P12_SUSIE` (hardcoded 5e-6)**, not the user's `p12` slider value.
**Processing:** SuSiE fine-mapping per dataset, then pairwise colocalization testing between credible sets.
**Output:** A list with `$summary` (one row per credible-set pair, each with its own PP.H0–H4).
**Role in Colocalisation:** Optional multiple-signal extension, only reachable with user-supplied LD.
**Scientific meaning:** Distinguishes multiple independent causal signals per locus, which `coloc.abf()` cannot.
**Audit assessment:** Correctly gated (no LD is ever inferred, estimated, or fabricated — `:65-66`, `:256`) and correctly guarded against a missing package or a failed/empty run (`:663-693`). **The hardcoded `p12` is a genuine deviation from the user-facing Priors panel and is flagged as finding A-1 in §26.**
**Potential issue:** See §26, A-1.

### `coloc::sensitivity()`

**Source:** package function (coloc 5.2.3).
**Where used:** `observeEvent(input$sens_prior_btn, ...)` (`:1044-1058`).
**Purpose:** Re-derives PP.H3/PP.H4 from the already-computed `abf_res` Bayes factors across a grid of alternative `p12` values, without re-reading data or re-filtering variants.
**Input:** The `abf_res` object from `coloc.abf()`, a `rule` string (e.g. `"H4 > 0.8"`), `npoints`.
**Processing:** Re-weights the existing per-SNP Bayes factors at each grid point's `p12`.
**Output:** A data.frame with `p12`, `PP.H3.abf`, `PP.H4.abf` (and others) per grid point.
**Role in Colocalisation:** Prior-robustness check, Tab 5.
**Scientific meaning:** Shows whether the H3/H4 conclusion depends sensitively on the somewhat-arbitrary `p12` prior.
**Audit assessment:** Correctly implemented; `doplot = FALSE` is passed so the package's own base-R plot is suppressed and the module's own ggplot is used instead (`:1047-1050`, `:1090-1105`).
**Potential issue:** None found — errors are caught and degrade to `sens <- NULL` (`:1053`), silently leaving the Prior sensitivity panel empty rather than crashing; this silent degradation (no `showNotification` on failure here, unlike the Validate/Run stages) is a minor UX gap, not a scientific correctness issue (see §26, finding C-2).

### `TwoSampleMR::format_data()`

**Source:** package function (TwoSampleMR 0.7.8).
**Where used:** `build_validate_state_upload()`, twice (`:486-501`, once per dataset).
**Purpose:** Standardises arbitrary user column names into TwoSampleMR's canonical schema (`beta.exposure`, `se.exposure`, `pval.exposure`, `effect_allele.exposure`, `other_allele.exposure`, `eaf.exposure`, `samplesize.exposure`, and the outcome-side equivalents).
**Input:** Raw data.frame + column-name arguments built dynamically via `do.call()` from the user's mapped `input$meth_*`/`input$gwas_*` selections; `chr_col`/`pos_col` supplied only if the user mapped them (`:490-491`).
**Processing:** Renames/recasts columns; drops rows unusable for MR.
**Output:** A `format_data()`-schema data.frame with a `type = "exposure"` or `"outcome"` label.
**Role in Colocalisation:** Prerequisite standardisation step before `harmonise_data()`.
**Scientific meaning:** Puts the mQTL and GWAS association statistics into a common, machine-readable shape for downstream allele harmonisation.
**Audit assessment:** Correctly implemented and wrapped in `tryCatch()`.
**Potential issue:** None found.

### `TwoSampleMR::harmonise_data()`

**Source:** package function (TwoSampleMR 0.7.8).
**Where used:** `build_validate_state_upload()` (`:506`).
**Purpose:** Aligns the effect allele/direction between the two independently-formatted datasets, flags/removes palindromic (A/T, G/C) and otherwise ambiguous or strand-inconsistent SNPs.
**Input:** The two `format_data()` outputs; `action = 2` (default; infer palindromic strand via EAF) or `action = 3` (drop all palindromic SNPs) depending on `input$remove_ambiguous` (`:505`).
**Processing:** Merges the two datasets by SNP ID, flips beta signs/alleles where needed, sets `mr_keep` and `palindromic`/`ambiguous` flags.
**Output:** One combined data.frame, one row per candidate SNP, with `mr_keep` marking which rows are safe to analyse.
**Role in Colocalisation:** This is the step that actually produces the shared, allele-aligned SNP set `coloc.abf()` will consume.
**Scientific meaning:** Without this step, a beta reported against the reference allele in one file and the alternate allele in the other would silently corrupt the comparison — this is the module's core allele-harmonisation safeguard.
**Audit assessment:** Correctly implemented, reuses the exact same function/mode `mod_methyl_mr.R`'s own upload route uses (per the file's own header comment, `:28-29`), and the module explicitly restricts to `mr_keep == TRUE` rows before analysis (`:518`).
**Potential issue:** None found in the call itself; see §23 for a data-provenance caveat specific to the Preloaded route's *validation-summary display*, not this function.

### `.mcol_verdict()`

**Source:** custom, `mod_methyl_coloc.R:104-108`.
**Where used:** Results tables (both routes), interpretation panels, CSV exports.
**Purpose:** Converts (PP.H3, PP.H4, threshold) into one of three text verdicts.
**Input:** `h3`, `h4` (numeric probabilities), `threshold` (default `MCOL_DEFAULT_PP_THRESHOLD = 0.8`).
**Processing:** `if (h4 >= threshold) "coloc-supported" else if (h3 >= threshold) "coloc-refuted" else "inconclusive"`.
**Output:** A character string.
**Role in Colocalisation:** Directly reproduces the upstream pipeline's own decision rule ("A CpG-locus pair was classified as coloc-supported when… PP.H4 reached ≥ 0.8, coloc-refuted when… PP.H3 reached ≥ 0.8, and inconclusive otherwise," `METHODS_mendelian_randomization.md` §2.FF.2).
**Scientific meaning:** A binary/ternary decision layered on top of continuous posterior probabilities, for reporting convenience.
**Audit assessment:** Correctly implemented and consistent with the upstream pipeline's own threshold logic. **Note:** the threshold is user-adjustable in the Upload route (`pp_threshold` input) but fixed at 0.8 in the Preloaded route (`.mcol_verdict()` called with its default argument at `:758`, `:853`, with no UI control to change it for that route) — this is by design (priors/threshold were "fixed by the upstream pipeline," `:769`), not an oversight.
**Potential issue:** None found.

### `.mcol_interpret()`

**Source:** custom, `mod_methyl_coloc.R:81-102`.
**Where used:** Results tab, both routes.
**Purpose:** Produces the cautious-language narrative panel from five posterior probabilities.
**Input:** `h0, h1, h2, h3, h4`.
**Processing:** Identifies the `which.max()` hypothesis, selects a `switch()`-based headline sentence conditioned on whether that hypothesis clears `MCOL_DEFAULT_PP_THRESHOLD`, and always appends a fixed causality caveat.
**Output:** A `tagList()` of `<p>` elements.
**Role in Colocalisation:** The module's primary defence against a user over-interpreting a coloc result as proof of causality (see §2, §26).
**Audit assessment:** Correctly implemented; language is appropriately hedged ("Moderate evidence," "not strongly resolved," "not, by itself, proof of biological causality").
**Potential issue:** None found.

### `.mcol_rerun_with()`

**Source:** custom, `mod_methyl_coloc.R:984-1005`.
**Where used:** Parameter-sensitivity handler (`:1060-1088`).
**Purpose:** Applies one perturbed filter set to the *original* harmonised data (`h_full`, i.e., pre-window/pre-association/pre-MAF filtering) and re-runs `coloc.abf()`.
**Input:** `h_full`, `gwas_type`, `case_frac`, `p1`, `p2`, `p12`, and the filter values to test.
**Processing:** Re-applies window/p-value/MAF/N filtering, rebuilds `d1`/`d2`, calls `coloc.abf()`, returns `NULL` on any early-exit condition (fewer than 6 SNPs, missing MAF, or a `coloc.abf()` error) rather than propagating an error.
**Output:** A named numeric vector `c(nsnps, PP.H3, PP.H4)`, or `NULL`.
**Role in Colocalisation:** Powers the parameter-sensitivity grid table.
**Audit assessment:** Correctly implemented; the `nrow(h) < 6` early-return (`:996`) is a stricter internal floor than the user-facing `f_min_shared` default of 10, applied silently (a perturbed row that fails this floor is simply dropped from the sensitivity table with no explanatory message — see §26, finding C-3).
**Potential issue:** See §26, C-3.

### `.mcol_prep_ld()`

**Source:** custom, `mod_methyl_coloc.R:67-76`.
**Purpose:** Aligns an uploaded, user-supplied LD (SNP×SNP correlation) matrix to a target SNP ID set.
**Input:** Raw uploaded table (first column = SNP IDs, header row = SNP IDs) and a vector of SNP IDs to match against.
**Processing:** Coerces to a numeric matrix, sets row names from the first column, intersects row/column names with the target SNP set, and requires ≥2 common SNPs.
**Output:** A square, subsetted numeric matrix, or `NULL`.
**Audit assessment:** Correctly implemented, with an explicit "never a guessed/simulated matrix" comment (`:65-66`) reflected faithfully in the code (returns `NULL`, never a fabricated identity/default matrix).
**Potential issue:** None found.

---

## 11. Statistical Method

**Exact implementation.** The module implements **Approximate Bayes Factor colocalization** (`coloc::coloc.abf()`, Giambartolomei et al. 2014) as its primary and only method for the Preloaded route, with an optional **SuSiE-based multi-signal colocalization** (`coloc::coloc.susie()`, Wallace 2021) for the Upload route when LD data is supplied.

**Required assumptions (coloc.abf):**
- Exactly one causal variant per trait in the tested region (the assumption `coloc.susie()` relaxes, at the cost of requiring LD).
- Summary statistics (beta, SE, and either MAF or a known phenotype SD/case fraction) are sufficient to approximate each SNP's Bayes factor without needing individual-level genotypes or an explicit LD matrix.
- The two datasets' SNP sets meaningfully overlap across the tested region.

**Region definition.** For the Upload route, the region is either the full harmonised SNP set (no window applied, `use_window = FALSE`, the default) or a user-defined ±`window_kb` window centred on either the lowest-mQTL-p-value SNP ("lead SNP") or a manually specified position (`:609-615`). For the Preloaded route, the region was fixed upstream at ±1 Mb per CpG and is not adjustable (`MCOL_DEFAULT_WINDOW_KB = 1000` kb, matching the code comment "08d_mr_coloc.R's own CIS_WINDOW_BP = 1e6," `:57`).

**Variant matching.** Purely by SNP-ID string match, first for the raw pre-overlap check (`:483-484`), then implicitly inside `TwoSampleMR::harmonise_data()`'s own merge logic. No positional (chr:pos) matching or liftover fallback exists.

**Allele harmonisation.** Delegated entirely to `TwoSampleMR::harmonise_data()` (§10); the module's own code never inspects or flips alleles directly.

**Effect-size handling.** Betas and SEs are passed through unmodified from the harmonised table into `coloc.abf()`'s `beta`/`varbeta` (`= se^2`) arguments (`:632-633`, `:637-638`) — no rescaling, standardisation, or unit conversion is applied.

**P-value handling.** P-values are used only for the optional pre-`coloc.abf()` filtering step (`f_pval_meth`/`f_pval_gwas`, default = 1, i.e., no filter) and for plotting (`-log10(p)`); they are never fed into `coloc.abf()` itself, which uses effect sizes and SEs, not p-values, to build its Bayes factors.

**Sample-size handling.** `N` is taken as the **median** of the harmonised subset's `samplesize.exposure`/`.outcome` column (`stats::median(..., na.rm = TRUE)`, `:633`, `:638`) — a single scalar per dataset, as `coloc.abf()` requires, rather than a per-SNP value (`coloc.abf()`'s API does not accept per-SNP N).

**MAF requirements.** Required for the methylation (quantitative) dataset — enforced by an explicit `validate(need(!is.null(d1$MAF), ...))` (`:635-636`). For the GWAS dataset, MAF is used only if `gwas_type == "quant"` and an EAF column was mapped (`:640`); for `gwas_type == "cc"` (the default), the case fraction (`s`) substitutes for MAF, matching `coloc.abf()`'s own API for case-control traits.

**LD requirements.** None for `coloc.abf()` (by design — this is precisely what the single-causal-variant approximation avoids needing). Required, and only ever user-supplied, for `coloc.susie()`.

**Priors.** `p1 = 1e-4`, `p2 = 1e-4`, `p12 = 1e-5` are `coloc`'s own package defaults (stated explicitly in the code comment, `:49-51`), pre-filled as the Upload route's defaults and echoed read-only for the Preloaded route (since they were fixed at cache-generation time). All three are user-adjustable sliders in the Upload route.

**Posterior probabilities / decision thresholds.** PP.H0–H4 sum to 1 by construction (a property of `coloc.abf()`'s output, not separately enforced by this module). The verdict threshold is `pp_threshold`, default 0.8, matching "this project's own script08d threshold" (`:270-271`).

**Interpretation of a significant result.** PP.H4 ≥ threshold → "coloc-supported (shared causal variant)"; PP.H3 ≥ threshold → "coloc-refuted (distinct causal variants)"; neither → "inconclusive." This module never claims a p-value-style "significance" — everything is reported as a posterior probability.

---

## 12. Data Preprocessing

**Upload route** (`build_run_state_upload()`, `:604-708`), applied to the harmonised data.frame `h`, in this exact order:

1. **Genomic window** (only if `use_window`): subset to `abs(pos - center) <= window_kb * 1000` (`:609-615`).
2. **Association filter**: `pval.exposure <= f_pval_meth & pval.outcome <= f_pval_gwas` (`:616`).
3. **MAF filter**: `maf_val` computed as `pmin(eaf, 1-eaf)` (folded to a true minor-allele frequency), kept if in `[f_maf_min, f_maf_max]` or missing (`:617-619`).
4. **Deduplication** (if `f_dedup`): `!duplicated(h$SNP)` (`:620`) — keeps the **first** occurrence encountered, not, e.g., the lowest-p-value one.
5. **Minimum sample size**: rows with `samplesize < f_min_n` on either dataset are dropped, unless `NA` (`:621-623`).
6. **Hard floor check**: `nrow(h) >= f_min_shared` or the run aborts with a validation message (`:626-628`).
7. **SE sanity check**: all retained SEs strictly positive on both datasets, or the run aborts (`:629`).

**Preloaded route.** No live preprocessing — only a post-hoc `nsnps >= pre_min_nsnps` row filter on the already-finished `coloc_results.csv` table (`:596-598`); the underlying instrument selection, LD clumping, and coloc run were performed entirely upstream, outside this application.

---

## 13. Variant/CpG Harmonization

Handled exclusively by `TwoSampleMR::harmonise_data()` (see §10), invoked once per Upload-route Validate click, with `action` set from `input$remove_ambiguous`:
- `action = 2` (default, `remove_ambiguous = FALSE`): palindromic SNPs are resolved via EAF where possible; unresolvable ones are dropped.
- `action = 3` (`remove_ambiguous = TRUE`): all palindromic SNPs are dropped outright, regardless of EAF availability.

Post-harmonisation, the module reports (and the user can inspect in the Data-validation summary box) counts of `aligned` (`mr_keep`), `palindromic`, `ambiguous`, and `removed` variants (`:514-516`). CpG-level harmonisation is not a distinct step — a single target CpG's rows are simply subset from the uploaded methylation file before formatting (`:469`); there is no cross-CpG matching or merging logic in this module.

---

## 14. Colocalization Computation

**Preloaded route:** no computation — `build_run_state_preloaded()` filters the cached `coloc_results.csv` table by `pre_min_nsnps` and returns it directly (`:593-602`); `coloc.abf()` is never called on this path.

**Upload route (`coloc.abf`):** After preprocessing (§12), `d1`/`d2` lists are built and passed to `coloc::coloc.abf(dataset1 = d1, dataset2 = d2, p1, p2, p12)` inside `tryCatch()`/`suppressWarnings()` (`:642`). The per-SNP results (`abf_res$results`) are joined back onto the harmonised table by `match(h$SNP, abf_res$results$snp)` to build `snp_df` (`:645-657`), which is then sorted by descending `snp_pp_h4` and ranked. The lead variant is `which.max(snp_df$snp_pp_h4)` (`:659`).

**Upload route (`coloc.susie`, optional):** Only reachable if `use_susie` is checked, both LD files are supplied, `susieR` is installed, and `.mcol_prep_ld()` finds ≥3 common SNPs on both matrices (`:662-694`). The credible-set-pair summary (`susie_res$summary`) is stored and, if empty/all-`NA`, treated as "no credible set found" rather than an error (`:689-692`).

**Region reporting.** A `region` object (`chr`, `start`, `end`) is derived from the harmonised table's own `chr.exposure`/`pos.exposure` columns when available (`:697-699`) — purely descriptive, not fed back into the analysis.

---

## 15. Output Generation

See §16–§18 for full detail. In summary: the Results tab renders value-box summaries + a hypothesis-probability table (+ SNP-level table and optional SuSiE table for Upload); the Visualisation tab renders up to three/four `ggplot2` plots; the Sensitivity tab renders one plot + two tables (Upload only); the Export tab exposes seven to ten `downloadHandler()`s depending on route and whether plotting/sensitivity have run.

---

## 16. Tables

| Table (output ID) | Route | Source | Content |
|---|---|---|---|
| `pre_preview_table` | Preloaded | `pre_coloc_tbl()` filtered to selected CpGs | `(cpg, nsnps)` only — descriptive preview, no PPs |
| `meth_preview_table` / `gwas_preview_table` | Upload | `meth_df_r()`/`gwas_df_r()` | First 20 raw rows of each uploaded file |
| `pre_results_table` | Preloaded | `run_state()$table` | Full per-CpG PP.H0–H4 + computed `verdict` |
| `snp_table` | Upload | `run_state()$snp_df` | Per-SNP chr/pos/alleles/MAF/betas/p-values/log-ABF/`snp_pp_h4`/rank |
| `susie_table` | Upload (if run) | `run_state()$susie_res$summary` | One row per credible-set pair, its own PP.H0–H4 |
| `sens_prior_table` | Upload | `sensitivity_state()$prior` | p12 grid × PP.H3/PP.H4 |
| `sens_param_table` | Upload | `sensitivity_state()$param` | One row per perturbed filter value × nsnps/PP.H3/PP.H4 |

All tables use `DT::datatable()`; numeric formatting via `DT::formatRound()` (PPs, 4 dp) or `DT::formatSignif()` (p-values/log-ABF/PPs, 4 significant figures).

---

## 17. Figures

| Plot (output ID) | Builder | X | Y | Grouping/colour | What it shows |
|---|---|---|---|---|---|
| `pre_pp_plot` | `build_pp_bar_plot()` | H0–H4 (factor) | posterior probability | H4 highlighted red | The five posterior probabilities for one selected CpG (Preloaded) |
| `region_plot` | `build_region_plot()` | SNP position (bp) | −log10(p) | faceted by track (methylation vs. GWAS), coloured by track | Regional association "Manhattan-style" comparison; a dashed vertical line marks the lead variant |
| `comparison_plot` | `build_comparison_plot()` | mQTL −log10(p) | GWAS −log10(p) | continuous colour = `snp_pp_h4` | Whether the SNPs most significant in one trait are also significant in the other, and whether those SNPs carry high H4 posterior weight |
| `posterior_plot` | `build_posterior_plot()` | SNP position (bp), or SNP ID (top 30) if no position mapped | `snp_pp_h4` | lead variant highlighted red | Which specific SNP(s) carry the posterior weight for a shared signal |
| `sens_prior_plot` | inline in `output$sens_prior_plot` | p12 (log scale) | posterior probability | line colour = H3 vs. H4 | How stable the H3/H4 split is as the shared-association prior varies |

**Statistical annotations.** Only the region/posterior plots' dashed/highlighted lead-variant marker; no confidence bands, LD-based point shading, or gene-track annotation is drawn on any plot (these require external gene/LD-reference data this module does not load).

**Filtering/thresholds shown in plots.** None of the plots apply an additional filter beyond what already went into `run_state()` — they visualise exactly the SNP set that was analysed, not a re-filtered subset.

**Misleading/incomplete visualization choices identified from the code.** The `posterior_plot`'s fallback branch (no position mapped) silently caps display at the top 30 SNPs by `snp_pp_h4` (`:922`) with no on-plot note of how many were omitted — a user comparing this plot's apparent variant count against the SNP-level table's full count could be misled about coverage (see §26, finding C-1).

---

## 18. Downloads

| Download button | Filename | Source object | Content | Reproducible? |
|---|---|---|---|---|
| `dl_hypotheses` | `coloc_hypothesis_probabilities.csv` | `run_state()` | Per-CpG table + verdict (Preloaded) or one-row PP summary + verdict (Upload) | Yes — deterministic from `run_state()` |
| `dl_validation` | `coloc_validation_summary.csv` | `validate_state()$summary` | Variant/dedup/missingness/harmonisation counts | Yes |
| `dl_params` | `coloc_analysis_parameters.csv` | `run_state()` | Mode, priors, method, sample sizes, threshold (Upload) or fixed priors/window/method (Preloaded) | Yes |
| `dl_snp` | `coloc_snp_results_<cpg>.csv` | `run_state()$snp_df` | Full per-SNP results | Yes (Upload only) |
| `dl_shared` | `coloc_shared_variants_<cpg>.csv` | `snp_df[snp_pp_h4 >= 0.01, ]` | High-posterior-weight variant subset | Yes, but the 0.01 cutoff is **hardcoded** and not shown in the UI (see §26, finding C-4) |
| `dl_sensitivity` | `coloc_sensitivity_analysis.csv` | `sensitivity_state()$prior` | Prior-sensitivity grid (empty data.frame if that analysis was never run) | Yes |
| `dl_plot_region`/`_comparison`/`_posterior` | `coloc_plot.<fmt>` | re-invokes the corresponding `build_*_plot()` | PNG/PDF/SVG, 9×6in, 300 dpi | Yes — regenerated from `run_state()`, not a cached image |

None of these downloads embed a session ID, timestamp, or package-version stamp — reproducing a downloaded CSV's exact provenance later requires cross-referencing it against the separately downloaded `dl_params` file (see §25).

---

## 19. Reactive Programming and Tab Connectivity

**Stage-gating mechanism.** `stage_flags <- reactiveValues(validate=FALSE, run=FALSE, plot=FALSE, sensitivity=FALSE)` (`:303`) plus `invalidate_from(stage)` (`:304-307`), which walks `.mcol_stage_order <- c("validate","run","plot","sensitivity")` (`:47`) from the given stage to the end, resetting each flag to `FALSE`. Every `renderUI()` that shows stage-gated content checks its flag with `req()`, and every one of those `renderUI()`s is wrapped in `outputOptions(output, ..., suspendWhenHidden = FALSE)` so it recomputes even while its tab isn't visually active (necessary because Tab 3's content must exist before the user switches to it via `updateTabsetPanel()` on a successful run, `:730`).

**Invalidation triggers, verified from the code:**
- Tab 1 defining inputs → `invalidate_from("validate")` (`:376-378`): `data_source`, `pre_cpgs`, `meth_file`, `gwas_file`, `meth_cpg`, `meth_target_cpg`, `gwas_label`, `gwas_type`, `case_frac`.
- Tab 2 defining inputs → `invalidate_from("run")` (`:581-586`): every filter/prior/method input.
- `validate_btn` click → `stage_flags$validate <- TRUE; invalidate_from("run")` (`:539-541`).
- `run_btn` click → `stage_flags$run <- TRUE; invalidate_from("plot")` (`:727-728`), plus a write into the shared cross-module `results$coloc` reactiveValues if a `results` argument was supplied to the module (`:732-743`) — see §23, finding B-6 for whether this write is ever consumed.
- `plot_btn` click → `stage_flags$plot <- TRUE` (`:967`) — this stage flag is set but **never checked anywhere downstream** (the actual plot gating uses a separate client-side `conditionalPanel(input[...] > 0)` check, `:943-944`, `:1024-1025`, `:1036-1037`, and each `render*` output re-checks `input$plot_btn > 0`/`input$sens_*_btn > 0` directly) — a redundant, unread piece of state, not a bug (see §26, finding C-5).
- `sens_prior_btn`/`sens_param_btn` clicks → set `stage_flags$sensitivity <- TRUE` (also unread elsewhere) and populate `sensitivity_state()`.

**Dependency chain in plain English.** Data & Setup → (Validate Data) → Filters & Parameters (route-specific body only appears post-validation) → (Run Colocalisation) → Results/Visualisation/Sensitivity/Export (all four unlock together on a successful run; Visualisation and Sensitivity additionally require their own in-tab button click before rendering plots/tables). Tabs are strictly **sequential**, not independent: running Tab 4 or Tab 5 without having completed Tab 1's Validate Data and Tab 3's Run Colocalisation is impossible, because their `renderUI()`s `req(stage_flags$run)`/`req(stage_flags$validate)` first.

**Objects shared between tabs.** `validate_state()`, `run_state()`, and `sensitivity_state()` — three `reactiveVal()`s, each written by exactly one stage's observer and read by every tab downstream of that stage. No tab writes to another tab's `reactiveVal()`.

**Unnecessary recomputation / stale outputs.** Not identified — every `render*` reads from the appropriate `reactiveVal()`/`stage_flags` combination, and `invalidate_from()` consistently clears every downstream stage on an upstream change, so a stale Results/Plot/Sensitivity view is not reachable through normal interaction.

**Missing dependencies / accidental triggering.** Not identified in the observed reactive graph — see §26 for the two informational (non-bug) observations about unread flags above.

---

## 20. End-to-End Data Pipeline

```text
User Input (Tab 1: choose Preloaded/Upload; upload files + map columns, OR pick CpG(s))
   ↓
Input Validation ("Validate Data" button → build_validate_state_upload()/_preloaded())
   [Upload]: file-readable check → CpG-subset check → raw SNP-ID overlap check
   [Preloaded]: CpG-selection check → coloc_results.csv row-exists check
   ↓
Data Loading
   [Upload]: read_uploaded_table() on both files
   [Preloaded]: load_default_meth_coloc_results() (+ load_default_mr_harmonised() for context only)
   ↓
Identifier/Variant Processing (Upload only)
   TwoSampleMR::format_data() ×2 → TwoSampleMR::harmonise_data() → mr_keep==TRUE subset
   ↓
Region Definition (Tab 2, Upload only)
   optional genomic window around lead SNP or a manual position
   ↓
Dataset Harmonisation
   already complete by this point (harmonise_data() ran during Validate Data)
   ↓
Statistical Preparation ("Run Colocalisation" button → build_run_state_upload()/_preloaded())
   [Upload]: association/MAF/N/dedup filters → min-shared-SNP floor → build d1/d2 lists
   [Preloaded]: filter cached table by pre_min_nsnps only
   ↓
Colocalisation Analysis
   [Upload]: coloc::coloc.abf(); optionally coloc::coloc.susie() if LD supplied
   [Preloaded]: none — cached result is looked up, not computed
   ↓
Posterior/Statistical Results
   PP.H0–H4 (+ per-SNP SNP.PP.H4, log-ABF) [Upload]; per-CpG PP.H0–H4 [Preloaded]
   ↓
Filtering/Classification
   .mcol_verdict() → coloc-supported / coloc-refuted / inconclusive
   ↓
Tables (Tab 3: Results)
   ↓
Plots (Tab 4: Visualisation — gated behind its own "Generate Regional Plot" click)
   ↓
Sensitivity (Tab 5, Upload only — gated behind its own two buttons)
   ↓
Downloads (Tab 6: Export)
```

Every stage above exists in the reviewed code exactly as shown; no stage was added for scientific completeness. The Preloaded route visibly skips "Statistical Preparation," "Colocalisation Analysis," and "Sensitivity" as live compute steps — this is a genuine data-availability constraint stated repeatedly in the code's own comments (`:14-25`, `:158`, `:939`, `:1010-1012`), not an implementation gap.

---

## 21. Code-to-Science Mapping

| Code operation | Computational meaning | Biological meaning |
|---|---|---|
| `read_uploaded_table()` / `load_default_meth_coloc_results()` | Parse a delimited file / cached CSV into a data.frame | Ingest a methylation-QTL or GWAS association study's summary statistics, or a previously completed colocalization run |
| `format_data()` ×2 | Standardise column names/types to a common schema | Recognise which numbers are effect sizes, SEs, alleles, etc., regardless of the source study's naming convention |
| `harmonise_data()` | Merge by SNP ID, flip signs/alleles, flag palindromic rows | Ensure "the effect of the A allele" means the same physical allele in both studies before comparing them |
| genomic-window / p-value / MAF / N filters | Row-subset a data.frame | Restrict to the biologically relevant local region and to variants informative and reliable enough to model |
| `coloc.abf(dataset1, dataset2, p1, p2, p12)` | Bayesian model comparison over 5 hypotheses from summary statistics | Test whether the methylation-associated and disease-associated signals in this region are compatible with one shared causal variant |
| `coloc.susie()` | SuSiE fine-mapping + pairwise credible-set colocalization | Test the same question while allowing for more than one independent causal variant per trait in the region |
| `PP.H4` / `SNP.PP.H4` | Posterior probability of the shared-signal hypothesis / per-SNP posterior weight under H4 | Overall confidence that the two traits share a genetic cause here / which specific variant most plausibly is that cause |
| `.mcol_verdict()` | Threshold PP.H3/PP.H4 into a label | Translate a continuous probability into the same "supported/refuted/inconclusive" categories the thesis chapter itself reports |
| `coloc::sensitivity()` / `.mcol_rerun_with()` | Re-derive PP under alternative priors/filters | Check whether the biological conclusion is a robust feature of the data or an artefact of one modelling choice |

---

## 22. Package and Dependency Audit

| Package | Version (renv.lock / DESCRIPTION) | Function(s) used | Why used |
|---|---|---|---|
| `coloc` | 5.2.3 (`DESCRIPTION:25`, `renv.lock:61-65`) | `coloc.abf()`, `coloc.susie()`, `sensitivity()` | Core colocalization engine |
| `TwoSampleMR` | 0.7.8 (`DESCRIPTION:101`, `renv.lock:421-424`, `Source: Unknown`) | `format_data()`, `harmonise_data()` | Standardisation/harmonisation shared with `mod_methyl_mr.R` |
| `susieR` | **Not present in `renv.lock`** | `coloc::coloc.susie()`'s internal fine-mapping dependency | Optional — only exercised if installed; module checks with `requireNamespace()` and degrades gracefully if absent |
| `data.table` | 1.18.2.1 | `fread()` (via `read_uploaded_table()`) | Fast, delimiter-flexible file reading |
| `DT` | (project-wide dependency, not individually re-checked here) | `datatable()`, `renderDataTable()`, `formatRound()`, `formatSignif()` | Interactive result tables |
| `ggplot2` | (project-wide dependency) | plotting + `ggsave()` | All four plots and their file export |
| `shiny`, `shinycssloaders` | (project-wide) | reactive framework, `withSpinner()` | UI/reactivity/loading indicator |

**Indirect dependencies.** `susieR` (via `coloc.susie()`) is the only indirect dependency this module can invoke; it is correctly treated as optional.

**Unused imports.** None identified — every package function listed in §9.2 has at least one live call site in this file.

**Potentially missing dependency (reproducibility concern).** Because `susieR` is not pinned in `renv.lock`, whichever `susieR` version happens to be installed on a given deployment (or none at all) determines whether, and with what numerical behaviour, `coloc.susie()` results can be reproduced across deployments — see §25.

---

## 23. Scientific Correctness Audit

Classified per the audit brief's scale (Correct / Minor / Moderate / Major / Critical); every item is backed by the code cited.

**A. Correct.**
- Allele harmonisation is fully delegated to `TwoSampleMR::harmonise_data()` and restricted to `mr_keep==TRUE` rows before analysis — the standard, published approach (`:505-518`).
- No LD is fabricated anywhere; `coloc.susie()` is unreachable without genuine user-supplied LD matrices, and `.mcol_prep_ld()` returns `NULL` rather than guessing (`:65-76`).
- The quantitative-trait MAF requirement for `coloc.abf()` is explicitly checked and blocks the run with an actionable message if missing (`:635-636`) — a real, non-hypothetical failure mode this code correctly anticipates.
- The verdict thresholds and priors mirror the upstream thesis pipeline's own published methodology exactly (cross-checked against `METHODS_mendelian_randomization.md` §2.FF.2 and `MCOL_DEFAULT_*` constants, `:53-59`).
- Posterior probabilities are always reported as probabilities, never re-labelled as p-values or "significance" — no multiple-testing correction is applied because none is claimed to be needed (each CpG's colocalization test is presented and interpreted independently, matching how the source thesis chapter itself reports these results, per `METHODS_mendelian_randomization.md` §2.FF.3).

**B. Minor/Moderate concerns.**
- **B-1 (Moderate).** In `build_validate_state_preloaded()` (`:436-456`), the "Dataset 1 variants" and "Dataset 2 variants" value boxes shown in the Data-validation summary are both computed as `length(unique(ctx$SNP))`, where `ctx` is a subset of `mr_harmonised_all_cpgs.csv` — the **already-LD-clumped MR instrument table** (`load_default_mr_harmonised()`, `global.R:562-567`: "already clumped and already `harmonise_data(action=2)`'d"). This is a **different, much smaller** SNP set than the one `08d_mr_coloc.R` actually fed into `coloc.abf()` for the cached result, which per `METHODS_mendelian_randomization.md` §2.FF.2 used "the same candidate-list rows" (the full, unclumped GoDMC candidate list) — the value the `nsnps` column in `coloc_results.csv` itself correctly reflects (used separately for "Shared variants"/"Final variants available" in the same summary row, `:452-453`). A user reading all four value boxes side by side could reasonably assume "Dataset 1 variants" describes the same analysis population as "Shared variants," when it does not.
- **B-2 (Minor).** Invalid p-values (outside [0,1]) are counted but not excluded before `format_data()`/`harmonise_data()`/`coloc.abf()` (`:481-482`, no corresponding filter step in `build_run_state_upload()`). `coloc.abf()` does not itself validate p-value range (it does not consume p-values at all, only beta/varbeta), so this has no direct numerical effect on the colocalization result, but a downstream **plot** (`build_region_plot()`, `build_comparison_plot()`) computing `-log10(p)` from an out-of-range p-value (e.g., negative) would produce a `NaN` point silently dropped by ggplot2, without any on-screen warning.
- **B-3 (Minor).** Duplicate SNPs are removed with `!duplicated(h$SNP)` (`:620`), which keeps the *first row encountered* rather than, e.g., the one with the smallest p-value or largest N — for a file with legitimate multi-row duplicates (different imputation batches, for example), which row survives is an artefact of input row order, not a principled choice. This is disabled entirely if the user unchecks `f_dedup`.
- **B-4 (Minor).** No genome-build check exists between the uploaded mQTL file, the uploaded GWAS file, and (if used) the CpG-position annotation implied by the "SNP position" column — a build mismatch (e.g., one file in GRCh37, the other GRCh38) would not be detected and would silently produce a scientifically meaningless "region" and window filter.
- **B-5 (Minor).** `coloc.susie()`'s `p12` argument is hardcoded to `MCOL_DEFAULT_P12_SUSIE = 5e-6` (`:685`) rather than following the user's `p12` slider (which does apply to `coloc.abf()`, via the shared `p1`/`p2`/`p12` variables computed at `:631`). The Priors panel's own tooltip for `p12` ("Prior probability a SNP is associated with BOTH traits. coloc default: 1e-5," `:248-249`) does not mention this exception, so a user who deliberately raises or lowers `p12` and then enables `coloc.susie()` would see the SuSiE-based result computed at a different, silently substituted prior.
- **B-6 (Minor/informational).** On a successful Upload-route run, the module writes a summary entry into the shared cross-module `results$coloc` reactiveValues list (`:732-743`), following the same pattern `mod_methyl_mr.R` and other Genetics-group modules use to publish results for potential downstream aggregation (e.g., a biomarker-card view). A repo-wide search (`grep -rn "results\$coloc" R/`) found **no other module currently reads this key** — the write is not incorrect, but is presently inert integration scaffolding rather than an active cross-module data flow.

No **Major** or **Critical** correctness issues were identified in the reviewed code.

---

## 24. Software Engineering Audit

- **Error handling.** Every fallible external call (`format_data()`, `harmonise_data()`, `coloc.abf()`, `coloc.susie()`, `coloc::sensitivity()`) is wrapped in `tryCatch()`, and user-facing `validate(need(...))`/`showNotification()` calls surface actionable messages rather than raw R errors, except `coloc::sensitivity()`'s failure path (`:1053`), which degrades silently to an empty panel with no notification — a minor, non-scientific UX gap (§19).
- **Hard-coded values.** `MCOL_DEFAULT_*` constants (`:53-59`) are the module's package/pipeline defaults, appropriately centralised rather than scattered as magic numbers; the one exception worth flagging is the `snp_pp_h4 >= 0.01` cutoff inside `dl_shared`'s `downloadHandler` (`:1176`), which is not a named constant, not shown anywhere in the UI, and not documented in the download's button label ("Shared/high-posterior variants (CSV)").
- **Hard-coded file paths.** None inside this file — all Preloaded-route paths are resolved via `METH_MR_DIR` (`data_paths.R:94`) inside the shared `global.R` loader functions, not duplicated here.
- **Unused variables/functions.** None identified — every helper and constant defined in this file has at least one call site.
- **Duplicate code / redundant computation.** The window-centring and filter-application logic in `build_run_state_upload()` (`:609-623`) and `.mcol_rerun_with()` (`:987-995`) is structurally similar but intentionally separate (the latter is a standalone, side-effect-free re-run helper reused across many sensitivity iterations) — not flagged as problematic duplication, since unifying them would require passing `input` reactives into a pure function, complicating the sensitivity loop.
- **Namespace conflicts.** None identified; every package call in this file is namespace-qualified (`coloc::`, `TwoSampleMR::`, `DT::`, `stats::`, `utils::`) or uses a base function with no ambiguity.
- **Deprecated functions.** None identified against the pinned package versions (`coloc` 5.2.3, `TwoSampleMR` 0.7.8).
- **Empty-result handling.** Explicitly handled at every stage: empty post-filter tables (`validate(need(nrow(...) > 0))`), an empty/absent `coloc.susie()` credible set (`:689-692`), and an empty sensitivity grid (`dl_sensitivity` falls back to `data.frame()`, `:1184`) are all anticipated rather than left to crash.
- **Performance.** `.mcol_rerun_with()` re-runs `coloc.abf()` up to roughly 8–10 times per "Run Sensitivity Analysis (Parameter)" click (baseline + up to 3 window multipliers + up to 3 MAF values + up to 3 p-value thresholds, each conditionally skipped if it equals the baseline, `:1076-1081`); for typical cis-window SNP counts (tens to low thousands, per `coloc_results.csv`'s own `nsnps` column, max 1,005) this is a cheap, sub-second-per-call operation, not a performance risk.
- **Reproducibility / randomness.** No random-number generation occurs anywhere in this module (`coloc.abf()`, `coloc.susie()`'s SuSiE step, and `sensitivity()` are all deterministic given their inputs) — no seed is needed and none is set.
- **Download reproducibility.** Every download regenerates its content from `run_state()`/`validate_state()`/`sensitivity_state()` at click time rather than serving a stale cached file; see §25 for the separate question of reproducing the *analysis* later.

---

## 25. Reproducibility Audit

| Criterion | Status | Evidence |
|---|---|---|
| Required inputs documented? | Yes | UI text explicitly states required columns/format per file (`:167`, `gwas_col_map_ui` labels) |
| Parameters visible? | Yes | All priors/filters are explicit numeric/checkbox inputs, not buried defaults |
| Defaults documented? | Yes | `MCOL_DEFAULT_*` constants + inline tooltips (`.mcol_tip()`) state each default and its rationale |
| Statistical method identifiable? | Yes | `coloc.abf`/`coloc.susie`/`coloc::sensitivity` named directly in UI text and outputs |
| Package versions controlled? | Partially | `coloc` 5.2.3 and `TwoSampleMR` 0.7.8 pinned via `renv.lock`; `susieR` is **not** pinned (§22) |
| Genome build specified? | Preloaded: yes (GRCh37, stated) — Upload: no, and not checked (§23, B-4) |
| Reference datasets documented? | Yes, for Preloaded (GoDMC/Ishigaki 2022, cited in UI and `METHODS_mendelian_randomization.md`) — N/A for Upload (user-supplied) |
| Randomness controlled? | N/A — no stochastic step exists |
| Results downloadable? | Yes — 7–10 CSV/plot downloads depending on route/stage (§18) |
| Intermediate data recoverable? | Partial — the harmonised, pre-`coloc.abf()` SNP set (`h`/`h_full`) is **not** independently downloadable; only the post-`coloc.abf()` `snp_df` is |
| Analysis assumptions documented? | Yes — single-causal-variant assumption, priors, and the causality caveat are all stated in-app (`.mcol_interpret()`, `:99-100`; Priors tooltips) |

**Reproducibility score: Moderate.**

**Justification.** The Upload route is highly transparent about its own parameters, defaults, and method, and every output is regenerable from `run_state()` at any time within a session — but a researcher wanting to reproduce a *specific past result outside the app* would need to separately retain the `dl_params` CSV (parameters), the `dl_snp`/`dl_hypotheses` CSV (results), and the original input files, since no single export bundles genome build, package versions, or a run timestamp together with the analysis parameters. The Preloaded route is fully documented at the pipeline-methods level (`METHODS_mendelian_randomization.md`) but is, by the module's own design, a **look-up of a fixed historical result**, not a live-reproducible computation within this application — `coloc.abf()` cannot be re-run for it here regardless of what the user does, because the per-SNP GoDMC/RA-GWAS rows are not bundled with the deployment.

---

## 26. Audit Findings

| ID | Finding | Severity | Code Location | Why It Matters | Recommended Action |
|---|---|---|---|---|---|
| A-1 | `coloc.susie()`'s `p12` is hardcoded to `5e-6`, ignoring the user's `p12` slider (which does apply to `coloc.abf()`) | Minor | `mod_methyl_coloc.R:685`, cf. `:631` | A user adjusting `p12` and enabling SuSiE would see an unexplained prior substitution | *Audit recommendation (not implemented): document this exception directly in the SuSiE `conditionalPanel` note (`:255-256`), or add a distinct, clearly-labelled `p12_susie` input.* |
| B-1 | Preloaded-route "Dataset 1/2 variants" value boxes are computed from the LD-clumped MR instrument table, not the unclumped candidate set actually used by the cached coloc run | Moderate | `mod_methyl_coloc.R:443-454` | Could mislead a reader comparing this figure against "Shared variants"/`nsnps` in the same summary row | *Audit recommendation (not implemented): relabel these two value boxes to make clear they describe the MR instrument context, not the coloc input set, or drop them from this summary.* |
| B-2 | Out-of-range p-values are counted but not filtered before `coloc.abf()`/plotting | Minor | `mod_methyl_coloc.R:481-482` | No numerical effect on `coloc.abf()` itself (which never reads p-values); could silently blank a plotted point | *Audit recommendation (not implemented): filter or flag invalid p-value rows before the region/comparison plots.* |
| B-3 | Deduplication keeps the first-encountered duplicate SNP row, an artefact of input order | Minor | `mod_methyl_coloc.R:620` | Result may differ across otherwise-identical files with reordered rows | *Audit recommendation (not implemented): document the tie-break rule in the UI tooltip.* |
| B-4 | No genome-build consistency check across the two uploaded files | Minor | `mod_methyl_coloc.R:604-708` (absence) | A build mismatch would silently corrupt the region/window filter | *Audit recommendation (not implemented): add an optional build-declaration input per file.* |
| C-1 | Fallback posterior plot silently caps display at the top 30 SNPs with no on-plot note | Minor | `mod_methyl_coloc.R:922-928` | Understates apparent variant coverage relative to the full SNP table | *Audit recommendation (not implemented): add a plot subtitle stating "(top 30 of N)".* |
| C-2 | `coloc::sensitivity()` failures degrade silently (no `showNotification`) | Minor | `mod_methyl_coloc.R:1044-1058` | User sees an empty panel with no explanation of why | *Audit recommendation (not implemented): add a `showNotification` on the `inherits(sens, "error")` branch.* |
| C-3 | `.mcol_rerun_with()`'s internal `nrow(h) < 6` floor is stricter than, and independent of, the user-facing `f_min_shared` default (10), and drops silently | Minor | `mod_methyl_coloc.R:996` | A perturbed sensitivity row can vanish from the table with no explanation | *Audit recommendation (not implemented): surface a "insufficient SNPs at this setting" row instead of omitting it.* |
| C-4 | The "Shared/high-posterior variants" download's `snp_pp_h4 >= 0.01` cutoff is hardcoded and not shown in the UI | Minor | `mod_methyl_coloc.R:1176` | A user cannot tell what threshold defines "shared" in this specific export without reading the source | *Audit recommendation (not implemented): state the threshold in the download button's own label.* |
| C-5 | `stage_flags$plot`/`stage_flags$sensitivity` are set but never read anywhere (actual gating uses `input$plot_btn > 0` etc. directly) | Informational | `mod_methyl_coloc.R:967`, `:1057`, `:1087` | Dead reactive state, no functional effect | *No action needed — informational only.* |
| B-6 | `results$coloc` is written on every successful Upload run but read by no other module in this codebase | Informational | `mod_methyl_coloc.R:732-743` | Currently inert integration scaffolding | *No action needed unless a downstream consumer is planned.* |

### Correctly Implemented
Allele harmonisation delegation; MAF pre-flight check for quantitative traits; no-fabricated-LD guarantee for `coloc.susie()`; verdict/threshold logic matching the upstream published methodology; deterministic, seed-free computation throughout; comprehensive `tryCatch()`/`validate()` coverage around every external statistical call except one (`coloc::sensitivity()`, see C-2).

### Minor Issues
A-1, B-2, B-3, B-4, C-1, C-2, C-3, C-4.

### Moderate Issues
B-1.

### Major Issues
No implementation issue of Major severity was identified for this component based on the reviewed code.

### Critical Issues
No implementation issue of Critical severity was identified for this component based on the reviewed code.

---

## 27. Educational Interpretation

1. **What biological question is being asked?** Does the same genetic variant plausibly explain both a methylation difference at a CpG and a disease-association signal nearby, or are these two separate genetic effects that merely sit close together on the chromosome?
2. **What data are needed?** Regional summary statistics (SNP ID, effect size, standard error, and either allele frequency or a case fraction) for both the methylation/mQTL signal and the GWAS trait, covering the same genomic window.
3. **Why are two association signals compared, not one?** Because overlap alone (both studies have "a hit near this CpG") cannot distinguish a shared cause from two neighbouring, LD-correlated causes — colocalization is the statistical machinery built specifically to make that distinction using each study's *full* regional pattern.
4. **What does LD have to do with the analysis?** LD is exactly what makes H3 (two distinct causal variants) look superficially like H4 (one shared variant) in the first place; `coloc.abf()` sidesteps needing an LD matrix by working entirely from marginal per-SNP effect estimates, while `coloc.susie()` uses LD directly to separate multiple independent signals.
5. **What does a shared signal (high PP.H4) mean?** The pattern of effect sizes across the region is best explained by one variant influencing both traits.
6. **What does it not mean?** It does not identify that variant's biological mechanism, does not prove a causal direction, and — per the module's own text — is "not, by itself, proof of biological causality" (`:100`).
7. **What does a high posterior probability mean?** Strong statistical support, under the stated priors and single-causal-variant (or SuSiE multi-signal) assumption, for the corresponding hypothesis — not certainty.
8. **What does a low posterior probability mean?** The data at this locus does not distinguish this hypothesis from the alternatives well, or actively favours a different one; it is not itself evidence of "no relationship," especially when instrument/variant coverage is sparse (as the Preloaded route's own GoDMC-candidate-list limitation illustrates, `:844`).
9. **How should the result be reported in a thesis?** As a posterior probability with its threshold and priors stated explicitly (e.g., "PP.H4 = 0.87 under coloc's default priors (p1=p2=1e-4, p12=1e-5), classified coloc-supported at the pre-specified 0.8 threshold"), always alongside the caveat that colocalization establishes statistical compatibility with a shared signal, not a proven causal mechanism.

---

## 28. Thesis Implementation Description

The Colocalization sub-module of the Methylomics section of XomicShiny implements Bayesian genetic colocalization analysis between a methylation-associated quantitative trait locus (mQTL) signal and a genome-wide association study (GWAS) trait signal at a single genomic locus, organised across six sequential tabs (Data & Setup, Filters & Parameters, Results, Visualisation, Sensitivity Analysis, and Export). Two data-acquisition routes are supported: a Preloaded route that reproduces, by direct lookup, a previously completed `coloc.abf()` analysis comparing the GoDMC cis-mQTL resource against the Ishigaki et al. (2022) rheumatoid-arthritis GWAS for each candidate CpG identified by the module's upstream Mendelian-randomization panel; and a fully live Upload route in which user-supplied methylation/mQTL and GWAS summary-statistic files are column-mapped, standardised and allele-harmonised via `TwoSampleMR::format_data()` and `harmonise_data()`, filtered by genomic window, association strength, minor allele frequency, and sample size, and analysed with the `coloc` package's `coloc.abf()` function under the Approximate Bayes Factor framework of Giambartolomei et al. (2014), with an optional SuSiE-based multi-signal extension (`coloc::coloc.susie()`) available when the user additionally supplies linkage-disequilibrium reference matrices for both datasets. The implementation reports the full posterior probability distribution across the five canonical colocalization hypotheses (H0–H4), applies a pre-specified 0.8 posterior threshold — matching the upstream thesis pipeline's own decision rule — to classify each locus as coloc-supported, coloc-refuted, or inconclusive, and provides per-SNP posterior weights, regional and comparison visualisations, and a two-part sensitivity analysis (prior-grid re-weighting via `coloc::sensitivity()`, and filter-parameter perturbation via a custom re-run helper) to assess the robustness of each conclusion. Throughout, the module explicitly and consistently distinguishes statistical colocalization from proof of biological causality, a scientifically important safeguard given colocalization's frequent misinterpretation in the applied genomics literature; the audit conducted alongside this documentation identified no major or critical implementation defects, one moderate concern regarding a validation-summary display sourced from a different SNP set than the one the cached analysis actually used, and several minor, well-contained gaps in threshold visibility and edge-case notification.

---

## 29. Short Thesis Version

This sub-module implements Bayesian genetic colocalization (`coloc::coloc.abf()`, with an optional SuSiE-based multi-signal extension) between methylation-QTL and GWAS trait summary statistics at a single genomic locus, reproducing a previously completed GoDMC-vs-Ishigaki-et-al.-2022 rheumatoid-arthritis analysis by lookup, or running the same method live on user-uploaded, harmonised (`TwoSampleMR`) summary statistics. It reports posterior probabilities for five colocalization hypotheses (shared vs. distinct causal variant, or no signal), classifies each result against a pre-specified 0.8 posterior threshold matching the source thesis pipeline's own criterion, and provides prior- and filter-sensitivity analyses to assess robustness. This addresses linkage-disequilibrium-driven pleiotropy, a confound that Mendelian-randomization sensitivity tests cannot detect when instrument counts are low — precisely the scenario affecting most of the candidate CpG panel this analysis was designed to complement.

---

## 30. Overall Assessment

The Methylomics Colocalization sub-module implements a scientifically coherent, code-consistent, single-purpose analysis: Bayesian colocalization (`coloc.abf`, optionally `coloc.susie`) between one methylation-associated genetic signal and one GWAS trait signal at one genomic region, across six sequential, gate-locked tabs. Two data routes are cleanly separated — a Preloaded lookup of a previously completed, fully documented pipeline result (7 CpGs, GoDMC vs. Ishigaki et al. 2022 RA GWAS), and a fully live Upload pipeline that performs its own column mapping, allele harmonisation, filtering, `coloc.abf()`/`coloc.susie()` computation, visualisation, and sensitivity analysis. Tabs are strictly sequential (Data & Setup → Filters & Parameters → Results/Visualisation/Sensitivity → Export), enforced by a consistent stage-flag/invalidation mechanism that leaves no stale-output pathway in the reviewed reactive graph. The dominant scientific strength is the module's explicit, repeated, and technically accurate refusal to overstate colocalization as proof of causality, paired with genuinely live sensitivity analysis rather than a cosmetic robustness claim. The main implementation weaknesses are narrow and well-contained: one moderate finding (a validation-summary display sourced from a different, smaller SNP set than the one the cached Preloaded result was actually computed from, §26 B-1) and a handful of minor UI-transparency gaps (a hardcoded, unlabelled SuSiE prior override; a hardcoded, unlabelled "high-posterior" export threshold; silent truncation in one plot's fallback branch; one silently-degrading sensitivity call). No major or critical scientific or software defect was found. Reproducibility is Moderate: parameters, defaults, and method are fully transparent and every output is regenerable within a session, but no single export bundles genome build, package versions, and run parameters together, and the optional `susieR` dependency is unpinned. Based on the code actually reviewed, the implementation is scientifically appropriate for its stated, narrow purpose — a haplotype-pleiotropy check for CpGs with sparse MR instruments — and should not, and per its own in-app language does not, be read as establishing causal claims on its own.
