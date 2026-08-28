# Differential Methylation Region (DMR) Module

**Scope of this document.** This document describes and audits exactly one submodule of the Methylomics section of the ArthOMix Shiny application: the "Differentially Methylated Regions (DMRs)" module. All statements are derived from reading the actual source code, the actual bundled data files, and the actual methods write-up shipped with the app. Nothing here describes the Dataset, Quality Control, Normalisation, Cell-Type Deconvolution, or Differential Methylation Position (DMP) submodules except where the DMR module directly calls into their code — and in those cases only the specific called function is described, not the rest of those modules.

**Primary source file:** [`ArthOMix/R/methylomics/mod_methyl_dmr.R`](../../ArthOMix/R/methylomics/mod_methyl_dmr.R) (1,119 lines) — this is the *entire* DMR module; there is no second file that defines DMR-specific UI or server logic.

**Directly reused helper files** (functions called by name, not sourced separately — this app loads every `R/**/*.R` file into one global environment):
- [`ArthOMix/R/methylomics/mod_methyl_dmp.R`](../../ArthOMix/R/methylomics/mod_methyl_dmp.R) — CpG-level plotting/filtering/model-fitting helpers reused directly by the DMR engine.
- [`ArthOMix/R/methylomics/qc.R`](../../ArthOMix/R/methylomics/qc.R) — probe-level QC filters (missingness, variance, SNP overlap) and sample-ID resolution.
- [`ArthOMix/R/methylomics/annotation.R`](../../ArthOMix/R/methylomics/annotation.R) — Illumina manifest annotation (chromosome/position/SNP).
- [`ArthOMix/R/methylomics/normalization.R`](../../ArthOMix/R/methylomics/normalization.R) — normalisation-status advisory message.
- [`ArthOMix/global.R`](../../ArthOMix/global.R) and [`ArthOMix/data_paths.R`](../../ArthOMix/data_paths.R) — bundled precomputed-table loaders and path constants.
- [`ArthOMix/R/submodules_registry.R`](../../ArthOMix/R/submodules_registry.R) — registers the module into the generic Methylomics tab grid.

**Bundled ground-truth data referenced throughout this document:**
- `data/preloaded/methylomics/tables/script04_dmr_sexstratified/tables/dmr_female_full.csv` (4,657 rows) and `dmr_male_full.csv` (3,414 rows)
- `data/preloaded/methylomics/tables/script04_dmr_sexstratified/METHODS_dmr_sexstratified.md` — the thesis chapter this module's "SVA" tab reproduces.
- `data/preloaded/methylomics/tables/script03_dmp_sva_sexstratified/tables/dmp_female_full.csv` — the per-CpG input the precomputed DMR pipeline was built from.

---

## 1. Module Overview

Differential Methylation Region (DMR) analysis is a genomics method that asks whether a contiguous stretch of genomic DNA — typically containing several neighbouring CpG sites measured on the same array — shows a coordinated, spatially consistent difference in methylation between two groups of samples. Rather than testing each CpG in isolation, a DMR method aggregates the individual per-CpG evidence within a genomic window and produces one test statistic per candidate region.

DMR analysis is useful because true methylation differences linked to a biological process (disease state, cell-type shift, genetic influence) are frequently distributed across several physically adjacent, correlated CpGs rather than concentrated in exactly one probe. Testing 400,000+ individual CpGs genome-wide requires a very strict per-test significance threshold (Bonferroni/BH correction over ~4×10^5 tests); a true but modest per-CpG effect, repeated consistently across several neighbouring CpGs, can fail every single-CpG test yet still represent a real, biologically coordinated signal. Aggregating first into regions reduces the effective number of independent tests and increases power to detect exactly this pattern.

**Biological question this module is designed to answer:** are there genomic regions where the methylation level differs, in a coordinated multi-CpG fashion, between two sample groups (in the bundled dataset: rheumatoid-arthritis cases vs. controls; in the live engine: whatever two levels of a user-chosen column the user selects)?

**What this specific implementation actually does (Code fact):** The module never derives DMRs by scanning for physical clusters of independently significant DMPs. Instead, in both of its two tabs, it calls Bioconductor's `DMRcate::dmrcate()` — a Gaussian-kernel-smoothing region-calling algorithm — either on precomputed per-CpG statistics (the "SVA" tab) or on statistics computed live from a limma fit run inside this app (the "DMR" tab). `DMRcate::extractRanges()` then converts DMRcate's internal candidate-region object into a genomic-coordinate table, and this module filters, visualises, annotates, and exports that table. This is the same algorithm, not merely a scientifically similar substitute, that the reference thesis chapter used to generate the bundled tables.

**Input data the module expects:**
- **SVA tab:** no upload is possible here at all; the tab reads two precomputed CSV tables (one per sex) already bundled with the deployment.
- **DMR tab:** a beta-value or M-value methylation matrix (probes × samples) plus a sample sheet/phenotype table with a group column with ≥2 levels, loaded via the Methylomics **Dataset** tab (out of scope for this document — this document treats `methyl_dataset` only as an input the DMR module reads).

**What the user provides:** in the SVA tab, only filter/threshold choices (sex, FDR, Δβ, direction, CpG-count/width display filters) and a "Run" click. In the DMR tab, additionally: which sample-sheet column is the sex column (auto-detected, not user-entered), which column defines the two groups being compared, which two levels of that column to compare, optional covariates, QC-filter thresholds, and the DMRcate tuning parameters (seeding p-value, lambda, C, min CpGs).

**What the module produces:** for each tab, a table of candidate genomic regions with region-level effect size (mean Δβ), a region-level FDR, an overlapping-gene annotation, several plots (volcano, Manhattan, and — DMR tab only — a top-N effect-size bar chart and a region×sample heatmap), and CSV downloads of the full/significant/filtered tables plus a machine-readable record of the analysis configuration.

**What statistical/biological interpretation the outputs support:** a region with region-level FDR below the user's threshold and a mean Δβ above their chosen minimum is evidence that several nearby CpGs, considered together, show a directionally consistent methylation difference between the two compared groups, at whatever quality-control and modelling assumptions were in force for that run (see §21 and §23 for exactly which assumptions the code does, and does not, satisfy).

---

## 2. Purpose of DMR Analysis

DMR analysis exists because single-CpG (DMP) testing and region-level testing answer related but distinct questions, and because region-level testing has specific power advantages that this project's own bundled methods document explicitly invokes as the reason this module exists at all.

**Code fact (from `METHODS_dmr_sexstratified.md`, §2.BB.1, bundled at `data/preloaded/methylomics/tables/script04_dmr_sexstratified/METHODS_dmr_sexstratified.md`):** "Single-CpG testing … established that the female stratum carries genuine, calibrated genome-wide-significant signal (eighteen CpGs at bacon-corrected FDR < 0.05), while the male stratum, at n=196, retained none at the same threshold despite the underlying model being properly calibrated. … Differentially methylated region (DMR) analysis addresses this directly: it aggregates spatially adjacent, correlated CpGs into a single region-level statistic … This both reduces the effective number of independent tests … and increases power specifically when a true effect is distributed across several neighbouring CpGs."

**Scientific interpretation:** this is the textbook rationale for region-level methylation analysis (Peters et al. 2015; Peters et al. 2021, both cited in the bundled methods document) — a diffuse, low-magnitude, multi-CpG signal can be statistically invisible at the single-CpG level yet detectable once neighbouring evidence is pooled.

DMR analysis in this module is also, per the same document, motivated by downstream compatibility: region-level results are described as "the more natural unit for the Mendelian randomisation and colocalisation stage this chapter works towards, since published methylation quantitative trait locus (mQTL) resources characterise genetic influence on methylation across local, multi-CpG regions." This module itself does not perform MR or colocalisation — that is out of scope here — but the design rationale documented in the bundled methods file explains why the app maintains a dedicated DMR submodule rather than treating DMP results as sufficient on their own.

---

## 3. DMP vs DMR

| | DMP (Differential Methylation Position) | DMR (this module) |
|---|---|---|
| Unit of analysis | One CpG probe | One genomic region spanning ≥2 CpGs |
| Statistic | Per-CpG t-statistic / p-value (limma) | Per-region Stouffer-combined statistic (DMRcate), then a **separate** region-level Benjamini–Hochberg FDR computed by this app on top of it |
| What "significant" means | One probe's methylation differs between groups | Several physically adjacent probes, in aggregate, show a coordinated methylation difference |
| Correction burden | One correction across ~400k+ tests | One correction across the (much smaller) number of *candidate regions* DMRcate proposes |

**Code fact — this module does not build DMRs by post-hoc clustering of significant DMPs.** Both DMR tabs feed a per-CpG statistics table into `DMRcate::dmrcate()`, which performs its own internal Gaussian-kernel smoothing and region-calling; the module then takes DMRcate's output as the DMR table. The relationship to DMP-level information is therefore **input-level reuse of statistics, not post-hoc grouping of already-declared-significant CpGs**:

- **SVA tab:** the per-CpG input is the *same, already-computed* SVA-adjusted, bacon-corrected DMP statistics from the Differential Methylation (DMPs) tab's own default analysis (`load_default_dmp("sva", sex)`), loaded a second time here (`load_default_dmp("sva", r$sex)`, [mod_methyl_dmr.R:416](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L416)) purely to let the user inspect, per selected region, the constituent CpGs' own per-CpG Δβ/FDR — the DMR table itself is a separately bundled, precomputed CSV (`dmr_{sex}_full.csv`), not recomputed from the DMP CSV at runtime.
- **DMR tab:** the per-CpG input is a **freshly computed** limma fit ([mod_methyl_dmr.R:743-748](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L743-L748)) on whatever matrix/groups/covariates the user configured — this fit is functionally the same kind of step the DMP tab's own live engine performs, but it is re-run independently inside the DMR module's own `eventReactive`, not read from the DMP tab's reactive state (Shiny modules do not share reactive state across module instances unless explicitly passed as an argument, and `mod_methyl_dmr_server()`'s only arguments are `methyl_dataset` and `methyl_results`, neither of which carries the DMP tab's live fit).

A CpG's inclusion in a called region requires that its raw (uncorrected) per-CpG p-value fall below a seeding threshold (`is.sig`, default nominal p<0.05, not the genome-wide FDR) — this is DMRcate's own internal mechanism for proposing candidate regions, and it is a documented DMRcate tuning parameter (see §11 of the bundled methods document, quoted in §12 below), not something this app added. The **final decision about which regions are "significant"** is a separate, explicit Benjamini–Hochberg correction this app/pipeline applies to DMRcate's own Stouffer combined-probability statistic across all candidate regions in a stratum (`dmr_fdr`) — i.e., the module does define a genuine region-level multiple-testing correction, distinct from and applied on top of the CpG-level correction already embedded in the input statistics.

---

## 4. DMR Module Structure

The DMR module is a single Shiny module (`mod_methyl_dmr_ui()` / `mod_methyl_dmr_server()`) registered into the Methylomics tab grid via `mod_methyl_dmr_config` ([`submodules_registry.R:44`](../../ArthOMix/R/submodules_registry.R#L44)):

```r
mod_methyl_dmr_config <- list(
  id = "dmr", title = "Differentially Methylated Regions (DMRs)",
  icon = "map-location-dot", group = "Data",
  description = "Finds region-level methylation differences using DMRcate. Uses the bundled
                  whole-blood analysis by default, or a live, configurable run on your own dataset."
)
```
([mod_methyl_dmr.R:40-43](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L40-L43))

Internally the module is organized into two independent code sections inside one `moduleServer()` call:
1. **"1. Default analysis (GSE42861)"** — reads and filters two precomputed CSVs; performs no new statistical computation.
2. **"2. DMR Analysis (configurable live engine)"** — runs a full limma + DMRcate pipeline on a live matrix.

These two sections are exposed to the user as two UI tabs, described next.

---

## 5. Number and Overview of Tabs

**Exact count, from the UI code ([mod_methyl_dmr.R:56-66](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L56-L66)): the DMR module has exactly two tabs**, both inside one `tabsetPanel(id = ns("dmr_subtabs"), type = "tabs")`:

```r
tabPanel("SVA", br(), withSpinner(uiOutput(ns("default_ui")), color = "#2563EB", type = 6)),
tabPanel("DMR", br(), withSpinner(uiOutput(ns("live_ui")), color = "#2563EB", type = 6))
```

| Exact tab label | Wraps output ID | Internal name in code comments |
|---|---|---|
| **"SVA"** | `default_ui` | "Default analysis (GSE42861)" |
| **"DMR"** | `live_ui` | "DMR Analysis" (configurable live engine) |

No other `tabPanel`/`tabsetPanel` call exists anywhere in the file — every other UI element in the module is `fluidRow`/`div`/`uiOutput` content nested inside one of these two tabs. A code comment at the top of the file ([mod_methyl_dmr.R:45-55](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L45-L55)) states this two-tab split is "purely a UI reorganization" of what used to be two stacked sections separated by an `<hr>`, with no output ID, reactive, or server logic changed — i.e. the tab split is presentation-only and does not itself gate any computation differently than before.

---

## 6. End-to-End DMR Workflow

The generic 14-stage template in the task instructions does not match this implementation one-to-one; several stages the template lists (e.g. a distinct "preprocessing/normalisation" stage, in the sense of background/dye-bias correction) are **not implemented inside this module** — normalisation is a separate Methylomics submodule and the DMR tab only ever reads whatever matrix is already in `methyl_dataset$beta`. The actual, implemented pipeline differs between the two tabs and is given separately below, with the object and function performing each transition.

### SVA tab (precomputed) — actual pipeline

| Stage | Object | Function / mechanism | What changes |
|---|---|---|---|
| Load precomputed tables | `default_data()` reactive | `load_default_dmr("female")`, `load_default_dmr("male")`, `load_default_meth_pheno()` | Reads two CSVs (one per sex) already containing DMRcate output with `dmr_fdr` pre-computed |
| User selects stratum + thresholds, clicks Run | `d_run()` eventReactive | Reads `input$d_sex/d_fdr/d_dbeta/d_direction/d_mincpgs/d_minwidth/d_maxwidth`, adds `direction` column (`ifelse(meandiff>0,"hyper","hypo")`) | Freezes the chosen table + thresholds at click time |
| Filter for display | `d_filtered()` reactive | `mod_methyl_dmr_filter()` → `mod_methyl_dmp_filter()` | Row-subsets by FDR/Δβ/direction/CpG-count/width |
| Visualise | `output$d_volcano`, `output$d_manhattan` | `mod_methyl_dmp_volcano()`, `mod_methyl_dmr_manhattan()` | Renders plots from `d_run()$df` (unfiltered) with threshold lines |
| Tabulate | `output$d_table` | `DT::datatable()` | Renders `d_filtered()` |
| Region-level CpG inspection | `d_region_selected()` reactive | `load_default_dmp("sva", sex)`, `methyl_champ_probe_positions()`, `merge()` | Recovers which per-CpG DMP rows fall inside the selected region's coordinates |
| Export | `d_download_*` handlers | `utils::write.csv()` | Writes complete / significant / filtered / configuration CSVs |

**No probe filtering, normalisation, model fitting, or region calling happens live in this tab.** Everything from "candidate region" onward was computed once, outside the app, by `script04_dmr_sexstratified/04_dmr_sexstratified.R` (documented in the bundled `METHODS_dmr_sexstratified.md`), and the app only ever filters/plots/exports that fixed table.

### DMR tab (live engine) — actual pipeline

| Stage | Object | Function | What changes |
|---|---|---|---|
| Load matrix + sheet | `methyl_dataset$beta`, `methyl_dataset$sample_sheet` | (from Dataset tab; out of scope) | — |
| Resolve sex column | `sex_col()` | `mod_methyl_dmp_sex_col()` | Auto-detects a `sex`/`Sex`/`gender`/`Gender` column |
| Resolve array annotation | `anno_result()` | `methyl_get_annotation()` | Loads chromosome/position/SNP manifest for 450K/EPIC only |
| Match samples | inside `live_result()` | `methyl_sheet_sample_ids()`, `intersect()` | Aligns matrix columns to sheet rows by sample ID |
| Subset by sex (optional) | `beta0`→subset | logical indexing on `sex_col` | Restricts to one sex if chosen |
| Subset by comparison groups | `beta1`, `grp` | logical indexing + `factor()` | Restricts to the two chosen group levels only |
| Subset by covariate completeness | `beta1`, `ph1` | `stats::complete.cases()` | Drops samples missing a selected covariate |
| Convert scale if needed | `beta_scale_full` | `2^beta1/(1+2^beta1)` (if M-value input) | Puts data on 0–1 beta scale for QC filters only |
| QC filter probes | `keep_probe` | `methyl_filter_missing()`, `methyl_filter_variance()`, `methyl_filter_snp()` (optional), position match | Drops probes failing missingness/variance/SNP/no-annotation criteria |
| Build M-values for the fit | `m` | `log2(beta/(1-beta))` (if beta input) | Converts filtered beta matrix to M-values for `limma` |
| Fit model | `fit`, `fit2`, `tt` | `stats::model.matrix()`, `methyl_chunked_lmfit()`, `limma::contrasts.fit()`, `limma::eBayes()`, `limma::topTable()` | Produces per-CpG t/p-value/logFC for the Group2-vs-Group1 contrast |
| Compute per-CpG Δβ | `dbeta` | `rowMeans()` on `beta_scale` per group | Region-level Δβ later derives from this |
| Build annotated GRanges | `gr`, `annot` | `GenomicRanges::GRanges()`, `stats::p.adjust(method="BH")`, `methods::new("CpGannotated", ...)` | Attaches chr/pos/stat/rawpval/diff/ind.fdr/is.sig per CpG |
| Call regions | `dmr_raw`, `ranges` | `DMRcate::dmrcate()`, `DMRcate::extractRanges()` | Produces one row per candidate genomic region |
| Region-level FDR | `dt$dmr_fdr` | `stats::p.adjust(dt$Stouffer, "BH")` | Genome-wide correction across all candidate regions in this run |
| Attach per-group means | `dt$ref_mean_beta`, `dt$comp_mean_beta` | `GenomicRanges::findOverlaps()`, `tapply()` | Region-level descriptive means, independent of DMRcate's own statistic |
| Filter for display | `live_filtered()` | `mod_methyl_dmr_filter()` | Adjustable post-hoc, without recalling regions |
| Visualise | `output$live_volcano/live_manhattan/live_effectplot/live_heatmap/live_qq` | `mod_methyl_dmp_volcano()`, `mod_methyl_dmr_manhattan()`, `mod_methyl_dmr_topplot()`, `mod_methyl_dmr_heatmap()`, `mod_methyl_qq_plot()` | Region- and CpG-level plots |
| Tabulate | `output$live_table` | `DT::datatable()` | Renders `live_filtered()` with a computed `significant` column |
| Region-level CpG inspection | `live_region_selected()` | logical range match, `mod_methyl_dmp_betadist()` | Per-sample β boxplot for the selected region's constituent CpGs |
| Export | `download_live_*` handlers | `utils::write.csv()` | Complete / significant / filtered / annotation-only / region-methylation-matrix / configuration CSVs |
| Publish summary | `methyl_results$dmr` | `observeEvent(live_result(), ...)` | Writes `comparison`, `n_regions`, `n_sig` for cross-module consumption (e.g. ArthOChat context) |

---

## 7. Tab 1: SVA

### Purpose
Reproduces, inside the app, the thesis's own published sex-stratified DMRcate region calling on the bundled GSE42861 dataset (Liu et al. 2013), letting the user filter and explore the fixed result rather than compute anything new.

### Inputs
- **Data source:** two bundled CSVs, `data/preloaded/methylomics/tables/script04_dmr_sexstratified/tables/dmr_{female,male}_full.csv`, read via `load_default_dmr(sex)` ([global.R:329-335](../../ArthOMix/global.R#L329-L335)). No upload path exists for this tab; it is gated entirely on `methyl_dataset$preloaded` being `TRUE` ([mod_methyl_dmr.R:230](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L230)).
- **Structure verified against the actual CSV file** (`dmr_female_full.csv`, 4,657 rows including header — matches `METHODS_dmr_sexstratified.md`'s reported 4,657 female candidate regions exactly): columns `seqnames, start, end, width, strand, no.cpgs, min_smoothed_fdr, Stouffer, HMFDR, Fisher, maxdiff, meandiff, overlapping.genes, sex, dmr_fdr` — one row per candidate genomic region, not per CpG.
- **User-supplied controls:** sex/stratum (radio), region-level FDR threshold, minimum |Δβ|, direction (any/hyper/hypo), minimum CpGs per region, minimum/maximum region width (bp), and a "Run DMR Analysis" button.
- **Required metadata:** none beyond the bundled CSVs; `load_default_meth_pheno()` is used only to report female/male sample counts in the descriptive banner text, not as an analysis input.

### Processing
No new statistics are computed. On clicking "Run DMR Analysis", `d_run()` selects the chosen sex's table, adds a `direction` column derived from the sign of `meandiff`, and freezes the current threshold inputs. `d_filtered()` then applies `mod_methyl_dmr_filter()` (FDR ≤ threshold, |Δβ| ≥ threshold, optional direction, CpG-count ≥ minimum, width within bounds) to produce the displayed/exported table.

### Functions Used
`load_default_dmr`, `load_default_meth_pheno` (data loading, global.R); `mod_methyl_dmr_filter`, `mod_methyl_dmp_filter` (filtering); `mod_methyl_dmp_volcano`, `mod_methyl_dmr_manhattan` (plotting); `DT::datatable`/`DT::formatSignif`/`DT::renderDataTable` (table rendering); `load_default_dmp`, `methyl_champ_probe_positions`, `merge()` (region-level CpG inspection); `utils::write.csv` (downloads).

### Outputs
- A volcano plot (mean Δβ vs. −log10 region FDR) and a Manhattan-style plot (genomic position vs. −log10 region FDR) of the *unfiltered* per-sex table, with threshold lines.
- Value boxes: candidate DMR count, count passing region FDR, count passing all current filters, sex label, hyper/hypomethylated counts (filtered).
- A sortable/searchable results table (`d_table`) of `d_filtered()`, showing `seqnames, start, end, width, no.cpgs, meandiff, maxdiff, Stouffer, dmr_fdr, overlapping.genes, direction`.
- A "Region-level inspection" panel: on selecting a table row, shows that region's constituent CpGs' own per-CpG `dbeta/p_bacon/fdr_bacon` from the DMP tab's default (SVA) analysis — not per-sample β values, because the precomputed pipeline does not bundle a per-sample matrix.
- Four CSV downloads: complete table, significant-only table, currently filtered table, and a one-row-per-parameter analysis-configuration record.

### Tab Connection
This tab is entirely self-contained: it neither reads from nor writes to the "DMR" tab. `methyl_results$dmr` (the one object shared with the rest of the app) is written **only** by the live "DMR" tab's `observeEvent(live_result(), ...)` ([mod_methyl_dmr.R:845-853](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L845-L853)) — running the SVA tab never updates `methyl_results$dmr`.

### Audit Notes
- **Confirmed implementation issue (reactivity):** after clicking "Run DMR Analysis" once, the FDR/Δβ/direction/CpG-count/width controls do **not** take effect until the button is clicked again — `d_filtered()` reads its thresholds from `d_run()$fdr/$dbeta/...`, which are captured only at click time ([mod_methyl_dmr.R:268-289](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L268-L289)), unlike the "DMR" tab where the equivalent controls are read live from `input$live_fdr` etc. every time ([mod_methyl_dmr.R:860-866](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L860-L866)). A user adjusting the FDR slider on the SVA tab and expecting the table/plots to update immediately will see no change until they press Run again.
- **Confirmed minor inconsistency (threshold boundary):** the value-box "significant" count ([mod_methyl_dmr.R:306](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L306)) and the "Significant DMRs" CSV download ([mod_methyl_dmr.R:377](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L377)) both use strict `dmr_fdr < r$fdr`, while the volcano/Manhattan plots' highlighted points and the main filtered table both use `dmr_fdr <= fdr_max` (via `mod_methyl_dmp_volcano`/`mod_methyl_dmr_manhattan`/`mod_methyl_dmp_filter`). A region whose `dmr_fdr` exactly equals the chosen threshold would be shown as significant on the plots/table but excluded from the "significant" count and download. This only manifests at an exact floating-point boundary and is unlikely to affect real analyses in practice.
- **No issue identified:** the tab is honestly scoped — it explicitly states in its own UI text that "nothing here recomputes anything," and the region-level CpG inspection panel correctly substitutes per-CpG statistics for the per-sample view it cannot provide, rather than fabricating or omitting the panel.

---

## 8. Tab 2: DMR

### Purpose
A fully configurable, live region-calling engine that runs the same `DMRcate::dmrcate()` algorithm used to build the bundled tables, but on whatever beta/M-value matrix and sample sheet is currently loaded (an upload, or — when this deployment has the raw matrix bundled — the preloaded dataset's own live matrix), for any two chosen group levels, optionally restricted to one sex, with the DMRcate tuning parameters exposed as controls.

### Inputs
- **Data source:** `methyl_dataset$beta` (probe × sample matrix) and `methyl_dataset$sample_sheet` (data frame), both populated by the Dataset tab (out of scope here). If `methyl_dataset$beta` is `NULL`, this tab shows an explanatory message instead of any controls ([mod_methyl_dmr.R:499-508](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L499-L508)).
- **Required sample-sheet content:** a group column with ≥2 distinct non-missing values (user-chosen, defaults to a `group`/`Group`/`disease`/`Disease` column if present, else the first column); optionally a sex column (auto-detected by name); optionally covariate columns (any column with ≥2 distinct values that isn't the group/sex/ID column, and isn't a unique-per-row character column).
- **Required genomic information:** chromosome + base-pair position for every tested CpG, from `methyl_get_annotation(methyl_dataset$array_type)` — **only available for 450K and EPIC(v1) arrays** ([annotation.R:18-21](../../ArthOMix/R/methylomics/annotation.R#L18-L21)); EPICv2/WGBS/RRBS/Custom array show an explicit "annotation unavailable" message and the analysis cannot run ([mod_methyl_dmr.R:529-530, 625-627](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L529-L530)).
- **User-supplied controls:** sex subset (if a sex column exists), group column, two group levels (Group 1/reference, Group 2/comparison), covariates (checkbox list), region-level FDR threshold, CpG seeding p-value, minimum |Δβ|, direction, minimum CpGs per region, minimum/maximum region width, minimum valid-sample percentage, minimum variance, SNP-probe removal (checkbox, default off), and — behind an "advanced" toggle — DMRcate's own `lambda`, `C`, and an optional manual candidate-region p-value cutoff override.
- **Accepted formats:** whatever the Dataset tab already parsed into `methyl_dataset$beta`/`$sample_sheet`/`$input_scale` (beta-scale or M-value-scale matrices are both handled, distinguished by `methyl_dataset$input_scale`).

### Processing
Full pipeline as detailed in §6 above: sample matching → optional sex subset → group subset → covariate-completeness subset → scale conversion for QC purposes only → missingness/variance/SNP/position probe filters → M-value conversion for the fit → design-matrix construction (with a rank check) → chunked `limma` fit → contrast + eBayes → per-CpG Δβ from group means → construction of a `CpGannotated` GRanges object → `DMRcate::dmrcate()` → `DMRcate::extractRanges()` → region-level BH-FDR on the Stouffer statistic → attaching per-group mean β per region.

### Functions Used
`mod_methyl_dmp_sex_col`, `mod_methyl_dmp_sex_choices`, `mod_methyl_dmp_covariate_cols`, `methyl_get_annotation`, `methyl_norm_status`, `methyl_sheet_sample_ids`, `methyl_filter_missing`, `methyl_filter_variance`, `methyl_filter_snp`, `stats::model.matrix`, `stats::complete.cases`, `methyl_chunked_lmfit` → `limma::lmFit`, `limma::makeContrasts`, `limma::contrasts.fit`, `limma::eBayes`, `limma::topTable`, `stats::p.adjust`, `GenomicRanges::GRanges`, `IRanges::IRanges`, `methods::new`, `DMRcate::dmrcate`, `DMRcate::extractRanges`, `GenomicRanges::findOverlaps`, `S4Vectors::queryHits`/`subjectHits`, `tapply`, `mod_methyl_lambda_gc`, `mod_methyl_qq_plot`, `mod_methyl_dmp_volcano`, `mod_methyl_dmr_manhattan`, `mod_methyl_dmr_topplot`, `mod_methyl_dmr_heatmap`, `mod_methyl_dmp_betadist`, `DT::datatable`, `utils::write.csv`.

### Outputs
- A configuration/sample-size summary card, a genomic-inflation-factor (λ) diagnostic with a QQ plot and an explicit warning when λ>1.1 that "this live engine does not apply SVA/bacon correction," value boxes (candidate/significant/hyper/hypo counts), a volcano plot, a Manhattan plot, a top-N Δβ bar chart (selectable N), a region×sample mean-methylation heatmap (capped at 50 rendered regions, all regions still in the table/exports), a results table with a computed `significant` Yes/No column, and a region-level inspection panel with a per-CpG boxplot (real per-sample β values, since a live matrix is available) plus a per-sample table.
- Six CSV downloads: complete, significant, filtered, annotation-only, a region×sample mean-methylation matrix, and a configuration record covering every parameter that reached the analysis (group/sex/covariates, seeding p-value, lambda, C, min.cpgs, cutoff, thresholds, sample counts, timestamp).
- `methyl_results$dmr <- list(comparison, n_regions, n_sig)`, written on every successful run, for consumption elsewhere in the app (e.g. the ArthOChat assistant's context).

### Tab Connection
Independent of the SVA tab in every respect except shared helper functions and a shared visual style. Does not read the SVA tab's `d_run()`/`d_filtered()` state, and the SVA tab does not read this tab's `live_result()`.

### Audit Notes
- **Confirmed implementation issue (mislabeled field):** the returned analysis object sets `n_cpgs_before_filter = nrow(beta1)` ([mod_methyl_dmr.R:832](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L832)), but `beta1` was already reassigned to the **post**-QC-filter subset two statements earlier (`beta1 <- beta1[keep_probe, , drop = FALSE]`, [mod_methyl_dmr.R:718](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L718)) before that final list is built. Because `m`/`tt` are derived from this same filtered `beta1`, `nrow(beta1)` at that point is numerically identical to `n_cpgs_tested`. The UI text "`%s CpGs tested after QC filters: %s of %s`" ([mod_methyl_dmr.R:894-897](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L894-L897)) is therefore always rendered as "X of X" (100%), never showing the true pre-filter probe count the label implies. This is a genuine reporting bug — it does not affect the statistical analysis itself (the correct filtered matrix is still used for the fit), only the displayed "before/after" QC summary.
- **Potential limitation (not a bug):** when "All samples" (`__all__`) is selected as the sex subset, sex is **not automatically added** as a covariate; the sex column simply remains available in the optional covariate checklist for the user to add manually ([mod_methyl_dmr.R:607-608](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L607-L608)). A combined-sex comparison run without manually adding sex as a covariate risks confounding a sex-linked methylation pattern with the group effect being tested. The code does not force this covariate, nor does it warn about its absence specifically (the genomic-inflation warning is a general, not sex-specific, diagnostic).
- **No issue identified:** the design-matrix rank check (`qr(design)$rank == ncol(design)`) correctly precedes the model fit and gives an informative message rather than letting `limma` fail opaquely or silently drop a coefficient.
- **No issue identified:** UI parameters that scientifically should reach `DMRcate::dmrcate()` do reach it verbatim — `lambda`, `C`, `min.cpgs`, and `pcutoff` are all read from their respective `input$live_*` values (with the "override" checkbox correctly gating whether a manual `pcutoff` or DMRcate's own `"fdr"` default is used) and passed directly into the `dmrcate()` call ([mod_methyl_dmr.R:783-789](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L783-L789)).

---

## 9. Complete Function Audit

Functions are grouped by category. "Source" identifies whether a function is custom (defined in this project), a base-R/`stats`/`utils` function, or a named package function.

### 9.1 Data input / validation

| Function | Source | Used in | Inputs | What it does | Audit |
|---|---|---|---|---|---|
| `load_default_dmr(sex)` | Custom, [global.R:329](../../ArthOMix/global.R#L329) | `default_data()` | `sex` ("female"/"male") | Reads the precomputed `dmr_{sex}_full.csv` via `data.table::fread`, returns `NULL` if the deployment lacks the bundled data folder | Correct: matches on-disk CSV exactly (verified: female file has 4,657 data rows, matching the methods document's reported count) |
| `load_default_meth_pheno()` | Custom, [global.R:340](../../ArthOMix/global.R#L340) | `default_data()` | none | Reads the QC-derived sample-metadata CSV, used only to compute descriptive female/male counts shown in a banner | Correct, non-critical (display only) |
| `load_default_dmp(stage, sex)` | Custom, [global.R:313](../../ArthOMix/global.R#L313) | `d_region_selected()` | `"sva"`, sex | Reads the SVA-adjusted per-CpG DMP table for a stratum | Correct; column names (`cpg,dbeta,p_bacon,fdr_bacon`) match what `d_region_table` selects |
| `methyl_get_annotation(array_type)` | Custom, [annotation.R:48](../../ArthOMix/R/methylomics/annotation.R#L48) | `anno_result()` | array type string | Loads/caches chromosome+position+SNP manifest from the array's Bioconductor annotation package, avoiding `minfi::getAnnotation()`'s namespace-attach side effect | Correct and deliberately defensive (documented reason: avoids `Biostrings` masking `strsplit()` app-wide); correctly returns `ok=FALSE` for array types with no bundled manifest package |
| `methyl_sheet_sample_ids(sheet, all_ids)` | Custom, [qc.R:456](../../ArthOMix/R/methylomics/qc.R#L456) | `live_result()` | sample sheet, matrix column names | Resolves the sheet's sample-ID column (or falls back to row order/rownames) so matrix columns can be matched to sheet rows | Correct; explicit id-column priority avoids blind positional matching except as a documented last resort |
| `stats::complete.cases(cc)` | base R (`stats`) | `live_result()` | covariate columns | Flags samples with no missing covariate values | Correctly used before model fitting |

### 9.2 QC / probe filtering (methylomics-shared, reused as-is)

| Function | Source | Used in | Inputs | What it does | Scientific significance | Audit |
|---|---|---|---|---|---|---|
| `methyl_filter_missing(mat, max_na_frac)` | Custom, [qc.R:30](../../ArthOMix/R/methylomics/qc.R#L30) | `live_result()` | beta-scale matrix, max allowed missing fraction | Drops probes exceeding a per-row missingness fraction | Removes unreliable probes before fitting | Correctly derives `max_na_frac` from the UI's "minimum valid %" ([mod_methyl_dmr.R:696](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L696)) |
| `methyl_filter_variance(mat, min_variance)` | Custom, [qc.R:36](../../ArthOMix/R/methylomics/qc.R#L36) | `live_result()` | beta-scale matrix, min variance | Drops near-invariant probes (default threshold 0, i.e. off unless the user sets one) | Near-zero-variance probes carry no group-discriminating information | Correct; default 0 means this filter is opt-in via the threshold rather than a checkbox — reasonable |
| `methyl_filter_snp(mat, anno_result)` | Custom, [qc.R:68](../../ArthOMix/R/methylomics/qc.R#L68) | `live_result()` (optional, checkbox default off) | beta-scale matrix, annotation object | Flags probes overlapping a known SNP via manifest `Probe_rs/CpG_rs/SBE_rs` columns | SNP-overlapping probes can show apparent "methylation" differences driven by genotype, not epigenetics | Correct; gracefully no-ops (keeps everything) when annotation is unavailable |
| Position match (`match(rownames(beta_scale_full), rownames(a))`) | Inline, [mod_methyl_dmr.R:712-714](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L712-L714) | `live_result()` | filtered matrix rownames, annotation data frame | Drops any probe absent from the manifest or missing chr/pos | Region calling is impossible without genomic coordinates | Correct and necessary — DMRcate/GRanges cannot place an unannotated probe |

### 9.3 Methylation / statistical processing

| Function | Source | Used in | Inputs | What it does | Scientific significance | Audit |
|---|---|---|---|---|---|---|
| `stats::model.matrix(~0+grp)` / `~covariates` | base R | `live_result()` | group factor, covariate data frame | Builds a no-intercept group design matrix, optionally column-bound with covariate terms (covariate names backtick-quoted for safety) | Standard `limma` design-matrix construction for a two-group moderated t-test with adjustment | Correct; explicitly checks `qr(design)$rank == ncol(design)` before fitting and gives an informative error on rank deficiency (e.g. a covariate collinear with group) |
| `methyl_chunked_lmfit(m, design, chunk_size=20000)` | Custom, [mod_methyl_dmp.R:143](../../ArthOMix/R/methylomics/mod_methyl_dmp.R#L143) | `live_result()` | M-value matrix, design matrix | Row-chunks `limma::lmFit()` calls and re-concatenates their per-gene summary fields into one `MArrayLM` object | Avoids exceeding R's vector-memory limit on a full ~400k-probe × hundreds-of-samples fit (documented as a reproduced real crash on the bundled 689-sample matrix at full size) | Code comment documents a bit-for-bit verification against a whole-matrix fit on synthetic data; **no issue identified** given that documented verification |
| `limma::makeContrasts`, `limma::contrasts.fit`, `limma::eBayes`, `limma::topTable` | Package (`limma`) | `live_result()` | fitted model, contrast string, design | Standard moderated-t-test contrast pipeline: builds the Group2-vs-Group1 contrast, refits, empirical-Bayes-moderates variance, extracts the full per-CpG results table (`sort.by="none"` preserves row order for downstream alignment) | This is the calibrated per-CpG statistical engine feeding region calling | Correct; `sort.by="none"` is the right choice here since row order must stay aligned with `chr`/`pos`/`dbeta` vectors built by position, not by name-matching afterwards |
| `rowMeans(beta_scale[, grp==level])` | base R | `live_result()` | filtered beta-scale matrix, group factor | Per-CpG mean β within each group | Basis for per-CpG and (after region aggregation) per-region Δβ | Correct |
| `stats::p.adjust(tt$P.Value, "BH")` (`ind.fdr`) | base R (`stats`) | `live_result()` | raw per-CpG p-values | Genome-wide BH correction across all tested CpGs | This is the individual-CpG FDR retained in the `CpGannotated` object, distinct from the region-level FDR computed later | Correct; matches the documented DMRcate convention of keeping `ind.fdr` genome-wide-corrected even though seeding uses the uncorrected p-value |
| `mod_methyl_lambda_gc(p)` | Custom, [mod_methyl_dmp.R:97](../../ArthOMix/R/methylomics/mod_methyl_dmp.R#L97) | `live_result()` | raw per-CpG p-values (`tt$P.Value`, pre-seeding) | Median-χ² genomic-inflation factor (λ) | Standard EWAS/GWAS calibration diagnostic; here it substitutes for the SVA/bacon correction the live engine deliberately does not apply | Correct; the code comment explicitly requires the *full, unfiltered* p-vector, and the call site does pass `tt$P.Value` before any region-seeding subset |
| `mod_methyl_qq_plot(p)` | Custom, [mod_methyl_dmp.R:106](../../ArthOMix/R/methylomics/mod_methyl_dmp.R#L106) | `output$live_qq` | same raw p-values | Observed-vs-expected −log10(p) QQ plot | Visual companion to λ | Correct |

### 9.4 Region identification (DMRcate / GenomicRanges core)

| Function | Source | Used in | Inputs | What it does | Scientific significance | Audit |
|---|---|---|---|---|---|---|
| `GenomicRanges::GRanges(...)` | Package (`GenomicRanges`) | `live_result()`, `mod_methyl_dmr_heatmap`, download handlers | chromosome, `IRanges::IRanges(pos,pos)`, plus per-CpG metadata columns (`stat, rawpval, diff, ind.fdr, is.sig`) | Builds a genomic-coordinate object annotated with the per-CpG limma statistics | This *is* the `CpGannotated`-shaped input DMRcate needs | Correct; constructed directly rather than via `DMRcate::cpg.annotate()`'s own internal refit, exactly matching the bundled pipeline's own documented approach (`METHODS_dmr_sexstratified.md` §2.BB.2) |
| `methods::new("CpGannotated", ranges = gr)` | base R (`methods`) + DMRcate's S4 class | `live_result()` | the `GRanges` object above | Wraps it in DMRcate's expected S4 class | Required so `DMRcate::dmrcate()` accepts the object as valid input | Correct; code comment explains the preceding `requireNamespace("DMRcate")` call is required first so the S4 class is registered before `methods::new()` is invoked |
| `DMRcate::dmrcate(annot, lambda, C, min.cpgs, pcutoff)` | Package (`DMRcate`) | `live_result()` | annotated object + tuning parameters | Gaussian-kernel-smooths the per-CpG statistics along the genome and proposes candidate regions from CpGs whose smoothed statistic crosses a threshold | This is the actual DMR-calling algorithm — the statistical core of the entire module | Correct: every exposed UI parameter (`lambda`, `C`, `min.cpgs`, `pcutoff`) is passed through unmodified; `min.cpgs` is additionally floored at 2 (`max(2, ...)`), a sensible guard since a "region" of one CpG is not a region |
| `DMRcate::extractRanges(dmr_raw, genome="hg19")` | Package (`DMRcate`) | `live_result()` | DMRcate's internal result object | Converts to a `GRanges` of final candidate regions with `no.cpgs`, `Stouffer`, `meandiff`, `maxdiff`, `overlapping.genes`, etc. | Produces the genomic-coordinate table the rest of the tab consumes | Correct; `genome="hg19"` matches the annotation build used elsewhere in the app (450K/EPIC v1 annotation packages are both hg19-based) |
| `GenomicRanges::findOverlaps(probe_gr, region_gr)` + `S4Vectors::queryHits`/`subjectHits` | Package | `live_result()`, `mod_methyl_dmr_heatmap`, `download_live_matrix` | per-probe GRanges, per-region GRanges | Vectorized join: which probes fall inside which called region | Needed to compute region-level per-group means and per-region-per-sample matrices without an explicit loop over base-pair ranges | Correct; the same idiom is documented as matching the bundled pipeline's own biomarker-panel construction |

### 9.5 Filtering, ranking, and multiple-testing correction (module-specific)

| Function | Source | Used in | Inputs | What it does | Audit |
|---|---|---|---|---|---|
| `mod_methyl_dmr_filter(df, fdr_col, effect_col, fdr_max, effect_min, direction, min_cpgs, min_width, max_width)` | Custom, [mod_methyl_dmr.R:75](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L75) | both tabs | region table + thresholds | Delegates FDR/Δβ/direction filtering to `mod_methyl_dmp_filter()`, then applies region-specific CpG-count and width bounds | Correct reuse; conditionally checks `"no.cpgs"`/`"width"` are present before filtering on them, so it degrades gracefully if a table lacks those columns |
| `mod_methyl_dmp_filter(df, fdr_col, effect_col, fdr_max, effect_min, direction)` | Custom, [mod_methyl_dmp.R:71](../../ArthOMix/R/methylomics/mod_methyl_dmp.R#L71) | via the above | any table with a named FDR + effect column | Generic `fdr_col <= fdr_max & abs(effect_col) >= effect_min`, optionally signed by direction | Correct, generic, reused verbatim from the DMP module rather than reimplemented |
| `stats::p.adjust(dt$Stouffer, method="BH")` (`dt$dmr_fdr`) | base R | `live_result()` | DMRcate's own Stouffer statistic per candidate region | The module's own region-level multiple-testing correction, treating each candidate region as one test | This is the correct unit of correction for a "how many regions are significant" question, and mirrors the bundled precomputed pipeline's own documented step exactly |

### 9.6 Annotation

| Function | Source | Used in | Inputs | What it does | Audit |
|---|---|---|---|---|---|
| `methyl_champ_probe_positions()` | Custom, [mod_methyl_dmr.R:189](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L189) | `d_region_selected()` (SVA tab only) | none (lazy-loads `ChAMPdata::probe.features`) | Builds a `cpg/chr/pos` lookup table from `ChAMPdata`, cached for the process lifetime | Correct: same annotation source the original precomputed pipeline used (per `METHODS_dmr_sexstratified.md`), so region-to-CpG coordinate matching for the SVA tab is internally consistent with how the regions were originally called |
| `overlapping.genes` (from `DMRcate::extractRanges`) | Package-derived field | both tabs | — | DMRcate's own nearest/overlapping-gene annotation, attached during region extraction | This is DMRcate's built-in annotation, not a separate annotation step this module performs |

### 9.7 Visualization

| Function | Source | Used in | Inputs | What it does | Audit |
|---|---|---|---|---|---|
| `mod_methyl_dmr_manhattan(dt, fdr_max)` | Custom, [mod_methyl_dmr.R:87](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L87) | both tabs | region table | Genomic-position-ordered −log10(region FDR) scatter, autosomes only (`chr_num` must parse to an integer, so X/Y/M regions are silently excluded from this plot only) | Correct; `validate(need(...))` gives an informative message rather than an opaque `ggplot` error if no autosomal region has both fields |
| `mod_methyl_dmr_topplot(dt, n)` | Custom, [mod_methyl_dmr.R:110](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L110) | DMR tab only | region table, N | Top-N-by-FDR horizontal Δβ bar chart, labeled by coordinates + gene | Correct; not present in the SVA tab (a deliberate scope difference, not a bug) |
| `mod_methyl_dmr_heatmap(sig_dt, beta_scale, probe_chr, probe_pos, grp, max_regions)` | Custom, [mod_methyl_dmr.R:134](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L134) | DMR tab only | significant regions, beta matrix, group factor | Region×sample mean-β heatmap, samples grouped/faceted, capped at 50 rendered regions | Correct; explicitly documents the cap is rendering-only (table/exports are uncapped) and the UI text states this to the user |
| `mod_methyl_dmp_volcano(df, effect_col, p_col, ..., fdr_max, effect_min)` | Custom, [mod_methyl_dmp.R:78](../../ArthOMix/R/methylomics/mod_methyl_dmp.R#L78) | both tabs | region table (reused unmodified from the DMP module, since it only needs a named effect + p/FDR column) | Effect-size vs. −log10(p) scatter with significance colouring | Correct reuse |
| `mod_methyl_dmp_betadist(beta_mat, cpgs, grp)` | Custom, [mod_methyl_dmp.R:245](../../ArthOMix/R/methylomics/mod_methyl_dmp.R#L245) | DMR tab, region inspection only | per-CpG beta submatrix, CpG IDs, group factor | Per-CpG boxplot by group, for the CpGs inside a selected region | Correct; only usable in the DMR tab because only there is a real per-sample matrix available |

### 9.8 Shiny / reactive framework

| Function | Source | Used in | Role |
|---|---|---|---|
| `moduleServer`, `NS`, `reactive`, `reactiveVal`, `eventReactive`, `observeEvent`, `renderUI`, `renderPlot`, `DT::renderDataTable`, `downloadHandler`, `validate`/`need`, `req`, `withProgress`, `withSpinner`, `outputOptions(..., suspendWhenHidden=FALSE)` | Shiny / `shinycssloaders` / `DT` | throughout | Standard Shiny reactive-programming primitives; see §17 for the specific reactive graph and gating logic |

### 9.9 File / download

| Function | Source | Used in | Purpose |
|---|---|---|---|
| `utils::write.csv(df, file, row.names=FALSE)` | base R (`utils`) | every `download*` handler | Writes the exact in-memory data frame (complete/significant/filtered/annotation/matrix/configuration) to the downloaded file, with no additional transformation — what is downloaded is exactly what was computed, not a re-derived summary |

---

## 10. Data Structures and Objects

| Object | Produced by | Structure | Consumed by |
|---|---|---|---|
| `default_data()` | SVA tab load | `list(pheno, dmr_f, dmr_m)` | `d_run()` |
| `d_run()` | SVA tab click | `list(df, sex, fdr, dbeta, min_cpgs, min_width, max_width, direction)` | `d_filtered()`, `d_valueboxes_ui`, `d_manhattan`, `d_volcano`, downloads |
| `d_filtered()` | SVA tab | region data frame, row-subset of `d_run()$df` | `d_table`, `d_region_selected()`, `d_download_filtered` |
| `d_region_selected()` | SVA tab, row click | `list(row, cpgs)` — `cpgs` is a merged DMP×position data frame restricted to the selected region's coordinates | `d_region_summary`, `d_region_table` |
| `live_result()` | DMR tab click | `list(dt, ref, comp, group_col, sex_label, sex_col, covariates, design_formula, n_ref, n_comp, n_cpgs_tested, n_cpgs_before_filter, n_seed_cpgs, seed_p, lambda, C, min_cpgs, pcutoff, min_valid_pct, min_variance, snp_filter, snp_note, missing_note, variance_note, norm_status, dataset_source, preloaded, beta_scale, grp, probe_chr, probe_pos, lambda_gc, cpg_p_raw, run_at)` | virtually every other reactive/output in the DMR tab |
| `live_result()$dt` | DMRcate → `extractRanges` → this app's own additions | columns `seqnames, start, end, width, strand, no.cpgs, min_smoothed_fdr, Stouffer, HMFDR, Fisher, maxdiff, meandiff, overlapping.genes, dmr_fdr, direction, dmr_id, ref_mean_beta, comp_mean_beta` (`dmr_fdr`, `direction`, `dmr_id`, `ref_mean_beta`, `comp_mean_beta` are added by this module; the rest come directly from DMRcate) | `live_filtered()`, `live_sig()`, all live-tab outputs |
| `live_filtered()` / `live_sig()` | DMR tab | row-subsets of `live_result()$dt` | table, heatmap, downloads |
| `live_region_selected()` | DMR tab, row click | `list(row, cpg, beta, grp, ref, comp)` | `live_region_summary`, `live_region_plot`, `live_region_table` |
| `methyl_dataset` | passed in (Dataset tab) | reactiveValues incl. `$preloaded, $beta, $sample_sheet, $array_type, $source, $input_scale, $rg_set` | read throughout the DMR tab; never written by this module |
| `methyl_results$dmr` | written by DMR tab only | `list(comparison, n_regions, n_sig)` | consumed elsewhere in the app (out of scope for this document) |

---

## 11. Input Data Audit

### SVA tab
- **Type:** precomputed region-level results table (not a beta/M-value matrix). No CpG-level or genomic-annotation upload is possible here.
- **Structure:** verified against the actual bundled files — one row per candidate region, 15 columns, ~4,600 (female) / ~3,400 (male) rows.
- **Validation implemented:** `req(isTRUE(methyl_dataset$preloaded))` gates rendering entirely; `validate(need(!is.null(d$dmr_f) || !is.null(d$dmr_m), ...))` gives an explicit message if the bundled data folder is missing; `d_run()`'s own `validate(need(!is.null(df) && nrow(df) > 0, ...))` checks the chosen stratum's table is non-empty.
- **No validation exists for:** malformed/corrupted CSV content (column presence is assumed, not checked) — if the bundled CSV were ever edited to drop a required column, the failure would surface downstream as a `NULL`-column subscript error inside `mod_methyl_dmr_manhattan`/`mod_methyl_dmr_filter` rather than a clear upstream message. This is a **potential limitation**, not a confirmed issue, since the files are bundled and versioned with the app rather than user-supplied.

### DMR tab
- **Type:** beta-value or M-value matrix (probes × samples), distinguished by `methyl_dataset$input_scale`; sample-sheet data frame with a group column and, optionally, sex/covariate columns.
- **Expected structure:** matrix column names are sample IDs; sample-sheet rows are matched to matrix columns via `methyl_sheet_sample_ids()` (by an ID column if present, else row order if counts match, else `rownames(sheet)`).
- **Validation implemented (in order):** dataset loaded (`req`/`validate`); sample sheet present; group column selected and present in the sheet; two distinct group levels chosen and different from each other; array annotation available for this array type; ≥6 samples match between matrix and sheet; if sex-subsetting, ≥6 samples remain in that sex; ≥6 samples remain after group restriction is implicit via the next check; if covariates chosen, ≥6 samples have complete covariate data; each of the two groups has ≥3 samples; design matrix is full rank; ≥20 CpGs remain after missingness/variance/SNP/position filters; ≥2 CpGs pass the seeding p-value; ≥1 candidate region returned by DMRcate.
- **What happens with missing/invalid input:** every one of the above is a `validate(need(...))` call, which — per Shiny's own contract — stops execution and shows the message in place of output, rather than emitting a partial or silently empty result. This is a **correct, comprehensive validation chain**; no silent-empty-result path was found for the primary analysis trigger.
- **Not validated:** the code does not check that the group column's two chosen levels are not simply relabelled duplicates of each other in disguise (e.g. `"RA"` vs `"ra"` would be treated as different groups since matching is exact-string, case-sensitive) — a **potential limitation**, not confirmed to cause any actual error, just a possible silent misconfiguration if a sample sheet has inconsistent casing.

---

## 12. Statistical Methodology Audit

| Parameter | Meaning | Source | Default | User-adjustable? | Effect of changing | Actually reaches the analysis function? |
|---|---|---|---|---|---|---|
| Region-level FDR threshold | Cutoff on `dmr_fdr` for "significant" | UI (`d_fdr`/`live_fdr`) | 0.05 | Yes | Changes which regions are flagged/exported as significant; does **not** change which regions DMRcate calls | Used only in display/filter/download logic, by design (it is a post-hoc threshold, not a DMRcate input) |
| Minimum \|Δβ\| | Effect-size floor | UI (`d_dbeta`/`live_dbeta`) | 0 (SVA) / 0.05 (DMR) | Yes | Same as above | Display/filter only, by design |
| Direction | Restrict to hyper/hypo | UI | "any" | Yes | Same as above | Display/filter only, by design |
| Minimum CpGs per region | Display floor on `no.cpgs` | UI (`d_mincpgs`/`live_mincpgs`) | 3 | Yes | SVA tab: display filter only. **DMR tab: this same value is also passed as DMRcate's own `min.cpgs` argument** ([mod_methyl_dmr.R:785, 789](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L785)), floored at 2 | Yes — DMR tab: changes which regions DMRcate calls at all (only takes effect on next "Run"); SVA tab: display-only, correctly labeled "(display filters)" in that tab's UI |
| Min/max region width (bp) | Display bounds on `width` | UI | 0/0 (no limit) | Yes | Display/filter only in both tabs — DMRcate itself has no width parameter | Correctly never passed to DMRcate; used only in `mod_methyl_dmr_filter` |
| CpG seeding p-value | Raw per-CpG p-value threshold that flags `is.sig` | UI (`live_seed_p`), DMR tab only | 0.05 | Yes | Controls which CpGs can seed a candidate region; lower = fewer, more stringent seeds | Reaches the analysis directly, via `gr$is.sig <- tt$P.Value < seed_p` before the `CpGannotated` object is built |
| Kernel bandwidth (`lambda`) | DMRcate Gaussian kernel width, nt | UI (`live_lambda`, behind "advanced" toggle) | 1000 | Yes | Wider bandwidth pools more distant CpGs into one region | Passed directly to `DMRcate::dmrcate(lambda=...)` |
| Bandwidth scaling (`C`) | DMRcate scaling factor | UI (`live_c`) | 2 | Yes | Controls kernel sharpness | Passed directly to `DMRcate::dmrcate(C=...)` |
| Candidate-region p cutoff (`pcutoff`) | DMRcate's own internal region-candidacy cutoff | UI (`live_pcutoff`, gated behind an explicit "override" checkbox) | DMRcate's own `"fdr"` default unless overridden | Yes (opt-in only) | Changes how permissive DMRcate is when nominating candidate regions | Passed directly; correctly defaults to DMRcate's own recommended `"fdr"` behavior unless the user explicitly opts into a manual value, and the UI itself labels the override "not recommended" |
| Minimum valid sample % | Probe missingness floor | UI (`live_min_valid_pct`) | 80 | Yes | Converted to `max_na_frac` and passed to `methyl_filter_missing()` | Yes |
| Minimum variance | Probe variance floor | UI (`live_min_variance`) | 0 (off) | Yes | Passed to `methyl_filter_variance()` | Yes |
| SNP-probe removal | Checkbox | UI (`live_snp_filter`) | FALSE | Yes | Removes probes overlapping a known SNP, via `methyl_filter_snp()` | Yes, but only when checked |
| Covariates | Additional design-matrix terms | UI (`live_covariates`) | none selected | Yes | Adjusts the limma fit for the selected phenotype columns | Yes, built into the design matrix with a rank check |

**Statistical test:** a moderated t-test per CpG (`limma::eBayes`) for a two-level group contrast, optionally covariate-adjusted, feeding DMRcate's Gaussian-kernel-smoothed region-calling algorithm. **Multiple-testing correction:** two distinct, correctly separated levels — (1) genome-wide BH correction of the raw per-CpG p-values into `ind.fdr`, retained in the `CpGannotated` object but not used to decide significance downstream, and (2) a second, independent BH correction of DMRcate's Stouffer combined-region statistic across all candidate regions in the run, producing `dmr_fdr`, which is what the module actually uses to decide "significant." **Region definition:** a DMRcate-called kernel-smoothed cluster of CpGs, not a fixed genomic window and not a post-hoc grouping of independently significant CpGs.

---

## 13. Region Identification Method

Both tabs use **`DMRcate::dmrcate()`** exclusively; no other region-calling method (bumphunter, comb-p, methylKit's region-level test, a fixed sliding window, etc.) is implemented anywhere in this module, despite those algorithm names appearing as **search aliases only** in the app's unrelated glossary/help feature (`global.R`, help registry entry `id="dmr"`, `aliases = c("dmr", "differentially methylated region", "dmrcate", "comb-p", "bumphunter")` — this list exists so a user typing "bumphunter" into the in-app help search is shown the DMRcate-based glossary entry; it is not a claim that bumphunter or comb-p are implemented).

DMRcate's algorithm: per-CpG test statistics are smoothed along genomic position using a Gaussian kernel of the configured bandwidth (`lambda`) and scaling (`C`); CpGs whose raw p-value passes the seeding threshold (`is.sig`) can start a candidate region; contiguous runs of CpGs meeting the (possibly smoothed) candidacy criterion, subject to `min.cpgs`, are merged into one region. `DMRcate::extractRanges()` then reports each region's boundary coordinates, CpG count, several combined-probability statistics (Stouffer, harmonic-mean FDR, Fisher), the maximum and mean per-CpG Δβ across the region's constituent CpGs, and any directly overlapping gene symbol(s).

The bundled SVA tab's regions were called this exact way, offline, with `lambda=1000, C=2` (per `METHODS_dmr_sexstratified.md` §2.BB.3) on hg19 coordinates. The live DMR tab exposes the same two parameters as user controls with the same defaults, so a user who does not change them and runs on comparable data is running the same algorithm configuration as the bundled analysis — though on their own matrix/groups/covariates, not a recomputation of the bundled result.

---

## 14. Filtering and Multiple-Testing Correction

Region-level FDR correction (`stats::p.adjust(Stouffer, "BH")`) is the only multiple-testing correction this module itself performs. It is applied once per run, across all candidate regions returned for that stratum/comparison — a design choice matching the bundled reference pipeline's own documented approach. There is no correction across multiple runs (e.g. running the DMR tab twice with different comparisons does not adjust across both runs' combined region counts — each run's `dmr_fdr` is computed independently within that run only, which is the scientifically expected behavior for one hypothesis test per run).

Non-statistical filters (minimum CpGs, minimum/maximum width) are display/export filters, not corrections — they change which rows are *shown*, not the p-values or FDR values themselves, except for the DMR tab's "minimum CpGs" control, which (as noted in §12) doubles as an actual DMRcate `min.cpgs` parameter and therefore *does* change which regions exist to be filtered, not merely which are displayed.

---

## 15. Annotation and Biological Interpretation

**Gene annotation:** `overlapping.genes`, produced entirely by `DMRcate::extractRanges()` — genes whose genomic span directly overlaps the called region's coordinates. This is a direct-overlap annotation, not a nearest-gene, regulatory-element, or functional annotation. A region with no `overlapping.genes` entry is not necessarily non-functional; it may simply fall in an intergenic interval.

**What a significant DMR legitimately supports, given this implementation:** a region-level statistical statement — that, under the fitted model (limma moderated t-test, optionally covariate-adjusted, on whichever matrix/probe-QC/groups were configured), several physically adjacent, DMRcate-selected CpGs show a directionally consistent methylation difference between the two compared groups, at a region-level FDR below the chosen threshold.

**What this implementation does not support, and the document does not claim:**
- **Causality.** No code path in this module performs or implies a causal test (that is the stated purpose of the app's separate Mendelian randomisation module, out of scope here). A significant DMR here is an association, not a claim about direction of causation between methylation and disease.
- **Cell-type-specific interpretation.** Nothing in this module adjusts for or reports cell-type composition (that is the Cell-Type Deconvolution submodule, out of scope). A DMR detected in a mixed cell population (e.g. whole blood, the bundled dataset) may reflect a shift in cell-type proportions between groups rather than a true within-cell-type methylation change, and the module's own outputs do not distinguish these possibilities.
- **Functional consequence.** An `overlapping.genes` hit is a coordinate overlap, not evidence that the gene's expression, splicing, or function is altered.
- **Independent validation of DMRcate's `min_smoothed_fdr`/`HMFDR`/`Fisher` columns.** These pass through unmodified from DMRcate; this module documents but does not itself re-derive or re-validate them.

---

## 16. Visualization and Outputs

| Output | Source function | Input | Filtered or raw? | Downloadable? |
|---|---|---|---|---|
| SVA volcano (`d_volcano`) | `mod_methyl_dmp_volcano` | `d_run()$df` | Raw (unfiltered) table, thresholds only shown as colour/lines | No (plot only) |
| SVA Manhattan (`d_manhattan`) | `mod_methyl_dmr_manhattan` | `d_run()$df` | Raw | No |
| SVA results table (`d_table`) | `DT::datatable` | `d_filtered()` | Filtered | Indirectly, via the separate download buttons on the same card |
| SVA region-CpG table (`d_region_table`) | `DT::datatable` | `d_region_selected()$cpgs` | Filtered to selected region's coordinates | No |
| SVA downloads | `write.csv` | `d_run()$df` / significant subset / `d_filtered()` / a config data frame | complete / significant / filtered / configuration | Yes, 4 files |
| DMR QQ plot (`live_qq`) | `mod_methyl_qq_plot` | `cpg_p_raw` (full, pre-seeding) | Raw | No |
| DMR volcano/Manhattan | `mod_methyl_dmp_volcano`/`mod_methyl_dmr_manhattan` | `live_result()$dt` | Raw | No |
| DMR top-N plot (`live_effectplot`) | `mod_methyl_dmr_topplot` | `live_result()$dt` | Raw, ranked | No |
| DMR heatmap (`live_heatmap`) | `mod_methyl_dmr_heatmap` | `live_sig()` (FDR+Δβ-filtered) | Filtered, capped at 50 rendered rows | Indirectly — full region×sample matrix is a separate download |
| DMR results table (`live_table`) | `DT::datatable` | `live_filtered()` + computed `significant` column | Filtered, with a boolean significance flag | Yes, via separate buttons |
| DMR region inspection plot/table | `mod_methyl_dmp_betadist` / `DT::datatable` | `live_region_selected()` | Restricted to selected region's CpGs, real per-sample values | No |
| DMR downloads | `write.csv` | `live_result()$dt` / `live_sig()` / `live_filtered()` / annotation columns / a region×sample matrix / a config data frame | complete / significant / filtered / annotation-only / region-methylation-matrix / configuration | Yes, 6 files |
| `methyl_results$dmr` | `observeEvent(live_result(), ...)` | `live_result()` | Summary only (`comparison, n_regions, n_sig`) | Not a file; a shared in-app object |

**Consistency check (plot vs. table vs. download):** Both tabs' plots are drawn from the *unfiltered* per-run table (so a user always sees the full effect-size/significance landscape, with the current threshold marked), while the *table* and most *downloads* reflect the currently filtered view — this is a deliberate, documented design (not a bug), but it means a region visible as a highlighted point on the volcano plot may not appear in the results table below it if it fails a display filter such as minimum CpGs or width. The SVA tab's boundary-condition inconsistency between plot-highlighting (`<=`) and the "significant" value-box/download (`<`) is noted as a confirmed minor issue in §7.

---

## 17. Reactive/Data Flow Architecture

**SVA tab:** `default_data()` (plain `reactive`, recomputed only if its own dependencies change — practically static) → `d_run()` (`eventReactive` on `input$d_run_btn`, `ignoreInit=TRUE` — computes nothing until the first click, and thereafter recomputes **only** on a click, not on any control change) → `d_filtered()` (plain `reactive` depending on `d_run()`, so it recomputes whenever `d_run()` changes, i.e., only after a click) → outputs. `d_has_run` (`reactiveVal`) is a one-way latch (`FALSE`→`TRUE` on first click, never reset) that gates whether `default_ui`/`d_valueboxes_ui`/`d_results_ui` show the "click Run" placeholder or the real content.

**DMR tab:** `sex_col()`, `id_cols()`, `anno_result()`, `dataset_norm_status()` are plain `reactive`s recomputed whenever `methyl_dataset$sample_sheet`/`$array_type`/`$beta` change. `live_result()` is an `eventReactive` on `input$live_run_btn` (`ignoreInit=TRUE`) wrapped in `withProgress()` — it re-executes the entire pipeline only on a click, never on a mere control change. `live_has_run` is a `reactiveVal` explicitly reset to `FALSE` whenever `methyl_dataset$beta` changes (`observeEvent(methyl_dataset$beta, live_has_run(FALSE), ignoreNULL=TRUE)`), forcing the user to re-run after loading a new dataset rather than silently continuing to show a stale result computed on the previous matrix. `live_filtered()`/`live_sig()` are plain `reactive`s that read `input$live_fdr/$live_dbeta/$live_direction/$live_mincpgs/$live_minwidth/$live_maxwidth` **directly**, not through `live_result()`, so they recompute immediately on every such control change without requiring a re-run — the opposite pattern from the SVA tab (see §7's audit note).

**Stale-result risk:** because `eventReactive` only ever recomputes on its triggering event, if a user changes the group column, sex, or covariates *without* clicking "Run" again, `live_result()` continues to reflect the previous configuration — this is standard, intended Shiny `eventReactive` behavior (and is what the "Run" button is for), not a defect; the UI's "Configuration & sample sizes" card does correctly display the configuration that was actually used for the currently shown `live_result()`, so a user can verify what is being shown even if they forgot to re-run after a change.

**All six `DT::renderDataTable` outputs** across the module (`d_table`, `d_region_table`, `live_table`, `live_region_table`, plus two more counted across both tabs) have `outputOptions(output, ..., suspendWhenHidden = FALSE)` explicitly set, with a code comment documenting this as a **confirmed, previously reproduced bug** in Shiny's own client-side visibility detection for DT tables nested two levels deep inside dynamically created `renderUI` blocks — without this override, the table's data would compute correctly server-side (verified via `testServer()`, per the comment) but never actually be sent to the browser. This is a **documented, applied fix**, not an open issue.

---

## 18. Tab-to-Tab Connections

The two DMR tabs are **fully independent** — neither reads the other's reactive state, and switching between them (a `tabsetPanel`) only suspends rendering of the hidden tab; per the file's own header comment, Shiny's tab suspension affects rendering only, not the underlying `eventReactive` results, so switching tabs and back never forces either `d_run()` or `live_result()` to recompute.

The only outward connection from the DMR module to the rest of the app is `methyl_results$dmr`, written exclusively by the DMR tab's `observeEvent(live_result(), ...)`. Running the SVA tab does not update this shared object at all — a limitation worth noting for anyone relying on `methyl_results$dmr` to reflect "the DMR module's most recent activity" generally, since it only ever reflects the live engine.

- **If a required input is missing** (no dataset, no sample sheet, no annotation for the array type): the DMR tab shows an explanatory `div` instead of the control panel; no error is thrown.
- **If a parameter is changed:** SVA tab — no visible effect until "Run" is clicked again (§7); DMR tab — FDR/Δβ/direction/CpG-count/width take effect immediately via `live_filtered()`/`live_sig()`; every other parameter (sex, group, covariates, seeding p, lambda, C, min.cpgs *as a DMRcate input*, probe-QC thresholds) requires clicking "Run DMR Analysis" again.
- **If the phenotype/group column is changed:** `live_level_ui` and `live_comparison_label_ui` update immediately to reflect the new column's levels (these are separate `renderUI` blocks reacting directly to `input$live_group_col`), but the actual analysis (`live_result()`) does not recompute until "Run" is clicked.
- **If the user re-runs:** the full pipeline re-executes from sample matching onward; `live_has_run` stays `TRUE`; any previous row selection in the results table is cleared by the table re-render (standard DT behavior on data replacement).
- **If the user downloads:** every download handler reads directly from the same reactive objects driving the on-screen table/plots at the moment of download (`live_result()`, `live_sig()`, `live_filtered()`), so a download always matches what is currently displayed, not a separately cached copy.

---

## 19. Packages and Dependencies

**Required for core DMR functionality (DMR tab):**
| Package | Functions used | Why required |
|---|---|---|
| `DMRcate` | `dmrcate()`, `extractRanges()`, `CpGannotated` S4 class | The region-calling algorithm itself — without it, no region can be called |
| `GenomicRanges` | `GRanges()`, `findOverlaps()` | Genomic-coordinate representation and region↔probe overlap joins |
| `IRanges` | `IRanges()` | Interval representation underlying `GRanges` |
| `S4Vectors` | `queryHits()`, `subjectHits()` | Extracting matched-pair indices from `findOverlaps()` results |
| `limma` | `lmFit()` (via `methyl_chunked_lmfit`), `makeContrasts()`, `contrasts.fit()`, `eBayes()`, `topTable()` | The per-CpG differential-methylation statistic feeding region calling; library-attached app-wide in `global.R` |
| `methods` (base) | `new("CpGannotated", ...)` | Constructs the S4 object DMRcate requires |

**Required for annotation (both tabs, indirectly):**
| Package | Used for | Notes |
|---|---|---|
| `IlluminaHumanMethylation450kanno.ilmn12.hg19` / `IlluminaHumanMethylationEPICanno.ilm10b4.hg19` | Chromosome/position/SNP manifest, via `methyl_get_annotation()` | Only these two array types are supported; loaded via `requireNamespace()` + `utils::data()`, never `library()`-attached, deliberately avoiding a documented `strsplit()`-masking risk from `Biostrings` |
| `ChAMPdata` | `probe.features` object, via `methyl_champ_probe_positions()` | SVA tab only — recovers chr/pos for the DMP tab's precomputed CpGs so the region-inspection panel can match them to a selected region |

**Used only for data manipulation / table display:**
| Package | Functions | Notes |
|---|---|---|
| `data.table` | `fread()` (in `global.R`'s loader functions) | Fast CSV reading for the bundled precomputed tables |
| `DT` | `datatable()`, `renderDataTable()`, `formatSignif()`, `dataTableOutput()` | All results tables |

**Used only for visualization:**
| Package | Functions | Notes |
|---|---|---|
| `ggplot2` | `ggplot()`, `geom_point()`, `geom_col()`, `geom_tile()`, `geom_boxplot()`, `geom_hline()`, `geom_abline()`, `scale_*`, `facet_grid()`, `theme()` | Every plot in the module; library-attached app-wide |
| `shinycssloaders` | `withSpinner()` | Loading-spinner wrapper around slow-rendering UI outputs |

**Base R / core statistics (no package attribution needed beyond `stats`/`utils`):** `stats::model.matrix`, `stats::complete.cases`, `stats::as.formula`, `stats::p.adjust`, `stats::aggregate`, `stats::na.omit`, `utils::head`, `utils::data`, `utils::write.csv`.

**Not used by this module** (despite appearing as unrelated search-alias text in the app's glossary feature): `bumphunter`, `comb-p`, `methylKit`, `minfi::getAnnotation()` (explicitly avoided by design, per `annotation.R`'s own comment).

---

## 20. Code and Scientific Audit

### Correctness
- **No issue identified:** sample alignment between the matrix and sheet is done by explicit ID matching (`methyl_sheet_sample_ids` + `intersect`), not positional assumption, and is re-validated (`>=6` common samples) before proceeding.
- **No issue identified:** phenotype/group labels are matched to samples via the same aligned `ph0`/`ph1` data frame the beta matrix subsetting uses, so group assignment cannot silently desynchronize from the matrix columns.
- **Confirmed implementation issue:** `n_cpgs_before_filter` mislabeling in the DMR tab (§8), causing the "CpGs tested after QC filters: X of Y" summary to always read "X of X."

### Statistical correctness
- **No issue identified:** the statistical test (moderated t-test via `limma`) is an appropriate, standard choice for a two-group methylation comparison, and is the same class of test the bundled reference pipeline used upstream of its own DMR step.
- **No issue identified:** multiple-testing correction is applied at the correct unit (per-region BH on the Stouffer statistic), separately from and in addition to the CpG-level correction already embedded in the input statistics — this correctly avoids double-dipping the same correction twice on the same axis.
- **Potential limitation:** the live engine explicitly does not apply SVA/bacon correction, and the module itself surfaces this as an in-app warning above λ>1.1 rather than hiding it — a transparent, code-confirmed caveat rather than a silent gap.

### Data integrity
- **No issue identified:** NA handling is explicit throughout (`is.na()` checks before FDR/Δβ comparisons, `na.rm=TRUE` in every group-mean/variance computation, `stats::complete.cases()` before the covariate-adjusted fit).
- **No issue identified:** duplicate CpG IDs are not explicitly deduplicated, but this was not observed to cause a specific failure mode in the traced code paths — flagged as a **potential limitation** only, since neither the precomputed CSVs nor typical Illumina manifests are expected to contain duplicate probe IDs.
- **No issue identified:** probes lacking genomic coordinates are explicitly excluded (`has_pos` check) before the fit, so DMRcate never receives an unplaceable probe.
- **No issue identified:** feature (probe) IDs are preserved as row names throughout the filtering/scale-conversion chain (`rownames(beta_scale_full)`, `rownames(tt)`), and match steps use these names rather than positional indices.

### Parameter integrity
- Covered in detail in §12 — every UI-exposed parameter either correctly reaches its intended analysis function, or is correctly restricted to display/filtering only, with one confirmed exception noted for the width/CpG-count controls' differing scope between tabs (documented, not miscoded).

### Reactive logic
- Covered in §17. **No issue identified** with dependency correctness (each analysis correctly re-triggers only on its intended events); the one asymmetry noted (SVA tab's non-live filter controls) is a **confirmed UI/implementation inconsistency**, not a crash risk.

### Output correctness
- **No issue identified:** downloads reflect exactly the reactive object driving the corresponding on-screen output at download time (verified by reading each `downloadHandler`'s `content` function against the reactive it calls).
- **Confirmed minor inconsistency:** SVA tab's plot-vs-valuebox/download significance-boundary mismatch (`<=` vs `<`), noted in §7.

### Error handling
- **No issue identified:** every major failure mode in the DMR tab's pipeline (missing data, insufficient samples, rank-deficient design, `limma`/`DMRcate` internal errors, zero seeded CpGs, zero returned regions) is caught via `validate(need(...))` or `tryCatch(..., error=function(e) validate(need(FALSE, ...)))`, producing a user-facing message rather than an uncaught exception or a silently empty table.

### Reproducibility
- **No issue identified:** no random-number generation was found anywhere in the traced DMR code path (`limma`'s moderated t-test and DMRcate's kernel smoothing are both deterministic given their inputs); no seed is set because none was found to be needed. The full parameter set (including `run_at`, a timestamp) is recorded in the downloadable configuration CSV, supporting reproducibility of a given run's settings.

---

## 21. Expected Methodology vs Actual Implementation

| Component | Expected/Scientific Concept | Actual Code Implementation | Assessment |
|---|---|---|---|
| Region calling | Any of several published DMR algorithms (DMRcate, bumphunter, comb-p, methylKit tiling) could plausibly be used | `DMRcate::dmrcate()` exclusively, both tabs | Matches the app's own stated scope and the bundled reference pipeline; no other algorithm is silently substituted or claimed |
| DMR-from-DMP grouping | A naive implementation might simply cluster already-significant single CpGs by proximity | Not implemented this way — DMRcate performs its own kernel-based candidate-region search over *all* tested CpGs' raw statistics, seeded by a nominal (not genome-wide) p-value | Scientifically more sound than naive clustering; matches §3's DMP-vs-DMR distinction |
| Multiple-testing correction | Region-level BH (or similar) FDR is expected for a defensible region-level significance claim | Implemented exactly this way, on the Stouffer statistic, per-run | Matches expectation |
| Confounder adjustment (SVA/surrogate variables) | A rigorous EWAS/region-level pipeline typically adjusts for unmeasured technical/batch confounding via SVA or similar, as the bundled reference pipeline itself does | The **SVA tab** does use SVA-corrected input (inherited from the DMP tab's precomputed statistics); the **DMR tab (live engine)** explicitly does **not** apply SVA/bacon correction and says so in-app | A confirmed, code-acknowledged limitation of the live engine specifically, not of the module as a whole |
| Cell-type adjustment | Whole-blood methylation studies often adjust for estimated cell-type proportions | Not performed anywhere in this module (a separate Cell-Type Deconvolution submodule exists but is not wired into the DMR fit) | Confirmed absence; not claimed as implemented anywhere in the DMR module's own UI text |
| Region annotation | Nearest-gene or regulatory-feature annotation is common in DMR pipelines | Direct-overlap gene annotation only, from `DMRcate::extractRanges()` | Narrower than "nearest gene," but this is what is actually computed and is accurately labeled `overlapping.genes`, not `nearest.gene` |
| Effect size | Region mean/max Δβ is standard | `meandiff`/`maxdiff` from DMRcate, plus this module's own added `ref_mean_beta`/`comp_mean_beta` from the actual filtered matrix | Matches expectation; the added per-group means are a genuine value-add over DMRcate's own output, independently computed and internally consistent with it |

---

## 22. Plain-English End-to-End Walkthrough

**SVA tab.** The module first reads two spreadsheets that were already computed outside the app — one row per genomic region, one file for female samples and one for male samples — from the app's own bundled data folder. It does not receive any per-sample methylation values here at all. When the user picks a sex and clicks "Run DMR Analysis," the module takes that one spreadsheet, tags each row as hyper- or hypomethylated based on the sign of its effect size, and remembers the current filter settings. It then produces a scatterplot and a genome-position plot straight from that full spreadsheet, and a separate, shorter table that only keeps the rows passing the user's chosen FDR/effect-size/CpG-count/width cutoffs. If the user clicks on a row of that table, the module goes back to a second bundled spreadsheet — this one holding individual CpG results, not regions — looks up each CpG's own chromosome and position from a third data source, and shows only the CpGs that physically fall inside the selected region's start/end coordinates, so the user can see which individual measurements made up that region-level result.

**DMR tab.** Here the module works from a real methylation matrix — one row per CpG probe, one column per sample — plus a spreadsheet describing each sample (which group it belongs to, its sex, any other covariates). The module first figures out which samples in the matrix have a matching row in the sample sheet, restricts to those, then further restricts to just the two groups (and, if chosen, just the one sex) the user wants to compare. It then removes probes that are missing too often, too invariant, sitting on a known SNP (if asked), or lacking a known chromosome position — because a region-calling algorithm cannot use a CpG it cannot place on the genome. On what remains, it fits one statistical model per CpG comparing the two groups (adjusting for any covariates the user picked), producing a test statistic and a raw p-value for every surviving CpG. It packages those per-CpG numbers, together with each CpG's chromosome and position, into the exact object format the DMRcate algorithm expects, and hands that to DMRcate, which slides a smoothing window along each chromosome and proposes candidate multi-CpG regions wherever enough nearby CpGs look different enough between the two groups. The module then applies its own separate correction across all of those candidate regions (treating each region as one statistical test) to decide which regions survive multiple-testing correction, attaches each region's own overlapping gene (from DMRcate) and its own directly recomputed group-wise mean methylation (from the actual matrix, not implied from DMRcate's summary statistic), and renders that as a table, several plots, and a heatmap. Clicking a region in the table shows the actual per-sample methylation values for that region's constituent CpGs, so the user can visually confirm the region reflects several CpGs moving together rather than one outlier probe.

---

## 23. Limitations and Important Caveats

1. **The live engine's statistics are not SVA/bacon-corrected.** The module surfaces this explicitly via a genomic-inflation-factor diagnostic and an in-app warning above λ>1.1, but does not itself apply the correction — a user comparing the live DMR tab's results to the SVA tab's numbers should expect the live tab's per-CpG (and therefore region-seeding) statistics to be less calibrated for datasets with unmeasured technical confounding, exactly as this dataset's own bundled methods document demonstrates was the case before SVA correction was applied to build the SVA tab's tables.
2. **No cell-type adjustment.** Neither tab adjusts for estimated cell-type composition; a DMR detected in a mixed-cell-type sample (e.g., whole blood) may partly or wholly reflect a shift in cell-type proportions between groups rather than a true within-cell-type methylation change.
3. **Annotation coverage is restricted to 450K/EPIC(v1).** The live engine cannot run region calling at all for EPICv2, WGBS, RRBS, or a generic "Custom array" dataset in this deployment, because no chromosome/position manifest is bundled for those array types.
4. **Region annotation is direct-overlap only,** not nearest-gene or regulatory-feature-aware; a biologically relevant but non-overlapping regulatory region near a DMR will not be reported.
5. **The SVA tab's controls do not update its display live** after the first "Run" click (§7); a user must click "Run DMR Analysis" again after changing any threshold to see the effect.
6. **The "CpGs tested after QC filters: X of Y" text in the DMR tab always reports X of X** due to the confirmed field-mislabeling bug in §8 — the true pre-filter probe count is not currently surfaced to the user anywhere in this module.
7. **No automatic sex covariate for combined-sex live comparisons** — left to the user's discretion (§8).
8. **DMRcate's candidate-region seeding uses a nominal, not genome-wide-corrected, p-value by design** (a documented DMRcate parameter, not a coding shortcut) — this is scientifically standard practice for this algorithm, but a reader unfamiliar with DMRcate's methodology could mistake "seed CpGs" for "significant CpGs," which they are not.

---

## 24. Thesis Implementation Paragraph

The Differentially Methylated Regions (DMR) module comprises two tabs within the application's Methylomics section. The first, "SVA," reproduces the sex-stratified region-level analysis of the reference thesis chapter by loading precomputed DMRcate (Peters et al., 2015, 2021) region tables — built with a Gaussian kernel bandwidth of 1000 nucleotides and a scaling factor of 2 on surrogate-variable-adjusted, bacon-corrected per-CpG statistics — and applying user-selected region-level FDR, effect-size, direction, CpG-count, and width filters without recomputing any statistic; selecting a called region additionally recovers its constituent CpGs' own differential-methylation statistics for inspection. The second tab, "DMR," implements a fully configurable, live region-calling pipeline: given a user-supplied or preloaded beta/M-value matrix and sample sheet, the module resolves sample identity, subsets by sex and by two user-chosen comparison groups, applies missingness/variance/SNP/genomic-position probe filters, fits a covariate-adjustable moderated t-test per CpG via `limma`, constructs a `DMRcate`-compatible annotated `GRanges` object from the resulting statistics, and calls candidate regions with `DMRcate::dmrcate()` under user-exposed tuning parameters (seeding p-value, kernel bandwidth, scaling factor, minimum CpGs per region). A region-level Benjamini–Hochberg correction is then applied to DMRcate's Stouffer combined-probability statistic across all candidate regions, and each surviving region is annotated with any directly overlapping gene and its own recomputed per-group mean methylation. Both tabs expose the resulting region table through volcano, genomic (Manhattan), and (live tab only) top-ranked effect-size and region×sample heatmap visualisations, alongside a region-level inspection view that recovers the constituent CpGs underlying any selected region, and both provide complete, significance-filtered, and analysis-configuration data exports. This design distinguishes DMR analysis from single-CpG (DMP) testing by aggregating spatially correlated methylation evidence into one region-level statistic, increasing statistical power to detect coordinated, biologically plausible multi-CpG effects that individual-CpG genome-wide correction can otherwise obscure, while making explicit, rather than concealing, the live engine's lack of surrogate-variable adjustment relative to the precomputed reference analysis it is built alongside.

## 25. Very Short Thesis Version

The DMR module takes either precomputed, surrogate-variable-adjusted per-CpG statistics (bundled sex-stratified thesis results) or a live, user-configured beta/M-value matrix and sample sheet, and calls genomic regions of coordinated differential methylation using `DMRcate::dmrcate()`, a Gaussian-kernel-smoothing region-calling algorithm, followed by an explicit Benjamini–Hochberg correction of the region-level Stouffer statistic. The output is a filterable, exportable table of candidate regions with region-level effect size, FDR, and overlapping-gene annotation, supported by volcano, Manhattan, and region-inspection visualisations. This region-level approach is used because true methylation differences are often distributed across several correlated, neighbouring CpGs rather than concentrated at one probe, a pattern single-CpG testing alone can fail to detect after genome-wide multiple-testing correction.

---

## 26. Final Audit Summary

| Area | Finding | Evidence from Code | Severity/Assessment |
|---|---|---|---|
| Input handling | Sample/sheet alignment by explicit ID matching, not position | `methyl_sheet_sample_ids()`, [qc.R:456](../../ArthOMix/R/methylomics/qc.R#L456) | Correct / appropriate |
| Input handling | SVA tab assumes bundled CSV column presence without an explicit schema check | [mod_methyl_dmr.R:218-224](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L218-L224) | Potential limitation (low risk — files are bundled, not user-supplied) |
| Preprocessing | Probe QC filters (missingness/variance/SNP/position) correctly gate region calling | [mod_methyl_dmr.R:693-720](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L693-L720) | Correct / appropriate |
| DMR calculation | `DMRcate::dmrcate()`/`extractRanges()` used exactly as documented by the bundled reference methodology | [mod_methyl_dmr.R:788-798](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L788-L798) | Correct / appropriate |
| DMR calculation | `n_cpgs_before_filter` mislabeled — equals post-filter count | [mod_methyl_dmr.R:718, 832](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L832) | Confirmed implementation issue (display text only, not the analysis itself) |
| Statistics | Moderated t-test (`limma`) per CpG, chunked for memory safety, documented as bit-for-bit verified against a whole-matrix fit | [mod_methyl_dmp.R:143](../../ArthOMix/R/methylomics/mod_methyl_dmp.R#L143) | Correct / appropriate |
| Statistics | Live engine has no SVA/bacon correction; disclosed via an in-app λ warning | [mod_methyl_dmr.R:903-905](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L903-L905) | Scientific limitation (disclosed, not hidden) |
| Multiple testing | Region-level BH on the Stouffer statistic, separate from CpG-level `ind.fdr` | [mod_methyl_dmr.R:767, 803](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L767) | Correct / appropriate |
| Multiple testing | SVA tab: plot-highlighting (`<=`) vs. value-box/download "significant" (`<`) boundary mismatch | [mod_methyl_dmr.R:98, 306, 377](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L306) | Confirmed minor issue (edge-case only) |
| Annotation | Direct-overlap gene annotation only, from DMRcate itself | [mod_methyl_dmr.R:800-802](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L800-L802) | Correct / appropriate (accurately labeled, not oversold) |
| Annotation | Region calling unsupported for EPICv2/WGBS/RRBS/Custom array | [annotation.R:18-21](../../ArthOMix/R/methylomics/annotation.R#L18-L21) | Scientific limitation (disclosed via in-app message) |
| Visualization | Plots use unfiltered data with threshold markers; tables/most downloads use filtered data — a documented, deliberate split | [mod_methyl_dmr.R:298-353](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L298-L353) | No issue identified (by design, and internally consistent) |
| Reactive flow | SVA tab: filter controls require re-clicking "Run" to take effect | [mod_methyl_dmr.R:268-289](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L268-L289) | Confirmed UI/implementation inconsistency (vs. the DMR tab's live filtering) |
| Reactive flow | `live_has_run` correctly resets on a new dataset load | [mod_methyl_dmr.R:616](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L616) | Correct / appropriate |
| Reactive flow | All DT tables forced to render regardless of Shiny's own (buggy, in this nested-`renderUI` layout) visibility detection | [mod_methyl_dmr.R:367, 456, 995, 1117](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L367) | Correct / appropriate (documented, applied fix) |
| Downloads | Every download reflects exactly the reactive object driving the matching on-screen output | e.g. [mod_methyl_dmr.R:997-1035](../../ArthOMix/R/methylomics/mod_methyl_dmr.R#L997-L1035) | Correct / appropriate |

### Overall Assessment

1. **What the DMR module currently does:** it provides two independent routes to a region-level differential-methylation result — a read-only exploration of a precomputed, thesis-grade, sex-stratified DMRcate analysis, and a fully live, user-configurable DMRcate pipeline operating on limma per-CpG statistics computed inside the app. Both routes call the same underlying algorithm (`DMRcate::dmrcate()`), and neither route implements a different or additional region-calling method than DMRcate.
2. **How the tabs work together:** they do not share analysis state; they are two views into related but separately triggered pipelines, connected only by shared helper functions and a shared visual style. Only the live tab publishes a summary object (`methyl_results$dmr`) for the rest of the app to read.
3. **What the major functions accomplish:** sample/probe QC and alignment (`methyl_filter_*`, `methyl_sheet_sample_ids`, `methyl_get_annotation`) prepare a genomically placeable, quality-controlled input; `limma` produces calibrated per-CpG statistics; `DMRcate` aggregates those statistics into genomic regions; this module's own code then applies a region-level FDR correction, attaches descriptive per-group means, and renders/exports the result.
4. **Whether the implementation matches its stated purpose:** yes — the module's own UI text ("Finds region-level methylation differences using DMRcate. Uses the bundled whole-blood analysis by default, or a live, configurable run on your own dataset") is an accurate description of what the code actually does; no functionality is claimed in the UI that the traced code does not implement.
5. **The most important confirmed limitations/issues:** (a) the mislabeled `n_cpgs_before_filter` field always reports 100% of probes retained regardless of actual filtering; (b) the SVA tab's filter controls are not live-reactive, requiring a second click of "Run" to take effect, unlike the DMR tab; (c) a minor `<=`/`<` boundary inconsistency in the SVA tab's "significant" definitions across plots vs. value-boxes/downloads; and (d), as a disclosed scientific limitation rather than a coding defect, the live engine's per-CpG statistics are not surrogate-variable/bacon-corrected, a gap the module itself surfaces via an inflation-factor warning rather than concealing.
