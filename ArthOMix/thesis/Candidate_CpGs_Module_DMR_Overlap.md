# Candidate CpGs (Module-DMR Overlap) — `mod_methyl_candidates.R`

**Source file:** `ArthOMix/R/methylomics/mod_methyl_candidates.R` (969 lines)
**Registration:** `mod_methyl_candidates_config` (id = `"candidates"`, group = `"Network"`, title = `"Candidate CpGs (Module-DMR Overlap)"`, icon = `"star"`), wired into the Methylomics vertical via `MX_MODULES` in `R/submodules_registry.R` and invoked as `m$server(paste0("mx_", m$config$id), methyl_dataset, methyl_results)` in `server.R:95`.
Prepared: 2026-08-26

This document is derived **exclusively** from the code in `mod_methyl_candidates.R`, cross-referenced only against the small number of external symbols it directly calls (`global.R`'s `load_default_wgcna_module_assignment()`, `load_default_dmr()`, `load_default_dmp()`, `METH_DATA_AVAILABLE`; `submodules_registry.R`'s `build_mx_context()`/`.format_results_block()`, which read this module's shared-results output). Nothing else in the application was inspected or modified. Anything not present in this code is not stated. Where the module's own in-code comments describe intent or scope decisions, they are quoted or paraphrased and attributed as such, not presented as independent claims. No code changes were made as part of this audit.

---

## 1. Module Overview

`mod_methyl_candidates.R`'s own header comment (lines 1–27) states the scope directly:

> "Methylation/WGCNA module CpGs -> DMR coordinates -> genomic overlap -> candidate CpGs -> filtering/prioritization -> results."

and:

> "This module owns its own data intake end to end (no reliance on the live WGCNA/DMR sub-modules' in-session reactive state), so this file is the only one touched to implement it."

Concretely: the module takes (a) a CpG-to-WGCNA-module assignment table and (b) a DMR (differentially methylated region) results table, resolves genomic coordinates for the module CpGs, computes which CpGs physically overlap which DMRs, tests whether each WGCNA module is statistically enriched for DMR-overlapping CpGs, and lets the user filter/rank the resulting CpG-DMR overlap pairs into a prioritized "candidate CpG" list. It is entirely a *downstream integration and filtering* module — it does not perform WGCNA clustering itself and does not call DMR region-calling; both are read as already-computed tables.

## 2. Scientific Purpose

**What a candidate CpG is (general concept).** In methylation biomarker discovery, a "candidate CpG" is a specific cytosine-guanine dinucleotide probe/site proposed as a plausible biological marker or mechanistic contributor, because independent lines of evidence converge on it — e.g., it sits inside a region of significant differential methylation *and* it belongs to a co-methylated network module implicated in the trait of interest.

**What a DMR is (general concept).** A differentially methylated region (DMR) is a genomic interval (not a single CpG) over which methylation differs, on aggregate, between two conditions/groups, typically identified by combining neighboring differentially-methylated-position (DMP) statistics (e.g., via DMRcate, comb-p, or bumphunter) into a region-level test with its own effect size and significance.

**What "module-DMR overlap" means in this application.** WGCNA (weighted gene co-expression network analysis, applied here to methylation) groups CpGs into co-methylation "modules" based on correlated methylation patterns across samples, independent of any DMR analysis. A CpG's module membership reflects a *network-level, unsupervised* clustering signal, while a DMR reflects a *supervised, region-level* differential-methylation signal from a separate statistical pipeline. "Module-DMR overlap," as implemented here, means: for every CpG that is both (i) assigned to a WGCNA module and (ii) has a resolvable genomic coordinate, test whether that coordinate falls inside a DMR (optionally with a flanking margin). This is a coordinate-based genomic intersection, not a co-membership or a Venn-diagram of gene lists.

**Why this overlap is scientifically useful.** A CpG that is both inside a significant DMR (independent statistical evidence of differential methylation) and inside a co-methylation module (independent evidence of coordinated methylation with other functionally related CpGs) is doubly supported by two orthogonal analyses. Concordance between the two increases confidence that the site is not a statistical artifact of a single pipeline, and module enrichment (are DMR-overlapping CpGs concentrated in a few specific modules, more than expected by chance?) can implicate an entire co-methylation program rather than an isolated site.

**Biological question this module is intended to answer.** "Which specific CpGs, and which WGCNA modules, show convergent evidence from both differential-region calling and co-methylation network structure — and can those CpGs be ranked into a short, prioritized candidate list for downstream validation or biomarker-panel construction?"

**What this module actually implements, precisely.** It performs (1) a coordinate-overlap join between module-assigned CpGs and filtered DMRs using `GenomicRanges::findOverlaps`, and (2) a per-module one-sided Fisher's exact test for enrichment of DMR-overlapping CpGs in each module against the tested CpG universe, with Benjamini-Hochberg FDR correction across modules. It does **not** implement any hypergeometric/permutation-based region-set enrichment, does not correct for spatial autocorrelation between neighboring CpGs (see §16), and does not perform any gene-level or pathway-level enrichment — those are explicitly out of scope for this file.

## 3. Scope of the Implementation

Per the header comment, this module is deliberately self-contained:

> "Upload path accepts a module-assignment table and a DMR-results table (a CpG annotation/coordinate table is optional if chromosome/position already live in the module table); column names are auto-detected the same way `mod_preprocessing.R`'s `pp_guess_col()` already does elsewhere in this app — exact names are never required."

The preloaded path reads the same on-disk tables that the (separately built) WGCNA and DMR sub-modules read (via `load_default_wgcna_module_assignment()`, `load_default_dmr()`, `load_default_dmp()` in `global.R`), but does **not** consume any in-session reactive state produced by `mod_methyl_wgcna.R` or `mod_methyl_dmr.R` — if a user has run a live WGCNA or DMR analysis elsewhere in the session with different parameters, this module is unaware of it and always re-reads the static preloaded files. Nor does it read `methyl_dataset` (the module's own `dataset` function argument, `mod_methyl_candidates.R:388`, is accepted but **never referenced anywhere in the function body** — a declared-but-unused parameter, consistent with the header's stated design, not a bug).

Normalisation, quality control, DMP calling, DMR calling, and WGCNA clustering are all performed **upstream**, outside this file. This submodule consumes their outputs only; it performs no normalisation and no re-clustering of its own.

## 4. Number of Tabs/Sub-tabs

**Total number of tabs/sub-tabs: 5**, defined in `mod_methyl_candidates_ui()` (`mod_methyl_candidates.R:371–384`) as a `tabsetPanel(id = ns("cd_subtabs"))`:

1. **Data & Filters**
2. **DMR-CpG Overlap**
3. **Module-DMR Overlap**
4. **Candidate CpGs**
5. **Visualization**

## 5. Input Data

### 5.1 Module assignment table (`module_assign` / `d$module_assign`)

- **Origin:** Preloaded — `load_default_wgcna_module_assignment(sex, merged = TRUE)` reading `METH_WGCNA_DIR/module_assignment_{sex}_merged10.csv` (`global.R`). Upload — a user-supplied CSV/TSV.
- **Standardized by:** `mcd_standardize_module_assign()` (lines 61–75).
- **Rows represent:** individual CpG probes.
- **Required columns (auto-detected):** a CpG/probe ID column (`MCD_CPG_ID_PATTERNS`, e.g. `cpg_id`, `probe_id`, `id`) and a module/module-color column (`MCD_MODULE_PATTERNS`, e.g. `module_color`, `module`, `color`).
- **Optional column:** a module-membership (kME) column (`MCD_KME_PATTERNS`).
- **Standardized structure:** `data.frame(cpg = character, module = character[, kme = numeric])`.
- **Missing values:** rows with `NA`/blank `cpg` or `module` are dropped silently (`out[!is.na(out$cpg) & nzchar(out$cpg) & ...]`).
- **Duplicates:** de-duplicated by CpG ID, keeping the **first** occurrence (`out[!duplicated(out$cpg), ]`) — silently, with no count of how many rows were dropped surfaced to the user.
- **Normalized upstream?** Yes — module assignment is the output of a completed WGCNA run; this module does not recompute it.

### 5.2 DMR results table (`dmr` / `d$dmr`)

- **Origin:** Preloaded — `load_default_dmr(sex)` reading `METH_DMR_DIR/dmr_{sex}_full.csv`. Upload — a user-supplied CSV/TSV.
- **Standardized by:** `mcd_standardize_dmr()` (lines 77–117).
- **Rows represent:** genomic regions (DMRs), not individual CpGs.
- **Required columns (auto-detected):** chromosome (`MCD_CHR_PATTERNS`), start (`^start$`/`dmr_start`/`region_start`), end (`^end$`/`dmr_end`/`region_end`).
- **Optional columns:** DMR ID (`MCD_DMR_ID_PATTERNS`, else auto-generated as `DMR0001`, `DMR0002`, ...), FDR (`MCD_FDR_PATTERNS`), raw p-value (`MCD_PVAL_PATTERNS`, only kept if distinct from the FDR column), Delta-Beta/effect size (`MCD_DBETA_PATTERNS`), CpG count per region (`MCD_NCPGS_PATTERNS`), gene (`MCD_GENE_PATTERNS`), and an explicit direction column (`MCD_DIRECTION_PATTERNS`) — if no explicit direction column exists but a Delta-Beta column does, direction is derived as `hyper`/`hypo`/`NA` from its sign.
- **Standardized structure:** `data.frame(chr, start, end, dmr_id[, dmr_fdr, dmr_pvalue, delta_beta, n_cpgs, gene, direction])`.
- **Missing values:** rows with missing chromosome/start/end are dropped.
- **Duplicates:** de-duplicated by `dmr_id`, first occurrence kept, silently.
- **Identifier format:** free-text chromosome strings, normalized later by `mcd_norm_chr()` (e.g. `"1"` → `"chr1"`, already-`chr`-prefixed values case-normalized; anything else left untouched).

### 5.3 CpG annotation / coordinate table (`annotation` / `d$annotation`)

- **Origin (upload path):** an optional user-uploaded file, standardized by `mcd_standardize_annotation()` (lines 119–141); if omitted, the module attempts to re-use the module-assignment file itself as a coordinate source (`ma_pos <- mcd_standardize_annotation(ma_raw)`) in case it already carries chromosome/position columns.
- **Origin (preloaded path):** `mcd_champ_full_annotation()` (lines 164–184), which loads `ChAMPdata::probe.features` (a Bioconductor reference annotation for Illumina methylation arrays) and extracts `cpg`, `chr` (prefixed `"chr"`), `pos` (`MAPINFO`), `gene`, `feature`, `island` (`cgi`). Cached in a package-level environment (`.mcd_champ_anno_cache`) so it is loaded from disk only once per R session.
- **Required columns (upload path, auto-detected):** CpG ID, chromosome, position.
- **Optional columns:** gene, CpG-island context (`MCD_ISLAND_PATTERNS`), genomic-region/feature (`MCD_FEATURE_PATTERNS`).
- **Missing values / duplicates:** rows with missing chromosome/position dropped; de-duplicated by CpG ID.
- **If unavailable:** genomic overlap and annotation-based filtering cannot run (surfaced to the user as an explicit note and, at analysis time, a `validate()` error).

### 5.4 CpG-level statistics table (`cpg_stats` / `d$cpg_stats`, optional)

- **Origin (preloaded path):** `load_default_dmp("sva", sex)` — the SVA-adjusted DMP table.
- **Origin (upload path):** re-uses the raw module-assignment file itself (`mcd_standardize_cpg_stats(ma_raw)`), i.e. if the uploaded module-assignment CSV also carries a p-value/FDR/Delta-Beta column, that column is picked up automatically; there is no separate upload control for a dedicated CpG-stats file.
- **Standardized by:** `mcd_standardize_cpg_stats()` (lines 147–162); requires a CpG ID column plus **at least one** of p-value, FDR, or Delta-Beta — if none present, returns "not ok" and the field is simply omitted.
- **Structure:** `data.frame(cpg[, p_value, fdr, delta_beta])`.
- **Consumed for:** merging CpG-level significance/effect size onto the CpG-DMR overlap table (subtab 2), and for the direction-consistency check between CpG-level and DMR-level methylation direction.

## 6. Complete End-to-End Pipeline

The actual sequence implemented in code (not the generic template in the audit brief):

```text
Data source selection (radio: Preloaded / Upload)
    |
Load & standardize  (mcd_standardize_module_assign / _dmr / _annotation / _cpg_stats)
    |  -> module_assign, dmr, annotation, cpg_stats  (loaded() eventReactive)
    v
Filter UI rendered from detected columns  (filters_ui: module / DMR significance / overlap definition)
    |
[Tab 2] DMR filtering  (mcd_filter_dmrs: FDR / p-value / |Delta-Beta| / min-CpGs / direction / chromosome)
    |
    + module filtering  (grey exclusion / specific modules / min size / min |kME|)
    v
Coordinate-based CpG<->DMR overlap  (mcd_compute_overlap: merge() by cpg -> GRanges -> findOverlaps)
    |  -> joined (one row per CpG-DMR pair) + cpg_universe (all coordinate-resolved module CpGs)
    v
Merge CpG-level stats onto joined pairs (optional); derive direction + direction_consistency;
apply direction-consistency filter
    |  -> overlap_run() result: "DMR-CpG Overlap" table + counts
    v
[Tab 3] Same DMR filters, but EVERY module (not just selected) x-ref'd against the SAME overlap join
    |  -> per-module n_overlap_cpgs / n_dmrs / one-sided Fisher's exact test / BH-FDR
    v
[Tab 4] Candidate CpGs: further filter (module / FDR / |Delta-Beta| / direction / consistency /
region / island) + score (significance | effect | membership | combined) the Tab-2 overlap table
    |  -> ranked "prioritized candidate CpGs" table
    v
[Tab 5] Visualization: six independent plots, each gated on Tab 2 or Tab 3 having been run
    |
Downloads (CSV for tables, PNG for plots) at every stage
    |
results$candidate_cpgs <- summary counts  (written for the Assistant/ArthOChat context — see §17.4)
```

## 7. Tab-to-Tab Data Flow

```text
Tab 1 (Data & Filters)
  produces: loaded()  [module_assign, dmr, annotation, cpg_stats, detected, notes]
  produces: filter inputs (f_module*, f_dmr_*, f_overlap_mode, f_flank_bp, f_chr, f_consistency)
      |
      v
Tab 2 (DMR-CpG Overlap)  --uses loaded() + ALL Tab-1 filters, including module selection--
  produces: overlap_run()  [table = CpG x DMR pairs, counts]
      |                                   |
      v                                   v
Tab 4 (Candidate CpGs)          Tab 5 (Visualization: viz_dmrbar, viz_direction, viz_annot, viz_volcano)
  further filters/scores               (viz_volcano also OPTIONALLY highlights Tab 4's
  Tab 2's overlap_run() table           candidate_result()$table$cpg if Tab 4 has been run)
  produces: candidate_result()

Tab 3 (Module-DMR Overlap)  --uses loaded() + Tab-1's DMR filters + grey-exclusion ONLY
                               (NOT module-selection/min-size/min-kME filters; every module
                               is scored, by design -- see header comment lines 675-678)--
  produces: module_overlap_run()  [per-module enrichment table]
      |
      v
Tab 5 (Visualization: viz_modbar, viz_heatmap)
```

Dependency structure, stated precisely:

- **Tab 1 -> Tab 2 and Tab 1 -> Tab 3** are both direct: both re-read `loaded()` and re-apply the DMR filters from Tab 1, but Tab 2 additionally applies the module-selection filters while Tab 3 deliberately does not (it scores every module for comparability).
- **Tab 2 -> Tab 4**: Tab 4 does not recompute the CpG-DMR overlap; it only filters and scores the table Tab 2 already produced (`ov <- overlap_run(); req(ov)`, line 833). Tab 4 cannot be opened meaningfully until Tab 2 has been run (`mod_methyl_candidates.R:775`: "Run 'DMR-CpG Overlap' first").
- **Tab 3 is independent of Tab 2's run** in the sense that it does not read `overlap_run()`'s result — it calls `mcd_compute_overlap()` a second time itself, with its own (unfiltered-by-module) CpG set. This means the overlap join is computed **twice** (once in Tab 2, once in Tab 3) whenever both are run, with no caching or sharing between them.
- **Tab 5** is a pure presentation layer over Tabs 2/3/4's already-computed reactives; it triggers no new analysis logic beyond the plotting functions themselves, and each of its six plots is gated independently on whichever upstream tab it depends on.
- Tabs are **not** sequential in a hard-enforced sense — Shiny does not prevent a user from clicking directly into Tab 4 or Tab 5 before running Tab 2/3, and the module handles this by rendering an explanatory `empty-note` instead of an error (see §14).

## 8. Tab 1 — Data & Filters

### Purpose
Select and validate a data source (preloaded or uploaded), report what was loaded and how columns were auto-detected, and expose every downstream filter control in one place.

### Inputs
- `data_source` (radio: `preloaded` / `upload`), `pre_sex` (radio: `female`/`male`), `load_btn`.
- Upload-only: `up_module_file`, `up_dmr_file`, `up_annot_file` (all `fileInput`, CSV/TSV/TXT).
- Filters (rendered only after `has_loaded()` is `TRUE`): `f_module_mode`, `f_module`, `f_module_minsize`, `f_module_min_kme`, `f_exclude_grey`, `f_dmr_fdr`, `f_dmr_p`, `f_dmr_dbeta`, `f_dmr_mincpgs`, `f_direction`, `f_overlap_mode`, `f_flank_bp`, `f_chr`, `f_consistency`.
- `goto_overlap_btn` (navigates to Tab 2).

### Input Data Structure
See §5 for the four possible input tables. On the preloaded path, only the `sex` stratum selector varies the data source; on the upload path, up to three files are combined.

### Functions Used
`mcd_standardize_module_assign()`, `mcd_standardize_dmr()`, `mcd_standardize_annotation()`, `mcd_standardize_cpg_stats()`, `mcd_champ_full_annotation()`, `mcd_detected_list()`, plus base/Shiny: `data.table::fread()`, `validate()`/`need()`, `conditionalPanel()`, `radioButtons()`, `fileInput()`, `numericInput()`, `selectizeInput()`, `checkboxInput()`.

### Processing
1. `has_loaded` (a `reactiveVal`) is set `TRUE` on `load_btn`, reset `FALSE` whenever `data_source` changes.
2. `loaded()` (`eventReactive` on `load_btn`) branches on `data_source`:
   - **Upload:** parses both required files with `data.table::fread()`, standardizes each, and — if no separate annotation file was supplied — attempts to reuse the module-assignment file as a coordinate source before giving up on genomic overlap.
   - **Preloaded:** requires `METH_DATA_AVAILABLE`; calls the three `load_default_*()` functions for the chosen sex, standardizes each, and reads the full `ChAMPdata` annotation.
3. `output$load_summary_ui` reports CpG/module/DMR counts, any warning notes (e.g. "no CpG-level statistics available"), and the detected column mapping.
4. `output$filters_ui` inspects which optional columns were actually detected (`has_kme`, `has_dmr_fdr`, `has_dmr_dbeta`, `has_dmr_dir`, `consistency_available`, ...) and renders only the controls that data supports — a control that has no backing column is replaced with an `empty-note` explaining why it is unavailable, rather than being hidden silently.

### Biological Interpretation
This tab has no biological output of its own; it establishes the module/DMR/CpG universe and analyst-chosen thresholds that every later biological interpretation depends on.

### Output Data
`loaded()` (the standardized four-table bundle), consumed by every other tab; no table/plot/download is produced directly by this tab.

### Output Data Structure
See §5.

### Connection to Other Tabs
Feeds `loaded()` and every filter input to Tabs 2, 3, and (indirectly, via Tab 2's table) Tab 4.

### Audit Findings
- Column auto-detection is transparent (the "Detected column mapping" `<details>` panel), which is good practice for a module whose correctness depends entirely on correct column identification.
- Silent row-dropping during standardization (duplicate CpG/DMR IDs, rows with missing coordinates) is not quantified anywhere in the UI — the user sees the *final* row count only, with no "N rows dropped due to missing/duplicate identifiers" figure. See §14 and §19 (MEDIUM).
- `dataset` (the shared `methyl_dataset` argument) is accepted by `mod_methyl_candidates_server()` but never used — see §17.3.

---

## 9. Tab 2 — DMR-CpG Overlap

### Purpose
Compute, for the currently selected module(s) and currently filtered DMRs, exactly which CpGs fall inside (or within a flanking distance of) which DMRs.

### Inputs
`overlap_run_btn`, plus every Tab-1 filter (module selection/size/kME, grey-exclusion, DMR FDR/p/Delta-Beta/min-CpGs/direction, overlap mode, flank distance, chromosome restriction, direction-consistency).

### Input Data Structure
`loaded()`'s `module_assign`, `annotation`, `dmr` tables (§5.1–5.3), plus optionally `cpg_stats` (§5.4).

### Functions Used
`mcd_filter_dmrs()`, `mcd_compute_overlap()` (which itself calls `merge()`, `GenomicRanges::GRanges()`, `IRanges::IRanges()`, `GenomicRanges::findOverlaps()`, `S4Vectors::queryHits()`/`subjectHits()`), `mcd_norm_chr()`, base `merge()` (for CpG stats), `validate()`/`need()`, `DT::renderDataTable()`/`DT::datatable()`/`DT::formatSignif()`, `downloadHandler()`.

### Processing (execution order, `overlap_run` eventReactive, lines 577–622)
1. Require `d$annotation` to be non-NULL (hard stop otherwise).
2. Filter DMRs via `mcd_filter_dmrs()` using every DMR-side Tab-1 filter; stop if zero DMRs remain.
3. Filter the module-assignment table: exclude grey/gray if checked, restrict to selected module(s) if in "specific" mode, drop modules below the minimum-size threshold, drop CpGs below the minimum `|kME|` if that column exists; stop if zero CpGs remain.
4. Compute `flank` (0 unless "flank" overlap mode is chosen).
5. Call `mcd_compute_overlap(ma, annotation, dmr_f, flank)`; stop if zero overlapping pairs.
6. If `cpg_stats` exists, left-merge it onto the overlap table by `cpg`.
7. Derive `direction` (from CpG-level `delta_beta` sign, falling back to the DMR's own direction if no CpG-level Delta-Beta exists) and, when both a CpG-level and DMR-level direction are available, `direction_consistency` ("Consistent"/"Inconsistent"), then apply the direction-consistency filter if requested; stop if zero rows remain.
8. Return the table plus summary counts (`n_dmr_tested`, `n_dmr_passing`, `n_cpg_tested`, `n_overlap_cpgs`, `n_overlap_dmrs`).

### Biological Interpretation
Each output row is one physical co-location event: a specific CpG, assigned to a specific WGCNA module, whose genomic coordinate falls inside (or near) a specific significant DMR. This is the base evidence unit for everything downstream — "which CpGs are supported by both analyses, and where."

### Output Data
A `valueBox` row (DMRs in table / DMRs passing filters / module CpGs tested / overlapping CpGs) and a full CpG-DMR pair `DT::datatable` with CSV download (`candidate_cpgs_dmr_overlap.csv`).

### Output Data Structure
Columns present depend on what was detected upstream, but typically include: `cpg`, `module`, `kme`, `chr`, `pos`, `dmr_id`, `dmr_chr`, `dmr_start`, `dmr_end`, `dmr_fdr`, `dmr_pvalue`, `dmr_delta_beta`, `dmr_n_cpgs`, `dmr_gene`/`gene`, `dmr_direction`, `p_value`, `fdr`, `delta_beta`, `direction`, `direction_consistency`. Displayed via `mcd_pretty()`, which renames a fixed whitelist of columns to human-readable labels (`MCD_PRETTY_MAP`) and silently drops any column not in that map.

### Connection to Other Tabs
Directly consumed, unmodified, by Tab 4 (`candidate_result()` filters this exact table) and by four of Tab 5's six plots (`viz_dmrbar`, `viz_direction`, `viz_annot`, `viz_volcano`).

### Audit Findings
- **`n_cpg_tested` denominator mismatch:** the "Module CpGs tested" `valueBox` reports `nrow(ma)` — the module-filtered CpG count **before** the coordinate join in `mcd_compute_overlap()` — not the count of CpGs that actually had a resolvable genomic coordinate. Any CpG whose ID fails to match the annotation table (case, whitespace, naming convention) is silently excluded from the overlap computation but still counted in this stat. See §19 (MEDIUM).
- A CpG overlapping multiple DMRs correctly produces multiple output rows (one per pair), and `n_overlap_cpgs` is explicitly computed as `length(unique(out$cpg))` rather than `nrow(out)` — the row-count-vs-unique-CpG-count distinction is handled correctly.
- `merge(module_assign, annotation, by = "cpg")` is an inner join; any module CpG absent from the annotation table (or vice versa) is dropped without an explicit count being surfaced (see previous point).

---

## 10. Tab 3 — Module-DMR Overlap

### Purpose
For every WGCNA module (not just a user-selected subset), test whether that module's CpGs are statistically enriched for overlap with the filtered DMRs, relative to the full tested CpG universe.

### Inputs
`modoverlap_run_btn`, plus Tab 1's DMR filters (FDR/p/Delta-Beta/min-CpGs/direction/chromosome) and grey-exclusion — **not** the module-selection, minimum-size, or minimum-kME filters (see header comment, lines 675–678, and §7).

### Input Data Structure
Same `module_assign`/`annotation`/`dmr` tables as Tab 2, but with the module set restricted only by grey-exclusion.

### Functions Used
`mcd_filter_dmrs()`, `mcd_compute_overlap()`, base `table()`/`unique()`/`length()`, `stats::fisher.test()`, `stats::p.adjust(method = "BH")`, `do.call(rbind, ...)`, `DT::renderDataTable()`.

### Processing (`module_overlap_run`, lines 680–722)
1. Filter DMRs identically to Tab 2 (same `mcd_filter_dmrs()` call with the same inputs).
2. Build `ma_all` = the full module assignment minus grey/gray (no module-selection filter applied).
3. Compute `flank` and call `mcd_compute_overlap(ma_all, annotation, dmr_f, flank)`, producing `universe` (all coordinate-resolved CpGs, `cpg_universe`) and `joined` (overlap pairs).
4. `total_n` = unique CpGs in the universe; `total_overlap` = unique CpGs overlapping any filtered DMR.
5. For each module `m`: `a` = CpGs in module `m` that overlap; `b` = CpGs in module `m` that don't; `c_` = overlapping CpGs outside module `m`; `dd` = the remainder. A 2x2 contingency table `[[a,b],[c_,dd]]` is tested with `stats::fisher.test(..., alternative = "greater")` (one-sided, testing enrichment only, not depletion). The Fisher test is skipped (p-value `NA`) if any cell would be negative or if `total_n == 0`, rather than erroring.
6. `res$fdr <- stats::p.adjust(res$p_value, method = "BH")` across all tested modules; results sorted by raw p-value.

### Biological Interpretation
A module with a low enrichment p-value/FDR contains disproportionately more DMR-overlapping CpGs than expected if overlap were distributed uniformly across modules — i.e., that co-methylation program is not just correlated internally, but its member CpGs are independently flagged by the DMR-calling pipeline more often than chance. `odds_ratio` quantifies the strength of that association.

### Output Data
An explanatory sentence (`%d modules tested against a background of %d CpGs...`), a per-module statistics table, and a CSV download (`module_dmr_overlap_statistics.csv`).

### Output Data Structure
`module, n_module_cpgs, n_overlap_cpgs, n_dmrs, pct_overlap, odds_ratio, p_value, fdr` — displayed with manually assigned column headers (not routed through `mcd_pretty()`/`MCD_PRETTY_MAP`, unlike Tabs 2 and 4).

### Connection to Other Tabs
Feeds two of Tab 5's six plots (`viz_modbar`, `viz_heatmap`). Does **not** feed Tab 4 (Candidate CpGs) — Tab 4 only reads Tab 2's table.

### Audit Findings
- **Only enrichment, not depletion, is tested** (`alternative = "greater"` is hard-coded). A module with significantly *fewer* DMR-overlapping CpGs than expected cannot be flagged by this test as currently implemented. See §19 (LOW/MEDIUM — depends on whether depletion is scientifically relevant to this application's intended use).
- The one-sided Fisher's exact test treats each CpG as an independent Bernoulli trial. Neighboring CpGs within the same DMR or the same module are not independent draws (methylation is spatially and co-methylation-network correlated) — this is a general, known caveat of naive Fisher enrichment on genomic co-localization data, not something this specific implementation corrects for or claims to correct for. See §16 (MEDIUM, scientific).
- Defensive guarding of the 2x2 table (`a>=0 && b>=0 && c_>=0 && dd>=0 && total_n>0`) before calling `fisher.test()` is correctly implemented and avoids a crash on a degenerate/empty contingency table.
- The overlap join (`mcd_compute_overlap()`) is recomputed independently from Tab 2's call, with a different (unfiltered-by-module) input — this is by design, but means the same `GenomicRanges::findOverlaps()` work happens twice per full walkthrough of the module, with no caching.
- `modoverlap_table`'s column headers are set manually (`colnames(df) <- c(...)`) rather than via `mcd_pretty()`/`MCD_PRETTY_MAP`, the mechanism Tabs 2 and 4 use — a minor internal inconsistency, not a correctness bug (§19, LOW).

---

## 11. Tab 4 — Candidate CpGs

### Purpose
Take Tab 2's CpG-DMR overlap table and narrow/rank it into a final "prioritized candidate CpGs" list, using a user-chosen combination of filters and a scoring formula.

### Inputs
Gated entirely on `ov_has_run()` (Tab 2 must have been run first; the tab otherwise shows: *"Run 'DMR-CpG Overlap' first — Candidate CpGs filters/prioritizes that overlap table."*). Controls, each rendered only if the backing column exists in Tab 2's output: `c4_module`, `c4_dmr_fdr`, `c4_cpg_fdr`, `c4_min_dbeta`, `c4_min_dmr_dbeta`, `c4_direction`, `c4_consistency`, `c4_feature`, `c4_island`, `cand_rank_mode` (`significance` / `effect` / `membership` / `combined`), `cand_apply_btn`.

### Input Data Structure
`overlap_run()$table` — the exact table produced by Tab 2 (§9's output structure), not re-derived from raw data.

### Functions Used
Base subsetting (`df[...]`), `validate()`/`need()`, `-log10()`, `order()`, `round()`, `DT::renderDataTable()`, `downloadHandler()`. `output$cand_formula_ui` dynamically renders which scoring terms are actually available before the user clicks "Apply," using the same column-presence checks as the scoring code itself.

### Processing (`candidate_result`, lines 832–866)
1. Sequential filtering by module selection, DMR FDR, CpG FDR, minimum CpG-level `|Delta-Beta|`, minimum DMR-level `|Delta-Beta|`, direction, consistency, genomic region (feature), and CpG-island context — each step conditional on the corresponding column actually being present; stop if zero rows remain.
2. **Score computation**, additive across whichever factors are available for the chosen `cand_rank_mode`:
   - *Significance*: `-log10(dmr_fdr)` + `-log10(fdr)` (CpG-level), each floored at `1e-300` before the log.
   - *Effect*: `10 * |delta_beta|` (CpG-level) + `10 * |dmr_delta_beta|` (DMR-level).
   - *Membership*: `|kme|`.
   - *Combined*: the sum of all of the above, plus `+1` if `direction_consistency == "Consistent"`.
3. `candidate_score` is attached and the table sorted descending by score.

### Biological Interpretation
The ranked output represents the module's best attempt at prioritizing which overlap-supported CpGs are most likely to be biologically meaningful, by whichever single axis (statistical significance, effect magnitude, network centrality/kME, or a combined additive score) the analyst chooses. This is explicitly a **heuristic ranking**, not a re-derived statistical test — the "Candidate Score" has no associated p-value, null distribution, or formal multiple-testing interpretation.

### Output Data
A summary sentence (candidate count + factors used in the score), a ranked `DT::datatable`, and a CSV download (`prioritized_candidate_cpgs.csv`).

### Output Data Structure
Same columns as Tab 2's output, plus `candidate_score`.

### Connection to Other Tabs
Consumes Tab 2's output exclusively. Optionally feeds back into Tab 5's `viz_volcano` plot (candidate CpGs are highlighted in red if Tab 4 has been run; the plot still renders, unhighlighted, if not).

### Audit Findings
- **The additive combined score mixes unlike scales without normalization**: `-log10(p)`-type terms (unbounded, easily reaching 10–300 for tiny p-values) are summed directly with `10*|Delta-Beta|` (bounded roughly 0–10 for a 0–1 methylation-proportion effect size) and `|kME|` (bounded 0–1). In `combined` mode, an extreme significance value can dominate the score almost completely, effectively reducing the "combined" ranking to a significance-only ranking in practice for most real data. This is a genuine scientific/statistical design concern, not a bug — the module does not claim the combined score is calibrated or weighted, but a reader could reasonably expect "combined" to balance the factors more evenly. See §19 (MEDIUM).
- The `1e-300` floor used here for the significance score differs from the `.Machine$double.xmin` (~4.9e-324) floor used in `mcd_plot_enrichment_heatmap()` (§13) — a minor internal inconsistency with negligible practical effect (both are far below any real FDR value), noted for completeness (§19, LOW).
- No independent statistical test (e.g., a combined-evidence significance measure) is computed for the ranked list — correctly, the module does not claim one exists.

---

## 12. Tab 5 — Visualization

### Purpose
Render one plot per pre-defined chart type, each independently gated on whichever upstream tab's result it depends on, via a single shared registration helper.

### Inputs
Six independent `*_btn` action buttons (one per chart); no additional filter controls of its own (all filtering happens upstream in Tabs 1/2/3/4).

### Input Data Structure
Whichever upstream reactive each plot depends on (see table below) — no new data is loaded or transformed at this tab beyond what the plotting functions do internally.

### Functions Used
`mcd_viz_block()` (UI card+button+placeholder), `register_viz()` (shared button/output/download wiring, lines 919–942), and the six plotting functions: `mcd_plot_module_bar()`, `mcd_plot_enrichment_heatmap()`, `mcd_plot_dmr_bar()`, `mcd_plot_direction()`, `mcd_plot_annotation_dist()`, `mcd_plot_candidate_volcano()` — all built on `ggplot2` (`ggplot()`, `geom_col()`, `geom_tile()`, `geom_point()`, `facet_wrap()`, `scale_fill_manual()`/`scale_fill_gradient()`/`scale_color_manual()`, `theme_arthomix()` — the app's shared ggplot theme, not defined in this file), plus `ggsave()` for PNG export.

### Processing / per-plot dependency and description text

| Plot (`id_prefix`) | Depends on | Description shown in UI | Plotting function |
|---|---|---|---|
| `viz_modbar` — "Overlapping CpGs by module" | Tab 3 (`modov_has_run()`) | "Requires 'Module-DMR Overlap' to have been run." | `mcd_plot_module_bar()` |
| `viz_heatmap` — "Module enrichment heatmap" | Tab 3 (`modov_has_run()`) | "-log10(p) per module, from the same Module-DMR Overlap run." | `mcd_plot_enrichment_heatmap()` |
| `viz_dmrbar` — "Candidate CpGs per DMR" | Tab 2 (`ov_has_run()`) | "Top 20 DMRs by number of overlapping candidate CpGs. Requires 'DMR-CpG Overlap' to have been run." | `mcd_plot_dmr_bar()` (top-20 hard-coded via `n = 20` default argument) |
| `viz_direction` — "Effect direction" | Tab 2 (`ov_has_run()`) | "Hyper- vs hypomethylated candidate CpGs." | `mcd_plot_direction()` |
| `viz_annot` — "Genomic annotation distribution" | Tab 2 (`ov_has_run()`) | "CpG island context and/or genomic region, when available in the loaded data." | `mcd_plot_annotation_dist()` |
| `viz_volcano` — "Candidate CpG plot" | Tab 2 (`ov_has_run()`); optionally Tab 4 (`cand_has_run()`) for highlighting | "Delta-Beta vs -log10(significance); prioritized candidates (Candidate CpGs tab) are highlighted when available." | `mcd_plot_candidate_volcano()` |

Each plot is wrapped in `validate(need(isTRUE(ok$ok), ok$reason))` inside `register_viz()`'s `plot_obj` reactive, so an unmet prerequisite produces an inline message rather than an error, and each has its own PNG download (`ggsave(..., width = 8, height = 5.5, dpi = 300, bg = "white")`).

### Biological Interpretation
- `viz_modbar` / `viz_heatmap`: which modules carry the most DMR-overlap signal, and how strongly (visual companions to Tab 3's table).
- `viz_dmrbar`: which specific DMRs are "hubs" accumulating the most candidate CpGs.
- `viz_direction`: whether candidate CpGs skew hyper- or hypomethylated overall.
- `viz_annot`: where candidate CpGs sit relative to CpG islands and genomic features (promoter/body/etc.), a classic methylation-biology context check.
- `viz_volcano`: the standard effect-size-vs-significance view, with prioritized candidates highlighted for at-a-glance triage.

### Output Data
Six `plotOutput`s and six PNG downloads; no new tables.

### Output Data Structure
N/A (plots only).

### Connection to Other Tabs
Purely a downstream consumer of Tabs 2/3/4; produces nothing consumed elsewhere.

### Audit Findings
- `mcd_plot_direction()`, `mcd_plot_annotation_dist()`, and `mcd_plot_candidate_volcano()` each independently `validate(need(...))` for the specific columns they require, correctly degrading to an informative message rather than crashing when, e.g., no direction information is available.
- The top-20-DMRs cap in `mcd_plot_dmr_bar()` is a function-default argument (`n = 20`), not exposed as a UI control — a hard-coded, non-configurable value (§13).
- `mcd_module_colors()` (used only by `mcd_plot_module_bar()`) recycles a fixed 7-color palette if there are more than 7 non-WGCNA-colour-named modules — see §13.

---

## 13. Hard-Coded Parameters

| Value | Location | Configurable? | Assessment |
|---|---|---|---|
| Default DMR FDR threshold `0.05` | `f_dmr_fdr` numericInput default | Yes (user can change) | Standard, reasonable default. |
| Default DMR raw-p threshold `1` (effectively off) | `f_dmr_p` numericInput default | Yes | Appropriate — off by default. |
| Default min `|Delta-Beta|`, min CpGs/DMR, min module size, min `|kME|` = `0` (off) | Several `numericInput`s | Yes | Appropriate — off by default. |
| Fisher's exact test `alternative = "greater"` | `mod_methyl_candidates.R:710` | **No** — not exposed in UI or as a function argument | One-sided only; enrichment-only, cannot detect depletion (§10, §19). |
| BH FDR correction method | `mod_methyl_candidates.R:719` | No | Standard, defensible default; not user-selectable but appropriate. |
| Significance-score p-value floor `1e-300` | `mod_methyl_candidates.R:851–852` | No | Inconsistent with the `.Machine$double.xmin` floor used in `mcd_plot_enrichment_heatmap()` (line 290) — cosmetic inconsistency only. |
| Top-20 DMR cap in `mcd_plot_dmr_bar()` | `mod_methyl_candidates.R:299` (`n = 20` default arg) | No (code-level only) | Reasonable for readability; not disclosed as a fixed cap anywhere in the UI text beyond the plot's own description sentence. |
| PNG export `width = 8, height = 5.5, dpi = 300` | `register_viz()`, line 940 | No | Minor; standard export defaults. |
| Categorical fallback palette (7 colors) in `mcd_module_colors()` | `mod_methyl_candidates.R:272–274` | No | Recycles (repeats colors) once module count exceeds 7 for non-standard (non-WGCNA-color) module names — could visually conflate distinct modules in the bar plot. |
| Chromosome-name normalization pattern (`"chr"` prefixing) | `mcd_norm_chr()`, lines 189–193 | No | No genome-build (hg19/hg38) awareness or validation — see §16. |
| Column-detection regex patterns (`MCD_*_PATTERNS`) | Lines 31–44 | No (code-level only) | Extensive and reasonably permissive; not user-configurable, but this is consistent with the app-wide `pp_guess_col()` convention the header comment cites. |
| Default `flank` distance `0` bp | `f_flank_bp` numericInput default | Yes | With flank `0`, "flank" mode and "inside" mode are functionally identical unless the user explicitly enters a positive value (§16, LOW). |

## 14. Edge-Case Audit

| Scenario | Actual code behaviour |
|---|---|
| No CpGs overlap the filtered DMRs | `overlap_run()` / `module_overlap_run()` each `validate(need(nrow(out) > 0, ...))` with an actionable message ("increasing the DMR FDR threshold, adding flanking distance, or changing the selected module"). |
| One CpG overlaps multiple DMRs | Produces multiple rows in the joined table (one per pair); `n_overlap_cpgs` correctly de-duplicates by unique CpG for the summary stat. |
| Multiple CpGs occur in one DMR | Handled identically — multiple rows, correctly aggregated in `mcd_plot_dmr_bar()` via `stats::aggregate(cpg ~ dmr_id, ..., function(x) length(unique(x)))`. |
| A module has no DMR-overlapping CpGs | In Tab 3, `a = 0` for that module; Fisher's test still runs (guarded, does not crash) and yields a non-significant/near-1 p-value; the module still appears in the results table with `n_overlap_cpgs = 0`. |
| DMR table is empty after filtering | `validate(need(nrow(dmr_f) > 0, "No significant DMRs remain after the current DMR filters..."))` in both Tab 2 and Tab 3. |
| Module-assignment table is empty after filtering | `validate(need(nrow(ma) > 0, "No CpGs remain in the selected module(s) after the module filters."))` (Tab 2); analogous grey-only check in Tab 3. |
| Missing CpG IDs in input | Dropped during standardization (`mcd_standardize_module_assign`/`_annotation`/`_cpg_stats`), silently, before the user ever sees the data. |
| Duplicate CpG IDs in input | Dropped (first occurrence kept) during standardization, silently, with no count surfaced. |
| Missing genomic coordinates | Rows dropped both at standardization (`mcd_standardize_annotation`/`_dmr`) and again inside `mcd_compute_overlap()` (`cpg_pos[!is.na(chr) & !is.na(pos)]`); if the *entire* CpG universe ends up with unresolved coordinates, `validate(need(nrow(universe) > 0, "None of the module-assigned CpGs have resolvable genomic coordinates."))` stops Tab 3 explicitly (Tab 2 would instead fail at the "no overlap" check, since `nrow(out)` would be 0). |
| Invalid/negative thresholds entered | No explicit server-side range validation beyond Shiny's own `numericInput(min=, max=)` soft constraints on some (not all) fields; a manually typed out-of-range value is not defensively re-clamped or rejected by this module's own code. |
| No rows pass a filter (any tab) | Every filtering `eventReactive` ends in a `validate(need(nrow(...) > 0, "..."))` with a specific, actionable message — consistently implemented across all four analysis tabs. |
| Only one module exists in the data | Tab 3 still runs the loop over `mods` (length 1); Fisher's test and FDR correction both function correctly on a single-module input (`p.adjust` on a length-1 vector is a no-op, which is mathematically correct). |
| A selected module has no candidates (Tab 4) | `validate(need(nrow(df) > 0, "No candidate CpGs remain after applying these filters. Try relaxing a threshold above."))`. |
| User changes filters repeatedly without re-running | Each analysis tab's "has run" flag (`ov_has_run`, `modov_has_run`, `cand_has_run`) is **not** reset when a filter input changes — only on `load_btn` (all three) or the tab's own run button (and, for `cand_has_run`, also on `overlap_run_btn`). A previously computed table therefore remains visible, unchanged, after the user edits a filter but before clicking "Run" again, with no visual "results are stale" indicator. This is a deliberate simplification, explicitly flagged to the user only on Tab 1 ("These filters take effect once you run an analysis — nothing here computes anything by itself.") but not repeated on Tabs 2–4 themselves. See §19 (MEDIUM). |
| User opens Tab 2/3/4/5 before loading data | Each tab's `renderUI` checks `has_loaded()` first and shows `"Load a data source on the 'Data & Filters' tab first."` instead of erroring. |
| User opens Tab 4 before running Tab 2 | Explicit guidance message, as noted in §11. |

## 15. Complete Function Inventory

### Shiny reactive/module functions (custom)
- `mod_methyl_candidates_ui(id)` — builds the 5-tab UI shell.
- `mod_methyl_candidates_server(id, dataset, results = NULL)` — the module server; `dataset` argument is accepted but unused (§17.3).
- `register_viz(id_prefix, requires, plot_fn, filename)` — shared factory that wires one plot's button, gating `reactive()`, `renderPlot()`, and `downloadHandler()` in a single call, replacing six near-duplicate blocks.

### Data intake / standardization functions (custom)
- `mcd_find_col(cols, patterns)` — first-match, case-insensitive regex column finder.
- `mcd_standardize_module_assign(df)` — validates/extracts `cpg`, `module`, optional `kme`.
- `mcd_standardize_dmr(df)` — validates/extracts `chr`, `start`, `end`, plus optional `dmr_id`/`dmr_fdr`/`dmr_pvalue`/`delta_beta`/`n_cpgs`/`gene`/`direction`.
- `mcd_standardize_annotation(df)` — validates/extracts `cpg`, `chr`, `pos`, plus optional `gene`/`island`/`feature`.
- `mcd_standardize_cpg_stats(df)` — validates/extracts `cpg` plus at least one of `p_value`/`fdr`/`delta_beta`.
- `mcd_champ_full_annotation()` — loads and caches `ChAMPdata::probe.features` as a coordinate/gene/island/feature reference.
- `mcd_norm_chr(x)` — normalizes chromosome-name strings to a consistent `"chr"`-prefixed form.

### Filtering / overlap computation functions (custom)
- `mcd_filter_dmrs(dmr, fdr_max, p_max, dbeta_min, mincpgs_min, direction, chr_restrict)` — applies every DMR-side filter in one pass.
- `mcd_compute_overlap(module_assign, annotation, dmr_filtered, flank)` — the core coordinate-overlap engine: inner-joins CpGs to coordinates, builds `GRanges` for both CpGs and (optionally flanked) DMRs, runs `findOverlaps()`, and returns both the joined pairs and the full coordinate-resolved CpG universe.

### Display / formatting functions (custom)
- `mcd_pretty(df)` — renames a whitelist of internal column names to display labels via `MCD_PRETTY_MAP`, dropping anything not in the map.
- `mcd_detected_list(detected)` — renders the "Detected column mapping" bullet list.
- `mcd_module_colors(mods)` — chooses either literal WGCNA colour names or a fallback categorical palette for module-colored plots.
- `mcd_viz_block(ns, id_prefix, title, icon_name, btn_label, desc)` — one visualization card's static UI shell.

### Plotting functions (custom, all `ggplot2`-based)
- `mcd_plot_module_bar(mod_tab)` — overlapping CpGs per module, horizontal bar.
- `mcd_plot_enrichment_heatmap(mod_tab)` — single-row `-log10(p)` heatmap across modules.
- `mcd_plot_dmr_bar(overlap_tab, n = 20)` — top-N DMRs by overlapping CpG count.
- `mcd_plot_direction(tab)` — hyper/hypo bar chart.
- `mcd_plot_annotation_dist(tab)` — faceted genomic-region / CpG-island distribution.
- `mcd_plot_candidate_volcano(tab, highlight_cpgs = NULL)` — Delta-Beta vs. `-log10(significance)` scatter, optional candidate highlighting.

### Important package/base R functions used
| Function | Package | Role in this module |
|---|---|---|
| `GenomicRanges::GRanges()`, `GenomicRanges::findOverlaps()` | GenomicRanges | Builds interval representations of CpGs and DMRs; computes the core coordinate overlap. |
| `IRanges::IRanges()` | IRanges | Builds the start/end interval component of each `GRanges` object. |
| `S4Vectors::queryHits()`, `S4Vectors::subjectHits()` | S4Vectors | Extracts matched-pair indices from the `findOverlaps()` result (`Hits` object). |
| `stats::fisher.test()` | stats | One-sided exact test of module enrichment for DMR overlap (Tab 3). |
| `stats::p.adjust(method = "BH")` | stats | Benjamini-Hochberg FDR correction across per-module p-values (Tab 3). |
| `stats::aggregate()` | stats | Counts unique CpGs per DMR for `mcd_plot_dmr_bar()`. |
| `data.table::fread()` | data.table | Fast delimited-file parsing for all uploaded/preloaded tables. |
| `merge()` (base) | base | Inner-joins module assignment to annotation (coordinates) and overlap table to CpG-level stats. |
| `validate()` / `need()` (shiny) | shiny | User-facing, non-crashing guard clauses at every stage where an empty/invalid intermediate result would otherwise propagate. |
| `req()` (shiny) | shiny | Used inside `candidate_result()`/`register_viz()`'s `plot_obj` to require an upstream reactive before proceeding. |
| `eventReactive()` / `reactive()` / `reactiveVal()` / `observeEvent()` (shiny) | shiny | The reactive graph itself — button-gated computation and has-run state tracking throughout. |
| `DT::renderDataTable()` / `DT::datatable()` / `DT::formatSignif()` (DT) | DT | All four result tables, with column-aware significant-figure formatting. |
| `downloadHandler()` (shiny) | shiny | Every CSV/PNG export. |
| `ggplot2::ggplot()` and geoms/scales/themes | ggplot2 | All six visualizations. |
| `ggplot2::ggsave()` | ggplot2 | PNG export for every plot. |
| `utils::write.csv()` | utils | CSV export for every table. |
| `grDevices::col2rgb()` | grDevices | Tests whether module names are valid R colour names, inside `mcd_module_colors()`. |
| `tools::toTitleCase()` | tools | Cosmetic capitalization of the sex stratum in the load-summary text. |

### File/download functions
- Four `downloadHandler()` instances for CSV export (`overlap_download`, `modoverlap_download`, `cand_download`, plus the implicit per-plot ones), and six PNG `downloadHandler()`s registered inside `register_viz()` — ten download handlers total.

## 16. Package Inventory

| Package | Functions used in this module | Purpose here |
|---|---|---|
| **GenomicRanges** | `GRanges()`, `findOverlaps()` | Core genomic-coordinate overlap engine — the single most important scientific dependency of this module. |
| **IRanges** | `IRanges()` | Interval construction for `GRanges`. |
| **S4Vectors** | `queryHits()`, `subjectHits()` | Reading matched pairs out of the `Hits` object `findOverlaps()` returns. |
| **stats** | `fisher.test()`, `p.adjust()`, `aggregate()`, `setNames()` | Module enrichment testing, multiple-testing correction, per-DMR CpG counting. |
| **ChAMPdata** | `data("probe.features")` | Illumina methylation array coordinate/gene/CpG-island/genomic-feature reference annotation (preloaded path only). |
| **data.table** | `fread()` | Fast table parsing for both uploaded and preloaded files. |
| **DT** | `renderDataTable()`, `datatable()`, `formatSignif()` | All interactive result tables. |
| **ggplot2** | `ggplot()`, `geom_col()`, `geom_tile()`, `geom_point()`, `facet_wrap()`, `scale_fill_manual()`, `scale_fill_gradient()`, `scale_color_manual()`, `labs()`, `element_text()`, `ggsave()` | All six visualizations and their PNG export. |
| **shiny** | `moduleServer()`, `NS()`, `tabsetPanel()`/`tabPanel()`, `radioButtons()`, `fileInput()`, `numericInput()`, `selectizeInput()`, `checkboxInput()`/`checkboxGroupInput()`, `conditionalPanel()`, `actionButton()`, `downloadButton()`/`downloadHandler()`, `validate()`/`need()`, `req()`, `reactive()`/`reactiveVal()`/`eventReactive()`/`observeEvent()`, `renderUI()`/`uiOutput()`, `updateTabsetPanel()` | The entire reactive/UI scaffold. |
| **shinycssloaders** (`withSpinner`, implied by app-wide convention; called directly in this file) | `withSpinner()` | Loading spinners on the load-summary panel and every plot. |
| **grDevices** | `col2rgb()` | Detecting whether a module name is a valid literal colour, for direct-colour plotting. |
| **tools** | `toTitleCase()` | Cosmetic text formatting only. |
| **utils** | `write.csv()`, `head()` | CSV export; truncating the top-N DMR list for the bar plot. |

Not listed: packages loaded globally elsewhere in the app (e.g. `shinydashboard` for `valueBox()`, used here but defined/loaded app-wide) that contribute only generic UI chrome rather than this module's specific scientific logic.

## 17. Code-to-Biology Mapping

| Code Operation | Computational Meaning | Biological Meaning |
|---|---|---|
| `mcd_standardize_module_assign()` / `_dmr()` / `_annotation()` | Detects and coerces heterogeneous input tables into a fixed internal schema | Ensures WGCNA-derived module labels and DMRcate/comb-p/bumphunter-derived region calls (whatever the upstream tool) can be compared on a common footing |
| `mcd_filter_dmrs()` | Row-subsets the DMR table on FDR/p/effect-size/CpG-count/direction/chromosome | Restricts the analysis to DMRs meeting the analyst's chosen statistical-confidence and biological-relevance bar |
| `GenomicRanges::findOverlaps()` (in `mcd_compute_overlap()`) | Finds all `(query, subject)` index pairs whose genomic intervals intersect | Identifies which specific CpGs physically lie inside (or near) which specific differentially methylated regions |
| Module-selection / min-size / min-`|kME|` filtering | Row-subsets the CpG-module table | Restricts to co-methylation modules (and CpGs with strong intramodular connectivity) considered biologically coherent enough to trust |
| Direction derivation + `direction_consistency` | String comparison of independently-derived `hyper`/`hypo` labels | Checks whether a CpG's own methylation change agrees in direction with its overlapping DMR's aggregate change — a basic internal-consistency sanity check |
| `stats::fisher.test(..., alternative = "greater")` per module | One-sided exact test on a 2x2 contingency table | Tests whether a WGCNA module's CpGs are enriched (not merely present) among DMR-overlapping CpGs, beyond what the module's size alone would predict |
| `stats::p.adjust(method = "BH")` | Multiple-testing correction across modules | Controls the false-discovery rate when scanning many modules for enrichment simultaneously |
| Candidate scoring (`candidate_score`) | Weighted/unweighted additive combination of `-log10(FDR)`, `10*|Delta-Beta|`, `|kME|`, consistency bonus | A heuristic prioritization of which overlap-supported CpGs are most likely to be biologically important, by the analyst's chosen emphasis |
| `mcd_plot_annotation_dist()` | Frequency tabulation of `feature`/`island` columns | Shows where candidate CpGs fall relative to CpG islands and gene-structural context (promoter, body, shore, etc.) — standard methylation-biology interpretive context |

## 18. Data Validation and Error Handling

- Every multi-step `eventReactive()` (`loaded`, `overlap_run`, `module_overlap_run`, `candidate_result`) uses `validate(need(...))` at each point where an empty or invalid intermediate result would otherwise propagate silently — this is applied consistently across all four analysis stages.
- File-parse failures (`data.table::fread()` throwing) are caught with `tryCatch(..., error = function(e) NULL)` and converted into a `validate()`-driven user message rather than an application error.
- Missing optional data (no kME, no DMR FDR, no CpG stats, no annotation) degrades gracefully at every stage: the relevant UI control is either omitted (with an `empty-note` explanation) or the relevant filter/scoring term is simply skipped, never crashing.
- What is **not** validated: numeric threshold ranges beyond Shiny's own soft `min=`/`max=` widget constraints (no explicit server-side re-check); the count of rows silently dropped during standardization (duplicates, missing coordinates) is not reported as a number anywhere in the UI, only as a pass/fail "notes" message in specific failure cases.

## 19. Scientific Audit

- **Genome-build assumption is implicit and unchecked (HIGH).** `mcd_norm_chr()` normalizes chromosome *names* but nothing in this module checks or asserts a genome build (e.g. hg19 vs. hg38) for the DMR table, the annotation table, or the module-assignment table. On the preloaded path, coordinates come from `ChAMPdata::probe.features` (a fixed build tied to the installed package version) merged against a DMR table produced by an upstream pipeline that is assumed, but not verified here, to use the same build. On the upload path, a user could supply a DMR table on a different genome build than the annotation table with no cross-check — silently producing biologically meaningless overlaps (coordinates that happen to numerically intersect on different builds). No liftover, no build-tag column, no warning.
- **Fisher's exact test independence assumption (MEDIUM).** As noted in §10, treating each CpG as an independent trial ignores known spatial/network correlation between neighboring CpGs and within co-methylation modules by construction — a standard limitation of this class of enrichment test, not remedied here (e.g., no permutation-based null, no LD/region-block correction).
- **One-sided-only enrichment test (LOW/MEDIUM).** Depletion cannot be detected; whether this matters depends on whether the downstream user cares about modules with *fewer* DMR-overlapping CpGs than expected.
- **Combined score is an unnormalized, unweighted heuristic sum (MEDIUM).** As detailed in §11, `-log10(p)`-scale terms can dominate `10*|Delta-Beta|`/`|kME|`-scale terms in `combined` mode; the module does not claim statistical calibration for this score, but the UI does not warn the user of the scale imbalance either.
- **No formal joint-evidence statistic.** The module correctly does not invent a combined p-value/statistic for "module-and-DMR-supported" CpGs; it reports the two lines of evidence (overlap membership, module enrichment) separately and lets the user combine them via the heuristic score. This is a defensible scope boundary, explicitly reflected in the code, not an oversight.

## 20. Technical / Reactive Audit

- **Reactive dependency structure is correct and consistently button-gated**: nothing recomputes on every keystroke; all four analysis stages are explicit `eventReactive()`s tied to their own action buttons, as intended and as stated in the Tab 1 UI text.
- **Stale-result risk (MEDIUM, restated from §14):** none of `ov_has_run`, `modov_has_run`, `cand_has_run` reset when a filter input changes (only on `load_btn`/relevant run button/`overlap_run_btn` as applicable) — a previously rendered table or plot can visually persist after filters change but before the tab is re-run, with no "results out of date" indicator anywhere except the one general sentence on Tab 1.
- **No circular dependencies** were found; the reactive graph is a strict DAG (Tab 1 -> {Tab 2, Tab 3} -> {Tab 4, Tab 5}).
- **`outputOptions(..., suspendWhenHidden = FALSE)`** is applied to all three `DT::renderDataTable()` outputs (`overlap_table`, `modoverlap_table`, `cand_table`), ensuring they render even while their containing tab isn't the active one — correct usage for a multi-tab layout where a user might switch tabs before a table has ever been visible.
- **The overlap join runs twice** (once each in `overlap_run()` and `module_overlap_run()`) whenever both Tab 2 and Tab 3 are used, with no shared cache — a performance, not correctness, concern; on the array sizes this app's preloaded data uses (450K/EPIC scale probe counts, one CSV per sex), this is unlikely to be a practical bottleneck, but would not scale gracefully to substantially larger uploaded CpG sets.

## 21. UI–Server Consistency Audit

- **Visualization gating text matches server-side gating exactly** for all six plots (§12) — no mismatch found.
- **The dynamic scoring-formula preview (`cand_formula_ui`) matches the actual scoring code** in `candidate_result()` — both independently check the same column-presence conditions per `cand_rank_mode`, and stay in sync because they're driven by the same table (`overlap_run()$table`), not duplicated hard-coded lists.
- **Module-selection filters are declared once, on Tab 1, but consumed differently by Tabs 2 and 3** (all filters in Tab 2; DMR-side filters + grey-exclusion only in Tab 3) — this is explicitly documented in an in-code comment (lines 675–678) and in Tab 3's own UI description text ("using the same DMR filters as the overlap tab, but every module rather than only a selected one"), so this is a *disclosed* design choice, not a silent mismatch.
- **Cross-module results-key mismatch (HIGH — the most significant finding in this audit).** At the end of the server function (lines 958–967), this module writes its summary to the shared `results` reactiveValues object as:
  ```r
  results$candidate_cpgs <- list(
    n_overlap_cpgs = r$n_overlap_cpgs, n_overlap_dmrs = r$n_overlap_dmrs,
    n_dmr_passing = r$n_dmr_passing, source = loaded()$source
  )
  ```
  But `mod_methyl_candidates_config$id` is `"candidates"` (line 367), and the Methylomics-wide Assistant/ArthOChat context builder, `build_mx_context()` in `submodules_registry.R:163–181`, reads this module's results by **that config id**:
  ```r
  ids <- focus_id %||% vapply(MX_MODULES, function(m) m$config$id, character(1))
  for (mid in ids) {
    lines <- c(lines, .format_results_block(MX_MODULES_BY_ID[[mid]]$config$title, methyl_results[[mid]]))
  }
  ```
  i.e. it looks up `methyl_results[["candidates"]]`, not `methyl_results[["candidate_cpgs"]]`. Because the key this module writes (`candidate_cpgs`) never matches the key the context builder reads (`candidates`), **`methyl_results[["candidates"]]` is always `NULL`**, and `.format_results_block()` always renders this sub-module's block as *"NOT YET RUN IN THIS SESSION — no numbers exist for this sub-module yet"* to the ArthOChat assistant — even after a user has fully run the DMR-CpG Overlap analysis and populated real counts. This is a genuine, code-verifiable UI/cross-module integration bug: the summary result is computed correctly but is unreachable by the one consumer that reads it by name. No other place in the codebase reads `results$candidate_cpgs` (confirmed by a full-repository search) — the value written is, as of this audit, **write-only, unused output**.

## 22. Implemented vs. Intended Functionality

**Actually implemented:**
- Column-auto-detected data intake for both preloaded and uploaded module-assignment/DMR/annotation/CpG-stats tables.
- Coordinate-based CpG-DMR overlap via `GenomicRanges::findOverlaps()`, with an "inside" vs. "flanked" overlap mode.
- Per-module one-sided Fisher's exact enrichment test with BH-FDR correction.
- Multi-criterion candidate filtering and a four-mode additive scoring/ranking system.
- Six gated visualizations with independent PNG export.
- CSV export at every table stage.

**Partially implemented:**
- Direction-consistency checking — implemented, but only ever a simple string-equality of two independently derived direction labels, not a formal concordance statistic.
- The "combined" candidate score — implemented, but unnormalized across heterogeneous scales (§19).

**UI-present but not implemented (or implemented elsewhere, not here):**
- Nothing was found in this file's UI that lacks a corresponding server-side implementation — every button, filter, and conditional panel traced to working reactive/render logic. (Contrast: the *cross-module* Assistant integration is implemented but broken, per §21 — that is a connectivity bug, not a missing feature.)

**Code-present but unused / unreachable:**
- `dataset` function argument to `mod_methyl_candidates_server()` — accepted, never referenced (by design, per the header comment; §3).
- `results$candidate_cpgs` — written, never read by name anywhere in the repository (bug; §21).

No other unused custom function was identified; every `mcd_*` helper defined in this file is called at least once (§15).

## 23. Reproducibility

To reproduce this module's output independently of the running Shiny app, one needs:

- **Packages:** GenomicRanges, IRanges, S4Vectors, ChAMPdata (preloaded path only), stats, data.table, DT, ggplot2, shiny (+ shinycssloaders, shinydashboard for UI chrome only).
- **Input tables:** a CpG-to-module assignment table (CpG ID + module/color columns, optional kME) and a DMR table (chromosome + start + end, optional ID/FDR/p/Delta-Beta/n-CpGs/gene/direction); a coordinate/annotation table if the module-assignment table lacks chromosome/position.
- **Preloaded-path specifics (not reproducible from this file alone):** the exact contents of `METH_WGCNA_DIR/module_assignment_{sex}_merged10.csv`, `METH_DMR_DIR/dmr_{sex}_full.csv`, and `load_default_dmp("sva", sex)`'s output live outside this repository's inspected scope (external `METH_DATA_ROOT`) — this module reads them as static files and applies no further transformation to them beyond standardization/filtering, so their upstream generation (WGCNA clustering parameters, DMR-calling method/thresholds) is **not** documented here and must be sourced from the DMR/WGCNA sub-modules' own methods documentation.
- **Parameters:** every filter/threshold described in §5–§13 must be recorded to reproduce a specific run's candidate list — the module itself keeps no run-history log beyond the current session's reactive state.
- **Genome build:** not recorded or asserted anywhere in this module (§19) — must be independently confirmed to match between the DMR source and the annotation source for results to be biologically valid.

## 24. Audit Findings and Recommendations Summary

**HIGH**
1. **Assistant/ArthOChat results-key mismatch** (§21): `results$candidate_cpgs` is written but `build_mx_context()` reads `methyl_results[["candidates"]]` — the assistant can never see this module's computed summary. *Recommendation (not applied, per scope): rename the written key to `results$candidates` to match `mod_methyl_candidates_config$id`.*
2. **No genome-build validation** (§19): DMR, annotation, and module-assignment coordinate sources are merged with no check that they share a genome build; a mismatch would silently produce meaningless overlaps. *Recommendation: surface a build assumption/label in the UI, or a spot-check heuristic.*

**MEDIUM**
1. Fisher's exact test treats CpGs as independent, ignoring known spatial/network correlation (§16, §19) — a standard, undisclosed caveat of this enrichment approach.
2. "Combined" candidate score sums unnormalized, differently-scaled terms, risking significance-term dominance (§11, §19).
3. Result tables can go stale relative to filter inputs without any in-tab visual warning once a filter changes after a run (§14, §20).
4. `n_cpg_tested` in Tab 2's summary counts pre-coordinate-resolution CpGs, not post-join CpGs, understating how many were actually excluded for unresolved coordinates (§9).
5. Silent row-dropping during table standardization (duplicate IDs, missing coordinates) is not quantified in the UI (§5, §14).
6. Only enrichment (not depletion) is tested by the module-level Fisher's exact test (§10).

**LOW**
1. `1e-300` vs. `.Machine$double.xmin` floor inconsistency between the candidate score and the enrichment heatmap (§11, §13).
2. `modoverlap_table` sets column headers manually instead of via `mcd_pretty()`/`MCD_PRETTY_MAP`, unlike Tabs 2 and 4 (§10).
3. Top-20 DMR cap in `mcd_plot_dmr_bar()` is a code-level default, not user-configurable (§13).
4. `mcd_module_colors()`'s 7-color fallback palette recycles for modules beyond 7 non-WGCNA-named modules (§13).
5. "Inside" and "flanked" overlap modes are functionally identical at the default `flank = 0` (§13).

---

## 25. Thesis Implementation Summary

The Candidate CpGs (Module-DMR Overlap) submodule comprises five sequential tabs — Data & Filters, DMR-CpG Overlap, Module-DMR Overlap, Candidate CpGs, and Visualization — that integrate a WGCNA co-methylation module assignment table with a differentially methylated region (DMR) results table, either from this application's own preloaded sex-stratified methylomics pipeline or from user-uploaded files with automatically detected columns. The implementation uses `GenomicRanges::findOverlaps()` to resolve which module-assigned CpGs physically lie inside (or near) which significant DMRs, `stats::fisher.test()` with Benjamini-Hochberg correction to test each WGCNA module for statistical enrichment of DMR-overlapping CpGs against the tested CpG universe, and a user-configurable additive scoring formula to filter and rank the resulting CpG-DMR overlap pairs into a prioritized candidate list, producing downloadable overlap tables, module-enrichment statistics, ranked candidate tables, and six diagnostic plots as its principal outputs. This analysis provides a coordinate- and network-integrated approach for identifying CpG sites supported by two independent lines of evidence — region-level differential methylation and co-methylation network structure — intended to narrow a genome-wide DMR/module output into a shorter, prioritized set of candidates for downstream biological interpretation or biomarker-panel construction; it does not itself perform WGCNA clustering, DMR calling, or any joint statistical test of combined evidence, and (per this audit) its summary output is not currently reaching the application's Assistant module due to a results-key mismatch between this file and the shared cross-module context builder.
