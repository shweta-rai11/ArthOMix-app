# `mod_coloc.R` — Colocalization Module: Teaching, Method, and Thesis-Documentation Notes

File: `ArthOMix/R/transcriptomics/mod_coloc.R` (285 lines)
Group: Genetics &nbsp;|&nbsp; Sibling module: `mod_mr.R` (Mendelian randomisation)
Prepared: 2026-08-25

This document is derived from the code in `mod_coloc.R`, the shared helpers it calls in
`global.R`, and the thesis methods reference `data/preloaded/transcriptomics/results/METHODS_2.7_colocalisation.md`
(the offline pipeline this live module is a simplified, interactive counterpart to).
Anything the live module does *not* do that the offline pipeline does is flagged
explicitly — the two are not the same analysis and should not be cited as if they were.

---

## 1. Why this module exists (the educational part)

### 1.1 The question colocalisation answers

`mod_mr.R` (the preceding module in the Genetics group) runs *cis*-Mendelian
randomisation: it uses a SNP that affects a gene's expression (an eQTL) as an
instrument to ask whether that gene's expression has a causal effect on rheumatoid
arthritis (RA) risk. This works only under an assumption MR does not itself test: that
the SNP driving the eQTL signal and the SNP driving the disease-risk signal at that
locus are **the same variant**.

That assumption fails whenever two nearby but distinct variants are correlated with
each other because they sit close together on the chromosome (linkage disequilibrium,
LD) — one variant genuinely alters the gene's expression, a different variant
genuinely alters disease risk, and LD makes them look associated with each other's
trait even though neither causes the other's effect. When this happens, MR returns a
confident, statistically significant, and **spurious** causal estimate. This is the
best-known weakness of *cis*-eQTL MR (Zhu et al., 2016; Wallace, 2020), and it is
exactly why a locus cannot be called "causal" on the strength of an MR p-value alone.

Colocalisation is the test that distinguishes the two scenarios:

- **One shared causal variant** — the eQTL and the disease association really do
  point at the same underlying variant. MR's assumption holds, and an MR estimate at
  this locus is trustworthy.
- **Two distinct causal variants in LD** — the eQTL and the disease association are
  driven by different variants that merely co-occur. MR's assumption is violated, and
  an MR estimate at this locus is confounded and should not be read as causal.

Overlap of two association peaks on a locus plot is **not** evidence for the first
scenario over the second — both scenarios can produce overlapping peaks. Colocalisation
replaces that visual judgement with a formal posterior probability.

### 1.2 The method: `coloc.abf` (Giambartolomei et al., 2014)

`coloc.abf` (R package `coloc`) takes GWAS-style **summary statistics** for two traits
over the same genomic region — no individual-level genotypes required — and evaluates
five mutually exclusive hypotheses about how many distinct causal variants exist in
that region and which trait(s) they affect:

| Hypothesis | Meaning |
|---|---|
| **H0** | No causal variant for either trait in this region |
| **H1** | A causal variant for trait 1 (the eQTL) only |
| **H2** | A causal variant for trait 2 (the GWAS/disease trait) only |
| **H3** | Both traits have a causal variant in the region, but **different** variants |
| **H4** | Both traits share **one** causal variant |

For each SNP, the method computes an approximate Bayes factor (Wakefield's ABF) from
its effect size and standard error alone (no LD matrix or individual genotypes are
needed for `coloc.abf` specifically — that is what makes it usable directly on GWAS
summary statistics). These per-SNP Bayes factors are combined, together with three
prior probabilities — p1 (a given SNP is causal for trait 1), p2 (causal for trait 2),
and p12 (causal for **both**) — into one posterior probability per hypothesis across the
whole region, `PP.H0` through `PP.H4`, which sum to 1.

**H3 is the failure mode that invalidates an MR estimate** (distinct causal variants
masquerading as one signal through LD); **H4 is the state of the world MR silently
assumes**. A high `PP.H4` is read as support for a single shared causal variant —
i.e., support for treating the gene's expression as genuinely on the causal path to
disease risk, not merely correlated with it through a nearby confounding variant.

**The method's core limitation**, inherited by every use of `coloc.abf` in this module:
it assumes **at most one causal variant per trait** in the analysed region. That
assumption is reasonable for a typical, narrow *cis*-eQTL window and becomes
unreliable in regions with multiple independent association signals — most notably the
MHC, which carries several independent RA signals (Raychaudhuri et al., 2012). Where
the assumption is violated, `PP.H3` is inflated (two traits each driven by several
variants can resemble "different causal variants" even when one of those variants is
truly shared), so a high `PP.H3` inside a multi-signal region like the MHC cannot be
read as clean evidence against colocalisation. The offline thesis pipeline addresses
this for MHC genes with a second, more complex method (`coloc.susie`, §7 below); **this
live module does not** — see §6.

---

## 2. Role in the app

- App: ArthOMix Explorer
  - Top-level omics area: **Transcriptomics** (`R/transcriptomics/`)
    - Group: **Genetics**
      - `mod_mr.R` — *cis*-MR (upstream: identifies candidate causal genes)
      - `mod_coloc.R` — **this file** (downstream: stress-tests whether a candidate
        gene's MR signal survives the shared-causal-variant test)

### Inputs it receives

- `id` — the Shiny module namespace string.
- `dataset` — the shared `reactiveValues` object passed to every transcriptomics
  module; **not read anywhere in this file** (colocalisation runs on its own bundled
  eQTL/GWAS summary statistics, not on the app's active expression matrix — the same
  pattern `mod_mr.R` follows).
- `results` — the shared `reactiveValues` object this module writes into
  (`results$coloc`, §5.3).
- On disk: `coloc_regions.rds` (`COLOC_REGIONS_RDS`, resolved in `data_paths.R:67` to
  `data/preloaded/transcriptomics/processed/new/coloc_regions.rds`) — a named list, one
  entry per gene, each holding a pre-extracted, pre-aligned pair of `eqtl` and `gwas`
  data frames for that gene's cis-window.
- User inputs: a data-source radio choice, a gene pick, a case-fraction slider, and
  (upload mode only) a GWAS file plus column-mapping dropdowns.

### Outputs it produces

- On-screen: a text summary, a posterior-probability bar chart, a two-track regional
  association plot, and a per-SNP results table with CSV download.
- `results$coloc$genes_tested[[gene]]` — one entry per completed run, written for
  consumption elsewhere in the app (the consumer is outside this file's scope, per the
  same documentation convention used for `mod_mr.R`).

---

## 3. What is bundled vs. what is live

The module header (`mod_coloc.R:1-21`) states this precisely, and it is the single
most important fact to get right when describing this module:

> The eQTL side is **always** the bundled cis-window instrument for the chosen gene, in
> both modes: `coloc_regions.rds` only has a prepared eQTL region for **33 genes**, and
> building an equivalent for an arbitrary gene would need the full eQTLGen dataset this
> app doesn't otherwise load. **Only the GWAS side is swappable.**

| | Bundled mode (`data_source == "project"`, default) | Upload mode (`data_source == "upload"`) |
|---|---|---|
| **eQTL side** | Bundled cis-window instrument for the picked gene (`coloc_regions.rds`) | *Same* bundled cis-window instrument — never changes |
| **GWAS side** | Bundled RA GWAS (Okada et al., 2014) | Any trait's GWAS summary statistics the user uploads |
| **Alignment needed?** | No — both sides were pre-aligned by the offline pipeline; a plain rsID intersection is safe | Yes — an arbitrary upload has no guaranteed allele coding relative to the bundled eQTL, so the module runs both sides through `TwoSampleMR::format_data()`/`harmonise_data()` (the same functions `mod_mr.R`'s own upload mode uses) to align alleles and flip β sign consistently before `coloc.abf()` |

This is why the gene picker (`selectInput(ns("gene"), ...)`, 33 choices) is shown and
required in **both** modes, and why only the GWAS-related controls sit inside the
`upload`-mode conditional panel.

---

## 4. Inputs (UI)

| Input | Widget | Choices / range | Default | Required in |
|---|---|---|---|---|
| `data_source` | `radioButtons` | `"project"` (Bundled RA GWAS), `"upload"` | `"project"` | Both |
| `gene` | `selectInput` (`selectize = FALSE`) | 33 genes = `sort(names(coloc_regions))` | first alphabetically | Both — this is always the eQTL side |
| `gwas_label` | `textInput` | free text, label only | `"Uploaded GWAS"` | Upload only |
| `gwas_file` | `fileInput` | `.csv`, `.tsv`, `.txt` | — | Upload only, **required** |
| Column-mapping (`gwas_snp`, `gwas_beta`, `gwas_se`, `gwas_pval`, `gwas_ea`, `gwas_oa`, `gwas_n`) | `selectInput`s generated by `gwas_col_map_ui(..., extra_fields = "n")` | column names of the uploaded file, auto-guessed by regex (`GWAS_COL_PATTERNS`) | best regex guess, or first column | Upload only, **all required** |
| `gwas_eaf` | `selectInput` | column names, or `"(none)"` | best guess or none | Upload only, **optional** |
| `case_frac` | `sliderInput` | 0.05 – 0.5, step 0.01 | 0.33 | Both — feeds `coloc.abf`'s case-control sample-proportion parameter (`s`), needed because the GWAS side is a case/control disease trait |

**Why `n` (sample size) is required here but not in `mod_mr.R`'s upload mode:**
`coloc.abf` needs a per-trait sample size `N` to compute its Bayes factors; `mod_mr.R`'s
Wald-ratio/IVW estimators do not. `gwas_col_map_ui()` is a shared helper — `mod_coloc.R`
is the module that actually requests the `n` field via `extra_fields = "n"`
(`global.R:1226-1228` documents this split explicitly).

**Validation:** `req()` guards on every required field before either run function
proceeds; `validate(need(...))` surfaces in-place (non-crashing) messages when fewer
than 10 SNPs are shared/complete after filtering, or when upload-mode harmonisation
finds no overlap at all.

---

## 5. Processing — how a run is performed

### 5.1 Dispatch

```r
coloc_result <- eventReactive(input$run_btn, {
  if (identical(input$data_source, "upload")) coloc_result_uploaded() else coloc_result_project()
}, ignoreInit = TRUE)
```
Fires only on "Run colocalisation"; the branch taken depends on `data_source` at click
time.

### 5.2 Bundled path — `coloc_result_project()` (`mod_coloc.R:120-147`)

```
r <- coloc_regions[[gene]]                     # pre-extracted eqtl + gwas data frames
common <- intersect(eqtl$rsid, gwas$rsid)        # rsID intersection — safe because both
                                                  # sides were pre-aligned by the offline pipeline
validate(need(length(common) >= 10, ...))
        ↓
ok <- complete.cases(beta, se, eaf, n) on both sides & se > 0 (both) & 0 < eaf < 1
validate(need(sum(ok) >= 10, ...))
        ↓
d1 <- list(beta, varbeta = se^2, N = round(median(n)), MAF = min(eaf, 1-eaf), type = "quant", snp)   # eQTL
d2 <- list(beta, varbeta = se^2, N = round(median(n)), type = "cc", s = case_frac, snp)                # GWAS
        ↓
res <- coloc::coloc.abf(dataset1 = d1, dataset2 = d2)     # package-default priors: p1=1e-4, p2=1e-4, p12=1e-5
        ↓
snp_df <- per-SNP table: snp, eqtl_beta, eqtl_p, gwas_beta, gwas_p, pos, snp_pp_h4 (SNP.PP.H4)
        ↓
list(gene, gwas_label = "Bundled RA GWAS (Okada 2014)", summary = res$summary, snp_df, n_snp, uploaded = FALSE)
```

Note `dataset1` is typed `"quant"` (expression is a quantitative trait) and `dataset2`
is typed `"cc"` (case/control — RA is a binary disease outcome), which is why `d2`
carries `s` (the assumed case fraction, from the `case_frac` slider) while `d1` carries
`MAF` instead.

### 5.3 Upload path — `coloc_result_uploaded()` (`mod_coloc.R:156-208`)

```
req(gene, gwas_file, all seven mapping inputs)
eqtl <- bundled cis-window region for `gene`  (unchanged — never swapped)
gwas_raw <- gwas_df_r()                        # cached fread() parse of the uploaded file
        ↓
exp_fmt <- TwoSampleMR::format_data(eqtl, type = "exposure",
             snp_col="rsid", beta_col="beta", se_col="se", pval_col="p",
             effect_allele_col="ea", other_allele_col="nea", eaf_col="eaf",
             samplesize_col="n", chr_col="chr", pos_col="position")
exp_fmt$exposure <- "<gene> eQTL (blood, bundled)"

out_fmt <- TwoSampleMR::format_data(gwas_raw, type = "outcome",
             snp_col/beta_col/se_col/pval_col/effect_allele_col/other_allele_col = user-mapped columns,
             eaf_col = mapped or "eaf", samplesize_col = gwas_n)
out_fmt$outcome <- gwas_label
        ↓
dat_up <- TwoSampleMR::harmonise_data(exp_fmt, out_fmt, action = 2)   # aligns alleles, flips β sign
validate(need(nrow(dat_up) > 0, "no overlapping SNPs..."))
dat_up <- dat_up[dat_up$mr_keep, ]                                    # drop unresolved palindromic SNPs etc.
        ↓
ok <- complete.cases(...) & se.exposure > 0 & se.outcome > 0 & 0 < eaf.exposure < 1
validate(need(sum(ok) >= 10, ...))
        ↓
d1 <- list(beta.exposure, varbeta.exposure, N, MAF, type="quant", snp)
d2 <- list(beta.outcome, varbeta.outcome, N, type="cc", s = case_frac, snp)
res <- coloc::coloc.abf(dataset1 = d1, dataset2 = d2)
        ↓
snp_df built in the SAME shape as the bundled path (pos = dat_up$pos.exposure), so every
downstream output (plots/table/download) needs no branching between the two modes.
        ↓
list(gene, gwas_label = user label, summary = res$summary, snp_df, n_snp = nrow(dat_up), uploaded = TRUE)
```

`action = 2` in `harmonise_data()` is `TwoSampleMR`'s standard palindromic-SNP handling:
attempt to infer strand from allele frequency for ambiguous A/T or C/G SNPs, and drop
what cannot be resolved — the same setting `mod_mr.R`'s own upload mode uses, so the
two modules' harmonisation behaviour cannot silently diverge.

### 5.4 What is common to both paths

- `coloc::coloc.abf()` is called with **package-default priors** in every run of this
  module: `p1 = 1e-4`, `p2 = 1e-4`, `p12 = 1e-5`. No prior-sensitivity re-run is
  performed here — see §6 for how this differs from the offline pipeline.
- `results$coloc$genes_tested[[gene]] <- list(n_snp, pp_h4 = round(PP.H4.abf, 3), gwas = gwas_label, uploaded)`
  is written on every completed run (`observeEvent(coloc_result(), ...)`,
  `mod_coloc.R:217-224`), keyed by gene, so repeat runs on the same gene overwrite that
  gene's entry rather than accumulating duplicates.

---

## 6. Outputs

| Output | Reactive | Content |
|---|---|---|
| **Summary text** | `output$summary_ui` | Upload-mode disclosure (if applicable); gene name + shared-SNP count + GWAS label; `PP.H4` as a percentage; a banded interpretation note — **PP.H4 > 0.8** "Strong support for colocalisation at this locus" (✓ icon), **0.5–0.8** "Moderate support" (i icon), **≤ 0.5** "Little support for a single shared causal variant" (⚠ icon) |
| **Posterior-probability plot** | `output$pp_plot` | Bar chart of `PP.H0`…`PP.H4`, H4 highlighted in red, others grey |
| **Regional association plot** | `output$region_plot` | Two stacked panels (`facet_wrap`, free y-scale), `-log10(p)` vs. genomic position, one panel for the eQTL track and one for the GWAS track — a visual complement to the numeric `PP.H4`, not a substitute for it |
| **SNP-level table** | `output$snp_table` (`DT`) | Per-SNP: `snp`, `eqtl_beta`, `eqtl_p`, `gwas_beta`, `gwas_p`, `pos`, `snp_pp_h4` (that SNP's own share of `PP.H4`) |
| **CSV download** | `output$download_coloc` | The same `snp_df`, filename `coloc_<gene>.csv` |
| **Cross-module** | `results$coloc$genes_tested` | `n_snp`, `pp_h4`, `gwas`, `uploaded`, keyed by gene |

---

## 7. Explicitly not implemented in this live module

The offline thesis pipeline (`METHODS_2.7_colocalisation.md`, §2.7.2–2.7.4) does
several things this interactive module deliberately does **not** reproduce, because
they require infrastructure (network queries to a summary-statistic server, an
LD-reference panel, per-region caching) beyond what a live, per-click Shiny run is
built for:

- **No dual-prior sensitivity check.** The offline pipeline runs every gene at both the
  package-default `p12 = 1e-5` and a conservative `p12 = 1e-6`, and only reports
  colocalisation support that holds at both. This module runs once, at the package
  defaults, and reports that single result.
- **No `coloc.susie` / multiple-causal-variant model.** The offline pipeline
  supplements `coloc.abf` with SuSiE-based colocalisation inside the MHC, where the
  single-causal-variant assumption is known to fail. This module always uses
  `coloc.abf` alone, for every gene, including any MHC gene among the 33 bundled
  regions — its `PP.H3`/`PP.H4` should be read with the same caution the offline
  methods document attaches to raw `coloc.abf` output in the MHC (§1.2 above).
- **No MHC flag or MHC-specific messaging** in the UI or summary text — unlike
  `mod_mr.R`, which explicitly flags and offers to exclude MHC instruments, this module
  is silent on whether the picked gene falls in the MHC.
- **No re-derivation of the eQTL region from raw summary statistics.** Both modes
  always read the same pre-extracted `coloc_regions.rds` entry for the chosen gene;
  there is no path to test an arbitrary, non-bundled gene's eQTL side.
- **No image/plot download** — only the SNP-level table has a `downloadButton`; the
  posterior-probability and regional-association plots are view-only.

---

## 8. Function reference

| Function | Location | Purpose |
|---|---|---|
| `mod_coloc_config` | `mod_coloc.R:23-28` | Static registration metadata (id, group, title, description, icon) |
| `mod_coloc_ui(id)` | `mod_coloc.R:30-93` | Builds the two-column layout: controls + result summary/PP plot (left), regional plot (right), SNP table (full width, post-run) |
| `mod_coloc_server(id, dataset, results)` | `mod_coloc.R:95-284` | All reactive/statistical logic (§5–§6) |
| `coloc_result_project()` | `mod_coloc.R:120-147` | Bundled-GWAS run (§5.2) |
| `coloc_result_uploaded()` | `mod_coloc.R:156-208` | Uploaded-GWAS run, with harmonisation (§5.3) |
| `gwas_df_r` (reactive) | `mod_coloc.R:114` | Cached parse of the uploaded GWAS file via `read_uploaded_table()` |
| `output$gwas_map_ui` | `mod_coloc.R:115` | Column-mapping UI via the shared `gwas_col_map_ui()` (`global.R:1229`), with `extra_fields = "n"` |

Shared helpers used (defined in `global.R`, also used by `mod_mr.R`):
`read_uploaded_table()`, `gwas_col_map_ui()`, `guess_gwas_col()`, `GWAS_COL_PATTERNS`.

---

## 9. Suggested paragraph for a thesis Methods/Results narrative

Use this as a starting point and adapt it to whichever section (interactive-tool
description vs. the offline §2.7 pipeline) you are writing into — they are not
interchangeable, and the paragraph below describes the **live app module**, not the
dual-prior/SuSiE offline pipeline.

> To determine whether a candidate gene's *cis*-eQTL association and its association
> with rheumatoid arthritis risk are likely driven by the same underlying variant —
> rather than by two distinct variants correlated through linkage disequilibrium, which
> would invalidate a Mendelian randomisation estimate at that locus — the ArthOMix
> Explorer application provides an interactive Bayesian colocalisation tool
> implementing `coloc.abf` (Giambartolomei et al., 2014). For a user-selected gene, the
> module always draws the eQTL side from a bundled, pre-extracted *cis*-window
> instrument (available for 33 candidate genes identified upstream by the Mendelian
> randomisation module), and tests it against either the bundled rheumatoid arthritis
> GWAS (Okada et al., 2014) or a GWAS of any other trait supplied by the user. Where a
> user-supplied GWAS is used, allele alignment and effect-direction harmonisation are
> performed live via `TwoSampleMR::harmonise_data()` before colocalisation, since an
> arbitrary upload carries no guaranteed allele coding relative to the bundled eQTL
> data. The tool reports posterior probabilities for the five standard colocalisation
> hypotheses (H0–H4), with particular emphasis on PP.H4 — the posterior probability
> that expression and disease risk share one causal variant — alongside per-SNP
> evidence and a regional association plot, allowing a user to assess, locus by locus,
> whether a gene prioritised by Mendelian randomisation withstands this
> shared-causal-variant test before it is described as a plausible causal candidate
> rather than merely a statistically associated one.

---

## References cited in this document

Giambartolomei, C., Vukcevic, D., Schadt, E.E., et al. (2014) 'Bayesian test for
colocalisation between pairs of genetic association studies using summary
statistics', *PLoS Genetics*, 10(5), e1004383.

Raychaudhuri, S., Sandor, C., Stahl, E.A., et al. (2012) 'Five amino acids in three
HLA proteins explain most of the association between MHC and seropositive rheumatoid
arthritis', *Nature Genetics*, 44(3), pp. 291–296.

Wallace, C. (2020) 'Eliciting priors and relaxing the single causal variant assumption
in colocalisation analyses', *PLoS Genetics*, 16(4), e1008720.

Zhu, Z., Zhang, F., Hu, H., et al. (2016) 'Integration of summary data from GWAS and
eQTL studies predicts complex trait gene targets', *Nature Genetics*, 48(5),
pp. 481–487.
