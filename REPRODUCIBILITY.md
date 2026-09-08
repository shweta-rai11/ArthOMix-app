# ArthOMix reproducibility status

This document tracks, honestly and per-pipeline-stage, which of ArthOMix's precomputed/headline
tables can be independently regenerated from data and code that actually exist in this repository.

**Background.** A 2026-09-07 hostile audit found that ~90-97% of the app's precomputed thesis
tables (methylomics pipeline stages, all cross-omics tables, all multiomics tables/DIABLO fits,
transcriptomics MR/coloc/WGCNA/FS objects) had no generating script anywhere in the repo or git
history — the original generating pipeline was deliberately kept in a separate, unlinked repo.
Only the transcriptomics blood DEG table had been proven reproducible (in a prior session).

**Reproducibility bar used here.** The original pipeline's exact scripts and hyperparameters are
gone; only its outputs and (for most stages) the live app's own processing functions remain. So
this work targets **statistically-equivalent reproduction** — running the app's own documented
methodology on the real raw/source data and reporting how closely the result matches the
precomputed table (effect-size correlation, direction agreement, overlap of significant
calls) — not byte-exact reproduction. Where a match is weak, that is reported as a finding, not
hidden.

Scripts live in `ArthOMix/reproduce/{methylomics,multiomics}/`. Each script's header documents
every parameter choice and its justification (a `preloaded`-mode default hardcoded in the live app
code, where one exists; otherwise the function's own generic default, flagged as lower-confidence).
Script outputs (comparison CSVs, downloaded raw GEO data) are gitignored — regenerate by re-running.

## Methylomics (GSE42861, 689-sample whole-blood cohort)

Raw data (`data/preloaded/methylomics/matrix/beta_raw.rds`, 2.0GB beta matrix + `pheno.rds`) is
genuinely bundled in the repo. Four of six flagged stages have `preloaded`-mode hyperparameter
defaults hardcoded directly in the live app's module code — meaning the live app itself asserts
what methodology was used, not just this write-up.

| Stage | Script | Result | Confidence |
|---|---|---|---|
| QC (PCA outliers, sex-check) | `01_qc.R` | PCA outlier flags 98.5% agree with precomputed. Sex-check (chrY_mean/sex_mismatch) is **not reproducible**: `beta_raw.rds` as bundled has zero chrX/chrY probes (already stripped upstream of what's bundled) — a genuine, disclosed gap, not a bug in this script. | High (for what could be tested) |
| DMP (SVA + limma + bacon) | `02_dmp.R` | Δβ correlation r=1.00, t-statistic correlation r=0.96 (F) / 0.98 (M). Zero false-positive significant CpGs in either sex; 11/18 (F) and 0/0 (M) significant calls recovered exactly. | **High** |
| DMR (DMRcate on DMP output) | `03_dmr.R` | 64.2% (F) / 62.9% (M) of precomputed regions overlap a regenerated region; 57% of precomputed significant DMRs recovered in both sexes, zero spurious significant calls. | **High** |
| WGCNA (co-methylation network) | `04_wgcna.R` | 27 modules found both sexes. Adjusted Rand Index vs precomputed module assignment: 0.573 (F), 0.831 (M). | **High** |
| ML Feature Selection (LASSO/RF) | `05_featureselection.R` | **Weak / honest negative.** LASSO selection-frequency correlation r=0.20 (F) / 0.05 (M). RF importance: not comparable — the precomputed table has only 9-12 rows, close to the DMP significant-CpG count, revealing the real candidate pool was almost certainly the DMP-significant panel, not a fresh top-variance selection like this script (and the live UI's generic default) assumed. | Low — see note below |
| Diagnostic Classifier | not attempted | Same structural problem as Feature Selection (no `preloaded`-mode hyperparameter defaults in the code — confirmed by direct grep, zero matches), so a generic-default attempt would carry the same low confidence. Skipped rather than repeat a known-weak pattern. | Not attempted |

**Why Feature Selection/Diagnostic Classifier are different from DMP/DMR/WGCNA:** those three have
literal `if (isTRUE(methyl_dataset$preloaded)) <value> else <other>` branches in the module code
that pin the exact original hyperparameters (SVA topn, WGCNA power vector, min module size, etc.).
Feature Selection and Diagnostic Classifier have **no such branches** — every hyperparameter (LASSO
alpha/nfolds, RF ntree/mtry, candidate-CpG universe) is a bare UI default with no code evidence it
matches what the original pipeline used. This is a structural fact about the codebase, not a
one-off judgment call, and the weak result above validates the concern rather than being a surprise.

## Transcriptomics

- Blood DEG (differential expression): previously verified byte-for-byte reproducible from raw
  counts via the app's own live limma path (prior session, not re-verified here).
- `combined_expr_batchcorrected.rds` (feeds WGCNA) and its downstream WGCNA: raw data (GSE93272 +
  GSE110169) and a live merge/ComBat code path both exist, but — like Feature Selection above —
  the batch-correction module has **zero** `preloaded`-conditional hyperparameter defaults (batch
  column, protect-columns, quantile-normalization method, ComBat vs limma vs SVA, prior type,
  reference batch are all bare UI defaults). Not attempted this session for the same reason.
- MR/coloc: exposure is bundled eQTLGen cis-eQTL data referenced by the live module but not
  actually bundled in `data/`; outcome is Okada et al. 2014 RA GWAS, partially fetchable via the
  app's existing `ieugwasr` dependency (only used for LD-clumping today, not full association
  fetch). Regenerable-with-external-fetch, not attempted this session.

## Cross-Omics

- The gene-level min-p-across-43-tests multiple-testing artifact (finding #1 of the 2026-09-07
  audit) was fixed this session in `crossomics_biomarkerconv_helpers.R` — see git history, not a
  reproduce/ script (it's a correctness fix, not a reproducibility exercise).
- The MR-stage module (`mod_cross_mr_stage.R`) is **confirmed genuinely not-regenerable**: its own
  header comment states it is a viewer for externally-run GoDMC (cis-mQTL) → Ishigaki et al. 2022
  RA GWAS (GCST90132223) MR results and deliberately never runs MR itself. GoDMC is access-gated
  consortium data, not a public download; this is an honest, permanent limitation, not a gap to close.

## Multi-Omics (Tao et al. 2021 anti-TNF cohort)

**The GEO accessions were not recorded anywhere in the repo or git history** (confirmed by
`git log --all -p`) — only the free-text label "Tao et al. 2021" existed. This session identified
and ground-truth-verified (via NCBI E-utilities, not just search-engine text) both halves of the
study:

- **GSE138746** — RNA-seq, 240 samples (80 patients × PBMC/monocyte/CD4+T cells)
- **GSE138653** — DNA methylation array (Illumina EPIC, GPL21145), ~80 samples

These are two independent GEO series (not a linked SuperSeries — that's why the app's own
`multi_geo_autosplit_fetch()` auto-detection finds nothing), and their public metadata carries no
shared subject-ID column, so exact patient-level linkage between the two series' *public* metadata
is not possible from GEO alone. Patient-matched data for 79 patients already existed in
`data/examples/multiomics_upload/` (`gse138746_pbmc_rnaseq_counts.csv`,
`gse138653_methylation_beta_top2000.csv`, `gse138747_sample_metadata.csv`) from a prior session
this one has no other record of; it was independently legitimacy-checked here (realistic RNA-seq
count and methylation beta-value ranges, drug counts 37 ADA/42 ETN vs. the paper's reported 38/42)
before use, but carries no verified sex labels, so a **pooled (non-sex-stratified)** test was run
instead of reproducing the exact sex-stratified `diablo_drugsex_*` fits.

`reproduce/multiomics/01_pooled_diablo_response.R` ran the app's own leakage-safe nested-CV DIABLO
engine (`mss_nested_cv()`, in its built-in pooled mode) on this freshly-fetched data, predicting
EULAR responder-vs-non-responder status from paired expression + methylation:

**Result: AUROC 0.483 [0.353, 0.613] — chance level.**

This independently confirms — from raw public data, through a different code path than the one
that produced the bundled tables — the audit's most consequential finding: the multi-omics
drug-response prediction does not recover meaningful signal, far below the published 0.79-0.89
benchmark. It is not an artifact of the app's own precomputed pipeline; it reproduces from scratch.

## Summary

| Area | Regenerable & verified | Honest negative / lower confidence | Confirmed not-regenerable |
|---|---|---|---|
| Methylomics | QC (partial), DMP, DMR, WGCNA | Feature Selection | Diagnostic Classifier (not attempted, same pattern) |
| Transcriptomics | DEG (prior session) | — | — |
| — not attempted this session | batch-correction/WGCNA, MR/coloc | | |
| Cross-Omics | min-p correctness fix (not a repro exercise) | — | MR-stage (GoDMC access-gated) |
| Multi-Omics | GEO accessions found + verified; pooled DIABLO independently confirms chance-level AUROC | sex-stratified fits (no verified sex linkage) | — |
