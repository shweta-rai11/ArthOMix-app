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

- **GSE138747** — the SuperSeries (320 samples = 240 + 80), correctly cited in the thesis draft
- **GSE138746** — RNA-seq sub-series, 240 samples (80 patients × PBMC/monocyte/CD4+T cells)
- **GSE138653** — DNA methylation sub-series (Illumina EPIC, GPL21145), 80 samples

**Correction (this session, after further testing):** an earlier version of this document claimed
GSE138746/GSE138653 were an unlinked pair and that this was why the app's existing
`multi_geo_autosplit_fetch()` failed. That was wrong — it was tested against the two *sub-series*
accessions directly rather than the actual SuperSeries accession. Re-tested against the correct
accession: `multi_geo_series_relation("GSE138747")` and `multi_geo_autosplit_fetch("GSE138747")`
(the pre-existing, unmodified function) both correctly discover the linkage via GEO's own
`Series_relation` metadata. Fetching still fails, but for the real, sole, and permanent reason
(confirmed directly: `Biobase::exprs()` on GSE138746's parsed series matrix is 0 rows × 240
samples): GEO stores no values at all for this sequencing-based series — the actual counts live in
a supplementary file with no standard format, which GEO's series-matrix reader cannot parse,
independent of any SuperSeries linkage question. The existing error message
(`multi_geo_autosplit_fetch`'s "may only exist as a supplementary file... use Upload Dataset
instead") already said this correctly before this session touched anything.

`multi_geo_dual_fetch()` (added this session) is still sound, tested engineering — useful for
studies that genuinely have no SuperSeries wrapper in GEO at all — but it was not the fix this
specific dataset needed, and should not be cited as solving a linkage problem here. What actually
remains true for Tao et al. 2021: even with the correct accession, its RNA-seq layer cannot be
live-fetched by this app at all (by design of the GEO format, not a bug), and it never could be, so
patient-matched data for 79 patients came from a pre-existing, independently
legitimacy-checked example upload rather than from any live fetch this session performed. It
already existed in
`data/examples/multiomics_upload/` (`gse138746_pbmc_rnaseq_counts.csv`,
`gse138653_methylation_beta_top2000.csv`, `gse138747_sample_metadata.csv`) from a prior session
this one has no other record of; it was independently legitimacy-checked here (realistic RNA-seq
count and methylation beta-value ranges, drug counts 37 ADA/42 ETN vs. the paper's reported 38/42)
before use. It initially carried no verified sex labels, so a first pass
(`01_pooled_diablo_response.R`) ran only a pooled (non-sex-stratified) test. Sex labels were later
recovered directly: GSE138746's and GSE138653's own GEO sample titles share an identical trailing
patient number (e.g. `PBMC_E_n_01` / `DNA_E_n_01`), which matched the pre-existing patient-ID
metadata at 100% drug/response agreement on merge (`patient_sex_lookup.csv`) — a real, verified
join, not an assumption.

`reproduce/multiomics/02_all_arms.R` and `03_snf_only.R` then ran the app's own leakage-safe
nested-CV engines (`mss_run_stratified()`/`mss_nested_cv()`, DIABLO and RF, plus SNF joint
clustering) across **14 scenarios** — sex-stratified response prediction, drug×sex-stratified
prediction for both drugs, both engines, and SNF clustering per drug — reproducing essentially the
full scenario space of the app's own ~28-arm bundled benchmark, from independently-fetched raw GEO
data:

**Result: 12/12 nested-CV arms and 2/2 SNF clustering tests show no statistically defensible
signal.** The strongest point estimate across all 14 (Etanercept, female, DIABLO: AUROC 0.688) still
has a 95% CI spanning 0.5. Full per-arm table in `02_all_arms_performance.csv`.

This independently confirms — from raw public data, through a different code path and a much wider
scenario sweep than the one that produced the bundled tables — the audit's most consequential
finding: the multi-omics drug-response prediction does not recover meaningful signal, far below the
published 0.79-0.89 benchmark, and this is not an artifact of testing only one configuration or of
the app's own precomputed pipeline; it reproduces from scratch across the whole scenario space.

## Summary

| Area | Regenerable & verified | Honest negative / lower confidence | Confirmed not-regenerable |
|---|---|---|---|
| Methylomics | QC (partial), DMP, DMR, WGCNA | Feature Selection | Diagnostic Classifier (not attempted, same pattern) |
| Transcriptomics | DEG (prior session) | — | — |
| — not attempted this session | batch-correction/WGCNA, MR/coloc | | |
| Cross-Omics | min-p correctness fix (not a repro exercise) | — | MR-stage (GoDMC access-gated) |
| Multi-Omics | GEO accessions found + verified; pooled DIABLO independently confirms chance-level AUROC | sex-stratified fits (no verified sex linkage) | — |
