# Methylomics precomputed-table reproducibility check (2026-09-08)

Standalone scripts that regenerate the Methylomics vertical's precomputed
tables (Section 2.5 of the thesis) end-to-end from the bundled raw matrix,
using the live application's own functions (loaded via `shiny::loadSupport()`,
the same bootstrap `runApp()` uses), and compare the result against the
tables shipped in `data/preloaded/methylomics/tables/`. Written in response
to the full-application audit finding that ~90% of the app's precomputed
tables had never been shown to reproduce from raw data.

Run from the `ArthOMix/` app directory, in order (each stage after the first
consumes the previous stage's own regenerated output, not the precomputed
one, so the chain measures true end-to-end divergence rather than being
reset to ground truth at every step):

```
Rscript reproduce/methylomics/01_qc.R
Rscript reproduce/methylomics/02_dmp.R
Rscript reproduce/methylomics/03_dmr.R
Rscript reproduce/methylomics/04_wgcna.R
```

Output (comparison CSVs, run logs, summary RDS files) is written to
`reproduce/methylomics/output/`.

## Results

**QC (689 samples).** The chromosome-Y sex-mismatch check is **not
reproducible from what is bundled**: `beta_raw.rds` already has zero
chrX/chrY probes (excluded upstream of what ships), so the sex-check columns
in the precomputed table cannot have come from this file. The PCA outlier
flag does not reproduce numerically either: regenerated PC1/PC2 are
essentially uncorrelated with the precomputed values (r = 0.05, -0.03), and
the two runs agree on zero of the samples they flag as outliers (3
precomputed, 7 regenerated, 0 in common).

**DMP (412,492 CpGs, both sexes).** The model reproduces closely at the
statistic level: t-statistic r = 0.955 (F) / 0.976 (M), Δβ r = 1.000 /
0.9997. The small FDR<0.05 counts are unstable at that correlation (F: 18
precomputed vs 11 regenerated, 11 overlap; M: 0 vs 0) — expected for a
faithfully reproduced model near a weak-signal significance boundary, not a
reproducibility failure.

**DMR.** Feeding the *regenerated* DMP statistics into DMRcate with the
app's own parameters recovers 64.2% of precomputed female DMRs and 62.9% of
precomputed male DMRs by genomic overlap (57.1% / 57.2% of the FDR<0.05
ones) — partial, not exact, reproduction.

**WGCNA.** Regeneration script written and (as of this writing) running;
results not yet recorded here. Compares module assignment against the
precomputed table via Adjusted Rand Index (module labels/colors are
arbitrary across runs, so ARI — not label matching — is the correct
agreement measure).

These findings are written into the thesis (`_MConverter.eu_ArthOMix.md_file.md`)
inline against each affected claim in Section 2.5, in a dedicated
"Reproducibility check" block before the Abstract, and as a Section 4.3
limitation.
