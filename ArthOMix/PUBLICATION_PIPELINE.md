# Publication Pipeline Map

Maps each analysis stage to the scientific workflow it implements, per vertical. Being filled in as the publication-ready reorganization proceeds — see `CODE_MAP.md` for the file-level map and `REFACTORING_NOTES.md` for known gaps/inconsistencies not addressed by this reorganization.

## Cross-Omics

```
Transcriptomics DEG results + Methylomics DMP/DMR results
 ↓
Data Loading & Harmonization           (01_Data/)
 ↓
Expression-Methylation Integration     (02_Expression_Methylation_Integration/)
   — classifies each gene/CpG pair by concordance direction
 ↓
Biomarker Convergence                  (03_Biomarker_Convergence/)
   — ranks/tiers candidates by cross-omics evidence strength
 ↓
Cross-Omics MR                         (04_Cross_Omics_MR/)
   — two-sample Mendelian Randomization on convergent candidates
```

Every stage's output (`cross_results$<stage>`) is summarized for ArthOChat by `build_cx_context()` in `R/submodules_registry.R`, and downstream biomarker-card modules (Transcriptomics, Multiomics) call back into `cx_classify_evidence()`/`cx_harmonize_gene_ids()` from `R/crossomics/functions/integration/crossomics_integration_helpers.R`.

## Methylomics

```
Preloaded GSE42861 / GEO fetch / upload (matrix+sample sheet or IDAT)
 ↓
Quality Control (02) → Normalisation (03) → Cell Type Deconvolution (04)
 ↓
Differential Methylation Position (05) → Differential Methylation Region (06) → Sex Interaction Analysis (07)
 ↓
WGCNA Co-Methylation Network (08) → Candidate CpGs (09)
 ↓
ML Feature Selection (10) → Mendelian Randomization (11) → Colocalization (12)
 ↓
Diagnostic Classifier (13) → Validation (14)
 ↓
Biomarker Analysis (15)
```

Every stage's output (`results$<stage>`) is summarized for ArthOChat by `build_mx_context()` in `R/submodules_registry.R`. See `R/methylomics/README.md` for the per-stage dependency table, including the DMP-module-owned helpers several downstream stages call directly.

## Multiomics, Transcriptomics

Not yet mapped in this document — pending reorganization of those verticals (see `CODE_MAP.md`). The stage order for each is already authoritative in `R/submodules_registry.R`'s `MULTI_MODULES`/`TX_MODULES` lists.
