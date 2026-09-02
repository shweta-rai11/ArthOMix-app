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

## Multiomics

```
Upload / GEO / preloaded fit → Active Multi-Omics Dataset (01_Data_Workspace, incl. optional live MOFA2)
 ↓
Cohort Harmonization (02)
 ↓
DIABLO & SNF Integration (03)
 ↓
SNF Clustering → patient subgroups (04)
 ↓
Biomarker Discovery, sex-stratified (05)
 ↓
Gene-CpG Concordance (06)
 ↓
Pathways (07)
 ↓
Biomarker Card — read-only interpretation (08)
 ↓
Results Summary & Reproducibility (09)
```

Every stage's output (`multi_results$<stage>`) is summarized for ArthOChat by `build_mo_context()` in `R/submodules_registry.R`. See `R/multiomics/README.md` for the per-stage dependency table, including how SNF Clustering delegates to Integration's own SNF fusion/clustering primitives.

## Transcriptomics

```
Upload / GEO / preloaded (01_Data)
 ↓
Overview and Datasets (02)
 ↓
Preprocessing and Batch Correction (03)
 ↓
Differential Expression (04) → WGCNA (05) → Candidate Gene Identification (06)
 ↓
Mendelian Randomization (07) → Colocalization (08)
 ↓
Feature Selection (09) → Diagnostic Model (10) → Sex Interaction Analysis (11)
 ↓
Cross-Tissue Validation (12) → Cross-Ancestral Validation (13)
 ↓
Functional Enrichment (14) → Immune Deconvolution (15)
 ↓
Nomogram (16)
 ↓
Biomarker Card (17)
```

Every stage's output (`results$<stage>`) is summarized for ArthOChat by `build_tx_context()` in `R/submodules_registry.R`. See `R/transcriptomics/README.md` for the per-stage dependency table.

---

All four verticals (Transcriptomics, Methylomics, Cross-Omics, Multiomics) are now reorganized. See `CODE_MAP.md` for the complete file-level index and `REFACTORING_NOTES.md` for known issues intentionally left unfixed.
