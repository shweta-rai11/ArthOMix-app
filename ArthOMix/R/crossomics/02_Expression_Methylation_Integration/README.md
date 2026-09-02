# 02_Expression_Methylation_Integration

`CX_MODULES` stage `id = "integration"`, titled **"Expression and Methylation"** (group "Data") in the Cross-Omics Sub-modules grid.

- **File**: `mod_cross_integration.R` (Shiny module trio).
- **Input**: the shared `cross_dataset` reactiveValues from `../01_Data/`, plus Transcriptomics `results$dge` and Methylomics `methyl_results$dmp`/`dmr` (passed in directly from `server.R` — this is the one stage in the vertical wired with extra arguments beyond the standard `CX_MODULES` loop signature).
- **Main analysis**: classifies each matched gene/CpG pair by direction-of-change concordance (e.g. Hyper-methylation + Down-regulation, Hypo + Up, and the two non-canonical combinations), using `cx_classify_evidence()`/`cx_classify_*()` from `../functions/integration/crossomics_integration_helpers.R`.
- **Output**: `cross_results$integration` — classified gene/CpG pairs plus summary counts per concordance category (`cx$summary$counts`), consumed by `build_cx_context()`'s bespoke integration formatter in `R/submodules_registry.R` and by Biomarker Convergence.
- **UI**: Cross-Omics → Sub-modules → "Expression and Methylation" card.
- **Dependencies**: `../functions/integration/crossomics_integration_helpers.R`, `../functions/integration/crossomics_integration_plots.R`.
