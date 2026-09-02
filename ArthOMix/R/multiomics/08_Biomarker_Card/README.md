# 08_Biomarker_Card

`MULTI_MODULES` stage `id = "biomarkercard"`, title **"Biomarker Card"**.

- **File**: `mod_multi_biomarkercard.R`.
- **Input**: `multi_results$concordance$df` exclusively — the one table already joining expression and methylation evidence.
- **Main analysis**: read-only integrated per-biomarker interpretation (does a candidate gene–CpG pair have Transcriptomics-only, Methylomics-only, or Multiomics-supported evidence?), using `cx_classify_evidence()` from `R/crossomics/functions/integration/crossomics_integration_helpers.R` — the exact same classifier Cross-Omics Integration uses, applied read-only for display.
- **Output**: none — never writes back to `multi_results`.
- **UI**: Multiomics → Sub-modules → Biomarker Card.
