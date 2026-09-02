# functions/biomarker_convergence

Shared helpers used by the Biomarker Convergence stage (`../../03_Biomarker_Convergence/`) and reused by the Cross-Omics MR stage (`../../04_Cross_Omics_MR/`).

- **`crossomics_biomarkerconv_helpers.R`** — `cx_bc_*` functions: evidence-tier ranking/scoring of candidate genes from `cross_results$integration`, plus `cx_bc_load_precomputed()`/`cx_bc_relabel()` which `../../04_Cross_Omics_MR/mod_cross_mr_stage.R` calls directly to carry ranked candidates forward into MR.
