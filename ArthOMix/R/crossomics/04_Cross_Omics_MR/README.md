# 04_Cross_Omics_MR

`CX_MODULES` stage `id = "mrstage"`, titled **"Cross-Omics MR"** (group "Genetics") in the Cross-Omics Sub-modules grid.

- **Files**: `mod_cross_mr_stage.R` (Shiny module trio), `crossomics_mrstage_helpers.R` (`cx_mr_*` — kept alongside the module here rather than under `functions/`, since it's used only by this one stage).
- **Input**: `cross_results$biomarkerconv` from `../03_Biomarker_Convergence/` (via `cx_bc_load_precomputed`/`cx_bc_relabel`, reused from `../functions/biomarker_convergence/`).
- **Main analysis**: two-sample Mendelian Randomization on the convergent gene/CpG candidates carried forward from Biomarker Convergence.
- **Output**: `cross_results$mrstage` — MR estimates (effect size, SE, p-value) per candidate.
- **UI**: Cross-Omics → Sub-modules → "Cross-Omics MR" card. `server.R` passes this stage's server call an extra `app_session = session` argument (for ArthOChat), unlike the other two `CX_MODULES` stages.
- **Dependencies**: `crossomics_mrstage_helpers.R` (same folder); `cx_bc_load_precomputed`/`cx_bc_relabel` from `../functions/biomarker_convergence/crossomics_biomarkerconv_helpers.R`; `cx_empty_state` from `../functions/integration/crossomics_integration_plots.R`.
