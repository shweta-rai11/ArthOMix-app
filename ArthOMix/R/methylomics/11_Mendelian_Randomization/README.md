# 11_Mendelian_Randomization

`MX_MODULES` stage `id = "mr"`.

- **File**: `mod_methyl_mr.R` (mostly self-contained mQTL/GWAS logic).
- **Input**: candidate CpGs + mQTL/GWAS summary statistics (upload or default instrument table).
- **Main analysis**: two-sample Mendelian Randomization on methylation instruments.
- **Output**: `results$mr`.
- **UI**: Methylomics → Sub-modules → Mendelian Randomization.
- **Dependencies**: `../functions/annotation.R`.
- **Note**: has no dedicated test file (`REFACTORING_NOTES.md`) — a pre-existing coverage gap, not touched by this reorganization.
