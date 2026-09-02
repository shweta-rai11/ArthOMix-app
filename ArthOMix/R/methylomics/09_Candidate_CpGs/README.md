# 09_Candidate_CpGs

`MX_MODULES` stage `id = "candidates"`.

- **File**: `mod_methyl_candidates.R` (self-contained — no `functions/` dependencies).
- **Input**: `results$dmr`, `results$wgcna`.
- **Main analysis**: intersects DMR hits with WGCNA co-methylation modules to shortlist candidate CpGs.
- **Output**: `results$candidates`.
- **UI**: Methylomics → Sub-modules → Candidate CpGs.
