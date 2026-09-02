# 06_Candidate_Gene_Identification

`TX_MODULES` stage `id = "candidates"`.

- **File**: `mod_candidates.R`.
- **Input**: `results$wgcna`, `results$dge`.
- **Main analysis**: intersects WGCNA module genes with DEGs, sex-stratified, with a Venn diagram.
- **Output**: `results$candidates`.
- **UI**: Transcriptomics → Sub-modules → Candidate Gene Identification.
