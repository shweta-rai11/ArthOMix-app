# functions

Shared helpers used across multiple Transcriptomics stages.

- **`expression_type.R`** — expression-matrix scale/type heuristics (`looks_like_raw_counts()`, `looks_like_normalized_totals()`) and the declare-then-verify upload validator. Promoted out of `04_Differential_Expression/mod_dge.R` (which used to keep these as local closures) so `01_Data/mod_dataset.R`'s upload path, `mod_dge.R`'s decoupled upload path, and `15_Immune_Deconvolution/mod_deconvolution.R`'s run gate all share one implementation instead of drifting duplicates. Mirrors `R/methylomics/functions/parse_upload.R`'s validation pattern.

See `../README.md` for the full per-stage dependency table.
