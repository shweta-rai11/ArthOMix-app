# functions

Shared helpers used across multiple Methylomics stages (kept flat, not split into `statistics/`/`plotting/`/etc. subfolders — the real code doesn't cleanly separate along those lines; `qc.R` alone mixes statistics, plotting, and report generation).

- **`qc.R`** — general-purpose statistics/plotting library for the vertical (filtering, batch correction, PCA/MDS, beta↔M conversion, sex-check, outlier detection, plots, report/HTML/zip generation). Used by 9 of 14 stages.
- **`annotation.R`** — array annotation lookup (`methyl_get_annotation`, `methyl_probe_is_cpg`). Used by 8 stages, plus called cross-vertical from `R/crossomics/functions/integration/crossomics_integration_helpers.R`.
- **`parse_upload.R`** — matrix/sample-sheet/probe-list/IDAT upload parsing and validation. Used by nearly every stage offering custom upload, plus mirrored by `R/transcriptomics/expression_type.R` and `R/crossomics/01_Data/crossomics_integration_upload.R`.
- **`normalization.R`** — normalization method engine (noob/funnorm/SWAN/dasen/BMIQ/PBC/quantile), diagnostics, status detection.
- **`idat_metrics.R`** — IDAT-derived QC metrics (bisulfite conversion, median intensity).

See `../README.md` for the full per-stage dependency table.
