# ArthOMix

**ArthOMix** is a Shiny application for multi-omics analysis of Rheumatoid Arthritis data — transcriptomics, methylomics, cross-omics integration, and multi-omics (MOFA2/DIABLO/SNF) workflows in one interactive tool.

[![R Tests](https://github.com/shweta-rai11/ArthOMix-app/actions/workflows/r-tests.yml/badge.svg)](https://github.com/shweta-rai11/ArthOMix-app/actions/workflows/r-tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

ArthOMix lets a researcher take Rheumatoid Arthritis omics data — uploaded files, a GEO accession, or one of the app's bundled reference datasets — through a full analysis pipeline entirely inside the browser: differential expression / methylation, WGCNA co-expression modules, cross-omics concordance, and multi-omics integration (MOFA2, DIABLO, similarity network fusion), with an AI assistant ("ArthOChat") that can summarize results from any stage.

It is a hand-rolled classic Shiny app (`ui.R` / `server.R` / `global.R` + `R/` submodules), **not** an R package — see [`ArthOMix/DESCRIPTION`](ArthOMix/DESCRIPTION) for the full dependency surface and [`ArthOMix/CODE_MAP.md`](ArthOMix/CODE_MAP.md) for how the code is organized.

## Features

- **Transcriptomics** — differential expression, WGCNA, enrichment, nomogram, candidate gene selection
- **Methylomics** — QC, DMP/DMR analysis, WGCNA, candidate CpG selection, external validation
- **Cross-Omics** — expression–methylation integration, biomarker convergence, Mendelian Randomization
- **Multi-Omics** — MOFA2, DIABLO, and Similarity Network Fusion (SNF) integration, with a Biomarker Card summarizing cross-omics evidence
- **Flexible data intake** — upload your own files, pull a dataset by GEO accession, or use bundled reference datasets, with each intake path isolated per pipeline (see [`ArthOMix/REFACTORING_NOTES.md`](ArthOMix/REFACTORING_NOTES.md))
- **ArthOChat** — an in-app AI assistant that summarizes live results from any of the four verticals
- **Provenance tracking** — every analysis run records a provenance manifest for reproducibility

## Pipeline architecture

```
Transcriptomics DEG results + Methylomics DMP/DMR results
 ↓
Cross-Omics: Data Loading & Harmonization
 ↓
Cross-Omics: Expression–Methylation Integration  (concordance classification)
 ↓
Cross-Omics: Biomarker Convergence                (ranks candidates by cross-omics evidence)
 ↓
Cross-Omics: Mendelian Randomization               (two-sample MR on convergent candidates)
```

See [`ArthOMix/PUBLICATION_PIPELINE.md`](ArthOMix/PUBLICATION_PIPELINE.md) for the full stage-by-stage map of every vertical, and [`ArthOMix/CODE_MAP.md`](ArthOMix/CODE_MAP.md) for where each analysis lives in the codebase.

## Installation & running

ArthOMix pins its exact dependency set via [`renv.lock`](ArthOMix/renv.lock) (R 4.4.2) and ships a working [`Dockerfile`](ArthOMix/Dockerfile). Docker is the recommended way to run it, since the dependency stack is Bioconductor-heavy (limma, minfi, mixOmics, WGCNA, SNFtool, EpiDISH, and more).

### Docker (recommended)

```sh
cd ArthOMix
docker build -t arthomix .
docker run --rm -p 3838:3838 arthomix
```

Then open `http://localhost:3838` in a browser.

### Local R / renv

```r
# from the ArthOMix/ directory
renv::restore()
shiny::runApp()
```

### Environment variables

Copy `ArthOMix/.Renviron.example` to `ArthOMix/.Renviron` and fill in the required credentials (e.g. Supabase auth) before starting the app.

## Sample data

The app supports three ways to bring in data for each vertical (Transcriptomics, Methylomics, Multi-Omics):

- **Upload** your own expression/methylation matrices
- **GEO accession** — pull a public dataset directly by GEO ID
- **Preloaded reference datasets** — bundled RA datasets for exploring the app without external data

## Testing

```sh
cd ArthOMix
Rscript -e 'testthat::test_dir("tests/testthat")'
```

Tests run automatically on push via [`.github/workflows/r-tests.yml`](.github/workflows/r-tests.yml).

## Contributing

Issues and pull requests are welcome. Please open an issue describing the change before submitting a large PR.

## License

MIT — see [LICENSE](LICENSE).
