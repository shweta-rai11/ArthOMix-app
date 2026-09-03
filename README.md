<p align="center">
  <img src="docs/logo.png" alt="ArthOMix logo" width="220">
</p>


**ArthOMix** is a Shiny application for both single and multi--omics analysis — transcriptomics, methylomics, cross-omics integration, and multi-omics ArthOChatin one interactive tool. 


## Overview

ArthOMix lets a researcher take omics data — to either uploaded files, fetch GEO accession, or use the pre-loaded reference datasets — through a full analysis pipeline entirely inside the browser: sex-specific, sex-pooled differential expression, methylation, WGCNA co-expression modules, cross-omics concordance, and multi-omics integration (MOFA2, DIABLO, similarity network fusion), with an AI assistant ("ArthOChat") that can summarize results from any stage.


## Features

- **Transcriptomics** — differential expression, WGCNA, enrichment, nomogram, candidate gene selection
- **Methylomics** — QC, DMP/DMR analysis, WGCNA, candidate CpG selection, external validation
- **Cross-Omics** — expression–methylation integration, biomarker convergence, Mendelian Randomization
- **Multi-Omics** — MOFA2, DIABLO, and Similarity Network Fusion (SNF) integration, with a Biomarker Card summarizing cross-omics evidence
- **Flexible data intake** — upload your own files, pull a dataset by GEO accession, or use pre-loaded datasets, with each intake path isolated per pipeline.
- **ArthOChat** — an in-app AI assistant that summarizes live results from any of the four sub-modules.
- **Result tracking** — every analysis run records a provenance manifest for reproducibility

## How to use ArthOmix?

- **User guide** — It consist of all the user guide for each modules with screenshots.

- **Video** — You can also see Video tutorial for ArthOMix.


See [`ArthOMix/PUBLICATION_PIPELINE.md`](ArthOMix/PUBLICATION_PIPELINE.md) for the full stage-by-stage map of every vertical, and [`ArthOMix/CODE_MAP.md`](ArthOMix/CODE_MAP.md) for where each analysis lives in the codebase.

## Installation & running

ArthOMix is a Shiny app web app, not an installable R package. To reproduce it, clone the repo and run it in place with one of the two routes below. It pins its exact dependency set via [`renv.lock`](ArthOMix/renv.lock) (R 4.4.2) and ships a working [`Dockerfile`](ArthOMix/Dockerfile). Docker is the recommended way to run it, since the dependency stack is Bioconductor-heavy (limma, minfi, mixOmics, WGCNA, SNFtool, EpiDISH, and others).

### Before you start

- **Git LFS** — `data/` (~3.2GB of `.rds`/`.RData`/`.rda`) is tracked via Git LFS, not plain git. Run `git lfs install` *before* cloning, or `git lfs pull` afterwards if you already cloned without it — otherwise those files are just small pointer stubs and the app will fail to start.
- **Run from the `ArthOMix/` app directory** — the app resolves its data paths off the current working directory, not an absolute or installed location. It isn't relocatable: always `cd ArthOMix` (the inner app folder) before building/running/restoring.
- **`GITHUB_PAT`** — several pinned dependencies are GitHub-only (see `Remotes:` in [`ArthOMix/DESCRIPTION`](ArthOMix/DESCRIPTION)). `renv::restore()` re-resolves each one against the GitHub API, which quickly exceeds the unauthenticated 60 requests/hour limit — set a `GITHUB_PAT` environment variable first (any token with public read access works).

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

This binds to `http://127.0.0.1:7788` by default (set in `.Rprofile`) — different from Docker's `3838`.

### Environment variables

Copy `ArthOMix/.Renviron.example` to `ArthOMix/.Renviron` and fill in `SUPABASE_URL`/`SUPABASE_ANON_KEY` before starting the app. This isn't optional for a working demo: every page is gated behind a Supabase-authenticated login screen, and without real credentials the app still boots and *looks* fine but no one can get past login.

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

## Links

- **Repository** 

## References

- Shiny dashboard best practices.

## License

The code in this project is licensed under MIT license.
