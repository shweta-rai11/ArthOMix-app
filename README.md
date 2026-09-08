<p align="center">
  <img src="docs/logo.png" alt="ArthOMix logo" width="220">
</p>


**ArthOMix** is a Shiny application for both single and multi-omics analysis - transcriptomics, methylomics, cross-omics integration, multi-omics and ArthOChat in one interactive tool. 

The app is not currently deployed behind a permanent, authenticated public URL. Run it locally - see [Installation & running](#installation--running) below.

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

- **Video** — You can also see Video tutorial for ArthOMix.


## Installation & running

ArthOMix is a Shiny web app, not an installable R package. To reproduce it, clone the repo and run it in place with one of the two routes below. It pins its exact dependency set via [`ArthOMix/renv.lock`](ArthOMix/renv.lock) (R 4.4.2) and ships a [`Dockerfile`](ArthOMix/Dockerfile); the dependency stack is Bioconductor-heavy (limma, minfi, mixOmics, WGCNA, SNFtool, MOFA2, and others).

### Before you start

- **Git LFS** — large precomputed files under `ArthOMix/data/` are tracked via Git LFS, not plain git. Run `git lfs install` *before* cloning, or `git lfs pull` afterwards if you already cloned without it. The app now fails at startup with a clear message naming the affected file(s) if any are still LFS pointer stubs, instead of failing deep inside a random module.
- **Run from the `ArthOMix/` app directory** — the app resolves its data paths off the current working directory, not an absolute or installed location. Always `cd ArthOMix` (the inner app folder) before building/running/restoring.
- **`GITHUB_PAT`** — several pinned dependencies are GitHub-only (see `Remotes:` in [`ArthOMix/DESCRIPTION`](ArthOMix/DESCRIPTION)). `renv::restore()` re-resolves each one against the GitHub API, which quickly exceeds the unauthenticated 60 requests/hour limit — set a `GITHUB_PAT` environment variable first (any token with public read access works).
- **No authentication** — the app does not currently gate access behind a login. Do not expose a running instance on a public, unauthenticated network address. For remote access, use an SSH port-forward rather than a public tunnel — see [`deploy/README.md`](deploy/README.md).

### Docker

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

This binds to `http://127.0.0.1:7788` by default (set in `ArthOMix/.Rprofile`) — different from Docker's `3838`.

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

- **Repository** — [github.com/shweta-rai11/ArthOMix-app](https://github.com/shweta-rai11/ArthOMix-app)

## References

- Shiny dashboard best practices.

## License

The code in this project is licensed under MIT license.
