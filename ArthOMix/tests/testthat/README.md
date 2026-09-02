# Automated Test Index

90 files, all flat in this directory — `testthat` (edition 3, see `DESCRIPTION`) only
discovers `test-*.R` files directly under `tests/testthat/`, it does not recurse into
subfolders, so this directory can't be split into per-vertical folders without breaking
`devtools::test()` / `R CMD check`. Organization instead comes from a strict filename
prefix, mirroring the numbered-stage folders under `R/` (see `../../CODE_MAP.md`).

| Prefix | Vertical / area | Maps to |
|---|---|---|
| `test-app-*` | Shared app infrastructure | auth, data loaders/paths, full mounted-app server, smoke e2e |
| `test-crossomics-*` | Cross-Omics | `R/crossomics/` |
| `test-methyl-*` | Methylomics | `R/methylomics/` |
| `test-multi-*` | Multiomics | `R/multiomics/` |
| `test-txn-*` | Transcriptomics | `R/transcriptomics/` |
| `test-provenance-*` | Shared provenance helpers | `R/provenance.R` |

Within a prefix, the rest of the filename follows
`test-<prefix>-<module-or-stage>-<functions\|server\|e2e>.R`, e.g.
`test-methyl-dmp-functions.R` (pure helpers for the DMP stage) vs.
`test-methyl-dmp-server.R` (`testServer()` reactive tests for the same module).

## Running a subset

`test_dir()`'s `filter` matches against the filename with `test-` and `.R` stripped, so
filter on the prefix to run just one vertical:

```r
# whole suite
testthat::test_dir("tests/testthat")

# one vertical
testthat::test_dir("tests/testthat", filter = "^methyl")
testthat::test_dir("tests/testthat", filter = "^txn")

# one module across both its functions/server files
testthat::test_dir("tests/testthat", filter = "^methyl-dmp")

# a single file
testthat::test_file("tests/testthat/test-methyl-dmp-functions.R")
```

See `tests/MANUAL_TESTING.md` for the manual end-to-end checklist that covers gaps
`testServer()` can't reach (full file-upload -> multi-tab reactive chains).
