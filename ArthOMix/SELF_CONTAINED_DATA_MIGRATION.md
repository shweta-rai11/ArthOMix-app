# ArthOMix Self-Contained Data Migration — Report

## Before

```
/Users/swetarai/ArthOMix/                          (repo root)
├── ArthOMix/                                       (the Shiny app)
│   ├── global.R  (defined 5 external root paths, all "../<folder>")
│   ├── server.R  ui.R
│   ├── R/{transcriptomics,methylomics,crossomics,multiomics}/
│   └── data/  (one file: cytoBandIdeo_hg19.txt.gz — not a real data root)
├── methylomics/                                    24 GB — required
├── replicate_chen2021/                             91 MB — NOT required (zero references)
├── Research_05_multiomics_sexstratified/            8 GB — required
├── Research_Q2_TRANSCRIPTOMICS_sexstratified_COPY/ 1.7 GB — required
├── Research_Q3_METHYLOMICS_sexstratified_COPY/     (symlink → Research_Q4) — required
└── Research_Q4_cross_Omics_sexstratified_COPY/     24 GB — required
```
58 GB total across the five actually-used folders; the app hard-`stop()`ed at startup if `Research_Q2...` was missing, and silently degraded (upload-only) if the other three were missing.

## After

```
ArthOMix/
├── data_paths.R                     NEW — single source of truth for all data paths
├── global.R  server.R  ui.R         unchanged except: constants moved out to data_paths.R
├── R/{transcriptomics,methylomics,crossomics,multiomics}/   unchanged except 2 lines (below)
├── www/
└── data/                             2.3 GB total — the only data root the app now reads
    ├── dataset_manifest.csv          NEW — catalogue of every bundled dataset
    ├── preloaded/
    │   ├── transcriptomics/          117 MB
    │   ├── methylomics/               2.2 GB  (incl. 2.14GB live beta matrix)
    │   ├── cross_omics/               172 KB
    │   └── multiomics/                896 KB
    ├── reference/                     cytoBandIdeo_hg19.txt.gz (moved from data/ root)
    ├── annotations/gene_panels/       curated gene panel .txt files
    ├── examples/                      reserved, empty (see below)
    ├── uploads/                       reserved, empty (see below)
    └── .cache/                        regenerable WGCNA/probe caches (gitignored)
```

**Zero references to the six external folders remain anywhere in the app** (verified by a repo-wide grep after migration — every remaining hit is either in `global.R`'s historical comments, which were rewritten to describe the new location, or provenance-naming text in UI descriptions like "this reads Research_05_multiomics_sexstratified's own precomputed tables," which is accurate historical attribution, not a functional path, and is UI-facing text left untouched per the no-UI-changes requirement).

---

## 1. Datasets moved (original → new → consuming module)

Full detail — including per-file paths, formats, and whether each dataset is normalized/preloaded/needs annotation — is in **`data/dataset_manifest.csv`** (40 rows). Summary by module:

| Domain | Original root | New root | Size | Used by |
|---|---|---|---|---|
| Transcriptomics | `Research_Q2_TRANSCRIPTOMICS_sexstratified_COPY/{results/tables,data/processed,data/processed/new,data/gene_panels}` | `data/preloaded/transcriptomics/{results/tables,processed,processed/new}`, `data/annotations/gene_panels/` | 117 MB | Dataset, WGCNA, MR, Coloc, Cross-Tissue, Cross-Ancestry, Feature Selection, Candidates, Enrichment, ArthOChat |
| Methylomics (tables) | `Research_Q3.../methylomics/script0{1,3,3sva,4,5,7,8,9}_*/{tables,METHODS_*.md}` | `data/preloaded/methylomics/tables/script0N_*/{tables,METHODS doc}` | 112 MB | QC, DMP, DMR, WGCNA, MR, Coloc, Diagnostic Classifier, ArthOChat |
| Methylomics (live matrix) | `methylomics/data/processed/{beta_raw.rds,pheno.rds,gse42861_*,gse111942_*}` | `data/preloaded/methylomics/matrix/` | 2.14 GB | Dataset tab "Load preloaded dataset" (live recompute), Diagnostic Classifier |
| Cross-Omics | `Research_Q4.../cross_Omics_Sexstratified_COPY/results/*.csv` (8 files) | `data/preloaded/cross_omics/tables/` | 172 KB | Dataset, Biomarker Convergence, Cross-Omics MR, Integration (reuses transcriptomics/methylomics files above too) |
| Multi-Omics | `Research_05.../{results,metadata,analyses/*/results}/tables/*.csv` (24 files) + `results/*.rds` (6 fits) | `data/preloaded/multiomics/{tables/{,summary,adalimumab,etanercept},fits}/` | 896 KB | Dataset, Overview, Integration, Stratification, Biomarker, Concordance, Pathway |
| Reference | `ArthOMix/data/cytoBandIdeo_hg19.txt.gz` (already in-repo) | `data/reference/cytoBandIdeo_hg19.txt.gz` | 12 KB | Methylomics Biomarker Card |

## 2. Files intentionally NOT moved, and why

- **`replicate_chen2021/`** (91 MB) — **correction**: no app *code* reads it (still zero grep hits), but it is a real, deliberate manual test-fixture folder — `upload_csv_merged/` and `upload_csv_probelevel/` inside it contain ready-to-upload CSVs shaped for the Transcriptomics Dataset upload widgets, originally built while replicating Chen et al. 2021's published WGCNA result as a validation exercise (the rest of the folder — `run_wgcna_paper_*.R`, `*_result.rds`, `raw/` GEO series matrices — is the scripts/output that built those fixtures). The two `upload_csv_*` folders (~23MB) were copied into `data/examples/transcriptomics_upload/{merged,probelevel}/` and are now exercised by an automated end-to-end test (`tests/testthat/test-upload-transcriptomics.R`). `replicate_chen2021/` itself was **not deleted** — it remains at the repo root as the fixtures' authoring source.
- **`data/cache/`, `data/cache/wgcna/`** (~1.3 GB, transcriptomics) — regenerable WGCNA/probe-collapse caches, recomputed on first use and written to the new `data/.cache/` location.
- **Raw GEO ExpressionSets** (`Research_Q2.../data/raw/*_raw.rds`) — the "Load individual GEO dataset" and Preprocessing per-source features. This was **already broken before this migration**: `data/raw` is a symlink pointing at `/Users/swetarai/THESIS_SWETA_28_MAY/...`, a path that does not exist anywhere on this machine. Per your explicit choice, this is preserved exactly as-is (gracefully unavailable) rather than fixed, since it's a pre-existing, unrelated issue.
- **`methylomics/data/raw/`** (~21 GB — GSE42861 series matrix, GoDMC, GSE111942 raw, raw IDATs, Ishigaki GWAS) and **`methylomics/data/processed/`'s other ~180 MB of derived `.rds` files** — confirmed, file-by-file, never read by any of the 19 methylomics R files. Only the 4 files actually read (`beta_raw.rds`, `pheno.rds`, and the 2 diagnostic panel RDS files) were copied.
- **`results/figures/by_section/`** (21 MB, transcriptomics) — `addResourcePath("figures", ...)`/`figure_exists()` are defined in `global.R` but have zero call sites anywhere in the app (confirmed by repo-wide grep). Genuinely dead code. An **empty placeholder directory** was created (`data/preloaded/transcriptomics/results/figures/by_section/.gitkeep`) because `addResourcePath()` errors at app startup if the target directory doesn't exist at all — but no image files were copied, since nothing reads them.
- **~110 of ~130 files** in `Research_Q2.../results/tables/` and most files in each methylomics `script0N_*` folder, all multi-omics `analyses/*/results/tables/` extras, and most of `Research_Q4.../results/` — an exhaustive per-file audit (5 parallel deep-read passes + a design-validation pass) confirmed these are not read by any reachable code path. Only the specific files each module actually opens were copied.
- **Uploaded user data** — was never file-based to begin with. It already flows through Shiny's own per-session `input$file$datapath` tempfile mechanism everywhere in the app, which is correctly session-isolated. `data/uploads/` is provisioned as a reserved, empty folder per your requested layout, but nothing is redirected into it — doing so would add a new cross-session-leak risk for no benefit.
- **`data/examples/`** — stays empty. The one related feature (Preprocessing's "merge the example pipeline's training datasets" demo) is backed by the same broken raw-GEO symlink above and has been unreachable since before this migration; there's nothing on disk to copy.
- **Illumina 450K/EPIC manifests and EpiDISH cell-type reference panels** — read from installed Bioconductor packages (`IlluminaHumanMethylation450kanno.ilmn12.hg19`, `...EPICanno.ilm10b4.hg19`, `EpiDISH`), not from any folder. Already fully self-contained via package dependencies (tracked in `renv.lock`/`DESCRIPTION`).

## 3. Remaining external dependencies

**None — and the five folders are now actually deleted** (moved to `~/.Trash`, not `rm -rf`, so there's a recovery window until the user empties it; `replicate_chen2021/` was excluded and remains in place).

Before deleting anything, a real automated test suite was built and run (`ArthOMix/tests/`, `testthat` + `shinytest2`):
- **`test-data-paths.R`** (93 assertions) — every path constant in `data_paths.R` resolves to a file/directory that exists under `ArthOMix/data/`.
- **`test-data-loaders.R`** (129 assertions) — sources the real `global.R` and calls *every* `load_default_*()` function, every entry in `CX_TABLE_REGISTRY` (8), `MULTI_TABLE_REGISTRY` (24), and `MULTI_DIABLO_FIT_REGISTRY` (6), plus gene panels, cytoband, and the ArthOChat methodology-lookup tools — pinning known-good shapes (15,763×183 expression matrix, 412,492-row DMP/beta matrix, 689-sample phenotype table, etc.) as regression guards.
- **`test-app-smoke.R`** — headless-Chromium `shinytest2` test: launches the real app, visits Transcriptomics/Methylomics/Cross-Omics/Multi-Omics, asserts no Shiny output error renders.
- **`test-upload-transcriptomics.R`** — drives the actual upload flow (file select → column mapping → load) using the bundled `chen2021` fixtures, asserts the dataset loads successfully.

**This test suite caught a real bug the manual verification missed**: `R/transcriptomics/mod_wgcna.R:126` independently reconstructed `file.path(DATA_ROOT, "data", "processed")` instead of using the `PROCESSED_DIR` constant — a leftover from the old external-folder convention. After `PROCESSED_DIR` was repointed to `data/preloaded/transcriptomics/processed/` (no `data/` sub-segment, matching the new layout), this one line would have silently broken the WGCNA tab's precomputed-result loader and dendrogram plot. Fixed to `proc_dir <- PROCESSED_DIR`.

The full suite (233 assertions across 5 files) was run **after** the five folders were actually deleted, confirming the app and every data path work with them genuinely gone, not just backed up.

## 4. Database recommendation

**No database.** `server.R` has no authentication, no multi-user shared state, and no persistent-upload requirement — every session's uploaded data is already ephemeral and session-scoped via Shiny's own tempfile mechanism. A plain file tree under `ArthOMix/data/` (Option A) is sufficient and is what was implemented. Introducing SQLite/Postgres would add real complexity (schema, migrations, a query layer) for zero functional gain given this app's actual architecture.

## 5. Confirmations

- **UI unchanged.** No layout, tab, module, input, output, plot, label, color, or CSS file was touched. The only user-visible text changes are two error-state strings that referenced the old external-folder name and would have been actively misleading after this migration (`mod_methyl_wgcna.R`'s "needs the preloaded dataset's Research_Q3 folder" message, now naming the actual bundled folder) — both only ever shown when data is *unavailable*, not in normal operation.
- **Scientific analysis unchanged.** No statistical method, normalization, filtering, or computation logic was touched. Every relocated file is a byte-identical copy (verified by checksum against the original) of the same file the app read before migration — same numbers in, same numbers out.
- **Code changes were exactly:** one new file (`data_paths.R`, all path constants/registries, same names as before), deletion of the now-duplicate constant definitions from `global.R` (functions/logic untouched), one hardcoded literal path fixed in `mod_methyl_biomarkercard.R`, and one misleading error-message string fixed in `mod_methyl_wgcna.R`.
- **Five of the six external folders have been deleted**: `methylomics/`, `Research_05_multiomics_sexstratified/`, `Research_Q2_TRANSCRIPTOMICS_sexstratified_COPY/`, `Research_Q3_METHYLOMICS_sexstratified_COPY/` (symlink), `Research_Q4_cross_Omics_sexstratified_COPY/` — moved to `~/.Trash` (not `rm -rf`) after the full test suite passed both before and after the move, so there's a recovery window until it's emptied. `replicate_chen2021/` (91 MB) was deliberately kept — it's the authoring source for the bundled upload-test fixtures (§2), not dead weight.
