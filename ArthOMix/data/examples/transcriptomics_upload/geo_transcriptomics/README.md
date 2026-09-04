# GEO transcriptomics test datasets (cancer + rare immune disease, group/treatment + sex)

Built 2026-09-03 from NCBI GEO for the **Transcriptomics** tab ("Upload your own data"). Every dataset has a
disease-vs-control and/or treated-vs-untreated grouping **and** both sexes recorded on GEO, so the same files
exercise Differential Expression (group contrast) and the Sex Interaction module (male vs female). Same layout as
`../../multiomics_upload/geo_multiomics/`:

- `<GSE>_exp.csv` - genes in rows (`gene` column = HGNC symbol), samples in columns (GSM ids). Array probes were
  collapsed to genes by the same MaxMean rule as the app's `collapse_probes_to_genes()`; RNA-seq files are the
  deposited raw counts (Ensembl ids mapped to symbols with `org.Hs.eg.db` where needed).
- `<GSE>_sample.csv` - one row per sample: `sample` (= GSM, the column header in the matrix), `gsm`, `title`,
  `group` (cleaned label used for the contrast), `sex` (F/M/NA), `group_geo` (the label exactly as on GEO), then
  **every** GEO sample characteristic (snake_case; a characteristic that collides with a base column gets `_geo`),
  `platform`, `dataset`.

The Upload tab auto-picks `sample`, `group` and `sex`; any other column (e.g. `glucocorticoids`, `treatment`,
`timepoint`, `chemotherapy_adjuvant_type`) can be chosen as the group or batch column instead. Orientation is
Feature x Sample ("Feature and Sample (first column = feature ID)").

| GSE | Disease / design | `group` values (n) | sex F/M (NA) | Platform, scale -> **Data type to pick** | Genes x samples, size | Other useful columns |
|---|---|---|---|---|---|---|
| GSE110169 | RA, SLE, healthy whole blood (PAXgene); steroid-treated vs not | Normal 77 / RA 84 / SLE 82 | 188/50 (5) | Affy HG-U219 (GPL13667), RMA log2 1.3-13.3 -> **Already log-transformed** | 19,041 x 243, 28 MB | `glucocorticoids` TRUE 72 / FALSE 89 / NA (healthy) 82, `batch` (4) |
| GSE8650 | Systemic-onset JIA, SLE, PAPA (rare autoinflammatory), healthy; PBMC | JIA Systemic Onset 58 / SLE 38 / PAPA 6 / Healthy 21 | 91/32 | Affy U133A (GPL96), MAS5 **converted to log2(x+1)** 0.1-18 -> **Already log-transformed** | 12,548 x 123, 13 MB | `treatment` (None 30 / Steroids PO / Methotrexate / Anakinra / ...), `age`, `symptoms` |
| GSE17114 | Behcet's disease PBMC vs healthy | BD 15 / Control 14 | 15/14 | Affy U133 Plus 2 (GPL570), log2 1.9-14.6 -> **Already log-transformed** | 21,655 x 29, 5 MB | `immunosuppressors` yes 9 / no 20, `geographical_origin` |
| GSE19314 | Sarcoidosis vs healthy whole blood | Sarcoidosis 38 / Healthy 20 / Hypersensitivity Pnemonitis 6 (+2 single-sample labels) | 41/25 | Affy U133 Plus 2 (GPL570), log2 1.1-14.2 -> **Already log-transformed** | 21,655 x 66, 12 MB | `dx` (sarc_low_lung_fx / sarc_high_lung_fx / health / hp), `array_batch` |
| GSE83456 | Sarcoidosis, pulmonary + extra-pulmonary TB, healthy; whole blood | HC 61 / Sarcoid 49 / PTB 45 / EPTB 47 | 85/117 | Illumina HT-12 v4 (GPL10558), **normalised log-ratios, negative values** (-8.8 to 7.0) -> **Already log-transformed** | 31,334 x 202, 53 MB | `ethnicity`, `age` |
| GSE39582 | Colon cancer, adjuvant chemotherapy vs none | AdjuvantChemo 240 / NoChemo 326 / NA 19 | 263/322 | Affy U133 Plus 2 (GPL570), RMA log2 0.8-16.1 -> **Already log-transformed** | 21,655 x 585, 101 MB | `chemotherapy_adjuvant_type` (5FU/FUFOL/FOLFOX/FOLFIRI), `tnm_stage`, `mmr_status`, `cimp_status`, `kras_mutation`, `braf_mutation`, `os_event`, `dataset_geo` (discovery/validation) |
| GSE68465 | Lung adenocarcinoma vs normal lung; adjuvant chemo / RT flags | Adenocarcinoma 443 / Normal 19 | 220/223 (19, all normals) | Affy U133A (GPL96), MAS5 **converted to log2(x+1)** 0.1-16.9 -> **Already log-transformed** | 12,548 x 462, 47 MB | `clinical_treatment_adjuvant_chemo` Yes 89 / No 341, `clinical_treatment_adjuvant_rt`, `disease_stage`, `race`, `vital_status` |
| GSE10846 | Diffuse large B-cell lymphoma, R-CHOP vs CHOP | R-CHOP 233 / CHOP 181 / NA 6 | 172/224 (24) | Affy U133 Plus 2 (GPL570), log2 0-17.9 -> **Already log-transformed** | 21,655 x 420, 56 MB | `final_microarray_diagnosis` (GCB/ABC), `stage`, `follow_up_status`, `ecog_performance_status` |
| GSE112057 | JIA subtypes + paediatric IBD vs control, whole blood RNA-seq | polyJIA 46 / oligoJIA 43 / sJIA 26 / CD 60 / UC 15 / Control 12 | 110/80 (12) | RNA-seq (GPL11154), **raw counts** -> **Raw counts** | 26,305 x 202, 16 MB | `age_at_diagnosis`, `ancestry` |
| GSE178388 | MIS-C, paediatric COVID-19, healthy; whole blood RNA-seq, before/after treatment | COVID 23 / MISC 8 / Healthy 8 | 19/20 | RNA-seq (GPL24676), **raw counts** -> **Raw counts** | 35,949 x 39, 4 MB | `timepoint` (T0 = before treatment, T1-T5 after), `immnosuppressive`, `comorbidities`, `subject_id` |

Notes:
- GSE110169 is the one file that covers **all three** designs the app tests: disease vs healthy (`group`),
  glucocorticoid-treated vs untreated (`glucocorticoids`) and male vs female (`sex`). The 5 samples with sex NA
  and the 82 healthy with glucocorticoids NA are dropped automatically when those columns are used.
- GSE39582 group is the treatment flag; to contrast tumour subgroups instead, pick `mmr_status` or `cimp_status`
  as the group column. The 19 `NA` group samples have no chemotherapy record on GEO.
- GSE68465 has no sex on its 19 normal-lung samples, so a sex analysis restricted to the Normal group is not possible.
- GSE8650 and GSE68465 were deposited as MAS5 signal (unlogged, up to ~2.6e5); they were converted to
  log2(x + 1) here so they pass as "Already log-transformed" like the other arrays. Everything else is exactly as deposited.
- GSE83456 values are centred (negative allowed); use "Already log-transformed", not "Normalized".
- GSE19314 keeps two single-sample GEO labels ("HCV/sarcoid", "sarcoid vs. BO") in `group`; ignore them or use `dx`.
- GSE112057 sex is NA on 12 samples; the Control group (n = 12) contains both sexes.
- GSE10846 characteristics on GEO are prefixed "Clinical info:"; that prefix was stripped before making columns.

## Using "Retrieve from GEO" in the app instead of these files

- All eight array series above load directly from their accession (their GEO series matrix has values). The
  Affymetrix ones (GPL570, GPL96, GPL13667) collapse to gene symbols on the fly. **GSE83456 (Illumina GPL10558)
  loads at probe level** on a live fetch because the app's `collapse_probes_to_genes()` regex does not match
  Illumina's `Symbol` column - use this CSV instead, or the Preprocessing collapse step afterwards.
- GSE8650 spans two platforms (GPL96/GPL97); the app shows a platform picker - choose GPL96.
- GSE112057 and GSE178388 (RNA-seq) have **empty** series matrices; the app reports this and asks for an upload.
  These CSVs were built from `GSE112057_RawCounts_dataset.txt.gz` and `GSE178388_PedBatch_all_counts_matrix.csv.gz`.

## Rebuilding / adding a dataset

`build_geo_transcriptomics.R` in this folder rebuilds any of these (one dataset per call, cached downloads):

```
Rscript data/examples/transcriptomics_upload/geo_transcriptomics/build_geo_transcriptomics.R GSE110169 <out_dir> <cache_dir>
```

To add another series, append a `specs` entry (group characteristic, sex characteristic, optional recode / log2 /
supplementary file). For a plain GEO array series the shorter `data/examples/geo_to_arthomix_upload.R` also works.

## Verification (2026-09-03)

Each `_exp.csv` was read back with `data.table::fread` (as the Upload tab does) and passed through the app's
`tx_validate_expr_upload()` with the data type listed above: all 10 return ok, 0 duplicated genes, 0 NA cells,
every matrix column present in `_sample.csv`, and both sexes present inside every main group (exceptions noted above).
Group and sex counts come from GEO sample characteristics, not from the publications.
