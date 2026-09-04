# GEO methylomics test datasets (treated vs control + sex)

Built 2026-09-03 from NCBI GEO for testing the **Methylomics** tab. Every series is a human Illumina 450K or EPIC
methylation array from cancer or a rare/immune-mediated disease, carries a treated-vs-control (or pre/post-treatment)
grouping **and** a sex column, and was verified against GEO's per-sample characteristics and the actual FTP files.

Files per series (same layout as `../../multiomics_upload/geo_multiomics`):

- `<GSE>_meth.csv` - CpGs in rows (`cpg` column), samples in columns (GSM ids). **Top 20,000 most-variable CpGs** only.
- `<GSE>_sample.csv` - one row per sample: `sample` (= matrix column header), `title`, `group`, `sex` (F/M), every GEO
  characteristic as an extra snake_case column, `platform`, `dataset`.
- `<GSE>.xlsx` - the overview row for that dataset + its sample sheet. `GEO_methylomics_verification.xlsx` /
  `GEO_methylomics_overview.csv` - all datasets, including the ones not converted.

Load in the app: Methylomics > Dataset > *Upload your own data* > matrix `<GSE>_meth.csv` (orientation Feature x Sample,
first column = feature ID), sample sheet `<GSE>_sample.csv`; map Sample ID = `sample`, Group = `group`, Sex = `sex`.
Set *Input scale* per the table (GSE263432 is M-values, everything else beta).

| GSE | Category | Disease / design | group values (as converted) | sex | Platform | CpGs x samples | Scale | GEO fetch works? |
|---|---|---|---|---|---|---|---|---|
| GSE136724 | cancer | Chronic lymphocytic leukemia (blood) | no=29; yes=43 | F=23; M=49 | GPL13534 450K | 20,000 x 72 | beta (0-1) | yes |
| GSE138653 | immunological (treated-only) | Rheumatoid arthritis PBMC before adalimumab/etanercept (BiOCURA) | good=22; moderate=17; no=41 | F=56; M=24 | GPL21145 EPIC | 20,000 x 80 | beta (0-1) | yes |
| GSE151278 | immunological (treated-only) | Psoriasis peripheral blood on anti-TNF | ER=49; PR=21 | F=27; M=43 | GPL13534 450K | 20,000 x 70 | beta (0-1) | yes |
| GSE263432 | rare immunological | Sarcoidosis bronchoalveolar lavage cells | control=16; sarcoidosis (non-progressive)=36; sarcoidosis (progressive)=25 | F=40; M=37 | GPL21145 EPIC | 20,000 x 77 | M-values | yes |
| GSE110776 | rare immunological | Sarcoidosis BAL: remitting vs progressive | Treated=5; Untreated=19 | F=7; M=17 | GPL13534 450K | 20,000 x 24 | beta (0-1) | yes |
| GSE110778 | rare immunological | Chronic beryllium disease / beryllium sensitisation / sarcoidosis BAL | BeS=8; CBD=8; Sarcoidosis=8 | F=4; M=20 | GPL13534 450K | 20,000 x 24 | beta (0-1) | yes |
| GSE189426 | rare immunological | Undifferentiated arthritis blood/synovial monocytes | HD=15; UA=56 | F=43; M=28 | GPL21145 EPIC | 20,000 x 71 | beta (0-1) | yes |
| GSE175758 | cancer | Acute myeloid leukemia blasts on decitabine (DECIDER trial) | 0=52; 15=28; 8=52 | F=56; M=76 | GPL13534 450K | 20,000 x 132 | beta (0-1) | no |
| GSE193879 | rare immunological | Multisystem inflammatory syndrome in children (MIS-C) whole blood | Control=69; COVID-19=15; MIS-C=43 | F=58; M=69 | GPL21145 EPIC | 20,000 x 127 | beta (0-1) | no |
| GSE89251 | rare immunological | Juvenile idiopathic arthritis CD4+ T cells anti-TNF withdrawal | *not converted* - T0 vs Tend; ID / NO ID / flare | F=98 / M=38 | GPL13534 450K |  |  | yes |
| GSE274821 | cancer | B-cell NHL T cells before/after chemotherapy (CAR-T) | *not converted* - pre-therapy=16 / post-therapy=15; T cells=22 / CAR T=9 | M=24 / F=7 | GPL21145 EPIC |  |  | no |
| GSE315367 | cancer | Longitudinal AML diagnosis/remission/relapse | *not converted* - Diagnosis=14 / Remission1=14 / Remission2=14 / Relapse=4 | M=36 / F=10 | GPL21145 EPIC |  |  | no |
| GSE279950 | cancer | IDH-mutant low-grade glioma on temozolomide | *not converted* - all TMZ-treated; MGMT class M=91 / U=18; batch1/2 | M=62 / F=47 | GPL21145 EPIC |  |  | no |
| GSE184500 | immunological | Psoriatic arthritis CD8+ T cells | *not converted* - control=9 / psoriasis=14 / PsA=8; before/after treatment 5/5 | M=18 / F=13 | GPL21145 EPIC |  |  | no |
| GSE176168 | immunological (treated) | Rheumatoid arthritis anti-TNF week 0 vs 12 | *not converted* - wk0=62 / wk12=51 paired; response R=94 / NR=19 | F=97 / M=16 | GPL21145 EPIC |  |  | no |
| GSE252021 | rare immunological | Giant cell arteritis CD4+ T cells | *not converted* - CTRL=28 / active untreated=17 / treated remission=26 / untreated remission=27 | F=65 / M=33 | GPL21145 EPIC |  |  | no |

Notes:

- GSE136724: GEO's declared `Sex` has 33 unknowns, so `sex` comes from the depositor's `predicted_sex_rnbeads_package_in_bioconductor`
  column (the declared column is kept as an extra column). `group` = `treated` (yes/no); `response_to_therapy` is an extra column.
- GSE175758: `group` = sampling day (day0 = before decitabine, day8 / day15 = on treatment); `cell_type` (blasts / T_cells),
  `course` and the patient id (in `title`) are extra columns for paired or cell-type-restricted designs.
- GSE110776: GEO codes treatment as `Untreated`/`untreated`/`Treated`; the converted `group` is normalised to Untreated/Treated.
  The remitting-vs-progressive sarcoidosis label (GEO `group`, which collides with the canonical column) is restored from the sample titles as the extra column `sarcoidosis_course`.
- GSE263432: series matrix is deposited as **M-values**; the converter keeps them as-is - choose *M-values* as input scale.
- GSE110778: sex is 20 M / 4 F - fine for group comparisons, too unbalanced for the sex-interaction module.
- Series with `geo_fetch = no` (GSE193879, GSE175758) have empty GEO series matrices; the CSVs were built from the processed
  supplementary files after dropping detection p-value columns and renaming columns to GSM ids (`data/examples/methylomics_upload/geo_methylomics/geo_supp_to_gsm_columns.py`).
- Not converted (listed in the overview only, all need a supplementary download): GSE252021 giant cell arteritis CD4 (4-arm
  treated/untreated), GSE274821 lymphoma pre/post chemo, GSE176168 RA anti-TNF wk0/wk12, GSE184500 psoriatic arthritis, GSE315367
  longitudinal AML, GSE279950 glioma on temozolomide. GSE89251 (JIA) is already covered in the multiomics folder.

## Rebuilding

```
Rscript data/examples/methylomics_upload/geo_methylomics/build_geo_methylomics.R --gse GSE136724 --group treated \
    --sex predicted_sex_rnbeads_package_in_bioconductor --out data/examples/methylomics_upload/geo_methylomics --verify
Rscript data/examples/methylomics_upload/geo_methylomics/build_geo_methylomics.R --gse GSE193879 --supp GSE193879_betas_gsm.csv.gz \
    --group group --sex biological_sex --out data/examples/methylomics_upload/geo_methylomics --verify
Rscript data/examples/methylomics_upload/geo_methylomics/build_geo_methylomics.R --check data/examples/methylomics_upload/geo_methylomics
```
