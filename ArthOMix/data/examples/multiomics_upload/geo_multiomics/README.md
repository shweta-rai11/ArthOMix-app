# GEO multi-omics test datasets (matched transcriptomics + DNA methylation, same patients)

Built 2026-09-03 from NCBI GEO. Each SuperSeries below has an expression layer and a DNA-methylation layer
from the **same individuals**, plus sex and a disease/treatment grouping. Files per SuperSeries:

- `<GSE>_exp.csv` - genes in rows (`gene` column), samples in columns. Gene symbols (probes collapsed by the
  same MaxMean rule as the app's `collapse_probes_to_genes()`), values exactly as deposited in GEO.
- `<GSE>_meth.csv` - CpGs in rows (`cpg` column), samples in columns. **Top 20,000 most-variable CpGs** only
  (full arrays are 27K-866K rows; the variance filter keeps files small, and does not exclude X/Y probes).
- `<GSE>_sample.csv` - one row per individual sample; `sample` is the **shared ID used as the column header
  in both layers**, plus `gsm_exp`, `gsm_meth`, `group`, `sex` (F/M/NA), `in_expression`, `in_methylation`,
  and dataset-specific extras. Rows present in only one layer have NA in the other layer's columns.

Because both layers use the same sample IDs, the multiomics Dataset tab's **"Exact Sample ID"** matching
works with no mapping file. Orientation is Feature x Sample (same as the bundled gse138746/gse138653 example).

| SuperSeries | Disease / design | group values | sex F/M | Expression (SubSeries, scale) | Methylation (SubSeries, scale) | Shared samples | Shared sample ID looks like |
|---|---|---|---|---|---|---|---|
| GSE89253 | Juvenile idiopathic arthritis CD4+ T cells, anti-TNF withdrawal; baseline T0 vs Tend | ID / NO ID / flare; also `timepoint` T0/Tend and `patient` | 95/30 | GSE89252, 15,209 genes x 125, log2 (4.3-14.6) | GSE89251, 450K, **M-values** (-14.6 to 15.0), 136 samples | 125 (11 methylation-only) | `13106_T0`, `13106_Tend` |
| GSE146093 | Systemic sclerosis CD4+ T cells, 48 SSc vs 16 healthy | HD / SSc; `subtype` lcSSc/dcSSc/ssSSc/Initial; `lung_disease` | 56/8 | GSE146088, 19,284 genes x 64, log2 (1.9-19.6) | GSE146092, EPIC, beta (0-1), 64 samples | 64 (all) | `H4` (healthy), `S7` (SSc) |
| GSE201754 | Giant cell arteritis CD14+ monocytes; healthy vs active vs remission on glucocorticoids vs remission untreated | HD / ACT / RT / RNT; `group_label` | 68/40 | GSE201753, 60,591 genes x 111, **raw RNA-seq counts** | GSE201752, EPIC, beta (0-1), 113 samples | 111 (2 methylation-only) | `V1M`, `V100M` |
| GSE32867 | Lung adenocarcinoma, tumour vs matched adjacent normal | Tumor / Normal; `subject`, `stage`, `smoking` | 90/26 | GSE32863, 38,275 genes x 116, log2 (6.3-16.1) | GSE32861, 27K, beta (0-1), 118 samples | 116 (2 methylation-only) | `05L12_N`, `05L12_T` |
| GSE117931 | Systemic sclerosis PBMC, 18 SSc vs 19 control | NC / SSc | 25/12 | GSE117928, 20,815 genes x 37, **quantile-normalised, NOT log** (-23 to 29,630) | GSE117929, 450K, beta (0-1), 37 samples | 37 (all) | `NC001`, `SSc012` |
| GSE65205 | Childhood atopic asthma, nasal epithelium | Asthma / Control | 35/34 | GSE65204, 13,865 genes x 69, log2 (3.0-16.0) | GSE65163, 450K subset (65,535 CpGs deposited), **M-values** | 69 (3 methylation-only) | `15-02-005-0` |
| GSE50387 | Seasonal allergic rhinitis CD4+ T cells, in vs outside pollen season | HC / Allergic; `season`, `patient` | 18/20 | GSE50101, 34,686 genes x 38, log2 (5.0-14.8) | GSE50222, 450K, beta (0-1), 32 samples | 30 (8 expression-only, 2 methylation-only) | `11H_during`, `02P_outside` |
| GSE87650 | Inflammatory bowel disease, whole blood + CD4/CD8/CD14 | UC / CD / HC; `cell_type` | 92/159 | GSE86434, 34,412 genes x 251, log2 | **not built** - GSE87640/GSE87648 use a different sample coding and a 944 MB matrix | expression only | GSM ids |

Notes on specific datasets:
- GSE117931: the expression SubSeries has no sex in GEO; `sex` was copied from the methylation SubSeries by
  sample title (`sex_source` column). Its expression matrix has negative values, so in the transcriptomics
  Upload tab only "Normalized" passes the validator.
- GSE146093: methylation samples V25T/V26T are healthy donors 25/26 (GEO characteristics say Healthy donor),
  so they were keyed H25/H26 to match the expression layer.
- GSE89253 and GSE65205 methylation are M-values, not betas; the app's value-type detector reports this.
- GSE201754: the expression layer is raw counts; 17,182 genes are all-zero.

## Using "Retrieve from GEO" in the app instead of these files

- Enter the **SubSeries** id, not the SuperSeries. A SuperSeries id also fetches, but GEO serves one
  series-matrix per platform, so the app shows a platform picker containing both the expression and the
  methylation platform. Picking the methylation one loads CpG betas as "expression".
- Empty series matrices (fetch reports it; upload required): GSE201753 (RNA-seq), GSE146092 and GSE201752
  (EPIC). These CSVs were built from their supplementary files.
- Illumina expression platforms GPL10558 / GPL6884 / GPL14951 name their symbol column `Symbol`, which the
  app's `collapse_probes_to_genes()` regex `^gene[ ._]?symbol$` does not match, so a live fetch of GSE89252,
  GSE32863, GSE86434, GSE117928 and GSE50101 lands at probe level without gene symbols. GSE65204 (Agilent,
  `GENE_SYMBOL`) collapses correctly. GSE146088 (Clariom S) has no symbols in GEO's platform table; these
  CSVs used Bioconductor `clariomshumantranscriptcluster.db`.
- Methylation series matrices with values (fetchable live by the methylomics tab): GSE89251, GSE32861,
  GSE117929, GSE65163, GSE50222.

## Same-individual verification (2026-09-03)

For every shared sample ID, the sex, age and group recorded by GEO on the expression SubSeries were compared
with those recorded independently on the methylation SubSeries. Counts are "agree / comparable".

| SuperSeries | Shared samples | Sex | Age | Group |
|---|---|---|---|---|
| GSE89253 | 125 | 125/125 | 125/125 | 125/125 |
| GSE146093 | 64 | 64/64 | not on GEO | 64/64 |
| GSE201754 | 111 | 108/108 (3 unannotated) | 111/111 | 111/111 |
| GSE32867 | 116 | 116/116 | 116/116 | 116/116 |
| GSE117931 | 37 | copied from methylation layer, so not an independent check | not on GEO | 37/37 (+ identical sample titles) |
| GSE65205 | 69 | 69/69 | 69/69 | 69/69 |
| GSE50387 | 30 | 30/30 | not on GEO | 30/30 |
| GSE87650 | 0 | no methylation layer built; expression SubSeries carries no patient id | | |

## Converting another GEO series yourself

`data/examples/geo_to_arthomix_upload.R` is a ~30-line script that turns any GEO series into this upload format:

```
Rscript data/examples/geo_to_arthomix_upload.R GSE89252 "clinical activity" sex          # expression layer
Rscript data/examples/geo_to_arthomix_upload.R GSE89251 "clinical activity" sex          # methylation layer (auto-detected)
```

Arguments: GSE id, group characteristic, sex characteristic (names as shown on the GEO sample page), optional
output folder (default `data/uploads/<GSE>`). It writes `<GSE>_exp.csv` or `<GSE>_meth.csv` plus `<GSE>_sample.csv`.
Series whose GEO matrix is empty (RNA-seq, EPIC) need their processed supplementary file instead.
