# 01. Overview — Methylomics Cell-Type Deconvolution

**Source files read in full:** `ArthOMix/R/methylomics/mod_methyl_celltype.R` (1,130 lines), `ArthOMix/R/methylomics/celltype.R` (679 lines).
**Dependencies read:** `ArthOMix/R/methylomics/qc.R`, `ArthOMix/R/methylomics/annotation.R`, `ArthOMix/R/methylomics/parse_upload.R`, `ArthOMix/R/methylomics/mod_methyl_featureselection.R` (relevant functions), `ArthOMix/R/methylomics/mod_methyl_dataset.R`, `ArthOMix/global.R` (chart-styling section), `ArthOMix/R/submodules_registry.R`, `ArthOMix/server.R`. The installed `EpiDISH` R package's own exported source (`epidish()`, `hepidish()`, and internal `DoCP`/`DoRPC`/`DoCBS`) was also inspected directly in this R environment to verify claims about determinism and constraint handling.

**Registration:** `mod_methyl_celltype_config` — `id = "celltype"`, `title = "Cell-Type Deconvolution"`, `icon = "people-group"`, `group = "Data"` (`mod_methyl_celltype.R:24-27`). Registered third in `MX_MODULES` (`submodules_registry.R:39-42`), between Normalization and DMP (Differential Methylation), invoked generically as `mod_methyl_celltype_server("mx_celltype", methyl_dataset, methyl_results)` via `lapply(MX_MODULES, function(m) m$server(paste0("mx_", m$config$id), methyl_dataset, methyl_results))` (`server.R:95`).

---

## 1.1 Biological motivation (general scientific background, not a claim about this code)

Bulk DNA methylation array or bisulfite-sequencing data measures the average methylation level at each CpG site across every cell in a tissue sample (e.g. whole blood, a solid-tissue biopsy). Because different cell types (e.g. B cells, T-cell subsets, monocytes, neutrophils in blood; epithelial vs. fibroblast vs. immune cells in a solid tissue) have distinct, stable methylation signatures at specific "marker" CpGs, an observed bulk methylation profile can be modeled as a linear mixture of that sample's constituent cell types' own methylation profiles, weighted by their proportions. Cell-type deconvolution estimates those proportions from the bulk signal alone, without physically sorting cells. This is important because differences in cell-type composition between disease and control groups (e.g. more neutrophils, fewer lymphocytes in an inflammatory condition) are a well-known confounder of differential-methylation analysis — a CpG that appears "differentially methylated between groups" may in fact simply track a shift in cell-type composition rather than any genuine within-cell-type epigenetic change. Estimated cell fractions are therefore used both as a QC/confounder-adjustment covariate for downstream DMP/DMR analysis and as a biological readout in their own right (e.g. immune-cell shifts as a disease biomarker).

The classical approach (Houseman et al. 2012) is a reference-based **constrained projection**: given a reference matrix of per-cell-type mean methylation at a curated set of marker CpGs, and the observed bulk methylation at those same CpGs, solve a constrained (non-negative, and optionally sum-to-one) least-squares/quadratic-programming problem for the mixing proportions. Newer methods relax or replace the plain least-squares fit: **RPC** (robust partial correlations, via robust linear regression) is less sensitive to outlier CpGs; **CBS** (a CIBERSORT-style formulation, via nu-support-vector regression) was adapted for methylation data by the `EpiDISH` package; **hepidish** performs a **two-stage hierarchical** decomposition (e.g. split a solid tissue into Epithelial/Fibroblast/Immune first, then further resolve the Immune fraction into blood-cell subtypes with a second, independent reference).

## 1.2 What this module implements vs. general theory

This submodule implements exactly one estimation engine — the Bioconductor `EpiDISH` package, which is installed in this deployment — exposing its three methods (`CP` = Houseman-style constrained projection, `RPC`, `CBS`) plus its two-stage `hepidish()` function. It does **not** implement Houseman's original algorithm from scratch; `EpiDISH::epidish(method = "CP")` is itself a reimplementation of the same constrained-projection idea via `quadprog::solve.QP()` (verified by reading `EpiDISH:::DoCP` directly — see `03_functions_and_code_audit.md` and `05_statistical_methodology.md`).

Explicitly **not implemented**, and explicitly surfaced as disabled UI choices with a stated reason rather than silently omitted or faked (`celltype.R:68-77`, rendered at `mod_methyl_celltype.R:177-178`):
- **MethylResolver** — package not installed in this deployment.
- **IDOL-optimized reference libraries** (`IDOLOptimizedCpGs` / `FlowSorted.Blood.EPIC`) — not installed.
- **True reference-free deconvolution** (e.g. `RefFreeEWAS`, `TOAST`) — no such package installed.
- **Differential-methylation-marker (DMC/DMR)-based CpG selection**, offered as a feature-selection method choice but marked "(unavailable)" in its own label (`mod_methyl_celltype.R:99`), with an explicit reason: it would need per-sample sorted-cell-type data, which the built-in EpiDISH reference panels do not provide (mean-beta centroids only, no per-sample replicates) (`mod_methyl_celltype.R:120-124`).

This "backend honesty" convention is stated as this file's own design contract in its header comment (`celltype.R:12-18`) and is verified against the code throughout this documentation set rather than merely restated.

## 1.3 Input data type, shape, and required metadata

**Primary input:** a numeric **beta-value matrix**, rows = CpG probe IDs, columns = samples, values nominally in `[0, 1]` (fraction methylated). This orientation is enforced end-to-end: `methyl_parse_matrix()` (`parse_upload.R:11-37`) sets `rownames(m) <- probe_ids` from the file's first column and treats every other column as one sample; every filter in `qc.R` used by this module operates row-wise for probe-level logic (`rowMeans`, e.g. `methyl_filter_missing()`) and column-wise for sample-level logic; `methyl_ct_reconstruct()` and `methyl_ct_validation_metrics()` (`celltype.R:356-386`) likewise assume CpG-rows/sample-columns throughout.

**Scale handling:** the module never assumes the input is already on the beta scale. `methyl_ct_detect_scale()` (`celltype.R:90-104`) inspects a sampled quantile range and classifies the matrix as `"beta"` (already `[0,1]`), `"percent"` (`[0,100]`), or `"m"` (unbounded M-values), and no transformation happens automatically — the UI requires an explicit "Apply Transformation" click (`mod_methyl_celltype.R:379-398`).

**Two data-source paths** (`ct_data_source` radio, `mod_methyl_celltype.R:40-43`):
1. **Shared dataset** — reads `methyl_dataset$beta`/`$input_scale`/`$sample_sheet`/`$array_type`/`$source`/`$detp` populated by the Dataset tab (`mod_methyl_dataset.R`); if `input_scale == "m"` it is converted to beta via `methyl_ct_m_to_beta()` before use (`mod_methyl_celltype.R:406-410`).
2. **Module-scoped upload** — an independent CSV/TSV matrix (+ optional sample sheet), parsed the same way, scoped only to this module and not written back to `methyl_dataset` (`mod_methyl_celltype.R:325-363`).

**Reference input:** either one of 7 built-in EpiDISH reference centroid matrices (a CpG-by-cell-type matrix of mean beta values per cell type; `methyl_ct_reference_registry()`, `celltype.R:27-54`) or a user-uploaded custom reference matrix in the same CpG-by-cell-type shape, validated to be finite, in `[0, 1]`, with ≥2 cell-type columns (`methyl_ct_parse_custom_reference()`, `celltype.R:135-143`).

**Optional metadata:** a sample sheet / phenotype table (any columns), used for (a) the per-sample "missing sample" QC filter's row alignment, (b) coloring the PCA/MDS composition plots, and (c) the Group Comparison tab's grouping variable. No sample-sheet column is mandatory for the module to run — deconvolution itself needs only the beta matrix and a reference.

## 1.4 What the module produces

Per-CpG marker rankings (with effect size and specificity, deliberately no p-value/FDR — see `05_statistical_methodology.md`), an estimated samples-by-cell-types fraction matrix from one of CP/RPC/CBS/hepidish, cell-composition visualizations (stacked bar, heatmap, box/violin, PCA/MDS, correlation matrix), group-comparison statistics against an uploaded phenotype column, a reconstruction-validation of the estimated fractions against the observed methylation, a cross-method agreement comparison, and 6 CSV + 1 text-report export. Downstream, only `results$celltype` (a small summary list: method, cell types, sample count, markers used, mean fraction per cell type) is written into the shared `methyl_results` reactiveValues store (`mod_methyl_celltype.R:730-739`); see `04_data_flow_and_pipeline.md` for the exact consumer.
