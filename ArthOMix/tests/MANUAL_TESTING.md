Manual testing for Shiny App using sample data
- why? shinytest/testServer coverage exists (see tests/testthat) but has hard-to-debug
  gaps for full end-to-end reactive chains (file upload -> validation -> multi-tab
  pipeline state); manual pass confirms real user-facing behavior across all four
  omics verticals.

Go to: run `Rscript -e "shiny::runApp('.', launch.browser = TRUE)"` from the ArthOMix/
       repo root and open the local URL Shiny prints in the console. Run
       `Rscript -e "renv::restore()"` first if packages aren't installed yet.

Manual test Log:
- DATE: __________, run by: __________ — mark each test PASS / FAIL, note observed values.


================================================================================
TRANSCRIPTOMICS  (default cohort: GSE93272 + GSE110169 merged, RA whole blood)
================================================================================

1) Dataset (01_Data)
objectives
- confirm all three load paths work and gate correctly
test 1 — Preloaded                                                    [ ] PASS  [ ] FAIL
- input: open Dataset tab, pick a cohort from dropdown (GSE93272 "Whole Blood Training
  Cohort A", GSE110169 "Whole Blood Training Cohort B", GSE15573 "PBMC Validation Cohort",
  or GSE89408 "Synovial Tissue Validation Cohort"), click "Load this dataset"
- output: dataset loads, downstream tabs unlock
- notes:

test 2 — GEO fetch                                                    [ ] PASS  [ ] FAIL
- input: enter a GSE accession in the textbox, click "Fetch", then "Load this dataset"
- output: dataset fetched from NCBI GEO and loads same as preloaded path
- notes:

test 3 — Upload                                                       [ ] PASS  [ ] FAIL
- input: use example files in ArthOMix/data/examples/transcriptomics_upload/merged/
  (or probelevel/ subfolder); Step 1 — pick declared_data_type radio (Raw counts /
  Normalized TPM-FPKM-CPM / Already log-transformed) matching the file, upload the
  expression matrix CSV and sample metadata CSV
- output: "Upload Data" button stays disabled until both files are uploaded
- notes:

test 4                                                                 [ ] PASS  [ ] FAIL
- input: Step 2 — confirm auto-guessed column mapping dropdowns (Sample ID, Group/
  diagnosis, Sex, Batch) populate; edit if wrong; click "Upload Data"
- output: dataset accepted only if >=4 matching sample IDs between expression and
  metadata; duplicate sample IDs block upload; duplicate gene/feature IDs produce a
  warning but are kept
- notes:

test 5 — mismatch guard                                               [ ] PASS  [ ] FAIL
- input: intentionally pick the wrong declared_data_type for the file's actual value range
- output: upload is hard-blocked with a validation error
- notes:

2) Overview and Datasets (02_Overview)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: load default cohort
- output: QC summary of loaded dataset renders (sample counts, group breakdown, etc.)
- notes:

3) Preprocessing / Batch Correction (03_Preprocessing_Batch_Correction)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Preprocess this dataset" / "Preprocess"
- output: normalization + batch-correction pipeline runs, diagnostic plots update
- notes:

test 2                                                                 [ ] PASS  [ ] FAIL
- input: click "Load and Preprocess Selected Cohorts"
- output: multi-cohort merge + correction executes
- notes:

test 3 — known broken path                                            [ ] PASS  [ ] FAIL
- input: look for a "merge example training datasets" demo/example-merge option
- output: KNOWN ISSUE — this demo path has a dead symlink and is expected to fail;
  confirm it still fails (or has since been fixed) and record which
- notes:

4) Differential Expression (04_Differential_Expression)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run differential expression"
- output: DEG table + volcano plot render; results should be sex-stratified per app design
- notes:

5) WGCNA (05_WGCNA)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run Soft Threshold Analysis"
- output: soft-threshold diagnostic plot renders
- notes:

test 2                                                                 [ ] PASS  [ ] FAIL
- input: click "Run WGCNA"
- output: modules detected, dendrogram/module-color plot renders
- notes:

test 3                                                                 [ ] PASS  [ ] FAIL
- input: click "Compute module-trait correlations", then "Compute hub genes"
- output: module-trait heatmap and hub-gene table render
- notes:

test 4                                                                 [ ] PASS  [ ] FAIL
- input: click "Run enrichment", then "Compute network heatmap"
- output: enrichment table and network heatmap render
- notes:

6) Candidate Gene Identification (06_Candidate_Gene_Identification)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: with DEG + WGCNA both run, open Candidates tab
- output: Venn/overlap of WGCNA hub genes and DEGs renders, candidate gene list populates
- notes:

7) Mendelian Randomization (07_Mendelian_Randomization)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run Clumping", then "Run MR", then "Run Sensitivity"
- output: MR effect-estimate table + sensitivity diagnostics (e.g. heterogeneity/pleiotropy
  tests) render; uses Okada 2014 RA GWAS instruments
- notes:

8) Colocalization (08_Colocalization)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run Colocalisation"
- output: coloc posterior probabilities table/plot renders per candidate region
- notes:

9) Feature Selection (09_Feature_Selection)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click each of "Run Univariate Selection", "Run LASSO / Elastic Net",
  "Run Random Forest", "Run RFE", "Run Stability Selection", "Run Consensus Selection"
  in turn, for both sexes
- output: each produces a selected-feature list; Consensus Selection combines across methods
- notes:

10) Diagnostic Model (10_Diagnostic_Model)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run Female", "Run Male", "Run All (pooled)"
- output: classifier (logistic/elastic net/RF/SVM) fits per stratum, performance metrics
  (ROC/AUC) render for each
- notes:

11) Sex Interaction Analysis (11_Sex_Interaction_Analysis)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run interaction model"
- output: sex-by-exposure interaction term results table renders
- notes:

12) Cross-Tissue Validation (12_Cross_Tissue_Validation)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run [Sex]" for the 4-classifier panel (uses GSE89408 synovial tissue as
  validation cohort)
- output: cross-tissue performance metrics render
- notes:

test 2                                                                 [ ] PASS  [ ] FAIL
- input: click "Run [Sex] Leakage-safe Validation"
- output: leakage-controlled validation metrics render, should differ from test 1 if
  leakage was present
- notes:

13) Cross-Ancestral Validation (13_Cross_Ancestral_Validation)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: open tab (uses cached cross-ancestry MR35 replication tables)
- output: replication concordance table/plot renders
- notes:

14) Functional Enrichment (14_Functional_Enrichment)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run enrichment" on a candidate gene list (GO/KEGG/Reactome ORA)
- output: enrichment results table + plot render
- notes:

15) Immune Deconvolution (15_Immune_Deconvolution)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Estimate cell composition"
- output: CIBERSORT/LM22 and MCP-counter cell-fraction estimates render per sample
- notes:

16) Nomogram (16_Nomogram)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Build nomogram"
- output: nomogram plot + decision curve analysis render
- notes:

17) Biomarker Card (17_Biomarker_Card)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: enter a gene from the candidate list (or the bundled ferroptosis panel,
  ArthOMix/data/annotations/gene_panels/ferroptosis_soorya_sundararajan_2025.txt),
  click each of "Run STRING/Open Targets/HPA/DGIdb/PubMed/GO/KEGG/Reactome/WikiPathways
  Query"
- output: each query populates its respective panel (interaction network, tractability,
  tissue expression, drug-gene interactions, literature, pathway memberships)
- known guard: WikiPathways panel shows an unavailable message if `msigdbr` isn't installed
- notes:


================================================================================
METHYLOMICS  (default cohort: GSE42861, Illumina 450K, whole blood, sex-stratified RA)
================================================================================

1) Dataset (01_Data)
test 1 — Preloaded                                                    [ ] PASS  [ ] FAIL
- input: pick GSE42861 preloaded cohort, click "Load this dataset"
- output: dataset loads
- notes:

test 2 — Upload (matrix)                                              [ ] PASS  [ ] FAIL
- input: use ArthOMix/data/examples/methylomics_upload/gse71841_expression_matrix.csv +
  gse71841_sample_metadata.csv (GSE71841, CD4+ T cells, 12 RA / 12 control, 485,577
  probes x 24 samples); pick array_type "EPIC" or "450K" as appropriate, upload_format
  "Beta/M-value matrix", input_scale "Beta values 0-1"
- output: dataset validates and loads
- notes:

test 3 — Upload (IDAT)                                                [ ] PASS  [ ] FAIL
- input: if raw .idat pairs available, use the "Raw IDAT files" upload_format instead
- output: beta values, detection p-values, bead counts derived automatically
- notes:

2) Quality Control (02_Quality_Control)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run Outlier/Sample/Probe/Sex/Batch QC" (buttons may be separate)
- output: QC plots/tables render for each check; sex-prediction-vs-reported-sex mismatch
  flagged if present
- notes:

test 2 — deployment guard                                             [ ] PASS  [ ] FAIL
- input: attempt to export the PDF QC report
- output: if `rmarkdown`/LaTeX not installed in this deployment, expect
  "not available in this deployment" message instead of a PDF
- notes:

3) Normalisation (03_Normalisation)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run Normalization"
- output: normalized beta/M-value distributions render (before/after density plots)
- notes:

4) Cell-Type Deconvolution (04_Cell_Type_Deconvolution)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run Cell-Type Deconvolution"
- output: estimated cell-type proportions per sample render
- notes:

test 2                                                                 [ ] PASS  [ ] FAIL
- input: click "Run Validation", then "Compare Methods"
- output: cross-method concordance comparison renders
- known guard: "Detection p-values require raw IDAT input - not available for this
  dataset" if loaded from a matrix (not IDAT) upload
- notes:

5) DMP (05_DMP)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run SVA-adjusted Analysis", then "Run DMP Analysis"
- output: differentially methylated position table + volcano plot render
- notes:

6) DMR (06_DMR)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run DMR Analysis"
- output: differentially methylated region table renders
- notes:

7) Sex Interaction Analysis (07_Sex_Interaction_Analysis)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run interaction model" equivalent for methylation
- output: sex-interaction results table renders
- notes:

8) WGCNA Co-Methylation Network (08_WGCNA_Co_Methylation_Network)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: run network construction, then click "Compute hub CpGs"
- output: co-methylation modules + hub CpG table render
- notes:

test 2                                                                 [ ] PASS  [ ] FAIL
- input: click "Run Functional Enrichment"
- output: enrichment results for module CpG-associated genes render
- notes:

9) Candidate CpGs (09_Candidate_CpGs)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run DMR-CpG Overlap", then "Run Module-DMR Overlap"
- output: overlap tables/Venn between DMPs, DMRs, and WGCNA modules render
- notes:

10) ML Feature Selection (10_ML_Feature_Selection)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: run each feature-selection method available (mirrors transcriptomics: LASSO/RF/
  stability/consensus)
- output: selected CpG feature lists render per method
- notes:

11) Mendelian Randomization (11_Mendelian_Randomization)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: run mQTL-MR pipeline
- output: MR results table renders
- known guard: "Needs exposure and outcome sample sizes - not available for this run" if
  sample size fields are unset
- notes:

12) Colocalization (12_Colocalization)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: run coloc on candidate CpG regions
- output: coloc posterior probability table/plot renders
- notes:

13) Diagnostic Classifier (13_Diagnostic_Classifier)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run Single-CpG Diagnostic Analysis"
- output: per-CpG diagnostic performance (ROC/AUC) renders
- notes:

test 2                                                                 [ ] PASS  [ ] FAIL
- input: click "Run Model" (multi-CpG panel classifier)
- output: panel-level classifier performance renders
- notes:

14) Validation (14_Validation)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run External Validation" (uses external GSE111942, 21-CpG panel)
- output: external validation performance metrics render
- known guard: methylomics results/candidate tabs show "The preloaded methylomics
  results folder is not available in this deployment - use Upload Data instead" if
  `METH_DATA_AVAILABLE` is false — confirm which applies
- notes:

15) Biomarker Card (15_Biomarker_Analysis)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: enter a candidate CpG/gene, click "Run EWAS Catalog/Atlas/KEGG/Reactome/
  WikiPathways Query"
- output: each query populates its panel
- known guard: panels needing `ChAMPdata`, `org.Hs.eg.db`/`httr2`, or `msigdbr` show an
  unavailable message if those packages aren't installed
- notes:


================================================================================
MULTI-OMICS  (default cohort: Tao et al. 2021 anti-TNF RA cohort, Adalimumab/Etanercept)
================================================================================

1) Dataset Workspace (01_Data_Workspace)
test 1 — Preloaded                                                    [ ] PASS  [ ] FAIL
- input: pick "Analysis cell (matched sex x drug/outcome subset)" dropdown, click
  "Load Reference Dataset"
- output: matched transcriptomics + methylomics subset loads
- notes:

test 2 — Upload                                                       [ ] PASS  [ ] FAIL
- input: use ArthOMix/data/examples/multiomics_upload/ files (gse138746_pbmc_rnaseq_
  counts.csv, gse138653_methylation_beta_top2000.csv, gse138747_sample_metadata.csv);
  set Data orientation and Table shape radios correctly for each block; use "+ Add
  Another Dataset" if needed; pick a sample-matching method
- output: click "Validate Datasets" — validation passes/fails with a clear message;
  click "Use Selected Datasets for Multi-Omics Analysis" to proceed
- notes:

test 3                                                                 [ ] PASS  [ ] FAIL
- input: click "Apply normalization, filtering, and scaling", then
  "Apply batch correction"
- output: processed data summary/plots update
- notes:

test 4 — MOFA2 (mounted in this tab)                                  [ ] PASS  [ ] FAIL
- input: click "Train MOFA2"
- output: MOFA2 factor model fits, variance-explained plot renders
- notes:

2) Cohort Harmonization (02_Cohort_Harmonization)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Cohort Harmonization", then "Analyze Cohort"
- output: harmonized cross-cohort summary renders
- notes:

3) DIABLO / SNF Integration (03_DIABLO_SNF_Integration)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run DIABLO"
- output: DIABLO block-integration fit + component plots render
- notes:

test 2                                                                 [ ] PASS  [ ] FAIL
- input: click "Run SNF"
- output: similarity-network-fusion result renders
- notes:

test 3                                                                 [ ] PASS  [ ] FAIL
- input: click "Run Comparison"
- output: DIABLO vs SNF comparison view renders
- known guard: "Feature-selection stability is not available for this configuration
  (needs more than one CV repeat)" if only 1 CV repeat configured
- notes:

4) SNF Clustering (04_SNF_Clustering)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run SNF Clustering"
- output: patient similarity clusters render
- notes:

test 2                                                                 [ ] PASS  [ ] FAIL
- input: click "Recompute stability", then "Run parameter sensitivity check"
- output: cluster stability metrics and sensitivity plot render
- notes:

5) Biomarker Discovery (05_Biomarker_Discovery)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run Sex-Stratified Analysis"
- output: sex-stratified multi-omic biomarker candidates table renders
- notes:

6) Gene-CpG Concordance (06_Gene_CpG_Concordance)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run Gene-CpG Analysis"
- output: gene-expression vs CpG-methylation concordance table/plot renders (feeds the
  read-only Biomarker Card evidence view)
- notes:

7) Pathways (07_Pathways)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run Pathway Analysis"
- output: multi-omic pathway enrichment table renders
- notes:

8) Biomarker Card (08_Biomarker_Card)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: look up a candidate biomarker from prior steps
- output: integrated evidence card renders
- known guard: "Patient-level values are not available for preloaded cohort tables
  (summary statistics only)" when using a preloaded (not uploaded) cohort
- notes:

9) Results Summary (09_Results_Summary)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: open tab after running several prior steps
- output: session summary renders with a downloadable results bundle; verify the
  download button actually produces a file
- notes:


================================================================================
CROSS-OMICS  (default: precomputed eQTL-MR x mQTL-MR x DEG x DMP/DMR convergence tables)
================================================================================

1) Dataset (01_Data)
test 1 — Example data                                                 [ ] PASS  [ ] FAIL
- input: source_mode "Example data", pick sex_stratum (ALL/FEMALE/MALE) and meth_level
  (CpG-level DMP / Region-level DMR), click "Load example data"
- output: dataset loads; on success, notification "Ready for Expression and
  Methylation." appears
- notes:

test 2 — Upload                                                       [ ] PASS  [ ] FAIL
- input: source_mode "Upload your own data"; use
  ArthOMix/data/examples/crossomics_upload/ files for transcriptomics (DEG: gene
  symbol + log2FC) and methylomics (DMP/DMR: gene symbol + delta-beta); columns
  auto-detected (no manual mapping step)
- output: click "Use this data" — loads successfully; "Clear" resets the tab
- notes:

2) Expression and Methylation Integration (02_Expression_Methylation_Integration)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: click "Run Integration"
- output: joint DEG/DMP-DMR integration table + concordance plot render
- notes:

3) Biomarker Convergence (03_Biomarker_Convergence)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: open tab (loads precomputed eQTL-MR/mQTL-MR evidence, no run button expected)
- output: convergence evidence table renders
- notes:

4) Cross-Omics MR (04_Cross_Omics_MR)
test 1                                                                 [ ] PASS  [ ] FAIL
- input: open tab, review classification of DEG/DMP/DMR/QTL evidence into the 5
  convergence categories
- output: classified convergence table renders
- known guard: "Cross-Omics MR source data is not available in this deployment" if
  `CX_MR_DATA_AVAILABLE` is false
- notes:


================================================================================
ARTHOCHAT (persistent slide-out drawer, not a tab)
================================================================================
test 1                                                                 [ ] PASS  [ ] FAIL
- input: open ArthOChat from outside a valid omics/module context (e.g. Home tab)
- output: expect error notification "ArthOChat is not available in this context."
- notes:

test 2                                                                 [ ] PASS  [ ] FAIL
- input: open ArthOChat from within a loaded module, ask a PubMed-grounded question
- output: response renders with cited sources
- notes:


================================================================================
NOTES FOR THE TESTER
================================================================================
- Every "not available" message listed above is conditional on deployment (missing
  optional R package, missing bundled results folder, or insufficient sample size) —
  not a permanent disabled feature. Confirm the actual on-screen message matches
  current environment state rather than assuming it will always fire.
- No numeric expected outputs (p-values, AUCs, etc.) are given here since they weren't
  captured from a live run — record actual observed values above on first PASS run;
  they can be pinned as regression baselines for future manual passes.
- The Preprocessing "merge example training datasets" demo option is flagged in the
  data manifest as a dead symlink as of the last codebase check — verify current status.
- Work top to bottom within each vertical: Dataset tab first, then submodules in
  numbered order, since later tabs (WGCNA, MR, Diagnostic Model, etc.) read state
  produced by earlier ones.
