# ARTHOMIX — Review annotated against the current code (2026-09-03)

**Legend.** <span style="color:red">**Red = corrected in the app code**</span> (verified against the code and by tests). <span style="color:#d17a00">**Orange = partly corrected**</span> (what is still open is stated). Plain text = not a code item (thesis/manuscript wording), or not addressed.

---

**ARTHOMIX**

**Scientific, Statistical, and Software Review**

Recommendations for a Shiny-based single- and multi-omics web application

| **Document reviewed**       | **Review focus**                                                                                                  |
| --------------------------- | ----------------------------------------------------------------------------------------------------------------- |
| ArthOMix_draft1_30JUNE.docx | Scientific validity, analytical design, software architecture, usability, reproducibility, and manuscript clarity |

Overall assessment. ArthOMix is a promising and potentially publishable application concept. Its strongest prospective differentiators are sex-aware analysis, combined transcriptomic–methylomic interpretation, blood–synovium comparison, integrated biomarker prioritization, and gene-level contextualization.

# Executive assessment

The central concept is strong, but the current description is more ambitious than the analytical design presently supports. The draft combines a software specification, a methods manuscript, and an RA case study. These should be separated conceptually before extensive interface development so that the application does not produce attractive visualizations without sufficient statistical safeguards.

**Recommended positioning:**

**_A modular, reproducible web application for transcriptomic and DNA-methylation analysis, sex-aware investigation, and matched-sample multi-omics integration, demonstrated using rheumatoid arthritis datasets._**

The team should decide whether ArthOMix is (1) a general omics platform using RA as a validation case study or (2) an RA-specific discovery platform. The first option is likely to have broader scientific and software value. If that direction is chosen, RA-specific annotations should be implemented as an optional disease-context layer rather than embedded throughout the core workflow.

*(Thesis decision, not a code item — still to be stated in the manuscript.)*

# Principal strengths

- An end-to-end workflow extending from quality control to biological interpretation.
- Explicit attention to sex, a scientifically important dimension in RA that could become a genuine differentiator.
- Complementary unsupervised and supervised integration through MOFA and DIABLO.
- Blood–synovium comparison, which connects biomarker accessibility with disease-relevant tissue biology.
- Gene cards that can make results accessible to researchers without extensive bioinformatics expertise.
- Potential for local or private-server deployment when unpublished or sensitive datasets are analyzed.

# Priority findings

| **Priority** | **Issue** | **Required response** | **Status in code** |
| --- | --- | --- | --- |
| Critical | <span style="color:red">Sex-stratified analysis does not test whether effects differ by sex.</span> | <span style="color:red">Add disease × sex interaction models and within-sex contrasts.</span> | <span style="color:red">**Corrected.** Sex Interaction Analysis (transcriptomics and methylomics) fits `outcome ~ group × sex + covariates` and reports the interaction, the disease effect within each sex, and the main sex effect from the same fit; multiple covariates; rank (estimability) check; per-cell sample counts; underpowered warning.</span> |
| Critical | <span style="color:red">Multi-omics methods are not clearly separated by sample-pairing structure.</span> | <span style="color:red">Use matched feature matrices for MOFA/DIABLO; use cohort-level synthesis for unpaired studies.</span> | <span style="color:red">**Corrected.** Multi-omics module takes matched sample × feature matrices with an explicit sample-matching step and overlap table feeding MOFA/DIABLO/SNF; Cross-omics module is the separate result-table (unpaired) pathway.</span> |
| Critical | <span style="color:red">Machine-learning workflow risks information leakage.</span> | <span style="color:red">Move preprocessing, feature selection, and tuning inside nested validation.</span> | <span style="color:red">**Corrected.** Held-out split taken before feature selection; CV-only tuning; train-only scaling; automatic nested-CV headline metric with on-screen disclosure when a run is not leakage-safe; trained models applied unchanged to an external cohort.</span> |
| High | <span style="color:red">Input matrices are treated too generically.</span> | <span style="color:red">Define typed input schemas and expose only statistically valid workflows.</span> | <span style="color:red">**Corrected.** Declared data type validated against the matrix (raw counts / normalised / log; beta / M / IDAT); result tables rejected as expression matrices; count models only for counts.</span> |
| High | <span style="color:red">Batch correction and study integration are oversimplified.</span> | <span style="color:red">Distinguish compatible pooled analysis from study-wise meta-analysis.</span> | <span style="color:red">**Corrected.** Batch × phenotype confounding check with a hard block and explicit override in transcriptomics, methylomics and multi-omics; new Study-wise Meta-analysis tab (per-study limma / limma-voom, DerSimonian–Laird random effects, τ², Q, I²).</span> |
| High | <span style="color:red">Causal language exceeds the proposed evidence.</span> | <span style="color:red">Replace causation with regulatory inference or add a defensible causal design.</span> | <span style="color:red">**Corrected.** Causal wording confined to MR and colocalisation; cross-omics and biomarker-card text state "associated with", never "causes".</span> |
| Medium | <span style="color:#d17a00">Software architecture and reproducibility requirements are not specified.</span> | <span style="color:#d17a00">Add session isolation, provenance, versioning, testing, security, and deployment details.</span> | <span style="color:#d17a00">**Partly.** Present: provenance manifests, fixed seeds, `renv.lock` + `DESCRIPTION` + `Dockerfile`, 94 test files, async jobs, upload-size cap, per-session state, disk caching limited to the bundled reference data. Open: job cancellation, concurrent-user/memory limits, per-run R-script export, accessibility markup.</span> |

# Scientific and analytical recommendations

## 1. Replace simple sex stratification with sex-aware modelling

Analyzing RA versus control separately in males and females does not demonstrate a sex difference. A finding can be significant in females and nonsignificant in males even when the underlying effects are statistically indistinguishable.

The primary model should include an interaction term:

<span style="color:red">**Molecular outcome ~ disease + sex + disease × sex + covariates**</span> — <span style="color:red">implemented in both omics.</span>

The application should report:

- <span style="color:red">overall disease effects;</span>
- <span style="color:red">overall sex effects;</span>
- <span style="color:red">disease-by-sex interaction effects;</span>
- <span style="color:red">within-sex disease contrasts as secondary analyses.</span>

<span style="color:red">The same framework should be available for expression and methylation.</span> <span style="color:red">User-selectable covariates should include age, batch, cell composition, medication, smoking, ancestry, tissue source, and relevant technical variables</span> (any metadata / sample-sheet column can be selected, several at once). <span style="color:red">The application should also check sample counts and model estimability</span> (per-cell counts, design-rank check, minimum-detectable-effect warning when a cell has fewer than 10 samples), particularly because some public RA datasets may contain too few male samples for reliable interaction analysis.

## 2. Make input type determine the analytical pathway

<span style="color:red">ArthOMix should not treat all uploaded matrices as equivalent. It should require or infer a declared data type and expose only valid preprocessing and modelling choices.</span>

| **Input type** | **Key handling requirement** | **Status** |
| --- | --- | --- |
| Raw RNA-seq counts | Count-aware filtering, normalization, and modelling; DESeq2, edgeR, or limma–voom. | <span style="color:red">Corrected (TMM / ComBat-seq / DESeq2 / limma-voom paths).</span> |
| Normalized RNA-seq expression | Do not apply count models; require provenance and transformation details. | <span style="color:red">Corrected (declared type checked against values).</span> |
| Microarray intensities | Platform-appropriate background correction, normalization, annotation, and QC. | <span style="color:#d17a00">Partly (GEO platform detection, quantile normalisation, probe collapse; no explicit microarray-vs-RNA-seq declared type).</span> |
| Processed microarray expression | Confirm scale, normalization status, probe mapping, and batch variables. | <span style="color:red">Corrected.</span> |
| Methylation IDAT files | Platform-specific minfi-style QC and preprocessing. | <span style="color:red">Corrected.</span> |
| Methylation beta values | Use for interpretation/visualization; model choice must be explicit. | <span style="color:red">Corrected.</span> |
| Methylation M-values | Preferred for many differential tests; retain beta values for interpretable effect display. | <span style="color:red">Corrected (models on M-values, Δβ reported).</span> |
| Result tables | Permit downstream interpretation only; do not present them as sample-level omics matrices. | <span style="color:red">Corrected (rejected as matrices; accepted only by Cross-omics).</span> |

## 3. Introduce a defensible batch and study-integration policy

Batch effects can occur within a single study, whereas indiscriminate correction after merging studies can remove genuine biological differences. ArthOMix should distinguish two strategies:

- <span style="color:red">Pooled analysis when datasets are genuinely compatible, with study and batch represented in the design model.</span>
- <span style="color:red">Study-wise analysis followed by effect-size meta-analysis when platforms, populations, tissues, or processing procedures differ.</span>

<span style="color:red">The interface should show pre- and post-adjustment PCA or similar diagnostics and warn when batch is confounded with disease or sex.</span> In a fully confounded design, statistical correction cannot reliably separate the technical and biological effects.

## 4. Redesign the multi-omics entry point around sample pairing

<span style="color:red">Differential-expression and differential-methylation result tables are not sufficient inputs for sample-level correlation, MOFA, or DIABLO. Paired integration requires appropriately processed feature-by-sample matrices, an explicit sample-matching table, and substantial overlap in measured individuals.</span>

| **Data relationship** | **Appropriate analyses** | **Status** |
| --- | --- | --- |
| Matched or substantially overlapping individuals | Sample-level expression–methylation correlation, MOFA, DIABLO, joint latent factors, and multiblock feature selection. | <span style="color:red">Corrected (Multi-omics module).</span> |
| Unpaired cohorts from the same condition | Independent differential analysis, pathway-level integration, rank aggregation, effect-direction concordance, meta-analysis, and network/knowledge-based integration. | <span style="color:red">Corrected (Cross-omics module; study-wise meta-analysis).</span> |
| Different tissues | Cross-tissue replication or projection with explicit handling of tissue and cell-composition differences. | <span style="color:red">Corrected (Cross-Tissue Validation is a separate replication module, not paired integration).</span> |

## 5. Add biological context to methylation–expression coupling

A simple overlap between differentially methylated genes and DEGs is biologically incomplete. The result object should retain CpG identifier and position, regulatory context, CpG-to-gene mapping source, methylation and expression directions, sample-level correlation where paired data exist, distance to the transcription start site, and optional chromatin or enhancer evidence.

<span style="color:#d17a00">Promoter hypermethylation accompanied by reduced expression should not be conflated with gene-body or enhancer methylation patterns. ArthOMix should report these categories separately</span> — <span style="color:#d17a00">**Partly.** Promoter / gene-body / enhancer context and TSS distance are shown on the Biomarker Card; sample-level correlation is used when paired matrices exist. Not independently verified at the base Cross-omics integration table.</span>

## 6. Replace unsupported causal terminology

<span style="color:red">Cross-sectional methylation and expression associations can generate regulatory hypotheses, but they do not ordinarily demonstrate causality. Rename "Dysregulation and causation" as "Regulatory concordance and causal-hypothesis prioritization."</span> <span style="color:red">Corrected — no such label exists; causal wording only in MR / colocalisation.</span>

## 7. Protect diagnostic modelling from information leakage

This is the highest-risk analytical component. Normalization, differential testing, LASSO/Boruta selection, imputation, and hyperparameter tuning must not use the complete dataset before validation.

1. <span style="color:red">Split data at the participant level; when feasible, reserve an independent study or cohort for external validation.</span>
2. <span style="color:red">Within each training fold, perform preprocessing, missing-value handling, feature filtering, feature selection, class balancing, and model tuning.</span>
3. <span style="color:red">Evaluate the locked workflow on the untouched validation fold.</span>
4. <span style="color:red">Use nested or properly repeated cross-validation for performance estimation and tuning.</span>
5. <span style="color:red">Report uncertainty and calibration in addition to discrimination.</span>

<span style="color:red">Outputs should include confidence intervals, sensitivity, specificity, predictive values, precision–recall AUC, calibration plots or metrics, and a clear distinction between exploratory biomarker discovery and a clinically validated diagnostic model.</span> TRIPOD+AI provides an appropriate reporting framework *(the metrics are reported; the framework is not named on screen by choice).*

*Remaining caveat: the bundled/pre-computed demo panel has no held-out split; the app labels that result "not leakage-safe" and substitutes the nested-CV estimate.*

# Software, reproducibility, and usability requirements

- <span style="color:red">Modular Shiny organization with user-interface code separated from tested analysis functions.</span>
- <span style="color:#d17a00">Asynchronous jobs, progress reporting, cancellation, and informative error recovery.</span> *(async + progress + error messages yes; cancellation no)*
- <span style="color:#d17a00">Maximum file size, expected dataset dimensions, memory safeguards, timeouts, and concurrent-user limits.</span> *(3 GB upload cap and HTTP timeouts yes; memory / concurrent-user limits no)*
- <span style="color:#d17a00">Session isolation, authentication options, encrypted transport, retention policy, and secure deletion for private data.</span> *(per-session state, authentication, no disk persistence of uploaded-data results yes; transport / retention are deployment-level)*
- <span style="color:#d17a00">API caching, rate-limit handling, versioned annotation sources, and graceful failure when external services are unavailable.</span> *(caching, retrieval dates, graceful failure yes; retry/back-off no)*
- <span style="color:red">Pinned R/Bioconductor package versions, containerized deployment, controlled random seeds, and automated dependency checks.</span>
- <span style="color:#d17a00">Downloadable parameters, warnings, analysis logs, session information, figures, tables, and reproducible R scripts.</span> *(parameters, provenance JSON, figures, tables yes; per-run R script no)*
- <span style="color:red">Unit, integration, regression, and known-answer tests using synthetic and public benchmark datasets.</span> *(94 test files with planted-effect fixtures)*
- <span style="color:#d17a00">Accessibility, clear help text, sensible defaults, and guardrails that disable invalid method combinations.</span> *(help text, defaults, guardrails yes; accessibility markup minimal)*

<span style="color:red">Every analysis should generate a machine-readable provenance manifest recording input checksums, selected parameters, software and annotation versions, warnings, random seeds, and timestamps.</span>

# Recommended gene-card design

- <span style="color:red">gene identifiers, synonyms, functional summary, and genomic location;</span>
- <span style="color:red">RA or other selected disease associations with source and provenance;</span>
- <span style="color:red">expression effect estimates, uncertainty, and sex interaction;</span>
- <span style="color:red">methylation results separated by promoter, gene body, enhancer, and other contexts;</span>
- <span style="color:red">co-expression or co-methylation module membership;</span>
- <span style="color:red">MOFA/DIABLO contribution and feature-stability evidence;</span>
- blood–synovium or other cross-tissue replication; *(not on the card — separate module)*
- an evidence/stability score whose formula is transparent; *(a documented evidence tier, not a numeric score)*
- <span style="color:#d17a00">external database links, versions, and retrieval dates.</span> *(retrieval date shown on the methylomics card; cached but not displayed on the transcriptomics card)*

<span style="color:#d17a00">External annotations should be cached with a database version or retrieval date.</span>

# Corrections needed in the current draft

*(Manuscript items — none are code. Several refer to the 30 June draft and no longer apply to the current thesis; the ones that still apply to the thesis are: submodule count (thesis says 41, code has 45), "Methyl omics" spelling, duplicated and mis-encoded references, and FDR vs Bonferroni wording.)*

- The manuscript states that there are 24 submodules and then says "out of 25." Reconcile the count and provide one definitive module inventory.
- Transcriptomics is said to contain nine submodules, but only eight are described.
- Methylomics is said to have five submodules, while earlier text implies that additional modules are shared with transcriptomics.
- Meta-analysis is named as a major module but is not described. <span style="color:red">*(Now exists in the app — Study-wise Meta-analysis tab — and must be described.)*</span>
- Gene set enrichment and functional enrichment are not consistently separated or defined.
- The methylation sex-stratified section incorrectly refers to differential gene expression.
- Do not define "top significant genes" simply as the overlap between DEGs and LASSO/Boruta features.
- Do not present Bonferroni correction and FDR as interchangeable defaults.
- Standardize capitalization and spelling: ArthOMix, Shiny, RNA-seq, multi-omics, methylomics, and rheumatoid arthritis.
- Correct unmatched parentheses, duplicated references, incomplete citations, and inconsistent reference formatting; verify all 2025–2026 citations.
- Replace embedded EMF figures with high-resolution PNG or SVG files.

# Recommended development roadmap

| **Step** | **Development priority** | **Status in code** |
| --- | --- | --- |
| 1 | Define the target user and decide whether the core platform is general or RA-specific. | Thesis decision |
| 2 | <span style="color:red">Create formal input schemas and automated validation for each supported data type.</span> | <span style="color:red">Done</span> |
| 3 | <span style="color:#d17a00">Implement a method-selection decision tree based on platform, data scale, sample pairing, tissue, and study design.</span> | <span style="color:#d17a00">Partly — enforced by guards (data type, pairing, confounding) rather than one explicit decision tree</span> |
| 4 | <span style="color:red">Build and validate transcriptomic and methylation single-omics workflows.</span> | <span style="color:red">Done</span> |
| 5 | <span style="color:red">Implement interaction-based sex-aware analysis with power and estimability warnings.</span> | <span style="color:red">Done</span> |
| 6 | <span style="color:red">Separate matched multi-omics integration from unpaired cohort-level synthesis.</span> | <span style="color:red">Done</span> |
| 7 | <span style="color:red">Build leakage-resistant predictive modelling with nested validation and external-cohort testing.</span> | <span style="color:red">Done</span> |
| 8 | <span style="color:#d17a00">Add provenance, versioning, downloadable reports, and reproducible code.</span> | <span style="color:#d17a00">Provenance, versioning, reports done; per-run reproducible script not</span> |
| 9 | <span style="color:#d17a00">Validate each module using simulated data and published benchmark datasets.</span> | <span style="color:#d17a00">Simulated known-effect tests yes; published benchmark validation not</span> |
| 10 | After the analytical core is stable, expand the interface, gene cards, and disease-context features. | Ongoing |

# Proposed core identity

ArthOMix should be distinguished by methodological guidance, not merely by the number of modules. Its strongest identity would be a platform that helps non-specialists choose analyses appropriate to their data while enforcing reproducibility, sex-aware modelling, and rigorous multi-omics integration.

# Selected methodological resources

- **Bioconductor limma RNA-seq workflow:** <https://www.bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html>
- **MOFA2: training and suitability guidance:** <https://bioconductor.org/packages/release/bioc/vignettes/MOFA2/inst/doc/getting_started_R.html>
- **MOFA2 FAQ:** <https://biofam.github.io/MOFA2/faq.html>
- **mixOmics DIABLO documentation:** <https://mixomics.org/mixdiablo/>
- **TRIPOD+AI reporting guideline:** <https://www.equator-network.org/reporting-guidelines/tripod-statement/>
