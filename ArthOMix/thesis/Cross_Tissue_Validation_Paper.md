# Cross-Tissue Validation: a sex-stratified module for evaluating blood-derived biomarker panels in an independent tissue cohort within ArthOMix

**Module:** `mod_crosstissue.R` (ArthOMix, Transcriptomics area, Section 2.11, "Validation" group)
Format follows the Implementation/Features structure of xOmicsShiny (Wang *et al.*, *Bioinformatics Advances*, 2025), adapted here to document one module of the ArthOMix web application rather than the whole application.

---

## Abstract

**Summary:** Biomarker panels selected in one tissue rarely stay informative, or even point in the same direction, when checked in the tissue where disease actually occurs. Cross-Tissue Validation is a module of the ArthOMix Shiny application that lets a user carry a blood-derived, sex-stratified gene panel into an independent synovial-tissue cohort and ask three questions of it, per sex: does each gene's effect size and significance in synovium support it; does its direction of association agree with blood; and does the panel, refit from scratch (never transported as a locked model) in synovium, still discriminate disease from control. The module ships with a bundled, pre-processed synovial RNA-seq cohort (GSE89408, *n* = 180) and, as of this update, also accepts a user-uploaded raw-count validation cohort with its own sample metadata, running the identical sex-stratified discovery and four-classifier evaluation pipeline on either data source. Results are reported with an apparent (resubstitution) fit reported explicitly as an optimistic upper bound and a pooled out-of-fold cross-validated estimate reported as the headline synovium performance number, never conflated with a locked-coefficient transfer of the blood model.

**Availability and implementation:** Cross-Tissue Validation is implemented in R (Shiny, `shinydashboard`) as `ArthOMix/R/transcriptomics/mod_crosstissue.R` within the ArthOMix application, and is distributed with the application's source code as part of this thesis project. It depends on `edgeR`, `limma`, `glmnet`, `randomForest`, `e1071`, `caret`, and `pROC`. It is not deployed at a public URL; it runs as one tab of the locally/institutionally hosted ArthOMix Shiny app.

---

## 1. Introduction

A gene panel selected from a blood transcriptome carries an implicit assumption: that the biology it captures is not specific to blood. Rheumatoid arthritis (RA), however, is a disease of the synovial joint, and a growing literature on tissue-specific gene regulation means this assumption cannot be taken for granted — a gene differentially expressed in blood may be silent, or even reversed in direction, at the site of disease. Testing that assumption requires an independent tissue cohort, and testing it honestly requires being explicit about what kind of test is actually being run: whether the blood-fitted model itself is being transported unchanged into new data (a strong, locked-coefficient generalizability test), or whether only the *identity* of the selected genes is being carried across, with all model parameters re-estimated inside the new tissue (a weaker, but still informative, portability test) (Justice, Covinsky and Berlin, 1999; Steyerberg and Harrell, 2016).

Cross-Tissue Validation is built around the second, weaker claim, and states this distinction directly in its own interface rather than leaving it to a methods appendix: every classifier fit inside this module is refit from scratch on the validation-tissue data, and only the panel's gene identities are inherited from blood. The module does this sex-stratified throughout, mirroring the sex-stratified design used earlier in the ArthOMix pipeline (Feature Selection, Diagnostic Model), because a panel's synovial behaviour is not guaranteed to be uniform across sexes.

Two extensions distinguish this module from a single-purpose validation script. First, it exposes both an apparent (resubstitution) estimate of panel-classifier discrimination and an out-of-fold, cross-validated pooled estimate, and refuses to let either one stand in for the other — the apparent estimate is labelled an optimistic upper bound at the point it is displayed (Harrell, Lee and Mark, 1996), not merely in a caveat elsewhere. Second — the extension documented in this paper — the module no longer requires the bundled GSE89408 cohort to be the only validation tissue available: a user can upload their own raw-count expression matrix and sample metadata for any independent tissue cohort, sex- and group-labelled, and the module runs the identical discovery-and-classification pipeline on it.

## 2. Methods

### 2.1 Overall architecture

Cross-Tissue Validation is registered as one config/UI/server triple (`mod_crosstissue_config`, `mod_crosstissue_ui`, `mod_crosstissue_server`) inside ArthOMix's transcriptomics module registry, appearing as Section 2.11 of the application, positioned after Sex Interaction Analysis (2.10) and before Cross-Ancestry Validation (2.12). The server function reads two pieces of already-computed session state from the wider application — a live Feature Selection consensus panel (`results$featureselection`) and a live or bundled blood differential-expression result (`results$dge_runs` / a bundled `dge_results.rds`) — but does not read the blood expression matrix itself; the whole point of the module is to operate on a compartment the rest of the application never touches.

A single reactive dispatch point, `val_active()`, resolves to one of two data sources with an identical internal shape (a per-gene differential-expression table, a log-CPM expression matrix, and per-sample sex/group labels): the bundled validation object, or a user-uploaded cohort processed live (§2.5). Every discovery and classification function downstream reads only from this dispatch point, so no analysis code is duplicated between the two data sources.

### 2.2 Sex-stratified gene-level discovery and blood-direction concordance

For a chosen sex, each requested panel gene present in the validation-tissue expression matrix is scored against three quantities: its tissue-level log₂ fold change and Benjamini–Hochberg-adjusted significance, both read from a single, sex-adjusted differential-expression fit shared across both sexes (§2.3); its concordance with blood, defined as agreement in the *sign* of its validation-tissue log₂ fold change and its blood log₂ fold change for the same sex; and a per-gene area under the receiver-operating-characteristic curve (AUC) computed within that sex's own validation-tissue samples, under two orientation conventions — a best-direction convention (AUC ≥ 0.5 by construction, "how much information does this gene carry") and a blood-train-fixed convention (oriented by blood's own direction, so an AUC below 0.5 signals a reversed association out of sample). A gene is flagged a validated cross-tissue biomarker only when all three criteria are jointly met: concordant direction, adjusted significance below a user-set threshold (default 0.05), and best-direction AUC at or above 0.70 (Hosmer, Lemeshow and Sturdivant's conventional "more than weak discrimination" cutoff) — one definition, read by every KPI tile, plot and table in the module, so the three views can never disagree about which genes qualify.

### 2.3 Differential expression in the validation tissue

The per-gene tissue log₂ fold change and adjusted significance used throughout §2.2 come from a single sex-adjusted limma-voom differential-expression fit, computed once per validation dataset and shared across both sexes' discovery views. Raw counts are assembled into an `edgeR::DGEList`; genes of insufficient expression are removed by `filterByExpr`, applied with disease/comparison group as the grouping factor (Chen, Lun and Smyth, 2016); library sizes are normalized by the trimmed mean of M-values (Robinson and Oshlack, 2010); and the design matrix contains both disease group and sex, so that the group coefficient is estimated adjusted for sex. Precision weights are estimated by `voom` (Law *et al.*, 2014), gene-wise linear models are fitted and moderated by empirical Bayes shrinkage (Smyth, 2004; Ritchie *et al.*, 2015), and p-values are corrected across genes by the Benjamini–Hochberg procedure (Benjamini and Hochberg, 1995). For the bundled dataset this fit is precomputed and bundled as part of the shipped validation object; for a user-uploaded cohort (§2.5) the identical pipeline is run live, so the two data sources are scored by the same method rather than two different ones.

### 2.4 Panel-classifier engine: apparent fit and cross-validated estimate

Independently of the per-gene discovery view, the module fits a multivariable panel classifier per sex, four ways: unpenalized logistic regression; elastic-net logistic regression (`glmnet`, alpha and lambda tuned by inner cross-validation; Friedman, Hastie and Tibshirani, 2010; Zou and Hastie, 2005); random forest (`randomForest`, `mtry` tuned by inner cross-validation; Breiman, 2001); and a support-vector machine (`e1071`, cost tuned by inner cross-validation, linear kernel by default; Cortes and Vapnik, 1995). Because the validation-tissue cohort has no natural held-out partition of its own — unlike the blood training cohort, it *is* the held-out compartment — each classifier is fit once on the entire validation-tissue sex-subset, gene-wise z-scored, and its resubstitution AUC is reported explicitly as an optimistic upper bound (Harrell, Lee and Mark, 1996), never as a performance estimate. A second, independent estimate is obtained by refitting each classifier across the folds of an outer, user-configurable cross-validation (stratified by disease status by default, or simple random to match the underlying pipeline's own published script), pooling every out-of-fold prediction into a single receiver-operating-characteristic curve before computing one pooled AUC and a DeLong confidence interval (DeLong, DeLong and Clarke-Pearson, 1988) — this pooled estimate, not the apparent one, is the module's headline synovium performance number. Hyperparameter tuning uses the same standard toolchain throughout (`caret`; Kuhn, 2008).

Only the identity of the panel's genes is inherited from blood; no blood-fitted coefficient is ever transported into a validation-tissue model, and this is stated directly in the module's own interface and in every exported model bundle's metadata (Ambroise and McLachlan, 2002, on the selection-bias risk of treating a fixed, externally derived panel as if its features were free).

### 2.5 User-uploaded validation data

The module accepts a second data source without altering the analysis described in §2.2–2.4: a raw RNA-seq count matrix and an accompanying sample-metadata table, uploaded as CSV or RDS. The user maps three metadata columns — sample identifier, sex, and disease/comparison group — through auto-populated dropdowns, and selects which group value is the reference (e.g. healthy/control) and which is the comparison (e.g. disease) class, the same "map your own columns" interaction already used elsewhere in the application's external-validation workflows. On submission, the module intersects the uploaded expression columns with the metadata's sample identifiers, resolves sex values by prefix match ("F…"/"M…"), restricts samples to the two chosen groups, and runs the identical filterByExpr/TMM/voom/limma/eBayes pipeline described in §2.3 to produce a differential-expression table and log-CPM matrix in the same shape as the bundled validation object. From that point forward — gene-level discovery, concordance, and all four panel classifiers — the uploaded cohort is processed by exactly the same code path as the bundled GSE89408 cohort; no separate analysis branch exists for uploaded data.

## 3. Features

### 3.1 Synovium Discovery & Concordance

Per sex: KPI tiles for the count of validated cross-tissue biomarkers, the fraction of panel genes present in the validation tissue, the fraction direction-concordant with blood, and the median AUC among validated biomarkers; a blood-vs-validation-tissue log₂ fold-change scatter with validated biomarkers highlighted and labelled; a ranked per-gene AUC plot under the user's chosen orientation convention, with reference lines at chance (0.5) and the biomarker cutoff (0.70); and a sortable, filterable per-gene results table (gene, biomarker flag, presence, tissue log₂FC, adjusted P, blood log₂FC, concordance, both AUC conventions), exportable as CSV.

### 3.2 Panel Classifier — Full Fit

Per sex, per classifier: KPI tiles for the apparent AUC, the mean cross-validated AUC (± SD), and the selected hyperparameter; an apparent ROC curve; a cross-validated-AUC-by-fold bar chart; a hyperparameter-search plot (or an explanatory note where no grid applies); a performance table with CSV export; and a downloadable trained-model object bundling the fitted model, its hyperparameters, the gene panel, and a written scoring recipe for applying it to new data. A per-sex comparison table lines up all four classifiers' apparent and cross-validated AUCs side by side.

### 3.3 Panel Classifier — Cross-Validated

Per sex, per classifier: KPI tiles for the pooled out-of-fold AUC with its 95% confidence interval (or an explicit "unavailable" state with a stated reason) and the number of samples contributing an out-of-fold prediction; the pooled ROC curve; and a performance table with CSV export. A per-sex comparison table lines up all four classifiers' pooled cross-validated AUCs.

### 3.4 Cross-Dataset Comparison

Read-only: for each sex, a note on whether a Diagnostic Model blood result is available this session for comparison; a grouped bar chart placing the validation-tissue apparent and pooled cross-validated AUCs alongside blood's own full-fit and cross-validated AUCs, per classifier; and a matching table. The module's own description text states explicitly that this is a side-by-side presentation, not a transfer of the blood model.

### 3.5 Preloaded and user-uploaded validation datasets

A single sidebar control switches the module's data source between the bundled synovial cohort (GSE89408, RA versus histologically normal synovium, *n* = 180: 120 female / 60 male, 28 normal / 152 RA) and an uploaded cohort (raw counts plus metadata, mapped through sample-ID/sex/group columns and a reference/comparison group choice, §2.5). Switching sources changes nothing about §3.1–3.4 beyond which dataset populates them; every plot, table, KPI tile and download in the module works identically against either source, and every exported artifact (trained-model bundle, per-gene table) records which dataset produced it.

## 4. Conclusions

Cross-Tissue Validation gives blood-derived, sex-stratified biomarker panels a single, disciplined route to being checked in the tissue where rheumatoid arthritis actually manifests — reporting tissue-specific effect size, statistical evidence, blood-direction concordance, and both an optimistic apparent classifier fit and a more conservative pooled cross-validated estimate, without ever presenting a re-estimated synovial model as a locked-coefficient transfer of the blood classifier. Extending the module to accept a user-uploaded validation cohort, processed by the identical limma-voom discovery pipeline and identical four-classifier panel engine as the bundled GSE89408 cohort, generalizes this check beyond one fixed synovial dataset without introducing a second analysis method to maintain or reason about.

---

## Availability and implementation

Cross-Tissue Validation is one module of the ArthOMix Shiny application, implemented in R at `ArthOMix/R/transcriptomics/mod_crosstissue.R` and registered via `ArthOMix/R/submodules_registry.R`. It depends on `edgeR` (≥ 4.4.2), `limma` (≥ 3.62.2), `glmnet`, `randomForest`, `e1071`, `caret`, and `pROC`, all already declared dependencies of the wider ArthOMix application. The bundled validation object (`val_synovium.rds`, GSE89408) ships with the application's preloaded data; user-uploaded validation data is processed in-session and is not persisted beyond the running Shiny session. The application is distributed as source code with this thesis project and is not deployed at a public URL.

## References

Ambroise, C. and McLachlan, G.J. (2002) 'Selection bias in gene extraction on the basis of microarray gene-expression data', *Proceedings of the National Academy of Sciences*, 99(10), pp. 6562–6566.

Benjamini, Y. and Hochberg, Y. (1995) 'Controlling the false discovery rate: a practical and powerful approach to multiple testing', *Journal of the Royal Statistical Society: Series B*, 57(1), pp. 289–300.

Breiman, L. (2001) 'Random forests', *Machine Learning*, 45(1), pp. 5–32.

Chen, Y., Lun, A.T.L. and Smyth, G.K. (2016) 'From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline', *F1000Research*, 5, 1438.

Cortes, C. and Vapnik, V. (1995) 'Support-vector networks', *Machine Learning*, 20(3), pp. 273–297.

DeLong, E.R., DeLong, D.M. and Clarke-Pearson, D.L. (1988) 'Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach', *Biometrics*, 44(3), pp. 837–845.

Friedman, J., Hastie, T. and Tibshirani, R. (2010) 'Regularization paths for generalized linear models via coordinate descent', *Journal of Statistical Software*, 33(1), pp. 1–22.

Harrell, F.E., Lee, K.L. and Mark, D.B. (1996) 'Multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors', *Statistics in Medicine*, 15(4), pp. 361–387.

Justice, A.C., Covinsky, K.E. and Berlin, J.A. (1999) 'Assessing the generalizability of prognostic information', *Annals of Internal Medicine*, 130(6), pp. 515–524.

Kuhn, M. (2008) 'Building predictive models in R using the caret package', *Journal of Statistical Software*, 28(5), pp. 1–26.

Law, C.W., Chen, Y., Shi, W. and Smyth, G.K. (2014) 'voom: precision weights unlock linear model analysis tools for RNA-seq read counts', *Genome Biology*, 15, R29.

Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W. and Smyth, G.K. (2015) 'limma powers differential expression analyses for RNA-sequencing and microarray studies', *Nucleic Acids Research*, 43(7), e47.

Robinson, M.D. and Oshlack, A. (2010) 'A scaling normalization method for differential expression analysis of RNA-seq data', *Genome Biology*, 11, R25.

Smyth, G.K. (2004) 'Linear models and empirical Bayes methods for assessing differential expression in microarray experiments', *Statistical Applications in Genetics and Molecular Biology*, 3, Article 3.

Steyerberg, E.W. and Harrell, F.E. (2016) 'Prediction models need appropriate internal, internal-external, and external validation', *Journal of Clinical Epidemiology*, 69, pp. 245–247.

Wang, J. *et al.* (2025) 'xOmicsShiny: an R Shiny application for cross-omics data analysis and pathway mapping', *Bioinformatics Advances*, 5(1), vbaf097.

Zou, H. and Hastie, T. (2005) 'Regularization and variable selection via the elastic net', *Journal of the Royal Statistical Society: Series B*, 67(2), pp. 301–320.

**Note on the dataset citation.** The primary publication for the synovial RNA-sequencing cohort (GSE89408) should be inserted before this document is submitted as part of the thesis; it is not fabricated here.
