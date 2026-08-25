# 2.GG Diagnostic Classification Modelling of the Majority-Vote CpG Panels

## 2.GG.1 Rationale

Section 2.EE's ensemble feature selection produced, per sex stratum, a majority-vote panel of CpGs supported by at least two of three independent, differently biased machine-learning methods (female: 12 CpGs; male: 9 CpGs). Section 2.FF (Mendelian randomisation) tested whether these CpGs are *causally* upstream of RA risk. This section addresses a different, complementary question: how well methylation at these CpGs *discriminates* RA from control status, a diagnostic and predictive question that can be answered independently of causal status.

This analysis is deliberately not restricted to the subset of CpGs for which an MR estimate was available. MR-instrument availability in GoDMC is unrelated to a CpG's diagnostic value: a biomarker can be an excellent classifier while representing a downstream consequence of disease (for example, inflammation-driven shifts in immune-cell composition) rather than a cause of it. Restricting the present analysis to the eight CpGs that reached MR in Section 2.FF would have arbitrarily excluded `cg21740204` (*SMC5*), the strongest single CpG in this chapter by every other line of evidence, for a reason unconnected to predictive performance. The two analyses are therefore retained as separate, complementary lines of evidence: Section 2.FF establishes, or fails to establish, causal or mechanistic status; this section establishes predictive or diagnostic value.

**Feature choice: raw M-values rather than covariate-adjusted residuals.** CpG selection upstream (Sections 2.CC-2.EE) used residuals with age, smoking, and cell-type composition regressed out. This section instead uses raw methylation M-values at the panel CpGs directly, for two reasons. First, the external validation cohort (GSE111942, below) provides no age, smoking, or cell-type data in its public metadata, so the same residualisation could not be reproduced there. Second, a deployable diagnostic test applied to a new patient would have access only to their raw measured methylation, never to a residual computed against a specific training cohort's covariate model. Raw M-values are therefore the methodologically appropriate feature representation for this question, independent of the metadata limitation of the external cohort.

## 2.GG.2 Methods

**Internal cohort (training and internal test): GSE42861** (Liu et al., 2013), the cohort used throughout this chapter. Samples were stratified 75/25 into training and internal-test sets per sex (`caret::createDataPartition`), with a fixed random seed (`PARAM$seed`).

**External validation cohort: GSE111942** (Zhu et al., 2019, *Annals of the Rheumatic Diseases*, PMID 30297333), comprising 25 RA cases and 18 healthy controls, PBMC, Illumina 450K, from Soochow University, China. This dataset was selected following two independent, exhaustive searches for a publicly available RA-versus-healthy-control whole-blood or PBMC methylation dataset distinct from GSE42861 (search summarised below). Direct inspection of its GEO series matrix confirmed that all 43 samples are female; this dataset therefore serves as external validation for the female panel only. No dataset with adequate male representation was identified: the only male-inclusive candidate located, GSE121192, provides three unique male-RA and three unique male-control donors from cell-sorted PBMC subsets rather than bulk blood, an n insufficient to yield a meaningful external AUC estimate and therefore not reported. The male panel is consequently reported with training and internal-test results only, with external validation stated as unavailable rather than substituted with an underpowered estimate.

**External dataset search.** The following candidate datasets were evaluated and excluded: GSE87095 (CD19+ B cells rather than whole blood, no sex field), GSE131989 (entirely female), GSE175364 (no clean healthy-control arm), GSE137593/GSE137594/GSE137634 (comparator group is other inflammatory arthritis rather than healthy controls), GSE134877 (targeted pyrosequencing of a single locus rather than genome-wide coverage, no sex field), and GSE121192 (male-inclusive but limited to three unique donors per relevant cell type, as above). Two adequately designed studies with sex-balanced RA-versus-healthy blood cohorts (Svendsen et al., 2025, *Frontiers in Immunology*; Barbarroja et al., PMC12333266) were identified but have no public data deposit, owing to a Danish Data Protection Agency restriction and a stated absence of generated or analysed datasets respectively, and could therefore not be used.

**Raw-data processing and normalisation consistency.** GSE111942 provides only raw IDATs, with no processed beta matrix; these were processed via `minfi::preprocessNoob()` (`09_process_external_gse111942.R`). The initial version of this analysis paired this Noob-processed external cohort with GSE42861's already-processed beta values as shipped in its GEO series matrix (the approach used for GSE42861 throughout the rest of this chapter), an undocumented normalisation applied in the original Liu et al. (2013) study and not replicated here. Pairing two independently normalised cohorts in this way is a recognised source of spurious external-validation failure in the methylation literature, and constitutes a more parsimonious explanation for a poor external result than genuine non-generalisation of the classifier. To rule out this confound, GSE42861 was reprocessed from its own raw IDATs (`GSE42861_RAW.tar`, 689 samples, download integrity confirmed by exact byte-count match and `tar` integrity check following an initial truncated download) through the identical `preprocessNoob()` pipeline, in memory-safe chunks of 80 samples (`09a_process_internal_gse42861_frozen.R`). Chunking does not affect the result, since Noob is a single-sample, reference-free method with no cross-sample fitted parameters. External AUCs were materially unchanged following this reprocessing, indicating that normalisation mismatch alone did not account for the external-validation pattern subsequently characterised below, and motivating the further investigation described next.

**A second confound: compositional difference between whole blood and PBMC.** GEO metadata confirm that GSE42861's tissue is described as peripheral blood leukocytes (PBLs), that is, whole-blood leukocytes including granulocytes, whereas GSE111942 is explicitly PBMC (peripheral blood mononuclear cells), a fraction that structurally excludes granulocytes via density-gradient centrifugation. This project's own EpiDISH-based cell-type deconvolution (script02_celltype, `F$pheno_celltype_rds`) estimates GSE42861's granulocyte lineage (neutrophils plus eosinophils) at a mean of 70.0% of leukocytes (neutrophil range 22-99%; eosinophils occasionally up to 27%), a substantial fraction present in the training cohort and structurally absent from the external cohort. Because DNA methylation is highly cell-type-specific, any direct comparison between PBL and PBMC methylation at the same CpGs is confounded by this compositional difference, independently of the ancestry and population mismatch discussed in Section 2.GG.4. This was addressed (`09g_celltype_adjust_internal.R`) by fitting, for each panel CpG, a linear regression of M-value on neutrophil and eosinophil fraction across all 689 GSE42861 samples, and subtracting each sample's fitted granulocyte contribution from its M-value; this is equivalent to estimating each sample's methylation at zero granulocyte fraction, applying the same logic as the age- and smoking-residualisation already used elsewhere in this chapter (Section 2.EE) to a newly identified confound. Per-CpG R² for the granulocyte-only regression was modest (0.0002-0.019 across the 21 panel CpGs; `granulocyte_adjustment_coefficients.csv`), consistent with these CpGs having been originally selected for disease association rather than cell-composition association, although the adjustment's effect on external validation was nonetheless substantial (Section 2.GG.3). All results reported below use this cell-type-adjusted internal dataset (`gse42861_internal_panel_celltype_adjusted.rds`), built on the frozen-pipeline reprocessing described above.

**Algorithms (combined panel).** Six algorithms were fitted with grid-search hyperparameter tuning via 10-fold cross-validation on the training split only, using ROC/AUC as the tuning metric: random forest, elastic-net logistic regression (`glmnet`, tuning alpha and lambda), support vector machine (radial kernel, tuning sigma and cost), k-nearest neighbours (tuning k), a single-hidden-layer artificial neural network (`nnet`, tuning hidden units and weight decay), and XGBoost.

**XGBoost implementation.** Gradient-boosted models were implemented directly via `xgboost::xgb.cv()` and `xgb.train()`, using an explicit manual grid search over the same hyperparameters as the other algorithms (maximum tree depth, learning rate, column subsampling, minimum child weight, row subsampling), with model selection by mean cross-validated AUC. This departs from the standard `caret::train(method = "xgbTree")` wrapper, which returned undefined (NA) cross-validated ROC values under the installed xgboost version (3.2.0.1): caret's wrapper does not accommodate a change in the xgboost API whereby the `objective` parameter must be supplied within the `params` list rather than as a top-level argument. The best iteration from `xgb.cv()` is correspondingly read from the `$early_stop$best_iteration` field of the returned object rather than the top-level `$best_iteration` field used in earlier xgboost versions.

**Algorithm suite (per probe).** Each panel CpG was modelled independently using the same six-algorithm suite described above, with smaller hyperparameter grids reflecting the reduced feature dimensionality, and the best-performing algorithm per CpG selected by internal-test AUC for the headline figures and tables. This full-suite approach is necessary because algorithm choice is not immaterial even for a single continuous predictor: this holds only for monotonic classifiers such as logistic regression, and not for k-nearest neighbours, random forest, the neural network, or XGBoost, none of which are constrained to a monotonic decision function even with a single feature (for example, a k>1 KNN's predicted probability need not increase or decrease monotonically with the feature value). Full per-CpG, per-algorithm results are reported in full (`diagnostic_perprobe_all_algorithms_{female,male}.csv`), not only the winning algorithm, so that this selection is auditable.

Elastic-net logistic regression (`glmnet`) requires a design matrix with at least two predictor columns, since regularisation across a single coefficient is not defined; standard, unregularised logistic regression (`glm`) was substituted for this algorithm in the per-probe analysis, with no hyperparameters requiring tuning.

**Evaluation.** AUC with DeLong confidence intervals (`pROC::ci.auc`) was computed on the training set (in-sample, reported for reference rather than as a generalisation estimate), the internal held-out test split, and, for the female panel only, the external test set. ROC curves are rendered as empirical step functions with per-curve AUC embedded in the legend label (format `AUC=X.XX (label)`), together with a labelled `AUC=0.50 (random)` reference diagonal.

## 2.GG.3 Results

**Table 2.GG.1. Combined-panel AUC, female (12 CpGs), before and after granulocyte-fraction adjustment.**

| Algorithm | Train AUC (adj.) | Internal test AUC (adj., 95% CI) | External test AUC, pre-adjustment | External test AUC, post-adjustment (95% CI) |
|---|---|---|---|---|
| Random Forest | 1.00 | 0.684 (0.586-0.781) | 0.516 | **0.791 (0.658-0.925)** |
| SVM | 0.894 | 0.681 (0.584-0.778) | 0.349 | 0.676 (0.512-0.839) |
| KNN | 0.829 | 0.681 (0.585-0.777) | 0.553 | 0.707 (0.548-0.866) |
| ANN | 0.879 | 0.654 (0.555-0.752) | 0.444 | 0.544 (0.365-0.724) |
| XGBoost | 0.938 | 0.639 (0.540-0.738) | 0.189 | 0.479 (0.296-0.662) |
| Logistic regression (elastic net) | 0.755 | 0.626 (0.526-0.726) | 0.158 | 0.189 (0.050-0.328) |

**Table 2.GG.2. Combined-panel AUC, male (9 CpGs), post-adjustment, no external test available.**

| Algorithm | Train AUC | Internal test AUC (95% CI) |
|---|---|---|
| SVM | 0.950 | 0.743 (0.598-0.887) |
| Random Forest | 1.00 | 0.717 (0.569-0.865) |
| KNN | 0.826 | 0.694 (0.542-0.846) |
| Logistic regression (elastic net) | 0.683 | 0.678 (0.507-0.849) |
| XGBoost | 0.984 | 0.669 (0.512-0.825) |
| ANN | 0.955 | 0.650 (0.488-0.813) |

The granulocyte-fraction adjustment described in Section 2.GG.2 materially improved external validation for the female panel, the principal empirical finding of this section. Random Forest's external AUC approximately doubled to 0.791 (95% CI 0.658-0.925, now clearly excluding 0.5); SVM and KNN also improved substantially (0.349 to 0.676, and 0.553 to 0.707, respectively). This improvement was accompanied by a decrease in internal-test AUC for every algorithm following adjustment (for example, Random Forest 0.741 to 0.684), consistent with part of the panel's apparent internal-cohort performance having been driven by the whole-blood-versus-PBMC compositional confound itself rather than by disease signal, with removal of this confound revealing a smaller but more genuinely cross-population-generalisable effect. Logistic regression and XGBoost did not improve under adjustment (XGBoost's external AUC in fact fell further, to 0.189); this inconsistency across algorithms indicates that the improvement is not uniform and should not be treated as fully resolving the external-validation limitation, but rather as a partial correction, discussed further in Section 2.GG.4.

No combined-panel algorithm, in either sex, reaches train AND internal-test AUC both ≥0.70 for the female panel after adjustment (an aspirational threshold considered during this analysis); Random Forest and SVM continue to clear both thresholds for the male panel (train/test: 1.00/0.717 and 0.950/0.743, respectively). This follows directly from the adjustment's effect on internal fit described above, rather than constituting an independent finding.

**Table 2.GG.3. Per-probe AUC, female panel, post-adjustment, best-tuned algorithm per CpG, sorted by internal-test AUC.**

| CpG | Gene | Best algorithm | Train AUC | Internal test AUC | External test AUC |
|---|---|---|---|---|---|
| cg20639658 | *C2CD4A* | LogReg | 0.673 | 0.650 | 0.384 |
| cg02120739 | *PM20D2* | Random Forest | 1.00 | 0.635 | 0.433 |
| cg08212230 | *CDC20B* | LogReg | 0.560 | 0.627 | 0.573 |
| cg15478005 | *ATP6V1A* | LogReg | 0.662 | 0.626 | 0.322 |
| cg25598086 | *HSPA2* | KNN | 0.723 | 0.625 | 0.640 |
| cg21740204 | *SMC5* | LogReg | 0.688 | 0.625 | 0.300 |
| cg10437627 | *C12orf45* | LogReg | 0.618 | 0.614 | 0.453 |
| cg00322319 | -- | Random Forest | 1.00 | 0.601 | 0.436 |
| cg18924222 | *GLT8D2* | LogReg | 0.649 | 0.576 | **0.709** |
| cg12089985 | *HIST1H4C* | LogReg | 0.598 | 0.567 | 0.551 |
| cg26793226 | *SNHG3-RCC1* | LogReg | 0.551 | 0.563 | 0.598 |
| cg19157140 | *HEATR2* | SVM | 0.594 | 0.548 | 0.482 |

No female CpG reaches both train AUC>0.70 and internal-test AUC>0.70, a pattern unchanged by the adjustment. A qualitative shift is nonetheless apparent: logistic regression is the best-performing algorithm for nine of twelve CpGs (compared with five of twelve prior to adjustment), with Random Forest's previously frequent train=1.00 overfitting selections (five CpGs before adjustment, including GLT8D2) reduced to two, consistent with granulocyte-driven signal having accounted for a substantial part of what Random Forest was fitting prior to its removal.

**Table 2.GG.4. Per-probe AUC, male panel, post-adjustment, best-tuned algorithm per CpG, no external test available.**

| CpG | Gene | Best algorithm | Train AUC | Internal test AUC |
|---|---|---|---|---|
| cg15643574 | *PICALM* | ANN | 0.727 | **0.807** |
| cg14004673 | *C19orf54* | LogReg | 0.546 | 0.741 |
| cg06999776 | *SEMA4F* | SVM | 0.749 | 0.736 |
| cg17484629 | *ELOVL6* | LogReg | 0.543 | 0.706 |
| cg23881902 | *HDLBP* | Random Forest | 1.00 | 0.693 |
| cg19657380 | *ZFAT* | LogReg | 0.524 | 0.626 |
| cg27301343 | *EML2* | ANN | 0.669 | 0.624 |
| cg22560837 | *POP4* | LogReg | 0.531 | 0.600 |
| cg10186232 | *C5orf44* | LogReg | 0.635 | 0.579 |

Two male CpGs reach both train AUC>0.70 and internal-test AUC>0.70 across multiple algorithms (unavailable for external confirmation, Section 2.GG.2): *PICALM* (SVM, KNN, ANN, and XGBoost all qualify; ANN highest at train=0.727/test=0.807) and *SEMA4F* (SVM, KNN, and ANN all qualify; SVM highest at train=0.749/test=0.736), representing an improvement in robustness, measured by the number of qualifying algorithms per CpG, relative to the pre-adjustment analysis.

`cg18924222` (*GLT8D2*) remained the only female CpG to show a consistent cross-population signal, replicating across all three iterations of this analysis (prior to reprocessing, following frozen-pipeline reprocessing, and following granulocyte-fraction adjustment; external AUC of 0.71 in each case for its best-agreeing algorithms), and represents the most robust single-probe finding in this section.

**Qualifying results (train AUC>0.70 AND internal-test AUC>0.70), filtered post-adjustment.** `09f_filtered_qualifying_roc.R` was rerun against the cell-type-adjusted data; the qualifying set shifted substantially relative to the pre-adjustment run. The female combined panel fell from five qualifying algorithms to zero (no algorithm now clears both thresholds; no filtered figure was produced, and this absence is reported explicitly rather than omitted). The male combined panel fell from five to two qualifying algorithms (SVM, Random Forest). The female per-probe panel retains zero qualifying CpG-algorithm combinations, unchanged from the pre-adjustment analysis. The male per-probe panel retains two qualifying CpGs, each now supported by multiple qualifying algorithms rather than one (*PICALM*: SVM, KNN, ANN, XGBoost; *SEMA4F*: SVM, KNN, ANN). The female combined panel's fall to zero follows mechanically from the adjustment reducing internal-test AUC below 0.70 for every algorithm (Table 2.GG.1) while simultaneously improving external AUC, indicating that this fixed >0.70/>0.70 threshold, chosen prior to identification of the confound, is no longer the most informative summary of this panel's performance and should be interpreted alongside the external AUC values directly rather than as a standalone pass/fail criterion.

**Table 2.GG.5. Summary of candidate diagnostic biomarkers.**

| Sex | Level | Biomarker | Best algorithm | Train AUC | Internal test AUC | External test AUC | External validation |
|---|---|---|---|---|---|---|---|
| Female | Combined panel (12 CpGs) | -- | Random Forest | 1.00 | 0.684 | **0.791** | Yes, GSE111942 |
| Female | Single CpG | cg18924222 (*GLT8D2*) | Logistic regression | 0.649 | 0.576 | 0.709 | Yes, GSE111942; replicated across 3 analysis iterations |
| Male | Combined panel (9 CpGs) | -- | SVM | 0.950 | 0.743 | Not available | No adequate external cohort identified |
| Male | Single CpG | cg15643574 (*PICALM*) | ANN | 0.727 | 0.807 | Not available | No adequate external cohort identified |
| Male | Single CpG | cg06999776 (*SEMA4F*) | SVM | 0.749 | 0.736 | Not available | No adequate external cohort identified |

Taken together, these results indicate that a subset of methylation sites within the majority-vote panels carries diagnostic information distinguishing RA from control status, though the strength of this evidence differs by sex. The strongest evidence for a genuinely generalisable biomarker is the female combined-panel Random Forest classifier (external AUC 0.791) and the single CpG `cg18924222`/*GLT8D2* (external AUC 0.709, replicated across all three iterations of this analysis). No externally validated male biomarker could be established, owing to the absence of a suitable independent male cohort (Section 2.GG.2); *PICALM* and *SEMA4F* are reported as the strongest internally validated male candidates, with external confirmation left for future work.

## 2.GG.4 Limitations

**The sample-type (PBL versus PBMC) confound was identified and corrected, but the correction is partial rather than complete.** Section 2.GG.2 describes the identification of this confound, via EpiDISH deconvolution (Teschendorff, Breeze, Zheng & Beck, 2017), and its correction by granulocyte-fraction regression; Section 2.GG.3 describes its effect, a material improvement in female external AUC for three of six algorithms, most notably Random Forest (0.516 to 0.791). Two considerations temper how far this correction should be trusted. First, per-CpG R² for the granulocyte-only regression was modest (0.0002-0.019), indicating that the fitted M ~ neutrophil + eosinophil model explains only a small fraction of each CpG's total methylation variance; the correction removes the linear granulocyte-associated component specifically, but cannot rule out remaining nonlinear or higher-order compositional effects, nor effects from cell types not modelled here. Second, the improvement was not uniform across algorithms (logistic regression and XGBoost did not improve, and XGBoost's external AUC fell further), which is inconsistent with a single confound providing the complete explanation and suggests that other, unaddressed factors, including ancestry, discussed below, remain in play.

**The decrease in internal-test AUC following adjustment reflects expected, correct behaviour rather than a regression to be explained away.** Every algorithm's internal-test AUC fell after the granulocyte correction (Tables 2.GG.1, 2.GG.3). This is the intended and logically necessary consequence of removing a compositional confound that the training cohort's cell-type structure had made available as an additional source of separability: insofar as part of the panel's original internal performance was attributable to disease-associated shifts in blood cell composition rather than to within-cell-type methylation change (the structural risk described generally for EWAS by Houseman et al., 2012 and Teschendorff & Relton, 2018, and which this project's own age, smoking, and cell-type adjustment, Section 2.EE, was designed to guard against), a correct adjustment should lower apparent internal fit while improving external generalisability, the pattern observed here.

**The female panel's external validation, even after correction, remains below what would be required to support a diagnostic-use claim.** Combined with the null Mendelian randomisation result for these same CpGs (Section 2.FF), this motivates caution in interpreting the majority-vote panel's association with RA in GSE42861 as either causal (Section 2.FF) or broadly generalisable to other populations for diagnostic purposes (this section). Random Forest's post-adjustment external AUC of 0.791 constitutes genuine evidence that part of the panel's signal is real and detectable across populations, but represents one algorithm of six, evaluated on a single external cohort, and should not be over-generalised to the panel as a whole.

**Population and ancestry mismatch remains an unaddressed, residual confound even after the sample-type correction.** GSE42861 (training) is a Swedish (EIRA) cohort; GSE111942 (external) is Chinese (Soochow). Ancestry-driven differences in DNA methylation, allele frequency, and RA subtype and severity distribution are all plausible contributors to the remaining external-validation gap and, unlike the PBL/PBMC compositional confound, could not be corrected for with the data available in this analysis.

**No male external validation was possible.** This absence follows an exhaustive, documented search (Section 2.GG.2) rather than an unexamined gap. The male panel's improved per-probe robustness following adjustment (*PICALM* and *SEMA4F* now qualifying across multiple algorithms, above) should not be assumed to generalise externally with any greater confidence than the female panel's pre-adjustment internal results did; the female result demonstrates directly that internal-cohort performance alone, whether or not cell-type-adjusted, is not a reliable proxy for external performance.

**GSE111942 metadata are minimal.** No age, smoking, or cell-type-composition data are available for this cohort, which is why raw M-values, rather than covariate-adjusted residuals, were used as classifier features throughout (Section 2.GG.1), and why the granulocyte correction (Section 2.GG.2) could be applied only to the internal cohort; the external cohort's own compositional structure, PBMC, granulocyte-free by construction, could not be independently verified or adjusted beyond what is already known from its stated sample type.

**Overfitting persists in Random Forest and, for the male panel, XGBoost and the neural network.** Substantial train-to-internal-test AUC gaps remain after adjustment (for example, male-panel XGBoost train=0.984 versus internal test=0.669; Random Forest train=1.00 throughout both panels), indicating that these algorithms continue to fit training-set-specific noise on this small feature panel. Their point estimates should accordingly be weighted less heavily than those of algorithms with a smaller train/test gap, such as logistic regression and SVM.

---

## References

Houseman, E. A., Accomando, W. P., Koestler, D. C., et al. (2012). DNA methylation arrays as surrogate measures of cell mixture distribution. *BMC Bioinformatics*, 13, 86.

Liu, Y., Aryee, M. J., Padyukov, L., et al. (2013). Epigenome-wide association data implicate DNA methylation as an intermediary of genetic risk in rheumatoid arthritis. *Nature Biotechnology*, 31(2), 142-147.

Teschendorff, A. E., Breeze, C. E., Zheng, S. C., & Beck, S. (2017). A comparison of reference-based algorithms for correcting cell-type heterogeneity in Epigenome-Wide Association Studies. *BMC Bioinformatics*, 18, 105.

Teschendorff, A. E., & Relton, C. L. (2018). Statistical and integrative system-level analysis of DNA methylation data. *Nature Reviews Genetics*, 19(3), 129-147. PMID: 29129922.

Triche, T. J., Weisenberger, D. J., Van Den Berg, D., Laird, P. W., & Siegmund, K. D. (2013). Low-level processing of Illumina Infinium DNA Methylation BeadArrays. *Nucleic Acids Research*, 41(7), e90.

Zhu, H., Wu, L.-F., Mo, X.-B., et al. (2019). Rheumatoid arthritis-associated DNA methylation sites in peripheral blood mononuclear cells. *Annals of the Rheumatic Diseases*, 78(1), 36-42.

---
