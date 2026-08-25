# 2.EE Ensemble Machine-Learning Feature Selection

## 2.EE.1 Rationale

Convergent module-membership and differential-methylation-region (DMR) evidence, as established in Section 2.DD, define a candidate pool of CpGs for each stratum (female: 1,490 CpGs; male: 557 CpGs), but do not themselves constitute a feature-selection procedure: no ranking, model, or inclusion threshold was applied in arriving at that overlap (Section 2.DD.1). Distinguishing which of these candidates carry the strongest association with rheumatoid arthritis (RA) status, from those retained only by the breadth of the two contributing criteria, requires a dedicated selection step. Three independent, methodologically distinct machine-learning algorithms were applied to each stratum's full candidate set for this purpose, and CpGs selected concordantly across methods were taken forward as a consensus panel.

The three algorithms selected, LASSO logistic regression, Boruta, and support vector machine recursive feature elimination (SVM-RFE), differ in their underlying inductive bias: LASSO is a regularised linear model that penalises the number of nonzero coefficients directly; Boruta is a random-forest-based all-relevant method that tests each feature's importance against randomised shadow copies of itself; and SVM-RFE is a margin-based discriminative method that iteratively removes the least-influential feature from a support vector classifier. A CpG selected by all three therefore reflects agreement across substantially different modelling assumptions, rather than an artefact particular to any one method. The three methods were run independently on the same candidate pool and outcome, and their outputs combined only at the final stage by set intersection, rather than in a sequential procedure in which one method's output constrains the next. This parallel-then-intersect design follows the approach of Su et al. (2025), who applied LASSO, Boruta, and SVM-RFE independently to a candidate gene set and derived a final biomarker panel from the intersection of the three outputs.

The same three-method combination is used in the companion transcriptomic analysis of this thesis and is applied identically here for consistency across omics layers: presenting different combination rules for different data types, without a stated justification, would weaken the overall analytical framework, even where, as in the male stratum (Section 2.EE.3), it yields a sparser consensus than an alternative combination might. Random forest recursive feature elimination was evaluated as a candidate substitute for SVM-RFE as the third method; this alternative produced an empty three-way consensus in both strata and was not adopted, in favour of retaining consistency with the transcriptomic chapter's combination rule. Random forest is retained in this analysis in a diagnostic capacity only, consistent with its role in the transcriptomic chapter's reference figures: fitted once on the final selected panel to characterise its out-of-bag classification error and variable importance, rather than contributing to which CpGs are selected (Section 2.EE.2, Method 4).

The panel derived here is intended as a candidate CpG (and nearby-gene) list for a subsequent cross-omics analysis, in which candidates will be screened against methylation-QTL (mQTL) and expression-QTL (eQTL) reference resources (e.g. GoDMC, ARIES) to identify loci with genetic instruments suitable for Mendelian randomisation, testing whether methylation, or co-regulated gene expression, at these loci is causally related to RA. For this reason, more than the single strictest tier of evidence is reported (Section 2.EE.3): a subsequent instrument-availability screen will itself remove candidates lacking usable genetic instruments, and an unduly narrow starting panel risks leaving no candidates to test, particularly in the male stratum.

## 2.EE.2 Methods

**Input data.** For each stratum, the feature matrix comprised covariate-adjusted methylation residuals (age, smoking status, and estimated cell-type composition regressed out by `limma::lmFit()`, with disease group excluded from the regression), restricted to that stratum's Section 2.DD candidate CpGs. Residuals were computed following the same specification as the null model of Section 2.CC.2, applied here to the smaller candidate set (1,490 and 557 CpGs, respectively, rather than the full 20,000-probe network). Class balance was verified before modelling and found approximately even in both strata (female: 238 Control, 253 RA; male: 95 Control, 101 RA); no class-weighting or resampling correction was applied.

**Method 1: LASSO** (Tibshirani, 1996; implemented via `glmnet`, Friedman et al., 2010). L1-penalised logistic regression was fitted by 10-fold cross-validation, repeated 100 times with independent fold assignments, retaining the coefficients at `lambda.1se`, the more conservative of glmnet's two standard rules, in each repeat. A CpG was retained if its coefficient was nonzero in at least 50% of the 100 repeats. This selection-frequency criterion provides robustness against the fold assignment of any single cross-validation split; it is distinct from, and should not be confused with, the formal stability-selection procedure of Meinshausen and Bühlmann (2010), which uses a specific complementary-pairs subsampling scheme with associated error-control guarantees and was not applied here.

**Method 2: Boruta** (Kursa & Rudnicki, 2010). Boruta is an all-relevant feature-selection wrapper built on random-forest variable importance: at each iteration, every feature's importance is compared against the maximum importance obtained by randomly permuted ("shadow") copies of all features, and features are confirmed, rejected, or left tentative according to a running binomial test across iterations. The algorithm was run for a maximum of 500 iterations (above the package default of 100, to obtain a more stable final classification), and `TentativeRoughFix()` was applied to resolve any remaining tentative features to a final Confirmed or Rejected status; only Confirmed features were carried forward.

**Method 3: SVM-RFE** (Guyon et al., 2002), implemented via `caret::rfe()` (Kuhn, 2008) with a linear-kernel support vector classifier, predictors centred and scaled, and 5-fold cross-validation repeated three times over a candidate subset-size grid ranging from 5 CpGs to the full candidate set. Cross-validated classification accuracy was used as the sole performance criterion for both the inner tuning of the classifier's cost parameter and the outer selection of subset size, so that a single consistent metric governed both stages.

Subset size was selected using a lenient tolerance rule rather than the single best-performing size. The cross-validated accuracy curve was maximal at the smallest tested subset size (5 CpGs) in both strata and declined thereafter, with several neighbouring sizes performing within approximately one cross-validated standard deviation of this maximum (Section 2.EE.3). A custom `selectSize` function, `pickSizeLenient`, was therefore used in place of `caret`'s default (`pickSizeBest`): this function selects the largest subset size whose cross-validated accuracy remains within a specified tolerance, set here to 5%, of the best-observed accuracy. This is the converse of `caret`'s built-in `pickSizeTolerance`, which selects the smallest such size in the interest of parsimony; the converse rule was adopted here because a larger, more inclusive candidate panel better serves the downstream mQTL/eQTL screening this analysis is intended to feed (Section 2.EE.1), and because the near-equivalence of accuracy across neighbouring sizes makes the single smallest-tested size an unreliable basis for exclusion.

**Method 4 (diagnostic, not a selection method): Random forest** (Breiman, 2001; `randomForest` package). A single random forest (500 trees, default `mtry`) was fitted to the majority-vote panel defined in Section 2.EE.3, rather than to the full candidate pool, and did not contribute to feature selection. Its role was diagnostic only: the resulting out-of-bag error and mean-decrease-in-Gini variable importance characterise the separability of RA from Control within the already-selected panel, in the same diagnostic role Random Forest occupies in the companion transcriptomic chapter.

**Reporting convention.** Feature-selection tables and figures report CpG probe identifiers rather than gene symbols, since a subset of candidate probes have no annotated gene in `ChAMPdata::probe.features` (e.g. intergenic probes); reporting probe identifiers throughout avoids mixing annotated and unannotated entries within the same table or figure axis. Gene names are given only in the accompanying narrative, where a specific CpG's nearest annotated gene is discussed.

**Consensus.** Three tiers of cross-method agreement are reported: union, comprising CpGs selected by at least one of the three methods; majority vote, comprising CpGs selected by at least two of the three; and strict (three-way) intersection, comprising CpGs selected by all three. The three-way intersection is treated as the primary result, consistent with the reporting convention of the companion transcriptomic chapter. Because this tier is empty in the male stratum (Section 2.EE.3), the majority-vote tier is additionally reported as the practical candidate set for downstream use in that stratum, since an empty panel cannot be carried forward to instrument screening.

**Reproducibility.** All analyses used a fixed random seed (`PARAM$seed = 42`, this thesis's standing seed); each of the 100 LASSO repeats used an additional offset seed (`PARAM$seed + i`) to ensure a distinct, reproducible fold assignment per repeat. All computations are implemented in `07_feature_selection_female.R` and `07b_feature_selection_male.R`, which are self-contained and rerunnable from the Section 2.DD candidate-CpG tables and the project's standard input files.

---

## 2.EE.3 Results

**Table 2.EE.1.** Features selected per method, and by the three agreement tiers (Section 2.EE.2), by stratum.

| Stratum | Candidate pool (Section 2.DD) | LASSO | Boruta | SVM-RFE (lenient size) | Union (≥1/3) | Majority (≥2/3) | Strict (3/3) |
|---|---|---|---|---|---|---|---|
| Female | 1,490 | 6 | 37 | 40 | 70 | **12** | **1** |
| Male | 557 | 1 | 27 | 25 | 44 | **9** | **0** |

In the female stratum, the majority-vote panel comprised 12 CpGs (Table 2.EE.2). `cg21740204` (*SMC5*) was selected by all three methods. A further five CpGs were selected jointly by LASSO and SVM-RFE: `cg10437627` (*C12orf45*), `cg26793226` (*SNHG3-RCC1*), `cg15478005` (*ATP6V1A*), `cg12089985` (*HIST1H4C*), and `cg08212230` (*CDC20B*). The remaining six were selected jointly by Boruta and SVM-RFE: `cg00322319` (unannotated), `cg02120739` (*PM20D2*), `cg19157140` (*HEATR2*), `cg25598086` (*HSPA2*), `cg20639658` (*C2CD4A*), and `cg18924222` (*GLT8D2*).

**Table 2.EE.2.** Female majority-vote panel (Majority, ≥2/3 methods).

| CpG | Gene | LASSO | Boruta | SVM-RFE | Votes |
|---|---|---|---|---|---|
| cg21740204 | *SMC5* | Yes | Yes | Yes | 3 |
| cg10437627 | *C12orf45* | Yes | No | Yes | 2 |
| cg26793226 | *SNHG3-RCC1* | Yes | No | Yes | 2 |
| cg15478005 | *ATP6V1A* | Yes | No | Yes | 2 |
| cg12089985 | *HIST1H4C* | Yes | No | Yes | 2 |
| cg08212230 | *CDC20B* | Yes | No | Yes | 2 |
| cg00322319 | unannotated | No | Yes | Yes | 2 |
| cg02120739 | *PM20D2* | No | Yes | Yes | 2 |
| cg19157140 | *HEATR2* | No | Yes | Yes | 2 |
| cg25598086 | *HSPA2* | No | Yes | Yes | 2 |
| cg20639658 | *C2CD4A* | No | Yes | Yes | 2 |
| cg18924222 | *GLT8D2* | No | Yes | Yes | 2 |

The single three-way consensus CpG, `cg21740204` (*SMC5*), was independently supported by four lines of evidence: identification as one of 18 genome-wide-significant female DMPs in the single-CpG differential methylation analysis (Section 2.BB); membership of the `black` WGCNA module, the most strongly RA-associated female co-methylation module (Section 2.CC.7); location within a genome-wide-significant DMR, by construction of the Section 2.DD candidate pool; and independent selection by all three feature-selection methods applied here. Within the wider 12-CpG majority-vote panel, this CpG ranked fifth by Gini importance in the Random Forest diagnostic (Section 2.EE.2, Method 4): a strong contributor to panel separability, though not the highest-ranked feature once the full panel is considered.

In the male stratum, the majority-vote panel comprised 9 CpGs (Table 2.EE.3), with no CpG reaching three-way agreement. `cg14004673` (*C19orf54*) was selected jointly by LASSO and SVM-RFE. The remaining eight CpGs were selected jointly by Boruta and SVM-RFE: `cg15643574` (*PICALM*), `cg22560837` (*POP4*), `cg23881902` (*HDLBP*), `cg10186232` (*C5orf44*), `cg19657380` (*ZFAT*), `cg06999776` (*SEMA4F*), `cg27301343` (*EML2*), and `cg17484629` (*ELOVL6*). `cg17484629`/*ELOVL6* also appears among this chapter's female module-level findings (Section 2.CC.7), although not within the female majority-vote panel itself.

**Table 2.EE.3.** Male majority-vote panel (Majority, ≥2/3 methods).

| CpG | Gene | LASSO | Boruta | SVM-RFE | Votes |
|---|---|---|---|---|---|
| cg14004673 | *C19orf54* | Yes | No | Yes | 2 |
| cg15643574 | *PICALM* | No | Yes | Yes | 2 |
| cg22560837 | *POP4* | No | Yes | Yes | 2 |
| cg23881902 | *HDLBP* | No | Yes | Yes | 2 |
| cg10186232 | *C5orf44* | No | Yes | Yes | 2 |
| cg19657380 | *ZFAT* | No | Yes | Yes | 2 |
| cg06999776 | *SEMA4F* | No | Yes | Yes | 2 |
| cg27301343 | *EML2* | No | Yes | Yes | 2 |
| cg17484629 | *ELOVL6* | No | Yes | Yes | 2 |

The absence of a three-way consensus CpG in the male stratum is consistent with the comparatively weaker male signal observed throughout this chapter, reflected in fewer significant DMPs, fewer significant DMRs, and the absence of a genome-wide-significant WGCNA module. This result is reported without substitution of a weaker criterion, as it is itself informative of the male stratum's more limited statistical power rather than indicative of a methodological failure; the majority-vote panel is reported as the practical male candidate set for downstream use, since a three-way requirement would leave no male candidates.

For downstream Mendelian randomisation screening, the majority-vote tier (Tables 2.EE.2 and 2.EE.3) is recommended as the primary screening set, with the three-way-consensus CpG it contains (`cg21740204`/*SMC5*, female stratum) prioritised should it retain a usable genetic instrument. The full union set (70 CpGs female, 44 CpGs male; `ensemble_votes_{female,male}.csv`) is retained as a secondary, more exploratory list should the majority-vote tier prove insufficiently covered by available instruments; CpGs supported by only one method within this wider set carry correspondingly weaker evidence.

Cross-validated accuracy for SVM-RFE was maximal at the smallest tested subset size in both strata (5 CpGs; female: 0.670, male: 0.602) and declined as subset size increased (female: falling to 0.378 at the full 1,490 CpGs; male: falling to a minimum of 0.411 near 300 CpGs before partially recovering to 0.498 at the full 557). Applying the lenient tolerance rule described in Section 2.EE.2 (5% tolerance) selected the largest subset within 5% of this maximum: 40 CpGs in the female stratum (accuracy 0.643) and 25 CpGs in the male stratum (accuracy 0.577). The decline in accuracy with increasing subset size is consistent with a low proportion of truly predictive CpGs among the wider candidate pool, and does not by itself indicate that the smaller subsets generalise well in an absolute sense, only that they cross-validated more favourably than larger subsets within this dataset.

The Random Forest diagnostic (Method 4), fitted to the 12-CpG female and 9-CpG male majority-vote panels, produced an out-of-bag error of 0.326 (female, `mtry` = 3) and 0.311 (male, `mtry` = 3), both well below the approximately 0.50 error expected under random classification for a near-balanced two-class problem, indicating that both panels carry genuine discriminative signal. By mean-decrease-in-Gini importance, the highest-ranked CpG in each panel was `cg00322319` (unannotated, female) and `cg22560837` (*POP4*, male); the female three-way consensus CpG, `cg21740204`/*SMC5*, ranked fifth of twelve.

The cross-validated accuracy values reported for SVM-RFE, and the selection frequencies reported for LASSO and Boruta, are model-selection criteria computed on the same data used to determine subset size and feature membership. They do not constitute an unbiased estimate of classification performance on independent data, which would require a held-out test set excluded from every selection step, and should not be interpreted as evidence that, for example, the female 5-CpG SVM-RFE subset would achieve 67% classification accuracy in an independent sample; they indicate only that, within this dataset's cross-validation, the smaller subsets cross-validated more favourably than larger ones.

<!-- Provenance note, not thesis text: rewritten 2026-08-03 to formal PhD-thesis register at the user's request ("write like a thesis PhD") -- the prior version narrated the analysis's development history in the main prose (an SVM-RFE inner/outer metric mismatch caught and fixed on first run; a log()/logmsg() name-shadowing bug in the LASSO diagnostic plotting code; the RF-as-third-selector detour) using informal bolded lead-ins ("A real error caught and fixed...", "A second real error caught...", "A caveat stated plainly, not implied:") and first-person development narrative. That register has been replaced throughout with a formal, declarative statement of the final adopted procedure. Substantive content preserved: the RF-as-third-selector detour and its rejection is retained as a formal methodological justification in 2.EE.1 (result: empty three-way consensus in both strata); the SVM-RFE metric issue is preserved implicitly as a positive methodological statement in Method 3 ("a single consistent metric governed both stages") rather than as a narrated bug-and-fix. Content removed from the thesis body as implementation-only detail with no bearing on any reported result: the log()/logmsg() base-R shadowing bug in the diagnostic-plotting code (affected only a representative LASSO shrinkage-path plot's internal log-lambda computation, not any selection result). The full bug narrative for both issues remains in the corresponding script comments (07_feature_selection_female.R Section 5, 07c_diagnostic_panels_and_venn.R) and in this project's conversation history if it needs to be recovered for an appendix or supplementary methods note. All numbers, gene annotations, table contents, and citations are unchanged from the previously finalised and independently re-verified version (see prior provenance entries below) -- this pass changed prose register and structure only, not content. Em-dash removal, "+"-to-"and" conversion, and in-text citation formatting (et al. for 3+ authors, "and" in narrative citations, "&" in parenthetical citations) from the preceding pass are preserved throughout this rewrite.

Prior provenance (content/verification history, retained for continuity): LASSO/Boruta/SVM-RFE selections and the LASSO diagnostic panels come from 07_feature_selection_female.R (rerun completed 2026-08-03 12:13:25-12:20:57) and 07b_feature_selection_male.R (rerun completed 2026-08-03 12:13:25-12:15:16), which added the pickSizeLenient custom selectSize function (5% tolerance, largest-within-tolerance) at the user's request to replace caret's default pickSizeBest, which had picked only the smallest tested size (5 CpGs) in both strata; LASSO and Boruta were unaffected by this change. The Venn diagram, SVM-RFE accuracy/error panels, and the Random Forest diagnostic were rebuilt by 07c_diagnostic_panels_and_venn.R (rerun completed 2026-08-03 12:47:52-12:48:20) from the above scripts' saved selections, without recomputation. The three-way Venn is rendered with eulerr::venn() rather than VennDiagram::draw.triple.venn, adopted because the latter's fixed non-proportional template misplaced the smallest set's (LASSO) category label for this topology (LASSO almost entirely nested inside SVM-RFE, several pairwise overlaps at or near zero); eulerr::venn() places every label and count inside its own region automatically, and required the polylabelr package as a transitive dependency, installed alongside it. No figure in this section carries a rendered title, and all use the same canvas dimensions (FIG$w x FIG$h), per this project's figure-standards convention (captions are written separately, outside the image). All counts, gene annotations, and cross-references were independently re-derived from the current saved tables/figures and matched exactly as of the most recent verification pass. Figures: figures/{venn_feature_selection,svmrfe_accuracy,svmrfe_error,boruta_importance,lasso_shrinkage_path,lasso_coefficient_direction,rf_oob_error,rf_importance}_{female,male}.png. Tables: tables/{lasso_selection_frequency,boruta_decisions,svmrfe_selected,svmrfe_size_performance,lasso_coefficient_direction,rf_importance,consensus_features,ensemble_votes,feature_selection_summary}_{female,male}.csv. Table number 2.EE.1 is a placeholder pending final chapter numbering. -->

---

## References

Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5-32.

Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. *Journal of Statistical Software*, 33(1), 1-22.

Guyon, I., Weston, J., Barnhill, S., & Vapnik, V. (2002). Gene selection for cancer classification using support vector machines. *Machine Learning*, 46, 389-422.

Kuhn, M. (2008). Building predictive models in R using the caret package. *Journal of Statistical Software*, 28(5), 1-26.

Kursa, M. B., & Rudnicki, W. R. (2010). Feature selection with the Boruta package. *Journal of Statistical Software*, 36(11), 1-13.

Meinshausen, N., & Bühlmann, P. (2010). Stability selection. *Journal of the Royal Statistical Society: Series B*, 72(4), 417-473.

Su, W., Huang, H., Ruan, Q., Liu, S., Dan, J., Zhao, Y., Zhang, H., & Huo, Q. (2025). Identification and validation of biomarkers associated with cellular senescence and demethylation in acute myocardial infarction. *Scientific Reports*, 15, 41385.

Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. *Journal of the Royal Statistical Society: Series B*, 58(1), 267-288.

---
