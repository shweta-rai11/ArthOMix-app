# Thesis §2.8 — Sex-stratified feature selection (METHODS ONLY)

Written against `scripts/goal2_sex_stratified/12_feature_selection.R` and
`12b_feature_selection_noMHC.R`, with upstream detail from
`scripts/00_shared/03_normalize_batch.R`, `04_apply_holdout.R` and `10_MR.R`.
Software versions from `ENVIRONMENT.txt`. Citation style: author–date (Harvard),
in-text only. Contains no results: no selected-gene counts, no tuned parameter
values, no panel membership.

Section numbering follows the canonical scheme of `METHODS_00_INDEX.md`.

---

## 2.8 Sex-stratified feature selection

### 2.8.1 Genes eligible for feature selection

Feature selection was restricted to genes carrying genetic evidence of a directional effect on
rheumatoid arthritis (RA) liability. A gene was eligible if its *cis*-eQTL-instrumented expression
(Võsa et al., 2021) was associated with RA risk in the two-sample Mendelian randomisation (MR)
analysis against the European-ancestry GWAS of Okada et al. (2014) at a Benjamini–Hochberg false
discovery rate below 0.05, with the correction applied within that sex's stratum (Benjamini and
Hochberg, 1995). Eligibility was independent of the direction of the MR estimate: both
risk-increasing (OR > 1) and protective (OR < 1) genes were carried forward.

The eligible set is therefore the product of two filters of different provenance, and the
distinction governs what the resampling of §2.9.6 is able to correct. The genes *submitted* to MR
were the candidate set of §2.5, formed by intersecting the disease-associated co-expression modules
of §2.4 with the sex-stratified differentially expressed genes of §2.3; both were computed on the
whole training partition and both used the training samples' diagnosis labels. The filter *applied*
to that set was estimated entirely outside this cohort, from external eQTL and genome-wide
association summary statistics, and could not have been influenced by any sample analysed here. The
eligible set is thus an externally filtered internal list rather than an external list, and it is
not described as free of selection bias on the strength of its genetic component alone. The
consequence for the performance estimates, and the reasons the residual optimism is bounded, are
set out in §2.9.2; the sealed holdout and the external cohorts are unaffected, having taken no part
in the differential-expression, network, prioritisation or selection stages.

Because both the eQTLGen exposure data and the Okada outcome GWAS are sex-combined, MR was
estimated once across the union of the female and male candidate sets and partitioned by stratum
afterwards (§2.6); a gene eligible in both strata therefore carries an identical estimate, and any
difference between the two eligible sets arises from upstream gene eligibility and from the
within-stratum correction denominator, not from sex-specific genetic effects. Since the
colocalisation analyses of §2.7 could not establish a shared causal variant for these genes, they
are referred to throughout as *genetically prioritised* rather than as causal genes.

### 2.8.2 Training data, partitioning and preprocessing

Feature selection was performed exclusively on the training partition. The female and male
analyses were run independently: every model was fitted within a single sex, with no exchange of
samples, cross-validation folds or fitted parameters between them.

Expression values from GSE93272 and GSE110169 were placed on the log₂ scale at import and merged
on their common gene set. The merged, pre-normalisation matrix was partitioned 70:30 into a
training set and a sealed internal-validation holdout by stratified random sampling on the
combination of dataset, diagnosis and sex (`createDataPartition`, caret v7.0.1; Kuhn, 2008).
Partitioning preceded every cross-sample operation by design: quantile normalisation learns a
reference distribution shared across samples, and ComBat estimates batch and covariate effects
using the outcome labels, so estimating either on the pooled data would let the holdout influence
the model's input representation and bias performance estimates optimistically (Ambroise and
McLachlan, 2002; Simon et al., 2003; Kaufman, Rosset and Perlich, 2012). Quantile normalisation
(Bolstad et al., 2003), applied with `normalizeBetweenArrays` (limma v3.62.2; Ritchie et al.,
2015), and parametric empirical-Bayes batch correction by ComBat with diagnosis and sex protected
in the model matrix (sva v3.54.0; Johnson, Li and Rabinovic, 2007; Leek et al., 2012) were
therefore estimated on the 70% training partition alone; the corresponding parameters were frozen
and applied to the holdout only at the evaluation stage (§2.9).

For each sex, a design matrix was assembled with samples as rows and that sex's genetically
prioritised genes as columns. The response was diagnosis, coded as RA versus healthy control with
healthy control as the reference level. Three supervised feature-selection algorithms — a
penalised regression, a tree ensemble and a margin-based selector — were then applied within each
sex, the primary hyperparameter of each tuned by 10-fold cross-validation internal to the training
partition. A single global seed (1234) was set before every stochastic step, so that fold
assignment, forest growth and cross-validated tuning are exactly reproducible. All analyses were
run in R v4.4.2.

### 2.8.3 LASSO logistic regression

LASSO logistic regression was fitted with glmnet v5.0 (Friedman, Hastie and Tibshirani, 2010),
specifying `family = "binomial"` and `alpha = 1` for a pure L₁ penalty (Tibshirani, 1996), which
performs embedded selection by shrinking uninformative coefficients to exactly zero. The penalty
parameter λ was tuned over the default glmnet λ sequence by 10-fold cross-validation with binomial
deviance as the loss function (`cv.glmnet`, `nfolds = 10`, `type.measure = "deviance"`). Genes
retaining a non-zero coefficient at λ_min, the value minimising cross-validated deviance,
constituted the LASSO selection; the sparser λ_1se solution was recorded alongside it.

### 2.8.4 Random forest importance

Random forests were grown with randomForest v4.7.1.2 (Breiman, 2001; Liaw and Wiener, 2002) using
1,000 trees, raised from the package default of 500 to stabilise the importance estimates. The
`mtry` parameter — the number of predictors sampled as split candidates at each node — was tuned
by 10-fold cross-validated grid search maximising the area under the receiver operating
characteristic curve (`train`, `method = "rf"`, `metric = "ROC"`, `trControl = trainControl(method
= "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)`; caret v7.0.1; Kuhn,
2008). For a candidate set of p genes the grid comprised 1, 2, ⌊√p⌋, ⌊p/3⌋, ⌊p/2⌋ and p, each
capped at p with duplicate values removed; ties were resolved in favour of the smaller `mtry`. The
final 1,000-tree forest was refitted at the tuned value with `importance = TRUE`. A gene was
selected when its mean decrease in Gini impurity exceeded the arithmetic mean of that statistic
across all candidate genes for that sex, retaining genes that contribute more than the average
candidate to node-impurity reduction.

### 2.8.5 Support vector machine recursive feature elimination

Recursive feature elimination was performed with a linear-kernel support vector machine (Cortes
and Vapnik, 1995) following Guyon et al. (2002), implemented through the e1071 v1.7.17 interface
to LIBSVM (Chang and Lin, 2011). Predictors were scaled internally by `svm`. The cost parameter C
was tuned first, by 10-fold cross-validated grid search (`tune`, `tune.control(sampling = "cross",
cross = 10)`) over C ∈ {0.01, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16}. At the tuned cost, a linear SVM was
trained on the full feature set, the single feature with the smallest squared weight w² was
removed, and the model refitted; this was iterated until one feature remained, producing a
complete ranking from most to least important. Panel size was then selected as the number of
top-ranked features minimising 10-fold cross-validated classification error, evaluated for every
k from 1 to the full candidate set.

### 2.8.6 Consensus panel definition

The primary sex-stratified panel was defined as the intersection of the three selections, so that
every member is independently identified by the penalised regression, the tree ensemble and the
margin-based selector. Three-way agreement was required in place of any single selector because
the male training stratum contains fewer samples than the female stratum relative to the number of
candidate features, a regime in which individual multivariable selectors are unstable. Because
selection used the training partition, no performance estimate is derived from these fits; panel
discrimination was assessed separately by nested cross-validation, on the sealed internal holdout
and in independent external cohorts (§2.9).

The identical procedure — same three algorithms, same tuning grids, same cross-validation scheme,
same seed and same consensus rule — was re-run on the candidate set remaining after exclusion of
the major histocompatibility complex, yielding the secondary MHC-free panel used in the
sensitivity analyses of §2.6 and §2.9.
