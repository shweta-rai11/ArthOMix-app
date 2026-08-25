# 2.9 Machine learning-based diagnostic model development

*Verified against `12_feature_selection.R`, `12b_feature_selection_noMHC.R`,
`14_model_training_nested_cv.R`, `15_model_training_elasticnet.R`,
`16_model_training_final_panel.R`, `16b_model_training_final_panel_noMHC.R`,
`16d_nested_cv_reconciliation.R`, `17_testing_blood_internal_external.R`,
`18_testing_blood_pergene_roc.R`, `19_testing_blood_clinical_utility.R`,
`00_shared/04_apply_holdout.R`, and against the output tables
`mr_fs_summary.csv`, `mr_fs_summary_noMHC.csv`, `NESTED_CV_AUTHORITATIVE.csv`,
`mr_nested_cv_summary.csv`, `PANEL_primary_vs_noMHC_{performance,delong}.csv`
and `renv.lock`.*

---

## Corrections applied in this revision (not for the thesis body)

These are substantive, not stylistic. Each was checked against the code.

1. **The ten-fold cross-validated training AUC must not be called the "honest"
   estimate.** An earlier draft described the training cohort as assessed by "two
   metrics: the apparent AUC and the honest training estimate, represented by the
   10-fold stratified cross-validated out-of-fold AUC". That ten-fold figure
   ([17_testing_blood_internal_external.R:131-149](scripts/goal2_sex_stratified/17_testing_blood_internal_external.R#L131-L149))
   cross-validates *only the logistic model* over a panel that was selected on
   the whole training set. It is the *flat* estimate, precisely the biased
   quantity the nested design exists to displace (Ambroise and McLachlan, 2002;
   Varma and Simon, 2006). The honest training figure is the nested estimate from
   `NESTED_CV_AUTHORITATIVE.csv`. The three quantities are now named apparent,
   flat and nested throughout, and the flat one is never described as honest.
2. **"log₂-transformed to a linear scale" was self-contradictory and the
   condition was omitted.** GSE15573 intensities are placed on the log₂ scale
   *only where they are not already logged*, the test being a maximum above 50
   ([17:84](scripts/goal2_sex_stratified/17_testing_blood_internal_external.R#L84);
   [16b:106](scripts/goal2_sex_stratified/16b_model_training_final_panel_noMHC.R#L106)).
3. **Elastic-net standardisation was misdescribed.** The frozen fold-train
   z-score applies to the consensus/logistic path. In the elastic-net path the
   raw fold-training matrix is passed to `cv.glmnet(..., standardize = TRUE)`, so
   standardisation is internal to glmnet and derived from the fold-training rows
   alone, coefficients being returned on the original scale
   ([16d:116-124, 157-161](scripts/goal2_sex_stratified/16d_nested_cv_reconciliation.R#L116-L161)).
   Leakage-free either way, but the mechanism differs and is now stated.
4. **Zero-variance handling differs between the two standardisation regimes.**
   In-fold, σ = 0 or NA is replaced by unity ([16d:150](scripts/goal2_sex_stratified/16d_nested_cv_reconciliation.R#L150));
   in cross-cohort transfer, `zrows` returns a vector of zeros for such a gene
   ([16b:122-124](scripts/goal2_sex_stratified/16b_model_training_final_panel_noMHC.R#L122-L124)).
   Both are now described, separately.
5. **Two stale code comments to fix, so that no reader quotes them.** The header
   of `14_model_training_nested_cv.R:18` still gives superseded candidate-set
   sizes ("14 female / 40 male") and the header of
   `18_testing_blood_pergene_roc.R:10-11` still names superseded panel genes
   (BNIP2, NMI; CLSTN1, GABBR1, HLA-DMA, SSRP1). The current sizes and panel
   membership are those written by the scripts into `mr_fs_summary.csv` and
   `mr_fs_summary_noMHC.csv`; quote the tables, never the comments. Candidate and
   panel sizes are results of §2.6 and §2.8 and are therefore not stated in this
   methods section.
6. **Per-gene orientation.** Fixed on training and carried unchanged across
   datasets in the cross-dataset table ([17:168-173, 223-235](scripts/goal2_sex_stratified/17_testing_blood_internal_external.R#L168-L235)).
   The train-only per-gene overlay figure instead orients each gene to its better
   direction within training ([18:34](scripts/goal2_sex_stratified/18_testing_blood_pergene_roc.R#L34)),
   which is legitimate for a within-training display but is *not* a transfer
   result; the prose now separates the two so the figure cannot be read as
   external validation.
7. **Frozen normalisation now carries its citations.** The holdout is projected
   using a frozen training quantile target and training-estimated ComBat
   parameters ([04_apply_holdout.R:59-69](scripts/00_shared/04_apply_holdout.R#L59-L69)),
   which requires Bolstad et al. (2003), Johnson, Li and Rabinovic (2007), Leek
   et al. (2012) and Nygaard, Rødland and Hovig (2016) rather than a bare
   assertion.
8. **Separation in the male stratum is now declared as a modelling condition**,
   together with the ridge stabilisation it forces in the utility analyses
   ([19:39-41](scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R#L39-L41))
   and the ≥ 0.999 flag it triggers in reporting.
9. **Seed policy made specific.** Global seed 1234; outer repeats seeded
   1000 + repeat index in the authoritative implementation. The superseded
   elastic-net loop in `15_model_training_elasticnet.R:105` uses 2000 + repeat
   index; no reported figure comes from it, and the text now says so.
10. **In-text citations added at every methodological claim**, and the reference
    list completed accordingly (39 entries, all cited in the text below).
11. **The claim that the fixed candidate universe introduces no leakage was
    wrong as written, and is the most exposed sentence the section contained.**
    The universe is an *externally filtered internal list*: the genes entering MR
    were `disease module ∩ sex-stratified DEG`, both computed on the whole
    training partition using the diagnosis labels, and only the *filter* applied
    to them was external. Because that universe is read once, outside the fold
    loop ([16d:76-80](scripts/goal2_sex_stratified/16d_nested_cv_reconciliation.R#L76-L80)),
    outer-fold test samples helped define the gene list their own predictions are
    built from. §2.9.2 now states this, bounds it, and scopes the nested
    estimate accordingly; §2.9.6 no longer says "the entire selection and fitting
    pipeline" was nested. The held-out and external estimates are untouched by
    the argument and are identified as the basis of the diagnostic claim.
12. **The per-gene orientation claim contradicted the code for the synovial set.**
    §2.9.7 asserted that training orientation was applied unchanged to the
    internal, external blood *and external synovial* sets, and that no gene was
    silently re-oriented. [20_testing_synovium_external.R:50-51](scripts/goal2_sex_stratified/20_testing_synovium_external.R#L50-L51)
    flips direction whenever AUC < 0.5, so the stored synovial per-gene values are
    best-direction. Both conventions are now defined, each is scoped to where it
    is genuinely used, and the cross-dataset re-expression performed in
    [24_crosstissue_pergene_auc.R:37-40](scripts/goal2_sex_stratified/24_crosstissue_pergene_auc.R#L37-L40)
    is named as the bridge between them.
13. **Three structurally different estimators shared one name.** "Panel AUC in
    dataset X" denoted a locked-transfer estimate in `PANEL_primary_vs_noMHC_performance.csv`,
    a within-dataset resampled estimate in `PANEL_auc_celladjusted.csv` and a
    nested estimate in `NESTED_CV_AUTHORITATIVE.csv`. §2.9.7 now opens by
    defining the three families, states that they answer different questions, and
    fixes the locked-transfer family as the basis of every transportability claim.

> **Duplication warning.** `METHODS_2.10_sexstratified_blood_biomarker.md`
> covers the same section. Sections 2.9.1–2.9.10 below supersede it; keep one
> file in the thesis or the two will drift apart again.

---

## 2.9.1 Overview and sex stratification

Sex-stratified consensus panels were used to construct and evaluate a diagnostic
classifier for rheumatoid arthritis (RA), with females and males modelled as
entirely separate analyses. No sample, cross-validation fold, tuning parameter or
fitted coefficient was exchanged between the two strata, so that the female panel
was evaluated only in females and the male panel only in males, and no estimate
reported in this chapter derives from a model that has seen data from both sexes.
Complete stratification, rather than the inclusion of sex as a covariate, was
adopted because the analytical question concerns whether the discriminating genes
themselves differ by sex, which a common-slope adjustment cannot address (Klein
and Flanagan, 2016; Oliva et al., 2020).

Two classifiers were developed for each sex. The first was a multivariable
logistic regression model fitted to the genes of the three-method consensus panel
of §2.8. The second was a penalised logistic regression fitted across the entire
within-sex set of genetically prioritised candidate genes, which permits the model
to retain correlated but informative genes that a strict three-way intersection
necessarily discards (Zou and Hastie, 2005). The two are distinct models, are
reported separately rather than pooled, and are never differenced.

## 2.9.2 Candidate predictors and data partitioning

The candidate gene universe was fixed before any expression-based modelling. It
comprised, for each sex, those genes whose *cis*-eQTL-instrumented expression was
associated with RA liability in the two-sample Mendelian randomisation analysis
of §2.6 at a Benjamini–Hochberg false discovery rate below 0.05 (Benjamini and
Hochberg, 1995). A second, restricted universe was defined in parallel by
excluding every candidate lying within the major histocompatibility complex, and
the whole of the development and evaluation procedure described below was run
independently on each universe as a pre-specified sensitivity analysis.

**The provenance of that universe must be stated precisely, because it determines
what the nested cross-validation of §2.9.6 does and does not correct for.** The
universe is the product of two stages of a different character. The genes
*entering* the Mendelian randomisation were the candidate set of §2.5, defined as
the intersection of the disease-associated co-expression modules of §2.4 with the
sex-stratified differentially expressed genes of §2.3; both of those analyses were
computed on the whole of the training partition and both used the diagnosis labels
of the training samples. The *filter applied to* that set was Mendelian
randomisation against external genome-wide association and eQTL summary statistics
(Okada et al., 2014; Võsa et al., 2021), which used no expression data from this
study whatever.

It follows that the universe is an externally filtered internal list, not an
external list, and the frequently made claim that a genetically defined candidate
set is immune to selection bias does not hold here without qualification. The
universe was held fixed across resampling folds — re-deriving the network,
differential expression and Mendelian randomisation inside every one of the fifty
outer folds was computationally out of reach — with the consequence that a sample
in an outer test fold contributed, through the training-wide differential
expression and module–trait analyses, to defining the gene list from which its own
prediction is subsequently built. The nested estimate of §2.9.6 therefore
quantifies selection bias arising in the three-selector stage and in the
classifier, and does not quantify selection bias arising in the upstream
differential-expression, network and prioritisation stages. It is accordingly an
upper bound on out-of-sample performance in that respect, and is reported as such
rather than as a fully leakage-free estimate (Ambroise and McLachlan, 2002; Simon
et al., 2003; Kaufman, Rosset and Perlich, 2012).

Three considerations bound the magnitude of the residual optimism without
eliminating it, and they are given rather than the assertion they replace. The
upstream stages are severe filters that are nonetheless indifferent to the
classifier: the module–trait rule of §2.4 selects on the correlation between a
module eigengene and diagnosis rather than on any gene's individual discriminative
power, and the Mendelian randomisation filter of §2.6 is estimated entirely
outside this cohort, so neither stage can select a gene *because* it separates
these particular samples. The candidate universe is small relative to the
transcriptome, comprising a few tens of genes drawn from 15,763, so the
combinatorial freedom available to the fold-level selectors — which is the
quantity that governs selection bias (Ambroise and McLachlan, 2002) — is
correspondingly small. And the two evaluations that carry the diagnostic claim,
namely the sealed internal holdout and the independent external cohorts of
§2.9.7, are unaffected by this argument entirely, because those samples took no
part in the differential expression, network construction, prioritisation or
selection at any stage and were projected onto the training scale using frozen
parameters. The residual optimism is confined to the training-cohort nested
estimate, and the held-out estimates are the ones on which the diagnostic claim
rests.

The magnitude attainable when selection is genuinely unconstrained is bounded
empirically by the supplementary transcriptome-wide benchmark described in
§2.9.9, in which differential expression and penalised selection are both
re-derived inside every outer training fold beginning from all 15,763 genes.

All model development was confined to the 70% training partition of §2.2. The
sealed 30% internal-validation holdout and the two independent external cohorts
were reserved for the evaluation of §2.9.7 and contributed to no selection,
tuning or fitting step. Quantile normalisation and empirical Bayes batch
correction had likewise been estimated on the training partition alone (Bolstad
et al., 2003; Johnson, Li and Rabinovic, 2007; Leek et al., 2012), the resulting
quantile target and batch parameters being frozen and applied to the holdout only
at the point of evaluation. A fresh joint normalisation would have allowed the
holdout to influence the input representation of the model, and ComBat in
particular uses the outcome design matrix when estimating its adjustment, so
re-estimating it on pooled data would import outcome information from the sealed
samples (Ambroise and McLachlan, 2002; Simon et al., 2003; Nygaard, Rødland and
Hovig, 2016).

## 2.9.3 Standardisation of expression values

To ensure applicability across microarray platforms, each gene was standardised
to a z-score before entering a model, under two regimes that are not
interchangeable and are distinguished here for that reason.

Within any resampling procedure conducted inside a single cohort, the mean and
standard deviation of each gene were computed from the training fold alone and
were then applied, unchanged, to the held-out fold, so that no test sample
influenced the standardisation; where the standard deviation of a gene was zero
or undefined it was replaced by unity, leaving that gene centred but unscaled
rather than yielding an undefined value. In the elastic-net path the raw
fold-training matrix was instead passed to `cv.glmnet` with internal
standardisation enabled, so that the centring and scaling constants were derived
from the fold-training rows by the fitting routine itself and coefficients were
returned on the original scale (Friedman, Hastie and Tibshirani, 2010); the
guarantee is the same, in that no held-out sample contributes to the
transformation.

When a locked model was transferred to a different cohort, each gene was
standardised within that target dataset, and within the relevant sex stratum of
it, independently. This was necessary because the internal holdout and the
external cohorts were generated on different array platforms, so that their
absolute intensities are not on a common scale, and it is legitimate because the
transformation is entirely unsupervised, using only the expression values of the
target dataset and never its diagnostic labels (Kaufman, Rosset and Perlich,
2012). A gene with zero or undefined variance in a target dataset was set to zero
across that dataset, and a panel gene absent from a target platform's annotation
was likewise assigned a z-score of zero, that is, the dataset mean, so that its
contribution to the linear predictor vanished rather than the affected samples
being discarded. The genes missing from each dataset were recorded alongside the
corresponding performance estimate rather than being suppressed.

## 2.9.4 Model specification

The consensus-panel classifier was a multivariable logistic regression model
(`glm`, family = binomial) in which the standardised consensus panel genes were
entered as predictors and diagnosis was the response, coded as RA versus healthy
control with healthy control as the reference level.

The penalised classifier was an elastic-net logistic regression (Zou and Hastie,
2005), fitted by cyclical coordinate descent in glmnet (Friedman, Hastie and
Tibshirani, 2010), a formulation that spans the ridge and least absolute
shrinkage and selection operator penalties through the mixing parameter α
(Hoerl and Kennard, 1970; Tibshirani, 1996). The mixing parameter was tuned over
the grid α ∈ {0.1, 0.3, 0.5, 0.7, 0.9, 1.0} and, for each value of α, the penalty
strength λ was tuned by inner five-fold cross-validation on binomial deviance;
the pair (α, λ) minimising the inner cross-validated deviance was selected
jointly and predictions were taken at λ_min. Because sparsity is an outcome of
tuning under this specification rather than a fixed design choice, the median
number of non-zero coefficients across folds is reported alongside every
elastic-net estimate.

## 2.9.5 Model training

For each sex, one model of each type was trained on the complete training
partition and then locked, such that no coefficient, tuning parameter or
standardisation constant was subsequently altered. The design matrix was
assembled with samples as rows and the relevant predictor genes as columns and
standardised as described in §2.9.3; all hyperparameters were tuned exclusively
by cross-validation internal to the training partition; and the model was
refitted at the selected values on the whole training partition, its coefficients
being recorded.

On the full training set the LASSO penalty was tuned by ten-fold
cross-validation on binomial deviance, yielding λ_min = 0.0434 in females and
0.0204 in males; the random forest was grown with 1,000 trees, raised from the
package default of 500 for stability of the importance ranking (Breiman, 2001;
Liaw and Wiener, 2002), with the number of predictors sampled at each split tuned
by ten-fold cross-validation over the grid {1, 2, ⌊√p⌋, ⌊p/3⌋, ⌊p/2⌋, p} using
cross-validated AUC as the tuning metric (Kuhn, 2008), yielding mtry = 2 in
females and mtry = 1 in males; and the linear support-vector cost was tuned by
ten-fold cross-validation over the grid {0.01, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16}
(Cortes and Vapnik, 1995; Chang and Lin, 2011), yielding C = 0.01 in both sexes,
with the panel size taken from a ten-fold cross-validated error curve (Guyon et
al., 2002).

Convergence of the unpenalised logistic models was monitored, and complete or
quasi-complete separation was treated as a diagnostic signal rather than as a
successful fit, since separation at small sample size produces unbounded
coefficients and degenerate standard errors (Albert and Anderson, 1984). Where
the panel separated the training data completely, as occurred in the male
stratum, a mild ridge penalty was applied for the clinical utility analyses of
§2.9.8 in order to stabilise the fit (Hoerl and Kennard, 1970; Harrell, 2015),
and the separation itself was reported as a condition of the fit rather than
absorbed into the estimate.

A single global seed of 1234 governed every stochastic step of training,
including cross-validation fold assignment, forest growth and bootstrap
resampling, and outer resampling repeats were seeded deterministically as
1000 + repeat index, so that all reported models and estimates are exactly
reproducible. Training used only the discovery cohorts, and at no point was a
model refitted after inspecting a validation result.

## 2.9.6 Nested cross-validation and the quantification of optimism

The consensus panel was selected using the entire training set. A conventional
cross-validated AUC computed on that same set would therefore be optimistically
biased, since each fold would be scoring genes that had already been chosen with
the assistance of the samples being held out. This phenomenon, termed
feature-selection bias, is sufficiently large in gene-expression data of this
dimensionality to generate apparently excellent classifiers from noise alone
(Ambroise and McLachlan, 2002; Simon et al., 2003). The three-selector procedure
and the classifier fitted on its output were therefore embedded within a nested
cross-validation framework, stratified by outcome within each sex (Varma and
Simon, 2006; Cawley and Talbot, 2010; Krstajic et al., 2014). What the fold loop
does and does not re-derive is stated exactly in §2.9.2: the candidate universe
is held fixed across folds, so the estimate below corrects the selection bias of
the panel-building stage and not that of the upstream differential-expression,
network and prioritisation stages.

The outer loop was a repeated stratified k-fold design: ten folds repeated five
times in females (n = 145) and five folds repeated ten times in males (n = 38).
Fewer folds were used in the male stratum because a ten-fold split at that sample
size leaves approximately four samples in each test fold, of which fewer still
are cases, so that a single sample changing side moves the fold-level AUC
substantially; the resulting instability was offset by increasing the number of
repeats, thereby averaging over fold assignment rather than conditioning on a
single arbitrary partition (Bengio and Grandvalet, 2004; Krstajic et al., 2014).
Folds were constructed with `caret` (Kuhn, 2008).

Within each outer training fold, and without any access to the outer test fold,
the complete three-selector procedure was re-derived. LASSO logistic regression
was refitted by coordinate descent (α = 1, family = binomial) with λ tuned by
inner five-fold cross-validation, and the genes retaining non-zero coefficients
at λ_min constituted its selection (Tibshirani, 1996; Friedman, Hastie and
Tibshirani, 2010). A random forest of 500 trees was grown, with the number of
predictors sampled at each split left at the package default of ⌊√p⌋, and those
genes whose mean decrease in Gini impurity exceeded the arithmetic mean of that
statistic across candidates constituted its selection (Breiman, 2001; Liaw and
Wiener, 2002). Recursive feature elimination was performed with a linear-kernel
support vector machine at fixed cost C = 1, eliminating at each iteration the
feature with the smallest squared weight and retaining the panel size that
minimised the inner five-fold cross-validated classification error (Cortes and
Vapnik, 1995; Guyon et al., 2002; Chang and Lin, 2011). The panel for that fold
was the intersection of the three selections wherever this contained at least two
genes; failing that, the procedure reverted, in a fixed and pre-specified order,
to the union of the three selections, then to the LASSO selection alone, and
finally to the complete candidate set, so that every fold yielded a usable panel
and no fold was discarded for want of one.

These in-fold selectors reproduce the structure of the selection procedure
applied to the full training set (§2.8) but deliberately not its hyperparameter
tuning: within the folds the forest was grown at 500 trees without tuning mtry
and the support vector machine used a fixed unit cost, only λ remaining
self-tuned, because re-tuning every hyperparameter within every fold of every
repeat was computationally prohibitive at this design size. The simplification is
conservative in the relevant direction, since an untuned selector is not more
likely than a tuned one to recover the panel obtained on the full training set,
and the resulting quantity is accordingly described as an estimate of the
pipeline under fixed selector hyperparameters rather than of the fully tuned
procedure (Cawley and Talbot, 2010).

A logistic regression classifier was subsequently fitted on each fold's panel,
standardised using the fold-training mean and standard deviation as described in
§2.9.3, and applied to predict the untouched outer test fold. In the
elastic-net variant the penalised model itself was refitted within each outer
training fold, with α and λ re-tuned by inner cross-validation on that fold's
data alone. The out-of-fold predicted probabilities were retained for every
sample. Since each sample is predicted exactly once within a repeat, an AUC was
computed for each repeat and summarised as a mean and standard deviation across
repeats; for the headline estimate, the out-of-fold probabilities of each sample
were averaged across repeats and a single receiver operating characteristic curve
was constructed from the averaged probabilities, with the case direction fixed a
priori rather than selected to maximise the AUC.

Three quantities were computed on the training data so that the influence of the
resampling design would be visible rather than assumed, and they are labelled
distinctly throughout because they are not comparable. The **apparent** AUC fits
the locked panel to all training samples and scores those same samples without
any resampling whatever; since the panel had been selected using those samples it
is an optimistic upper anchor by construction and is never presented as a
validation estimate. The **flat** cross-validated AUC selects the panel once on
the whole training set and cross-validates only the logistic model, with
z-scoring re-estimated inside each training fold; it is an explicit
representation of the biased procedure that the nested design exists to displace,
and it is reported for that purpose alone rather than as an honest estimate of
performance. The **nested** AUC re-derives selection within every fold and is the
figure quoted for the training cohort. Optimism was reported as the apparent
estimate minus the nested estimate, following established practice in the
internal validation of prediction models (Harrell, Lee and Mark, 1996;
Steyerberg et al., 2001).

The difference between the flat and nested estimates was deliberately not
reported as the magnitude of selection bias, because two effects act upon it in
opposing directions: the flat estimate is inflated by the panel having been
chosen using all of the data, whereas the nested estimate pools out-of-fold
probabilities across repeats in which each fold selected its own panel, so that
averaging those probabilities constitutes a form of ensembling and tends to raise
the AUC (Dietterich, 2000). The two estimates therefore do not form a clean
decomposition of bias, and both are tabulated with this caveat rather than
differenced.

Because the panel is re-derived within every fold, the frequency with which each
gene of the reported panel re-entered the in-fold three-way consensus provides a
direct measure of the reproducibility of the selection under resampling, and was
reported as a percentage of fitted folds for every panel gene (Meinshausen and
Bühlmann, 2010). A gene selected on the full training set but seldom recovered
within folds was treated as an unstable selection rather than as a validated
biomarker.

Finally, several analyses in the pipeline reported a quantity that could
reasonably be taken to be a single nested cross-validated AUC, and these did not
agree, because they in fact measured different procedures. Three axes of
variation were involved: the candidate set, comprising either the primary
genetically prioritised set or the set remaining after exclusion of the major
histocompatibility complex; the classifier, comprising either the three-selector
consensus with a logistic model or the tuned elastic net; and the resampling
scheme, which differs between the sexes as described above. All eight
combinations were therefore recomputed within a single process under a single
seed policy, with the three axes recorded as explicit columns of one
authoritative table, and the earlier duplicate implementations were retired
rather than reconciled line by line. The variants were not constrained to agree,
since they estimate the performance of genuinely different models; one row was
designated as the headline estimate and the remainder are reported as a
sensitivity grid. The dispersion observed across procedures in the male stratum
is interpreted as a consequence of that stratum's sample size rather than as
evidence favouring one procedure over another, and the male arm is reported
throughout as a power-limited exploratory analysis.

## 2.9.7 Model evaluation

### The three estimator families, and why they are named separately

Several analyses in this chapter report a quantity that a reader would reasonably
take to be "the area under the curve of the panel in dataset X", and they do not
agree, because they estimate three structurally different things. They are named
distinctly here, every table states which family it belongs to, and values from
different families are never compared with one another or averaged.

A **locked-transfer** estimate applies the model fitted and frozen on the training
partition — its panel, its coefficients and its standardisation rule — to a dataset
that contributed to none of them, and scores it once. This is the estimator that
answers whether the diagnostic model works on new patients, and it is the estimator
on which every claim of transportability in this chapter rests.

A **within-dataset resampled** estimate takes only the *identity* of the panel
genes from training and re-estimates the classifier's coefficients inside the
evaluation dataset under cross-validation. It answers whether the gene set carries
discriminative information in that dataset, which is a weaker and different
question, and it is the appropriate estimator when the evaluation data are on a
scale or platform that makes transported coefficients meaningless, as in the
synovial analysis of §2.11, or when the panel must be compared like-for-like
against a competing model fitted in the same dataset, as in the
cell-composition benchmark. A within-dataset resampled estimate is never a test of
transfer, because the coefficients were fitted on the samples being scored, and it
is never quoted as out-of-sample performance of the locked model.

A **nested** estimate re-derives the panel-building procedure inside every outer
training fold of the training partition, as described in §2.9.6, and estimates the
performance of the *procedure* rather than of any one fitted model.

The consequence is stated so that it cannot be lost: a within-dataset resampled
estimate and a locked-transfer estimate computed on the same samples answer
different questions and will differ, sometimes substantially, and neither is a
correction of the other. Where both are available for a dataset, both are reported,
labelled by family, and the locked-transfer value is the one carried into any
statement about diagnostic performance.

### Settings of evaluation

Each locked model was evaluated in three settings of increasing independence from
the data on which it was trained (Justice, Covinsky and Berlin, 1999; Steyerberg
and Harrell, 2016).

The first was the training cohort itself, where the apparent and nested estimates
defined in §2.9.6 are reported, the former explicitly as an optimistic upper
bound and the latter as the figure quoted for the training cohort.

The second was the sealed internal-validation holdout, comprising the 30% of the
discovery cohorts withheld before any cross-sample operation and projected onto
the training scale using the frozen quantile target and batch-correction
parameters of §2.9.2, so that no holdout sample influenced any training
parameter. The locked model was applied to these samples once.

The third was an independent external blood cohort profiled on a different
platform, which constitutes a genuine external test of transportability rather
than a further internal split (Justice, Covinsky and Berlin, 1999; Steyerberg and
Harrell, 2016). Intensities were placed on the log₂ scale where they were not
already, the test being whether the maximum observed value exceeded 50, and probe
sets were collapsed to gene symbols by retaining, for each symbol, the probe with
the highest mean expression across samples, a rule that outperforms averaging
across probes when probes for one symbol differ in specificity (Miller et al.,
2011). The locked model was applied to these samples once.

No element of the evaluation was tuned to its result. The direction of the
receiver operating characteristic comparison was fixed a priori for every
dataset, models were never scored on the samples used to fit them except where
explicitly labelled as apparent performance, and no random seed was chosen on the
basis of the estimate it produced.

In a complementary univariate analysis, each panel gene was evaluated alone as a
classifier in every dataset. Because the area under the curve of a single gene
depends on which direction of expression is taken to mark disease, and because
that choice can be made either before or after seeing the evaluation data, two
orientation conventions arise and are distinguished throughout; they yield
different values for the same gene and must never be read interchangeably.

Under the **training-fixed** convention the direction of expression taken to mark
RA is determined on the training data and applied unchanged to the evaluation
dataset, so that a gene whose association reverses out of sample yields an area
under the curve below 0.5. This is the convention used for the internal holdout
and the external blood cohort, and it is the convention used wherever per-gene
estimates from several datasets are tabulated or plotted together, including the
external synovial cohort. It is the only convention under which a failure of
transfer remains visible, and it is therefore the one on which any statement about
transportability is based.

Under the **best-direction** convention the orientation is chosen within the
evaluation dataset itself, so that every reported value is at least 0.5 by
construction. This convention answers a different and narrower question, namely
how much information a gene carries within that dataset irrespective of the
direction of its association, and it cannot support a claim of transfer. It is
used in two places, both confined to a single dataset: the within-training overlay
figure, and the synovial discovery table and its accompanying synovium-only
overlay, in which the direction of effect between blood and synovium is carried
separately by an explicit concordance flag rather than by the sign of the
statistic. Every table and figure states which convention it uses, and the
cross-dataset comparisons of §2.11 re-express the synovial values on the
training-fixed convention before placing them alongside the blood estimates.

## 2.9.8 Evaluation of clinical utility

Discrimination alone is an insufficient basis on which to recommend a diagnostic
model, since a model may discriminate well yet be poorly calibrated or offer no
net benefit at clinically relevant decision thresholds (Vickers and Elkin, 2006;
Steyerberg et al., 2010). Three further analyses were therefore conducted on the
training cohort for each sex, using logistic models fitted with the rms package
(Harrell, 2015).

Calibration, the agreement between predicted and observed risk, was assessed by
plotting observed against predicted probability, both apparently and after bias
correction by bootstrap resampling with 200 repetitions, which corrects the
calibration curve for the optimism induced by fitting and evaluating on the same
data (Harrell, Lee and Mark, 1996; Steyerberg et al., 2001).

Net benefit was assessed by decision curve analysis across threshold
probabilities from 0.01 to 0.99, comparing the panel model against the two
default strategies of treating all patients and treating none (Vickers and Elkin,
2006). Net benefit was computed as the proportion of true positives minus the
proportion of false positives weighted by the odds of the threshold probability,
and the threshold axis was additionally annotated with the corresponding
cost-to-benefit ratio.

Clinical impact was expressed, per 1,000 hypothetical individuals, as the number
classified as high risk and the number classified as high risk who genuinely have
RA, across the same range of thresholds, with 95% confidence bands obtained from
500 bootstrap resamples of the training data (Kerr et al., 2016). A nomogram was
constructed from the fitted logistic model to render the score transparent and
hand-computable at the point of use (Harrell, 2015).

These analyses were conducted on the training cohort and are therefore internal;
they describe the potential utility of the panel under the fitted model rather
than demonstrating utility in an independent population, and are reported as
such.

## 2.9.9 Statistical inference, small-sample safeguards and evidence tiering

Areas under the receiver operating characteristic curve and their 95% confidence
intervals were computed with pROC (Robin et al., 2011), intervals being obtained
by the method of DeLong, DeLong and Clarke-Pearson (1988). Where a stratum
comprised fewer than twenty samples, for which the DeLong interval is unreliable
and may report spurious precision, a stratified bootstrap interval with 2,000
resamples was substituted (Carpenter and Bithell, 2000). Paired comparisons
between panels were made on the same samples using the DeLong test for
correlated curves (DeLong, DeLong and Clarke-Pearson, 1988).

Two further safeguards were applied because of the limited size of the male
stratum. First, any AUC at or above 0.999 was explicitly flagged as indicating
separation. At small sample sizes the number of case–control pairs available to be
discordant is itself small, and the AUC is by construction the probability that a
randomly selected case is ranked above a randomly selected control (Hanley and
McNeil, 1982); a perfect value under such conditions is an unremarkable
consequence of sample size rather than evidence of a perfect classifier, and the
degenerate interval accompanying it must not be interpreted as precision.
Second, each analysis was assigned an evidence tier, either primary or
exploratory, on the basis of sample size alone and never on the magnitude of the
estimate, using pre-specified thresholds applied to the training, internal and
external sample counts. The male stratum falls below these thresholds, its case
and control counts yielding a ratio of events to candidate predictors far below
conventional minima for a stable multivariable model (Peduzzi et al., 1996;
Riley et al., 2019). It is consequently reported throughout as a power-limited
exploratory analysis, is not presented alongside the female estimates as though
the two were comparable, and supports no diagnostic claim. Reporting follows the
recommendations of the Transparent Reporting of a multivariable prediction model
for Individual Prognosis Or Diagnosis statement (Collins et al., 2015).

A supplementary benchmark analysis was performed in which differential
expression and LASSO selection were both conducted within each outer training
fold, beginning from the transcriptome rather than from the genetically
prioritised candidate set. This analysis omits the co-expression, Mendelian
randomisation and three-selector consensus stages from the fold loop and is
therefore not a leakage-free validation of the reported panels. It is presented
for the single purpose of bounding the magnitude of feature-selection bias
attainable on these data when selection is unconstrained (Ambroise and McLachlan,
2002), and is identified as a benchmark wherever it appears.

## 2.9.10 Software and reproducibility

Analyses were conducted in R v4.4.2 using glmnet v5.0 (Friedman, Hastie and
Tibshirani, 2010), randomForest v4.7.1.2 (Liaw and Wiener, 2002), e1071 v1.7.17
interfacing LIBSVM (Chang and Lin, 2011), caret v7.0.1 for stratified fold
construction and hyperparameter grid search (Kuhn, 2008), pROC v1.19.0.1 (Robin
et al., 2011), rms v8.1.0 (Harrell, 2015), limma v3.62.2 (Ritchie et al., 2015),
sva v3.54.0 (Leek et al., 2012), preprocessCore v1.68.0 for frozen quantile
normalisation (Bolstad et al., 2003) and data.table v1.18.2.1. Package versions
were captured at run time from the session information and recorded in a
lockfile. Outer fold assignment was seeded as a fixed base value of 1000
incremented by the repeat index and a single global seed of 1234 governed every
other stochastic step, so that every repeat is exactly reproducible and every
variant contributing to the reconciled table was generated under an identical
seed policy.

---

## References cited in this section

Albert, A. and Anderson, J.A. (1984) 'On the existence of maximum likelihood estimates in logistic regression models', *Biometrika*, 71(1), pp. 1–10.

Ambroise, C. and McLachlan, G.J. (2002) 'Selection bias in gene extraction on the basis of microarray gene-expression data', *Proceedings of the National Academy of Sciences*, 99(10), pp. 6562–6566.

Bengio, Y. and Grandvalet, Y. (2004) 'No unbiased estimator of the variance of K-fold cross-validation', *Journal of Machine Learning Research*, 5, pp. 1089–1105.

Benjamini, Y. and Hochberg, Y. (1995) 'Controlling the false discovery rate: a practical and powerful approach to multiple testing', *Journal of the Royal Statistical Society: Series B*, 57(1), pp. 289–300.

Bolstad, B.M., Irizarry, R.A., Åstrand, M. and Speed, T.P. (2003) 'A comparison of normalization methods for high density oligonucleotide array data based on variance and bias', *Bioinformatics*, 19(2), pp. 185–193.

Breiman, L. (2001) 'Random forests', *Machine Learning*, 45(1), pp. 5–32.

Carpenter, J. and Bithell, J. (2000) 'Bootstrap confidence intervals: when, which, what? A practical guide for medical statisticians', *Statistics in Medicine*, 19(9), pp. 1141–1164.

Cawley, G.C. and Talbot, N.L.C. (2010) 'On over-fitting in model selection and subsequent selection bias in performance evaluation', *Journal of Machine Learning Research*, 11, pp. 2079–2107.

Chang, C.-C. and Lin, C.-J. (2011) 'LIBSVM: a library for support vector machines', *ACM Transactions on Intelligent Systems and Technology*, 2(3), 27.

Collins, G.S., Reitsma, J.B., Altman, D.G. and Moons, K.G.M. (2015) 'Transparent reporting of a multivariable prediction model for individual prognosis or diagnosis (TRIPOD): the TRIPOD statement', *Annals of Internal Medicine*, 162(1), pp. 55–63.

Cortes, C. and Vapnik, V. (1995) 'Support-vector networks', *Machine Learning*, 20(3), pp. 273–297.

DeLong, E.R., DeLong, D.M. and Clarke-Pearson, D.L. (1988) 'Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach', *Biometrics*, 44(3), pp. 837–845.

Dietterich, T.G. (2000) 'Ensemble methods in machine learning', in *Multiple Classifier Systems*. Lecture Notes in Computer Science 1857. Berlin: Springer, pp. 1–15.

Friedman, J., Hastie, T. and Tibshirani, R. (2010) 'Regularization paths for generalized linear models via coordinate descent', *Journal of Statistical Software*, 33(1), pp. 1–22.

Guyon, I., Weston, J., Barnhill, S. and Vapnik, V. (2002) 'Gene selection for cancer classification using support vector machines', *Machine Learning*, 46(1–3), pp. 389–422.

Hanley, J.A. and McNeil, B.J. (1982) 'The meaning and use of the area under a receiver operating characteristic (ROC) curve', *Radiology*, 143(1), pp. 29–36.

Harrell, F.E. (2015) *Regression Modeling Strategies: With Applications to Linear Models, Logistic and Ordinal Regression, and Survival Analysis*. 2nd edn. New York: Springer.

Harrell, F.E., Lee, K.L. and Mark, D.B. (1996) 'Multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors', *Statistics in Medicine*, 15(4), pp. 361–387.

Hoerl, A.E. and Kennard, R.W. (1970) 'Ridge regression: biased estimation for nonorthogonal problems', *Technometrics*, 12(1), pp. 55–67.

Johnson, W.E., Li, C. and Rabinovic, A. (2007) 'Adjusting batch effects in microarray expression data using empirical Bayes methods', *Biostatistics*, 8(1), pp. 118–127.

Justice, A.C., Covinsky, K.E. and Berlin, J.A. (1999) 'Assessing the generalizability of prognostic information', *Annals of Internal Medicine*, 130(6), pp. 515–524.

Kaufman, S., Rosset, S. and Perlich, C. (2012) 'Leakage in data mining: formulation, detection, and avoidance', *ACM Transactions on Knowledge Discovery from Data*, 6(4), 15.

Kerr, K.F., Brown, M.D., Zhu, K. and Janes, H. (2016) 'Assessing the clinical impact of risk prediction models with decision curves: guidance for correct interpretation and appropriate use', *Journal of Clinical Oncology*, 34(21), pp. 2534–2540.

Klein, S.L. and Flanagan, K.L. (2016) 'Sex differences in immune responses', *Nature Reviews Immunology*, 16(10), pp. 626–638.

Krstajic, D., Buturovic, L.J., Leahy, D.E. and Thomas, S. (2014) 'Cross-validation pitfalls when selecting and assessing regression and classification models', *Journal of Cheminformatics*, 6, 10.

Kuhn, M. (2008) 'Building predictive models in R using the caret package', *Journal of Statistical Software*, 28(5), pp. 1–26.

Leek, J.T., Johnson, W.E., Parker, H.S., Jaffe, A.E. and Storey, J.D. (2012) 'The sva package for removing batch effects and other unwanted variation in high-throughput experiments', *Bioinformatics*, 28(6), pp. 882–883.

Liaw, A. and Wiener, M. (2002) 'Classification and regression by randomForest', *R News*, 2(3), pp. 18–22.

Meinshausen, N. and Bühlmann, P. (2010) 'Stability selection', *Journal of the Royal Statistical Society: Series B*, 72(4), pp. 417–473.

Miller, J.A., Cai, C., Langfelder, P., Geschwind, D.H., Kurian, S.M., Salomon, D.R. and Horvath, S. (2011) 'Strategies for aggregating gene expression data: the collapseRows R function', *BMC Bioinformatics*, 12, 322.

Nygaard, V., Rødland, E.A. and Hovig, E. (2016) 'Methodological variations in ComBat and their impact on the analysis of gene expression data', *Briefings in Bioinformatics*, 17(1), pp. 29–39.

Okada, Y., Wu, D., Trynka, G., Raj, T., Terao, C., Ikari, K., et al. (2014) 'Genetics of rheumatoid arthritis contributes to biology and drug discovery', *Nature*, 506(7488), pp. 376–381.

Oliva, M., Muñoz-Aguirre, M., Kim-Hellmuth, S., Wucher, V., Gewirtz, A.D.H., Cotter, D.J., et al. (2020) 'The impact of sex on gene expression across human tissues', *Science*, 369(6509), eaba3066.

Peduzzi, P., Concato, J., Kemper, E., Holford, T.R. and Feinstein, A.R. (1996) 'A simulation study of the number of events per variable in logistic regression analysis', *Journal of Clinical Epidemiology*, 49(12), pp. 1373–1379.

Riley, R.D., Snell, K.I.E., Ensor, J., Burke, D.L., Harrell, F.E., Moons, K.G.M. and Collins, G.S. (2019) 'Minimum sample size for developing a multivariable prediction model: Part II, binary and time-to-event outcomes', *Statistics in Medicine*, 38(7), pp. 1276–1296.

Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W. and Smyth, G.K. (2015) 'limma powers differential expression analyses for RNA-sequencing and microarray studies', *Nucleic Acids Research*, 43(7), e47.

Robin, X., Turck, N., Hainard, A., Tiberti, N., Lisacek, F., Sanchez, J.-C. and Müller, M. (2011) 'pROC: an open-source package for R and S+ to analyze and compare ROC curves', *BMC Bioinformatics*, 12, 77.

Simon, R., Radmacher, M.D., Dobbin, K. and McShane, L.M. (2003) 'Pitfalls in the use of DNA microarray data for diagnostic and prognostic classification', *Journal of the National Cancer Institute*, 95(1), pp. 14–18.

Steyerberg, E.W. and Harrell, F.E. (2016) 'Prediction models need appropriate internal, internal-external, and external validation', *Journal of Clinical Epidemiology*, 69, pp. 245–247.

Steyerberg, E.W., Harrell, F.E., Borsboom, G.J.J.M., Eijkemans, M.J.C., Vergouwe, Y. and Habbema, J.D.F. (2001) 'Internal validation of predictive models: efficiency of some procedures for logistic regression analysis', *Journal of Clinical Epidemiology*, 54(8), pp. 774–781.

Steyerberg, E.W., Vickers, A.J., Cook, N.R., Gerds, T., Gonen, M., Obuchowski, N., Pencina, M.J. and Kattan, M.W. (2010) 'Assessing the performance of prediction models: a framework for traditional and novel measures', *Epidemiology*, 21(1), pp. 128–138.

Tibshirani, R. (1996) 'Regression shrinkage and selection via the lasso', *Journal of the Royal Statistical Society: Series B*, 58(1), pp. 267–288.

Varma, S. and Simon, R. (2006) 'Bias in error estimation when using cross-validation for model selection', *BMC Bioinformatics*, 7, 91.

Vickers, A.J. and Elkin, E.B. (2006) 'Decision curve analysis: a novel method for evaluating prediction models', *Medical Decision Making*, 26(6), pp. 565–574.

Võsa, U., Claringbould, A., Westra, H.-J., Bonder, M.J., Deelen, P., Zeng, B., et al. (2021) 'Large-scale cis- and trans-eQTL analyses identify thousands of genetic loci and polygenic scores that regulate blood gene expression', *Nature Genetics*, 53(9), pp. 1300–1310.

Zou, H. and Hastie, T. (2005) 'Regularization and variable selection via the elastic net', *Journal of the Royal Statistical Society: Series B*, 67(2), pp. 301–320.

---

## Parameter appendix (verified against code and output tables; not for the thesis body)

| Step | Parameter | Value | Source |
|---|---|---|---|
| Global seed | `GLOBAL_SEED` | 1234 | 12, 12b, 16, 16b, 16d |
| Outer fold seed | per repeat | 1000 + repeat index | 14:97, 16d:139 |
| Superseded seed offset | elastic-net loop in 15 | 2000 + repeat index (no reported figure) | 15:105 |
| Outer resampling (F) | k × repeats | 10 × 5, n = 145 | 16d:82 |
| Outer resampling (M) | k × repeats | 5 × 10, n = 38 | 16d:82 |
| Candidate universes | primary / MHC-free | run independently through the whole procedure | 16d:76-80 |
| In-fold LASSO | α, nfolds, s | 1, 5, `lambda.min` | 16d:103-107 |
| In-fold RF | ntree, mtry, cutoff | 500, default ⌊√p⌋, Gini > mean | 16d:108-109 |
| In-fold SVM-RFE | kernel, cost, size CV | linear, C = 1, 5-fold error | 16d:88-102 |
| Fold panel rule | intersection ≥ 2 genes | else union ≥ 1 → LASSO → all candidates | 16d:111-114 |
| Full-set LASSO | nfolds, loss, λ_min | 10, deviance; 0.0434 (F) / 0.0204 (M) | 12:228-232, `mr_fs_summary.csv` |
| Full-set RF | ntree, mtry grid, metric | 1000; {1,2,⌊√p⌋,⌊p/3⌋,⌊p/2⌋,p}; ROC, 10-fold | 12:194-204 |
| Full-set RF tuned | mtry | 2 (F), 1 (M) | `mr_fs_summary.csv` |
| Full-set SVM | cost grid, CV, tuned | {0.01…16}, 10-fold; C = 0.01 both sexes | 12:183-189, `mr_fs_summary.csv` |
| Elastic net | α grid, λ, standardise | {0.1,0.3,0.5,0.7,0.9,1.0}, inner 5-fold, λ.min, `standardize = TRUE` | 16:56,86-88; 16d:116-124 |
| Standardisation (in-cohort) | frozen fold-train μ, σ | σ = 0 or NA → 1 | 14:107-110, 16d:149-152, 17:139-142 |
| Standardisation (transfer) | per-dataset, per-sex z | σ = 0 or NA → all zeros; missing gene → z = 0 | 16b:122-124,256; 17:106-108,186 |
| Flat CV (script 17) | single stratified 10-fold | z re-estimated per fold; **not** an honest estimate | 17:131-149 |
| Flat CV (script 14) | same repeated scheme as nested | fixed full-training panel | 14:152-182 |
| Apparent AUC | resampling | none; locked panel fitted and scored on the same samples | 14:186-197, 16b:262-264 |
| Optimism | definition | apparent − nested | 14:207-208, 245 |
| Reconciliation grid | axes × combinations | candidate set × selector × sex, 8 rows, one process, one seed policy | 16d:183-195 |
| MHC sensitivity test | statistic | paired DeLong test on the nested out-of-fold predictions | 16b:206-227 |
| Separation flag | threshold | AUC ≥ 0.999 labelled SEPARATION | 16b:133-137 |
| Evidence tiering | basis | sample size only (train < 50 or internal/external n < 20 → exploratory) | 16:135,149 |
| External blood | log₂ condition, probe collapse | log₂(x+1) if max > 50; MaxMean per symbol | 17:84-88, 16b:106-111 |
| ROC | direction, levels | fixed `"<"`, `c("HC","RA")` | all |
| CI | DeLong; bootstrap if n < 20 | `boot.n = 2000` | 16b:126-132, 16d:171-173, 17:112-118 |
| Calibration | bootstrap reps | B = 200 (`rms::calibrate`) | 19:58 |
| Clinical impact | bootstrap reps | B = 500 | 19:85 |
| DCA | threshold range | 0.01 to 0.99 by 0.01 | 19:43 |
| Male utility fit | penalty | mild ridge (separation at n = 38) | 19:39-41 |
| Software | R, key packages | R 4.4.2; glmnet 5.0, randomForest 4.7.1.2, e1071 1.7.17, caret 7.0.1, pROC 1.19.0.1, rms 8.1.0, limma 3.62.2, sva 3.54.0, preprocessCore 1.68.0, data.table 1.18.2.1 | `renv.lock` |
