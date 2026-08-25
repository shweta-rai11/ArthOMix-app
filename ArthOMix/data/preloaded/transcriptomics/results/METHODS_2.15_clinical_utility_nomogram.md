# 2.15 Nomogram construction and clinical evaluation

## 2.15.1 Definitions of the four instruments

A **nomogram** is a graphical representation of a fitted regression model in which
each predictor is mapped onto a points scale, so that the linear predictor, and
hence the predicted probability, can be evaluated by hand for an individual
patient without recourse to software (Harrell, 2015; Iasonos et al., 2008). It
adds no information to the model from which it is drawn; its purpose is to make
that model transparent and computable at the point of use.

**Calibration** is the agreement between the probabilities a model predicts and
the frequencies with which the outcome actually occurs, and it is a property
distinct from discrimination: a model may rank cases above controls well, and so
achieve a high area under the receiver operating characteristic curve, while
systematically over- or under-stating absolute risk (Van Calster et al., 2019;
Alba et al., 2017). Because a curve obtained by fitting and evaluating on the same
samples lies optimistically close to the diagonal by construction, calibration is
reported both as the **apparent** curve and as a **bias-corrected** curve, the
latter obtained by resampling to estimate and remove that optimism (Harrell, Lee
and Mark, 1996; Steyerberg et al., 2001).

**Decision-curve analysis** evaluates whether acting on a model would do more good
than harm, by placing true and false positives on a common scale. Its central
quantity, the **net benefit**, is the proportion of true positives minus the
proportion of false positives weighted by the odds of the **threshold
probability** p_t, where p_t is the predicted risk at which a clinician would
choose to act and its odds, p_t/(1 − p_t), express the relative harm of an
unnecessary intervention against a missed case (Vickers and Elkin, 2006; Vickers,
van Calster and Steyerberg, 2016). A model is worth using at a given threshold
only if its net benefit exceeds that of both default strategies, treating every
patient and treating none.

A **clinical-impact curve** translates the same model into the units in which
consequences are counted, reporting, for a hypothetical population at each
threshold, how many individuals would be classified as high risk and how many of
those would genuinely have the disease (Kerr et al., 2016).

All four were computed separately within each sex, as described below.

## 2.15.2 Sex stratification

The procedure was carried out entirely within sex. Distinct instruments were
constructed for females and males, each based on a model fitted to the training
samples of that sex, using that sex's consensus panel genes, so that no sample,
coefficient, predictor axis or scaling constant was shared between them (Harrell,
2015). Sex was not included as a predictor; rather, it served as the partition
within which all subsequent computations were performed, resulting in instruments
devoid of a sex axis and inapplicable across sexes.

## 2.15.3 Model fitting

For each sex, a binary logistic regression model was fitted using `rms::lrm`, with
rheumatoid arthritis (RA) status regressed on the panel genes. The predictors were
incorporated as additive linear terms on the logit scale, without spline
expansions, interaction terms or categorisation of expression, thereby encoding an
assumption of linearity between each gene's log₂ expression and the log-odds of RA
across the observed range of that gene. The genes were entered as unstandardised
log₂ expression values, as they appear in the processed training matrix, rather
than as z-scores. A `datadist` object summarising the distribution of each
predictor within that sex was computed and registered prior to fitting, as rms
requires it to determine the plotted range of each predictor axis. For the male
panel a mild ridge penalty (penalty = 5) was applied to stabilise a fit that
almost completely separates the training data at that stratum's sample size,
whereas the female model was fitted without penalty (Hoerl and Kennard, 1970;
Albert and Anderson, 1984; Harrell, 2015). The design and response matrices were
retained (`x = TRUE`, `y = TRUE`) so that the resampling procedures below could be
applied to the fitted object.

## 2.15.4 Nomogram construction

The nomogram was generated from the fitted model using `rms::nomogram`. Each panel
gene was represented by its own axis, which spanned the range of that gene's
observed expression within the corresponding sex as recorded in the `datadist`
object, so that the axes were empirical and stratum-specific rather than nominal.
Points were allocated according to the standard nomogram scaling, in which the
predictor whose fitted contribution to the linear predictor spans the widest
interval was assigned a range of 0 to 100 points and every remaining predictor was
scaled in proportion to the width of its own contribution. The relative length of
the gene axes therefore reflects the product of each coefficient and the spread of
that gene's expression in that sex, not the coefficient alone, and is not a
ranking of effect size (Balachandran et al., 2015). The instrument was read by
identifying a patient's value on each gene axis, noting the points indicated above
it, summing those points across the panel genes, and projecting the total onto a
predicted-probability axis. That axis was derived by applying the logistic
transformation to the corresponding value of the linear predictor (`fun = plogis`),
was labelled "Risk of RA", and was annotated at predicted probabilities of 0.05,
0.1, 0.3, 0.5, 0.7, 0.9 and 0.99. Points were therefore expressed per unit of log₂
expression, and a reading is valid only for expression values placed on the same
processed scale as the training matrix from which the model was derived.

## 2.15.5 Calibration

Calibration was assessed within each sex by plotting observed against predicted
probability, reporting the apparent curve and a bootstrap bias-corrected curve
obtained from 200 resamples (`rms::calibrate`), so that the optimism induced by
fitting and evaluating on the same samples was estimated rather than assumed
absent (Harrell, Lee and Mark, 1996; Steyerberg et al., 2001; Van Calster et al.,
2019). Where the resampling procedure failed to return a calibration object, the
corresponding panel was left blank rather than replaced by an uncorrected curve.

## 2.15.6 Decision-curve analysis

Within each sex, net benefit was computed as

NB = TP/n − FP/n × p_t/(1 − p_t)

across threshold probabilities from 0.01 to 0.99 in steps of 0.01, and was compared
against the two default strategies of treating all patients, whose net benefit is
determined by the observed prevalence in that sex stratum alone, and treating none,
whose net benefit is zero by definition (Vickers and Elkin, 2006). Because the
treat-all reference is prevalence-dependent, the female and male decision curves
are each interpretable only against their own stratum's reference and cannot be
overlaid (Vickers, van Calster and Steyerberg, 2016). A secondary cost-to-benefit
ratio axis, p_t/(1 − p_t), was added for interpretation, annotated at ratios of
1:100, 1:4, 2:3, 3:2, 4:1 and 100:1. The probabilities entering this curve were the
apparent, resubstitution predictions of the fitted model, so the curve was not
corrected for optimism.

## 2.15.7 Clinical-impact analysis

Within each sex, the number of individuals per 1,000 classified as high risk, and
the number of those who genuinely have RA, were plotted across the same range of
thresholds (Kerr et al., 2016). Uncertainty was quantified by 500 bootstrap
resamples drawn with replacement from within that sex's training samples only, each
seeded by its own replicate index for reproducibility, with the model refitted in
each resample under the same penalty as the corresponding primary model, and with
the 2.5th and 97.5th percentiles across replicates providing 95% percentile bands
for both curves (Carpenter and Bithell, 2000). Within each replicate the refitted
model was evaluated on that replicate's own resampled observations, so the bands
describe the sampling variability of the apparent curves and do not correct them
for optimism.

## 2.15.8 Scope and limitations

Three limitations follow from this construction and are stated rather than left to
inference. First, each instrument is specific to the stratum in which it was built,
in its predictors, its axis ranges and its point scaling alike, so the female
nomogram is not applicable to a male patient or the reverse, and the point totals
of the two instruments are not comparable quantities. Second, the male instruments
derive from penalised coefficients while the female ones do not, so they rest on
different degrees of shrinkage, and the male arm is reported throughout as a
power-limited exploratory analysis. Third, all four analyses were conducted on the
training cohort, and only the calibration curve is bias-corrected; they therefore
describe the potential utility of the panel under the fitted model rather than
demonstrating utility in an independent population (Steyerberg and Harrell, 2016).

---

## Parameters used

| Item | Value |
|---|---|
| Stratification | females and males analysed independently; nothing shared between the two sets of instruments |
| Cohort | training cohort, that sex's samples only |
| Predictors | that sex's consensus panel genes, as unstandardised log₂ expression |
| Model | `rms::lrm`, RA (0/1) ~ panel genes, additive linear terms |
| Model form | no splines, no interactions, no categorisation, no sex term |
| Penalty | female 0; male 5 (ridge, for near-complete separation) |
| Stored components | `x = TRUE`, `y = TRUE` |
| Predictor metadata | `datadist` computed within sex, registered via `options(datadist = ...)` |
| Nomogram call | `nomogram(fit, fun = plogis, funlabel = "Risk of RA", fun.at = c(0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99))` |
| Nomogram axes | one per gene, spanning that gene's observed range within the sex |
| Points scaling | widest fitted contribution scaled 0 to 100; remaining predictors in proportion |
| Calibration | `calibrate(fit, B = 200)`; apparent plus bootstrap bias-corrected; blank panel on failure |
| Net benefit | NB = TP/n − FP/n × p_t/(1 − p_t) |
| Threshold grid | 0.01 to 0.99, step 0.01 |
| Reference strategies | treat all (prevalence-determined) and treat none (zero) |
| Cost:benefit axis | p_t/(1 − p_t), ticks at 1:100, 1:4, 2:3, 3:2, 4:1, 100:1 |
| Clinical impact | per 1,000: number high risk, and number high risk with RA |
| Bootstrap | B = 500, seeded per replicate, refitted with the same penalty, evaluated within the replicate |
| Interval type | 2.5th and 97.5th percentile bands across replicates |
| Optimism correction | calibration only; nomogram, decision and clinical-impact curves are apparent |
| Output | one composite four-panel figure per sex, plus one decision-curve table per sex |
| Software | R with rms, data.table and magick |

---

## References cited in this section

Alba, A.C., Agoritsas, T., Walsh, M., Hanna, S., Iorio, A., Devereaux, P.J., McGinn, T. and Guyatt, G. (2017) 'Discrimination and calibration of clinical prediction models: users' guides to the medical literature', *JAMA*, 318(14), pp. 1377–1384.

Albert, A. and Anderson, J.A. (1984) 'On the existence of maximum likelihood estimates in logistic regression models', *Biometrika*, 71(1), pp. 1–10.

Balachandran, V.P., Gonen, M., Smith, J.J. and DeMatteo, R.P. (2015) 'Nomograms in oncology: more than meets the eye', *The Lancet Oncology*, 16(4), pp. e173–e180.

Carpenter, J. and Bithell, J. (2000) 'Bootstrap confidence intervals: when, which, what? A practical guide for medical statisticians', *Statistics in Medicine*, 19(9), pp. 1141–1164.

Harrell, F.E. (2015) *Regression Modeling Strategies: With Applications to Linear Models, Logistic and Ordinal Regression, and Survival Analysis*. 2nd edn. New York: Springer.

Harrell, F.E., Lee, K.L. and Mark, D.B. (1996) 'Multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors', *Statistics in Medicine*, 15(4), pp. 361–387.

Hoerl, A.E. and Kennard, R.W. (1970) 'Ridge regression: biased estimation for nonorthogonal problems', *Technometrics*, 12(1), pp. 55–67.

Iasonos, A., Schrag, D., Raj, G.V. and Panageas, K.S. (2008) 'How to build and interpret a nomogram for cancer prognosis', *Journal of Clinical Oncology*, 26(8), pp. 1364–1370.

Kerr, K.F., Brown, M.D., Zhu, K. and Janes, H. (2016) 'Assessing the clinical impact of risk prediction models with decision curves: guidance for correct interpretation and appropriate use', *Journal of Clinical Oncology*, 34(21), pp. 2534–2540.

Steyerberg, E.W. and Harrell, F.E. (2016) 'Prediction models need appropriate internal, internal-external, and external validation', *Journal of Clinical Epidemiology*, 69, pp. 245–247.

Steyerberg, E.W., Harrell, F.E., Borsboom, G.J.J.M., Eijkemans, M.J.C., Vergouwe, Y. and Habbema, J.D.F. (2001) 'Internal validation of predictive models: efficiency of some procedures for logistic regression analysis', *Journal of Clinical Epidemiology*, 54(8), pp. 774–781.

Van Calster, B., McLernon, D.J., van Smeden, M., Wynants, L. and Steyerberg, E.W. (2019) 'Calibration: the Achilles heel of predictive analytics', *BMC Medicine*, 17, 230.

Vickers, A.J. and Elkin, E.B. (2006) 'Decision curve analysis: a novel method for evaluating prediction models', *Medical Decision Making*, 26(6), pp. 565–574.

Vickers, A.J., van Calster, B. and Steyerberg, E.W. (2016) 'Net benefit approaches to the evaluation of prediction models, molecular markers, and diagnostic tests', *BMJ*, 352, i6.
