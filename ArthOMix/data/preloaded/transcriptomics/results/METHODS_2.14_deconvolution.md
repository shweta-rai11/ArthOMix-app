# 2.14 Immune deconvolution

To characterise the immune-cell composition underlying the transcriptomic signal,
the whole-blood expression matrices were deconvolved using CIBERSORT with the
LM22 signature matrix, which resolves 22 human immune-cell subsets, as
implemented in the IOBR package (`deconvo_tme`, method = "cibersort") (Newman et
al., 2015; Zeng et al., 2021). LM22 was derived from microarray data and is
therefore the appropriate reference for the Affymetrix cohorts analysed here.
Deconvolution was applied independently to the training cohort, the sealed
internal-validation holdout and the external blood cohort, and was deliberately
not fitted jointly across them: the holdout is a sealed test set, and a jointly
estimated composition would allow its samples to influence a quantity
subsequently used to adjust the training model. Since CIBERSORT requires a
non-logarithmic mixture, the ComBat batch-corrected log₂ expression matrix was
first converted back to the linear scale (2^expression) prior to deconvolution,
any non-finite value being set to zero, and the analysis was conducted in
microarray mode (arrays = TRUE), which applies the quantile-normalised handling of
the mixture appropriate for array data (Bolstad et al., 2003). The goodness-of-fit
for deconvolution was evaluated through 100 permutations (perm = 100), providing a
per-sample empirical p-value along with the CIBERSORT correlation and
root-mean-square-error diagnostics, and a fixed seed (`set.seed(1234)`) was
established for reproducibility. These three diagnostics were retained for quality
assessment only and were excluded from all downstream modelling, in which only the
22 LM22 fraction columns were used. As a methodologically independent check,
MCP-counter was additionally run on the training and external cohorts (Becht et
al., 2016); because it returns relative abundance scores rather than fractions, it
was used solely to corroborate the direction of compositional differences and
never as the basis of any adjustment, on the principle that agreement between two
estimators that fail in different ways is more informative than either alone
(Shen-Orr and Gaujoux, 2013). The output consisted of a table of the 22 estimated
cell-type fractions for each sample in each dataset, which was integrated with the
sample-level sex and disease-group annotations for subsequent analysis.

To examine the relationship between immune composition and the disease signature,
the estimated cell-type fractions in the training cohort were compared between
rheumatoid arthritis (RA) cases and healthy controls within each sex, the two
strata being tested separately so that no comparison pooled female and male
samples. For each of the 22 cell types, the difference in fraction between RA and
control samples was assessed using a Wilcoxon rank-sum test (Wilcoxon, 1945; Mann
and Whitney, 1947), a distribution-free test being preferred because cell
fractions are bounded, frequently zero-inflated and not plausibly Gaussian. Cell
types with zero variance within a stratum were skipped rather than tested. The
resulting p-values were corrected for multiple testing across the 22 subsets by
the Benjamini–Hochberg procedure, applied separately within each sex, and a subset
was declared to differ by disease status at a false discovery rate below 0.05
(Benjamini and Hochberg, 1995). This test determines whether composition is a live
confounder of, or mediator on, the differential-expression signal: were the
fractions not to differ, composition could not be driving that signal and any
adjustment for it would be a formality, whereas a difference requires the
unadjusted expression results to be reported as confounded with, or mediated by,
composition.

Because fractions are compositional and sum to one, they were not entered
individually into any adjustment model. Subsets whose mean fraction fell below
0.001 were dropped, remaining zeros were replaced by half the smallest observed
positive fraction, and the retained fractions were transformed to centred
log-ratio coordinates (Aitchison, 1982). Principal component analysis was applied
to the centred log-ratio matrix without further scaling, and the first three
components were retained as a low-dimensional summary of composition for use as
covariates in the composition-adjusted analyses.

---

## Removed from the draft as not present in the code

- The comparison of cell-type fractions **between female and male RA patients**.
  No script performs it; the only Wilcoxon test on LM22 fractions is RA versus HC
  within each sex ([05c_deconvolution.R:227-244](scripts/00_shared/05c_deconvolution.R#L227-L244)).
- The **uncorrected significance stars** (p ≤ 0.05, 0.01, 0.001, 0.0001). That
  scheme is defined in [30_figure_fig4A_expression_groups.R:11,46-47](scripts/goal2_sex_stratified/30_figure_fig4A_expression_groups.R#L11-L47)
  and applies to gene expression, not cell fractions. The fraction tests use
  Benjamini–Hochberg correction within sex.
- The **per-sample stacked bar charts** of the 22 fractions and the **per-cell-type
  female-versus-male boxplots**. Both were deliberately deleted from the pipeline:
  [33_figure_fig4_composite.R:8](scripts/goal2_sex_stratified/33_figure_fig4_composite.R#L8)
  records "former panels D and E (CIBERSORT stacked bar + boxplot) were REMOVED".

Two code-level notes, unrelated to the text above: this section and the
cross-ancestry section are both numbered 2.14, and 05c decides whether to log the
external cohort on the 99th percentile > 100
([05c:126](scripts/00_shared/05c_deconvolution.R#L126)) while scripts 16b and 17
use maximum > 50 ([17:84](scripts/goal2_sex_stratified/17_testing_blood_internal_external.R#L84)).

---

## Parameters used

| Item | Value |
|---|---|
| Method | CIBERSORT, LM22 signature (22 immune subsets), via IOBR `deconvo_tme` |
| Mode | `arrays = TRUE` (microarray mixture, quantile-normalised) |
| Input scale | linear, obtained as 2^(ComBat-corrected log₂ matrix); non-finite → 0 |
| Datasets | training, sealed internal holdout, external blood; each deconvolved independently |
| Permutations | `perm = 100` (per-sample empirical p-value) |
| Retained diagnostics | CIBERSORT p-value, correlation, RMSE; quality assessment only, never modelled |
| Modelled columns | the 22 LM22 fraction columns only |
| Second estimator | MCP-counter (training and external), direction corroboration only |
| Seed | `set.seed(1234)` |
| Group comparison | Wilcoxon rank-sum, RA vs HC, within each sex separately |
| Zero-variance subsets | skipped, not tested |
| Multiple testing | Benjamini–Hochberg across the 22 subsets, within each sex; FDR < 0.05 |
| Compositional handling | mean fraction < 0.001 dropped; zeros → half the smallest positive value; centred log-ratio; `prcomp` with `center = TRUE`, `scale. = FALSE`; first 3 PCs retained |
| Software | R with IOBR, limma and data.table |

---

## References cited in this section

Aitchison, J. (1982) 'The statistical analysis of compositional data', *Journal of the Royal Statistical Society: Series B*, 44(2), pp. 139–177.

Becht, E., Giraldo, N.A., Lacroix, L., Buttard, B., Elarouci, N., Petitprez, F., et al. (2016) 'Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression', *Genome Biology*, 17, 218.

Benjamini, Y. and Hochberg, Y. (1995) 'Controlling the false discovery rate: a practical and powerful approach to multiple testing', *Journal of the Royal Statistical Society: Series B*, 57(1), pp. 289–300.

Bolstad, B.M., Irizarry, R.A., Åstrand, M. and Speed, T.P. (2003) 'A comparison of normalization methods for high density oligonucleotide array data based on variance and bias', *Bioinformatics*, 19(2), pp. 185–193.

Mann, H.B. and Whitney, D.R. (1947) 'On a test of whether one of two random variables is stochastically larger than the other', *Annals of Mathematical Statistics*, 18(1), pp. 50–60.

Newman, A.M., Liu, C.L., Green, M.R., Gentles, A.J., Feng, W., Xu, Y., Hoang, C.D., Diehn, M. and Alizadeh, A.A. (2015) 'Robust enumeration of cell subsets from tissue expression profiles', *Nature Methods*, 12(5), pp. 453–457.

Shen-Orr, S.S. and Gaujoux, R. (2013) 'Computational deconvolution: extracting cell type-specific information from heterogeneous samples', *Current Opinion in Immunology*, 25(5), pp. 571–578.

Wilcoxon, F. (1945) 'Individual comparisons by ranking methods', *Biometrics Bulletin*, 1(6), pp. 80–83.

Zeng, D., Ye, Z., Shen, R., Yu, G., Wu, J., Xiong, Y., et al. (2021) 'IOBR: multi-omics immuno-oncology biological research to decode tumor microenvironment and signatures', *Frontiers in Immunology*, 12, 687975.
