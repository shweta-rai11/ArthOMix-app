# 2.Z Sex-Stratified Differential Methylation Analysis

## 2.Z.1 Rationale and design

The objective of this stage is to identify CpG sites differentially methylated between rheumatoid arthritis (RA) cases and controls independently within each self-reported sex, so that a female-specific and a male-specific candidate set can each be resolved on their own terms rather than derived from a single pooled model. Sex-stratified EWAS designs of this kind are an established approach for surfacing sex-specific methylation associations that a sex-adjusted, pooled analysis can obscure (Tesfaye et al., 2024). The two strata established in script01 (female n=492, male n=197; Table 2.X.1) are fitted as fully independent models in this stage: no parameter, covariate estimate, or significance threshold is shared between them beyond the modelling choices described below, which are applied identically to preserve comparability.

## 2.Z.2 Data preparation

**Cohort.** The QC-filtered, probe-filtered beta matrix (`beta_raw.rds`, Section 2.X) and the cell-type-augmented phenotype table (`pheno_celltype.rds`, Section 2.Y) were used as input, joined and order-checked by sample identifier as in prior stages.

**Smoking-status exclusion.** Two of the 689 samples carry an unresolved smoking-status value, as already noted in Section 2.X.3's Results. Direct inspection of `pheno_celltype.rds` established that this value is encoded as the literal string `"na"` (lower-case, a residue of the raw GEO field text) rather than as an R `NA` or an empty string. This distinction is recorded here because an initial implementation of the exclusion filter, written as "keep rows where the value is not `NA` and not empty," silently retained both of these samples instead of excluding them, since neither condition matches a literal `"na"` string; the discrepancy was caught only by cross-checking the number of samples the filter actually excluded against the number already established in Section 2.X.3 (two), before this script was run to produce any official output, and the filter was corrected to match the string explicitly, with a hard assertion that exactly two samples are excluded, so that a future re-run silently reverting to the earlier behaviour would fail immediately rather than pass unnoticed. The two affected samples were excluded from the modelled dataset (689 to 687 samples: 492 to 491 female, 197 to 196 male) rather than assigned an invented "unknown" factor level, since a near-singleton category risks rank deficiency or near-perfect separation once the data are further split by sex.

**Cell-type covariates.** Of the seven estimated cell-type fractions from Section 2.Y, six were used as model covariates, with the seventh omitted as an implicit reference category to avoid the collinearity inherent in a compositional variable set (Section 2.Y.6). The omitted category was chosen, programmatically rather than by a fixed a priori choice, as whichever fraction had the largest cohort mean; on this cohort that was `Neutro` (mean fraction 0.699, Table 2.Y.1), so `B`, `NK`, `CD4T`, `CD8T`, `Mono`, and `Eosino` were entered as covariates, identically in both sex-stratified models.

**M-value transformation.** Beta values were converted to M-values, M = log2(beta / (1 - beta)), for use as the model's response variable, following the recommendation of Du et al. (2010) that M-values have statistical properties (approximate homoscedasticity, an unbounded range) better suited to linear-model-based testing than the bounded beta scale; effect sizes are reported separately on the beta scale (as Δbeta, Section 2.Z.3) for biological interpretability, consistent with the same source's recommendation. Because M is undefined at beta = 0 or beta = 1 exactly, beta values were clipped to [1e-6, 1 - 1e-6] immediately beforehand; the quality-controlled, probe-mean-imputed matrix from Section 2.X is not expected to contain exact 0/1 values, but this clip guards against that possibility rather than assuming it.

## 2.Z.3 Sex-stratified models

Within each sex stratum independently, a linear model was fitted to the M-value matrix with `limma::lmFit()`, followed by empirical Bayes moderation (`limma::eBayes()`), using the design `~ group + age + smoking + B + NK + CD4T + CD8T + Mono + Eosino` (`group`: Control as reference, RA as the tested level). This lmFit/eBayes/topTable pattern is the standard approach for probe-wise differential methylation testing on array data (Maksimovic, Phipson & Oshlack, 2016). The moderated t-statistic for the `group` coefficient was extracted for every retained probe and passed to `bacon::bacon()` (van Iterson, van Zwet & BIOS Consortium, 2017) to obtain bias- and inflation-corrected p-values (`bacon::pval(..., corrected = TRUE)`) ahead of Benjamini-Hochberg FDR correction, applied independently within each stratum. Each CpG's descriptive effect size, Δbeta (mean beta in RA minus mean beta in Control, computed directly on the beta scale within that stratum, not from the M-value coefficient), was retained alongside the corrected statistics for panel construction (Section 2.Z.4).

## 2.Z.4 Gene panel construction

Each sex stratum's own bacon-corrected, genome-wide Benjamini-Hochberg FDR-significant (FDR < 0.05) CpGs constitute that sex's panel directly, mapped to an annotated gene via `ChAMPdata::probe.features$gene` (a single nearest/overlapping RefGene symbol per probe, or no gene for an intergenic CpG, consistent with the gene annotation resource already established for probe-level filtering in Section 2.X.6), with each CpG additionally classified as hyper- or hypomethylated in RA relative to Control by the sign of its Δbeta.

## 2.Z.5 Limitations: residual calibration problems in both stratified models

Neither sex-stratified model calibrated cleanly, and this is reported here as an explicit limitation rather than omitted. Bacon's own inflation estimate was high for the female-stratum model and moderately elevated for the male-stratum model, against an ideal value of 1 (Section 2.Z.6); visual inspection of the bacon-corrected QQ plot for each stratified model shows the observed distribution of -log10 p-values still bowing visibly above the expected null line across most of the distribution after correction, more severely for the female stratum. A bacon-corrected QQ plot that remains bowed above the diagonal indicates that bacon's mean/variance empirical-null correction (van Iterson et al., 2017) did not fully resolve the departure from the expected null distribution for either stratified model. Both stratified volcano plots additionally show an unusual bimodal ("double-peak") shape, with the great majority of CpGs distributed into two symmetric plumes at small negative and small positive Δbeta rather than concentrated near Δbeta = 0, consistent with (though not, on its own, proof of) a genome-wide, low-magnitude, systematically directional effect rather than the sparse, near-zero-centred pattern expected if only a small minority of CpGs carried a true RA-associated effect. This chapter does not adjudicate here between the candidate explanations this pattern is consistent with, namely a genuinely widespread, low-magnitude methylation response to a systemic inflammatory disease process, or an unmodelled residual technical confound (for instance array batch, which Section 2.X.1 already establishes cannot be recovered from the processed series-matrix data used throughout this chapter). This calibration problem, and an approach to resolving it, is addressed directly in the following stage of this chapter.

---

## 2.Z.6 Results

**Model fitting.** Both models fitted successfully on the 687-sample smoking-resolved cohort (491 female, 196 male). Table 2.Z.1 summarises the bacon inflation estimate and the number of Benjamini-Hochberg FDR-significant CpGs (bacon-corrected, FDR < 0.05) recovered by each stratum.

**Table 2.Z.1.** Model calibration and genome-wide-significant CpG counts, by sex.

| Stratum | n | bacon lambda | FDR<0.05 CpGs |
|---|---|---|---|
| Female (RA vs Control) | 491 | 2.285 | 0 |
| Male (RA vs Control) | 196 | 1.189 | 0 |

Neither stratum retained any CpG at genome-wide FDR-significance once bacon correction was applied: the resulting female and male gene panels are both empty at this stage. Given the residual calibration problem documented in Section 2.Z.5, an empty panel here is not interpreted as evidence of an absent sex-specific methylation effect; it reflects unresolved genomic inflation that a subsequent, methodologically distinct stage of this chapter addresses directly.

---

## References

Du, P., Zhang, X., Huang, C. C., Jafari, N., Kibbe, W. A., Hou, L., & Lin, S. M. (2010). Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis. *BMC Bioinformatics*, 11, 587.

Maksimovic, J., Phipson, B., & Oshlack, A. (2016). A cross-package Bioconductor workflow for analysing methylation array data. *F1000Research*, 5, 1281.

Tesfaye, M., Spindola, L. M., Stavrum, A. K., Shadrin, A., Melle, I., Andreassen, O. A., & Le Hellard, S. (2024). Sex effects on DNA methylation affect discovery in epigenome-wide association study of schizophrenia. *Molecular Psychiatry*, 29, 2467-2477. PMID: 38503926.

van Iterson, M., van Zwet, E. W., & BIOS Consortium (2017). Controlling bias and inflation in epigenome- and transcriptome-wide association studies using the empirical null distribution. *Genome Biology*, 18, 19.

---
