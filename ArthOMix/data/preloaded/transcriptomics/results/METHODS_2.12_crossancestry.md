# Thesis §2.12 — Cross-ancestry evaluation of the prioritised biomarker genes (METHODS ONLY)

**One of three parallel biomarker methods sections** (§2.9 sex-stratified blood;
§2.11 cross-tissue; §2.12 cross-ancestry, *this file*). Written against
`scripts/goal2_sex_stratified/26_crossancestry_biomarker_mr.R` and
`27_crossancestry_biomarker_figure.R`, with instrument construction inherited from
`scripts/00_shared/10_MR.R`. Software versions from `ENVIRONMENT.txt`. Citation
style: author–date (Harvard), in-text. Contains **no results**: no gene counts, no
odds ratios, no replication tallies.

---

## 2.12 Cross-ancestry evaluation by comparative Mendelian randomisation

### 2.12.1 Question, and the level at which it is answered

The panels of §2.9 rest on genetic evidence obtained in a European-ancestry
population: candidate genes were admitted only if their *cis*-eQTL-instrumented
expression was associated with RA liability in a European GWAS. This section asks
whether that underlying genetic evidence is specific to the discovery cohort, is
reproducible in a second European cohort, or extends to a population of different
ancestry.

The question is answered at the level of the **genetic evidence for the genes**,
not at the level of the classifier. No expression data from any non-European
cohort were available, so the diagnostic model itself is not re-evaluated across
ancestries and no cross-ancestry performance estimate is produced. What is compared
is the Mendelian randomisation (MR) estimate for each prioritised gene against RA
outcome GWAS from three cohorts. This distinction is stated first because it
determines the strongest claim the section can support: evidence that a gene's
inferred effect on RA liability is or is not shared across populations, not
evidence that a blood-expression panel would classify patients of another ancestry.

A second constraint applies throughout. Both the eQTL exposure resource and all
three RA outcome GWAS are sex-combined. The MR estimates are therefore not
sex-specific: the per-sex tables reported here **partition genes** according to the
stratum in which each was prioritised, and a gene prioritised in both strata
carries an identical estimate in both. Any difference between the female and male
cross-ancestry summaries arises from which genes entered each stratum upstream, not
from sex-specific genetic effects.

### 2.12.2 Genes tested and genetic instruments

The genes carried into this analysis were the within-sex prioritised sets defined
in §2.6 — those surviving Benjamini–Hochberg correction at a false discovery rate
below 0.05 in the European discovery MR (Benjamini and Hochberg, 1995) — combined
into their union for extraction and partitioned back by stratum for reporting.

Instruments were the European whole-blood *cis*-eQTL variants used in the discovery
MR and were reused unchanged rather than reselected, so that the three outcome
cohorts are compared on an identical exposure definition. Their construction is
described in full in §2.6 and is summarised here for readability: variants
associated with the expression of the candidate gene at genome-wide significance
in eQTLGen (Võsa et al., 2021), LD-clumped to near independence against a European
reference panel (1000 Genomes Project Consortium, 2015), restricted to variants
lying within one megabase of the gene body on the same chromosome so that the
instrument is genuinely *cis*-acting, and retained only where the instrument
strength statistic exceeded the conventional threshold of ten (Burgess and
Thompson, 2011).

### 2.12.3 Outcome cohorts

Three RA case–control GWAS were used as outcomes, accessed through the MRC IEU
OpenGWAS infrastructure (Elsworth et al., 2020) via the TwoSampleMR interface
(Hemani et al., 2018). The first is the European discovery GWAS used to prioritise
the genes in the first place (Okada et al., 2014), whose estimates were reused from
the discovery pipeline rather than refetched, so that the comparator column is
identical to the one on which selection was based. The second is an independent
European RA GWAS used as a same-ancestry replication cohort (Stahl et al., 2010).
The third is an East Asian RA GWAS from the Biobank Japan resource (Ishigaki et
al., 2020), used as the cross-ancestry comparator.

For the two cohorts not already harmonised in the discovery pipeline, outcome
summary statistics were extracted for the instrument variants and harmonised
against the exposure data with the strand of palindromic variants inferred from
allele frequency, palindromic variants of intermediate frequency being dropped
rather than guessed. Gene labels were re-attached to variants after harmonisation
by explicit matching, since harmonisation may reorder or drop rows.

### 2.12.4 Estimation

For each gene and each outcome, the primary causal estimate was obtained by the
same method ladder used in the discovery analysis, so that no gene is advantaged by
being estimated differently in one cohort than in another. Where three or more
instrument variants survived harmonisation, inverse-variance-weighted MR (Burgess,
Butterworth and Thompson, 2013) was computed together with MR-Egger regression
(Bowden, Davey Smith and Burgess, 2015) and the weighted median estimator (Bowden
et al., 2016). Where exactly two variants survived, only the
inverse-variance-weighted estimate was computed, and where a single variant
survived, the Wald ratio was used. The reported estimate for each gene is taken
from a fixed method priority — inverse-variance weighted, then Wald ratio, then
weighted median, then MR-Egger — applied identically in every cohort, so the choice
of estimator is a function of how many instruments survived and never of which
estimator gave the more favourable answer.

Effect estimates were exponentiated to odds ratios per unit increase in genetically
predicted expression, with ninety-five per cent confidence intervals from the
normal approximation on the log-odds scale, and each gene was classified in each
cohort as risk-increasing, protective, or non-significant.

### 2.12.5 Cross-cohort classification

Two comparisons were made against the European discovery estimate. A gene was
counted as replicated in Europeans when its estimate in the second European cohort
lay on the same side of the null and reached nominal significance. It was counted
as transferable to East Asians when its estimate in the East Asian cohort likewise
agreed in direction and reached nominal significance. Both criteria require
directional agreement as well as a p-value, so that a nominally significant
estimate of opposite sign is not miscounted as support.

The significance threshold used for these two comparisons is nominal rather than
false-discovery-rate corrected, and this differs deliberately from the discovery
stage of §2.7, where correction was applied. The rationale is that replication and
transfer are being assessed for a small, pre-specified set of genes chosen before
these cohorts were examined, rather than screened from the transcriptome; the
threshold is therefore not a discovery threshold. The distinction is stated
explicitly because the two stages use different rules and a reader comparing gene
counts between sections would otherwise be misled.

Genes were assigned to one of four mutually exclusive classes: untestable in East
Asians, where no instrument variant survived harmonisation against that cohort;
shared between the two ancestries; replicated in Europeans but not transferring to
East Asians; and European-discovery only.

### 2.12.6 Instrument transferability diagnostics

The class "untestable" exists because a null result in a cross-ancestry MR has two
quite different explanations that must not be conflated: the gene may have no
causal effect in that population, or the European instrument may simply not be
usable there. Diagnostics were therefore computed to separate them.

For every gene and every outcome cohort, the number of instrument variants
surviving harmonisation was recorded, so that a gene reaching the East Asian
analysis with no usable instrument is reported as untestable rather than as null.
In addition, the mean absolute difference in effect-allele frequency between the
European exposure data and each outcome cohort was computed per gene, giving a
direct measure of how far the instrument's allele frequencies diverge between
populations.

### 2.12.7 The central caveat: ancestry mismatch on the exposure side

The eQTL resource used to instrument gene expression is European. No comparable
East Asian whole-blood *cis*-eQTL resource was available through the same
infrastructure, so the East Asian arm necessarily tests **European-derived
instruments against an East Asian outcome**. The two European arms are
ancestry-matched on both sides; the East Asian arm is matched on neither the linkage
disequilibrium structure nor the allele-frequency spectrum of its exposure.

The consequence is that the East Asian analysis is exploratory and is not a
like-for-like replication, and it is labelled as such wherever it appears. A variant
tagging the causal variant in Europeans may tag it poorly or not at all in East
Asians, and allele frequencies may differ enough that the instrument is
underpowered even where it is present; both mechanisms attenuate an estimate
towards the null without any change in the underlying biology. This is the same
portability problem documented for genetic prediction more generally across
ancestries (Martin et al., 2019). A failure to transfer must therefore be read as
an absence of evidence under mismatched instruments, not as evidence of absence of
a causal effect, and the transferability diagnostics of §2.12.6 are reported
alongside every such gene so that a reader can see which of the two situations
applies. A properly ancestry-matched analysis would require East Asian *cis*-eQTL
data and, ideally, a multi-ancestry RA GWAS with per-ancestry summary statistics
(Ishigaki et al., 2022); this is identified as a limitation and as the natural
extension of the analysis rather than presented as having been done.

### 2.12.8 Software and reproducibility

Analyses were conducted in R v4.4.2 using TwoSampleMR v0.7.8 (Hemani et al., 2018),
dplyr and data.table v1.18.2.1, with figures produced in ggplot2 v4.0.2 (Wickham,
2016). A fixed seed was set at the start of the script. Outputs were written under
a dedicated file prefix so that no result of the discovery MR pipeline is
overwritten by this analysis.

---

## References cited in this section

1000 Genomes Project Consortium (2015) 'A global reference for human genetic
variation', *Nature*, 526(7571), pp. 68–74.

Benjamini, Y. and Hochberg, Y. (1995) 'Controlling the false discovery rate: a
practical and powerful approach to multiple testing', *Journal of the Royal
Statistical Society: Series B*, 57(1), pp. 289–300.

Bowden, J., Davey Smith, G. and Burgess, S. (2015) 'Mendelian randomization with
invalid instruments: effect estimation and bias detection through Egger
regression', *International Journal of Epidemiology*, 44(2), pp. 512–525.

Bowden, J., Davey Smith, G., Haycock, P.C. and Burgess, S. (2016) 'Consistent
estimation in Mendelian randomization with some invalid instruments using a
weighted median estimator', *Genetic Epidemiology*, 40(4), pp. 304–314.

Burgess, S., Butterworth, A. and Thompson, S.G. (2013) 'Mendelian randomization
analysis with multiple genetic variants using summarized data', *Genetic
Epidemiology*, 37(7), pp. 658–665.

Burgess, S. and Thompson, S.G. (2011) 'Avoiding bias from weak instruments in
Mendelian randomization studies', *International Journal of Epidemiology*, 40(3),
pp. 755–764.

Elsworth, B., Lyon, M., Alexander, T., Liu, Y., Matthews, P., Hallett, J., Bates,
P., Palmer, T., Haberland, V., Davey Smith, G., Zheng, J., Haycock, P., Gaunt, T.R.
and Hemani, G. (2020) 'The MRC IEU OpenGWAS data infrastructure', *bioRxiv*.
doi:10.1101/2020.08.10.244293.

Hemani, G., Zheng, J., Elsworth, B., Wade, K.H., Haberland, V., Baird, D.,
Laurin, C., Burgess, S., Bowden, J., Langdon, R., Tan, V.Y., Yarmolinsky, J.,
Shihab, H.A., Timpson, N.J., Evans, D.M., Relton, C., Martin, R.M., Davey Smith,
G., Gaunt, T.R. and Haycock, P.C. (2018) 'The MR-Base platform supports systematic
causal inference across the human phenome', *eLife*, 7, e34408.

Ishigaki, K., Akiyama, M., Kanai, M., Takahashi, A., Kawakami, E., Sugishita, H.,
et al. (2020) 'Large-scale genome-wide association study in a Japanese population
identifies novel susceptibility loci across multiple human diseases', *Nature
Genetics*, 52(7), pp. 669–679.

Ishigaki, K., Sakaue, S., Terao, C., Luo, Y., Sonehara, K., Yamaguchi, K., et al.
(2022) 'Multi-ancestry genome-wide association analyses identify novel genetic
mechanisms in rheumatoid arthritis', *Nature Genetics*, 54(11), pp. 1640–1651.

Martin, A.R., Kanai, M., Kamatani, Y., Okada, Y., Neale, B.M. and Daly, M.J.
(2019) 'Clinical use of current polygenic risk scores may exacerbate health
disparities', *Nature Genetics*, 51(4), pp. 584–591.

Okada, Y., Wu, D., Trynka, G., Raj, T., Terao, C., Ikari, K., et al. (2014)
'Genetics of rheumatoid arthritis contributes to biology and drug discovery',
*Nature*, 506(7488), pp. 376–381.

Stahl, E.A., Raychaudhuri, S., Remmers, E.F., Xie, G., Eyre, S., Thomson, B.P., et
al. (2010) 'Genome-wide association study meta-analysis identifies seven new
rheumatoid arthritis risk loci', *Nature Genetics*, 42(6), pp. 508–514.

Võsa, U., Claringbould, A., Westra, H.-J., Bonder, M.J., Deelen, P., Zeng, B., et
al. (2021) 'Large-scale *cis*- and *trans*-eQTL analyses identify thousands of
genetic loci and polygenic scores that regulate blood gene expression', *Nature
Genetics*, 53(9), pp. 1300–1310.

Wickham, H. (2016) *ggplot2: Elegant Graphics for Data Analysis*. 2nd edn. New
York: Springer.

---

## One citation you must verify before submission

`26_crossancestry_biomarker_mr.R` is internally inconsistent about the Biobank
Japan source: line 10 of its header calls the dataset "BBJ 2019", the closing
reference line calls it "Ishigaki K Nat Genet 2022", and
`scripts/00_shared/REFERENCES.md` lists Ishigaki et al. (2020) *Nature Genetics*
52:669–679 with a "⚠️ confirm" flag. These are different papers: the 2020 paper is
the Biobank Japan multi-disease GWAS, and the 2022 paper is the multi-ancestry RA
GWAS. The section above cites the 2020 paper as the source of the outcome dataset,
because that is the publication conventionally attributed to the `bbj-a-151`
accession, and cites the 2022 paper only where a properly ancestry-matched design
is discussed as future work. **Confirm the attribution against the OpenGWAS record
for that accession before submitting, and correct the script header so it names one
paper rather than three years.**

---

## Points a reader or examiner will press on, and where the text answers them

Line references are to `26_crossancestry_biomarker_mr.R` unless stated.

1. **"Is the cross-ancestry analysis sex-specific?"** No. Both the eQTL exposures
   and all three outcome GWAS are sex-combined, so the per-sex tables partition
   genes rather than estimating sex-specific effects (lines 43-47 partition the
   prioritised sets; the estimation itself is run once on the union). §2.12.1 says
   this in the second paragraph, consistent with the same statement in §2.6.
2. **"Is the East Asian arm a replication?"** No, and the script's own header says
   so (lines 15-25). The exposure instruments are European because no East Asian
   *cis*-eQTL resource was available; §2.12.7 gives the mechanism by which this
   attenuates estimates and why a null is not evidence of absence.
3. **"How is 'untestable' distinguished from 'not causal'?"** By the instrument
   transferability diagnostics — surviving instrument count per cohort and
   EUR-versus-EAS allele-frequency divergence (lines 99-105) — and by making
   `untestable in EAS` its own class rather than folding it into the null class
   (lines 122-126). §2.12.6 describes both.
4. **"Why nominal p < 0.05 here when discovery used BH-FDR?"** Because the gene set
   was fixed before these cohorts were examined, so this is not a screening step
   (lines 118-121 apply `p < 0.05` with a directional requirement). §2.12.5 states
   the difference explicitly rather than leaving two thresholds unexplained across
   sections.
5. **"Was the estimator chosen per gene?"** Only by instrument count, under a fixed
   priority applied identically in every cohort (lines 56-70). §2.12.4 says so.
6. **"Were the discovery estimates refetched?"** No — the European discovery column
   is reused from the primary MR objects (lines 92-95), so the comparator is
   identical to the one selection was based on.
