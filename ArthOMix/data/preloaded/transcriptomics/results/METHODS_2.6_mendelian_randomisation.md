# 2.6 Instrument selection and Mendelian randomisation design

*(Replacement prose for thesis §2.6. Every number is taken from the current saved
analysis objects, not from the run logs — `results/logs/10_MR.log` pre-dates the cis
filter and reports superseded counts. Companion document: `results/METHODS_MR.md`.)*

---

## 2.6.1 Rationale and instrumental-variable assumptions

Differential expression analysis identifies genes whose transcript abundance differs between
cases and controls, but it cannot separate a driver of disease from a downstream consequence of
chronic inflammation or of immunosuppressive treatment. Two-sample cis-Mendelian randomisation
(MR) was therefore used as a causal-anchoring step. Germline variants that regulate a gene's
expression in cis (expression quantitative trait loci, eQTLs) serve as instrumental variables for
that gene's expression level, and their association with rheumatoid arthritis (RA) risk is tested.
Because genotype is fixed at conception and assorted randomly at meiosis, an estimate obtained
this way is not subject to reverse causation, nor to the environmental, disease-activity and
treatment confounders that bias observational expression–disease associations.

MR was positioned between differential-expression screening and machine-learning feature
selection: the disease-associated differentially expressed genes were the exposures tested, and
only genes surviving MR at a false-discovery rate below 0.05 were passed forward as the
feature-selection input.

The analysis rests on the three standard instrumental-variable assumptions: **relevance**, that
each instrument is robustly associated with the exposure; **independence**, that no instrument
shares a common cause with the outcome; and the **exclusion restriction**, that each instrument
influences RA only through the expression of the gene it instruments. Section 2.6.6 states which
of these assumptions were formally tested, and for how many genes.

---

## 2.6.2 Data sources

Instruments were drawn from the eQTLGen Consortium whole-blood cis-eQTL resource (Võsa et al.,
2021), accessed through the IEU OpenGWAS database as datasets `eqtl-a-<ENSG>`. Whole blood was
chosen deliberately so that the instrumented exposure is the same quantity measured by the
transcriptomic discovery cohorts and by the resulting biomarker panels.

The primary outcome was the European RA genome-wide association study of Okada et al. (2014),
OpenGWAS `ieu-a-832`. Two further outcomes were used for replication and cross-ancestry
assessment (§2.6.9): Stahl et al. (2010), European, `ieu-a-834`; and the Biobank Japan RA GWAS,
East Asian, `bbj-a-151`.

HGNC gene symbols were mapped to Ensembl gene identifiers using `org.Hs.eg.db`. Only publicly
available summary statistics were used; no individual-level genotype data were accessed, and the
work therefore required no ethical approval beyond that held by the contributing studies.

---

## 2.6.3 Instrument selection

Candidate exposures were the disease-associated differentially expressed genes carried forward
from the sex-stratified screen: 2,045 female and 2,079 male candidates, 2,351 unique. MR was
executed once on the union of both lists and the results partitioned by stratum afterwards
(§2.6.7), so that the two strata differ only in which genes are eligible and never in how a gene
is estimated.

Instruments were extracted with `TwoSampleMR::extract_instruments()` under four criteria.

**Relevance.** Association with expression at p < 5 × 10⁻⁸.

**Independence.** Linkage-disequilibrium clumping at r² < 0.001 within a 10,000 kb window, using
the European 1000 Genomes reference panel, to retain approximately independent signals.

**Instrument strength.** The F-statistic, F = (β_exposure / SE_exposure)², was computed for every
instrument. It is important to be precise about the role of this quantity. At p < 5 × 10⁻⁸ the
absolute Z-score is at least 5.451, so F = Z² is at least 29.7 for every retained instrument **by
construction**; the conventional F ≥ 10 criterion removed no variants and cannot remove any at
this significance threshold. F is therefore reported as a **verification** that instrument
strength comfortably exceeds the weak-instrument boundary, not as a selection step. The observed
distribution confirms this: minimum F = 29.7, median F = 119.0, maximum F = 40,302.

**Cis restriction.** Each instrument was required to lie on the same chromosome as the gene it
instruments and within 1 Mb of the gene body, using GRCh37 coordinates from `EnsDb.Hsapiens.v75`.
The GRCh37 annotation was chosen to match the genome build of the eQTLGen and OpenGWAS variant
positions; a GRCh38 annotation would render the coordinates non-comparable and corrupt the filter
silently.

This restriction requires explicit statement because it is not supplied by the data source.
eQTLGen datasets held in OpenGWAS contain trans as well as cis associations, and
`extract_instruments()` returns both; a cis-MR design that does not enforce the window is
therefore not, in fact, a cis design. Applying the filter reduced the instrument set from 4,932
variants across 1,980 genes to **3,897 variants across 1,947 genes**. Of the 1,035 variants
removed, 996 — twenty per cent of the original set — were located on a *different chromosome* from
the gene they nominally instrumented; the remainder were same-chromosome variants lying more than
1 Mb from the gene body, or belonged to genes without a GRCh37 coordinate.

The filter was consequential rather than cosmetic. Before it was applied, the two largest effect
estimates in the entire analysis arose from trans instruments situated on the two strongest RA
loci in the genome: *HNRNPM* on chromosome 19, instrumented by chr6:32,431,962 in the MHC class II
region (OR 287.3), and *FOXP3* on chromosome X, instrumented by chr1:114,303,808 at the *PTPN22*
locus (OR 263.9). In both cases the variant reaches RA through HLA-DRB1 or PTPN22 rather than
through the expression of the instrumented gene, so the exclusion restriction is violated by
construction rather than merely untested.

Instruments falling within the extended MHC (chromosome 6, 25–34 Mb, GRCh37) were flagged but
retained in the primary analysis — 84 variants across 33 genes — and removed entirely in a full
parallel sensitivity analysis described in §2.6.8.

---

## 2.6.4 Outcome retrieval and harmonisation

Outcome associations were retrieved with `extract_outcome_data()` with LD proxy substitution
enabled, so that instruments absent from the outcome GWAS could be represented by variants in
strong linkage disequilibrium rather than discarded. Requests were issued in batches to avoid
query truncation.

Exposure and outcome datasets were harmonised with `harmonise_data(action = 2)`. For palindromic
variants, whose alleles are uninformative about strand, this procedure infers the strand from
allele frequency and discards those whose frequencies remain too close to 0.5 for the inference to
be reliable; non-palindromic variants are resolved directly from their alleles and require no
frequency-based inference. Harmonisation produced 3,081 exposure–outcome records across 1,844
genes, of which 2,703 were retained for estimation.

Gene labels were re-attached to harmonised records by matching on the pair of variant identifier
**and** exposure dataset identifier, rather than on the variant alone. Because a single variant may
instrument more than one gene, matching on the variant alone assigns every record for that variant
to whichever gene happens to appear first, silently transferring instruments between genes; the
exposure dataset identifier `eqtl-a-<ENSG>` is what actually identifies the gene.

---

## 2.6.5 Effect estimation

MR was estimated separately for each gene, with the set of estimators determined by the number of
instruments surviving harmonisation. Genes with three or more instruments were analysed by
inverse-variance weighted (IVW) regression together with MR-Egger and the weighted median
estimator; genes with two instruments by IVW alone; and genes with a single instrument by the Wald
ratio.

Where more than one estimator was fitted, a single primary estimate per gene was selected by the
**pre-specified** hierarchy IVW > Wald ratio > weighted median > MR-Egger. Fixing this order in
advance ensures that the reported estimate is never selected on the basis of its p-value. Effects
are reported as odds ratios for RA per unit increase in normalised gene expression — not per risk
allele — with 95% confidence intervals given by exp(β ± 1.96 SE).

Estimates were obtained for 1,694 genes. The distribution of instrument counts is strongly skewed
towards single-variant support: 1,000 genes had one instrument, 458 had two, 179 had three, 49 had
four and 8 had five or more. Consequently the primary estimate is a Wald ratio for 1,000 genes and
an IVW estimate for 694.

---

## 2.6.6 Sensitivity analyses and assumption testing

For the 236 genes instrumented by three or more variants, two diagnostics were computed:
**Cochran's Q**, quantifying heterogeneity among instrument-specific estimates as an indirect
indicator of pleiotropy; and the **MR-Egger intercept**, for which a non-zero value indicates
directional horizontal pleiotropy. Funnel plots, leave-one-out IVW re-estimation and per-variant
Wald ratios were generated for the same subset, together with an explicit availability table
recording, for each prioritised gene, which diagnostics were possible.

This coverage must be stated plainly rather than implied. Neither diagnostic is defined for a gene
instrumented by one or two variants, and that describes the majority of genes analysed here. For
those genes **the exclusion restriction is untested, not satisfied**. Every results table
accordingly carries both the instrument count and a `sensitivity_testable` flag, so that no
estimate can be read without knowing whether any pleiotropy assessment was possible for it.

| Assumption | How addressed | Coverage |
|---|---|---|
| Relevance | p < 5 × 10⁻⁸; F verified | All 3,897 instruments (min F 29.7) |
| Independence | Clumping r² < 0.001 / 10 Mb; two-sample design | All genes |
| Exclusion restriction | Cis restriction ±1 Mb; MR-Egger intercept; Cochran's Q; MHC re-run; colocalisation (§2.7) | Cis filter: all genes. Formal tests: 236 of 1,694 |
| Population stratification | European exposure matched to European outcome | Primary and Stahl arms; the East Asian arm is deliberately mismatched and reported as exploratory |

A further limitation of cis-eQTL MR is that it establishes association between an eQTL and disease
but cannot, on its own, establish that the eQTL signal and the disease signal arise from the
**same** causal variant rather than from two distinct variants in linkage disequilibrium. This was
addressed by Bayesian colocalisation, described in §2.7.

---

## 2.6.7 Multiple-testing correction and the sex strata

P-values were corrected using the Benjamini–Hochberg procedure at a false-discovery rate of 0.05,
applied **within each sex stratum**: the correction denominator is the set of genes actually tested
in that stratum — 1,477 female and 1,478 male — rather than the union of both. In a sex-stratified
design the relevant family of hypotheses is the set of genes eligible in that stratum, and pooling
would inflate the denominator with genes that were never tested in that sex. A pooled correction is
retained in the output files as an unreported legacy column.

The FDR-surviving set, and not the nominal screen, was passed to feature selection. This
distinction is material to everything downstream: the nominal p < 0.05 screen returned 113 female
and 115 male genes against approximately 74.5 and 74.7 expected false positives respectively —
roughly two-thirds noise — and had that list been carried forward it would have driven all
subsequent model building. After correction, **32 female and 25 male genes** were prioritised.

### A note on sex

Both eQTLGen and all three RA GWAS used here are **sex-combined**, and a survey of all 37 RA
datasets available in OpenGWAS confirmed that no sex-stratified RA GWAS exists. The MR estimate for
a given gene is therefore numerically **identical** in the female and male outputs; the two strata
differ only in which genes were eligible for testing and in the correction denominator, which are
near-identical at 1,477 versus 1,478.

The design is consequently one of **sex-stratified discovery followed by sex-combined genetic
validation**, and no sex-differential genetic claim is supported by the MR. It follows that the
difference between the two prioritised lists is not evidence of a sex difference: it reflects only
upstream eligibility, which is itself **sex-stratified output rather than a sex-differential
finding** — the within-sex differential-expression contrasts were computed separately in each
stratum and contain no cross-sex comparison. Lists are accordingly labelled "female-list-only" and
"male-list-only" rather than "female-specific" and "male-specific".

The point is worth stating plainly because it is easily inverted. **Nothing in the genetic arm of
this work is sex-differential**, and the stratification here is organisational, not evidential. The
only analysis in the thesis that estimates a difference between the sexes is the diagnosis-by-sex
interaction test of §2.10, and it is that analysis — not this one — which supplies the estimand
comparable to the interaction design used in the methylomic chapter.

---

## 2.6.8 MHC sensitivity analysis

The strongest objection to any causal claim arising from these data concerns the extended MHC, for
two reasons. HLA-DRB1 is by a wide margin the dominant RA susceptibility locus; and the extended
MHC carries the most extensive long-range linkage disequilibrium in the human genome. A cis-eQTL
for *any* gene in that region is therefore correlated with the HLA-DRB1 risk haplotype whether or
not the gene is causal, and a cis-MR estimate for an MHC gene cannot be distinguished from the
HLA-DRB1 signal read through a proxy. This is not a marginal concern here: 14 of the 32 female
(44%) and 10 of the 25 male (40%) FDR-surviving genes are MHC-flagged, most of them estimated by a
single-variant Wald ratio for which no pleiotropy test is possible.

The entire MR was therefore re-run with MHC instruments removed and reported as a full parallel
column against the primary analysis. The instrument set, the estimator hierarchy, the
heterogeneity and pleiotropy tests, and the per-stratum correction were all re-computed — the
estimator hierarchy in particular must be re-applied, because removing variants can demote a gene
from IVW to a Wald ratio. Instrument extraction, clumping, cis filtering and harmonisation were not
repeated; the sensitivity analysis reads the primary run's own cached harmonised data, so it is
fully deterministic and cannot drift from the analysis it is testing.

The correction denominator legitimately shrinks under this analysis, from 1,477 to 1,448 in
females and from 1,478 to 1,456 in males. A gene whose only instrument lay in the MHC has not been
tested and found null; it has become **untestable**, and retaining it in the denominator would
penalise testable genes for an absence of evidence about untestable ones. The shrinkage is reported
explicitly rather than absorbed silently.

Each gene received one of five verdicts. The taxonomy distinguishes two mechanisms that would
otherwise be conflated: genuine MHC dependence, in which the gene's own instrument set changed and
its estimate is therefore different; and Benjamini–Hochberg re-ranking, in which the estimate is
bit-identical but removing the strongly significant MHC genes ranked above it raises its adjusted
p-value. The first is evidence of confounding; the second is multiplicity bookkeeping, and treating
them alike would overstate the damage.

| Verdict | Meaning |
|---|---|
| ROBUST | FDR < 0.05 both with and without MHC instruments |
| UNTESTABLE without MHC | No non-MHC instrument remains |
| MHC-DEPENDENT | Instrument set changed and significance was lost |
| FDR-RANK ONLY | Estimate identical; significance lost solely to re-ranking |
| ns in both | Not significant under either analysis |

Female prioritised genes fell from 32 to 14 surviving (14 ROBUST, 14 UNTESTABLE, 0 MHC-DEPENDENT,
4 FDR-RANK ONLY); male genes from 25 to 14 (14 ROBUST, 10 UNTESTABLE, 0 MHC-DEPENDENT, 1 FDR-RANK
ONLY). Effect direction agreed between the two analyses for every gene estimable in both
(1,448 of 1,448 female; 1,456 of 1,456 male). Of the six genes in each final biomarker panel, three
are ROBUST, two UNTESTABLE and one FDR-RANK ONLY, in both sexes.

**The reporting rule adopted throughout this thesis follows directly.** A gene may carry a *causal*
claim only where its verdict is ROBUST. Genes classified MHC-DEPENDENT or UNTESTABLE are described
as **associated, not causal**, and are reported as MHC-confounded. No gene is silently discarded;
the fate of every prioritised gene and every panel gene is tabulated.

---

## 2.6.9 Replication and cross-ancestry transferability

Each prioritised gene was re-tested against two further RA GWAS using the identical estimator
hierarchy: Stahl et al. (2010, European) as an ancestry-matched replication, and Biobank Japan
(East Asian) as an exploratory cross-ancestry arm. A gene was counted as replicated or transferable
where the outcome-specific estimate was nominally significant (p < 0.05) and lay in the same
direction relative to the null as the primary estimate.

The East Asian arm carries an important caveat. OpenGWAS holds European eQTL data only, and no East
Asian cis-eQTL resource is available there, so the *same European instruments* are necessarily
tested against the Japanese outcome. The exposure side is thus ancestry-mismatched, and this arm is
exploratory rather than a like-for-like replication. So that a null result remains interpretable,
instrument-transferability diagnostics are reported alongside the estimates — the number of European
instruments surviving harmonisation in each outcome, and the mean European-versus-East-Asian allele
frequency divergence per gene — allowing a gene to be classified as *untestable in East Asians*,
because its instruments are absent or non-polymorphic, as distinct from *not causal in East Asians*.

Of the 32 female prioritised genes, 23 replicated in the European replication cohort, 7 were
transferable to the East Asian cohort, 2 were untestable there and 5 were European-discovery only.
Of the 25 male genes, 18 replicated, 6 were transferable, 2 untestable and 3 European-discovery
only.

---

## 2.6.10 Gene flow through the analysis

2,351 unique candidate genes → 1,947 with at least one cis instrument → 1,694 with an MR estimate →
32 female and 25 male surviving FDR < 0.05 → 14 female and 14 male classified ROBUST after MHC
exclusion.

---

## 2.6.11 Software

R 4.4.2 with `TwoSampleMR` 0.7.8, `ieugwasr` 1.1.0.9000, `dplyr` 1.2.0, `data.table` 1.18.2.1,
`EnsDb.Hsapiens.v75` 2.99.0 and `org.Hs.eg.db` 3.20.0. A fixed seed (2024) was set in every script;
note that two-sample MR involves no stochastic component, so this is a reproducibility convention
rather than a component of the design. Instruments and harmonised data were cached, and the MHC
sensitivity and cross-ancestry analyses read from those caches, so both reproduce exactly without
re-querying OpenGWAS.

---

## Limitations to state explicitly

1. No sex-stratified RA GWAS exists, so the causal arm is sex-combined and contributes no
   sex-specific evidence.
2. Most genes are instrumented by a single cis-eQTL, so pleiotropy and heterogeneity are formally
   untestable for them; the exclusion restriction is untested rather than satisfied.
3. MHC confounding cannot be resolved by MR alone; approximately 40% of prioritised genes lose
   testability once MHC instruments are removed and are reported as associated rather than causal.
4. Whole-blood eQTLs may not represent regulation in synovium, the tissue of primary pathology.
5. European-only exposure data restricts the cross-ancestry analysis to an exploratory
   interpretation.
6. Estimates reflect lifelong, genetically determined differences in expression, which are
   typically small and not equivalent in magnitude to a therapeutic intervention on the same gene.
7. Reverse MR was not performed; the design is directional by construction (cis-eQTL → disease).
8. MR-PRESSO was not applied because it requires at least four instruments to detect outliers,
   a threshold met by only 57 of 1,694 genes. Multivariable MR was not used for the MHC because the
   region's long-range linkage disequilibrium renders the instruments non-independent, violating
   that method's own assumptions; MHC confounding was instead addressed by instrument excision
   (§2.6.8) and by multiple-causal-variant colocalisation (§2.7).
