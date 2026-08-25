# Thesis §2.7 — Bayesian colocalisation (METHODS ONLY)

Written against `scripts/00_shared/10d_coloc_panel_genes.R` (single causal variant)
and `scripts/00_shared/10e_coloc_susie_mhc.R` (multiple causal variants). Software
versions from `ENVIRONMENT.txt`. Citation style: author–date (Harvard), in-text
only. Contains **no results**: no posterior probabilities, no per-gene verdicts, no
counts.

Section numbering follows the canonical scheme of `METHODS_00_INDEX.md`. This
section was previously carried as Section 10 of the combined document
`METHODS_MR.md`, and as "§2.8" in one superseded draft; both are archived.

---

## 2.7 Bayesian colocalisation

### 2.7.1 Why colocalisation is required

The *cis*-Mendelian randomisation of §2.6 establishes that the expression quantitative
trait locus (eQTL) for a gene is *associated* with rheumatoid arthritis (RA). It does
not establish that the eQTL signal and the RA risk signal are driven by the **same**
causal variant. The two questions come apart whenever the eQTL variant is merely in
linkage disequilibrium with a distinct disease-causing variant, and in that situation
*cis*-MR returns a confident, highly significant and entirely spurious estimate. This
is the principal recognised weakness of *cis*-eQTL MR (Zhu et al., 2016; Wallace,
2020), and it is the reason MR estimates within the major histocompatibility complex
(MHC) cannot be taken at face value.

Colocalisation is therefore not an optional addition to the causal argument but the
step that converts the assumption *cis*-MR makes silently into a quantity that is
measured. It is positioned here, between the Mendelian randomisation of §2.6 and the
feature selection of §2.8, because its verdict governs the language in which every
prioritised gene may subsequently be described.

### 2.7.2 `coloc.abf` — the single-causal-variant model

Bayesian colocalisation was performed with `coloc::coloc.abf` (Giambartolomei et al.,
2014) for every gene surviving false-discovery-rate correction in §2.6 and for every
gene entering the biomarker panels of §2.8. The method evaluates five mutually
exclusive hypotheses over a genomic region: H0, no causal variant for either trait;
H1, a causal variant for expression only; H2, a causal variant for RA only; H3, both
traits have a causal variant but **different** variants; and H4, both traits share
**one** causal variant. H3 is the failure mode that invalidates a *cis*-MR estimate,
and H4 is the state of the world that the MR estimate assumes.

**Input data.** Full *regional* summary statistics were used rather than the clumped
instruments of §2.6, because colocalisation requires all variants in the window; this
step therefore queries the summary-statistic server directly rather than reading the
cached instrument file. The window was the gene body ± 250 kb in GRCh37 coordinates
taken from `EnsDb.Hsapiens.v75`, the same coordinate source as the *cis* filter of
§2.6, so that the two analyses cannot disagree about where a gene lies. A window
narrower than the 1 Mb *cis* boundary used for instrument selection was chosen
deliberately: `coloc.abf` assumes at most one causal variant per trait within the
region analysed, and that assumption becomes less tenable as the window widens to
encompass additional independent signals (Giambartolomei et al., 2014; Wallace,
2021). All regional extracts are cached, so re-runs are deterministic and reproduce
without network access.

**Harmonisation**, deliberately conservative. Variants were merged on rsID and only
those present in both datasets were retained. The genome-wide association effect
allele was aligned to the eQTL effect allele, with sign inversion of β where the
alleles were swapped. Strand-ambiguous variants (A/T and C/G) were **dropped**
outright, because the outcome extract carries no allele frequency, strand therefore
cannot be resolved from minor-allele frequency, and an incorrect flip would be
silent. Variants whose allele pairs did not match after flipping — multi-allelic
sites or annotation mismatches — were dropped. Minor-allele frequency was taken from
the eQTL dataset, both studies being of European ancestry. Regions containing fewer
than 50 shared variants were not analysed at all, because with few variants the
posterior is dominated by the prior rather than by the data and the resulting
probabilities are not interpretable (Wallace, 2020).

**Prior sensitivity.** The p12 prior — the probability that a variant is causal for
both traits — is the principal subjective input to the method and the usual point of
attack upon it (Wallace, 2020). Every gene was therefore analysed twice, at the
package default p12 = 1 × 10⁻⁵ and at a conservative p12 = 1 × 10⁻⁶, and both
posteriors are reported. A colocalisation appearing only at the permissive prior is
**not** reported as support.

**Interpretation thresholds**, conventional (Giambartolomei et al., 2014):

| Criterion | Reading |
|---|---|
| PP.H4 ≥ 0.80 | Strong support for a shared causal variant |
| 0.50 ≤ PP.H4 < 0.80 | Suggestive |
| PP.H3 ≥ 0.80 | Positive evidence **against** colocalisation: distinct causal variants, so the MR estimate is confounded by linkage disequilibrium |
| PP.H4 / (PP.H3 + PP.H4) ≥ 0.80 | Conditional support given that both traits are associated; informative where the region is underpowered and PP.H0–H2 absorb most of the posterior mass |

These thresholds are conventions and not error rates. A gene marginally below a
threshold is not meaningfully different from one marginally above it, so the full
posterior distribution is reported for every gene and no reader is confined to the
dichotomy.

**Where the single-causal-variant assumption breaks, and what follows.**
`coloc.abf` assumes at most one causal variant per trait per region. That is
defensible for a typical *cis*-eQTL window and is **false in the MHC**, where RA has
multiple independent, well-mapped signals — *HLA-DRB1* positions 11, 71 and 74,
*HLA-B* position 9 and *HLA-DPB1* position 9 (Raychaudhuri et al., 2012) — and where
regions carry large numbers of variants in extended linkage disequilibrium (de Bakker
et al., 2006). Violating the assumption **inflates PP.H3**: two traits each driven by
several variants resemble "different causal variants" even when they share one. A
high PP.H3 inside the MHC therefore cannot be read as clean evidence of distinct
causal variants, and — this is the point that must not be softened — the method is
unreliable there in **both** directions, so no causal claim is supportable inside the
MHC on the basis of `coloc.abf` whether the posterior favours H3 or H4. MHC genes are
flagged in an `assumption_valid` column and are reclassified accordingly, rather than
being reported under the same verdict vocabulary as genes outside the region. Outside
the MHC the assumption is defensible and PP.H3 carries its usual meaning.

### 2.7.3 `coloc.susie` — the multiple-causal-variant model

To address the MHC on its own terms, colocalisation was repeated under a model that
permits multiple causal signals per trait. SuSiE, the "Sum of Single Effects" model
(Wang et al., 2020), decomposes each trait's regional association into a set of
credible sets, each representing one independent causal signal. `coloc.susie`
(Wallace, 2021) then tests colocalisation between every pair of credible sets — eQTL
signal *i* against genome-wide association signal *j* — returning a posterior for
each pair. A gene may therefore be shown to colocalise with **one** of several RA
signals while remaining independent of the others, which is precisely the situation
`coloc.abf` cannot represent and precisely the situation the MHC presents.

SuSiE was fitted with a maximum of ten credible sets per trait (*L* = 10) at 95%
coverage, the package defaults (Wang et al., 2020). *L* = 10 comfortably exceeds the
number of independent RA association signals mapped within the MHC (Raychaudhuri et
al., 2012), so the limit is not binding on the result. Regional statistics were read
from the cache written in §2.7.2, so no re-download was required and the two
colocalisation analyses are guaranteed to operate on identical data.

**What this analysis requires that `coloc.abf` did not, and the limitations that
follow.** `coloc.susie` needs a regional linkage-disequilibrium matrix in the same
allele orientation as the summary statistics. That matrix was obtained from the 1000
Genomes European reference panel (1000 Genomes Project Consortium, 2015) and
re-oriented to the effect alleles of the summary statistics. Two consequences are
stated rather than buried.

1. **Reference-panel mismatch.** The linkage-disequilibrium matrix is out-of-sample,
   derived from approximately 500 reference individuals rather than from the eQTL or
   genome-wide association samples themselves. In the MHC, where haplotype structure
   is extreme, an out-of-sample reference is a genuine source of error and can both
   create and destroy credible sets (Zou et al., 2022). In-sample linkage
   disequilibrium would be required to settle the region, and neither source dataset
   releases it.
2. **Variant cap.** The server-side matrix service is capped, so each region was
   reduced to its 450 most strongly associated shared variants, ranked by the smaller
   of the two traits' *p*-values. Truncating a region can drop a causal variant.

This analysis therefore **improves upon `coloc.abf` within the MHC without settling
it**, and its output is reported as *indicative* rather than definitive. It supersedes
the `coloc.abf` verdict for MHC genes only; outside the MHC the single-causal-variant
assumption is tenable and §2.7.2 stands. Where SuSiE identifies no credible set in
one or both traits, the gene is reported as **not resolvable** — a statement about
statistical power and reference-panel quality, and **never** as evidence of absence.

### 2.7.4 The reading rule adopted throughout this thesis

The following rule was fixed before the colocalisation results were examined and is
applied without exception.

A gene may carry a **causal** claim only if it is **both** (i) classified ROBUST in
the MHC sensitivity analysis of §2.6, meaning it survives false-discovery-rate
correction with MHC instruments excluded, **and** (ii) supported at PP.H4 ≥ 0.80 at
**both** priors in §2.7.2. A gene meeting only the first condition is described as
*associated with*, and never as *causal for*, RA.

Three corollaries follow and are stated so that they cannot be quietly set aside.
PP.H3 ≥ 0.80 outside the MHC is positive evidence that the MR estimate is confounded
by linkage disequilibrium, not merely an absence of evidence for colocalisation, and
is reported as such. Genes classified as inconclusive are underpowered for
colocalisation — the discovery outcome dataset predates contemporary imputation
panels, so regional variant density rather than the method is the limiting factor —
and are reported as *untested*, never as negative. And where no gene in a panel meets
the rule, the Mendelian randomisation stage retains its standing as a
**genetically informed filter on the candidate space**, which is a legitimate and
unusual way to construct a candidate set, but it forfeits its standing as evidence of
causality, and the vocabulary used throughout the thesis changes accordingly.

---

## References cited in this section

1000 Genomes Project Consortium (2015) 'A global reference for human genetic variation', *Nature*, 526(7571), pp. 68–74.

de Bakker, P.I.W., McVean, G., Sabeti, P.C., et al. (2006) 'A high-resolution HLA and SNP haplotype map for disease association studies in the extended human MHC', *Nature Genetics*, 38(10), pp. 1166–1172.

Giambartolomei, C., Vukcevic, D., Schadt, E.E., et al. (2014) 'Bayesian test for colocalisation between pairs of genetic association studies using summary statistics', *PLoS Genetics*, 10(5), e1004383.

Raychaudhuri, S., Sandor, C., Stahl, E.A., et al. (2012) 'Five amino acids in three HLA proteins explain most of the association between MHC and seropositive rheumatoid arthritis', *Nature Genetics*, 44(3), pp. 291–296.

Wallace, C. (2020) 'Eliciting priors and relaxing the single causal variant assumption in colocalisation analyses', *PLoS Genetics*, 16(4), e1008720.

Wallace, C. (2021) 'A more accurate method for colocalisation analysis allowing for multiple causal variants', *PLoS Genetics*, 17(9), e1009440.

Wang, G., Sarkar, A., Carbonetto, P. and Stephens, M. (2020) 'A simple new approach to variable selection in regression, with application to genetic fine mapping', *Journal of the Royal Statistical Society: Series B*, 82(5), pp. 1273–1300.

Zhu, Z., Zhang, F., Hu, H., et al. (2016) 'Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets', *Nature Genetics*, 48(5), pp. 481–487.

Zou, Y., Carbonetto, P., Wang, G. and Stephens, M. (2022) 'Fine-mapping from summary data with the "Sum of Single Effects" model', *PLoS Genetics*, 18(7), e1010299.
