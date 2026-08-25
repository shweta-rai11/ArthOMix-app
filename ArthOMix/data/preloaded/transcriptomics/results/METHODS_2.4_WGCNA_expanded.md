# 2.4 Co-expression network analysis — WGCNA

*Complete replacement for the current Section 2.4, with in-text citations. Every
parameter stated here matches `scripts/00_shared/06_WGCNA.R` as executed
(`results/logs/06_WGCNA_run3.log`). Full reference list at the end of the document.*

---

## 2.4.1 Rationale: what weighted gene co-expression network analysis does

The differential-expression analysis of Section 2.3 tests each transcript in isolation,
treating the 15,763 genes as 15,763 independent hypotheses. This is a deliberate
simplification, and a biologically implausible one: genes act through co-regulated
programmes, and a coordinated shift distributed across hundreds of functionally related
transcripts may leave every individual gene short of significance while the programme
itself is strongly disease-associated. A gene-by-gene analysis is also blind to
*context* — it cannot distinguish a transcript that is dysregulated as part of a
coherent inflammatory programme from one that moves in isolation, although the two carry
very different evidential weight as candidate biomarkers.

Weighted gene co-expression network analysis (WGCNA) addresses this by reversing the
unit of analysis (Zhang & Horvath, 2005; Langfelder & Horvath, 2008). Rather than
testing genes, it constructs a network in which every gene is a node and the weight of
the edge between two genes is a function of their expression correlation across samples,
then partitions that network into **modules** — clusters of densely interconnected genes
whose expression varies together across individuals. Each module is summarised by a
single representative variable, its **module eigengene**, defined as the first principal
component of the expression of its member genes (Langfelder & Horvath, 2007).
Association with disease is then tested at the level of the module eigengene rather than
the gene, which reduces the multiple testing burden from thousands of tests to a few
dozen and asks a question about coordinated programmes rather than isolated transcripts.

Three properties make it the appropriate method for this thesis specifically.

First, **it is unsupervised**. Modules are defined from the correlation structure of the
expression data alone, with no reference to diagnosis. Disease status enters only
afterwards, when the resulting eigengenes are correlated with the RA indicator. A module
that emerges from the data without knowledge of the phenotype and *then* proves strongly
associated with it constitutes stronger evidence than a gene set assembled by selecting
on that same phenotype.

Second, **it supplies a second, independent line of evidence for candidate selection**.
The candidate genes carried forward to Mendelian randomisation are defined as
`disease-module genes ∩ sex-stratified DEGs` (Section 2.5). Module membership is a
network-level, sample-correlation-based property; differential expression is a
gene-level, group-mean-based property. Requiring both means a candidate must be
individually dysregulated *and* embedded in a coordinated disease programme, and the two
criteria fail in different ways, so their intersection is more robust than either alone.

Third, **it permits a formal, sample-size-controlled test of whether network structure
differs between the sexes**. The male stratum is far smaller than the female stratum, so
comparing DEG lists between sexes conflates a genuine biological difference with
differential statistical power. Module preservation statistics (Section 2.4.14) are
constructed by permutation specifically to control for the size of the test set
(Langfelder et al., 2011), and so address the sex question in a way that a comparison of
gene lists cannot.

The overall workflow — network construction, module detection, module preservation
between two sample groups, and functional enrichment of the resulting modules —
corresponds to the published protocol of Nguyen & Zeng (2025), against which the
parameter choices below are benchmarked and, where they differ, explicitly justified.

**Limitation stated at the outset.** Co-expression is a correlational, undirected
construct. Modules identify genes that vary together; they do not establish regulatory
direction, causal ordering, or physical interaction. All causal inference in this thesis
is carried by the Mendelian randomisation of Chapter *n*, for which the network analysis
supplies candidates rather than conclusions.

---

## 2.4.2 Input data, software and reproducibility

The network was constructed on the **70 % training split only** (183 samples × 15,763
genes common to both platforms), using the quantile-normalised (Bolstad et al., 2003),
ComBat-corrected (Johnson et al., 2007) matrix of Section 2.1. The held-out 30 % was
never loaded by the network script, so no information from the evaluation set can enter
module definition, disease-module selection or candidate identification.

Analyses used R 4.4.2 with WGCNA 1.74 (Langfelder & Horvath, 2008) and clusterProfiler
4.14.6 (Yu et al., 2012; Wu et al., 2021). A fixed seed (`set.seed(1234)`, and
`randomSeed = 1234` passed to `blockwiseModules`) was applied throughout. All tunable
parameters are declared in a single configuration block, and the analysis cache is keyed
to a hash of every network-relevant parameter together with a fingerprint of the input
matrix, so that any change to a parameter or to the data forces recomputation and a
stale network cannot be silently reused.

---

## 2.4.3 Data-integrity screening

Genes and samples were first screened with `goodSamplesGenes` (Langfelder & Horvath,
2008), which removes genes of zero variance or excessive missingness and samples with
too many missing entries. This is a precondition rather than a filter: a zero-variance
gene has an undefined correlation and would propagate missing values through the entire
adjacency matrix. **All 15,763 genes and all 183 samples passed**, so the step removed
nothing; it is reported as evidence that the preprocessing of Section 2.1 left no
degenerate features.

---

## 2.4.4 Sample-outlier detection and removal

Following the procedure of the WGCNA authors' tutorial (Langfelder & Horvath, 2008),
samples were clustered by average-linkage hierarchical clustering on Euclidean distance,
and a data-driven outlier cut height was derived as the mean plus three standard
deviations of the merge heights. Arrays falling outside the principal cluster at that
height were **removed**, not merely reported.

**Why removal rather than inspection.** A weighted co-expression network is constructed
entirely from between-sample correlations (Zhang & Horvath, 2005), so a single aberrant
array perturbs every one of the ~124 million gene pairs simultaneously and generates
co-expression structure that reflects the artefact rather than the biology.
Co-expression analysis is considerably more sensitive to this failure than differential
expression, in which a poor array inflates the residual variance of one gene at a time
and is further absorbed by the empirical array quality weights of Section 2.3 (Ritchie
et al., 2006). Detecting outliers and then retaining them is therefore not the
conservative option but an uncorrected source of spurious modules.

**Why average linkage on Euclidean distance.** Average linkage avoids the chaining
behaviour of single linkage, which absorbs isolated arrays into the principal cluster,
while not imposing the equal-sized, spherical clusters assumed by Ward's method (Ward,
1963). Neither property is appropriate when the target is a small number of isolated
arrays standing against one large group.

**Why a computed cut height rather than a visual one.** The tutorial selects the cut
height by eye from the dendrogram, which is unreproducible and leaves the criterion
unstated. Deriving the cut from the merge-height distribution makes it a computed
property of the data that transfers unchanged to a new dataset. Three standard
deviations is the conventional choice and is the conservative direction, removing fewer
arrays than a two-standard-deviation rule.

**Alternatives considered and not adopted.** Standardised connectivity (Z.k < −2.5;
Oldham et al., 2012) measures the same construct through a different metric and was not
run; principal-component or MDS inspection supplies no threshold and uses only the
leading components; and raw-intensity diagnostics such as `arrayQualityMetrics`
(Kauffmann et al., 2009) or the RLE and NUSE statistics (Brettschneider et al., 2008)
operate on probe-level intensities and no longer apply to a matrix already
quantile-normalised and ComBat-corrected.

**Independent corroboration.** The removed set was cross-checked against a criterion
computed from an entirely different quantity — the limma REML array quality weights of
Section 2.3 (Ritchie et al., 2006), which derive from gene-wise residual variance rather
than between-sample distance. Six of the ten removed arrays are also among the ten
lowest-weighted arrays under that criterion (hypergeometric *p* = 8.5 × 10⁻⁷). The
convergence of two unrelated quality criteria on the same arrays establishes that the
removal identifies genuinely poor-quality data rather than an arbitrary trim.

The network was accordingly built on **173 samples**.

---

## 2.4.5 Gene filtering: none applied

All 15,763 genes surviving the integrity screen entered the network; **a variance filter
was not applied.**

**This is a deliberate departure from common practice, including from the protocol
followed elsewhere in this section.** Nguyen & Zeng (2025) describe filtering out
low-variance genes as "a key preprocessing step in WGCNA" and retain only the top 5 % of
genes by variance, and the WGCNA package FAQ permits filtering by mean or variance
(Langfelder & Horvath, 2008). The departure is justified by a structural feature of the
present design that does not apply to those settings.

**The structural argument.** Where the module is itself the endpoint of the analysis —
as in Nguyen & Zeng (2025), in which modules are characterised and compared between
tumour and normal tissue — a variance filter discards genes that would have contributed
little to module structure and costs the result little. Here the network is not the
endpoint. Its sole function is to supply the disease-module background that is
intersected with each sex's DEG list to define the candidate genes carried into
Mendelian randomisation (Section 2.5). A gene removed before module detection can never
be assigned to a disease module, can never enter that intersection, and can therefore
never be tested for a causal relationship with RA. The filter would operate as a silent
and permanent veto over the study's final output, exercised on marginal variability — a
criterion bearing no relationship to the question being asked.

**The cost was measured rather than assumed.** An earlier configuration of this pipeline
applied a 40 % variance filter, which removed 6,305 of the 15,763 genes and with them
**14.9 % of the 5,131 female DEGs and 20.9 % of the 5,820 male DEGs**. The loss is
asymmetric, and it falls in the direction least affordable to this study: the male
stratum comprises 38 of the 183 training samples and 17 RA cases, and is already the
weaker of the two in statistical power (Section 3.1), so a filter removing a fifth of
its DEGs erodes precisely the stratum this chapter exists to characterise. A
95th-percentile filter of the kind used in the protocol literature would have been more
severe still. No computational consideration required the filter in compensation: all
15,763 genes across 173 samples were accommodated in a single block, so modules were
detected globally rather than block-wise.

**Why not the other filtering strategies.** A mean-expression filter addresses the
count-driven mean–variance dependence characteristic of RNA-seq data (Law et al., 2014),
which is absent from a quantile-normalised, batch-corrected microarray matrix. Retaining
a fixed number of the most variable genes replaces a percentile with an equally
arbitrary round number while inheriting the same structural objection. Filtering to
differentially expressed genes is explicitly forbidden by the WGCNA FAQ (Langfelder &
Horvath, 2008) and would render the analysis circular, since genes selected for their
association with RA would then be used to form modules whose correlation with RA is
reported as a finding — the form of selection bias described by Ambroise & McLachlan
(2002) and Simon et al. (2003). Filtering for computational tractability was
unnecessary: all 15,763 genes were accommodated in a single block, so modules were
detected globally rather than block-wise.

**Why this does not admit noise.** Low-information genes are not ignored; they are
excluded *softly*, by the network itself. Soft-thresholding at β = 12 drives weak
correlations toward zero (Zhang & Horvath, 2005), so near-invariant genes acquire
negligible connectivity and are assigned to the unassigned (grey) module rather than
distorting genuine ones — 6,019 genes were so assigned. The distinction is that this
exclusion is applied *after* the correlation structure is known and on the basis of
connectivity, whereas a variance filter is applied *before* it and on the basis of
marginal variability alone: a gene of modest variance that genuinely co-expresses with a
disease module survives the former but cannot survive the latter.

---

## 2.4.6 Correlation coefficient

Pearson correlation was used (`corType = "pearson"`).

The WGCNA FAQ recommends the biweight midcorrelation (`bicor`, `maxPOutliers = 0.05`) in
preference to Pearson, because it is resistant to outlying observations *within* a gene
(Langfelder & Horvath, 2008). A small number of extreme expression values can inflate a
Pearson correlation and thereby create edges — and ultimately module membership — that
reflect a handful of arrays rather than a consistent co-expression relationship.

`bicor` was implemented and evaluated, and was reverted on computational grounds.
Whereas Pearson correlation has a closed-form solution, `bicor` estimates its weights
iteratively, requiring additional working copies of the data at roughly three times the
computational cost. This interacts poorly with the decision of Section 2.4.5 to retain
all 15,763 genes: a single 15,763 × 15,763 double-precision matrix occupies
approximately 2 GB, and the topological-overlap calculation holds the adjacency matrix
and the TOM in memory simultaneously together with intermediates. On the analysis
machine (16 GB) the `bicor` run entered swap and did not complete, whereas the Pearson
network completed in a single block.

This is a computational constraint and is not advanced as a claim that Pearson is
preferable to `bicor`. It is also, and should be acknowledged as, a consequence of the
gene-filtering decision: on the available hardware, retaining every gene and using the
FAQ-recommended correlation were mutually exclusive. Of the two FAQ-permitted
departures, the one preserving every gene's eligibility for the candidate intersection
was chosen, because a gene removed by a variance filter is lost permanently from the
study's output (Section 2.5), whereas the robustness forgone by using Pearson is partly
recoverable by other means.

Robustness to outliers was accordingly addressed at the **sample** level instead. The
ten outlying arrays identified in Section 2.4.4 were removed, and that removal was
independently corroborated against the limma array quality weights (Ritchie et al.,
2006). Because an aberrant array contributes an extreme value to every gene
simultaneously, deleting it removes the most common source of within-gene outliers at
source: an array-driven artefact cannot survive the deletion of the array that produced
it.

**Limitation.** Sample-level removal does not address within-gene outliers arising from
causes other than global array quality — for example a probe behaving erratically in a
subset of individuals. Such values remain capable of inflating the Pearson correlations
involving that gene and could in principle strengthen or create its module membership.
Three considerations bound this risk without eliminating it: quantile normalisation
(Bolstad et al., 2003) constrains the range of extreme values; the disease modules are
large, so no single gene's correlations determine the module eigengene; and hub genes
are defined by |kME| and |GS| jointly (Section 2.4.12), so a gene whose apparent
centrality rested on a few extreme observations would still require a genuine
correlation with disease status to be reported. The residual risk is nonetheless real
and is recorded rather than discounted.

---

## 2.4.7 Network type and topological overlap

A **signed** network was constructed (`networkType = "signed"`, `TOMType = "signed"`),
consistent with the protocol of Nguyen & Zeng (2025).

**Why signed rather than unsigned.** An unsigned network uses |cor|, so a gene pair that
is strongly *anti*-correlated is treated as strongly connected and the two genes may be
placed in the same module. For a disease-programme analysis this is undesirable: a
transcript up-regulated in RA and one down-regulated in RA are biologically opposed, and
a module containing both has no coherent direction, no interpretable eigengene, and
cannot be described as "up in RA" or "down in RA" — which is precisely the description
required by the directional consistency check of Section 2.5. A signed network preserves
the sign of the correlation, so modules are directionally homogeneous and each carries an
unambiguous relationship to disease status.

**Why not signed-hybrid.** Signed-hybrid sets negative correlations to zero rather than
mapping them to low adjacency. Fully signed was preferred for its monotone, continuous
adjacency transformation; the two behave similarly for module detection, and the
distinction is not material to any conclusion drawn here.

**Why topological overlap rather than raw correlation.** Modules were defined from the
topological overlap matrix (TOM) rather than the adjacency matrix directly (Zhang &
Horvath, 2005). The TOM measures the extent to which two genes share their network
neighbours in addition to being connected to each other, so it reflects shared
regulatory context rather than pairwise association alone. This substantially suppresses
spurious edges arising from noise in individual correlations, since a spurious pair is
unlikely to share a neighbourhood, and yields more robust and reproducible modules than
clustering on correlation directly.

---

## 2.4.8 Soft-thresholding power

Scale-free topology fit and mean connectivity were evaluated across candidate powers
(1–10, then 12–20 in steps of two) with `pickSoftThreshold` for a signed network.

**The scale-free criterion, and what it does and does not license.** Soft-thresholding
raises the correlation to a power β, so that strong correlations are preserved while weak
ones are driven toward zero — a continuous alternative to imposing a hard correlation
cut-off, which would discard the magnitude information that distinguishes a correlation
of 0.9 from one of 0.5 (Zhang & Horvath, 2005). β is chosen so that the resulting
network approximates scale-free topology, the connectivity distribution characteristic
of biological networks.

`pickSoftThreshold` returned an estimate of **β = 3**, the first power reaching
R² ≥ 0.85 (R² = 0.873). **This power was rejected on connectivity grounds.** The
estimate is defined as the *first* power clearing the R² threshold, not the best one, and
at β = 3 the mean connectivity is 2,240 — a network in which the average gene is
connected to one in seven of all others, far too dense for meaningful module resolution.
The signed R² curve plateaus at approximately 0.90 from β = 9 onwards while mean
connectivity falls steeply, so powers in the plateau achieve equivalent scale-free fit at
a far better conditioned connectivity.

**β = 12 was used**, where R² = 0.901, slope = −1.78 and mean connectivity = 47.4. The
power was therefore selected on the **joint criterion of scale-free fit and network
connectivity**, with the full fit-versus-connectivity table reported at every candidate
power (Section 3.5) so that the choice can be audited.

This joint criterion is the one specified by Nguyen & Zeng (2025), who select β by
requiring scale-free fit R² ≥ 0.9 while maintaining mean connectivity below 300. The
value used here satisfies both of those published thresholds (R² = 0.901; mean
connectivity 47.4), whereas the automatic estimate of β = 3 satisfies neither
(R² = 0.873; mean connectivity 2,240). The rejection of the automatic estimate is
therefore not an ad hoc preference but the application of a criterion set out
independently in the published protocol literature.

**On the WGCNA FAQ sample-size table.** The FAQ recommends β = 12 for a signed network
with more than 40 samples, which coincides with the value used. This coincidence is
reported as *corroboration only* and is deliberately not offered as the justification,
because the FAQ table is formally **conditional on the scale-free fit failing**. Here the
fit did not fail — an estimate was returned — so the table does not formally apply, and
citing it as the primary reason would misapply the authors' own guidance. The
connectivity argument is the justification; the agreement of the FAQ value is a
supporting observation.

---

## 2.4.9 Module detection

Modules were detected in a single block with `blockwiseModules` at β = 12, using the
dynamic tree cut algorithm (Langfelder, Zhang & Horvath, 2008). The parameters, and the
reason for each, are set out in Table 2.x.

**Table 2.x — Module-detection parameters, with justification.**

| Parameter | Value | Why this value, and why not another |
|---|---|---|
| `maxBlockSize` | 20000 | Exceeds the 15,763 genes, so the network is built in **one block**. The default (5,000) would split the genes into four blocks clustered independently, and modules could not span blocks — an arbitrary partition of the gene space imposed by memory management rather than biology. Raised from 16,000 when the variance filter was removed. |
| `minModuleSize` | 30 | The package default, and the value used by Nguyen & Zeng (2025). Smaller values (10–20) fragment the dendrogram into many small modules whose eigengenes are unstable and whose enrichment is uninterpretable; larger values (50+) merge distinct programmes. Left at the default rather than tuned, so that module count was not adjusted toward a preferred outcome. |
| `mergeCutHeight` | 0.25 | Merges modules whose eigengenes correlate above 0.75; again the default and the protocol value (Nguyen & Zeng, 2025). Two modules at that correlation are largely redundant summaries of the same variation, and retaining both would inflate the number of module–trait tests without adding information. A stricter value (0.15, r > 0.85) leaves near-duplicate modules; a looser one (0.35, r > 0.65) merges genuinely distinct programmes. |
| `deepSplit` | 2 (default) | Controls dynamic tree cut sensitivity on a 0–4 scale; the value used by Nguyen & Zeng (2025). Left at the default; not tuned. |
| `reassignThreshold` | 0 | Disables post-hoc reassignment of genes between modules on the basis of kME significance. Genes therefore remain in the module assigned by the dendrogram cut, so module membership is determined by one criterion (topological overlap) rather than two, and assignment cannot drift under a second, correlation-based rule. |
| `pamRespectsDendro` | FALSE | Allows PAM-stage assignment of otherwise unassigned genes on topological overlap without being constrained by dendrogram branch boundaries. This is the setting in the authors' tutorial. |
| `numericLabels` | TRUE | Modules are produced as integers and converted to colours only for display (see below). |
| `randomSeed` | 1234 | Fixes the stochastic component so module assignment is exactly reproducible. |

Detection at these settings yielded **12 modules**, with 6,019 genes unassigned (grey).
The grey module is not a module — it collects genes not co-expressed with any coherent
group — and it is excluded from disease-module selection and from all downstream
analysis.

Each module was summarised by its **module eigengene**, the first principal component of
its member genes' expression, and eigengenes were ordered by hierarchical clustering
(`orderMEs`; Langfelder & Horvath, 2007).

**On module colour names.** `labels2colors` assigns colours by module *size rank*, so
the same biology is renamed whenever the gene filter, the soft power or the sample set
changes; in this project the same disease programme was named "green" under the filtered
configuration and "yellow" under the unfiltered one. Colour names are therefore recorded
as an **output** of the analysis and are never used as an **input**: no downstream step
refers to a colour literal, and every module is selected by its correlation with RA.

---

## 2.4.10 Module–trait association

A binary trait matrix encoded disease status (RA, HC), sex, the four group-by-sex
combinations, and age. Each trait was correlated with every module eigengene by Pearson
correlation, with Student asymptotic *p* values (`corPvalueStudent`; Langfelder &
Horvath, 2008).

The same association was additionally recomputed **within each sex stratum separately**
(all samples, female only, male only), so that the disease association of every module
can be read inside each sex. The female and male columns are separate within-sex
contrasts and are **not compared with one another**: a module correlating more strongly
with RA in men than in women is not evidence of a sex difference, because the male
stratum is far smaller and its correlations correspondingly noisier and more extreme
(Gelman & Stern, 2006). The formal sex comparison is the preservation analysis of
Section 2.4.14.

*Multiple testing.* The reported *p* values are unadjusted Student asymptotic values.
No correction was applied because the disease-module threshold used for selection
(*p* < 1 × 10⁻⁸, Section 2.4.11) is several orders of magnitude more stringent than a
Bonferroni correction over the tests performed, so the selection decision is unaffected
by the choice of adjustment.

---

## 2.4.11 Disease-module selection: a pre-specified, data-driven rule

Disease-associated modules were selected from the overall, all-sample module–trait
analysis by a pre-specified rule: **|cor(ME, RA)| ≥ 0.5 and association p < 1 × 10⁻⁸**,
with the grey module excluded.

**Why a rule rather than a ranking.** Describing modules as "ranked by their correlation
with RA" and then selecting the top ones leaves the number selected undetermined and the
threshold implicit, so the selection can be adjusted after seeing the result. An
explicit, pre-declared threshold fixes the decision before the outcome is known and can
be checked by a reader against the reported correlation table.

**Why these thresholds.** The correlation cut-off of 0.5 requires a module to explain a
substantial share of disease-related expression variance rather than merely reaching
significance; the *p* threshold guards against a spurious module–trait relationship. The
two act on different failure modes — effect size and sampling noise — and both must be
met.

**Why both directions are retained.** Modules meeting the rule in either direction were
kept, so that the up-regulated and down-regulated disease programmes are both captured.
Restricting to positive correlations would discard half of the disease signal and would
bias the candidate set toward genes up-regulated in RA.

**Why the same modules are used for both sexes.** The identical module set is carried
forward for both sexes, so that no sex is assigned a different co-expression background
and all sex resolution enters through the sex-stratified expression contrast alone. This
corrects an earlier design in which one module was assigned to each sex — a procedure
that guarantees a difference between the sexes by construction.

---

## 2.4.12 Gene-level importance and hub genes

Within each selected disease module, gene importance was quantified by three
complementary measures (Zhang & Horvath, 2005; Langfelder & Horvath, 2008):

- **module membership (kME)** — the correlation between a gene's expression and its
  module eigengene, i.e. how representative the gene is of the module's programme;
- **gene significance (GS)** — the correlation between a gene's expression and RA
  status, i.e. its individual disease association;
- **intramodular connectivity** — the row sums of the signed adjacency at β = 12
  computed *within* the module, i.e. how central the gene is to the module's topology.

Hub genes were defined as genes with **|kME| > 0.8 and |GS| > 0.2**. Requiring both
means a hub must be central to its module *and* disease-associated; kME alone would
select genes central to a module irrespective of their relationship to RA, and GS alone
would select disease-associated genes irrespective of network position, which is simply
a weaker restatement of the differential-expression analysis. The kME threshold of 0.8
is the conventional value and identifies genes tracking the eigengene closely; the GS
threshold of 0.2 is deliberately permissive, because GS is a whole-cohort correlation
with a binary trait and its magnitude is bounded well below that of kME.

Intramodular connectivity was computed within the module rather than across the whole
network, so that hub status reflects centrality in the disease programme rather than
global connectivity. These analyses were performed separately for each selected disease
module, so that the up-regulated and down-regulated programmes receive equivalent
treatment.

Hub identification here is numerical rather than visual. Nguyen & Zeng (2025) export
edge and node lists for interactive topological exploration in Cytoscape; no such export
was performed, because hub status in this study is defined by explicit numeric
thresholds that are reported in full and require no visual inspection to reproduce.

---

## 2.4.13 Functional enrichment

Gene Ontology enrichment across all three ontologies (`enrichGO`, `ont = "ALL"`) and
KEGG pathway enrichment (`enrichKEGG`, `organism = "hsa"`, following identifier
conversion with `bitr`) were applied to the pooled disease-module gene set with
clusterProfiler (Yu et al., 2012; Wu et al., 2021), using Benjamini–Hochberg correction
(Benjamini & Hochberg, 1995) and p- and q-value cut-offs of 0.05. This is the enrichment
procedure specified by Nguyen & Zeng (2025).

**The background universe was restricted to all genes present in the network** rather
than to the whole genome or the whole annotation. This is the correct comparator: the
disease-module genes were drawn from the genes that entered the network, so those are
the genes that *could* have been selected. Using a genome-wide universe would test the
disease modules against genes that were never eligible, and would report as enrichment
any bias in the platform's probe coverage or in the two-platform gene intersection.

---

## 2.4.14 Sex-stratified networks and module preservation

Two further networks were constructed independently on the female and the male samples
at the same soft-threshold power. Whether co-expression structure differs between sexes
was then tested with `modulePreservation` (Langfelder et al., 2011), using the female
data as the reference network and the male data as the test network, with **200
permutations**.

**Why a permutation statistic rather than a direct comparison.** Preservation Z-summary
statistics are constructed by permuting sample labels within the test set, so the null
distribution is generated at the test set's own size. This controls for the fact that
the male stratum is much smaller than the female stratum. A direct comparison of male
and female module structure, or of male and female DEG lists, does not control for this,
and would report reduced power in the male stratum as a biological sex difference.

**Why the reference assignment comes from the combined network.** The reference module
assignment was taken from the **combined-network** modules rather than from the
female-network modules. This is essential: female-network modules are a different gene
partition from the disease modules even where colour names coincide, so using them would
answer the question "are female-network modules preserved in men?" rather than the
question actually required, "are the *disease* modules preserved between sexes?".

**Interpretation thresholds.** A preservation Z-summary below 2 indicates non-preserved
structure, between 2 and 10 moderate preservation, and above 10 strong preservation
(Langfelder et al., 2011; Nguyen & Zeng, 2025).

**Number of permutations.** Two hundred permutations were run, double the 100 used by
Nguyen & Zeng (2025), to stabilise the permutation null in the presence of a small test
set; the male stratum is considerably smaller than either group in that protocol, so the
Z-summary is correspondingly more sensitive to permutation noise.

**Direction of the comparison.** Preservation was assessed in one direction only, with
the female network as reference and the male network as test. Nguyen & Zeng (2025)
additionally reverse the reference and test assignments to assess preservation
bidirectionally. The single direction was chosen here because the male stratum is too
small to serve as a stable reference network, and the question of interest is whether
structure established in the larger stratum recurs in the smaller one. That the
comparison is therefore not symmetric is recorded as a limitation.

**What this statistic does not license.** It is a network *stability* statistic. It is
not a per-gene test of whether the RA effect differs between the sexes, and no
diagnosis-by-sex interaction test is performed anywhere in this chapter. No gene is
described as sex-specific on the basis of appearing in one stratum's list and not the
other's.

---

## 2.4.15 Relationship to the differential-expression analysis

The differential-expression contrasts of Section 2.3 were fitted on all 183 training
samples, whereas the network was constructed on the 173 samples surviving array-outlier
removal. Because the candidate set of Section 2.5 intersects a gene list derived from
the first analysis with one derived from the second, the discrepancy requires
justification rather than mere disclosure.

**The two analyses did not ignore array quality differently; they controlled it by the
mechanism each method provides.** limma incorporates array quality *inside* the
estimator, through the REML empirical array weights of Section 2.3 (Ritchie et al.,
2006), which down-weight poor arrays continuously in proportion to their residual
variance rather than deleting them. WGCNA offers no equivalent: `blockwiseModules`
admits no per-sample weight, and every sample contributes equally to every pairwise
correlation, so deletion is the only available lever. Applying weighting where weighting
is available and deletion where it is not is a method-appropriate choice, not an
inconsistency.

**The difference between the two treatments is of degree, not of kind.** All ten arrays
removed from the network carry a limma array weight below 1.0, and six of them fall
below 0.5 — six of only eight such arrays in the entire training set. The median weight
of the removed set is 0.482 against 1.079 for the retained samples, a 2.2-fold
difference. These arrays were therefore already heavily discounted in the
differential-expression analysis; they were not contributing to the DEG results as
full-weight observations. The network deletes what limma down-weights.

**Deletion was not adopted for the differential-expression analysis because its cost
runs in the opposite direction.** Statistical power in a linear model depends directly
on the number of samples per group, and the male stratum comprises only 38 samples and
17 RA cases (Section 3.1); removing ten arrays would erode the very contrast this
chapter is least able to afford to weaken. Weighting retains the information in a
marginal array while discounting its unreliable component. A co-expression network, by
contrast, loses little between 173 and 183 samples — correlation estimates are stable at
either size, and the WGCNA FAQ's minimum is 15 to 20 samples — while being acutely
sensitive to a single aberrant array, which perturbs every gene pair simultaneously. The
trade-off between sample count and array quality therefore resolves differently for the
two methods, and it was resolved separately for each.

**The removal does not alter cohort composition materially.** The proportion of RA cases
changes from 56.3 % (103/183) to 54.9 % (95/173) and the proportion of male participants
from 20.8 % to 20.2 %. The removed set contains eight RA and two HC samples, an
imbalance that does not depart significantly from the cohort ratio (hypergeometric
*p* = 0.108). The network is therefore built on the same population as the
differential-expression analysis, not on a differently constituted subset.

**What the intersection requires.** The operation of Section 2.5 is a set intersection
over gene identities within a shared feature space of 15,763 common genes. Its validity
requires that each list be individually valid on the samples from which it was derived,
not that the two sample sets coincide. What the discrepancy affects is the evidential
basis of each list, and the number of samples entering each analysis is therefore
reported explicitly wherever results are given.

---

## Appendix to 2.3 — Consolidated parameter table

Every network-relevant parameter, its value, and the rejected alternative. Useful as a
viva preparation sheet. The final column records agreement or divergence with the
published protocol of Nguyen & Zeng (2025).

| Decision | Chosen | Rejected alternative | One-line reason | vs. protocol |
|---|---|---|---|---|
| Input set | 70 % training split | Full cohort | Prevents leakage into the held-out evaluation set | n/a |
| Integrity screen | `goodSamplesGenes` | — | Removed nothing; reported to show clean preprocessing | — |
| Sample outliers | Removed (10 arrays, 183 → 173) | Examined only | Correlation-based method; one bad array perturbs all gene pairs | — |
| Outlier rule | mean + 3 SD of merge heights | Visual cut height | Computed, reproducible, transfers to new data | — |
| Linkage | Average, Euclidean | Single / Ward | No chaining; no spherical-cluster assumption | — |
| Gene filter | **None** (15,763 genes) | 40 % variance filter | Filtered genes can never reach `module ∩ DEG`; cost 14.9 % F / 20.9 % M DEGs | **Diverges** (protocol keeps top 5 %) |
| — | | Mean-expression filter | RNA-seq remedy; irrelevant to a normalised microarray matrix | — |
| — | | Filter to DEGs | Forbidden by FAQ; circular | — |
| Correlation | Pearson | `bicor` | FAQ-preferred but exceeded memory on the unfiltered network; sample-level removal substitutes | — |
| Network type | Signed | Unsigned | Unsigned merges up- and down-regulated genes; no interpretable module direction | Agrees |
| Similarity | TOM | Raw adjacency | Shared neighbourhood suppresses spurious single-pair edges | Agrees |
| Soft power | β = 12 (R² = 0.901, k̄ = 47.4) | β = 3 (auto estimate) | Auto estimate is the *first* power clearing R², k̄ = 2,240 — fails both protocol thresholds | Agrees (R² ≥ 0.9, k̄ < 300) |
| — | | FAQ table value | Coincides at 12, but the table applies only when the fit fails; corroboration only | — |
| Block size | 20000 (one block) | 5000 (default) | Default splits genes into 4 blocks; modules could not span them | — |
| `minModuleSize` | 30 | 10–20 / 50+ | Default; not tuned | Agrees |
| `mergeCutHeight` | 0.25 (r > 0.75) | 0.15 / 0.35 | Default; merges redundant eigengenes without merging distinct programmes | Agrees |
| `deepSplit` | 2 | 0–4 | Default; not tuned | Agrees |
| `reassignThreshold` | 0 | Non-zero | Membership determined by topological overlap alone, not a second kME rule | — |
| Module identity | correlation with RA | colour name | `labels2colors` assigns by size rank; colours are unstable across runs | — |
| Disease-module rule | \|cor\| ≥ 0.5 **and** p < 1×10⁻⁸ | Top-*n* ranking | Threshold pre-specified; ranking leaves the cut-off implicit and post-hoc | — |
| Direction | Both retained | Positive only | Discarding down-regulated modules halves the disease signal | — |
| Hub rule | \|kME\| > 0.8 **and** \|GS\| > 0.2 | Either alone | kME alone ignores disease; GS alone restates the DEG analysis | — |
| Hub exploration | Numeric thresholds | Cytoscape edge/node export | Thresholds fully reported; no visual step needed to reproduce | Diverges (protocol exports) |
| Enrichment | `enrichGO` + `enrichKEGG`, BH, 0.05 | — | Standard | Agrees |
| Enrichment universe | Genes in the network | Whole genome | Only network genes could have been selected; avoids platform-coverage bias | — |
| Preservation | 200 permutations | 100 permutations | Small male stratum; stabilises the permutation null | Diverges (protocol uses 100) |
| Preservation direction | F → M only | Bidirectional swap | Male stratum too small to serve as a stable reference | Diverges (protocol reverses) |
| Preservation reference | Combined-network modules | Female-network modules | Otherwise tests a different gene partition than the disease modules | — |
| Z-summary bands | < 2 / 2–10 / > 10 | — | Standard interpretation | Agrees |

---

## References cited in Section 2.3

Ambroise C, McLachlan GJ. Selection bias in gene extraction on the basis of microarray
gene-expression data. *Proc Natl Acad Sci USA* 2002;99(10):6562–6.

Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful
approach to multiple testing. *J R Stat Soc Series B* 1995;57(1):289–300.

Bolstad BM, Irizarry RA, Åstrand M, Speed TP. A comparison of normalization methods for
high density oligonucleotide array data based on variance and bias. *Bioinformatics*
2003;19(2):185–93.

Brettschneider J, Collin F, Bolstad BM, Speed TP. Quality assessment for short
oligonucleotide microarray data. *Technometrics* 2008;50(3):241–64.

Gelman A, Stern H. The difference between "significant" and "not significant" is not
itself statistically significant. *Am Stat* 2006;60(4):328–31.

Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data
using empirical Bayes methods. *Biostatistics* 2007;8(1):118–27.

Kauffmann A, Gentleman R, Huber W. arrayQualityMetrics — a bioconductor package for
quality assessment of microarray data. *Bioinformatics* 2009;25(3):415–6.

Langfelder P, Horvath S. Eigengene networks for studying the relationships between
co-expression modules. *BMC Syst Biol* 2007;1:54.

Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis.
*BMC Bioinformatics* 2008;9:559.

Langfelder P, Zhang B, Horvath S. Defining clusters from a hierarchical cluster tree:
the Dynamic Tree Cut package for R. *Bioinformatics* 2008;24(5):719–20.

Langfelder P, Luo R, Oldham MC, Horvath S. Is my network module preserved and
reproducible? *PLoS Comput Biol* 2011;7(1):e1001057.

Law CW, Chen Y, Shi W, Smyth GK. voom: precision weights unlock linear model analysis
tools for RNA-seq read counts. *Genome Biol* 2014;15:R29.

**Nguyen P, Zeng E. A protocol for weighted gene co-expression network analysis with
module preservation and functional enrichment analysis for tumor and normal
transcriptomic data. *Bio Protoc* 2025;15(18):e5447. doi:10.21769/BioProtoc.5447.
PMID: 41000162; PMCID: PMC12457846.**

Oldham MC, Langfelder P, Horvath S. Network methods for describing sample relationships
in genomic datasets: application to Huntington's disease. *BMC Syst Biol* 2012;6:63.

Ritchie ME, Diyagama D, Neilson J, et al. Empirical array quality weights in the
analysis of microarray data. *BMC Bioinformatics* 2006;7:261.

Simon R, Radmacher MD, Dobbin K, McShane LM. Pitfalls in the use of DNA microarray data
for diagnostic and prognostic classification. *J Natl Cancer Inst* 2003;95(1):14–8.

Ward JH Jr. Hierarchical grouping to optimize an objective function. *J Am Stat Assoc*
1963;58(301):236–44.

Wu T, Hu E, Xu S, et al. clusterProfiler 4.0: a universal enrichment tool for
interpreting omics data. *Innovation (Camb)* 2021;2(3):100141.

Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological
themes among gene clusters. *OMICS* 2012;16(5):284–7.

Zhang B, Horvath S. A general framework for weighted gene co-expression network
analysis. *Stat Appl Genet Mol Biol* 2005;4:Article 17.
