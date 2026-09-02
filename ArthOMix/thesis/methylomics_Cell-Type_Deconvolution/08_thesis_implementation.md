# Cell-Type Deconvolution — Thesis Implementation

Seven sequential tabs comprise the Cell-Type Deconvolution module: Data & QC, CpG Feature
Selection, Reference & Method, Deconvolution, Cell Composition, Validation, and Export. Using a
bulk DNA methylation beta-value matrix (CpG probes by samples), they estimate the underlying
cell-type proportions in each sample. An optional SNP-associated probe filter, an optional
cross-reactive probe filter, an optional detection p-value filter (only when raw IDAT data is
available), a missing-CpG threshold, chromosome-scope filtering (removal of sex chromosomes), and
finally a missing-sample threshold are applied to the working matrix before estimation. Probe and
sample counts are reported through this process in a retention cascade table on the Data & QC tab;
however, the per-step statistics following the first filter should be interpreted as approximate
rather than precise.

The EpiDISH package's (Teschendorff et al., 2017) Houseman-style constrained projection (Houseman
et al., 2012), robust partial correlations, and CIBERSORT-style support-vector methods (Newman et
al., 2015), along with a two-stage hierarchical option (hepidish; Zheng et al., 2018), are used to
execute deconvolution itself. These methods are run against either a user-supplied reference or
one of seven published reference panels, including the Reinius/Houseman blood panel (Reinius et
al., 2012) and the Salas 12-cell-type blood panel (Salas et al., 2022). By reconstructing observed
methylation from the fitted mixture, the resulting cell-type fractions can be internally confirmed,
statistically compared among phenotype groups, and visualised.

The biological validity of any downstream differential methylation or regional analysis is
strengthened by quantifying and accounting for cell-type heterogeneity, a well-researched
confounder in bulk methylation studies (Houseman et al., 2012). Differences attributed to disease
or treatment may instead reflect changes in underlying cell composition.

## References

Houseman EA, Accomando WP, Koestler DC, Christensen BC, Marsit CJ, Nelson HH, Wiencke JK, Kelsey
KT. DNA methylation arrays as surrogate measures of cell mixture distribution. *BMC
Bioinformatics*. 2012;13:86. doi:10.1186/1471-2105-13-86

Newman AM, Liu CL, Green MR, Gentles AJ, Feng W, Xu Y, Hoang CD, Diehn M, Alizadeh AA. Robust
enumeration of cell subsets from tissue expression profiles. *Nat Methods*. 2015;12(5):453-457.
doi:10.1038/nmeth.3337

Reinius LE, Acevedo N, Joerink M, Pershagen G, Dahlén SE, Greco D, Söderhäll C, Scheynius A, Kere
J. Differential DNA methylation in purified human blood cells: implications for cell lineage and
studies on disease susceptibility. *PLoS One*. 2012;7(7):e41361. doi:10.1371/journal.pone.0041361

Salas LA, Zhang Z, Koestler DC, Butler RA, Hansen HM, Molinaro AM, Wiencke JK, Kelsey KT,
Christensen BC. Enhanced cell deconvolution of peripheral blood using DNA methylation for
high-resolution immune profiling. *Nat Commun*. 2022;13:761. doi:10.1038/s41467-021-27864-7

Teschendorff AE, Breeze CE, Zheng SC, Beck S. A comparison of reference-based algorithms for
correcting cell-type heterogeneity in Epigenome-Wide Association Studies. *BMC Bioinformatics*.
2017;18:105. doi:10.1186/s12859-017-1511-5

Zheng SC, Breeze CE, Beck S, Teschendorff AE. A novel cell-type deconvolution algorithm reveals
substantial contamination by immune cells in saliva, buccal and cervix. *Epigenomics*.
2018;10(7):925-940. doi:10.2217/epi-2018-0037
