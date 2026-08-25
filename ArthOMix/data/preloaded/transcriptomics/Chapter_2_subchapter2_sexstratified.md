Sub-chapter 2 - Sex-Stratified transcriptomic biomarker for RA Diagnosis.

Chapter aim and objectives

Research question:

Does biological sex define distinct whole-blood transcriptomic diagnostic signatures in rheumatoid arthritis, and do these sex-stratified signatures replicate in an independent blood cohort and in synovial tissue? The question is answered in two parts, which are kept separate throughout: whether panels built separately within each sex transfer to independent data (a sex-stratified question), and whether the RA effect on a gene genuinely differs between the sexes (a sex-differential question, answerable only by the interaction test of §2.10)

Sub-question:

- Do sex-stratified signatures discriminate RA from healthy controls within each sex? Female and male performance are reported separately and are not compared as though the comparison were a test of a sex difference, given the male sample-size limit
- Do the blood biomarkers generalise to an independent PBMC cohort (GSE15573) and to synovium (GSE89408), and which genes are concordant across both
- Which pathways and immune-cell composition differences underlie the signatures, and do sex-chromosome genes contribute to classification rather than acting only as sex markers

Objectives:

- To discover RA diagnostic transcriptomic signatures separately in females and males from whole blood
- To test whether the resulting sex-distinct signatures replicate in independent blood and in synovial tissue.
- To identify the tissue-specific biomarker and cross-tissue biomarker in RA female and RA male samples.
- To validate the gene panel against the baseline models.
- To identify the diagnostic performance of the gene panels in males and females.
- To identify cross-ancestral biomarkers in RA
- To validate the biomarkers through Mendelian Randomisation

**General methods and analytical framework - write methodology and dataset used in the study**

**2\. Methodology**

**Terminology used throughout this chapter, and why it is fixed in advance**

Three terms are used precisely and are not interchangeable. They are defined here
because the distinction between them is the difference between a design and a
finding, and because the multi-omics integration of a later chapter requires that
each omics layer state which of the three it is reporting.

**Sex-stratified** describes an *analysis design*: the analysis is run separately
within each sex, with no sample, fold, tuning parameter or fitted coefficient
shared between strata. It is a statement about how the computation was organised
and carries **no claim whatever** that anything differs between the sexes. The
differential-expression contrasts of §2.3, the candidate sets of §2.5, the feature
selection of §2.8 and the diagnostic panels of §2.9 are all sex-stratified in this
sense.

**Sex-differential** describes a *finding*: the effect of RA versus healthy control
on a gene differs between the sexes, established by a formal diagnosis-by-sex
interaction test whose null hypothesis is that the two within-sex effects are
equal (§2.10). This is the only term in this chapter that asserts a difference
between the sexes, and it is applied only to genes that pass that test.

**Sex-specific** — meaning present in one sex and absent in the other — is
**avoided**, and no gene in this chapter is described that way. A significant
interaction establishes that an effect *differs* in magnitude or direction; it does
not establish that the effect is *absent* in one sex, which at these sample sizes
would be an inference from a non-significant result in the smaller stratum. Where
the phrase appears below it appears only in the negative, to state what is *not*
being claimed.

The practical consequence, stated once here and not repeated: a gene appearing in
one stratum's list and not the other's is **sex-stratified output**, not a
sex-differential finding. Only §2.10 produces sex-differential findings.

**2.1 Dataset**

The datasets used for this study are GSE93272 and GSE110169, which are blood samples and are considered training data, while GSE15573 (blood) and GSE89408 (synovium) are the external validation data.

**2.2 Data preprocessing and batch correction**

The datasets were loaded using the R package (GEOquery), which retrieved the biological data such as metadata, annotation data, and expression data. The datasets used were from whole blood having both healthy controls and Rheumatoid Arthritis (RA). The dataset had both male and female samples.

Different algorithms were applied to standardise the two series. In GSE93272, the frozen robust multiarray analysis (fRMA) was used to prepare the transcriptomics data (Tasaki et al., 2018), and for this dataset, GSE110169, the Robust Multi-array Average algorithm (RMA) was used to normalise the data using the "affy" package in R version 3.2.1. (Hu et al., 2018). In this study, the prepared transcriptomics data were used for further analysis.

In the next step, the training datasets were analysed to find the overlapping genes between the two datasets. The overlapping genes were then taken into consideration for further analysis.

Microarray data measure probes, but for our study we are interested in the genes, so one of the important steps is collapsing probes to genes using the WGCNA package (R version 4.4.2) and the method "MaxMean". This method was chosen to retain the single-probe information using the highest mean expression across samples (Langfelder et al., 2008). In this study, the "MaxMean" method is used over averaging because averaging these genes would lead to blending the high-signal probes; moreover, there are multiple probes for the same genes, and these same genes are not replicates; instead, it targets different regions of the transcript and have different GC content, whereas "MaxMean" picks the probe with the highest average intensity which identifies the probes that has the best dynamic range. This step is important as the two different datasets were collected from different platforms and machines.

After collapsing the gene to probe, merging of the data was done using "cbind" as the platform is Microarray for both the datasets, followed by splitting the merged data into training and internal testing data in the ratio of 70:30. The entire analysis further takes up the training data and the internal testing data are kept untouched and only used for internal validation in the later analysis.

**Normalisation**

The "normalizeBetweenArrays" function is used to align the expression values between the samples using methods such as "quantile". (Bolstad et al., 2003). This method aligns the per-sample intensity distributions to a common reference distribution that learned from the training samples.

**Batch correction**

The ComBat empirical-Bayes technique (sva package; Leek et al., 2012) was then used to eliminate residual technical variance owing to study or platform origin from the training set. ComBat was applied with parametric priors (par.prior = TRUE) on the quantile-normalised training matrix, using a batch variable defined at the study-and-internal-batch level. Crucially, the biological signal of interest was protected by providing a model matrix of covariates, that is (mod = ~ group + sex), such that sex-associated expression differences and RA-versus-HC were kept instead of being eliminated as batch effects (Nygaard et al., 2016). The 70% training set was used to estimate all ComBat and normalisation parameters.

**Hold-out internal**

At evaluation, the hold-out was placed on the training scale during validation without re-estimating any parameters. The training ComBat standardisation and empirical-Bayes batch parameters were reapplied out-of-sample after quantile normalisation and were applied to the hold-out against the fixed training distribution (normalise.quantiles.use.target, preprocessCore). By recreating the training matrix from the fitted parameters and comparing it with the saved ComBat output, the fidelity of this frozen approach was verified. This ensured that the internal-validation cohort was on the same scale as the training data.

**2.3 Differential gene expression analysis (DGE)**

Differential gene expression analysis is performed in this study to identify which genes are expressed at different levels between two different conditions, such as healthy and control in female samples and male samples (McDermaid et al., 2019). About 70% of the training set (183 samples × 15,763 common genes) was subjected to differential expression analysis (Miñoza et al., 2022). The training set used was quantile-normalised and batch-corrected during the preprocessing step. The "limma" package's linear models with empirical-Bayes moderation of the variance were used to find differentially expressed genes (DEGs) in rheumatoid arthritis (RA) and healthy control (HC) samples. This method was performed in this study as the number of male samples is very limited (Umezawa et al., 1995). The limma fits the RA vs HC for each gene and computes the log2FC values, whereas empirical Bayes fixes the unstable variance and computes the t-value and p-value. In this study. the limma with eBayes is used.

Moreover, as the datasets were from different arrays such as \[HG-U219\] Affymetrix Human Genome U219 Array and \[HG-U133_Plus_2\] Affymetrix Human Genome U133 Plus 2.0 Array

"arrayWeights" were also applied (Ritchie et al., 2006).

Across the different arrays in the microarray data, a poor-quality array is a result of poor hybridisation or degraded RNA, which gives false results due to noise. In order to prevent noisy results, in this study "arrayWeights" are applied, which estimate one weight per array across all the ~15000 genes. It identifies the array residual larger than the model's expectation. It is assigned a value below 1, then limma uses those weights and flags the noisy array. (Ritchie et al., 2006) This ultimately leads to less contribution of those flagged arrays to estimate the log fold change.

Three reasons to perform arrayWeights in this study are - the two datasets are from two different platforms, GPL570 and GPL13667 arrays. Quantile normalisation and ComBat align the distribution of the individual gene sets but are unable to solve the issue related to noisy arrays. Another reason is the limited number of male samples in this study.

For the RA coefficient (groupRA) to estimate the log₂ fold-change of RA in relation to healthy control (HC), a design matrix was built as (~ group) for each comparison, with HC defined as the reference level. "lmFit" was used to fit gene-wise linear models, and "topTable" (ranked by p-value, genes) was used to extract moderated t-statistics, p-values, and log₂ fold-changes from "eBayes". According to Benjamini and Hochberg (1995), a gene was deemed significantly differentially expressed if it met both an effect-size and a statistical threshold: |log₂ fold-change| > 0.1 and a Benjamini-Hochberg false-discovery-rate (FDR)-adjusted p-value &lt; 0.05. The term "up in RA" (log₂FC &gt; 0) or "down in RA" (log₂FC < 0) was used to describe the direction of change (Rosati et al., 2024). Three configurations of the analysis were performed at identical thresholds: a female-only contrast, a male-only contrast, and a pooled contrast using all training samples. Each is an RA-versus-HC comparison computed *within* the stratum concerned, and no gene is described as sex-specific on the basis of appearing in one stratum's list and not the other's (Gelman & Stern, 2006); that question is addressed only by the interaction test of §2.10.

**Fold-change threshold: what it does and does not achieve.** The rule above applies the |log₂FC| filter *after* testing against the null hypothesis log₂FC = 0, which is the conventional but not the statistically correct way to demand a minimum effect size. As a sensitivity analysis the same contrasts were therefore refitted with limma's `treat`, which places the fold-change threshold inside the null hypothesis (H₀: |true log₂FC| ≤ τ) and is the formally correct and deliberately conservative procedure (McCarthy & Smyth, 2009). Thresholds of τ = 0.1, log₂(1.2) and 0.5 were evaluated and the resulting counts compared with the post-hoc-filtered counts.

The comparison must be reported rather than relegated, because it bears on how much prioritisation this stage performs. A |log₂FC| threshold of 0.1 corresponds to a 1.07-fold change, which is a permissive effect-size gate, and formal `treat` testing at a 1.2-fold threshold reduces the significant-gene counts by more than an order of magnitude. The practical consequence is that the intersection with the disease modules in §2.5 retains a large majority of the disease-module genes, so **this stage is better described as a coherence filter than as a prioritisation step**, and the prioritisation burden falls on the Mendelian randomisation of §2.6 and the feature selection of §2.8. The permissive threshold was retained for the primary analysis so that genes of modest effect remained eligible for the causal screen — the male stratum is small and its power correspondingly limited — but the `treat` counts are reported alongside, and no claim is made that the DEG stage substantially narrowed the candidate space.

**2.4 Co-expression network analysis — WGCNA**

The DEG analysis is used to identify changes in the expression pattern in different sample groups (HC vs RA), and it only considers a limited number of pairwise comparisons. The major limitation is that it does not capture the co-expression pattern between the genes. So, to overcome this limitation, a co-expression network analysis package called as WGCNA was developed, which allows us to identify modules; these are the group of genes which has similar behaviour and are co-expressed. (Zhang et al., 2005)

The weighted gene co-expression network analysis (WGCNA) framework is developed based on a network of genes. In this package, every genes are nodes and the two nodes are connected to each other by edges. Modules are the highly interconnected genes (nodes) whose expression pattern vary together across individuals. In this study, after performing DEG, WGCNA was performed to find modules of co-regulated genes linked to rheumatoid arthritis (RA) on the 70% training set (quantile-normalised, ComBat-corrected expression), and the 30% testing data was not loaded in the network as it is considered internal validation data. The analyses were performed with R version 4.4.2 with "WGCNA" and "clusterProfiler" 4.14.6. A fixed seed randomSeed = 1234 was passed to blockwiseModules.

Initially, the training data were screened using "goodSamplesGenes", which removes zero variance or too many missing values. Then, using a data-driven cut height mean and three standard deviations of the merge heights, sample clustering (average-linkage hierarchical clustering on Euclidean distance) was examined for outlier detection. Arrays beyond the clusters were removed.

A weighted co-expression network is entirely derived from correlations between samples. Consequently, a single noisy array can simultaneously disrupt gene pairs, resulting in a co-expression structure that mirrors the artefact rather than the underlying biology. Co-expression analysis is significantly more susceptible to this issue compared to differential expression analysis, where a noisy array increases the residual variance of individual genes and is further mitigated by the array quality weights. Therefore, detecting and retaining outliers is not a conservative approach but rather an uncorrected source of erroneous modules. The average linkage method mitigates the chaining effect observed in single linkage, which tends to incorporate isolated arrays into the main cluster. Unlike Ward's method, it does not assume clusters of equal size and spherical shape. These characteristics are unsuitable when the objective is to identify a few isolated arrays in contrast to a large group. Moreover, the conventional approach is to use three standard deviations, which is more conservative and results in fewer arrays being removed compared to using two standard deviations.

In this study, all the genes were considered, and no variance pre-filter was performed. This is because the variance filter discards the genes that would contribute to the module. Here, the sole purpose of performing WGCNA is to get the disease-module that is intersected within each sex's DEG list to define the candidate genes carried forward to the Mendelian randomisation. Therefore, a gene removed before the module detection can never contribute to the intersection and can never be tested for a causal relationship with RA.

Initially, the WGCNA pipeline was applied in order to perform the gene variance filtering, and it resulted in the removal of many genes. Moreover, the male sample are also very limited. Therefore, in order to compensate this loss, all the genes across 173 samples were taken as a single block (maxBlockSize = 20000) and the modules were detected globally rather than in blocks. The genes with lower information are not ignored; instead, they are excluded by the network. A soft-thresholding at the chosen power drives the weak correlation towards zero, and these are considered as unassigned modules (grey).

**Correlation**

A Pearson correlation was used in this study (corType = "pearson"). According to Langfelder & Horvath (2008), bicor was suggested instead of Pearson. The implementation and evaluation of bicor were conducted. However, it was subsequently reverted due to computational constraints. Unlike Pearson correlation, which possesses a closed-form solution, bicor determines its weights through an iterative process, necessitating additional data copies and incurring approximately three times the computational cost. This poses a challenge when combined with the decision above to retain all 15,763 genes. Specifically, a single 15,763 × 15,763 matrix requires approximately 2 GB of memory, and the topological-overlap calculation necessitates simultaneous storage of the adjacency matrix, the TOM, and intermediate data. On the analysis machine, which has 16 GB of memory, the bicor process entered swap and failed to complete, whereas the Pearson network was successfully processed in a single block. This is a computational constraint and not a claim that Pearson is better than bicor.

**Network type**

A network of signed co-expression was constructed. An unsigned network employs |cor|, resulting in a scenario where a gene pair exhibiting strong anti-correlation is considered strongly connected, potentially leading to their inclusion in the same module. This is problematic for disease-programme analysis, as a transcript that is up-regulated in rheumatoid arthritis (RA) and another that is down-regulated in RA are biologically opposed. Consequently, a module containing both lacks coherent direction, an interpretable eigengene, and cannot be characterised as either "up in RA" or "down in RA," which is essential for the directional consistency check outlined in §2.5. In contrast, a signed network retains the correlation sign, ensuring that modules are directionally homogeneous and each maintains a clear relationship to disease status.

Modules were defined from the topological overlap matrix (TOM) rather than directly from the adjacency matrix (Zhang & Horvath, 2005). The TOM measures the extent to which two genes share their networks in addition to being connected. This approach substantially suppresses edges arising from noise in individual correlations, as unlikely pairs do not typically share modules, and yields more robust and reproducible modules than clustering based solely on correlation.

**Soft-thresholding power**

Soft-thresholding raises the correlation to a power β so that strong correlations are preserved while weak ones are driven towards zero, providing a continuous alternative to a hard correlation cut-off, which would discard the magnitude information distinguishing a correlation of 0.9 from one of 0.5 (Zhang & Horvath, 2005). β is chosen so that the resulting network approximates scale-free topology. Scale-free topology fit and mean connectivity were evaluated across candidate powers (1–10, then 12–20 in steps of two) with pickSoftThreshold for a signed network.

**β = 12 was selected on the joint criterion of scale-free fit and network connectivity**, and the reasoning is given in full because the automatic estimate was rejected. pickSoftThreshold returns the *first* power clearing the R² threshold rather than the best one, and the power it returned produced a network in which the average gene was connected to a large fraction of all others — a mean connectivity of the order of thousands, far too dense for meaningful module resolution, and accompanied by a high proportion of unassigned genes and negligible pathway enrichment. The signed R² curve plateaus at approximately 0.90 from β = 9 onwards while mean connectivity falls steeply, so powers on the plateau achieve equivalent scale-free fit at a far better conditioned connectivity. β = 12 sits on that plateau. The full fit-versus-connectivity table at every candidate power is reported in the Results so that the choice can be audited.

The WGCNA FAQ recommends β = 12 for a signed network with more than 40 samples, which coincides with the value used. **This coincidence is reported as corroboration only and is deliberately not offered as the justification**, because the FAQ table is formally conditional on the scale-free fit *failing*. Here the fit did not fail — an estimate was returned — so the table does not formally apply, and citing it as the primary reason would misapply the authors' own guidance (Langfelder & Horvath, 2008).

**Module-trait association**

The four group-by-sex combinations, age, sex, and disease status (RA, HC), were all encoded in a binary trait matrix. Each characteristic was correlated with module eigengenes (corPvalueStudent, Pearson correlation with student asymptotic p-values). The same modules were recomputed within each sex separately (all samples, female, male) so that the disease association of every module can be studied in each sex.

**Module detection and parameters**

Modules were identified within a single block using blockwiseModules at β = 12, employing the dynamic tree cut algorithm as described by Langfelder, Zhang, and Horvath (2008). The parameters and their respective justifications are detailed in the **table** below

**Gene-level importance and hub genes**

Module membership (kME, the correlation between a gene and its module eigengene), gene significance (the correlation between a gene and RA status), and intramodular connectivity (signed adjacency, power = 12) were used to measure gene-level importance within the RA-associated modules (Zhang & Horvath, 2005).

Hub genes were identified as genes satisfying both a module-membership and a trait-significance criterion, |kME| > 0.8 and |GS| > 0.2, and were additionally characterised by their intramodular connectivity. Requiring both criteria means that a hub must be central to its module *and* associated with disease: kME alone would select genes central to a module irrespective of any relationship to RA, and GS alone would select disease-associated genes irrespective of network position, which is merely a weaker restatement of the differential-expression analysis. Connectivity was computed within the module rather than across the whole network, so that hub status reflects centrality within the disease programme rather than global connectivity. No protein–protein interaction network and no betweenness-centrality analysis were performed; hub status in this chapter rests on kME, GS and intramodular connectivity alone.

**Functional enrichment**

Gene Ontology (enrichGO, all ontologies) and KEGG pathway enrichment (enrichKEGG, organism = "hsa") were applied to the genes in the disease-associated modules using clusterProfiler, with p- and q-value cutoffs of 0.05 and Benjamini-Hochberg correction. Ultimately, the intersection of the RA-associated module genes and the significant DEGs was kept as the high-confidence "DEG", "WGCNA" gene set and independently annotated, serving as the input for the Mendelian-randomisation and subsequent candidate-gene phases. (Nguyen et al., 2025)

**2.5 Candidate gene identification (disease module ∩ sex-stratified DEG)**

The next step combined this network-level evidence with the sex-stratified differential expression results to identify candidate genes for rheumatoid arthritis (RA), building directly on the co-expression analysis presented in the previous section, where weighted gene co-expression network analysis (WGCNA) divided the training transcriptome into modules of tightly co-regulated genes and linked each module to disease status.

The goal of this integration was to maintain the sex resolution that is essential to this work while keeping only the genes that were concurrently supported by two independent lines of evidence: membership in a disease-associated co-expression module, which indicates participation in a coordinated, systemically dysregulated program, and significant differential expression within a given sex, which indicates an individual disease association.

This design introduces the sex-stratified expression contrast, using the identical disease-associated module set for both sexes. In order to attribute any downstream sex stratification to expression biology rather than an analytical asymmetry, female and male candidates are selected from a shared, disease-relevant gene background.

Using Pearson correlation with Student asymptotic p-values, each module eigengene the first principal component summarising the expression of the genes in a module was correlated with a binary indicator of disease status (RA versus control) to identify disease-associated modules from the overall (all-sample) module-trait association analysis. When a module met both a strict effect-size threshold and a strict significance condition an absolute correlation with RA status of at least 0.5 and an association p-value below 1 × 10⁻<sup>8</sup> was kept as disease-associated. While the conservative p-value threshold prevented false module-trait relationships, the correlation cutoff ensured that retained modules captured a significant portion of disease-related expression variance. The resulting backdrop included both up-regulated and down-regulated disease programs since modules that changed in either direction were eligible. The disease-module gene background utilised for both sexes was created by combining the genes from the chosen modules.

The list of significantly differentially expressed genes for each sex from the sex-stratified limma analysis was then intersected with this disease-module backdrop, separately for each sex. The intersection of the disease-module genes and the female DEGs was identified as the female candidate set, while the intersection of the same disease-module genes with the male DEGs was identified as the male candidate set. The sex-stratified expression contrast is the only factor causing divergence between the two lists, because the module background is identical for both sexes; any downstream difference is therefore attributable to expression biology rather than to an analytical asymmetry in which sex was assigned which module.

**A gene appearing in one sex's list but not the other is not thereby shown to differ between the sexes.** A gene may fail to reach significance in one stratum through reduced statistical power alone, and the male stratum is far the smaller of the two. No claim of sex-specificity is made anywhere in this chapter on the basis of such a set difference; the only evidence that can speak to that question is the formal diagnosis-by-sex interaction test of §2.10, whose null hypothesis is that the RA-versus-control effect on a gene is the same in both sexes (Gelman & Stern, 2006).

The co-expression module, the direction of dysregulation (up- or down-regulated in RA), the log₂ fold-change and false-discovery-rate-adjusted significance in the relevant sex-stratified contrast, and the module membership value (kME, the correlation between the gene and its module eigengene), which measures how central the gene is within its module, were all annotated for each retained candidate. A two-set Venn diagram for each sex was also used to visualise the overlap within each candidate list.

This resulted in female (disease-module genes intersection female DEGs) and male candidate genes (disease-module genes intersection male DEGs).

**Sex stratified network and module preservation**

Two networks were constructed at the soft-thresholding parameter (beta =12). The potential differences in co-expression structure between sexes were examined using modulePreservation (Langfelder et al., 2011). In this analysis, the female data served as the reference network, while the male data functioned as the test network, with 200 permutations conducted.

Preservation Z-summary statistics are constructed by permuting sample labels within the test set, generating the null distribution at the test set's own size. This controls for the fact that the male stratum is much smaller than the female stratum. A direct comparison of male and female module structures, or of male and female DEG lists, does not account for this and would report reduced power in the male stratum as a biological sex difference.

The reference module assignment was derived from the combined-network modules rather than the female-network modules. This distinction is crucial as female-network modules represent a different gene partition from the disease modules, even when colour names coincide. Thus, using female-network modules would address the question, "Are female-network modules preserved in men?" rather than the intended question, "Are the disease modules preserved between sexes?"

Interpretation thresholds are as follows: a preservation Z-summary below 2 indicates a non-preserved structure, between 2 and 10 suggests moderate preservation, and above 10 denotes strong preservation (Langfelder et al., 2011; Nguyen & Zeng, 2025).

Regarding the number of permutations, two hundred permutations were conducted, which is double the 100 used by Nguyen & Zeng (2025). This increase was implemented to stabilise the permutation null in the context of a small test set, as the male stratum is significantly smaller than either group in that protocol. Consequently, the Z-summary is more susceptible to permutation noise.

The direction of the comparison was unidirectional, with the female network serving as the reference and the male network as the test. Nguyen & Zeng (2025) also reverse the reference and test assignments to evaluate preservation bidirectionally. The unidirectional approach was selected here because the male stratum is too small to function as a stable reference network. The asymmetry of this comparison is acknowledged as a limitation.

**2.6 Mendelian randomisation**

**Rationale and design**

Differential gene expression identifies the genes that differ between control and treated samples, but the limitation of DEG is that it cannot distinguish the driver of the disease. In this study, Two-sample Mendelian randomisation (MR) was used for the causal anchoring step. These are germline genetic variants that regulate gene expression, known as expression quantitative trait loci (eQTLs). An eQTL is a locus that accounts for a portion of the genetic variance in a gene expression phenotype. Standard eQTL analysis involves a direct association test between markers of genetic variation and gene expression levels, typically measured in tens or hundreds of individuals. This association analysis can be conducted either proximally or distally to the gene. The regulatory variants can be studied as cis and trans. Variants within 1 Mb (megabase) on either side of the TSS site are called cis, while those at least 5 Mb downstream or upstream of TSS or on different chromosomes are considered trans. In this study, cis is performed as cis variant lies within the gene's own regulatory mechanism such as the promoter, gene body, enhancer, splice sites, etc and it is used as instrumental variables for that gene expression and association with Rheumatoid arthritis.

MR was strategically placed between differential-expression screening and machine-learning feature selection. The differentially expressed candidate genes served as the exposures tested, and only those genes that passed MR with a false-discovery rate below 0.05 were advanced as inputs for feature selection.

Associations between exposure and outcome were estimated using distinct samples from published GWAS and eQTL resources that surpass the power of any single cohort (Hartwig et al., 2016). This design presupposes that the two samples originate from the same underlying population without participant overlap; significant overlap biases the estimates towards the confounded observational association, with the bias magnitude increasing as instrument strength decreases (Burgess et al., 2016). The exposure data (eQTLGen, blood expression in population cohorts) and the outcome data (Okada et al., 2014, a case-control RA meta-analysis) are derived from separate European resources with no known participant overlap, and the instrument-selection threshold adopted below guarantees strong instruments, so overlap bias is expected to be negligible.

Whole blood was chosen as the eQTL tissue to match the tissue of the transcriptomic discovery datasets, ensuring that the instrumented exposure is the same quantity measured by the biomarker panels. This choice has limitations. Genetic regulatory effects are substantially tissue-specific (GTEx Consortium, 2020), and whole blood may not represent transcriptional regulation in synovium, the tissue of primary pathology in RA. A gene with a genuine causal effect exerted only in synovium could be missed, and a blood eQTL effect may not correspond in magnitude or direction to the synovial effect. The trade-off was accepted because no synovial eQTL resource of comparable size exists, and tissue-matching to the discovery data is what makes the MR estimate relevant to the biomarker panels being validated.

All GWAS and eQTL summary statistics were accessed programmatically through the MRC IEU OpenGWAS database (Elsworth et al., 2020) via TwoSampleMR (Hemani et al., 2018); no individual-level genotype data were used, and the study therefore required no additional ethical approval beyond that of the contributing studies.

Instrument selection

The exposures tested were the disease-associated differentially expressed genes carried forward from the sex-stratified screen. MR was run once on the union of the female and male candidate sets, and the results were then partitioned by stratum. Running the estimation once ensures that a gene appearing in both strata carries an identical instrument set and an identical estimate, so any difference between the two prioritized lists arises from gene eligibility alone rather than from estimation.

For each candidate gene, HGNC symbols were mapped to Ensembl gene identifiers (org.Hs.eg.db 3.20.0), and instruments were extracted with TwoSampleMR::extract_instruments() under the following criteria:

Genome-wide significance for the exposure: p < 5 × 10⁻⁸. This threshold satisfies the relevance assumption while controlling the family-wise error rate across the genome under the standard correction for the effective number of independent common variants in European populations (Pe'er et al., 2008). A more permissive threshold would admit weak instruments and, in a two-sample design, bias estimates towards the null (Burgess and Thompson, 2011).

Independence: LD clumping at r² < 0.001 within a 10,000 kb window (TwoSampleMR defaults, European 1000 Genomes Phase 3 reference panel; Auton et al., 2015). Correlated instruments violate the independence assumption of the inverse-variance weighted estimator, which treats instrument-specific estimates as independent observations. Clumping to near-zero r² is the standard remedy (Burgess et al., 2013; Hemani et al., 2018).

Instrument strength: F = (β_exposure / SE_exposure)², evaluated against the conventional weak-instrument threshold F ≥ 10 (Burgess and Thompson, 2011).

Cis restriction: the variant must lie on the same chromosome as the gene it instruments and within 1 Mb of the gene body (GRCh37). The 1 Mb window was chosen because it is the cis boundary used by eQTLGen itself when generating these summary statistics (Võsa et al., 2021), so the instrument definition matches the definition under which the exposure data were generated. Cis instruments are preferred over trans because a variant acting locally on its own transcript has a mechanistically constrained route to the outcome, whereas trans variants are, by definition, candidates for horizontal pleiotropy (Swerdlow et al., 2016; Zhu et al., 2016).

Regarding the F statistic, at a significance level of p < 5 × 10⁻⁸, the absolute value of Z is at least 5.451, resulting in F = Z² being at least 29.7 for each instrument by design. Consequently, the criterion of F ≥ 10 does not exclude any variant selected under the initial criterion. This criterion is reported to confirm that the instrument strength surpasses the conventional weak-instrument threshold, rather than serving as a selection step. Therefore, weak-instrument bias, which in a two-sample design tends to attenuate estimates towards the null (Burgess and Thompson, 2011; Hartwig et al., 2016), is not a significant concern.

Concerning the cis restriction, eQTLGen datasets available in OpenGWAS include both trans and cis associations, and the function extract_instruments() retrieves both types. Thus, the restriction is not inherently implied by the dataset choice and must be explicitly enforced. This is crucial because a trans variant located on a strong RA locus can instrument a distant gene while affecting the outcome through the locus it resides on, thereby violating the exclusion restriction by design rather than leaving it untested. Gene coordinates were sourced from EnsDb.Hsapiens.v75 (Ensembl 75 = GRCh37) to align with the genome build of the eQTLGen/OpenGWAS variant positions; a GRCh38 annotation would not be positionally comparable.

Regarding the major histocompatibility complex, instruments within the extended MHC (chr6:25-34 Mb, GRCh37; Horton et al., 2004) were identified but retained in the primary analysis and excluded in a comprehensive parallel sensitivity analysis, described below. Retention is deemed appropriate for a primary estimate. However, given that the extended MHC exhibits the most extensive long-range linkage disequilibrium in the human genome (de Bakker et al., 2006) and contains the primary RA susceptibility locus HLA-DRB1 (Raychaudhuri et al., 2012), MHC status is included as a column in every results table to ensure that no estimate is interpreted without this context.

Outcome extraction and harmonization

Outcome associations were obtained using the extract_outcome_data() function, processed in 250-SNP segments with LD proxy substitution enabled, allowing up to three retries per segment. This approach was necessary as a single unsegmented query for the complete instrument set exceeded the server timeout.

The harmonisation of exposure and outcome datasets was conducted using harmonise_data(action = 2). This method deduces the strand of palindromic variants based on allele frequency and excludes those that remain ambiguous at intermediate frequencies (Hemani et al., 2018).

Gene labels were reassigned to harmonised rows by matching on the (SNP, exposure dataset ID) pair, rather than solely on SNP. In cases where a variant serves as an instrument for multiple genes, it is assigned the label of the first gene listed, which may result in the silent reassignment of instruments between genes. The exposure dataset ID (eqtl-a-) is the definitive identifier for the gene under investigation.

Effect estimation

The estimator applied to each gene was determined by the number of instruments surviving harmonisation. For three or more instruments, Inverse-variance weighted (IVW; Burgess et al., 2013), MR-Egger (Bowden et al., 2015), weighted median (Bowden et al., 2016), for 2 IVW, for 1 Wald ratio.

The three estimators serve complementary roles, each relaxing distinct assumptions. The Inverse Variance Weighted (IVW) estimator is the most statistically efficient when all instruments are valid; however, it exhibits bias if any instrument is pleiotropic (Burgess et al., 2013). The weighted median estimator maintains consistency as long as at least 50% of the weight is derived from valid instruments (Bowden et al., 2016). The MR-Egger estimator permits all instruments to be pleiotropic, provided that pleiotropic effects are independent of instrument strength, as per the InSIDE assumption, though this results in significantly reduced precision (Bowden et al., 2015). Thus, concordance among the three estimators suggests that the estimate is not contingent upon the specific assumption employed.

In instances where multiple estimators were applied, a single primary estimate per gene was selected according to a predetermined hierarchy: IVW > Wald ratio > weighted median > MR-Egger. This hierarchy was established in advance to ensure that the reported estimate is not selected based on its p-value. IVW is prioritized because it is the most efficient estimator under the null hypothesis of no pleiotropy and is the standard primary analysis in two-sample Mendelian Randomization (MR), with the other estimators functioning as sensitivity analyses rather than as competing primary estimates (Burgess et al., 2013; Hemani et al., 2018; Sanderson et al., 2022).

Effects are expressed as odds ratios for rheumatoid arthritis (RA) per unit increase in normalized gene expression, accompanied by 95% confidence intervals calculated as exp(β ± 1.96 SE). Genes were categorized as risk-increasing (OR > 1) or protective (OR < 1) based on the false discovery rate (FDR)-adjusted p-value, rather than the nominal p-value.

Sensitivity analyses for pleiotropy and heterogeneity

For genes instrumented by three or more variants, Cochran's Q test and MR-Egger intercept test were used. In Cochran's Q test, the variation among estimates specific to each instrument serves as an indirect sign of pleiotropy (Bowden et al., 2015), and in the MR-Egger intercept test, the intercept that is not zero suggests the presence of directional (unbalanced) horizontal pleiotropy, as noted by Bowden et al. (2015).

Multiple-testing correction and sex stratified

P-values were adjusted using the Benjamini-Hochberg procedure at a false discovery rate (FDR) of less than 0.05 (Benjamini and Hochberg, 1995), applied separately within each sex stratum. The denominator for this adjustment is the set of genes tested within each stratum, rather than the combined set from both strata. In a sex-stratified design, using a combined denominator would transfer the multiplicity burden from one stratum to the other. A pooled correction is included in the output as a legacy column but is not reported.

The set of genes surviving the FDR threshold, rather than those meeting a nominal p-value of less than 0.05, was used for feature selection. At a nominal threshold, the expected false positive rate is 5% of the genes tested, which, given the size of the screen, constitutes a significant portion of the genes identified. Relying on this list would have allowed predominantly spurious associations to influence all subsequent model development.

The eQTLGen and all three rheumatoid arthritis (RA) genome-wide association studies (GWAS) utilized here are sex-combined. A review of all 37 RA datasets in OpenGWAS confirmed the absence of sex-stratified RA GWAS. Consequently, the Mendelian randomisation (MR) estimate for a given gene is numerically identical in both female and male outputs; the two strata differ only in the genes eligible for testing and in the FDR denominator.

The design thus involves sex-stratified discovery followed by sex-combined causal validation.

The most significant challenge to any causal assertion presented here pertains to the extended MHC. Two key factors render this region particularly inhospitable to cis-MR: firstly, HLA-DRB1 is overwhelmingly the primary RA susceptibility locus (Raychaudhuri et al., 2012); secondly, the extended MHC exhibits the most extensive long-range linkage disequilibrium within the human genome (de Bakker et al., 2006). Consequently, a cis-eQTL for any gene within this region is inherently correlated with the HLA-DRB1 risk haplotype, irrespective of the gene's causal involvement, thereby directly contravening the exclusion restriction. This is because the instrument influences RA through a pathway that bypasses the exposure. As a result, a cis-MR estimate for an MHC gene is indistinguishable from the HLA-DRB1 signal when interpreted via a proxy. This issue is further exacerbated for genes that are instrumented by a single variant, as there is no available pleiotropy or heterogeneity test to identify the violation.

A gene that is solely represented by an MHC variant has not been tested and deemed null; rather, it has become untestable. Including it in the denominator would unfairly disadvantage testable genes due to the lack of evidence regarding untestable ones.

A gene may fall below the false discovery rate (FDR) threshold in the sensitivity analysis for reasons that possess distinct evidential significance. Firstly, if the gene's own instrument set is altered, resulting in a change in its estimate, this indicates genuine major histocompatibility complex (MHC) dependence. Conversely, if the gene's instruments are entirely outside the MHC and its estimate remains numerically unchanged, the Benjamini-Hochberg procedure, being a rank-based method, may increase its adjusted p-value by removing the highly significant MHC genes ranked above it, despite no change in the gene's own evidence. This latter scenario pertains to multiplicity bookkeeping rather than confounding, and conflating these two scenarios would exaggerate the impact. These scenarios can be differentiated by testing whether the p-value and the instrument count remain unchanged within a numerical tolerance of 1 × 10⁻¹².

**2.7 Colocalisation**

A cis-MR estimate indicates that the eQTL for a gene is linked to RA. However, it does not confirm that the eQTL and the RA risk signal are influenced by the same causal variant. These two issues diverge when the eQTL variant is merely in linkage disequilibrium with a separate disease-causing variant. In such cases, cis-MR yields a confident, highly significant, yet entirely spurious causal estimate. This is the primary acknowledged limitation of cis-eQTL MR (Zhu et al., 2016; Wallace, 2020), which is why MR estimates within the MHC should not be accepted without scrutiny. Colocalisation transforms the implicit assumption made by MR into a measurable quantity.

Single Causal Variant

Bayesian colocalisation was conducted using coloc::coloc.abf (Giambartolomei et al., 2014) for each FDR-surviving prioritised gene and every gene in the final biomarker panels (10d_coloc_panel_genes.R). This method assesses five mutually exclusive hypotheses over a genomic region: H0 (no causal variant for either trait), H1 (causal variant for expression only), H2 (causal variant for RA only), H3 (both traits have a causal variant, but different variants), and H4 (both traits share one causal variant). H3 represents the failure mode that invalidates a cis-MR estimate, while H4 represents the scenario assumed by the MR estimate.

Full regional summary statistics were used instead of the clumped instruments, as colocalisation requires all variants in the window. Therefore, this step queries OpenGWAS rather than reading the cached instrument file. The window was defined as the gene body ± 250 kb using GRCh37 coordinates from EnsDb.Hsapiens.v75, the same coordinate source as the Section 3 cis filter, ensuring that the two analyses cannot disagree on gene location. A window narrower than the 1 Mb MR cis boundary was deliberately chosen because coloc.abf assumes at most one causal variant per trait within the analysed region, and this assumption becomes less tenable as the window widens to include additional independent signals (Giambartolomei et al., 2014; Wallace, 2021).

Harmonisation was conducted conservatively: variants were consolidated based on rsID, retaining only those present in both datasets. The GWAS effect allele was aligned with the eQTL effect allele, applying sign inversion of β when alleles were swapped. Strand-ambiguous variants (A/T, C/G) were excluded, as the Okada extract lacks allele frequency data, preventing strand resolution by MAF, and an incorrect flip would remain undetected. Variants with mismatched allele pairs after flipping (due to multi-allelic or annotation discrepancies) were also excluded. MAF was sourced from the eQTL dataset, with both studies being European. Regions with fewer than 50 shared variants were not analysed, as the posterior in such cases is dominated by the prior rather than the data, rendering the resulting probabilities uninterpretable (Wallace, 2020).

Regarding prior sensitivity, the p12 prior is the probability that a variant is causal for both traits that serves as the primary subjective input and is typically scrutinised (Wallace, 2020). Consequently, each gene was analysed using the default p12 = 1 × 10⁻⁵ and a conservative p12 = 1 × 10⁻⁶, with both results reported.

The assumption of a single causal variant per trait per region is challenged in certain contexts. The coloc.abf method presumes a maximum of one causal variant per trait within a given region, which is a reasonable assumption for a typical cis-eQTL window. However, this assumption does not hold true in the Major Histocompatibility Complex (MHC), where rheumatoid arthritis (RA) exhibits multiple independent and well-mapped signals, such as those at HLA-DRB1 positions 11/71/74, HLA-B position 9, and HLA-DPB1 position 9 (Raychaudhuri et al., 2012). In the MHC, regions often contain numerous variants in extended linkage disequilibrium (LD). When this assumption is violated, it results in an inflated posterior probability (PP.H3), as two traits influenced by several variants may appear to have "different causal variants" even if they share one. Consequently, a high PP.H3 within the MHC should not be interpreted as definitive evidence of distinct causal variants. The appropriate and more cautious conclusion is that within the MHC, both cis-MR and coloc.abf are unreliable, rendering any causal claims unsupportable in either direction. MHC genes are identified in an assumption_valid column and reclassified from "DISTINCT VARIANTS" to "MHC - unreliable both ways." Outside the MHC, the assumption remains valid, and PP.H3 retains its conventional interpretation.

Multiple Causal Variants in the MHC

To investigate the MHC, colocalisation was re-evaluated, allowing for multiple causal signals per trait (10e_coloc_susie_mhc.R). The "Sum of Single Effects" model, SuSiE (Wang et al., 2020), decomposes each trait's regional association into several credible sets, each representing an independent causal signal. coloc.susie (Wallace, 2021) subsequently assesses colocalisation between each pair of credible sets, comparing eQTL signal i with GWAS signal j, and assigns a posterior probability to each pair. Consequently, a gene may colocalise with one of several RA signals while remaining independent of others-a scenario that coloc.abf cannot depict, as exemplified by the MHC. SuSiE was implemented with a maximum of 10 credible sets per trait (L = 10) at 95% coverage, which are the package defaults (Wang et al., 2020); L = 10 comfortably exceeds the number of independent RA association signals identified within the MHC (Raychaudhuri et al., 2012), ensuring that this limit does not constrain the outcome. Regional statistics were retrieved from the cache created in the single-causal-variant analysis above, obviating the need for re-downloading.

Requirements and Limitations of coloc.susie Compared to coloc.abf. coloc.susie requires a regional LD matrix aligned with the allele orientation of the summary statistics, sourced here from the 1000 Genomes European reference panel using ieugwasr::ld_matrix() and adjusted to match the effect alleles of the summary statistics. Two specific consequences are noted:

Reference-panel mismatch. The LD matrix is out-of-sample (1000 Genomes EUR, n ≈ 500) rather than derived from the eQTLGen or Okada samples. In the MHC, where haplotype structure is highly complex, an out-of-sample LD reference can introduce errors and alter credible sets (Zou et al., 2022). An in-sample LD would be necessary to accurately resolve the region, but neither eQTLGen nor Okada provide it.

Variant cap. The ld_matrix is limited server-side, so each region was narrowed down to its 450 most strongly associated shared variants (based on the smaller p-values of the two traits). Reducing a region can result in the omission of a causal variant.

A gene can be considered to carry a causal claim only if it satisfies two conditions: (i) it is classified ROBUST in the MHC sensitivity analysis of §2.6, and (ii) it is supported by a PP.H4 value of at least 0.80 at both priors in §2.7. If a gene meets only the first condition, it is reported as associated with, rather than causal for, rheumatoid arthritis (RA).

A PP.H3 value of at least 0.80 outside the MHC provides positive evidence that the MR estimate is LD-confounded, rather than merely indicating an absence of evidence, and is reported accordingly. Genes classified as inconclusive lack sufficient power for colocalisation. This is because Okada (2014) predates contemporary imputation panels, and the limiting factor is regional variant density rather than the method itself.

**2.8 Feature selection — LASSO, SVM-RFE, Random Forest**

2.8 Sex-stratified feature selection and biomarker panel derivation

2.8.1 Candidate genes entering feature selection

To refine the sex-stratified candidate genes into compact diagnostic biomarker panels, supervised machine-learning feature selection was employed. Only genes with genetic evidence indicating a directional effect on rheumatoid arthritis (RA) liability were considered eligible. The candidates entering this phase were those whose cis-eQTL-instrumented expression (Võsa et al., 2021) was linked to RA risk in the Mendelian randomisation (MR) analysis against the European-ancestry GWAS of Okada et al. (2014), with a Benjamini-Hochberg false discovery rate below 0.05 within each sex's stratum (Benjamini and Hochberg, 1995). Genes were retained regardless of the direction of the MR estimate, allowing both risk-increasing (OR > 1) and protective (OR < 1) genes to be included. Of the 2,045 female and 2,079 male disease-associated candidates, 1,477 and 1,478, respectively, yielded usable cis instruments and were tested, with 32 female and 25 male genes surviving correction (§2.6).

Three qualifications are explicitly noted. First, since both the eQTLGen exposure data and the Okada outcome GWAS are sex-combined, MR was estimated once across the union of the two candidate sets and subsequently partitioned by stratum. Consequently, a gene present in both lists carries an identical estimate, and the difference between the female and male candidate sets reflects upstream gene eligibility and the within-stratum correction denominator, rather than sex-specific genetic effects. Second, colocalisation analysis (§2.7) later demonstrated that several of these genes are influenced by variants distinct from the RA association signal. They are therefore described throughout as genetically prioritised rather than causal genes, and no member of the final panels is attributed with a causal claim.

Third, and bearing directly on the resampling estimates of §2.9, the eligible set is not an external gene list. The genes *submitted* to Mendelian randomisation were the candidates of §2.5, derived from the training expression matrix using the training diagnosis labels; only the *filter* applied to them was external. The eligible set is thus an externally filtered internal list, and it is not described as free of selection bias merely because its final filter was genetic. The consequences for the performance estimates, and the reasons the residual optimism is bounded, are set out in §2.9.

2.8.2 Training data and design

Feature selection was conducted solely on the training cohort, with analyses for females and males maintained entirely separate. Each model was fitted within a single sex, ensuring no exchange of samples, folds, or fitted parameters between them.

Expression values from GSE93272 and GSE110169 were log₂-transformed upon import and merged based on the common gene set. The merged, pre-normalisation matrix was then divided into a 70:30 ratio, forming a training set and a sealed internal-validation holdout. This partitioning employed stratified random sampling based on the combination of dataset, diagnosis, and sex (caret::createDataPartition; Kuhn, 2008). The decision to split prior to any cross-sample operation was intentional: quantile normalisation requires a shared reference distribution across samples, and ComBat estimates batch and covariate effects using outcome labels. Estimating these on pooled data could allow the holdout to influence the model's input representation, leading to an optimistic bias in performance (Ambroise and McLachlan, 2002; Simon et al., 2003; Kaufman, Rosset and Perlich, 2012). Consequently, quantile normalisation (Bolstad et al., 2003, as implemented in limma; Ritchie et al., 2015) and empirical-Bayes batch correction by ComBat, which protects diagnosis and sex in the model matrix (Johnson, Li and Rabinovic, 2007; Leek et al., 2012), were estimated solely on the 70% training partition and subsequently applied to the holdout as frozen parameters. The resulting training cohort consisted of 183 samples: 145 female (86 RA, 59 healthy control) and 38 male (17 RA, 21 healthy control).

For each sex, a design matrix was constructed with samples as rows and genetically prioritised candidate genes for that sex as columns. The response variable was diagnosis, coded as RA versus healthy control, with healthy control serving as the reference level. Three complementary supervised feature-selection algorithms were applied within each sex, with the primary hyperparameter of each tuned through internal 10-fold cross-validation on the training data. A single global seed (1234) was set before every stochastic step to ensure that fold assignment, forest growth, and cross-validated tuning are exactly reproducible.

2.8.3 LASSO logistic regression

The LASSO logistic regression model was implemented using the glmnet package (Friedman, Hastie, and Tibshirani, 2010), with the family set to "binomial" and alpha set to 1 to apply a pure L₁ penalty (Tibshirani, 1996). This L₁ penalty facilitates embedded selection by reducing non-informative coefficients to zero. The penalty parameter λ was optimized through 10-fold cross-validation, employing binomial deviance as the loss function (cv.glmnet, type.measure = "deviance"). Genes with non-zero coefficients at λ_min, the value that minimizes cross-validated deviance, were included in the LASSO selection; the more parsimonious λ_1se solution was also documented. The tuning process yielded λ_min values of 0.0434 for females and 0.0204 for males, resulting in the selection of 7 and 10 genes, respectively.

2.8.4 Random forest importance

A random forest model was developed using the randomForest package (Breiman, 2001; Liaw and Wiener, 2002), comprising 1,000 trees, an increase from the default 500 to enhance the stability of importance estimates. The mtry parameter, representing the number of variables sampled as candidates at each split, was optimized through a 10-fold cross-validated grid search aimed at maximizing the area under the receiver operating characteristic curve (caret::train, method = "rf", metric = "ROC"; Kuhn, 2008). The grid included 1, 2, ⌊√p⌋, ⌊p/3⌋, ⌊p/2⌋, and p candidate genes, with duplicate values removed, resulting in {1, 2, 5, 10, 16, 32} for females and {1, 2, 5, 8, 12, 25} for males. The final 1,000-tree forest was refitted at the optimized mtry value (mtry = 2 for females, 1 for males). A gene was selected if its mean decrease in Gini impurity exceeded the average mean-decrease-in-Gini across all candidate genes for that sex, retaining only those genes that contributed more than the average candidate to node-impurity reduction; this process identified 10 genes for each sex.

2.8.5 Support vector machine recursive feature elimination

Recursive feature elimination was conducted using a linear-kernel support vector machine (Cortes and Vapnik, 1995) via the e1071 interface to LIBSVM (Chang and Lin, 2011), following the methodology of Guyon et al. (2002). The regularization parameter C was initially optimized through a 10-fold cross-validated grid search over the set {0.01, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16}, resulting in C = 0.01 for both sexes. At this optimized cost, a linear SVM was trained on all candidate features, and the feature with the smallest squared weight was iteratively removed, with the model being refitted each time, until only a single feature remained, thereby producing a complete importance ranking. The panel size was then objectively determined as the number of top-ranked features that minimized the 10-fold cross-validated classification error, resulting in the selection of 8 genes for females and 11 for males.

The primary sex-stratified panel was established as the intersection of the three selections. Consequently, each panel member was independently identified through a penalised regression, a tree ensemble, and a margin-based selector.

**2.9 Machine-learning-based diagnostic model development**

Sex-stratified consensus panels were used to construct and evaluate a diagnostic classifier for rheumatoid arthritis (RA), with females and males treated as entirely separate analyses. No sample, cross-validation fold, tuning parameter, or fitted coefficient was exchanged between the two strata, ensuring that the female panel was evaluated only in females and the male panel only in males. No estimate reported in this chapter derives from a model that has seen data from both sexes.

A multivariable logistic regression model (glm, family = binomial) was developed for each sex to distinguish RA from healthy controls, with healthy control as the reference level and the standardized consensus panel genes entered as predictors. Additionally, a second penalized classifier was fitted over the entire within-sex set of genetically prioritized candidate genes, allowing the model to retain correlated but informative genes that a strict three-way intersection would necessarily discard (Zou and Hastie, 2005). This was an elastic-net logistic regression in which the mixing parameter was tuned over the grid α ∈ {0.1, 0.3, 0.5, 0.7, 0.9, 1.0}. For each value of α, the penalty strength λ was tuned by inner five-fold cross-validation on binomial deviance, with the pair (α, λ) minimizing the inner cross-validated deviance being selected jointly and predictions taken at λ_min (Tibshirani, 1996; Friedman, Hastie, and Tibshirani, 2010). The two classifiers are distinct models and are reported separately rather than pooled.

To ensure compatibility across microarray platforms, each gene was standardized to a z-score prior to model integration, under two distinct regimes that are not interchangeable and are delineated here for clarity. Within any resampling procedure conducted within a single cohort, the mean and standard deviation of each gene were calculated solely from the training fold and subsequently applied, unchanged, to the held-out fold, ensuring that no test sample influenced the standardization. In instances where a gene's standard deviation was zero or undefined, it was substituted with unity, resulting in the gene being centered but unscaled, thereby avoiding an undefined value. When a locked model was transferred to a different cohort, each gene was independently standardized within that target dataset. This approach was necessary due to the different array platforms used for the internal validation holdout and the external cohorts, which resulted in their absolute intensities not being on a common scale. This method is legitimate as the transformation is entirely unsupervised, relying solely on the expression values of the target dataset and not its diagnostic labels (Kaufman, Rosset, and Perlich, 2012). A panel gene absent from a target platform's annotation was assigned a z-score of zero, equivalent to the dataset mean, ensuring its contribution to the linear predictor was nullified rather than discarding the affected samples. The genes missing from each dataset were documented alongside the corresponding performance estimate.

Two quantities were reported on the training partition, distinctly labeled throughout due to their non-comparability. The apparent area under the receiver operating characteristic curve (AUC) was derived by fitting the locked panel to all training samples and scoring those same samples without any resampling; since the panel was selected using those samples, it serves as an optimistic upper anchor by design and is never presented as a validation estimate. The cross-validated figure quoted for the training cohort was obtained under the nested design described below, with z-scoring recalculated within each cross-validation training fold.

The consensus panel was selected utilizing the entire training set. Consequently, a conventional cross-validated AUC computed on this same set would exhibit optimistic bias, as each fold would evaluate genes previously selected with the aid of the samples being withheld. This phenomenon, known as feature-selection bias, is sufficiently pronounced in gene-expression data of this dimensionality to produce seemingly excellent classifiers from noise alone (Ambroise and McLachlan, 2002; Simon et al., 2003). To obtain an unbiased estimate, the entire selection and fitting process was incorporated within a nested cross-validation framework, stratified by outcome within each sex (Varma and Simon, 2006; Cawley and Talbot, 2010; Krstajic et al., 2014).

The outer loop employed a repeated stratified k-fold design: ten folds repeated five times for females (n = 145) and five folds repeated ten times for males (n = 38). Fewer folds were utilized in the male stratum because a ten-fold split at that sample size results in approximately four samples per test fold, with even fewer being cases, leading to substantial changes in fold-level AUC with the movement of a single sample; this instability was mitigated by increasing the number of repeats, thereby averaging over fold assignment rather than relying on a single arbitrary partition (Bengio and Grandvalet, 2004; Krstajic et al., 2014). Outer fold assignment was deterministically seeded with a fixed base value of 1000, incremented by the repeat index, and all other stochastic steps were controlled by a single global seed of 1234, ensuring that all reported estimates are precisely reproducible.

Within each outer training fold, without access to the outer test fold, the complete three-selector procedure was re-derived. LASSO logistic regression was refitted using coordinate descent (cv.glmnet, α = 1, family = binomial), with λ tuned through inner five-fold cross-validation. Genes retaining non-zero coefficients at λ_min were selected (Tibshirani, 1996; Friedman, Hastie, and Tibshirani, 2010). A random forest comprising 500 trees was constructed, with the number of predictors sampled at each split set to the package default of ⌊√p⌋. Genes whose mean decrease in Gini impurity exceeded the arithmetic mean of that statistic across candidates were selected (Breiman, 2001; Liaw and Wiener, 2002). Recursive feature elimination was conducted using a linear-kernel support vector machine with a fixed cost C = 1, eliminating at each iteration the feature with the smallest squared weight and retaining the panel size that minimized the inner five-fold cross-validated classification error (Cortes and Vapnik, 1995; Guyon et al., 2002; Chang and Lin, 2011). The panel for that fold was the intersection of the three selections, provided it contained at least two genes. If not, the procedure reverted, in a fixed and pre-specified order, to the union of the three selections, then to the LASSO selection alone, and finally to the complete candidate set, ensuring every fold yielded a usable panel and none was discarded.

These in-fold selectors replicate the structure of the selection procedure applied to the full training set (§2.8) but intentionally exclude its hyperparameter tuning. On the full training set, the LASSO penalty was tuned by ten-fold cross-validation on binomial deviance, resulting in λ_min = 0.0434 for females and 0.0204 for males. The forest was grown with 1,000 trees, and the number of predictors per split was tuned by ten-fold cross-validation over the grid {1, 2, ⌊√p⌋, ⌊p/3⌋, ⌊p/2⌋, p}, using cross-validated AUC as the tuning metric, yielding mtry = 2 for females and mtry = 1 for males. The support-vector cost was tuned by ten-fold cross-validation over the grid {0.01, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16}, resulting in C = 0.01 for both sexes, with the panel size derived from a ten-fold cross-validated error curve (Kuhn, 2008). Within the folds, however, the forest was grown with 500 trees without tuning mtry, and the support vector machine used a fixed unit cost, with only λ remaining self-tuned, as re-tuning every hyperparameter within every fold of every repeat was computationally prohibitive at this design size. This simplification is conservative, as an untuned selector is less likely than a tuned one to recover the panel obtained on the full training set. Consequently, the resulting quantity is described as an estimate of the pipeline under fixed selector hyperparameters rather than of the fully tuned procedure.

A logistic regression classifier was fitted to each fold's panel, standardized using the fold-training mean and standard deviation, as previously described, and subsequently applied to predict the untouched outer test fold. In the elastic-net variant, the penalized model was refitted within each outer training fold, with α and λ re-tuned through inner cross-validation on that fold's data alone. The out-of-fold predicted probabilities were retained for each sample. Since each sample is predicted exactly once within a repeat, an AUC was computed for each repeat and summarized as a mean and standard deviation across repeats. For the headline estimate, the out-of-fold probabilities of each sample were averaged across repeats, and a single receiver operating characteristic curve was constructed from these averaged probabilities. The case direction was predetermined rather than selected to maximize the AUC. Areas under the curve and their 95% confidence intervals were computed using pROC (Robin et al., 2011). Intervals were obtained by the method of DeLong, DeLong, and Clarke-Pearson (1988) and replaced by a stratified bootstrap interval with 2,000 resamples wherever a stratum comprised fewer than twenty samples, as the DeLong interval is unreliable and may report spurious precision in such cases (Carpenter and Bithell, 2000). Paired comparisons between panels were conducted on the same samples using the DeLong test for correlated curves.

Three quantities were computed on the training data to make the influence of the resampling design visible rather than assumed. The apparent AUC, as defined above, serves as the unambiguously optimistic upper anchor. The flat cross-validated AUC selects the panel once on the entire training set and cross-validates only the logistic model, thus explicitly representing the biased procedure that the nested design aims to displace. The nested AUC re-derives selection within every fold. Optimism was reported as the apparent estimate minus the nested estimate, following established practice in the internal validation of prediction models (Harrell, Lee, and Mark, 1996; Steyerberg et al., 2001). The difference between the flat and nested estimates was deliberately not reported as the magnitude of selection bias, because two effects act upon it in opposing directions: the flat estimate is inflated by the panel having been chosen using all of the data, whereas the nested estimate pools out-of-fold probabilities across repeats in which each fold selected its own panel, so that averaging those probabilities constitutes a form of ensembling and tends to raise the AUC (Dietterich, 2000). Therefore, the two estimates do not form a clean decomposition of bias, and both are tabulated with this caveat rather than differenced.

As the panel is re-derived within every fold, the frequency with which each gene of the reported panel re-entered the in-fold three-way consensus provides a direct measure of the reproducibility of the selection under resampling, and was reported as a percentage of fitted folds for each panel gene (Meinshausen and Bühlmann, 2010). A gene selected on the full training set but seldom recovered within folds was considered an unstable selection rather than a validated biomarker.

**Model evaluation**

**The three estimator families, and why they are named separately.** Several analyses in this chapter report a quantity that a reader would reasonably take to be "the area under the curve of the panel in dataset X", and they do not agree, because they estimate three structurally different things. They are named distinctly, every table states which family it belongs to, and values from different families are never compared with one another or averaged.

A **locked-transfer** estimate applies the model fitted and frozen on the training partition — its panel, its coefficients and its standardisation rule — to a dataset that contributed to none of them, and scores it once. This is the estimator that answers whether the diagnostic model works on new patients, and every claim of transportability in this chapter rests on it.

A **within-dataset resampled** estimate takes only the *identity* of the panel genes from training and re-estimates the classifier's coefficients inside the evaluation dataset under cross-validation. It answers whether the gene set carries discriminative information in that dataset, which is a weaker and different question. It is the appropriate estimator where transported coefficients would be meaningless because the target data are on a different scale or platform — as in the synovial analysis of §2.11 — or where the panel must be compared like-for-like against a competing model fitted in the same dataset, as in the cell-composition benchmark. It is never a test of transfer, because the coefficients were fitted on the samples being scored, and it is never quoted as out-of-sample performance of the locked model.

A **nested** estimate re-derives the panel-building procedure inside every outer training fold of the training partition and estimates the performance of the *procedure* rather than of any one fitted model.

The consequence is stated so that it cannot be lost: a within-dataset resampled estimate and a locked-transfer estimate computed on the same samples answer different questions and will differ, sometimes substantially, and neither is a correction of the other. Where both are available, both are reported and labelled by family, and the locked-transfer value is the one carried into any statement about diagnostic performance.

**Orientation conventions for per-gene evaluation.** Because the area under the curve of a single gene depends on which direction of expression is taken to mark disease, and because that choice may be made either before or after seeing the evaluation data, two conventions arise and are distinguished throughout. Under the **training-fixed** convention the direction is determined on the training data and applied unchanged, so that a gene whose association reverses out of sample yields a value below 0.5; this is the convention used for the internal holdout, the external blood cohort, and every table or figure that places several datasets side by side, and it is the only convention under which a failure of transfer remains visible. Under the **best-direction** convention the orientation is chosen within the evaluation dataset itself, so every value is at least 0.5 by construction; this answers only how much information a gene carries within that dataset, cannot support a claim of transfer, and is confined to single-dataset displays, where the direction of effect is instead carried by an explicit concordance flag. Every table and figure states which convention it uses.

In the training cohort, the panel model was evaluated using two metrics, neither of which provides an unbiased estimate of performance: the apparent AUC, derived through resubstitution, and the flat cross-validated AUC, obtained from a stratified ten-fold split where the panel remained fixed and only the logistic model underwent cross-validation. During this procedure, z-scoring was recalculated within each cross-validation training fold to prevent test samples from influencing their own standardization. Despite this, since the panel was selected using all training samples, each fold still evaluates genes chosen with the aid of the samples being held out. Consequently, the flat estimate is reported as an explicit representation of feature-selection bias rather than as a validation figure (Ambroise and McLachlan, 2002; Simon et al., 2003). The reliable training estimate is the nested cross-validated AUC described below.

Given that the panel was selected using the entire training set, a standard cross-validated AUC would exhibit optimistic bias. Consequently, the entire pipeline was re-evaluated through nested cross-validation, stratified by outcome within each sex (Varma and Simon, 2006; Cawley and Talbot, 2010; Krstajic et al., 2014). The outer resampling employed a repeated stratified k-fold scheme: 10 folds repeated 5 times for females and 5 folds repeated 10 times for males. This approach was chosen because, in the smaller male stratum, a ten-fold split results in only a few samples per test fold, with even fewer cases, necessitating additional repetitions to average over fold assignment rather than relying on a single arbitrary partition (Bengio and Grandvalet, 2004). Each repetition was reproducibly seeded (seed = 1000 + repeat index) under a single global seed of 1234, and folds were constructed using the caret package (Kuhn, 2008). Within each outer training fold, and without accessing the outer test fold, feature selection was conducted anew using fixed hyperparameters: LASSO (cv.glmnet, alpha = 1, inner 5-fold cross-validation, genes with non-zero coefficients at λ.min); random forest (ntree = 500, mtry at the default ⌊√p⌋, genes with above-mean Gini importance); and SVM-RFE (linear kernel, cost = 1, panel size determined by inner 5-fold cross-validated error) (Tibshirani, 1996; Breiman, 2001; Guyon et al., 2002; Liaw and Wiener, 2002; Friedman, Hastie and Tibshirani, 2010; Chang and Lin, 2011).

A consensus was reached among the three methods, incorporating a defined fallback strategy that prioritizes the union, followed by LASSO, and finally the full candidate set, to ensure each fold yields a viable panel. These in-fold selectors replicate the structure of the full-training-set procedure, excluding hyperparameter tuning due to computational constraints at this design scale. This simplification is conservative, as an untuned selector is not more likely than a tuned one to recover the reported panel, thus the resulting quantity is considered an estimate of the pipeline under fixed selector hyperparameters (Cawley and Talbot, 2010). Subsequently, a logistic-regression classifier was applied to the fold's panel, standardized using the fold-training mean and standard deviation, with any zero or undefined standard deviation replaced by unity, and used to predict the untouched outer-test fold. Out-of-fold predicted probabilities were aggregated across all folds and repetitions to generate a pooled AUC with a 95% DeLong confidence interval, which was substituted by a stratified bootstrap interval with 2,000 resamples when a stratum contained fewer than twenty samples (DeLong, DeLong, and Clarke-Pearson, 1988; Carpenter and Bithell, 2000; Robin et al., 2011). This was accompanied by the per-repeat AUC distribution and the re-selection frequency of each panel gene across folds, indicating panel stability (Meinshausen and Bühlmann, 2010). Optimism was calculated as the apparent estimate minus the nested estimate (Harrell, Lee, and Mark, 1996; Steyerberg et al., 2001); the flat and nested estimates were not differenced, as pooling out-of-fold probabilities across repeats, where each fold selects its own panel, constitutes a form of ensembling that increases the AUC, thus preventing a clear decomposition of bias (Dietterich, 2000).

**The candidate gene set was held fixed across folds, and the consequence of that decision is stated precisely rather than presented as an absence of leakage.** The candidate set is the product of two stages of different character. The genes *submitted* to Mendelian randomisation were the candidates of §2.5, formed by intersecting the disease-associated co-expression modules of §2.4 with the sex-stratified differentially expressed genes of §2.3; both of those analyses were computed on the whole training partition and both used the training samples' diagnosis labels. The *filter applied* to that set was estimated entirely outside this cohort, from external genome-wide association and eQTL summary statistics (Okada et al., 2014; Võsa et al., 2021), and could not have been influenced by any sample analysed here. The candidate set is therefore an externally filtered internal list rather than an external list, and it is not claimed to be free of selection bias on the strength of its genetic component alone.

Because re-deriving the network, the differential expression and the Mendelian randomisation inside every outer fold of every repeat was computationally out of reach, a sample in an outer test fold contributed, through the training-wide differential-expression and module–trait analyses, to defining the gene list from which its own prediction is subsequently built. The nested estimate accordingly corrects the selection bias arising in the three-selector panel-building stage and in the classifier, and does **not** correct selection bias arising in the upstream differential-expression, network and prioritisation stages; it is an upper bound on out-of-sample performance in that respect (Ambroise and McLachlan, 2002; Simon et al., 2003; Kaufman, Rosset, and Perlich, 2012).

Three considerations bound the residual optimism without eliminating it. The upstream stages are severe filters that are nonetheless indifferent to the classifier: the module–trait rule of §2.4 selects on the correlation between a module eigengene and diagnosis rather than on any individual gene's discriminative power, and the Mendelian randomisation filter of §2.6 is estimated entirely outside this cohort, so neither stage can select a gene *because* it separates these particular samples. The candidate set is small relative to the transcriptome, so the combinatorial freedom available to the fold-level selectors — the quantity that governs selection bias — is correspondingly small. And the two evaluations on which the diagnostic claim actually rests, the sealed internal holdout and the independent external cohorts, are unaffected by this argument entirely, those samples having taken no part in differential expression, network construction, prioritisation or selection at any stage.

Subsequently, the training model was applied to two distinct datasets: an external blood cohort (GSE15573) and a sealed 30% internal-validation hold-out. For the external cohort, probe intensities were converted to the log₂ scale if not already in that format, with a maximum observed intensity above 50 serving as the test criterion. These intensities were aggregated to genes using the MaxMean rule, retaining the probe with the highest mean expression for each symbol (Miller et al., 2011). The internal-validation hold-out was projected onto the training scale using frozen quantile-normalisation and ComBat parameters, ensuring that no hold-out sample influenced any training parameter (Bolstad et al., 2003; Johnson, Li, and Rabinovic, 2007; Leek et al., 2012; Nygaard, Rødland, and Hovig, 2016). The locked model was applied once to each dataset without refitting or tuning. Within each dataset, genes were standardised independently within the relevant sex stratum. This standardisation was necessary due to the lack of a common intensity scale across platforms and was legitimate because the transformation was unsupervised, relying solely on expression values without using diagnostic labels (Kaufman, Rosset, and Perlich, 2012). A panel gene absent from a target platform or exhibiting zero variance was assigned a z-score of zero, equivalent to the dataset mean, ensuring its contribution to the linear predictor was nullified rather than discarding affected samples. Genes missing from each dataset were recorded alongside the corresponding estimate. The external cohort, profiled on a different platform, represents a genuine test of transportability rather than an additional internal split (Justice, Covinsky, and Berlin, 1999; Steyerberg and Harrell, 2016).

No aspect of the evaluation was adjusted based on its outcome: the direction for comparing the receiver operating characteristic was predetermined for each dataset, rather than being selected to optimize the area under the curve (AUC). No model was refitted after examining a validation estimate, and no random seed was chosen based on the estimate it generated. An AUC of 0.999 or higher was flagged as indicating separation, as the AUC inherently represents the probability that a randomly selected case is ranked above a randomly selected control (Hanley and McNeil, 1982). At small sample sizes, a perfect AUC value reflects the limited number of case-control pairs available to be discordant, rather than indicating a perfect classifier. Each analysis was categorized into an evidence tier, either primary or exploratory, based solely on sample size and not on the magnitude of the estimate. The male stratum falls below the pre-specified thresholds, with its events per candidate predictor significantly below the conventional minimums for a stable multivariable model (Peduzzi et al., 1996; Riley et al., 2019), and is consistently reported as a power-limited exploratory analysis. Reporting adheres to the Transparent Reporting of a multivariable prediction model for Individual Prognosis Or Diagnosis (TRIPOD) statement (Collins et al., 2015).

**2.10 Diagnosis-by-sex interaction testing**

The stratified design of §2.3 to §2.9 cannot, by construction, establish that anything differs between the sexes: two lists thresholded separately in strata of unequal size will differ even when the underlying effects are identical. A formal diagnosis-by-sex interaction test was therefore fitted across the transcriptome, and it is the **only** analysis in this chapter capable of supporting a sex-differential claim. Its null hypothesis is that the RA-versus-healthy-control effect on a gene is the same in both sexes; a gene rejecting that null is termed **sex-differential**. A gene may be strongly associated with RA within a single sex without its effect differing between sexes, and such a gene is not sex-differential no matter how significant it is within its own stratum.

**This test is also the point of contact with the other omics layers of this thesis, and it is placed here for that reason.** The methylomic analysis of the companion chapter estimates a diagnosis-by-sex interaction directly, one coefficient per probe. The estimand tested in this section — the `groupRA:sexM` coefficient of a `~ group * sex + batch` model — is the transcriptomic counterpart of that quantity, and the two are therefore directly comparable at the level of the gene. The sex-stratified panels of §2.9 are **not** comparable to it: they are the output of a different design, they carry no cross-sex contrast, and integrating them with methylomic interaction hits as though the two were the same kind of object would compare a design with a finding. Any multi-omics integration on the sex axis must therefore be built on the interaction results of this section, with the stratified panels entering only as a separate, diagnostic-performance layer.

The primary test was fitted on the **pre-ComBat, quantile-normalised** training matrix with batch modelled explicitly, using a limma linear model of the form expression ~ group \* sex + batch_full, with healthy control and female as reference levels. The coefficient of interest was the group × sex interaction term (groupRA:sexM), which quantifies the difference between males and females in the RA-versus-healthy-control expression change. Gene-wise models were fitted using lmFit and moderated by empirical Bayes (eBayes), and the interaction coefficient was extracted for all genes using topTable.

**The choice of matrix is deliberate and is the reason the ComBat-corrected model is not the primary analysis.** ComBat was run with a covariate model matrix of the form ~ group + sex, which protects the two main effects but **not their interaction**. Any group-by-sex structure correlated with batch is therefore partly absorbed and partly redistributed by the correction, and a model fitted on the ComBat output is additionally never debited the degrees of freedom that ComBat itself consumed (Nygaard, Rødland and Hovig, 2016). Fitting the interaction on the ComBat matrix consequently returns a substantially larger and **inflated** gene count. Both matrices were fitted and the counts compared; the pre-ComBat model is the one reported, and the ComBat count is retained solely as a sensitivity figure. A third specification, adding leukocyte-composition principal components to the pre-ComBat design (§2.14), was fitted to establish how much of the interaction signal survives adjustment for blood composition.

Two limits on interpretation are recorded here rather than left to the Results. Interaction tests carry roughly a quarter of the power of main-effect tests at the same sample size, and the male stratum contributes the smaller number of cases, so the number of genes recovered is a **floor** rather than an estimate of how many sex-differential genes exist. And an interaction analysis conducted on the transcriptome does not retrospectively confer sex-specificity on genes selected by a different procedure: the diagnostic panels of §2.8 are sex-**stratified**, meaning fitted separately within each sex, and are described as such throughout.

Multiple-testing correction was implemented via the Benjamini-Hochberg procedure across the entire transcriptome, and deliberately not restricted to the panel genes, so that false-discovery-rate control reflects the full set of tests conducted rather than a subset chosen after the fact. Interaction statistics were then looked up for the consensus panel genes of §2.9, a gene being classified **sex-differential** where its interaction reached transcriptome-wide FDR below 0.05 and **sex-shared** otherwise.

The direction of this procedure must be stated carefully, because it is a lookup and not a second test. The panels were selected by a procedure that never examined a cross-sex contrast, so this step cannot confer a sex-differential property on them retrospectively; it can only report whether genes selected on other grounds happen also to reject the interaction null. A panel consisting entirely of sex-shared genes is a coherent and unsurprising result, and would mean that the two panels differ in membership because the stratified selection ran on different samples, not because the underlying biology differs by sex. That is the expected outcome of a stratified design and is reported as such rather than as a failure.

**2.11 Cross-tissue sex-stratified biomarker validation**

To assess the discriminatory efficacy of blood-derived sex-stratified biomarker panels within disease-target tissue, these panels were validated using an independent RNA-sequencing dataset of rheumatoid arthritis synovium (GEO accession GSE89408). This validation serves as a cross-tissue generalisation test, given that the panels were initially developed using peripheral blood microarray data, whereas synovium represents a distinct tissue evaluated through a different technological platform. The transferable element between these settings is the panel gene set; however, the coefficients of the blood-trained classifier were not applied to the synovial data due to the differing scales of the two platforms. Consequently, the multivariable model was re-estimated within the synovium, as detailed below. Successful performance of the gene set in this context would suggest a tissue-transferable disease signal rather than a blood-specific artefact. It is important to note that this does not constitute an out-of-sample validation of the blood model itself and is not reported as such (Harrell, Lee, and Mark, 1996; Steyerberg and Harrell, 2016).

Initially, raw synovial RNA-seq counts were loaded, with samples limited to the rheumatoid arthritis (RA) and normal (control) groups. Disease status was determined from the sample-name prefix of the count matrix. The sex of each sample was assigned based on the sex field in the series metadata and aligned with the count matrix by column position, ensuring complete alignment across all columns prior to analysis. The counts were processed using edgeR: a DGEList was constructed, and genes with low expression were excluded using filterByExpr, categorized by disease status (Robinson, McCarthy, and Smyth, 2010; Chen, Lun, and Smyth, 2016). Library sizes were normalized using the trimmed mean of M-values (TMM) method via calcNormFactors, and expression levels were represented as log₂ counts-per-million (cpm, log = TRUE, prior.count = 1) (Robinson and Oshlack, 2010). Differential expression between RA and normal samples was assessed using limma-voom, employing a design matrix of the form ~ sex + grp to adjust for the RA effect (coefficient grpRA) concerning sex. This design estimates a single RA coefficient applicable to both sexes, thus the differential-expression step is sex-adjusted rather than sex-stratified, whereas subsequent discrimination analyses are conducted within each sex. The mean-variance relationship was modeled using voom, gene-wise linear models were fitted using lmFit, and moderated with eBayes (Smyth, 2004; Law et al., 2014; Ritchie et al., 2015). Per-gene log₂ fold-changes and Benjamini-Hochberg-adjusted p-values were extracted using topTable without re-sorting, preserving the gene order of the count matrix (Benjamini and Hochberg, 1995). For each panel gene, the direction of its synovial RA-versus-normal change was compared with the direction of its blood training differential expression, specifically the sign of the training log₂ fold-change in the corresponding sex. A gene was flagged as concordant when the two directions matched, facilitating the explicit identification of genes reversing direction across tissues rather than concealing them.

Each gene within the panel was assessed as a univariate classifier to differentiate rheumatoid arthritis (RA) from normal synovium, using its log₂-CPM values. The area under the ROC curve (AUC) was computed both overall, encompassing pooled RA and normal samples from both sexes, and separately for each sex (Hanley and McNeil, 1982; Robin et al., 2011). In both the tabulated per-gene AUCs and the per-gene ROC overlay, each gene was oriented towards its optimal direction, ensuring that every reported AUC is at least 0.5 by design; thus, direction reversal between tissues is not visible in these values and is instead indicated by the concordance flag defined earlier. For the cross-tissue comparison figure, the same AUCs were re-expressed based on the blood training orientation, with a discordant gene reported as one minus its best-direction AUC, so that a reversal is indicated by an AUC below 0.5; these two conventions are labeled wherever they appear and should not be interpreted interchangeably. The entire panel was analyzed as a combined multivariable logistic regression classifier (glm, family = binomial), fitted on the z-scored log₂-CPM values of the panel genes within the synovial samples of the respective sex, with undefined z-scores set to zero. Since this model was fitted and scored on the same synovial samples, its AUC is an apparent, resubstitution estimate and is labeled accordingly; it is neither cross-validated nor bootstrap-corrected, and it lacks a confidence interval, with the ROC direction predetermined using normal as the reference level (Harrell, Lee, and Mark, 1996; Robin et al., 2011).

**2.12 Cross-ancestry sex-stratified biomarker evaluation**

To evaluate the consistency of genes identified as putatively causal in the European discovery analysis across different ancestral groups, the European MR-prioritised gene sets for each sex were re-assessed using two-sample Mendelian randomisation against rheumatoid arthritis (RA) GWAS from three cohorts representing two ancestries (Davey Smith and Hemani, 2014; Hemani et al., 2018). The Okada 2014 European GWAS (ieu-a-832) served as the discovery outcome, the Stahl 2010 European GWAS (ieu-a-834) was employed as an ancestry-matched replication, and the Biobank Japan East-Asian GWAS (bbj-a-151) functioned as the cross-ancestry test (Stahl et al., 2010; Okada et al., 2014; Ishigaki et al., 2020), all accessed through the IEU OpenGWAS infrastructure (Elsworth et al., 2020). Given that OpenGWAS provides only European cis-eQTL data (eQTLGen) and lacks an East-Asian cis-eQTL resource, the same European eQTLGen instruments were applied to all three outcomes (Võsa et al., 2021). Consequently, the two European arms are ancestry-matched, whereas the East-Asian arm is ancestry-mismatched on the exposure side and was considered exploratory cross-ancestry support rather than a direct replication, due to differences in allele frequencies and linkage-disequilibrium structure at the instrument loci between the two populations, and an instrument valid in Europeans may not tag the same regulatory variation in East Asians (Martin et al., 2019). A fixed seed (\`set.seed(2024)\`) was consistently used. It is important to note that all three outcome GWAS are sex-combined, so the Mendelian randomisation estimates themselves are not sex-specific; what is stratified by sex is the gene set under test, drawn from the female and male MR-prioritised sets, respectively. The estimation was thus performed once over the union of the two sets, and the estimates were then partitioned by sex for reporting and classification, ensuring that a gene belonging to both panels carries an identical estimate in each. This section therefore investigates whether the sex-stratified candidate gene sets are supported in other ancestries, and not whether the genetic effects themselves differ by sex, which these data cannot address.

The instruments used were the European eQTLGen cis-eQTL SNPs, which had been previously selected in the primary Mendelian randomisation (MR) analysis of §2.6 based on an instrument-strength F-statistic exceeding 10, a conventional threshold below which weak-instrument bias becomes significant (Burgess and Thompson, 2011). These instruments were reused from the stored primary MR objects and were restricted to the prioritised genes, ensuring consistency in the exposure side across all outcomes. For the two newly tested outcomes, the instrument SNPs were retrieved using the \`extract_outcome_data\` function and harmonised with the exposure data through the \`harmonise_data\` function, employing action = 2, which infers strand orientation from allele frequencies and excludes ambiguous palindromic SNPs (Hartwig et al., 2016; Hemani et al., 2018). The Okada results were reused from the primary pipeline without re-fetching, ensuring that the discovery arm is numerically identical to §2.6 rather than a re-estimation. MR was conducted for each gene using the same SNP-count-adaptive method ladder as in the primary analysis: inverse-variance weighted, MR-Egger, and weighted median methods were applied when three or more instruments were available; inverse-variance weighted was used for two instruments; and the Wald ratio was applied for a single instrument (Burgess, Butterworth and Thompson, 2013; Bowden, Davey Smith and Burgess, 2015; Bowden et al., 2016). A single primary estimate per gene was selected following a fixed hierarchy: inverse-variance weighted > Wald ratio > weighted median > MR-Egger. The estimates were expressed as odds ratios (OR = exp(β)) with 95% confidence intervals obtained on the log-odds scale by the normal approximation, exp(β ± 1.96 × SE). Each gene was then classified within each cohort by evaluating significance first: a gene with p ≥ 0.05 was labelled as having no significant effect, and only otherwise was it classified as conferring risk (OR > 1) or protection (OR < 1). These thresholds are nominal and uncorrected, in contrast to the Benjamini-Hochberg false discovery rate applied when the candidate set was first defined (Benjamini and Hochberg, 1995); the replication and cross-ancestry arms are therefore interpreted as directional corroboration of a pre-specified gene set rather than as independent hypothesis tests.

In the East-Asian cohort, a null result may indicate either a true absence of effect or the unavailability of European instruments. To address this, two diagnostics for instrument transferability were calculated alongside each estimate. The first diagnostic involved determining the number of instrument SNPs per gene that remained after harmonisation for each outcome. This allows for the classification of a gene as untestable in East Asians if no instrument survives, rather than concluding it is non-causal in this population. The second diagnostic was the mean absolute difference in effect-allele frequencies between exposure and outcome across a gene's retained instruments, serving as an index of allele-frequency divergence between the European exposure panel and the outcome cohort (Martin et al., 2019). Directional agreement between cohorts was evaluated by examining the sign of the odds ratio relative to the null, sign(OR - 1). Based on these criteria, each gene was categorized into one of four mutually exclusive classes, assessed in a fixed sequence: untestable in East Asians, where no instrument survived; shared European and East Asian, where the East-Asian estimate aligned directionally with the discovery estimate at p < 0.05; European-replicated but not East Asian, where the Stahl estimate aligned directionally at p < 0.05; and European-discovery only, where neither condition was satisfied. This classification reflects instrument availability and directional concordance, and a gene classified in the final category is not necessarily non-causal in other ancestries.

**2.13 Functional enrichment**

Functional interpretation was conducted at two distinct stages within the analytical pipeline. The first stage involved the analysis of differentially expressed genes (DEGs) following expression analysis, while the second stage focused on Mendelian-randomisation prioritised genes after the causal-inference phase. Both stages utilised the clusterProfiler framework (Yu et al., 2012; Wu et al., 2021), employing human annotation from org.Hs.eg.db (Carlson, 2019). Gene symbols were consistently mapped to Entrez identifiers using bitr; symbols lacking an Entrez mapping were excluded at this step, resulting in both the tested gene sets and the background universes being the mapped subsets rather than the nominal ones. Two additional over-representation analyses, concerning the WGCNA disease-module genes and the sex-by-diagnosis interaction genes, are detailed in separate sections and are not included in the two stages discussed here. Over-representation analysis and ranked gene-set enrichment address different questions: the former assesses whether a fixed gene list is enriched relative to a background, while the latter evaluates whether a pathway is systematically displaced in a ranking of all genes. Consequently, both methods were retained rather than considered alternatives (Goeman and Bühlmann, 2007).

During the differential-expression phase, significant DEG sets were determined for the pooled cohort and separately for females and males prior to conducting over-representation analysis. Gene Ontology enrichment was performed using enrichGO across the three ontologies: biological process, cellular component, and molecular function (Ashburner et al., 2000; The Gene Ontology Consortium, 2021). Concurrently, KEGG pathway enrichment was executed with enrichKEGG for Homo sapiens (organism = "hsa") (Kanehisa and Goto, 2000). Both analyses employed the Benjamini-Hochberg correction, applying a cutoff of 0.05 to both raw and adjusted p-values, and a q-value threshold of 0.2 (Benjamini and Hochberg, 1995; Storey and Tibshirani, 2003). Gene Ontology results were converted to readable gene symbols, while KEGG identifiers were maintained in Entrez form at this stage. As no background universe was provided, comparisons were made against the complete set of annotated human genes rather than the genes detectable on the arrays analyzed here. This approach is more permissive than the universe-restricted analysis described subsequently, as gene classes disproportionately represented among array-detectable transcripts may be enriched for that reason alone, rendering the two stages non-interchangeable (Timmons, Szkop, and Gallagher, 2015; Wijesooriya et al., 2022). Additionally, a ranked Gene Set Enrichment Analysis of KEGG pathways was conducted using gseKEGG (organism = "hsa", Benjamini-Hochberg correction, p-value cutoff 0.05, seed = TRUE to ensure reproducibility of the permutation procedure) (Subramanian et al., 2005). For each comparison, genes were ranked by their limma moderated t-statistic for the RA-versus-healthy-control contrast, retaining a single ranking value per gene, specifically the largest absolute t-statistic when a gene mapped to multiple identifiers (Smyth, 2004; Ritchie et al., 2015). A positive normalised enrichment score indicated a pathway up-regulated in RA, whereas a negative score indicated down-regulation in RA.

During the causal-gene phase, Mendelian randomisation was employed to prioritise gene sets for each sex, which were then annotated separately to clarify the biological characteristics of genes with genetic support for a causal role in rheumatoid arthritis (RA). The primary candidate sets, which retained the major histocompatibility complex (MHC), were utilised, indicating that genes from the MHC contribute to any enriched immune terms. The MHC-free sensitivity sets, defined elsewhere in this chapter, were not annotated separately, and thus, no enrichment result reported here can be considered independent of the HLA region. The background universe was limited to genes expressed in the training cohort, specifically the rows of the training expression matrix mapped to Entrez identifiers, rather than the entire genome. This ensured that enrichment was assessed against genes measurable in the study (Timmons, Szkop, and Gallagher, 2015). Gene Ontology biological-process enrichment was conducted using enrichGO (ont = "BP") with this expressed-gene universe, applying a p-value cutoff of 0.05, a q-value cutoff of 0.2, and readable output. KEGG pathway enrichment was performed with enrichKEGG (organism = "hsa"), using the same universe and a p-value cutoff of 0.05, with the default q-value cutoff of 0.2, both under the default Benjamini-Hochberg adjustment, with results mapped back to gene symbols. The most significant enriched terms for each sex were visualised as gene-ratio dotplots, which combined the top ten Gene Ontology biological-process terms and the top ten KEGG terms ranked by adjusted p-value. These were faceted by source, with the gene ratio on the horizontal axis, point size proportional to the number of genes in the term, and point colour encoding the adjusted p-value.

**2.14 Immune deconvolution and composition-adjusted expression**

To characterize the immune-cell composition underlying the transcriptomic signal, whole-blood expression matrices were deconvolved using CIBERSORT with the LM22 signature matrix, which resolves 22 human immune-cell subsets, as implemented in the IOBR package (deconvo_tme, method = "cibersort") (Newman et al., 2015; Zeng et al., 2021). LM22, derived from microarray data, serves as the appropriate reference for the Affymetrix cohorts analyzed here. Deconvolution was applied independently to the training cohort, the sealed internal-validation holdout, and the external blood cohort, deliberately avoiding joint fitting across them. This approach ensures that the holdout remains a sealed test set, preventing its samples from influencing a quantity subsequently used to adjust the training model. Since CIBERSORT requires a non-logarithmic mixture, the ComBat batch-corrected log₂ expression matrix was first converted back to the linear scale (2^expression) prior to deconvolution, with any non-finite value set to zero. The analysis was conducted in microarray mode (arrays = TRUE), applying the quantile-normalized handling of the mixture appropriate for array data (Bolstad et al., 2003). The goodness-of-fit for deconvolution was evaluated through 100 permutations (perm = 100), providing a per-sample empirical p-value along with the CIBERSORT correlation and root-mean-square-error diagnostics. A fixed seed (set.seed(1234)) was established for reproducibility. These three diagnostics were retained solely for quality assessment and were explicitly excluded from all downstream modeling, where only the 22 LM22 fraction columns were utilized. As a methodologically independent check, MCP-counter was additionally run on the training and external cohorts (Becht et al., 2016). Since it returns relative abundance scores rather than fractions, it was used solely to corroborate the direction of compositional differences and never as the basis of any adjustment, based on the principle that agreement between two estimators that fail in different ways is more informative than either alone (Shen-Orr and Gaujoux, 2013). The output consisted of a table of the 22 estimated cell-type fractions for each sample in each dataset, which was integrated with the sample-level sex and disease-group annotations for subsequent analysis.

To investigate the association between immune composition and disease signature, the estimated cell-type fractions were compared between rheumatoid arthritis (RA) cases and healthy controls, with analyses conducted separately for each sex to avoid pooling female and male samples. For each of the 22 cell types, the difference in fractions between RA and control samples was evaluated using the Wilcoxon rank-sum test (Wilcoxon, 1945; Mann and Whitney, 1947). This non-parametric test was chosen due to the bounded nature of cell fractions, their frequent zero inflation, and their deviation from a Gaussian distribution. Cell types exhibiting zero variance within a stratum were excluded from testing. The resulting p-values were adjusted for multiple comparisons across the 22 subsets using the Benjamini-Hochberg procedure, applied separately within each sex. A subset was identified as differing by disease status at a false discovery rate below 0.05 (Benjamini and Hochberg, 1995). This analysis determines whether composition acts as a confounder or mediator of the differential-expression signal. If the fractions do not differ, composition cannot drive the signal, rendering any adjustment a formality. Conversely, a difference necessitates reporting unadjusted expression results as confounded with or mediated by composition. Given that fractions are compositional and sum to one, they were subsequently summarized in centered log-ratio space using principal components for inclusion as covariates in the composition-adjusted analyses, detailed in a separate section (Aitchison, 1982).

**2.15 Nomogram construction and clinical evaluation**

A nomogram is a graphical tool that represents a fitted regression model by mapping each predictor onto a points scale. This allows for the manual evaluation of the linear predictor and the predicted probability for an individual patient without the need for software (Harrell, 2015; Iasonos et al., 2008). The nomogram does not add new information to the model from which it is derived; its function is to render the model transparent and computable at the point of use.

Calibration refers to the concordance between the probabilities predicted by a model and the actual frequencies of the outcome. This property is distinct from discrimination, which pertains to a model's ability to rank cases above controls effectively, thereby achieving a high area under the receiver operating characteristic curve, even if it systematically over- or underestimates absolute risk (Van Calster et al., 2019; Alba et al., 2017). Since a curve fitted and evaluated on the same samples is optimistically close to the diagonal by design, calibration is reported as both the apparent curve and a bias-corrected curve. The latter is obtained through resampling to estimate and eliminate this optimism (Harrell, Lee, and Mark, 1996; Steyerberg et al., 2001).

Decision-curve analysis assesses whether acting on a model yields more benefit than harm by placing true and false positives on a common scale. The central metric, net benefit, is calculated as the proportion of true positives minus the proportion of false positives, weighted by the odds of the threshold probability \\( p_t \\). Here, \\( p_t \\) represents the predicted risk at which a clinician would decide to act, and its odds, \\( p_t/(1 - p_t) \\), reflect the relative harm of an unnecessary intervention versus a missed case (Vickers and Elkin, 2006; Vickers, van Calster, and Steyerberg, 2016). A model is deemed useful at a given threshold only if its net benefit surpasses that of both default strategies: treating every patient and treating none.

A clinical-impact curve translates the model into units that quantify consequences, indicating, for a hypothetical population at each threshold, how many individuals would be classified as high risk and how many of those would genuinely have the disease (Kerr et al., 2016).

All four analyses were conducted separately for each sex, as detailed below.

The procedure was conducted separately for each sex. Distinct instruments were developed for females and males, each based on a model tailored to the training samples of the respective sex, utilizing the consensus panel genes specific to that sex. Consequently, no sample, coefficient, predictor axis, or scaling constant was shared between the sexes (Harrell, 2015). Sex was not included as a predictor; instead, it functioned as the partition within which all subsequent computations were executed, resulting in instruments that lacked a sex axis and were not applicable across sexes.

2.15.3 Model fitting

For each sex, a binary logistic regression model was constructed using rms::lrm, with rheumatoid arthritis (RA) status as the dependent variable and the panel genes as predictors. The predictors were incorporated as additive linear terms on the logit scale, without employing spline expansions, interaction terms, or categorisation of expression. This approach assumes a linear relationship between each gene's log₂ expression and the log-odds of RA across the observed range of that gene. The genes were included as unstandardised log₂ expression values, as presented in the processed training matrix, rather than as z-scores. A datadist object summarising the distribution of each predictor within that sex was computed and registered prior to model fitting, as rms requires it to determine the plotted range of each predictor axis. For the male panel, a mild ridge penalty (penalty = 5) was applied to stabilise a fit that nearly completely separates the training data at that stratum's sample size, whereas the female model was fitted without penalty (Hoerl and Kennard, 1970; Albert and Anderson, 1984; Harrell, 2015). The design and response matrices were retained (x = TRUE, y = TRUE) to facilitate the application of resampling procedures to the fitted object.

2.15.4 Nomogram construction

The nomogram was derived from the fitted model using rms::nomogram. Each panel gene was represented by its own axis, which spanned the range of that gene's observed expression within the corresponding sex as recorded in the datadist object, ensuring that the axes were empirical and stratum-specific rather than nominal. Points were allocated according to the standard nomogram scaling, where the predictor with the widest contribution to the linear predictor was assigned a range of 0 to 100 points, and all other predictors were scaled proportionally to the width of their contributions. Consequently, the relative length of the gene axes reflects the product of each coefficient and the spread of that gene's expression in that sex, rather than the coefficient alone, and does not represent a ranking of effect size (Balachandran et al., 2015). The instrument was interpreted by identifying a patient's value on each gene axis, noting the points indicated above it, summing those points across the panel genes, and projecting the total onto a predicted-probability axis. This axis was derived by applying the logistic transformation to the corresponding value of the linear predictor (fun = plogis), labelled "Risk of RA", and annotated at predicted probabilities of 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, and 0.99. Points were thus expressed per unit of log₂ expression, and a reading is valid only for expression values placed on the same processed scale as the training matrix from which the model was derived.

2.15.5 Calibration

Calibration was evaluated separately for each sex by plotting observed probabilities against predicted probabilities. This involved reporting both the apparent curve and a bootstrap bias-corrected curve derived from 200 resamples (rms::calibrate), allowing for the estimation of optimism induced by fitting and evaluating on the same samples, rather than assuming its absence (Harrell, Lee, and Mark, 1996; Steyerberg et al., 2001; Van Calster et al., 2019). In instances where the resampling procedure did not yield a calibration object, the corresponding panel was left blank instead of being replaced by an uncorrected curve.

2.15.6 Decision-curve analysis

For each sex, net benefit was calculated using the formula:

NB = TP/n − FP/n × p_t/(1 − p_t)

This calculation was performed across threshold probabilities ranging from 0.01 to 0.99 in increments of 0.01. The net benefit was then compared to two default strategies: treating all patients, where net benefit is determined solely by the observed prevalence in that sex stratum, and treating none, where net benefit is zero by definition (Vickers and Elkin, 2006). Given that the treat-all reference is dependent on prevalence, the decision curves for females and males are interpretable only within their respective strata and cannot be overlaid (Vickers, van Calster, and Steyerberg, 2016). An additional cost-to-benefit ratio axis, p_t/(1 − p_t), was included for interpretation, annotated at ratios of 1:100, 1:4, 2:3, 3:2, 4:1, and 100:1. The probabilities used in this curve were the apparent, resubstitution predictions of the fitted model, and thus the curve was not corrected for optimism.

2.15.7 Clinical-impact analysis

For each sex, the number of individuals per 1,000 classified as high risk, along with the number of those who genuinely have RA, was plotted across the same range of thresholds (Kerr et al., 2016). Uncertainty was quantified using 500 bootstrap resamples drawn with replacement from the training samples of each sex, each seeded by its own replicate index for reproducibility. The model was refitted in each resample under the same penalty as the primary model, with the 2.5th and 97.5th percentiles across replicates providing 95% percentile bands for both curves (Carpenter and Bithell, 2000). Within each replicate, the refitted model was evaluated on that replicate's own resampled observations, so the bands describe the sampling variability of the apparent curves and do not correct them for optimism.

**3. Results**

**Sex-stratified transcriptomic signatures for RA**

Results are reported in the same order as the methodology of §2, so that every
number below can be traced to the analytical step that produced it. Every
figure quoted in this section is machine-derived from the corresponding
output table in `results/tables/` and is reproducible by re-running the
script named in the corresponding methods subsection; none is transcribed by
hand from an intermediate note.

**3.1 Whole-blood cohort description**

After probe-to-gene collapse, GSE93272 contributed 20,848 genes and
GSE110169 contributed 19,041 genes; the intersection of the two platforms'
feature spaces, carried forward through the remainder of the chapter, was
**15,763 genes**. The stratified 70:30 partition assigned **183 samples** to
training — 103 RA and 80 healthy control, of whom 145 were female (86 RA, 59
HC) and 38 were male (17 RA, 21 HC) — and **74 samples** to the sealed
internal hold-out, comprising 36 female RA, 25 female HC, 6 male RA and 7
male HC. The male stratum is therefore a fixed constraint on every result
reported below, at roughly a fifth of the training cohort; this reflects the
sex composition of the two source series, both recruited from predominantly
female RA populations, rather than a property of the sampling procedure.

**3.2 Normalisation and batch correction**

Before quantile normalisation, the merged training matrix showed a standard
deviation of per-sample medians of 0.574 and of per-sample interquartile
ranges of 0.799, indicating distributional heterogeneity between the two
source series sufficient to require correction. Quantile normalisation
reduced these to 3.3 × 10⁻⁴ and 4.6 × 10⁻⁴ respectively — a reduction of more
than three orders of magnitude — confirming that the merged distributions
were aligned. ComBat, applied at the study-and-internal-batch level with
diagnosis and sex protected in the model matrix, resolved **six batches**.

**3.3 Whole-blood differential gene expression**

All three contrasts (pooled, female-only, male-only) were fitted on the 183
training samples at |log₂FC| > 0.1 and FDR < 0.05. The pooled contrast
returned **6,422 significant genes** (2,760 up in RA, 3,662 down; 40.7% of
the tested transcriptome).

*Female RA versus female control DEGs.* The female contrast identified
**5,131 DEGs** (2,238 up in RA, 2,893 down; 32.6% of the tested
transcriptome). The most significant genes were *HMGB2* (log₂FC +0.67, FDR
2.7 × 10⁻¹²), *MAGED1* (−0.39), *KDM1A* (−0.29), *C1GALT1C1* (+0.65) and
*S100A8* (+0.50, FDR 5.1 × 10⁻¹¹); the largest effects among significant
genes were *ARG1* (+1.36), *CLEC4D* (+1.17), *BCL2A1* (+1.16), *COX7B*
(+1.02) and *DEFA4* (+1.00), with a median |log₂FC| among significant genes
of 0.192 and a minimum attainable FDR of 2.7 × 10⁻¹².

*Male RA versus male control DEGs.* The male contrast identified **5,820
DEGs** (2,510 up, 3,310 down; 36.9%). Leading genes by significance were
*STARD3NL* (+0.70), *TBC1D15* (+1.01), *MIER1* (+0.65), *TTC33* (+1.30) and
*TRIM23* (+1.21), all around FDR 3.3–3.4 × 10⁻⁸; the largest effects were
*BCL2A1* (+1.84), *COX7B* (+1.70), *RPL22L1* (+1.62), *COMMD6* (+1.60) and
*SCOC* (+1.58), with a median |log₂FC| of 0.337. **This larger male DEG count
is a statement about estimator variance, not about a stronger male signal**:
the median male effect size is 1.8-fold the female value and the minimum
attainable male FDR (3.3 × 10⁻⁸) is four orders of magnitude weaker than the
female minimum, exactly the pattern produced by noisier per-gene estimates
at n = 38 checked against a fixed effect-size gate rather than by a larger
true effect.

*Shared female and male blood DEGs.* Comparing the two significant-gene
lists, **3,857 genes are significant in both sexes**, 1,274 in women only and
1,963 in men only (union 7,094). Among the shared genes, concordance in
direction is near-complete: **3,854 of 3,857 (99.9%)** move the same way in
both sexes.

*Genes reaching significance in one stratum only.* A gene appearing on one
sex's list and not the other's is, by the design stated in §2.3, not thereby
shown to differ between the sexes: the female and male strata are
thresholded separately at markedly different sample sizes, and a set
difference between two independently thresholded lists conflates unequal
power with unequal effect. No gene is described as sex-specific on this
basis anywhere in this chapter. The only test that speaks directly to
whether the RA effect on a gene differs by sex is the diagnosis-by-sex
interaction test, reported in §3.10.

**3.4 Co-expression network analysis (WGCNA)**

Average-linkage clustering at a data-driven cut height of **69.32** (mean + 3
SD of merge heights) identified **10 outlying arrays**, which were removed,
leaving a network matrix of **173 samples × 15,763 genes**. `pickSoftThreshold`
returned an estimate of β = 3 (signed R² = 0.873), rejected on connectivity
grounds — mean connectivity 2,240.3 is far too dense for meaningful module
resolution. The signed-R² curve plateaus from approximately β = 9 onward
while mean connectivity continues to fall (R² = 0.897, mean k = 125.0 at
β = 9; R² = 0.908, mean k = 28.2 at β = 14; R² = 0.914, mean k = 8.5 at
β = 20). **β = 12 was used** (R² = 0.901, mean connectivity 47.4).

Blockwise construction at β = 12 detected **12 modules**, with **6,019 genes**
left unassigned (grey); the largest were turquoise (2,584 genes), blue
(2,554), brown (1,878) and yellow (708). Correlating each module eigengene
with RA status across all 173 samples, two modules met the pre-specified
selection rule (|cor| ≥ 0.5, p < 1 × 10⁻⁸): the **yellow module** (708 genes,
r = +0.556, p = 1.91 × 10⁻¹⁵, up in RA) and the **brown module** (1,878 genes,
r = −0.590, p = 1.29 × 10⁻¹⁷, down in RA), giving a pooled disease-module
background of **2,586 genes**. Neither disease module is a sex module: their
correlation with the male indicator is −0.176 (yellow) and +0.190 (brown),
an order of magnitude weaker than their disease correlations. Recomputed
within each sex, both modules retain their sign and significance — yellow is
+0.485 (p = 1.7 × 10⁻⁹) in women and +0.772 (p = 5.5 × 10⁻⁸) in men; brown is
−0.515 (p = 1.0 × 10⁻¹⁰) in women and −0.813 (p = 2.9 × 10⁻⁹) in men — though
the female and male coefficients are not compared directly against one
another, for the reasons given in §2.4.

Applying the hub definition (|kME| > 0.8 and |GS| > 0.2) identified **217 hub
genes in the yellow module** and **144 in the brown module**. The
highest-connectivity yellow hubs are *ZNF267*, *ZFYVE16*, *SP3*, *FAR1* and
*SNX13* (kME 0.94–0.95, connectivity 111–118); the highest-connectivity brown
hubs are *KHSRP*, *GPI*, *SF3B3*, *CLSTN1* and *XRCC6* (kME 0.88–0.91,
connectivity 118–130).

Preservation of the combined-network modules from a female reference to a
male test set (200 permutations) shows the random-gene benchmark module
(`gold`) already highly preserved at Zsummary = 34.19. Against that
reference, the yellow disease module is clearly preserved (Zsummary 48.14,
+13.95 over gold), while the **brown disease module exceeds the random
benchmark only marginally** (Zsummary 36.95, +2.76 over gold) — the
conventional rule of thumb that Zsummary > 10 denotes strong preservation
overstates the evidence for brown once it is compared against the correct
random-gene reference rather than a fixed cut-off.

**3.5 Candidate gene identification**

Intersecting the 2,586-gene disease-module background with each sex's
significant-DEG list gave **2,045 female candidates** (595 from yellow, 1,450
from brown) and **2,079 male candidates** (585 from yellow, 1,494 from
brown). Because no variance filter was applied before network construction,
**zero DEGs were lost to filtering in either sex**, and every candidate's
direction of differential expression agrees with its module's direction of
association with RA in both sexes. Of the two sex-specific candidate lists,
**1,773 candidates are shared**, 272 are female-list-only and 306 are
male-list-only, giving a **union of 2,351 genes** carried forward into the
Mendelian randomisation screen.

**3.6 Instrument selection and Mendelian randomisation results**

An MR estimate was obtained for **1,477 of the 2,351 union candidates in the
female stratum and 1,478 in the male stratum**. Within-stratum
Benjamini–Hochberg correction identified genes provisionally described as
MR-prioritised: **32 in the female stratum and 25 in the male stratum**, of
which 22/32 (female) and 19/25 (male) rest on a single instrument SNP. The
strongest associations by p-value were *HLA-DRB1* (OR 0.337, FDR
1.5 × 10⁻²⁴⁷, 1 SNP, Wald ratio), *WDR46* (OR 4.531, FDR 3.8 × 10⁻³⁹, 1 SNP),
*AIF1* (OR 0.556, FDR 1.2 × 10⁻²¹, 2 SNP, IVW), *VPS52* (OR 0.414, FDR
1.5 × 10⁻²¹, 1 SNP) and *GNL1* (OR 4.629, FDR 2.7 × 10⁻¹⁵, 1 SNP) in women,
with the same leading MHC signal reproduced in men alongside *AP4B1* (OR
1.424, FDR 3.2 × 10⁻¹⁰, 1 SNP). Comparing the two stratum-level
MR-prioritised lists, **24 genes are shared**, 8 are female-only (*FCGR2B*,
*GNL1*, *HLA-DPA1*, *IKZF4*, *LAX1*, *LSM2*, *RETSAT*, *ZBTB9*) and 1 is
male-only (*TAB1*). As stated in §2.6, because both the eQTLGen exposure and
the Okada outcome are sex-combined resources, the MR estimate for a given
gene is numerically identical regardless of which stratum's table it is read
from; sex-specificity at this stage arises entirely from which genes each
stratum's candidate list submitted to the screen.

**3.7 Sensitivity analyses: MHC exclusion and colocalisation**

Repeating the MR screen with every instrument inside the extended MHC
removed reclassified the female MR-prioritised set as **14 robust to MHC
exclusion, 14 untestable without an MHC instrument, 0 MHC-dependent and 4
surviving only through FDR re-ranking**; the male picture is nearly
identical (**14 robust, 10 untestable, 0 MHC-dependent, 1 FDR-rank-only**).
Restricted to the six primary panel genes reported in §3.8 for each sex, the
same pattern holds in both sexes: **3 of 6 robust, 2 of 6 untestable, 1 of 6
FDR-rank-only** — three of the six genes in each panel sit inside the MHC and
cannot be validated by cis-MR at all.

Bayesian colocalisation (`coloc.abf`) was run for every testable
MR-prioritised gene. Across the full set of **33 causal-screen genes**: 0
robustly colocalised, 2 prior-fragile, 4 suggestive, 9 with distinct causal
variants (all non-MHC and therefore interpretable), 13 inside the MHC where
the method's single-causal-variant assumption is known to be violated, and 5
inconclusive. At the level of individual panel genes, *IKZF3* reaches
PP.H4 = 0.774 at the default prior but collapses to 0.255 at the
conservative prior — the strongest panel gene, and still short of the 0.8
threshold under either prior; *INPP5B* (PP.H3 0.929) and *ESYT1* (PP.H3
0.912) show distinct causal variants, positive evidence that their cis-MR
estimates are LD-confounded rather than causal; and the four MHC panel genes
— *GNL1*, *C6orf136*, *VPS52*, *HLA-DMA* — each return PP.H3 ≈ 1.000 under an
assumption known to be violated in that region, and are reported as
MHC-unreliable rather than as evidence of distinct variants. **No panel gene
colocalises with the RA association under either prior.**

`coloc.susie`, which permits multiple causal variants and is therefore the
appropriate tool inside the MHC, was attempted on all **14 MHC genes and
resolved 3**. *WDR46* (1 eQTL credible set, 9 RA credible sets) reaches best
PP.H4 = 0.918 and is the only MHC gene with positive colocalisation support
under the multi-variant model; *GNL1* (1 eQTL set, 10 RA sets) returns best
PP.H4 = 0.000 and is now demonstrably driven by distinct variants rather
than merely unresolved; *HLA-DMA* remains inconclusive (PP.H4 = 0.323).
Across the MHC regions tested, SuSiE consistently found **9–10 independent
RA credible sets**, direct confirmation that `coloc.abf`'s
single-causal-variant assumption was badly violated there. Consistent with
the criterion set out in §2.7, no gene in either final panel satisfies both
conditions required for a causal claim, and every occurrence of "causal" in
earlier drafts of this chapter has accordingly been replaced by
"MR-prioritised".

**3.8 Feature selection**

In the female stratum the three selectors returned 7 (LASSO), 10 (RF) and 8
(SVM-RFE) genes from the 32 MR-prioritised inputs; their three-way consensus
was **6 genes: *C6orf136*, *ESYT1*, *GNL1*, *IKZF3*, *MED1*, *SMARCC2***. In
the male stratum the selectors returned 10 (LASSO), 10 (RF) and 11 (SVM-RFE)
genes from 25 inputs; consensus was **6 genes: *ESYT1*, *HLA-DMA*, *INPP5B*,
*MED1*, *SMARCC2*, *VPS52***. Repeating the identical procedure on the
MHC-free MR-prioritised input gave a **4-gene MHC-free female panel** (*CDC37*,
*IKZF3*, *MED1*, *SMARCC2*) and a **5-gene MHC-free male panel** (*MED1*,
*NCOA5*, *PHF19*, *SMARCC2*, *TAB1*). Comparing panel membership directly:
in the female arm, *C6orf136*, *ESYT1* and *GNL1* are dropped and *CDC37* is
newly added, while *IKZF3*, *MED1* and *SMARCC2* are retained; in the male
arm, *ESYT1*, *HLA-DMA*, *INPP5B* and *VPS52* are dropped and *NCOA5*,
*PHF19* and *TAB1* are newly added, while *MED1* and *SMARCC2* are retained.
*MED1* and *SMARCC2* are therefore the two genes retained in both sexes
under both the primary and MHC-free procedures, and are the most
selection-robust members of either panel.

**3.9 Diagnostic model development and validation**

*Model tuning and cross-validation performance.* Under the reconciled nested
cross-validation implementation, the **female, primary-panel, consensus**
configuration gives nested AUC **0.816 (95% CI 0.747–0.884)**, per-repeat SD
0.010, median 5 genes used, n = 145; the MHC-free consensus configuration is
statistically indistinguishable (AUC 0.798, 95% CI 0.725–0.871, median 3
genes) and is the configuration recommended in the authoritative table, on
the strength of equivalent performance from fewer, MHC-independent genes.
For the **male, primary-panel, consensus** configuration — reported as
exploratory throughout, at n = 38 — nested AUC is 0.896 (95% CI
0.800–0.993); the MHC-free consensus configuration gives 0.924 (95% CI
0.841–1.000) and is likewise the recommended male configuration. A DeLong
comparison of the primary versus MHC-free nested AUC found no significant
difference in either sex (female 0.816 vs 0.798, p = 0.322; male 0.896 vs
0.924, p = 0.582): **removing the MHC costs no detectable performance**, and
in the male arm the point estimate rises. Decomposing apparent, flat-CV and
nested-CV AUC quantifies the optimism of the panel-selection procedure
directly: female apparent 0.869 vs nested 0.816 (optimism +0.054); male
apparent 1.000 vs nested 0.896 (optimism +0.104).

*ROC and locked-transfer performance, including the independent external
blood cohort.* Transferring the locked, frozen female primary panel — its
coefficients unchanged — to unseen data gives apparent-train AUC 0.869
(0.813–0.925, n = 145), internal-hold-out AUC 0.823 (0.719–0.928, n = 61) and
external-blood AUC 0.957 (0.886–1.000, n = 24, GSE15573). The locked female
MHC-free panel gives apparent-train AUC 0.845 (0.782–0.908), internal AUC
0.806 (0.692–0.919) and external-blood AUC 0.743 (0.535–0.951) — lower on
external blood than the primary panel, but without the primary panel's
perfect-separation artefact. The locked male primary panel shows perfect
separation on training and internal-hold-out data (AUC 1.000 in both,
flagged `SEPARATION` and not read as a performance estimate, per the
criterion of §2.9) and AUC 0.750 (0.399–1.000, n = 9) on external blood; the
locked male MHC-free panel shows perfect train separation, internal AUC
0.917 (0.750–1.000, n = 13) and external-blood AUC 0.700 (0.400–1.000, n =
9).

Benchmarked under within-dataset resampling against a composition-only model
fitted from leukocyte fractions alone (no gene expression; see §3.14), the
**female primary panel** adds signal on training data (panel 0.831 vs
composition 0.662, Δ +0.169, p = 1.2 × 10⁻¹⁰) and on external blood (1.000 vs
0.429, Δ +0.571, p = 2.4 × 10⁻⁵, flagged separation), but is only marginally
better than composition on the internal test (0.721 vs 0.714, Δ +0.007, p =
0.0014). The **female MHC-free panel** beats composition on all three
datasets — training 0.812 vs 0.662 (Δ +0.150), internal test 0.789 vs 0.714
(Δ +0.074, p = 0.032) and external blood 0.843 vs 0.429 (Δ +0.414, no
separation flag) — making it the only female configuration to clear the
composition benchmark on every dataset without a separation artefact. The
**male primary panel** beats composition on training (0.933 vs 0.857, Δ
+0.076) and external blood (0.725 vs 0.575, Δ +0.150, not significant) and
matches it exactly on the degenerate, perfectly-separated internal test
(1.000 vs 0.798, Δ +0.202); the **male MHC-free panel** shows a smaller
training margin (0.888 vs 0.857, Δ +0.031) and a smaller external margin
(0.700 vs 0.575, Δ +0.125). Every male comparison against composition is
flagged `EXPLORATORY`, and the internal and external comparisons in both
panels carry a `SEPARATION` warning that limits how much weight the point
estimates can bear.

**3.10 Diagnosis-by-sex interaction testing**

The diagnosis-by-sex interaction model (`~ group * sex + batch_full`, fitted
on the pre-ComBat quantile-normalised matrix as specified in §2.10) identifies
**53 genes at FDR < 0.05**; the ComBat-matrix sensitivity model, whose
covariate specification does not protect the interaction term, returns an
inflated 270 genes and is not the reported figure. Of the 53, the pattern
breakdown is **35 male-restricted**, **12 magnitude-difference**, **4
opposite-direction**, **1 female-restricted** and **1** significant only as
an interaction with neither within-sex effect individually significant. The
leading genes by interaction FDR are *ZNF800*, *BCLAF1*, *MAP4K5*, *UBXN4*
and *RIF1* (interaction FDR 0.0024–0.0056, all male-restricted). *PTPN22*,
an established RA susceptibility gene, is up-regulated in both sexes but
roughly five-fold more strongly in men (female +0.171, FDR 0.046; male
+0.893, FDR 1.2 × 10⁻⁶; interaction FDR 0.038). Among the 53, only one gene
is sex-chromosome-linked (*KDM5D*), and **no gene from either diagnostic
panel appears among the 53** — consistent with the panels being
sex-stratified rather than sex-differential, as stated in §2.10. Refitting
the interaction model with leukocyte-composition principal components added
leaves only **7 of the 53** genes significant (*ZNF800*, *BCLAF1*, *MAP4K5*,
*UBXN4*, *RIF1*, *RAB2A*, *KDM5D*), indicating that most of the apparent
sex-differential signal at this stage is attributable to differences in
blood composition rather than to cell-intrinsic regulation. Given that
interaction tests carry a fraction of the power of the corresponding
main-effect tests and the male stratum contributes only 17 RA cases, the
53-gene set is read as a floor on the number of sex-differential genes
rather than a stable enumeration.

**3.11 Cross-tissue evaluation (synovium)**

Of the female panel's six genes, four retained expression and were tested in
synovial tissue (GSE89408): *C6orf136*, *GNL1* and *SMARCC2* replicate
direction from blood to synovium (synovial FDR 2.1 × 10⁻²¹ to 1.1 × 10⁻³²,
single-gene AUC 0.906–0.996), while *ESYT1* is directionally concordant but
weakly discriminating (AUC 0.511), and *IKZF3* and *MED1* reverse direction
between tissues. Overall, **4 of 6 female panel genes are directionally
concordant** between blood and synovium. In the male panel, *ESYT1*,
*INPP5B*, *SMARCC2* and *VPS52* replicate direction (AUC 0.575–0.983) while
*HLA-DMA* and *MED1* reverse; **4 of 6 male panel genes are directionally
concordant**. Refitting each sex's panel from scratch within the synovium
dataset, rather than transferring the blood-trained coefficients, gives
apparent AUC 1.000 and 10-fold CV AUC 0.986 (0.969–1.000) for the female
refit (n = 120: 106 RA, 14 normal), and apparent AUC 0.994 with 10-fold CV
AUC 0.776 (0.573–0.978) for the male refit (n = 60: 46 RA, 14 normal); as
elsewhere, the male synovium result is read as exploratory given the small
comparator group.

**3.12 Cross-ancestry evaluation**

Of the **32 female MR-prioritised genes tested** against the East Asian
(BioBank Japan) and ancestry-matched European replication GWAS, 23 replicate
in the European stratum, 7 are transferable to the East Asian dataset, 2 are
untestable in the East Asian data for lack of a surviving instrument, and 5
are European-only. Of the **25 male MR-prioritised genes tested**, 18
replicate in European data, 6 are East-Asian-transferable, 2 are
East-Asian-untestable and 3 are European-only.

**3.13 Functional enrichment**

*Differentially expressed genes.* Over-representation analysis of the
sex-stratified DEG lists returned **1,408 enriched terms for the female
DEGs**, **1,513 for the male DEGs** and **1,659 for the pooled DEGs**. Ranked
KEGG gene-set enrichment gave a strongly directional result in every
contrast: **19 significant pathways pooled (18 suppressed, 1 activated)**,
**17 in the female analysis (16 suppressed, 1 activated)** and **22 in the
male analysis (21 suppressed, 1 activated)**. The single activated pathway
in every contrast was **ribosome** (NES +1.82 to +1.89); the most suppressed
pathways included ATP-dependent chromatin remodelling and primary
immunodeficiency (pooled and female) and carbon metabolism and the citrate
cycle (male).

*Disease-module genes.* Over-representation analysis of the 2,586
disease-module genes (yellow + brown) against the 15,763-gene network
background returned **218 GO terms** and **16 KEGG pathways** at FDR < 0.05,
led by mRNA processing, ribonucleoprotein complex biogenesis and RNA
splicing (GO) and spliceosome and the citrate/TCA cycle (KEGG). Spliceosome,
ATP-dependent chromatin remodelling and carbon metabolism/TCA cycle are
therefore suppressed at the transcriptome-wide GSEA level *and* enriched
among the down-in-RA brown module genes — two methodologically independent
analyses converging on the same suppressed biology.

*MR-prioritised genes.* Annotating the 32 female and 25 male MR-prioritised
genes separately, GO:BP over-representation returned **46 significant terms
in the female set** and **45 in the male set**, and KEGG returned **25
pathways in the female set** and **14 in the male set**. The female
MR-prioritised gene set is dominated by antigen-presentation biology, headed
by antigen processing and presentation of exogenous peptide antigen via MHC
class II (4/31 genes, FDR 1.2 × 10⁻⁴) and *Staphylococcus aureus infection*
(4/19 genes, FDR 1.4 × 10⁻³) — a direct reflection of the panel's MHC-heavy
composition at this stage (§3.6–§3.7). The male MR-prioritised gene set is
instead dominated by lymphocyte-activation biology, headed by positive
regulation of T-cell activation (5/24 genes, FDR 0.018) and *Herpes simplex
virus 1 infection* (4/16 genes, FDR 0.016). Because these annotations were
run on the MHC-retained primary MR-prioritised sets, as specified in §2.13,
no enrichment result reported here can be considered independent of the HLA
region, and the female result in particular should be read with that
caveat.

**3.14 Immune deconvolution and composition-adjusted expression**

Composition differs by disease status within each sex, following a
lymphopenia pattern: in women, CD8 T cells fall from a mean fraction of
0.162 in controls to 0.093 in RA (FDR 1.6 × 10⁻⁸), with smaller shifts in γδ
T cells and naive B cells; in men, CD8 T cells fall similarly (0.154 →
0.091, FDR 0.018) alongside increases in γδ T cells and eosinophils. Across
both sexes, **3 female and 4 male cell subsets differ by disease status at
FDR < 0.05**. Refitting the within-sex DEG contrasts with composition
principal components added reduces the DEG count substantially — the female
contrast falls from 5,131 to **2,709 DEGs (49.4% retained**, ρ = 0.889
between adjusted and unadjusted log-fold-changes**)**; the male contrast
falls from 5,820 to **1,450 DEGs (23.8% retained**, ρ = 0.842**)**. As stated
in §2.14, because leukocyte composition in RA is plausibly a disease
consequence rather than a confounder, this adjustment is read as removing
part of the true composition-mediated disease signal rather than as a bias
correction; the unadjusted model is retained as primary for the diagnostic
claim, and the adjusted model is the only basis for a cell-intrinsic
interpretation of a specific gene. Applied to the primary panel genes, **10
of 12 gene–sex pairs retain their signal after composition adjustment**; in
the male panel, *VPS52* (69.6% of logFC retained) and *INPP5B* (56.7%
retained) are **lost to composition adjustment**, indicating that their
apparent RA association is substantially explained by the shift in
leukocyte composition rather than by cell-intrinsic regulation.

**3.15 Nomogram construction and clinical evaluation**

A penalised logistic model was fitted on the consensus panel genes for each
sex on the full training cohort and converted to a nomogram, a
bootstrap-corrected calibration curve and a Vickers decision curve. The
female model (6 genes, n = 145: 86 RA, 59 HC) achieved a bootstrap-corrected
mean absolute calibration error of 0.036, mean squared error 0.00155 and a
90th-percentile absolute error of 0.053. The male model (6 genes, n = 38: 17
RA, 21 HC), fitted with the stabilising ridge penalty specified in §2.15,
achieved a mean absolute calibration error of 0.112, mean squared error
0.01286 and a 90th-percentile absolute error of 0.127 — a considerably worse
calibration than the female model, consistent with the small-sample caveats
attached to every other male-stratum result in this chapter. At a 20% risk
threshold, the female panel's net benefit is 0.5069 against 0.4914 for
treating every patient; at 50% it is 0.3655 against 0.1862. In the male
model, net benefit at 20% is 0.4013 against 0.3092 for treat-all, and at 50%
is 0.4211 against a *negative* treat-all net benefit of −0.1053 — at high
risk thresholds, treating every male patient becomes actively net-harmful
under this framework, while the panel-guided decision remains
net-beneficial. Across the tested threshold range (0.01–0.99), the panel's
net benefit exceeds both the treat-all strategy and zero in both sexes.

**Discussion**

**Chapter summary**