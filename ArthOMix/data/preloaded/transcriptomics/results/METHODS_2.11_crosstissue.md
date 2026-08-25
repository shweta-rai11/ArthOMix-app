## 2.11 Cross-tissue evaluation of the sex-stratified biomarker panels

### 2.11.1 Rationale and scope of the analysis

The panels described in §2.9 were derived in whole blood, a tissue that is
readily accessible but is not the site at which rheumatoid arthritis (RA)
pathology occurs. A cross-tissue analysis was therefore undertaken in RA synovial
tissue in order to establish whether the same genes remain informative at the site
of disease, and whether the direction of their association with disease is
preserved between the two compartments.

The scope of the resulting claim is deliberately narrow, and is stated at the
outset because it governs how every estimate in this section may be read. The
analysis is not an external validation of the blood classifier. The locked blood
model, with its blood-derived coefficients, was not transported into synovium;
both the panel-level and the gene-level synovial analyses re-estimate their
parameters within the synovial data themselves. What is tested is therefore the
portability of the gene set and the consistency of the direction of effect across
tissues, rather than the transferability of a fitted diagnostic model. An estimate
obtained in this manner cannot be quoted as out-of-sample performance of the blood
panel, and is not so quoted here. This distinction follows the conventional
hierarchy of prediction-model validation, in which the transport of a fully
specified model to a new setting is a stronger test than the re-estimation of a
model using a previously specified set of predictors (Justice, Covinsky and
Berlin, 1999; Steyerberg and Harrell, 2016).

### 2.11.2 Synovial cohort, sample selection and sex assignment

Cross-tissue evaluation used a publicly available synovial tissue RNA-sequencing
dataset (GSE89408). Samples were restricted to the RA and histologically normal
groups, with osteoarthritis, arthralgia and undifferentiated arthritis samples
excluded, so that the contrast tested in synovium corresponds to the
disease-versus-control contrast tested in blood. Disease group was determined from
the sample identifier prefix and sex from the series metadata field, and the
correspondence between the columns of the count matrix and the rows of the
metadata was verified for every sample before use rather than assumed. Analyses
were then conducted within each sex, with the female panel evaluated only in
female samples and the male panel only in male samples, in accordance with the
stratified design of §2.9.

The disease and control groups are markedly unbalanced in this cohort, control
samples forming a small minority of those retained. This imbalance is carried into
every synovial estimate reported below, is stated alongside those estimates, and
constrains both the precision of the sex-stratified analyses and the behaviour of
the cross-validation described in §2.11.5.

### 2.11.3 Expression processing and differential expression

Raw counts were assembled into a `DGEList` object and genes of insufficient
expression were removed using the filtering rule implemented in `filterByExpr`,
applied with disease group as the grouping factor so that a gene expressed in only
one group is retained (Chen, Lun and Smyth, 2016). Library sizes were normalised
by the trimmed mean of M-values method (Robinson and Oshlack, 2010) as implemented
in edgeR v4.4.2 (Robinson, McCarthy and Smyth, 2010), and counts were transformed
to log₂ counts per million with a prior count of one in order to stabilise the
variance of genes observed at low abundance.

Differential expression between RA and normal synovium was assessed within the
limma-voom framework (Law et al., 2014; Ritchie et al., 2015). The design matrix
contained sex and disease group, so that the disease coefficient is estimated
adjusted for sex. Precision weights were estimated by `voom`, gene-wise linear
models were fitted and their residual variances moderated by empirical Bayes
shrinkage (Smyth, 2004), and p-values were corrected across genes by the
Benjamini-Hochberg procedure (Benjamini and Hochberg, 1995). For each panel gene,
the synovial log₂ fold change, its adjusted p-value and the sign of the fold
change were extracted.

### 2.11.4 Concordance of effect direction between blood and synovium

Direction concordance was defined by comparing the sign of a gene's synovial log₂
fold change with the sign of its fold change in the blood training data for the
same sex, so that a gene shared between the two panels is assessed against the
training direction of the stratum in which it was selected. A gene was classified
as concordant when the two signs agreed, and as significantly concordant when it
was additionally differentially expressed in synovium at a Benjamini-Hochberg
false discovery rate below 0.05. Concordance was displayed as a scatter of each
gene's log₂ fold change in the two tissues, the quadrants of which distinguish
concordant from discordant genes.

### 2.11.5 Model training and discrimination in synovial tissue

For each sex, the panel genes present in the synovial expression matrix were
standardised to gene-wise z-scores within the synovial data and entered as
predictors into a multivariable logistic regression model of RA versus normal
tissue. Where standardisation produced an undefined value, the corresponding entry
was set to zero, that is, to the gene mean. The coefficients of this model were
estimated within the synovial data; only the identity of the predictors was
inherited from blood.

Two estimates of discrimination were computed. The apparent area under the
receiver operating characteristic curve (AUC) scores the model on the samples used
to fit it and is reported explicitly as an optimistic upper bound rather than as
performance (Harrell, Lee and Mark, 1996). The cross-validated AUC was obtained by
ten-fold cross-validation, the model being refitted on each training fold and
out-of-fold predicted probabilities collected across folds before a single
receiver operating characteristic curve was constructed. Confidence intervals of
95% were computed by the method of DeLong, DeLong and Clarke-Pearson (1988) using
pROC v1.19.0.1 (Robin et al., 2011), and the case direction was fixed a priori
rather than selected to maximise the AUC.

Two limitations of this cross-validation are stated rather than left implicit.
First, fold assignment was by simple random allocation rather than allocation
stratified on disease status, in contrast to the stratified schemes employed in
§2.9; in combination with the group imbalance noted in §2.11.2, individual folds
may therefore contain very few control samples, which inflates the variance of the
fold-level estimates. Second, only the logistic coefficients are re-estimated
across folds, the panel itself having been fixed in blood before these data were
examined. The latter is not a source of optimism within synovium, since no
synovial sample contributed to selecting the panel, but it does mean that this
estimate incorporates none of the selection-bias correction described in §2.9 and
is not comparable with a nested estimate (Ambroise and McLachlan, 2002; Varma and
Simon, 2006).

### 2.11.6 Univariate evaluation and the orientation convention

Each panel gene was additionally evaluated as a univariate classifier in synovium,
both across all retained samples and within each sex. The unit of evaluation
throughout is the gene-panel pair rather than the gene alone, because the female
and male panels share genes and a shared gene constitutes a different classifier
in each stratum, carrying a different estimate of discrimination; treating the
gene alone as the unit would collapse two distinct stratum-specific estimates
into one.

Two orientation conventions were used, and because they yield different values for
the same gene the distinction is recorded explicitly. In the synovial discovery
table and in the receiver operating characteristic overlay produced for synovium
alone, each gene was oriented so that its AUC is at least 0.5, the data themselves
determining which direction of expression marks disease; this convention addresses
how much information the gene carries within synovium irrespective of direction.
In the cross-dataset comparisons, by contrast, every gene was held to the
orientation fixed in the blood training data, a gene found concordant retaining
its synovial AUC and a gene found discordant having its AUC replaced by one minus
that value. Under the second convention an AUC below 0.5 is informative rather
than anomalous, being the signature of a gene whose direction of association
reverses between blood and synovium, and the corresponding figures are annotated
accordingly. Only the second, training-fixed convention was used when synovium was
placed alongside the training, internal holdout and external blood estimates, so
that the four datasets are compared under a single consistent rule and a reversal
of direction is made visible rather than concealed by silent re-orientation.

### 2.11.7 Presentation across datasets

Panel-level and gene-level results were assembled across all four datasets,
comprising the training cohort, the sealed internal holdout, the external blood
cohort and the external synovial cohort, using the receiver operating
characteristic coordinates stored during the blood analysis of §2.9 together with
the synovial curves computed here. Figures were produced with ggplot2 v4.0.2
(Wickham, 2016) and composed with patchwork.

### 2.11.8 Software

Analyses were conducted in R v4.4.2 using edgeR v4.4.2 (Robinson, McCarthy and
Smyth, 2010), limma v3.62.2 (Ritchie et al., 2015), Biobase v2.66.0, pROC
v1.19.0.1 (Robin et al., 2011) and data.table v1.18.2.1. A single global seed was
set before fold assignment, so that the cross-validated synovial estimates are
exactly reproducible.

---

### References cited in this section

Ambroise, C. and McLachlan, G.J. (2002) 'Selection bias in gene extraction on the basis of microarray gene-expression data', *Proceedings of the National Academy of Sciences*, 99(10), pp. 6562–6566.

Benjamini, Y. and Hochberg, Y. (1995) 'Controlling the false discovery rate: a practical and powerful approach to multiple testing', *Journal of the Royal Statistical Society: Series B*, 57(1), pp. 289–300.

Chen, Y., Lun, A.T.L. and Smyth, G.K. (2016) 'From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline', *F1000Research*, 5, 1438.

DeLong, E.R., DeLong, D.M. and Clarke-Pearson, D.L. (1988) 'Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach', *Biometrics*, 44(3), pp. 837–845.

Harrell, F.E., Lee, K.L. and Mark, D.B. (1996) 'Multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors', *Statistics in Medicine*, 15(4), pp. 361–387.

Justice, A.C., Covinsky, K.E. and Berlin, J.A. (1999) 'Assessing the generalizability of prognostic information', *Annals of Internal Medicine*, 130(6), pp. 515–524.

Law, C.W., Chen, Y., Shi, W. and Smyth, G.K. (2014) 'voom: precision weights unlock linear model analysis tools for RNA-seq read counts', *Genome Biology*, 15, R29.

Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W. and Smyth, G.K. (2015) 'limma powers differential expression analyses for RNA-sequencing and microarray studies', *Nucleic Acids Research*, 43(7), e47.

Robin, X., Turck, N., Hainard, A., Tiberti, N., Lisacek, F., Sanchez, J.-C. and Müller, M. (2011) 'pROC: an open-source package for R and S+ to analyze and compare ROC curves', *BMC Bioinformatics*, 12, 77.

Robinson, M.D., McCarthy, D.J. and Smyth, G.K. (2010) 'edgeR: a Bioconductor package for differential expression analysis of digital gene expression data', *Bioinformatics*, 26(1), pp. 139–140.

Robinson, M.D. and Oshlack, A. (2010) 'A scaling normalization method for differential expression analysis of RNA-seq data', *Genome Biology*, 11, R25.

Smyth, G.K. (2004) 'Linear models and empirical Bayes methods for assessing differential expression in microarray experiments', *Statistical Applications in Genetics and Molecular Biology*, 3, Article 3.

Steyerberg, E.W. and Harrell, F.E. (2016) 'Prediction models need appropriate internal, internal-external, and external validation', *Journal of Clinical Epidemiology*, 69, pp. 245–247.

Varma, S. and Simon, R. (2006) 'Bias in error estimation when using cross-validation for model selection', *BMC Bioinformatics*, 7, 91.

Wickham, H. (2016) *ggplot2: Elegant Graphics for Data Analysis*. 2nd edn. New York: Springer.

**Note on the dataset citation.** The primary publication for the synovial
RNA-sequencing cohort (GSE89408) must be inserted before submission; this
accession is flagged as unverified in the project reference list and no citation
for it has been fabricated here.
