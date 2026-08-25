# 2.13 Functional interpretation

Functional interpretation was executed at two distinct stages within the
analytical pipeline: initially on the differentially expressed genes (DEGs)
following the expression analysis, and subsequently on the
Mendelian-randomisation prioritised genes after the causal-inference phase. Both
stages employed the clusterProfiler framework (Yu et al., 2012; Wu et al., 2021),
utilising human annotation from org.Hs.eg.db (Carlson, 2019), with gene symbols
consistently mapped to Entrez identifiers using `bitr`; symbols without an Entrez
mapping are dropped at this step, so both the tested gene sets and the background
universes are the mapped subsets rather than the nominal ones. Two further
over-representation analyses, of the WGCNA disease-module genes and of the
sex-by-diagnosis interaction genes, are described in their own sections and are
not part of the two stages reported here. Over-representation analysis and ranked
gene-set enrichment answer different questions, the former testing whether a
fixed gene list is enriched relative to a background and the latter testing
whether a pathway is systematically displaced in a ranking of all genes, and both
were therefore retained rather than treated as alternatives (Goeman and Bühlmann,
2007).

In the differential-expression phase, each significant DEG set was calculated for
the pooled cohort and separately for females and males before
over-representation analysis. Gene Ontology enrichment was conducted using
`enrichGO` across all three ontologies (biological process, cellular component
and molecular function) (Ashburner et al., 2000; The Gene Ontology Consortium,
2021), while KEGG pathway enrichment was performed with `enrichKEGG` for *Homo
sapiens* (organism = "hsa") (Kanehisa and Goto, 2000). Both analyses applied the
Benjamini–Hochberg correction with a cutoff of 0.05 imposed on the raw and
adjusted p-values alike and a q-value threshold of 0.2 (Benjamini and Hochberg,
1995; Storey and Tibshirani, 2003), with Gene Ontology results converted to
readable gene symbols while the KEGG identifiers were retained in Entrez form at
this stage. No background universe was supplied at this stage, so the comparison
is made against the full set of annotated human genes rather than against the
genes detectable on the arrays analysed here; this is more permissive than the
universe-restricted analysis described below, since gene classes that are
disproportionately represented among array-detectable transcripts can be enriched
for that reason alone, and the two stages are consequently not interchangeable
(Timmons, Szkop and Gallagher, 2015; Wijesooriya et al., 2022). Additionally, a
ranked Gene Set Enrichment Analysis of KEGG pathways was conducted using
`gseKEGG` (organism = "hsa", Benjamini–Hochberg correction, p-value cutoff 0.05,
`seed = TRUE` so that the permutation procedure is reproducible) (Subramanian et
al., 2005). For each comparison, genes were ranked by their limma moderated
t-statistic for the RA-versus-healthy-control contrast, retaining a single ranking
value per gene, namely the largest absolute t-statistic where a gene mapped to
multiple identifiers (Smyth, 2004; Ritchie et al., 2015). A positive normalised
enrichment score indicated a pathway up-regulated in RA, while a negative score
indicated down-regulation in RA.

In the causal-gene phase, the Mendelian-randomisation prioritised gene sets for
each sex were annotated separately to elucidate the biological characteristics of
genes with genetic support for a causal role in RA. The primary,
MHC-retained candidate sets were used, so genes of the major histocompatibility
complex contribute to any enriched immune terms; the MHC-free sensitivity sets
defined elsewhere in this chapter were not separately annotated, and no enrichment
result reported here can therefore be assumed to be independent of the HLA region.
The background universe was confined to genes expressed in the training cohort,
that is, the rows of the training expression matrix mapped to Entrez identifiers,
rather than the entire genome, ensuring that enrichment was assessed against
genes measurable in the study (Timmons, Szkop and Gallagher, 2015). Gene Ontology
biological-process enrichment was performed using `enrichGO` (ont = "BP") with
this expressed-gene universe, a p-value cutoff of 0.05, a q-value cutoff of 0.2
and readable output, and KEGG pathway enrichment was conducted with `enrichKEGG`
(organism = "hsa", same universe, p-value cutoff 0.05, the default q-value cutoff
of 0.2 applying), both under the default Benjamini–Hochberg adjustment, with
results mapped back to gene symbols. The most significant enriched terms per sex
were visualised as gene-ratio dotplots, combining the top ten Gene Ontology
biological-process terms and the top ten KEGG terms ranked by adjusted p-value,
faceted by source, with the gene ratio on the horizontal axis, point size
proportional to the number of genes in the term and point colour encoding the
adjusted p-value.

---

## Details corrected against the scripts

1. **The DEG-phase enrichment uses no background universe.** `enrichGO` and
   `enrichKEGG` at [05_dge.R:176-181](scripts/00_shared/05_dge.R#L176-L181) are
   called without a `universe` argument, so the background is the whole annotated
   genome, whereas the causal-gene phase explicitly restricts the universe to
   expressed genes ([34_pathway_enrichment.R:20-21, 29-32](scripts/goal2_sex_stratified/34_pathway_enrichment.R#L20-L32)).
   The passage described the restriction only for the second stage, which is
   correct as far as it goes, but the asymmetry needs declaring: two stages of the
   same thesis test enrichment against two different backgrounds, and the
   unrestricted one is the more permissive of the two.
2. **The 0.05 cutoff is not simply "a p-value threshold".** In clusterProfiler
   `pvalueCutoff` is applied to the raw *and* the adjusted p-value, with
   `qvalueCutoff` applied on top; `PADJ <- 0.05` at
   [05_dge.R:39](scripts/00_shared/05_dge.R#L39) is passed to that argument.
   Phrased accordingly, so the filter is not understood as uncorrected.
3. **GSEA reproducibility.** `gseKEGG(..., seed = TRUE)`
   ([05_dge.R:229-231](scripts/00_shared/05_dge.R#L229-L231)). The permutation
   test is stochastic and the seed is what makes the reported pathways
   reproducible; the original passage omitted it.
4. **`readable` applies to GO only in stage 1.** KEGG gene identifiers stay in
   Entrez form there ([05_dge.R:179-181](scripts/00_shared/05_dge.R#L179-L181)),
   whereas stage 2 wraps `enrichKEGG` in `setReadable`
   ([34:31-32](scripts/goal2_sex_stratified/34_pathway_enrichment.R#L31-L32)).
   Now stated for each stage rather than once.
5. **KEGG's q-value cutoff in stage 2 is not switched off.** Only
   `pvalueCutoff = 0.05` is passed, so the package default `qvalueCutoff = 0.2`
   still applies. Reporting the p cutoff alone understates the filtering.
6. **Unmapped symbols are silently dropped** by `bitr` under `suppressWarnings`,
   in both stages, so the effective gene sets and universes are smaller than the
   nominal ones. Worth one clause, because it is the first place gene counts stop
   agreeing with the DE tables.
7. **The causal-gene phase annotates the MHC-retained sets only.** Script 34
   reads `FS_input_{female,male}.csv`, not the `_noMHC` variants
   ([34:25](scripts/goal2_sex_stratified/34_pathway_enrichment.R#L25)). Given that
   this chapter's headline panel is the MHC-free one, an enrichment result driven
   by HLA-region genes cannot be quoted in support of it without this caveat.
8. **"Two distinct stages" is not true of the pipeline as a whole.** Enrichment is
   also run on the WGCNA disease-module genes
   ([06_WGCNA.R:423-427](scripts/00_shared/06_WGCNA.R#L423-L427), ont = "ALL",
   universe-restricted, q < 0.05) and on the sex-interaction genes
   ([05d_interaction_report.R:275-281](scripts/00_shared/05d_interaction_report.R#L275-L281),
   ont = "BP", universe-restricted, q < 0.2). One clause added so the claim is
   accurate; those two analyses belong to their own sections.
9. **Two things to action before this section is submitted.**
   - The stage-2 outputs are **not present on disk**: no
     `enrich_GO_BP_{female,male}.csv`, no `enrich_KEGG_{female,male}.csv` and no
     `fig_enrich_*` figures. Either script 34 has not been run in this tree, or it
     ran and found no significant terms, in which case it writes nothing and
     prints "no significant terms". The paragraph above currently describes an
     analysis with no artefacts behind it. Re-run it and confirm which case
     applies.
   - `34_pathway_enrichment.R:6` states "per sex (Female 74, Male 55)"; the
     candidate sets it actually reads are much smaller. Stale comment, same class
     of defect as in scripts 14 and 18.
10. Gene counts, term counts and enriched pathways are results and are
    deliberately absent above.

---

## Parameters used

| Item | Value |
|---|---|
| Framework | clusterProfiler with org.Hs.eg.db; symbols → Entrez via `bitr` (unmapped dropped) |
| **Stage 1: DEG over-representation** | |
| Gene sets tested | significant DEGs for pooled cohort, females, males (separately) |
| GO | `enrichGO`, ont = "ALL" (BP + CC + MF), `readable = TRUE` |
| KEGG | `enrichKEGG`, organism = "hsa", identifiers left as Entrez |
| Background universe | none supplied → all annotated human genes |
| Multiple testing | `pAdjustMethod = "BH"`, `pvalueCutoff = 0.05`, `qvalueCutoff = 0.2` |
| **Stage 1: ranked GSEA** | |
| Method | `gseKEGG`, organism = "hsa", `pAdjustMethod = "BH"`, `pvalueCutoff = 0.05`, `seed = TRUE` |
| Ranking statistic | limma moderated t for RA vs HC, descending |
| Duplicate handling | one rank per Entrez ID, keeping the largest \|t\| |
| Sign convention | NES > 0 = up in RA; NES < 0 = down in RA |
| **Stage 2: MR-prioritised gene enrichment** | |
| Gene sets tested | `FS_input_{female,male}.csv` (MHC-retained candidate sets) |
| Background universe | Entrez-mapped rows of the training expression matrix |
| GO | `enrichGO`, ont = "BP", `pvalueCutoff = 0.05`, `qvalueCutoff = 0.2`, `readable = TRUE` |
| KEGG | `enrichKEGG`, organism = "hsa", same universe, `pvalueCutoff = 0.05`, default `qvalueCutoff = 0.2`, `setReadable` applied |
| Multiple testing | package default `pAdjustMethod = "BH"` |
| Figure | dotplot, top 10 GO:BP + top 10 KEGG by adjusted p; x = gene ratio, size = gene count, colour = adjusted p; faceted by source |

---

## References cited in this section

Ashburner, M., Ball, C.A., Blake, J.A., Botstein, D., Butler, H., Cherry, J.M., et al. (2000) 'Gene Ontology: tool for the unification of biology', *Nature Genetics*, 25(1), pp. 25–29.

Benjamini, Y. and Hochberg, Y. (1995) 'Controlling the false discovery rate: a practical and powerful approach to multiple testing', *Journal of the Royal Statistical Society: Series B*, 57(1), pp. 289–300.

Carlson, M. (2019) *org.Hs.eg.db: Genome Wide Annotation for Human*. R package version 3.8.2. Bioconductor.

Goeman, J.J. and Bühlmann, P. (2007) 'Analyzing gene expression data in terms of gene sets: methodological issues', *Bioinformatics*, 23(8), pp. 980–987.

Kanehisa, M. and Goto, S. (2000) 'KEGG: Kyoto Encyclopedia of Genes and Genomes', *Nucleic Acids Research*, 28(1), pp. 27–30.

Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W. and Smyth, G.K. (2015) 'limma powers differential expression analyses for RNA-sequencing and microarray studies', *Nucleic Acids Research*, 43(7), e47.

Smyth, G.K. (2004) 'Linear models and empirical Bayes methods for assessing differential expression in microarray experiments', *Statistical Applications in Genetics and Molecular Biology*, 3(1), 3.

Storey, J.D. and Tibshirani, R. (2003) 'Statistical significance for genomewide studies', *Proceedings of the National Academy of Sciences*, 100(16), pp. 9440–9445.

Subramanian, A., Tamayo, P., Mootha, V.K., Mukherjee, S., Ebert, B.L., Gillette, M.A., et al. (2005) 'Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles', *Proceedings of the National Academy of Sciences*, 102(43), pp. 15545–15550.

The Gene Ontology Consortium (2021) 'The Gene Ontology resource: enriching a GOld mine', *Nucleic Acids Research*, 49(D1), pp. D325–D334.

Timmons, J.A., Szkop, K.J. and Gallagher, I.J. (2015) 'Multiple sources of bias confound functional enrichment analysis of globally truncated genes', *Genome Biology*, 16, 186.

Wijesooriya, K., Jadaan, S.A., Perera, K.L., Kaur, T. and Ziemann, M. (2022) 'Urgent need for consistent standards in functional enrichment analysis', *PLoS Computational Biology*, 18(3), e1009935.

Wu, T., Hu, E., Xu, S., Chen, M., Guo, P., Dai, Z., et al. (2021) 'clusterProfiler 4.0: a universal enrichment tool for interpreting omics data', *The Innovation*, 2(3), 100141.

Yu, G., Wang, L.-G., Han, Y. and He, Q.-Y. (2012) 'clusterProfiler: an R package for comparing biological themes among gene clusters', *OMICS: A Journal of Integrative Biology*, 16(5), pp. 284–287.
