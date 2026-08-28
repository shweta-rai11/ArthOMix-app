# Cross-Omics — Thesis Implementation Paragraphs

Both paragraphs below are constrained to statements traceable to the inspected source code, cited
in `Cross-Omics_Dataset_Tab_Code_Audit.md` and `Cross-Omics_End_to_End_Function_Documentation.md`.
No function, dataset, statistical method, or output not located in the codebase is claimed here.
Style is informed by how this codebase's own METHODS documentation
(`data/preloaded/methylomics/tables/script0N_*/METHODS_*.md`) and its Cross-Omics module's own
extensive code comments describe implementation choices — plain, method-first, and explicit about
what is and is not statistically claimed — rather than by copying any external framework's wording.

## Implementation Paragraph

The Cross-Omics module (`R/crossomics/`) is organized as a Dataset tab and three registered
sub-modules — Expression and Methylation, Biomarker Convergence, and Cross-Omics MR — presented
across 21 `tabPanel` elements in total. The Dataset tab (`mod_cross_dataset.R`) is the entry point
for the Expression-and-Methylation analysis: it loads either this application's own bundled,
sex-stratified Transcriptomics differential-expression (DEG) and Methylomics differential-
methylation (DMP or DMR) results, or a user-uploaded file of the same shape, and standardizes both
onto one common schema (`gene, log2fc, pvalue, fdr` for expression; `cpg, gene, dbeta, pvalue, fdr`
for methylation) via `cx_standardize_expression()`/`cx_standardize_methylation()`
(`crossomics_integration_helpers.R`), publishing the result into a shared `cross_dataset` store
that only the Expression-and-Methylation sub-module reads. For the bundled data, the upstream
quality control, probe filtering, surrogate-variable adjustment, and bacon bias correction that
produced the underlying DMP/DMR statistics were performed by an external pipeline documented in
this repository's own methods notes (`METHODS_load_qc.md`, `METHODS_dmp_sva_sexstratified.md`,
`METHODS_dmr_sexstratified.md`) but not re-executed, or re-verifiable as source code, inside this
application; the Dataset tab itself performs no normalization or QC of its own. Gene-CpG
annotation for the CpG-level path is resolved directly from the Illumina 450K Bioconductor
manifest (`cx_get_region_annotation()`), taking the first listed gene for any multiply-annotated
probe, a limitation disclosed to the user in the run's own provenance text. On "Run Integration,"
gene identifiers are harmonized across both layers by exact symbol/Entrez/Ensembl match or
single-candidate alias resolution only (`cx_harmonize_gene_ids()` — no fuzzy matching), multiple
CpGs per gene are collapsed by one of seven user-selectable methods
(`cx_aggregate_methylation()`), with the mean/median methods combining each gene's constituent
CpG p-values via a directional Stouffer's Z statistic — so the reported gene-level effect size and
its significance are always derived from the same set of CpGs — followed by gene-level
Benjamini-Hochberg FDR correction. Genes are classified into four regulatory-direction categories
plus an Evidence Level tier, using association language throughout ("potential
methylation-associated repression/activation," never a causal claim), with the "Strong candidate"
tier additionally requiring a significant negative sample-level correlation computed only when
genuinely paired per-sample data is detected (`cx_detect_sample_pairing()`) — never inferred. The
Biomarker Convergence and Cross-Omics MR sub-modules operate independently of the Dataset tab,
relabeling — never re-deriving — the pipeline's own precomputed eQTL-MR/mQTL-MR/DEG/DMP/DMR join
and single-instrument Mendelian-randomization results at user-adjustable significance thresholds.
Outputs across the module include the integrated gene-level table, CpG-level detail, quadrant and
heatmap visualizations, a static gene-CpG network, evidence-combination category tables, and
downloadable CSV/TSV/XLSX results with an accompanying Markdown provenance report.

## Dataset Tab Paragraph

The Dataset tab receives either of two forms of the same standardized input: this application's
own bundled sex-stratified DEG and DMP/DMR result tables (read via `cx_load_default_deg()` and
`cx_load_default_methylation()`/`cx_load_default_dmr()`), or a user-uploaded CSV/TSV/TXT/XLSX file
for either omics layer, read and column-mapped automatically by `cx_read_and_detect()` and
`cx_detect_columns()`. Both paths are passed through the same standardization functions,
`cx_standardize_expression()` and `cx_standardize_methylation()`, which require a resolvable gene
identifier and a numeric effect-size column (log2FC or Δβ), coerce the remaining fields to numeric,
and drop rows with a missing gene value; expression-side duplicate gene symbols are collapsed to
the most significant row, while methylation-side duplicate CpGs per gene are deliberately left
uncollapsed at this stage. The only quality checks performed by the tab itself are these structural
presence/parseability checks — no beta-value plausibility range, missingness threshold, or
probe-quality filter is applied here; those checks, where they exist, were performed upstream of
this application for the bundled example data, and are documented, not re-executed, in this
repository. The object handed to the downstream Expression-and-Methylation sub-module is the
standardized `gene, log2fc, pvalue, fdr` and `cpg, gene, dbeta, pvalue, fdr, chr, pos, region,
region_fine, island_context` pair, published into the shared `cross_dataset` reactiveValues store
together with its provenance label and, for uploads only, the original wide table and its detected
per-sample columns — the latter being the sole basis on which any later sample-level correlation
can be computed.
