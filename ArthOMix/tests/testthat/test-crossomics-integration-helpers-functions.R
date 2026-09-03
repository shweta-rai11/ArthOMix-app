## Module 4 (Cross-omics) - Expression-and-Methylation Integration's core
## pure engine (crossomics_integration_helpers.R): column auto-detection,
## standardization, region bucketing, gene-level methylation aggregation

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "functions", "integration", "crossomics_integration_helpers.R"))

test_that("cx_match_column() matches by regex pattern, excluding already-claimed columns", {
  cols <- c("Gene_Symbol", "log2FC", "P.Value")
  expect_equal(cx_match_column(cols, CX_FIELD_PATTERNS$gene), "Gene_Symbol")
  expect_true(is.na(cx_match_column(cols, CX_FIELD_PATTERNS$gene, exclude = "Gene_Symbol")))
})

test_that("cx_detect_columns() greedily/order-independently assigns each canonical field to a distinct column", {
  df <- data.frame(gene_symbol = 1, logFC = 1, adj.P.Val = 1, P.Value = 1)
  out <- cx_detect_columns(df, "expression")
  expect_equal(unname(out["gene"]), "gene_symbol")
  expect_equal(unname(out["log2fc"]), "logFC")
  expect_equal(unname(out["fdr"]), "adj.P.Val")
  expect_equal(unname(out["pvalue"]), "P.Value")
})

test_that("cx_detect_sample_columns() finds numeric non-metadata columns, excluding common summary-statistic names", {
  df <- data.frame(gene = c("A", "B"), logFC = c(1, 2), AveExpr = c(3, 4), S1 = c(5, 6), S2 = c(7, 8), stringsAsFactors = FALSE)
  mapping <- c(gene = "gene", log2fc = "logFC")
  out <- cx_detect_sample_columns(df, mapping)
  expect_setequal(out, c("S1", "S2"))
})

test_that("cx_standardize_expression() builds the canonical gene/log2fc/pvalue/fdr table and keeps unmapped extra columns", {
  df <- data.frame(sym = c("TP53", "BRCA1"), lfc = c(2, -1), pv = c(0.01, 0.02), notes = c("x", "y"), stringsAsFactors = FALSE)
  mapping <- c(gene = "sym", log2fc = "lfc", pvalue = "pv", fdr = NA_character_)
  out <- cx_standardize_expression(df, mapping)
  expect_true(out$ok)
  expect_equal(out$df$gene, c("TP53", "BRCA1"))
  expect_true("notes" %in% colnames(out$df))
})

test_that("cx_standardize_expression() refuses without a gene or log2fc column", {
  expect_false(cx_standardize_expression(data.frame(x = 1), c(gene = NA_character_, log2fc = "x"))$ok)
  expect_false(cx_standardize_expression(data.frame(x = 1), c(gene = "x", log2fc = NA_character_))$ok)
})

test_that("cx_dedup_by_gene() keeps the duplicate row with the smallest FDR, falling back to smallest p-value", {
  df <- data.frame(gene = c("A", "A", "B"), fdr = c(0.2, 0.01, NA), pvalue = c(0.5, 0.5, 0.3))
  out <- cx_dedup_by_gene(df)
  expect_equal(nrow(out), 2L)
  expect_equal(out$fdr[out$gene == "A"], 0.01)
})

test_that("cx_region_bucket() classifies Illumina UCSC_RefGene_Group and free-text region strings into Promoter/Gene body/Other", {
  expect_equal(cx_region_bucket(c("TSS200", "Body", "3'UTR", "intergenic", NA)), c("Promoter", "Gene body", "Gene body", "Other", NA))
  expect_equal(cx_region_bucket("promoter region"), "Promoter")
  expect_equal(cx_region_bucket("gene body"), "Gene body")
})

test_that("cx_region_fine() takes the first semicolon/comma-separated token, trimmed", {
  expect_equal(cx_region_fine(c("TSS200;Body", " 5'UTR ,1stExon", NA, "")), c("TSS200", "5'UTR", NA, NA))
})

test_that("cx_standardize_methylation() requires gene + (dbeta or beta), computes region/region_fine from region_raw", {
  df <- data.frame(cpg = "cg1", gene = "TP53", db = -0.3, reg = "TSS200", stringsAsFactors = FALSE)
  mapping <- c(cpg = "cpg", gene = "gene", dbeta = "db", beta = NA_character_, pvalue = NA_character_,
                fdr = NA_character_, chr = NA_character_, pos = NA_character_, end = NA_character_,
                n_cpgs = NA_character_, region_id = NA_character_, region = "reg", island = NA_character_)
  out <- cx_standardize_methylation(df, mapping)
  expect_true(out$ok)
  expect_equal(out$df$region, "Promoter")
  expect_equal(out$df$region_fine, "TSS200")
})

test_that("cx_standardize_methylation() refuses without a gene column, or without either dbeta or beta", {
  mapping_no_gene <- c(gene = NA_character_, cpg = "cpg", dbeta = "db", beta = NA_character_)
  expect_false(cx_standardize_methylation(data.frame(cpg = "cg1", db = 1), mapping_no_gene)$ok)
  mapping_no_val <- c(gene = "gene", cpg = "cpg", dbeta = NA_character_, beta = NA_character_)
  expect_false(cx_standardize_methylation(data.frame(cpg = "cg1", gene = "TP53"), mapping_no_val)$ok)
})

test_that("cx_cpg_level_table() flags sig_cpg correctly and labels methylation_direction from dbeta's sign", {
  meth_std <- data.frame(gene = c("A", "A"), cpg = c("cg1", "cg2"), dbeta = c(0.2, -0.05),
                           pvalue = c(0.001, 0.5), fdr = c(0.001, 0.5), region_raw = NA, region = NA, region_fine = NA,
                           island_context = NA, chr = NA, pos = NA)
  out <- cx_cpg_level_table(meth_std, meth_thresh = 0.1, meth_fdr_thresh = 0.05)
  expect_equal(out$methylation_direction, c("Hyper", "Hypo"))
  expect_equal(out$sig_cpg, c(TRUE, FALSE))
})

test_that("cx_cpg_counts_per_gene() correctly tallies per-gene CpG counts by significance/direction/region/island", {
  cpg_level <- data.frame(
    gene = c("A", "A", "A", "B"),
    cpg = c("cg1", "cg2", "cg3", "cg4"),
    region_fine = c("TSS200", "Body", "Body", "TSS1500"),
    island_context = c("Island", "OpenSea", "OpenSea", "Island"),
    dbeta = c(0.3, -0.2, 0.1, 0.4),
    sig_cpg = c(TRUE, TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  out <- cx_cpg_counts_per_gene(cpg_level)
  a_row <- out[out$gene == "A", ]
  expect_equal(a_row$n_cpg_total, 3L)
  expect_equal(a_row$n_cpg_significant, 2L)
  expect_equal(a_row$n_cpg_hyper_sig, 1L)
  expect_equal(a_row$n_body, 2L)
  expect_equal(a_row$primary_region, "Body")
})

test_that("cx_aggregate_methylation() 'mean' aggregation averages dbeta per gene and combines p-values via real Stouffer's Z (reduces to the single CpG's own p-value for a 1-CpG gene)", {
  meth_std <- data.frame(
    gene = c("A", "A", "B"), cpg = c("cg1", "cg2", "cg3"),
    dbeta = c(0.3, 0.1, -0.2), pvalue = c(0.001, 0.01, 0.02),
    fdr = c(0.001, 0.01, 0.02), chr = NA, pos = NA, region_raw = NA, region = NA, region_fine = NA, island_context = NA,
    stringsAsFactors = FALSE
  )
  out <- cx_aggregate_methylation(meth_std, method = "mean")
  expect_true(out$ok)
  a_row <- out$df[out$df$gene == "A", ]
  expect_equal(a_row$dbeta, mean(c(0.3, 0.1)))
  b_row <- out$df[out$df$gene == "B", ]
  expect_equal(b_row$pvalue, 0.02, tolerance = 1e-8)
  expect_true(b_row$n_probes == 1L)
  expect_true(a_row$n_probes == 2L)
})

test_that("cx_aggregate_methylation() 'min_fdr'/'max_abs_dbeta' pick the correct representative CpG per gene", {
  meth_std <- data.frame(
    gene = c("A", "A"), cpg = c("cg1", "cg2"), dbeta = c(0.1, -0.5),
    pvalue = c(0.5, 0.001), fdr = c(0.5, 0.001), chr = NA, pos = NA, region_raw = NA, region = NA, region_fine = NA, island_context = NA,
    stringsAsFactors = FALSE
  )
  out_fdr <- cx_aggregate_methylation(meth_std, method = "min_fdr")
  expect_equal(out_fdr$df$cpg, "cg2")
  out_dbeta <- cx_aggregate_methylation(meth_std, method = "max_abs_dbeta")
  expect_equal(out_dbeta$df$cpg, "cg2")
})

test_that("cx_aggregate_methylation() region-restricted methods ('promoter_only' etc.) filter CpGs before aggregating", {
  meth_std <- data.frame(
    gene = c("A", "A"), cpg = c("cg1", "cg2"), dbeta = c(0.3, 0.1),
    pvalue = c(0.01, 0.01), fdr = c(0.01, 0.01), chr = NA, pos = NA,
    region_raw = c("TSS200", "Body"), region = c("Promoter", "Gene body"), region_fine = c("TSS200", "Body"), island_context = NA,
    stringsAsFactors = FALSE
  )
  out <- cx_aggregate_methylation(meth_std, method = "promoter_only")
  expect_true(out$ok)
  expect_equal(out$df$n_probes, 1L)
  expect_equal(out$df$dbeta, 0.3)
})

test_that("cx_aggregate_methylation() refuses when no CpGs remain after a region filter", {
  meth_std <- data.frame(gene = "A", cpg = "cg1", dbeta = 0.3, pvalue = 0.01, fdr = 0.01, chr = NA, pos = NA,
                           region_raw = "Body", region = "Gene body", region_fine = "Body", island_context = NA, stringsAsFactors = FALSE)
  out <- cx_aggregate_methylation(meth_std, method = "promoter_only")
  expect_false(out$ok)
})

test_that("cx_classify() assigns the correct 4-way Hyper/Hypo x Up/Down category only when BOTH layers clear significance", {
  df <- data.frame(log2fc = c(2, -2, 2, -2, 0.1), expr_fdr = c(0.01, 0.01, 0.01, 0.01, 0.01),
                     dbeta = c(-0.3, 0.3, 0.3, -0.3, 0.3), meth_fdr = c(0.01, 0.01, 0.01, 0.01, 0.01))
  out <- cx_classify(df, expr_thresh = 1, expr_fdr_thresh = 0.05, meth_thresh = 0.1, meth_fdr_thresh = 0.05)
  expect_equal(as.character(out$category), c("Hypo + Up", "Hyper + Down", "Hyper + Up", "Hypo + Down", "Not significant"))
  expect_equal(out$category_label[1], CX_CATEGORY_LABELS[["Hypo + Up"]])
})

test_that("cx_classify_evidence() requires has_correlation=TRUE AND a significant negative correlation for 'Strong candidate'", {
  df <- data.frame(sig_expression = TRUE, sig_methylation = TRUE, log2fc = 2, dbeta = -0.3, correlation_r = -0.8, correlation_fdr = 0.001)
  expect_equal(as.character(cx_classify_evidence(df, has_correlation = TRUE)), "Strong candidate")
  expect_equal(as.character(cx_classify_evidence(df, has_correlation = FALSE)), "Moderate candidate")
})

test_that("cx_filter_by_category()/cx_filter_by_region()/cx_filter_by_island()/cx_filter_by_evidence()/cx_filter_by_min_cpg()/cx_filter_by_correlation_direction() each apply correctly and no-op safely when their column is absent", {
  df <- data.frame(
    gene = c("A", "B", "C"), sig_expression = c(TRUE, TRUE, FALSE), sig_methylation = c(TRUE, FALSE, FALSE),
    category = c("Hyper + Down", "Not significant", "Not significant"),
    primary_region = c("Body", "TSS200", NA), n_island = c(1, 0, 0), n_shore = c(0, 0, 0), n_shelf = c(0, 0, 0), n_open_sea = c(0, 1, 0),
    evidence_level = c("Strong candidate", "Expression-only", "Insufficient evidence"),
    n_cpg_total = c(5, 1, 0), correlation_r = c(-0.6, 0.2, NA), correlation_fdr = c(0.01, 0.2, NA),
    stringsAsFactors = FALSE
  )
  expect_equal(cx_filter_by_category(df, "sig_both")$gene, "A")
  expect_equal(cx_filter_by_category(df, "sig_expr_only")$gene, "B")
  expect_equal(cx_filter_by_region(df, "Body")$gene, "A")
  expect_equal(cx_filter_by_region(df, NULL), df)
  expect_equal(cx_filter_by_island(df, "Island")$gene, "A")
  expect_equal(cx_filter_by_evidence(df, "Strong candidate")$gene, "A")
  expect_equal(cx_filter_by_min_cpg(df, 3)$gene, "A")
  expect_equal(cx_filter_by_correlation_direction(df, "neg")$gene, "A")
  expect_equal(nrow(cx_filter_by_region(data.frame(gene = "A"), "Body")), 1L)
})

test_that("cx_build_gene_sample_matrix() averages duplicate-id rows into one, preserves unique-id rows unchanged", {
  df <- data.frame(id = c("A", "A", "B"), S1 = c(1, 3, 5), S2 = c(2, 4, 6))
  out <- cx_build_gene_sample_matrix(df, "id", c("S1", "S2"))
  expect_equal(unname(out["A", "S1"]), 2)
  expect_equal(unname(out["B", "S2"]), 6)
})

test_that("cx_detect_sample_pairing() is conservative: paired=FALSE without both sides having per-sample columns or enough overlap", {
  expect_false(cx_detect_sample_pairing(character(0), c("S1", "S2", "S3"))$paired)
  out <- cx_detect_sample_pairing(c("S1", "S2", "S3", "S4"), c("S3", "S4", "S5", "S6"), min_overlap = 3)
  expect_false(out$paired)
  out2 <- cx_detect_sample_pairing(c("S1", "S2", "S3"), c("S2", "S3", "S4"), min_overlap = 2)
  expect_true(out2$paired)
})

test_that("cx_gene_correlation() computes a real per-gene correlation matching stats::cor.test() directly", {
  set.seed(2000)
  samples <- paste0("S", 1:10)
  expr <- matrix(rnorm(20), 2, 10, dimnames = list(c("G1", "G2"), samples))
  meth <- matrix(rnorm(20), 2, 10, dimnames = list(c("G1", "G2"), samples))
  meth["G1", ] <- -expr["G1", ] + rnorm(10, sd = 0.05)
  out <- cx_gene_correlation(expr, meth, samples, method = "pearson")
  expect_true(out$ok)
  g1 <- out$df[out$df$gene == "G1", ]
  expect_true(g1$r < -0.9)
  ref <- stats::cor.test(expr["G1", samples], meth["G1", samples])
  expect_equal(g1$r, unname(ref$estimate))
  expect_equal(g1$p, ref$p.value)
})

test_that("cx_adjust_p() applies real BH/bonferroni adjustment matching stats::p.adjust() directly", {
  p <- c(0.001, 0.01, 0.02, 0.5)
  expect_equal(cx_adjust_p(p, "BH"), stats::p.adjust(p, "BH"))
  expect_equal(cx_adjust_p(p, "bonferroni"), stats::p.adjust(p, "bonferroni"))
})

test_that("cx_get_region_annotation() loads real chr/pos/gene/region/island columns from the real 450K manifest and caches on repeat calls", {
  skip_if_not(requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE), "450K annotation not installed")
  out1 <- cx_get_region_annotation("450K")
  expect_true(out1$ok)
  expect_true(all(c("chr", "pos", "gene", "region_raw", "island_context") %in% colnames(out1$anno)))
  out2 <- cx_get_region_annotation("450K")
  expect_identical(out1$anno, out2$anno)
})

test_that("cx_get_region_annotation() fails soft for an unconfigured array type", {
  out <- cx_get_region_annotation("not_a_real_array")
  expect_false(out$ok)
})

test_that("cx_load_default_deg() reads a real preloaded DEG table with real gene/logFC columns", {
  skip_if_not(METH_DATA_AVAILABLE, "preloaded data not available")
  d <- cx_load_default_deg("female")
  skip_if(is.null(d), "no preloaded DEG table for this deployment")
  expect_true(nrow(d) > 0)
  expect_true("gene" %in% colnames(d))
})

test_that("cx_load_default_methylation() reads and annotates a real preloaded DMP table into the standardized methylation schema", {
  skip_if_not(METH_DATA_AVAILABLE, "preloaded data not available")
  out <- cx_load_default_methylation("female")
  skip_if_not(out$ok, "preloaded methylation table not available")
  expect_true(nrow(out$df) > 0)
  expect_true(all(c("gene", "cpg", "dbeta", "region", "region_fine") %in% colnames(out$df)))
})

test_that("cx_load_default_dmr() reads a real preloaded DMR table with real region-level dbeta/gene annotation", {
  skip_if_not(METH_DATA_AVAILABLE, "preloaded data not available")
  out <- cx_load_default_dmr("female")
  skip_if_not(out$ok, "preloaded DMR table not available")
  expect_true(nrow(out$df) > 0)
  expect_true(all(c("gene", "dbeta") %in% colnames(out$df)))
})

test_that("cx_build_provenance() reports every real parameter used, '(not set)' for anything genuinely missing", {
  lines <- cx_build_provenance(list(sex_stratum = "female", expr_thresh = 1, expr_fdr_thresh = 0.05, meth_thresh = 0.1, meth_fdr_thresh = 0.05, agg_method = "mean"))
  text <- paste(lines, collapse = "\n")
  expect_true(grepl("Sex stratum: female", text))
  expect_true(grepl("Input mode: \\(not set\\)", text))
  expect_true(grepl("Mean .* across CpGs", text))
})

test_that("cx_harmonize_gene_ids() resolves a real Entrez ID and flags a genuinely ambiguous alias as 'ambiguous', never guessed", {
  skip_if_not(requireNamespace("org.Hs.eg.db", quietly = TRUE), "org.Hs.eg.db not installed")
  out <- cx_harmonize_gene_ids(c("TP53", "7157", "not_a_real_gene_xyz123"))
  expect_true(out$ok)
  expect_equal(out$df$canonical_symbol[out$df$input_id == "7157"], "TP53")
  expect_equal(out$df$match_type[out$df$input_id == "not_a_real_gene_xyz123"], "unmatched")
})

test_that("cx_apply_harmonization() rewrites only unambiguous matches, leaving ambiguous/unmatched entries as their original text", {
  harm_df <- data.frame(input_id = c("7157", "weird_id"), canonical_symbol = c("TP53", "A;B"),
                          match_type = c("exact_entrez", "ambiguous"), stringsAsFactors = FALSE)
  out <- cx_apply_harmonization(c("7157", "weird_id", "unrelated"), harm_df)
  expect_equal(out, c("TP53", "weird_id", "unrelated"))
})

test_that("cx_validate_dataset() reports real per-check pass/fail and readiness from whatever's actually loaded", {
  expr <- data.frame(gene = c("A", "B"), log2fc = c(1, NA), fdr = c(0.01, NA))
  meth <- data.frame(gene = c("A", "C"), dbeta = c(0.2, NA))
  out <- cx_validate_dataset(expr, meth)
  expect_true(out$ready)
  expect_true(out$transcriptomics[[1]]$ok)
  expect_true(out$compatibility[[2]]$ok)

  out_none <- cx_validate_dataset(NULL, NULL)
  expect_false(out_none$ready)
  expect_false(out_none$transcriptomics[[1]]$ok)
})
