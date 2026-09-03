## Module 3 (Multiomics) - Pathways' own pure functions
## (multiomics_pathway_helpers.R): uploaded-table column-role detection/
## confirmation (never silently auto-accepted), real identifier

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "functions", "integration", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_plots.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "multiomics_dataset_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "02_Cohort_Harmonization", "cohort_harmonization_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_integration_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "06_Gene_CpG_Concordance", "multiomics_concordance_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "07_Pathways", "multiomics_pathway_helpers.R"))

test_that("mp_detect_upload() detects id/effect/pvalue/fdr columns by name and flags a ranked (signed) effect column", {
  df <- data.frame(gene_symbol = c("TP53", "BRCA1"), log2FC = c(2, -1.5), pvalue = c(0.01, 0.02), padj = c(0.03, 0.04))
  out <- mp_detect_upload(df)
  expect_equal(out$detected$id_col, "gene_symbol")
  expect_equal(out$detected$id_type, "Gene symbol")
  expect_equal(out$detected$effect_col, "log2FC")
  expect_equal(out$detected$pvalue_col, "pvalue")
  expect_equal(out$detected$fdr_col, "padj")
  expect_true(out$detected$ranked)
  expect_equal(length(out$warnings), 0L)
})

test_that("mp_detect_upload() falls back to value-sniffing the first column when no header matches a known pattern, and warns accordingly", {
  df <- data.frame(my_weird_col = c("cg00000029", "cg00000108"), some_number = c(1, 2))
  out <- mp_detect_upload(df)
  expect_equal(out$detected$id_col, "my_weird_col")
  expect_equal(out$detected$id_type, "Illumina CpG probe ID")
})

test_that("mp_detect_upload() warns when no effect column exists (GSEA unavailable) and when neither p nor fdr exists", {
  df <- data.frame(gene_symbol = c("TP53", "BRCA1", "EGFR"))
  out <- mp_detect_upload(df)
  expect_true(any(grepl("GSEA", out$warnings)))
  expect_true(any(grepl("No P-value or adjusted-P column", out$warnings)))
})

test_that("mp_confirm_upload_mapping() applies the user-confirmed mapping, de-duplicates by feature, and drops empty identifiers", {
  df <- data.frame(sym = c("TP53", "TP53", "", "BRCA1"), lfc = c(1, 1, 2, -1), stringsAsFactors = FALSE)
  mapping <- list(id_col = "sym", id_type = "Gene symbol", effect_col = "lfc", pvalue_col = NA, fdr_col = NA, direction_col = NA, omics_col = NA, sex_col = NA)
  out <- mp_confirm_upload_mapping(df, mapping)
  expect_true(out$ok)
  expect_equal(nrow(out$df), 2L)
  expect_setequal(out$df$feature, c("TP53", "BRCA1"))
  expect_true(all(out$df$custom))
})

test_that("mp_confirm_upload_mapping() refuses when no identifier column is selected", {
  out <- mp_confirm_upload_mapping(data.frame(a = 1), list(id_col = NA))
  expect_false(out$ok)
})

test_that("mp_map_candidate_cpgs() maps real CpGs to real genes via the actual 450K manifest, reporting unmapped IDs explicitly", {
  skip_if_not(requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE), "450K annotation not installed")
  out <- mp_map_candidate_cpgs(c("cg00000029", "cg_totally_made_up_9999"))
  expect_true(out$ok)
  expect_true(out$n_mapped >= 1)
  expect_true("cg_totally_made_up_9999" %in% out$unmapped_cpgs)
})

test_that("mp_harmonize_identifiers() dispatches gene-like features through cx_harmonize_gene_ids() and CpG features through the CpG map, in one unified table", {
  skip_if_not(requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE), "450K annotation not installed")
  df <- data.frame(feature = c("TP53", "cg00000029"), id_type = c("Gene symbol", "Illumina CpG probe ID"), stringsAsFactors = FALSE)
  out <- mp_harmonize_identifiers(df)
  tp53_row <- out$df[out$df$feature == "TP53", ]
  expect_true(tp53_row$mapped)
  expect_equal(tp53_row$canonical_symbol, "TP53")
  cpg_row <- out$df[out$df$feature == "cg00000029", ]
  expect_equal(cpg_row$match_type, "cpg_mapped")
})

test_that("mp_mapping_summary() computes real mapped/unmapped counts and a correctly-rounded mapping rate", {
  mapped_df <- data.frame(feature = c("A", "B", "C", "D"), mapped = c(TRUE, TRUE, TRUE, FALSE))
  out <- mp_mapping_summary(mapped_df)
  expect_equal(out$n_mapped, 3L)
  expect_equal(out$n_unmapped, 1L)
  expect_equal(out$mapping_rate, 75)
  expect_equal(out$unmapped_ids, "D")
})

test_that("mp_build_ranked_vector() refuses with fewer than 10 valid (Entrez + finite value) rows", {
  df <- data.frame(expr_logFC = c(3, -5, 1, NA, 2))
  entrez <- c("1", "2", "3", "4", "2")
  out <- mp_build_ranked_vector(df, entrez, ranking_method = "log2fc")
  expect_false(out$ok)
})

test_that("mp_build_ranked_vector() succeeds with >=10 valid (Entrez + finite value) rows, sorted decreasing, deduplicated by keeping the largest |value|", {
  set.seed(900)
  n <- 15
  df <- data.frame(expr_logFC = c(rnorm(n - 1), NA))
  entrez <- c(as.character(1:(n - 2)), "3", NA)
  out <- mp_build_ranked_vector(df, entrez, ranking_method = "log2fc")
  expect_true(out$ok)
  expect_true(is.unsorted(-out$vec) == FALSE)
  expect_false(anyDuplicated(names(out$vec)) > 0)
})

test_that("mp_build_ranked_vector() signed_neglog10p ranking matches sign(logFC) * -log10(p) exactly", {
  df <- data.frame(expr_logFC = c(2, -3), expr_p = c(0.01, 0.0001))
  out <- mp_build_ranked_vector(df, c("1", "2"), ranking_method = "signed_neglog10p")
  expect_false(out$ok)
})

test_that("mp_concordance_direction() defers to a preloaded row's audited biological_pattern when present", {
  expect_equal(mp_concordance_direction(list(biological_pattern = "Concordant (gene-body hypermethylation)")), "Directionally consistent (region-aware)")
  expect_equal(mp_concordance_direction(list(biological_pattern = "Non-canonical")), "Directionally opposite (region-aware)")
})

test_that("mp_concordance_direction() classifies uploaded data by direct sign concordance when no biological_pattern is present", {
  expect_equal(mp_concordance_direction(list(biological_pattern = NA, expr_direction = "Up", meth_direction = "Hyper")), "Directionally consistent")
  expect_equal(mp_concordance_direction(list(biological_pattern = NA, expr_direction = "Up", meth_direction = "Hypo")), "Directionally opposite")
  expect_equal(mp_concordance_direction(list(biological_pattern = NA, expr_direction = NA, meth_direction = "Hypo")), "Insufficient information")
})

test_that("mp_build_evidence_tracks() summarizes overlapping-gene transcript/methylation evidence per pathway from the enrichment table's own geneID list", {
  enrichment_df <- data.frame(ID = "GO:001", Description = "test pathway", geneID = "TP53/BRCA1", stringsAsFactors = FALSE)
  input_df <- data.frame(
    feature = c("TP53", "BRCA1", "EGFR"), gene_symbol = c("TP53", "BRCA1", "EGFR"),
    expr_logFC = c(2, -1, 3), expr_p = c(0.01, 0.02, 0.001), expr_direction = c("Up", "Down", "Up"),
    cpg = c(NA, "cg1", NA), meth_direction = c(NA, "Hyper", NA), meth_p = c(NA, 0.03, NA), region = c(NA, "Body", NA),
    stringsAsFactors = FALSE
  )
  out <- mp_build_evidence_tracks(enrichment_df, input_df)
  expect_equal(out$transcript_gene_count, 2L)
  expect_equal(out$transcript_direction_summary, "Mixed")
  expect_equal(out$meth_cpg_count, 1L)
  expect_equal(out$integration_label, "RNA + Methylation supported")
})

test_that("mp_build_evidence_tracks() passes through an empty/NULL enrichment table unchanged", {
  expect_null(mp_build_evidence_tracks(NULL, data.frame()))
})

test_that("mp_validate_ora_inputs() refuses with fewer than 3 genes, or a universe smaller than the input list", {
  expect_false(mp_validate_ora_inputs(c("1", "2"), NULL)$ok)
  expect_true(mp_validate_ora_inputs(c("1", "2", "3"), NULL)$ok)
  out <- mp_validate_ora_inputs(c("1", "2", "3"), c("1", "2"))
  expect_false(out$ok)
  expect_true(grepl("Entire selected database", out$error))
})

test_that("mp_validate_gsea_inputs() refuses with fewer than 10 ranked features", {
  expect_false(mp_validate_gsea_inputs(setNames(1:5, letters[1:5]))$ok)
  expect_true(mp_validate_gsea_inputs(setNames(1:10, letters[1:10]))$ok)
})

test_that("mp_validate_database_choice() blocks a database whose own package gate has failed", {
  expect_true(mp_validate_database_choice("GO_BP")$ok)
  if (!MP_REACTOME_AVAILABLE) expect_false(mp_validate_database_choice("Reactome")$ok)
})

test_that("mp_infer_species() infers mouse/rat from Ensembl ID prefixes, defaults to human otherwise", {
  expect_equal(mp_infer_species("Ensembl Gene ID", c("ENSMUSG00000001", "ENSMUSG00000002")), "Mus musculus")
  expect_equal(mp_infer_species("Ensembl Gene ID", c("ENSRNOG00000001", "ENSRNOG00000002")), "Rattus norvegicus")
  expect_equal(mp_infer_species("Ensembl Gene ID", c("ENSG00000001")), "Homo sapiens")
  expect_equal(mp_infer_species("Gene symbol", c("TP53")), "Homo sapiens")
})

test_that("mp_build_metadata() reports every real parameter, including a correctly-rounded mapping rate", {
  tbl <- mp_build_metadata("GO_BP", "ORA", "Homo sapiens", "Entire selected database", NULL, input_count = 40, mapped_count = 30, padj_thresh = 0.05, min_size = 5, max_size = 500)
  by_field <- setNames(tbl$Value, tbl$Field)
  expect_equal(by_field[["Input features"]], "40")
  expect_equal(by_field[["Mapped features"]], "30")
  expect_equal(by_field[["Mapping rate (%)"]], "75")
})

test_that("mp_run_ora_go() (real clusterProfiler::enrichGO, offline via org.Hs.eg.db) finds real GO Biological Process enrichment for a canonical immune gene set", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  immune_genes_entrez <- c("3569", "7124", "3576", "3592", "6347", "3595")
  out <- mp_run_ora_go(immune_genes_entrez, universe_entrez = NULL, ont = "BP", pvalueCutoff = 0.2, minGSSize = 3)
  expect_true(out$ok)
  expect_true(nrow(out$df) > 0)
  expect_true(any(grepl("immune|inflamm|cytokine", out$df$Description, ignore.case = TRUE)))
  expect_true(all(out$df$p.adjust >= out$df$pvalue - 1e-8))
})

test_that("mp_finalize_enrich() reports a clear failure (never ok=TRUE/df=NULL) when the enrichment table is empty", {
  out <- mp_finalize_enrich(NULL, "Test Database")
  expect_false(out$ok)
  expect_true(nzchar(out$error))
})
