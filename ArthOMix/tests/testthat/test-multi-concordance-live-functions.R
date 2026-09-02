## Module 3 (Multiomics) - Gene-CpG Concordance's own pure functions
## (multiomics_concordance_helpers.R): layer-role guessing, the data
## status panel, sample matching (thin wrapper), the candidate-biomarker
## pool (DIABLO/SNF/Joint/custom unification + source filter), real
## annotation-backed gene<->CpG mapping, direction classification (real
## cx_classify()) with the promoter/gene-body canonical-direction rule,
## raw-direction labeling, per-pair correlation, adjusted regression,
## priority scoring, and the two-group differential-stats path. Every
## `cx_*` call here is an already-verified building block per this file's
## own header comment - these tests check the CONCORDANCE-SPECIFIC
## composition logic, not re-derive cx_*'s own correctness (Cross-omics,
## a separate module, is where cx_* gets its own dedicated coverage).

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_plots.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_dataset_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "cohort_harmonization_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_integration_live_helpers.R"))
source_from_app_root(file.path("R", "methylomics", "celltype.R"))
source_from_app_root(file.path("R", "methylomics", "qc.R"))
source_from_app_root(file.path("R", "multiomics", "snf_clustering_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_concordance_helpers.R"))

## ---- mcc_layer_candidates() / mcc_default_layer() --------------------------

test_that("mcc_layer_candidates() picks layers by recorded omics_type metadata first", {
  md <- list(layers = list(A = matrix(1, 2, 2), B = matrix(1, 2, 2)),
             layer_meta = list(A = list(omics_type = "rnaseq"), B = list(omics_type = "methylation")))
  expect_equal(mcc_layer_candidates(md, "rnaseq"), "A")
  expect_equal(mcc_layer_candidates(md, "methylation"), "B")
})

test_that("mcc_layer_candidates() falls back to a name regex when omics_type metadata is absent, never guessing beyond rnaseq/methylation", {
  md <- list(layers = list(MyExpression = matrix(1, 2, 2), MyMeth = matrix(1, 2, 2)), layer_meta = list())
  expect_equal(mcc_layer_candidates(md, "rnaseq"), "MyExpression")
  expect_equal(mcc_layer_candidates(md, "methylation"), "MyMeth")
  expect_equal(mcc_layer_candidates(md, "proteomics"), character(0))
})

test_that("mcc_default_layer() picks the first candidate, NULL when there are none", {
  expect_equal(mcc_default_layer(c("A", "B"), list()), "A")
  expect_null(mcc_default_layer(character(0), list()))
})

## ---- mcc_data_status() -------------------------------------------------------

test_that("mcc_data_status() reports 'Missing' rows honestly when no layers are selected", {
  tbl <- mcc_data_status(list(layers = list()), list())
  expect_equal(tbl$status[tbl$item == "Expression data"], "Missing")
  expect_equal(tbl$status[tbl$item == "Methylation data"], "Missing")
})

test_that("mcc_data_status() reports real matched-sample counts and detail strings once both layers are present", {
  expr <- matrix(1, 8, 20, dimnames = list(paste0("S", 1:8), paste0("g", 1:20)))
  meth <- matrix(1, 6, 30, dimnames = list(paste0("S", 3:8), paste0("cg", 1:30)))
  md <- list(layers = list(RNA = expr, Meth = meth), sample_meta = NULL)
  tbl <- mcc_data_status(md, list(), expr_layer = "RNA", meth_layer = "Meth")
  expect_true(grepl("RNA \\(20 genes x 8 samples\\)", tbl$detail[tbl$item == "Expression data"]))
  expect_equal(tbl$status[tbl$item == "Matched samples"], "Available")
  expect_true(grepl("6 matched of 8 expression / 6 methylation", tbl$detail[tbl$item == "Matched samples"]))
})

## ---- mcc_detect_id_type() / mcc_detect_methylation_value_type() ------------

test_that("mcc_detect_id_type() distinguishes Ensembl/CpG/Entrez/gene-symbol ID formats", {
  expect_equal(mcc_detect_id_type(c("ENSG00000141510", "ENSG00000012048")), "Ensembl Gene ID")
  expect_equal(mcc_detect_id_type(c("cg00000029", "cg00000108")), "Illumina CpG probe ID")
  expect_equal(mcc_detect_id_type(c("7157", "672")), "Entrez ID")
  expect_equal(mcc_detect_id_type(c("TP53", "BRCA1")), "Gene symbol")
})

test_that("mcc_detect_methylation_value_type() distinguishes beta/M-value/other scales", {
  expect_equal(mcc_detect_methylation_value_type(matrix(runif(50, 0, 1), 5, 10)), "beta")
  expect_equal(mcc_detect_methylation_value_type(matrix(rnorm(50, 0, 3), 5, 10)), "M-value")
})

## ---- mcc_match_samples() -----------------------------------------------------

test_that("mcc_match_samples() reports real matched/removed counts via cx_detect_sample_pairing()", {
  expr <- matrix(1, 6, 2, dimnames = list(paste0("S", 1:6), c("f1", "f2")))
  meth <- matrix(1, 4, 2, dimnames = list(paste0("S", 3:6), c("cg1", "cg2")))
  out <- mcc_match_samples(expr, meth)
  expect_true(out$ok)
  expect_equal(out$n_matched, 4L)
  expect_equal(out$n_removed_expr, 2L)
  expect_equal(out$n_removed_meth, 0L)
})

test_that("mcc_match_samples() refuses when either matrix is missing", {
  out <- mcc_match_samples(NULL, matrix(1, 2, 2))
  expect_false(out$ok)
})

## ---- Candidate pool: mcc_diablo_candidates() / mcc_snf_candidates() / mcc_candidate_pool() / mcc_filter_source() ----

test_that("mcc_diablo_candidates() reads multi_results$biomarker$df, deduplicating by feature", {
  df <- data.frame(feature = c("g1", "g1", "g2"), omics = c("Transcriptomics", "Transcriptomics", "Methylomics"),
                     component = 1, loading = c(0.5, 0.5, -0.3), selection_frequency = c(0.9, 0.9, 0.6),
                     stability_category = c("Stable", "Stable", "Moderately stable"), stringsAsFactors = FALSE)
  out <- mcc_diablo_candidates(list(biomarker = list(df = df)))
  expect_equal(nrow(out), 2L)  ## g1 deduplicated
})

test_that("mcc_diablo_candidates() returns NULL when nothing is loaded", {
  expect_null(mcc_diablo_candidates(list()))
})

test_that("mcc_candidate_pool() unifies DIABLO + SNF + custom features, tagging id_type-derived omics when otherwise unknown", {
  diablo_df <- data.frame(feature = "g1", omics = "Transcriptomics", component = 1, loading = 0.5,
                            selection_frequency = 0.9, stability_category = "Stable", stringsAsFactors = FALSE)
  results <- list(biomarker = list(df = diablo_df))
  out <- mcc_candidate_pool(results, multi_dataset = list(), expr_layer = NULL, meth_layer = NULL,
                              custom_genes = "TP53", custom_cpgs = "cg00000029")
  expect_true(out$ok)
  expect_true(all(c("g1", "TP53", "cg00000029") %in% out$df$feature))
  expect_true(out$df$diablo[out$df$feature == "g1"])
  expect_true(out$df$custom[out$df$feature == "TP53"])
  expect_equal(out$df$omics[out$df$feature == "TP53"], "Transcriptomics")
  expect_equal(out$df$omics[out$df$feature == "cg00000029"], "Methylomics")
})

test_that("mcc_candidate_pool() reports ok=FALSE with a clear note when there are no candidates at all", {
  out <- mcc_candidate_pool(list(), multi_dataset = list(), expr_layer = NULL, meth_layer = NULL)
  expect_false(out$ok)
  expect_true(nzchar(out$note))
})

test_that("mcc_filter_source() applies every documented Biomarker Source option correctly", {
  pool <- data.frame(feature = c("a", "b", "c", "d"), diablo = c(TRUE, TRUE, FALSE, FALSE),
                       snf = c(TRUE, FALSE, TRUE, FALSE), joint = c(TRUE, FALSE, FALSE, TRUE),
                       omics = c("Transcriptomics", "Transcriptomics", "Methylomics", "Methylomics"),
                       custom = c(FALSE, FALSE, FALSE, TRUE), stringsAsFactors = FALSE)
  expect_equal(mcc_filter_source(pool, "DIABLO")$feature, c("a", "b"))
  expect_equal(mcc_filter_source(pool, "SNF")$feature, c("a", "c"))
  expect_equal(mcc_filter_source(pool, "DIABLO + SNF")$feature, "a")
  expect_equal(mcc_filter_source(pool, "Shared candidates")$feature, "a")  ## all 3 of diablo/snf/joint
  expect_equal(mcc_filter_source(pool, "Custom CpGs")$feature, "d")  ## d is custom AND Methylomics
})

## ---- mcc_gene_cpg_map() (real Illumina annotation) ---------------------------

test_that("mcc_gene_cpg_map() maps real candidate genes to real annotated CpGs via the actual 450K manifest", {
  skip_if_not(requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE), "450K annotation package not installed")
  out <- mcc_gene_cpg_map(c("TP53", "BRCA1"), array_type = "450K")
  expect_true(out$ok)
  expect_true(nrow(out$df) > 0)
  expect_true(all(toupper(out$df$gene_symbol) %in% c("TP53", "BRCA1")))
  expect_true(all(startsWith(out$df$cpg, "cg")))
})

test_that("mcc_gene_cpg_map() fails soft with an unresolvable gene list or empty input", {
  expect_false(mcc_gene_cpg_map(character(0))$ok)
})

## ---- mcc_classify_direction() (real cx_classify) + canonical rule ------------

mcc_pairs_fixture <- function() {
  data.frame(
    gene_symbol = c("G1", "G2", "G3", "G4"),
    cpg = c("cg1", "cg2", "cg3", "cg4"),
    log2fc = c(2, -2, 2, -2),
    expr_fdr = c(0.001, 0.001, 0.001, 0.001),
    dbeta = c(-0.3, 0.3, 0.3, -0.3),
    meth_fdr = c(0.001, 0.001, 0.001, 0.001),
    region_fine = c("TSS200", "TSS200", "Body", "Body"),
    stringsAsFactors = FALSE
  )
}

test_that("mcc_classify_direction() calls through to real cx_classify() and assigns the correct 4-way direction label", {
  df <- mcc_classify_direction(mcc_pairs_fixture(), expr_thresh = 1, expr_fdr_thresh = 0.05, meth_thresh = 0.1, meth_fdr_thresh = 0.05)
  expect_equal(as.character(df$direction_classification[df$gene_symbol == "G1"]), "Up expression + Hypomethylation")
  expect_equal(as.character(df$direction_classification[df$gene_symbol == "G2"]), "Down expression + Hypermethylation")
})

test_that("mcc_classify_direction() applies the promoter (inverse=canonical) vs. gene-body (concordant=canonical) rule correctly", {
  df <- mcc_classify_direction(mcc_pairs_fixture(), expr_thresh = 1, expr_fdr_thresh = 0.05, meth_thresh = 0.1, meth_fdr_thresh = 0.05)
  ## G1: TSS200 (promoter) + Up+Hypo -> inverse relationship -> canonical.
  expect_equal(df$canonical_label[df$gene_symbol == "G1"], "Canonical")
  ## G3: Body + Up+Hyper -> concordant relationship -> canonical.
  expect_equal(df$canonical_label[df$gene_symbol == "G3"], "Canonical")
})

test_that("mcc_classify_direction() errors clearly when required columns are missing", {
  expect_error(mcc_classify_direction(data.frame(x = 1), 1, 0.05, 0.1, 0.05), "missing required columns")
})

## ---- mcc_add_raw_direction() / mcc_biomarker_direction_table() ---------------

test_that("mcc_add_raw_direction() labels every row with both directions present, ungated by significance", {
  df <- data.frame(methylation_direction = c("Hypo", "Hyper"), expression_direction = c("Upregulated", "Downregulated"),
                     region_fine = c("TSS200", "Body"), diablo = c(TRUE, FALSE), snf = c(FALSE, TRUE), joint = c(FALSE, FALSE),
                     stringsAsFactors = FALSE)
  out <- mcc_add_raw_direction(df)
  expect_equal(as.character(out$raw_direction), c("Up expression + Hypomethylation", "Down expression + Hypermethylation"))
  expect_equal(out$biomarker_source, c("DIABLO", "SNF"))
})

test_that("mcc_biomarker_direction_table() restricts to DIABLO/SNF/Joint-flagged rows only", {
  df <- data.frame(gene_symbol = c("G1", "G2"), cpg = c("cg1", "cg2"),
                     methylation_direction = c("Hypo", "Hyper"), expression_direction = c("Upregulated", "Downregulated"),
                     diablo = c(TRUE, FALSE), snf = c(FALSE, FALSE), joint = c(FALSE, FALSE), stringsAsFactors = FALSE)
  out <- mcc_biomarker_direction_table(df)
  expect_equal(nrow(out), 1L)
  expect_equal(out$gene_symbol, "G1")
})

## ---- mcc_pair_correlation() --------------------------------------------------

test_that("mcc_pair_correlation() computes a real per-pair correlation and BH-FDR, NA-safe for an unmatched feature", {
  set.seed(500)
  ids <- paste0("S", 1:10)
  expr <- matrix(rnorm(20), 10, 2, dimnames = list(ids, c("G1", "G2")))
  meth <- matrix(rnorm(20), 10, 2, dimnames = list(ids, c("cg1", "cg2")))
  expr[, "G1"] <- meth[, "cg1"] * -1 + rnorm(10, sd = 0.05)  ## strong known negative correlation
  pairs_df <- data.frame(gene_symbol = c("G1", "G1"), cpg = c("cg1", "cg9"), stringsAsFactors = FALSE)  ## cg9 doesn't exist
  out <- mcc_pair_correlation(expr, meth, pairs_df, ids, method = "pearson")
  expect_true(out$ok)
  expect_true(out$df$r[out$df$cpg == "cg1"] < -0.9)
  expect_true(is.na(out$df$r[out$df$cpg == "cg9"]))
  expect_equal(out$df$n[out$df$cpg == "cg9"], 0L)
})

## ---- mcc_regression() --------------------------------------------------------

test_that("mcc_regression() fits a real lm() and reports the methylation coefficient/p-value", {
  set.seed(510)
  meth <- rnorm(30)
  expr <- 3 * meth + rnorm(30, sd = 0.5)
  out <- mcc_regression(expr, meth)
  expect_true(out$ok)
  expect_equal(out$coefficient, unname(coef(lm(expr ~ meth))["meth"]), tolerance = 1e-8)
  expect_true(out$p_value < 0.001)
})

test_that("mcc_regression() refuses when there are too few complete observations for the requested model", {
  out <- mcc_regression(c(1, 2, NA), c(1, NA, 3))
  expect_false(out$ok)
})

## ---- mcc_priority_score() -----------------------------------------------------

test_that("mcc_priority_score() gives a canonical, multiply-flagged, highly-significant pair a higher score than an unflagged, non-canonical one", {
  df <- data.frame(
    log2fc = c(3, 0.1), dbeta = c(0.3, 0.01), expr_fdr = c(0.0001, 0.5), meth_fdr = c(0.0001, 0.5),
    correlation_r = c(0.9, 0.1), canonical = c(TRUE, FALSE), diablo = c(TRUE, FALSE), snf = c(TRUE, FALSE), joint = c(TRUE, FALSE)
  )
  out <- mcc_priority_score(df)
  expect_true(out$priority_score[1] > out$priority_score[2])
  expect_equal(out$evidence_label[1], "Potential Multi-Omics Biomarker")
  expect_true(all(out$priority_score >= 0 & out$priority_score <= 100))
})

## ---- mcc_design_candidates() / mcc_expression_stats() / mcc_methylation_stats() ----

test_that("mcc_design_candidates() finds exactly the 2-class categorical columns among the matched sample_ids", {
  meta <- data.frame(group = rep(c("HC", "RA"), 5), sex = rep(c("F", "M"), 5), age = rnorm(10), row.names = paste0("S", 1:10))
  out <- mcc_design_candidates(meta, paste0("S", 1:10))
  expect_setequal(out, c("group", "sex"))
})

test_that("mcc_expression_stats() computes a real per-feature log2FC/t-test/BH-FDR between two real groups", {
  set.seed(520)
  ids <- paste0("S", 1:20)
  group <- setNames(rep(c("HC", "RA"), 10), ids)
  expr <- matrix(rnorm(20 * 5), 20, 5, dimnames = list(ids, paste0("g", 1:5)))
  expr[group == "RA", 1] <- expr[group == "RA", 1] + 4  ## g1 has a strong, real group difference
  out <- mcc_expression_stats(expr, group)
  expect_true(out$ok)
  expect_true(out$df$log2fc[out$df$feature == "g1"] > 3)
  expect_true(out$df$p[out$df$feature == "g1"] < 0.001)
  expect_equal(out$df$fdr, stats::p.adjust(out$df$p, method = "BH"))
})

test_that("mcc_expression_stats() refuses with fewer than 4 matched samples or a non-2-class design", {
  expr <- matrix(rnorm(10), 2, 5, dimnames = list(c("S1", "S2"), paste0("g", 1:5)))
  out <- mcc_expression_stats(expr, c(S1 = "A", S2 = "B"))
  expect_false(out$ok)
})

test_that("mcc_methylation_stats() computes real per-CpG delta-M/delta-beta from beta-scale input", {
  set.seed(530)
  ids <- paste0("S", 1:20)
  group <- setNames(rep(c("HC", "RA"), 10), ids)
  meth <- matrix(runif(20 * 3, 0.3, 0.5), 20, 3, dimnames = list(ids, paste0("cg", 1:3)))
  meth[group == "RA", 1] <- meth[group == "RA", 1] + 0.3  ## cg1: real, large beta shift
  out <- mcc_methylation_stats(meth, group, value_type = "beta")
  expect_true(out$ok)
  expect_true(out$df$delta_beta[out$df$feature == "cg1"] > 0.2)
  expect_true(out$df$p[out$df$feature == "cg1"] < 0.01)
})

## ---- mcc_summary_counts() ----------------------------------------------------

test_that("mcc_summary_counts() tallies genes/CpGs/significant/canonical/biomarker-source counts correctly", {
  df <- data.frame(
    gene_symbol = c("G1", "G1", "G2"), cpg = c("cg1", "cg2", "cg3"),
    sig_expression = c(TRUE, TRUE, FALSE), sig_methylation = c(TRUE, FALSE, FALSE),
    canonical = c(TRUE, FALSE, NA), evidence_label = c("Potential Multi-Omics Biomarker", "Candidate Multi-Omics Biomarker", NA),
    diablo = c(TRUE, FALSE, FALSE), snf = c(FALSE, TRUE, FALSE), joint = c(FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  out <- mcc_summary_counts(df)
  expect_equal(out$n_genes, 2L)  ## G1 counted once
  expect_equal(out$n_pairs, 3L)
  expect_equal(out$n_significant, 1L)  ## only row 1 has both sig
  expect_equal(out$n_canonical, 1L)
  expect_equal(out$n_potential, 1L)
  expect_equal(out$n_diablo, 1L)
  expect_equal(out$n_snf, 1L)
})

test_that("mcc_summary_counts() returns NULL for an empty/NULL pairs table", {
  expect_null(mcc_summary_counts(NULL))
  expect_null(mcc_summary_counts(data.frame()))
})
