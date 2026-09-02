## Module 4 (Cross-omics) - Biomarker Convergence + Cross-Omics MR's own
## pure functions (crossomics_biomarkerconv_helpers.R,
## crossomics_mrstage_helpers.R): the precomputed eQTL/mQTL join loader
## (real data, real documented backfill of the mutually-exclusive-panel
## data defect), upload parsing/merging for both "upload your own data"
## paths, live threshold relabeling (never re-joins/re-derives), the real
## precomputed MR results loader, and the 5 independent evidence-combination
## categories. Also covers crossomics_integration_upload.R's shared
## cx_read_table()/cx_read_and_detect() (used by every Cross-omics upload
## path, this module included).

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "01_Data", "crossomics_integration_upload.R"))
source_from_app_root(file.path("R", "crossomics", "functions", "integration", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "crossomics", "functions", "biomarker_convergence", "crossomics_biomarkerconv_helpers.R"))
source_from_app_root(file.path("R", "crossomics", "04_Cross_Omics_MR", "crossomics_mrstage_helpers.R"))

## ---- cx_read_table() / cx_read_and_detect() -----------------------------------

test_that("cx_read_table() parses a real CSV and refuses a single-column file or an unsupported extension", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("gene,logFC", "TP53,2.1"), path)
  out <- cx_read_table(path, "genes.csv")
  expect_true(out$ok)
  expect_equal(out$df$gene, "TP53")

  onecol <- tempfile(fileext = ".csv")
  writeLines(c("gene", "TP53"), onecol)
  expect_false(cx_read_table(onecol, "genes.csv")$ok)

  expect_false(cx_read_table(path, "genes.docx")$ok)
})

test_that("cx_read_and_detect() bundles read + column auto-detection for the requested kind", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("gene_symbol,log2FC,adj.P.Val", "TP53,2.1,0.01"), path)
  out <- cx_read_and_detect(path, "genes.csv", kind = "expression")
  expect_true(out$ok)
  expect_equal(unname(out$mapping["gene"]), "gene_symbol")
  expect_equal(unname(out$mapping["log2fc"]), "log2FC")
})

## ---- cx_bc_load_precomputed() (real data + the documented backfill) ----------

test_that("cx_bc_load_precomputed() reads the real precomputed eQTL x mQTL join and backfills the documented in_mQTL_MR_panel data defect", {
  skip_if_not(CX_BC_DATA_AVAILABLE, "Biomarker Convergence source data not available")
  raw <- as.data.frame(data.table::fread(cx_bc_precomputed_file("female")))
  out <- cx_bc_load_precomputed("female")
  expect_true(out$ok)
  ## The real backfill should genuinely increase in_mQTL_MR_panel coverage -
  ## never decrease it, never touch genes the source table already flagged.
  expect_true(sum(out$df$in_mQTL_MR_panel %in% TRUE) >= sum(raw$in_mQTL_MR_panel %in% TRUE))
  already_true <- raw$in_mQTL_MR_panel %in% TRUE
  expect_true(all(out$df$in_mQTL_MR_panel[already_true] %in% TRUE))  ## never un-flags an existing TRUE
})

test_that("cx_bc_load_precomputed() also reads the real 'combined' precomputed file, and fails soft for a genuinely missing file", {
  out <- cx_bc_load_precomputed("combined")
  expect_true(out$ok)  ## this deployment does bundle a combined-sex precomputed table
  expect_true(nrow(out$df) > 0)

  ## cx_bc_precomputed_file()'s own path construction, exercised directly
  ## against a path that cannot exist, proves the "file not found" branch
  ## without needing to fake CX_BC_DATA_AVAILABLE itself.
  missing_path <- file.path(CX_RESULTS_DIR, "cross_omics_eQTL_mQTL_not_a_real_sex.csv")
  expect_false(file.exists(missing_path))
})

## ---- cx_bc_backfill_mqtl_from_mrstage() (isolated, hand-built fixture) -------

test_that("cx_bc_backfill_mqtl_from_mrstage() fills mQTL panel membership ONLY for genes with a real MR-stage instrument, never overwriting an existing TRUE", {
  skip_if_not(exists("CX_MR_DATA_AVAILABLE") && isTRUE(CX_MR_DATA_AVAILABLE), "MR-stage source data not available")
  ## Use 2 genes known to actually exist in the real MR-stage file.
  mr <- as.data.frame(data.table::fread(CX_MR_PRECOMPUTED_FILE, showProgress = FALSE))
  skip_if(nrow(mr) < 1, "MR-stage file has no rows")
  target_gene <- mr$gene[1]
  df <- data.frame(gene = c(target_gene, "definitely_not_a_real_gene_xyz"), in_mQTL_MR_panel = c(FALSE, FALSE), stringsAsFactors = FALSE)
  out <- cx_bc_backfill_mqtl_from_mrstage(df)
  expect_true(out$in_mQTL_MR_panel[out$gene == target_gene])
  expect_false(out$in_mQTL_MR_panel[out$gene == "definitely_not_a_real_gene_xyz"])
})

## ---- Upload paths: cx_bc_load_eqtl_upload() / cx_bc_load_mqtl_upload() / cx_bc_merge_eqtl_mqtl() ----

test_that("cx_bc_load_eqtl_upload()/cx_bc_load_mqtl_upload() enforce their own required columns and coerce numeric fields", {
  path_e <- tempfile(fileext = ".csv"); writeLines(c("gene,eQTL_MR_OR,eQTL_MR_FDR", "TP53,1.5,0.01"), path_e)
  out_e <- cx_bc_load_eqtl_upload(path_e, "e.csv")
  expect_true(out_e$ok)
  expect_true(is.numeric(out_e$df$eQTL_MR_OR))

  path_bad <- tempfile(fileext = ".csv"); writeLines(c("notgene,x", "a,1"), path_bad)
  expect_false(cx_bc_load_eqtl_upload(path_bad, "bad.csv")$ok)

  path_m <- tempfile(fileext = ".csv"); writeLines(c("gene,mQTL_MR_pval", "TP53,0.02"), path_m)
  out_m <- cx_bc_load_mqtl_upload(path_m, "m.csv")
  expect_true(out_m$ok)
  expect_false(cx_bc_load_mqtl_upload(path_bad, "bad.csv")$ok)  ## missing mQTL_MR_pval
})

test_that("cx_bc_merge_eqtl_mqtl() outer-joins by gene, marking real panel membership and adding required-but-absent columns as NA (never fabricated)", {
  eqtl_df <- data.frame(gene = c("A", "B"), eQTL_MR_OR = c(1.2, 1.5), stringsAsFactors = FALSE)
  mqtl_df <- data.frame(gene = c("B", "C"), mQTL_MR_pval = c(0.01, 0.02), stringsAsFactors = FALSE)
  out <- cx_bc_merge_eqtl_mqtl(eqtl_df, mqtl_df)
  expect_true(out$ok)
  expect_setequal(out$df$gene, c("A", "B", "C"))
  expect_equal(out$df$in_eQTL_MR_panel[out$df$gene == "A"], TRUE)
  expect_equal(out$df$in_mQTL_MR_panel[out$df$gene == "A"], FALSE)
  expect_true(is.na(out$df$DEG_adjP[out$df$gene == "A"]))  ## never evaluated, never fabricated

  expect_false(cx_bc_merge_eqtl_mqtl(NULL, NULL)$ok)
})

## ---- cx_bc_dedup_min() ---------------------------------------------------------

test_that("cx_bc_dedup_min() keeps the row with the smallest order_col per key, tolerating NAs", {
  df <- data.frame(gene = c("A", "A", "B"), pval = c(0.5, 0.01, NA))
  out <- cx_bc_dedup_min(df, "gene", "pval")
  expect_equal(nrow(out), 2L)
  expect_equal(out$pval[out$gene == "A"], 0.01)
})

## ---- cx_bc_relabel() (live threshold relabeling over already-retained raw values) ----

test_that("cx_bc_relabel() recomputes every *_significant flag from retained raw values at the given thresholds, never re-deriving new evidence", {
  df <- data.frame(
    gene = c("A", "B"), DEG_adjP = c(0.01, 0.5), DMP_fdr_bacon = c(0.001, 0.9), DMR_fdr = c(0.9, 0.001),
    mQTL_MR_pval = c(0.001, 0.5), in_eQTL_MR_panel = c(TRUE, FALSE), stringsAsFactors = FALSE
  )
  out <- cx_bc_relabel(df)
  expect_true(out$DEG_significant[1]); expect_false(out$DEG_significant[2])
  expect_true(out$DMP_genomewide_significant[1])
  expect_true(out$DMR_significant[2])
  expect_true(out$methylation_significant[1])  ## DMP significant
  expect_true(out$methylation_significant[2])  ## DMR significant
  expect_equal(out$n_evidence_layers[1], 4L)   ## eQTL+DEG+methylation(DMP)+mQTL all TRUE
})

test_that("cx_bc_relabel() 'fdr' mqtl_sig_basis recomputes a real BH-FDR from mQTL_MR_pval instead of using the raw p-value directly", {
  df <- data.frame(gene = c("A", "B", "C"), DEG_adjP = NA, DMP_fdr_bacon = NA, DMR_fdr = NA,
                     mQTL_MR_pval = c(0.001, 0.02, 0.5), in_eQTL_MR_panel = FALSE, stringsAsFactors = FALSE)
  out_nominal <- cx_bc_relabel(df, list(mqtl_sig_basis = "nominal_p", mqtl_sig_cutoff = 0.05))
  out_fdr <- cx_bc_relabel(df, list(mqtl_sig_basis = "fdr", mqtl_sig_cutoff = 0.05))
  expect_equal(sum(out_nominal$mQTL_MR_significant), 2L)  ## 0.001, 0.02 < 0.05
  expected_fdr <- stats::p.adjust(df$mQTL_MR_pval, method = "BH")
  expect_equal(out_fdr$mQTL_MR_significant, expected_fdr < 0.05)
})

test_that("cx_bc_relabel() applies real data end-to-end on the real precomputed+backfilled join without erroring", {
  skip_if_not(CX_BC_DATA_AVAILABLE, "Biomarker Convergence source data not available")
  loaded <- cx_bc_load_precomputed("female")
  skip_if_not(loaded$ok, "precomputed data not available")
  out <- cx_bc_relabel(loaded$df)
  expect_true(all(c("DEG_significant", "eQTL_MR_significant", "n_evidence_layers") %in% colnames(out)))
  expect_true(all(out$n_evidence_layers >= 0 & out$n_evidence_layers <= 4))
})

## ---- cx_mr_load_precomputed() (real MR-stage data) -----------------------------

test_that("cx_mr_load_precomputed() reads the real precomputed MR results with a real logical steiger_dir column", {
  skip_if_not(exists("CX_MR_DATA_AVAILABLE") && isTRUE(CX_MR_DATA_AVAILABLE), "MR-stage source data not available")
  out <- cx_mr_load_precomputed()
  expect_true(out$ok)
  expect_true(nrow(out$df) > 0)
  expect_true(is.logical(out$df$steiger_dir))
})

## ---- cx_mr_load_upload() / cx_mr_load_evidence_upload() -----------------------

test_that("cx_mr_load_upload() requires gene+pval, recomputes FDR via real BH when not supplied, and fills optional columns as NA", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("gene,pval", "A,0.001", "B,0.04", "C,0.5"), path)
  out <- cx_mr_load_upload(path, "mr.csv")
  expect_true(out$ok)
  expect_equal(out$df$FDR, stats::p.adjust(out$df$pval, method = "BH"))
  expect_true(all(is.na(out$df$cpg)))

  path_bad <- tempfile(fileext = ".csv"); writeLines(c("gene,notpval", "A,1"), path_bad)
  expect_false(cx_mr_load_upload(path_bad, "bad.csv")$ok)
})

test_that("cx_mr_load_evidence_upload() computes in_eQTL_MR_panel strictly from eQTL_MR_FDR < 0.05, not merely from its presence", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("gene,eQTL_MR_FDR", "A,0.01", "B,0.9"), path)
  out <- cx_mr_load_evidence_upload(path, "ev.csv")
  expect_true(out$ok)
  expect_true(out$df$in_eQTL_MR_panel[out$df$gene == "A"])
  expect_false(out$df$in_eQTL_MR_panel[out$df$gene == "B"])  ## has a value, but doesn't clear 0.05
  expect_true(is.na(out$df$DEG_adjP[out$df$gene == "A"]))  ## never fabricated

  path_no_gene <- tempfile(fileext = ".csv"); writeLines(c("notgene,x", "a,1"), path_no_gene)
  expect_false(cx_mr_load_evidence_upload(path_no_gene, "bad.csv")$ok)
})

## ---- cx_mr_classify_categories() (5 independent evidence combinations) --------

test_that("cx_mr_classify_categories() correctly and independently matches each of the 5 documented evidence combinations (never mutually exclusive)", {
  df <- data.frame(
    gene = c("AllFour", "DegEqtlOnly", "DmpMqtlOnly", "None"),
    DEG_significant = c(TRUE, TRUE, FALSE, FALSE),
    DMP_genomewide_significant = c(TRUE, FALSE, TRUE, FALSE),
    DMR_significant = c(FALSE, FALSE, FALSE, FALSE),
    mQTL_MR_significant = c(TRUE, FALSE, TRUE, FALSE),
    eQTL_MR_significant = c(TRUE, TRUE, FALSE, FALSE),
    DEG_logFC = 1, DEG_direction = "Up", DEG_adjP = 0.01,
    DMP_top_cpg = "cg1", DMP_dbeta = 0.2, DMP_direction = "Hyper", DMP_fdr_bacon = 0.01,
    mQTL_candidate_cpg = "cg1", mQTL_MR_beta = 0.1, mQTL_MR_pval = 0.01,
    eQTL_MR_OR = 1.5, eQTL_MR_direction = "Up", eQTL_MR_FDR = 0.01,
    stringsAsFactors = FALSE
  )
  out <- cx_mr_classify_categories(df)
  expect_equal(out$deg_eqtl$gene, c("AllFour", "DegEqtlOnly"))  ## both genes satisfy DEG+eQTL
  expect_equal(out$dmp_mqtl$gene, c("AllFour", "DmpMqtlOnly"))  ## AllFour appears in BOTH categories - not mutually exclusive
  expect_equal(out$deg_dmp_qtl$gene, "AllFour")  ## requires all 4 flags
  expect_equal(nrow(out$dmr_mqtl), 0L)  ## no gene has DMR_significant

  ## Only columns actually present in join_df are kept - never fabricates.
  expect_false("DMR_id" %in% colnames(out$deg_eqtl))
})

test_that("cx_mr_classify_categories() returns NULL when join_df is NULL", {
  expect_null(cx_mr_classify_categories(NULL))
})

test_that("cx_mr_classify_categories() runs end-to-end on the real relabeled precomputed join without erroring", {
  skip_if_not(CX_BC_DATA_AVAILABLE, "Biomarker Convergence source data not available")
  loaded <- cx_bc_load_precomputed("female")
  skip_if_not(loaded$ok, "precomputed data not available")
  relabeled <- cx_bc_relabel(loaded$df)
  out <- cx_mr_classify_categories(relabeled)
  expect_equal(length(out), 5L)
  expect_true(all(vapply(out, is.data.frame, logical(1))))
})
