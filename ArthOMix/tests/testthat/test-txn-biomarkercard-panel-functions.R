## Regression guard for the Gene Panel mode added to the Biomarker Card
## (2026-08-26): identifier resolution must never silently drop a submitted
## gene, the per-database evidence-summary status must reflect real

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "functions", "integration", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "transcriptomics", "17_Biomarker_Card", "mod_biomarkercard.R"))

test_that("tbc_panel_identity resolves real gene symbols and reports counts", {
  res <- tbc_panel_identity(c("TNF", "IL6", "STAT3"))
  expect_true(res$ok)
  expect_equal(res$n_submitted, 3L)
  expect_equal(res$n_resolved, 3L)
  expect_equal(res$n_unresolved, 0L)
  expect_true(all(res$df$resolved))
  expect_true(all(grepl("^Resolved", res$df$status_label)))
})

test_that("tbc_panel_identity never silently drops an unrecognized identifier", {
  res <- tbc_panel_identity(c("TNF", "NOTAREALGENEXYZ123"))
  expect_true(res$ok)
  expect_equal(res$n_submitted, 2L)
  expect_equal(res$n_resolved, 1L)
  expect_equal(res$n_unresolved, 1L)
  expect_setequal(res$df$input_id, c("TNF", "NOTAREALGENEXYZ123"))
  bad_row <- res$df[res$df$input_id == "NOTAREALGENEXYZ123", ]
  expect_false(bad_row$resolved)
  expect_equal(bad_row$match_type, "unmatched")
  expect_equal(bad_row$status_label, "Unresolved - identifier not recognized")
})

test_that("tbc_panel_identity resolves duplicate and blank input safely", {
  res <- tbc_panel_identity(c("TNF", "tnf", "  ", "", "TNF"))
  expect_true(res$ok)
  expect_equal(res$n_submitted, 2L)
  expect_equal(res$n_resolved, 2L)
})

test_that("tbc_panel_identity reports failure cleanly on empty input", {
  res <- tbc_panel_identity(character(0))
  expect_false(res$ok)
  expect_match(res$reason, "No gene identifiers")
})

test_that("tbc_evidence_status classifies not-run/failed/empty/found correctly", {
  expect_equal(tbc_evidence_status(NULL, "pathways"), "Not yet run")
  expect_equal(tbc_evidence_status(list(ok = FALSE, reason = "boom"), "pathways"), "Failed")
  expect_equal(tbc_evidence_status(list(ok = TRUE, pathways = data.frame()), "pathways"), "No results")
  expect_equal(tbc_evidence_status(list(ok = TRUE, pathways = data.frame(id = "hsa04010")), "pathways"), "Results found")
})

test_that("tbc_panel_evidence_status mirrors tbc_evidence_status for panel-shaped results", {
  expect_equal(tbc_panel_evidence_status(NULL, "df"), "Not yet run")
  expect_equal(tbc_panel_evidence_status(list(ok = FALSE, error = "no genes"), "df"), "Failed")
  expect_equal(tbc_panel_evidence_status(list(ok = TRUE, df = NULL), "df"), "No results")
  expect_equal(tbc_panel_evidence_status(list(ok = TRUE, df = data.frame(ID = "GO:1")), "df"), "Results found")
})

test_that("tbc_aggregate_convergence groups by the shared item and lists all genes", {
  long_df <- data.frame(
    Gene = c("TNF", "IL6", "TNF", "STAT3"),
    Disease = c("Rheumatoid arthritis", "Rheumatoid arthritis", "Sepsis", "Rheumatoid arthritis"),
    stringsAsFactors = FALSE
  )
  agg <- tbc_aggregate_convergence(long_df, "Disease")
  expect_setequal(agg$Disease, c("Rheumatoid arthritis", "Sepsis"))
  ra_row <- agg[agg$Disease == "Rheumatoid arthritis", ]
  expect_equal(ra_row$`Gene count`, 3L)
  expect_equal(ra_row$Genes, "IL6, STAT3, TNF")
  sepsis_row <- agg[agg$Disease == "Sepsis", ]
  expect_equal(sepsis_row$`Gene count`, 1L)
  expect_equal(agg$Disease[1], "Rheumatoid arthritis")
})

test_that("tbc_aggregate_convergence returns NULL for empty input rather than erroring", {
  expect_null(tbc_aggregate_convergence(NULL, "Disease"))
  expect_null(tbc_aggregate_convergence(data.frame(Gene = character(0), Disease = character(0)), "Disease"))
})

test_that("tbc_split_gene_text tokenizes on comma/newline/tab/space and dedupes", {
  toks <- tbc_split_gene_text("TNF, IL6\nSTAT3\tTNF  FOXP3")
  expect_equal(toks, c("TNF", "IL6", "STAT3", "FOXP3"))
})

test_that("tbc_single_gene_roc reports near-perfect performance for a clearly separable gene", {
  set.seed(1)
  x <- c(rnorm(15, mean = 5, sd = 0.3), rnorm(15, mean = 2, sd = 0.3))
  y <- c(rep("Case", 15), rep("Control", 15))
  res <- tbc_single_gene_roc(x, y, "Case", "Control")
  expect_true(res$ok)
  expect_gt(res$auc, 0.9)
  expect_gt(res$sensitivity, 0.8)
  expect_gt(res$specificity, 0.8)
  expect_equal(res$n_case, 15L)
  expect_equal(res$n_control, 15L)
  expect_equal(sum(res$confusion), 30L)
})

test_that("tbc_single_gene_roc fails honestly (not a fabricated result) with too few samples per group", {
  res <- tbc_single_gene_roc(c(1, 2, 5, 6), c("Case", "Case", "Control", "Control"), "Case", "Control")
  expect_false(res$ok)
  expect_match(res$reason, "Fewer than 3")
})

test_that("tbc_single_gene_cv produces a pooled internal-validation AUC for a separable gene", {
  set.seed(1)
  x <- c(rnorm(20, mean = 5, sd = 0.3), rnorm(20, mean = 2, sd = 0.3))
  y <- c(rep("Case", 20), rep("Control", 20))
  res <- tbc_single_gene_cv(x, y, "Case", "Control", k = 5)
  expect_true(res$ok)
  expect_equal(res$k, 5L)
  expect_gt(res$auc, 0.8)
  expect_true(res$n_used <= 40L)
})

test_that("tbc_single_gene_cv fails honestly with too few samples per group", {
  res <- tbc_single_gene_cv(c(1, 2, 5, 6), c("Case", "Case", "Control", "Control"), "Case", "Control", k = 5)
  expect_false(res$ok)
  expect_match(res$reason, "Not enough samples")
})

test_that("tbc_de_rank computes rank by FDR ordering without fabricating a value for a missing gene", {
  tab <- data.frame(gene = c("A", "B", "C"), adj.P.Val = c(0.2, 0.001, 0.05), stringsAsFactors = FALSE)
  expect_equal(tbc_de_rank("B", tab)$rank, 1L)
  expect_equal(tbc_de_rank("A", tab)$rank, 3L)
  expect_true(is.na(tbc_de_rank("ZZZ", tab)$rank))
  expect_true(is.na(tbc_de_rank("A", NULL)$rank))
})

test_that("tbc_parse_geo_source extracts accession/platform only when the exact pattern is present", {
  res <- tbc_parse_geo_source("NCBI GEO: GSE12345 (Affymetrix HG-U133, collapsed)")
  expect_equal(res$accession, "GSE12345")
  expect_equal(res$platform, "Affymetrix HG-U133")
  res2 <- tbc_parse_geo_source("Example dataset: sex-stratified RA blood cohort")
  expect_true(is.na(res2$accession))
  expect_true(is.na(res2$platform))
})

test_that("tbc_evidence_classification never calls a gene with no significant DE anything but Insufficient evidence", {
  d <- list(live = list(ok = FALSE), dge_hits = NULL, diagnostic_match = list())
  cl <- tbc_evidence_classification(d, ext = NULL, sgd = list(ok = FALSE), sgcv = list(ok = FALSE))
  expect_equal(cl$tier, "Insufficient evidence")
  expect_false(Filter(function(x) x$label == "Significant differential expression", cl$checklist)[[1]]$met)
})

test_that("tbc_evidence_classification stops at Candidate biomarker on DE alone (never calls it validated)", {
  d <- list(live = list(ok = FALSE),
            dge_hits = data.frame(direction = "Up", adj.P.Val = 0.001, stringsAsFactors = FALSE),
            diagnostic_match = list())
  cl <- tbc_evidence_classification(d, ext = NULL, sgd = list(ok = FALSE), sgcv = list(ok = FALSE))
  expect_equal(cl$tier, "Candidate biomarker")
})

test_that("tbc_evidence_classification reaches Supported candidate with diagnostic + internal validation evidence", {
  d <- list(live = list(ok = FALSE),
            dge_hits = data.frame(direction = "Up", adj.P.Val = 0.001, stringsAsFactors = FALSE),
            diagnostic_match = list())
  sgd <- list(ok = TRUE, auc = 0.9)
  sgcv <- list(ok = TRUE, auc = 0.85, k = 5)
  cl <- tbc_evidence_classification(d, ext = NULL, sgd = sgd, sgcv = sgcv)
  expect_equal(cl$tier, "Supported candidate")
})

test_that("tbc_evidence_classification tops out at Supported candidate on training + internal + biological evidence alone - external validation is never inferred from those", {
  d <- list(live = list(ok = FALSE),
            dge_hits = data.frame(direction = "Up", adj.P.Val = 0.001, stringsAsFactors = FALSE),
            diagnostic_match = list())
  sgd <- list(ok = TRUE, auc = 0.95)
  sgcv <- list(ok = TRUE, auc = 0.93, k = 5)
  ext <- list(genetics = list(ok = TRUE, n_diseases = 5), drugs = list(ok = TRUE, drugs = data.frame(Drug = "X")))
  cl <- tbc_evidence_classification(d, ext = ext, sgd = sgd, sgcv = sgcv)
  expect_equal(cl$tier, "Supported candidate")
  expect_false(Filter(function(x) startsWith(x$label, "External validation"), cl$checklist)[[1]]$met)
})

## ---- External-validation evidence tier (audit 2026-09-05): reachable only from a persisted
## frozen-model scoring written by the Diagnostic Model's External Validation tab ----

bmc_min_d <- function(external = NULL) {
  list(
    live = list(ok = FALSE), dge_hits = data.frame(direction = "Up", adj.P.Val = 0.001, stringsAsFactors = FALSE),
    diagnostic_match = list(female = list(in_panel = TRUE, panel_size = 5, n_samples = 100,
                                          lr_auc = 0.9, lr_cv_auc = 0.85, enet_auc = 0.9, enet_cv_auc = 0.86,
                                          rf_auc = 0.95, rf_cv_auc = 0.84, svm_auc = 0.9, svm_cv_auc = 0.83,
                                          genes = c("TNF", "IL6"), external = external))
  )
}

bmc_external_fixture <- function(available = TRUE) {
  list(cohort = "GSE15573 (bundled external blood cohort, PBMC, 33 samples, log2-transformed on load)", cohort_source = "bundled",
       n_ref = 10, n_comp = 14, ref_group = "HC", comp_group = "RA", sex_restricted = TRUE,
       n_genes_present = 2, n_genes_panel = 2, genes = c("TNF", "IL6"), models_scored = available,
       models = list(
         list(key = "lr", label = "Logistic Regression", available = available, auc = if (available) 0.81 else NA_real_,
              ci_lo = 0.7, ci_hi = 0.92, sensitivity = 0.8, specificity = 0.7, reason = if (available) "" else "Could not score"),
         list(key = "rf", label = "Random Forest", available = available, auc = if (available) 0.77 else NA_real_,
              ci_lo = 0.65, ci_hi = 0.89, sensitivity = 0.75, specificity = 0.7, reason = "")))
}

test_that("tbc_external_validation_best()/rows() read the persisted frozen-model scoring and pick the best model", {
  expect_null(tbc_external_validation_rows(NULL))
  expect_null(tbc_external_validation_rows(bmc_min_d()$diagnostic_match))
  expect_true(is.na(tbc_external_validation_best(bmc_min_d()$diagnostic_match)$auc))

  dm <- bmc_min_d(external = bmc_external_fixture())$diagnostic_match
  best <- tbc_external_validation_best(dm)
  expect_equal(best$auc, 0.81)
  expect_match(best$label, "female panel, Logistic Regression, on GSE15573")
  rows <- tbc_external_validation_rows(dm)
  expect_equal(nrow(rows), 2)
  expect_equal(rows$Stratum, c("female", "female"))
  expect_equal(rows$`External AUC (95% CI)`[1], "0.810 (0.700-0.920)")

  dm_unavail <- bmc_min_d(external = bmc_external_fixture(available = FALSE))$diagnostic_match
  expect_true(is.na(tbc_external_validation_best(dm_unavail)$auc))
  expect_null(tbc_external_validation_rows(dm_unavail))
})

test_that("tbc_evidence_classification() reaches 'Strong candidate' only when a frozen-model external scoring is persisted", {
  cl_no_ext <- tbc_evidence_classification(bmc_min_d(), ext = NULL)
  expect_equal(cl_no_ext$tier, "Supported candidate")
  ext_item <- Filter(function(it) grepl("External validation", it$label), cl_no_ext$checklist)[[1]]
  expect_false(ext_item$met)

  cl_ext <- tbc_evidence_classification(bmc_min_d(external = bmc_external_fixture()), ext = NULL)
  expect_equal(cl_ext$tier, "Strong candidate")
  ext_item <- Filter(function(it) grepl("External validation", it$label), cl_ext$checklist)[[1]]
  expect_true(ext_item$met)
  expect_true("Strong candidate" %in% names(TBC_TIER_CLASS))

  cl_unavail <- tbc_evidence_classification(bmc_min_d(external = bmc_external_fixture(available = FALSE)), ext = NULL)
  expect_equal(cl_unavail$tier, "Supported candidate")
})
