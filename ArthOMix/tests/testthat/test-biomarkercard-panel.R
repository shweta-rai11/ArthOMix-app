## Regression guard for the Gene Panel mode added to the Biomarker Card
## (2026-08-26): identifier resolution must never silently drop a submitted
## gene, the per-database evidence-summary status must reflect real
## retrieval outcomes (not just "the database exists"), and the disease/drug
## convergence aggregator must group by the shared item (disease/drug) and
## list every gene that hit it - all pure logic, tested here without any
## live network call.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "transcriptomics", "mod_biomarkercard.R"))

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
  ## Both inputs must still appear in the resolution table, in original form.
  expect_setequal(res$df$input_id, c("TNF", "NOTAREALGENEXYZ123"))
  bad_row <- res$df[res$df$input_id == "NOTAREALGENEXYZ123", ]
  expect_false(bad_row$resolved)
  expect_equal(bad_row$match_type, "unmatched")
  expect_equal(bad_row$status_label, "Unresolved - identifier not recognized")
})

test_that("tbc_panel_identity resolves duplicate and blank input safely", {
  res <- tbc_panel_identity(c("TNF", "tnf", "  ", "", "TNF"))
  expect_true(res$ok)
  ## Blanks dropped before harmonization; case-insensitive dup collapses to one row.
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
  ## Sorted by descending gene count - the most-convergent item comes first.
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

## ---- Single-gene diagnostic performance (2026-08-27 evidence-dossier -------
## refactor): pure local computation on synthetic expression vectors, no
## network call, no dependency on mod_diagnostic.R's private reactives.

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

## ---- Evidence tier classification (spec: DE alone is never "validated") ----

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

test_that("tbc_evidence_classification never reaches Strong candidate - external validation is never persisted in this deployment", {
  d <- list(live = list(ok = FALSE),
            dge_hits = data.frame(direction = "Up", adj.P.Val = 0.001, stringsAsFactors = FALSE),
            diagnostic_match = list())
  sgd <- list(ok = TRUE, auc = 0.95)
  sgcv <- list(ok = TRUE, auc = 0.93, k = 5)
  ext <- list(genetics = list(ok = TRUE, n_diseases = 5), drugs = list(ok = TRUE, drugs = data.frame(Drug = "X")))
  cl <- tbc_evidence_classification(d, ext = ext, sgd = sgd, sgcv = sgcv)
  expect_equal(cl$tier, "Supported candidate")
  expect_false(Filter(function(x) x$label == "External validation", cl$checklist)[[1]]$met)
})
