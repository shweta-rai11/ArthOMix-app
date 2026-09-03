## Regression guard for the Methylomics Biomarker Card's API-driven Evidence
## Explorer extension (2026-08-26): identifier detection/resolution must
## never silently drop a submitted gene or CpG, the source/API provenance

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "functions", "annotation.R"))
source_from_app_root(file.path("R", "crossomics", "functions", "integration", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "methylomics", "15_Biomarker_Analysis", "mod_methyl_biomarkercard.R"))

test_that("bc_detect_identifier_type routes CpG/Ensembl/Entrez/gene-symbol tokens correctly", {
  expect_equal(bc_detect_identifier_type("cg00000029"), "cpg")
  expect_equal(bc_detect_identifier_type("CG00000029"), "cpg")
  expect_equal(bc_detect_identifier_type("ENSG00000012048"), "ensembl")
  expect_equal(bc_detect_identifier_type("ENSG00000012048.15"), "ensembl")
  expect_equal(bc_detect_identifier_type("672"), "entrez")
  expect_equal(bc_detect_identifier_type("BRCA1"), "gene_symbol")
  expect_equal(bc_detect_identifier_type("NOTAREALGENEXYZ123"), "gene_symbol")
})

test_that("bc_split_tokens tokenizes on comma/newline/tab/space and dedupes", {
  toks <- bc_split_tokens("BRCA1, cg00000029\nTP53\tBRCA1  cg00000029")
  expect_equal(toks, c("BRCA1", "cg00000029", "TP53"))
})

test_that("bc_split_tokens returns an empty character vector for blank/NULL input", {
  expect_equal(bc_split_tokens(""), character(0))
  expect_equal(bc_split_tokens(NULL), character(0))
  expect_equal(bc_split_tokens("   "), character(0))
})

test_that("bc_resolve_identifiers resolves a mixed gene+CpG list and never silently drops an entry", {
  res <- bc_resolve_identifiers(c("BRCA1", "cg00000029", "NOTAREALGENEXYZ123"), array_type = "450K")
  expect_true(res$ok)
  expect_equal(res$n_submitted, 3L)
  expect_setequal(res$df$input_id, c("BRCA1", "cg00000029", "NOTAREALGENEXYZ123"))
  expect_true(res$df$resolved[res$df$input_id == "BRCA1"])
  expect_true(res$df$resolved[res$df$input_id == "cg00000029"])
  expect_false(res$df$resolved[res$df$input_id == "NOTAREALGENEXYZ123"])
  expect_equal(res$df$detected_type[res$df$input_id == "cg00000029"], "CpG probe")
  expect_equal(res$df$detected_type[res$df$input_id == "BRCA1"], "Gene identifier")
})

test_that("bc_resolve_identifiers handles duplicate and blank input safely", {
  res <- bc_resolve_identifiers(c("BRCA1", "brca1", "  ", "", "BRCA1"))
  expect_true(res$ok)
  expect_equal(res$n_submitted, 2L)
  expect_equal(res$n_resolved, 2L)
})

test_that("bc_resolve_identifiers reports failure cleanly on empty input", {
  res <- bc_resolve_identifiers(character(0))
  expect_false(res$ok)
  expect_match(res$reason, "No identifiers")
  expect_null(res$df)
})

test_that("an invalid CpG probe ID is reported unresolved, not dropped", {
  res <- bc_resolve_identifiers("cg99999999999", array_type = "450K")
  expect_true(res$ok)
  expect_equal(res$n_resolved, 0L)
  expect_false(res$df$resolved[1])
  expect_match(res$df$status_label[1], "Unresolved")
})

test_that("bc_meta builds a complete provenance envelope with a live timestamp", {
  before <- Sys.time()
  m <- bc_meta("Test DB", "BRCA1", "example.org/api", "Success", 5L)
  after <- Sys.time()
  expect_equal(m$source, "Test DB")
  expect_equal(m$query, "BRCA1")
  expect_equal(m$endpoint, "example.org/api")
  expect_equal(m$status, "Success")
  expect_equal(m$n_records, 5L)
  expect_true(m$retrieved_at >= before && m$retrieved_at <= after)
})

test_that("bc_evidence_status classifies not-run/failed/no-results/found correctly", {
  expect_equal(bc_evidence_status(NULL, "pathways"), "Not yet run")
  expect_equal(bc_evidence_status(list(ok = FALSE, reason = "boom"), "pathways"), "Failed")
  expect_equal(bc_evidence_status(list(ok = TRUE, pathways = data.frame()), "pathways"), "No results")
  expect_equal(bc_evidence_status(list(ok = TRUE, pathways = data.frame(id = "hsa04010")), "pathways"), "Results found")
  expect_equal(bc_evidence_status(list(ok = TRUE, name = "BRCA1"), "name"), "Results found")
})

test_that("BC_EVIDENCE_DBS covers every new live-API client with a matching key", {
  keys <- vapply(BC_EVIDENCE_DBS, `[[`, character(1), "key")
  expect_true(all(c("ncbi_gene", "ensembl", "regulatory", "genetics", "gwas_catalog", "hpa", "encode", "geo", "biostudies", "literature") %in% keys))
})

test_that("bc_aggregate_convergence groups by the shared item and lists all biomarkers", {
  long_df <- data.frame(
    ID = c("BRCA1", "TP53", "BRCA1", "STAT3"),
    Disease = c("Breast cancer", "Breast cancer", "Ovarian cancer", "Breast cancer"),
    stringsAsFactors = FALSE
  )
  agg <- bc_aggregate_convergence(long_df, "Disease")
  expect_setequal(agg$Disease, c("Breast cancer", "Ovarian cancer"))
  bc_row <- agg[agg$Disease == "Breast cancer", ]
  expect_equal(bc_row$`Biomarker count`, 3L)
  expect_equal(bc_row$Biomarkers, "BRCA1, STAT3, TP53")
  expect_equal(agg$Disease[1], "Breast cancer")
})

test_that("bc_aggregate_convergence returns NULL for empty input rather than erroring", {
  expect_null(bc_aggregate_convergence(NULL, "Disease"))
  expect_null(bc_aggregate_convergence(data.frame(ID = character(0), Disease = character(0)), "Disease"))
})

test_that("bc_literature_classify tags plausible categories from title text", {
  cls <- bc_literature_classify("An epigenome-wide association study (EWAS) of BRCA1 promoter methylation as a biomarker for breast cancer risk")
  expect_true("EWAS" %in% cls)
  expect_true("DNA methylation" %in% cls)
  expect_true("Biomarker study" %in% cls)
  expect_true("Disease association" %in% cls)
})

test_that("bc_literature_classify returns an empty vector for blank input, never errors", {
  expect_equal(bc_literature_classify(""), character(0))
  expect_equal(bc_literature_classify(NULL), character(0))
})

test_that("bc_literature_classify flags a review article", {
  cls <- bc_literature_classify("Systematic review of DNA methylation biomarkers in autoimmune disease")
  expect_true("Review" %in% cls)
})

test_that("bc_literature_query builds gene/CpG-flavored preset queries", {
  expect_equal(bc_literature_query("BRCA1", "Gene + methylation"), "BRCA1 methylation")
  expect_equal(bc_literature_query("BRCA1", "Gene + CpG"), "BRCA1 CpG methylation")
  expect_equal(bc_literature_query("BRCA1", "disease", "rheumatoid arthritis"), "BRCA1 rheumatoid arthritis")
  expect_equal(bc_literature_query("BRCA1", "disease"), "BRCA1 disease")
})

test_that("bc_panel_gene_rows returns one row per gene, NA-filled for an unresolvable symbol", {
  rows <- bc_panel_gene_rows(c("BRCA1", "NOTAREALGENEXYZ123"))
  expect_equal(nrow(rows), 2L)
  expect_true(all(c("Gene", "NCBI Entrez ID", "Ensembl Gene ID", "Chromosome") %in% colnames(rows)))
  expect_true(is.na(rows$`NCBI Entrez ID`[rows$Gene == "NOTAREALGENEXYZ123"]))
})

test_that("bc_panel_gene_rows / bc_panel_cpg_rows return NULL for empty input", {
  expect_null(bc_panel_gene_rows(character(0)))
  expect_null(bc_panel_cpg_rows(character(0), "450K"))
})
