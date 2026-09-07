## Regression coverage for the data-governance guarantee ArthOChat's context
## builder makes: no sample-level identifiers or raw data values are ever
## transmitted to the LLM backend, at any nesting depth of a stored result.
##
## Background (2026-09-06 audit): .format_results_block() previously checked
## is.matrix()/is.data.frame() only on a field's *immediate* value. Any
## data.frame, matrix, or sample-ID vector wrapped one level inside list() -
## a natural R idiom used throughout the submodules - fell through to
## as.character() and was deparsed verbatim with no size cap. Three real
## leaks existed: multi_results$integration$snf (raw per-sample metadata +
## raw omics matrices + clustering object), multi_results$overview$
## harmonization$matched_ids (matched sample IDs), and results$
## featureselection$<sex>$holdout_sample_ids (held-out sample IDs). Matrix
## column names were a fourth vector: genes x samples and sample x sample
## matrices carry sample IDs in their dimnames.
##
## These tests exercise the pure formatter with synthetic data shaped like
## those real result objects. No live LLM is needed.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
all_r_files <- list.files(file.path(app_dir, "R"), pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
all_r_files <- all_r_files[!grepl("modules_index[.]R$", all_r_files)]
for (f in all_r_files) suppressWarnings(suppressMessages(source(f)))
source_from_app_root(file.path("R", "modules_index.R"))

sample_ids <- sprintf("GSM%07d", 1:30)
emits <- function(lines, needle) any(grepl(needle, lines, fixed = TRUE))

## ---------------------------------------------------------------------------
## Nested data.frames / matrices / model objects (the SNF-bundle shape)
## ---------------------------------------------------------------------------

test_that("a nested list of raw matrices, per-sample metadata and a clustering object emits no sample IDs or values", {
  set.seed(1)
  expr <- matrix(rnorm(30 * 5), 30, 5, dimnames = list(sample_ids, paste0("g", 1:5)))
  meta <- data.frame(sex = rep(c("F", "M"), 15), age = 40:69, row.names = sample_ids)
  snf_res <- structure(list(clusters = setNames(sample(1:3, 30, TRUE), sample_ids),
                            W = expr %*% t(expr)), class = "snf_fit")
  res <- list(n_clusters = 3L,
              snf = list(sample_meta = meta, layers = list(mrna = expr, meth = expr), result = snf_res))
  out <- .format_results_block("Integration", res)

  expect_false(emits(out, "GSM"))
  expect_false(emits(out, sprintf("%.4f", expr[1, 1])))
  expect_true(emits(out, "snf.sample_meta: 30 x 2 table (columns: sex, age)"))
  expect_true(emits(out, "snf.layers.mrna: 30 x 5 matrix"))
  expect_true(emits(out, "snf.result.W: 30 x 30 matrix"))
  expect_true(emits(out, "- n_clusters: 3"))
})

test_that("matrix dimnames are never emitted, even for a top-level matrix (genes x samples orientation)", {
  expr <- matrix(0, 4, 3, dimnames = list(paste0("gene", 1:4), sample_ids[1:3]))
  out <- .format_results_block("X", list(expr = expr))
  expect_true(emits(out, "- expr: 4 x 3 matrix"))
  expect_false(emits(out, "GSM"))
  expect_false(emits(out, "gene1"))
})

test_that("data.frame column names (schema) are still emitted, capped at 10", {
  df <- as.data.frame(setNames(replicate(12, 1:2, simplify = FALSE), paste0("col", 1:12)))
  out <- .format_results_block("X", list(tbl = df))
  expect_true(emits(out, "- tbl: 2 x 12 table (columns: col1, col2, col3, col4, col5, col6, col7, col8, col9, col10)"))
  expect_false(emits(out, "col11"))
})

## ---------------------------------------------------------------------------
## Sample-identifier vectors are withheld by field name, at any depth
## ---------------------------------------------------------------------------

test_that("matched_ids at any depth is withheld and only counted", {
  res <- list(harmonization = list(n_matched = 30L, matched_ids = sample_ids))
  out <- .format_results_block("Overview", res)
  expect_false(emits(out, "GSM"))
  expect_true(emits(out, "harmonization.matched_ids: 30 identifiers (withheld)"))
  expect_true(emits(out, "harmonization.n_matched: 30"))
})

test_that("holdout_sample_ids nested per sex is withheld while the feature-level gene panel is still grounded", {
  res <- list(female = list(holdout_sample_ids = sample_ids[1:6],
                            consensus_genes = c("TNF", "IL6", "STAT1"), auc = 0.81))
  out <- .format_results_block("Feature Selection", res)
  expect_false(emits(out, "GSM"))
  expect_true(emits(out, "female.holdout_sample_ids: 6 identifiers (withheld)"))
  expect_true(emits(out, "female.consensus_genes: TNF, IL6, STAT1"))
  expect_true(emits(out, "female.auc: 0.81"))
})

test_that("the identifier-name pattern covers the identifier field names used across the codebase", {
  for (nm in c("sample_ids", "sampleIds", "patient_id", "subject_ids", "train_sample_ids",
               "holdout_sample_ids", "matched_ids", "reserved_ids", "ids", "rownames")) {
    out <- .format_results_block("X", setNames(list(sample_ids), nm))
    expect_false(emits(out, "GSM"), info = nm)
    expect_true(emits(out, "identifiers (withheld)"), info = nm)
  }
})

test_that("feature-level identifier fields are NOT withheld (they are the intended grounding signal)", {
  for (nm in c("top_hits", "genes", "consensus_genes", "selected_cpgs", "biomarker_genes", "unmapped_ids_note")) {
    out <- .format_results_block("X", setNames(list(c("TNF", "IL6")), nm))
    expect_true(emits(out, "TNF, IL6"), info = nm)
  }
})

## ---------------------------------------------------------------------------
## Size and depth caps
## ---------------------------------------------------------------------------

test_that("a long feature vector is capped at 20 values with a total count", {
  out <- .format_results_block("X", list(genes = paste0("G", 1:100)))
  expect_true(emits(out, "G20 ... (100 total)"))
  expect_false(emits(out, "G21"))
})

test_that("names on an atomic vector are dropped before emission", {
  clusters <- setNames(c(1L, 2L, 1L), sample_ids[1:3])
  out <- .format_results_block("X", list(clusters = clusters))
  expect_true(emits(out, "- clusters: 1, 2, 1"))
  expect_false(emits(out, "GSM"))
})

test_that("nesting beyond the depth cap collapses to an element count rather than expanding", {
  deep <- list(a = list(b = list(c = list(d = list(e = sample_ids)))))
  out <- .format_results_block("Deep", deep)
  expect_false(emits(out, "GSM"))
  expect_true(emits(out, "- a.b.c: nested structure with 1 elements (not expanded)"))
})

test_that("non-list, non-atomic objects show their class only", {
  out <- .format_results_block("X", list(fit = new.env()))
  expect_true(emits(out, "- fit: <environment object, not shown>"))
})

test_that("a result block is capped in total line count", {
  wide <- setNames(as.list(seq_len(200)), paste0("k", seq_len(200)))
  out <- .format_results_block("Wide", wide)
  expect_lte(length(out), 1L + .ARTHOCHAT_MAX_LINES_PER_BLOCK)
})

## ---------------------------------------------------------------------------
## Pre-existing behaviour preserved
## ---------------------------------------------------------------------------

test_that("the NOT YET RUN block for a NULL result is unchanged", {
  out <- .format_results_block("Differential Expression", NULL)
  expect_true(emits(out, "### Differential Expression"))
  expect_true(emits(out, "NOT YET RUN IN THIS SESSION"))
})

test_that("top-level scalars, top_hits and data.frames render exactly as before", {
  res <- list(n_sig = 142L, top_hits = c("TNF", "IL6"), table = data.frame(gene = "TNF", p = 0.01))
  out <- .format_results_block("DGE", res)
  expect_identical(out[[1]], "### DGE")
  expect_true(emits(out, "- n_sig: 142"))
  expect_true(emits(out, "- top_hits: TNF, IL6"))
  expect_true(emits(out, "- table: 1 x 2 table (columns: gene, p)"))
})
