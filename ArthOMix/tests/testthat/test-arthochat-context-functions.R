## Regression coverage for ArthOChat fixes made in the forensic audit
## (2026-09-04): (B) build_cx_context()'s dataset-scope line, which was
## previously entirely absent for the Cross-Omics vertical, causing the
## model to fabricate cohort claims for a sub-module it had no context for
## (see tests/arthochat_verification/README.md); (C) a new, code-level
## fabrication guard that doesn't require a live Ollama server to test;
## (D) case-insensitive contrast/group matching in the DGE tool, closing an
## asymmetry with the other 3 ArthOChat-invocable analysis tools. None of
## these tests require a live LLM - they exercise pure, LLM-independent
## helper functions with synthetic strings/data.

## build_cx_context()/CX_MODULES (R/modules_index.R) construct every
## vertical's module registry at source-time from each individual module
## file's own *_config object (e.g. mod_overview_config, mod_cross_*_config)
## - exactly like the real app's own R/ autoload does - so modules_index.R
## cannot be sourced standalone without every module file already loaded
## first. Mirror that by sourcing the whole R/ tree (order-independent for
## function *definitions*; only modules_index.R's own top-level registry
## construction needs everything else to exist first) before modules_index.R
## and mod_arthochat.R themselves.
suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
all_r_files <- list.files(file.path(app_dir, "R"), pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
all_r_files <- all_r_files[!grepl("modules_index[.]R$", all_r_files)]
for (f in all_r_files) suppressWarnings(suppressMessages(source(f)))
source_from_app_root(file.path("R", "modules_index.R"))
source_from_app_root(file.path("R", "shared", "mod_arthochat.R"))

## ---------------------------------------------------------------------------
## (B) build_cx_context() dataset-scope line regression tests
## ---------------------------------------------------------------------------

test_that("build_cx_context() emits a dataset-scope line for 'integration' even with entirely empty results", {
  ctx <- build_cx_context(list(), list(), focus_id = "integration")
  expect_match(ctx, "Cross-Omics: currently loaded dataset", fixed = TRUE)
  expect_match(ctx, "NOT YET LOADED IN THIS SESSION", fixed = TRUE)
})

test_that("build_cx_context() emits a dataset-scope line for 'biomarkerconv' even with entirely empty results", {
  ctx <- build_cx_context(list(), list(), focus_id = "biomarkerconv")
  expect_match(ctx, "Cross-Omics: currently loaded dataset", fixed = TRUE)
  expect_match(ctx, "NOT YET LOADED IN THIS SESSION", fixed = TRUE)
})

test_that("build_cx_context() emits a dataset-scope line for 'mrstage' even with entirely empty results", {
  ctx <- build_cx_context(list(), list(), focus_id = "mrstage")
  expect_match(ctx, "Cross-Omics: currently loaded dataset", fixed = TRUE)
  expect_match(ctx, "NOT YET LOADED IN THIS SESSION", fixed = TRUE)
})

test_that("build_cx_context() reports the real loaded sources once a live Cross-Omics dataset exists", {
  cross_dataset <- list(user_expr_source = "Live Transcriptomics DGE: HC vs RA",
                         user_meth_source = "Live Methylomics DMP: HC vs RA")
  ctx <- build_cx_context(cross_dataset, list(), focus_id = "biomarkerconv")
  expect_match(ctx, "Live Transcriptomics DGE: HC vs RA", fixed = TRUE)
  expect_match(ctx, "Live Methylomics DMP: HC vs RA", fixed = TRUE)
  expect_false(grepl("NOT YET LOADED IN THIS SESSION", ctx, fixed = TRUE))
})

## ---------------------------------------------------------------------------
## (C) arthochat_detect_ungrounded_reference() - the fabrication guard
## ---------------------------------------------------------------------------

## A real module title, not a hand-invented "WGCNA" label - every actual
## WGCNA title (R/modules_index.R) carries a trailing/parenthetical
## qualifier, so a bare "WGCNA" header never occurs in production context
## text and a test built on one would silently stop exercising the guard.
WGCNA_TITLE <- Find(function(m) grepl("^WGCNA", m), .arthochat_known_modules())

mk_context <- function(populated_module = NULL, not_run_module = NULL) {
  lines <- character(0)
  if (!is.null(populated_module)) {
    lines <- c(lines, sprintf("## %s", populated_module),
               "n_significant: 42", "top_hits: GENE1, GENE2, GENE3", "")
  }
  if (!is.null(not_run_module)) {
    lines <- c(lines, sprintf("## %s", not_run_module),
               "NOT YET RUN IN THIS SESSION.", "")
  }
  paste(lines, collapse = "\n")
}

test_that("a response mentioning a module WITH real populated context is not flagged", {
  ctx <- mk_context(populated_module = "Differential Expression")
  resp <- "Differential Expression found 42 significant genes, including GENE1 and GENE2."
  chk <- arthochat_detect_ungrounded_reference(resp, ctx)
  expect_false(chk$flagged)
  expect_length(chk$modules, 0)
})

test_that("a response mentioning a module ABSENT from context (marked not-yet-run) IS flagged", {
  ctx <- mk_context(not_run_module = WGCNA_TITLE)
  resp <- "WGCNA identified 6 co-expression modules, with the turquoise module most strongly associated with disease status."
  chk <- arthochat_detect_ungrounded_reference(resp, ctx)
  expect_true(chk$flagged)
  expect_true(WGCNA_TITLE %in% chk$modules)
})

test_that("a response that mentions no specific known sub-module at all is not flagged", {
  ctx <- mk_context(not_run_module = WGCNA_TITLE)
  resp <- "Rheumatoid arthritis is an autoimmune disease that primarily affects the joints."
  chk <- arthochat_detect_ungrounded_reference(resp, ctx)
  expect_false(chk$flagged)
})

test_that("a response that hedges appropriately when discussing an unrun module is not flagged", {
  ctx <- mk_context(not_run_module = WGCNA_TITLE)
  resp <- "WGCNA hasn't been run in this session, so I can't tell you about any co-expression modules yet."
  chk <- arthochat_detect_ungrounded_reference(resp, ctx)
  expect_false(chk$flagged)
})

test_that("a response naming the correct, genuinely populated module is not falsely flagged, even alongside an unrelated unrun module", {
  ctx <- mk_context(populated_module = "Differential Expression", not_run_module = WGCNA_TITLE)
  resp <- "Differential Expression found 42 significant genes, including GENE1 and GENE2."
  chk <- arthochat_detect_ungrounded_reference(resp, ctx)
  expect_false(chk$flagged)
  expect_false(WGCNA_TITLE %in% chk$modules)
})

test_that(".arthochat_classify_context_modules() resolves a header to its most specific matching title, not a shorter title it happens to contain", {
  ## "ML Feature Selection" (Methylomics) contains "Feature Selection"
  ## (Transcriptomics) as a substring - the classifier must not misresolve
  ## the header to the shorter, wrong-module title.
  ml_title <- Find(function(m) grepl("^ML Feature Selection", m), .arthochat_known_modules())
  ctx <- mk_context(not_run_module = ml_title)
  cls <- .arthochat_classify_context_modules(ctx)
  expect_true(ml_title %in% cls$not_run)
  expect_false("Feature Selection" %in% cls$not_run)
})

test_that("a hedge about one unrun module does not suppress flagging of a DIFFERENT unrun module in the same response", {
  ## Regression guard for the 2026-09-07 defense audit finding: hedge
  ## suppression was response-wide, so hedging about module A disabled
  ## ungrounded-reference flagging for every other not-run module named in
  ## that same response, not just A.
  ml_title <- Find(function(m) grepl("^ML Feature Selection", m), .arthochat_known_modules())
  ctx <- mk_context(not_run_module = WGCNA_TITLE)
  ctx <- paste(ctx, sprintf("## %s\nNOT YET RUN IN THIS SESSION.\n", ml_title), sep = "\n")
  resp <- paste0(
    ml_title, " hasn't been run in this session yet. ",
    "WGCNA identified 6 co-expression modules, with the turquoise module most strongly associated with disease status."
  )
  chk <- arthochat_detect_ungrounded_reference(resp, ctx)
  expect_true(chk$flagged)
  expect_true(WGCNA_TITLE %in% chk$modules)
  expect_false(ml_title %in% chk$modules)
})

test_that("an unrelated word containing an unrun module's acronym as a substring does not false-flag it", {
  ## Regression guard: the acronym matcher used fixed substring matching, so
  ## "ML" (from "ML Feature Selection") matched any response containing
  ## "html", "normally", "warmly", etc.
  ml_title <- Find(function(m) grepl("^ML Feature Selection", m), .arthochat_known_modules())
  ctx <- mk_context(not_run_module = ml_title)
  resp <- "You can export this table as html or csv, and it renders normally in most browsers."
  chk <- arthochat_detect_ungrounded_reference(resp, ctx)
  expect_false(chk$flagged)
  expect_false(ml_title %in% chk$modules)
})

test_that("the acronym itself, on a real word boundary, still correctly flags the unrun module", {
  ml_title <- Find(function(m) grepl("^ML Feature Selection", m), .arthochat_known_modules())
  ctx <- mk_context(not_run_module = ml_title)
  resp <- "ML found the top 20 CpGs most predictive of RA status."
  chk <- arthochat_detect_ungrounded_reference(resp, ctx)
  expect_true(chk$flagged)
  expect_true(ml_title %in% chk$modules)
})

## ---------------------------------------------------------------------------
## (C2) arthochat_affirms_pending_run() - code-level consent check for the
## DGE-execution tool (2026-09-07 defense audit finding: consent for the one
## write-capable ArthoChat tool was entirely model-judged, never checked
## against the user's actual last message in code).
## ---------------------------------------------------------------------------

test_that("clear affirmative replies are recognized as consent", {
  for (msg in c("Yes", "yes, go ahead", "sure, run it", "ok do it", "Yeah please proceed",
                "sounds good", "Let's do it", "please run it", "confirmed")) {
    expect_true(arthochat_affirms_pending_run(msg), info = msg)
  }
})

test_that("declines and off-topic messages are NOT recognized as consent", {
  for (msg in c("no", "no thanks", "don't run it", "wait, not yet", "cancel that",
                "actually what does WGCNA do?", "", "   ")) {
    expect_false(arthochat_affirms_pending_run(msg), info = msg)
  }
})

test_that("NULL input is not treated as consent", {
  expect_false(arthochat_affirms_pending_run(NULL))
})

test_that("a vague reply that isn't a clear yes/no is not treated as consent", {
  expect_false(arthochat_affirms_pending_run("hmm, tell me more about the covariate option first"))
})

test_that("an empty or NULL response is never flagged", {
  ctx <- mk_context(not_run_module = WGCNA_TITLE)
  expect_false(arthochat_detect_ungrounded_reference("", ctx)$flagged)
  expect_false(arthochat_detect_ungrounded_reference(NULL, ctx)$flagged)
  expect_false(arthochat_detect_ungrounded_reference("   ", ctx)$flagged)
})

## ---------------------------------------------------------------------------
## arthochat_grounded_modules_label() - the transparency footer helper
## ---------------------------------------------------------------------------

test_that("arthochat_grounded_modules_label() lists only modules with real (not-yet-run-free) context", {
  ctx <- mk_context(populated_module = "Differential Expression", not_run_module = WGCNA_TITLE)
  lbl <- arthochat_grounded_modules_label(ctx)
  expect_match(lbl, "Differential Expression", fixed = TRUE)
  expect_false(grepl(WGCNA_TITLE, lbl, fixed = TRUE))
})

test_that("arthochat_grounded_modules_label() returns an empty string when nothing is grounded", {
  ctx <- mk_context(not_run_module = WGCNA_TITLE)
  expect_equal(arthochat_grounded_modules_label(ctx), "")
})

## ---------------------------------------------------------------------------
## End-to-end: the guard must fire against REAL context text, not just the
## hand-rolled mk_context() fixture above. This is the exact gap that let the
## "## " vs "### " header-regex bug (and the hardcoded-title-drift bug) ship
## silently: every existing test built context by hand instead of calling the
## real build_tx_context()/build_mx_context()/build_cx_context().
## ---------------------------------------------------------------------------

test_that("the guard fires against real build_tx_context() output for an unrun sub-module", {
  dataset <- list(source = "demo", expr = matrix(1, 2, 2), meta = data.frame(group = c("Case", "Control")))
  ctx <- build_tx_context(dataset, results = list(), focus_id = "dge")
  resp <- "Differential Expression found 42 significant genes, including GENE1 and GENE2."
  chk <- arthochat_detect_ungrounded_reference(resp, ctx)
  expect_true(chk$flagged)
  expect_true("Differential Expression" %in% chk$modules)
})

test_that("the guard does not fire against real build_tx_context() output once the sub-module has results", {
  dataset <- list(source = "demo", expr = matrix(1, 2, 2), meta = data.frame(group = c("Case", "Control")))
  ctx <- build_tx_context(dataset, results = list(dge = list(n_sig = 42, contrast = "Case vs Control")), focus_id = "dge")
  resp <- "Differential Expression found 42 significant genes."
  chk <- arthochat_detect_ungrounded_reference(resp, ctx)
  expect_false(chk$flagged)
  expect_match(arthochat_grounded_modules_label(ctx), "Differential Expression", fixed = TRUE)
})

## ---------------------------------------------------------------------------
## (D) case-insensitive contrast/group matching in the DGE tool
## ---------------------------------------------------------------------------

dge_ci_fixture <- function(n_per_group = 6, seed = 91) {
  set.seed(seed)
  n <- n_per_group * 2
  genes <- paste0("GENE", 1:40)
  samples <- paste0("S", 1:n)
  grp <- rep(c("HC", "RA"), each = n_per_group)
  m <- matrix(rnorm(40 * n, mean = 8, sd = 1.2), 40, n, dimnames = list(genes, samples))
  m[1:5, grp == "RA"] <- m[1:5, grp == "RA"] + 2
  meta <- data.frame(sample = samples, group = grp, stringsAsFactors = FALSE)
  shiny::reactiveValues(expr = m, meta = meta, source = "ci test cohort", source_type = "uploaded",
                          is_bundled_reference = FALSE, geo_ids = character(0))
}

test_that("a differently-cased contrast column and group values still resolve and fit correctly", {
  dataset <- dge_ci_fixture()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    ## Real metadata has column "group" with values "HC"/"RA" - request the
    ## contrast using different casing throughout.
    fit <- compute_dge_fit(contrast_col = "GROUP", ref_group = "hc", comp_group = "Ra",
                            covariate_col = NULL, covariate_mode = NULL, covariate_level = NULL,
                            method = "limma")
    expect_true(is.data.frame(fit$table) || is.data.frame(fit))
  })
})

test_that("a genuinely nonexistent group value is still rejected (case-insensitivity does not weaken real validation)", {
  dataset <- dge_ci_fixture()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    err <- tryCatch(
      compute_dge_fit(contrast_col = "group", ref_group = "hc", comp_group = "nonexistent_level",
                       covariate_col = NULL, covariate_mode = NULL, covariate_level = NULL, method = "limma"),
      error = function(e) e
    )
    expect_true(inherits(err, "error"))
  })
})

test_that("an exact-case contrast (the button-driven path's normal case) behaves exactly as before", {
  dataset <- dge_ci_fixture()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    fit <- compute_dge_fit(contrast_col = "group", ref_group = "HC", comp_group = "RA",
                            covariate_col = NULL, covariate_mode = NULL, covariate_level = NULL,
                            method = "limma")
    expect_true(is.data.frame(fit$table) || is.data.frame(fit))
  })
})

test_that("requesting the same group as both reference and comparison (via different casing) is still rejected as 'must be different'", {
  dataset <- dge_ci_fixture()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    err <- tryCatch(
      compute_dge_fit(contrast_col = "group", ref_group = "hc", comp_group = "HC",
                       covariate_col = NULL, covariate_mode = NULL, covariate_level = NULL, method = "limma"),
      error = function(e) e
    )
    expect_true(inherits(err, "shiny.silent.error") || inherits(err, "error"))
  })
})
