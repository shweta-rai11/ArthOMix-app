## Module 1 (Transcriptomics) - Candidate Gene Identification: WGCNA module
## background intersected with DEG lists (sex-stratified when supported),
## optional MR/Colocalization causal-evidence refinement, and the
## results$candidates hand-off other tabs read from - via testServer().

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "mod_candidates.R"))

## Prime helper for the ignoreInit=TRUE bare-actionButton gotcha (see
## feedback_shiny_testserver_ignoreinit_actionbutton_priming memory) -
## sex_candidates()'s eventReactive uses this exact pattern.
click <- function(session, input_id) {
  do.call(session$setInputs, setNames(list(0), input_id))
  do.call(session$setInputs, setNames(list(1), input_id))
}

cand_wgcna_fixture <- function() {
  set.seed(100)
  turquoise <- paste0("GENE", 1:30)
  blue <- paste0("GENE", 31:50)
  list(module_genes = list(turquoise = turquoise, blue = blue, grey = paste0("GENE", 51:60)),
       significant_trait_modules = "turquoise")
}

## A DGE run's saved table shape (matches results$dge_runs[[id]]$table from mod_dge.R).
cand_dge_run <- function(contrast, sig_genes, all_genes) {
  n <- length(all_genes)
  df <- data.frame(gene = all_genes, logFC = rnorm(n), adj.P.Val = runif(n, 0.1, 0.9),
                    direction = "Not significant", stringsAsFactors = FALSE)
  df$adj.P.Val[df$gene %in% sig_genes] <- 0.01
  df$direction[df$gene %in% sig_genes] <- ifelse(df$logFC[df$gene %in% sig_genes] > 0, "Up", "Down")
  list(contrast = contrast, table = df)
}

cand_dataset <- function(with_sex = TRUE, n_samples = 20) {
  set.seed(101)
  genes <- paste0("GENE", 1:60)
  samples <- paste0("S", 1:n_samples)
  expr <- matrix(rnorm(60 * n_samples, mean = 8), 60, n_samples, dimnames = list(genes, samples))
  meta <- data.frame(sample = samples, group = rep(c("HC", "RA"), length.out = n_samples), stringsAsFactors = FALSE)
  if (with_sex) meta$sex <- rep(c("F", "M"), length.out = n_samples)
  shiny::reactiveValues(expr = expr, meta = meta, source = "candidates test",
                          source_type = "uploaded", is_bundled_reference = FALSE)
}

test_that("prereqs() correctly gates on both WGCNA and DGE having been run this session", {
  dataset <- cand_dataset()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_candidates_server, args = list(id = "cand", dataset = dataset, results = results), {
    p <- prereqs()
    expect_false(p$wgcna_ok)
    expect_false(p$dge_ok)

    results$wgcna <- cand_wgcna_fixture()
    p2 <- prereqs()
    expect_true(p2$wgcna_ok)
    expect_false(p2$dge_ok)

    results$dge_runs <- list(run1 = cand_dge_run("RA vs HC", paste0("GENE", 1:10), paste0("GENE", 1:60)))
    p3 <- prereqs()
    expect_true(p3$wgcna_ok)
    expect_true(p3$dge_ok)
  })
})

test_that("module_choices() reports the disease-associated module as pre-selected among the non-grey modules", {
  dataset <- cand_dataset()
  results <- shiny::reactiveValues(wgcna = cand_wgcna_fixture())
  shiny::testServer(mod_candidates_server, args = list(id = "cand", dataset = dataset, results = results), {
    mc <- module_choices()
    expect_setequal(mc$mods, c("turquoise", "blue"))
    expect_equal(mc$disease, "turquoise")
  })
})

test_that("Case 1 (M+F): sex_available() is TRUE with both sexes represented", {
  dataset <- cand_dataset(with_sex = TRUE)
  shiny::testServer(mod_candidates_server, args = list(id = "cand", dataset = dataset, results = shiny::reactiveValues()), {
    expect_true(sex_available())
  })
})

test_that("Case 2/3 (single-sex-only data): sex_available() is FALSE, falling back to the pooled panel", {
  dataset <- cand_dataset(with_sex = TRUE)
  shiny::isolate(dataset$meta$sex <- "F")  ## collapse to a single represented sex
  shiny::testServer(mod_candidates_server, args = list(id = "cand", dataset = dataset, results = shiny::reactiveValues()), {
    expect_false(sex_available())
  })
})

test_that("Case 4 (no sex metadata column at all): sex_available() is FALSE (the documented pooled/fallback behavior)", {
  dataset <- cand_dataset(with_sex = FALSE)
  shiny::testServer(mod_candidates_server, args = list(id = "cand", dataset = dataset, results = shiny::reactiveValues()), {
    expect_false(sex_available())
  })
})

test_that("the female panel computes the correct module-background ∩ significant-DEG overlap with a valid hypergeometric p-value", {
  dataset <- cand_dataset()
  wg <- cand_wgcna_fixture()
  ## 5 significant DEGs, all inside the turquoise module (GENE1-30) - a clean, known overlap.
  sig <- c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5")
  results <- shiny::reactiveValues(wgcna = wg, dge_runs = list(f1 = cand_dge_run("RA vs HC (female)", sig, paste0("GENE", 1:60))))
  shiny::testServer(mod_candidates_server, args = list(id = "cand", dataset = dataset, results = results), {
    session$setInputs(wgcna_module_choice = "turquoise")
    session$setInputs(gene_panel_choice = "")
    session$setInputs(female_deg_run = "f1")
    click(session, "female_run_btn")

    r <- female_result()
    expect_setequal(r$overlap, sig)
    expect_equal(r$n_bg, 30L)
    expect_equal(r$n_deg, 5L)
    expect_true(r$p_value >= 0 && r$p_value <= 1)
    expect_true(all(c("gene", "logFC", "adj.P.Val", "direction", "wgcna_module") %in% colnames(r$stats)))
    expect_true(all(r$stats$wgcna_module == "turquoise"))
  })
})

test_that("no overlap between the module background and the DEG list is rejected with a clear validation error", {
  dataset <- cand_dataset()
  wg <- cand_wgcna_fixture()
  ## Significant genes entirely outside turquoise/blue (in the "grey"/excluded pool).
  sig <- c("GENE51", "GENE52")
  results <- shiny::reactiveValues(wgcna = wg, dge_runs = list(f1 = cand_dge_run("RA vs HC (female)", sig, paste0("GENE", 1:60))))
  shiny::testServer(mod_candidates_server, args = list(id = "cand", dataset = dataset, results = results), {
    session$setInputs(wgcna_module_choice = "turquoise")
    session$setInputs(gene_panel_choice = "")
    session$setInputs(female_deg_run = "f1")
    click(session, "female_run_btn")
    err <- tryCatch(female_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("No genes are shared between the module background and this DEG list", conditionMessage(err)))
  })
})

test_that("a DEG contrast with no significant genes at all is rejected before any overlap is attempted", {
  dataset <- cand_dataset()
  wg <- cand_wgcna_fixture()
  results <- shiny::reactiveValues(wgcna = wg, dge_runs = list(f1 = cand_dge_run("RA vs HC (female)", character(0), paste0("GENE", 1:60))))
  shiny::testServer(mod_candidates_server, args = list(id = "cand", dataset = dataset, results = results), {
    session$setInputs(wgcna_module_choice = "turquoise")
    session$setInputs(gene_panel_choice = "")
    session$setInputs(female_deg_run = "f1")
    click(session, "female_run_btn")
    err <- tryCatch(female_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("no significant genes", conditionMessage(err)))
  })
})

test_that("union/intersection/female-only/male-only final-set selection all compute the correct gene sets", {
  dataset <- cand_dataset()
  wg <- cand_wgcna_fixture()
  sig_f <- c("GENE1", "GENE2", "GENE3")
  sig_m <- c("GENE2", "GENE3", "GENE4")
  results <- shiny::reactiveValues(
    wgcna = wg,
    dge_runs = list(
      f1 = cand_dge_run("RA vs HC (female)", sig_f, paste0("GENE", 1:60)),
      m1 = cand_dge_run("RA vs HC (male)", sig_m, paste0("GENE", 1:60))
    )
  )
  shiny::testServer(mod_candidates_server, args = list(id = "cand", dataset = dataset, results = results), {
    session$setInputs(wgcna_module_choice = "turquoise")
    session$setInputs(gene_panel_choice = "")
    session$setInputs(female_deg_run = "f1", male_deg_run = "m1")
    click(session, "female_run_btn")
    click(session, "male_run_btn")

    session$setInputs(final_candidate_set = "union")
    expect_setequal(final_candidates()$genes, c("GENE1", "GENE2", "GENE3", "GENE4"))

    session$setInputs(final_candidate_set = "intersection")
    expect_setequal(final_candidates()$genes, c("GENE2", "GENE3"))

    session$setInputs(final_candidate_set = "female")
    expect_setequal(final_candidates()$genes, sig_f)

    session$setInputs(final_candidate_set = "male")
    expect_setequal(final_candidates()$genes, sig_m)
  })
})

test_that("results$candidates is published with the correct female/male/final structure after both panels run", {
  dataset <- cand_dataset()
  wg <- cand_wgcna_fixture()
  sig_f <- c("GENE1", "GENE2")
  sig_m <- c("GENE2", "GENE3")
  results <- shiny::reactiveValues(
    wgcna = wg,
    dge_runs = list(
      f1 = cand_dge_run("RA vs HC (female)", sig_f, paste0("GENE", 1:60)),
      m1 = cand_dge_run("RA vs HC (male)", sig_m, paste0("GENE", 1:60))
    )
  )
  shiny::testServer(mod_candidates_server, args = list(id = "cand", dataset = dataset, results = results), {
    session$setInputs(wgcna_module_choice = "turquoise")
    session$setInputs(gene_panel_choice = "")
    session$setInputs(female_deg_run = "f1", male_deg_run = "m1")
    click(session, "female_run_btn")
    click(session, "male_run_btn")
    session$setInputs(final_candidate_set = "union")
    session$flushReact()

    expect_false(is.null(results$candidates))
    expect_setequal(results$candidates$female$genes, sig_f)
    expect_setequal(results$candidates$male$genes, sig_m)
    expect_setequal(results$candidates$final$genes, union(sig_f, sig_m))
    expect_equal(results$candidates$final$selection, "union")
  })
})

test_that("requiring MR support filters the final candidate set down to only MR-supported genes", {
  dataset <- cand_dataset()
  wg <- cand_wgcna_fixture()
  sig <- c("GENE1", "GENE2", "GENE3")
  results <- shiny::reactiveValues(
    wgcna = wg,
    dge_runs = list(f1 = cand_dge_run("RA vs HC (female)", sig, paste0("GENE", 1:60))),
    ## mr_supported_genes() keeps genes with p < 0.05 from results$mr$genes_tested.
    mr = list(genes_tested = list(GENE1 = list(p = 0.01), GENE2 = list(p = 0.5), GENE3 = list(p = 0.02)))
  )
  shiny::testServer(mod_candidates_server, args = list(id = "cand", dataset = dataset, results = results), {
    session$setInputs(wgcna_module_choice = "turquoise")
    session$setInputs(gene_panel_choice = "")
    session$setInputs(female_deg_run = "f1")
    click(session, "female_run_btn")
    session$setInputs(final_candidate_set = "female")  ## only the female panel was run

    expect_setequal(mr_supported_genes(), c("GENE1", "GENE3"))

    session$setInputs(require_mr = TRUE)
    rf <- refined_final()
    expect_setequal(rf$genes, c("GENE1", "GENE3"))
  })
})
