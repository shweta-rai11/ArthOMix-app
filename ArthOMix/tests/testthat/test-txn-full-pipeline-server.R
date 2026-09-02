## Module 1 (Transcriptomics) - FULL CROSS-MODULE PIPELINE integration test:
## Dataset -> Differential Expression -> WGCNA -> Candidate Gene
## Identification -> Feature Selection, chained through the REAL shared
## `dataset`/`results` reactiveValues exactly as the running app does (each
## stage's testServer() call receives the SAME reactiveValues objects the
## previous stage populated - nothing is hand-copied between stages).
##
## This is the thing per-module tests (which each construct their own
## plausible-looking `results$dge_runs`/`results$wgcna` fixture by hand)
## cannot prove on their own: that DGE's REAL output shape is something
## WGCNA-independent Candidates can REALLY consume, and that Candidates'
## REAL output shape is something Feature Selection can REALLY consume -
## i.e. genuine cross-module communication, not just each module's own
## internal correctness in isolation.
##
## Confirmed by reading each module's source first (not assumed): DGE
## publishes to results$dge_runs; WGCNA publishes to results$wgcna
## (module_genes); Candidates reads both of those and publishes
## results$candidates$final; Feature Selection reads results$candidates$final/
## female/male. mod_preprocessing.R was checked too and deliberately does NOT
## write back into `dataset`/`results` at all (its own "Download merged
## data" buttons are an explicit export-and-manually-reupload checkpoint,
## not an automatic pipeline stage) - so it is intentionally not part of
## this chain; that is the module's real, by-design behavior, not a gap.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "01_Data", "crossomics_integration_upload.R"))
source_from_app_root(file.path("R", "provenance.R"))
source_from_app_root(file.path("R", "transcriptomics", "functions", "expression_type.R"))
source_from_app_root(file.path("R", "transcriptomics", "04_Differential_Expression", "mod_dge.R"))
source_from_app_root(file.path("R", "transcriptomics", "05_WGCNA", "mod_wgcna.R"))
source_from_app_root(file.path("R", "transcriptomics", "06_Candidate_Gene_Identification", "mod_candidates.R"))
source_from_app_root(file.path("R", "transcriptomics", "09_Feature_Selection", "mod_featureselection.R"))

test_that("a real dataset flows through DGE -> WGCNA -> Candidates -> Feature Selection via the shared dataset/results reactiveValues, with real computation at every stage", {
  skip_if_not_installed("WGCNA")
  ## Restore R's default RNG kind before seeding - a package loaded earlier
  ## in the same suite run (parallel-backend registration prints "Allowing
  ## parallel execution..." repeatedly throughout this suite) can leave
  ## RNGkind() switched to "L'Ecuyer-CMRG" globally, which makes set.seed()
  ## produce a DIFFERENT random stream than in an isolated run even with an
  ## identical seed value - silently changing this fixture's planted
  ## signal and making the WGCNA/DGE assertions below order-dependent
  ## rather than genuinely deterministic. Pin it explicitly so this test
  ## behaves identically alone or anywhere in the full suite.
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
  set.seed(500)
  n_genes <- 80; n_per_group <- 12; n <- n_per_group * 2
  genes <- paste0("gene", 1:n_genes)
  samples <- paste0("S", 1:n)
  grp <- rep(c("HC", "RA"), each = n_per_group)

  ## Genes 1-20 form one real co-expressed module (a shared per-sample latent
  ## factor drives them together, giving WGCNA real correlation structure to
  ## detect) - genes 1-10 of that module ALSO carry a real group mean-shift,
  ## so they are simultaneously real WGCNA-module members AND real DGE hits,
  ## which is exactly the overlap Candidate Gene Identification is built to
  ## find. Genes 21-80 are independent background noise.
  latent <- rnorm(n)
  expr <- matrix(rnorm(n_genes * n, mean = 8, sd = 1), n_genes, n, dimnames = list(genes, samples))
  loadings <- runif(20, 0.8, 1.2)
  for (i in 1:20) expr[i, ] <- expr[i, ] + loadings[i] * latent
  expr[1:10, grp == "RA"] <- expr[1:10, grp == "RA"] + 3  ## real, planted DEG signal

  meta <- data.frame(sample = samples, group = grp, stringsAsFactors = FALSE)
  dataset <- shiny::reactiveValues(expr = expr, meta = meta, source = "full-pipeline test cohort",
                                     source_type = "uploaded", is_bundled_reference = FALSE, geo_ids = character(0))
  results <- shiny::reactiveValues()

  ## ---- Stage 1: Differential Expression (real limma fit) --------------------
  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    session$setInputs(data_source = "pipeline", method = "limma", contrast_col = "group",
                        ref_group = "HC", comp_group = "RA", padj_cut = 0.05, lfc_cut = 0.1)
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)
    invisible(fit_result())
  })
  ## Reads of a shared reactiveValues object OUTSIDE any reactive
  ## consumer (i.e. between testServer() calls, at this test's own top
  ## level) must be wrapped in shiny::isolate() - the same rule that
  ## applies inside a real moduleServer body (see
  ## feedback_shiny_reactivevalues_setup_time_read).
  dge_id <- shiny::isolate({
    expect_false(is.null(results$dge_runs))
    id <- names(results$dge_runs)[1]
    dge_tab <- results$dge_runs[[id]]$table
    expect_true(sum(dge_tab$direction[dge_tab$gene %in% paste0("gene", 1:10)] != "Not significant") >= 5)
    id
  })

  ## ---- Stage 2: WGCNA (real blockwiseModules(), manual power to skip the ----
  ## separate soft-threshold-search step - a real, supported configuration)
  shiny::testServer(mod_wgcna_server, args = list(id = "wg", dataset = dataset, results = results), {
    session$setInputs(gene_filter_method = "all", exclude_pattern = "", remove_outliers = FALSE,
                        power_mode = "manual", manual_power = 6,
                        network_type = "signed", tom_type = "signed", cor_method = "pearson",
                        deep_split = 2, min_module_size = 5, merge_cut_height = 0.25, pam_respects_dendro = TRUE)
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)
    invisible(net_result())
  })
  real_module_genes <- paste0("gene", 1:20)
  shiny::isolate({
    expect_false(is.null(results$wgcna))
    expect_true(length(results$wgcna$module_genes) > 0)
    ## Real network structure recovered: at least one detected (non-grey)
    ## module should substantially overlap the real 20-gene co-expressed set.
    best_overlap <- max(vapply(setdiff(names(results$wgcna$module_genes), "grey"), function(m) {
      length(intersect(results$wgcna$module_genes[[m]], real_module_genes))
    }, integer(1)))
    expect_true(best_overlap >= 10)
  })

  ## ---- Stage 3: Candidate Gene Identification (real hypergeometric overlap, ----
  ## reading BOTH results$dge_runs and results$wgcna published just above)
  shiny::testServer(mod_candidates_server, args = list(id = "cand", dataset = dataset, results = results), {
    best_module <- names(which.max(vapply(setdiff(names(results$wgcna$module_genes), "grey"), function(m) {
      length(intersect(results$wgcna$module_genes[[m]], paste0("gene", 1:20)))
    }, integer(1))))
    session$setInputs(wgcna_module_choice = best_module, gene_panel_choice = "")
    session$setInputs(pooled_deg_run = dge_id)
    session$setInputs(pooled_run_btn = 0)
    session$setInputs(pooled_run_btn = 1)
    r <- pooled_result()
    expect_true(length(r$overlap) > 0)
    ## The real overlap should be drawn from the real planted intersection
    ## (module genes 1-20) AND (DEG-significant genes, concentrated in 1-10).
    expect_true(all(r$overlap %in% real_module_genes))

    ## results$candidates$final publishes automatically via an observe() on
    ## pooled_result()/refined_final() (no separate "publish" button exists,
    ## and no sex column in this fixture's meta means the pooled path is the
    ## active one) - force a reactive flush so that observer has run.
    session$flushReact()
  })
  shiny::isolate({
    expect_false(is.null(results$candidates$final))
    expect_true(length(results$candidates$final$genes) > 0)
    expect_true(all(results$candidates$final$genes %in% paste0("gene", 1:20)))
  })

  ## ---- Stage 4: Feature Selection (reads results$candidates$final, ----------
  ## published two stages upstream, through the SAME shared `results`)
  shiny::testServer(mod_featureselection_server, args = list(id = "fs", dataset = dataset, results = results), {
    ## Calls Feature Selection's own project_candidate_genes("pooled") - the
    ## exact function its "project" data-source branch uses - proving ITS OWN
    ## reactive wiring genuinely reads the real upstream candidate set
    ## published by Candidate Gene Identification two stages earlier, not a
    ## hand-fabricated stand-in.
    live <- project_candidate_genes("pooled")
    expect_true(isTRUE(live$is_live))
    expect_true(length(live$genes) > 0)
    expect_setequal(live$genes, results$candidates$final$genes)
  })
})
