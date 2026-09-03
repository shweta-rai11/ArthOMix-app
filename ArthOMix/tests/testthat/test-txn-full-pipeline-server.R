## Module 1 (Transcriptomics) - FULL CROSS-MODULE PIPELINE integration test:
## Dataset -> Differential Expression -> WGCNA -> Candidate Gene
## Identification -> Feature Selection, chained through the REAL shared

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
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
  set.seed(500)
  n_genes <- 80; n_per_group <- 12; n <- n_per_group * 2
  genes <- paste0("gene", 1:n_genes)
  samples <- paste0("S", 1:n)
  grp <- rep(c("HC", "RA"), each = n_per_group)

  latent <- rnorm(n)
  expr <- matrix(rnorm(n_genes * n, mean = 8, sd = 1), n_genes, n, dimnames = list(genes, samples))
  loadings <- runif(20, 0.8, 1.2)
  for (i in 1:20) expr[i, ] <- expr[i, ] + loadings[i] * latent
  expr[1:10, grp == "RA"] <- expr[1:10, grp == "RA"] + 3

  meta <- data.frame(sample = samples, group = grp, stringsAsFactors = FALSE)
  dataset <- shiny::reactiveValues(expr = expr, meta = meta, source = "full-pipeline test cohort",
                                     source_type = "uploaded", is_bundled_reference = FALSE, geo_ids = character(0))
  results <- shiny::reactiveValues()

  shiny::testServer(mod_dge_server, args = list(id = "dge", dataset = dataset, results = results), {
    session$setInputs(data_source = "pipeline", method = "limma", contrast_col = "group",
                        ref_group = "HC", comp_group = "RA", padj_cut = 0.05, lfc_cut = 0.1)
    session$setInputs(run_btn = 0)
    session$setInputs(run_btn = 1)
    invisible(fit_result())
  })
  dge_id <- shiny::isolate({
    expect_false(is.null(results$dge_runs))
    id <- names(results$dge_runs)[1]
    dge_tab <- results$dge_runs[[id]]$table
    expect_true(sum(dge_tab$direction[dge_tab$gene %in% paste0("gene", 1:10)] != "Not significant") >= 5)
    id
  })

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
    best_overlap <- max(vapply(setdiff(names(results$wgcna$module_genes), "grey"), function(m) {
      length(intersect(results$wgcna$module_genes[[m]], real_module_genes))
    }, integer(1)))
    expect_true(best_overlap >= 10)
  })

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
    expect_true(all(r$overlap %in% real_module_genes))

    session$flushReact()
  })
  shiny::isolate({
    expect_false(is.null(results$candidates$final))
    expect_true(length(results$candidates$final$genes) > 0)
    expect_true(all(results$candidates$final$genes %in% paste0("gene", 1:20)))
  })

  shiny::testServer(mod_featureselection_server, args = list(id = "fs", dataset = dataset, results = results), {
    live <- project_candidate_genes("pooled")
    expect_true(isTRUE(live$is_live))
    expect_true(length(live$genes) > 0)
    expect_setequal(live$genes, results$candidates$final$genes)
  })
})
