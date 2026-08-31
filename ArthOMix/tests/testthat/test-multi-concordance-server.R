## Module 3 (Multiomics) - Gene-CpG Concordance sub-module, via testServer():
## the "Run Gene-CpG Analysis" pipeline is SYNCHRONOUS (no ExtendedTask,
## unlike Integration/Stratification/Biomarker Discovery) - so this drives
## a REAL, full click-through run: real gene->CpG annotation mapping (real
## 450K manifest), real per-feature expression/methylation stats (real
## Welch t-tests), and real per-pair correlation, using genuine CpG probe
## IDs looked up from the real annotation for TP53 rather than invented ones.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "crossomics", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_plots.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_live_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "cohort_harmonization_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_integration_live_helpers.R"))
source_from_app_root(file.path("R", "methylomics", "qc.R"))
source_from_app_root(file.path("R", "methylomics", "celltype.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_concordance_live_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "multiomics_concordance_plots.R"))
source_from_app_root(file.path("R", "multiomics", "snf_clustering_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "mod_multi_concordance.R"))

test_that("clicking 'Run Gene-CpG Analysis' on a real live dataset computes a real end-to-end pooled result (real annotation, real stats, real correlation)", {
  skip_if_not(requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE), "450K annotation not installed")
  ar <- cx_get_region_annotation("450K")
  skip_if_not(isTRUE(ar$ok), "450K annotation failed to load")
  real_cpgs <- rownames(ar$anno)[toupper(ar$anno$gene) %in% "TP53"]
  skip_if(length(real_cpgs) < 2, "fewer than 2 real CpGs annotated to TP53 in this deployment's manifest")
  real_cpgs <- utils::head(real_cpgs, 3)

  set.seed(1400)
  n <- 20
  ids <- paste0("S", seq_len(n))
  grp <- rep(c("HC", "RA"), each = n / 2)
  expr <- matrix(rnorm(n * 3), n, 3, dimnames = list(ids, c("TP53", "BRCA1", "EGFR")))
  expr[grp == "RA", "TP53"] <- expr[grp == "RA", "TP53"] + 3  ## real, planted group difference
  meth <- matrix(runif(n * length(real_cpgs), 0.3, 0.5), n, length(real_cpgs), dimnames = list(ids, real_cpgs))
  meth[grp == "RA", ] <- meth[grp == "RA", ] + 0.15
  meta <- data.frame(design = grp, row.names = ids, stringsAsFactors = FALSE)

  multi_dataset <- shiny::reactiveValues(
    active = TRUE, source = "upload", layers = list(Transcriptomics = expr, Methylomics = meth),
    layer_meta = list(Transcriptomics = list(omics_type = "rnaseq"), Methylomics = list(omics_type = "methylation")),
    sample_meta = meta, table_label = NULL
  )
  multi_results <- shiny::reactiveValues()

  shiny::testServer(mod_multi_concordance_server, args = list(id = "cc", multi_dataset = multi_dataset, multi_results = multi_results), {
    session$setInputs(
      data_source = "active", expr_layer = "Transcriptomics", meth_layer = "Methylomics",
      source = "All candidates", sex = "all", array_type = "450K",
      custom_genes = "TP53", custom_cpgs = character(0), design_col = "design",
      expr_thresh = 1, expr_fdr_thresh = 0.05, meth_thresh = 0.1, meth_fdr_thresh = 0.05,
      cor_method = "pearson"
    )
    session$setInputs(run_btn = 1)

    r <- state$result
    expect_true(isTRUE(r$ok))
    expect_true(nrow(r$pairs_df) > 0)
    expect_true(all(r$pairs_df$gene_symbol == "TP53"))
    expect_setequal(r$pairs_df$sex, "pooled")
    ## The planted TP53 group difference should show up as a real,
    ## significant expression effect - not a fabricated statistic.
    expect_true(all(r$pairs_df$log2fc > 1))
    expect_true(all(r$pairs_df$expr_p < 0.05))
  })
})

test_that("the Run handler reports a clear, non-crashing error when no 2-class design column is selected", {
  set.seed(1410)
  n <- 10
  ids <- paste0("S", seq_len(n))
  expr <- matrix(rnorm(n * 3), n, 3, dimnames = list(ids, c("TP53", "BRCA1", "EGFR")))
  meth <- matrix(runif(n * 3, 0.3, 0.5), n, 3, dimnames = list(ids, c("cg1", "cg2", "cg3")))
  multi_dataset <- shiny::reactiveValues(
    active = TRUE, source = "upload", layers = list(Transcriptomics = expr, Methylomics = meth),
    layer_meta = list(Transcriptomics = list(omics_type = "rnaseq"), Methylomics = list(omics_type = "methylation")),
    sample_meta = NULL
  )
  multi_results <- shiny::reactiveValues()

  shiny::testServer(mod_multi_concordance_server, args = list(id = "cc", multi_dataset = multi_dataset, multi_results = multi_results), {
    session$setInputs(
      data_source = "active", expr_layer = "Transcriptomics", meth_layer = "Methylomics",
      source = "All candidates", sex = "all", array_type = "450K",
      custom_genes = "TP53", custom_cpgs = character(0)
    )
    session$setInputs(run_btn = 1)
    r <- state$result
    expect_false(isTRUE(r$ok))
    expect_true(nzchar(r$error))
  })
})
