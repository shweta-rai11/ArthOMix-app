## Module 2 (Methylomics) - DMR tab's "SVA" tab, live SVA-adjusted engine
## (mod_methyl_dmr_svalive_panel_ui() / the "1b." server block in
## mod_methyl_dmr.R), via testServer(). Verifies this mirrors
## mod_methyl_dmp.R's own SVA tab exactly: when the active dataset isn't
## the preloaded GSE42861 cohort, the "SVA" tab itself (not the separate
## "DMR" tab) swaps from the precomputed display to a live sva::sva() +
## bacon::bacon() engine feeding DMRcate region calling, run against
## Upload/GEO-fetched data - and that the precomputed preloaded path is
## unaffected.
##
## Unlike test-methyl-dmr-functions.R's documented coverage gap (a small
## synthetic fixture can't give DMRcate anything to find), this fixture uses
## real EPIC manifest CpG IDs/positions so the kernel-smoothed region caller
## has genuine genomic adjacency to work with: one real, dense ~3kb cluster
## of CpGs on chr1 carries the injected group effect, and a scattered set of
## real background CpGs elsewhere carries a plate/batch confound for SVA to
## soak up (same fixture design as test-methyl-dmp-sva-live-server.R).

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "qc.R"))
source_from_app_root(file.path("R", "methylomics", "annotation.R"))
source_from_app_root(file.path("R", "methylomics", "normalization.R"))
source_from_app_root(file.path("R", "provenance.R"))
source_from_app_root(file.path("R", "methylomics", "mod_methyl_dmp.R"))
source_from_app_root(file.path("R", "methylomics", "mod_methyl_dmr.R"))

## Locates the densest real EPIC CpG cluster on chr1 within a 3kb window
## (skipped if the EPIC annotation package isn't installed in this
## deployment - same guard methyl_dmr_engine_pkgs_ok() exists for).
dmr_sva_find_cluster <- function(anno, window_bp = 3000) {
  a1 <- anno[!is.na(anno$chr) & !is.na(anno$pos) & anno$chr == "chr1", ]
  a1 <- a1[order(a1$pos), ]
  pos <- a1$pos
  n <- length(pos)
  best_count <- 0; best_i <- 1; best_j <- 1
  for (i in seq_len(n)) {
    j <- i
    while (j <= n && pos[j] - pos[i] <= window_bp) j <- j + 1
    if (j - i > best_count) { best_count <- j - i; best_i <- i; best_j <- j - 1 }
  }
  list(ids = rownames(a1)[best_i:best_j], range = c(pos[best_i], pos[best_j]))
}

## Real-CpG-ID fixture: `cluster$ids` (real, genomically adjacent) get a
## true methylation shift in the RA group (the signal DMRcate should
## recover as one region); a block of real, scattered background CpGs
## carries a plate/batch confound uncorrelated with group or sex, for SVA
## to estimate away (mirrors sva_fixture_dataset() in
## test-methyl-dmp-sva-live-server.R).
dmr_sva_fixture_dataset <- function(n_per_group = 20, seed = 411, source_type = "uploaded") {
  ar <- methyl_get_annotation("EPIC")
  testthat::skip_if_not(isTRUE(ar$ok), "EPIC annotation package not available in this deployment")
  cluster <- dmr_sva_find_cluster(ar$anno)
  testthat::skip_if(length(cluster$ids) < 20, "Could not locate a large enough real CpG cluster to seed a DMR")

  bg <- ar$anno[!is.na(ar$anno$chr) & !is.na(ar$anno$pos) & ar$anno$chr %in% paste0("chr", 2:5), ]
  set.seed(seed)
  bg_ids <- sample(rownames(bg), 300)
  all_ids <- c(cluster$ids, bg_ids)

  n <- n_per_group * 2
  m <- matrix(stats::runif(length(all_ids) * n, 0.2, 0.8), length(all_ids), n,
              dimnames = list(all_ids, paste0("S", seq_len(n))))
  ## Real signal on the clustered CpGs for the RA group.
  m[cluster$ids, (n_per_group + 1):n] <- pmin(m[cluster$ids, (n_per_group + 1):n] + 0.35, 0.99)
  ## Batch confound spread across a block of background probes.
  batch <- rep(c("plate1", "plate2"), length.out = n)[sample.int(n)]
  batch_shift <- ifelse(batch == "plate2", 0.15, 0)
  bg_block <- bg_ids[1:150]
  m[bg_block, ] <- pmin(pmax(sweep(m[bg_block, ], 2, batch_shift, `+`), 0.01), 0.99)

  sheet <- data.frame(sample = colnames(m), group = rep(c("HC", "RA"), each = n_per_group),
                        sex = rep(c("F", "M"), length.out = n), batch = batch, stringsAsFactors = FALSE)
  list(
    cluster_range = cluster$range,
    methyl_dataset = shiny::reactiveValues(
      beta = m, sample_sheet = sheet, input_scale = "beta", array_type = "EPIC",
      rg_set = NULL, mset = NULL, detp = NULL, beadcount = NULL,
      preloaded = FALSE, source_type = source_type, source = "dmr-sva-live test"
    )
  )
}

test_that("the SVA tab's default_ui swaps to the live SVA-adjusted DMR panel for non-preloaded data (mirrors mod_methyl_dmp.R's own SVA tab)", {
  fx <- dmr_sva_fixture_dataset()
  methyl_results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_dmr_server, args = list(id = "dmr", methyl_dataset = fx$methyl_dataset, methyl_results = methyl_results), {
    html <- output$default_ui$html
    expect_true(grepl("SVA-adjusted DMR Analysis", html))
    expect_true(grepl("svalive_run_btn", html))
    expect_false(grepl("DMR configuration", html))  ## the precomputed panel's own header must NOT appear
  })
})

test_that("the SVA tab's default_ui still shows the precomputed panel when the dataset IS preloaded (regression check)", {
  methyl_dataset <- shiny::reactiveValues(beta = NULL, sample_sheet = NULL, input_scale = "beta", array_type = "450K",
                                            rg_set = NULL, mset = NULL, detp = NULL, beadcount = NULL,
                                            preloaded = TRUE, source_type = "preloaded", source = "GSE42861")
  methyl_results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_dmr_server, args = list(id = "dmr", methyl_dataset = methyl_dataset, methyl_results = methyl_results), {
    html <- output$default_ui$html
    expect_true(grepl("DMR configuration", html))
    expect_false(grepl("SVA-adjusted DMR Analysis", html))
  })
})

test_that("the live SVA-adjusted DMR engine estimates SVs, corrects inflation, and recovers the injected region", {
  fx <- dmr_sva_fixture_dataset()
  methyl_results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_dmr_server, args = list(id = "dmr", methyl_dataset = fx$methyl_dataset, methyl_results = methyl_results), {
    session$setInputs(svalive_sex = "__all__", svalive_group_col = "group", svalive_ref = "HC", svalive_comp = "RA",
                        svalive_min_valid_pct = 80, svalive_min_variance = 0, svalive_snp_filter = FALSE,
                        svalive_covariates = character(0), svalive_seed_p = 0.05, svalive_dbeta = 0.05, svalive_mincpgs = 3)
    session$setInputs(svalive_run_btn = 1)

    r <- svalive_result()
    expect_true(is.numeric(r$n_sv) && r$n_sv >= 0)
    expect_true(grepl("surrogate variable", r$design_formula))
    expect_true(is.finite(r$lambda_bacon) && r$lambda_bacon > 0)
    expect_true(is.finite(r$bias_bacon))
    expect_true(nrow(r$dt) >= 1)
    expect_true(all(c("dmr_fdr", "meandiff", "no.cpgs") %in% colnames(r$dt)))

    ## The region DMRcate calls should overlap the real cluster the group
    ## effect was actually injected into.
    hit <- r$dt$seqnames == "chr1" & r$dt$start <= fx$cluster_range[2] & r$dt$end >= fx$cluster_range[1]
    expect_true(any(hit))
    expect_true(r$dt$dmr_fdr[which(hit)[1]] < 0.05)

    ## Publishes to methyl_results for downstream Candidate CpGs pickup,
    ## same convention as the plain "DMR" tab's own live engine.
    expect_equal(methyl_results$dmr_table, r$dt)
  })
})

test_that("the live SVA-adjusted DMR engine also runs for a live GEO-fetched dataset (not just Upload)", {
  fx <- dmr_sva_fixture_dataset(seed = 99, source_type = "geo")
  methyl_results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_dmr_server, args = list(id = "dmr", methyl_dataset = fx$methyl_dataset, methyl_results = methyl_results), {
    session$setInputs(svalive_sex = "__all__", svalive_group_col = "group", svalive_ref = "HC", svalive_comp = "RA",
                        svalive_min_valid_pct = 80, svalive_min_variance = 0, svalive_snp_filter = FALSE,
                        svalive_covariates = character(0), svalive_seed_p = 0.05, svalive_dbeta = 0.05, svalive_mincpgs = 3)
    session$setInputs(svalive_run_btn = 1)

    r <- svalive_result()
    expect_true(nrow(r$dt) >= 1)
    hit <- r$dt$seqnames == "chr1" & r$dt$start <= fx$cluster_range[2] & r$dt$end >= fx$cluster_range[1]
    expect_true(any(hit))
  })
})

test_that("identical reference and comparison groups are rejected by the live SVA-adjusted DMR engine", {
  fx <- dmr_sva_fixture_dataset()
  methyl_results <- shiny::reactiveValues()
  shiny::testServer(mod_methyl_dmr_server, args = list(id = "dmr", methyl_dataset = fx$methyl_dataset, methyl_results = methyl_results), {
    session$setInputs(svalive_sex = "__all__", svalive_group_col = "group", svalive_ref = "HC", svalive_comp = "HC")
    session$setInputs(svalive_run_btn = 1)
    err <- tryCatch(svalive_result(), error = function(e) e)
    expect_s3_class(err, "validation")
    expect_true(grepl("must be different", conditionMessage(err)))
  })
})
