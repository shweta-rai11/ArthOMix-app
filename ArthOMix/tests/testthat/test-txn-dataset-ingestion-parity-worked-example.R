## Worked example (2026-09-07 defense audit, item 5): a single real, public GEO
## accession (GSE110169, whole-blood RA/SLE/HC microarray cohort) loaded through
## two of the app's three ingestion paths and cross-checked for logical
## consistency - not byte-identical matrices (probe-level vs. gene-symbol-
## collapsed representations differ legitimately), but the same sample set with
## the same group/sex labels reaching the same shared dataset structure.
##
## Path A (preloaded): mod_dataset.R's individual_dataset_entry("GSE110169"),
## the app's own bundled raw probe-level Affymetrix data for this accession.
## Path B (upload): data/examples/transcriptomics_upload/geo_transcriptomics/
## GSE110169_{exp,sample}.csv, a gene-symbol-collapsed export of the same GEO
## series, uploaded exactly the way a reviewer would upload their own data.
## Path C (GEO fetch, live/optional): fetching GSE110169 directly from NCBI -
## network-dependent, gated behind ARTHOMIX_RUN_LIVE_GEO_TESTS like the other
## live GEO test in test-txn-dataset-geo-server.R.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "functions", "expression_type.R"))
source_from_app_root(file.path("R", "transcriptomics", "01_Data", "mod_dataset.R"))

GSE110169_EXP_CSV  <- file.path(app_dir, "data", "examples", "transcriptomics_upload", "geo_transcriptomics", "GSE110169_exp.csv")
GSE110169_META_CSV <- file.path(app_dir, "data", "examples", "transcriptomics_upload", "geo_transcriptomics", "GSE110169_sample.csv")

test_that("Path A (preloaded) and Path B (upload) of the SAME GEO accession converge to the same sample set with consistent group/sex labels", {
  skip_if_not(file.exists(GSE110169_EXP_CSV) && file.exists(GSE110169_META_CSV),
              "GSE110169 example-upload fixture not present in this checkout.")

  ## --- Path A: preloaded, via the real Dataset-tab server code ---
  dataset_a <- shiny::reactiveValues()
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset_a), {
    session$setInputs(preloaded_choice = "GSE110169")
    session$setInputs(load_preloaded_btn = 1)
  })
  expect_equal(shiny::isolate(dataset_a$source_type), "preloaded")
  expect_true(shiny::isolate(nrow(dataset_a$meta)) > 0)

  ## --- Path B: upload, via the real Dataset-tab server code, using a public,
  ## versioned export of the same accession ---
  dataset_b <- shiny::reactiveValues()
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset_b), {
    session$setInputs(expr_file = fx_mkfile(GSE110169_EXP_CSV))
    session$setInputs(meta_file = fx_mkfile(GSE110169_META_CSV))
    session$setInputs(map_id = "sample", map_group = "group", map_sex = "sex", map_batch = "batch")
    session$setInputs(load_btn = 1)
  })
  expect_equal(shiny::isolate(dataset_b$source_type), "uploaded")

  meta_a <- shiny::isolate(dataset_a$meta)
  meta_b <- shiny::isolate(dataset_b$meta)

  ## Same conceptual data type and orientation: both are sample-level matrices
  ## keyed by the same real GEO GSM accessions (not probe/row identity, which
  ## legitimately differs - Path A is raw probe-level, Path B is gene-symbol-
  ## collapsed, exactly the two representations the app itself supports).
  common <- intersect(meta_a$sample, meta_b$sample)
  expect_gt(length(common), 200)                 # both cover ~all 243 GSM samples
  expect_equal(length(common), nrow(meta_a))     # every preloaded sample is present in the upload
  expect_equal(length(common), nrow(meta_b))     # ... and vice versa (exact same sample set)

  pm <- meta_a[match(common, meta_a$sample), ]
  um <- meta_b[match(common, meta_b$sample), ]

  ## Group labels differ in spelling (preloaded: HC/RA/SLE; upload's raw GEO
  ## export: Normal/RA/SLE) but must be a CONSISTENT relabeling - every "Normal"
  ## upload sample must be "HC" in the preloaded path and nowhere else, and
  ## RA/SLE must line up exactly with no cross-contamination between disease
  ## groups. This is the actual "logically consistent results" check the
  ## worked example is for - a real, disease-label ingestion-parity bug would
  ## show up as off-diagonal entries here.
  group_map <- c(Normal = "HC", RA = "RA", SLE = "SLE")
  expect_equal(unname(group_map[um$group]), pm$group)

  ## Sex must agree exactly between the two ingestion paths for every shared sample.
  sex_common <- !is.na(pm$sex) & !is.na(um$sex) & nzchar(um$sex)
  expect_gt(sum(sex_common), 200)
  expect_equal(pm$sex[sex_common], toupper(substr(um$sex[sex_common], 1, 1)))

  ## Same conceptual scale: both are log2 microarray intensities in a plausible range.
  expr_a <- shiny::isolate(dataset_a$expr)
  expr_b <- shiny::isolate(dataset_b$expr)
  expect_true(all(is.finite(range(expr_a))) && max(expr_a) < 20 && min(expr_a) > -5)
  expect_true(all(is.finite(range(expr_b))) && max(expr_b) < 20 && min(expr_b) > -5)
})

test_that("live: Path C (GEO fetch) of the same accession also reaches the same sample set", {
  skip_if_not(
    identical(Sys.getenv("ARTHOMIX_RUN_LIVE_GEO_TESTS"), "1"),
    "Set ARTHOMIX_RUN_LIVE_GEO_TESTS=1 to run this network-dependent test."
  )
  skip_if_not(file.exists(GSE110169_META_CSV), "GSE110169 example-upload fixture not present in this checkout.")

  dataset_c <- shiny::reactiveValues()
  shiny::testServer(mod_dataset_server, args = list(id = "ds", dataset = dataset_c), {
    session$setInputs(geo_accession = "GSE110169")
    session$setInputs(geo_fetch_btn = 1)
    em <- geo_expr_meta()
    req(em, !inherits(em, "error"))
    session$setInputs(geo_map_group = "disease state:ch1", geo_map_sex = "Sex:ch1",
                       geo_map_batch = "(none)", geo_declared_data_type = "normalized")
    session$setInputs(geo_load_btn = 1)
  })
  expect_equal(shiny::isolate(dataset_c$source_type), "geo")
  up_meta <- data.table::fread(GSE110169_META_CSV, data.table = FALSE)
  common <- intersect(shiny::isolate(dataset_c$meta)$sample, up_meta$sample)
  expect_gt(length(common), 200)
})
