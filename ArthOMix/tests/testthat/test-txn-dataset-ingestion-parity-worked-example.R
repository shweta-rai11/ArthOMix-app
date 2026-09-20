## Worked example: GSE110169 via the preloaded and upload paths must give the same samples and consistent group/sex labels.
## Path C (live GEO fetch) is optional and gated behind ARTHOMIX_RUN_LIVE_GEO_TESTS.

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

  ## Both are sample-level matrices keyed by the same GSM accessions; rows legitimately differ (probe vs gene symbol).
  common <- intersect(meta_a$sample, meta_b$sample)
  expect_gt(length(common), 200)                 # both cover ~all 243 GSM samples
  expect_equal(length(common), nrow(meta_a))     # every preloaded sample is present in the upload
  expect_equal(length(common), nrow(meta_b))     # ... and vice versa (exact same sample set)

  pm <- meta_a[match(common, meta_a$sample), ]
  um <- meta_b[match(common, meta_b$sample), ]

  ## Group labels differ in spelling (HC vs Normal) but must relabel consistently, with no off-diagonal entries.
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
