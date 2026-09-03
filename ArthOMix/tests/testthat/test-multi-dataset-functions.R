## Module 3 (Multiomics) - Dataset Workspace's pure functions
## (mod_multi_dataset.R): status-badge/summary-table/provenance rendering,
## the preloaded-branch dataset blocks (real RNA-seq/methylation QC

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "multiomics_dataset_helpers.R"))
source_from_app_root(file.path("R", "multiomics", "functions", "multiomics_plots.R"))
source_from_app_root(file.path("R", "multiomics", "01_Data_Workspace", "mod_multi_dataset.R"))

test_that("mo_status_badge() renders the label and lists every reason for a non-ready status", {
  html <- htmltools::doRenderTags(mo_status_badge(list(level = "review", label = "Review Required", reasons = c("Low sample overlap", "Unrecognized omics type"))))
  expect_true(grepl("Review Required", html))
  expect_true(grepl("Low sample overlap", html))
  expect_true(grepl("Unrecognized omics type", html))
})

test_that("mo_status_badge() renders no reason list when reasons is empty (a ready status)", {
  html <- htmltools::doRenderTags(mo_status_badge(list(level = "ready", label = "Ready", reasons = character(0))))
  expect_true(grepl("Ready", html))
  expect_false(grepl("<ul", html, fixed = TRUE))
})

test_that("mo_load_first_msg() phrases the message for preloaded vs. upload/GEO sources", {
  expect_true(grepl("Load Reference Dataset", mo_load_first_msg("preloaded")))
  expect_true(grepl("Validate Datasets", mo_load_first_msg("upload")))
  expect_true(grepl("Validate Datasets", mo_load_first_msg("geo")))
})

test_that("mo_dataset_block_card() shows the label and formatted sample/feature counts", {
  html <- htmltools::doRenderTags(mo_dataset_block_card("Transcriptomics", 1234, 56789, list(level = "ready", label = "Ready", reasons = character(0))))
  expect_true(grepl("Transcriptomics", html))
  expect_true(grepl("1,234", html))
  expect_true(grepl("56,789", html))
})

test_that("mo_summary_table() builds one row per layer with real sample/feature counts from validation, 'Not processed'/'Unknown' when absent", {
  layer_meta <- list(
    expression = list(omics_type = "rnaseq", validation = list(n_samples = 40, n_features = 500), processing = "Normalized", status = list(label = "Ready")),
    methylation = list(omics_type = "methylation", validation = NULL, processing = NULL, status = NULL)
  )
  tbl <- mo_summary_table(layer_meta)
  expect_equal(nrow(tbl), 2L)
  expr_row <- tbl[tbl$Dataset == "expression", ]
  meth_row <- tbl[tbl$Dataset == "methylation", ]
  expect_equal(expr_row$Samples, 40)
  expect_equal(expr_row$Features, 500)
  expect_equal(expr_row$Processing, "Normalized")
  expect_true(is.na(meth_row$Samples))
  expect_equal(meth_row$Processing, "Not processed")
  expect_equal(meth_row$Status, "Unknown")
})

test_that("mo_summary_table() returns NULL for an empty layer_meta", {
  expect_null(mo_summary_table(list()))
})

test_that("mo_provenance_ui() shows an empty-state note with no layers, and each layer's provenance fields otherwise", {
  html_empty <- htmltools::doRenderTags(mo_provenance_ui(list()))
  expect_true(grepl("No datasets selected", html_empty))

  layer_meta <- list(expression = list(omics_type = "rnaseq", provenance = list(source = "Upload", detail = "my_expr.csv", imported_at = "2026-08-30 10:00:00")))
  html <- htmltools::doRenderTags(mo_provenance_ui(layer_meta))
  expect_true(grepl("Upload", html))
  expect_true(grepl("my_expr.csv", html))
  expect_true(grepl("2026-08-30 10:00:00", html))
})

test_that("mo_preloaded_blocks_ui() shows real RNA-seq/methylation QC-derived sample and feature counts, never hardcoded", {
  skip_if_not(MULTI_DATA_AVAILABLE, "multi-omics preloaded data not available in this deployment")
  rna <- multi_read_registry_table("RNA-seq QC summary")
  meth <- multi_read_registry_table("Methylation QC summary")
  skip_if_not(rna$ok && meth$ok, "RNA-seq/Methylation QC summary tables not available")
  rna_pbmc <- rna$df[rna$df$cell_type == "PBMC", , drop = FALSE]

  html <- htmltools::doRenderTags(mo_preloaded_blocks_ui())
  expect_true(grepl(format(rna_pbmc$n_samples[1], big.mark = ","), html, fixed = TRUE))
  expect_true(grepl(format(rna_pbmc$n_genes_retained[1], big.mark = ","), html, fixed = TRUE))
  expect_true(grepl(format(meth$df$n_samples_retained[1], big.mark = ","), html, fixed = TRUE))
  expect_true(grepl("Tao et al. 2021", html))
})

test_that("mo_block_id() prefixes by mode so upload and GEO block ids never collide", {
  expect_equal(mo_block_id(1, "upload"), "ublock1")
  expect_equal(mo_block_id(1, "geo"), "gblock1")
  expect_false(mo_block_id(3, "upload") == mo_block_id(3, "geo"))
})

test_that("mo_file_input_id()/mo_meta_file_input_id() suffix by generation so a pipeline switch always yields a fresh, never-before-used input id", {
  id_gen0 <- mo_file_input_id("ublock1", 0)
  id_gen1 <- mo_file_input_id("ublock1", 1)
  expect_false(id_gen0 == id_gen1)
  expect_true(grepl("_file_g0$", id_gen0))
  expect_true(grepl("_meta_file_g2$", mo_meta_file_input_id("ublock1", 2)))
})

test_that("mo_label_omics_type() maps the two fixed preloaded labels, 'other' for anything else", {
  expect_equal(mo_label_omics_type("Transcriptomics", NULL, 0, 0, "preloaded"), "rnaseq")
  expect_equal(mo_label_omics_type("Methylomics", NULL, 0, 0, "preloaded"), "methylation")
  expect_equal(mo_label_omics_type("Something else", NULL, 0, 0, "preloaded"), "other")
})

test_that("mo_label_omics_type() looks up an upload/GEO block's own omics-type input by matching its display label", {
  input <- list(ublock1_label = "My Expression Data", ublock1_type = "rnaseq",
                 ublock2_label = "My Methylation Data", ublock2_type = "methylation")
  expect_equal(mo_label_omics_type("My Expression Data", input, n_upload = 2, n_geo = 0, mode = "upload"), "rnaseq")
  expect_equal(mo_label_omics_type("My Methylation Data", input, n_upload = 2, n_geo = 0, mode = "upload"), "methylation")
  expect_equal(mo_label_omics_type("No such label", input, n_upload = 2, n_geo = 0, mode = "upload"), "other")
})

test_that("mo_apply_matching() with method='patient_id' remaps rownames via a metadata Patient ID column and drops unmapped samples", {
  m1 <- matrix(1:6, 3, 2, dimnames = list(c("S1", "S2", "S3"), c("f1", "f2")))
  meta <- data.frame(patient_id = c("P1", "P2"), row.names = c("S1", "S2"))
  out <- mo_apply_matching(list(rna = m1), method = "patient_id", meta = meta, patient_col = "patient_id")
  expect_setequal(rownames(out$mats$rna), c("P1", "P2"))
  expect_equal(nrow(out$mats$rna), 2L)
})

test_that("mo_apply_matching() with method='mapping' remaps rownames via a per-label mapping file column and counts dropped samples", {
  m1 <- matrix(1:6, 3, 2, dimnames = list(c("a1", "a2", "a3"), c("f1", "f2")))
  mapping_df <- data.frame(canonical = c("P1", "P2"), RNA = c("a1", "a2"), stringsAsFactors = FALSE)
  out <- mo_apply_matching(list(RNA = m1), method = "mapping", mapping_df = mapping_df)
  expect_setequal(rownames(out$mats$RNA), c("P1", "P2"))
  expect_equal(out$dropped$RNA, 1L)
})

test_that("mo_apply_matching() with method='exact' (or missing prerequisites) leaves matrices untouched", {
  m1 <- matrix(1:4, 2, 2, dimnames = list(c("S1", "S2"), c("f1", "f2")))
  out <- mo_apply_matching(list(rna = m1), method = "exact")
  expect_identical(out$mats$rna, m1)
  expect_equal(out$dropped, list())
})

test_that("mo_merge_sample_meta() unions sample IDs from both frames and lets 'b' win on a column-name collision", {
  a <- data.frame(sex = c("F", "M"), group = c("HC", "HC"), row.names = c("S1", "S2"))
  b <- data.frame(sex = c("Female", "Male", "Female"), row.names = c("S1", "S2", "S3"))
  out <- mo_merge_sample_meta(a, b)
  expect_setequal(rownames(out), c("S1", "S2", "S3"))
  expect_equal(out["S1", "sex"], "Female")
  expect_equal(out["S1", "group"], "HC")
  expect_true(is.na(out["S3", "group"]))
})

test_that("mo_merge_sample_meta() passes through when either side is NULL/empty", {
  a <- data.frame(sex = "F", row.names = "S1")
  expect_identical(mo_merge_sample_meta(NULL, a), a)
  expect_identical(mo_merge_sample_meta(a, NULL), a)
  expect_identical(mo_merge_sample_meta(a, data.frame()), a)
})

test_that("mo_read_meta_file() reads a metadata CSV keyed by its first column as rownames", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("sample_id,sex,age", "S1,F,40", "S2,M,50"), path)
  out <- mo_read_meta_file(list(datapath = path))
  expect_setequal(rownames(out), c("S1", "S2"))
  expect_equal(out["S1", "sex"], "F")
})

test_that("mo_read_meta_file() returns NULL for a NULL fileInput value or an unreadable/empty file", {
  expect_null(mo_read_meta_file(NULL))
  path <- tempfile(fileext = ".csv")
  writeLines(character(0), path)
  expect_null(mo_read_meta_file(list(datapath = path)))
})
