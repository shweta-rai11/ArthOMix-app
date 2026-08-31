## Regression guard for a gap found in the transcriptomics audit (2026-08-26):
## mod_pp_source_ui()/mod_pp_source_server() - a fully-written per-dataset
## upload panel (fileInput + column mapping + filters), instantiated
## MAX_PP_SOURCES times as `pp_sources` - was never mounted by any renderUI(),
## so genuine multi-file upload through the Preprocessing tab itself was
## unreachable; the only way a user's own data reached Preprocessing was
## indirectly, via the single upload staged once on the Dataset tab. It is
## now mounted as an "Upload Your Own Data" box, and its results feed into
## merge_inputs() alongside the bundled-cohort picker.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "mod_dataset.R"))
source_from_app_root(file.path("R", "transcriptomics", "mod_preprocessing_explore.R"))
source_from_app_root(file.path("R", "transcriptomics", "mod_preprocessing.R"))

pp_write_fixture <- function(dir, seed, offset, sample_prefix) {
  set.seed(seed)
  genes <- paste0("GENE", 1:80)
  m <- matrix(rnorm(80 * 10, mean = 8 + offset, sd = 1), 80, 10,
              dimnames = list(genes, paste0(sample_prefix, 1:10)))
  expr_path <- file.path(dir, paste0("expr_", sample_prefix, ".csv"))
  meta_path <- file.path(dir, paste0("meta_", sample_prefix, ".csv"))
  write.csv(data.frame(gene = rownames(m), m, check.names = FALSE), expr_path, row.names = FALSE)
  write.csv(data.frame(sample = colnames(m), group = rep(c("HC", "RA"), 5)), meta_path, row.names = FALSE)
  list(expr_path = expr_path, meta_path = meta_path)
}

pp_mkfile <- function(path) {
  data.frame(name = basename(path), size = file.info(path)$size, type = "text/csv",
             datapath = path, stringsAsFactors = FALSE)
}

test_that("two independently uploaded datasets combine through the Preprocessing upload panel", {
  dir <- withr::local_tempdir()
  fx1 <- pp_write_fixture(dir, seed = 11, offset = 0, sample_prefix = "A")
  fx2 <- pp_write_fixture(dir, seed = 12, offset = 0.2, sample_prefix = "B")

  d0 <- load_default_dataset()
  dataset <- shiny::reactiveValues(expr = d0$expr, meta = d0$meta, source = d0$source, source_type = "preloaded")

  shiny::testServer(mod_preprocessing_server, args = list(id = "pp", dataset = dataset, results = shiny::reactiveValues()), {
    session$setInputs(n_sources = 2)
    session$setInputs(`src1-source_type` = "upload")
    session$setInputs(`src1-expr_file` = pp_mkfile(fx1$expr_path))
    session$setInputs(`src1-meta_file` = pp_mkfile(fx1$meta_path))
    session$setInputs(`src1-map_id` = "sample", `src1-map_group` = "group")
    session$setInputs(`src1-log2` = "skip", `src1-max_na_pct` = 0)
    session$setInputs(`src1-run` = 1)

    session$setInputs(`src2-source_type` = "upload")
    session$setInputs(`src2-expr_file` = pp_mkfile(fx2$expr_path))
    session$setInputs(`src2-meta_file` = pp_mkfile(fx2$meta_path))
    session$setInputs(`src2-map_id` = "sample", `src2-map_group` = "group")
    session$setInputs(`src2-log2` = "skip", `src2-max_na_pct` = 0)
    session$setInputs(`src2-run` = 1)

    own <- own_upload_results()
    expect_length(own, 2)

    session$setInputs(merge_mode = "own")
    session$setInputs(merge_selected = vapply(own, function(x) x$label, character(1)))
    session$setInputs(merge_btn = 1)
    m <- merged()
    expect_equal(ncol(m$expr), 20L)
    expect_true(grepl("Uploaded dataset:", m$sources))

    session$setInputs(color_by = "group", batch_col = "group", norm_method = "skip",
                       skip_combat = TRUE, mad_k = 3, min_pct = 0, variance_pct = 0)
    session$setInputs(run_btn = 1)
    session$setInputs(activate_btn = 1)
    expect_true(grepl("^Uploaded dataset", dataset$source))
  })
})
