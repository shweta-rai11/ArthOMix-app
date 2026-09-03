## Regression guard for a gap found in the transcriptomics audit (2026-08-26):
## mod_dataset.R's upload message promises that duplicated feature identifiers
## keep "only the first occurrence" downstream, but the single-dataset branch

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "functions", "expression_type.R"))
source_from_app_root(file.path("R", "transcriptomics", "01_Data", "mod_dataset.R"))
source_from_app_root(file.path("R", "transcriptomics", "03_Preprocessing_Batch_Correction", "mod_preprocessing_explore.R"))
source_from_app_root(file.path("R", "transcriptomics", "03_Preprocessing_Batch_Correction", "mod_preprocessing.R"))

test_that("a single loaded dataset with a duplicate feature ID is deduplicated, keeping the first occurrence", {
  d0 <- load_default_dataset()
  expr <- d0$expr[1:50, 1:20]
  rownames(expr)[10] <- rownames(expr)[3]
  expr[10, ] <- expr[3, ] + 0.5
  meta <- d0$meta[1:20, c("sample", "group")]
  dataset <- shiny::reactiveValues(expr = expr, meta = meta, source = "Uploaded dataset: dup_test.csv", source_type = "uploaded")

  shiny::testServer(mod_preprocessing_server, args = list(id = "pp", dataset = dataset, results = shiny::reactiveValues()), {
    session$setInputs(preloaded_selected = "__current__", preloaded_log2 = "auto")
    session$setInputs(preloaded_run = 1)
    session$setInputs(merge_mode = "own")
    session$setInputs(merge_btn = 1)
    m <- merged()

    expect_equal(nrow(m$expr), nrow(expr) - 1L)
    expect_equal(as.numeric(m$expr[rownames(expr)[3], ]), as.numeric(expr[3, ]))
  })
})
