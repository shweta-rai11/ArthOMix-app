## Module 1 (Transcriptomics) - Cross-Ancestry MR Replication: the upload
## path's uploaded_outcome() called TwoSampleMR::format_data() unguarded.
## format_data() usually degrades gracefully (coerces bad columns to NA with
## a warning), but when every mapped column collapses onto the same
## non-SNP column - the realistic outcome of uploading a file with none of
## the expected column names, so every dropdown guess falls back to the
## file's first column - it throws a hard "SNP column not found" error that,
## unguarded, would crash the reactive chain instead of showing the app's
## normal red validate() banner. These tests cover the fix.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "transcriptomics", "13_Cross_Ancestral_Validation", "mod_crossancestry.R"))

ca_mkfile <- function(path) {
  data.frame(name = basename(path), size = file.info(path)$size, type = "text/csv",
             datapath = path, stringsAsFactors = FALSE)
}

test_that("a file with none of the expected GWAS columns (every field mapped to the same non-SNP column) fails with an actionable message, not a raw format_data crash", {
  bad_path <- tempfile(fileext = ".csv")
  write.csv(
    data.frame(sample_id = c("x1", "x2", "x3"), expression = c(1.1, 2.2, 3.3),
               log2fc = c(0.5, 0.6, 0.7), group = c("a", "b", "c")),
    bad_path, row.names = FALSE
  )

  dataset <- shiny::reactiveValues()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_crossancestry_server, args = list(id = "ca", dataset = dataset, results = results), {
    session$setInputs(data_source = "upload", arm_replaced = "eas",
                       gwas_file = ca_mkfile(bad_path), gwas_label = "Bad file")
    ## Every dropdown's auto-guess falls back to the file's first column,
    ## exactly as gwas_col_map_ui() does when none of GWAS_COL_PATTERNS match.
    session$setInputs(gwas_snp = "sample_id", gwas_beta = "sample_id", gwas_se = "sample_id",
                       gwas_pval = "sample_id", gwas_ea = "sample_id", gwas_oa = "sample_id", gwas_eaf = "")
    session$setInputs(run_female_btn = 1)

    err <- tryCatch(result_female(), error = function(e) e)
    expect_true(inherits(err, "shiny.silent.error"))
    expect_match(conditionMessage(err), "Could not parse the uploaded file", fixed = TRUE)
  })
})

test_that("an unreadable/empty upload is rejected before reaching format_data, with a distinct message", {
  empty_path <- tempfile(fileext = ".csv")
  writeLines(character(0), empty_path)

  dataset <- shiny::reactiveValues()
  results <- shiny::reactiveValues()
  shiny::testServer(mod_crossancestry_server, args = list(id = "ca", dataset = dataset, results = results), {
    session$setInputs(data_source = "upload", arm_replaced = "eas",
                       gwas_file = ca_mkfile(empty_path), gwas_label = "Empty file")
    session$setInputs(gwas_snp = "V1", gwas_beta = "V1", gwas_se = "V1",
                       gwas_pval = "V1", gwas_ea = "V1", gwas_oa = "V1", gwas_eaf = "")
    session$setInputs(run_female_btn = 1)

    err <- tryCatch(result_female(), error = function(e) e)
    expect_true(inherits(err, "shiny.silent.error"))
    expect_match(conditionMessage(err), "Could not read the uploaded GWAS file", fixed = TRUE)
  })
})
