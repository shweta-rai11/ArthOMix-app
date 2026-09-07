## Module 3 (Multiomics) - same-patient ID validation in the Dataset Workspace.
## Matched-sample Multi-Omics integration must only proceed when every layer
## carries exactly the same set of patient IDs (order-insensitive); mismatches,
## missing IDs and duplicates block Preprocessing and Activation with the
## required message instead of silently intersecting the datasets.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
for (f in list.files(file.path(app_dir, "R", "multiomics"), pattern = "[.]R$", recursive = TRUE, full.names = TRUE)) {
  suppressWarnings(suppressMessages(source(f)))
}

sp_genes <- c("TNF", "IL6", "CXCL8", "STAT1", "IRF7", "MX1", "ISG15", "IFI44L", "OAS1", "JAK2")
sp_cpgs <- sprintf("cg%08d", 1:10)

sp_write_expr <- function(ids, file, id_col = TRUE) {
  set.seed(11)
  df <- data.frame(matrix(round(rnorm(length(ids) * 10, 8, 2), 3), nrow = length(ids)))
  colnames(df) <- sp_genes
  if (id_col) df <- cbind(Patient_ID = ids, df)
  utils::write.csv(df, file, row.names = FALSE)
  file
}
sp_write_meth <- function(ids, file) {
  set.seed(12)
  df <- data.frame(Patient_ID = ids, matrix(round(runif(length(ids) * 10), 3), nrow = length(ids)), check.names = FALSE)
  colnames(df) <- c("Patient_ID", sp_cpgs)
  utils::write.csv(df, file, row.names = FALSE)
  file
}

## Drives the real Dataset Workspace server: upload both files, Validate, then
## attempt Preprocessing and Activation exactly as a user would.
sp_drive <- function(expr_ids, meth_ids, expr_has_id_col = TRUE, meta = NULL, method = "exact", patient_col = NULL) {
  ef <- sp_write_expr(expr_ids, tempfile(fileext = ".csv"), id_col = expr_has_id_col)
  mf <- sp_write_meth(meth_ids, tempfile(fileext = ".csv"))
  meta_file <- if (!is.null(meta)) { p <- tempfile(fileext = ".csv"); utils::write.csv(meta, p, row.names = FALSE); p } else NULL
  multi_dataset <- shiny::reactiveValues()
  multi_results <- shiny::reactiveValues()
  result <- NULL
  shiny::testServer(mod_multi_dataset_server, args = list(id = "ds", multi_dataset = multi_dataset, multi_results = multi_results), {
    session$setInputs(dataset_source = "upload")
    gen <- block_reset_gen()
    session$setInputs(ublock1_label = "Transcriptomics", ublock1_type = "rnaseq", ublock1_table_shape = "wide", ublock1_orient = "samples_rows",
                      ublock2_label = "Methylomics", ublock2_type = "methylation", ublock2_table_shape = "wide", ublock2_orient = "samples_rows")
    do.call(session$setInputs, stats::setNames(list(data.frame(name = basename(ef), datapath = ef)), paste0("ublock1_file_g", gen)))
    do.call(session$setInputs, stats::setNames(list(data.frame(name = basename(mf), datapath = mf)), paste0("ublock2_file_g", gen)))
    if (!is.null(meta_file)) do.call(session$setInputs, stats::setNames(list(data.frame(name = basename(meta_file), datapath = meta_file)), paste0("ublock1_meta_file_g", gen)))
    session$setInputs(validate_btn = 1)
    session$setInputs(matching_method = method)
    if (!is.null(patient_col)) session$setInputs(patient_id_col = patient_col)
    sp <- same_patient()
    cmp <- compat()
    session$setInputs(impute_method = "none", max_sample_missing = 100, max_feature_missing = 100,
                      norm_Transcriptomics = "none", norm_Methylomics = "none", preprocess_btn = 1)
    session$setInputs(active_layers = c("Transcriptomics", "Methylomics"), activate_btn = 1)
    result <<- list(
      layers = names(raw$mats), check = sp, overall = cmp$overall_label, sample_matching_ok = cmp$sample_matching_ok,
      preprocessed = !is.null(proc$scaled_mats),
      activated = isTRUE(multi_dataset$active) && length(multi_dataset$layers) == 2,
      activated_ids = if (length(multi_dataset$layers) > 0) rownames(multi_dataset$layers[[1]]) else NULL
    )
  })
  result
}

expect_blocked <- function(r, category) {
  expect_setequal(r$layers, c("Transcriptomics", "Methylomics"))
  expect_false(isTRUE(r$check$ok))
  expect_equal(r$check$category, category)
  expect_false(r$sample_matching_ok)
  expect_equal(r$overall, "NOT READY - patient IDs are not the same across datasets")
  expect_false(r$preprocessed)
  expect_false(r$activated)
}
expect_allowed <- function(r, ids) {
  expect_setequal(r$layers, c("Transcriptomics", "Methylomics"))
  expect_true(r$check$ok)
  expect_equal(r$overall, "READY")
  expect_true(r$preprocessed)
  expect_true(r$activated)
  expect_setequal(r$activated_ids, ids)
}

test_that("multi_live_same_patient_check(): identical ID sets pass regardless of order, numeric vs character, case or whitespace", {
  mk <- function(ids) { m <- matrix(1, nrow = length(ids), ncol = 2); rownames(m) <- as.character(ids); m }
  ov <- multi_live_sample_overlap(list(Transcriptomics = mk(c("P001", "P002", "P003")), Methylomics = mk(c("P003", "P001", "P002"))))
  expect_true(multi_live_same_patient_check(ov)$ok)
  expect_equal(multi_live_same_patient_check(ov)$n_patients, 3)
  ov <- multi_live_sample_overlap(list(Transcriptomics = mk(c(1, 2, 3)), Methylomics = mk(c("1", "2", "3"))))
  expect_true(multi_live_same_patient_check(ov)$ok)
  ov <- multi_live_sample_overlap(list(Transcriptomics = mk(c("P001", "P002", "P003")), Methylomics = mk(c("p001 ", "P002", "p003"))))
  expect_true(multi_live_same_patient_check(ov)$ok)
})

test_that("multi_live_same_patient_check(): completely different, partially overlapping and extra-patient sets fail with the exact required message", {
  mk <- function(ids) { m <- matrix(1, nrow = length(ids), ncol = 2); rownames(m) <- as.character(ids); m }
  for (meth in list(c("P004", "P005", "P006"), c("P001", "P002", "P004"), c("P001", "P002", "P003", "P004"), c("P001", "P002"))) {
    ov <- multi_live_sample_overlap(list(Transcriptomics = mk(c("P001", "P002", "P003")), Methylomics = mk(meth)))
    chk <- multi_live_same_patient_check(ov)
    expect_false(chk$ok)
    expect_equal(chk$category, "mismatch")
    expect_equal(chk$message, MULTI_SAME_PATIENT_MESSAGE)
  }
  expect_equal(MULTI_SAME_PATIENT_MESSAGE,
               "Patient ID should be same, this requires expression and methylation data to be from the same patient ID, you can explore cross-omics if your Patient ID are different.")
})

test_that("multi_live_same_patient_check(): duplicate, blank and unassignable patient IDs and empty datasets are detected rather than treated as valid", {
  mk <- function(ids) { m <- matrix(1, nrow = length(ids), ncol = 2); rownames(m) <- as.character(ids); m }
  ov <- multi_live_sample_overlap(list(Transcriptomics = mk(c("P001", "P002", "P003")), Methylomics = mk(c("P001", "P001", "P003"))))
  chk <- multi_live_same_patient_check(ov)
  expect_false(chk$ok); expect_equal(chk$category, "duplicate"); expect_match(chk$detail, "Methylomics has 1 duplicated patient ID")
  ov <- multi_live_sample_overlap(list(Transcriptomics = mk(c("P001", "P002", "P003")), Methylomics = mk(c("P001", "", "P003"))))
  chk <- multi_live_same_patient_check(ov)
  expect_false(chk$ok); expect_equal(chk$category, "missing")
  ov <- multi_live_sample_overlap(list(Transcriptomics = mk(c("P001", "P002", "P003")), Methylomics = mk(c("P001", "P002", "P003"))))
  chk <- multi_live_same_patient_check(ov, dropped = list(Transcriptomics = 0L, Methylomics = 1L))
  expect_false(chk$ok); expect_equal(chk$category, "missing")
  ov <- multi_live_sample_overlap(list(Transcriptomics = mk(c("P001", "P002", "P003")), Methylomics = mk(character(0))))
  chk <- multi_live_same_patient_check(ov)
  expect_false(chk$ok); expect_equal(chk$category, "missing")
  expect_true(multi_live_same_patient_check(NULL)$pending)
})

test_that("multi_dataset_compatibility() reports patient-ID mismatch as NOT READY and clears sample_matching_ok", {
  v <- list(ok = TRUE, layer = "A", n_samples = 5, n_features = 10, pct_missing = 0, n_duplicate_samples = 0, n_duplicate_features = 0, n_zero_variance = 0)
  ov <- list(ok = TRUE, n_shared = 4, per_layer = list(A = 5, B = 4))
  fail <- list(ok = FALSE, pending = FALSE, category = "mismatch")
  cmp <- multi_dataset_compatibility(list(A = v, B = v), overlap = ov, same_patient = fail)
  expect_false(cmp$sample_matching_ok)
  expect_true(cmp$same_patient_failed)
  expect_equal(cmp$overall_label, "NOT READY - patient IDs are not the same across datasets")
  cmp_ok <- multi_dataset_compatibility(list(A = v, B = v), overlap = ov, same_patient = list(ok = TRUE, pending = FALSE))
  expect_true(cmp_ok$sample_matching_ok)
  expect_equal(cmp_ok$overall_label, "READY")
  expect_equal(multi_dataset_compatibility(list(A = v, B = v), overlap = ov)$overall_label, "READY")
})

test_that("Dataset Workspace: same patients in the same order pass and the existing workflow proceeds to activation", {
  expect_allowed(sp_drive(c("P001", "P002", "P003"), c("P001", "P002", "P003")), c("P001", "P002", "P003"))
})

test_that("Dataset Workspace: same patients in a different order pass (set-based, not row-position matching)", {
  expect_allowed(sp_drive(c("P001", "P002", "P003"), c("P003", "P001", "P002")), c("P001", "P002", "P003"))
})

test_that("Dataset Workspace: completely different patient IDs are blocked with the required message", {
  r <- sp_drive(c("P001", "P002", "P003"), c("P004", "P005", "P006"))
  expect_blocked(r, "mismatch")
  expect_equal(r$check$message, MULTI_SAME_PATIENT_MESSAGE)
})

test_that("Dataset Workspace: partially overlapping and extra-patient ID sets are blocked, never silently intersected", {
  r <- sp_drive(c("P001", "P002", "P003"), c("P001", "P002", "P004"))
  expect_blocked(r, "mismatch"); expect_equal(r$check$message, MULTI_SAME_PATIENT_MESSAGE)
  r <- sp_drive(c("P001", "P002", "P003", "P004"), c("P001", "P002", "P005", "P006"))
  expect_blocked(r, "mismatch"); expect_equal(r$check$message, MULTI_SAME_PATIENT_MESSAGE)
  r <- sp_drive(c("P001", "P002", "P003", "P004"), c("P001", "P002", "P003"))
  expect_blocked(r, "mismatch")
})

test_that("Dataset Workspace: a missing patient ID (blank ID row, or no ID column at all) blocks integration", {
  expect_blocked(sp_drive(c("P001", "P002", "P003"), c("P001", "", "P003")), "missing")
  expect_blocked(sp_drive(c("P001", "P002", "P003"), c("P001", "P002", "P003"), expr_has_id_col = FALSE), "mismatch")
})

test_that("Dataset Workspace: duplicate patient IDs in either dataset block integration", {
  expect_blocked(sp_drive(c("P001", "P002", "P003"), c("P001", "P001", "P003")), "duplicate")
  expect_blocked(sp_drive(c("P001", "P002", "P002"), c("P001", "P002", "P003")), "duplicate")
})

test_that("Dataset Workspace: numeric IDs in one file and character IDs in the other are the same patients", {
  expect_allowed(sp_drive(c(1, 2, 3), c("1", "2", "3")), c("1", "2", "3"))
})

test_that("Dataset Workspace: the metadata patient-ID matching route (GEO-style differing sample IDs) still passes when every sample maps to the same patients, and fails when one cannot be assigned", {
  meta <- data.frame(sample_id = c("GSM1", "GSM2", "GSM3", "GSM4", "GSM5", "GSM6"), patient = c("P001", "P002", "P003", "P003", "P001", "P002"))
  expect_allowed(sp_drive(c("GSM1", "GSM2", "GSM3"), c("GSM4", "GSM5", "GSM6"), meta = meta, method = "patient_id", patient_col = "patient"), c("P001", "P002", "P003"))
  expect_blocked(sp_drive(c("GSM1", "GSM2", "GSM3"), c("GSM4", "GSM5", "GSM6"), meta = meta[meta$sample_id != "GSM6", ], method = "patient_id", patient_col = "patient"), "missing")
})

test_that("Dataset Workspace: the preloaded reference cell (identical sample sets) is unaffected by the gate", {
  skip_if_not(MULTI_DATA_AVAILABLE, "multi-omics preloaded data not available in this deployment")
  multi_dataset <- shiny::reactiveValues()
  multi_results <- shiny::reactiveValues()
  shiny::testServer(mod_multi_dataset_server, args = list(id = "ds", multi_dataset = multi_dataset, multi_results = multi_results), {
    session$setInputs(dataset_source = "preloaded")
    session$setInputs(preloaded_pick = "ra_antitnf", preloaded_cell = "female_Etanercept")
    session$setInputs(load_preloaded_btn = 1)
    expect_true(same_patient()$ok)
    expect_equal(compat()$overall_label, "READY")
    session$setInputs(impute_method = "none", max_sample_missing = 100, max_feature_missing = 100, norm_Transcriptomics = "none", norm_Methylomics = "none")
    session$setInputs(preprocess_btn = 1)
    expect_false(is.null(proc$scaled_mats))
    session$setInputs(active_layers = c("Transcriptomics", "Methylomics"), activate_btn = 1)
    expect_equal(length(multi_dataset$layers), 2)
  })
})
