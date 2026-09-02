## Module 1 (Transcriptomics) - Cross-Tissue Validation's metadata-aware
## tissue classification (the "Tissue-type validation" section of
## mod_crosstissue.R). Pure functions only, tested independently of the
## server body: tissue_normalize/tissue_is_blood_derived, tissue_extract_
## field/classify_tissue_from_metadata/classify_tissue_string, ct_classify_
## training_tissue, and validate_cross_tissue's accept/reject rules.

source_from_app_root(file.path("R", "transcriptomics", "12_Cross_Tissue_Validation", "mod_crosstissue.R"))

## ---- tissue_normalize() -----------------------------------------------------

test_that("tissue_normalize() collapses common formatting variants to the same string", {
  expect_equal(tissue_normalize("Whole Blood"), "whole blood")
  expect_equal(tissue_normalize("whole_blood"), "whole blood")
  expect_equal(tissue_normalize("whole-blood"), "whole blood")
  expect_equal(tissue_normalize("  Whole   Blood  "), "whole blood")
})

test_that("tissue_normalize() strips a leading GEO-style 'label: value' prefix", {
  expect_equal(tissue_normalize("tissue: whole blood"), "whole blood")
  expect_equal(tissue_normalize("Sample type: Whole Blood"), "whole blood")
  expect_equal(tissue_normalize("source: PBMC"), "pbmc")
})

test_that("tissue_normalize() returns NA for NULL/NA/empty input", {
  expect_true(is.na(tissue_normalize(NULL)))
  expect_true(is.na(tissue_normalize(NA)))
  expect_true(is.na(tissue_normalize("")))
  expect_true(is.na(tissue_normalize("   ")))
})

## ---- tissue_is_blood_derived() ----------------------------------------------

test_that("tissue_is_blood_derived() recognizes every documented blood-derived term", {
  terms <- c("whole blood", "peripheral blood", "pbmc", "peripheral blood mononuclear cells",
             "buffy coat", "plasma", "serum", "venous blood", "arterial blood", "cord blood",
             "circulating blood", "blood derived", "blood")
  for (t in terms) expect_true(tissue_is_blood_derived(t), info = t)
})

test_that("tissue_is_blood_derived() recognizes blood terminology embedded in a longer phrase", {
  expect_true(tissue_is_blood_derived(tissue_normalize("sample obtained from peripheral blood")))
  expect_true(tissue_is_blood_derived(tissue_normalize("PBMC isolated from blood")))
  expect_true(tissue_is_blood_derived(tissue_normalize("whole blood RNA")))
  expect_true(tissue_is_blood_derived(tissue_normalize("blood-derived mononuclear cells")))
})

test_that("tissue_is_blood_derived() does not flag unrelated tissue names", {
  expect_false(tissue_is_blood_derived(tissue_normalize("synovium")))
  expect_false(tissue_is_blood_derived(tissue_normalize("synovial tissue")))
  expect_false(tissue_is_blood_derived(tissue_normalize("muscle")))
  expect_false(tissue_is_blood_derived(tissue_normalize("liver biopsy")))
})

## ---- classify_tissue_string() -----------------------------------------------

test_that("classify_tissue_string() classifies a blood-derived value", {
  out <- classify_tissue_string("Whole Blood", field = "tissue")
  expect_equal(out$classification, "blood")
  expect_equal(out$field, "tissue")
})

test_that("classify_tissue_string() classifies a non-blood value", {
  out <- classify_tissue_string("Synovium", field = "tissue")
  expect_equal(out$classification, "non-blood")
})

test_that("classify_tissue_string() returns unknown for NA/empty input", {
  expect_equal(classify_tissue_string(NA)$classification, "unknown")
  expect_equal(classify_tissue_string("")$classification, "unknown")
  expect_equal(classify_tissue_string(NULL)$classification, "unknown")
})

## ---- classify_tissue_from_metadata() / tissue_extract_field() --------------

test_that("classify_tissue_from_metadata() reads an exact 'tissue' column", {
  meta <- data.frame(sample = c("S1", "S2"), tissue = c("Whole Blood", "Whole Blood"))
  out <- classify_tissue_from_metadata(meta)
  expect_equal(out$field, "tissue")
  expect_equal(out$classification, "blood")
})

test_that("classify_tissue_from_metadata() reads a 'tissue_type' column", {
  meta <- data.frame(sample = c("S1", "S2"), tissue_type = c("Synovium", "Synovium"))
  out <- classify_tissue_from_metadata(meta)
  expect_equal(out$classification, "non-blood")
})

test_that("classify_tissue_from_metadata() reads a GEO-style 'characteristics_ch1' column", {
  meta <- data.frame(sample = c("S1", "S2"), characteristics_ch1 = c("tissue: whole blood", "tissue: whole blood"))
  out <- classify_tissue_from_metadata(meta)
  expect_equal(out$field, "characteristics_ch1")
  expect_equal(out$classification, "blood")
})

test_that("classify_tissue_from_metadata() reads a 'source_name_ch1' column", {
  meta <- data.frame(sample = c("S1", "S2"), source_name_ch1 = c("Peripheral blood", "Peripheral blood"))
  out <- classify_tissue_from_metadata(meta)
  expect_equal(out$classification, "blood")
})

test_that("classify_tissue_from_metadata() skips a same-pattern column whose value isn't tissue-related", {
  meta <- data.frame(sample = c("S1", "S2"), characteristics_ch1.3 = c("age: 45", "age: 52"))
  out <- classify_tissue_from_metadata(meta)
  expect_equal(out$classification, "unknown")
})

test_that("classify_tissue_from_metadata() returns unknown when no tissue-like column exists", {
  meta <- data.frame(sample = c("S1", "S2"), group = c("RA", "HC"), sex = c("F", "M"))
  out <- classify_tissue_from_metadata(meta)
  expect_equal(out$classification, "unknown")
  expect_true(is.na(out$field))
})

test_that("classify_tissue_from_metadata() returns unknown for a NULL/empty metadata table", {
  expect_equal(classify_tissue_from_metadata(NULL)$classification, "unknown")
  expect_equal(classify_tissue_from_metadata(data.frame())$classification, "unknown")
})

## ---- ct_classify_training_tissue() ------------------------------------------

fake_geo_sources <- list(
  list(gse = "GSE93272",  role = "Training (whole blood)"),
  list(gse = "GSE110169", role = "Training (whole blood)"),
  list(gse = "GSE89408",  role = "Validation (synovium)")
)

test_that("ct_classify_training_tissue() prefers a real metadata tissue column over the source label", {
  meta <- data.frame(sample = "S1", tissue = "Synovium")
  out <- ct_classify_training_tissue(meta, source = "Example dataset: sex-stratified RA blood cohort",
                                      source_type = "preloaded", geo_ids = c("GSE93272", "GSE110169"),
                                      geo_sources = fake_geo_sources)
  expect_equal(out$classification, "non-blood")
  expect_equal(out$field, "tissue")
})

test_that("ct_classify_training_tissue() falls back to the dataset source label for a preloaded pick with no tissue column", {
  meta <- data.frame(sample = "S1", group = "RA", sex = "F")
  out <- ct_classify_training_tissue(meta, source = "Example dataset: sex-stratified RA blood cohort (GSE93272 + GSE110169)",
                                      source_type = "preloaded", geo_ids = c("GSE93272", "GSE110169"),
                                      geo_sources = fake_geo_sources)
  expect_equal(out$classification, "blood")
})

test_that("ct_classify_training_tissue() falls back to the GEO_SOURCES registry when neither metadata nor source label carries tissue text", {
  meta <- data.frame(sample = "S1", group = "RA", sex = "F")
  out <- ct_classify_training_tissue(meta, source = "Individual dataset: GSE89408",
                                      source_type = "preloaded", geo_ids = "GSE89408",
                                      geo_sources = fake_geo_sources)
  expect_equal(out$classification, "non-blood")
})

test_that("ct_classify_training_tissue() never derives tissue from an uploaded dataset's filename-based source label", {
  meta <- data.frame(sample = "S1", group = "RA", sex = "F")  ## no tissue column
  out <- ct_classify_training_tissue(meta, source = "Uploaded dataset: expr.csv + meta.csv",
                                      source_type = "uploaded", geo_ids = character(0),
                                      geo_sources = fake_geo_sources)
  expect_equal(out$classification, "unknown")
})

test_that("ct_classify_training_tissue() reads a real metadata column for an uploaded/GEO source", {
  meta <- data.frame(sample = "S1", tissue = "Whole Blood")
  out <- ct_classify_training_tissue(meta, source = "Uploaded dataset: expr.csv + meta.csv",
                                      source_type = "uploaded", geo_ids = character(0),
                                      geo_sources = fake_geo_sources)
  expect_equal(out$classification, "blood")
})

## ---- validate_cross_tissue() ------------------------------------------------

blood <- classify_tissue_string("Whole Blood", field = "tissue")
pbmc <- classify_tissue_string("PBMC", field = "tissue")
plasma <- classify_tissue_string("Plasma", field = "tissue")
serum <- classify_tissue_string("Serum", field = "tissue")
peripheral_blood <- classify_tissue_string("Peripheral Blood", field = "tissue")
synovium <- classify_tissue_string("Synovium", field = "tissue")
synovial_tissue <- classify_tissue_string("Synovial tissue", field = "tissue")
muscle <- classify_tissue_string("Muscle", field = "tissue")
unknown_val <- classify_tissue_string(NA)

test_that("validate_cross_tissue() accepts blood training vs. a non-blood validation tissue", {
  expect_true(validate_cross_tissue(blood, synovium)$valid)
  expect_true(validate_cross_tissue(blood, synovial_tissue)$valid)
})

test_that("validate_cross_tissue() rejects blood training paired with any blood-derived validation sample", {
  expect_false(validate_cross_tissue(blood, peripheral_blood)$valid)
  expect_false(validate_cross_tissue(blood, pbmc)$valid)
  expect_false(validate_cross_tissue(blood, plasma)$valid)
  expect_false(validate_cross_tissue(blood, serum)$valid)
})

test_that("validate_cross_tissue() rejects identical non-blood tissue on both sides", {
  out <- validate_cross_tissue(synovium, synovium)
  expect_false(out$valid)
  expect_match(out$status, "not cross-tissue")
})

test_that("validate_cross_tissue() accepts two distinct non-blood tissues", {
  out <- validate_cross_tissue(synovium, muscle)
  expect_true(out$valid)
})

test_that("validate_cross_tissue() accepts a non-blood training tissue validated against blood", {
  out <- validate_cross_tissue(synovium, blood)
  expect_true(out$valid)
})

test_that("validate_cross_tissue() rejects (fails safely) when either side is unknown", {
  expect_false(validate_cross_tissue(unknown_val, synovium)$valid)
  expect_false(validate_cross_tissue(blood, unknown_val)$valid)
  expect_false(validate_cross_tissue(unknown_val, unknown_val)$valid)
})

test_that("validate_cross_tissue() reason text names both datasets on rejection", {
  out <- validate_cross_tissue(blood, pbmc)
  expect_match(out$reason, "blood-derived")
  expect_match(out$reason, "Whole Blood")
  expect_match(out$reason, "PBMC")
})
