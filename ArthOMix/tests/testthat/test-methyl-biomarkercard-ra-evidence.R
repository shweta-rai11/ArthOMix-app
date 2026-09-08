## Regression guard for a RED finding (2026-09-07 defense audit): the
## Methylomics Biomarker Card's "RA-associated biomarker" green banner (1) fired
## on the bare regex "rheumatoid|arthritis", matching any arthritis subtype
## (osteoarthritis, psoriatic, juvenile idiopathic, reactive) as if it were
## rheumatoid arthritis specifically, and (2) had no significance/effect-size/
## replication gate - a single unfiltered EWAS Catalog/Atlas hit was enough.

suppressWarnings(suppressMessages(
  source_from_app_root("global.R")
))
source_from_app_root(file.path("R", "methylomics", "functions", "annotation.R"))
source_from_app_root(file.path("R", "crossomics", "functions", "integration", "crossomics_integration_helpers.R"))
source_from_app_root(file.path("R", "methylomics", "15_Biomarker_Analysis", "mod_methyl_biomarkercard.R"))

test_that("bc_is_ra_trait matches rheumatoid arthritis but not other arthritis subtypes", {
  expect_true(bc_is_ra_trait("Rheumatoid arthritis"))
  expect_true(bc_is_ra_trait("rheumatoid arthritis (seropositive)"))
  expect_false(bc_is_ra_trait("Osteoarthritis"))
  expect_false(bc_is_ra_trait("Psoriatic arthritis"))
  expect_false(bc_is_ra_trait("Juvenile idiopathic arthritis"))
  expect_false(bc_is_ra_trait("Reactive arthritis"))
  expect_false(bc_is_ra_trait("Gouty arthritis"))
  expect_false(bc_is_ra_trait("Generic inflammatory disease"))
  expect_false(bc_is_ra_trait(NA_character_))
})

ra_ext <- function(ra_rows) list(ra_rows = ra_rows)

test_that("no rheumatoid match: neutral 'no evidence' banner, no table", {
  html <- fx_html_text(bc_section_ra_evidence(ra_ext(NULL)))
  expect_true(grepl("No rheumatoid-arthritis-specific evidence", html))
  expect_false(grepl("circle-check", html))
})

test_that("a single unreplicated hit with no p-value gets the weak/unreplicated banner, not the green 'associated' banner", {
  ra <- data.frame(source = "EWAS Atlas", trait = "Rheumatoid arthritis", effect = 0.1,
                    p = NA_character_, pmid = "12345678", stringsAsFactors = FALSE)
  html <- fx_html_text(bc_section_ra_evidence(ra_ext(ra)))
  expect_true(grepl("unreplicated", html, ignore.case = TRUE))
  expect_false(grepl("Rheumatoid arthritis-associated CpG", html))
})

test_that("a single hit with a significant p-value passes the gate and shows the associated banner", {
  ra <- data.frame(source = "EWAS Catalog", trait = "Rheumatoid arthritis", effect = 0.2,
                    p = "0.001", pmid = "12345678", stringsAsFactors = FALSE)
  html <- fx_html_text(bc_section_ra_evidence(ra_ext(ra)))
  expect_true(grepl("Rheumatoid arthritis-associated CpG", html))
  expect_true(grepl("p below 0.05", html))
})

test_that("two independent PMIDs with no p-value reach the gate via replication", {
  ra <- data.frame(source = c("EWAS Atlas", "EWAS Atlas"), trait = c("Rheumatoid arthritis", "Rheumatoid arthritis"),
                    effect = c(0.1, 0.15), p = c(NA_character_, NA_character_), pmid = c("111", "222"), stringsAsFactors = FALSE)
  html <- fx_html_text(bc_section_ra_evidence(ra_ext(ra)))
  expect_true(grepl("Rheumatoid arthritis-associated CpG", html))
  expect_true(grepl("2 independent studies", html))
})

test_that("an osteoarthritis-only row (would have matched the old 'arthritis' regex) never reaches this section as RA evidence", {
  ## ewas_combined() builds ra_rows via bc_is_ra_trait() before this function ever
  ## sees the data - reproduce that filtering step directly against a mixed trait set.
  combined <- data.frame(
    source = c("EWAS Catalog", "EWAS Catalog"),
    trait = c("Osteoarthritis", "Psoriatic arthritis"),
    effect = c(0.3, 0.4), p = c("0.001", "0.002"), pmid = c("999", "888"),
    stringsAsFactors = FALSE
  )
  ra_rows <- combined[bc_is_ra_trait(combined$trait), , drop = FALSE]
  expect_equal(nrow(ra_rows), 0L)
  html <- fx_html_text(bc_section_ra_evidence(ra_ext(if (nrow(ra_rows) > 0) ra_rows else NULL)))
  expect_true(grepl("No rheumatoid-arthritis-specific evidence", html))
})
