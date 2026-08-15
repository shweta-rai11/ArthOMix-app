## R/methylomics/mod_methyl_wgcna.R
## Submodule: WGCNA / Co-Methylation Network (script05_wgcna_sexstratified,
## METHODS 2.CC). Not built yet - scaffold only, so the tile exists in the
## Methylomics > Sub-modules grid. Will reproduce the sex-stratified signed
## weighted co-methylation network (WGCNA on the broad, correlation-
## structure-preserving CpG set, not the DMP/DMR panels) from the default
## GSE42861 analysis, plus a live network build on the uploaded dataset.

mod_methyl_wgcna_config <- list(
  id = "wgcna", title = "WGCNA (Co-Methylation Network)", icon = "circle-nodes", group = "Network",
  description = "Sex-stratified signed weighted co-methylation network analysis: the default modules from the preloaded whole-blood dataset reproduced from their own results, plus a live network build on your own uploaded/QC'd/normalized dataset."
)

mod_methyl_wgcna_ui <- function(id) {
  ns <- NS(id)
  box(
    width = 12, status = "primary", solidHeader = FALSE,
    div(
      class = "coming-soon",
      icon("hammer", class = "coming-soon-icon"),
      h4("This module is not built yet"),
      p(mod_methyl_wgcna_config$description)
    )
  )
}

mod_methyl_wgcna_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    NULL
  })
}
