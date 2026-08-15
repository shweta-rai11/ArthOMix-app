## R/methylomics/mod_methyl_celltype.R
## Submodule: Cell-Type Deconvolution (script02_celltype, METHODS 2.Y)
## Not built yet - scaffold only, so the tile exists in the Methylomics >
## Sub-modules grid. Will reproduce EpiDISH (Houseman constrained-projection,
## Reinius 7-cell-type reference) whole-blood cell-type composition
## estimates from the sex-stratified GSE42861 analysis, plus a live
## estimate on the uploaded/QC'd dataset.

mod_methyl_celltype_config <- list(
  id = "celltype", title = "Cell-Type Deconvolution", icon = "people-group", group = "Data",
  description = "Whole-blood cell-type composition (EpiDISH, Houseman constrained-projection, Reinius 7-cell-type reference): the default sex-stratified estimates from the preloaded whole-blood dataset, plus a live estimate on your own uploaded/QC'd dataset."
)

mod_methyl_celltype_ui <- function(id) {
  ns <- NS(id)
  box(
    width = 12, status = "primary", solidHeader = FALSE,
    div(
      class = "coming-soon",
      icon("hammer", class = "coming-soon-icon"),
      h4("This module is not built yet"),
      p(mod_methyl_celltype_config$description)
    )
  )
}

mod_methyl_celltype_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    NULL
  })
}
