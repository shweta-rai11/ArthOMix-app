## R/methylomics/mod_methyl_candidates.R
## Submodule: Candidate CpGs / Module-DMR Overlap (script06_module_dmr_venn,
## METHODS 2.DD). Not built yet - scaffold only, so the tile exists in the
## Methylomics > Sub-modules grid. Will reproduce the sex-stratified
## convergence check between WGCNA module membership and the DMR panel
## (Venn overlap) from the default GSE42861 analysis, plus a live overlap
## on the uploaded dataset's own WGCNA/DMR results.

mod_methyl_candidates_config <- list(
  id = "candidates", title = "Candidate CpGs (Module-DMR Overlap)", icon = "star", group = "Network",
  description = "Sex-stratified convergence check between WGCNA module membership and the DMR panel (Venn overlap): the default candidate CpGs from the preloaded whole-blood dataset, reproduced from their own results, plus a live overlap on your own uploaded dataset's WGCNA/DMR results."
)

mod_methyl_candidates_ui <- function(id) {
  ns <- NS(id)
  box(
    width = 12, status = "primary", solidHeader = FALSE,
    div(
      class = "coming-soon",
      icon("hammer", class = "coming-soon-icon"),
      h4("This module is not built yet"),
      p(mod_methyl_candidates_config$description)
    )
  )
}

mod_methyl_candidates_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    NULL
  })
}
