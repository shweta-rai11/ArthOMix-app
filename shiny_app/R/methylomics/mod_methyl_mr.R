## R/methylomics/mod_methyl_mr.R
## Submodule: Mendelian Randomization (script08_mendelian_randomization,
## 08b/08c, METHODS 2.FF). Not built yet - scaffold only, so the tile
## exists in the Methylomics > Sub-modules grid. Will reproduce the
## sex-stratified MR of the majority-vote CpG panel against RA (method
## tiered by surviving cis-mQTL instrument count: Wald ratio / IVW-only /
## full IVW+Egger+weighted-median+mode suite) from the default GoDMC/RA-GWAS
## analysis.

mod_methyl_mr_config <- list(
  id = "mr", title = "Mendelian Randomization", icon = "route", group = "Genetics",
  description = "Sex-stratified Mendelian randomization of the majority-vote CpG panel against RA, using harmonised cis-mQTL (GoDMC) / RA-GWAS instruments: method tiered by surviving instrument count (Wald ratio, IVW-only, or the full IVW+Egger+weighted-median+mode suite with sensitivity analyses)."
)

mod_methyl_mr_ui <- function(id) {
  ns <- NS(id)
  box(
    width = 12, status = "primary", solidHeader = FALSE,
    div(
      class = "coming-soon",
      icon("hammer", class = "coming-soon-icon"),
      h4("This module is not built yet"),
      p(mod_methyl_mr_config$description)
    )
  )
}

mod_methyl_mr_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    NULL
  })
}
