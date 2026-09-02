## R/0_load_omics_modules.R
## Transcriptomics, Methylomics, Cross-Omics, and Multiomics each keep their
## mod_*.R files in their own subfolder (R/transcriptomics/, R/methylomics/,
## R/crossomics/, R/multiomics/) rather than one flat R/ directory - and,
## within a vertical, in publication-oriented numbered stage folders plus a
## functions/ folder for helpers shared across stages (e.g.
## R/crossomics/02_Expression_Methylation_Integration/,
## R/crossomics/functions/integration/ - see each vertical's README.md).
## shiny's own auto-loader (shiny:::loadSupport()) only scans R/*.R
## non-recursively - list.files(helpersDir, pattern = "\\.[rR]$",
## recursive = FALSE) - so any file below R/ that isn't directly in R/ itself
## is invisible to it unless sourced explicitly, which is all this file does.
## `recursive = TRUE` below walks the numbered/functions subfolders the same
## way it already walked the one-level vertical split.
##
## Named with a leading digit so C-locale sort_c() (loadSupport's own sort,
## same as every other R/*.R file) runs this before submodules_registry.R,
## which needs every mod_*_config/_ui/_server already defined, and before
## ui_shell.R.
##
## Sourcing order *within* a vertical (numbered-folder order vs. functions/
## vs. alphabetical) does not matter: every file here only *defines*
## functions/lists at its own top level (module config/ui/server, helpers).
## The one place that *reads* those bindings at source time rather than at
## Shiny runtime is submodules_registry.R's TX_MODULES/MX_MODULES/etc lists -
## and that file is guaranteed to load after every file below regardless of
## this loop's internal order, since "submodules_registry.R" sorts after
## "0_load_omics_modules.R" at the R/ top level. Cross-file *calls* (e.g. a
## dmr module calling a helper physically defined in the dmp module) happen
## inside moduleServer() bodies, which don't run until server.R instantiates
## them - long after all of R/ has been sourced - so they're order-independent
## too.
##
## source(..., local = TRUE) at top level of a script itself being sourced
## with envir = renv (loadSupport's shared child environment - see
## submodules_registry.R's own header comment) evaluates each subfolder
## file in that same renv, exactly where it would land if this were still
## one flat directory - ui.R/server.R/submodules_registry.R see no
## difference from before the split.
for (.omics_dir in c("transcriptomics", "methylomics", "crossomics", "multiomics")) {
  .omics_files <- sort(list.files(file.path("R", .omics_dir), pattern = "\\.[rR]$", full.names = TRUE, recursive = TRUE))
  for (.omics_file in .omics_files) source(.omics_file, local = TRUE)
}
rm(.omics_dir, .omics_files, .omics_file)
