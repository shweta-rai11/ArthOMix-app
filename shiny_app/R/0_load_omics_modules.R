## R/0_load_omics_modules.R
## Transcriptomics, Methylomics, and Cross-Omics each keep their mod_*.R
## files in their own subfolder (R/transcriptomics/, R/methylomics/,
## R/crossomics/) rather than one flat R/ directory. shiny's own auto-loader
## (shiny:::loadSupport()) only
## scans R/*.R non-recursively - list.files(helpersDir, pattern = "\\.[rR]$",
## recursive = FALSE) - so files in either subfolder are invisible to it
## unless sourced explicitly, which is all this file does.
##
## Named with a leading digit so C-locale sort_c() (loadSupport's own sort,
## same as every other R/*.R file) runs this before submodules_registry.R,
## which needs every mod_*_config/_ui/_server already defined, and before
## ui_shell.R.
##
## source(..., local = TRUE) at top level of a script itself being sourced
## with envir = renv (loadSupport's shared child environment - see
## submodules_registry.R's own header comment) evaluates each subfolder
## file in that same renv, exactly where it would land if this were still
## one flat directory - ui.R/server.R/submodules_registry.R see no
## difference from before the split.
for (.omics_dir in c("transcriptomics", "methylomics", "crossomics")) {
  .omics_files <- sort(list.files(file.path("R", .omics_dir), pattern = "\\.[rR]$", full.names = TRUE))
  for (.omics_file in .omics_files) source(.omics_file, local = TRUE)
}
rm(.omics_dir, .omics_files, .omics_file)
