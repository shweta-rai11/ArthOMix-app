## R/0_load_omics_modules.R
## Transcriptomics, Methylomics, Cross-Omics, and Multiomics each keep their
## mod_*.R files in their own subfolder (R/transcriptomics/, R/methylomics/,
for (.omics_dir in c("transcriptomics", "methylomics", "crossomics", "multiomics")) {
  .omics_files <- sort(list.files(file.path("R", .omics_dir), pattern = "\\.[rR]$", full.names = TRUE, recursive = TRUE))
  for (.omics_file in .omics_files) source(.omics_file, local = TRUE)
}
rm(.omics_dir, .omics_files, .omics_file)
