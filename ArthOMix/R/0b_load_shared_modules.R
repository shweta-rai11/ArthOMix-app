## R/0b_load_shared_modules.R
## Sources R/shared/*.R, same pattern as 0_load_omics_modules.R.
.shared_files <- sort(list.files(file.path("R", "shared"), pattern = "\\.[rR]$", full.names = TRUE))
for (.shared_file in .shared_files) source(.shared_file, local = TRUE)
rm(.shared_files, .shared_file)
