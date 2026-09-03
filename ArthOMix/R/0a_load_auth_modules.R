## R/0a_load_auth_modules.R
## Sources R/auth/*.R the same way 0_load_omics_modules.R sources
## transcriptomics/methylomics/crossomics/multiomics - shiny's own
.auth_files <- sort(list.files(file.path("R", "auth"), pattern = "\\.[rR]$", full.names = TRUE))
for (.auth_file in .auth_files) source(.auth_file, local = TRUE)
rm(.auth_files, .auth_file)
