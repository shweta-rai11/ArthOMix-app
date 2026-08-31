## R/0a_load_auth_modules.R
## Sources R/auth/*.R the same way 0_load_omics_modules.R sources
## transcriptomics/methylomics/crossomics/multiomics - shiny's own
## auto-loader (shiny:::loadSupport()) only scans R/*.R non-recursively, so
## R/auth/ is invisible to it unless sourced explicitly here. Named "0a_" so
## C-locale sort_c() runs it right after 0_load_omics_modules.R, before
## ui_shell.R/ui.R/server.R need mod_auth_ui()/mod_auth_server() defined.
.auth_files <- sort(list.files(file.path("R", "auth"), pattern = "\\.[rR]$", full.names = TRUE))
for (.auth_file in .auth_files) source(.auth_file, local = TRUE)
rm(.auth_files, .auth_file)
