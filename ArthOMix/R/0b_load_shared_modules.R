## R/0b_load_shared_modules.R
## Sources R/shared/*.R the same way 0_load_omics_modules.R sources
## transcriptomics/methylomics/crossomics/multiomics and 0a_load_auth_modules.R
## sources R/auth/ - shiny's own auto-loader (shiny:::loadSupport()) only
## scans R/*.R non-recursively, so R/shared/ is invisible to it unless
## sourced explicitly here. R/shared/ holds modules that are genuinely used
## across more than one of the four omics verticals (currently just
## mod_arthochat.R, the AI assistant, whose context builders span
## transcriptomics/methylomics/crossomics/multiomics) rather than belonging
## to any single one of them.
##
## Named "0b_" so C-locale sort_c() runs it right after
## 0a_load_auth_modules.R, before ui_shell.R/ui.R/server.R need
## mod_arthochat_ui()/mod_arthochat_server() defined.
.shared_files <- sort(list.files(file.path("R", "shared"), pattern = "\\.[rR]$", full.names = TRUE))
for (.shared_file in .shared_files) source(.shared_file, local = TRUE)
rm(.shared_files, .shared_file)
