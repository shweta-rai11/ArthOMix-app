## R/ui_shell.R
## Shared app-shell UI components for the SaaS-dashboard redesign: header,
## per-module left sidebar, and the right-column pipeline-summary timeline.

app_header <- function(user_email = NULL) {
  tagList(
    tags$div(
      class = "app-header",
      tags$div(
        class = "app-header-brand",
        "ArthOMix"
      ),
      tags$div(
        class = "app-header-actions",
        tags$a(
          "Ask ArthOChat", href = "#", class = "btn btn-primary btn-sm",
          onclick = ARTHOCHAT_DRAWER_OPEN_JS
        ),
        actionButton("theme_toggle_btn", label = "", icon = icon("moon"),
                     class = "app-header-icon-btn", title = "Toggle light / dark mode"),
        tags$div(
          class = "app-header-account",
          tags$span(class = "app-header-avatar", title = "Account", icon("user")),
          tags$span(class = "app-header-account-email", user_email),
          actionButton("logout_btn", "Log Out", icon = icon("right-from-bracket"),
                       class = "btn btn-default btn-sm", title = "Log out")
        )
      )
    ),
    tags$script(HTML(
      "$(function(){
         $(document).on('shown.bs.tab', 'a[data-toggle=\"tab\"]', function(e){
           var shownText = $(e.target).text().trim();
           // Every module's sidebar stays in the DOM at once (Shiny keeps
           // every tabsetPanel pane mounted, just display:none-d when not
           // the active one) - data-match values like \"Sub-modules\"/
           // \"Dataset\" repeat across Transcriptomics/Methylomics/
           // Cross-Omics/Multi-Omics, so an unscoped match highlighted the
           // same-named item in every other module's sidebar too, not just
           // the module actually on screen. :visible (which correctly
           // accounts for a display:none ancestor) restricts both the
           // reset and the highlight to whichever sidebar is actually shown.
           $('.sidebar-nav-item:visible').removeClass('active');
           $('.sidebar-nav-item[data-match=\"' + shownText.replace(/\"/g, '') + '\"]:visible').addClass('active');
         });
       });"
    ))
  )
}

omics_sidebar <- function(module_id, module_label, nav_items, extra_sidebar_content = NULL, dynamic_nav_output_id = NULL) {
  tags$div(
    class = "omics-sidebar",
    tags$div(class = "omics-sidebar-heading", toupper(module_label)),
    tags$ul(
      class = "sidebar-nav",
      lapply(nav_items, function(it) {
        tags$li(
          tags$a(
            id = paste0("sidebar_nav_", module_id, "_", it$id), href = "#",
            class = "sidebar-nav-item action-button",
            `data-match` = it$match %||% it$label,
            icon(it$icon), it$label
          )
        )
      })
    ),
    if (!is.null(dynamic_nav_output_id)) {
      uiOutput(dynamic_nav_output_id, container = function(...) tags$ul(class = "sidebar-nav", ...))
    },
    tags$div(class = "omics-sidebar-heading", "QUICK LINKS"),
    tags$div(
      class = "sidebar-quicklinks",
      tags$a(href = "#", icon("book"), "Documentation"),
      tags$a(href = "#", icon("graduation-cap"), "Tutorials")
    ),
    extra_sidebar_content
  )
}

pipeline_summary_ui <- function(steps) {
  tags$div(
    class = "card",
    tags$div(class = "card-title", icon("list-check"), "Pipeline summary"),
    tags$ul(
      class = "pipeline-summary-list",
      lapply(steps, function(s) {
        tags$li(
          class = paste("pipeline-summary-step", paste0("state-", s$state)),
          tags$div(
            class = "pipeline-summary-dot",
            if (identical(s$state, "done")) icon("check") else as.character(s$number)
          ),
          tags$div(
            tags$div(class = "pipeline-summary-title", s$label),
            tags$div(class = "pipeline-summary-sub", s$sublabel)
          )
        )
      })
    ),
    tags$div(
      class = "pipeline-summary-footer",
      tags$a(href = "#", "View full pipeline docs ", icon("arrow-right"))
    )
  )
}
