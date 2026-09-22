## R/ui_shell.R
## Shared app-shell UI: header, per-module sidebar, pipeline-summary timeline.

app_header <- function() {
  tagList(
    tags$div(
      class = "app-header",
      ## The mark sits in the brand slot, immediately left of the Home tab on the same row.
      tags$div(
        class = "app-header-brand",
        "ArthOMix"
      ),
      tags$div(
        class = "app-header-actions",
        if (ARTHOMIX_CHAT_ENABLED) tags$a(
          "Ask ArthOChat", href = "#", class = "btn btn-primary btn-sm",
          onclick = ARTHOCHAT_DRAWER_OPEN_JS
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
           // A sidebar marked data-server-nav decides its own active row (see omics_sidebar).
           $('.omics-sidebar:not([data-server-nav]) .sidebar-nav-item:visible').removeClass('active');
           $('.omics-sidebar:not([data-server-nav]) .sidebar-nav-item[data-match=\"' + shownText.replace(/\"/g, '') + '\"]:visible').addClass('active');
         });
       });"
    ))
  )
}

## server_nav = TRUE marks a sidebar whose active row is decided on the server (the step-grouped
## Transcriptomics one). app_header()'s shown.bs.tab handler leaves those alone, so it cannot strip
## the active class off a row the server just set - two rows there share one tab (Diagnostic Model
## appears in step 7 and, on its External Validation sub-tab, in step 8), which data-match cannot tell
## apart. A sidebar with no static nav_items renders only its dynamic list.
omics_sidebar <- function(module_id, module_label, nav_items, extra_sidebar_content = NULL,
                          dynamic_nav_output_id = NULL, server_nav = FALSE) {
  tags$div(
    class = "omics-sidebar",
    `data-server-nav` = if (isTRUE(server_nav)) "true" else NULL,
    tags$div(class = "omics-sidebar-heading", toupper(module_label)),
    if (length(nav_items)) tags$ul(
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
      uiOutput(dynamic_nav_output_id, container = function(...) tags$div(class = "sidebar-nav-dynamic", ...))
    },
    extra_sidebar_content
  )
}

## One row of a step-grouped sidebar. `key` is a TX_STEP_MAP key, so the two Diagnostic Model rows
## stay distinguishable; clicks go through one event input rather than per-row action buttons, which
## would fire spuriously every time the list re-renders with a new status.
omics_sidebar_step_row <- function(key, label, icon_name, status, reason, active = FALSE,
                                   input_id = "tx_nav_go") {
  tags$li(
    tags$a(
      href = "#",
      class = paste("sidebar-nav-item sidebar-step-item", paste0("sm-state-", status),
                    if (isTRUE(active)) "active"),
      title = reason,
      onclick = sprintf("Shiny.setInputValue('%s', '%s', {priority: 'event'}); return false;", input_id, key),
      icon(icon_name), tags$span(class = "sidebar-step-label", label),
      ## "Ready" is the resting state of every row, so saying it adds nothing. A tag appears only
      ## when there is something to report: a stored result, or one that needs re-running, or a
      ## dependency that is still missing.
      if (!identical(status, "ready")) {
        tags$span(class = paste("sm-state-pill", paste0("sm-state-", status)),
                  WF_STATUS_LABEL[[status]] %||% status)
      }
    )
  )
}

## Header for one workflow step in the sidebar, carrying that step's own rolled-up status.
omics_sidebar_step_header <- function(n, label, status) {
  tags$li(
    class = paste("sidebar-step-header", paste0("sm-state-", status)),
    tags$span(class = "sidebar-step-num", if (identical(status, "done")) icon("check") else n),
    tags$span(class = "sidebar-step-name", label)
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
    )
  )
}
