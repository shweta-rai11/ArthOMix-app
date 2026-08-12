## R/ui_shell.R
## Shared app-shell UI components for the SaaS-dashboard redesign: header,
## per-module left sidebar, and the right-column pipeline-summary timeline.
## These are presentational only - they call into existing inputs/outputs
## (arthochat_shortcut_ui() from global.R, the "header_search_submit" and
## "sidebar_nav_*" inputs wired up as new, additive observers in server.R)
## rather than defining any new reactive/analysis logic. Written so
## omics_sidebar() can be reused for Methylomics/Cross-Omics/Multi-Omics
## once those modules exist, not just Transcriptomics.

## ---------------------------------------------------------------------------
## Header: brand, search (jumps to a module/sub-module on Enter via the
## "header_search_submit" input, handled in server.R), Ask ArthOChat
## (same client-side tab-show trick arthochat_shortcut_ui() already uses),
## a light-mode-only theme toggle (server.R just shows a notification), and
## a decorative avatar (no auth system exists in this app).
## ---------------------------------------------------------------------------

app_header <- function() {
  tagList(
    tags$div(
      class = "app-header",
      tags$div(
        class = "app-header-brand",
        tags$span(class = "app-header-brand-icon", icon("share-nodes")),
        "ArthOMix Explorer"
      ),
      tags$div(
        class = "app-header-search",
        tags$span(class = "app-header-search-icon", icon("magnifying-glass")),
        tags$input(
          type = "text", id = "header_search", class = "form-control",
          placeholder = "Search modules, datasets, genes, etc...", autocomplete = "off"
        ),
        tags$span(class = "app-header-search-kbd", "Enter")
      ),
      tags$div(
        class = "app-header-actions",
        ## Always visible, on every page - which dataset every module is
        ## currently reading from (dataset$source, set once at app start by
        ## global.R::load_default_dataset(), changed only by mod_dataset.R's
        ## upload or mod_preprocessing.R's "Use this as the active dataset").
        ## Exists specifically so switching between the app's own preloaded
        ## example and your own uploaded/merged/corrected data is never
        ## ambiguous - no need to go check the Dataset tab to know which one
        ## every other module is about to use.
        uiOutput("active_dataset_badge", inline = TRUE),
        tags$a(
          "Ask ArthOChat", href = "#", class = "btn btn-primary btn-sm",
          onclick = ARTHOCHAT_DRAWER_OPEN_JS
        ),
        actionButton("theme_toggle_btn", label = "", icon = icon("moon"),
                     class = "app-header-icon-btn", title = "Light mode only, for now"),
        tags$span(class = "app-header-avatar", title = "Account", icon("user"))
      )
    ),
    tags$script(HTML(
      "$(function(){
         $('#header_search').on('keydown', function(e){
           if (e.which === 13) {
             Shiny.setInputValue('header_search_submit', $(this).val(), {priority: 'event'});
           }
         });
         $(document).on('shown.bs.tab', 'a[data-toggle=\"tab\"]', function(e){
           var shownText = $(e.target).text().trim();
           $('.sidebar-nav-item').removeClass('active');
           $('.sidebar-nav-item[data-match=\"' + shownText.replace(/\"/g, '') + '\"]').addClass('active');
         });
       });"
    ))
  )
}

## ---------------------------------------------------------------------------
## Left sidebar: nav items are plain <a class="action-button"> elements (the
## same minimal pattern shiny::actionLink() itself renders down to), each
## with an id server.R can observeEvent() on, and a data-match attribute
## used only client-side (see the JS above) to highlight whichever item's
## destination tab is currently showing. `nav_items` is a list of
## list(id=, label=, icon=, match=) - `match` is the exact visible tab title
## the item navigates to. No ArthOChat card here - the header's persistent
## "Ask ArthOChat" button already covers every page; a second, contextual
## card lives inside the Preprocessing step's own content instead (see
## mod_preprocessing.R), rather than duplicating it here on every page.
## ---------------------------------------------------------------------------

omics_sidebar <- function(module_id, module_label, nav_items) {
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
    tags$div(class = "omics-sidebar-heading", "QUICK LINKS"),
    tags$div(
      class = "sidebar-quicklinks",
      tags$a(href = "#", icon("book"), "Documentation"),
      tags$a(href = "#", icon("graduation-cap"), "Tutorials"),
      tags$a(href = "#", icon("code"), "API Reference"),
      tags$a(href = "#", icon("clipboard-list"), "Release Notes")
    )
  )
}

## ---------------------------------------------------------------------------
## Right-column "Pipeline summary" vertical timeline. `steps` is a plain list
## of list(number=, label=, sublabel=, state=) with state one of
## "done"/"current"/"future" - computed by the caller from reactive signals
## that already exist (e.g. mod_preprocessing.R's own n_ready/merged_ok/
## batch_ok logic and the shared `results` reactiveValues), not by this
## function.
## ---------------------------------------------------------------------------

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
