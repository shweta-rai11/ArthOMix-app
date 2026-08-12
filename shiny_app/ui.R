## ui.R
## ArthOMix Explorer
## Layout follows xOmicsShiny's pattern: no marketing landing page, the app
## opens directly on Modules (pick an omics layer) -> Transcriptomics. Inside
## Transcriptomics: a Dataset tab (what's currently loaded), a Sub-modules
## tab (the grouped analysis picker), and one tab per analysis the user has
## added. Methylomics, Cross-Omics and Multi-Omics are placeholders for
## later work. ArthOChat is its own top-level tab, after Multi-Omics - one
## shared assistant for the whole app (sees whatever dataset/results are
## current, in every module), not a chat instance duplicated per sub-module.

## ---------------------------------------------------------------------------
## Home page
## First tab in the navbar, selected on load. Purely descriptive - every
## control on it hands off to machinery that already exists elsewhere
## (updateTabsetPanel(session, "sidebar_tabs", ...) in server.R, the same
## ARTHOCHAT_DRAWER_OPEN_JS every other "Ask ArthOChat" trigger uses) rather
## than introducing any new reactive logic of its own. Content mirrors the
## real pipeline: the six SUBMODULE_GROUP_ORDER stages defined below and the
## MODULE_REGISTRY omics layers from global.R, so it can't drift out of sync
## with what the app actually does. The hero's search box and quick-jump
## pills are genuinely functional, not decorative: both just set
## `header_search_submit`, the exact input R/ui_shell.R's own header search
## already sets and server.R's existing observeEvent(input$header_search_submit,
## ...) already handles - no new server-side logic needed. The stat-counter
## animation and FAQ accordion are pure client-side chrome (an
## IntersectionObserver count-up and native <details>/<summary>), scoped to
## this page's own script tag.
## ---------------------------------------------------------------------------

home_stat_card <- function(value, label, icn, count_to = NULL) {
  div(
    class = "home-stat-card",
    div(class = "home-stat-icon", icon(icn)),
    div(
      div(
        class = "home-stat-value",
        if (is.null(count_to)) value else tags$span(class = "home-stat-counter", `data-count-to` = count_to, "0")
      ),
      div(class = "home-stat-label", label)
    )
  )
}

home_pipeline_step <- function(n, title, blurb) {
  div(
    class = "home-pipeline-step",
    div(class = "home-pipeline-num", n),
    div(class = "home-pipeline-title", title),
    p(class = "home-pipeline-blurb", blurb)
  )
}

home_feature_card <- function(icn, title, desc) {
  div(
    class = "info-card home-feature-card",
    div(class = "info-icon", icon(icn)),
    h4(title),
    p(desc)
  )
}

## Quick-jump pill: sets the same `header_search_submit` input the hero
## search box and R/ui_shell.R's header search box set on Enter, so it
## reuses server.R's existing observeEvent(input$header_search_submit, ...)
## instead of adding a second navigation path.
home_hero_pill <- function(label, query = label) {
  tags$a(
    label, href = "#", class = "home-hero-pill",
    onclick = sprintf(
      "Shiny.setInputValue('header_search_submit', %s, {priority: 'event'}); return false;",
      jsonlite::toJSON(query, auto_unbox = TRUE)
    )
  )
}

home_faq_item <- function(question, answer) {
  tags$details(
    class = "home-faq-item",
    tags$summary(tags$span(question), tags$span(class = "home-faq-icon", icon("plus"))),
    div(class = "home-faq-answer", p(answer))
  )
}

## Decorative only (aria-hidden): five icons echoing the pipeline stages,
## floating over the hero's empty right-hand side in place of the static
## browser-chrome screenshot this used to be.
HOME_HERO_FLOAT_ICONS <- c("dna", "chart-column", "circle-nodes", "route", "stethoscope")

homeUI <- function() {
  tagList(
    div(
      class = "home-hero",
      div(class = "home-hero-pattern"),
      div(class = "home-hero-orb home-hero-orb-1", `aria-hidden` = "true"),
      div(class = "home-hero-orb home-hero-orb-2", `aria-hidden` = "true"),
      div(
        class = "home-hero-floating-icons", `aria-hidden` = "true",
        lapply(seq_along(HOME_HERO_FLOAT_ICONS), function(i) {
          div(class = paste0("home-hero-float-icon fi-", i), icon(HOME_HERO_FLOAT_ICONS[i]))
        })
      ),
      div(
        class = "home-hero-content",
        span(
          class = "home-hero-eyebrow",
          tags$span(class = "home-hero-live-dot"),
          "Sex-stratified multi-omics platform · Transcriptomics live"
        ),
        h1(class = "home-hero-title", "Find robust RA biomarkers, layer by layer."),
        p(
          class = "home-hero-subtitle",
          "ArthOMix Explorer takes a rheumatoid arthritis blood cohort from a raw expression matrix to a ",
          "validated diagnostic gene panel - sex-stratified differential expression, WGCNA co-expression, ",
          "Mendelian randomisation and colocalisation, and cross-tissue / cross-ancestry validation, all in one place."
        ),
        div(
          class = "home-hero-actions",
          actionButton("home_open_transcriptomics", tagList(icon("play"), "Open Transcriptomics"), class = "btn btn-primary"),
          actionButton("home_browse_modules", tagList(icon("table-cells-large"), "Browse all modules"), class = "btn btn-default"),
          tags$a(icon("comments"), " Ask ArthOChat", href = "#", class = "btn btn-default", onclick = ARTHOCHAT_DRAWER_OPEN_JS)
        ),
        div(class = "home-hero-search-label", "Or jump straight to a module, gene or dataset"),
        div(
          class = "home-hero-search",
          tags$span(class = "home-hero-search-icon", icon("magnifying-glass")),
          tags$input(
            type = "text", id = "home_hero_search", class = "form-control",
            placeholder = "Try “WGCNA”, “TNF”, “diagnostic model”…", autocomplete = "off"
          ),
          tags$span(class = "home-hero-search-kbd", "Enter")
        ),
        div(
          class = "home-hero-pills",
          home_hero_pill("Preprocessing"),
          home_hero_pill("Differential expression"),
          home_hero_pill("WGCNA"),
          home_hero_pill("Mendelian randomization"),
          home_hero_pill("Diagnostic model"),
          home_hero_pill("ArthOChat")
        ),
        div(
          class = "home-hero-trust",
          icon("database"),
          "Built on GSE93272 + GSE110169, harmonised against the Okada et al. 2014 RA GWAS - 1,701-gene MR instrument set"
        )
      ),
      tags$script(HTML(
        "$(function(){
           $('#home_hero_search').on('keydown', function(e){
             if (e.which === 13) {
               Shiny.setInputValue('header_search_submit', $(this).val(), {priority: 'event'});
             }
           });
         });"
      ))
    ),
    div(
      class = "home-stats-row",
      home_stat_card("15", "Analysis sub-modules", "layer-group", count_to = 15),
      home_stat_card("6", "Pipeline stages", "route", count_to = 6),
      home_stat_card("2", "Cohorts merged", "database", count_to = 2),
      home_stat_card("1 of 4", "Omics layers live", "circle-nodes")
    ),
    tags$script(HTML(
      "$(function(){
         var counters = document.querySelectorAll('.home-stat-counter');
         if (!counters.length) return;
         function animateCounter(el){
           var target = parseInt(el.getAttribute('data-count-to'), 10) || 0;
           var duration = 900, start = null;
           function step(ts){
             if (!start) start = ts;
             var progress = Math.min((ts - start) / duration, 1);
             var eased = 1 - Math.pow(1 - progress, 3);
             el.textContent = Math.round(eased * target);
             if (progress < 1) window.requestAnimationFrame(step);
           }
           window.requestAnimationFrame(step);
         }
         if ('IntersectionObserver' in window) {
           var obs = new IntersectionObserver(function(entries){
             entries.forEach(function(entry){
               if (entry.isIntersecting) { animateCounter(entry.target); obs.unobserve(entry.target); }
             });
           }, {threshold: 0.4});
           counters.forEach(function(c){ obs.observe(c); });
         } else {
           counters.forEach(animateCounter);
         }
       });"
    )),
    div(
      class = "page-header",
      h2("Why ArthOMix Explorer"),
      p("Built around how the underlying thesis pipeline actually works, not a generic omics dashboard.")
    ),
    div(
      class = "home-features-grid",
      home_feature_card("venus-mars", "Sex-stratified by design",
        "Differential expression, WGCNA, MR and the diagnostic model all run male/female contrasts side by side, not as an afterthought."),
      home_feature_card("route", "One pipeline, start to finish",
        "From a raw expression matrix to a validated gene panel - 15 sub-modules across 6 stages, without swivel-chairing between tools."),
      home_feature_card("comments", "ArthOChat, grounded in your data",
        "Ask about your dataset or results and get answers grounded in a live PubMed search with citations, not a generic chatbot."),
      home_feature_card("flask", "Validated, not just discovered",
        "Candidate panels are checked against synovium (cross-tissue) and independent ancestry cohorts before you trust them.")
    ),
    div(
      class = "page-header",
      h2("How the pipeline works"),
      p("Six stages, in the order the thesis pipeline actually runs. Open Transcriptomics and add any sub-module from its Sub-modules tab.")
    ),
    div(
      class = "home-pipeline-strip",
      lapply(seq_along(SUBMODULE_GROUP_ORDER), function(i) {
        g <- SUBMODULE_GROUP_ORDER[i]
        home_pipeline_step(i, g, SUBMODULE_GROUP_BLURB[[g]])
      })
    ),
    div(
      class = "page-header",
      h2("Omics modules"),
      p("Transcriptomics is live end to end. The rest are reserved for later layers of the same cohort.")
    ),
    div(
      class = "module-grid",
      ## "home_card_", not "home_" - the hero CTA button above is already
      ## actionButton("home_open_transcriptomics", ...); id_prefix = "home_"
      ## would produce that exact same id a second time, on this same page.
      lapply(Filter(function(m) m$kind %in% c("Single-omics", "Multi-omics"), MODULE_REGISTRY),
             function(m) moduleCardUI(m, id_prefix = "home_card_"))
    ),
    div(
      class = "page-header",
      h2("Frequently asked"),
      p("What people ask before loading their own cohort.")
    ),
    div(
      class = "home-faq-list",
      home_faq_item(
        "Do I need to write code to train the diagnostic model?",
        paste(
          "No. Every step through the diagnostic gene panel and its machine-learning model - feature",
          "selection, training, evaluation - is point-and-click inside the Biomarker modeling sub-modules.",
          "No R or Python required."
        )
      ),
      home_faq_item(
        "What data does ArthOMix Explorer ship with?",
        paste(
          "A merged rheumatoid arthritis blood cohort (GSE93272 + GSE110169), harmonised against the",
          "Okada et al. 2014 RA GWAS with a 1,701-gene Mendelian randomisation instrument set. It loads",
          "automatically, so you can explore the full pipeline before bringing your own data."
        )
      ),
      home_faq_item(
        "Can I load my own cohort instead?",
        "Yes - swap in your own expression matrix from the Dataset tab inside Transcriptomics; every downstream sub-module runs against whatever is currently loaded."
      ),
      home_faq_item(
        "Why run everything sex-stratified?",
        "Differential expression, WGCNA, Mendelian randomisation and the diagnostic model all run male/female contrasts side by side by design, not as an optional toggle - RA biomarkers can behave differently by sex, and averaging over that can hide a real signal."
      ),
      home_faq_item(
        "What about methylomics, proteomics or other omics layers?",
        "Transcriptomics is live end to end today. Methylomics, Cross-Omics and Multi-Omics are reserved on the Modules tab for later layers of the same cohort."
      )
    ),
    div(
      class = "home-cta-banner",
      div(class = "home-cta-pattern"),
      div(
        class = "home-cta-content",
        h3("Ready to work with your own cohort?"),
        p("Load a dataset on the Transcriptomics Dataset tab, or explore the example RA blood cohort that's already wired up."),
        actionButton("home_cta_open_transcriptomics", tagList(icon("arrow-right"), "Open Transcriptomics"), class = "btn btn-primary")
      )
    ),
    tags$footer(
      class = "home-footer",
      div(
        class = "home-footer-brand",
        tags$span(class = "app-header-brand-icon", icon("share-nodes")),
        div(
          div(class = "home-footer-title", "ArthOMix Explorer"),
          div(class = "home-footer-tagline", "Sex-stratified multi-omics biomarker discovery for rheumatoid arthritis.")
        )
      ),
      div(
        class = "home-footer-links",
        tags$a(href = "#", icon("book"), "Documentation"),
        tags$a(href = "#", icon("graduation-cap"), "Tutorials"),
        tags$a(href = "#", icon("code"), "API Reference"),
        tags$a(href = "#", icon("clipboard-list"), "Release Notes")
      )
    )
  )
}

## ---------------------------------------------------------------------------
## Modules landing page
## ---------------------------------------------------------------------------

## id_prefix distinguishes the "Open module" actionButton's inputId when
## this card is rendered on more than one tab at once - navbarPage mounts
## every tabPanel's content into the DOM upfront (not lazily on first
## visit), so the Modules tab's and Home tab's cards for the same module
## both exist simultaneously; without distinct ids they'd collide as a
## duplicate Shiny input id. Left at "" (unprefixed, same id as before) on
## the Modules tab itself so nothing there changes; homeUI() passes
## "home_" for its own copy, wired to the same updateTabsetPanel target in
## server.R.
moduleCardUI <- function(m, id_prefix = "") {
  is_available <- identical(m$status, "available")
  div(
    class = paste("module-card", if (!is_available) "module-card-disabled"),
    div(class = "module-card-icon", icon(m$icon)),
    div(
      class = "module-card-body",
      div(
        class = "module-card-title-row",
        h4(m$title),
        span(
          class = if (is_available) "badge-available" else "badge-soon",
          if (is_available) "Available" else "Coming soon"
        )
      ),
      p(class = "module-card-tagline", m$tagline),
      if (identical(m$id, "arthochat")) {
        ## ArthOChat is a slide-out drawer (see ui.R's app shell), not a
        ## navbarPage tab - open it the same way every other "Ask
        ## ArthOChat" trigger in the app does, not a server round-trip.
        tags$a("Open module", href = "#", class = "btn btn-primary btn-sm", onclick = ARTHOCHAT_DRAWER_OPEN_JS)
      } else if (is_available) {
        actionButton(paste0(id_prefix, "open_", m$id), "Open module", class = "btn-primary btn-sm")
      } else {
        tags$button("Open module", class = "btn btn-sm btn-default", disabled = "disabled")
      }
    )
  )
}

modulesLandingUI <- function() {
  single <- Filter(function(m) m$kind == "Single-omics", MODULE_REGISTRY)
  multi  <- Filter(function(m) m$kind == "Multi-omics", MODULE_REGISTRY)
  assist <- Filter(function(m) m$kind == "Assistant", MODULE_REGISTRY)
  tagList(
    div(
      class = "page-header",
      h2("Modules"),
      p("Each card is one omics layer. Transcriptomics is the only one running so far. Open it to load a dataset and add sub-modules; the rest are reserved for later.")
    ),
    h4(class = "module-group-title", "Single-omics modules"),
    div(class = "module-grid", lapply(single, moduleCardUI)),
    h4(class = "module-group-title", "Multi-omics modules"),
    div(class = "module-grid", lapply(multi, moduleCardUI)),
    if (length(assist) > 0) tagList(
      h4(class = "module-group-title", "Assistant"),
      div(class = "module-grid", lapply(assist, moduleCardUI))
    )
  )
}

## ---------------------------------------------------------------------------
## Transcriptomics module: Dataset tab, Sub-modules tab (grouped picker),
## dynamic analysis tabs inserted/removed by server.R
## ---------------------------------------------------------------------------

## Order controls the section order on the Sub-modules tab; it follows the
## order the thesis pipeline actually runs in.
SUBMODULE_GROUP_ORDER <- c(
  "Data", "Network", "Genetics", "Biomarker modeling", "Validation", "Interpretation"
)
SUBMODULE_GROUP_BLURB <- c(
  "Data" = "Load, inspect and prepare the expression matrix.",
  "Network" = "Co-expression structure and the candidate genes it points to.",
  "Genetics" = "Causal evidence from GWAS summary statistics.",
  "Biomarker modeling" = "Turn candidate genes into a panel and a diagnostic model.",
  "Validation" = "Check the panel holds up outside the discovery cohort.",
  "Interpretation" = "What the panel means biologically and clinically."
)

submoduleCardUI <- function(cfg) {
  div(
    id = paste0("smcard_wrap_", cfg$id), class = "sm-card-wrap",
    `data-sm-title` = tolower(cfg$title),
    div(
      class = "sm-card", id = paste0("smcard_", cfg$id),
      div(class = "sm-card-icon", icon(cfg$icon)),
      div(
        class = "sm-card-body",
        div(class = "sm-card-title-row", h4(cfg$title)),
        p(class = "sm-card-desc", cfg$description)
      ),
      tags$button(
        id = paste0("sm_toggle_", cfg$id), type = "button",
        class = "btn sm-toggle-btn action-button",
        span(id = paste0("smstate_", cfg$id), "Add")
      )
    )
  )
}

build_submodule_grid <- function() {
  by_group <- split(TX_MODULES, vapply(TX_MODULES, function(m) m$config$group %||% "Data", character(1)))
  groups <- intersect(SUBMODULE_GROUP_ORDER, names(by_group))
  tagList(
    lapply(groups, function(g) {
      tagList(
        div(
          class = "sm-group-header",
          h4(g),
          p(class = "sm-group-blurb", SUBMODULE_GROUP_BLURB[[g]])
        ),
        div(class = "sm-grid", lapply(by_group[[g]], function(m) submoduleCardUI(m$config)))
      )
    })
  )
}

## Left sidebar nav items for the Transcriptomics module (see
## R/ui_shell.R::omics_sidebar()). `match` is the exact visible tab title
## each item navigates to / highlights against - Overview and Datasets is
## one outer tx_menu tab; Preprocessing/Batch Correction/Merge Datasets are
## the three inner tabs of the "Preprocessing and Batch Correction"
## sub-module (see mod_preprocessing_ui's workflow stepper). The
## corresponding "sidebar_nav_transcriptomics_<id>" click observers live in
## server.R.
TRANSCRIPTOMICS_SIDEBAR_NAV <- list(
  list(id = "overview", label = "Overview", icon = "table-cells", match = "Overview and Datasets"),
  list(id = "dataset", label = "Datasets", icon = "database", match = "Overview and Datasets"),
  list(id = "preprocessing", label = "Preprocessing", icon = "broom", match = "Preprocessing"),
  list(id = "batchcorrection", label = "Batch Correction", icon = "wand-magic-sparkles", match = "Batch correction"),
  list(id = "mergedatasets", label = "Merge Datasets", icon = "code-merge", match = "Merge datasets"),
  list(id = "submodules", label = "Sub-modules", icon = "layer-group", match = "Sub-modules")
)

transcriptomicsUI <- function() {
  fluidRow(
    column(3, div(class = "omics-sidebar-col", omics_sidebar(
      "transcriptomics", "Transcriptomics", TRANSCRIPTOMICS_SIDEBAR_NAV
    ))),
    column(
      9,
      div(
        class = "page-header page-header-tight",
        div(class = "page-header-pattern"),
        h2(icon("dna"), " Transcriptomics"),
        uiOutput("tx_page_subtitle")
      ),
      div(
        class = "tx-menu-wrap",
        tabsetPanel(
          id = "tx_menu", type = "tabs",
          tabPanel("Dataset", br(), mod_dataset_ui("tx_dataset")),
          tabPanel(
            "Sub-modules", br(),
            div(
              class = "sm-toolbar",
              div(class = "sm-toolbar-count", uiOutput("sm_active_count")),
              div(
                class = "sm-toolbar-search",
                tags$span(icon("magnifying-glass")),
                textInput("sm_search", NULL, placeholder = "Filter sub-modules by name...", width = "260px")
              )
            ),
            build_submodule_grid()
          )
        )
      )
    )
  )
}

## ---------------------------------------------------------------------------
## Placeholder modules
## ---------------------------------------------------------------------------

comingSoonUI <- function(title, blurb) {
  tagList(
    div(class = "page-header", h2(title)),
    box(
      width = 12, status = "primary", solidHeader = FALSE,
      div(
        class = "coming-soon",
        icon("hammer", class = "coming-soon-icon"),
        h4("This module is not built yet"),
        p(blurb)
      )
    )
  )
}

## ---------------------------------------------------------------------------
## App shell
## Matches xOmicsShiny's actual page structure (github.com/interactivereport/
## xOmicsShiny, app.R): fluidPage + shinytheme("cerulean") + a plain
## titlePanel, with navbarPage nested inside for the tab bar - not
## shinydashboard's AdminLTE chrome. box()/valueBox() (shinydashboard
## components used throughout the sub-modules) still need AdminLTE's CSS to
## render correctly; shinydashboard normally attaches that only via
## dashboardPage(), so it's attached explicitly below. CSS only, not the
## AdminLTE/shinydashboard JS - that JS expects dashboardHeader/Sidebar DOM
## (sidebar toggle, box collapse) that doesn't exist here and throws on load;
## nothing in this app uses box(collapsible = TRUE) or any other bit of that
## JS, so it's safe to drop.
## ---------------------------------------------------------------------------

addCssDepsOnly <- function(tag) {
  deps <- htmltools::findDependencies(shinydashboard:::addDeps(tags$div()))
  js_free <- lapply(deps, function(d) { d$script <- character(0); d })
  htmltools::attachDependencies(tag, js_free, append = TRUE)
}

addCssDepsOnly(
  fluidPage(
    theme = shinythemes::shinytheme("cerulean"),
    useShinyjs(),
    tags$head(
      tags$title("ArthOMix Explorer"),
      ## ?v=<mtime> cache-busts these on every app restart - www/ assets have
      ## no Cache-Control header, so browsers are free to serve a stale copy
      ## from disk cache on a plain reload; a changing query string forces a
      ## real refetch without requiring a hard reload from whoever's testing.
      tags$link(rel = "stylesheet", type = "text/css",
                href = paste0("custom.css?v=", as.integer(file.mtime("www/custom.css")))),
      tags$link(rel = "stylesheet", type = "text/css",
                href = paste0("menuhex.css?v=", as.integer(file.mtime("www/menuhex.css"))))
    ),
    app_header(),
    navbarPage(
      title = "", id = "sidebar_tabs", selected = "home",
      tabPanel(tagList(icon("house"), "Home"), value = "home", homeUI()),
      tabPanel(tagList(icon("table-cells-large"), "Modules"), value = "modules", modulesLandingUI()),
      tabPanel(tagList(icon("dna"), "Transcriptomics"), value = "transcriptomics", transcriptomicsUI()),
      tabPanel(tagList(icon("circle-nodes"), "Methylomics (soon)"), value = "methylomics", comingSoonUI(
        "Methylomics",
        "DNA methylation analysis for this cohort will be added as a second single-omics module, following the same Dataset and Sub-modules layout as transcriptomics."
      )),
      tabPanel(tagList(icon("arrows-left-right"), "Cross-Omics (soon)"), value = "crossomics", comingSoonUI(
        "Cross-Omics",
        "Paired comparison between two omics layers, for example correlating methylation and expression at matched genes, will be added once methylomics is available."
      )),
      tabPanel(tagList(icon("layer-group"), "Multi-Omics (soon)"), value = "multiomics", comingSoonUI(
        "Multi-Omics",
        "Joint integration across all available omics layers into a single biomarker or disease model will be added last, once each single-omics module is complete."
      ))
      ## No "ArthOChat" tabPanel here - see the drawer below. Opening chat
      ## used to replace whatever module you were looking at with a
      ## separate full-page tab; it's now a slide-out panel mounted once,
      ## outside the navbarPage, so it overlays on top of the current page
      ## instead of navigating away from it.
    ),
    div(
      id = "arthochat_drawer", class = "chat-drawer",
      div(
        class = "chat-drawer-header",
        div(icon("comments"), strong(" ArthOChat")),
        tags$button(
          icon("xmark"), class = "chat-drawer-close", title = "Close",
          onclick = "document.getElementById('arthochat_drawer').classList.remove('open')"
        )
      ),
      div(class = "chat-drawer-body", mod_arthochat_ui("arthochat"))
    ),
    tags$div(
      id = "arthochat_drawer_backdrop", class = "chat-drawer-backdrop",
      onclick = "document.getElementById('arthochat_drawer').classList.remove('open')"
    )
  )
)
