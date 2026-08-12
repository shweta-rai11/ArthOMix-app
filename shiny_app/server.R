## server.R
## ArthOMix Explorer

function(input, output, session) {

  ## ---- Shared dataset, read by every Transcriptomics submodule -----------
  dataset <- local({
    d <- load_default_dataset()
    do.call(reactiveValues, d)
  })

  ## ---- Shared computed-results store, read by the Assistant sub-module ---
  ## Populated by each analysis module right after its main reactive fires;
  ## NULL until that module has been run at least once this session.
  results <- reactiveValues()

  ## Reset whenever the loaded dataset changes, so stale results from a
  ## previous dataset aren't shown to the assistant as if still current.
  observeEvent(dataset$source, {
    for (nm in names(results)) results[[nm]] <- NULL
  }, ignoreInit = TRUE)

  ## Always-visible header badge (R/ui_shell.R::app_header()) - which
  ## dataset every Transcriptomics submodule is currently reading from.
  ## Neutral grey for this app's own preloaded example; blue, plus sample
  ## count, for anyone's own uploaded/merged/corrected data - so switching
  ## between the two is never ambiguous no matter which page you're on.
  output$active_dataset_badge <- renderUI({
    src <- dataset$source %||% ""
    is_default <- grepl("^Example dataset:", src)
    n <- tryCatch(ncol(dataset$expr), error = function(e) NA)
    tags$span(
      class = paste("header-dataset-badge", if (is_default) "is-default" else "is-custom"),
      icon(if (is_default) "flask" else "table"),
      if (is_default) {
        sprintf("ArthOMix default%s", if (!is.na(n)) sprintf(" (%s samples)", n) else "")
      } else {
        sprintf("Your own data%s", if (!is.na(n)) sprintf(" (%s samples)", n) else "")
      }
    )
  })

  ## ---- Home page -> jump into the app -------------------------------------
  observeEvent(input$home_open_transcriptomics, {
    updateTabsetPanel(session, "sidebar_tabs", selected = "transcriptomics")
  }, ignoreInit = TRUE)
  observeEvent(input$home_cta_open_transcriptomics, {
    updateTabsetPanel(session, "sidebar_tabs", selected = "transcriptomics")
  }, ignoreInit = TRUE)
  observeEvent(input$home_browse_modules, {
    updateTabsetPanel(session, "sidebar_tabs", selected = "modules")
  }, ignoreInit = TRUE)
  ## Home page's own "Omics modules" grid (homeUI() in ui.R, moduleCardUI()
  ## called there with id_prefix = "home_card_" so its buttons don't
  ## collide with the Modules tab's identical cards, mounted in the DOM at
  ## the same time - see moduleCardUI()'s own comment in ui.R).
  lapply(MODULE_REGISTRY, function(m) {
    if (identical(m$status, "available") && !identical(m$id, "arthochat")) {
      observeEvent(input[[paste0("home_card_open_", m$id)]], {
        updateTabsetPanel(session, "sidebar_tabs", selected = m$tab)
      }, ignoreInit = TRUE)
    }
  })

  ## ---- Modules landing page -> open a specific module ---------------------
  ## "arthochat" is excluded: it's a slide-out drawer now (see ui.R), opened
  ## by a plain client-side onclick (ARTHOCHAT_DRAWER_OPEN_JS) rather than
  ## an actionButton/updateTabsetPanel round-trip - see moduleCardUI() in
  ## ui.R for its card's "Open module" link.
  lapply(MODULE_REGISTRY, function(m) {
    if (identical(m$status, "available") && !identical(m$id, "arthochat")) {
      observeEvent(input[[paste0("open_", m$id)]], {
        updateTabsetPanel(session, "sidebar_tabs", selected = m$tab)
      }, ignoreInit = TRUE)
    }
  })

  ## ---- Dataset tab ---------------------------------------------------------
  mod_dataset_server("tx_dataset", dataset)

  ## ---- ArthOChat: one shared assistant session for the whole app, its own
  ## top-level tab (see ui.R) rather than nested inside a sub-module ---------
  mod_arthochat_server("arthochat", dataset, results)

  ## ---- Instantiate every submodule server once; visibility is controlled --
  ## by insertTab()/removeTab() below, not by when the server is created.
  lapply(TX_MODULES, function(m) m$server(paste0("tx_", m$config$id), dataset, results))

  ## ---- Sub-modules tab: card Add/Remove toggles ----------------------------
  added <- reactiveValues(ids = character(0))

  lapply(TX_MODULES, function(m) {
    hid <- m$config$id
    toggle_input <- paste0("sm_toggle_", hid)

    observeEvent(input[[toggle_input]], {
      if (hid %in% added$ids) {
        removeTab(session = session, inputId = "tx_menu", target = m$config$title)
        added$ids <- setdiff(added$ids, hid)
        shinyjs::removeClass(id = paste0("smcard_", hid), class = "sm-card-active")
        shinyjs::html(id = paste0("smstate_", hid), html = "Add")
      } else {
        insertTab(
          session = session, inputId = "tx_menu",
          tabPanel(m$config$title, br(), m$ui(paste0("tx_", hid))),
          target = "Sub-modules", position = "before", select = TRUE
        )
        added$ids <- union(added$ids, hid)
        shinyjs::addClass(id = paste0("smcard_", hid), class = "sm-card-active")
        shinyjs::html(id = paste0("smstate_", hid), html = "Added")
      }
    }, ignoreInit = TRUE)
  })

  output$sm_active_count <- renderUI({
    n <- length(added$ids)
    span(
      class = "sm-count-badge",
      sprintf("%d of %d sub-modules added", n, length(TX_MODULES))
    )
  })

  ## ---- Sub-modules tab: live filter by name --------------------------------
  observeEvent(input$sm_search, {
    q <- tolower(trimws(input$sm_search %||% ""))
    lapply(TX_MODULES, function(m) {
      hid <- m$config$id
      match <- q == "" || grepl(q, tolower(m$config$title), fixed = TRUE)
      shinyjs::toggle(id = paste0("smcard_wrap_", hid), condition = match)
    })
  }, ignoreNULL = FALSE)

  ## ===========================================================================
  ## App-shell navigation (header search, left sidebar, theme toggle, page
  ## subtitle) - additive UI glue only. Every call below is either a plain
  ## updateTabsetPanel()/updateTextInput() (the same functions the "Modules"
  ## landing page's open_ observers and the Sub-modules search box already
  ## use) or a read of state that already exists (`added$ids`, `input$tx_menu`,
  ## `TX_MODULES`/`TX_MODULES_BY_ID`, `MODULE_REGISTRY`). No reactive
  ## computation, analysis logic, or existing observer is touched.
  ## ===========================================================================

  ## Jump into one Transcriptomics sub-module tab, and (if it's been added)
  ## optionally straight to one of its own inner tabs - e.g. "preprocessing"'s
  ## "Merge datasets"/"Batch correction" steps. If the sub-module hasn't been
  ## added yet, falls back to the Sub-modules picker with the search box
  ## pre-filled, reusing the existing sm_search filter observer above instead
  ## of duplicating its logic.
  jump_to_submodule <- function(mod_id, inner_tab = NULL, sm_filter = NULL) {
    cfg <- TX_MODULES_BY_ID[[mod_id]]$config
    if (mod_id %in% added$ids) {
      updateTabsetPanel(session, "tx_menu", selected = cfg$title)
      if (!is.null(inner_tab)) {
        updateTabsetPanel(session, paste0("tx_", mod_id, "-tabs"), selected = inner_tab)
      }
    } else {
      updateTabsetPanel(session, "tx_menu", selected = "Sub-modules")
      updateTextInput(session, "sm_search", value = sm_filter %||% cfg$title)
    }
  }

  observeEvent(input$sidebar_nav_transcriptomics_overview, {
    jump_to_submodule("overview", sm_filter = "Overview")
  }, ignoreInit = TRUE)
  observeEvent(input$sidebar_nav_transcriptomics_dataset, {
    jump_to_submodule("overview", inner_tab = "Datasets", sm_filter = "Overview")
  }, ignoreInit = TRUE)
  observeEvent(input$sidebar_nav_transcriptomics_preprocessing, {
    jump_to_submodule("preprocessing", inner_tab = "Preprocessing", sm_filter = "Preprocessing")
  }, ignoreInit = TRUE)
  observeEvent(input$sidebar_nav_transcriptomics_batchcorrection, {
    jump_to_submodule("preprocessing", inner_tab = "Batch correction", sm_filter = "Preprocessing")
  }, ignoreInit = TRUE)
  observeEvent(input$sidebar_nav_transcriptomics_mergedatasets, {
    jump_to_submodule("preprocessing", inner_tab = "Merge datasets", sm_filter = "Preprocessing")
  }, ignoreInit = TRUE)
  observeEvent(input$sidebar_nav_transcriptomics_submodules, {
    updateTabsetPanel(session, "tx_menu", selected = "Sub-modules")
  }, ignoreInit = TRUE)

  ## Header search: "Enter" in the search box (see R/ui_shell.R::app_header())
  ## sets this input; matched first against top-level modules, then against
  ## Transcriptomics sub-module titles, via the same jump helpers above.
  observeEvent(input$header_search_submit, {
    q <- tolower(trimws(input$header_search_submit %||% ""))
    req(nzchar(q))

    mod_hit <- Find(function(m) identical(m$status, "available") && grepl(q, tolower(m$title), fixed = TRUE), MODULE_REGISTRY)
    if (!is.null(mod_hit)) {
      if (identical(mod_hit$id, "arthochat")) {
        shinyjs::runjs(ARTHOCHAT_DRAWER_OPEN_JS_STATEMENT)
      } else {
        updateTabsetPanel(session, "sidebar_tabs", selected = mod_hit$tab)
      }
      return()
    }

    sm_hit <- Find(function(m) grepl(q, tolower(m$config$title), fixed = TRUE), TX_MODULES)
    if (!is.null(sm_hit)) {
      updateTabsetPanel(session, "sidebar_tabs", selected = "transcriptomics")
      jump_to_submodule(sm_hit$config$id, sm_filter = sm_hit$config$title)
      return()
    }

    showNotification(sprintf('No module or sub-module matched "%s".', input$header_search_submit), type = "warning")
  }, ignoreInit = TRUE)

  ## Theme toggle: light-mode-only for now (see R/ui_shell.R::app_header()).
  observeEvent(input$theme_toggle_btn, {
    showNotification("Dark mode is coming soon.", type = "message", duration = 3)
  }, ignoreInit = TRUE)

  ## Transcriptomics page subtitle, directly under the "Transcriptomics"
  ## title - mirrors whichever tx_menu tab (Dataset / Overview and Datasets /
  ## Preprocessing and Batch Correction / Sub-modules / ...) is currently
  ## selected. Purely presentational; reads the existing tx_menu input.
  output$tx_page_subtitle <- renderUI({
    sel <- input$tx_menu %||% "Dataset"
    txt <- switch(sel,
      "Dataset" = "Load and manage the working dataset",
      "Sub-modules" = "Add or remove pipeline analyses",
      sel
    )
    p(txt)
  })
}
