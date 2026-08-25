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
  observeEvent(input$home_browse_modules, {
    updateTabsetPanel(session, "sidebar_tabs", selected = "modules")
  }, ignoreInit = TRUE)
  observeEvent(input$home_cta_browse_modules, {
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

  ## ---- Methylomics: shared dataset + computed-results store, separate from
  ## the transcriptomics `dataset`/`results` above - no default/preloaded
  ## methylation dataset exists, so every field starts NULL until the
  ## Methylomics Dataset tab (mod_methyl_dataset.R) loads something.
  methyl_dataset <- reactiveValues(
    beta = NULL, input_scale = NULL, array_type = NULL, sample_sheet = NULL,
    rg_set = NULL, mset = NULL, detp = NULL, beadcount = NULL, source = NULL,
    ## TRUE once "Load preloaded dataset" is clicked on the Dataset tab -
    ## every sub-module checks this first to decide whether to show the
    ## default GSE42861 analysis reproduction or the live upload-driven
    ## tool. Deliberately independent of `beta` (still NULL in preloaded
    ## mode - see METH_DATA_ROOT's comment in global.R for why the ~2.1GB
    ## QC'd matrix isn't bundled).
    preloaded = FALSE
  )
  methyl_results <- reactiveValues()
  mod_methyl_dataset_server("mx_dataset", methyl_dataset)
  lapply(MX_MODULES, function(m) m$server(paste0("mx_", m$config$id), methyl_dataset, methyl_results))

  ## ---- Cross-Omics: shared dataset + computed-results store, separate
  ## from `dataset`/`methyl_dataset` above - starts NULL until the
  ## Cross-Omics Dataset tab (mod_cross_dataset.R) loads one of the
  ## pipeline's own precomputed convergence/MR tables.
  cross_dataset <- reactiveValues(table_label = NULL, df = NULL, source = NULL)
  cross_results <- reactiveValues()
  mod_cross_dataset_server("cx_dataset", cross_dataset)
  ## The "integration" sub-module additionally receives the already-in-scope
  ## `dataset`/`results`/`methyl_dataset`/`methyl_results` handles (read-only)
  ## so its "Load from Transcriptomics/Methylomics" buttons can reuse results
  ## already generated in those tabs, without those modules being edited or
  ## their own tabs being affected. The other Cross-Omics sub-modules keep
  ## their original (id, cross_dataset, cross_results) call signature.
  lapply(CX_MODULES, function(m) {
    if (identical(m$config$id, "integration")) {
      m$server(paste0("cx_", m$config$id), cross_dataset, cross_results, dataset, results, methyl_dataset, methyl_results)
    } else {
      m$server(paste0("cx_", m$config$id), cross_dataset, cross_results)
    }
  })

  ## ---- Multi-Omics: shared dataset + computed-results store, separate from
  ## every other module's own reactiveValues. The Dataset Workspace tab
  ## (mod_multi_dataset.R) is the only place that writes to this - it
  ## publishes one validated Active Multi-Omics Dataset (source =
  ## "preloaded"/"upload"/"geo", `active` = TRUE once selected) that every
  ## other Multi-Omics sub-module reads from, via multi_active_dataset_banner()
  ## and, for Live Analysis, `layers`/`sample_meta` directly. `table_label`/
  ## `df` are kept for the Dataset Workspace's own precomputed-table browser.
  multi_dataset <- reactiveValues(
    table_label = NULL, df = NULL, source = NULL,
    layers = list(), layer_meta = list(), sample_meta = NULL,
    overlap = NULL, active = FALSE, loaded_at = NULL
  )
  multi_results <- reactiveValues()
  mod_multi_dataset_server("mo_dataset", multi_dataset, multi_results)
  lapply(MULTI_MODULES, function(m) m$server(paste0("mo_", m$config$id), multi_dataset, multi_results))

  ## ---- ArthOChat: one shared assistant session for the whole app, its own
  ## top-level tab (see ui.R) rather than nested inside a sub-module ---------
  mod_arthochat_server("arthochat", dataset, results, cross_results)

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

  ## ---- Methylomics Sub-modules tab: card Add/Remove toggles, count, and
  ## live filter - identical wiring to the Transcriptomics block above,
  ## against MX_MODULES/"mx_menu"/methyl_dataset+methyl_results instead of
  ## TX_MODULES/"tx_menu"/dataset+results, and the "mx_" id_prefix
  ## submoduleCardUI() was given (see ui.R) so DOM/input ids never collide
  ## with the Transcriptomics cards mounted in the same page.
  mx_added <- reactiveValues(ids = character(0))

  lapply(MX_MODULES, function(m) {
    hid <- m$config$id
    toggle_input <- paste0("mx_sm_toggle_", hid)

    observeEvent(input[[toggle_input]], {
      if (hid %in% mx_added$ids) {
        removeTab(session = session, inputId = "mx_menu", target = m$config$title)
        mx_added$ids <- setdiff(mx_added$ids, hid)
        shinyjs::removeClass(id = paste0("mx_smcard_", hid), class = "sm-card-active")
        shinyjs::html(id = paste0("mx_smstate_", hid), html = "Add")
      } else {
        insertTab(
          session = session, inputId = "mx_menu",
          tabPanel(m$config$title, br(), m$ui(paste0("mx_", hid))),
          target = "Sub-modules", position = "before", select = TRUE
        )
        mx_added$ids <- union(mx_added$ids, hid)
        shinyjs::addClass(id = paste0("mx_smcard_", hid), class = "sm-card-active")
        shinyjs::html(id = paste0("mx_smstate_", hid), html = "Added")
      }
    }, ignoreInit = TRUE)
  })

  output$mx_sm_active_count <- renderUI({
    n <- length(mx_added$ids)
    span(class = "sm-count-badge", sprintf("%d of %d sub-modules added", n, length(MX_MODULES)))
  })

  observeEvent(input$mx_sm_search, {
    q <- tolower(trimws(input$mx_sm_search %||% ""))
    lapply(MX_MODULES, function(m) {
      hid <- m$config$id
      match <- q == "" || grepl(q, tolower(m$config$title), fixed = TRUE)
      shinyjs::toggle(id = paste0("mx_smcard_wrap_", hid), condition = match)
    })
  }, ignoreNULL = FALSE)

  ## ---- Cross-Omics Sub-modules tab: card Add/Remove toggles, count, and
  ## live filter - identical wiring to the Transcriptomics/Methylomics
  ## blocks above, against CX_MODULES/"cx_menu"/cross_dataset+cross_results
  ## and the "cx_" id_prefix submoduleCardUI() was given (see ui.R).
  cx_added <- reactiveValues(ids = character(0))

  lapply(CX_MODULES, function(m) {
    hid <- m$config$id
    toggle_input <- paste0("cx_sm_toggle_", hid)

    observeEvent(input[[toggle_input]], {
      if (hid %in% cx_added$ids) {
        removeTab(session = session, inputId = "cx_menu", target = m$config$title)
        cx_added$ids <- setdiff(cx_added$ids, hid)
        shinyjs::removeClass(id = paste0("cx_smcard_", hid), class = "sm-card-active")
        shinyjs::html(id = paste0("cx_smstate_", hid), html = "Add")
      } else {
        insertTab(
          session = session, inputId = "cx_menu",
          tabPanel(m$config$title, br(), m$ui(paste0("cx_", hid))),
          target = "Sub-modules", position = "before", select = TRUE
        )
        cx_added$ids <- union(cx_added$ids, hid)
        shinyjs::addClass(id = paste0("cx_smcard_", hid), class = "sm-card-active")
        shinyjs::html(id = paste0("cx_smstate_", hid), html = "Added")
      }
    }, ignoreInit = TRUE)
  })

  output$cx_sm_active_count <- renderUI({
    n <- length(cx_added$ids)
    span(class = "sm-count-badge", sprintf("%d of %d sub-modules added", n, length(CX_MODULES)))
  })

  observeEvent(input$cx_sm_search, {
    q <- tolower(trimws(input$cx_sm_search %||% ""))
    lapply(CX_MODULES, function(m) {
      hid <- m$config$id
      match <- q == "" || grepl(q, tolower(m$config$title), fixed = TRUE)
      shinyjs::toggle(id = paste0("cx_smcard_wrap_", hid), condition = match)
    })
  }, ignoreNULL = FALSE)

  ## ---- Multi-Omics Sub-modules tab: card Add/Remove toggles, count, and
  ## live filter - identical wiring to the Transcriptomics/Methylomics/
  ## Cross-Omics blocks above, against MULTI_MODULES/"mo_menu"/
  ## multi_dataset+multi_results and the "mo_" id_prefix submoduleCardUI()
  ## was given (see ui.R).
  mo_added <- reactiveValues(ids = character(0))

  lapply(MULTI_MODULES, function(m) {
    hid <- m$config$id
    toggle_input <- paste0("mo_sm_toggle_", hid)

    observeEvent(input[[toggle_input]], {
      if (hid %in% mo_added$ids) {
        removeTab(session = session, inputId = "mo_menu", target = m$config$title)
        mo_added$ids <- setdiff(mo_added$ids, hid)
        shinyjs::removeClass(id = paste0("mo_smcard_", hid), class = "sm-card-active")
        shinyjs::html(id = paste0("mo_smstate_", hid), html = "Add")
      } else {
        insertTab(
          session = session, inputId = "mo_menu",
          tabPanel(m$config$title, br(), m$ui(paste0("mo_", hid))),
          target = "Sub-modules", position = "before", select = TRUE
        )
        mo_added$ids <- union(mo_added$ids, hid)
        shinyjs::addClass(id = paste0("mo_smcard_", hid), class = "sm-card-active")
        shinyjs::html(id = paste0("mo_smstate_", hid), html = "Added")
      }
    }, ignoreInit = TRUE)
  })

  output$mo_sm_active_count <- renderUI({
    n <- length(mo_added$ids)
    span(class = "sm-count-badge", sprintf("%d of %d sub-modules added", n, length(MULTI_MODULES)))
  })

  observeEvent(input$mo_sm_search, {
    q <- tolower(trimws(input$mo_sm_search %||% ""))
    lapply(MULTI_MODULES, function(m) {
      hid <- m$config$id
      match <- q == "" || grepl(q, tolower(m$config$title), fixed = TRUE)
      shinyjs::toggle(id = paste0("mo_smcard_wrap_", hid), condition = match)
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

  ## Methylomics sidebar nav - same jump-or-fall-back-to-picker pattern as
  ## jump_to_submodule() above, against MX_MODULES/"mx_menu"/mx_added$ids.
  jump_to_mx_submodule <- function(mod_id, sm_filter = NULL) {
    cfg <- MX_MODULES_BY_ID[[mod_id]]$config
    if (mod_id %in% mx_added$ids) {
      updateTabsetPanel(session, "mx_menu", selected = cfg$title)
    } else {
      updateTabsetPanel(session, "mx_menu", selected = "Sub-modules")
      updateTextInput(session, "mx_sm_search", value = sm_filter %||% cfg$title)
    }
  }

  observeEvent(input$sidebar_nav_methylomics_dataset, {
    updateTabsetPanel(session, "mx_menu", selected = "Dataset")
  }, ignoreInit = TRUE)
  observeEvent(input$sidebar_nav_methylomics_qc, {
    jump_to_mx_submodule("qc", sm_filter = "Quality Control")
  }, ignoreInit = TRUE)
  observeEvent(input$sidebar_nav_methylomics_normalization, {
    jump_to_mx_submodule("normalization", sm_filter = "Normalization")
  }, ignoreInit = TRUE)
  observeEvent(input$sidebar_nav_methylomics_dmp, {
    jump_to_mx_submodule("dmp", sm_filter = "Differential Methylation")
  }, ignoreInit = TRUE)
  ## Cell-Type Deconvolution has no dedicated sidebar entry of its own (it's
  ## still an unbuilt stub, only reachable via the Sub-modules grid) - this
  ## lets Quality Control's Cell Composition tab link out to it with the
  ## same jump-or-fall-back-to-picker mechanism as every sidebar link above,
  ## via a plain (non-namespaced) actionLink placed inside that tab's UI.
  observeEvent(input$sidebar_nav_methylomics_celltype, {
    jump_to_mx_submodule("celltype", sm_filter = "Cell-Type Deconvolution")
  }, ignoreInit = TRUE)
  observeEvent(input$sidebar_nav_methylomics_submodules, {
    updateTabsetPanel(session, "mx_menu", selected = "Sub-modules")
  }, ignoreInit = TRUE)

  ## Cross-Omics sidebar nav - same jump-or-fall-back-to-picker pattern as
  ## jump_to_mx_submodule() above, against CX_MODULES/"cx_menu"/cx_added$ids.
  jump_to_cx_submodule <- function(mod_id, sm_filter = NULL) {
    cfg <- CX_MODULES_BY_ID[[mod_id]]$config
    if (mod_id %in% cx_added$ids) {
      updateTabsetPanel(session, "cx_menu", selected = cfg$title)
    } else {
      updateTabsetPanel(session, "cx_menu", selected = "Sub-modules")
      updateTextInput(session, "cx_sm_search", value = sm_filter %||% cfg$title)
    }
  }

  observeEvent(input$sidebar_nav_crossomics_dataset, {
    updateTabsetPanel(session, "cx_menu", selected = "Dataset")
  }, ignoreInit = TRUE)
  observeEvent(input$sidebar_nav_crossomics_submodules, {
    updateTabsetPanel(session, "cx_menu", selected = "Sub-modules")
  }, ignoreInit = TRUE)

  ## Multi-Omics sidebar nav - same jump-or-fall-back-to-picker pattern as
  ## jump_to_cx_submodule() above, against MULTI_MODULES/"mo_menu"/mo_added$ids.
  jump_to_mo_submodule <- function(mod_id, sm_filter = NULL) {
    cfg <- MULTI_MODULES_BY_ID[[mod_id]]$config
    if (mod_id %in% mo_added$ids) {
      updateTabsetPanel(session, "mo_menu", selected = cfg$title)
    } else {
      updateTabsetPanel(session, "mo_menu", selected = "Sub-modules")
      updateTextInput(session, "mo_sm_search", value = sm_filter %||% cfg$title)
    }
  }

  observeEvent(input$sidebar_nav_multiomics_dataset, {
    updateTabsetPanel(session, "mo_menu", selected = "Dataset")
  }, ignoreInit = TRUE)
  observeEvent(input$sidebar_nav_multiomics_submodules, {
    updateTabsetPanel(session, "mo_menu", selected = "Sub-modules")
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

    mx_hit <- Find(function(m) grepl(q, tolower(m$config$title), fixed = TRUE), MX_MODULES)
    if (!is.null(mx_hit)) {
      updateTabsetPanel(session, "sidebar_tabs", selected = "methylomics")
      jump_to_mx_submodule(mx_hit$config$id, sm_filter = mx_hit$config$title)
      return()
    }

    cx_hit <- Find(function(m) grepl(q, tolower(m$config$title), fixed = TRUE), CX_MODULES)
    if (!is.null(cx_hit)) {
      updateTabsetPanel(session, "sidebar_tabs", selected = "crossomics")
      jump_to_cx_submodule(cx_hit$config$id, sm_filter = cx_hit$config$title)
      return()
    }

    mo_hit <- Find(function(m) grepl(q, tolower(m$config$title), fixed = TRUE), MULTI_MODULES)
    if (!is.null(mo_hit)) {
      updateTabsetPanel(session, "sidebar_tabs", selected = "multiomics")
      jump_to_mo_submodule(mo_hit$config$id, sm_filter = mo_hit$config$title)
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

  ## Methylomics page subtitle - mirrors tx_page_subtitle above, against mx_menu.
  output$mx_page_subtitle <- renderUI({
    sel <- input$mx_menu %||% "Dataset"
    txt <- switch(sel,
      "Dataset" = "Load a methylation dataset",
      "Sub-modules" = "Add or remove pipeline analyses",
      sel
    )
    p(txt)
  })

  ## Cross-Omics page subtitle - mirrors mx_page_subtitle above, against cx_menu.
  output$cx_page_subtitle <- renderUI({
    sel <- input$cx_menu %||% "Dataset"
    txt <- switch(sel,
      "Dataset" = "Browse the pipeline's biomarker-convergence tables",
      "Sub-modules" = "Add or remove pipeline analyses",
      sel
    )
    p(txt)
  })

  ## Multi-Omics page subtitle - mirrors cx_page_subtitle above, against mo_menu.
  ## No subtitle for the Dataset tab itself (its own header covers that).
  output$mo_page_subtitle <- renderUI({
    sel <- input$mo_menu %||% "Dataset"
    txt <- switch(sel,
      "Dataset" = NULL,
      "Sub-modules" = "Add or remove pipeline analyses",
      sel
    )
    if (is.null(txt)) return(NULL)
    p(txt)
  })
}
