## R/transcriptomics/01_Data/mod_dataset.R
## Dataset tab: pick a preloaded dataset, upload your own, or fetch from GEO.

mod_dataset_config <- list(
  id = "dataset", 
  title = "Dataset", 
  icon = "database",
  description = "Pick a preloaded dataset or upload your own - either way it's what every sub-module below reads from."
)

INDIVIDUAL_DATASET_LABELS <- c(
  "GSE93272"  = "Whole Blood Training Cohort A",
  "GSE110169" = "Whole Blood Training Cohort B",
  "GSE15573"  = "PBMC Validation Cohort",
  "GSE89408"  = "Synovial Tissue Validation Cohort"
)

individual_dataset_entry <- function(gse_id) {
  list(
    id = gse_id,
    label = if (gse_id %in% names(INDIVIDUAL_DATASET_LABELS)) INDIVIDUAL_DATASET_LABELS[[gse_id]] else gse_id,
    load = function() {
      d <- load_individual_dataset(gse_id)
      validate(need(!is.null(d), paste("Raw data for", gse_id, "was not found on disk.")))
      list(expr = d$expr, meta = d$meta, source = paste0("Individual dataset: ", d$label))
    }
  )
}

default_dataset_entry <- list(
  id = "__default_merged__",
  label = "Merged Data",
  load = function() {
    d <- load_default_dataset()
    list(expr = d$expr, meta = d$meta, source = d$source)
  }
)

PRELOADED_DATASETS <- c(
  list(default_dataset_entry),
  lapply(vapply(GEO_SOURCES, `[[`, character(1), "gse"), individual_dataset_entry)
)

preloaded_choices <- function() {
  setNames(vapply(PRELOADED_DATASETS, `[[`, character(1), "id"),
           vapply(PRELOADED_DATASETS, `[[`, character(1), "label"))
}

MERGED_DEFAULT_GEO_IDS <- c("GSE93272", "GSE110169")

entry_geo_ids <- function(entry) {
  if (identical(entry$id, default_dataset_entry$id)) MERGED_DEFAULT_GEO_IDS else entry$id
}

mod_dataset_ui <- function(id) {
  ns <- NS(id) # NS means namespace
  tagList(  # returns multiple UI components together
    fluidRow(  # for bootstrap grid layout
      column(
        6,
        box( # Box is a UI container in shiny dashboard
          width = NULL, title = "Preloaded Datasets", status = "primary", solidHeader = FALSE,
          selectInput(ns("preloaded_choice"), "Individual dataset",
                      choices = preloaded_choices(), selected = character(0), width = "100%"),
          uiOutput(ns("preloaded_note")),
          uiOutput(ns("preloaded_geo_card_ui")),
          div(style = "display: flex; align-items: center; gap: 10px; flex-wrap: wrap;",
              actionButton(ns("load_preloaded_btn"), "Load this dataset", icon = icon("rotate-left"), class = "btn-primary btn-sm"),
              uiOutput(ns("preloaded_load_message"), inline = TRUE))
        ),
        box(
          width = NULL, title = "Fetch from NCBI GEO", status = "primary", solidHeader = FALSE,
          p(strong("GEO Series accession"), " (e.g. ", code("GSE12345"), "). Fetches the expression data and metadata from NCBI. For your own data, use \"Upload your own data\" instead."),
          fluidRow(
            column(4, textInput(ns("geo_accession"), NULL, placeholder = "GSE12345", width = "100%")),
            column(3, actionButton(ns("geo_fetch_btn"), "Fetch", icon = icon("download"), class = "btn-primary btn-sm", width = "100%"))
          ),

          tags$script(HTML(sprintf(
            "document.addEventListener('click', function(e) {
               var btn = e.target.closest('#%s');
               if (!btn) return;
               var statusDiv = document.getElementById('%s');
               if (!statusDiv) return;
               btn.disabled = true;
               statusDiv.innerHTML = '<div class=\"empty-note\"><i class=\"fas fa-spinner fa-spin\" role=\"presentation\"></i> Fetching from NCBI GEO - this can take a few seconds to a couple of minutes. Please don\\'t navigate away or click Fetch again.</div>';
               var observer = new MutationObserver(function() {
                 var liveBtn = document.getElementById('%s');
                 if (liveBtn) liveBtn.disabled = false;
                 observer.disconnect();
               });
               observer.observe(statusDiv, { childList: true, subtree: true });
             });",
            ns("geo_fetch_btn"), ns("geo_fetch_status"), ns("geo_fetch_btn")
          ))),
          uiOutput(ns("geo_fetch_status")),
          uiOutput(ns("geo_platform_ui")),
          uiOutput(ns("geo_column_mapping")),
          actionButton(ns("geo_load_btn"), "Load this dataset", icon = icon("upload"), class = "btn-primary btn-sm")
        )
      ),
      column(
        6,
        box(
          width = NULL, title = "Upload your own data", status = "primary", solidHeader = FALSE,
          div(class = "upload-step-label", "STEP 1 · Choose your files"),
          p(strong("Expression matrix"), " - CSV or RDS. Genes in rows, samples in columns; for CSV, the first column is the gene ID."),
          radioButtons(ns("declared_data_type"), "Data type", inline = TRUE,
                       choices = c("Raw counts" = "raw", "Normalized (TPM/FPKM/CPM)" = "normalized", "Already log-transformed" = "logtransformed"),
                       selected = "normalized"),
          fileInput(ns("expr_file"), "Expression matrix", accept = c(".csv", ".rds", ".Rds")),
          p(strong("Sample metadata"), " - CSV or RDS data frame, one row per sample."),
          fileInput(ns("meta_file"), "Sample metadata", accept = c(".csv", ".rds", ".Rds")),
          uiOutput(ns("upload_preview_ui")),
          uiOutput(ns("upload_preview_tables_ui")),
          tags$hr(),
          div(class = "upload-step-label", "STEP 2 · Map the columns"),
          uiOutput(ns("column_mapping")),
          tags$hr(),
          div(class = "upload-step-label", "STEP 3 · Confirm"),
          actionButton(ns("load_btn"), "Upload Data", icon = icon("upload"), class = "btn-primary btn-sm"),
          div(class = "empty-note", style = "margin-top:6px;", icon("circle-info"),
              "Loads exactly what you provide, as-is - no merging, normalising, or batch correction.")
        )
      )
    ),
    uiOutput(ns("load_message")),
    box(
      width = NULL, title = "Reserve validation samples (sealed hold-out)", status = "primary", solidHeader = FALSE,
      p(class = "submodule-desc",
        "Set aside samples before any analysis runs. Reserved samples are hidden from every step (DE, WGCNA, Candidate Genes, Feature Selection, Model Training). They're only scored later in Diagnostic Model → External Validation, giving the one result that's leakage-safe end-to-end. Note: this turns off the fast precomputed results, since those were built using all samples."),
      fluidRow(
        column(3, numericInput(ns("reserve_frac_pct"), "Reserved fraction (%)", value = 30, min = 10, max = 50, step = 5)),
        column(3, numericInput(ns("reserve_seed"), "Split seed", value = 1234, min = 1, step = 1)),
        column(3, checkboxInput(ns("reserve_stratify_sex"), "Stratify by sex as well as group", value = TRUE)),
        column(3, div(style = "margin-top: 26px; display: flex; gap: 8px; flex-wrap: wrap;",
                      actionButton(ns("reserve_btn"), "Reserve", icon = icon("lock"), class = "btn-primary btn-sm"),
                      actionButton(ns("release_btn"), "Release", icon = icon("lock-open"), class = "btn-default btn-sm")))
      ),
      uiOutput(ns("reserve_status"))
    )
  )
}

DATASET_RESERVE_MIN_TOTAL <- 20L
DATASET_RESERVE_MIN_PER_GROUP <- 4L

## Stratified split into discovery/reserved samples. Pure function, testable outside Shiny.
tx_reserve_validation_ids <- function(meta, frac, seed = 1234, stratify_sex = TRUE) {
  stopifnot(is.data.frame(meta), "sample" %in% colnames(meta), "group" %in% colnames(meta))
  frac <- min(max(as.numeric(frac), 0.1), 0.5)
  ok <- !is.na(meta$group)
  meta <- meta[ok, , drop = FALSE]
  if (nrow(meta) < DATASET_RESERVE_MIN_TOTAL) stop(sprintf("At least %d samples with a group label are needed to reserve a validation set (this dataset has %d).", DATASET_RESERVE_MIN_TOTAL, nrow(meta)))
  grp <- as.character(meta$group)
  if (length(unique(grp)) < 2) stop("The group column needs at least two distinct values to reserve a stratified validation set.")
  strata <- grp
  if (isTRUE(stratify_sex) && "sex" %in% colnames(meta) && any(!is.na(meta$sex))) {
    sx <- as.character(meta$sex); sx[is.na(sx)] <- "unknown"
    strata <- paste(grp, sx, sep = "|")
  }
  set.seed(as.integer(seed))
  strata_f <- factor(strata)
  train_idx <- as.integer(caret::createDataPartition(strata_f, p = 1 - frac, list = FALSE)[, 1])
  reserved <- meta$sample[setdiff(seq_len(nrow(meta)), train_idx)]
  disc_tab <- table(grp[train_idx]); res_tab <- table(grp[-train_idx])
  if (any(disc_tab < DATASET_RESERVE_MIN_PER_GROUP) || length(res_tab) < 2 || any(res_tab < DATASET_RESERVE_MIN_PER_GROUP)) {
    stop(sprintf("Reserving %d%% would leave fewer than %d samples in a group on one side of the split (discovery: %s; reserved: %s). Lower the fraction or use a larger dataset.",
                 round(frac * 100), DATASET_RESERVE_MIN_PER_GROUP,
                 paste(sprintf("%s = %d", names(disc_tab), as.integer(disc_tab)), collapse = ", "),
                 paste(sprintf("%s = %d", names(res_tab), as.integer(res_tab)), collapse = ", ")))
  }
  reserved
}

DATASET_RESERVED_SUFFIX_RE <- " \\[[0-9]+ validation samples reserved\\]$"

mod_dataset_server <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$preloaded_note <- renderUI({
      req(input$preloaded_choice)
      if (identical(input$preloaded_choice, default_dataset_entry$id)) {
        return(p(class = "empty-note", icon("circle-info"),
          "The same merged, batch-corrected training cohort the app loads by default on startup - pick this to switch back to it after loading something else."))
      }
      raw_unavailable <- input$preloaded_choice %in% c("GSE93272", "GSE110169") &&
        {
          eset <- get_raw_eset(input$preloaded_choice)
          is.null(eset) || nrow(Biobase::exprs(eset)) == 0
        }
      if (raw_unavailable) {
        p(class = "empty-note", icon("circle-info"),
          "This source's raw probe-level file isn't available here. This loads only its samples, filtered from the merged, batch-corrected training cohort - not raw single-platform data. Pick \"Merged Data\" above to see it merged with the other source.")
      } else {
        p(class = "empty-note", icon("triangle-exclamation"),
          "Raw, single-platform data - probe-level, not merged or normalised. Sub-modules can run on it directly, but most assume the merged cohort, so results may differ. To inspect it without changing what sub-modules run on, use the QC tab on Overview and Datasets instead.")
      }
    })

    output$preloaded_geo_card_ui <- renderUI({
      req(input$preloaded_choice)
      if (identical(input$preloaded_choice, default_dataset_entry$id)) return(NULL)
      gse_id <- input$preloaded_choice
      src <- Find(function(s) identical(s$gse, gse_id), GEO_SOURCES)
      req(src)
      eset <- get_raw_eset(gse_id)
      div(
        class = "info-card",
        div(
          class = "module-card-title-row",
          h4(gse_id),
          tags$a(href = geo_link(gse_id), target = "_blank", rel = "noopener",
                  icon("up-right-from-square"), " NCBI GEO")
        ),
        if (!is.null(eset)) {
          tagList(
            p(class = "module-card-tagline",
              tryCatch(Biobase::experimentData(eset)@title, error = arthomix_null_on_error)),
            p(strong("Role: "), src$role, br(), strong("Used for: "), src$used_in),
            p(strong("Platform: "), Biobase::annotation(eset), br(),
              strong("Samples: "), ncol(eset), ", ", strong("Probes: "), format(nrow(eset), big.mark = ","))
          )
        } else {
          tagList(
            p(strong("Role: "), src$role, br(), strong("Used for: "), src$used_in),
            div(class = "empty-note", icon("triangle-exclamation"), "Raw file not found on disk.")
          )
        }
      )
    })

    meta_raw <- reactive({
      req(input$meta_file)
      path <- input$meta_file$datapath
      if (grepl("\\.rds$", input$meta_file$name, ignore.case = TRUE)) {
        loaded <- safe_read_rds(path)
        validate(need(isTRUE(loaded$ok), loaded$error %||% "Could not read this .rds file."))
        d <- loaded$value
        validate(need(is.data.frame(d), "The uploaded metadata RDS file must contain a data frame."))
        as.data.frame(d)
      } else {
        as.data.frame(data.table::fread(path, showProgress = FALSE))
      }
    })

    expr_raw <- reactive({
      req(input$expr_file)
      path <- input$expr_file$datapath
      if (grepl("\\.rds$", input$expr_file$name, ignore.case = TRUE)) {
        loaded <- safe_read_rds(path)
        validate(need(isTRUE(loaded$ok), loaded$error %||% "Could not read this .rds file."))
        loaded$value
      } else {
        m <- as.data.frame(data.table::fread(path, showProgress = FALSE))
        rn <- as.character(m[[1]])
        m <- as.matrix(m[, -1, drop = FALSE])
        if (!is.numeric(m)) {
          storage.mode(m) <- "character"
          m_num <- arthomix_quiet(matrix(as.numeric(m), nrow = nrow(m), ncol = ncol(m), dimnames = dimnames(m)))
          validate(need(!any(!is.na(m) & m != "" & is.na(m_num)),
            "This file has non-numeric values outside the first (gene ID) column. Everything after the ID column must be numeric - check for stray text, footnotes, or thousands separators."))
          m <- m_num
        }
        rownames(m) <- rn
        m
      }
    })

    upload_preview_data <- reactive({
      req(input$expr_file, input$meta_file)
      tryCatch(list(expr = expr_raw(), meta = meta_raw()), error = function(e) e)
    })

    output$upload_preview_ui <- renderUI({
      preview <- upload_preview_data()
      if (inherits(preview, "error")) {
        return(div(class = "empty-note", icon("triangle-exclamation"),
                    paste("Could not read the uploaded file(s):", conditionMessage(preview))))
      }
      div(class = "empty-note", icon("circle-info"),
          sprintf("Read %s: %s features x %s samples. Read %s: %s rows. Map the columns below, then click \"Load dataset\".",
                  input$expr_file$name, format(nrow(preview$expr), big.mark = ","), ncol(preview$expr),
                  input$meta_file$name, nrow(preview$meta)))
    })

    output$upload_preview_tables_ui <- renderUI({
      preview <- upload_preview_data()
      req(!inherits(preview, "error"))
      tagList(
        div(class = "upload-step-label", "Preview of what you uploaded"),
        p(class = "submodule-desc", "First 5 rows/columns, exactly as read - before any column mapping below."),
        strong("Sample metadata (first 5 rows)"),
        DT::dataTableOutput(ns("upload_preview_meta_table")),
        br(),
        strong("Expression matrix (first 5 features x 5 samples)"),
        DT::dataTableOutput(ns("upload_preview_expr_table"))
      )
    })

    output$upload_preview_meta_table <- DT::renderDataTable({
      preview <- upload_preview_data()
      req(!inherits(preview, "error"))
      DT::datatable(head(preview$meta, 5), rownames = FALSE,
                     options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    output$upload_preview_expr_table <- DT::renderDataTable({
      preview <- upload_preview_data()
      req(!inherits(preview, "error"))
      m <- preview$expr
      cols <- seq_len(min(5, ncol(m)))
      df <- data.frame(feature = rownames(m), round(m[, cols, drop = FALSE], 3), check.names = FALSE)
      DT::datatable(head(df, 5), rownames = FALSE,
                     options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    guess_col <- function(cols, exact, contains = exact, fallback = cols[1]) {
      for (term in exact) {
        hit <- cols[tolower(cols) == tolower(term)]
        if (length(hit) > 0) return(structure(hit[1], matched = TRUE))
      }
      for (term in contains) {
        hit <- cols[grepl(term, cols, ignore.case = TRUE)]
        if (length(hit) > 0) return(structure(hit[1], matched = TRUE))
      }
      structure(fallback, matched = FALSE)
    }

    unconfident_guess_note <- function(guess, label) {
      if (isTRUE(attr(guess, "matched"))) return(NULL)
      div(class = "empty-note", icon("triangle-exclamation"),
          sprintf("Couldn't confidently guess the %s from its column name - double-check the selection above.", label))
    }

    output$column_mapping <- renderUI({
      req(input$meta_file)
      cols <- colnames(meta_raw())
      id_guess <- guess_col(cols, c("sample", "sample_id", "id", "geo_accession", "accession"))
      group_guess <- guess_col(cols, c("group", "diagnosis", "disease", "condition", "status", "phenotype"))
      tagList(
        selectInput(ns("map_id"), "Sample ID column", choices = cols, selected = id_guess, selectize = FALSE),
        unconfident_guess_note(id_guess, "sample-ID column"),
        selectInput(ns("map_group"), "Group / diagnosis column", choices = cols, selected = group_guess, selectize = FALSE),
        unconfident_guess_note(group_guess, "group/diagnosis column"),
        selectInput(ns("map_sex"), "Sex column (optional)", choices = c("(none)", cols),
                    selected = guess_col(cols, c("sex", "gender"), fallback = "(none)"),
                    selectize = FALSE),
        selectInput(ns("map_batch"), "Batch column (optional)", choices = c("(none)", cols),
                    selected = guess_col(cols, c("batch", "cohort", "platform", "dataset"), fallback = "(none)"),
                    selectize = FALSE)
      )
    })

    observe({
      ready <- !is.null(input$expr_file) && !is.null(input$meta_file) &&
        !is.null(input$map_id) && !is.null(input$map_group)
      if (isTRUE(ready)) shinyjs::enable("load_btn") else shinyjs::disable("load_btn")
    })

    observe({
      shinyjs::toggleState("load_preloaded_btn", condition = !is.null(input$preloaded_choice))
    })

    activate_dataset <- function(expr, meta, source, source_type, is_bundled_reference,
                                  geo_ids, declared_data_type = NA_character_) {
      dataset$staged_expr <- expr
      dataset$staged_meta <- meta
      dataset$staged_source <- source
      dataset$expr <- expr
      dataset$meta <- meta
      dataset$source <- source
      dataset$source_type <- source_type
      dataset$is_bundled_reference <- is_bundled_reference
      dataset$geo_ids <- geo_ids
      dataset$declared_data_type <- declared_data_type
      ## a newly loaded dataset never inherits the previous dataset's sealed hold-out
      dataset$reserved_ids <- character(0)
      dataset$reserved_expr <- NULL
      dataset$reserved_meta <- NULL
      dataset$reserved_info <- NULL
      ## Shiny reactiveValues skip invalidating dependent observers when a field is
      ## set to a value identical() to its current one. dataset$source is built from
      ## filenames alone, so two different re-uploads sharing a filename (a corrected
      ## re-upload saved under the same name; two GEO exports both saved as
      ## "expression.csv") produce an identical source string and silently never
      ## reset stale WGCNA/Candidate-gene results (RED finding, 2026-09-07 defense
      ## audit). load_id is a plain monotonic counter - guaranteed to differ from its
      ## previous value on every single activation regardless of content - so
      ## observers keyed on it always fire, independent of what dataset$source says.
      dataset$load_id <- (dataset$load_id %||% 0L) + 1L
      sum(duplicated(rownames(expr)))
    }

    ## ---- Sealed validation hold-out --------------------------------------------------------------
    reserve_msg <- reactiveVal(NULL)

    seal_samples <- function(ids, info) {
      expr <- dataset$expr; meta <- dataset$meta
      ids <- intersect(ids, colnames(expr))
      keep <- setdiff(colnames(expr), ids)
      dataset$reserved_expr <- expr[, ids, drop = FALSE]
      dataset$reserved_meta <- meta[match(ids, meta$sample), , drop = FALSE]
      dataset$reserved_ids <- ids
      dataset$reserved_info <- info
      dataset$expr <- expr[, keep, drop = FALSE]
      dataset$meta <- meta[match(keep, meta$sample), , drop = FALSE]
      ## dataset$staged_expr/staged_meta are what Preprocessing's "currently loaded
      ## dataset" fallback reads (`use_expr <- dataset$staged_expr %||% dataset$expr`,
      ## mod_preprocessing.R) and were never updated here before this fix - so
      ## ComBat/SVA/TMM/quantile-normalisation parameters could be estimated on the
      ## pre-seal matrix, INCLUDING the reserved hold-out samples, even though only
      ## the preprocessing OUTPUT was re-sealed afterward (RED finding, 2026-09-07
      ## defense audit: this directly contradicted the in-UI claim that reserved
      ## samples are "hidden from every step"). Keep staged_expr/staged_meta in
      ## lock-step with expr/meta so every downstream fallback also only ever sees
      ## the discovery subset.
      if (!is.null(dataset$staged_expr)) {
        staged_ids <- intersect(colnames(dataset$staged_expr), ids)
        staged_keep <- setdiff(colnames(dataset$staged_expr), staged_ids)
        dataset$staged_expr <- dataset$staged_expr[, staged_keep, drop = FALSE]
        if (!is.null(dataset$staged_meta)) {
          dataset$staged_meta <- dataset$staged_meta[match(staged_keep, dataset$staged_meta$sample), , drop = FALSE]
        }
      }
    }

    observeEvent(input$reserve_btn, {
      req(dataset$expr, dataset$meta)
      if (length(dataset$reserved_ids %||% character(0)) > 0) {
        reserve_msg(list(kind = "warn", text = "Validation samples are already reserved for this dataset. Release them first if you want a different split."))
        return()
      }
      frac <- (input$reserve_frac_pct %||% 30) / 100
      ids <- tryCatch(tx_reserve_validation_ids(dataset$meta, frac, seed = input$reserve_seed %||% 1234,
                                                stratify_sex = isTRUE(input$reserve_stratify_sex)),
                      error = function(e) e)
      if (inherits(ids, "error")) { reserve_msg(list(kind = "warn", text = conditionMessage(ids))); return() }
      info <- list(frac = frac, seed = input$reserve_seed %||% 1234, stratify_sex = isTRUE(input$reserve_stratify_sex),
                   reserved_at = Sys.time(), source_at_reservation = dataset$source,
                   was_bundled_reference = isTRUE(dataset$is_bundled_reference), resealed = 0L)
      seal_samples(ids, info)
      ## bundled shortcuts used all samples, so disable them for a reserved run
      dataset$is_bundled_reference <- FALSE
      base <- sub(DATASET_RESERVED_SUFFIX_RE, "", dataset$source %||% "Currently loaded dataset")
      dataset$source <- sprintf("%s [%d validation samples reserved]", base, length(ids))
      arthomix_provenance_push(arthomix_provenance_record(
        module = "mod_dataset_reserve_validation",
        checksum_input = list(reserved_ids = sort(ids), discovery_ids = sort(colnames(dataset$expr))),
        params = list(fraction = frac, seed = info$seed, stratify_sex = info$stratify_sex,
                      n_reserved = length(ids), n_discovery = ncol(dataset$expr),
                      reserved_groups = paste(sprintf("%s = %d", names(table(dataset$reserved_meta$group)), as.integer(table(dataset$reserved_meta$group))), collapse = ", "),
                      bundled_shortcuts_disabled = info$was_bundled_reference),
        seed = info$seed, packages = "caret"
      ))
      reserve_msg(list(kind = "ok", text = sprintf(
        "Reserved %d of %d samples (%s) as a sealed validation set, stratified by %s. Every sub-module now sees only the %d discovery samples. Score them later in Diagnostic Model → External Validation → \"Sealed validation samples\".%s",
        length(ids), length(ids) + ncol(dataset$expr),
        paste(sprintf("%s = %d", names(table(dataset$reserved_meta$group)), as.integer(table(dataset$reserved_meta$group))), collapse = ", "),
        if (info$stratify_sex) "group and sex" else "group", ncol(dataset$expr),
        if (info$was_bundled_reference) " Precomputed bundled-reference shortcuts are disabled for this run because they were computed on all samples." else "")))
    }, ignoreInit = TRUE)

    observeEvent(input$release_btn, {
      ids <- dataset$reserved_ids %||% character(0)
      if (length(ids) == 0) { reserve_msg(list(kind = "warn", text = "No validation samples are reserved.")); return() }
      info <- dataset$reserved_info
      rexpr <- dataset$reserved_expr; rmeta <- dataset$reserved_meta
      common_genes <- intersect(rownames(dataset$expr), rownames(rexpr))
      dataset$expr <- cbind(dataset$expr[common_genes, , drop = FALSE], rexpr[common_genes, , drop = FALSE])
      dataset$meta <- rbind(dataset$meta, rmeta[, colnames(dataset$meta), drop = FALSE])
      ## Mirror the restore into staged_expr/staged_meta - see the matching comment in
      ## seal_samples() for why these must always stay in lock-step with expr/meta.
      if (!is.null(dataset$staged_expr)) {
        common_genes_staged <- intersect(rownames(dataset$staged_expr), rownames(rexpr))
        dataset$staged_expr <- cbind(dataset$staged_expr[common_genes_staged, , drop = FALSE], rexpr[common_genes_staged, , drop = FALSE])
        if (!is.null(dataset$staged_meta)) {
          dataset$staged_meta <- rbind(dataset$staged_meta, rmeta[, colnames(dataset$staged_meta), drop = FALSE])
        }
      }
      dataset$reserved_ids <- character(0); dataset$reserved_expr <- NULL; dataset$reserved_meta <- NULL; dataset$reserved_info <- NULL
      if (isTRUE(info$was_bundled_reference)) dataset$is_bundled_reference <- TRUE
      dataset$source <- sub(DATASET_RESERVED_SUFFIX_RE, "", dataset$source %||% "Currently loaded dataset")
      reserve_msg(list(kind = "ok", text = sprintf("Released %d reserved samples back into the active dataset (%d samples). Results computed on the discovery subset have been cleared.", length(ids), ncol(dataset$expr))))
    }, ignoreInit = TRUE)

    ## Re-seal reserved samples if another sub-module reintroduces them into the active matrix.
    observeEvent(dataset$expr, {
      ids <- dataset$reserved_ids %||% character(0)
      if (length(ids) == 0 || is.null(dataset$expr)) return()
      back <- intersect(ids, colnames(dataset$expr))
      if (length(back) == 0) return()
      info <- dataset$reserved_info %||% list()
      info$resealed <- (info$resealed %||% 0L) + 1L
      info$resealed_from <- dataset$source
      seal_samples(ids, info)
      if (isTRUE(dataset$is_bundled_reference)) dataset$is_bundled_reference <- FALSE
      if (!grepl(DATASET_RESERVED_SUFFIX_RE, dataset$source %||% "")) {
        dataset$source <- sprintf("%s [%d validation samples reserved]", dataset$source %||% "Currently loaded dataset", length(dataset$reserved_ids))
      }
      reserve_msg(list(kind = "ok", text = sprintf(
        "The active matrix was replaced (%s) and contained %d reserved samples. They were sealed away again, now with the same processing as the discovery samples.",
        dataset$source, length(back))))
    })

    output$reserve_status <- renderUI({
      ids <- dataset$reserved_ids %||% character(0)
      msg <- reserve_msg()
      tagList(
        if (!is.null(msg)) div(class = "empty-note", icon(if (identical(msg$kind, "ok")) "check" else "triangle-exclamation"), msg$text),
        if (length(ids) > 0) {
          info <- dataset$reserved_info %||% list()
          rt <- table(dataset$reserved_meta$group)
          div(class = "empty-note", icon("lock"),
              sprintf("Sealed hold-out active: %d reserved samples (%s), seed %s, %d%% of the dataset%s. Reserved on: %s.%s",
                      length(ids), paste(sprintf("%s = %d", names(rt), as.integer(rt)), collapse = ", "),
                      format(info$seed %||% NA), round((info$frac %||% NA) * 100),
                      if (isTRUE(info$stratify_sex)) ", stratified by group and sex" else ", stratified by group",
                      info$source_at_reservation %||% "(unknown)",
                      if ((info$resealed %||% 0L) > 0) sprintf(" Re-sealed %d time(s) after the active matrix was reprocessed.", info$resealed) else ""))
        } else div(class = "empty-note", icon("circle-info"), "No validation samples reserved - every loaded sample is visible to every sub-module.")
      )
    })

    duplicate_feature_note <- function(n_dup) {
      if (n_dup == 0) return(NULL)
      div(class = "empty-note", icon("triangle-exclamation"),
          sprintf("%d duplicated feature identifier(s) found. All rows are kept here, but downstream steps keyed by row name (e.g. the Preprocessing merge tab) keep only the first occurrence of each. Rename duplicates in your source file if unintended.", n_dup))
    }

    observeEvent(input$load_preloaded_btn, {
      req(input$preloaded_choice)
      entry <- Find(function(d) d$id == input$preloaded_choice, PRELOADED_DATASETS)
      req(entry)
      d <- entry$load()
      n_dup <- activate_dataset(
        expr = d$expr, meta = d$meta, source = d$source,
        source_type = "preloaded",
        is_bundled_reference = identical(entry$id, default_dataset_entry$id),
        geo_ids = entry_geo_ids(entry)
      )
      output$preloaded_load_message <- renderUI(tagList(
        span(style = "color: var(--color-success); font-size: 13px; font-weight: 600;", icon("check"), " ",
             sprintf("Loaded %s: %s genes across %s samples. Every sub-module now runs on this dataset. Optionally go to Preprocessing and pick \"Currently loaded dataset\" to merge, normalise, or batch-correct it first.",
                      entry$label, format(nrow(d$expr), big.mark = ","), ncol(d$expr))),
        duplicate_feature_note(n_dup)
      ))
    })

    observeEvent(input$load_btn, {
      req(input$expr_file, input$meta_file, input$map_id, input$map_group)

      result <- tryCatch({
        expr <- expr_raw()
        checked <- tx_validate_expr_upload(expr, input$declared_data_type)
        validate(need(isTRUE(checked$ok), checked$error))
        expr <- checked$mat
        meta <- meta_raw()
        meta$sample <- as.character(meta[[input$map_id]])
        meta$group  <- as.character(meta[[input$map_group]])
        meta$sex    <- if (!identical(input$map_sex, "(none)")) as.character(meta[[input$map_sex]]) else NA_character_
        meta$batch  <- if (!identical(input$map_batch, "(none)")) as.character(meta[[input$map_batch]]) else NA_character_

        common <- intersect(colnames(expr), meta$sample)
        validate(need(
          length(common) >= 4,
          "Fewer than 4 sample IDs in the expression matrix match the metadata sample-ID column. Check the column mapping."
        ))
        dup_ids <- unique(meta$sample[duplicated(meta$sample) & meta$sample %in% common])
        validate(need(
          length(dup_ids) == 0,
          sprintf("The metadata sample-ID column has duplicate value(s) matching expression-matrix columns: %s. Make sample IDs unique before loading.", paste(dup_ids, collapse = ", "))
        ))

        expr <- expr[, common, drop = FALSE]
        meta <- meta[match(common, meta$sample), , drop = FALSE]
        list(expr = expr, meta = meta, type_note = checked$note)
      }, error = function(e) e)

      if (inherits(result, "error")) {
        output$load_message <- renderUI(
          div(class = "empty-note", icon("triangle-exclamation"), paste("Could not load this dataset:", conditionMessage(result)))
        )
      } else {
        source_label <- paste0("Uploaded dataset: ", input$expr_file$name, " + ", input$meta_file$name)
        n_dup <- activate_dataset(
          expr = result$expr, meta = result$meta, source = source_label,
          source_type = "uploaded", is_bundled_reference = FALSE,
          geo_ids = character(0), declared_data_type = input$declared_data_type
        )
        output$load_message <- renderUI(
          tagList(
            div(class = "empty-note", icon("check"),
                sprintf("Loaded %s genes across %s samples, as-is. Every sub-module now runs on this dataset.",
                        format(nrow(result$expr), big.mark = ","), ncol(result$expr))),
            duplicate_feature_note(n_dup),
            if (!is.null(result$type_note) && nzchar(result$type_note)) div(class = "empty-note", icon("triangle-exclamation"), result$type_note)
          )
        )
      }
    })

    geo_fetch_result <- eventReactive(input$geo_fetch_btn, {
      if (!requireNamespace("GEOquery", quietly = TRUE)) {
        return(simpleError("The GEOquery package is not installed in this deployment. Install it with BiocManager::install(\"GEOquery\") to enable fetching by GEO accession, or use \"Upload your own data\" instead."))
      }
      acc <- toupper(trimws(input$geo_accession %||% ""))
      if (!grepl("^GSE[0-9]+$", acc)) {
        return(simpleError("Enter a valid GEO Series accession, e.g. GSE12345."))
      }
      tryCatch({
        gpl_skipped <- FALSE
        gset <- tryCatch(
          suppressMessages(GEOquery::getGEO(acc, GSEMatrix = TRUE)),
          error = function(e) {
            if (!grepl("series_data_table_begin", conditionMessage(e), fixed = TRUE)) stop(e)
            gpl_skipped <<- TRUE
            suppressMessages(GEOquery::getGEO(acc, GSEMatrix = TRUE, getGPL = FALSE))
          }
        )
        if (!is.list(gset) || length(gset) == 0) {
          stop(paste(acc, "returned no series matrix from GEO - check the accession is a Series (GSExxxxx), not a Sample (GSM) or Platform (GPL) ID."))
        }
        list(acc = acc, platforms = gset, gpl_skipped = gpl_skipped)
      }, error = function(e) e)
    })

    output$geo_platform_ui <- renderUI({
      res <- geo_fetch_result()
      req(res); req(!inherits(res, "error"))
      if (length(res$platforms) <= 1) return(NULL)
      selectInput(ns("geo_platform_choice"), "This series spans multiple platforms - pick one",
                  choices = names(res$platforms), width = "100%")
    })

    geo_eset <- reactive({
      res <- geo_fetch_result()
      req(res); req(!inherits(res, "error"))
      if (length(res$platforms) > 1) {
        req(input$geo_platform_choice)
        res$platforms[[input$geo_platform_choice]]
      } else {
        res$platforms[[1]]
      }
    })

    geo_expr_meta <- reactive({
      res <- geo_fetch_result()
      req(res); req(!inherits(res, "error"))
      if (length(res$platforms) > 1 && is.null(input$geo_platform_choice)) return(NULL)
      tryCatch({
        eset <- geo_eset()
        ex <- Biobase::exprs(eset)
        if (nrow(ex) == 0 || ncol(ex) == 0) {
          stop("This GEO series has no expression matrix in its series matrix file. This is common for RNA-seq series that only deposit raw counts as supplementary files. Download that file from GEO and use \"Upload your own data\" instead.")
        }
        collapsed <- tryCatch(collapse_probes_to_genes(eset), error = arthomix_null_on_error)
        used_collapse <- !is.null(collapsed) && nrow(collapsed) > 0 && isTRUE(attr(collapsed, "collapsed"))
        expr <- if (used_collapse) collapsed else ex
        list(expr = expr, meta = as.data.frame(Biobase::pData(eset)),
             collapsed = used_collapse, platform = Biobase::annotation(eset))
      }, error = function(e) e)
    })

    output$geo_fetch_status <- renderUI({
      req(input$geo_fetch_btn)
      res <- geo_fetch_result()
      if (inherits(res, "error")) {
        return(div(class = "empty-note", icon("triangle-exclamation"), paste("Could not fetch from GEO:", conditionMessage(res))))
      }
      em <- geo_expr_meta()
      req(em) # NULL = still waiting on a platform pick above; render nothing yet, not an error
      if (inherits(em, "error")) {
        return(div(class = "empty-note", icon("triangle-exclamation"), conditionMessage(em)))
      }
      tagList(
        div(class = "empty-note", icon("circle-info"),
            sprintf("Fetched %s (%s): %s samples x %s %s. Map the columns below, then click \"Load this dataset\".",
                    res$acc, em$platform, ncol(em$expr), format(nrow(em$expr), big.mark = ","),
                    if (em$collapsed) "genes" else "probes/features")),
        if (!em$collapsed) {
          div(class = "empty-note", icon("triangle-exclamation"),
              if (isTRUE(res$gpl_skipped))
                "GEO's platform-annotation file couldn't be fetched (NCBI may be rate-limiting requests). Loaded at probe/feature-ID level without gene symbols. Try again later, or use Preprocessing's probe-collapse step afterward."
              else
                "No gene-symbol annotation found for this platform - left at probe/feature-ID level. You can still load it as-is, or use Preprocessing's own probe-collapse step afterward.")
        }
      )
    })

    output$geo_column_mapping <- renderUI({
      em <- geo_expr_meta()
      req(em); req(!inherits(em, "error"))
      cols <- colnames(em$meta)
      geo_group_guess <- guess_col(
        cols,
        exact = c("group", "diagnosis", "disease", "condition", "phenotype"),
        contains = c("group", "diagnosis", "disease", "condition", "phenotype", "characteristics_ch1")
      )
      tagList(
        radioButtons(ns("geo_declared_data_type"), "Data type", inline = TRUE,
                     choices = c("Raw counts" = "raw", "Normalized (TPM/FPKM/CPM)" = "normalized", "Already log-transformed" = "logtransformed"),
                     selected = "normalized"),
        div(class = "empty-note", icon("circle-info"),
            "GEO series matrices are usually already normalized (microarray) or summarized counts - check the series' own \"Data processing\" description on GEO if unsure."),
        selectInput(ns("geo_map_group"), "Group / diagnosis column", choices = cols,
                    selected = geo_group_guess, selectize = FALSE),
        unconfident_guess_note(geo_group_guess, "group/diagnosis column"),
        selectInput(ns("geo_map_sex"), "Sex column (optional)", choices = c("(none)", cols),
                    selected = guess_col(cols, c("sex", "gender"), fallback = "(none)"),
                    selectize = FALSE),
        selectInput(ns("geo_map_batch"), "Batch column (optional)", choices = c("(none)", cols),
                    selected = guess_col(cols, c("batch", "cohort", "platform", "dataset"), fallback = "(none)"),
                    selectize = FALSE)
      )
    })

    observe({
      em <- tryCatch(geo_expr_meta(), error = function(e) e)
      ready <- !is.null(em) && !inherits(em, "error") && !is.null(input$geo_map_group)
      if (isTRUE(ready)) shinyjs::enable("geo_load_btn") else shinyjs::disable("geo_load_btn")
    })

    observeEvent(input$geo_load_btn, {
      req(input$geo_map_group)

      result <- tryCatch({
        em <- geo_expr_meta()
        validate(need(!inherits(em, "error"), "No GEO data fetched yet."))
        checked <- tx_validate_expr_upload(em$expr, input$geo_declared_data_type)
        validate(need(isTRUE(checked$ok), checked$error))
        expr <- checked$mat
        meta <- em$meta
        meta$sample <- rownames(meta)
        meta$group  <- as.character(meta[[input$geo_map_group]])
        meta$sex    <- if (!identical(input$geo_map_sex, "(none)")) as.character(meta[[input$geo_map_sex]]) else NA_character_
        meta$batch  <- if (!identical(input$geo_map_batch, "(none)")) as.character(meta[[input$geo_map_batch]]) else NA_character_

        common <- intersect(colnames(expr), meta$sample)
        validate(need(
          length(common) >= 4,
          "Fewer than 4 sample IDs matched between the fetched expression matrix and metadata."
        ))
        dup_ids <- unique(meta$sample[duplicated(meta$sample) & meta$sample %in% common])
        validate(need(
          length(dup_ids) == 0,
          sprintf("The metadata sample-ID column has duplicate value(s) matching expression-matrix columns: %s. Make sample IDs unique before loading.", paste(dup_ids, collapse = ", "))
        ))
        expr <- expr[, common, drop = FALSE]
        meta <- meta[match(common, meta$sample), , drop = FALSE]
        acc <- geo_fetch_result()$acc
        label <- sprintf("%s (%s, %s)", acc, em$platform,
                          if (em$collapsed) "collapsed to genes" else "probe-level, raw")
        list(expr = expr, meta = meta, label = label, acc = acc, type_note = checked$note)
      }, error = function(e) e)

      if (inherits(result, "error")) {
        output$load_message <- renderUI(
          div(class = "empty-note", icon("triangle-exclamation"), paste("Could not load this GEO dataset:", conditionMessage(result)))
        )
      } else {
        source_label <- paste0("NCBI GEO: ", result$label)
        n_dup <- activate_dataset(
          expr = result$expr, meta = result$meta, source = source_label,
          source_type = "geo", is_bundled_reference = FALSE, geo_ids = result$acc,
          declared_data_type = input$geo_declared_data_type
        )
        n_samples <- ncol(result$expr)
        wgcna_note <- if (n_samples < 15) {
          " Note: WGCNA needs at least 15 samples to detect modules reliably. Its tab stays blank until a larger dataset is loaded (e.g. merged with another GEO series or the preloaded cohort)."
        } else {
          ""
        }
        probe_note <- if (!isTRUE(em$collapsed)) " Note: no gene-symbol annotation was found for this platform, so rows are raw probe/feature IDs, not gene symbols." else ""
        output$load_message <- renderUI(tagList(
          div(class = "empty-note", icon("check"),
              sprintf("Loaded %s genes across %s samples. Every sub-module now runs on this GEO dataset.%s%s",
                      format(nrow(result$expr), big.mark = ","), n_samples, wgcna_note, probe_note)),
          duplicate_feature_note(n_dup),
          if (!is.null(result$type_note) && nzchar(result$type_note)) div(class = "empty-note", icon("triangle-exclamation"), result$type_note)
        ))
      }
    })

  })
}
