## R/multiomics/mod_multi_pathway.R
## Submodule: Pathways - a live GO/KEGG/Reactome/WikiPathways pathway-
## enrichment engine (ORA + GSEA) over either the app's preloaded multi-omics
## candidate panel (DIABLO/SNF-selected genes, methylation-linked genes via
## CpG->gene mapping, real per-feature effect sizes from Table42/45) or a
## user-uploaded gene/CpG-level table whose structure is auto-detected, never
## assumed. Nothing below is computed - no enrichment table, plot, or pathway
## map - until the blue "Run Pathway Analysis" button is clicked; only
## filters/settings and the Detected Data summary are visible before that.
## All heavy lifting is in multiomics_pathway_helpers.R/_plots.R; this file
## is UI wiring only (same split as mod_multi_concordance.R).

mod_multi_pathway_config <- list(
  id = "pathway", title = "Pathways", icon = "sitemap", group = "Interpretation",
  description = "GO/KEGG/Reactome/WikiPathways/Hallmark enrichment on candidate biomarkers or your own upload, with real pathway maps."
)

MP_ORA_FIELDS <- c("id_col", "id_type", "effect_col", "pvalue_col", "fdr_col", "direction_col", "omics_col")

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------

mod_multi_pathway_ui <- function(id) {
  ns <- NS(id)
  tagList(uiOutput(ns("active_dataset_banner")), fluidRow(
    column(
      4,
      box(
        width = NULL, title = "1. Data source", status = "primary", solidHeader = FALSE,
        radioButtons(ns("data_source"), NULL, choices = c("Preloaded candidate panel" = "preloaded", "Upload my own data" = "upload"), selected = "preloaded"),
        selectInput(ns("sex"), "Sex", choices = c("Both (pooled)" = "", "Female" = "female", "Male" = "male")),
        p(class = "submodule-desc", "For an upload, this filters on a Sex column you map below (Detected data) - optional; only applies if your file actually carries one."),
        conditionalPanel(condition = sprintf("input['%s'] == 'preloaded'", ns("data_source")),
                          selectInput(ns("cohort"), "Cohort", choices = MP_PRELOADED_COHORTS, width = "100%"),
                          uiOutput(ns("layer_pick_ui")),
                          selectInput(ns("source"), "Candidate source", choices = MCC_BIOMARKER_SOURCES, selected = "All candidates"),
                          selectizeInput(ns("custom_genes"), "Custom genes (optional)", choices = NULL, multiple = TRUE, options = list(create = TRUE, placeholder = "Type a gene symbol and press Enter")),
                          selectizeInput(ns("custom_cpgs"), "Custom CpGs (optional)", choices = NULL, multiple = TRUE, options = list(create = TRUE, placeholder = "Type a CpG ID and press Enter"))),
        conditionalPanel(condition = sprintf("input['%s'] == 'upload'", ns("data_source")),
                          fileInput(ns("upload_file"), "File (CSV/TSV/TXT/XLSX)", accept = c(".csv", ".tsv", ".txt", ".xlsx")))
      ),
      conditionalPanel(condition = sprintf("input['%s'] == 'upload'", ns("data_source")),
                        box(width = NULL, title = "2. Detected data", status = "primary", solidHeader = FALSE, uiOutput(ns("detected_ui")))),
      box(
        width = NULL, title = "3. Filters", status = "primary", solidHeader = FALSE, collapsible = TRUE,
        checkboxGroupInput(ns("database"), "Database", choices = stats::setNames(names(MP_DATABASES), vapply(MP_DATABASES, function(d) d$label, character(1))), selected = "GO_BP"),
        uiOutput(ns("database_availability_note")),
        radioButtons(ns("method"), "Analysis", choices = c("Over-representation (ORA)" = "ORA", "Gene-set enrichment (GSEA)" = "GSEA", "Topology (Reactome only)" = "Topology"), selected = "ORA"),
        selectInput(ns("species"), "Species", choices = c("Homo sapiens" = "human"), selected = "human"),
        div(class = "empty-note", icon("circle-info"), "Human (Homo sapiens) only - other species are rejected."),
        conditionalPanel(condition = sprintf("input['%s'] != 'GSEA'", ns("method")),
                          fluidRow(column(6, numericInput(ns("pval_cut"), "P-value cutoff", value = 1, min = 0, max = 1, step = 0.05)),
                                   column(6, numericInput(ns("fdr_cut"), "FDR cutoff", value = 0.25, min = 0, max = 1, step = 0.01)))),
        conditionalPanel(condition = sprintf("input['%s'] == 'GSEA'", ns("method")),
                          selectInput(ns("ranking_method"), "Ranking", choices = c("log2 fold change" = "log2fc", "Signed -log10(P)" = "signed_neglog10p")),
                          numericInput(ns("gsea_fdr_cut"), "FDR cutoff", value = 0.25, min = 0, max = 1, step = 0.01)),
        fluidRow(column(6, numericInput(ns("min_size"), "Min gene-set size", value = 5, min = 1, step = 1)),
                 column(6, numericInput(ns("max_size"), "Max gene-set size", value = 500, min = 5, step = 10))),
        selectInput(ns("background"), "Background / universe", choices = c(
          "Auto (measured features in active dataset)" = "auto_experimental",
          "Preloaded cohort candidate-gene list (small - not genome-wide)" = "preloaded_universe",
          "Uploaded file's own identifier list" = "uploaded_background",
          "Entire selected database (no experimental universe)" = "entire_database"
        ), selected = "auto_experimental"),
        conditionalPanel(condition = sprintf("input['%s'] == 'entire_database'", ns("background")),
                          div(class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"),
                              "No experimental universe supplied - results will be labeled accordingly."))
      ),
      box(width = NULL, title = "4. Run", status = "primary", solidHeader = FALSE,
          uiOutput(ns("readiness_ui")),
          actionButton(ns("run_btn"), "Run Pathway Analysis", icon = icon("play"), class = "btn-primary btn-lg", width = "100%")),
      br(), br()
    ),
    column(
      8,
      tabsetPanel(
        id = ns("tabs"), type = "tabs",
        tabPanel("Enrichment", br(), uiOutput(ns("enrichment_ui"))),
        tabPanel("Pathway Map", br(), uiOutput(ns("map_ui"))),
        tabPanel("Gene Sets", br(), uiOutput(ns("genesets_ui"))),
        tabPanel("Results", br(), uiOutput(ns("results_ui")))
      )
    )
  ))
}

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------

mod_multi_pathway_server <- function(id, multi_dataset = NULL, multi_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    output$active_dataset_banner <- renderUI(multi_active_dataset_banner(multi_dataset))

    ## =========================================================================
    ## Preloaded-path layer pickers (spec 2/3) - populated only from what
    ## multi_dataset actually contains, same helpers mod_multi_concordance.R uses.
    ## =========================================================================

    output$layer_pick_ui <- renderUI({
      layers <- multi_dataset$layers %||% list()
      if (length(layers) == 0) return(div(class = "empty-note", icon("circle-info"), "No Active Multi-Omics Dataset loaded - the preloaded candidate panel (DIABLO/Biomarker Discovery/Patient Stratification) is still usable without one."))
      expr_cand <- mcc_layer_candidates(multi_dataset, "rnaseq")
      meth_cand <- mcc_layer_candidates(multi_dataset, "methylation")
      tagList(
        selectInput(ns("expr_layer"), "Expression layer (for SNF candidates/background)", choices = c("(none)", names(layers)), selected = mcc_default_layer(expr_cand, layers) %||% "(none)"),
        selectInput(ns("meth_layer"), "Methylation layer (for SNF candidates/background)", choices = c("(none)", names(layers)), selected = mcc_default_layer(meth_cand, layers) %||% "(none)")
      )
    })
    expr_layer_val <- reactive(if (identical(input$expr_layer, "(none)")) NULL else input$expr_layer)
    meth_layer_val <- reactive(if (identical(input$meth_layer, "(none)")) NULL else input$meth_layer)

    output$database_availability_note <- renderUI({
      unavail <- Filter(function(nm) !isTRUE(MP_DATABASES[[nm]]$available), names(MP_DATABASES))
      if (length(unavail) == 0) return(NULL)
      div(class = "empty-note", icon("triangle-exclamation"), sprintf("Unavailable in this deployment: %s.", paste(vapply(unavail, function(nm) MP_DATABASES[[nm]]$label, character(1)), collapse = ", ")))
    })

    ## =========================================================================
    ## Upload path (spec 2, 10, 15) - detect, report, never silently accept.
    ## =========================================================================

    upload_raw <- reactiveVal(NULL)
    upload_detected <- reactiveVal(NULL)

    observeEvent(input$upload_file, {
      res <- mp_read_upload_table(input$upload_file$datapath, input$upload_file$name)
      if (!isTRUE(res$ok)) { showNotification(res$error, type = "error"); upload_raw(NULL); upload_detected(NULL); return() }
      upload_raw(res$df)
      upload_detected(mp_detect_upload(res$df))
    })

    output$detected_ui <- renderUI({
      df <- upload_raw(); det <- upload_detected()
      if (is.null(df) || is.null(det)) return(multi_empty_state("Upload a CSV/TSV/TXT/XLSX file to see what was detected."))
      d <- det$detected
      col_choices <- c("(none)", colnames(df))
      pick <- function(field_id, label, default) selectInput(ns(field_id), label, choices = col_choices, selected = default %||% "(none)")
      tagList(
        p(class = "submodule-desc", sprintf("%s - %s rows, %s columns.", input$upload_file$name, format(d$n_rows, big.mark = ","), d$n_cols)),
        tags$ul(
          tags$li(sprintf("Identifier type: %s", d$id_type)),
          tags$li(sprintf("Ranked list: %s", if (d$ranked) "Yes (GSEA available)" else "No (ORA only)")),
          tags$li(sprintf("Species (guess): %s", d$species_guess)),
          tags$li(sprintf("Omics (guess): %s", d$omics_guess))
        ),
        if (length(det$warnings) > 0) div(class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"), tags$ul(lapply(det$warnings, tags$li))) else NULL,
        pick("map_id_col", "Identifier column", d$id_col),
        pick("map_effect_col", "Effect size / log2FC / statistic column", d$effect_col),
        pick("map_pvalue_col", "P-value column", d$pvalue_col),
        pick("map_fdr_col", "Adjusted P / FDR column", d$fdr_col),
        pick("map_direction_col", "Direction column", d$direction_col),
        pick("map_omics_col", "Omics-type column", d$omics_col),
        pick("map_sex_col", "Sex column (optional)", d$sex_col),
        actionButton(ns("confirm_mapping"), "Confirm mapping", icon = icon("check"), class = "btn-sm"),
        uiOutput(ns("confirm_status_ui"))
      )
    })

    upload_confirmed <- reactiveVal(NULL)
    observeEvent(input$confirm_mapping, {
      df <- req(upload_raw()); det <- req(upload_detected())
      mapping <- list(id_col = input$map_id_col, id_type = det$detected$id_type, effect_col = input$map_effect_col,
                       pvalue_col = input$map_pvalue_col, fdr_col = input$map_fdr_col, direction_col = input$map_direction_col,
                       omics_col = input$map_omics_col, omics_guess = det$detected$omics_guess, sex_col = input$map_sex_col)
      res <- mp_confirm_upload_mapping(df, mapping)
      if (!isTRUE(res$ok)) { showNotification(res$error, type = "error"); upload_confirmed(NULL); return() }
      upload_confirmed(res)
      showNotification(res$source_detail, type = "message")
    })
    output$confirm_status_ui <- renderUI({
      if (is.null(upload_confirmed())) return(div(class = "empty-note", icon("circle-info"), "Not confirmed yet."))
      div(class = "empty-note", icon("circle-check"), upload_confirmed()$source_detail)
    })

    ## =========================================================================
    ## Built input (spec 2/3) - dispatches by data source; drives readiness
    ## and the GSEA ranking-column choices, but triggers NO enrichment itself.
    ## =========================================================================

    built_input <- reactive({
      if (identical(input$data_source, "upload")) {
        bi <- upload_confirmed()
        ## Sex is preloaded-only in mp_build_preloaded_input() (it filters a
        ## precomputed table's own `sex` column) - for an upload, the same
        ## filter is applied here, post-hoc, against the optional `sex`
        ## column mp_confirm_upload_mapping() now carries through when the
        ## user's own file has one (map_sex_col above). Never invents a sex
        ## value when the uploaded file didn't provide one.
        if (!is.null(bi) && isTRUE(bi$ok) && nzchar(input$sex %||% "") && "sex" %in% colnames(bi$df) && any(!is.na(bi$df$sex))) {
          bi$df <- bi$df[!is.na(bi$df$sex) & tolower(trimws(bi$df$sex)) == input$sex, , drop = FALSE]
          bi$source_detail <- sprintf("%s Filtered to Sex = %s.", bi$source_detail, input$sex)
        }
        return(bi)
      }
      mp_build_preloaded_input(multi_results, multi_dataset, expr_layer_val(), meth_layer_val(), input$cohort,
                                if (nzchar(input$sex %||% "")) input$sex else character(0),
                                input$custom_genes %||% character(0), input$custom_cpgs %||% character(0))
    })

    output$readiness_ui <- renderUI({
      bi <- tryCatch(built_input(), error = function(e) NULL)
      row <- function(label, ok, detail) tags$tr(tags$td(label), tags$td(tags$strong(style = sprintf("color:%s;", if (isTRUE(ok)) ARTHOMIX_STATUS$good else ARTHOMIX_STATUS$warning), if (isTRUE(ok)) "Ready" else "Not ready")), tags$td(style = "color:var(--color-ink-muted,#898781);", detail))
      pool_filtered <- if (!is.null(bi) && isTRUE(bi$ok) && identical(input$data_source, "preloaded")) mcc_filter_source(bi$df, input$source) else if (!is.null(bi) && isTRUE(bi$ok)) bi$df else NULL
      tags$table(class = "table table-condensed", style = "font-size:0.85em;", tags$tbody(
        row("Candidate data", !is.null(bi) && isTRUE(bi$ok), if (!is.null(bi) && isTRUE(bi$ok)) sprintf("%d feature(s) - %s", nrow(pool_filtered %||% bi$df), bi$source_detail) else "Not yet available."),
        row("Database selected", length(input$database %||% character(0)) > 0, paste(input$database, collapse = ", ") %||% "None selected."),
        row("Ranked list (GSEA)", identical(input$method, "GSEA") == FALSE || (!is.null(bi) && isTRUE(bi$ok) && sum(!is.na(pool_filtered$expr_logFC %||% bi$df$expr_logFC)) >= 10), if (identical(input$method, "GSEA")) sprintf("%d features with a numeric effect size.", sum(!is.na(pool_filtered$expr_logFC %||% bi$df$expr_logFC %||% NA))) else "Not required for ORA/Topology.")
      ))
    })

    ## =========================================================================
    ## RUN (spec 1, 19-20) - the only place computation happens.
    ## =========================================================================

    has_run <- reactiveVal(FALSE)
    result <- reactiveVal(NULL)

    observeEvent(input$run_btn, {
      bi <- built_input()
      if (is.null(bi) || !isTRUE(bi$ok)) { showNotification("No candidate data available to analyze - see the readiness panel.", type = "error"); return() }
      if (length(input$database %||% character(0)) == 0) { showNotification("Select at least one database.", type = "error"); return() }

      withProgress(message = "Running pathway analysis...", value = 0.1, {
        df <- if (identical(input$data_source, "preloaded")) mcc_filter_source(bi$df, input$source) else bi$df
        if (is.null(df) || nrow(df) == 0) df <- bi$df
        mapped <- mp_harmonize_identifiers(df)
        summ <- mp_mapping_summary(mapped$df)

        species <- mp_infer_species(df$id_type[1] %||% "unknown", df$feature)
        if (!identical(species, "Homo sapiens")) {
          showNotification(sprintf("Detected species: %s - Homo sapiens only in this build.", species), type = "error")
          result(list(ok = FALSE, error = sprintf("Detected species %s - this deployment supports Homo sapiens only.", species)))
          has_run(TRUE); return()
        }

        genes_entrez <- unique(stats::na.omit(mapped$df$entrez_id))
        univ <- mp_resolve_universe(input$background, multi_dataset, expr_layer_val(), meth_layer_val(),
                                     uploaded_universe_ids = if (identical(input$data_source, "upload")) bi$df$feature else NULL)
        if (!isTRUE(univ$ok) && !identical(input$background, "entire_database")) {
          showNotification(univ$error %||% "Could not resolve a background/universe.", type = "error")
          result(list(ok = FALSE, error = univ$error)); has_run(TRUE); return()
        }
        ov <- mp_validate_ora_inputs(genes_entrez, univ$universe_entrez)
        if (identical(input$method, "GSEA")) {
          rk <- mp_build_ranked_vector(mapped$df, mapped$df$entrez_id, input$ranking_method)
        }
        params <- list(pvalueCutoff = input$pval_cut %||% 1, minGSSize = input$min_size %||% 5, maxGSSize = input$max_size %||% 500,
                        fdrCutoff = if (identical(input$method, "GSEA")) input$gsea_fdr_cut %||% 0.25 else input$fdr_cut %||% 0.25)

        rows <- list(); errs <- character(0)
        for (db in input$database) {
          incProgress(1 / max(1, length(input$database)))
          dv <- mp_validate_database_choice(db)
          if (!isTRUE(dv$ok)) { errs <- c(errs, dv$error); next }
          r <- if (identical(input$method, "GSEA")) {
            if (!isTRUE(rk$ok)) { errs <- c(errs, rk$error); next }
            mp_run_gsea(db, rk$vec, params)
          } else if (identical(input$method, "Topology")) {
            if (!identical(db, "Reactome")) { errs <- c(errs, sprintf("Topology is only available for Reactome - skipped %s.", db)); next }
            if (!isTRUE(ov$ok)) { errs <- c(errs, ov$error); next }
            mp_run_topology(genes_entrez, univ$universe_entrez, params)
          } else {
            if (!isTRUE(ov$ok)) { errs <- c(errs, ov$error); next }
            mp_run_ora(db, genes_entrez, univ$universe_entrez, params)
          }
          if (isTRUE(r$ok) && !is.null(r$df)) rows[[db]] <- r$df else errs <- c(errs, sprintf("%s: %s", db, r$error %||% "no result"))
        }

        if (length(rows) == 0) {
          result(list(ok = FALSE, error = paste(unique(errs), collapse = " | ")))
          has_run(TRUE); return()
        }
        tab <- do.call(rbind, rows)
        fdr_thresh <- if (identical(input$method, "GSEA")) input$gsea_fdr_cut %||% 0.25 else input$fdr_cut %||% 0.25
        ## Every tested term is kept and shown, ranked by p.adjust - a small
        ## input gene list (e.g. a 21-gene candidate panel) naturally clears
        ## FDR for only a handful of terms; previously this dropped every
        ## other tested term from view entirely, which is why KEGG/etc. could
        ## show just one pathway even though the database itself was queried
        ## in full. `significant` flags which rows actually clear FDR - nothing
        ## is hidden, nothing is relabeled as significant that isn't.
        tab <- tab[order(tab$p.adjust, tab$pvalue), , drop = FALSE]
        tab$significant <- !is.na(tab$p.adjust) & tab$p.adjust < fdr_thresh
        n_sig <- sum(tab$significant)
        tab_show <- mp_build_evidence_tracks(tab, df)
        tab_show$concordance <- vapply(seq_len(nrow(tab_show)), function(i) {
          overlap <- unique(trimws(unlist(strsplit(tab_show$geneID[i] %||% "", "/"))))
          sub <- df[toupper(df$gene_symbol) %in% toupper(overlap) | toupper(df$feature) %in% toupper(overlap), , drop = FALSE]
          if (nrow(sub) == 0) return("Insufficient information")
          mp_concordance_direction(as.list(sub[1, , drop = FALSE]))
        }, character(1))

        meta <- mp_build_metadata(input$database, input$method, "Homo sapiens", univ$universe_label, if (identical(input$method, "GSEA")) input$ranking_method else NULL,
                                   summ$n_input, summ$n_mapped, fdr_thresh, input$min_size, input$max_size)

        result(list(ok = TRUE, table = tab_show, mapping = summ, mapping_df = mapped$df, cpg_map = mapped$cpg_map,
                     input_df = df, metadata = meta, warnings = errs, fdr_thresh = fdr_thresh,
                     method = input$method, background_label = univ$universe_label))
        has_run(TRUE)
      })
    })

    res_ok <- function() { r <- result(); if (!is.null(r) && isTRUE(r$ok)) r else NULL }

    ## =========================================================================
    ## Enrichment tab
    ## =========================================================================

    output$enrichment_ui <- renderUI({
      if (!has_run()) return(box(width = NULL, title = "Enrichment", status = "primary", solidHeader = FALSE, div(class = "empty-note", icon("circle-info"), "Not run yet. Set filters on the left, then click \"Run Pathway Analysis\".")))
      r <- result()
      if (is.null(r) || !isTRUE(r$ok)) return(box(width = NULL, title = "Enrichment", status = "primary", solidHeader = FALSE, div(class = "empty-note", style = "border-color: var(--color-danger, #e34948);", icon("circle-xmark"), r$error %||% "Analysis failed.")))
      tagList(
        div(class = "empty-note", icon("circle-check"), sprintf("%s term(s) tested across %s - %s significant at FDR < %s.", format(nrow(r$table), big.mark = ","), paste(unique(r$table$source), collapse = ", "), sum(r$table$significant %||% FALSE), r$fdr_thresh)),
        if (length(r$warnings) > 0) div(class = "empty-note", style = "border-color: var(--color-warning, #eda100);", icon("triangle-exclamation"), tags$ul(lapply(unique(r$warnings), tags$li))) else NULL,
        fluidRow(
          column(6, box(width = NULL, title = "Dot plot", status = "primary", solidHeader = FALSE,
                         sliderInput(ns("dot_top_n"), "Top N", min = 5, max = 40, value = 20, step = 5),
                         multi_plot_or_empty(function() mp_dot_plot(r$table, if (identical(r$method, "GSEA")) "GSEA" else "ORA", input$dot_top_n %||% 20), ns("dot_plot"), height = "420px"),
                         downloadButton(ns("dl_dot_png"), "Download (PNG)", class = "btn-sm"))),
          column(6, box(width = NULL, title = "Bar plot", status = "primary", solidHeader = FALSE,
                         selectInput(ns("bar_sort"), "Sort by", choices = c("FDR", "p.adjust", "NES", "Count")),
                         multi_plot_or_empty(function() mp_bar_plot(r$table, 20, input$bar_sort %||% "FDR"), ns("bar_plot"), height = "420px"),
                         downloadButton(ns("dl_bar_png"), "Download (PNG)", class = "btn-sm")))
        ),
        box(width = NULL, title = "Pathway x Omics evidence heatmap", status = "primary", solidHeader = FALSE,
            multi_plot_or_empty(function() mp_omics_heatmap(r$table), ns("heatmap_plot"), "Not enough evidence-scored pathways for a heatmap.", height = "460px")),
        box(width = NULL, title = "Gene-Pathway network", status = "primary", solidHeader = FALSE,
            uiOutput(ns("network_pathway_pick_ui")),
            multi_plot_or_empty(function() mp_gene_pathway_network(r$table, r$input_df, input$network_pathway_pick), ns("network_plot"), "Select pathway(s), or too few overlapping genes for a readable network.", height = "460px"),
            downloadButton(ns("dl_network_png"), "Download (PNG)", class = "btn-sm")),
        box(width = NULL, title = "Enrichment table", status = "primary", solidHeader = FALSE,
            div(class = "table-toolbar", downloadButton(ns("dl_enrich_csv"), "Download table (CSV)", class = "btn-sm")),
            DT::dataTableOutput(ns("enrich_table")))
      )
    })
    output$dot_plot <- renderPlot(mp_dot_plot(req(res_ok())$table, if (identical(req(res_ok())$method, "GSEA")) "GSEA" else "ORA", input$dot_top_n %||% 20))
    output$dl_dot_png <- multi_png_download(function() mp_dot_plot(req(res_ok())$table, if (identical(req(res_ok())$method, "GSEA")) "GSEA" else "ORA", input$dot_top_n %||% 20), function() "pathway_dotplot.png")
    output$bar_plot <- renderPlot(mp_bar_plot(req(res_ok())$table, 20, input$bar_sort %||% "FDR"))
    output$dl_bar_png <- multi_png_download(function() mp_bar_plot(req(res_ok())$table, 20, input$bar_sort %||% "FDR"), function() "pathway_barplot.png")
    output$heatmap_plot <- renderPlot(mp_omics_heatmap(req(res_ok())$table))

    output$network_pathway_pick_ui <- renderUI({
      r <- req(res_ok())
      selectizeInput(ns("network_pathway_pick"), "Restrict network to pathway(s)", choices = stats::setNames(r$table$ID, r$table$Description), multiple = TRUE, options = list(placeholder = "Top 8 by adjusted P (default)"))
    })
    output$network_plot <- renderPlot(mp_gene_pathway_network(req(res_ok())$table, req(res_ok())$input_df, input$network_pathway_pick))
    output$dl_network_png <- multi_png_download(function() mp_gene_pathway_network(req(res_ok())$table, req(res_ok())$input_df, input$network_pathway_pick), function() "pathway_gene_network.png")

    output$enrich_table <- DT::renderDataTable({
      r <- req(res_ok())
      cols <- intersect(c("source", "ID", "Description", "significant", "GeneRatio", "Count", "pvalue", "p.adjust", "qvalue", "NES", "ES", "geneID",
                           "transcript_gene_count", "meth_cpg_count", "integration_label", "concordance"), colnames(r$table))
      DT::datatable(r$table[, cols, drop = FALSE], rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    output$dl_enrich_csv <- downloadHandler(function() "pathway_enrichment.csv", function(file) utils::write.csv(req(res_ok())$table, file, row.names = FALSE))

    ## =========================================================================
    ## Pathway Map tab (spec 12, Plots 5-6)
    ## =========================================================================

    output$map_ui <- renderUI({
      if (!has_run()) return(box(width = NULL, title = "Pathway Map", status = "primary", solidHeader = FALSE, div(class = "empty-note", icon("circle-info"), "Not run yet.")))
      r <- res_ok()
      if (is.null(r)) return(box(width = NULL, title = "Pathway Map", status = "primary", solidHeader = FALSE, div(class = "empty-note", icon("circle-info"), "Analysis did not complete successfully - see the Enrichment tab.")))
      kegg_rows <- r$table[r$table$source == "KEGG", , drop = FALSE]
      reactome_rows <- r$table[r$table$source == "Reactome", , drop = FALSE]
      tagList(
        if (nrow(kegg_rows) > 0) box(width = NULL, title = "KEGG pathway map", status = "primary", solidHeader = FALSE,
            selectInput(ns("kegg_pick"), "Pathway", choices = stats::setNames(kegg_rows$ID, kegg_rows$Description)),
            withSpinner(imageOutput(ns("kegg_image"), height = 500), color = "#2c6fbb", type = 6)) else NULL,
        if (nrow(reactome_rows) > 0) box(width = NULL, title = "Reactome pathway visualization", status = "primary", solidHeader = FALSE,
            selectInput(ns("reactome_pick"), "Pathway", choices = stats::setNames(reactome_rows$ID, reactome_rows$Description)),
            withSpinner(imageOutput(ns("reactome_image"), height = 500), color = "#2c6fbb", type = 6)) else NULL,
        if (nrow(kegg_rows) == 0 && nrow(reactome_rows) == 0) multi_empty_state("Select KEGG and/or Reactome as a database and re-run to see pathway maps here.") else NULL
      )
    })
    output$kegg_image <- renderImage({
      r <- req(res_ok()); req(input$kegg_pick)
      genes_entrez <- stats::setNames(r$mapping_df$entrez_id, r$mapping_df$feature)
      effect <- r$input_df$expr_logFC[match(r$mapping_df$feature, r$input_df$feature)]
      names(effect) <- r$mapping_df$entrez_id
      effect <- effect[!is.na(names(effect)) & !is.na(effect)]
      out <- mp_kegg_pathway_map(input$kegg_pick, effect, out_dir = tempdir())
      validate(need(isTRUE(out$ok), out$error %||% "Could not render this KEGG pathway map."))
      list(src = out$path, contentType = "image/png", width = "100%", alt = "KEGG pathway map")
    }, deleteFile = FALSE)
    output$reactome_image <- renderImage({
      req(input$reactome_pick)
      out <- mp_fetch_reactome_diagram_png(input$reactome_pick, out_dir = tempdir())
      validate(need(isTRUE(out$ok), out$error %||% "No diagram available for this pathway."))
      list(src = out$path, contentType = "image/png", width = "100%", alt = "Reactome pathway diagram")
    }, deleteFile = FALSE)

    ## =========================================================================
    ## Gene Sets tab (spec 13)
    ## =========================================================================

    output$genesets_ui <- renderUI({
      if (!has_run()) return(box(width = NULL, title = "Gene Sets", status = "primary", solidHeader = FALSE, div(class = "empty-note", icon("circle-info"), "Not run yet.")))
      r <- res_ok()
      if (is.null(r)) return(multi_empty_state("Analysis did not complete successfully - see the Enrichment tab."))
      tagList(
        selectInput(ns("geneset_pick"), "Pathway", choices = stats::setNames(r$table$ID, r$table$Description), width = "100%"),
        div(class = "table-toolbar", downloadButton(ns("dl_geneset_csv"), "Download pathway -> genes (CSV)", class = "btn-sm")),
        DT::dataTableOutput(ns("geneset_table"))
      )
    })
    geneset_rows <- reactive({
      r <- req(res_ok()); req(input$geneset_pick)
      row <- r$table[r$table$ID == input$geneset_pick, , drop = FALSE][1, ]
      genes <- unique(trimws(unlist(strsplit(row$geneID %||% "", "/"))))
      sub <- r$input_df[toupper(r$input_df$gene_symbol) %in% toupper(genes) | toupper(r$input_df$feature) %in% toupper(genes), , drop = FALSE]
      data.frame(pathway_id = row$ID, pathway_name = row$Description, database = row$source, gene_set_size = row$BgRatio,
                 overlap_count = row$Count, gene = if (nrow(sub) > 0) sub$gene_symbol %||% sub$feature else genes,
                 omics_source = if (nrow(sub) > 0) sub$evidence_source else NA_character_,
                 expr_logFC = if (nrow(sub) > 0) sub$expr_logFC else NA_real_,
                 meth_delta_M = if (nrow(sub) > 0) sub$meth_delta_M else NA_real_, stringsAsFactors = FALSE)
    })
    output$geneset_table <- DT::renderDataTable(DT::datatable(geneset_rows(), rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact"))
    output$dl_geneset_csv <- downloadHandler(function() "pathway_genes.csv", function(file) utils::write.csv(geneset_rows(), file, row.names = FALSE))

    ## =========================================================================
    ## Results tab (spec 11, 14)
    ## =========================================================================

    output$results_ui <- renderUI({
      if (!has_run()) return(box(width = NULL, title = "Results", status = "primary", solidHeader = FALSE, div(class = "empty-note", icon("circle-info"), "Not run yet.")))
      r <- res_ok()
      if (is.null(r)) return(multi_empty_state("Analysis did not complete successfully - see the Enrichment tab."))
      card <- function(label, value) div(class = "card", style = "flex:1 1 130px; text-align:center; padding:10px;",
                                          div(style = sprintf("font-size:1.25em; font-weight:600; color:%s;", ARTHOMIX_COLORS$blue), value %||% "NA"),
                                          div(style = "font-size:0.8em; color:var(--color-ink-muted,#898781);", label))
      tagList(
        div(style = "display:flex; flex-wrap:wrap; gap:8px;",
            card("Database(s)", paste(unique(r$table$source), collapse = ", ")), card("Method", r$method),
            card("Input features", r$mapping$n_input), card("Mapped features", r$mapping$n_mapped),
            card("Mapping rate", sprintf("%s%%", r$mapping$mapping_rate)), card("Terms tested", nrow(r$table)),
            card("Significant pathways", sum(r$table$significant %||% FALSE)), card("FDR threshold", r$fdr_thresh)),
        br(),
        box(width = NULL, title = "Analysis metadata (reproducibility)", status = "primary", solidHeader = FALSE, DT::dataTableOutput(ns("meta_table"))),
        box(width = NULL, title = "Mapping results (input -> mapped identifier)", status = "primary", solidHeader = FALSE,
            div(class = "table-toolbar", downloadButton(ns("dl_mapping_csv"), "Download (CSV)", class = "btn-sm")),
            DT::dataTableOutput(ns("mapping_table"))),
        box(width = NULL, title = "Pathway genes (full table)", status = "primary", solidHeader = FALSE,
            div(class = "table-toolbar", downloadButton(ns("dl_all_genesets_csv"), "Download (CSV)", class = "btn-sm")),
            DT::dataTableOutput(ns("all_genesets_table"))),
        div(class = "table-toolbar", downloadButton(ns("dl_enrich_csv2"), "Enrichment results (CSV)", class = "btn-sm"),
            downloadButton(ns("dl_metadata_csv"), "Analysis metadata (CSV)", class = "btn-sm"))
      )
    })
    output$meta_table <- DT::renderDataTable(DT::datatable(req(res_ok())$metadata, rownames = FALSE, options = list(dom = "t", pageLength = 20), class = "stripe hover compact"))
    output$mapping_table <- DT::renderDataTable({
      r <- req(res_ok())
      cols <- intersect(c("feature", "id_type", "canonical_symbol", "entrez_id", "ensembl_id", "match_type", "mapped"), colnames(r$mapping_df))
      DT::datatable(r$mapping_df[, cols, drop = FALSE], rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    output$dl_mapping_csv <- downloadHandler(function() "pathway_mapping_results.csv", function(file) utils::write.csv(req(res_ok())$mapping_df, file, row.names = FALSE))
    all_genesets <- reactive({
      r <- req(res_ok())
      do.call(rbind, lapply(seq_len(nrow(r$table)), function(i) {
        genes <- unique(trimws(unlist(strsplit(r$table$geneID[i] %||% "", "/"))))
        if (length(genes) == 0 || all(!nzchar(genes))) return(NULL)
        data.frame(pathway_id = r$table$ID[i], pathway_name = r$table$Description[i], database = r$table$source[i], gene = genes, stringsAsFactors = FALSE)
      }))
    })
    output$all_genesets_table <- DT::renderDataTable(DT::datatable(req(all_genesets()), rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact"))
    output$dl_all_genesets_csv <- downloadHandler(function() "pathway_genes_full.csv", function(file) utils::write.csv(req(all_genesets()), file, row.names = FALSE))
    output$dl_enrich_csv2 <- downloadHandler(function() "pathway_enrichment.csv", function(file) utils::write.csv(req(res_ok())$table, file, row.names = FALSE))
    output$dl_metadata_csv <- downloadHandler(function() "pathway_analysis_metadata.csv", function(file) utils::write.csv(req(res_ok())$metadata, file, row.names = FALSE))

    ## =========================================================================
    ## Publish - keeps the multi_results$pathway$df contract other code
    ## reads (multi_qc_scorecard()/multi_analysis_summary_table()).
    ## =========================================================================

    observe({
      r <- res_ok()
      if (is.null(r) || is.null(multi_results)) return()
      multi_results$pathway <- list(df = r$table, mapping = r$mapping, metadata = r$metadata)
    })
  })
}
