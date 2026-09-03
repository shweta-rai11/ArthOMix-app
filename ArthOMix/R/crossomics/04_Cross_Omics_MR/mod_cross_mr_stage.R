## R/crossomics/04_Cross_Omics_MR/mod_cross_mr_stage.R
## Submodule: Cross-Omics MR - loads the pipeline's own already-run
## single-instrument mQTL-MR results (Wald ratio, GoDMC exposure -> Ishigaki

mod_cross_mr_stage_config <- list(
  id = "mrstage", title = "Cross-Omics MR", icon = "arrow-right-arrow-left", group = "Genetics",
  description = "Classifies already-computed DEG/DMP/DMR/QTL evidence into 5 named convergence categories (DEG-DMP-QTL, DEG-DMR-QTL, DEG-eQTL, DMP-mQTL, DMR-mQTL)."
)

mod_cross_mr_stage_ui <- function(id) {
  ns <- NS(id)
  category_tabs <- lapply(CX_MR_CATEGORIES, function(cat) {
    tabPanel(cat$tab, br(), uiOutput(ns(paste0("cat_", cat$id, "_ui"))))
  })
  fluidRow(
    column(
      4,
      box(
        width = NULL, title = "1. MR Data", status = "primary", solidHeader = FALSE,
        radioButtons(ns("mr_source"), "Data source",
                     choices = c("Preloaded (this project's pipeline)" = "preloaded", "Upload your own data" = "upload"),
                     selected = "preloaded"),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'preloaded'", ns("mr_source")),
          p(class = "empty-note", icon("circle-info"),
            "Precomputed mQTL-MR results (GoDMC exposure, Ishigaki 2022 RA outcome) - already run."),
          actionButton(ns("load_mr"), "Load MR Results", icon = icon("database"), class = "btn-primary btn-sm", width = "100%")
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'upload'", ns("mr_source")),
          fileInput(ns("upload_mr_file"), "MR instrument results", accept = c(".csv", ".tsv", ".txt", ".xlsx"), placeholder = "CSV / TSV / TXT / XLSX"),
          p(class = "empty-note", icon("circle-info"),
            "One row per CpG-instrument. Required: gene, pval. Optional: cpg, SNP, nsnp, b, se, OR, OR_lo, OR_hi, FDR (recomputed via Benjamini-Hochberg from pval if omitted), steiger_dir, steiger_pval."),
          actionButton(ns("load_mr_upload"), "Load Uploaded MR Results", icon = icon("upload"), class = "btn-primary btn-sm", width = "100%"),
          tags$hr(),
          p(class = "submodule-desc", icon("robot"), " Not sure which exposure/outcome dataset fits your own trait? Ask ArthOChat - it can search GWAS Catalog/OpenGWAS and PubMed for real candidate datasets (it will not fetch raw summary statistics or run MR itself; you'd still run that externally and upload the results above)."),
          textInput(ns("ai_trait_query"), NULL, placeholder = "Describe your trait/disease, e.g. \"type 2 diabetes\""),
          actionButton(ns("ai_suggest_mr"), "Ask ArthOChat for a suggested dataset", icon = icon("comments"), class = "btn-outline-primary btn-sm", width = "100%")
        ),
        uiOutput(ns("panel_info_ui"))
      ),
      box(
        width = NULL, title = "2. Evidence Source", status = "primary", solidHeader = FALSE,
        radioButtons(ns("evidence_source"), "Data source",
                     choices = c("Preloaded (this project's pipeline)" = "preloaded", "Upload your own data" = "upload"),
                     selected = "preloaded"),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'preloaded'", ns("evidence_source")),
          radioButtons(ns("sex"), "DEG/DMP/DMR/mQTL-MR/eQTL-MR evidence from", choices = c("Female" = "female", "Male" = "male", "Both (female and male)" = "combined"), selected = "female", inline = TRUE),
          p(class = "submodule-desc", "Thresholds are set on the Biomarker Convergence tab.")
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'upload'", ns("evidence_source")),
          fileInput(ns("upload_evidence_file"), "Gene-level evidence table", accept = c(".csv", ".tsv", ".txt", ".xlsx"), placeholder = "CSV / TSV / TXT / XLSX"),
          p(class = "empty-note", icon("circle-info"),
            "One row per gene. Required: gene. Optional (each defaults to \"not significant/not evaluated\" if omitted, never fabricated): DEG_adjP, DEG_logFC, DEG_direction, DMP_fdr_bacon, DMP_dbeta, DMP_direction, DMP_top_cpg, DMR_fdr, DMR_meandiff, DMR_direction, DMR_id, mQTL_MR_pval, mQTL_MR_beta, mQTL_candidate_cpg, eQTL_MR_FDR, eQTL_MR_OR, eQTL_MR_direction. Significance is relabeled at the same default thresholds (FDR < 0.05) as the Biomarker Convergence tab's own preloaded data."),
          actionButton(ns("load_evidence_upload"), "Load Uploaded Evidence", icon = icon("upload"), class = "btn-primary btn-sm", width = "100%")
        )
      )
    ),
    column(
      8,
      uiOutput(ns("summary_cards")),
      do.call(tabsetPanel, c(
        list(id = ns("result_tabs"), type = "tabs"),
        category_tabs,
        list(
          tabPanel("Results Table", br(), uiOutput(ns("results_table_ui"))),
          tabPanel("Volcanoplot", br(), uiOutput(ns("volcano_ui"))),
          tabPanel("Downloads", br(), uiOutput(ns("downloads_ui")))
        )
      ))
    )
  )
}

mod_cross_mr_stage_server <- function(id, cross_dataset, cross_results = NULL, app_session = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    mrs <- reactiveValues(df = NULL, loaded_at = NULL)

    observeEvent(input$mr_source, {
      mrs$df <- NULL; mrs$loaded_at <- NULL
    }, ignoreInit = TRUE)

    observeEvent(input$load_mr, {
      res <- cx_mr_load_precomputed()
      if (!res$ok) { showNotification(res$error, type = "error"); return() }
      mrs$df <- res$df
      mrs$loaded_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      showNotification(sprintf("Loaded %s MR instrument results.", format(nrow(res$df), big.mark = ",")), type = "message")
    })

    observeEvent(input$load_mr_upload, {
      validate(need(!is.null(input$upload_mr_file), "Upload an MR results file first."))
      res <- cx_mr_load_upload(input$upload_mr_file$datapath, input$upload_mr_file$name)
      if (!res$ok) { showNotification(res$error, type = "error", duration = 12); return() }
      mrs$df <- res$df
      mrs$loaded_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      showNotification(sprintf("Loaded %s MR instrument results from \"%s\".", format(nrow(res$df), big.mark = ","), input$upload_mr_file$name), type = "message")
    })

    output$panel_info_ui <- renderUI({
      if (identical(input$mr_source, "preloaded") && !CX_MR_DATA_AVAILABLE) return(div(class = "empty-note", icon("triangle-exclamation"), "Cross-Omics MR source data is not available in this deployment."))
      if (is.null(mrs$df)) return(NULL)
      div(class = "empty-note", icon("dna"), sprintf(
        "%s CpG-instrument results covering %s genes.", format(nrow(mrs$df), big.mark = ","), format(length(unique(mrs$df$gene)), big.mark = ",")))
    })

    observeEvent(input$ai_suggest_mr, {
      trait <- trimws(input$ai_trait_query %||% "")
      validate(need(nzchar(trait), "Describe a trait or disease first."))
      if (is.null(app_session)) { showNotification("ArthOChat is not available in this context.", type = "error"); return() }
      prompt <- sprintf(
        "Suggest Mendelian Randomization exposure and outcome datasets for this trait/disease: \"%s\". Use your gwas_catalog_search and pubmed_search tools to find real candidate GWAS Catalog / OpenGWAS accession IDs and PMIDs - the same way this project's own precomputed Cross-Omics MR data used GoDMC cis-mQTL summary statistics as the exposure and Ishigaki et al. 2022 rheumatoid arthritis GWAS (GCST90132223) as the outcome. List 2-3 concrete candidate outcome GWAS (accession ID, trait, sample size) and a plausible exposure QTL source. Note clearly that this app does not fetch raw summary statistics or run new MR itself - any dataset found would need to be run through MR externally, with the results then uploaded here via \"Upload your own data\".",
        trait
      )
      shinyjs::runjs(ARTHOCHAT_DRAWER_OPEN_JS_STATEMENT)
      shinychat::update_chat_user_input("arthochat-chat", value = prompt, submit = TRUE, session = app_session)
    })

    uploaded_evidence <- reactiveValues(df = NULL)

    observeEvent(input$load_evidence_upload, {
      validate(need(!is.null(input$upload_evidence_file), "Upload an evidence file first."))
      res <- cx_mr_load_evidence_upload(input$upload_evidence_file$datapath, input$upload_evidence_file$name)
      if (!res$ok) { showNotification(res$error, type = "error", duration = 12); return() }
      uploaded_evidence$df <- cx_bc_relabel(res$df, CX_BC_DEFAULT_PARAMS)
      showNotification(sprintf("Loaded evidence for %s genes from \"%s\".", format(nrow(res$df), big.mark = ","), input$upload_evidence_file$name), type = "message")
    })

    join_df <- reactive({
      if (identical(input$evidence_source, "upload")) return(uploaded_evidence$df)
      if (!is.null(cross_results) && !is.null(cross_results$biomarkerconv) && identical(cross_results$biomarkerconv$sex, input$sex)) {
        return(cross_results$biomarkerconv$df)
      }
      j <- cx_bc_load_precomputed(input$sex)
      if (!isTRUE(j$ok)) return(NULL)
      cx_bc_relabel(j$df)
    })

    categories <- reactive({
      req(mrs$df)
      cx_mr_classify_categories(join_df())
    })

    mrs_with_mhc <- reactive({
      df <- mrs$df; req(df)
      jd <- tryCatch(join_df(), error = function(e) NULL)
      if (!is.null(jd) && all(c("gene", "eQTL_MHC_region") %in% colnames(jd))) {
        lut <- stats::setNames(jd$eQTL_MHC_region, jd$gene)
        df$MHC_region <- unname(lut[df$gene])
      } else {
        df$MHC_region <- NA
      }
      df
    })

    mhc_warning <- function(d, mhc_col = "eQTL_MHC_region", gene_col = "gene") {
      if (is.null(d) || !mhc_col %in% colnames(d)) return(NULL)
      n_mhc <- sum(d[[mhc_col]] %in% TRUE)
      if (n_mhc == 0) return(NULL)
      n_genes <- length(unique(d[[gene_col]][d[[mhc_col]] %in% TRUE]))
      p(class = "empty-note", icon("triangle-exclamation"), style = "border-color: var(--color-warning, #e0a800);",
        sprintf(
          "%d of %d row(s) here (%d gene(s), including HLA/MHC-region genes) fall in the MHC region (chr6, ~25-34Mb) - the single most notorious horizontal-pleiotropy hotspot in autoimmune-disease genetics. Extreme regional LD means an MR estimate here is more likely to reflect the region than the gene itself; treat these rows as considerably less reliable than non-MHC hits, regardless of how significant they look.",
          n_mhc, nrow(d), n_genes
        ))
    }

    build_volcano_plot <- function(d, title, cutoff = 0.05) {
      validate(need(nrow(d) > 0, "No MR instrument data to plot here."))
      d$neglog10fdr <- -log10(pmax(d$FDR, 1e-300))
      d$sig <- d$FDR < cutoff
      d$direction <- factor(
        ifelse(!d$sig, "Not significant", ifelse(d$b > 0, "Up (risk, b > 0)", "Down (protective, b < 0)")),
        levels = CX_MR_VOLCANO_LEVELS
      )
      ggplot2::ggplot(d, ggplot2::aes(x = b, y = neglog10fdr, color = direction)) +
        ggplot2::geom_point(alpha = 0.75, size = 1.8) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = ARTHOMIX_COLORS$axis) +
        ggplot2::geom_hline(yintercept = -log10(cutoff), linetype = "dashed", color = ARTHOMIX_COLORS$axis) +
        ggrepel::geom_text_repel(data = d[d$sig %in% TRUE, , drop = FALSE], ggplot2::aes(label = gene), color = ARTHOMIX_COLORS$ink, size = 3, max.overlaps = 30) +
        ggplot2::scale_color_manual(
          values = setNames(c(ARTHOMIX_COLORS$red, ARTHOMIX_COLORS$blue, ARTHOMIX_COLORS$ink_muted), CX_MR_VOLCANO_LEVELS),
          drop = FALSE, name = NULL
        ) +
        ggplot2::labs(x = "MR estimate (b, log-odds scale)", y = "-log10(FDR)", title = title) +
        theme_arthomix()
    }

    output$summary_cards <- renderUI({
      if (is.null(mrs$df)) return(NULL)
      df <- mrs$df; cats <- categories()
      card <- function(label, value, color = "blue") {
        div(class = "card", style = "flex: 1 1 150px; text-align:center; padding:10px;",
            div(style = sprintf("font-size:1.4em; font-weight:600; color:%s;", ARTHOMIX_COLORS[[color]] %||% ARTHOMIX_COLORS$ink), format(value, big.mark = ",")),
            div(style = "font-size:0.82em; color:var(--color-ink-muted, #898781);", label))
      }
      cat_cards <- if (!is.null(cats)) lapply(CX_MR_CATEGORIES, function(cat) card(cat$tab, nrow(cats[[cat$id]]), "violet")) else list()
      div(style = "display:flex; gap:10px; flex-wrap:wrap; margin-bottom:12px;",
          c(list(card("Instruments", nrow(df)), card("Genes", length(unique(df$gene)), "aqua")), cat_cards)
      )
    })

    lapply(CX_MR_CATEGORIES, function(cat) {
      local({
        this_cat <- cat
        output[[paste0("cat_", this_cat$id, "_ui")]] <- renderUI({
          if (is.null(mrs$df)) return(cx_empty_state("Load MR results on the \"1. MR Data\" panel to see results here."))
          cats <- categories()
          if (is.null(cats)) return(div(class = "empty-note", icon("triangle-exclamation"), sprintf("DEG/DMP/DMR/mQTL-MR/eQTL-MR evidence was not available (%s).", if (identical(input$evidence_source, "upload")) "upload an evidence file above" else sprintf("%s not available", toupper(input$sex)))))
          d <- cats[[this_cat$id]]
          tagList(
            p(class = "submodule-desc", this_cat$rule_text),
            p(class = "submodule-desc", sprintf("%s gene(s) match.", format(nrow(d), big.mark = ","))),
            mhc_warning(d),
            div(class = "table-toolbar", downloadButton(ns(paste0("dl_cat_", this_cat$id, "_csv")), "CSV", class = "btn-sm")),
            DT::dataTableOutput(ns(paste0("cat_", this_cat$id, "_table")))
          )
        })
        output[[paste0("cat_", this_cat$id, "_table")]] <- DT::renderDataTable({
          cats <- categories()
          req(cats)
          DT::datatable(cats[[this_cat$id]], rownames = FALSE, filter = "top", options = list(scrollX = TRUE, pageLength = 15), class = "stripe hover compact")
        })
        output[[paste0("dl_cat_", this_cat$id, "_csv")]] <- downloadHandler(
          filename = function() sprintf("cross_omics_mr_%s_%s.csv", this_cat$id, Sys.Date()),
          content = function(file) {
            cats <- categories()
            utils::write.csv(if (!is.null(cats)) cats[[this_cat$id]] else data.frame(), file, row.names = FALSE)
          }
        )
      })
    })

    output$results_table_ui <- renderUI({
      if (is.null(mrs$df)) return(cx_empty_state("Load MR results on the \"1. MR Data\" panel to see results here."))
      tagList(
        p(class = "submodule-desc", "Mendelian randomization estimates - valid under the standard instrumental-variable assumptions. A single-instrument gene cannot be tested for validity via heterogeneity; a gene with several independent CpG-instruments (shown as separate rows here) could be, but this module does not currently aggregate per gene or run that test - treat every row as an association consistent with a causal effect, not as proof of one. No instrument-strength (F-statistic) figure is available from the precomputed file."),
        mhc_warning(mrs_with_mhc(), mhc_col = "MHC_region"),
        if (any(mrs$df$steiger_dir %in% FALSE)) p(class = "empty-note", icon("triangle-exclamation"),
          sprintf("%d row(s) fail the Steiger directionality test (steiger_dir = FALSE) - possible reverse causation. Filter the \"steiger_dir\" column below to inspect them.", sum(mrs$df$steiger_dir %in% FALSE))),
        div(class = "table-toolbar", downloadButton(ns("dl_table_csv"), "CSV", class = "btn-sm"), downloadButton(ns("dl_table_xlsx"), "XLSX", class = "btn-sm")),
        DT::dataTableOutput(ns("results_table"))
      )
    })
    output$results_table <- DT::renderDataTable({
      d <- mrs_with_mhc()
      cols <- intersect(c("gene", "cpg", "SNP", "nsnp", "b", "se", "pval", "OR", "OR_lo", "OR_hi", "FDR", "mr_significant", "steiger_dir", "steiger_pval", "MHC_region"), colnames(d))
      dt <- DT::datatable(d[, cols, drop = FALSE], rownames = FALSE, filter = "top", options = list(scrollX = TRUE, pageLength = 15), class = "stripe hover compact")
      if ("steiger_dir" %in% cols) {
        dt <- DT::formatStyle(dt, "steiger_dir", target = "row",
                               backgroundColor = DT::styleEqual(c(FALSE), c("#fbe4e0")))
      }
      dt
    })

    volcano_cat <- reactiveVal(NULL)
    lapply(CX_MR_CATEGORIES, function(cat) {
      local({
        this_id <- cat$id
        observeEvent(input[[paste0("volc_btn_", this_id)]], {
          volcano_cat(if (identical(volcano_cat(), this_id)) NULL else this_id)
        })
      })
    })

    output$volcano_ui <- renderUI({
      if (is.null(mrs$df)) return(cx_empty_state("Load MR results on the \"1. MR Data\" panel to see results here."))
      sel <- volcano_cat()
      tagList(
        div(class = "table-toolbar",
            lapply(CX_MR_CATEGORIES, function(cat) actionButton(
              ns(paste0("volc_btn_", cat$id)), cat$tab,
              class = if (identical(sel, cat$id)) "btn-sm btn-primary" else "btn-sm btn-outline-primary"
            )),
            downloadButton(ns("dl_volcano_png"), "PNG", class = "btn-sm")),
        p(class = "submodule-desc",
          if (is.null(sel)) "Showing every loaded instrument. Click a category above to restrict the plot to just its matching genes; click it again to go back to all instruments."
          else sprintf("Showing only %s's matching genes. Click it again to go back to all instruments.", Find(function(c) identical(c$id, sel), CX_MR_CATEGORIES)$tab)),
        plotOutput(ns("volcano_plot"), height = "460px")
      )
    })
    current_volcano_plot <- reactive({
      df <- mrs$df
      validate(need(!is.null(df) && nrow(df) > 0, "Load MR results first."))
      d <- df[!is.na(df$b) & !is.na(df$FDR), , drop = FALSE]
      sel <- volcano_cat()
      if (is.null(sel)) {
        build_volcano_plot(d, "Cross-Omics MR (Wald ratio) volcano - every loaded instrument")
      } else {
        cats <- categories()
        validate(need(!is.null(cats), "DEG/DMP/DMR/mQTL-MR/eQTL-MR evidence is not available - load it on the Evidence Source panel first."))
        cat_tab <- Find(function(c) identical(c$id, sel), CX_MR_CATEGORIES)$tab
        genes <- unique(cats[[sel]]$gene)
        d <- d[d$gene %in% genes, , drop = FALSE]
        build_volcano_plot(d, sprintf("%s volcano (%s matching gene(s) with an MR instrument)", cat_tab, format(length(unique(d$gene)), big.mark = ",")))
      }
    })
    output$volcano_plot <- renderPlot({ current_volcano_plot() })
    output$dl_volcano_png <- downloadHandler(
      filename = function() sprintf("cross_omics_mr_volcano_%s_%s.png", volcano_cat() %||% "all", Sys.Date()),
      content = function(file) ggplot2::ggsave(file, plot = current_volcano_plot(), width = 9, height = 6.5, dpi = 200, bg = "white")
    )

    output$downloads_ui <- renderUI({
      if (is.null(mrs$df)) return(cx_empty_state("Load MR results on the \"1. MR Data\" panel to see results here."))
      tagList(
        downloadButton(ns("dl_table_csv"), "MR results (CSV)", class = "btn-sm"),
        downloadButton(ns("dl_table_xlsx"), "MR results (XLSX)", class = "btn-sm"),
        tags$hr(),
        p(class = "submodule-desc", "Evidence-combination category tables (same 5 tables shown in their own tabs above):"),
        tagList(lapply(CX_MR_CATEGORIES, function(cat) downloadButton(ns(paste0("dl_cat_", cat$id, "_csv")), sprintf("%s (CSV)", cat$tab), class = "btn-sm")))
      )
    })
    output$dl_table_csv <- downloadHandler(function() paste0("cross_omics_mr_", Sys.Date(), ".csv"), function(file) utils::write.csv(mrs$df, file, row.names = FALSE))
    output$dl_table_xlsx <- downloadHandler(function() paste0("cross_omics_mr_", Sys.Date(), ".xlsx"), function(file) {
      if (!requireNamespace("openxlsx", quietly = TRUE)) { utils::write.csv(mrs$df, file, row.names = FALSE); return() }
      openxlsx::write.xlsx(mrs$df, file)
    })

    observe({
      if (is.null(mrs$df) || is.null(cross_results)) return()
      cross_results$mrstage <- list(df = mrs$df, categories = categories(), run_at = mrs$loaded_at)
    })
  })
}
