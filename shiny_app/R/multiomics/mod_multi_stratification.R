## R/multiomics/mod_multi_stratification.R
## Submodule: Patient Stratification (SNF Clusters) - the pipeline's own
## unsupervised SNF joint-biomarker-discovery workflow (script 07c, Wang et
## al. 2014 method) fused one similarity network per omics layer across the
## FULL cohort (not response-supervised) and spectral-clustered the result;
## this tab shows those real per-patient cluster assignments and the
## pipeline's own Fisher-test association with clinical response - it does
## not re-run SNF or re-cluster.

mod_multi_stratification_config <- list(
  id = "stratification", title = "Patient Stratification (SNF Clusters)", icon = "diagram-project", group = "Data",
  description = "Real per-patient clusters from unsupervised similarity-network fusion across expression + methylation, and whether cluster membership associates with clinical response."
)

MULTI_STRAT_DRUG_CHOICES <- c("Adalimumab" = "Adalimumab", "Etanercept" = "Etanercept")

mod_multi_stratification_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      4,
      box(
        width = NULL, title = "1. Drug cohort", status = "primary", solidHeader = FALSE,
        selectInput(ns("drug"), "Drug", choices = MULTI_STRAT_DRUG_CHOICES, width = "100%"),
        p(class = "submodule-desc", "SNF joint clustering was run on the full drug cohort (both sexes together), unsupervised - cluster labels are not derived from the response outcome."),
        actionButton(ns("load_btn"), "Load clusters", icon = icon("play"), class = "btn-primary btn-sm", width = "100%")
      )
    ),
    column(
      8,
      uiOutput(ns("assoc_banner")),
      conditionalPanel(
        condition = sprintf("input['%s'] > 0", ns("load_btn")),
        tabsetPanel(
          id = ns("tabs"), type = "tabs",
          tabPanel("Cluster composition", br(), uiOutput(ns("composition_ui"))),
          tabPanel("Omics contribution (NMI)", br(), uiOutput(ns("nmi_ui"))),
          tabPanel("Patient assignments", br(), uiOutput(ns("assignments_ui")))
        )
      )
    )
  )
}

mod_multi_stratification_server <- function(id, multi_dataset = NULL, multi_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    dat <- reactiveValues(assign = NULL, assoc = NULL, drug = NULL, nmi = NULL)

    observeEvent(input$load_btn, {
      assign_label <- sprintf("SNF patient clusters — %s", input$drug)
      assoc_label  <- sprintf("SNF cluster-response association — %s", input$drug)
      nmi_label    <- sprintf("SNF concordance NMI — %s", input$drug)
      a <- multi_read_registry_table(assign_label)
      s <- multi_read_registry_table(assoc_label)
      n <- multi_read_registry_table(nmi_label)
      if (!a$ok) { showNotification(a$error, type = "error"); return() }
      dat$assign <- a$df
      dat$assoc <- if (s$ok) s$df else NULL
      dat$nmi <- if (n$ok) n$df else NULL
      dat$drug <- input$drug
      showNotification(sprintf("Loaded SNF clusters for %s (%d patients).", input$drug, nrow(a$df)), type = "message")
    })

    output$assoc_banner <- renderUI({
      if (is.null(dat$assoc) || nrow(dat$assoc) == 0) return(if (is.null(dat$assign)) multi_empty_state("Pick a drug cohort and click \"Load clusters\".") else NULL)
      rows <- dat$assoc
      tagList(lapply(seq_len(nrow(rows)), function(i) {
        r <- rows[i, ]
        sig <- !is.na(r$fisher_p_cluster_vs_response) && r$fisher_p_cluster_vs_response < 0.05
        div(class = "empty-note", style = if (sig) "border-color: var(--color-warning, #eda100);" else NULL,
            icon(if (sig) "circle-check" else "circle-info"),
            sprintf(" %s, %s: %d fused clusters, n=%d, Fisher p (cluster vs. response) = %.3f%s",
                    r$sex, dat$drug %||% "", r$n_fused_clusters, r$n, r$fisher_p_cluster_vs_response,
                    if (sig) " (nominally significant)" else " (not significant at p<0.05)"))
      }))
    })

    output$composition_ui <- renderUI({
      if (is.null(dat$assign)) return(multi_empty_state())
      tagList(
        selectInput(ns("sex_filter"), "Sex", choices = c("Both" = "both", "Female" = "female", "Male" = "male"), selected = "both"),
        multi_plot_or_empty(comp_plot_fn, ns("composition_plot"), height = "360px"),
        div(class = "table-toolbar", downloadButton(ns("dl_comp_png"), "Download plot (PNG)", class = "btn-sm"))
      )
    })
    filtered_assign <- reactive({
      req(dat$assign)
      if (identical(input$sex_filter %||% "both", "both")) return(dat$assign)
      multi_filter_cell(dat$assign, sex = input$sex_filter)
    })
    comp_plot_fn <- reactive(multi_cluster_composition_plot(filtered_assign()))
    output$composition_plot <- renderPlot(comp_plot_fn())
    output$dl_comp_png <- multi_png_download(comp_plot_fn, function() sprintf("multiomics_snf_clusters_%s.png", dat$drug %||% "cohort"))

    output$nmi_ui <- renderUI({
      if (is.null(dat$nmi)) return(multi_empty_state("No concordance-NMI table loaded for this drug cohort."))
      tagList(
        p(class = "submodule-desc", "Normalized Mutual Information between each per-omics similarity network and the fused network - the real, already-computed measure of how much each omics layer agrees with the final clustering (Wang et al. 2014 SNF method), not a MOFA-style invented number."),
        multi_plot_or_empty(nmi_plot_fn, ns("nmi_plot"), height = "320px"),
        div(class = "table-toolbar", downloadButton(ns("dl_nmi_png"), "Download plot (PNG)", class = "btn-sm"),
            downloadButton(ns("dl_nmi_csv"), "Download data (CSV)", class = "btn-sm"))
      )
    })
    nmi_plot_fn <- reactive(multi_nmi_heatmap(dat$nmi))
    output$nmi_plot <- renderPlot(nmi_plot_fn())
    output$dl_nmi_png <- multi_png_download(nmi_plot_fn, function() sprintf("multiomics_snf_nmi_%s.png", dat$drug %||% "cohort"))
    output$dl_nmi_csv <- downloadHandler(function() sprintf("multiomics_snf_nmi_%s.csv", dat$drug %||% "cohort"), function(file) utils::write.csv(dat$nmi, file, row.names = FALSE))

    output$assignments_ui <- renderUI({
      if (is.null(dat$assign)) return(multi_empty_state())
      tagList(
        div(class = "table-toolbar", downloadButton(ns("dl_assign_csv"), "Download assignments (CSV)", class = "btn-sm")),
        DT::dataTableOutput(ns("assign_table"))
      )
    })
    output$assign_table <- DT::renderDataTable({
      req(dat$assign)
      DT::datatable(dat$assign, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    output$dl_assign_csv <- downloadHandler(function() sprintf("multiomics_snf_cluster_assignments_%s.csv", dat$drug %||% "cohort"), function(file) utils::write.csv(dat$assign, file, row.names = FALSE))

    observe({
      if (is.null(dat$assign) || is.null(multi_results)) return()
      multi_results$stratification <- list(drug = dat$drug, assignments = dat$assign, association = dat$assoc, nmi = dat$nmi)
    })
  })
}
