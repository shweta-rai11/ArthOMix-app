## Sex Interaction Analysis submodule (Section 2.10): live limma
## group*sex interaction model on the currently loaded dataset.

mod_interaction_config <- list(
  id = "interaction", group = "Biomarker modeling",
  title = "Sex Interaction Analysis",
  description = "Diagnosis-by-sex interaction model on the preloaded or uploaded data, showing which genes respond to the group difference differently in each sex.",
  icon = "venus-mars"
)

mod_interaction_ui <- function(id) {
  ns <- NS(id)
  tagList(
      fluidRow(
        column(
          4,
          box(
            width = NULL, title = "Model", status = "primary", solidHeader = FALSE,
            uiOutput(ns("controls")),
            numericInput(ns("padj_cut"), "Adjusted p-value cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
            actionButton(ns("run_btn"), "Run interaction model", icon = icon("play"), class = "btn-primary btn-sm")
          )
        ),
        column(
          8,
          box(
            width = NULL, title = "Result", status = "primary", solidHeader = FALSE,
            withSpinner(uiOutput(ns("summary_ui")), color = "#2c6fbb", type = 6)
          )
        )
      ),
      box(
        width = 12, title = "Interaction result table", status = "primary", solidHeader = FALSE,
        div(class = "table-toolbar", downloadButton(ns("download_int"), "Download CSV", class = "btn-sm")),
        DT::dataTableOutput(ns("int_table"))
      )
  )
}

mod_interaction_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$controls <- renderUI({
      groups <- sort(unique(na.omit(dataset$meta$group)))
      sexes  <- sort(unique(na.omit(dataset$meta$sex)))
      validate(need(length(groups) >= 2, "Needs at least two group values."))
      validate(need(length(sexes) >= 2, "Needs a sex column with at least two values to test an interaction."))
      tagList(
        selectInput(ns("ref_group"), "Reference group", choices = groups, selected = groups[1], selectize = FALSE),
        selectInput(ns("comp_group"), "Comparison group", choices = groups, selected = groups[min(2, length(groups))], selectize = FALSE),
        selectInput(ns("ref_sex"), "Reference sex", choices = sexes, selected = sexes[1], selectize = FALSE),
        selectInput(ns("comp_sex"), "Comparison sex", choices = sexes, selected = sexes[min(2, length(sexes))], selectize = FALSE)
      )
    })

    int_has_run <- reactiveVal(FALSE)
    observeEvent(input$run_btn, int_has_run(TRUE), ignoreInit = TRUE)

    ## fit_result() is an eventReactive, so it keeps the previous dataset's fit
    ## until Run is clicked again; clear the gate so nothing stale stays on screen.
    observeEvent(dataset$source, {
      int_has_run(FALSE)
    }, ignoreInit = TRUE)

    fit_result <- eventReactive(input$run_btn, {
      req(input$ref_group, input$comp_group, input$ref_sex, input$comp_sex)
      validate(need(input$ref_group != input$comp_group, "Reference and comparison group must be different."))
      validate(need(input$ref_sex != input$comp_sex, "Reference and comparison sex must be different."))

      meta <- dataset$meta
      meta <- meta[meta$group %in% c(input$ref_group, input$comp_group) &
                     !is.na(meta$sex) & meta$sex %in% c(input$ref_sex, input$comp_sex), , drop = FALSE]
      common <- intersect(colnames(dataset$expr), meta$sample)
      validate(need(length(common) >= 12, "Fewer than 12 samples match this selection across both sexes and groups."))
      meta <- meta[match(common, meta$sample), , drop = FALSE]
      expr <- dataset$expr[, common, drop = FALSE]

      grp <- factor(meta$group, levels = c(input$ref_group, input$comp_group))
      sx  <- factor(meta$sex, levels = c(input$ref_sex, input$comp_sex))
      validate(need(all(table(grp, sx) >= 2), "Each group-by-sex cell needs at least 2 samples."))

      design <- model.matrix(~ grp * sx)
      fit <- limma::eBayes(limma::lmFit(expr, design))
      coef_name <- tail(colnames(design), 1)
      tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
      tt$gene <- rownames(tt)
      rownames(tt) <- NULL
      tt <- tt[, c("gene", setdiff(colnames(tt), "gene"))]
      list(table = tt, coef_name = coef_name)
    }, ignoreInit = TRUE)

    sig_table <- reactive({
      res <- fit_result()
      req(res)
      res$table %>% mutate(significant = adj.P.Val < input$padj_cut)
    })

    ## Publishes the fitted interaction model into shared results$interaction so
    ## ArthOChat and the biomarker card can see it, mirroring mod_dge.R's
    ## results$dge contract (contrast label, n tested/significant, top hits) -
    ## silently skipped if the fit failed validation, same as mod_dge.R.
    observeEvent(input$run_btn, {
      df <- tryCatch(sig_table(), error = function(e) NULL)
      req(df)
      res <- fit_result()
      results$interaction <- list(
        contrast = sprintf("%s vs %s interaction with sex (%s vs %s)",
                            input$comp_group, input$ref_group, input$comp_sex, input$ref_sex),
        coef_name = res$coef_name,
        n_tested = nrow(df),
        n_significant = sum(df$significant),
        top_hits = head(df$gene[order(df$adj.P.Val)], 10),
        timestamp = Sys.time()
      )
    }, ignoreInit = TRUE)

    output$summary_ui <- renderUI({
      if (!int_has_run()) {
        return(div(class = "empty-note", icon("circle-info"),
                    "Not run yet. Set the groups and sexes on the left, then click \"Run interaction model\"."))
      }
      res <- fit_result()
      df <- sig_table()
      tagList(
        p(strong("Interaction term: "), res$coef_name),
        p(strong(format(nrow(df), big.mark = ",")), " genes tested, ",
          strong(format(sum(df$significant), big.mark = ",")), " with a significant group-by-sex interaction at the current cutoff.")
      )
    })

    output$int_table <- DT::renderDataTable({
      ## req(), not `if (!int_has_run()) return(NULL)`: an explicit NULL still
      ## gets sent to the client as a real value, and DT's own JS binding
      ## throws "Cannot read properties of null (reading 'lazyRender')" the
      ## first time it's asked to render one - which happens as soon as this
      ## tab is added (insertTab(..., select = TRUE) in server.R makes it
      ## visible immediately, before Run is ever clicked). Same fix as
      ## mod_dge.R's dge_table.
      req(int_has_run())
      DT::datatable(sig_table(), rownames = FALSE, filter = "top",
                     options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    outputOptions(output, "int_table", suspendWhenHidden = FALSE)

    output$download_int <- downloadHandler(
      filename = function() "sex_interaction.csv",
      content = function(file) {
        if (!int_has_run()) stop("No interaction model run yet in this session - click \"Run interaction model\" first.")
        write.csv(sig_table(), file, row.names = FALSE)
      }
    )
  })
}
