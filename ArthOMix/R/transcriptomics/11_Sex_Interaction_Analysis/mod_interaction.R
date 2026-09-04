## Sex Interaction Analysis submodule (Section 2.10): live limma
## group*sex interaction model on the currently loaded dataset.

mod_interaction_config <- list(
  id = "interaction", group = "Biomarker modeling",
  title = "Sex Interaction Analysis",
  description = "Diagnosis-by-sex interaction model on the preloaded or uploaded data, showing which genes respond to the group difference differently in each sex, with the within-sex disease effects from the same fit.",
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
        width = 12, title = "Result table", status = "primary", solidHeader = FALSE,
        selectInput(ns("view_coef"), "Show coefficient", choices = INTERACTION_VIEW_CHOICES,
                    selected = "interaction", selectize = FALSE),
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
      covar_choices <- setdiff(colnames(dataset$meta), c("sample", "group", "sex"))
      tagList(
        selectInput(ns("ref_group"), "Reference group", choices = groups, selected = groups[1], selectize = FALSE),
        selectInput(ns("comp_group"), "Comparison group", choices = groups, selected = groups[min(2, length(groups))], selectize = FALSE),
        selectInput(ns("ref_sex"), "Reference sex", choices = sexes, selected = sexes[1], selectize = FALSE),
        selectInput(ns("comp_sex"), "Comparison sex", choices = sexes, selected = sexes[min(2, length(sexes))], selectize = FALSE),
        selectizeInput(ns("covariate_cols"), "Covariates to adjust for (optional, e.g. age, batch, cell-type fractions)",
                       choices = covar_choices, selected = NULL, multiple = TRUE,
                       options = list(placeholder = "None - add any metadata column")),
        p(class = "submodule-desc", style = "font-size: 12px;",
          "Covariates enter the model additively (~ group * sex + covariates). Samples with a missing value in any chosen covariate are dropped for this run.")
      )
    })

    int_has_run <- reactiveVal(FALSE)
    observeEvent(input$run_btn, int_has_run(TRUE), ignoreInit = TRUE)

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

      covariate_cols <- setdiff(intersect(input$covariate_cols %||% character(0), colnames(meta)), c("sample", "group", "sex"))
      for (cc in covariate_cols) meta <- meta[!is.na(meta[[cc]]), , drop = FALSE]

      common <- intersect(colnames(dataset$expr), meta$sample)
      validate(need(length(common) >= 12, "Fewer than 12 samples match this selection across both sexes and groups."))
      meta <- meta[match(common, meta$sample), , drop = FALSE]
      expr <- dataset$expr[, common, drop = FALSE]

      grp <- factor(meta$group, levels = c(input$ref_group, input$comp_group))
      sx  <- factor(meta$sex, levels = c(input$ref_sex, input$comp_sex))
      cell_table <- table(group = grp, sex = sx)
      validate(need(all(cell_table >= 2), "Each group-by-sex cell needs at least 2 samples."))

      covariate_df <- if (length(covariate_cols)) meta[, covariate_cols, drop = FALSE] else NULL
      built <- interaction_build_design(grp, sx, covariate_df)

      fit0 <- limma::lmFit(expr, built$design)
      fit  <- limma::eBayes(fit0)
      fit_wc <- limma::eBayes(limma::contrasts.fit(fit0, built$contrast))

      add_gene_col <- function(tt) {
        tt$gene <- rownames(tt)
        rownames(tt) <- NULL
        tt[, c("gene", setdiff(colnames(tt), "gene"))]
      }
      cn <- built$coef_names
      list(
        tables = list(
          interaction     = add_gene_col(limma::topTable(fit, coef = cn$interaction, number = Inf, sort.by = "P")),
          group           = add_gene_col(limma::topTable(fit, coef = cn$group, number = Inf, sort.by = "P")),
          within_comp_sex = add_gene_col(limma::topTable(fit_wc, coef = "within_comp_sex", number = Inf, sort.by = "P")),
          sex             = add_gene_col(limma::topTable(fit, coef = cn$sex, number = Inf, sort.by = "P"))
        ),
        coef_names = cn,
        coef_name = cn$interaction,
        covariates = covariate_cols,
        cell_table = cell_table,
        min_cell_n = min(cell_table),
        min_detectable_effect = interaction_min_detectable_effect(cell_table),
        ref_sex = input$ref_sex, comp_sex = input$comp_sex,
        ref_group = input$ref_group, comp_group = input$comp_group
      )
    }, ignoreInit = TRUE)

    sig_table <- reactive({
      res <- fit_result()
      req(res)
      view <- input$view_coef %||% "interaction"
      tbl <- res$tables[[view]]
      req(tbl)
      tbl %>% mutate(significant = adj.P.Val < input$padj_cut)
    })

    observeEvent(input$run_btn, {
      res <- tryCatch(fit_result(), error = function(e) NULL)
      req(res)
      padj_cut <- input$padj_cut
      sig_counts <- lapply(res$tables, function(tt) sum(tt$adj.P.Val < padj_cut, na.rm = TRUE))
      int_tbl <- res$tables$interaction
      results$interaction <- list(
        contrast = sprintf("%s vs %s interaction with sex (%s vs %s)",
                            input$comp_group, input$ref_group, input$comp_sex, input$ref_sex),
        coef_name = res$coef_names$interaction,
        covariates = res$covariates,
        n_tested = nrow(int_tbl),
        n_significant = sig_counts$interaction,
        n_significant_within_ref_sex = sig_counts$group,
        n_significant_within_comp_sex = sig_counts$within_comp_sex,
        n_significant_main_sex = sig_counts$sex,
        min_cell_n = res$min_cell_n,
        min_detectable_effect = res$min_detectable_effect,
        top_hits = head(int_tbl$gene[order(int_tbl$adj.P.Val)], 10),
        timestamp = Sys.time()
      )
    }, ignoreInit = TRUE)

    output$summary_ui <- renderUI({
      if (!int_has_run()) {
        return(div(class = "empty-note", icon("circle-info"),
                    "Not run yet. Set the groups and sexes on the left, then click \"Run interaction model\"."))
      }
      res <- fit_result()
      padj_cut <- input$padj_cut
      effect_row <- function(label, coef_name, tt) {
        n_sig <- sum(tt$adj.P.Val < padj_cut, na.rm = TRUE)
        p(strong(paste0(label, ": ")), tags$code(coef_name), " — ",
          strong(format(nrow(tt), big.mark = ",")), " genes tested, ",
          strong(format(n_sig, big.mark = ",")), " significant at the current cutoff.")
      }
      ct <- res$cell_table
      cell_tbl <- tags$table(class = "table table-condensed", style = "width:auto; margin-bottom:6px;",
        tags$thead(tags$tr(tags$th("Samples per cell"), lapply(colnames(ct), tags$th))),
        tags$tbody(lapply(rownames(ct), function(r) tags$tr(tags$td(strong(r)), lapply(colnames(ct), function(cl) tags$td(ct[r, cl])))))
      )
      power_note <- if (is.finite(res$min_detectable_effect)) {
        sprintf("Smallest group-by-sex cell: n = %d. Minimum detectable interaction effect at 80%% power (alpha = 0.05, per gene, before multiple-testing adjustment) is about %.2f SD of expression - effects smaller than this will usually be missed.",
                res$min_cell_n, res$min_detectable_effect)
      } else NULL
      tagList(
        cell_tbl,
        if (res$min_cell_n < 10)
          div(class = "empty-note", style = "border-color: var(--color-warning, #f0ad4e);", icon("triangle-exclamation"),
              " Underpowered interaction test: at least one group-by-sex cell has fewer than 10 samples. ", power_note)
        else if (!is.null(power_note)) p(class = "submodule-desc", icon("circle-info"), " ", power_note),
        effect_row("Group x sex interaction (primary)", res$coef_names$interaction, res$tables$interaction),
        effect_row(sprintf("Disease effect within %s (reference sex)", res$ref_sex), res$coef_names$group, res$tables$group),
        effect_row(sprintf("Disease effect within %s (comparison sex)", res$comp_sex), res$coef_names$within_comp_sex, res$tables$within_comp_sex),
        effect_row(sprintf("Main sex effect within %s (reference group)", res$ref_group), res$coef_names$sex, res$tables$sex),
        if (length(res$covariates)) p(class = "empty-note", icon("circle-info"), sprintf("Adjusted for: %s", paste(res$covariates, collapse = ", "))) else NULL,
        p(class = "empty-note", icon("circle-info"),
          "All four results come from the same fit. The interaction term is the primary result - it tests whether the disease effect itself differs between sexes. The two within-sex disease effects are the secondary, stratified view (a gene can be significant in one sex and not the other without a significant interaction, which is why the interaction term is reported first).")
      )
    })

    output$int_table <- DT::renderDataTable({
      req(int_has_run())
      DT::datatable(sig_table(), rownames = FALSE, filter = "top",
                     options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    outputOptions(output, "int_table", suspendWhenHidden = FALSE)

    output$download_int <- downloadHandler(
      filename = function() sprintf("sex_interaction_%s.csv", input$view_coef %||% "interaction"),
      content = function(file) {
        if (!int_has_run()) stop("No interaction model run yet in this session - click \"Run interaction model\" first.")
        write.csv(sig_table(), file, row.names = FALSE)
      }
    )
  })
}
