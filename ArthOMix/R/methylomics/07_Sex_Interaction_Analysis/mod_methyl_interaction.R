## R/methylomics/07_Sex_Interaction_Analysis/mod_methyl_interaction.R
## Methylomics sub-module: Sex Interaction Analysis - a genuine disease*sex
## interaction limma model, the methylomics port of

mod_methyl_interaction_config <- list(
  id = "interaction", group = "Biomarker modeling",
  title = "Sex Interaction Analysis",
  description = "Diagnosis-by-sex interaction model (fit on M-values) on the currently loaded methylation data, showing which CpGs respond to the group difference differently in each sex, with the within-sex disease effects from the same fit.",
  icon = "venus-mars"
)

mod_methyl_interaction_ui <- function(id) {
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
            withSpinner(uiOutput(ns("summary_ui")), color = "#2563EB", type = 6)
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

mod_methyl_interaction_server <- function(id, methyl_dataset, methyl_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$controls <- renderUI({
      validate(need(!is.null(methyl_dataset$beta), "Load a dataset with a beta/M-value matrix on the Methylomics Dataset tab first."))
      sheet <- methyl_dataset$sample_sheet
      validate(need(!is.null(sheet), "Load a dataset with a sample sheet on the Methylomics Dataset tab first."))
      sc <- mod_methyl_dmp_sex_col(sheet)
      sexes <- if (!is.null(sc)) sort(unique(as.character(stats::na.omit(sheet[[sc]])))) else character(0)
      validate(need(length(sexes) >= 2, "Needs a sex column with at least two values to test an interaction."))
      cols <- colnames(sheet)
      tagList(
        selectInput(ns("group_col"), "Group column", choices = cols,
                    selected = intersect(c("group", "Group", "disease", "Disease"), cols)[1] %||% cols[1]),
        uiOutput(ns("group_level_ui")),
        selectInput(ns("ref_sex"), "Reference sex", choices = sexes, selected = sexes[1], selectize = FALSE),
        selectInput(ns("comp_sex"), "Comparison sex", choices = sexes, selected = sexes[min(2, length(sexes))], selectize = FALSE),
        selectizeInput(ns("covariate_cols"), "Covariates to adjust for (optional, e.g. age, smoking, chip, cell-type fractions)",
                       choices = setdiff(cols, sc), selected = NULL, multiple = TRUE,
                       options = list(placeholder = "None - add any sample-sheet column")),
        p(class = "submodule-desc", style = "font-size: 12px;",
          "Covariates enter the model additively (~ group * sex + covariates). Samples with a missing value in any chosen covariate are dropped for this run.")
      )
    })

    output$group_level_ui <- renderUI({
      req(input$group_col)
      sheet <- methyl_dataset$sample_sheet
      req(sheet, input$group_col %in% colnames(sheet))
      levels_available <- sort(unique(as.character(stats::na.omit(sheet[[input$group_col]]))))
      validate(need(length(levels_available) >= 2, "Needs at least two group values."))
      tagList(
        selectInput(ns("ref_group"), "Reference group", choices = levels_available, selected = levels_available[1], selectize = FALSE),
        selectInput(ns("comp_group"), "Comparison group", choices = levels_available, selected = levels_available[min(2, length(levels_available))], selectize = FALSE)
      )
    })

    int_has_run <- reactiveVal(FALSE)
    observeEvent(input$run_btn, int_has_run(TRUE), ignoreInit = TRUE)

    observeEvent(methyl_dataset$source, {
      int_has_run(FALSE)
    }, ignoreInit = TRUE)

    fit_result <- eventReactive(input$run_btn, {
      req(input$group_col, input$ref_group, input$comp_group, input$ref_sex, input$comp_sex)
      validate(need(input$ref_group != input$comp_group, "Reference and comparison group must be different."))
      validate(need(input$ref_sex != input$comp_sex, "Reference and comparison sex must be different."))

      validate(need(!is.null(methyl_dataset$beta), "Load a dataset first."))
      sheet <- methyl_dataset$sample_sheet
      validate(need(!is.null(sheet), "No sample sheet loaded."))
      sc <- mod_methyl_dmp_sex_col(sheet)
      validate(need(!is.null(sc), "Needs a sex column with at least two values to test an interaction."))

      sample_ids <- methyl_sheet_sample_ids(sheet, colnames(methyl_dataset$beta))
      common_all <- intersect(colnames(methyl_dataset$beta), sample_ids)
      sheet_all <- as.data.frame(sheet)[match(common_all, sample_ids), , drop = FALSE]
      beta_all <- methyl_dataset$beta[, common_all, drop = FALSE]

      covariate_cols <- setdiff(intersect(input$covariate_cols %||% character(0), colnames(sheet_all)), c(input$group_col, sc))

      keep <- !is.na(sheet_all[[input$group_col]]) & sheet_all[[input$group_col]] %in% c(input$ref_group, input$comp_group) &
        !is.na(sheet_all[[sc]]) & sheet_all[[sc]] %in% c(input$ref_sex, input$comp_sex)
      for (cc in covariate_cols) keep <- keep & !is.na(sheet_all[[cc]])
      validate(need(sum(keep) >= 12, "Fewer than 12 samples match this selection across both sexes and groups."))
      sheet1 <- sheet_all[keep, , drop = FALSE]
      beta1 <- beta_all[, keep, drop = FALSE]

      grp <- factor(as.character(sheet1[[input$group_col]]), levels = c(input$ref_group, input$comp_group))
      sx  <- factor(as.character(sheet1[[sc]]), levels = c(input$ref_sex, input$comp_sex))
      cell_table <- table(group = grp, sex = sx)
      validate(need(all(cell_table >= 2), "Each group-by-sex cell needs at least 2 samples."))

      is_m_scale <- identical(methyl_dataset$input_scale, "m")
      beta_scale <- if (is_m_scale) 2^beta1 / (1 + 2^beta1) else beta1
      m <- if (is_m_scale) {
        beta1
      } else {
        b <- pmin(pmax(beta1, 1e-6), 1 - 1e-6)
        log2(b / (1 - b))
      }

      covariate_df <- if (length(covariate_cols)) sheet1[, covariate_cols, drop = FALSE] else NULL
      built <- interaction_build_design(grp, sx, covariate_df)

      fit0 <- methyl_chunked_lmfit(m, built$design)
      fit  <- limma::eBayes(fit0)
      fit_wc <- limma::eBayes(limma::contrasts.fit(fit0, built$contrast))

      add_cpg_col <- function(tt) {
        tt$cpg <- rownames(tt)
        rownames(tt) <- NULL
        tt[, c("cpg", setdiff(colnames(tt), "cpg"))]
      }
      cn <- built$coef_names
      tt_interaction <- add_cpg_col(limma::topTable(fit, coef = cn$interaction, number = Inf, sort.by = "P"))
      tt_group       <- add_cpg_col(limma::topTable(fit, coef = cn$group, number = Inf, sort.by = "P"))
      tt_within_comp <- add_cpg_col(limma::topTable(fit_wc, coef = "within_comp_sex", number = Inf, sort.by = "P"))
      tt_sex         <- add_cpg_col(limma::topTable(fit, coef = cn$sex, number = Inf, sort.by = "P"))

      cell_ref_refsex   <- grp == levels(grp)[1] & sx == levels(sx)[1]
      cell_comp_refsex  <- grp == levels(grp)[2] & sx == levels(sx)[1]
      cell_ref_compsex  <- grp == levels(grp)[1] & sx == levels(sx)[2]
      cell_comp_compsex <- grp == levels(grp)[2] & sx == levels(sx)[2]
      mean_beta <- function(mask) rowMeans(beta_scale[, mask, drop = FALSE], na.rm = TRUE)
      dbeta_group_refsex  <- mean_beta(cell_comp_refsex) - mean_beta(cell_ref_refsex)
      dbeta_group_compsex <- mean_beta(cell_comp_compsex) - mean_beta(cell_ref_compsex)
      dbeta_interaction   <- dbeta_group_compsex - dbeta_group_refsex
      dbeta_sex           <- mean_beta(cell_ref_compsex) - mean_beta(cell_ref_refsex)
      tt_interaction$dbeta_interaction <- dbeta_interaction[tt_interaction$cpg]
      tt_group$dbeta_group <- dbeta_group_refsex[tt_group$cpg]
      tt_within_comp$dbeta_group <- dbeta_group_compsex[tt_within_comp$cpg]
      tt_sex$dbeta_sex <- dbeta_sex[tt_sex$cpg]

      list(
        tables = list(interaction = tt_interaction, group = tt_group, within_comp_sex = tt_within_comp, sex = tt_sex),
        coef_names = cn,
        coef_name = cn$interaction,
        covariates = covariate_cols,
        cell_table = cell_table,
        min_cell_n = min(cell_table),
        min_detectable_effect = interaction_min_detectable_effect(cell_table),
        ref_sex = input$ref_sex, comp_sex = input$comp_sex,
        ref_group = input$ref_group, comp_group = input$comp_group,
        is_m_scale = is_m_scale
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
      methyl_results$interaction <- list(
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
        top_hits = head(int_tbl$cpg[order(int_tbl$adj.P.Val)], 10),
        timestamp = Sys.time()
      )
    }, ignoreInit = TRUE)

    output$summary_ui <- renderUI({
      if (!int_has_run()) {
        return(div(class = "empty-note", icon("circle-info"),
                    "Not run yet. Set the group column, groups, and sexes on the left, then click \"Run interaction model\"."))
      }
      res <- fit_result()
      padj_cut <- input$padj_cut
      effect_row <- function(label, coef_name, tt) {
        n_sig <- sum(tt$adj.P.Val < padj_cut, na.rm = TRUE)
        p(strong(paste0(label, ": ")), tags$code(coef_name), " — ",
          strong(format(nrow(tt), big.mark = ",")), " CpGs tested, ",
          strong(format(n_sig, big.mark = ",")), " significant at the current cutoff.")
      }
      ct <- res$cell_table
      cell_tbl <- tags$table(class = "table table-condensed", style = "width:auto; margin-bottom:6px;",
        tags$thead(tags$tr(tags$th("Samples per cell"), lapply(colnames(ct), tags$th))),
        tags$tbody(lapply(rownames(ct), function(r) tags$tr(tags$td(strong(r)), lapply(colnames(ct), function(cl) tags$td(ct[r, cl])))))
      )
      power_note <- if (is.finite(res$min_detectable_effect)) {
        sprintf("Smallest group-by-sex cell: n = %d. Minimum detectable interaction effect at 80%% power (alpha = 0.05, per CpG, before multiple-testing adjustment) is about %.2f SD of the M-value - effects smaller than this will usually be missed.",
                res$min_cell_n, res$min_detectable_effect)
      } else NULL
      tagList(
        p(class = "submodule-desc", sprintf("Model fit on %s.", if (isTRUE(res$is_m_scale)) "M-values" else "M-values (logit-transformed from beta)")),
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
          "All four results come from the same fit. The interaction term is the primary result - it tests whether the disease effect itself differs between sexes. The two within-sex disease effects are the secondary, stratified view (a CpG can be significant in one sex and not the other without a significant interaction, which is why the interaction term is reported first). Delta-beta columns are on the beta scale for interpretability; the test statistics are on the M-value scale.")
      )
    })

    output$int_table <- DT::renderDataTable({
      if (!int_has_run()) return(NULL)
      df <- sig_table()
      num_cols <- intersect(c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "dbeta_interaction", "dbeta_group", "dbeta_sex"), colnames(df))
      DT::datatable(df, rownames = FALSE, filter = "top",
                     options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact") %>%
        DT::formatSignif(columns = num_cols, digits = 4)
    })

    output$download_int <- downloadHandler(
      filename = function() sprintf("methylation_sex_interaction_%s.csv", input$view_coef %||% "interaction"),
      content = function(file) write.csv(sig_table(), file, row.names = FALSE)
    )
  })
}
