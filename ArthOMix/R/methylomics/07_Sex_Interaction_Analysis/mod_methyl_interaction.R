## R/methylomics/07_Sex_Interaction_Analysis/mod_methyl_interaction.R
## Methylomics sub-module: Sex Interaction Analysis - a genuine disease*sex
## interaction limma model, the methylomics port of
## R/transcriptomics/11_Sex_Interaction_Analysis/mod_interaction.R (see that file's header comment for
## why this is needed: a "sex-stratified" analysis that just re-runs a
## disease-only model separately per sex never actually tests whether the
## disease effect differs by sex - only a group*sex interaction term does).
##
## Methylomics-specific differences from the transcriptomics version:
##  - no fixed "group"/"sex" column names - the group column is picked from
##    the loaded sample sheet (same input$live_group_col-style pattern as
##    mod_methyl_dmp.R/mod_methyl_dmr.R), and the sex column is
##    auto-detected via mod_methyl_dmp_sex_col() (mod_methyl_dmp.R).
##  - sample IDs are reconciled between the beta matrix and the sample sheet
##    via methyl_sheet_sample_ids() (R/methylomics/functions/qc.R), not a fixed
##    "sample" column.
##  - the model is fit on M-values (logit-transformed beta, or the matrix
##    as-is when already M-scale), same convention as mod_methyl_dmp.R's
##    live engine, via the shared memory-safe methyl_chunked_lmfit()
##    wrapper (mod_methyl_dmp.R) instead of calling limma::lmFit()
##    directly - a full genome-wide array can be 400k+ probes. A
##    beta-scale delta-beta for the interaction cells is reported alongside
##    the M-value-scale logFC for interpretability, mirroring how
##    mod_methyl_dmp.R reports both scales.

mod_methyl_interaction_config <- list(
  id = "interaction", group = "Biomarker modeling",
  title = "Sex Interaction Analysis",
  description = "Diagnosis-by-sex interaction model (fit on M-values) on the currently loaded methylation data, showing which CpGs respond to the group difference differently in each sex.",
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
        width = 12, title = "Interaction result table", status = "primary", solidHeader = FALSE,
        div(class = "table-toolbar", downloadButton(ns("download_int"), "Download CSV", class = "btn-sm")),
        DT::dataTableOutput(ns("int_table"))
      )
  )
}

mod_methyl_interaction_server <- function(id, methyl_dataset, methyl_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## ---- Stage 1: sex-column detection + group-column picker ---------------

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
        selectInput(ns("comp_sex"), "Comparison sex", choices = sexes, selected = sexes[min(2, length(sexes))], selectize = FALSE)
      )
    })

    ## ---- Stage 2: group levels, dependent on the chosen group column --------

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

    ## fit_result() is an eventReactive, so it keeps the previous dataset's fit
    ## until Run is clicked again; clear the gate so nothing stale stays on screen.
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

      keep <- !is.na(sheet_all[[input$group_col]]) & sheet_all[[input$group_col]] %in% c(input$ref_group, input$comp_group) &
        !is.na(sheet_all[[sc]]) & sheet_all[[sc]] %in% c(input$ref_sex, input$comp_sex)
      validate(need(sum(keep) >= 12, "Fewer than 12 samples match this selection across both sexes and groups."))
      sheet1 <- sheet_all[keep, , drop = FALSE]
      beta1 <- beta_all[, keep, drop = FALSE]

      grp <- factor(as.character(sheet1[[input$group_col]]), levels = c(input$ref_group, input$comp_group))
      sx  <- factor(as.character(sheet1[[sc]]), levels = c(input$ref_sex, input$comp_sex))
      validate(need(all(table(grp, sx) >= 2), "Each group-by-sex cell needs at least 2 samples."))

      ## ---- M-value conversion (same convention as mod_methyl_dmp.R's live
      ## engine): the model is fit on M-values; a beta-scale matrix is kept
      ## alongside for a delta-beta effect-size column.
      is_m_scale <- identical(methyl_dataset$input_scale, "m")
      beta_scale <- if (is_m_scale) 2^beta1 / (1 + 2^beta1) else beta1
      m <- if (is_m_scale) {
        beta1
      } else {
        b <- pmin(pmax(beta1, 1e-6), 1 - 1e-6)
        log2(b / (1 - b))
      }

      design <- stats::model.matrix(~ grp * sx)
      fit <- limma::eBayes(methyl_chunked_lmfit(m, design))
      coef_name <- utils::tail(colnames(design), 1)
      tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
      tt$cpg <- rownames(tt)
      rownames(tt) <- NULL
      tt <- tt[, c("cpg", setdiff(colnames(tt), "cpg"))]

      ## Delta-beta companion: the interaction contrast in beta-scale terms -
      ## (compGroup-refGroup difference in compSex) minus (compGroup-refGroup
      ## difference in refSex) - for interpretability alongside the
      ## M-value-scale logFC that the actual statistical test used.
      cell_ref_refsex  <- grp == levels(grp)[1] & sx == levels(sx)[1]
      cell_comp_refsex <- grp == levels(grp)[2] & sx == levels(sx)[1]
      cell_ref_compsex  <- grp == levels(grp)[1] & sx == levels(sx)[2]
      cell_comp_compsex <- grp == levels(grp)[2] & sx == levels(sx)[2]
      mean_beta <- function(mask) rowMeans(beta_scale[, mask, drop = FALSE], na.rm = TRUE)
      dbeta_interaction <- (mean_beta(cell_comp_compsex) - mean_beta(cell_ref_compsex)) -
        (mean_beta(cell_comp_refsex) - mean_beta(cell_ref_refsex))
      tt$dbeta_interaction <- dbeta_interaction[tt$cpg]

      list(table = tt, coef_name = coef_name, is_m_scale = is_m_scale)
    }, ignoreInit = TRUE)

    sig_table <- reactive({
      res <- fit_result()
      req(res)
      res$table %>% mutate(significant = adj.P.Val < input$padj_cut)
    })

    ## Publishes the fitted interaction model into shared methyl_results$interaction
    ## so ArthOChat and the biomarker card can see it, same field shape as the
    ## transcriptomics module's results$interaction contract (contrast label,
    ## n tested/significant, top hits) - silently skipped if the fit failed
    ## validation, same as mod_interaction.R.
    observeEvent(input$run_btn, {
      df <- tryCatch(sig_table(), error = function(e) NULL)
      req(df)
      res <- fit_result()
      methyl_results$interaction <- list(
        contrast = sprintf("%s vs %s interaction with sex (%s vs %s)",
                            input$comp_group, input$ref_group, input$comp_sex, input$ref_sex),
        coef_name = res$coef_name,
        n_tested = nrow(df),
        n_significant = sum(df$significant),
        top_hits = head(df$cpg[order(df$adj.P.Val)], 10),
        timestamp = Sys.time()
      )
    }, ignoreInit = TRUE)

    output$summary_ui <- renderUI({
      if (!int_has_run()) {
        return(div(class = "empty-note", icon("circle-info"),
                    "Not run yet. Set the group column, groups, and sexes on the left, then click \"Run interaction model\"."))
      }
      res <- fit_result()
      df <- sig_table()
      tagList(
        p(strong("Interaction term: "), res$coef_name),
        p(class = "submodule-desc", sprintf("Model fit on %s.", if (isTRUE(res$is_m_scale)) "M-values" else "M-values (logit-transformed from beta)")),
        p(strong(format(nrow(df), big.mark = ",")), " CpGs tested, ",
          strong(format(sum(df$significant), big.mark = ",")), " with a significant group-by-sex interaction at the current cutoff.")
      )
    })

    output$int_table <- DT::renderDataTable({
      if (!int_has_run()) return(NULL)
      DT::datatable(sig_table(), rownames = FALSE, filter = "top",
                     options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact") %>%
        DT::formatSignif(columns = c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "dbeta_interaction"), digits = 4)
    })

    output$download_int <- downloadHandler(
      filename = function() "methylation_sex_interaction.csv",
      content = function(file) write.csv(sig_table(), file, row.names = FALSE)
    )
  })
}
