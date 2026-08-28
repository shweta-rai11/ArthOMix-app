## R/mod_deconvolution.R
## Immune Deconvolution: CIBERSORT (LM22 via IOBR) as the primary fraction
## estimator, MCP-counter as an independent corroborating check, then a
## group-composition test (Wilcoxon/Kruskal-Wallis, BH FDR).

mod_deconvolution_config <- list(
  id = "deconvolution", group = "Interpretation",
  title = "Immune Deconvolution",
  description = "Estimate immune cell composition (CIBERSORT/LM22, corroborated with MCP-counter) for pre-loaded expression matrix or uploaded data.",
  icon = "chart-pie"
)

mod_deconvolution_ui <- function(id) {
  ns <- NS(id)
  tagList(
      p(class = "submodule-desc", "Runs CIBERSORT (LM22 signature, via IOBR) as the primary estimator and MCP-counter as an independent corroborating check, on whatever dataset is currently loaded; gene symbols must be in the expression matrix rownames."),
      uiOutput(ns("advanced_ui")),
      actionButton(ns("run_btn"), "Estimate cell composition", icon = icon("play"), class = "btn-primary btn-sm"),
      div(class = "submodule-desc", style = "margin-top:6px;", "CIBERSORT's permutation testing can take up to a minute or two for larger cohorts."),
      uiOutput(ns("not_run_note")),
      uiOutput(ns("id_overlap_warning_ui")),
      br(),
      withSpinner(plotOutput(ns("cell_plot"), height = 420), color = "#2c6fbb", type = 6),
      box(
        width = 12, title = "Estimated cell fractions (CIBERSORT, LM22)", status = "primary", solidHeader = FALSE,
        div(class = "table-toolbar", downloadButton(ns("download_cell"), "Download CSV", class = "btn-sm")),
        DT::dataTableOutput(ns("cell_table"))
      ),
      box(
        width = 12, title = "Group comparison (CIBERSORT)", status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "Per cell type: does the CIBERSORT fraction differ across the comparison groups? Benjamini-Hochberg FDR-adjusted across all cell types tested."),
        div(class = "table-toolbar", downloadButton(ns("download_stats"), "Download CSV", class = "btn-sm")),
        DT::dataTableOutput(ns("stats_table"))
      ),
      tags$details(
        class = "box box-primary",
        tags$summary(class = "box-header", style = "cursor: pointer;",
                      tags$h3(class = "box-title", icon("magnifying-glass-chart"), "Corroborating check: MCP-counter")),
        div(class = "box-body",
          p(class = "submodule-desc", "MCP-counter is an independent, marker-gene-based estimator. It gives relative abundance scores, not fractions, so it's shown here to corroborate direction - agreement between two methods that fail in different ways is worth more than either alone - never as a replacement for the CIBERSORT fractions above."),
          withSpinner(plotOutput(ns("mcp_plot"), height = 380), color = "#2c6fbb", type = 6),
          div(class = "table-toolbar", downloadButton(ns("download_mcp"), "Download scores CSV", class = "btn-sm")),
          DT::dataTableOutput(ns("mcp_table")),
          tags$hr(),
          p(class = "submodule-desc", style = "margin-bottom:4px;", "Group comparison (MCP-counter), same test and FDR as above:"),
          div(class = "table-toolbar", downloadButton(ns("download_mcp_stats"), "Download CSV", class = "btn-sm")),
          DT::dataTableOutput(ns("mcp_stats_table"))
        )
      )
  )
}

## % of gene_ids that are known human HUGO symbols (org.Hs.eg.db), the
## identifier type CIBERSORT (LM22) and MCP-counter both require. Neither
## package validates this itself - they silently intersect their own marker
## gene lists against whatever names they're given, so a matrix keyed by
## Ensembl IDs, probe IDs, or an unrecognised casing previously produced a
## low-confidence or degenerate result with no signal to the user short of a
## hard failure in the rare case IOBR returns zero fraction columns. A
## top-level function (not nested in mod_deconvolution_server) so it can be
## unit-tested without running the full CIBERSORT/MCP-counter pipeline.
deconv_gene_id_overlap_pct <- function(gene_ids) {
  known_symbols <- tryCatch(AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL"),
                             error = function(e) NULL)
  if (is.null(known_symbols)) return(NA_real_)
  round(100 * mean(gene_ids %in% known_symbols), 1)
}

mod_deconvolution_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Non-negative values with 99th percentile > 100 flag linear-scale
    ## expression (same scale-detection rule as mod_dge.R/mod_preprocessing.R).
    is_linear_counts <- function(expr) {
      vals <- as.numeric(expr)
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0 || any(vals < 0)) return(FALSE)
      q99 <- suppressWarnings(stats::quantile(vals[vals > 0], 0.99, na.rm = TRUE))
      isTRUE(!is.na(q99) && q99 > 100)
    }

    ## Extract LM22 fraction columns from IOBR's CIBERSORT output, dropping
    ## the P-value/Correlation/RMSE diagnostic columns.
    cib_frac_cols <- function(cib_raw) {
      all_cib <- grep("_CIBERSORT$", colnames(cib_raw), value = TRUE)
      diag_cols <- grep("^(P-value|Correlation|RMSE)_CIBERSORT$", all_cib, value = TRUE)
      setdiff(all_cib, diag_cols)
    }
    clean_cib_name <- function(x) gsub("_", " ", sub("_CIBERSORT$", "", x))

    ## Metadata columns with 2-20 distinct values, i.e. plausible grouping
    ## columns (same rule as mod_dge.R's candidate_columns()).
    candidate_columns <- reactive({
      meta <- dataset$meta
      req(meta)
      cols <- setdiff(colnames(meta), "sample")
      keep <- vapply(cols, function(cl) {
        nu <- length(unique(stats::na.omit(meta[[cl]])))
        nu >= 2 && nu <= 20
      }, logical(1))
      cols[keep]
    })

    ## Collapsed by default (native <details>, since box(collapsible=TRUE)
    ## fails once inserted via renderUI); all controls are optional and
    ## defaults reproduce the standard one-click comparison.
    output$advanced_ui <- renderUI({
      cols <- candidate_columns()
      validate(need(length(cols) > 0, "The loaded metadata has no column with 2-20 distinct values to compare cell composition by."))
      default_group <- if ("group" %in% cols) "group" else cols[1]

      tags$details(
        class = "box box-primary",
        tags$summary(class = "box-header", style = "cursor: pointer;",
                      tags$h3(class = "box-title", icon("sliders"), "Advanced: comparison group, sample filter & test")),
        div(class = "box-body",
          p(class = "submodule-desc",
            "Defaults to comparing by \"group\" with Wilcoxon rank-sum (2 groups) or Kruskal-Wallis (3+ groups), FDR-adjusted - the same test this project's own deconvolution pipeline uses. Open this only if you want to compare by a different column, restrict to a sample subset first, force a specific test, or tune CIBERSORT's permutation count."),
          selectInput(ns("group_col"), "Compare cell composition across", choices = cols, selected = default_group, selectize = FALSE),
          radioButtons(
            ns("test_method"), "Statistical test",
            choiceNames = list("Auto (recommended - Wilcoxon for 2 groups, Kruskal-Wallis for 3+)",
                                "Wilcoxon rank-sum (2 groups only)",
                                "Kruskal-Wallis (2 or more groups)"),
            choiceValues = list("auto", "wilcox", "kruskal"), selected = "auto"
          ),
          ## 100 matches CFG$perm in 05c_deconvolution.R; 0 skips permutation p-values for a faster run.
          numericInput(ns("cib_perm"), "CIBERSORT permutations (higher = slower, more precise deconvolution p-values; 0 = skip)",
                       value = 100, min = 0, max = 1000, step = 50),
          tags$hr(),
          div(class = "filter-section-header", icon("users"), "Sample filter (optional)"),
          selectInput(ns("filter_col"), "Restrict to samples where", choices = c("(no filter)", cols), selected = "(no filter)", selectize = FALSE),
          conditionalPanel(
            condition = sprintf("input['%s'] != '(no filter)'", ns("filter_col")),
            uiOutput(ns("filter_val_ui"))
          )
        )
      )
    })

    ## All levels pre-checked by default, so opening the filter picker keeps
    ## every sample until something is unchecked.
    output$filter_val_ui <- renderUI({
      req(input$filter_col, !identical(input$filter_col, "(no filter)"))
      meta <- dataset$meta
      req(input$filter_col %in% colnames(meta))
      vals <- sort(unique(stats::na.omit(as.character(meta[[input$filter_col]]))))
      checkboxGroupInput(ns("filter_vals"), "...equals one of", choices = vals, selected = vals, inline = TRUE)
    })

    result <- eventReactive(input$run_btn, {
      expr <- dataset$expr
      validate(need(nrow(expr) >= 100, "Too few genes in the current expression matrix to run a deconvolution meaningfully."))

      meta <- dataset$meta
      group_col <- input$group_col %||% "group"
      validate(need(group_col %in% colnames(meta), sprintf("Column \"%s\" was not found in the loaded metadata.", group_col)))

      if (!is.null(input$filter_col) && !identical(input$filter_col, "(no filter)") && length(input$filter_vals) > 0) {
        meta <- meta[as.character(meta[[input$filter_col]]) %in% input$filter_vals, , drop = FALSE]
      }
      meta <- meta[!is.na(meta[[group_col]]), , drop = FALSE]

      common <- intersect(colnames(expr), meta$sample)
      validate(need(length(common) >= 4, "Fewer than 4 samples remain after the sample filter (and have a value for the comparison column); loosen the filter, or pick a different comparison column."))
      expr <- expr[, common, drop = FALSE]
      meta <- meta[match(common, meta$sample), , drop = FALSE]
      keep_cols <- intersect(unique(c("sample", "group", "sex", group_col)), colnames(meta))
      meta_sub <- meta[, keep_cols, drop = FALSE]

      ## Neither CIBERSORT (LM22) nor MCP-counter validate that rownames(expr)
      ## are actually HUGO gene symbols - both packages silently intersect
      ## their own marker gene lists against whatever names are given, so a
      ## matrix keyed by Ensembl IDs, probe IDs, or an unrecognised casing
      ## would previously produce a low-confidence or degenerate result with
      ## no signal to the user beyond a hard failure in the rare case IOBR
      ## returns literally zero fraction columns. Check the overlap against
      ## the known HUGO symbol universe up front instead.
      id_overlap_pct <- deconv_gene_id_overlap_pct(rownames(expr))
      validate(need(is.na(id_overlap_pct) || id_overlap_pct >= 10,
                    sprintf("Only %s%% of this expression matrix's row names match known human gene symbols - CIBERSORT/MCP-counter need HUGO gene symbols as row names (not Ensembl IDs, probe IDs, or another identifier type). Check the loaded dataset's feature IDs.",
                            id_overlap_pct)))

      ## --- CIBERSORT (LM22, primary) ---------------------------------
      is_counts <- is_linear_counts(expr)
      lin_expr <- if (is_counts) as.matrix(expr) else 2^as.matrix(expr)
      lin_expr[!is.finite(lin_expr)] <- 0
      perm <- suppressWarnings(as.integer(input$cib_perm %||% 100))
      if (is.na(perm) || perm < 0) perm <- 100

      cib_raw <- tryCatch(
        suppressMessages(suppressWarnings(
          IOBR::deconvo_tme(eset = lin_expr, method = "cibersort", arrays = !is_counts, perm = perm)
        )),
        error = function(e) validate(need(FALSE, paste("CIBERSORT could not run on this expression matrix:", conditionMessage(e))))
      )
      cib_raw <- as.data.frame(cib_raw)
      colnames(cib_raw)[colnames(cib_raw) == "ID"] <- "sample"
      cib_frac <- cib_frac_cols(cib_raw)
      validate(need(length(cib_frac) > 0, "CIBERSORT returned no LM22 fraction columns for this expression matrix."))
      cib_df <- cib_raw[, c("sample", cib_frac), drop = FALSE]
      colnames(cib_df)[colnames(cib_df) != "sample"] <- clean_cib_name(cib_frac)
      cib_cell_cols <- setdiff(colnames(cib_df), "sample")
      cib_table <- merge(meta_sub, cib_df, by = "sample")

      ## --- MCP-counter (independent check) ---------------------------
      ## Each score is a mean of marker genes (MCPcounter:::appendSignatures),
      ## which needs log-scale expression; log2-transform if linear.
      mcp_input <- if (is_counts) log2(as.matrix(expr) + 1) else as.matrix(expr)
      mcp_est <- tryCatch(
        MCPcounter::MCPcounter.estimate(as.data.frame(mcp_input), featuresType = "HUGO_symbols"),
        error = function(e) validate(need(FALSE, paste("MCP-counter could not run on this expression matrix:", conditionMessage(e))))
      )
      mcp_df <- as.data.frame(t(mcp_est))
      mcp_cell_cols <- colnames(mcp_df)
      mcp_df$sample <- rownames(mcp_df)
      mcp_table <- merge(meta_sub, mcp_df, by = "sample")

      list(
        group_col = group_col,
        cib = list(table = cib_table, cell_cols = cib_cell_cols),
        mcp = list(table = mcp_table, cell_cols = mcp_cell_cols),
        id_overlap_pct = id_overlap_pct
      )
    }, ignoreInit = TRUE)

    output$id_overlap_warning_ui <- renderUI({
      req(deconv_has_run())
      pct <- tryCatch(result()$id_overlap_pct, error = function(e) NA_real_)
      req(!is.na(pct), pct < 70)
      div(class = "empty-note", style = "margin-top: 8px;", icon("triangle-exclamation"),
          sprintf("Only %s%% of this expression matrix's row names matched known human gene symbols. CIBERSORT and MCP-counter can only score genes they recognise, so the fractions below may be based on a small, possibly unrepresentative subset of the LM22/MCP-counter marker panels - check that the loaded dataset's feature IDs are gene symbols.", pct))
    })

    deconv_has_run <- reactiveVal(FALSE)
    observeEvent(input$run_btn, deconv_has_run(TRUE), ignoreInit = TRUE)

    output$not_run_note <- renderUI({
      if (deconv_has_run()) return(NULL)
      div(class = "empty-note", style = "margin-top: 8px;", icon("circle-info"),
          "Not run yet. Click \"Estimate cell composition\" above.")
    })

    ## Per-cell-type group-difference test, shared by both estimators so
    ## their stats tables/plots and the results summary agree on what's significant.
    compute_group_stats <- function(df, cell_cols, group_col, method_input) {
      lvls <- sort(unique(as.character(stats::na.omit(df[[group_col]]))))
      validate(need(length(lvls) >= 2, sprintf("\"%s\" has fewer than two distinct values among the samples analysed, so cell composition can't be compared across it.", group_col)))

      if (identical(method_input, "wilcox")) {
        validate(need(length(lvls) == 2, "Wilcoxon rank-sum needs exactly two groups. Pick Kruskal-Wallis instead, or restrict the sample filter to two levels of this column."))
      }
      used_method <- if (identical(method_input, "auto")) (if (length(lvls) == 2) "wilcox" else "kruskal") else method_input
      test_label <- if (identical(used_method, "wilcox")) "Wilcoxon rank-sum" else "Kruskal-Wallis"

      grp_full <- factor(as.character(df[[group_col]]), levels = lvls)
      rows <- lapply(cell_cols, function(ct) {
        x <- df[[ct]]
        ok <- is.finite(x) & !is.na(grp_full)
        g <- droplevels(grp_full[ok]); xv <- x[ok]
        p <- NA_real_; direction <- NA_character_
        if (nlevels(g) >= 2 && all(table(g) >= 2) && stats::sd(xv) > 0) {
          p <- tryCatch(
            suppressWarnings(
              if (identical(used_method, "wilcox") && nlevels(g) == 2) {
                stats::wilcox.test(xv ~ g)$p.value
              } else {
                stats::kruskal.test(xv, g)$p.value
              }
            ),
            error = function(e) NA_real_
          )
          if (nlevels(g) == 2) {
            m <- tapply(xv, g, mean)
            direction <- if (m[[2]] > m[[1]]) sprintf("Higher in %s", levels(g)[2]) else sprintf("Higher in %s", levels(g)[1])
          }
        }
        data.frame(cell_type = ct, p.value = p, direction = direction, stringsAsFactors = FALSE)
      })
      stats_df <- do.call(rbind, rows)
      stats_df$p.adj <- stats::p.adjust(stats_df$p.value, method = "BH")
      stats_df$significant <- !is.na(stats_df$p.adj) & stats_df$p.adj < 0.05
      stats_df <- stats_df[order(stats_df$p.value), , drop = FALSE]

      list(table = stats_df, test_label = test_label, group_col = group_col)
    }

    cib_stats <- reactive({
      res <- result()
      compute_group_stats(res$cib$table, res$cib$cell_cols, res$group_col, input$test_method %||% "auto")
    })
    mcp_stats <- reactive({
      res <- result()
      compute_group_stats(res$mcp$table, res$mcp$cell_cols, res$group_col, input$test_method %||% "auto")
    })

    ## Boxplot with significance stars; shared rendering so the CIBERSORT and
    ## MCP-counter panels stay visually consistent.
    render_group_boxplot <- function(df, cell_cols, group_col, gs, y_lab) {
      long <- tidyr::pivot_longer(df, cols = all_of(cell_cols), names_to = "cell_type", values_to = "score")
      p <- ggplot(long, aes(x = cell_type, y = score, fill = .data[[group_col]])) +
        geom_boxplot(outlier.size = 0.6) +
        labs(x = NULL, y = y_lab, fill = tools::toTitleCase(group_col)) +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 40, hjust = 1))

      if (!is.null(gs)) {
        st <- gs$table
        st$label <- ifelse(is.na(st$p.adj), "",
                     ifelse(st$p.adj < 0.001, "***",
                      ifelse(st$p.adj < 0.01, "**",
                       ifelse(st$p.adj < 0.05, "*", ""))))
        st <- st[nzchar(st$label), , drop = FALSE]
        if (nrow(st) > 0) {
          ymax <- stats::aggregate(score ~ cell_type, data = long, FUN = max, na.rm = TRUE)
          st <- merge(st, ymax, by = "cell_type")
          p <- p +
            ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.12))) +
            ggplot2::geom_text(data = st, aes(x = cell_type, y = score * 1.05, label = label),
                                inherit.aes = FALSE, size = 5, fontface = "bold", color = "#B3261E")
        }
      }
      p
    }

    ## Publish a run summary to results$deconvolution for cross-tab use
    ## (ArthOChat's session context in submodules_registry.R).
    observeEvent(input$run_btn, {
      if (is.null(results)) return(invisible(NULL))
      cgs <- tryCatch(cib_stats(), error = function(e) NULL)
      mgs <- tryCatch(mcp_stats(), error = function(e) NULL)
      req(cgs)
      results$deconvolution <- list(
        compared_across = cgs$group_col,
        cibersort_test = cgs$test_label,
        cibersort_n_cell_types = nrow(cgs$table),
        cibersort_n_significant = sum(cgs$table$significant, na.rm = TRUE),
        cibersort_significant_cell_types = cgs$table$cell_type[cgs$table$significant],
        mcpcounter_n_significant = if (!is.null(mgs)) sum(mgs$table$significant, na.rm = TRUE) else NA_integer_,
        mcpcounter_significant_cell_types = if (!is.null(mgs)) mgs$table$cell_type[mgs$table$significant] else character(0)
      )
    })

    ## --- CIBERSORT outputs (primary) -----------------------------------

    output$cell_table <- DT::renderDataTable({
      req(deconv_has_run())
      DT::datatable(result()$cib$table, rownames = FALSE, filter = "top",
                     options = list(pageLength = 12, scrollX = TRUE), class = "stripe hover compact")
    })

    output$cell_plot <- renderPlot({
      req(deconv_has_run())
      res <- result()
      gs <- tryCatch(cib_stats(), error = function(e) NULL)
      render_group_boxplot(res$cib$table, res$cib$cell_cols, res$group_col, gs, "CIBERSORT fraction")
    })

    output$stats_table <- DT::renderDataTable({
      req(deconv_has_run())
      gs <- tryCatch(cib_stats(), error = function(e) NULL)
      req(gs)
      st <- gs$table
      out <- data.frame(
        cell_type = st$cell_type, test = gs$test_label,
        direction = ifelse(is.na(st$direction), "", st$direction),
        p.value = signif(st$p.value, 3), FDR = signif(st$p.adj, 3),
        significant = ifelse(st$significant, "✓", ""),
        stringsAsFactors = FALSE
      )
      DT::datatable(out, rownames = FALSE, filter = "top",
                     options = list(pageLength = 12, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_cell <- downloadHandler(
      filename = function() "immune_deconvolution_cibersort.csv",
      content = function(file) write.csv(result()$cib$table, file, row.names = FALSE)
    )

    output$download_stats <- downloadHandler(
      filename = function() "immune_deconvolution_cibersort_group_tests.csv",
      content = function(file) write.csv(cib_stats()$table, file, row.names = FALSE)
    )

    ## --- MCP-counter outputs (corroborating check) ---------------------

    output$mcp_table <- DT::renderDataTable({
      req(deconv_has_run())
      DT::datatable(result()$mcp$table, rownames = FALSE, filter = "top",
                     options = list(pageLength = 12, scrollX = TRUE), class = "stripe hover compact")
    })

    output$mcp_plot <- renderPlot({
      req(deconv_has_run())
      res <- result()
      gs <- tryCatch(mcp_stats(), error = function(e) NULL)
      render_group_boxplot(res$mcp$table, res$mcp$cell_cols, res$group_col, gs, "MCP-counter score")
    })

    output$mcp_stats_table <- DT::renderDataTable({
      req(deconv_has_run())
      gs <- tryCatch(mcp_stats(), error = function(e) NULL)
      req(gs)
      st <- gs$table
      out <- data.frame(
        cell_type = st$cell_type, test = gs$test_label,
        direction = ifelse(is.na(st$direction), "", st$direction),
        p.value = signif(st$p.value, 3), FDR = signif(st$p.adj, 3),
        significant = ifelse(st$significant, "✓", ""),
        stringsAsFactors = FALSE
      )
      DT::datatable(out, rownames = FALSE, filter = "top",
                     options = list(pageLength = 12, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_mcp <- downloadHandler(
      filename = function() "immune_deconvolution_mcpcounter.csv",
      content = function(file) write.csv(result()$mcp$table, file, row.names = FALSE)
    )

    output$download_mcp_stats <- downloadHandler(
      filename = function() "immune_deconvolution_mcpcounter_group_tests.csv",
      content = function(file) write.csv(mcp_stats()$table, file, row.names = FALSE)
    )
  })
}
