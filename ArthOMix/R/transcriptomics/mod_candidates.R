## R/mod_candidates.R
## Candidate Gene Identification (Section 2.5): intersects a shared disease-associated
## WGCNA module background with the female and male DEG lists separately, one Venn
## diagram per sex. Reads WGCNA modules and DGE runs live from the shared `results` store.

mod_candidates_config <- list(
  id = "candidates", group = "Network",
  title = "Candidate Gene Identification",
  description = "Intersect the disease-associated WGCNA module with the sex-based DEG list, and visualise and export the Venn diagram and candidate table.",
  icon = "star"
)

## One sex panel (DEG picker, Venn vs. WGCNA background, candidate table), shared
## template for both "female" and "male" prefixes, gated behind its own run button.
mod_candidates_sex_panel_ui <- function(ns, prefix, title) {
  box(
    width = NULL, title = title, status = "primary", solidHeader = FALSE,
    uiOutput(ns(paste0(prefix, "_deg_picker_ui"))),
    actionButton(ns(paste0(prefix, "_run_btn")), sprintf("Compute %s candidates", prefix), icon = icon("play"), class = "btn-primary btn-sm"),
    div(style = "margin-top:10px;",
        withSpinner(uiOutput(ns(paste0(prefix, "_summary_ui"))), color = "#2c6fbb", type = 6),
        withSpinner(plotOutput(ns(paste0(prefix, "_venn")), height = 340), color = "#2c6fbb", type = 6),
        div(class = "table-toolbar",
            downloadButton(ns(paste0(prefix, "_download_venn")), "Diagram (PNG)", class = "btn-sm"),
            downloadButton(ns(paste0(prefix, "_download_table")), "Candidates (CSV)", class = "btn-sm")),
        DT::dataTableOutput(ns(paste0(prefix, "_table")))
    )
  )
}

## Body is server-rendered so it can be gated on WGCNA + DGE prereqs (see server below).
mod_candidates_ui <- function(id) {
  ns <- NS(id)
  uiOutput(ns("body_ui"))
}

mod_candidates_server <- function(id, dataset, results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Hard gate: requires both an upstream WGCNA run and at least one DGE run.
    prereqs <- reactive({
      list(
        wgcna_ok = !is.null(results$wgcna) && length(results$wgcna$module_genes) > 0,
        dge_ok = !is.null(results$dge_runs) && length(results$dge_runs) > 0
      )
    })

    output$body_ui <- renderUI({
      pr <- prereqs()
      if (!pr$wgcna_ok || !pr$dge_ok) {
        missing <- c(
          if (!pr$wgcna_ok) "WGCNA - run Step 3 (Modules) and Step 4 (Module-Trait) on the WGCNA Co-expression Network tab",
          if (!pr$dge_ok) "Differential Expression - run at least one contrast on the Differential Expression tab"
        )
        return(
          box(
            width = 12, title = "Run WGCNA and Differential Expression first", status = "warning", solidHeader = FALSE,
            p(class = "submodule-desc", "This tab intersects a WGCNA disease-module background with sex-stratified DEG lists - it cannot compute anything until both have actually been run this session. Still needed:"),
            tags$ul(lapply(missing, tags$li)),
            p(class = "submodule-desc", "Once both are run, the module picker and \"Compute ... candidates\" buttons below become available - no need to reload.")
          )
        )
      }
      tagList(
        box(
          width = 12, title = "WGCNA module background", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "The DEG data is taken from results$dge_runs and WGCNA data from results$wgcna in order to perform this sub-module."),
          uiOutput(ns("module_picker_ui"))
        ),
        uiOutput(ns("gene_panel_box_ui")),
        fluidRow(
          column(6, mod_candidates_sex_panel_ui(ns, "female", "Female candidates (module background ∩ female DEGs)")),
          column(6, mod_candidates_sex_panel_ui(ns, "male", "Male candidates (module background ∩ male DEGs)"))
        ),
        box(
          width = 12, title = "Final candidate gene set", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Pick which of the two panels above (or their overlap) is \"the\" candidate gene set going forward - this choice is what other tabs in the app read from results$candidates$final."),
          uiOutput(ns("final_set_picker_ui")),
          withSpinner(uiOutput(ns("final_summary_ui")), color = "#2c6fbb", type = 6),
          div(class = "table-toolbar", downloadButton(ns("final_download_table"), "Final set (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns("final_table"))
        )
      )
    })

    ## ---- Shared WGCNA module background ----
    module_choices <- reactive({
      wg <- results$wgcna
      req(wg, length(wg$module_genes) > 0)
      mods <- setdiff(names(wg$module_genes), "grey")
      validate(need(length(mods) > 0, "No non-grey modules were detected in the WGCNA run."))
      list(mods = mods, sizes = lengths(wg$module_genes[mods]),
           disease = intersect(wg$significant_trait_modules %||% character(0), mods))
    })

    ## Multi-select module picker, pre-selected to the disease-associated modules
    ## from WGCNA Step 4 (|cor| >= 0.5, p < 1e-8), but any combination is selectable.
    output$module_picker_ui <- renderUI({
      mc <- tryCatch(module_choices(), error = function(e) NULL)
      if (is.null(mc)) {
        return(div(class = "empty-note", icon("circle-info"),
          "Run WGCNA first (Step 3, Modules, then Step 4, Module-Trait) - detected co-expression modules will appear here."))
      }
      choices <- stats::setNames(
        mc$mods,
        sprintf("%s (n=%d)%s", mc$mods, mc$sizes, ifelse(mc$mods %in% mc$disease, " - disease-associated", ""))
      )
      default_sel <- if (length(mc$disease) > 0) mc$disease else mc$mods[1]
      tagList(
        selectInput(ns("wgcna_module_choice"), NULL, choices = choices, selected = default_sel, multiple = TRUE),
        if (length(mc$disease) > 0) p(class = "submodule-desc",
          sprintf("Pre-selected: %s - the disease-associated module(s) from WGCNA Step 4 (|cor| ≥ 0.5, p < 1e-8). Add or remove modules above to use a different combination.",
                  paste(mc$disease, collapse = " + ")))
      )
    })

    module_background <- reactive({
      wg <- results$wgcna
      req(wg, length(input$wgcna_module_choice) > 0)
      chosen <- intersect(input$wgcna_module_choice, names(wg$module_genes))
      req(length(chosen) > 0)
      unique(unlist(wg$module_genes[chosen], use.names = FALSE))
    })

    ## Optional third gene set (bundled panel or pasted custom list); NULL/"(none)" is a
    ## no-op. Upload-only: hidden unless dataset$source starts with "Uploaded dataset".
    output$gene_panel_box_ui <- renderUI({
      req(grepl("^Uploaded dataset", dataset$source %||% ""))
      box(
        width = 12, title = "Narrow further with a gene panel (optional)", status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "Off by default - the module/DEG overlap above is unchanged unless you pick a panel here. When set, candidates are additionally intersected with this gene list, e.g. narrowing a general disease-module overlap down to genes in a specific biological process (ferroptosis, autophagy, whatever the panel covers) - the same 3-way intersection a published study's own curated gene list would apply."),
        uiOutput(ns("gene_panel_picker_ui"))
      )
    })

    output$gene_panel_picker_ui <- renderUI({
      panels <- list_gene_panels()
      choices <- c("(none)" = "", panels, "Paste my own list" = "__custom__")
      tagList(
        selectInput(ns("gene_panel_choice"), "Gene panel", choices = choices, selected = "", selectize = FALSE),
        conditionalPanel(
          condition = sprintf("input['%s'] == '__custom__'", ns("gene_panel_choice")),
          textAreaInput(ns("gene_panel_custom"), "Gene symbols (one per line, or comma/space-separated)", rows = 5,
                        placeholder = "GPX4\nSLC7A11\nACSL4\n...")
        )
      )
    })

    gene_panel <- reactive({
      choice <- input$gene_panel_choice %||% ""
      if (!nzchar(choice)) return(NULL)
      if (identical(choice, "__custom__")) {
        raw <- input$gene_panel_custom %||% ""
        genes <- trimws(strsplit(raw, "[,\\s\\n]+", perl = TRUE)[[1]])
        genes <- unique(genes[nzchar(genes)])
        validate(need(length(genes) > 0, "Paste at least one gene symbol, or pick \"(none)\" to turn this off."))
        return(genes)
      }
      validate(need(file.exists(choice), "That bundled gene panel file could not be found on disk."))
      load_gene_panel(choice)
    })

    ## DEG contrast pickers, one per sex, listing every DGE run and defaulting to
    ## whichever contrast label matches that sex (word or "F"/"M" code, word-boundary).
    deg_run_choices <- reactive({
      runs <- results$dge_runs
      req(length(runs) > 0)
      ids <- names(runs)
      stats::setNames(ids, vapply(ids, function(i) runs[[i]]$contrast, character(1)))
    })

    guess_run <- function(ch, pattern) {
      hit <- ch[grepl(pattern, names(ch), ignore.case = TRUE, perl = TRUE)]
      unname(utils::tail(if (length(hit) > 0) hit else ch, 1))
    }

    output$female_deg_picker_ui <- renderUI({
      ch <- tryCatch(deg_run_choices(), error = function(e) NULL)
      if (is.null(ch)) {
        return(div(class = "empty-note", icon("circle-info"),
          "Run Differential Expression for the female stratum first (e.g. contrast column \"group\", second column \"sex\" filtered to F)."))
      }
      selectInput(ns("female_deg_run"), "Female DEG contrast", choices = ch, selected = guess_run(ch, "\\bfemale\\b|\\bF\\b"), selectize = FALSE)
    })

    output$male_deg_picker_ui <- renderUI({
      ch <- tryCatch(deg_run_choices(), error = function(e) NULL)
      if (is.null(ch)) {
        return(div(class = "empty-note", icon("circle-info"),
          "Run Differential Expression for the male stratum first (e.g. contrast column \"group\", second column \"sex\" filtered to M)."))
      }
      selectInput(ns("male_deg_run"), "Male DEG contrast", choices = ch, selected = guess_run(ch, "\\bmale\\b|\\bM\\b"), selectize = FALSE)
    })

    ## Module background ∩ one sex's significant DEGs (+ optional gene panel), gated
    ## behind that sex's run button so changing inputs doesn't recompute automatically.
    sex_candidates <- function(deg_input_name, run_btn_name) {
      eventReactive(input[[run_btn_name]], {
        req(input[[deg_input_name]])
        bg <- module_background()
        validate(need(length(bg) > 0, "Pick at least one WGCNA module above."))
        run <- results$dge_runs[[input[[deg_input_name]]]]
        req(run)
        tab <- run$table
        sig <- tab[!is.na(tab$adj.P.Val) & tab$direction != "Not significant", , drop = FALSE]
        deg_genes <- sig$gene
        validate(need(length(deg_genes) > 0, "This DEG contrast has no significant genes at its current cutoffs."))

        overlap <- intersect(bg, deg_genes)
        validate(need(length(overlap) > 0, "No genes are shared between the module background and this DEG list."))

        panel <- tryCatch(gene_panel(), error = function(e) NULL)
        if (!is.null(panel)) {
          overlap <- intersect(overlap, panel)
          validate(need(length(overlap) > 0, "No genes are shared between the module/DEG overlap and the selected gene panel."))
        }

        expr <- dataset$expr
        in_expr <- overlap %in% rownames(expr)
        cand <- data.frame(gene = overlap, stringsAsFactors = FALSE)
        m <- match(cand$gene, tab$gene)
        cand$logFC <- round(tab$logFC[m], 3)
        cand$adj.P.Val <- signif(tab$adj.P.Val[m], 3)
        cand$direction <- tab$direction[m]
        wg <- results$wgcna
        mod_lookup <- stats::setNames(rep(names(wg$module_genes), lengths(wg$module_genes)), unlist(wg$module_genes, use.names = FALSE))
        cand$wgcna_module <- unname(mod_lookup[cand$gene])
        cand$in_current_dataset <- in_expr
        cand$mean_expr <- NA_real_
        if (any(in_expr)) cand$mean_expr[in_expr] <- rowMeans(expr[overlap[in_expr], , drop = FALSE])

        ## One-sided hypergeometric test: is the overlap bigger than chance predicts.
        N <- length(union(rownames(expr), union(bg, deg_genes)))
        p_value <- stats::phyper(length(overlap) - 1, length(bg), N - length(bg), length(deg_genes), lower.tail = FALSE)

        sets <- list(`WGCNA module(s)` = bg, `Significant DEGs` = deg_genes)
        if (!is.null(panel)) sets[["Gene panel"]] <- panel

        list(
          sets = sets,
          overlap = overlap, stats = cand, contrast = run$contrast,
          n_bg = length(bg), n_deg = length(deg_genes), n_panel = if (!is.null(panel)) length(panel) else NA_integer_,
          p_value = p_value
        )
      }, ignoreInit = TRUE)
    }

    female_result <- sex_candidates("female_deg_run", "female_run_btn")
    male_result <- sex_candidates("male_deg_run", "male_run_btn")

    ## Registers the summary/Venn/table/download outputs for one sex panel.
    register_panel <- function(prefix, res) {
      output[[paste0(prefix, "_summary_ui")]] <- renderUI({
        r <- tryCatch(res(), error = function(e) NULL)
        if (is.null(r)) {
          return(div(class = "empty-note", icon("circle-info"),
            sprintf("Not run yet. Click \"Compute %s candidates\" above.", prefix)))
        }
        tagList(
          p(strong(r$n_bg), " module-background genes, ", strong(r$n_deg), " significant DEGs (", r$contrast, ")",
            if (!is.na(r$n_panel)) tagList(", ", strong(r$n_panel), " genes in the selected gene panel") else NULL, "."),
          p(strong(length(r$overlap)), " candidate genes in the overlap - hypergeometric enrichment ",
            HTML("<em>p</em>"), " = ", signif(r$p_value, 3), ".")
        )
      })

      ## Green for Female, brown for Male, matching this project's reference figures.
      venn_fill_high <- c(female = "#1a7a3c", male = "#7a4a26")[[prefix]]

      venn_obj <- reactive({
        r <- tryCatch(res(), error = function(e) NULL)
        req(r)
        draw_overlap_venn(r$sets, title = sprintf("%s: %d candidates", tools::toTitleCase(prefix), length(r$overlap)),
                           fill_high = venn_fill_high)
      })

      output[[paste0(prefix, "_venn")]] <- renderPlot({ venn_obj() })

      output[[paste0(prefix, "_download_venn")]] <- downloadHandler(
        filename = function() sprintf("%s_candidate_venn.png", prefix),
        content = function(file) ggsave(file, plot = venn_obj(), width = 6.5, height = 5.5, dpi = 300, bg = "white")
      )

      output[[paste0(prefix, "_table")]] <- DT::renderDataTable({
        r <- tryCatch(res(), error = function(e) NULL)
        req(r)
        DT::datatable(r$stats, rownames = FALSE, filter = "top",
                       options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
      })

      output[[paste0(prefix, "_download_table")]] <- downloadHandler(
        filename = function() sprintf("%s_candidate_genes.csv", prefix),
        content = function(file) write.csv(res()$stats, file, row.names = FALSE)
      )
    }

    register_panel("female", female_result)
    register_panel("male", male_result)

    ## Final candidate gene set: lets the user pick female/male/union/intersection as
    ## "the" set for downstream tabs. Defaults to Union, matching this project's own
    ## MR script (00_shared/10_MR.R), which runs MR once on union(fem, mal).
    output$final_set_picker_ui <- renderUI({
      radioButtons(
        ns("final_candidate_set"), NULL,
        choices = c("Female candidates only" = "female",
                    "Male candidates only" = "male",
                    "Union (candidate in either sex)" = "union",
                    "Intersection (candidate in both sexes)" = "intersection"),
        selected = "union", inline = TRUE
      )
    })

    final_candidates <- reactive({
      req(input$final_candidate_set)
      fr <- tryCatch(female_result(), error = function(e) NULL)
      mr <- tryCatch(male_result(), error = function(e) NULL)
      switch(input$final_candidate_set,
        female = {
          validate(need(!is.null(fr), "Female candidates are not available yet - see the Female panel above."))
          list(genes = fr$overlap, stats = fr$stats)
        },
        male = {
          validate(need(!is.null(mr), "Male candidates are not available yet - see the Male panel above."))
          list(genes = mr$overlap, stats = mr$stats)
        },
        union = {
          validate(need(!is.null(fr) && !is.null(mr), "Both the Female and Male panels above must have results before their union can be computed."))
          genes <- union(fr$overlap, mr$overlap)
          combined_stats <- rbind(fr$stats, mr$stats[!mr$stats$gene %in% fr$stats$gene, , drop = FALSE])
          list(genes = genes, stats = combined_stats[combined_stats$gene %in% genes, , drop = FALSE])
        },
        intersection = {
          validate(need(!is.null(fr) && !is.null(mr), "Both the Female and Male panels above must have results before their intersection can be computed."))
          genes <- intersect(fr$overlap, mr$overlap)
          validate(need(length(genes) > 0, "No candidate genes are shared between the female and male lists."))
          list(genes = genes, stats = fr$stats[fr$stats$gene %in% genes, , drop = FALSE])
        }
      )
    })

    output$final_summary_ui <- renderUI({
      fc <- final_candidates()
      p(strong(length(fc$genes)), " genes in the final candidate set (",
        switch(input$final_candidate_set, female = "female only", male = "male only", union = "found in either sex", intersection = "shared by both sexes"),
        ").")
    })

    output$final_table <- DT::renderDataTable({
      DT::datatable(final_candidates()$stats, rownames = FALSE, filter = "top",
                     options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })

    output$final_download_table <- downloadHandler(
      filename = function() sprintf("final_candidate_genes_%s.csv", input$final_candidate_set),
      content = function(file) write.csv(final_candidates()$stats, file, row.names = FALSE)
    )

    ## Read by the Assistant sub-module, same as every other analysis tab.
    observe({
      fr <- tryCatch(female_result(), error = function(e) NULL)
      mr <- tryCatch(male_result(), error = function(e) NULL)
      fc <- tryCatch(final_candidates(), error = function(e) NULL)
      if (is.null(fr) && is.null(mr)) return()
      results$candidates <- list(
        female = if (!is.null(fr)) list(n_candidates = length(fr$overlap), genes = fr$overlap, contrast = fr$contrast) else NULL,
        male = if (!is.null(mr)) list(n_candidates = length(mr$overlap), genes = mr$overlap, contrast = mr$contrast) else NULL,
        final = if (!is.null(fc)) list(selection = input$final_candidate_set, n_candidates = length(fc$genes), genes = fc$genes) else NULL
      )
    })
  })
}
