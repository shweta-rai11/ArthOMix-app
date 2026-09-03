## R/transcriptomics/06_Candidate_Gene_Identification/mod_candidates.R
## Candidate Gene Identification (Section 2.5): intersects a shared disease-associated
## WGCNA module background with a DEG list, and visualises/exports the Venn diagram and

mod_candidates_config <- list(
  id = "candidates", group = "Network",
  title = "Candidate Gene Identification",
  description = "Intersect the disease-associated WGCNA module with the DEG list (sex-stratified when the data supports it), optionally refined by MR/Colocalization support, and visualise and export the Venn diagram and candidate table.",
  icon = "star"
)

mod_candidates_sex_panel_ui <- function(ns, prefix, title, btn_label = sprintf("Compute %s candidates", prefix)) {
  box(
    width = NULL, title = title, status = "primary", solidHeader = FALSE,
    uiOutput(ns(paste0(prefix, "_deg_picker_ui"))),
    actionButton(ns(paste0(prefix, "_run_btn")), btn_label, icon = icon("play"), class = "btn-primary btn-sm"),
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

mod_candidates_ui <- function(id) {
  ns <- NS(id)
  uiOutput(ns("body_ui"))
}

mod_candidates_server <- function(id, dataset, results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    prereqs <- reactive({
      list(
        wgcna_ok = !is.null(results$wgcna) && length(results$wgcna$module_genes) > 0,
        dge_ok = !is.null(results$dge_runs) && length(results$dge_runs) > 0
      )
    })

    sex_available <- reactive({
      m <- dataset$meta
      !is.null(m) && "sex" %in% names(m) && length(unique(stats::na.omit(m$sex))) >= 2
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
            p(class = "submodule-desc", "This tab intersects a WGCNA disease-module background with DEG lists - it cannot compute anything until both have actually been run this session. Still needed:"),
            tags$ul(lapply(missing, tags$li)),
            p(class = "submodule-desc", "Once both are run, the module picker and \"Compute ... candidates\" buttons below become available - no need to reload.")
          )
        )
      }
      shared_boxes <- tagList(
        box(
          width = 12, title = "WGCNA module background", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "The DEG data is taken from results$dge_runs and WGCNA data from results$wgcna in order to perform this sub-module."),
          uiOutput(ns("module_picker_ui"))
        ),
        uiOutput(ns("gene_panel_box_ui"))
      )
      if (!isTRUE(sex_available())) {
        return(tagList(
          shared_boxes,
          div(class = "empty-note", icon("circle-info"),
              "No male/female metadata was found for the currently loaded dataset (a GEO fetch or upload without a mapped sex column leaves this unset) - showing one pooled candidate panel instead of separate female/male ones."),
          mod_candidates_sex_panel_ui(ns, "pooled", "Candidate biomarkers (module background ∩ DEGs)", btn_label = "Compute candidates"),
          box(
            width = 12, title = "Final candidate gene set", status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "This is \"the\" candidate gene set going forward - what other tabs in the app (Feature Selection, Diagnostic Model, ...) read from results$candidates$final."),
            uiOutput(ns("causal_refine_ui")),
            withSpinner(uiOutput(ns("final_summary_ui")), color = "#2c6fbb", type = 6),
            div(class = "table-toolbar", downloadButton(ns("final_download_table"), "Final set (CSV)", class = "btn-sm")),
            DT::dataTableOutput(ns("final_table"))
          )
        ))
      }
      tagList(
        shared_boxes,
        fluidRow(
          column(6, mod_candidates_sex_panel_ui(ns, "female", "Female candidates (module background ∩ female DEGs)")),
          column(6, mod_candidates_sex_panel_ui(ns, "male", "Male candidates (module background ∩ male DEGs)"))
        ),
        box(
          width = 12, title = "Final candidate gene set", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Pick which of the two panels above (or their overlap) is \"the\" candidate gene set going forward - this choice is what other tabs in the app read from results$candidates$final."),
          uiOutput(ns("final_set_picker_ui")),
          uiOutput(ns("causal_refine_ui")),
          withSpinner(uiOutput(ns("final_summary_ui")), color = "#2c6fbb", type = 6),
          div(class = "table-toolbar", downloadButton(ns("final_download_table"), "Final set (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns("final_table"))
        )
      )
    })

    module_choices <- reactive({
      wg <- results$wgcna
      req(wg, length(wg$module_genes) > 0)
      mods <- setdiff(names(wg$module_genes), "grey")
      validate(need(length(mods) > 0, "No non-grey modules were detected in the WGCNA run."))
      list(mods = mods, sizes = lengths(wg$module_genes[mods]),
           disease = intersect(wg$significant_trait_modules %||% character(0), mods))
    })

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

    output$gene_panel_box_ui <- renderUI({
      req(!isTRUE(dataset$is_bundled_reference))
      box(
        width = 12, title = "Narrow further with a gene panel (optional)", status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "Off by default - the module/DEG overlap above is unchanged unless you pick a panel here. When set, it narrows the overlap to genes in a specific biological process or panel."),
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

    output$pooled_deg_picker_ui <- renderUI({
      ch <- tryCatch(deg_run_choices(), error = function(e) NULL)
      if (is.null(ch)) {
        return(div(class = "empty-note", icon("circle-info"),
          "Run Differential Expression first (Differential Expression tab)."))
      }
      selectInput(ns("pooled_deg_run"), "DEG contrast", choices = ch, selected = unname(utils::tail(ch, 1)), selectize = FALSE)
    })

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
    pooled_result <- sex_candidates("pooled_deg_run", "pooled_run_btn")

    female_has_run <- reactiveVal(FALSE)
    male_has_run <- reactiveVal(FALSE)
    pooled_has_run <- reactiveVal(FALSE)
    observeEvent(input$female_run_btn, female_has_run(TRUE), ignoreInit = TRUE)
    observeEvent(input$male_run_btn, male_has_run(TRUE), ignoreInit = TRUE)
    observeEvent(input$pooled_run_btn, pooled_has_run(TRUE), ignoreInit = TRUE)
    observeEvent(dataset$source, {
      female_has_run(FALSE); male_has_run(FALSE); pooled_has_run(FALSE)
    }, ignoreInit = TRUE)

    female_safe <- function() { req(female_has_run()); female_result() }
    male_safe <- function() { req(male_has_run()); male_result() }
    pooled_safe <- function() { req(pooled_has_run()); pooled_result() }

    register_panel <- function(prefix, res, btn_label = sprintf("Compute %s candidates", prefix)) {
      output[[paste0(prefix, "_summary_ui")]] <- renderUI({
        r <- tryCatch(res(), error = function(e) NULL)
        if (is.null(r)) {
          return(div(class = "empty-note", icon("circle-info"),
            sprintf("Not run yet. Click \"%s\" above.", btn_label)))
        }
        tagList(
          p(strong(r$n_bg), " module-background genes, ", strong(r$n_deg), " significant DEGs (", r$contrast, ")",
            if (!is.na(r$n_panel)) tagList(", ", strong(r$n_panel), " genes in the selected gene panel") else NULL, "."),
          p(strong(length(r$overlap)), " candidate genes in the overlap - hypergeometric enrichment ",
            HTML("<em>p</em>"), " = ", signif(r$p_value, 3), ".")
        )
      })

      venn_fill_high <- c(female = "#1a7a3c", male = "#7a4a26", pooled = "#2c6fbb")[[prefix]]

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

    register_panel("female", female_safe)
    register_panel("male", male_safe)
    register_panel("pooled", pooled_safe, btn_label = "Compute candidates")

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
      fr <- tryCatch(female_safe(), error = function(e) NULL)
      mr <- tryCatch(male_safe(), error = function(e) NULL)
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

    MR_P_CUT <- 0.05
    COLOC_PP4_CUT <- 0.8

    mr_supported_genes <- reactive({
      gt <- results$mr$genes_tested
      if (is.null(gt) || length(gt) == 0) return(character(0))
      ps <- vapply(gt, function(e) e$p %||% NA_real_, numeric(1))
      names(ps)[!is.na(ps) & ps < MR_P_CUT]
    })

    coloc_supported_genes <- reactive({
      gt <- results$coloc$genes_tested
      if (is.null(gt) || length(gt) == 0) return(character(0))
      pp <- vapply(gt, function(e) e$pp_h4 %||% NA_real_, numeric(1))
      names(pp)[!is.na(pp) & pp >= COLOC_PP4_CUT]
    })

    base_final <- reactive({
      if (isTRUE(sex_available())) {
        fc <- final_candidates()
        list(genes = fc$genes, stats = fc$stats)
      } else {
        pr <- pooled_safe()
        list(genes = pr$overlap, stats = pr$stats)
      }
    })

    output$causal_refine_ui <- renderUI({
      base <- tryCatch(base_final(), error = function(e) NULL)
      if (is.null(base) || length(base$genes) == 0) return(NULL)
      mr_hits <- intersect(base$genes, mr_supported_genes())
      coloc_hits <- intersect(base$genes, coloc_supported_genes())
      if (length(mr_hits) == 0 && length(coloc_hits) == 0) return(NULL)
      tagList(
        tags$hr(),
        p(strong("Optional: refine by causal genetic evidence."),
          " Mendelian Randomization and Colocalization are separate, optional analyses (Genetics section) - skip this if you don't have GWAS/eQTL data for your trait; everything downstream works fine on the set above alone."),
        if (length(mr_hits) > 0) checkboxInput(ns("require_mr"),
          sprintf("Require Mendelian Randomization support (p < %.2f) - %d of %d candidates qualify", MR_P_CUT, length(mr_hits), length(base$genes)), value = FALSE),
        if (length(coloc_hits) > 0) checkboxInput(ns("require_coloc"),
          sprintf("Require Colocalization support (PP.H4 ≥ %.1f) - %d of %d candidates qualify", COLOC_PP4_CUT, length(coloc_hits), length(base$genes)), value = FALSE)
      )
    })

    refined_final <- reactive({
      base <- base_final()
      genes <- base$genes
      if (isTRUE(input$require_mr)) genes <- intersect(genes, mr_supported_genes())
      if (isTRUE(input$require_coloc)) genes <- intersect(genes, coloc_supported_genes())
      validate(need(length(genes) > 0, "No candidate genes pass the selected causal-evidence filter(s) above - uncheck one to continue."))
      list(genes = genes, stats = base$stats[base$stats$gene %in% genes, , drop = FALSE])
    })

    output$final_summary_ui <- renderUI({
      fc <- refined_final()
      basis <- if (isTRUE(sex_available())) {
        switch(input$final_candidate_set %||% "union",
               female = "female only", male = "male only",
               union = "found in either sex", intersection = "shared by both sexes")
      } else "module ∩ DEG overlap"
      refined_note <- if (isTRUE(input$require_mr) || isTRUE(input$require_coloc)) ", after the causal-evidence filter(s) above" else ""
      p(strong(length(fc$genes)), " genes in the final candidate set (", basis, refined_note, ").")
    })

    output$final_table <- DT::renderDataTable({
      DT::datatable(refined_final()$stats, rownames = FALSE, filter = "top",
                     options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })

    output$final_download_table <- downloadHandler(
      filename = function() sprintf("final_candidate_genes_%s.csv", if (isTRUE(sex_available())) input$final_candidate_set %||% "union" else "pooled"),
      content = function(file) write.csv(refined_final()$stats, file, row.names = FALSE)
    )

    observe({
      pr <- prereqs()
      if (!isTRUE(pr$wgcna_ok) || !isTRUE(pr$dge_ok)) {
        results$candidates <- NULL
        return()
      }
      fc <- tryCatch(refined_final(), error = function(e) NULL)
      if (isTRUE(sex_available())) {
        fr <- tryCatch(female_safe(), error = function(e) NULL)
        mr <- tryCatch(male_safe(), error = function(e) NULL)
        if (is.null(fr) && is.null(mr)) return()
        results$candidates <- list(
          female = if (!is.null(fr)) list(n_candidates = length(fr$overlap), genes = fr$overlap, contrast = fr$contrast) else NULL,
          male = if (!is.null(mr)) list(n_candidates = length(mr$overlap), genes = mr$overlap, contrast = mr$contrast) else NULL,
          final = if (!is.null(fc)) list(selection = input$final_candidate_set, n_candidates = length(fc$genes), genes = fc$genes) else NULL
        )
      } else {
        pr <- tryCatch(pooled_safe(), error = function(e) NULL)
        if (is.null(pr)) return()
        results$candidates <- list(
          female = NULL, male = NULL,
          final = if (!is.null(fc)) list(selection = "pooled", n_candidates = length(fc$genes), genes = fc$genes, contrast = pr$contrast) else NULL
        )
      }
    })
  })
}
