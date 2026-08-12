## R/mod_dge.R
## Submodule: Differential Expression (Section 2.3)
## "Your analysis" fits a live model comparing two levels of a user-chosen
## metadata column in the currently loaded dataset (not just "group" - any
## column with 2-20 distinct values, e.g. group, sex, dataset, batch, tx),
## with an optional second column either used as a filter (restrict to one
## of its levels) or as a covariate the model adjusts for. Two methods are
## available: limma (continuous/log-scale data - microarray, or already-
## normalised RNA-seq) and DESeq2 (raw, non-negative integer RNA-seq
## counts) - "Auto-detect" picks between them by inspecting the loaded
## expression matrix's values.

mod_dge_config <- list(
  id = "dge", group = "Data",
  title = "Differential Expression",
  description = "Fit a live limma or DESeq2 model comparing two levels of any metadata column in the currently loaded dataset, with an optional second column to filter or adjust for.",
  icon = "chart-column"
)

mod_dge_ui <- function(id) {
  ns <- NS(id)
  tagList(
      fluidRow(
        column(
          4,
          arthochat_shortcut_ui(
            "New to differential expression? Ask ArthOChat.",
            compact = TRUE
          ),
          box(
            width = NULL, title = "Contrast", status = "primary", solidHeader = FALSE,
            uiOutput(ns("method_and_columns_ui")),
            uiOutput(ns("ref_comp_ui")),
            uiOutput(ns("covariate_controls_ui")),
            numericInput(ns("padj_cut"), "Adjusted p-value cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
            ## 0.1, not a rounder-looking 0.5: this project's own primary
            ## sex-stratified analysis used |log2FC| > 0.1 (a permissive
            ## 1.07-fold gate deliberately used as a coherence filter rather
            ## than a prioritisation step - see Chapter_2_subchapter2_
            ## sexstratified.md Section 2.3), not a stricter fold-change cut.
            numericInput(ns("lfc_cut"), "Absolute log2 fold-change cutoff", value = 0.1, min = 0, step = 0.1),
            actionButton(ns("run_btn"), "Run differential expression", icon = icon("play"), class = "btn-primary btn-sm"),
            div(style = "margin-top:10px;", uiOutput(ns("saved_runs_ui")))
          )
        ),
        column(
          8,
          box(
            width = NULL, title = "Result", status = "primary", solidHeader = FALSE,
            withSpinner(uiOutput(ns("summary_ui")), color = "#2563EB", type = 6),
            withSpinner(plotOutput(ns("volcano"), height = 460), color = "#2563EB", type = 6),
            div(class = "table-toolbar", downloadButton(ns("download_volcano_png"), "Download volcano plot (PNG, 7×6in @ 300dpi)", class = "btn-sm"))
          )
        )
      ),
      box(
        width = 12, title = "Heatmap of top differentially expressed genes", status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "Per-gene z-scored expression, clustered by gene and sample, for the most significant genes in this contrast."),
        fluidRow(
          column(4, numericInput(ns("heatmap_n"), "Top significant genes to show", value = 30, min = 2, max = 200, step = 5)),
          column(8, div(style = "padding-top: 25px;", downloadButton(ns("download_heatmap_png"), "Download heatmap (PNG, 300dpi)", class = "btn-sm")))
        ),
        withSpinner(plotOutput(ns("heatmap"), height = 520), color = "#2563EB", type = 6)
      ),
      box(
        width = 12, title = "Result table", status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "Every tested gene, with an up/down/not-significant “direction” column. Use the column filter below the header to show only Up or only Down genes."),
        div(class = "table-toolbar", downloadButton(ns("download_dge"), "Download CSV", class = "btn-sm")),
        DT::dataTableOutput(ns("dge_table"))
      ),
      uiOutput(ns("references_box_ui"))
  )
}

mod_dge_server <- function(id, dataset, results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Two pipelines, one module: the References box below cites "the
    ## papers this project's own methodology cites" - true for the default/
    ## preloaded cohort, not for an arbitrary dataset a user brought in
    ## themselves. dataset$source starts with "Uploaded dataset" for BOTH
    ## upload paths: mod_dataset.R's plain "Uploaded dataset: X + Y" and
    ## mod_preprocessing.R's activate button, "Uploaded dataset
    ## (preprocessed + batch-corrected): ..." - match the prefix without
    ## requiring an immediate colon so both are caught.
    output$references_box_ui <- renderUI({
      req(!grepl("^Uploaded dataset", dataset$source %||% ""))
      box(
        width = 12, title = "References", status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "Background reading for the methods this tab implements, and the papers this project's own methodology cites. Verified against live PubMed lookups, not from memory."),
        tags$ul(
          class = "dge-ref-list",
          tags$li(strong("limma (moderated t-test): "), "Ritchie ME et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. ",
                  tags$em("Nucleic Acids Research"), ". ", tags$a(href = "https://pubmed.ncbi.nlm.nih.gov/25605792", target = "_blank", "PMID: 25605792")),
          tags$li(strong("DESeq2 (negative binomial GLM): "), "Love MI, Huber W, Anders S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. ",
                  tags$em("Genome Biology"), ". ", tags$a(href = "https://pubmed.ncbi.nlm.nih.gov/25516281", target = "_blank", "PMID: 25516281")),
          tags$li(strong("Benjamini-Hochberg FDR (behind adj.P.Val/padj in both methods): "), "Benjamini Y, Hochberg Y. (1995). Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. ",
                  tags$em("Journal of the Royal Statistical Society: Series B"), ", 57(1), 289-300 (not in PubMed; a statistics journal)."),
          tags$li(strong("arrayWeights (used in this project's own limma pipeline): "), "Ritchie ME et al. (2006). Empirical array quality weights in the analysis of microarray data. ",
                  tags$em("BMC Bioinformatics"), ". ", tags$a(href = "https://pubmed.ncbi.nlm.nih.gov/16712727", target = "_blank", "PMID: 16712727")),
          tags$li(strong("New to interpreting DGE results? "), "McDermaid A et al. (2019). Interpretation of differential gene expression results of RNA-seq data: review and integration. ",
                  tags$em("Briefings in Bioinformatics"), ". ", tags$a(href = "https://pubmed.ncbi.nlm.nih.gov/30099484", target = "_blank", "PMID: 30099484"))
        ),
        p(class = "submodule-desc", "Have a more specific question? ", strong("Ask ArthOChat"), " above. It can search PubMed live for anything not covered here.")
      )
    })

    ## Any metadata column with 2-20 distinct non-missing values is a
    ## plausible contrast/covariate variable - not hardcoded to "group"/
    ## "sex", so this works on whatever columns the currently loaded
    ## dataset (default, preloaded, or the user's own upload) actually has.
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

    output$method_and_columns_ui <- renderUI({
      cols <- candidate_columns()
      validate(need(length(cols) > 0, "The loaded metadata has no column with 2-20 distinct values to contrast on."))
      default_contrast <- if ("group" %in% cols) "group" else cols[1]
      tagList(
        radioButtons(
          ns("method"), "Method",
          choiceNames = list(
            "Auto-detect (recommended)",
            "limma (microarray, log-scale, or already-normalised data)",
            "DESeq2 (RNA-seq raw counts)"
          ),
          choiceValues = list("auto", "limma", "deseq2"), selected = "auto"
        ),
        selectInput(ns("contrast_col"), "Contrast column", choices = cols, selected = default_contrast, selectize = FALSE),
        selectInput(ns("covariate_col"), "Second column (optional)", choices = c("(none)", cols), selected = "(none)", selectize = FALSE)
      )
    })

    ## Reference/comparison level choices depend on which column is picked
    ## above, so this is a separate output from method_and_columns_ui
    ## (which defines contrast_col itself) - same "depends on the previous
    ## select's value" split already used elsewhere in this app, e.g.
    ## mod_preprocessing.R's group_filter_ui/filter_val_ui.
    output$ref_comp_ui <- renderUI({
      req(input$contrast_col)
      meta <- dataset$meta
      req(input$contrast_col %in% colnames(meta))
      lvls <- sort(unique(na.omit(as.character(meta[[input$contrast_col]]))))
      validate(need(length(lvls) >= 2, "This column has fewer than two distinct values."))
      tagList(
        selectInput(ns("ref_group"), "Reference level", choices = lvls, selected = lvls[1], selectize = FALSE),
        selectInput(ns("comp_group"), "Comparison level", choices = lvls, selected = lvls[min(2, length(lvls))], selectize = FALSE)
      )
    })

    output$covariate_controls_ui <- renderUI({
      req(input$covariate_col)
      if (identical(input$covariate_col, "(none)")) return(NULL)
      tagList(
        radioButtons(
          ns("covariate_mode"), NULL,
          choiceNames = list("Filter to one level of this column", "Adjust for this column as a covariate"),
          choiceValues = list("filter", "adjust"), selected = "filter"
        ),
        uiOutput(ns("covariate_level_ui"))
      )
    })

    output$covariate_level_ui <- renderUI({
      req(input$covariate_col, input$covariate_mode)
      if (identical(input$covariate_col, "(none)") || !identical(input$covariate_mode, "filter")) return(NULL)
      meta <- dataset$meta
      req(input$covariate_col %in% colnames(meta))
      lvls <- sort(unique(na.omit(as.character(meta[[input$covariate_col]]))))
      selectInput(ns("covariate_level"), "Restrict to", choices = lvls, selected = lvls[1], selectize = FALSE)
    })

    ## Auto-detect heuristic: same non-negative / "99th percentile > 100
    ## means linear scale" rule already used throughout this app for scale
    ## detection (see mod_preprocessing.R's Scale check, and
    ## needs_quantile_norm() in global.R) - raw counts are non-negative and
    ## span a wide linear range; continuous/log-scale data (microarray, or
    ## RNA-seq already turned into log-CPM) doesn't. Deliberately NOT an
    ## "are these exact integers" check - real-world count matrices
    ## (RSEM-style estimated counts, e.g. 472.59) are often fractional.
    looks_like_raw_counts <- function(m) {
      vals <- as.numeric(m)
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0 || any(vals < 0)) return(FALSE)
      q99 <- suppressWarnings(stats::quantile(vals[vals > 0], 0.99, na.rm = TRUE))
      isTRUE(!is.na(q99) && q99 > 100)
    }

    fit_result <- eventReactive(input$run_btn, {
      req(input$contrast_col, input$ref_group, input$comp_group)
      validate(need(input$ref_group != input$comp_group, "Reference and comparison level must be different."))

      contrast_col <- input$contrast_col
      meta <- dataset$meta
      meta <- meta[!is.na(meta[[contrast_col]]) & as.character(meta[[contrast_col]]) %in% c(input$ref_group, input$comp_group), , drop = FALSE]

      use_covariate <- !is.null(input$covariate_col) && !identical(input$covariate_col, "(none)")
      covariate_mode <- input$covariate_mode %||% "filter"
      covariate_label <- ""
      if (use_covariate) {
        meta <- meta[!is.na(meta[[input$covariate_col]]), , drop = FALSE]
        if (identical(covariate_mode, "filter")) {
          req(input$covariate_level)
          meta <- meta[as.character(meta[[input$covariate_col]]) == input$covariate_level, , drop = FALSE]
          covariate_label <- sprintf(" (%s = %s)", input$covariate_col, input$covariate_level)
        }
      }

      common <- intersect(colnames(dataset$expr), meta$sample)
      validate(need(length(common) >= 6, "Fewer than 6 samples match this contrast (and filter, if set); pick a different combination."))
      meta <- meta[match(common, meta$sample), , drop = FALSE]
      expr <- dataset$expr[, common, drop = FALSE]

      grp <- factor(as.character(meta[[contrast_col]]), levels = c(input$ref_group, input$comp_group))
      validate(need(all(table(grp) >= 3), "Each level needs at least 3 samples in this contrast to fit a model."))

      adjust_for_covariate <- use_covariate && identical(covariate_mode, "adjust")
      covar <- NULL
      if (adjust_for_covariate) {
        covar <- factor(as.character(meta[[input$covariate_col]]))
        validate(need(length(unique(covar)) >= 2, "The covariate column has only one level left in this contrast. Pick \"Filter\" instead, or a different column."))
        covariate_label <- sprintf(" adjusted for %s", input$covariate_col)
      }

      method <- input$method %||% "auto"
      is_counts <- looks_like_raw_counts(expr)
      used_method <- if (identical(method, "auto")) (if (is_counts) "deseq2" else "limma") else method
      if (identical(used_method, "deseq2")) {
        validate(need(is_counts, "DESeq2 needs raw, non-negative integer counts, but this data has negative or non-integer values (it looks already normalised/log-scale). Pick limma instead, or load raw counts directly via Dataset → Upload your own data (Preprocessing → Batch Correction always outputs normalised, log-scale data, even with Preprocessing's own log2 set to \"Skip\")."))
      }

      tt <- if (identical(used_method, "deseq2")) {
        counts <- round(as.matrix(expr))
        storage.mode(counts) <- "integer"
        col_data <- data.frame(grp = grp, row.names = colnames(counts))
        design <- ~grp
        if (adjust_for_covariate) { col_data$covar <- covar; design <- ~covar + grp }
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = design)
        dds <- tryCatch(suppressMessages(DESeq2::DESeq(dds)),
                         error = function(e) validate(need(FALSE, paste("DESeq2 could not fit this model:", conditionMessage(e)))))
        res <- DESeq2::results(dds, contrast = c("grp", levels(grp)[2], levels(grp)[1]))
        out <- as.data.frame(res)
        out$gene <- rownames(out)
        rownames(out) <- NULL
        out <- dplyr::rename(out, logFC = log2FoldChange, P.Value = pvalue, adj.P.Val = padj)
        out <- out[!is.na(out$P.Value), , drop = FALSE]
        out[order(out$P.Value), c("gene", setdiff(colnames(out), "gene"))]
      } else {
        design <- tryCatch(
          if (adjust_for_covariate) model.matrix(~0 + grp + covar) else model.matrix(~0 + grp),
          error = function(e) validate(need(FALSE, paste("Could not build a design matrix for this contrast/covariate combination:", conditionMessage(e))))
        )
        colnames(design)[seq_len(nlevels(grp))] <- levels(grp)
        ## Per-array quality weights (Ritchie et al. 2006, cited in
        ## References below) - this project's own methodology
        ## (Chapter_2_subchapter2_sexstratified.md, Section 2.2) applies
        ## these because the training cohort combines two different
        ## microarray platforms (HG-U219/GPL570 and HG-U133_Plus_2/
        ## GPL13667): a noisy or poor-quality array gets an estimated
        ## weight below 1 across all ~15,000 genes, so lmFit gives it
        ## proportionally less influence on the fit instead of every array
        ## counting equally regardless of quality. The References box below
        ## already cited arrayWeights as used in this pipeline; this is
        ## that citation actually being applied, not just referenced.
        aw <- limma::arrayWeights(expr, design)
        fit <- limma::lmFit(expr, design, weights = aw)
        cm <- limma::makeContrasts(contrasts = paste0(levels(grp)[2], "-", levels(grp)[1]), levels = design)
        fit2 <- tryCatch(limma::eBayes(limma::contrasts.fit(fit, cm)),
                          error = function(e) validate(need(FALSE, paste("limma could not fit this contrast:", conditionMessage(e)))))
        out <- limma::topTable(fit2, number = Inf, sort.by = "P")
        out$gene <- rownames(out)
        rownames(out) <- NULL
        out[, c("gene", setdiff(colnames(out), "gene"))]
      }

      ## Shown verbatim in summary_ui below, for reproducibility - the exact
      ## design formula and test this specific run used, not a generic
      ## textbook description.
      design_formula <- if (identical(used_method, "deseq2")) {
        if (adjust_for_covariate) sprintf("~ %s + %s", input$covariate_col, contrast_col) else sprintf("~ %s", contrast_col)
      } else {
        if (adjust_for_covariate) sprintf("~0 + %s + %s", contrast_col, input$covariate_col) else sprintf("~0 + %s", contrast_col)
      }
      test_label <- switch(used_method,
        limma = "moderated t-test (limma eBayes, array-quality-weighted)",
        deseq2 = "Wald test (DESeq2)",
        used_method
      )

      list(
        table = tt, method = used_method,
        ## expr/grp are the exact sample subset and group assignment this
        ## contrast was fit on - kept here so the heatmap below can reuse
        ## them directly instead of re-deriving the same filtering logic.
        expr = expr, grp = grp,
        design_formula = design_formula, test_label = test_label,
        n_ref = sum(grp == levels(grp)[1]), n_comp = sum(grp == levels(grp)[2]),
        label = paste0(input$comp_group, " vs ", input$ref_group, " (", contrast_col, ")", covariate_label)
      )
    }, ignoreInit = TRUE)

    dge_has_run <- reactiveVal(FALSE)
    observeEvent(input$run_btn, dge_has_run(TRUE), ignoreInit = TRUE)

    sig_table <- reactive({
      res <- fit_result()
      req(res)
      res$table %>%
        mutate(
          significant = !is.na(adj.P.Val) & adj.P.Val < input$padj_cut & abs(logFC) > input$lfc_cut,
          direction = case_when(
            !significant ~ "Not significant",
            logFC > 0 ~ "Up",
            TRUE ~ "Down"
          )
        )
    })

    ## Triggered on the same button fit_result() itself listens to, and
    ## silently does nothing if the fit failed validation - not on
    ## sig_table()/fit_result() directly, which would re-throw that same
    ## validation error into this observer too (see the tryCatch()+req()
    ## note above output$volcano below for why that matters).
    observeEvent(input$run_btn, {
      df <- tryCatch(sig_table(), error = function(e) NULL)
      req(df)
      res <- fit_result()
      results$dge <- list(
        contrast = res$label,
        method = res$method,
        n_samples = res$n_ref + res$n_comp,
        n_tested = nrow(df),
        n_significant = sum(df$significant),
        n_up = sum(df$direction == "Up"),
        n_down = sum(df$direction == "Down"),
        top_hits = head(df$gene[order(df$adj.P.Val)], 10)
      )

      ## Appended (not overwritten) so the Candidate Gene Identification tab
      ## can offer a click-based pick from EVERY contrast run this session -
      ## e.g. a female-vs-control run and a later male-vs-control run both
      ## stay available to compare, instead of the second silently replacing
      ## the first the way results$dge above does. Capped at the 8 most
      ## recent runs so a long session doesn't grow this without bound; a
      ## run is unlikely to still be wanted many contrasts later.
      runs <- results$dge_runs %||% list()
      run_id <- paste0("run", length(runs) + 1L)
      runs[[run_id]] <- list(
        contrast = res$label, method = res$method,
        n_samples = res$n_ref + res$n_comp,
        table = df[, c("gene", "logFC", "adj.P.Val", "direction")]
      )
      if (length(runs) > 8) runs <- utils::tail(runs, 8)
      results$dge_runs <- runs

      ## Saving into results$dge_runs above is the ONLY thing that makes a
      ## run pickable later (e.g. Candidate Gene Identification's "Female
      ## DEG contrast"/"Male DEG contrast" pickers, which read this exact
      ## list) - but it happened silently in a background observer with
      ## nothing on screen to show it. This notification plus
      ## saved_runs_ui below are that visible confirmation.
      showNotification(
        sprintf("Saved \"%s\" as run %d of %d this session - pick it up in Candidate Gene Identification's DEG contrast pickers.",
                res$label, length(runs), length(runs)),
        type = "message", duration = 6
      )
    })

    ## Female DEG / Male DEG completion status, directly under the Run
    ## button - by default this app expects the sex-stratified workflow
    ## (run once filtered/adjusted to female, once to male; see
    ## mod_candidates.R's header comment), but nothing on this tab ever
    ## showed whether that had actually happened. Every run saved into
    ## results$dge_runs above is silent otherwise, so without this a user
    ## who's run both female and male contrasts has no on-screen sign of
    ## it here - only in Candidate Gene Identification, a different tab.
    ## Detection mirrors mod_candidates.R's guess_run() exactly (same
    ## \\bfemale\\b|\\bF\\b / \\bmale\\b|\\bM\\b pattern against the
    ## contrast label, most recent match wins) so "completed" here means
    ## precisely what Candidate Gene Identification will actually pick up
    ## - not a separate, potentially-inconsistent guess.
    latest_run_matching <- function(pattern) {
      runs <- results$dge_runs %||% list()
      if (length(runs) == 0) return(NULL)
      labels <- vapply(runs, function(r) r$contrast, character(1))
      hit <- which(grepl(pattern, labels, ignore.case = TRUE, perl = TRUE))
      if (length(hit) == 0) return(NULL)
      labels[[utils::tail(hit, 1)]]
    }

    output$saved_runs_ui <- renderUI({
      female_label <- latest_run_matching("\\bfemale\\b|\\bF\\b")
      male_label <- latest_run_matching("\\bmale\\b|\\bM\\b")

      status_row <- function(sex, label) {
        if (is.null(label)) {
          tags$li(icon("circle-minus", style = "color: #8A929C;"), sprintf(" %s DEG - not run yet", sex))
        } else {
          tags$li(icon("check", style = "color: #1a9c5f;"), strong(sprintf(" %s DEG completed: ", sex)), label)
        }
      }

      tagList(
        p(class = "submodule-desc", style = "margin-bottom: 4px;",
          "Sex-stratified status (used by Candidate Gene Identification):"),
        tags$ul(style = "padding-left: 18px; margin-bottom: 0; list-style: none;",
                status_row("Female", female_label), status_row("Male", male_label))
      )
    })

    output$summary_ui <- renderUI({
      if (!dge_has_run()) {
        return(div(class = "empty-note", icon("circle-info"),
                    "Not run yet. Set a contrast on the left, then click \"Run differential expression\"."))
      }
      res <- fit_result()
      df <- sig_table()
      n_sig <- sum(df$significant)
      n_up <- sum(df$direction == "Up")
      n_down <- sum(df$direction == "Down")
      method_label <- switch(res$method, limma = "limma", deseq2 = "DESeq2", res$method)
      tagList(
        p(strong("Contrast: "), res$label, sprintf(" (n = %d vs %d, %s)", res$n_comp, res$n_ref, method_label)),
        p(strong("Model: "), tags$code(res$design_formula), sprintf(", %s.", res$test_label)),
        p(strong(format(nrow(df), big.mark = ",")), " genes tested, ",
          strong(format(n_sig, big.mark = ",")), " significant at the current cutoffs."),
        div(class = "pipeline-status-strip",
            span(class = "badge-up", sprintf("%s upregulated", format(n_up, big.mark = ","))),
            span(class = "badge-down", sprintf("%s downregulated", format(n_down, big.mark = ","))))
      )
    })

    ## Publication-style volcano: theme_classic() base, bold axis titles with
    ## proper log subscript notation, up/down counts annotated directly on
    ## the plot (common journal convention), gene labels for the top hits,
    ## and a horizontal reference line at the largest raw p-value among
    ## currently-significant points - the actual data-driven boundary the
    ## adjusted-p cutoff maps to here, rather than a fixed but misleading
    ## transform of padj_cut plotted against the raw-p y-axis. Built once as
    ## its own reactive so the on-screen plot and the PNG download always
    ## show exactly the same figure (same pattern as the Venn diagrams in
    ## mod_preprocessing.R).
    volcano_plot_obj <- reactive({
      df <- sig_table()
      top_labels <- df %>% dplyr::filter(significant) %>% dplyr::arrange(adj.P.Val) %>% head(15)
      n_up <- sum(df$direction == "Up")
      n_down <- sum(df$direction == "Down")
      sig_p <- suppressWarnings(max(df$P.Value[df$significant], na.rm = TRUE))
      xr <- max(abs(df$logFC), na.rm = TRUE) * 1.08

      p <- ggplot(df, aes(x = logFC, y = -log10(P.Value), color = direction, size = direction)) +
        geom_point(alpha = 0.7) +
        scale_color_manual(
          values = c(Up = ARTHOMIX_COLORS$red, Down = ARTHOMIX_COLORS$blue, `Not significant` = "#B8BEC8"),
          name = NULL
        ) +
        scale_size_manual(values = c(Up = 1.7, Down = 1.7, `Not significant` = 1.0), guide = "none") +
        geom_vline(xintercept = c(-input$lfc_cut, input$lfc_cut), linetype = "dashed", linewidth = 0.4, color = "#8A929C") +
        scale_x_continuous(limits = c(-xr, xr)) +
        labs(x = expression(log[2] ~ "fold-change"), y = expression(-log[10] ~ "(" * italic(P) * "-value)")) +
        theme_classic(base_size = 14) +
        theme(
          axis.title = element_text(face = "bold", color = "#1A1A1A"),
          axis.text = element_text(color = "#1A1A1A"),
          axis.line = element_line(color = "#1A1A1A", linewidth = 0.4),
          axis.ticks = element_line(color = "#1A1A1A", linewidth = 0.4),
          legend.position = "top", legend.text = element_text(size = 11),
          plot.margin = ggplot2::margin(16, 18, 10, 10)
        )
      if (is.finite(sig_p)) {
        p <- p + geom_hline(yintercept = -log10(sig_p), linetype = "dashed", linewidth = 0.4, color = "#8A929C")
      }
      p +
        ggrepel::geom_text_repel(data = top_labels, mapping = aes(x = logFC, y = -log10(P.Value), label = gene),
                                  inherit.aes = FALSE, size = 3.2, fontface = "italic", color = "#1A1A1A",
                                  segment.size = 0.3, max.overlaps = 20, show.legend = FALSE) +
        annotate("text", x = -xr, y = Inf, label = sprintf("%s down", format(n_down, big.mark = ",")),
                 hjust = 0, vjust = 1.6, size = 4.2, fontface = "bold", color = ARTHOMIX_COLORS$blue) +
        annotate("text", x = xr, y = Inf, label = sprintf("%s up", format(n_up, big.mark = ",")),
                 hjust = 1, vjust = 1.6, size = 4.2, fontface = "bold", color = ARTHOMIX_COLORS$red)
    })

    ## summary_ui below is the one place a failed fit_result() validation
    ## message is shown - every other output that depends on the same
    ## fit_result()/sig_table() chain (volcano, heatmap, dge_table) would
    ## otherwise redisplay that identical message in its own box, since
    ## each output independently re-evaluates (and re-throws) the same
    ## cached validation error. tryCatch() + req() here swallows it and
    ## leaves this panel blank instead of repeating the text.
    output$volcano <- renderPlot({
      p <- tryCatch(volcano_plot_obj(), error = function(e) NULL)
      req(p)
      p
    })

    ## width/height chosen to match a standard single-to-1.5-column journal
    ## figure (~7in @ 300dpi is comfortably above the ~85-183mm range most
    ## journals specify for a single panel) - deliberately fixed regardless
    ## of the on-screen preview size, since the download is for actual use
    ## in a manuscript/poster, not just a bigger version of the preview.
    output$download_volcano_png <- downloadHandler(
      filename = function() "volcano_plot.png",
      content = function(file) ggsave(file, plot = volcano_plot_obj(), width = 7, height = 6, dpi = 300, bg = "white")
    )

    ## ---------------------------------------------------------------------
    ## Heatmap of the top significant genes from this contrast, z-scored
    ## per gene across the exact same sample subset fit_result() modeled,
    ## with a top color bar for the contrast group. DESeq2's raw counts are
    ## log2(x+1)-transformed first (limma's expr is already continuous/log
    ## scale). pheatmap (already installed) gives clustering + an
    ## annotation bar in one call - the standard tool for this figure type.
    ## ---------------------------------------------------------------------

    top_de_genes <- reactive({
      df <- sig_table()
      n <- input$heatmap_n %||% 30
      df %>% dplyr::filter(significant) %>% dplyr::arrange(adj.P.Val) %>% head(n) %>% dplyr::pull(gene)
    })

    heatmap_matrix <- reactive({
      res <- fit_result()
      genes <- top_de_genes()
      validate(need(length(genes) >= 2, "Fewer than 2 significant genes at the current cutoffs to draw a heatmap. Lower the cutoffs, or pick a contrast with more signal."))
      m <- res$expr[genes, , drop = FALSE]
      if (identical(res$method, "deseq2")) m <- log2(m + 1)
      m
    })

    heatmap_args <- reactive({
      res <- fit_result()
      mat <- heatmap_matrix()
      mat_z <- t(scale(t(mat)))
      mat_z[!is.finite(mat_z)] <- 0
      ann_col <- data.frame(Group = res$grp, row.names = colnames(mat_z))
      grp_colors <- setNames(c(ARTHOMIX_COLORS$blue, ARTHOMIX_COLORS$red), levels(res$grp))
      list(
        mat = mat_z, annotation_col = ann_col, annotation_colors = list(Group = grp_colors),
        color = grDevices::colorRampPalette(c(ARTHOMIX_COLORS$blue, "#FFFFFF", ARTHOMIX_COLORS$red))(100),
        show_rownames = nrow(mat_z) <= 60, show_colnames = FALSE,
        clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
        clustering_method = "complete", fontsize = 10, fontsize_row = 7, border_color = NA
      )
    })

    output$heatmap <- renderPlot({
      a <- tryCatch(heatmap_args(), error = function(e) NULL)
      req(a)
      ph <- pheatmap::pheatmap(
        a$mat, annotation_col = a$annotation_col, annotation_colors = a$annotation_colors,
        color = a$color, show_rownames = a$show_rownames, show_colnames = a$show_colnames,
        clustering_distance_rows = a$clustering_distance_rows, clustering_distance_cols = a$clustering_distance_cols,
        clustering_method = a$clustering_method, fontsize = a$fontsize, fontsize_row = a$fontsize_row,
        border_color = a$border_color, silent = TRUE
      )
      grid::grid.newpage()
      grid::grid.draw(ph$gtable)
    })

    output$download_heatmap_png <- downloadHandler(
      filename = function() "dge_heatmap.png",
      content = function(file) {
        a <- heatmap_args()
        pheatmap::pheatmap(
          a$mat, annotation_col = a$annotation_col, annotation_colors = a$annotation_colors,
          color = a$color, show_rownames = a$show_rownames, show_colnames = a$show_colnames,
          clustering_distance_rows = a$clustering_distance_rows, clustering_distance_cols = a$clustering_distance_cols,
          clustering_method = a$clustering_method, fontsize = a$fontsize, fontsize_row = a$fontsize_row,
          border_color = a$border_color,
          filename = file, width = 7, height = max(5, 0.15 * nrow(a$mat) + 3), res = 300
        )
      }
    )

    output$dge_table <- DT::renderDataTable({
      df <- tryCatch(sig_table(), error = function(e) NULL)
      req(df)
      DT::datatable(df, rownames = FALSE, filter = "top",
                     options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_dge <- downloadHandler(
      filename = function() "differential_expression.csv",
      content = function(file) write.csv(sig_table(), file, row.names = FALSE)
    )
  })
}
