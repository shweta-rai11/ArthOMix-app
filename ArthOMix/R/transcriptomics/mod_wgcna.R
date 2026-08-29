## R/mod_wgcna.R
## WGCNA co-expression network submodule (Section 2.4): a 6-step wizard
## (Filter & QC, Soft Power, Modules, Module-Trait, Hub Genes, Enrichment)
## following the standard WGCNA workflow (Langfelder & Horvath 2008;
## Zhang & Horvath 2005), parameterized similarly to xOmicsShiny's WGCNA
## module but reimplemented against this app's own data/UI plumbing and
## extended with outlier QC, a custom gene panel, scale-free topology check,
## module enrichment, Cytoscape export and STRING-DB lookup.

mod_wgcna_config <- list(
  id = "wgcna", group = "Network",
  title = "WGCNA Co-expression Network", section = "Section 2.4",
  description = "Co-expression network analysis performed by gene or sample filtering, soft-threshold power, module detection, module-trait correlation, hub genes, functional enrichment and Cytoscape export.",
  icon = "circle-nodes"
)

## Recodes a non-numeric trait column to numeric (factor levels sorted
## alphabetically); levels_keep optionally restricts to specific levels,
## setting the rest to NA so downstream cor(..., use = "p") ignores them.
wgcna_encode_trait <- function(meta, col, levels_keep = NULL) {
  v <- meta[[col]]
  if (is.numeric(v)) return(as.numeric(v))
  v <- as.character(v)
  if (!is.null(levels_keep) && length(levels_keep) > 0) v[!(v %in% levels_keep)] <- NA
  as.numeric(factor(v, levels = sort(unique(stats::na.omit(v)))))
}

## Widget-id suffix for a trait column's "which levels?" checkbox group, built
## from the column's position among eligible_cols (this module's own
## eligible_traits(), the same full/stable list every time for a fixed active
## dataset) rather than the raw column name - GEO's own raw pData column
## names routinely contain spaces/colons (e.g. "disease state:ch1", confirmed
## live off a real GSE93272 fetch), and Shiny's client dispatches messages
## keyed as "inputType:inputId" - a colon (or a space, which breaks the
## jQuery selector/DOM id) inside the id corrupts that routing and throws an
## uncaught "No handler registered for type ..." error that kills the whole
## session (confirmed live: the entire page goes grey/unresponsive - a Shiny
## disconnection, not a rendering glitch). Same fix as mod_overview.R's
## filter_id() for the identical risk there.
wgcna_trait_widget_id <- function(col, eligible_cols) paste0("trait_levels_", match(col, eligible_cols))

## Reads the Step 4 per-trait level-filter checkbox; NULL for numeric traits.
wgcna_trait_levels <- function(input, meta, col, eligible_cols) {
  if (is.numeric(meta[[col]])) return(NULL)
  input[[wgcna_trait_widget_id(col, eligible_cols)]]
}

## Picks the correlation function (bicor vs Pearson) matching Step 2's choice.
wgcna_cor_fnc <- function(cor_method) {
  if (identical(cor_method, "bicor")) WGCNA::bicor else WGCNA::cor
}
wgcna_cor_fnc_name <- function(cor_method) if (identical(cor_method, "bicor")) "bicor" else "cor"

## Builds a STRING-DB multi-gene network lookup URL (species fixed to human).
wgcna_string_url <- function(genes, species = 9606) {
  ids <- paste(vapply(genes, utils::URLencode, character(1), reserved = TRUE), collapse = "%0d")
  sprintf("https://string-db.org/cgi/network?identifiers=%s&species=%d", ids, species)
}

## Builds the STRING-DB REST image URL that renders a PNG network for a gene list.
wgcna_string_image_url <- function(genes, species = 9606) {
  ids <- paste(vapply(genes, utils::URLencode, character(1), reserved = TRUE), collapse = "%0d")
  sprintf("https://string-db.org/api/image/network?identifiers=%s&species=%d", ids, species)
}

## Loads the precomputed WGCNA result (scripts/00_shared/06_WGCNA.R output)
## from data/processed/wgcna_results.rds, plus the raw blockwiseModules
## object for the dendrogram if a wgcna_net_*.rds file is present.
load_precomputed_wgcna_result <- function() {
  proc_dir <- PROCESSED_DIR
  results_path <- file.path(proc_dir, "wgcna_results.rds")
  validate(need(file.exists(results_path),
                sprintf("No precomputed result found at %s.", results_path)))
  res <- readRDS(results_path)

  texpr <- res$datExpr
  module_colors <- res$moduleColors
  gene_module <- data.frame(gene = colnames(texpr), module = module_colors, stringsAsFactors = FALSE)
  gene_module <- gene_module[order(gene_module$module), ]

  net_files <- sort(list.files(proc_dir, pattern = "^wgcna_net_.*\\.rds$", full.names = TRUE), decreasing = TRUE)
  dendro <- NULL
  block_colors <- module_colors
  if (length(net_files) > 0) {
    net <- tryCatch(readRDS(net_files[1]), error = function(e) NULL)
    if (!is.null(net) && length(net$dendrograms) > 0) {
      dendro <- net$dendrograms[[1]]
      block_colors <- module_colors[net$blockGenes[[1]]]
    }
  }

  list(
    texpr = texpr, meta = res$meta, power = res$soft_power,
    network_type = res$config$network_type %||% "signed",
    tom_type = res$config$tom_type %||% "signed",
    cor_method = res$config$cor_type %||% "pearson",
    module_colors = module_colors, MEs = res$MEs,
    dendro = dendro, block_colors = block_colors,
    gene_module = gene_module,
    module_sizes = as.data.frame(table(module = module_colors)),
    n_genes = ncol(texpr), n_samples = nrow(texpr),
    n_modules = length(unique(module_colors[module_colors != "grey"])),
    ## Precomputed hub-gene table (disease modules, group trait); Step 5
    ## falls back to live computation for any other module/trait selection.
    hub_table_precomputed = res$hub_table
  )
}

## Step 2's counterpart to load_precomputed_wgcna_result() above: reads the same
## data/processed/wgcna_results.rds for just the soft-power diagnostics (the
## pickSoftThreshold fit table and the power actually used), so the app's own
## bundled default dataset doesn't have to re-run pickSoftThreshold/softConnectivity
## over its full ~15.7k-gene matrix live (a multi-minute computation - see the
## comment on get_or_compute_wgcna_blocks() in global.R). connectivity is left NULL
## since softConnectivity()'s raw per-gene output isn't saved by the offline script;
## the scale-free check panel shows an explanatory message instead of that one plot.
load_precomputed_wgcna_sft <- function() {
  proc_dir <- PROCESSED_DIR
  results_path <- file.path(proc_dir, "wgcna_results.rds")
  validate(need(file.exists(results_path),
                sprintf("No precomputed result found at %s.", results_path)))
  res <- readRDS(results_path)
  list(sft_df = res$sft_fit, power = res$soft_power, connectivity = NULL)
}

## Icon+label+sublabel markup for a workflow-stepper tab (no step number, unlike mod_preprocessing.R's pp_step_title()).
wgcna_step_title <- function(ic, label, sublabel) {
  tagList(
    icon(ic),
    span(class = "step-text",
         span(class = "step-label", label),
         span(class = "step-sublabel", sublabel))
  )
}

## Boxless results section: heading + divider instead of a bordered box() card (settings panels keep box()).
wgcna_result <- function(icon_name, title, ..., desc = NULL) {
  tagList(
    div(class = "wgcna-result-heading", icon(icon_name), title),
    if (!is.null(desc)) p(class = "submodule-desc wgcna-result-desc", desc),
    ...
  )
}

## Wraps a fluidRow of result columns with a divider class instead of per-column box borders.
wgcna_result_row <- function(...) {
  div(class = "wgcna-result-row", fluidRow(...))
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------

## Scoped CSS: equal-height boxes per row, flex-based nav-tabs (Bootstrap 3's
## default float layout reflows unpredictably for 6 icon+label tabs), and the
## boxless results-section styling above. Scoped to .wgcna-wrap only.
wgcna_layout_css <- function() {
  tags$style(HTML("
    .wgcna-wrap .row { display: flex; flex-wrap: wrap; }
    .wgcna-wrap .row > div[class*='col-'] { display: flex; flex-direction: column; }
    .wgcna-wrap .row > div[class*='col-'] > .box { flex: 1 1 auto; display: flex; flex-direction: column; margin-bottom: 0; }
    .wgcna-wrap .row > div[class*='col-'] > .box > .box-body { flex: 1 1 auto; }

    .wgcna-wrap .nav-tabs {
      display: flex !important;
      flex-wrap: wrap;
      border-bottom: 1px solid var(--color-border);
    }
    .wgcna-wrap .nav-tabs > li:not(:last-child)::after { display: none; }
    .wgcna-wrap .nav-tabs > li { min-width: 140px; flex: 1 1 30%; float: none; }
    .wgcna-wrap .nav-tabs > li > a {
      display: flex !important;
      align-items: center !important;
      border: none !important;
      background: transparent !important;
      border-radius: 0 !important;
      border-bottom: 2px solid transparent !important;
      padding: 7px 8px !important;
      gap: 6px !important;
      min-width: 0;
      color: var(--color-ink-secondary);
      outline: none;
    }
    .wgcna-wrap .nav-tabs > li > a:hover { color: var(--color-ink); background: transparent !important; }
    .wgcna-wrap .nav-tabs > li.active > a,
    .wgcna-wrap .nav-tabs > li.active > a:hover,
    .wgcna-wrap .nav-tabs > li.active > a:focus {
      color: var(--color-primary);
      background: transparent !important;
      border-bottom: 2px solid var(--color-primary) !important;
    }
    .wgcna-wrap .step-text { min-width: 0; display: flex; flex-direction: column; }
    /* wrap only at word boundaries, never mid-word */
    .wgcna-wrap .step-label, .wgcna-wrap .step-sublabel { white-space: normal; overflow-wrap: normal; word-break: normal; }
    .wgcna-wrap .step-label { font-size: 10.5px; line-height: 1.25; color: inherit; }
    .wgcna-wrap .nav-tabs > li.active .step-label { font-weight: 600; }
    .wgcna-wrap .step-sublabel { font-size: 9px; line-height: 1.25; color: var(--color-ink-secondary); }

    /* ---- Results section (boxless): heading + dividers instead of box() cards ---- */
    .wgcna-wrap .wgcna-result-heading {
      display: flex; align-items: center; gap: 8px;
      font-size: 15px; font-weight: 600; color: var(--color-ink);
      margin: 28px 0 6px 0; padding-top: 20px;
      border-top: 1px solid var(--color-border);
      clear: both; /* clears the floated box() wrapper above it */
    }
    .wgcna-wrap .box + .wgcna-result-heading,
    .wgcna-wrap .box + div > .wgcna-result-heading:first-child {
      margin-top: 20px;
    }
    .wgcna-wrap .wgcna-result-heading svg, .wgcna-wrap .wgcna-result-heading .fa { color: var(--color-primary); }
    .wgcna-wrap .wgcna-result-desc { margin-bottom: 12px; }
    .wgcna-wrap .wgcna-result-row { margin-bottom: 4px; clear: both; }
    .wgcna-wrap .wgcna-result-row .row { row-gap: 20px; }
    @media (min-width: 900px) {
      .wgcna-wrap .wgcna-result-row .row > div[class*='col-']:not(:first-child) {
        border-left: 1px solid var(--color-border);
      }
    }
    .wgcna-wrap .wgcna-result-subtitle {
      font-size: 13px; font-weight: 600; color: var(--color-ink);
      margin-bottom: 8px;
    }
  "))
}

## Left rail (ArthOChat shortcut) in column(3) beside the 6-step wizard in column(9).
mod_wgcna_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      3,
      div(
        class = "pipeline-rail",
        uiOutput(ns("pipeline_summary"))
      )
    ),
    column(
      9,
      div(
        class = "workflow-stepper-wrap wgcna-wrap",
        wgcna_layout_css(),
        tabsetPanel(
          id = ns("tabs"), type = "tabs",
          header = tagList(
            p(class = "submodule-desc", style = "margin-bottom: 8px;",
              "Run WGCNA on the currently loaded dataset: filter genes and samples, pick a network power, detect co-expression modules, relate them to your traits, then pull out and annotate hub genes."),
            tags$hr(style = "margin: 10px 0;")
          ),
          tabPanel(
            value = "Filter", title = wgcna_step_title("filter", "Filter & QC", "Genes & samples"),
            br(), uiOutput(ns("step1_ui"))
          ),
          tabPanel(
            value = "Power", title = wgcna_step_title("wave-square", "Soft Power", "Network power"),
            br(), uiOutput(ns("step2_controls_ui")), uiOutput(ns("step2_ui"))
          ),
          tabPanel(
            value = "Modules", title = wgcna_step_title("diagram-project", "Modules", "Detect & cluster"),
            br(), uiOutput(ns("step3_ui"))
          ),
          tabPanel(
            value = "Traits", title = wgcna_step_title("table-cells", "Module-Trait", "Correlate traits"),
            br(), uiOutput(ns("step4_ui"))
          ),
          tabPanel(
            value = "Hubs", title = wgcna_step_title("star", "Hub Genes", "Rank & export"),
            br(), uiOutput(ns("step5_ui"))
          ),
          tabPanel(
            value = "Enrichment", title = wgcna_step_title("flask", "Enrichment", "GO & KEGG"),
            br(), uiOutput(ns("step6_controls_ui")), uiOutput(ns("step6_ui"))
          )
        )
      )
    )
  )
}

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------

mod_wgcna_server <- function(id, dataset, results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## True whenever the active dataset is NOT the exact default merged training
    ## cohort - i.e. uploaded, GEO-fetched, or an individual raw preloaded GSE
    ## picked on its own; drives paper-methodology defaults in steps 1-3.
    ## Despite the name (kept to avoid touching every call site below),
    ## this is NOT "uploaded specifically" - it's dataset$is_bundled_reference
    ## inverted, the same flag the rest of the app uses to gate the project's
    ## own bundled/precomputed results. Originally checked only
    ## source_type=="uploaded"/a "^Uploaded dataset" source-string prefix,
    ## which meant a GEO-fetched dataset (e.g. a completely unrelated disease
    ## fetched live from NCBI) still got shown "12, this project's own
    ## choice" - a soft-power value calibrated specifically for the default
    ## merged whole-blood RA cohort - instead of the generic "uploaded data"
    ## defaults meant for exactly this situation. Confirmed live off a real
    ## GSE21942 (multiple sclerosis) fetch.
    dataset_is_uploaded <- reactive({
      !isTRUE(dataset$is_bundled_reference)
    })

    ## =====================================================================
    ## Step 1: Filter & QC
    ## =====================================================================

    qc <- reactive({
      expr <- dataset$expr
      meta <- dataset$meta
      common <- intersect(colnames(expr), meta$sample)
      meta <- meta[match(common, meta$sample), , drop = FALSE]
      expr <- expr[, common, drop = FALSE]
      expr <- expr[stats::complete.cases(expr), , drop = FALSE]
      validate(need(nrow(meta) >= 15, "Fewer than 15 samples in the current dataset; WGCNA needs more samples to detect modules reliably."))
      list(expr = expr, meta = meta)
    })

    ## Applies the chosen gene-selection strategy (top-N variance, top-median,
    ## variance/MAD percentile, mean floor, custom gene list, or none), then a regex exclude filter.
    gene_selection <- reactive({
      d <- qc()
      expr <- d$expr
      method <- input$gene_filter_method %||% "topvar"
      keep <- switch(method,
        topvar = {
          v <- apply(expr, 1, stats::var)
          n <- min(input$n_genes %||% 2000, nrow(expr))
          order(v, decreasing = TRUE)[seq_len(n)]
        },
        ## Top-N by median expression, per the paper-methodology default (see note below).
        topmedian = {
          m <- apply(expr, 1, stats::median)
          n <- min(input$n_genes_median %||% 5000, nrow(expr))
          order(m, decreasing = TRUE)[seq_len(n)]
        },
        var_pct = {
          v <- apply(expr, 1, stats::var)
          thr <- stats::quantile(v, (input$var_pct %||% 75) / 100, na.rm = TRUE)
          which(v >= thr)
        },
        mad_pct = {
          m <- apply(expr, 1, stats::mad)
          thr <- stats::quantile(m, (input$mad_pct %||% 75) / 100, na.rm = TRUE)
          which(m >= thr)
        },
        mean_cutoff = which(rowMeans(expr) >= (input$mean_cutoff %||% 0)),
        custom_list = {
          genes <- unique(trimws(unlist(strsplit(input$custom_genes %||% "", "[,\n\t ]+"))))
          genes <- genes[nzchar(genes)]
          which(rownames(expr) %in% genes)
        },
        all = seq_len(nrow(expr))
      )

      pattern <- trimws(input$exclude_pattern %||% "")
      if (nzchar(pattern) && length(keep) > 0) {
        matched <- tryCatch(grepl(pattern, rownames(expr)[keep], perl = TRUE), error = function(e) NULL)
        validate(need(!is.null(matched), paste("Invalid regex pattern:", pattern)))
        keep <- keep[!matched]
      }

      min_genes <- if (identical(method, "custom_list")) 20 else 50
      validate(need(length(keep) >= min_genes,
                    sprintf("Fewer than %d genes remain after this filter. Loosen the threshold and try again.", min_genes)))
      list(expr = expr[keep, , drop = FALSE], meta = d$meta, n_kept = length(keep), n_total = nrow(expr))
    })

    ## Sample clustering on the gene-filtered matrix, for outlier detection.
    sample_tree <- reactive({
      gs <- gene_selection()
      stats::hclust(stats::dist(t(gs$expr)), method = "average")
    })

    output$sample_tree_plot <- renderPlot({
      tree <- sample_tree()
      graphics::plot(tree, main = "Sample clustering", sub = "", xlab = "",
                      cex = 0.75, cex.main = 1, col = "#2c6fbb")
      if (isTRUE(input$remove_outliers) && !is.null(input$outlier_height)) {
        graphics::abline(h = input$outlier_height, col = "#c0392b", lty = 2, lwd = 1.5)
      }
    })

    output$outlier_height_ui <- renderUI({
      tree <- sample_tree()
      rng <- range(tree$height)
      sliderInput(ns("outlier_height"), "Cut height (samples below this line are removed)",
                  min = 0, max = ceiling(rng[2]), value = round(rng[2] * 0.9, 1),
                  step = max(round(rng[2] / 100, 2), 0.01))
    })

    ## Cuts the sample tree at the chosen height and keeps only the largest cluster (Zhang & Horvath 2005 outlier rule).
    final_input <- reactive({
      gs <- gene_selection()
      texpr_full <- t(gs$expr)
      removed <- character(0)
      if (isTRUE(input$remove_outliers) && !is.null(input$outlier_height)) {
        tree <- sample_tree()
        clust <- WGCNA::cutreeStatic(tree, cutHeight = input$outlier_height, minSize = 2)
        tab <- table(clust)
        main_label <- as.integer(names(tab)[which.max(tab)])
        keep <- clust == main_label
        removed <- rownames(texpr_full)[!keep]
        texpr_full <- texpr_full[keep, , drop = FALSE]
      }
      validate(need(nrow(texpr_full) >= 15,
                    "Fewer than 15 samples remain after outlier removal; raise the cut height or uncheck \"Remove outlier samples\"."))

      ## WGCNA's own recommended QC gate (Langfelder & Horvath tutorial), run
      ## unconditionally before soft-thresholding/module detection - flags and
      ## drops zero-variance/excess-missing-data genes and samples. Every
      ## gene_filter_method other than variance/MAD-percentile filters
      ## (including "all genes", the default for this project's own dataset)
      ## has no incidental variance floor of its own, so constant genes could
      ## otherwise reach pickSoftThreshold()/blockwiseModules(), which the
      ## WGCNA documentation warns produces undefined correlations.
      removed_genes <- character(0)
      gsg <- WGCNA::goodSamplesGenes(texpr_full, verbose = 0)
      if (!isTRUE(gsg$allOK)) {
        removed_genes <- colnames(texpr_full)[!gsg$goodGenes]
        removed <- union(removed, rownames(texpr_full)[!gsg$goodSamples])
        texpr_full <- texpr_full[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
      }
      validate(need(nrow(texpr_full) >= 15,
                    "Fewer than 15 samples remain after WGCNA's own sample/gene quality check (goodSamplesGenes); raise the cut height, uncheck \"Remove outlier samples\", or loosen the gene filter."))
      validate(need(ncol(texpr_full) >= 20,
                    "Fewer than 20 genes remain after WGCNA's own sample/gene quality check (goodSamplesGenes) removed zero-variance or excess-missing-data genes; loosen the gene filter and try again."))

      meta <- gs$meta[match(rownames(texpr_full), gs$meta$sample), , drop = FALSE]
      list(texpr = texpr_full, meta = meta,
           n_samples_before = ncol(gs$expr), n_samples_after = nrow(texpr_full),
           removed_samples = removed, removed_genes = removed_genes)
    })

    output$gene_filter_summary <- renderUI({
      gs <- gene_selection()
      div(class = "empty-note", icon("check"),
          sprintf("%s of %s genes kept.", format(gs$n_kept, big.mark = ","), format(gs$n_total, big.mark = ",")))
    })

    output$sample_qc_summary <- renderUI({
      fi <- final_input()
      tagList(
        if (length(fi$removed_samples) > 0) {
          div(class = "empty-note", icon("triangle-exclamation"),
              sprintf("%d of %d samples kept. Removed as outliers or by WGCNA's own quality check: %s.",
                      fi$n_samples_after, fi$n_samples_before, paste(fi$removed_samples, collapse = ", ")))
        } else {
          div(class = "empty-note", icon("check"),
              sprintf("%d of %d samples kept.", fi$n_samples_after, fi$n_samples_before))
        },
        if (length(fi$removed_genes %||% character(0)) > 0) {
          gene_list_txt <- paste0(paste(utils::head(fi$removed_genes, 10), collapse = ", "),
                                   if (length(fi$removed_genes) > 10) sprintf(" (+%d more)", length(fi$removed_genes) - 10) else "")
          div(class = "empty-note", icon("triangle-exclamation"),
              sprintf("WGCNA's own quality check (goodSamplesGenes) removed %d gene(s) as zero-variance or excess-missing-data before network construction: %s.",
                      length(fi$removed_genes), gene_list_txt))
        }
      )
    })

    output$step1_ui <- renderUI({
      uploaded <- dataset_is_uploaded()
      tagList(
        box(
          width = 12, title = tagList(icon("dna"), " Gene filtering"), status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Choose which genes go into the network - fewer, more variable genes run faster and add less noise."),
          if (uploaded) div(
            class = "empty-note", style = "margin-bottom: 10px;", icon("book"),
            "This isn't the app's default reference cohort, so the settings below (and in Soft Power / Modules) are pre-filled to match the published methodology's own general-purpose defaults: top 5,000 genes by highest median expression, power β = 7, minimum module size 66, merge cut height 0.3. Adjust any of them before you run."
          ),
          ## Defaults to "all genes" for this project's own dataset (matches Section 2.4.1's no-prefilter choice); uploaded data defaults to the paper's median top-N instead.
          radioButtons(
            ns("gene_filter_method"), "Gene selection method",
            choiceNames = list(
              "Use all genes (recommended: this project's own choice)",
              "Most highly-expressed genes by median (top N) (paper default)",
              "Most variable genes (top N)",
              "Variance percentile cutoff",
              "MAD percentile cutoff",
              "Mean expression cutoff",
              "Custom gene panel (paste your own list)"
            ),
            choiceValues = list("all", "topmedian", "topvar", "var_pct", "mad_pct", "mean_cutoff", "custom_list"),
            selected = if (uploaded) "topmedian" else "all"
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'topmedian'", ns("gene_filter_method")),
            numericInput(ns("n_genes_median"), "Highest-median-expression genes to use", value = 5000, min = 200, max = 20000, step = 100)
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'topvar'", ns("gene_filter_method")),
            numericInput(ns("n_genes"), "Most variable genes to use", value = 2000, min = 200, max = 10000, step = 100)
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'var_pct'", ns("gene_filter_method")),
            sliderInput(ns("var_pct"), "Keep genes at or above this variance percentile", min = 50, max = 99, value = 75, step = 1)
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'mad_pct'", ns("gene_filter_method")),
            sliderInput(ns("mad_pct"), "Keep genes at or above this MAD percentile", min = 50, max = 99, value = 75, step = 1)
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'mean_cutoff'", ns("gene_filter_method")),
            numericInput(ns("mean_cutoff"), "Minimum mean expression", value = 0, step = 0.5)
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'custom_list'", ns("gene_filter_method")),
            textAreaInput(ns("custom_genes"), "Gene symbols (one per line or comma-separated)", rows = 5, placeholder = "TNF\nIL6\nSTAT3\n...")
          ),
          textInput(ns("exclude_pattern"), "Exclude genes matching a pattern (regex, optional)", value = "", placeholder = "e.g. ^HLA for HLA genes"),
          uiOutput(ns("gene_filter_summary"))
        ),
        wgcna_result("users", "Sample quality control",
          desc = "A sample that sits alone on its own branch is usually a technical outlier, not real biology - remove it before building the network.",
          withSpinner(plotOutput(ns("sample_tree_plot"), height = 300), color = "#2c6fbb", type = 6),
          checkboxInput(ns("remove_outliers"), "Remove outlier samples", value = FALSE),
          conditionalPanel(condition = sprintf("input['%s']", ns("remove_outliers")), uiOutput(ns("outlier_height_ui"))),
          uiOutput(ns("sample_qc_summary"))
        )
      )
    })

    ## =====================================================================
    ## Step 2: Soft Power
    ## =====================================================================

    ## reactiveValues + observeEvent (matching Step 3's run_btn/net_store below) rather
    ## than eventReactive: the compute here (pickSoftThreshold + softConnectivity, "the
    ## single slowest step in this app" per global.R) previously ran lazily inside
    ## step2_ui's own renderUI the first time something pulled sft_result() - which
    ## meant clicking "Compute power" invalidated and re-executed that ENTIRE renderUI
    ## (box, radio/slider inputs, the button itself) synchronously around the multi-
    ## minute blocking call, so every click looked like nothing happened (no spinner
    ## reachable until the whole box had already re-rendered) and any real error from
    ## pickSoftThreshold/softConnectivity was swallowed by the `tryCatch(..., error =
    ## function(e) NULL)` guards used everywhere sft_result() was read, indistinguishable
    ## from "button not clicked yet". Running the compute in an observer instead keeps
    ## it out of the render path entirely, and errors are captured into sft_store$error
    ## instead of vanishing.
    sft_store <- reactiveValues(result = NULL, error = NULL)

    observeEvent(input$compute_power_btn, {
      sft_store$error <- NULL
      result <- tryCatch({
        ## Loads the offline script's own soft-power fit instead of re-running
        ## pickSoftThreshold/softConnectivity live over the full ~15.7k-gene default
        ## matrix - mirrors net_store's identical short-circuit for Step 3 below.
        ## Assigns sft_store$result directly (rather than through the `result <-`/
        ## error-handling path below) since return() here exits this whole observer.
        if (isTRUE(dataset$is_bundled_reference)) {
          sft_store$result <- load_precomputed_wgcna_sft()
          return(invisible(NULL))
        }
        texpr <- final_input()$texpr
        powers <- c(1:10, seq(12, 20, 2))
        corName <- wgcna_cor_fnc_name(input$cor_method)
        ## Cached like net_result() below - pickSoftThreshold tests 12 power values, expensive to repeat.
        get_or_compute_wgcna_blocks(
          list(texpr = texpr, powers = powers, network_type = input$network_type,
               cor_method = input$cor_method, r_sq_cutoff = input$r_sq_cutoff, step = "sft"),
          function() {
            sft <- WGCNA::pickSoftThreshold(
              texpr, powerVector = powers, networkType = input$network_type,
              corFnc = wgcna_cor_fnc(input$cor_method), corOptions = list(use = "p"),
              RsquaredCut = input$r_sq_cutoff, verbose = 0
            )
            power <- sft$powerEstimate
            if (is.na(power)) power <- sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)]
            ## Degree-distribution check at the chosen power (the WGCNA paper's own scale-free diagnostic).
            k <- WGCNA::softConnectivity(
              texpr, corFnc = corName, corOptions = "use = 'p'",
              type = input$network_type, power = power, verbose = 0
            )
            list(sft_df = sft$fitIndices, power = power, connectivity = k)
          }
        )
      }, error = function(e) e)
      if (inherits(result, "error")) {
        sft_store$error <- conditionMessage(result)
      } else {
        sft_store$result <- result
      }
    }, ignoreInit = TRUE)

    sft_result <- reactive({
      req(sft_store$result)
      sft_store$result
    })

    effective_power <- reactive({
      if (identical(input$power_mode, "manual")) {
        req(input$manual_power)
        as.numeric(input$manual_power)
      } else {
        sft_result()$power
      }
    })

    output$power_status_ui <- renderUI({
      if (!is.null(sft_store$error)) {
        return(div(class = "empty-note", icon("triangle-exclamation"),
                    paste("Could not compute soft power:", sft_store$error)))
      }
      res <- tryCatch(sft_result(), error = function(e) NULL)
      if (is.null(res)) {
        return(div(class = "empty-note", icon("circle-info"), "Not computed yet. Click \"Compute power\"."))
      }
      eff <- effective_power()
      div(class = "empty-note", icon("check"),
          sprintf("Estimated power: %s. Using power = %s (%s).", res$power, eff,
                  if (identical(input$power_mode, "manual")) "manual override" else "automatic"))
    })

    ## Scale-free fit (R^2) vs power; expression(R^2) for a proper superscript.
    output$sft_r2_plot <- renderPlot({
      res <- tryCatch(sft_result(), error = function(e) NULL)
      req(res)
      df <- res$sft_df
      ggplot(df, aes(x = Power, y = -sign(slope) * SFT.R.sq)) +
        geom_point(color = "#2c6fbb", size = 2) +
        geom_line(color = "#2c6fbb") +
        geom_hline(yintercept = input$r_sq_cutoff, linetype = "dashed", color = "#c0392b") +
        geom_vline(xintercept = effective_power(), linetype = "dotted", color = "#5b6470") +
        labs(x = "Soft-threshold power", y = expression(paste("Scale-free topology fit (", R^2, ")"))) +
        theme_minimal(base_size = 12)
    })

    output$sft_k_plot <- renderPlot({
      res <- tryCatch(sft_result(), error = function(e) NULL)
      req(res)
      df <- res$sft_df
      ggplot(df, aes(x = Power, y = mean.k.)) +
        geom_point(color = "#2c6fbb", size = 2) +
        geom_line(color = "#2c6fbb") +
        geom_vline(xintercept = effective_power(), linetype = "dotted", color = "#5b6470") +
        labs(x = "Soft-threshold power", y = "Mean connectivity") +
        theme_minimal(base_size = 12)
    })

    ## Truncated-exponential fit R^2 vs power (pickSoftThreshold's alternative goodness-of-fit).
    output$sft_trunc_plot <- renderPlot({
      res <- tryCatch(sft_result(), error = function(e) NULL)
      req(res)
      df <- res$sft_df
      ggplot(df, aes(x = Power, y = truncated.R.sq)) +
        geom_point(color = "#2c6fbb", size = 2) +
        geom_line(color = "#2c6fbb") +
        geom_vline(xintercept = effective_power(), linetype = "dotted", color = "#5b6470") +
        labs(x = "Soft-threshold power", y = expression(paste("Exponential fit (", R^2, ")"))) +
        theme_minimal(base_size = 12)
    })

    ## Reimplements WGCNA::scaleFreePlot()'s fit (same binning/models) manually
    ## since scaleFreePlot always forces its own title text with no way to suppress it.
    scale_free_check_stats <- reactive({
      res <- tryCatch(sft_result(), error = function(e) NULL)
      req(res)
      validate(need(!is.null(res$connectivity),
                    "Not available for the precomputed default dataset (this diagnostic needs a live softConnectivity() run)."))
      k <- res$connectivity
      k <- k[is.finite(k)]
      validate(need(length(k) >= 20, "Too few genes have a valid connectivity value to check the scale-free fit."))

      n_breaks <- 10
      discretized_k <- cut(k, n_breaks)
      dk <- tapply(k, discretized_k, mean)
      p_dk <- as.vector(tapply(k, discretized_k, length) / length(k))
      breaks1 <- seq(from = min(k), to = max(k), length = n_breaks + 1)
      hist1 <- suppressWarnings(graphics::hist(k, breaks = breaks1, equidist = FALSE, plot = FALSE, right = TRUE))
      dk2 <- hist1$mids
      dk <- ifelse(is.na(dk), dk2, dk)
      dk <- ifelse(dk == 0, dk2, dk)
      p_dk <- ifelse(is.na(p_dk), 0, p_dk)
      log_dk <- as.vector(log10(dk))
      log_pdk <- as.numeric(log10(p_dk + 1e-09))

      lm_linear <- stats::lm(log_pdk ~ log_dk)
      lm_exp <- tryCatch(stats::lm(log_pdk ~ log_dk + I(10^log_dk)), error = function(e) NULL)

      list(
        log_dk = log_dk, log_pdk = log_pdk,
        fit_linear = as.numeric(stats::predict(lm_linear)),
        fit_exp = if (!is.null(lm_exp)) as.numeric(stats::predict(lm_exp)) else NULL,
        r2_linear = round(summary(lm_linear)$adj.r.squared, 2),
        slope = round(lm_linear$coefficients[[2]], 2),
        r2_exp = if (!is.null(lm_exp)) round(summary(lm_exp)$adj.r.squared, 2) else NA_real_
      )
    })

    draw_scale_free_check <- function() {
      st <- scale_free_check_stats()
      graphics::plot(st$log_dk, st$log_pdk, xlab = "log10(k)", ylab = "log10(p(k))", main = "", pch = 1, col = "black")
      graphics::lines(st$log_dk, st$fit_linear, col = "black")
      if (!is.null(st$fit_exp)) graphics::lines(st$log_dk, st$fit_exp, col = "red")
    }

    ## Catches draw-time errors and turns them into a validate() message instead of a raw error.
    output$scale_free_check_plot <- renderPlot({
      ok <- tryCatch({ draw_scale_free_check(); TRUE },
                      error = function(e) { message("WGCNA scale_free_check_plot failed: ", conditionMessage(e)); FALSE })
      validate(need(ok, "Could not draw the scale-free check plot for the current power."))
    })

    output$scale_free_check_annotation <- renderUI({
      st <- scale_free_check_stats()
      p(class = "submodule-desc", HTML(sprintf(
        "Fit quality: R<sup>2</sup> = %.2f (slope %.2f) for the black line &bull; R<sup>2</sup> = %s for the red (exponential) line.",
        st$r2_linear, st$slope, if (is.na(st$r2_exp)) "n/a" else sprintf("%.2f", st$r2_exp)
      )))
    })

    output$download_scale_free_png <- downloadHandler(
      filename = function() "wgcna_scale_free_check.png",
      content = function(file) {
        grDevices::png(file, width = 1400, height = 1000, res = 150)
        draw_scale_free_check()
        grDevices::dev.off()
      }
    )

    ## Controls only - deliberately NOT dependent on sft_result()/sft_store, so this box
    ## (and the compute_power_btn inside it) doesn't get torn down and rebuilt - resetting
    ## every input back to its hardcoded default - each time a computation finishes. See
    ## the sft_store observer above for why the compute itself was pulled out of the render path.
    output$step2_controls_ui <- renderUI({
      uploaded <- dataset_is_uploaded()
      box(
        width = 12, title = tagList(icon("wave-square"), " Network settings"), status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "WGCNA raises gene-gene correlations to a power so the network is approximately scale-free. Pick the lowest power that reaches a good fit below."),
        fluidRow(
          column(
            6,
            ## Defaults to signed for this project's own data (Section 2.4.7:
            ## avoids merging anti-correlated genes); unsigned (the WGCNA
            ## package default) for uploaded data, since networkType is unstated in the paper.
            radioButtons(ns("network_type"), "Network type",
                         choiceNames = list(
                           if (uploaded) "Signed" else "Signed (recommended)",
                           if (uploaded) "Unsigned (recommended: paper default - unstated in the methodology, so left at the WGCNA package default)" else "Unsigned",
                           "Signed hybrid"
                         ),
                         choiceValues = list("signed", "unsigned", "signed hybrid"),
                         selected = if (uploaded) "unsigned" else "signed"),
            radioButtons(ns("cor_method"), "Correlation method",
                         choiceNames = list("Pearson (standard)", "Bicor (robust to outliers)"),
                         choiceValues = list("pearson", "bicor"), selected = "pearson")
          ),
          column(
            6,
            ## R^2 >= 0.9 cutoff, per this project's own joint fit/connectivity criterion (Section 2.4.8).
            sliderInput(ns("r_sq_cutoff"), HTML("Target scale-free fit (R<sup>2</sup>)"), min = 0.7, max = 0.95, value = 0.9, step = 0.01),
            ## Defaults to manual power 12 for this project's data (Section 2.4.5, rejecting the auto-estimate's runaway connectivity); manual 7 (the paper's value) for uploaded data.
            radioButtons(ns("power_mode"), "Power",
                         choiceNames = list(
                           "Automatic",
                           if (uploaded) "Manual (recommended: 7, paper default)" else "Manual (recommended: 12, this project's own choice)"
                         ),
                         choiceValues = list("auto", "manual"), selected = "manual"),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'manual'", ns("power_mode")),
              numericInput(ns("manual_power"), "Power value", value = if (uploaded) 7 else 12, min = 1, max = 30, step = 1)
            )
          )
        ),
        actionButton(ns("compute_power_btn"), "Compute power", icon = icon("play"), class = "btn-primary btn-sm"),
        if (!isTRUE(dataset$is_bundled_reference)) p(class = "submodule-desc", style = "margin-top:6px;",
          "Can take several minutes depending on gene count and dataset size - narrowing Step 1's gene filter (e.g. the top 2,000-4,000 most variable genes instead of all genes) runs much faster."),
        div(style = "margin-top:8px;", uiOutput(ns("power_status_ui")))
      )
    })

    ## Gate + results only - re-renders when sft_result()/sft_store changes, but the
    ## controls box above (with the button itself) is a separate output and stays put.
    output$step2_ui <- renderUI({
      ## Gate: no result boxes until "Compute power" has been clicked.
      if (!is.null(sft_store$error) || is.null(tryCatch(sft_result(), error = function(e) NULL))) {
        div(class = "empty-note", icon("circle-info"),
            "Not computed yet. Set your options above, then click \"Compute power\" to see the diagnostics below.")
      } else tagList(
          wgcna_result("chart-line", "Power diagnostics",
            wgcna_result_row(
              column(4, div(class = "wgcna-result-subtitle", "Scale-free fit"),
                     withSpinner(plotOutput(ns("sft_r2_plot"), height = 280), color = "#2c6fbb", type = 6)),
              column(4, div(class = "wgcna-result-subtitle", "Mean connectivity"),
                     withSpinner(plotOutput(ns("sft_k_plot"), height = 280), color = "#2c6fbb", type = 6)),
              column(4, div(class = "wgcna-result-subtitle", "Exponential fit"),
                     withSpinner(plotOutput(ns("sft_trunc_plot"), height = 280), color = "#2c6fbb", type = 6))
            )
          ),
          wgcna_result("magnifying-glass-chart", "Scale-free topology check (at the power in use)",
            desc = "The actual degree-distribution diagnostic from the WGCNA paper: log10(p(k)) vs log10(k) should be roughly linear.",
            withSpinner(plotOutput(ns("scale_free_check_plot"), height = 420), color = "#2c6fbb", type = 6),
            uiOutput(ns("scale_free_check_annotation")),
            div(class = "table-toolbar", downloadButton(ns("download_scale_free_png"), "Download figure (PNG)", class = "btn-sm"))
          )
        )
    })

    ## =====================================================================
    ## Step 3: Modules
    ## =====================================================================

    ## "Run" loads the precomputed result for the default dataset, or computes live for uploaded/preprocessed data.
    net_store <- reactiveValues(result = NULL, source = NULL, error = NULL)

    ## Runs blockwiseModules() (or loads the precomputed result for the example dataset); errors are caught and stored since this isn't a render context.
    observeEvent(input$run_btn, {
      net_store$error <- NULL
      result <- tryCatch({
      if (isTRUE(dataset$is_bundled_reference)) {
        net_store$result <- load_precomputed_wgcna_result()
        net_store$source <- "loaded"
        return(invisible(NULL))
      }

      fi <- final_input()
      texpr <- fi$texpr
      power <- tryCatch(effective_power(), error = function(e) NA)
      validate(need(!is.null(power) && length(power) == 1 && is.finite(power),
                    "Compute the network power in Step 2 first, or switch to a manual power there."))

      ## Cached by digest of matrix + settings (get_or_compute_wgcna_blocks in global.R).
      ## reassignThreshold: 0 for this project's data (Table 2.x, membership by topological overlap alone); package default 1e-06 for uploaded data, since the paper doesn't mention it.
      reassign_threshold <- if (dataset_is_uploaded()) 1e-06 else 0
      net <- get_or_compute_wgcna_blocks(
        list(texpr = texpr, power = power, network_type = input$network_type,
             tom_type = input$tom_type, cor_method = input$cor_method,
             deep_split = input$deep_split, min_module_size = input$min_module_size,
             merge_cut_height = input$merge_cut_height,
             pam_respects_dendro = isTRUE(input$pam_respects_dendro),
             reassign_threshold = reassign_threshold),
        function() WGCNA::blockwiseModules(
          texpr,
          power = power, networkType = input$network_type, TOMType = input$tom_type,
          corType = input$cor_method,
          deepSplit = input$deep_split, minModuleSize = input$min_module_size,
          reassignThreshold = reassign_threshold, mergeCutHeight = input$merge_cut_height,
          pamRespectsDendro = isTRUE(input$pam_respects_dendro),
          ## Fixed seed (matches 06_WGCNA.R's CFG$seed = 1234) - PAM-stage reassignment is stochastic, and module colors are size-rank based, so an unpinned seed could relabel modules between runs.
          randomSeed = 1234,
          numericLabels = FALSE, maxBlockSize = ncol(texpr) + 1, verbose = 0
        )
      )
      module_colors <- net$colors
      MEs <- net$MEs

      gene_module <- data.frame(gene = colnames(texpr), module = module_colors, stringsAsFactors = FALSE)
      gene_module <- gene_module[order(gene_module$module), ]

      net_store$result <- list(
        texpr = texpr, meta = fi$meta, power = power,
        network_type = input$network_type, tom_type = input$tom_type, cor_method = input$cor_method,
        module_colors = module_colors, MEs = MEs,
        dendro = net$dendrograms[[1]], block_colors = module_colors[net$blockGenes[[1]]],
        gene_module = gene_module,
        module_sizes = as.data.frame(table(module = module_colors)),
        hub_table_precomputed = NULL,
        n_genes = ncol(texpr), n_samples = nrow(texpr),
        n_modules = length(unique(module_colors[module_colors != "grey"]))
      )
      net_store$source <- "computed"
      }, error = function(e) e)
      if (inherits(result, "error")) net_store$error <- conditionMessage(result)
    })

    net_result <- reactive({
      req(net_store$result)
      net_store$result
    })

    ## Publishes module results into shared results$wgcna, preserving prior Step 4 fields (accumulate pattern), for the Candidate Gene tab's module-genes picker.
    observeEvent(net_result(), {
      res <- net_result()
      base <- results$wgcna %||% list()
      base$n_genes <- res$n_genes
      base$n_samples <- res$n_samples
      base$soft_power <- res$power
      base$network_type <- res$network_type
      base$n_modules <- res$n_modules
      base$module_genes <- split(res$gene_module$gene, res$gene_module$module)
      results$wgcna <- base
    })

    output$module_run_status_ui <- renderUI({
      if (!is.null(net_store$error)) {
        return(div(class = "empty-note", icon("triangle-exclamation"), net_store$error))
      }
      res <- tryCatch(net_result(), error = function(e) NULL)
      if (is.null(res)) {
        return(div(class = "empty-note", icon("circle-info"), "Not run yet. Click \"Run\" above."))
      }
      div(class = "empty-note", icon("check"),
          sprintf("%s genes, %s samples, %s modules detected (grey = unassigned), power %s.",
                  res$n_genes, res$n_samples, res$n_modules, res$power))
    })

    draw_dendro <- function() {
      net <- tryCatch(net_result(), error = function(e) NULL)
      req(net)
      WGCNA::plotDendroAndColors(
        net$dendro, net$block_colors, "Module colors",
        dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05
      )
    }

    output$dendro_plot <- renderPlot({ draw_dendro() })

    output$download_dendro_png <- downloadHandler(
      filename = function() "wgcna_dendrogram.png",
      content = function(file) {
        grDevices::png(file, width = 1400, height = 900, res = 150)
        draw_dendro()
        grDevices::dev.off()
      }
    )

    ## Gene-network heatmap (TOM, per Langfelder & Horvath tutorial's TOMplot()):
    ## computes the full TOM once, then subsamples genes for a renderable plot.
    ## Gated behind its own button since TOMsimilarityFromExpr is as expensive as blockwiseModules.
    tom_result <- eventReactive(input$compute_tom_btn, {
      net <- net_result()
      full_tom <- get_or_compute_wgcna_blocks(
        list(texpr = net$texpr, power = net$power, network_type = net$network_type,
             tom_type = net$tom_type, cor_method = net$cor_method, step = "tom_full"),
        function() WGCNA::TOMsimilarityFromExpr(
          net$texpr, power = net$power, networkType = net$network_type,
          TOMType = net$tom_type, corType = net$cor_method, verbose = 0
        )
      )
      n_genes <- ncol(net$texpr)
      n_select <- min(input$tom_n_select %||% 400, n_genes)
      set.seed(1234) ## fixed so the subsample doesn't reshuffle on every re-render
      select <- sample(seq_len(n_genes), size = n_select)
      diss_tom <- 1 - full_tom
      select_tom <- diss_tom[select, select]
      select_tree <- stats::hclust(stats::as.dist(select_tom), method = "average")
      select_colors <- net$module_colors[select]
      ## ^7 and NA diagonal match TOMplot()'s recipe: raises contrast, blanks the self-similar diagonal.
      plot_diss <- select_tom^7
      diag(plot_diss) <- NA
      list(plot_diss = plot_diss, select_tree = select_tree, select_colors = select_colors,
           n_select = n_select, n_total = n_genes)
    })

    draw_tom_plot <- function() {
      tr <- tom_result()
      WGCNA::TOMplot(tr$plot_diss, tr$select_tree, tr$select_colors, main = "")
    }

    output$tom_plot <- renderPlot({
      ok <- tryCatch({ draw_tom_plot(); TRUE },
                      error = function(e) { message("WGCNA tom_plot failed: ", conditionMessage(e)); FALSE })
      validate(need(ok, "Could not draw the network heatmap for the current settings."))
    })

    output$tom_status_ui <- renderUI({
      res <- tryCatch(tom_result(), error = function(e) NULL)
      if (is.null(res)) {
        return(div(class = "empty-note", icon("circle-info"),
                    "Not computed yet. Click \"Compute network heatmap\" - this recomputes the full topological overlap matrix and can take as long as Step 3's own \"Run\"."))
      }
      div(class = "empty-note", icon("check"),
          sprintf("Showing a random subsample of %d genes (of %d in the network).", res$n_select, res$n_total))
    })

    output$download_tom_png <- downloadHandler(
      filename = function() "wgcna_network_heatmap.png",
      content = function(file) {
        grDevices::png(file, width = 1400, height = 1400, res = 150)
        draw_tom_plot()
        grDevices::dev.off()
      }
    )

    ## Module size, colored by the module's own WGCNA color - module names
    ## in this app ARE valid R color names, so fill can map directly to them.
    output$module_size_plot <- renderPlot({
      net <- tryCatch(net_result(), error = function(e) NULL)
      req(net)
      df <- net$module_sizes
      ggplot(df, aes(x = stats::reorder(module, Freq), y = Freq, fill = module)) +
        geom_col() + coord_flip() +
        scale_fill_manual(values = stats::setNames(as.character(df$module), df$module), guide = "none") +
        labs(x = NULL, y = "Genes") +
        theme_minimal(base_size = 12)
    })

    ## Eigengene network: how correlated modules are with each other -
    ## highly correlated modules describe overlapping biology.
    output$eigengene_network_plot <- renderPlot({
      net <- tryCatch(net_result(), error = function(e) NULL)
      req(net)
      me_cor <- wgcna_cor_fnc(net$cor_method)(as.matrix(net$MEs), use = "p")
      rownames(me_cor) <- colnames(me_cor) <- sub("^ME", "", colnames(net$MEs))
      df <- as.data.frame(as.table(me_cor)); colnames(df) <- c("module1", "module2", "cor")
      ## Cell/axis text scales down as module count grows, to avoid label overlap.
      n_mod <- length(unique(df$module1))
      cell_text_size <- max(1.6, min(2.7, 30 / n_mod))
      axis_text_size <- max(6, min(9, 90 / n_mod))
      ggplot(df, aes(module1, module2, fill = cor)) +
        geom_tile(color = "white") +
        geom_text(aes(label = sprintf("%.2f", cor)), size = cell_text_size) +
        scale_fill_gradient2(low = "#2c6fbb", mid = "white", high = "#c0392b", midpoint = 0, limits = c(-1, 1)) +
        labs(x = NULL, y = NULL, fill = "r") +
        coord_fixed() +
        theme_minimal(base_size = 11) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
              axis.text.y = element_text(size = axis_text_size))
    })

    output$module_table <- DT::renderDataTable({
      net <- tryCatch(net_result(), error = function(e) NULL)
      req(net)
      DT::datatable(net$gene_module, rownames = FALSE, filter = "top",
                     options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_modules <- downloadHandler(
      filename = function() "wgcna_gene_modules.csv",
      content = function(file) write.csv(net_result()$gene_module, file, row.names = FALSE)
    )

    ## Exports the full net_result() list (colors, MEs, dendrogram, settings) as .rds, for downstream use in R.
    output$download_wgcna_rds <- downloadHandler(
      filename = function() "wgcna_module_detection_result.rds",
      content = function(file) saveRDS(net_result(), file)
    )

    output$step3_ui <- renderUI({
      ## minModuleSize/mergeCutHeight default to the paper's values (66/0.3) for
      ## uploaded data vs this project's own choices (30/0.25); pamRespectsDendro
      ## defaults to the package default (TRUE) for uploaded data since the paper doesn't state it.
      uploaded <- dataset_is_uploaded()
      tagList(
        box(
          width = 12, title = tagList(icon("diagram-project"), " Module detection settings"), status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Genes are clustered by co-expression and cut into modules; similar modules are then merged into one."),
          fluidRow(
            column(
              6,
              numericInput(ns("min_module_size"), "Minimum module size", value = if (uploaded) 66 else 30, min = 5, max = 1000, step = 5),
              sliderInput(ns("deep_split"), "Deep split (0 = coarse, 4 = fine)", min = 0, max = 4, value = 2, step = 1)
            ),
            column(
              6,
              sliderInput(ns("merge_cut_height"), "Merge cut height", min = 0, max = 1, value = if (uploaded) 0.3 else 0.25, step = 0.01),
              selectInput(ns("tom_type"), "TOM type", choices = c("Signed" = "signed", "Unsigned" = "unsigned", "Signed Nowick" = "signed Nowick"),
                          selected = "signed", selectize = FALSE),
              checkboxInput(ns("pam_respects_dendro"), "PAM stage respects dendrogram", value = uploaded)
            )
          ),
          ## Instant for the default dataset (loads precomputed result); live compute for uploaded data.
          actionButton(ns("run_btn"), "Run", icon = icon("play"), class = "btn-primary btn-sm"),
          if (!isTRUE(dataset$is_bundled_reference)) p(class = "submodule-desc", style = "margin-top:6px;",
            "Can take several minutes to tens of minutes depending on gene count and dataset size - narrowing Step 1's gene filter (e.g. the top 2,000-4,000 most variable genes instead of all genes) runs much faster for quick exploration."),
          div(style = "margin-top:8px;", uiOutput(ns("module_run_status_ui")))
        ),
        ## Gate: no result boxes until "Run" has been clicked.
        if (is.null(tryCatch(net_result(), error = function(e) NULL))) {
          div(class = "empty-note", icon("circle-info"),
              "Not run yet. Set your module-detection options above, then click \"Run\" to see the results below.")
        } else tagList(
          wgcna_result("diagram-project", "Gene dendrogram & modules",
            withSpinner(plotOutput(ns("dendro_plot"), height = 400), color = "#2c6fbb", type = 6),
            div(class = "table-toolbar", downloadButton(ns("download_dendro_png"), "Download dendrogram (PNG)", class = "btn-sm"))
          ),
          ## Uploaded dataset only; opt-in button since recomputing the full TOM is expensive.
          if (uploaded) wgcna_result("table-cells", "Gene network heatmap (TOM)",
            desc = "Topological overlap between genes, on a random subsample (full topological overlap at thousands of genes can't be rendered directly). Darker = higher shared connectivity; the color strip along each edge marks module membership.",
            fluidRow(
              column(4, numericInput(ns("tom_n_select"), "Genes to sample for this plot", value = 400, min = 50, max = 1000, step = 50)),
              column(8, div(style = "margin-top:24px;",
                            actionButton(ns("compute_tom_btn"), "Compute network heatmap", icon = icon("play"), class = "btn-primary btn-sm")))
            ),
            div(style = "margin-bottom:8px;", uiOutput(ns("tom_status_ui"))),
            withSpinner(plotOutput(ns("tom_plot"), height = 420), color = "#2c6fbb", type = 6),
            div(class = "table-toolbar", downloadButton(ns("download_tom_png"), "Download figure (PNG)", class = "btn-sm"))
          ),
          wgcna_result_row(
            column(6, div(class = "wgcna-result-subtitle", "Module sizes"),
                   withSpinner(plotOutput(ns("module_size_plot"), height = 280), color = "#2c6fbb", type = 6)),
            column(6, div(class = "wgcna-result-subtitle", "Eigengene network"),
                   withSpinner(plotOutput(ns("eigengene_network_plot"), height = 420), color = "#2c6fbb", type = 6))
          ),
          wgcna_result("table-list", "Gene-module assignment",
            div(class = "table-toolbar",
                downloadButton(ns("download_modules"), "Download CSV", class = "btn-sm"),
                downloadButton(ns("download_wgcna_rds"), "Full result (RDS)", class = "btn-sm")),
            DT::dataTableOutput(ns("module_table"))
          )
        )
      )
    })

    ## =====================================================================
    ## Step 4: Module-Trait
    ## =====================================================================

    eligible_traits <- reactive({
      meta <- qc()$meta
      cols <- setdiff(colnames(meta), "sample")
      Filter(function(cl) {
        v <- meta[[cl]]
        !all(is.na(v)) && length(unique(stats::na.omit(v))) >= 2
      }, cols)
    })

    output$trait_picker_ui <- renderUI({
      cols <- eligible_traits()
      validate(need(length(cols) > 0, "No usable trait columns were found in the current dataset's metadata."))
      default_sel <- intersect(c("group", "sex"), cols)
      if (length(default_sel) == 0) default_sel <- cols[1]
      checkboxGroupInput(ns("trait_cols"), NULL, choices = cols, selected = default_sel, inline = TRUE)
    })

    ## Per-categorical-trait "which levels?" checkboxes, so e.g. a 3-level group column can be narrowed to just RA vs HC.
    output$trait_level_filters_ui <- renderUI({
      req(input$trait_cols)
      meta <- qc()$meta
      cat_traits <- Filter(function(cl) !is.numeric(meta[[cl]]), input$trait_cols)
      if (length(cat_traits) == 0) return(NULL)
      tagList(
        p(class = "submodule-desc", "Each level you keep below becomes its own column in the heatmap (e.g. group -> \"RA\" and \"HC\" columns) - untick a level (e.g. \"other\") to drop it from the comparison."),
        lapply(cat_traits, function(cl) {
          lv <- sort(unique(stats::na.omit(as.character(meta[[cl]]))))
          checkboxGroupInput(ns(wgcna_trait_widget_id(cl, eligible_traits())), sprintf("\"%s\" levels to include", cl),
                              choices = lv, selected = lv, inline = TRUE)
        })
      )
    })

    ## Builds the trait matrix: each selected level of a categorical trait becomes its own binary column
    ## (samples outside the chosen levels are NA, not 0); numeric traits pass through as-is.
    module_trait <- eventReactive(input$run_traits_btn, {
      net <- net_result()
      validate(need(length(input$trait_cols) > 0, "Select at least one trait above."))
      traits <- data.frame(row.names = rownames(net$texpr))
      used_names <- character(0)
      uniquify <- function(nm, cl) if (nm %in% used_names) paste0(nm, " (", cl, ")") else nm

      for (cl in input$trait_cols) {
        v <- net$meta[[cl]]
        if (is.numeric(v)) {
          nm <- uniquify(cl, cl)
          traits[[nm]] <- as.numeric(v)
          used_names <- c(used_names, nm)
        } else {
          v <- as.character(v)
          lv <- wgcna_trait_levels(input, net$meta, cl, eligible_traits())
          if (is.null(lv) || length(lv) == 0) lv <- sort(unique(stats::na.omit(v)))
          validate(need(length(lv) >= 2, sprintf("Select at least two levels for \"%s\" above.", cl)))
          in_set <- v %in% lv
          for (level in lv) {
            nm <- uniquify(level, cl)
            ind <- rep(NA_real_, length(v))
            ind[in_set] <- as.numeric(v[in_set] == level)
            traits[[nm]] <- ind
            used_names <- c(used_names, nm)
          }
        }
      }

      ## Drop grey (WGCNA's unassigned-gene bin) - not a real co-expression module.
      MEs_use <- net$MEs[, colnames(net$MEs) != "MEgrey", drop = FALSE]
      ## Uploaded data: cluster-order eigengenes (WGCNA::orderMEs) so correlated modules sit together.
      if (dataset_is_uploaded()) MEs_use <- WGCNA::orderMEs(MEs_use)
      cor_mat <- wgcna_cor_fnc(net$cor_method)(MEs_use, traits, use = "p")
      ## Per-column sample size, since a filtered categorical trait excludes some samples via NA.
      n_per_trait <- vapply(traits, function(x) sum(!is.na(x)), numeric(1))
      n_mat <- matrix(rep(n_per_trait, each = nrow(cor_mat)), nrow = nrow(cor_mat), dimnames = dimnames(cor_mat))
      p_mat <- WGCNA::corPvalueStudent(cor_mat, n_mat)
      rownames(cor_mat) <- sub("^ME", "", rownames(cor_mat))
      rownames(p_mat) <- sub("^ME", "", rownames(p_mat))
      list(cor = cor_mat, p = p_mat)
    }, ignoreInit = TRUE)

    wgcna_traits_has_run <- reactiveVal(FALSE)
    observeEvent(input$run_traits_btn, wgcna_traits_has_run(TRUE), ignoreInit = TRUE)

    observeEvent(module_trait(), {
      mt <- module_trait()
      ## Disease-module cutoff: |r| >= 0.5 AND p < 1e-8 (Section 2.4.9), not just p < 0.05.
      hit <- (abs(mt$cor) >= 0.5) & (mt$p < 1e-8)
      sig <- rownames(mt$p)[apply(hit, 1, any, na.rm = TRUE)]
      base <- results$wgcna %||% list()
      base$significant_trait_modules <- sig
      base$traits_tested <- colnames(mt$p)
      results$wgcna <- base
    })

    mt_plot_obj <- reactive({
      mt <- module_trait()
      cor_long <- as.data.frame(as.table(mt$cor)); colnames(cor_long) <- c("module", "trait", "cor")
      p_long <- as.data.frame(as.table(mt$p)); colnames(p_long) <- c("module", "trait", "p")
      df <- merge(cor_long, p_long, by = c("module", "trait"))
      ## Uploaded data: restore eigengene-clustered row order (merge() alphabetises it) so the heatmap reads top-to-bottom in clustered order.
      uploaded <- dataset_is_uploaded()
      if (uploaded) df$module <- factor(df$module, levels = rev(rownames(mt$cor)))
      ## Uploaded data: two-line "r / (P)" cell label matching the published figure's style.
      df$cell_label <- if (uploaded) {
        sprintf("%.2f\n(%s)", df$cor, formatC(df$p, format = "e", digits = 0))
      } else {
        sprintf("%.2f", df$cor)
      }
      ggplot(df, aes(x = trait, y = module, fill = cor)) +
        geom_tile(color = "white") +
        geom_text(aes(label = cell_label), size = if (uploaded) 2.6 else 3.2, lineheight = 0.85) +
        scale_fill_gradient2(low = "#2c6fbb", mid = "white", high = "#c0392b", midpoint = 0, limits = c(-1, 1)) +
        labs(x = NULL, y = "Module", fill = "r") +
        theme_minimal(base_size = 12)
    })

    output$mt_plot <- renderPlot({ mt_plot_obj() })

    output$download_mt_png <- downloadHandler(
      filename = function() "wgcna_module_trait_heatmap.png",
      content = function(file) ggsave(file, plot = mt_plot_obj(), width = 8, height = 6, dpi = 300, bg = "white")
    )

    output$download_mt_csv <- downloadHandler(
      filename = function() "wgcna_module_trait_correlation.csv",
      content = function(file) {
        mt <- module_trait()
        cor_long <- as.data.frame(as.table(mt$cor)); colnames(cor_long) <- c("module", "trait", "cor")
        p_long <- as.data.frame(as.table(mt$p)); colnames(p_long) <- c("module", "trait", "p")
        write.csv(merge(cor_long, p_long, by = c("module", "trait")), file, row.names = FALSE)
      }
    )

    output$step4_ui <- renderUI({
      tagList(
        box(
          width = 12, title = tagList(icon("table-cells"), " Traits"), status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Any metadata column with at least two values can be used as a trait, not just group and sex."),
          uiOutput(ns("trait_picker_ui")),
          uiOutput(ns("trait_level_filters_ui")),
          actionButton(ns("run_traits_btn"), "Compute module-trait correlations", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        if (!wgcna_traits_has_run()) {
          div(class = "empty-note", icon("circle-info"),
              "Not run yet. Pick traits above, then click \"Compute module-trait correlations\".")
        } else {
          wgcna_result("table-cells", "Module-trait correlation",
            desc = "Strongly red or blue cells mark modules whose overall expression tracks that trait.",
            withSpinner(plotOutput(ns("mt_plot"), height = 360), color = "#2c6fbb", type = 6),
            div(class = "table-toolbar",
                downloadButton(ns("download_mt_png"), "Download heatmap (PNG)", class = "btn-sm"),
                downloadButton(ns("download_mt_csv"), "Download table (CSV)", class = "btn-sm"))
          )
        }
      )
    })

    ## =====================================================================
    ## Step 5: Hub Genes
    ## =====================================================================

    ## Intramodular connectivity (kWithin), computed once per module-detection run and reused across module/trait picks.
    intramod_conn <- reactive({
      net <- net_result()
      ik <- WGCNA::intramodularConnectivity.fromExpr(
        net$texpr, colors = net$module_colors,
        corFnc = wgcna_cor_fnc_name(net$cor_method), corOptions = "use = 'p'",
        networkType = net$network_type, power = net$power
      )
      rownames(ik) <- colnames(net$texpr)
      ik
    })

    output$hub_picker_ui <- renderUI({
      net <- net_result()
      mods <- sort(setdiff(unique(net$module_colors), "grey"))
      validate(need(length(mods) > 0, "No non-grey modules were detected. Try a smaller minimum module size in Step 3."))
      traits <- eligible_traits()
      validate(need(length(traits) > 0, "No usable trait columns were found in the current dataset's metadata."))
      fluidRow(
        column(6, selectInput(ns("hub_module"), "Module", choices = mods, selectize = FALSE)),
        column(6, selectInput(ns("hub_trait"), "Trait", choices = traits, selectize = FALSE))
      )
    })

    ## Same "which levels?" idea as Step 4, for the single trait picked here; hidden for numeric traits.
    output$hub_trait_levels_ui <- renderUI({
      req(input$hub_trait)
      meta <- net_result()$meta
      if (is.numeric(meta[[input$hub_trait]])) return(NULL)
      lv <- sort(unique(stats::na.omit(as.character(meta[[input$hub_trait]]))))
      checkboxGroupInput(ns("hub_trait_levels"), sprintf("\"%s\" levels to include", input$hub_trait),
                          choices = lv, selected = lv, inline = TRUE)
    })

    hub_data <- eventReactive(input$run_hubs_btn, {
      net <- net_result()
      req(input$hub_module, input$hub_trait)
      me_col <- paste0("ME", input$hub_module)
      validate(need(me_col %in% colnames(net$MEs), "Selected module not found - re-run Step 3."))

      ## Use the precomputed hub table only when module+trait matches exactly what it was built for
      ## (disease module, trait = "group", no level filter); otherwise fall through to live computation.
      pre <- net$hub_table_precomputed
      lv <- if (!is.numeric(net$meta[[input$hub_trait]])) input$hub_trait_levels else NULL
      use_precomputed <- !is.null(pre) && identical(input$hub_trait, "group") &&
        input$hub_module %in% unique(pre$module) &&
        (is.null(lv) || setequal(lv, sort(unique(stats::na.omit(as.character(net$meta$group))))))

      if (use_precomputed) {
        sub <- pre[pre$module == input$hub_module, , drop = FALSE]
        ## is_hub kept as-is from the pipeline (kME/GS_RA here are rounded for display, not full precision).
        df <- data.frame(gene = sub$gene, kME = sub$kME, GS = sub$GS_RA, kWithin = sub$connectivity,
                          is_hub = sub$is_hub, stringsAsFactors = FALSE)
      } else {
        corFnc <- wgcna_cor_fnc(net$cor_method)
        kme <- as.numeric(corFnc(net$texpr, net$MEs[[me_col]], use = "p"))
        names(kme) <- colnames(net$texpr)
        trait_vec <- wgcna_encode_trait(net$meta, input$hub_trait, levels_keep = lv)
        validate(need(length(unique(stats::na.omit(trait_vec))) >= 2,
                      "Select at least two levels for this trait above."))
        gs <- as.numeric(corFnc(net$texpr, trait_vec, use = "p"))
        names(gs) <- colnames(net$texpr)
        ik <- intramod_conn()
        module_genes <- net$gene_module$gene[net$gene_module$module == input$hub_module]
        df <- data.frame(
          gene = module_genes, kME = kme[module_genes], GS = gs[module_genes],
          kWithin = ik[module_genes, "kWithin"], is_hub = NA, stringsAsFactors = FALSE
        )
      }
      df$abs_kME <- abs(df$kME); df$abs_GS <- abs(df$GS)
      df
    }, ignoreInit = TRUE)

    wgcna_hubs_has_run <- reactiveVal(FALSE)
    observeEvent(input$run_hubs_btn, wgcna_hubs_has_run(TRUE), ignoreInit = TRUE)

    ## Filters hub_data() by the kME/GS thresholds; explicit !is.na(...) guards
    ## avoid NA rows surviving subsetting and breaking the network plot downstream.
    hub_filtered <- reactive({
      df <- hub_data()
      ## At the documented default threshold, use the pipeline's own is_hub flag (exact match) rather than re-deriving from rounded display columns.
      at_default <- isTRUE(all.equal(input$kme_thr, 0.8)) && isTRUE(all.equal(input$gs_thr, 0.2))
      if (at_default && !is.null(df$is_hub) && !any(is.na(df$is_hub))) {
        keep <- df$is_hub
      } else {
        keep <- !is.na(df$abs_kME) & !is.na(df$abs_GS) & df$abs_kME >= input$kme_thr & df$abs_GS >= input$gs_thr
      }
      df <- df[keep, , drop = FALSE]
      df[order(-df$abs_kME, -df$abs_GS), c("gene", "kME", "GS", "kWithin")]
    })

    output$kme_gs_plot <- renderPlot({
      df <- hub_data()
      ggplot(df, aes(x = kME, y = GS)) +
        geom_point(color = "#2c6fbb", alpha = 0.6, size = 1.6) +
        geom_vline(xintercept = c(-1, 1) * input$kme_thr, linetype = "dashed", color = "#c0392b") +
        geom_hline(yintercept = c(-1, 1) * input$gs_thr, linetype = "dashed", color = "#c0392b") +
        labs(x = sprintf("Module membership (kME) in %s", input$hub_module),
             y = sprintf("Gene significance for %s", input$hub_trait)) +
        theme_minimal(base_size = 12)
    })

    output$hub_table <- DT::renderDataTable({
      DT::datatable(hub_filtered(), rownames = FALSE, filter = "top",
                     options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_hubs <- downloadHandler(
      filename = function() sprintf("wgcna_hub_genes_%s.csv", input$hub_module),
      content = function(file) write.csv(hub_filtered(), file, row.names = FALSE)
    )

    ## Module membership vs intramodular connectivity - genes high on both are the most defensible hub calls.
    output$connectivity_plot <- renderPlot({
      df <- hub_data()
      ggplot(df, aes(x = kME, y = kWithin)) +
        geom_point(color = "#2c6fbb", alpha = 0.6, size = 1.6) +
        labs(x = sprintf("Module membership (kME) in %s", input$hub_module),
             y = "Intramodular connectivity (kWithin)") +
        theme_minimal(base_size = 12)
    })

    ## Module eigengene vs trait at sample level (Step 4's heatmap shown as a single r).
    output$me_trait_plot <- renderPlot({
      net <- net_result()
      req(input$hub_module, input$hub_trait)
      me_col <- paste0("ME", input$hub_module)
      raw <- net$meta[[input$hub_trait]]
      df <- data.frame(ME = net$MEs[[me_col]], trait = raw)
      if (!is.numeric(raw) && length(input$hub_trait_levels %||% character(0)) > 0) {
        df <- df[df$trait %in% input$hub_trait_levels, , drop = FALSE]
      }
      if (is.numeric(raw)) {
        ggplot(df, aes(x = trait, y = ME)) +
          geom_point(color = "#2c6fbb", alpha = 0.7) +
          geom_smooth(method = "lm", se = TRUE, color = "#c0392b", linewidth = 0.6) +
          labs(x = input$hub_trait, y = sprintf("Eigengene (%s)", input$hub_module)) +
          theme_minimal(base_size = 12)
      } else {
        df$trait <- factor(df$trait)
        ggplot(df, aes(x = trait, y = ME, fill = trait)) +
          geom_boxplot(outlier.size = 0.6) +
          scale_fill_manual(values = arthomix_pair(df$trait), guide = "none") +
          labs(x = input$hub_trait, y = sprintf("Eigengene (%s)", input$hub_module)) +
          theme_minimal(base_size = 12)
      }
    })

    ## Hub genes if >= 3 pass the thresholds, else the top 30 by |kME| in the module (so tight sliders don't leave the network blank).
    network_genes <- reactive({
      hf <- hub_filtered()$gene
      if (length(hf) >= 3) return(hf)
      df <- hub_data()
      df[order(-df$abs_kME), "gene"][seq_len(min(30, nrow(df)))]
    })

    output$network_status_ui <- renderUI({
      hf <- hub_filtered()$gene
      ng <- network_genes()
      if (length(hf) >= 3) {
        div(class = "empty-note", icon("check"), sprintf("Showing the %d hub genes that pass the kME/GS thresholds above.", length(hf)))
      } else {
        div(class = "empty-note", icon("circle-info"),
            sprintf("Fewer than 3 genes pass the kME/GS thresholds above, so showing the top %d genes by module membership (kME) in \"%s\" instead. Lower the thresholds to use your own filtered list.",
                    length(ng), input$hub_module))
      }
    })

    output$string_link_ui <- renderUI({
      genes <- network_genes()
      if (length(genes) == 0) {
        return(div(class = "empty-note", icon("circle-info"), "No genes to look up yet."))
      }
      tags$a(href = wgcna_string_url(genes), target = "_blank", rel = "noopener",
             class = "btn btn-default btn-sm", icon("diagram-project"), " Open in STRING-DB")
    })

    ## Live STRING-DB network image for these genes (capped at 60; hidden via onerror if the request fails).
    output$string_image_ui <- renderUI({
      genes <- network_genes()
      if (length(genes) < 2) {
        return(div(class = "empty-note", icon("circle-info"), "Need at least 2 genes to request a STRING-DB image."))
      }
      url <- wgcna_string_image_url(utils::head(genes, 60))
      tagList(
        tags$img(
          src = url, style = "max-width:100%; display:block; border:1px solid var(--color-border); border-radius:4px;",
          onerror = "this.style.display='none'; this.nextElementSibling.style.display='block';"
        ),
        div(class = "empty-note", style = "display:none;", icon("triangle-exclamation"),
            "Could not load the STRING-DB image right now - the external service may be unreachable. Try \"Open in STRING-DB\" above instead.")
      )
    })

    ## Builds Cytoscape edge/node tables for network_genes() via exportNetworkToCytoscape (in-memory, no file write).
    cyto_export <- reactive({
      net <- net_result()
      hub_genes <- network_genes()
      validate(need(length(hub_genes) >= 3, "This module has too few genes to build a network."))
      result <- tryCatch({
        sub_expr <- net$texpr[, hub_genes, drop = FALSE]
        adj <- WGCNA::adjacency(sub_expr, power = net$power, type = net$network_type,
                                 corFnc = wgcna_cor_fnc_name(net$cor_method))
        tom <- WGCNA::TOMsimilarity(adj, TOMType = net$tom_type, verbose = 0)
        dimnames(tom) <- dimnames(adj)
        WGCNA::exportNetworkToCytoscape(
          tom, edgeFile = NULL, nodeFile = NULL, weighted = TRUE,
          threshold = input$cyto_threshold, nodeNames = hub_genes,
          nodeAttr = unname(net$module_colors[hub_genes])
        )
      }, error = function(e) { message("WGCNA cyto_export failed: ", conditionMessage(e)); NULL })
      validate(need(!is.null(result), "Could not build the network for this gene set - try a different module or fewer genes."))
      result
    })

    output$download_cyto_edges <- downloadHandler(
      filename = function() sprintf("wgcna_cytoscape_edges_%s.csv", input$hub_module),
      content = function(file) write.csv(cyto_export()$edgeData, file, row.names = FALSE)
    )
    output$download_cyto_nodes <- downloadHandler(
      filename = function() sprintf("wgcna_cytoscape_nodes_%s.csv", input$hub_module),
      content = function(file) write.csv(cyto_export()$nodeData, file, row.names = FALSE)
    )

    ## In-app network figure from the same edge/node data as the Cytoscape export (nodes colored by
    ## module, edge width = TOM weight, force-directed layout); failures become a validate() message.
    network_plot_obj <- reactive({
      cyto <- cyto_export()
      edges <- cyto$edgeData
      validate(need(nrow(edges) > 0, "No edges pass the current weight threshold below - lower it to see a network."))
      nodes <- cyto$nodeData
      colnames(nodes)[3] <- "module"

      g <- tryCatch(
        igraph::graph_from_data_frame(
          edges[, c("fromNode", "toNode", "weight")], directed = FALSE,
          vertices = nodes[, c("nodeName", "module")]
        ),
        error = function(e) { message("WGCNA network_plot_obj: graph build failed: ", conditionMessage(e)); NULL }
      )
      validate(need(!is.null(g), "Could not build a graph from the current edges - try a different module or threshold."))

      ## geom_edge_link0() (not geom_edge_link()): avoids a namespace collision
      ## with Bioconductor packages loaded elsewhere in the app (IRanges/S4Vectors/BiocGenerics).
      p <- tryCatch(
        ggraph::ggraph(g, layout = "fr") +
          ggraph::geom_edge_link0(aes(width = weight), alpha = 0.3, colour = "#9aa2ac") +
          ggraph::geom_node_point(aes(color = module), size = 5) +
          ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 3, color = "#33393f") +
          ggplot2::scale_color_identity() +
          ggraph::scale_edge_width(range = c(0.3, 2), guide = "none") +
          ggplot2::theme_void() +
          ggplot2::labs(title = sprintf("%s module - %d genes, %d edges", input$hub_module, igraph::vcount(g), igraph::ecount(g))),
        error = function(e) { message("WGCNA network_plot_obj: plot build failed: ", conditionMessage(e)); NULL }
      )
      validate(need(!is.null(p), "Could not draw the network figure - try a different module or adjust the edge-weight slider."))
      p
    })

    ## Explicit print() catches draw-time failures too (ggraph's layout only runs when printed, not when built).
    output$network_plot <- renderPlot({
      p <- network_plot_obj()
      ok <- tryCatch({ print(p); TRUE },
                      error = function(e) { message("WGCNA network_plot draw failed: ", conditionMessage(e)); FALSE })
      validate(need(ok, "Could not draw the network figure for this gene set - try a different module or a higher edge-weight threshold."))
    })

    output$download_network_png <- downloadHandler(
      filename = function() sprintf("wgcna_network_%s.png", input$hub_module),
      content = function(file) {
        p <- network_plot_obj()
        tryCatch(
          ggsave(file, plot = p, width = 8, height = 7, dpi = 300, bg = "white"),
          error = function(e) {
            message("WGCNA download_network_png failed: ", conditionMessage(e))
            grDevices::png(file, width = 800, height = 700)
            graphics::plot.new()
            graphics::text(0.5, 0.5, "Could not render this network figure.")
            grDevices::dev.off()
          }
        )
      }
    )

    output$step5_ui <- renderUI({
      tagList(
        box(
          width = 12, title = tagList(icon("star"), " Module membership & gene significance"), status = "primary", solidHeader = FALSE,
          ## Defaults |kME|>0.8 AND |GS|>0.2 (Section 2.4.12): a hub must be both central to its module and disease-associated.
          p(class = "submodule-desc", "Hub genes are central to their module (high kME) and associated with your trait (high GS)."),
          uiOutput(ns("hub_picker_ui")),
          uiOutput(ns("hub_trait_levels_ui")),
          fluidRow(
            column(6, sliderInput(ns("kme_thr"), "Minimum |kME|", min = 0, max = 1, value = 0.8, step = 0.05)),
            column(6, sliderInput(ns("gs_thr"), "Minimum |GS|", min = 0, max = 1, value = 0.2, step = 0.05))
          ),
          actionButton(ns("run_hubs_btn"), "Compute hub genes", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        if (!wgcna_hubs_has_run()) {
          div(class = "empty-note", icon("circle-info"),
              "Not run yet. Pick a module and trait above, then click \"Compute hub genes\".")
        } else tagList(
          wgcna_result("star", "Hub gene diagnostics",
            wgcna_result_row(
              column(6, div(class = "wgcna-result-subtitle", "kME vs GS"),
                     withSpinner(plotOutput(ns("kme_gs_plot"), height = 300), color = "#2c6fbb", type = 6)),
              column(6, div(class = "wgcna-result-subtitle", "Hub genes"),
                     div(class = "table-toolbar", downloadButton(ns("download_hubs"), "Download CSV", class = "btn-sm")),
                     DT::dataTableOutput(ns("hub_table")))
            )
          ),
          wgcna_result_row(
            column(6, div(class = "wgcna-result-subtitle", "Intramodular connectivity"),
                   withSpinner(plotOutput(ns("connectivity_plot"), height = 280), color = "#2c6fbb", type = 6)),
            column(6, div(class = "wgcna-result-subtitle", "Eigengene vs trait"),
                   withSpinner(plotOutput(ns("me_trait_plot"), height = 280), color = "#2c6fbb", type = 6))
          ),
          wgcna_result("share-nodes", "Network figure & export",
            desc = "A co-expression network for the selected module - nodes colored by module, edges weighted by topological overlap. Adjust the slider to show fewer/more edges.",
            uiOutput(ns("network_status_ui")),
            sliderInput(ns("cyto_threshold"), "Minimum edge weight to include", min = 0, max = 0.5, value = 0.02, step = 0.01),
            withSpinner(plotOutput(ns("network_plot"), height = 460), color = "#2c6fbb", type = 6),
            div(class = "table-toolbar",
                downloadButton(ns("download_network_png"), "Network figure (PNG)", class = "btn-primary btn-sm"),
                downloadButton(ns("download_cyto_edges"), "Cytoscape edges (CSV)", class = "btn-sm"),
                downloadButton(ns("download_cyto_nodes"), "Cytoscape nodes (CSV)", class = "btn-sm"),
                uiOutput(ns("string_link_ui"), inline = TRUE))
          ),
          wgcna_result("image", "STRING-DB network image (live)",
            desc = "A real image fetched live from the STRING protein-protein interaction database for the same genes - independent evidence (known/predicted interactions), not this dataset's own co-expression.",
            withSpinner(uiOutput(ns("string_image_ui")), color = "#2c6fbb", type = 6)
          )
        )
      )
    })

    ## =====================================================================
    ## Step 6: Enrichment (GO/KEGG via clusterProfiler + org.Hs.eg.db, same as mod_enrichment.R, on a module's genes)
    ## =====================================================================

    output$enrich_hub_note <- renderUI({
      hf <- tryCatch(hub_filtered(), error = function(e) NULL)
      if (is.null(hf) || nrow(hf) == 0) {
        div(class = "empty-note", icon("circle-info"), "No hub genes yet - set thresholds in Step 5 first.")
      } else {
        div(class = "empty-note", icon("check"), sprintf("Using %d hub genes from module \"%s\" (Step 5).", nrow(hf), input$hub_module))
      }
    })

    output$enrich_picker_ui <- renderUI({
      net <- net_result()
      mods <- sort(setdiff(unique(net$module_colors), "grey"))
      validate(need(length(mods) > 0, "No non-grey modules were detected."))
      tagList(
        radioButtons(ns("enrich_gene_set"), "Gene set",
                     choiceNames = list("A whole module", "Hub genes from Step 5"),
                     choiceValues = list("module", "hub"), selected = "module"),
        conditionalPanel(condition = sprintf("input['%s'] == 'module'", ns("enrich_gene_set")),
                          selectInput(ns("enrich_module"), "Module", choices = mods, selectize = FALSE)),
        conditionalPanel(condition = sprintf("input['%s'] == 'hub'", ns("enrich_gene_set")),
                          uiOutput(ns("enrich_hub_note"))),
        radioButtons(ns("enrich_db"), "Database", choices = c("GO Biological Process" = "GO_BP", "KEGG pathways" = "KEGG"), selected = "GO_BP"),
        sliderInput(ns("enrich_qcut"), "Q-value cutoff", min = 0.01, max = 0.25, value = 0.05, step = 0.01),
        actionButton(ns("run_enrich_btn"), "Run enrichment", icon = icon("play"), class = "btn-primary btn-sm")
      )
    })

    ## reactiveValues + observeEvent (matching net_store/sft_store above) rather than
    ## eventReactive - enrichKEGG() in particular hits the live KEGG REST API, so a slow
    ## or unreachable network makes this run for a while; running it in an observer
    ## keeps that off the render path (so the controls box below isn't torn down/reset
    ## on every click) and lets any real error - network failure, a term-mapping issue -
    ## reach the UI instead of being swallowed by the `tryCatch(..., error = function(e)
    ## NULL)` guards used everywhere module_enrich() is read.
    enrich_store <- reactiveValues(result = NULL, error = NULL)

    observeEvent(input$run_enrich_btn, {
      enrich_store$error <- NULL
      result <- tryCatch({
        net <- net_result()
        use_hub <- identical(input$enrich_gene_set, "hub")
        gene_pool <- if (use_hub) {
          hf <- hub_filtered()
          validate(need(nrow(hf) > 0, "No hub genes available - set thresholds in Step 5 first."))
          hf$gene
        } else {
          net$gene_module$gene[net$gene_module$module == input$enrich_module]
        }
        validate(need(length(gene_pool) >= 5, "Need at least 5 genes to test enrichment. Pick a larger module or lower the Step 5 thresholds."))

        ## Background is the network's own gene universe (Step 1 output), not the whole dataset.
        universe_symbols <- colnames(net$texpr)
        map <- suppressMessages(AnnotationDbi::select(
          org.Hs.eg.db, keys = universe_symbols, keytype = "SYMBOL", columns = "ENTREZID"
        ))
        map <- map[!is.na(map$ENTREZID) & !duplicated(map$SYMBOL), ]
        gene_entrez <- map$ENTREZID[map$SYMBOL %in% gene_pool]
        validate(need(length(gene_entrez) >= 3, "Fewer than 3 of these genes could be mapped to Entrez IDs in the analyzed gene set."))

        ego <- if (identical(input$enrich_db, "GO_BP")) {
          clusterProfiler::enrichGO(
            gene = gene_entrez, universe = map$ENTREZID, OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
            pvalueCutoff = 1, qvalueCutoff = 1
          )
        } else {
          clusterProfiler::enrichKEGG(
            gene = gene_entrez, universe = map$ENTREZID, organism = "hsa",
            pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1
          )
        }
        df <- as.data.frame(ego)
        validate(need(nrow(df) > 0, "No enriched terms were found for this gene set."))
        label <- if (use_hub) sprintf("hub genes (%s)", input$hub_module) else input$enrich_module
        list(table = df, n_input = length(gene_pool), n_mapped = length(gene_entrez), label = label)
      }, error = function(e) e)
      if (inherits(result, "error")) {
        enrich_store$error <- conditionMessage(result)
      } else {
        enrich_store$result <- result
      }
    }, ignoreInit = TRUE)

    module_enrich <- reactive({
      req(enrich_store$result)
      enrich_store$result
    })

    output$enrich_summary_ui <- renderUI({
      ## step6_ui above only renders this once module_enrich() has a result (its own
      ## gate handles the error/not-run-yet states), so no need to re-check those here.
      res <- module_enrich()
      df <- res$table %>% dplyr::filter(qvalue < input$enrich_qcut)
      tagList(
        p(strong(res$n_mapped), " of ", res$n_input, " genes (", res$label, ") mapped to Entrez IDs and used in the test."),
        p(strong(nrow(df)), " terms significant at q < ", input$enrich_qcut, " (", nrow(res$table), " tested in total).")
      )
    })

    output$enrich_bar_plot <- renderPlot({
      res <- tryCatch(module_enrich(), error = function(e) NULL)
      req(res)
      df <- res$table %>% dplyr::filter(qvalue < input$enrich_qcut) %>% dplyr::arrange(qvalue) %>% utils::head(15)
      validate(need(nrow(df) > 0, "No terms pass the current q-value cutoff."))
      ggplot(df, aes(x = stats::reorder(Description, -log10(qvalue)), y = -log10(qvalue))) +
        geom_col(fill = "#2c6fbb") +
        coord_flip() +
        labs(x = NULL, y = "-log10 q-value") +
        theme_minimal(base_size = 12)
    })

    output$enrich_table <- DT::renderDataTable({
      res <- tryCatch(module_enrich(), error = function(e) NULL)
      req(res)
      DT::datatable(res$table, rownames = FALSE, filter = "top",
                     options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_enrich <- downloadHandler(
      filename = function() "wgcna_module_enrichment.csv",
      content = function(file) write.csv(module_enrich()$table, file, row.names = FALSE)
    )

    ## Controls only - not dependent on module_enrich()/enrich_store, so this box
    ## (including run_enrich_btn) doesn't get torn down and reset every time a run
    ## finishes. See the enrich_store observer above for why the actual enrichment
    ## call was pulled out of the render path.
    output$step6_controls_ui <- renderUI({
      box(
        width = 12, title = tagList(icon("flask"), " Module enrichment"), status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "Tests whether a module (or its hub genes) is enriched for GO terms or KEGG pathways, against the genes analyzed in Step 1 as background."),
        uiOutput(ns("enrich_picker_ui"))
      )
    })

    ## Gate + results only - re-renders when module_enrich()/enrich_store changes, but
    ## the controls box above (with the button itself) is a separate output and stays put.
    output$step6_ui <- renderUI({
      ## Gate: no result boxes until "Run enrichment" has been clicked - and a real
      ## failure (e.g. enrichKEGG's live KEGG API call failing) shows its own message
      ## rather than falling back to the same "Not run yet" text a fresh tab would show.
      if (!is.null(enrich_store$error)) {
        return(div(class = "empty-note", icon("triangle-exclamation"),
                    paste("Could not run enrichment:", enrich_store$error)))
      }
      if (is.null(tryCatch(module_enrich(), error = function(e) NULL))) {
        return(div(class = "empty-note", icon("circle-info"),
                    "Not run yet. Pick a module above, then click \"Run enrichment\" to see the results below."))
      }
      tagList(
        wgcna_result("chart-column", "Result",
          withSpinner(uiOutput(ns("enrich_summary_ui")), color = "#2c6fbb", type = 6),
          withSpinner(plotOutput(ns("enrich_bar_plot"), height = 380), color = "#2c6fbb", type = 6)
        ),
        wgcna_result("table-list", "Enriched terms",
          div(class = "table-toolbar", downloadButton(ns("download_enrich"), "Download CSV", class = "btn-sm")),
          DT::dataTableOutput(ns("enrich_table"))
        )
      )
    })

    ## =====================================================================
    ## Pipeline summary (shared vertical-timeline component, R/ui_shell.R)
    ## =====================================================================

    wgcna_progress <- reactive({
      power_ok <- !is.null(tryCatch(sft_result(), error = function(e) NULL))
      modules_ok <- !is.null(tryCatch(net_result(), error = function(e) NULL))
      traits_ok <- !is.null(tryCatch(module_trait(), error = function(e) NULL))
      hubs_ok <- !is.null(tryCatch(hub_filtered(), error = function(e) NULL))
      enrich_ok <- !is.null(tryCatch(module_enrich(), error = function(e) NULL))
      list(power_ok = power_ok, modules_ok = modules_ok, traits_ok = traits_ok, hubs_ok = hubs_ok, enrich_ok = enrich_ok)
    })

    output$pipeline_summary <- renderUI({
      pr <- wgcna_progress()
      step_state <- function(done, current) if (done) "done" else if (current) "current" else "future"
      steps <- list(
        list(number = 1, label = "Filter & QC", sublabel = "Genes & samples", state = "done"),
        list(number = 2, label = "Soft Power", sublabel = "Network power", state = step_state(pr$power_ok, TRUE)),
        list(number = 3, label = "Modules", sublabel = "Detect & cluster", state = step_state(pr$modules_ok, pr$power_ok && !pr$modules_ok)),
        list(number = 4, label = "Module-Trait", sublabel = "Correlate traits", state = step_state(pr$traits_ok, pr$modules_ok && !pr$traits_ok)),
        list(number = 5, label = "Hub Genes", sublabel = "Rank & export", state = step_state(pr$hubs_ok, pr$traits_ok && !pr$hubs_ok)),
        list(number = 6, label = "Enrichment", sublabel = "GO & KEGG", state = step_state(pr$enrich_ok, pr$hubs_ok && !pr$enrich_ok))
      )
      pipeline_summary_ui(steps)
    })

    ## Forces every step's UI to render on mount instead of only when its tab is first viewed
    ## (Shiny's default suspendWhenHidden = TRUE would leave later steps' inputs NULL until visited).
    for (step in paste0("step", 1:6, "_ui")) outputOptions(output, step, suspendWhenHidden = FALSE)
  })
}
