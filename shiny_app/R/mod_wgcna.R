## R/mod_wgcna.R
## Submodule: WGCNA Co-expression Network (Section 2.4)
## "Your analysis" runs a live, end-to-end WGCNA pipeline on the currently
## loaded dataset, as a 6-step wizard (same workflow-stepper pattern as
## Preprocessing, see mod_preprocessing.R):
##   1. Filter & QC     - which genes/samples go in (incl. a custom gene
##                         panel and a regex exclude filter), outlier removal
##   2. Soft Power       - network type, correlation method, power, and a
##                         scale-free topology check at the chosen power
##   3. Modules           - TOM/tree-cut parameters, dendrogram, module sizes,
##                         eigengene network (how modules relate to each other)
##   4. Module-Trait      - correlate module eigengenes with any trait
##   5. Hub Genes          - module membership, gene significance,
##                         intramodular connectivity, eigengene-vs-trait plot,
##                         Cytoscape export, STRING-DB lookup
##   6. Enrichment         - live GO/KEGG over-representation for a module or
##                         its hub genes (org.Hs.eg.db + clusterProfiler -
##                         the same databases mod_enrichment.R already uses)
##
## Wizard sits in an 8-wide column with a persistent 4-wide left rail
## alongside it - same outer split, same left/right order as
## mod_dge.R/mod_preprocessing.R (rail column first, main content second),
## not full page width. The rail holds only the "Ask ArthOChat" shortcut
## (see arthochat_shortcut_ui() in global.R), not a pipeline-summary panel
## like Preprocessing, so it stays short and doesn't compete for space with
## the step content. Every box/plot-grid width inside the 6 steps below is
## already relative (box(width = 12), column(4)/column(6) grids), so
## narrowing the wizard's outer column from 12 to 8 just makes them
## proportionally narrower - nothing clips, same shared drawer/session as
## everywhere else in the app, just a WGCNA-specific hint.
##
## The parameter set (gene-filtering method, network type/TOMType, deepSplit,
## mergeCutHeight, PAM stage, kME/GS-based hub-gene filtering, Cytoscape
## export) follows the standard WGCNA tutorial workflow (Langfelder P,
## Horvath S. WGCNA: an R package for weighted correlation network analysis.
## BMC Bioinformatics 2008;9:559; Zhang B, Horvath S. A general framework for
## weighted gene co-expression network analysis. Stat Appl Genet Mol Biol
## 2005;4:Article17) and mirrors the parameters exposed by xOmicsShiny's own
## WGCNA module (https://github.com/interactivereport/xOmicsShiny,
## wgcnamodule.R: top-N variable genes, mergeCutHeight, deepSplit = 2,
## pamRespectsDendro = FALSE, signed network type by default) - reimplemented
## here against this app's own dataset/results plumbing and UI system rather
## than reused directly, and extended with sample-outlier QC, a custom gene
## panel, a manual power override, a scale-free topology check, intramodular
## connectivity, module enrichment, a Cytoscape export and a STRING-DB
## lookup that xOmicsShiny's module doesn't expose.

mod_wgcna_config <- list(
  id = "wgcna", group = "Network",
  title = "WGCNA Co-expression Network", section = "Section 2.4",
  description = "Weighted gene co-expression network analysis, end to end: gene/sample filtering, soft-threshold power, module detection, module-trait correlation, hub genes, GO/KEGG enrichment and Cytoscape export.",
  icon = "circle-nodes"
)

## Trait values that aren't already numeric are recoded as.numeric(factor(.))
## - the same approach the original module used for `group`/`sex`. The sign
## of the resulting correlation is therefore arbitrary (whichever level sorts
## first alphabetically becomes 1).
##
## `levels_keep` (categorical traits only - ignored for numeric columns)
## restricts the comparison to specific levels, e.g. "RA" vs "Control" out of
## a group column that also has "other": samples whose value isn't in
## levels_keep become NA for this trait, and every cor(..., use = "p") call
## in this module already drops NAs pairwise, so that trait's correlation is
## computed only over the chosen levels without touching any other trait's
## sample set.
wgcna_encode_trait <- function(meta, col, levels_keep = NULL) {
  v <- meta[[col]]
  if (is.numeric(v)) return(as.numeric(v))
  v <- as.character(v)
  if (!is.null(levels_keep) && length(levels_keep) > 0) v[!(v %in% levels_keep)] <- NA
  as.numeric(factor(v, levels = sort(unique(stats::na.omit(v)))))
}

## The dynamic per-trait level-filter checkbox input is named
## "trait_levels_<col>" (Step 4) - NULL for a numeric column (no filter UI is
## shown for those) or before the checkboxes have rendered.
wgcna_trait_levels <- function(input, meta, col) {
  if (is.numeric(meta[[col]])) return(NULL)
  input[[paste0("trait_levels_", col)]]
}

## Same corFnc-selection idiom pickSoftThreshold/blockwiseModules use
## internally, so every correlation computed in this module (soft-threshold
## fitting, module detection, module membership, gene significance) respects
## the same "Pearson vs bicor" choice made in Step 2.
wgcna_cor_fnc <- function(cor_method) {
  if (identical(cor_method, "bicor")) WGCNA::bicor else WGCNA::cor
}
wgcna_cor_fnc_name <- function(cor_method) if (identical(cor_method, "bicor")) "bicor" else "cor"

## STRING-DB (string-db.org) multi-gene network lookup URL - same idea as
## global.R's own geo_link() for NCBI GEO accessions: a plain link builder
## for a well-known external bioinformatics resource, not a guessed URL.
## Species is fixed to human (9606), matching this app's RA blood cohort.
wgcna_string_url <- function(genes, species = 9606) {
  ids <- paste(vapply(genes, utils::URLencode, character(1), reserved = TRUE), collapse = "%0d")
  sprintf("https://string-db.org/cgi/network?identifiers=%s&species=%d", ids, species)
}

## STRING's REST image endpoint - a real, documented API that returns an
## actual rendered PNG network image for a gene list (confirmed directly:
## curl against this exact URL pattern returns a valid PNG). This is the one
## external tool in this app's toolchain that genuinely hosts a "figure from
## a database" the way it was asked for - Cytoscape itself is a desktop
## application with no such hosted API, so the in-app network_plot below
## (built from this run's own WGCNA edges) is what stands in for it.
wgcna_string_image_url <- function(genes, species = 9606) {
  ids <- paste(vapply(genes, utils::URLencode, character(1), reserved = TRUE), collapse = "%0d")
  sprintf("https://string-db.org/api/image/network?identifiers=%s&species=%d", ids, species)
}

## Loads this project's own already-run WGCNA result straight from disk -
## the actual output of scripts/00_shared/06_WGCNA.R on the real 173-
## sample training network (data/processed/wgcna_results.rds: datExpr,
## meta, moduleColors, MEs, soft_power, config - already-labelled colours,
## not numeric), not a re-derivation inside the app. Instant, in place of
## the several-minutes-to-tens-of-minutes live "Detect modules" run on the
## full 15,763-gene network. wgcna_results.rds alone is enough for every
## downstream step except the dendrogram plot, which needs the raw
## blockwiseModules() object (wgcna_net_<hash>.rds, hash changes if the
## script is ever rerun with different settings, hence the glob instead of
## a hardcoded filename) - falls back to no dendrogram rather than failing
## the whole load if that file isn't present.
load_precomputed_wgcna_result <- function() {
  proc_dir <- file.path(DATA_ROOT, "data", "processed")
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
    ## Precomputed kME/GS_RA/connectivity/is_hub, disease modules only
    ## (yellow+brown - see hub_data() in Step 5 below) - the actual output
    ## of the real pipeline's hub-gene rule, not a live re-derivation. Used
    ## whenever Step 5's module+trait selection exactly matches what this
    ## was computed against (disease module, "group" trait, RA vs HC, no
    ## narrower level filter); falls back to live computation otherwise,
    ## e.g. for any other module or trait.
    hub_table_precomputed = res$hub_table
  )
}

## Step-title markup for the workflow-stepper nav - identical pattern to
## mod_preprocessing.R's pp_step_title(), duplicated (not shared) so each
## module file stays self-contained the way this codebase already keeps
## mod_<id>.R files independent of one another. Labels are kept short (one
## or two words) since this stepper has 6 items, not 3.
wgcna_step_title <- function(number, ic, label, sublabel) {
  tagList(
    span(class = "step-number", as.character(number)),
    icon(ic),
    span(class = "step-text",
         span(class = "step-label", label),
         span(class = "step-sublabel", sublabel))
  )
}

## ---------------------------------------------------------------------------
## Results section - a boxless alternative to box() for RESULTS/PLOTS only
## (settings/config panels in each step keep their existing box() treatment
## unchanged, per the brief this replaces: "remove the individual boxes/cards
## around results and plots", "keep the existing layout everywhere else").
## A subtle top divider + heading instead of a bordered card; wgcna_result_row()
## below adds a light vertical divider between side-by-side result columns so
## multi-plot rows still read as grouped without a box around each one.
## ---------------------------------------------------------------------------

wgcna_result <- function(icon_name, title, ..., desc = NULL) {
  tagList(
    div(class = "wgcna-result-heading", icon(icon_name), title),
    if (!is.null(desc)) p(class = "submodule-desc wgcna-result-desc", desc),
    ...
  )
}

## Wraps a fluidRow of same-level result columns so the shared CSS below can
## draw a subtle vertical divider between them (desktop widths only - see
## wgcna_layout_css()) instead of each column getting its own box border.
wgcna_result_row <- function(...) {
  div(class = "wgcna-result-row", fluidRow(...))
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------

## Scoped, WGCNA-only equal-height fix: shinydashboard box()es in the same
## fluidRow don't stretch to match each other's height by default (a plain
## Bootstrap row/column doesn't do that on its own), so a row pairing a
## short box next to a taller one - e.g. Step 2's three power plots, Step
## 3's dendrogram/eigengene-network pair, Step 5's plot/table pairs - looks
## visibly misaligned along the bottom edge. This makes every box in a
## `.wgcna-wrap` row stretch to the row's tallest box instead. Scoped to
## `.wgcna-wrap` (not the shared `.workflow-stepper-wrap` class Preprocessing
## also uses) so it can't change any other module's layout.
wgcna_layout_css <- function() {
  tags$style(HTML("
    .wgcna-wrap .row { display: flex; flex-wrap: wrap; }
    .wgcna-wrap .row > div[class*='col-'] { display: flex; flex-direction: column; }
    .wgcna-wrap .row > div[class*='col-'] > .box { flex: 1 1 auto; display: flex; flex-direction: column; margin-bottom: 0; }
    .wgcna-wrap .row > div[class*='col-'] > .box > .box-body { flex: 1 1 auto; }

    /* The shared workflow-stepper tab bar was designed for Preprocessing's
       3 steps; 6 steps in the same space only gave the ACTIVE tab a visible
       border/box while the other 5 floated as bare text, reading as
       unbalanced. Give every tab the same card treatment, tighten spacing
       to fit 6 comfortably, and drop the connector dashes (tuned for 3
       evenly-spaced items, not 6 narrower ones).
       `min-width: 0` below is load-bearing, not decorative: flex items
       default to `min-width: auto`, meaning a flex-shrink item still
       refuses to shrink below its own content's natural width - adding the
       border/padding above made that natural width bigger, so the row
       overflowed past the box's right edge instead of shrinking to fit.
       `min-width: 0` removes that floor, letting each tab shrink to the
       row's available width instead of overflowing it. Text wraps onto a
       second line at that narrower width (`white-space: normal`) rather
       than being truncated with an ellipsis - full label/sublabel text
       stays readable; the tab item just gets a little taller instead of
       cutting words off. `align-items: flex-start` keeps the icon/number
       lined up with the first line of text once it wraps, instead of
       re-centering against the now-taller two-line block.
       flex-wrap: wrap (the shared .workflow-stepper-wrap rule is nowrap,
       tuned for 3 items) plus a real min-width floor on each li: below
       roughly 3 tabs' worth of row width, shrinking all 6 to fit one row
       stops being legible no matter how small the font gets - text sits
       flush against the icon/number with no breathing room, reading as
       overlapping even though nothing literally overlaps in the DOM.
       Wrapping to 2 rows of 3 once that floor is hit keeps every tab at a
       comfortable, constant width instead. */
    .wgcna-wrap .nav-tabs { flex-wrap: wrap; }
    .wgcna-wrap .nav-tabs > li:not(:last-child)::after { display: none; }
    .wgcna-wrap .nav-tabs > li { min-width: 140px; flex: 1 1 30%; }
    .wgcna-wrap .nav-tabs > li > a {
      border: 1px solid transparent !important;
      border-radius: 8px !important;
      padding: 7px 8px !important;
      gap: 6px !important;
      min-width: 0;
      align-items: flex-start !important;
    }
    .wgcna-wrap .nav-tabs > li > a:hover { background: #F5F7FA !important; }
    .wgcna-wrap .nav-tabs > li.active > a { border-color: var(--color-border) !important; }
    .wgcna-wrap .step-number { width: 20px; height: 20px; font-size: 10.5px; flex-shrink: 0; margin-top: 1px; }
    .wgcna-wrap .step-text { min-width: 0; }
    /* overflow-wrap: normal (not break-word) - wrapping is only ever
       allowed at a real word boundary (a space, or the hyphen in
       Module-Trait), never mid-word, so a single word like Modules
       never gets split into two broken fragments. Smaller text below is
       what actually buys the room to keep whole words on one line. */
    .wgcna-wrap .step-label, .wgcna-wrap .step-sublabel { white-space: normal; overflow-wrap: normal; word-break: normal; }
    .wgcna-wrap .step-label { font-size: 10.5px; line-height: 1.25; }
    .wgcna-wrap .step-sublabel { font-size: 9px; line-height: 1.25; }

    /* ---------- Results section (boxless): consistent grid, equal margins,
       subtle dividers instead of box() cards - settings panels above each
       results section keep their existing box() look, untouched. ---------- */
    .wgcna-wrap .wgcna-result-heading {
      display: flex; align-items: center; gap: 8px;
      font-size: 15px; font-weight: 600; color: var(--color-ink);
      margin: 28px 0 6px 0; padding-top: 20px;
      border-top: 1px solid var(--color-border);
      /* box() (shinydashboard) wraps itself in a floated <div class=
         'col-sm-12'> even outside a .row - without a .row's own clearfix
         after it, that float is never cleared, so a plain flex/block
         sibling right after a settings box() tries to render beside it
         instead of below it and gets squeezed into the (zero) remaining
         width. clear: both forces this heading below any such float. */
      clear: both;
    }
    /* First results section right under a settings box needs no extra
       top gap of its own - the box already ends with clear margin below it. */
    .wgcna-wrap .box + .wgcna-result-heading,
    .wgcna-wrap .box + div > .wgcna-result-heading:first-child {
      margin-top: 20px;
    }
    .wgcna-wrap .wgcna-result-heading svg, .wgcna-wrap .wgcna-result-heading .fa { color: var(--color-primary); }
    .wgcna-wrap .wgcna-result-desc { margin-bottom: 12px; }
    /* Equal margins/gutters for every result plot/table/output - same
       principle box-body padding gave each card, applied consistently
       instead of per-box. */
    .wgcna-wrap .wgcna-result-row { margin-bottom: 4px; clear: both; }
    .wgcna-wrap .wgcna-result-row .row { row-gap: 20px; }
    @media (min-width: 900px) {
      .wgcna-wrap .wgcna-result-row .row > div[class*='col-']:not(:first-child) {
        border-left: 1px solid var(--color-border);
      }
    }
    /* Small column-level label inside a wgcna-result-row (was a box()
       title before) - lighter weight than the section heading above it,
       since it's identifying one plot within an already-titled group. */
    .wgcna-wrap .wgcna-result-subtitle {
      font-size: 13px; font-weight: 600; color: var(--color-ink);
      margin-bottom: 8px;
    }
  "))
}

## Same left-rail pattern mod_dge.R uses: ArthOChat stacked with a settings/
## status card in column(3), main content in column(9) beside it (Preprocessing's
## column(4)/column(8) split is too narrow here - by the time this nests inside
## the app's own outer sidebar column, an 8/12 share of that for a 6-step
## stepper plus its result cards left barely 425px to work with; 9/12 buys
## back some of that without touching Preprocessing or any other module,
## since this split is local to mod_wgcna_ui).
mod_wgcna_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      3,
      div(
        class = "pipeline-rail",
        arthochat_shortcut_ui(
          "Questions about WGCNA? Ask ArthOChat.",
          compact = TRUE
        ),
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
            value = "Filter", title = wgcna_step_title(1, "filter", "Filter & QC", "Genes & samples"),
            br(), uiOutput(ns("step1_ui"))
          ),
          tabPanel(
            value = "Power", title = wgcna_step_title(2, "wave-square", "Soft Power", "Network power"),
            br(), uiOutput(ns("step2_ui"))
          ),
          tabPanel(
            value = "Modules", title = wgcna_step_title(3, "diagram-project", "Modules", "Detect & cluster"),
            br(), uiOutput(ns("step3_ui"))
          ),
          tabPanel(
            value = "Traits", title = wgcna_step_title(4, "table-cells", "Module-Trait", "Correlate traits"),
            br(), uiOutput(ns("step4_ui"))
          ),
          tabPanel(
            value = "Hubs", title = wgcna_step_title(5, "star", "Hub Genes", "Rank & export"),
            br(), uiOutput(ns("step5_ui"))
          ),
          tabPanel(
            value = "Enrichment", title = wgcna_step_title(6, "flask", "Enrichment", "GO & KEGG"),
            br(), uiOutput(ns("step6_ui"))
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

    ## True as soon as the active dataset came from the user's own data,
    ## whether loaded straight off the "Upload your own data" form on the
    ## Dataset tab (mod_dataset.R sets dataset$source to "Uploaded dataset:
    ## <files>", as-is, no normalisation/batch correction) or run through
    ## Preprocessing and Batch Correction first and then activated
    ## (mod_preprocessing.R's activate_btn sets "Uploaded dataset
    ## (preprocessed + normalised/batch-corrected): <files>") - both start
    ## with "Uploaded dataset", just without a colon in the second case, so
    ## match on the shared prefix rather than the punctuation. Drives the
    ## paper-methodology defaults below. Referenced (not just called) from
    ## inside step1_ui/step2_ui/step3_ui's renderUI so those steps
    ## re-render, and so re-default, every time a new dataset is loaded -
    ## including switching from one upload to another, or from an upload
    ## back to a preloaded/example dataset.
    dataset_is_uploaded <- reactive({
      isTRUE(grepl("^Uploaded dataset", dataset$source %||% ""))
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

    ## Six gene-selection strategies (the original module only offered the
    ## first): a hard top-N cutoff by variance (matches xOmicsShiny's
    ## "Select Top N Genes"), variance/MAD percentile cutoffs that scale with
    ## the dataset instead of a fixed count, a mean-expression floor, a
    ## user-pasted candidate gene panel (targeted WGCNA), or no filtering at
    ## all. A regex exclude filter then applies on top of any of these.
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
        ## Matches the published methodology this module quotes verbatim in
        ## its "paper defaults" note below: "The 5000 genes exhibiting the
        ## highest median expression values were then incorporated into a
        ## WGCNA analysis." Median, not mean/variance - a deliberate
        ## different ranking from every other method here.
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

    ## Sample clustering runs on the gene-filtered matrix (the same input
    ## that will go into the network), so the dendrogram the user sees is
    ## exactly what outlier removal will act on.
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

    ## Standard WGCNA tutorial outlier rule (Zhang & Horvath 2005 / the
    ## Langfelder-Horvath tutorial's own sample-clustering step): cut the
    ## sample tree at the chosen height and keep only the largest resulting
    ## cluster; everything else is flagged as an outlier and dropped.
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
      meta <- gs$meta[match(rownames(texpr_full), gs$meta$sample), , drop = FALSE]
      list(texpr = texpr_full, meta = meta,
           n_samples_before = ncol(gs$expr), n_samples_after = nrow(texpr_full),
           removed_samples = removed)
    })

    output$gene_filter_summary <- renderUI({
      gs <- gene_selection()
      div(class = "empty-note", icon("check"),
          sprintf("%s of %s genes kept.", format(gs$n_kept, big.mark = ","), format(gs$n_total, big.mark = ",")))
    })

    output$sample_qc_summary <- renderUI({
      fi <- final_input()
      if (length(fi$removed_samples) > 0) {
        div(class = "empty-note", icon("triangle-exclamation"),
            sprintf("%d of %d samples kept. Removed as outliers: %s.",
                    fi$n_samples_after, fi$n_samples_before, paste(fi$removed_samples, collapse = ", ")))
      } else {
        div(class = "empty-note", icon("check"),
            sprintf("%d of %d samples kept.", fi$n_samples_after, fi$n_samples_before))
      }
    })

    output$step1_ui <- renderUI({
      ## Re-evaluated (and so re-defaulted) every time the active dataset
      ## changes - see dataset_is_uploaded() above.
      uploaded <- dataset_is_uploaded()
      tagList(
        box(
          width = 12, title = tagList(icon("dna"), " Gene filtering"), status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Choose which genes go into the network - fewer, more variable genes run faster and add less noise."),
          if (uploaded) div(
            class = "empty-note", style = "margin-bottom: 10px;", icon("book"),
            "You uploaded your own dataset, so the settings below (and in Soft Power / Modules) are pre-filled to match the published methodology: top 5,000 genes by highest median expression, power β = 7, minimum module size 66, merge cut height 0.3. Adjust any of them before you run."
          ),
          ## Defaults to "Use all genes" for this project's own dataset (it
          ## deliberately skipped any variance pre-filter - Chapter_2_
          ## subchapter2_sexstratified.md Section 2.4.1 - so a gene removed
          ## before module detection can never be assigned to a module and
          ## so can never enter the disease-module intersection candidate
          ## genes are drawn from later, Section 2.5). For an uploaded
          ## dataset it instead defaults to the median-expression top-N the
          ## published methodology above describes. Every other filter
          ## method here is still fully available either way.
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

    sft_result <- eventReactive(input$compute_power_btn, {
      texpr <- final_input()$texpr
      powers <- c(1:10, seq(12, 20, 2))
      corName <- wgcna_cor_fnc_name(input$cor_method)
      ## Cached the same way net_result() below caches blockwiseModules -
      ## pickSoftThreshold tests 12 separate power values across the full
      ## gene set, comparably expensive at this project's own "all genes"
      ## default, so a repeat call with the same matrix/settings is instant
      ## instead of re-paying that cost.
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
          ## Actual degree-distribution check at the chosen power (the
          ## diagnostic the WGCNA paper itself uses to define
          ## "scale-free"), not just the fitted-curve summary above.
          k <- WGCNA::softConnectivity(
            texpr, corFnc = corName, corOptions = "use = 'p'",
            type = input$network_type, power = power, verbose = 0
          )
          list(sft_df = sft$fitIndices, power = power, connectivity = k)
        }
      )
    }, ignoreInit = TRUE)

    effective_power <- reactive({
      if (identical(input$power_mode, "manual")) {
        req(input$manual_power)
        as.numeric(input$manual_power)
      } else {
        sft_result()$power
      }
    })

    output$power_status_ui <- renderUI({
      res <- tryCatch(sft_result(), error = function(e) NULL)
      if (is.null(res)) {
        return(div(class = "empty-note", icon("circle-info"), "Not computed yet. Click \"Compute power\"."))
      }
      eff <- effective_power()
      div(class = "empty-note", icon("check"),
          sprintf("Estimated power: %s. Using power = %s (%s).", res$power, eff,
                  if (identical(input$power_mode, "manual")) "manual override" else "automatic"))
    })

    ## expression(R^2) draws a true mathematical superscript (matches base R
    ## plotmath rendering) instead of the Unicode "²" glyph, which some
    ## fonts/graphics devices render as a small inline "2" rather than a
    ## proper superscript - "correct format" per feedback on this label.
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

    ## Third power diagnostic: the truncated-exponential fit's own R² per
    ## power (WGCNA::pickSoftThreshold's fitIndices$truncated.R.sq column -
    ## an alternative goodness-of-fit already computed alongside the linear
    ## SFT.R.sq above, per WGCNA's own documentation, not a new statistic).
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

    ## Same underlying fit as WGCNA::scaleFreePlot() (same binning, same
    ## linear + exponential models - see WGCNA::scaleFreePlot source), but
    ## drawn with a direct plot()/lines() call instead of calling that
    ## function: scaleFreePlot ALWAYS appends its own "scale free R^2=...,
    ## slope=..., trunc.R^2=..." text to the plot's title, even when
    ## main = "" is passed in - there's no argument that removes it. Doing
    ## the drawing here instead keeps the same classic look (points, a black
    ## fit line, a red fit line) with a genuinely empty title, and the R²
    ## numbers are shown once, properly, in the write-up below. `k` is
    ## filtered to finite values first: WGCNA::softConnectivity() can return
    ## NA for genes with too few valid pairwise samples, and an NA in
    ## min()/max() (no na.rm) would otherwise silently break the binning.
    scale_free_check_stats <- reactive({
      res <- tryCatch(sft_result(), error = function(e) NULL)
      req(res)
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

    ## Draw-time errors (e.g. from a plotting library's own layout/repel
    ## step) only surface when the plot is actually printed, not when it's
    ## built - print() is called explicitly here, inside this tryCatch, so
    ## any such failure still becomes a plain validate() message instead of
    ## an unreadable raw error object.
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

    output$step2_ui <- renderUI({
      ## Re-evaluated (and so re-defaulted) every time the active dataset
      ## changes - see dataset_is_uploaded() above.
      uploaded <- dataset_is_uploaded()
      tagList(
        box(
          width = 12, title = tagList(icon("wave-square"), " Network settings"), status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "WGCNA raises gene-gene correlations to a power so the network is approximately scale-free. Pick the lowest power that reaches a good fit below."),
          fluidRow(
            column(
              6,
              ## Defaults to Unsigned for an uploaded dataset, Signed for
              ## this project's own: the published methodology's WGCNA
              ## paragraph states only pickSoftThreshold, minModuleSize = 66
              ## and mergeCutHeight = 0.3 - networkType is never mentioned,
              ## so nothing in the text overrides blockwiseModules()'s own
              ## package default of "unsigned". This project's own analysis
              ## chose signed deliberately and documents why at length
              ## (METHODS_2.4_WGCNA_expanded.md Section 2.4.7 - an unsigned
              ## network merges anti-correlated genes into one module with
              ## no coherent direction); an unstated parameter in someone
              ## else's paper is not evidence they made that same choice, so
              ## it isn't carried over to the "paper defaults" preset.
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
              ## 0.9, not a looser 0.8: this project's own analysis required
              ## R^2 >= 0.9 jointly with a connectivity ceiling before
              ## accepting a power (METHODS_2.4_WGCNA_expanded.md, Section
              ## 2.4.8), and rejected pickSoftThreshold's own auto-estimate
              ## for landing on a lower-R^2 power with runaway connectivity.
              sliderInput(ns("r_sq_cutoff"), HTML("Target scale-free fit (R<sup>2</sup>)"), min = 0.7, max = 0.95, value = 0.9, step = 0.01),
              ## Defaults to Manual, power 12, not Automatic, for this
              ## project's own dataset: that analysis explicitly rejected
              ## pickSoftThreshold's own auto-estimate (it landed on a
              ## lower-R^2 power with runaway connectivity) and chose beta =
              ## 12 by hand on the joint scale-free-fit/connectivity
              ## criterion (METHODS_2.4_WGCNA_expanded.md; Chapter_2_
              ## subchapter2_sexstratified.md Section 2.4.5). For an
              ## uploaded dataset it instead defaults to manual power 7, the
              ## published methodology's own selected soft-thresholding
              ## power. Still a free choice either way - switch back to
              ## Automatic, or type a different manual value, any time.
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
          div(style = "margin-top:8px;", uiOutput(ns("power_status_ui")))
        ),
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

    ## One button, one action, routed internally by which dataset is
    ## active - not exposed to the user as a "precomputed vs. live"
    ## choice. dataset$source is set once at app start to this project's
    ## own default (load_default_dataset() in global.R) and changes only
    ## when the user explicitly activates their own uploaded/preprocessed
    ## data (mod_preprocessing.R's "Use this as the active dataset"), so
    ## it's the same signal "Select data for this pipeline" (example vs.
    ## your own data) already drives elsewhere in the app. On the default
    ## dataset, "Run" loads this project's own real analysis instantly;
    ## on anyone's own data, "Run" computes live instead, with every
    ## filter/setting above (gene selection method, power mode, network
    ## type, min module size, etc.) fully available either way.
net_store <- reactiveValues(result = NULL, source = NULL, error = NULL)

    ## The whole body runs inside tryCatch so a validate()/need() failure
    ## (e.g. an unavailable power) or any other error is CAUGHT AND STORED
    ## instead of silently vanishing - this is a plain observeEvent, not a
    ## render context, so an uncaught validate() here previously showed no
    ## message anywhere at all (see the comment on suspendWhenHidden above,
    ## which fixes the actual root cause; this is the belt to that braces,
    ## for any other way "Run" could fail this same silent way).
    observeEvent(input$run_btn, {
      net_store$error <- NULL
      result <- tryCatch({
      if (isTRUE(grepl("^Example dataset:", dataset$source %||% ""))) {
        net_store$result <- load_precomputed_wgcna_result()
        net_store$source <- "loaded"
        return(invisible(NULL))
      }

      fi <- final_input()
      texpr <- fi$texpr
      power <- tryCatch(effective_power(), error = function(e) NA)
      validate(need(!is.null(power) && length(power) == 1 && is.finite(power),
                    "Compute the network power in Step 2 first, or switch to a manual power there."))

      ## Cached (global.R's get_or_compute_wgcna_blocks()) on a digest of
      ## the exact expression matrix plus every setting below - the same
      ## combination of dataset/filter/power/network settings only ever
      ## pays this cost once, on this machine, for good; a different
      ## combination just computes and caches separately. Nothing about
      ## which settings are choosable changes.
      ## reassignThreshold: 0 (disabled) for this project's own dataset, a
      ## deliberate choice (METHODS_2.4_WGCNA_expanded.md Table 2.x - keeps
      ## module membership determined by topological overlap alone). For an
      ## uploaded dataset, the published methodology's WGCNA paragraph never
      ## mentions this parameter either, so it's left at blockwiseModules()'s
      ## own package default (1e-06) instead of silently inheriting this
      ## project's own unrelated choice.
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
          ## Fixed, not user-exposed - matches this project's own
          ## 06_WGCNA.R (CFG$seed = 1234, passed through as
          ## blockwiseModules(randomSeed = CFG$seed)). The PAM-stage
          ## reassignment (active here since pamRespectsDendro = FALSE)
          ## has a stochastic tie-break; without pinning the same seed,
          ## an otherwise-identical run can land on different module
          ## boundaries, and because colour names are assigned by module
          ## SIZE RANK (not a stable identity), even a small boundary
          ## shift can relabel which colour is "yellow" vs "brown" vs
          ## anything else entirely.
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

    ## base <- results$wgcna %||% list() preserves significant_trait_modules/
    ## traits_tested from a PRIOR run of Step 4 below rather than wiping them
    ## on every re-run of Step 3 - the same accumulate-into-one-list pattern
    ## Step 4's own observer already uses. module_genes (module color -> gene
    ## vector, from this run's gene_module assignment) is what lets the
    ## Candidate Gene Identification tab offer "this module's genes" as a
    ## click-based pick instead of a pasted list.
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

    ## ---------------------------------------------------------------------
    ## Gene-network heatmap (TOM) - uploaded dataset only. The classic WGCNA
    ## network heatmap ("interactions between genes in coexpression
    ## modules" - Langfelder & Horvath tutorial, TOMplot()), not the
    ## module-level eigengene-network heatmap already shown below (that's a
    ## different, much smaller N-module x N-module summary). A literal
    ## gene x gene TOM at thousands of genes is tens of millions of cells -
    ## infeasible to render directly - so this follows the tutorial's own
    ## standard approach: compute the TOM once on the full filtered gene
    ## set, then randomly subsample a fixed number of genes for the actual
    ## plot. Gated behind its own button (like Step 2's "Compute power")
    ## rather than running on every Step 3 "Run", since
    ## TOMsimilarityFromExpr() on thousands of genes is comparably
    ## expensive to blockwiseModules() itself.
    ## ---------------------------------------------------------------------

    tom_result <- eventReactive(input$compute_tom_btn, {
      net <- net_result()
      ## Cached the same way net_result()/sft_result() are (global.R's
      ## get_or_compute_wgcna_blocks()) - a repeat call on the same
      ## expression matrix/power/network settings is instant.
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
      ## Fixed seed (matches this module's other fixed-seed choices, e.g.
      ## blockwiseModules's randomSeed = 1234 above) - which genes land in
      ## the illustrative subsample doesn't affect any reported statistic,
      ## but pinning it means the plot doesn't reshuffle on every re-render.
      set.seed(1234)
      select <- sample(seq_len(n_genes), size = n_select)
      diss_tom <- 1 - full_tom
      select_tom <- diss_tom[select, select]
      select_tree <- stats::hclust(stats::as.dist(select_tom), method = "average")
      select_colors <- net$module_colors[select]
      ## ^7 and the NA diagonal match the tutorial's own TOMplot() recipe -
      ## raises contrast between strong/weak overlap and blanks the (always
      ## self-similar) diagonal so it doesn't dominate the color scale.
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
      ## Cell text/axis size scales down as module count grows (this
      ## project's own default run has 12+grey = 13) - at the previous
      ## fixed size = 2.7, cells this small overlapped both each other's
      ## correlation numbers and the axis labels. coord_fixed() keeps
      ## cells square instead of stretching to fill the panel.
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

    ## Full net_result() list (module colors, MEs, dendrogram, gene-module
    ## table, power/settings used) as a single .rds - everything this run
    ## produced, not just the flattened CSV above, for anyone who wants to
    ## keep working on it directly in R (e.g. re-drawing the dendrogram
    ## with different graphical parameters, or feeding MEs into their own
    ## downstream script) instead of only what's exported as a table/plot
    ## in this tab. This is the same object get_or_compute_wgcna_blocks()
    ## already caches to data/cache/wgcna/ on disk for reuse inside the
    ## app - this button just lets the user take a copy of it with them.
    output$download_wgcna_rds <- downloadHandler(
      filename = function() "wgcna_module_detection_result.rds",
      content = function(file) saveRDS(net_result(), file)
    )

    output$step3_ui <- renderUI({
      ## Re-evaluated (and so re-defaulted) every time the active dataset
      ## changes - see dataset_is_uploaded() above. minModuleSize = 66 and
      ## mergeCutHeight = 0.3 are the published methodology's own values
      ## ("minmodulesize = 66; mergecutheight = 0.3"); 30/0.25 remain this
      ## project's own choice for its own dataset. pamRespectsDendro is
      ## never mentioned in that same methodology paragraph, so for an
      ## uploaded dataset this defaults to TRUE, blockwiseModules()'s own
      ## package default - not this project's own deliberate FALSE (chosen
      ## to match the WGCNA authors' tutorial, METHODS_2.4_WGCNA_expanded.md
      ## Table 2.x), which has no basis in someone else's unstated methods.
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
          ## One button. What it does is routed internally by which
          ## dataset is active (see the run_btn observer above) - on this
          ## project's own default dataset it's instant, on anyone's own
          ## uploaded/preprocessed data it computes live, using every
          ## filter above (including "Use all genes" or a faster top-N
          ## variable-gene subset in Step 1). The time-estimate note only
          ## shows for the live-compute case, since the default dataset
          ## never needs it.
          actionButton(ns("run_btn"), "Run", icon = icon("play"), class = "btn-primary btn-sm"),
          if (!isTRUE(grepl("^Example dataset:", dataset$source %||% ""))) p(class = "submodule-desc", style = "margin-top:6px;",
            "Can take several minutes to tens of minutes depending on gene count and dataset size - narrowing Step 1's gene filter (e.g. the top 2,000-4,000 most variable genes instead of all genes) runs much faster for quick exploration."),
          div(style = "margin-top:8px;", uiOutput(ns("module_run_status_ui")))
        ),
        wgcna_result("diagram-project", "Gene dendrogram & modules",
          withSpinner(plotOutput(ns("dendro_plot"), height = 400), color = "#2c6fbb", type = 6),
          div(class = "table-toolbar", downloadButton(ns("download_dendro_png"), "Download dendrogram (PNG)", class = "btn-sm"))
        ),
        ## Uploaded dataset only - the classic WGCNA network heatmap
        ## ("interactions between genes in coexpression modules"), distinct
        ## from the module-level "Eigengene network" plot below. Its own
        ## opt-in button (not run automatically with the "Run" button above)
        ## since recomputing the full TOM is comparably expensive to module
        ## detection itself and not every user needs this panel.
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

    ## One "which levels?" checkbox row per selected categorical trait (not
    ## shown for numeric traits, e.g. age) - lets a group column with 3+
    ## levels (e.g. RA / HC / other) be narrowed to just the two you want to
    ## compare, e.g. RA vs HC, instead of forcing every level into one
    ## arbitrary numeric scale.
    output$trait_level_filters_ui <- renderUI({
      req(input$trait_cols)
      meta <- qc()$meta
      cat_traits <- Filter(function(cl) !is.numeric(meta[[cl]]), input$trait_cols)
      if (length(cat_traits) == 0) return(NULL)
      tagList(
        p(class = "submodule-desc", "Each level you keep below becomes its own column in the heatmap (e.g. group -> \"RA\" and \"HC\" columns) - untick a level (e.g. \"other\") to drop it from the comparison."),
        lapply(cat_traits, function(cl) {
          lv <- sort(unique(stats::na.omit(as.character(meta[[cl]]))))
          checkboxGroupInput(ns(paste0("trait_levels_", cl)), sprintf("\"%s\" levels to include", cl),
                              choices = lv, selected = lv, inline = TRUE)
        })
      )
    })

    ## Each SELECTED LEVEL of a categorical trait becomes its own binary
    ## (0/1) column - e.g. picking "group" shows as "RA" and "HC" columns in
    ## the heatmap, not one arbitrarily-numbered "group" column - so the sign
    ## of every correlation has an unambiguous, labeled meaning ("this module
    ## goes up in RA" vs "goes up in HC"), matching how the WGCNA tutorial
    ## itself recommends handling categorical traits. Samples outside the
    ## levels chosen in the "levels to include" checkboxes above are NA in
    ## every one of that trait's columns (excluded, not coded as 0), so e.g.
    ## an unchecked "other" group doesn't get silently folded into either
    ## remaining indicator. Numeric traits (e.g. age) pass through unchanged
    ## as a single column, exactly as before.
    module_trait <- reactive({
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
          lv <- wgcna_trait_levels(input, net$meta, cl)
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

      ## Drop the grey pseudo-module's eigengene - "grey" is WGCNA's bin for
      ## genes that didn't cluster into any real module (see Step 3's own
      ## "grey = unassigned" note and the grey exclusion in Steps 5/6), so
      ## correlating it against traits and surfacing it here would mix an
      ## artifact of leftover genes in with actual co-expression modules.
      MEs_use <- net$MEs[, colnames(net$MEs) != "MEgrey", drop = FALSE]
      ## Uploaded-dataset only: order by hierarchical clustering of the
      ## eigengenes themselves (WGCNA::orderMEs), so correlated modules sit
      ## next to each other the way a published module-trait figure reads,
      ## instead of merge()'s default alphabetical-by-colour sort below.
      ## Scoped to dataset_is_uploaded() so the preloaded/example dataset's
      ## existing module-trait heatmap ordering is untouched either way -
      ## this only changes what an uploaded dataset's heatmap looks like.
      if (dataset_is_uploaded()) MEs_use <- WGCNA::orderMEs(MEs_use)
      cor_mat <- wgcna_cor_fnc(net$cor_method)(MEs_use, traits, use = "p")
      ## Per-column sample size (not one shared nrow(texpr)): a filtered
      ## categorical trait excludes some samples via NA above, so its
      ## p-value must be computed on however many samples that column
      ## actually has, not the full cohort - corPvalueStudent(cor, nSamples)
      ## accepts a matrix of per-cell sample sizes for exactly this case.
      n_per_trait <- vapply(traits, function(x) sum(!is.na(x)), numeric(1))
      n_mat <- matrix(rep(n_per_trait, each = nrow(cor_mat)), nrow = nrow(cor_mat), dimnames = dimnames(cor_mat))
      p_mat <- WGCNA::corPvalueStudent(cor_mat, n_mat)
      rownames(cor_mat) <- sub("^ME", "", rownames(cor_mat))
      rownames(p_mat) <- sub("^ME", "", rownames(p_mat))
      list(cor = cor_mat, p = p_mat)
    })

    observeEvent(module_trait(), {
      mt <- module_trait()
      ## Disease-module cutoff per this project's own methodology - |r| >=
      ## 0.5 AND p < 1e-8 (METHODS_2.4_WGCNA_expanded.md, Section 2.4.9),
      ## not a generic p < 0.05 significance screen: the effect-size
      ## condition is what keeps a module that's merely large-sample-
      ## significant but weakly correlated from counting as disease-
      ## associated.
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
      ## Uploaded-dataset only, matching the equally-scoped orderMEs() call
      ## in module_trait() above: merge(sort = TRUE by default) alphabetises
      ## "module", discarding the eigengene-clustered row order that call
      ## produced - reapplying it explicitly here (reversed, since ggplot's
      ## discrete y axis draws its first factor level at the BOTTOM) is what
      ## makes the heatmap read top-to-bottom in clustered order instead of
      ## A-Z. Left as merge()'s own alphabetical order for the preloaded/
      ## example dataset, exactly as before.
      uploaded <- dataset_is_uploaded()
      if (uploaded) df$module <- factor(df$module, levels = rev(rownames(mt$cor)))
      ## Uploaded-dataset only: two-line "r / (P)" cell label, matching the
      ## published figure's own module-trait panel ("Correlations and P
      ## values are shown in each cell as the first and second lines,
      ## respectively"). formatC(..., format = "e", digits = 0) gives a
      ## single-significant-figure mantissa (e.g. "8e-07"), the same style
      ## as that figure. Left as the plain r-only label for the preloaded/
      ## example dataset, exactly as before.
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
          uiOutput(ns("trait_level_filters_ui"))
        ),
        wgcna_result("table-cells", "Module-trait correlation",
          desc = "Strongly red or blue cells mark modules whose overall expression tracks that trait.",
          withSpinner(plotOutput(ns("mt_plot"), height = 360), color = "#2c6fbb", type = 6),
          div(class = "table-toolbar",
              downloadButton(ns("download_mt_png"), "Download heatmap (PNG)", class = "btn-sm"),
              downloadButton(ns("download_mt_csv"), "Download table (CSV)", class = "btn-sm"))
        )
      )
    })

    ## =====================================================================
    ## Step 5: Hub Genes
    ## =====================================================================

    ## Intramodular connectivity (kWithin: how connected a gene is to other
    ## genes in its OWN module) - depends only on net_result(), so it's
    ## computed once per module-detection run and reused across every module/
    ## trait combination picked below rather than recomputed each time.
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

    ## Same "which levels?" idea as Step 4, for the single trait picked here
    ## - hidden entirely for a numeric trait like age.
    output$hub_trait_levels_ui <- renderUI({
      req(input$hub_trait)
      meta <- net_result()$meta
      if (is.numeric(meta[[input$hub_trait]])) return(NULL)
      lv <- sort(unique(stats::na.omit(as.character(meta[[input$hub_trait]]))))
      checkboxGroupInput(ns("hub_trait_levels"), sprintf("\"%s\" levels to include", input$hub_trait),
                          choices = lv, selected = lv, inline = TRUE)
    })

    hub_data <- reactive({
      net <- net_result()
      req(input$hub_module, input$hub_trait)
      me_col <- paste0("ME", input$hub_module)
      validate(need(me_col %in% colnames(net$MEs), "Selected module not found - re-run Step 3."))

      ## Use this project's own precomputed kME/GS_RA/connectivity (exact
      ## values from the real pipeline, not a live re-derivation) whenever
      ## the current module+trait selection is EXACTLY what that table was
      ## built for: a disease module, trait = "group", and every level
      ## (HC + RA) still included, i.e. no narrower level filter applied
      ## above. Any other selection (a different module, a different
      ## trait, or a level filter) falls through to live computation below
      ## instead, since the precomputed table wouldn't match that choice.
      pre <- net$hub_table_precomputed
      lv <- if (!is.numeric(net$meta[[input$hub_trait]])) input$hub_trait_levels else NULL
      use_precomputed <- !is.null(pre) && identical(input$hub_trait, "group") &&
        input$hub_module %in% unique(pre$module) &&
        (is.null(lv) || setequal(lv, sort(unique(stats::na.omit(as.character(net$meta$group))))))

      if (use_precomputed) {
        sub <- pre[pre$module == input$hub_module, , drop = FALSE]
        ## is_hub carried through as-is (not re-derived from the stored
        ## kME/GS_RA columns, which are rounded for display and no longer
        ## reproduce the exact >0.8/>0.2 boundary the real pipeline used
        ## at full precision) - hub_filtered() below uses it directly at
        ## the documented default threshold so the hub count matches
        ## exactly (217 yellow / 144 brown), and falls back to a live
        ## threshold check on these same kME/GS_RA values for any other
        ## threshold the user sets.
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
    })

    ## Root cause of the network figure's earlier draw failures: a gene with
    ## NA kME or GS (e.g. near-zero variance among the samples actually used
    ## - more likely now that Step 4/5's trait-level filter can narrow the
    ## sample set) makes the `>=` comparison below evaluate to NA, not
    ## FALSE. Subsetting a data.frame on a logical vector containing NA
    ## does NOT drop that row - it inserts a row of NAs (confirmed
    ## directly: df[c(TRUE, NA, FALSE), ] keeps 2 rows, one of them all-NA).
    ## That phantom row's gene name is NA, which flows into
    ## network_genes() -> cyto_export() -> the graph's node names, breaking
    ## ggraph's label layer. `!is.na(...) &` makes NA rows drop out cleanly
    ## instead, the same way FALSE would.
    hub_filtered <- reactive({
      df <- hub_data()
      ## At exactly the documented default threshold (|kME|>0.8, |GS|>0.2)
      ## on precomputed data, use the pipeline's own is_hub flag rather
      ## than re-deriving it from the rounded-for-display kME/GS_RA
      ## columns, which don't reproduce the exact boundary at full
      ## precision (see hub_data() above) - this is what makes the count
      ## match exactly (217 yellow / 144 brown) instead of landing a
      ## handful of genes off. Any other threshold recomputes from the
      ## same kME/GS values as always.
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

    ## Module membership vs intramodular connectivity - both measure how
    ## "central" a gene is to its module, from different angles; genes high
    ## on both are the most defensible hub-gene calls.
    output$connectivity_plot <- renderPlot({
      df <- hub_data()
      ggplot(df, aes(x = kME, y = kWithin)) +
        geom_point(color = "#2c6fbb", alpha = 0.6, size = 1.6) +
        labs(x = sprintf("Module membership (kME) in %s", input$hub_module),
             y = "Intramodular connectivity (kWithin)") +
        theme_minimal(base_size = 12)
    })

    ## Module eigengene vs the chosen trait, one point/box per sample -
    ## the same relationship Step 4's heatmap summarises as a single r,
    ## shown here at the sample level.
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

    ## The network figure/export always has something to show: it uses the
    ## kME/GS-thresholded hub genes when there are at least 3 of them, and
    ## otherwise falls back to the top-connected genes (by |kME|) in the
    ## selected module, capped at 30 for a readable figure - so tightening
    ## the sliders above to "nothing passes" doesn't leave the network blank.
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

    ## An actual <img> fetched live from STRING's own REST image API
    ## (string-db.org/api/image/network) - a real database-hosted figure for
    ## these genes, distinct from the in-app network_plot above which is
    ## built from this run's own WGCNA edges. Capped at 60 genes: STRING's
    ## own image rendering gets unreadable well before that, and the browser
    ## simply hides the <img> (via onerror below) if the request ever fails.
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

    ## On-demand Cytoscape export for network_genes() above (WGCNA::
    ## exportNetworkToCytoscape with edgeFile/nodeFile = NULL returns the
    ## edge/node tables instead of writing to disk, so the CSVs below are
    ## written directly inside the downloadHandler).
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

    ## An actual network figure, in the app, built from the exact same
    ## edge/node data the Cytoscape CSV export above uses - so you can see
    ## the network without opening the separate Cytoscape desktop app.
    ## Nodes are colored by their own WGCNA module color (a valid R color
    ## name already), edge thickness is TOM weight; layout is force-directed
    ## (Fruchterman-Reingold, igraph's default "fr" layout).
    ## Every step below is wrapped so a failure turns into a plain validate()
    ## message in the plot panel - never a raw, unreadable error object -
    ## and the real cause still lands in the server log via message() for
    ## debugging.
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

      ## geom_edge_link0(), not geom_edge_link(): found by reproducing this
      ## exact failure outside Shiny with this app's full ~30-package stack
      ## loaded (it never failed with a smaller package set, which is why
      ## earlier fixes here didn't stick). The real error, unwrapped via
      ## rlang::cnd_message(e, inherit = TRUE): "Caused by error in
      ## pathAttr(): SET_VECTOR_ELT() can only be applied to a 'list', not a
      ## 'integer'" - geom_edge_link()'s internals hit a generic-function
      ## namespace collision with one of the Bioconductor packages this app
      ## loads (IRanges/S4Vectors/BiocGenerics, pulled in by clusterProfiler/
      ## org.Hs.eg.db/coloc) - the same class of collision global.R's own
      ## `cor <- WGCNA::cor` / `validate <- shiny::validate` lines already
      ## document and work around for other functions. geom_edge_link0()
      ## draws the same straight, width-mapped edges through simpler
      ## grid internals (segmentsGrob) that don't hit the broken code path -
      ## confirmed by reproducing this exact module/gene set with both,
      ## outside Shiny, using ragg::agg_png (Shiny's own renderPlot device).
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

    ## ggplot/ggraph objects are lazily evaluated: building one with `+`
    ## (inside network_plot_obj() above) doesn't run the force-directed
    ## layout algorithm or ggrepel's label-placement step - that only
    ## happens when the plot is actually printed/drawn. A failure there
    ## would previously happen OUTSIDE any tryCatch (network_plot_obj() had
    ## already returned successfully), surfacing as a raw, unreadable error.
    ## print() is called explicitly here so a draw-time failure is caught
    ## too, not just a build-time one.
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
          ## Defaults 0.8/0.2, not 0.7/0.2: this project's own hub-gene rule
          ## (Chapter_2_subchapter2_sexstratified.md Section 2.4.12,
          ## METHODS_2.4_WGCNA_expanded.md line 445) is |kME| > 0.8 AND
          ## |GS| > 0.2 - requiring both means a hub must be central to its
          ## module AND disease-associated; kME alone (the old 0.7 default)
          ## would just restate module membership. Both sliders remain
          ## freely adjustable.
          p(class = "submodule-desc", "Hub genes are central to their module (high kME) and associated with your trait (high GS)."),
          uiOutput(ns("hub_picker_ui")),
          uiOutput(ns("hub_trait_levels_ui")),
          fluidRow(
            column(6, sliderInput(ns("kme_thr"), "Minimum |kME|", min = 0, max = 1, value = 0.8, step = 0.05)),
            column(6, sliderInput(ns("gs_thr"), "Minimum |GS|", min = 0, max = 1, value = 0.2, step = 0.05))
          )
        ),
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
    })

    ## =====================================================================
    ## Step 6: Enrichment (GO / KEGG, via clusterProfiler + org.Hs.eg.db -
    ## the same databases mod_enrichment.R already uses for user-pasted gene
    ## lists; here the gene list comes from a detected module instead).
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

    module_enrich <- eventReactive(input$run_enrich_btn, {
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

      ## Background = the genes actually analyzed in Step 1 (the network's
      ## own universe), not the whole loaded dataset - the statistically
      ## correct background for a WGCNA module enrichment test.
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
    }, ignoreInit = TRUE)

    output$enrich_summary_ui <- renderUI({
      res <- tryCatch(module_enrich(), error = function(e) NULL)
      if (is.null(res)) {
        return(div(class = "empty-note", icon("circle-info"), "Not run yet. Click \"Run enrichment\" above."))
      }
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

    output$step6_ui <- renderUI({
      tagList(
        box(
          width = 12, title = tagList(icon("flask"), " Module enrichment"), status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Tests whether a module (or its hub genes) is enriched for GO terms or KEGG pathways, against the genes analyzed in Step 1 as background."),
          uiOutput(ns("enrich_picker_ui"))
        ),
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
    ## Pipeline summary - same vertical-timeline card component
    ## (pipeline_summary_ui(), R/ui_shell.R) Preprocessing uses, not a
    ## separate custom status display.
    ## =====================================================================

    wgcna_progress <- reactive({
      power_ok <- !is.null(tryCatch(sft_result(), error = function(e) NULL))
      modules_ok <- !is.null(tryCatch(net_result(), error = function(e) NULL))
      traits_ok <- modules_ok && length(input$trait_cols %||% character(0)) > 0
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

    ## Every step's uiOutput must render as soon as this module mounts, not
    ## only once its own tab is actually clicked into view - Shiny suspends
    ## renderUI/plotOutput bindings in hidden tabs by default
    ## (suspendWhenHidden = TRUE), so without this, a step's own inputs
    ## (e.g. Step 2's manual power numericInput) simply don't exist
    ## server-side until that tab has been visited at least once. That
    ## silently broke Step 3's "Run": effective_power() reads
    ## input$manual_power, which was NULL (never rendered) for anyone who
    ## went straight from Step 1 to Step 3 - req() inside it failed, caught
    ## as NA by the run_btn observer's tryCatch, then blocked by a
    ## validate() that (being inside a plain observeEvent, not a render
    ## context) showed NO message anywhere: "Run" just silently did
    ## nothing, forever, with zero visible feedback of why. Must come after
    ## every output$stepN_ui <- renderUI(...) above - outputOptions() errors
    ## if the named output doesn't exist in `output` yet.
    for (step in paste0("step", 1:6, "_ui")) outputOptions(output, step, suspendWhenHidden = FALSE)
  })
}
