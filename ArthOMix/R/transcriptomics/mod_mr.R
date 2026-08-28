## R/mod_mr.R
## Submodule: Mendelian Randomization
## A general-purpose two-sample MR tool. "Your analysis" runs a live MR test
## for any gene against a bundled RA dataset, or against your own uploaded
## exposure/outcome GWAS summary statistics for any trait - unlike the rest
## of this app (which reads a currently loaded expression matrix), GWAS
## summary statistics have no single "your dataset" equivalent, so both
## paths live side by side behind one switch.
##
## Bundled dataset: data/processed/new/MR_primary_objects.rds - cis-
## restricted (+/-1Mb of the gene body, GRCh37), MHC-flagged, and harmonised
## (action=2) cis-eQTL instruments against the Okada 2014 RA GWAS (ieu-a-832)
## for 1,701 genes with at least one usable (mr_keep) harmonised instrument.
## The advanced filters below (instrument p-value, F-statistic, MHC
## handling) subset THIS cached, already-cis-filtered object live, rather
## than re-deriving cis restriction from scratch (that needs a GRCh37
## annotation this app does not otherwise load - 20% of unfiltered
## instruments would sit on the wrong chromosome entirely without it).
##
## TWO layers, not one - "Your analysis" (the interactive tool, ArthOChat
## shortcut top-left, matching every other sub-module's layout) comes first
## on the page; the bundled RA reference results follow as their own
## full-width box, gated behind a click like every other analysis in this
## app (nothing computes or renders until asked):
##  1. "Your analysis" runs the estimator hierarchy LIVE, gene by gene or in
##     a batch screen, starting from recommended default filters so a user
##     who changes nothing reproduces layer 2's numbers for that gene.
##  2. The bundled RA reference results read pre-computed saved CSVs
##     (MR_MHC_sensitivity_{sex}.csv etc.) directly - the exact, already-
##     FDR-corrected, already-MHC-sensitivity-tested numbers for this
##     dataset. Nothing here is recomputed.
## A live data-integrity check (`relabel_check` below) verified that
## `$dat$gene` in the cached RDS is stale for 21 SNPs / 37 genes that share
## an eQTL with a neighbouring gene - the pair-matched relabelling fix was
## written into the upstream pipeline but the cached object predates it.
## This module recomputes the correct label from `$inst` (SNP + exposure-
## dataset ID pair) at load time rather than trusting the cached column, so
## layer 1 is never affected by that staleness; layer 2 is unaffected too,
## because none of the 37 genes are in either sex's FDR-surviving list.

## Small hover tooltip (native browser title attribute - no extra JS/CSS
## dependency) for a filter's rationale, so the control row stays one line
## instead of a paragraph underneath it.
mr_info_tip <- function(text) tags$span(icon("circle-info", style = "color:#8A929C; cursor: help; margin-left: 4px;"), title = text)

mod_mr_config <- list(
  id = "mr", group = "Genetics",
  title = "Mendelian Randomization",
  description = "Two-sample Mendelian randomisation on preloaded or uploaded data. GWAS summary statistics for any trait are needed.",
  icon = "route"
)

mod_mr_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tabsetPanel(
      id = ns("mr_top_tabs"), type = "tabs",
      tabPanel("MR overall", br(), mod_mr_overall_ui(ns)),
      tabPanel("MR female", br(), mod_mr_sex_tab_ui(ns, "female")),
      tabPanel("MR male", br(), mod_mr_sex_tab_ui(ns, "male"))
    )
  )
}

## "MR female"/"MR male": the bundled dataset's own reference results for
## that sex - but, like every other analysis in this app (Run MR, batch
## screen, ...), nothing renders - not the stats, not the plot, not the
## table - until the user explicitly clicks "Show". Opening the tab alone
## isn't a click; the button inside it is.
mod_mr_sex_tab_ui <- function(ns, sx) {
  uiOutput(ns(paste0(sx, "_tab_body")))
}

mod_mr_overall_ui <- function(ns) {
  tagList(
      fluidRow(
        column(
          4,
          box(
            width = NULL, title = "Exposure & outcome data", status = "primary", solidHeader = FALSE,
            radioButtons(
              ns("data_source"), NULL,
              choiceNames = list(
                tagList(icon("database"), " Bundled RA dataset (default)"),
                tagList(icon("upload"), " Upload your own GWAS summary statistics")
              ),
              choiceValues = list("project", "upload"), selected = "project"
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'project'", ns("data_source")),
              uiOutput(ns("gene_ui")),
              div(class = "empty-note", style = "margin-top: 8px; font-size: 12.5px;",
                  icon("circle-info"), strong(" GWAS used: "),
                  "outcome = ", tags$a(href = "https://gwas.mrcieu.ac.uk/datasets/ieu-a-832/", target = "_blank", "Okada et al. 2014 RA GWAS (OpenGWAS ieu-a-832)"),
                  "; exposure = ", tags$a(href = "https://pubmed.ncbi.nlm.nih.gov/34475573", target = "_blank", "eQTLGen whole-blood cis-eQTL (Vösa et al. 2021)"), ".")
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'upload'", ns("data_source")),
              p(class = "submodule-desc",
                "Two delimited files (CSV/TSV), one row per SNP: your own eQTL/exposure summary statistics and an outcome GWAS's. Not sure which outcome dataset fits your trait? ",
                strong("Ask ArthOChat"), " to search the OpenGWAS catalogue for it."),
              textInput(ns("upload_label"), "Exposure name (for labelling only)", value = "Uploaded exposure", width = "100%"),
              fileInput(ns("exp_file"), "Exposure file", accept = c(".csv", ".tsv", ".txt")),
              uiOutput(ns("exp_map_ui")),
              fileInput(ns("out_file"), "Outcome file", accept = c(".csv", ".tsv", ".txt")),
              uiOutput(ns("out_map_ui")),
              uiOutput(ns("upload_data_summary")),
              div(
                style = "display:flex; align-items:center; gap:4px; margin-top: 8px;",
                checkboxInput(ns("do_clump"), "LD-clump exposure SNPs before harmonising (needs network access)", value = FALSE),
                mr_info_tip("Runs ieugwasr::ld_clump() live against the OpenGWAS LD reference to drop correlated SNPs, the same independence step the bundled dataset's own instruments were built with. Off by default since it needs network access; if it fails, filtering proceeds on the unclumped set.")
              ),
              conditionalPanel(
                condition = sprintf("input['%s']", ns("do_clump")),
                fluidRow(
                  column(6, selectInput(ns("clump_r2"), "Clump r² threshold", c("0.1" = "0.1", "0.01" = "0.01", "0.001 (default)" = "0.001"), selected = "0.001", width = "100%")),
                  column(6, selectInput(ns("clump_kb"), "Clump window (kb)", c("250" = "250", "1,000" = "1000", "10,000 (default)" = "10000"), selected = "10000", width = "100%"))
                ),
                selectInput(ns("clump_pop"), "LD reference population", c("European (EUR)" = "EUR", "East Asian (EAS)" = "EAS", "South Asian (SAS)" = "SAS", "African (AFR)" = "AFR", "Admixed American (AMR)" = "AMR"), selected = "EUR", width = "100%")
              )
            )
          ),
          box(
            width = NULL, title = "Instrument & analysis filters", status = "primary", solidHeader = FALSE,
            div(class = "empty-note", icon("circle-check"),
                "Pre-selected to recommended defaults - leave as-is and click “Run MR” to match the bundled reference estimate for whichever gene you pick (Result panel shows both side by side)."),
            div(style = "display:flex; align-items:center; gap:4px;",
              selectInput(
                ns("pval_cut"), "eQTL p-value threshold",
                choices = c(
                  "Default: p < 5×10⁻⁸" = "5e-8", "Stricter: p < 5×10⁻⁹" = "5e-9",
                  "Stricter: p < 1×10⁻¹⁰" = "1e-10", "Stricter: p < 1×10⁻¹²" = "1e-12"
                ),
                selected = "5e-8", selectize = FALSE, width = "100%"
              ),
              mr_info_tip("Genome-wide significance. Every cached instrument already clears this (min F = 29.7).")
            ),
            div(style = "display:flex; align-items:center; gap:4px;",
              selectInput(
                ns("fstat_cut"), "Minimum F-statistic",
                choices = c(
                  "Default: F ≥ 10" = "10", "No additional filter: F ≥ 0" = "0",
                  "Stricter: F ≥ 30 (~median)" = "30", "Stricter: F ≥ 100" = "100"
                ),
                selected = "10", selectize = FALSE, width = "100%"
              ),
              mr_info_tip("Conventional weak-instrument threshold. Non-binding at the default p-value cutoff - p<5×10⁻⁸ already implies F≥29.7.")
            ),
            div(style = "display:flex; align-items:center; gap:4px;",
              selectInput(
                ns("max_snps"), "Cap instruments used per gene (strongest by p-value)",
                choices = c("No limit (default)" = "Inf", "Top 20" = "20", "Top 10" = "10", "Top 5" = "5", "Top 3" = "3"),
                selected = "Inf", selectize = FALSE, width = "100%"
              ),
              mr_info_tip("Explore how the estimate degrades with fewer, stronger instruments - e.g. to stress-test sensitivity to instrument count.")
            ),
            div(style = "display:flex; align-items:center; gap:4px;",
              selectInput(
                ns("ci_level"), "Confidence interval level",
                choices = c("90%" = "0.90", "95% (default)" = "0.95", "99%" = "0.99"),
                selected = "0.95", selectize = FALSE, width = "100%"
              ),
              mr_info_tip("Widens or narrows every reported confidence interval; does not change the point estimate or p-value.")
            ),
            div(style = "display:flex; align-items:center; gap:4px;",
              numericInput(ns("alpha_cut"), "Significance threshold (α) for flagging results", value = 0.05, min = 0.0001, max = 0.5, step = 0.005, width = "100%"),
              mr_info_tip("Used to flag a significant result in the Result panel and to count/highlight significant genes in the batch screen below.")
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'project'", ns("data_source")),
              div(style = "display:flex; align-items:center; gap:4px;",
                radioButtons(
                  ns("mhc_mode"), "Extended MHC region (chr6:25-34Mb)",
                  choiceNames = list("Include (default)", "Exclude (sensitivity)"),
                  choiceValues = list("include", "exclude"), selected = "include", inline = TRUE
                ),
                mr_info_tip("HLA-DRB1 is the dominant RA locus with the genome's longest-range LD, so any MHC gene's cis-eQTL is confounded with it. MHC-flagged genes are best read as “associated, not causal” unless they also survive with MHC instruments excluded.")
              )
            ),
            div(style = "display:flex; align-items:center; gap:4px;",
              checkboxInput(ns("include_mode"), "Also compute weighted-mode (exploratory, ≥3 instruments)", value = FALSE)
            ),
            div(style = "display:flex; align-items:center; gap:4px;",
              checkboxInput(ns("run_presso"), "Also run MR-PRESSO outlier test (≥4 instruments)", value = FALSE),
              mr_info_tip("Simulation-based test (MRPRESSO::mr_presso) for horizontal-pleiotropy outlier instruments, with an outlier-corrected estimate when any are found. Off by default since it's simulation-based and adds a few seconds per run.")
            ),
            div(style = "display: flex; gap: 8px; margin-top: 6px;",
                actionButton(ns("run_btn"), "Run MR", icon = icon("play"), class = "btn-primary btn-sm"),
                actionButton(ns("reset_btn"), "Reset to defaults", icon = icon("rotate-left"), class = "btn-default btn-sm"))
          )
        ),
        column(
          8,
          ## No box/title here until Run MR is clicked (mr_has_run) - the
          ## trigger lives in the "Instrument filters" box in the left
          ## column, so hiding these entirely first-load costs nothing for
          ## discoverability. Result comes last, below the plots it explains.
          ## Pre-click, uiOutput(ns("scatter_diag_ui")) itself fills this
          ## space with a plain empty-state card (see server) instead of
          ## leaving the whole column blank next to a much taller filters
          ## panel - still no result content, box title or number shown.
          uiOutput(ns("scatter_diag_ui")),
          uiOutput(ns("result_box_ui"))
        )
      ),
      ## Nothing below - title included - exists until its own trigger is
      ## clicked. Each uiOutput returns NULL (no box, no chrome) pre-click.
      uiOutput(ns("instrument_mr_tables_ui")),
      uiOutput(ns("batch_screen_section"))
  )
}

mod_mr_server <- function(id, dataset, results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## mr_obj$dat: one row per (gene, SNP), already cis-restricted and
    ## MHC-flagged, harmonised against Okada 2014 - see the header comment
    ## above for why this file (not the raw MR_instruments.rds this module
    ## used previously) is the right base for live filtering. Loading,
    ## mr_keep filtering, and the SNP->gene relabel fix (see global.R's
    ## load_mr_instrument_table() for the full rationale) live in global.R
    ## now, shared with mod_crossancestry.R's live-upload arm, which
    ## re-harmonises these SAME exposure-side rows against a different
    ## outcome GWAS.
    mr_loaded <- load_mr_instrument_table()
    mr_dat_all <- mr_loaded$dat
    mr_primary_all <- mr_loaded$primary
    available_genes <- sort(unique(mr_dat_all$gene))

    ## Transparency record of exactly what the relabelling changed, surfaced
    ## to the user via data_quality_note below rather than fixed silently -
    ## a correctness fix a user can't see isn't verifiable.
    relabel_check <- mr_loaded$relabel_check

    ## -----------------------------------------------------------------------
    ## "This project's own MR results" (Section 2.6): read straight from the
    ## pipeline's saved output, never recomputed. MR_MHC_sensitivity_{sex}.csv
    ## covers every gene tested in that stratum (1,477 female / 1,478 male)
    ## with the primary estimate, its 95% CI, and the MHC-exclusion verdict
    ## already computed - filtering to FDR_primary < 0.05 reproduces
    ## FS_input_{sex}.csv / MR_causal_FDR_{sex}.csv exactly (verified: 32
    ## female / 25 male rows) while also carrying OR_lo/OR_hi for the forest
    ## plot and the ROBUST/UNTESTABLE/MHC-DEPENDENT/FDR-RANK-ONLY verdict
    ## from the full MHC sensitivity re-run (Section 2.6.8).
    ## -----------------------------------------------------------------------
    verdict_short <- function(v) {
      dplyr::case_when(
        startsWith(v, "ROBUST")      ~ "ROBUST",
        startsWith(v, "UNTESTABLE")  ~ "UNTESTABLE w/o MHC",
        startsWith(v, "MHC-DEPENDENT") ~ "MHC-DEPENDENT",
        startsWith(v, "FDR-RANK")    ~ "FDR-rank only",
        TRUE ~ "ns in both"
      )
    }
    verdict_color <- c(
      "ROBUST" = ARTHOMIX_STATUS$good, "UNTESTABLE w/o MHC" = ARTHOMIX_STATUS$warning,
      "MHC-DEPENDENT" = ARTHOMIX_STATUS$critical, "FDR-rank only" = "#8A929C", "ns in both" = "#8A929C"
    )
    verdict_badge <- function(v) sprintf(
      '<span style="background:%s;color:#fff;font-size:11px;font-weight:600;padding:2px 8px;border-radius:10px;white-space:nowrap;">%s</span>',
      verdict_color[v], v
    )

    project_survivors <- lapply(list(female = "female", male = "male"), function(sx) {
      d <- read_table_safe(sprintf("MR_MHC_sensitivity_%s.csv", sx))
      req_cols <- c("gene", "MHC_gene", "nSNP_primary", "method_primary", "OR_primary",
                    "OR_lo_primary", "OR_hi_primary", "p_primary", "FDR_primary", "verdict")
      if (is.null(d) || !all(req_cols %in% colnames(d))) return(NULL)
      d <- d[d$FDR_primary < 0.05, req_cols, drop = FALSE]
      d$direction <- ifelse(d$OR_primary > 1, "risk (OR>1)", "protective (OR<1)")
      d$verdict_short <- verdict_short(d$verdict)
      d[order(d$p_primary), , drop = FALSE]
    })

    project_n_tested <- vapply(c(female = "female", male = "male"), function(sx) {
      d <- read_table_safe(sprintf("MR_MHC_sensitivity_%s.csv", sx))
      if (is.null(d)) NA_integer_ else nrow(d)
    }, integer(1))
    project_n_candidates <- vapply(c(female = "female", male = "male"), function(sx) {
      d <- read_table_safe(sprintf("candidates_%s_disease.csv", sx))
      if (is.null(d)) NA_integer_ else length(unique(d$gene))
    }, integer(1))

    ## -----------------------------------------------------------------------
    ## Upload-your-own-GWAS mode: parse a user's exposure/outcome summary
    ## statistics files with TwoSampleMR::format_data() + harmonise_data() -
    ## the SAME functions scripts/00_shared/10_MR.R itself calls (`::`-
    ## qualified so this file never attaches the package and risks the
    ## dplyr::slice-style namespace collisions that motivated 10_MR.R's own
    ## D1 fix) - rather than hand-rolling allele-flipping/strand logic. Both
    ## functions are pure data-frame operations, no OpenGWAS network access
    ## needed, so this works entirely offline on whatever the user uploads.
    ##
    ## Column-guessing/mapping (GWAS_COL_PATTERNS, guess_gwas_col(),
    ## read_uploaded_table(), gwas_col_map_ui()) lives in global.R now,
    ## shared with mod_coloc.R and mod_crossancestry.R's own upload modes.
    ## -----------------------------------------------------------------------

    ## Each uploaded file is parsed from disk exactly ONCE, cached by Shiny's
    ## own reactive graph (re-reads only when input$exp_file/out_file itself
    ## changes) - three different consumers (the column-mapping UI below,
    ## the "Your data" summary, and mr_result_uploaded() when Run MR is
    ## clicked) previously each called read_uploaded_table() independently,
    ## paying the parse cost up to 3x per file for no reason, which is what
    ## made switching between the mapping UI and running MR feel slow on a
    ## large outcome GWAS file.
    exp_df_r <- reactive({ req(input$exp_file); read_uploaded_table(input$exp_file$datapath) })
    out_df_r <- reactive({ req(input$out_file); read_uploaded_table(input$out_file$datapath) })

    output$exp_map_ui <- gwas_col_map_ui(ns, reactive(input$exp_file), exp_df_r, "exp", "Exposure file")
    output$out_map_ui <- gwas_col_map_ui(ns, reactive(input$out_file), out_df_r, "out", "Outcome file")

    ## "Your data" reference summary - the upload-mode counterpart to the
    ## "GWAS used" note shown for the bundled dataset. Updates the instant a
    ## file is selected (no "Run MR" click needed), reusing exp_df_r/
    ## out_df_r above rather than re-reading either file.
    output$upload_data_summary <- renderUI({
      req(identical(input$data_source, "upload"))
      if (is.null(input$exp_file) && is.null(input$out_file)) return(NULL)
      exp_part <- if (!is.null(input$exp_file)) {
        df <- tryCatch(exp_df_r(), error = function(e) NULL)
        label <- if (nzchar(trimws(input$upload_label %||% ""))) trimws(input$upload_label) else "Uploaded exposure"
        sprintf("%s - %s (%s SNPs)", label, input$exp_file$name, if (!is.null(df)) format(nrow(df), big.mark = ",") else "reading...")
      } else "not yet uploaded"
      out_part <- if (!is.null(input$out_file)) {
        df <- tryCatch(out_df_r(), error = function(e) NULL)
        sprintf("%s (%s SNPs)", input$out_file$name, if (!is.null(df)) format(nrow(df), big.mark = ",") else "reading...")
      } else "not yet uploaded"
      div(class = "empty-note", style = "margin-top: 8px; font-size: 12.5px;",
          icon("circle-info"), strong(" Your data: "), "exposure = ", exp_part, "; outcome = ", out_part, ".")
    })

    output$gene_ui <- renderUI({
      ## selectize = TRUE here (unlike the plain <select> this app uses for
      ## every other single-pick dropdown, e.g. mod_dge.R's contrast_col) is
      ## deliberate: this list has ~1,700 genes, not 2-20, and needs to be
      ## searchable rather than scrolled.
      default_gene <- if ("LPCAT2" %in% available_genes) "LPCAT2" else available_genes[1]
      tagList(
        selectInput(ns("gene"), "Candidate gene", choices = available_genes, selected = default_gene, selectize = TRUE),
        p(class = "submodule-desc", style = "margin-bottom: 0; font-size: 12px;",
          sprintf("%s genes have a cached, harmonised cis-eQTL instrument available.", format(length(available_genes), big.mark = ",")))
      )
    })

    ## -----------------------------------------------------------------------
    ## Project-results outputs (stats, forest plot, table, downloads) per
    ## sex - written once as factories over `sx` so the female/male panels
    ## can never drift apart in how they're built.
    ## -----------------------------------------------------------------------
    ## One factory shared by the "MR female"/"MR male" top-level tabs so
    ## they can never drift apart. Each tab gets its OWN output id
    ## (data_quality_note_female / _male) - two uiOutputs bound to the same
    ## Shiny output id would collide in the DOM, since both tabs exist in
    ## the page simultaneously (tabsetPanel keeps inactive tab content in
    ## the DOM, just hidden) rather than one being inserted only on demand.
    render_data_quality_note <- function() renderUI({
      rc <- relabel_check
      if (rc$n_snp == 0) return(NULL)
      surv_genes <- unique(unlist(lapply(project_survivors, function(d) if (!is.null(d)) d$gene)))
      overlap <- intersect(rc$genes, surv_genes)
      tags$span(
        mr_info_tip(sprintf(
          "Data-integrity check: %d shared instrument SNPs across %d genes (e.g. %s) were mislabelled to the wrong gene in the cached file - auto-corrected here before every estimate. %s",
          rc$n_snp, length(rc$genes), paste(utils::head(rc$genes, 3), collapse = ", "),
          if (length(overlap) == 0) "None affect the survivor lists below."
          else sprintf("Affected genes in a survivor list: %s.", paste(overlap, collapse = ", "))
        )),
        " ", em(style = "color:#8A929C; font-size:12px;", sprintf("Auto-corrected %d relabelled SNP%s (hover for detail)", rc$n_snp, if (rc$n_snp != 1) "s" else ""))
      )
    })
    output$data_quality_note_female <- render_data_quality_note()
    output$data_quality_note_male   <- render_data_quality_note()

    render_proj_stats <- function(sx) renderUI({
      surv <- project_survivors[[sx]]
      n_robust <- if (!is.null(surv)) sum(surv$verdict_short == "ROBUST") else 0L
      fluidRow(
        valueBox(format(project_n_candidates[[sx]], big.mark = ","), "Candidate genes (DE screen)", icon = icon("dna"), color = "light-blue", width = 3),
        valueBox(format(project_n_tested[[sx]], big.mark = ","), "MR-tested (cis instrument found)", icon = icon("route"), color = "purple", width = 3),
        valueBox(if (!is.null(surv)) nrow(surv) else 0, "Survive FDR < 0.05", icon = icon("check"), color = "yellow", width = 3),
        valueBox(n_robust, "ROBUST after MHC exclusion", icon = icon("shield-halved"), color = "green", width = 3)
      )
    })
    output$proj_stats_female <- render_proj_stats("female")
    output$proj_stats_male   <- render_proj_stats("male")

    ## Shared by the on-screen plot and its PNG download, so the two can
    ## never show different figures. Title/subtitle kept short and wrapped
    ## across two short lines (rather than one long technical sentence) so
    ## nothing overflows the plot panel and gets clipped at narrow widths;
    ## the method-hierarchy detail lives in the "Result" panel and page text
    ## elsewhere instead of being crammed into the figure itself.
    build_proj_forest <- function(sx) {
      surv <- project_survivors[[sx]]
      validate(need(!is.null(surv) && nrow(surv) > 0, "No reference results available for this sex."))
      d <- surv[order(surv$OR_primary), , drop = FALSE]
      d$gene <- factor(d$gene, levels = d$gene)
      ggplot(d, aes(OR_primary, gene, colour = verdict_short)) +
        geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
        geom_errorbarh(aes(xmin = OR_lo_primary, xmax = OR_hi_primary), height = 0, linewidth = 0.45) +
        geom_point(size = 2) +
        scale_x_log10() +
        scale_colour_manual(values = verdict_color, name = "MHC-sensitivity verdict",
                             breaks = c("ROBUST", "UNTESTABLE w/o MHC", "MHC-DEPENDENT", "FDR-rank only")) +
        labs(title = sprintf("%s: %d genes, FDR < 0.05 vs Okada 2014 RA GWAS", tools::toTitleCase(sx), nrow(d)),
             subtitle = "Colour = verdict under MHC-instrument-exclusion sensitivity re-run",
             x = "Odds ratio for RA (log scale)", y = NULL) +
        theme_arthomix(base_size = 11) +
        theme(axis.text.y = element_text(size = 8), legend.position = "top",
              plot.title = element_text(size = 12.5), plot.subtitle = element_text(size = 10))
    }
    render_proj_forest <- function(sx) renderPlot(build_proj_forest(sx))
    output$proj_forest_female <- render_proj_forest("female")
    output$proj_forest_male   <- render_proj_forest("male")

    make_dl_forest_png <- function(sx) downloadHandler(
      filename = function() sprintf("MR_%s_forest_plot.png", sx),
      content = function(file) {
        n_genes <- max(1, nrow(project_survivors[[sx]]))
        ggsave(file, plot = build_proj_forest(sx), width = 9, height = max(5, n_genes * 0.16), dpi = 300, bg = "white", limitsize = FALSE)
      }
    )
    output$dl_proj_forest_female <- make_dl_forest_png("female")
    output$dl_proj_forest_male   <- make_dl_forest_png("male")

    render_proj_table <- function(sx) DT::renderDataTable({
      surv <- project_survivors[[sx]]
      validate(need(!is.null(surv), "Reference results not found."))
      out <- data.frame(
        Gene = surv$gene, Direction = surv$direction, OR = round(surv$OR_primary, 3),
        `95% CI` = sprintf("%.3f-%.3f", surv$OR_lo_primary, surv$OR_hi_primary),
        p = signif(surv$p_primary, 3), FDR = signif(surv$FDR_primary, 3),
        nSNP = surv$nSNP_primary, Method = surv$method_primary,
        `In MHC` = ifelse(surv$MHC_gene, "Yes", "No"),
        Verdict = vapply(surv$verdict_short, verdict_badge, character(1)),
        check.names = FALSE
      )
      DT::datatable(out, rownames = FALSE, escape = FALSE, filter = "top",
                     options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    output$proj_table_female <- render_proj_table("female")
    output$proj_table_male   <- render_proj_table("male")

    make_dl_survivors <- function(sx) downloadHandler(
      filename = function() sprintf("MR_%s_RA_reference_survivors.csv", sx),
      content = function(file) write.csv(project_survivors[[sx]], file, row.names = FALSE)
    )
    output$dl_proj_survivors_female <- make_dl_survivors("female")
    output$dl_proj_survivors_male   <- make_dl_survivors("male")

    make_dl_xlsx <- function(sx) downloadHandler(
      filename = function() sprintf("MR_%s_all_tables.xlsx", sx),
      content = function(file) file.copy(file.path(TABLES_DIR, sprintf("MR_%s_all_tables.xlsx", sx)), file, overwrite = TRUE)
    )
    output$dl_proj_xlsx_female <- make_dl_xlsx("female")
    output$dl_proj_xlsx_male   <- make_dl_xlsx("male")

    ## "MR female"/"MR male" tab bodies: nothing (no stats, no plot, no
    ## table) until the user clicks "Show ... MR results" - same click-gated
    ## pattern as "Run MR" and the batch screen, applied here too so the
    ## whole module behaves consistently: opening a tab is not, itself, a
    ## request to compute or reveal anything.
    sex_shown <- reactiveValues(female = FALSE, male = FALSE)
    observeEvent(input$show_female_btn, sex_shown$female <- TRUE)
    observeEvent(input$show_male_btn, sex_shown$male <- TRUE)

    render_sex_tab_body <- function(sx) renderUI({
      if (!isTRUE(sex_shown[[sx]])) {
        return(div(
          style = "margin: 16px 0;",
          p(class = "submodule-desc",
            sprintf("Pre-computed MR results for %s differential-expression candidate genes against the bundled RA GWAS (Okada et al. 2014, OpenGWAS ieu-a-832).", sx)),
          actionButton(ns(paste0("show_", sx, "_btn")), sprintf("Show %s MR results", sx), icon = icon("chart-bar"), class = "btn-primary btn-sm")
        ))
      }
      tagList(
        uiOutput(ns(paste0("data_quality_note_", sx))),
        uiOutput(ns(paste0("proj_stats_", sx))),
        div(class = "table-toolbar", downloadButton(ns(paste0("dl_proj_forest_", sx)), "Forest plot (PNG)", class = "btn-sm")),
        withSpinner(plotOutput(ns(paste0("proj_forest_", sx)), height = 560), color = "#2c6fbb", type = 6),
        div(style = "margin-top: 16px;",
          div(class = "table-toolbar",
              downloadButton(ns(paste0("dl_proj_survivors_", sx)), "Survivors (CSV)", class = "btn-sm"),
              downloadButton(ns(paste0("dl_proj_xlsx_", sx)), "Full TABLE1-4 (XLSX)", class = "btn-sm")),
          DT::dataTableOutput(ns(paste0("proj_table_", sx)))
        )
      )
    })
    output$female_tab_body <- render_sex_tab_body("female")
    output$male_tab_body   <- render_sex_tab_body("male")

    observeEvent(input$reset_btn, {
      updateSelectInput(session, "pval_cut", selected = "5e-8")
      updateSelectInput(session, "fstat_cut", selected = "10")
      updateSelectInput(session, "max_snps", selected = "Inf")
      updateSelectInput(session, "ci_level", selected = "0.95")
      updateNumericInput(session, "alpha_cut", value = 0.05)
      updateRadioButtons(session, "mhc_mode", selected = "include")
      updateCheckboxInput(session, "include_mode", value = FALSE)
      updateCheckboxInput(session, "run_presso", value = FALSE)
    })

    ## estimate_mr_set() - the project's own estimator hierarchy (Section
    ## 2.6.5): >=3 SNPs get IVW + MR-Egger + weighted median (+ weighted mode
    ## if requested), 2 SNPs get IVW alone, 1 SNP gets the Wald ratio - now
    ## lives in global.R, shared with mod_crossancestry.R's live-upload arm
    ## so the two can never silently drift apart on how a gene is estimated.

    ## Shared by both branches below: current confidence level as a
    ## fraction, and an optional cap on how many instruments (strongest by
    ## exposure p-value) actually reach estimate_mr_set() - lets a user
    ## stress-test how an estimate degrades with fewer, stronger instruments.
    current_ci_level <- function() as.numeric(input$ci_level %||% "0.95")
    apply_max_snps_cap <- function(d) {
      cap <- suppressWarnings(as.numeric(input$max_snps %||% "Inf"))
      if (is.na(cap) || !is.finite(cap) || nrow(d) <= cap) return(d)
      d[order(d$pval.exposure), , drop = FALSE][seq_len(cap), , drop = FALSE]
    }

    ## Uploaded-data branch: format + harmonise the user's own two files with
    ## TwoSampleMR (same functions, same defaults 10_MR.R itself uses), then
    ## hand off to the SAME estimate_mr_set() the project-data branch below
    ## uses, so the estimator hierarchy, its output shape and every
    ## downstream render (scatter, diagnostics, tables) is identical either
    ## way - only how `d` gets built differs.
    mr_result_uploaded <- function() {
      req(input$exp_file, input$out_file)
      req(input$exp_snp, input$exp_beta, input$exp_se, input$exp_pval, input$exp_ea, input$exp_oa)
      req(input$out_snp, input$out_beta, input$out_se, input$out_pval, input$out_ea, input$out_oa)

      exp_raw <- exp_df_r()
      out_raw <- out_df_r()
      validate(need(!is.null(exp_raw) && !is.null(out_raw), "Could not read one of the uploaded files."))

      label <- if (nzchar(trimws(input$upload_label %||% ""))) trimws(input$upload_label) else "Uploaded exposure"
      exp_fmt <- TwoSampleMR::format_data(
        exp_raw, type = "exposure", snp_col = input$exp_snp, beta_col = input$exp_beta,
        se_col = input$exp_se, pval_col = input$exp_pval, effect_allele_col = input$exp_ea,
        other_allele_col = input$exp_oa, eaf_col = if (nzchar(input$exp_eaf %||% "")) input$exp_eaf else "eaf"
      )
      out_fmt <- TwoSampleMR::format_data(
        out_raw, type = "outcome", snp_col = input$out_snp, beta_col = input$out_beta,
        se_col = input$out_se, pval_col = input$out_pval, effect_allele_col = input$out_ea,
        other_allele_col = input$out_oa, eaf_col = if (nzchar(input$out_eaf %||% "")) input$out_eaf else "eaf"
      )
      exp_fmt$exposure <- label
      out_fmt$outcome <- "Uploaded outcome"

      dat_up <- tryCatch(TwoSampleMR::harmonise_data(exp_fmt, out_fmt, action = 2), error = function(e) NULL)
      validate(need(!is.null(dat_up) && nrow(dat_up) > 0,
        "Harmonisation found no overlapping SNPs between the two files - check that both use the same SNP identifiers (e.g. rsIDs) and that allele columns are mapped correctly."))
      dat_up$gene <- label
      dat_up$Fstat <- (dat_up$beta.exposure / dat_up$se.exposure)^2
      dat_up$MHC <- FALSE
      dat_up <- dat_up[dat_up$mr_keep, , drop = FALSE]
      validate(need(nrow(dat_up) > 0, "No SNPs survived harmonisation (all were dropped as ambiguous/unresolvable, e.g. unresolved palindromic variants)."))

      n_before <- nrow(dat_up)
      pval_cut <- as.numeric(input$pval_cut %||% "5e-8")
      fstat_cut <- as.numeric(input$fstat_cut %||% "10")
      d <- dat_up[dat_up$pval.exposure <= pval_cut & dat_up$Fstat >= fstat_cut, , drop = FALSE]
      validate(need(nrow(d) >= 1, sprintf(
        "No SNPs remain after filtering (%d harmonised before filtering). Loosen the p-value or F-statistic threshold - uploaded data doesn't necessarily reach genome-wide significance the way the bundled eQTLGen instruments do.",
        n_before
      )))

      ## Optional live LD clumping (independence): off by default (needs
      ## network access to the OpenGWAS LD reference). A failure here
      ## (no token, unreachable, or no LD panel match) degrades to
      ## proceeding on the unclumped set rather than blocking the run - the
      ## clump_note flag records which happened so the Result panel can say
      ## so honestly instead of silently pretending it always succeeds.
      ## Distinguishes three outcomes explicitly rather than lumping them
      ## into one "failed" state: a genuine API/network/token failure (the
      ## call itself errors) is a different problem from the call
      ## succeeding but returning nothing usable (e.g. the uploaded rsIDs
      ## aren't genotyped in the chosen population's LD reference panel) -
      ## confirmed against this app's own test fixtures, where synthetic
      ## rsIDs correctly clump to zero survivors while real cached rsIDs
      ## clump successfully. Either way, the unclumped set is kept rather
      ## than silently discarding every instrument.
      clump_status <- "off"
      n_before_clump <- nrow(d)
      if (isTRUE(input$do_clump) && nrow(d) > 1) {
        clumped <- tryCatch(
          ieugwasr::ld_clump(
            data.frame(rsid = d$SNP, pval = d$pval.exposure, id = label),
            clump_kb = as.numeric(input$clump_kb %||% "10000"),
            clump_r2 = as.numeric(input$clump_r2 %||% "0.001"),
            clump_p = 1, pop = input$clump_pop %||% "EUR"
          ),
          error = function(e) NULL
        )
        if (is.null(clumped)) {
          clump_status <- "api_error"
        } else if (nrow(clumped) == 0) {
          clump_status <- "no_match"
        } else {
          d <- d[d$SNP %in% clumped$rsid, , drop = FALSE]
          clump_status <- if (nrow(d) < n_before_clump) "applied" else "applied_no_change"
        }
      }
      n_after_clump <- nrow(d)
      d <- apply_max_snps_cap(d)

      est <- estimate_mr_set(d, include_mode = isTRUE(input$include_mode),
                              ci_level = current_ci_level(), run_presso = isTRUE(input$run_presso))
      list(
        gene = label, d = d, n_before = n_before, n_after = nrow(d),
        n_mhc_dropped = 0L, n_mhc_retained = 0L, mhc_mode = "n/a",
        pval_cut = pval_cut, fstat_cut = fstat_cut, est = est,
        clump_status = clump_status, n_before_clump = n_before_clump, n_after_clump = n_after_clump,
        project_reference = mr_primary_all[0, ], uploaded = TRUE
      )
    }

    ## Parameterised on `gene` (not always input$gene) so the female/male
    ## tab pickers below can reuse the exact same bundled-dataset logic -
    ## same filters, same estimator hierarchy, same reference comparison -
    ## instead of a second, drift-prone copy.
    mr_result_for_gene <- function(gene) {
      req(gene)
      pval_cut <- as.numeric(input$pval_cut %||% "5e-8")
      fstat_cut <- as.numeric(input$fstat_cut %||% "10")
      mhc_mode <- input$mhc_mode %||% "include"

      d0 <- mr_dat_all[mr_dat_all$gene == gene, , drop = FALSE]
      validate(need(nrow(d0) > 0, "No harmonised instrument rows are cached for this gene."))
      n_before <- nrow(d0)

      d <- d0[d0$pval.exposure <= pval_cut & d0$Fstat >= fstat_cut, , drop = FALSE]
      n_mhc_in_filtered <- sum(d$MHC)
      if (identical(mhc_mode, "exclude")) d <- d[!d$MHC, , drop = FALSE]

      validate(need(nrow(d) >= 1, sprintf(
        "No instruments remain for %s after these filters (%d were available before filtering). Loosen the p-value or F-statistic threshold, or switch MHC handling back to \"Include\".",
        gene, n_before
      )))
      d <- apply_max_snps_cap(d)

      est <- estimate_mr_set(d, include_mode = isTRUE(input$include_mode),
                              ci_level = current_ci_level(), run_presso = isTRUE(input$run_presso))

      list(
        gene = gene, d = d, n_before = n_before, n_after = nrow(d),
        n_mhc_dropped = if (identical(mhc_mode, "exclude")) n_mhc_in_filtered else 0L,
        n_mhc_retained = if (identical(mhc_mode, "include")) n_mhc_in_filtered else 0L,
        mhc_mode = mhc_mode, pval_cut = pval_cut, fstat_cut = fstat_cut,
        est = est,
        project_reference = mr_primary_all[mr_primary_all$gene == gene, ], uploaded = FALSE
      )
    }
    mr_result_project <- function() mr_result_for_gene(input$gene)

    mr_result <- eventReactive(input$run_btn, {
      if (identical(input$data_source, "upload")) mr_result_uploaded() else mr_result_project()
    })

    ## Nothing about a run's results - Result/scatter/diagnostics/instrument
    ## table/MR-estimates, box titles included - exists on the page until
    ## "Run MR" has actually been clicked at least once. A validate() error
    ## inside mr_result() (e.g. no instruments left after filtering) still
    ## surfaces normally once revealed, same as every other tab in this app.
    mr_has_run <- reactiveVal(FALSE)
    observeEvent(input$run_btn, mr_has_run(TRUE))

    output$result_box_ui <- renderUI({
      if (!mr_has_run()) return(NULL)
      box(
        width = NULL, title = "Result", status = "primary", solidHeader = FALSE,
        withSpinner(uiOutput(ns("summary_ui")), color = "#2c6fbb", type = 6)
      )
    })

    output$scatter_diag_ui <- renderUI({
      if (!mr_has_run()) {
        ## Empty state, not a result: no box title, no numbers, nothing
        ## this app has actually computed - just fills the column next to
        ## the (now taller) filters panel so the page doesn't read as
        ## broken while still showing nothing result-related pre-click.
        return(box(
          width = NULL, solidHeader = FALSE,
          div(class = "coming-soon", style = "padding: 90px 20px;",
              icon("route", class = "coming-soon-icon"),
              h4("Nothing run yet"),
              p("Pick a gene (or upload your own data) on the left, adjust filters if you like, and click “Run MR” - the SNP effect scatter, sensitivity diagnostics and result will appear here."))
        ))
      }
      tagList(
        box(
          width = NULL, title = "SNP effect scatter", status = "primary", solidHeader = FALSE,
          withSpinner(plotOutput(ns("scatter_plot"), height = 380), color = "#2c6fbb", type = 6)
        ),
        box(
          width = NULL, title = "Sensitivity diagnostics", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Funnel plot (asymmetry suggests directional pleiotropy) and leave-one-out IVW (whether one SNP drives the whole estimate) - both only defined with ≥3 instruments, same as Cochran's Q and the MR-Egger intercept shown in the Result panel."),
          uiOutput(ns("diagnostics_ui"))
        )
      )
    })

    output$instrument_mr_tables_ui <- renderUI({
      if (!mr_has_run()) return(NULL)
      tagList(
        box(
          width = 12, title = "Instruments used in this run", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "SNPs surviving your current filters, after harmonisation against the outcome."),
          DT::dataTableOutput(ns("instrument_table"))
        ),
        box(
          width = 12, title = "MR estimates", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Every estimator this gene's instrument count qualifies for. “Project primary” marks the pre-specified hierarchy's pick (IVW > Wald ratio > weighted median > MR-Egger) - fixed in advance, never chosen by which p-value looks best."),
          div(class = "table-toolbar", downloadButton(ns("download_mr"), "Download CSV", class = "btn-sm")),
          DT::dataTableOutput(ns("mr_table"))
        )
      )
    })

    observeEvent(mr_result(), {
      res <- mr_result()
      primary_row <- res$est$res_table[res$est$res_table$primary, , drop = FALSE][1, ]
      entry <- list(
        gene = res$gene, method = primary_row$method, n_snp = res$est$n_snp,
        estimate = round(primary_row$estimate, 4), p = signif(primary_row$p, 3),
        filters = sprintf("eQTL p<=%.0e, F>=%.0f, MHC %s", res$pval_cut, res$fstat_cut,
                           if (identical(res$mhc_mode, "exclude")) "excluded" else "included")
      )
      existing <- results$mr %||% list(genes_tested = list())
      existing$genes_tested[[res$gene]] <- entry
      results$mr <- existing
    })

    output$summary_ui <- renderUI({
      res <- mr_result()
      pr <- res$project_reference
      het <- res$est$heterogeneity
      pleio <- res$est$pleiotropy
      presso <- res$est$presso
      alpha <- as.numeric(input$alpha_cut %||% 0.05)
      primary_row <- res$est$res_table[res$est$res_table$primary, , drop = FALSE][1, ]
      is_sig <- !is.na(primary_row$p) && primary_row$p < alpha

      tagList(
        if (isTRUE(res$uploaded)) p(class = "empty-note", icon("upload"),
          "Result from your own uploaded exposure/outcome files - no MHC flag or bundled reference estimate applies to custom data."),
        if (isTRUE(res$uploaded) && identical(res$clump_status, "applied"))
          p(class = "empty-note", icon("circle-check"), sprintf("LD clumping removed %d correlated SNP%s (%d -> %d instruments).", res$n_before_clump - res$n_after_clump, if ((res$n_before_clump - res$n_after_clump) != 1) "s" else "", res$n_before_clump, res$n_after_clump)),
        if (isTRUE(res$uploaded) && identical(res$clump_status, "applied_no_change"))
          p(class = "empty-note", icon("circle-check"), "LD clumping ran and found no correlated SNPs to remove."),
        if (isTRUE(res$uploaded) && identical(res$clump_status, "api_error"))
          p(class = "empty-note", icon("triangle-exclamation"), "LD clumping could not reach the OpenGWAS LD reference (network or token issue) - proceeded on the unclumped SNP set."),
        if (isTRUE(res$uploaded) && identical(res$clump_status, "no_match"))
          p(class = "empty-note", icon("triangle-exclamation"), "LD clumping ran but found none of these SNP IDs in the chosen population's LD reference panel (check they're rsIDs, and that the population matches your data's ancestry) - proceeded on the unclumped SNP set."),
        p(strong(res$gene), ": ", strong(res$est$n_snp),
          paste0(" instrument SNP", if (res$est$n_snp != 1) "s" else "", " after your current filters (",
                 format(res$n_before, big.mark = ","), if (isTRUE(res$uploaded)) " harmonised" else " cached", " before filtering).")),
        if (identical(res$mhc_mode, "exclude") && res$n_mhc_dropped > 0)
          p(class = "empty-note", icon("circle-info"), sprintf("%d MHC-region instrument%s excluded for this gene (sensitivity analysis).", res$n_mhc_dropped, if (res$n_mhc_dropped != 1) "s" else "")),
        if (identical(res$mhc_mode, "include") && res$n_mhc_retained > 0)
          p(class = "empty-note", icon("triangle-exclamation"), sprintf("%d of these instrument%s sit in the extended MHC - interpret this estimate as associated, not necessarily causal, unless it also survives with MHC instruments excluded.", res$n_mhc_retained, if (res$n_mhc_retained != 1) "s" else "")),
        p(strong("Primary estimate: "), primary_row$method, sprintf(
          ", b = %.4f (%.0f%% CI %.4f to %.4f), p = %s",
          primary_row$estimate, res$est$ci_level * 100, primary_row$ci_low, primary_row$ci_high, format(signif(primary_row$p, 3), scientific = TRUE)
        ), " ", tags$span(style = sprintf("color:%s; font-weight:600;", if (is_sig) ARTHOMIX_STATUS$good else "#8A929C"),
                            sprintf("(%s at α = %.3g)", if (is_sig) "significant" else "not significant", alpha))),
        if (!is.null(het)) p(strong("Cochran's Q: "), sprintf("%.2f on %d df, p = %s", het$Q, het$Q_df, format(signif(het$Q_pval, 3), scientific = TRUE)),
                              " (heterogeneity among instruments - a low p suggests some instrument disagrees with the rest)."),
        if (!is.null(pleio)) p(strong("MR-Egger intercept: "), sprintf("%.4f, p = %s, I² = %.3f", pleio$intercept, format(signif(pleio$p, 3), scientific = TRUE), pleio$I2),
                                " (a p < 0.05 intercept suggests directional pleiotropy; I² well below 0.9 means the Egger slope itself is less reliable - NOME assumption)."),
        if (!is.null(presso) && is.null(presso$error)) tagList(
          p(strong("MR-PRESSO global test: "), sprintf("p = %s", format(presso$global_p, scientific = TRUE)),
            sprintf(" (%d outlier%s flagged)", presso$n_outliers, if (presso$n_outliers != 1) "s" else ""),
            " - a significant global test indicates horizontal pleiotropy somewhere in the instrument set."),
          if (!is.na(presso$corrected_estimate))
            p(strong("MR-PRESSO outlier-corrected estimate: "), sprintf("b = %.4f, p = %s", presso$corrected_estimate, format(signif(presso$corrected_p, 3), scientific = TRUE)),
              " (re-estimated with the flagged outlier instrument(s) removed).")
        ),
        if (!is.null(presso) && !is.null(presso$error))
          p(class = "empty-note", icon("triangle-exclamation"), sprintf("MR-PRESSO could not run: %s", presso$error)),
        if (res$est$n_snp < 3) p(class = "empty-note", icon("circle-info"),
                              "Fewer than 3 SNPs: only the IVW (or Wald ratio) estimate is shown - MR-Egger, weighted median, Cochran's Q and the MR-Egger intercept all need at least 3 for a meaningful pleiotropy check."),
        if (res$est$n_snp >= 3 && res$est$n_snp < 4 && isTRUE(input$run_presso))
          p(class = "empty-note", icon("circle-info"), "MR-PRESSO needs at least 4 instruments; not run for this gene."),
        if (res$gene %in% relabel_check$genes) p(class = "empty-note", icon("circle-info"),
                              "This gene shares at least one instrument SNP with a neighbouring gene, and the cached reference estimate shown below was computed before this app's own SNP-to-gene relabelling fix - your live estimate above may therefore use a slightly different instrument set than the cached one."),
        if (nrow(pr) > 0) tagList(
          hr(),
          p(class = "submodule-desc", style = "margin-bottom: 0;",
            sprintf("The bundled reference pipeline (before separate female/male FDR correction) estimated this gene at b = %.4f, p = %s, using %d instrument%s%s.",
                    pr$b[1], format(signif(pr$pval[1], 3), scientific = TRUE), pr$nSNP[1], if (pr$nSNP[1] != 1) "s" else "",
                    if (isTRUE(pr$MHC_gene[1])) " (MHC-flagged)" else ""))
        )
      )
    })

    ## ---------------------------------------------------------------------
    ## Instrument table: exactly the rows estimate_mr_set() above was fit
    ## on for this gene, so what a user sees here is never out of step with
    ## the estimate below it.
    ## ---------------------------------------------------------------------
    output$instrument_table <- DT::renderDataTable({
      res <- mr_result()
      ## Uploaded data has no chr/pos (format_data() wasn't given chr_col/
      ## pos_col, since this app doesn't ask the user to map them) and no
      ## real MHC flag (always FALSE, set in mr_result_uploaded) - build the
      ## column list from what's actually present instead of assuming the
      ## project-cached shape.
      cols <- c(SNP = "SNP", Chr = "chr.exposure", `Position (GRCh37)` = "pos.exposure",
                `β exposure` = "beta.exposure", `SE exposure` = "se.exposure",
                `eQTL/exposure p-value` = "pval.exposure", `F-statistic` = "Fstat", `In MHC` = "MHC")
      cols <- cols[cols %in% colnames(res$d)]
      d <- res$d[, cols, drop = FALSE]
      colnames(d) <- names(cols)
      if ("eQTL/exposure p-value" %in% colnames(d)) d[["eQTL/exposure p-value"]] <- signif(d[["eQTL/exposure p-value"]], 3)
      if ("β exposure" %in% colnames(d)) d[["β exposure"]] <- round(d[["β exposure"]], 4)
      if ("SE exposure" %in% colnames(d)) d[["SE exposure"]] <- round(d[["SE exposure"]], 4)
      if ("F-statistic" %in% colnames(d)) d[["F-statistic"]] <- round(d[["F-statistic"]], 1)
      if ("In MHC" %in% colnames(d)) d[["In MHC"]] <- ifelse(d[["In MHC"]], "Yes", "No")
      DT::datatable(d, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE, dom = "tip"), class = "stripe hover compact")
    })

    output$mr_table <- DT::renderDataTable({
      df <- mr_result()$est$res_table
      df$estimate <- round(df$estimate, 4)
      df$se <- round(df$se, 4)
      df$ci_low <- round(df$ci_low, 4)
      df$ci_high <- round(df$ci_high, 4)
      df$p <- signif(df$p, 3)
      df$primary <- ifelse(df$primary, "✓", "")
      colnames(df) <- c("Method", "Estimate (b)", "SE", "95% CI low", "95% CI high", "p-value", "Project primary")
      DT::datatable(df, rownames = FALSE, options = list(pageLength = 5, scrollX = TRUE, dom = "t"), class = "stripe hover compact")
    })

    scatter_plot_obj <- reactive({
      res <- mr_result()
      d <- res$d
      primary_row <- res$est$res_table[res$est$res_table$primary, , drop = FALSE][1, ]
      d$mhc_label <- ifelse(d$MHC, "MHC region", "Non-MHC")
      p <- ggplot(d, aes(x = beta.exposure, y = beta.outcome, color = mhc_label)) +
        geom_hline(yintercept = 0, color = "#d5d9de") +
        geom_vline(xintercept = 0, color = "#d5d9de") +
        geom_errorbar(aes(ymin = beta.outcome - 1.96 * se.outcome, ymax = beta.outcome + 1.96 * se.outcome), width = 0, alpha = 0.6) +
        geom_errorbarh(aes(xmin = beta.exposure - 1.96 * se.exposure, xmax = beta.exposure + 1.96 * se.exposure), height = 0, alpha = 0.6) +
        geom_point(size = 2.6) +
        geom_abline(intercept = 0, slope = primary_row$estimate, color = "#c0392b", linetype = "dashed") +
        scale_color_manual(values = c("Non-MHC" = ARTHOMIX_COLORS$blue, "MHC region" = ARTHOMIX_COLORS$red), name = NULL) +
        labs(x = "SNP effect on expression (cis-eQTL beta)", y = "SNP effect on RA risk (GWAS beta)",
             title = paste0(res$gene, " (", res$est$n_snp, " SNP", if (res$est$n_snp != 1) "s" else "", ", ", primary_row$method, " shown)")) +
        theme_minimal(base_size = 12) + theme(legend.position = "top")
      p
    })

    output$scatter_plot <- renderPlot({
      p <- tryCatch(scatter_plot_obj(), error = function(e) NULL)
      req(p)
      p
    })

    ## ---------------------------------------------------------------------
    ## Leave-one-out: re-fit IVW with each SNP dropped in turn, plus the
    ## all-SNP estimate as the reference row - the same diagnostic
    ## scripts/goal2_sex_stratified/35_figure_mr_diagnostics.R draws for
    ## the >=3-instrument subset of the thesis's own prioritised genes.
    ## ---------------------------------------------------------------------
    loo_result <- reactive({
      res <- mr_result()
      d <- res$d
      if (nrow(d) < 3) return(NULL)
      fit_ivw <- function(dd) {
        mrobj <- MendelianRandomization::mr_input(bx = dd$beta.exposure, bxse = dd$se.exposure, by = dd$beta.outcome, byse = dd$se.outcome, snps = dd$SNP)
        MendelianRandomization::mr_ivw(mrobj, model = "random")
      }
      rows <- lapply(seq_len(nrow(d)), function(i) {
        ivw_i <- fit_ivw(d[-i, , drop = FALSE])
        data.frame(snp = paste("Excl.", d$SNP[i]), estimate = ivw_i@Estimate, ci_low = ivw_i@CILower, ci_high = ivw_i@CIUpper)
      })
      ivw_all <- fit_ivw(d)
      out <- do.call(rbind, rows)
      out <- rbind(out, data.frame(snp = "All SNPs (IVW)", estimate = ivw_all@Estimate, ci_low = ivw_all@CILower, ci_high = ivw_all@CIUpper))
      out$snp <- factor(out$snp, levels = rev(out$snp))
      out
    })

    ## Per-SNP Wald ratio (beta.outcome/beta.exposure) against its own
    ## precision (1/SE of that ratio) - the standard MR funnel-plot
    ## construction: symmetric scatter around the pooled IVW estimate is
    ## consistent with no directional pleiotropy, a lopsided funnel is not.
    funnel_data <- reactive({
      res <- mr_result()
      d <- res$d
      if (nrow(d) < 3) return(NULL)
      ratio <- d$beta.outcome / d$beta.exposure
      ratio_se <- abs(d$se.outcome / d$beta.exposure)
      data.frame(SNP = d$SNP, ratio = ratio, precision = 1 / ratio_se)
    })

    output$diagnostics_ui <- renderUI({
      res <- mr_result()
      if (res$est$n_snp < 3) {
        return(p(class = "empty-note", icon("circle-info"),
                  sprintf("%s has %d instrument%s after your current filters - fewer than the 3 needed for a funnel plot or leave-one-out.",
                          res$gene, res$est$n_snp, if (res$est$n_snp != 1) "s" else "")))
      }
      fluidRow(
        column(6, withSpinner(plotOutput(ns("funnel_plot"), height = 320), color = "#2c6fbb", type = 6)),
        column(6, withSpinner(plotOutput(ns("loo_plot"), height = 320), color = "#2c6fbb", type = 6))
      )
    })

    output$funnel_plot <- renderPlot({
      fd <- funnel_data()
      req(fd)
      primary_row <- mr_result()$est$res_table[mr_result()$est$res_table$primary, , drop = FALSE][1, ]
      ggplot(fd, aes(x = ratio, y = precision)) +
        geom_vline(xintercept = primary_row$estimate, linetype = "dashed", color = "#c0392b") +
        geom_point(color = ARTHOMIX_COLORS$blue, size = 2.4, alpha = 0.8) +
        labs(x = "Per-SNP causal estimate (Wald ratio)", y = "Precision (1 / SE)", title = "Funnel plot") +
        theme_minimal(base_size = 12)
    })

    output$loo_plot <- renderPlot({
      out <- loo_result()
      req(out)
      ggplot(out, aes(x = estimate, y = snp)) +
        geom_vline(xintercept = 0, color = "#d5d9de") +
        geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.15, color = "#a9b1bb") +
        geom_point(aes(color = snp == "All SNPs (IVW)"), size = 2.6, show.legend = FALSE) +
        scale_color_manual(values = c(`TRUE` = ARTHOMIX_COLORS$red, `FALSE` = ARTHOMIX_COLORS$blue)) +
        labs(x = "IVW estimate (95% CI)", y = NULL, title = "Leave-one-out") +
        theme_minimal(base_size = 12)
    })

    output$download_mr <- downloadHandler(
      filename = function() paste0("mr_", mr_result()$gene, ".csv"),
      content = function(file) write.csv(mr_result()$est$res_table, file, row.names = FALSE)
    )

    ## No box, no title, until "Run batch screen" is clicked - this box's
    ## own button is its only entry point (unlike Result/scatter/
    ## diagnostics, which are triggered by the separate "Run MR" button in
    ## the "Instrument filters" box).
    batch_has_run <- reactiveVal(FALSE)
    observeEvent(input$run_batch_btn, batch_has_run(TRUE))

    output$batch_screen_section <- renderUI({
      ## Bundled-dataset only (it loops over the ~1,700 cached genes), so it
      ## is entirely absent - not just disabled - when "MR overall" is set
      ## to upload mode, matching that tab's own data-source switch.
      if (!identical(input$data_source, "project")) return(NULL)
      if (!batch_has_run()) {
        return(div(
          style = "margin: 4px 0 16px 0;",
          p(class = "submodule-desc", style = "display:inline; margin-right: 8px;",
            "Want FDR-corrected estimates across every available gene under your current filters, not just one (advanced, bundled-dataset mode only)?"),
          actionButton(ns("run_batch_btn"), "Run batch screen", icon = icon("bolt"), class = "btn-default btn-sm")
        ))
      }
      box(
        width = 12, title = "Batch screen: MR across every available gene", status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "Applies the filters above to every gene with a cached instrument (up to 1,701), then BH-FDR corrects across genes tested under these filters - this session's own screen, not a separate, precomputed per-sex FDR reference, which used a different denominator."),
        withSpinner(uiOutput(ns("batch_summary_ui")), color = "#2c6fbb", type = 6),
        div(class = "table-toolbar", downloadButton(ns("download_batch"), "Download CSV", class = "btn-sm")),
        DT::dataTableOutput(ns("batch_table"))
      )
    })

    ## ---------------------------------------------------------------------
    ## Batch screen: the same filters and the same estimate_mr_set() used
    ## for one gene above, run across every gene with a cached instrument,
    ## then BH-FDR corrected across whatever set of genes the current
    ## filters actually leave testable - this session's own screen, kept
    ## explicitly distinct from the thesis's separate, sex-stratified
    ## denominators (Section 2.6.7) in every label this produces.
    ## ---------------------------------------------------------------------
    batch_result <- eventReactive(input$run_batch_btn, {
      validate(need(identical(input$data_source, "project"),
        "The batch screen runs across the bundled dataset's genes and isn't available in \"Upload your own GWAS\" mode - switch the data source above to the bundled dataset first."))
      pval_cut <- as.numeric(input$pval_cut %||% "5e-8")
      fstat_cut <- as.numeric(input$fstat_cut %||% "10")
      mhc_mode <- input$mhc_mode %||% "include"

      d_all <- mr_dat_all[mr_dat_all$pval.exposure <= pval_cut & mr_dat_all$Fstat >= fstat_cut, , drop = FALSE]
      if (identical(mhc_mode, "exclude")) d_all <- d_all[!d_all$MHC, , drop = FALSE]
      genes <- sort(unique(d_all$gene))
      validate(need(length(genes) > 0, "No genes have any instrument left after these filters."))

      split_d <- split(d_all, d_all$gene)
      rows <- vector("list", length(genes))
      withProgress(message = "Running MR across all available genes", value = 0, {
        for (i in seq_along(genes)) {
          dg <- apply_max_snps_cap(split_d[[genes[i]]])
          est <- tryCatch(estimate_mr_set(dg, include_mode = FALSE, full = FALSE, ci_level = current_ci_level()), error = function(e) NULL)
          if (!is.null(est)) {
            prim <- est$res_table[est$res_table$primary, , drop = FALSE][1, ]
            rows[[i]] <- data.frame(
              gene = genes[i], method = prim$method, nSNP = est$n_snp,
              estimate = prim$estimate, se = prim$se, p = prim$p,
              OR = exp(prim$estimate), OR_lo = exp(prim$ci_low), OR_hi = exp(prim$ci_high),
              MHC_gene = any(dg$MHC), sensitivity_testable = est$n_snp >= 3,
              stringsAsFactors = FALSE
            )
          }
          if (i %% 25 == 0 || i == length(genes)) setProgress(value = i / length(genes), detail = sprintf("%d / %d genes", i, length(genes)))
        }
      })
      out <- do.call(rbind, rows)
      out$FDR <- stats::p.adjust(out$p, method = "BH")
      out <- out[order(out$p), , drop = FALSE]
      list(table = out, n_genes = length(genes), pval_cut = pval_cut, fstat_cut = fstat_cut, mhc_mode = mhc_mode)
    })

    observeEvent(batch_result(), {
      br <- batch_result()
      alpha <- as.numeric(input$alpha_cut %||% 0.05)
      existing <- results$mr %||% list(genes_tested = list())
      existing$batch_screen <- list(
        n_genes_tested = br$n_genes,
        n_fdr_significant = sum(br$table$FDR < alpha, na.rm = TRUE),
        filters = sprintf("eQTL p<=%.0e, F>=%.0f, MHC %s, alpha=%.3g", br$pval_cut, br$fstat_cut,
                           if (identical(br$mhc_mode, "exclude")) "excluded" else "included", alpha)
      )
      results$mr <- existing
    })

    output$batch_summary_ui <- renderUI({
      br <- batch_result()
      alpha <- as.numeric(input$alpha_cut %||% 0.05)
      n_sig <- sum(br$table$FDR < alpha, na.rm = TRUE)
      p(strong(format(br$n_genes, big.mark = ",")), " genes tested under the current filters, ",
        strong(n_sig), sprintf(" significant at FDR < %.3g.", alpha))
    })

    output$batch_table <- DT::renderDataTable({
      df <- batch_result()$table
      df$estimate <- round(df$estimate, 4); df$se <- round(df$se, 4)
      df$p <- signif(df$p, 3); df$FDR <- signif(df$FDR, 3)
      df$OR <- round(df$OR, 3); df$OR_lo <- round(df$OR_lo, 3); df$OR_hi <- round(df$OR_hi, 3)
      df$MHC_gene <- ifelse(df$MHC_gene, "Yes", "No")
      df$sensitivity_testable <- ifelse(df$sensitivity_testable, "Yes", "No")
      colnames(df) <- c("Gene", "Method", "nSNP", "b", "SE", "p", "OR", "OR low", "OR high", "In MHC", "≥3 SNP (pleiotropy testable)", "FDR")
      DT::datatable(df, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_batch <- downloadHandler(
      filename = function() "mr_batch_screen.csv",
      content = function(file) write.csv(batch_result()$table, file, row.names = FALSE)
    )
  })
}
