## R/methylomics/11_Mendelian_Randomization/mod_methyl_mr.R
## Two-sample MR for methylation exposures: cis/trans mQTL instruments for
## a CpG -> GWAS outcome. Separate from R/transcriptomics/mod_mr.R (eQTL MR).
##
## Two data routes: Preloaded reproduces script08_mendelian_randomization's
## completed run (GoDMC cis-mQTL, already clumped/harmonised against the
## Ishigaki et al. 2022 RA GWAS) - MR estimation/sensitivity/plots run live
## from the cached harmonised table. Upload Dataset runs the full live
## pipeline (instrument selection, ld_clump() via OpenGWAS API,
## harmonisation, MR, sensitivity) on user-supplied summary stats.
##
## Stage-gated UI: Data -> Filters & Instruments -> LD Clumping ->
## Harmonisation -> MR Analysis -> Sensitivity -> Results -> Plots; changing
## a stage's inputs invalidates everything downstream (has-run reactiveVal
## pattern, same as mod_methyl_dmp.R/mod_methyl_normalization.R).

## ---------------------------------------------------------------------------
## Small local helpers (not shared with R/transcriptomics/mod_mr.R)
## ---------------------------------------------------------------------------

.mmr_tip <- function(text) tags$span(icon("circle-info", style = "color:#8A929C; cursor: help; margin-left: 4px;"), title = text)

.mmr_stage_order <- c("data", "instruments", "clump", "harmonise", "mr", "sensitivity")

## GoDMC/script08's cis-mQTL definition (METHODS_mendelian_randomization.md
## Section 3); pre-filled default for the Upload route's cis window,
## echoed read-only for the Preloaded route.
MMR_DEFAULT_CIS_WINDOW_BP <- 1e6
MMR_DEFAULT_CIS_PVAL <- 5e-8
MMR_DEFAULT_MIN_F <- 10
MMR_DEFAULT_CLUMP_R2 <- 0.001
MMR_DEFAULT_CLUMP_KB <- 10000

mod_methyl_mr_config <- list(
  id = "mr", title = "Mendelian Randomization", icon = "route", group = "Genetics",
  description = "Two-sample Mendelian randomisation of mQTL instruments against a GWAS outcome."
)

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------

mod_methyl_mr_ui <- function(id) {
  ns <- NS(id)
  div(
    class = "tx-menu-wrap",
    tabsetPanel(
      id = ns("mr_tabs"), type = "tabs",
      tabPanel("1. Data", br(), mmr_data_ui(ns)),
      tabPanel("2. Filters & Instruments", br(), mmr_filters_ui(ns)),
      tabPanel("3. LD Clumping", br(), mmr_clump_ui(ns)),
      tabPanel("4. Harmonisation", br(), mmr_harmonise_ui(ns)),
      tabPanel("5. MR Analysis", br(), mmr_analysis_ui(ns)),
      tabPanel("6. Sensitivity", br(), mmr_sensitivity_ui(ns)),
      tabPanel("7. Results", br(), mmr_results_ui(ns)),
      tabPanel("8. Plots", br(), mmr_plots_ui(ns))
    )
  )
}

## ---- Tab 1: Data -----------------------------------------------------------

mmr_data_ui <- function(ns) {
  fluidRow(
    column(
      4,
      box(
        width = NULL, title = "Data source", status = "primary", solidHeader = FALSE,
        radioButtons(
          ns("data_source"), NULL,
          choiceNames = list(
            tagList(icon("upload"), " Upload Dataset"),
            tagList(icon("database"), " Use Preloaded Data")
          ),
          choiceValues = list("upload", "preloaded"), selected = "preloaded"
        )
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'preloaded'", ns("data_source")),
        box(
          width = NULL, title = "Preloaded methylomics MR dataset", status = "primary", solidHeader = FALSE,
          if (METH_DATA_AVAILABLE) tagList(
            p(class = "submodule-desc", "Reproduces ", strong("script08_mendelian_randomization"), "."),
            radioButtons(ns("pre_sex"), "Sex stratum", inline = TRUE,
                         choices = c("Female" = "female", "Male" = "male", "Both (union)" = "both"), selected = "female"),
            uiOutput(ns("pre_cpg_ui")),
            checkboxInput(ns("pre_binary_outcome"), "Outcome is binary (rheumatoid arthritis case/control - show odds ratios)", value = TRUE)
          ) else p(class = "empty-note", icon("triangle-exclamation"),
                   "The preloaded methylomics MR dataset isn't available in this deployment - use Upload Dataset instead.")
        )
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'upload'", ns("data_source")),
        box(
          width = NULL, title = "Upload exposure (mQTL) data", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "One row per SNP-CpG association. CSV/TSV."),
          fileInput(ns("exp_file"), "Exposure / methylation-QTL file", accept = c(".csv", ".tsv", ".txt")),
          uiOutput(ns("exp_map_ui")),
          uiOutput(ns("exp_extra_map_ui"))
        ),
        box(
          width = NULL, title = "Upload outcome GWAS data", status = "primary", solidHeader = FALSE,
          textInput(ns("outcome_label"), "Outcome name (for labelling only)", value = "Uploaded outcome", width = "100%"),
          checkboxInput(ns("up_binary_outcome"), "Outcome is a binary trait (show odds ratios)", value = FALSE),
          fileInput(ns("out_file"), "Outcome GWAS file", accept = c(".csv", ".tsv", ".txt")),
          uiOutput(ns("out_map_ui"))
        )
      )
    ),
    column(
      8,
      uiOutput(ns("data_validation_ui")),
      div(style = "margin-top: 10px;",
          actionButton(ns("validate_btn"), "Validate Data", icon = icon("check"), class = "btn-primary btn-sm"))
    )
  )
}

## ---- Tab 2: Filters & Instruments ------------------------------------------

mmr_filters_ui <- function(ns) {
  uiOutput(ns("filters_tab_body"))
}

mmr_filters_controls <- function(ns) {
  tagList(
    div(style = "display:flex; align-items:center; gap:4px;",
        numericInput(ns("f_pval"), "Instrument p-value threshold", value = MMR_DEFAULT_CIS_PVAL, min = 0, max = 1, step = 1e-9, width = "100%"),
        .mmr_tip("Default 5x10^-8 (conventional genome-wide significance), matching this project's own GoDMC instrument threshold. Editable, not hardcoded.")),
    radioButtons(ns("f_region_mode"), "Instrument region", inline = TRUE,
                 choices = c("All instruments" = "all", "Cis-mQTL only (default)" = "cis", "Trans allowed" = "trans"), selected = "cis"),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'cis'", ns("f_region_mode")),
      div(style = "display:flex; align-items:center; gap:4px;",
          numericInput(ns("f_cis_window"), "Cis window (+/- kb from CpG)", value = MMR_DEFAULT_CIS_WINDOW_BP / 1000, min = 1, step = 100, width = "100%"),
          .mmr_tip("Default +/-1,000kb, this project's own cis-mQTL window (script08_mendelian_randomization). SNPs are cis if within this distance of the CpG on the same chromosome."))
    ),
    div(style = "display:flex; align-items:center; gap:4px;",
        numericInput(ns("f_maf"), "MAF / EAF threshold", value = 0, min = 0, max = 0.5, step = 0.01, width = "100%"),
        .mmr_tip("Minimum minor allele frequency (min(EAF, 1-EAF)). Only applied to instruments with an EAF value - rows missing EAF are flagged, not silently dropped, unless you also require EAF below.")),
    checkboxInput(ns("f_require_eaf"), "Require EAF to be present", value = FALSE),
    div(style = "display:flex; align-items:center; gap:4px;",
        numericInput(ns("f_min_instruments"), "Minimum instruments per CpG", value = 1, min = 1, step = 1, width = "100%"),
        numericInput(ns("f_max_instruments"), "Maximum instruments per CpG (strongest by p-value)", value = NA, min = 1, step = 1, width = "100%")),
    div(style = "display:flex; align-items:center; gap:4px;",
        numericInput(ns("f_min_f"), "Minimum F-statistic", value = MMR_DEFAULT_MIN_F, min = 0, step = 1, width = "100%"),
        .mmr_tip("Default F >= 10, the conventional weak-instrument threshold and this project's own MIN_F_STAT. Weak instruments are flagged, not silently removed - see the retained/excluded breakdown below.")),
    checkboxInput(ns("f_exclude_weak"), "Exclude weak instruments (F below threshold) from downstream MR", value = TRUE),
    div(style = "margin-top: 6px;",
        actionButton(ns("instruments_btn"), "Select Instruments", icon = icon("filter"), class = "btn-primary btn-sm"))
  )
}

## ---- Tab 3: LD Clumping -----------------------------------------------------

mmr_clump_ui <- function(ns) {
  uiOutput(ns("clump_tab_body"))
}

mmr_clump_controls <- function(ns) {
  tagList(
    div(class = "empty-note", icon("circle-info"),
        "Only OpenGWAS-API-based clumping (ieugwasr::ld_clump()) is available in this deployment - no local PLINK binary/LD reference is bundled, so \"Local PLINK reference\" is shown disabled below."),
    fluidRow(
      column(6, selectInput(ns("clump_pop"), "LD population", c("European (EUR)" = "EUR", "African (AFR)" = "AFR", "Admixed American (AMR)" = "AMR", "East Asian (EAS)" = "EAS", "South Asian (SAS)" = "SAS"), selected = "EUR", width = "100%")),
      column(6, div(style = "opacity:0.5;", checkboxInput(ns("clump_local_plink"), "Use local PLINK reference (unavailable)", value = FALSE)))
    ),
    fluidRow(
      column(4, numericInput(ns("clump_r2"), "r^2 threshold", value = MMR_DEFAULT_CLUMP_R2, min = 0, max = 1, step = 0.001, width = "100%")),
      column(4, numericInput(ns("clump_kb"), "Window (kb)", value = MMR_DEFAULT_CLUMP_KB, min = 1, step = 100, width = "100%")),
      column(4, numericInput(ns("clump_p"), "p-value for index SNP", value = 1, min = 0, max = 1, step = 0.01, width = "100%"))
    ),
    div(style = "margin-top: 6px;",
        actionButton(ns("clump_btn"), "Run Clumping", icon = icon("scissors"), class = "btn-primary btn-sm"))
  )
}

## ---- Tab 4: Harmonisation ----------------------------------------------------

mmr_harmonise_ui <- function(ns) {
  uiOutput(ns("harmonise_tab_body"))
}

mmr_harmonise_controls <- function(ns, editable = TRUE) {
  tagList(
    radioButtons(
      ns("harmonise_action"), "Harmonisation strategy",
      choiceNames = list(
        "1. Assume forward-strand alleles",
        "2. Infer strand using allele frequency (default)",
        "3. Correct non-palindromic variants and remove palindromic variants"
      ),
      choiceValues = list("1", "2", "3"), selected = "2"
    ),
    if (!editable) div(class = "empty-note", icon("circle-info"),
      "The Preloaded route's cached instrument set was harmonised with strategy 2 (script08_mendelian_randomization's default) - switch to Upload Dataset to harmonise your own data with a different strategy."),
    div(style = "margin-top: 6px;",
        actionButton(ns("harmonise_btn"), "Harmonise", icon = icon("link"), class = "btn-primary btn-sm"))
  )
}

## ---- Tab 5: MR Analysis ------------------------------------------------------

mmr_analysis_ui <- function(ns) {
  uiOutput(ns("analysis_tab_body"))
}

## MR-RAPS / penalised weighted median are offered only when exists()
## confirms the installed TwoSampleMR version provides them. Contamination
## mixture isn't provided by this version, and Radial MR (RadialMR) uses
## its own API rather than mr()'s method_list - neither is offered here.
mmr_robust_method_choices <- function() {
  extra <- list()
  if (exists("mr_raps", where = asNamespace("TwoSampleMR"), inherits = FALSE)) extra[["MR-RAPS (robust)"]] <- "mr_raps"
  if (exists("mr_penalised_weighted_median", where = asNamespace("TwoSampleMR"), inherits = FALSE)) extra[["Penalised weighted median (robust)"]] <- "mr_penalised_weighted_median"
  extra
}

## Picks exactly one, pre-specified row per CpG from a multi-method mr() results data
## frame - IVW for CpGs with >=2 instruments, Wald ratio for single-instrument CpGs -
## matching the tiering rule already documented in the UI ("1 instrument -> Wald ratio
## only, 2 -> IVW only, >=3 -> your selection above"). This selection never looks at
## p-values: choosing the row with the smallest p-value across several methods run for
## the same CpG, then correcting only across CpGs, understates the true number of tests
## performed (a method is silently "selected" by its own p-value) and inflates the
## apparent number of significant CpGs. Every other method run for a CpG stays visible
## in the full Results table as a disclosed sensitivity check; it just isn't the one
## used for headline significance or for the Manhattan plot.
mmr_primary_row_per_cpg <- function(results_df) {
  ivw_name <- {
    ml <- tryCatch(TwoSampleMR::mr_method_list(), error = function(e) NULL)
    nm <- if (!is.null(ml)) ml$name[ml$obj == "mr_ivw"] else character(0)
    if (length(nm)) nm[1] else "Inverse variance weighted"
  }
  do.call(rbind, lapply(split(results_df, results_df$exposure), function(x) {
    primary <- x[x$method == ivw_name, , drop = FALSE]
    if (nrow(primary) == 0) primary <- x[x$method == "Wald ratio", , drop = FALSE]
    if (nrow(primary) == 0) primary <- x[1, , drop = FALSE]
    primary[1, , drop = FALSE]
  }))
}

mmr_analysis_controls <- function(ns) {
  extra <- mmr_robust_method_choices()
  base_names <- list("Inverse variance weighted (IVW)", "MR-Egger", "Weighted median", "Weighted mode", "Simple mode")
  base_values <- list("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_simple_mode")
  tagList(
    checkboxGroupInput(
      ns("mr_methods"), "MR methods",
      choiceNames = c(base_names, names(extra)),
      choiceValues = c(base_values, unname(extra)),
      selected = c(unlist(base_values), unname(extra))
    ),
    p(class = "submodule-desc", "Applied per CpG per this project's own tiered rule: 1 instrument -> Wald ratio only, 2 -> IVW only, >=3 -> your selection above (intersected with what's estimable)."),
    p(class = "empty-note", style = "font-size:12px;", icon("circle-info"),
      "Contamination mixture and Radial MR are not wired into method selection in this deployment (contamination mixture isn't provided by the installed TwoSampleMR version; Radial MR uses a different interface than the other methods here)."),
    div(style = "display:flex; align-items:center; gap:4px;",
        checkboxInput(ns("mr_run_presso"), "Also run MR-PRESSO outlier test (>=4 instruments)", value = FALSE),
        .mmr_tip("Simulation-based horizontal-pleiotropy/outlier test (MRPRESSO::mr_presso). Off by default - adds runtime per CpG.")),
    div(style = "display:flex; align-items:center; gap:4px;",
        numericInput(ns("mr_ci_level"), "Confidence interval level", value = 0.95, min = 0.5, max = 0.999, step = 0.01, width = "100%")),
    div(style = "margin-top: 6px;",
        actionButton(ns("mr_run_btn"), "Run MR", icon = icon("play"), class = "btn-primary btn-sm"))
  )
}

## ---- Tab 6: Sensitivity -------------------------------------------------------

mmr_sensitivity_ui <- function(ns) {
  uiOutput(ns("sensitivity_tab_body"))
}

## ---- Tab 7: Results ------------------------------------------------------------

mmr_results_ui <- function(ns) {
  uiOutput(ns("results_tab_body"))
}

## ---- Tab 8: Plots ---------------------------------------------------------------

mmr_plots_ui <- function(ns) {
  uiOutput(ns("plots_tab_body"))
}

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------

mod_methyl_mr_server <- function(id, methyl_dataset, methyl_results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## ------------------------------------------------------------------
    ## Stage has-run flags + generic downstream invalidation
    ## ------------------------------------------------------------------
    stage_flags <- reactiveValues(data = FALSE, instruments = FALSE, clump = FALSE, harmonise = FALSE, mr = FALSE, sensitivity = FALSE)
    invalidate_from <- function(stage) {
      idx <- match(stage, .mmr_stage_order)
      for (s in .mmr_stage_order[idx:length(.mmr_stage_order)]) stage_flags[[s]] <- FALSE
    }

    ## CpG panel picker for the Preloaded route, from the cached
    ## mr_estimates_{sex}.csv CpG list (script08's majority-vote panel).
    pre_panel_cpgs <- reactive({
      switch(input$pre_sex %||% "female",
        female = { d <- load_default_mr_estimates("female"); if (is.null(d)) character(0) else sort(unique(d$exposure)) },
        male   = { d <- load_default_mr_estimates("male");   if (is.null(d)) character(0) else sort(unique(d$exposure)) },
        both   = {
          f <- load_default_mr_estimates("female"); m <- load_default_mr_estimates("male")
          sort(unique(c(if (!is.null(f)) f$exposure else character(0), if (!is.null(m)) m$exposure else character(0))))
        }
      )
    })

    output$pre_cpg_ui <- renderUI({
      req(METH_DATA_AVAILABLE)
      cpgs <- pre_panel_cpgs()
      validate(need(length(cpgs) > 0, "No cached MR results found for this stratum."))
      tagList(
        selectInput(ns("pre_cpgs"), "CpG panel", choices = cpgs, selected = cpgs, multiple = TRUE, width = "100%"),
        p(class = "submodule-desc", style = "font-size:12px; margin-bottom:0;",
          sprintf("%d CpG(s) - the script07 majority-vote panel (n_votes >= 2) for this stratum.", length(cpgs)))
      )
    })

    ## Any Data-tab defining input change invalidates everything downstream.
    observeEvent(list(input$data_source, input$pre_sex, input$pre_cpgs, input$pre_binary_outcome,
                       input$exp_file, input$out_file, input$outcome_label, input$up_binary_outcome),
                 invalidate_from("instruments"), ignoreInit = TRUE, ignoreNULL = FALSE)

    ## ------------------------------------------------------------------
    ## Upload route: shared column-mapping utilities (global.R; also used
    ## by mod_mr.R and mod_coloc.R's upload modes)
    ## ------------------------------------------------------------------
    exp_df_r <- reactive({ req(input$exp_file); read_uploaded_table(input$exp_file$datapath) })
    out_df_r <- reactive({ req(input$out_file); read_uploaded_table(input$out_file$datapath) })

    output$exp_map_ui <- gwas_col_map_ui(ns, reactive(input$exp_file), exp_df_r, "exp", "Exposure file", extra_fields = "n")
    output$out_map_ui <- gwas_col_map_ui(ns, reactive(input$out_file), out_df_r, "out", "Outcome file", extra_fields = "n")

    ## CpG ID + optional chr/pos/gene columns - mQTL-specific, outside the
    ## shared gwas_col_map_ui (generic GWAS fields only).
    output$exp_extra_map_ui <- renderUI({
      req(input$exp_file)
      df <- exp_df_r()
      validate(need(!is.null(df), "Could not read the exposure file."))
      cols <- colnames(df)
      guess <- function(patterns) { g <- guess_gwas_col(cols, patterns); if (is.na(g)) NA_character_ else g }
      cpg_guess <- guess(c("^cpg", "^probe", "^cg_?id"))
      tagList(
        selectInput(ns("exp_cpg"), "CpG ID column", cols, selected = if (!is.na(cpg_guess)) cpg_guess else cols[1]),
        p(class = "submodule-desc", style = "font-size:12px;", "Optional - only needed if you want cis/trans restriction or a CpG-level Manhattan overview:"),
        fluidRow(
          column(3, selectInput(ns("exp_snp_chr"), "SNP chr", c("(none)" = "", cols), selected = { g <- guess(c("^chr", "^chromosome")); if (is.na(g)) "" else g })),
          column(3, selectInput(ns("exp_snp_pos"), "SNP position", c("(none)" = "", cols), selected = { g <- guess(c("^pos", "^bp$", "base_pair")); if (is.na(g)) "" else g })),
          column(3, selectInput(ns("exp_cpg_chr"), "CpG chr", c("(none)" = "", cols), selected = "")),
          column(3, selectInput(ns("exp_cpg_pos"), "CpG position", c("(none)" = "", cols), selected = ""))
        ),
        selectInput(ns("exp_gene"), "Gene (optional annotation)", c("(none)" = "", cols), selected = { g <- guess(c("^gene")); if (is.na(g)) "" else g })
      )
    })

    ## ------------------------------------------------------------------
    ## Stage 1: Data -> validated raw inputs + Data Validation panel
    ## ------------------------------------------------------------------
    data_state <- reactiveVal(NULL)

    build_data_state <- function() {
      if (identical(input$data_source, "preloaded")) {
        req(METH_DATA_AVAILABLE)
        cpgs <- input$pre_cpgs
        validate(need(length(cpgs) > 0, "Pick at least one CpG."))
        harm <- load_default_mr_harmonised()
        validate(need(!is.null(harm), "Preloaded harmonised MR dataset not found."))
        harm <- harm[harm$exposure %in% cpgs & harm$mr_keep %in% TRUE, , drop = FALSE]
        list(mode = "preloaded", cpgs = cpgs, harmonised = harm,
             binary_outcome = isTRUE(input$pre_binary_outcome),
             outcome_label = if (nrow(harm) > 0) harm$outcome[1] else "Rheumatoid arthritis")
      } else {
        req(input$exp_file, input$out_file, input$exp_cpg,
            input$exp_snp, input$exp_beta, input$exp_se, input$exp_pval, input$exp_ea, input$exp_oa,
            input$out_snp, input$out_beta, input$out_se, input$out_pval, input$out_ea, input$out_oa)
        exp_raw <- exp_df_r(); out_raw <- out_df_r()
        validate(need(!is.null(exp_raw) && !is.null(out_raw), "Could not read one of the uploaded files."))
        list(mode = "upload", exp_raw = exp_raw, out_raw = out_raw,
             binary_outcome = isTRUE(input$up_binary_outcome),
             outcome_label = if (nzchar(trimws(input$outcome_label %||% ""))) trimws(input$outcome_label) else "Uploaded outcome")
      }
    }

    observeEvent(input$validate_btn, {
      st <- tryCatch(build_data_state(), error = function(e) e)
      if (inherits(st, "error") || inherits(st, "shiny.silent.error")) return()
      data_state(st)
      stage_flags$data <- TRUE
      invalidate_from("instruments")
    })

    output$data_validation_ui <- renderUI({
      req(stage_flags$data, data_state())
      st <- data_state()
      if (identical(st$mode, "preloaded")) {
        h <- st$harmonised
        box(
          width = NULL, title = "Data validation", status = "primary", solidHeader = FALSE,
          div(class = "methyl-stats-row", fluidRow(
            valueBox(length(st$cpgs), "CpGs selected", icon = icon("dna"), color = "light-blue", width = 3),
            valueBox(length(unique(h$SNP)), "Unique instrument SNPs", icon = icon("code-branch"), color = "purple", width = 3),
            valueBox(sum(h$mr_keep %in% TRUE), "Harmonised instrument rows", icon = icon("link"), color = "green", width = 3),
            valueBox("GRCh37", "Genome build", icon = icon("dna"), color = "yellow", width = 3)
          )),
          p(class = "empty-note", icon("circle-info"),
            sprintf("Outcome: %s. EAF available for exposure: %s; outcome: %s. Duplicate SNP-CpG rows: %d.",
                    st$outcome_label,
                    sprintf("%.0f%%", 100 * mean(!is.na(h$eaf.exposure))),
                    sprintf("%.0f%%", 100 * mean(!is.na(h$eaf.outcome))),
                    sum(duplicated(h[, c("SNP", "exposure")]))))
        )
      } else {
        exp_raw <- st$exp_raw; out_raw <- st$out_raw
        cpg_col <- input$exp_cpg
        req(cpg_col %in% colnames(exp_raw))
        n_cpg <- length(unique(exp_raw[[cpg_col]]))
        n_snp_exp <- length(unique(exp_raw[[input$exp_snp]]))
        n_snp_out <- length(unique(out_raw[[input$out_snp]]))
        overlap <- length(intersect(exp_raw[[input$exp_snp]], out_raw[[input$out_snp]]))
        miss_key <- function(df, cols) sum(!stats::complete.cases(df[, cols[cols %in% colnames(df)], drop = FALSE]))
        invalid_alleles <- sum(!grepl("^[ACGTDI]+$", toupper(as.character(exp_raw[[input$exp_ea]]))) |
                                  !grepl("^[ACGTDI]+$", toupper(as.character(exp_raw[[input$exp_oa]]))))
        dup_snp <- sum(duplicated(exp_raw[[input$exp_snp]]))
        eaf_avail_exp <- if (nzchar(input$exp_eaf %||% "")) sprintf("%.0f%%", 100 * mean(!is.na(exp_raw[[input$exp_eaf]]))) else "no EAF column mapped"
        eaf_avail_out <- if (nzchar(input$out_eaf %||% "")) sprintf("%.0f%%", 100 * mean(!is.na(out_raw[[input$out_eaf]]))) else "no EAF column mapped"
        chr_pos_avail <- nzchar(input$exp_snp_chr %||% "") && nzchar(input$exp_snp_pos %||% "")
        pos_range <- if (chr_pos_avail) suppressWarnings(range(as.numeric(exp_raw[[input$exp_snp_pos]]), na.rm = TRUE)) else c(NA, NA)
        build_guess <- if (chr_pos_avail && is.finite(pos_range[2])) {
          if (pos_range[2] < 3e8) "Not distinguishable from position alone (GRCh37/GRCh38 ranges overlap) - preloaded pipeline uses GRCh37." else "Could not estimate."
        } else "Not determinable (no chr/position mapped)."

        box(
          width = NULL, title = "Data validation", status = "primary", solidHeader = FALSE,
          div(class = "methyl-stats-row", fluidRow(
            valueBox(n_cpg, "CpGs (exposure file)", icon = icon("dna"), color = "light-blue", width = 3),
            valueBox(format(n_snp_exp, big.mark = ","), "Unique SNPs (exposure)", icon = icon("code-branch"), color = "purple", width = 3),
            valueBox(format(n_snp_out, big.mark = ","), "Unique SNPs (outcome)", icon = icon("code-branch"), color = "aqua", width = 3),
            valueBox(format(overlap, big.mark = ","), "Overlapping SNPs", icon = icon("circle-check"), color = "green", width = 3)
          )),
          tags$table(class = "table table-sm",
            tags$tbody(
              tags$tr(tags$td("Missing beta/SE/p-value (exposure)"), tags$td(miss_key(exp_raw, c(input$exp_beta, input$exp_se, input$exp_pval)))),
              tags$tr(tags$td("Missing beta/SE/p-value (outcome)"), tags$td(miss_key(out_raw, c(input$out_beta, input$out_se, input$out_pval)))),
              tags$tr(tags$td("Duplicate SNP rows (exposure)"), tags$td(dup_snp)),
              tags$tr(tags$td("Invalid allele codes (exposure)"), tags$td(invalid_alleles)),
              tags$tr(tags$td("EAF available - exposure / outcome"), tags$td(sprintf("%s / %s", eaf_avail_exp, eaf_avail_out))),
              tags$tr(tags$td("Chromosome/position available (exposure SNPs)"), tags$td(if (chr_pos_avail) "Yes" else "No")),
              tags$tr(tags$td("Exposure / outcome sample size columns mapped"), tags$td(sprintf("%s / %s", if (nzchar(input$exp_n %||% "")) "Yes" else "No", if (nzchar(input$out_n %||% "")) "Yes" else "No"))),
              tags$tr(tags$td("Detected genome build"), tags$td(build_guess))
            )
          )
        )
      }
    })
    outputOptions(output, "data_validation_ui", suspendWhenHidden = FALSE)

    ## ------------------------------------------------------------------
    ## Stage 2: Filters & Instruments
    ## ------------------------------------------------------------------
    instruments_state <- reactiveVal(NULL)

    build_instruments_state <- function() {
      st <- data_state(); req(st)
      pval_cut <- input$f_pval %||% MMR_DEFAULT_CIS_PVAL
      min_f <- input$f_min_f %||% MMR_DEFAULT_MIN_F
      maf_cut <- input$f_maf %||% 0
      region_mode <- input$f_region_mode %||% "cis"
      cis_window <- (input$f_cis_window %||% (MMR_DEFAULT_CIS_WINDOW_BP / 1000)) * 1000

      if (identical(st$mode, "preloaded")) {
        h <- st$harmonised
        h$F_stat_recomputed <- (h$beta.exposure / h$se.exposure)^2
        d <- h
        excl_pval <- d$pval.exposure > pval_cut
        excl_maf <- rep(FALSE, nrow(d))
        if (isTRUE(input$f_require_eaf)) excl_maf <- excl_maf | is.na(d$eaf.exposure)
        maf_val <- pmin(d$eaf.exposure, 1 - d$eaf.exposure)
        excl_maf <- excl_maf | (!is.na(maf_val) & maf_val < maf_cut)
        weak <- d$F_stat_recomputed < min_f
        ## Preloaded pipeline already restricted to cis (+/-1Mb) at build
        ## time - this control is informational only for this route.
        excl_region <- rep(FALSE, nrow(d))

        d$excluded_pval <- excl_pval
        d$excluded_region <- excl_region
        d$excluded_maf <- excl_maf
        d$weak_instrument <- weak
        d$retained <- !excl_pval & !excl_region & !excl_maf & (!weak | !isTRUE(input$f_exclude_weak))

        counts_before <- table(d$exposure[d$retained | d$weak_instrument | d$excluded_pval | d$excluded_maf])
        list(mode = "preloaded", binary_outcome = st$binary_outcome, outcome_label = st$outcome_label,
             region_mode_note = "Already restricted to cis (+/-1Mb) when this cached instrument set was built.",
             d = d)
      } else {
        exp_raw <- st$exp_raw
        cpg_col <- input$exp_cpg
        fmt_args <- list(
          exp_raw, type = "exposure", snp_col = input$exp_snp, beta_col = input$exp_beta,
          se_col = input$exp_se, pval_col = input$exp_pval, effect_allele_col = input$exp_ea,
          other_allele_col = input$exp_oa, phenotype_col = cpg_col, id_col = cpg_col
        )
        if (nzchar(input$exp_eaf %||% "")) fmt_args$eaf_col <- input$exp_eaf
        if (nzchar(input$exp_n %||% "")) fmt_args$samplesize_col <- input$exp_n
        if (nzchar(input$exp_snp_chr %||% "")) fmt_args$chr_col <- input$exp_snp_chr
        if (nzchar(input$exp_snp_pos %||% "")) fmt_args$pos_col <- input$exp_snp_pos
        exp_fmt <- do.call(TwoSampleMR::format_data, fmt_args)
        validate(need(nrow(exp_fmt) > 0, "No usable rows after formatting the exposure file - check column mapping."))

        exp_fmt$F_stat <- (exp_fmt$beta.exposure / exp_fmt$se.exposure)^2
        excl_pval <- exp_fmt$pval.exposure > pval_cut
        maf_val <- if (!is.null(exp_fmt$eaf.exposure)) pmin(exp_fmt$eaf.exposure, 1 - exp_fmt$eaf.exposure) else rep(NA_real_, nrow(exp_fmt))
        excl_maf <- (isTRUE(input$f_require_eaf) & is.na(exp_fmt$eaf.exposure)) | (!is.na(maf_val) & maf_val < maf_cut)
        weak <- exp_fmt$F_stat < min_f

        has_region_cols <- all(c("chr.exposure", "pos.exposure") %in% colnames(exp_fmt)) &&
          nzchar(input$exp_cpg_chr %||% "") && nzchar(input$exp_cpg_pos %||% "")
        excl_region <- rep(FALSE, nrow(exp_fmt))
        region_note <- NULL
        if (identical(region_mode, "cis")) {
          if (has_region_cols) {
            cpg_chr <- as.character(exp_raw[[input$exp_cpg_chr]])
            cpg_pos <- suppressWarnings(as.numeric(exp_raw[[input$exp_cpg_pos]]))
            ## Rejoin by SNP+phenotype since format_data() may reorder rows.
            key <- match(paste(exp_fmt$SNP, exp_fmt$exposure), paste(exp_raw[[input$exp_snp]], exp_raw[[cpg_col]]))
            row_cpg_chr <- cpg_chr[key]; row_cpg_pos <- cpg_pos[key]
            same_chr <- as.character(exp_fmt$chr.exposure) == row_cpg_chr
            within <- abs(exp_fmt$pos.exposure - row_cpg_pos) <= cis_window
            excl_region <- !(same_chr & within & !is.na(same_chr) & !is.na(within))
          } else {
            region_note <- "Cis/trans restriction requested but SNP and/or CpG chromosome/position columns weren't mapped - region filter not applied."
          }
        }

        exp_fmt$excluded_pval <- excl_pval
        exp_fmt$excluded_region <- excl_region
        exp_fmt$excluded_maf <- excl_maf
        exp_fmt$weak_instrument <- weak
        exp_fmt$retained <- !excl_pval & !excl_region & !excl_maf & (!weak | !isTRUE(input$f_exclude_weak))

        list(mode = "upload", binary_outcome = st$binary_outcome, outcome_label = st$outcome_label,
             out_raw = st$out_raw, region_mode_note = region_note, d = exp_fmt)
      }
    }

    observeEvent(input$instruments_btn, {
      ist <- tryCatch(build_instruments_state(), error = function(e) e)
      if (inherits(ist, "error") || inherits(ist, "shiny.silent.error")) return()
      ## Enforce min/max instruments per CpG (retained set), keeping the
      ## strongest by p-value; excluded rows are flagged (excluded_maxcap),
      ## not dropped.
      d <- ist$d
      d$excluded_maxcap <- FALSE
      max_n <- input$f_max_instruments
      if (!is.null(max_n) && !is.na(max_n) && is.finite(max_n)) {
        for (cpg_id in unique(d$exposure[d$retained])) {
          idx <- which(d$exposure == cpg_id & d$retained)
          if (length(idx) > max_n) {
            keep <- idx[order(d$pval.exposure[idx])][seq_len(max_n)]
            drop <- setdiff(idx, keep)
            d$excluded_maxcap[drop] <- TRUE
            d$retained[drop] <- FALSE
          }
        }
      }
      min_n <- input$f_min_instruments %||% 1
      n_per_cpg <- table(d$exposure[d$retained])
      too_few_cpgs <- names(n_per_cpg)[n_per_cpg < min_n]
      if (length(too_few_cpgs) > 0) d$retained[d$exposure %in% too_few_cpgs] <- FALSE
      ist$d <- d
      instruments_state(ist)
      stage_flags$instruments <- TRUE
      invalidate_from("clump")
    })
    observeEvent(list(input$f_pval, input$f_region_mode, input$f_cis_window, input$f_maf, input$f_require_eaf,
                       input$f_min_instruments, input$f_max_instruments, input$f_min_f, input$f_exclude_weak),
                 invalidate_from("clump"), ignoreInit = TRUE)

    output$filters_tab_body <- renderUI({
      req(stage_flags$data)
      tagList(
        fluidRow(column(6, mmr_filters_controls(ns)), column(6, uiOutput(ns("instruments_summary_ui"))))
      )
    })
    outputOptions(output, "filters_tab_body", suspendWhenHidden = FALSE)

    output$instruments_summary_ui <- renderUI({
      req(stage_flags$instruments, instruments_state())
      ist <- instruments_state()
      d <- ist$d
      n_retained <- sum(d$retained)
      n_weak <- sum(d$weak_instrument & !d$excluded_pval & !d$excluded_region & !d$excluded_maf)
      n_missing_invalid <- sum(d$excluded_pval | d$excluded_region | d$excluded_maf)
      tagList(
        box(
          width = NULL, title = "Instrument selection summary", status = "primary", solidHeader = FALSE,
          div(class = "methyl-stats-row", fluidRow(
            valueBox(n_retained, "Retained instruments", icon = icon("check"), color = "green", width = 4),
            valueBox(n_weak, "Excluded: weak (low F)", icon = icon("triangle-exclamation"), color = "yellow", width = 4),
            valueBox(n_missing_invalid, "Excluded: filtered", icon = icon("xmark"), color = "red", width = 4)
          )),
          if (!is.null(ist$region_mode_note)) p(class = "empty-note", icon("circle-info"), ist$region_mode_note),
          p(class = "empty-note", style = "font-size:12px;", sprintf("%d CpG(s) have at least one retained instrument.", length(unique(d$exposure[d$retained]))))
        )
      )
    })
    outputOptions(output, "instruments_summary_ui", suspendWhenHidden = FALSE)

    ## ------------------------------------------------------------------
    ## Stage 3: LD Clumping
    ## ------------------------------------------------------------------
    clump_state <- reactiveVal(NULL)

    observeEvent(input$clump_btn, {
      ist <- instruments_state(); req(ist)
      d <- ist$d[ist$d$retained, , drop = FALSE]
      validate(need(nrow(d) > 0, "No retained instruments to clump."))

      if (identical(ist$mode, "preloaded")) {
        counts <- load_default_mr_instrument_counts()
        clump_state(list(mode = "preloaded", d = d, counts = counts,
                          n_before = NA_integer_, n_after = if (!is.null(counts)) sum(counts$n_snp[counts$cpg %in% unique(d$exposure)]) else nrow(d),
                          note = "Preloaded route: already LD-clumped when this cached instrument set was built (r2<0.001, 10,000kb, GoDMC's own reference) - not re-run live. Counts below are the recorded per-CpG instrument counts after that clumping."))
      } else {
        n_before <- nrow(d)
        clumped <- tryCatch(
          ieugwasr::ld_clump(
            dplyr::tibble(rsid = d$SNP, pval = d$pval.exposure, id = d$exposure),
            clump_kb = input$clump_kb %||% MMR_DEFAULT_CLUMP_KB,
            clump_r2 = input$clump_r2 %||% MMR_DEFAULT_CLUMP_R2,
            clump_p = input$clump_p %||% 1,
            pop = input$clump_pop %||% "EUR"
          ),
          error = function(e) NULL
        )
        if (is.null(clumped)) {
          clump_state(list(mode = "upload", d = d, n_before = n_before, n_after = n_before, status = "api_error",
                            note = "Could not reach the OpenGWAS LD reference (network/token issue) - proceeding on the unclumped instrument set."))
        } else if (nrow(clumped) == 0) {
          clump_state(list(mode = "upload", d = d, n_before = n_before, n_after = n_before, status = "no_match",
                            note = "Clumping ran but matched none of these SNP IDs in the chosen population's LD reference - proceeding on the unclumped instrument set."))
        } else {
          d2 <- d[d$SNP %in% clumped$rsid, , drop = FALSE]
          removed <- d[!(d$SNP %in% clumped$rsid), , drop = FALSE]
          clump_state(list(mode = "upload", d = d2, n_before = n_before, n_after = nrow(d2), status = "applied",
                            removed = removed,
                            note = sprintf("Removed %d SNP(s) as LD-correlated (r^2 >= %.3g within %g kb) with a retained index SNP.",
                                            nrow(removed), input$clump_r2 %||% MMR_DEFAULT_CLUMP_R2, input$clump_kb %||% MMR_DEFAULT_CLUMP_KB)))
        }
      }
      stage_flags$clump <- TRUE
      invalidate_from("harmonise")
    })
    observeEvent(list(input$clump_pop, input$clump_r2, input$clump_kb, input$clump_p), invalidate_from("harmonise"), ignoreInit = TRUE)

    output$clump_tab_body <- renderUI({
      req(stage_flags$instruments)
      st_mode <- (data_state())$mode
      tagList(
        fluidRow(
          column(6, mmr_clump_controls(ns)),
          column(6, uiOutput(ns("clump_summary_ui")))
        )
      )
    })
    outputOptions(output, "clump_tab_body", suspendWhenHidden = FALSE)

    output$clump_summary_ui <- renderUI({
      req(stage_flags$clump, clump_state())
      cs <- clump_state()
      box(
        width = NULL, title = "Clumping summary", status = "primary", solidHeader = FALSE,
        div(class = "methyl-stats-row", fluidRow(
          valueBox(if (is.na(cs$n_before)) "n/a" else cs$n_before, "Before clumping", icon = icon("list"), color = "light-blue", width = 6),
          valueBox(cs$n_after, "After clumping", icon = icon("check"), color = "green", width = 6)
        )),
        p(class = "empty-note", icon("circle-info"), cs$note)
      )
    })
    outputOptions(output, "clump_summary_ui", suspendWhenHidden = FALSE)

    ## ------------------------------------------------------------------
    ## Stage 4: Harmonisation
    ## ------------------------------------------------------------------
    harmonise_state <- reactiveVal(NULL)

    observeEvent(input$harmonise_btn, {
      cs <- clump_state(); req(cs)
      if (identical(cs$mode, "preloaded")) {
        d <- cs$d
        summary_counts <- list(
          before = nrow(d), aligned = sum(d$mr_keep %in% TRUE), flipped = NA_integer_,
          palindromic = sum(d$palindromic %in% TRUE), ambiguous = sum(d$ambiguous %in% TRUE),
          mismatches = sum(d$remove %in% TRUE, na.rm = TRUE), removed = sum(!(d$mr_keep %in% TRUE)),
          final = sum(d$mr_keep %in% TRUE)
        )
        harmonise_state(list(mode = "preloaded", harmonised = d[d$mr_keep %in% TRUE, , drop = FALSE], summary = summary_counts,
                              binary_outcome = (data_state())$binary_outcome, outcome_label = (data_state())$outcome_label))
      } else {
        st <- data_state()
        out_raw <- st$out_raw
        fmt_args <- list(
          out_raw, type = "outcome", snp_col = input$out_snp, beta_col = input$out_beta,
          se_col = input$out_se, pval_col = input$out_pval, effect_allele_col = input$out_ea, other_allele_col = input$out_oa
        )
        if (nzchar(input$out_eaf %||% "")) fmt_args$eaf_col <- input$out_eaf
        if (nzchar(input$out_n %||% "")) fmt_args$samplesize_col <- input$out_n
        out_fmt <- do.call(TwoSampleMR::format_data, fmt_args)
        validate(need(nrow(out_fmt) > 0, "No usable rows after formatting the outcome file - check column mapping."))
        out_fmt$outcome <- st$outcome_label

        harmonised <- tryCatch(TwoSampleMR::harmonise_data(cs$d, out_fmt, action = as.integer(input$harmonise_action %||% "2")), error = function(e) NULL)
        validate(need(!is.null(harmonised) && nrow(harmonised) > 0,
          "Harmonisation found no overlapping SNPs between exposure and outcome - check that both use the same SNP identifiers (rsIDs) and that allele columns are mapped correctly."))

        summary_counts <- list(
          before = nrow(harmonised), aligned = sum(harmonised$mr_keep %in% TRUE),
          flipped = sum(harmonised$effect_allele.exposure != harmonised$effect_allele.outcome & harmonised$mr_keep %in% TRUE, na.rm = TRUE),
          palindromic = sum(harmonised$palindromic %in% TRUE), ambiguous = sum(harmonised$ambiguous %in% TRUE),
          mismatches = sum(harmonised$remove %in% TRUE, na.rm = TRUE), removed = sum(!(harmonised$mr_keep %in% TRUE)),
          final = sum(harmonised$mr_keep %in% TRUE)
        )
        harmonised <- harmonised[harmonised$mr_keep %in% TRUE, , drop = FALSE]
        validate(need(nrow(harmonised) > 0, "No instruments survived harmonisation (all ambiguous/unresolvable, e.g. unresolved palindromic SNPs)."))
        harmonise_state(list(mode = "upload", harmonised = harmonised, summary = summary_counts,
                              binary_outcome = st$binary_outcome, outcome_label = st$outcome_label))
      }
      stage_flags$harmonise <- TRUE
      invalidate_from("mr")
    })
    observeEvent(input$harmonise_action, invalidate_from("mr"), ignoreInit = TRUE)

    output$harmonise_tab_body <- renderUI({
      req(stage_flags$clump)
      editable <- identical((data_state())$mode, "upload")
      fluidRow(
        column(6, mmr_harmonise_controls(ns, editable = editable)),
        column(6, uiOutput(ns("harmonise_summary_ui")))
      )
    })
    outputOptions(output, "harmonise_tab_body", suspendWhenHidden = FALSE)

    output$harmonise_summary_ui <- renderUI({
      req(stage_flags$harmonise, harmonise_state())
      s <- harmonise_state()$summary
      box(
        width = NULL, title = "Harmonisation QC summary", status = "primary", solidHeader = FALSE,
        tags$table(class = "table table-sm", tags$tbody(
          tags$tr(tags$td("SNPs before harmonisation"), tags$td(s$before)),
          tags$tr(tags$td("Aligned"), tags$td(s$aligned)),
          tags$tr(tags$td("Flipped (allele swap)"), tags$td(if (is.na(s$flipped)) "not recorded (preloaded)" else s$flipped)),
          tags$tr(tags$td("Palindromic"), tags$td(s$palindromic)),
          tags$tr(tags$td("Ambiguous"), tags$td(s$ambiguous)),
          tags$tr(tags$td("Allele mismatches"), tags$td(s$mismatches)),
          tags$tr(tags$td("Removed"), tags$td(s$removed)),
          tags$tr(tags$td(strong("Final instruments")), tags$td(strong(s$final)))
        ))
      )
    })
    outputOptions(output, "harmonise_summary_ui", suspendWhenHidden = FALSE)

    ## ------------------------------------------------------------------
    ## Stage 5: MR Analysis
    ## ------------------------------------------------------------------
    mr_state <- reactiveVal(NULL)

    tier_methods <- function(n_snp, selected) {
      if (n_snp == 1) return("mr_wald_ratio")
      if (n_snp == 2) return("mr_ivw")
      estimable <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_simple_mode",
                      unname(unlist(mmr_robust_method_choices())))
      intersect(estimable, selected)
    }

    observeEvent(input$mr_run_btn, {
      hs <- harmonise_state(); req(hs)
      dat <- hs$harmonised
      selected <- input$mr_methods %||% c("mr_ivw")
      cpgs <- unique(dat$exposure)
      validate(need(length(cpgs) > 0, "No harmonised instruments to run MR on."))

      results <- withProgress(message = "Running MR", value = 0, {
        out <- lapply(seq_along(cpgs), function(i) {
          incProgress(1 / length(cpgs), detail = cpgs[i])
          sub <- dat[dat$exposure == cpgs[i], , drop = FALSE]
          n_snp <- nrow(sub)
          methods <- tier_methods(n_snp, selected)
          if (length(methods) == 0) methods <- if (n_snp == 1) "mr_wald_ratio" else "mr_ivw"
          res <- tryCatch(TwoSampleMR::mr(sub, method_list = methods), error = function(e) NULL)
          if (is.null(res)) return(NULL)
          res$n_snp_tier <- n_snp
          res
        })
        do.call(rbind, out[!vapply(out, is.null, logical(1))])
      })
      validate(need(!is.null(results) && nrow(results) > 0, "MR produced no estimable results for the current selection."))

      alpha <- 1 - (input$mr_ci_level %||% 0.95)
      z <- stats::qnorm(1 - alpha / 2)
      results$ci_low <- results$b - z * results$se
      results$ci_high <- results$b + z * results$se
      if (isTRUE(hs$binary_outcome)) {
        results$or <- exp(results$b); results$or_low <- exp(results$ci_low); results$or_high <- exp(results$ci_high)
      }

      presso <- NULL
      if (isTRUE(input$mr_run_presso) && requireNamespace("MRPRESSO", quietly = TRUE)) {
        presso <- lapply(cpgs, function(cpg_id) {
          sub <- dat[dat$exposure == cpg_id, , drop = FALSE]
          if (nrow(sub) < 4) return(NULL)
          tryCatch({
            pr <- MRPRESSO::mr_presso(
              BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
              SdOutcome = "se.outcome", SdExposure = "se.exposure", data = as.data.frame(sub),
              OUTLIERtest = TRUE, DISTORTIONtest = TRUE, NbDistribution = 1000, SignifThreshold = 0.05, seed = 2024
            )
            list(cpg = cpg_id, global_p = pr$`MR-PRESSO results`$`Global Test`$Pvalue)
          }, error = function(e) NULL)
        })
        names(presso) <- cpgs
      }

      mr_state(list(cpgs = cpgs, results = results, dat = dat, binary_outcome = isTRUE(hs$binary_outcome),
                     outcome_label = hs$outcome_label, ci_level = input$mr_ci_level %||% 0.95, presso = presso))
      stage_flags$mr <- TRUE
      invalidate_from("sensitivity")
    })
    observeEvent(list(input$mr_methods, input$mr_ci_level, input$mr_run_presso), invalidate_from("sensitivity"), ignoreInit = TRUE)

    output$analysis_tab_body <- renderUI({
      req(stage_flags$harmonise)
      fluidRow(column(6, mmr_analysis_controls(ns)), column(6, uiOutput(ns("mr_summary_ui"))))
    })
    outputOptions(output, "analysis_tab_body", suspendWhenHidden = FALSE)

    output$mr_summary_ui <- renderUI({
      req(stage_flags$mr, mr_state())
      ms <- mr_state()
      box(width = NULL, title = "MR run summary", status = "primary", solidHeader = FALSE,
          p(sprintf("%d CpG(s) tested, %d method-level estimate(s) produced. See Results and Plots tabs.", length(ms$cpgs), nrow(ms$results))))
    })
    outputOptions(output, "mr_summary_ui", suspendWhenHidden = FALSE)

    ## ------------------------------------------------------------------
    ## Stage 6: Sensitivity
    ## ------------------------------------------------------------------
    sensitivity_state <- reactiveVal(NULL)

    observeEvent(input$sensitivity_run_btn, {
      ms <- mr_state(); req(ms)
      dat <- ms$dat
      n_per_cpg <- table(dat$exposure)
      cpgs_ge3 <- names(n_per_cpg)[n_per_cpg >= 3]

      het <- pleio <- loo <- single_snp <- NULL
      if (length(cpgs_ge3) > 0) {
        sub <- dat[dat$exposure %in% cpgs_ge3, , drop = FALSE]
        het <- tryCatch(TwoSampleMR::mr_heterogeneity(sub, method_list = c("mr_ivw", "mr_egger_regression")), error = function(e) NULL)
        pleio <- tryCatch(TwoSampleMR::mr_pleiotropy_test(sub), error = function(e) NULL)
        loo <- tryCatch(TwoSampleMR::mr_leaveoneout(sub), error = function(e) NULL)
      }
      single_snp <- tryCatch(TwoSampleMR::mr_singlesnp(dat), error = function(e) NULL)

      directionality <- NULL
      has_n <- all(c("samplesize.exposure", "samplesize.outcome") %in% colnames(dat)) &&
        !all(is.na(dat$samplesize.exposure)) && !all(is.na(dat$samplesize.outcome))
      if (has_n) directionality <- tryCatch(TwoSampleMR::directionality_test(dat), error = function(e) NULL)

      f_tab <- data.frame(
        cpg = dat$exposure, SNP = dat$SNP,
        F_stat = (dat$beta.exposure / dat$se.exposure)^2
      )

      sensitivity_state(list(het = het, pleio = pleio, loo = loo, single_snp = single_snp,
                              directionality = directionality, f_tab = f_tab, cpgs_ge3 = cpgs_ge3,
                              n_per_cpg = n_per_cpg))
      stage_flags$sensitivity <- TRUE
    })

    output$sensitivity_tab_body <- renderUI({
      req(stage_flags$mr)
      tagList(
        p(class = "submodule-desc",
          "Heterogeneity (Cochran's Q), pleiotropy (MR-Egger intercept), and leave-one-out are only computable for a CpG with >=3 instruments. Single-SNP estimates and F-statistics are always shown."),
        actionButton(ns("sensitivity_run_btn"), "Run Sensitivity", icon = icon("magnifying-glass-chart"), class = "btn-primary btn-sm"),
        uiOutput(ns("sensitivity_results_ui"))
      )
    })
    outputOptions(output, "sensitivity_tab_body", suspendWhenHidden = FALSE)

    output$sensitivity_results_ui <- renderUI({
      req(stage_flags$sensitivity, sensitivity_state())
      ss <- sensitivity_state()
      tagList(
        br(),
        box(width = 12, title = "Weak-instrument / F-statistic diagnostics", status = "primary", solidHeader = FALSE,
            DT::dataTableOutput(ns("f_stat_table"))),
        box(width = 12, title = "Heterogeneity (Cochran's Q)", status = "primary", solidHeader = FALSE,
            if (is.null(ss$het) || nrow(ss$het) == 0) p(class = "empty-note", icon("circle-info"), "No CpG in this run has >=3 instruments - heterogeneity is not computable.")
            else DT::dataTableOutput(ns("het_table"))),
        box(width = 12, title = "Horizontal pleiotropy (MR-Egger intercept)", status = "primary", solidHeader = FALSE,
            if (is.null(ss$pleio) || nrow(ss$pleio) == 0) p(class = "empty-note", icon("circle-info"), "Not computable (needs >=3 instruments).")
            else DT::dataTableOutput(ns("pleio_table"))),
        box(width = 12, title = "Leave-one-out", status = "primary", solidHeader = FALSE,
            if (is.null(ss$loo) || nrow(ss$loo) == 0) p(class = "empty-note", icon("circle-info"), "Not computable (needs >=3 instruments).")
            else DT::dataTableOutput(ns("loo_table"))),
        box(width = 12, title = "Single-SNP estimates", status = "primary", solidHeader = FALSE,
            if (is.null(ss$single_snp)) p(class = "empty-note", icon("circle-info"), "Not available.") else DT::dataTableOutput(ns("single_snp_table"))),
        box(width = 12, title = "Steiger directionality", status = "primary", solidHeader = FALSE,
            if (is.null(ss$directionality)) p(class = "empty-note", icon("circle-info"), "Needs exposure and outcome sample sizes - not available for this run.")
            else DT::dataTableOutput(ns("directionality_table")))
      )
    })
    outputOptions(output, "sensitivity_results_ui", suspendWhenHidden = FALSE)

    output$f_stat_table <- DT::renderDataTable({
      req(sensitivity_state())
      d <- sensitivity_state()$f_tab
      d$F_stat <- round(d$F_stat, 1)
      d$weak <- d$F_stat < (input$f_min_f %||% MMR_DEFAULT_MIN_F)
      colnames(d) <- c("CpG", "SNP", "F-statistic", "Weak (< threshold)")
      DT::datatable(d, rownames = FALSE, filter = "top", options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })
    outputOptions(output, "f_stat_table", suspendWhenHidden = FALSE)
    output$het_table <- DT::renderDataTable({
      req(sensitivity_state()$het)
      DT::datatable(sensitivity_state()$het, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })
    outputOptions(output, "het_table", suspendWhenHidden = FALSE)
    output$pleio_table <- DT::renderDataTable({
      req(sensitivity_state()$pleio)
      DT::datatable(sensitivity_state()$pleio, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })
    outputOptions(output, "pleio_table", suspendWhenHidden = FALSE)
    output$loo_table <- DT::renderDataTable({
      req(sensitivity_state()$loo)
      DT::datatable(sensitivity_state()$loo, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })
    outputOptions(output, "loo_table", suspendWhenHidden = FALSE)
    output$single_snp_table <- DT::renderDataTable({
      req(sensitivity_state()$single_snp)
      DT::datatable(sensitivity_state()$single_snp, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })
    outputOptions(output, "single_snp_table", suspendWhenHidden = FALSE)
    output$directionality_table <- DT::renderDataTable({
      req(sensitivity_state()$directionality)
      DT::datatable(sensitivity_state()$directionality, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })
    outputOptions(output, "directionality_table", suspendWhenHidden = FALSE)

    ## ------------------------------------------------------------------
    ## CpG annotation (chr/pos/gene) via 450K manifest lookup (annotation.R).
    ## ------------------------------------------------------------------
    cpg_annotation <- function(cpgs) {
      anno <- tryCatch(methyl_get_annotation("450K"), error = function(e) list(ok = FALSE))
      if (isTRUE(anno$ok)) {
        hit <- anno$anno[cpgs, c("chr", "pos", "gene"), drop = FALSE]
        rownames(hit) <- cpgs
        hit
      } else {
        data.frame(chr = NA_character_, pos = NA_integer_, gene = NA_character_, row.names = cpgs)
      }
    }

    ## ------------------------------------------------------------------
    ## Stage 7: Results
    ## ------------------------------------------------------------------
    output$results_tab_body <- renderUI({
      req(stage_flags$mr)
      ms <- mr_state()
      anno <- cpg_annotation(ms$cpgs)
      tagList(
        box(width = 12, title = "Primary MR results", status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "One row per CpG per MR method actually estimable given its instrument count. A single nominal p < 0.05 is not, on its own, evidence of a causal effect - check the Sensitivity tab and the adjusted p-value below."),
            div(class = "table-toolbar", downloadButton(ns("dl_mr_results"), "MR results (CSV)", class = "btn-sm")),
            DT::dataTableOutput(ns("primary_results_table"))),
        box(width = 12, title = "Adjusted results (multiple-testing correction)", status = "primary", solidHeader = FALSE,
            if (length(unique(ms$results$exposure)) < 2)
              p(class = "empty-note", icon("circle-info"), "Only one CpG was tested - multiple-testing correction is not meaningful here.")
            else tagList(
              selectInput(ns("padj_method"), "Correction method", c("Benjamini-Hochberg (FDR)" = "BH", "Bonferroni" = "bonferroni"), selected = "BH", width = "260px"),
              DT::dataTableOutput(ns("adjusted_results_table"))
            )),
        box(width = 12, title = "QC summary", status = "primary", solidHeader = FALSE,
            uiOutput(ns("qc_summary_ui"))),
        div(class = "table-toolbar",
            downloadButton(ns("dl_harmonised"), "Harmonised instruments (CSV)", class = "btn-sm"),
            if (stage_flags$sensitivity) downloadButton(ns("dl_sensitivity"), "Sensitivity results (CSV)", class = "btn-sm"),
            downloadButton(ns("dl_qc_report"), "QC report (TXT)", class = "btn-sm"),
            downloadButton(ns("dl_full_report"), "Complete MR report (TXT)", class = "btn-sm"))
      )
    })
    outputOptions(output, "results_tab_body", suspendWhenHidden = FALSE)

    build_primary_results_df <- function() {
      ms <- mr_state(); req(ms)
      anno <- cpg_annotation(ms$cpgs)
      r <- ms$results
      r$gene <- anno[r$exposure, "gene"]
      r$chr <- anno[r$exposure, "chr"]
      r$pos <- anno[r$exposure, "pos"]
      r
    }

    output$primary_results_table <- DT::renderDataTable({
      r <- build_primary_results_df()
      cols <- c(CpG = "exposure", Gene = "gene", Chr = "chr", Position = "pos", Outcome = "outcome",
                Method = "method", nSNP = "nsnp", beta = "b", SE = "se", p = "pval", `CI low` = "ci_low", `CI high` = "ci_high")
      if ("or" %in% colnames(r)) cols <- c(cols, OR = "or", `OR CI low` = "or_low", `OR CI high` = "or_high")
      cols <- cols[cols %in% colnames(r)]
      d <- r[, cols, drop = FALSE]; colnames(d) <- names(cols)
      for (nm in intersect(c("beta", "SE", "CI low", "CI high", "OR", "OR CI low", "OR CI high"), colnames(d))) d[[nm]] <- round(d[[nm]], 4)
      if ("p" %in% colnames(d)) d[["p"]] <- signif(d[["p"]], 3)
      DT::datatable(d, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    outputOptions(output, "primary_results_table", suspendWhenHidden = FALSE)

    adjusted_results_df <- reactive({
      r <- build_primary_results_df()
      primary_per_cpg <- mmr_primary_row_per_cpg(r)
      method_r <- switch(input$padj_method %||% "BH", BH = "BH", bonferroni = "bonferroni")
      primary_per_cpg$p_adj <- stats::p.adjust(primary_per_cpg$pval, method = method_r)
      primary_per_cpg$significant <- primary_per_cpg$p_adj < 0.05
      primary_per_cpg
    })

    output$adjusted_results_table <- DT::renderDataTable({
      d <- adjusted_results_df()
      out <- data.frame(CpG = d$exposure, Gene = d$gene, Method = d$method, `Raw p` = signif(d$pval, 3),
                         `Adjusted p` = signif(d$p_adj, 3), `N tested` = nrow(d),
                         Significant = ifelse(d$significant, "Yes", "No"), check.names = FALSE)
      DT::datatable(out, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    outputOptions(output, "adjusted_results_table", suspendWhenHidden = FALSE)

    output$qc_summary_ui <- renderUI({
      cs <- clump_state(); hs <- harmonise_state(); ss <- sensitivity_state()
      tags$table(class = "table table-sm", tags$tbody(
        tags$tr(tags$td("Instruments before clumping"), tags$td(if (is.null(cs) || is.na(cs$n_before)) "n/a (preloaded)" else cs$n_before)),
        tags$tr(tags$td("Instruments after clumping"), tags$td(if (is.null(cs)) "n/a" else cs$n_after)),
        tags$tr(tags$td("Harmonisation exclusions"), tags$td(if (is.null(hs)) "n/a" else hs$summary$removed)),
        tags$tr(tags$td("Weak instruments flagged"), tags$td(sum((instruments_state()$d$weak_instrument) %||% FALSE))),
        tags$tr(tags$td("CpGs with heterogeneity/pleiotropy computed (>=3 instruments)"), tags$td(if (is.null(ss)) "not run yet" else length(ss$cpgs_ge3)))
      ))
    })
    outputOptions(output, "qc_summary_ui", suspendWhenHidden = FALSE)

    output$dl_mr_results <- downloadHandler(
      filename = function() "methylomics_mr_results.csv",
      content = function(file) write.csv(build_primary_results_df(), file, row.names = FALSE)
    )
    output$dl_harmonised <- downloadHandler(
      filename = function() "methylomics_mr_harmonised_instruments.csv",
      content = function(file) write.csv(mr_state()$dat, file, row.names = FALSE)
    )
    output$dl_sensitivity <- downloadHandler(
      filename = function() "methylomics_mr_sensitivity.csv",
      content = function(file) {
        ss <- sensitivity_state()
        write.csv(if (!is.null(ss$het)) ss$het else data.frame(), file, row.names = FALSE)
      }
    )
    output$dl_qc_report <- downloadHandler(
      filename = function() "methylomics_mr_qc_report.txt",
      content = function(file) {
        cs <- clump_state(); hs <- harmonise_state()
        lines <- c(
          "Methylomics Mendelian Randomization - QC report",
          sprintf("Instruments before clumping: %s", if (is.null(cs) || is.na(cs$n_before)) "n/a (preloaded)" else cs$n_before),
          sprintf("Instruments after clumping: %s", if (is.null(cs)) "n/a" else cs$n_after),
          sprintf("Harmonisation - final instruments: %s", if (is.null(hs)) "n/a" else hs$summary$final)
        )
        writeLines(lines, file)
      }
    )
    output$dl_full_report <- downloadHandler(
      filename = function() "methylomics_mr_full_report.txt",
      content = function(file) {
        r <- build_primary_results_df()
        lines <- c(
          "Methylomics Mendelian Randomization - complete report",
          sprintf("Generated for %d CpG(s).", length(unique(r$exposure))),
          "", "Primary results:", utils::capture.output(print(r))
        )
        writeLines(lines, file)
      }
    )

    ## ------------------------------------------------------------------
    ## Stage 8: Plots
    ## ------------------------------------------------------------------
    output$plots_tab_body <- renderUI({
      req(stage_flags$mr)
      ms <- mr_state()
      tagList(
        fluidRow(
          column(6, selectInput(ns("plot_cpg"), "CpG (for single-CpG plots)", choices = ms$cpgs, selected = ms$cpgs[1], width = "100%")),
          column(6, selectInput(ns("plot_format"), "Download format", c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"), selected = "png", width = "100%"))
        ),
        tabsetPanel(
          type = "tabs",
          tabPanel("Scatter", br(), actionButton(ns("show_scatter"), "Show plot", class = "btn-sm btn-primary"),
                   uiOutput(ns("scatter_plot_ui"))),
          tabPanel("Forest", br(), actionButton(ns("show_forest"), "Show plot", class = "btn-sm btn-primary"),
                   uiOutput(ns("forest_plot_ui"))),
          tabPanel("Funnel", br(), actionButton(ns("show_funnel"), "Show plot", class = "btn-sm btn-primary"),
                   uiOutput(ns("funnel_plot_ui"))),
          tabPanel("Leave-one-out", br(), actionButton(ns("show_loo"), "Show plot", class = "btn-sm btn-primary"),
                   uiOutput(ns("loo_plot_ui"))),
          tabPanel("Single-SNP forest", br(), actionButton(ns("show_snp_forest"), "Show plot", class = "btn-sm btn-primary"),
                   uiOutput(ns("snp_forest_plot_ui"))),
          tabPanel("Manhattan overview", br(), actionButton(ns("show_manhattan"), "Show plot", class = "btn-sm btn-primary"),
                   uiOutput(ns("manhattan_plot_ui")))
        )
      )
    })
    outputOptions(output, "plots_tab_body", suspendWhenHidden = FALSE)

    show_flags <- reactiveValues(scatter = FALSE, forest = FALSE, funnel = FALSE, loo = FALSE, snp_forest = FALSE, manhattan = FALSE)
    observeEvent(input$show_scatter, show_flags$scatter <- TRUE)
    observeEvent(input$show_forest, show_flags$forest <- TRUE)
    observeEvent(input$show_funnel, show_flags$funnel <- TRUE)
    observeEvent(input$show_loo, show_flags$loo <- TRUE)
    observeEvent(input$show_snp_forest, show_flags$snp_forest <- TRUE)
    observeEvent(input$show_manhattan, show_flags$manhattan <- TRUE)

    cpg_sub <- reactive({ req(mr_state(), input$plot_cpg); mr_state()$dat[mr_state()$dat$exposure == input$plot_cpg, , drop = FALSE] })
    cpg_res <- reactive({ req(mr_state(), input$plot_cpg); mr_state()$results[mr_state()$results$exposure == input$plot_cpg, , drop = FALSE] })

    build_scatter <- function() { validate(need(nrow(cpg_sub()) > 0, "No data.")); TwoSampleMR::mr_scatter_plot(cpg_res(), cpg_sub())[[1]] + theme_arthomix(11) }
    build_forest <- function() { ss <- tryCatch(TwoSampleMR::mr_singlesnp(cpg_sub()), error = function(e) NULL); validate(need(!is.null(ss), "Not available.")); TwoSampleMR::mr_forest_plot(ss)[[1]] + theme_arthomix(11) }
    build_funnel <- function() { ss <- tryCatch(TwoSampleMR::mr_singlesnp(cpg_sub()), error = function(e) NULL); validate(need(!is.null(ss), "Not available.")); TwoSampleMR::mr_funnel_plot(ss)[[1]] + theme_arthomix(11) }
    build_loo <- function() { validate(need(nrow(cpg_sub()) >= 3, "Needs >=3 instruments.")); lo <- tryCatch(TwoSampleMR::mr_leaveoneout(cpg_sub()), error = function(e) NULL); validate(need(!is.null(lo), "Not available.")); TwoSampleMR::mr_leaveoneout_plot(lo)[[1]] + theme_arthomix(11) }
    build_manhattan <- function() {
      ms <- mr_state(); anno <- cpg_annotation(ms$cpgs)
      d <- mmr_primary_row_per_cpg(ms$results)
      d$pos <- suppressWarnings(as.numeric(anno[d$exposure, "pos"]))
      validate(need(length(unique(d$exposure)) > 1, "Needs more than one CpG."))
      validate(need(any(!is.na(d$pos)), "No genomic positions available for these CpGs."))
      ggplot(d, aes(x = pos, y = -log10(pval), label = exposure)) +
        geom_point(color = ARTHOMIX_COLORS$blue, size = 2) +
        ggrepel::geom_text_repel(size = 3) +
        geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey50") +
        labs(title = "CpG MR results overview", x = "Genomic position", y = "-log10(p), primary method") +
        theme_arthomix(11)
    }

    render_plot_slot <- function(show_flag, build_fn, dl_name) {
      renderUI({
        if (!isTRUE(show_flags[[show_flag]])) return(p(class = "empty-note", icon("circle-info"), "Click \"Show plot\" to render."))
        tagList(
          div(class = "table-toolbar", downloadButton(ns(paste0("dl_", show_flag)), "Download", class = "btn-sm")),
          withSpinner(plotOutput(ns(paste0("plot_", show_flag)), height = 420), color = "#2563EB", type = 6)
        )
      })
    }
    output$scatter_plot_ui <- render_plot_slot("scatter", build_scatter, "scatter"); outputOptions(output, "scatter_plot_ui", suspendWhenHidden = FALSE)
    output$forest_plot_ui <- render_plot_slot("forest", build_forest, "forest"); outputOptions(output, "forest_plot_ui", suspendWhenHidden = FALSE)
    output$funnel_plot_ui <- render_plot_slot("funnel", build_funnel, "funnel"); outputOptions(output, "funnel_plot_ui", suspendWhenHidden = FALSE)
    output$loo_plot_ui <- render_plot_slot("loo", build_loo, "loo"); outputOptions(output, "loo_plot_ui", suspendWhenHidden = FALSE)
    output$snp_forest_plot_ui <- render_plot_slot("snp_forest", build_forest, "snp_forest"); outputOptions(output, "snp_forest_plot_ui", suspendWhenHidden = FALSE)
    output$manhattan_plot_ui <- render_plot_slot("manhattan", build_manhattan, "manhattan"); outputOptions(output, "manhattan_plot_ui", suspendWhenHidden = FALSE)

    output$plot_scatter <- renderPlot(build_scatter())
    output$plot_forest <- renderPlot(build_forest())
    output$plot_funnel <- renderPlot(build_funnel())
    output$plot_loo <- renderPlot(build_loo())
    output$plot_snp_forest <- renderPlot(build_forest())
    output$plot_manhattan <- renderPlot(build_manhattan())

    make_plot_dl <- function(build_fn, base_name) downloadHandler(
      filename = function() sprintf("%s.%s", base_name, input$plot_format %||% "png"),
      content = function(file) {
        ggplot2::ggsave(file, plot = build_fn(), width = 9, height = 6, dpi = 300, device = input$plot_format %||% "png")
      }
    )
    output$dl_scatter <- make_plot_dl(build_scatter, "mr_scatter")
    output$dl_forest <- make_plot_dl(build_forest, "mr_forest")
    output$dl_funnel <- make_plot_dl(build_funnel, "mr_funnel")
    output$dl_loo <- make_plot_dl(build_loo, "mr_leaveoneout")
    output$dl_snp_forest <- make_plot_dl(build_forest, "mr_single_snp_forest")
    output$dl_manhattan <- make_plot_dl(build_manhattan, "mr_manhattan")

    NULL
  })
}
