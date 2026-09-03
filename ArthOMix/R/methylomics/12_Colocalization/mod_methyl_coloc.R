## R/methylomics/12_Colocalization/mod_methyl_coloc.R
## Bayesian colocalisation (coloc.abf, optionally coloc.susie) between an
## mQTL/CpG signal and a GWAS trait signal in one region: shared causal

.mcol_tip <- function(text) tags$span(icon("circle-info", style = "color:#8A929C; cursor: help; margin-left: 4px;"), title = text)

.mcol_stage_order <- c("validate", "run", "plot", "sensitivity")

MCOL_DEFAULT_P1 <- 1e-4
MCOL_DEFAULT_P2 <- 1e-4
MCOL_DEFAULT_P12 <- 1e-5
MCOL_DEFAULT_P12_SUSIE <- 5e-6
MCOL_DEFAULT_WINDOW_KB <- 1000
MCOL_DEFAULT_MIN_SHARED_SNPS <- 10
MCOL_DEFAULT_PP_THRESHOLD <- 0.8

.mcol_prep_ld <- function(raw_df, snp_ids) {
  if (is.null(raw_df) || nrow(raw_df) == 0 || ncol(raw_df) < 3) return(NULL)
  df <- as.data.frame(raw_df)
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- as.character(df[[1]])
  common <- intersect(as.character(snp_ids), intersect(rownames(mat), colnames(mat)))
  if (length(common) < 2) return(NULL)
  mat[common, common, drop = FALSE]
}

.mcol_interpret <- function(h0, h1, h2, h3, h4) {
  hs <- c(H0 = h0, H1 = h1, H2 = h2, H3 = h3, H4 = h4)
  strongest <- names(hs)[which.max(hs)]
  lead_p <- unname(hs[strongest])
  headline <- switch(strongest,
    H4 = if (h4 >= MCOL_DEFAULT_PP_THRESHOLD)
      "Strong evidence supporting a shared genetic signal between the methylation-associated signal and the GWAS trait in this region."
    else "Moderate evidence supporting a shared genetic signal between the two traits in this region - below the project's own 0.8 posterior threshold for a firm \"coloc-supported\" verdict.",
    H3 = if (h3 >= MCOL_DEFAULT_PP_THRESHOLD)
      "Both traits show association in this region, but the evidence favours two distinct, LD-linked causal variants rather than one shared signal."
    else "Both traits show some association in this region; the evidence leans toward distinct causal variants, though this is not strongly resolved.",
    H0 = "Little evidence of association with either trait in this region.",
    H1 = "Evidence primarily for the methylation-associated signal in this region, with little evidence for the GWAS trait.",
    H2 = "Evidence primarily for the GWAS trait in this region, with little evidence for the methylation-associated signal."
  )
  tagList(
    p(strong(sprintf("Strongest-supported hypothesis: %s", strongest)), sprintf(" (posterior probability %.1f%%).", 100 * lead_p)),
    p(headline),
    p(class = "empty-note", icon("circle-info"),
      "Colocalisation identifies statistical compatibility with a shared causal signal - it is not, by itself, proof of biological causality. A high PP.H4 does not mean a specific SNP causes the disease; it means the methylation- and disease-associated signals in this region are consistent with arising from the same underlying genetic variant.")
  )
}

.mcol_verdict <- function(h3, h4, threshold = MCOL_DEFAULT_PP_THRESHOLD) {
  if (h4 >= threshold) "coloc-supported (shared causal variant)"
  else if (h3 >= threshold) "coloc-refuted (distinct causal variants)"
  else sprintf("inconclusive (<%.2g posterior for either H3 or H4)", threshold)
}

mod_methyl_coloc_config <- list(
  id = "coloc", title = "Colocalisation", icon = "bullseye", group = "Genetics",
  description = "Tests whether an mQTL and a GWAS signal at a CpG's region share a causal variant"
)

mod_methyl_coloc_ui <- function(id) {
  ns <- NS(id)
  div(
    class = "tx-menu-wrap",
    tabsetPanel(
      id = ns("coloc_tabs"), type = "tabs",
      tabPanel("1. Data & Setup", br(), mcol_data_ui(ns)),
      tabPanel("2. Filters & Parameters", br(), mcol_filters_ui(ns)),
      tabPanel("3. Results", br(), mcol_results_ui(ns)),
      tabPanel("4. Visualisation", br(), mcol_plots_ui(ns)),
      tabPanel("5. Sensitivity Analysis", br(), mcol_sensitivity_ui(ns)),
      tabPanel("6. Export", br(), mcol_export_ui(ns))
    )
  )
}

mcol_data_ui <- function(ns) {
  fluidRow(
    column(
      4,
      box(
        width = NULL, title = "Data source", status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "Colocalisation between a methylation-associated genetic signal (mQTL) and a GWAS trait's signal at one genomic region - not eQTL/eGene colocalisation."),
        radioButtons(
          ns("data_source"), NULL,
          choiceNames = list(tagList(icon("database"), " Preloaded Data"), tagList(icon("upload"), " Upload Data")),
          choiceValues = list("preloaded", "upload"), selected = "preloaded"
        )
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'preloaded'", ns("data_source")),
        box(
          width = NULL, title = "Preloaded methylomics colocalisation dataset", status = "primary", solidHeader = FALSE,
          if (METH_DATA_AVAILABLE) tagList(
            p(class = "submodule-desc", "Reproduces ", strong("script08_mendelian_randomization"), "'s coloc.abf() run: the GoDMC cis-mQTL signal vs the Ishigaki et al. (2022) rheumatoid arthritis GWAS signal, at each CpG carried into MR."),
            uiOutput(ns("pre_cpg_ui")),
            p(class = "empty-note", icon("circle-info"),
              "Only the completed run's PP.H0-H4 summary per CpG is bundled with this deployment - the underlying per-SNP GoDMC/RA-GWAS region data isn't, so coloc.abf() can't be re-run live here. Results below are looked up, not recomputed, and SNP-level output, regional plots, and prior-sensitivity analysis are unavailable for this route. Use Upload Data for a fully live analysis.")
          ) else p(class = "empty-note", icon("triangle-exclamation"),
                   "The preloaded methylomics colocalisation dataset isn't available in this deployment - use Upload Data instead.")
        )
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'upload'", ns("data_source")),
        box(
          width = NULL, title = "Dataset 1: Methylation / mQTL data", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "One row per SNP, optionally x CpG (long format). CSV/TSV. Methylation is treated as a quantitative trait - map the effect-allele-frequency column below (coloc.abf needs it to analyse a quantitative trait without a directly-known phenotype SD)."),
          fileInput(ns("meth_file"), "Methylation/mQTL file", accept = c(".csv", ".tsv", ".txt")),
          uiOutput(ns("meth_map_ui")),
          uiOutput(ns("meth_extra_map_ui"))
        ),
        box(
          width = NULL, title = "Dataset 2: GWAS / trait data", status = "primary", solidHeader = FALSE,
          textInput(ns("gwas_label"), "GWAS trait name (for labelling only)", value = "Uploaded GWAS trait", width = "100%"),
          radioButtons(ns("gwas_type"), "Trait type", inline = TRUE,
                       choices = c("Binary (case/control)" = "cc", "Quantitative" = "quant"), selected = "cc"),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'cc'", ns("gwas_type")),
            sliderInput(ns("case_frac"), "Case fraction in the GWAS", value = 0.33, min = 0.01, max = 0.99, step = 0.01)
          ),
          fileInput(ns("gwas_file"), "GWAS file", accept = c(".csv", ".tsv", ".txt")),
          uiOutput(ns("gwas_map_ui"))
        )
      )
    ),
    column(
      8,
      uiOutput(ns("preview_ui")),
      uiOutput(ns("validation_ui")),
      div(style = "margin-top: 10px;",
          actionButton(ns("validate_btn"), "Validate Data", icon = icon("check"), class = "btn-primary btn-sm"))
    )
  )
}

mcol_filters_ui <- function(ns) uiOutput(ns("filters_tab_body"))

mcol_filters_controls_preloaded <- function(ns) {
  tagList(
    p(class = "empty-note", icon("circle-info"),
      sprintf("The Preloaded route reproduces the completed coloc.abf() run per CpG - priors (p1=%.0e, p2=%.0e, p12=%.0e, coloc's own defaults) and the +/-1Mb GoDMC cis window were fixed when this cached result was produced and are not re-run live for this route. Only a minimum-instrument-count filter and a focus-CpG choice apply here.",
              MCOL_DEFAULT_P1, MCOL_DEFAULT_P2, MCOL_DEFAULT_P12)),
    numericInput(ns("pre_min_nsnps"), "Minimum SNPs used per CpG", value = 0, min = 0, step = 10, width = "100%"),
    div(style = "margin-top: 6px;",
        actionButton(ns("run_btn"), "Run Colocalisation", icon = icon("play"), class = "btn-primary btn-sm"))
  )
}

mcol_filters_controls_upload <- function(ns) {
  tagList(
    tags$h5("Genomic-region filters"),
    checkboxInput(ns("use_window"), "Restrict to a genomic window", value = FALSE),
    conditionalPanel(
      condition = sprintf("input['%s']", ns("use_window")),
      radioButtons(ns("window_center"), "Window centred on", inline = TRUE,
                   choices = c("Lead SNP (smallest mQTL p-value)" = "lead", "Specified position" = "manual"), selected = "lead"),
      conditionalPanel(condition = sprintf("input['%s'] == 'manual'", ns("window_center")),
                       numericInput(ns("window_pos"), "Centre position (bp)", value = NA, width = "100%")),
      numericInput(ns("window_kb"), "Window (+/- kb)", value = MCOL_DEFAULT_WINDOW_KB, min = 1, step = 100, width = "100%")
    ),
    div(style = "display:flex; align-items:center; gap:4px;",
        numericInput(ns("f_min_shared"), "Minimum shared SNPs required", value = MCOL_DEFAULT_MIN_SHARED_SNPS, min = 3, step = 1, width = "100%"),
        .mcol_tip("coloc needs a minimally informative set of SNPs tested across the region to distinguish PP.H3 from PP.H4 - this project's own pipeline required >=10 GoDMC SNPs per CpG.")),
    tags$h5("Association filters"),
    fluidRow(
      column(6, numericInput(ns("f_pval_meth"), "mQTL p-value threshold", value = 1, min = 0, max = 1, step = 0.0001, width = "100%")),
      column(6, numericInput(ns("f_pval_gwas"), "GWAS p-value threshold", value = 1, min = 0, max = 1, step = 0.0001, width = "100%"))
    ),
    fluidRow(
      column(6, numericInput(ns("f_maf_min"), "Minimum MAF", value = 0, min = 0, max = 0.5, step = 0.01, width = "100%")),
      column(6, numericInput(ns("f_maf_max"), "Maximum MAF", value = 0.5, min = 0, max = 0.5, step = 0.01, width = "100%"))
    ),
    numericInput(ns("f_min_n"), "Minimum sample size", value = 0, min = 0, step = 100, width = "100%"),
    tags$h5("Variant filtering"),
    checkboxInput(ns("f_dedup"), "Remove duplicated SNPs", value = TRUE),
    checkboxInput(ns("remove_ambiguous"), "Remove ambiguous/palindromic alleles (harmonisation strategy 3) instead of inferring via EAF (strategy 2)", value = FALSE),
    tags$h5("Priors"),
    fluidRow(
      column(4, div(style = "display:flex; align-items:center; gap:4px;",
                    numericInput(ns("p1"), "p1", value = MCOL_DEFAULT_P1, min = 0, max = 1, step = 1e-6, width = "100%"),
                    .mcol_tip("Prior probability a SNP is associated with trait 1 (methylation) alone. coloc default: 1e-4."))),
      column(4, div(style = "display:flex; align-items:center; gap:4px;",
                    numericInput(ns("p2"), "p2", value = MCOL_DEFAULT_P2, min = 0, max = 1, step = 1e-6, width = "100%"),
                    .mcol_tip("Prior probability a SNP is associated with trait 2 (GWAS) alone. coloc default: 1e-4."))),
      column(4, div(style = "display:flex; align-items:center; gap:4px;",
                    numericInput(ns("p12"), "p12", value = MCOL_DEFAULT_P12, min = 0, max = 1, step = 1e-7, width = "100%"),
                    .mcol_tip("Prior probability a SNP is associated with BOTH traits. coloc default: 1e-5.")))
    ),
    tags$h5("Method"),
    checkboxInput(ns("use_susie"), "Also run multiple-signal colocalisation (coloc.susie)", value = FALSE),
    conditionalPanel(
      condition = sprintf("input['%s']", ns("use_susie")),
      div(class = "empty-note", icon("circle-info"),
          "coloc.susie() needs a SNP-by-SNP LD (correlation) matrix for BOTH datasets - upload a delimited file with SNP IDs as the first column and as the header row, matching the SNP ID column mapped above. No LD is inferred, estimated, or simulated if not provided; without both matrices this option stays unavailable."),
      fileInput(ns("ld1_file"), "LD matrix - methylation/mQTL SNPs", accept = c(".csv", ".tsv")),
      fileInput(ns("ld2_file"), "LD matrix - GWAS SNPs", accept = c(".csv", ".tsv")),
      fluidRow(
        column(6, numericInput(ns("susie_coverage"), "Credible-set coverage", value = 0.95, min = 0.5, max = 0.99, step = 0.01, width = "100%")),
        column(6, numericInput(ns("susie_maxit"), "Max iterations", value = 100, min = 10, step = 10, width = "100%"))
      )
    ),
    tags$details(
      class = "box box-primary", style = "margin-top: 10px;",
      tags$summary(class = "box-header", style = "cursor: pointer;", tags$h3(class = "box-title", "Advanced parameters")),
      div(class = "box-body",
        numericInput(ns("pp_threshold"), "Posterior-probability threshold for a supported/refuted verdict",
                     value = MCOL_DEFAULT_PP_THRESHOLD, min = 0.5, max = 0.99, step = 0.01, width = "100%"),
        p(class = "submodule-desc", style = "font-size:12px;",
          "Matches this project's own script08d threshold: PP.H4 >= this value -> \"coloc-supported\"; PP.H3 >= this value -> \"coloc-refuted\"; otherwise \"inconclusive\".")
      )
    ),
    div(style = "margin-top: 10px;",
        actionButton(ns("run_btn"), "Run Colocalisation", icon = icon("play"), class = "btn-primary btn-sm"))
  )
}

mcol_results_ui <- function(ns) uiOutput(ns("results_tab_body"))

mcol_plots_ui <- function(ns) uiOutput(ns("plots_tab_body"))

mcol_sensitivity_ui <- function(ns) uiOutput(ns("sensitivity_tab_body"))

mcol_export_ui <- function(ns) uiOutput(ns("export_tab_body"))

mod_methyl_coloc_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    stage_flags <- reactiveValues(validate = FALSE, run = FALSE, plot = FALSE, sensitivity = FALSE)
    invalidate_from <- function(stage) {
      idx <- match(stage, .mcol_stage_order)
      for (s in .mcol_stage_order[idx:length(.mcol_stage_order)]) stage_flags[[s]] <- FALSE
    }

    pre_coloc_tbl <- reactive({
      req(METH_DATA_AVAILABLE)
      tbl <- load_default_meth_coloc_results()
      validate(need(!is.null(tbl), "Preloaded coloc results not found."))
      tbl
    })

    output$pre_cpg_ui <- renderUI({
      req(METH_DATA_AVAILABLE)
      tbl <- pre_coloc_tbl()
      cpgs <- sort(unique(tbl$cpg))
      tagList(
        selectInput(ns("pre_cpgs"), "CpG(s)", choices = cpgs, selected = cpgs, multiple = TRUE, width = "100%"),
        p(class = "submodule-desc", style = "font-size:12px; margin-bottom:0;",
          sprintf("%d CpG(s) have a bundled coloc.abf() result (>=10 GoDMC candidate SNPs in the cis window).", length(cpgs)))
      )
    })

    meth_df_r <- reactive({ req(input$meth_file); read_uploaded_table(input$meth_file$datapath) })
    gwas_df_r <- reactive({ req(input$gwas_file); read_uploaded_table(input$gwas_file$datapath) })

    output$meth_map_ui <- gwas_col_map_ui(ns, reactive(input$meth_file), meth_df_r, "meth", "Methylation/mQTL file", extra_fields = "n")
    output$gwas_map_ui <- gwas_col_map_ui(ns, reactive(input$gwas_file), gwas_df_r, "gwas", "GWAS file", extra_fields = "n")

    output$meth_extra_map_ui <- renderUI({
      req(input$meth_file)
      df <- meth_df_r()
      validate(need(!is.null(df), "Could not read the methylation/mQTL file."))
      cols <- colnames(df)
      guess <- function(patterns) { g <- guess_gwas_col(cols, patterns); if (is.na(g)) NA_character_ else g }
      cpg_guess <- guess(c("^cpg", "^probe", "^cg_?id"))
      tagList(
        selectInput(ns("meth_cpg"), "CpG ID column", cols, selected = if (!is.na(cpg_guess)) cpg_guess else cols[1]),
        p(class = "submodule-desc", style = "font-size:12px;", "Optional - needed for the genomic-window region filter and the regional/comparison plots:"),
        fluidRow(
          column(6, selectInput(ns("meth_snp_chr"), "SNP chromosome", c("(none)" = "", cols), selected = { g <- guess(c("^chr", "^chromosome")); if (is.na(g)) "" else g })),
          column(6, selectInput(ns("meth_snp_pos"), "SNP position", c("(none)" = "", cols), selected = { g <- guess(c("^pos", "^bp$", "base_pair")); if (is.na(g)) "" else g }))
        ),
        selectInput(ns("meth_gene"), "Gene (optional annotation)", c("(none)" = "", cols), selected = { g <- guess(c("^gene")); if (is.na(g)) "" else g }),
        uiOutput(ns("meth_target_cpg_ui"))
      )
    })

    output$meth_target_cpg_ui <- renderUI({
      req(input$meth_file, input$meth_cpg)
      df <- meth_df_r(); req(df)
      cpg_col <- input$meth_cpg
      req(cpg_col %in% colnames(df))
      cpgs <- sort(unique(as.character(df[[cpg_col]])))
      if (length(cpgs) <= 1) {
        return(p(class = "empty-note", icon("circle-info"),
                  sprintf("Single CpG detected: %s.", if (length(cpgs) == 1) cpgs[1] else "(none)")))
      }
      tagList(
        selectInput(ns("meth_target_cpg"), "Target CpG (this file has multiple)", cpgs, selected = cpgs[1]),
        p(class = "submodule-desc", style = "font-size:12px;",
          sprintf("%d CpGs found in this file - colocalisation is run for one CpG's region at a time.", length(cpgs)))
      )
    })

    observeEvent(list(input$data_source, input$pre_cpgs, input$meth_file, input$gwas_file,
                       input$meth_cpg, input$meth_target_cpg, input$gwas_label, input$gwas_type, input$case_frac),
                 invalidate_from("validate"), ignoreInit = TRUE, ignoreNULL = FALSE)

    output$preview_ui <- renderUI({
      if (identical(input$data_source, "upload")) {
        tagList(
          if (!is.null(input$meth_file)) {
            df <- meth_df_r()
            if (!is.null(df)) box(width = NULL, title = "Methylation/mQTL file preview", status = "primary", solidHeader = FALSE,
                                   DT::dataTableOutput(ns("meth_preview_table")))
          },
          if (!is.null(input$gwas_file)) {
            df <- gwas_df_r()
            if (!is.null(df)) box(width = NULL, title = "GWAS file preview", status = "primary", solidHeader = FALSE,
                                   DT::dataTableOutput(ns("gwas_preview_table")))
          }
        )
      } else {
        req(METH_DATA_AVAILABLE)
        cpgs <- input$pre_cpgs
        req(length(cpgs) > 0)
        tbl <- pre_coloc_tbl()
        tbl <- tbl[tbl$cpg %in% cpgs, c("cpg", "nsnps"), drop = FALSE]
        box(width = NULL, title = "Available preloaded data (CpG, SNPs used)", status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "Descriptive only - hypothesis probabilities remain hidden until Run Colocalisation."),
            DT::dataTableOutput(ns("pre_preview_table")))
      }
    })
    outputOptions(output, "preview_ui", suspendWhenHidden = FALSE)

    output$pre_preview_table <- DT::renderDataTable({
      req(METH_DATA_AVAILABLE)
      cpgs <- input$pre_cpgs
      req(length(cpgs) > 0)
      tbl <- pre_coloc_tbl()
      tbl <- tbl[tbl$cpg %in% cpgs, c("cpg", "nsnps"), drop = FALSE]
      DT::datatable(tbl, rownames = FALSE, options = list(pageLength = 8, dom = "tp"), class = "stripe hover compact")
    })

    output$meth_preview_table <- DT::renderDataTable({
      req(input$meth_file)
      df <- meth_df_r(); validate(need(!is.null(df), "Could not read the methylation/mQTL file."))
      DT::datatable(utils::head(df, 20), rownames = FALSE, options = list(scrollX = TRUE, dom = "tp"), class = "stripe hover compact")
    })
    output$gwas_preview_table <- DT::renderDataTable({
      req(input$gwas_file)
      df <- gwas_df_r(); validate(need(!is.null(df), "Could not read the GWAS file."))
      DT::datatable(utils::head(df, 20), rownames = FALSE, options = list(scrollX = TRUE, dom = "tp"), class = "stripe hover compact")
    })

    validate_state <- reactiveVal(NULL)

    build_validate_state_preloaded <- function() {
      req(METH_DATA_AVAILABLE)
      cpgs <- input$pre_cpgs
      validate(need(length(cpgs) > 0, "Pick at least one CpG."))
      tbl <- pre_coloc_tbl()
      tbl <- tbl[tbl$cpg %in% cpgs, , drop = FALSE]
      validate(need(nrow(tbl) > 0, "No preloaded coloc result for the selected CpG(s)."))
      harm <- load_default_mr_harmonised()
      ctx <- if (!is.null(harm)) harm[harm$exposure %in% cpgs & harm$mr_keep %in% TRUE, , drop = FALSE] else NULL

      list(
        mode = "preloaded", cpgs = cpgs, tbl = tbl, ctx = ctx,
        outcome_label = if (!is.null(ctx) && nrow(ctx) > 0) ctx$outcome[1] else "Rheumatoid arthritis (Ishigaki et al. 2022)",
        summary = list(
          dataset1_variants = if (!is.null(ctx)) length(unique(ctx$SNP)) else NA_integer_,
          dataset2_variants = if (!is.null(ctx)) length(unique(ctx$SNP)) else NA_integer_,
          shared_variants = sum(tbl$nsnps), duplicate_variants = 0L, missing_values = 0L,
          final_variants = sum(tbl$nsnps)
        )
      )
    }

    build_validate_state_upload <- function() {
      req(input$meth_file, input$gwas_file, input$meth_cpg,
          input$meth_snp, input$meth_beta, input$meth_se, input$meth_pval, input$meth_ea, input$meth_oa, input$meth_n,
          input$gwas_snp, input$gwas_beta, input$gwas_se, input$gwas_pval, input$gwas_ea, input$gwas_oa, input$gwas_n)
      meth_raw <- meth_df_r(); gwas_raw <- gwas_df_r()
      validate(need(!is.null(meth_raw), "Could not read the methylation/mQTL file."))
      validate(need(!is.null(gwas_raw), "Could not read the GWAS file."))

      cpg_col <- input$meth_cpg
      req(cpg_col %in% colnames(meth_raw))
      target_cpg <- input$meth_target_cpg %||% NA_character_
      meth_sub <- if (!is.na(target_cpg) && nzchar(target_cpg)) meth_raw[as.character(meth_raw[[cpg_col]]) == target_cpg, , drop = FALSE] else meth_raw
      validate(need(nrow(meth_sub) > 0, "No rows for the selected CpG in the methylation/mQTL file."))
      cpg_label <- if (!is.na(target_cpg) && nzchar(target_cpg)) target_cpg else as.character(meth_sub[[cpg_col]][1])

      n_dataset1_variants <- length(unique(meth_sub[[input$meth_snp]]))
      n_dataset2_variants <- length(unique(gwas_raw[[input$gwas_snp]]))
      dup1 <- sum(duplicated(meth_sub[[input$meth_snp]]))
      dup2 <- sum(duplicated(gwas_raw[[input$gwas_snp]]))
      req_cols1 <- c(input$meth_snp, input$meth_beta, input$meth_se, input$meth_pval)
      req_cols2 <- c(input$gwas_snp, input$gwas_beta, input$gwas_se, input$gwas_pval)
      miss1 <- sum(!stats::complete.cases(meth_sub[, req_cols1[req_cols1 %in% colnames(meth_sub)], drop = FALSE]))
      miss2 <- sum(!stats::complete.cases(gwas_raw[, req_cols2[req_cols2 %in% colnames(gwas_raw)], drop = FALSE]))
      invalid_p1 <- sum(!is.na(meth_sub[[input$meth_pval]]) & (meth_sub[[input$meth_pval]] < 0 | meth_sub[[input$meth_pval]] > 1))
      invalid_p2 <- sum(!is.na(gwas_raw[[input$gwas_pval]]) & (gwas_raw[[input$gwas_pval]] < 0 | gwas_raw[[input$gwas_pval]] > 1))
      overlap_raw <- length(intersect(meth_sub[[input$meth_snp]], gwas_raw[[input$gwas_snp]]))
      validate(need(overlap_raw > 0, "No overlapping SNP IDs between the methylation/mQTL and GWAS files - check that both use the same SNP identifier convention (e.g. rsIDs)."))

      fmt1_args <- list(meth_sub, type = "exposure", snp_col = input$meth_snp, beta_col = input$meth_beta, se_col = input$meth_se,
                         pval_col = input$meth_pval, effect_allele_col = input$meth_ea, other_allele_col = input$meth_oa,
                         samplesize_col = input$meth_n, phenotype_col = cpg_col, id_col = cpg_col)
      if (nzchar(input$meth_eaf %||% "")) fmt1_args$eaf_col <- input$meth_eaf
      if (nzchar(input$meth_snp_chr %||% "")) fmt1_args$chr_col <- input$meth_snp_chr
      if (nzchar(input$meth_snp_pos %||% "")) fmt1_args$pos_col <- input$meth_snp_pos
      meth_fmt <- tryCatch(do.call(TwoSampleMR::format_data, fmt1_args), error = function(e) NULL)
      validate(need(!is.null(meth_fmt) && nrow(meth_fmt) > 0, "No usable rows after formatting the methylation/mQTL file - check column mapping."))
      meth_fmt$exposure <- cpg_label

      fmt2_args <- list(gwas_raw, type = "outcome", snp_col = input$gwas_snp, beta_col = input$gwas_beta, se_col = input$gwas_se,
                         pval_col = input$gwas_pval, effect_allele_col = input$gwas_ea, other_allele_col = input$gwas_oa,
                         samplesize_col = input$gwas_n)
      if (nzchar(input$gwas_eaf %||% "")) fmt2_args$eaf_col <- input$gwas_eaf
      gwas_fmt <- tryCatch(do.call(TwoSampleMR::format_data, fmt2_args), error = function(e) NULL)
      validate(need(!is.null(gwas_fmt) && nrow(gwas_fmt) > 0, "No usable rows after formatting the GWAS file - check column mapping."))
      gwas_label <- if (nzchar(trimws(input$gwas_label %||% ""))) trimws(input$gwas_label) else "Uploaded GWAS trait"
      gwas_fmt$outcome <- gwas_label

      harm_action <- if (isTRUE(input$remove_ambiguous)) 3L else 2L
      harmonised <- tryCatch(TwoSampleMR::harmonise_data(meth_fmt, gwas_fmt, action = harm_action), error = function(e) NULL)
      validate(need(!is.null(harmonised) && nrow(harmonised) > 0,
        "Harmonisation found no overlapping, alignable SNPs between the methylation/mQTL and GWAS datasets - check that both use the same SNP identifiers (rsIDs) and that allele columns are mapped correctly."))

      summary_counts <- list(
        dataset1_variants = n_dataset1_variants, dataset2_variants = n_dataset2_variants,
        shared_variants = overlap_raw, duplicate_variants = dup1 + dup2, missing_values = miss1 + miss2,
        invalid_pvalues = invalid_p1 + invalid_p2,
        aligned = sum(harmonised$mr_keep %in% TRUE), palindromic = sum(harmonised$palindromic %in% TRUE),
        ambiguous = sum(harmonised$ambiguous %in% TRUE), removed = sum(!(harmonised$mr_keep %in% TRUE)),
        final_variants = sum(harmonised$mr_keep %in% TRUE)
      )
      harmonised <- harmonised[harmonised$mr_keep %in% TRUE, , drop = FALSE]
      validate(need(nrow(harmonised) > 0, "No variants survived harmonisation (all ambiguous/unresolvable, e.g. unresolved palindromic SNPs)."))

      list(mode = "upload", harmonised = harmonised, summary = summary_counts, cpg = cpg_label,
           gwas_label = gwas_label, gwas_type = input$gwas_type %||% "cc", case_frac = input$case_frac %||% 0.33,
           has_chr_pos = all(c("chr.exposure", "pos.exposure") %in% colnames(harmonised)) && any(!is.na(harmonised$pos.exposure)))
    }

    observeEvent(input$validate_btn, {
      st <- tryCatch(
        if (identical(input$data_source, "upload")) build_validate_state_upload() else build_validate_state_preloaded(),
        error = function(e) e
      )
      if (inherits(st, "shiny.silent.error")) {
        showNotification(conditionMessage(st), type = "warning", duration = 10)
        return()
      }
      if (inherits(st, "error")) {
        showNotification(sprintf("Validate Data failed: %s", conditionMessage(st)), type = "error", duration = 10)
        return()
      }
      validate_state(st)
      stage_flags$validate <- TRUE
      invalidate_from("run")
      showNotification("Validation complete - see the summary below, then continue to Filters & Parameters.", type = "message", duration = 6)
    })

    output$validation_ui <- renderUI({
      req(stage_flags$validate, validate_state())
      vs <- validate_state(); s <- vs$summary
      box(
        width = NULL, title = "Data validation / harmonisation summary", status = "primary", solidHeader = FALSE,
        div(class = "methyl-stats-row", fluidRow(
          valueBox(format(s$dataset1_variants %||% NA, big.mark = ","), "Dataset 1 variants", icon = icon("dna"), color = "light-blue", width = 3),
          valueBox(format(s$dataset2_variants %||% NA, big.mark = ","), "Dataset 2 variants", icon = icon("code-branch"), color = "purple", width = 3),
          valueBox(format(s$shared_variants %||% NA, big.mark = ","), "Shared variants", icon = icon("circle-check"), color = "green", width = 3),
          valueBox(format(s$final_variants %||% NA, big.mark = ","), "Final variants available", icon = icon("check-double"), color = "aqua", width = 3)
        )),
        tags$table(class = "table table-sm", tags$tbody(
          tags$tr(tags$td("Duplicate variants"), tags$td(s$duplicate_variants %||% "n/a")),
          tags$tr(tags$td("Missing required values"), tags$td(s$missing_values %||% "n/a")),
          if (!is.null(s$invalid_pvalues)) tags$tr(tags$td("Invalid p-values (outside [0,1])"), tags$td(s$invalid_pvalues)),
          if (!is.null(s$aligned)) tags$tr(tags$td("Allele-aligned (harmonised)"), tags$td(s$aligned)),
          if (!is.null(s$palindromic)) tags$tr(tags$td("Palindromic"), tags$td(s$palindromic)),
          if (!is.null(s$ambiguous)) tags$tr(tags$td("Ambiguous"), tags$td(s$ambiguous)),
          if (!is.null(s$removed)) tags$tr(tags$td("Removed at harmonisation"), tags$td(s$removed))
        )),
        if (identical(vs$mode, "preloaded")) p(class = "empty-note", icon("circle-info"),
          "Preloaded route: this reflects the already-completed upstream pipeline's own inputs (GoDMC cis-mQTL candidate rows, Ishigaki 2022 RA GWAS) - no re-validation is performed live, since the underlying per-SNP rows aren't bundled with this deployment.")
      )
    })
    outputOptions(output, "validation_ui", suspendWhenHidden = FALSE)

    output$filters_tab_body <- renderUI({
      req(stage_flags$validate)
      vs <- validate_state()
      if (identical(vs$mode, "preloaded")) mcol_filters_controls_preloaded(ns) else mcol_filters_controls_upload(ns)
    })
    outputOptions(output, "filters_tab_body", suspendWhenHidden = FALSE)

    observeEvent(list(input$use_window, input$window_center, input$window_pos, input$window_kb, input$f_min_shared,
                       input$f_pval_meth, input$f_pval_gwas, input$f_maf_min, input$f_maf_max, input$f_min_n,
                       input$f_dedup, input$remove_ambiguous, input$p1, input$p2, input$p12, input$use_susie,
                       input$ld1_file, input$ld2_file, input$susie_coverage, input$susie_maxit, input$pp_threshold,
                       input$pre_min_nsnps),
                 invalidate_from("run"), ignoreInit = TRUE, ignoreNULL = FALSE)

    run_state <- reactiveVal(NULL)

    build_run_state_preloaded <- function() {
      vs <- validate_state(); req(vs)
      tbl <- vs$tbl
      min_nsnps <- input$pre_min_nsnps %||% 0
      tbl <- tbl[tbl$nsnps >= min_nsnps, , drop = FALSE]
      validate(need(nrow(tbl) > 0, "No CpGs remain after the minimum-SNP-count filter."))
      focus <- tbl$cpg[1]
      list(mode = "preloaded", table = tbl, focus_cpg = focus, outcome_label = vs$outcome_label,
           n_variants = sum(tbl$nsnps), cpgs = tbl$cpg)
    }

    build_run_state_upload <- function() {
      vs <- validate_state(); req(vs)
      h_full <- vs$harmonised
      h <- h_full

      if (isTRUE(input$use_window) && vs$has_chr_pos) {
        pos <- h$pos.exposure
        center <- if (identical(input$window_center, "manual") && !is.na(input$window_pos)) input$window_pos
                  else pos[which.min(h$pval.exposure)]
        win_bp <- (input$window_kb %||% MCOL_DEFAULT_WINDOW_KB) * 1000
        h <- h[!is.na(pos) & abs(pos - center) <= win_bp, , drop = FALSE]
      }
      h <- h[h$pval.exposure <= (input$f_pval_meth %||% 1) & h$pval.outcome <= (input$f_pval_gwas %||% 1), , drop = FALSE]
      maf_val <- if (!is.null(h$eaf.exposure)) pmin(h$eaf.exposure, 1 - h$eaf.exposure) else rep(NA_real_, nrow(h))
      maf_lo <- input$f_maf_min %||% 0; maf_hi <- input$f_maf_max %||% 0.5
      h <- h[is.na(maf_val) | (maf_val >= maf_lo & maf_val <= maf_hi), , drop = FALSE]
      if (isTRUE(input$f_dedup)) h <- h[!duplicated(h$SNP), , drop = FALSE]
      min_n <- input$f_min_n %||% 0
      h <- h[(is.na(h$samplesize.exposure) | h$samplesize.exposure >= min_n) &
             (is.na(h$samplesize.outcome) | h$samplesize.outcome >= min_n), , drop = FALSE]

      min_shared <- input$f_min_shared %||% MCOL_DEFAULT_MIN_SHARED_SNPS
      validate(need(nrow(h) >= min_shared,
        sprintf("Fewer than %d shared SNPs remain after filtering (%d available) - colocalisation needs a minimally informative set of SNPs across the region. Loosen the filters or widen the genomic window.",
                min_shared, nrow(h))))
      validate(need(all(h$se.exposure > 0) && all(h$se.outcome > 0), "Some retained variants have a zero or negative standard error - cannot compute a Bayes factor for them."))

      p1 <- input$p1 %||% MCOL_DEFAULT_P1; p2 <- input$p2 %||% MCOL_DEFAULT_P2; p12 <- input$p12 %||% MCOL_DEFAULT_P12
      validate(need(is.numeric(p1) && length(p1) == 1 && !is.na(p1) && p1 > 0 && p1 < 1,
                    "p1 must be a probability strictly between 0 and 1."))
      validate(need(is.numeric(p2) && length(p2) == 1 && !is.na(p2) && p2 > 0 && p2 < 1,
                    "p2 must be a probability strictly between 0 and 1."))
      validate(need(is.numeric(p12) && length(p12) == 1 && !is.na(p12) && p12 > 0 && p12 < 1,
                    "p12 must be a probability strictly between 0 and 1."))
      validate(need(p12 <= min(p1, p2),
                    "p12 (probability a SNP affects both the methylation signal and the GWAS trait) cannot exceed p1 or p2 (probability it affects only one) - this is a coloc sanity convention. Lower p12, or raise p1/p2."))
      d1 <- list(beta = h$beta.exposure, varbeta = h$se.exposure^2,
                 N = round(stats::median(h$samplesize.exposure, na.rm = TRUE)), type = "quant", snp = h$SNP)
      if (!is.null(h$eaf.exposure)) d1$MAF <- pmin(h$eaf.exposure, 1 - h$eaf.exposure)
      validate(need(!is.null(d1$MAF),
        "coloc.abf needs a minor allele frequency to analyse a quantitative trait (methylation) without a directly-supplied phenotype SD - map an effect-allele-frequency column for the methylation/mQTL file in Data & Setup and re-validate."))
      d2 <- list(beta = h$beta.outcome, varbeta = h$se.outcome^2,
                 N = round(stats::median(h$samplesize.outcome, na.rm = TRUE)), type = vs$gwas_type, snp = h$SNP)
      if (identical(vs$gwas_type, "cc")) d2$s <- vs$case_frac
      else if (!is.null(h$eaf.outcome)) d2$MAF <- pmin(h$eaf.outcome, 1 - h$eaf.outcome)

      abf_res <- tryCatch(suppressWarnings(coloc::coloc.abf(dataset1 = d1, dataset2 = d2, p1 = p1, p2 = p2, p12 = p12)), error = function(e) e)
      validate(need(!inherits(abf_res, "error"), sprintf("coloc.abf() failed: %s", if (inherits(abf_res, "error")) conditionMessage(abf_res) else "unknown error")))

      snp_df <- data.frame(
        snp = h$SNP,
        chr = if (!is.null(h$chr.exposure)) h$chr.exposure else NA_character_,
        pos = if (!is.null(h$pos.exposure)) h$pos.exposure else NA_real_,
        effect_allele = h$effect_allele.exposure, other_allele = h$other_allele.exposure,
        maf = if (!is.null(h$eaf.exposure)) round(pmin(h$eaf.exposure, 1 - h$eaf.exposure), 4) else NA_real_,
        gwas_beta = h$beta.outcome, gwas_se = h$se.outcome, gwas_pval = h$pval.outcome,
        meth_beta = h$beta.exposure, meth_se = h$se.exposure, meth_pval = h$pval.exposure,
        log_abf = abf_res$results$internal.sum.lABF[match(h$SNP, abf_res$results$snp)],
        snp_pp_h4 = abf_res$results$SNP.PP.H4[match(h$SNP, abf_res$results$snp)],
        stringsAsFactors = FALSE
      )
      snp_df <- snp_df[order(-snp_df$snp_pp_h4), , drop = FALSE]
      snp_df$rank <- seq_len(nrow(snp_df))
      lead_variant <- snp_df$snp[which.max(snp_df$snp_pp_h4)]

      susie_res <- NULL; susie_note <- NULL
      if (isTRUE(input$use_susie)) {
        if (!requireNamespace("susieR", quietly = TRUE)) {
          susie_note <- "susieR is not installed in this deployment - multiple-signal colocalisation is unavailable."
        } else if (is.null(input$ld1_file) || is.null(input$ld2_file)) {
          susie_note <- "Multiple-signal colocalisation (coloc.susie) needs an LD (SNP-by-SNP correlation) matrix for both datasets - upload both to enable it. No LD was inferred or fabricated, so this method was not run."
        } else {
          ld1_raw <- read_uploaded_table(input$ld1_file$datapath)
          ld2_raw <- read_uploaded_table(input$ld2_file$datapath)
          ld1 <- .mcol_prep_ld(ld1_raw, h$SNP); ld2 <- .mcol_prep_ld(ld2_raw, h$SNP)
          common_ld <- if (!is.null(ld1) && !is.null(ld2)) intersect(colnames(ld1), colnames(ld2)) else character(0)
          if (length(common_ld) < 3) {
            susie_note <- "Could not align the uploaded LD matrices to the harmonised SNP set (row/column names must be SNP IDs matching the mapped SNP column, with >=3 SNPs in common) - multiple-signal colocalisation was not run."
          } else {
            h_ld <- h[h$SNP %in% common_ld, , drop = FALSE]
            d1s <- list(beta = h_ld$beta.exposure, varbeta = h_ld$se.exposure^2, N = d1$N, type = "quant",
                        snp = h_ld$SNP, LD = ld1[h_ld$SNP, h_ld$SNP, drop = FALSE])
            d2s <- list(beta = h_ld$beta.outcome, varbeta = h_ld$se.outcome^2, N = d2$N, type = vs$gwas_type,
                        snp = h_ld$SNP, LD = ld2[h_ld$SNP, h_ld$SNP, drop = FALSE])
            if (identical(vs$gwas_type, "cc")) d2s$s <- vs$case_frac
            susie_res <- tryCatch(
              suppressWarnings(suppressMessages(coloc::coloc.susie(
                dataset1 = d1s, dataset2 = d2s,
                susie.args = list(maxit = input$susie_maxit %||% 100, coverage = input$susie_coverage %||% 0.95),
                p1 = p1, p2 = p2, p12 = MCOL_DEFAULT_P12_SUSIE
              ))), error = function(e) e)
            if (inherits(susie_res, "error")) {
              susie_note <- sprintf("coloc.susie() failed: %s", conditionMessage(susie_res)); susie_res <- NULL
            } else if (is.null(susie_res$summary) || nrow(susie_res$summary) == 0 || all(is.na(susie_res$summary$nsnps))) {
              susie_note <- "coloc.susie() ran but found no credible set of putative causal variants in one or both datasets at the requested coverage - no multi-signal result to report."
              susie_res <- NULL
            }
          }
        }
      }

      region <- if (!is.null(h$chr.exposure) && any(!is.na(h$pos.exposure)))
        list(chr = h$chr.exposure[1], start = min(h$pos.exposure, na.rm = TRUE), end = max(h$pos.exposure, na.rm = TRUE))
      else NULL

      list(mode = "upload", cpg = vs$cpg, gwas_label = vs$gwas_label, gwas_type = vs$gwas_type,
           n_variants_pre_filter = nrow(h_full), n_variants = nrow(h), region = region,
           n_meth = d1$N, n_gwas = d2$N, priors = list(p1 = p1, p2 = p2, p12 = p12),
           pp_threshold = input$pp_threshold %||% MCOL_DEFAULT_PP_THRESHOLD,
           abf_res = abf_res, snp_df = snp_df, lead_variant = lead_variant,
           susie_res = susie_res, susie_note = susie_note,
           d1 = d1, d2 = d2, h = h, h_full = h_full)
    }

    observeEvent(input$run_btn, {
      rs <- tryCatch(
        if (identical((validate_state())$mode, "preloaded")) build_run_state_preloaded() else build_run_state_upload(),
        error = function(e) e
      )
      if (inherits(rs, "shiny.silent.error")) {
        showNotification(conditionMessage(rs), type = "warning", duration = 10)
        return()
      }
      if (inherits(rs, "error")) {
        showNotification(sprintf("Run Colocalisation failed: %s", conditionMessage(rs)), type = "error", duration = 10)
        return()
      }
      run_state(rs)
      stage_flags$run <- TRUE
      invalidate_from("plot")
      showNotification("Colocalisation complete - see the Results tab.", type = "message", duration = 6)
      updateTabsetPanel(session, "coloc_tabs", selected = "3. Results")

      if (!is.null(results)) {
        entry <- if (identical(rs$mode, "upload")) {
          pp <- rs$abf_res$summary
          list(mode = "upload", cpg = rs$cpg, gwas = rs$gwas_label, n_snp = rs$n_variants,
               pp_h4 = round(unname(pp["PP.H4.abf"]), 3), pp_h3 = round(unname(pp["PP.H3.abf"]), 3))
        } else {
          list(mode = "preloaded", cpgs_tested = rs$cpgs, n_cpgs = length(rs$cpgs))
        }
        existing <- results$coloc %||% list()
        existing[[if (identical(rs$mode, "upload")) rs$cpg else "preloaded"]] <- entry
        results$coloc <- existing
      }
    }, ignoreInit = TRUE)

    output$results_tab_body <- renderUI({
      if (!stage_flags$run) {
        return(p(class = "empty-note", icon("circle-info"),
                  "Not run yet. Complete Data & Setup and Filters & Parameters, then click \"Run Colocalisation\"."))
      }
      rs <- run_state()

      if (identical(rs$mode, "preloaded")) {
        tbl <- rs$table
        tbl$verdict <- mapply(.mcol_verdict, tbl$PP.H3, tbl$PP.H4)
        focus_row <- tbl[tbl$cpg == rs$focus_cpg, , drop = FALSE][1, ]
        tagList(
          box(
            width = NULL, title = "Colocalisation summary", status = "primary", solidHeader = FALSE,
            div(class = "methyl-stats-row", fluidRow(
              valueBox(rs$outcome_label, "GWAS trait", icon = icon("dna"), color = "light-blue", width = 3),
              valueBox(length(rs$cpgs), "CpGs tested", icon = icon("dna"), color = "purple", width = 3),
              valueBox(format(rs$n_variants, big.mark = ","), "Total SNPs across CpGs", icon = icon("code-branch"), color = "aqua", width = 3),
              valueBox("coloc.abf (preloaded)", "Analysis method", icon = icon("calculator"), color = "green", width = 3)
            )),
            p(class = "submodule-desc", "Genomic build: GRCh37. Priors: p1=1e-4, p2=1e-4, p12=1e-5 (coloc defaults, fixed by the upstream pipeline)."),
            selectInput(ns("pre_focus_select"), "Focus CpG (for the interpretation panel below)", choices = tbl$cpg, selected = rs$focus_cpg, width = "300px")
          ),
          uiOutput(ns("pre_focus_ui")),
          box(
            width = 12, title = "Per-CpG results", status = "primary", solidHeader = FALSE,
            p(class = "empty-note", icon("circle-info"),
              "SNP-level results aren't available for the Preloaded route - only the completed run's per-CpG PP.H0-H4 summary is bundled with this deployment."),
            div(class = "table-toolbar", downloadButton(ns("dl_hypotheses"), "Download CSV", class = "btn-sm")),
            DT::dataTableOutput(ns("pre_results_table"))
          )
        )
      } else {
        pp <- rs$abf_res$summary
        h0 <- unname(pp["PP.H0.abf"]); h1 <- unname(pp["PP.H1.abf"]); h2 <- unname(pp["PP.H2.abf"])
        h3 <- unname(pp["PP.H3.abf"]); h4 <- unname(pp["PP.H4.abf"])
        tagList(
          box(
            width = NULL, title = "Colocalisation summary", status = "primary", solidHeader = FALSE,
            div(class = "methyl-stats-row", fluidRow(
              valueBox(rs$gwas_label, "GWAS trait", icon = icon("dna"), color = "light-blue", width = 3),
              valueBox(rs$cpg, "Methylation phenotype (CpG)", icon = icon("dna"), color = "purple", width = 3),
              valueBox(if (!is.null(rs$region)) sprintf("chr%s:%s-%s", rs$region$chr, format(round(rs$region$start), big.mark = ","), format(round(rs$region$end), big.mark = ",")) else "n/a (no chr/pos mapped)",
                       "Genomic region", icon = icon("map"), color = "aqua", width = 3),
              valueBox(rs$n_variants, "Variants analysed", icon = icon("code-branch"), color = "green", width = 3)
            )),
            div(class = "methyl-stats-row", fluidRow(
              valueBox(unname(pp["nsnps"]), "Shared SNPs used", icon = icon("circle-check"), color = "yellow", width = 3),
              valueBox(if (!is.null(rs$susie_res)) "coloc.abf + coloc.susie" else "coloc.abf", "Analysis method", icon = icon("calculator"), color = "light-blue", width = 3),
              valueBox(sprintf("%s / %s", format(rs$n_meth, big.mark = ","), format(rs$n_gwas, big.mark = ",")), "Sample sizes (mQTL / GWAS)", icon = icon("users"), color = "purple", width = 3),
              valueBox(sprintf("%.0e / %.0e / %.0e", rs$priors$p1, rs$priors$p2, rs$priors$p12), "Priors (p1/p2/p12)", icon = icon("sliders"), color = "aqua", width = 3)
            )),
            p(strong("Lead variant (highest SNP.PP.H4): "), rs$lead_variant),
            tags$table(class = "table table-sm", tags$tbody(
              tags$tr(tags$td("H0 - no association with either trait"), tags$td(sprintf("%.4f", h0))),
              tags$tr(tags$td("H1 - association with methylation only"), tags$td(sprintf("%.4f", h1))),
              tags$tr(tags$td("H2 - association with GWAS trait only"), tags$td(sprintf("%.4f", h2))),
              tags$tr(tags$td("H3 - both associated, different causal variants"), tags$td(sprintf("%.4f", h3))),
              tags$tr(tags$td(strong("H4 - both associated, shared causal variant")), tags$td(strong(sprintf("%.4f", h4))))
            )),
            p(class = "empty-note", if (h4 >= rs$pp_threshold) icon("circle-check") else if (h3 >= rs$pp_threshold) icon("circle-exclamation") else icon("circle-info"),
              strong(" Verdict: "), .mcol_verdict(h3, h4, rs$pp_threshold))
          ),
          box(width = NULL, title = "Interpretation", status = "primary", solidHeader = FALSE, .mcol_interpret(h0, h1, h2, h3, h4)),
          if (!is.null(rs$susie_res)) box(
            width = NULL, title = "Multiple-signal colocalisation (coloc.susie)", status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "Each row is one detected pair of credible sets (one per dataset) - relaxes the single-causal-variant assumption coloc.abf makes."),
            DT::dataTableOutput(ns("susie_table"))
          ) else if (!is.null(rs$susie_note)) p(class = "empty-note", icon("circle-info"), rs$susie_note),
          box(
            width = 12, title = "SNP-level results", status = "primary", solidHeader = FALSE,
            div(class = "table-toolbar", downloadButton(ns("dl_snp"), "Download CSV", class = "btn-sm")),
            DT::dataTableOutput(ns("snp_table"))
          )
        )
      }
    })
    outputOptions(output, "results_tab_body", suspendWhenHidden = FALSE)

    output$pre_focus_ui <- renderUI({
      req(stage_flags$run, run_state())
      rs <- run_state()
      focus <- input$pre_focus_select %||% rs$focus_cpg
      row <- rs$table[rs$table$cpg == focus, , drop = FALSE][1, ]
      req(nrow(row) == 1 && !is.na(row$cpg))
      box(
        width = NULL, title = sprintf("Interpretation - %s", focus), status = "primary", solidHeader = FALSE,
        p(sprintf("%d GoDMC cis-mQTL SNPs tested against %s in this CpG's region.", row$nsnps, rs$outcome_label)),
        tags$table(class = "table table-sm", tags$tbody(
          tags$tr(tags$td("H0"), tags$td(sprintf("%.4f", row$PP.H0))), tags$tr(tags$td("H1"), tags$td(sprintf("%.4f", row$PP.H1))),
          tags$tr(tags$td("H2"), tags$td(sprintf("%.4f", row$PP.H2))), tags$tr(tags$td("H3"), tags$td(sprintf("%.4f", row$PP.H3))),
          tags$tr(tags$td(strong("H4")), tags$td(strong(sprintf("%.4f", row$PP.H4))))
        )),
        .mcol_interpret(row$PP.H0, row$PP.H1, row$PP.H2, row$PP.H3, row$PP.H4),
        p(class = "empty-note", icon("triangle-exclamation"),
          "GoDMC's cis-mQTL rows are a pre-filtered candidate list (each contributing cohort thresholded at p<1e-5 before meta-analysis), not an exhaustive dense scan of the region - coloc's power to distinguish PP.H3 from PP.H4 is bounded by this, a property of the source data, not of this analysis.")
      )
    })
    outputOptions(output, "pre_focus_ui", suspendWhenHidden = FALSE)

    output$pre_results_table <- DT::renderDataTable({
      req(stage_flags$run, run_state())
      rs <- run_state()
      tbl <- rs$table
      tbl$verdict <- mapply(.mcol_verdict, tbl$PP.H3, tbl$PP.H4)
      DT::datatable(tbl, rownames = FALSE, filter = "top", options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact") |>
        DT::formatRound(columns = c("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"), digits = 4)
    })
    outputOptions(output, "pre_results_table", suspendWhenHidden = FALSE)

    output$snp_table <- DT::renderDataTable({
      req(stage_flags$run, run_state())
      rs <- run_state(); req(identical(rs$mode, "upload"))
      DT::datatable(rs$snp_df, rownames = FALSE, filter = "top", options = list(pageLength = 12, scrollX = TRUE), class = "stripe hover compact") |>
        DT::formatSignif(columns = c("meth_pval", "gwas_pval", "log_abf", "snp_pp_h4"), digits = 4)
    })
    outputOptions(output, "snp_table", suspendWhenHidden = FALSE)

    output$susie_table <- DT::renderDataTable({
      req(stage_flags$run, run_state())
      rs <- run_state(); req(!is.null(rs$susie_res))
      DT::datatable(as.data.frame(rs$susie_res$summary), rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })
    outputOptions(output, "susie_table", suspendWhenHidden = FALSE)

    build_pp_bar_plot <- function(h0, h1, h2, h3, h4, title = NULL) {
      df <- data.frame(hypothesis = factor(c("H0", "H1", "H2", "H3", "H4"), levels = c("H0", "H1", "H2", "H3", "H4")),
                        pp = c(h0, h1, h2, h3, h4))
      ggplot(df, aes(x = hypothesis, y = pp, fill = hypothesis == "H4")) +
        geom_col() +
        scale_fill_manual(values = c(`TRUE` = ARTHOMIX_COLORS$red, `FALSE` = ARTHOMIX_COLORS$ink_muted), guide = "none") +
        labs(x = NULL, y = "Posterior probability", title = title) +
        theme_arthomix(12)
    }

    build_region_plot <- function(rs) {
      df <- rs$snp_df
      req(all(!is.na(df$pos)))
      meth_df <- data.frame(pos = df$pos, nlp = -log10(pmax(df$meth_pval, .Machine$double.xmin)), track = sprintf("Methylation/mQTL (%s)", rs$cpg))
      gwas_df <- data.frame(pos = df$pos, nlp = -log10(pmax(df$gwas_pval, .Machine$double.xmin)), track = rs$gwas_label)
      plot_df <- rbind(meth_df, gwas_df)
      cols <- stats::setNames(c(ARTHOMIX_COLORS$blue, ARTHOMIX_COLORS$red), c(unique(meth_df$track), unique(gwas_df$track)))
      ggplot(plot_df, aes(x = pos, y = nlp, color = track)) +
        geom_point(alpha = 0.7, size = 1.6) +
        geom_vline(xintercept = df$pos[df$snp == rs$lead_variant][1], linetype = "dashed", color = ARTHOMIX_COLORS$ink_muted) +
        facet_wrap(~track, ncol = 1, scales = "free_y") +
        scale_color_manual(values = cols, guide = "none") +
        labs(x = "Position (bp)", y = "-log10(p-value)", title = "Regional association") +
        theme_arthomix(12)
    }

    build_comparison_plot <- function(rs) {
      df <- rs$snp_df
      ggplot(df, aes(x = -log10(pmax(meth_pval, .Machine$double.xmin)), y = -log10(pmax(gwas_pval, .Machine$double.xmin)), color = snp_pp_h4)) +
        geom_point(alpha = 0.8, size = 2) +
        scale_color_gradient(low = ARTHOMIX_COLORS$grid, high = ARTHOMIX_COLORS$red, name = "SNP.PP.H4") +
        labs(x = "Methylation/mQTL -log10(p)", y = sprintf("%s -log10(p)", rs$gwas_label), title = "Comparison of association signals") +
        theme_arthomix(12)
    }

    build_posterior_plot <- function(rs) {
      df <- rs$snp_df
      if (all(!is.na(df$pos))) {
        ggplot(df, aes(x = pos, y = snp_pp_h4)) +
          geom_col(fill = ARTHOMIX_COLORS$blue) +
          geom_point(data = df[df$snp == rs$lead_variant, , drop = FALSE], aes(x = pos, y = snp_pp_h4), color = ARTHOMIX_COLORS$red, size = 3) +
          labs(x = "Position (bp)", y = "SNP-level posterior probability of a shared signal (SNP.PP.H4)", title = "Posterior support by variant") +
          theme_arthomix(12)
      } else {
        df2 <- df[order(-df$snp_pp_h4), , drop = FALSE][seq_len(min(30, nrow(df))), , drop = FALSE]
        df2$snp <- factor(df2$snp, levels = rev(df2$snp))
        ggplot(df2, aes(x = snp, y = snp_pp_h4)) +
          geom_col(fill = ARTHOMIX_COLORS$blue) +
          coord_flip() +
          labs(x = NULL, y = "SNP.PP.H4", title = "Posterior support by variant (top 30 - no position mapped)") +
          theme_arthomix(12)
      }
    }

    output$plots_tab_body <- renderUI({
      req(stage_flags$run)
      rs <- run_state()
      tagList(
        box(
          width = NULL, title = "Generate visualisation", status = "primary", solidHeader = FALSE,
          if (identical(rs$mode, "preloaded"))
            p(class = "submodule-desc", "Preloaded route: only the per-CpG posterior-probability summary can be plotted - no per-SNP GoDMC/RA-GWAS region data is bundled to build a regional or posterior plot from (see Data & Setup).")
          else p(class = "submodule-desc", "Builds the regional association, comparison, and posterior plots from the current Run Colocalisation result."),
          actionButton(ns("plot_btn"), "Generate Regional Plot", icon = icon("chart-area"), class = "btn-primary btn-sm")
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] > 0", ns("plot_btn")),
          if (identical(rs$mode, "preloaded")) tagList(
            box(width = 12, title = "PP.H0-H4 by CpG", status = "primary", solidHeader = FALSE,
                selectInput(ns("plot_pre_cpg"), "CpG", choices = rs$cpgs, selected = rs$cpgs[1], width = "260px"),
                withSpinner(plotOutput(ns("pre_pp_plot"), height = 320), color = "#2c6fbb", type = 6))
          ) else tagList(
            fluidRow(
              column(6, box(width = NULL, title = "Regional association", status = "primary", solidHeader = FALSE,
                            if (!is.null(rs$region))
                              withSpinner(plotOutput(ns("region_plot"), height = 420), color = "#2c6fbb", type = 6)
                            else p(class = "empty-note", icon("circle-info"),
                              "Regional plot needs a SNP position - map a \"SNP position\" column for the methylation/mQTL file in Data & Setup to enable it. No position was inferred."))),
              column(6, box(width = NULL, title = "Comparison of signals", status = "primary", solidHeader = FALSE,
                            withSpinner(plotOutput(ns("comparison_plot"), height = 420), color = "#2c6fbb", type = 6)))
            ),
            box(width = 12, title = "Posterior support (per-SNP shared-signal probability)", status = "primary", solidHeader = FALSE,
                withSpinner(plotOutput(ns("posterior_plot"), height = 360), color = "#2c6fbb", type = 6))
          )
        )
      )
    })
    outputOptions(output, "plots_tab_body", suspendWhenHidden = FALSE)

    observeEvent(input$plot_btn, { stage_flags$plot <- TRUE }, ignoreInit = TRUE)

    output$pre_pp_plot <- renderPlot({
      req(input$plot_btn > 0, run_state())
      rs <- run_state(); row <- rs$table[rs$table$cpg == (input$plot_pre_cpg %||% rs$cpgs[1]), , drop = FALSE][1, ]
      req(nrow(row) == 1)
      build_pp_bar_plot(row$PP.H0, row$PP.H1, row$PP.H2, row$PP.H3, row$PP.H4, title = row$cpg)
    })
    output$region_plot <- renderPlot({ req(input$plot_btn > 0, run_state(), (run_state())$region); build_region_plot(run_state()) })
    output$comparison_plot <- renderPlot({ req(input$plot_btn > 0, run_state()); build_comparison_plot(run_state()) })
    output$posterior_plot <- renderPlot({ req(input$plot_btn > 0, run_state()); build_posterior_plot(run_state()) })

    sensitivity_state <- reactiveVal(NULL)

    .mcol_rerun_with <- function(h_full, gwas_type, case_frac, p1, p2, p12,
                                  pval_meth = 1, pval_gwas = 1, maf_min = 0, maf_max = 0.5,
                                  window_kb = NULL, window_center = NULL, min_n = 0) {
      h <- h_full
      if (!is.null(window_kb) && !is.null(window_center) && "pos.exposure" %in% colnames(h)) {
        h <- h[!is.na(h$pos.exposure) & abs(h$pos.exposure - window_center) <= window_kb * 1000, , drop = FALSE]
      }
      h <- h[h$pval.exposure <= pval_meth & h$pval.outcome <= pval_gwas, , drop = FALSE]
      maf_val <- if (!is.null(h$eaf.exposure)) pmin(h$eaf.exposure, 1 - h$eaf.exposure) else rep(NA_real_, nrow(h))
      h <- h[is.na(maf_val) | (maf_val >= maf_min & maf_val <= maf_max), , drop = FALSE]
      h <- h[(is.na(h$samplesize.exposure) | h$samplesize.exposure >= min_n) &
             (is.na(h$samplesize.outcome) | h$samplesize.outcome >= min_n), , drop = FALSE]
      if (nrow(h) < 6) return(NULL)
      d1 <- list(beta = h$beta.exposure, varbeta = h$se.exposure^2, N = round(stats::median(h$samplesize.exposure, na.rm = TRUE)), type = "quant", snp = h$SNP)
      if (!is.null(h$eaf.exposure)) d1$MAF <- pmin(h$eaf.exposure, 1 - h$eaf.exposure)
      if (is.null(d1$MAF)) return(NULL)
      d2 <- list(beta = h$beta.outcome, varbeta = h$se.outcome^2, N = round(stats::median(h$samplesize.outcome, na.rm = TRUE)), type = gwas_type, snp = h$SNP)
      if (identical(gwas_type, "cc")) d2$s <- case_frac
      res <- tryCatch(suppressWarnings(coloc::coloc.abf(dataset1 = d1, dataset2 = d2, p1 = p1, p2 = p2, p12 = p12)), error = function(e) NULL)
      if (is.null(res)) return(NULL)
      c(nsnps = nrow(h), PP.H3 = unname(res$summary["PP.H3.abf"]), PP.H4 = unname(res$summary["PP.H4.abf"]))
    }

    output$sensitivity_tab_body <- renderUI({
      req(stage_flags$run)
      rs <- run_state()
      if (identical(rs$mode, "preloaded")) {
        return(p(class = "empty-note", icon("circle-info"),
                  "Sensitivity re-analysis needs per-SNP Bayes factors, which aren't bundled with the preloaded results - only summary posterior probabilities are available for this route. Use Upload Data with your own summary statistics to run sensitivity analysis."))
      }
      tagList(
        box(
          width = NULL, title = "Prior sensitivity", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Re-weights the current result's Bayes factors across a grid of p12 values (coloc::sensitivity()) to test how robust the H3/H4 conclusion is to the chosen prior."),
          fluidRow(
            column(6, numericInput(ns("sens_threshold"), "Decision rule: H4 above", value = MCOL_DEFAULT_PP_THRESHOLD, min = 0.1, max = 0.99, step = 0.05, width = "100%")),
            column(6, numericInput(ns("sens_npoints"), "Grid points", value = 100, min = 10, max = 1000, step = 10, width = "100%"))
          ),
          actionButton(ns("sens_prior_btn"), "Run Sensitivity Analysis", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] > 0", ns("sens_prior_btn")),
          box(width = 12, title = "PP.H3 / PP.H4 across the prior grid", status = "primary", solidHeader = FALSE,
              withSpinner(plotOutput(ns("sens_prior_plot"), height = 340), color = "#2c6fbb", type = 6),
              div(class = "table-toolbar", downloadButton(ns("dl_sensitivity"), "Download CSV", class = "btn-sm")),
              DT::dataTableOutput(ns("sens_prior_table")))
        ),
        box(
          width = NULL, title = "Parameter sensitivity", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Re-runs coloc.abf() varying one filter at a time from the current Run Colocalisation settings (genomic window, p-value thresholds, MAF, minimum sample size), holding the others fixed, to show how many shared SNPs and PP.H3/PP.H4 change."),
          actionButton(ns("sens_param_btn"), "Run Sensitivity Analysis", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] > 0", ns("sens_param_btn")),
          box(width = 12, title = "Parameter sensitivity table", status = "primary", solidHeader = FALSE, DT::dataTableOutput(ns("sens_param_table")))
        )
      )
    })
    outputOptions(output, "sensitivity_tab_body", suspendWhenHidden = FALSE)

    observeEvent(input$sens_prior_btn, {
      rs <- run_state(); req(rs, identical(rs$mode, "upload"))
      sens <- tryCatch(
        suppressWarnings(suppressMessages(coloc::sensitivity(
          rs$abf_res, rule = sprintf("H4 > %g", input$sens_threshold %||% MCOL_DEFAULT_PP_THRESHOLD),
          npoints = input$sens_npoints %||% 100, doplot = FALSE
        ))),
        error = function(e) e
      )
      if (inherits(sens, "error")) { sens <- NULL }
      st <- sensitivity_state() %||% list()
      st$prior <- sens
      sensitivity_state(st)
      stage_flags$sensitivity <- TRUE
    }, ignoreInit = TRUE)

    observeEvent(input$sens_param_btn, {
      rs <- run_state(); req(rs, identical(rs$mode, "upload"))
      base_window_kb <- if (isTRUE(input$use_window)) (input$window_kb %||% MCOL_DEFAULT_WINDOW_KB) else NULL
      center <- if (!is.null(base_window_kb) && "pos.exposure" %in% colnames(rs$h_full))
        rs$h_full$pos.exposure[which.min(rs$h_full$pval.exposure)] else NULL
      base <- list(pval_meth = input$f_pval_meth %||% 1, pval_gwas = input$f_pval_gwas %||% 1,
                   maf_min = input$f_maf_min %||% 0, maf_max = input$f_maf_max %||% 0.5,
                   window_kb = base_window_kb, window_center = center, min_n = input$f_min_n %||% 0)
      p1 <- rs$priors$p1; p2 <- rs$priors$p2; p12 <- rs$priors$p12

      rows <- list()
      add_row <- function(label, value, ...) {
        args <- utils::modifyList(base, list(...))
        out <- do.call(.mcol_rerun_with, c(list(h_full = rs$h_full, gwas_type = rs$gwas_type, case_frac = (validate_state())$case_frac %||% 0.33, p1 = p1, p2 = p2, p12 = p12), args))
        if (!is.null(out)) rows[[length(rows) + 1]] <<- data.frame(parameter = label, value = value, nsnps = out["nsnps"], PP.H3 = out["PP.H3"], PP.H4 = out["PP.H4"], row.names = NULL)
      }
      add_row("baseline", "current settings")
      if (!is.null(base_window_kb)) {
        for (mult in c(0.5, 2, 4)) add_row("window (kb)", base_window_kb * mult, window_kb = base_window_kb * mult)
      }
      for (mv in unique(c(base$maf_min, 0, 0.01, 0.05))) if (mv != base$maf_min) add_row("min MAF", mv, maf_min = mv)
      for (pv in unique(c(base$pval_meth, 1, 0.05, 5e-8))) if (pv != base$pval_meth) add_row("mQTL p-value threshold", pv, pval_meth = pv)

      tbl <- if (length(rows) > 0) do.call(rbind, rows) else NULL
      st <- sensitivity_state() %||% list()
      st$param <- tbl
      sensitivity_state(st)
      stage_flags$sensitivity <- TRUE
    }, ignoreInit = TRUE)

    output$sens_prior_plot <- renderPlot({
      req(input$sens_prior_btn > 0)
      st <- sensitivity_state(); req(st, !is.null(st$prior))
      df <- st$prior
      validate(need(all(c("p12", "PP.H3.abf", "PP.H4.abf") %in% colnames(df)), "Unexpected output shape from coloc::sensitivity()."))
      long <- rbind(
        data.frame(p12 = df$p12, pp = df$PP.H3.abf, hypothesis = "H3"),
        data.frame(p12 = df$p12, pp = df$PP.H4.abf, hypothesis = "H4")
      )
      ggplot(long, aes(x = p12, y = pp, color = hypothesis)) +
        geom_line(linewidth = 0.9) +
        scale_x_log10() +
        scale_color_manual(values = stats::setNames(c(ARTHOMIX_COLORS$ink_muted, ARTHOMIX_COLORS$red), c("H3", "H4"))) +
        labs(x = "p12 (log scale)", y = "Posterior probability", color = NULL, title = "Robustness of H3/H4 to the shared-association prior") +
        theme_arthomix(12)
    })

    output$sens_prior_table <- DT::renderDataTable({
      req(input$sens_prior_btn > 0)
      st <- sensitivity_state(); req(st, !is.null(st$prior))
      DT::datatable(st$prior, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })

    output$sens_param_table <- DT::renderDataTable({
      req(input$sens_param_btn > 0)
      st <- sensitivity_state(); req(st, !is.null(st$param))
      DT::datatable(st$param, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact") |>
        DT::formatSignif(columns = c("PP.H3", "PP.H4"), digits = 4)
    })

    output$export_tab_body <- renderUI({
      if (!stage_flags$run) {
        return(p(class = "empty-note", icon("circle-info"), "Export becomes available once Run Colocalisation has completed successfully."))
      }
      rs <- run_state()
      tagList(
        box(
          width = NULL, title = "Export results", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Every download below reflects the most recent Run Colocalisation result."),
          fluidRow(
            column(4, downloadButton(ns("dl_hypotheses"), "Hypothesis probabilities (CSV)", class = "btn-sm", style = "width:100%; margin-bottom:8px;")),
            column(4, downloadButton(ns("dl_validation"), "Validation/harmonisation summary (CSV)", class = "btn-sm", style = "width:100%; margin-bottom:8px;")),
            column(4, downloadButton(ns("dl_params"), "Analysis parameters (CSV)", class = "btn-sm", style = "width:100%; margin-bottom:8px;"))
          ),
          if (identical(rs$mode, "upload")) fluidRow(
            column(4, downloadButton(ns("dl_snp"), "SNP-level results (CSV)", class = "btn-sm", style = "width:100%; margin-bottom:8px;")),
            column(4, downloadButton(ns("dl_shared"), "Shared/high-posterior variants (CSV)", class = "btn-sm", style = "width:100%; margin-bottom:8px;")),
            column(4, downloadButton(ns("dl_sensitivity"), "Sensitivity-analysis results (CSV)", class = "btn-sm", style = "width:100%; margin-bottom:8px;"))
          ),
          if (identical(rs$mode, "upload") && stage_flags$plot) tagList(
            selectInput(ns("plot_format"), "Plot format", choices = c(PNG = "png", PDF = "pdf", SVG = "svg"), selected = "png", width = "160px"),
            fluidRow(
              if (!is.null(rs$region)) column(4, downloadButton(ns("dl_plot_region"), "Regional plot", class = "btn-sm", style = "width:100%;")),
              column(4, downloadButton(ns("dl_plot_comparison"), "Comparison plot", class = "btn-sm", style = "width:100%;")),
              column(4, downloadButton(ns("dl_plot_posterior"), "Posterior plot", class = "btn-sm", style = "width:100%;"))
            )
          )
        )
      )
    })
    outputOptions(output, "export_tab_body", suspendWhenHidden = FALSE)

    output$dl_hypotheses <- downloadHandler(
      filename = function() "coloc_hypothesis_probabilities.csv",
      content = function(file) {
        rs <- run_state()
        df <- if (identical(rs$mode, "preloaded")) {
          tbl <- rs$table; tbl$verdict <- mapply(.mcol_verdict, tbl$PP.H3, tbl$PP.H4); tbl
        } else {
          pp <- rs$abf_res$summary
          data.frame(cpg = rs$cpg, gwas = rs$gwas_label, as.list(pp),
                     verdict = .mcol_verdict(unname(pp["PP.H3.abf"]), unname(pp["PP.H4.abf"]), rs$pp_threshold))
        }
        write.csv(df, file, row.names = FALSE)
      }
    )
    output$dl_snp <- downloadHandler(
      filename = function() sprintf("coloc_snp_results_%s.csv", (run_state())$cpg %||% "result"),
      content = function(file) write.csv((run_state())$snp_df, file, row.names = FALSE)
    )
    output$dl_shared <- downloadHandler(
      filename = function() sprintf("coloc_shared_variants_%s.csv", (run_state())$cpg %||% "result"),
      content = function(file) {
        rs <- run_state(); df <- rs$snp_df[rs$snp_df$snp_pp_h4 >= 0.01, , drop = FALSE]
        write.csv(df, file, row.names = FALSE)
      }
    )
    output$dl_sensitivity <- downloadHandler(
      filename = function() "coloc_sensitivity_analysis.csv",
      content = function(file) {
        st <- sensitivity_state()
        df <- if (!is.null(st) && !is.null(st$prior)) st$prior else data.frame()
        write.csv(df, file, row.names = FALSE)
      }
    )
    output$dl_params <- downloadHandler(
      filename = function() "coloc_analysis_parameters.csv",
      content = function(file) {
        rs <- run_state()
        df <- if (identical(rs$mode, "preloaded")) {
          data.frame(parameter = c("mode", "p1", "p2", "p12", "window_bp", "method"),
                     value = c("preloaded", MCOL_DEFAULT_P1, MCOL_DEFAULT_P2, MCOL_DEFAULT_P12, 1e6, "coloc.abf"))
        } else {
          data.frame(parameter = c("mode", "cpg", "gwas_trait", "gwas_type", "n_variants", "n_meth", "n_gwas",
                                     "p1", "p2", "p12", "pp_threshold", "method"),
                     value = c("upload", rs$cpg, rs$gwas_label, rs$gwas_type, rs$n_variants, rs$n_meth, rs$n_gwas,
                               rs$priors$p1, rs$priors$p2, rs$priors$p12, rs$pp_threshold,
                               if (!is.null(rs$susie_res)) "coloc.abf + coloc.susie" else "coloc.abf"))
        }
        write.csv(df, file, row.names = FALSE)
      }
    )
    output$dl_validation <- downloadHandler(
      filename = function() "coloc_validation_summary.csv",
      content = function(file) {
        s <- (validate_state())$summary
        write.csv(data.frame(metric = names(s), value = unlist(lapply(s, function(x) if (is.null(x)) NA else x))), file, row.names = FALSE)
      }
    )

    make_plot_dl <- function(build_fn) downloadHandler(
      filename = function() sprintf("coloc_plot.%s", input$plot_format %||% "png"),
      content = function(file) ggplot2::ggsave(file, plot = build_fn(run_state()), width = 9, height = 6, dpi = 300, device = input$plot_format %||% "png")
    )
    output$dl_plot_region <- make_plot_dl(build_region_plot)
    output$dl_plot_comparison <- make_plot_dl(build_comparison_plot)
    output$dl_plot_posterior <- make_plot_dl(build_posterior_plot)
  })
}
