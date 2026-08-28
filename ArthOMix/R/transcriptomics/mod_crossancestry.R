## R/mod_crossancestry.R
## Submodule: Cross-Ancestry Validation (Section 2.12)
## "Your analysis" recomputes, on a per-sex Run click, which prioritised
## genes replicate in an independent European RA GWAS (Stahl et al. 2010)
## and transfer to an East Asian RA GWAS (BioBank Japan, Ishigaki et al.
## 2020), at whatever significance threshold and direction-of-effect rule
## the user sets - instead of the fixed 0.05/direction-required columns
## bundled in MR35_crossancestry_{female,male}.csv. Same house style as
## Cross-Tissue Validation (mod_crosstissue.R): a headline "validated
## biomarker" count driven by one shared classification rule, a concordance
## scatter, a ranked evidence plot, and a per-gene table - adapted to what
## this module's OWN data actually is (per-SNP MR odds ratios across three
## GWAS cohorts, not expression samples to fit a classifier on). Adjusting a
## slider/radio does NOT recompute anything until Run is clicked again -
## same click-to-lock-in-results convention as Diagnostic Model and
## Cross-Tissue Validation, even though this module's own computation is
## cheap enough that it could have stayed fully live.
##
## WHAT "VALIDATED CROSS-ANCESTRY BIOMARKER" MEANS HERE: a gene whose MR
## estimate (a) replicates in the second European cohort (same direction as
## the European discovery GWAS, nominal p < threshold) AND (b) transfers to
## the East Asian cohort (same direction, nominal p < threshold) - i.e. this
## module's "shared EUR+EAS" class. This is genetic evidence for a causal
## effect on RA liability, not a re-evaluated diagnostic classifier: no
## expression data from any non-European cohort exists, so the panel itself
## is never re-fit here (Section 2.12.1).
##
## WHY THE EAST ASIAN ARM ISN'T A LIKE-FOR-LIKE REPLICATION: the genetic
## instruments are European cis-eQTL variants - no comparable East Asian
## whole-blood eQTL resource exists - so that arm tests European-derived
## instruments against an East Asian outcome GWAS (Section 2.12.7). A gene
## failing to transfer may reflect the instrument tagging the causal variant
## poorly in East Asians, not an absence of a real effect - this is why
## "untestable in East Asian" (no instrument variant survived harmonisation)
## is its own KPI tile and table state, never folded into "not transferable".
##
## NOT SEX-SPECIFIC ESTIMATES: both the eQTL exposure and all three outcome
## GWAS are sex-combined (Section 2.12.1) - the Female/Male tabs here
## partition genes by which sex-stratified panel prioritised them, not by a
## different MR model per sex, so a gene appearing in both tabs carries an
## identical estimate in both.

mod_crossancestry_config <- list(
  id = "crossancestry", group = "Validation",
  title = "Cross-Ancestry Validation",
  description = "Validation of the diagnostic model based on different ancestral type and sex. Either work on preloaded data or upload your own data.",
  icon = "earth-americas"
)

## Single source of truth for "did this gene replicate in Europeans and
## transfer to East Asians, at the user's current thresholds" - every KPI
## tile, plot and table column reads from these three columns so they can
## never disagree with one another.
ca_classify <- function(df, p_eur, p_eas, require_dir) {
  dir_eur <- if (require_dir) df$dir_okada_eq_stahl else TRUE
  dir_eas <- if (require_dir) df$dir_okada_eq_bbj else TRUE
  ## A gene with no computable estimate in a cohort (e.g. zero surviving
  ## instrument SNPs) yields NA here, not FALSE - coerced explicitly so it
  ## reads as "not replicated"/"not transferable" rather than silently
  ## poisoning every downstream sum()/table with NA.
  replicated <- dir_eur & df$p_stahl < p_eur
  transferable <- df$testable_EAS & dir_eas & df$p_bbj < p_eas
  df$replicated_EUR <- !is.na(replicated) & replicated
  df$transferable_EAS <- !is.na(transferable) & transferable
  df$biomarker <- df$replicated_EUR & df$transferable_EAS
  df$ancestry_class <- ifelse(
    df$biomarker, "shared EUR+EAS",
    ifelse(df$replicated_EUR & !df$testable_EAS, "EUR-replicated, untestable in EAS",
           ifelse(df$replicated_EUR & !df$transferable_EAS, "EUR-replicated, not EAS",
                  "not EUR-replicated"))
  )
  df
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------

mod_crossancestry_sex_panel <- function(ns, sex_label) {
  run_id <- paste0("run_", sex_label, "_btn")
  sex_title <- tools::toTitleCase(sex_label)
  tagList(
    actionButton(ns(run_id), paste("Run", sex_title), icon = icon("play"), class = "btn-primary btn-sm"),
    div(style = "height:10px;"),
    conditionalPanel(
      condition = sprintf("input['%s'] > 0", ns(run_id)),
      withSpinner(uiOutput(ns(paste0(sex_label, "_stats"))), color = "#2c6fbb", type = 6),
      fluidRow(
        column(6, withSpinner(plotOutput(ns(paste0(sex_label, "_concordance_plot")), height = 380), color = "#2c6fbb", type = 6)),
        column(6, withSpinner(plotOutput(ns(paste0(sex_label, "_rank_plot")), height = "auto"), color = "#2c6fbb", type = 6))
      ),
      div(class = "table-toolbar", downloadButton(ns(paste0(sex_label, "_download")), "Per-gene table (CSV)", class = "btn-sm")),
      DT::dataTableOutput(ns(paste0(sex_label, "_table")))
    )
  )
}

mod_crossancestry_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$style(HTML(".ca-module .small-box h3 { font-size: 20px; white-space: normal; overflow-wrap: break-word; line-height: 1.2; }")),
    div(class = "ca-module",
    fluidRow(
      column(
        3,
        box(
          width = NULL, title = "Ancestry GWAS", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Applies to both Female and Male tabs (sex-combined GWAS)."),
          radioButtons(
            ns("data_source"), NULL,
            choiceNames = list(
              tagList(icon("database"), " Bundled Stahl (EUR) + BioBank Japan (EAS) (default)"),
              tagList(icon("upload"), " Upload your own replication/transfer GWAS")
            ),
            choiceValues = list("project", "upload"), selected = "project"
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'upload'", ns("data_source")),
            selectInput(
              ns("arm_replaced"), "Which arm does this replace?",
              choices = c("European replication (like Stahl 2010)" = "eur", "Transfer ancestry (like BioBank Japan)" = "eas"),
              selected = "eas", width = "100%"
            ),
            p(class = "submodule-desc",
              "A delimited file (CSV/TSV), one row per SNP: any second GWAS, tested against this project's own cis-eQTL instruments (the same ones the European/Okada discovery arm uses) for whichever genes are in the panel below. SNP IDs must be rsIDs to match those instruments."),
            textInput(ns("gwas_label"), "GWAS label (for display only)", value = "Uploaded GWAS", width = "100%"),
            fileInput(ns("gwas_file"), "GWAS file", accept = c(".csv", ".tsv", ".txt")),
            uiOutput(ns("gwas_map_ui"))
          )
        ),
        box(
          width = NULL, title = "Significance thresholds", status = "primary", solidHeader = FALSE,
          sliderInput(ns("p_eur"), "European replication (Stahl et al.)", value = 0.05, min = 0.001, max = 0.2, step = 0.001),
          sliderInput(ns("p_eas"), "East Asian transfer (BioBank Japan)", value = 0.05, min = 0.001, max = 0.2, step = 0.001),
          radioButtons(
            ns("require_dir"), "Direction of effect",
            choiceNames = list("Must agree with European discovery (recommended)", "Ignore direction (p-value only)"),
            choiceValues = list("yes", "no"), selected = "yes"
          ),
          div(class = "empty-note", style = "font-size: 12.5px;", icon("circle-info"),
              "Adjust, then click Run (Female/Male tab) - recomputed from per-SNP MR estimates, not the fixed 0.05 in the bundled columns.")
        )
      ),
      column(
        9,
        tabsetPanel(
          id = ns("sex_tabs"), type = "tabs",
          tabPanel("Female", br(), mod_crossancestry_sex_panel(ns, "female")),
          tabPanel("Male", br(), mod_crossancestry_sex_panel(ns, "male"))
        )
      )
    ),
    uiOutput(ns("references_box_ui"))
    )
  )
}

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------

mod_crossancestry_server <- function(id, dataset, results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    tables <- list(
      female = read_table_safe("MR35_crossancestry_female.csv"),
      male   = read_table_safe("MR35_crossancestry_male.csv")
    )

    ## ---------------------------------------------------------------------
    ## Upload-a-replacement-GWAS mode: parse + map the uploaded file with
    ## the SAME shared helpers mod_mr.R's and mod_coloc.R's own upload modes
    ## use (global.R). Loaded lazily and cached: MR_primary_objects.rds is
    ## the same file mod_mr.R already loads, cheap either way, but there's
    ## no reason to pay for it when this module's upload mode is never used
    ## this session.
    ## ---------------------------------------------------------------------
    gwas_df_r <- reactive({ req(input$gwas_file); read_uploaded_table(input$gwas_file$datapath) })
    output$gwas_map_ui <- gwas_col_map_ui(ns, reactive(input$gwas_file), gwas_df_r, "gwas", "GWAS file")

    mr_instrument_cache <- NULL
    get_mr_instruments <- function() {
      if (is.null(mr_instrument_cache)) mr_instrument_cache <<- load_mr_instrument_table()
      mr_instrument_cache
    }

    ## Formatted once per upload (Shiny's own reactive cache), reused by
    ## whichever sex's Run is clicked - the uploaded GWAS itself doesn't
    ## depend on sex.
    uploaded_outcome <- reactive({
      req(input$gwas_file, input$gwas_snp, input$gwas_beta, input$gwas_se, input$gwas_pval, input$gwas_ea, input$gwas_oa)
      raw <- gwas_df_r()
      validate(need(!is.null(raw), "Could not read the uploaded GWAS file."))
      label <- if (nzchar(trimws(input$gwas_label %||% ""))) trimws(input$gwas_label) else "Uploaded GWAS"
      out_fmt <- TwoSampleMR::format_data(
        raw, type = "outcome", snp_col = input$gwas_snp, beta_col = input$gwas_beta,
        se_col = input$gwas_se, pval_col = input$gwas_pval, effect_allele_col = input$gwas_ea,
        other_allele_col = input$gwas_oa, eaf_col = if (nzchar(input$gwas_eaf %||% "")) input$gwas_eaf else "eaf"
      )
      out_fmt$outcome <- label
      list(dat = out_fmt, label = label)
    })

    ## Re-harmonises this project's own Okada-discovery cis-eQTL instruments
    ## (the SAME exposure-side rows mod_mr.R's bundled path reads, via the
    ## shared global.R loader) against the uploaded outcome instead of
    ## Okada, restricted to whichever genes are in the current sex's panel -
    ## the exact TwoSampleMR harmonise_data() + estimate_mr_set() machinery
    ## mod_mr.R's own upload/batch modes already use, applied here to a
    ## small, fixed gene list instead of one gene or all ~1,700.
    live_arm_for_genes <- function(genes) {
      out <- uploaded_outcome()
      instr_all <- get_mr_instruments()$dat
      exp_cols <- c("SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure",
                    "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure",
                    "chr.exposure", "pos.exposure", "id.exposure", "exposure",
                    "mr_keep.exposure", "pval_origin.exposure", "data_source.exposure")
      instr <- instr_all[instr_all$gene %in% genes & instr_all$pval.exposure <= 5e-8 & instr_all$Fstat >= 10, , drop = FALSE]
      if (nrow(instr) == 0) return(list(table = NULL, n_input_genes = length(unique(genes)), n_output_genes = 0L, label = out$label))

      dat <- tryCatch(TwoSampleMR::harmonise_data(instr[, exp_cols, drop = FALSE], out$dat, action = 2), error = function(e) NULL)
      if (is.null(dat) || nrow(dat) == 0) return(list(table = NULL, n_input_genes = length(unique(instr$gene)), n_output_genes = 0L, label = out$label))
      dat <- dat[dat$mr_keep, , drop = FALSE]
      if (nrow(dat) == 0) return(list(table = NULL, n_input_genes = length(unique(instr$gene)), n_output_genes = 0L, label = out$label))
      ## Re-derive gene from (SNP, id.exposure) rather than trust a
      ## passed-through column - the same defensive pattern global.R's
      ## load_mr_instrument_table() uses for the SAME reason: harmonise_data
      ## isn't guaranteed to preserve columns it doesn't itself recognise.
      dat$gene <- instr$gene[match(paste(dat$SNP, dat$id.exposure), paste(instr$SNP, instr$id.exposure))]

      rows <- lapply(split(dat, dat$gene), function(dg) {
        est <- tryCatch(estimate_mr_set(dg, full = FALSE), error = function(e) NULL)
        if (is.null(est)) return(NULL)
        prim <- est$res_table[est$res_table$primary, , drop = FALSE][1, ]
        data.frame(gene = dg$gene[1], OR = exp(prim$estimate), p = prim$p, nSNP = est$n_snp, stringsAsFactors = FALSE)
      })
      tbl <- do.call(rbind, rows)
      list(table = tbl, n_input_genes = length(unique(instr$gene)),
           n_output_genes = if (is.null(tbl)) 0L else nrow(tbl), label = out$label)
    }

    result_for <- function(sex_label) {
      df <- tables[[sex_label]]
      validate(need(!is.null(df), "The cross-ancestry results table was not found on disk."))
      ## Every submodule server is instantiated at app startup regardless of
      ## whether its UI has been added yet (server.R inserts UI lazily on
      ## demand), so these inputs can still be NULL the first time this runs
      ## - fall back to the sliders'/radio's own UI defaults rather than
      ## erroring, matching mod_crosstissue.R's ct_advanced_params() pattern.
      p_eur <- input$p_eur %||% 0.05
      p_eas <- input$p_eas %||% 0.05
      require_dir <- !identical(input$require_dir %||% "yes", "no")

      live_note <- NULL
      if (identical(input$data_source %||% "project", "upload")) {
        arm <- input$arm_replaced %||% "eas"
        live <- live_arm_for_genes(df$gene)
        validate(need(!is.null(live$table) && nrow(live$table) > 0, sprintf(
          "None of the %d genes in this panel produced a usable harmonised instrument against your uploaded GWAS - check that SNP IDs are rsIDs and allele columns are mapped correctly.",
          live$n_input_genes
        )))
        m <- match(df$gene, live$table$gene)
        if (identical(arm, "eur")) {
          df$OR_stahl <- live$table$OR[m]
          df$p_stahl  <- live$table$p[m]
          df$nSNP_stahl <- live$table$nSNP[m]
          df$dir_okada_eq_stahl <- (df$OR_okada > 1) == (df$OR_stahl > 1)
        } else {
          df$OR_bbj <- live$table$OR[m]
          df$p_bbj  <- live$table$p[m]
          df$nSNP_bbj <- live$table$nSNP[m]
          df$dir_okada_eq_bbj <- (df$OR_okada > 1) == (df$OR_bbj > 1)
          df$testable_EAS <- !is.na(m)
          df$eafgap_bbj <- NA_real_  ## BBJ-specific diagnostic, meaningless for an arbitrary uploaded GWAS
        }
        live_note <- list(arm = arm, label = live$label, n_input_genes = live$n_input_genes, n_output_genes = live$n_output_genes)
      }

      list(df = ca_classify(df, p_eur, p_eas, require_dir), live_note = live_note)
    }
    ## Locked to whatever the sliders/radio/upload said at the moment Run was
    ## last clicked, for that sex - changing a threshold, or re-uploading,
    ## afterwards does nothing until Run is clicked again (ignoreInit = TRUE:
    ## skip the phantom "click" every actionButton starts at value 0, not a
    ## real click).
    result_female <- eventReactive(input$run_female_btn, result_for("female"), ignoreInit = TRUE)
    result_male <- eventReactive(input$run_male_btn, result_for("male"), ignoreInit = TRUE)
    res_sex <- function(sex_label) if (identical(sex_label, "female")) result_female else result_male

    ## Saved into results$crossancestry as each sex's Run completes - same
    ## results$<id> convention every other submodule uses, so ArthOChat sees
    ## this module's headline numbers as soon as they exist.
    save_result <- function(sex_label, res) {
      df <- res$df
      results$crossancestry <- utils::modifyList(
        results$crossancestry %||% list(),
        setNames(list(list(
          n_tested = nrow(df), n_testable_EAS = sum(df$testable_EAS, na.rm = TRUE),
          n_biomarkers = sum(df$biomarker, na.rm = TRUE), biomarker_genes = df$gene[df$biomarker],
          uploaded = !is.null(res$live_note), arm_replaced = res$live_note$arm %||% NA_character_
        )), sex_label)
      )
    }
    observeEvent(result_female(), save_result("female", result_female()))
    observeEvent(result_male(), save_result("male", result_male()))

    register_sex_outputs <- function(sex_label, res) {
      output[[paste0(sex_label, "_stats")]] <- renderUI({
        rr <- res(); df <- rr$df
        n_total <- nrow(df); n_testable <- sum(df$testable_EAS, na.rm = TRUE)
        n_eur <- sum(df$replicated_EUR, na.rm = TRUE)
        n_eas <- sum(df$transferable_EAS, na.rm = TRUE)
        n_bio <- sum(df$biomarker, na.rm = TRUE)
        n_untestable <- sum(!df$testable_EAS, na.rm = TRUE)
        tagList(
          if (!is.null(rr$live_note)) p(class = "empty-note", icon("upload"), sprintf(
            "%s arm replaced with your uploaded GWAS (%s) - %d of %d panel genes produced a usable harmonised instrument.",
            if (identical(rr$live_note$arm, "eur")) "European replication" else "Transfer ancestry",
            rr$live_note$label, rr$live_note$n_output_genes, rr$live_note$n_input_genes
          )),
          fluidRow(
            valueBox(n_bio, "Validated cross-ancestry biomarkers", icon = icon("award"), color = "green", width = 3),
            valueBox(sprintf("%d / %d", n_eur, n_total), "Replicated in Europeans", icon = icon("check-double"), color = "light-blue", width = 3),
            valueBox(sprintf("%d / %d", n_eas, n_testable), "Transferable to East Asian", icon = icon("earth-asia"), color = "purple", width = 3),
            valueBox(n_untestable, "Untestable in East Asian", icon = icon("triangle-exclamation"), color = "yellow", width = 3)
          )
        )
      })

      ## Ancestry concordance: European replication OR vs East Asian transfer
      ## OR (log2 scale, so OR=1 sits at 0 on both axes) - the direct
      ## cross-ancestry analog of Cross-Tissue Validation's blood-vs-synovium
      ## log2FC scatter. Untestable genes (no East Asian instrument) are
      ## excluded rather than plotted at a misleading position.
      output[[paste0(sex_label, "_concordance_plot")]] <- renderPlot({
        df <- res()$df
        d <- df[df$testable_EAS & !is.na(df$OR_stahl) & !is.na(df$OR_bbj) & df$OR_stahl > 0 & df$OR_bbj > 0, , drop = FALSE]
        req(nrow(d) > 0)
        d$Status <- factor(ifelse(d$biomarker, "Validated biomarker", "Not validated"), levels = c("Validated biomarker", "Not validated"))
        d$log_or_stahl <- log2(d$OR_stahl); d$log_or_bbj <- log2(d$OR_bbj)
        ggplot(d, aes(x = log_or_stahl, y = log_or_bbj, color = Status)) +
          geom_hline(yintercept = 0, color = ARTHOMIX_COLORS$axis, linewidth = 0.3) +
          geom_vline(xintercept = 0, color = ARTHOMIX_COLORS$axis, linewidth = 0.3) +
          geom_point(aes(size = Status), alpha = 0.9) +
          ggrepel::geom_text_repel(data = d[d$biomarker, , drop = FALSE], aes(label = gene), size = 3.4,
                                    color = ARTHOMIX_COLORS$ink, fontface = "bold", show.legend = FALSE, max.overlaps = 30) +
          scale_color_manual(values = c(`Validated biomarker` = ARTHOMIX_STATUS$good, `Not validated` = ARTHOMIX_COLORS$ink_muted)) +
          scale_size_manual(values = c(`Validated biomarker` = 3.6, `Not validated` = 2.2), guide = "none") +
          labs(title = "Ancestry concordance", subtitle = "European (Stahl) vs. East Asian (BioBank Japan) MR estimate",
               x = "log2 OR - European replication", y = "log2 OR - East Asian transfer", color = NULL) +
          theme_arthomix(base_size = 12)
      }, alt = sprintf("Scatter plot comparing each testable %s gene's European replication odds ratio against its East Asian transfer odds ratio, with validated cross-ancestry biomarkers highlighted and labelled.", sex_label))

      ## Cross-ancestry biomarker identification: genes testable in East
      ## Asian, ranked by transfer significance - a lollipop, not a bar
      ## chart, with a live threshold line at the current p_eas cutoff (this
      ## line moves as the slider moves) and only validated genes labelled
      ## with their East Asian OR.
      ## Height grows with the number of testable genes (one row each) so
      ## labels never overlap regardless of panel size - matched to
      ## plotOutput(height = "auto") in the UI, which reads this pixel value.
      rank_plot_height <- function() {
        df <- res()$df
        n <- sum(df$testable_EAS & !is.na(df$p_bbj))
        max(320, 60 + 20 * n)
      }
      output[[paste0(sex_label, "_rank_plot")]] <- renderPlot({
        df <- res()$df
        d <- df[df$testable_EAS & !is.na(df$p_bbj), , drop = FALSE]
        req(nrow(d) > 0)
        d$Status <- factor(ifelse(d$biomarker, "Validated biomarker", "Not validated"), levels = c("Validated biomarker", "Not validated"))
        d$neglogp <- -log10(pmax(d$p_bbj, 1e-300))
        d$gene <- factor(d$gene, levels = d$gene[order(d$neglogp)])
        ## Extreme hits (e.g. HLA associations at p < 1e-50) otherwise
        ## squash every other gene into an unreadable cluster near zero on a
        ## linear scale - capped for DISPLAY only; ranking, classification
        ## and the labelled OR value all still use the real p-value.
        CAP <- 20
        d$neglogp_plot <- pmin(d$neglogp, CAP)
        subtitle <- if (any(d$neglogp > CAP)) sprintf("East Asian transfer significance (capped at -log10 p = %d)", CAP) else "East Asian transfer significance"
        cutoff <- -log10(input$p_eas %||% 0.05)
        ggplot(d, aes(x = neglogp_plot, y = gene, color = Status)) +
          geom_vline(xintercept = cutoff, linetype = "dashed", color = ARTHOMIX_STATUS$good, linewidth = 0.4) +
          geom_segment(aes(x = 0, xend = neglogp_plot, yend = gene), linewidth = 0.7) +
          geom_point(aes(size = Status)) +
          geom_text(data = d[d$biomarker, , drop = FALSE], aes(x = neglogp_plot, label = sprintf("OR=%.2f", OR_bbj)),
                    hjust = -0.15, size = 3.2, color = ARTHOMIX_COLORS$ink, show.legend = FALSE) +
          scale_color_manual(values = c(`Validated biomarker` = ARTHOMIX_STATUS$good, `Not validated` = ARTHOMIX_COLORS$ink_muted)) +
          scale_size_manual(values = c(`Validated biomarker` = 3.4, `Not validated` = 2.4), guide = "none") +
          scale_x_continuous(expand = expansion(mult = c(0.02, 0.16))) +
          labs(title = "Cross-ancestry biomarkers", subtitle = subtitle,
               x = "-log10(p), BioBank Japan", y = NULL, color = NULL) +
          theme_arthomix(base_size = 12) + theme(legend.position = "bottom")
      }, height = rank_plot_height, alt = sprintf("Ranked dot plot of each testable %s gene's East Asian transfer significance, with validated cross-ancestry biomarkers highlighted and labelled against the current p-value cutoff.", sex_label))

      output[[paste0(sex_label, "_table")]] <- DT::renderDataTable({
        df <- res()$df
        d2 <- data.frame(
          gene = df$gene, cross_ancestry_biomarker = ifelse(df$biomarker, "✓", ""), ancestry_class = df$ancestry_class,
          OR_okada = round(df$OR_okada, 3), p_okada = signif(df$p_okada, 3),
          OR_stahl = round(df$OR_stahl, 3), p_stahl = signif(df$p_stahl, 3),
          OR_bbj = round(df$OR_bbj, 3), p_bbj = signif(df$p_bbj, 3),
          nSNP_bbj = df$nSNP_bbj, eafgap_bbj = round(df$eafgap_bbj, 3),
          stringsAsFactors = FALSE
        )
        d2 <- d2[order(-df$biomarker, df$p_bbj), ]
        DT::datatable(d2, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
      })
      output[[paste0(sex_label, "_download")]] <- downloadHandler(
        filename = function() sprintf("%s_crossancestry_recomputed.csv", sex_label),
        content = function(file) write.csv(res()$df, file, row.names = FALSE)
      )
    }

    register_sex_outputs("female", result_female)
    register_sex_outputs("male", result_male)

    ## A plain shinydashboard::box(collapsible = TRUE) never actually opens
    ## here: AdminLTE's box-widget click handler only binds to boxes present
    ## at initial page load, not ones inserted later via renderUI/uiOutput
    ## (this box's entire reason for existing) - so its "+" toggle looks like
    ## a working collapse but silently does nothing. A native <details>
    ## disclosure needs no JS at all and always works, dynamically inserted
    ## or not; box/box-header/box-title/box-body classes are reused purely
    ## for the visual match to every other box on this page, not their JS.
    output$references_box_ui <- renderUI({
      tags$details(
        class = "box box-primary",
        tags$summary(class = "box-header", style = "cursor: pointer;", tags$h3(class = "box-title", "References")),
        div(class = "box-body",
        tags$ul(
          class = "dge-ref-list",
          tags$li(strong("Discovery / replication / transfer GWAS: "), "Okada Y, et al. (2014). Genetics of rheumatoid arthritis contributes to biology and drug discovery. ", tags$em("Nature"), ", 506(7488), 376-381; Stahl EA, et al. (2010). ", tags$em("Nat Genet"), ", 42(6), 508-514; Ishigaki K, et al. (2020). ", tags$em("Nat Genet"), ", 52(7), 669-679."),
          tags$li(strong("Instrument construction & strength: "), "Võsa U, et al. (2021). ", tags$em("Nat Genet"), ", 53(9), 1300-1310; Burgess S, Thompson SG (2011). ", tags$em("Int J Epidemiol"), ", 40(3), 755-764."),
          tags$li(strong("MR estimators: "), "Burgess S, Butterworth A, Thompson SG (2013). ", tags$em("Genet Epidemiol"), ", 37(7), 658-665; Bowden J, Davey Smith G, Burgess S (2015). ", tags$em("Int J Epidemiol"), ", 44(2), 512-525; Bowden J, et al. (2016). ", tags$em("Genet Epidemiol"), ", 40(4), 304-314."),
          tags$li(strong("MR-Base / TwoSampleMR / OpenGWAS: "), "Hemani G, et al. (2018). ", tags$em("eLife"), ", 7, e34408; Elsworth B, et al. (2020). ", tags$em("bioRxiv"), "."),
          tags$li(strong("Cross-ancestry portability caveat: "), "Martin AR, et al. (2019). Clinical use of current polygenic risk scores may exacerbate health disparities. ", tags$em("Nat Genet"), ", 51(4), 584-591."),
          tags$li(strong("Multiple-testing correction (gene set entry, Section 2.6): "), "Benjamini Y, Hochberg Y (1995). ", tags$em("J R Stat Soc B"), ", 57(1), 289-300.")
        ),
        p(class = "submodule-desc", strong("Ask ArthOChat"), " for a plain-language walkthrough or a live citation.")
        )
      )
    })
  })
}
