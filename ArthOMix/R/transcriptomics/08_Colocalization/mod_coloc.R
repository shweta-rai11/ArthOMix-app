## R/transcriptomics/08_Colocalization/mod_coloc.R
## Submodule: Colocalization (Section 2.7)
## Runs coloc.abf between the bundled eQTL cis-window instrument (33 genes) and

mod_coloc_config <- list(
  id = "coloc", group = "Genetics",
  title = "Colocalization",
  description = "Bayesian colocalisation testing whether the eQTL and a GWAS's disease-risk association at each candidate locus share a single causal variant.",
  icon = "map-location-dot"
)

mod_coloc_ui <- function(id) {
  ns <- NS(id)
  tagList(
      fluidRow(
        column(
          4,
          box(
            width = NULL, title = "Candidate region & GWAS", status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "Overlapping eQTL and GWAS peaks can mean two distinct causal variants in LD, not one variant driving both. coloc.abf tests that directly: for the chosen gene's cis-window it compares eQTL vs GWAS association patterns and returns posterior probabilities for five hypotheses, plus per-SNP evidence and regional/posterior-probability plots. The eQTL side is always this project's bundled cis-window instrument (33 genes); only the GWAS side can be swapped between the bundled RA GWAS and an uploaded GWAS for another trait."),
            radioButtons(
              ns("data_source"), NULL,
              choiceNames = list(
                tagList(icon("database"), " Bundled RA GWAS (default)"),
                tagList(icon("upload"), " Upload your own GWAS summary statistics")
              ),
              choiceValues = list("project", "upload"), selected = "project"
            ),
            uiOutput(ns("controls")),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'upload'", ns("data_source")),
              p(class = "submodule-desc",
                "A delimited file (CSV/TSV), one row per SNP: any trait's GWAS summary statistics, tested against the bundled eQTL region for the gene picked above. SNP IDs must be rsIDs (dbSNP rs#) to match the bundled eQTL instrument, and a sample-size column is required (coloc.abf needs N for both sides)."),
              textInput(ns("gwas_label"), "GWAS trait name (for labelling only)", value = "Uploaded GWAS", width = "100%"),
              radioButtons(
                ns("gwas_type"), "GWAS trait type", inline = TRUE,
                choices = c("Binary (case/control)" = "cc", "Quantitative" = "quant"), selected = "cc"
              ),
              fileInput(ns("gwas_file"), "GWAS file", accept = c(".csv", ".tsv", ".txt")),
              uiOutput(ns("gwas_map_ui"))
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'project' || input['%s'] == 'cc'", ns("data_source"), ns("gwas_type")),
              sliderInput(ns("case_frac"), "Assumed case fraction in the GWAS", value = 0.33, min = 0.05, max = 0.5, step = 0.01)
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'upload' && input['%s'] == 'quant'", ns("data_source"), ns("gwas_type")),
              p(class = "submodule-desc",
                "Quantitative trait selected: coloc.abf needs an effect-allele-frequency column for this GWAS to analyse it without a directly-known phenotype SD - map one below (the optional \"Effect allele frequency\" field).")
            ),
            tags$details(
              class = "box box-primary", style = "margin-top: 10px;",
              tags$summary(class = "box-header", style = "cursor: pointer;", tags$h3(class = "box-title", "Advanced: colocalisation priors")),
              div(class = "box-body",
                p(class = "submodule-desc",
                  "coloc.abf's Bayesian priors: the probability a given SNP is associated with the eQTL only (p1), the GWAS trait only (p2), or both (p12). These drive the PP.H4 (shared causal variant) conclusion, so it's worth sensitivity-testing them rather than trusting the package defaults blindly. Defaults below are coloc's own documented conventions and leave existing results unchanged unless you edit them."),
                fluidRow(
                  column(4, numericInput(ns("p1"), "p1 (eQTL only)", value = 1e-4, min = 1e-8, max = 1e-2, step = 1e-6, width = "100%")),
                  column(4, numericInput(ns("p2"), "p2 (GWAS only)", value = 1e-4, min = 1e-8, max = 1e-2, step = 1e-6, width = "100%")),
                  column(4, numericInput(ns("p12"), "p12 (both)", value = 1e-5, min = 1e-8, max = 1e-2, step = 1e-7, width = "100%"))
                ),
                helpText("coloc's conventional defaults: p1 = 1e-4, p2 = 1e-4, p12 = 1e-5. p12 should not exceed either p1 or p2 (a SNP can't be more likely to affect both traits than it is to affect either one alone).")
              )
            ),
            actionButton(ns("run_btn"), "Run colocalisation", icon = icon("play"), class = "btn-primary btn-sm")
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] > 0", ns("run_btn")),
            box(
              width = NULL, title = "Result", status = "primary", solidHeader = FALSE,
              withSpinner(uiOutput(ns("summary_ui")), color = "#2c6fbb", type = 6),
              withSpinner(plotOutput(ns("pp_plot"), height = 260), color = "#2c6fbb", type = 6)
            )
          )
        ),
        column(
          8,
          conditionalPanel(
            condition = sprintf("input['%s'] > 0", ns("run_btn")),
            box(
              width = NULL, title = "Regional association", status = "primary", solidHeader = FALSE,
              withSpinner(plotOutput(ns("region_plot"), height = 460), color = "#2c6fbb", type = 6)
            )
          )
        )
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] > 0", ns("run_btn")),
        box(
          width = 12, title = "SNP-level results", status = "primary", solidHeader = FALSE,
          div(class = "table-toolbar", downloadButton(ns("download_coloc"), "Download CSV", class = "btn-sm")),
          DT::dataTableOutput(ns("snp_table"))
        )
      )
  )
}

mod_coloc_server <- function(id, dataset, results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    coloc_regions <- readRDS(COLOC_REGIONS_RDS)
    available_genes <- sort(names(coloc_regions))

    output$controls <- renderUI({
      tagList(
        selectInput(ns("gene"), "Candidate gene", choices = available_genes, selectize = FALSE),
        p(class = "submodule-desc", style = "margin-bottom: 0; font-size: 12px;",
          sprintf("%d genes have a bundled cis-window eQTL instrument available - this is always the eQTL side, in both modes above.", length(available_genes)))
      )
    })

    gwas_df_r <- reactive({ req(input$gwas_file); read_uploaded_table(input$gwas_file$datapath) })
    output$gwas_map_ui <- gwas_col_map_ui(ns, reactive(input$gwas_file), gwas_df_r, "gwas", "GWAS file", extra_fields = "n")

    validated_priors <- function() {
      p1 <- input$p1 %||% 1e-4
      p2 <- input$p2 %||% 1e-4
      p12 <- input$p12 %||% 1e-5
      validate(need(is.numeric(p1) && length(p1) == 1 && !is.na(p1) && p1 > 0 && p1 < 1,
                    "p1 must be a probability strictly between 0 and 1."))
      validate(need(is.numeric(p2) && length(p2) == 1 && !is.na(p2) && p2 > 0 && p2 < 1,
                    "p2 must be a probability strictly between 0 and 1."))
      validate(need(is.numeric(p12) && length(p12) == 1 && !is.na(p12) && p12 > 0 && p12 < 1,
                    "p12 must be a probability strictly between 0 and 1."))
      validate(need(p12 <= min(p1, p2),
                    "p12 (probability a SNP affects both the eQTL and the GWAS trait) cannot exceed p1 or p2 (probability it affects only one) - this is a coloc sanity convention. Lower p12, or raise p1/p2."))
      list(p1 = p1, p2 = p2, p12 = p12)
    }

    coloc_result_project <- function() {
      req(input$gene)
      r <- coloc_regions[[input$gene]]
      eqtl <- as.data.frame(r$eqtl)
      gwas <- as.data.frame(r$gwas)

      common <- intersect(eqtl$rsid, gwas$rsid)
      validate(need(length(common) >= 10, "Fewer than 10 SNPs are shared between the eQTL and GWAS summary statistics for this region."))
      e <- eqtl[match(common, eqtl$rsid), ]
      g <- gwas[match(common, gwas$rsid), ]

      ok <- complete.cases(e[, c("beta", "se", "eaf", "n")]) & complete.cases(g[, c("beta", "se", "n")]) &
        e$se > 0 & g$se > 0 & e$eaf > 0 & e$eaf < 1
      validate(need(sum(ok) >= 10, "Fewer than 10 SNPs have complete summary statistics for this region."))
      common <- common[ok]; e <- e[ok, ]; g <- g[ok, ]

      priors <- validated_priors()
      d1 <- list(beta = e$beta, varbeta = e$se^2, N = round(median(e$n)), MAF = pmin(e$eaf, 1 - e$eaf), type = "quant", snp = common)
      d2 <- list(beta = g$beta, varbeta = g$se^2, N = round(median(g$n)), type = "cc", s = input$case_frac, snp = common)
      res <- tryCatch(suppressWarnings(coloc::coloc.abf(dataset1 = d1, dataset2 = d2, p1 = priors$p1, p2 = priors$p2, p12 = priors$p12)), error = function(e) e)
      validate(need(!inherits(res, "error"), sprintf("coloc.abf() failed: %s", if (inherits(res, "error")) conditionMessage(res) else "unknown error")))

      snp_df <- data.frame(
        snp = common, eqtl_beta = e$beta, eqtl_p = e$p, gwas_beta = g$beta, gwas_p = g$p, pos = e$position,
        snp_pp_h4 = res$results$SNP.PP.H4[match(common, res$results$snp)]
      )

      list(gene = input$gene, gwas_label = "Bundled RA GWAS (Okada 2014)", gwas_type = "cc", summary = res$summary,
           snp_df = snp_df, n_snp = length(common), uploaded = FALSE, priors = priors)
    }

    coloc_result_uploaded <- function() {
      req(input$gene, input$gwas_file)
      req(input$gwas_snp, input$gwas_beta, input$gwas_se, input$gwas_pval, input$gwas_ea, input$gwas_oa, input$gwas_n)

      gwas_type <- input$gwas_type %||% "cc"
      if (identical(gwas_type, "quant")) {
        validate(need(nzchar(input$gwas_eaf %||% ""),
          "Quantitative-trait colocalisation needs an effect-allele-frequency column for the GWAS (coloc.abf needs it to analyse a quantitative trait without a directly-known phenotype SD) - map one in \"Effect allele frequency\" below, or switch GWAS trait type to Binary if this GWAS has no allele-frequency column."))
      }

      r <- coloc_regions[[input$gene]]
      eqtl <- as.data.frame(r$eqtl)
      gwas_raw <- gwas_df_r()
      validate(need(!is.null(gwas_raw), "Could not read the uploaded GWAS file."))

      label <- if (nzchar(trimws(input$gwas_label %||% ""))) trimws(input$gwas_label) else "Uploaded GWAS"

      exp_fmt <- TwoSampleMR::format_data(
        eqtl, type = "exposure", snp_col = "rsid", beta_col = "beta", se_col = "se",
        pval_col = "p", effect_allele_col = "ea", other_allele_col = "nea",
        eaf_col = "eaf", samplesize_col = "n", chr_col = "chr", pos_col = "position"
      )
      exp_fmt$exposure <- sprintf("%s eQTL (blood, bundled)", input$gene)
      out_fmt <- TwoSampleMR::format_data(
        gwas_raw, type = "outcome", snp_col = input$gwas_snp, beta_col = input$gwas_beta,
        se_col = input$gwas_se, pval_col = input$gwas_pval, effect_allele_col = input$gwas_ea,
        other_allele_col = input$gwas_oa, eaf_col = if (nzchar(input$gwas_eaf %||% "")) input$gwas_eaf else "eaf",
        samplesize_col = input$gwas_n
      )
      out_fmt$outcome <- label

      dat_up <- tryCatch(TwoSampleMR::harmonise_data(exp_fmt, out_fmt, action = 2), error = function(e) NULL)
      validate(need(!is.null(dat_up) && nrow(dat_up) > 0,
        "Harmonisation found no overlapping SNPs between the bundled eQTL region and your uploaded GWAS - check that both use rsIDs and that allele columns are mapped correctly."))
      dat_up <- dat_up[dat_up$mr_keep, , drop = FALSE]
      dat_up <- dat_up[!duplicated(dat_up$SNP), , drop = FALSE]

      out_cols <- if (identical(gwas_type, "quant")) c("beta.outcome", "se.outcome", "samplesize.outcome", "eaf.outcome")
                  else c("beta.outcome", "se.outcome", "samplesize.outcome")
      ok <- complete.cases(dat_up[, c("beta.exposure", "se.exposure", "eaf.exposure", "samplesize.exposure")]) &
        complete.cases(dat_up[, out_cols]) &
        dat_up$se.exposure > 0 & dat_up$se.outcome > 0 & dat_up$eaf.exposure > 0 & dat_up$eaf.exposure < 1
      if (identical(gwas_type, "quant")) ok <- ok & dat_up$eaf.outcome > 0 & dat_up$eaf.outcome < 1
      validate(need(sum(ok) >= 10,
        "Fewer than 10 SNPs have complete, harmonised summary statistics (with a valid sample size on both sides) for this region."))
      dat_up <- dat_up[ok, , drop = FALSE]

      priors <- validated_priors()
      d1 <- list(beta = dat_up$beta.exposure, varbeta = dat_up$se.exposure^2,
                 N = round(median(dat_up$samplesize.exposure)), MAF = pmin(dat_up$eaf.exposure, 1 - dat_up$eaf.exposure),
                 type = "quant", snp = dat_up$SNP)
      d2 <- if (identical(gwas_type, "quant")) {
        list(beta = dat_up$beta.outcome, varbeta = dat_up$se.outcome^2, N = round(median(dat_up$samplesize.outcome)),
             MAF = pmin(dat_up$eaf.outcome, 1 - dat_up$eaf.outcome), type = "quant", snp = dat_up$SNP)
      } else {
        list(beta = dat_up$beta.outcome, varbeta = dat_up$se.outcome^2,
             N = round(median(dat_up$samplesize.outcome)), type = "cc", s = input$case_frac, snp = dat_up$SNP)
      }
      res <- tryCatch(suppressWarnings(coloc::coloc.abf(dataset1 = d1, dataset2 = d2, p1 = priors$p1, p2 = priors$p2, p12 = priors$p12)), error = function(e) e)
      validate(need(!inherits(res, "error"), sprintf("coloc.abf() failed: %s", if (inherits(res, "error")) conditionMessage(res) else "unknown error")))

      snp_df <- data.frame(
        snp = dat_up$SNP, eqtl_beta = dat_up$beta.exposure, eqtl_p = dat_up$pval.exposure,
        gwas_beta = dat_up$beta.outcome, gwas_p = dat_up$pval.outcome, pos = dat_up$pos.exposure,
        snp_pp_h4 = res$results$SNP.PP.H4[match(dat_up$SNP, res$results$snp)]
      )

      list(gene = input$gene, gwas_label = label, gwas_type = gwas_type, summary = res$summary,
           snp_df = snp_df, n_snp = nrow(dat_up), uploaded = TRUE, priors = priors)
    }

    coloc_result <- eventReactive(input$run_btn, {
      if (identical(input$data_source, "upload")) coloc_result_uploaded() else coloc_result_project()
    }, ignoreInit = TRUE)

    coloc_has_run <- reactiveVal(FALSE)
    observeEvent(input$run_btn, coloc_has_run(TRUE), ignoreInit = TRUE)

    observeEvent(coloc_result(), {
      res <- coloc_result()
      entry <- list(n_snp = res$n_snp, pp_h4 = round(unname(res$summary["PP.H4.abf"]), 3),
                    gwas = res$gwas_label, uploaded = res$uploaded)
      existing <- results$coloc %||% list(genes_tested = list())
      existing$genes_tested[[res$gene]] <- entry
      results$coloc <- existing
    })

    output$summary_ui <- renderUI({
      if (!coloc_has_run()) {
        return(div(class = "empty-note", icon("circle-info"),
                    "Not run yet. Pick a gene on the left, then click \"Run colocalisation\"."))
      }
      res <- coloc_result()
      pp4 <- res$summary["PP.H4.abf"]
      tagList(
        if (isTRUE(res$uploaded)) p(class = "empty-note", icon("upload"),
          sprintf("Tested against your uploaded %s GWAS (%s), harmonised against the bundled eQTL region - no bundled reference estimate applies to custom data.",
                  if (identical(res$gwas_type, "quant")) "quantitative-trait" else "case/control", res$gwas_label)),
        p(strong(res$gene), ": ", strong(res$n_snp), " shared SNPs tested against ", res$gwas_label, "."),
        p(class = "submodule-desc", style = "font-size:12px;",
          sprintf("Priors used: p1 = %.3g, p2 = %.3g, p12 = %.3g.", res$priors$p1, res$priors$p2, res$priors$p12)),
        p("Posterior probability of a shared causal variant (PP.H4): ", strong(sprintf("%.1f%%", 100 * pp4))),
        if (pp4 > 0.8) p(class = "empty-note", icon("circle-check"), "Strong support for colocalisation at this locus.")
        else if (pp4 > 0.5) p(class = "empty-note", icon("circle-info"), "Moderate support for colocalisation at this locus.")
        else p(class = "empty-note", icon("circle-exclamation"), "Little support for a single shared causal variant at this locus.")
      )
    })

    output$pp_plot <- renderPlot({
      req(coloc_has_run())
      res <- coloc_result()
      df <- data.frame(
        hypothesis = factor(c("H0", "H1", "H2", "H3", "H4"), levels = c("H0", "H1", "H2", "H3", "H4")),
        pp = as.numeric(res$summary[c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")])
      )
      ggplot(df, aes(x = hypothesis, y = pp, fill = hypothesis == "H4")) +
        geom_col() +
        scale_fill_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "#8a929c"), guide = "none") +
        labs(x = NULL, y = "Posterior probability") +
        theme_minimal(base_size = 12)
    })

    output$region_plot <- renderPlot({
      req(coloc_has_run())
      res <- coloc_result()
      df <- res$snp_df
      eqtl_df <- data.frame(pos = df$pos, nlp = -log10(df$eqtl_p), track = "eQTL (blood)")
      gwas_df <- data.frame(pos = df$pos, nlp = -log10(df$gwas_p), track = res$gwas_label)
      plot_df <- rbind(eqtl_df, gwas_df)
      ggplot(plot_df, aes(x = pos, y = nlp, color = track)) +
        geom_point(alpha = 0.7, size = 1.6) +
        facet_wrap(~track, ncol = 1, scales = "free_y") +
        scale_color_manual(values = stats::setNames(c("#2c6fbb", "#c0392b"), c("eQTL (blood)", res$gwas_label)), guide = "none") +
        labs(x = "Position", y = "-log10 p-value") +
        theme_minimal(base_size = 12)
    })

    output$snp_table <- DT::renderDataTable({
      req(coloc_has_run())
      DT::datatable(coloc_result()$snp_df, rownames = FALSE, filter = "top",
                     options = list(pageLength = 12, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_coloc <- downloadHandler(
      filename = function() paste0("coloc_", coloc_result()$gene, ".csv"),
      content = function(file) {
        res <- coloc_result()
        df <- res$snp_df
        df$prior_p1 <- res$priors$p1
        df$prior_p2 <- res$priors$p2
        df$prior_p12 <- res$priors$p12
        write.csv(df, file, row.names = FALSE)
      }
    )
  })
}
