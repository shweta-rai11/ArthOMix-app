## R/mod_coloc.R
## Submodule: Colocalization (Section 2.7)
## "Your analysis" runs a live coloc.abf test for a chosen candidate gene's
## region, using the same eQTL and RA GWAS summary statistics
## (coloc_regions.rds) as the thesis pipeline - or, in "Upload your own GWAS
## summary statistics" mode, the SAME bundled eQTL region tested against any
## other trait's GWAS a user brings.
##
## The eQTL side is ALWAYS the bundled cis-window instrument for the chosen
## gene, in both modes: coloc_regions.rds only has a prepared eQTL region for
## 33 genes, and building an equivalent for an arbitrary gene needs the full
## eQTLGen dataset this app doesn't otherwise load. Only the GWAS side is
## swappable.
##
## Unlike the bundled path (both sides pre-aligned by the pipeline, so a
## plain rsid intersection is safe), an uploaded GWAS has no guaranteed
## allele coding relative to the bundled eQTL - the upload branch below runs
## both through TwoSampleMR::format_data()/harmonise_data() (the SAME
## functions mod_mr.R's own upload mode uses) to align alleles and flip beta
## sign consistently before handing off to coloc.abf(), rather than assuming
## an arbitrary upload's beta sign already matches.

mod_coloc_config <- list(
  id = "coloc", group = "Genetics",
  title = "Colocalization",
  description = "Bayesian colocalisation (coloc.abf) testing whether the eQTL and a GWAS's disease-risk association at each candidate locus share a single causal variant - the bundled RA GWAS by default, or upload your own GWAS summary statistics for any trait.",
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
            p(class = "submodule-desc", "Runs coloc.abf testing whether the eQTL association and a GWAS's disease-risk association at one gene's region share a single causal variant. The eQTL side always uses this project's own bundled cis-window instrument for the gene picked below - only the GWAS side is swappable."),
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
              fileInput(ns("gwas_file"), "GWAS file", accept = c(".csv", ".tsv", ".txt")),
              uiOutput(ns("gwas_map_ui"))
            ),
            sliderInput(ns("case_frac"), "Assumed case fraction in the GWAS", value = 0.33, min = 0.05, max = 0.5, step = 0.01),
            actionButton(ns("run_btn"), "Run colocalisation", icon = icon("play"), class = "btn-primary btn-sm")
          ),
          box(
            width = NULL, title = "Result", status = "primary", solidHeader = FALSE,
            withSpinner(uiOutput(ns("summary_ui")), color = "#2c6fbb", type = 6),
            withSpinner(plotOutput(ns("pp_plot"), height = 260), color = "#2c6fbb", type = 6)
          )
        ),
        column(
          8,
          box(
            width = NULL, title = "Regional association", status = "primary", solidHeader = FALSE,
            withSpinner(plotOutput(ns("region_plot"), height = 460), color = "#2c6fbb", type = 6)
          )
        )
      ),
      box(
        width = 12, title = "SNP-level results", status = "primary", solidHeader = FALSE,
        div(class = "table-toolbar", downloadButton(ns("download_coloc"), "Download CSV", class = "btn-sm")),
        DT::dataTableOutput(ns("snp_table"))
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

    ## Upload-your-own-GWAS mode: parse + map the uploaded file with the
    ## SAME shared helpers mod_mr.R's own upload mode uses (global.R), plus
    ## a required "n" (sample size) column mod_mr.R's estimator never needed
    ## but coloc.abf does.
    gwas_df_r <- reactive({ req(input$gwas_file); read_uploaded_table(input$gwas_file$datapath) })
    output$gwas_map_ui <- gwas_col_map_ui(ns, reactive(input$gwas_file), gwas_df_r, "gwas", "GWAS file", extra_fields = "n")

    ## Bundled path - unchanged from before upload mode existed. Both sides
    ## were pre-aligned by the pipeline, so a plain rsid intersection (no
    ## allele harmonisation) is safe here.
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

      d1 <- list(beta = e$beta, varbeta = e$se^2, N = round(median(e$n)), MAF = pmin(e$eaf, 1 - e$eaf), type = "quant", snp = common)
      d2 <- list(beta = g$beta, varbeta = g$se^2, N = round(median(g$n)), type = "cc", s = input$case_frac, snp = common)
      res <- suppressWarnings(coloc::coloc.abf(dataset1 = d1, dataset2 = d2))

      snp_df <- data.frame(
        snp = common, eqtl_beta = e$beta, eqtl_p = e$p, gwas_beta = g$beta, gwas_p = g$p, pos = e$position,
        snp_pp_h4 = res$results$SNP.PP.H4[match(common, res$results$snp)]
      )

      list(gene = input$gene, gwas_label = "Bundled RA GWAS (Okada 2014)", summary = res$summary,
           snp_df = snp_df, n_snp = length(common), uploaded = FALSE)
    }

    ## Upload path: format the bundled eQTL region and the uploaded GWAS as
    ## TwoSampleMR exposure/outcome objects and harmonise them (allele
    ## alignment + beta-sign flipping) before building coloc.abf's d1/d2 -
    ## the harmonisation step the bundled path gets for free but an
    ## arbitrary upload needs done live. `snp_df`'s shape matches the
    ## bundled path exactly, so pp_plot/region_plot/snp_table/download below
    ## need no branching of their own.
    coloc_result_uploaded <- function() {
      req(input$gene, input$gwas_file)
      req(input$gwas_snp, input$gwas_beta, input$gwas_se, input$gwas_pval, input$gwas_ea, input$gwas_oa, input$gwas_n)

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

      ok <- complete.cases(dat_up[, c("beta.exposure", "se.exposure", "eaf.exposure", "samplesize.exposure")]) &
        complete.cases(dat_up[, c("beta.outcome", "se.outcome", "samplesize.outcome")]) &
        dat_up$se.exposure > 0 & dat_up$se.outcome > 0 & dat_up$eaf.exposure > 0 & dat_up$eaf.exposure < 1
      validate(need(sum(ok) >= 10,
        "Fewer than 10 SNPs have complete, harmonised summary statistics (with a valid sample size on both sides) for this region."))
      dat_up <- dat_up[ok, , drop = FALSE]

      d1 <- list(beta = dat_up$beta.exposure, varbeta = dat_up$se.exposure^2,
                 N = round(median(dat_up$samplesize.exposure)), MAF = pmin(dat_up$eaf.exposure, 1 - dat_up$eaf.exposure),
                 type = "quant", snp = dat_up$SNP)
      d2 <- list(beta = dat_up$beta.outcome, varbeta = dat_up$se.outcome^2,
                 N = round(median(dat_up$samplesize.outcome)), type = "cc", s = input$case_frac, snp = dat_up$SNP)
      res <- suppressWarnings(coloc::coloc.abf(dataset1 = d1, dataset2 = d2))

      snp_df <- data.frame(
        snp = dat_up$SNP, eqtl_beta = dat_up$beta.exposure, eqtl_p = dat_up$pval.exposure,
        gwas_beta = dat_up$beta.outcome, gwas_p = dat_up$pval.outcome, pos = dat_up$pos.exposure,
        snp_pp_h4 = res$results$SNP.PP.H4[match(dat_up$SNP, res$results$snp)]
      )

      list(gene = input$gene, gwas_label = label, summary = res$summary,
           snp_df = snp_df, n_snp = nrow(dat_up), uploaded = TRUE)
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
          sprintf("Tested against your uploaded GWAS (%s), harmonised against the bundled eQTL region - no bundled reference estimate applies to custom data.", res$gwas_label)),
        p(strong(res$gene), ": ", strong(res$n_snp), " shared SNPs tested against ", res$gwas_label, "."),
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
      content = function(file) write.csv(coloc_result()$snp_df, file, row.names = FALSE)
    )
  })
}
