## R/mod_nomogram.R
## Submodule: Clinical Utility Nomogram (Section 2.15)
## "Your analysis" fits a per-sex rms::lrm logistic model on a chosen gene
## panel (unstandardised log2 expression, no covariate standardisation - as
## unstandardised log2 expression, matching METHODS_2.15.3) and evaluates it
## the same FOUR ways this project's own script
## (scripts/goal2_sex_stratified/19_testing_blood_clinical_utility.R) does:
##   (A) NOMOGRAM       - rms::nomogram, points per gene -> total -> risk
##   (B) CALIBRATION    - apparent + bootstrap bias-corrected (rms::calibrate)
##   (C) DECISION CURVE - net benefit of the panel vs treat-all vs treat-none
##       (Vickers & Elkin 2006), with the secondary cost:benefit ratio axis
##   (D) CLINICAL IMPACT - per-N patients flagged high-risk, and how many of
##       those are genuine cases, with bootstrap 95% bands
## exactly as METHODS_2.15_clinical_utility_nomogram.md documents: no splines,
## no interactions, no standardisation, sex never a predictor (it is the
## partition), and a project-default ridge penalty of 5 for the male fit
## (near-complete separation at this cohort's male sample size) vs 0 for
## female - overridable below as an advanced filter.
##
## FOUR gene-panel sources, one radio choice (mirrors mod_enrichment.R's own
## "Your own / Same-tissue / Cross-tissue / Cross-ancestry" picker), each
## coupled to whichever expression cohort that source actually means:
##   "Your own"       - pasted gene list; user also picks which cohort to
##                       fit on (blood or synovium).
##   "Same-tissue"     - this session's live Feature Selection consensus
##                       panel (or this project's bundled one), fit on the
##                       currently loaded blood dataset - the project's
##                       primary instrument (Section 2.15).
##   "Cross-tissue"    - the SAME blood-derived panel, refit from scratch in
##       the independent synovial cohort (GSE89408 / val_synovium.rds) -
##       never a transfer of the blood model's coefficients, same principle
##       mod_crosstissue.R documents at length in its own header.
##   "Cross-ancestry"  - genes whose Mendelian-randomisation causal signal
##       replicates in European AND transfers to East-Asian ancestry GWAS
##       (results$crossancestry / bundled MR35_crossancestry_*.csv), fit on
##       the blood dataset - this project has no ancestry-stratified
##       expression cohort, so (exactly as mod_crossancestry.R's own header
##       explains for its MR panel) the ancestry evidence is genetic, not a
##       re-fit expression model in another ancestry.
##
## Fitting is per-sex (Female/Male radio, always shown - Section 2.15.2: sex
## is the partition, never a predictor), triggered by one "Build nomogram"
## button - adjusting a filter afterwards does nothing until clicked again,
## the same click-to-lock-in-results convention as every other Validation
## submodule in this app.

mod_nomogram_config <- list(
  id = "nomogram", group = "Interpretation",
  title = "Clinical Utility Nomogram",
  description = "A clinical nomogram built from the final gene panel and evaluated for net clinical benefit using decision curve analysis.",
  icon = "file-medical"
)

## ---------------------------------------------------------------------------
## Shared fitting helpers (pure functions - no input/results/session; every
## mod_*.R file is sourced into one shared environment, so these names only
## need to stay unique within this app, hence the "nom_" prefix - same
## convention as mod_crosstissue.R's "ct_" and mod_crossancestry.R's "ca_").
## ---------------------------------------------------------------------------

NOM_DEFAULT_PARAMS <- list(
  penalty_mode = "project", penalty_manual = 0,
  calibrate_B = 200, dca_step = 0.01,
  impact_N = 1000, impact_B = 500, seed = 1234
)

## Vickers & Elkin (2006) net benefit: NB = TP/n - FP/n * pt/(1-pt).
## `all` is the prevalence-determined "treat everyone" reference; "treat
## no one" is zero by definition and not computed. Matches this project's
## own dca() in 19_testing_blood_clinical_utility.R exactly.
nom_dca <- function(y, p, th) {
  n <- length(y); ev <- mean(y == 1)
  model <- vapply(th, function(pt) {
    pos <- p >= pt
    sum(pos & y == 1) / n - (sum(pos & y == 0) / n) * (pt / (1 - pt))
  }, numeric(1))
  all <- vapply(th, function(pt) ev - (1 - ev) * (pt / (1 - pt)), numeric(1))
  list(model = model, all = all)
}

## Clinical-impact curve (Kerr et al. 2016): per N patients, how many are
## flagged high-risk at each threshold, and how many of those are genuine
## cases - point estimate from the apparent fit, with 95% bands from B
## bootstrap resamples, each refit under the SAME penalty as the primary
## fit and evaluated on its own resampled observations (so the bands
## describe sampling variability, not optimism-corrected performance) -
## matches this project's own script section (D) exactly.
nom_clinical_impact <- function(df, form, penalty, p, y, th, N, B, seed) {
  nhigh  <- vapply(th, function(pt) mean(p >= pt) * N, numeric(1))
  nevent <- vapply(th, function(pt) mean(p >= pt & y == 1) * N, numeric(1))
  bh <- matrix(NA_real_, length(th), B)
  be <- matrix(NA_real_, length(th), B)
  for (b in seq_len(B)) {
    set.seed(seed + b)
    idx <- sample(nrow(df), replace = TRUE)
    dfb <- df[idx, , drop = FALSE]
    fb <- tryCatch(rms::lrm(form, data = dfb, penalty = penalty), error = function(e) NULL)
    if (is.null(fb)) next
    pb <- tryCatch(as.numeric(plogis(stats::predict(fb, newdata = dfb))), error = function(e) NULL)
    if (is.null(pb)) next
    yb <- dfb$y
    bh[, b] <- vapply(th, function(pt) mean(pb >= pt) * N, numeric(1))
    be[, b] <- vapply(th, function(pt) mean(pb >= pt & yb == 1) * N, numeric(1))
  }
  q <- function(M, pr) apply(M, 1, stats::quantile, probs = pr, na.rm = TRUE)
  list(threshold = th, nhigh = nhigh, nevent = nevent,
       nhigh_lo = q(bh, .025), nhigh_hi = q(bh, .975),
       nevent_lo = q(be, .025), nevent_hi = q(be, .975), N = N, B = B)
}

## One sex's full instrument: fit + nomogram + calibration + decision curve
## + clinical impact + coefficients, from a data.frame `df` already
## restricted to one sex, with unstandardised predictor columns plus a
## 0/1 `y` outcome column. `event_label` names the "1" outcome (e.g. "RA").
nom_fit_core <- function(df, predictors, penalty, event_label, params) {
  form <- as.formula(paste("y ~", paste(predictors, collapse = " + ")))

  ## datadist must be a GLOBAL option (rms's own requirement) - registered
  ## under a private name and torn down on exit, same pattern the original
  ## version of this module used, so a nomogram built here never leaks a
  ## stale datadist into another session's fit.
  assign(".arthomix_nomogram_dd", rms::datadist(df), envir = .GlobalEnv)
  old_opts <- options(datadist = ".arthomix_nomogram_dd")
  on.exit({
    options(old_opts)
    if (exists(".arthomix_nomogram_dd", envir = .GlobalEnv)) rm(".arthomix_nomogram_dd", envir = .GlobalEnv)
  }, add = TRUE)

  fit <- rms::lrm(form, data = df, x = TRUE, y = TRUE, penalty = penalty)
  nom <- rms::nomogram(fit, fun = plogis, funlabel = paste0("Risk of ", event_label),
                        fun.at = c(0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99), lp = FALSE)

  set.seed(params$seed)
  cal <- tryCatch(rms::calibrate(fit, B = params$calibrate_B), error = function(e) NULL)

  p <- as.numeric(plogis(stats::predict(fit)))
  y <- df$y
  th <- seq(0.01, 0.99, by = params$dca_step)
  dca_raw <- nom_dca(y, p, th)
  dca_df <- data.frame(threshold = th, NB_panel = dca_raw$model, NB_all = dca_raw$all, NB_none = 0)

  impact <- nom_clinical_impact(df, form, penalty, p, y, th, params$impact_N, params$impact_B, params$seed)

  coef_df <- data.frame(term = names(coef(fit)), coefficient = round(as.numeric(coef(fit)), 4), row.names = NULL)

  list(
    fit = fit, nom = nom, cal = cal, dca_df = dca_df, impact = impact, coef_df = coef_df,
    c_stat = unname(fit$stats["C"]), n_samples = nrow(df),
    n_pos = sum(y == 1), n_neg = sum(y == 0), penalty = penalty, event_label = event_label
  )
}

## Decision-curve plot with the secondary Cost:Benefit ratio axis (rmda
## convention, ticks at 1:100 .. 100:1) - matches this project's own script
## section (C) exactly, base graphics so it composes the same way the
## nomogram/calibration plots already do in this app.
nom_render_dca_plot <- function(dca_df, event_label) {
  th <- dca_df$threshold
  ymax <- suppressWarnings(max(c(dca_df$NB_panel, dca_df$NB_all), na.rm = TRUE))
  if (!is.finite(ymax) || ymax <= 0) ymax <- 0.1
  op <- par(mar = c(7.5, 4, 3, 1)); on.exit(par(op))
  plot(th, dca_df$NB_panel, type = "l", col = "#C0392B", lwd = 2.5, ylim = c(-0.05, ymax * 1.1),
       xlab = "", ylab = "Net benefit", main = "")
  mtext("Threshold probability", side = 1, line = 2.2)
  lines(th, dca_df$NB_all, col = "grey55", lwd = 1.5)
  abline(h = 0, col = "black", lwd = 1.2)
  axis(1, at = c(1 / 101, 0.2, 0.4, 0.6, 0.8, 100 / 101),
       labels = c("1:100", "1:4", "2:3", "3:2", "4:1", "100:1"), line = 3.2, cex.axis = 0.85)
  mtext("Cost:Benefit Ratio", side = 1, line = 5.4)
  legend("topright", c(paste("Panel genes -", event_label), "Treat all", "Treat none"),
         col = c("#C0392B", "grey55", "black"), lwd = c(2.5, 1.5, 1.2), bty = "n")
}

## Clinical-impact plot - matches this project's own script section (D).
nom_render_impact_plot <- function(impact, event_label) {
  th <- impact$threshold; N <- impact$N
  op <- par(mar = c(7.5, 4, 3, 1)); on.exit(par(op))
  plot(th, impact$nhigh, type = "l", col = "#C0392B", lwd = 2.5, ylim = c(0, N),
       xlab = "", ylab = sprintf("Number per %s", format(N, big.mark = ",")), main = "")
  mtext("High-risk threshold", side = 1, line = 2.2)
  lines(th, impact$nhigh_lo, col = "#C0392B", lwd = 0.8); lines(th, impact$nhigh_hi, col = "#C0392B", lwd = 0.8)
  lines(th, impact$nevent, col = "#2E86C1", lwd = 2.5, lty = 2)
  lines(th, impact$nevent_lo, col = "#2E86C1", lwd = 0.8, lty = 3); lines(th, impact$nevent_hi, col = "#2E86C1", lwd = 0.8, lty = 3)
  axis(1, at = c(1 / 101, 0.2, 0.4, 0.6, 0.8, 100 / 101),
       labels = c("1:100", "1:4", "2:3", "3:2", "4:1", "100:1"), line = 3.2, cex.axis = 0.85)
  mtext("Cost:Benefit Ratio", side = 1, line = 5.4)
  legend("topright", c("Number high risk", sprintf("Number high risk with %s", event_label), "95% CI (bootstrap)"),
         col = c("#C0392B", "#2E86C1", "grey40"), lwd = c(2.5, 2.5, 0.8), lty = c(1, 2, 1), bty = "n")
}

nom_impact_df <- function(impact) {
  data.frame(
    threshold = impact$threshold,
    n_high_risk = round(impact$nhigh, 1),
    n_high_risk_lo95 = round(impact$nhigh_lo, 1),
    n_high_risk_hi95 = round(impact$nhigh_hi, 1),
    n_high_risk_true_case = round(impact$nevent, 1),
    n_high_risk_true_case_lo95 = round(impact$nevent_lo, 1),
    n_high_risk_true_case_hi95 = round(impact$nevent_hi, 1)
  )
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------

mod_nomogram_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$style(HTML("
      .nom-module .small-box h3 { font-size: 20px; white-space: normal; overflow-wrap: break-word; line-height: 1.2; }
      .nom-module .source-picker-grid .btn-group-justified { display: grid; grid-template-columns: 1fr 1fr; gap: 6px; }
      .nom-module .source-picker-grid .btn-group-justified > .btn-group { display: block; width: 100%; }
      .nom-module .source-picker-grid .btn-group-justified .btn { width: 100%; white-space: normal; font-size: 12px; line-height: 1.3; padding: 8px 6px; }
      .nom-module .preset-panel { background: #F8FAFC; border: 1px solid var(--color-border); border-radius: var(--radius-card); padding: 14px 16px; margin-bottom: 8px; }
      .nom-module .gene-chip-row { display: flex; flex-wrap: wrap; gap: 6px; margin: 10px 0; }
      .nom-module .gene-chip { background: var(--color-surface-blue); color: var(--color-primary-dark); font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 12px; font-weight: 600; padding: 4px 9px; border-radius: 20px; line-height: 1.4; }
    ")),
    div(class = "nom-module",
      fluidRow(
        column(
          4,
          arthochat_shortcut_ui(
            "Questions about this nomogram? Ask ArthOChat.",
            compact = TRUE
          ),
          box(
            width = NULL, title = "Gene panel & cohort", status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "Fits a per-sex model, then builds a nomogram, calibration curve, decision curve and clinical-impact curve."),
            div(class = "filter-section-header", "SEX"),
            shinyWidgets::radioGroupButtons(
              ns("fit_sex"), NULL, choices = c("Female" = "female", "Male" = "male"),
              selected = "female", status = "default", size = "sm", justified = TRUE
            ),
            div(class = "filter-section-header", style = "margin-top:10px;", "GENE SOURCE"),
            div(class = "source-picker-grid",
              shinyWidgets::radioGroupButtons(
                ns("gene_source"), NULL,
                choices = c("Your own" = "own", "Same-tissue" = "same_tissue",
                            "Cross-tissue" = "cross_tissue", "Cross-ancestry" = "cross_ancestry"),
                selected = "same_tissue", status = "default", size = "sm", justified = TRUE
              )
            ),
            uiOutput(ns("preset_panel_ui")),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'own'", ns("gene_source")),
              radioButtons(
                ns("own_cohort"), "Fit & evaluate on",
                choiceNames = list("Same-tissue (blood, current dataset)", "Cross-tissue (synovium, GSE89408)"),
                choiceValues = list("blood", "synovium"), selected = "blood"
              ),
              textAreaInput(ns("gene_list"), "Gene panel (2-12 genes; this project's own panels use 4-7)", rows = 5, placeholder = "TNF\nIL6\nSTAT3\n...")
            ),
            tags$hr(),
            uiOutput(ns("cohort_contrast_ui")),
            actionButton(ns("run_btn"), "Build nomogram", icon = icon("play"), class = "btn-primary btn-sm")
          ),
          box(
            width = NULL, title = "Advanced filters", status = "primary", solidHeader = FALSE,
            strong("Ridge penalty"),
            radioButtons(
              ns("penalty_mode"), NULL,
              choiceNames = list(
                tagList("Project default ", tags$small("(0 female / 5 male - stabilises the smaller, near-separated male fit)")),
                "Custom"
              ),
              choiceValues = list("project", "manual"), selected = "project"
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'manual'", ns("penalty_mode")),
              numericInput(ns("penalty_manual"), "Penalty value", value = 0, min = 0, max = 30, step = 0.5)
            ),
            tags$hr(),
            numericInput(ns("calibrate_B"), "Calibration bootstrap reps", value = 200, min = 50, max = 1000, step = 50),
            numericInput(ns("dca_step"), "Decision-curve threshold step", value = 0.01, min = 0.01, max = 0.05, step = 0.01),
            selectInput(ns("impact_N"), "Clinical impact: per how many patients", choices = c(100, 1000, 10000), selected = 1000),
            numericInput(ns("impact_B"), "Clinical impact bootstrap reps", value = 500, min = 50, max = 1000, step = 50),
            numericInput(ns("seed"), "Random seed", value = 1234, min = 1, step = 1),
            div(class = "empty-note", style = "font-size: 12.5px;", icon("circle-info"),
                "Nothing recomputes until “Build nomogram” is clicked again.")
          )
        ),
        column(
          8,
          ## Nothing in this column is rendered client-side until "Build
          ## nomogram" has been clicked at least once - `nom_result()` is an
          ## eventReactive keyed on that button (ignoreInit = TRUE), and
          ## reading it beforehand from a bare plotOutput() would otherwise
          ## surface as a raw dispatch error rather than a clean empty
          ## state (the same reason mod_crosstissue.R's per-sex panels are
          ## wrapped in an identical `input['...'] > 0` conditionalPanel).
          conditionalPanel(
            condition = sprintf("input['%s'] > 0", ns("run_btn")),
            tabsetPanel(
              id = ns("result_tabs"), type = "tabs",
              tabPanel(
                "Nomogram", br(),
                box(width = NULL, status = "primary", solidHeader = FALSE,
                    p(class = "submodule-desc", "Read a patient's value on each gene axis, sum the points, and project the total onto Risk."),
                    withSpinner(plotOutput(ns("nomogram_plot"), height = 460), color = "#2c6fbb", type = 6))
              ),
              tabPanel(
                "Calibration", br(),
                box(width = NULL, status = "primary", solidHeader = FALSE,
                    p(class = "submodule-desc", "Apparent vs bootstrap bias-corrected calibration - agreement between predicted and observed probability."),
                    withSpinner(plotOutput(ns("calibration_plot"), height = 460), color = "#2c6fbb", type = 6))
              ),
              tabPanel(
                "Decision Curve", br(),
                box(width = NULL, status = "primary", solidHeader = FALSE,
                    p(class = "submodule-desc", "Net clinical benefit of using the panel vs treating everyone vs treating no one, across threshold probabilities (Vickers & Elkin, 2006)."),
                    withSpinner(plotOutput(ns("dca_plot"), height = 440), color = "#2c6fbb", type = 6),
                    div(class = "table-toolbar", downloadButton(ns("download_dca"), "Download CSV", class = "btn-sm")),
                    DT::dataTableOutput(ns("dca_table")))
              ),
              tabPanel(
                "Clinical Impact", br(),
                box(width = NULL, status = "primary", solidHeader = FALSE,
                    p(class = "submodule-desc", "Per-N patients: how many would be flagged high-risk, and how many of those genuinely have the disease, with 95% bootstrap bands (Kerr et al., 2016)."),
                    withSpinner(plotOutput(ns("impact_plot"), height = 440), color = "#2c6fbb", type = 6),
                    div(class = "table-toolbar", downloadButton(ns("download_impact"), "Download CSV", class = "btn-sm")),
                    DT::dataTableOutput(ns("impact_table")))
              )
            ),
            box(width = NULL, title = "Result", status = "primary", solidHeader = FALSE,
                withSpinner(uiOutput(ns("summary_ui")), color = "#2c6fbb", type = 6))
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 0", ns("run_btn")),
            box(width = NULL, status = "primary", solidHeader = FALSE,
                div(class = "empty-note", icon("circle-info"),
                    "Pick a sex and gene source on the left, then click “Build nomogram”."))
          )
        )
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] > 0", ns("run_btn")),
        box(
          width = 12, title = "Model coefficients", status = "primary", solidHeader = FALSE,
          div(class = "table-toolbar",
              downloadButton(ns("download_nom"), "Coefficients (CSV)", class = "btn-sm"),
              downloadButton(ns("download_model"), "Fitted model (.rds)", class = "btn-sm")),
          DT::dataTableOutput(ns("coef_table"))
        )
      )
    )
  )
}

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------

mod_nomogram_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Independent synovium cohort for the Cross-tissue source/cohort -
    ## loaded once, read-only, same object mod_crosstissue.R uses.
    val <- tryCatch(readRDS(VAL_SYNOVIUM_RDS), error = function(e) NULL)
    bundled_venn <- read_table_safe("FS_venn_membership.csv")

    panel_labels <- c(
      same_tissue = "Same-tissue panel (blood)",
      cross_tissue = "Cross-tissue panel (blood → synovium)",
      cross_ancestry = "Cross-ancestry panel (MR-replicated)"
    )
    panel_blurbs <- c(
      same_tissue = "This session's Feature Selection consensus panel (or this project's bundled one) - fit and evaluated in blood.",
      cross_tissue = "The blood-derived panel, fit and evaluated fresh in the independent synovial cohort (GSE89408) - never a transfer of blood-fitted coefficients.",
      cross_ancestry = "Genes whose MR causal signal replicates in European AND transfers to East-Asian ancestry GWAS - fit in blood, since this project has no ancestry-stratified expression cohort."
    )

    ## -----------------------------------------------------------------
    ## Gene-panel sources - same live-session-first, bundled-fallback
    ## convention as mod_enrichment.R's same_tissue_panel()/cross_tissue_
    ## panel()/cross_ancestry_panel() (and mod_crosstissue.R's project
    ## panel getter), reused here rather than re-derived.
    ## -----------------------------------------------------------------

    nom_same_tissue_panel <- function(sex_label) {
      live <- results$featureselection[[sex_label]]$consensus_genes
      if (!is.null(live) && length(live) >= 2) {
        return(list(genes = live, is_live = TRUE,
          note = sprintf("%d genes - this session's live Feature Selection %s consensus panel.", length(live), sex_label)))
      }
      if (!is.null(bundled_venn)) {
        genes <- unique(bundled_venn$gene[
          bundled_venn$candidate_set == "noMHC" & bundled_venn$sex == sex_label &
          as.character(bundled_venn$in_consensus) == "TRUE"
        ])
        if (length(genes) >= 2) {
          return(list(genes = genes, is_live = FALSE,
            note = sprintf("%d genes - this project's bundled %s diagnostic consensus panel (MHC-excluded).", length(genes), sex_label)))
        }
      }
      list(genes = character(0), is_live = FALSE, note = sprintf("No %s panel available.", sex_label))
    }

    nom_cross_tissue_panel <- function(sex_label) {
      live <- results$featureselection[[sex_label]]$consensus_genes
      if (!is.null(live) && length(live) >= 2) {
        return(list(genes = live, is_live = TRUE,
          note = sprintf("%d genes - this session's live Feature Selection %s consensus panel, evaluated in synovium.", length(live), sex_label)))
      }
      bundled <- if (!is.null(val)) (if (identical(sex_label, "female")) val$fsig else val$msig) else NULL
      bundled <- unique(as.character(bundled))
      if (length(bundled) >= 2) {
        return(list(genes = bundled, is_live = FALSE,
          note = sprintf("%d genes - this project's bundled %s panel, the one its own synovial validation used (MHC-retained).", length(bundled), sex_label)))
      }
      list(genes = character(0), is_live = FALSE, note = sprintf("No %s panel available.", sex_label))
    }

    nom_cross_ancestry_panel <- function(sex_label) {
      live <- results$crossancestry[[sex_label]]$biomarker_genes
      if (!is.null(live) && length(live) >= 2) {
        return(list(genes = live, is_live = TRUE,
          note = sprintf("%d genes - this session's live Cross-Ancestry replicated-biomarker list.", length(live))))
      }
      bundled <- read_table_safe(sprintf("MR35_crossancestry_%s.csv", sex_label))
      if (!is.null(bundled)) {
        genes <- unique(bundled$gene[bundled$ancestry_class == "shared EUR+EAS"])
        if (length(genes) >= 2) {
          return(list(genes = genes, is_live = FALSE,
            note = sprintf("%d genes - replicated in European AND transferable to East-Asian ancestry GWAS (this project's bundled %s cross-ancestry table).", length(genes), sex_label)))
        }
      }
      list(genes = character(0), is_live = FALSE, note = sprintf("No %s panel available.", sex_label))
    }

    panel_fns <- list(same_tissue = nom_same_tissue_panel, cross_tissue = nom_cross_tissue_panel, cross_ancestry = nom_cross_ancestry_panel)

    output$preset_panel_ui <- renderUI({
      src <- input$gene_source %||% "own"
      if (identical(src, "own")) return(NULL)
      sex_label <- input$fit_sex %||% "female"
      panel <- panel_fns[[src]](sex_label)
      div(class = "preset-panel",
        div(class = "card-title", style = "font-size: 13.5px;", panel_labels[[src]]),
        div(class = "card-subtitle", panel_blurbs[[src]]),
        if (!length(panel$genes)) {
          div(class = "empty-note", style = "margin-top: 10px;", icon("circle-info"), panel$note)
        } else {
          tagList(
            div(class = "gene-chip-row", lapply(panel$genes, function(g) span(class = "gene-chip", g))),
            div(class = "empty-note", style = "margin-top: 0;",
                icon(if (isTRUE(panel$is_live)) "check" else "circle-info"), panel$note)
          )
        }
      )
    })

    ## Which expression cohort a click on "Build nomogram" will actually
    ## fit on: cross_tissue always means synovium; same_tissue and
    ## cross_ancestry always mean blood (see module header); own lets the
    ## user pick either.
    current_cohort <- reactive({
      src <- input$gene_source %||% "own"
      if (identical(src, "cross_tissue")) "synovium"
      else if (identical(src, "own")) (input$own_cohort %||% "blood")
      else "blood"
    })

    output$cohort_contrast_ui <- renderUI({
      if (identical(current_cohort(), "synovium")) {
        return(div(class = "empty-note", icon("circle-info"),
          "Synovium contrast is fixed: Normal (reference) vs RA - GSE89408 has only these two groups."))
      }
      sex_label <- input$fit_sex %||% "female"
      sex_code <- if (identical(sex_label, "female")) "F" else "M"
      meta <- dataset$meta
      if (!("sex" %in% colnames(meta))) {
        return(div(class = "empty-note", icon("triangle-exclamation"), "The loaded metadata has no 'sex' column - required for this per-sex nomogram."))
      }
      groups <- sort(unique(na.omit(meta$group[meta$sex == sex_code])))
      validate(need(length(groups) >= 2, sprintf("The %s subset of the loaded dataset needs at least two group values.", sex_label)))
      ref_default <- if ("HC" %in% groups) "HC" else groups[1]
      comp_default <- if ("RA" %in% groups) "RA" else groups[groups != ref_default][1]
      tagList(
        selectInput(ns("ref_group"), "Reference group (0)", choices = groups, selected = ref_default, selectize = FALSE),
        selectInput(ns("comp_group"), "Comparison group (1)", choices = groups, selected = comp_default, selectize = FALSE),
        if ("age" %in% colnames(meta)) checkboxInput(ns("include_age"), "Include age as a covariate", value = FALSE)
      )
    })

    nom_advanced_params <- function() {
      list(
        penalty_mode = input$penalty_mode %||% NOM_DEFAULT_PARAMS$penalty_mode,
        penalty_manual = input$penalty_manual %||% NOM_DEFAULT_PARAMS$penalty_manual,
        calibrate_B = input$calibrate_B %||% NOM_DEFAULT_PARAMS$calibrate_B,
        dca_step = input$dca_step %||% NOM_DEFAULT_PARAMS$dca_step,
        impact_N = suppressWarnings(as.numeric(input$impact_N)) %||% NOM_DEFAULT_PARAMS$impact_N,
        impact_B = input$impact_B %||% NOM_DEFAULT_PARAMS$impact_B,
        seed = input$seed %||% NOM_DEFAULT_PARAMS$seed
      )
    }

    ## -----------------------------------------------------------------
    ## Run: locks in the current sex/source/cohort/filters and builds the
    ## full instrument. Nothing above recomputes until this fires again.
    ## -----------------------------------------------------------------

    nom_result <- eventReactive(input$run_btn, {
      sex_label <- input$fit_sex %||% "female"
      sex_code <- if (identical(sex_label, "female")) "F" else "M"
      src <- input$gene_source %||% "own"
      cohort <- current_cohort()
      params <- nom_advanced_params()

      genes_req <- if (identical(src, "own")) {
        g <- unique(trimws(unlist(strsplit(input$gene_list %||% "", "[,\n\t ]+"))))
        g[nzchar(g)]
      } else {
        panel_fns[[src]](sex_label)$genes
      }
      validate(need(length(genes_req) >= 2, "At least 2 genes are needed to build a nomogram."))
      validate(need(length(genes_req) <= 12, "Use 12 or fewer genes; a nomogram with many predictors becomes unreadable (this project's own panels use 4-7)."))

      if (identical(cohort, "synovium")) {
        validate(need(!is.null(val), "The bundled synovium validation object could not be loaded."))
        idx_sex <- which(val$sex == sex_code)
        validate(need(length(idx_sex) > 0, sprintf("No %s samples in the synovium dataset.", sex_label)))
        genes_present <- intersect(genes_req, rownames(val$logcpm))
        validate(need(length(genes_present) >= 2, sprintf("Fewer than 2 panel genes are present in the synovium dataset for %s.", sex_label)))
        y_full <- droplevels(val$grp[idx_sex])
        validate(need(length(unique(y_full)) == 2 && all(table(y_full) >= 10),
                      sprintf("The %s synovium subset needs at least 10 samples in each group (RA/Normal).", sex_label)))
        expr_sub <- val$logcpm[genes_present, idx_sex, drop = FALSE]
        df <- as.data.frame(t(expr_sub)); names(df) <- make.names(names(df))
        df$y <- as.numeric(y_full) - 1  ## factor levels "Normal","RA" -> Normal = 0, RA = 1
        event_label <- "RA"
        tissue_label <- "Cross-tissue — synovium (GSE89408)"
      } else {
        meta <- dataset$meta
        validate(need("sex" %in% colnames(meta), "The loaded dataset's metadata has no 'sex' column - required for this per-sex nomogram."))
        meta_sex <- meta[!is.na(meta$sex) & meta$sex == sex_code, , drop = FALSE]
        validate(need(nrow(meta_sex) > 0, sprintf("No %s samples in the currently loaded dataset.", sex_label)))
        groups <- sort(unique(na.omit(meta_sex$group)))
        validate(need(length(groups) >= 2, sprintf("The %s subset needs at least two group values.", sex_label)))
        ref_group <- if (!is.null(input$ref_group) && input$ref_group %in% groups) input$ref_group else if ("HC" %in% groups) "HC" else groups[1]
        comp_group <- if (!is.null(input$comp_group) && input$comp_group %in% groups) input$comp_group else if ("RA" %in% groups) "RA" else groups[groups != ref_group][1]
        validate(need(ref_group != comp_group, "Reference and comparison group must be different."))
        meta_sub <- meta_sex[meta_sex$group %in% c(ref_group, comp_group), , drop = FALSE]
        genes_present <- intersect(genes_req, rownames(dataset$expr))
        validate(need(length(genes_present) >= 2, sprintf("Fewer than 2 panel genes are present in the current expression matrix for %s.", sex_label)))
        common <- intersect(colnames(dataset$expr), meta_sub$sample)
        meta_sub <- meta_sub[match(common, meta_sub$sample), , drop = FALSE]
        expr_sub <- dataset$expr[genes_present, common, drop = FALSE]
        df <- as.data.frame(t(expr_sub)); names(df) <- make.names(names(df))
        df$y <- as.numeric(factor(meta_sub$group, levels = c(ref_group, comp_group))) - 1
        if (isTRUE(input$include_age) && "age" %in% colnames(meta_sub)) df$age <- meta_sub$age
        df <- df[stats::complete.cases(df), , drop = FALSE]
        validate(need(nrow(df) >= 30, sprintf("Fewer than 30 complete %s samples for this contrast (after dropping missing age, if included).", sex_label)))
        validate(need(all(table(df$y) >= 10), sprintf("Each group needs at least 10 complete %s samples.", sex_label)))
        event_label <- comp_group
        tissue_label <- if (identical(src, "cross_ancestry")) {
          "Cross-ancestry panel, fit on same-tissue blood (this project has no ancestry-stratified expression cohort)"
        } else {
          "Same-tissue — blood (current dataset)"
        }
      }

      predictors <- setdiff(names(df), "y")
      penalty <- if (identical(params$penalty_mode, "manual")) params$penalty_manual else if (identical(sex_label, "male")) 5 else 0

      res <- nom_fit_core(df, predictors, penalty, event_label, params)
      res$sex_label <- sex_label; res$src <- src; res$cohort <- cohort; res$tissue_label <- tissue_label
      res$n_input <- length(genes_req); res$n_present <- length(genes_present)
      res$age_included <- "age" %in% predictors
      res
    }, ignoreInit = TRUE)

    output$summary_ui <- renderUI({
      res <- nom_result()
      tagList(
        p(strong(res$n_present), " of ", res$n_input, " requested genes",
          if (isTRUE(res$age_included)) " + age" else "", " used — ",
          tools::toTitleCase(res$sex_label), ", ", res$tissue_label, "."),
        p(strong(res$n_samples), " complete samples (", strong(res$n_pos), " ", res$event_label,
          " / ", strong(res$n_neg), " reference)."),
        p("Model C-statistic (discrimination): ", strong(sprintf("%.3f", res$c_stat)),
          " · Ridge penalty used: ", strong(res$penalty)),
        if (identical(res$src, "cross_ancestry")) {
          div(class = "empty-note", icon("circle-info"),
              "Ancestry replication here is genetic (Mendelian randomisation across European + East-Asian GWAS); the model above is still fit on this project's blood cohort - see Cross-Ancestry Validation.")
        }
      )
    })

    output$nomogram_plot <- renderPlot({
      plot(nom_result()$nom)
    })

    output$calibration_plot <- renderPlot({
      res <- nom_result()
      if (is.null(res$cal)) {
        plot.new()
        text(0.5, 0.5, "Calibration bootstrap did not converge for this fit.\nTry fewer genes, more samples, or a larger penalty.", cex = 1.05)
      } else {
        plot(res$cal, xlab = "Predicted probability", ylab = "Actual probability", subtitles = FALSE, main = "")
      }
    })

    output$dca_plot <- renderPlot({
      res <- nom_result()
      nom_render_dca_plot(res$dca_df, res$event_label)
    })

    output$impact_plot <- renderPlot({
      res <- nom_result()
      nom_render_impact_plot(res$impact, res$event_label)
    })

    output$dca_table <- DT::renderDataTable({
      d <- nom_result()$dca_df
      d[c("NB_panel", "NB_all", "NB_none")] <- lapply(d[c("NB_panel", "NB_all", "NB_none")], round, 4)
      DT::datatable(d, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })

    output$impact_table <- DT::renderDataTable({
      DT::datatable(nom_impact_df(nom_result()$impact), rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
    })

    output$coef_table <- DT::renderDataTable({
      DT::datatable(nom_result()$coef_df, rownames = FALSE, options = list(pageLength = 10, dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_nom <- downloadHandler(
      filename = function() sprintf("nomogram_%s_coefficients.csv", nom_result()$sex_label),
      content = function(file) write.csv(nom_result()$coef_df, file, row.names = FALSE)
    )
    output$download_dca <- downloadHandler(
      filename = function() sprintf("nomogram_%s_decision_curve.csv", nom_result()$sex_label),
      content = function(file) write.csv(nom_result()$dca_df, file, row.names = FALSE)
    )
    output$download_impact <- downloadHandler(
      filename = function() sprintf("nomogram_%s_clinical_impact.csv", nom_result()$sex_label),
      content = function(file) write.csv(nom_impact_df(nom_result()$impact), file, row.names = FALSE)
    )
    output$download_model <- downloadHandler(
      filename = function() sprintf("nomogram_%s_model.rds", nom_result()$sex_label),
      content = function(file) saveRDS(nom_result()$fit, file)
    )

    ## Published for ArthOChat / build_assistant_context(), same
    ## per-sex-merge convention every other Validation submodule uses.
    observeEvent(nom_result(), {
      res <- nom_result()
      results$nomogram <- utils::modifyList(
        results$nomogram %||% list(),
        setNames(list(list(
          n_genes = res$n_present, n_samples = res$n_samples, c_stat = round(res$c_stat, 4),
          penalty = res$penalty, tissue = res$tissue_label, gene_source = res$src
        )), res$sex_label)
      )
    })
  })
}
