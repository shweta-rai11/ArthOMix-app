## Feature Selection module: fits LASSO, random forest, and SVM-RFE
## independently per sex on a candidate gene panel, then reports the
## consensus overlap. Female/male are always modeled separately, never
## pooled with sex as a covariate. Candidate genes come from live Candidate
## Gene Identification output, an uploaded expression matrix, or an
## uploaded DEG/gene list (radioButtons "data_source": project/expr/deg).

## SVM-RFE: recursive elimination by smallest squared linear-SVM weight,
## one feature per round, then a CV error curve over the ranking to pick
## panel size objectively.

## CSV formula-injection guard: gene/feature identifiers written into
## downloadHandler CSVs can come from a user upload or GEO metadata, not
## just this app's own curated panels. A symbol starting with =, +, -, or @
## is interpreted as a formula by Excel/Sheets on open; prefixing a single
## quote neutralizes the leading character while leaving the value's own
## text intact everywhere else (including when re-read by read.csv()).
tx_csv_safe <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^[=+\\-@]", x), paste0("'", x), x)
}

## Linear kernel only - the per-round elimination weight (squared SVM
## coefficient) is only meaningful for a linear decision boundary.
fs_svm_rfe_rank <- function(X, y, cost = 1, tolerance = 0.001, class_weights = NULL) {
  ## Zero-variance-column guard (mirrors mod_diagnostic.R's diag_zrows()
  ## convention): e1071::svm(..., scale = TRUE) internally z-scores each
  ## column via scale(), which divides by the column SD - a constant
  ## candidate column (SD 0 or NA, e.g. a gene with identical expression
  ## across every sample after upstream filtering) turns into NaN for every
  ## row, silently corrupting the fit and the resulting elimination weights.
  ## Such columns carry no discriminative signal for a linear SVM anyway, so
  ## they're dropped before ranking starts, then appended back at the end of
  ## the ranking (least-informative first) so every input feature is still
  ## accounted for in the returned order.
  zero_var <- vapply(as.data.frame(X), function(col) {
    s <- stats::sd(col, na.rm = TRUE)
    is.na(s) || s == 0
  }, logical(1))
  dropped <- colnames(X)[zero_var]
  X <- X[, !zero_var, drop = FALSE]

  feats <- colnames(X)
  ranking <- character(0)
  while (length(feats) > 1) {
    m <- e1071::svm(X[, feats, drop = FALSE], y, kernel = "linear", scale = TRUE, cost = cost,
                     tolerance = tolerance, class.weights = class_weights)
    w2 <- ((t(m$coefs) %*% m$SV)[1, ])^2
    drop <- names(sort(w2))[1]
    ranking <- c(drop, ranking)
    feats <- setdiff(feats, drop)
  }
  c(feats, ranking, dropped)
}

fs_svm_rfe_curve <- function(X, y, rank, cost = 1, seed = 1234, folds = 10, tolerance = 0.001, class_weights = NULL) {
  ks <- seq_along(rank)
  err <- vapply(ks, function(k) {
    set.seed(seed)
    Xk <- X[, rank[seq_len(k)], drop = FALSE]
    ## Same zero-variance-column guard as fs_svm_rfe_rank() above: rank's
    ## own zero-variance columns are appended at its tail, so a k large
    ## enough to reach them would otherwise feed a constant column straight
    ## into e1071::svm(..., scale = TRUE) here too, producing NaN.
    keep <- vapply(as.data.frame(Xk), function(col) {
      s <- stats::sd(col, na.rm = TRUE)
      !(is.na(s) || s == 0)
    }, logical(1))
    if (!any(keep)) return(NA_real_)
    acc <- e1071::svm(Xk[, keep, drop = FALSE], y, kernel = "linear",
                       scale = TRUE, cost = cost, cross = folds,
                       tolerance = tolerance, class.weights = class_weights)$tot.accuracy
    1 - acc / 100
  }, numeric(1))
  list(k = ks, err = err, best = ks[which.min(err)], besterr = min(err))
}

FS_SVM_COST_GRID <- c(0.01, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16)

## SVM-RFE refits once per eliminated feature (~2p fits for p genes) - fine
## for a 25-40 gene MR panel but hours-long on a raw WGCNA module (1000s of
## genes). Candidate sets above this cap get reduced to their most variable
## genes first, in fs_build_sex() below, with a note shown on-screen.
FS_MAX_CANDIDATE_GENES <- 200

## Minimum samples required in EACH group (e.g. HC and RA, within one sex) for
## LASSO/RF/SVM-RFE to fit at all - checked in fs_build_sex() below. fs_fit_sex()
## already scales its own CV fold count down to whatever's actually available
## (`nf <- max(2, min(cv_folds, min_class))`), so this floor exists only to stop
## a fit degenerating to a near-meaningless 2-fold split (1 sample per class per
## fold) - at 4, every method still gets at least a real 4-fold CV.
FS_MIN_GROUP_SAMPLES <- 4

## A second, higher bar purely for the UI caveat below FS_MIN_GROUP_SAMPLES - a
## fit between the two thresholds runs (and its CV fold count is still valid),
## but with fewer than this many samples in its smallest group, results should
## be read as exploratory rather than a fit anyone should treat as final.
FS_RELIABLE_GROUP_SAMPLES <- 6

## Default hyperparameters, matching this project's own analysis script;
## overridable per-method from the UI (see fs_fit_sex()'s `params`).
FS_DEFAULT_PARAMS <- list(
  # equal = unweighted (default); balanced = inverse-frequency; manual = fixed ratio
  class_weight_mode = "equal", class_weight_ratio = 1,
  lasso_cv_folds = 10, lasso_alpha = 1, lasso_lambda_choice = "lambda.min",
  lasso_nlambda = 100, lasso_type_measure = "deviance",
  rf_cv_folds = 10, rf_ntree = 1000, rf_mtry_mode = "auto", rf_mtry_manual = NULL,
  rf_selection_rule = "above_mean", rf_top_n = 10,
  rf_nodesize = 1, rf_maxnodes = NULL,
  svm_cv_folds = 10, svm_cost_mode = "auto", svm_cost_manual = 1, svm_cost_grid = FS_SVM_COST_GRID,
  svm_panel_mode = "auto", svm_manual_k = 10, svm_tolerance = 0.001,
  # methods intersected for consensus; e.g. drop RandomForest to replicate Chen et al. 2021/2022's LASSO∩SVM-RFE panel
  consensus_methods = c("LASSO", "RandomForest", "SVM_RFE")
)

## Class weights by level for imbalanced groups; feeds randomForest's classwt
## and e1071's class.weights directly, and glmnet's weights= via fs_obs_weights().
fs_class_weight_levels <- function(y, mode, ratio) {
  lv <- levels(y)
  if (identical(mode, "balanced")) {
    n <- table(y)
    w <- max(n) / n
    stats::setNames(as.numeric(w[lv]), lv)
  } else if (identical(mode, "manual")) {
    stats::setNames(c(1, ratio %||% 1), lv)
  } else {
    stats::setNames(c(1, 1), lv)
  }
}

fs_obs_weights <- function(y, mode, ratio) {
  wl <- fs_class_weight_levels(y, mode, ratio)
  unname(wl[as.character(y)])
}

## Fits LASSO + tuned random forest + tuned SVM-RFE on one sex's data.
## Gene symbols are made syntactically safe (e.g. HLA-A) for the fits and
## translated back on the way out.
fs_fit_sex <- function(X, y, params = list()) {
  ## caret::train(classProbs = TRUE) below requires factor levels that are
  ## valid R variable names - it make.names()s them internally to build its
  ## own predicted-probability column names, so a raw group label with a
  ## space (e.g. "multiple sclerosis") desyncs from any levels(y)-based
  ## lookup once caret has already renamed its own columns to
  ## "multiple.sclerosis". Sanitized once here, up front, matching
  ## mod_diagnostic.R::diag_fit_sex()'s identical fix; callers keep the real
  ## group names for their own display text, so nothing user-visible changes.
  levels(y) <- make.names(levels(y), unique = TRUE)
  params <- utils::modifyList(FS_DEFAULT_PARAMS, params)
  GLOBAL_SEED <- ARTHOMIX_TX_ML_SEED
  p <- ncol(X)
  min_class <- min(table(y))
  nf_lasso <- max(2, min(params$lasso_cv_folds, min_class))
  nf_rf <- max(2, min(params$rf_cv_folds, min_class))
  nf_svm <- max(2, min(params$svm_cv_folds, min_class))

  safe <- make.names(colnames(X), unique = TRUE)
  lk <- stats::setNames(colnames(X), safe)
  colnames(X) <- safe
  back <- function(v) unname(lk[v])

  # shared class weights for all three methods below; "equal" mode is a no-op (all weights 1)
  cw_levels <- fs_class_weight_levels(y, params$class_weight_mode, params$class_weight_ratio)
  obs_w <- fs_obs_weights(y, params$class_weight_mode, params$class_weight_ratio)

  # (1) LASSO logistic regression - non-zero coefficients at the chosen lambda
  set.seed(GLOBAL_SEED)
  cv <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = params$lasso_alpha, nfolds = nf_lasso,
                           type.measure = params$lasso_type_measure %||% "deviance", nlambda = params$lasso_nlambda,
                           weights = obs_w)
  lambda_s <- if (identical(params$lasso_lambda_choice, "lambda.1se")) "lambda.1se" else "lambda.min"
  co <- coef(cv, s = lambda_s)[-1, 1, drop = TRUE]
  lasso_genes <- back(names(co)[co != 0])

  # (2) Random forest importance - mtry CV-tuned (or manual), genes kept above mean Gini (or top-N)
  ntree <- max(100, round(params$rf_ntree))
  rf_nodesize <- max(1, round(params$rf_nodesize %||% 1))
  rf_maxnodes <- if (!is.null(params$rf_maxnodes) && is.finite(params$rf_maxnodes)) max(2, round(params$rf_maxnodes)) else NULL
  if (identical(params$rf_mtry_mode, "manual") && !is.null(params$rf_mtry_manual)) {
    best_mtry <- min(p, max(1, round(params$rf_mtry_manual)))
  } else {
    mtry_grid <- sort(unique(pmin(p, c(1, 2, floor(sqrt(p)), floor(p / 3), floor(p / 2), p))))
    ctrl <- caret::trainControl(method = "cv", number = nf_rf, classProbs = TRUE, summaryFunction = caret::twoClassSummary)
    set.seed(GLOBAL_SEED)
    rf_tune <- tryCatch(
      caret::train(x = X, y = y, method = "rf", metric = "ROC", trControl = ctrl,
                   tuneGrid = expand.grid(mtry = mtry_grid), ntree = ntree, importance = TRUE,
                   nodesize = rf_nodesize, maxnodes = rf_maxnodes, classwt = cw_levels),
      error = function(e) NULL
    )
    best_mtry <- if (!is.null(rf_tune)) rf_tune$bestTune$mtry else max(1, floor(sqrt(p)))
  }
  set.seed(GLOBAL_SEED)
  rf <- randomForest::randomForest(X, y, importance = TRUE, ntree = ntree, mtry = best_mtry,
                                    nodesize = rf_nodesize, maxnodes = rf_maxnodes, classwt = cw_levels)
  gini <- sort(rf$importance[, "MeanDecreaseGini"], decreasing = TRUE)
  gini_thr <- mean(gini)
  rf_genes <- if (identical(params$rf_selection_rule, "top_n")) {
    back(head(names(gini), min(p, max(1, round(params$rf_top_n)))))
  } else {
    back(names(gini)[gini > gini_thr])
  }
  names(gini) <- back(names(gini))

  # (3) SVM-RFE, linear kernel - cost CV-tuned (or manual), panel size = CV-optimal k (or manual)
  svm_tolerance <- params$svm_tolerance %||% 0.001
  if (identical(params$svm_cost_mode, "manual") && !is.null(params$svm_cost_manual)) {
    best_cost <- params$svm_cost_manual
  } else {
    grid <- params$svm_cost_grid
    if (!is.numeric(grid) || length(grid) == 0) grid <- FS_SVM_COST_GRID
    set.seed(GLOBAL_SEED)
    svm_tune <- tryCatch(
      e1071::tune(e1071::svm, train.x = X, train.y = y, kernel = "linear", scale = TRUE,
                  tolerance = svm_tolerance, class.weights = cw_levels,
                  ranges = list(cost = grid),
                  tunecontrol = e1071::tune.control(sampling = "cross", cross = nf_svm)),
      error = function(e) NULL
    )
    best_cost <- if (!is.null(svm_tune)) svm_tune$best.parameters$cost else 1
  }
  set.seed(GLOBAL_SEED)
  rank <- fs_svm_rfe_rank(X, y, cost = best_cost, tolerance = svm_tolerance, class_weights = cw_levels)
  curve <- fs_svm_rfe_curve(X, y, rank, cost = best_cost, seed = GLOBAL_SEED, folds = nf_svm, tolerance = svm_tolerance, class_weights = cw_levels)
  svm_rank <- back(rank)
  svm_genes <- if (identical(params$svm_panel_mode, "manual") && !is.null(params$svm_manual_k)) {
    head(svm_rank, min(length(svm_rank), max(1, round(params$svm_manual_k))))
  } else {
    svm_rank[seq_len(curve$best)]
  }

  sets <- list(LASSO = lasso_genes, RandomForest = rf_genes, SVM_RFE = svm_genes)

  # consensus = intersection of the selected methods; falls back to all three if none selected
  consensus_methods <- intersect(params$consensus_methods %||% names(sets), names(sets))
  if (length(consensus_methods) == 0) consensus_methods <- names(sets)
  consensus_genes <- Reduce(intersect, sets[consensus_methods])

  list(
    cv = cv, lasso_genes = lasso_genes, lasso_alpha = params$lasso_alpha, lasso_lambda_choice = lambda_s,
    gini = gini, gini_thr = gini_thr, rf_mtry = best_mtry, rf_ntree = ntree, rf_genes = rf_genes,
    rf_selection_rule = params$rf_selection_rule, rf_top_n = params$rf_top_n,
    svm_rank = svm_rank, svm_curve = curve, svm_cost = best_cost, svm_genes = svm_genes,
    svm_panel_mode = params$svm_panel_mode, svm_manual_k = params$svm_manual_k,
    sets = sets, consensus = consensus_genes, consensus_methods = consensus_methods,
    n_input = p, n_samples = nrow(X), fast_path = FALSE
  )
}

mod_featureselection_config <- list(
  id = "featureselection", group = "Biomarker modeling",
  title = "Feature Selection",
  description = "Select features based on LASSO, random forest and SVM-RFE by sex.",
  icon = "sliders"
)

## Reveals a sex's technique box only after that sex's own Run button has been clicked.
fs_sex_panel <- function(ns, run_btn_id, ...) {
  conditionalPanel(condition = sprintf("input['%s'] > 0", ns(run_btn_id)), ...)
}

## One technique's result box (summary/plot/table/download); instantiated per method x sex.
## clickable = TRUE (used only for the Overlap/consensus panel, which draws a Venn
## diagram whose regions correspond to real gene sets) wires the plot's click event
## through to the server, so clicking a region can filter the table below to just
## that region's genes - see the "_plot_click" observer in register_sex_technique_outputs().
mod_featureselection_technique_panel <- function(ns, prefix, title, plot_height = 300, clickable = FALSE) {
  box(
    width = NULL, title = title, status = "primary", solidHeader = FALSE,
    withSpinner(uiOutput(ns(paste0(prefix, "_summary"))), color = "#2563EB", type = 6),
    withSpinner(
      plotOutput(ns(paste0(prefix, "_plot")), height = plot_height,
                 click = if (clickable) ns(paste0(prefix, "_plot_click")) else NULL),
      color = "#2563EB", type = 6
    ),
    if (clickable) div(class = "empty-note", style = "margin-top: -4px; margin-bottom: 8px;", icon("hand-pointer"),
                        "Click a region of the diagram to filter the table below to just that region's genes."),
    div(class = "table-toolbar", downloadButton(ns(paste0(prefix, "_download")), "Genes (CSV)", class = "btn-sm")),
    DT::dataTableOutput(ns(paste0(prefix, "_table")))
  )
}

## One method's parameter box; always visible/live-editable, one independent instance per method.
mod_featureselection_params_box <- function(ns, prefix, method_label, defaults_desc, ...) {
  box(
    width = 12, title = sprintf("%s parameters", method_label), status = "primary", solidHeader = FALSE,
    p(class = "submodule-desc", defaults_desc),
    ...
  )
}

mod_featureselection_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        4,
        box(
          width = NULL, title = "Candidate genes & samples", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Female and male fit separately by default; \"Run All\" pools every sample."),
          # run buttons pinned at the top via .fs-run-bar (custom.css) so they stay reachable while scrolling settings
          div(
            class = "fs-run-bar",
            div(style = "display:flex; gap:8px; flex-wrap: wrap;",
                actionButton(ns("run_female_btn"), "Run Female", icon = icon("play"), class = "btn-primary btn-sm"),
                actionButton(ns("run_male_btn"), "Run Male", icon = icon("play"), class = "btn-primary btn-sm"),
                actionButton(ns("run_pooled_btn"), "Run All (pooled)", icon = icon("play"), class = "btn-primary btn-sm")
            ),
            div(style = "margin-top: 6px;", uiOutput(ns("speed_hint_ui")))
          ),
          tags$hr(),
          radioButtons(
            ns("data_source"), NULL,
            choiceNames = list(
              tagList(icon("diagram-project"), " Follow this project's pipeline (recommended)"),
              tagList(icon("table"), " Upload my own expression data"),
              tagList(icon("file-arrow-up"), " Upload my own DEG / candidate gene list")
            ),
            choiceValues = list("project", "expr", "deg"), selected = "project"
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'project'", ns("data_source")),
            uiOutput(ns("project_source_ui"))
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'expr'", ns("data_source")),
            p(strong("Expression matrix"), " - CSV or RDS. Genes in rows, samples in columns; for CSV, the first column is the gene ID."),
            fileInput(ns("expr_file"), "Expression matrix", accept = c(".csv", ".rds", ".Rds")),
            p(strong("Sample metadata"), " - CSV or RDS data frame, one row per sample, with a group and a sex column."),
            fileInput(ns("meta_file"), "Sample metadata", accept = c(".csv", ".rds", ".Rds")),
            uiOutput(ns("expr_column_mapping")),
            radioButtons(ns("gene_source"), "Candidate genes", inline = TRUE,
                         choices = c("Most variable genes" = "variable", "Paste my own list" = "custom", "A WGCNA module from this session" = "wgcna_module"),
                         selected = "variable"),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'variable'", ns("gene_source")),
              numericInput(ns("n_genes"), "Number of most variable genes (per sex)", value = 50, min = 10, max = 150, step = 10),
              div(class = "empty-note", style = "font-size: 12.5px;", icon("circle-info"),
                  "Keep this modest - SVM-RFE gets slow above 100 genes.")
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'custom'", ns("gene_source")),
              textAreaInput(ns("gene_list"), NULL, rows = 5, placeholder = "TNF\nIL6\nSTAT3\n...")
            ),
            # separate picker (wgcna_module_pick_expr) from the "deg" source's own, to avoid a duplicate DOM id
            conditionalPanel(
              condition = sprintf("input['%s'] == 'wgcna_module'", ns("gene_source")),
              uiOutput(ns("wgcna_module_pick_expr_ui"))
            )
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'deg'", ns("data_source")),
            p(class = "submodule-desc", "Expression values for these genes come from whatever dataset is currently loaded in the app (see the Dataset tab)."),
            radioButtons(ns("deg_source_mode"), NULL,
                         choices = c("Uploaded file" = "file", "A WGCNA module from this session" = "wgcna"),
                         selected = "file", inline = TRUE),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'file'", ns("deg_source_mode")),
              p(class = "submodule-desc", "One CSV per button, with a \"gene\" column. An optional \"direction\" or \"adj.P.Val\" column restricts to significant genes."),
              fileInput(ns("female_deg_file"), "Female DEG / candidate gene file", accept = c(".csv", ".tsv", ".txt")),
              fileInput(ns("male_deg_file"), "Male DEG / candidate gene file", accept = c(".csv", ".tsv", ".txt")),
              fileInput(ns("pooled_deg_file"), "Pooled (all-sample) DEG / candidate gene file - for \"Run All\"", accept = c(".csv", ".tsv", ".txt"))
            ),
            # reads results$wgcna$module_genes directly (from mod_wgcna.R Step 3); same module used for all sexes
            conditionalPanel(
              condition = sprintf("input['%s'] == 'wgcna'", ns("deg_source_mode")),
              uiOutput(ns("wgcna_module_pick_ui"))
            )
          ),
          tags$hr(),
          uiOutput(ns("group_controls_ui")),
          div(style = "margin-top:10px;", uiOutput(ns("saved_runs_ui")))
        )
      ),
      column(
        8,
        # hidden (not absent) via shinyjs so outputs stay bound in the DOM before the first Run click
        shinyjs::hidden(div(
        id = ns("fs_results_wrap"),
        tabsetPanel(
          id = ns("technique_tabs"), type = "tabs",
          # static UI (not req()-gated) so shinycssloaders spinners have something to attach to on first click
          tabPanel(
            "LASSO", br(),
            uiOutput(ns("lasso_params_ui")),
            div(style = "margin-bottom:10px;",
                actionButton(ns("lasso_show_btn"), "Show LASSO Results", icon = icon("eye"), class = "btn-outline-primary btn-sm")),
            conditionalPanel(
              condition = sprintf("input['%s'] > 0", ns("lasso_show_btn")),
              div(class = "fs-pair-row",
                fs_sex_panel(ns, "run_female_btn", mod_featureselection_technique_panel(ns, "female_lasso", "Female - LASSO")),
                fs_sex_panel(ns, "run_male_btn", mod_featureselection_technique_panel(ns, "male_lasso", "Male - LASSO"))
              ),
              fluidRow(column(12, fs_sex_panel(ns, "run_pooled_btn", mod_featureselection_technique_panel(ns, "pooled_lasso", "Pooled (all) - LASSO"))))
            )
          ),
          tabPanel(
            "Random Forest", br(),
            uiOutput(ns("rf_params_ui")),
            div(style = "margin-bottom:10px;",
                actionButton(ns("rf_show_btn"), "Show Random Forest Results", icon = icon("eye"), class = "btn-outline-primary btn-sm")),
            conditionalPanel(
              condition = sprintf("input['%s'] > 0", ns("rf_show_btn")),
              div(class = "fs-pair-row",
                fs_sex_panel(ns, "run_female_btn", mod_featureselection_technique_panel(ns, "female_rf", "Female - Random Forest")),
                fs_sex_panel(ns, "run_male_btn", mod_featureselection_technique_panel(ns, "male_rf", "Male - Random Forest"))
              ),
              fluidRow(column(12, fs_sex_panel(ns, "run_pooled_btn", mod_featureselection_technique_panel(ns, "pooled_rf", "Pooled (all) - Random Forest"))))
            )
          ),
          tabPanel(
            "SVM-RFE", br(),
            uiOutput(ns("svm_params_ui")),
            div(style = "margin-bottom:10px;",
                actionButton(ns("svm_show_btn"), "Show SVM-RFE Results", icon = icon("eye"), class = "btn-outline-primary btn-sm")),
            conditionalPanel(
              condition = sprintf("input['%s'] > 0", ns("svm_show_btn")),
              div(class = "fs-pair-row",
                fs_sex_panel(ns, "run_female_btn", mod_featureselection_technique_panel(ns, "female_svm", "Female - SVM-RFE")),
                fs_sex_panel(ns, "run_male_btn", mod_featureselection_technique_panel(ns, "male_svm", "Male - SVM-RFE"))
              ),
              fluidRow(column(12, fs_sex_panel(ns, "run_pooled_btn", mod_featureselection_technique_panel(ns, "pooled_svm", "Pooled (all) - SVM-RFE"))))
            )
          ),
          tabPanel(
            "Overlap", br(),
            uiOutput(ns("consensus_params_ui")),
            div(class = "fs-pair-row",
              fs_sex_panel(ns, "run_female_btn", mod_featureselection_technique_panel(ns, "female_consensus", "Female - Overlap", plot_height = 340, clickable = TRUE)),
              fs_sex_panel(ns, "run_male_btn", mod_featureselection_technique_panel(ns, "male_consensus", "Male - Overlap", plot_height = 340, clickable = TRUE))
            ),
            fluidRow(column(12, fs_sex_panel(ns, "run_pooled_btn", mod_featureselection_technique_panel(ns, "pooled_consensus", "Pooled (all) - Overlap", plot_height = 340, clickable = TRUE))))
          )
        ),
        box(
          width = NULL, title = "Result", status = "primary", solidHeader = FALSE,
          withSpinner(uiOutput(ns("summary_ui")), color = "#2563EB", type = 6)
        )
        ))
      )
    )
  )
}

mod_featureselection_server <- function(id, dataset, results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # own expression/metadata upload, separate from the Dataset tab, since sex mapping is mandatory here
    own_meta_raw <- reactive({
      req(input$meta_file)
      path <- input$meta_file$datapath
      if (grepl("\\.rds$", input$meta_file$name, ignore.case = TRUE)) {
        loaded <- safe_read_rds(path)
        validate(need(isTRUE(loaded$ok), loaded$error %||% "Could not read this .rds file."))
        d <- loaded$value
        validate(need(is.data.frame(d), "The uploaded metadata RDS file must contain a data frame."))
        as.data.frame(d)
      } else {
        as.data.frame(data.table::fread(path, showProgress = FALSE))
      }
    })

    output$expr_column_mapping <- renderUI({
      req(input$meta_file)
      cols <- colnames(own_meta_raw())
      tagList(
        selectInput(ns("map_id"), "Sample ID column", choices = cols, selected = cols[1], selectize = FALSE),
        selectInput(ns("map_group"), "Group / diagnosis column", choices = cols, selectize = FALSE),
        selectInput(ns("map_sex"), "Sex column", choices = cols, selectize = FALSE)
      )
    })

    # unified expr/meta source: "project"/"deg" read the app-wide dataset; "expr" reads the own upload above
    source_expr_meta <- reactive({
      if (identical(input$data_source, "expr")) {
        req(input$expr_file, input$meta_file, input$map_id, input$map_group, input$map_sex)
        expr <- if (grepl("\\.rds$", input$expr_file$name, ignore.case = TRUE)) {
          res <- tx_parse_expr_matrix_rds(input$expr_file$datapath)
          validate(need(res$ok, res$error))
          res$mat
        } else {
          m <- as.data.frame(data.table::fread(input$expr_file$datapath, showProgress = FALSE))
          rn <- as.character(m[[1]])
          m <- as.matrix(m[, -1, drop = FALSE])
          rownames(m) <- rn
          m
        }
        meta <- own_meta_raw()
        meta$sample <- as.character(meta[[input$map_id]])
        meta$group  <- as.character(meta[[input$map_group]])
        meta$sex    <- as.character(meta[[input$map_sex]])
        common <- intersect(colnames(expr), meta$sample)
        validate(need(length(common) >= 20, "Fewer than 20 sample IDs in the expression matrix match the metadata sample-ID column. Check the column mapping."))
        list(expr = expr[, common, drop = FALSE], meta = meta[match(common, meta$sample), , drop = FALSE])
      } else {
        req(dataset$expr, dataset$meta)
        # sex column is only required later, in fs_build_sex(), so "Run All" still works without one
        list(expr = dataset$expr, meta = dataset$meta)
      }
    })

    # maps the sex column's actual values to female/male (usually "F"/"M", but an upload might spell it out)
    sex_levels <- reactive({
      lv <- unique(stats::na.omit(as.character(source_expr_meta()$meta$sex)))
      validate(need(length(lv) >= 2, "The sex column needs at least two distinct values (e.g. F and M)."))
      f <- lv[grepl("^f", lv, ignore.case = TRUE)]
      m <- lv[grepl("^m", lv, ignore.case = TRUE)]
      lv_sorted <- sort(lv)
      list(female = if (length(f) > 0) f[1] else lv_sorted[1],
           male   = if (length(m) > 0) m[1] else lv_sorted[min(2, length(lv_sorted))])
    })

    output$group_controls_ui <- renderUI({
      meta <- tryCatch(source_expr_meta()$meta, error = function(e) NULL)
      if (is.null(meta) || !("group" %in% colnames(meta))) {
        return(div(class = "empty-note", icon("circle-info"), "Pick a candidate gene source above (with a group column) first."))
      }
      groups <- sort(unique(stats::na.omit(as.character(meta$group))))
      validate(need(length(groups) >= 2, "The group column needs at least two distinct values."))
      tagList(
        selectInput(ns("ref_group"), "Reference group", choices = groups, selected = groups[1], selectize = FALSE),
        selectInput(ns("comp_group"), "Comparison group", choices = groups, selected = groups[min(2, length(groups))], selectize = FALSE),
        radioButtons(ns("class_weight_mode"), "Class weighting (imbalanced groups)",
                     choices = c("Equal - this project's own methodology (default)" = "equal",
                                 "Balanced - auto inverse-frequency" = "balanced",
                                 "Manual ratio" = "manual"),
                     selected = "equal"),
        conditionalPanel(condition = sprintf("input['%s'] == 'manual'", ns("class_weight_mode")),
                          numericInput(ns("class_weight_ratio"), "Weight ratio (comparison : reference)", value = 1, min = 0.05, max = 20, step = 0.05)),
        div(class = "empty-note", style = "font-size: 12.5px; margin-top: -6px;", icon("circle-info"),
            "Applied to all three methods (LASSO, Random Forest, SVM-RFE) for both sexes.")
      )
    })

    # candidate gene panel per sex, one function per data source

    # "project" candidates: live Candidate Gene Identification output, else bundled MR-prioritised lists.
    # "pooled" has no candidate list of its own (candidate discovery is sex-stratified), so it uses the
    # union of the female/male panels instead. The bundled lists were only ever computed from the app's
    # own default merged cohort - only offer them when that exact dataset is still active
    # (dataset$is_bundled_reference), never for an uploaded/GEO/individual-preloaded dataset.
    project_candidate_genes <- function(sex_label) {
      if (identical(sex_label, "pooled")) {
        ## Candidate Gene Identification writes a genuinely non-sex-stratified
        ## panel to results$candidates$final (selection == "pooled") when the
        ## loaded dataset has no usable male/female metadata (see
        ## mod_candidates.R's sex_available()) - a live GEO fetch mapped to
        ## "(none)" for sex, most commonly. That IS the pooled candidate set
        ## for this dataset, so it's checked first, ahead of the
        ## union(female, male) derivation below (which only applies when the
        ## dataset *does* support sex-stratified discovery and both panels
        ## have been run).
        cand_final <- results$candidates$final
        if (!is.null(cand_final) && identical(cand_final$selection, "pooled") && length(cand_final$genes) >= 2) {
          return(list(genes = cand_final$genes, is_live = TRUE,
                      note = sprintf("%d genes from this session's live pooled candidate panel (Candidate Gene Identification - no sex-stratified metadata on this dataset).",
                                     length(cand_final$genes))))
        }
        live_f <- results$candidates$female$genes
        live_m <- results$candidates$male$genes
        live <- unique(c(live_f, live_m))
        if (length(live) >= 2) {
          return(list(genes = live, is_live = TRUE,
                      note = sprintf("%d genes - union of this session's live candidate panels (%d female, %d male).",
                                     length(live), length(unique(live_f)), length(unique(live_m)))))
        }
        if (isTRUE(dataset$is_bundled_reference)) {
          mhc <- if (isTRUE(input$mhc_exclude)) "_noMHC" else ""
          bf <- read_table_safe(sprintf("FS_input_female%s.csv", mhc))
          bm <- read_table_safe(sprintf("FS_input_male%s.csv", mhc))
          genes_f <- if (!is.null(bf) && "gene" %in% colnames(bf)) unique(as.character(bf$gene)) else character(0)
          genes_m <- if (!is.null(bm) && "gene" %in% colnames(bm)) unique(as.character(bm$gene)) else character(0)
          genes <- unique(c(genes_f, genes_m))
          if (length(genes) >= 2) {
            return(list(genes = genes, is_live = FALSE,
                        note = sprintf("%d genes - union of the bundled female + male candidate lists (%d female, %d male).",
                                       length(genes), length(genes_f), length(genes_m))))
          }
        }
        return(list(genes = character(0), is_live = FALSE,
                    note = "No pooled candidate genes available - run Candidate Gene Identification on the currently loaded dataset first."))
      }
      live <- results$candidates[[sex_label]]$genes
      if (!is.null(live) && length(live) >= 2) {
        return(list(genes = live, is_live = TRUE,
                    note = sprintf("%d genes from this session's live %s candidate panel.",
                                   length(live), sex_label)))
      }
      if (isTRUE(dataset$is_bundled_reference)) {
        fname <- sprintf("FS_input_%s%s.csv", sex_label, if (isTRUE(input$mhc_exclude)) "_noMHC" else "")
        bundled <- read_table_safe(fname)
        if (!is.null(bundled) && nrow(bundled) >= 2 && "gene" %in% colnames(bundled)) {
          return(list(genes = unique(as.character(bundled$gene)), is_live = FALSE,
                      note = sprintf("%d genes from the bundled %s candidate list (%s).",
                                     nrow(bundled), sex_label, fname)))
        }
      }
      list(genes = character(0), is_live = FALSE,
           note = sprintf("No live %s candidate genes yet - run Candidate Gene Identification on the currently loaded dataset first.", sex_label))
    }

    # fast path: precomputed LASSO/RF/SVM-RFE run (ml_features.rds / ml_features_noMHC.rds), served
    # instead of recomputing when the default pipeline is used uncustomized (see use_fast_path below).
    load_precomputed_fs <- function(sex_label, mhc_exclude) {
      fname <- if (isTRUE(mhc_exclude)) "ml_features_noMHC.rds" else "ml_features.rds"
      path <- file.path(PROCESSED_NEW_DIR, fname)
      if (!file.exists(path)) return(NULL)
      obj <- tryCatch(readRDS(path), error = function(e) NULL)
      sx <- obj[[sex_label]]
      if (is.null(sx)) return(NULL)
      grp <- obj$meta$group[match(sx$samples, obj$meta$sample)]
      list(
        cv = sx$cv, lasso_genes = sx$sets$LASSO, lasso_alpha = 1, lasso_lambda_choice = "lambda.min",
        gini = sx$gini, gini_thr = sx$gini_thr, rf_mtry = sx$tuned$rf_mtry, rf_ntree = 1000,
        rf_genes = sx$sets$RandomForest, rf_selection_rule = "above_mean", rf_top_n = NA,
        svm_rank = sx$svm_rank, svm_curve = sx$svm_curve, svm_cost = sx$tuned$svm_cost,
        svm_genes = sx$sets$SVM_RFE, svm_panel_mode = "auto", svm_manual_k = NA,
        sets = sx$sets, consensus = sx$consensus,
        n_input = length(sx$mr_genes), n_samples = length(sx$samples),
        n_ref = sum(grp == "HC", na.rm = TRUE), n_comp = sum(grp == "RA", na.rm = TRUE),
        fast_path = TRUE,
        candidate_note = sprintf("%d genes from the precomputed %s run (%s) - shown instantly. Change any parameter to run live instead.",
                                  length(sx$mr_genes), sex_label, fname)
      )
    }

    # "deg" candidates: uploaded gene list/DEG table per sex, or a WGCNA module's gene list (same module for all sexes).
    deg_candidate_genes <- function(sex_label) {
      if (identical(input$deg_source_mode, "wgcna")) {
        req(input$wgcna_module_pick)
        mg <- results$wgcna$module_genes
        validate(need(!is.null(mg) && input$wgcna_module_pick %in% names(mg),
                      "Run WGCNA Step 3 (Modules) first, then pick a module above."))
        genes <- unique(as.character(mg[[input$wgcna_module_pick]]))
        return(list(genes = genes,
                    note = sprintf("%d genes from WGCNA module \"%s\" (this session).",
                                    length(genes), input$wgcna_module_pick)))
      }
      fi <- switch(sex_label,
        female = input$female_deg_file, male = input$male_deg_file, pooled = input$pooled_deg_file)
      req(fi)
      d <- as.data.frame(data.table::fread(fi$datapath, showProgress = FALSE))
      colnames(d) <- tolower(colnames(d))
      validate(need("gene" %in% colnames(d), sprintf("The uploaded %s file needs a \"gene\" column.", sex_label)))
      if ("direction" %in% colnames(d)) {
        d <- d[!is.na(d$direction) & d$direction != "Not significant", , drop = FALSE]
      } else if ("adj.p.val" %in% colnames(d)) {
        d <- d[!is.na(d[["adj.p.val"]]) & d[["adj.p.val"]] < 0.05, , drop = FALSE]
      }
      genes <- unique(as.character(d$gene))
      list(genes = genes, note = sprintf("%d genes from your uploaded %s file (%s).", length(genes), sex_label, fi$name))
    }

    ## "expr": most-variable-in-this-sex-subset, a pasted list, or a WGCNA
    ## module's gene list, computed/intersected against the uploaded
    ## expression matrix.
    expr_candidate_genes <- function(sex_label, expr_sub) {
      if (identical(input$gene_source, "wgcna_module")) {
        req(input$wgcna_module_pick_expr)
        mg <- results$wgcna$module_genes
        validate(need(!is.null(mg) && input$wgcna_module_pick_expr %in% names(mg),
                      "Run WGCNA Step 3 (Modules) first, then pick a module above."))
        genes <- unique(as.character(mg[[input$wgcna_module_pick_expr]]))
        genes <- intersect(genes, rownames(expr_sub))
        list(genes = genes, note = sprintf("%d genes from WGCNA module \"%s\" present in the %s expression matrix.",
                                            length(genes), input$wgcna_module_pick_expr, sex_label))
      } else if (identical(input$gene_source, "custom")) {
        genes <- unique(trimws(unlist(strsplit(input$gene_list %||% "", "[,\n\t ]+"))))
        genes <- genes[nzchar(genes)]
        genes <- intersect(genes, rownames(expr_sub))
        list(genes = genes, note = sprintf("%d pasted genes present in the %s expression matrix.", length(genes), sex_label))
      } else {
        v <- apply(expr_sub, 1, stats::var)
        n <- min(input$n_genes %||% 50, length(v))
        genes <- names(sort(v, decreasing = TRUE))[seq_len(n)]
        list(genes = genes, note = sprintf("Top %d most variable genes in the %s subset.", n, sex_label))
      }
    }

    # shared module-picker builder for both WGCNA-module entry points; "grey" (unassigned) is excluded.
    wgcna_module_pick_ui_builder <- function(input_id, expr_note) {
      mg <- results$wgcna$module_genes
      if (is.null(mg) || length(mg) == 0) {
        return(div(class = "empty-note", icon("circle-info"),
                    "No WGCNA modules yet - run Step 3 (Modules) in the WGCNA Co-expression Network tab first."))
      }
      sizes <- vapply(mg, length, integer(1))
      choices <- setNames(names(mg), sprintf("%s (%d genes)", names(mg), sizes))
      choices <- choices[names(mg) != "grey"]
      tagList(
        selectInput(ns(input_id), "Module", choices = choices, selectize = FALSE),
        div(class = "empty-note", icon("circle-info"), expr_note)
      )
    }
    output$wgcna_module_pick_ui <- renderUI({
      wgcna_module_pick_ui_builder("wgcna_module_pick",
        "Uses expression values from the Dataset tab - use the same dataset WGCNA was run on.")
    })
    output$wgcna_module_pick_expr_ui <- renderUI({
      wgcna_module_pick_ui_builder("wgcna_module_pick_expr",
        "Genes are intersected with your uploaded expression matrix above.")
    })

    output$project_source_ui <- renderUI({
      f_live <- results$candidates$female$genes
      m_live <- results$candidates$male$genes
      has_sex_live <- length(f_live) >= 2 && length(m_live) >= 2
      ## Pooled (no sex-stratified metadata) live panel from Candidate Gene
      ## Identification - see project_candidate_genes("pooled") above for why
      ## this is checked separately from the female/male panels.
      cand_final <- results$candidates$final
      has_pooled_live <- is.null(results$candidates$female) && is.null(results$candidates$male) &&
        !is.null(cand_final) && identical(cand_final$selection, "pooled") && length(cand_final$genes) >= 2
      has_live <- has_sex_live || has_pooled_live
      tagList(
        if (has_sex_live) {
          div(class = "empty-note", icon("check"),
              sprintf("Using live candidates: %d female / %d male genes.", length(f_live), length(m_live)))
        } else if (has_pooled_live) {
          div(class = "empty-note", icon("check"),
              sprintf("Using %d live pooled candidates (this dataset has no sex-stratified metadata).", length(cand_final$genes)))
        } else {
          div(class = "empty-note", icon("circle-info"),
              "No live candidates yet - run Candidate Gene Identification first, or use the bundled example lists below.")
        },
        checkboxInput(ns("mhc_exclude"), "Exclude MHC region (chr6:25-34Mb) - sensitivity panel", value = FALSE),
        if (has_live) {
          div(class = "empty-note", style = "font-size: 12.5px;", icon("triangle-exclamation"),
              "MHC exclusion only applies to the bundled lists, not live candidates.")
        } else NULL
      )
    })

    val_eq <- function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-8))
    # true if any control differs from FS_DEFAULT_PARAMS; decides fast-path eligibility in use_fast_path below
    fs_any_customized <- function() {
      p <- fs_advanced_params()
      d <- FS_DEFAULT_PARAMS
      !(
        identical(p$class_weight_mode, d$class_weight_mode) && val_eq(p$class_weight_ratio, d$class_weight_ratio) &&
        val_eq(p$lasso_cv_folds, d$lasso_cv_folds) && val_eq(p$lasso_alpha, d$lasso_alpha) &&
        identical(p$lasso_lambda_choice, d$lasso_lambda_choice) &&
        val_eq(p$lasso_nlambda, d$lasso_nlambda) && identical(p$lasso_type_measure, d$lasso_type_measure) &&
        val_eq(p$rf_cv_folds, d$rf_cv_folds) && val_eq(p$rf_ntree, d$rf_ntree) &&
        identical(p$rf_mtry_mode, d$rf_mtry_mode) && identical(p$rf_selection_rule, d$rf_selection_rule) &&
        val_eq(p$rf_nodesize, d$rf_nodesize) && is.null(p$rf_maxnodes) == is.null(d$rf_maxnodes) &&
        val_eq(p$svm_cv_folds, d$svm_cv_folds) && identical(p$svm_cost_mode, d$svm_cost_mode) &&
        identical(p$svm_panel_mode, d$svm_panel_mode) && val_eq(p$svm_tolerance, d$svm_tolerance) &&
        isTRUE(all.equal(sort(as.numeric(p$svm_cost_grid)), sort(as.numeric(d$svm_cost_grid)), tolerance = 1e-8)) &&
        setequal(p$consensus_methods, d$consensus_methods)
      )
    }

    # tells the user whether the next Run will be instant or live; mirrors use_fast_path below exactly
    output$speed_hint_ui <- renderUI({
      note <- function(icon_name, txt) div(class = "empty-note", style = "font-size: 12.5px;", icon(icon_name), txt)

      if (!identical(input$data_source, "project")) {
        return(note("clock", "Uploaded data always runs live - can take a few seconds to a couple of minutes."))
      }
      f_live <- results$candidates$female$genes
      m_live <- results$candidates$male$genes
      if (length(f_live) >= 2 && length(m_live) >= 2) {
        big <- max(length(f_live), length(m_live)) > FS_MAX_CANDIDATE_GENES
        return(note("clock", sprintf(
          "Using live candidates (%d female / %d male genes) - always runs live.%s",
          length(f_live), length(m_live),
          if (big) sprintf(" Reduced to the top %d most variable genes per sex first.", FS_MAX_CANDIDATE_GENES) else ""
        )))
      }
      if (!isTRUE(dataset$is_bundled_reference)) {
        return(note("clock", "A non-default dataset is loaded, so this runs live."))
      }
      if (fs_any_customized()) {
        return(note("clock", "Parameters are customized, so this run will be live - can take a little while."))
      }
      if (!(identical(input$ref_group, "HC") && identical(input$comp_group, "RA"))) {
        return(note("clock", "Only the standard HC vs RA contrast is precomputed - this pick runs live."))
      }
      NULL
    })

    # reads current values of the LASSO/RF/SVM-RFE controls, falling back to FS_DEFAULT_PARAMS
    fs_advanced_params <- function() {
      cost_grid <- suppressWarnings(as.numeric(trimws(strsplit(input$svm_cost_grid %||% "", ",")[[1]])))
      cost_grid <- cost_grid[!is.na(cost_grid) & cost_grid > 0]
      list(
        class_weight_mode = input$class_weight_mode %||% FS_DEFAULT_PARAMS$class_weight_mode,
        class_weight_ratio = input$class_weight_ratio %||% FS_DEFAULT_PARAMS$class_weight_ratio,
        lasso_cv_folds = input$lasso_cv_folds %||% FS_DEFAULT_PARAMS$lasso_cv_folds,
        lasso_alpha = input$lasso_alpha %||% FS_DEFAULT_PARAMS$lasso_alpha,
        lasso_lambda_choice = input$lasso_lambda_choice %||% FS_DEFAULT_PARAMS$lasso_lambda_choice,
        lasso_nlambda = input$lasso_nlambda %||% FS_DEFAULT_PARAMS$lasso_nlambda,
        lasso_type_measure = input$lasso_type_measure %||% FS_DEFAULT_PARAMS$lasso_type_measure,
        rf_cv_folds = input$rf_cv_folds %||% FS_DEFAULT_PARAMS$rf_cv_folds,
        rf_ntree = input$rf_ntree %||% FS_DEFAULT_PARAMS$rf_ntree,
        rf_mtry_mode = input$rf_mtry_mode %||% FS_DEFAULT_PARAMS$rf_mtry_mode,
        rf_mtry_manual = input$rf_mtry_manual,
        rf_selection_rule = input$rf_selection_rule %||% FS_DEFAULT_PARAMS$rf_selection_rule,
        rf_top_n = input$rf_top_n %||% FS_DEFAULT_PARAMS$rf_top_n,
        rf_nodesize = input$rf_nodesize %||% FS_DEFAULT_PARAMS$rf_nodesize,
        rf_maxnodes = if (isTRUE(input$rf_maxnodes_unlimited)) NULL else (input$rf_maxnodes %||% FS_DEFAULT_PARAMS$rf_maxnodes),
        svm_cv_folds = input$svm_cv_folds %||% FS_DEFAULT_PARAMS$svm_cv_folds,
        svm_cost_mode = input$svm_cost_mode %||% FS_DEFAULT_PARAMS$svm_cost_mode,
        svm_cost_manual = input$svm_cost_manual %||% FS_DEFAULT_PARAMS$svm_cost_manual,
        svm_cost_grid = if (length(cost_grid) > 0) cost_grid else FS_SVM_COST_GRID,
        svm_panel_mode = input$svm_panel_mode %||% FS_DEFAULT_PARAMS$svm_panel_mode,
        svm_manual_k = input$svm_manual_k %||% FS_DEFAULT_PARAMS$svm_manual_k,
        svm_tolerance = input$svm_tolerance %||% FS_DEFAULT_PARAMS$svm_tolerance,
        consensus_methods = input$consensus_methods %||% FS_DEFAULT_PARAMS$consensus_methods
      )
    }

    # fits LASSO/RF/SVM-RFE for one sex, given the current data source, parameters, and group contrast
    fs_build_sex <- function(sex_label, sex_value) {
      req(input$ref_group, input$comp_group)
      validate(need(input$ref_group != input$comp_group, "Reference and comparison group must be different."))

      sem <- source_expr_meta()
      # only Female/Male need a usable sex column; "Run All" (sex_value NULL) pools regardless
      if (!is.null(sex_value)) {
        validate(need("sex" %in% colnames(sem$meta) && length(unique(stats::na.omit(sem$meta$sex))) >= 2,
                      "No usable sex column in this dataset - load one with sex data, or use \"Run All (pooled)\" instead."))
      }
      adv_params <- fs_advanced_params()
      any_customized <- fs_any_customized()

      # candidate lookup is cheap; resolve it first to decide fast-path eligibility before the expensive fit
      cand_project <- if (identical(input$data_source, "project")) project_candidate_genes(sex_label) else NULL

      # fast path only when: default pipeline, unmodified bundled candidates, no customized params,
      # standard HC-vs-RA contrast, and the project's own example dataset still loaded
      use_fast_path <- identical(input$data_source, "project") &&
        !isTRUE(cand_project$is_live) && !any_customized &&
        identical(input$ref_group, "HC") && identical(input$comp_group, "RA") &&
        isTRUE(dataset$is_bundled_reference)

      if (use_fast_path) {
        fit <- load_precomputed_fs(sex_label, isTRUE(input$mhc_exclude))
        if (!is.null(fit)) {
          fit$ref_group <- input$ref_group; fit$comp_group <- input$comp_group
          fit$mhc_mode <- if (isTRUE(input$mhc_exclude)) "exclude" else "include"
          return(fit)
        }
      }

      meta <- sem$meta
      sex_ok <- if (is.null(sex_value)) rep(TRUE, nrow(meta)) else (!is.na(meta$sex) & as.character(meta$sex) == sex_value)
      meta <- meta[sex_ok & !is.na(meta$group) & as.character(meta$group) %in% c(input$ref_group, input$comp_group), , drop = FALSE]
      common <- intersect(colnames(sem$expr), meta$sample)
      validate(need(length(common) >= 10, sprintf("Fewer than 10 %s samples match this contrast; feature selection needs more samples to be meaningful.", sex_label)))
      meta <- meta[match(common, meta$sample), , drop = FALSE]
      expr_sub <- sem$expr[, common, drop = FALSE]

      cand <- switch(input$data_source,
        project = cand_project,
        deg = deg_candidate_genes(sex_label),
        expr = expr_candidate_genes(sex_label, expr_sub)
      )
      genes <- intersect(cand$genes, rownames(expr_sub))
      validate(need(length(genes) >= 3, sprintf("Fewer than 3 %s candidate genes are present in the currently loaded expression matrix.", sex_label)))

      # cap candidate count (see FS_MAX_CANDIDATE_GENES) by keeping the most variable genes in this sex's samples
      n_before_cap <- length(genes)
      if (n_before_cap > FS_MAX_CANDIDATE_GENES) {
        v <- apply(expr_sub[genes, , drop = FALSE], 1, stats::var)
        genes <- names(sort(v, decreasing = TRUE))[seq_len(FS_MAX_CANDIDATE_GENES)]
        cand$note <- sprintf("%s Reduced from %d to the %d most variable genes before fitting.",
                              cand$note, n_before_cap, FS_MAX_CANDIDATE_GENES)
      }

      X <- t(expr_sub[genes, , drop = FALSE])
      y <- factor(as.character(meta$group), levels = c(input$ref_group, input$comp_group))
      grp_counts <- table(y)
      validate(need(all(grp_counts >= FS_MIN_GROUP_SAMPLES), sprintf(
        "Each group needs at least %d %s samples for feature selection (this contrast has %s). This is about patients, not candidate genes - a larger candidate gene panel doesn't add patients to a group.",
        FS_MIN_GROUP_SAMPLES, sex_label, paste(sprintf("%s = %d", names(grp_counts), as.integer(grp_counts)), collapse = ", ")
      )))

      fit <- fs_fit_sex(X, y, params = adv_params)
      fit$candidate_note <- cand$note
      fit$n_ref <- sum(y == input$ref_group)
      fit$n_comp <- sum(y == input$comp_group)
      fit$min_group_n <- min(grp_counts)
      fit$ref_group <- input$ref_group; fit$comp_group <- input$comp_group
      fit$mhc_mode <- if (identical(input$data_source, "project")) (if (isTRUE(input$mhc_exclude)) "exclude" else "include") else "n/a"
      fit
    }

    # ignoreInit avoids an early req() halt before ref_group/comp_group have round-tripped to the server
    fs_result_female <- eventReactive(input$run_female_btn, {
      fs_build_sex("female", sex_levels()$female)
    }, ignoreInit = TRUE)
    fs_result_male <- eventReactive(input$run_male_btn, {
      fs_build_sex("male", sex_levels()$male)
    }, ignoreInit = TRUE)
    # sex_value = NULL skips the sex filter entirely, pooling every sample regardless of sex
    fs_result_pooled <- eventReactive(input$run_pooled_btn, {
      fs_build_sex("pooled", NULL)
    }, ignoreInit = TRUE)

    # Each fs_result_*() is an eventReactive, so it holds the fit from whichever dataset was loaded
    # when its Run button was last clicked - switching datasets would otherwise leave that fit on
    # screen labelled as the current one. Marked stale per sex here, cleared by that sex's own Run
    # button below, and read by fs_result_or_null()/fs_result_error_msg() so every downstream
    # summary/plot/table/Venn falls back to the "not run yet" empty state.
    fs_stale <- reactiveValues(female = FALSE, male = FALSE, pooled = FALSE)
    observeEvent(dataset$source, {
      fs_stale$female <- TRUE; fs_stale$male <- TRUE; fs_stale$pooled <- TRUE
    }, ignoreInit = TRUE)
    fs_is_stale <- function(sex_label) isTRUE(fs_stale[[sex_label]])

    # Reads one sex's fs_result_*() and returns its real validate()/need() failure message, or
    # NULL if it hasn't failed. Distinguishes a genuine failure from "hasn't been run yet": both
    # are eventReactive halts of the same shiny "validation" class, but eventReactive's own
    # pre-first-click halt (ignoreInit's req(FALSE)) always carries an EMPTY message, while a
    # validate(need(...)) failure inside fs_build_sex() always carries the real one - so
    # `nzchar()` on the caught message tells them apart. Every other read of fs_result_*()
    # elsewhere in this file (`tryCatch(..., error = function(e) NULL)`) collapsed both cases to
    # NULL, which made a real failure (e.g. too few female samples for this contrast) look
    # exactly like the button never having been clicked - the bug behind "Female result is not
    # showing" with no indication why.
    fs_result_error_msg <- function(sex_label) {
      if (fs_is_stale(sex_label)) return(NULL)
      fr <- switch(sex_label, female = fs_result_female, male = fs_result_male, pooled = fs_result_pooled)
      tryCatch({ fr(); NULL }, error = function(e) {
        msg <- conditionMessage(e)
        if (nzchar(msg)) msg else NULL
      })
    }

    # Single stale-aware read of one sex's fit: NULL both when it has never been run and when the
    # last run belongs to a dataset that is no longer loaded.
    fs_result_or_null <- function(sex_label) {
      if (fs_is_stale(sex_label)) return(NULL)
      fr <- switch(sex_label, female = fs_result_female, male = fs_result_male, pooled = fs_result_pooled)
      tryCatch(fr(), error = function(e) NULL)
    }

    # reveals the results area on the first Run click (stays visible after); raw id since shinyjs auto-namespaces
    fs_has_run <- reactiveVal(FALSE)
    observeEvent(input$run_female_btn, {
      fs_has_run(TRUE)
      fs_stale$female <- FALSE
      shinyjs::show(id = "fs_results_wrap")
      err <- fs_result_error_msg("female")
      if (!is.null(err)) showNotification(paste("Female feature selection failed:", err), type = "error", duration = 10)
    }, ignoreInit = TRUE)
    observeEvent(input$run_male_btn, {
      fs_has_run(TRUE)
      fs_stale$male <- FALSE
      shinyjs::show(id = "fs_results_wrap")
      err <- fs_result_error_msg("male")
      if (!is.null(err)) showNotification(paste("Male feature selection failed:", err), type = "error", duration = 10)
    }, ignoreInit = TRUE)
    observeEvent(input$run_pooled_btn, {
      fs_has_run(TRUE)
      fs_stale$pooled <- FALSE
      shinyjs::show(id = "fs_results_wrap")
      err <- fs_result_error_msg("pooled")
      if (!is.null(err)) showNotification(paste("Pooled feature selection failed:", err), type = "error", duration = 10)
    }, ignoreInit = TRUE)

    output$lasso_params_ui <- renderUI({
      req(fs_has_run())
      mod_featureselection_params_box(
        ns, "lasso", "LASSO",
        "Defaults: alpha = 1, lambda.min, 10-fold CV. Change and re-run to use different settings.",
        fluidRow(
          column(4, numericInput(ns("lasso_cv_folds"), "Cross-validation folds", value = 10, min = 3, max = 10, step = 1)),
          column(4, sliderInput(ns("lasso_alpha"), "Alpha (1 = LASSO, 0 = ridge)", min = 0, max = 1, value = 1, step = 0.05)),
          column(4, radioButtons(ns("lasso_lambda_choice"), "Lambda", choices = c("lambda.min (default)" = "lambda.min", "lambda.1se (sparser)" = "lambda.1se"), selected = "lambda.min"))
        ),
        h5("Advanced", style = "margin-top: 6px;"),
        fluidRow(
          column(4, numericInput(ns("lasso_nlambda"), "Number of lambda values searched (nlambda)", value = FS_DEFAULT_PARAMS$lasso_nlambda, min = 20, max = 300, step = 10)),
          column(4, radioButtons(ns("lasso_type_measure"), "CV metric to optimize",
                                  choices = c("Deviance (default)" = "deviance", "AUC" = "auc", "Misclassification error" = "class"),
                                  selected = "deviance"))
        )
      )
    })

    output$rf_params_ui <- renderUI({
      req(fs_has_run())
      mod_featureselection_params_box(
        ns, "rf", "Random Forest",
        "Defaults: ntree = 1000, mtry tuned by CV, genes kept above mean Gini. Change and re-run to use different settings.",
        fluidRow(
          column(3, numericInput(ns("rf_cv_folds"), "Cross-validation folds", value = 10, min = 3, max = 10, step = 1)),
          column(3, numericInput(ns("rf_ntree"), "Number of trees", value = 1000, min = 100, max = 5000, step = 100)),
          column(3,
            radioButtons(ns("rf_mtry_mode"), "mtry (per split)", choices = c("Auto-tune (default)" = "auto", "Manual" = "manual"), selected = "auto"),
            conditionalPanel(condition = sprintf("input['%s'] == 'manual'", ns("rf_mtry_mode")),
                              numericInput(ns("rf_mtry_manual"), "mtry value", value = 5, min = 1, max = 500, step = 1))
          ),
          column(3,
            radioButtons(ns("rf_selection_rule"), "Selection rule", choices = c("Above mean Gini (default)" = "above_mean", "Top N" = "top_n"), selected = "above_mean"),
            conditionalPanel(condition = sprintf("input['%s'] == 'top_n'", ns("rf_selection_rule")),
                              numericInput(ns("rf_top_n"), "N genes to keep", value = 10, min = 1, max = 500, step = 1))
          )
        ),
        h5("Advanced", style = "margin-top: 6px;"),
        fluidRow(
          column(4, numericInput(ns("rf_nodesize"), "Minimum terminal node size (nodesize)", value = FS_DEFAULT_PARAMS$rf_nodesize, min = 1, max = 50, step = 1)),
          column(4,
            checkboxInput(ns("rf_maxnodes_unlimited"), "Unlimited tree depth (maxnodes, default)", value = TRUE),
            conditionalPanel(condition = sprintf("!input['%s']", ns("rf_maxnodes_unlimited")),
                              numericInput(ns("rf_maxnodes"), "Max terminal nodes per tree", value = 20, min = 2, max = 2000, step = 1))
          )
        )
      )
    })

    output$svm_params_ui <- renderUI({
      req(fs_has_run())
      mod_featureselection_params_box(
        ns, "svm", "SVM-RFE",
        "Defaults: linear kernel, cost tuned by CV, panel size = top-k minimising CV error. Change and re-run to use different settings.",
        fluidRow(
          column(3, numericInput(ns("svm_cv_folds"), "Cross-validation folds", value = 10, min = 3, max = 10, step = 1)),
          column(4,
            radioButtons(ns("svm_cost_mode"), "Cost (C)", choices = c("Auto-tune via CV grid (default)" = "auto", "Manual" = "manual"), selected = "auto"),
            conditionalPanel(condition = sprintf("input['%s'] == 'auto'", ns("svm_cost_mode")),
                              textInput(ns("svm_cost_grid"), "Cost grid (comma-separated)", value = paste(FS_SVM_COST_GRID, collapse = ", "))),
            conditionalPanel(condition = sprintf("input['%s'] == 'manual'", ns("svm_cost_mode")),
                              numericInput(ns("svm_cost_manual"), "Cost value", value = 1, min = 0.001, step = 0.1))
          ),
          column(5,
            radioButtons(ns("svm_panel_mode"), "Panel size (top-ranked genes)", choices = c("Auto: minimise CV error (default)" = "auto", "Manual" = "manual"), selected = "auto"),
            conditionalPanel(condition = sprintf("input['%s'] == 'manual'", ns("svm_panel_mode")),
                              numericInput(ns("svm_manual_k"), "k (top-ranked genes to keep)", value = 10, min = 1, max = 500, step = 1))
          )
        ),
        h5("Advanced", style = "margin-top: 6px;"),
        fluidRow(
          column(4, numericInput(ns("svm_tolerance"), "Convergence tolerance", value = FS_DEFAULT_PARAMS$svm_tolerance, min = 0.00001, max = 0.1, step = 0.0001))
        )
      )
    })

    # methods to intersect for consensus; all three by default, untick to replicate e.g. Chen et al. 2021/2022's LASSO∩SVM-RFE panel
    output$consensus_params_ui <- renderUI({
      req(fs_has_run())
      mod_featureselection_params_box(
        ns, "consensus", "Overlap",
        "All three methods are intersected by default. Pick LASSO only, any two, or all three - the overlap below updates instantly, no re-run needed.",
        checkboxGroupInput(ns("consensus_methods"), "Methods to intersect",
                            choices = c("LASSO" = "LASSO", "Random Forest" = "RandomForest", "SVM-RFE" = "SVM_RFE"),
                            selected = c("LASSO", "RandomForest", "SVM_RFE"), inline = TRUE)
      )
    })

    # each sex saves into results$featureselection independently via modifyList, without clobbering the other
    observeEvent(fs_result_female(), {
      r <- fs_result_female()
      results$featureselection <- utils::modifyList(
        results$featureselection %||% list(),
        list(
          data_source = input$data_source,
          female = list(n_input = r$n_input, n_samples = r$n_samples,
                        n_lasso = length(r$lasso_genes), n_rf = length(r$rf_genes),
                        n_svm = length(r$svm_genes), n_consensus = length(r$consensus),
                        consensus_genes = r$consensus)
        )
      )

      runs <- results$featureselection_runs %||% list()
      run_id <- paste0("run", length(runs) + 1L)
      runs[[run_id]] <- list(
        sex = "female", contrast = sprintf("%s vs %s", r$comp_group, r$ref_group), data_source = input$data_source,
        female_consensus = r$consensus, male_consensus = NULL
      )
      if (length(runs) > 8) runs <- utils::tail(runs, 8)
      results$featureselection_runs <- runs

      showNotification(
        sprintf("Female feature selection saved: %d consensus genes%s.",
                length(r$consensus), if (isTRUE(r$fast_path)) " (instant, precomputed)" else ""),
        type = "message", duration = 6
      )
    })

    observeEvent(fs_result_male(), {
      r <- fs_result_male()
      results$featureselection <- utils::modifyList(
        results$featureselection %||% list(),
        list(
          data_source = input$data_source,
          male = list(n_input = r$n_input, n_samples = r$n_samples,
                      n_lasso = length(r$lasso_genes), n_rf = length(r$rf_genes),
                      n_svm = length(r$svm_genes), n_consensus = length(r$consensus),
                      consensus_genes = r$consensus)
        )
      )

      runs <- results$featureselection_runs %||% list()
      run_id <- paste0("run", length(runs) + 1L)
      runs[[run_id]] <- list(
        sex = "male", contrast = sprintf("%s vs %s", r$comp_group, r$ref_group), data_source = input$data_source,
        female_consensus = NULL, male_consensus = r$consensus
      )
      if (length(runs) > 8) runs <- utils::tail(runs, 8)
      results$featureselection_runs <- runs

      showNotification(
        sprintf("Male feature selection saved: %d consensus genes%s.",
                length(r$consensus), if (isTRUE(r$fast_path)) " (instant, precomputed)" else ""),
        type = "message", duration = 6
      )
    })

    observeEvent(fs_result_pooled(), {
      r <- fs_result_pooled()
      results$featureselection <- utils::modifyList(
        results$featureselection %||% list(),
        list(
          data_source = input$data_source,
          pooled = list(n_input = r$n_input, n_samples = r$n_samples,
                        n_lasso = length(r$lasso_genes), n_rf = length(r$rf_genes),
                        n_svm = length(r$svm_genes), n_consensus = length(r$consensus),
                        consensus_genes = r$consensus)
        )
      )

      runs <- results$featureselection_runs %||% list()
      run_id <- paste0("run", length(runs) + 1L)
      runs[[run_id]] <- list(
        sex = "pooled", contrast = sprintf("%s vs %s", r$comp_group, r$ref_group), data_source = input$data_source,
        female_consensus = NULL, male_consensus = NULL, pooled_consensus = r$consensus
      )
      if (length(runs) > 8) runs <- utils::tail(runs, 8)
      results$featureselection_runs <- runs

      showNotification(
        sprintf("Pooled (all) feature selection saved: %d consensus genes%s.",
                length(r$consensus), if (isTRUE(r$fast_path)) " (instant, precomputed)" else ""),
        type = "message", duration = 6
      )
    })

    output$saved_runs_ui <- renderUI({
      res_f <- fs_result_or_null("female")
      res_m <- fs_result_or_null("male")
      res_p <- fs_result_or_null("pooled")
      status_row <- function(sex, sex_label, r) {
        if (!is.null(r)) {
          ## Ran (cleared FS_MIN_GROUP_SAMPLES), but still below FS_RELIABLE_GROUP_SAMPLES -
          ## flagged, not blocked, so a genuinely small group (e.g. 4-5 patients) still
          ## produces a result, just clearly marked as exploratory.
          low_n <- !isTRUE(r$fast_path) && !is.null(r$min_group_n) && r$min_group_n < FS_RELIABLE_GROUP_SAMPLES
          return(tags$li(icon("check", style = "color: #1a9c5f;"), strong(sprintf(" %s feature selection completed: ", sex)),
                          sprintf("%d consensus genes%s", length(r$consensus), if (isTRUE(r$fast_path)) " (instant)" else " (live)"),
                          if (low_n) span(style = "color: #b8860b;", sprintf(" - smallest group has only %d samples, treat as exploratory", r$min_group_n)) else NULL))
        }
        ## Real validate() failure (e.g. too few female samples) vs. genuinely never clicked -
        ## see fs_result_error_msg() above; both used to render as this same "not run yet" row.
        err <- fs_result_error_msg(sex_label)
        if (!is.null(err)) {
          tags$li(icon("triangle-exclamation", style = "color: #c0392b;"), strong(sprintf(" %s feature selection failed: ", sex)), err)
        } else {
          tags$li(icon("circle-minus", style = "color: #8A929C;"), sprintf(" %s feature selection - not run yet", sex))
        }
      }
      tagList(
        p(class = "submodule-desc", style = "margin-bottom: 4px;", "Status:"),
        tags$ul(style = "padding-left: 18px; margin-bottom: 0; list-style: none;",
                status_row("Female", "female", res_f),
                status_row("Male", "male", res_m),
                status_row("Pooled (all)", "pooled", res_p))
      )
    })

    # renders each sex's summary line independently, showing "not run yet" for the others
    output$summary_ui <- renderUI({
      res_f <- fs_result_or_null("female")
      res_m <- fs_result_or_null("male")
      res_p <- fs_result_or_null("pooled")
      mhc_mode <- if (identical(input$data_source, "project")) (if (isTRUE(input$mhc_exclude)) "exclude" else "include") else "n/a"
      sex_line <- function(sex, sex_label, r) {
        if (!is.null(r)) {
          low_n <- !isTRUE(r$fast_path) && !is.null(r$min_group_n) && r$min_group_n < FS_RELIABLE_GROUP_SAMPLES
          return(p(strong(sprintf("%s: ", sex)),
            if (isTRUE(r$fast_path)) icon("bolt") else NULL,
            sprintf(" %d candidate genes, %d samples (%d vs %d) → %d LASSO / %d random forest / %d SVM-RFE → %d consensus%s.",
                    r$n_input, r$n_samples, r$n_comp, r$n_ref,
                    length(r$lasso_genes), length(r$rf_genes), length(r$svm_genes), length(r$consensus),
                    if (isTRUE(r$fast_path)) " (instant, precomputed)" else " (live)"),
            if (low_n) span(style = "color: #b8860b;",
              sprintf(" Smallest group has only %d samples (below the recommended %d) - treat this fit as exploratory, not a reliable panel.",
                      r$min_group_n, FS_RELIABLE_GROUP_SAMPLES)) else NULL))
        }
        ## Same never-run vs. actually-failed distinction as status_row() above.
        err <- fs_result_error_msg(sex_label)
        if (!is.null(err)) return(p(strong(sprintf("%s: ", sex)), span(style = "color: #c0392b;", sprintf("failed - %s", err))))
        p(strong(sprintf("%s: ", sex)), "not run yet.")
      }
      tagList(
        p(strong("Contrast: "), sprintf("%s vs %s", input$comp_group %||% "?", input$ref_group %||% "?"),
          if (!identical(mhc_mode, "n/a")) sprintf(" (MHC region %s)", if (identical(mhc_mode, "exclude")) "excluded" else "included") else NULL, "."),
        sex_line("Female", "female", res_f), sex_line("Male", "male", res_m), sex_line("Pooled (all)", "pooled", res_p)
      )
    })

    # returns NULL (not a req()-halt) when a sex hasn't been run yet, so outputs can show "not run yet"
    res_sex <- function(sex_label) reactive(fs_result_or_null(sex_label))

    # binds all summary/plot/table/download outputs for one sex's LASSO/RF/SVM-RFE/consensus boxes
    register_sex_technique_outputs <- function(sex_label, res) {
      sex_color <- switch(sex_label, female = "#1a7a3c", male = "#7a4a26", pooled = "#2563EB")
      sex_title <- tools::toTitleCase(sex_label)
      # Shows the real validate()/need() failure (e.g. too few female samples for this contrast)
      # instead of the generic "not run yet" note when this sex's Run button WAS clicked but
      # fs_build_sex() failed - see fs_result_error_msg() above for how the two are told apart.
      # plot/table outputs below use validate(need(...)) rather than req(): a req() halt is silent,
      # so the previous dataset's already-rendered content would stay in the DOM.
      not_yet_msg <- function() {
        err <- fs_result_error_msg(sex_label)
        if (!is.null(err)) return(sprintf("%s feature selection failed: %s", sex_title, err))
        if (fs_is_stale(sex_label)) {
          return("Not run for the currently loaded dataset yet. Click Run Female, Run Male, or Run All (pooled) on the left.")
        }
        "No result yet. Click Run Female, Run Male, or Run All (pooled) on the left."
      }
      not_yet_note <- function() {
        failed <- !is.null(fs_result_error_msg(sex_label))
        div(class = "empty-note", icon(if (failed) "triangle-exclamation" else "circle-info"), not_yet_msg())
      }

      # LASSO
      output[[paste0(sex_label, "_lasso_summary")]] <- renderUI({
        r <- res()
        if (is.null(r)) return(not_yet_note())
        lambda_used <- if (identical(r$lasso_lambda_choice, "lambda.1se")) r$cv$lambda.1se else r$cv$lambda.min
        p(strong(length(r$lasso_genes)),
          sprintf(" of %d candidate genes selected (alpha = %.2f, %s = %.4f).", r$n_input, r$lasso_alpha, r$lasso_lambda_choice, lambda_used))
      })
      # r <- res() first (not inline plot(res()$cv)) so req()'s halt doesn't hit plot()'s S3 dispatch
      output[[paste0(sex_label, "_lasso_plot")]] <- renderPlot({
        r <- res()
        validate(need(!is.null(r), not_yet_msg()))
        plot(r$cv)
      })
      output[[paste0(sex_label, "_lasso_table")]] <- DT::renderDataTable({
        r <- res()
        validate(need(!is.null(r), not_yet_msg()))
        DT::datatable(data.frame(gene = r$lasso_genes, stringsAsFactors = FALSE), rownames = FALSE,
                      options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
      })
      output[[paste0(sex_label, "_lasso_download")]] <- downloadHandler(
        filename = function() sprintf("%s_lasso_genes.csv", sex_label),
        content = function(file) write.csv(data.frame(gene = tx_csv_safe(res()$lasso_genes)), file, row.names = FALSE)
      )

      # Random Forest
      output[[paste0(sex_label, "_rf_summary")]] <- renderUI({
        r <- res()
        if (is.null(r)) return(not_yet_note())
        rule_txt <- if (identical(r$rf_selection_rule, "top_n")) {
          sprintf("top %d by importance", length(r$rf_genes))
        } else {
          sprintf("Gini > mean = %.3f", r$gini_thr)
        }
        p(strong(length(r$rf_genes)), sprintf(" of %d candidate genes selected (ntree = %d, mtry = %d, %s).", r$n_input, r$rf_ntree, r$rf_mtry, rule_txt))
      })
      output[[paste0(sex_label, "_rf_plot")]] <- renderPlot({
        r <- res()
        validate(need(!is.null(r), not_yet_msg()))
        df <- head(data.frame(gene = names(r$gini), importance = as.numeric(r$gini), stringsAsFactors = FALSE), 20)
        p <- ggplot(df, aes(x = reorder(gene, importance), y = importance)) +
          geom_col(fill = sex_color) +
          coord_flip() + labs(x = NULL, y = "Mean decrease in Gini") + theme_arthomix(base_size = 12)
        # only draw the mean-Gini cutoff when it's the active selection rule
        if (!identical(r$rf_selection_rule, "top_n")) {
          p <- p + geom_hline(yintercept = r$gini_thr, linetype = "dashed", color = "#8A929C")
        }
        p
      })
      output[[paste0(sex_label, "_rf_table")]] <- DT::renderDataTable({
        r <- res()
        validate(need(!is.null(r), not_yet_msg()))
        df <- data.frame(gene = names(r$gini), gini_importance = round(as.numeric(r$gini), 4), stringsAsFactors = FALSE)
        df$selected <- df$gene %in% r$rf_genes
        DT::datatable(df, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
      })
      output[[paste0(sex_label, "_rf_download")]] <- downloadHandler(
        filename = function() sprintf("%s_random_forest_genes.csv", sex_label),
        content = function(file) {
          r <- res()
          df <- data.frame(gene = names(r$gini), gini_importance = as.numeric(r$gini), selected = names(r$gini) %in% r$rf_genes)
          df$gene <- tx_csv_safe(df$gene)
          write.csv(df, file, row.names = FALSE)
        }
      )

      # SVM-RFE
      output[[paste0(sex_label, "_svm_summary")]] <- renderUI({
        r <- res()
        if (is.null(r)) return(not_yet_note())
        k_txt <- if (identical(r$svm_panel_mode, "manual")) {
          sprintf("manual k = %d", length(r$svm_genes))
        } else {
          sprintf("CV-optimal k = %d, CV error = %.3f", r$svm_curve$best, r$svm_curve$besterr)
        }
        p(strong(length(r$svm_genes)), sprintf(" of %d candidate genes selected (cost = %s, %s).", r$n_input, r$svm_cost, k_txt))
      })
      output[[paste0(sex_label, "_svm_plot")]] <- renderPlot({
        r <- res()
        validate(need(!is.null(r), not_yet_msg()))
        df <- data.frame(k = r$svm_curve$k, error = r$svm_curve$err)
        # marks whichever panel size was actually used - CV-optimal k or manual k
        marker_k <- if (identical(r$svm_panel_mode, "manual")) length(r$svm_genes) else r$svm_curve$best
        ggplot(df, aes(x = k, y = error)) +
          geom_line(color = sex_color) + geom_point(color = sex_color) +
          geom_vline(xintercept = marker_k, linetype = "dashed", color = "#8A929C") +
          labs(x = "Top-k ranked genes", y = "10-fold CV classification error") + theme_arthomix(base_size = 12)
      })
      output[[paste0(sex_label, "_svm_table")]] <- DT::renderDataTable({
        r <- res()
        validate(need(!is.null(r), not_yet_msg()))
        df <- data.frame(gene = r$svm_rank, rank = seq_along(r$svm_rank), stringsAsFactors = FALSE)
        df$selected <- df$gene %in% r$svm_genes
        DT::datatable(df, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
      })
      output[[paste0(sex_label, "_svm_download")]] <- downloadHandler(
        filename = function() sprintf("%s_svm_rfe_genes.csv", sex_label),
        content = function(file) {
          r <- res()
          df <- data.frame(gene = r$svm_rank, rank = seq_along(r$svm_rank), selected = r$svm_rank %in% r$svm_genes)
          df$gene <- tx_csv_safe(df$gene)
          write.csv(df, file, row.names = FALSE)
        }
      )

      # Consensus/Overlap - live-reactive to the "Methods to intersect" checkboxes
      # (input$consensus_methods) so picking LASSO only, any two, or all three
      # re-intersects the already-fitted r$sets instantly, with no re-run needed.
      method_labels <- c(LASSO = "LASSO", RandomForest = "Random Forest", SVM_RFE = "SVM-RFE")
      consensus_used_methods <- function(r) {
        chosen <- intersect(input$consensus_methods %||% names(r$sets), names(r$sets))
        if (length(chosen) == 0) names(r$sets) else chosen
      }
      consensus_used_genes <- function(r, used) Reduce(intersect, r$sets[used])
      consensus_venn_obj <- reactive({
        r <- res()
        validate(need(!is.null(r), not_yet_msg()))
        used <- consensus_used_methods(r)
        genes <- consensus_used_genes(r, used)
        draw_overlap_venn(r$sets[used], title = sprintf("%s: %d-gene consensus", tools::toTitleCase(sex_label), length(genes)), fill_high = sex_color)
      })
      ## Which Venn region (if any) the user's last click on this plot landed in -
      ## a pure reactive (not a stored reactiveVal), so it's always consistent with
      ## whatever's currently drawn: it re-resolves automatically if the "Methods to
      ## intersect" checkboxes change the diagram, and reading input$..._plot_click
      ## before any click throws (caught below) rather than resolving to a stale value.
      consensus_clicked_region <- reactive({
        click <- input[[paste0(sex_label, "_consensus_plot_click")]]
        req(click)
        r <- res(); req(r)
        used <- consensus_used_methods(r)
        venn_region_at_point(r$sets[used], click$x, click$y)
      })
      output[[paste0(sex_label, "_consensus_summary")]] <- renderUI({
        r <- res()
        if (is.null(r)) return(not_yet_note())
        used <- consensus_used_methods(r)
        genes <- consensus_used_genes(r, used)
        clicked <- tryCatch(consensus_clicked_region(), error = function(e) NULL)
        tagList(
          p(strong(length(genes)), sprintf(" genes selected by %s: ", paste(unname(method_labels[used]), collapse = " ∩ ")),
            if (length(genes) > 0) paste(genes, collapse = ", ") else "none", "."),
          if (!is.null(clicked)) {
            div(class = "empty-note", icon("filter"),
                sprintf("Table filtered to the \"%s\" region (%d gene%s): %s. Click elsewhere on the diagram to show all genes again.",
                        gsub("/", " ∩ ", clicked$name, fixed = TRUE), clicked$count, if (clicked$count == 1) "" else "s",
                        if (clicked$count > 0) paste(clicked$item, collapse = ", ") else "none"))
          },
          p(class = "submodule-desc", r$candidate_note)
        )
      })
      output[[paste0(sex_label, "_consensus_plot")]] <- renderPlot({ consensus_venn_obj() })
      output[[paste0(sex_label, "_consensus_table")]] <- DT::renderDataTable({
        r <- res()
        validate(need(!is.null(r), not_yet_msg()))
        used <- consensus_used_methods(r)
        genes <- consensus_used_genes(r, used)
        df <- data.frame(gene = union(union(r$lasso_genes, r$rf_genes), r$svm_genes), stringsAsFactors = FALSE)
        df$lasso <- df$gene %in% r$lasso_genes
        df$random_forest <- df$gene %in% r$rf_genes
        df$svm_rfe <- df$gene %in% r$svm_genes
        df$consensus <- df$gene %in% genes
        clicked <- tryCatch(consensus_clicked_region(), error = function(e) NULL)
        if (!is.null(clicked)) df <- df[df$gene %in% clicked$item, , drop = FALSE]
        df <- df[order(-df$consensus, df$gene), ]
        DT::datatable(df, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
      })
      output[[paste0(sex_label, "_consensus_download")]] <- downloadHandler(
        filename = function() sprintf("%s_consensus_genes.csv", sex_label),
        content = function(file) {
          r <- res()
          used <- consensus_used_methods(r)
          genes <- consensus_used_genes(r, used)
          df <- data.frame(gene = union(union(r$lasso_genes, r$rf_genes), r$svm_genes), stringsAsFactors = FALSE)
          df$lasso <- df$gene %in% r$lasso_genes
          df$random_forest <- df$gene %in% r$rf_genes
          df$svm_rfe <- df$gene %in% r$svm_genes
          df$consensus <- df$gene %in% genes
          df$gene <- tx_csv_safe(df$gene)
          write.csv(df, file, row.names = FALSE)
        }
      )
    }

    register_sex_technique_outputs("female", res_sex("female"))
    register_sex_technique_outputs("male", res_sex("male"))
    register_sex_technique_outputs("pooled", res_sex("pooled"))
  })
}
