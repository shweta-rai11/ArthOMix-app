## R/mod_featureselection.R
## Submodule: Feature Selection (Section 2.8)
## "Your analysis" fits three independent, sex-stratified feature-selection
## algorithms - LASSO logistic regression, random forest importance, and
## SVM-RFE - on a candidate gene panel, and reports the 3-way consensus.
## Female and male are ALWAYS modeled completely separately (two independent
## fits per method, never a combined model with sex as a covariate) - this
## mirrors this project's own methodology (Chapter_2_subchapter2_
## sexstratified.md Section 2.8) exactly, including its hyperparameter
## tuning (10-fold CV for LASSO's lambda, random forest's mtry, and SVM-RFE's
## cost) and its selection rules (LASSO: non-zero coefficient at lambda.min;
## RF: Mean Decrease in Gini above the mean; SVM-RFE: the top-k ranked
## features minimising 10-fold CV error), and its single global seed (1234)
## reset before every stochastic step.
##
## Three ways to get to a candidate gene panel (radioButtons "data_source"):
##   "project" (default) - use this session's own live Candidate Gene
##       Identification output (results$candidates$female/$male, i.e. WGCNA
##       module background intersected with sex-stratified DEGs) if it's
##       been run, otherwise fall back to this project's own bundled,
##       MR-prioritised candidate lists (FS_input_{female,male}.csv) as a
##       worked example - so this tab is never empty on first use.
##   "expr" - upload your own expression matrix + metadata (own fileInputs,
##       not the app-wide Dataset tab, since the sex column mapping here is
##       mandatory rather than optional) and pick candidate genes from it
##       directly (most variable, or a pasted list).
##   "deg"  - upload your own female/male DEG (or any ranked/significant
##       gene list) CSVs directly; expression values for those genes still
##       come from whatever dataset is currently loaded in the app.

## -----------------------------------------------------------------------
## SVM-RFE helpers - ports of this project's own script (goal2_sex_
## stratified/12_feature_selection.R) verbatim: recursive elimination by
## smallest squared linear-SVM weight, one feature per round, giving a full
## most-to-least-important ranking; then a 10-fold CV error curve over the
## top-k ranked features to pick the panel size objectively.
## -----------------------------------------------------------------------

## kernel is deliberately fixed to "linear" throughout SVM-RFE (both here and
## in fs_svm_rfe_curve() below) - the per-round elimination rule (dropping
## the feature with the smallest squared linear weight, w2 <- ((t(m$coefs)
## %*% m$SV)[1, ])^2) is only mathematically valid for a linear decision
## boundary; a radial/polynomial kernel has no single feature-wise weight to
## rank by. tolerance and class.weights (see fs_class_weight_levels() above)
## are still meaningful for a linear-kernel SVM and are exposed as advanced
## parameters.
fs_svm_rfe_rank <- function(X, y, cost = 1, tolerance = 0.001, class_weights = NULL) {
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
  c(feats, ranking)
}

fs_svm_rfe_curve <- function(X, y, rank, cost = 1, seed = 1234, folds = 10, tolerance = 0.001, class_weights = NULL) {
  ks <- seq_along(rank)
  err <- vapply(ks, function(k) {
    set.seed(seed)
    acc <- e1071::svm(X[, rank[seq_len(k)], drop = FALSE], y, kernel = "linear",
                       scale = TRUE, cost = cost, cross = folds,
                       tolerance = tolerance, class.weights = class_weights)$tot.accuracy
    1 - acc / 100
  }, numeric(1))
  list(k = ks, err = err, best = ks[which.min(err)], besterr = min(err))
}

FS_SVM_COST_GRID <- c(0.01, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16)

## SVM-RFE eliminates one feature per round via a full model refit
## (fs_svm_rfe_rank), then re-fits again per round to build the CV-error
## curve (fs_svm_rfe_curve) - roughly 2p linear-SVM fits for p candidate
## genes. This project's own MR-prioritised panels are 25-40 genes (a few
## minutes at most); this project's Candidate Gene Identification output,
## by contrast, can be a raw WGCNA-module-sized list of 1,000-2,000+ genes
## - at that size this is thousands of SVM fits and would run for hours.
## Candidate sets above this cap are reduced to their most variable genes
## (within that sex's samples) before fitting, in fs_build_sex() below -
## noted on-screen (candidate_note) so it's never a silent change.
FS_MAX_CANDIDATE_GENES <- 200

## Every entry here is the exact value this project's own script
## (goal2_sex_stratified/12_feature_selection.R) uses - fs_fit_sex()'s
## `params` argument is this list, with any subset overridden by the
## per-tab "Customize parameters" UI (each of LASSO/Random Forest/SVM-RFE
## has its own independent cv_folds and settings). Leaving every tab
## uncustomized reproduces the thesis panels exactly - and, on the default
## "project pipeline" data source, is served instantly from this project's
## own precomputed run rather than recomputed (see load_precomputed_fs()
## in the server below); nothing here changes what "the default pipeline"
## means.
FS_DEFAULT_PARAMS <- list(
  ## Class weighting - shared across all three methods (see fs_class_weight_
  ## levels()/fs_obs_weights() below), one control for the whole per-sex run
  ## rather than three duplicated ones. "equal" = this project's own
  ## methodology (no reweighting); "balanced" = auto inverse-frequency;
  ## "manual" = a fixed comparison:reference weight ratio.
  class_weight_mode = "equal", class_weight_ratio = 1,
  lasso_cv_folds = 10, lasso_alpha = 1, lasso_lambda_choice = "lambda.min",
  lasso_nlambda = 100, lasso_type_measure = "deviance",
  rf_cv_folds = 10, rf_ntree = 1000, rf_mtry_mode = "auto", rf_mtry_manual = NULL,
  rf_selection_rule = "above_mean", rf_top_n = 10,
  rf_nodesize = 1, rf_maxnodes = NULL,
  svm_cv_folds = 10, svm_cost_mode = "auto", svm_cost_manual = 1, svm_cost_grid = FS_SVM_COST_GRID,
  svm_panel_mode = "auto", svm_manual_k = 10, svm_tolerance = 0.001,
  ## Which of the three methods count toward "consensus" - all three
  ## (this project's own methodology) by default, but overridable to any
  ## non-empty subset, e.g. c("LASSO", "SVM_RFE") to replicate a published
  ## two-method overlap (Chen et al. 2021/2022's own Fig. 3c intersects only
  ## LASSO and SVM-RFE, never a random forest step).
  consensus_methods = c("LASSO", "RandomForest", "SVM_RFE")
)

## Named-by-level (reference, comparison) weights for imbalanced groups -
## a duplicate of mod_diagnostic.R's own diag_class_weight_levels()/diag_
## obs_weights(), per this codebase's per-module self-containment
## convention (see wgcna_module_pick_ui_builder()'s own note on this).
## Shared by randomForest's classwt and e1071's class.weights directly, and
## (converted to a per-observation vector by fs_obs_weights()) glmnet's
## weights=. "equal" (default) returns all-1s, i.e. this project's own
## unweighted methodology.
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
## Column names are made syntactically safe first (some gene symbols contain
## "-", e.g. HLA-A) and translated back to symbols on the way out - same
## precaution this project's own script takes. `params` (see
## FS_DEFAULT_PARAMS above) lets each method's hyperparameters and
## selection rule be overridden from the UI; any field left out falls back
## to the thesis default.
fs_fit_sex <- function(X, y, params = list()) {
  params <- utils::modifyList(FS_DEFAULT_PARAMS, params)
  GLOBAL_SEED <- 1234
  p <- ncol(X)
  min_class <- min(table(y))
  nf_lasso <- max(2, min(params$lasso_cv_folds, min_class))
  nf_rf <- max(2, min(params$rf_cv_folds, min_class))
  nf_svm <- max(2, min(params$svm_cv_folds, min_class))

  safe <- make.names(colnames(X), unique = TRUE)
  lk <- stats::setNames(colnames(X), safe)
  colnames(X) <- safe
  back <- function(v) unname(lk[v])

  ## Class weighting (see fs_class_weight_levels() above) - shared by all
  ## three methods below. class_weight_mode = "equal" (the default) makes
  ## obs_w all 1s and cw_levels c(1, 1), numerically identical to each
  ## method's unweighted fit, so this project's own numeric replication is
  ## unaffected unless class weighting is actively turned on.
  cw_levels <- fs_class_weight_levels(y, params$class_weight_mode, params$class_weight_ratio)
  obs_w <- fs_obs_weights(y, params$class_weight_mode, params$class_weight_ratio)

  ## (1) LASSO logistic regression. alpha (1 = pure LASSO, 0 = ridge),
  ## which lambda to read coefficients at, and the CV fold count are all
  ## user-overridable independently of the other two methods. nlambda and
  ## the CV metric (deviance/AUC/misclassification error) are additional
  ## advanced controls.
  set.seed(GLOBAL_SEED)
  cv <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = params$lasso_alpha, nfolds = nf_lasso,
                           type.measure = params$lasso_type_measure %||% "deviance", nlambda = params$lasso_nlambda,
                           weights = obs_w)
  lambda_s <- if (identical(params$lasso_lambda_choice, "lambda.1se")) "lambda.1se" else "lambda.min"
  co <- coef(cv, s = lambda_s)[-1, 1, drop = TRUE]
  lasso_genes <- back(names(co)[co != 0])

  ## (2) Random forest importance. mtry is either tuned by its own nf-fold
  ## CV grid search (caret, ROC metric, default) or a fixed manual value;
  ## ntree and the selection rule (above-mean Gini, default, or a fixed
  ## top-N) are both user-overridable. nodesize/maxnodes (tree complexity)
  ## and classwt (class weighting) are additional advanced controls.
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

  ## (3) SVM-RFE, linear kernel (see fs_svm_rfe_rank()'s own comment on why
  ## the kernel itself isn't user-selectable here). Cost is either tuned by
  ## its own nf-fold CV grid search (e1071::tune, default, over a
  ## user-editable grid) or a fixed manual value; the final panel size is
  ## either the top-k minimising nf-fold CV error (default) or a fixed
  ## manual k. tolerance and class.weights are additional advanced controls.
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

  ## Consensus = intersection of only the user-selected subset of methods
  ## (params$consensus_methods, default all three) - falls back to all
  ## three if somehow left empty, so "no methods ticked" can't silently
  ## zero out the panel.
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
  description = "LASSO, random forest and SVM-RFE consensus feature selection, by sex.",
  icon = "sliders"
)

## Wraps one sex's technique box so it only renders once THAT sex's OWN Run
## button has been clicked - not "any of the three", the way the outer
## fs_results_wrap toggle (shown on the first click of any button) already
## works. Applied to every Female/Male/Pooled box across all four tabs
## (LASSO/Random Forest/SVM-RFE/Overlap), so e.g. clicking "Run Female"
## reveals only the Female boxes; Male and Pooled stay hidden until their
## own buttons are clicked, instead of all three appearing together.
fs_sex_panel <- function(ns, run_btn_id, ...) {
  conditionalPanel(condition = sprintf("input['%s'] > 0", ns(run_btn_id)), ...)
}

## One technique's Female/Male pair of result boxes. `body` is a tagList of
## the sex-specific outputs (plot/summary/table/download) for that
## technique, built by the caller - written once and instantiated 4 times
## (LASSO / Random Forest / SVM-RFE / Consensus) rather than laid out
## separately for each.
mod_featureselection_technique_panel <- function(ns, prefix, title, plot_height = 300) {
  box(
    width = NULL, title = title, status = "primary", solidHeader = FALSE,
    withSpinner(uiOutput(ns(paste0(prefix, "_summary"))), color = "#2563EB", type = 6),
    withSpinner(plotOutput(ns(paste0(prefix, "_plot")), height = plot_height), color = "#2563EB", type = 6),
    div(class = "table-toolbar", downloadButton(ns(paste0(prefix, "_download")), "Genes (CSV)", class = "btn-sm")),
    DT::dataTableOutput(ns(paste0(prefix, "_table")))
  )
}

## One method's own filter controls, always visible and live-editable -
## no "customize" checkbox to find first. Every control's own default
## value (set on the numericInput/sliderInput/radioButtons themselves,
## below) is exactly this project's own methodology; change any of them
## to run with different settings. fs_any_customized() in the server
## compares the current values of these inputs against those same
## defaults to decide whether a run still qualifies for the instant
## precomputed path (see load_precomputed_fs()) or needs to be live.
## Each of the three methods gets its own instance - not one shared box -
## so folds/settings for one method never entangle with another's.
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
        arthochat_shortcut_ui(
          "New to feature selection? Ask ArthOChat.",
          compact = TRUE
        ),
        box(
          width = NULL, title = "Candidate genes & samples", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Female and male fit separately by default; \"Run All\" pools every sample."),
          ## Run controls pinned at the top of the box (position: sticky, see
          ## .fs-run-bar in custom.css) rather than after every data-source/
          ## group setting below - previously they sat at the very end of a
          ## long form, off-screen until scrolled past everything. Sticky
          ## keeps them reachable the whole time you're adjusting settings
          ## further down, without duplicating the actionButtons (which
          ## would need a second input id and a merged trigger).
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
            ## Same results$wgcna$module_genes source as the "deg" data
            ## source's own WGCNA-module option below, just intersected with
            ## THIS uploaded expression matrix's own genes instead of the
            ## app-wide Dataset tab's - a separate selectInput/output
            ## (wgcna_module_pick_expr) rather than reusing the "deg" one,
            ## so both can sit in the DOM at once without a duplicate
            ## element ID (only one is ever visible at a time, per
            ## data_source, but Shiny still renders both conditionalPanels).
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
            ## Reads results$wgcna$module_genes directly - the same
            ## color -> gene-vector split mod_wgcna.R's Step 3 already
            ## populates the shared results object with (see the
            ## observeEvent(net_result(), ...) block there) the moment WGCNA
            ## has been run this session, on whatever dataset is currently
            ## active - no CSV export/re-upload round-trip needed. The SAME
            ## module selection is used for Female, Male and "Run All
            ## (pooled)" alike (unlike the per-sex uploaded-file mode above,
            ## a WGCNA module isn't sex-specific by construction).
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
        ## Hidden (not absent - see shinyjs::show(ns("fs_results_wrap")) in
        ## the server, fired the same moment fs_has_run() first flips TRUE)
        ## rather than gated behind a req(fs_has_run())-style renderUI: the
        ## withSpinner()-wrapped outputs inside still need to already exist
        ## in the DOM, bound, the moment "Run Female"/"Run Male" is first
        ## clicked (see the comment on tabsetPanel below) - shinyjs::hidden()
        ## keeps that binding intact while still showing nothing at all
        ## until that first click.
        shinyjs::hidden(div(
        id = ns("fs_results_wrap"),
        tabsetPanel(
          id = ns("technique_tabs"), type = "tabs",
          ## Every technique/consensus box below is written directly into
          ## this static UI - NOT behind a server-side req(res())-gated
          ## renderUI - so its withSpinner()-wrapped outputs already exist
          ## in the DOM the moment "Run Female"/"Run Male" is clicked. Only
          ## the outputs THEMSELVES (bound below in the server, e.g.
          ## female_lasso_summary) are invalidated by that click; the boxes
          ## around them are not, so shinycssloaders' busy-detection (which
          ## only has a visible effect on already-existing output elements)
          ## actually gets something to attach the spinner to. This mirrors
          ## the working pattern already used elsewhere in this app (e.g.
          ## mod_wgcna.R's Step 3 box, which is never gated on its own Run
          ## button either). Each output's own render function still starts
          ## from "not run yet" and only shows real content once its sex has
          ## a result - see not_yet_note() in the server below.
          tabPanel(
            "LASSO", br(),
            uiOutput(ns("lasso_params_ui")),
            fluidRow(
              column(6, fs_sex_panel(ns, "run_female_btn", mod_featureselection_technique_panel(ns, "female_lasso", "Female - LASSO"))),
              column(6, fs_sex_panel(ns, "run_male_btn", mod_featureselection_technique_panel(ns, "male_lasso", "Male - LASSO")))
            ),
            fluidRow(column(12, fs_sex_panel(ns, "run_pooled_btn", mod_featureselection_technique_panel(ns, "pooled_lasso", "Pooled (all) - LASSO"))))
          ),
          tabPanel(
            "Random Forest", br(),
            uiOutput(ns("rf_params_ui")),
            fluidRow(
              column(6, fs_sex_panel(ns, "run_female_btn", mod_featureselection_technique_panel(ns, "female_rf", "Female - Random Forest"))),
              column(6, fs_sex_panel(ns, "run_male_btn", mod_featureselection_technique_panel(ns, "male_rf", "Male - Random Forest")))
            ),
            fluidRow(column(12, fs_sex_panel(ns, "run_pooled_btn", mod_featureselection_technique_panel(ns, "pooled_rf", "Pooled (all) - Random Forest"))))
          ),
          tabPanel(
            "SVM-RFE", br(),
            uiOutput(ns("svm_params_ui")),
            fluidRow(
              column(6, fs_sex_panel(ns, "run_female_btn", mod_featureselection_technique_panel(ns, "female_svm", "Female - SVM-RFE"))),
              column(6, fs_sex_panel(ns, "run_male_btn", mod_featureselection_technique_panel(ns, "male_svm", "Male - SVM-RFE")))
            ),
            fluidRow(column(12, fs_sex_panel(ns, "run_pooled_btn", mod_featureselection_technique_panel(ns, "pooled_svm", "Pooled (all) - SVM-RFE"))))
          ),
          tabPanel(
            "Overlap", br(),
            uiOutput(ns("consensus_params_ui")),
            fluidRow(
              column(6, fs_sex_panel(ns, "run_female_btn", mod_featureselection_technique_panel(ns, "female_consensus", "Female - Overlap", plot_height = 340))),
              column(6, fs_sex_panel(ns, "run_male_btn", mod_featureselection_technique_panel(ns, "male_consensus", "Male - Overlap", plot_height = 340)))
            ),
            fluidRow(column(12, fs_sex_panel(ns, "run_pooled_btn", mod_featureselection_technique_panel(ns, "pooled_consensus", "Pooled (all) - Overlap", plot_height = 340))))
          )
        ),
        box(
          width = NULL, title = "Result", status = "primary", solidHeader = FALSE,
          withSpinner(uiOutput(ns("summary_ui")), color = "#2563EB", type = 6)
        )
        ))
      )
    ),
    uiOutput(ns("references_box_ui"))
  )
}

mod_featureselection_server <- function(id, dataset, results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## -----------------------------------------------------------------
    ## "Upload my own expression data" - a self-contained upload, separate
    ## from the app-wide Dataset tab, because the sex column mapping here is
    ## mandatory (sex-stratification is this whole module's point) rather
    ## than optional the way it is on the Dataset tab.
    ## -----------------------------------------------------------------

    own_meta_raw <- reactive({
      req(input$meta_file)
      path <- input$meta_file$datapath
      if (grepl("\\.rds$", input$meta_file$name, ignore.case = TRUE)) {
        d <- readRDS(path)
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

    ## Unified expr/meta source: "project" and "deg" both read the app-wide
    ## `dataset` (whatever's loaded on the Dataset tab); "expr" reads this
    ## module's own upload above.
    source_expr_meta <- reactive({
      if (identical(input$data_source, "expr")) {
        req(input$expr_file, input$meta_file, input$map_id, input$map_group, input$map_sex)
        expr <- if (grepl("\\.rds$", input$expr_file$name, ignore.case = TRUE)) {
          readRDS(input$expr_file$datapath)
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
        ## No sex-column check here any more - it only matters for the
        ## Female/Male buttons (checked in fs_build_sex() itself, right
        ## before it's actually needed) and would otherwise block "Run All
        ## (pooled)" on a dataset that legitimately has no sex data, e.g.
        ## one merged from sources where only some contributed a sex column.
        list(expr = dataset$expr, meta = dataset$meta)
      }
    })

    ## Which value in the sex column means "female" / "male" - matches
    ## whatever's actually in the loaded metadata (this app's own datasets
    ## harmonise to single-letter "F"/"M", but an upload might spell it out).
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

    ## -----------------------------------------------------------------
    ## Candidate gene panel per sex, one function per data source.
    ## -----------------------------------------------------------------

    ## "project": this session's own live Candidate Gene Identification
    ## output if available, else this project's own bundled, MR-prioritised
    ## candidate lists (with an MHC-excluded variant when the checkbox below
    ## is ticked) - so the tab is never empty on first use.
    project_candidate_genes <- function(sex_label) {
      live <- results$candidates[[sex_label]]$genes
      if (!is.null(live) && length(live) >= 2) {
        return(list(genes = live, is_live = TRUE,
                    note = sprintf("%d genes from this session's live %s candidate panel.",
                                   length(live), sex_label)))
      }
      fname <- sprintf("FS_input_%s%s.csv", sex_label, if (isTRUE(input$mhc_exclude)) "_noMHC" else "")
      bundled <- read_table_safe(fname)
      if (!is.null(bundled) && nrow(bundled) >= 2 && "gene" %in% colnames(bundled)) {
        return(list(genes = unique(as.character(bundled$gene)), is_live = FALSE,
                    note = sprintf("%d genes from the bundled %s candidate list (%s).",
                                   nrow(bundled), sex_label, fname)))
      }
      list(genes = character(0), is_live = FALSE,
           note = sprintf("No %s candidate genes available.", sex_label))
    }

    ## -----------------------------------------------------------------
    ## Fast path: this project's own precomputed LASSO/RF/SVM-RFE run
    ## (data/processed/new/ml_features.rds, or ml_features_noMHC.rds),
    ## already verified to reproduce fs_fit_sex()'s live output exactly on
    ## this project's own bundled candidates. Used instead of recomputing
    ## whenever the default pipeline is selected uncustomized (see
    ## use_fast_path in build_sex() below) - the same "instant on the
    ## project's own data, live otherwise" pattern mod_wgcna.R and
    ## mod_mr.R already use elsewhere in this app.
    ## -----------------------------------------------------------------

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

    ## "deg": either a user-uploaded gene list/DEG table (one file per sex,
    ## plus a dedicated pooled file for "Run All"), or - deg_source_mode ==
    ## "wgcna" - the gene list straight from a WGCNA module run this session,
    ## the SAME module regardless of sex_label since a co-expression module
    ## isn't sex-specific by construction (unlike the per-sex uploaded files).
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

    ## Module-color dropdown, shared builder for both entry points into
    ## results$wgcna$module_genes (set by mod_wgcna.R's Step 3 the moment a
    ## run completes) - the "deg" data source's own WGCNA-module mode
    ## (input_id = "wgcna_module_pick", reads expression from the app-wide
    ## Dataset tab) and the "expr" data source's matching gene_source choice
    ## (input_id = "wgcna_module_pick_expr", reads expression from this
    ## module's own uploaded matrix instead). "grey" is excluded from both,
    ## matching WGCNA's own Step 4/Step 5 treatment of it as an unassigned
    ## bin rather than a real co-expression module.
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
      has_live <- length(f_live) >= 2 && length(m_live) >= 2
      tagList(
        if (has_live) {
          div(class = "empty-note", icon("check"),
              sprintf("Using live candidates: %d female / %d male genes.", length(f_live), length(m_live)))
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

    ## -----------------------------------------------------------------
    ## Advanced method parameters - each of the three methods has its own
    ## "Customize" checkbox (in its own tab); only a checked method's
    ## fields are actually read, so customizing e.g. Random Forest alone
    ## leaves LASSO and SVM-RFE at this project's own defaults untouched.
    ## fs_fit_sex()'s utils::modifyList(FS_DEFAULT_PARAMS, ...) fills in
    ## defaults for anything not present here.
    ## -----------------------------------------------------------------

    ## Every LASSO/Random Forest/SVM-RFE control is always visible and
    ## live-editable now (no "customize" checkbox gating them) - so
    ## "customized" is decided by comparing the current value of each
    ## control against its own default (FS_DEFAULT_PARAMS), not by
    ## whether a checkbox is checked. use_fast_path in build_sex() below
    ## needs exactly this to decide whether the precomputed result is
    ## still valid to serve.
    val_eq <- function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-8))
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

    ## Tells the user, before they click Run, whether THIS click will be
    ## instant (precomputed) or live - and if live, exactly which of the
    ## fast-path conditions isn't met. Mirrors use_fast_path in build_sex()
    ## below exactly, so this is never wrong about what's about to happen.
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
      if (!isTRUE(grepl("^Example dataset:", dataset$source %||% ""))) {
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

    ## Reads every LASSO/Random Forest/SVM-RFE control's current value
    ## directly - no checkbox to gate them, so a control left untouched
    ## just reads back its own default (set on the widget itself in the
    ## UI) and a control the user changes takes effect on the next Run,
    ## same as any other filter in this app.
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

    ## -----------------------------------------------------------------
    ## Run: fit LASSO + random forest + SVM-RFE for both sexes on one
    ## click - whichever data source, per-tab parameters, and group
    ## contrast are currently set. Female and male are still always fit
    ## completely separately internally (fs_build_sex is called once per
    ## sex, independently), just triggered by the same button.
    ## -----------------------------------------------------------------

    fs_build_sex <- function(sex_label, sex_value) {
      req(input$ref_group, input$comp_group)
      validate(need(input$ref_group != input$comp_group, "Reference and comparison group must be different."))

      sem <- source_expr_meta()
      ## Only Female/Male need a usable sex column - "Run All" (sex_value
      ## NULL) pools every sample regardless, so it works even on a dataset
      ## that never had (or lost) sex information, e.g. one merged from
      ## sources where only some contributed a sex column.
      if (!is.null(sex_value)) {
        validate(need("sex" %in% colnames(sem$meta) && length(unique(stats::na.omit(sem$meta$sex))) >= 2,
                      "No usable sex column in this dataset - load one with sex data, or use \"Run All (pooled)\" instead."))
      }
      adv_params <- fs_advanced_params()
      any_customized <- fs_any_customized()

      ## Candidate gene lookup is cheap (a CSV read or an existing
      ## results$candidates list) regardless of path - only the model fit
      ## itself (fs_fit_sex, below) is the expensive part, so it's always
      ## resolved first to decide fast-path eligibility.
      cand_project <- if (identical(input$data_source, "project")) project_candidate_genes(sex_label) else NULL

      ## Fast path: default pipeline, this project's own unmodified bundled
      ## candidates (not live Candidate Gene Identification output, which
      ## wouldn't match the precomputed run), every method's parameters
      ## left uncustomized, the standard HC-vs-RA contrast this precomputed
      ## run was itself built on, AND the project's own default example
      ## dataset still actually loaded (same guard mod_wgcna.R uses for its
      ## own precomputed shortcut) - without this last check, someone who
      ## switched to a different dataset on the Dataset tab but still had
      ## HC/RA-labelled groups would silently get this project's numbers
      ## instead of their own.
      use_fast_path <- identical(input$data_source, "project") &&
        !isTRUE(cand_project$is_live) && !any_customized &&
        identical(input$ref_group, "HC") && identical(input$comp_group, "RA") &&
        isTRUE(grepl("^Example dataset:", dataset$source %||% ""))

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

      ## See FS_MAX_CANDIDATE_GENES above: a raw WGCNA-module-sized
      ## candidate set (thousands of genes) makes SVM-RFE's per-round
      ## refitting take hours - reduce to the most variable genes within
      ## this sex's own samples first, and say so on-screen.
      n_before_cap <- length(genes)
      if (n_before_cap > FS_MAX_CANDIDATE_GENES) {
        v <- apply(expr_sub[genes, , drop = FALSE], 1, stats::var)
        genes <- names(sort(v, decreasing = TRUE))[seq_len(FS_MAX_CANDIDATE_GENES)]
        cand$note <- sprintf("%s Reduced from %d to the %d most variable genes before fitting.",
                              cand$note, n_before_cap, FS_MAX_CANDIDATE_GENES)
      }

      X <- t(expr_sub[genes, , drop = FALSE])
      y <- factor(as.character(meta$group), levels = c(input$ref_group, input$comp_group))
      validate(need(all(table(y) >= 6), sprintf("Each group needs at least 6 %s samples for feature selection.", sex_label)))

      fit <- fs_fit_sex(X, y, params = adv_params)
      fit$candidate_note <- cand$note
      fit$n_ref <- sum(y == input$ref_group)
      fit$n_comp <- sum(y == input$comp_group)
      fit$ref_group <- input$ref_group; fit$comp_group <- input$comp_group
      fit$mhc_mode <- if (identical(input$data_source, "project")) (if (isTRUE(input$mhc_exclude)) "exclude" else "include") else "n/a"
      fit
    }

    ## ignoreInit = TRUE: without it, this evaluates once at session start
    ## (using run_female_btn/run_male_btn's initial value of 0), before the
    ## renderUI-driven ref_group/comp_group selectInputs have round-tripped
    ## to the server - req() inside fs_build_sex() then halts with a bare,
    ## unstyled error that briefly flashes in every technique panel until
    ## the real click.
    ##
    ## Female and male each have their own button and their own
    ## eventReactive, so one sex can be run (and re-run, with different
    ## per-tab parameters) without touching the other - fs_build_sex()
    ## itself is unchanged and still fits each sex completely independently
    ## internally; only the trigger is now split.
    fs_result_female <- eventReactive(input$run_female_btn, {
      fs_build_sex("female", sex_levels()$female)
    }, ignoreInit = TRUE)
    fs_result_male <- eventReactive(input$run_male_btn, {
      fs_build_sex("male", sex_levels()$male)
    }, ignoreInit = TRUE)
    ## "Run All": pools every sample regardless of sex - sex_value = NULL
    ## tells fs_build_sex() to skip the sex filter (and its sex-column
    ## requirement) entirely, so this works even when sex isn't available.
    ## For replicating a published method that never stratified by sex.
    fs_result_pooled <- eventReactive(input$run_pooled_btn, {
      fs_build_sex("pooled", NULL)
    }, ignoreInit = TRUE)

    ## Everything to the right of the "Candidate genes & samples" box -
    ## each method's own parameter/filter box, every Female/Male result
    ## panel (hidden via shinyjs, see "fs_results_wrap" in the UI above),
    ## the bottom Result box, and the References box - stays completely
    ## invisible until EITHER "Run Female" or "Run Male" has been clicked
    ## once. Stays TRUE for the rest of the session after that first click,
    ## so the parameter boxes remain visible/editable to adjust and re-run -
    ## only the very first reveal is gated.
    ## shinyjs functions called from inside moduleServer() already apply
    ## session$ns() to id= internally (shinyjs:::jsFuncHelper checks for a
    ## session_proxy) - passing an already-ns()-wrapped id here would get
    ## double-namespaced ("tx_featureselection-tx_featureselection-...")
    ## and silently match nothing client-side, so this has to be the RAW
    ## "fs_results_wrap", same as the raw id ns() wrapped once in the UI.
    fs_has_run <- reactiveVal(FALSE)
    observeEvent(input$run_female_btn, {
      fs_has_run(TRUE)
      shinyjs::show(id = "fs_results_wrap")
    }, ignoreInit = TRUE)
    observeEvent(input$run_male_btn, {
      fs_has_run(TRUE)
      shinyjs::show(id = "fs_results_wrap")
    }, ignoreInit = TRUE)
    observeEvent(input$run_pooled_btn, {
      fs_has_run(TRUE)
      shinyjs::show(id = "fs_results_wrap")
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

    ## Which methods count toward "consensus" - all three ticked by default
    ## (this project's own methodology: LASSO ∩ Random Forest ∩ SVM-RFE),
    ## but untick any one (or two) to intersect a smaller subset instead -
    ## e.g. LASSO + SVM-RFE only, to replicate a published two-method
    ## overlap (this exact choice is what Chen et al. 2021/2022's own Fig.
    ## 3c does: their hub-gene panel is LASSO ∩ SVM-RFE, no random forest
    ## step at all). At least one method must stay ticked; unticking the
    ## last one silently falls back to all three rather than erroring (see
    ## fs_fit_sex()'s own fallback), so this only ever narrows the panel,
    ## never breaks it. Changing this always forces a live re-run (see
    ## fs_any_customized()) since the fast-path precomputed result was only
    ## ever built as the 3-way intersection.
    output$consensus_params_ui <- renderUI({
      req(fs_has_run())
      mod_featureselection_params_box(
        ns, "consensus", "Overlap",
        "All three methods are intersected by default. Untick any to use a smaller subset, then re-run.",
        checkboxGroupInput(ns("consensus_methods"), "Methods to intersect",
                            choices = c("LASSO" = "LASSO", "Random Forest" = "RandomForest", "SVM-RFE" = "SVM_RFE"),
                            selected = c("LASSO", "RandomForest", "SVM_RFE"), inline = TRUE)
      )
    })

    output$references_box_ui <- renderUI({
      req(fs_has_run())
      box(
        width = 12, title = "References", status = "primary", solidHeader = FALSE,
        p(class = "submodule-desc", "Background reading for the methods used on this tab."),
        tags$ul(
          class = "dge-ref-list",
          tags$li(strong("LASSO (glmnet): "), "Friedman J, Hastie T, Tibshirani R (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. ", tags$em("Journal of Statistical Software"), ", 33(1)."),
          tags$li(strong("Random forests: "), "Breiman L (2001). Random Forests. ", tags$em("Machine Learning"), ", 45, 5-32; and Liaw A, Wiener M (2002). Classification and Regression by randomForest. ", tags$em("R News"), ", 2(3)."),
          tags$li(strong("SVM-RFE: "), "Guyon I, Weston J, Barnhill S, Vapnik V (2002). Gene Selection for Cancer Classification using Support Vector Machines. ", tags$em("Machine Learning"), ", 46, 389-422."),
          tags$li(strong("caret (hyperparameter tuning): "), "Kuhn M (2008). Building Predictive Models in R Using the caret Package. ", tags$em("Journal of Statistical Software"), ", 28(5).")
        ),
        p(class = "submodule-desc", strong("Ask ArthOChat"), " above for a plain-language walkthrough of any method on this page.")
      )
    })

    ## Female and male are saved into results$featureselection independently
    ## as each finishes - utils::modifyList() only ever overwrites its own
    ## sex's entry, so running Male alone (or re-running it later) never
    ## clobbers an already-saved Female result, and vice versa.
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
      res_f <- tryCatch(fs_result_female(), error = function(e) NULL)
      res_m <- tryCatch(fs_result_male(), error = function(e) NULL)
      res_p <- tryCatch(fs_result_pooled(), error = function(e) NULL)
      status_row <- function(sex, r) {
        if (is.null(r)) {
          tags$li(icon("circle-minus", style = "color: #8A929C;"), sprintf(" %s feature selection - not run yet", sex))
        } else {
          tags$li(icon("check", style = "color: #1a9c5f;"), strong(sprintf(" %s feature selection completed: ", sex)),
                  sprintf("%d consensus genes%s", length(r$consensus), if (isTRUE(r$fast_path)) " (instant)" else " (live)"))
        }
      }
      tagList(
        p(class = "submodule-desc", style = "margin-bottom: 4px;", "Status:"),
        tags$ul(style = "padding-left: 18px; margin-bottom: 0; list-style: none;",
                status_row("Female", res_f),
                status_row("Male", res_m),
                status_row("Pooled (all)", res_p))
      )
    })

    ## Female and male render their line in the Result box independently -
    ## whichever has actually been run, showing "not run yet" for the
    ## other rather than waiting on both. Contrast/MHC mode are read
    ## straight off the shared input controls (not off either sex's fit
    ## object) since both sexes are always fit against the same
    ## reference/comparison group and MHC setting.
    output$summary_ui <- renderUI({
      res_f <- tryCatch(fs_result_female(), error = function(e) NULL)
      res_m <- tryCatch(fs_result_male(), error = function(e) NULL)
      res_p <- tryCatch(fs_result_pooled(), error = function(e) NULL)
      mhc_mode <- if (identical(input$data_source, "project")) (if (isTRUE(input$mhc_exclude)) "exclude" else "include") else "n/a"
      sex_line <- function(sex, r) {
        if (is.null(r)) return(p(strong(sprintf("%s: ", sex)), "not run yet."))
        p(strong(sprintf("%s: ", sex)),
          if (isTRUE(r$fast_path)) icon("bolt") else NULL,
          sprintf(" %d candidate genes, %d samples (%d vs %d) → %d LASSO / %d random forest / %d SVM-RFE → %d consensus%s.",
                  r$n_input, r$n_samples, r$n_comp, r$n_ref,
                  length(r$lasso_genes), length(r$rf_genes), length(r$svm_genes), length(r$consensus),
                  if (isTRUE(r$fast_path)) " (instant, precomputed)" else " (live)"))
      }
      tagList(
        p(strong("Contrast: "), sprintf("%s vs %s", input$comp_group %||% "?", input$ref_group %||% "?"),
          if (!identical(mhc_mode, "n/a")) sprintf(" (MHC region %s)", if (identical(mhc_mode, "exclude")) "excluded" else "included") else NULL, "."),
        sex_line("Female", res_f), sex_line("Male", res_m), sex_line("Pooled (all)", res_p)
      )
    })

    ## Every per-technique output below reads through this: NULL (not a
    ## req()-halt) when that sex hasn't been run yet, so each output can
    ## explicitly show "not run yet" rather than sitting there blank.
    res_sex <- function(sex_label) reactive({
      fr <- switch(sex_label, female = fs_result_female, male = fs_result_male, pooled = fs_result_pooled)
      tryCatch(fr(), error = function(e) NULL)
    })

    register_sex_technique_outputs <- function(sex_label, res) {
      sex_color <- switch(sex_label, female = "#1a7a3c", male = "#7a4a26", pooled = "#2563EB")
      ## Each technique/consensus box now lives directly in the static UI
      ## (see mod_featureselection_ui above) rather than behind a
      ## req(res())-gated renderUI, so this note covers BOTH "never clicked
      ## yet" and "clicked but the run failed validation" - there's no
      ## separate state to distinguish them by once the box always exists.
      not_yet_note <- function() {
        div(class = "empty-note", icon("circle-info"),
            "No result yet. Click Run Female, Run Male, or Run All (pooled) on the left.")
      }
      sex_title <- tools::toTitleCase(sex_label)

      ## LASSO. Every renderUI summary below explicitly branches on
      ## is.null(r) to show not_yet_note() - the whole reason res_sex()
      ## returns a plain NULL rather than req()-halting. Plot/table/
      ## download outputs instead `req(r)` right after `r <- res()` (a
      ## plain statement, not nested in a generic-dispatch argument - see
      ## the plot() note below) so they just stay blank under the note.
      output[[paste0(sex_label, "_lasso_summary")]] <- renderUI({
        r <- res()
        if (is.null(r)) return(not_yet_note())
        lambda_used <- if (identical(r$lasso_lambda_choice, "lambda.1se")) r$cv$lambda.1se else r$cv$lambda.min
        p(strong(length(r$lasso_genes)),
          sprintf(" of %d candidate genes selected (alpha = %.2f, %s = %.4f).", r$n_input, r$lasso_alpha, r$lasso_lambda_choice, lambda_used))
      })
      ## `r <- res()` as its own statement, not `plot(res()$cv)` inline:
      ## plot() is an S3 generic, so passing res()$cv directly forces the
      ## halt to be evaluated INSIDE plot()'s own method-dispatch
      ## machinery, which surfaces an ugly error banner instead of a clean
      ## blank panel. Forcing res() first avoids that entirely.
      output[[paste0(sex_label, "_lasso_plot")]] <- renderPlot({
        r <- res()
        req(r)
        plot(r$cv)
      })
      output[[paste0(sex_label, "_lasso_table")]] <- DT::renderDataTable({
        r <- res()
        req(r)
        DT::datatable(data.frame(gene = r$lasso_genes, stringsAsFactors = FALSE), rownames = FALSE,
                      options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
      })
      output[[paste0(sex_label, "_lasso_download")]] <- downloadHandler(
        filename = function() sprintf("%s_lasso_genes.csv", sex_label),
        content = function(file) write.csv(data.frame(gene = res()$lasso_genes), file, row.names = FALSE)
      )

      ## Random Forest
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
        req(r)
        df <- head(data.frame(gene = names(r$gini), importance = as.numeric(r$gini), stringsAsFactors = FALSE), 20)
        p <- ggplot(df, aes(x = reorder(gene, importance), y = importance)) +
          geom_col(fill = sex_color) +
          coord_flip() + labs(x = NULL, y = "Mean decrease in Gini") + theme_arthomix(base_size = 12)
        ## Only draw the mean-Gini cutoff line when it's actually the
        ## selection rule in effect - misleading otherwise (top-N mode
        ## picks genes by rank, not by this threshold).
        if (!identical(r$rf_selection_rule, "top_n")) {
          p <- p + geom_hline(yintercept = r$gini_thr, linetype = "dashed", color = "#8A929C")
        }
        p
      })
      output[[paste0(sex_label, "_rf_table")]] <- DT::renderDataTable({
        r <- res()
        req(r)
        df <- data.frame(gene = names(r$gini), gini_importance = round(as.numeric(r$gini), 4), stringsAsFactors = FALSE)
        df$selected <- df$gene %in% r$rf_genes
        DT::datatable(df, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
      })
      output[[paste0(sex_label, "_rf_download")]] <- downloadHandler(
        filename = function() sprintf("%s_random_forest_genes.csv", sex_label),
        content = function(file) {
          r <- res()
          df <- data.frame(gene = names(r$gini), gini_importance = as.numeric(r$gini), selected = names(r$gini) %in% r$rf_genes)
          write.csv(df, file, row.names = FALSE)
        }
      )

      ## SVM-RFE
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
        req(r)
        df <- data.frame(k = r$svm_curve$k, error = r$svm_curve$err)
        ## Marks wherever the panel size actually came from - the
        ## CV-error-minimising k (default) or the manual k, if that mode
        ## was chosen instead.
        marker_k <- if (identical(r$svm_panel_mode, "manual")) length(r$svm_genes) else r$svm_curve$best
        ggplot(df, aes(x = k, y = error)) +
          geom_line(color = sex_color) + geom_point(color = sex_color) +
          geom_vline(xintercept = marker_k, linetype = "dashed", color = "#8A929C") +
          labs(x = "Top-k ranked genes", y = "10-fold CV classification error") + theme_arthomix(base_size = 12)
      })
      output[[paste0(sex_label, "_svm_table")]] <- DT::renderDataTable({
        r <- res()
        req(r)
        df <- data.frame(gene = r$svm_rank, rank = seq_along(r$svm_rank), stringsAsFactors = FALSE)
        df$selected <- df$gene %in% r$svm_genes
        DT::datatable(df, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
      })
      output[[paste0(sex_label, "_svm_download")]] <- downloadHandler(
        filename = function() sprintf("%s_svm_rfe_genes.csv", sex_label),
        content = function(file) {
          r <- res()
          df <- data.frame(gene = r$svm_rank, rank = seq_along(r$svm_rank), selected = r$svm_rank %in% r$svm_genes)
          write.csv(df, file, row.names = FALSE)
        }
      )

      ## Consensus / Overlap - r$consensus_methods (default all three, see
      ## FS_DEFAULT_PARAMS$consensus_methods / the "Overlap" tab's checkbox
      ## group) is the authoritative set this whole box reflects; the Venn
      ## plot, "selected by" wording, and the table/download's own
      ## "consensus" column all key off it directly rather than re-deriving
      ## a hardcoded 3-way AND, so a 2-method selection (e.g. LASSO ∩
      ## SVM-RFE only) renders a genuine 2-circle Venn and a 2-method
      ## consensus column, not a 3-circle one with an unused RF ring.
      method_labels <- c(LASSO = "LASSO", RandomForest = "Random Forest", SVM_RFE = "SVM-RFE")
      consensus_venn_obj <- reactive({
        r <- res()
        req(r)
        used <- r$consensus_methods %||% names(r$sets)
        draw_overlap_venn(r$sets[used], title = sprintf("%s: %d-gene consensus", tools::toTitleCase(sex_label), length(r$consensus)), fill_high = sex_color)
      })
      output[[paste0(sex_label, "_consensus_summary")]] <- renderUI({
        r <- res()
        if (is.null(r)) return(not_yet_note())
        used <- r$consensus_methods %||% names(r$sets)
        tagList(
          p(strong(length(r$consensus)), sprintf(" genes selected by %s: ", paste(unname(method_labels[used]), collapse = " ∩ ")),
            if (length(r$consensus) > 0) paste(r$consensus, collapse = ", ") else "none", "."),
          p(class = "submodule-desc", r$candidate_note)
        )
      })
      output[[paste0(sex_label, "_consensus_plot")]] <- renderPlot({ consensus_venn_obj() })
      output[[paste0(sex_label, "_consensus_table")]] <- DT::renderDataTable({
        r <- res()
        req(r)
        df <- data.frame(gene = union(union(r$lasso_genes, r$rf_genes), r$svm_genes), stringsAsFactors = FALSE)
        df$lasso <- df$gene %in% r$lasso_genes
        df$random_forest <- df$gene %in% r$rf_genes
        df$svm_rfe <- df$gene %in% r$svm_genes
        df$consensus <- df$gene %in% r$consensus
        df <- df[order(-df$consensus, df$gene), ]
        DT::datatable(df, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
      })
      output[[paste0(sex_label, "_consensus_download")]] <- downloadHandler(
        filename = function() sprintf("%s_consensus_genes.csv", sex_label),
        content = function(file) {
          r <- res()
          df <- data.frame(gene = union(union(r$lasso_genes, r$rf_genes), r$svm_genes), stringsAsFactors = FALSE)
          df$lasso <- df$gene %in% r$lasso_genes
          df$random_forest <- df$gene %in% r$rf_genes
          df$svm_rfe <- df$gene %in% r$svm_genes
          df$consensus <- df$gene %in% r$consensus
          write.csv(df, file, row.names = FALSE)
        }
      )
    }

    register_sex_technique_outputs("female", res_sex("female"))
    register_sex_technique_outputs("male", res_sex("male"))
    register_sex_technique_outputs("pooled", res_sex("pooled"))
  })
}
