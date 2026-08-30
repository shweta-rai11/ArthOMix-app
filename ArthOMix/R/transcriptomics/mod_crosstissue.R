## R/mod_crosstissue.R
## Submodule: Cross-Tissue Validation (Section 2.11)
## "Your analysis" evaluates a user-chosen gene panel in the independent RA
## synovium dataset (GSE89408, val_synovium.rds): sex-stratified discovery
## (synovium log2FC/significance, direction concordance with blood) plus a
## full four-classifier panel model (logistic regression, elastic net,
## random forest, SVM) - the SAME four techniques and the SAME box/tab/KPI
## layout as Diagnostic Model (mod_diagnostic.R), so the two "feel" like one
## family of tools, but fit on this module's OWN data (synovium, not blood)
## and against this module's OWN methodology, not a copy of Diagnostic
## Model's train/test split.
##
## WHY THIS ISN'T A TRAIN/TEST SPLIT LIKE DIAGNOSTIC MODEL: this cohort has
## no natural held-out partition the way the blood training cohort does (see
## mod_diagnostic.R's header) - it is itself the held-out compartment. This
## project's own Section 2.11.5 methodology therefore fits each classifier
## ONCE on the full synovium sex-subset (an "apparent"/resubstitution AUC,
## reported explicitly as an optimistic upper bound, never as performance -
## Harrell, Lee and Mark 1996) and separately estimates discrimination by
## refitting across outer folds, pooling every out-of-fold prediction into
## ONE ROC curve before computing AUC + a DeLong confidence interval. This
## module reports BOTH that pooled estimate (this project's own headline
## number) AND, in the same "CV AUC by fold" bar chart Diagnostic Model
## uses, the per-fold AUC values themselves - format parity with Diagnostic
## Model without inventing a train/test split this cohort doesn't have.
##
## WHAT'S ALSO NOT COPIED: only the identity of the panel genes transfers
## from blood - never a blood-fitted model's coefficients (Section 2.11.1).
## Every classifier here is refit from scratch within synovium. The "Panel
## genes present in synovium" gene set, the synovium DE table, and the panel
## consensus lists (val$fsig / val$msig) all come from this module's own
## bundled val_synovium.rds - not from the blood expression matrix Diagnostic
## Model reads from `dataset`.
##
## THREE analyses, one Run button per sex (mirroring Diagnostic Model's
## per-sex Run, shared across its two tabs):
##   "Synovium Discovery & Concordance" - per-gene synovium log2FC/adjP
##       (sex-adjusted limma-voom DE, already computed and stored in
##       val_synovium.rds$tt), direction concordance against blood, and
##       per-gene AUC under BOTH orientation conventions this project's own
##       script emits explicitly (best-direction: AUC>=0.5 by construction,
##       "how much information does this gene carry"; train-fixed: oriented
##       by the blood direction, so AUC<0.5 means the association reversed
##       out of sample - the only convention valid for cross-dataset
##       comparison). See METHODS_2.11_crosstissue.md Section 2.11.6.
##   "Panel Classifier - Full Fit" - the four classifiers, apparent AUC.
##   "Panel Classifier - Cross-Validated" - the same four classifiers,
##       scored out-of-fold (pooled ROC + per-fold bar chart).
## A fourth, read-only tab, "Cross-Dataset Comparison", lines up this
## module's synovium AUCs against Diagnostic Model's OWN saved blood AUCs
## for the same sex (results$diagnostic, written by mod_diagnostic.R) when
## available this session - Section 2.11.7's "presentation across datasets",
## built from already-computed shared app state rather than any duplicated
## blood data or re-run blood model.
##
## HYPERPARAMETERS: elastic net alpha, random forest mtry and SVM cost are
## each tuned by an INNER cross-validation grid search on the full synovium
## sex-subset (identical tuning code to Diagnostic Model's, minus any
## train/test split) - each model type gets its own advanced-parameters box,
## switched to whichever model pill was clicked most recently in either
## sex's Full Fit or Cross-Validated tab. The OUTER cross-validation used for
## the pooled/per-fold performance estimate is a single, shared setting
## (fold count + stratification) for all four models at once, set in the
## "Advanced filters" box on the left - not a per-model choice, because it
## is what "performance in synovium" means here, independent of which
## classifier is being scored.
##
## FOLD STRATIFICATION: this project's own script (23_crosstissue_biomarker_
## discovery.R) assigns outer folds by simple random allocation, not
## stratified by disease status - a documented limitation (Section 2.11.5),
## since this cohort's disease/control imbalance (Normal is a small minority
## in both sexes) means an unstratified fold can end up with very few
## control samples. This module defaults to STRATIFIED folds (a genuine
## improvement over the published script) but keeps "simple random - matches
## this project's own script exactly" selectable, so the published numbers
## remain reproducible on demand.

mod_crosstissue_config <- list(
  id = "crosstissue", group = "Validation",
  title = "Cross-Tissue Validation",
  description = "Validation of the diagnostic model based on different tissue type and sex. Four-classifier panel model (logistic regression, elastic net, random forest, SVM). Performs analysis on both preloaded or uploaded data, based on sex.",
  icon = "shuffle"
)

## ---------------------------------------------------------------------------
## Shared fitting helpers (pure functions - no `input`/`results`/`session`,
## same split as mod_diagnostic.R's diag_* helpers vs. its server body).
## Reuses mod_diagnostic.R's diag_zrows()/diag_auc_ci()/diag_perf_at_cutoff()/
## diag_hyperparam_value()/diag_separation_note() directly - every mod_*.R
## file is sourced into the same shared environment (see submodules_
## registry.R's own header note on this), and none of the calls below happen
## at source time, only inside function bodies invoked once the whole app
## has finished loading, so file sourcing order does not matter here.
## ---------------------------------------------------------------------------

CT_SVM_COST_GRID <- c(0.01, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16)

CT_DEFAULT_PARAMS <- list(
  cv_folds = 10, stratified_folds = TRUE,
  enet_cv_folds = 5, enet_alpha_grid = c(0.1, 0.3, 0.5, 0.7, 0.9, 1.0), enet_lambda_choice = "lambda.min",
  rf_cv_folds = 5, rf_ntree = 1000, rf_mtry_mode = "auto", rf_mtry_manual = NULL,
  svm_cv_folds = 5, svm_kernel = "linear", svm_cost_mode = "auto", svm_cost_manual = 1, svm_cost_grid = CT_SVM_COST_GRID
)

CT_TECHNIQUES <- list(
  list(key = "lr", label = "Logistic Regression"),
  list(key = "enet", label = "Elastic Net"),
  list(key = "rf", label = "Random Forest"),
  list(key = "svm", label = "SVM")
)

ct_youden <- function(roc_obj) pROC::coords(roc_obj, "best", best.method = "youden",
                                             ret = c("threshold", "sensitivity", "specificity", "accuracy"), transpose = FALSE)

## A panel gene counts as a validated cross-tissue biomarker when all three
## hold: its direction of effect agrees between blood and synovium, its
## synovium differential expression is significant at the chosen FDR, and its
## synovium AUC (best-direction) clears a "more than weak" discrimination bar
## (Hosmer & Lemeshow 2013's conventional 0.70 cutoff). This is the single
## definition every KPI tile, plot and table column below reads from, so the
## three views can never disagree about which genes qualify.
CT_BIOMARKER_AUC_MIN <- 0.70

ct_biomarker_flag <- function(d, sig_cut) {
  !is.na(d$concordant) & d$concordant &
    !is.na(d$syn_adjP) & d$syn_adjP < sig_cut &
    !is.na(d$auc_bestdir) & d$auc_bestdir >= CT_BIOMARKER_AUC_MIN
}

## One gene's synovium AUC, best-direction convention (>= 0.5 by
## construction) - matches 20_testing_synovium_external.R's gene_stat().
ct_gene_auc <- function(values, y) {
  r <- tryCatch(pROC::roc(y, as.numeric(values), direction = "<", levels = levels(y), quiet = TRUE), error = function(e) NULL)
  if (is.null(r)) return(NA_real_)
  a <- as.numeric(pROC::auc(r))
  if (is.na(a)) return(NA_real_)
  if (a < 0.5) 1 - a else a
}

## Per-gene synovium discovery table for one sex: log2FC/adjP from the
## already-computed, sex-adjusted limma-voom DE (val$tt - every synovium
## gene, computed once in val_synovium.rds, not recomputed here), direction
## concordance against `blood_dir`'s logFC for the same gene, and AUC under
## both orientation conventions (Section 2.11.6). `genes` is the REQUESTED
## panel (may include genes absent from synovium - flagged via `present`).
ct_discovery_table <- function(genes, sex_code, val, blood_dir) {
  idx_sex <- if (is.null(sex_code)) seq_along(val$sex) else which(val$sex == sex_code)
  y_sex <- droplevels(val$grp[idx_sex])
  tt <- val$tt
  rows <- lapply(genes, function(g) {
    present <- g %in% tt$gene && g %in% rownames(val$logcpm)
    if (!present) {
      return(data.frame(gene = g, present = FALSE, syn_log2FC = NA_real_, syn_adjP = NA_real_,
                         blood_log2FC = NA_real_, concordant = NA, auc_bestdir = NA_real_, auc_trainorient = NA_real_,
                         stringsAsFactors = FALSE))
    }
    row <- tt[tt$gene == g, , drop = FALSE][1, ]
    blood_fc <- unname(blood_dir$logfc[g])
    conc <- if (is.na(blood_fc)) NA else (sign(row$logFC) == sign(blood_fc))
    a_best <- ct_gene_auc(val$logcpm[g, idx_sex], y_sex)
    a_train <- if (is.na(conc) || is.na(a_best)) NA_real_ else if (isTRUE(conc)) a_best else 1 - a_best
    data.frame(gene = g, present = TRUE, syn_log2FC = row$logFC, syn_adjP = row$adj.P.Val,
               blood_log2FC = blood_fc, concordant = conc, auc_bestdir = a_best, auc_trainorient = a_train,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

## Outer cross-validation for the panel classifier: ONE fold assignment
## drives both (a) a per-fold AUC vector, for the same "CV AUC by fold" bar
## chart Diagnostic Model uses, and (b) a single pooled ROC/AUC/CI built from
## every out-of-fold prediction across all folds - this project's own
## Section 2.11.5 headline estimate (matches 23_crosstissue_biomarker_
## discovery.R's panel_auc(), generalised from logistic regression to all
## four classifiers and made fold-stratification-selectable). Re-standardises
## with the FOLD-TRAIN mean/SD only (leakage-free), same as Diagnostic
## Model's diag_cv_auc(). `pooled$available/ci_lo/n` deliberately match
## mod_diagnostic.R's `ev` shape so diag_separation_note() can be reused
## as-is on the pooled result.
ct_cv_eval <- function(Xraw, y, n_folds, refit_fn, predict_fn, stratified = TRUE, seed = 1234) {
  nf <- max(2, min(n_folds, min(table(y))))
  set.seed(seed)
  folds <- if (stratified) {
    caret::createFolds(y, k = nf)
  } else {
    fold_id <- sample(rep_len(seq_len(nf), length(y)))
    split(seq_along(y), fold_id)
  }
  pooled <- rep(NA_real_, length(y))
  fold_auc <- rep(NA_real_, length(folds))
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_along(y), te)
    if (length(unique(y[tr])) < 2) next
    mu <- colMeans(Xraw[tr, , drop = FALSE]); sg <- apply(Xraw[tr, , drop = FALSE], 2, stats::sd)
    sg[is.na(sg) | sg == 0] <- 1
    Ztr <- scale(Xraw[tr, , drop = FALSE], center = mu, scale = sg)
    Zte <- scale(Xraw[te, , drop = FALSE], center = mu, scale = sg)
    fit_i <- tryCatch(refit_fn(Ztr, y[tr]), error = function(e) NULL)
    if (is.null(fit_i)) next
    p <- tryCatch(predict_fn(fit_i, Zte), error = function(e) NULL)
    if (is.null(p)) next
    pooled[te] <- p
    if (length(unique(y[te])) == 2) {
      roc_i <- tryCatch(pROC::roc(y[te], p, quiet = TRUE, levels = levels(y), direction = "<"), error = function(e) NULL)
      fold_auc[i] <- if (is.null(roc_i)) NA_real_ else as.numeric(pROC::auc(roc_i))
    }
  }
  ok <- !is.na(pooled)
  pooled_res <- list(available = FALSE, reason = "Fewer than 4 samples ended up with an out-of-fold prediction covering both groups - try fewer folds, or simple random folds.")
  if (sum(ok) >= 4 && length(unique(y[ok])) == 2) {
    r <- tryCatch(pROC::roc(y[ok], pooled[ok], quiet = TRUE, levels = levels(y), direction = "<"), error = function(e) NULL)
    if (!is.null(r)) {
      ci <- diag_auc_ci(r)
      pooled_res <- list(available = TRUE, roc = r, auc = unname(ci["auc"]), ci_lo = unname(ci["lo"]), ci_hi = unname(ci["hi"]),
                          n = sum(ok), n_pos = sum(y[ok] == levels(y)[2]))
    } else {
      pooled_res <- list(available = FALSE, reason = "ROC could not be computed from the pooled out-of-fold predictions.")
    }
  }
  list(fold_auc = fold_auc, n_folds = length(folds), pooled = pooled_res)
}

## Fits all four classifiers on one sex's FULL synovium sample pool -
## `expr_full` is genes x samples for this sex, raw log-CPM (not yet
## z-scored). No train/test split (see module header): apparent/
## resubstitution fit + tuning both use the full sex-subset; ct_cv_eval()
## above supplies the out-of-fold estimate.
ct_fit_sex <- function(expr_full, y_full, params = list()) {
  ## caret::train(classProbs = TRUE) below requires factor levels that are
  ## valid R variable names - it make.names()s them internally to build its
  ## own predicted-probability column names, so a raw group label with a
  ## space (e.g. "multiple sclerosis") desyncs from any levels(y)-based
  ## lookup once caret has already renamed its own columns to
  ## "multiple.sclerosis". Sanitized once here, up front, matching
  ## mod_diagnostic.R::diag_fit_sex()'s identical fix; callers keep the real
  ## group names for their own display text, so nothing user-visible changes.
  levels(y_full) <- make.names(levels(y_full), unique = TRUE)
  params <- utils::modifyList(CT_DEFAULT_PARAMS, params)
  GLOBAL_SEED <- 1234
  genes <- rownames(expr_full)
  safe <- make.names(genes, unique = TRUE)
  y <- y_full

  validate(need(length(unique(y)) == 2 && all(table(y) >= 4),
                "Each group needs at least 4 samples in this sex's synovium subset to fit a panel model."))

  Zfull <- diag_zrows(expr_full); rownames(Zfull) <- safe
  Xfull <- t(Zfull); colnames(Xfull) <- safe          # z-scored, sample x gene - apparent fit + hyperparameter tuning
  Xraw <- t(expr_full); colnames(Xraw) <- safe        # raw - re-scaled per fold inside ct_cv_eval() (leakage-free)

  ## -------------------------------------------------------------------
  ## (1) Elastic net - alpha tuned by minimising CV deviance over a grid,
  ## lambda by glmnet's own inner CV, exactly as Diagnostic Model's enet fit.
  ## -------------------------------------------------------------------
  nf_a <- max(2, min(params$enet_cv_folds, min(table(y))))
  best <- NULL; bcv <- Inf
  alpha_search <- data.frame(alpha = numeric(0), cv_deviance = numeric(0))
  set.seed(GLOBAL_SEED)
  for (a in params$enet_alpha_grid) {
    cv <- tryCatch(glmnet::cv.glmnet(Xfull, y, family = "binomial", alpha = a, nfolds = nf_a, standardize = TRUE), error = function(e) NULL)
    if (!is.null(cv)) {
      alpha_search <- rbind(alpha_search, data.frame(alpha = a, cv_deviance = min(cv$cvm)))
      if (min(cv$cvm) < bcv) { bcv <- min(cv$cvm); best <- list(cv = cv, alpha = a) }
    }
  }
  validate(need(!is.null(best), "Elastic net fitting failed for every alpha in the grid - check the gene panel and sample sizes."))
  alpha_search$chosen <- alpha_search$alpha == best$alpha
  lambda_s <- if (identical(params$enet_lambda_choice, "lambda.1se")) "lambda.1se" else "lambda.min"
  enet_pred_full <- as.numeric(predict(best$cv, newx = Xfull, s = lambda_s, type = "response"))
  enet_roc_full <- pROC::roc(y, enet_pred_full, quiet = TRUE, levels = levels(y), direction = "<")
  enet_best <- ct_youden(enet_roc_full)
  enet_cv <- ct_cv_eval(
    Xraw, y, params$cv_folds,
    refit_fn = function(Ztr, ytr) glmnet::cv.glmnet(Ztr, ytr, family = "binomial", alpha = best$alpha,
                                                     nfolds = max(2, min(params$enet_cv_folds, min(table(ytr)))), standardize = TRUE),
    predict_fn = function(m, Zte) as.numeric(predict(m, newx = Zte, s = lambda_s, type = "response")),
    stratified = params$stratified_folds, seed = GLOBAL_SEED
  )
  enet <- list(model = best$cv, model_type = "enet", label = "Elastic Net", alpha = best$alpha,
               lambda_choice = lambda_s, lambda_used = best$cv[[lambda_s]], tuning_search = alpha_search,
               pred_full = enet_pred_full, roc_full = enet_roc_full, full_auc = as.numeric(pROC::auc(enet_roc_full)),
               best = enet_best, perf_full = diag_perf_at_cutoff(enet_pred_full, y, enet_best$threshold, levels(y)[2]), cv = enet_cv)

  ## -------------------------------------------------------------------
  ## (2) Random forest - mtry tuned by CV, ntree fixed, both user-overridable.
  ## -------------------------------------------------------------------
  p <- ncol(Xfull); ntree <- max(100, round(params$rf_ntree))
  mtry_search <- NULL
  if (identical(params$rf_mtry_mode, "manual") && !is.null(params$rf_mtry_manual)) {
    rf_mtry <- min(p, max(1, round(params$rf_mtry_manual)))
  } else {
    nf_rf <- max(2, min(params$rf_cv_folds, min(table(y))))
    mtry_grid <- sort(unique(pmin(p, c(1, 2, floor(sqrt(p)), floor(p / 3), floor(p / 2), p))))
    ctrl <- caret::trainControl(method = "cv", number = nf_rf, classProbs = TRUE, summaryFunction = caret::twoClassSummary)
    set.seed(GLOBAL_SEED)
    rf_tune <- tryCatch(caret::train(x = Xfull, y = y, method = "rf", metric = "ROC", trControl = ctrl,
                                      tuneGrid = expand.grid(mtry = mtry_grid), ntree = ntree), error = function(e) NULL)
    rf_mtry <- if (!is.null(rf_tune)) rf_tune$bestTune$mtry else max(1, floor(sqrt(p)))
    if (!is.null(rf_tune)) { mtry_search <- rf_tune$results[, c("mtry", "ROC")]; mtry_search$chosen <- mtry_search$mtry == rf_mtry }
  }
  set.seed(GLOBAL_SEED)
  rf_model <- randomForest::randomForest(Xfull, y, ntree = ntree, mtry = rf_mtry)
  rf_pred_full <- predict(rf_model, Xfull, type = "prob")[, levels(y)[2]]
  rf_roc_full <- pROC::roc(y, rf_pred_full, quiet = TRUE, levels = levels(y), direction = "<")
  rf_best <- ct_youden(rf_roc_full)
  rf_cv <- ct_cv_eval(
    Xraw, y, params$cv_folds,
    refit_fn = function(Ztr, ytr) randomForest::randomForest(Ztr, ytr, ntree = ntree, mtry = min(rf_mtry, ncol(Ztr))),
    predict_fn = function(m, Zte) predict(m, Zte, type = "prob")[, levels(y)[2]],
    stratified = params$stratified_folds, seed = GLOBAL_SEED
  )
  rf <- list(model = rf_model, model_type = "rf", label = "Random Forest", ntree = ntree, mtry = rf_mtry, tuning_search = mtry_search,
             pred_full = rf_pred_full, roc_full = rf_roc_full, full_auc = as.numeric(pROC::auc(rf_roc_full)),
             best = rf_best, perf_full = diag_perf_at_cutoff(rf_pred_full, y, rf_best$threshold, levels(y)[2]), cv = rf_cv)

  ## -------------------------------------------------------------------
  ## (3) SVM - cost tuned by CV over a grid, kernel user-selectable.
  ## scale = FALSE: the data is already gene-wise z-scored, so e1071's own
  ## internal re-scaling is deliberately skipped rather than double-applied.
  ## -------------------------------------------------------------------
  kernel <- params$svm_kernel %||% "linear"
  cost_search <- NULL
  if (identical(params$svm_cost_mode, "manual") && !is.null(params$svm_cost_manual)) {
    svm_cost <- params$svm_cost_manual
  } else {
    nf_svm <- max(2, min(params$svm_cv_folds, min(table(y))))
    grid <- params$svm_cost_grid
    if (!is.numeric(grid) || length(grid) == 0) grid <- CT_SVM_COST_GRID
    set.seed(GLOBAL_SEED)
    svm_tune <- tryCatch(e1071::tune(e1071::svm, train.x = Xfull, train.y = y, kernel = kernel, scale = FALSE,
                                      ranges = list(cost = grid),
                                      tunecontrol = e1071::tune.control(sampling = "cross", cross = nf_svm)), error = function(e) NULL)
    svm_cost <- if (!is.null(svm_tune)) svm_tune$best.parameters$cost else 1
    if (!is.null(svm_tune)) { cost_search <- svm_tune$performances[, c("cost", "error")]; cost_search$chosen <- cost_search$cost == svm_cost }
  }
  set.seed(GLOBAL_SEED)
  svm_model <- e1071::svm(Xfull, y, kernel = kernel, cost = svm_cost, scale = FALSE, probability = TRUE)
  svm_pred_full <- attr(predict(svm_model, Xfull, probability = TRUE), "probabilities")[, levels(y)[2]]
  svm_roc_full <- pROC::roc(y, svm_pred_full, quiet = TRUE, levels = levels(y), direction = "<")
  svm_best <- ct_youden(svm_roc_full)
  svm_cv <- ct_cv_eval(
    Xraw, y, params$cv_folds,
    refit_fn = function(Ztr, ytr) e1071::svm(Ztr, ytr, kernel = kernel, cost = svm_cost, scale = FALSE, probability = TRUE),
    predict_fn = function(m, Zte) attr(predict(m, Zte, probability = TRUE), "probabilities")[, levels(y)[2]],
    stratified = params$stratified_folds, seed = GLOBAL_SEED
  )
  svm_fit <- list(model = svm_model, model_type = "svm", label = "SVM", kernel = kernel, cost = svm_cost, tuning_search = cost_search,
                   pred_full = svm_pred_full, roc_full = svm_roc_full, full_auc = as.numeric(pROC::auc(svm_roc_full)),
                   best = svm_best, perf_full = diag_perf_at_cutoff(svm_pred_full, y, svm_best$threshold, levels(y)[2]), cv = svm_cv)

  ## -------------------------------------------------------------------
  ## (4) Logistic regression - plain, unpenalized glm(y ~ ., family =
  ## binomial) on every gene in the panel; the classifier this project's own
  ## Section 2.11.5 methodology actually fits. No hyperparameters to tune.
  ## -------------------------------------------------------------------
  lr_predict <- function(m, Znew) as.numeric(predict(m, newdata = data.frame(Znew, check.names = FALSE), type = "response"))
  lr_model <- suppressWarnings(stats::glm(y ~ ., data = data.frame(y, Xfull, check.names = FALSE), family = stats::binomial))
  lr_pred_full <- as.numeric(predict(lr_model, type = "response"))
  lr_roc_full <- pROC::roc(y, lr_pred_full, quiet = TRUE, levels = levels(y), direction = "<")
  lr_best <- ct_youden(lr_roc_full)
  lr_cv <- ct_cv_eval(
    Xraw, y, params$cv_folds,
    refit_fn = function(Ztr, ytr) suppressWarnings(stats::glm(ytr ~ ., data = data.frame(ytr, Ztr, check.names = FALSE), family = stats::binomial)),
    predict_fn = lr_predict, stratified = params$stratified_folds, seed = GLOBAL_SEED
  )
  lr <- list(model = lr_model, model_type = "lr", label = "Logistic Regression", tuning_search = NULL,
             pred_full = lr_pred_full, roc_full = lr_roc_full, full_auc = as.numeric(pROC::auc(lr_roc_full)),
             best = lr_best, perf_full = diag_perf_at_cutoff(lr_pred_full, y, lr_best$threshold, levels(y)[2]), cv = lr_cv)

  list(lr = lr, enet = enet, rf = rf, svm = svm_fit, genes = genes,
       n_samples = nrow(Xfull), n_pos = sum(y == levels(y)[2]), n_neg = sum(y == levels(y)[1]))
}

## ---------------------------------------------------------------------------
## User-uploaded validation dataset - Option B (Section 5). Builds a
## val_synovium.rds-shaped object (tt/logcpm/sex/grp; fsig/msig left empty,
## since that bundled consensus panel is specific to the project's own
## synovium script) from an uploaded raw-count expression matrix and sample
## metadata, so ct_build_sex() below can run the IDENTICAL sex-stratified
## discovery (ct_discovery_table) and panel-classifier (ct_fit_sex) code on
## either data source - no separate analysis path for uploaded data.
## ---------------------------------------------------------------------------

## Sex-adjusted limma-voom DE (filterByExpr, TMM, voom, limma, eBayes) - the
## SAME pipeline already cited in this module's own References box as the
## provenance of the bundled val_synovium.rds$tt, applied live to an
## uploaded raw-count matrix rather than invoking a new statistical method.
## `grp` must be a 2-level factor (reference = level 1, comparison = level
## 2, matching val_synovium.rds's own Normal/RA convention); `sex` a
## "F"/"M" character vector, same length as ncol(counts).
ct_voom_de_table <- function(counts, grp, sex) {
  validate(need(all(counts >= 0, na.rm = TRUE) && all(counts == round(counts), na.rm = TRUE),
                "The uploaded validation expression matrix must be raw (non-negative, whole-number) RNA-seq counts, not already-normalised or log-scale values."))
  counts <- round(as.matrix(counts)); storage.mode(counts) <- "integer"
  dge0 <- edgeR::DGEList(counts = counts)
  keepg <- edgeR::filterByExpr(dge0, group = grp)
  validate(need(sum(keepg) >= 50, "Fewer than 50 genes pass edgeR's expression filter (filterByExpr) for this group split - check the uploaded counts and group mapping."))
  dge <- edgeR::calcNormFactors(dge0[keepg, ], method = "TMM")
  design_df <- data.frame(grp = grp, sex = factor(sex))
  design <- if (length(unique(sex)) == 2) stats::model.matrix(~ grp + sex, data = design_df) else stats::model.matrix(~ grp, data = design_df)
  v <- limma::voom(dge, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  coef_name <- paste0("grp", levels(grp)[2])
  tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  list(tt = tt, logcpm = v$E)
}

## Assembles the val_synovium.rds-shaped object for an uploaded cohort.
## `meta` is the raw uploaded metadata data.frame; `id_col`/`sex_col`/
## `group_col` are the user's chosen column mappings (same "map your own
## columns" pattern as Diagnostic Model's External Validation tab);
## `ref_group`/`comp_group` are the two group values to keep (comp = the
## positive/disease class, matching val_synovium.rds's own RA-as-level-2
## convention).
ct_build_uploaded_val <- function(expr, meta, id_col, sex_col, group_col, ref_group, comp_group) {
  sample_id <- as.character(meta[[id_col]])
  common <- intersect(colnames(expr), sample_id)
  validate(need(length(common) >= 12, "Fewer than 12 sample IDs in the uploaded validation expression matrix match the metadata sample-ID column. Check the column mapping."))
  expr <- expr[, common, drop = FALSE]
  meta <- meta[match(common, sample_id), , drop = FALSE]

  sex_raw <- as.character(meta[[sex_col]])
  sex <- ifelse(grepl("^f", sex_raw, ignore.case = TRUE), "F",
                ifelse(grepl("^m", sex_raw, ignore.case = TRUE), "M", NA_character_))
  validate(need(!anyNA(sex), "The sex column must contain values starting with \"F\"/\"f\" (female) or \"M\"/\"m\" (male) for every matched sample."))

  grp_raw <- as.character(meta[[group_col]])
  keep <- grp_raw %in% c(ref_group, comp_group)
  validate(need(sum(keep) >= 12, "Fewer than 12 samples match the chosen reference/comparison groups."))
  expr <- expr[, keep, drop = FALSE]; sex <- sex[keep]
  grp <- factor(grp_raw[keep], levels = c(ref_group, comp_group))
  validate(need(all(table(sex) >= 4), "Each sex needs at least 4 samples (after the reference/comparison group filter) in the uploaded validation cohort."))

  de <- ct_voom_de_table(expr, grp, sex)
  list(logcpm = de$logcpm, grp = grp, sex = sex, tt = de$tt, fsig = character(0), msig = character(0))
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------

mod_crosstissue_params_box <- function(ns, prefix, method_label, defaults_desc, body = NULL) {
  box(
    width = 12, title = sprintf("%s parameters", method_label), status = "primary", solidHeader = FALSE,
    p(class = "submodule-desc", defaults_desc),
    body
  )
}

## One box, reused for either the "Full Fit" (apparent) or "Cross-Validated"
## (pooled out-of-fold) view of one model x sex - full width, matching
## Diagnostic Model's "utilise the entire page" layout.
mod_crosstissue_fullfit_panel <- function(ns, prefix, roc_height = 250) {
  box(
    width = NULL, status = "primary", solidHeader = FALSE,
    withSpinner(uiOutput(ns(paste0(prefix, "_full_stats"))), color = "#2c6fbb", type = 6),
    fluidRow(
      column(4, h5("ROC (apparent / resubstitution fit)"), withSpinner(plotOutput(ns(paste0(prefix, "_full_roc_plot")), height = roc_height), color = "#2c6fbb", type = 6)),
      column(4, h5("Cross-validated AUC by fold"), withSpinner(plotOutput(ns(paste0(prefix, "_cv_fold_plot")), height = roc_height), color = "#2c6fbb", type = 6)),
      column(4, h5("Hyperparameter tuning - explore the grid"), withSpinner(plotOutput(ns(paste0(prefix, "_tuning_plot")), height = roc_height), color = "#2c6fbb", type = 6))
    ),
    div(class = "table-toolbar",
        downloadButton(ns(paste0(prefix, "_full_download")), "Performance (CSV)", class = "btn-sm"),
        downloadButton(ns(paste0(prefix, "_model_download")), "Model (.rds)", class = "btn-sm")),
    DT::dataTableOutput(ns(paste0(prefix, "_full_table")))
  )
}

mod_crosstissue_cv_panel <- function(ns, prefix, roc_height = 300) {
  box(
    width = NULL, status = "primary", solidHeader = FALSE,
    withSpinner(uiOutput(ns(paste0(prefix, "_cv_summary"))), color = "#2c6fbb", type = 6),
    withSpinner(plotOutput(ns(paste0(prefix, "_cv_roc_plot")), height = roc_height), color = "#2c6fbb", type = 6),
    div(class = "table-toolbar", downloadButton(ns(paste0(prefix, "_cv_download")), "Performance (CSV)", class = "btn-sm")),
    DT::dataTableOutput(ns(paste0(prefix, "_cv_table")))
  )
}

mod_crosstissue_discovery_sex_panel <- function(ns, sex_label) {
  run_id <- paste0("run_", sex_label, "_btn_disc")
  cond <- sprintf("input['%s'] > 0 || input['%s'] > 0 || input['%s'] > 0",
                   ns(run_id), ns(paste0("run_", sex_label, "_btn_full")), ns(paste0("run_", sex_label, "_btn_cv")))
  sex_title <- tools::toTitleCase(sex_label)
  tagList(
    actionButton(ns(run_id), paste("Run", sex_title), icon = icon("play"), class = "btn-primary btn-sm"),
    div(style = "height:10px;"),
    conditionalPanel(
      condition = cond,
      withSpinner(uiOutput(ns(paste0(sex_label, "_disc_stats"))), color = "#2c6fbb", type = 6),
      fluidRow(
        column(6, withSpinner(plotOutput(ns(paste0(sex_label, "_concordance_plot")), height = 380), color = "#2c6fbb", type = 6)),
        column(6, withSpinner(plotOutput(ns(paste0(sex_label, "_geneauc_plot")), height = 380), color = "#2c6fbb", type = 6))
      ),
      div(class = "table-toolbar", downloadButton(ns(paste0(sex_label, "_disc_download")), "Per-gene table (CSV)", class = "btn-sm")),
      DT::dataTableOutput(ns(paste0(sex_label, "_disc_table")))
    )
  )
}

mod_crosstissue_fullfit_sex_panel <- function(ns, sex_label) {
  run_id <- paste0("run_", sex_label, "_btn_full")
  cond <- sprintf("input['%s'] > 0 || input['%s'] > 0 || input['%s'] > 0",
                   ns(run_id), ns(paste0("run_", sex_label, "_btn_disc")), ns(paste0("run_", sex_label, "_btn_cv")))
  sex_title <- tools::toTitleCase(sex_label)
  tagList(
    actionButton(ns(run_id), paste("Run", sex_title), icon = icon("play"), class = "btn-primary btn-sm"),
    div(style = "height:10px;"),
    conditionalPanel(
      condition = cond,
      tabsetPanel(
        id = ns(paste0(sex_label, "_full_model_pills")), type = "pills",
        tabPanel("Logistic Regression", br(), mod_crosstissue_fullfit_panel(ns, paste0(sex_label, "_lr"))),
        tabPanel("Elastic Net", br(), mod_crosstissue_fullfit_panel(ns, paste0(sex_label, "_enet"))),
        tabPanel("Random Forest", br(), mod_crosstissue_fullfit_panel(ns, paste0(sex_label, "_rf"))),
        tabPanel("SVM", br(), mod_crosstissue_fullfit_panel(ns, paste0(sex_label, "_svm")))
      ),
      box(width = NULL, title = sprintf("Full-fit comparison - %s", sex_title), status = "primary", solidHeader = FALSE,
          DT::dataTableOutput(ns(paste0(sex_label, "_full_compare_table")))),
      uiOutput(ns(paste0(sex_label, "_full_result_line")))
    )
  )
}

mod_crosstissue_cv_sex_panel <- function(ns, sex_label) {
  run_id <- paste0("run_", sex_label, "_btn_cv")
  cond <- sprintf("input['%s'] > 0 || input['%s'] > 0 || input['%s'] > 0",
                   ns(run_id), ns(paste0("run_", sex_label, "_btn_disc")), ns(paste0("run_", sex_label, "_btn_full")))
  sex_title <- tools::toTitleCase(sex_label)
  tagList(
    actionButton(ns(run_id), paste("Run", sex_title), icon = icon("play"), class = "btn-primary btn-sm"),
    div(style = "height:10px;"),
    conditionalPanel(
      condition = cond,
      tabsetPanel(
        id = ns(paste0(sex_label, "_cv_model_pills")), type = "pills",
        tabPanel("Logistic Regression", br(), mod_crosstissue_cv_panel(ns, paste0(sex_label, "_lr"))),
        tabPanel("Elastic Net", br(), mod_crosstissue_cv_panel(ns, paste0(sex_label, "_enet"))),
        tabPanel("Random Forest", br(), mod_crosstissue_cv_panel(ns, paste0(sex_label, "_rf"))),
        tabPanel("SVM", br(), mod_crosstissue_cv_panel(ns, paste0(sex_label, "_svm")))
      ),
      box(width = NULL, title = sprintf("Cross-validated comparison - %s", sex_title), status = "primary", solidHeader = FALSE,
          DT::dataTableOutput(ns(paste0(sex_label, "_cv_compare_table"))))
    )
  )
}

mod_crosstissue_crossdata_sex_panel <- function(ns, sex_label) {
  sex_title <- tools::toTitleCase(sex_label)
  tagList(
    h4(sex_title),
    withSpinner(uiOutput(ns(paste0(sex_label, "_crossdata_note"))), color = "#2c6fbb", type = 6),
    withSpinner(plotOutput(ns(paste0(sex_label, "_crossdata_plot")), height = 320), color = "#2c6fbb", type = 6),
    DT::dataTableOutput(ns(paste0(sex_label, "_crossdata_table"))),
    div(style = "height:18px;")
  )
}

mod_crosstissue_ui <- function(id) {
  ns <- NS(id)
  tagList(
    ## Scoped to this module only (.ct-module), same defensive floor/wrap as
    ## Diagnostic Model's own .diag-module rule, so a longer valueBox() value
    ## (e.g. an AUC ± CI) never overflows its tile.
    tags$style(HTML(".ct-module .small-box h3 { font-size: 20px; white-space: normal; overflow-wrap: break-word; line-height: 1.2; }")),
    div(class = "ct-module",
    fluidRow(
      column(
        3,
        box(
          width = NULL, title = "Validation dataset", status = "primary", solidHeader = FALSE,
          radioButtons(
            ns("val_source"), NULL,
            choiceNames = list(
              tagList(icon("database"), " Use preloaded validation dataset (Synovium, GSE89408)"),
              tagList(icon("file-arrow-up"), " Upload my own validation dataset")
            ),
            choiceValues = list("preloaded", "upload"), selected = "preloaded"
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'upload'", ns("val_source")),
            p(class = "submodule-desc", "Provide a raw RNA-seq count matrix and sample metadata for an independent validation-tissue cohort. The same sex-stratified discovery and panel-classifier workflow below then runs on this cohort instead of the bundled synovium dataset."),
            fileInput(ns("val_expr_file"), "Validation expression matrix (raw counts)", accept = c(".csv", ".rds", ".Rds")),
            div(class = "empty-note", style = "font-size: 12.5px; margin-top: -8px;", icon("circle-info"),
                "CSV or RDS. Genes in rows, samples in columns; for CSV, the first column is the gene ID."),
            fileInput(ns("val_meta_file"), "Validation sample metadata", accept = c(".csv", ".rds", ".Rds")),
            uiOutput(ns("val_column_mapping"))
          )
        ),
        box(
          width = NULL, title = "Gene panel & synovium contrast", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "RA vs Normal synovium (GSE89408, or your uploaded validation cohort), sex-stratified."),
          radioButtons(
            ns("panel_source"), NULL,
            choiceNames = list(
              tagList(icon("diagram-project"), " Follow this project's pipeline (recommended)"),
              tagList(icon("file-arrow-up"), " Paste my own gene list")
            ),
            choiceValues = list("project", "own"), selected = "project"
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'project'", ns("panel_source")),
            uiOutput(ns("project_source_ui"))
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'own'", ns("panel_source")),
            p(class = "submodule-desc", "Same list for both sexes."),
            textAreaInput(ns("gene_list"), NULL, rows = 5, placeholder = "TNF\nIL6\nSTAT3\n...")
          ),
          tags$hr(),
          uiOutput(ns("blood_direction_ui")),
          div(style = "margin-top:10px;", uiOutput(ns("saved_runs_ui")))
        ),
        box(
          width = NULL, title = "Advanced filters", status = "primary", solidHeader = FALSE,
          numericInput(ns("sig_cutoff"), "Significance threshold (BH-adjusted P)", value = 0.05, min = 0.001, max = 0.5, step = 0.005),
          radioButtons(
            ns("orient_view"), "Gene AUC orientation",
            choiceNames = list(
              tagList("Best-direction ", tags$small("(discovery)")),
              tagList("Train-fixed ", tags$small("(cross-dataset; AUC<0.5 = direction reversed)"))
            ),
            choiceValues = list("bestdir", "trainorient"), selected = "bestdir"
          ),
          tags$hr(),
          strong("Panel classifier - outer cross-validation"),
          numericInput(ns("cv_folds"), "Outer folds", value = 10, min = 3, max = 20, step = 1),
          radioButtons(
            ns("stratified_folds"), NULL,
            choiceNames = list("Stratified by disease status (recommended)",
                                "Simple random (matches this project's own script)"),
            choiceValues = list("stratified", "random"), selected = "stratified"
          ),
          div(class = "empty-note", style = "font-size: 12.5px;", icon("circle-info"),
              "Shared by all four models; each model's own tuning folds are set once you Run.")
        )
      ),
      column(
        9,
        tabsetPanel(
          id = ns("main_tabs"), type = "tabs",
          tabPanel(
            "Synovium Discovery & Concordance", br(),
            p(class = "submodule-desc", "Pick Pooled (all), Female or Male, then Run."),
            tabsetPanel(
              id = ns("disc_sex_tabs"), type = "tabs",
              tabPanel("Pooled (All)", br(), mod_crosstissue_discovery_sex_panel(ns, "pooled")),
              tabPanel("Female", br(), mod_crosstissue_discovery_sex_panel(ns, "female")),
              tabPanel("Male", br(), mod_crosstissue_discovery_sex_panel(ns, "male"))
            )
          ),
          tabPanel(
            "Panel Classifier - Full Fit", br(),
            p(class = "submodule-desc", "Apparent fit on the full synovium sample - an optimistic upper bound, not held-out performance."),
            fluidRow(
              column(
                9,
                tabsetPanel(
                  id = ns("full_sex_tabs"), type = "tabs",
                  tabPanel("Pooled (All)", br(), mod_crosstissue_fullfit_sex_panel(ns, "pooled")),
                  tabPanel("Female", br(), mod_crosstissue_fullfit_sex_panel(ns, "female")),
                  tabPanel("Male", br(), mod_crosstissue_fullfit_sex_panel(ns, "male"))
                )
              ),
              ## Small, right-hand box at the same level as the Female/Male
              ## tabs above - whichever model pill was clicked most recently
              ## in either sex's Full Fit or Cross-Validated tab (see
              ## active_model_pill below), rather than a full-width box
              ## stacked above the tabs.
              column(3, uiOutput(ns("model_params_ui")))
            )
          ),
          tabPanel(
            "Panel Classifier - Cross-Validated", br(),
            p(class = "submodule-desc", "Out-of-fold performance - this project's headline synovium estimate."),
            tabsetPanel(
              id = ns("cv_sex_tabs"), type = "tabs",
              tabPanel("Pooled (All)", br(), mod_crosstissue_cv_sex_panel(ns, "pooled")),
              tabPanel("Female", br(), mod_crosstissue_cv_sex_panel(ns, "female")),
              tabPanel("Male", br(), mod_crosstissue_cv_sex_panel(ns, "male"))
            )
          ),
          tabPanel(
            "Cross-Dataset Comparison", br(),
            p(class = "submodule-desc", "Synovium AUC alongside Diagnostic Model's saved blood AUC - not a transfer of the blood model."),
            mod_crosstissue_crossdata_sex_panel(ns, "pooled"),
            tags$hr(),
            mod_crosstissue_crossdata_sex_panel(ns, "female"),
            tags$hr(),
            mod_crosstissue_crossdata_sex_panel(ns, "male")
          )
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

mod_crosstissue_server <- function(id, dataset, results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Independent, read-only validation cohort - loaded once, not from the
    ## shared `dataset` reactiveValues (unlike every blood-facing submodule),
    ## since this module's whole point is a compartment `dataset` never
    ## holds. `tt` converted to a plain data.frame once here so every
    ## gene == g row lookup below doesn't repeat data.table dispatch.
    val_bundled <- { v <- readRDS(VAL_SYNOVIUM_RDS); v$tt <- as.data.frame(v$tt); v }
    bundled_dge <- tryCatch(readRDS(DGE_RESULTS_RDS), error = function(e) NULL)

    ## -----------------------------------------------------------------
    ## Validation dataset source - Option A (bundled val_synovium.rds,
    ## read once above) vs Option B (user upload, Section 5/6). Both sides
    ## of `val_active()` return the same shape (tt/logcpm/sex/grp/fsig/
    ## msig), so every downstream function below (ct_project_panel_genes,
    ## ct_build_sex, ct_discovery_table, ct_fit_sex) is unchanged by which
    ## source is active.
    ## -----------------------------------------------------------------

    val_meta_raw <- reactive({
      req(input$val_meta_file)
      path <- input$val_meta_file$datapath
      if (grepl("\\.rds$", input$val_meta_file$name, ignore.case = TRUE)) {
        d <- readRDS(path)
        validate(need(is.data.frame(d), "The uploaded validation metadata RDS file must contain a data frame."))
        as.data.frame(d)
      } else {
        as.data.frame(data.table::fread(path, showProgress = FALSE))
      }
    })

    val_expr_raw <- reactive({
      req(input$val_expr_file)
      if (grepl("\\.rds$", input$val_expr_file$name, ignore.case = TRUE)) {
        res <- tx_parse_expr_matrix_rds(input$val_expr_file$datapath)
        validate(need(res$ok, res$error))
        res$mat
      } else {
        m <- as.data.frame(data.table::fread(input$val_expr_file$datapath, showProgress = FALSE))
        rn <- as.character(m[[1]]); m <- as.matrix(m[, -1, drop = FALSE]); rownames(m) <- rn
        m
      }
    })

    output$val_column_mapping <- renderUI({
      req(input$val_expr_file, val_meta_raw())
      cols <- colnames(val_meta_raw())
      tagList(
        fluidRow(
          column(4, selectInput(ns("val_map_id"), "Sample ID column", choices = cols, selected = cols[1], selectize = FALSE)),
          column(4, selectInput(ns("val_map_sex"), "Sex column", choices = cols, selectize = FALSE)),
          column(4, selectInput(ns("val_map_group"), "Group column", choices = cols, selectize = FALSE))
        ),
        uiOutput(ns("val_group_pick_ui"))
      )
    })

    output$val_group_pick_ui <- renderUI({
      req(input$val_map_group)
      groups <- sort(unique(stats::na.omit(as.character(val_meta_raw()[[input$val_map_group]]))))
      validate(need(length(groups) >= 2, "The chosen group column needs at least two distinct values."))
      fluidRow(
        column(6, selectInput(ns("val_ref_group"), "Reference group (e.g. healthy / control)", choices = groups, selected = groups[1], selectize = FALSE)),
        column(6, selectInput(ns("val_comp_group"), "Comparison group (e.g. disease)", choices = groups, selected = groups[min(2, length(groups))], selectize = FALSE))
      )
    })

    val_uploaded <- reactive({
      req(input$val_map_id, input$val_map_sex, input$val_map_group, input$val_ref_group, input$val_comp_group)
      validate(need(input$val_ref_group != input$val_comp_group, "Reference and comparison group must be different."))
      ct_build_uploaded_val(val_expr_raw(), val_meta_raw(), input$val_map_id, input$val_map_sex, input$val_map_group,
                            input$val_ref_group, input$val_comp_group)
    })

    val_active <- reactive({
      if (identical(input$val_source %||% "preloaded", "upload")) val_uploaded() else val_bundled
    })

    ## -----------------------------------------------------------------
    ## Gene panel sources - same two-radio-button pattern as Diagnostic
    ## Model, but the "project pipeline" fallback is va$fsig/va$msig (this
    ## project's own consensus panels, already the ones the bundled synovium
    ## validation itself used - empty for an uploaded cohort, which has no
    ## such bundled panel), not Diagnostic Model's blood FS_input CSVs.
    ## -----------------------------------------------------------------

    ct_project_panel_genes <- function(sex_label) {
      live <- results$featureselection[[sex_label]]$consensus_genes
      if (!is.null(live) && length(live) >= 2) {
        return(list(genes = live, is_live = TRUE,
                    note = sprintf("%d genes from this session's live Feature Selection %s consensus panel.", length(live), sex_label)))
      }
      va <- val_active()
      bundled <- if (isTRUE(dataset$is_bundled_reference)) switch(sex_label, female = va$fsig, male = va$msig, NULL) else NULL
      bundled <- unique(as.character(bundled))
      if (length(bundled) >= 2) {
        return(list(genes = bundled, is_live = FALSE,
                    note = sprintf("%d genes from this project's own bundled %s consensus panel (the one its own synovium validation used).", length(bundled), sex_label)))
      }
      list(genes = character(0), is_live = FALSE,
           note = sprintf("No %s candidate genes available%s.", sex_label,
                           if (identical(input$val_source %||% "preloaded", "upload")) " - run Feature Selection or paste a gene list; an uploaded validation cohort has no bundled consensus panel" else ""))
    }

    ct_own_panel_genes <- function(sex_label) {
      genes <- unique(trimws(unlist(strsplit(input$gene_list %||% "", "[,\n\t ]+"))))
      genes <- genes[nzchar(genes)]
      list(genes = genes, is_live = FALSE, note = sprintf("%d pasted genes.", length(genes)))
    }

    output$project_source_ui <- renderUI({
      f_live <- results$featureselection$female$consensus_genes
      m_live <- results$featureselection$male$consensus_genes
      has_live <- length(f_live) >= 2 && length(m_live) >= 2
      if (has_live) {
        div(class = "empty-note", icon("check"), sprintf("Live panel: %d F / %d M genes.", length(f_live), length(m_live)))
      } else if (isTRUE(dataset$is_bundled_reference) && !identical(input$val_source %||% "preloaded", "upload")) {
        div(class = "empty-note", icon("circle-info"), "Bundled consensus panel (no live Feature Selection yet).")
      } else {
        div(class = "empty-note", icon("triangle-exclamation"), "No consensus panel available for the active dataset - run Feature Selection first, or paste a gene list instead.")
      }
    })

    ## -----------------------------------------------------------------
    ## Blood direction of effect (for concordance) - live if this session
    ## has already run Differential Expression for that sex
    ## (results$dge_runs, same guess_run() pattern mod_candidates.R uses to
    ## match a contrast label to a sex), else this project's own bundled
    ## dge_results.rds. Full per-gene logFC (not just sign), so the
    ## concordance scatter plot can show magnitude, not just direction.
    ## -----------------------------------------------------------------

    ct_blood_direction <- function(sex_label) {
      runs <- results$dge_runs %||% list()
      sex_word_pattern <- "\\b(female|male)\\b|\\bF\\b|\\bM\\b"
      pattern <- switch(sex_label, female = "\\bfemale\\b|\\bF\\b", male = "\\bmale\\b|\\bM\\b", NULL)
      if (length(runs) > 0) {
        labels <- vapply(runs, function(r) r$contrast, character(1))
        matches <- if (is.null(pattern)) !grepl(sex_word_pattern, labels, ignore.case = TRUE, perl = TRUE)
                   else grepl(pattern, labels, ignore.case = TRUE, perl = TRUE)
        hit_ids <- names(runs)[matches]
        if (length(hit_ids) > 0) {
          r <- runs[[utils::tail(hit_ids, 1)]]
          return(list(logfc = stats::setNames(r$table$logFC, r$table$gene),
                      note = "live DE run", is_live = TRUE))
        }
      }
      if (isTRUE(dataset$is_bundled_reference) && !is.null(bundled_dge)) {
        key <- switch(sex_label, female = "Female", male = "Male", pooled = "All", NULL)
        d <- if (!is.null(key)) bundled_dge$res[[key]] else NULL
        if (!is.null(d) && nrow(d) > 0) {
          return(list(logfc = stats::setNames(d$logFC, d$gene),
                      note = "bundled DE", is_live = FALSE))
        }
      }
      list(logfc = stats::setNames(numeric(0), character(0)), note = "unavailable", is_live = FALSE)
    }

    output$blood_direction_ui <- renderUI({
      p <- ct_blood_direction("pooled"); f <- ct_blood_direction("female"); m <- ct_blood_direction("male")
      div(class = "empty-note", title = "Reference used to check each gene's direction of effect against blood.",
          icon("arrows-turn-to-dots"), sprintf("Blood direction: Pooled = %s, F = %s, M = %s.", p$note, f$note, m$note))
    })

    ## -----------------------------------------------------------------
    ## Advanced parameters - live-editable, one box per model type (mirrors
    ## Diagnostic Model / Feature Selection), plus the shared outer-CV
    ## controls in the left "Advanced filters" box (fold count/stratification
    ## apply to all four models at once - see the module header).
    ## -----------------------------------------------------------------

    ct_advanced_params <- function() {
      alpha_grid <- suppressWarnings(as.numeric(trimws(strsplit(input$enet_alpha_grid %||% "", ",")[[1]])))
      alpha_grid <- alpha_grid[!is.na(alpha_grid) & alpha_grid >= 0 & alpha_grid <= 1]
      cost_grid <- suppressWarnings(as.numeric(trimws(strsplit(input$svm_cost_grid %||% "", ",")[[1]])))
      cost_grid <- cost_grid[!is.na(cost_grid) & cost_grid > 0]
      list(
        cv_folds = input$cv_folds %||% CT_DEFAULT_PARAMS$cv_folds,
        stratified_folds = !identical(input$stratified_folds, "random"),
        enet_cv_folds = input$enet_cv_folds %||% CT_DEFAULT_PARAMS$enet_cv_folds,
        enet_alpha_grid = if (length(alpha_grid) > 0) alpha_grid else CT_DEFAULT_PARAMS$enet_alpha_grid,
        enet_lambda_choice = input$enet_lambda_choice %||% CT_DEFAULT_PARAMS$enet_lambda_choice,
        rf_cv_folds = input$rf_cv_folds %||% CT_DEFAULT_PARAMS$rf_cv_folds,
        rf_ntree = input$rf_ntree %||% CT_DEFAULT_PARAMS$rf_ntree,
        rf_mtry_mode = input$rf_mtry_mode %||% CT_DEFAULT_PARAMS$rf_mtry_mode,
        rf_mtry_manual = input$rf_mtry_manual,
        svm_cv_folds = input$svm_cv_folds %||% CT_DEFAULT_PARAMS$svm_cv_folds,
        svm_kernel = input$svm_kernel %||% CT_DEFAULT_PARAMS$svm_kernel,
        svm_cost_mode = input$svm_cost_mode %||% CT_DEFAULT_PARAMS$svm_cost_mode,
        svm_cost_manual = input$svm_cost_manual %||% CT_DEFAULT_PARAMS$svm_cost_manual,
        svm_cost_grid = if (length(cost_grid) > 0) cost_grid else CT_SVM_COST_GRID
      )
    }

    ## -----------------------------------------------------------------
    ## Run: one Run per sex computes discovery + all four classifiers at
    ## once (shared trigger across the three Run buttons that sex has, one
    ## per tab - same fan-in pattern as Diagnostic Model's two buttons).
    ## -----------------------------------------------------------------

    ct_build_sex <- function(sex_label) {
      va <- val_active()
      is_upload <- identical(input$val_source %||% "preloaded", "upload")
      sex_code <- switch(sex_label, female = "F", male = "M", NULL)
      idx_sex <- if (is.null(sex_code)) seq_along(va$sex) else which(va$sex == sex_code)
      validate(need(length(idx_sex) > 0, sprintf("No %s samples in the validation dataset.", sex_label)))
      y_full <- droplevels(va$grp[idx_sex])
      validate(need(length(unique(y_full)) == 2 && all(table(y_full) >= 4),
                    sprintf("The %s validation subset needs at least 4 samples in each group (%s).", sex_label, paste(levels(va$grp), collapse = " / "))))

      cand <- if (identical(input$panel_source, "project")) ct_project_panel_genes(sex_label) else ct_own_panel_genes(sex_label)
      genes_req <- cand$genes
      validate(need(length(genes_req) >= 2, sprintf("Fewer than 2 %s genes from the chosen panel.", sex_label)))
      genes_present <- intersect(genes_req, rownames(va$logcpm))
      validate(need(length(genes_present) >= 2, sprintf("Fewer than 2 %s panel genes are present in the validation dataset.", sex_label)))

      bd <- ct_blood_direction(sex_label)
      discovery <- ct_discovery_table(genes_req, sex_code, va, bd)

      expr_sub <- va$logcpm[genes_present, idx_sex, drop = FALSE]
      fit <- ct_fit_sex(expr_sub, y_full, params = ct_advanced_params())
      fit$discovery <- discovery
      fit$candidate_note <- cand$note
      fit$blood_note <- bd$note
      fit$n_input <- length(genes_req); fit$n_present <- length(genes_present)
      fit$sex_label <- sex_label
      fit$grp_levels <- levels(va$grp)
      fit$dataset_label <- if (is_upload) "user-uploaded validation cohort" else "synovium, GSE89408"
      fit
    }

    pooled_run_trigger <- reactiveVal(0)
    female_run_trigger <- reactiveVal(0)
    male_run_trigger <- reactiveVal(0)
    lapply(c("run_pooled_btn_disc", "run_pooled_btn_full", "run_pooled_btn_cv"), function(bid) {
      observeEvent(input[[bid]], { pooled_run_trigger(isolate(pooled_run_trigger()) + 1) }, ignoreInit = TRUE)
    })
    lapply(c("run_female_btn_disc", "run_female_btn_full", "run_female_btn_cv"), function(bid) {
      observeEvent(input[[bid]], { female_run_trigger(isolate(female_run_trigger()) + 1) }, ignoreInit = TRUE)
    })
    lapply(c("run_male_btn_disc", "run_male_btn_full", "run_male_btn_cv"), function(bid) {
      observeEvent(input[[bid]], { male_run_trigger(isolate(male_run_trigger()) + 1) }, ignoreInit = TRUE)
    })

    ct_result_pooled <- eventReactive(pooled_run_trigger(), ct_build_sex("pooled"), ignoreInit = TRUE)
    ct_result_female <- eventReactive(female_run_trigger(), ct_build_sex("female"), ignoreInit = TRUE)
    ct_result_male <- eventReactive(male_run_trigger(), ct_build_sex("male"), ignoreInit = TRUE)

    ct_has_run <- reactiveVal(FALSE)
    observeEvent(pooled_run_trigger(), ct_has_run(TRUE), ignoreInit = TRUE)
    observeEvent(female_run_trigger(), ct_has_run(TRUE), ignoreInit = TRUE)
    observeEvent(male_run_trigger(), ct_has_run(TRUE), ignoreInit = TRUE)

    ## Which model's params box to show: whichever pill was clicked most
    ## recently, in any stratum's Full Fit or Cross-Validated tab.
    active_model_pill <- reactiveVal("Logistic Regression")
    lapply(c("pooled_full_model_pills", "female_full_model_pills", "male_full_model_pills",
             "pooled_cv_model_pills", "female_cv_model_pills", "male_cv_model_pills"), function(iid) {
      observeEvent(input[[iid]], { active_model_pill(input[[iid]]) }, ignoreInit = TRUE)
    })

    lr_params_box <- function() {
      mod_crosstissue_params_box(
        ns, "lr", "Logistic Regression",
        "Plain logistic regression over all panel genes, with no shrinkage (this project's own methodology, Section 2.11.5). Nothing to tune here - fold count and splitting are set in \"Advanced filters\" on the left."
      )
    }

    ## Bodies are stacked (no fluidRow/column split) rather than side-by-side
    ## fields: this box now sits in the narrow right-hand column next to the
    ## Female/Male tabs (see mod_crosstissue_ui's Full Fit tabPanel), not a
    ## full-width box, so side-by-side columns would just cramp into slivers.
    enet_params_box <- function() {
      mod_crosstissue_params_box(
        ns, "enet", "Elastic Net",
        "Elastic Net shrinks less-useful genes toward zero, which helps when panel genes are correlated. Alpha is auto-tuned over the grid below; lambda.min is used by default.",
        tagList(
          numericInput(ns("enet_cv_folds"), "Inner tuning folds (alpha search)", value = CT_DEFAULT_PARAMS$enet_cv_folds, min = 3, max = 10, step = 1),
          textInput(ns("enet_alpha_grid"), "Alpha grid (0 = ridge … 1 = LASSO, comma-separated)", value = paste(CT_DEFAULT_PARAMS$enet_alpha_grid, collapse = ", ")),
          radioButtons(ns("enet_lambda_choice"), "Lambda", choices = c("lambda.min (default)" = "lambda.min", "lambda.1se (sparser)" = "lambda.1se"), selected = "lambda.min")
        )
      )
    }

    rf_params_box <- function() {
      mod_crosstissue_params_box(
        ns, "rf", "Random Forest",
        "Random Forest builds many decision trees, each looking at a random subset of genes, and averages their votes. How many genes each tree considers (\"mtry\") is picked automatically from the grid below; the number of trees is fixed.",
        tagList(
          numericInput(ns("rf_cv_folds"), "Inner tuning folds (mtry search)", value = CT_DEFAULT_PARAMS$rf_cv_folds, min = 3, max = 10, step = 1),
          numericInput(ns("rf_ntree"), "Number of trees", value = CT_DEFAULT_PARAMS$rf_ntree, min = 100, max = 5000, step = 100),
          radioButtons(ns("rf_mtry_mode"), "mtry (per split)", choices = c("Auto-tune (default)" = "auto", "Manual" = "manual"), selected = "auto"),
          conditionalPanel(condition = sprintf("input['%s'] == 'manual'", ns("rf_mtry_mode")),
                            numericInput(ns("rf_mtry_manual"), "mtry value", value = 5, min = 1, max = 500, step = 1))
        )
      )
    }

    svm_params_box <- function() {
      mod_crosstissue_params_box(
        ns, "svm", "SVM",
        "SVM finds the boundary that best separates RA from Normal samples. \"Cost\" (how strictly it fits the training data) is auto-tuned from the grid below; the boundary is linear by default.",
        tagList(
          numericInput(ns("svm_cv_folds"), "Inner tuning folds (cost search)", value = CT_DEFAULT_PARAMS$svm_cv_folds, min = 3, max = 10, step = 1),
          radioButtons(ns("svm_kernel"), "Kernel", choices = c("Linear (default)" = "linear", "Radial" = "radial"), selected = "linear"),
          radioButtons(ns("svm_cost_mode"), "Cost (C)", choices = c("Auto-tune via CV grid (default)" = "auto", "Manual" = "manual"), selected = "auto"),
          conditionalPanel(condition = sprintf("input['%s'] == 'auto'", ns("svm_cost_mode")),
                            textInput(ns("svm_cost_grid"), "Cost grid (comma-separated)", value = paste(CT_SVM_COST_GRID, collapse = ", "))),
          conditionalPanel(condition = sprintf("input['%s'] == 'manual'", ns("svm_cost_mode")),
                            numericInput(ns("svm_cost_manual"), "Cost value", value = 1, min = 0.001, step = 0.1))
        )
      )
    }

    output$model_params_ui <- renderUI({
      req(ct_has_run())
      switch(active_model_pill(),
        "Elastic Net" = enet_params_box(),
        "Random Forest" = rf_params_box(),
        "SVM" = svm_params_box(),
        lr_params_box()
      )
    })

    ## A plain shinydashboard::box(collapsible = TRUE) never actually opens
    ## here: AdminLTE's box-widget click handler only binds to boxes present
    ## at initial page load, not ones inserted later via renderUI/uiOutput
    ## (this box's entire reason for existing) - so its "+" toggle looks like
    ## a working collapse but silently does nothing. A native <details>
    ## disclosure needs no JS at all and always works, dynamically inserted
    ## or not; box/box-header/box-title/box-body classes are reused purely
    ## for the visual match to every other box on this page, not their JS.
    output$references_box_ui <- renderUI({
      req(ct_has_run())
      tags$details(
        class = "box box-primary",
        tags$summary(class = "box-header", style = "cursor: pointer;", tags$h3(class = "box-title", "References")),
        div(class = "box-body",
        tags$ul(
          class = "dge-ref-list",
          tags$li(strong("Scope of a gene-set, non-transferred validation: "), "Justice AC, Covinsky KE, Berlin JA (1999). Assessing the Generalizability of Prognostic Information. ", tags$em("Ann Intern Med"), ", 130(6), 515-524; Steyerberg EW, Harrell FE (2016). Prediction models need appropriate internal, internal-external, and external validation. ", tags$em("J Clin Epidemiol"), ", 69, 245-247."),
          tags$li(strong("Apparent AUC as an optimistic upper bound: "), "Harrell FE, Lee KL, Mark DB (1996). Multivariable prognostic models. ", tags$em("Stat Med"), ", 15(4), 361-387."),
          tags$li(strong("Synovium DE (filterByExpr, TMM, voom, limma, eBayes): "), "Chen Y, Lun ATL, Smyth GK (2016). ", tags$em("F1000Research"), ", 5, 1438; Robinson MD, Oshlack A (2010). ", tags$em("Genome Biology"), ", 11, R25; Law CW, et al. (2014). voom. ", tags$em("Genome Biology"), ", 15, R29; Ritchie ME, et al. (2015). limma. ", tags$em("Nucleic Acids Research"), ", 43(7), e47; Smyth GK (2004). ", tags$em("Stat Appl Genet Mol Biol"), ", 3, Article 3."),
          tags$li(strong("Multiple-testing correction: "), "Benjamini Y, Hochberg Y (1995). ", tags$em("J R Stat Soc B"), ", 57(1), 289-300."),
          tags$li(strong("ROC / AUC confidence intervals: "), "DeLong ER, et al. (1988). ", tags$em("Biometrics"), ", 44, 837-845 (n ≥ 20); Carpenter J, Bithell J (2000). ", tags$em("Stat Med"), ", 19, 1141-1164 (n < 20)."),
          tags$li(strong("Elastic net (glmnet), random forest, SVM, caret: "), "same as Diagnostic Model's References box - Friedman, Hastie & Tibshirani (2010); Zou & Hastie (2005); Breiman (2001); Cortes & Vapnik (1995); Kuhn (2008)."),
          tags$li(strong("Selection-bias caveat on any fixed, blood-derived panel: "), "Ambroise C, McLachlan GJ (2002). ", tags$em("PNAS"), ", 99(10), 6562-6566.")
        ),
        p(class = "submodule-desc", strong("Ask ArthOChat"), " for a plain-language walkthrough or a live citation.")
        )
      )
    })

    ## Saved into results$crosstissue independently per sex, as each
    ## finishes - same save pattern as Diagnostic Model's save_result().
    save_result <- function(sex_label, r) {
      d <- r$discovery
      sig_cut <- input$sig_cutoff %||% 0.05
      is_bio <- ct_biomarker_flag(d, sig_cut)
      results$crosstissue <- utils::modifyList(
        results$crosstissue %||% list(),
        setNames(list(list(
          n_input = r$n_input, n_present = r$n_present, n_samples = r$n_samples, n_pos = r$n_pos, n_neg = r$n_neg,
          n_biomarkers = sum(is_bio), biomarker_genes = d$gene[is_bio],
          n_concordant = sum(d$concordant, na.rm = TRUE),
          n_significant_concordant = sum(d$concordant & d$syn_adjP < sig_cut, na.rm = TRUE),
          lr_apparent_auc = round(r$lr$full_auc, 3),
          lr_pooled_cv_auc = if (isTRUE(r$lr$cv$pooled$available)) round(r$lr$cv$pooled$auc, 3) else NA_real_,
          enet_apparent_auc = round(r$enet$full_auc, 3),
          enet_pooled_cv_auc = if (isTRUE(r$enet$cv$pooled$available)) round(r$enet$cv$pooled$auc, 3) else NA_real_,
          rf_apparent_auc = round(r$rf$full_auc, 3),
          rf_pooled_cv_auc = if (isTRUE(r$rf$cv$pooled$available)) round(r$rf$cv$pooled$auc, 3) else NA_real_,
          svm_apparent_auc = round(r$svm$full_auc, 3),
          svm_pooled_cv_auc = if (isTRUE(r$svm$cv$pooled$available)) round(r$svm$cv$pooled$auc, 3) else NA_real_,
          genes = r$genes, dataset_source = r$dataset_label
        )), sex_label)
      )
      showNotification(
        sprintf("%s: %d validated cross-tissue biomarker%s found (of %d genes present in synovium).",
                tools::toTitleCase(sex_label), sum(is_bio), if (sum(is_bio) == 1) "" else "s", r$n_present),
        type = "message", duration = 6
      )
    }
    observeEvent(ct_result_pooled(), save_result("pooled", ct_result_pooled()))
    observeEvent(ct_result_female(), save_result("female", ct_result_female()))
    observeEvent(ct_result_male(), save_result("male", ct_result_male()))

    output$saved_runs_ui <- renderUI({
      res_p <- tryCatch(ct_result_pooled(), error = function(e) NULL)
      res_f <- tryCatch(ct_result_female(), error = function(e) NULL)
      res_m <- tryCatch(ct_result_male(), error = function(e) NULL)
      sig_cut <- input$sig_cutoff %||% 0.05
      status_row <- function(sex, r) {
        if (is.null(r)) {
          tags$li(icon("circle-minus", style = "color: #8A929C;"), sprintf(" %s - not run yet", sex))
        } else {
          n_bio <- sum(ct_biomarker_flag(r$discovery, sig_cut))
          tags$li(icon("award", style = sprintf("color: %s;", ARTHOMIX_STATUS$good)), strong(sprintf(" %s: ", sex)),
                  sprintf("%d biomarker%s", n_bio, if (n_bio == 1) "" else "s"))
        }
      }
      tags$ul(style = "padding-left: 18px; margin-bottom: 0; list-style: none;",
              status_row("Pooled (All)", res_p), status_row("Female", res_f), status_row("Male", res_m))
    })

    res_sex <- function(sex_label) reactive({
      fr <- switch(sex_label, female = ct_result_female, male = ct_result_male, ct_result_pooled)
      tryCatch(fr(), error = function(e) NULL)
    })

    ## -----------------------------------------------------------------
    ## Discovery & Concordance outputs, per sex.
    ## -----------------------------------------------------------------

    register_discovery_outputs <- function(sex_label, res) {
      not_yet_note <- function() div(class = "empty-note", icon("circle-info"), "Not run yet, or the last run failed validation - check above.")

      ## Single source of truth for "is this gene a cross-tissue biomarker"
      ## (ct_biomarker_flag(), module header) - the KPI tiles, both plots and
      ## the table all mark the exact same genes, so a gene highlighted in
      ## one view is never absent from another.
      disc_marked <- reactive({
        r <- res(); req(r)
        d <- r$discovery
        d$biomarker <- ct_biomarker_flag(d, input$sig_cutoff %||% 0.05)
        d
      })

      output[[paste0(sex_label, "_disc_stats")]] <- renderUI({
        r <- res(); if (is.null(r)) return(not_yet_note())
        d <- disc_marked()
        n_present <- sum(d$present); n_total <- nrow(d)
        n_bio <- sum(d$biomarker)
        n_conc <- sum(d$concordant, na.rm = TRUE)
        med_auc <- stats::median(d$auc_bestdir[d$biomarker], na.rm = TRUE)
        fluidRow(
          valueBox(n_bio, "Validated cross-tissue biomarkers", icon = icon("award"), color = "green", width = 3),
          valueBox(sprintf("%d / %d", n_present, n_total), "Panel genes present in synovium", icon = icon("dna"), color = "light-blue", width = 3),
          valueBox(sprintf("%d / %d", n_conc, n_present), "Direction-concordant with blood", icon = icon("arrows-turn-to-dots"), color = "purple", width = 3),
          valueBox(if (is.na(med_auc)) "N/A" else sprintf("%.3f", med_auc), "Median AUC (biomarkers)", icon = icon("chart-line"), color = "aqua", width = 3)
        )
      })

      ## Direction concordance: every present gene plotted, but only
      ## validated biomarkers get a label and the highlight color - the
      ## story is which genes cleared the bar, not a wall of gene names.
      output[[paste0(sex_label, "_concordance_plot")]] <- renderPlot({
        d <- disc_marked(); d <- d[d$present & !is.na(d$blood_log2FC), , drop = FALSE]
        req(nrow(d) > 0)
        d$Status <- factor(ifelse(d$biomarker, "Validated biomarker", "Not validated"), levels = c("Validated biomarker", "Not validated"))
        ggplot(d, aes(x = blood_log2FC, y = syn_log2FC, color = Status)) +
          geom_hline(yintercept = 0, color = ARTHOMIX_COLORS$axis, linewidth = 0.3) +
          geom_vline(xintercept = 0, color = ARTHOMIX_COLORS$axis, linewidth = 0.3) +
          geom_point(aes(size = Status), alpha = 0.9) +
          ggrepel::geom_text_repel(data = d[d$biomarker, , drop = FALSE], aes(label = gene), size = 3.4,
                                    color = ARTHOMIX_COLORS$ink, fontface = "bold", show.legend = FALSE, max.overlaps = 30) +
          scale_color_manual(values = c(`Validated biomarker` = ARTHOMIX_STATUS$good, `Not validated` = ARTHOMIX_COLORS$ink_muted)) +
          scale_size_manual(values = c(`Validated biomarker` = 3.6, `Not validated` = 2.2), guide = "none") +
          labs(title = "Direction concordance", subtitle = "Blood vs. synovium log2 fold change",
               x = "Blood log2FC (RA vs HC)", y = "Synovium log2FC (RA vs Normal)", color = NULL) +
          theme_arthomix(base_size = 12)
      }, alt = sprintf("Scatter plot comparing each %s panel gene's blood log2 fold change against its synovium log2 fold change, with validated cross-tissue biomarkers highlighted and labelled.", sex_label))

      ## Cross-tissue biomarker identification: a ranked lollipop instead of
      ## a bar chart - thinner marks, only biomarkers labelled with their
      ## AUC, a reference line at chance (0.5) and at the biomarker cutoff
      ## (CT_BIOMARKER_AUC_MIN).
      output[[paste0(sex_label, "_geneauc_plot")]] <- renderPlot({
        d <- disc_marked(); d <- d[d$present, , drop = FALSE]
        req(nrow(d) > 0)
        use_train <- identical(input$orient_view, "trainorient")
        d$auc_show <- if (use_train) d$auc_trainorient else d$auc_bestdir
        d <- d[!is.na(d$auc_show), , drop = FALSE]
        req(nrow(d) > 0)
        d$Status <- factor(ifelse(d$biomarker, "Validated biomarker", "Not validated"), levels = c("Validated biomarker", "Not validated"))
        d$gene <- factor(d$gene, levels = d$gene[order(d$auc_show)])
        ggplot(d, aes(x = auc_show, y = gene, color = Status)) +
          geom_vline(xintercept = 0.5, linetype = "dashed", color = ARTHOMIX_COLORS$axis, linewidth = 0.4) +
          geom_vline(xintercept = CT_BIOMARKER_AUC_MIN, linetype = "dashed", color = ARTHOMIX_STATUS$good, linewidth = 0.4) +
          geom_segment(aes(x = 0.5, xend = auc_show, yend = gene), linewidth = 0.7) +
          geom_point(aes(size = Status)) +
          geom_text(data = d[d$biomarker, , drop = FALSE], aes(label = sprintf("%.2f", auc_show)),
                    hjust = -0.4, size = 3.2, color = ARTHOMIX_COLORS$ink, show.legend = FALSE) +
          scale_color_manual(values = c(`Validated biomarker` = ARTHOMIX_STATUS$good, `Not validated` = ARTHOMIX_COLORS$ink_muted)) +
          scale_size_manual(values = c(`Validated biomarker` = 3.4, `Not validated` = 2.4), guide = "none") +
          scale_x_continuous(limits = c(0.3, 1.05), breaks = c(0.5, CT_BIOMARKER_AUC_MIN, 1)) +
          labs(title = "Cross-tissue biomarkers", subtitle = sprintf("Synovium AUC (%s)", if (use_train) "train-fixed" else "best-direction"),
               x = NULL, y = NULL, color = NULL) +
          theme_arthomix(base_size = 12) + theme(legend.position = "bottom")
      }, alt = sprintf("Ranked dot plot of each %s panel gene's synovium AUC, with validated cross-tissue biomarkers highlighted and labelled against a 0.70 validation cutoff.", sex_label))

      output[[paste0(sex_label, "_disc_table")]] <- DT::renderDataTable({
        d <- disc_marked()
        d2 <- data.frame(
          gene = d$gene, cross_tissue_biomarker = ifelse(d$biomarker, "✓", ""), present = d$present,
          synovium_log2FC = round(d$syn_log2FC, 3), synovium_adjP = signif(d$syn_adjP, 3),
          blood_log2FC = round(d$blood_log2FC, 3), concordant = d$concordant,
          auc_bestdir = round(d$auc_bestdir, 3), auc_trainorient = round(d$auc_trainorient, 3)
        )
        d2 <- d2[order(-d$biomarker, -d$auc_bestdir), ]
        DT::datatable(d2, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
      })
      output[[paste0(sex_label, "_disc_download")]] <- downloadHandler(
        filename = function() sprintf("%s_synovium_discovery.csv", sex_label),
        content = function(file) write.csv(res()$discovery, file, row.names = FALSE)
      )
    }

    ## -----------------------------------------------------------------
    ## Panel classifier outputs (Full Fit + Cross-Validated), per sex x
    ## technique - same output-name/box structure as Diagnostic Model,
    ## driven off ct_fit_sex()'s result shape instead of diag_fit_sex()'s.
    ## -----------------------------------------------------------------

    build_full_perf_table <- function(rr) {
      do.call(rbind, list(
        data.frame(dataset = "Full synovium fit (apparent)", metric = "AUC", value = sprintf("%.3f", rr$full_auc), stringsAsFactors = FALSE),
        data.frame(dataset = "Full synovium fit (apparent)", metric = "Threshold (prob., Youden)", value = sprintf("%.3f", rr$best$threshold), stringsAsFactors = FALSE),
        data.frame(dataset = "Full synovium fit (apparent)", metric = "Sensitivity", value = sprintf("%.3f", rr$perf_full$sensitivity), stringsAsFactors = FALSE),
        data.frame(dataset = "Full synovium fit (apparent)", metric = "Specificity", value = sprintf("%.3f", rr$perf_full$specificity), stringsAsFactors = FALSE),
        data.frame(dataset = "Full synovium fit (apparent)", metric = "Accuracy", value = sprintf("%.3f", rr$perf_full$accuracy), stringsAsFactors = FALSE)
      ))
    }

    build_cv_perf_table <- function(rr) {
      fa <- rr$cv$fold_auc
      base_rows <- list(
        data.frame(dataset = sprintf("%d-fold CV (per-fold)", length(fa)), metric = "Mean AUC", value = sprintf("%.3f", mean(fa, na.rm = TRUE)), stringsAsFactors = FALSE),
        data.frame(dataset = sprintf("%d-fold CV (per-fold)", length(fa)), metric = "SD AUC", value = sprintf("%.3f", stats::sd(fa, na.rm = TRUE)), stringsAsFactors = FALSE)
      )
      pooled <- rr$cv$pooled
      pooled_rows <- if (isTRUE(pooled$available)) {
        auc_str <- paste0(sprintf("%.3f (%.3f-%.3f)", pooled$auc, pooled$ci_lo, pooled$ci_hi), diag_separation_note(pooled))
        list(
          data.frame(dataset = "Pooled out-of-fold (headline)", metric = "AUC (95% CI)", value = auc_str, stringsAsFactors = FALSE),
          data.frame(dataset = "Pooled out-of-fold (headline)", metric = "n samples (n positive)", value = sprintf("%d (%d)", pooled$n, pooled$n_pos), stringsAsFactors = FALSE)
        )
      } else {
        list(data.frame(dataset = "Pooled out-of-fold (headline)", metric = "Status", value = pooled$reason %||% "Unavailable.", stringsAsFactors = FALSE))
      }
      do.call(rbind, c(base_rows, pooled_rows))
    }

    build_model_bundle <- function(r, key, sex_label) {
      rr <- r[[key]]
      pooled <- rr$cv$pooled
      list(
        model_type = rr$model_type, model_label = rr$label, model = rr$model,
        sex = sex_label, genes = r$genes,
        contrast = sprintf("%s vs %s (%s)", r$grp_levels[2], r$grp_levels[1], r$dataset_label),
        hyperparams = switch(rr$model_type,
          lr = list(),
          enet = list(alpha = rr$alpha, lambda_choice = rr$lambda_choice, lambda_used = rr$lambda_used),
          rf = list(ntree = rr$ntree, mtry = rr$mtry),
          svm = list(kernel = rr$kernel, cost = rr$cost)
        ),
        full_fit_auc = rr$full_auc,
        cv_pooled_auc = if (isTRUE(pooled$available)) pooled$auc else NA_real_,
        cv_mean_fold_auc = mean(rr$cv$fold_auc, na.rm = TRUE),
        scoring_note = paste(
          "Fit on the FULL synovium sex-subset (apparent/resubstitution) - not a held-out estimate; see cv_pooled_auc for the",
          "out-of-fold estimate this project's own methodology reports. To score new samples with `model`: for each gene in",
          "`genes`, z-score it across the NEW dataset's own samples (independent of this training set's own mean/SD); a gene",
          "in `genes` absent from the new data becomes a column of 0. Align columns to `genes` order (via make.names if",
          "needed), then: logistic regression - predict(model, newdata = as.data.frame(X), type = \"response\"); elastic net -",
          "predict(model, newx = X, s = hyperparams$lambda_choice, type = \"response\"); random forest -",
          "predict(model, X, type = \"prob\")[, \"RA\"]; SVM -",
          "attr(predict(model, X, probability = TRUE), \"probabilities\")[, \"RA\"]."
        ),
        trained_at = as.character(Sys.time())
      )
    }

    register_sex_model_outputs <- function(sex_label, res) {
      sex_color <- switch(sex_label, female = "#1a7a3c", male = "#7a4a26", "#2c6fbb")
      cv_color <- ARTHOMIX_COLORS$orange
      not_yet_note <- function() div(class = "empty-note", icon("circle-info"), "Not run yet, or the last run failed validation - check above.")

      lapply(CT_TECHNIQUES, function(tech) {
        key <- tech$key; prefix <- paste0(sex_label, "_", key)

        ## ---- Full Fit tab ----
        output[[paste0(prefix, "_full_stats")]] <- renderUI({
          r <- res(); if (is.null(r)) return(not_yet_note())
          rr <- r[[key]]
          fluidRow(
            valueBox(sprintf("%.3f", rr$full_auc), "Apparent AUC (optimistic upper bound)", icon = icon("triangle-exclamation"), color = "light-blue", width = 4),
            valueBox(sprintf("%.3f", mean(rr$cv$fold_auc, na.rm = TRUE)),
                     sprintf("%d-fold CV AUC (± %.3f)", length(rr$cv$fold_auc), stats::sd(rr$cv$fold_auc, na.rm = TRUE)),
                     icon = icon("layer-group"), color = "purple", width = 4),
            valueBox(diag_hyperparam_value(rr), "Hyperparameter", icon = icon("sliders"), color = "yellow", width = 4)
          )
        })
        output[[paste0(prefix, "_full_roc_plot")]] <- renderPlot({
          r <- res(); req(r); rr <- r[[key]]
          pROC::ggroc(rr$roc_full, color = sex_color, size = 1) +
            geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), linetype = "dashed", color = "#a9b1bb") +
            labs(x = "Specificity", y = "Sensitivity", title = sprintf("Apparent AUC = %.3f", rr$full_auc)) +
            theme_arthomix(base_size = 12)
        }, alt = sprintf("ROC curve for the %s %s apparent (resubstitution) fit on the full synovium sex-subset.", sex_label, tech$label))
        output[[paste0(prefix, "_cv_fold_plot")]] <- renderPlot({
          r <- res(); req(r); rr <- r[[key]]
          df <- data.frame(fold = factor(seq_along(rr$cv$fold_auc)), auc = rr$cv$fold_auc)
          ggplot(df, aes(x = fold, y = auc)) +
            geom_col(fill = sex_color) +
            geom_hline(yintercept = mean(rr$cv$fold_auc, na.rm = TRUE), linetype = "dashed", color = ARTHOMIX_COLORS$red) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(x = "Fold", y = "CV AUC (out-of-fold)") + theme_arthomix(base_size = 12)
        }, alt = sprintf("Bar chart of the %s %s cross-validated AUC in each outer fold, with the mean marked.", sex_label, tech$label))
        output[[paste0(prefix, "_tuning_plot")]] <- renderPlot({
          r <- res(); req(r); rr <- r[[key]]
          ts <- rr$tuning_search
          if (is.null(ts) || nrow(ts) < 2) {
            msg <- if (identical(rr$model_type, "lr")) {
              "Plain logistic regression has no\nregularisation path - nothing to tune."
            } else {
              "Manual value used -\nno tuning grid was searched."
            }
            return(ggplot() + annotate("text", x = 0, y = 0, label = msg, size = 4.2, color = "#64748B") + theme_void())
          }
          switch(rr$model_type,
            enet = ggplot(ts, aes(x = factor(alpha), y = cv_deviance, fill = chosen)) +
              geom_col() +
              scale_fill_manual(values = c(`TRUE` = sex_color, `FALSE` = "#CBD5E1"), guide = "none") +
              labs(x = "Alpha", y = "CV deviance (lower is better)", title = sprintf("Selected alpha = %.2f", rr$alpha)) +
              theme_arthomix(base_size = 12),
            rf = ggplot(ts, aes(x = mtry, y = ROC)) +
              geom_line(color = sex_color) + geom_point(color = sex_color, size = 2) +
              geom_vline(xintercept = rr$mtry, linetype = "dashed", color = ARTHOMIX_COLORS$red) +
              labs(x = "mtry", y = "CV ROC AUC", title = sprintf("Selected mtry = %d", rr$mtry)) +
              theme_arthomix(base_size = 12),
            svm = ggplot(ts, aes(x = cost, y = error)) +
              geom_line(color = sex_color) + geom_point(color = sex_color, size = 2) +
              scale_x_log10() +
              geom_vline(xintercept = rr$cost, linetype = "dashed", color = ARTHOMIX_COLORS$red) +
              labs(x = "Cost (log scale)", y = "CV classification error", title = sprintf("Selected cost = %s", format(rr$cost, trim = TRUE))) +
              theme_arthomix(base_size = 12)
          )
        }, alt = sprintf("Hyperparameter tuning grid for the %s %s model, with the selected value marked, or a note explaining why no grid was searched.", sex_label, tech$label))
        output[[paste0(prefix, "_full_table")]] <- DT::renderDataTable({
          r <- res(); req(r)
          DT::datatable(build_full_perf_table(r[[key]]), rownames = FALSE, width = "100%",
                        options = list(pageLength = 15, dom = "t", scrollX = TRUE, autoWidth = FALSE), class = "stripe hover compact")
        })
        output[[paste0(prefix, "_full_download")]] <- downloadHandler(
          filename = function() sprintf("%s_%s_crosstissue_fullfit_performance.csv", sex_label, key),
          content = function(file) write.csv(build_full_perf_table(res()[[key]]), file, row.names = FALSE)
        )
        output[[paste0(prefix, "_model_download")]] <- downloadHandler(
          filename = function() sprintf("%s_%s_crosstissue_model.rds", sex_label, key),
          content = function(file) saveRDS(build_model_bundle(res(), key, sex_label), file)
        )

        ## ---- Cross-Validated tab ----
        output[[paste0(prefix, "_cv_summary")]] <- renderUI({
          r <- res(); if (is.null(r)) return(not_yet_note())
          pooled <- r[[key]]$cv$pooled
          pooled_tile <- if (isTRUE(pooled$available)) {
            valueBox(sprintf("%.3f", pooled$auc), sprintf("Pooled out-of-fold AUC (95%% CI %.3f-%.3f)%s", pooled$ci_lo, pooled$ci_hi, diag_separation_note(pooled)),
                     icon = icon("flask"), color = "light-blue", width = 6)
          } else {
            valueBox("N/A", "Pooled out-of-fold AUC", icon = icon("flask"), color = "red", width = 6)
          }
          n_tile <- valueBox(if (isTRUE(pooled$available)) pooled$n else "-", "Samples with an out-of-fold prediction", icon = icon("users"), color = "purple", width = 6)
          tagList(
            fluidRow(pooled_tile, n_tile),
            if (!isTRUE(pooled$available)) p(class = "submodule-desc", style = "font-size: 12.5px;", icon("circle-info"), pooled$reason %||% "Pooled estimate unavailable.") else NULL
          )
        })
        output[[paste0(prefix, "_cv_roc_plot")]] <- renderPlot({
          r <- res(); req(r); rr <- r[[key]]
          req(isTRUE(rr$cv$pooled$available))
          pROC::ggroc(rr$cv$pooled$roc, color = cv_color, size = 1) +
            geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), linetype = "dashed", color = "#a9b1bb") +
            labs(x = "Specificity", y = "Sensitivity", title = sprintf("Pooled out-of-fold AUC = %.3f", rr$cv$pooled$auc)) +
            theme_arthomix(base_size = 12)
        }, alt = sprintf("ROC curve built from every %s %s out-of-fold prediction pooled across the outer cross-validation folds.", sex_label, tech$label))
        output[[paste0(prefix, "_cv_table")]] <- DT::renderDataTable({
          r <- res(); req(r)
          DT::datatable(build_cv_perf_table(r[[key]]), rownames = FALSE, width = "100%",
                        options = list(pageLength = 15, dom = "t", scrollX = TRUE, autoWidth = FALSE), class = "stripe hover compact")
        })
        output[[paste0(prefix, "_cv_download")]] <- downloadHandler(
          filename = function() sprintf("%s_%s_crosstissue_cv_performance.csv", sex_label, key),
          content = function(file) write.csv(build_cv_perf_table(res()[[key]]), file, row.names = FALSE)
        )
      })

      output[[paste0(sex_label, "_full_compare_table")]] <- DT::renderDataTable({
        r <- res(); req(r)
        df <- do.call(rbind, lapply(CT_TECHNIQUES, function(tech) {
          rr <- r[[tech$key]]
          data.frame(model = tech$label, apparent_auc = round(rr$full_auc, 3),
                     cv_mean_fold_auc = round(mean(rr$cv$fold_auc, na.rm = TRUE), 3),
                     cv_sd_fold_auc = round(stats::sd(rr$cv$fold_auc, na.rm = TRUE), 3),
                     stringsAsFactors = FALSE)
        }))
        DT::datatable(df, rownames = FALSE, width = "100%", options = list(pageLength = 5, dom = "t", scrollX = TRUE, autoWidth = FALSE), class = "stripe hover compact")
      })
      output[[paste0(sex_label, "_cv_compare_table")]] <- DT::renderDataTable({
        r <- res(); req(r)
        df <- do.call(rbind, lapply(CT_TECHNIQUES, function(tech) {
          rr <- r[[tech$key]]; pooled <- rr$cv$pooled
          data.frame(
            model = tech$label,
            pooled_cv_auc = if (isTRUE(pooled$available)) round(pooled$auc, 3) else NA_real_,
            ci_lo = if (isTRUE(pooled$available)) round(pooled$ci_lo, 3) else NA_real_,
            ci_hi = if (isTRUE(pooled$available)) round(pooled$ci_hi, 3) else NA_real_,
            note = trimws(diag_separation_note(pooled)),
            stringsAsFactors = FALSE
          )
        }))
        DT::datatable(df, rownames = FALSE, width = "100%", options = list(pageLength = 5, dom = "t", scrollX = TRUE, autoWidth = FALSE), class = "stripe hover compact")
      })
      output[[paste0(sex_label, "_full_result_line")]] <- renderUI({
        r <- res(); if (is.null(r)) return(NULL)
        p(strong("Result: "),
          sprintf("%d genes (%d present in synovium), %d samples (%d RA vs %d Normal) → apparent AUC logistic regression %.3f / elastic net %.3f / random forest %.3f / SVM %.3f.",
                  r$n_input, r$n_present, r$n_samples, r$n_pos, r$n_neg, r$lr$full_auc, r$enet$full_auc, r$rf$full_auc, r$svm$full_auc))
      })
    }

    ## -----------------------------------------------------------------
    ## Cross-Dataset Comparison outputs, per sex - reads results$diagnostic
    ## (written by mod_diagnostic.R) rather than recomputing or storing any
    ## copy of the blood data; absent if Diagnostic Model hasn't been run
    ## this session for that sex.
    ## -----------------------------------------------------------------

    CROSSDATA_ORDER <- c("Blood train (full fit)", "Blood train (k-fold CV)", "Synovium (apparent)", "Synovium (pooled CV)")
    CROSSDATA_PAL <- c(`Blood train (full fit)` = "#1b6ca8", `Blood train (k-fold CV)` = "#7fb2dd",
                        `Synovium (apparent)` = "#c0392b", `Synovium (pooled CV)` = "#e8a598")

    register_crossdata_outputs <- function(sex_label, res) {
      output[[paste0(sex_label, "_crossdata_note")]] <- renderUI({
        r <- res()
        if (is.null(r)) {
          return(div(class = "empty-note", icon("circle-info"),
                      sprintf("Run %s in Discovery, Full Fit or Cross-Validated first.", tools::toTitleCase(sex_label))))
        }
        diag <- results$diagnostic[[sex_label]]
        if (is.null(diag)) {
          return(div(class = "empty-note", icon("circle-info"),
                      sprintf("No Diagnostic Model result for %s yet - only synovium is shown below.", tools::toTitleCase(sex_label))))
        }
        overlap <- length(intersect(diag$genes, r$genes))
        tagList(
          div(class = "empty-note",
              title = "Coefficients are re-estimated within synovium using the same gene identities - not a transfer of the blood-fitted model (Section 2.11.1).",
              icon("circle-info"),
              sprintf("%d of %d synovium-panel genes also appear in the %d-gene blood panel.", overlap, length(r$genes), length(diag$genes)))
        )
      })
      output[[paste0(sex_label, "_crossdata_plot")]] <- renderPlot({
        r <- res(); req(r)
        diag <- results$diagnostic[[sex_label]]; req(diag)
        df <- do.call(rbind, lapply(CT_TECHNIQUES, function(tech) {
          rr <- r[[tech$key]]
          rbind(
            data.frame(model = tech$label, dataset = "Blood train (full fit)", auc = diag[[paste0(tech$key, "_auc")]]),
            data.frame(model = tech$label, dataset = "Blood train (k-fold CV)", auc = diag[[paste0(tech$key, "_cv_auc")]]),
            data.frame(model = tech$label, dataset = "Synovium (apparent)", auc = rr$full_auc),
            data.frame(model = tech$label, dataset = "Synovium (pooled CV)", auc = if (isTRUE(rr$cv$pooled$available)) rr$cv$pooled$auc else NA_real_)
          )
        }))
        df$dataset <- factor(df$dataset, levels = CROSSDATA_ORDER)
        ggplot(df, aes(x = model, y = auc, fill = dataset)) +
          geom_col(position = position_dodge(0.75), width = 0.7, na.rm = TRUE) +
          geom_hline(yintercept = 0.5, linetype = "dashed", color = "#8a929c") +
          scale_fill_manual(values = CROSSDATA_PAL, name = NULL) +
          coord_cartesian(ylim = c(0, 1)) +
          labs(x = NULL, y = "AUC") + theme_arthomix(base_size = 12) + theme(legend.position = "bottom")
      }, alt = sprintf("Grouped bar chart comparing, per model, %s's blood-train full-fit and cross-validated AUC against this synovium panel's apparent and pooled cross-validated AUC.", tools::toTitleCase(sex_label)))
      output[[paste0(sex_label, "_crossdata_table")]] <- DT::renderDataTable({
        r <- res(); req(r)
        diag <- results$diagnostic[[sex_label]]; req(diag)
        df <- do.call(rbind, lapply(CT_TECHNIQUES, function(tech) {
          rr <- r[[tech$key]]
          data.frame(model = tech$label,
                     blood_train_full = round(diag[[paste0(tech$key, "_auc")]], 3),
                     blood_train_cv = round(diag[[paste0(tech$key, "_cv_auc")]], 3),
                     synovium_apparent = round(rr$full_auc, 3),
                     synovium_pooled_cv = if (isTRUE(rr$cv$pooled$available)) round(rr$cv$pooled$auc, 3) else NA_real_,
                     stringsAsFactors = FALSE)
        }))
        DT::datatable(df, rownames = FALSE, width = "100%", options = list(pageLength = 5, dom = "t", scrollX = TRUE, autoWidth = FALSE), class = "stripe hover compact")
      })
    }

    register_discovery_outputs("pooled", res_sex("pooled"))
    register_discovery_outputs("female", res_sex("female"))
    register_discovery_outputs("male", res_sex("male"))
    register_sex_model_outputs("pooled", res_sex("pooled"))
    register_sex_model_outputs("female", res_sex("female"))
    register_sex_model_outputs("male", res_sex("male"))
    register_crossdata_outputs("pooled", res_sex("pooled"))
    register_crossdata_outputs("female", res_sex("female"))
    register_crossdata_outputs("male", res_sex("male"))
  })
}
