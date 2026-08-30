## R/mod_diagnostic.R
## Diagnostic Model submodule: fits logistic regression, elastic net, random
## forest and SVM classifiers on a user-chosen gene panel, sex-stratified.
## Model Training splits each sex into Train/Test once and tunes on Train
## only; Model Testing scores that same fit once on the held-out Test split.

mod_diagnostic_config <- list(
  id = "diagnostic", group = "Biomarker modeling",
  title = "Diagnostic Model",
  description = "Diagnostic model using logistic regression, elastic net, random forest and SVM diagnostic models, by sex.",
  icon = "stethoscope"
)

## ---------------------------------------------------------------------------
## Shared fitting helpers
## ---------------------------------------------------------------------------

## Per-gene (row) z-score, independent of any other dataset's mean/SD.
diag_zrows <- function(M) {
  t(apply(M, 1, function(v) {
    s <- stats::sd(v, na.rm = TRUE)
    if (is.na(s) || s == 0) rep(0, length(v)) else (v - mean(v, na.rm = TRUE)) / s
  }))
}

DIAG_SVM_COST_GRID <- c(0.01, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16)

DIAG_DEFAULT_PARAMS <- list(
  test_frac = 0.3,
  ## Class weighting shared across all four models: "equal" = unweighted,
  ## "balanced" = inverse-frequency, "manual" = fixed comparison:reference ratio.
  class_weight_mode = "equal", class_weight_ratio = 1,
  lr_cv_folds = 5,
  enet_cv_folds = 5, enet_alpha_grid = c(0.1, 0.3, 0.5, 0.7, 0.9, 1.0), enet_lambda_choice = "lambda.min",
  enet_nlambda = 100, enet_type_measure = "deviance",
  rf_cv_folds = 5, rf_ntree = 1000, rf_mtry_mode = "auto", rf_mtry_manual = NULL,
  rf_nodesize = 1, rf_maxnodes = NULL,
  svm_cv_folds = 5, svm_kernel = "linear", svm_cost_mode = "auto", svm_cost_manual = 1, svm_cost_grid = DIAG_SVM_COST_GRID,
  svm_gamma_mode = "auto", svm_gamma_manual = 1, svm_degree = 3, svm_tolerance = 0.001
)

## Named (reference, comparison) class weights for randomForest/e1071/glm weights=.
diag_class_weight_levels <- function(y, mode, ratio) {
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

diag_obs_weights <- function(y, mode, ratio) {
  wl <- diag_class_weight_levels(y, mode, ratio)
  unname(wl[as.character(y)])
}

## Stratified train/test split of one sex's pool, fixed single split (seed 1234).
diag_split_train_test <- function(y_full, test_frac, seed = 1234) {
  set.seed(seed)
  train_idx <- as.integer(caret::createDataPartition(y_full, p = 1 - test_frac, list = FALSE))
  list(train = train_idx, test = setdiff(seq_along(y_full), train_idx))
}

## k-fold CV AUC for an already-tuned classifier; rescales per fold with
## fold-train mean/SD only (leakage-free), refitting via refit_fn/predict_fn.
diag_cv_auc <- function(Xraw, y, n_folds, refit_fn, predict_fn, seed = 1234) {
  nf <- max(2, min(n_folds, min(table(y))))
  set.seed(seed)
  folds <- caret::createFolds(y, k = nf)
  vapply(folds, function(te) {
    tr <- setdiff(seq_along(y), te)
    if (length(unique(y[tr])) < 2) return(NA_real_)
    mu <- colMeans(Xraw[tr, , drop = FALSE])
    sg <- apply(Xraw[tr, , drop = FALSE], 2, stats::sd)
    sg[is.na(sg) | sg == 0] <- 1
    Ztr <- scale(Xraw[tr, , drop = FALSE], center = mu, scale = sg)
    Zte <- scale(Xraw[te, , drop = FALSE], center = mu, scale = sg)
    fit_i <- tryCatch(refit_fn(Ztr, y[tr]), error = function(e) NULL)
    if (is.null(fit_i)) return(NA_real_)
    p <- tryCatch(predict_fn(fit_i, Zte), error = function(e) NULL)
    if (is.null(p)) return(NA_real_)
    roc_i <- tryCatch(pROC::roc(y[te], p, quiet = TRUE, levels = levels(y), direction = "<"), error = function(e) NULL)
    if (is.null(roc_i)) NA_real_ else as.numeric(pROC::auc(roc_i))
  }, numeric(1))
}

## AUC + CI: DeLong when n >= 20, else stratified bootstrap (seed 1234, 2000 reps).
diag_auc_ci <- function(r) {
  n <- length(r$cases) + length(r$controls)
  ci <- if (n < 20) {
    set.seed(1234)
    suppressWarnings(tryCatch(as.numeric(pROC::ci.auc(r, method = "bootstrap", boot.n = 2000)), error = function(e) c(NA, NA, NA)))
  } else {
    suppressWarnings(tryCatch(as.numeric(pROC::ci.auc(r)), error = function(e) c(NA, NA, NA)))
  }
  c(auc = as.numeric(pROC::auc(r)), lo = ci[1], hi = ci[3])
}

## Publication-style ROC plot (1-specificity vs sensitivity, journal convention
## per Chen et al. 2021/2022 Fig. 4a). Matches diag_roc_plot_traintest()'s
## boxed bottom-right legend (color+linetype swatch next to the AUC/CI text,
## Chance drawn as a labelled series) instead of a floating text annotation.
diag_roc_plot_pub <- function(roc_obj, color, title = NULL, ci = NULL) {
  co <- pROC::coords(roc_obj, "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
  df <- data.frame(fpr = 1 - co$specificity, tpr = co$sensitivity)
  df <- df[order(df$fpr, df$tpr), ]
  auc_val <- as.numeric(pROC::auc(roc_obj))
  ci_val <- if (!is.null(ci) && all(is.finite(ci))) ci else diag_auc_ci(roc_obj)[c("lo", "auc", "hi")]
  curve_label <- if (all(is.finite(ci_val))) {
    sprintf("AUC = %.2f (95%% CI %.2f–%.2f)", auc_val, ci_val[1], ci_val[3])
  } else {
    sprintf("AUC = %.2f", auc_val)
  }
  chance_label <- "Chance (AUC = 0.50)"
  series_levels <- c(curve_label, chance_label)
  df$series <- factor(curve_label, levels = series_levels)
  chance_df <- data.frame(fpr = c(0, 1), tpr = c(0, 1), series = factor(chance_label, levels = series_levels))
  colors <- setNames(c(color, "#9CA3AF"), series_levels)
  linetypes <- setNames(c("solid", "dashed"), series_levels)

  ggplot() +
    geom_line(data = chance_df, aes(x = fpr, y = tpr, color = series, linetype = series), linewidth = 0.7) +
    geom_line(data = df, aes(x = fpr, y = tpr, color = series, linetype = series), linewidth = 1.15) +
    scale_color_manual(name = NULL, values = colors, breaks = series_levels) +
    scale_linetype_manual(name = NULL, values = linetypes, breaks = series_levels) +
    scale_x_continuous(name = "1 − Specificity", limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
    scale_y_continuous(name = "Sensitivity", limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
    coord_equal() +
    labs(title = title) +
    guides(color = guide_legend(override.aes = list(linewidth = 1))) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 11, color = "black"),
      panel.border = element_rect(color = "black", linewidth = 0.7, fill = NA),
      legend.position = c(0.98, 0.03), legend.justification = c(1, 0),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
      legend.margin = ggplot2::margin(6, 10, 6, 8),
      legend.text = element_text(size = 10.5, color = "black"),
      legend.key = element_rect(fill = "white"),
      legend.key.width = unit(1.4, "lines")
    )
}

## Train + Test ROC curves overlaid, styled as a journal figure: a single
## bottom-right legend box combines each series' color/linetype swatch with
## its AUC (95% CI) text (and, for Test, a DeLong p-value against Train) -
## rather than a floating color legend plus a separate stack of annotate()
## text labels. Chance (AUC = 0.50) is drawn as a labelled series in the same
## legend instead of an unlabelled diagonal.
diag_roc_plot_traintest <- function(train_roc, test_info, cv_mean = NA_real_, cv_sd = NA_real_, cv_n = NA_integer_,
                                     color_train = ARTHOMIX_COLORS$blue, color_test = ARTHOMIX_COLORS$orange, title = NULL) {
  curve_df <- function(roc_obj, series) {
    co <- pROC::coords(roc_obj, "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
    df <- data.frame(fpr = 1 - co$specificity, tpr = co$sensitivity, series = series, stringsAsFactors = FALSE)
    df[order(df$fpr, df$tpr), ]
  }
  has_test <- isTRUE(test_info$available)

  train_ci <- diag_auc_ci(train_roc)
  train_label <- sprintf("Train: AUC = %.2f (95%% CI %.2f–%.2f)", train_ci["auc"], train_ci["lo"], train_ci["hi"])
  df <- curve_df(train_roc, train_label)

  test_label <- NA_character_
  if (has_test) {
    p_val <- tryCatch(pROC::roc.test(train_roc, test_info$roc, quiet = TRUE)$p.value, error = function(e) NA_real_)
    test_label <- sprintf(
      "Test: AUC = %.2f (95%% CI %.2f–%.2f)%s", test_info$auc, test_info$ci_lo, test_info$ci_hi,
      if (is.finite(p_val)) sprintf(", p = %s", format.pval(p_val, digits = 2, eps = 0.001)) else ""
    )
    df <- rbind(df, curve_df(test_info$roc, test_label))
  }
  chance_label <- "Chance (AUC = 0.50)"
  series_levels <- c(train_label, if (has_test) test_label, chance_label)
  df$series <- factor(df$series, levels = series_levels)
  chance_df <- data.frame(fpr = c(0, 1), tpr = c(0, 1), series = factor(chance_label, levels = series_levels))

  colors <- setNames(c(color_train, if (has_test) color_test, "#9CA3AF"), series_levels)
  ## Train dashed, Test solid (journal convention) so the two curves stay
  ## visually distinct even when they nearly overlap at perfect separation.
  linetypes <- setNames(c("dashed", if (has_test) "solid", "dotted"), series_levels)

  p <- ggplot() +
    geom_line(data = chance_df, aes(x = fpr, y = tpr, color = series, linetype = series), linewidth = 0.7) +
    geom_line(data = df, aes(x = fpr, y = tpr, color = series, linetype = series), linewidth = 1.15) +
    scale_color_manual(name = NULL, values = colors, breaks = series_levels) +
    scale_linetype_manual(name = NULL, values = linetypes, breaks = series_levels) +
    scale_x_continuous(name = "1 − Specificity", limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
    scale_y_continuous(name = "Sensitivity", limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
    coord_equal() +
    labs(title = title, caption = if (is.finite(cv_mean)) sprintf("%d-fold CV AUC: %.2f ± %.2f", cv_n, cv_mean, cv_sd) else NULL) +
    guides(color = guide_legend(override.aes = list(linewidth = 1))) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      plot.caption = element_text(size = 9.5, color = "#4B5563", hjust = 0),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 11, color = "black"),
      panel.border = element_rect(color = "black", linewidth = 0.7, fill = NA),
      legend.position = c(0.98, 0.03), legend.justification = c(1, 0),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
      legend.margin = ggplot2::margin(6, 10, 6, 8),
      legend.text = element_text(size = 10.5, color = "black"),
      legend.key = element_rect(fill = "white"),
      legend.key.width = unit(1.4, "lines")
    )
  p
}

## Flags a CI collapsed against 1.000 at small n as a perfect-separation artefact.
diag_separation_note <- function(ev) {
  if (isTRUE(ev$available) && !is.na(ev$ci_lo) && ev$ci_lo >= 0.999) {
    sprintf(" [separation, n=%d]", ev$n)
  } else ""
}

## Human-readable selected hyperparameter for the Training tab's KPI tile.
diag_hyperparam_value <- function(rr) {
  switch(rr$model_type,
    lr = "none (unpenalized)",
    enet = sprintf("α = %.2f", rr$alpha),
    rf = sprintf("mtry = %d", rr$mtry),
    svm = sprintf("cost = %s", format(rr$cost, trim = TRUE))
  )
}

## Per-gene univariate ROC/AUC on raw expression, plus a Wilcoxon P (Chen et
## al. 2021/2022's hub-gene rule: AUC > 0.85 and P < 0.05).
diag_gene_roc <- function(expr_sub, y) {
  genes <- rownames(expr_sub)
  rocs <- vector("list", length(genes)); names(rocs) <- genes
  aucs <- setNames(numeric(length(genes)), genes)
  pvals <- setNames(numeric(length(genes)), genes)
  for (g in genes) {
    ## direction = "auto": reports whichever of down/up-regulation gives AUC >= 0.5.
    r <- tryCatch(pROC::roc(y, as.numeric(expr_sub[g, ]), quiet = TRUE, levels = levels(y), direction = "auto"), error = function(e) NULL)
    rocs[[g]] <- r
    aucs[g] <- if (is.null(r)) NA_real_ else as.numeric(pROC::auc(r))
    pvals[g] <- tryCatch(stats::wilcox.test(as.numeric(expr_sub[g, ]) ~ y)$p.value, error = function(e) NA_real_)
  }
  list(genes = genes, rocs = rocs, auc = aucs, p = pvals)
}

diag_perf_at_cutoff <- function(prob, y, threshold, positive_level) {
  pred_pos <- prob >= threshold
  obs_pos <- y == positive_level
  list(
    sensitivity = if (sum(obs_pos) > 0) sum(pred_pos & obs_pos) / sum(obs_pos) else NA_real_,
    specificity = if (sum(!obs_pos) > 0) sum(!pred_pos & !obs_pos) / sum(!obs_pos) else NA_real_,
    accuracy = mean(pred_pos == obs_pos)
  )
}

## Fits logistic regression + elastic net + random forest + SVM on one sex's
## full sample pool, splitting into Train/Test first; Test is scored once at
## the end. `expr_full` is genes x samples, raw (not yet z-scored).
diag_fit_sex <- function(expr_full, y_full, params = list()) {
  ## caret::train(classProbs = TRUE) requires factor levels that are valid R
  ## variable names - it make.names()s them internally to build its own
  ## predicted-probability column names, so a raw group label with a space
  ## (e.g. "multiple sclerosis") desyncs from every levels(y)[2]-style lookup
  ## below once caret has already renamed its own columns to
  ## "multiple.sclerosis" - exactly the subscript-out-of-bounds error caret's
  ## own startup warning predicts. Sanitized once here, up front, so every
  ## lookup in this function stays consistent with what caret/randomForest/
  ## e1071 actually produce; callers keep the real group names for their own
  ## display text (stored separately, e.g. fit$ref_group/$comp_group), so
  ## nothing user-visible changes.
  levels(y_full) <- make.names(levels(y_full), unique = TRUE)
  params <- utils::modifyList(DIAG_DEFAULT_PARAMS, params)
  GLOBAL_SEED <- 1234
  genes <- rownames(expr_full)
  safe <- make.names(genes, unique = TRUE)

  split <- diag_split_train_test(y_full, params$test_frac, seed = GLOBAL_SEED)
  validate(need(length(split$train) >= 10 && length(split$test) >= 4,
                "Not enough samples for a train/test split at this ratio - lower the test-set size or provide more samples for this sex."))
  y <- y_full[split$train]
  validate(need(length(unique(y)) == 2 && all(table(y) >= 3),
                "The training split ended up with fewer than 3 samples in one group - lower the test-set size or provide more samples for this sex."))
  ytest <- y_full[split$test]
  validate(need(length(unique(ytest)) == 2,
                "The test split ended up with only one group present - lower the test-set size or provide more samples for this sex."))

  ## Class weights computed once on the full training split (folds recompute their own).
  cw_levels <- diag_class_weight_levels(y, params$class_weight_mode, params$class_weight_ratio)
  obs_w <- diag_obs_weights(y, params$class_weight_mode, params$class_weight_ratio)

  expr_train_sub <- expr_full[, split$train, drop = FALSE]
  expr_test_sub <- expr_full[, split$test, drop = FALSE]

  Ztr <- diag_zrows(expr_train_sub); rownames(Ztr) <- safe
  Xtr_full <- t(Ztr)                                    # z-scored, sample x gene(safe) - used for tuning + full fit
  Xraw <- t(expr_train_sub); colnames(Xraw) <- safe      # raw, sample x gene(safe) - re-scaled per fold inside diag_cv_auc

  ## Test split standardised with Train's own per-gene mean/SD (leakage-free).
  mu_tr <- colMeans(Xraw); sg_tr <- apply(Xraw, 2, stats::sd); sg_tr[is.na(sg_tr) | sg_tr == 0] <- 1
  Xtest_full <- scale(t(expr_test_sub), center = mu_tr, scale = sg_tr)
  colnames(Xtest_full) <- safe

  youden <- function(roc_obj) pROC::coords(roc_obj, "best", best.method = "youden",
                                            ret = c("threshold", "sensitivity", "specificity", "accuracy"), transpose = FALSE)

  ## Scores the Test split with the locked model, never re-optimising
  ## anything on it.
  score_eval <- function(pred_eval, y_eval, avail, reason) {
    if (!avail) return(list(available = FALSE, reason = reason))
    roc_e <- tryCatch(pROC::roc(y_eval, pred_eval, quiet = TRUE, levels = levels(y), direction = "<"), error = function(e) NULL)
    if (is.null(roc_e)) return(list(available = FALSE, reason = "ROC could not be computed for this split/contrast."))
    ci <- diag_auc_ci(roc_e)
    list(available = TRUE, roc = roc_e, auc = unname(ci["auc"]), ci_lo = unname(ci["lo"]), ci_hi = unname(ci["hi"]),
         n = length(y_eval), n_pos = sum(y_eval == levels(y)[2]))
  }

  ## (1) Elastic net: alpha tuned by minimising CV deviance/AUC/error over a
  ## grid, lambda by glmnet's own inner CV.
  nf_a <- max(2, min(params$enet_cv_folds, min(table(y))))
  type_measure <- params$enet_type_measure %||% "deviance"
  bigger_is_better <- identical(type_measure, "auc")
  best <- NULL; bcv <- if (bigger_is_better) -Inf else Inf
  ## Every alpha tried and its best score, for the hyperparameter search plot.
  alpha_search <- data.frame(alpha = numeric(0), cv_metric = numeric(0))
  set.seed(GLOBAL_SEED)
  for (a in params$enet_alpha_grid) {
    cv <- tryCatch(glmnet::cv.glmnet(Xtr_full, y, family = "binomial", alpha = a, nfolds = nf_a, standardize = TRUE,
                                      weights = obs_w, nlambda = params$enet_nlambda, type.measure = type_measure), error = function(e) NULL)
    if (!is.null(cv)) {
      m <- if (bigger_is_better) max(cv$cvm) else min(cv$cvm)
      alpha_search <- rbind(alpha_search, data.frame(alpha = a, cv_metric = m))
      if ((bigger_is_better && m > bcv) || (!bigger_is_better && m < bcv)) { bcv <- m; best <- list(cv = cv, alpha = a) }
    }
  }
  validate(need(!is.null(best), "Elastic net fitting failed for every alpha in the grid - check the gene panel and sample sizes."))
  alpha_search$chosen <- alpha_search$alpha == best$alpha
  lambda_s <- if (identical(params$enet_lambda_choice, "lambda.1se")) "lambda.1se" else "lambda.min"
  enet_pred_full <- as.numeric(predict(best$cv, newx = Xtr_full, s = lambda_s, type = "response"))
  enet_roc_full <- pROC::roc(y, enet_pred_full, quiet = TRUE, levels = levels(y), direction = "<")
  enet_best <- youden(enet_roc_full)
  enet_cv_auc <- diag_cv_auc(
    Xraw, y, params$enet_cv_folds,
    refit_fn = function(Ztr_i, ytr_i) glmnet::cv.glmnet(Ztr_i, ytr_i, family = "binomial", alpha = best$alpha,
                                                          nfolds = max(2, min(params$enet_cv_folds, min(table(ytr_i)))), standardize = TRUE,
                                                          weights = diag_obs_weights(ytr_i, params$class_weight_mode, params$class_weight_ratio),
                                                          nlambda = params$enet_nlambda, type.measure = type_measure),
    predict_fn = function(m, Zte) as.numeric(predict(m, newx = Zte, s = lambda_s, type = "response")),
    seed = GLOBAL_SEED
  )
  enet_pred_test <- as.numeric(predict(best$cv, newx = Xtest_full, s = lambda_s, type = "response"))
  enet_test <- score_eval(enet_pred_test, ytest, TRUE, NULL)
  if (isTRUE(enet_test$available)) enet_test$perf <- diag_perf_at_cutoff(enet_pred_test, ytest, enet_best$threshold, levels(y)[2])
  enet <- list(model = best$cv, model_type = "enet", label = "Elastic Net", alpha = best$alpha,
               lambda_choice = lambda_s, lambda_used = best$cv[[lambda_s]], tuning_search = alpha_search,
               type_measure = type_measure,
               pred_full = enet_pred_full, roc_full = enet_roc_full, full_auc = as.numeric(pROC::auc(enet_roc_full)),
               best = enet_best, cv_auc = enet_cv_auc, test = enet_test)

  ## (2) Random forest: mtry tuned by CV, ntree fixed, both user-overridable.
  p <- ncol(Xtr_full)
  ntree <- max(100, round(params$rf_ntree))
  rf_nodesize <- max(1, round(params$rf_nodesize %||% 1))
  rf_maxnodes <- if (!is.null(params$rf_maxnodes) && is.finite(params$rf_maxnodes)) max(2, round(params$rf_maxnodes)) else NULL
  mtry_search <- NULL
  if (identical(params$rf_mtry_mode, "manual") && !is.null(params$rf_mtry_manual)) {
    rf_mtry <- min(p, max(1, round(params$rf_mtry_manual)))
  } else {
    nf_rf <- max(2, min(params$rf_cv_folds, min(table(y))))
    mtry_grid <- sort(unique(pmin(p, c(1, 2, floor(sqrt(p)), floor(p / 3), floor(p / 2), p))))
    ctrl <- caret::trainControl(method = "cv", number = nf_rf, classProbs = TRUE, summaryFunction = caret::twoClassSummary)
    set.seed(GLOBAL_SEED)
    rf_tune <- tryCatch(caret::train(x = Xtr_full, y = y, method = "rf", metric = "ROC", trControl = ctrl,
                                      tuneGrid = expand.grid(mtry = mtry_grid), ntree = ntree,
                                      nodesize = rf_nodesize, maxnodes = rf_maxnodes, classwt = cw_levels), error = function(e) NULL)
    rf_mtry <- if (!is.null(rf_tune)) rf_tune$bestTune$mtry else max(1, floor(sqrt(p)))
    ## Every mtry tried and its CV ROC, for the hyperparameter search plot.
    if (!is.null(rf_tune)) { mtry_search <- rf_tune$results[, c("mtry", "ROC")]; mtry_search$chosen <- mtry_search$mtry == rf_mtry }
  }
  set.seed(GLOBAL_SEED)
  rf_model <- randomForest::randomForest(Xtr_full, y, ntree = ntree, mtry = rf_mtry,
                                          nodesize = rf_nodesize, maxnodes = rf_maxnodes, classwt = cw_levels)
  rf_pred_full <- predict(rf_model, Xtr_full, type = "prob")[, levels(y)[2]]
  rf_roc_full <- pROC::roc(y, rf_pred_full, quiet = TRUE, levels = levels(y), direction = "<")
  rf_best <- youden(rf_roc_full)
  rf_cv_auc <- diag_cv_auc(
    Xraw, y, params$rf_cv_folds,
    refit_fn = function(Ztr_i, ytr_i) randomForest::randomForest(Ztr_i, ytr_i, ntree = ntree, mtry = min(rf_mtry, ncol(Ztr_i)),
                                                                   nodesize = rf_nodesize, maxnodes = rf_maxnodes,
                                                                   classwt = diag_class_weight_levels(ytr_i, params$class_weight_mode, params$class_weight_ratio)),
    predict_fn = function(m, Zte) predict(m, Zte, type = "prob")[, levels(y)[2]],
    seed = GLOBAL_SEED
  )
  rf_pred_test <- predict(rf_model, Xtest_full, type = "prob")[, levels(y)[2]]
  rf_test <- score_eval(rf_pred_test, ytest, TRUE, NULL)
  if (isTRUE(rf_test$available)) rf_test$perf <- diag_perf_at_cutoff(rf_pred_test, ytest, rf_best$threshold, levels(y)[2])
  rf <- list(model = rf_model, model_type = "rf", label = "Random Forest", ntree = ntree, mtry = rf_mtry, tuning_search = mtry_search,
             pred_full = rf_pred_full, roc_full = rf_roc_full, full_auc = as.numeric(pROC::auc(rf_roc_full)),
             best = rf_best, cv_auc = rf_cv_auc, test = rf_test)

  ## (3) SVM: cost tuned by CV over a grid, kernel user-selectable (linear
  ## default). scale = FALSE since data is already gene-wise z-scored.
  kernel <- params$svm_kernel %||% "linear"
  svm_degree <- max(1, round(params$svm_degree %||% 3))
  svm_tolerance <- params$svm_tolerance %||% 0.001
  svm_gamma <- if (identical(params$svm_gamma_mode, "manual") && !is.null(params$svm_gamma_manual)) {
    params$svm_gamma_manual
  } else {
    1 / ncol(Xtr_full)
  }
  cost_search <- NULL
  if (identical(params$svm_cost_mode, "manual") && !is.null(params$svm_cost_manual)) {
    svm_cost <- params$svm_cost_manual
  } else {
    nf_svm <- max(2, min(params$svm_cv_folds, min(table(y))))
    grid <- params$svm_cost_grid
    if (!is.numeric(grid) || length(grid) == 0) grid <- DIAG_SVM_COST_GRID
    set.seed(GLOBAL_SEED)
    svm_tune <- tryCatch(e1071::tune(e1071::svm, train.x = Xtr_full, train.y = y, kernel = kernel, scale = FALSE,
                                      gamma = svm_gamma, degree = svm_degree, tolerance = svm_tolerance, class.weights = cw_levels,
                                      ranges = list(cost = grid),
                                      tunecontrol = e1071::tune.control(sampling = "cross", cross = nf_svm)), error = function(e) NULL)
    svm_cost <- if (!is.null(svm_tune)) svm_tune$best.parameters$cost else 1
    ## Every cost tried and its CV error, for the hyperparameter search plot.
    if (!is.null(svm_tune)) { cost_search <- svm_tune$performances[, c("cost", "error")]; cost_search$chosen <- cost_search$cost == svm_cost }
  }
  set.seed(GLOBAL_SEED)
  svm_model <- e1071::svm(Xtr_full, y, kernel = kernel, cost = svm_cost, scale = FALSE, probability = TRUE,
                           gamma = svm_gamma, degree = svm_degree, tolerance = svm_tolerance, class.weights = cw_levels)
  svm_pred_full <- attr(predict(svm_model, Xtr_full, probability = TRUE), "probabilities")[, levels(y)[2]]
  svm_roc_full <- pROC::roc(y, svm_pred_full, quiet = TRUE, levels = levels(y), direction = "<")
  svm_best <- youden(svm_roc_full)
  svm_cv_auc <- diag_cv_auc(
    Xraw, y, params$svm_cv_folds,
    refit_fn = function(Ztr_i, ytr_i) e1071::svm(Ztr_i, ytr_i, kernel = kernel, cost = svm_cost, scale = FALSE, probability = TRUE,
                                                   gamma = if (identical(params$svm_gamma_mode, "manual")) svm_gamma else 1 / ncol(Ztr_i),
                                                   degree = svm_degree, tolerance = svm_tolerance,
                                                   class.weights = diag_class_weight_levels(ytr_i, params$class_weight_mode, params$class_weight_ratio)),
    predict_fn = function(m, Zte) attr(predict(m, Zte, probability = TRUE), "probabilities")[, levels(y)[2]],
    seed = GLOBAL_SEED
  )
  svm_pred_test <- attr(predict(svm_model, Xtest_full, probability = TRUE), "probabilities")[, levels(y)[2]]
  svm_test <- score_eval(svm_pred_test, ytest, TRUE, NULL)
  if (isTRUE(svm_test$available)) svm_test$perf <- diag_perf_at_cutoff(svm_pred_test, ytest, svm_best$threshold, levels(y)[2])
  svm_fit <- list(model = svm_model, model_type = "svm", label = "SVM", kernel = kernel, cost = svm_cost, tuning_search = cost_search,
                   pred_full = svm_pred_full, roc_full = svm_roc_full, full_auc = as.numeric(pROC::auc(svm_roc_full)),
                   best = svm_best, cv_auc = svm_cv_auc, test = svm_test)

  ## (4) Logistic regression: plain, unpenalized glm on every gene in the
  ## panel. No hyperparameters, so tuning_search stays NULL.
  lr_predict <- function(m, Znew) as.numeric(predict(m, newdata = data.frame(Znew, check.names = FALSE), type = "response"))
  lr_model <- suppressWarnings(stats::glm(y ~ ., data = data.frame(y, Xtr_full, check.names = FALSE), family = stats::binomial, weights = obs_w))
  lr_pred_full <- as.numeric(predict(lr_model, type = "response"))
  lr_roc_full <- pROC::roc(y, lr_pred_full, quiet = TRUE, levels = levels(y), direction = "<")
  lr_best <- youden(lr_roc_full)
  lr_cv_auc <- diag_cv_auc(
    Xraw, y, params$lr_cv_folds,
    refit_fn = function(Ztr_i, ytr_i) suppressWarnings(stats::glm(ytr_i ~ ., data = data.frame(ytr_i, Ztr_i, check.names = FALSE), family = stats::binomial,
                                                                     weights = diag_obs_weights(ytr_i, params$class_weight_mode, params$class_weight_ratio))),
    predict_fn = lr_predict,
    seed = GLOBAL_SEED
  )
  lr_pred_test <- lr_predict(lr_model, Xtest_full)
  lr_test <- score_eval(lr_pred_test, ytest, TRUE, NULL)
  if (isTRUE(lr_test$available)) lr_test$perf <- diag_perf_at_cutoff(lr_pred_test, ytest, lr_best$threshold, levels(y)[2])
  lr <- list(model = lr_model, model_type = "lr", label = "Logistic Regression", tuning_search = NULL,
             pred_full = lr_pred_full, roc_full = lr_roc_full, full_auc = as.numeric(pROC::auc(lr_roc_full)),
             best = lr_best, cv_auc = lr_cv_auc, test = lr_test)

  list(lr = lr, enet = enet, rf = rf, svm = svm_fit, genes = genes, n_input = length(genes),
       n_samples = nrow(Xtr_full), n_test = nrow(Xtest_full), test_frac = params$test_frac,
       gene_roc_train = diag_gene_roc(expr_train_sub, y),
       gene_roc_test = diag_gene_roc(expr_test_sub, ytest))
}

DIAG_TECHNIQUES <- list(
  list(key = "lr", label = "Logistic Regression"),
  list(key = "enet", label = "Elastic Net"),
  list(key = "rf", label = "Random Forest"),
  list(key = "svm", label = "SVM")
)

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------

## Model Training tab box: KPI tiles plus ROC / CV-by-fold / tuning plots, full width.
## The ROC plot gets its own full-width row at a large height - coord_equal()
## needs real room to keep its title, axis labels and in-plot legend box
## legible instead of being squeezed into a one-third-width column.
mod_diagnostic_training_panel <- function(ns, prefix, title, roc_height = 440, side_height = 260) {
  box(
    width = NULL, title = title, status = "primary", solidHeader = FALSE,
    withSpinner(uiOutput(ns(paste0(prefix, "_train_stats"))), color = "#2563EB", type = 6),
    fluidRow(
      column(12, h5("ROC (Train vs Test)"), withSpinner(plotOutput(ns(paste0(prefix, "_train_roc_plot")), height = roc_height), color = "#2563EB", type = 6))
    ),
    fluidRow(
      column(6, h5("Cross-validated AUC by fold"), withSpinner(plotOutput(ns(paste0(prefix, "_train_cv_plot")), height = side_height), color = "#2563EB", type = 6)),
      column(6, h5("Hyperparameter tuning - explore the grid"), withSpinner(plotOutput(ns(paste0(prefix, "_train_tuning_plot")), height = side_height), color = "#2563EB", type = 6))
    ),
    div(class = "table-toolbar",
        downloadButton(ns(paste0(prefix, "_train_download")), "Performance (CSV)", class = "btn-sm"),
        downloadButton(ns(paste0(prefix, "_model_download")), "Model (.rds)", class = "btn-sm")),
    DT::dataTableOutput(ns(paste0(prefix, "_train_table")))
  )
}

## Model Testing (Internal) tab box: scores the held-out Test split for this sex.
mod_diagnostic_testing_panel <- function(ns, prefix, title, roc_height = 300) {
  box(
    width = NULL, title = title, status = "primary", solidHeader = FALSE,
    withSpinner(uiOutput(ns(paste0(prefix, "_test_summary"))), color = "#2563EB", type = 6),
    withSpinner(plotOutput(ns(paste0(prefix, "_test_roc_plot")), height = roc_height), color = "#2563EB", type = 6),
    div(class = "table-toolbar", downloadButton(ns(paste0(prefix, "_test_download")), "Performance (CSV)", class = "btn-sm")),
    DT::dataTableOutput(ns(paste0(prefix, "_test_table")))
  )
}

mod_diagnostic_params_box <- function(ns, prefix, method_label, defaults_desc, ...) {
  box(
    width = 12, title = sprintf("%s parameters", method_label), status = "primary", solidHeader = FALSE,
    p(class = "submodule-desc", defaults_desc),
    ...
  )
}

## One sex's single-gene ROC box, with adjustable hub-gene thresholds (AUC/P
## defaults follow Chen et al. 2021/2022's WGCNA panel rule).
mod_diagnostic_generoc_box_sex <- function(ns, sex_label, mode, title) {
  box(
    width = NULL, title = title, status = "primary", solidHeader = FALSE,
    withSpinner(plotOutput(ns(paste0(sex_label, "_", mode, "_generoc_plot")), height = "auto"), color = "#2563EB", type = 6),
    fluidRow(
      column(6, numericInput(ns(paste0(sex_label, "_", mode, "_hub_auc_thr")), "Hub gene rule: AUC ≥",
                              value = 0.85, min = 0.5, max = 1, step = 0.01)),
      column(6, numericInput(ns(paste0(sex_label, "_", mode, "_hub_p_thr")), "and P <",
                              value = 0.05, min = 0.0001, max = 1, step = 0.005))
    ),
    div(class = "table-toolbar",
        downloadButton(ns(paste0(sex_label, "_", mode, "_generoc_download")), "Gene AUCs (CSV)", class = "btn-sm"),
        downloadButton(ns(paste0(sex_label, "_", mode, "_hub_download")), "Hub genes only (CSV)", class = "btn-sm")),
    DT::dataTableOutput(ns(paste0(sex_label, "_", mode, "_generoc_table")))
  )
}

## One sex's full "Model Training" drill-down: Run button plus model pills,
## per-gene ROC and comparison table, hidden via conditionalPanel until run.
mod_diagnostic_training_sex_panel <- function(ns, sex_label) {
  run_id <- paste0("run_", sex_label, "_btn")
  cond <- sprintf("input['%s'] > 0 || input['%s'] > 0", ns(run_id), ns(paste0(run_id, "_test")))
  sex_title <- tools::toTitleCase(sex_label)
  tagList(
    actionButton(ns(run_id), paste("Run", sex_title), icon = icon("play"), class = "btn-primary btn-sm"),
    div(style = "height:10px;"),
    conditionalPanel(
      condition = cond,
      tabsetPanel(
        id = ns(paste0(sex_label, "_train_model_pills")), type = "pills",
        tabPanel("Logistic Regression", br(), mod_diagnostic_training_panel(ns, paste0(sex_label, "_lr"), NULL)),
        tabPanel("Elastic Net", br(), mod_diagnostic_training_panel(ns, paste0(sex_label, "_enet"), NULL)),
        tabPanel("Random Forest", br(), mod_diagnostic_training_panel(ns, paste0(sex_label, "_rf"), NULL)),
        tabPanel("SVM", br(), mod_diagnostic_training_panel(ns, paste0(sex_label, "_svm"), NULL))
      ),
      mod_diagnostic_generoc_box_sex(ns, sex_label, "train", sprintf("Per-gene ROC/AUC - Train vs Test (%s)", sex_title)),
      box(width = NULL, title = sprintf("Training comparison - %s", sex_title), status = "primary", solidHeader = FALSE,
          DT::dataTableOutput(ns(paste0(sex_label, "_train_compare_table")))),
      uiOutput(ns(paste0(sex_label, "_result_line")))
    )
  )
}

mod_diagnostic_testing_sex_panel <- function(ns, sex_label) {
  run_id <- paste0("run_", sex_label, "_btn_test")
  cond <- sprintf("input['%s'] > 0 || input['%s'] > 0", ns(paste0("run_", sex_label, "_btn")), ns(run_id))
  sex_title <- tools::toTitleCase(sex_label)
  tagList(
    actionButton(ns(run_id), paste("Run", sex_title), icon = icon("play"), class = "btn-primary btn-sm"),
    div(style = "height:10px;"),
    conditionalPanel(
      condition = cond,
      tabsetPanel(
        type = "pills",
        tabPanel("Logistic Regression", br(), mod_diagnostic_testing_panel(ns, paste0(sex_label, "_lr"), NULL)),
        tabPanel("Elastic Net", br(), mod_diagnostic_testing_panel(ns, paste0(sex_label, "_enet"), NULL)),
        tabPanel("Random Forest", br(), mod_diagnostic_testing_panel(ns, paste0(sex_label, "_rf"), NULL)),
        tabPanel("SVM", br(), mod_diagnostic_testing_panel(ns, paste0(sex_label, "_svm"), NULL))
      ),
      box(width = NULL, title = sprintf("Testing comparison - %s", sex_title), status = "primary", solidHeader = FALSE,
          DT::dataTableOutput(ns(paste0(sex_label, "_test_compare_table"))))
    )
  )
}

mod_diagnostic_ui <- function(id) {
  ns <- NS(id)
  tagList(
    ## Scoped to this module only (.diag-module) - a defensive floor/wrap on
    ## valueBox() headline numbers so a longer value (e.g. an AUC ± SD) never
    ## overflows its tile, whatever text ends up in one later.
    tags$style(HTML(".diag-module .small-box h3 { font-size: 22px; white-space: normal; overflow-wrap: break-word; line-height: 1.2; }")),
    div(class = "diag-module",
    fluidRow(
      column(
        3,
        box(
          width = NULL, title = "Gene panel & samples", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", "Pick a gene panel and contrast. Run each sex from its own tab under Model Training."),
          radioButtons(
            ns("panel_source"), NULL,
            choiceNames = list(
              tagList(icon("diagram-project"), " Follow this project's pipeline (recommended)"),
              tagList(icon("file-arrow-up"), " Paste my own gene list"),
              tagList(icon("circle-nodes"), " A WGCNA module from this session")
            ),
            choiceValues = list("project", "own", "wgcna"), selected = "project"
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'project'", ns("panel_source")),
            uiOutput(ns("project_source_ui"))
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'own'", ns("panel_source")),
            p(class = "submodule-desc", "Same list for Female, Male and Pooled."),
            textAreaInput(ns("gene_list"), NULL, rows = 5, placeholder = "TNF\nIL6\nSTAT3\n...")
          ),
          ## WGCNA module list, shared across Female/Male/Pooled, from mod_wgcna.R's Step 3.
          conditionalPanel(
            condition = sprintf("input['%s'] == 'wgcna'", ns("panel_source")),
            p(class = "submodule-desc", "Same module for Female, Male and Pooled."),
            uiOutput(ns("wgcna_module_pick_ui"))
          ),
          tags$hr(),
          uiOutput(ns("contrast_controls")),
          div(style = "margin-top:10px;", uiOutput(ns("saved_runs_ui")))
        )
      ),
      column(
        9,
        tabsetPanel(
          id = ns("main_tabs"), type = "tabs",
          tabPanel(
            "Model Training", br(),
            p(class = "submodule-desc", "Pick Female or Male, then Run - nothing below renders until then."),
            ## Params box tracks whichever model pill was most recently clicked
            ## in either sex's tab, via an "active pill" reactiveVal.
            uiOutput(ns("model_params_ui")),
            tabsetPanel(
              id = ns("train_sex_tabs"), type = "tabs",
              tabPanel("Female", br(), mod_diagnostic_training_sex_panel(ns, "female")),
              tabPanel("Male", br(), mod_diagnostic_training_sex_panel(ns, "male")),
              ## Pools all samples regardless of sex into one model.
              tabPanel("Pooled (all)", br(), mod_diagnostic_training_sex_panel(ns, "pooled"))
            )
          ),
          tabPanel(
            "Model Testing (Internal)", br(),
            p(class = "submodule-desc", "Each Train-fit model scored once on its held-out Test split."),
            tabsetPanel(
              id = ns("test_sex_tabs"), type = "tabs",
              tabPanel("Female", br(), mod_diagnostic_testing_sex_panel(ns, "female")),
              tabPanel("Male", br(), mod_diagnostic_testing_sex_panel(ns, "male")),
              tabPanel("Pooled (all)", br(), mod_diagnostic_testing_sex_panel(ns, "pooled"))
            )
          ),
          ## A genuinely separate uploaded cohort (Chen et al. 2021/2022's
          ## GSE77298 role): per-gene AUC/P validation, not a refit model.
          tabPanel(
            "External Validation", br(),
            mod_diagnostic_external_panel(ns)
          )
        )
      )
    ),
    uiOutput(ns("references_box_ui"))
    )
  )
}

## External-validation tab: separate file upload, its own group-column mapping,
## and a panel selector reusing the same gene panel set up in the sidebar.
mod_diagnostic_external_panel <- function(ns) {
  tagList(
    box(
      width = NULL, title = "External validation dataset", status = "primary", solidHeader = FALSE,
      p(class = "submodule-desc", "Upload a separate cohort to check whether the gene panel below holds up outside the training data. This dataset is used for validation only - it is never used to train or refit any model."),
      fluidRow(
        column(6,
          fileInput(ns("ext_expr_file"), "External validation expression matrix", accept = c(".csv", ".rds", ".Rds")),
          div(class = "empty-note", style = "font-size: 12.5px; margin-top: -8px;", icon("circle-info"),
              "CSV or RDS. Genes in rows, samples in columns; for CSV, the first column is the gene ID.")
        ),
        column(6, fileInput(ns("ext_meta_file"), "External validation sample metadata", accept = c(".csv", ".rds", ".Rds")))
      ),
      uiOutput(ns("ext_column_mapping")),
      selectInput(ns("ext_panel_choice"), "Gene panel to validate (uses the panel source set on the left)",
                  choices = c("Pooled" = "pooled", "Female" = "female", "Male" = "male"), selected = "pooled", selectize = FALSE),
      fluidRow(
        column(6, numericInput(ns("ext_hub_auc_thr"), "Hub gene rule: AUC ≥", value = 0.85, min = 0.5, max = 1, step = 0.01)),
        column(6, numericInput(ns("ext_hub_p_thr"), "and P <", value = 0.05, min = 0.0001, max = 1, step = 0.005))
      ),
      actionButton(ns("run_ext_btn"), "Run external validation", icon = icon("play"), class = "btn-primary btn-sm")
    ),
    uiOutput(ns("ext_status_ui")),
    withSpinner(uiOutput(ns("ext_results_ui")), color = "#2563EB", type = 6)
  )
}

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------

mod_diagnostic_server <- function(id, dataset, results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Detects which value in the loaded sex column means "female"/"male".
    sex_levels <- reactive({
      lv <- unique(stats::na.omit(as.character(dataset$meta$sex)))
      validate(need(length(lv) >= 2, "The loaded metadata needs a \"sex\" column with at least two distinct values."))
      f <- lv[grepl("^f", lv, ignore.case = TRUE)]
      m <- lv[grepl("^m", lv, ignore.case = TRUE)]
      lv_sorted <- sort(lv)
      list(female = if (length(f) > 0) f[1] else lv_sorted[1],
           male   = if (length(m) > 0) m[1] else lv_sorted[min(2, length(lv_sorted))])
    })

    output$contrast_controls <- renderUI({
      groups <- sort(unique(stats::na.omit(dataset$meta$group)))
      validate(need(length(groups) >= 2, "The loaded metadata needs at least two group values."))
      tagList(
        selectInput(ns("ref_group"), "Reference group (negative class)", choices = groups, selected = groups[1], selectize = FALSE),
        selectInput(ns("comp_group"), "Comparison group (positive class)", choices = groups, selected = groups[min(2, length(groups))], selectize = FALSE),
        sliderInput(ns("test_frac_pct"), "Train : Test ratio (default 70:30)", min = 10, max = 50, value = 30, step = 5, post = "% test"),
        div(style = "font-size: 12.5px; color: #64748B; margin-top: -8px; margin-bottom: 8px;", strong(textOutput(ns("ratio_caption"), inline = TRUE))),
        div(class = "empty-note", style = "font-size: 12.5px; margin-top: -2px;", icon("circle-info"),
            "Same split used for all four models per sex."),
        radioButtons(ns("class_weight_mode"), "Class weighting (imbalanced groups)",
                     choices = c("Equal - this project's own methodology (default)" = "equal",
                                 "Balanced - auto inverse-frequency" = "balanced",
                                 "Manual ratio" = "manual"),
                     selected = "equal"),
        conditionalPanel(condition = sprintf("input['%s'] == 'manual'", ns("class_weight_mode")),
                          numericInput(ns("class_weight_ratio"), "Weight ratio (comparison : reference)", value = 1, min = 0.05, max = 20, step = 0.05)),
        div(class = "empty-note", style = "font-size: 12.5px; margin-top: -6px;", icon("circle-info"),
            "Applied to all four models, both sexes.")
      )
    })

    ## Gene panel sources. The bundled/precomputed panel was only ever computed
    ## from the app's own default merged cohort - only offer it as a fallback
    ## when that exact dataset is still active (dataset$is_bundled_reference),
    ## never for an uploaded, GEO-fetched, or individual raw preloaded dataset,
    ## where it would silently mix results from a different dataset entirely.
    project_panel_genes <- function(sex_label) {
      live <- results$featureselection[[sex_label]]$consensus_genes
      if (!is.null(live) && length(live) >= 2) {
        return(list(genes = live, is_live = TRUE,
                    note = sprintf("%d genes from this session's live %s consensus panel.", length(live), sex_label)))
      }
      if (isTRUE(dataset$is_bundled_reference)) {
        bundled <- read_table_safe(sprintf("FS_input_%s.csv", sex_label))
        if (!is.null(bundled) && nrow(bundled) >= 2 && "gene" %in% colnames(bundled)) {
          return(list(genes = unique(as.character(bundled$gene)), is_live = FALSE,
                      note = sprintf("%d genes from the bundled %s panel.", nrow(bundled), sex_label)))
        }
      }
      list(genes = character(0), is_live = FALSE,
           note = sprintf("No live %s panel yet - run Feature Selection on the currently loaded dataset first.", sex_label))
    }

    own_panel_genes <- function(sex_label) {
      genes <- unique(trimws(unlist(strsplit(input$gene_list %||% "", "[,\n\t ]+"))))
      genes <- genes[nzchar(genes)]
      list(genes = genes, is_live = FALSE, note = sprintf("%d pasted genes.", length(genes)))
    }

    ## Same WGCNA module regardless of sex_label (not sex-specific by construction).
    wgcna_panel_genes <- function(sex_label) {
      req(input$wgcna_module_pick)
      mg <- results$wgcna$module_genes
      validate(need(!is.null(mg) && input$wgcna_module_pick %in% names(mg),
                    "Run WGCNA Step 3 (Modules) first, then pick a module above."))
      genes <- unique(as.character(mg[[input$wgcna_module_pick]]))
      list(genes = genes, is_live = TRUE,
           note = sprintf("%d genes from WGCNA module \"%s\" (this session).", length(genes), input$wgcna_module_pick))
    }

    ## Module-color dropdown (grey excluded), mirroring mod_featureselection.R's own picker.
    output$wgcna_module_pick_ui <- renderUI({
      mg <- results$wgcna$module_genes
      if (is.null(mg) || length(mg) == 0) {
        return(div(class = "empty-note", icon("circle-info"),
                    "No WGCNA modules yet - run Step 3 (Modules) in the WGCNA Co-expression Network tab first."))
      }
      sizes <- vapply(mg, length, integer(1))
      choices <- setNames(names(mg), sprintf("%s (%d genes)", names(mg), sizes))
      choices <- choices[names(mg) != "grey"]
      tagList(
        selectInput(ns("wgcna_module_pick"), "Module", choices = choices, selectize = FALSE),
        div(class = "empty-note", icon("circle-info"),
            "Uses expression values from the Dataset tab - use the same dataset WGCNA was run on.")
      )
    })

    output$project_source_ui <- renderUI({
      f_live <- results$featureselection$female$consensus_genes
      m_live <- results$featureselection$male$consensus_genes
      p_live <- results$featureselection$pooled$consensus_genes
      bits <- c(
        if (length(f_live) >= 2) sprintf("%d female", length(f_live)),
        if (length(m_live) >= 2) sprintf("%d male", length(m_live)),
        if (length(p_live) >= 2) sprintf("%d pooled", length(p_live))
      )
      if (length(bits) > 0) {
        div(class = "empty-note", icon("check"),
            sprintf("Live Feature Selection panel: %s genes.", paste(bits, collapse = " / ")))
      } else if (isTRUE(dataset$is_bundled_reference)) {
        div(class = "empty-note", icon("circle-info"),
            "No live Feature Selection panel yet - using the bundled female/male panel. Pooled needs a live run first.")
      } else {
        div(class = "empty-note", icon("triangle-exclamation"),
            "No live Feature Selection panel yet for this dataset - the bundled panel only applies to the app's default reference cohort. Run Feature Selection on the currently loaded dataset first.")
      }
    })

    ## External Validation: a separately uploaded dataset, scored with
    ## diag_gene_roc()'s per-gene AUC/Wilcoxon-P (no multivariate refit).
    ext_meta_raw <- reactive({
      req(input$ext_meta_file)
      path <- input$ext_meta_file$datapath
      if (grepl("\\.rds$", input$ext_meta_file$name, ignore.case = TRUE)) {
        d <- readRDS(path)
        validate(need(is.data.frame(d), "The uploaded metadata RDS file must contain a data frame."))
        as.data.frame(d)
      } else {
        as.data.frame(data.table::fread(path, showProgress = FALSE))
      }
    })

    output$ext_column_mapping <- renderUI({
      req(input$ext_expr_file, ext_meta_raw())
      cols <- colnames(ext_meta_raw())
      tagList(
        fluidRow(
          column(6, selectInput(ns("ext_map_id"), "Sample ID column", choices = cols, selected = cols[1], selectize = FALSE)),
          column(6, selectInput(ns("ext_map_group"), "Group / diagnosis column", choices = cols, selectize = FALSE))
        ),
        uiOutput(ns("ext_group_pick_ui"))
      )
    })

    ext_data <- reactive({
      req(input$ext_expr_file, input$ext_meta_file, input$ext_map_id, input$ext_map_group)
      expr <- if (grepl("\\.rds$", input$ext_expr_file$name, ignore.case = TRUE)) {
        res <- tx_parse_expr_matrix_rds(input$ext_expr_file$datapath)
        validate(need(res$ok, res$error))
        res$mat
      } else {
        m <- as.data.frame(data.table::fread(input$ext_expr_file$datapath, showProgress = FALSE))
        rn <- as.character(m[[1]])
        m <- as.matrix(m[, -1, drop = FALSE])
        rownames(m) <- rn
        m
      }
      meta <- ext_meta_raw()
      meta$sample <- as.character(meta[[input$ext_map_id]])
      meta$group <- as.character(meta[[input$ext_map_group]])
      common <- intersect(colnames(expr), meta$sample)
      validate(need(length(common) >= 6, "Fewer than 6 sample IDs in the external expression matrix match the metadata sample-ID column. Check the column mapping."))
      list(expr = expr[, common, drop = FALSE], meta = meta[match(common, meta$sample), , drop = FALSE])
    })

    output$ext_group_pick_ui <- renderUI({
      groups <- sort(unique(stats::na.omit(ext_data()$meta$group)))
      validate(need(length(groups) >= 2, "The external metadata's group column needs at least two distinct values."))
      fluidRow(
        column(6, selectInput(ns("ext_ref_group"), "Reference group (negative class)", choices = groups, selected = groups[1], selectize = FALSE)),
        column(6, selectInput(ns("ext_comp_group"), "Comparison group (positive class)", choices = groups, selected = groups[min(2, length(groups))], selectize = FALSE))
      )
    })

    ext_result <- eventReactive(input$run_ext_btn, {
      d <- ext_data()
      req(input$ext_ref_group, input$ext_comp_group)
      validate(need(input$ext_ref_group != input$ext_comp_group, "Reference and comparison group must be different."))

      panel_sex <- input$ext_panel_choice %||% "pooled"
      cand <- switch(input$panel_source,
        project = project_panel_genes(panel_sex),
        wgcna = wgcna_panel_genes(panel_sex),
        own_panel_genes(panel_sex)
      )
      genes <- intersect(cand$genes, rownames(d$expr))
      validate(need(length(genes) >= 1, "None of the chosen gene panel's genes are present in the uploaded external dataset."))

      keep <- as.character(d$meta$group) %in% c(input$ext_ref_group, input$ext_comp_group)
      meta_sub <- d$meta[keep, , drop = FALSE]
      y <- factor(as.character(meta_sub$group), levels = c(input$ext_ref_group, input$ext_comp_group))
      validate(need(all(table(y) >= 3), "Each group needs at least 3 samples in the external dataset."))
      expr_sub <- d$expr[genes, meta_sub$sample, drop = FALSE]

      gr <- diag_gene_roc(expr_sub, y)
      list(gr = gr, expr_sub = expr_sub, y = y, genes = genes, panel_note = cand$note,
           n_ref = sum(y == input$ext_ref_group), n_comp = sum(y == input$ext_comp_group),
           ref_group = input$ext_ref_group, comp_group = input$ext_comp_group)
    }, ignoreInit = TRUE)

    ## Same never-run vs. actually-failed distinction as diag_result_error_msg() above,
    ## for the External Validation upload pipeline - a separate, easy-to-miss failure
    ## mode here is the chosen gene panel being empty (e.g. no live Feature Selection
    ## panel yet for an uploaded reference dataset), which used to show the exact same
    ## "Not run yet" text as never having clicked Run at all.
    ext_result_error_msg <- function() {
      tryCatch({ ext_result(); NULL }, error = function(e) {
        msg <- conditionMessage(e)
        if (nzchar(msg)) msg else NULL
      })
    }
    observeEvent(input$run_ext_btn, {
      err <- ext_result_error_msg()
      if (!is.null(err)) showNotification(paste("External validation could not be completed:", err), type = "error", duration = 10)
    }, ignoreInit = TRUE)

    output$ext_status_ui <- renderUI({
      r <- tryCatch(ext_result(), error = function(e) NULL)
      if (!is.null(r)) {
        return(div(class = "empty-note", icon("check"),
            sprintf("%d panel genes present in the external dataset, %d samples (%d %s vs %d %s). %s",
                    length(r$genes), length(r$y), r$n_comp, r$comp_group, r$n_ref, r$ref_group, r$panel_note)))
      }
      err <- ext_result_error_msg()
      if (!is.null(err)) {
        div(class = "empty-note", icon("triangle-exclamation"), sprintf("External validation failed: %s", err))
      } else {
        div(class = "empty-note", icon("circle-info"), "Not run yet. Upload files, map columns, then click Run.")
      }
    })

    ext_gene_df <- reactive({
      r <- ext_result()
      p <- if (!is.null(r$gr$p)) unname(r$gr$p) else rep(NA_real_, length(r$gr$genes))
      df <- data.frame(gene = r$gr$genes, auc = round(unname(r$gr$auc), 3), p = signif(p, 3), stringsAsFactors = FALSE)
      auc_thr <- input$ext_hub_auc_thr %||% 0.85
      p_thr <- input$ext_hub_p_thr %||% 0.05
      df$hub <- !is.na(df$auc) & !is.na(df$p) & df$auc >= auc_thr & df$p < p_thr
      df[order(-df$auc), ]
    })

    output$ext_gene_table <- DT::renderDataTable({
      req(ext_result())
      DT::datatable(ext_gene_df(), rownames = FALSE, width = "100%",
                    options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact") |>
        DT::formatStyle("hub", target = "row", backgroundColor = DT::styleEqual(c(TRUE, FALSE), c("#e6f4ea", "")))
    })

    output$ext_gene_download <- downloadHandler(
      filename = function() "external_validation_gene_auc.csv",
      content = function(file) write.csv(ext_gene_df(), file, row.names = FALSE)
    )

    ## Per-gene expression boxplot (Fig. 3d style), capped to top 24 genes by AUC.
    EXT_BOXPLOT_MAX_GENES <- 24
    output$ext_boxplot <- renderPlot({
      r <- ext_result(); req(r)
      top_genes <- head(ext_gene_df()$gene, EXT_BOXPLOT_MAX_GENES)
      df <- do.call(rbind, lapply(top_genes, function(g) {
        data.frame(gene = g, expr = as.numeric(r$expr_sub[g, ]), group = as.character(r$y), stringsAsFactors = FALSE)
      }))
      ggplot(df, aes(x = group, y = expr, fill = group)) +
        geom_boxplot(outlier.size = 0.6) +
        facet_wrap(~gene, scales = "free_y") +
        scale_fill_manual(values = arthomix_pair(factor(df$group)), guide = "none") +
        labs(x = NULL, y = "Expression (external dataset)") +
        theme_arthomix(base_size = 11)
    })

    output$ext_results_ui <- renderUI({
      r <- tryCatch(ext_result(), error = function(e) NULL)
      if (is.null(r)) return(NULL)
      tagList(
        box(
          width = NULL, title = "Per-gene ROC/AUC - External validation", status = "primary", solidHeader = FALSE,
          div(class = "table-toolbar", downloadButton(ns("ext_gene_download"), "Gene AUCs (CSV)", class = "btn-sm")),
          DT::dataTableOutput(ns("ext_gene_table"))
        ),
        box(
          width = NULL, title = "Expression by group - External validation", status = "primary", solidHeader = FALSE,
          p(class = "submodule-desc", sprintf("Top %d genes by AUC shown; download the CSV above for the full panel.", EXT_BOXPLOT_MAX_GENES)),
          plotOutput(ns("ext_boxplot"), height = 500)
        )
      )
    })

    ## Reads per-model advanced parameters from the UI, falling back to defaults.
    diag_advanced_params <- function() {
      alpha_grid <- suppressWarnings(as.numeric(trimws(strsplit(input$enet_alpha_grid %||% "", ",")[[1]])))
      alpha_grid <- alpha_grid[!is.na(alpha_grid) & alpha_grid >= 0 & alpha_grid <= 1]
      cost_grid <- suppressWarnings(as.numeric(trimws(strsplit(input$svm_cost_grid %||% "", ",")[[1]])))
      cost_grid <- cost_grid[!is.na(cost_grid) & cost_grid > 0]
      list(
        test_frac = (input$test_frac_pct %||% (DIAG_DEFAULT_PARAMS$test_frac * 100)) / 100,
        class_weight_mode = input$class_weight_mode %||% DIAG_DEFAULT_PARAMS$class_weight_mode,
        class_weight_ratio = input$class_weight_ratio %||% DIAG_DEFAULT_PARAMS$class_weight_ratio,
        lr_cv_folds = input$lr_cv_folds %||% DIAG_DEFAULT_PARAMS$lr_cv_folds,
        enet_cv_folds = input$enet_cv_folds %||% DIAG_DEFAULT_PARAMS$enet_cv_folds,
        enet_alpha_grid = if (length(alpha_grid) > 0) alpha_grid else DIAG_DEFAULT_PARAMS$enet_alpha_grid,
        enet_lambda_choice = input$enet_lambda_choice %||% DIAG_DEFAULT_PARAMS$enet_lambda_choice,
        enet_nlambda = input$enet_nlambda %||% DIAG_DEFAULT_PARAMS$enet_nlambda,
        enet_type_measure = input$enet_type_measure %||% DIAG_DEFAULT_PARAMS$enet_type_measure,
        rf_cv_folds = input$rf_cv_folds %||% DIAG_DEFAULT_PARAMS$rf_cv_folds,
        rf_ntree = input$rf_ntree %||% DIAG_DEFAULT_PARAMS$rf_ntree,
        rf_mtry_mode = input$rf_mtry_mode %||% DIAG_DEFAULT_PARAMS$rf_mtry_mode,
        rf_mtry_manual = input$rf_mtry_manual,
        rf_nodesize = input$rf_nodesize %||% DIAG_DEFAULT_PARAMS$rf_nodesize,
        rf_maxnodes = if (isTRUE(input$rf_maxnodes_unlimited)) NULL else (input$rf_maxnodes %||% DIAG_DEFAULT_PARAMS$rf_maxnodes),
        svm_cv_folds = input$svm_cv_folds %||% DIAG_DEFAULT_PARAMS$svm_cv_folds,
        svm_kernel = input$svm_kernel %||% DIAG_DEFAULT_PARAMS$svm_kernel,
        svm_cost_mode = input$svm_cost_mode %||% DIAG_DEFAULT_PARAMS$svm_cost_mode,
        svm_cost_manual = input$svm_cost_manual %||% DIAG_DEFAULT_PARAMS$svm_cost_manual,
        svm_cost_grid = if (length(cost_grid) > 0) cost_grid else DIAG_SVM_COST_GRID,
        svm_gamma_mode = input$svm_gamma_mode %||% DIAG_DEFAULT_PARAMS$svm_gamma_mode,
        svm_gamma_manual = input$svm_gamma_manual %||% DIAG_DEFAULT_PARAMS$svm_gamma_manual,
        svm_degree = input$svm_degree %||% DIAG_DEFAULT_PARAMS$svm_degree,
        svm_tolerance = input$svm_tolerance %||% DIAG_DEFAULT_PARAMS$svm_tolerance
      )
    }

    output$ratio_caption <- renderText({
      t <- input$test_frac_pct %||% 30
      sprintf("Currently %d : %d (train : test)", 100 - t, t)
    })

    ## Builds and fits all four models for one sex; sex_value = NULL means pooled.
    diag_build_sex <- function(sex_label, sex_value) {
      req(input$ref_group, input$comp_group)
      validate(need(input$ref_group != input$comp_group, "Reference and comparison group must be different."))

      meta <- dataset$meta
      sex_ok <- if (is.null(sex_value)) rep(TRUE, nrow(meta)) else (!is.na(meta$sex) & as.character(meta$sex) == sex_value)
      meta <- meta[sex_ok &
                     !is.na(meta$group) & as.character(meta$group) %in% c(input$ref_group, input$comp_group), , drop = FALSE]
      common <- intersect(colnames(dataset$expr), meta$sample)
      validate(need(length(common) >= 10, sprintf("Fewer than 10 %s samples match this contrast.", sex_label)))
      meta <- meta[match(common, meta$sample), , drop = FALSE]
      y <- factor(as.character(meta$group), levels = c(input$ref_group, input$comp_group))
      validate(need(all(table(y) >= 6), sprintf("Each group needs at least 6 %s samples.", sex_label)))

      cand <- switch(input$panel_source,
        project = project_panel_genes(sex_label),
        wgcna = wgcna_panel_genes(sex_label),
        own_panel_genes(sex_label)
      )
      ## Checked separately from the "present in the matrix" validate() below - an
      ## empty panel (most commonly: an uploaded/GEO dataset with no live Feature
      ## Selection run yet, since the bundled panel only applies to the default
      ## reference cohort) needs its own message. Both used to collapse into the
      ## same generic "fewer than 2 genes ... present in the matrix" text, which
      ## reads as a data problem when the real issue is "run Feature Selection
      ## first" - cand$note already has that exact, specific explanation.
      validate(need(length(cand$genes) >= 1, sprintf("No gene panel available for %s: %s", sex_label, cand$note)))
      genes <- intersect(cand$genes, rownames(dataset$expr))
      validate(need(length(genes) >= 2, sprintf(
        "Fewer than 2 %s genes from the chosen panel are present in the currently loaded expression matrix (panel: %s).",
        sex_label, cand$note
      )))

      expr_sub <- dataset$expr[genes, common, drop = FALSE]

      ## withProgress (same mechanism already used in mod_wgcna.R/mod_featureselection.R)
      ## - fitting 4 CV-tuned models per sex is genuinely slow, and previously gave no
      ## visible feedback at all while running.
      fit <- withProgress(
        message = sprintf("Fitting %s diagnostic models (logistic regression, elastic net, random forest, SVM)...", sex_label),
        value = 0.3,
        diag_fit_sex(expr_sub, y, params = diag_advanced_params())
      )
      fit$candidate_note <- cand$note
      fit$n_ref <- sum(y == input$ref_group); fit$n_comp <- sum(y == input$comp_group)
      fit$ref_group <- input$ref_group; fit$comp_group <- input$comp_group
      fit
    }

    ## Training and Testing tabs each have their own Run button per sex; both feed
    ## the same shared trigger below.
    ## An eventReactive can't be reset to its unfired state, so this tracks whether
    ## each sex's cached diag_result_*() still belongs to the active dataset; cleared
    ## when dataset$source changes, set again when that sex's Run button is clicked.
    diag_valid <- reactiveValues(female = FALSE, male = FALSE, pooled = FALSE)

    female_run_trigger <- reactiveVal(0)
    male_run_trigger <- reactiveVal(0)
    lapply(c("run_female_btn", "run_female_btn_test"), function(bid) {
      observeEvent(input[[bid]], { diag_valid$female <- TRUE; female_run_trigger(isolate(female_run_trigger()) + 1) }, ignoreInit = TRUE)
    })
    lapply(c("run_male_btn", "run_male_btn_test"), function(bid) {
      observeEvent(input[[bid]], { diag_valid$male <- TRUE; male_run_trigger(isolate(male_run_trigger()) + 1) }, ignoreInit = TRUE)
    })
    pooled_run_trigger <- reactiveVal(0)
    lapply(c("run_pooled_btn", "run_pooled_btn_test"), function(bid) {
      observeEvent(input[[bid]], { diag_valid$pooled <- TRUE; pooled_run_trigger(isolate(pooled_run_trigger()) + 1) }, ignoreInit = TRUE)
    })

    diag_result_female <- eventReactive(female_run_trigger(), {
      diag_build_sex("female", sex_levels()$female)
    }, ignoreInit = TRUE)
    diag_result_male <- eventReactive(male_run_trigger(), {
      diag_build_sex("male", sex_levels()$male)
    }, ignoreInit = TRUE)
    diag_result_pooled <- eventReactive(pooled_run_trigger(), {
      diag_build_sex("pooled", NULL)
    }, ignoreInit = TRUE)

    ## Reads one sex's diag_result_*() and returns its real validate()/need() failure
    ## message, or NULL if it hasn't failed. Distinguishes a genuine failure (e.g. no
    ## live gene panel for an uploaded dataset, too few samples in a group) from
    ## "hasn't been run yet": both are eventReactive halts of the same shiny
    ## "validation" class, but eventReactive's own pre-first-click halt always carries
    ## an EMPTY message, while a validate(need(...)) failure inside diag_build_sex()
    ## always carries the real one - so nzchar() on the caught message tells them
    ## apart. Every existing read of diag_result_*() elsewhere in this file
    ## (`tryCatch(..., error = function(e) NULL)`) collapsed both cases to NULL,
    ## which made clicking Run on an uploaded dataset with no live Feature Selection
    ## panel yet look exactly like the button doing nothing - no error, no result,
    ## just a generic "not run yet" note with no explanation.
    diag_result_error_msg <- function(sex_label) {
      if (!isTRUE(diag_valid[[sex_label]])) return(NULL)
      fr <- switch(sex_label, female = diag_result_female, male = diag_result_male, pooled = diag_result_pooled)
      tryCatch({ fr(); NULL }, error = function(e) {
        msg <- conditionMessage(e)
        if (nzchar(msg)) msg else NULL
      })
    }

    ## Single read path for a sex's cached result: NULL both when it never ran and
    ## when the cached fit belongs to a dataset the user has since switched away from.
    diag_result_value <- function(sex_label) {
      if (!isTRUE(diag_valid[[sex_label]])) return(NULL)
      fr <- switch(sex_label, female = diag_result_female, male = diag_result_male, pooled = diag_result_pooled)
      tryCatch(fr(), error = function(e) NULL)
    }
    lapply(c("female", "male", "pooled"), function(sex_label) {
      trigger <- switch(sex_label, female = female_run_trigger, male = male_run_trigger, pooled = pooled_run_trigger)
      observeEvent(trigger(), {
        err <- diag_result_error_msg(sex_label)
        if (!is.null(err)) {
          showNotification(sprintf("%s diagnostic models could not be completed: %s", tools::toTitleCase(sex_label), err),
                            type = "error", duration = 10)
        }
      }, ignoreInit = TRUE)
    })

    diag_has_run <- reactiveVal(FALSE)
    observeEvent(female_run_trigger(), diag_has_run(TRUE), ignoreInit = TRUE)
    observeEvent(male_run_trigger(), diag_has_run(TRUE), ignoreInit = TRUE)
    observeEvent(pooled_run_trigger(), diag_has_run(TRUE), ignoreInit = TRUE)

    observeEvent(dataset$source, {
      diag_valid$female <- FALSE; diag_valid$male <- FALSE; diag_valid$pooled <- FALSE
      diag_has_run(FALSE)
    }, ignoreInit = TRUE)

    ## Which model's params box to show: whichever pill was clicked most
    ## recently, in any sex's own Model Training tab.
    active_model_pill <- reactiveVal("Logistic Regression")
    observeEvent(input$female_train_model_pills, active_model_pill(input$female_train_model_pills), ignoreInit = TRUE)
    observeEvent(input$male_train_model_pills, active_model_pill(input$male_train_model_pills), ignoreInit = TRUE)
    observeEvent(input$pooled_train_model_pills, active_model_pill(input$pooled_train_model_pills), ignoreInit = TRUE)

    lr_params_box <- function() {
      mod_diagnostic_params_box(
        ns, "lr", "Logistic Regression",
        "Plain, unpenalized logistic regression on the full panel - nothing to tune except CV folds.",
        fluidRow(
          column(4, numericInput(ns("lr_cv_folds"), "Cross-validation folds", value = 5, min = 3, max = 10, step = 1))
        )
      )
    }

    enet_params_box <- function() {
      mod_diagnostic_params_box(
        ns, "enet", "Elastic Net",
        "Alpha tuned by CV deviance; lambda.min by default.",
        fluidRow(
          column(4, numericInput(ns("enet_cv_folds"), "Cross-validation folds", value = 5, min = 3, max = 10, step = 1)),
          column(4, textInput(ns("enet_alpha_grid"), "Alpha grid (0 = ridge … 1 = LASSO, comma-separated)", value = paste(DIAG_DEFAULT_PARAMS$enet_alpha_grid, collapse = ", "))),
          column(4, radioButtons(ns("enet_lambda_choice"), "Lambda", choices = c("lambda.min (default)" = "lambda.min", "lambda.1se (sparser)" = "lambda.1se"), selected = "lambda.min"))
        ),
        h5("Advanced", style = "margin-top: 6px;"),
        fluidRow(
          column(4, numericInput(ns("enet_nlambda"), "Number of lambda values searched (nlambda)", value = DIAG_DEFAULT_PARAMS$enet_nlambda, min = 20, max = 300, step = 10)),
          column(4, radioButtons(ns("enet_type_measure"), "CV metric to optimize",
                                  choices = c("Deviance (default)" = "deviance", "AUC" = "auc", "Misclassification error" = "class"),
                                  selected = "deviance"))
        )
      )
    }

    rf_params_box <- function() {
      mod_diagnostic_params_box(
        ns, "rf", "Random Forest",
        "mtry tuned by CV; ntree fixed.",
        fluidRow(
          column(4, numericInput(ns("rf_cv_folds"), "Cross-validation folds", value = 5, min = 3, max = 10, step = 1)),
          column(4, numericInput(ns("rf_ntree"), "Number of trees", value = 1000, min = 100, max = 5000, step = 100)),
          column(4,
            radioButtons(ns("rf_mtry_mode"), "mtry (per split)", choices = c("Auto-tune (default)" = "auto", "Manual" = "manual"), selected = "auto"),
            conditionalPanel(condition = sprintf("input['%s'] == 'manual'", ns("rf_mtry_mode")),
                              numericInput(ns("rf_mtry_manual"), "mtry value", value = 5, min = 1, max = 500, step = 1))
          )
        ),
        h5("Advanced", style = "margin-top: 6px;"),
        fluidRow(
          column(4, numericInput(ns("rf_nodesize"), "Minimum terminal node size (nodesize)", value = DIAG_DEFAULT_PARAMS$rf_nodesize, min = 1, max = 50, step = 1)),
          column(4,
            checkboxInput(ns("rf_maxnodes_unlimited"), "Unlimited tree depth (maxnodes, default)", value = TRUE),
            conditionalPanel(condition = sprintf("!input['%s']", ns("rf_maxnodes_unlimited")),
                              numericInput(ns("rf_maxnodes"), "Max terminal nodes per tree", value = 20, min = 2, max = 2000, step = 1))
          )
        )
      )
    }

    svm_params_box <- function() {
      mod_diagnostic_params_box(
        ns, "svm", "SVM",
        "Cost tuned by CV; linear kernel by default.",
        fluidRow(
          column(3, numericInput(ns("svm_cv_folds"), "Cross-validation folds", value = 5, min = 3, max = 10, step = 1)),
          column(3, radioButtons(ns("svm_kernel"), "Kernel", choices = c("Linear (default)" = "linear", "Radial" = "radial", "Polynomial" = "polynomial"), selected = "linear")),
          column(6,
            radioButtons(ns("svm_cost_mode"), "Cost (C)", choices = c("Auto-tune via CV grid (default)" = "auto", "Manual" = "manual"), selected = "auto"),
            conditionalPanel(condition = sprintf("input['%s'] == 'auto'", ns("svm_cost_mode")),
                              textInput(ns("svm_cost_grid"), "Cost grid (comma-separated)", value = paste(DIAG_SVM_COST_GRID, collapse = ", "))),
            conditionalPanel(condition = sprintf("input['%s'] == 'manual'", ns("svm_cost_mode")),
                              numericInput(ns("svm_cost_manual"), "Cost value", value = 1, min = 0.001, step = 0.1))
          )
        ),
        h5("Advanced", style = "margin-top: 6px;"),
        p(class = "submodule-desc", "Gamma/degree only apply to Radial/Polynomial kernels."),
        fluidRow(
          column(3,
            radioButtons(ns("svm_gamma_mode"), "Gamma (kernel coefficient)", choices = c("Auto: 1 / number of genes (default)" = "auto", "Manual" = "manual"), selected = "auto"),
            conditionalPanel(condition = sprintf("input['%s'] == 'manual'", ns("svm_gamma_mode")),
                              numericInput(ns("svm_gamma_manual"), "Gamma value", value = 1, min = 0.0001, step = 0.01))
          ),
          column(3, conditionalPanel(condition = sprintf("input['%s'] == 'polynomial'", ns("svm_kernel")),
                                      numericInput(ns("svm_degree"), "Polynomial degree", value = DIAG_DEFAULT_PARAMS$svm_degree, min = 1, max = 10, step = 1))),
          column(3, numericInput(ns("svm_tolerance"), "Convergence tolerance", value = DIAG_DEFAULT_PARAMS$svm_tolerance, min = 0.00001, max = 0.1, step = 0.0001))
        )
      )
    }

    output$model_params_ui <- renderUI({
      if (!diag_has_run()) return(NULL)
      switch(active_model_pill(),
        "Elastic Net" = enet_params_box(),
        "Random Forest" = rf_params_box(),
        "SVM" = svm_params_box(),
        lr_params_box()
      )
    })

    output$references_box_ui <- renderUI({
      if (!diag_has_run()) return(NULL)
      box(
        width = 12, title = "References", status = "primary", solidHeader = FALSE,
        tags$ul(
          class = "dge-ref-list",
          tags$li(strong("Logistic regression: "), "Hosmer DW, Lemeshow S, Sturdivant RX (2013). ", tags$em("Applied Logistic Regression"), ", 3rd ed. Wiley."),
          tags$li(strong("Elastic net (glmnet): "), "Friedman J, Hastie T, Tibshirani R (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. ", tags$em("Journal of Statistical Software"), ", 33(1); Zou H, Hastie T (2005). Regularization and Variable Selection via the Elastic Net. ", tags$em("J R Stat Soc B"), ", 67(2)."),
          tags$li(strong("Random forests: "), "Breiman L (2001). Random Forests. ", tags$em("Machine Learning"), ", 45, 5-32."),
          tags$li(strong("SVM: "), "Cortes C, Vapnik V (1995). Support-Vector Networks. ", tags$em("Machine Learning"), ", 20, 273-297."),
          tags$li(strong("ROC / AUC confidence intervals: "), "DeLong ER, et al. (1988). Biometrics, 44, 837-845 (n ≥ 20); Carpenter J, Bithell J (2000). Bootstrap CIs. ", tags$em("Stat Med"), ", 19, 1141-1164 (n < 20)."),
          tags$li(strong("caret (hyperparameter tuning): "), "Kuhn M (2008). Building Predictive Models in R Using the caret Package. ", tags$em("Journal of Statistical Software"), ", 28(5).")
        ),
        p(class = "submodule-desc", strong("Ask ArthOChat"), " for a plain-language walkthrough or a live citation.")
      )
    })

    ## Saves each sex's fit into results$diagnostic independently as it finishes.
    save_result <- function(sex_label, r) {
      results$diagnostic <- utils::modifyList(
        results$diagnostic %||% list(),
        setNames(list(list(
          n_input = r$n_input, n_samples = r$n_samples,
          lr_auc = round(r$lr$full_auc, 3), lr_cv_auc = round(mean(r$lr$cv_auc, na.rm = TRUE), 3),
          enet_auc = round(r$enet$full_auc, 3), enet_cv_auc = round(mean(r$enet$cv_auc, na.rm = TRUE), 3),
          rf_auc = round(r$rf$full_auc, 3), rf_cv_auc = round(mean(r$rf$cv_auc, na.rm = TRUE), 3),
          svm_auc = round(r$svm$full_auc, 3), svm_cv_auc = round(mean(r$svm$cv_auc, na.rm = TRUE), 3),
          genes = r$genes
        )), sex_label)
      )
      showNotification(
        sprintf("%s diagnostic models saved: logistic regression CV-AUC %.3f, elastic net %.3f, random forest %.3f, SVM %.3f.",
                tools::toTitleCase(sex_label), mean(r$lr$cv_auc, na.rm = TRUE), mean(r$enet$cv_auc, na.rm = TRUE), mean(r$rf$cv_auc, na.rm = TRUE), mean(r$svm$cv_auc, na.rm = TRUE)),
        type = "message", duration = 6
      )
    }
    observeEvent(diag_result_female(), save_result("female", diag_result_female()))
    observeEvent(diag_result_male(), save_result("male", diag_result_male()))
    observeEvent(diag_result_pooled(), save_result("pooled", diag_result_pooled()))

    output$saved_runs_ui <- renderUI({
      res_f <- diag_result_value("female")
      res_m <- diag_result_value("male")
      res_p <- diag_result_value("pooled")
      status_row <- function(sex, sex_label, r) {
        if (!is.null(r)) {
          best_label <- c("Logistic Regression", "Elastic Net", "Random Forest", "SVM")[which.max(c(mean(r$lr$cv_auc, na.rm = TRUE), mean(r$enet$cv_auc, na.rm = TRUE), mean(r$rf$cv_auc, na.rm = TRUE), mean(r$svm$cv_auc, na.rm = TRUE)))]
          return(tags$li(icon("check", style = "color: #1a9c5f;"), strong(sprintf(" %s completed: ", sex)),
                          sprintf("best CV-AUC = %s", best_label)))
        }
        ## Real validate() failure vs. genuinely never clicked - see diag_result_error_msg() above.
        err <- diag_result_error_msg(sex_label)
        if (!is.null(err)) {
          tags$li(icon("triangle-exclamation", style = "color: #c0392b;"), strong(sprintf(" %s failed: ", sex)), err)
        } else {
          tags$li(icon("circle-minus", style = "color: #8A929C;"), sprintf(" %s diagnostic models - not run yet", sex))
        }
      }
      tagList(
        p(class = "submodule-desc", style = "margin-bottom: 4px;", "Status:"),
        tags$ul(style = "padding-left: 18px; margin-bottom: 0; list-style: none;",
                status_row("Female", "female", res_f), status_row("Male", "male", res_m), status_row("Pooled (all)", "pooled", res_p))
      )
    })

    ## One-line results summary shown inside that sex's own tab.
    diag_result_line <- function(sex_label, r) {
      if (!is.null(r)) {
        return(p(strong("Result: "),
          sprintf("%d genes, %d samples (%d vs %d), %s vs %s → CV-AUC logistic regression %.3f / elastic net %.3f / random forest %.3f / SVM %.3f.",
                  r$n_input, r$n_samples, r$n_comp, r$n_ref, r$comp_group, r$ref_group,
                  mean(r$lr$cv_auc, na.rm = TRUE), mean(r$enet$cv_auc, na.rm = TRUE), mean(r$rf$cv_auc, na.rm = TRUE), mean(r$svm$cv_auc, na.rm = TRUE))))
      }
      err <- diag_result_error_msg(sex_label)
      if (!is.null(err)) p(strong("Result: "), span(style = "color: #c0392b;", sprintf("failed - %s", err))) else NULL
    }
    output$female_result_line <- renderUI({ diag_result_line("female", diag_result_value("female")) })
    output$male_result_line <- renderUI({ diag_result_line("male", diag_result_value("male")) })
    output$pooled_result_line <- renderUI({ diag_result_line("pooled", diag_result_value("pooled")) })

    ## Per-technique, per-sex (or pooled) outputs.
    res_sex <- function(sex_label) reactive({ diag_result_value(sex_label) })

    build_train_perf_table <- function(rr) {
      rows <- list(
        data.frame(dataset = "Training (full fit)", metric = "AUC", value = sprintf("%.3f", rr$full_auc), stringsAsFactors = FALSE),
        data.frame(dataset = "Training (full fit)", metric = "Threshold (prob., Youden)", value = sprintf("%.3f", rr$best$threshold), stringsAsFactors = FALSE),
        data.frame(dataset = "Training (full fit)", metric = "Sensitivity", value = sprintf("%.3f", rr$best$sensitivity), stringsAsFactors = FALSE),
        data.frame(dataset = "Training (full fit)", metric = "Specificity", value = sprintf("%.3f", rr$best$specificity), stringsAsFactors = FALSE),
        data.frame(dataset = "Training (full fit)", metric = "Accuracy", value = sprintf("%.3f", rr$best$accuracy), stringsAsFactors = FALSE),
        data.frame(dataset = sprintf("Training (%d-fold CV)", length(rr$cv_auc)), metric = "Mean AUC", value = sprintf("%.3f", mean(rr$cv_auc, na.rm = TRUE)), stringsAsFactors = FALSE),
        data.frame(dataset = sprintf("Training (%d-fold CV)", length(rr$cv_auc)), metric = "SD AUC", value = sprintf("%.3f", stats::sd(rr$cv_auc, na.rm = TRUE)), stringsAsFactors = FALSE)
      )
      do.call(rbind, rows)
    }

    ## Test-split eval block as table rows.
    build_eval_rows <- function(ev, label) {
      if (!isTRUE(ev$available)) {
        return(list(data.frame(dataset = label, metric = "Status", value = ev$reason, stringsAsFactors = FALSE)))
      }
      auc_str <- paste0(sprintf("%.3f (%.3f-%.3f)", ev$auc, ev$ci_lo, ev$ci_hi), diag_separation_note(ev))
      list(
        data.frame(dataset = label, metric = "AUC (95% CI)", value = auc_str, stringsAsFactors = FALSE),
        data.frame(dataset = label, metric = "n samples (n positive)", value = sprintf("%d (%d)", ev$n, ev$n_pos), stringsAsFactors = FALSE),
        data.frame(dataset = label, metric = "Sensitivity @ training cutoff", value = sprintf("%.3f", ev$perf$sensitivity), stringsAsFactors = FALSE),
        data.frame(dataset = label, metric = "Specificity @ training cutoff", value = sprintf("%.3f", ev$perf$specificity), stringsAsFactors = FALSE),
        data.frame(dataset = label, metric = "Accuracy @ training cutoff", value = sprintf("%.3f", ev$perf$accuracy), stringsAsFactors = FALSE)
      )
    }

    build_test_perf_table <- function(rr) {
      do.call(rbind, build_eval_rows(rr$test, "Test split (held out from Train)"))
    }

    build_model_bundle <- function(r, key, sex_label) {
      rr <- r[[key]]
      list(
        model_type = rr$model_type, model_label = rr$label, model = rr$model,
        sex = sex_label, genes = r$genes,
        ref_group = r$ref_group, comp_group = r$comp_group,
        hyperparams = switch(rr$model_type,
          lr = list(),
          enet = list(alpha = rr$alpha, lambda_choice = rr$lambda_choice, lambda_used = rr$lambda_used),
          rf = list(ntree = rr$ntree, mtry = rr$mtry),
          svm = list(kernel = rr$kernel, cost = rr$cost)
        ),
        training_auc = rr$full_auc, training_cv_mean_auc = mean(rr$cv_auc, na.rm = TRUE),
        test_split_auc = if (isTRUE(rr$test$available)) rr$test$auc else NA_real_,
        scoring_note = paste(
          "To score new samples with `model`: for each gene in `genes`, z-score it across the NEW dataset's own samples",
          "(subtract that gene's mean, divide by its SD within the new data - independent of this training set's own mean/SD,",
          "matching how this project's own script transfers across cohorts); a gene in `genes` absent from the new data becomes",
          "a column of 0. Align columns to `genes` order (via make.names if needed), then: logistic regression -",
          "predict(model, newdata = as.data.frame(X), type = \"response\"); elastic net -",
          "predict(model, newx = X, s = hyperparams$lambda_choice, type = \"response\"); random forest -",
          "predict(model, X, type = \"prob\")[, comp_group]; SVM -",
          "attr(predict(model, X, probability = TRUE), \"probabilities\")[, comp_group]."
        ),
        trained_at = as.character(Sys.time())
      )
    }

    register_sex_model_outputs <- function(sex_label, res) {
      sex_color <- switch(sex_label, female = "#1a7a3c", male = "#7a4a26", pooled = "#2563EB")
      test_color <- ARTHOMIX_COLORS$orange
      ## Shows the real validate()/need() failure (e.g. "no live gene panel - run
      ## Feature Selection first") instead of the generic "not run yet" note when
      ## this sex's Run button WAS clicked but diag_build_sex() failed - see
      ## diag_result_error_msg() above for how the two are told apart.
      not_yet_note <- function() {
        err <- diag_result_error_msg(sex_label)
        if (!is.null(err)) {
          div(class = "empty-note", icon("triangle-exclamation"), sprintf("%s diagnostic models failed: %s", tools::toTitleCase(sex_label), err))
        } else {
          div(class = "empty-note", icon("circle-info"), "Not run yet. Click Run above.")
        }
      }

      lapply(DIAG_TECHNIQUES, function(tech) {
        key <- tech$key; label <- tech$label
        prefix <- paste0(sex_label, "_", key)

        ## ---- Model Training tab ----
        output[[paste0(prefix, "_train_stats")]] <- renderUI({
          r <- res()
          if (is.null(r)) return(not_yet_note())
          rr <- r[[key]]
          cv_mean <- mean(rr$cv_auc, na.rm = TRUE)
          tagList(
            fluidRow(
              valueBox(sprintf("%.3f", rr$full_auc), "Train AUC", icon = icon("chart-line"), color = "light-blue", width = 4),
              valueBox(sprintf("%.3f", cv_mean),
                       sprintf("%d-fold CV AUC (± %.3f)", length(rr$cv_auc), stats::sd(rr$cv_auc, na.rm = TRUE)),
                       icon = icon("layer-group"), color = "purple", width = 4),
              valueBox(diag_hyperparam_value(rr), "Hyperparameter", icon = icon("sliders"), color = "yellow", width = 4)
            ),
            ## Overfitting flag: near-perfect train AUC with a much weaker CV AUC
            ## suggests the model memorised the training samples.
            if (isTRUE(rr$full_auc >= 0.95) && isTRUE((rr$full_auc - cv_mean) >= 0.25)) {
              div(class = "empty-note", style = "border-left: 3px solid #c0392b; padding-left: 8px;",
                  icon("triangle-exclamation"),
                  sprintf(
                    "Likely overfitting: Train AUC (%.3f) far exceeds %d-fold CV AUC (%.3f) with %d genes on %d samples%s. Trust the CV-AUC.",
                    rr$full_auc, length(rr$cv_auc), cv_mean, r$n_input, r$n_samples,
                    if (identical(rr$model_type, "lr")) " and no regularisation" else ""
                  ))
            }
          )
        })
        output[[paste0(prefix, "_train_roc_plot")]] <- renderPlot({
          r <- res(); if (is.null(r)) return(NULL)
          rr <- r[[key]]
          diag_roc_plot_traintest(rr$roc_full, rr$test,
                                   cv_mean = mean(rr$cv_auc, na.rm = TRUE), cv_sd = stats::sd(rr$cv_auc, na.rm = TRUE), cv_n = length(rr$cv_auc),
                                   title = sprintf("ROC - %s (%s)", label, tools::toTitleCase(sex_label)))
        })
        output[[paste0(prefix, "_train_cv_plot")]] <- renderPlot({
          r <- res(); if (is.null(r)) return(NULL)
          rr <- r[[key]]
          df <- data.frame(fold = factor(seq_along(rr$cv_auc)), auc = rr$cv_auc)
          ggplot(df, aes(x = fold, y = auc)) +
            geom_col(fill = sex_color) +
            geom_hline(yintercept = mean(rr$cv_auc, na.rm = TRUE), linetype = "dashed", color = ARTHOMIX_COLORS$red) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(x = "Fold", y = "CV AUC (training)") + theme_arthomix(base_size = 12)
        })
        ## Plots the full hyperparameter search (alpha/mtry/cost vs CV score), chosen value marked.
        output[[paste0(prefix, "_train_tuning_plot")]] <- renderPlot({
          r <- res(); if (is.null(r)) return(NULL)
          rr <- r[[key]]
          ts <- rr$tuning_search
          if (is.null(ts) || nrow(ts) < 2) {
            msg <- if (identical(rr$model_type, "lr")) {
              "Plain logistic regression has no\nregularisation path - nothing to tune."
            } else {
              "Manual value used -\nno tuning grid was searched."
            }
            return(
              ggplot() +
                annotate("text", x = 0, y = 0, label = msg, size = 4.2, color = "#64748B") +
                theme_void()
            )
          }
          switch(rr$model_type,
            enet = {
              metric_lab <- switch(rr$type_measure %||% "deviance",
                auc = "CV AUC (higher is better)",
                class = "CV misclassification error (lower is better)",
                "CV deviance (lower is better)"
              )
              ggplot(ts, aes(x = factor(alpha), y = cv_metric, fill = chosen)) +
                geom_col() +
                scale_fill_manual(values = c(`TRUE` = sex_color, `FALSE` = "#CBD5E1"), guide = "none") +
                labs(x = "Alpha", y = metric_lab, title = sprintf("Selected alpha = %.2f", rr$alpha)) +
                theme_arthomix(base_size = 12)
            },
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
        })
        output[[paste0(prefix, "_train_table")]] <- DT::renderDataTable({
          r <- res(); if (is.null(r)) return(NULL)
          DT::datatable(build_train_perf_table(r[[key]]), rownames = FALSE, width = "100%",
                        options = list(pageLength = 15, dom = "t", scrollX = TRUE, autoWidth = FALSE), class = "stripe hover compact")
        })
        output[[paste0(prefix, "_train_download")]] <- downloadHandler(
          filename = function() sprintf("%s_%s_training_performance.csv", sex_label, key),
          content = function(file) write.csv(build_train_perf_table(res()[[key]]), file, row.names = FALSE)
        )
        output[[paste0(prefix, "_model_download")]] <- downloadHandler(
          filename = function() sprintf("%s_%s_model.rds", sex_label, key),
          content = function(file) saveRDS(build_model_bundle(res(), key, sex_label), file)
        )

        ## ---- Model Testing (Internal) tab ----
        output[[paste0(prefix, "_test_summary")]] <- renderUI({
          r <- res()
          if (is.null(r)) return(not_yet_note())
          rr <- r[[key]]
          test_tile <- if (isTRUE(rr$test$available)) {
            valueBox(sprintf("%.3f", rr$test$auc), sprintf("Test-split AUC (n=%d)%s", rr$test$n, diag_separation_note(rr$test)),
                     icon = icon("flask"), color = "light-blue", width = 6)
          } else {
            valueBox("N/A", "Test-split AUC", icon = icon("flask"), color = "red", width = 6)
          }
          n_tile <- valueBox(if (isTRUE(rr$test$available)) rr$test$n else "-", "Test samples held out", icon = icon("users"), color = "purple", width = 6)
          tagList(
            fluidRow(test_tile, n_tile),
            if (!isTRUE(rr$test$available)) p(class = "submodule-desc", style = "font-size: 12.5px;", icon("circle-info"), rr$test$reason %||% "Test split unavailable.") else NULL
          )
        })
        output[[paste0(prefix, "_test_roc_plot")]] <- renderPlot({
          r <- res(); if (is.null(r)) return(NULL)
          rr <- r[[key]]
          req(isTRUE(rr$test$available))
          diag_roc_plot_pub(rr$test$roc, color = test_color,
                             title = sprintf("Test-split ROC - %s (%s)", label, tools::toTitleCase(sex_label)),
                             ci = c(rr$test$ci_lo, rr$test$auc, rr$test$ci_hi))
        })
        output[[paste0(prefix, "_test_table")]] <- DT::renderDataTable({
          r <- res(); if (is.null(r)) return(NULL)
          DT::datatable(build_test_perf_table(r[[key]]), rownames = FALSE, width = "100%",
                        options = list(pageLength = 15, dom = "t", scrollX = TRUE, autoWidth = FALSE), class = "stripe hover compact")
        })
        output[[paste0(prefix, "_test_download")]] <- downloadHandler(
          filename = function() sprintf("%s_%s_testing_performance.csv", sex_label, key),
          content = function(file) write.csv(build_test_perf_table(res()[[key]]), file, row.names = FALSE)
        )
      })

      output[[paste0(sex_label, "_train_compare_table")]] <- DT::renderDataTable({
        r <- res(); if (is.null(r)) return(NULL)
        df <- do.call(rbind, lapply(DIAG_TECHNIQUES, function(tech) {
          rr <- r[[tech$key]]
          data.frame(model = tech$label, training_auc = round(rr$full_auc, 3),
                     cv_mean_auc = round(mean(rr$cv_auc, na.rm = TRUE), 3), cv_sd_auc = round(stats::sd(rr$cv_auc, na.rm = TRUE), 3),
                     stringsAsFactors = FALSE)
        }))
        DT::datatable(df, rownames = FALSE, width = "100%", options = list(pageLength = 5, dom = "t", scrollX = TRUE, autoWidth = FALSE), class = "stripe hover compact")
      })
      output[[paste0(sex_label, "_test_compare_table")]] <- DT::renderDataTable({
        r <- res(); if (is.null(r)) return(NULL)
        df <- do.call(rbind, lapply(DIAG_TECHNIQUES, function(tech) {
          rr <- r[[tech$key]]
          data.frame(
            model = tech$label,
            test_split_auc = if (isTRUE(rr$test$available)) round(rr$test$auc, 3) else NA_real_,
            note = trimws(diag_separation_note(rr$test)),
            stringsAsFactors = FALSE
          )
        }))
        DT::datatable(df, rownames = FALSE, width = "100%", options = list(pageLength = 5, dom = "t", scrollX = TRUE, autoWidth = FALSE), class = "stripe hover compact")
      })

      ## Per-gene ROC/AUC: Train and Test curves overlaid per gene panel.
      register_generoc_outputs <- function(mode) {
        plot_id <- paste0(sex_label, "_", mode, "_generoc_plot")
        table_id <- paste0(sex_label, "_", mode, "_generoc_table")
        download_id <- paste0(sex_label, "_", mode, "_generoc_download")
        hub_download_id <- paste0(sex_label, "_", mode, "_hub_download")
        auc_thr_id <- paste0(sex_label, "_", mode, "_hub_auc_thr")
        p_thr_id <- paste0(sex_label, "_", mode, "_hub_p_thr")

        ## "hub" is decided on training stats only; `gr$p` falls back to NA for older cached runs.
        gene_auc_df <- function(gr_tr, gr_te) {
          p_tr <- if (!is.null(gr_tr$p)) unname(gr_tr$p) else rep(NA_real_, length(gr_tr$genes))
          df <- data.frame(gene = gr_tr$genes, train_auc = round(unname(gr_tr$auc), 3), train_p = signif(p_tr, 3), stringsAsFactors = FALSE)
          idx <- match(df$gene, gr_te$genes)
          p_te <- if (!is.null(gr_te$p)) unname(gr_te$p) else rep(NA_real_, length(gr_te$genes))
          df$test_auc <- round(unname(gr_te$auc)[idx], 3)
          df$test_p <- signif(p_te[idx], 3)
          auc_thr <- input[[auc_thr_id]] %||% 0.85
          p_thr <- input[[p_thr_id]] %||% 0.05
          df$hub <- !is.na(df$train_auc) & !is.na(df$train_p) & df$train_auc >= auc_thr & df$train_p < p_thr
          df[order(-df$train_auc), ]
        }

        ## One small ROC panel per gene (Train+Test overlaid), capped to top genes by training AUC.
        GENEROC_MAX_FACETS <- 24
        output[[plot_id]] <- renderPlot({
          r <- res(); if (is.null(r)) return(NULL)
          gr_tr <- r$gene_roc_train; gr_te <- r$gene_roc_test
          auc_df <- gene_auc_df(gr_tr, gr_te)
          req(nrow(auc_df) > 0)
          top_genes <- head(auc_df$gene, GENEROC_MAX_FACETS)
          gene_curve <- function(gr, g, ds) {
            rc <- gr$rocs[[g]]
            if (is.null(rc)) return(NULL)
            data.frame(gene = g, fpr = 1 - rc$specificities, tpr = rc$sensitivities, Dataset = ds, stringsAsFactors = FALSE)
          }
          df <- do.call(rbind, c(lapply(top_genes, function(g) gene_curve(gr_tr, g, "Train")),
                                  lapply(top_genes, function(g) gene_curve(gr_te, g, "Test"))))
          req(nrow(df) > 0)
          df$gene <- factor(df$gene, levels = top_genes)
          df$Dataset <- factor(df$Dataset, levels = c("Test", "Train"))
          auc_lab <- auc_df[match(top_genes, auc_df$gene), c("gene", "train_auc", "test_auc")]
          auc_lab$gene <- factor(auc_lab$gene, levels = top_genes)
          auc_lab$label <- ifelse(!is.na(auc_lab$test_auc),
                                   sprintf("Train AUC=%.3f\nTest AUC=%.3f", auc_lab$train_auc, auc_lab$test_auc),
                                   sprintf("Train AUC=%.3f", auc_lab$train_auc))
          n_col <- min(4, length(top_genes))
          plot_title <- if (nrow(auc_df) > GENEROC_MAX_FACETS) {
            sprintf("Top %d of %d genes by training AUC", GENEROC_MAX_FACETS, nrow(auc_df))
          } else NULL
          ggplot(df, aes(x = fpr, y = tpr, color = Dataset)) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#9CA3AF", linewidth = 0.5) +
            geom_line(linewidth = 1) +
            scale_color_manual(name = "Dataset", values = c(Train = ARTHOMIX_COLORS$yellow, Test = ARTHOMIX_COLORS$blue), breaks = c("Test", "Train")) +
            geom_text(data = auc_lab, aes(x = 0.98, y = 0.04, label = label), inherit.aes = FALSE,
                      hjust = 1, vjust = 0, size = 2.7, color = "#1F2937", lineheight = 0.9) +
            facet_wrap(~gene, ncol = n_col) +
            scale_x_continuous(name = "1 − Specificity (FPR)", limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0.02, 0.02)) +
            scale_y_continuous(name = "Sensitivity (TPR)", limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0.02, 0.02)) +
            coord_equal() +
            labs(title = plot_title) +
            theme_bw(base_size = 11) +
            theme(
              panel.grid.minor = element_blank(),
              panel.grid.major = element_line(color = "#EEF0F3", linewidth = 0.3),
              strip.background = element_rect(fill = "#F1F5F9", color = NA),
              strip.text = element_text(face = "bold", size = 10),
              plot.title = element_text(size = 10.5, color = "#6B7280"),
              axis.title = element_text(face = "bold", size = 11),
              axis.text = element_text(size = 9, color = "black"),
              panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA),
              legend.position = "bottom"
            )
        }, height = function() {
          r <- tryCatch(res(), error = function(e) NULL)
          if (is.null(r)) return(320)
          n <- min(nrow(gene_auc_df(r$gene_roc_train, r$gene_roc_test)), GENEROC_MAX_FACETS)
          if (n <= 0) return(320)
          n_col <- min(4, n)
          n_row <- ceiling(n / n_col)
          max(320, n_row * 230)
        })
        output[[table_id]] <- DT::renderDataTable({
          r <- res(); if (is.null(r)) return(NULL)
          DT::datatable(gene_auc_df(r$gene_roc_train, r$gene_roc_test), rownames = FALSE, width = "100%",
                        options = list(pageLength = 8, scrollX = TRUE, autoWidth = FALSE), class = "stripe hover compact") |>
            DT::formatStyle("hub", target = "row", backgroundColor = DT::styleEqual(c(TRUE, FALSE), c("#e6f4ea", "")))
        })
        output[[download_id]] <- downloadHandler(
          filename = function() sprintf("%s_gene_auc.csv", sex_label),
          content = function(file) write.csv(gene_auc_df(res()$gene_roc_train, res()$gene_roc_test), file, row.names = FALSE)
        )
        ## Just the rows currently passing the AUC/P thresholds set above.
        output[[hub_download_id]] <- downloadHandler(
          filename = function() sprintf("%s_hub_genes.csv", sex_label),
          content = function(file) {
            df <- gene_auc_df(res()$gene_roc_train, res()$gene_roc_test)
            write.csv(df[df$hub, c("gene", "train_auc", "train_p", "test_auc", "test_p")], file, row.names = FALSE)
          }
        )
      }
      register_generoc_outputs("train")
    }

    register_sex_model_outputs("female", res_sex("female"))
    register_sex_model_outputs("male", res_sex("male"))
    register_sex_model_outputs("pooled", res_sex("pooled"))

  })
}
