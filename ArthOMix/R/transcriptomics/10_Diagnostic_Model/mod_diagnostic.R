## R/transcriptomics/10_Diagnostic_Model/mod_diagnostic.R
## Diagnostic Model submodule: fits logistic regression, elastic net, random
## forest and SVM classifiers on a user-chosen gene panel, sex-stratified.

mod_diagnostic_config <- list(
  id = "diagnostic", group = "Biomarker modeling",
  title = "Diagnostic Model",
  description = "Diagnostic model using logistic regression, elastic net, random forest and SVM diagnostic models, by sex.",
  icon = "stethoscope"
)

diag_zrows <- function(M) {
  t(apply(M, 1, function(v) {
    s <- stats::sd(v, na.rm = TRUE)
    if (is.na(s) || s == 0) rep(0, length(v)) else (v - mean(v, na.rm = TRUE)) / s
  }))
}

DIAG_SVM_COST_GRID <- c(0.01, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16)

DIAG_DEFAULT_PARAMS <- list(
  test_frac = 0.3,
  class_weight_mode = "equal", class_weight_ratio = 1,
  lr_cv_folds = 5,
  enet_cv_folds = 5, enet_alpha_grid = c(0.1, 0.3, 0.5, 0.7, 0.9, 1.0), enet_lambda_choice = "lambda.min",
  enet_nlambda = 100, enet_type_measure = "deviance",
  rf_cv_folds = 5, rf_ntree = 1000, rf_mtry_mode = "auto", rf_mtry_manual = NULL,
  rf_nodesize = 1, rf_maxnodes = NULL,
  svm_cv_folds = 5, svm_kernel = "linear", svm_cost_mode = "auto", svm_cost_manual = 1, svm_cost_grid = DIAG_SVM_COST_GRID,
  svm_gamma_mode = "auto", svm_gamma_manual = 1, svm_degree = 3, svm_tolerance = 0.001
)

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

diag_split_train_test <- function(y_full, test_frac, seed = 1234, sample_ids = NULL, holdout_ids = NULL) {
  if (!is.null(sample_ids) && length(holdout_ids) > 0) {
    test_idx <- which(sample_ids %in% holdout_ids)
    train_idx <- setdiff(seq_along(y_full), test_idx)
    if (length(test_idx) >= 4 && length(train_idx) >= 10 &&
        length(unique(y_full[train_idx])) == 2 && length(unique(y_full[test_idx])) == 2) {
      return(list(train = train_idx, test = test_idx, leakage_safe = TRUE))
    }
  }
  set.seed(seed)
  train_idx <- as.integer(caret::createDataPartition(y_full, p = 1 - test_frac, list = FALSE))
  list(train = train_idx, test = setdiff(seq_along(y_full), train_idx), leakage_safe = FALSE)
}

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

diag_separation_note <- function(ev) {
  if (isTRUE(ev$available) && !is.na(ev$ci_lo) && ev$ci_lo >= 0.999) {
    sprintf(" [separation, n=%d]", ev$n)
  } else ""
}

diag_hyperparam_value <- function(rr) {
  switch(rr$model_type,
    lr = "none (unpenalized)",
    enet = sprintf("α = %.2f", rr$alpha),
    rf = sprintf("mtry = %d", rr$mtry),
    svm = sprintf("cost = %s", format(rr$cost, trim = TRUE))
  )
}

diag_gene_roc <- function(expr_sub, y) {
  genes <- rownames(expr_sub)
  rocs <- vector("list", length(genes)); names(rocs) <- genes
  aucs <- setNames(numeric(length(genes)), genes)
  pvals <- setNames(numeric(length(genes)), genes)
  for (g in genes) {
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
  TP <- sum(pred_pos & obs_pos); TN <- sum(!pred_pos & !obs_pos)
  FP <- sum(pred_pos & !obs_pos); FN <- sum(!pred_pos & obs_pos)
  sens <- if (sum(obs_pos) > 0) TP / sum(obs_pos) else NA_real_
  spec <- if (sum(!obs_pos) > 0) TN / sum(!obs_pos) else NA_real_
  prec <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  npv  <- if ((TN + FN) > 0) TN / (TN + FN) else NA_real_
  f1   <- if (!is.na(prec) && !is.na(sens) && (prec + sens) > 0) 2 * prec * sens / (prec + sens) else NA_real_
  bal_acc <- if (!is.na(sens) && !is.na(spec)) (sens + spec) / 2 else NA_real_
  mcc_den <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- if (is.na(mcc_den) || mcc_den == 0) NA_real_ else (TP * TN - FP * FN) / mcc_den
  y_bin <- as.integer(obs_pos)
  brier <- mean((prob - y_bin)^2)
  pr_auc <- if (sum(y_bin == 1) > 0 && sum(y_bin == 0) > 0) {
    tryCatch(PRROC::pr.curve(scores.class0 = prob[y_bin == 1], scores.class1 = prob[y_bin == 0], curve = FALSE)$auc.integral,
              error = function(e) NA_real_)
  } else NA_real_
  list(
    sensitivity = sens, specificity = spec, accuracy = mean(pred_pos == obs_pos),
    precision = prec, npv = npv, f1 = f1, balanced_accuracy = bal_acc, mcc = mcc,
    brier = brier, pr_auc = pr_auc
  )
}

diag_calibration <- function(y, prob, positive_level, bins = 10) {
  y_bin <- as.integer(y == positive_level)
  brks <- unique(seq(0, 1, length.out = bins + 1))
  bin_id <- cut(prob, brks, include.lowest = TRUE)
  df <- data.frame(prob = prob, y = y_bin, bin = bin_id)
  agg <- stats::aggregate(cbind(mean_pred = prob, mean_obs = y) ~ bin, data = df, FUN = mean)
  n_bin <- as.data.frame(table(bin = df$bin)); names(n_bin) <- c("bin", "n")
  agg <- merge(agg, n_bin, by = "bin", all.x = TRUE)
  df$logit_prob <- stats::qlogis(pmin(pmax(prob, 1e-6), 1 - 1e-6))
  fit <- tryCatch(stats::glm(y ~ logit_prob, family = stats::binomial(), data = df), error = function(e) NULL)
  list(table = agg, brier = mean((prob - y_bin)^2),
       slope = if (!is.null(fit)) unname(stats::coef(fit)[2]) else NA_real_,
       intercept = if (!is.null(fit)) unname(stats::coef(fit)[1]) else NA_real_)
}

diag_plot_calibration <- function(cal, title) {
  validate(need(!is.null(cal) && nrow(cal$table) > 0, "Not enough held-out predictions to draw a calibration curve."))
  ggplot(cal$table, aes(x = mean_pred, y = mean_obs)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_line(color = ARTHOMIX_COLORS$blue) + geom_point(aes(size = n), color = ARTHOMIX_COLORS$blue) +
    scale_x_continuous(name = "Mean predicted probability", limits = c(0, 1)) +
    scale_y_continuous(name = "Observed proportion", limits = c(0, 1)) +
    labs(title = title, subtitle = sprintf("Brier score = %.3f", cal$brier), size = "Bin n") +
    theme_arthomix(base_size = 11)
}

diag_fit_sex <- function(expr_full, y_full, params = list(), holdout_ids = character(0)) {
  levels(y_full) <- make.names(levels(y_full), unique = TRUE)
  params <- utils::modifyList(DIAG_DEFAULT_PARAMS, params)
  GLOBAL_SEED <- ARTHOMIX_TX_ML_SEED
  genes <- rownames(expr_full)
  safe <- make.names(genes, unique = TRUE)

  split <- diag_split_train_test(y_full, params$test_frac, seed = GLOBAL_SEED,
                                  sample_ids = colnames(expr_full), holdout_ids = holdout_ids)
  validate(need(length(split$train) >= 10 && length(split$test) >= 4,
                "Not enough samples for a train/test split at this ratio - lower the test-set size or provide more samples for this sex."))
  y <- y_full[split$train]
  validate(need(length(unique(y)) == 2 && all(table(y) >= 3),
                "The training split ended up with fewer than 3 samples in one group - lower the test-set size or provide more samples for this sex."))
  ytest <- y_full[split$test]
  validate(need(length(unique(ytest)) == 2,
                "The test split ended up with only one group present - lower the test-set size or provide more samples for this sex."))

  cw_levels <- diag_class_weight_levels(y, params$class_weight_mode, params$class_weight_ratio)
  obs_w <- diag_obs_weights(y, params$class_weight_mode, params$class_weight_ratio)

  expr_train_sub <- expr_full[, split$train, drop = FALSE]
  expr_test_sub <- expr_full[, split$test, drop = FALSE]

  Ztr <- diag_zrows(expr_train_sub); rownames(Ztr) <- safe
  Xtr_full <- t(Ztr)                                    # z-scored, sample x gene(safe) - used for tuning + full fit
  Xraw <- t(expr_train_sub); colnames(Xraw) <- safe      # raw, sample x gene(safe) - re-scaled per fold inside diag_cv_auc

  mu_tr <- colMeans(Xraw); sg_tr <- apply(Xraw, 2, stats::sd); sg_tr[is.na(sg_tr) | sg_tr == 0] <- 1
  Xtest_full <- scale(t(expr_test_sub), center = mu_tr, scale = sg_tr)
  colnames(Xtest_full) <- safe

  youden <- function(roc_obj) {
    b <- pROC::coords(roc_obj, "best", best.method = "youden",
                       ret = c("threshold", "sensitivity", "specificity", "accuracy"), transpose = FALSE)
    if (nrow(b) > 1) b <- b[1, , drop = FALSE]
    b
  }

  score_eval <- function(pred_eval, y_eval, avail, reason) {
    if (!avail) return(list(available = FALSE, reason = reason))
    roc_e <- tryCatch(pROC::roc(y_eval, pred_eval, quiet = TRUE, levels = levels(y), direction = "<"), error = function(e) NULL)
    if (is.null(roc_e)) return(list(available = FALSE, reason = "ROC could not be computed for this split/contrast."))
    ci <- diag_auc_ci(roc_e)
    list(available = TRUE, roc = roc_e, auc = unname(ci["auc"]), ci_lo = unname(ci["lo"]), ci_hi = unname(ci["hi"]),
         n = length(y_eval), n_pos = sum(y_eval == levels(y)[2]), pred = pred_eval, y_eval = y_eval)
  }

  nf_a <- max(2, min(params$enet_cv_folds, min(table(y))))
  type_measure <- params$enet_type_measure %||% "deviance"
  bigger_is_better <- identical(type_measure, "auc")
  best <- NULL; bcv <- if (bigger_is_better) -Inf else Inf
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

  p <- ncol(Xtr_full)
  ntree <- max(100, round(params$rf_ntree))
  rf_nodesize <- max(1, round(params$rf_nodesize %||% 1))
  rf_maxnodes <- if (!is.null(params$rf_maxnodes) && is.finite(params$rf_maxnodes)) max(2, round(params$rf_maxnodes)) else NULL
  mtry_search <- NULL
  if (identical(params$rf_mtry_mode, "manual") && !is.null(params$rf_mtry_manual)) {
    rf_mtry <- min(p, max(1, round(params$rf_mtry_manual)))
  } else {
    nf_rf <- max(2, min(params$rf_cv_folds, min(table(y))))
    mtry_grid <- sort(unique(pmax(1, pmin(p, c(1, 2, floor(sqrt(p)), floor(p / 3), floor(p / 2), p)))))
    ctrl <- caret::trainControl(method = "cv", number = nf_rf, classProbs = TRUE, summaryFunction = caret::twoClassSummary)
    set.seed(GLOBAL_SEED)
    rf_tune <- tryCatch(caret::train(x = Xtr_full, y = y, method = "rf", metric = "ROC", trControl = ctrl,
                                      tuneGrid = expand.grid(mtry = mtry_grid), ntree = ntree,
                                      nodesize = rf_nodesize, maxnodes = rf_maxnodes, classwt = cw_levels), error = function(e) NULL)
    rf_mtry <- if (!is.null(rf_tune)) rf_tune$bestTune$mtry else max(1, floor(sqrt(p)))
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
       gene_roc_test = diag_gene_roc(expr_test_sub, ytest),
       leakage_safe = isTRUE(split$leakage_safe))
}

DIAG_TECHNIQUES <- list(
  list(key = "lr", label = "Logistic Regression"),
  list(key = "enet", label = "Elastic Net"),
  list(key = "rf", label = "Random Forest"),
  list(key = "svm", label = "SVM")
)

diag_apply_models_external <- function(r, expr_ext, y_ext) {
  genes <- r$genes
  safe <- make.names(genes, unique = TRUE)
  present <- genes %in% rownames(expr_ext)
  X <- matrix(0, nrow = ncol(expr_ext), ncol = length(genes), dimnames = list(colnames(expr_ext), safe))
  if (any(present)) {
    sub <- t(expr_ext[genes[present], , drop = FALSE])
    mu <- colMeans(sub, na.rm = TRUE)
    sg <- apply(sub, 2, stats::sd, na.rm = TRUE); sg[is.na(sg) | sg == 0] <- 1
    X[, safe[present]] <- scale(sub, center = mu, scale = sg)
  }
  X[is.na(X)] <- 0
  pos <- make.names(r$comp_group %||% levels(y_ext)[2])
  pick_pos <- function(pm) { cn <- colnames(pm); if (!is.null(cn) && pos %in% cn) pm[, pos] else pm[, ncol(pm)] }
  models <- lapply(DIAG_TECHNIQUES, function(tech) {
    rr <- r[[tech$key]]
    if (is.null(rr) || is.null(rr$model)) return(list(available = FALSE, label = tech$label, key = tech$key, reason = "Model not trained."))
    prob <- tryCatch(switch(rr$model_type,
      lr   = as.numeric(predict(rr$model, newdata = as.data.frame(X), type = "response")),
      enet = as.numeric(predict(rr$model, newx = X, s = rr$lambda_choice, type = "response")),
      rf   = as.numeric(pick_pos(predict(rr$model, X, type = "prob"))),
      svm  = as.numeric(pick_pos(attr(predict(rr$model, X, probability = TRUE), "probabilities")))
    ), error = function(e) NULL)
    if (is.null(prob) || length(prob) != length(y_ext) || all(is.na(prob))) {
      return(list(available = FALSE, label = tech$label, key = tech$key, reason = "Could not score the external samples with this model."))
    }
    roc_e <- tryCatch(pROC::roc(y_ext, prob, quiet = TRUE, levels = levels(y_ext), direction = "<"), error = function(e) NULL)
    if (is.null(roc_e)) return(list(available = FALSE, label = tech$label, key = tech$key, reason = "ROC could not be computed."))
    ci <- diag_auc_ci(roc_e)
    thr <- suppressWarnings(as.numeric(rr$best$threshold))
    if (!is.finite(thr)) thr <- 0.5
    perf <- diag_perf_at_cutoff(prob, y_ext, thr, levels(y_ext)[2])
    list(available = TRUE, label = tech$label, key = tech$key, roc = roc_e,
         auc = unname(ci["auc"]), ci_lo = unname(ci["lo"]), ci_hi = unname(ci["hi"]),
         threshold = thr, perf = perf, prob = prob, y = y_ext)
  })
  names(models) <- vapply(DIAG_TECHNIQUES, `[[`, character(1), "key")
  list(models = models, n_genes_present = sum(present), n_genes_panel = length(genes))
}

diag_validate_nested <- function(expr_candidates, y_full, outer_k = 5, uni_top_n = 100,
                                  lasso_alpha = 1, seed = ARTHOMIX_TX_ML_SEED) {
  y_full <- droplevels(factor(as.character(y_full)))
  validate(need(length(levels(y_full)) == 2, "Leakage-safe validation needs exactly two groups."))
  genes <- rownames(expr_candidates)
  safe <- make.names(genes, unique = TRUE)
  rownames(expr_candidates) <- safe

  nf <- max(2, min(outer_k, min(table(y_full))))
  set.seed(seed)
  folds <- caret::createFolds(y_full, k = nf, list = TRUE)

  pred_oof <- rep(NA_real_, length(y_full))
  per_fold <- list()

  for (fi in seq_along(folds)) {
    te <- folds[[fi]]; tr <- setdiff(seq_along(y_full), te)
    y_tr <- droplevels(y_full[tr]); y_te <- droplevels(y_full[te])
    if (length(unique(y_tr)) < 2 || length(unique(y_te)) < 2) next
    expr_tr <- expr_candidates[, tr, drop = FALSE]
    expr_te <- expr_candidates[, te, drop = FALSE]

    design <- stats::model.matrix(~y_tr)
    uni_fit <- tryCatch(limma::eBayes(limma::lmFit(expr_tr, design)), error = function(e) NULL)
    if (is.null(uni_fit)) next
    tt <- tryCatch(limma::topTable(uni_fit, coef = 2, number = Inf, sort.by = "P"), error = function(e) NULL)
    if (is.null(tt) || nrow(tt) < 2) next
    uni_genes <- rownames(tt)[seq_len(min(uni_top_n, nrow(tt)))]

    Xtr_raw <- t(expr_tr[uni_genes, , drop = FALSE])
    nf_lasso <- max(2, min(5, min(table(y_tr))))
    cv <- tryCatch(glmnet::cv.glmnet(Xtr_raw, y_tr, family = "binomial", alpha = lasso_alpha, nfolds = nf_lasso),
                   error = function(e) NULL)
    panel <- character(0)
    if (!is.null(cv)) {
      co <- coef(cv, s = "lambda.min")[-1, 1, drop = TRUE]
      panel <- names(co)[co != 0]
    }
    if (length(panel) < 1) panel <- utils::head(uni_genes, min(5, length(uni_genes)))

    mu <- rowMeans(expr_tr[panel, , drop = FALSE], na.rm = TRUE)
    sg <- apply(expr_tr[panel, , drop = FALSE], 1, stats::sd)
    sg[is.na(sg) | sg == 0] <- 1
    Ztr <- (expr_tr[panel, , drop = FALSE] - mu) / sg
    Zte <- (expr_te[panel, , drop = FALSE] - mu) / sg
    fit_df <- data.frame(y = y_tr, t(Ztr), check.names = FALSE)
    model <- tryCatch(suppressWarnings(stats::glm(y ~ ., data = fit_df, family = stats::binomial)), error = function(e) NULL)
    if (is.null(model)) next
    pred <- tryCatch(as.numeric(predict(model, newdata = data.frame(t(Zte), check.names = FALSE), type = "response")),
                      error = function(e) NULL)
    if (is.null(pred)) next
    pred_oof[te] <- pred

    roc_i <- tryCatch(pROC::roc(y_te, pred, quiet = TRUE, levels = levels(y_full), direction = "<"), error = function(e) NULL)
    auc_i <- if (!is.null(roc_i)) as.numeric(pROC::auc(roc_i)) else NA_real_
    per_fold[[length(per_fold) + 1]] <- data.frame(fold = fi, n_panel = length(panel), auc = round(auc_i, 3))
  }

  per_fold_df <- if (length(per_fold) > 0) do.call(rbind, per_fold) else data.frame(fold = integer(0), n_panel = integer(0), auc = numeric(0))
  have_oof <- !is.na(pred_oof)
  pooled <- list(available = FALSE, reason = "Not enough folds completed to score a pooled AUC.")
  if (sum(have_oof) >= 10 && length(unique(y_full[have_oof])) == 2) {
    roc_pooled <- tryCatch(pROC::roc(y_full[have_oof], pred_oof[have_oof], quiet = TRUE, levels = levels(y_full), direction = "<"),
                            error = function(e) NULL)
    if (!is.null(roc_pooled)) {
      ci <- diag_auc_ci(roc_pooled)
      pooled <- list(available = TRUE, auc = unname(ci["auc"]), ci_lo = unname(ci["lo"]), ci_hi = unname(ci["hi"]),
                      n = sum(have_oof), roc = roc_pooled)
    }
  }
  list(pooled = pooled, per_fold = per_fold_df, n_folds_completed = nrow(per_fold_df), outer_k = nf)
}

## Decides which AUC is this sex-run's headline (primary) metric. A genuine
## leakage-safe held-out split (leakage_safe == TRUE) keeps the naive
## Test-split AUC as primary - it's a valid estimate in that case. Otherwise
## (the default/bundled-panel path) the automatically-computed nested-CV AUC
## becomes primary whenever it's available; the naive Test-split AUC is
## always still attached/shown, just demoted in the UI. Pure and testable
## independent of the nested-CV computation itself.
diag_attach_headline <- function(fit, nested_cv = NULL) {
  if (isTRUE(fit$leakage_safe)) {
    fit$nested_cv <- NULL
    fit$headline_metric <- "test_split"
    return(fit)
  }
  fit$nested_cv <- nested_cv
  fit$headline_metric <- if (!is.null(nested_cv) && isTRUE(nested_cv$pooled$available)) "nested_cv" else "test_split"
  fit
}

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

mod_diagnostic_testing_panel <- function(ns, prefix, title, roc_height = 300) {
  box(
    width = NULL, title = title, status = "primary", solidHeader = FALSE,
    withSpinner(uiOutput(ns(paste0(prefix, "_test_summary"))), color = "#2563EB", type = 6),
    fluidRow(
      column(6, withSpinner(plotOutput(ns(paste0(prefix, "_test_roc_plot")), height = roc_height), color = "#2563EB", type = 6)),
      column(6, withSpinner(plotOutput(ns(paste0(prefix, "_test_calibration_plot")), height = roc_height), color = "#2563EB", type = 6))
    ),
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
      uiOutput(ns(paste0(sex_label, "_headline_ui"))),
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

mod_diagnostic_leakagesafe_sex_panel <- function(ns, sex_label) {
  run_id <- paste0("run_", sex_label, "_leakagesafe_btn")
  sex_title <- tools::toTitleCase(sex_label)
  tagList(
    p(class = "empty-note", icon("triangle-exclamation"),
      "The Test-split AUC in Model Testing (Internal) evaluates a gene panel that was already chosen using this same data - Candidate Gene Identification and Feature Selection's LASSO/RF/SVM-RFE consensus both run on the full sample pool, with no held-out split - so that AUC is optimistic. This mode reselects the panel with the Univariate + LASSO steps inside every outer fold instead, for a more honest estimate."),
    p(class = "empty-note", icon("circle-info"),
      "Leakage-safe mode reselects the panel from this sex's WGCNA-candidate gene set (Candidate Gene Identification's output) using Univariate ranking + LASSO inside every outer fold - not literally the same panel shown in Feature Selection - and it does not rerun WGCNA, Random Forest, or SVM-RFE per fold, since a full discovery-pipeline refit per fold is impractical in a live session. The candidate gene list itself, and its cap above 200 genes (kept by full-pool variance), are still computed once before this outer cross-validation runs - so \"leakage-safe\" here means the supervised feature-selection step is redone per fold, not that the estimate is completely free of leakage."),
    fluidRow(
      column(4, numericInput(ns(paste0(sex_label, "_leakagesafe_k")), "Outer folds (k)", value = 5, min = 3, max = 10, step = 1)),
      column(4, div(style = "margin-top: 25px;",
                    actionButton(ns(run_id), paste("Run", sex_title, "Leakage-safe Validation"), icon = icon("play"), class = "btn-primary btn-sm")))
    ),
    withSpinner(uiOutput(ns(paste0(sex_label, "_leakagesafe_result_ui"))), color = "#2563EB", type = 6)
  )
}

mod_diagnostic_ui <- function(id) {
  ns <- NS(id)
  tagList(
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
            uiOutput(ns("model_params_ui")),
            tabsetPanel(
              id = ns("train_sex_tabs"), type = "tabs",
              tabPanel("Female", br(), mod_diagnostic_training_sex_panel(ns, "female")),
              tabPanel("Male", br(), mod_diagnostic_training_sex_panel(ns, "male")),
              tabPanel("Pooled (all)", br(), mod_diagnostic_training_sex_panel(ns, "pooled"))
            )
          ),
          tabPanel(
            "Model Testing (Internal)", br(),
            p(class = "submodule-desc", "Each Train-fit model scored once on its held-out Test split. This scores the panel as-is - see \"Leakage-safe Validation\" for an estimate that also accounts for how that panel was chosen."),
            tabsetPanel(
              id = ns("test_sex_tabs"), type = "tabs",
              tabPanel("Female", br(), mod_diagnostic_testing_sex_panel(ns, "female")),
              tabPanel("Male", br(), mod_diagnostic_testing_sex_panel(ns, "male")),
              tabPanel("Pooled (all)", br(), mod_diagnostic_testing_sex_panel(ns, "pooled"))
            )
          ),
          tabPanel(
            "Leakage-safe Validation", br(),
            p(class = "submodule-desc", "Nested cross-validation that reselects the gene panel inside each outer fold instead of scoring the panel already chosen using this same data - see the disclosure under each sex's Run button. Needs a live Candidate Gene Identification run (\"Follow this project's pipeline\" panel source) on this dataset."),
            tabsetPanel(
              id = ns("leakagesafe_sex_tabs"), type = "tabs",
              tabPanel("Female", br(), mod_diagnostic_leakagesafe_sex_panel(ns, "female")),
              tabPanel("Male", br(), mod_diagnostic_leakagesafe_sex_panel(ns, "male")),
              tabPanel("Pooled (all)", br(), mod_diagnostic_leakagesafe_sex_panel(ns, "pooled"))
            )
          ),
          tabPanel(
            "External Validation", br(),
            mod_diagnostic_external_panel(ns)
          )
        )
      )
    )
    )
  )
}

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

mod_diagnostic_server <- function(id, dataset, results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

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

    project_panel_genes <- function(sex_label) {
      live <- results$featureselection[[sex_label]]$consensus_genes
      if (!is.null(live) && length(live) >= 2) {
        return(list(genes = live, is_live = TRUE,
                    note = sprintf("%d genes from this session's live %s consensus panel.", length(live), sex_label),
                    holdout_sample_ids = results$featureselection[[sex_label]]$holdout_sample_ids %||% character(0)))
      }
      if (isTRUE(dataset$is_bundled_reference)) {
        bundled <- read_table_safe(sprintf("FS_input_%s.csv", sex_label))
        if (!is.null(bundled) && nrow(bundled) >= 2 && "gene" %in% colnames(bundled)) {
          return(list(genes = unique(as.character(bundled$gene)), is_live = FALSE,
                      note = sprintf("%d genes from the bundled %s panel.", nrow(bundled), sex_label),
                      holdout_sample_ids = character(0)))
        }
      }
      list(genes = character(0), is_live = FALSE,
           note = sprintf("No live %s panel yet - run Feature Selection on the currently loaded dataset first.", sex_label),
           holdout_sample_ids = character(0))
    }

    own_panel_genes <- function(sex_label) {
      genes <- unique(trimws(unlist(strsplit(input$gene_list %||% "", "[,\n\t ]+"))))
      genes <- genes[nzchar(genes)]
      list(genes = genes, is_live = FALSE, note = sprintf("%d pasted genes.", length(genes)), holdout_sample_ids = character(0))
    }

    wgcna_panel_genes <- function(sex_label) {
      req(input$wgcna_module_pick)
      mg <- results$wgcna$module_genes
      validate(need(!is.null(mg) && input$wgcna_module_pick %in% names(mg),
                    "Run WGCNA Step 3 (Modules) first, then pick a module above."))
      genes <- unique(as.character(mg[[input$wgcna_module_pick]]))
      list(genes = genes, is_live = TRUE,
           note = sprintf("%d genes from WGCNA module \"%s\" (this session).", length(genes), input$wgcna_module_pick),
           holdout_sample_ids = character(0))
    }

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

    ext_meta_raw <- reactive({
      req(input$ext_meta_file)
      path <- input$ext_meta_file$datapath
      if (grepl("\\.rds$", input$ext_meta_file$name, ignore.case = TRUE)) {
        loaded <- safe_read_rds(path)
        validate(need(isTRUE(loaded$ok), loaded$error %||% "Could not read this .rds file."))
        d <- loaded$value
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
      r_fit <- diag_result_value(panel_sex)
      ext_models <- if (!is.null(r_fit) && length(r_fit$genes)) {
        tryCatch(diag_apply_models_external(r_fit, d$expr[, meta_sub$sample, drop = FALSE], y), error = function(e) NULL)
      } else NULL
      list(gr = gr, expr_sub = expr_sub, y = y, genes = genes, panel_note = cand$note,
           n_ref = sum(y == input$ext_ref_group), n_comp = sum(y == input$ext_comp_group),
           ref_group = input$ext_ref_group, comp_group = input$ext_comp_group,
           ext_models = ext_models, panel_sex = panel_sex)
    }, ignoreInit = TRUE)

    output$ext_model_ui <- renderUI({
      r <- ext_result(); req(r)
      if (is.null(r$ext_models)) {
        return(div(class = "empty-note", icon("circle-info"),
                   sprintf("Run Model Training for the %s panel first. The trained logistic regression, elastic net, random forest and SVM models are then scored on this external cohort exactly as trained: training Youden threshold reused, external samples z-scored within the external cohort, no refitting or re-tuning.", r$panel_sex)))
      }
      tagList(
        p(class = "submodule-desc", sprintf(
          "%d of %d panel genes are present in the external cohort (a missing gene enters as 0 after z-scoring). Models were fit on the %s training split and are applied unchanged; the decision threshold is the training Youden cutoff.",
          r$ext_models$n_genes_present, r$ext_models$n_genes_panel, r$panel_sex)),
        fluidRow(
          column(6, plotOutput(ns("ext_model_roc"), height = 360)),
          column(6, DT::dataTableOutput(ns("ext_model_table")))
        )
      )
    })

    output$ext_model_table <- DT::renderDataTable({
      r <- ext_result(); req(r$ext_models)
      rows <- lapply(r$ext_models$models, function(mm) {
        if (!isTRUE(mm$available)) {
          return(data.frame(model = mm$label, `AUC (95% CI)` = NA_character_, sensitivity = NA_real_, specificity = NA_real_,
                            PPV = NA_real_, NPV = NA_real_, accuracy = NA_real_, `PR-AUC` = NA_real_, Brier = NA_real_,
                            note = mm$reason, check.names = FALSE, stringsAsFactors = FALSE))
        }
        data.frame(model = mm$label, `AUC (95% CI)` = sprintf("%.3f (%.3f-%.3f)", mm$auc, mm$ci_lo, mm$ci_hi),
                   sensitivity = round(mm$perf$sensitivity, 3), specificity = round(mm$perf$specificity, 3),
                   PPV = round(mm$perf$precision, 3), NPV = round(mm$perf$npv, 3), accuracy = round(mm$perf$accuracy, 3),
                   `PR-AUC` = round(mm$perf$pr_auc, 3), Brier = round(mm$perf$brier, 3), note = "",
                   check.names = FALSE, stringsAsFactors = FALSE)
      })
      DT::datatable(do.call(rbind, rows), rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    output$ext_model_roc <- renderPlot({
      r <- ext_result(); req(r$ext_models)
      curves <- Filter(function(mm) isTRUE(mm$available), r$ext_models$models)
      validate(need(length(curves) > 0, "No trained model could be applied to this external cohort."))
      df <- do.call(rbind, lapply(curves, function(mm) {
        co <- pROC::coords(mm$roc, "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
        d <- data.frame(fpr = 1 - co$specificity, tpr = co$sensitivity, model = sprintf("%s (AUC %.2f)", mm$label, mm$auc))
        d[order(d$fpr, d$tpr), ]
      }))
      ggplot(df, aes(x = fpr, y = tpr, color = model)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
        geom_path(linewidth = 1) + coord_equal() +
        labs(x = "1 - Specificity", y = "Sensitivity", title = "External cohort ROC - trained models, no refitting", color = NULL) +
        theme_arthomix(base_size = 11)
    })

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
      df <- df[order(-df$auc), ]
      df
    })

    ext_gene_df_csv_safe <- reactive({
      df <- ext_gene_df()
      df$gene <- tx_csv_safe(df$gene)
      df
    })

    output$ext_gene_table <- DT::renderDataTable({
      req(ext_result())
      DT::datatable(ext_gene_df(), rownames = FALSE, width = "100%",
                    options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact") |>
        DT::formatStyle("hub", target = "row", backgroundColor = DT::styleEqual(c(TRUE, FALSE), c("#e6f4ea", "")))
    })

    output$ext_gene_download <- downloadHandler(
      filename = function() "external_validation_gene_auc.csv",
      content = function(file) write.csv(ext_gene_df_csv_safe(), file, row.names = FALSE)
    )

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
          width = NULL, title = "Trained classifiers applied to the external cohort (frozen models, no refitting)", status = "primary", solidHeader = FALSE,
          uiOutput(ns("ext_model_ui"))
        ),
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
      validate(need(length(cand$genes) >= 1, sprintf("No gene panel available for %s: %s", sex_label, cand$note)))
      genes <- intersect(cand$genes, rownames(dataset$expr))
      validate(need(length(genes) >= 2, sprintf(
        "Fewer than 2 %s genes from the chosen panel are present in the currently loaded expression matrix (panel: %s).",
        sex_label, cand$note
      )))

      expr_sub <- dataset$expr[genes, common, drop = FALSE]

      fit <- withProgress(
        message = sprintf("Fitting %s diagnostic models (logistic regression, elastic net, random forest, SVM)...", sex_label),
        value = 0.3,
        diag_fit_sex(expr_sub, y, params = diag_advanced_params(), holdout_ids = cand$holdout_sample_ids %||% character(0))
      )
      fit$candidate_note <- cand$note
      fit$n_ref <- sum(y == input$ref_group); fit$n_comp <- sum(y == input$comp_group)
      fit$ref_group <- input$ref_group; fit$comp_group <- input$comp_group

      # Automatic leakage-safe headline metric: when this run's Test-split AUC is NOT
      # leakage-safe (the default/bundled-panel path), compute the nested-CV AUC right
      # here - so it's ready the moment results render, no separate manual tab/click
      # needed - and make it the primary metric instead of the naive Test-split AUC.
      nested <- NULL
      if (!isTRUE(fit$leakage_safe)) {
        compute_nested_cv <- function() {
          cand_genes <- tryCatch(diag_leakagesafe_candidate_genes(sex_label), error = function(e) character(0))
          cand_genes <- intersect(cand_genes, rownames(dataset$expr))
          if (length(cand_genes) < 5) {
            return(list(
              pooled = list(available = FALSE, reason = sprintf(
                "Fewer than 5 %s WGCNA-candidate genes are available for automatic leakage-safe validation - run Candidate Gene Identification on this dataset for a leakage-safe headline metric.",
                sex_label)),
              per_fold = data.frame(fold = integer(0), n_panel = integer(0), auc = numeric(0)),
              outer_k = NA_integer_, n_folds_completed = 0
            ))
          }
          if (length(cand_genes) > FS_MAX_CANDIDATE_GENES) {
            v <- apply(dataset$expr[cand_genes, common, drop = FALSE], 1, stats::var)
            cand_genes <- names(sort(v, decreasing = TRUE))[seq_len(FS_MAX_CANDIDATE_GENES)]
          }
          expr_candidates <- dataset$expr[cand_genes, common, drop = FALSE]
          diag_validate_nested(expr_candidates, y, outer_k = 5)
        }
        nested <- tryCatch(
          withProgress(
            message = sprintf("Computing automatic leakage-safe nested-CV AUC for %s (this run's headline metric)...", sex_label),
            value = 0.7,
            compute_nested_cv()
          ),
          error = function(e) NULL
        )
      }
      fit <- diag_attach_headline(fit, nested)
      fit
    }

    diag_leakagesafe_candidate_genes <- function(sex_label) {
      if (identical(sex_label, "pooled")) {
        cand_final <- results$candidates$final
        if (!is.null(cand_final) && identical(cand_final$selection, "pooled") && length(cand_final$genes) >= 3) {
          return(unique(as.character(cand_final$genes)))
        }
        return(unique(c(results$candidates$female$genes, results$candidates$male$genes)))
      }
      unique(as.character(results$candidates[[sex_label]]$genes))
    }

    diag_build_leakagesafe_sex <- function(sex_label, sex_value, outer_k) {
      req(input$ref_group, input$comp_group)
      validate(need(input$ref_group != input$comp_group, "Reference and comparison group must be different."))

      cand_genes <- diag_leakagesafe_candidate_genes(sex_label)
      validate(need(length(cand_genes) >= 5, sprintf(
        "Leakage-safe validation needs this session's live Candidate Gene Identification output for %s (the WGCNA-candidate list Feature Selection started from) - run Candidate Gene Identification (and Feature Selection) on this dataset first.",
        sex_label
      )))

      meta <- dataset$meta
      sex_ok <- if (is.null(sex_value)) rep(TRUE, nrow(meta)) else (!is.na(meta$sex) & as.character(meta$sex) == sex_value)
      meta <- meta[sex_ok &
                     !is.na(meta$group) & as.character(meta$group) %in% c(input$ref_group, input$comp_group), , drop = FALSE]
      common <- intersect(colnames(dataset$expr), meta$sample)
      validate(need(length(common) >= 10, sprintf("Fewer than 10 %s samples match this contrast.", sex_label)))
      meta <- meta[match(common, meta$sample), , drop = FALSE]
      y_full <- factor(as.character(meta$group), levels = c(input$ref_group, input$comp_group))
      validate(need(all(table(y_full) >= 6), sprintf("Each group needs at least 6 %s samples.", sex_label)))

      genes <- intersect(cand_genes, rownames(dataset$expr))
      validate(need(length(genes) >= 5, sprintf(
        "Fewer than 5 %s WGCNA-candidate genes are present in the currently loaded expression matrix.", sex_label
      )))
      if (length(genes) > FS_MAX_CANDIDATE_GENES) {
        v <- apply(dataset$expr[genes, common, drop = FALSE], 1, stats::var)
        genes <- names(sort(v, decreasing = TRUE))[seq_len(FS_MAX_CANDIDATE_GENES)]
      }
      expr_candidates <- dataset$expr[genes, common, drop = FALSE]

      withProgress(
        message = sprintf("Running %s leakage-safe validation (Univariate + LASSO reselected per outer fold)...", sex_label),
        value = 0.3,
        diag_validate_nested(expr_candidates, y_full, outer_k = outer_k)
      )
    }

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

    diag_leakagesafe_female <- eventReactive(input$run_female_leakagesafe_btn, {
      diag_build_leakagesafe_sex("female", sex_levels()$female, input$female_leakagesafe_k %||% 5)
    }, ignoreInit = TRUE)
    diag_leakagesafe_male <- eventReactive(input$run_male_leakagesafe_btn, {
      diag_build_leakagesafe_sex("male", sex_levels()$male, input$male_leakagesafe_k %||% 5)
    }, ignoreInit = TRUE)
    diag_leakagesafe_pooled <- eventReactive(input$run_pooled_leakagesafe_btn, {
      diag_build_leakagesafe_sex("pooled", NULL, input$pooled_leakagesafe_k %||% 5)
    }, ignoreInit = TRUE)

    diag_leakagesafe_has_run <- reactiveValues(female = FALSE, male = FALSE, pooled = FALSE)
    observeEvent(input$run_female_leakagesafe_btn, diag_leakagesafe_has_run$female <- TRUE, ignoreInit = TRUE)
    observeEvent(input$run_male_leakagesafe_btn, diag_leakagesafe_has_run$male <- TRUE, ignoreInit = TRUE)
    observeEvent(input$run_pooled_leakagesafe_btn, diag_leakagesafe_has_run$pooled <- TRUE, ignoreInit = TRUE)
    observeEvent(dataset$source, {
      diag_leakagesafe_has_run$female <- FALSE; diag_leakagesafe_has_run$male <- FALSE; diag_leakagesafe_has_run$pooled <- FALSE
    }, ignoreInit = TRUE)

    lapply(c("female", "male", "pooled"), function(sex_label) {
      res_fn <- switch(sex_label, female = diag_leakagesafe_female, male = diag_leakagesafe_male, pooled = diag_leakagesafe_pooled)
      output[[paste0(sex_label, "_leakagesafe_result_ui")]] <- renderUI({
        if (!isTRUE(diag_leakagesafe_has_run[[sex_label]])) {
          return(div(class = "empty-note", icon("circle-info"), "Not run yet. Click Run above."))
        }
        r <- tryCatch(res_fn(), error = function(e) {
          msg <- conditionMessage(e)
          if (nzchar(msg)) msg else NULL
        })
        if (is.null(r)) return(NULL)
        if (is.character(r)) return(div(class = "empty-note", icon("triangle-exclamation"), r))
        pooled <- r$pooled
        if (!isTRUE(pooled$available)) {
          return(div(class = "empty-note", icon("triangle-exclamation"), pooled$reason %||% "Leakage-safe AUC unavailable."))
        }
        tagList(
          fluidRow(
            valueBox(sprintf("%.3f", pooled$auc),
                     sprintf("Leakage-safe pooled AUC (95%% CI %.3f-%.3f, n=%d)", pooled$ci_lo, pooled$ci_hi, pooled$n),
                     icon = icon("shield-halved"), color = "light-blue", width = 6),
            valueBox(sprintf("%d / %d", r$n_folds_completed, r$outer_k), "Outer folds completed",
                     icon = icon("layer-group"), color = "purple", width = 6)
          ),
          DT::dataTableOutput(ns(paste0(sex_label, "_leakagesafe_table")))
        )
      })
      output[[paste0(sex_label, "_leakagesafe_table")]] <- DT::renderDataTable({
        req(isTRUE(diag_leakagesafe_has_run[[sex_label]]))
        r <- tryCatch(res_fn(), error = function(e) NULL)
        req(r, is.list(r), !is.null(r$per_fold))
        DT::datatable(r$per_fold, rownames = FALSE, width = "100%",
                      options = list(pageLength = 10, dom = "t", scrollX = TRUE), class = "stripe hover compact")
      })
      outputOptions(output, paste0(sex_label, "_leakagesafe_table"), suspendWhenHidden = FALSE)
    })

    diag_result_error_msg <- function(sex_label) {
      if (!isTRUE(diag_valid[[sex_label]])) return(NULL)
      fr <- switch(sex_label, female = diag_result_female, male = diag_result_male, pooled = diag_result_pooled)
      tryCatch({ fr(); NULL }, error = function(e) {
        msg <- conditionMessage(e)
        if (nzchar(msg)) msg else NULL
      })
    }

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
          cv_aucs <- c(mean(r$lr$cv_auc, na.rm = TRUE), mean(r$enet$cv_auc, na.rm = TRUE), mean(r$rf$cv_auc, na.rm = TRUE), mean(r$svm$cv_auc, na.rm = TRUE))
          best_label <- c("Logistic Regression", "Elastic Net", "Random Forest", "SVM")[which.max(cv_aucs)]
          return(tags$li(icon("check", style = "color: #1a9c5f;"), strong(sprintf(" %s completed: ", sex)),
                          sprintf("best model %s (CV-AUC = %.3f)", best_label, max(cv_aucs, na.rm = TRUE))))
        }
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

    diag_result_line <- function(sex_label, r) {
      if (!is.null(r)) {
        return(tagList(
          p(strong("Result: "),
            sprintf("%d genes, %d samples (%d vs %d), %s vs %s → CV-AUC logistic regression %.3f / elastic net %.3f / random forest %.3f / SVM %.3f.",
                    r$n_input, r$n_samples, r$n_comp, r$n_ref, r$comp_group, r$ref_group,
                    mean(r$lr$cv_auc, na.rm = TRUE), mean(r$enet$cv_auc, na.rm = TRUE), mean(r$rf$cv_auc, na.rm = TRUE), mean(r$svm$cv_auc, na.rm = TRUE))),
          if (isTRUE(r$leakage_safe))
            p(class = "empty-note", icon("shield-halved"), "Leakage-safe: the Test split above is the held-out sample set this gene panel was never selected against.")
          else if (identical(r$headline_metric, "nested_cv"))
            p(class = "empty-note", icon("shield-halved"),
              sprintf(
                "Not leakage-safe via the Test split above (this gene panel was chosen using the full sample pool) - the leakage-safe headline metric is instead the automatic nested-CV AUC = %.3f (95%% CI %.3f-%.3f), computed for every run. See Model Testing (Internal) for details.",
                r$nested_cv$pooled$auc, r$nested_cv$pooled$ci_lo, r$nested_cv$pooled$ci_hi
              ))
          else
            p(class = "empty-note", icon("triangle-exclamation"), "Not leakage-safe: this gene panel's Test-set AUC is exploratory, not a validated estimate - either it came from the bundled/precomputed panel, or Feature Selection's held-out split was disabled or didn't overlap this cohort.")
        ))
      }
      err <- diag_result_error_msg(sex_label)
      if (!is.null(err)) p(strong("Result: "), span(style = "color: #c0392b;", sprintf("failed - %s", err))) else NULL
    }
    output$female_result_line <- renderUI({ diag_result_line("female", diag_result_value("female")) })
    output$male_result_line <- renderUI({ diag_result_line("male", diag_result_value("male")) })
    output$pooled_result_line <- renderUI({ diag_result_line("pooled", diag_result_value("pooled")) })

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

    build_eval_rows <- function(ev, label) {
      if (!isTRUE(ev$available)) {
        return(list(data.frame(dataset = label, metric = "Status", value = ev$reason, stringsAsFactors = FALSE)))
      }
      auc_str <- paste0(sprintf("%.3f (%.3f-%.3f)", ev$auc, ev$ci_lo, ev$ci_hi), diag_separation_note(ev))
      list(
        data.frame(dataset = label, metric = "AUC (95% CI)", value = auc_str, stringsAsFactors = FALSE),
        data.frame(dataset = label, metric = "PR-AUC", value = sprintf("%.3f", ev$perf$pr_auc), stringsAsFactors = FALSE),
        data.frame(dataset = label, metric = "n samples (n positive)", value = sprintf("%d (%d)", ev$n, ev$n_pos), stringsAsFactors = FALSE),
        data.frame(dataset = label, metric = "Sensitivity @ training cutoff", value = sprintf("%.3f", ev$perf$sensitivity), stringsAsFactors = FALSE),
        data.frame(dataset = label, metric = "Specificity @ training cutoff", value = sprintf("%.3f", ev$perf$specificity), stringsAsFactors = FALSE),
        data.frame(dataset = label, metric = "Precision (PPV) @ training cutoff", value = sprintf("%.3f", ev$perf$precision), stringsAsFactors = FALSE),
        data.frame(dataset = label, metric = "NPV @ training cutoff", value = sprintf("%.3f", ev$perf$npv), stringsAsFactors = FALSE),
        data.frame(dataset = label, metric = "F1 @ training cutoff", value = sprintf("%.3f", ev$perf$f1), stringsAsFactors = FALSE),
        data.frame(dataset = label, metric = "Balanced accuracy @ training cutoff", value = sprintf("%.3f", ev$perf$balanced_accuracy), stringsAsFactors = FALSE),
        data.frame(dataset = label, metric = "MCC @ training cutoff", value = sprintf("%.3f", ev$perf$mcc), stringsAsFactors = FALSE),
        data.frame(dataset = label, metric = "Accuracy @ training cutoff", value = sprintf("%.3f", ev$perf$accuracy), stringsAsFactors = FALSE),
        data.frame(dataset = label, metric = "Brier score", value = sprintf("%.3f", ev$perf$brier), stringsAsFactors = FALSE)
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

      output[[paste0(sex_label, "_headline_ui")]] <- renderUI({
        r <- res()
        if (is.null(r) || isTRUE(r$leakage_safe)) return(NULL)
        if (identical(r$headline_metric, "nested_cv")) {
          pooled <- r$nested_cv$pooled
          box(
            width = NULL, status = "primary", solidHeader = FALSE,
            title = sprintf("Leakage-safe headline metric - %s (computed automatically)", tools::toTitleCase(sex_label)),
            fluidRow(
              valueBox(sprintf("%.3f", pooled$auc),
                       sprintf("Nested-CV AUC (95%% CI %.3f-%.3f, n=%d) - PRIMARY, leakage-safe", pooled$ci_lo, pooled$ci_hi, pooled$n),
                       icon = icon("shield-halved"), color = "light-blue", width = 6),
              valueBox(sprintf("%d / %d", r$nested_cv$n_folds_completed, r$nested_cv$outer_k), "Outer folds completed (auto-run)",
                       icon = icon("layer-group"), color = "purple", width = 6)
            ),
            p(class = "empty-note", icon("circle-info"),
              "Computed automatically because this gene panel was chosen using the full sample pool (no live held-out split) - see the disclosure under Leakage-safe Validation for how this reselects the panel per outer fold. Each technique tab below still shows its own Test-split AUC, now marked exploratory/not leakage-safe.")
          )
        } else {
          div(class = "empty-note", style = "border-left: 3px solid #c0392b; padding-left: 8px;",
              icon("triangle-exclamation"),
              sprintf("Not leakage-safe, and the automatic nested-CV headline metric is unavailable: %s Each technique tab's Test-split AUC below is exploratory only.",
                      r$nested_cv$pooled$reason %||% "Reselect a live gene panel or provide more candidate genes."))
        }
      })

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

        output[[paste0(prefix, "_test_summary")]] <- renderUI({
          r <- res()
          if (is.null(r)) return(not_yet_note())
          rr <- r[[key]]
          leak_safe <- isTRUE(r$leakage_safe)
          test_title <- if (leak_safe) sprintf("Test-split AUC (n=%d)%s", rr$test$n, diag_separation_note(rr$test)) else
            sprintf("Test-split AUC (n=%d)%s - NOT leakage-safe, exploratory only", rr$test$n, diag_separation_note(rr$test))
          test_tile <- if (isTRUE(rr$test$available)) {
            valueBox(sprintf("%.3f", rr$test$auc), test_title,
                     icon = icon("flask"), color = if (leak_safe) "light-blue" else "orange", width = 6)
          } else {
            valueBox("N/A", "Test-split AUC", icon = icon("flask"), color = "red", width = 6)
          }
          n_tile <- valueBox(if (isTRUE(rr$test$available)) rr$test$n else "-", "Test samples held out", icon = icon("users"), color = "purple", width = 6)
          tagList(
            if (!leak_safe)
              div(class = "empty-note", style = "border-left: 3px solid #c0392b; padding-left: 8px; margin-bottom: 8px;",
                  icon("triangle-exclamation"),
                  "This model's Test-split AUC below evaluates a gene panel chosen using the full sample pool - exploratory only, not a validated estimate. See the leakage-safe headline metric above.")
            else NULL,
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
        output[[paste0(prefix, "_test_calibration_plot")]] <- renderPlot({
          r <- res(); if (is.null(r)) return(NULL)
          rr <- r[[key]]
          req(isTRUE(rr$test$available))
          cal <- diag_calibration(rr$test$y_eval, rr$test$pred, levels(rr$test$y_eval)[2])
          diag_plot_calibration(cal, sprintf("Test-split calibration - %s (%s)", label, tools::toTitleCase(sex_label)))
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

      register_generoc_outputs <- function(mode) {
        plot_id <- paste0(sex_label, "_", mode, "_generoc_plot")
        table_id <- paste0(sex_label, "_", mode, "_generoc_table")
        download_id <- paste0(sex_label, "_", mode, "_generoc_download")
        hub_download_id <- paste0(sex_label, "_", mode, "_hub_download")
        auc_thr_id <- paste0(sex_label, "_", mode, "_hub_auc_thr")
        p_thr_id <- paste0(sex_label, "_", mode, "_hub_p_thr")

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
          content = function(file) {
            df <- gene_auc_df(res()$gene_roc_train, res()$gene_roc_test)
            df$gene <- tx_csv_safe(df$gene)
            write.csv(df, file, row.names = FALSE)
          }
        )
        output[[hub_download_id]] <- downloadHandler(
          filename = function() sprintf("%s_hub_genes.csv", sex_label),
          content = function(file) {
            df <- gene_auc_df(res()$gene_roc_train, res()$gene_roc_test)
            df <- df[df$hub, c("gene", "train_auc", "train_p", "test_auc", "test_p")]
            df$gene <- tx_csv_safe(df$gene)
            write.csv(df, file, row.names = FALSE)
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
