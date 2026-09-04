## R/methylomics/13_Diagnostic_Classifier/mod_methyl_diagnostic.R
## Diagnostic Classifier submodule (script09_diagnostic_classifier). Methylomics only -
## transcriptomics' mod_diagnostic.R is a separate, unrelated tab.

mod_methyl_diagnostic_config <- list(
  id = "diagnostic", title = "Diagnostic Classifier", icon = "stethoscope", group = "Biomarker modeling",
  description = "Diagnostic classifiers (logistic regression, elastic net, SVM, random forest, XGBoost, kNN) on single CpGs or a panel, with cross-validated tuning."
)

dxm_beta_to_m <- function(beta) {
  b <- pmin(pmax(beta, 1e-6), 1 - 1e-6)
  log2(b / (1 - b))
}

dxm_parse_num_list <- function(txt, default) {
  if (is.null(txt) || !nzchar(trimws(txt %||% ""))) return(default)
  v <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",")[[1]])))
  v <- v[!is.na(v)]
  if (length(v) == 0) default else v
}

dxm_sex_label <- function(sex_sel) switch(sex_sel %||% "female", all = "All samples", female = "Female", male = "Male", tools::toTitleCase(as.character(sex_sel)))

dxm_normalize_sex <- function(x) {
  v <- toupper(substr(trimws(as.character(x)), 1, 1))
  ifelse(v %in% c("F", "M"), v, NA_character_)
}

dxm_load_wgcna_for_sex <- function(sex_sel) {
  if (!identical(sex_sel, "all")) return(tryCatch(load_default_wgcna_module_assignment(sex_sel), error = function(e) NULL))
  parts <- Filter(Negate(is.null), lapply(c("female", "male"), function(s) tryCatch(load_default_wgcna_module_assignment(s), error = function(e) NULL)))
  if (length(parts) == 0) NULL else do.call(rbind, parts)
}

dxm_load_fs_votes_for_sex <- function(sex_sel) {
  if (!identical(sex_sel, "all")) return(tryCatch(load_default_diagnostic_ensemble_votes(sex_sel), error = function(e) NULL))
  parts <- Filter(Negate(is.null), lapply(c("female", "male"), function(s) tryCatch(load_default_diagnostic_ensemble_votes(s), error = function(e) NULL)))
  if (length(parts) == 0) return(NULL)
  tbl <- do.call(rbind, parts)
  tbl[order(-tbl$n_votes), ][!duplicated(tbl$cpg), ]
}

DXM_POS <- "Class1"
DXM_NEG <- "Class0"

dxm_validate_checklist <- function(dxm) {
  rows <- list()
  add <- function(check, status, detail) rows[[length(rows) + 1]] <<- data.frame(Check = check, Status = status, Detail = detail, stringsAsFactors = FALSE)

  n_train <- nrow(dxm$train_X); n_test <- nrow(dxm$test_internal_X)
  add("Sample counts", "OK", sprintf("%d training samples, %d test samples", n_train, n_test))
  dup_train <- sum(duplicated(rownames(dxm$train_X)))
  add("Duplicated samples", if (dup_train == 0) "OK" else "WARN", sprintf("%d duplicate sample ID(s) in training data", dup_train))
  dup_cpg <- sum(duplicated(colnames(dxm$train_X)))
  add("Duplicated CpGs", if (dup_cpg == 0) "OK" else "WARN", sprintf("%d duplicate CpG column(s)", dup_cpg))
  n_missing <- sum(is.na(as.matrix(dxm$train_X)))
  add("Missing methylation values", if (n_missing == 0) "OK" else "WARN", sprintf("%d missing value(s) in training matrix", n_missing))
  non_numeric <- !all(vapply(dxm$train_X, is.numeric, logical(1)))
  add("Non-numeric methylation values", if (!non_numeric) "OK" else "FAIL", if (non_numeric) "One or more feature columns are not numeric" else "All feature columns numeric")
  const_cols <- vapply(dxm$train_X, function(x) length(unique(x[!is.na(x)])) <= 1, logical(1))
  add("Constant features", if (!any(const_cols)) "OK" else "WARN", sprintf("%d constant feature(s) in training data", sum(const_cols)))
  nzv <- tryCatch(caret::nearZeroVar(dxm$train_X), error = function(e) integer(0))
  add("Near-zero-variance features", if (length(nzv) == 0) "OK" else "WARN", sprintf("%d near-zero-variance feature(s)", length(nzv)))
  cls_tbl <- table(dxm$train_y)
  imb_ratio <- if (length(cls_tbl) == 2 && min(cls_tbl) > 0) max(cls_tbl) / min(cls_tbl) else NA_real_
  add("Class balance (training)", if (!is.na(imb_ratio) && imb_ratio <= 3) "OK" else "WARN",
      sprintf("%s=%d, %s=%d", dxm$ref_level, cls_tbl[[DXM_NEG]] %||% 0, dxm$comp_level, cls_tbl[[DXM_POS]] %||% 0))
  na_pheno <- sum(is.na(dxm$train_y)) + sum(is.na(dxm$test_internal_y))
  add("Missing phenotype labels", if (na_pheno == 0) "OK" else "FAIL", sprintf("%d sample(s) with an unmapped class label", na_pheno))
  shared <- intersect(colnames(dxm$train_X), colnames(dxm$test_internal_X))
  add("Train/test feature compatibility", if (length(shared) == ncol(dxm$train_X)) "OK" else "WARN",
      sprintf("%d of %d training features present in the test data", length(shared), ncol(dxm$train_X)))
  add("Minimum class size for cross-validation", if (min(cls_tbl) >= 10) "OK" else "WARN",
      sprintf("Smallest training class has %d sample(s)", min(cls_tbl)))
  do.call(rbind, rows)
}

dxm_intersect_features <- function(train_ids, test_ids) {
  shared <- train_ids[train_ids %in% test_ids]
  list(train = train_ids, test = test_ids, shared = shared, unmatched = setdiff(train_ids, test_ids))
}

dxm_smote_fold <- function(x, y) {
  x <- as.data.frame(x)
  tab <- table(y)
  if (length(tab) != 2 || min(tab) < 6) return(list(x = x, y = y))
  res <- tryCatch({
    sm <- smotefamily::SMOTE(X = x, target = as.character(y), K = min(5, min(tab) - 1))
    out <- sm$data
    list(x = out[, setdiff(colnames(out), "class"), drop = FALSE], y = factor(out$class, levels = levels(y)))
  }, error = function(e) list(x = x, y = y))
  res
}

dxm_cv_control <- function(input) {
  folds <- input$cv_folds %||% 10
  repeats <- input$cv_repeats %||% 1
  search <- if (identical(input$search_method, "random")) "random" else "grid"
  sampling <- switch(input$imbalance_mode %||% "none",
                      "weighted" = "up",
                      "smote" = list(name = "smote (train-fold only)", func = dxm_smote_fold, first = TRUE),
                      NULL)
  caret::trainControl(method = "repeatedcv", number = folds, repeats = repeats,
                       classProbs = TRUE, summaryFunction = caret::twoClassSummary,
                       savePredictions = "final", search = search, sampling = sampling)
}

dxm_fit_caret <- function(method, X, y, tune_grid, tune_length, ctrl, seed, preProcess = NULL, extra = list()) {
  set.seed(seed)
  args <- c(list(x = X, y = y, method = method, trControl = ctrl, metric = "ROC", preProcess = preProcess), extra)
  if (!is.null(tune_grid)) args$tuneGrid <- tune_grid else args$tuneLength <- tune_length
  list(model = do.call(caret::train, args), kind = "caret")
}

dxm_xgb_grid <- function(input, mid) {
  expand.grid(
    max_depth = dxm_parse_num_list(input[[paste0(mid, "_max_depth")]], c(2, 3, 4)),
    eta = dxm_parse_num_list(input[[paste0(mid, "_eta")]], c(0.05, 0.1, 0.3)),
    min_child_weight = dxm_parse_num_list(input[[paste0(mid, "_min_child_weight")]], c(1, 3)),
    subsample = dxm_parse_num_list(input[[paste0(mid, "_subsample")]], c(0.8, 1)),
    colsample_bytree = dxm_parse_num_list(input[[paste0(mid, "_colsample")]], c(0.8, 1)),
    gamma = dxm_parse_num_list(input[[paste0(mid, "_gamma")]], c(0))
  )
}

dxm_fit_xgb_native <- function(X, y, input, mid, ctrl, seed) {
  validate(need(!identical(input$imbalance_mode, "smote"),
                "SMOTE imbalance handling isn't available for Gradient Boosting/XGBoost (this model trains via its own native cross-validation path, not caret's fold-safe sampling hook) - choose \"None\" or \"Class weighting\" on the Filters & Parameters tab instead."))
  grid <- dxm_xgb_grid(input, mid)
  if (identical(ctrl$search, "random") && nrow(grid) > 12) grid <- grid[sample(nrow(grid), 12), , drop = FALSE]
  nrounds_max <- input[[paste0(mid, "_nrounds")]] %||% 200
  early_stop <- input[[paste0(mid, "_early_stop")]] %||% 20
  folds_n <- ctrl$number %||% 10
  scale_pos_weight <- if (identical(input$imbalance_mode, "weighted")) {
    tab <- table(y); as.numeric(tab[[DXM_NEG]] / tab[[DXM_POS]])
  } else NULL

  Xm <- data.matrix(X)
  yb <- as.integer(y == DXM_POS)
  dtrain <- xgboost::xgb.DMatrix(data = Xm, label = yb)
  best <- NULL
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    params <- list(objective = "binary:logistic", eval_metric = "auc",
                    max_depth = g$max_depth, eta = g$eta, min_child_weight = g$min_child_weight,
                    subsample = g$subsample, colsample_bytree = g$colsample_bytree, gamma = g$gamma)
    if (!is.null(scale_pos_weight)) params$scale_pos_weight <- scale_pos_weight
    set.seed(seed)
    cv <- tryCatch(xgboost::xgb.cv(params = params, data = dtrain, nrounds = nrounds_max, nfold = folds_n,
                                    stratified = TRUE, early_stopping_rounds = early_stop, verbose = 0, showsd = FALSE),
                    error = function(e) NULL)
    if (is.null(cv)) next
    best_iter <- cv$early_stop$best_iteration %||% cv$niter
    mean_auc <- cv$evaluation_log$test_auc_mean[best_iter]
    if (length(mean_auc) == 1 && !is.na(mean_auc) && (is.null(best) || mean_auc > best$mean_auc)) {
      best <- list(params = params, nrounds = best_iter, mean_auc = mean_auc, grid_row = g)
    }
  }
  validate(need(!is.null(best), "XGBoost cross-validated grid search failed for every combination tried - widen the grid or increase CV folds."))
  set.seed(seed)
  model <- xgboost::xgb.train(params = best$params, data = dtrain, nrounds = best$nrounds, verbose = 0)
  list(model = model, kind = "xgb", params = best$params, nrounds = best$nrounds, best = best)
}

dxm_predict_prob <- function(fit, X) {
  if (inherits(fit, "train")) {
    p <- stats::predict(fit, newdata = X, type = "prob")
    p[[DXM_POS]]
  } else if (inherits(fit, "xgb.Booster")) {
    stats::predict(fit, xgboost::xgb.DMatrix(data.matrix(X)))
  } else stop("Unsupported fit object")
}

dxm_roc_bundle <- function(y, prob) {
  r <- tryCatch(pROC::roc(response = y, predictor = as.numeric(prob), levels = c(DXM_NEG, DXM_POS), direction = "auto", quiet = TRUE),
                error = function(e) NULL)
  if (is.null(r)) return(NULL)
  ci <- tryCatch(as.numeric(pROC::ci.auc(r, method = "delong")),
                 error = function(e) tryCatch(as.numeric(pROC::ci.auc(r, method = "bootstrap", boot.n = 1000)), error = function(e) c(NA, NA, NA)))
  co <- pROC::coords(r, "all", ret = c("threshold", "specificity", "sensitivity"), transpose = FALSE)
  list(roc = r, auc = as.numeric(r$auc), ci_lo = ci[1], ci_hi = ci[3], n = length(y), coords = co)
}

dxm_cv_roc_from_fit <- function(fit) {
  pred <- fit$pred
  if (is.null(pred)) return(NULL)
  bt <- fit$bestTune
  for (nm in names(bt)) pred <- pred[pred[[nm]] == bt[[nm]], , drop = FALSE]
  fold_ids <- unique(pred$Resample)
  curves <- lapply(fold_ids, function(f) dxm_roc_bundle(pred$obs[pred$Resample == f], pred[[DXM_POS]][pred$Resample == f]))
  curves <- Filter(Negate(is.null), curves)
  overall <- dxm_roc_bundle(pred$obs, pred[[DXM_POS]])
  aucs <- vapply(curves, function(x) x$auc, numeric(1))
  list(folds = curves, overall = overall, pooled = data.frame(obs = pred$obs, prob = pred[[DXM_POS]]),
       mean_auc = mean(aucs), sd_auc = stats::sd(aucs), n_folds = length(curves))
}

dxm_xgb_cv_roc <- function(X, y, params, nrounds, folds_n, seed) {
  set.seed(seed)
  folds <- caret::createFolds(y, k = folds_n, list = TRUE)
  preds <- lapply(seq_along(folds), function(i) {
    te_idx <- folds[[i]]; tr_idx <- setdiff(seq_len(nrow(X)), te_idx)
    if (length(unique(y[tr_idx])) < 2) return(NULL)
    dtr <- xgboost::xgb.DMatrix(data.matrix(X[tr_idx, , drop = FALSE]), label = as.integer(y[tr_idx] == DXM_POS))
    dte <- xgboost::xgb.DMatrix(data.matrix(X[te_idx, , drop = FALSE]))
    m <- xgboost::xgb.train(params = params, data = dtr, nrounds = nrounds, verbose = 0)
    data.frame(obs = y[te_idx], prob = stats::predict(m, dte), Resample = paste0("Fold", i))
  })
  preds <- Filter(Negate(is.null), preds)
  pred <- do.call(rbind, preds)
  curves <- lapply(unique(pred$Resample), function(f) dxm_roc_bundle(pred$obs[pred$Resample == f], pred$prob[pred$Resample == f]))
  curves <- Filter(Negate(is.null), curves)
  overall <- dxm_roc_bundle(pred$obs, pred$prob)
  aucs <- vapply(curves, function(x) x$auc, numeric(1))
  list(folds = curves, overall = overall, pooled = data.frame(obs = pred$obs, prob = pred$prob),
       mean_auc = mean(aucs), sd_auc = stats::sd(aucs), n_folds = length(curves))
}

dxm_pick_threshold <- function(strategy, roc_bundle) {
  if (is.null(roc_bundle)) return(0.5)
  co <- roc_bundle$coords
  pick_first_meeting <- function(cond) {
    ok <- cond & !is.infinite(co$threshold)
    if (!any(ok)) return(0.5)
    as.numeric(co$threshold[which.max(ok)])
  }
  switch(strategy,
    "youden" = {
      j <- co$sensitivity + co$specificity - 1
      as.numeric(co$threshold[which.max(j)])
    },
    "sensitivity" = pick_first_meeting(co$sensitivity >= 0.9),
    "specificity" = pick_first_meeting(co$specificity >= 0.9),
    0.5
  )
}

dxm_confusion <- function(y, prob, threshold) {
  pred <- factor(ifelse(prob >= threshold, DXM_POS, DXM_NEG), levels = levels(y))
  tab <- table(Predicted = pred, Actual = y)
  TP <- tab[DXM_POS, DXM_POS]; TN <- tab[DXM_NEG, DXM_NEG]; FP <- tab[DXM_POS, DXM_NEG]; FN <- tab[DXM_NEG, DXM_POS]
  sens <- TP / (TP + FN); spec <- TN / (TN + FP); prec <- TP / (TP + FP); npv <- TN / (TN + FN)
  f1 <- 2 * prec * sens / (prec + sens); acc <- (TP + TN) / sum(tab); bal_acc <- (sens + spec) / 2
  mcc_den <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- if (is.na(mcc_den) || mcc_den == 0) NA_real_ else (TP * TN - FP * FN) / mcc_den
  list(table = tab, TP = TP, TN = TN, FP = FP, FN = FN, sensitivity = sens, specificity = spec,
       precision = prec, npv = npv, f1 = f1, accuracy = acc, balanced_accuracy = bal_acc, mcc = mcc)
}

dxm_metrics_bundle <- function(y, prob, threshold, roc_bundle) {
  conf <- dxm_confusion(y, prob, threshold)
  y_bin <- as.integer(y == DXM_POS)
  brier <- mean((prob - y_bin)^2)
  pr_auc <- tryCatch(PRROC::pr.curve(scores.class0 = prob[y_bin == 1], scores.class1 = prob[y_bin == 0], curve = FALSE)$auc.integral,
                      error = function(e) NA_real_)
  c(conf, list(auc = roc_bundle$auc %||% NA_real_, auc_ci_lo = roc_bundle$ci_lo %||% NA_real_, auc_ci_hi = roc_bundle$ci_hi %||% NA_real_,
               brier = brier, pr_auc = pr_auc, n = length(y)))
}

dxm_calibration <- function(y, prob, bins = 10) {
  y_bin <- as.integer(y == DXM_POS)
  brks <- unique(seq(0, 1, length.out = bins + 1))
  bin_id <- cut(prob, brks, include.lowest = TRUE)
  df <- data.frame(prob = prob, y = y_bin, bin = bin_id)
  agg <- stats::aggregate(cbind(mean_pred = prob, mean_obs = y) ~ bin, data = df, FUN = mean)
  n_bin <- as.data.frame(table(bin = df$bin)); names(n_bin) <- c("bin", "n")
  agg <- merge(agg, n_bin, by = "bin", all.x = TRUE)
  fit <- tryCatch(stats::glm(y ~ prob, family = stats::binomial(), data = df), error = function(e) NULL)
  list(table = agg, brier = mean((prob - y_bin)^2),
       slope = if (!is.null(fit)) unname(stats::coef(fit)[2]) else NA_real_,
       intercept = if (!is.null(fit)) unname(stats::coef(fit)[1]) else NA_real_)
}

dxm_overfitting_note <- function(train_auc, cv_auc, test_auc, test_label = "independent test") {
  if (is.na(train_auc) || is.na(test_auc)) return("Not enough completed evaluations yet to assess overfitting.")
  gap <- train_auc - test_auc
  if (gap > 0.15) {
    sprintf("Training discrimination (AUC = %.3f) is substantially higher than %s discrimination (AUC = %.3f), suggesting possible overfitting or limited generalization to new samples.", train_auc, test_label, test_auc)
  } else if (!is.na(cv_auc) && (train_auc - cv_auc) > 0.15) {
    sprintf("Training discrimination (AUC = %.3f) is substantially higher than cross-validated discrimination (AUC = %.3f), suggesting the model may be fitting noise in the training split.", train_auc, cv_auc)
  } else {
    sprintf("Training (AUC = %.3f) and %s (AUC = %.3f) discrimination are broadly consistent, with no strong evidence of overfitting in this comparison.", train_auc, test_label, test_auc)
  }
}

DXM_MAX_CANDIDATE_CPGS <- 200

## Leakage-safe nested-CV validator (no held-out split): reselects the panel inside every outer fold, training-fold labels only.
dxm_validate_nested <- function(X_full, y_full, outer_k = 5, uni_top_n = 100, lasso_alpha = 1, seed = 42) {
  y_full <- droplevels(factor(as.character(y_full), levels = c(DXM_NEG, DXM_POS)))
  validate(need(nlevels(y_full) == 2, "Leakage-safe nested-CV validation needs exactly two classes."))
  n <- nrow(X_full)

  nf <- max(2, min(outer_k, min(table(y_full))))
  set.seed(seed)
  folds <- caret::createFolds(y_full, k = nf, list = TRUE)

  pred_oof <- rep(NA_real_, n)
  per_fold <- list()

  for (fi in seq_along(folds)) {
    te <- folds[[fi]]; tr <- setdiff(seq_len(n), te)
    y_tr <- droplevels(y_full[tr]); y_te <- droplevels(y_full[te])
    if (nlevels(y_tr) < 2 || nlevels(y_te) < 2) next

    m_tr <- t(as.matrix(X_full[tr, , drop = FALSE]))  # CpG rows x sample cols, training fold only

    design <- stats::model.matrix(~y_tr)
    uni_fit <- tryCatch(limma::eBayes(limma::lmFit(m_tr, design)), error = function(e) NULL)
    if (is.null(uni_fit)) next
    tt <- tryCatch(limma::topTable(uni_fit, coef = 2, number = Inf, sort.by = "P"), error = function(e) NULL)
    if (is.null(tt) || nrow(tt) < 2) next
    uni_cpgs <- rownames(tt)[seq_len(min(uni_top_n, nrow(tt)))]

    Xtr_raw <- as.matrix(X_full[tr, uni_cpgs, drop = FALSE])
    nf_lasso <- max(2, min(5, min(table(y_tr))))
    cv <- tryCatch(glmnet::cv.glmnet(Xtr_raw, y_tr, family = "binomial", alpha = lasso_alpha, nfolds = nf_lasso),
                   error = function(e) NULL)
    panel <- character(0)
    if (!is.null(cv)) {
      co <- stats::coef(cv, s = "lambda.min")[-1, 1, drop = TRUE]
      panel <- names(co)[co != 0]
    }
    if (length(panel) < 1) panel <- utils::head(uni_cpgs, min(5, length(uni_cpgs)))

    fit_df <- data.frame(y = y_tr, X_full[tr, panel, drop = FALSE], check.names = FALSE)
    model <- tryCatch(suppressWarnings(stats::glm(y ~ ., data = fit_df, family = stats::binomial)), error = function(e) NULL)
    if (is.null(model)) next
    newdata_te <- data.frame(X_full[te, panel, drop = FALSE], check.names = FALSE)
    pred <- tryCatch(as.numeric(stats::predict(model, newdata = newdata_te, type = "response")), error = function(e) NULL)
    if (is.null(pred)) next
    pred_oof[te] <- pred

    rb <- dxm_roc_bundle(y_te, pred)
    auc_i <- if (!is.null(rb)) rb$auc else NA_real_
    per_fold[[length(per_fold) + 1]] <- data.frame(fold = fi, n_panel = length(panel), auc = round(auc_i, 3))
  }

  per_fold_df <- if (length(per_fold) > 0) do.call(rbind, per_fold) else data.frame(fold = integer(0), n_panel = integer(0), auc = numeric(0))
  have_oof <- !is.na(pred_oof)
  pooled <- list(available = FALSE, reason = "Not enough folds completed to score a pooled AUC.")
  if (sum(have_oof) >= 10 && length(unique(y_full[have_oof])) == 2) {
    rb_pooled <- dxm_roc_bundle(y_full[have_oof], pred_oof[have_oof])
    if (!is.null(rb_pooled)) {
      pooled <- list(available = TRUE, auc = rb_pooled$auc, ci_lo = rb_pooled$ci_lo, ci_hi = rb_pooled$ci_hi, n = sum(have_oof))
    }
  }
  list(pooled = pooled, per_fold = per_fold_df, n_folds_completed = nrow(per_fold_df), outer_k = nf)
}

## Headline AUC: naive Test AUC if leakage_safe, else nested-CV AUC when available (Test AUC still shown, just demoted).
dxm_attach_headline <- function(leakage_safe, nested_cv = NULL) {
  if (isTRUE(leakage_safe)) {
    return(list(headline_metric = "test_split", nested_cv = NULL))
  }
  metric <- if (!is.null(nested_cv) && isTRUE(nested_cv$pooled$available)) "nested_cv" else "test_split"
  list(headline_metric = metric, nested_cv = nested_cv)
}

dxm_learning_curve <- function(X, y, fit_one, fracs, seed) {
  out <- lapply(fracs, function(fr) {
    set.seed(seed)
    idx <- tryCatch(caret::createDataPartition(y, p = fr, list = FALSE)[, 1], error = function(e) NULL)
    if (is.null(idx) || length(idx) < 10) return(NULL)
    Xs <- X[idx, , drop = FALSE]; ys <- y[idx]
    if (length(unique(ys)) < 2 || min(table(ys)) < 3) return(NULL)
    res <- tryCatch(fit_one(Xs, ys), error = function(e) NULL)
    if (is.null(res)) return(NULL)
    data.frame(frac = fr, n = length(ys), train_auc = res$train_auc, cv_auc = res$cv_auc)
  })
  do.call(rbind, Filter(Negate(is.null), out))
}

dxm_theme <- function() ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom")

dxm_plot_roc <- function(bundles_named, title, subtitle = NULL) {
  dfs <- lapply(names(bundles_named), function(nm) {
    b <- bundles_named[[nm]]; if (is.null(b)) return(NULL)
    co <- b$coords
    data.frame(fpr = 1 - co$specificity, tpr = co$sensitivity, series = sprintf("%s (AUC = %.3f, n = %d)", nm, b$auc, b$n))
  })
  df <- do.call(rbind, Filter(Negate(is.null), dfs))
  validate(need(!is.null(df) && nrow(df) > 0, "Not enough data to draw an ROC curve yet - run the model and evaluate test data first."))
  df <- df[order(df$series, df$fpr), ]
  ggplot2::ggplot(df, ggplot2::aes(x = fpr, y = tpr, color = series)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_step(linewidth = 1) + ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)",
                  title = title, subtitle = subtitle, color = NULL) + dxm_theme() +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 8)) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 3, byrow = TRUE))
}

dxm_plot_cv_roc <- function(cv_roc, title) {
  validate(need(!is.null(cv_roc) && length(cv_roc$folds) > 0, "Not enough data to draw a cross-validated ROC curve yet."))
  fold_df <- do.call(rbind, lapply(seq_along(cv_roc$folds), function(i) {
    co <- cv_roc$folds[[i]]$coords
    data.frame(fpr = 1 - co$specificity, tpr = co$sensitivity, fold = paste0("Fold ", i))
  }))
  mean_co <- cv_roc$overall$coords
  mean_df <- data.frame(fpr = 1 - mean_co$specificity, tpr = mean_co$sensitivity)
  ggplot2::ggplot() +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_step(data = fold_df[order(fold_df$fold, fold_df$fpr), ], ggplot2::aes(x = fpr, y = tpr, group = fold), color = "grey70", alpha = 0.6) +
    ggplot2::geom_step(data = mean_df[order(mean_df$fpr), ], ggplot2::aes(x = fpr, y = tpr), color = "#2563EB", linewidth = 1.2) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)", title = title,
                  subtitle = sprintf("Mean CV AUC = %.3f +/- %.3f (SD) across %d fold(s); pooled out-of-fold AUC = %.3f",
                                      cv_roc$mean_auc, cv_roc$sd_auc, cv_roc$n_folds, cv_roc$overall$auc)) +
    dxm_theme() + ggplot2::theme(legend.position = "none")
}

dxm_plot_calibration <- function(cal, title) {
  validate(need(!is.null(cal), "Not enough data to draw a calibration curve yet."))
  ggplot2::ggplot(cal$table, ggplot2::aes(x = mean_pred, y = mean_obs)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_line(color = "#2563EB") + ggplot2::geom_point(ggplot2::aes(size = n), color = "#2563EB") +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::labs(x = "Mean predicted probability", y = "Observed proportion", title = title,
                  subtitle = sprintf("Brier score = %.3f", cal$brier), size = "Bin n") + dxm_theme()
}

dxm_plot_learning_curve <- function(lc, title) {
  validate(need(!is.null(lc) && nrow(lc) > 0, "Not enough data to draw a learning curve yet."))
  df <- rbind(data.frame(n = lc$n, score = lc$train_auc, series = "Training AUC"),
              data.frame(n = lc$n, score = lc$cv_auc, series = "Cross-validated AUC"))
  ggplot2::ggplot(df, ggplot2::aes(x = n, y = score, color = series)) + ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::ylim(0, 1) + ggplot2::labs(x = "Training sample size", y = "ROC-AUC", title = title, color = NULL) + dxm_theme()
}

dxm_plot_roc_compare <- function(bundles_named, title) {
  dfs <- lapply(names(bundles_named), function(nm) {
    b <- bundles_named[[nm]]; if (is.null(b)) return(NULL)
    co <- b$coords
    data.frame(fpr = 1 - co$specificity, tpr = co$sensitivity, series = sprintf("%s (AUC = %.3f)", nm, b$auc))
  })
  df <- do.call(rbind, Filter(Negate(is.null), dfs))
  validate(need(!is.null(df) && nrow(df) > 0, "Select at least one completed run to compare."))
  df <- df[order(df$series, df$fpr), ]
  ggplot2::ggplot(df, ggplot2::aes(x = fpr, y = tpr, color = series)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_step(linewidth = 1) + ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)", title = title, color = NULL) +
    dxm_theme() +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 8)) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 3, byrow = TRUE))
}

DXM_MODEL_SPECS <- list(
  lr = list(id = "lr", label = "Logistic Regression", icon = "chart-line", kind = "caret",
    params_ui = function(ns, mid) tagList(
      radioButtons(ns(paste0(mid, "_penalty")), "Penalty", choices = c("L2 (ridge)" = "0", "L1 (lasso)" = "1"), selected = "0", inline = TRUE),
      textInput(ns(paste0(mid, "_C")), "C (inverse regularization strength, comma-separated)", value = "0.01, 0.1, 1, 10"),
      numericInput(ns(paste0(mid, "_maxit")), "Maximum iterations", value = 1e5, min = 100, step = 1000),
      numericInput(ns(paste0(mid, "_tol")), "Tolerance", value = 1e-7, min = 1e-10, step = 1e-7),
      helpText("Single CpG mode fits a plain (unregularized) logistic regression instead - glmnet requires at least two predictor columns.")
    ),
    fit = function(X, y, input, mid, ctrl, seed) {
      if (ncol(X) < 2) return(dxm_fit_caret("glm", X, y, NULL, NULL, ctrl, seed))
      alpha_fixed <- as.numeric(input[[paste0(mid, "_penalty")]] %||% "0")
      C <- dxm_parse_num_list(input[[paste0(mid, "_C")]], c(0.01, 0.1, 1, 10))
      grid <- data.frame(alpha = alpha_fixed, lambda = sort(1 / C))
      dxm_fit_caret("glmnet", X, y, grid, NULL, ctrl, seed,
                     extra = list(thresh = input[[paste0(mid, "_tol")]] %||% 1e-7, maxit = input[[paste0(mid, "_maxit")]] %||% 1e5))
    }),
  enet = list(id = "enet", label = "Elastic Net", icon = "sliders-h", kind = "caret",
    params_ui = function(ns, mid) tagList(
      checkboxInput(ns(paste0(mid, "_auto")), "Automatic hyperparameter search", value = TRUE),
      conditionalPanel(sprintf("!input['%s']", ns(paste0(mid, "_auto"))),
        textInput(ns(paste0(mid, "_alpha")), "alpha / L1 ratio (comma-separated)", value = "0, 0.2, 0.4, 0.6, 0.8, 1"),
        textInput(ns(paste0(mid, "_lambda")), "lambda (comma-separated)", value = "0.0001, 0.001, 0.01, 0.1, 1")
      ),
      numericInput(ns(paste0(mid, "_tune_length")), "Search size (when automatic)", value = 10, min = 3, max = 30),
      numericInput(ns(paste0(mid, "_maxit")), "Maximum iterations", value = 1e5, min = 100, step = 1000),
      numericInput(ns(paste0(mid, "_tol")), "Tolerance", value = 1e-7, min = 1e-10, step = 1e-7),
      helpText("Single CpG mode fits a plain (unregularized) logistic regression instead - glmnet requires at least two predictor columns.")
    ),
    fit = function(X, y, input, mid, ctrl, seed) {
      if (ncol(X) < 2) return(dxm_fit_caret("glm", X, y, NULL, NULL, ctrl, seed))
      extra <- list(thresh = input[[paste0(mid, "_tol")]] %||% 1e-7, maxit = input[[paste0(mid, "_maxit")]] %||% 1e5)
      if (isTRUE(input[[paste0(mid, "_auto")]])) {
        dxm_fit_caret("glmnet", X, y, NULL, input[[paste0(mid, "_tune_length")]] %||% 10, ctrl, seed, extra = extra)
      } else {
        grid <- expand.grid(alpha = dxm_parse_num_list(input[[paste0(mid, "_alpha")]], seq(0, 1, 0.2)),
                             lambda = dxm_parse_num_list(input[[paste0(mid, "_lambda")]], 10^seq(-4, 0, length.out = 5)))
        dxm_fit_caret("glmnet", X, y, grid, NULL, ctrl, seed, extra = extra)
      }
    }),
  svm = list(id = "svm", label = "Support Vector Machine", icon = "vector-square", kind = "caret",
    params_ui = function(ns, mid) tagList(
      radioButtons(ns(paste0(mid, "_kernel")), "Kernel", choices = c("Linear" = "linear", "RBF (radial)" = "radial", "Polynomial" = "poly"), selected = "radial", inline = TRUE),
      uiOutput(ns(paste0(mid, "_kernel_params")))
    ),
    fit = function(X, y, input, mid, ctrl, seed) {
      k <- input[[paste0(mid, "_kernel")]] %||% "radial"
      C <- dxm_parse_num_list(input[[paste0(mid, "_C")]], c(0.25, 0.5, 1, 2, 4))
      method <- switch(k, linear = "svmLinear", radial = "svmRadial", poly = "svmPoly")
      grid <- switch(k,
        linear = data.frame(C = C),
        radial = expand.grid(sigma = dxm_parse_num_list(input[[paste0(mid, "_sigma")]], c(0.01, 0.05, 0.1)), C = C),
        poly = expand.grid(degree = dxm_parse_num_list(input[[paste0(mid, "_degree")]], c(2, 3)),
                            scale = dxm_parse_num_list(input[[paste0(mid, "_scale")]], c(0.01, 0.1)), C = C)
      )
      dxm_fit_caret(method, X, y, grid, NULL, ctrl, seed, preProcess = c("center", "scale"))
    }),
  rf = list(id = "rf", label = "Random Forest", icon = "tree", kind = "caret",
    params_ui = function(ns, mid) tagList(
      numericInput(ns(paste0(mid, "_ntree")), "Number of trees", value = 500, min = 50, max = 5000, step = 50),
      textInput(ns(paste0(mid, "_mtry")), "mtry (comma-separated; blank = 1..p)", value = ""),
      numericInput(ns(paste0(mid, "_nodesize")), "Minimum terminal node size (0 = default)", value = 0, min = 0, step = 1),
      numericInput(ns(paste0(mid, "_maxnodes")), "Maximum terminal nodes (0 = unlimited)", value = 0, min = 0, step = 1)
    ),
    fit = function(X, y, input, mid, ctrl, seed) {
      mtry_vals <- dxm_parse_num_list(input[[paste0(mid, "_mtry")]], seq_len(ncol(X)))
      grid <- data.frame(mtry = round(mtry_vals))
      extra <- list(ntree = input[[paste0(mid, "_ntree")]] %||% 500)
      ns_val <- input[[paste0(mid, "_nodesize")]] %||% 0; if (isTRUE(ns_val > 0)) extra$nodesize <- ns_val
      mn_val <- input[[paste0(mid, "_maxnodes")]] %||% 0; if (isTRUE(mn_val > 0)) extra$maxnodes <- mn_val
      dxm_fit_caret("rf", X, y, grid, NULL, ctrl, seed, extra = extra)
    }),
  gbm = list(id = "gbm", label = "Gradient Boosting / XGBoost", icon = "bolt", kind = "xgb",
    params_ui = function(ns, mid) tagList(
      numericInput(ns(paste0(mid, "_nrounds")), "Max boosting rounds (n_estimators)", value = 200, min = 10, max = 2000, step = 10),
      numericInput(ns(paste0(mid, "_early_stop")), "Early stopping rounds", value = 20, min = 0, max = 200, step = 5),
      textInput(ns(paste0(mid, "_eta")), "learning_rate / eta (comma-separated)", value = "0.05, 0.1, 0.3"),
      textInput(ns(paste0(mid, "_max_depth")), "max_depth (comma-separated)", value = "2, 3, 4"),
      textInput(ns(paste0(mid, "_min_child_weight")), "min_child_weight (comma-separated)", value = "1, 3"),
      textInput(ns(paste0(mid, "_subsample")), "subsample (comma-separated)", value = "0.8, 1"),
      textInput(ns(paste0(mid, "_colsample")), "colsample_bytree (comma-separated)", value = "0.8, 1"),
      textInput(ns(paste0(mid, "_gamma")), "gamma (comma-separated)", value = "0")
    ),
    fit = dxm_fit_xgb_native),
  knn = list(id = "knn", label = "k-Nearest Neighbors", icon = "circle-nodes", kind = "caret",
    params_ui = function(ns, mid) tagList(
      textInput(ns(paste0(mid, "_k")), "Number of neighbors k (comma-separated)", value = "3, 5, 7, 9, 11, 15, 21"),
      helpText("Weights fixed to uniform and distance metric fixed to Euclidean: the optional kknn package (needed for distance-weighted kNN / Minkowski p / leaf size) isn't installed in this deployment.")
    ),
    fit = function(X, y, input, mid, ctrl, seed) {
      grid <- data.frame(k = round(dxm_parse_num_list(input[[paste0(mid, "_k")]], seq(3, 21, 2))))
      dxm_fit_caret("knn", X, y, grid, NULL, ctrl, seed, preProcess = c("center", "scale"))
    })
)

dxm_metrics_display <- function(m) {
  data.frame(
    Metric = c("Accuracy", "Balanced accuracy", "Sensitivity (recall)", "Specificity", "Precision", "NPV", "F1",
               "MCC", "ROC-AUC", "AUC 95% CI", "PR-AUC", "Brier score", "N"),
    Value = c(sprintf("%.3f", m$accuracy), sprintf("%.3f", m$balanced_accuracy), sprintf("%.3f", m$sensitivity),
              sprintf("%.3f", m$specificity), sprintf("%.3f", m$precision), sprintf("%.3f", m$npv), sprintf("%.3f", m$f1),
              sprintf("%.3f", m$mcc), sprintf("%.3f", m$auc), sprintf("%.3f - %.3f", m$auc_ci_lo, m$auc_ci_hi),
              sprintf("%.3f", m$pr_auc), sprintf("%.3f", m$brier), as.character(m$n)),
    stringsAsFactors = FALSE
  )
}

dxm_render_model_panel <- function(mid, spec, ns, input, dxm, feat, ms, headline = NULL) {
  if (!isTRUE(dxm$validated)) return(p(class = "text-muted", "Validate your data on the Datasets tab first."))
  single <- identical(input$analysis_type, "single")

  setup_box <- box(width = 12, status = "primary", solidHeader = TRUE, title = "Setup",
    tags$ul(
      tags$li(sprintf("Feature source: %s", switch(feat$source %||% "none",
        wgcna = "Methylomics WGCNA", featureselection = "Methylomics Feature Selection",
        uploaded = "Uploaded panel", manual = "Manual selection", "not loaded yet"))),
      tags$li(sprintf("Analysis type: %s", if (single) sprintf("Single CpG (%s)", input$single_cpg %||% "none selected")
                       else sprintf("Combined panel (%d CpG%s)", length(feat$selected), if (length(feat$selected) == 1) "" else "s"))),
      tags$li(sprintf("Phenotype: %s (reference) vs %s (comparison/positive class)", dxm$ref_level, dxm$comp_level)),
      tags$li(sprintf("Test data: %s", if (!is.null(ms$test_internal_prob)) "evaluated" else "not yet evaluated"))
    ))

  params_box <- box(width = 12, status = "primary", solidHeader = TRUE, title = "Parameters",
    spec$params_ui(ns, mid), br(),
    actionButton(ns(paste0(mid, "_run_btn")), if (single) "Run Single-CpG Diagnostic Analysis" else "Run Model", icon = icon("play"), class = "btn-primary"))

  if (is.null(ms$fit)) return(tagList(setup_box, params_box))

  results_box <- box(width = 12, status = "success", solidHeader = TRUE, title = "Results: Training & Cross-Validation",
    fluidRow(
      valueBox(sprintf("%.3f", ms$train_roc$auc), "Training AUC", icon = icon("chart-line"), color = "blue", width = 3),
      valueBox(sprintf("%.3f", ms$cv_roc$mean_auc), sprintf("Mean CV AUC (+/- %.3f SD, %d folds)", ms$cv_roc$sd_auc, ms$cv_roc$n_folds), icon = icon("layer-group"), color = "blue", width = 3),
      valueBox(sprintf("%.3f", ms$threshold), "Classification threshold", icon = icon("ruler"), color = "light-blue", width = 3),
      valueBox(length(ms$feature_ids), "Feature(s) used", icon = icon("dna"), color = "light-blue", width = 3)
    ),
    DT::dataTableOutput(ns(paste0(mid, "_train_metrics_table"))),
    br(), actionButton(ns(paste0(mid, "_test_btn")), "Run Test Evaluation", icon = icon("play"), class = "btn-primary"))

  out <- list(setup_box, params_box, results_box)

  if (!is.null(ms$test_internal_prob)) {
    leak_safe <- isTRUE(dxm$leakage_safe)
    headline <- headline %||% list(headline_metric = "test_split", nested_cv = NULL)
    is_primary_nested <- !leak_safe && identical(headline$headline_metric, "nested_cv")

    if (!leak_safe) {
      out <- c(out, list(
        if (is_primary_nested) {
          pooled <- headline$nested_cv$pooled
          box(width = 12, status = "primary", solidHeader = TRUE, title = "Leakage-safe headline metric (computed automatically)",
            fluidRow(
              valueBox(sprintf("%.3f", pooled$auc),
                       sprintf("Nested-CV AUC (95%% CI %.3f-%.3f, n=%d) - PRIMARY, leakage-safe", pooled$ci_lo, pooled$ci_hi, pooled$n),
                       icon = icon("shield-halved"), color = "light-blue", width = 6),
              valueBox(sprintf("%d / %d", headline$nested_cv$n_folds_completed, headline$nested_cv$outer_k),
                       "Outer folds completed (auto-run)", icon = icon("layer-group"), color = "purple", width = 6)
            ),
            p(class = "submodule-desc", icon("circle-info"),
              "Computed automatically because this run has no confirmed leakage-safe held-out split. Reselects the CpG panel with limma (moderated t) + LASSO inside every outer fold, on this cohort's candidate CpGs - see \"Results: Test Internal Data\" below for the exploratory, not-leakage-safe Test AUC."))
        } else {
          div(class = "empty-note", style = "border-left: 3px solid #c0392b; padding-left: 8px;", icon("triangle-exclamation"),
              sprintf("Not leakage-safe, and the automatic nested-CV headline metric is unavailable: %s The Test AUC below is exploratory only.",
                      headline$nested_cv$pooled$reason %||% "Not enough candidate CpGs or completed folds."))
        }
      ))
    }

    test_box_title <- if (leak_safe) "Results: Test Internal Data" else "Results: Test Internal Data (NOT leakage-safe - exploratory only)"
    test_auc_label <- if (leak_safe) sprintf("Test AUC (95%% CI %.3f-%.3f)", ms$test_internal_metrics$auc_ci_lo, ms$test_internal_metrics$auc_ci_hi) else
      sprintf("Test AUC (95%% CI %.3f-%.3f) - not leakage-safe", ms$test_internal_metrics$auc_ci_lo, ms$test_internal_metrics$auc_ci_hi)
    out <- c(out, list(box(width = 12, status = if (leak_safe) "success" else "warning", solidHeader = TRUE, title = test_box_title,
      if (!leak_safe)
        div(class = "empty-note", style = "border-left: 3px solid #c0392b; padding-left: 8px; margin-bottom: 8px;", icon("triangle-exclamation"),
            "This panel's Test AUC below was not confirmed to be selected without seeing these test samples - exploratory only, not a validated estimate. See the leakage-safe headline metric above.")
      else NULL,
      fluidRow(
        valueBox(sprintf("%.3f", ms$test_internal_metrics$auc), test_auc_label, icon = icon("chart-area"), color = if (leak_safe) "blue" else "orange", width = 3),
        valueBox(sprintf("%.3f", ms$test_internal_metrics$sensitivity), "Sensitivity", icon = icon("check"), color = "light-blue", width = 3),
        valueBox(sprintf("%.3f", ms$test_internal_metrics$specificity), "Specificity", icon = icon("shield"), color = "light-blue", width = 3),
        valueBox(ms$test_internal_metrics$n, "Test N", icon = icon("users"), color = "teal", width = 3)
      ),
      h5("Test Internal Data"), DT::dataTableOutput(ns(paste0(mid, "_test_metrics_table"))),
      hr(), p(class = "submodule-desc", dxm_overfitting_note(ms$train_metrics$auc, ms$cv_roc$mean_auc, ms$test_internal_metrics$auc, test_label = "internal validation")))))

    out <- c(out, list(box(width = 12, status = "primary", solidHeader = TRUE, title = "ROC / AUC",
      actionButton(ns(paste0(mid, "_roc_btn")), "Generate ROC/AUC", icon = icon("chart-area"), class = "btn-primary"),
      downloadButton(ns(paste0(mid, "_roc_png")), "Download ROC plot (PNG)", class = "btn-sm"),
      br(), br(),
      if (isTRUE(ms$roc_generated)) fluidRow(
        column(6, plotOutput(ns(paste0(mid, "_roc_traincv_plot")))),
        column(6, plotOutput(ns(paste0(mid, "_roc_cv_plot"))))
      ),
      if (isTRUE(ms$roc_generated)) fluidRow(
        column(6, plotOutput(ns(paste0(mid, "_roc_test_plot"))))
      ))))

    out <- c(out, list(box(width = 12, status = "primary", solidHeader = TRUE, title = "Diagnostics",
      fluidRow(
        column(6, h5("Training Confusion Matrix"), DT::dataTableOutput(ns(paste0(mid, "_confusion_train_table")))),
        column(6, h5("Test Confusion Matrix"), DT::dataTableOutput(ns(paste0(mid, "_confusion_test_table"))))
      ),
      hr(),
      actionButton(ns(paste0(mid, "_calib_btn")), "Generate Calibration", icon = icon("chart-line"), class = "btn-primary btn-sm"),
      if (isTRUE(ms$calib_generated)) plotOutput(ns(paste0(mid, "_calib_plot"))),
      hr(),
      actionButton(ns(paste0(mid, "_lc_btn")), "Generate Learning Curve", icon = icon("chart-line"), class = "btn-primary btn-sm"),
      if (isTRUE(ms$lc_generated)) plotOutput(ns(paste0(mid, "_lc_plot"))))))
  }
  tagList(out)
}

dxm_do_run_model <- function(mid, spec, input, dxm, feat, ms) {
  ids <- intersect(dxm_active_ids(mid, input, feat), colnames(dxm$train_X))
  validate(need(length(ids) > 0, "No selected CpG(s) are present in the training data - pick a CpG (single mode) or load a feature panel (combined mode) on the Filters & Parameters / Feature Source tabs first."))
  Xtr <- dxm$train_X[, ids, drop = FALSE]
  ytr <- dxm$train_y
  ctrl <- dxm_cv_control(input)
  seed <- input$dxm_seed %||% 42

  fit <- tryCatch(spec$fit(Xtr, ytr, input, mid, ctrl, seed), error = function(e) e)
  if (inherits(fit, "error")) {
    showNotification(paste0(spec$label, ": model fit failed - ", conditionMessage(fit)), type = "error", duration = 10)
    return(invisible(FALSE))
  }

  ms$fit <- fit; ms$analysis_type <- input$analysis_type; ms$feature_ids <- ids
  ms$ran_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  prob_train <- dxm_predict_prob(fit$model, Xtr)
  ms$train_prob <- prob_train
  ms$train_roc <- dxm_roc_bundle(ytr, prob_train)
  ms$cv_roc <- if (identical(fit$kind, "caret")) dxm_cv_roc_from_fit(fit$model) else dxm_xgb_cv_roc(Xtr, ytr, fit$params, fit$nrounds, ctrl$number %||% 10, seed)
  ms$threshold <- dxm_pick_threshold(input$threshold_strategy %||% "default", ms$train_roc)
  ms$train_metrics <- dxm_metrics_bundle(ytr, prob_train, ms$threshold, ms$train_roc)
  ms$confusion_train <- dxm_confusion(ytr, prob_train, ms$threshold)

  ms$test_internal_prob <- NULL; ms$test_internal_roc <- NULL; ms$test_internal_metrics <- NULL; ms$confusion_test <- NULL
  ms$roc_generated <- FALSE; ms$calib <- NULL; ms$calib_generated <- FALSE; ms$lc <- NULL; ms$lc_generated <- FALSE

  showNotification(sprintf("%s: trained on %d feature(s), %d training samples (train AUC = %.3f, mean CV AUC = %.3f). See the %s tab for full results.",
                            spec$label, length(ids), nrow(Xtr), ms$train_roc$auc, ms$cv_roc$mean_auc, spec$label), type = "message", duration = 8)
  invisible(TRUE)
}

dxm_active_ids <- function(mid, input, feat) {
  if (identical(input$analysis_type, "single")) {
    if (is.null(input$single_cpg) || !nzchar(input$single_cpg)) return(character(0))
    input$single_cpg
  } else feat$selected
}

dxm_register_model_server <- function(mid, spec, input, output, session, ns, dxm, feat, ms, runs, results, dxm_headline = NULL) {

  if (identical(mid, "svm")) {
    output[[paste0(mid, "_kernel_params")]] <- renderUI({
      k <- input[[paste0(mid, "_kernel")]] %||% "radial"
      tagList(
        textInput(ns(paste0(mid, "_C")), "C (comma-separated)", value = "0.25, 0.5, 1, 2, 4"),
        if (identical(k, "radial")) textInput(ns(paste0(mid, "_sigma")), "Gamma / sigma (comma-separated, RBF only)", value = "0.01, 0.05, 0.1"),
        if (identical(k, "poly")) tagList(
          textInput(ns(paste0(mid, "_degree")), "Degree (comma-separated, polynomial only)", value = "2, 3"),
          textInput(ns(paste0(mid, "_scale")), "Scale (comma-separated, polynomial only)", value = "0.01, 0.1")
        )
      )
    })
  }

  observeEvent(input[[paste0(mid, "_run_btn")]], {
    req(dxm$validated)
    dxm_do_run_model(mid, spec, input, dxm, feat, ms)
  })

  observeEvent(input[[paste0(mid, "_test_btn")]], {
    req(ms$fit)
    ids <- ms$feature_ids
    Xte <- dxm$test_internal_X[, ids, drop = FALSE]; yte <- dxm$test_internal_y
    ms$test_internal_prob <- dxm_predict_prob(ms$fit$model, Xte)
    ms$test_internal_roc <- dxm_roc_bundle(yte, ms$test_internal_prob)
    ms$test_internal_metrics <- dxm_metrics_bundle(yte, ms$test_internal_prob, ms$threshold, ms$test_internal_roc)
    ms$confusion_test <- dxm_confusion(yte, ms$test_internal_prob, ms$threshold)

    key <- paste(mid, ms$analysis_type, paste(ids, collapse = "|"))
    runs[[key]] <- list(model_id = mid, label = spec$label, analysis_type = ms$analysis_type, feature_ids = ids,
                         threshold = ms$threshold, train_roc = ms$train_roc, cv_roc = ms$cv_roc,
                         test_internal_roc = ms$test_internal_roc, ran_at = ms$ran_at)

    if (!is.null(results)) {
      results$diagnostic <- list(last_model = spec$label, analysis_type = ms$analysis_type, n_features = length(ids),
                                  stratum = dxm$sex, mode = dxm$mode,
                                  train_auc = round(ms$train_roc$auc, 3), cv_auc = round(ms$cv_roc$mean_auc, 3),
                                  test_auc = round(ms$test_internal_roc$auc, 3))

      model_entry <- list(
        model_id = mid, label = spec$label, kind = ms$fit$kind,
        fit = ms$fit, feature_ids = ids, threshold = ms$threshold,
        analysis_type = ms$analysis_type,
        ref_level = dxm$ref_level, comp_level = dxm$comp_level,
        sex_stratum = dxm$sex, mode = dxm$mode, seed = input$dxm_seed %||% 42,
        train_sample_ids = rownames(dxm$train_X),
        train_X = dxm$train_X[, ids, drop = FALSE],
        train_n = nrow(dxm$train_X), train_class_table = table(dxm$train_y),
        train_metrics = ms$train_metrics, train_roc = ms$train_roc, confusion_train = ms$confusion_train,
        cv_roc = ms$cv_roc,
        test_internal_metrics = ms$test_internal_metrics,
        ran_at = ms$ran_at
      )
      all_models <- results$diagnostic_models %||% list()
      all_models[[key]] <- model_entry
      results$diagnostic_models <- all_models
    }
    showNotification(sprintf("%s: test evaluation complete (n=%d, test AUC = %.3f).", spec$label, nrow(Xte), ms$test_internal_roc$auc), type = "message")
  })

  observeEvent(input[[paste0(mid, "_roc_btn")]], { req(ms$test_internal_roc); ms$roc_generated <- TRUE })

  observeEvent(input[[paste0(mid, "_calib_btn")]], {
    req(ms$cv_roc$pooled)
    ms$calib <- dxm_calibration(ms$cv_roc$pooled$obs, ms$cv_roc$pooled$prob)
    ms$calib_generated <- TRUE
  })

  observeEvent(input[[paste0(mid, "_lc_btn")]], {
    req(ms$fit)
    ids <- ms$feature_ids
    Xtr <- dxm$train_X[, ids, drop = FALSE]; ytr <- dxm$train_y
    ctrl5 <- caret::trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = caret::twoClassSummary, savePredictions = "final")
    seed <- input$dxm_seed %||% 42
    fit_one <- function(Xs, ys) {
      if (length(unique(ys)) < 2) return(NULL)
      if (identical(spec$kind, "xgb")) {
        m <- xgboost::xgb.train(params = ms$fit$params, data = xgboost::xgb.DMatrix(data.matrix(Xs), label = as.integer(ys == DXM_POS)), nrounds = ms$fit$nrounds, verbose = 0)
        p_tr <- dxm_predict_prob(m, Xs)
        cvr <- dxm_xgb_cv_roc(Xs, ys, ms$fit$params, ms$fit$nrounds, 5, seed)
      } else {
        f <- spec$fit(Xs, ys, input, mid, ctrl5, seed)
        p_tr <- dxm_predict_prob(f$model, Xs)
        cvr <- dxm_cv_roc_from_fit(f$model)
      }
      list(train_auc = dxm_roc_bundle(ys, p_tr)$auc, cv_auc = cvr$mean_auc)
    }
    ms$lc <- dxm_learning_curve(Xtr, ytr, fit_one, c(0.4, 0.6, 0.8, 1.0), seed)
    ms$lc_generated <- TRUE
  })

  output[[paste0(mid, "_train_metrics_table")]] <- DT::renderDataTable({
    req(ms$train_metrics); DT::datatable(dxm_metrics_display(ms$train_metrics), rownames = FALSE, options = list(dom = "t", paging = FALSE))
  })
  output[[paste0(mid, "_test_metrics_table")]] <- DT::renderDataTable({
    req(ms$test_internal_metrics); DT::datatable(dxm_metrics_display(ms$test_internal_metrics), rownames = FALSE, options = list(dom = "t", paging = FALSE))
  })
  output[[paste0(mid, "_confusion_train_table")]] <- DT::renderDataTable({
    req(ms$confusion_train); DT::datatable(as.data.frame.matrix(ms$confusion_train$table), options = list(dom = "t", paging = FALSE))
  })
  output[[paste0(mid, "_confusion_test_table")]] <- DT::renderDataTable({
    req(ms$confusion_test); DT::datatable(as.data.frame.matrix(ms$confusion_test$table), options = list(dom = "t", paging = FALSE))
  })

  output[[paste0(mid, "_roc_traincv_plot")]] <- renderPlot({
    req(ms$roc_generated); p <- dxm_plot_roc(list(Training = ms$train_roc), sprintf("%s - Training ROC", spec$label)); ms$last_roc_plot <- p; p
  })
  output[[paste0(mid, "_roc_cv_plot")]] <- renderPlot({
    req(ms$roc_generated, ms$cv_roc); dxm_plot_cv_roc(ms$cv_roc, sprintf("%s - Cross-Validated ROC", spec$label))
  })
  output[[paste0(mid, "_roc_test_plot")]] <- renderPlot({
    req(ms$roc_generated, ms$test_internal_roc); dxm_plot_roc(list(Test = ms$test_internal_roc), sprintf("%s - Test ROC", spec$label))
  })
  output[[paste0(mid, "_calib_plot")]] <- renderPlot({
    req(ms$calib_generated, ms$calib); p <- dxm_plot_calibration(ms$calib, sprintf("%s - Calibration (cross-validated predictions)", spec$label)); ms$last_calib_plot <- p; p
  })
  output[[paste0(mid, "_lc_plot")]] <- renderPlot({
    req(ms$lc_generated, ms$lc); p <- dxm_plot_learning_curve(ms$lc, sprintf("%s - Learning Curve", spec$label)); ms$last_lc_plot <- p; p
  })

  output[[paste0(mid, "_roc_png")]] <- downloadHandler(
    filename = function() sprintf("methylomics_diagnostic_%s_roc.png", mid),
    content = function(file) { req(ms$last_roc_plot); ggplot2::ggsave(file, ms$last_roc_plot, width = 7, height = 6, dpi = 150) }
  )

  output[[paste0("panel_", mid)]] <- renderUI({
    dxm_render_model_panel(mid, spec, ns, input, dxm, feat, ms,
                            headline = if (!is.null(dxm_headline)) dxm_headline() else NULL)
  })

  invisible(NULL)
}

mod_methyl_diagnostic_ui <- function(id) {
  ns <- NS(id)
  model_tabs <- lapply(DXM_MODEL_SPECS, function(spec) {
    tabPanel(spec$label, br(), withSpinner(uiOutput(ns(paste0("panel_", spec$id))), color = "#2563EB", type = 6))
  })
  do.call(tabsetPanel, c(
    list(id = ns("main_tabs"), type = "tabs",
      tabPanel("Datasets", br(), withSpinner(uiOutput(ns("setup_ui")), color = "#2563EB", type = 6)),
      tabPanel("Feature Source", br(), withSpinner(uiOutput(ns("feature_ui")), color = "#2563EB", type = 6)),
      tabPanel("Filters & Parameters", br(), withSpinner(uiOutput(ns("filters_ui")), color = "#2563EB", type = 6))),
    unname(model_tabs),
    list(
      tabPanel("Model Comparison", br(), withSpinner(uiOutput(ns("compare_ui")), color = "#2563EB", type = 6)),
      tabPanel("Test Internal Data", br(), withSpinner(uiOutput(ns("testdata_ui")), color = "#2563EB", type = 6)),
      tabPanel("Export", br(), withSpinner(uiOutput(ns("export_ui")), color = "#2563EB", type = 6))
    )
  ))
}

mod_methyl_diagnostic_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    dxm <- reactiveValues(train_X = NULL, train_y = NULL, test_internal_X = NULL, test_internal_y = NULL,
                           full_X = NULL, full_y = NULL, all_cpgs = character(0),
                           ref_level = "Control", comp_level = "RA", sex = "female", mode = "preloaded",
                           validated = FALSE, validation_report = NULL, leakage_safe = FALSE)
    feat <- reactiveValues(source = NULL, table = NULL, selected = character(0))

    dxm_apply_holdout_split <- function(holdout_ids) {
      req(dxm$full_X, dxm$full_y)
      avail <- rownames(dxm$full_X)
      test_ids <- intersect(as.character(holdout_ids), avail)
      if (length(test_ids) < 4 || (length(avail) - length(test_ids)) < 4) return(FALSE)
      test_rows <- avail %in% test_ids
      if (length(unique(dxm$full_y[!test_rows])) < 2 || length(unique(dxm$full_y[test_rows])) < 2) return(FALSE)
      dxm$train_X <- dxm$full_X[!test_rows, , drop = FALSE]; dxm$train_y <- dxm$full_y[!test_rows]
      dxm$test_internal_X <- dxm$full_X[test_rows, , drop = FALSE]; dxm$test_internal_y <- dxm$full_y[test_rows]
      dxm$leakage_safe <- TRUE
      TRUE
    }
    runs <- reactiveValues()
    compare_state <- reactiveValues(generated = FALSE, bundles = NULL)

    output$setup_ui <- renderUI({
      tagList(
        box(width = 12, status = "primary", solidHeader = TRUE, title = "Data source",
          radioButtons(ns("data_mode"), NULL,
            choices = c("Preloaded whole-blood dataset (sex-stratified RA)" = "preloaded", "Upload your own dataset" = "upload"),
            selected = if (METH_DIAG_DATA_AVAILABLE) "preloaded" else "upload"),
          conditionalPanel(sprintf("input['%s'] == 'preloaded'", ns("data_mode")),
            if (METH_DIAG_DATA_AVAILABLE) tagList(
              radioButtons(ns("sex_stratum"), "Sex stratum", choices = c("All samples" = "all", "Female" = "female", "Male" = "male"), selected = "female", inline = TRUE),
              helpText("Female/Male match script09's own sex-stratified panels; All samples combines both.")
            ) else p(class = "text-muted", "The preloaded diagnostic train/test data isn't available in this deployment - switch to Upload Data instead.")
          ),
          conditionalPanel(sprintf("input['%s'] == 'upload'", ns("data_mode")),
            fileInput(ns("upload_matrix"), "Methylation matrix (CSV/TSV: CpG rows x sample columns)"),
            radioButtons(ns("upload_scale"), "Input scale", choices = c("Beta-value" = "beta", "M-value" = "m"), selected = "beta", inline = TRUE),
            fileInput(ns("upload_sheet"), "Sample sheet (CSV/TSV: one row per sample)"),
            uiOutput(ns("upload_col_ui")),
            radioButtons(ns("upload_sex_stratum"), "Sex stratum", choices = c("All samples" = "all", "Female" = "female", "Male" = "male"), selected = "all", inline = TRUE),
            uiOutput(ns("upload_sex_col_ui"))
          ),
          fluidRow(
            column(4, textInput(ns("ref_level"), "Reference (control) class label", value = "Control")),
            column(4, textInput(ns("comp_level"), "Comparison (disease) class label", value = "RA")),
            column(4, numericInput(ns("train_frac"), "Training fraction", value = 0.75, min = 0.5, max = 0.9, step = 0.05))
          ),
          numericInput(ns("dxm_seed"), "Random seed", value = 42, min = 1, step = 1),
          actionButton(ns("validate_btn"), "Validate Data", icon = icon("vial"), class = "btn-primary")
        ),
        if (isTRUE(dxm$validated)) box(width = 12, status = if (any(dxm$validation_report$Status == "FAIL")) "danger" else "primary",
          solidHeader = TRUE, title = "Validation summary",
          DT::dataTableOutput(ns("validation_table")), hr(),
          p(class = "submodule-desc", sprintf("Active dataset: %s. %d training samples, %d test samples, %d candidate CpG(s).",
              if (identical(dxm$mode, "preloaded")) sprintf("Preloaded whole-blood dataset (%s)", dxm_sex_label(dxm$sex)) else sprintf("Uploaded dataset (%s)", dxm_sex_label(dxm$sex)),
              nrow(dxm$train_X), nrow(dxm$test_internal_X),
              length(dxm$all_cpgs))))
      )
    })

    output$upload_col_ui <- renderUI({
      req(input$upload_sheet)
      ps <- methyl_parse_sample_sheet(input$upload_sheet$datapath, input$upload_sheet$name)
      if (!isTRUE(ps$ok)) return(p(class = "text-danger", ps$error))
      cols <- colnames(ps$df)
      guess <- cols[which(grepl("class|phenotype|group|status|disease", cols, ignore.case = TRUE))[1]]
      selectInput(ns("upload_pheno_col"), "Phenotype / class column", choices = cols, selected = guess %||% cols[1])
    })

    output$upload_sex_col_ui <- renderUI({
      req(input$upload_sheet)
      ps <- methyl_parse_sample_sheet(input$upload_sheet$datapath, input$upload_sheet$name)
      if (!isTRUE(ps$ok)) return(NULL)
      cols <- colnames(ps$df)
      guess <- cols[which(grepl("sex|gender", cols, ignore.case = TRUE))[1]]
      selectInput(ns("upload_sex_col"), "Sex column (used when Sex stratum is Female/Male)", choices = cols, selected = guess %||% cols[1])
    })

    observeEvent(input$validate_btn, {
      mode <- input$data_mode %||% "preloaded"
      ref_lab <- trimws(input$ref_level %||% "Control"); comp_lab <- trimws(input$comp_level %||% "RA")
      validate(need(nzchar(ref_lab) && nzchar(comp_lab) && !identical(ref_lab, comp_lab), "Enter two distinct, non-empty reference and comparison class labels."))
      seed <- input$dxm_seed %||% 42
      train_frac <- input$train_frac %||% 0.75

      if (identical(mode, "preloaded")) {
        validate(need(METH_DIAG_DATA_AVAILABLE, "The preloaded diagnostic train/test data isn't available in this deployment - switch to Upload Data instead."))
        dd <- load_default_diagnostic_train_test()
        validate(need(!is.null(dd), "Could not load the preloaded diagnostic train/test data."))
        sex_sel <- input$sex_stratum %||% "female"
        sex_code <- switch(sex_sel, male = "M", female = "F", NA_character_)
        internal <- dd$internal
        sex_keep <- if (is.na(sex_code)) rep(TRUE, nrow(internal$pheno)) else internal$pheno$sex == sex_code
        keep <- sex_keep & internal$pheno$group %in% c(ref_lab, comp_lab)
        validate(need(sum(keep) > 20, "Fewer than 20 samples match this sex stratum and class labels - check the reference/comparison class labels (default: Control / RA)."))
        beta_sub <- internal$beta[, keep, drop = FALSE]
        pheno_sub <- internal$pheno[keep, , drop = FALSE]
        y_all <- factor(ifelse(pheno_sub$group == comp_lab, DXM_POS, DXM_NEG), levels = c(DXM_NEG, DXM_POS))
        set.seed(seed)
        train_idx <- caret::createDataPartition(y_all, p = train_frac, list = FALSE)[, 1]
        Xm <- as.data.frame(t(dxm_beta_to_m(beta_sub))); rownames(Xm) <- pheno_sub$gsm
        dxm$train_X <- Xm[train_idx, , drop = FALSE]; dxm$train_y <- y_all[train_idx]
        dxm$test_internal_X <- Xm[-train_idx, , drop = FALSE]; dxm$test_internal_y <- y_all[-train_idx]
        dxm$full_X <- Xm; dxm$full_y <- y_all
        dxm$all_cpgs <- colnames(Xm)
        dxm$mode <- "preloaded"; dxm$sex <- sex_sel
      } else {
        req(input$upload_matrix, input$upload_sheet)
        pm <- methyl_parse_matrix(input$upload_matrix$datapath, input$upload_matrix$name)
        validate(need(isTRUE(pm$ok), pm$error %||% "Could not parse the uploaded methylation matrix."))
        ps <- methyl_parse_sample_sheet(input$upload_sheet$datapath, input$upload_sheet$name)
        validate(need(isTRUE(ps$ok), ps$error %||% "Could not parse the uploaded sample sheet."))
        mat <- pm$mat; sheet <- ps$df
        sample_ids <- methyl_sheet_sample_ids(sheet, colnames(mat))
        common <- intersect(colnames(mat), sample_ids)
        validate(need(length(common) >= 10, "Fewer than 10 samples matched between the uploaded matrix and sample sheet."))
        mat <- mat[, common, drop = FALSE]; sheet <- sheet[match(common, sample_ids), , drop = FALSE]
        pheno_col <- input$upload_pheno_col
        validate(need(!is.null(pheno_col) && nzchar(pheno_col) && pheno_col %in% colnames(sheet), "Select a phenotype/class column from the uploaded sample sheet."))
        grp_raw <- trimws(as.character(sheet[[pheno_col]]))
        validate(need(all(c(ref_lab, comp_lab) %in% grp_raw), "The chosen reference/comparison labels were not both found in the phenotype column - check spelling/case."))
        sex_sel <- input$upload_sex_stratum %||% "all"
        sex_keep <- rep(TRUE, nrow(sheet))
        if (!identical(sex_sel, "all")) {
          sex_col <- input$upload_sex_col
          validate(need(!is.null(sex_col) && nzchar(sex_col) && sex_col %in% colnames(sheet),
                        "Select a sex column from the uploaded sample sheet to filter by Female/Male, or choose \"All samples\"."))
          sex_norm <- dxm_normalize_sex(sheet[[sex_col]])
          target <- if (identical(sex_sel, "male")) "M" else "F"
          sex_keep <- !is.na(sex_norm) & sex_norm == target
          validate(need(sum(sex_keep) > 0, sprintf("No samples matched sex = %s in the selected sex column.", dxm_sex_label(sex_sel))))
        }
        keep <- sex_keep & (grp_raw %in% c(ref_lab, comp_lab))
        mat <- mat[, keep, drop = FALSE]; grp_raw <- grp_raw[keep]
        validate(need(sum(duplicated(rownames(mat))) == 0, "Uploaded matrix has duplicated CpG IDs."))
        m_vals <- if (identical(input$upload_scale, "m")) mat else dxm_beta_to_m(mat)
        y_all <- factor(ifelse(grp_raw == comp_lab, DXM_POS, DXM_NEG), levels = c(DXM_NEG, DXM_POS))
        validate(need(min(table(y_all)) >= 6, "Fewer than 6 samples in the smaller class - not enough to fit or validate a classifier."))
        set.seed(seed)
        train_idx <- caret::createDataPartition(y_all, p = train_frac, list = FALSE)[, 1]
        Xm <- as.data.frame(t(m_vals))
        dxm$train_X <- Xm[train_idx, , drop = FALSE]; dxm$train_y <- y_all[train_idx]
        dxm$test_internal_X <- Xm[-train_idx, , drop = FALSE]; dxm$test_internal_y <- y_all[-train_idx]
        dxm$full_X <- Xm; dxm$full_y <- y_all
        dxm$all_cpgs <- colnames(Xm)
        dxm$mode <- "upload"; dxm$sex <- sex_sel
      }

      dxm$ref_level <- ref_lab; dxm$comp_level <- comp_lab
      dxm$validation_report <- dxm_validate_checklist(dxm)
      dxm$validated <- TRUE
      dxm$leakage_safe <- FALSE
      feat$source <- NULL; feat$table <- NULL; feat$selected <- character(0)
      showNotification("Data validated - see the Validation summary below, then continue to Feature Source.", type = "message")
    })

    output$validation_table <- DT::renderDataTable({
      req(dxm$validation_report)
      DT::datatable(dxm$validation_report, rownames = FALSE, options = list(dom = "t", paging = FALSE))
    })

    output$feature_ui <- renderUI({
      req(dxm$validated)
      tagList(
        box(width = 12, status = "primary", solidHeader = TRUE, title = "Feature source",
          radioButtons(ns("feature_source"), NULL,
            choices = c("Methylomics WGCNA" = "wgcna", "Methylomics Feature Selection" = "featureselection",
                        "Uploaded/Preloaded Feature Set" = "uploaded", "Manual CpG Selection" = "manual"),
            selected = "featureselection", inline = TRUE),
          conditionalPanel(sprintf("input['%s'] == 'wgcna'", ns("feature_source")),
            if (identical(dxm$mode, "preloaded")) {
              mt <- dxm_load_wgcna_for_sex(dxm$sex)
              mod_col <- if (!is.null(mt)) intersect(c("module", "Module", "color"), colnames(mt))[1] else NA_character_
              mod_choices <- if (!is.null(mt) && !is.na(mod_col)) c("All modules" = "__all__", stats::setNames(sort(unique(as.character(mt[[mod_col]]))), sort(unique(as.character(mt[[mod_col]]))))) else c("All modules" = "__all__")
              tagList(selectInput(ns("wgcna_module"), "Module", choices = mod_choices, selected = "__all__"),
                      numericInput(ns("wgcna_top_n"), "Number of top CpGs (0 = all)", value = 0, min = 0, step = 1))
            } else fileInput(ns("wgcna_upload"), "WGCNA module-assignment/hub-CpG CSV (from the WGCNA module's own download)"),
            actionButton(ns("wgcna_load_btn"), "Load WGCNA Features", icon = icon("download"), class = "btn-primary btn-sm")
          ),
          conditionalPanel(sprintf("input['%s'] == 'featureselection'", ns("feature_source")),
            if (identical(dxm$mode, "preloaded")) tagList(
              numericInput(ns("fs_min_votes"), "Minimum selection-method votes (of 3)", value = 2, min = 1, max = 3, step = 1),
              numericInput(ns("fs_top_n"), "Number of top CpGs (0 = all)", value = 0, min = 0, step = 1)
            ) else fileInput(ns("fs_upload"), "Feature Selection RDS export (from that module's own \"Save Model as RDS\" download)"),
            actionButton(ns("fs_load_btn"), "Load Feature-Selection Features", icon = icon("download"), class = "btn-primary btn-sm")
          ),
          conditionalPanel(sprintf("input['%s'] == 'uploaded'", ns("feature_source")),
            fileInput(ns("panel_upload"), "CpG panel list (CSV/TXT with a cpg/CpG/probe/ID column, or one ID per line)"),
            actionButton(ns("panel_load_btn"), "Load Uploaded Panel", icon = icon("download"), class = "btn-primary btn-sm")
          ),
          conditionalPanel(sprintf("input['%s'] == 'manual'", ns("feature_source")),
            selectizeInput(ns("manual_cpg_select"), "Search / select CpGs", choices = dxm$all_cpgs, multiple = TRUE,
                           options = list(placeholder = "Type to search available CpGs...")),
            actionButton(ns("manual_load_btn"), "Load Selected CpGs", icon = icon("download"), class = "btn-primary btn-sm")
          )
        ),
        if (!is.null(feat$selected) && length(feat$selected) > 0) box(width = 12, status = "primary", solidHeader = TRUE,
          title = sprintf("Selected features (%d CpG%s)", length(feat$selected), if (length(feat$selected) == 1) "" else "s"),
          DT::dataTableOutput(ns("feature_table")))
      )
    })

    output$feature_table <- DT::renderDataTable({
      req(feat$table); DT::datatable(feat$table, rownames = FALSE, options = list(pageLength = 10))
    })

    observeEvent(input$wgcna_load_btn, {
      req(dxm$validated)
      dxm$leakage_safe <- FALSE
      if (identical(dxm$mode, "preloaded")) {
        mt <- dxm_load_wgcna_for_sex(dxm$sex)
        validate(need(!is.null(mt), "Published WGCNA module assignment isn't available in this deployment."))
        cpg_col <- intersect(c("cpg", "CpG", "probe", "ID"), colnames(mt))[1]
        mod_col <- intersect(c("module", "Module", "color"), colnames(mt))[1]
        validate(need(!is.na(cpg_col), "Unexpected WGCNA module-assignment table format (no cpg/CpG/probe/ID column)."))
        sub <- mt
        if (!is.na(mod_col) && !is.null(input$wgcna_module) && !identical(input$wgcna_module, "__all__")) sub <- mt[mt[[mod_col]] == input$wgcna_module, , drop = FALSE]
        cpgs <- intersect(as.character(sub[[cpg_col]]), dxm$all_cpgs)
        src_tbl <- sub[as.character(sub[[cpg_col]]) %in% cpgs, , drop = FALSE]
        names(src_tbl)[names(src_tbl) == cpg_col] <- "cpg"
      } else {
        req(input$wgcna_upload)
        up <- tryCatch(as.data.frame(data.table::fread(input$wgcna_upload$datapath, showProgress = FALSE)), error = function(e) NULL)
        validate(need(!is.null(up), "Could not parse the uploaded WGCNA CSV."))
        cpg_col <- intersect(c("cpg", "CpG", "probe", "ID"), colnames(up))[1]
        validate(need(!is.na(cpg_col), "Uploaded WGCNA CSV must have a cpg/CpG/probe/ID column."))
        cpgs <- intersect(as.character(up[[cpg_col]]), dxm$all_cpgs)
        src_tbl <- up[as.character(up[[cpg_col]]) %in% cpgs, , drop = FALSE]
        names(src_tbl)[names(src_tbl) == cpg_col] <- "cpg"
      }
      validate(need(length(cpgs) > 0, "None of the WGCNA-derived CpGs are present in the active training data."))
      n_top <- input$wgcna_top_n %||% 0
      if (isTRUE(n_top > 0)) cpgs <- utils::head(cpgs, n_top)
      feat$source <- "wgcna"; feat$table <- src_tbl[match(cpgs, src_tbl$cpg), , drop = FALSE]; feat$selected <- cpgs
      showNotification(sprintf("Loaded %d WGCNA-derived CpG(s).", length(cpgs)), type = "message")
    })

    observeEvent(input$fs_load_btn, {
      req(dxm$validated)
      dxm$leakage_safe <- FALSE
      if (identical(dxm$mode, "preloaded")) {
        tbl <- dxm_load_fs_votes_for_sex(dxm$sex)
        validate(need(!is.null(tbl), "Published Feature Selection ensemble-vote table isn't available in this deployment."))
        min_votes <- input$fs_min_votes %||% 2
        sub <- tbl[tbl$n_votes >= min_votes, , drop = FALSE]
        cpgs <- intersect(as.character(sub$cpg), dxm$all_cpgs)
        src_tbl <- sub[as.character(sub$cpg) %in% cpgs, , drop = FALSE]
      } else {
        req(input$fs_upload)
        loaded <- safe_read_rds(input$fs_upload$datapath)
        obj <- if (isTRUE(loaded$ok)) loaded$value else NULL
        validate(need(!is.null(obj) && identical(obj$module, "mod_methyl_featureselection"), loaded$error %||% "Upload a Feature Selection RDS export (from that module's own \"Save Model as RDS\" download)."))
        panel_ids <- as.character(obj$final_panel$cpg_ids %||% character(0))
        cpgs <- intersect(panel_ids, dxm$all_cpgs)
        src_tbl <- data.frame(cpg = cpgs)
        holdout_ids <- as.character(obj$holdout_sample_ids %||% character(0))
        if (length(holdout_ids) > 0) {
          applied <- dxm_apply_holdout_split(holdout_ids)
          if (isTRUE(applied)) {
            showNotification(sprintf("Leakage-safe: using the %d sample(s) this panel reserved as held-out as the internal test set.", length(dxm$test_internal_y)), type = "message")
          } else {
            showNotification("This panel's held-out samples don't overlap enough with the currently loaded cohort - falling back to the Validate & Split partition (not leakage-safe against this panel).", type = "warning")
          }
        }
      }
      validate(need(length(cpgs) > 0, "None of the Feature-Selection-derived CpGs are present in the active training data."))
      n_top <- input$fs_top_n %||% 0
      if (isTRUE(n_top > 0)) cpgs <- utils::head(cpgs, n_top)
      feat$source <- "featureselection"; feat$table <- src_tbl[match(cpgs, src_tbl$cpg), , drop = FALSE]; feat$selected <- cpgs
      showNotification(sprintf("Loaded %d Feature-Selection-derived CpG(s).", length(cpgs)), type = "message")
    })

    observeEvent(input$panel_load_btn, {
      req(dxm$validated, input$panel_upload)
      dxm$leakage_safe <- FALSE
      is_txt <- grepl("\\.txt$", input$panel_upload$name, ignore.case = TRUE)
      if (is_txt) {
        pl <- methyl_parse_probe_list(input$panel_upload$datapath, input$panel_upload$name)
        validate(need(isTRUE(pl$ok), pl$error %||% "Could not parse the uploaded panel file."))
        ids <- pl$ids
      } else {
        up <- tryCatch(as.data.frame(data.table::fread(input$panel_upload$datapath, showProgress = FALSE)), error = function(e) NULL)
        validate(need(!is.null(up), "Could not parse the uploaded panel file."))
        cpg_col <- intersect(c("cpg", "CpG", "probe", "probe_id", "feature", "ID"), colnames(up))[1]
        ids <- if (!is.na(cpg_col)) as.character(up[[cpg_col]]) else as.character(up[[1]])
      }
      cpgs <- intersect(ids, dxm$all_cpgs)
      validate(need(length(cpgs) > 0, "None of the uploaded panel's CpG IDs are present in the active training data."))
      feat$source <- "uploaded"; feat$table <- data.frame(cpg = cpgs); feat$selected <- cpgs
      showNotification(sprintf("Loaded %d CpG(s) from the uploaded panel.", length(cpgs)), type = "message")
    })

    observeEvent(input$manual_load_btn, {
      req(dxm$validated)
      dxm$leakage_safe <- FALSE
      cpgs <- intersect(input$manual_cpg_select, dxm$all_cpgs)
      validate(need(length(cpgs) > 0, "Select at least one CpG present in the active training data."))
      feat$source <- "manual"; feat$table <- data.frame(cpg = cpgs); feat$selected <- cpgs
      showNotification(sprintf("Loaded %d manually selected CpG(s).", length(cpgs)), type = "message")
    })

    output$filters_ui <- renderUI({
      req(dxm$validated)
      tagList(
        box(width = 12, status = "primary", solidHeader = TRUE, title = "Analysis type & feature panel",
          if (length(feat$selected) > 0)
            (if (isTRUE(dxm$leakage_safe))
               p(class = "empty-note", icon("shield-halved"), "Leakage-safe: the internal test set below is the held-out sample set this panel's CpGs were never selected against.")
             else
               p(class = "empty-note", icon("triangle-exclamation"), "Not leakage-safe: this panel's CpGs were not confirmed to be selected without seeing the internal test-set samples - treat any internal-test AUC as exploratory, not a validated estimate."))
          else NULL,
          radioButtons(ns("analysis_type"), "Analysis Type", choices = c("Single CpG" = "single", "Combined CpG Panel" = "combined"), selected = "combined", inline = TRUE),
          conditionalPanel(sprintf("input['%s'] == 'single'", ns("analysis_type")),
            selectInput(ns("single_cpg"), "Select CpG", choices = if (length(feat$selected) > 0) feat$selected else dxm$all_cpgs)),
          conditionalPanel(sprintf("input['%s'] == 'combined'", ns("analysis_type")),
            if (length(feat$selected) > 0) p(sprintf("Combined panel: %d CpG(s) from %s.", length(feat$selected),
                switch(feat$source %||% "", wgcna = "Methylomics WGCNA", featureselection = "Methylomics Feature Selection",
                       uploaded = "an uploaded panel", manual = "manual selection", "no source yet")))
            else p(class = "text-muted", "Load a feature set on the Feature Source tab first."))
        ),
        box(width = 12, status = "primary", solidHeader = TRUE, title = "Class balance",
          DT::dataTableOutput(ns("class_balance_table")),
          radioButtons(ns("imbalance_mode"), "Imbalance handling", inline = TRUE,
            choices = c("None" = "none", "Class weighting (fold-safe up-sampling)" = "weighted", "SMOTE (fold-safe)" = "smote"))
        ),
        box(width = 12, status = "primary", solidHeader = TRUE, title = "Cross-validation & search",
          fluidRow(
            column(3, numericInput(ns("cv_folds"), "CV folds", value = 10, min = 3, max = 20, step = 1)),
            column(3, numericInput(ns("cv_repeats"), "CV repeats", value = 1, min = 1, max = 10, step = 1)),
            column(3, radioButtons(ns("search_method"), "Hyperparameter search", choices = c("Grid Search" = "grid", "Randomized Search" = "random"), selected = "grid")),
            column(3, selectInput(ns("threshold_strategy"), "Classification threshold",
                     choices = c("Default (0.50)" = "default", "Youden's J" = "youden", "Sensitivity-focused (>=0.90 sens)" = "sensitivity", "Specificity-focused (>=0.90 spec)" = "specificity")))
          )
        ),
        box(width = 12, status = "success", solidHeader = TRUE, title = "Run",
          p("Runs a model using the panel and settings configured above (identical to clicking \"Run Model\" on that model's own tab), then jumps you there to see the results."),
          fluidRow(
            column(6, selectInput(ns("run_model_choice"), "Model",
                     choices = stats::setNames(names(DXM_MODEL_SPECS), vapply(DXM_MODEL_SPECS, function(s) s$label, character(1))))),
            column(6, br(), actionButton(ns("run_selected_model_btn"), "Run Model", icon = icon("play"), class = "btn-primary"))
          )
        )
      )
    })

    output$class_balance_table <- DT::renderDataTable({
      req(dxm$validated)
      tbl <- table(dxm$train_y)
      df <- data.frame(Class = c(dxm$ref_level, dxm$comp_level), N = c(tbl[[DXM_NEG]] %||% 0, tbl[[DXM_POS]] %||% 0))
      df$Percent <- round(100 * df$N / sum(df$N), 1)
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE))
    })

    # Automatic leakage-safe headline metric: whenever this session has no confirmed
    # genuine held-out split (dxm$leakage_safe == FALSE - the default/preloaded path,
    # or any live feature-source load that didn't apply a real holdout), compute the
    # nested-CV AUC right here so it's ready by the time any model's results render -
    # no separate manual step needed. Memoized as ordinary reactives: only recomputes
    # when the underlying data/leakage state actually changes.
    dxm_nested_result <- reactive({
      req(dxm$validated)
      if (isTRUE(dxm$leakage_safe)) return(NULL)
      cpgs <- intersect(dxm$all_cpgs, colnames(dxm$full_X))
      if (length(cpgs) < 5) {
        return(list(
          pooled = list(available = FALSE, reason = "Fewer than 5 candidate CpGs are available for automatic leakage-safe validation."),
          per_fold = data.frame(fold = integer(0), n_panel = integer(0), auc = numeric(0)),
          outer_k = NA_integer_, n_folds_completed = 0
        ))
      }
      if (length(cpgs) > DXM_MAX_CANDIDATE_CPGS) {
        v <- apply(dxm$full_X[, cpgs, drop = FALSE], 2, stats::var)
        cpgs <- names(sort(v, decreasing = TRUE))[seq_len(DXM_MAX_CANDIDATE_CPGS)]
      }
      tryCatch(
        withProgress(
          message = "Computing automatic leakage-safe nested-CV AUC (headline metric)...",
          value = 0.6,
          dxm_validate_nested(dxm$full_X[, cpgs, drop = FALSE], dxm$full_y, outer_k = 5, seed = input$dxm_seed %||% 42)
        ),
        error = function(e) NULL
      )
    })
    dxm_headline <- reactive({
      dxm_attach_headline(isTRUE(dxm$leakage_safe), dxm_nested_result())
    })

    model_states <- lapply(DXM_MODEL_SPECS, function(spec) {
      ms <- reactiveValues(fit = NULL, analysis_type = NULL, feature_ids = NULL, ran_at = NULL,
                            train_prob = NULL, train_roc = NULL, train_metrics = NULL, confusion_train = NULL,
                            cv_roc = NULL, threshold = NULL,
                            test_internal_prob = NULL, test_internal_roc = NULL, test_internal_metrics = NULL, confusion_test = NULL,
                            roc_generated = FALSE, calib = NULL, calib_generated = FALSE, lc = NULL, lc_generated = FALSE,
                            last_roc_plot = NULL, last_calib_plot = NULL, last_lc_plot = NULL)
      dxm_register_model_server(spec$id, spec, input, output, session, ns, dxm, feat, ms, runs, results, dxm_headline = dxm_headline)
      ms
    })
    names(model_states) <- names(DXM_MODEL_SPECS)

    observeEvent(input$run_selected_model_btn, {
      req(dxm$validated)
      mid <- input$run_model_choice %||% names(DXM_MODEL_SPECS)[1]
      spec <- DXM_MODEL_SPECS[[mid]]
      validate(need(!is.null(spec), "Unknown model selected."))
      ok <- dxm_do_run_model(mid, spec, input, dxm, feat, model_states[[mid]])
      if (isTRUE(ok)) updateTabsetPanel(session, "main_tabs", selected = spec$label)
    })

    output$compare_ui <- renderUI({
      keys <- names(shiny::reactiveValuesToList(runs))
      if (length(keys) == 0) return(p(class = "text-muted", "No completed model runs yet - run a model and evaluate test data on any model tab first."))
      tagList(
        box(width = 12, status = "primary", solidHeader = TRUE, title = "Model Comparison",
          p(class = "submodule-desc", icon("triangle-exclamation"),
            " \"Test AUC\" is drawn from the same cohort the WGCNA- and Feature-Selection-derived CpG panels were originally selected on, so it can be optimistically biased for those two feature sources - it is internal-test performance only. For a fully independent held-out cohort, evaluate the trained model in the Validation sub-module (External Validation)."),
          DT::dataTableOutput(ns("compare_table")),
          downloadButton(ns("compare_download"), "Download comparison (CSV)", class = "btn-sm")
        ),
        box(width = 12, status = "primary", solidHeader = TRUE, title = "Compare ROC Curves",
          selectizeInput(ns("compare_select"), "Select runs to compare",
            choices = stats::setNames(keys, vapply(keys, function(k) sprintf("%s - %s", runs[[k]]$label, runs[[k]]$analysis_type), character(1))), multiple = TRUE),
          selectInput(ns("compare_curve"), "Curve to compare", choices = c("Test" = "test", "Training" = "train", "Cross-Validated" = "cv")),
          actionButton(ns("compare_roc_btn"), "Generate ROC Comparison", icon = icon("chart-area"), class = "btn-primary"),
          br(), br(),
          if (isTRUE(compare_state$generated)) plotOutput(ns("compare_roc_plot"))
        ),
        box(width = 12, status = "primary", solidHeader = TRUE, title = "Single CpG vs Combined Panel",
          selectInput(ns("svp_model"), "Model", choices = stats::setNames(names(DXM_MODEL_SPECS), vapply(DXM_MODEL_SPECS, function(s) s$label, character(1)))),
          p(class = "submodule-desc", "Shows this model's single-CpG run alongside its combined-panel run, when both have been executed and test-evaluated."),
          DT::dataTableOutput(ns("svp_table"))
        )
      )
    })

    output$compare_table <- DT::renderDataTable({
      all_runs <- shiny::reactiveValuesToList(runs)
      req(length(all_runs) > 0)
      tbl <- do.call(rbind, lapply(all_runs, function(r) data.frame(
        Model = r$label, `Feature set` = sprintf("%s (%d)", r$analysis_type, length(r$feature_ids)),
        `Train AUC` = round(r$train_roc$auc, 3), `CV AUC` = round(r$cv_roc$mean_auc, 3),
        `Test AUC` = round(r$test_internal_roc$auc, 3),
        Threshold = round(r$threshold, 3), `Ran at` = r$ran_at, check.names = FALSE)))
      DT::datatable(tbl, rownames = FALSE, filter = "top", options = list(pageLength = 10))
    })

    observeEvent(input$compare_roc_btn, {
      sel <- input$compare_select
      validate(need(length(sel) > 0, "Select at least one run to compare."))
      all_runs <- shiny::reactiveValuesToList(runs)
      bundles <- stats::setNames(lapply(sel, function(k) switch(input$compare_curve,
                    test = all_runs[[k]]$test_internal_roc, train = all_runs[[k]]$train_roc, cv = all_runs[[k]]$cv_roc$overall)),
                    vapply(sel, function(k) sprintf("%s (%s)", all_runs[[k]]$label, all_runs[[k]]$analysis_type), character(1)))
      compare_state$bundles <- bundles; compare_state$generated <- TRUE
    })
    output$compare_roc_plot <- renderPlot({
      req(compare_state$generated)
      dxm_plot_roc_compare(compare_state$bundles, sprintf("ROC Comparison (%s)",
        switch(input$compare_curve, test = "Test Internal Data", train = "Training", cv = "Cross-Validated")))
    })
    output$compare_download <- downloadHandler(
      filename = function() "methylomics_diagnostic_classifier_comparison.csv",
      content = function(file) {
        all_runs <- shiny::reactiveValuesToList(runs)
        tbl <- do.call(rbind, lapply(all_runs, function(r) data.frame(
          model = r$label, analysis_type = r$analysis_type, n_features = length(r$feature_ids),
          features = paste(r$feature_ids, collapse = ";"), threshold = r$threshold,
          train_auc = r$train_roc$auc, cv_auc = r$cv_roc$mean_auc, cv_auc_sd = r$cv_roc$sd_auc,
          test_auc = r$test_internal_roc$auc, ran_at = r$ran_at)))
        utils::write.csv(tbl, file, row.names = FALSE)
      }
    )

    output$svp_table <- DT::renderDataTable({
      req(input$svp_model)
      all_runs <- shiny::reactiveValuesToList(runs)
      rel <- Filter(function(r) identical(r$model_id, input$svp_model), all_runs)
      req(length(rel) > 0)
      tbl <- do.call(rbind, lapply(rel, function(r) data.frame(
        `Analysis type` = r$analysis_type, `Feature(s)` = paste(r$feature_ids, collapse = ", "),
        `Train AUC` = round(r$train_roc$auc, 3), `CV AUC` = round(r$cv_roc$mean_auc, 3),
        `Test AUC` = round(r$test_internal_roc$auc, 3), check.names = FALSE)))
      DT::datatable(tbl, rownames = FALSE, options = list(dom = "t", paging = FALSE))
    })

    output$testdata_ui <- renderUI({
      req(dxm$validated)
      intr <- dxm_intersect_features(colnames(dxm$train_X), colnames(dxm$test_internal_X))
      out <- list(box(width = 12, status = "primary", solidHeader = TRUE, title = "Test Internal Data",
        p(sprintf("%d training samples; %d test samples.", nrow(dxm$train_X), nrow(dxm$test_internal_X))),
        tags$ul(tags$li(sprintf("Training features: %d", length(intr$train))), tags$li(sprintf("Test features: %d", length(intr$test))),
                tags$li(sprintf("Shared features: %d", length(intr$shared))), tags$li(sprintf("Unmatched (dropped) features: %d", length(intr$unmatched))))))
      tagList(out)
    })

    output$export_ui <- renderUI({
      keys <- names(shiny::reactiveValuesToList(runs))
      if (length(keys) == 0) return(p(class = "text-muted", "No completed model runs yet - run a model and evaluate test data first."))
      box(width = 12, status = "primary", solidHeader = TRUE, title = "Export",
        downloadButton(ns("export_metrics_csv"), "Download all metrics (CSV)"),
        downloadButton(ns("export_panel_csv"), "Download selected feature panel (CSV)"))
    })
    output$export_metrics_csv <- downloadHandler(
      filename = function() "methylomics_diagnostic_classifier_metrics.csv",
      content = function(file) {
        all_runs <- shiny::reactiveValuesToList(runs)
        tbl <- do.call(rbind, lapply(all_runs, function(r) data.frame(
          model = r$label, analysis_type = r$analysis_type, n_features = length(r$feature_ids),
          features = paste(r$feature_ids, collapse = ";"), threshold = r$threshold,
          train_auc = r$train_roc$auc, cv_auc = r$cv_roc$mean_auc, cv_auc_sd = r$cv_roc$sd_auc,
          test_auc = r$test_internal_roc$auc, ran_at = r$ran_at)))
        utils::write.csv(tbl, file, row.names = FALSE)
      }
    )
    output$export_panel_csv <- downloadHandler(
      filename = function() "methylomics_diagnostic_classifier_feature_panel.csv",
      content = function(file) utils::write.csv(feat$table %||% data.frame(cpg = feat$selected), file, row.names = FALSE)
    )

    invisible(NULL)
  })
}
