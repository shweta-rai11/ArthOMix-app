## R/mod_diagnostic.R
## Submodule: Diagnostic Model (Section 2.9)
## "Your analysis" fits FOUR live classifiers - plain logistic regression,
## elastic-net logistic regression, random forest and SVM - on a user-chosen
## gene panel and two-group contrast, sex-stratified exactly like Feature
## Selection (female and male are ALWAYS fit completely separately, never
## combined with sex as a covariate). Two tabs, matching the two tiers this
## is evaluated at:
##   "Model Training" - each sex's own available samples are split ONCE
##       (stratified, seed 1234) into Train / Test at a user-chosen ratio
##       (default 70:30 - see the "Test-set size" slider in the left panel).
##       This tab fits and tunes on Train ONLY, and reports Train's own
##       full-fit ROC/AUC plus a k-fold cross-validated AUC by fold (not just
##       in-sample fit) - Test is never touched here.
##   "Model Testing (Internal)" - the SAME Train-fit model, applied ONCE to
##       the Test split held out above ("internal" because it's held out
##       from the same loaded cohort, not a separate dataset) - genuine
##       held-out performance, not in-sample. This tab deliberately stops at
##       that: it does NOT score this project's own separately bundled
##       external validation holdout or blood cohort (two different,
##       cross-platform datasets this project's own scripts also validate
##       against) - out of scope for a live, point-and-click training tool
##       that only ever has the one dataset currently loaded on the Dataset
##       tab to work with.
## This tab also does NOT reproduce this project's own nested (nested =
## feature-selection-redone-per-fold) CV - deliberately out of scope for the
## same reason. See the per-method parameter boxes below for exactly which
## simplification that is.
##
## TWO of these four are this project's own methodology, not just one:
##   - Elastic net (15_model_training_elasticnet.R / 16_model_training_
##     final_panel.R): alpha tuned by minimising CV deviance over a grid,
##     lambda by glmnet's own inner CV. This is the FINAL LOCKED model
##     16_model_training_final_panel.R reports (full-train fit, scored once
##     on the internal + external holdouts).
##   - Logistic regression, plain and unpenalized (14_model_training_
##     nested_cv.R / 16b_model_training_final_panel_noMHC.R / 16d_nested_cv_
##     reconciliation.R): glm(y ~ ., family = binomial) fit directly on the
##     gene panel, no regularisation or tuning. This is the classifier every
##     one of this project's OWN nested-CV and downstream validation scripts
##     (14, 16b, 16d, and the blood/synovium/cross-tissue testing scripts,
##     17/17b/21/22/23) actually scores - 16d's own "authoritative"
##     reconciliation table marks its noMHC + logistic-regression row
##     `recommended`, standing on equal footing with elastic net rather than
##     as a discarded baseline. Report both; neither is "the" winner here.
## Random forest and SVM are NOT part of this project's own diagnostic-model
## methodology (this project only ever used them as feature *selectors*, in
## Feature Selection) - they're offered here as additional classifier choices
## on the same gene panel, so treat their defaults as reasonable off-the-shelf
## settings rather than a methodology citation.
##
## SCALING: the Test split is the SAME loaded dataset/platform as Train,
## just samples Train never saw, so it's standardised with TRAIN's own
## per-gene mean/SD (leakage-free supervised transfer - standard practice
## for an in-dataset split). The k-fold CV-within-Train step instead uses
## the FOLD-TRAIN mean/SD only (leakage-free, matching this project's own
## nested-CV scripts) - both are correct for what each is estimating.
##
## HYPERPARAMETERS are tuned ONCE per sex on Train (alpha / mtry / cost) and
## then reused, unchanged, inside every k-fold CV refit - only the model's
## coefficients are refit per fold. This project's own nested-CV scripts
## instead re-tune inside every fold; skipping that here is a deliberate
## simplification to keep a live "Run Female" / "Run Male" click fast on
## this UI. The result is therefore an upper bound, same as this project's
## own flat-vs-nested CV documentation already says about its own flat
## comparison - not a new caveat invented for this tab.

mod_diagnostic_config <- list(
  id = "diagnostic", group = "Biomarker modeling",
  title = "Diagnostic Model",
  description = "Sex-stratified logistic regression, elastic-net, random forest and SVM models, trained on a Train:Test split and evaluated on the held-out Test set.",
  icon = "stethoscope"
)

## ---------------------------------------------------------------------------
## Shared fitting helpers
## ---------------------------------------------------------------------------

## Direct port of this project's own zrows() (16_model_training_final_panel.R
## / 15_model_training_elasticnet.R): per-GENE (row) z-score, independent of
## any other dataset's mean/SD.
diag_zrows <- function(M) {
  t(apply(M, 1, function(v) {
    s <- stats::sd(v, na.rm = TRUE)
    if (is.na(s) || s == 0) rep(0, length(v)) else (v - mean(v, na.rm = TRUE)) / s
  }))
}

DIAG_SVM_COST_GRID <- c(0.01, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16)

DIAG_DEFAULT_PARAMS <- list(
  test_frac = 0.3,
  ## Class weighting - shared across all four models (see diag_class_weight_
  ## levels()/diag_obs_weights() below), one control for the whole per-sex
  ## run rather than four duplicated ones, same "one shared setting" pattern
  ## as test_frac above. "equal" = this project's own methodology (no
  ## reweighting); "balanced" = auto inverse-frequency (larger class
  ## down-weighted so the effective sample size is unchanged); "manual" =
  ## a fixed comparison:reference weight ratio.
  class_weight_mode = "equal", class_weight_ratio = 1,
  lr_cv_folds = 5,
  enet_cv_folds = 5, enet_alpha_grid = c(0.1, 0.3, 0.5, 0.7, 0.9, 1.0), enet_lambda_choice = "lambda.min",
  enet_nlambda = 100, enet_type_measure = "deviance",
  rf_cv_folds = 5, rf_ntree = 1000, rf_mtry_mode = "auto", rf_mtry_manual = NULL,
  rf_nodesize = 1, rf_maxnodes = NULL,
  svm_cv_folds = 5, svm_kernel = "linear", svm_cost_mode = "auto", svm_cost_manual = 1, svm_cost_grid = DIAG_SVM_COST_GRID,
  svm_gamma_mode = "auto", svm_gamma_manual = 1, svm_degree = 3, svm_tolerance = 0.001
)

## Named-by-level (reference, comparison) weights for imbalanced groups -
## shared by randomForest's classwt, e1071's class.weights, and (converted to
## a per-observation vector by diag_obs_weights()) glm/glmnet's weights=.
## "balanced" divides the larger class's count into itself so the smaller
## class gets weight > 1 (standard inverse-frequency reweighting); "manual"
## fixes the comparison-class weight at class_weight_ratio (reference always
## 1); "equal" (default) returns all-1s, i.e. this project's own unweighted
## methodology.
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

## Stratified train/test split of one sex's full available pool - a single
## fixed split (seed 1234), not repeated/resampled, so "Test split" below
## means the same held-out samples every run unless test_frac itself changes.
diag_split_train_test <- function(y_full, test_frac, seed = 1234) {
  set.seed(seed)
  train_idx <- as.integer(caret::createDataPartition(y_full, p = 1 - test_frac, list = FALSE))
  list(train = train_idx, test = setdiff(seq_along(y_full), train_idx))
}

## k-fold CV AUC for one already-tuned classifier: re-standardises with the
## FOLD-TRAIN mean/SD only (leakage-free) and refits `refit_fn`'s model type
## on each training fold, predicting the untouched fold with `predict_fn`.
## `Xraw` is samples x genes, NOT yet z-scored (scaling happens per fold,
## inside this loop).
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

## AUC + CI, matching this project's own auc_ci() exactly: DeLong when
## n >= 20, stratified bootstrap (seed 1234, 2000 reps) when n < 20 - the
## same small-n handling 16_model_training_final_panel.R uses for its own
## internal/external scoring (Carpenter & Bithell 2000).
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

## Publication-style ROC curve, used for every ROC plot in this module
## (Training, Test-split, and External Validation could all reuse it) -
## deliberately NOT pROC::ggroc()'s own default rendering:
##  - x-axis is 1 - specificity (false positive rate) increasing left to
##    right, y-axis is sensitivity (true positive rate) increasing bottom
##    to top - the convention most biomedical journals use (and the one
##    Chen et al. 2021/2022's own Fig. 4a uses), rather than ggroc()'s
##    default of specificity decreasing 1 -> 0 on the x-axis. Purely a
##    re-orientation of the same curve/AUC - mathematically identical,
##    just the more immediately readable convention for a manuscript figure.
##  - AUC (95% CI) reported both as a subtitle and as in-plot text, the two
##    places a reviewer looks for it first.
##  - 1:1 aspect ratio (coord_equal) so the diagonal chance line reads as a
##    true 45 degrees, a solid panel border, and larger, bold axis titles -
##    all standard journal-figure conventions a dashboard-style ggplot
##    theme (theme_arthomix, used everywhere else in this app) doesn't aim
##    for on its own.
diag_roc_plot_pub <- function(roc_obj, color, title = NULL, ci = NULL) {
  co <- pROC::coords(roc_obj, "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
  df <- data.frame(fpr = 1 - co$specificity, tpr = co$sensitivity)
  df <- df[order(df$fpr, df$tpr), ]
  auc_val <- as.numeric(pROC::auc(roc_obj))
  ci_val <- if (!is.null(ci) && all(is.finite(ci))) ci else diag_auc_ci(roc_obj)[c("lo", "auc", "hi")]
  auc_label <- if (all(is.finite(ci_val))) {
    sprintf("AUC = %.3f\n(95%% CI %.3f-%.3f)", auc_val, ci_val[1], ci_val[3])
  } else {
    sprintf("AUC = %.3f", auc_val)
  }
  ggplot(df, aes(x = fpr, y = tpr)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#9CA3AF", linewidth = 0.6) +
    geom_line(color = color, linewidth = 1.1) +
    annotate("text", x = 0.97, y = 0.05, label = auc_label, size = 3.8, hjust = 1, vjust = 0, color = "#1F2937", lineheight = 0.95) +
    scale_x_continuous(name = "1 − Specificity (False Positive Rate)", limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0.01, 0.01)) +
    scale_y_continuous(name = "Sensitivity (True Positive Rate)", limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0.01, 0.01)) +
    coord_equal() +
    labs(title = title, subtitle = sub("\n", ", ", auc_label)) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#EEF0F3", linewidth = 0.4),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 11, color = "#374151"),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 11, color = "black"),
      panel.border = element_rect(color = "black", linewidth = 0.6, fill = NA)
    )
}

## Sensitivity/specificity/accuracy at a FIXED, already-locked threshold
## (the training Youden cutoff) - never re-optimised on the Test split, since
## re-optimising a cutoff on the same data used to report performance would
## leak exactly the way this project's own methodology warns against.
## Matches this project's own fmt() (16_model_training_final_panel.R): a CI
## collapsed against 1.000 at small n is a perfect-separation artefact of few
## discordant pairs, not evidence of a strong classifier.
diag_separation_note <- function(ev) {
  if (isTRUE(ev$available) && !is.na(ev$ci_lo) && ev$ci_lo >= 0.999) {
    sprintf(" [separation, n=%d]", ev$n)
  } else ""
}

## Human-readable selected-hyperparameter value for the Training tab's KPI
## tile - one model type's chosen setting, whether it came from CV tuning or
## a manual override.
diag_hyperparam_value <- function(rr) {
  switch(rr$model_type,
    lr = "none (unpenalized)",
    enet = sprintf("α = %.2f", rr$alpha),
    rf = sprintf("mtry = %d", rr$mtry),
    svm = sprintf("cost = %s", format(rr$cost, trim = TRUE))
  )
}

## One ROC curve per gene, using that gene's own raw expression as a
## single-variable predictor - model-agnostic (a gene's own univariate
## diagnostic value doesn't depend on which multivariate classifier is later
## fit on top of it, and AUC is invariant to monotonic transforms like
## z-scoring, so this runs directly on raw expression - no re-scaling
## needed). Same shape whether given the Train or the Test split.
## p-value alongside AUC: a two-sided Wilcoxon rank-sum test of this gene's
## raw expression between the two groups - the standard companion statistic
## to a univariate AUC in this kind of biomarker screen (e.g. Chen et al.
## 2021/2022's own hub-gene rule: "AUC > 85% and a P value < 0.05"), and
## rank-based so it needs no normality assumption on raw expression, unlike
## a t-test. Kept independent of the AUC computation (its own tryCatch) so
## one failing gene's test doesn't blank out that gene's otherwise-valid AUC.
diag_gene_roc <- function(expr_sub, y) {
  genes <- rownames(expr_sub)
  rocs <- vector("list", length(genes)); names(rocs) <- genes
  aucs <- setNames(numeric(length(genes)), genes)
  pvals <- setNames(numeric(length(genes)), genes)
  for (g in genes) {
    ## direction = "auto" (not a hardcoded "<") - pROC picks whichever of
    ## "<"/">" gives AUC >= 0.5. A fixed "<" would report a gene that's
    ## DOWN-regulated in the comparison group as e.g. AUC = 0.05, making a
    ## strong down-regulated marker look like noise under an "AUC >=
    ## threshold" hub-gene filter, when 1 - 0.05 = 0.95 is its real
    ## discriminative power just in the other direction.
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
## FULL available sample pool, which this function itself splits into a
## Train / Test partition (see diag_split_train_test() - a single fixed
## stratified split, not repeated). Every fit, every hyperparameter, and the
## k-fold CV are all computed on Train only; Test is genuinely held out and
## scored once at the end. `expr_full` is genes x all available samples for
## this sex (raw values, not yet z-scored).
diag_fit_sex <- function(expr_full, y_full, params = list()) {
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

  ## Class weighting (see diag_class_weight_levels() above) - computed once
  ## on the full training split; diag_cv_auc()'s refit_fn closures below each
  ## recompute their own fold-train weights from fold_y directly (never
  ## these outer ones), same leakage-free convention as the fold-local
  ## z-scoring already documented on diag_cv_auc() itself.
  cw_levels <- diag_class_weight_levels(y, params$class_weight_mode, params$class_weight_ratio)
  obs_w <- diag_obs_weights(y, params$class_weight_mode, params$class_weight_ratio)

  expr_train_sub <- expr_full[, split$train, drop = FALSE]
  expr_test_sub <- expr_full[, split$test, drop = FALSE]

  Ztr <- diag_zrows(expr_train_sub); rownames(Ztr) <- safe
  Xtr_full <- t(Ztr)                                    # z-scored, sample x gene(safe) - used for tuning + full fit
  Xraw <- t(expr_train_sub); colnames(Xraw) <- safe      # raw, sample x gene(safe) - re-scaled per fold inside diag_cv_auc

  ## Test split - the SAME loaded dataset/platform as Train, just samples
  ## Train never saw, so it's correct and standard to standardise it with
  ## TRAIN's own per-gene mean/SD (leakage-free: those stats never touch the
  ## test rows).
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

  ## -------------------------------------------------------------------
  ## (1) Elastic net - this project's own methodology: alpha tuned by
  ## minimising CV deviance over a grid, lambda by glmnet's inner CV.
  ## enet_type_measure lets the CV metric itself be swapped (deviance/
  ## AUC/misclassification error); for "auc", bigger is better, so the
  ## alpha-selection direction below is metric-dependent rather than
  ## always "smallest wins".
  ## -------------------------------------------------------------------
  nf_a <- max(2, min(params$enet_cv_folds, min(table(y))))
  type_measure <- params$enet_type_measure %||% "deviance"
  bigger_is_better <- identical(type_measure, "auc")
  best <- NULL; bcv <- if (bigger_is_better) -Inf else Inf
  ## Captured for the "explore hyperparameter tuning" plot below - every
  ## alpha actually tried and its own best (inner-CV-tuned lambda) score
  ## under the chosen CV metric, not just the winner.
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

  ## -------------------------------------------------------------------
  ## (2) Random forest classifier. Not this project's own methodology (RF
  ## was only ever a feature *selector* here) - mtry tuned by CV, ntree
  ## fixed, both user-overridable. nodesize/maxnodes (tree complexity) and
  ## classwt (class weighting, see diag_class_weight_levels() above) are
  ## additional advanced, user-overridable controls.
  ## -------------------------------------------------------------------
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
    ## Captured for the "explore hyperparameter tuning" plot below - every
    ## mtry actually tried and its own CV ROC (twoClassSummary), not just
    ## the winner.
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

  ## -------------------------------------------------------------------
  ## (3) SVM classifier. Not this project's own methodology either (SVM was
  ## only ever SVM-RFE, a feature selector, here) - cost tuned by CV over a
  ## grid, kernel user-selectable (linear, matching this project's own
  ## SVM-RFE convention, by default; radial and polynomial also available).
  ## scale = FALSE throughout: the data is already gene-wise z-scored above,
  ## so e1071's own internal re-scaling is deliberately skipped rather than
  ## double-applied. gamma (radial/polynomial kernel coefficient, either
  ## e1071's own 1/ncol(x) default or a manual value) and degree (polynomial
  ## only) are ignored by e1071::svm() for a linear kernel, so it's always
  ## safe to pass them regardless of which kernel is selected. tolerance is
  ## the SMO termination tolerance; class.weights is the shared class-
  ## weighting control (see diag_class_weight_levels() above) - epsilon is
  ## deliberately NOT exposed here, since it only affects epsilon-regression
  ## and has no effect on the C-classification this module always fits.
  ## -------------------------------------------------------------------
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
    ## Captured for the "explore hyperparameter tuning" plot below - every
    ## cost actually tried and its own CV classification error, not just
    ## the winner.
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

  ## -------------------------------------------------------------------
  ## (4) Logistic regression - plain, unpenalized multivariable glm on
  ## every gene in the panel. This IS this project's own methodology, just
  ## as much as elastic net is (see the module header): 14_model_training_
  ## nested_cv.R / 16b_model_training_final_panel_noMHC.R / 16d_nested_cv_
  ## reconciliation.R all fit exactly this - glm(y ~ ., family = binomial)
  ## on the z-scored panel - as the classifier scored for nested-CV AUC, and
  ## it's the classifier every downstream validation script (17/17b/21/22/
  ## 23) scores against the internal/external/cross-tissue cohorts. No
  ## hyperparameters to tune - unlike elastic net's regularisation path,
  ## there's no grid here, so tuning_search stays NULL throughout.
  ## -------------------------------------------------------------------
  ## class_weight_mode = "equal" (the default) means obs_w is all 1s, which
  ## makes weights = obs_w numerically identical to the unweighted glm() fit
  ## this project's own scripts use - so this project's own numeric
  ## replication is unaffected unless class weighting is actively turned on.
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

## Model Training tab: one box per model x sex, full width (see mod_diagnostic_ui
## below - "utilise the entire page, results are very narrow" - so the KPI
## tiles and all three plots (ROC / CV-by-fold / hyperparameter tuning) sit
## across the full width of ONE box rather than each panel sharing a
## half-width column with the other sex).
mod_diagnostic_training_panel <- function(ns, prefix, title, roc_height = 250) {
  box(
    width = NULL, title = title, status = "primary", solidHeader = FALSE,
    withSpinner(uiOutput(ns(paste0(prefix, "_train_stats"))), color = "#2563EB", type = 6),
    fluidRow(
      column(4, h5("ROC (training fit)"), withSpinner(plotOutput(ns(paste0(prefix, "_train_roc_plot")), height = roc_height), color = "#2563EB", type = 6)),
      column(4, h5("Cross-validated AUC by fold"), withSpinner(plotOutput(ns(paste0(prefix, "_train_cv_plot")), height = roc_height), color = "#2563EB", type = 6)),
      column(4, h5("Hyperparameter tuning - explore the grid"), withSpinner(plotOutput(ns(paste0(prefix, "_train_tuning_plot")), height = roc_height), color = "#2563EB", type = 6))
    ),
    div(class = "table-toolbar",
        downloadButton(ns(paste0(prefix, "_train_download")), "Performance (CSV)", class = "btn-sm"),
        downloadButton(ns(paste0(prefix, "_model_download")), "Model (.rds)", class = "btn-sm")),
    DT::dataTableOutput(ns(paste0(prefix, "_train_table")))
  )
}

## Model Testing (Internal) tab: same full-width treatment. "Internal" means
## the Test split held out from this sex's OWN loaded samples (per the
## Train:Test ratio control) - this tab only ever scores that split, nothing
## from outside the currently loaded dataset.
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

## One sex's own single-gene ROC curve box (Female and Male get entirely
## separate boxes now, each only ever placed inside that sex's own tab - see
## mod_diagnostic_training_sex_panel()/mod_diagnostic_testing_sex_panel()
## below - so there is no shared box where the other sex's column would show
## up while you're looking at this one).
## Hub-gene rule thresholds: adjustable, not hardcoded to any one paper's
## values - defaults (AUC >= 0.85, P < 0.05) match Chen et al. 2021/2022's
## own rule for this exact WGCNA-module gene panel, but this box is reused
## for every dataset/panel run through this app, so both fields are plain
## live-editable numericInputs rather than baked-in constants. "hub" in the
## table/download below means "passes both thresholds currently set here",
## recomputed reactively - change either number and the flagged set updates
## without re-running the models themselves (AUC/p per gene don't depend on
## the threshold, only which rows count as "hub").
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

## One sex's ENTIRE "Model Training" drill-down: its own Run button, then -
## hidden until that button (or its Model Testing counterpart) is clicked -
## the model pills, per-gene ROC and comparison table, all Female-only or
## Male-only (never both, since Female and Male are now separate OUTER tabs
## with their own independent copy of every box below, not columns sharing
## one box). conditionalPanel (a client-side CSS toggle), not a server-side
## req()-gated renderUI, is what hides this: the boxes still exist in the DOM
## - and therefore keep working with shinycssloaders' spinners - just
## visually hidden until the condition is true, the same way Bootstrap's own
## tab-pane switching already hides the inactive model pill without breaking
## its spinners.
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
      mod_diagnostic_generoc_box_sex(ns, sex_label, "train", sprintf("Per-gene ROC/AUC - Training (%s)", sex_title)),
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
      mod_diagnostic_generoc_box_sex(ns, sex_label, "test", sprintf("Per-gene ROC/AUC - Test (%s)", sex_title)),
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
        arthochat_shortcut_ui(
          "Questions about these models? Ask ArthOChat.",
          compact = TRUE
        ),
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
          ## Same results$wgcna$module_genes source Feature Selection's own
          ## WGCNA-module option reads from (set by mod_wgcna.R's Step 3) -
          ## same module list used for Female, Male and Pooled alike, same
          ## as the "own" pasted-list option above.
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
            ## One params box, switched server-side to match whichever model
            ## pill was clicked most recently in EITHER sex's tab (params are
            ## shared, not sex-specific). A single output driven by an
            ## "active pill" reactiveVal, rather than a conditionalPanel OR-ing
            ## both sexes' raw pill inputs - the OR approach left the
            ## Logistic Regression box stuck visible forever, since the other
            ## sex's pill (never clicked) always sits on its "Logistic
            ## Regression" default and kept satisfying that box's condition.
            uiOutput(ns("model_params_ui")),
            tabsetPanel(
              id = ns("train_sex_tabs"), type = "tabs",
              tabPanel("Female", br(), mod_diagnostic_training_sex_panel(ns, "female")),
              tabPanel("Male", br(), mod_diagnostic_training_sex_panel(ns, "male")),
              ## Pools every sample regardless of sex (or when sex isn't
              ## available) into one model - for replicating a published
              ## method that never stratified by sex, same convention as
              ## Feature Selection's own "Run All (pooled)".
              tabPanel("Pooled (all)", br(), mod_diagnostic_training_sex_panel(ns, "pooled"))
            )
          ),
          tabPanel(
            "Model Testing (Internal)", br(),
            p(class = "submodule-desc", "Each Train-fit model scored once on its held-out Test split. Pick Female, Male or Pooled, then Run."),
            tabsetPanel(
              id = ns("test_sex_tabs"), type = "tabs",
              tabPanel("Female", br(), mod_diagnostic_testing_sex_panel(ns, "female")),
              tabPanel("Male", br(), mod_diagnostic_testing_sex_panel(ns, "male")),
              tabPanel("Pooled (all)", br(), mod_diagnostic_testing_sex_panel(ns, "pooled"))
            )
          ),
          ## A genuinely SEPARATE cohort, not a train/test split of the same
          ## loaded dataset (that's what "Model Testing (Internal)" above
          ## already covers) - the same role Chen et al. 2021/2022's own
          ## GSE77298 plays: upload a different dataset entirely and check
          ## whether the gene panel's univariate AUC/P still holds up on it.
          ## Per-gene, not a refit multivariate model - this is exactly what
          ## that paper's own Fig. 3d/4a validation step is (expression
          ## boxplots + per-gene ROC/AUC on the external cohort), not a new
          ## kind of analysis invented for this app.
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

## External-validation tab: its own file upload (a genuinely different
## dataset than whatever's loaded on the app-wide Dataset tab), its own
## group-column mapping and reference/comparison pick, and a "which panel"
## selector reusing the SAME panel_source config (project/own/wgcna) as
## Model Training/Testing - so whatever gene panel you've already set up on
## the left sidebar is exactly what gets validated here, just against new
## samples instead of a held-out split of the training samples.
mod_diagnostic_external_panel <- function(ns) {
  tagList(
    box(
      width = NULL, title = "External validation dataset", status = "primary", solidHeader = FALSE,
      p(class = "submodule-desc", "Upload a genuinely separate cohort - not a train/test split of the currently loaded dataset, a different one entirely (e.g. the paper's own GSE77298-style external validation set) - to check whether the gene panel below holds up outside the data it was chosen on."),
      fluidRow(
        column(6, fileInput(ns("ext_expr_file"), "Expression matrix", accept = c(".csv", ".rds", ".Rds"))),
        column(6, fileInput(ns("ext_meta_file"), "Sample metadata", accept = c(".csv", ".rds", ".Rds")))
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

    ## Which value in the app's currently loaded sex column means "female" /
    ## "male" - same detection FS uses, since a user upload might not use
    ## this project's own single-letter "F"/"M" convention.
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
            "Applied to all four models (Logistic Regression, Elastic Net, Random Forest, SVM) for both sexes.")
      )
    })

    ## -----------------------------------------------------------------
    ## Gene panel sources.
    ## -----------------------------------------------------------------

    project_panel_genes <- function(sex_label) {
      live <- results$featureselection[[sex_label]]$consensus_genes
      if (!is.null(live) && length(live) >= 2) {
        return(list(genes = live, is_live = TRUE,
                    note = sprintf("%d genes from this session's live Feature Selection %s consensus panel (LASSO ∩ Random Forest ∩ SVM-RFE).", length(live), sex_label)))
      }
      bundled <- read_table_safe(sprintf("FS_input_%s.csv", sex_label))
      if (!is.null(bundled) && nrow(bundled) >= 2 && "gene" %in% colnames(bundled)) {
        return(list(genes = unique(as.character(bundled$gene)), is_live = FALSE,
                    note = sprintf("%d genes from this project's bundled %s panel (FS_input_%s.csv).", nrow(bundled), sex_label, sex_label)))
      }
      list(genes = character(0), is_live = FALSE, note = sprintf("No %s candidate genes available.", sex_label))
    }

    own_panel_genes <- function(sex_label) {
      genes <- unique(trimws(unlist(strsplit(input$gene_list %||% "", "[,\n\t ]+"))))
      genes <- genes[nzchar(genes)]
      list(genes = genes, is_live = FALSE, note = sprintf("%d pasted genes.", length(genes)))
    }

    ## Same module for every sex_label - a WGCNA module isn't sex-specific
    ## by construction, matching own_panel_genes()'s "same list for both
    ## sexes" behavior above.
    wgcna_panel_genes <- function(sex_label) {
      req(input$wgcna_module_pick)
      mg <- results$wgcna$module_genes
      validate(need(!is.null(mg) && input$wgcna_module_pick %in% names(mg),
                    "Run WGCNA Step 3 (Modules) first, then pick a module above."))
      genes <- unique(as.character(mg[[input$wgcna_module_pick]]))
      list(genes = genes, is_live = TRUE,
           note = sprintf("%d genes from WGCNA module \"%s\" (this session).", length(genes), input$wgcna_module_pick))
    }

    ## Module-color dropdown, same results$wgcna$module_genes source and
    ## "grey excluded" convention as mod_featureselection.R's own
    ## wgcna_module_pick_ui_builder() - duplicated rather than shared across
    ## module files, per this codebase's existing per-module
    ## self-containment convention (see mod_wgcna.R's own note on
    ## wgcna_step_title()).
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
            "Uses whichever gene set is currently loaded on the Dataset tab for expression values - make sure that's the same dataset WGCNA was run on.")
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
      } else {
        div(class = "empty-note", icon("circle-info"),
            "No live Feature Selection panel yet - using this project's bundled female/male panel (no bundled pooled panel exists, so Pooled here needs a live Feature Selection run first).")
      }
    })

    ## -----------------------------------------------------------------
    ## External Validation - a genuinely separate uploaded dataset, not a
    ## split of whatever's on the app-wide Dataset tab. Reuses
    ## diag_gene_roc() (the same per-gene AUC/Wilcoxon-P this module's own
    ## hub-gene boxes use) rather than refitting the multivariate models on
    ## new data - this is deliberately the same scope as the paper's own
    ## external-validation step (per-gene ROC/AUC + expression boxplots),
    ## not a new kind of analysis.
    ## -----------------------------------------------------------------

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
        readRDS(input$ext_expr_file$datapath)
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

    output$ext_status_ui <- renderUI({
      r <- tryCatch(ext_result(), error = function(e) NULL)
      if (is.null(r)) {
        return(div(class = "empty-note", icon("circle-info"), "Not run yet. Upload both files, map the columns, pick reference/comparison groups, then click \"Run external validation\"."))
      }
      div(class = "empty-note", icon("check"),
          sprintf("%d panel genes present in the external dataset, %d samples (%d %s vs %d %s). %s",
                  length(r$genes), length(r$y), r$n_comp, r$comp_group, r$n_ref, r$ref_group, r$panel_note))
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

    ## Expression boxplot, one panel per gene, RA vs HC (matching the
    ## paper's own Fig. 3d style) - capped to the top 24 genes by AUC so a
    ## whole-module-sized panel doesn't render an unreadable wall of facets;
    ## the full per-gene table/download above has every gene regardless.
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

    ## -----------------------------------------------------------------
    ## Advanced parameters - always visible and live-editable, one
    ## independent box per model type (mirrors Feature Selection).
    ## -----------------------------------------------------------------

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

    ## -----------------------------------------------------------------
    ## Run: fit elastic net + random forest + SVM for one sex on one click.
    ## -----------------------------------------------------------------

    ## sex_value = NULL means "pooled" - every sample regardless of sex,
    ## skipping the sex filter (and its column requirement) entirely, the
    ## same convention Feature Selection's fs_build_sex() uses. For
    ## replicating a published method (like Chen et al. 2021/2022's own
    ## LASSO/SVM-RFE/AUC pipeline) that never stratified by sex at all.
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
      genes <- intersect(cand$genes, rownames(dataset$expr))
      validate(need(length(genes) >= 2, sprintf("Fewer than 2 %s genes from the chosen panel are present in the currently loaded expression matrix.", sex_label)))

      expr_sub <- dataset$expr[genes, common, drop = FALSE]

      fit <- diag_fit_sex(expr_sub, y, params = diag_advanced_params())
      fit$candidate_note <- cand$note
      fit$n_ref <- sum(y == input$ref_group); fit$n_comp <- sum(y == input$comp_group)
      fit$ref_group <- input$ref_group; fit$comp_group <- input$comp_group
      fit
    }

    ## One Run button per sex on the Training tab, one more on the Testing
    ## tab (mod_diagnostic_training_sex_panel() / mod_diagnostic_testing_sex_
    ## panel() in the UI) so a run can be started from wherever you are -
    ## both feed the same shared per-sex trigger below, and the SAME two
    ## button IDs are what each tab's conditionalPanel condition checks to
    ## decide whether to reveal itself.
    female_run_trigger <- reactiveVal(0)
    male_run_trigger <- reactiveVal(0)
    lapply(c("run_female_btn", "run_female_btn_test"), function(bid) {
      observeEvent(input[[bid]], { female_run_trigger(isolate(female_run_trigger()) + 1) }, ignoreInit = TRUE)
    })
    lapply(c("run_male_btn", "run_male_btn_test"), function(bid) {
      observeEvent(input[[bid]], { male_run_trigger(isolate(male_run_trigger()) + 1) }, ignoreInit = TRUE)
    })
    pooled_run_trigger <- reactiveVal(0)
    lapply(c("run_pooled_btn", "run_pooled_btn_test"), function(bid) {
      observeEvent(input[[bid]], { pooled_run_trigger(isolate(pooled_run_trigger()) + 1) }, ignoreInit = TRUE)
    })

    diag_result_female <- eventReactive(female_run_trigger(), {
      diag_build_sex("female", sex_levels()$female)
    }, ignoreInit = TRUE)
    diag_result_male <- eventReactive(male_run_trigger(), {
      diag_build_sex("male", sex_levels()$male)
    }, ignoreInit = TRUE)
    ## sex_value = NULL: every sample regardless of sex, works even without
    ## a usable sex column - see diag_build_sex()'s own comment above.
    diag_result_pooled <- eventReactive(pooled_run_trigger(), {
      diag_build_sex("pooled", NULL)
    }, ignoreInit = TRUE)

    diag_has_run <- reactiveVal(FALSE)
    observeEvent(female_run_trigger(), diag_has_run(TRUE), ignoreInit = TRUE)
    observeEvent(male_run_trigger(), diag_has_run(TRUE), ignoreInit = TRUE)
    observeEvent(pooled_run_trigger(), diag_has_run(TRUE), ignoreInit = TRUE)

    ## Which model's params box to show: whichever pill was clicked most
    ## recently, in any sex's own Model Training tab.
    active_model_pill <- reactiveVal("Logistic Regression")
    observeEvent(input$female_train_model_pills, active_model_pill(input$female_train_model_pills), ignoreInit = TRUE)
    observeEvent(input$male_train_model_pills, active_model_pill(input$male_train_model_pills), ignoreInit = TRUE)
    observeEvent(input$pooled_train_model_pills, active_model_pill(input$pooled_train_model_pills), ignoreInit = TRUE)

    lr_params_box <- function() {
      mod_diagnostic_params_box(
        ns, "lr", "Logistic Regression",
        "This IS this project's own methodology, same as elastic net (see 14_model_training_nested_cv.R / 16b_model_training_final_panel_noMHC.R / 16d_nested_cv_reconciliation.R): plain, unpenalized glm(y ~ ., family = binomial) on every gene in the panel - no regularisation, so nothing to tune except how many folds the cross-validated AUC below uses. Class weighting (Reference/Comparison group settings above) applies here too.",
        fluidRow(
          column(4, numericInput(ns("lr_cv_folds"), "Cross-validation folds", value = 5, min = 3, max = 10, step = 1))
        )
      )
    }

    enet_params_box <- function() {
      mod_diagnostic_params_box(
        ns, "enet", "Elastic Net",
        "Alpha tuned by CV deviance; lambda.min used by default. Class weighting (Reference/Comparison group settings above) applies here too.",
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
        "mtry tuned by CV; ntree fixed. Class weighting (Reference/Comparison group settings above) applies here too.",
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
        "Cost tuned by CV; linear kernel by default. Class weighting (Reference/Comparison group settings above) applies here too.",
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
        p(class = "submodule-desc", "Gamma/degree only matter for Radial/Polynomial kernels above - both are ignored (no effect) while Kernel is set to Linear."),
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
      req(diag_has_run())
      switch(active_model_pill(),
        "Elastic Net" = enet_params_box(),
        "Random Forest" = rf_params_box(),
        "SVM" = svm_params_box(),
        lr_params_box()
      )
    })

    output$references_box_ui <- renderUI({
      req(diag_has_run())
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

    ## Female and male are saved into results$diagnostic independently as
    ## each finishes, mirroring Feature Selection's own save pattern.
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
      res_f <- tryCatch(diag_result_female(), error = function(e) NULL)
      res_m <- tryCatch(diag_result_male(), error = function(e) NULL)
      res_p <- tryCatch(diag_result_pooled(), error = function(e) NULL)
      status_row <- function(sex, r) {
        if (is.null(r)) {
          tags$li(icon("circle-minus", style = "color: #8A929C;"), sprintf(" %s diagnostic models - not run yet", sex))
        } else {
          best_label <- c("Logistic Regression", "Elastic Net", "Random Forest", "SVM")[which.max(c(mean(r$lr$cv_auc, na.rm = TRUE), mean(r$enet$cv_auc, na.rm = TRUE), mean(r$rf$cv_auc, na.rm = TRUE), mean(r$svm$cv_auc, na.rm = TRUE)))]
          tags$li(icon("check", style = "color: #1a9c5f;"), strong(sprintf(" %s completed: ", sex)),
                  sprintf("best CV-AUC = %s", best_label))
        }
      }
      tagList(
        p(class = "submodule-desc", style = "margin-bottom: 4px;", "Status:"),
        tags$ul(style = "padding-left: 18px; margin-bottom: 0; list-style: none;",
                status_row("Female", res_f), status_row("Male", res_m), status_row("Pooled (all)", res_p))
      )
    })

    ## One result line per sex, each shown ONLY inside that sex's own tab
    ## (mod_diagnostic_training_sex_panel() in the UI) - no shared box where
    ## the other sex's line would appear while you're looking at this one.
    diag_result_line <- function(sex_label, r) {
      if (is.null(r)) return(NULL)
      p(strong("Result: "),
        sprintf("%d genes, %d samples (%d vs %d), %s vs %s → CV-AUC logistic regression %.3f / elastic net %.3f / random forest %.3f / SVM %.3f.",
                r$n_input, r$n_samples, r$n_comp, r$n_ref, r$comp_group, r$ref_group,
                mean(r$lr$cv_auc, na.rm = TRUE), mean(r$enet$cv_auc, na.rm = TRUE), mean(r$rf$cv_auc, na.rm = TRUE), mean(r$svm$cv_auc, na.rm = TRUE)))
    }
    output$female_result_line <- renderUI({ diag_result_line("female", tryCatch(diag_result_female(), error = function(e) NULL)) })
    output$male_result_line <- renderUI({ diag_result_line("male", tryCatch(diag_result_male(), error = function(e) NULL)) })
    output$pooled_result_line <- renderUI({ diag_result_line("pooled", tryCatch(diag_result_pooled(), error = function(e) NULL)) })

    ## -----------------------------------------------------------------
    ## Per-technique, per-sex (or pooled) outputs.
    ## -----------------------------------------------------------------

    res_sex <- function(sex_label) reactive({
      fr <- switch(sex_label, female = diag_result_female, male = diag_result_male, pooled = diag_result_pooled)
      tryCatch(fr(), error = function(e) NULL)
    })

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
      ## Every box these outputs feed is written directly into the static UI
      ## (mod_diagnostic_ui above) rather than behind a req(res())-gated
      ## renderUI, so this note covers BOTH "never clicked yet" and "clicked
      ## but the run failed validation" - there's no separate state to tell
      ## them apart by once the box always exists.
      not_yet_note <- function() {
        div(class = "empty-note", icon("circle-info"), "Not run yet, or the last run failed validation - check above.")
      }

      lapply(DIAG_TECHNIQUES, function(tech) {
        key <- tech$key; label <- tech$label
        prefix <- paste0(sex_label, "_", key)

        ## ---- Model Training tab ----
        ## KPI tiles (shinydashboard::valueBox - this app's own established
        ## "SaaS stat card" pattern, already styled by www/custom.css and
        ## used the same way in mod_overview.R/mod_mr.R/mod_preprocessing.R)
        ## instead of a plain paragraph.
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
            ## Overfitting flag: a near-perfect training fit alongside a much
            ## weaker cross-validated AUC means the model memorised the
            ## training samples rather than learning real signal - classic
            ## p >= n behavior with a gene panel this large relative to the
            ## sample count, and Logistic Regression (no regularisation
            ## path at all) is the most exposed of the four techniques here.
            ## Purely informational and driven by this run's own AUC gap -
            ## no default/setting changes anywhere, so a well-sized panel
            ## (e.g. this project's own 25-40 gene MR-prioritised lists)
            ## never triggers it.
            if (isTRUE(rr$full_auc >= 0.95) && isTRUE((rr$full_auc - cv_mean) >= 0.25)) {
              div(class = "empty-note", style = "border-left: 3px solid #c0392b; padding-left: 8px;",
                  icon("triangle-exclamation"),
                  sprintf(
                    "Likely overfitting: Train AUC (%.3f) far exceeds the %d-fold CV AUC (%.3f) - with %d candidate genes on %d samples%s, trust the CV-AUC, not the training AUC. Consider a smaller, pre-selected panel (e.g. a Feature Selection consensus set) instead of a whole WGCNA module.",
                    rr$full_auc, length(rr$cv_auc), cv_mean, r$n_input, r$n_samples,
                    if (identical(rr$model_type, "lr")) " and no regularisation (plain logistic regression)" else ""
                  ))
            }
          )
        })
        output[[paste0(prefix, "_train_roc_plot")]] <- renderPlot({
          r <- res(); req(r)
          rr <- r[[key]]
          diag_roc_plot_pub(rr$roc_full, color = sex_color, title = sprintf("Training ROC — %s (%s)", label, tools::toTitleCase(sex_label)))
        })
        output[[paste0(prefix, "_train_cv_plot")]] <- renderPlot({
          r <- res(); req(r)
          rr <- r[[key]]
          df <- data.frame(fold = factor(seq_along(rr$cv_auc)), auc = rr$cv_auc)
          ggplot(df, aes(x = fold, y = auc)) +
            geom_col(fill = sex_color) +
            geom_hline(yintercept = mean(rr$cv_auc, na.rm = TRUE), linetype = "dashed", color = ARTHOMIX_COLORS$red) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(x = "Fold", y = "CV AUC (training)") + theme_arthomix(base_size = 12)
        })
        ## "Explore hyperparameter tuning": every value in the grid this
        ## technique actually searched (diag_fit_sex's alpha_search /
        ## mtry_search / cost_search, captured alongside the winner) plotted
        ## against its own CV score, with the value this run actually chose
        ## marked - not just the final number in the KPI tile above.
        output[[paste0(prefix, "_train_tuning_plot")]] <- renderPlot({
          r <- res(); req(r)
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
          r <- res(); req(r)
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
          r <- res(); req(r)
          rr <- r[[key]]
          req(isTRUE(rr$test$available))
          diag_roc_plot_pub(rr$test$roc, color = test_color,
                             title = sprintf("Test-split ROC — %s (%s)", label, tools::toTitleCase(sex_label)),
                             ci = c(rr$test$ci_lo, rr$test$auc, rr$test$ci_hi))
        })
        output[[paste0(prefix, "_test_table")]] <- DT::renderDataTable({
          r <- res(); req(r)
          DT::datatable(build_test_perf_table(r[[key]]), rownames = FALSE, width = "100%",
                        options = list(pageLength = 15, dom = "t", scrollX = TRUE, autoWidth = FALSE), class = "stripe hover compact")
        })
        output[[paste0(prefix, "_test_download")]] <- downloadHandler(
          filename = function() sprintf("%s_%s_testing_performance.csv", sex_label, key),
          content = function(file) write.csv(build_test_perf_table(res()[[key]]), file, row.names = FALSE)
        )
      })

      output[[paste0(sex_label, "_train_compare_table")]] <- DT::renderDataTable({
        r <- res(); req(r)
        df <- do.call(rbind, lapply(DIAG_TECHNIQUES, function(tech) {
          rr <- r[[tech$key]]
          data.frame(model = tech$label, training_auc = round(rr$full_auc, 3),
                     cv_mean_auc = round(mean(rr$cv_auc, na.rm = TRUE), 3), cv_sd_auc = round(stats::sd(rr$cv_auc, na.rm = TRUE), 3),
                     stringsAsFactors = FALSE)
        }))
        DT::datatable(df, rownames = FALSE, width = "100%", options = list(pageLength = 5, dom = "t", scrollX = TRUE, autoWidth = FALSE), class = "stripe hover compact")
      })
      output[[paste0(sex_label, "_test_compare_table")]] <- DT::renderDataTable({
        r <- res(); req(r)
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

      ## Per-gene ROC/AUC - same rendering for Training (gene_roc_train) and
      ## Testing (gene_roc_test), just a different field name.
      register_generoc_outputs <- function(mode, field) {
        plot_id <- paste0(sex_label, "_", mode, "_generoc_plot")
        table_id <- paste0(sex_label, "_", mode, "_generoc_table")
        download_id <- paste0(sex_label, "_", mode, "_generoc_download")
        hub_download_id <- paste0(sex_label, "_", mode, "_hub_download")
        auc_thr_id <- paste0(sex_label, "_", mode, "_hub_auc_thr")
        p_thr_id <- paste0(sex_label, "_", mode, "_hub_p_thr")

        ## `gr$p` is absent from older cached results$diagnostic runs (a
        ## session resumed from before this column existed) - NA rather than
        ## a hard error, so the AUC column still renders and only "hub"
        ## (which needs p) comes out empty instead of the whole table.
        gene_auc_df <- function(gr) {
          p <- if (!is.null(gr$p)) unname(gr$p) else rep(NA_real_, length(gr$genes))
          df <- data.frame(gene = gr$genes, auc = round(unname(gr$auc), 3), p = signif(p, 3), stringsAsFactors = FALSE)
          auc_thr <- input[[auc_thr_id]] %||% 0.85
          p_thr <- input[[p_thr_id]] %||% 0.05
          df$hub <- !is.na(df$auc) & !is.na(df$p) & df$auc >= auc_thr & df$p < p_thr
          df[order(-df$auc), ]
        }

        ## One small ROC panel PER GENE (not all genes overlaid on one axes,
        ## which becomes unreadable spaghetti past 3-4 genes and hid each
        ## gene's own AUC entirely) - same FPR/TPR orientation and 1:1
        ## aspect ratio as diag_roc_plot_pub()'s model-level ROC, with that
        ## gene's own AUC printed inside its panel. Capped to the top genes
        ## by AUC so a whole-module-sized panel doesn't render an
        ## unreadable wall of facets - the full per-gene table/download
        ## below always has every gene regardless of this cap.
        GENEROC_MAX_FACETS <- 24
        output[[plot_id]] <- renderPlot({
          r <- res(); req(r)
          gr <- r[[field]]
          auc_df <- gene_auc_df(gr)
          req(nrow(auc_df) > 0)
          top_genes <- head(auc_df$gene, GENEROC_MAX_FACETS)
          df <- do.call(rbind, lapply(top_genes, function(g) {
            rc <- gr$rocs[[g]]
            if (is.null(rc)) return(NULL)
            data.frame(gene = g, fpr = 1 - rc$specificities, tpr = rc$sensitivities, stringsAsFactors = FALSE)
          }))
          req(nrow(df) > 0)
          df$gene <- factor(df$gene, levels = top_genes)
          auc_lab <- auc_df[match(top_genes, auc_df$gene), c("gene", "auc")]
          auc_lab$gene <- factor(auc_lab$gene, levels = top_genes)
          auc_lab$label <- sprintf("AUC = %.3f", auc_lab$auc)
          n_col <- min(4, length(top_genes))
          plot_title <- if (nrow(auc_df) > GENEROC_MAX_FACETS) {
            sprintf("Top %d of %d genes by AUC - see the table below for the full panel", GENEROC_MAX_FACETS, nrow(auc_df))
          } else NULL
          ggplot(df, aes(x = fpr, y = tpr)) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#9CA3AF", linewidth = 0.5) +
            geom_line(color = sex_color, linewidth = 1) +
            geom_text(data = auc_lab, aes(x = 0.98, y = 0.04, label = label), inherit.aes = FALSE,
                      hjust = 1, vjust = 0, size = 3.1, color = "#1F2937") +
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
              panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA)
            )
        }, height = function() {
          r <- tryCatch(res(), error = function(e) NULL)
          if (is.null(r)) return(320)
          n <- min(nrow(gene_auc_df(r[[field]])), GENEROC_MAX_FACETS)
          if (n <= 0) return(320)
          n_col <- min(4, n)
          n_row <- ceiling(n / n_col)
          max(320, n_row * 230)
        })
        output[[table_id]] <- DT::renderDataTable({
          r <- res(); req(r)
          DT::datatable(gene_auc_df(r[[field]]), rownames = FALSE, width = "100%",
                        options = list(pageLength = 8, scrollX = TRUE, autoWidth = FALSE), class = "stripe hover compact") |>
            DT::formatStyle("hub", target = "row", backgroundColor = DT::styleEqual(c(TRUE, FALSE), c("#e6f4ea", "")))
        })
        output[[download_id]] <- downloadHandler(
          filename = function() sprintf("%s_%s_gene_auc.csv", sex_label, mode),
          content = function(file) write.csv(gene_auc_df(res()[[field]]), file, row.names = FALSE)
        )
        ## Just the rows currently passing this box's own AUC/P thresholds -
        ## the direct candidate list to carry into the next pipeline step
        ## (e.g. back into Feature Selection, or reported as this run's own
        ## "hub genes"), without hand-filtering the full CSV yourself.
        output[[hub_download_id]] <- downloadHandler(
          filename = function() sprintf("%s_%s_hub_genes.csv", sex_label, mode),
          content = function(file) {
            df <- gene_auc_df(res()[[field]])
            write.csv(df[df$hub, c("gene", "auc", "p")], file, row.names = FALSE)
          }
        )
      }
      register_generoc_outputs("train", "gene_roc_train")
      register_generoc_outputs("test", "gene_roc_test")
    }

    register_sex_model_outputs("female", res_sex("female"))
    register_sex_model_outputs("male", res_sex("male"))
    register_sex_model_outputs("pooled", res_sex("pooled"))
  })
}
