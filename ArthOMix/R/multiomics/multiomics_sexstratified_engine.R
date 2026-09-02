## R/multiomics/multiomics_sexstratified_engine.R
## Live sex-stratified DIABLO/Random-Forest engine - a parameterized,
## general-dataset port of Research_05_multiomics_sexstratified's own
## nested-CV pipeline (analyses/05_female_response/scripts/
## 11_diablo_response_sexstratified.R, .../13_random_forest_secondary_model.R,
## .../07_cross_analysis_summary/scripts/14_diablo_drugtype_sexstratified.R),
## so a user's own uploaded data (or a preloaded cell's own stored blocks)
## gets the SAME leakage-safe, in-fold-covariate-adjusted, nested 5x5-CV
## analysis that produced Table34/35/37/39 - not a generic freeform DIABLO
## run. Every default constant below (top-K per block, ncomp, folds,
## repeats, keepX bounds, ntree, seeds) is copied verbatim from those
## scripts; callers may override them via `params`, but the defaults are
## what reproduces the pipeline's own numbers when pointed at the pipeline's
## own cohorts.
##
## Kept deliberately separate from mi_diablo_run()
## (multiomics_integration_helpers.R) and mb_cv_roc()
## (multiomics_biomarker_helpers.R) - those are the app's existing freeform
## DIABLO engines (single outcome, no covariate adjustment, mixOmics's own
## tune.block.splsda()/perf() machinery) and are untouched by this file.
##
## Matrix convention throughout this file matches the rest of R/multiomics/:
## samples in rows, features in columns (`multi_dataset$layers`,
## mi_preloaded_cell_dataset()'s `fit$X`) - the reverse of the pipeline
## scripts' own features-in-rows convention. Every port below is transposed
## accordingly; the underlying arithmetic (which values get averaged/tested/
## imputed) is unchanged.

## ---------------------------------------------------------------------------
## Pipeline-derived default constants (verbatim from the scripts named above)
## ---------------------------------------------------------------------------

MSS_DEFAULTS <- list(
  top_expr = 50L, top_meth = 100L,            ## DIABLO per-block top-K (script 11/14)
  top_expr_rf = 50L, top_meth_rf = 100L,      ## RF integrated-blocks top-K (script 13, items 5-6)
  ncomp = 2L,
  folds = 5L, repeats = 5L,
  keepx_min = 5L, keepx_max = 15L,
  ntree = 500L,
  min_selected_diablo = 5L, min_selected_rf = 3L,
  min_train_rows_rf = 4L
)

## Below this many matched samples, a sex stratum is too small to run a
## nested-CV model at all - mirrors mcc_build_live()'s own `< 6` guard
## (mod_multi_concordance.R) so the same "too small, silently skip" rule
## applies everywhere a live sex split happens in this module.
MSS_MIN_STRATUM_N <- 6L

MSS_SEX_MODE_CHOICES <- c(
  "All (pooled)" = "pooled",
  "Female only" = "female",
  "Male only" = "male",
  "Female and Male separately" = "both"
)

MSS_ENGINE_CHOICES <- c("DIABLO (mixOmics::block.splsda)" = "diablo", "Random Forest (randomForest, ntree=500)" = "rf")

## ---------------------------------------------------------------------------
## 1. Cohort assembly - generic analogue of the pipeline's own get_cohort():
## intersects sample IDs present in BOTH chosen live matrices with the
## requested stratum's sample set. Never filters by outcome/covariate here -
## those enter only inside the CV loop as label/adjustment, exactly as in
## the pipeline.
## ---------------------------------------------------------------------------

mss_cohort_for_stratum <- function(expr_mat, meth_mat, sample_meta, sample_ids) {
  common <- intersect(rownames(expr_mat), rownames(meth_mat))
  ids <- intersect(common, sample_ids)
  if (length(ids) < MSS_MIN_STRATUM_N) {
    return(list(ok = FALSE, error = sprintf(
      "Only %d matched sample(s) in this stratum - at least %d are needed for a sex-stratified nested-CV run.",
      length(ids), MSS_MIN_STRATUM_N)))
  }
  list(ok = TRUE, ids = ids, expr = expr_mat[ids, , drop = FALSE], meth = meth_mat[ids, , drop = FALSE],
       meta = if (!is.null(sample_meta)) sample_meta[ids, , drop = FALSE] else NULL)
}

## ---------------------------------------------------------------------------
## 2. Leakage-safe per-fold imputation - direct port of the pipeline's
## impute_fold_matrix(): per-feature TRAIN-only mean, applied to both train
## and test. A feature whose training column is entirely NA is left as-is
## (matches the pipeline's own `if (is.na(train_mean[i])) next` guard).
## ---------------------------------------------------------------------------

mss_impute_fold_matrix <- function(X, train_idx, test_idx) {
  X_train <- X[train_idx, , drop = FALSE]
  X_test  <- X[test_idx, , drop = FALSE]
  train_mean <- colMeans(X_train, na.rm = TRUE)
  na_train <- which(colSums(is.na(X_train)) > 0)
  for (j in na_train) { if (is.na(train_mean[j])) next; idx <- is.na(X_train[, j]); X_train[idx, j] <- train_mean[j] }
  na_test <- which(colSums(is.na(X_test)) > 0)
  for (j in na_test) { if (is.na(train_mean[j])) next; idx <- is.na(X_test[, j]); X_test[idx, j] <- train_mean[j] }
  list(train = X_train, test = X_test)
}

## ---------------------------------------------------------------------------
## 3. In-fold feature selection - training-fold-only limma moderated-t
## ranking, covariate-adjusted (`~outcome + covariate`, coef=2 tests the
## outcome term). DIABLO uses topTable()'s own p-ranking (script 11); RF
## uses the pipeline's separate `rank_by_pvalue()` idiom (sorts the raw
## eBayes p-value column directly) - both are in-fold-only, never touch
## test_idx.
## ---------------------------------------------------------------------------

mss_limma_design <- function(outcome_train, covariate_train) {
  if (is.null(covariate_train)) stats::model.matrix(~outcome_train) else stats::model.matrix(~outcome_train + covariate_train)
}

mss_select_features_diablo <- function(expr_train, meth_train, outcome_train, covariate_train, top_expr, top_meth) {
  design <- mss_limma_design(outcome_train, covariate_train)
  sel_e <- tryCatch({
    fit <- limma::eBayes(limma::lmFit(t(expr_train), design))
    head(rownames(limma::topTable(fit, coef = 2, number = Inf, sort.by = "P")), top_expr)
  }, error = function(e) NULL)
  sel_m <- tryCatch({
    fit <- limma::eBayes(limma::lmFit(t(meth_train), design))
    head(rownames(limma::topTable(fit, coef = 2, number = Inf, sort.by = "P")), top_meth)
  }, error = function(e) NULL)
  list(expr = sel_e, meth = sel_m)
}

mss_rank_features_pvalue <- function(mat_train, outcome_train, covariate_train, top_k) {
  design <- mss_limma_design(outcome_train, covariate_train)
  tryCatch({
    fit <- limma::eBayes(limma::lmFit(t(mat_train), design))
    head(names(sort(fit$p.value[, 2])), top_k)
  }, error = function(e) NULL)
}

## ---------------------------------------------------------------------------
## 4. Per-fold model fits. Both return NULL to signal "skip this fold"
## (mirrors the pipeline's own `tryCatch(..., error=function(e) NULL)` /
## `next` guards) or `list(score, selected)` - `selected` (the in-fold chosen
## feature names) feeds the selection-frequency/stability summary below.
## ---------------------------------------------------------------------------

## DIABLO: select -> impute -> drop near-zero-variance training columns
## (pipeline lines 72-73/75-76) -> correlation-derived 2x2 design matrix ->
## keepX clamped to [keepx_min, keepx_max] -> block.splsda -> the held-out
## fold's WeightedVote max.dist score at the final component.
mss_diablo_fold <- function(expr, meth, outcome_full, covariate_full, train_idx, test_idx, params) {
  outcome_train <- outcome_full[train_idx]
  if (length(unique(as.character(outcome_train))) < 2) return(NULL)
  covariate_train <- if (!is.null(covariate_full)) covariate_full[train_idx] else NULL

  sel <- mss_select_features_diablo(expr[train_idx, , drop = FALSE], meth[train_idx, , drop = FALSE], outcome_train, covariate_train, params$top_expr, params$top_meth)
  if (is.null(sel$expr) || is.null(sel$meth) || length(sel$expr) < params$min_selected_diablo || length(sel$meth) < params$min_selected_diablo) return(NULL)

  imp_e <- mss_impute_fold_matrix(expr[, sel$expr, drop = FALSE], train_idx, test_idx)
  imp_m <- mss_impute_fold_matrix(meth[, sel$meth, drop = FALSE], train_idx, test_idx)

  keep_e <- apply(imp_e$train, 2, function(col) stats::var(col, na.rm = TRUE) > 0)
  keep_m <- apply(imp_m$train, 2, function(col) stats::var(col, na.rm = TRUE) > 0)
  Xtr <- list(expression = imp_e$train[, keep_e, drop = FALSE], methylation = imp_m$train[, keep_m, drop = FALSE])
  Xte <- list(expression = imp_e$test[, keep_e, drop = FALSE], methylation = imp_m$test[, keep_m, drop = FALSE])
  if (ncol(Xtr$expression) < 2 || ncol(Xtr$methylation) < 2) return(NULL)

  block_cor <- suppressWarnings(stats::cor(rowMeans(scale(Xtr$expression)), rowMeans(scale(Xtr$methylation))))
  if (is.na(block_cor)) block_cor <- 0
  design_mat <- matrix(abs(block_cor), 2, 2, dimnames = list(names(Xtr), names(Xtr)))
  diag(design_mat) <- 0

  keepX <- lapply(Xtr, function(x) max(params$keepx_min, min(params$keepx_max, ncol(x) - 1)))

  fit <- tryCatch(mixOmics::block.splsda(X = Xtr, Y = outcome_train, ncomp = params$ncomp, keepX = keepX, design = design_mat), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  pred <- tryCatch(stats::predict(fit, newdata = Xte), error = function(e) NULL)
  if (is.null(pred)) return(NULL)
  vote <- tryCatch(pred$WeightedVote$max.dist[, params$ncomp], error = function(e) NULL)
  if (is.null(vote)) return(NULL)

  positive_class <- levels(outcome_full)[2]
  list(score = as.numeric(vote == positive_class), selected = list(expr = sel$expr, meth = sel$meth))
}

## Random Forest: select (own top-K/ranking) -> impute -> concatenate
## expression+methylation features into one matrix, plus an explicit
## covariate column (pipeline lines ~ script 13:169-172 - DIABLO instead
## only adjusts for the covariate at the limma design-matrix stage, RF adds
## it as a model input) -> randomForest(ntree) -> predicted probability of
## the positive class.
mss_rf_fold <- function(expr, meth, outcome_full, covariate_full, train_idx, test_idx, params) {
  outcome_train <- outcome_full[train_idx]
  if (length(unique(as.character(outcome_train))) < 2) return(NULL)
  covariate_train <- if (!is.null(covariate_full)) covariate_full[train_idx] else NULL

  sel_e <- mss_rank_features_pvalue(expr[train_idx, , drop = FALSE], outcome_train, covariate_train, params$top_expr_rf)
  sel_m <- mss_rank_features_pvalue(meth[train_idx, , drop = FALSE], outcome_train, covariate_train, params$top_meth_rf)
  if (is.null(sel_e) || is.null(sel_m) || length(sel_e) < params$min_selected_rf || length(sel_m) < params$min_selected_rf) return(NULL)

  imp_e <- mss_impute_fold_matrix(expr[, sel_e, drop = FALSE], train_idx, test_idx)
  imp_m <- mss_impute_fold_matrix(meth[, sel_m, drop = FALSE], train_idx, test_idx)

  Xtr <- cbind(imp_e$train, imp_m$train)
  Xte <- cbind(imp_e$test, imp_m$test)
  if (!is.null(covariate_full)) {
    Xtr <- cbind(Xtr, .covariate = as.integer(covariate_train))
    Xte <- cbind(Xte, .covariate = as.integer(covariate_full[test_idx]))
  }
  if (nrow(Xtr) < params$min_train_rows_rf) return(NULL)

  fit <- tryCatch(randomForest::randomForest(x = Xtr, y = outcome_train, ntree = params$ntree), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  positive_class <- levels(outcome_full)[2]
  prob <- tryCatch(stats::predict(fit, Xte, type = "prob")[, positive_class], error = function(e) NULL)
  if (is.null(prob)) return(NULL)
  list(score = as.numeric(prob), selected = list(expr = sel_e, meth = sel_m))
}

## ---------------------------------------------------------------------------
## 5. Nested CV orchestrator - caret::createFolds(outcome, k), OUTER_REPEATS
## repeats, DIABLO's own seed scheme (3000+rep_i, script 11/14) reused for
## both engines so a Female/Male "both separately" run draws the same fold
## partitions regardless of which engine is selected. OOF predictions are
## averaged per patient (rowMeans(na.rm=TRUE)), then scored with
## pROC::roc(..., direction="<") + ci.auc() exactly as the pipeline does.
## Requires a strictly two-class outcome - the pipeline's own scope
## (responder/non-responder, or drug A/B) never modeled more than two
## classes with this design.
## ---------------------------------------------------------------------------

mss_nested_cv <- function(expr, meth, outcome, covariate = NULL, engine = c("diablo", "rf"), params = list()) {
  engine <- match.arg(engine)
  p <- utils::modifyList(MSS_DEFAULTS, params %||% list())
  ids <- rownames(expr)
  n <- length(ids)

  outcome_levels <- levels(droplevels(factor(outcome[ids])))
  if (length(outcome_levels) != 2) return(list(ok = FALSE, error = "This engine requires a two-class outcome (e.g. responder/non-responder)."))
  outcome_full <- factor(outcome[ids], levels = outcome_levels)

  covariate_full <- NULL
  if (!is.null(covariate)) {
    cov_try <- factor(covariate[ids])
    if (nlevels(cov_try) >= 2) covariate_full <- cov_try
  }

  min_class_n <- min(table(outcome_full))
  if (min_class_n < 3) return(list(ok = FALSE, error = sprintf("The smaller outcome class has only %d sample(s) - at least 3 are needed per class.", min_class_n)))

  folds_req <- max(2L, min(as.integer(p$folds), min_class_n))
  repeats_req <- max(1L, as.integer(p$repeats))

  oof <- matrix(NA_real_, nrow = n, ncol = repeats_req, dimnames = list(ids, NULL))
  fold_selected <- list()

  for (rep_i in seq_len(repeats_req)) {
    set.seed(3000 + rep_i)
    fold_test_ids <- caret::createFolds(outcome_full, k = folds_req, list = TRUE)
    for (f in fold_test_ids) {
      test_idx <- ids[f]
      train_idx <- setdiff(ids, test_idx)
      if (length(unique(as.character(outcome_full[train_idx]))) < 2) next
      if (!is.null(covariate_full) && length(unique(as.character(covariate_full[train_idx]))) < 2) next

      res <- if (identical(engine, "diablo")) {
        mss_diablo_fold(expr, meth, outcome_full, covariate_full, train_idx, test_idx, p)
      } else {
        mss_rf_fold(expr, meth, outcome_full, covariate_full, train_idx, test_idx, p)
      }
      if (!is.null(res)) {
        oof[test_idx, rep_i] <- res$score
        fold_selected[[length(fold_selected) + 1]] <- res$selected
      }
    }
  }

  patient_score <- rowMeans(oof, na.rm = TRUE)
  ok_idx <- !is.na(patient_score)
  if (sum(ok_idx) < 6) return(list(ok = FALSE, error = "Too few patients received an out-of-fold score to compute performance - this configuration may be too small or too unstable for nested CV."))

  y_num <- as.numeric(outcome_full[ok_idx] == outcome_levels[2])
  roc_obj <- tryCatch(pROC::roc(y_num, patient_score[ok_idx], quiet = TRUE, direction = "<"), error = function(e) NULL)
  if (is.null(roc_obj)) return(list(ok = FALSE, error = "AUROC could not be computed for this configuration."))
  ci <- tryCatch(as.numeric(pROC::ci.auc(roc_obj)), error = function(e) c(NA_real_, NA_real_, NA_real_))

  list(
    ok = TRUE,
    performance = data.frame(
      n = sum(ok_idx), auroc = as.numeric(roc_obj$auc), ci_lo = ci[1], ci_hi = ci[3],
      excludes_chance = isTRUE(!is.na(ci[1]) && !is.na(ci[3]) && (ci[1] > 0.5 || ci[3] < 0.5)),
      engine = engine, positive_class = outcome_levels[2], stringsAsFactors = FALSE
    ),
    scores = data.frame(patient_id = ids, score = as.numeric(patient_score), outcome = as.character(outcome_full), stringsAsFactors = FALSE),
    fold_selected_features = fold_selected,
    params = p, outcome_levels = outcome_levels
  )
}

## Per-feature selection frequency across every fold that actually fit a
## model (out of the up-to- folds*repeats attempts) - the stability metric
## Biomarker Discovery's "Sex-Stratified" tab surfaces instead of
## re-deriving one from mixOmics::perf()'s own stability structure (which
## this engine, unlike mi_diablo_run(), never calls).
mss_selection_frequency <- function(fold_selected_features, block = c("expr", "meth")) {
  block <- match.arg(block)
  n_folds <- length(fold_selected_features)
  if (n_folds == 0) return(NULL)
  feats <- unlist(lapply(fold_selected_features, function(s) s[[block]]))
  if (length(feats) == 0) return(NULL)
  tab <- table(feats)
  df <- data.frame(feature = names(tab), times_selected = as.integer(tab), n_folds = n_folds,
                    selection_frequency = as.numeric(tab) / n_folds, stringsAsFactors = FALSE)
  df[order(-df$selection_frequency), , drop = FALSE]
}

## ---------------------------------------------------------------------------
## 6. Full-cohort descriptive panel - no CV, purely for the panel/loadings
## artifact (mirrors the pipeline's own separate full-cohort refit, script
## 11 lines 122-144). DIABLO's panel is `selectVar()` component-1 loadings,
## exactly as the pipeline reports it (Table35's own shape). Random Forest
## has no equivalent pipeline artifact - its panel here (MeanDecreaseGini
## importance from a full-cohort refit) is a live extension, flagged as such
## in the returned `note` for the UI to display verbatim.
## ---------------------------------------------------------------------------

mss_full_cohort_panel <- function(expr, meth, outcome, covariate = NULL, engine = c("diablo", "rf"), params = list()) {
  engine <- match.arg(engine)
  p <- utils::modifyList(MSS_DEFAULTS, params %||% list())
  ids <- rownames(expr)

  outcome_levels <- levels(droplevels(factor(outcome[ids])))
  if (length(outcome_levels) != 2) return(list(ok = FALSE, error = "This engine requires a two-class outcome."))
  outcome_full <- factor(outcome[ids], levels = outcome_levels)
  covariate_full <- NULL
  if (!is.null(covariate)) {
    cov_try <- factor(covariate[ids])
    if (nlevels(cov_try) >= 2) covariate_full <- cov_try
  }

  design <- mss_limma_design(outcome_full, covariate_full)

  impute_full <- function(mat) {
    means <- colMeans(mat, na.rm = TRUE)
    for (j in which(colSums(is.na(mat)) > 0)) { if (is.na(means[j])) next; mat[is.na(mat[, j]), j] <- means[j] }
    mat
  }
  ## Both blocks get the same train-only-style mean imputation the per-fold
  ## path already applies (mss_impute_fold_matrix) - expr was previously left
  ## raw here while meth was imputed, so any NA in the expression block made
  ## this descriptive panel fail (randomForest rejects NA predictors outright)
  ## even when the CV performance path succeeded via its own fold-level
  ## imputation, an inconsistent and confusing "performance shown, panel
  ## silently missing" outcome for the same dataset.
  expr <- impute_full(expr)
  meth_full <- impute_full(meth)

  if (identical(engine, "diablo")) {
    sel_e <- tryCatch(head(rownames(limma::topTable(limma::eBayes(limma::lmFit(t(expr), design)), coef = 2, number = Inf, sort.by = "P")), p$top_expr), error = function(e) NULL)
    sel_m <- tryCatch(head(rownames(limma::topTable(limma::eBayes(limma::lmFit(t(meth_full), design)), coef = 2, number = Inf, sort.by = "P")), p$top_meth), error = function(e) NULL)
    if (is.null(sel_e) || is.null(sel_m) || length(sel_e) < 2 || length(sel_m) < 2) return(list(ok = FALSE, error = "Feature selection failed on the full cohort."))

    Xfull <- list(expression = expr[, sel_e, drop = FALSE], methylation = meth_full[, sel_m, drop = FALSE])
    block_cor <- suppressWarnings(stats::cor(rowMeans(scale(Xfull$expression)), rowMeans(scale(Xfull$methylation))))
    if (is.na(block_cor)) block_cor <- 0
    design_mat <- matrix(abs(block_cor), 2, 2, dimnames = list(names(Xfull), names(Xfull)))
    diag(design_mat) <- 0
    keepX <- lapply(Xfull, function(x) max(p$keepx_min, min(p$keepx_max, ncol(x) - 1)))

    fit <- tryCatch(mixOmics::block.splsda(X = Xfull, Y = outcome_full, ncomp = p$ncomp, keepX = keepX, design = design_mat), error = function(e) NULL)
    if (is.null(fit)) return(list(ok = FALSE, error = "DIABLO full-cohort refit failed."))
    panel <- do.call(rbind, lapply(names(Xfull), function(blk) {
      sv <- tryCatch(mixOmics::selectVar(fit, comp = 1, block = blk), error = function(e) NULL)
      w <- sv[[blk]]$value
      if (is.null(w)) return(NULL)
      data.frame(view = blk, feature = rownames(w), loading = w$value.var, stringsAsFactors = FALSE)
    }))
    return(list(ok = TRUE, panel = panel, fit = fit, note = NULL))
  }

  ## engine == "rf"
  sel_e <- tryCatch(head(names(sort(limma::eBayes(limma::lmFit(t(expr), design))$p.value[, 2])), p$top_expr_rf), error = function(e) NULL)
  sel_m <- tryCatch(head(names(sort(limma::eBayes(limma::lmFit(t(meth_full), design))$p.value[, 2])), p$top_meth_rf), error = function(e) NULL)
  if (is.null(sel_e) || is.null(sel_m) || length(sel_e) < 2 || length(sel_m) < 2) return(list(ok = FALSE, error = "Feature selection failed on the full cohort."))

  Xfull <- cbind(expr[, sel_e, drop = FALSE], meth_full[, sel_m, drop = FALSE])
  if (!is.null(covariate_full)) Xfull <- cbind(Xfull, .covariate = as.integer(covariate_full))
  fit <- tryCatch(randomForest::randomForest(x = Xfull, y = outcome_full, ntree = p$ntree, importance = TRUE), error = function(e) NULL)
  if (is.null(fit)) return(list(ok = FALSE, error = "Random Forest full-cohort refit failed."))
  imp <- randomForest::importance(fit)
  panel <- data.frame(
    view = ifelse(rownames(imp) %in% sel_e, "expression", ifelse(rownames(imp) %in% sel_m, "methylation", "covariate")),
    feature = rownames(imp), importance = imp[, "MeanDecreaseGini"], stringsAsFactors = FALSE
  )
  panel <- panel[order(-panel$importance), , drop = FALSE]
  list(ok = TRUE, panel = panel, fit = fit,
       note = "Random Forest panel = MeanDecreaseGini importance (full-cohort refit) - a live extension, not in the original pipeline.")
}

## ---------------------------------------------------------------------------
## 7. Top-level orchestrator - runs one or more sex strata end to end
## (nested-CV performance + full-cohort descriptive panel). This is the one
## function the Integration/Biomarker "Sex-Stratified" tabs call directly;
## everything above is an implementation detail of this function.
## ---------------------------------------------------------------------------

mss_run_stratified <- function(expr_mat, meth_mat, sample_meta, outcome_col, covariate_col = NULL,
                                sex_mode = c("pooled", "female", "male", "both"), engine = c("diablo", "rf"), params = list()) {
  engine <- match.arg(engine)
  sex_mode <- match.arg(sex_mode)
  if (is.null(sample_meta) || !outcome_col %in% colnames(sample_meta)) return(list(ok = FALSE, error = "Select a valid outcome column."))

  outcome_vec <- stats::setNames(sample_meta[[outcome_col]], rownames(sample_meta))
  covariate_vec <- if (!is.null(covariate_col) && nzchar(covariate_col %||% "") && covariate_col %in% colnames(sample_meta)) {
    stats::setNames(sample_meta[[covariate_col]], rownames(sample_meta))
  } else NULL

  common <- intersect(rownames(expr_mat), rownames(meth_mat))
  sex_cands <- multi_sex_candidates(sample_meta)
  sex_col <- if (length(sex_cands) == 0) NULL else sex_cands[1]
  groups <- if (!is.null(sex_col)) multi_sex_groups(sample_meta, sex_col, common) else NULL

  strata <- if (identical(sex_mode, "both")) {
    if (is.null(groups) || length(groups) < 2) return(list(ok = FALSE, error = "No sex/gender column detected (or fewer than two sexes represented among matched samples) - cannot run Female and Male separately."))
    groups
  } else if (sex_mode %in% c("female", "male")) {
    if (is.null(groups) || !sex_mode %in% tolower(names(groups))) return(list(ok = FALSE, error = sprintf("No %s samples found via the detected sex/gender column.", sex_mode)))
    stats::setNames(list(groups[[names(groups)[tolower(names(groups)) == sex_mode][1]]]), sex_mode)
  } else {
    stats::setNames(list(common), "pooled")
  }

  results <- list()
  for (label in names(strata)) {
    cohort <- mss_cohort_for_stratum(expr_mat, meth_mat, sample_meta, intersect(common, strata[[label]]))
    if (!isTRUE(cohort$ok)) { results[[label]] <- list(ok = FALSE, error = cohort$error, n = length(intersect(common, strata[[label]]))); next }
    cv <- mss_nested_cv(cohort$expr, cohort$meth, outcome_vec, covariate_vec, engine = engine, params = params)
    panel <- if (isTRUE(cv$ok)) mss_full_cohort_panel(cohort$expr, cohort$meth, outcome_vec, covariate_vec, engine = engine, params = params) else NULL
    results[[label]] <- list(ok = isTRUE(cv$ok), error = cv$error, cv = cv, panel = panel, n = length(cohort$ids))
  }

  ok_strata <- Filter(function(r) isTRUE(r$ok), results)
  if (length(ok_strata) == 0) {
    first_error <- results[[1]]$error %||% "No stratum produced a usable result."
    return(list(ok = FALSE, error = first_error, strata = results))
  }

  performance <- do.call(rbind, lapply(names(ok_strata), function(label) cbind(stratum = label, ok_strata[[label]]$cv$performance, stringsAsFactors = FALSE)))
  panels <- do.call(rbind, lapply(names(ok_strata), function(label) {
    pr <- ok_strata[[label]]$panel
    if (is.null(pr) || !isTRUE(pr$ok)) return(NULL)
    cbind(stratum = label, pr$panel, stringsAsFactors = FALSE)
  }))
  panel_note <- Find(Negate(is.null), lapply(ok_strata, function(r) r$panel$note))

  list(ok = TRUE, engine = engine, sex_col = sex_col, strata = results, performance = performance, panels = panels,
       panels_wide = mss_panel_wide_by_sex(panels), panel_note = panel_note)
}

## ---------------------------------------------------------------------------
## 8. Side-by-side sex comparison - pivots the combined `panels` table (long
## format: stratum, view, feature, loading/importance) into one row per
## feature with a separate column per stratum actually run (e.g. female,
## male, pooled), so a user can see at a glance which sex(es) a given
## feature was a candidate biomarker in, rather than scrolling one merged
## long table sorted by stratum.
## ---------------------------------------------------------------------------

mss_panel_wide_by_sex <- function(panels) {
  if (is.null(panels) || nrow(panels) == 0) return(NULL)
  value_col <- if ("loading" %in% colnames(panels)) "loading" else if ("importance" %in% colnames(panels)) "importance" else NULL
  if (is.null(value_col)) return(NULL)
  long <- panels[, c("view", "feature", "stratum", value_col)]
  wide <- tryCatch(
    as.data.frame(tidyr::pivot_wider(long, id_cols = c("view", "feature"), names_from = "stratum", values_from = dplyr::all_of(value_col))),
    error = function(e) NULL
  )
  if (is.null(wide)) return(NULL)
  wide[order(wide$view, wide$feature), , drop = FALSE]
}
