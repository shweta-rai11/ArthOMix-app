## R/transcriptomics/03_Preprocessing_Batch_Correction/mod_preprocessing_explore.R
## Data Exploration tab: standalone EDA/QC module (own raw-data upload, own
## namespace "eda") - independent of the shared `dataset` reactiveValues,

EDA_MAX_POOLED_VALUES   <- 200000  # pooled histogram/density/Q-Q sampling cap
EDA_MAX_SHAPIRO_N       <- 5000    # stats::shapiro.test()'s own hard cap is 5000
EDA_MAX_VARIANCE_FEATURES <- 2000  # top-variance features used for PCA/correlation/clustering
EDA_MAX_MEANVAR_POINTS  <- 5000    # points drawn on the mean-variance scatter
EDA_MAX_VIOLIN_FEATURES <- 3000    # features sampled per violin plot
EDA_MAX_DENSITY_SAMPLES <- 200     # per-sample density overlay / violin sample cap

eda_parse_upload <- function(datapath, filename) {
  df <- tryCatch(
    as.data.frame(data.table::fread(datapath, showProgress = FALSE,
                                     na.strings = c("NA", "", "NaN", "null", "NULL", "#N/A"))),
    error = function(e) NULL
  )
  if (is.null(df)) {
    return(list(ok = FALSE, error = "The file could not be parsed as a delimited table. Please upload a CSV, TSV, or TXT file with a consistent delimiter."))
  }
  if (nrow(df) == 0 || ncol(df) < 2) {
    return(list(ok = FALSE, error = "The file has no data rows, or fewer than two columns. Expecting a feature-by-sample matrix: one identifier column plus at least one numeric sample column."))
  }

  first_col <- df[[1]]
  first_num <- suppressWarnings(as.numeric(as.character(first_col)))
  id_is_char <- mean(is.na(first_num)) > 0.5
  if (id_is_char) {
    ids <- as.character(first_col)
    rest <- df[, -1, drop = FALSE]
    id_col_name <- colnames(df)[1]
  } else {
    ids <- paste0("row_", seq_len(nrow(df)))
    rest <- df
    id_col_name <- NULL
  }
  if (ncol(rest) == 0) {
    return(list(ok = FALSE, error = "No columns remain after the identifier column - at least one numeric sample column is required."))
  }

  is_num_col <- vapply(rest, function(col) {
    if (is.numeric(col)) return(TRUE)
    v <- suppressWarnings(as.numeric(as.character(col)))
    blank <- is.na(col) | (is.character(col) & trimws(as.character(col)) == "")
    mean(is.na(v) & !blank) < 0.2
  }, logical(1))

  if (!any(is_num_col)) {
    return(list(ok = FALSE, error = "No numeric sample columns were found. Expecting a feature-by-sample matrix with numeric expression/intensity values."))
  }

  numeric_part <- rest[, is_num_col, drop = FALSE]
  expr <- vapply(numeric_part, function(col) suppressWarnings(as.numeric(as.character(col))), numeric(nrow(numeric_part)))
  if (is.null(dim(expr))) expr <- matrix(expr, nrow = 1, dimnames = list(NULL, names(expr)))
  rownames(expr) <- ids
  colnames(expr) <- colnames(numeric_part)

  if (!any(is.finite(expr))) {
    return(list(ok = FALSE, error = "The detected numeric columns contain no usable (finite) values."))
  }

  list(
    ok = TRUE, expr = expr, id_col_name = id_col_name,
    nonnumeric_cols = colnames(rest)[!is_num_col],
    n_rows_orig = nrow(df), n_cols_orig = ncol(df), filename = filename
  )
}

eda_skewness <- function(x) {
  x <- x[is.finite(x)]; n <- length(x)
  if (n < 3) return(NA_real_)
  m <- mean(x); s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(0)
  (sum((x - m)^3) / n) / s^3
}

eda_kurtosis <- function(x) {
  x <- x[is.finite(x)]; n <- length(x)
  if (n < 4) return(NA_real_)
  m <- mean(x); s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(0)
  (sum((x - m)^4) / n) / s^4 - 3
}

eda_robust_z <- function(x) {
  med <- stats::median(x, na.rm = TRUE)
  mad <- stats::mad(x, na.rm = TRUE)
  if (!is.finite(mad) || mad == 0) return(rep(0, length(x)))
  0.6745 * (x - med) / mad
}

eda_skew_label <- function(skew) {
  if (is.na(skew)) return("Unknown (not enough data)")
  a <- abs(skew)
  if (a < 0.5) "Approximately symmetric"
  else if (a < 1) "Moderately skewed"
  else "Strongly skewed"
}

eda_overview <- function(parsed) {
  m <- parsed$expr
  vals <- as.numeric(m)
  finite_vals <- vals[is.finite(vals)]
  ids <- rownames(m)

  fp <- apply(m, 2, function(col) {
    col <- col[is.finite(col)]
    if (length(col) == 0) return("empty")
    paste(round(c(mean(col), stats::sd(col), stats::quantile(col, probs = c(0, .25, .5, .75, 1))), 6), collapse = "_")
  })

  row_var <- apply(m, 1, stats::var, na.rm = TRUE)
  row_var_valid <- row_var[is.finite(row_var)]
  n_constant <- sum(row_var_valid == 0)
  positive_var <- row_var_valid[row_var_valid > 0]
  nzv_cutoff <- if (length(positive_var) > 0) stats::quantile(positive_var, 0.01) else NA_real_
  n_near_zero_var <- if (is.na(nzv_cutoff)) 0L else sum(row_var_valid > 0 & row_var_valid <= nzv_cutoff)

  list(
    n_features = nrow(m), n_samples = ncol(m),
    n_numeric_cols = ncol(m), n_nonnumeric_cols = length(parsed$nonnumeric_cols),
    n_missing = sum(is.na(vals)), pct_missing = 100 * sum(is.na(vals)) / length(vals),
    n_infinite = sum(is.infinite(vals)),
    n_duplicated_features = sum(duplicated(ids)),
    n_duplicated_samples = sum(duplicated(fp)),
    n_constant_features = n_constant,
    n_near_zero_var_features = n_near_zero_var,
    mean = mean(finite_vals), median = stats::median(finite_vals), sd = stats::sd(finite_vals),
    var = stats::var(finite_vals),
    min = min(finite_vals), max = max(finite_vals),
    q1 = unname(stats::quantile(finite_vals, 0.25)), q3 = unname(stats::quantile(finite_vals, 0.75)),
    iqr = stats::IQR(finite_vals)
  )
}

eda_descriptive_stats <- function(m, margin) {
  agg <- function(v) {
    v <- v[is.finite(v)]
    n <- length(v)
    if (n < 2) {
      return(c(n = n, mean = NA_real_, median = NA_real_, sd = NA_real_, var = NA_real_,
                min = NA_real_, max = NA_real_, q1 = NA_real_, q3 = NA_real_, iqr = NA_real_,
                range = NA_real_, skewness = NA_real_, kurtosis = NA_real_, cv = NA_real_))
    }
    mn <- mean(v); sdv <- stats::sd(v)
    c(n = n, mean = mn, median = stats::median(v), sd = sdv, var = stats::var(v),
      min = min(v), max = max(v), q1 = unname(stats::quantile(v, .25)), q3 = unname(stats::quantile(v, .75)),
      iqr = stats::IQR(v), range = max(v) - min(v),
      skewness = eda_skewness(v), kurtosis = eda_kurtosis(v),
      cv = if (is.finite(mn) && mn != 0) sdv / abs(mn) else NA_real_)
  }
  out <- t(apply(m, margin, agg))
  df <- data.frame(id = if (margin == 1) rownames(m) else colnames(m), out, row.names = NULL, check.names = FALSE)
  colnames(df)[1] <- if (margin == 1) "feature" else "sample"
  df
}

eda_normality_summary <- function(m) {
  vals <- as.numeric(m); vals <- vals[is.finite(vals)]
  skew <- eda_skewness(vals); kurt <- eda_kurtosis(vals)
  shapiro <- NULL
  shapiro_n <- NA_integer_
  if (length(vals) >= 20) {
    shapiro_n <- min(length(vals), EDA_MAX_SHAPIRO_N)
    shapiro <- tryCatch(stats::shapiro.test(sample(vals, shapiro_n)), error = function(e) NULL)
  }
  list(
    skewness = skew, kurtosis = kurt, label = eda_skew_label(skew),
    shapiro_W = if (!is.null(shapiro)) unname(shapiro$statistic) else NA_real_,
    shapiro_p = if (!is.null(shapiro)) shapiro$p.value else NA_real_,
    shapiro_n = shapiro_n,
    qq_values = if (length(vals) > EDA_MAX_POOLED_VALUES) sample(vals, EDA_MAX_POOLED_VALUES) else vals
  )
}

eda_normalization_assessment <- function(expr) {
  m <- as.matrix(expr)
  finite_vals <- m[is.finite(m)]
  dt <- detect_expr_data_type(m)
  diag <- summarize_norm_diagnostics(m)
  differs <- needs_quantile_norm(diag)
  frac_integer <- mean(abs(finite_vals - round(finite_vals)) < 1e-6)
  has_negative <- any(finite_vals < 0)
  skew <- eda_skewness(finite_vals)
  evidence <- character(0)

  if (identical(dt, "counts")) {
    verdict <- "not_normalized"; label <- "Likely not normalized (raw counts)"
    evidence <- c(evidence,
      sprintf("%.0f%% of finite values are at or near integers and non-negative, with a maximum value of %s - the signature of raw sequencing read/count data rather than a continuous, transformed measurement.",
              frac_integer * 100, format(round(max(finite_vals)), big.mark = ",")),
      "Raw counts of this kind have not yet been adjusted for library size or composition (e.g. via TMM, DESeq2 size factors, or CPM).")
  } else if (has_negative) {
    verdict <- "normalized"; label <- "Likely normalized / transformed"
    evidence <- c(evidence, "Negative values are present, which is only possible after a log-ratio, z-score, or similarly centered transform - raw intensity or count data cannot be negative.")
    if (!differs) {
      evidence <- c(evidence, "Per-sample medians and interquartile ranges also agree closely across samples, consistent with a normalized matrix.")
    } else {
      evidence <- c(evidence, "However, per-sample medians and/or interquartile ranges still show some spread - normalization may be incomplete, or a residual batch-like effect may remain.")
    }
  } else if (!differs) {
    verdict <- "normalized"; label <- "Likely normalized"
    evidence <- c(evidence,
      "Per-sample medians and interquartile ranges are closely comparable across samples (their between-sample spread is small) - the expected signature of a normalized expression matrix.",
      sprintf("The overall value range (up to %.2f) and the %s pooled distribution shape (skewness = %.2f) are also consistent with continuous, transformed expression data rather than raw counts.",
              max(finite_vals), tolower(eda_skew_label(skew)), skew))
  } else {
    verdict <- "inconclusive"; label <- "Inconclusive"
    evidence <- c(evidence,
      "Per-sample medians and/or interquartile ranges differ noticeably between samples. This can indicate data that has not yet been normalized, but it can also reflect genuine biological or batch variation in an already-normalized dataset.",
      sprintf("Values are continuous rather than count-like (%.0f%% near-integer) and the overall range (%.2f to %.2f) does not clearly indicate raw counts either.",
              frac_integer * 100, min(finite_vals), max(finite_vals)))
  }

  list(verdict = verdict, label = label, evidence = evidence, detected_type = dt,
       between_sample_differs = differs, frac_integer = frac_integer, has_negative = has_negative, diag = diag)
}

eda_impute_median <- function(m) {
  mm <- m
  mm[!is.finite(mm)] <- NA
  if (anyNA(mm)) {
    row_med <- apply(mm, 1, stats::median, na.rm = TRUE)
    na_idx <- which(is.na(mm), arr.ind = TRUE)
    if (nrow(na_idx) > 0) mm[na_idx] <- row_med[na_idx[, 1]]
  }
  if (anyNA(mm)) mm[is.na(mm)] <- stats::median(mm, na.rm = TRUE)
  mm
}

eda_prep_for_structure <- function(m, max_features = EDA_MAX_VARIANCE_FEATURES) {
  mm <- eda_impute_median(m)
  row_var <- apply(mm, 1, stats::var, na.rm = TRUE)
  row_var[!is.finite(row_var)] <- 0
  n_keep <- min(max_features, sum(row_var > 0))
  if (n_keep < 3) return(NULL)
  keep <- order(row_var, decreasing = TRUE)[seq_len(n_keep)]
  list(sub = mm[keep, , drop = FALSE], n_features_used = n_keep, n_features_total = nrow(m))
}

eda_pca <- function(m) {
  prep <- eda_prep_for_structure(m)
  if (is.null(prep) || ncol(prep$sub) < 3) return(NULL)
  pr <- tryCatch(stats::prcomp(t(prep$sub), scale. = TRUE, center = TRUE), error = function(e) NULL)
  if (is.null(pr) || ncol(pr$x) < 2) return(NULL)
  var_exp <- (pr$sdev^2) / sum(pr$sdev^2) * 100
  list(scores = pr$x, var_exp = var_exp, n_features_used = prep$n_features_used, n_features_total = prep$n_features_total)
}

eda_sample_correlation <- function(m) {
  prep <- eda_prep_for_structure(m)
  if (is.null(prep) || ncol(prep$sub) < 2) return(NULL)
  list(cor = stats::cor(prep$sub, use = "pairwise.complete.obs"),
       n_features_used = prep$n_features_used, n_features_total = prep$n_features_total)
}

eda_hclust <- function(cor_mat) stats::hclust(stats::as.dist(1 - cor_mat), method = "average")

eda_sample_outliers <- function(m, pca) {
  qc <- tryCatch(compute_sample_qc(eda_impute_median(m)), error = function(e) NULL)
  if (is.null(qc)) {
    qc <- data.frame(sample = colnames(m), signal = NA_real_, detected = NA_real_, mean_cor = NA_real_,
                      flag_signal = FALSE, flag_detected = FALSE, flag_cor = FALSE, stringsAsFactors = FALSE)
  }
  qc$pca_distance <- NA_real_; qc$flag_pca <- FALSE
  if (!is.null(pca) && ncol(pca$scores) >= 2) {
    scores <- pca$scores[, 1:2, drop = FALSE]
    center <- colMeans(scores)
    dist <- sqrt(rowSums((scores - matrix(center, nrow(scores), ncol(scores), byrow = TRUE))^2))
    z <- eda_robust_z(dist)
    idx <- match(qc$sample, rownames(scores))
    qc$pca_distance <- dist[idx]
    qc$flag_pca <- ifelse(is.na(idx), FALSE, abs(z[idx]) > 3.5)
  }
  qc$n_flags <- rowSums(qc[, c("flag_signal", "flag_detected", "flag_cor", "flag_pca")], na.rm = TRUE)
  qc
}

eda_feature_outliers <- function(m, desc_df) {
  miss_pct <- rowMeans(is.na(m)) * 100
  var_z <- eda_robust_z(desc_df$var)
  data.frame(
    feature = desc_df$feature, variance = desc_df$var, skewness = desc_df$skewness, pct_missing = miss_pct,
    flag_extreme_variance = is.finite(var_z) & abs(var_z) > 3.5,
    flag_extreme_skew = is.finite(desc_df$skewness) & abs(desc_df$skewness) > 2,
    flag_high_missing = miss_pct > 20,
    stringsAsFactors = FALSE
  )
}

eda_missingness <- function(m) {
  by_sample <- data.frame(sample = colnames(m), n_missing = colSums(is.na(m)),
                           pct_missing = 100 * colMeans(is.na(m)), stringsAsFactors = FALSE)
  fm <- rowMeans(is.na(m)) * 100
  by_feature_summary <- data.frame(
    bucket = c("0%", ">0-5%", ">5-20%", ">20-50%", ">50%"),
    n_features = c(sum(fm == 0), sum(fm > 0 & fm <= 5), sum(fm > 5 & fm <= 20), sum(fm > 20 & fm <= 50), sum(fm > 50)),
    stringsAsFactors = FALSE
  )
  list(by_sample = by_sample, by_feature_summary = by_feature_summary)
}

eda_mean_variance_df <- function(m) {
  data.frame(mean = rowMeans(m, na.rm = TRUE), variance = apply(m, 1, stats::var, na.rm = TRUE))
}

eda_transform_diagnostic <- function(m) {
  vals <- as.numeric(m); vals <- vals[is.finite(vals)]
  can_log <- mean(vals <= 0) < 0.01
  log_vals <- NULL
  if (can_log) {
    v <- vals[vals > 0]
    log_vals <- log2(v)
  }
  raw_sample <- if (length(vals) > EDA_MAX_POOLED_VALUES) sample(vals, EDA_MAX_POOLED_VALUES) else vals
  log_sample <- if (!is.null(log_vals)) {
    if (length(log_vals) > EDA_MAX_POOLED_VALUES) sample(log_vals, EDA_MAX_POOLED_VALUES) else log_vals
  } else NULL
  list(can_log = can_log, raw_sample = raw_sample, log_sample = log_sample,
       skew_raw = eda_skewness(vals), skew_log = if (!is.null(log_vals)) eda_skewness(log_vals) else NA_real_)
}

eda_final_summary <- function(overview, norm_assess, normality, samp_outliers, feat_outliers) {
  n_bad_samples <- sum(samp_outliers$n_flags >= 2)
  pct_nzv <- 100 * (overview$n_constant_features + overview$n_near_zero_var_features) / max(1, overview$n_features)

  quality <- if (overview$n_infinite > 0 || overview$pct_missing > 20 ||
                 overview$n_duplicated_features > 0.05 * overview$n_features) {
    "Needs attention"
  } else if (overview$pct_missing > 5 || n_bad_samples > 0) {
    "Moderate"
  } else "Good"

  outlier_state <- if (n_bad_samples == 0) {
    "No major potential outlier samples detected"
  } else if (n_bad_samples == 1) {
    "One potential outlier sample detected - investigate before removal"
  } else {
    sprintf("%d samples show multiple outlier signals - investigate before removal", n_bad_samples)
  }

  missing_state <- if (overview$pct_missing == 0) "None"
    else if (overview$pct_missing < 1) "Low"
    else if (overview$pct_missing < 10) "Moderate"
    else "High"

  variance_state <- if (pct_nzv > 20) "Many low-variance features"
    else if (pct_nzv > 5) "Some low-variance features present"
    else "Appropriate"

  next_steps <- character(0)
  if (identical(norm_assess$verdict, "not_normalized")) {
    next_steps <- c(next_steps, "Consider normalization (Preprocessing and Batch Correction tabs) before downstream analysis.")
  } else if (identical(norm_assess$verdict, "inconclusive")) {
    next_steps <- c(next_steps, "Review the normalization evidence above; if in doubt, treat this dataset as not yet normalized.")
  }
  if (n_bad_samples > 0) next_steps <- c(next_steps, "Investigate the flagged sample(s) before deciding whether to exclude them.")
  if (missing_state %in% c("Moderate", "High")) next_steps <- c(next_steps, "Review missing-data patterns; consider imputation or feature filtering during preprocessing.")
  if (overview$n_duplicated_features > 0) next_steps <- c(next_steps, "Resolve duplicated feature identifiers before any row-name-keyed downstream step.")
  if (length(next_steps) == 0) next_steps <- "Proceed to downstream preprocessing, normalization, and analysis."

  list(quality = quality, normalization = norm_assess$label, distribution = normality$label,
       outliers = outlier_state, missing = missing_state, variance = variance_state, next_steps = next_steps)
}

eda_value_axis_label <- function(m) {
  v <- m[is.finite(m)]
  if (length(v) == 0) return("Value")
  if (any(v < 0) || max(v) <= 30) "Log-scale value" else "Value"
}

eda_hist_plot <- function(vals, x_label = "Value") {
  bins <- tryCatch({
    b <- grDevices::nclass.FD(vals)
    if (!is.finite(b) || b < 10) grDevices::nclass.Sturges(vals) else b
  }, error = function(e) grDevices::nclass.Sturges(vals))
  bins <- max(20, min(100, bins))
  ggplot(data.frame(value = vals), aes(x = value)) +
    geom_histogram(bins = bins, fill = ARTHOMIX_COLORS$blue, color = "white", linewidth = 0.15, alpha = 0.9) +
    labs(x = x_label, y = "Frequency") + theme_arthomix()
}

eda_density_plot <- function(vals, x_label = "Value") {
  ggplot(data.frame(value = vals), aes(x = value)) +
    geom_density(fill = ARTHOMIX_COLORS$blue, alpha = 0.35, color = ARTHOMIX_COLORS$blue) +
    labs(x = x_label, y = "Density") + theme_arthomix()
}

eda_box_plot <- function(m) {
  qs <- apply(m, 2, stats::quantile, probs = c(0.05, .25, .5, .75, .95), na.rm = TRUE)
  df <- data.frame(sample = colnames(m), ymin = qs[1, ], lower = qs[2, ], middle = qs[3, ], upper = qs[4, ], ymax = qs[5, ])
  ggplot(df, aes(x = reorder(sample, middle), ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax)) +
    geom_boxplot(stat = "identity", width = 0.7, linewidth = 0.15, fill = ARTHOMIX_COLORS$blue, alpha = 0.5) +
    labs(x = NULL, y = "Value") + theme_arthomix() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank())
}

eda_box_plot_interactive <- function(m, outliers_df) {
  qs <- apply(m, 2, stats::quantile, probs = c(0.05, .25, .5, .75, .95), na.rm = TRUE)
  df <- data.frame(sample = colnames(m), ymin = qs[1, ], lower = qs[2, ], middle = qs[3, ], upper = qs[4, ], ymax = qs[5, ],
                    mean_val = colMeans(m, na.rm = TRUE), stringsAsFactors = FALSE)
  df <- merge(df, outliers_df[, c("sample", "n_flags")], by = "sample", all.x = TRUE)
  df$n_flags[is.na(df$n_flags)] <- 0
  df$outlier <- ifelse(df$n_flags >= 2, "Potential outlier", "Sample")
  df$hover_text <- sprintf("Sample: %s<br>Median: %.2f<br>Mean: %.2f<br>Outlier signals: %d", df$sample, df$middle, df$mean_val, df$n_flags)
  p <- ggplot(df, aes(x = reorder(sample, middle), ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax,
                        fill = outlier, text = hover_text)) +
    geom_boxplot(stat = "identity", width = 0.7, linewidth = 0.15) +
    scale_fill_manual(values = c("Sample" = ARTHOMIX_COLORS$blue, "Potential outlier" = ARTHOMIX_COLORS$red)) +
    labs(x = NULL, y = "Value", fill = NULL) + theme_arthomix() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank())
  plotly::ggplotly(p, tooltip = "text") %>% plotly::layout(legend = list(orientation = "h", y = -0.2))
}

eda_violin_plot <- function(m, max_features = EDA_MAX_VIOLIN_FEATURES, max_samples = EDA_MAX_DENSITY_SAMPLES) {
  samples <- colnames(m)
  if (length(samples) > max_samples) samples <- samples[round(seq(1, length(samples), length.out = max_samples))]
  sub <- m[, samples, drop = FALSE]
  if (nrow(sub) > max_features) sub <- sub[sample(nrow(sub), max_features), , drop = FALSE]
  df <- as.data.frame(sub, check.names = FALSE)
  df$feature <- rownames(df)
  long <- tidyr::pivot_longer(df, -feature, names_to = "sample", values_to = "value")
  long <- long[is.finite(long$value), ]
  med_order <- tapply(long$value, long$sample, stats::median, na.rm = TRUE)
  long$sample <- factor(long$sample, levels = names(sort(med_order)))
  ggplot(long, aes(x = sample, y = value)) +
    geom_violin(fill = ARTHOMIX_COLORS$blue, alpha = 0.35, color = ARTHOMIX_COLORS$blue, linewidth = 0.2, scale = "width") +
    labs(x = NULL, y = "Value") + theme_arthomix() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

eda_sample_density_plot <- function(m, max_samples = EDA_MAX_DENSITY_SAMPLES, x_label = "Value") {
  samples <- colnames(m)
  if (length(samples) > max_samples) samples <- samples[round(seq(1, length(samples), length.out = max_samples))]
  dens_list <- lapply(samples, function(s) {
    v <- m[, s]; v <- v[is.finite(v)]
    if (length(v) < 2 || stats::sd(v) == 0) return(NULL)
    d <- stats::density(v)
    data.frame(sample = s, x = d$x, y = d$y)
  })
  df <- do.call(rbind, dens_list)
  validate(need(!is.null(df) && nrow(df) > 0, "Not enough finite, non-constant values per sample to estimate a distribution."))
  ggplot(df, aes(x = x, y = y, group = sample)) +
    geom_line(alpha = 0.3, linewidth = 0.25, color = ARTHOMIX_COLORS$blue) +
    labs(x = x_label, y = "Density") + theme_arthomix() + theme(legend.position = "none")
}

eda_qq_plot <- function(vals) {
  ggplot(data.frame(value = vals), aes(sample = value)) +
    geom_qq(color = ARTHOMIX_COLORS$blue, alpha = 0.4, size = 0.8) +
    geom_qq_line(color = ARTHOMIX_COLORS$red, linewidth = 0.6) +
    labs(x = "Theoretical quantiles", y = "Sample quantiles") + theme_arthomix()
}

eda_meanvar_plot <- function(df, max_points = EDA_MAX_MEANVAR_POINTS) {
  df <- df[is.finite(df$mean) & is.finite(df$variance), ]
  if (nrow(df) > max_points) df <- df[sample(nrow(df), max_points), ]
  ggplot(df, aes(x = mean, y = variance)) +
    geom_point(alpha = 0.25, size = 0.8, color = ARTHOMIX_COLORS$blue) +
    geom_smooth(method = "loess", se = FALSE, color = ARTHOMIX_COLORS$red, linewidth = 0.6, formula = y ~ x) +
    labs(x = "Mean expression (per feature)", y = "Variance (per feature)") + theme_arthomix()
}

eda_missing_bar_plot <- function(by_sample_df) {
  ggplot(by_sample_df, aes(x = reorder(sample, pct_missing), y = pct_missing)) +
    geom_col(fill = ARTHOMIX_COLORS$blue, alpha = 0.85) +
    labs(x = NULL, y = "% missing") + theme_arthomix() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

eda_transform_diag_plot <- function(diag) {
  df1 <- data.frame(value = diag$raw_sample, stage = "Raw scale (diagnostic view)")
  if (is.null(diag$log_sample)) {
    return(ggplot(df1, aes(x = value)) +
             geom_density(fill = ARTHOMIX_COLORS$blue, alpha = 0.35, color = ARTHOMIX_COLORS$blue) +
             labs(x = "Value", y = "Density",
                  title = "Log2 view unavailable - too many non-positive values") + theme_arthomix())
  }
  df2 <- data.frame(value = diag$log_sample, stage = "Log2 (diagnostic only)")
  df <- rbind(df1, df2)
  df$stage <- factor(df$stage, levels = c("Raw scale (diagnostic view)", "Log2 (diagnostic only)"))
  ggplot(df, aes(x = value, fill = stage)) +
    geom_density(alpha = 0.4, color = NA) +
    scale_fill_manual(values = stats::setNames(c(ARTHOMIX_COLORS$blue, ARTHOMIX_COLORS$orange), levels(df$stage))) +
    facet_wrap(~stage, scales = "free", ncol = 2) +
    labs(x = "Value", y = "Density", fill = NULL) + theme_arthomix() + theme(legend.position = "none")
}

eda_corr_heatmap_plot <- function(cor_mat) {
  ord <- eda_hclust(cor_mat)$order
  cm <- cor_mat[ord, ord]
  df <- as.data.frame(as.table(cm))
  colnames(df) <- c("sample_x", "sample_y", "correlation")
  df$sample_x <- factor(df$sample_x, levels = rownames(cm))
  df$sample_y <- factor(df$sample_y, levels = rev(rownames(cm)))
  rng <- range(df$correlation, na.rm = TRUE)
  ggplot(df, aes(x = sample_x, y = sample_y, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = ARTHOMIX_COLORS$red, mid = "white", high = ARTHOMIX_COLORS$blue,
                          midpoint = mean(rng), limits = rng) +
    labs(x = NULL, y = NULL, fill = "r") + theme_arthomix() +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
}

eda_dendro_plot <- function(cor_mat) {
  hc <- eda_hclust(cor_mat)
  graphics::plot(hc, main = "", xlab = "", sub = "", ylab = "Distance (1 - correlation)",
                  cex = if (length(hc$order) > 60) 0.01 else 0.75)
}

eda_pca_plot <- function(pca, outliers_df) {
  df <- as.data.frame(pca$scores[, 1:2, drop = FALSE])
  colnames(df) <- c("PC1", "PC2")
  df$sample <- rownames(df)
  df <- merge(df, outliers_df[, c("sample", "n_flags")], by = "sample", all.x = TRUE)
  df$n_flags[is.na(df$n_flags)] <- 0
  df$outlier <- ifelse(df$n_flags >= 2, "Potential outlier", "Sample")
  df$hover_text <- sprintf("Sample: %s<br>PC1: %.2f<br>PC2: %.2f", df$sample, df$PC1, df$PC2)
  p <- ggplot(df, aes(x = PC1, y = PC2, color = outlier, text = hover_text)) +
    geom_point(size = 2.2, alpha = 0.85) +
    scale_color_manual(values = c("Sample" = ARTHOMIX_COLORS$blue, "Potential outlier" = ARTHOMIX_COLORS$red)) +
    labs(x = sprintf("PC1 (%.1f%% variance)", pca$var_exp[1]), y = sprintf("PC2 (%.1f%% variance)", pca$var_exp[2]), color = NULL) +
    theme_arthomix()
  plotly::ggplotly(p, tooltip = "text") %>% plotly::layout(legend = list(orientation = "h", y = -0.2))
}

eda_scree_plot <- function(pca, max_pcs = 10) {
  n <- min(max_pcs, length(pca$var_exp))
  df <- data.frame(pc = factor(paste0("PC", seq_len(n)), levels = paste0("PC", seq_len(n))), var_exp = pca$var_exp[seq_len(n)])
  ggplot(df, aes(x = pc, y = var_exp)) +
    geom_col(fill = ARTHOMIX_COLORS$blue, alpha = 0.85) +
    labs(x = NULL, y = "Variance explained (%)") + theme_arthomix()
}

eda_status_panel_ui <- function(norm_assess, overview) {
  state_class <- switch(norm_assess$verdict, normalized = "explore-status-good",
                          not_normalized = "explore-status-info", "explore-status-unknown")
  state_icon <- switch(norm_assess$verdict, normalized = "circle-check",
                         not_normalized = "circle-info", "circle-question")
  div(class = paste("card explore-status-panel", state_class),
      div(class = "explore-status-head",
          icon(state_icon, class = "explore-status-icon"),
          div(span(class = "explore-status-eyebrow", "Normalization assessment"),
              div(class = "explore-status-headline", norm_assess$label))),
      tags$ul(lapply(norm_assess$evidence, tags$li)),
      p(class = "explore-status-detail",
        "This assessment is diagnostic/inferential, based on the observed value distribution - it is not proof, since the original processing history of this file is unknown to the application."),
      div(class = "explore-status-facts",
          div(class = "explore-status-fact", span(class = "explore-status-fact-label", "Detected type"),
              span(class = "explore-status-fact-value", switch(norm_assess$detected_type,
                    counts = "Raw counts", already_normalised = "Normalized / transformed", "Linear- or log-scale expression"))),
          div(class = "explore-status-fact", span(class = "explore-status-fact-label", "Near-integer values"),
              span(class = "explore-status-fact-value", sprintf("%.0f%%", norm_assess$frac_integer * 100))),
          div(class = "explore-status-fact", span(class = "explore-status-fact-label", "Negative values present"),
              span(class = "explore-status-fact-value", if (norm_assess$has_negative) "Yes" else "No")),
          div(class = "explore-status-fact", span(class = "explore-status-fact-label", "Between-sample agreement"),
              span(class = "explore-status-fact-value", if (norm_assess$between_sample_differs) "Differences detected" else "Comparable"))
      ))
}

eda_summary_card_ui <- function(summ) {
  fact <- function(label, value) div(class = "explore-status-fact",
    span(class = "explore-status-fact-label", label), span(class = "explore-status-fact-value", value))
  div(class = "card explore-status-panel explore-status-good",
      div(class = "explore-status-head", icon("flag-checkered", class = "explore-status-icon"),
          div(span(class = "explore-status-eyebrow", "EDA summary"), div(class = "explore-status-headline", "Exploratory data analysis complete"))),
      div(class = "explore-status-facts",
          fact("Dataset quality", summ$quality), fact("Normalization", summ$normalization),
          fact("Distribution", summ$distribution), fact("Outliers", summ$outliers),
          fact("Missing data", summ$missing), fact("Variance", summ$variance)),
      tags$hr(),
      div(class = "explore-status-eyebrow", "Recommended next step"),
      tags$ul(lapply(summ$next_steps, tags$li))
  )
}

eda_upload_info_ui <- function(ns, parsed) {
  h <- expr_raw_health(parsed$expr)
  div(class = "card",
      div(class = "card-title", icon("table"), "Step 2 - Review dataset structure"),
      div(class = "explore-stats-row", fluidRow(
        valueBox(format(h$n_features, big.mark = ","), "Rows (features)", icon = icon("dna"), color = "light-blue", width = 3),
        valueBox(format(h$n_samples, big.mark = ","), "Numeric (sample) columns", icon = icon("vial"), color = "green", width = 3),
        valueBox(format(length(parsed$nonnumeric_cols), big.mark = ","), "Non-numeric columns", icon = icon("font"), color = "purple", width = 3),
        valueBox(format(h$n_missing, big.mark = ","), "Missing values", icon = icon("circle-question"), color = if (h$n_missing > 0) "yellow" else "green", width = 3)
      )),
      p(class = "submodule-desc", paste(
        sprintf("Detected \"%s\" as the feature-identifier column.", parsed$id_col_name %||% "(none found - using row position)"),
        if (length(parsed$nonnumeric_cols) > 0) {
          sprintf("%d additional non-numeric column(s) (%s) were set aside and excluded from numeric analysis.",
                  length(parsed$nonnumeric_cols), paste(utils::head(parsed$nonnumeric_cols, 5), collapse = ", "))
        } else "Every remaining column was detected as numeric.",
        if (h$n_infinite > 0) sprintf("%d infinite value(s) were detected and will be treated as missing throughout this tab.", h$n_infinite) else NULL,
        if (h$n_duplicated_features > 0) sprintf("%d duplicated feature identifier(s) were detected.", h$n_duplicated_features) else NULL
      )),
      h5("First rows / columns"),
      DT::dataTableOutput(ns("head_preview_table"))
  )
}

mod_data_exploration_ui <- function(id) {
  ns <- NS(id)
  div(
    div(class = "empty-note", icon("circle-info"),
        "Upload an individual raw molecular dataset here to assess its structure and quality before preprocessing, normalization, or downstream analysis. This tool is fully independent of the app's shared active dataset (set on the Dataset tab) - nothing you do here changes it, and nothing you upload here is written back anywhere."),
    withSpinner(uiOutput(ns("body_ui")), color = "#2563EB", type = 6)
  )
}

mod_data_exploration_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    raw_data <- reactiveVal(NULL)
    raw_version <- reactiveVal(0L)

    observeEvent(input$raw_file, {
      fi <- input$raw_file
      parsed <- tryCatch(eda_parse_upload(fi$datapath, fi$name),
                          error = function(e) list(ok = FALSE, error = "The file could not be processed. Please check its format and try again."))
      raw_data(parsed)
      raw_version(isolate(raw_version()) + 1L)
    }, ignoreInit = TRUE)

    output$body_ui <- renderUI({
      parsed <- raw_data()
      tagList(
        div(class = "card",
            div(class = "card-title", icon("upload"), "Step 1 - Upload raw data"),
            p(class = "submodule-desc",
              "One feature (gene / probe / CpG site / protein) per row, one sample per column, with a feature-identifier column first - the same layout every raw expression/intensity matrix in this app uses. Accepted formats: CSV, TSV, or TXT. Excel files are not currently supported anywhere in this application; please export to CSV first."),
            fileInput(ns("raw_file"), "Raw data file", accept = c(".csv", ".tsv", ".txt"), width = "100%"),
            if (!is.null(parsed) && !isTRUE(parsed$ok)) div(class = "empty-note", icon("triangle-exclamation"), parsed$error)
        ),
        if (!is.null(parsed) && isTRUE(parsed$ok)) tagList(
          eda_upload_info_ui(ns, parsed),
          div(class = "card",
              div(class = "card-title", icon("play"), "Step 3 - Run Exploratory Data Analysis"),
              p(class = "submodule-desc",
                "Runs the full EDA pipeline in one step: descriptive statistics, distribution and normality diagnostics, a normalization-status assessment, outlier detection, PCA and sample structure, sample correlation, missing-data and low-variance feature analysis, the mean-variance relationship, and a final plain-language summary with a recommended next step. Nothing here modifies the file you uploaded."),
              actionButton(ns("run_btn"), "Run Exploratory Data Analysis", icon = icon("flask"), class = "btn-primary")
          ),
          withSpinner(uiOutput(ns("results_ui")), color = "#2563EB", type = 6)
        )
      )
    })

    output$head_preview_table <- DT::renderDataTable({
      req(raw_data(), isTRUE(raw_data()$ok))
      m <- raw_data()$expr
      m <- m[seq_len(min(8, nrow(m))), seq_len(min(10, ncol(m))), drop = FALSE]
      df <- data.frame(feature = rownames(m), round(m, 3), check.names = FALSE, stringsAsFactors = FALSE)
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact")
    })

    eda_result <- eventReactive(input$run_btn, {
      parsed <- raw_data()
      validate(need(isTRUE(parsed$ok), "Upload a valid raw data file first."))
      m <- parsed$expr

      withProgress(message = "Running exploratory data analysis", value = 0, {
        incProgress(0.05, detail = "Dataset overview")
        overview <- eda_overview(parsed)

        incProgress(0.15, detail = "Descriptive statistics")
        feat_stats <- eda_descriptive_stats(m, 1)
        samp_stats <- eda_descriptive_stats(m, 2)

        incProgress(0.1, detail = "Distribution diagnostics")
        vals_finite <- as.numeric(m); vals_finite <- vals_finite[is.finite(vals_finite)]
        validate(need(length(vals_finite) > 0, "No finite numeric values were found to analyze."))
        pooled_sample <- if (length(vals_finite) > EDA_MAX_POOLED_VALUES) sample(vals_finite, EDA_MAX_POOLED_VALUES) else vals_finite
        x_label <- eda_value_axis_label(m)
        normality <- eda_normality_summary(m)

        incProgress(0.1, detail = "Normalization assessment")
        norm_assess <- eda_normalization_assessment(m)

        incProgress(0.15, detail = "PCA and sample structure")
        pca <- tryCatch(eda_pca(m), error = function(e) NULL)
        corr <- tryCatch(eda_sample_correlation(m), error = function(e) NULL)

        incProgress(0.15, detail = "Outlier detection")
        feat_outliers <- eda_feature_outliers(m, feat_stats)
        samp_outliers <- eda_sample_outliers(m, pca)

        incProgress(0.15, detail = "Missing data, variance & transformation")
        missingness <- eda_missingness(m)
        meanvar_df <- eda_mean_variance_df(m)
        transform_diag <- eda_transform_diagnostic(m)

        incProgress(0.15, detail = "Summary")
        summ <- eda_final_summary(overview, norm_assess, normality, samp_outliers, feat_outliers)
      })

      list(m = m, parsed = parsed, overview = overview, feat_stats = feat_stats, samp_stats = samp_stats,
           pooled_sample = pooled_sample, x_label = x_label, normality = normality, norm_assess = norm_assess,
           pca = pca, corr = corr, feat_outliers = feat_outliers, samp_outliers = samp_outliers,
           missingness = missingness, meanvar_df = meanvar_df, transform_diag = transform_diag, summary = summ,
           run_version = raw_version())
    })

    output$overview_ui <- renderUI({
      res <- eda_result(); o <- res$overview
      tagList(
        div(class = "explore-stats-row", fluidRow(
          valueBox(format(o$n_samples, big.mark = ","), "Samples (numeric columns)", icon = icon("vial"), color = "light-blue", width = 3),
          valueBox(format(o$n_features, big.mark = ","), "Features (rows)", icon = icon("dna"), color = "green", width = 3),
          valueBox(sprintf("%.2f%%", o$pct_missing), "Missing values", icon = icon("circle-question"), color = if (o$pct_missing > 0) "yellow" else "green", width = 3),
          valueBox(format(o$n_infinite, big.mark = ","), "Infinite values", icon = icon("triangle-exclamation"), color = if (o$n_infinite > 0) "red" else "green", width = 3)
        ), fluidRow(
          valueBox(format(o$n_duplicated_features, big.mark = ","), "Duplicated feature IDs", icon = icon("clone"), color = if (o$n_duplicated_features > 0) "yellow" else "green", width = 3),
          valueBox(format(o$n_duplicated_samples, big.mark = ","), "Duplicated sample columns", icon = icon("copy"), color = if (o$n_duplicated_samples > 0) "yellow" else "green", width = 3),
          valueBox(format(o$n_constant_features, big.mark = ","), "Constant features", icon = icon("minus"), color = if (o$n_constant_features > 0) "yellow" else "green", width = 3),
          valueBox(format(o$n_near_zero_var_features, big.mark = ","), "Near-zero-variance features", icon = icon("chart-line"), color = "purple", width = 3)
        )),
        DT::dataTableOutput(ns("overview_stats_table"))
      )
    })

    output$overview_stats_table <- DT::renderDataTable({
      o <- eda_result()$overview
      df <- data.frame(statistic = c("Mean", "Median", "SD", "Variance", "Min", "Max", "Q1", "Q3", "IQR"),
                        value = c(o$mean, o$median, o$sd, o$var, o$min, o$max, o$q1, o$q3, o$iqr))
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE), class = "stripe hover compact") %>%
        DT::formatRound(columns = "value", digits = 3)
    })

    output$feat_stats_table <- DT::renderDataTable({
      df <- eda_result()$feat_stats
      DT::datatable(df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatRound(columns = setdiff(colnames(df), c("feature", "n")), digits = 3)
    })
    output$samp_stats_table <- DT::renderDataTable({
      df <- eda_result()$samp_stats
      DT::datatable(df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatRound(columns = setdiff(colnames(df), c("sample", "n")), digits = 3)
    })

    output$hist_plot <- renderPlot(eda_hist_plot(eda_result()$pooled_sample, eda_result()$x_label))
    output$density_plot <- renderPlot(eda_density_plot(eda_result()$pooled_sample, eda_result()$x_label))
    output$box_plot <- renderPlot(eda_box_plot(eda_result()$m))
    output$violin_plot <- renderPlot(eda_violin_plot(eda_result()$m))
    output$sample_density_plot <- renderPlot(eda_sample_density_plot(eda_result()$m, x_label = eda_result()$x_label))

    output$normality_ui <- renderUI({
      n <- eda_result()$normality
      tagList(
        fluidRow(
          column(6, withSpinner(plotOutput(ns("qq_plot"), height = 300), color = "#2563EB", type = 6)),
          column(6,
            div(class = "card explore-diag-card",
                div(class = "card-title", icon("clipboard-check"), "Distribution diagnostics"),
                div(class = "explore-diag-row", span(class = "explore-diag-label", "Skewness"), span(class = "explore-diag-value", sprintf("%.3f", n$skewness))),
                div(class = "explore-diag-row", span(class = "explore-diag-label", "Excess kurtosis"), span(class = "explore-diag-value", sprintf("%.3f", n$kurtosis))),
                div(class = "explore-diag-row", span(class = "explore-diag-label", "Assessment"), span(class = "explore-diag-value", n$label)),
                div(class = "explore-diag-row", span(class = "explore-diag-label", "Shapiro-Wilk W"),
                    span(class = "explore-diag-value", if (is.na(n$shapiro_W)) "Not computed" else sprintf("%.4f", n$shapiro_W))),
                div(class = "explore-diag-row", span(class = "explore-diag-label", "Shapiro-Wilk p"),
                    span(class = "explore-diag-value", if (is.na(n$shapiro_p)) "Not computed" else format.pval(n$shapiro_p, digits = 3))),
                p(class = "empty-note", icon("circle-info"),
                  if (is.na(n$shapiro_p)) {
                    "Not enough finite values to run a normality test."
                  } else {
                    sprintf("Computed on a random sample of %s values (Shapiro-Wilk is capped at 5,000 and becomes extremely sensitive at large n - a significant p-value here is common for real data and is not, by itself, evidence of a meaningful problem). Treat skewness/kurtosis and the Q-Q plot as the primary signal.", format(n$shapiro_n, big.mark = ","))
                  })
            ))
        )
      )
    })
    output$qq_plot <- renderPlot(eda_qq_plot(eda_result()$normality$qq_values))

    output$normalization_ui <- renderUI(eda_status_panel_ui(eda_result()$norm_assess, eda_result()$overview))

    output$outliers_ui <- renderUI({
      res <- eda_result()
      tagList(
        p(class = "submodule-desc",
          "Sample-level flags combine robust (median/MAD-based) signal, detection-rate, and correlation-to-cohort checks with a PCA-distance-from-centroid check. Feature-level flags identify extreme variance, extreme skewness, or excessive missingness. Nothing is removed automatically - flagged rows/columns are for investigation, not automatic exclusion."),
        withSpinner(plotly::plotlyOutput(ns("outlier_box_plot"), height = 340), color = "#2563EB", type = 6),
        h5("Sample-level flags"), DT::dataTableOutput(ns("samp_outlier_table")),
        h5("Feature-level flags (most-flagged first)"), DT::dataTableOutput(ns("feat_outlier_table"))
      )
    })
    output$outlier_box_plot <- plotly::renderPlotly(eda_box_plot_interactive(eda_result()$m, eda_result()$samp_outliers))
    output$samp_outlier_table <- DT::renderDataTable({
      df <- eda_result()$samp_outliers
      df <- df[order(-df$n_flags), ]
      DT::datatable(df, rownames = FALSE, options = list(dom = "tp", scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatRound(columns = intersect(c("signal", "detected", "mean_cor", "pca_distance"), colnames(df)), digits = 3)
    })
    output$feat_outlier_table <- DT::renderDataTable({
      df <- eda_result()$feat_outliers
      df$n_flags <- rowSums(df[, c("flag_extreme_variance", "flag_extreme_skew", "flag_high_missing")])
      df <- df[order(-df$n_flags, -df$pct_missing), ]
      DT::datatable(df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatRound(columns = c("variance", "skewness", "pct_missing"), digits = 3)
    })

    output$pca_ui <- renderUI({
      res <- eda_result()
      if (is.null(res$pca)) {
        return(div(class = "empty-note", icon("circle-info"), "PCA requires at least 3 samples and 3 informative features - not enough usable data to compute it for this dataset."))
      }
      tagList(
        p(class = "submodule-desc", sprintf(
          "Computed on the %s most-variable features (of %s total) after median-imputing missing values for this diagnostic view only - the uploaded data itself is unmodified.",
          format(res$pca$n_features_used, big.mark = ","), format(res$pca$n_features_total, big.mark = ","))),
        fluidRow(
          column(8, withSpinner(plotly::plotlyOutput(ns("pca_plot"), height = 380), color = "#2563EB", type = 6)),
          column(4, withSpinner(plotOutput(ns("scree_plot"), height = 380), color = "#2563EB", type = 6))
        )
      )
    })
    output$pca_plot <- plotly::renderPlotly(eda_pca_plot(eda_result()$pca, eda_result()$samp_outliers))
    output$scree_plot <- renderPlot(eda_scree_plot(eda_result()$pca))

    output$correlation_ui <- renderUI({
      res <- eda_result()
      if (is.null(res$corr)) {
        return(div(class = "empty-note", icon("circle-info"), "Sample correlation requires at least 2 samples - not enough usable data to compute it for this dataset."))
      }
      tagList(
        p(class = "submodule-desc", sprintf(
          "Pairwise sample correlation and hierarchical clustering (average linkage, 1 - correlation distance) on the %s most-variable features (of %s total).",
          format(res$corr$n_features_used, big.mark = ","), format(res$corr$n_features_total, big.mark = ","))),
        fluidRow(
          column(6, withSpinner(plotOutput(ns("corr_heatmap"), height = 360), color = "#2563EB", type = 6)),
          column(6, withSpinner(plotOutput(ns("dendro_plot"), height = 360), color = "#2563EB", type = 6))
        )
      )
    })
    output$corr_heatmap <- renderPlot(eda_corr_heatmap_plot(eda_result()$corr$cor))
    output$dendro_plot <- renderPlot(eda_dendro_plot(eda_result()$corr$cor))

    output$missing_ui <- renderUI({
      res <- eda_result()
      tagList(
        fluidRow(
          column(7, withSpinner(plotOutput(ns("missing_bar_plot"), height = 300), color = "#2563EB", type = 6)),
          column(5, DT::dataTableOutput(ns("missing_feature_table")))
        )
      )
    })
    output$missing_bar_plot <- renderPlot(eda_missing_bar_plot(eda_result()$missingness$by_sample))
    output$missing_feature_table <- DT::renderDataTable({
      DT::datatable(eda_result()$missingness$by_feature_summary, rownames = FALSE,
                    colnames = c("% missing (per feature)", "Number of features"),
                    options = list(dom = "t"), class = "stripe hover compact")
    })

    output$meanvar_plot <- renderPlot(eda_meanvar_plot(eda_result()$meanvar_df))

    output$transform_ui <- renderUI({
      d <- eda_result()$transform_diag
      tagList(
        div(class = "empty-note", icon("circle-info"), "Diagnostic only - the uploaded raw data has not been modified."),
        withSpinner(plotOutput(ns("transform_plot"), height = 280), color = "#2563EB", type = 6),
        p(class = "submodule-desc", sprintf(
          "Pooled skewness: %.2f raw vs. %s log2. %s",
          d$skew_raw, if (is.na(d$skew_log)) "n/a" else sprintf("%.2f", d$skew_log),
          if (!d$can_log) "Log2 view unavailable: too many non-positive values in this matrix (log2 is undefined for values <= 0)."
          else if (is.na(d$skew_log) || abs(d$skew_log) < abs(d$skew_raw)) "The log2 view is less skewed than the raw scale, which is common for count-like or heavy-tailed expression data."
          else "The log2 view is not less skewed than the raw scale here - this dataset may already be on a log or otherwise transformed scale."
        ))
      )
    })
    output$transform_plot <- renderPlot(eda_transform_diag_plot(eda_result()$transform_diag))

    output$summary_ui <- renderUI(eda_summary_card_ui(eda_result()$summary))

    output$download_summary <- downloadHandler(
      filename = function() "data_exploration_summary.csv",
      content = function(file) {
        res <- eda_result(); o <- res$overview; s <- res$summary; na_ass <- res$norm_assess
        df <- data.frame(
          field = c("Samples", "Features", "Missing (%)", "Infinite values", "Duplicated feature IDs",
                      "Duplicated sample columns", "Constant features", "Near-zero-variance features",
                      "Normalization assessment", "Distribution", "Dataset quality", "Outliers", "Missing data", "Variance",
                      "Recommended next step"),
          value = c(o$n_samples, o$n_features, sprintf("%.2f", o$pct_missing), o$n_infinite, o$n_duplicated_features,
                      o$n_duplicated_samples, o$n_constant_features, o$n_near_zero_var_features,
                      s$normalization, s$distribution, s$quality, s$outliers, s$missing, s$variance,
                      paste(s$next_steps, collapse = " | ")),
          stringsAsFactors = FALSE
        )
        write.csv(df, file, row.names = FALSE)
      }
    )

    output$results_ui <- renderUI({
      cur_version <- raw_version()
      res <- tryCatch(eda_result(), error = function(e) NULL)
      if (is.null(res)) return(NULL)
      if (!identical(res$run_version, cur_version)) {
        return(div(class = "empty-note", icon("circle-info"),
                   "A new file has been uploaded. Click \"Run Exploratory Data Analysis\" above to analyze it - the results below are still from the previous file."))
      }
      tagList(
        box(width = 12, title = "A. Dataset overview", status = "primary", solidHeader = FALSE, uiOutput(ns("overview_ui"))),
        box(width = 12, title = "B. Descriptive statistics", status = "primary", solidHeader = FALSE, collapsible = TRUE,
            tabsetPanel(
              tabPanel("Feature-level", br(), DT::dataTableOutput(ns("feat_stats_table"))),
              tabPanel("Sample-level", br(), DT::dataTableOutput(ns("samp_stats_table")))
            )),
        box(width = 12, title = "C. Distribution analysis", status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "Global (pooled, sampled where the matrix is very large) and per-sample views of the value distribution."),
            tabsetPanel(
              tabPanel("Histogram", br(), withSpinner(plotOutput(ns("hist_plot"), height = 340), color = "#2563EB", type = 6)),
              tabPanel("Density", br(), withSpinner(plotOutput(ns("density_plot"), height = 340), color = "#2563EB", type = 6)),
              tabPanel("Per-sample density", br(), withSpinner(plotOutput(ns("sample_density_plot"), height = 340), color = "#2563EB", type = 6)),
              tabPanel("Boxplot (per sample)", br(), withSpinner(plotOutput(ns("box_plot"), height = 340), color = "#2563EB", type = 6)),
              tabPanel("Violin (per sample)", br(), withSpinner(plotOutput(ns("violin_plot"), height = 340), color = "#2563EB", type = 6))
            )),
        box(width = 12, title = "D. Normality / distribution assessment", status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "High-dimensional molecular data is not expected to be normally distributed at the whole-matrix level; this section reports the measured shape rather than asserting normality either way."),
            uiOutput(ns("normality_ui"))),
        box(width = 12, title = "E. Normalization status assessment", status = "primary", solidHeader = FALSE, uiOutput(ns("normalization_ui"))),
        box(width = 12, title = "F. Outlier detection", status = "primary", solidHeader = FALSE, collapsible = TRUE, uiOutput(ns("outliers_ui"))),
        box(width = 12, title = "G. PCA / sample structure", status = "primary", solidHeader = FALSE, uiOutput(ns("pca_ui"))),
        box(width = 12, title = "H. Sample correlation / distance", status = "primary", solidHeader = FALSE, uiOutput(ns("correlation_ui"))),
        box(width = 12, title = "I. Missing data analysis", status = "primary", solidHeader = FALSE, collapsible = TRUE, uiOutput(ns("missing_ui"))),
        box(width = 12, title = "J. Low-variance features", status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", sprintf(
              "%s constant feature(s) (zero variance) and %s additional near-zero-variance feature(s) (bottom 1%% of the non-zero variance distribution) were detected out of %s total. These are candidates to consider filtering downstream - they are not removed here.",
              format(res$overview$n_constant_features, big.mark = ","), format(res$overview$n_near_zero_var_features, big.mark = ","),
              format(res$overview$n_features, big.mark = ",")))),
        box(width = 12, title = "K. Mean-variance relationship", status = "primary", solidHeader = FALSE,
            p(class = "submodule-desc", "For count-like molecular data, variance typically increases with mean expression; a strong such trend is one signal that a variance-stabilizing transformation (e.g. log2 or a count-appropriate normalization) may be worth considering."),
            withSpinner(plotOutput(ns("meanvar_plot"), height = 340), color = "#2563EB", type = 6)),
        box(width = 12, title = "L. Before / after transformation diagnostic", status = "primary", solidHeader = FALSE, uiOutput(ns("transform_ui"))),
        box(width = 12, title = "M. EDA summary", status = "primary", solidHeader = FALSE,
            div(class = "table-toolbar", downloadButton(ns("download_summary"), "Download summary (CSV)", class = "btn-sm")),
            uiOutput(ns("summary_ui")))
      )
    })
  })
}
