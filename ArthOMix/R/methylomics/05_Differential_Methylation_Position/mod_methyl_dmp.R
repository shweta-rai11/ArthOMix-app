## R/methylomics/05_Differential_Methylation_Position/mod_methyl_dmp.R
## Methylomics sub-module: Differential Methylation (DMPs).
##
## Two tabs, each gated on data availability:
##   "SVA" tab -
##     1a. "Default analysis (GSE42861)" - reproduces the published
##         sex-stratified limma + bacon-correction model (plain and
##         SVA-adjusted) from precomputed tables; nothing recomputes here.
##         Shown only when methyl_dataset$preloaded is TRUE.
##     1b. Live SVA-adjusted engine - same statistical method
##         (sva::sva() surrogate-variable estimation, appended as model
##         covariates, then bacon::bacon() bias/inflation correction ahead
##         of Benjamini-Hochberg FDR; see
##         data/preloaded/methylomics/tables/script03_dmp_sva_sexstratified/
##         METHODS_dmp_sva_sexstratified.md for the method this reproduces)
##         run live against whatever beta/M-value matrix is loaded via
##         Upload or a live GEO fetch (i.e. methyl_dataset$beta is set and
##         methyl_dataset$preloaded is not TRUE). Shares its sample/probe
##         subsetting logic (mod_methyl_dmp_prepare_subset()) with the DMP
##         tab below so selection can't silently diverge between the two.
##   "DMP" tab -
##     2. "DMP Analysis" - full configurable live EWAS engine (sex/group/
##        covariate selection, filters, limma model, plots, export) that
##        runs against whatever beta/M-value matrix is currently loaded,
##        upload or preloaded. Shown whenever methyl_dataset$beta is set.
##        Deliberately does not apply SVA/bacon correction - it's the
##        fast, unadjusted counterpart to the SVA tab's live engine above.

mod_methyl_dmp_config <- list(
  id = "dmp", title = "Differential Methylation (DMPs)", icon = "chart-scatter", group = "Data",
  description = "Performs differentially methylated positions analysis."
)

## "SVA" and "DMP" subtabs just wrap the existing default_ui/live_ui
## outputs; switching tabs doesn't force sva_run()/live_result() to
## recompute since Shiny only suspends rendering of the hidden tab.
mod_methyl_dmp_ui <- function(id) {
  ns <- NS(id)
  div(
    class = "tx-menu-wrap",
    tabsetPanel(
      id = ns("dmp_subtabs"), type = "tabs",
      tabPanel("SVA", br(), withSpinner(uiOutput(ns("default_ui")), color = "#2563EB", type = 6)),
      tabPanel("DMP", br(), withSpinner(uiOutput(ns("live_ui")), color = "#2563EB", type = 6))
    )
  )
}

## ---- Shared: filter + summarize + volcano, parameterized by which table --

mod_methyl_dmp_filter <- function(df, fdr_col, effect_col, fdr_max, effect_min, direction) {
  keep <- !is.na(df[[fdr_col]]) & df[[fdr_col]] <= fdr_max & abs(df[[effect_col]]) >= effect_min
  if (identical(direction, "hyper")) keep <- keep & df[[effect_col]] > 0
  if (identical(direction, "hypo"))  keep <- keep & df[[effect_col]] < 0
  df[keep, , drop = FALSE]
}

mod_methyl_dmp_volcano <- function(df, effect_col, p_col, effect_label, fdr_max, effect_min) {
  df$neglog10p <- -log10(pmax(df[[p_col]], .Machine$double.xmin))
  df$sig <- !is.na(df[[p_col]]) & df[[p_col]] <= fdr_max & abs(df[[effect_col]]) >= effect_min
  gg <- ggplot(df, aes(x = .data[[effect_col]], y = neglog10p, color = sig)) +
    geom_point(alpha = 0.4, size = 0.8) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = ARTHOMIX_STATUS$critical)) +
    labs(x = effect_label, y = expression(-log[10](p)), color = "Significant") +
    theme_arthomix()
  gg
}

## Genomic inflation factor (lambda_gc), standard median-chi-square formula.
## Must be fed the full raw per-CpG p-value vector (before FDR/effect-size
## filtering). Shared with mod_methyl_dmr.R's live engine.
mod_methyl_lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (length(p) == 0) return(NA_real_)
  stats::median(stats::qchisq(1 - p, df = 1), na.rm = TRUE) / stats::qchisq(0.5, df = 1)
}

## Observed-vs-expected -log10(p) QQ plot; companion to lambda_gc.
mod_methyl_qq_plot <- function(p) {
  p <- sort(p[is.finite(p) & p > 0 & p <= 1])
  n <- length(p)
  if (n == 0) return(NULL)
  df <- data.frame(observed = -log10(p), expected = -log10(stats::ppoints(n)))
  ggplot(df, aes(x = expected, y = observed)) +
    geom_point(color = ARTHOMIX_COLORS$blue, alpha = 0.4, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = expression(-log[10] ~ "(expected p)"), y = expression(-log[10] ~ "(observed p)")) +
    theme_arthomix()
}

## ---- Live-engine helpers ---------------------------------------------------

## Row-chunked limma::lmFit(), reused by mod_methyl_dmr.R's live engine too.
## A full genome-wide fit (400k+ probes) can exceed R's vector memory limit,
## so rows are fit in chunks and the per-row fit summaries (coefficients,
## stdev.unscaled, sigma, df.residual, Amean) are concatenated back into one
## MArrayLM object. Verified bit-for-bit identical topTable() output vs. a
## whole-matrix fit on synthetic data.
methyl_chunked_lmfit <- function(m, design, chunk_size = 20000) {
  n <- nrow(m)
  if (n <= chunk_size) return(limma::lmFit(m, design))
  starts <- seq(1, n, by = chunk_size)
  fits <- lapply(starts, function(s) {
    e <- min(s + chunk_size - 1, n)
    limma::lmFit(m[s:e, , drop = FALSE], design)
  })
  per_row_1d <- c("df.residual", "sigma", "Amean")
  per_row_2d <- c("coefficients", "stdev.unscaled")
  out <- fits[[1]]
  for (nm in per_row_1d) if (!is.null(out[[nm]])) out[[nm]] <- do.call(c, lapply(fits, `[[`, nm))
  for (nm in per_row_2d) if (!is.null(out[[nm]])) out[[nm]] <- do.call(rbind, lapply(fits, `[[`, nm))
  out
}

## Sample + probe subsetting shared by the DMP tab's plain live engine and
## the SVA tab's live SVA-adjusted engine, so the two can't silently diverge
## in which samples/probes they analyze - only the model fit itself differs
## downstream. Mirrors the validation/subsetting steps previously inlined
## in the DMP tab's own live_result() one-for-one (sex subset -> group/ref/
## comp subset -> covariate complete-cases -> missingness/variance/SNP probe
## filters -> beta/M-value matrix construction).
mod_methyl_dmp_prepare_subset <- function(methyl_dataset, sex_choice, sex_col_name,
                                            group_col, ref, comp, covariate_cols,
                                            min_valid_pct, min_variance, snp_filter, anno) {
  sheet <- methyl_dataset$sample_sheet
  validate(need(!is.null(methyl_dataset$beta), "Load a dataset first."))
  validate(need(!is.null(sheet), "No sample sheet loaded."))
  validate(need(!is.null(group_col) && group_col %in% colnames(sheet), "Pick a group column."))
  validate(need(!is.null(ref) && !is.null(comp), "Pick a reference and a comparison group."))
  validate(need(ref != comp, "Reference and comparison groups must be different."))

  sample_ids <- methyl_sheet_sample_ids(sheet, colnames(methyl_dataset$beta))
  common <- intersect(colnames(methyl_dataset$beta), sample_ids)
  validate(need(length(common) >= 6, "Fewer than 6 samples match between the matrix and the sample sheet."))
  beta0 <- methyl_dataset$beta[, common, drop = FALSE]
  ## Preloaded sample sheet is a data.table - coerce to data.frame so
  ## single-bracket column selection by name behaves consistently below.
  ph0 <- as.data.frame(sheet)[match(common, sample_ids), , drop = FALSE]

  ## ---- sex subset --------------------------------------------------
  sex_label <- "All samples"
  if (!identical(sex_choice, "__all__")) {
    validate(need(!is.null(sex_col_name), "No sex column available to subset on."))
    keep_sex <- !is.na(ph0[[sex_col_name]]) & as.character(ph0[[sex_col_name]]) == sex_choice
    validate(need(sum(keep_sex) >= 6, sprintf("Fewer than 6 samples remain after restricting to sex = \"%s\".", sex_choice)))
    beta0 <- beta0[, keep_sex, drop = FALSE]
    ph0 <- ph0[keep_sex, , drop = FALSE]
    sex_label <- sex_choice
  }
  n_after_sex <- ncol(beta0)

  ## ---- group/reference/comparison subset ----------------------------
  grp_raw <- as.character(ph0[[group_col]])
  keep_grp <- !is.na(grp_raw) & grp_raw %in% c(ref, comp)
  beta1 <- beta0[, keep_grp, drop = FALSE]
  ph1 <- ph0[keep_grp, , drop = FALSE]
  grp <- factor(grp_raw[keep_grp], levels = c(ref, comp))
  rm(beta0)  ## free the full-size matrix now that beta1 exists (memory pressure on full arrays)

  ## ---- optional covariates (complete cases only) --------------------
  cov_df <- NULL
  if (length(covariate_cols) > 0) {
    cc <- ph1[, covariate_cols, drop = FALSE]
    complete <- stats::complete.cases(cc)
    validate(need(sum(complete) >= 6, "Fewer than 6 samples have no missing values in the selected covariates. Deselect a covariate, or pick different ones."))
    beta1 <- beta1[, complete, drop = FALSE]
    ph1 <- ph1[complete, , drop = FALSE]
    grp <- grp[complete]
    cov_df <- ph1[, covariate_cols, drop = FALSE]
    for (cl in covariate_cols) if (!is.numeric(cov_df[[cl]])) cov_df[[cl]] <- factor(as.character(cov_df[[cl]]))
  }

  n_ref <- sum(grp == ref, na.rm = TRUE)
  n_comp <- sum(grp == comp, na.rm = TRUE)
  validate(need(n_ref >= 3 && n_comp >= 3,
    sprintf("Each group needs at least 3 samples to fit a model (reference \"%s\": %d, comparison \"%s\": %d).",
            ref, n_ref, comp, n_comp)))

  ## ---- probe filters (missingness / variance / SNP), before the M-value transform ----
  is_m_scale <- identical(methyl_dataset$input_scale, "m")
  beta_scale_full <- if (is_m_scale) 2^beta1 / (1 + 2^beta1) else beta1

  max_na_frac <- 1 - (min_valid_pct %||% 80) / 100
  f_miss <- methyl_filter_missing(beta_scale_full, max_na_frac)
  keep_probe <- f_miss$keep
  f_var <- methyl_filter_variance(beta_scale_full, min_variance %||% 0)
  keep_probe <- keep_probe & f_var$keep
  snp_note <- NULL
  if (isTRUE(snp_filter)) {
    if (isTRUE(anno$ok)) {
      f_snp <- methyl_filter_snp(beta_scale_full, anno)
      keep_probe <- keep_probe & f_snp$keep
      snp_note <- f_snp$note
    } else snp_note <- anno$reason
  }
  validate(need(sum(keep_probe) >= 10,
    "Fewer than 10 CpGs remain after the missingness/variance/SNP filters. Relax the filters and try again."))

  beta1 <- beta1[keep_probe, , drop = FALSE]
  if (is_m_scale) rm(beta_scale_full)
  gc(FALSE)

  ## ---- beta/M-value matrices (filtered subset only) ---------------------
  m <- if (is_m_scale) beta1 else log2(pmin(pmax(beta1, 1e-6), 1 - 1e-6) / (1 - pmin(pmax(beta1, 1e-6), 1 - 1e-6)))
  beta_scale <- if (is_m_scale) 2^m / (1 + 2^m) else beta1

  list(
    m = m, beta_scale = beta_scale, grp = grp, cov_df = cov_df,
    sex_label = sex_label, n_after_sex = n_after_sex, n_ref = n_ref, n_comp = n_comp,
    n_probes_before_filter = nrow(beta1),
    missing_note = f_miss$note, variance_note = f_var$note, snp_note = snp_note
  )
}

## Live SVA + bacon-corrected model, mirroring the precomputed GSE42861
## pipeline (script03_dmp_sva_sexstratified/METHODS_dmp_sva_sexstratified.md)
## so the same statistical method is available live for Upload/GEO-fetched
## data, not just the preloaded reference cohort.
##
## Surrogate variables are estimated with sva::sva(), contrasting a full
## model (~ group [+ covariates]) against the same model with `group`
## removed - exactly the full-vs-null contrast the precomputed pipeline
## uses, so that variation attributable to the group comparison itself
## isn't absorbed into the surrogate variables. Estimation is restricted to
## the `sv_topn` most variable probes (same top-variable-probe convention
## as the precomputed pipeline) since sva::sva() is expensive at
## genome-wide scale; the final limma fit still uses the full filtered
## probe set, with the estimated SVs appended as ordinary covariates.
mod_methyl_sva_fit <- function(m, grp, cov_df, sv_topn = 20000) {
  cov_names <- if (!is.null(cov_df)) colnames(cov_df) else character(0)
  safe_cov <- if (length(cov_names) > 0) sprintf("`%s`", cov_names) else character(0)

  mod_df <- data.frame(grp = grp)
  if (!is.null(cov_df)) mod_df <- cbind(mod_df, cov_df)
  full_form <- stats::as.formula(paste("~ grp", if (length(safe_cov) > 0) paste("+", paste(safe_cov, collapse = " + ")) else ""))
  null_form <- if (length(safe_cov) > 0) stats::as.formula(paste("~", paste(safe_cov, collapse = " + "))) else ~1
  mod1 <- stats::model.matrix(full_form, data = mod_df)
  mod0 <- stats::model.matrix(null_form, data = mod_df)

  n <- nrow(m)
  top_idx <- if (n > sv_topn) order(matrixStats::rowVars(m), decreasing = TRUE)[seq_len(sv_topn)] else seq_len(n)

  sv_obj <- tryCatch(
    sva::sva(dat = as.matrix(m[top_idx, , drop = FALSE]), mod = mod1, mod0 = mod0, method = "irw"),
    error = function(e) NULL
  )

  sv_mat <- NULL
  if (!is.null(sv_obj) && !is.null(sv_obj$n.sv) && sv_obj$n.sv > 0 && !is.null(sv_obj$sv) && NCOL(sv_obj$sv) > 0) {
    sv_mat <- as.matrix(sv_obj$sv)[, seq_len(sv_obj$n.sv), drop = FALSE]
    colnames(sv_mat) <- paste0("SV", seq_len(ncol(sv_mat)))
  }

  design_grp <- stats::model.matrix(~0 + grp)
  colnames(design_grp) <- levels(grp)
  design_cov <- NULL
  if (length(cov_names) > 0) {
    dc <- stats::model.matrix(stats::as.formula(paste("~", paste(safe_cov, collapse = " + "))), data = cov_df)
    design_cov <- dc[, setdiff(colnames(dc), "(Intercept)"), drop = FALSE]
  }

  design <- design_grp
  if (!is.null(design_cov)) design <- cbind(design, design_cov)
  if (!is.null(sv_mat)) design <- cbind(design, sv_mat)

  ## SVA's own component count is a heuristic (Buja-Eyuboglu "be" method) and
  ## can occasionally overshoot on a small live cohort; drop trailing SVs
  ## (the last-added columns) rather than failing the model outright if that
  ## makes the design matrix rank-deficient.
  n_sv <- if (!is.null(sv_mat)) ncol(sv_mat) else 0L
  while (n_sv > 0 && qr(design)$rank < ncol(design)) {
    n_sv <- n_sv - 1L
    design <- design[, seq_len(ncol(design) - 1), drop = FALSE]
  }

  list(design = design, n_sv = n_sv)
}

## Auto-detects a sex/gender column by name; same candidates as
## mod_methyl_qc.R's sex-check panel.
mod_methyl_dmp_sex_col <- function(sheet) {
  if (is.null(sheet)) return(NULL)
  cols <- intersect(c("sex", "Sex", "SEX", "gender", "Gender"), colnames(sheet))
  if (length(cols) == 0) return(NULL)
  cols[1]
}

## Builds "All samples / Female / Male" radio choices from the sex column's
## actual values - maps to Female/Male labels only when exactly one level
## matches each pattern (F/Female, M/Male), else falls back to raw values.
mod_methyl_dmp_sex_choices <- function(sheet, sex_col) {
  base <- c("All samples" = "__all__")
  if (is.null(sex_col) || !sex_col %in% colnames(sheet)) return(base)
  lvls <- unique(as.character(stats::na.omit(sheet[[sex_col]])))
  if (length(lvls) == 0) return(base)
  fem <- lvls[grepl("^f(emale)?$", lvls, ignore.case = TRUE)]
  mal <- lvls[grepl("^m(ale)?$", lvls, ignore.case = TRUE)]
  if (length(fem) == 1 && length(mal) == 1) {
    return(c(base, stats::setNames(fem, "Female"), stats::setNames(mal, "Male")))
  }
  c(base, stats::setNames(lvls, lvls))
}

## Candidate covariate columns: >=2 distinct non-missing values, excluding
## free-text/ID-like character columns (every value unique); numeric columns
## are kept regardless since continuous covariates can be all-unique too.
mod_methyl_dmp_covariate_cols <- function(sheet, exclude) {
  cols <- setdiff(colnames(sheet), exclude)
  keep <- vapply(cols, function(cl) {
    v <- sheet[[cl]]
    nu <- length(unique(stats::na.omit(v)))
    if (nu < 2) return(FALSE)
    if (!is.numeric(v) && nu == nrow(sheet)) return(FALSE)
    TRUE
  }, logical(1))
  cols[keep]
}

mod_methyl_dmp_manhattan <- function(df, fdr_max = 0.05) {
  d <- df[!is.na(df$chr) & !is.na(df$pos) & !is.na(df$fdr), , drop = FALSE]
  d$chr_num <- suppressWarnings(as.integer(gsub("^chr", "", d$chr, ignore.case = TRUE)))
  d <- d[!is.na(d$chr_num), , drop = FALSE]
  validate(need(nrow(d) > 0, "No tested CpG has both a chromosome and a genomic position from the available annotation - a Manhattan plot needs both."))
  d <- d[order(d$chr_num, d$pos), , drop = FALSE]
  d$idx <- seq_len(nrow(d))
  d$sig <- d$fdr <= fdr_max
  axis_df <- stats::aggregate(idx ~ chr_num, d, mean)
  ggplot(d, aes(x = idx, y = -log10(pmax(fdr, .Machine$double.xmin)))) +
    geom_point(aes(color = factor(chr_num %% 2)), size = 0.6, alpha = 0.5, show.legend = FALSE) +
    { if (any(d$sig)) geom_point(data = d[d$sig, , drop = FALSE], color = ARTHOMIX_STATUS$critical, size = 1.1) } +
    scale_color_manual(values = c("0" = ARTHOMIX_COLORS$blue, "1" = ARTHOMIX_COLORS$ink_muted)) +
    scale_x_continuous(breaks = axis_df$idx, labels = axis_df$chr_num, expand = expansion(mult = 0.01)) +
    geom_hline(yintercept = -log10(fdr_max), linetype = "dashed", color = "grey50") +
    labs(x = "Chromosome", y = expression(-log[10](FDR))) +
    theme_arthomix() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))
}

mod_methyl_dmp_topplot <- function(df, rank_by = c("fdr", "dbeta", "combined"), n = 20) {
  rank_by <- match.arg(rank_by)
  d <- df[!is.na(df$fdr) & !is.na(df$dbeta), , drop = FALSE]
  validate(need(nrow(d) > 0, "No tested CpG has both an FDR and a Δβ value to rank."))
  ord <- switch(rank_by,
    fdr      = order(d$fdr, -abs(d$dbeta)),
    dbeta    = order(-abs(d$dbeta), d$fdr),
    combined = order(d$fdr / pmax(abs(d$dbeta), 1e-6))
  )
  d <- utils::head(d[ord, , drop = FALSE], n)
  label <- ifelse(!is.na(d$gene) & nzchar(d$gene), sprintf("%s (%s)", d$cpg, d$gene), d$cpg)
  d$label <- factor(label, levels = rev(label))
  gg <- ggplot(d, aes(x = dbeta, y = label, fill = dbeta > 0)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = ARTHOMIX_STATUS$critical, "FALSE" = ARTHOMIX_COLORS$blue), guide = "none") +
    labs(x = "Δβ", y = NULL) +
    theme_arthomix()
  list(plot = gg, cpgs = as.character(d$cpg))
}

mod_methyl_dmp_betadist <- function(beta_mat, cpgs, grp) {
  cpgs <- intersect(cpgs, rownames(beta_mat))
  validate(need(length(cpgs) > 0, "None of the top-ranked CpGs are present in the analyzed beta matrix."))
  m <- beta_mat[cpgs, , drop = FALSE]
  long <- data.frame(
    cpg = rep(rownames(m), times = ncol(m)),
    beta = as.vector(m),
    group = rep(as.character(grp), each = nrow(m))
  )
  long$cpg <- factor(long$cpg, levels = cpgs)
  ggplot(long, aes(x = cpg, y = beta, fill = group)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.75)) +
    scale_fill_manual(values = arthomix_pair(levels(grp))) +
    labs(x = NULL, y = "β value", fill = NULL) +
    theme_arthomix() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

## UI for the SVA tab's live engine (Upload/GEO-fetched data). Deliberately
## mirrors the DMP tab's own live_ui card/control layout (same fluidRow(4,8)
## split, same control types) so the two tabs read as one consistent
## module, differing only in which inputs/outputs it targets (svalive_*
## instead of live_*) and that it also runs SVA + bacon correction.
mod_methyl_svalive_panel_ui <- function(ns, methyl_dataset, sc, anno) {
  sheet <- methyl_dataset$sample_sheet
  cols <- colnames(sheet)
  arr_ok <- isTRUE(anno$ok)
  tagList(
    div(class = "card",
        div(class = "card-title", icon("flask"), "SVA-adjusted Analysis (live)"),
        p(class = "empty-note", icon("circle-info"),
          sprintf("Dataset: %s. %s probes x %s samples. Estimates surrogate variables (sva::sva, full-vs-null contrast) and applies bacon bias/inflation correction on top of a limma model - the same method used for the preloaded cohort's reproduced analysis above, run live here against your own data.",
                  methyl_dataset$source %||% "(unnamed)",
                  format(nrow(methyl_dataset$beta), big.mark = ","), ncol(methyl_dataset$beta))),
        fluidRow(
          column(4,
            tags$h5("Sex"),
            radioButtons(ns("svalive_sex"), NULL, inline = TRUE, choices = mod_methyl_dmp_sex_choices(sheet, sc), selected = "__all__"),
            if (length(mod_methyl_dmp_sex_choices(sheet, sc)) <= 1) p(class = "empty-note", icon("circle-info"), "No usable sex information was found for this dataset - showing pooled analysis only."),

            tags$h5("Comparison"),
            selectInput(ns("svalive_group_col"), "Group column", choices = cols,
                        selected = intersect(c("group", "Group", "disease", "Disease"), cols)[1] %||% cols[1]),
            uiOutput(ns("svalive_level_ui")),

            tags$h5("Filters"),
            numericInput(ns("svalive_fdr"), "Adjusted p-value (FDR) threshold", value = 0.05, min = 0, max = 1, step = 0.01),
            numericInput(ns("svalive_dbeta"), "Absolute Δβ threshold", value = 0, min = 0, max = 1, step = 0.01),
            radioButtons(ns("svalive_direction"), "Direction", inline = TRUE,
                         choices = c("All DMPs" = "any", "Hypermethylated" = "hyper", "Hypomethylated" = "hypo"), selected = "any"),
            numericInput(ns("svalive_min_valid_pct"), "Minimum valid (non-missing) sample %", value = 80, min = 0, max = 100, step = 5),
            numericInput(ns("svalive_min_variance"), "Minimum methylation variance (optional)", value = 0, min = 0, step = 0.001),
            if (arr_ok) checkboxInput(ns("svalive_snp_filter"), "Remove SNP-associated probes (manifest Probe_rs/CpG_rs/SBE_rs)", value = FALSE)
            else p(class = "empty-note", icon("circle-info"), anno$reason),

            tags$h5("Covariates (optional)"),
            uiOutput(ns("svalive_covariate_ui")),

            actionButton(ns("svalive_run_btn"), "Run SVA-adjusted Analysis", icon = icon("play"), class = "btn-primary")
          ),
          column(8, withSpinner(uiOutput(ns("svalive_results_ui")), color = "#2563EB", type = 6))
        )
    )
  )
}

mod_methyl_dmp_server <- function(id, methyl_dataset, methyl_results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## ================= 1. Default analysis (GSE42861) =====================

    default_data <- reactive({
      if (!METH_DATA_AVAILABLE) return(NULL)
      list(
        pheno = load_default_meth_pheno(),
        plain_f = load_default_dmp("plain", "female"), plain_m = load_default_dmp("plain", "male"),
        sva_f   = load_default_dmp("sva", "female"),   sva_m   = load_default_dmp("sva", "male")
      )
    })

    output$default_ui <- renderUI({
      ## Three data sources this tab has to cover: preloaded (reproduced,
      ## precomputed - below), Upload/GEO with no data yet, and Upload/GEO
      ## with data loaded (live SVA engine, further down this file).
      if (!isTRUE(methyl_dataset$preloaded)) {
        if (is.null(methyl_dataset$beta)) {
          return(div(class = "card",
            div(class = "card-title", icon("upload"), "SVA-adjusted Analysis"),
            p(class = "submodule-desc",
              "Upload a beta/M-value matrix or fetch a dataset from GEO on the Methylomics Dataset tab to run a live surrogate-variable-adjusted, bacon-corrected differential methylation model - the same statistical method used for the preloaded reference cohort's reproduced analysis.")
          ))
        }
        if (is.null(methyl_dataset$sample_sheet)) {
          return(div(class = "card",
            div(class = "card-title", icon("triangle-exclamation"), "No sample sheet"),
            p(class = "submodule-desc", "An SVA-adjusted model needs a group variable - re-load on the Dataset tab with a sample sheet/phenotype file included.")
          ))
        }
        return(mod_methyl_svalive_panel_ui(ns, methyl_dataset, sex_col(), anno_result()))
      }
      d <- default_data()
      req(d)
      n_ra <- sum(d$pheno$group == "RA", na.rm = TRUE); n_ctrl <- sum(d$pheno$group == "Control", na.rm = TRUE)
      n_f <- sum(d$pheno$sex == "F", na.rm = TRUE); n_m <- sum(d$pheno$sex == "M", na.rm = TRUE)
      n_sig_f <- sum(!is.na(d$sva_f$fdr_bacon) & d$sva_f$fdr_bacon < 0.05)
      n_sig_m <- sum(!is.na(d$sva_m$fdr_bacon) & d$sva_m$fdr_bacon < 0.05)
      n_sig_plain_f <- sum(!is.na(d$plain_f$fdr_bacon) & d$plain_f$fdr_bacon < 0.05)
      n_sig_plain_m <- sum(!is.na(d$plain_m$fdr_bacon) & d$plain_m$fdr_bacon < 0.05)
      tagList(
        div(class = "empty-note", icon("flask"),
            sprintf("Default analysis: preloaded whole-blood dataset, %d samples (%d RA / %d Control; %d female / %d male). Sex-stratified limma on M-values (~group + age + smoking + cell-type estimates), bacon-corrected, Benjamini-Hochberg FDR.",
                    nrow(d$pheno), n_ra, n_ctrl, n_f, n_m)),
        div(class = "card",
            div(class = "card-title", icon("check"), "SVA-adjusted model: the usable panel"),
            fluidRow(
              column(4,
                radioButtons(ns("sva_sex"), "Stratum", inline = TRUE, choices = c("Female" = "female", "Male" = "male"), selected = "female"),
                numericInput(ns("sva_fdr"), "FDR threshold", value = 0.05, min = 0, max = 1, step = 0.01),
                numericInput(ns("sva_dbeta"), "Min |Δβ|", value = 0, min = 0, max = 1, step = 0.01),
                radioButtons(ns("sva_direction"), "Direction", inline = TRUE,
                             choices = c("Any" = "any", "Hypermethylated" = "hyper", "Hypomethylated" = "hypo"), selected = "any"),
                actionButton(ns("sva_run_btn"), "View results", icon = icon("play"), class = "btn-primary btn-sm")
              ),
              column(8, withSpinner(uiOutput(ns("sva_volcano_ui")), color = "#2563EB", type = 6))
            ),
            uiOutput(ns("sva_valueboxes_ui")),
            uiOutput(ns("sva_table_ui"))
        ),
        div(class = "card",
            div(class = "card-title", icon("circle-info"), "Why an SVA-adjusted stage exists"),
            p(class = "submodule-desc",
              sprintf("The plain (unadjusted) sex-stratified model found %d female-stratum and %d male-stratum genome-wide-significant CpGs - both strata showed residual genomic inflation after bacon-correction, a real documented calibration problem, not a bug or a null biological result. The SVA-adjusted model above resolves it and is the panel actually used downstream (DMR calling, biomarker panels).",
                      n_sig_plain_f, n_sig_plain_m))
        )
      )
    })

    ## eventReactive on the button click - controls don't update results until "View results" is clicked again.
    sva_run <- eventReactive(input$sva_run_btn, {
      d <- default_data(); req(d)
      df <- if (identical(input$sva_sex, "male")) d$sva_m else d$sva_f
      req(df)
      list(
        df = df, sex = input$sva_sex,
        fdr = input$sva_fdr %||% 0.05, dbeta = input$sva_dbeta %||% 0,
        direction = if (identical(input$sva_direction, "any")) NULL else input$sva_direction
      )
    }, ignoreInit = TRUE)

    sva_has_run <- reactiveVal(FALSE)
    observeEvent(input$sva_run_btn, sva_has_run(TRUE), ignoreInit = TRUE)

    default_filtered <- reactive({
      r <- sva_run()
      mod_methyl_dmp_filter(r$df, "fdr_bacon", "dbeta", r$fdr, r$dbeta, r$direction)
    })

    output$sva_volcano_ui <- renderUI({
      if (!sva_has_run()) {
        return(p(class = "empty-note", icon("circle-info"), "Pick a stratum and filters on the left, then click \"View results\"."))
      }
      withSpinner(plotOutput(ns("default_volcano"), height = 320), color = "#2563EB", type = 6)
    })

    output$default_volcano <- renderPlot({
      r <- sva_run()
      mod_methyl_dmp_volcano(r$df, "dbeta", "fdr_bacon", "Δβ (RA - Control)", r$fdr, r$dbeta)
    })

    output$sva_valueboxes_ui <- renderUI({
      req(sva_has_run())
      r <- sva_run()
      n_sig <- sum(!is.na(r$df$fdr_bacon) & r$df$fdr_bacon < 0.05)
      sex_label <- if (identical(r$sex, "male")) "Male" else "Female"
      div(class = "methyl-stats-row", fluidRow(
        valueBox(format(n_sig), sprintf("%s FDR<0.05", sex_label),
                 icon = icon(if (identical(r$sex, "male")) "mars" else "venus"),
                 color = if (n_sig > 0) "green" else "light-blue", width = 4),
        valueBox(format(nrow(r$df), big.mark = ","), sprintf("CpGs tested (%s)", tolower(sex_label)), icon = icon("dna"), color = "purple", width = 4),
        valueBox(format(nrow(default_filtered()), big.mark = ","), "Passing current filters", icon = icon("filter"), color = "light-blue", width = 4)
      ))
    })

    output$sva_table_ui <- renderUI({
      req(sva_has_run())
      tagList(
        downloadButton(ns("download_default"), "Download filtered CSV", class = "btn-default btn-sm", style = "margin: 8px 0;"),
        DT::dataTableOutput(ns("default_table"))
      )
    })

    output$default_table <- DT::renderDataTable({
      DT::datatable(default_filtered(), rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatSignif(columns = c("logFC_M", "dbeta", "t", "p_raw", "p_bacon", "fdr_bacon"), digits = 4)
    })
    ## Nested two levels deep in dynamic renderUI; Shiny's suspendWhenHidden
    ## visibility tracking misdetects it, so the table never reaches the
    ## browser unless forced to always compute/send.
    outputOptions(output, "default_table", suspendWhenHidden = FALSE)

    output$download_default <- downloadHandler(
      filename = function() { r <- sva_run(); sprintf("gse42861_dmp_sva_%s_fdr%s.csv", r$sex, r$fdr) },
      content = function(file) utils::write.csv(default_filtered(), file, row.names = FALSE)
    )

    ## ================= 1b. Live SVA-adjusted engine (Upload / GEO) ========
    ## Same "SVA" tab, shown instead of the block above whenever the active
    ## dataset isn't the preloaded reference cohort. Reuses
    ## mod_methyl_dmp_prepare_subset() (shared with the DMP tab's plain live
    ## engine, section 2 below) for sample/probe selection, then estimates
    ## surrogate variables and applies bacon correction - see
    ## mod_methyl_sva_fit() above for the method itself.

    output$svalive_level_ui <- renderUI({
      req(input$svalive_group_col)
      sheet <- methyl_dataset$sample_sheet
      req(input$svalive_group_col %in% colnames(sheet))
      levels_available <- sort(unique(as.character(stats::na.omit(sheet[[input$svalive_group_col]]))))
      validate(need(length(levels_available) >= 2, "This column has fewer than two distinct values - pick a different group column."))
      tagList(
        selectInput(ns("svalive_ref"), "Reference group", choices = levels_available, selected = levels_available[1]),
        selectInput(ns("svalive_comp"), "Comparison group", choices = levels_available, selected = levels_available[min(2, length(levels_available))])
      )
    })

    output$svalive_covariate_ui <- renderUI({
      sheet <- methyl_dataset$sample_sheet
      req(sheet, input$svalive_group_col)
      exclude <- c(id_cols(), input$svalive_group_col)
      if (!identical(input$svalive_sex %||% "__all__", "__all__") && !is.null(sex_col())) exclude <- c(exclude, sex_col())
      cand <- mod_methyl_dmp_covariate_cols(sheet, exclude)
      if (length(cand) == 0) {
        return(p(class = "empty-note", icon("circle-info"), "No additional phenotype columns are available to use as covariates."))
      }
      checkboxGroupInput(ns("svalive_covariates"), NULL, choices = cand, selected = character(0))
    })

    svalive_has_run <- reactiveVal(FALSE)
    observeEvent(methyl_dataset$beta, svalive_has_run(FALSE), ignoreNULL = TRUE)

    svalive_result <- eventReactive(input$svalive_run_btn, withProgress(
      message = "Estimating surrogate variables and fitting the bacon-corrected model - slower than the plain DMP model, and can take several minutes on a full genome-wide array...",
      value = 0.15, {
        sub <- mod_methyl_dmp_prepare_subset(
          methyl_dataset, input$svalive_sex %||% "__all__", sex_col(),
          input$svalive_group_col, input$svalive_ref, input$svalive_comp,
          input$svalive_covariates %||% character(0),
          input$svalive_min_valid_pct, input$svalive_min_variance,
          isTRUE(input$svalive_snp_filter), anno_result()
        )
        m <- sub$m; beta_scale <- sub$beta_scale; grp <- sub$grp; cov_df <- sub$cov_df

        incProgress(0.3, detail = "Estimating surrogate variables (sva::sva)")
        sv_fit <- mod_methyl_sva_fit(m, grp, cov_df)
        design <- sv_fit$design
        validate(need(qr(design)$rank == ncol(design),
          "The selected group/covariate combination produces a rank-deficient design even after dropping surrogate variables. Remove a covariate or change the comparison."))

        incProgress(0.3, detail = "Fitting limma model")
        fit <- methyl_chunked_lmfit(m, design)
        cm <- tryCatch(limma::makeContrasts(contrasts = paste0(input$svalive_comp, "-", input$svalive_ref), levels = design),
                        error = function(e) validate(need(FALSE, "Could not build the reference/comparison contrast - group names may contain characters limma can't use directly (try renaming the group levels in your sample sheet).")))
        fit2 <- tryCatch(limma::eBayes(limma::contrasts.fit(fit, cm)),
                          error = function(e) validate(need(FALSE, paste("limma could not fit this model:", conditionMessage(e)))))
        tt <- limma::topTable(fit2, number = Inf, sort.by = "P")

        incProgress(0.2, detail = "bacon bias/inflation correction")
        bc <- tryCatch(bacon::bacon(teststatistics = tt$t, verbose = FALSE),
                        error = function(e) validate(need(FALSE, paste("bacon could not fit a bias/inflation model on these test statistics:", conditionMessage(e)))))
        p_bacon <- as.numeric(bacon::pval(bc))
        fdr <- stats::p.adjust(p_bacon, method = "BH")

        beta_ref <- rowMeans(beta_scale[, grp == input$svalive_ref, drop = FALSE], na.rm = TRUE)
        beta_comp <- rowMeans(beta_scale[, grp == input$svalive_comp, drop = FALSE], na.rm = TRUE)
        dbeta <- (beta_comp - beta_ref)[rownames(tt)]

        df <- data.frame(
          cpg = rownames(tt), t = tt$t, p_raw = tt$P.Value,
          p_bacon = p_bacon, fdr = fdr,
          dbeta = dbeta, ref_mean_beta = beta_ref[rownames(tt)], comp_mean_beta = beta_comp[rownames(tt)],
          row.names = NULL, stringsAsFactors = FALSE
        )
        ar <- anno_result()
        df$chr <- NA_character_; df$pos <- NA_real_; df$gene <- NA_character_
        if (isTRUE(ar$ok)) {
          a <- ar$anno
          hit <- df$cpg %in% rownames(a)
          df$chr[hit] <- a[df$cpg[hit], "chr"]
          df$pos[hit] <- a[df$cpg[hit], "pos"]
          df$gene[hit] <- a[df$cpg[hit], "gene"]
        }
        df$direction <- ifelse(df$dbeta > 0, "hyper", "hypo")

        cov_names <- if (!is.null(cov_df)) colnames(cov_df) else character(0)
        design_formula <- sprintf("Methylation ~ %s%s + %d surrogate variable(s)", input$svalive_group_col,
                                   if (length(cov_names) > 0) paste0(" + ", paste(cov_names, collapse = " + ")) else "",
                                   sv_fit$n_sv)

        list(
          df = df, ref = input$svalive_ref, comp = input$svalive_comp,
          group_col = input$svalive_group_col, sex_label = sub$sex_label, sex_col = sex_col(),
          covariates = cov_names, n_sv = sv_fit$n_sv, design_formula = design_formula,
          n_after_sex = sub$n_after_sex, n_ref = sub$n_ref, n_comp = sub$n_comp,
          n_probes_tested = nrow(df), n_probes_before_filter = sub$n_probes_before_filter,
          min_valid_pct = input$svalive_min_valid_pct %||% 80, min_variance = input$svalive_min_variance %||% 0,
          snp_filter = isTRUE(input$svalive_snp_filter), snp_note = sub$snp_note,
          missing_note = sub$missing_note, variance_note = sub$variance_note,
          anno_ok = isTRUE(ar$ok), norm_status = dataset_norm_status(),
          dataset_source = methyl_dataset$source %||% "(unnamed)",
          beta_scale = beta_scale, grp = grp,
          lambda_gc = mod_methyl_lambda_gc(df$p_raw),
          lambda_bacon = bacon::inflation(bc), bias_bacon = bacon::bias(bc),
          run_at = Sys.time()
        )
      }
    ), ignoreInit = TRUE)

    observeEvent(svalive_result(), svalive_has_run(TRUE))

    svalive_filtered <- reactive({
      r <- svalive_result()
      direction <- if (identical(input$svalive_direction, "any")) NULL else input$svalive_direction
      mod_methyl_dmp_filter(r$df, "fdr", "dbeta", input$svalive_fdr %||% 0.05, input$svalive_dbeta %||% 0, direction)
    })

    svalive_sig_count <- reactive({
      r <- svalive_result()
      sum(!is.na(r$df$fdr) & r$df$fdr <= (input$svalive_fdr %||% 0.05) & abs(r$df$dbeta) >= (input$svalive_dbeta %||% 0), na.rm = TRUE)
    })

    output$svalive_results_ui <- renderUI({
      req(svalive_has_run(), svalive_result())
      r <- svalive_result()
      n_sig <- svalive_sig_count()
      n_hyper <- sum(!is.na(r$df$fdr) & r$df$fdr <= (input$svalive_fdr %||% 0.05) & r$df$dbeta >= (input$svalive_dbeta %||% 0), na.rm = TRUE)
      n_hypo  <- sum(!is.na(r$df$fdr) & r$df$fdr <= (input$svalive_fdr %||% 0.05) & -r$df$dbeta >= (input$svalive_dbeta %||% 0), na.rm = TRUE)

      tagList(
        div(class = "card",
            div(class = "card-title", icon("circle-check"), "Configuration & sample sizes"),
            p(strong("Model: "), tags$code(r$design_formula), " (limma moderated t-test, eBayes; bacon-corrected)."),
            p(strong("Sex: "), r$sex_label, " | ", strong("Reference: "), r$ref, sprintf(" (n=%d)", r$n_ref),
              " | ", strong("Comparison: "), r$comp, sprintf(" (n=%d)", r$n_comp),
              " | ", strong("Total analyzed: "), r$n_ref + r$n_comp),
            if (length(r$covariates) > 0) p(strong("Covariates: "), paste(r$covariates, collapse = ", ")) else p(class = "submodule-desc", "No covariates selected."),
            p(strong("Surrogate variables: "), r$n_sv,
              if (r$n_sv == 0) " - none were estimated for this cohort; results below are bacon-corrected but not SVA-adjusted." else ""),
            if ((r$n_ref < 10 || r$n_comp < 10))
              p(class = "empty-note", icon("triangle-exclamation"), "One or both groups have fewer than 10 samples - results may be underpowered, and surrogate-variable estimation is less reliable at this scale."),
            p(class = "submodule-desc", sprintf("%s Probes tested after QC filters: %s of %s.",
                                                  paste(c(r$missing_note, r$variance_note, r$snp_note), collapse = " "),
                                                  format(r$n_probes_tested, big.mark = ","), format(r$n_probes_before_filter, big.mark = ",")))
        ),
        div(class = "card",
            div(class = "card-title", icon("gauge-high"), "Genomic inflation diagnostic"),
            p(strong("Genomic inflation factor before bacon (λ): "), if (is.na(r$lambda_gc)) "not available" else sprintf("%.2f", r$lambda_gc)),
            p(strong("bacon inflation: "), sprintf("%.2f", r$lambda_bacon), " | ", strong("bacon bias: "), sprintf("%.3f", r$bias_bacon),
              " - bacon's own bias/inflation estimates on the corrected test statistics; 1.0 inflation and 0 bias indicate a well-calibrated model."),
            withSpinner(plotOutput(ns("svalive_qq"), height = 320), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("chart-simple"), "Summary"),
            div(class = "methyl-stats-row",
              fluidRow(
                valueBox(format(r$n_probes_tested, big.mark = ","), "CpGs tested", icon = icon("dna"), color = "purple", width = 3),
                valueBox(format(n_sig, big.mark = ","), "Significant (FDR + Δβ)", icon = icon("star"), color = if (n_sig > 0) "green" else "light-blue", width = 3),
                valueBox(format(n_hyper, big.mark = ","), "Hypermethylated", icon = icon("arrow-up"), color = "red", width = 3),
                valueBox(format(n_hypo, big.mark = ","), "Hypomethylated", icon = icon("arrow-down"), color = "blue", width = 3)
              )
            ),
            if (n_sig == 0) p(class = "empty-note", icon("circle-info"),
              "No DMPs passed the selected FDR and Δβ thresholds. Consider relaxing the filtering thresholds or reviewing the sample/group configuration.")
        ),
        div(class = "card",
            div(class = "card-title", icon("chart-scatter"), "Volcano plot"),
            withSpinner(plotOutput(ns("svalive_volcano"), height = 340), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("chart-column"), "Manhattan plot"),
            if (r$anno_ok) withSpinner(plotOutput(ns("svalive_manhattan"), height = 320), color = "#2563EB", type = 6)
            else p(class = "empty-note", icon("circle-info"), "No chromosome/position annotation is available for this array type in this deployment.")
        ),
        div(class = "card",
            div(class = "card-title", icon("ranking-star"), "Top DMPs"),
            fluidRow(
              column(6, selectInput(ns("svalive_rank_by"), "Rank by", choices = c("FDR" = "fdr", "Absolute Δβ" = "dbeta", "Combined (FDR/|Δβ|)" = "combined"), selected = "fdr")),
              column(6, selectInput(ns("svalive_top_n"), "Show top", choices = c(10, 20, 50, 100), selected = 20))
            ),
            withSpinner(plotOutput(ns("svalive_topplot"), height = 380), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("chart-simple"), "β-value distribution (top-ranked CpGs)"),
            withSpinner(plotOutput(ns("svalive_betadist"), height = 340), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("table"), "Results table"),
            div(class = "table-toolbar",
                downloadButton(ns("download_svalive"), "Download filtered CSV", class = "btn-default btn-sm"),
                downloadButton(ns("download_svalive_config"), "Download analysis configuration", class = "btn-default btn-sm")),
            DT::dataTableOutput(ns("svalive_table"))
        )
      )
    })

    output$svalive_qq <- renderPlot({
      r <- svalive_result()
      mod_methyl_qq_plot(r$df$p_raw)
    })

    output$svalive_volcano <- renderPlot({
      r <- svalive_result()
      mod_methyl_dmp_volcano(r$df, "dbeta", "fdr", "Δβ", input$svalive_fdr %||% 0.05, input$svalive_dbeta %||% 0)
    })

    output$svalive_manhattan <- renderPlot({
      r <- svalive_result()
      req(r$anno_ok)
      mod_methyl_dmp_manhattan(r$df, input$svalive_fdr %||% 0.05)
    })

    svalive_top <- reactive({
      r <- svalive_result()
      mod_methyl_dmp_topplot(r$df, rank_by = input$svalive_rank_by %||% "fdr", n = as.integer(input$svalive_top_n %||% 20))
    })

    output$svalive_topplot <- renderPlot({ svalive_top()$plot })

    output$svalive_betadist <- renderPlot({
      r <- svalive_result()
      mod_methyl_dmp_betadist(r$beta_scale, svalive_top()$cpgs, r$grp)
    })

    output$svalive_table <- DT::renderDataTable({
      df <- svalive_filtered()
      show_cols <- c("cpg", "gene", "chr", "pos", "ref_mean_beta", "comp_mean_beta", "dbeta", "p_raw", "p_bacon", "fdr", "direction")
      df$significant <- ifelse(!is.na(df$fdr) & df$fdr <= (input$svalive_fdr %||% 0.05) & abs(df$dbeta) >= (input$svalive_dbeta %||% 0), "Yes", "No")
      DT::datatable(df[, c(show_cols, "significant")], rownames = FALSE, filter = "top",
                    options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatSignif(columns = c("ref_mean_beta", "comp_mean_beta", "dbeta", "p_raw", "p_bacon", "fdr"), digits = 4)
    })
    outputOptions(output, "svalive_table", suspendWhenHidden = FALSE)  ## see default_table's own comment above

    output$download_svalive <- downloadHandler(
      filename = function() { r <- svalive_result(); sprintf("sva_dmp_%s_vs_%s_%s.csv", r$comp, r$ref, r$sex_label) },
      content = function(file) utils::write.csv(svalive_filtered(), file, row.names = FALSE)
    )

    output$download_svalive_config <- downloadHandler(
      filename = function() "sva_dmp_analysis_configuration.csv",
      content = function(file) {
        r <- svalive_result()
        cfg <- data.frame(
          parameter = c("dataset", "sex", "reference_group", "comparison_group", "n_reference", "n_comparison",
                        "n_total_analyzed", "covariates", "n_surrogate_variables", "statistical_method", "design_formula",
                        "fdr_threshold", "dbeta_threshold", "direction_filter", "min_valid_sample_pct",
                        "min_variance", "snp_filter_applied", "n_cpgs_tested", "n_significant",
                        "bacon_inflation", "bacon_bias", "run_at"),
          value = c(r$dataset_source, r$sex_label, r$ref, r$comp, r$n_ref, r$n_comp, r$n_ref + r$n_comp,
                    if (length(r$covariates) > 0) paste(r$covariates, collapse = ";") else "(none)",
                    r$n_sv, "sva::sva() + limma (eBayes) + bacon::bacon()", r$design_formula,
                    input$svalive_fdr %||% 0.05, input$svalive_dbeta %||% 0, input$svalive_direction %||% "any",
                    r$min_valid_pct, r$min_variance, r$snp_filter, r$n_probes_tested, svalive_sig_count(),
                    r$lambda_bacon, r$bias_bacon, format(r$run_at)),
          stringsAsFactors = FALSE
        )
        utils::write.csv(cfg, file, row.names = FALSE)
      }
    )

    ## ================= 2. DMP Analysis (configurable live engine) =========
    ## Runs against whatever methyl_dataset$beta currently holds (upload or
    ## preloaded). Controls appear once a matrix + sample sheet are loaded,
    ## but the model only fits on "Run DMP Analysis".

    sex_col <- reactive({ mod_methyl_dmp_sex_col(methyl_dataset$sample_sheet) })

    id_cols <- reactive({
      sheet <- methyl_dataset$sample_sheet
      if (is.null(sheet)) return(character(0))
      intersect(c("sample", "Sample", "sample_id", "Sample_ID", "gsm", "GSM"), colnames(sheet))
    })

    anno_result <- reactive({
      req(methyl_dataset$array_type)
      methyl_get_annotation(methyl_dataset$array_type)
    })

    dataset_norm_status <- reactive({
      req(methyl_dataset$beta)
      methyl_norm_status(methyl_dataset$beta, methyl_dataset, anno_result())
    })

    output$live_ui <- renderUI({
      ## Available whenever a real matrix is loaded (upload or preloaded's live matrix).
      if (is.null(methyl_dataset$beta)) {
        msg <- if (isTRUE(methyl_dataset$preloaded))
          "The preloaded dataset's live beta matrix isn't available in this deployment, so only the sex-stratified reproduced analysis above is available (no live All-Samples/Female-only/Male-only option)."
        else
          "Upload a beta/M-value matrix or IDAT files on the Methylomics Dataset tab to configure and run a live differential methylation model - including an All-Samples (combined) option."
        return(div(class = "card",
          div(class = "card-title", icon("upload"), "DMP Analysis"),
          p(class = "submodule-desc", msg)
        ))
      }
      if (is.null(methyl_dataset$sample_sheet)) {
        return(div(class = "card",
          div(class = "card-title", icon("triangle-exclamation"), "No sample sheet"),
          p(class = "submodule-desc", "A DMP model needs a group variable - re-load on the Dataset tab with a sample sheet/phenotype file included.")
        ))
      }
      sheet <- methyl_dataset$sample_sheet
      cols <- colnames(sheet)
      sc <- sex_col()
      ns_ <- dataset_norm_status()
      arr_ok <- isTRUE(anno_result()$ok)
      tagList(
        div(class = "card",
            div(class = "card-title", icon("flask"), "DMP Analysis"),
            p(class = "empty-note", icon("circle-info"),
              sprintf("Dataset: %s. %s probes × %s samples. %s",
                      methyl_dataset$source %||% "(unnamed)",
                      format(nrow(methyl_dataset$beta), big.mark = ","), ncol(methyl_dataset$beta),
                      ns_$message)),
            fluidRow(
              column(4,
                tags$h5("Sex"),
                radioButtons(ns("live_sex"), NULL, inline = TRUE, choices = mod_methyl_dmp_sex_choices(sheet, sc), selected = "__all__"),
                if (length(mod_methyl_dmp_sex_choices(sheet, sc)) <= 1) p(class = "empty-note", icon("circle-info"), "No usable sex information was found for this dataset - showing pooled analysis only."),

                tags$h5("Comparison"),
                selectInput(ns("live_group_col"), "Group column", choices = cols,
                            selected = intersect(c("group", "Group", "disease", "Disease"), cols)[1] %||% cols[1]),
                uiOutput(ns("live_level_ui")),

                tags$h5("Filters"),
                numericInput(ns("live_fdr"), "Adjusted p-value (FDR) threshold", value = 0.05, min = 0, max = 1, step = 0.01),
                numericInput(ns("live_dbeta"), "Absolute Δβ threshold", value = 0, min = 0, max = 1, step = 0.01),
                radioButtons(ns("live_direction"), "Direction", inline = TRUE,
                             choices = c("All DMPs" = "any", "Hypermethylated" = "hyper", "Hypomethylated" = "hypo"), selected = "any"),
                numericInput(ns("live_min_valid_pct"), "Minimum valid (non-missing) sample %", value = 80, min = 0, max = 100, step = 5),
                numericInput(ns("live_min_variance"), "Minimum methylation variance (optional)", value = 0, min = 0, step = 0.001),
                if (arr_ok) checkboxInput(ns("live_snp_filter"), "Remove SNP-associated probes (manifest Probe_rs/CpG_rs/SBE_rs)", value = FALSE)
                else p(class = "empty-note", icon("circle-info"), anno_result()$reason),

                tags$h5("Covariates (optional)"),
                uiOutput(ns("live_covariate_ui")),

                actionButton(ns("live_run_btn"), "Run DMP Analysis", icon = icon("play"), class = "btn-primary")
              ),
              column(8, withSpinner(uiOutput(ns("live_results_ui")), color = "#2563EB", type = 6))
            )
        )
      )
    })

    output$live_level_ui <- renderUI({
      req(input$live_group_col)
      sheet <- methyl_dataset$sample_sheet
      req(input$live_group_col %in% colnames(sheet))
      levels_available <- sort(unique(as.character(stats::na.omit(sheet[[input$live_group_col]]))))
      validate(need(length(levels_available) >= 2, "This column has fewer than two distinct values - pick a different group column."))
      tagList(
        selectInput(ns("live_ref"), "Reference group", choices = levels_available, selected = levels_available[1]),
        selectInput(ns("live_comp"), "Comparison group", choices = levels_available, selected = levels_available[min(2, length(levels_available))])
      )
    })

    output$live_covariate_ui <- renderUI({
      sheet <- methyl_dataset$sample_sheet
      req(sheet, input$live_group_col)
      exclude <- c(id_cols(), input$live_group_col)
      ## Sex is constant within a single-sex run, so only offer it as a covariate for "All samples".
      if (!identical(input$live_sex %||% "__all__", "__all__") && !is.null(sex_col())) exclude <- c(exclude, sex_col())
      cand <- mod_methyl_dmp_covariate_cols(sheet, exclude)
      if (length(cand) == 0) {
        return(p(class = "empty-note", icon("circle-info"), "No additional phenotype columns are available to use as covariates."))
      }
      checkboxGroupInput(ns("live_covariates"), NULL, choices = cand, selected = character(0))
    })

    ## Reset when the underlying dataset changes so a stale result isn't left showing.
    live_has_run <- reactiveVal(FALSE)
    observeEvent(methyl_dataset$beta, live_has_run(FALSE), ignoreNULL = TRUE)

    live_result <- eventReactive(input$live_run_btn, withProgress(message = "Running DMP Analysis - a full genome-wide array (400k+ CpGs) can take a minute or more...", value = 0.2, {
      sub <- mod_methyl_dmp_prepare_subset(
        methyl_dataset, input$live_sex %||% "__all__", sex_col(),
        input$live_group_col, input$live_ref, input$live_comp,
        input$live_covariates %||% character(0),
        input$live_min_valid_pct, input$live_min_variance,
        isTRUE(input$live_snp_filter), anno_result()
      )
      m <- sub$m; beta_scale <- sub$beta_scale; grp <- sub$grp; cov_df <- sub$cov_df

      ## ---- design + limma fit ---------------------------------------------
      design_grp <- stats::model.matrix(~0 + grp)
      colnames(design_grp) <- levels(grp)
      cov_names <- if (!is.null(cov_df)) colnames(cov_df) else character(0)
      if (length(cov_names) > 0) {
        safe <- sprintf("`%s`", cov_names)
        cov_form <- stats::as.formula(paste("~", paste(safe, collapse = " + ")))
        design_cov <- tryCatch(stats::model.matrix(cov_form, data = cov_df),
                                error = function(e) validate(need(FALSE, paste("Could not build a design matrix for the selected covariates:", conditionMessage(e)))))
        design_cov <- design_cov[, setdiff(colnames(design_cov), "(Intercept)"), drop = FALSE]
        design <- cbind(design_grp, design_cov)
      } else {
        design <- design_grp
      }
      validate(need(qr(design)$rank == ncol(design),
        "The selected group/covariate combination produces a rank-deficient design (e.g. a covariate that perfectly predicts the group, or is constant). Remove a covariate or change the comparison."))

      fit <- methyl_chunked_lmfit(m, design)
      cm <- tryCatch(limma::makeContrasts(contrasts = paste0(input$live_comp, "-", input$live_ref), levels = design),
                      error = function(e) validate(need(FALSE, "Could not build the reference/comparison contrast - group names may contain characters limma can't use directly (try renaming the group levels in your sample sheet).")))
      fit2 <- tryCatch(limma::eBayes(limma::contrasts.fit(fit, cm)),
                        error = function(e) validate(need(FALSE, paste("limma could not fit this model:", conditionMessage(e)))))
      tt <- limma::topTable(fit2, number = Inf, sort.by = "P")

      beta_ref <- rowMeans(beta_scale[, grp == input$live_ref, drop = FALSE], na.rm = TRUE)
      beta_comp <- rowMeans(beta_scale[, grp == input$live_comp, drop = FALSE], na.rm = TRUE)
      dbeta <- (beta_comp - beta_ref)[rownames(tt)]

      df <- data.frame(
        cpg = rownames(tt), t = tt$t, p_raw = tt$P.Value, fdr = tt$adj.P.Val,
        dbeta = dbeta, ref_mean_beta = beta_ref[rownames(tt)], comp_mean_beta = beta_comp[rownames(tt)],
        row.names = NULL, stringsAsFactors = FALSE
      )
      ar <- anno_result()
      df$chr <- NA_character_; df$pos <- NA_real_; df$gene <- NA_character_
      if (isTRUE(ar$ok)) {
        a <- ar$anno
        hit <- df$cpg %in% rownames(a)
        df$chr[hit] <- a[df$cpg[hit], "chr"]
        df$pos[hit] <- a[df$cpg[hit], "pos"]
        df$gene[hit] <- a[df$cpg[hit], "gene"]
      }
      df$direction <- ifelse(df$dbeta > 0, "hyper", "hypo")

      design_formula <- sprintf("Methylation ~ %s%s", input$live_group_col,
                                 if (length(cov_names) > 0) paste0(" + ", paste(cov_names, collapse = " + ")) else "")

      list(
        df = df, ref = input$live_ref, comp = input$live_comp,
        group_col = input$live_group_col, sex_label = sub$sex_label, sex_col = sex_col(),
        covariates = cov_names, design_formula = design_formula,
        n_after_sex = sub$n_after_sex, n_ref = sub$n_ref, n_comp = sub$n_comp,
        n_probes_tested = nrow(df), n_probes_before_filter = sub$n_probes_before_filter,
        min_valid_pct = input$live_min_valid_pct %||% 80, min_variance = input$live_min_variance %||% 0,
        snp_filter = isTRUE(input$live_snp_filter), snp_note = sub$snp_note,
        missing_note = sub$missing_note, variance_note = sub$variance_note,
        anno_ok = isTRUE(ar$ok), norm_status = dataset_norm_status(),
        dataset_source = methyl_dataset$source %||% "(unnamed)", preloaded = isTRUE(methyl_dataset$preloaded),
        beta_scale = beta_scale, grp = grp,
        lambda_gc = mod_methyl_lambda_gc(df$p_raw),
        run_at = Sys.time()
      )
    }), ignoreInit = TRUE)

    observeEvent(live_result(), {
      r <- live_result()
      methyl_results$dmp <- list(comparison = sprintf("%s vs %s (%s)", r$comp, r$ref, r$sex_label),
                                  n_probes = r$n_probes_tested,
                                  n_sig = sum(!is.na(r$df$fdr) & r$df$fdr < 0.05, na.rm = TRUE))
      live_has_run(TRUE)
    })

    ## FDR/Δβ/direction filters apply without refitting; only sex/group/covariate/probe-filter changes need Run again.
    live_filtered <- reactive({
      r <- live_result()
      direction <- if (identical(input$live_direction, "any")) NULL else input$live_direction
      mod_methyl_dmp_filter(r$df, "fdr", "dbeta", input$live_fdr %||% 0.05, input$live_dbeta %||% 0, direction)
    })

    live_sig_count <- reactive({
      r <- live_result()
      sum(!is.na(r$df$fdr) & r$df$fdr <= (input$live_fdr %||% 0.05) & abs(r$df$dbeta) >= (input$live_dbeta %||% 0), na.rm = TRUE)
    })

    output$live_results_ui <- renderUI({
      req(live_has_run(), live_result())
      r <- live_result()
      n_sig <- live_sig_count()
      n_hyper <- sum(!is.na(r$df$fdr) & r$df$fdr <= (input$live_fdr %||% 0.05) & r$df$dbeta >= (input$live_dbeta %||% 0), na.rm = TRUE)
      n_hypo  <- sum(!is.na(r$df$fdr) & r$df$fdr <= (input$live_fdr %||% 0.05) & -r$df$dbeta >= (input$live_dbeta %||% 0), na.rm = TRUE)
      n_pass_fdr <- sum(!is.na(r$df$fdr) & r$df$fdr <= (input$live_fdr %||% 0.05), na.rm = TRUE)
      n_pass_dbeta <- sum(!is.na(r$df$dbeta) & abs(r$df$dbeta) >= (input$live_dbeta %||% 0), na.rm = TRUE)

      tagList(
        div(class = "card",
            div(class = "card-title", icon("circle-check"), "Configuration & sample sizes"),
            p(strong("Model: "), tags$code(r$design_formula), " (limma moderated t-test, eBayes)."),
            p(strong("Sex: "), r$sex_label, " | ", strong("Reference: "), r$ref, sprintf(" (n=%d)", r$n_ref),
              " | ", strong("Comparison: "), r$comp, sprintf(" (n=%d)", r$n_comp),
              " | ", strong("Total analyzed: "), r$n_ref + r$n_comp),
            if (length(r$covariates) > 0) p(strong("Covariates: "), paste(r$covariates, collapse = ", ")) else p(class = "submodule-desc", "No covariates selected."),
            p(class = "submodule-desc", r$norm_status$message),
            if ((r$n_ref < 10 || r$n_comp < 10))
              p(class = "empty-note", icon("triangle-exclamation"), "One or both groups have fewer than 10 samples - results may be underpowered."),
            p(class = "submodule-desc", sprintf("%s Probes tested after QC filters: %s of %s.",
                                                  paste(c(r$missing_note, r$variance_note, r$snp_note), collapse = " "),
                                                  format(r$n_probes_tested, big.mark = ","), format(r$n_probes_before_filter, big.mark = ",")))
        ),
        div(class = "card",
            div(class = "card-title", icon("gauge-high"), "Genomic inflation diagnostic"),
            p(strong("Genomic inflation factor (λ): "), if (is.na(r$lambda_gc)) "not available" else sprintf("%.2f", r$lambda_gc),
              " - the ratio of observed to expected median test statistic across all tested CpGs; 1.0 indicates no inflation."),
            if (!is.na(r$lambda_gc) && r$lambda_gc > 1.1)
              p(class = "empty-note", icon("triangle-exclamation"),
                "Genomic inflation is elevated (λ > 1.1) - this live engine doesn't apply SVA/bacon correction, so treat these p-values/FDR as potentially optimistic."),
            withSpinner(plotOutput(ns("live_qq"), height = 320), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("chart-simple"), "Summary"),
            div(class = "methyl-stats-row",
              fluidRow(
                valueBox(format(r$n_probes_tested, big.mark = ","), "CpGs tested", icon = icon("dna"), color = "purple", width = 3),
                valueBox(format(n_sig, big.mark = ","), "Significant (FDR + Δβ)", icon = icon("star"), color = if (n_sig > 0) "green" else "light-blue", width = 3),
                valueBox(format(n_hyper, big.mark = ","), "Hypermethylated", icon = icon("arrow-up"), color = "red", width = 3),
                valueBox(format(n_hypo, big.mark = ","), "Hypomethylated", icon = icon("arrow-down"), color = "blue", width = 3)
              ),
              fluidRow(
                valueBox(format(n_pass_fdr, big.mark = ","), sprintf("Pass FDR < %s", input$live_fdr %||% 0.05), icon = icon("filter"), color = "light-blue", width = 4),
                valueBox(format(n_pass_dbeta, big.mark = ","), sprintf("Pass |Δβ| ≥ %s", input$live_dbeta %||% 0), icon = icon("filter"), color = "light-blue", width = 4),
                valueBox(r$sex_label, "Sex subset", icon = icon("venus-mars"), color = "black", width = 4)
              )
            ),
            if (n_sig == 0) p(class = "empty-note", icon("circle-info"),
              "No DMPs passed the selected FDR and Δβ thresholds. Consider relaxing the filtering thresholds or reviewing the sample/group configuration.")
        ),
        div(class = "card",
            div(class = "card-title", icon("chart-scatter"), "Volcano plot"),
            withSpinner(plotOutput(ns("live_volcano"), height = 340), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("chart-column"), "Manhattan plot"),
            if (r$anno_ok) withSpinner(plotOutput(ns("live_manhattan"), height = 320), color = "#2563EB", type = 6)
            else p(class = "empty-note", icon("circle-info"), "No chromosome/position annotation is available for this array type in this deployment.")
        ),
        div(class = "card",
            div(class = "card-title", icon("ranking-star"), "Top DMPs"),
            fluidRow(
              column(6, selectInput(ns("live_rank_by"), "Rank by", choices = c("FDR" = "fdr", "Absolute Δβ" = "dbeta", "Combined (FDR/|Δβ|)" = "combined"), selected = "fdr")),
              column(6, selectInput(ns("live_top_n"), "Show top", choices = c(10, 20, 50, 100), selected = 20))
            ),
            withSpinner(plotOutput(ns("live_topplot"), height = 380), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("chart-simple"), "β-value distribution (top-ranked CpGs)"),
            withSpinner(plotOutput(ns("live_betadist"), height = 340), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("table"), "Results table"),
            div(class = "table-toolbar",
                downloadButton(ns("download_live"), "Download filtered CSV", class = "btn-default btn-sm"),
                downloadButton(ns("download_live_config"), "Download analysis configuration", class = "btn-default btn-sm"),
                downloadButton(ns("download_provenance"), "Download analysis record (.json)", class = "btn-default btn-sm")),
            DT::dataTableOutput(ns("live_table"))
        )
      )
    })

    output$live_qq <- renderPlot({
      r <- live_result()
      mod_methyl_qq_plot(r$df$p_raw)
    })

    output$live_volcano <- renderPlot({
      r <- live_result()
      mod_methyl_dmp_volcano(r$df, "dbeta", "fdr", "Δβ", input$live_fdr %||% 0.05, input$live_dbeta %||% 0)
    })

    output$live_manhattan <- renderPlot({
      r <- live_result()
      req(r$anno_ok)
      mod_methyl_dmp_manhattan(r$df, input$live_fdr %||% 0.05)
    })

    live_top <- reactive({
      r <- live_result()
      mod_methyl_dmp_topplot(r$df, rank_by = input$live_rank_by %||% "fdr", n = as.integer(input$live_top_n %||% 20))
    })

    output$live_topplot <- renderPlot({ live_top()$plot })

    output$live_betadist <- renderPlot({
      r <- live_result()
      mod_methyl_dmp_betadist(r$beta_scale, live_top()$cpgs, r$grp)
    })

    output$live_table <- DT::renderDataTable({
      df <- live_filtered()
      show_cols <- c("cpg", "gene", "chr", "pos", "ref_mean_beta", "comp_mean_beta", "dbeta", "p_raw", "fdr", "direction")
      df$significant <- ifelse(!is.na(df$fdr) & df$fdr <= (input$live_fdr %||% 0.05) & abs(df$dbeta) >= (input$live_dbeta %||% 0), "Yes", "No")
      DT::datatable(df[, c(show_cols, "significant")], rownames = FALSE, filter = "top",
                    options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact") %>%
        DT::formatSignif(columns = c("ref_mean_beta", "comp_mean_beta", "dbeta", "p_raw", "fdr"), digits = 4)
    })
    outputOptions(output, "live_table", suspendWhenHidden = FALSE)  ## see default_table's own comment above

    output$download_live <- downloadHandler(
      filename = function() { r <- live_result(); sprintf("dmp_%s_vs_%s_%s.csv", r$comp, r$ref, r$sex_label) },
      content = function(file) utils::write.csv(live_filtered(), file, row.names = FALSE)
    )

    output$download_live_config <- downloadHandler(
      filename = function() "dmp_analysis_configuration.csv",
      content = function(file) {
        r <- live_result()
        cfg <- data.frame(
          parameter = c("dataset", "sex", "reference_group", "comparison_group", "n_reference", "n_comparison",
                        "n_total_analyzed", "covariates", "statistical_method", "design_formula",
                        "fdr_threshold", "dbeta_threshold", "direction_filter", "min_valid_sample_pct",
                        "min_variance", "snp_filter_applied", "n_cpgs_tested", "n_significant", "run_at"),
          value = c(r$dataset_source, r$sex_label, r$ref, r$comp, r$n_ref, r$n_comp, r$n_ref + r$n_comp,
                    if (length(r$covariates) > 0) paste(r$covariates, collapse = ";") else "(none)",
                    "limma (moderated t-test, eBayes)", r$design_formula,
                    input$live_fdr %||% 0.05, input$live_dbeta %||% 0, input$live_direction %||% "any",
                    r$min_valid_pct, r$min_variance, r$snp_filter, r$n_probes_tested, live_sig_count(),
                    format(r$run_at)),
          stringsAsFactors = FALSE
        )
        utils::write.csv(cfg, file, row.names = FALSE)
      }
    )

    ## Provenance manifest (R/provenance.R): checksum of the exact filtered
    ## beta/M-value matrix + group assignment that went into the limma fit
    ## (m/grp, captured inside live_result() above), plus the group/
    ## covariate/threshold choices already surfaced in download_live_config's
    ## own CSV above. No seed - this live DMP engine's only statistical
    ## step is limma::eBayes(), which is deterministic (no internal RNG).
    dmp_provenance_record <- reactive({
      r <- live_result()
      arthomix_provenance_record(
        module = "mod_methyl_dmp",
        checksum_input = list(beta_scale = r$beta_scale, grp = as.character(r$grp)),
        params = list(
          dataset_source = r$dataset_source, sex = r$sex_label,
          group_col = r$group_col, reference_level = r$ref, comparison_level = r$comp,
          covariates = r$covariates, design_formula = r$design_formula,
          fdr_threshold = input$live_fdr %||% 0.05, dbeta_threshold = input$live_dbeta %||% 0,
          direction_filter = input$live_direction %||% "any",
          min_valid_sample_pct = r$min_valid_pct, min_variance = r$min_variance,
          snp_filter_applied = r$snp_filter,
          n_reference = r$n_ref, n_comparison = r$n_comp,
          n_probes_tested = r$n_probes_tested, n_probes_before_filter = r$n_probes_before_filter
        ),
        seed = NULL,
        packages = "limma",
        extra = list(lambda_gc = r$lambda_gc, n_significant = live_sig_count())
      )
    })

    output$download_provenance <- arthomix_provenance_download_handler(dmp_provenance_record, "mod_methyl_dmp_provenance")
  })
}
