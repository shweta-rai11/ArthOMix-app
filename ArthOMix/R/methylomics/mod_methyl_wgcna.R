## R/methylomics/mod_methyl_wgcna.R
## Submodule: WGCNA (Co-Methylation Network). Self-contained methylomics
## implementation - CpGs as network nodes, samples as observations. Every
## helper below is mx_wgcna_*-prefixed and lives only in this file; nothing
## here is shared with, or can collide with, R/transcriptomics/mod_wgcna.R's
## own wgcna_*-prefixed helpers (a different, unrelated gene-expression
## implementation this file never reads from or writes to).
##
## Two data pathways feed the SAME live pipeline (Data & Filtering -> Sample
## QC -> Soft Threshold -> Network & Modules -> Module-Trait -> Hub CpGs ->
## Results & Export); only the starting matrix and each control's DEFAULT
## value differ - every control stays user-editable in both paths:
##  - Preloaded: methyl_dataset$beta/$sample_sheet from the preloaded
##    whole-blood cohort (GSE42861, Liu et al. 2013), when the live raw
##    matrix is available in this deployment. Defaults reproduce
##    script05_wgcna_sexstratified/05_wgcna_sexstratified.R's actual
##    published analysis: per-sex residualization against age+smoking+
##    cell-type composition (minus the highest-mean cell-type column, kept
##    as an implicit reference), top 20,000 CpGs by residual MAD, signed
##    network, Pearson correlation, automatic power selection via the
##    script's own "first power reaching the R-squared cutoff, else the
##    best-observed power" rule over its exact tested power vector,
##    minModuleSize=20, mergeCutHeight=0.25, deepSplit=2, maxBlockSize=5000.
##  - Uploaded: whatever methyl_dataset$beta/$sample_sheet an upload
##    produced. Defaults follow general WGCNA-on-methylation best practice
##    (bicor, signed, auto power 1-20, minModuleSize=30) instead.
##
## Sex stratification (Female / Male / All samples combined) is a
## first-class, always-visible control above the sub-tabs, read at the top
## of Data & Filtering's own filter step (not a separate stage) - the
## published script residualizes and MAD-ranks CpGs SEPARATELY per sex, so
## reproducing it faithfully means subsetting to one sex BEFORE
## residualization/variability-ranking, not filtering globally and
## subsetting afterward.
##
## Results-visibility: every stage after Data & Filtering is gated behind
## its own action button via the same idiom R/transcriptomics/mod_wgcna.R
## uses - a stage's renderUI does tryCatch(<eventReactive>(), error =
## function(e) NULL) and shows only an empty-note "not run yet" message
## when NULL. Nothing computes or renders speculatively.

mod_methyl_wgcna_config <- list(
  id = "wgcna", title = "WGCNA (Co-Methylation Network)", icon = "circle-nodes", group = "Network",
  description = "Weighted co-methylation network analysis: filtering, soft-threshold selection, module detection, module-trait correlation, hub CpGs, and enrichment. Works on the preloaded dataset or your own uploaded, normalized data."
)

## ---- small pure helpers (all mx_wgcna_*-prefixed) --------------------------

mx_wgcna_tab_title <- function(ic, label) tagList(icon(ic), " ", label)

## Row-wise MAD via matrixStats when available - same performance reasoning
## methyl_row_vars() in qc.R already documents for rowVars() at 400k+ probes.
mx_wgcna_row_mads <- function(m) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowMads(m, na.rm = TRUE)
  } else {
    apply(m, 1, function(x) stats::mad(x, na.rm = TRUE))
  }
}

## Ranks probes by the chosen variability statistic and keeps the top N -
## qc.R's own filters return a keep vector for a fixed threshold; this is
## the top-N-by-rank companion Data & Filtering needs. `mat` is probes
## (rows) x samples (cols).
mx_wgcna_top_variable <- function(mat, method = c("mad", "variance", "sd", "iqr"), top_n) {
  method <- match.arg(method)
  stat <- switch(method,
    mad      = mx_wgcna_row_mads(mat),
    variance = methyl_row_vars(mat),
    sd       = sqrt(pmax(methyl_row_vars(mat), 0)),
    iqr      = apply(mat, 1, function(x) stats::IQR(x, na.rm = TRUE))
  )
  stat[!is.finite(stat)] <- 0
  keep_n <- max(1, min(top_n, sum(stat > 0)))
  ord <- order(stat, decreasing = TRUE)[seq_len(keep_n)]
  list(mat = mat[ord, , drop = FALSE], method = method, top_n = top_n, n_kept = keep_n)
}

## Encodes a trait column as a numeric vector for module-eigengene
## correlation: numeric columns pass through; a two-level column becomes
## 0/1; anything else is unsupported (module-trait correlation is a linear
## step, not a multi-level test). Independent from (does not call)
## mod_wgcna.R's own wgcna_encode_trait() - see file header.
mx_wgcna_encode_trait <- function(sheet, col) {
  v <- sheet[[col]]
  if (is.numeric(v)) return(list(ok = TRUE, vec = as.numeric(v)))
  lv <- sort(unique(as.character(stats::na.omit(v))))
  if (length(lv) != 2) {
    return(list(ok = FALSE, reason = sprintf(
      "\"%s\" has %d non-missing level(s) - module-trait correlation needs a numeric column or exactly two levels.", col, length(lv))))
  }
  list(ok = TRUE, vec = as.numeric(factor(as.character(v), levels = lv)) - 1, levels = lv)
}

## Detects cell-type-composition columns among candidate covariates (>=3
## numeric columns whose per-row values sum to ~1, the compositional
## signature Houseman-style cell estimates have) and returns the name of
## the one with the highest mean fraction - the published script's own
## implicit-reference choice ("cell-type covariates: %s used as implicit
## reference"), reproduced generically (never hardcoding literal column
## names) so an uploaded dataset with its own differently-named cell-type
## estimates is handled the same way.
mx_wgcna_celltype_reference <- function(sheet, candidate_cols) {
  if (length(candidate_cols) < 3) return(NULL)
  is_num <- vapply(candidate_cols, function(cl) is.numeric(sheet[[cl]]), logical(1))
  num_cols <- candidate_cols[is_num]
  if (length(num_cols) < 3) return(NULL)
  mat <- vapply(num_cols, function(cl) as.numeric(sheet[[cl]]), numeric(nrow(sheet)))
  row_sums <- rowSums(mat, na.rm = TRUE)
  if (mean(abs(row_sums - 1) < 0.1, na.rm = TRUE) < 0.8) return(NULL)
  means <- colMeans(mat, na.rm = TRUE)
  names(which.max(means))
}

## Scientific-guardrail booleans (spec: warn, don't silently force a
## meaningless analysis through). Any argument left NULL is skipped.
mx_wgcna_guardrails <- function(n_samples = NULL, n_probes = NULL, max_r_sq = NULL, module_colors = NULL) {
  list(
    low_n_samples = !is.null(n_samples) && n_samples < 15,
    low_n_probes  = !is.null(n_probes) && n_probes < 500,
    poor_sft_fit  = !is.null(max_r_sq) && is.finite(max_r_sq) && max_r_sq < 0.7,
    all_grey      = !is.null(module_colors) && length(setdiff(unique(module_colors), "grey")) == 0,
    single_module = !is.null(module_colors) && length(setdiff(unique(module_colors), "grey")) == 1
  )
}

mod_methyl_wgcna_ui <- function(id) {
  ns <- NS(id)
  tagList(
    shinyjs::useShinyjs(),
    withSpinner(uiOutput(ns("body_ui")), color = "#2563EB", type = 6)
  )
}

mod_methyl_wgcna_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    methyl_dataset <- dataset

    ## ==== 0. Shared context reactives ====================================

    sex_col <- reactive(mod_methyl_dmp_sex_col(methyl_dataset$sample_sheet))
    sex_choices_r <- reactive(mod_methyl_dmp_sex_choices(methyl_dataset$sample_sheet, sex_col()))
    mx_default_sex <- reactive({
      ch <- sex_choices_r()
      if (isTRUE(methyl_dataset$preloaded) && "Female" %in% names(ch)) return(unname(ch["Female"]))
      "__all__"
    })

    id_cols <- reactive({
      sheet <- methyl_dataset$sample_sheet
      if (is.null(sheet)) return(character(0))
      intersect(c("sample", "Sample", "sample_id", "Sample_ID", "gsm", "GSM"), colnames(sheet))
    })
    trait_col_default <- reactive({
      sheet <- methyl_dataset$sample_sheet
      if (is.null(sheet)) return(NULL)
      cols <- intersect(c("group", "Group", "disease", "Disease"), colnames(sheet))
      if (length(cols) == 0) return(NULL)
      cols[1]
    })
    covariate_choices <- reactive({
      sheet <- methyl_dataset$sample_sheet
      if (is.null(sheet)) return(character(0))
      exclude <- c(id_cols(), trait_col_default())
      if (!identical(input$sex_stratum %||% "__all__", "__all__") && !is.null(sex_col())) exclude <- c(exclude, sex_col())
      mod_methyl_dmp_covariate_cols(sheet, exclude)
    })
    ## Preloaded default residualization set: age + smoking + cell-type
    ## columns minus the auto-detected reference cell type - reproduces the
    ## published script's own choice (see mx_wgcna_celltype_reference()).
    ## Uploaded default: nothing pre-ticked.
    covariate_default_selected <- reactive({
      cand <- covariate_choices()
      if (!isTRUE(methyl_dataset$preloaded)) return(character(0))
      sheet <- methyl_dataset$sample_sheet
      ref_ct <- mx_wgcna_celltype_reference(sheet, cand)
      age_smoke <- cand[grepl("^age$", cand, ignore.case = TRUE) | grepl("smok", cand, ignore.case = TRUE)]
      ct_cols <- cand[grepl("^(B|NK|CD4T|CD8T|Mono|Neutro|Eosino)$", cand, ignore.case = TRUE)]
      sel <- unique(c(age_smoke, ct_cols))
      if (!is.null(ref_ct)) sel <- setdiff(sel, ref_ct)
      sel
    })

    ## ==== Status strip (always visible, no button) ========================

    status_ui <- function() {
      if (!is.null(methyl_dataset$beta)) {
        div(class = "empty-note", icon("flask"),
            sprintf("%s: %s CpGs × %s samples.",
                    if (isTRUE(methyl_dataset$preloaded)) "Preloaded whole-blood dataset (live matrix)" else "Uploaded dataset",
                    format(nrow(methyl_dataset$beta), big.mark = ","), ncol(methyl_dataset$beta)))
      } else {
        div(class = "empty-note", icon("circle-info"),
            "Preloaded whole-blood dataset (metadata only) - the live beta matrix isn't available in this deployment, so live computation is disabled here. The Results & Export tab's \"Compare with published results\" panel still works from static reference tables.")
      }
    }

    ## ==== 1. Data & Filtering =============================================

    output$orientation_check_ui <- renderUI({
      req(methyl_dataset$beta)
      mat <- methyl_dataset$beta
      if (nrow(mat) < ncol(mat)) {
        p(class = "empty-note", icon("triangle-exclamation"), sprintf(
          "This matrix has more columns (%d) than rows (%d) - CpG methylation matrices normally have far more probes than samples. If this is actually samples-in-rows, tick \"transpose\" below.",
          ncol(mat), nrow(mat)))
      } else NULL
    })

    output$covariate_ui <- renderUI({
      sheet <- methyl_dataset$sample_sheet
      if (is.null(sheet)) return(p(class = "empty-note", icon("circle-info"), "No sample sheet available - residualization is unavailable without phenotype metadata."))
      cand <- covariate_choices()
      if (length(cand) == 0) return(p(class = "empty-note", icon("circle-info"), "No usable covariate columns detected in the sample sheet."))
      ref_ct <- mx_wgcna_celltype_reference(sheet, cand)
      tagList(
        checkboxGroupInput(ns("resid_covariates"), NULL, choices = cand, selected = covariate_default_selected()),
        if (!is.null(ref_ct)) p(class = "empty-note", icon("circle-info"), sprintf(
          "\"%s\" (the cell-type column with the highest mean fraction) is left unticked by default and used as an implicit reference, avoiding a rank-deficient design when the other cell-type columns are included - matches the published script's own approach.", ref_ct))
      )
    })

    tab_data_ui <- function() {
      tagList(
        div(class = "card",
            div(class = "card-title", icon("filter"), "Data & Filtering"),
            p(class = "submodule-desc", "CpGs are the network features (nodes), samples are observations. Nothing below is applied until you click \"Build filtered matrix\"."),
            uiOutput(ns("orientation_check_ui")),
            fluidRow(
              column(6,
                h5("Missingness"),
                numericInput(ns("max_probe_missing"), "Max missingness per CpG (%)", value = 5, min = 0, max = 100, step = 1),
                checkboxInput(ns("force_transpose"), "Matrix looks transposed - transpose before use", value = FALSE)
              ),
              column(6,
                h5("Variability filtering"),
                selectInput(ns("var_method"), "Method", choices = c("MAD" = "mad", "Variance" = "variance", "Standard deviation" = "sd", "IQR" = "iqr"), selected = "mad"),
                numericInput(ns("top_n"), "Number of most variable CpGs", value = if (isTRUE(methyl_dataset$preloaded)) 20000 else 5000, min = 100, step = 500)
              )
            ),
            h5("Optional covariate residualization (limma)"),
            p(class = "empty-note", icon("circle-info"), "Regresses out the selected covariates (on the M-value scale) before ranking CpGs by variability - the group/disease-status column is never offered here, so any group-associated signal stays in the residuals for Module-Trait Analysis to detect later."),
            uiOutput(ns("covariate_ui")),
            actionButton(ns("filter_btn"), "Build filtered matrix", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        withSpinner(uiOutput(ns("filter_result_ui")), color = "#2563EB", type = 6)
      )
    }

    ## Sex-stratum subsetting happens FIRST, inside this stage, not as a
    ## separate step afterward - the published script residualizes and
    ## MAD-ranks CpGs SEPARATELY per sex (build_residuals(), called once per
    ## stratum), so filtering globally then splitting by sex afterward would
    ## select a different (less faithful) CpG set per stratum than the
    ## script's own per-sex ranking.
    ## Memory discipline below mirrors mod_methyl_dmp.R's own live_result()
    ## engine (see its header comment on methyl_chunked_lmfit) - at the full
    ## ~412K-probe preloaded matrix, holding more than one or two full-size
    ## (probes x samples) copies alive at once is enough to exceed a 16GB
    ## machine's vector memory limit. Every intermediate is rm()'d + gc()'d
    ## as soon as it's no longer needed, in-place clipping is used instead of
    ## pmin()/pmax() copies, and the residual subtraction (the single most
    ## memory-hungry step - a naive `m - coefs %*% t(design)` briefly holds
    ## THREE full-size copies at once) is done in row-chunks instead.
    mx_wgcna_filtered <- eventReactive(input$filter_btn, {
      validate(need(!is.null(methyl_dataset$beta), "No live beta matrix available - see the status message above."))
      n_probes_total <- nrow(methyl_dataset$beta)
      withProgress(message = "Building the filtered matrix - a full genome-wide probe set can take a minute or more...", value = 0.1, {
        mat_full <- methyl_dataset$beta
        if (isTRUE(input$force_transpose)) mat_full <- t(mat_full)
        validate(need(nrow(mat_full) >= ncol(mat_full), "This matrix has more samples than probes even after transposing - check the orientation setting."))

        sheet <- methyl_dataset$sample_sheet
        sub <- methyl_qc_subgroup_filter(mat_full, sheet, sex_col(), input$sex_stratum %||% "__all__", min_n = 15)
        rm(mat_full)
        mat <- sub$mat; stratum_label <- sub$label
        rm(sub); gc(FALSE)
        validate(need(ncol(mat) >= 6, sprintf("Only %d sample(s) in the \"%s\" stratum - fewer than the minimum of 6 needed to fit anything.", ncol(mat), stratum_label)))
        incProgress(0.1, detail = "checking scale")

        is_m_scale <- identical(methyl_dataset$input_scale, "m")
        if (!is_m_scale) {
          rng <- range(mat, na.rm = TRUE)
          if (rng[1] < -0.1 || rng[2] > 1.1) {
            if (rng[2] <= 100 && rng[2] > 1.1) {
              validate(need(FALSE, sprintf(
                "Values range from %.2f to %.2f - this looks like a 0-100 percentage scale, not 0-1 beta values. Rescale to 0-1 before uploading, or mark the dataset's scale correctly on the Dataset tab.", rng[1], rng[2])))
            } else {
              validate(need(FALSE, sprintf(
                "Values range from %.2f to %.2f, well outside the expected 0-1 beta-value range. If this is M-value or log-transformed data, mark it as such on the Dataset tab; otherwise this doesn't look like methylation beta-value data.", rng[1], rng[2])))
            }
          }
          mat[mat < 1e-6] <- 1e-6
          mat[mat > 1 - 1e-6] <- 1 - 1e-6
        }
        incProgress(0.1, detail = "missingness filter")

        f_miss <- methyl_filter_missing(mat, (input$max_probe_missing %||% 5) / 100)
        mat <- mat[f_miss$keep, , drop = FALSE]
        validate(need(nrow(mat) >= 50, "Fewer than 50 CpGs remain after the missingness filter - relax the threshold."))

        m <- if (is_m_scale) mat else log2(mat / (1 - mat))
        rm(mat); gc(FALSE)
        incProgress(0.15, detail = "residualization")

        resid_cols <- input$resid_covariates %||% character(0)
        resid_note <- "No residualization applied."
        if (length(resid_cols) > 0) {
          ids <- methyl_sheet_sample_ids(sheet, colnames(m))
          common <- intersect(colnames(m), ids)
          validate(need(length(common) >= 10, "Fewer than 10 samples match between the matrix and the sample sheet for residualization."))
          m <- m[, common, drop = FALSE]
          ph <- as.data.frame(sheet)[match(common, ids), , drop = FALSE]
          cc <- ph[, resid_cols, drop = FALSE]
          complete <- stats::complete.cases(cc)
          validate(need(sum(complete) >= 10, "Fewer than 10 samples have no missing values in the selected covariates - untick one, or leave a covariate deselected."))
          m <- m[, complete, drop = FALSE]
          ph <- ph[complete, , drop = FALSE]
          design <- tryCatch(
            stats::model.matrix(stats::as.formula(paste("~", paste(sprintf("`%s`", resid_cols), collapse = " + "))), data = ph),
            error = function(e) validate(need(FALSE, paste("Could not build a design matrix for the selected covariates:", conditionMessage(e)))))
          validate(need(qr(design)$rank == ncol(design),
            "The selected covariates produce a rank-deficient design (e.g. cell-type fractions that sum to a constant) - untick one, or leave the auto-suggested reference column unticked."))
          fit <- methyl_chunked_lmfit(m, design)
          design_t <- t(design)
          chunk <- 20000L
          starts <- seq(1L, nrow(m), by = chunk)
          for (s in starts) {
            e <- min(s + chunk - 1L, nrow(m))
            m[s:e, ] <- m[s:e, ] - fit$coefficients[s:e, , drop = FALSE] %*% design_t
          }
          rm(fit, design, design_t, ph, cc); gc(FALSE)
          resid_note <- sprintf("Residualized against: %s (n=%d samples).", paste(resid_cols, collapse = ", "), ncol(m))
        }
        incProgress(0.2, detail = "ranking by variability")

        vt <- mx_wgcna_top_variable(m, method = input$var_method %||% "mad", top_n = input$top_n %||% 5000)
        rm(m); gc(FALSE)
        m <- vt$mat
        incProgress(0.1, detail = "sample/probe QC")

        gsg <- tryCatch(WGCNA::goodSamplesGenes(t(m), verbose = 0), error = function(e) NULL)
        gsg_note <- NULL
        if (!is.null(gsg) && !isTRUE(gsg$allOK)) {
          n_bad_g <- sum(!gsg$goodGenes); n_bad_s <- sum(!gsg$goodSamples)
          m <- m[gsg$goodGenes, gsg$goodSamples, drop = FALSE]
          gsg_note <- sprintf("WGCNA::goodSamplesGenes() flagged %d probe(s) and %d sample(s) as unusable (near-zero variance / excess missingness); removed automatically.", n_bad_g, n_bad_s)
        }
        validate(need(ncol(m) >= 6, "Fewer than 6 samples remain in the filtered matrix."))
        validate(need(nrow(m) >= 20, "Fewer than 20 CpGs remain in the filtered matrix."))

        list(mat = m, stratum_label = stratum_label, stratum_low_n = ncol(m) < 15,
             n_probes_in = n_probes_total, n_probes_kept = nrow(m), n_samples_kept = ncol(m),
             missing_note = f_miss$note, resid_note = resid_note, resid_covariates = resid_cols,
             var_method = vt$method, top_n = vt$top_n, gsg_note = gsg_note)
      })
    })

    output$filter_result_ui <- renderUI({
      f <- tryCatch(mx_wgcna_filtered(), error = function(e) NULL)
      if (is.null(f)) return(div(class = "empty-note", icon("circle-info"), "Not built yet. Set filters above, then click \"Build filtered matrix\"."))
      tagList(
        div(class = "card",
            div(class = "card-title", icon("check"), "Filtering Summary"),
            fluidRow(
              valueBox(format(f$n_probes_in, big.mark = ","), "CpGs before filtering", icon = icon("dna"), color = "aqua", width = 3),
              valueBox(format(f$n_probes_kept, big.mark = ","), "CpGs after filtering", icon = icon("dna"), color = "green", width = 3),
              valueBox(f$n_samples_kept, sprintf("Samples (%s)", f$stratum_label), icon = icon("vial"), color = "purple", width = 3),
              valueBox(toupper(f$var_method), "Variability method", icon = icon("chart-line"), color = "light-blue", width = 3)
            ),
            p(class = "empty-note", icon("circle-info"), f$missing_note),
            p(class = "empty-note", icon("circle-info"), f$resid_note),
            if (!is.null(f$gsg_note)) p(class = "empty-note", icon("triangle-exclamation"), f$gsg_note),
            if (isTRUE(f$stratum_low_n)) p(class = "empty-note", icon("triangle-exclamation"), sprintf(
              "Only %d samples in this stratum - WGCNA module detection is not reliable below ~15 samples. Consider \"All samples combined\".", f$n_samples_kept))
        )
      )
    })

    ## ==== 2. Sample QC (read-only, no button - reacts to the filtered matrix) ==

    tab_sampleqc_ui <- function() {
      f <- tryCatch(mx_wgcna_filtered(), error = function(e) NULL)
      if (is.null(f)) return(div(class = "empty-note", icon("circle-info"), "Build the filtered matrix on \"Data & Filtering\" first."))
      gr <- mx_wgcna_guardrails(n_samples = f$n_samples_kept, n_probes = f$n_probes_kept)
      div(class = "card",
          div(class = "card-title", icon("magnifying-glass-chart"), "Sample QC"),
          p(class = "empty-note", icon("circle-info"), sprintf("Stratum: %s. %s CpGs × %s samples.", f$stratum_label, format(f$n_probes_kept, big.mark = ","), f$n_samples_kept)),
          if (gr$low_n_samples) p(class = "empty-note", icon("triangle-exclamation"), sprintf("Only %d samples - module detection is not reliable below ~15 samples.", f$n_samples_kept)),
          fluidRow(
            column(6, withSpinner(plotOutput(ns("sampleqc_dendro_plot"), height = 320), color = "#2563EB", type = 6)),
            column(6, withSpinner(plotOutput(ns("sampleqc_pca_plot"), height = 320), color = "#2563EB", type = 6))
          ),
          withSpinner(plotOutput(ns("sampleqc_missing_plot"), height = 240), color = "#2563EB", type = 6)
      )
    }

    output$sampleqc_dendro_plot <- renderPlot({
      f <- tryCatch(mx_wgcna_filtered(), error = function(e) NULL); req(f)
      tr <- stats::hclust(stats::dist(t(f$mat)), method = "average")
      graphics::plot(tr, main = "", xlab = "", sub = "", cex = 0.6)
    })
    output$sampleqc_pca_plot <- renderPlot({
      f <- tryCatch(mx_wgcna_filtered(), error = function(e) NULL); req(f)
      validate(need(ncol(f$mat) >= 3 && nrow(f$mat) >= 3, "Not enough samples/probes for a PCA."))
      pca <- stats::prcomp(t(f$mat), scale. = TRUE)
      df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
      ggplot(df, aes(x = PC1, y = PC2)) + geom_point(color = ARTHOMIX_COLORS$blue, size = 2) +
        labs(x = "PC1", y = "PC2") + theme_arthomix()
    })
    output$sampleqc_missing_plot <- renderPlot({
      f <- tryCatch(mx_wgcna_filtered(), error = function(e) NULL); req(f)
      miss <- colMeans(is.na(f$mat)) * 100
      df <- data.frame(sample = names(miss), pct = miss)
      df$sample <- factor(df$sample, levels = df$sample[order(-df$pct)])
      ggplot(df, aes(x = sample, y = pct)) + geom_col(fill = ARTHOMIX_COLORS$blue) +
        labs(x = NULL, y = "Missing (%)") + theme_arthomix() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    })

    ## ==== 3. Soft Threshold ================================================

    tab_power_ui <- function() {
      f <- tryCatch(mx_wgcna_filtered(), error = function(e) NULL)
      tagList(
        div(class = "card",
            div(class = "card-title", icon("wave-square"), "Soft Threshold Selection"),
            if (is.null(f)) p(class = "empty-note", icon("circle-info"), "Complete Data & Filtering before continuing to Soft Threshold.") else tagList(
              fluidRow(
                column(4, radioButtons(ns("network_type"), "Network type",
                  choices = c("Signed" = "signed", "Signed hybrid" = "signed hybrid", "Unsigned" = "unsigned"), selected = "signed")),
                column(4, radioButtons(ns("cor_method"), "Correlation method",
                  choices = c("Biweight midcorrelation (bicor)" = "bicor", "Pearson" = "pearson"),
                  selected = if (isTRUE(methyl_dataset$preloaded)) "pearson" else "bicor")),
                column(4, radioButtons(ns("power_mode"), "Power selection", choices = c("Automatic" = "auto", "Manual" = "manual"), selected = "auto"))
              ),
              fluidRow(
                column(4, checkboxInput(ns("use_custom_powers"), "Custom power vector", value = isTRUE(methyl_dataset$preloaded))),
                column(4, numericInput(ns("r_sq_cutoff"), "Target scale-free R² cutoff", value = 0.85, min = 0, max = 1, step = 0.01)),
                column(4, conditionalPanel(condition = sprintf("input['%s'] == 'manual'", ns("power_mode")),
                  numericInput(ns("manual_power"), "Manual power", value = 6, min = 1, max = 30, step = 1)))
              ),
              conditionalPanel(condition = sprintf("input['%s']", ns("use_custom_powers")),
                textInput(ns("custom_powers"), "Powers (comma-separated)",
                  value = if (isTRUE(methyl_dataset$preloaded)) "1,2,3,4,5,6,7,8,9,10,12,14,16,18,20" else "")),
              conditionalPanel(condition = sprintf("!input['%s']", ns("use_custom_powers")),
                fluidRow(
                  column(4, numericInput(ns("power_min"), "Minimum power", value = 1, min = 1, max = 30, step = 1)),
                  column(4, numericInput(ns("power_max"), "Maximum power", value = 20, min = 2, max = 30, step = 1)),
                  column(4, numericInput(ns("power_step"), "Step", value = 1, min = 1, max = 5, step = 1))
                )),
              actionButton(ns("power_btn"), "Run Soft Threshold Analysis", icon = icon("play"), class = "btn-primary btn-sm")
            )
        ),
        withSpinner(uiOutput(ns("power_result_ui")), color = "#2563EB", type = 6)
      )
    }

    mx_wgcna_sft <- eventReactive(input$power_btn, {
      f <- tryCatch(mx_wgcna_filtered(), error = function(e) NULL)
      validate(need(!is.null(f), "Build the filtered matrix (Data & Filtering) first."))
      texpr <- t(f$mat)
      powers <- if (isTRUE(input$use_custom_powers)) {
        p <- suppressWarnings(as.numeric(trimws(strsplit(input$custom_powers %||% "", ",")[[1]])))
        p <- sort(unique(p[is.finite(p) & p > 0]))
        validate(need(length(p) >= 2, "Enter at least two valid powers, comma-separated."))
        p
      } else {
        seq(input$power_min %||% 1, input$power_max %||% 20, by = input$power_step %||% 1)
      }
      cor_fnc_name <- if (identical(input$cor_method, "bicor")) "bicor" else "cor"
      sft <- WGCNA::pickSoftThreshold(texpr, powerVector = powers, networkType = input$network_type %||% "signed",
                                       corFnc = cor_fnc_name, RsquaredCut = input$r_sq_cutoff %||% 0.85, verbose = 0)
      fi <- sft$fitIndices
      max_r_sq <- suppressWarnings(max(fi$SFT.R.sq, na.rm = TRUE))
      reached <- is.finite(max_r_sq) && max_r_sq >= (input$r_sq_cutoff %||% 0.85)
      ## Same rule the published script uses: first power reaching the
      ## cutoff, else the single best-observed fit across the tested range -
      ## never a hardcoded literal power.
      auto_power <- if (reached) min(fi$Power[fi$SFT.R.sq >= (input$r_sq_cutoff %||% 0.85)]) else fi$Power[which.max(fi$SFT.R.sq)]
      power <- if (identical(input$power_mode, "manual")) (input$manual_power %||% 6) else auto_power
      list(fit_indices = fi, power = power, auto_power = auto_power, reached_cutoff = reached,
           network_type = input$network_type %||% "signed",
           cor_method = if (identical(input$cor_method, "bicor")) "bicor" else "pearson",
           r_sq_cutoff = input$r_sq_cutoff %||% 0.85, powers = powers, max_r_sq = max_r_sq)
    })

    output$power_result_ui <- renderUI({
      res <- tryCatch(mx_wgcna_sft(), error = function(e) NULL)
      if (is.null(res)) return(NULL)
      gr <- mx_wgcna_guardrails(max_r_sq = res$max_r_sq)
      row_at_power <- res$fit_indices[res$fit_indices$Power == res$power, , drop = FALSE][1, ]
      tagList(
        div(class = "card",
            div(class = "card-title", icon("chart-line"), "Soft Threshold Results"),
            p(class = "empty-note", icon(if (gr$poor_sft_fit) "triangle-exclamation" else "check"), sprintf(
              "Recommended power: %s (%s). Scale-free fit at this power: R²=%.3f. Mean connectivity: %.1f.",
              res$power,
              if (identical(input$power_mode, "manual")) "manual override" else if (res$reached_cutoff) "first power reaching the target cutoff" else "best-observed fit - cutoff not reached by any tested power",
              row_at_power$SFT.R.sq, row_at_power$mean.k.)),
            if (gr$poor_sft_fit) p(class = "empty-note", icon("triangle-exclamation"),
              "No strong scale-free topology fit was detected across the tested powers. Methylation networks may not reach conventional transcriptomic-style fit thresholds; consider widening the power range, reducing highly redundant CpGs, or evaluating the network on connectivity/interpretability rather than forcing a threshold."),
            fluidRow(
              column(6, withSpinner(plotOutput(ns("sft_r2_plot"), height = 300), color = "#2563EB", type = 6)),
              column(6, withSpinner(plotOutput(ns("sft_k_plot"), height = 300), color = "#2563EB", type = 6))
            ),
            DT::dataTableOutput(ns("sft_table")),
            downloadButton(ns("download_sft_csv"), "Soft-threshold table (CSV)", class = "btn-default btn-sm")
        )
      )
    })
    output$sft_r2_plot <- renderPlot({
      res <- tryCatch(mx_wgcna_sft(), error = function(e) NULL); req(res)
      df <- res$fit_indices
      ggplot(df, aes(x = Power, y = -sign(slope) * SFT.R.sq)) +
        geom_point(color = ARTHOMIX_COLORS$blue, size = 2) + geom_line(color = ARTHOMIX_COLORS$blue) +
        geom_hline(yintercept = res$r_sq_cutoff, linetype = "dashed", color = ARTHOMIX_STATUS$critical) +
        geom_vline(xintercept = res$power, linetype = "dotted", color = ARTHOMIX_COLORS$ink_muted) +
        labs(x = "Soft-thresholding power", y = expression(paste("Scale free topology (", R^2, ")"))) + theme_arthomix()
    })
    output$sft_k_plot <- renderPlot({
      res <- tryCatch(mx_wgcna_sft(), error = function(e) NULL); req(res)
      df <- res$fit_indices
      ggplot(df, aes(x = Power, y = mean.k.)) +
        geom_point(color = ARTHOMIX_COLORS$blue, size = 2) + geom_line(color = ARTHOMIX_COLORS$blue) +
        geom_vline(xintercept = res$power, linetype = "dotted", color = ARTHOMIX_COLORS$ink_muted) +
        labs(x = "Soft-thresholding power", y = "Mean connectivity") + theme_arthomix()
    })
    output$sft_table <- DT::renderDataTable({
      res <- tryCatch(mx_wgcna_sft(), error = function(e) NULL); req(res)
      df <- res$fit_indices[, c("Power", "SFT.R.sq", "slope", "truncated.R.sq", "mean.k.", "median.k.", "max.k.")]
      colnames(df) <- c("Power", "SFT.R.sq", "slope", "truncated R²", "Mean connectivity", "Median connectivity", "Max connectivity")
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", pageLength = 30), class = "stripe hover compact") %>% DT::formatRound(columns = 2:7, digits = 3)
    })
    output$download_sft_csv <- downloadHandler(
      filename = function() "methylomics_wgcna_soft_threshold.csv",
      content = function(file) utils::write.csv(mx_wgcna_sft()$fit_indices, file, row.names = FALSE)
    )

    ## ==== 4. Network & Modules =============================================

    tab_modules_ui <- function() {
      sft <- tryCatch(mx_wgcna_sft(), error = function(e) NULL)
      tagList(
        div(class = "card",
            div(class = "card-title", icon("diagram-project"), "Network & Modules"),
            if (is.null(sft)) p(class = "empty-note", icon("circle-info"), "Compute the soft-threshold power first.") else tagList(
              fluidRow(
                column(4, numericInput(ns("net_power"), "Soft-threshold power", value = sft$power, min = 1, max = 30, step = 1)),
                column(4, selectInput(ns("tom_type"), "TOM type", choices = c(
                  "Signed" = "signed", "Unsigned" = "unsigned", "Signed Nowick" = "signed Nowick",
                  "Unsigned 2" = "unsigned 2", "Signed 2" = "signed 2", "Signed Nowick 2" = "signed Nowick 2"), selected = "signed")),
                column(4, numericInput(ns("max_block_size"), "Maximum block size", value = 5000, min = 500, step = 500))
              ),
              fluidRow(
                column(3, numericInput(ns("min_module_size"), "Minimum module size", value = if (isTRUE(methyl_dataset$preloaded)) 20 else 30, min = 2, step = 1)),
                column(3, selectInput(ns("deep_split"), "Deep split", choices = stats::setNames(0:4, 0:4), selected = "2")),
                column(3, numericInput(ns("merge_cut_height"), "Merge cut height", value = 0.25, min = 0, max = 1, step = 0.01)),
                column(3, checkboxInput(ns("pam_stage"), "Use PAM refinement (pamStage)", value = TRUE))
              ),
              fluidRow(
                column(3, checkboxInput(ns("pam_respects_dendro"), "PAM respects dendrogram", value = TRUE)),
                column(3, numericInput(ns("reassign_threshold"), "Module reassignment threshold", value = 1e-6, min = 0, max = 1, step = 1e-6)),
                column(3, numericInput(ns("min_kme_to_stay"), "Minimum KME to stay", value = 0.3, min = 0, max = 1, step = 0.01)),
                column(3, numericInput(ns("min_core_kme"), "Minimum core KME", value = 0.5, min = 0, max = 1, step = 0.01))
              ),
              p(class = "empty-note", icon("circle-info"), "Maximum block size controls memory use on large probe sets - the published pipeline lowered this from 20,000 to 5,000 after a single 20,000×20,000 signed TOM exceeded available RAM; raise it only if this deployment has memory to spare."),
              actionButton(ns("modules_btn"), "Run WGCNA", icon = icon("play"), class = "btn-primary")
            )
        ),
        withSpinner(uiOutput(ns("modules_result_ui")), color = "#2563EB", type = 6)
      )
    }

    mx_wgcna_net <- eventReactive(input$modules_btn, {
      f <- tryCatch(mx_wgcna_filtered(), error = function(e) NULL)
      sft <- tryCatch(mx_wgcna_sft(), error = function(e) NULL)
      validate(need(!is.null(f), "Build the filtered matrix (Data & Filtering) first."))
      validate(need(!is.null(sft), "Compute the soft-threshold power first."))
      texpr <- t(f$mat)
      cor_type <- if (identical(sft$cor_method, "bicor")) "bicor" else "pearson"
      key_parts <- list(texpr = texpr, power = input$net_power, network_type = sft$network_type, tom_type = input$tom_type,
                         cor_type = cor_type, deep_split = as.integer(input$deep_split), min_module_size = input$min_module_size,
                         merge_cut_height = input$merge_cut_height, max_block_size = input$max_block_size,
                         pam_stage = isTRUE(input$pam_stage), pam_respects_dendro = isTRUE(input$pam_respects_dendro),
                         reassign_threshold = input$reassign_threshold, min_kme_to_stay = input$min_kme_to_stay,
                         min_core_kme = input$min_core_kme)
      net <- withProgress(message = "Running WGCNA (blockwiseModules) - this can take a while on large probe sets...", value = 0.3, {
        get_or_compute_meth_wgcna_blocks(key_parts, function() {
          WGCNA::blockwiseModules(
            texpr, power = input$net_power, networkType = sft$network_type, TOMType = input$tom_type,
            corType = cor_type, deepSplit = as.integer(input$deep_split), minModuleSize = input$min_module_size,
            mergeCutHeight = input$merge_cut_height, maxBlockSize = input$max_block_size,
            pamStage = isTRUE(input$pam_stage), pamRespectsDendro = isTRUE(input$pam_respects_dendro),
            reassignThreshold = input$reassign_threshold, minKMEtoStay = input$min_kme_to_stay,
            minCoreKME = input$min_core_kme, numericLabels = FALSE, saveTOMs = FALSE, randomSeed = 1234, verbose = 0
          )
        })
      })
      module_colors <- net$colors
      names(module_colors) <- colnames(texpr)
      module_sizes <- as.data.frame(table(module = module_colors), stringsAsFactors = FALSE)
      colnames(module_sizes) <- c("module", "n_cpgs")
      list(net = net, texpr = texpr, module_colors = module_colors, MEs = net$MEs,
           dendro = net$dendrograms[[1]], block_colors = net$colors[net$blockGenes[[1]]], module_sizes = module_sizes,
           power = input$net_power, network_type = sft$network_type, tom_type = input$tom_type, cor_type = cor_type,
           min_module_size = input$min_module_size, deep_split = as.integer(input$deep_split),
           merge_cut_height = input$merge_cut_height, max_block_size = input$max_block_size,
           stratum_label = f$stratum_label)
    })

    output$modules_result_ui <- renderUI({
      res <- tryCatch(mx_wgcna_net(), error = function(e) NULL)
      if (is.null(res)) return(NULL)
      gr <- mx_wgcna_guardrails(module_colors = res$module_colors)
      n_real <- length(setdiff(unique(res$module_colors), "grey"))
      tagList(
        div(class = "card",
            div(class = "card-title", icon("diagram-project"), "Network Results"),
            p(class = "empty-note", icon("check"), sprintf(
              "%s: %d real module(s) detected (+ grey/unassigned), %d CpGs, %d samples, power %s.",
              res$stratum_label, n_real, length(res$module_colors), nrow(res$texpr), res$power)),
            if (gr$all_grey) p(class = "empty-note", icon("triangle-exclamation"), "No real modules were detected - everything fell into the grey (unassigned) module. Try a lower minimum module size, a different power, or less aggressive variability filtering."),
            if (!gr$all_grey && gr$single_module) p(class = "empty-note", icon("triangle-exclamation"), "Only one real module was detected - the network may be under-resolved at these settings."),
            fluidRow(
              column(6, withSpinner(plotOutput(ns("dendro_plot"), height = 340), color = "#2563EB", type = 6)),
              column(6, withSpinner(plotOutput(ns("module_size_plot"), height = 340), color = "#2563EB", type = 6))
            ),
            DT::dataTableOutput(ns("module_size_table")),
            downloadButton(ns("download_module_assignment"), "Module assignment (CSV)", class = "btn-default btn-sm")
        )
      )
    })
    output$dendro_plot <- renderPlot({
      res <- tryCatch(mx_wgcna_net(), error = function(e) NULL); req(res)
      ok <- tryCatch({
        WGCNA::plotDendroAndColors(res$dendro, res$block_colors, "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
        TRUE
      }, error = function(e) FALSE)
      validate(need(ok, "Could not draw the dendrogram for the current settings."))
    })
    output$module_size_plot <- renderPlot({
      res <- tryCatch(mx_wgcna_net(), error = function(e) NULL); req(res)
      df <- res$module_sizes
      ggplot(df, aes(x = stats::reorder(module, n_cpgs), y = n_cpgs, fill = module)) +
        geom_col() + coord_flip() +
        scale_fill_manual(values = stats::setNames(as.character(df$module), df$module), guide = "none") +
        labs(x = NULL, y = "CpGs") + theme_arthomix()
    })
    output$module_size_table <- DT::renderDataTable({
      res <- tryCatch(mx_wgcna_net(), error = function(e) NULL); req(res)
      DT::datatable(res$module_sizes[order(-res$module_sizes$n_cpgs), ], rownames = FALSE, options = list(dom = "t", pageLength = 30), class = "stripe hover compact")
    })
    output$download_module_assignment <- downloadHandler(
      filename = function() "methylomics_wgcna_module_assignment.csv",
      content = function(file) {
        res <- mx_wgcna_net()
        utils::write.csv(data.frame(cpg = names(res$module_colors), module = as.character(res$module_colors)), file, row.names = FALSE)
      }
    )

    ## ==== 5. Module-Trait Analysis =========================================

    tab_traits_ui <- function() {
      net <- tryCatch(mx_wgcna_net(), error = function(e) NULL)
      sheet <- methyl_dataset$sample_sheet
      tagList(
        div(class = "card",
            div(class = "card-title", icon("table-cells"), "Module-Trait Analysis"),
            if (is.null(net)) p(class = "empty-note", icon("circle-info"), "Run WGCNA (Network & Modules) before continuing.")
            else if (is.null(sheet)) p(class = "empty-note", icon("circle-info"), "No sample sheet available - module-trait correlation needs phenotype metadata.")
            else tagList(
              fluidRow(
                column(4, selectInput(ns("trait_col"), "Trait column", choices = colnames(sheet), selected = trait_col_default() %||% colnames(sheet)[1])),
                column(4, radioButtons(ns("mt_cor_method"), "Correlation method", choices = c("Pearson" = "pearson", "Spearman" = "spearman"), selected = "pearson")),
                column(4, radioButtons(ns("mt_correction"), "Multiple testing correction", choices = c("Benjamini-Hochberg (FDR)" = "BH", "Bonferroni" = "bonferroni"), selected = "BH"))
              ),
              numericInput(ns("mt_sig_thr"), "Significance threshold (FDR)", value = 0.05, min = 0, max = 1, step = 0.01),
              actionButton(ns("traits_btn"), "Correlate modules with trait", icon = icon("play"), class = "btn-primary btn-sm")
            )
        ),
        withSpinner(uiOutput(ns("traits_result_ui")), color = "#2563EB", type = 6)
      )
    }

    mx_wgcna_module_trait <- eventReactive(input$traits_btn, {
      net <- tryCatch(mx_wgcna_net(), error = function(e) NULL)
      validate(need(!is.null(net), "Run WGCNA (Network & Modules) first."))
      sheet <- methyl_dataset$sample_sheet
      validate(need(!is.null(sheet), "No sample sheet available."))
      ids <- methyl_sheet_sample_ids(sheet, rownames(net$texpr))
      common <- intersect(rownames(net$texpr), ids)
      validate(need(length(common) >= 3, "Fewer than 3 samples match between the network and the sample sheet."))
      ph <- as.data.frame(sheet)[match(common, ids), , drop = FALSE]
      MEs_all <- WGCNA::orderMEs(net$MEs[common, , drop = FALSE])
      enc <- mx_wgcna_encode_trait(ph, input$trait_col)
      validate(need(isTRUE(enc$ok), enc$reason %||% sprintf("Cannot use \"%s\" as a trait.", input$trait_col)))
      trait_vec <- enc$vec
      validate(need(sum(!is.na(trait_vec)) >= 3, "Fewer than 3 samples have a non-missing trait value."))
      corfn <- if (identical(input$mt_cor_method, "spearman")) function(x, y, use) stats::cor(x, y, use = use, method = "spearman") else WGCNA::cor
      module_cor <- corfn(as.matrix(MEs_all), trait_vec, use = "p")
      n_used <- sum(!is.na(trait_vec))
      module_p <- WGCNA::corPvalueStudent(module_cor, n_used)
      module_names <- sub("^ME", "", colnames(MEs_all))
      real <- module_names != "grey"
      fdr <- rep(NA_real_, length(module_p))
      method_p_adjust <- if (identical(input$mt_correction, "bonferroni")) "bonferroni" else "BH"
      fdr[real] <- stats::p.adjust(module_p[real], method = method_p_adjust)
      df <- data.frame(module = module_names, n_cpgs = as.integer(table(net$module_colors)[module_names]),
                        cor = as.numeric(module_cor), p_value = as.numeric(module_p), fdr = fdr, stringsAsFactors = FALSE)
      df <- df[order(df$p_value), ]
      list(table = df, trait_col = input$trait_col, trait_levels = enc$levels, cor_method = input$mt_cor_method,
           correction = method_p_adjust, sig_thr = input$mt_sig_thr %||% 0.05, MEs = MEs_all, n_used = n_used)
    })

    output$traits_result_ui <- renderUI({
      res <- tryCatch(mx_wgcna_module_trait(), error = function(e) NULL)
      if (is.null(res)) return(NULL)
      n_sig <- sum(res$table$fdr < res$sig_thr, na.rm = TRUE)
      tagList(
        div(class = "card",
            div(class = "card-title", icon("table-cells"), "Module-Trait Results"),
            p(class = "empty-note", icon("check"), sprintf(
              "%d of %d real module(s) significant at FDR<%.2f for trait \"%s\" (n=%d).",
              n_sig, sum(res$table$module != "grey"), res$sig_thr, res$trait_col, res$n_used)),
            withSpinner(plotOutput(ns("mt_heatmap_plot"), height = "auto"), color = "#2563EB", type = 6),
            DT::dataTableOutput(ns("mt_table")),
            downloadButton(ns("download_module_trait"), "Module-trait table (CSV)", class = "btn-default btn-sm")
        )
      )
    })
    ## Two-column Control/RA-style heatmap when the trait is a two-level
    ## factor (matching the published script's own moduletrait_plot(): one
    ## column per level, the "other" level's column is the exact negative
    ## correlation - a two-level indicator and its complement are perfectly
    ## anti-correlated, so cor(x, 1-vec) = -cor(x, vec) exactly, same
    ## p-value/FDR). Falls back to a single column labeled with the trait's
    ## own column name (not left blank) for a numeric trait. Label text is
    ## single-line and both its size and the plot's own height scale with
    ## the module count, so rows never overlap regardless of how many real
    ## modules were detected.
    mt_heatmap_df <- reactive({
      res <- tryCatch(mx_wgcna_module_trait(), error = function(e) NULL); req(res)
      df <- res$table[res$table$module != "grey", , drop = FALSE]
      validate(need(nrow(df) > 0, "No real (non-grey) modules to display."))
      mod_order <- df$module[order(df$cor)]
      if (!is.null(res$trait_levels) && length(res$trait_levels) == 2) {
        lo <- res$trait_levels[1]; hi <- res$trait_levels[2]
        long <- rbind(
          data.frame(module = df$module, group = lo, r = -df$cor, p = df$p_value, fdr = df$fdr, stringsAsFactors = FALSE),
          data.frame(module = df$module, group = hi, r = df$cor, p = df$p_value, fdr = df$fdr, stringsAsFactors = FALSE)
        )
        long$group <- factor(long$group, levels = c(lo, hi))
      } else {
        long <- data.frame(module = df$module, group = res$trait_col, r = df$cor, p = df$p_value, fdr = df$fdr, stringsAsFactors = FALSE)
        long$group <- factor(long$group, levels = res$trait_col)
      }
      long$module <- factor(long$module, levels = rev(mod_order))
      long$label <- sprintf("r=%.2f (p=%s)", long$r, formatC(long$p, format = "e", digits = 0))
      list(df = long, n_mod = length(mod_order))
    })
    output$mt_heatmap_plot <- renderPlot({
      d <- mt_heatmap_df()
      cex_text <- max(2.0, min(3.2, 60 / d$n_mod))
      ggplot(d$df, aes(x = group, y = module, fill = r)) +
        geom_tile(color = "white") + geom_text(aes(label = label), size = cex_text, lineheight = 0.9) +
        scale_fill_gradient2(low = ARTHOMIX_COLORS$blue, mid = "white", high = ARTHOMIX_STATUS$critical, midpoint = 0, limits = c(-1, 1)) +
        labs(x = NULL, y = "Module", fill = "r") + theme_arthomix()
    }, height = function() {
      d <- tryCatch(mt_heatmap_df(), error = function(e) NULL)
      if (is.null(d)) return(380)
      max(380, 60 + d$n_mod * 24)
    })
    output$mt_table <- DT::renderDataTable({
      res <- tryCatch(mx_wgcna_module_trait(), error = function(e) NULL); req(res)
      DT::datatable(res$table, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact") %>% DT::formatRound(columns = c("cor", "p_value", "fdr"), digits = 4)
    })
    output$download_module_trait <- downloadHandler(
      filename = function() "methylomics_wgcna_module_trait.csv",
      content = function(file) utils::write.csv(mx_wgcna_module_trait()$table, file, row.names = FALSE)
    )

    ## ==== 6. Hub CpGs =======================================================

    tab_hubs_ui <- function() {
      net <- tryCatch(mx_wgcna_net(), error = function(e) NULL)
      tagList(
        div(class = "card",
            div(class = "card-title", icon("star"), "Hub CpGs"),
            if (is.null(net)) p(class = "empty-note", icon("circle-info"), "Run WGCNA (Network & Modules) before continuing.")
            else tagList(
              fluidRow(
                column(4, selectInput(ns("hub_module"), "Module", choices = sort(setdiff(unique(net$module_colors), "grey")))),
                column(4, numericInput(ns("hub_kme_thr"), "Minimum |kME| (module membership)", value = 0.7, min = 0, max = 1, step = 0.05)),
                column(4, selectInput(ns("hub_top_n"), "Top N hub CpGs", choices = c("10", "20", "50", "100", "Custom"), selected = "20"))
              ),
              conditionalPanel(condition = sprintf("input['%s'] == 'Custom'", ns("hub_top_n")), numericInput(ns("hub_top_n_custom"), NULL, value = 20, min = 1, step = 1)),
              actionButton(ns("hubs_btn"), "Compute hub CpGs", icon = icon("play"), class = "btn-primary btn-sm")
            )
        ),
        withSpinner(uiOutput(ns("hubs_result_ui")), color = "#2563EB", type = 6)
      )
    }

    ## Intramodular connectivity needs the WHOLE network's TOM/adjacency (not
    ## just the selected module), so it's one of the most expensive single
    ## calls in this file - a plain reactive() (not part of the eventReactive
    ## body below) so it's computed once per network and reused across every
    ## module switch, instead of recomputed on every "Compute hub CpGs"
    ## click. Same optimization R/transcriptomics/mod_wgcna.R's own
    ## intramod_conn reactive already uses for the identical reason.
    mx_wgcna_intramod_conn <- reactive({
      net <- tryCatch(mx_wgcna_net(), error = function(e) NULL)
      req(net)
      ik <- tryCatch(
        WGCNA::intramodularConnectivity.fromExpr(net$texpr, colors = net$module_colors,
          corFnc = if (identical(net$cor_type, "bicor")) "bicor" else "cor", networkType = net$network_type, power = net$power),
        error = function(e) NULL)
      if (!is.null(ik)) rownames(ik) <- colnames(net$texpr)
      ik
    })

    mx_wgcna_hubs <- eventReactive(input$hubs_btn, {
      net <- tryCatch(mx_wgcna_net(), error = function(e) NULL)
      validate(need(!is.null(net), "Run WGCNA (Network & Modules) first."))
      req(input$hub_module)
      me_col <- paste0("ME", input$hub_module)
      validate(need(me_col %in% colnames(net$MEs), "Selected module not found - re-run Network & Modules."))
      cor_fnc_r <- if (identical(net$cor_type, "bicor")) WGCNA::bicor else WGCNA::cor
      kme <- as.numeric(cor_fnc_r(net$texpr, net$MEs[[me_col]], use = "p"))
      names(kme) <- colnames(net$texpr)
      module_cpgs <- names(net$module_colors)[net$module_colors == input$hub_module]
      ik <- mx_wgcna_intramod_conn()
      df <- data.frame(cpg = module_cpgs, kME = kme[module_cpgs], stringsAsFactors = FALSE)
      if (!is.null(ik)) {
        df$kWithin <- ik[module_cpgs, "kWithin"]
      } else df$kWithin <- NA_real_
      df$abs_kME <- abs(df$kME)

      anno <- tryCatch(methyl_get_norm_annotation(methyl_dataset$array_type), error = function(e) list(ok = FALSE))
      anno_ok <- isTRUE(anno$ok) && !is.null(anno$anno)
      if (anno_ok) {
        a <- anno$anno[match(df$cpg, rownames(anno$anno)), , drop = FALSE]
        for (cl in intersect(c("chr", "pos", "gene", "island_relation", "gene_region"), colnames(a))) df[[cl]] <- a[[cl]]
      }

      thr <- input$hub_kme_thr %||% 0.7
      filtered <- df[!is.na(df$abs_kME) & df$abs_kME >= thr, , drop = FALSE]
      filtered <- filtered[order(-filtered$abs_kME), ]
      top_n <- if (identical(input$hub_top_n, "Custom")) (input$hub_top_n_custom %||% 20) else as.integer(input$hub_top_n %||% 20)
      filtered <- utils::head(filtered, top_n)

      list(table = filtered, module = input$hub_module, threshold = thr, top_n = top_n, annotation_available = anno_ok)
    })

    output$hubs_result_ui <- renderUI({
      res <- tryCatch(mx_wgcna_hubs(), error = function(e) NULL)
      if (is.null(res)) return(NULL)
      tagList(
        div(class = "card",
            div(class = "card-title", icon("star"), "Hub CpGs"),
            p(class = "empty-note", icon("check"), sprintf("Top %d hub CpG(s) in module \"%s\" (|kME|≥%.2f).", nrow(res$table), res$module, res$threshold)),
            if (!res$annotation_available) p(class = "empty-note", icon("circle-info"), "CpG annotation is unavailable for this dataset - upload a compatible annotation, or use a 450K/EPIC array, to enable genomic interpretation."),
            DT::dataTableOutput(ns("hubs_table")),
            downloadButton(ns("download_hubs"), "Hub CpGs (CSV)", class = "btn-default btn-sm")
        )
      )
    })
    output$hubs_table <- DT::renderDataTable({
      res <- tryCatch(mx_wgcna_hubs(), error = function(e) NULL); req(res)
      DT::datatable(res$table, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })
    output$download_hubs <- downloadHandler(
      filename = function() sprintf("methylomics_wgcna_hub_cpgs_%s.csv", mx_wgcna_hubs()$module),
      content = function(file) utils::write.csv(mx_wgcna_hubs()$table, file, row.names = FALSE)
    )

    ## ==== 7. Results & Export ===============================================

    output$summary_table_ui <- DT::renderDataTable({
      f <- tryCatch(mx_wgcna_filtered(), error = function(e) NULL)
      req(f)
      net <- tryCatch(mx_wgcna_net(), error = function(e) NULL)
      rows <- list(
        c("Dataset source", if (isTRUE(methyl_dataset$preloaded)) "Preloaded whole-blood dataset (GSE42861)" else (methyl_dataset$source %||% "Uploaded dataset")),
        c("Sex stratum", f$stratum_label),
        c("Initial CpGs", format(f$n_probes_in, big.mark = ",")),
        c("Filtered CpGs", format(f$n_probes_kept, big.mark = ",")),
        c("Samples (post-filter)", f$n_samples_kept),
        c("Variability filter", sprintf("%s, top %s", f$var_method, f$top_n)),
        c("Residualization", f$resid_note)
      )
      if (!is.null(net)) rows <- c(rows, list(
        c("Correlation method", net$cor_type),
        c("Network type", net$network_type),
        c("Soft power", net$power),
        c("TOM type", net$tom_type),
        c("Minimum module size", net$min_module_size),
        c("Deep split", net$deep_split),
        c("Merge cut height", net$merge_cut_height),
        c("Max block size", net$max_block_size),
        c("Modules detected", length(setdiff(unique(net$module_colors), "grey")))
      ))
      df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
      colnames(df) <- c("Parameter", "Value")
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", pageLength = 30), class = "stripe hover compact")
    })

    output$download_filtered_cpgs <- downloadHandler(
      filename = function() "methylomics_wgcna_filtered_cpgs.csv",
      content = function(file) utils::write.csv(data.frame(cpg = rownames(mx_wgcna_filtered()$mat)), file, row.names = FALSE)
    )
    output$download_module_assignment2 <- downloadHandler(
      filename = function() "methylomics_wgcna_module_assignment.csv",
      content = function(file) {
        net <- mx_wgcna_net()
        utils::write.csv(data.frame(cpg = names(net$module_colors), module = as.character(net$module_colors)), file, row.names = FALSE)
      }
    )
    output$download_params <- downloadHandler(
      filename = function() "methylomics_wgcna_parameters.csv",
      content = function(file) {
        f <- mx_wgcna_filtered(); net <- tryCatch(mx_wgcna_net(), error = function(e) NULL)
        rows <- list(
          c("dataset_source", if (isTRUE(methyl_dataset$preloaded)) "preloaded" else "uploaded"),
          c("sex_stratum", f$stratum_label),
          c("max_probe_missing_pct", input$max_probe_missing),
          c("variability_method", f$var_method), c("top_n_cpgs", f$top_n),
          c("residualization_covariates", paste(f$resid_covariates, collapse = "; "))
        )
        if (!is.null(net)) rows <- c(rows, list(
          c("network_type", net$network_type), c("cor_type", net$cor_type), c("power", net$power),
          c("tom_type", net$tom_type), c("min_module_size", net$min_module_size), c("deep_split", net$deep_split),
          c("merge_cut_height", net$merge_cut_height), c("max_block_size", net$max_block_size)
        ))
        df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
        colnames(df) <- c("parameter", "value")
        utils::write.csv(df, file, row.names = FALSE)
      }
    )

    ## ---- Compare with published results (preloaded dataset only) --------

    reference_panel_ui <- function() {
      div(class = "card",
          div(class = "card-title", icon("book"), "Compare with published results"),
          p(class = "submodule-desc", "Static reference tables from the published sex-stratified analysis (script05_wgcna_sexstratified) - not this run's live output."),
          radioButtons(ns("ref_sex"), "Sex", inline = TRUE, choices = c("Female" = "female", "Male" = "male"),
                       selected = if (identical(input$sex_stratum, unname(sex_choices_r()["Male"]))) "male" else "female"),
          DT::dataTableOutput(ns("ref_module_trait_table"))
      )
    }
    output$ref_module_trait_table <- DT::renderDataTable({
      sex <- input$ref_sex %||% "female"
      mt <- load_default_wgcna_module_trait(sex)
      validate(need(!is.null(mt), "Published reference tables aren't available in this deployment."))
      DT::datatable(mt, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })

    ## ---- Optional Functional Enrichment (Fisher's exact test against the --
    ## DMP/DMR biomarker panel - the published pipeline's own convergent-
    ## evidence check, reused here rather than inventing a separate GO/KEGG
    ## path from a CpG->gene mapping)

    enrichment_ui <- function() {
      mt <- tryCatch(mx_wgcna_module_trait(), error = function(e) NULL)
      div(class = "card",
          div(class = "card-title", icon("flask"), "Functional Enrichment"),
          if (is.null(mt)) p(class = "empty-note", icon("circle-info"), "Run Module-Trait Analysis before continuing.")
          else tagList(
            p(class = "submodule-desc", "Fisher's exact test of each significant module's CpGs against the DMP/DMR biomarker panel (script04_dmr_sexstratified) - the same convergent-evidence check the published pipeline itself performs."),
            actionButton(ns("enrich_btn"), "Run Functional Enrichment", icon = icon("play"), class = "btn-primary btn-sm"),
            withSpinner(uiOutput(ns("enrich_result_ui")), color = "#2563EB", type = 6)
          )
      )
    }

    mx_wgcna_enrichment <- eventReactive(input$enrich_btn, {
      mt <- tryCatch(mx_wgcna_module_trait(), error = function(e) NULL)
      net <- tryCatch(mx_wgcna_net(), error = function(e) NULL)
      validate(need(!is.null(mt) && !is.null(net), "Run Module-Trait Analysis first."))
      sex_choice <- input$sex_stratum %||% "__all__"
      panel_sex <- if (identical(sex_choice, unname(sex_choices_r()["Male"]))) "male" else "female"
      panel <- load_default_dmr_biomarker_panel(panel_sex)
      validate(need(!is.null(panel), "The DMP/DMR biomarker panel reference isn't available in this deployment (needs the bundled data/preloaded/methylomics/ folder)."))
      sig_modules <- mt$table$module[!is.na(mt$table$fdr) & mt$table$fdr < mt$sig_thr & mt$table$module != "grey"]
      validate(need(length(sig_modules) > 0, "No real (non-grey) module reached the significance threshold on the Module-Trait Analysis tab - enrichment is not applicable."))
      background <- names(net$module_colors)
      in_panel <- background %in% panel$cpg
      out <- do.call(rbind, lapply(sig_modules, function(m) {
        in_module <- net$module_colors == m
        tab <- matrix(c(sum(in_module & in_panel), sum(in_module & !in_panel),
                         sum(!in_module & in_panel), sum(!in_module & !in_panel)), nrow = 2)
        ft <- stats::fisher.test(tab)
        data.frame(module = m, n_cpgs = sum(in_module), n_panel_in_module = sum(in_module & in_panel),
                   n_panel_total = sum(in_panel), odds_ratio = unname(ft$estimate), p_value = ft$p.value, stringsAsFactors = FALSE)
      }))
      out$fdr <- stats::p.adjust(out$p_value, method = "BH")
      out <- out[order(out$p_value), ]
      list(table = out, panel_sex = panel_sex, n_panel = nrow(panel))
    })

    output$enrich_result_ui <- renderUI({
      res <- tryCatch(mx_wgcna_enrichment(), error = function(e) NULL)
      if (is.null(res)) return(NULL)
      tagList(
        p(class = "empty-note", icon("check"), sprintf("Tested against the %s biomarker panel (%d CpGs).", res$panel_sex, res$n_panel)),
        DT::dataTableOutput(ns("enrich_table")),
        downloadButton(ns("download_enrichment"), "Enrichment table (CSV)", class = "btn-default btn-sm")
      )
    })
    output$enrich_table <- DT::renderDataTable({
      res <- tryCatch(mx_wgcna_enrichment(), error = function(e) NULL); req(res)
      DT::datatable(res$table, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact") %>% DT::formatRound(columns = c("odds_ratio", "p_value", "fdr"), digits = 4)
    })
    output$download_enrichment <- downloadHandler(
      filename = function() "methylomics_wgcna_enrichment.csv",
      content = function(file) utils::write.csv(mx_wgcna_enrichment()$table, file, row.names = FALSE)
    )

    tab_export_ui <- function() {
      f <- tryCatch(mx_wgcna_filtered(), error = function(e) NULL)
      net <- tryCatch(mx_wgcna_net(), error = function(e) NULL)
      tagList(
        div(class = "card",
            div(class = "card-title", icon("chart-simple"), "Analysis Summary"),
            if (is.null(f)) p(class = "empty-note", icon("circle-info"), "Nothing has been run yet.") else DT::dataTableOutput(ns("summary_table_ui"))
        ),
        if (!is.null(net)) div(class = "card",
          div(class = "card-title", icon("download"), "Export"),
          div(style = "display:flex; gap:8px; flex-wrap:wrap;",
              downloadButton(ns("download_filtered_cpgs"), "Filtered CpG list (CSV)", class = "btn-default btn-sm"),
              downloadButton(ns("download_module_assignment2"), "Module assignment (CSV)", class = "btn-default btn-sm"),
              downloadButton(ns("download_params"), "Analysis parameters (CSV)", class = "btn-default btn-sm")
          )
        ),
        if (isTRUE(methyl_dataset$preloaded)) reference_panel_ui(),
        if (!is.null(net)) enrichment_ui()
      )
    }

    ## ==== Top-level assembly ================================================

    main_ui <- function() {
      tagList(
        status_ui(),
        div(class = "tx-menu-wrap",
            ## Each tab body is its own uiOutput/renderUI pair (not a plain
            ## function called inline here) so that triggering one stage's
            ## eventReactive - which tabsetPanel's eager rendering would
            ## otherwise evaluate for every tab in a single pass - only
            ## invalidates/spins the ONE Shiny output that actually reads
            ## it, instead of blocking the whole tabsetPanel under one
            ## spinner every time any stage is (re)computed.
            tabsetPanel(
              id = ns("subtabs"), type = "tabs", header = tagList(tags$hr()),
              tabPanel(value = "data", title = mx_wgcna_tab_title("filter", "Data & Filtering"), br(), withSpinner(uiOutput(ns("tabout_data")), color = "#2563EB", type = 6)),
              tabPanel(value = "sampleqc", title = mx_wgcna_tab_title("magnifying-glass-chart", "Sample QC"), br(), withSpinner(uiOutput(ns("tabout_sampleqc")), color = "#2563EB", type = 6)),
              tabPanel(value = "power", title = mx_wgcna_tab_title("wave-square", "Soft Threshold"), br(), withSpinner(uiOutput(ns("tabout_power")), color = "#2563EB", type = 6)),
              tabPanel(value = "modules", title = mx_wgcna_tab_title("diagram-project", "Network & Modules"), br(), withSpinner(uiOutput(ns("tabout_modules")), color = "#2563EB", type = 6)),
              tabPanel(value = "traits", title = mx_wgcna_tab_title("table-cells", "Module-Trait Analysis"), br(), withSpinner(uiOutput(ns("tabout_traits")), color = "#2563EB", type = 6)),
              tabPanel(value = "hubs", title = mx_wgcna_tab_title("star", "Hub CpGs"), br(), withSpinner(uiOutput(ns("tabout_hubs")), color = "#2563EB", type = 6)),
              tabPanel(value = "export", title = mx_wgcna_tab_title("download", "Results & Export"), br(), withSpinner(uiOutput(ns("tabout_export")), color = "#2563EB", type = 6))
            )
        ),
        div(class = "card",
            div(class = "card-title", icon("venus-mars"), "Sex Stratum"),
            radioButtons(ns("sex_stratum"), NULL, inline = TRUE, choices = sex_choices_r(), selected = mx_default_sex()),
            if (isTRUE(methyl_dataset$preloaded) && !is.null(sex_col()))
              p(class = "empty-note", icon("circle-info"), "The published analysis is sex-stratified - pick Female or Male to reproduce it per sex. \"All samples\" runs a live network across both sexes together, which was not part of the published methodology.")
            else if (is.null(sex_col()))
              p(class = "empty-note", icon("circle-info"), "No sex/gender column detected in the sample sheet - only \"All samples\" is available.")
        )
      )
    }

    output$tabout_data <- renderUI(tab_data_ui())
    output$tabout_sampleqc <- renderUI(tab_sampleqc_ui())
    output$tabout_power <- renderUI(tab_power_ui())
    output$tabout_modules <- renderUI(tab_modules_ui())
    output$tabout_traits <- renderUI(tab_traits_ui())
    output$tabout_hubs <- renderUI(tab_hubs_ui())
    output$tabout_export <- renderUI(tab_export_ui())

    output$body_ui <- renderUI({
      if (!isTRUE(methyl_dataset$preloaded) && is.null(methyl_dataset$beta)) {
        return(div(class = "card",
          div(class = "card-title", icon("triangle-exclamation"), "No methylation dataset loaded"),
          p(class = "submodule-desc", "Load the preloaded whole-blood dataset or upload your own beta/M-value matrix on the Methylomics Dataset tab first.")
        ))
      }
      main_ui()
    })
  })
}
