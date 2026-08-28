## R/methylomics/mod_methyl_candidates.R
## Submodule: Candidate CpGs (Module-DMR Overlap).
##
## Module CpGs -> DMR coordinates -> genomic overlap -> candidate CpGs ->
## filtering/prioritization -> results. Five sub-tabs:
##   1. Data & Filters      - pick Preloaded/Upload, load, set filters.
##   2. DMR-CpG Overlap     - coordinate overlap (GenomicRanges::findOverlaps)
##                            between module CpGs and filtered DMRs.
##   3. Module-DMR Overlap  - per-module overlap counts + one-sided Fisher's
##                            exact enrichment test vs. the tested CpG universe.
##   4. Candidate CpGs      - filter/rank the subtab-2 overlap table.
##   5. Visualization       - one Generate-plot button per chart, each gated
##                            on its own upstream analysis having run.
##
## Preloaded path reads global.R's loaders (load_default_wgcna_module_
## assignment/load_default_dmr/load_default_dmp) plus ChAMPdata::probe.features
## for coordinates/gene/island/region - nothing here reruns WGCNA or DMR
## calling. Upload path accepts a module-assignment table and a DMR-results
## table (annotation table optional if coordinates are already in the module
## table); columns are auto-detected the same way mod_preprocessing.R's
## pp_guess_col() does.

## ---- Column-name detection --------------------------------------------

MCD_CPG_ID_PATTERNS   <- c("^cpg_id$", "^cpg$", "^probe_id$", "^probeid$", "^probe$", "^illumina_id$", "^id$", "cpg", "probe")
MCD_MODULE_PATTERNS   <- c("^module_color$", "^modulecolor$", "^module$", "^color$", "module", "color")
MCD_KME_PATTERNS      <- c("^kme$", "module_membership", "^mm$", "membership", "kme")
MCD_CHR_PATTERNS      <- c("^chr$", "^chromosome$", "^seqnames$", "chr")
MCD_POS_PATTERNS      <- c("^pos$", "^position$", "^mapinfo$", "^coordinate$", "^start_pos$", "pos")
MCD_DMR_ID_PATTERNS   <- c("^dmr_id$", "^dmrid$", "^dmr$", "region_id", "^id$")
MCD_FDR_PATTERNS      <- c("^dmr_fdr$", "^fdr_bacon$", "^fdr$", "^adj\\.?p\\.?val$", "^q_?value$", "^padj$", "fdr")
MCD_PVAL_PATTERNS     <- c("^dmr_p(value)?$", "^p_bacon$", "^p_?value$", "^pvalue$", "^p\\.value$", "^pval$", "^stouffer$", "^fisher$", "^p$")
MCD_DBETA_PATTERNS    <- c("^meandiff$", "^mean_diff$", "^delta_beta$", "^deltabeta$", "^dbeta$", "^effect_size$", "^logfc$", "delta", "diff")
MCD_NCPGS_PATTERNS    <- c("^no\\.cpgs$", "^n_cpgs$", "^ncpgs$", "^num_cpgs$", "^no_cpgs$", "cpgs")
MCD_GENE_PATTERNS     <- c("^overlapping\\.genes$", "^gene_symbol$", "^genesymbol$", "^gene$", "^genes$", "^nearest_gene$", "gene")
MCD_ISLAND_PATTERNS   <- c("^cgi$", "^island$", "^cpg_island$", "^relation_to_island$", "^feat\\.cgi$", "island")
MCD_FEATURE_PATTERNS  <- c("^feature$", "^genomic_region$", "^region$", "^annotation$", "feature")
MCD_DIRECTION_PATTERNS <- c("^direction$")

## First column matching any pattern (case-insensitive regex), tried in
## order; NULL if none match. Callers treat that as required (validation
## error) or optional (column omitted downstream).
mcd_find_col <- function(cols, patterns) {
  for (p in patterns) {
    hit <- cols[grepl(p, cols, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  NULL
}

mcd_standardize_module_assign <- function(df) {
  cols <- colnames(df)
  cpg_col <- mcd_find_col(cols, MCD_CPG_ID_PATTERNS)
  mod_col <- mcd_find_col(cols, MCD_MODULE_PATTERNS)
  if (is.null(cpg_col) || is.null(mod_col)) {
    return(list(ok = FALSE, reason = "Could not detect a CpG/probe ID column and a module/module-color column in the module assignment table."))
  }
  kme_col <- mcd_find_col(cols, MCD_KME_PATTERNS)
  out <- data.frame(cpg = as.character(df[[cpg_col]]), module = as.character(df[[mod_col]]), stringsAsFactors = FALSE)
  if (!is.null(kme_col)) out$kme <- suppressWarnings(as.numeric(df[[kme_col]]))
  out <- out[!is.na(out$cpg) & nzchar(out$cpg) & !is.na(out$module) & nzchar(out$module), , drop = FALSE]
  out <- out[!duplicated(out$cpg), , drop = FALSE]
  if (nrow(out) == 0) return(list(ok = FALSE, reason = "No usable rows remained in the module assignment table after removing blank CpG IDs/modules."))
  list(ok = TRUE, df = out, detected = list(cpg_id = cpg_col, module = mod_col, kme = kme_col))
}

mcd_standardize_dmr <- function(df) {
  cols <- colnames(df)
  chr_col <- mcd_find_col(cols, MCD_CHR_PATTERNS)
  start_col <- mcd_find_col(cols, c("^start$", "dmr_start", "region_start"))
  end_col <- mcd_find_col(cols, c("^end$", "dmr_end", "region_end"))
  if (is.null(chr_col) || is.null(start_col) || is.null(end_col)) {
    return(list(ok = FALSE, reason = "Could not detect chromosome/start/end columns in the DMR results table."))
  }
  id_col <- mcd_find_col(cols, MCD_DMR_ID_PATTERNS)
  fdr_col <- mcd_find_col(cols, MCD_FDR_PATTERNS)
  p_col <- mcd_find_col(cols, MCD_PVAL_PATTERNS)
  if (!is.null(fdr_col) && identical(fdr_col, p_col)) p_col <- NULL
  dbeta_col <- mcd_find_col(cols, MCD_DBETA_PATTERNS)
  ncpg_col <- mcd_find_col(cols, MCD_NCPGS_PATTERNS)
  gene_col <- mcd_find_col(cols, MCD_GENE_PATTERNS)
  dir_col <- mcd_find_col(cols, MCD_DIRECTION_PATTERNS)

  out <- data.frame(
    chr = as.character(df[[chr_col]]),
    start = suppressWarnings(as.numeric(df[[start_col]])),
    end = suppressWarnings(as.numeric(df[[end_col]])),
    stringsAsFactors = FALSE
  )
  out$dmr_id <- if (!is.null(id_col)) as.character(df[[id_col]]) else sprintf("DMR%04d", seq_len(nrow(df)))
  if (!is.null(fdr_col)) out$dmr_fdr <- suppressWarnings(as.numeric(df[[fdr_col]]))
  if (!is.null(p_col)) out$dmr_pvalue <- suppressWarnings(as.numeric(df[[p_col]]))
  if (!is.null(dbeta_col)) out$delta_beta <- suppressWarnings(as.numeric(df[[dbeta_col]]))
  if (!is.null(ncpg_col)) out$n_cpgs <- suppressWarnings(as.numeric(df[[ncpg_col]]))
  if (!is.null(gene_col)) out$gene <- as.character(df[[gene_col]])
  if (!is.null(dir_col)) {
    out$direction <- tolower(as.character(df[[dir_col]]))
  } else if (!is.null(dbeta_col)) {
    out$direction <- ifelse(out$delta_beta > 0, "hyper", ifelse(out$delta_beta < 0, "hypo", NA_character_))
  }
  out <- out[!is.na(out$chr) & nzchar(out$chr) & !is.na(out$start) & !is.na(out$end), , drop = FALSE]
  out <- out[!duplicated(out$dmr_id), , drop = FALSE]
  if (nrow(out) == 0) return(list(ok = FALSE, reason = "No usable rows remained in the DMR results table after removing rows with missing chromosome/start/end."))
  list(ok = TRUE, df = out, detected = list(chr = chr_col, start = start_col, end = end_col, dmr_id = id_col,
                                             fdr = fdr_col, pvalue = p_col, delta_beta = dbeta_col, n_cpgs = ncpg_col,
                                             gene = gene_col, direction = dir_col))
}

mcd_standardize_annotation <- function(df) {
  cols <- colnames(df)
  cpg_col <- mcd_find_col(cols, MCD_CPG_ID_PATTERNS)
  chr_col <- mcd_find_col(cols, MCD_CHR_PATTERNS)
  pos_col <- mcd_find_col(cols, MCD_POS_PATTERNS)
  if (is.null(cpg_col) || is.null(chr_col) || is.null(pos_col)) {
    return(list(ok = FALSE, reason = "Could not detect CpG ID, chromosome, and position columns in the annotation/coordinate table."))
  }
  gene_col <- mcd_find_col(cols, MCD_GENE_PATTERNS)
  island_col <- mcd_find_col(cols, MCD_ISLAND_PATTERNS)
  feature_col <- mcd_find_col(cols, MCD_FEATURE_PATTERNS)
  out <- data.frame(
    cpg = as.character(df[[cpg_col]]), chr = as.character(df[[chr_col]]),
    pos = suppressWarnings(as.numeric(df[[pos_col]])), stringsAsFactors = FALSE
  )
  if (!is.null(gene_col)) out$gene <- as.character(df[[gene_col]])
  if (!is.null(island_col)) out$island <- as.character(df[[island_col]])
  if (!is.null(feature_col)) out$feature <- as.character(df[[feature_col]])
  out <- out[!is.na(out$cpg) & nzchar(out$cpg) & !is.na(out$chr) & !is.na(out$pos), , drop = FALSE]
  out <- out[!duplicated(out$cpg), , drop = FALSE]
  if (nrow(out) == 0) return(list(ok = FALSE, reason = "No usable rows remained in the annotation/coordinate table after removing rows with missing chromosome/position."))
  list(ok = TRUE, df = out, detected = list(cpg_id = cpg_col, chr = chr_col, pos = pos_col, gene = gene_col, island = island_col, feature = feature_col))
}

## Shared by the optional CpG-level stats upload and the preloaded
## load_default_dmp("sva", sex) table.
mcd_standardize_cpg_stats <- function(df) {
  cols <- colnames(df)
  cpg_col <- mcd_find_col(cols, MCD_CPG_ID_PATTERNS)
  if (is.null(cpg_col)) return(list(ok = FALSE, reason = "Could not detect a CpG/probe ID column."))
  p_col <- mcd_find_col(cols, MCD_PVAL_PATTERNS)
  fdr_col <- mcd_find_col(cols, MCD_FDR_PATTERNS)
  if (!is.null(fdr_col) && identical(fdr_col, p_col)) p_col <- NULL
  dbeta_col <- mcd_find_col(cols, MCD_DBETA_PATTERNS)
  out <- data.frame(cpg = as.character(df[[cpg_col]]), stringsAsFactors = FALSE)
  if (!is.null(p_col)) out$p_value <- suppressWarnings(as.numeric(df[[p_col]]))
  if (!is.null(fdr_col)) out$fdr <- suppressWarnings(as.numeric(df[[fdr_col]]))
  if (!is.null(dbeta_col)) out$delta_beta <- suppressWarnings(as.numeric(df[[dbeta_col]]))
  if (ncol(out) == 1) return(list(ok = FALSE, reason = "No p-value/FDR/delta-beta column was detected."))
  out <- out[!duplicated(out$cpg), , drop = FALSE]
  list(ok = TRUE, df = out)
}

## ChAMPdata::probe.features: chr/pos plus gene/genomic-region/CpG-island
## columns not exposed by mod_methyl_dmr.R's own extraction. Cached separately.
.mcd_champ_anno_cache <- new.env(parent = emptyenv())
mcd_champ_full_annotation <- function() {
  if (!is.null(.mcd_champ_anno_cache$anno)) return(.mcd_champ_anno_cache$anno)
  if (!requireNamespace("ChAMPdata", quietly = TRUE)) return(NULL)
  e <- new.env(parent = emptyenv())
  ok <- tryCatch({ utils::data("probe.features", package = "ChAMPdata", envir = e); TRUE }, error = function(e) FALSE)
  if (!ok || is.null(e$probe.features)) return(NULL)
  pf <- e$probe.features
  out <- data.frame(
    cpg = rownames(pf), chr = paste0("chr", pf$CHR), pos = pf$MAPINFO,
    gene = pf$gene, feature = pf$feature, island = pf$cgi,
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$chr) & !is.na(out$pos), , drop = FALSE]
  .mcd_champ_anno_cache$anno <- out
  out
}

## "1"/"X"/"MT" -> "chr1"/"chrX"/"chrMT"; already-prefixed values are just
## case-normalized; anything else is left as-is.
mcd_norm_chr <- function(x) {
  x <- trimws(as.character(x))
  ifelse(grepl("^chr", x, ignore.case = TRUE), paste0("chr", sub("^chr", "", x, ignore.case = TRUE)),
         ifelse(grepl("^([0-9]{1,2}|[XYM]|MT)$", x, ignore.case = TRUE), paste0("chr", toupper(x)), x))
}

## ---- DMR filtering / genomic overlap ----------------------------------

mcd_filter_dmrs <- function(dmr, fdr_max = NULL, p_max = NULL, dbeta_min = 0, mincpgs_min = 0,
                             direction = "all", chr_restrict = character(0)) {
  keep <- rep(TRUE, nrow(dmr))
  if (!is.null(fdr_max) && "dmr_fdr" %in% names(dmr)) keep <- keep & !is.na(dmr$dmr_fdr) & dmr$dmr_fdr <= fdr_max
  if (!is.null(p_max) && "dmr_pvalue" %in% names(dmr)) keep <- keep & (is.na(dmr$dmr_pvalue) | dmr$dmr_pvalue <= p_max)
  if (isTRUE(dbeta_min > 0) && "delta_beta" %in% names(dmr)) keep <- keep & (is.na(dmr$delta_beta) | abs(dmr$delta_beta) >= dbeta_min)
  if (isTRUE(mincpgs_min > 0) && "n_cpgs" %in% names(dmr)) keep <- keep & (is.na(dmr$n_cpgs) | dmr$n_cpgs >= mincpgs_min)
  if (identical(direction, "hyper") && "direction" %in% names(dmr)) keep <- keep & !is.na(dmr$direction) & dmr$direction == "hyper"
  if (identical(direction, "hypo") && "direction" %in% names(dmr)) keep <- keep & !is.na(dmr$direction) & dmr$direction == "hypo"
  if (length(chr_restrict) > 0) keep <- keep & dmr$chr %in% chr_restrict
  dmr[keep, , drop = FALSE]
}

## Joins module_assign to annotation for coordinates, builds CpG/DMR GRanges
## (DMR optionally flanked), and returns every CpG x DMR overlap pair
## (GenomicRanges::findOverlaps) plus the CpG universe with resolvable
## coordinates.
mcd_compute_overlap <- function(module_assign, annotation, dmr_filtered, flank = 0) {
  cpg_pos <- merge(module_assign, annotation, by = "cpg")
  cpg_pos <- cpg_pos[!is.na(cpg_pos$chr) & !is.na(cpg_pos$pos), , drop = FALSE]
  if (nrow(cpg_pos) == 0 || nrow(dmr_filtered) == 0) {
    return(list(joined = cpg_pos[0, , drop = FALSE], cpg_universe = cpg_pos))
  }
  cpg_gr <- GenomicRanges::GRanges(mcd_norm_chr(cpg_pos$chr), IRanges::IRanges(cpg_pos$pos, cpg_pos$pos))
  dmr_gr <- GenomicRanges::GRanges(mcd_norm_chr(dmr_filtered$chr),
                                    IRanges::IRanges(pmax(1, dmr_filtered$start - flank), dmr_filtered$end + flank))
  hits <- GenomicRanges::findOverlaps(cpg_gr, dmr_gr)
  if (length(hits) == 0) return(list(joined = cpg_pos[0, , drop = FALSE], cpg_universe = cpg_pos))
  qh <- S4Vectors::queryHits(hits); sh <- S4Vectors::subjectHits(hits)
  left <- cpg_pos[qh, , drop = FALSE]
  right <- dmr_filtered[sh, , drop = FALSE]
  ren <- c(chr = "dmr_chr", start = "dmr_start", end = "dmr_end", delta_beta = "dmr_delta_beta",
           n_cpgs = "dmr_n_cpgs", gene = "dmr_gene", direction = "dmr_direction")
  for (old in names(ren)) if (old %in% names(right)) names(right)[names(right) == old] <- ren[[old]]
  joined <- cbind(left, right)
  if (!"gene" %in% names(joined) && "dmr_gene" %in% names(joined)) joined$gene <- joined$dmr_gene
  rownames(joined) <- NULL
  list(joined = joined, cpg_universe = cpg_pos)
}

## ---- Display helpers ----------------------------------------------------

MCD_PRETTY_MAP <- c(
  cpg = "CpG ID", module = "Module", kme = "kME", dmr_id = "DMR ID", chr = "Chr", pos = "Position",
  dmr_start = "DMR Start", dmr_end = "DMR End", dmr_fdr = "DMR FDR", dmr_pvalue = "DMR p-value",
  p_value = "CpG p-value", fdr = "CpG FDR", delta_beta = "Delta Beta", dmr_delta_beta = "DMR Delta Beta",
  direction = "Direction", direction_consistency = "Consistency",
  gene = "Gene", island = "CpG Island", feature = "Genomic Region", candidate_score = "Candidate Score"
)
mcd_pretty <- function(df) {
  cols <- intersect(names(MCD_PRETTY_MAP), colnames(df))
  out <- df[, cols, drop = FALSE]
  colnames(out) <- unname(MCD_PRETTY_MAP[cols])
  out
}

mcd_detected_list <- function(detected) {
  rows <- unlist(lapply(names(detected), function(section) {
    d <- detected[[section]]
    if (is.null(d)) return(NULL)
    vapply(names(d), function(field) sprintf("%s / %s: %s", section, field, if (is.null(d[[field]])) "(not found)" else d[[field]]),
           character(1))
  }))
  tags$ul(style = "font-size:12.5px; color:var(--color-ink-secondary); margin-bottom:0;", lapply(rows, tags$li))
}

## WGCNA module names are usually literal colours - used directly as fill
## values when all valid R colour names, else falls back to the app palette
## (e.g. uploaded "M1"/"M2" labels).
mcd_module_colors <- function(mods) {
  ok <- vapply(mods, function(m) !inherits(tryCatch(grDevices::col2rgb(m), error = function(e) e), "error"), logical(1))
  if (length(mods) > 0 && all(ok)) return(stats::setNames(as.character(mods), as.character(mods)))
  pal <- c(ARTHOMIX_COLORS$blue, ARTHOMIX_COLORS$orange, ARTHOMIX_COLORS$aqua, ARTHOMIX_COLORS$violet,
           ARTHOMIX_COLORS$magenta, ARTHOMIX_COLORS$yellow, ARTHOMIX_COLORS$red)
  stats::setNames(rep(pal, length.out = length(mods)), as.character(mods))
}

## ---- Plots ---------------------------------------------------------------

mcd_plot_module_bar <- function(mod_tab) {
  d <- mod_tab[order(-mod_tab$n_overlap_cpgs), , drop = FALSE]
  d$module <- factor(d$module, levels = rev(d$module))
  ggplot(d, aes(x = n_overlap_cpgs, y = module, fill = module)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = mcd_module_colors(levels(d$module))) +
    labs(x = "Overlapping CpGs", y = NULL) + theme_arthomix()
}

mcd_plot_enrichment_heatmap <- function(mod_tab) {
  d <- mod_tab
  d$neglog10p <- -log10(pmax(d$p_value, .Machine$double.xmin))
  d$module <- factor(d$module, levels = d$module[order(-d$neglog10p)])
  ggplot(d, aes(x = module, y = "Module-DMR enrichment", fill = neglog10p)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "#EAF3FB", high = ARTHOMIX_COLORS$red, name = expression(-log[10](p))) +
    labs(x = NULL, y = NULL) + theme_arthomix() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

mcd_plot_dmr_bar <- function(overlap_tab, n = 20) {
  agg <- stats::aggregate(cpg ~ dmr_id, overlap_tab, function(x) length(unique(x)))
  colnames(agg) <- c("dmr_id", "n_cpgs")
  agg <- utils::head(agg[order(-agg$n_cpgs), , drop = FALSE], n)
  agg$dmr_id <- factor(agg$dmr_id, levels = rev(agg$dmr_id))
  ggplot(agg, aes(x = n_cpgs, y = dmr_id)) +
    geom_col(fill = ARTHOMIX_COLORS$blue) +
    labs(x = "Candidate CpGs", y = "DMR") + theme_arthomix()
}

mcd_plot_direction <- function(tab) {
  validate(need("direction" %in% names(tab), "No direction information is available to plot."))
  d <- tab[!is.na(tab$direction), , drop = FALSE]
  validate(need(nrow(d) > 0, "No direction information is available to plot."))
  agg <- as.data.frame(table(direction = d$direction))
  ggplot(agg, aes(x = direction, y = Freq, fill = direction)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = c(hyper = ARTHOMIX_STATUS$critical, hypo = ARTHOMIX_COLORS$blue)) +
    labs(x = NULL, y = "Candidate CpGs") + theme_arthomix()
}

mcd_plot_annotation_dist <- function(tab) {
  parts <- list()
  if ("feature" %in% names(tab)) {
    t1 <- as.data.frame(table(value = tab$feature)); t1$type <- "Genomic region"; parts[[length(parts) + 1]] <- t1
  }
  if ("island" %in% names(tab)) {
    t2 <- as.data.frame(table(value = tab$island)); t2$type <- "CpG island context"; parts[[length(parts) + 1]] <- t2
  }
  validate(need(length(parts) > 0, "No genomic-region or CpG-island annotation is available to plot."))
  d <- do.call(rbind, parts)
  ggplot(d, aes(x = value, y = Freq, fill = type)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~type, scales = "free", ncol = 1) +
    labs(x = NULL, y = "Candidate CpGs") + theme_arthomix() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

mcd_plot_candidate_volcano <- function(tab, highlight_cpgs = NULL) {
  dbeta_col <- if ("delta_beta" %in% names(tab)) "delta_beta" else if ("dmr_delta_beta" %in% names(tab)) "dmr_delta_beta" else NA
  sig_col <- if ("fdr" %in% names(tab)) "fdr" else if ("dmr_fdr" %in% names(tab)) "dmr_fdr" else NA
  validate(need(!is.na(dbeta_col) && !is.na(sig_col), "Neither a Delta Beta nor an FDR/p-value column is available for a volcano plot."))
  d <- tab
  d$x <- d[[dbeta_col]]
  d$sig <- -log10(pmax(d[[sig_col]], .Machine$double.xmin))
  d$is_candidate <- if (!is.null(highlight_cpgs)) d$cpg %in% highlight_cpgs else FALSE
  ggplot(d, aes(x = x, y = sig, color = is_candidate)) +
    geom_point(alpha = 0.5, size = 1.2) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = ARTHOMIX_STATUS$critical),
                        guide = if (!is.null(highlight_cpgs)) "legend" else "none",
                        labels = c("FALSE" = "Overlap CpG", "TRUE" = "Prioritized candidate")) +
    labs(x = sprintf("Delta Beta (%s)", if (identical(dbeta_col, "delta_beta")) "CpG-level" else "DMR-level"),
         y = sprintf("-log10(%s)", if (identical(sig_col, "fdr")) "CpG FDR" else "DMR FDR"), color = NULL) +
    theme_arthomix()
}

mcd_viz_block <- function(ns, id_prefix, title, icon_name, btn_label, desc) {
  div(class = "card",
      div(class = "card-title", icon(icon_name), title),
      p(class = "submodule-desc", desc),
      actionButton(ns(paste0(id_prefix, "_btn")), btn_label, icon = icon("chart-column"), class = "btn-primary btn-sm"),
      uiOutput(ns(paste0(id_prefix, "_ui")))
  )
}

## ---- Config / UI -----------------------------------------------------

mod_methyl_candidates_config <- list(
  id = "candidates", title = "Candidate CpGs (Module-DMR Overlap)", icon = "star", group = "Network",
  description = "Finds CpGs that overlap significant DMRs and ranks them by WGCNA module enrichment."
)

mod_methyl_candidates_ui <- function(id) {
  ns <- NS(id)
  div(
    class = "tx-menu-wrap",
    tabsetPanel(
      id = ns("cd_subtabs"), type = "tabs",
      tabPanel("Data & Filters", br(), uiOutput(ns("data_source_ui"))),
      tabPanel("DMR-CpG Overlap", br(), uiOutput(ns("overlap_tab_ui"))),
      tabPanel("Module-DMR Overlap", br(), uiOutput(ns("modoverlap_tab_ui"))),
      tabPanel("Candidate CpGs", br(), uiOutput(ns("candidates_tab_ui"))),
      tabPanel("Visualization", br(), uiOutput(ns("viz_tab_ui")))
    )
  )
}

## ---- Server ------------------------------------------------------------

mod_methyl_candidates_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## ===================== 1. Data & Filters ============================

    has_loaded <- reactiveVal(FALSE)
    observeEvent(input$load_btn, has_loaded(TRUE), ignoreInit = TRUE)
    observeEvent(input$data_source, has_loaded(FALSE), ignoreInit = TRUE)

    loaded <- eventReactive(input$load_btn, {
      if (identical(input$data_source, "upload")) {
        validate(need(!is.null(input$up_module_file), "Upload a CpG / WGCNA module assignment table before clicking \"Validate & Load Data\"."))
        validate(need(!is.null(input$up_dmr_file), "Upload a DMR results table before clicking \"Validate & Load Data\"."))

        ma_raw <- tryCatch(as.data.frame(data.table::fread(input$up_module_file$datapath, showProgress = FALSE)), error = function(e) NULL)
        validate(need(!is.null(ma_raw) && nrow(ma_raw) > 0, "Could not parse the module assignment file as a delimited table (CSV/TSV)."))
        dmr_raw <- tryCatch(as.data.frame(data.table::fread(input$up_dmr_file$datapath, showProgress = FALSE)), error = function(e) NULL)
        validate(need(!is.null(dmr_raw) && nrow(dmr_raw) > 0, "Could not parse the DMR results file as a delimited table (CSV/TSV)."))

        ma <- mcd_standardize_module_assign(ma_raw)
        validate(need(isTRUE(ma$ok), ma$reason))
        dmr <- mcd_standardize_dmr(dmr_raw)
        validate(need(isTRUE(dmr$ok), dmr$reason))

        annot <- NULL; annot_note <- NULL
        if (!is.null(input$up_annot_file)) {
          annot_raw <- tryCatch(as.data.frame(data.table::fread(input$up_annot_file$datapath, showProgress = FALSE)), error = function(e) NULL)
          if (!is.null(annot_raw)) {
            a <- mcd_standardize_annotation(annot_raw)
            if (isTRUE(a$ok)) annot <- a$df else annot_note <- a$reason
          } else {
            annot_note <- "Could not parse the annotation/coordinate file as a delimited table."
          }
        }
        if (is.null(annot)) {
          ma_pos <- mcd_standardize_annotation(ma_raw)
          if (isTRUE(ma_pos$ok)) { annot <- ma_pos$df; annot_note <- NULL }
        }
        if (is.null(annot)) {
          annot_note <- annot_note %||% "No CpG chromosome/position information was found - upload a CpG annotation/coordinate file, or include chromosome/position columns in the module assignment file. Genomic overlap cannot run without it."
        }

        cs <- mcd_standardize_cpg_stats(ma_raw)
        cpg_stats <- if (isTRUE(cs$ok)) cs$df else NULL

        list(source = "upload", sex = NA_character_,
             module_assign = ma$df, dmr = dmr$df, annotation = annot, cpg_stats = cpg_stats,
             detected = list(module_assignment = ma$detected, dmr_results = dmr$detected),
             notes = c(annot_note,
                       if (is.null(cpg_stats)) "No CpG-level p-value/FDR/Delta-Beta columns were detected in the module assignment file - CpG-level significance/effect-size filters and columns will be omitted."))
      } else {
        validate(need(METH_DATA_AVAILABLE, "The preloaded methylomics results folder is not available in this deployment - use \"Upload Data\" instead."))
        sex <- input$pre_sex %||% "female"
        ma_raw <- load_default_wgcna_module_assignment(sex, merged = TRUE)
        dmr_raw <- load_default_dmr(sex)
        validate(need(!is.null(ma_raw), "No preloaded WGCNA module assignment table is available for this stratum in this deployment."))
        validate(need(!is.null(dmr_raw), "No preloaded DMR results table is available for this stratum in this deployment."))

        ma <- mcd_standardize_module_assign(ma_raw)
        validate(need(isTRUE(ma$ok), ma$reason %||% "Could not read the preloaded module assignment table."))
        dmr <- mcd_standardize_dmr(dmr_raw)
        validate(need(isTRUE(dmr$ok), dmr$reason %||% "Could not read the preloaded DMR results table."))

        annot <- tryCatch(mcd_champ_full_annotation(), error = function(e) NULL)
        annot_note <- if (is.null(annot)) "Genomic coordinate/annotation data (ChAMPdata::probe.features) is not available in this deployment - overlap and annotation-based filters cannot run." else NULL

        cpg_stats_raw <- load_default_dmp("sva", sex)
        cpg_stats <- NULL
        if (!is.null(cpg_stats_raw)) {
          cs <- mcd_standardize_cpg_stats(cpg_stats_raw)
          if (isTRUE(cs$ok)) cpg_stats <- cs$df
        }

        list(source = "preloaded", sex = sex,
             module_assign = ma$df, dmr = dmr$df, annotation = annot, cpg_stats = cpg_stats,
             detected = list(module_assignment = ma$detected, dmr_results = dmr$detected),
             notes = c(annot_note,
                       if (is.null(cpg_stats)) sprintf("No preloaded CpG-level statistics are available for the %s stratum - CpG-level significance/effect-size filters and columns will be omitted.", sex)))
      }
    }, ignoreInit = TRUE)

    output$data_source_ui <- renderUI({
      tagList(
        div(class = "card",
            div(class = "card-title", icon("database"), "Data source"),
            radioButtons(ns("data_source"), NULL, inline = TRUE,
                         choiceNames = list(tagList(icon("database"), " Use Preloaded Data"), tagList(icon("upload"), " Upload Data")),
                         choiceValues = list("preloaded", "upload"), selected = "preloaded"),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'preloaded'", ns("data_source")),
              p(class = "submodule-desc", "Reproduces the sex-stratified WGCNA module assignments and DMR results already computed by this app's preloaded methylomics pipeline (script05_wgcna_sexstratified / script04_dmr_sexstratified) - nothing here reruns WGCNA or DMR calling."),
              radioButtons(ns("pre_sex"), "Sex / stratum", inline = TRUE, choices = c("Female" = "female", "Male" = "male"), selected = "female"),
              actionButton(ns("load_btn"), "Load Preloaded Data", icon = icon("play"), class = "btn-primary btn-sm")
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'upload'", ns("data_source")),
              p(class = "submodule-desc", "Upload your own module-assignment and DMR-results tables (from any WGCNA/DMR pipeline). Column names are auto-detected - exact matches aren't required."),
              fileInput(ns("up_module_file"), "CpG / WGCNA module assignment table (required)", accept = c(".csv", ".tsv", ".txt")),
              fileInput(ns("up_dmr_file"), "DMR results table (required)", accept = c(".csv", ".tsv", ".txt")),
              fileInput(ns("up_annot_file"), "CpG annotation / coordinate table (chromosome + position; optional if already in the module table)", accept = c(".csv", ".tsv", ".txt")),
              actionButton(ns("load_btn"), "Validate & Load Data", icon = icon("play"), class = "btn-primary btn-sm")
            )
        ),
        withSpinner(uiOutput(ns("load_summary_ui")), color = "#2563EB", type = 6),
        uiOutput(ns("filters_ui"))
      )
    })

    output$load_summary_ui <- renderUI({
      if (!has_loaded()) return(NULL)
      d <- loaded()
      tagList(
        div(class = "empty-note", icon("circle-check"),
            sprintf("Loaded: %s CpGs across %s modules, %s DMRs.%s",
                    format(nrow(d$module_assign), big.mark = ","), length(unique(d$module_assign$module)),
                    format(nrow(d$dmr), big.mark = ","),
                    if (identical(d$source, "preloaded")) sprintf(" Stratum: %s.", tools::toTitleCase(d$sex)) else "")),
        if (length(d$notes) > 0) div(class = "empty-note", style = "border-left-color:#fab219;", icon("triangle-exclamation"),
                                      tagList(lapply(d$notes, function(n) p(style = "margin:2px 0;", n)))),
        tags$details(tags$summary("Detected column mapping"), mcd_detected_list(d$detected))
      )
    })

    output$filters_ui <- renderUI({
      if (!has_loaded()) return(NULL)
      d <- loaded()
      mods <- sort(unique(d$module_assign$module))
      chr_choices <- sort(unique(d$dmr$chr))
      has_kme <- "kme" %in% colnames(d$module_assign)
      has_dmr_fdr <- "dmr_fdr" %in% colnames(d$dmr)
      has_dmr_p <- "dmr_pvalue" %in% colnames(d$dmr)
      has_dmr_dbeta <- "delta_beta" %in% colnames(d$dmr)
      has_dmr_ncpgs <- "n_cpgs" %in% colnames(d$dmr)
      has_dmr_dir <- has_dmr_dbeta || "direction" %in% colnames(d$dmr)
      has_cpg_dbeta <- !is.null(d$cpg_stats) && "delta_beta" %in% colnames(d$cpg_stats)
      consistency_available <- has_cpg_dbeta && has_dmr_dir

      tagList(
        div(class = "card",
            div(class = "card-title", icon("layer-group"), "Module selection"),
            radioButtons(ns("f_module_mode"), NULL, inline = TRUE,
                         choices = c("All modules" = "all", "Select specific module(s)" = "specific"), selected = "all"),
            conditionalPanel(condition = sprintf("input['%s'] == 'specific'", ns("f_module_mode")),
                             selectizeInput(ns("f_module"), "Module(s)", choices = mods, multiple = TRUE)),
            fluidRow(
              column(6, numericInput(ns("f_module_minsize"), "Module size threshold (min CpGs)", value = 0, min = 0, step = 1)),
              column(6, if (has_kme) numericInput(ns("f_module_min_kme"), "Minimum |module membership (kME)|", value = 0, min = 0, max = 1, step = 0.05)
                     else div(class = "empty-note", icon("circle-info"), "Module membership (kME) not available in this data."))
            ),
            checkboxInput(ns("f_exclude_grey"), "Exclude the grey (unassigned) module", value = TRUE)
        ),
        div(class = "card",
            div(class = "card-title", icon("map-location-dot"), "DMR significance & effect size"),
            fluidRow(
              column(6, if (has_dmr_fdr) numericInput(ns("f_dmr_fdr"), "Maximum DMR FDR / adjusted p-value", value = 0.05, min = 0, max = 1, step = 0.01)
                     else div(class = "empty-note", icon("circle-info"), "No DMR FDR column detected.")),
              column(6, if (has_dmr_p) numericInput(ns("f_dmr_p"), "Maximum DMR raw p-value", value = 1, min = 0, max = 1, step = 0.01))
            ),
            fluidRow(
              column(6, if (has_dmr_dbeta) numericInput(ns("f_dmr_dbeta"), "Minimum absolute DMR Delta-Beta", value = 0, min = 0, max = 1, step = 0.01)),
              column(6, if (has_dmr_ncpgs) numericInput(ns("f_dmr_mincpgs"), "Minimum CpGs per DMR", value = 0, min = 0, step = 1))
            ),
            if (has_dmr_dir) radioButtons(ns("f_direction"), "DMR direction", inline = TRUE,
                                           choices = c("All" = "all", "Hypermethylated" = "hyper", "Hypomethylated" = "hypo"), selected = "all")
        ),
        div(class = "card",
            div(class = "card-title", icon("sliders"), "Genomic overlap"),
            radioButtons(ns("f_overlap_mode"), "Overlap definition",
                         choices = c("CpG strictly inside the DMR" = "inside", "CpG within DMR +/- flanking distance" = "flank"), selected = "inside"),
            conditionalPanel(condition = sprintf("input['%s'] == 'flank'", ns("f_overlap_mode")),
                             numericInput(ns("f_flank_bp"), "Flanking distance (bp)", value = 0, min = 0, step = 100)),
            if (length(chr_choices) > 0) selectizeInput(ns("f_chr"), "Restrict to chromosome(s) (optional)", choices = chr_choices, multiple = TRUE),
            if (consistency_available) radioButtons(ns("f_consistency"), "Direction consistency", inline = TRUE,
                                                      choices = c("All" = "all", "Consistent only" = "consistent", "Inconsistent only" = "inconsistent"), selected = "all")
            else div(class = "empty-note", icon("circle-info"), "CpG-level and DMR-level methylation direction are not both available - direction-consistency filtering is disabled.")
        ),
        div(class = "empty-note", style = "display:flex; align-items:center; justify-content:space-between; gap:12px; flex-wrap:wrap;",
            tagList(icon("circle-info"), "These filters take effect once you run an analysis - nothing here computes anything by itself."),
            actionButton(ns("goto_overlap_btn"), "Continue to DMR-CpG Overlap", icon = icon("arrow-right"), class = "btn-primary btn-sm"))
      )
    })

    observeEvent(input$goto_overlap_btn, {
      updateTabsetPanel(session, "cd_subtabs", selected = "DMR-CpG Overlap")
    })

    ## ===================== 2. DMR-CpG Overlap ============================

    overlap_run <- eventReactive(input$overlap_run_btn, {
      d <- loaded(); req(d)
      validate(need(!is.null(d$annotation), "Genomic coordinates are required for CpG-DMR overlap, but none were detected in the loaded data. Provide a CpG annotation/coordinate file with chromosome and position columns."))

      dmr_f <- mcd_filter_dmrs(d$dmr, fdr_max = input$f_dmr_fdr, p_max = input$f_dmr_p,
                                dbeta_min = input$f_dmr_dbeta %||% 0, mincpgs_min = input$f_dmr_mincpgs %||% 0,
                                direction = input$f_direction %||% "all", chr_restrict = input$f_chr %||% character(0))
      validate(need(nrow(dmr_f) > 0, "No significant DMRs remain after the current DMR filters. Try increasing the FDR threshold or lowering the minimum |Delta-Beta|."))

      ma <- d$module_assign
      if (isTRUE(input$f_exclude_grey)) ma <- ma[!tolower(ma$module) %in% c("grey", "gray"), , drop = FALSE]
      if (identical(input$f_module_mode, "specific") && length(input$f_module) > 0) ma <- ma[ma$module %in% input$f_module, , drop = FALSE]
      if (isTRUE((input$f_module_minsize %||% 0) > 0)) {
        sizes <- table(ma$module)
        ma <- ma[ma$module %in% names(sizes)[sizes >= input$f_module_minsize], , drop = FALSE]
      }
      if (!is.null(ma$kme) && isTRUE((input$f_module_min_kme %||% 0) > 0)) {
        ma <- ma[!is.na(ma$kme) & abs(ma$kme) >= input$f_module_min_kme, , drop = FALSE]
      }
      validate(need(nrow(ma) > 0, "No CpGs remain in the selected module(s) after the module filters."))

      flank <- if (identical(input$f_overlap_mode, "flank")) (input$f_flank_bp %||% 0) else 0
      ov <- mcd_compute_overlap(ma, d$annotation, dmr_f, flank)
      out <- ov$joined
      validate(need(nrow(out) > 0, "No CpGs overlap the selected DMRs under the current filters. Try increasing the DMR FDR threshold, adding flanking distance, or changing the selected module."))

      if (!is.null(d$cpg_stats)) out <- merge(out, d$cpg_stats, by = "cpg", all.x = TRUE)

      if ("delta_beta" %in% names(out)) {
        out$direction <- ifelse(out$delta_beta > 0, "hyper", ifelse(out$delta_beta < 0, "hypo", NA_character_))
      } else if ("dmr_direction" %in% names(out)) {
        out$direction <- out$dmr_direction
      }
      if ("delta_beta" %in% names(out) && "dmr_direction" %in% names(out)) {
        out$direction_consistency <- ifelse(is.na(out$direction) | is.na(out$dmr_direction), NA_character_,
                                             ifelse(out$direction == out$dmr_direction, "Consistent", "Inconsistent"))
        cf <- input$f_consistency %||% "all"
        if (identical(cf, "consistent")) out <- out[!is.na(out$direction_consistency) & out$direction_consistency == "Consistent", , drop = FALSE]
        if (identical(cf, "inconsistent")) out <- out[!is.na(out$direction_consistency) & out$direction_consistency == "Inconsistent", , drop = FALSE]
        validate(need(nrow(out) > 0, "No CpG-DMR pairs remain after the direction-consistency filter."))
      }

      list(table = out, flank = flank,
           n_dmr_tested = nrow(d$dmr), n_dmr_passing = nrow(dmr_f), n_cpg_tested = nrow(ma),
           n_overlap_cpgs = length(unique(out$cpg)), n_overlap_dmrs = length(unique(out$dmr_id)))
    }, ignoreInit = TRUE)

    ov_has_run <- reactiveVal(FALSE)
    observeEvent(input$overlap_run_btn, ov_has_run(TRUE), ignoreInit = TRUE)
    observeEvent(input$load_btn, ov_has_run(FALSE), ignoreInit = TRUE)

    output$overlap_tab_ui <- renderUI({
      if (!has_loaded()) return(div(class = "empty-note", icon("circle-info"), "Select and load a data source on the \"Data & Filters\" tab first."))
      tagList(
        div(class = "card",
            div(class = "card-title", icon("dna"), "DMR-CpG Overlap"),
            p(class = "submodule-desc", "Coordinate-based overlap (GenomicRanges::findOverlaps): assigns each module CpG to the DMR(s) it falls inside, using the module and DMR filters set on the Data & Filters tab."),
            actionButton(ns("overlap_run_btn"), "Run DMR-CpG Overlap", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        uiOutput(ns("overlap_results_ui"))
      )
    })

    output$overlap_results_ui <- renderUI({
      if (!ov_has_run()) return(NULL)
      r <- overlap_run()
      tagList(
        div(class = "methyl-stats-row",
            fluidRow(
              valueBox(format(r$n_dmr_tested, big.mark = ","), "DMRs in table", icon = icon("map-location-dot"), color = "light-blue", width = 3),
              valueBox(format(r$n_dmr_passing, big.mark = ","), "DMRs passing filters", icon = icon("filter"), color = "purple", width = 3),
              valueBox(format(r$n_cpg_tested, big.mark = ","), "Module CpGs tested", icon = icon("dna"), color = "light-blue", width = 3),
              valueBox(format(r$n_overlap_cpgs, big.mark = ","), "Overlapping CpGs", icon = icon("star"), color = if (r$n_overlap_cpgs > 0) "green" else "light-blue", width = 3)
            )
        ),
        div(class = "card",
            div(class = "card-title", icon("table"), "DMR-CpG overlap table"),
            div(class = "table-toolbar", downloadButton(ns("overlap_download"), "Overlap table (CSV)", class = "btn-sm")),
            DT::dataTableOutput(ns("overlap_table"))
        )
      )
    })

    output$overlap_table <- DT::renderDataTable({
      req(ov_has_run())
      df <- mcd_pretty(overlap_run()$table)
      numcols <- intersect(c("DMR FDR", "DMR p-value", "CpG p-value", "CpG FDR", "Delta Beta", "DMR Delta Beta"), colnames(df))
      dt <- DT::datatable(df, rownames = FALSE, filter = "top", options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact")
      if (length(numcols) > 0) dt <- DT::formatSignif(dt, columns = numcols, digits = 4)
      dt
    })
    outputOptions(output, "overlap_table", suspendWhenHidden = FALSE)

    output$overlap_download <- downloadHandler(
      filename = function() "candidate_cpgs_dmr_overlap.csv",
      content = function(file) utils::write.csv(overlap_run()$table, file, row.names = FALSE)
    )

    ## ===================== 3. Module-DMR Overlap =========================
    ## Independent of the module filter above (every module is scored, not
    ## just the currently-selected one) so modules can be compared - but
    ## reuses the SAME DMR filters/overlap definition as subtab 2.

    module_overlap_run <- eventReactive(input$modoverlap_run_btn, {
      d <- loaded(); req(d)
      validate(need(!is.null(d$annotation), "Genomic coordinates are required for this analysis - none were detected."))

      dmr_f <- mcd_filter_dmrs(d$dmr, fdr_max = input$f_dmr_fdr, p_max = input$f_dmr_p,
                                dbeta_min = input$f_dmr_dbeta %||% 0, mincpgs_min = input$f_dmr_mincpgs %||% 0,
                                direction = input$f_direction %||% "all", chr_restrict = input$f_chr %||% character(0))
      validate(need(nrow(dmr_f) > 0, "No significant DMRs remain after the current DMR filters."))

      ma_all <- d$module_assign
      if (isTRUE(input$f_exclude_grey)) ma_all <- ma_all[!tolower(ma_all$module) %in% c("grey", "gray"), , drop = FALSE]
      validate(need(nrow(ma_all) > 0, "No CpGs remain after excluding the grey/unassigned module."))

      flank <- if (identical(input$f_overlap_mode, "flank")) (input$f_flank_bp %||% 0) else 0
      ov <- mcd_compute_overlap(ma_all, d$annotation, dmr_f, flank)
      universe <- ov$cpg_universe
      joined <- ov$joined
      validate(need(nrow(universe) > 0, "None of the module-assigned CpGs have resolvable genomic coordinates."))

      total_n <- length(unique(universe$cpg))
      total_overlap <- length(unique(joined$cpg))
      mods <- sort(unique(universe$module))

      rows <- lapply(mods, function(m) {
        mod_cpgs <- unique(universe$cpg[universe$module == m])
        n_mod <- length(mod_cpgs)
        overlap_cpgs <- unique(joined$cpg[joined$module == m])
        a <- length(overlap_cpgs); b <- n_mod - a
        c_ <- total_overlap - a; dd <- total_n - a - b - c_
        ft <- if (a >= 0 && b >= 0 && c_ >= 0 && dd >= 0 && total_n > 0) {
          tryCatch(stats::fisher.test(matrix(c(a, b, c_, dd), nrow = 2, byrow = TRUE), alternative = "greater"), error = function(e) NULL)
        } else NULL
        n_dmrs <- length(unique(joined$dmr_id[joined$module == m]))
        data.frame(module = m, n_module_cpgs = n_mod, n_overlap_cpgs = a, n_dmrs = n_dmrs,
                   pct_overlap = if (n_mod > 0) round(100 * a / n_mod, 2) else NA_real_,
                   odds_ratio = if (!is.null(ft)) round(unname(ft$estimate), 3) else NA_real_,
                   p_value = if (!is.null(ft)) ft$p.value else NA_real_, stringsAsFactors = FALSE)
      })
      res <- do.call(rbind, rows)
      res$fdr <- stats::p.adjust(res$p_value, method = "BH")
      res <- res[order(res$p_value), , drop = FALSE]
      list(table = res, total_universe = total_n, total_overlap = total_overlap, n_modules = nrow(res))
    }, ignoreInit = TRUE)

    modov_has_run <- reactiveVal(FALSE)
    observeEvent(input$modoverlap_run_btn, modov_has_run(TRUE), ignoreInit = TRUE)
    observeEvent(input$load_btn, modov_has_run(FALSE), ignoreInit = TRUE)

    output$modoverlap_tab_ui <- renderUI({
      if (!has_loaded()) return(div(class = "empty-note", icon("circle-info"), "Load a data source on the \"Data & Filters\" tab first."))
      tagList(
        div(class = "card",
            div(class = "card-title", icon("layer-group"), "Module-DMR Overlap"),
            p(class = "submodule-desc", "For every module (using the same DMR filters as the overlap tab, but every module rather than only a selected one): the number of CpGs overlapping the filtered DMRs, and a one-sided Fisher's exact test for enrichment against the tested CpG universe."),
            actionButton(ns("modoverlap_run_btn"), "Run Module-DMR Overlap", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        uiOutput(ns("modoverlap_results_ui"))
      )
    })

    output$modoverlap_results_ui <- renderUI({
      if (!modov_has_run()) return(NULL)
      r <- module_overlap_run()
      tagList(
        p(class = "empty-note", icon("circle-info"),
          sprintf("%d modules tested against a background of %d CpGs with resolvable coordinates; %d CpGs overlap the filtered DMRs overall. One-sided Fisher's exact test (enrichment).",
                  r$n_modules, r$total_universe, r$total_overlap)),
        div(class = "card",
            div(class = "card-title", icon("table"), "Module-DMR overlap statistics"),
            div(class = "table-toolbar", downloadButton(ns("modoverlap_download"), "Module statistics (CSV)", class = "btn-sm")),
            DT::dataTableOutput(ns("modoverlap_table"))
        )
      )
    })

    output$modoverlap_table <- DT::renderDataTable({
      req(modov_has_run())
      df <- module_overlap_run()$table
      colnames(df) <- c("Module", "Module CpGs", "Overlapping CpGs", "DMRs", "% Overlap", "Odds Ratio", "P-value", "FDR")
      DT::datatable(df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15), class = "stripe hover compact") %>%
        DT::formatSignif(columns = c("Odds Ratio", "P-value", "FDR"), digits = 4)
    })
    outputOptions(output, "modoverlap_table", suspendWhenHidden = FALSE)

    output$modoverlap_download <- downloadHandler(
      filename = function() "module_dmr_overlap_statistics.csv",
      content = function(file) utils::write.csv(module_overlap_run()$table, file, row.names = FALSE)
    )

    ## ===================== 4. Candidate CpGs =============================
    ## Filters/prioritizes the overlap table already computed in subtab 2 -
    ## no recomputation of the overlap itself, only narrowing + scoring it.

    output$candidates_tab_ui <- renderUI({
      if (!has_loaded()) return(div(class = "empty-note", icon("circle-info"), "Load a data source on the \"Data & Filters\" tab first."))
      if (!ov_has_run()) return(div(class = "empty-note", icon("circle-info"), "Run \"DMR-CpG Overlap\" first - Candidate CpGs filters/prioritizes that overlap table."))
      df <- overlap_run()$table
      mods <- sort(unique(df$module))
      has_dmr_fdr <- "dmr_fdr" %in% names(df); has_cpg_fdr <- "fdr" %in% names(df)
      has_dbeta <- "delta_beta" %in% names(df); has_dmr_dbeta <- "dmr_delta_beta" %in% names(df)
      has_direction <- "direction" %in% names(df); has_consistency <- "direction_consistency" %in% names(df)
      has_feature <- "feature" %in% names(df); has_island <- "island" %in% names(df)

      tagList(
        div(class = "card",
            div(class = "card-title", icon("list-check"), "Candidate CpG filters"),
            fluidRow(
              column(4, selectizeInput(ns("c4_module"), "Module", choices = c("All" = "__all__", mods), selected = "__all__", multiple = TRUE)),
              column(4, if (has_dmr_fdr) numericInput(ns("c4_dmr_fdr"), "DMR FDR <=", value = 1, min = 0, max = 1, step = 0.01)),
              column(4, if (has_cpg_fdr) numericInput(ns("c4_cpg_fdr"), "CpG FDR <=", value = 1, min = 0, max = 1, step = 0.01))
            ),
            fluidRow(
              column(4, if (has_dbeta) numericInput(ns("c4_min_dbeta"), "Min |Delta-Beta| (CpG)", value = 0, min = 0, max = 1, step = 0.01)),
              column(4, if (has_dmr_dbeta) numericInput(ns("c4_min_dmr_dbeta"), "Min |Delta-Beta| (DMR)", value = 0, min = 0, max = 1, step = 0.01)),
              column(4, if (has_direction) radioButtons(ns("c4_direction"), "Direction", inline = TRUE, choices = c("All" = "all", "Hyper" = "hyper", "Hypo" = "hypo"), selected = "all"))
            ),
            fluidRow(
              column(4, if (has_consistency) radioButtons(ns("c4_consistency"), "Consistency", inline = TRUE, choices = c("All" = "all", "Consistent" = "consistent"), selected = "all")),
              column(4, if (has_feature) checkboxGroupInput(ns("c4_feature"), "Genomic region", choices = sort(unique(stats::na.omit(df$feature))))),
              column(4, if (has_island) checkboxGroupInput(ns("c4_island"), "CpG island context", choices = sort(unique(stats::na.omit(df$island)))))
            ),
            tags$h5("Candidate ranking mode"),
            radioButtons(ns("cand_rank_mode"), NULL, inline = TRUE,
                         choices = c("Statistical significance" = "significance", "Effect size" = "effect",
                                     "Module membership" = "membership", "Combined score" = "combined"),
                         selected = "significance"),
            uiOutput(ns("cand_formula_ui")),
            actionButton(ns("cand_apply_btn"), "Apply Candidate Filters", icon = icon("filter"), class = "btn-primary btn-sm")
        ),
        uiOutput(ns("cand_results_ui"))
      )
    })

    output$cand_formula_ui <- renderUI({
      req(ov_has_run())
      mode <- input$cand_rank_mode %||% "significance"
      df <- overlap_run()$table
      parts <- character(0)
      if (mode %in% c("significance", "combined")) {
        if ("dmr_fdr" %in% names(df)) parts <- c(parts, "-log10(DMR FDR)")
        if ("fdr" %in% names(df)) parts <- c(parts, "-log10(CpG FDR)")
      }
      if (mode %in% c("effect", "combined")) {
        if ("delta_beta" %in% names(df)) parts <- c(parts, "10x|Delta-Beta CpG|")
        if ("dmr_delta_beta" %in% names(df)) parts <- c(parts, "10x|Delta-Beta DMR|")
      }
      if (mode %in% c("membership", "combined") && "kme" %in% names(df)) parts <- c(parts, "|kME|")
      if (identical(mode, "combined") && "direction_consistency" %in% names(df)) parts <- c(parts, "+1 if direction-consistent")
      if (length(parts) == 0) return(div(class = "empty-note", icon("triangle-exclamation"), "None of this ranking mode's factors are available in the loaded data - candidates will be shown unranked (score = 0)."))
      div(class = "empty-note", icon("ranking-star"), strong("Score = "), paste(parts, collapse = " + "))
    })

    candidate_result <- eventReactive(input$cand_apply_btn, {
      ov <- overlap_run(); req(ov)
      df <- ov$table

      sel <- input$c4_module %||% "__all__"
      if (!("__all__" %in% sel)) df <- df[df$module %in% sel, , drop = FALSE]
      if ("dmr_fdr" %in% names(df) && !is.null(input$c4_dmr_fdr)) df <- df[is.na(df$dmr_fdr) | df$dmr_fdr <= input$c4_dmr_fdr, , drop = FALSE]
      if ("fdr" %in% names(df) && !is.null(input$c4_cpg_fdr)) df <- df[is.na(df$fdr) | df$fdr <= input$c4_cpg_fdr, , drop = FALSE]
      if ("delta_beta" %in% names(df) && isTRUE((input$c4_min_dbeta %||% 0) > 0)) df <- df[is.na(df$delta_beta) | abs(df$delta_beta) >= input$c4_min_dbeta, , drop = FALSE]
      if ("dmr_delta_beta" %in% names(df) && isTRUE((input$c4_min_dmr_dbeta %||% 0) > 0)) df <- df[is.na(df$dmr_delta_beta) | abs(df$dmr_delta_beta) >= input$c4_min_dmr_dbeta, , drop = FALSE]
      if (!identical(input$c4_direction %||% "all", "all") && "direction" %in% names(df)) df <- df[!is.na(df$direction) & df$direction == input$c4_direction, , drop = FALSE]
      if (identical(input$c4_consistency %||% "all", "consistent") && "direction_consistency" %in% names(df)) df <- df[!is.na(df$direction_consistency) & df$direction_consistency == "Consistent", , drop = FALSE]
      if (length(input$c4_feature %||% character(0)) > 0 && "feature" %in% names(df)) df <- df[df$feature %in% input$c4_feature, , drop = FALSE]
      if (length(input$c4_island %||% character(0)) > 0 && "island" %in% names(df)) df <- df[df$island %in% input$c4_island, , drop = FALSE]
      validate(need(nrow(df) > 0, "No candidate CpGs remain after applying these filters. Try relaxing a threshold above."))

      mode <- input$cand_rank_mode %||% "significance"
      score <- rep(0, nrow(df)); factors_used <- character(0)
      if (mode %in% c("significance", "combined")) {
        if ("dmr_fdr" %in% names(df)) { score <- score + -log10(pmax(df$dmr_fdr, 1e-300)); factors_used <- c(factors_used, "-log10(DMR FDR)") }
        if ("fdr" %in% names(df)) { score <- score + -log10(pmax(df$fdr, 1e-300)); factors_used <- c(factors_used, "-log10(CpG FDR)") }
      }
      if (mode %in% c("effect", "combined")) {
        if ("delta_beta" %in% names(df)) { score <- score + 10 * abs(df$delta_beta); factors_used <- c(factors_used, "10x|Delta-Beta CpG|") }
        if ("dmr_delta_beta" %in% names(df)) { score <- score + 10 * abs(df$dmr_delta_beta); factors_used <- c(factors_used, "10x|Delta-Beta DMR|") }
      }
      if (mode %in% c("membership", "combined") && "kme" %in% names(df)) { score <- score + abs(df$kme); factors_used <- c(factors_used, "|kME|") }
      if (identical(mode, "combined") && "direction_consistency" %in% names(df)) {
        score <- score + ifelse(!is.na(df$direction_consistency) & df$direction_consistency == "Consistent", 1, 0)
        factors_used <- c(factors_used, "+1 if direction-consistent")
      }
      df$candidate_score <- round(score, 4)
      df <- df[order(-df$candidate_score), , drop = FALSE]
      list(table = df, factors_used = unique(factors_used), mode = mode, n = nrow(df))
    }, ignoreInit = TRUE)

    cand_has_run <- reactiveVal(FALSE)
    observeEvent(input$cand_apply_btn, cand_has_run(TRUE), ignoreInit = TRUE)
    observeEvent(input$load_btn, cand_has_run(FALSE), ignoreInit = TRUE)
    observeEvent(input$overlap_run_btn, cand_has_run(FALSE), ignoreInit = TRUE)

    output$cand_results_ui <- renderUI({
      if (!cand_has_run()) return(NULL)
      r <- candidate_result()
      tagList(
        p(class = "empty-note", icon("ranking-star"), strong(format(r$n, big.mark = ",")), " candidate CpGs after filtering, ranked by: ",
          if (length(r$factors_used) > 0) paste(r$factors_used, collapse = " + ") else "(no ranking factors available; unranked)", "."),
        div(class = "card",
            div(class = "card-title", icon("star"), "Prioritized candidate CpGs"),
            div(class = "table-toolbar", downloadButton(ns("cand_download"), "Candidates (CSV)", class = "btn-sm")),
            DT::dataTableOutput(ns("cand_table"))
        )
      )
    })

    output$cand_table <- DT::renderDataTable({
      req(cand_has_run())
      df <- mcd_pretty(candidate_result()$table)
      numcols <- intersect(c("DMR FDR", "DMR p-value", "CpG p-value", "CpG FDR", "Delta Beta", "DMR Delta Beta", "Candidate Score", "kME"), colnames(df))
      dt <- DT::datatable(df, rownames = FALSE, filter = "top", options = list(scrollX = TRUE, pageLength = 10), class = "stripe hover compact")
      if (length(numcols) > 0) dt <- DT::formatSignif(dt, columns = numcols, digits = 4)
      dt
    })
    outputOptions(output, "cand_table", suspendWhenHidden = FALSE)

    output$cand_download <- downloadHandler(
      filename = function() "prioritized_candidate_cpgs.csv",
      content = function(file) utils::write.csv(candidate_result()$table, file, row.names = FALSE)
    )

    ## ===================== 5. Visualization ===============================

    output$viz_tab_ui <- renderUI({
      if (!has_loaded()) return(div(class = "empty-note", icon("circle-info"), "Load a data source on the \"Data & Filters\" tab first."))
      tagList(
        mcd_viz_block(ns, "viz_modbar", "Overlapping CpGs by module", "chart-column", "Generate Module Overlap Bar Plot", "Requires \"Module-DMR Overlap\" to have been run."),
        mcd_viz_block(ns, "viz_heatmap", "Module enrichment heatmap", "table-cells", "Generate Enrichment Heatmap", "-log10(p) per module, from the same Module-DMR Overlap run."),
        mcd_viz_block(ns, "viz_dmrbar", "Candidate CpGs per DMR", "map-location-dot", "Generate DMR Overlap Plot", "Top 20 DMRs by number of overlapping candidate CpGs. Requires \"DMR-CpG Overlap\" to have been run."),
        mcd_viz_block(ns, "viz_direction", "Effect direction", "scale-balanced", "Generate Direction Plot", "Hyper- vs hypomethylated candidate CpGs."),
        mcd_viz_block(ns, "viz_annot", "Genomic annotation distribution", "dna", "Generate Annotation Distribution Plot", "CpG island context and/or genomic region, when available in the loaded data."),
        mcd_viz_block(ns, "viz_volcano", "Candidate CpG plot", "magnifying-glass-chart", "Generate Candidate Volcano", "Delta-Beta vs -log10(significance); prioritized candidates (Candidate CpGs tab) are highlighted when available.")
      )
    })

    ## Registers a plot's button/output pair: gates on its own upstream analysis (`requires()`), gets a has-run flag, renderPlot, and PNG download.
    register_viz <- function(id_prefix, requires, plot_fn, filename) {
      has_run <- reactiveVal(FALSE)
      observeEvent(input[[paste0(id_prefix, "_btn")]], has_run(TRUE), ignoreInit = TRUE)
      observeEvent(input$load_btn, has_run(FALSE), ignoreInit = TRUE)

      output[[paste0(id_prefix, "_ui")]] <- renderUI({
        if (!has_run()) return(NULL)
        tagList(
          div(class = "table-toolbar", downloadButton(ns(paste0(id_prefix, "_download")), "Plot (PNG)", class = "btn-sm")),
          withSpinner(plotOutput(ns(paste0(id_prefix, "_plot")), height = 380), color = "#2563EB", type = 6)
        )
      })
      plot_obj <- reactive({
        req(has_run())
        ok <- requires()
        validate(need(isTRUE(ok$ok), ok$reason %||% "Prerequisite analysis has not been run yet."))
        plot_fn()
      })
      output[[paste0(id_prefix, "_plot")]] <- renderPlot({ plot_obj() })
      output[[paste0(id_prefix, "_download")]] <- downloadHandler(
        filename = function() filename,
        content = function(file) ggsave(file, plot = plot_obj(), width = 8, height = 5.5, dpi = 300, bg = "white")
      )
    }

    register_viz("viz_modbar", function() list(ok = modov_has_run(), reason = "Run \"Module-DMR Overlap\" first."),
                 function() mcd_plot_module_bar(module_overlap_run()$table), "module_dmr_overlap_bar.png")
    register_viz("viz_heatmap", function() list(ok = modov_has_run(), reason = "Run \"Module-DMR Overlap\" first."),
                 function() mcd_plot_enrichment_heatmap(module_overlap_run()$table), "module_enrichment_heatmap.png")
    register_viz("viz_dmrbar", function() list(ok = ov_has_run(), reason = "Run \"DMR-CpG Overlap\" first."),
                 function() mcd_plot_dmr_bar(overlap_run()$table), "candidates_per_dmr.png")
    register_viz("viz_direction", function() list(ok = ov_has_run(), reason = "Run \"DMR-CpG Overlap\" first."),
                 function() mcd_plot_direction(overlap_run()$table), "direction_distribution.png")
    register_viz("viz_annot", function() list(ok = ov_has_run(), reason = "Run \"DMR-CpG Overlap\" first."),
                 function() mcd_plot_annotation_dist(overlap_run()$table), "annotation_distribution.png")
    register_viz("viz_volcano", function() list(ok = ov_has_run(), reason = "Run \"DMR-CpG Overlap\" first."),
                 function() mcd_plot_candidate_volcano(overlap_run()$table, if (isTRUE(cand_has_run())) candidate_result()$table$cpg else NULL),
                 "candidate_volcano.png")

    ## Read by the Assistant sub-module, same as every other analysis tab.
    observe({
      if (!isTRUE(ov_has_run())) return()
      r <- tryCatch(overlap_run(), error = function(e) NULL)
      if (is.null(r) || is.null(results)) return()
      results$candidate_cpgs <- list(
        n_overlap_cpgs = r$n_overlap_cpgs, n_overlap_dmrs = r$n_overlap_dmrs,
        n_dmr_passing = r$n_dmr_passing, source = loaded()$source
      )
    })
  })
}
