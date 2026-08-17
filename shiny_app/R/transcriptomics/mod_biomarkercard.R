## R/transcriptomics/mod_biomarkercard.R
## Submodule: Biomarker Card - a single-gene transcriptomic biomarker
## profile.
##
## Gene/transcript -> expression evidence in the loaded dataset -> saved
## Differential Expression evidence -> Candidate Gene Identification /
## ML Feature Selection membership -> Diagnostic Classifier performance
## (only when it was actually run on a matching gene panel) -> a compact
## pathway/GO interpretation -> a transparent, un-scored evidence
## checklist -> sources -> downloadable report.
##
## Deliberately mirrors the UI/UX of R/methylomics/mod_methyl_biomarkercard.R
## (same card/tab/spinner/table conventions) but NOT its biology: there is
## no methylation value, DMP/DMR, CpG island, or array manifest anywhere in
## this file. The transcriptomics equivalent of "CpG -> beta -> DMP" is
## "gene -> expression value -> DEG" (mod_dge.R), and the equivalent of the
## methylation card's own preloaded SVA/bacon pipeline is this app's shared
## `results` store, written live by mod_dge.R / mod_candidates.R /
## mod_featureselection.R / mod_diagnostic.R as those tabs are run this
## session - there is no separate offline transcriptomics pipeline output
## to read the way the methylation card reads preloaded DMP/DMR files.
##
## Every performance number shown (AUC, CV-AUC) is read directly from
## results$diagnostic, never recomputed or invented here - if no Diagnostic
## Classifier run matches the gene/signature in view, the card says
## "Not calculated", not a fabricated figure. Same rule for pathway
## evidence: KEGG/GO/Reactome lookups are live, opt-in, and fail soft to an
## explicit reason rather than a silent guess.

## ---- Gene identity (org.Hs.eg.db - already library()'d in global.R) -------

.tbc_identity_cache <- new.env(parent = emptyenv())

tbc_gene_identity <- function(symbol) {
  if (is.null(symbol) || is.na(symbol) || !nzchar(symbol)) return(list(ok = FALSE, reason = "No gene symbol given."))
  cached <- .tbc_identity_cache[[symbol]]
  if (!is.null(cached)) return(cached)
  res <- tryCatch({
    map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = symbol, keytype = "SYMBOL",
                                                    columns = c("ENTREZID", "ENSEMBL", "GENENAME")))
    map <- map[!is.na(map$ENTREZID), , drop = FALSE]
    if (nrow(map) == 0) return(list(ok = FALSE, reason = sprintf("No NCBI Entrez Gene entry found for symbol \"%s\".", symbol)))
    list(ok = TRUE, symbol = symbol, entrez = map$ENTREZID[1], ensembl = map$ENSEMBL[1], genename = map$GENENAME[1])
  }, error = function(e) e)
  if (inherits(res, "error")) res <- list(ok = FALSE, reason = sprintf("Could not resolve gene identity: %s", conditionMessage(res)))
  .tbc_identity_cache[[symbol]] <- res
  res
}

tbc_go_terms <- function(entrez, n = 6) {
  if (is.null(entrez) || is.na(entrez) || !nzchar(entrez)) return(NULL)
  if (!requireNamespace("GO.db", quietly = TRUE)) return(NULL)
  res <- tryCatch({
    go <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = entrez, keytype = "ENTREZID", columns = c("GOALL", "ONTOLOGYALL")))
    go <- go[!is.na(go$GOALL) & go$ONTOLOGYALL == "BP", , drop = FALSE]
    ids <- unique(go$GOALL)
    if (length(ids) == 0) return(NULL)
    terms <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = ids, keytype = "GOID", columns = "TERM"))
    utils::head(unique(stats::na.omit(terms$TERM)), n)
  }, error = function(e) NULL)
  if (inherits(res, "error")) return(NULL)
  res
}

## keggList()'s ~370-row pathway-name table, cached process-wide - own copy,
## not shared with the methylation card's identical-purpose cache, matching
## this app's "each module owns its own data intake" convention.
.tbc_kegg_names_cache <- new.env(parent = emptyenv())
tbc_kegg_pathway_names <- function() {
  if (!is.null(.tbc_kegg_names_cache$v)) return(.tbc_kegg_names_cache$v)
  if (!requireNamespace("KEGGREST", quietly = TRUE)) return(NULL)
  v <- tryCatch(KEGGREST::keggList("pathway", "hsa"), error = function(e) NULL)
  .tbc_kegg_names_cache$v <- v
  v
}
tbc_kegg_pathways_for_gene <- function(entrez) {
  if (is.null(entrez) || is.na(entrez) || !nzchar(entrez)) return(list(ok = FALSE, pathways = NULL, reason = "No NCBI Gene ID available for KEGG lookup."))
  if (!requireNamespace("KEGGREST", quietly = TRUE)) return(list(ok = FALSE, pathways = NULL, reason = "KEGGREST is not installed in this deployment."))
  res <- tryCatch({
    links <- KEGGREST::keggLink("pathway", sprintf("hsa:%s", entrez))
    key <- sub("^path:", "", unname(links))
    names_map <- tbc_kegg_pathway_names()
    nm <- if (!is.null(names_map)) unname(names_map[key]) else rep(NA_character_, length(key))
    nm <- sub(" - Homo sapiens.*$", "", nm)
    data.frame(id = key, name = nm, stringsAsFactors = FALSE)
  }, error = function(e) e)
  if (inherits(res, "error")) return(list(ok = FALSE, pathways = NULL, reason = sprintf("KEGG lookup failed: %s", conditionMessage(res))))
  list(ok = TRUE, pathways = res, reason = NULL)
}

## Same UniProt-fallback-loop rationale as the methylation card's own
## bc_reactome_pathways_for_gene(): SYMBOL -> UNIPROT is 1:many, tries each
## until one has Reactome pathway data.
tbc_reactome_pathways_for_gene <- function(symbol) {
  if (is.null(symbol) || is.na(symbol) || !nzchar(symbol)) return(list(ok = FALSE, pathways = NULL, reason = "No gene symbol available for Reactome lookup."))
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, pathways = NULL, reason = "httr2 is not installed in this deployment."))
  uniprot_ids <- tryCatch({
    m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = symbol, keytype = "SYMBOL", columns = "UNIPROT"))
    unique(stats::na.omit(m$UNIPROT))
  }, error = function(e) character(0))
  if (length(uniprot_ids) == 0) return(list(ok = FALSE, pathways = NULL, reason = sprintf("No UniProt ID found for gene symbol \"%s\".", symbol)))
  for (uid in uniprot_ids) {
    res <- tryCatch({
      url <- sprintf("https://reactome.org/ContentService/data/mapping/UniProt/%s/pathways?species=9606", uid)
      resp <- httr2::request(url) %>% httr2::req_timeout(15) %>%
        httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
      if (httr2::resp_status(resp) != 200) NULL else httr2::resp_body_json(resp, simplifyVector = TRUE)
    }, error = function(e) NULL)
    if (!is.null(res) && is.data.frame(res) && nrow(res) > 0) {
      out <- res[, intersect(c("stId", "displayName"), colnames(res)), drop = FALSE]
      return(list(ok = TRUE, pathways = out, reason = NULL, uniprot = uid))
    }
  }
  list(ok = TRUE, pathways = data.frame(stId = character(0), displayName = character(0)), reason = NULL)
}

tbc_pathway_evidence <- function(symbol, entrez) {
  kegg <- tbc_kegg_pathways_for_gene(entrez)
  reactome <- tbc_reactome_pathways_for_gene(symbol)
  go <- tbc_go_terms(entrez)
  list(kegg = kegg, reactome = reactome, go = go)
}

## ---- Sample-sheet column detection (own copy, same rationale as the
## methylation card's bc_find_col()/*_PATTERNS - each module owns its own
## data intake rather than sharing helpers across mod_*.R files). ----------

TBC_GROUP_PATTERNS <- c("^group$", "^status$", "^phenotype$", "^disease$", "^condition$", "^diagnosis$", "group")
TBC_SEX_PATTERNS <- c("^sex$", "^gender$", "sex")

tbc_find_col <- function(cols, patterns) {
  for (p in patterns) {
    hit <- cols[grepl(p, cols, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  NULL
}

tbc_pick_case_control <- function(levels) {
  lv <- levels
  case_idx <- which(grepl("case|patient|disease|^ra$|affected|positive", lv, ignore.case = TRUE))
  ctrl_idx <- which(grepl("control|healthy|normal|^hc$|negative", lv, ignore.case = TRUE))
  if (length(case_idx) >= 1 && length(ctrl_idx) >= 1 && case_idx[1] != ctrl_idx[1]) {
    return(list(case = lv[case_idx[1]], control = lv[ctrl_idx[1]]))
  }
  if (length(ctrl_idx) >= 1) {
    other <- setdiff(seq_along(lv), ctrl_idx[1])
    if (length(other) >= 1) return(list(case = lv[other[1]], control = lv[ctrl_idx[1]]))
  }
  list(case = lv[2], control = lv[1])
}

## ---- Expression-scale-aware live evidence (computed directly on the
## loaded dataset$expr - never assumes raw counts behave like continuous
## measurements; see detect_expr_data_type() in global.R, already used by
## the Data Exploration / Preprocessing tabs). ------------------------------

## A quick, transparent single-gene preview only - NOT a substitute for the
## Differential Expression tab's own count-aware limma/DESeq2 model. Raw
## counts are log2(CPM+1)-transformed first (edgeR::cpm(), already
## library()'d) so a gene's expression is never treated as ordinary
## continuous data on its native count scale; already-normalised/log-scale
## data is used as-is. The card is explicit about which of the two
## happened wherever this value is shown.
tbc_gene_expr_values <- function(gene, expr) {
  data_type <- detect_expr_data_type(expr)
  row <- expr[gene, ]
  transformed <- FALSE
  if (identical(data_type, "counts")) {
    cpm_mat <- tryCatch(edgeR::cpm(expr, log = TRUE, prior.count = 1), error = function(e) NULL)
    if (!is.null(cpm_mat)) { row <- cpm_mat[gene, ]; transformed <- TRUE }
  }
  list(values = row, data_type = data_type, transformed = transformed)
}

tbc_live_stats <- function(expr_row, group_vec, case_label, control_label) {
  keep <- !is.na(expr_row) & !is.na(group_vec) & group_vec %in% c(case_label, control_label)
  expr_row <- expr_row[keep]; grp <- group_vec[keep]
  case_vals <- expr_row[grp == case_label]; ctrl_vals <- expr_row[grp == control_label]
  if (length(case_vals) < 2 || length(ctrl_vals) < 2) {
    return(list(ok = FALSE, reason = "Fewer than 2 non-missing samples in one of the two groups - cannot compute a t-test.",
                case_label = case_label, control_label = control_label, n_case = length(case_vals), n_control = length(ctrl_vals)))
  }
  tt <- tryCatch(stats::t.test(case_vals, ctrl_vals), error = function(e) NULL)
  mc <- mean(case_vals); mo <- mean(ctrl_vals)
  list(ok = TRUE, mean_case = mc, mean_control = mo, log2fc = mc - mo,
       p_value = if (!is.null(tt)) tt$p.value else NA_real_,
       ci = if (!is.null(tt)) as.numeric(tt$conf.int) else c(NA_real_, NA_real_),
       n_case = length(case_vals), n_control = length(ctrl_vals),
       direction = if (mc > mo) "Up in case group" else if (mc < mo) "Down in case group" else "No change",
       case_label = case_label, control_label = control_label)
}

## Works identically for the preloaded example cohort or an uploaded one -
## dataset$expr/$meta is the same shape either way (see mod_dataset.R).
tbc_dataset_evidence <- function(gene, dataset) {
  expr <- dataset$expr; meta <- dataset$meta
  if (is.null(expr) || is.null(meta)) return(list(ok = FALSE, reason = "No expression matrix / sample metadata is currently loaded on the Dataset tab."))
  if (!gene %in% rownames(expr)) return(list(ok = FALSE, reason = "This gene is not present in the currently loaded expression matrix."))
  samp_ids <- colnames(expr)
  if (!"sample" %in% colnames(meta) || !all(samp_ids %in% as.character(meta$sample))) {
    return(list(ok = FALSE, reason = "Could not align the sample metadata to the expression matrix columns (no matching \"sample\" column)."))
  }
  meta <- meta[match(samp_ids, as.character(meta$sample)), , drop = FALSE]
  gv <- tbc_gene_expr_values(gene, expr)
  group_col <- tbc_find_col(colnames(meta), TBC_GROUP_PATTERNS)
  if (is.null(group_col)) return(list(ok = FALSE, reason = "Could not auto-detect a case/control grouping column in the loaded sample metadata.", values = gv$values, data_type = gv$data_type))
  group_vec <- as.character(meta[[group_col]])
  levels_present <- unique(group_vec[!is.na(group_vec)])
  if (length(levels_present) != 2) {
    return(list(ok = FALSE, reason = sprintf("The detected grouping column (\"%s\") does not have exactly two levels (found: %s).",
                                              group_col, paste(levels_present, collapse = ", "))))
  }
  cc <- tbc_pick_case_control(levels_present)
  overall <- tbc_live_stats(gv$values, group_vec, cc$case, cc$control)
  sex_col <- tbc_find_col(colnames(meta), TBC_SEX_PATTERNS)
  sex_vec <- if (!is.null(sex_col)) as.character(meta[[sex_col]]) else NULL
  by_sex <- NULL
  if (!is.null(sex_vec)) {
    sex_levels <- unique(sex_vec[!is.na(sex_vec)])
    if (length(sex_levels) > 0) {
      by_sex <- stats::setNames(lapply(sex_levels, function(s) {
        idx <- which(sex_vec == s)
        tbc_live_stats(gv$values[idx], group_vec[idx], cc$case, cc$control)
      }), sex_levels)
    }
  }
  list(ok = TRUE, values = gv$values, data_type = gv$data_type, transformed = gv$transformed,
       group_vec = group_vec, sex_vec = sex_vec, group_col = group_col, sex_col = sex_col,
       case_label = cc$case, control_label = cc$control, overall = overall, by_sex = by_sex)
}

## ---- Saved Differential Expression evidence (results$dge_runs, written
## live by mod_dge.R this session - the transcriptomics equivalent of the
## methylation card's preloaded SVA/bacon DMP pipeline, except this one is
## a live in-session analysis, not a static offline file). -----------------

tbc_dge_matches <- function(gene, results) {
  runs <- results$dge_runs %||% list()
  if (length(runs) == 0) return(NULL)
  hits <- lapply(names(runs), function(rid) {
    r <- runs[[rid]]
    row <- r$table[r$table$gene == gene, , drop = FALSE]
    if (nrow(row) == 0) return(NULL)
    data.frame(run = rid, contrast = r$contrast, method = r$method, n_samples = r$n_samples,
               logFC = row$logFC[1], adj.P.Val = row$adj.P.Val[1], direction = row$direction[1],
               stringsAsFactors = FALSE)
  })
  hits <- hits[!vapply(hits, is.null, logical(1))]
  if (length(hits) == 0) return(NULL)
  do.call(rbind, hits)
}

## ---- Candidate Gene Identification / ML Feature Selection membership
## (results$candidates, results$featureselection - both written live by
## their own tabs this session). --------------------------------------------

tbc_candidate_status <- function(gene, results) {
  cand <- results$candidates
  if (is.null(cand)) return(list(any = FALSE, in_female = FALSE, in_male = FALSE, in_final = FALSE, final_selection = NA_character_))
  list(
    any = (!is.null(cand$female) && gene %in% cand$female$genes) || (!is.null(cand$male) && gene %in% cand$male$genes),
    in_female = !is.null(cand$female) && gene %in% cand$female$genes,
    in_male = !is.null(cand$male) && gene %in% cand$male$genes,
    in_final = !is.null(cand$final) && gene %in% cand$final$genes,
    final_selection = if (!is.null(cand$final)) cand$final$selection else NA_character_
  )
}

tbc_signature_membership <- function(gene, results) {
  fs <- results$featureselection
  if (is.null(fs)) return(list())
  sexes <- intersect(c("female", "male", "pooled"), names(fs))
  stats::setNames(lapply(sexes, function(s) {
    genes <- fs[[s]]$consensus_genes
    list(in_signature = !is.null(genes) && gene %in% genes, signature_size = length(genes), genes = genes)
  }), sexes)
}

## ---- Diagnostic Classifier performance (results$diagnostic - only
## populated once that tab has actually been run on a gene panel; every
## number here is read verbatim from that store, never recomputed or
## estimated here). NULL/"Not calculated" propagates all the way to the UI
## whenever no saved run's panel contains this gene. -------------------------

tbc_diagnostic_lookup <- function(gene, results) {
  diag <- results$diagnostic
  if (is.null(diag)) return(list())
  sexes <- names(diag)
  out <- stats::setNames(lapply(sexes, function(s) {
    r <- diag[[s]]
    list(in_panel = !is.null(r$genes) && gene %in% r$genes, panel_size = r$n_input, n_samples = r$n_samples,
         lr_auc = r$lr_auc, lr_cv_auc = r$lr_cv_auc, enet_auc = r$enet_auc, enet_cv_auc = r$enet_cv_auc,
         rf_auc = r$rf_auc, rf_cv_auc = r$rf_cv_auc, svm_auc = r$svm_auc, svm_cv_auc = r$svm_cv_auc, genes = r$genes)
  }), sexes)
  Filter(function(x) isTRUE(x$in_panel), out)
}

## ---- Plots -----------------------------------------------------------------

tbc_plot_expression_dist <- function(df, y_label, facet_sex = FALSE) {
  if (isTRUE(facet_sex) && "sex" %in% names(df)) df <- df[!is.na(df$sex), , drop = FALSE]
  df <- df[!is.na(df$expr) & !is.na(df$group), , drop = FALSE]
  validate(need(nrow(df) > 0, "No non-missing expression values are available for this gene in the current groups."))
  p <- ggplot(df, aes(x = group, y = expr, fill = group)) +
    geom_violin(alpha = 0.35, color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.08, alpha = 0.5, size = 1.1) +
    scale_fill_manual(values = arthomix_pair(sort(unique(df$group)))) +
    labs(x = NULL, y = y_label) + theme_arthomix() + theme(legend.position = "none")
  if (isTRUE(facet_sex) && "sex" %in% names(df)) p <- p + facet_wrap(~sex)
  p
}

## Reuses the same run table mod_dge.R saved (results$dge_runs[[run]]$table:
## gene/logFC/adj.P.Val/direction) - the card never refits its own volcano
## model, only highlights the selected gene within an already-computed run.
tbc_plot_volcano_highlight <- function(tab, gene) {
  validate(need(!is.null(tab) && nrow(tab) > 0, "No Differential Expression run is available to plot."))
  tab$neglog10p <- -log10(pmax(tab$adj.P.Val, 1e-300))
  cols <- c("Up" = ARTHOMIX_COLORS$red, "Down" = ARTHOMIX_COLORS$blue, "Not significant" = ARTHOMIX_COLORS$ink_muted)
  p <- ggplot(tab, aes(x = logFC, y = neglog10p, color = direction)) +
    geom_point(alpha = 0.55, size = 1.4) +
    scale_color_manual(values = cols, name = NULL) +
    labs(x = "log2 fold-change", y = expression(-log[10]~"(adjusted p-value)")) +
    theme_arthomix()
  hit <- tab[tab$gene == gene, , drop = FALSE]
  if (nrow(hit) > 0) {
    p <- p + geom_point(data = hit, aes(x = logFC, y = neglog10p), color = ARTHOMIX_STATUS$critical, size = 3.4, inherit.aes = FALSE) +
      ggrepel::geom_text_repel(data = hit, aes(x = logFC, y = neglog10p, label = gene), color = ARTHOMIX_STATUS$critical, fontface = "bold", inherit.aes = FALSE)
  }
  p
}

tbc_ggsave_datauri <- function(plot, width = 8, height = 3.4, dpi = 110) {
  if (is.null(plot) || !requireNamespace("base64enc", quietly = TRUE)) return(NULL)
  tf <- tempfile(fileext = ".png")
  ok <- tryCatch({ ggplot2::ggsave(tf, plot = plot, width = width, height = height, dpi = dpi, bg = "white"); TRUE }, error = function(e) FALSE)
  if (!ok || !file.exists(tf)) return(NULL)
  uri <- base64enc::dataURI(file = tf, mime = "image/png")
  unlink(tf)
  uri
}

## ---- Display helpers -------------------------------------------------------

tbc_fmt_field <- function(x) {
  if (is.null(x) || length(x) == 0) return("Not available")
  if (length(x) == 1 && is.na(x)) return("Not available")
  if (is.character(x) && !nzchar(trimws(x))) return("Not available")
  as.character(x)
}
tbc_kv_table <- function(pairs) {
  df <- data.frame(Field = names(pairs), Value = vapply(pairs, tbc_fmt_field, character(1)), stringsAsFactors = FALSE)
  DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE, ordering = FALSE), class = "stripe hover compact")
}

## ---- Section builders (plain tags - reused identically on-screen and in
## the downloadable report). ------------------------------------------------

## Session-wide headline numbers (Section 12 of the spec: "Transcriptomic
## Biomarker Summary") - NOT gene-specific, computed from whatever this
## session's shared `results` store already holds; any figure that hasn't
## actually been computed this session reads "Not calculated", never a
## placeholder number.
tbc_section_session_summary <- function(results, ext = NULL) {
  dge <- results$dge
  cand <- results$candidates
  fs <- results$featureselection
  diag <- results$diagnostic

  n_significant <- if (!is.null(dge)) dge$n_significant else NA_integer_
  top_candidate <- if (!is.null(dge) && length(dge$top_hits) > 0) dge$top_hits[1] else NA_character_
  n_candidates <- if (!is.null(cand$final)) cand$final$n_candidates else NA_integer_
  sig_size <- if (!is.null(fs)) {
    sizes <- vapply(intersect(c("female", "male", "pooled"), names(fs)), function(s) fs[[s]]$n_consensus %||% NA_integer_, numeric(1))
    if (length(sizes) > 0 && any(!is.na(sizes))) max(sizes, na.rm = TRUE) else NA_integer_
  } else NA_integer_
  best_auc <- if (!is.null(diag)) {
    aucs <- unlist(lapply(diag, function(r) c(r$lr_cv_auc, r$enet_cv_auc, r$rf_cv_auc, r$svm_cv_auc)))
    if (length(aucs) > 0 && any(!is.na(aucs))) max(aucs, na.rm = TRUE) else NA_real_
  } else NA_real_
  n_pathways <- if (!is.null(ext)) {
    n_kegg <- if (isTRUE(ext$kegg$ok)) nrow(ext$kegg$pathways) else 0
    n_reactome <- if (isTRUE(ext$reactome$ok)) nrow(ext$reactome$pathways) else 0
    n_kegg + n_reactome
  } else NA_integer_

  pairs <- list(
    "Significant genes (latest Differential Expression run)" = n_significant,
    "Top candidate (latest Differential Expression run)" = top_candidate,
    "Candidate biomarkers identified (Candidate Gene Identification, final set)" = n_candidates,
    "Largest multi-gene signature size (ML Feature Selection consensus)" = if (!is.na(sig_size)) sig_size else NA,
    "Best classification performance seen (CV-AUC, Diagnostic Classifier)" = if (!is.na(best_auc)) sprintf("%.3f", best_auc) else NA,
    "External validation" = "Not available in this deployment",
    "Biological pathways found for the current gene" = if (!is.na(n_pathways)) sprintf("%d pathway(s)", n_pathways) else "Not yet looked up"
  )
  div(class = "card",
      div(class = "card-title", icon("clipboard-list"), "Transcriptomic Biomarker Summary"),
      p(class = "submodule-desc", "Session-wide headline numbers, read directly from whatever Differential Expression / Candidate Gene Identification / ML Feature Selection / Diagnostic Classifier runs have actually been computed this session - not gene-specific, and never fabricated when a figure hasn't been computed yet."),
      tbc_kv_table(pairs)
  )
}

tbc_section_identity <- function(d) {
  gi <- d$gene_identity
  pairs <- list(
    "Gene symbol" = d$gene,
    "Full name" = if (isTRUE(gi$ok)) gi$genename else NA,
    "NCBI Entrez Gene ID" = if (isTRUE(gi$ok)) gi$entrez else NA,
    "Ensembl Gene ID" = if (isTRUE(gi$ok)) gi$ensembl else NA,
    "Present in currently loaded expression matrix" = if (d$in_dataset) "Yes" else "No"
  )
  div(class = "card",
      div(class = "card-title", icon("dna"), "Gene Identity"),
      if (!isTRUE(gi$ok)) div(class = "empty-note", icon("triangle-exclamation"), gi$reason) else NULL,
      tbc_kv_table(pairs)
  )
}

tbc_section_expression_data <- function(dataset, live) {
  n_genes <- tryCatch(nrow(dataset$expr), error = function(e) NA_integer_)
  n_samples <- tryCatch(ncol(dataset$expr), error = function(e) NA_integer_)
  data_type <- if (isTRUE(live$ok) || !is.null(live$data_type)) live$data_type else tryCatch(detect_expr_data_type(dataset$expr), error = function(e) NA_character_)
  type_label <- switch(data_type %||% "unknown",
    counts = "Raw / count-scale data (integer-like values, no negatives) - log2(CPM+1)-transformed here before any statistical comparison, exactly as count-aware RNA-seq workflows require; the Differential Expression tab instead fits a proper count model (DESeq2).",
    already_normalised = "Already-normalised expression data (e.g. batch-corrected, quantile-normalised, or log-scale microarray) - used as-is, with no further normalisation applied here.",
    expression = "Linear- or log-scale expression that has not been shown to already agree across samples - the Preprocessing / Data Exploration tabs offer quantile normalisation for this case.",
    "Could not be determined for the currently loaded expression matrix."
  )
  pairs <- list(
    "Dataset source" = dataset$source,
    "Genes x samples" = if (!is.na(n_genes)) sprintf("%s x %s", format(n_genes, big.mark = ","), n_samples) else NA,
    "Detected data type" = type_label
  )
  div(class = "card",
      div(class = "card-title", icon("table-cells"), "Expression Data Input"),
      p(class = "submodule-desc", "Detection follows the same heuristic used by the Data Exploration and Preprocessing tabs (detect_expr_data_type()), so a gene's evidence below is never computed as if raw counts were ordinary continuous measurements."),
      tbc_kv_table(pairs)
  )
}

tbc_section_discovery_context <- function(dataset, live) {
  group_label <- if (isTRUE(live$ok) || !is.null(live$group_col)) live$group_col else "Not auto-detected in the loaded sample metadata"
  case_label <- if (isTRUE(live$ok)) sprintf("%s (comparison) vs %s (reference)", live$case_label, live$control_label) else "Not available"
  sex_label <- if (isTRUE(live$ok) && !is.null(live$sex_col)) live$sex_col else "Not recorded in the loaded sample metadata"
  n_case <- if (isTRUE(live$ok) && isTRUE(live$overall$ok)) live$overall$n_case else NA
  n_ctrl <- if (isTRUE(live$ok) && isTRUE(live$overall$ok)) live$overall$n_control else NA
  pairs <- list(
    "Dataset / cohort" = dataset$source,
    "Comparison group vs reference group" = case_label,
    "Grouping column used" = group_label,
    "Sample sizes (comparison / reference)" = if (!is.na(n_case)) sprintf("%s / %s", n_case, n_ctrl) else NA,
    "Sex column" = sex_label,
    "Tissue / sample type" = "Not recorded in the loaded sample metadata",
    "Other clinical covariates" = "Not recorded in the loaded sample metadata beyond sex/batch (see the Dataset tab)"
  )
  div(class = "card",
      div(class = "card-title", icon("magnifying-glass-chart"), "Transcriptomic Biomarker Discovery"),
      p(class = "submodule-desc", "The biological condition being distinguished with gene-expression data in the currently loaded dataset - fields read directly from the loaded sample metadata, never inferred."),
      tbc_kv_table(pairs)
  )
}

tbc_section_differential_expression <- function(d) {
  live <- d$live
  dge_hits <- d$dge_hits
  rows <- list()
  if (isTRUE(live$ok) && isTRUE(live$overall$ok)) {
    ov <- live$overall
    rows[[length(rows) + 1]] <- data.frame(
      Source = "Your loaded dataset (quick preview, uncorrected)",
      Gene = d$gene, `Log2 fold-change` = round(ov$log2fc, 3), `P-value` = signif(ov$p_value, 4), `Adjusted p-value` = NA_real_,
      Direction = ov$direction, `Biomarker candidate status` = if (!is.na(ov$p_value) && ov$p_value <= 0.05) "Candidate (uncorrected p <= 0.05)" else "Not significant",
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (!is.null(dge_hits) && nrow(dge_hits) > 0) {
    for (i in seq_len(nrow(dge_hits))) {
      r <- dge_hits[i, ]
      rows[[length(rows) + 1]] <- data.frame(
        Source = sprintf("Saved Differential Expression run (%s, %s)", r$contrast, r$method),
        Gene = d$gene, `Log2 fold-change` = round(r$logFC, 3), `P-value` = NA_real_, `Adjusted p-value` = signif(r$adj.P.Val, 4),
        Direction = r$direction, `Biomarker candidate status` = if (identical(r$direction, "Not significant")) "Not significant" else "Candidate transcriptomic biomarker",
        check.names = FALSE, stringsAsFactors = FALSE)
    }
  }
  body <- if (length(rows) == 0) {
    div(class = "empty-note", icon("circle-info"),
        "No differential-expression evidence is available for this gene yet - it is either absent from the loaded expression matrix, the grouping column could not be auto-detected, or no Differential Expression run this session included it. Run the Differential Expression tab, or check the Dataset tab.")
  } else {
    DT::datatable(do.call(rbind, rows), rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")
  }
  div(class = "card",
      div(class = "card-title", icon("chart-column"), "Differential Expression"),
      p(class = "submodule-desc", "Thresholds are whatever each source actually used - the quick dataset preview uses an uncorrected p <= 0.05 (labelled as such, never presented as FDR-controlled); saved Differential Expression runs use that run's own adjusted-p-value/log2FC cutoffs from the Differential Expression tab. \"Candidate\" here means statistically flagged, not clinically validated."),
      body
  )
}

tbc_section_prioritization <- function(d) {
  cand <- d$candidate_status
  sig <- d$signature_membership
  diag <- d$diagnostic_match
  live_sig <- isTRUE(d$live$ok) && isTRUE(d$live$overall$ok) && !is.na(d$live$overall$p_value) && d$live$overall$p_value <= 0.05
  dge_sig <- !is.null(d$dge_hits) && any(d$dge_hits$direction != "Not significant")
  effect_size <- if (!is.null(d$dge_hits) && nrow(d$dge_hits) > 0) sprintf("|log2FC| = %.2f (%s)", abs(d$dge_hits$logFC[1]), d$dge_hits$contrast[1])
                 else if (isTRUE(d$live$ok) && isTRUE(d$live$overall$ok)) sprintf("|log2FC| = %.2f (quick preview)", abs(d$live$overall$log2fc))
                 else "Not calculated"
  cv_stable <- if (length(diag) > 0) "Yes (see Biomarker Performance below)" else "Not calculated"

  rows <- list(
    c("Statistical significance (dataset or saved Differential Expression run)", if (live_sig || dge_sig) "Yes" else "No"),
    c("Effect size", effect_size),
    c("Candidate Gene Identification membership (WGCNA module x DEG overlap)", if (isTRUE(cand$any)) "Yes" else "No"),
    c("ML Feature Selection consensus membership", if (length(sig) > 0 && any(vapply(sig, function(x) isTRUE(x$in_signature), logical(1)))) "Yes" else "No"),
    c("Classification performance calculated", if (length(diag) > 0) "Yes" else "Not calculated"),
    c("Cross-validation stability", cv_stable),
    c("Independent cohort replication", "Not calculated - no independent external cohort is wired into this app"),
    c("Biological / pathway relevance", "See Biological Interpretation below")
  )
  df <- data.frame(Evidence = vapply(rows, `[`, character(1), 1), Status = vapply(rows, `[`, character(1), 2), stringsAsFactors = FALSE)
  div(class = "card",
      div(class = "card-title", icon("ranking-star"), "Biomarker Candidate Prioritization"),
      p(class = "submodule-desc", "Each row is a transparent, individually-checked evidence dimension, not a combined score - this build does not compute a composite priority score, so none is shown."),
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE, ordering = FALSE), class = "stripe hover compact")
  )
}

tbc_section_signature <- function(d) {
  sig <- d$signature_membership
  if (length(sig) == 0) {
    return(div(class = "card",
        div(class = "card-title", icon("layer-group"), "Single-Gene vs Multi-Gene Signature"),
        div(class = "empty-note", icon("circle-info"), "ML Feature Selection has not been run this session, so no multi-gene signature membership can be shown yet - this gene is only viewable as a single-gene candidate. Run the ML Feature Selection tab to build a consensus signature.")
    ))
  }
  rows <- lapply(names(sig), function(s) {
    x <- sig[[s]]
    data.frame(Stratum = tools::toTitleCase(s),
               `Signature size` = x$signature_size,
               `This gene included` = if (isTRUE(x$in_signature)) "Yes" else "No",
               `Mode` = if (x$signature_size <= 1) "Single-gene candidate" else "Multi-gene transcriptomic signature",
               check.names = FALSE, stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  in_any <- any(vapply(sig, function(x) isTRUE(x$in_signature), logical(1)))
  full_lists <- lapply(names(sig), function(s) {
    x <- sig[[s]]
    if (!isTRUE(x$in_signature) || is.null(x$genes)) return(NULL)
    div(p(strong(sprintf("%s consensus signature (%d genes): ", tools::toTitleCase(s), length(x$genes))), paste(x$genes, collapse = ", ")))
  })
  div(class = "card",
      div(class = "card-title", icon("layer-group"), "Single-Gene vs Multi-Gene Signature"),
      p(class = "submodule-desc", "Membership in each sex-stratified consensus signature from the ML Feature Selection tab (LASSO/Random Forest/SVM-RFE agreement) - correlated expression features can discriminate a phenotype collectively even when no single gene does so alone."),
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE), class = "stripe hover compact"),
      if (in_any) tagList(full_lists) else div(class = "empty-note", icon("circle-info"), "This gene is not part of any consensus signature computed so far this session.")
  )
}

tbc_section_performance <- function(d) {
  diag <- d$diagnostic_match
  if (length(diag) == 0) {
    return(div(class = "card",
        div(class = "card-title", icon("gauge-high"), "Biomarker Performance"),
        div(class = "empty-note", icon("circle-info"), "Validation not performed - no Diagnostic Classifier run this session used a gene panel containing this gene. Run the Diagnostic Classifier tab on a panel that includes it (e.g. its consensus signature above) to populate this section."),
        tbc_kv_table(list("ROC / AUC" = NA, "Sensitivity" = NA, "Specificity" = NA, "Accuracy" = NA, "Cross-validation performance" = NA))
    ))
  }
  rows <- do.call(rbind, lapply(names(diag), function(s) {
    x <- diag[[s]]
    data.frame(
      Stratum = tools::toTitleCase(s), `Panel size` = x$panel_size, `n samples` = x$n_samples,
      `Logistic regression (CV-AUC)` = tbc_fmt_field(x$lr_cv_auc), `Elastic net (CV-AUC)` = tbc_fmt_field(x$enet_cv_auc),
      `Random forest (CV-AUC)` = tbc_fmt_field(x$rf_cv_auc), `SVM (CV-AUC)` = tbc_fmt_field(x$svm_cv_auc),
      check.names = FALSE, stringsAsFactors = FALSE)
  }))
  div(class = "card",
      div(class = "card-title", icon("gauge-high"), "Biomarker Performance"),
      p(class = "submodule-desc", "Cross-validated AUC per model, read directly from the Diagnostic Classifier tab's saved results for whichever gene panel(s) this gene is part of - sensitivity/specificity/ROC curves are only shown on that tab itself (not persisted to the shared results store used here), so they are not duplicated below; open the Diagnostic Classifier tab for the full curve."),
      DT::datatable(rows, rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact"),
      div(class = "empty-note", icon("circle-info"), "Sensitivity, specificity, and full ROC curves - Not calculated here; see the Diagnostic Classifier tab's own \"Result\" panel for the matching run.")
  )
}

tbc_section_discovery_validation <- function(d, dataset) {
  cv_available <- length(d$diagnostic_match) > 0
  pairs <- list(
    "Discovery cohort" = dataset$source,
    "Internally evaluated (held-out test split, Diagnostic Classifier)" = if (cv_available) "Yes - see Biomarker Performance above" else "Not calculated",
    "Cross-validated (k-fold CV, Diagnostic Classifier)" = if (cv_available) "Yes - see Biomarker Performance above" else "Not calculated",
    "Independently replicated in a second in-app dataset" = "Not applicable - only one dataset is loaded at a time in this build",
    "External validation cohort" = "Not available in this deployment - no independent external cohort is currently wired into this app"
  )
  div(class = "card",
      div(class = "card-title", icon("flask-vial"), "Discovery vs Validation"),
      p(class = "submodule-desc", "Passing a p-value/FDR cutoff in the Discovery cohort is not the same as validation - this section states plainly which of Discovery / Internal evaluation / Cross-validation / External validation has actually happened for this gene, this session."),
      tbc_kv_table(pairs)
  )
}

tbc_section_pathways <- function(ext) {
  if (is.null(ext)) {
    return(div(class = "card",
        div(class = "card-title", icon("diagram-project"), "Biological Interpretation"),
        p(class = "submodule-desc", "Pathway, Gene Ontology, and disease-gene-network evidence for this gene - live lookups, opt-in since they're external network calls."),
        div(class = "empty-note", icon("circle-info"), "Not yet looked up - click \"Look Up Biological Interpretation\" below.")
    ))
  }
  go_body <- if (!is.null(ext$go) && length(ext$go) > 0) tags$ul(lapply(ext$go, tags$li)) else div(class = "empty-note", "No Gene Ontology (biological process) terms found, or GO.db is not installed in this deployment.")
  kegg_body <- if (isTRUE(ext$kegg$ok) && nrow(ext$kegg$pathways) > 0) {
    DT::datatable(ext$kegg$pathways, colnames = c("KEGG ID", "Pathway"), rownames = FALSE, options = list(dom = "t", paging = FALSE), class = "stripe hover compact")
  } else div(class = "empty-note", icon("circle-info"), if (isTRUE(ext$kegg$ok)) "No KEGG pathways found for this gene." else (ext$kegg$reason %||% "KEGG lookup unavailable."))
  reactome_body <- if (isTRUE(ext$reactome$ok) && nrow(ext$reactome$pathways) > 0) {
    DT::datatable(ext$reactome$pathways, colnames = c("Reactome ID", "Pathway"), rownames = FALSE, options = list(dom = "t", paging = FALSE), class = "stripe hover compact")
  } else div(class = "empty-note", icon("circle-info"), if (isTRUE(ext$reactome$ok)) "No Reactome pathways found for this gene." else (ext$reactome$reason %||% "Reactome lookup unavailable."))

  div(class = "card",
      div(class = "card-title", icon("diagram-project"), "Biological Interpretation"),
      p(class = "submodule-desc", "Live lookups: Gene Ontology biological-process terms (org.Hs.eg.db / GO.db), KEGG pathway membership (KEGGREST), and Reactome pathway membership (Reactome ContentService)."),
      tags$b("Gene Ontology (biological process)"), go_body,
      tags$div(style = "margin-top:10px;", tags$b("KEGG pathways")), kegg_body,
      tags$div(style = "margin-top:10px;", tags$b("Reactome pathways")), reactome_body
  )
}

tbc_section_interpretation <- function(d, ext = NULL) {
  cand <- d$candidate_status
  sig_any <- length(d$signature_membership) > 0 && any(vapply(d$signature_membership, function(x) isTRUE(x$in_signature), logical(1)))
  dge_sig <- !is.null(d$dge_hits) && any(d$dge_hits$direction != "Not significant")
  perf_done <- length(d$diagnostic_match) > 0
  pathway_done <- FALSE; pathway_label <- "Not yet looked up"
  if (!is.null(ext)) {
    n_kegg <- if (isTRUE(ext$kegg$ok)) nrow(ext$kegg$pathways) else 0
    n_reactome <- if (isTRUE(ext$reactome$ok)) nrow(ext$reactome$pathways) else 0
    pathway_done <- (n_kegg + n_reactome) > 0
    pathway_label <- if (pathway_done) sprintf("%d pathway(s)", n_kegg + n_reactome) else "None found"
  }
  nodes <- list(
    list(label = "Gene", value = d$gene, done = TRUE),
    list(label = "Expression change", value = if (dge_sig) "Differentially expressed" else "Not significant / not tested", done = dge_sig),
    list(label = "Candidate status", value = if (isTRUE(cand$any)) "Candidate biomarker" else "Not a candidate yet", done = isTRUE(cand$any)),
    list(label = "Signature", value = if (sig_any) "In multi-gene signature" else "Single-gene only", done = sig_any),
    list(label = "Performance", value = if (perf_done) "Calculated" else "Not calculated", done = perf_done),
    list(label = "Pathway", value = pathway_label, done = pathway_done)
  )
  div(class = "card",
      div(class = "card-title", icon("route"), "Functional Interpretation"),
      p(class = "submodule-desc", "Cautious, evidence-only chain - each arrow means \"supported by\", never a causal or diagnostic claim."),
      div(style = "display:flex; flex-wrap:wrap; align-items:center; gap:6px;",
          lapply(seq_along(nodes), function(i) {
            n <- nodes[[i]]
            tagList(
              div(style = sprintf("border:1px solid %s; border-radius:8px; padding:8px 12px; min-width:110px; text-align:center; background:%s;",
                                   if (n$done) ARTHOMIX_COLORS$blue else "#cccccc", if (n$done) "#EAF3FB" else "#F5F5F5"),
                  tags$div(style = "font-size:11px; text-transform:uppercase; color:var(--color-ink-secondary);", n$label),
                  tags$div(style = "font-weight:600; font-size:13px;", n$value)
              ),
              if (i < length(nodes)) icon("arrow-right", style = "color:#999;")
            )
          })
      )
  )
}

tbc_section_evidence_summary <- function(d, ext = NULL) {
  live_sig <- isTRUE(d$live$ok) && isTRUE(d$live$overall$ok) && !is.na(d$live$overall$p_value) && d$live$overall$p_value <= 0.05
  dge_sig <- !is.null(d$dge_hits) && any(d$dge_hits$direction != "Not significant")
  pathway_found <- FALSE
  if (!is.null(ext)) {
    n_kegg <- if (isTRUE(ext$kegg$ok)) nrow(ext$kegg$pathways) else 0
    n_reactome <- if (isTRUE(ext$reactome$ok)) nrow(ext$reactome$pathways) else 0
    pathway_found <- (n_kegg + n_reactome) > 0
  }
  rows <- list(
    c("Significant in loaded dataset (quick preview, uncorrected p <= 0.05)", if (live_sig) "Yes" else "No"),
    c("Significant in a saved Differential Expression run", if (dge_sig) "Yes" else "No"),
    c("Candidate Gene Identification membership", if (isTRUE(d$candidate_status$any)) "Yes" else "No"),
    c("ML Feature Selection consensus membership", if (length(d$signature_membership) > 0 && any(vapply(d$signature_membership, function(x) isTRUE(x$in_signature), logical(1)))) "Yes" else "No"),
    c("Diagnostic Classifier performance calculated", if (length(d$diagnostic_match) > 0) "Yes" else "Not calculated"),
    c("Cross-validated performance available", if (length(d$diagnostic_match) > 0) "Yes" else "Not calculated"),
    c("Pathway / GO annotation found", if (is.null(ext)) "Not yet checked" else if (pathway_found) "Yes" else "No"),
    c("Independent external validation", "Not available in this deployment")
  )
  df <- data.frame(Evidence = vapply(rows, `[`, character(1), 1), Status = vapply(rows, `[`, character(1), 2), stringsAsFactors = FALSE)
  div(class = "card",
      div(class = "card-title", icon("clipboard-check"), "Evidence Summary"),
      p(class = "submodule-desc", "Each row is a transparent, individually-computed check - not a combined score."),
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE, ordering = FALSE), class = "stripe hover compact")
  )
}

tbc_section_sources <- function(ext = NULL) {
  used <- c(
    "This app's own loaded expression dataset (Dataset tab) - live per-gene evidence, computed here",
    "org.Hs.eg.db - NCBI Gene ID, Ensembl Gene ID, gene name mapping",
    "This session's shared results store: Differential Expression (mod_dge.R), Candidate Gene Identification (mod_candidates.R), ML Feature Selection (mod_featureselection.R), Diagnostic Classifier (mod_diagnostic.R) - read verbatim, nothing recomputed or re-derived here"
  )
  external_available <- c(
    "Gene Ontology (GO.db) - biological process terms",
    "KEGG - live query via KEGGREST, https://rest.kegg.jp/",
    "Reactome - live query, https://reactome.org/ContentService/"
  )
  external_label <- if (is.null(ext)) "Available (click \"Look Up Biological Interpretation\" to query):" else "Used in this Biomarker Card (external, live-queried):"
  div(class = "card",
      div(class = "card-title", icon("database"), "Database Sources"),
      tags$b("Used in this Biomarker Card (local/in-session):"), tags$ul(lapply(used, tags$li)),
      tags$b(external_label), tags$ul(lapply(external_available, tags$li))
  )
}

tbc_report_css <- function() {
  "body{font-family:-apple-system,Helvetica,Arial,sans-serif; max-width:900px; margin:24px auto; color:#222;}
   .card{border:1px solid #ddd; border-radius:10px; padding:14px 18px; margin-bottom:16px;}
   .card-title{font-weight:700; font-size:15px; margin-bottom:8px;}
   .submodule-desc{color:#666; font-size:12.5px;}
   .empty-note{background:#f6f6f6; border-left:3px solid #999; padding:8px 12px; border-radius:4px; font-size:13px;}
   table{border-collapse:collapse; width:100%;} td,th{border:1px solid #eee; padding:4px 8px; font-size:13px; text-align:left;}"
}

tbc_build_report_tags <- function(d, dataset, results, ext = NULL) {
  dist_plot <- if (isTRUE(d$live$ok) && isTRUE(d$live$overall$ok))
    tryCatch(tbc_plot_expression_dist(data.frame(expr = d$live$values, group = d$live$group_vec), "Expression (analysis scale)"), error = function(e) NULL) else NULL
  volcano_plot <- if (!is.null(d$selected_run_table))
    tryCatch(tbc_plot_volcano_highlight(d$selected_run_table, d$gene), error = function(e) NULL) else NULL

  img_tag <- function(p) {
    uri <- tbc_ggsave_datauri(p)
    if (is.null(uri)) div(class = "empty-note", "Plot unavailable.") else tags$img(src = uri, style = "max-width:100%; height:auto;")
  }

  tagList(
    tags$h2(sprintf("Transcriptomic Biomarker Card: %s", d$gene)),
    tbc_section_session_summary(results, ext),
    tbc_section_identity(d),
    tbc_section_discovery_context(dataset, d$live),
    tbc_section_expression_data(dataset, d$live),
    tbc_section_differential_expression(d),
    if (!is.null(volcano_plot)) div(class = "card", div(class = "card-title", "Differential Expression (Volcano)"), img_tag(volcano_plot)) else NULL,
    if (!is.null(dist_plot)) div(class = "card", div(class = "card-title", "Expression Distribution (Your Dataset)"), img_tag(dist_plot)) else NULL,
    tbc_section_prioritization(d), tbc_section_signature(d), tbc_section_performance(d),
    tbc_section_discovery_validation(d, dataset), tbc_section_pathways(ext),
    tbc_section_interpretation(d, ext), tbc_section_evidence_summary(d, ext), tbc_section_sources(ext),
    tags$p(style = "color:#888; font-size:12px; margin-top:16px;",
           if (is.null(ext)) "Biological Interpretation (GO/KEGG/Reactome) was not looked up before this report was generated - go back to the card and click \"Look Up Biological Interpretation\" first if you want it included."
           else "Sensitivity/specificity/ROC curves are not duplicated here - see the Diagnostic Classifier tab's own Result panel for the matching run.")
  )
}

## ---- Config / UI -----------------------------------------------------

mod_biomarkercard_config <- list(
  id = "biomarkercard", group = "Biomarker modeling",
  title = "Biomarker Card",
  description = "A single-gene transcriptomic biomarker profile: expression evidence in your loaded dataset, saved Differential Expression results, candidate/signature membership, classification performance (only when calculated), and a compact pathway interpretation.",
  icon = "id-card"
)

mod_biomarkercard_ui <- function(id) {
  ns <- NS(id)
  div(
    class = "tx-menu-wrap",
    tabsetPanel(
      id = ns("bmc_subtabs"), type = "tabs",
      tabPanel("Select Biomarker", br(), uiOutput(ns("select_ui"))),
      tabPanel("Biomarker Card", br(), uiOutput(ns("bmc_card_ui")))
    )
  )
}

## ---- Server ------------------------------------------------------------

mod_biomarkercard_server <- function(id, dataset, results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$select_ui <- renderUI({
      dge_runs <- results$dge_runs %||% list()
      fs <- results$featureselection
      cand <- results$candidates
      tagList(
        div(class = "card",
            div(class = "card-title", icon("list-check"), "Find a biomarker"),
            radioButtons(ns("bmc_search_mode"), NULL, inline = TRUE,
                         choices = c("Type a gene symbol" = "gene",
                                     "Browse Differential Expression results" = "dge",
                                     "Browse candidate genes" = "candidates",
                                     "Browse a feature-selected signature" = "signature",
                                     "Upload a gene list" = "upload")),
            conditionalPanel(condition = sprintf("input['%s'] == 'gene'", ns("bmc_search_mode")),
                              textInput(ns("bmc_gene_input"), "Gene symbol", placeholder = "TNF")),
            conditionalPanel(condition = sprintf("input['%s'] == 'dge'", ns("bmc_search_mode")),
                              if (length(dge_runs) == 0) div(class = "empty-note", icon("circle-info"), "No Differential Expression run yet this session - run the Differential Expression tab first, or type a gene symbol directly.")
                              else tagList(
                                selectInput(ns("bmc_dge_run"), "Differential Expression run", choices = stats::setNames(names(dge_runs), vapply(dge_runs, function(r) r$contrast, character(1)))),
                                uiOutput(ns("bmc_dge_results_ui"))
                              )),
            conditionalPanel(condition = sprintf("input['%s' ] == 'candidates'", ns("bmc_search_mode")),
                              if (is.null(cand)) div(class = "empty-note", icon("circle-info"), "Candidate Gene Identification has not been run this session - run that tab first, or type a gene symbol directly.")
                              else tagList(
                                radioButtons(ns("bmc_cand_set"), "Candidate set", inline = TRUE,
                                             choices = Filter(function(v) !is.null(cand[[v]]),
                                                              c("Female" = "female", "Male" = "male", "Final (combined)" = "final"))),
                                uiOutput(ns("bmc_cand_results_ui"))
                              )),
            conditionalPanel(condition = sprintf("input['%s'] == 'signature'", ns("bmc_search_mode")),
                              if (is.null(fs)) div(class = "empty-note", icon("circle-info"), "ML Feature Selection has not been run this session - run that tab first, or type a gene symbol directly.")
                              else tagList(
                                radioButtons(ns("bmc_sig_sex"), "Stratum", inline = TRUE,
                                             choices = stats::setNames(intersect(c("female", "male", "pooled"), names(fs)), tools::toTitleCase(intersect(c("female", "male", "pooled"), names(fs))))),
                                uiOutput(ns("bmc_sig_results_ui"))
                              )),
            conditionalPanel(condition = sprintf("input['%s'] == 'upload'", ns("bmc_search_mode")),
                              p(class = "submodule-desc", "Upload a plain gene-symbol list (one per line, or the first column of a CSV/TSV), or a Diagnostic Classifier RDS export (that tab's own \"Save trained model\" download) - auto-detected by file extension."),
                              fileInput(ns("bmc_upload_file"), "Biomarker list (.csv, .txt, or .rds)", accept = c(".csv", ".txt", ".rds")),
                              actionButton(ns("bmc_upload_load_btn"), "Load Uploaded List", icon = icon("play"), class = "btn-sm"),
                              uiOutput(ns("bmc_upload_results_ui")))
        ),
        div(class = "empty-note", style = "display:flex; align-items:center; justify-content:space-between; gap:12px; flex-wrap:wrap;",
            uiOutput(ns("bmc_selection_status_ui"), inline = TRUE),
            actionButton(ns("bmc_generate_btn"), "Generate Biomarker Card", icon = icon("id-card"), class = "btn-primary btn-sm"))
      )
    })

    bmc_picked_gene <- reactiveVal(NULL)
    has_card <- reactiveVal(FALSE)
    observeEvent(input$bmc_search_mode, bmc_picked_gene(NULL), ignoreInit = TRUE)

    output$bmc_selection_status_ui <- renderUI({
      mode <- input$bmc_search_mode %||% "gene"
      gene <- if (identical(mode, "gene")) trimws(input$bmc_gene_input %||% "") else bmc_picked_gene()
      if (!is.null(gene) && nzchar(gene)) {
        tagList(icon("circle-check", style = "color:#0ca30c;"), tags$b(sprintf("Selected: %s", gene)), " - click Generate to build the card.")
      } else {
        tagList(icon("circle-info"), "Select a biomarker above (click a table row, or type a gene symbol) before clicking Generate.")
      }
    })

    ## ---- Browse a saved Differential Expression run ----
    output$bmc_dge_results_ui <- renderUI({
      req(input$bmc_dge_run)
      tagList(p(class = "submodule-desc", "Click a row to select that gene, then click \"Generate Biomarker Card\"."),
              DT::dataTableOutput(ns("bmc_dge_table")))
    })
    output$bmc_dge_table <- DT::renderDataTable({
      req(input$bmc_dge_run)
      run <- (results$dge_runs %||% list())[[input$bmc_dge_run]]
      req(run)
      df <- run$table[order(run$table$adj.P.Val), , drop = FALSE]
      DT::datatable(utils::head(df, 500), colnames = c("Gene", "Log2 fold-change", "Adjusted p-value", "Direction"),
                    rownames = FALSE, selection = "single", options = list(pageLength = 10, scrollX = TRUE))
    })
    outputOptions(output, "bmc_dge_table", suspendWhenHidden = FALSE)
    observeEvent(input$bmc_dge_table_rows_selected, {
      run <- (results$dge_runs %||% list())[[input$bmc_dge_run]]; req(run)
      df <- run$table[order(run$table$adj.P.Val), , drop = FALSE]
      idx <- input$bmc_dge_table_rows_selected
      if (length(idx) == 1) bmc_picked_gene(utils::head(df, 500)$gene[idx])
    })

    ## ---- Browse Candidate Gene Identification output ----
    output$bmc_cand_results_ui <- renderUI({
      req(input$bmc_cand_set)
      tagList(p(class = "submodule-desc", "Click a row to select that gene, then click \"Generate Biomarker Card\"."),
              DT::dataTableOutput(ns("bmc_cand_table")))
    })
    output$bmc_cand_table <- DT::renderDataTable({
      cand <- results$candidates; req(cand); set <- cand[[input$bmc_cand_set %||% "final"]]; req(set)
      DT::datatable(data.frame(Gene = set$genes, stringsAsFactors = FALSE), rownames = FALSE, selection = "single", options = list(pageLength = 10, scrollX = TRUE))
    })
    outputOptions(output, "bmc_cand_table", suspendWhenHidden = FALSE)
    observeEvent(input$bmc_cand_table_rows_selected, {
      cand <- results$candidates; req(cand); set <- cand[[input$bmc_cand_set %||% "final"]]; req(set)
      idx <- input$bmc_cand_table_rows_selected
      if (length(idx) == 1) bmc_picked_gene(set$genes[idx])
    })

    ## ---- Browse a feature-selected consensus signature ----
    output$bmc_sig_results_ui <- renderUI({
      req(input$bmc_sig_sex)
      tagList(p(class = "submodule-desc", "Click a row to select that gene, then click \"Generate Biomarker Card\"."),
              DT::dataTableOutput(ns("bmc_sig_table")))
    })
    output$bmc_sig_table <- DT::renderDataTable({
      fs <- results$featureselection; req(fs); s <- fs[[input$bmc_sig_sex %||% ""]]; req(s)
      DT::datatable(data.frame(Gene = s$consensus_genes, stringsAsFactors = FALSE), rownames = FALSE, selection = "single", options = list(pageLength = 10, scrollX = TRUE))
    })
    outputOptions(output, "bmc_sig_table", suspendWhenHidden = FALSE)
    observeEvent(input$bmc_sig_table_rows_selected, {
      fs <- results$featureselection; req(fs); s <- fs[[input$bmc_sig_sex %||% ""]]; req(s)
      idx <- input$bmc_sig_table_rows_selected
      if (length(idx) == 1) bmc_picked_gene(s$consensus_genes[idx])
    })

    ## ---- Upload a gene list, or a Diagnostic Classifier RDS export (that
    ## tab's build_model_bundle() output: has $genes and $model_type) ----
    upload_table <- eventReactive(input$bmc_upload_load_btn, {
      validate(need(!is.null(input$bmc_upload_file), "Upload a .csv, .txt, or .rds file first."))
      path <- input$bmc_upload_file$datapath; name <- input$bmc_upload_file$name
      if (grepl("\\.rds$", name, ignore.case = TRUE)) {
        obj <- tryCatch(readRDS(path), error = function(e) NULL)
        validate(need(!is.null(obj) && !is.null(obj$genes) && !is.null(obj$model_type),
                      "Upload a Diagnostic Classifier RDS export (from that tab's own \"Save trained model\" download), or a plain .csv/.txt gene-symbol list."))
        df <- data.frame(gene = as.character(obj$genes), stringsAsFactors = FALSE)
      } else if (grepl("\\.txt$", name, ignore.case = TRUE)) {
        ids <- tryCatch(trimws(readLines(path)), error = function(e) character(0))
        ids <- unique(ids[nzchar(ids)])
        validate(need(length(ids) > 0, "No gene symbols were found in the uploaded file."))
        df <- data.frame(gene = ids, stringsAsFactors = FALSE)
      } else {
        up <- tryCatch(as.data.frame(data.table::fread(path, showProgress = FALSE)), error = function(e) NULL)
        validate(need(!is.null(up) && nrow(up) > 0, "Could not parse the uploaded file as a delimited table (CSV/TSV)."))
        gene_col <- intersect(c("gene", "Gene", "symbol", "Symbol", "GENE_SYMBOL", "ID"), colnames(up))[1]
        ids <- if (!is.na(gene_col)) as.character(up[[gene_col]]) else as.character(up[[1]])
        df <- data.frame(gene = unique(ids[nzchar(ids)]), stringsAsFactors = FALSE)
      }
      validate(need(nrow(df) > 0, "No gene symbols were found in the uploaded file."))
      df
    }, ignoreInit = TRUE)

    output$bmc_upload_results_ui <- renderUI({
      req(input$bmc_upload_load_btn)
      tagList(p(class = "submodule-desc", sprintf("%d gene(s) loaded from the uploaded file. Click a row to select it, then click \"Generate Biomarker Card\".", nrow(upload_table()))),
              DT::dataTableOutput(ns("bmc_upload_table")))
    })
    output$bmc_upload_table <- DT::renderDataTable({
      df <- upload_table(); req(df)
      DT::datatable(df, rownames = FALSE, selection = "single", options = list(pageLength = 10, scrollX = TRUE))
    })
    outputOptions(output, "bmc_upload_table", suspendWhenHidden = FALSE)
    observeEvent(input$bmc_upload_table_rows_selected, {
      df <- upload_table(); req(df)
      idx <- input$bmc_upload_table_rows_selected
      if (length(idx) == 1) bmc_picked_gene(df$gene[idx])
    })

    ## ---- Generate ----
    observeEvent(input$bmc_generate_btn, {
      has_card(TRUE)
      updateTabsetPanel(session, "bmc_subtabs", selected = "Biomarker Card")
    }, ignoreInit = TRUE)

    card_data <- eventReactive(input$bmc_generate_btn, {
      gene <- switch(input$bmc_search_mode,
        gene = trimws(input$bmc_gene_input %||% ""),
        dge = bmc_picked_gene(),
        candidates = bmc_picked_gene(),
        signature = bmc_picked_gene(),
        upload = bmc_picked_gene(),
        NULL
      )
      validate(need(!is.null(gene) && nzchar(gene),
                    "Select or enter a gene symbol before generating the Biomarker Card - type a gene symbol, or pick one from a saved Differential Expression run / the candidate list / a feature-selected signature / your uploaded list."))

      in_dataset <- tryCatch(gene %in% rownames(dataset$expr), error = function(e) FALSE)
      gene_identity <- tbc_gene_identity(gene)
      validate(need(in_dataset || isTRUE(gene_identity$ok),
                    sprintf("\"%s\" is neither present in the currently loaded expression matrix nor resolvable as a known human gene symbol - check the spelling.", gene)))

      live <- tbc_dataset_evidence(gene, dataset)
      dge_hits <- tbc_dge_matches(gene, results)
      selected_run_table <- if (identical(isolate(input$bmc_search_mode), "dge") && !is.null(isolate(input$bmc_dge_run))) {
        r <- (results$dge_runs %||% list())[[isolate(input$bmc_dge_run)]]
        if (!is.null(r)) r$table else NULL
      } else if (!is.null(results$dge_runs) && length(results$dge_runs) > 0) {
        utils::tail(results$dge_runs, 1)[[1]]$table
      } else NULL

      candidate_status <- tbc_candidate_status(gene, results)
      signature_membership <- tbc_signature_membership(gene, results)
      diagnostic_match <- tbc_diagnostic_lookup(gene, results)

      list(gene = gene, in_dataset = in_dataset, gene_identity = gene_identity, live = live,
           dge_hits = dge_hits, selected_run_table = selected_run_table,
           candidate_status = candidate_status, signature_membership = signature_membership,
           diagnostic_match = diagnostic_match)
    }, ignoreInit = TRUE)

    observeEvent(card_data(), {
      d <- card_data()
      best_adjp <- if (!is.null(d$dge_hits) && nrow(d$dge_hits) > 0) min(d$dge_hits$adj.P.Val, na.rm = TRUE) else NA_real_
      results$biomarkercard <- list(
        gene = d$gene, in_dataset = d$in_dataset,
        candidate = isTRUE(d$candidate_status$any),
        in_signature = length(d$signature_membership) > 0 && any(vapply(d$signature_membership, function(x) isTRUE(x$in_signature), logical(1))),
        best_adj_p = best_adjp
      )
    }, ignoreInit = TRUE)

    ## Biological Interpretation (GO/KEGG/Reactome) is opt-in (own
    ## button/reactive, not part of card_data()) so a slow/down external
    ## service never blocks the fully local, already-computed core card
    ## above - same rationale as the methylation card's own external
    ## evidence button.
    ext_data <- reactiveVal(NULL)
    observeEvent(card_data(), ext_data(NULL), ignoreInit = TRUE)
    observeEvent(input$bmc_ext_lookup_btn, {
      d <- card_data(); req(d)
      ext_data(tbc_pathway_evidence(d$gene, d$gene_identity$entrez %||% NA_character_))
    }, ignoreInit = TRUE)

    output$bmc_card_ui <- renderUI({
      if (!has_card()) return(div(class = "empty-note", icon("circle-info"), "Select a biomarker and click \"Generate Biomarker Card\" on the \"Select Biomarker\" tab."))
      d <- card_data(); req(d)
      tagList(
        tbc_section_session_summary(results, ext_data()),
        tbc_section_identity(d),
        tbc_section_discovery_context(dataset, d$live),
        tbc_section_expression_data(dataset, d$live),
        tbc_section_differential_expression(d),
        div(class = "card",
            div(class = "card-title", icon("chart-column"), "Differential Expression (Volcano)"),
            p(class = "submodule-desc", "This gene highlighted within the most recently selected/used Differential Expression run this session."),
            if (is.null(d$selected_run_table)) div(class = "empty-note", icon("circle-info"), "No Differential Expression run is available - run the Differential Expression tab first.")
            else withSpinner(plotOutput(ns("bmc_volcano_plot"), height = "380px"), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("chart-simple"), "Expression Distribution (Your Dataset)"),
            p(class = "submodule-desc", "Source: the currently loaded expression matrix, grouped by the detected case/control column."),
            withSpinner(plotly::plotlyOutput(ns("bmc_dist_plot"), height = "300px"), color = "#2563EB", type = 6)
        ),
        div(class = "card",
            div(class = "card-title", icon("venus-mars"), "Sex-Specific Expression"),
            p(class = "submodule-desc", "Same comparison, split by the detected sex column (if present in the loaded sample metadata)."),
            withSpinner(plotly::plotlyOutput(ns("bmc_sex_dist_plot"), height = "320px"), color = "#2563EB", type = 6)
        ),
        tbc_section_prioritization(d),
        tbc_section_signature(d),
        tbc_section_performance(d),
        tbc_section_discovery_validation(d, dataset),
        div(class = "card",
            div(class = "card-title", icon("globe"), "Biological Interpretation"),
            p(class = "submodule-desc", "Live queries to Gene Ontology (GO.db), KEGG (KEGGREST), and Reactome - opt-in, since these are external network calls independent of everything above."),
            actionButton(ns("bmc_ext_lookup_btn"), "Look Up Biological Interpretation", icon = icon("magnifying-glass-location"), class = "btn-primary btn-sm")
        ),
        withSpinner(uiOutput(ns("bmc_ext_sections_ui")), color = "#2563EB", type = 6),
        tbc_section_interpretation(d, ext_data()),
        tbc_section_evidence_summary(d, ext_data()),
        tbc_section_sources(ext_data()),
        div(class = "card", div(class = "card-title", icon("download"), "Download"),
            p(class = "submodule-desc", "A self-contained HTML report covering every section above (including Biological Interpretation, if you've looked it up)."),
            downloadButton(ns("bmc_download_report"), "Download Biomarker Report (HTML)", class = "btn-primary btn-sm"))
      )
    })

    output$bmc_ext_sections_ui <- renderUI({ tbc_section_pathways(ext_data()) })

    output$bmc_volcano_plot <- renderPlot({
      d <- card_data(); req(d)
      validate(need(!is.null(d$selected_run_table), "No Differential Expression run is available to plot."))
      tbc_plot_volcano_highlight(d$selected_run_table, d$gene)
    })

    output$bmc_dist_plot <- plotly::renderPlotly({
      d <- card_data(); req(d)
      validate(need(isTRUE(d$live$ok) && isTRUE(d$live$overall$ok), "No live dataset evidence is available to plot for this gene."))
      df <- data.frame(expr = d$live$values, group = d$live$group_vec, stringsAsFactors = FALSE)
      y_lab <- if (isTRUE(d$live$transformed)) "log2(CPM + 1)" else "Expression (analysis scale)"
      plotly::ggplotly(tbc_plot_expression_dist(df, y_lab))
    })

    output$bmc_sex_dist_plot <- plotly::renderPlotly({
      d <- card_data(); req(d)
      validate(need(isTRUE(d$live$ok) && !is.null(d$live$sex_vec), "No sex information is available in the loaded sample metadata to compare female vs male expression."))
      df <- data.frame(expr = d$live$values, group = d$live$group_vec, sex = d$live$sex_vec, stringsAsFactors = FALSE)
      y_lab <- if (isTRUE(d$live$transformed)) "log2(CPM + 1)" else "Expression (analysis scale)"
      plotly::ggplotly(tbc_plot_expression_dist(df, y_lab, facet_sex = TRUE))
    })

    output$bmc_download_report <- downloadHandler(
      filename = function() sprintf("transcriptomic_biomarker_card_%s.html", card_data()$gene),
      content = function(file) {
        d <- card_data(); req(d)
        body <- tbc_build_report_tags(d, dataset, results, ext_data())
        page <- tags$html(
          tags$head(tags$meta(charset = "utf-8"), tags$title(sprintf("Transcriptomic Biomarker Card %s", d$gene)), tags$style(tbc_report_css())),
          tags$body(body)
        )
        htmltools::save_html(page, file = file)
      }
    )
  })
}
