## R/methylomics/15_Biomarker_Analysis/mod_methyl_biomarkercard.R
## Submodule: Biomarker Card - single-CpG epigenomic profile page.
##

.bc_illumina_anno_cache <- new.env(parent = emptyenv())

bc_illumina_full_annotation <- function(array_type) {
  pkg <- METHYL_ANNOTATION_PACKAGES[[array_type]]
  if (is.null(pkg)) {
    return(list(ok = FALSE, anno = NULL, reason = sprintf(
      "No Bioconductor manifest annotation is available for %s in this deployment.", array_type)))
  }
  cached <- .bc_illumina_anno_cache[[pkg]]
  if (!is.null(cached)) return(list(ok = TRUE, anno = cached, reason = NULL))
  if (!requireNamespace(pkg, quietly = TRUE)) {
    return(list(ok = FALSE, anno = NULL, reason = sprintf("The %s annotation package is not installed in this deployment.", pkg)))
  }
  anno <- tryCatch({
    e <- new.env(parent = emptyenv())
    utils::data(list = c("Locations", "Other", "Islands.UCSC"), package = pkg, envir = e)
    loc <- as.data.frame(e$Locations)
    oth <- as.data.frame(e$Other)
    isl <- as.data.frame(e$Islands.UCSC)
    ids <- rownames(loc)
    enhancer <- if ("Enhancer" %in% colnames(oth)) ifelse(oth[ids, "Enhancer"] == "TRUE", "TRUE", NA_character_) else rep(NA_character_, length(ids))
    dhs <- if ("DHS" %in% colnames(oth)) ifelse(oth[ids, "DHS"] == "TRUE", "TRUE", NA_character_) else rep(NA_character_, length(ids))
    data.frame(
      chr = loc[ids, "chr"], pos = loc[ids, "pos"], strand = loc[ids, "strand"],
      gene_names = oth[ids, "UCSC_RefGene_Name"], gene_group = oth[ids, "UCSC_RefGene_Group"],
      regulatory_group = oth[ids, "Regulatory_Feature_Group"], regulatory_name = oth[ids, "Regulatory_Feature_Name"],
      enhancer = enhancer, dhs = dhs,
      island_name = isl[ids, "Islands_Name"], island_relation = isl[ids, "Relation_to_Island"],
      row.names = ids, stringsAsFactors = FALSE
    )
  }, error = function(e) e)
  if (inherits(anno, "error")) {
    return(list(ok = FALSE, anno = NULL, reason = sprintf("Could not extract annotation from %s: %s", pkg, conditionMessage(anno))))
  }
  .bc_illumina_anno_cache[[pkg]] <- anno
  list(ok = TRUE, anno = anno, reason = NULL)
}

.bc_champ_anno_cache <- new.env(parent = emptyenv())
bc_champ_annotation <- function() {
  if (!is.null(.bc_champ_anno_cache$anno)) return(.bc_champ_anno_cache$anno)
  if (!requireNamespace("ChAMPdata", quietly = TRUE)) return(NULL)
  e <- new.env(parent = emptyenv())
  ok <- tryCatch({ utils::data("probe.features", package = "ChAMPdata", envir = e); TRUE }, error = function(e) FALSE)
  if (!ok || is.null(e$probe.features)) return(NULL)
  pf <- e$probe.features
  out <- data.frame(
    cpg = rownames(pf), chr = paste0("chr", pf$CHR), pos = pf$MAPINFO,
    gene = as.character(pf$gene), feature = as.character(pf$feature), island = as.character(pf$cgi),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$chr) & !is.na(out$pos), , drop = FALSE]
  .bc_champ_anno_cache$anno <- out
  out
}

BC_GROUP_PATTERNS <- c("^group$", "^status$", "^phenotype$", "^disease$", "^condition$", "^diagnosis$", "group")
BC_SEX_PATTERNS <- c("^sex$", "^gender$", "sex")
BC_ID_PATTERNS <- c("^sample$", "^gsm$", "^sample_id$", "^id$")

bc_find_col <- function(cols, patterns) {
  for (p in patterns) {
    hit <- cols[grepl(p, cols, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  NULL
}

bc_island_parse <- function(islands_name) {
  if (is.null(islands_name) || is.na(islands_name) || !nzchar(islands_name)) return(NULL)
  m <- regmatches(islands_name, regexec("^(chr[^:]+):([0-9]+)-([0-9]+)$", islands_name))[[1]]
  if (length(m) != 4) return(NULL)
  start <- as.numeric(m[3]); end <- as.numeric(m[4])
  list(chr = m[2], start = start, end = end, length = end - start + 1)
}

BC_FEATURE_MAP <- c(
  "TSS200" = "Promoter/TSS", "TSS1500" = "Promoter/TSS", "5'UTR" = "5' UTR",
  "1stExon" = "1st Exon", "Body" = "Gene body", "3'UTR" = "3' UTR", "ExonBnd" = "Exon boundary"
)
bc_feature_label <- function(gene_group) {
  if (is.null(gene_group) || is.na(gene_group) || !nzchar(gene_group)) return("Intergenic (no annotated RefGene overlap)")
  toks <- unique(trimws(strsplit(gene_group, ";")[[1]]))
  labs <- unique(stats::na.omit(unname(BC_FEATURE_MAP[toks])))
  if (length(labs) == 0) return(gene_group)
  paste(labs, collapse = " / ")
}

bc_island_context_label <- function(relation) {
  if (is.null(relation) || is.na(relation) || !nzchar(relation)) return("Open Sea (no annotated CpG island within range)")
  switch(relation,
    "Island" = "CpG Island", "N_Shore" = "Shore (North)", "S_Shore" = "Shore (South)",
    "N_Shelf" = "Shelf (North)", "S_Shelf" = "Shelf (South)", "OpenSea" = "Open Sea",
    relation)
}

bc_enhancer_flag <- function(r, array_type) {
  if (!identical(array_type, "450K")) return("Not available for this array type in this build")
  if (is.na(r$enhancer)) return("No annotated enhancer overlap")
  "Yes - annotated enhancer region (450K manifest)"
}
bc_dhs_flag <- function(r, array_type) {
  if (!identical(array_type, "450K")) return("Not available for this array type in this build")
  if (is.na(r$dhs)) return("No annotated DHS overlap")
  "Yes - annotated DNase hypersensitive site (450K manifest)"
}

bc_resolve_cpg <- function(cpg_id, array_type) {
  a <- bc_illumina_full_annotation(array_type)
  champ <- bc_champ_annotation()
  illumina_row <- if (isTRUE(a$ok) && cpg_id %in% rownames(a$anno)) a$anno[cpg_id, ] else NULL
  champ_row <- if (!is.null(champ) && cpg_id %in% champ$cpg) champ[champ$cpg == cpg_id, , drop = FALSE][1, ] else NULL
  list(
    cpg = cpg_id, found_in_manifest = !is.null(illumina_row), found_in_champ = !is.null(champ_row),
    chr = if (!is.null(illumina_row)) illumina_row$chr else if (!is.null(champ_row)) champ_row$chr else NA_character_,
    pos = if (!is.null(illumina_row)) illumina_row$pos else if (!is.null(champ_row)) champ_row$pos else NA_real_,
    strand = if (!is.null(illumina_row)) illumina_row$strand else NA_character_,
    gene_names = if (!is.null(illumina_row)) illumina_row$gene_names else NA_character_,
    gene_group = if (!is.null(illumina_row)) illumina_row$gene_group else NA_character_,
    regulatory_group = if (!is.null(illumina_row)) illumina_row$regulatory_group else NA_character_,
    regulatory_name = if (!is.null(illumina_row)) illumina_row$regulatory_name else NA_character_,
    enhancer = if (!is.null(illumina_row)) illumina_row$enhancer else NA_character_,
    dhs = if (!is.null(illumina_row)) illumina_row$dhs else NA_character_,
    island_name = if (!is.null(illumina_row)) illumina_row$island_name else NA_character_,
    island_relation = if (!is.null(illumina_row)) illumina_row$island_relation else if (!is.null(champ_row)) champ_row$island else NA_character_,
    champ_gene = if (!is.null(champ_row)) champ_row$gene else NA_character_,
    champ_feature = if (!is.null(champ_row)) champ_row$feature else NA_character_
  )
}

bc_search_gene <- function(symbol, array_type, max_n = 300) {
  a <- bc_illumina_full_annotation(array_type)
  if (!isTRUE(a$ok)) return(list(ok = FALSE, reason = a$reason, df = NULL))
  anno <- a$anno
  esc <- gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", symbol)
  pattern <- sprintf("(^|;)%s(;|$)", esc)
  hit <- !is.na(anno$gene_names) & grepl(pattern, anno$gene_names, ignore.case = TRUE)
  ids <- rownames(anno)[hit]
  if (length(ids) == 0) return(list(ok = TRUE, df = data.frame(cpg = character(0), chr = character(0), pos = numeric(0),
                                                                 gene_names = character(0), gene_group = character(0),
                                                                 island_relation = character(0), stringsAsFactors = FALSE)))
  out <- anno[ids, c("chr", "pos", "gene_names", "gene_group", "island_relation"), drop = FALSE]
  out$cpg <- ids
  out <- out[order(out$pos), , drop = FALSE]
  if (nrow(out) > max_n) out <- out[seq_len(max_n), , drop = FALSE]
  list(ok = TRUE, df = out)
}

bc_gene_relationship_table <- function(resolved) {
  names_str <- resolved$gene_names; group_str <- resolved$gene_group
  if (is.null(names_str) || is.na(names_str) || !nzchar(names_str)) return(NULL)
  nm <- trimws(strsplit(names_str, ";")[[1]])
  gp <- if (!is.null(group_str) && !is.na(group_str) && nzchar(group_str)) trimws(strsplit(group_str, ";")[[1]]) else character(0)
  n <- length(nm)
  if (length(gp) < n) gp <- c(gp, rep(NA_character_, n - length(gp)))
  df <- data.frame(
    Gene = nm,
    Relationship = vapply(gp[seq_len(n)], function(x) if (is.na(x) || !nzchar(x)) "Overlapping (unspecified)" else bc_feature_label(x), character(1)),
    stringsAsFactors = FALSE
  )
  df[!duplicated(df$Gene), , drop = FALSE]
}

.bc_txdb_genes_cache <- new.env(parent = emptyenv())
bc_txdb_genes <- function() {
  if (!is.null(.bc_txdb_genes_cache$g)) return(.bc_txdb_genes_cache$g)
  if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE) || !requireNamespace("GenomicFeatures", quietly = TRUE)) return(NULL)
  g <- tryCatch(suppressWarnings(suppressMessages(
    GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
  )), error = function(e) NULL)
  .bc_txdb_genes_cache$g <- g
  g
}

.bc_txdb_exons_cache <- new.env(parent = emptyenv())
bc_txdb_exons_by_gene <- function() {
  if (!is.null(.bc_txdb_exons_cache$ex)) return(.bc_txdb_exons_cache$ex)
  if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE) || !requireNamespace("GenomicFeatures", quietly = TRUE)) return(NULL)
  ex <- tryCatch(GenomicFeatures::exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene"), error = function(e) NULL)
  .bc_txdb_exons_cache$ex <- ex
  ex
}

bc_gene_structure <- function(symbol) {
  if (is.null(symbol) || is.na(symbol) || !nzchar(symbol)) return(list(ok = FALSE, reason = "No associated gene symbol to resolve."))
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) return(list(ok = FALSE, reason = "org.Hs.eg.db is not installed in this deployment."))
  res <- tryCatch({
    map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = symbol, keytype = "SYMBOL",
                                                    columns = c("ENTREZID", "ENSEMBL", "GENENAME")))
    map <- map[!is.na(map$ENTREZID), , drop = FALSE]
    if (nrow(map) == 0) return(list(ok = FALSE, reason = sprintf("No NCBI Entrez Gene entry found for symbol \"%s\".", symbol)))
    entrez <- map$ENTREZID[1]
    g <- bc_txdb_genes()
    gene_span <- if (!is.null(g) && entrez %in% names(g)) g[entrez] else NULL
    ex_all <- bc_txdb_exons_by_gene()
    ex <- if (!is.null(ex_all) && entrez %in% names(ex_all)) ex_all[[entrez]] else NULL
    list(ok = TRUE, symbol = symbol, entrez = entrez, ensembl = map$ENSEMBL[1], genename = map$GENENAME[1],
         chr = if (!is.null(gene_span)) as.character(GenomicRanges::seqnames(gene_span))[1] else NA_character_,
         start = if (!is.null(gene_span)) BiocGenerics::start(gene_span)[1] else NA_real_,
         end = if (!is.null(gene_span)) BiocGenerics::end(gene_span)[1] else NA_real_,
         strand = if (!is.null(gene_span)) as.character(BiocGenerics::strand(gene_span))[1] else NA_character_,
         exons = ex)
  }, error = function(e) e)
  if (inherits(res, "error")) return(list(ok = FALSE, reason = sprintf("Could not resolve gene structure: %s", conditionMessage(res))))
  res
}

bc_go_terms <- function(entrez, n = 8) {
  if (is.null(entrez) || is.na(entrez) || !nzchar(entrez)) return(NULL)
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) || !requireNamespace("GO.db", quietly = TRUE)) return(NULL)
  res <- tryCatch({
    go <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = entrez, keytype = "ENTREZID", columns = c("GOALL", "ONTOLOGYALL")))
    go <- go[!is.na(go$GOALL) & go$ONTOLOGYALL == "BP", , drop = FALSE]
    ids <- unique(go$GOALL)
    if (length(ids) == 0) return(NULL)
    terms <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = ids, keytype = "GOID", columns = "TERM"))
    utils::head(unique(stats::na.omit(terms[, c("GOID", "TERM")])), n)
  }, error = function(e) NULL)
  if (inherits(res, "error")) return(NULL)
  res
}

BC_CYTOBAND_PATH <- get_reference_path("cytoBandIdeo_hg19.txt.gz")
.bc_cytoband_cache <- new.env(parent = emptyenv())
bc_cytoband_hg19 <- function() {
  if (!is.null(.bc_cytoband_cache$df)) return(.bc_cytoband_cache$df)
  if (!file.exists(BC_CYTOBAND_PATH)) return(NULL)
  df <- tryCatch(
    utils::read.delim(gzfile(BC_CYTOBAND_PATH), header = FALSE,
                       col.names = c("chrom", "chromStart", "chromEnd", "name", "gieStain"),
                       stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  .bc_cytoband_cache$df <- df
  df
}

BC_GIESTAIN_COLORS <- c(
  gneg = "#FFFFFF", gpos25 = "#C8C8C8", gpos50 = "#969696", gpos75 = "#5A5A5A", gpos100 = "#141414",
  acen = "#B33A3A", gvar = "#A9A9A9", stalk = "#6E8FB0"
)

.bc_dmp_cache <- new.env(parent = emptyenv())
bc_dmp_table <- function(sex, stage = "sva") {
  key <- paste(stage, sex, sep = "_")
  cached <- .bc_dmp_cache[[key]]
  if (!is.null(cached)) return(cached)
  df <- load_default_dmp(stage, sex)
  if (!is.null(df)) .bc_dmp_cache[[key]] <- df
  df
}
.bc_dmr_cache <- new.env(parent = emptyenv())
bc_dmr_table <- function(sex) {
  cached <- .bc_dmr_cache[[sex]]
  if (!is.null(cached)) return(cached)
  df <- load_default_dmr(sex)
  if (!is.null(df)) .bc_dmr_cache[[sex]] <- df
  df
}
bc_dmp_lookup <- function(cpg, sex, stage = "sva") {
  df <- bc_dmp_table(sex, stage)
  if (is.null(df)) return(NULL)
  row <- df[df$cpg == cpg, , drop = FALSE]
  if (nrow(row) == 0) return(NULL)
  row[1, ]
}
bc_dmr_lookup <- function(chr, pos, sex) {
  df <- bc_dmr_table(sex)
  if (is.null(df) || is.null(chr) || is.na(chr) || is.null(pos) || is.na(pos)) return(NULL)
  chr_norm <- sub("^chr", "", chr, ignore.case = TRUE)
  seq_norm <- sub("^chr", "", as.character(df$seqnames), ignore.case = TRUE)
  hit <- df[seq_norm == chr_norm & df$start <= pos & df$end >= pos, , drop = FALSE]
  if (nrow(hit) == 0) return(NULL)
  hit
}

bc_pick_case_control <- function(levels) {
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

bc_live_stats <- function(beta_row, group_vec, case_label, control_label) {
  keep <- !is.na(beta_row) & !is.na(group_vec) & group_vec %in% c(case_label, control_label)
  beta_row <- beta_row[keep]; grp <- group_vec[keep]
  case_vals <- beta_row[grp == case_label]; ctrl_vals <- beta_row[grp == control_label]
  if (length(case_vals) < 2 || length(ctrl_vals) < 2) {
    return(list(ok = FALSE, reason = "Fewer than 2 non-missing samples in one of the two groups - cannot compute a t-test.",
                case_label = case_label, control_label = control_label, n_case = length(case_vals), n_control = length(ctrl_vals)))
  }
  tt <- tryCatch(stats::t.test(case_vals, ctrl_vals), error = function(e) NULL)
  mc <- mean(case_vals); mo <- mean(ctrl_vals)
  list(ok = TRUE, mean_case = mc, mean_control = mo, dbeta = mc - mo,
       p_value = if (!is.null(tt)) tt$p.value else NA_real_,
       ci = if (!is.null(tt)) as.numeric(tt$conf.int) else c(NA_real_, NA_real_),
       n_case = length(case_vals), n_control = length(ctrl_vals),
       direction = if (mc > mo) "Hypermethylated" else if (mc < mo) "Hypomethylated" else "No change",
       case_label = case_label, control_label = control_label)
}

bc_dataset_evidence <- function(cpg, dataset) {
  beta <- dataset$beta; sheet <- dataset$sample_sheet
  if (is.null(beta) || is.null(sheet)) return(list(ok = FALSE, reason = "No beta matrix / sample sheet is currently loaded on the Dataset tab."))
  if (!cpg %in% rownames(beta)) return(list(ok = FALSE, reason = "This CpG is not present in the currently loaded beta matrix."))
  samp_ids <- colnames(beta)
  id_col <- bc_find_col(colnames(sheet), BC_ID_PATTERNS)
  if (!is.null(id_col) && all(samp_ids %in% as.character(sheet[[id_col]]))) {
    sheet <- sheet[match(samp_ids, as.character(sheet[[id_col]])), , drop = FALSE]
  } else if (nrow(sheet) != length(samp_ids)) {
    return(list(ok = FALSE, reason = "Could not align the sample sheet to the beta matrix columns (no shared sample-ID column and row counts differ)."))
  }
  beta_row <- beta[cpg, ]
  group_col <- bc_find_col(colnames(sheet), BC_GROUP_PATTERNS)
  if (is.null(group_col)) return(list(ok = FALSE, reason = "Could not auto-detect a case/control grouping column in the loaded sample sheet.", beta_row = beta_row))
  group_vec <- as.character(sheet[[group_col]])
  levels_present <- unique(group_vec[!is.na(group_vec)])
  if (length(levels_present) != 2) {
    return(list(ok = FALSE, reason = sprintf("The detected grouping column (\"%s\") does not have exactly two levels (found: %s).",
                                              group_col, paste(levels_present, collapse = ", ")), beta_row = beta_row))
  }
  cc <- bc_pick_case_control(levels_present)
  overall <- bc_live_stats(beta_row, group_vec, cc$case, cc$control)
  sex_col <- bc_find_col(colnames(sheet), BC_SEX_PATTERNS)
  sex_vec <- if (!is.null(sex_col)) as.character(sheet[[sex_col]]) else NULL
  by_sex <- NULL
  if (!is.null(sex_vec)) {
    sex_levels <- unique(sex_vec[!is.na(sex_vec)])
    if (length(sex_levels) > 0) {
      by_sex <- stats::setNames(lapply(sex_levels, function(s) {
        idx <- which(sex_vec == s)
        bc_live_stats(beta_row[idx], group_vec[idx], cc$case, cc$control)
      }), sex_levels)
    }
  }
  list(ok = TRUE, beta_row = beta_row, group_vec = group_vec, sex_vec = sex_vec,
       group_col = group_col, sex_col = sex_col, case_label = cc$case, control_label = cc$control,
       overall = overall, by_sex = by_sex)
}

bc_meta <- function(source, query, endpoint, status, n_records = NA_integer_) {
  list(source = source, query = query, endpoint = endpoint, retrieved_at = Sys.time(), status = status, n_records = n_records)
}

bc_ewascatalog_query <- function(type, value) {
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, df = NULL, reason = "httr2 is not installed.", meta = NULL))
  url <- sprintf("https://www.ewascatalog.org/api/?%s=%s", type, utils::URLencode(value, reserved = TRUE))
  meta_base <- bc_meta("EWAS Catalog", sprintf("%s=%s", type, value), "www.ewascatalog.org/api/", NA_character_)
  res <- tryCatch({
    resp <- httr2::request(url) %>% httr2::req_timeout(15) %>% httr2::req_perform()
    httr2::resp_body_json(resp, simplifyVector = FALSE)
  }, error = function(e) e)
  if (inherits(res, "error")) {
    meta_base$status <- "Failed"
    return(list(ok = FALSE, df = NULL, reason = sprintf("EWAS Catalog lookup failed: %s", conditionMessage(res)), meta = meta_base))
  }
  fields <- unlist(res$fields)
  rows <- res$results
  if (is.null(fields)) {
    meta_base$status <- "Failed (no field schema)"
    return(list(ok = FALSE, df = NULL, reason = "EWAS Catalog response did not include a field schema.", meta = meta_base))
  }
  if (length(rows) == 0) {
    empty <- as.data.frame(matrix(character(0), ncol = length(fields)), stringsAsFactors = FALSE)
    colnames(empty) <- fields
    meta_base$status <- "Success"; meta_base$n_records <- 0L
    return(list(ok = TRUE, df = empty, reason = NULL, meta = meta_base))
  }
  mat <- do.call(rbind, lapply(rows, function(r) {
    r <- lapply(r, function(x) if (is.null(x)) NA_character_ else as.character(x))
    unlist(r)
  }))
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(df) <- fields
  rownames(df) <- NULL
  meta_base$status <- "Success"; meta_base$n_records <- nrow(df)
  list(ok = TRUE, df = df, reason = NULL, meta = meta_base)
}

bc_ewasatlas_query <- function(endpoint, param, value) {
  meta_base <- bc_meta("EWAS Atlas", sprintf("%s=%s", param, value), sprintf("ngdc.cncb.ac.cn/ewas/rest/%s", endpoint), NA_character_)
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, data = NULL, reason = "httr2 is not installed.", meta = NULL))
  url <- sprintf("https://ngdc.cncb.ac.cn/ewas/rest/%s?%s=%s", endpoint, param, utils::URLencode(value, reserved = TRUE))
  res <- tryCatch({
    resp <- httr2::request(url) %>% httr2::req_timeout(15) %>% httr2::req_perform()
    httr2::resp_body_json(resp, simplifyVector = TRUE)
  }, error = function(e) e)
  if (inherits(res, "error")) {
    meta_base$status <- "Failed"
    return(list(ok = FALSE, data = NULL, reason = sprintf("EWAS Atlas lookup failed: %s", conditionMessage(res)), meta = meta_base))
  }
  code <- res$code
  if (!is.null(code) && !identical(as.integer(code), 0L)) {
    meta_base$status <- "Failed (API error)"
    return(list(ok = FALSE, data = NULL, reason = res$msg %||% "EWAS Atlas returned an error.", meta = meta_base))
  }
  meta_base$status <- "Success"
  meta_base$n_records <- tryCatch(if (is.data.frame(res$data)) nrow(res$data) else length(res$data), error = function(e) NA_integer_)
  list(ok = TRUE, data = res$data, reason = NULL, meta = meta_base)
}
bc_ewasatlas_probe <- function(probe_id) bc_ewasatlas_query("probe", "probeId", probe_id)

.bc_kegg_pathway_names_cache <- new.env(parent = emptyenv())
bc_kegg_pathway_names <- function() {
  if (!is.null(.bc_kegg_pathway_names_cache$v)) return(.bc_kegg_pathway_names_cache$v)
  if (!requireNamespace("KEGGREST", quietly = TRUE)) return(NULL)
  v <- tryCatch(KEGGREST::keggList("pathway", "hsa"), error = function(e) NULL)
  .bc_kegg_pathway_names_cache$v <- v
  v
}
bc_kegg_pathways_for_gene <- function(entrez) {
  meta_base <- bc_meta("KEGG", sprintf("hsa:%s", entrez %||% NA), "rest.kegg.jp/link/pathway/hsa:{entrez} (via KEGGREST)", NA_character_)
  if (is.null(entrez) || is.na(entrez) || !nzchar(entrez)) return(list(ok = FALSE, pathways = NULL, reason = "No NCBI Gene ID available for KEGG lookup.", meta = meta_base))
  if (!requireNamespace("KEGGREST", quietly = TRUE)) return(list(ok = FALSE, pathways = NULL, reason = "KEGGREST is not installed in this deployment.", meta = meta_base))
  res <- tryCatch({
    links <- KEGGREST::keggLink("pathway", sprintf("hsa:%s", entrez))
    key <- sub("^path:", "", unname(links))
    names_map <- bc_kegg_pathway_names()
    nm <- if (!is.null(names_map)) unname(names_map[key]) else rep(NA_character_, length(key))
    nm <- sub(" - Homo sapiens.*$", "", nm)
    data.frame(id = key, name = nm, stringsAsFactors = FALSE)
  }, error = function(e) e)
  if (inherits(res, "error")) {
    meta_base$status <- "Failed"
    return(list(ok = FALSE, pathways = NULL, reason = sprintf("KEGG lookup failed: %s", conditionMessage(res)), meta = meta_base))
  }
  meta_base$status <- "Success"; meta_base$n_records <- nrow(res)
  list(ok = TRUE, pathways = res, reason = NULL, meta = meta_base)
}

.bc_kegg_ra_pathway_cache <- new.env(parent = emptyenv())
bc_kegg_ra_pathway_genes <- function() {
  if (!is.null(.bc_kegg_ra_pathway_cache$obj)) return(.bc_kegg_ra_pathway_cache$obj)
  if (!requireNamespace("KEGGREST", quietly = TRUE)) return(NULL)
  res <- tryCatch({
    p <- KEGGREST::keggGet("hsa05323")[[1]]
    genes <- p$GENE
    ids <- genes[seq(1, length(genes), by = 2)]
    descs <- genes[seq(2, length(genes), by = 2)]
    syms <- sub(";.*", "", descs)
    list(df = data.frame(entrez = ids, symbol = syms, stringsAsFactors = FALSE),
         name = paste(p$NAME, collapse = " "), description = paste(p$DESCRIPTION, collapse = " "))
  }, error = function(e) e)
  if (inherits(res, "error")) return(NULL)
  .bc_kegg_ra_pathway_cache$obj <- res
  res
}
bc_kegg_ra_pathway_check <- function(entrez, symbol) {
  ra <- bc_kegg_ra_pathway_genes()
  if (is.null(ra)) return(list(ok = FALSE, reason = "KEGG RA pathway lookup failed or is unavailable in this deployment."))
  in_by_id <- !is.null(entrez) && !is.na(entrez) && entrez %in% ra$df$entrez
  in_by_symbol <- !is.null(symbol) && !is.na(symbol) && toupper(symbol) %in% toupper(ra$df$symbol)
  list(ok = TRUE, in_pathway = in_by_id || in_by_symbol, pathway_name = ra$name, pathway_description = ra$description, genes = ra$df)
}

bc_reactome_pathways_for_gene <- function(symbol) {
  meta_base <- bc_meta("Reactome", symbol %||% NA, "reactome.org/ContentService/data/mapping/UniProt/{id}/pathways", NA_character_)
  if (is.null(symbol) || is.na(symbol) || !nzchar(symbol)) return(list(ok = FALSE, pathways = NULL, reason = "No gene symbol available for Reactome lookup.", meta = meta_base))
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) || !requireNamespace("httr2", quietly = TRUE)) {
    return(list(ok = FALSE, pathways = NULL, reason = "org.Hs.eg.db / httr2 not available in this deployment.", meta = meta_base))
  }
  uniprot_ids <- tryCatch({
    m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = symbol, keytype = "SYMBOL", columns = "UNIPROT"))
    unique(stats::na.omit(m$UNIPROT))
  }, error = function(e) character(0))
  if (length(uniprot_ids) == 0) {
    meta_base$status <- "Failed (no UniProt ID)"
    return(list(ok = FALSE, pathways = NULL, reason = sprintf("No UniProt ID found for gene symbol \"%s\".", symbol), meta = meta_base))
  }
  for (uid in uniprot_ids) {
    res <- tryCatch({
      url <- sprintf("https://reactome.org/ContentService/data/mapping/UniProt/%s/pathways?species=9606", uid)
      resp <- httr2::request(url) %>% httr2::req_timeout(15) %>%
        httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
      if (httr2::resp_status(resp) != 200) NULL else httr2::resp_body_json(resp, simplifyVector = TRUE)
    }, error = function(e) NULL)
    if (!is.null(res) && is.data.frame(res) && nrow(res) > 0) {
      out <- res[, intersect(c("stId", "displayName"), colnames(res)), drop = FALSE]
      meta_base$status <- "Success"; meta_base$n_records <- nrow(out)
      return(list(ok = TRUE, pathways = out, reason = NULL, uniprot = uid, meta = meta_base))
    }
  }
  meta_base$status <- "Success"; meta_base$n_records <- 0L
  list(ok = TRUE, pathways = data.frame(stId = character(0), displayName = character(0)), reason = NULL, meta = meta_base)
}

bc_pubmed_summaries <- function(pmids) {
  pmids <- unique(pmids[!is.na(pmids) & nzchar(pmids)])
  empty <- data.frame(pmid = character(0), title = character(0), journal = character(0), year = character(0), stringsAsFactors = FALSE)
  meta_base <- bc_meta("PubMed (NCBI E-utilities)", paste(pmids, collapse = ","), "eutils.ncbi.nlm.nih.gov/esummary.fcgi?db=pubmed", NA_character_)
  if (length(pmids) == 0) { meta_base$status <- "Success"; meta_base$n_records <- 0L; return(list(ok = TRUE, table = empty, reason = NULL, meta = meta_base)) }
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, table = NULL, reason = "httr2 is not installed.", meta = NULL))
  pmids <- utils::head(pmids, 20)
  res <- tryCatch({
    httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi") %>%
      httr2::req_url_query(db = "pubmed", id = paste(pmids, collapse = ","), retmode = "json", tool = "arthomix-explorer") %>%
      httr2::req_timeout(15) %>%
      httr2::req_perform() %>%
      httr2::resp_body_json()
  }, error = function(e) e)
  if (inherits(res, "error")) {
    meta_base$status <- "Failed"
    return(list(ok = FALSE, table = NULL, reason = sprintf("PubMed lookup failed: %s", conditionMessage(res)), meta = meta_base))
  }
  ids <- unlist(res$result$uids)
  if (is.null(ids)) { meta_base$status <- "Success"; meta_base$n_records <- 0L; return(list(ok = TRUE, table = empty, reason = NULL, meta = meta_base)) }
  rows <- lapply(ids, function(id) {
    rec <- res$result[[id]]
    data.frame(pmid = id, title = rec$title %||% NA_character_,
               journal = rec$fulljournalname %||% rec$source %||% NA_character_,
               year = substr(rec$pubdate %||% "", 1, 4), stringsAsFactors = FALSE)
  })
  tbl <- do.call(rbind, rows)
  meta_base$status <- "Success"; meta_base$n_records <- nrow(tbl)
  list(ok = TRUE, table = tbl, reason = NULL, meta = meta_base)
}

bc_methbank_link <- function() "https://ngdc.cncb.ac.cn/methbank/"

BC_DISEASE_CATEGORIES <- list(
  "Rheumatoid Arthritis" = "rheumatoid arthritis|\\bra\\b",
  "Osteoarthritis" = "osteoarthritis|\\boa\\b",
  "SLE / Lupus" = "lupus|\\bsle\\b",
  "Psoriasis" = "psoriasi",
  "Inflammatory Bowel Disease" = "inflammatory bowel|crohn|ulcerative colitis|\\bibd\\b",
  "Type 1 Diabetes" = "type 1 diabetes|\\bt1d\\b",
  "Cancer" = "cancer|carcinoma|tumou?r|neoplasm|leukemia|lymphoma",
  "Cardiovascular Disease" = "cardiovascular|coronary|myocardial|atherosclero|heart failure"
)

bc_ncbi_gene_summary <- function(entrez) {
  meta_base <- bc_meta("NCBI Gene", sprintf("db=gene&id=%s", entrez %||% NA), "eutils.ncbi.nlm.nih.gov/esummary.fcgi?db=gene", NA_character_)
  if (is.null(entrez) || is.na(entrez) || !nzchar(entrez)) return(list(ok = FALSE, reason = "No NCBI Gene ID available.", meta = meta_base))
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, reason = "httr2 is not installed.", meta = meta_base))
  res <- tryCatch({
    resp <- httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi") %>%
      httr2::req_url_query(db = "gene", id = entrez, retmode = "json", tool = "arthomix-explorer") %>%
      httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
    if (httr2::resp_status(resp) != 200) NULL else httr2::resp_body_json(resp)
  }, error = function(e) NULL)
  if (is.null(res)) { meta_base$status <- "Failed"; return(list(ok = FALSE, reason = "NCBI Gene lookup failed (network error or timeout).", meta = meta_base)) }
  rec <- res$result[[as.character(entrez)]]
  if (is.null(rec)) { meta_base$status <- "No results"; meta_base$n_records <- 0L; return(list(ok = FALSE, reason = sprintf("No NCBI Gene record found for Entrez ID %s.", entrez), meta = meta_base)) }
  meta_base$status <- "Success"; meta_base$n_records <- 1L
  list(ok = TRUE, reason = NULL,
       name = rec$name %||% NA_character_, description = rec$description %||% NA_character_, summary = rec$summary %||% NA_character_,
       aliases = if (!is.null(rec$otheraliases) && nzchar(rec$otheraliases)) rec$otheraliases else NA_character_,
       other_designations = rec$otherdesignations %||% NA_character_,
       chromosome = rec$chromosome %||% NA_character_, map_location = rec$maplocation %||% NA_character_, meta = meta_base)
}

bc_ensembl_gene_lookup <- function(symbol) {
  meta_base <- bc_meta("Ensembl (GRCh37)", symbol %||% NA, "grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}", NA_character_)
  if (is.null(symbol) || is.na(symbol) || !nzchar(symbol)) return(list(ok = FALSE, reason = "No gene symbol given.", meta = meta_base))
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, reason = "httr2 is not installed.", meta = meta_base))
  url <- sprintf("https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/%s", utils::URLencode(symbol, reserved = TRUE))
  res <- tryCatch({
    resp <- httr2::request(url) %>% httr2::req_url_query(`content-type` = "application/json") %>%
      httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
    if (httr2::resp_status(resp) != 200) NULL else httr2::resp_body_json(resp, simplifyVector = FALSE)
  }, error = function(e) NULL)
  if (is.null(res)) {
    meta_base$status <- "No results"; meta_base$n_records <- 0L
    return(list(ok = FALSE, reason = sprintf("No Ensembl (GRCh37) record found for gene symbol \"%s\".", symbol), meta = meta_base))
  }
  meta_base$status <- "Success"; meta_base$n_records <- 1L
  list(ok = TRUE, reason = NULL, ensembl_id = res$id %||% NA_character_, biotype = res$biotype %||% NA_character_,
       chr = res$seq_region_name %||% NA_character_, start = res$start %||% NA_real_, end = res[["end"]] %||% NA_real_,
       strand = res$strand %||% NA_integer_, description = sub("\\s*\\[Source.*$", "", res$description %||% NA_character_), meta = meta_base)
}

bc_ensembl_regulatory_overlap <- function(chr, start, end = NULL) {
  region_label <- if (!is.null(chr) && !is.na(chr) && !is.null(start) && !is.na(start)) sprintf("%s:%s-%s", chr, start, end %||% start) else NA_character_
  meta_base <- bc_meta("Ensembl Regulatory Build (GRCh37)", region_label, "grch37.rest.ensembl.org/overlap/region/human/{region}?feature=regulatory", NA_character_)
  if (is.null(chr) || is.na(chr) || is.null(start) || is.na(start)) return(list(ok = FALSE, reason = "No genomic coordinate available for this CpG.", meta = meta_base))
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, reason = "httr2 is not installed.", meta = meta_base))
  chr_n <- sub("^chr", "", chr, ignore.case = TRUE)
  win_start <- max(1, as.integer(start) - 1000); win_end <- as.integer(end %||% start) + 1000
  url <- sprintf("https://grch37.rest.ensembl.org/overlap/region/human/%s:%s-%s", chr_n, win_start, win_end)
  res <- tryCatch({
    resp <- httr2::request(url) %>% httr2::req_url_query(feature = "regulatory", `content-type` = "application/json") %>%
      httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
    if (httr2::resp_status(resp) != 200) NULL else httr2::resp_body_json(resp, simplifyVector = TRUE)
  }, error = function(e) NULL)
  if (is.null(res)) { meta_base$status <- "Failed"; return(list(ok = FALSE, reason = "Ensembl regulatory-overlap lookup failed (network error or timeout).", meta = meta_base)) }
  if (!is.data.frame(res) || nrow(res) == 0) { meta_base$status <- "No results"; meta_base$n_records <- 0L; return(list(ok = TRUE, reason = NULL, features = NULL, window = c(win_start, win_end), meta = meta_base)) }
  keep <- intersect(c("feature_type", "description", "start", "end", "id"), colnames(res))
  df <- res[, keep, drop = FALSE]
  df <- df[order(df$start), , drop = FALSE]
  meta_base$status <- "Success"; meta_base$n_records <- nrow(df)
  list(ok = TRUE, reason = NULL, features = df, window = c(win_start, win_end), meta = meta_base)
}

bc_wikipathways_pathways_for_gene <- function(entrez) {
  meta_base <- bc_meta("WikiPathways (msigdbr)", entrez %||% NA, "local msigdbr C2:CP:WIKIPATHWAYS cache (mp_get_wikipathways_termgene)", NA_character_)
  if (is.null(entrez) || is.na(entrez) || !nzchar(entrez)) return(list(ok = FALSE, pathways = NULL, reason = "No NCBI Gene ID available for WikiPathways lookup.", meta = meta_base))
  t2g <- tryCatch(mp_get_wikipathways_termgene(), error = function(e) NULL)
  if (is.null(t2g)) { meta_base$status <- "Failed"; return(list(ok = FALSE, pathways = NULL, reason = "WikiPathways gene sets (msigdbr) are not available in this deployment.", meta = meta_base)) }
  hit <- t2g$TERM2GENE[!is.na(t2g$TERM2GENE$ncbi_gene) & as.character(t2g$TERM2GENE$ncbi_gene) == as.character(entrez), , drop = FALSE]
  if (nrow(hit) == 0) { meta_base$status <- "Success"; meta_base$n_records <- 0L; return(list(ok = TRUE, pathways = data.frame(id = character(0), name = character(0)), reason = NULL, meta = meta_base)) }
  ids <- unique(hit$gs_exact_source)
  nm <- t2g$TERM2NAME$gs_name[match(ids, t2g$TERM2NAME$gs_exact_source)]
  meta_base$status <- "Success"; meta_base$n_records <- length(ids)
  list(ok = TRUE, pathways = data.frame(id = ids, name = nm, stringsAsFactors = FALSE), reason = NULL, meta = meta_base)
}

BC_OT_TRACTABILITY_TIERS <- c("Approved Drug", "Advanced Clinical", "Phase 1 Clinical")
BC_OT_MODALITY_LABELS <- c(SM = "Small molecule", AB = "Antibody / biologic", PR = "PROTAC-type degrader", OC = "Other clinical modality")

bc_opentargets_evidence_for_gene <- function(ensembl, top_n_diseases = 8) {
  meta_base <- bc_meta("Open Targets", ensembl %||% NA, "api.platform.opentargets.org/api/v4/graphql", NA_character_)
  if (is.null(ensembl) || is.na(ensembl) || !nzchar(ensembl)) return(list(ok = FALSE, reason = "No Ensembl Gene ID available for Open Targets lookup.", meta = meta_base))
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, reason = "httr2 is not installed.", meta = meta_base))
  query <- 'query($id:String!, $size:Int!){ target(ensemblId:$id){ id approvedSymbol
    tractability { label modality value }
    associatedDiseases(page:{index:0,size:$size}){ count rows { score disease { id name } datatypeScores { id score } } } } }'
  res <- tryCatch({
    resp <- httr2::request("https://api.platform.opentargets.org/api/v4/graphql") %>%
      httr2::req_headers(`Content-Type` = "application/json") %>%
      httr2::req_body_json(list(query = query, variables = list(id = ensembl, size = as.integer(top_n_diseases)))) %>%
      httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
    if (httr2::resp_status(resp) != 200) NULL else httr2::resp_body_json(resp, simplifyVector = FALSE)
  }, error = function(e) NULL)
  tgt <- res$data$target
  if (is.null(tgt)) { meta_base$status <- "No results"; meta_base$n_records <- 0L; return(list(ok = FALSE, reason = sprintf("No Open Targets entry found for Ensembl ID \"%s\".", ensembl), meta = meta_base)) }

  tract <- tgt$tractability
  modality_summary <- stats::setNames(lapply(names(BC_OT_MODALITY_LABELS), function(mod) {
    rows <- Filter(function(r) identical(r$modality, mod), tract)
    approved <- any(vapply(rows, function(r) identical(r$label, "Approved Drug") && isTRUE(r$value), logical(1)))
    tier <- Find(function(t) any(vapply(rows, function(r) identical(r$label, t) && isTRUE(r$value), logical(1))), BC_OT_TRACTABILITY_TIERS)
    feasible <- any(vapply(rows, function(r) isTRUE(r$value), logical(1)))
    if (approved) "Approved drug exists" else if (!is.null(tier)) tier else if (feasible) "Structurally feasible" else "No tractability evidence"
  }), names(BC_OT_MODALITY_LABELS))

  dis_rows <- tgt$associatedDiseases$rows %||% list()
  all_datatypes <- unique(unlist(lapply(dis_rows, function(r) vapply(r$datatypeScores, function(x) x$id, character(1)))))
  diseases <- if (length(dis_rows) > 0) {
    do.call(rbind, lapply(dis_rows, function(r) {
      dt <- stats::setNames(vapply(r$datatypeScores, function(x) x$score, numeric(1)), vapply(r$datatypeScores, function(x) x$id, character(1)))
      row <- data.frame(Disease = r$disease$name, `Overall score` = round(r$score, 3), check.names = FALSE, stringsAsFactors = FALSE)
      for (nm in all_datatypes) row[[nm]] <- if (nm %in% names(dt)) round(dt[[nm]], 3) else NA_real_
      row
    }))
  } else NULL
  if (!is.null(diseases)) colnames(diseases) <- gsub("\\bRna\\b", "RNA", gsub("_", " ", tools::toTitleCase(colnames(diseases))))

  meta_base$status <- "Success"; meta_base$n_records <- tgt$associatedDiseases$count %||% 0L
  list(ok = TRUE, n_diseases = tgt$associatedDiseases$count %||% 0, diseases = diseases, tractability = modality_summary, meta = meta_base)
}

bc_hpa_evidence_for_gene <- function(ensembl) {
  meta_base <- bc_meta("Human Protein Atlas", ensembl %||% NA, "www.proteinatlas.org/{ensembl}.json", NA_character_)
  if (is.null(ensembl) || is.na(ensembl) || !nzchar(ensembl)) return(list(ok = FALSE, reason = "No Ensembl Gene ID available for Human Protein Atlas lookup.", meta = meta_base))
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, reason = "httr2 is not installed.", meta = meta_base))
  res <- tryCatch({
    resp <- httr2::request(sprintf("https://www.proteinatlas.org/%s.json", ensembl)) %>%
      httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
    if (httr2::resp_status(resp) != 200) NULL else httr2::resp_body_json(resp, simplifyVector = FALSE)
  }, error = function(e) NULL)
  if (is.null(res)) { meta_base$status <- "No results"; meta_base$n_records <- 0L; return(list(ok = FALSE, reason = sprintf("No Human Protein Atlas entry found for Ensembl ID \"%s\".", ensembl), meta = meta_base)) }

  named_list_to_df <- function(x, value_label) {
    if (is.null(x) || length(x) == 0) return(NULL)
    out <- data.frame(Tissue = names(x), Value = vapply(x, function(v) as.character(v %||% NA), character(1)), stringsAsFactors = FALSE)
    stats::setNames(out, c("Tissue / cell type", value_label))
  }
  meta_base$status <- "Success"; meta_base$n_records <- 1L
  list(ok = TRUE,
       tissue_specificity = res[["RNA tissue specificity"]] %||% "Not available",
       tissue_top = named_list_to_df(res[["RNA tissue specific nTPM"]], "nTPM"),
       blood_specificity = res[["RNA blood lineage specificity"]] %||% "Not available",
       blood_top = named_list_to_df(res[["RNA blood lineage specific nTPM"]], "nTPM"),
       blood_cluster = res[["Blood expression cluster"]] %||% NA_character_,
       secretome = res[["Secretome location"]] %||% NA_character_,
       protein_class = if (!is.null(res[["Protein class"]])) paste(unlist(res[["Protein class"]]), collapse = ", ") else NA_character_,
       meta = meta_base)
}

bc_gwas_catalog_by_gene <- function(symbol, max_results = 15) {
  meta_base <- bc_meta("GWAS Catalog (NHGRI-EBI)", symbol %||% NA, "www.ebi.ac.uk/gwas/api/search?q={symbol}", NA_character_)
  if (is.null(symbol) || is.na(symbol) || !nzchar(symbol)) return(list(ok = FALSE, traits = NULL, reason = "No gene symbol given.", meta = meta_base))
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, traits = NULL, reason = "httr2 is not installed.", meta = meta_base))
  res <- tryCatch({
    resp <- httr2::request("https://www.ebi.ac.uk/gwas/api/search") %>%
      httr2::req_url_query(q = symbol, facet = "false", rows = 200) %>%
      httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
    if (httr2::resp_status(resp) != 200) NULL else httr2::resp_body_json(resp, check_type = FALSE, simplifyVector = FALSE)
  }, error = function(e) NULL)
  if (is.null(res)) { meta_base$status <- "Failed"; return(list(ok = FALSE, traits = NULL, reason = "GWAS Catalog lookup failed (network error or timeout).", meta = meta_base)) }
  docs <- Filter(function(d) identical(d$resourcename, "trait"), res$response$docs %||% list())
  if (length(docs) == 0) { meta_base$status <- "No results"; meta_base$n_records <- 0L; return(list(ok = TRUE, traits = NULL, reason = NULL, meta = meta_base)) }
  df <- do.call(rbind, lapply(utils::head(docs, max_results), function(d) {
    data.frame(Trait = d$mappedTrait %||% NA_character_, `EFO ID` = d$shortForm[[1]] %||% NA_character_,
               Studies = d$studyCount %||% NA_integer_, Associations = d$associationCount %||% NA_integer_,
               `Reported as` = if (!is.null(d$reportedTrait)) paste(utils::head(unlist(d$reportedTrait), 2), collapse = "; ") else NA_character_,
               check.names = FALSE, stringsAsFactors = FALSE)
  }))
  df <- df[order(-df$Associations), , drop = FALSE]
  meta_base$status <- "Success"; meta_base$n_records <- length(docs)
  list(ok = TRUE, traits = df, reason = NULL, n_total = length(docs), meta = meta_base)
}

bc_encode_search <- function(query, max_results = 15) {
  meta_base <- bc_meta("ENCODE", query %||% NA, "www.encodeproject.org/search/?type=Experiment", NA_character_)
  if (is.null(query) || !nzchar(trimws(query %||% ""))) return(list(ok = FALSE, experiments = NULL, reason = "No search term given.", meta = meta_base))
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, experiments = NULL, reason = "httr2 is not installed.", meta = meta_base))
  res <- tryCatch({
    resp <- httr2::request("https://www.encodeproject.org/search/") %>%
      httr2::req_url_query(searchTerm = query, type = "Experiment", status = "released", format = "json", limit = max_results) %>%
      httr2::req_headers(Accept = "application/json") %>%
      httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
    if (httr2::resp_status(resp) != 200) NULL else httr2::resp_body_json(resp, simplifyVector = FALSE)
  }, error = function(e) NULL)
  if (is.null(res)) { meta_base$status <- "Failed"; return(list(ok = FALSE, experiments = NULL, reason = "ENCODE lookup failed (network error or timeout).", meta = meta_base)) }
  rows <- res[["@graph"]] %||% list()
  if (length(rows) == 0) { meta_base$status <- "No results"; meta_base$n_records <- 0L; return(list(ok = TRUE, experiments = NULL, reason = NULL, meta = meta_base)) }
  df <- do.call(rbind, lapply(rows, function(r) {
    data.frame(Accession = r$accession %||% NA_character_, Assay = r$assay_title %||% NA_character_,
               Biosample = r$biosample_summary %||% NA_character_, Status = r$status %||% NA_character_,
               `Link path` = r[["@id"]] %||% NA_character_, check.names = FALSE, stringsAsFactors = FALSE)
  }))
  meta_base$status <- "Success"; meta_base$n_records <- res$total %||% nrow(df)
  list(ok = TRUE, experiments = df, reason = NULL, n_total = res$total %||% nrow(df), meta = meta_base)
}

bc_geo_search <- function(query, max_results = 10) {
  meta_base <- bc_meta("GEO (NCBI E-utilities, db=gds)", query %||% NA, "eutils.ncbi.nlm.nih.gov/esearch.fcgi?db=gds", NA_character_)
  if (is.null(query) || !nzchar(trimws(query %||% ""))) return(list(ok = FALSE, series = NULL, reason = "No search query given.", meta = meta_base))
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, series = NULL, reason = "httr2 is not installed.", meta = meta_base))
  max_results <- min(max(as.integer(max_results %||% 10), 1L), 25L)
  res <- tryCatch({
    esearch <- httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi") %>%
      httr2::req_url_query(db = "gds", term = query, retmode = "json", retmax = max_results, tool = "arthomix-explorer") %>%
      httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
    if (httr2::resp_status(esearch) != 200) return(NULL)
    ids <- unlist(httr2::resp_body_json(esearch)$esearchresult$idlist)
    if (length(ids) == 0) return(list(ids = character(0), recs = NULL))
    esummary <- httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi") %>%
      httr2::req_url_query(db = "gds", id = paste(ids, collapse = ","), retmode = "json", tool = "arthomix-explorer") %>%
      httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
    if (httr2::resp_status(esummary) != 200) return(NULL)
    list(ids = ids, recs = httr2::resp_body_json(esummary)$result)
  }, error = function(e) NULL)
  if (is.null(res)) { meta_base$status <- "Failed"; return(list(ok = FALSE, series = NULL, reason = "GEO lookup failed (network error, timeout, or malformed response).", meta = meta_base)) }
  if (length(res$ids) == 0) { meta_base$status <- "No results"; meta_base$n_records <- 0L; return(list(ok = TRUE, series = NULL, reason = NULL, meta = meta_base)) }
  rows <- lapply(res$ids, function(id) {
    r <- res$recs[[id]]
    if (is.null(r)) return(NULL)
    data.frame(Accession = r$accession %||% NA_character_, Title = r$title %||% NA_character_,
               `Data type` = r$gdstype %||% NA_character_, Organism = r$taxon %||% NA_character_,
               `Samples` = r$n_samples %||% NA_integer_, Date = r$pdat %||% NA_character_,
               check.names = FALSE, stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, Filter(Negate(is.null), rows))
  meta_base$status <- "Success"; meta_base$n_records <- if (!is.null(df)) nrow(df) else 0L
  list(ok = TRUE, series = df, reason = NULL, meta = meta_base)
}

bc_biostudies_search <- function(query, max_results = 10) {
  meta_base <- bc_meta("BioStudies / ArrayExpress (EBI)", query %||% NA, "www.ebi.ac.uk/biostudies/api/v1/search", NA_character_)
  if (is.null(query) || !nzchar(trimws(query %||% ""))) return(list(ok = FALSE, studies = NULL, reason = "No search query given.", meta = meta_base))
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, studies = NULL, reason = "httr2 is not installed.", meta = meta_base))
  res <- tryCatch({
    resp <- httr2::request("https://www.ebi.ac.uk/biostudies/api/v1/search") %>%
      httr2::req_url_query(query = query, pageSize = max_results) %>%
      httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
    if (httr2::resp_status(resp) != 200) NULL else httr2::resp_body_json(resp, simplifyVector = FALSE)
  }, error = function(e) NULL)
  if (is.null(res)) { meta_base$status <- "Failed"; return(list(ok = FALSE, studies = NULL, reason = "BioStudies lookup failed (network error or timeout).", meta = meta_base)) }
  hits <- res$hits %||% list()
  if (length(hits) == 0) { meta_base$status <- "No results"; meta_base$n_records <- 0L; return(list(ok = TRUE, studies = NULL, reason = NULL, meta = meta_base)) }
  df <- do.call(rbind, lapply(hits, function(h) {
    data.frame(Accession = h$accession %||% NA_character_, Title = h$title %||% NA_character_,
               Type = h$type %||% NA_character_, `Release date` = h$release_date %||% NA_character_,
               check.names = FALSE, stringsAsFactors = FALSE)
  }))
  meta_base$status <- "Success"; meta_base$n_records <- res$totalHits %||% nrow(df)
  list(ok = TRUE, studies = df, reason = NULL, n_total = res$totalHits %||% nrow(df), meta = meta_base)
}

bc_literature_search <- function(query, max_results = 12) {
  meta_base <- bc_meta("PubMed (NCBI E-utilities)", query %||% NA, "eutils.ncbi.nlm.nih.gov/esearch.fcgi?db=pubmed", NA_character_)
  if (is.null(query) || !nzchar(trimws(query %||% ""))) return(list(ok = FALSE, papers = NULL, reason = "No search query given.", meta = meta_base))
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, papers = NULL, reason = "httr2 is not installed.", meta = meta_base))
  max_results <- min(max(as.integer(max_results %||% 12), 1L), 30L)
  res <- tryCatch({
    esearch <- httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi") %>%
      httr2::req_url_query(db = "pubmed", term = query, retmode = "json", retmax = max_results, sort = "relevance", tool = "arthomix-biomarkercard") %>%
      httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
    if (httr2::resp_status(esearch) != 200) return(NULL)
    ids <- unlist(httr2::resp_body_json(esearch)$esearchresult$idlist)
    if (length(ids) == 0) return(data.frame())
    esummary <- httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi") %>%
      httr2::req_url_query(db = "pubmed", id = paste(ids, collapse = ","), retmode = "json", tool = "arthomix-biomarkercard") %>%
      httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
    if (httr2::resp_status(esummary) != 200) return(NULL)
    recs <- httr2::resp_body_json(esummary)$result
    rows <- lapply(ids, function(pmid) {
      r <- recs[[pmid]]
      if (is.null(r)) return(NULL)
      authors <- vapply(r$authors, function(a) a$name %||% "?", character(1))
      author_str <- if (length(authors) > 1) paste0(authors[1], " et al.") else if (length(authors) == 1) authors[1] else "Unknown author"
      data.frame(Title = r$title %||% "Untitled", Authors = author_str, Journal = r$fulljournalname %||% r$source %||% NA_character_,
                 Year = substr(r$pubdate %||% "", 1, 4), PMID = pmid, Type = if (length(r$pubtype) > 0) paste(unlist(r$pubtype), collapse = ", ") else NA_character_,
                 stringsAsFactors = FALSE)
    })
    do.call(rbind, Filter(Negate(is.null), rows))
  }, error = function(e) NULL)
  if (is.null(res)) { meta_base$status <- "Failed"; return(list(ok = FALSE, papers = NULL, reason = "PubMed lookup failed (network error, timeout, or malformed response).", meta = meta_base)) }
  meta_base$status <- "Success"; meta_base$n_records <- if (is.data.frame(res)) nrow(res) else 0L
  list(ok = TRUE, papers = if (is.data.frame(res) && nrow(res) > 0) res else NULL, reason = NULL, meta = meta_base)
}

BC_LITERATURE_PRESETS <- c(
  "Gene only" = "%s", "Gene + methylation" = "%s methylation", "Gene + CpG" = "%s CpG methylation",
  "Gene + epigenetic biomarker" = "%s epigenetic biomarker", "Gene + rheumatoid arthritis" = "%s rheumatoid arthritis"
)
bc_literature_query <- function(identifier, preset, custom_disease = NULL) {
  if (identical(preset, "disease")) {
    ctx <- if (!is.null(custom_disease) && nzchar(trimws(custom_disease))) trimws(custom_disease) else "disease"
    return(sprintf("%s %s", identifier, ctx))
  }
  tmpl <- BC_LITERATURE_PRESETS[[preset]] %||% "%s"
  sprintf(tmpl, identifier)
}

BC_LIT_CLASS_PATTERNS <- list(
  "EWAS" = "epigenome-wide|\\bewas\\b",
  "DNA methylation" = "methylat",
  "Biomarker study" = "biomarker",
  "Disease association" = "associat|risk of|susceptibility",
  "Functional study" = "function|regulat|expression|knockdown|knockout",
  "Mechanistic study" = "mechanism|pathway|signal(l)?ing",
  "Review" = "^review|systematic review|meta-analysis",
  "Clinical study" = "clinical trial|cohort|patients? with",
  "Validation study" = "validat|replicat"
)
bc_literature_classify <- function(title, abstract = NULL) {
  text <- tolower(paste(title %||% "", abstract %||% ""))
  if (!nzchar(trimws(text))) return(character(0))
  hits <- names(BC_LIT_CLASS_PATTERNS)[vapply(BC_LIT_CLASS_PATTERNS, function(p) grepl(p, text, ignore.case = TRUE), logical(1))]
  hits
}

BC_ID_TYPE_PATTERNS <- list(cpg = "^cg[0-9]{6,}$", ensembl = "^ENSG[0-9]{6,}(\\.[0-9]+)?$", entrez = "^[0-9]+$")
bc_detect_identifier_type <- function(x) {
  x <- trimws(x)
  if (grepl(BC_ID_TYPE_PATTERNS$cpg, x, ignore.case = TRUE)) return("cpg")
  if (grepl(BC_ID_TYPE_PATTERNS$ensembl, x, ignore.case = TRUE)) return("ensembl")
  if (grepl(BC_ID_TYPE_PATTERNS$entrez, x)) return("entrez")
  "gene_symbol"
}

bc_split_tokens <- function(text) {
  toks <- trimws(unlist(strsplit(text %||% "", "[,\n\t ]+")))
  unique(toks[nzchar(toks)])
}

bc_resolve_identifiers <- function(raw_ids, array_type = "450K") {
  raw_ids <- unique(trimws(as.character(raw_ids)))
  raw_ids <- raw_ids[nzchar(raw_ids)]
  if (length(raw_ids) == 0) return(list(ok = FALSE, reason = "No identifiers were provided.", df = NULL, n_submitted = 0L, n_resolved = 0L, n_unresolved = 0L))
  types <- vapply(raw_ids, bc_detect_identifier_type, character(1))
  cpg_ids <- raw_ids[types == "cpg"]
  gene_ids <- raw_ids[types != "cpg"]

  cpg_rows <- if (length(cpg_ids) > 0) do.call(rbind, lapply(cpg_ids, function(id) {
    r <- bc_resolve_cpg(id, array_type)
    resolved <- isTRUE(r$found_in_manifest) || isTRUE(r$found_in_champ)
    gene <- if (!is.na(r$champ_gene) && nzchar(r$champ_gene)) r$champ_gene
            else if (!is.na(r$gene_names) && nzchar(r$gene_names)) trimws(strsplit(r$gene_names, ";")[[1]][1]) else NA_character_
    data.frame(input_id = id, detected_type = "CpG probe", resolved = resolved,
               resolved_id = if (resolved) id else NA_character_,
               status_label = if (resolved) "Resolved (found in manifest annotation)" else sprintf("Unresolved - not found in the %s manifest", array_type),
               annotated_gene = gene, stringsAsFactors = FALSE)
  })) else NULL

  gene_rows <- NULL
  if (length(gene_ids) > 0) {
    harm <- tryCatch(cx_harmonize_gene_ids(gene_ids), error = function(e) list(ok = FALSE))
    if (isTRUE(harm$ok)) {
      hdf <- harm$df
      status_labels <- c(exact_symbol = "Resolved (exact gene symbol)", exact_entrez = "Resolved (NCBI Entrez ID)",
                          exact_ensembl = "Resolved (Ensembl Gene ID)", alias_resolved = "Resolved (via known alias)",
                          ambiguous = "Ambiguous - multiple candidate symbols, not guessed", unmatched = "Unresolved - identifier not recognized")
      gene_rows <- data.frame(input_id = hdf$input_id, detected_type = "Gene identifier",
                               resolved = hdf$match_type %in% c("exact_symbol", "exact_entrez", "exact_ensembl", "alias_resolved"),
                               resolved_id = hdf$canonical_symbol, status_label = unname(status_labels[hdf$match_type]),
                               annotated_gene = hdf$canonical_symbol, stringsAsFactors = FALSE)
    } else {
      gene_rows <- data.frame(input_id = gene_ids, detected_type = "Gene identifier", resolved = FALSE,
                               resolved_id = NA_character_, status_label = "Unresolved - gene identifier harmonization is unavailable in this deployment",
                               annotated_gene = NA_character_, stringsAsFactors = FALSE)
    }
  }
  df <- rbind(cpg_rows, gene_rows)
  df <- df[match(raw_ids, df$input_id), , drop = FALSE]
  list(ok = TRUE, reason = NULL, df = df, n_submitted = length(raw_ids), n_resolved = sum(df$resolved), n_unresolved = sum(!df$resolved))
}

bc_aggregate_convergence <- function(long_df, item_col, id_col = "ID") {
  if (is.null(long_df) || nrow(long_df) == 0) return(NULL)
  agg <- stats::aggregate(stats::reformulate(item_col, response = id_col), long_df, function(x) paste(sort(unique(x)), collapse = ", "))
  agg$`Biomarker count` <- vapply(strsplit(agg[[id_col]], ", "), length, integer(1))
  agg <- agg[order(-agg$`Biomarker count`), , drop = FALSE]
  colnames(agg)[colnames(agg) == id_col] <- "Biomarkers"
  agg
}

bc_panel_gene_rows <- function(gene_symbols) {
  if (length(gene_symbols) == 0) return(NULL)
  do.call(rbind, lapply(gene_symbols, function(g) {
    gs <- bc_gene_structure(g)
    data.frame(Gene = g, `NCBI Entrez ID` = if (isTRUE(gs$ok)) gs$entrez else NA_character_,
               `Ensembl Gene ID` = if (isTRUE(gs$ok)) gs$ensembl else NA_character_,
               Chromosome = if (isTRUE(gs$ok)) gs$chr else NA_character_,
               check.names = FALSE, stringsAsFactors = FALSE)
  }))
}

bc_panel_cpg_rows <- function(cpg_ids, array_type) {
  if (length(cpg_ids) == 0) return(NULL)
  do.call(rbind, lapply(cpg_ids, function(cpg) {
    r <- bc_resolve_cpg(cpg, array_type)
    dmp_f <- bc_dmp_lookup(cpg, "female"); dmp_m <- bc_dmp_lookup(cpg, "male")
    gene <- if (!is.na(r$champ_gene) && nzchar(r$champ_gene)) r$champ_gene
            else if (!is.na(r$gene_names) && nzchar(r$gene_names)) trimws(strsplit(r$gene_names, ";")[[1]][1]) else NA_character_
    data.frame(CpG = cpg, Chromosome = r$chr, Position = r$pos, Gene = gene,
               `Island relation` = bc_island_context_label(r$island_relation),
               `Female FDR (pipeline)` = if (!is.null(dmp_f)) signif(dmp_f$fdr_bacon, 3) else NA_real_,
               `Male FDR (pipeline)` = if (!is.null(dmp_m)) signif(dmp_m$fdr_bacon, 3) else NA_real_,
               check.names = FALSE, stringsAsFactors = FALSE)
  }))
}

bc_panel_disease_convergence_genes <- function(gene_symbols, max_genes = 15) {
  genes <- utils::head(unique(gene_symbols), max_genes)
  rows <- lapply(genes, function(g) {
    gs <- bc_gene_structure(g)
    if (!isTRUE(gs$ok) || is.na(gs$ensembl)) return(NULL)
    r <- bc_opentargets_evidence_for_gene(gs$ensembl, top_n_diseases = 10)
    if (!isTRUE(r$ok) || is.null(r$diseases) || nrow(r$diseases) == 0) return(NULL)
    data.frame(ID = g, Disease = r$diseases$Disease, stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(list(ok = FALSE, reason = "No Open Targets disease associations were found for any resolved gene in this list.", table = NULL, n_queried = length(genes)))
  list(ok = TRUE, reason = NULL, table = bc_aggregate_convergence(do.call(rbind, rows), "Disease"), n_queried = length(genes))
}

bc_panel_ewas_trait_convergence_cpgs <- function(cpg_ids, max_cpgs = 15) {
  cpgs <- utils::head(unique(cpg_ids), max_cpgs)
  rows <- lapply(cpgs, function(cpg) {
    r <- bc_ewascatalog_query("cpg", cpg)
    if (!isTRUE(r$ok) || nrow(r$df) == 0) return(NULL)
    data.frame(ID = cpg, Trait = r$df$trait, stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(list(ok = FALSE, reason = "No EWAS Catalog trait associations were found for any CpG in this list.", table = NULL, n_queried = length(cpgs)))
  list(ok = TRUE, reason = NULL, table = bc_aggregate_convergence(do.call(rbind, rows), "Trait"), n_queried = length(cpgs))
}

bc_panel_pathway_convergence_genes <- function(gene_symbols, max_genes = 15) {
  genes <- utils::head(unique(gene_symbols), max_genes)
  rows <- lapply(genes, function(g) {
    gs <- bc_gene_structure(g)
    if (!isTRUE(gs$ok)) return(NULL)
    k <- bc_kegg_pathways_for_gene(gs$entrez)
    r <- bc_reactome_pathways_for_gene(g)
    out <- NULL
    if (isTRUE(k$ok) && !is.null(k$pathways) && nrow(k$pathways) > 0) out <- rbind(out, data.frame(ID = g, Pathway = k$pathways$name, stringsAsFactors = FALSE))
    if (isTRUE(r$ok) && !is.null(r$pathways) && nrow(r$pathways) > 0) out <- rbind(out, data.frame(ID = g, Pathway = r$pathways$displayName, stringsAsFactors = FALSE))
    out
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(list(ok = FALSE, reason = "No KEGG/Reactome pathways were found for any resolved gene in this list.", table = NULL, n_queried = length(genes)))
  list(ok = TRUE, reason = NULL, table = bc_aggregate_convergence(do.call(rbind, rows), "Pathway"), n_queried = length(genes))
}

bc_plot_region <- function(chr, pos, island, gene_struct, flank_bp) {
  validate(need(!is.null(chr) && !is.na(chr) && !is.null(pos) && !is.na(pos), "No genomic coordinate is available to plot for this CpG."))
  xmin <- pos - flank_bp; xmax <- pos + flank_bp
  p <- ggplot() +
    coord_cartesian(xlim = c(xmin, xmax), ylim = c(0, 1)) +
    theme_arthomix() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    labs(x = sprintf("%s position (bp) - hg19", chr), y = NULL)

  if (!is.null(island) && identical(island$chr, chr)) {
    isl_df <- data.frame(xmin = island$start, xmax = island$end, ymin = 0.55, ymax = 0.78)
    p <- p + geom_rect(data = isl_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                        fill = ARTHOMIX_COLORS$aqua, alpha = 0.55, inherit.aes = FALSE) +
      annotate("text", x = (island$start + island$end) / 2, y = 0.88, label = "CpG Island", size = 3.1, color = ARTHOMIX_COLORS$aqua)
  }

  if (isTRUE(gene_struct$ok) && !is.na(gene_struct$chr) && identical(gene_struct$chr, chr)) {
    gy <- 0.22
    gene_df <- data.frame(x = gene_struct$start, xend = gene_struct$end, y = gy, yend = gy)
    p <- p + geom_segment(data = gene_df, aes(x = x, xend = xend, y = y, yend = yend),
                           color = ARTHOMIX_COLORS$ink_muted, linewidth = 1, inherit.aes = FALSE)
    if (!is.null(gene_struct$exons) && length(gene_struct$exons) > 0) {
      ex_df <- as.data.frame(gene_struct$exons)
      ex_df <- ex_df[ex_df$end >= xmin & ex_df$start <= xmax, , drop = FALSE]
      if (nrow(ex_df) > 0) {
        ex_df$ymin <- gy - 0.05; ex_df$ymax <- gy + 0.05
        p <- p + geom_rect(data = ex_df, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax),
                            fill = ARTHOMIX_COLORS$blue, inherit.aes = FALSE)
      }
    }
    p <- p + annotate("text", x = xmin, y = gy + 0.14, label = gene_struct$symbol, hjust = 0, size = 3.1,
                       fontface = "italic", color = ARTHOMIX_COLORS$ink_muted)
  }

  cpg_df <- data.frame(x = pos, y = 0.42)
  p <- p + geom_point(data = cpg_df, aes(x = x, y = y), color = ARTHOMIX_STATUS$critical, size = 3.2, inherit.aes = FALSE) +
    annotate("text", x = pos, y = 0.32, label = "CpG", size = 3, color = ARTHOMIX_STATUS$critical)
  p
}

bc_plot_ideogram <- function(chr, pos, cytoband_df) {
  validate(need(!is.null(cytoband_df), "The bundled UCSC cytoband reference (hg19) is not available in this deployment."))
  validate(need(!is.null(chr) && !is.na(chr), "No chromosome is available to plot."))
  d <- cytoband_df[cytoband_df$chrom == chr, , drop = FALSE]
  validate(need(nrow(d) > 0, sprintf("No cytoband data available for %s.", chr)))
  d$fillcol <- BC_GIESTAIN_COLORS[d$gieStain]
  d$fillcol[is.na(d$fillcol)] <- "#DDDDDD"
  mk_df <- data.frame(x = pos, y = 0)
  ggplot() +
    geom_rect(data = d, aes(xmin = chromStart, xmax = chromEnd, ymin = -0.35, ymax = 0.35),
              fill = d$fillcol, color = "grey45", linewidth = 0.15) +
    geom_point(data = mk_df, aes(x = x, y = y), color = ARTHOMIX_STATUS$critical, size = 4, shape = 24, fill = ARTHOMIX_STATUS$critical) +
    annotate("text", x = pos, y = 0.62, label = sprintf("%s:%s", chr, format(round(pos), big.mark = ",")),
             size = 3.1, color = ARTHOMIX_STATUS$critical, fontface = "bold") +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(x = sprintf("%s position (bp) - hg19 cytoBandIdeo", chr), y = NULL, title = sprintf("Chromosome %s", sub("^chr", "", chr))) +
    theme_arthomix() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())
}

bc_plot_methylation_dist <- function(df, facet_sex = FALSE) {
  if (isTRUE(facet_sex) && "sex" %in% names(df)) df <- df[!is.na(df$sex), , drop = FALSE]
  df <- df[!is.na(df$beta) & !is.na(df$group), , drop = FALSE]
  validate(need(nrow(df) > 0, "No non-missing beta values are available for this CpG in the current groups."))
  p <- ggplot(df, aes(x = group, y = beta, fill = group)) +
    geom_violin(alpha = 0.35, color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.08, alpha = 0.5, size = 1.1) +
    scale_fill_manual(values = arthomix_pair(sort(unique(df$group)))) +
    labs(x = NULL, y = "Beta value (0-1)") + theme_arthomix() + theme(legend.position = "none")
  if (isTRUE(facet_sex) && "sex" %in% names(df)) p <- p + facet_wrap(~sex)
  p
}

bc_ggsave_datauri <- function(plot, width = 8, height = 3.2, dpi = 110) {
  if (is.null(plot) || !requireNamespace("base64enc", quietly = TRUE)) return(NULL)
  tf <- tempfile(fileext = ".png")
  ok <- tryCatch({ ggplot2::ggsave(tf, plot = plot, width = width, height = height, dpi = dpi, bg = "white"); TRUE }, error = function(e) FALSE)
  if (!ok || !file.exists(tf)) return(NULL)
  uri <- base64enc::dataURI(file = tf, mime = "image/png")
  unlink(tf)
  uri
}

bc_fmt_field <- function(x) {
  if (is.null(x) || length(x) == 0) return("Not available")
  if (length(x) == 1 && is.na(x)) return("Not available")
  if (is.character(x) && !nzchar(trimws(x))) return("Not available")
  as.character(x)
}
bc_kv_table <- function(pairs) {
  df <- data.frame(Field = names(pairs), Value = vapply(pairs, bc_fmt_field, character(1)), stringsAsFactors = FALSE)
  DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE, ordering = FALSE), class = "stripe hover compact")
}

bc_section_summary <- function(d) {
  r <- d$resolved
  cand <- list(Female = d$dmp_female, Male = d$dmp_male)
  cand <- cand[!vapply(cand, is.null, logical(1))]
  headline <- NULL; headline_label <- NULL
  if (length(cand) > 0) {
    fdrs <- vapply(cand, function(x) if (is.na(x$fdr_bacon)) Inf else x$fdr_bacon, numeric(1))
    headline_label <- names(cand)[which.min(fdrs)]
    headline <- cand[[headline_label]]
  }
  dbeta <- if (!is.null(headline)) headline$dbeta else NA_real_
  fdr <- if (!is.null(headline)) headline$fdr_bacon else NA_real_
  direction <- if (!is.null(headline) && !is.na(dbeta)) {
    if (dbeta > 0) "Hypermethylated" else if (dbeta < 0) "Hypomethylated" else "No change"
  } else NA_character_
  sig <- if (!is.na(fdr)) (fdr <= 0.05) else NA
  gene_relationship <- if (!is.na(r$gene_names) && nzchar(r$gene_names)) paste(unique(trimws(strsplit(r$gene_names, ";")[[1]])), collapse = ", ") else NA_character_

  pairs <- list(
    "Gene" = d$primary_gene,
    "Chromosome" = r$chr,
    "Genomic position" = if (!is.na(r$pos)) sprintf("%s:%s", r$chr, format(r$pos, big.mark = ",")) else NA_character_,
    "Genome build" = "GRCh37 / hg19",
    "Methylation direction" = direction,
    "Delta Beta" = if (!is.na(dbeta)) sprintf("%.4f (%s stratum, SVA/bacon pipeline)", dbeta, headline_label) else NA_character_,
    "P-value" = if (!is.null(headline) && !is.na(headline$p_bacon)) signif(headline$p_bacon, 4) else NA_real_,
    "Adjusted P-value (FDR)" = if (!is.na(fdr)) signif(fdr, 4) else NA_real_,
    "Effect size (|Delta Beta|)" = if (!is.na(dbeta)) round(abs(dbeta), 4) else NA_real_,
    "CpG context" = bc_island_context_label(r$island_relation),
    "Genomic feature" = bc_feature_label(r$gene_group),
    "Gene relationship" = gene_relationship,
    "Biomarker status" = if (is.na(sig)) NA_character_ else if (sig) "Significant (FDR <= 0.05)" else "Not significant"
  )
  div(class = "card",
      div(class = "card-title", icon("id-card"), sprintf("Biomarker: %s", d$cpg)),
      p(class = "submodule-desc", "Source: SVA/bacon-corrected sex-stratified DMP pipeline (preloaded) + Illumina manifest / ChAMPdata annotation (hg19)."),
      bc_kv_table(pairs)
  )
}

bc_section_genomic_context <- function(d) {
  r <- d$resolved; isl <- d$island
  dist_txt <- if (!is.null(isl) && !is.na(r$pos)) {
    if (r$pos >= isl$start && r$pos <= isl$end) "0 (CpG lies within the island)"
    else if (r$pos < isl$start) sprintf("%s bp upstream of island start", format(isl$start - r$pos, big.mark = ","))
    else sprintf("%s bp downstream of island end", format(r$pos - isl$end, big.mark = ","))
  } else NA_character_

  pairs <- list(
    "CpG ID" = d$cpg, "Chromosome" = r$chr,
    "Exact genomic coordinate" = if (!is.na(r$pos)) format(r$pos, big.mark = ",") else NA_character_,
    "Genome assembly / build" = "GRCh37 / hg19", "Strand" = r$strand,
    "CpG island name / coordinates" = if (!is.null(isl)) sprintf("%s:%s-%s", isl$chr, format(isl$start, big.mark = ","), format(isl$end, big.mark = ",")) else NA_character_,
    "CpG island length (bp)" = if (!is.null(isl)) format(isl$length, big.mark = ",") else NA_character_,
    "Distance from CpG to island boundary" = dist_txt,
    "CpG island relation" = bc_island_context_label(r$island_relation),
    "Promoter / TSS relation" = bc_feature_label(r$gene_group),
    "Gene body relationship" = bc_feature_label(r$gene_group),
    "Enhancer relationship" = bc_enhancer_flag(r, d$array_type),
    "Regulatory feature" = if (!is.na(r$regulatory_group) && nzchar(r$regulatory_group))
      sprintf("%s%s", r$regulatory_group, if (!is.na(r$regulatory_name) && nzchar(r$regulatory_name)) sprintf(" (%s)", r$regulatory_name) else "")
      else "No annotated regulatory feature",
    "DNase hypersensitivity (DHS)" = bc_dhs_flag(r, d$array_type),
    "Nearby / overlapping gene(s)" = if (!is.na(r$gene_names) && nzchar(r$gene_names)) paste(unique(trimws(strsplit(r$gene_names, ";")[[1]])), collapse = ", ") else "None annotated"
  )
  div(class = "card",
      div(class = "card-title", icon("map-location-dot"), "CpG Genomic Context"),
      p(class = "submodule-desc", sprintf("Source: %s manifest annotation (Illumina/UCSC, hg19).", d$array_type)),
      bc_kv_table(pairs)
  )
}

bc_primary_gene_detail <- function(d) {
  gs <- d$gene_struct
  if (!isTRUE(gs$ok)) return(div(class = "empty-note", icon("circle-info"), gs$reason %||% "Gene structure could not be resolved."))
  tss <- if (identical(gs$strand, "-")) gs$end else gs$start
  pos <- d$resolved$pos
  dist <- if (!is.na(pos) && !is.na(gs$start) && !is.na(gs$end)) {
    if (pos >= gs$start && pos <= gs$end) "0 (CpG falls within the gene body)" else format(abs(pos - tss), big.mark = ",")
  } else NA_character_
  pairs <- list(
    "Gene symbol" = gs$symbol, "Gene name" = gs$genename, "NCBI Gene ID (Entrez)" = gs$entrez, "Ensembl Gene ID" = gs$ensembl,
    "Chromosome" = gs$chr,
    "Genomic location" = if (!is.na(gs$start) && !is.na(gs$end)) sprintf("%s:%s-%s", gs$chr, format(gs$start, big.mark = ","), format(gs$end, big.mark = ",")) else NA_character_,
    "Strand / transcription direction" = if (identical(gs$strand, "+")) "+ (forward)" else if (identical(gs$strand, "-")) "- (reverse)" else gs$strand,
    "Distance from CpG to TSS (bp)" = dist,
    "Top GO Biological Process terms" = if (!is.null(d$go_terms) && nrow(d$go_terms) > 0) paste(d$go_terms$TERM, collapse = "; ") else "Not available"
  )
  tagList(
    bc_kv_table(pairs),
    p(style = "font-size:11px; color:var(--color-ink-muted); margin-top:6px;",
      "Source: org.Hs.eg.db (NCBI/Ensembl mapping, GO annotation), TxDb.Hsapiens.UCSC.hg19.knownGene (gene structure), GO.db (GO term names).")
  )
}

bc_section_genes <- function(d) {
  rel_tab <- bc_gene_relationship_table(d$resolved)
  div(class = "card",
      div(class = "card-title", icon("dna"), "Associated Gene(s)"),
      p(class = "submodule-desc", "Source: Illumina manifest UCSC RefGene annotation (hg19) for the relationship table; NCBI/Ensembl/GO detail below is for the primary associated gene only."),
      if (!is.null(rel_tab) && nrow(rel_tab) > 0) DT::datatable(rel_tab, rownames = FALSE, options = list(dom = "t", paging = FALSE), class = "stripe hover compact")
      else div(class = "empty-note", icon("circle-info"), "No RefGene-annotated gene overlaps this CpG (intergenic)."),
      tags$hr(), tags$b("Primary gene detail"),
      bc_primary_gene_detail(d)
  )
}

bc_ext_not_fetched <- function(title, icon_name) {
  div(class = "card", div(class = "card-title", icon(icon_name), title),
      div(class = "empty-note", icon("circle-info"), "Not yet run - click Run below."))
}

bc_section_ewascatalog <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("EWAS Catalog", "table-list"))
  if (!isTRUE(ext$ok)) return(div(class = "card", div(class = "card-title", icon("table-list"), "EWAS Catalog"),
                                   div(class = "empty-note", icon("triangle-exclamation"), ext$reason %||% "EWAS Catalog lookup unavailable.")))
  body <- if (is.null(ext$df) || nrow(ext$df) == 0) div(class = "empty-note", icon("circle-info"), "No EWAS Catalog trait associations found for this CpG.")
          else DT::datatable(ext$df, rownames = FALSE, filter = "top", options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
  div(class = "card", div(class = "card-title", icon("table-list"), "EWAS Catalog"),
      p(class = "submodule-desc", "MRC-IEU EWAS Catalog trait associations for this CpG."),
      body)
}

bc_section_ewasatlas <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("EWAS Atlas", "table-list"))
  if (!isTRUE(ext$ok)) return(div(class = "card", div(class = "card-title", icon("table-list"), "EWAS Atlas"),
                                   div(class = "empty-note", icon("triangle-exclamation"), ext$reason %||% "EWAS Atlas lookup unavailable.")))
  assoc <- ext$data$associationList
  body <- if (is.null(assoc) || !is.data.frame(assoc) || nrow(assoc) == 0) div(class = "empty-note", icon("circle-info"), "No EWAS Atlas trait associations found for this CpG.")
          else DT::datatable(assoc, rownames = FALSE, filter = "top", options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
  div(class = "card", div(class = "card-title", icon("table-list"), "EWAS Atlas"),
      p(class = "submodule-desc", "NGDC EWAS Atlas trait associations for this CpG."),
      body)
}

bc_section_disease_evidence <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("Disease Evidence", "notes-medical"))
  cmb <- ext$combined_disease
  body <- if (is.null(cmb) || nrow(cmb) == 0) {
    div(class = "empty-note", icon("circle-info"), "No evidence identified in the connected databases (EWAS Catalog, EWAS Atlas) for this CpG.")
  } else {
    disp <- cmb[order(cmb$source, cmb$trait), , drop = FALSE]
    DT::datatable(disp, colnames = c("Source", "Trait", "Tissue", "Effect (beta/correlation)", "P-value", "PMID"),
                  rownames = FALSE, filter = "top", options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
  }
  div(class = "card", div(class = "card-title", icon("notes-medical"), "Disease Evidence"),
      p(class = "submodule-desc", "Source: MRC-IEU EWAS Catalog + EWAS Atlas (this CpG)."),
      body)
}

bc_section_ra_evidence <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("Rheumatoid Arthritis Evidence", "hand-dots"))
  ra <- ext$ra_rows
  found <- !is.null(ra) && nrow(ra) > 0
  banner <- if (found) div(class = "empty-note", style = "border-left-color:#0ca30c;", icon("circle-check"), tags$b("RA-associated biomarker"))
            else div(class = "empty-note", icon("circle-info"), tags$b("No RA-specific evidence identified in the connected databases."))
  body <- if (found) DT::datatable(ra[, c("source", "trait", "effect", "p", "pmid")], colnames = c("Source", "Trait", "Effect", "P-value", "PMID"),
                                    rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact") else NULL
  div(class = "card", div(class = "card-title", icon("hand-dots"), "Rheumatoid Arthritis Evidence"),
      p(class = "submodule-desc", "Source: MRC-IEU EWAS Catalog + EWAS Atlas, filtered to rheumatoid-arthritis-related traits (this CpG)."),
      banner, body)
}

bc_section_disease_comparison <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("Disease Comparison", "chart-simple"))
  dc <- ext$disease_counts
  if (is.null(dc) || nrow(dc) == 0) {
    return(div(class = "card", div(class = "card-title", icon("chart-simple"), "Disease Comparison"),
               div(class = "empty-note", icon("circle-info"), "No disease/trait associations were found to compare.")))
  }
  dc <- dc[order(-dc$n), , drop = FALSE]
  dc$category <- factor(dc$category, levels = rev(dc$category))
  p <- ggplot(dc, aes(x = n, y = category)) + geom_col(fill = ARTHOMIX_COLORS$blue) +
    labs(x = "Real associations found (EWAS Catalog + EWAS Atlas)", y = NULL) + theme_arthomix()
  div(class = "card", div(class = "card-title", icon("chart-simple"), "Disease Comparison"),
      p(class = "submodule-desc", "Source: MRC-IEU EWAS Catalog + EWAS Atlas trait associations for this CpG, bucketed by disease category. Bar length is a real association count, not a synthetic score."),
      renderPlot_static(p))
}

bc_section_tissue_evidence <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("Tissue / Cell-Type Evidence", "microscope"))
  tc <- ext$tissue_counts
  body <- if (is.null(tc) || nrow(tc) == 0) div(class = "empty-note", icon("circle-info"), "No tissue information was returned by the connected databases for this CpG.")
          else DT::datatable(tc[order(-tc$Freq), , drop = FALSE], colnames = c("Tissue", "Studies"), rownames = FALSE,
                              options = list(dom = "t", paging = FALSE), class = "stripe hover compact")
  div(class = "card", div(class = "card-title", icon("microscope"), "Tissue / Cell-Type Evidence"),
      p(class = "submodule-desc", "Source: MRC-IEU EWAS Catalog (tissue field, this CpG)."),
      body)
}

bc_section_kegg <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("KEGG Pathways", "diagram-project"))
  kegg <- ext$kegg
  body <- if (isTRUE(kegg$ok) && nrow(kegg$pathways) > 0) DT::datatable(kegg$pathways, colnames = c("KEGG ID", "Pathway name"), rownames = FALSE,
                                                                         options = list(pageLength = 5, scrollX = TRUE), class = "stripe hover compact")
          else div(class = "empty-note", icon("circle-info"), if (isTRUE(kegg$ok)) "No KEGG pathway annotation found for this gene." else (kegg$reason %||% "KEGG lookup unavailable."))
  div(class = "card", div(class = "card-title", icon("diagram-project"), "KEGG Pathways"),
      p(class = "submodule-desc", "Database annotation (membership lookup for the associated gene) - not a computed enrichment test, which needs a gene set and a background to test against."),
      body)
}

bc_section_reactome <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("Reactome Pathways", "route"))
  reactome <- ext$reactome
  body <- if (isTRUE(reactome$ok) && nrow(reactome$pathways) > 0) DT::datatable(reactome$pathways, colnames = c("Reactome Stable ID", "Pathway name"), rownames = FALSE,
                                                                                 options = list(pageLength = 5, scrollX = TRUE), class = "stripe hover compact")
          else div(class = "empty-note", icon("circle-info"), if (isTRUE(reactome$ok)) "No Reactome pathway annotation found for this gene." else (reactome$reason %||% "Reactome lookup unavailable."))
  div(class = "card", div(class = "card-title", icon("route"), "Reactome Pathways"),
      p(class = "submodule-desc", "Database annotation (membership lookup for the associated gene) - not a computed enrichment test, which needs a gene set and a background to test against."),
      body)
}

bc_section_wikipathways <- function(ext) {
  if (is.null(ext)) return(div(class = "card", div(class = "card-title", icon("map"), "WikiPathways"),
                                div(class = "empty-note", icon("circle-info"), "Not yet looked up - click Run below.")))
  wp <- ext$wikipathways
  body <- if (!isTRUE(wp$ok)) div(class = "empty-note", icon("circle-info"), wp$reason %||% "WikiPathways lookup unavailable.")
          else if (is.null(wp$pathways) || nrow(wp$pathways) == 0) div(class = "empty-note", icon("circle-info"), "No WikiPathways pathways found for this gene.")
          else {
            wdf <- wp$pathways
            wdf$Link <- sprintf('<a href="https://www.wikipathways.org/instance/%s" target="_blank" rel="noopener">%s</a>', wdf$id, wdf$id)
            DT::datatable(wdf[, c("Link", "name")], colnames = c("WikiPathways ID", "Pathway"), rownames = FALSE, escape = 1,
                          options = list(dom = "t", paging = FALSE), class = "stripe hover compact")
          }
  div(class = "card", div(class = "card-title", icon("map"), "WikiPathways"),
      p(class = "submodule-desc", "Community-curated pathway diagrams, via the same gene-set source (msigdbr) the Multi-Omics Pathways module uses."),
      body)
}

bc_section_kegg_ra_pathway <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("KEGG RA Pathway (hsa05323)", "network-wired"))
  kr <- ext$kegg_ra
  if (!isTRUE(kr$ok)) {
    return(div(class = "card", div(class = "card-title", icon("network-wired"), "KEGG RA Pathway (hsa05323)"),
               div(class = "empty-note", icon("triangle-exclamation"), kr$reason %||% "KEGG RA pathway lookup unavailable.")))
  }
  banner <- if (isTRUE(kr$in_pathway)) div(class = "empty-note", style = "border-left-color:#0ca30c;", icon("circle-check"),
                                            tags$b(sprintf("%s is a member gene of the KEGG Rheumatoid Arthritis pathway (hsa05323).", ext$gene_symbol %||% "This gene")))
             else div(class = "empty-note", icon("circle-info"), tags$b(sprintf("%s is not listed as a member gene of the KEGG Rheumatoid Arthritis pathway (hsa05323).", ext$gene_symbol %||% "This gene")))
  div(class = "card", div(class = "card-title", icon("network-wired"), "KEGG RA Pathway (hsa05323)"),
      p(class = "submodule-desc", "Source: KEGG (pathway hsa05323 gene list)."),
      banner,
      p(style = "margin-top:8px;", tags$b(kr$pathway_name)),
      p(style = "font-size:12.5px; color:var(--color-ink-secondary);", kr$pathway_description),
      tags$details(tags$summary(sprintf("All %d genes on this pathway map", nrow(kr$genes))),
                   DT::datatable(kr$genes, colnames = c("NCBI Gene ID", "Symbol"), rownames = FALSE, filter = "top",
                                 options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact"))
  )
}

bc_section_publications <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("Supporting Publications", "book"))
  pubs <- ext$publications
  body <- if (!isTRUE(pubs$ok)) div(class = "empty-note", icon("triangle-exclamation"), pubs$reason %||% "PubMed lookup unavailable.")
          else if (is.null(pubs$table) || nrow(pubs$table) == 0) div(class = "empty-note", icon("circle-info"), "No publication PMIDs were returned by the connected databases for this CpG.")
          else tags$ul(lapply(seq_len(nrow(pubs$table)), function(i) {
            r <- pubs$table[i, ]
            tags$li(sprintf("%s (%s). %s. ", r$title %||% "Untitled", r$year %||% "n.d.", r$journal %||% ""),
                    tags$a(href = sprintf("https://pubmed.ncbi.nlm.nih.gov/%s/", r$pmid), target = "_blank", sprintf("PMID: %s", r$pmid)))
          }))
  div(class = "card", div(class = "card-title", icon("book"), "Supporting Publications"),
      p(class = "submodule-desc", "Source: PMIDs from EWAS Catalog / EWAS Atlas results, resolved via NCBI PubMed."),
      body)
}

bc_section_methbank <- function(ext) {
  if (is.null(ext)) return(bc_ext_not_fetched("MethBank", "link"))
  div(class = "card", div(class = "card-title", icon("link"), "MethBank"),
      div(class = "empty-note", icon("circle-info"),
          "MethBank has no public API in this deployment - it is referenced here as an outbound link only.",
          tags$br(), tags$a(href = ext$methbank_link, target = "_blank", "Search MethBank manually")))
}

renderPlot_static <- function(p) {
  tf <- tempfile(fileext = ".png")
  ok <- tryCatch({ ggplot2::ggsave(tf, plot = p, width = 7, height = 3, dpi = 110, bg = "white"); TRUE }, error = function(e) FALSE)
  if (!ok || !requireNamespace("base64enc", quietly = TRUE)) return(div(class = "empty-note", "Plot unavailable."))
  uri <- base64enc::dataURI(file = tf, mime = "image/png")
  unlink(tf)
  tags$img(src = uri, style = "max-width:100%; height:auto;")
}

bc_section_sex_specific <- function(d, plot_tag) {
  mkrow <- function(label, r) {
    if (is.null(r)) return(data.frame(Analysis = label, `Delta Beta` = NA_real_, `P-value` = NA_real_, FDR = NA_real_, Direction = NA_character_, Significant = NA_character_, check.names = FALSE, stringsAsFactors = FALSE))
    data.frame(Analysis = label, `Delta Beta` = round(r$dbeta, 4), `P-value` = signif(r$p_bacon, 4), FDR = signif(r$fdr_bacon, 4),
               Direction = if (is.na(r$dbeta)) NA_character_ else if (r$dbeta > 0) "Hyper" else if (r$dbeta < 0) "Hypo" else "None",
               Significant = if (is.na(r$fdr_bacon)) NA_character_ else if (r$fdr_bacon <= 0.05) "Yes" else "No",
               check.names = FALSE, stringsAsFactors = FALSE)
  }
  pipeline_df <- rbind(mkrow("Female (pipeline)", d$dmp_female), mkrow("Male (pipeline)", d$dmp_male))

  ov <- if (isTRUE(d$live$ok)) d$live$overall else NULL
  overall_row <- if (!is.null(ov) && isTRUE(ov$ok)) {
    data.frame(Analysis = "Overall (your data, computed)", `Delta Beta` = round(ov$dbeta, 4), `P-value` = signif(ov$p_value, 4), FDR = NA_real_,
               Direction = ov$direction, Significant = if (!is.na(ov$p_value) && ov$p_value <= 0.05) "Yes (uncorrected p)" else "No",
               check.names = FALSE, stringsAsFactors = FALSE)
  } else data.frame(Analysis = "Overall (your data, computed)", `Delta Beta` = NA_real_, `P-value` = NA_real_, FDR = NA_real_, Direction = NA_character_, Significant = NA_character_, check.names = FALSE, stringsAsFactors = FALSE)

  live_rows <- NULL
  if (isTRUE(d$live$ok) && !is.null(d$live$by_sex)) {
    live_rows <- do.call(rbind, lapply(names(d$live$by_sex), function(s) {
      r <- d$live$by_sex[[s]]
      if (!isTRUE(r$ok)) return(data.frame(Analysis = sprintf("%s (your data, computed)", s), `Delta Beta` = NA_real_, `P-value` = NA_real_, FDR = NA_real_, Direction = NA_character_, Significant = NA_character_, check.names = FALSE, stringsAsFactors = FALSE))
      data.frame(Analysis = sprintf("%s (your data, computed)", s), `Delta Beta` = round(r$dbeta, 4), `P-value` = signif(r$p_value, 4), FDR = NA_real_,
                 Direction = r$direction, Significant = if (!is.na(r$p_value) && r$p_value <= 0.05) "Yes (uncorrected p)" else "No",
                 check.names = FALSE, stringsAsFactors = FALSE)
    }))
  }
  all_df <- rbind(pipeline_df, overall_row, if (!is.null(live_rows)) live_rows)

  f_sig <- !is.null(d$dmp_female) && !is.na(d$dmp_female$fdr_bacon) && d$dmp_female$fdr_bacon <= 0.05
  m_sig <- !is.null(d$dmp_male) && !is.na(d$dmp_male$fdr_bacon) && d$dmp_male$fdr_bacon <= 0.05
  verdict <- if (f_sig && m_sig) "Significant in both sexes (pipeline, FDR <= 0.05)."
             else if (f_sig && !m_sig) "Female-specific signal (pipeline, FDR <= 0.05 in females only)."
             else if (m_sig && !f_sig) "Male-specific signal (pipeline, FDR <= 0.05 in males only)."
             else if (!is.null(d$dmp_female) && !is.null(d$dmp_male)) "Not significant in either sex-stratified pipeline analysis (FDR <= 0.05)."
             else "Sex-stratified pipeline results are not both available for this CpG."

  div(class = "card",
      div(class = "card-title", icon("venus-mars"), "Sex-Specific Evidence"),
      p(class = "submodule-desc", "\"Your data\" rows are unadjusted - not comparable to pipeline FDR."),
      div(class = "empty-note", icon("circle-info"), tags$b(verdict)),
      DT::datatable(all_df, rownames = FALSE, options = list(dom = "t", paging = FALSE), class = "stripe hover compact"),
      tags$div(style = "margin-top:10px;", plot_tag)
  )
}

bc_section_interpretation <- function(d, ext = NULL) {
  gene_label <- d$primary_gene %||% "Unknown gene"
  region_label <- paste0(bc_feature_label(d$resolved$gene_group), if (!is.null(d$island)) " / CpG island" else "")
  cand <- list(Female = d$dmp_female, Male = d$dmp_male); cand <- cand[!vapply(cand, is.null, logical(1))]
  dir_label <- "Not available"
  if (length(cand) > 0) {
    fdrs <- vapply(cand, function(x) if (is.na(x$fdr_bacon)) Inf else x$fdr_bacon, numeric(1))
    best <- cand[[names(cand)[which.min(fdrs)]]]
    if (!is.na(best$dbeta)) dir_label <- if (best$dbeta > 0) "Hypermethylated" else if (best$dbeta < 0) "Hypomethylated" else "No change"
  }
  pathway_label <- "Not yet connected"; pathway_done <- FALSE
  disease_label <- "Not yet connected"; disease_done <- FALSE
  ra_label <- "Not yet connected"; ra_done <- FALSE
  if (!is.null(ext)) {
    n_kegg <- if (isTRUE(ext$kegg$ok)) nrow(ext$kegg$pathways) else 0
    n_reactome <- if (isTRUE(ext$reactome$ok)) nrow(ext$reactome$pathways) else 0
    pathway_label <- if (n_kegg + n_reactome > 0) sprintf("%d pathway(s) found", n_kegg + n_reactome) else "None found"
    pathway_done <- (n_kegg + n_reactome) > 0
    n_disease <- if (!is.null(ext$combined_disease)) nrow(ext$combined_disease) else 0
    disease_label <- if (n_disease > 0) sprintf("%d association(s) found", n_disease) else "None found"
    disease_done <- n_disease > 0
    ra_terms <- if (!is.null(ext$ra_rows) && nrow(ext$ra_rows) > 0) unique(ext$ra_rows$trait) else character(0)
    ra_in_kegg <- isTRUE(ext$kegg_ra$ok) && isTRUE(ext$kegg_ra$in_pathway)
    ra_done <- length(ra_terms) > 0 || ra_in_kegg
    ra_label <- if (length(ra_terms) > 0 && ra_in_kegg) sprintf("%s + KEGG RA pathway", ra_terms[1])
                else if (length(ra_terms) > 0) ra_terms[1]
                else if (ra_in_kegg) "In KEGG RA pathway (hsa05323)"
                else "No RA-specific evidence found"
  }
  nodes <- list(
    list(label = "CpG", value = d$cpg, done = TRUE),
    list(label = "Region", value = region_label, done = TRUE),
    list(label = "Gene", value = gene_label, done = !is.na(d$primary_gene %||% NA)),
    list(label = "Methylation", value = dir_label, done = !identical(dir_label, "Not available")),
    list(label = "Pathway", value = pathway_label, done = pathway_done),
    list(label = "Disease", value = disease_label, done = disease_done),
    list(label = "RA Evidence", value = ra_label, done = ra_done)
  )
  div(class = "card",
      div(class = "card-title", icon("diagram-project"), "Functional Interpretation"),
      p(class = "submodule-desc", "Cautious, evidence-only chain - each arrow means \"associated with\" / \"located within\", never a causal claim."),
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

bc_section_evidence_summary <- function(d, ext = NULL) {
  user_sig <- FALSE
  if (!is.null(d$dmp_female) && !is.na(d$dmp_female$fdr_bacon) && d$dmp_female$fdr_bacon <= 0.05) user_sig <- TRUE
  if (!is.null(d$dmp_male) && !is.na(d$dmp_male$fdr_bacon) && d$dmp_male$fdr_bacon <= 0.05) user_sig <- TRUE
  if (isTRUE(d$live$ok) && isTRUE(d$live$overall$ok) && !is.na(d$live$overall$p_value) && d$live$overall$p_value <= 0.05) user_sig <- TRUE
  gene_assoc <- !is.na(d$primary_gene %||% NA)
  island_reg <- !is.null(d$island) || (!is.na(d$resolved$regulatory_group) && nzchar(d$resolved$regulatory_group))
  sex_specific <- !is.null(d$dmp_female) && !is.null(d$dmp_male)

  ewas_repl <- "Not yet checked"; ra_evid <- "Not yet checked"; pathway_assoc <- "Not yet checked"
  if (!is.null(ext)) {
    n_disease <- if (!is.null(ext$combined_disease)) nrow(ext$combined_disease) else 0
    ewas_repl <- if (n_disease > 0) "Yes" else "No"
    ra_evid <- if (!is.null(ext$ra_rows) && nrow(ext$ra_rows) > 0) "Yes" else "No"
    n_kegg <- if (isTRUE(ext$kegg$ok)) nrow(ext$kegg$pathways) else 0
    n_reactome <- if (isTRUE(ext$reactome$ok)) nrow(ext$reactome$pathways) else 0
    pathway_assoc <- if ((n_kegg + n_reactome) > 0) "Yes" else "No"
  }

  rows <- list(
    c("User dataset significant (FDR/P <= 0.05)", if (user_sig) "Yes" else "No"),
    c("Gene association", if (gene_assoc) "Yes" else "No"),
    c("CpG island / regulatory annotation", if (island_reg) "Yes" else "No"),
    c("Sex-specific evidence (both sexes tested)", if (sex_specific) "Yes" else "No"),
    c("Independent EWAS replication (EWAS Catalog / EWAS Atlas)", ewas_repl),
    c("RA evidence", ra_evid),
    c("Pathway association (KEGG / Reactome)", pathway_assoc),
    c("Methylation-expression evidence", "Not yet integrated")
  )
  df <- data.frame(Evidence = vapply(rows, `[`, character(1), 1), Status = vapply(rows, `[`, character(1), 2), stringsAsFactors = FALSE)
  div(class = "card",
      div(class = "card-title", icon("clipboard-check"), "Summary"),
      p(class = "submodule-desc", "Each row is a transparent, individually-computed check - not a combined score. Rows marked \"Not yet integrated\" are planned for a later phase, not evaluated as absent."),
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE, ordering = FALSE), class = "stripe hover compact")
  )
}

BC_NOT_INTEGRATED_DBS <- c(
  "DisGeNET - requires API registration/licensing not configured in this deployment (same call already made for the sibling transcriptomics Biomarker Card)",
  "GTEx - would duplicate the baseline tissue/blood expression already shown via Human Protein Atlas above",
  "OMIM, Human Phenotype Ontology - require registration/licensing not verified for this deployment"
)

bc_single_cpg_roc <- function(beta_row, group_vec, case_label, control_label) {
  if (!requireNamespace("pROC", quietly = TRUE)) return(list(ok = FALSE, reason = "pROC is not installed in this deployment."))
  keep <- !is.na(beta_row) & !is.na(group_vec) & group_vec %in% c(case_label, control_label)
  x <- beta_row[keep]; y <- group_vec[keep]
  n_case <- sum(y == case_label); n_control <- sum(y == control_label)
  if (n_case < 3 || n_control < 3) {
    return(list(ok = FALSE, reason = sprintf("Fewer than 3 samples in one of the two groups (case=%d, control=%d) - cannot compute a reliable single-CpG ROC curve.", n_case, n_control)))
  }
  roc_obj <- tryCatch(pROC::roc(y, x, levels = c(control_label, case_label), direction = "auto", quiet = TRUE), error = function(e) NULL)
  if (is.null(roc_obj)) return(list(ok = FALSE, reason = "pROC could not fit a ROC curve for this CpG (e.g. constant methylation)."))
  ci <- tryCatch(as.numeric(pROC::ci.auc(roc_obj, quiet = TRUE)), error = function(e) c(NA_real_, NA_real_, NA_real_))
  best <- tryCatch(pROC::coords(roc_obj, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = FALSE), error = function(e) NULL)
  if (is.null(best) || nrow(best) == 0) return(list(ok = FALSE, reason = "Could not determine an optimal threshold for this CpG."))
  best <- best[1, ]
  pred_case <- if (identical(roc_obj$direction, "<")) x >= best$threshold else x <= best$threshold
  obs_case <- y == case_label
  tp <- sum(pred_case & obs_case); fn <- sum(!pred_case & obs_case)
  tn <- sum(!pred_case & !obs_case); fp <- sum(pred_case & !obs_case)
  list(ok = TRUE, auc = as.numeric(pROC::auc(roc_obj)), ci_lo = ci[1], ci_hi = ci[3],
       threshold = best$threshold, sensitivity = best$sensitivity, specificity = best$specificity,
       accuracy = best$accuracy, ppv = best$ppv, npv = best$npv, balanced_accuracy = (best$sensitivity + best$specificity) / 2,
       n_case = n_case, n_control = n_control,
       confusion = matrix(c(tp, fp, fn, tn), 2, 2, dimnames = list(Predicted = c("Case", "Control"), Actual = c("Case", "Control"))),
       roc_obj = roc_obj, case_label = case_label, control_label = control_label)
}

bc_single_cpg_cv <- function(beta_row, group_vec, case_label, control_label, k = 5, seed = 1234) {
  if (!requireNamespace("pROC", quietly = TRUE)) return(list(ok = FALSE, reason = "pROC is not installed in this deployment."))
  keep <- !is.na(beta_row) & !is.na(group_vec) & group_vec %in% c(case_label, control_label)
  x <- beta_row[keep]; y <- group_vec[keep]
  n_case <- sum(y == case_label); n_control <- sum(y == control_label)
  k_eff <- min(k, n_case, n_control)
  if (k_eff < 3) return(list(ok = FALSE, reason = sprintf("Not enough samples per group for cross-validation (case=%d, control=%d; need >=3 in each).", n_case, n_control)))
  set.seed(seed)
  fold <- integer(length(y))
  fold[y == case_label] <- sample(rep(seq_len(k_eff), length.out = n_case))
  fold[y == control_label] <- sample(rep(seq_len(k_eff), length.out = n_control))
  oof_score <- rep(NA_real_, length(y)); oof_pred <- rep(NA_character_, length(y))
  for (f in seq_len(k_eff)) {
    te <- fold == f; tr <- !te
    if (sum(tr & y == case_label) < 2 || sum(tr & y == control_label) < 2) next
    roc_tr <- tryCatch(pROC::roc(y[tr], x[tr], levels = c(control_label, case_label), direction = "auto", quiet = TRUE), error = function(e) NULL)
    if (is.null(roc_tr)) next
    best_tr <- tryCatch(pROC::coords(roc_tr, "best", best.method = "youden", ret = "threshold", transpose = FALSE), error = function(e) NULL)
    if (is.null(best_tr) || nrow(best_tr) == 0) next
    thr <- best_tr$threshold[1]
    oof_score[te] <- x[te]
    oof_pred[te] <- if (identical(roc_tr$direction, "<")) ifelse(x[te] >= thr, case_label, control_label) else ifelse(x[te] <= thr, case_label, control_label)
  }
  usable <- !is.na(oof_pred)
  if (sum(usable) < 6) return(list(ok = FALSE, reason = "Cross-validation could not produce enough out-of-fold predictions (folds too small, or this CpG wasn't separable in any training fold)."))
  y_used <- y[usable]; pred_used <- oof_pred[usable]; score_used <- oof_score[usable]
  roc_pooled <- tryCatch(pROC::roc(y_used, score_used, levels = c(control_label, case_label), direction = "auto", quiet = TRUE), error = function(e) NULL)
  auc_pooled <- if (!is.null(roc_pooled)) as.numeric(pROC::auc(roc_pooled)) else NA_real_
  obs_case <- y_used == case_label; pred_case <- pred_used == case_label
  tp <- sum(pred_case & obs_case); fn <- sum(!pred_case & obs_case)
  tn <- sum(!pred_case & !obs_case); fp <- sum(pred_case & !obs_case)
  sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  list(ok = TRUE, k = k_eff, n_used = sum(usable), auc = auc_pooled, sensitivity = sens, specificity = spec,
       accuracy = (tp + tn) / (tp + tn + fp + fn), balanced_accuracy = (sens + spec) / 2,
       confusion = matrix(c(tp, fp, fn, tn), 2, 2, dimnames = list(Predicted = c("Case", "Control"), Actual = c("Case", "Control"))),
       roc_obj = roc_pooled, case_label = case_label, control_label = control_label)
}

bc_pr_from_roc <- function(roc_obj) {
  if (is.null(roc_obj)) return(list(ok = FALSE, reason = "No ROC curve available."))
  co <- tryCatch(pROC::coords(roc_obj, "all", ret = c("recall", "precision"), transpose = FALSE), error = function(e) NULL)
  if (is.null(co) || nrow(co) == 0) return(list(ok = FALSE, reason = "Could not derive a precision-recall curve from this ROC curve."))
  df <- data.frame(recall = co$recall, precision = co$precision)
  df <- df[order(df$recall), , drop = FALSE]
  keep <- !is.na(df$recall) & !is.na(df$precision)
  df <- df[keep, , drop = FALSE]
  if (nrow(df) < 2) return(list(ok = FALSE, reason = "Not enough non-missing precision/recall points to draw a curve."))
  pr_auc <- sum(diff(df$recall) * (utils::head(df$precision, -1) + utils::tail(df$precision, -1)) / 2)
  list(ok = TRUE, df = df, pr_auc = pr_auc)
}

bc_plot_roc <- function(roc_obj, label = "ROC curve", auc = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_) {
  if (is.null(roc_obj)) return(NULL)
  co <- pROC::coords(roc_obj, "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
  df <- data.frame(fpr = 1 - co$specificity, tpr = co$sensitivity)
  sub <- if (!is.na(ci_lo) && !is.na(ci_hi)) sprintf("AUC %.3f (95%% CI %.3f-%.3f)", auc, ci_lo, ci_hi) else sprintf("AUC %.3f", auc)
  ggplot(df, aes(x = fpr, y = tpr)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = ARTHOMIX_COLORS$grid) +
    geom_line(color = ARTHOMIX_COLORS$blue, linewidth = 1) +
    labs(title = label, subtitle = sub, x = "1 - specificity", y = "Sensitivity") +
    coord_equal() + theme_arthomix()
}

bc_plot_pr <- function(pr, label = "Precision-recall curve") {
  if (is.null(pr) || !isTRUE(pr$ok)) return(NULL)
  ggplot(pr$df, aes(x = recall, y = precision)) +
    geom_line(color = ARTHOMIX_COLORS$violet, linewidth = 1) +
    labs(title = label, subtitle = sprintf("PR-AUC %.3f", pr$pr_auc), x = "Recall (sensitivity)", y = "Precision (PPV)") +
    ylim(0, 1) + theme_arthomix()
}

bc_fmt_num <- function(x, digits = 3) if (is.null(x) || length(x) == 0 || is.na(x)) "Not available" else sprintf(sprintf("%%.%df", digits), x)

bc_confusion_table <- function(cm) {
  if (is.null(cm)) return(NULL)
  df <- as.data.frame.matrix(cm)
  df <- cbind(Predicted = rownames(df), df)
  DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE, ordering = FALSE), class = "stripe hover compact")
}

.bc_ensemble_votes_cache <- new.env(parent = emptyenv())
bc_ensemble_votes_table <- function(sex) {
  cached <- .bc_ensemble_votes_cache[[sex]]
  if (!is.null(cached)) return(cached)
  df <- load_default_diagnostic_ensemble_votes(sex)
  if (!is.null(df)) .bc_ensemble_votes_cache[[sex]] <- df
  df
}
.bc_panel_auc_cache <- new.env(parent = emptyenv())
bc_panel_auc_table <- function(sex) {
  cached <- .bc_panel_auc_cache[[sex]]
  if (!is.null(cached)) return(cached)
  df <- load_default_diagnostic_panel_auc(sex)
  if (!is.null(df)) .bc_panel_auc_cache[[sex]] <- df
  df
}
.bc_perprobe_auc_cache <- new.env(parent = emptyenv())
bc_perprobe_auc_table <- function(sex) {
  cached <- .bc_perprobe_auc_cache[[sex]]
  if (!is.null(cached)) return(cached)
  df <- load_default_diagnostic_perprobe_auc(sex)
  if (!is.null(df)) .bc_perprobe_auc_cache[[sex]] <- df
  df
}

bc_panel_membership_lookup <- function(cpg) {
  stats::setNames(lapply(c("female", "male"), function(sex) {
    votes <- bc_ensemble_votes_table(sex)
    row <- if (!is.null(votes)) votes[votes$cpg == cpg, , drop = FALSE] else NULL
    list(in_panel = !is.null(row) && nrow(row) > 0, row = if (!is.null(row) && nrow(row) > 0) row[1, ] else NULL)
  }), c("female", "male"))
}
bc_panel_perprobe_lookup <- function(cpg) {
  stats::setNames(lapply(c("female", "male"), function(sex) {
    pp <- bc_perprobe_auc_table(sex)
    row <- if (!is.null(pp)) pp[pp$cpg == cpg, , drop = FALSE] else NULL
    if (is.null(row) || nrow(row) == 0) NULL else row
  }), c("female", "male"))
}

bc_multi_cpg_perf_table <- function(membership) {
  strata <- Filter(function(s) isTRUE(membership[[s]]$in_panel), names(membership))
  if (length(strata) == 0) return(NULL)
  rows <- lapply(strata, function(sex) {
    auc <- bc_panel_auc_table(sex)
    if (is.null(auc) || nrow(auc) == 0) return(NULL)
    best <- auc[which.max(auc$internal_test_auc), ]
    votes <- bc_ensemble_votes_table(sex)
    data.frame(Stratum = tools::toTitleCase(sex), `Panel size` = if (!is.null(votes)) nrow(votes) else NA_integer_,
               `Best algorithm (by internal AUC)` = toupper(best$algorithm),
               `Training AUC` = bc_fmt_num(best$train_auc), `Internal validation AUC` = bc_fmt_num(best$internal_test_auc),
               `External validation AUC` = bc_fmt_num(best$external_test_auc),
               check.names = FALSE, stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

bc_section_signature_comparison <- function(d, sgd, roc_widget = NULL) {
  single_body <- if (isTRUE(sgd$ok)) {
    bc_kv_table(list(
      AUC = bc_fmt_num(sgd$auc), "95% CI" = sprintf("%s - %s", bc_fmt_num(sgd$ci_lo), bc_fmt_num(sgd$ci_hi)),
      Threshold = bc_fmt_num(sgd$threshold), Sensitivity = bc_fmt_num(sgd$sensitivity), Specificity = bc_fmt_num(sgd$specificity),
      Accuracy = bc_fmt_num(sgd$accuracy), `Balanced accuracy` = bc_fmt_num(sgd$balanced_accuracy),
      `n (case / control)` = sprintf("%d / %d", sgd$n_case, sgd$n_control)
    ))
  } else div(class = "empty-note", icon("circle-info"), sgd$reason %||% "Single-CpG performance unavailable.")

  membership <- d$panel_membership
  perprobe <- d$panel_perprobe
  membership_rows <- do.call(rbind, lapply(names(membership), function(sex) {
    m <- membership[[sex]]
    data.frame(Stratum = tools::toTitleCase(sex),
               `In majority-vote panel` = if (isTRUE(m$in_panel)) "Yes" else "No",
               `In LASSO` = if (isTRUE(m$in_panel)) as.character(m$row$in_LASSO) else NA_character_,
               `In Boruta` = if (isTRUE(m$in_panel)) as.character(m$row$in_Boruta) else NA_character_,
               `In SVM-RFE` = if (isTRUE(m$in_panel)) as.character(m$row$in_SVMRFE) else NA_character_,
               Votes = if (isTRUE(m$in_panel)) m$row$n_votes else NA_integer_,
               check.names = FALSE, stringsAsFactors = FALSE)
  }))
  perprobe_rows <- do.call(rbind, Filter(Negate(is.null), lapply(names(perprobe), function(sex) {
    r <- perprobe[[sex]]
    if (is.null(r)) return(NULL)
    data.frame(Stratum = tools::toTitleCase(sex), Algorithm = toupper(r$algorithm[1]),
               `Training AUC` = bc_fmt_num(r$train_auc[1]), `Internal validation AUC` = bc_fmt_num(r$internal_test_auc[1]),
               `External validation AUC` = bc_fmt_num(r$external_test_auc[1]),
               check.names = FALSE, stringsAsFactors = FALSE)
  })))

  multi_table <- bc_multi_cpg_perf_table(membership)
  multi_body <- if (is.null(multi_table)) {
    NULL
  } else DT::datatable(multi_table, rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")

  comparison <- NULL
  scv <- d$single_cpg_cv
  cv_available <- !is.null(scv) && isTRUE(scv$ok) && is.finite(scv$auc %||% NA_real_)
  single_auc <- if (cv_available) scv$auc else if (isTRUE(sgd$ok)) sgd$auc else NA_real_
  single_auc_label <- if (cv_available) "this CpG alone, cross-validated" else "this CpG alone, training AUC (CV unavailable)"
  if (!is.null(multi_table) && is.finite(single_auc)) {
    best_multi_internal <- suppressWarnings(max(as.numeric(gsub("[^0-9.]", "", multi_table$`Internal validation AUC`)), na.rm = TRUE))
    if (is.finite(best_multi_internal)) {
      diff <- best_multi_internal - single_auc
      comparison <- div(class = "empty-note", icon("scale-balanced"),
        if (abs(diff) < 0.01) sprintf("The multi-CpG panel's internal-validation AUC (%.3f) is essentially the same as %s (%.3f) - no meaningful signature benefit is shown by the numbers here.", best_multi_internal, single_auc_label, single_auc)
        else if (diff > 0) sprintf("The multi-CpG panel's internal-validation AUC (%.3f) is %.3f points higher than %s (%.3f) - the multi-CpG signature adds real, measured value here.", best_multi_internal, diff, single_auc_label, single_auc)
        else sprintf("This CpG alone (%s, AUC %.3f) actually scores %.3f points higher than the multi-CpG panel's internal-validation AUC (%.3f).", if (cv_available) "cross-validated" else "training AUC, CV unavailable", single_auc, -diff, best_multi_internal))
    }
  }

  tagList(
    div(class = "card",
        div(class = "card-title", icon("layer-group"), "Single-CpG Performance"),
        p(class = "submodule-desc", "Single-CpG classifier, full-fit on this dataset."),
        single_body,
        if (!is.null(roc_widget)) tagList(tags$div(style = "margin-top:10px;"), roc_widget) else NULL
    ),
    div(class = "card",
        div(class = "card-title", icon("layer-group"), "Multi-CpG Panel Membership"),
        DT::datatable(membership_rows, rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")
    ),
    div(class = "card",
        div(class = "card-title", icon("layer-group"), "Multi-CpG Panel Performance"),
        multi_body,
        if (!is.null(perprobe_rows)) tagList(tags$div(style = "margin-top:10px;", tags$b("This CpG's own performance within the trained panel model")),
          DT::datatable(perprobe_rows, rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")) else NULL
    ),
    if (!is.null(comparison)) div(class = "card", div(class = "card-title", icon("scale-balanced"), "Comparison"), comparison) else NULL
  )
}

bc_section_biomarker_performance <- function(d, sgd, sgcv, train_roc_widget = NULL, train_pr_widget = NULL, internal_roc_widget = NULL) {
  train_body <- if (isTRUE(sgd$ok)) {
    tagList(
      bc_kv_table(list(AUC = bc_fmt_num(sgd$auc), "95% CI" = sprintf("%s - %s", bc_fmt_num(sgd$ci_lo), bc_fmt_num(sgd$ci_hi)),
                        Threshold = bc_fmt_num(sgd$threshold), Sensitivity = bc_fmt_num(sgd$sensitivity), Specificity = bc_fmt_num(sgd$specificity),
                        Accuracy = bc_fmt_num(sgd$accuracy), `Balanced accuracy` = bc_fmt_num(sgd$balanced_accuracy),
                        PPV = bc_fmt_num(sgd$ppv), NPV = bc_fmt_num(sgd$npv), `n (case / control)` = sprintf("%d / %d", sgd$n_case, sgd$n_control))),
      tags$div(style = "margin-top:10px;", tags$b("Confusion matrix (single-CpG, at optimal threshold)")), bc_confusion_table(sgd$confusion),
      if (!is.null(train_roc_widget)) tagList(tags$div(style = "margin-top:10px;", tags$b("ROC curve")), train_roc_widget) else NULL,
      if (!is.null(train_pr_widget)) tagList(tags$div(style = "margin-top:10px;", tags$b("Precision-recall curve")), train_pr_widget) else NULL
    )
  } else div(class = "empty-note", icon("circle-info"), sgd$reason %||% "Not available.")

  internal_body <- if (isTRUE(sgcv$ok)) {
    tagList(
      bc_kv_table(list(AUC = bc_fmt_num(sgcv$auc), Sensitivity = bc_fmt_num(sgcv$sensitivity), Specificity = bc_fmt_num(sgcv$specificity),
                        Accuracy = bc_fmt_num(sgcv$accuracy), `Balanced accuracy` = bc_fmt_num(sgcv$balanced_accuracy),
                        `Folds (k)` = sgcv$k, `n used` = sgcv$n_used)),
      tags$div(style = "margin-top:10px;", tags$b("Confusion matrix (single-CpG, pooled out-of-fold predictions)")), bc_confusion_table(sgcv$confusion),
      if (!is.null(internal_roc_widget)) tagList(tags$div(style = "margin-top:10px;", tags$b("ROC curve (pooled out-of-fold)")), internal_roc_widget) else NULL
    )
  } else div(class = "empty-note", icon("circle-info"), sgcv$reason %||% "Not available.")

  membership <- d$panel_membership
  multi_table <- bc_multi_cpg_perf_table(membership)
  has_external <- !is.null(multi_table) && any(!is.na(suppressWarnings(as.numeric(gsub("[^0-9.]", "", multi_table$`External validation AUC`)))))
  external_body <- if (has_external) {
    tagList(
      p(class = "submodule-desc", "External cohort GSE111942 (preloaded, held-out)."),
      DT::datatable(multi_table[, c("Stratum", "Best algorithm (by internal AUC)", "External validation AUC")], rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")
    )
  } else div(class = "empty-note", icon("triangle-exclamation"),
    "Not in the external-validation panel. See the Diagnostic Classifier tab.")

  div(class = "card",
      div(class = "card-title", icon("vial-circle-check"), "Biomarker Performance"),
      p(class = "submodule-desc", "Single-CpG: computed here. Multi-CpG panel: preloaded, read verbatim."),
      tags$h5("Training"), train_body,
      tags$hr(),
      tags$h5("Internal Validation (cross-validation)"), internal_body,
      tags$hr(),
      tags$h5("External Validation"), external_body,
      tags$hr(),
      tags$b("Calibration"), div(class = "empty-note", icon("circle-info"), "Not applicable (single raw value)."),
      if (!is.null(multi_table)) tagList(tags$hr(), tags$b("Multi-CpG panel context (this CpG is part of the evaluated panel)"),
        DT::datatable(multi_table, rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")) else NULL
  )
}

bc_section_dataset_cohort <- function(dataset, live) {
  n_probes <- if (is.null(dataset$beta)) NA_integer_ else nrow(dataset$beta)
  n_samples <- if (is.null(dataset$beta)) NA_integer_ else ncol(dataset$beta)
  scale_label <- switch(dataset$input_scale %||% "unknown",
    beta = "Beta values (0-1 scale)", m = "M-values (logit-transformed)", "Not determined")
  group_label <- if (isTRUE(live$ok) || !is.null(live$group_col)) live$group_col else "Not auto-detected"
  case_label <- if (isTRUE(live$ok)) sprintf("%s vs %s", live$case_label, live$control_label) else "Not available"
  sex_label <- if (isTRUE(live$ok) && !is.null(live$sex_col)) live$sex_col else "Not recorded"
  n_case <- if (isTRUE(live$ok) && isTRUE(live$overall$ok)) live$overall$n_case else NA
  n_ctrl <- if (isTRUE(live$ok) && isTRUE(live$overall$ok)) live$overall$n_control else NA
  sex_breakdown <- if (isTRUE(live$ok) && !is.null(live$sex_vec)) {
    tb <- table(live$sex_vec[!is.na(live$sex_vec)])
    paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = ", ")
  } else NA_character_
  pairs <- list(
    "Dataset / cohort label" = dataset$source,
    "Array type" = dataset$array_type,
    "Probes x samples" = if (!is.na(n_probes)) sprintf("%s x %s", format(n_probes, big.mark = ","), n_samples) else NA,
    "Number of samples" = n_samples,
    "Data type / preprocessing status" = scale_label,
    "Species" = "Not available (not recorded for this dataset)",
    "Tissue / sample type" = "Not available (not recorded for this dataset)",
    "Analysis group / contrast" = sprintf("%s (grouping column: %s)", case_label, group_label),
    "Cases / Controls" = if (!is.na(n_case)) sprintf("%s / %s", n_case, n_ctrl) else NA,
    "Sex / group breakdown" = sex_breakdown,
    "Dataset designation" = "Training / discovery dataset (this session's loaded dataset). Internal validation = cross-validation folds within it; external validation = the preloaded panel's own external cohort (see Biomarker Performance tab), not a separate cohort loaded this session.",
    "Feature identifier" = "CpG probe ID (rows of the loaded beta matrix)",
    "Gene mapping" = "Illumina manifest annotation / ChAMPdata - see CpG description tab"
  )
  div(class = "card",
      div(class = "card-title", icon("table-cells"), "Dataset & Cohort"),
      p(class = "submodule-desc", "Where this biomarker's evidence comes from. Fields with no equivalent recorded anywhere in this app's dataset metadata are shown as \"Not available\", never guessed."),
      bc_kv_table(pairs)
  )
}

bc_section_sources <- function(ext = NULL) {
  used <- c(
    "Illumina manifest annotation (IlluminaHumanMethylation450kanno.ilmn12.hg19 / IlluminaHumanMethylationEPICanno.ilm10b4.hg19) - UCSC hg19-based probe location, CpG island, RefGene, and regulatory-feature tracks",
    "ChAMPdata::probe.features - UCSC hg19-based gene/feature/CpG-island annotation",
    "TxDb.Hsapiens.UCSC.hg19.knownGene - UCSC knownGene transcript/exon structure (hg19)",
    "org.Hs.eg.db - NCBI Gene ID, Ensembl Gene ID, gene name, GO term mapping",
    "GO.db - Gene Ontology term names",
    "UCSC cytoBandIdeo (hg19) - chromosome ideogram/cytoband track, bundled from hgdownload.soe.ucsc.edu",
    "This app's own preloaded sex-stratified SVA/bacon-corrected DMP pipeline results, and/or your uploaded dataset"
  )
  external_available <- c(
    "MRC-IEU EWAS Catalog - https://www.ewascatalog.org/api/",
    "EWAS Atlas - https://ngdc.cncb.ac.cn/ewas/rest/",
    "NCBI Gene - https://eutils.ncbi.nlm.nih.gov/ (esummary, db=gene)",
    "Ensembl (GRCh37/hg19) - https://grch37.rest.ensembl.org/",
    "Ensembl Regulatory Build - query by CpG coordinate, https://grch37.rest.ensembl.org/overlap/region/",
    "KEGG - via KEGGREST, https://rest.kegg.jp/",
    "Reactome - https://reactome.org/ContentService/",
    "WikiPathways - community-curated pathways (msigdbr C2:CP:WIKIPATHWAYS), https://www.wikipathways.org/",
    "Open Targets - genetic association (aggregates GWAS Catalog) + druggability, https://api.platform.opentargets.org/",
    "GWAS Catalog (NHGRI-EBI) - trait-name search, https://www.ebi.ac.uk/gwas/",
    "Human Protein Atlas - baseline tissue / blood-lineage RNA expression, https://www.proteinatlas.org/",
    "ENCODE - experiment/dataset search, https://www.encodeproject.org/",
    "GEO (NCBI E-utilities, db=gds) - search/summary only, no automatic download, https://www.ncbi.nlm.nih.gov/gds/",
    "BioStudies / ArrayExpress (EBI) - search, https://www.ebi.ac.uk/biostudies/",
    "PubMed (NCBI E-utilities) - https://eutils.ncbi.nlm.nih.gov/"
  )
  external_label <- if (is.null(ext)) "Available (pick a subtab and click its lookup button to query):" else "Used in this Biomarker Card (external, queried):"
  div(class = "card",
      div(class = "card-title", icon("database"), "Database Sources"),
      tags$b("Used in this Biomarker Card (local/offline):"), tags$ul(lapply(used, tags$li)),
      tags$b(external_label), tags$ul(lapply(external_available, tags$li)),
      tags$b("Referenced, no public API (outbound link/browser view only):"),
      tags$ul(tags$li("MethBank"), tags$li("UCSC Genome Browser (real track view, deep-linked at the CpG's coordinate)")),
      tags$div(style = "margin-top:10px;", tags$b("Considered but not integrated (avoided rather than faked - see reasons):")),
      tags$ul(lapply(BC_NOT_INTEGRATED_DBS, tags$li))
  )
}

bc_db_provenance <- function(api_domain, live_url, live_label) {
  div(class = "empty-note", style = "display:flex; align-items:center; justify-content:space-between; gap:10px; flex-wrap:wrap;",
      tagList(icon("satellite-dish"), sprintf("Query to %s - nothing here is precomputed or cached in this app.", api_domain)),
      tags$a(href = live_url, target = "_blank", rel = "noopener", class = "btn btn-sm btn-default",
             icon("arrow-up-right-from-square"), sprintf(" Open %s ", live_label))
  )
}

bc_ucsc_browser_link <- function(chr, pos, flank = 200) {
  if (is.null(chr) || is.na(chr) || is.null(pos) || is.na(pos)) return(NULL)
  sprintf("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=%s:%s-%s", chr,
          format(round(pos - flank), scientific = FALSE, trim = TRUE), format(round(pos + flank), scientific = FALSE, trim = TRUE))
}

bc_section_gene_genome <- function(d, ncbi, ensembl) {
  ucsc_link <- bc_ucsc_browser_link(d$resolved$chr, d$resolved$pos)
  ncbi_prov <- bc_db_provenance("eutils.ncbi.nlm.nih.gov", sprintf("https://www.ncbi.nlm.nih.gov/gene/%s", d$gene_struct$entrez %||% ""), "on NCBI Gene")
  ncbi_body <- if (is.null(ncbi)) div(class = "empty-note", icon("circle-info"), "Not yet looked up - click \"Look Up Gene Study Evidence\" below.")
               else if (!isTRUE(ncbi$ok)) div(class = "empty-note", icon("triangle-exclamation"), ncbi$reason %||% "NCBI Gene lookup unavailable.")
               else bc_kv_table(list("Official full name" = ncbi$description, "Aliases" = ncbi$aliases, "Chromosome" = ncbi$chromosome,
                                      "Cytogenetic map location" = ncbi$map_location, "Gene summary" = ncbi$summary))
  ens_prov <- bc_db_provenance("grch37.rest.ensembl.org", sprintf("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=%s", ensembl$ensembl_id %||% d$gene_struct$ensembl %||% ""), "on Ensembl (GRCh37)")
  ens_body <- if (is.null(ensembl)) div(class = "empty-note", icon("circle-info"), "Not yet looked up - click \"Look Up Gene Study Evidence\" below.")
              else if (!isTRUE(ensembl$ok)) div(class = "empty-note", icon("triangle-exclamation"), ensembl$reason %||% "Ensembl lookup unavailable.")
              else bc_kv_table(list(
                "Ensembl Gene ID" = ensembl$ensembl_id, "Biotype" = ensembl$biotype,
                "Genomic location (GRCh37/hg19)" = if (!is.na(ensembl$chr)) sprintf("chr%s:%s-%s (%s strand)", ensembl$chr, format(ensembl$start, big.mark = ","), format(ensembl$end, big.mark = ","), if (identical(ensembl$strand, -1L)) "reverse" else "forward") else NA,
                "Description" = ensembl$description
              ))
  tagList(
    div(class = "card", div(class = "card-title", icon("circle-nodes"), "Structure Annotation"),
        bc_primary_gene_detail(d)),
    div(class = "card", div(class = "card-title", icon("dna"), "NCBI Gene"), ncbi_prov, ncbi_body),
    div(class = "card", div(class = "card-title", icon("dna"), "Ensembl (GRCh37/hg19 - matches this app's genome build)"), ens_prov, ens_body),
    div(class = "card", div(class = "card-title", icon("map-location-dot"), "UCSC Genome Browser"),
        if (!is.null(ucsc_link)) tagList(
          p(class = "submodule-desc", "The real UCSC track view at this CpG's exact coordinate - not redrawn locally (spec: link to the original visualization)."),
          tags$a(href = ucsc_link, target = "_blank", rel = "noopener", class = "btn btn-sm btn-default", icon("arrow-up-right-from-square"), " Open in UCSC Genome Browser ")
        ) else div(class = "empty-note", icon("circle-info"), "No genomic coordinate available to link to UCSC."))
  )
}

bc_section_go <- function(d) {
  go <- d$go_terms
  body <- if (is.null(go) || nrow(go) == 0) div(class = "empty-note", icon("circle-info"), "No Gene Ontology (biological process) terms found for this gene, or GO.db is not installed in this deployment.")
          else {
            go2 <- go
            go2$Link <- sprintf('<a href="https://www.ebi.ac.uk/QuickGO/term/%s" target="_blank" rel="noopener">%s</a>', go2$GOID, go2$GOID)
            DT::datatable(go2[, c("Link", "TERM")], colnames = c("GO ID", "Biological process"), rownames = FALSE, escape = 1,
                          options = list(dom = "t", paging = FALSE), class = "stripe hover compact")
          }
  div(class = "card", div(class = "card-title", icon("circle-nodes"), "Gene Ontology (Biological Process)"),
      p(class = "submodule-desc", "Local functional annotation (GO.db/org.Hs.eg.db) for this CpG's primary associated gene, with deep links to the original QuickGO record."),
      body)
}

bc_section_disease_associations <- function(ext_all, gene_symbol, ensembl_id) {
  gen <- ext_all$genetics
  ot_prov <- bc_db_provenance("api.platform.opentargets.org", sprintf("https://platform.opentargets.org/target/%s", ensembl_id %||% ""), "on Open Targets")
  ot_body <- if (is.null(gen)) div(class = "empty-note", icon("circle-info"), "Not yet looked up - click \"Look Up Disease/Genetics Evidence\" below.")
             else if (!isTRUE(gen$ok)) div(class = "empty-note", icon("triangle-exclamation"), gen$reason %||% "Open Targets lookup unavailable.")
             else if (is.null(gen$diseases) || nrow(gen$diseases) == 0) div(class = "empty-note", icon("circle-info"), "No disease associations found in Open Targets for this gene.")
             else DT::datatable(gen$diseases, rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")
  gwas <- ext_all$gwas_catalog
  gwas_prov <- bc_db_provenance("www.ebi.ac.uk/gwas", sprintf("https://www.ebi.ac.uk/gwas/genes/%s", gene_symbol %||% ""), "on GWAS Catalog")
  gwas_body <- if (is.null(gwas)) div(class = "empty-note", icon("circle-info"), "Not yet looked up - click \"Look Up Disease/Genetics Evidence\" below.")
               else if (!isTRUE(gwas$ok)) div(class = "empty-note", icon("triangle-exclamation"), gwas$reason %||% "GWAS Catalog lookup unavailable.")
               else if (is.null(gwas$traits)) div(class = "empty-note", icon("circle-info"), "No GWAS Catalog traits matched this gene symbol.")
               else tagList(
                 p(class = "submodule-desc", sprintf("Trait-name text-index match (%d total matched trait record(s), showing top %d) - a supplementary cross-check, not a curated gene->trait join; Open Targets above is the primary curated disease-association source.", gwas$n_total %||% nrow(gwas$traits), nrow(gwas$traits))),
                 DT::datatable(gwas$traits, rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")
               )
  tagList(
    div(class = "empty-note", icon("circle-info"), "EWAS-derived trait associations (EWAS Catalog / EWAS Atlas) are shown on the \"EWAS Catalog\" / \"EWAS Atlas\" options of the External Databases tab. This subtab covers curated gene-level disease/genetic-association databases."),
    div(class = "card", div(class = "card-title", icon("dna"), "Genetic Association & Tractability (Open Targets)"), ot_prov,
        p(class = "submodule-desc", "Aggregates GWAS Catalog and other genetic evidence via a curated gene-to-disease join."), ot_body),
    div(class = "card", div(class = "card-title", icon("magnifying-glass"), "GWAS Catalog (trait search)"), gwas_prov, gwas_body),
    div(class = "card", div(class = "card-title", icon("ban"), "DisGeNET"),
        div(class = "empty-note", icon("circle-info"), "External database - API requires registration/licensing not configured in this deployment.",
            tags$br(), tags$a(href = "https://www.disgenet.org/", target = "_blank", rel = "noopener", "Search DisGeNET manually")))
  )
}

bc_section_regulatory <- function(ext_all, d) {
  reg <- ext_all$regulatory
  reg_link <- if (!is.na(d$resolved$chr) && !is.na(d$resolved$pos))
    sprintf("https://grch37.ensembl.org/Homo_sapiens/Location/View?r=%s:%s-%s", sub("^chr", "", d$resolved$chr),
            round(d$resolved$pos - 1000), round(d$resolved$pos + 1000)) else NA_character_
  reg_prov <- bc_db_provenance("grch37.rest.ensembl.org", reg_link, "on Ensembl (regulatory track, GRCh37)")
  reg_body <- if (is.null(reg)) div(class = "empty-note", icon("circle-info"), "Not yet looked up - click \"Look Up Regulatory Evidence\" below.")
              else if (!isTRUE(reg$ok)) div(class = "empty-note", icon("triangle-exclamation"), reg$reason %||% "Ensembl regulatory-overlap lookup unavailable.")
              else if (is.null(reg$features)) div(class = "empty-note", icon("circle-info"), sprintf("No Ensembl-annotated regulatory features within the +/-1kb window around this CpG (chr%s:%s-%s).", d$resolved$chr, format(reg$window[1], big.mark = ","), format(reg$window[2], big.mark = ",")))
              else DT::datatable(reg$features, colnames = c("Feature type", "Description", "Start", "End", "Ensembl Regulatory ID"),
                                  rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")
  enc <- ext_all$encode
  enc_prov <- bc_db_provenance("www.encodeproject.org", sprintf("https://www.encodeproject.org/search/?searchTerm=%s&type=Experiment", utils::URLencode(d$primary_gene %||% "", reserved = TRUE)), "on ENCODE")
  enc_body <- if (is.null(enc)) div(class = "empty-note", icon("circle-info"), "Not yet looked up - click \"Look Up Regulatory Evidence\" below.")
              else if (!isTRUE(enc$ok)) div(class = "empty-note", icon("triangle-exclamation"), enc$reason %||% "ENCODE lookup unavailable.")
              else if (is.null(enc$experiments)) div(class = "empty-note", icon("circle-info"), "No ENCODE experiments found for this gene.")
              else {
                df <- enc$experiments
                df$Link <- sprintf('<a href="https://www.encodeproject.org%s" target="_blank" rel="noopener">Open</a>', df[["Link path"]])
                tagList(
                  p(class = "submodule-desc", sprintf("%d matching ENCODE experiment(s) total (showing top %d).", enc$n_total %||% nrow(df), nrow(df))),
                  DT::datatable(df[, c("Accession", "Assay", "Biosample", "Status", "Link")], rownames = FALSE, escape = 5,
                                options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")
                )
              }
  tagList(
    div(class = "card", div(class = "card-title", icon("layer-group"), "Regulatory Region (Ensembl Regulatory Build, by CpG coordinate)"),
        p(class = "submodule-desc", "Regulatory features within 1kb of this CpG."),
        reg_prov, reg_body),
    div(class = "card", div(class = "card-title", icon("chart-simple"), "ENCODE Experiments (gene-level)"),
        p(class = "submodule-desc", "Regulatory/epigenomic assays (ChIP-seq, DNase-seq, ATAC-seq, etc.) associated with this gene."),
        enc_prov, enc_body)
  )
}

bc_section_expression <- function(ext_all, gene, ensembl_id) {
  hpa <- ext_all$hpa
  prov <- bc_db_provenance("www.proteinatlas.org", sprintf("https://www.proteinatlas.org/%s-%s", ensembl_id %||% "", gene %||% ""), "on Human Protein Atlas")
  if (is.null(hpa)) return(div(class = "card", div(class = "card-title", icon("layer-group"), "Baseline Tissue Expression (Human Protein Atlas)"),
                                div(class = "empty-note", icon("circle-info"), "Not yet looked up - click \"Look Up Expression Evidence\" below.")))
  if (!isTRUE(hpa$ok)) return(div(class = "card", div(class = "card-title", icon("layer-group"), "Baseline Tissue Expression (Human Protein Atlas)"), prov,
                                   div(class = "empty-note", icon("triangle-exclamation"), hpa$reason %||% "Human Protein Atlas lookup unavailable.")))
  pairs <- list("Tissue specificity" = hpa$tissue_specificity, "Blood lineage specificity" = hpa$blood_specificity,
                "Blood expression cluster" = hpa$blood_cluster, "Secretome" = hpa$secretome, "Protein class" = hpa$protein_class)
  div(class = "card", div(class = "card-title", icon("layer-group"), "Baseline Tissue Expression (Human Protein Atlas)"), prov,
      p(class = "submodule-desc", "Sanity check: does this gene's baseline expression pattern plausibly match the tissue/cell type your methylation data comes from?"),
      bc_kv_table(pairs),
      if (!is.null(hpa$tissue_top)) tagList(tags$div(style = "margin-top:10px;", tags$b("Top tissues (nTPM)")), DT::datatable(hpa$tissue_top, rownames = FALSE, options = list(dom = "t", paging = FALSE), class = "stripe hover compact")) else NULL,
      if (!is.null(hpa$blood_top)) tagList(tags$div(style = "margin-top:10px;", tags$b("Top blood lineages (nTPM)")), DT::datatable(hpa$blood_top, rownames = FALSE, options = list(dom = "t", paging = FALSE), class = "stripe hover compact")) else NULL,
      div(class = "empty-note", style = "margin-top:8px;", icon("circle-info"), "GTEx is not integrated separately in this deployment - it would duplicate the baseline tissue expression already shown above via Human Protein Atlas.")
  )
}

bc_section_external_datasets <- function(ext_all, query_label) {
  geo <- ext_all$geo; bio <- ext_all$biostudies; enc <- ext_all$encode
  geo_prov <- bc_db_provenance("eutils.ncbi.nlm.nih.gov (db=gds)", sprintf("https://www.ncbi.nlm.nih.gov/gds/?term=%s", utils::URLencode(query_label %||% "", reserved = TRUE)), "on GEO")
  geo_body <- if (is.null(geo)) div(class = "empty-note", icon("circle-info"), "Not yet searched - click \"Search External Datasets\" below.")
              else if (!isTRUE(geo$ok)) div(class = "empty-note", icon("triangle-exclamation"), geo$reason %||% "GEO lookup unavailable.")
              else if (is.null(geo$series)) div(class = "empty-note", icon("circle-info"), "No matching GEO series/datasets found.")
              else {
                df <- geo$series
                df$Link <- sprintf('<a href="%s" target="_blank" rel="noopener">Open</a>', vapply(df$Accession, geo_link, character(1)))
                DT::datatable(df[, c("Accession", "Title", "Data type", "Organism", "Samples", "Date", "Link")], rownames = FALSE, escape = 7,
                              options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")
              }
  bio_prov <- bc_db_provenance("www.ebi.ac.uk/biostudies", sprintf("https://www.ebi.ac.uk/biostudies/studies?search=%s", utils::URLencode(query_label %||% "", reserved = TRUE)), "on BioStudies")
  bio_body <- if (is.null(bio)) div(class = "empty-note", icon("circle-info"), "Not yet searched - click \"Search External Datasets\" below.")
              else if (!isTRUE(bio$ok)) div(class = "empty-note", icon("triangle-exclamation"), bio$reason %||% "BioStudies lookup unavailable.")
              else if (is.null(bio$studies)) div(class = "empty-note", icon("circle-info"), "No matching BioStudies/ArrayExpress records found.")
              else {
                df <- bio$studies
                df$Link <- sprintf('<a href="https://www.ebi.ac.uk/biostudies/studies/%s" target="_blank" rel="noopener">Open</a>', df$Accession)
                tagList(
                  p(class = "submodule-desc", sprintf("%d total matching record(s) (showing top %d).", bio$n_total %||% nrow(df), nrow(df))),
                  DT::datatable(df[, c("Accession", "Title", "Type", "Release date", "Link")], rownames = FALSE, escape = 5,
                                options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")
                )
              }
  enc_body2 <- if (is.null(enc)) div(class = "empty-note", icon("circle-info"), "Not yet searched - the \"Look Up Regulatory Evidence\" button on the Regulatory/Epigenomics subtab queries the same ENCODE index; its result is shown there.")
               else if (!isTRUE(enc$ok)) div(class = "empty-note", icon("triangle-exclamation"), enc$reason %||% "ENCODE lookup unavailable.")
               else if (is.null(enc$experiments)) div(class = "empty-note", icon("circle-info"), "No ENCODE datasets found.")
               else p(class = "submodule-desc", sprintf("%d ENCODE experiment dataset(s) found for this gene - see the Regulatory/Epigenomics subtab for the full table.", enc$n_total %||% nrow(enc$experiments)))
  tagList(
    div(class = "card", div(class = "card-title", icon("database"), "GEO"), geo_body),
    div(class = "card", div(class = "card-title", icon("database"), "BioStudies / ArrayExpress"), bio_prov, bio_body),
    div(class = "card", div(class = "card-title", icon("dna"), "ENCODE Datasets"), enc_body2)
  )
}

bc_section_pubmed_literature <- function(ext_all, identifier) {
  lit <- ext_all$literature; q <- ext_all$literature_query
  prov <- bc_db_provenance("eutils.ncbi.nlm.nih.gov", sprintf("https://pubmed.ncbi.nlm.nih.gov/?term=%s", utils::URLencode(q %||% identifier %||% "", reserved = TRUE)), "on PubMed")
  body <- if (is.null(lit)) div(class = "empty-note", icon("circle-info"), "Not yet searched - pick a query preset below and click \"Search Literature\".")
          else if (!isTRUE(lit$ok)) div(class = "empty-note", icon("triangle-exclamation"), lit$reason %||% "PubMed lookup unavailable.")
          else if (is.null(lit$papers)) div(class = "empty-note", icon("circle-info"), "No matching PubMed records were returned for this query.")
          else {
            df <- lit$papers
            df$Classification <- vapply(seq_len(nrow(df)), function(i) {
              cls <- bc_literature_classify(df$Title[i])
              if (length(cls) == 0) "Unclassified" else paste(cls, collapse = ", ")
            }, character(1))
            df$Link <- sprintf('<a href="https://pubmed.ncbi.nlm.nih.gov/%s/" target="_blank" rel="noopener">%s</a>', df$PMID, df$PMID)
            DT::datatable(df[, c("Title", "Authors", "Journal", "Year", "Classification", "Link")],
                          colnames = c("Title", "Authors", "Journal", "Year", "Classification (automated)", "PMID"),
                          rownames = FALSE, escape = 6, filter = "top", options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
          }
  div(class = "card", div(class = "card-title", icon("book-open"), "PubMed & Literature"), prov,
      if (!is.null(q)) p(class = "submodule-desc", sprintf('Query: "%s"', q)) else NULL,
      if (!is.null(lit) && isTRUE(lit$ok) && !is.null(lit$papers)) div(class = "empty-note", style = "margin-bottom:8px;", icon("circle-info"),
        "\"Classification\" is an automated keyword heuristic over the title only, not a database-provided fact - always verify by reading the abstract.") else NULL,
      body)
}

bc_collect_meta_rows <- function(ext_all) {
  if (is.null(ext_all)) return(NULL)
  metas <- Filter(Negate(is.null), lapply(ext_all, function(x) if (is.list(x)) x$meta else NULL))
  if (length(metas) == 0) return(NULL)
  do.call(rbind, lapply(metas, function(m) {
    data.frame(Source = m$source %||% NA_character_, Query = as.character(m$query %||% NA_character_)[1],
               Endpoint = m$endpoint %||% NA_character_,
               `Retrieved (UTC)` = if (!is.null(m$retrieved_at)) format(m$retrieved_at, tz = "UTC", usetz = TRUE) else NA_character_,
               Status = m$status %||% NA_character_, `Records returned` = m$n_records %||% NA_integer_,
               check.names = FALSE, stringsAsFactors = FALSE)
  }))
}

bc_section_source_info <- function(ext_all) {
  df <- bc_collect_meta_rows(ext_all)
  body <- if (is.null(df) || nrow(df) == 0) div(class = "empty-note", icon("circle-info"), "No external database has been queried yet this session for this biomarker - use the subtabs above to look up evidence.")
          else DT::datatable(df[order(df$Source), , drop = FALSE], rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
  div(class = "card", div(class = "card-title", icon("clipboard-list"), "Source & API Information"),
      p(class = "submodule-desc", "One row per external query actually made this session for this biomarker: source, query submitted, endpoint, retrieval timestamp, response status, and record count. Every number shown elsewhere on this card traces back either to a row here, or to a locally-bundled annotation package listed on the Database Sources card below."),
      body)
}

bc_evidence_status <- function(res, result_field) {
  if (is.null(res)) return("Not yet run")
  if (!isTRUE(res$ok)) return("Failed")
  val <- res
  for (part in strsplit(result_field, "\\$")[[1]]) val <- val[[part]]
  has_rows <- !is.null(val) && (!is.data.frame(val) || nrow(val) > 0)
  if (isTRUE(has_rows)) "Results found" else "No results"
}

BC_EVIDENCE_DBS <- list(
  list(key = "ewascatalog_cpg", label = "EWAS Catalog", field = "df"),
  list(key = "ewasatlas", label = "EWAS Atlas", field = "data$associationList"),
  list(key = "ncbi_gene", label = "NCBI Gene", field = "name"),
  list(key = "ensembl", label = "Ensembl", field = "ensembl_id"),
  list(key = "regulatory", label = "Ensembl Regulatory Build", field = "features"),
  list(key = "kegg", label = "KEGG", field = "pathways"),
  list(key = "reactome", label = "Reactome", field = "pathways"),
  list(key = "wikipathways", label = "WikiPathways", field = "pathways"),
  list(key = "genetics", label = "Open Targets (disease/genetics)", field = "diseases"),
  list(key = "gwas_catalog", label = "GWAS Catalog", field = "traits"),
  list(key = "hpa", label = "Human Protein Atlas", field = "tissue_top"),
  list(key = "encode", label = "ENCODE", field = "experiments"),
  list(key = "geo", label = "GEO", field = "series"),
  list(key = "biostudies", label = "BioStudies / ArrayExpress", field = "studies"),
  list(key = "literature", label = "PubMed literature search", field = "papers")
)

bc_status_chip <- function(label, state) {
  cls <- switch(state, "Results found" = "status-done", "Failed" = "status-pending", "status-neutral")
  ic <- switch(state, "Results found" = "circle-check", "No results" = "circle-minus", "Failed" = "triangle-exclamation", "circle-info")
  span(class = paste("pipeline-status-chip", cls), icon(ic), sprintf("%s: %s", label, state))
}

bc_section_db_comparison <- function(ext_all) {
  if (is.null(ext_all)) return(NULL)
  rows <- list()
  add <- function(resource, evidence_type, results_n, significance, source) {
    rows[[length(rows) + 1]] <<- data.frame(Resource = resource, `Evidence type` = evidence_type, Results = results_n,
                                             Significance = significance, Source = source, check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (isTRUE(ext_all$ewascatalog_cpg$ok)) add("EWAS Catalog", "Epigenome-wide association (this CpG)", nrow(ext_all$ewascatalog_cpg$df), "Reported P-value / effect size (own study, not adjusted here)", "ewascatalog.org")
  if (isTRUE(ext_all$ewasatlas$ok)) add("EWAS Atlas", "Epigenome-wide association (this CpG)", if (!is.null(ext_all$ewasatlas$data$associationList)) nrow(ext_all$ewasatlas$data$associationList) else 0L, "Reported correlation (own study, not adjusted here)", "ngdc.cncb.ac.cn/ewas")
  if (isTRUE(ext_all$kegg$ok)) add("KEGG", "Curated pathway membership", nrow(ext_all$kegg$pathways), "Not applicable (membership, not enrichment)", "KEGGREST")
  if (isTRUE(ext_all$reactome$ok)) add("Reactome", "Curated pathway membership", nrow(ext_all$reactome$pathways), "Not applicable (membership, not enrichment)", "Reactome ContentService")
  if (isTRUE(ext_all$wikipathways$ok)) add("WikiPathways", "Curated pathway membership", nrow(ext_all$wikipathways$pathways), "Not applicable (membership, not enrichment)", "msigdbr")
  if (isTRUE(ext_all$genetics$ok)) add("Open Targets", "Genetic association", ext_all$genetics$n_diseases %||% 0L, "Association score (not a p-value)", "Open Targets Platform")
  if (isTRUE(ext_all$gwas_catalog$ok)) add("GWAS Catalog", "Trait-name text search", ext_all$gwas_catalog$n_total %||% 0L, "Not applicable (search index, not a curated join)", "GWAS Catalog")
  if (isTRUE(ext_all$hpa$ok)) add("Human Protein Atlas", "Baseline tissue/blood expression", if (!is.null(ext_all$hpa$tissue_top)) nrow(ext_all$hpa$tissue_top) else 0L, "Not applicable", "proteinatlas.org")
  if (isTRUE(ext_all$encode$ok)) add("ENCODE", "Regulatory/epigenomic experiments", ext_all$encode$n_total %||% 0L, "Not applicable", "encodeproject.org")
  if (isTRUE(ext_all$geo$ok)) add("GEO", "Independent dataset search", if (!is.null(ext_all$geo$series)) nrow(ext_all$geo$series) else 0L, "Not applicable", "NCBI GEO")
  if (isTRUE(ext_all$biostudies$ok)) add("BioStudies/ArrayExpress", "Independent dataset search", ext_all$biostudies$n_total %||% 0L, "Not applicable", "EBI BioStudies")
  if (isTRUE(ext_all$literature$ok)) add("PubMed", "Literature", if (!is.null(ext_all$literature$papers)) nrow(ext_all$literature$papers) else 0L, "Not applicable", "NCBI E-utilities")
  if (length(rows) == 0) {
    return(div(class = "card", div(class = "card-title", icon("scale-balanced"), "Database Comparison"),
               div(class = "empty-note", icon("circle-info"), "Run at least one database above to populate this comparison.")))
  }
  df <- do.call(rbind, rows)
  div(class = "card", div(class = "card-title", icon("scale-balanced"), "Database Comparison"),
      p(class = "submodule-desc", "Side-by-side view of what each already-run database returned, so you can judge whether independent resources converge on the same biology."),
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact"))
}

bc_section_panel_identity <- function(res) {
  if (is.null(res) || !isTRUE(res$ok) || is.null(res$df)) {
    return(div(class = "card", div(class = "card-title", icon("list-check"), "Identifier Resolution"),
               div(class = "empty-note", icon("triangle-exclamation"), res$reason %||% "No identifiers to resolve.")))
  }
  df <- res$df[, c("input_id", "detected_type", "status_label", "resolved_id", "annotated_gene")]
  colnames(df) <- c("Submitted identifier", "Detected type", "Resolution status", "Resolved ID", "Annotated gene")
  div(class = "card", div(class = "card-title", icon("list-check"), "Identifier Resolution"),
      p(class = "submodule-desc", sprintf("%d identifier(s) submitted, %d resolved, %d unresolved. Every submitted identifier is listed below - none are silently dropped.", res$n_submitted, res$n_resolved, res$n_unresolved)),
      DT::datatable(df, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact"))
}

bc_section_panel_overview <- function(gene_rows, cpg_rows) {
  tagList(
    if (!is.null(gene_rows)) tagList(tags$b(sprintf("Resolved genes (%d)", nrow(gene_rows))),
      DT::datatable(gene_rows, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")) else NULL,
    if (!is.null(cpg_rows)) tagList(tags$div(style = "margin-top:10px;"), tags$b(sprintf("Resolved CpGs (%d)", nrow(cpg_rows))),
      DT::datatable(cpg_rows, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")) else NULL,
    if (is.null(gene_rows) && is.null(cpg_rows)) div(class = "empty-note", icon("circle-info"), "No resolved genes or CpGs to summarize.") else NULL
  )
}

bc_section_panel_convergence <- function(conv, item_label, btn_label) {
  if (is.null(conv)) return(div(class = "card", div(class = "card-title", icon("circle-nodes"), sprintf("%s Convergence", item_label)),
                                 div(class = "empty-note", icon("circle-info"), sprintf("Not yet looked up - click \"%s\" below.", btn_label))))
  if (!isTRUE(conv$ok) || is.null(conv$table)) return(div(class = "card", div(class = "card-title", icon("circle-nodes"), sprintf("%s Convergence", item_label)),
                                    div(class = "empty-note", icon("circle-info"), conv$reason %||% "No convergence found.")))
  div(class = "card", div(class = "card-title", icon("circle-nodes"), sprintf("%s Convergence", item_label)),
      p(class = "submodule-desc", sprintf("Which %s(s) are shared across multiple biomarkers in this list (queried %d biomarker(s), capped to avoid overloading the external service) - a real overlap count, not a synthetic score.", tolower(item_label), conv$n_queried)),
      DT::datatable(conv$table, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact"))
}

bc_section_panel_evidence_summary <- function(disease_conv, pathway_conv) {
  chips <- list(
    bc_status_chip("Disease/Trait convergence", bc_evidence_status(disease_conv, "table")),
    bc_status_chip("Pathway convergence", bc_evidence_status(pathway_conv, "table"))
  )
  div(class = "card", div(class = "card-title", icon("clipboard-check"), "Evidence Summary (Panel)"),
      p(class = "submodule-desc", "Actual retrieval status - a checkmark only appears once real convergence results were returned."),
      div(class = "pipeline-status-strip", chips))
}

bc_build_panel_report_tags <- function(pcd, disease_conv, pathway_conv) {
  tagList(
    tags$h2(sprintf("Biomarker Panel Report (%d identifier(s))", pcd$resolution$n_submitted)),
    bc_section_panel_evidence_summary(disease_conv, pathway_conv),
    bc_section_panel_identity(pcd$resolution),
    div(class = "card", div(class = "card-title", icon("table"), "Panel Overview"), bc_section_panel_overview(pcd$gene_rows, pcd$cpg_rows)),
    bc_section_panel_convergence(disease_conv, if (length(pcd$gene_ids) > 0) "Disease" else "Trait", "Look Up Disease/Trait Convergence"),
    bc_section_panel_convergence(pathway_conv, "Pathway", "Look Up Pathway Convergence"),
    bc_section_sources(NULL),
    tags$p(style = "color:#888; font-size:12px; margin-top:16px;",
           "Per-biomarker detail (genomic context, EWAS evidence, regulatory context, etc.): select an individual biomarker in single-identifier mode for its full evidence dashboard.")
  )
}

bc_report_css <- function() {
  "body{font-family:-apple-system,Helvetica,Arial,sans-serif; max-width:900px; margin:24px auto; color:#222;}
   .card{border:1px solid #ddd; border-radius:10px; padding:14px 18px; margin-bottom:16px;}
   .card-title{font-weight:700; font-size:15px; margin-bottom:8px;}
   .submodule-desc{color:#666; font-size:12.5px;}
   .empty-note{background:#f6f6f6; border-left:3px solid #999; padding:8px 12px; border-radius:4px; font-size:13px;}
   table{border-collapse:collapse; width:100%;} td,th{border:1px solid #eee; padding:4px 8px; font-size:13px; text-align:left;}"
}

bc_build_report_tags <- function(d, cytoband_df, ext = NULL) {
  region_plot <- tryCatch(bc_plot_region(d$resolved$chr, d$resolved$pos, d$island, d$gene_struct, 5000), error = function(e) NULL)
  ideogram_plot <- tryCatch(bc_plot_ideogram(d$resolved$chr, d$resolved$pos, cytoband_df), error = function(e) NULL)
  dist_plot <- if (isTRUE(d$live$ok) && isTRUE(d$live$overall$ok))
    tryCatch(bc_plot_methylation_dist(data.frame(beta = d$live$beta_row, group = d$live$group_vec)), error = function(e) NULL) else NULL
  sex_plot <- if (isTRUE(d$live$ok) && !is.null(d$live$sex_vec))
    tryCatch(bc_plot_methylation_dist(data.frame(beta = d$live$beta_row, group = d$live$group_vec, sex = d$live$sex_vec), facet_sex = TRUE), error = function(e) NULL) else NULL

  img_tag <- function(p) {
    uri <- bc_ggsave_datauri(p)
    if (is.null(uri)) div(class = "empty-note", "Plot unavailable.") else tags$img(src = uri, style = "max-width:100%; height:auto;")
  }
  ens_id <- if (isTRUE(ext$ensembl$ok)) ext$ensembl$ensembl_id else d$gene_struct$ensembl
  identifier <- if (!is.na(d$primary_gene)) d$primary_gene else d$cpg
  res_ids <- bc_resolve_identifiers(if (!is.na(d$primary_gene)) c(d$cpg, d$primary_gene) else d$cpg, d$array_type)
  res_body <- if (isTRUE(res_ids$ok) && !is.null(res_ids$df)) {
    rdf <- res_ids$df[, c("input_id", "detected_type", "status_label")]
    colnames(rdf) <- c("Submitted identifier", "Detected type", "Resolution status")
    DT::datatable(rdf, rownames = FALSE, options = list(dom = "t", paging = FALSE), class = "stripe hover compact")
  } else NULL

  tagList(
    tags$h2(sprintf("Biomarker Card: %s", d$cpg)),
    tags$p(sprintf("Genome build: GRCh37 / hg19. Array: %s.", d$array_type)),
    div(class = "card", div(class = "card-title", "Identifier Resolution"), res_body),
    bc_section_summary(d), bc_section_interpretation(d, ext), bc_section_evidence_summary(d, ext),
    bc_section_genomic_context(d),
    div(class = "card", div(class = "card-title", "Local Region View"), img_tag(region_plot)),
    div(class = "card", div(class = "card-title", "Chromosome View"), img_tag(ideogram_plot)),
    bc_section_genes(d),
    bc_section_gene_genome(d, ext$ncbi_gene, ext$ensembl),
    bc_section_go(d),
    bc_section_disease_evidence(ext), bc_section_ra_evidence(ext), bc_section_disease_comparison(ext), bc_section_tissue_evidence(ext),
    bc_section_disease_associations(ext, d$primary_gene, ens_id),
    bc_section_kegg(ext), bc_section_kegg_ra_pathway(ext), bc_section_reactome(ext), bc_section_wikipathways(ext),
    bc_section_regulatory(ext, d),
    bc_section_expression(ext, d$primary_gene, ens_id),
    bc_section_external_datasets(ext, if (!is.na(d$primary_gene)) sprintf("%s methylation", d$primary_gene) else sprintf("%s methylation", d$cpg)),
    bc_section_publications(ext), bc_section_methbank(ext),
    bc_section_pubmed_literature(ext, identifier),
    bc_section_sex_specific(d, img_tag(sex_plot)),
    bc_section_db_comparison(ext), bc_section_source_info(ext), bc_section_sources(ext),
    tags$p(style = "color:#888; font-size:12px; margin-top:16px;",
           if (is.null(ext)) "No external database evidence was looked up before this report was generated - go back to the card and use each subtab's lookup button first if you want it included."
           else "Not yet available in this build: Methylation<->Expression correlation evidence, and a rendered KEGG/Reactome pathway diagram image (a real membership/gene-list lookup is shown instead).")
  )
}

mod_methyl_biomarkercard_config <- list(
  id = "biomarkercard", title = "Biomarker Card", icon = "id-card", group = "Interpretation",
  description = "Shows the potential methylomics biomarker profile."
)

mod_methyl_biomarkercard_ui <- function(id) {
  ns <- NS(id)
  div(
    class = "tx-menu-wrap",
    tabsetPanel(
      id = ns("bc_subtabs"), type = "tabs",
      tabPanel("Select Biomarker", br(), uiOutput(ns("select_ui"))),
      tabPanel("Biomarker Card", br(), uiOutput(ns("bc_card_ui")))
    )
  )
}

mod_methyl_biomarkercard_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$select_ui <- renderUI({
      tagList(
        div(class = "card",
            div(class = "card-title", icon("magnifying-glass"), "Data source"),
            radioButtons(ns("bc_source"), NULL, inline = TRUE,
                         choiceNames = list(
                           tagList(icon("database"), " Preloaded methylation data"),
                           tagList(icon("upload"), " Your loaded dataset (Dataset tab)")
                         ),
                         choiceValues = list("preloaded", "live"),
                         selected = if (METH_DATA_AVAILABLE) "preloaded" else "live"),
            if (!METH_DATA_AVAILABLE) div(class = "empty-note", icon("triangle-exclamation"), "The preloaded methylomics results folder is not available in this deployment - only your loaded dataset can be used.")
        ),
        div(class = "card",
            div(class = "card-title", icon("list-check"), "Find a biomarker"),
            radioButtons(ns("bc_search_mode"), NULL, inline = TRUE,
                         choices = c("Paste/type a CpG ID" = "cpg", "Search by gene symbol" = "gene", "Browse preloaded DMP results" = "browse",
                                     "Browse diagnostic model panel" = "panel", "Upload a biomarker list" = "upload",
                                     "Gene/CpG list (multiple biomarkers)" = "list")),
            conditionalPanel(condition = sprintf("input['%s'] == 'cpg'", ns("bc_search_mode")),
                              textInput(ns("bc_cpg_input"), "CpG ID", placeholder = "cg12277888")),
            conditionalPanel(condition = sprintf("input['%s'] == 'list'", ns("bc_search_mode")),
                              p(class = "submodule-desc", "Paste a mixed list of gene symbols and/or CpG IDs (one per line, or comma/space-separated) to explore them together as a panel - resolution status is shown for every entry, never silently dropped."),
                              textAreaInput(ns("bc_list_input"), "Gene symbols and/or CpG IDs", rows = 5, placeholder = "BRCA1\nBRCA2\nTP53\ncg00000029\ncg27665659"),
                              fileInput(ns("bc_list_upload_file"), "...or upload a list (.csv, .txt)", accept = c(".csv", ".txt")),
                              actionButton(ns("bc_list_load_btn"), "Load List", icon = icon("play"), class = "btn-sm")),
            conditionalPanel(condition = sprintf("input['%s'] == 'gene'", ns("bc_search_mode")),
                              fluidRow(
                                column(6, textInput(ns("bc_gene_input"), "Gene symbol", placeholder = "TRIM44")),
                                column(6, br(), actionButton(ns("bc_gene_search_btn"), "Search Gene", icon = icon("magnifying-glass"), class = "btn-sm"))
                              ),
                              uiOutput(ns("bc_gene_results_ui"))),
            conditionalPanel(condition = sprintf("input['%s'] == 'browse'", ns("bc_search_mode")),
                              if (!METH_DATA_AVAILABLE) div(class = "empty-note", icon("circle-info"), "Not available - the preloaded results folder is not present in this deployment.")
                              else tagList(
                                fluidRow(
                                  column(6, radioButtons(ns("bc_browse_sex"), "Sex stratum", inline = TRUE, choices = c("All samples" = "all", "Female" = "female", "Male" = "male"), selected = "female")),
                                  column(6, br(), actionButton(ns("bc_browse_load_btn"), "Load Top Results", icon = icon("play"), class = "btn-sm"))
                                ),
                                uiOutput(ns("bc_browse_results_ui"))
                              )),
            conditionalPanel(condition = sprintf("input['%s'] == 'panel'", ns("bc_search_mode")),
                              if (!METH_DATA_AVAILABLE) div(class = "empty-note", icon("circle-info"), "Not available - the preloaded results folder is not present in this deployment.")
                              else tagList(
                                p(class = "submodule-desc", "The Diagnostic Classifier's own \"potential biomarkers\": the script07 majority-vote CpG panel (in_LASSO/in_Boruta/in_SVMRFE, n_votes) - the same table that module's own \"Feature Source\" tab reads."),
                                fluidRow(
                                  column(6, radioButtons(ns("bc_panel_sex"), "Sex stratum", inline = TRUE, choices = c("Female" = "female", "Male" = "male"), selected = "female")),
                                  column(6, br(), actionButton(ns("bc_panel_load_btn"), "Load Panel", icon = icon("play"), class = "btn-sm"))
                                ),
                                uiOutput(ns("bc_panel_results_ui"))
                              )),
            conditionalPanel(condition = sprintf("input['%s'] == 'upload'", ns("bc_search_mode")),
                              p(class = "submodule-desc", "Upload a plain CpG ID list (one per line, or the first column of a CSV/TSV) or a Feature Selection RDS export (that module's own \"Save Model as RDS\" download) - auto-detected by file extension."),
                              fileInput(ns("bc_upload_file"), "Biomarker list (.csv, .txt, or .rds)", accept = c(".csv", ".txt", ".rds")),
                              actionButton(ns("bc_upload_load_btn"), "Load Uploaded List", icon = icon("play"), class = "btn-sm"),
                              uiOutput(ns("bc_upload_results_ui")))
        ),
        div(class = "empty-note", style = "display:flex; align-items:center; justify-content:space-between; gap:12px; flex-wrap:wrap;",
            uiOutput(ns("bc_selection_status_ui"), inline = TRUE),
            actionButton(ns("bc_generate_btn"), "Generate Biomarker Card", icon = icon("id-card"), class = "btn-primary btn-sm"))
      )
    })

    bc_picked_cpg <- reactiveVal(NULL)
    has_card <- reactiveVal(FALSE)
    observeEvent(input$bc_source, { bc_picked_cpg(NULL); has_card(FALSE) }, ignoreInit = TRUE)
    observeEvent(input$bc_search_mode, bc_picked_cpg(NULL), ignoreInit = TRUE)

    output$bc_selection_status_ui <- renderUI({
      mode <- input$bc_search_mode %||% "cpg"
      if (identical(mode, "list")) {
        n <- length(bc_panel_raw_ids())
        if (n > 0) return(tagList(icon("circle-check", style = "color:#0ca30c;"), tags$b(sprintf("%d identifier(s) loaded as a panel", n)), " - click Generate to build the panel report."))
        return(tagList(icon("circle-info"), "Paste or upload a gene/CpG list and click \"Load List\" (or \"Use Full Table as Panel\" from the upload tab) before clicking Generate."))
      }
      cpg <- if (identical(mode, "cpg")) trimws(input$bc_cpg_input %||% "") else bc_picked_cpg()
      if (!is.null(cpg) && nzchar(cpg)) {
        tagList(icon("circle-check", style = "color:#0ca30c;"), tags$b(sprintf("Selected: %s", cpg)), " - click Generate to build the card.")
      } else {
        tagList(icon("circle-info"), "Select a biomarker above (click a table row, or type a CpG ID) before clicking Generate.")
      }
    })

    gene_search_result <- eventReactive(input$bc_gene_search_btn, {
      req(nzchar(trimws(input$bc_gene_input %||% "")))
      array_type <- if (identical(input$bc_source, "preloaded")) "450K" else (dataset$array_type %||% "450K")
      bc_search_gene(trimws(input$bc_gene_input), array_type)
    }, ignoreInit = TRUE)

    output$bc_gene_results_ui <- renderUI({
      r <- gene_search_result(); req(r)
      if (!isTRUE(r$ok)) return(div(class = "empty-note", icon("triangle-exclamation"), r$reason))
      if (nrow(r$df) == 0) return(div(class = "empty-note", icon("circle-info"), "No CpG probes were found annotated to this gene symbol on this array."))
      tagList(
        p(class = "submodule-desc", sprintf("%d CpG probe(s) found. Click a row to select it, then click \"Generate Biomarker Card\".", nrow(r$df))),
        DT::dataTableOutput(ns("bc_gene_table"))
      )
    })
    output$bc_gene_table <- DT::renderDataTable({
      r <- gene_search_result(); req(r); req(isTRUE(r$ok)); req(nrow(r$df) > 0)
      DT::datatable(r$df[, c("cpg", "chr", "pos", "gene_names", "gene_group", "island_relation")],
                    colnames = c("CpG", "Chr", "Position", "Gene(s)", "Region", "Island relation"),
                    rownames = FALSE, selection = "single", options = list(pageLength = 10, scrollX = TRUE))
    })
    outputOptions(output, "bc_gene_table", suspendWhenHidden = FALSE)
    observeEvent(input$bc_gene_table_rows_selected, {
      r <- gene_search_result(); req(r)
      idx <- input$bc_gene_table_rows_selected
      if (length(idx) == 1) bc_picked_cpg(r$df$cpg[idx])
    })

    browse_table <- eventReactive(input$bc_browse_load_btn, {
      validate(need(METH_DATA_AVAILABLE, "The preloaded methylomics results folder is not available in this deployment."))
      df <- bc_dmp_table(input$bc_browse_sex %||% "female")
      validate(need(!is.null(df), "The preloaded DMP results table is not available in this deployment."))
      df <- df[order(df$fdr_bacon), , drop = FALSE]
      utils::head(df, 500)
    }, ignoreInit = TRUE)

    output$bc_browse_results_ui <- renderUI({
      req(input$bc_browse_load_btn)
      tagList(
        p(class = "submodule-desc", "Top 500 CpGs by FDR (SVA/bacon pipeline). Click a row to select it, then click \"Generate Biomarker Card\"."),
        DT::dataTableOutput(ns("bc_browse_table"))
      )
    })
    output$bc_browse_table <- DT::renderDataTable({
      df <- browse_table(); req(df)
      disp <- df[, c("cpg", "dbeta", "p_bacon", "fdr_bacon")]
      dt <- DT::datatable(disp, colnames = c("CpG", "Delta Beta", "P-value", "FDR"), rownames = FALSE, selection = "single",
                           options = list(pageLength = 10, scrollX = TRUE))
      DT::formatSignif(dt, columns = c("dbeta", "p_bacon", "fdr_bacon"), digits = 4)
    })
    outputOptions(output, "bc_browse_table", suspendWhenHidden = FALSE)
    observeEvent(input$bc_browse_table_rows_selected, {
      df <- browse_table(); req(df)
      idx <- input$bc_browse_table_rows_selected
      if (length(idx) == 1) bc_picked_cpg(df$cpg[idx])
    })

    panel_table <- eventReactive(input$bc_panel_load_btn, {
      validate(need(METH_DATA_AVAILABLE, "The preloaded methylomics results folder is not available in this deployment."))
      df <- load_default_diagnostic_ensemble_votes(input$bc_panel_sex %||% "female")
      validate(need(!is.null(df), "The preloaded diagnostic ensemble-vote panel is not available in this deployment."))
      df[order(-df$n_votes), , drop = FALSE]
    }, ignoreInit = TRUE)

    output$bc_panel_results_ui <- renderUI({
      req(input$bc_panel_load_btn)
      tagList(
        p(class = "submodule-desc", sprintf("%d CpG(s) in the majority-vote panel. Click a row to select it, then click \"Generate Biomarker Card\".", nrow(panel_table()))),
        DT::dataTableOutput(ns("bc_panel_table"))
      )
    })
    output$bc_panel_table <- DT::renderDataTable({
      df <- panel_table(); req(df)
      DT::datatable(df, colnames = c("CpG", "in LASSO", "in Boruta", "in SVM-RFE", "Votes"), rownames = FALSE, selection = "single",
                    options = list(pageLength = 10, scrollX = TRUE))
    })
    outputOptions(output, "bc_panel_table", suspendWhenHidden = FALSE)
    observeEvent(input$bc_panel_table_rows_selected, {
      df <- panel_table(); req(df)
      idx <- input$bc_panel_table_rows_selected
      if (length(idx) == 1) bc_picked_cpg(df$cpg[idx])
    })

    upload_table <- eventReactive(input$bc_upload_load_btn, {
      validate(need(!is.null(input$bc_upload_file), "Upload a .csv, .txt, or .rds file first."))
      path <- input$bc_upload_file$datapath; name <- input$bc_upload_file$name
      if (grepl("\\.rds$", name, ignore.case = TRUE)) {
        loaded <- safe_read_rds(path)
        obj <- if (isTRUE(loaded$ok)) loaded$value else NULL
        validate(need(!is.null(obj) && identical(obj$module, "mod_methyl_featureselection"),
                      loaded$error %||% "Upload a Feature Selection RDS export (from that module's own \"Save Model as RDS\" download), or a plain .csv/.txt CpG ID list."))
        ids <- as.character(obj$final_panel$cpg_ids %||% character(0))
        validate(need(length(ids) > 0, "This RDS export's final_panel$cpg_ids is empty."))
        df <- data.frame(cpg = ids, stringsAsFactors = FALSE)
        ann <- obj$final_panel$annotation
        if (!is.null(ann) && "cpg" %in% colnames(ann)) df <- merge(df, ann, by = "cpg", all.x = TRUE)
      } else if (grepl("\\.txt$", name, ignore.case = TRUE)) {
        pl <- methyl_parse_probe_list(path, name)
        validate(need(isTRUE(pl$ok), pl$error %||% "Could not parse the uploaded file."))
        df <- data.frame(cpg = pl$ids, stringsAsFactors = FALSE)
      } else {
        up <- tryCatch(as.data.frame(data.table::fread(path, showProgress = FALSE)), error = function(e) NULL)
        validate(need(!is.null(up) && nrow(up) > 0, "Could not parse the uploaded file as a delimited table (CSV/TSV)."))
        cpg_col <- intersect(c("cpg", "CpG", "probe", "probe_id", "feature", "ID"), colnames(up))[1]
        ids <- if (!is.na(cpg_col)) as.character(up[[cpg_col]]) else as.character(up[[1]])
        df <- data.frame(cpg = unique(ids[nzchar(ids)]), stringsAsFactors = FALSE)
      }
      validate(need(nrow(df) > 0, "No CpG IDs were found in the uploaded file."))
      df
    }, ignoreInit = TRUE)

    output$bc_upload_results_ui <- renderUI({
      req(input$bc_upload_load_btn)
      tagList(
        p(class = "submodule-desc", sprintf("%d CpG(s) loaded from the uploaded file. Click a row to select it and \"Generate Biomarker Card\" for a single-CpG profile, or send the whole table below as a multi-biomarker panel - this is the intended handoff point for an existing candidate-biomarker table (CpG/Gene/Delta-Beta/P/FDR) from an upstream Methylomics analysis (spec: use it as the starting point, never recompute it).", nrow(upload_table()))),
        actionButton(ns("bc_upload_as_panel_btn"), "Use Full Table as Panel", icon = icon("layer-group"), class = "btn-sm"),
        DT::dataTableOutput(ns("bc_upload_table"))
      )
    })
    output$bc_upload_table <- DT::renderDataTable({
      df <- upload_table(); req(df)
      DT::datatable(df, rownames = FALSE, selection = "single", options = list(pageLength = 10, scrollX = TRUE))
    })
    outputOptions(output, "bc_upload_table", suspendWhenHidden = FALSE)
    observeEvent(input$bc_upload_table_rows_selected, {
      df <- upload_table(); req(df)
      idx <- input$bc_upload_table_rows_selected
      if (length(idx) == 1) bc_picked_cpg(df$cpg[idx])
    })

    bc_panel_raw_ids <- reactiveVal(NULL)
    observeEvent(input$bc_list_load_btn, {
      from_text <- bc_split_tokens(input$bc_list_input %||% "")
      from_file <- if (!is.null(input$bc_list_upload_file)) {
        tryCatch({
          up <- as.data.frame(data.table::fread(input$bc_list_upload_file$datapath, header = FALSE, showProgress = FALSE))
          as.character(up[[1]])
        }, error = function(e) character(0))
      } else character(0)
      bc_panel_raw_ids(unique(c(from_text, from_file)))
    }, ignoreInit = TRUE)
    observeEvent(input$bc_upload_as_panel_btn, {
      df <- upload_table(); req(df)
      id_col <- if ("cpg" %in% colnames(df)) "cpg" else colnames(df)[1]
      bc_panel_raw_ids(unique(as.character(df[[id_col]])))
      updateRadioButtons(session, "bc_search_mode", selected = "list")
    }, ignoreInit = TRUE)

    panel_mode_active <- reactiveVal(FALSE)
    observeEvent(input$bc_generate_btn, {
      panel_mode_active(identical(isolate(input$bc_search_mode), "list"))
      has_card(TRUE)
      updateTabsetPanel(session, "bc_subtabs", selected = "Biomarker Card")
    }, ignoreInit = TRUE)

    panel_card_data <- eventReactive(input$bc_generate_btn, {
      req(identical(isolate(input$bc_search_mode), "list"))
      raw <- bc_panel_raw_ids()
      validate(need(!is.null(raw) && length(raw) > 0, "No gene/CpG list loaded - paste or upload a list (or click \"Use Full Table as Panel\" on the upload tab), then click \"Load List\" before Generate."))
      array_type <- if (identical(isolate(input$bc_source), "preloaded")) "450K" else (dataset$array_type %||% "450K")
      res <- bc_resolve_identifiers(raw, array_type)
      validate(need(isTRUE(res$ok) && res$n_resolved > 0, "None of the submitted identifiers could be resolved - check spelling/format (gene symbols, or CpG IDs like cg00000029)."))
      gene_ids <- unique(res$df$resolved_id[res$df$detected_type == "Gene identifier" & res$df$resolved])
      cpg_ids <- unique(res$df$input_id[res$df$detected_type == "CpG probe" & res$df$resolved])
      list(raw = raw, array_type = array_type, resolution = res, gene_ids = gene_ids, cpg_ids = cpg_ids,
           gene_rows = bc_panel_gene_rows(gene_ids), cpg_rows = bc_panel_cpg_rows(cpg_ids, array_type))
    }, ignoreInit = TRUE)

    panel_disease_conv <- reactiveVal(NULL); panel_pathway_conv <- reactiveVal(NULL)
    observeEvent(panel_card_data(), { panel_disease_conv(NULL); panel_pathway_conv(NULL) }, ignoreInit = TRUE)
    observeEvent(input$bc_panel_disease_btn, {
      pcd <- panel_card_data(); req(pcd)
      if (length(pcd$gene_ids) > 0) panel_disease_conv(bc_panel_disease_convergence_genes(pcd$gene_ids))
      else if (length(pcd$cpg_ids) > 0) panel_disease_conv(bc_panel_ewas_trait_convergence_cpgs(pcd$cpg_ids))
      else panel_disease_conv(list(ok = FALSE, reason = "No resolved genes or CpGs to query.", table = NULL, n_queried = 0L))
    }, ignoreInit = TRUE)
    observeEvent(input$bc_panel_pathway_btn, {
      pcd <- panel_card_data(); req(pcd)
      panel_pathway_conv(bc_panel_pathway_convergence_genes(pcd$gene_ids))
    }, ignoreInit = TRUE)

    output$bc_panel_tab_overview <- renderUI({
      pcd <- panel_card_data(); req(pcd)
      tagList(bc_section_panel_identity(pcd$resolution), div(class = "card", div(class = "card-title", icon("table"), "Panel Overview"), bc_section_panel_overview(pcd$gene_rows, pcd$cpg_rows)))
    })
    output$bc_panel_tab_convergence <- renderUI({
      pcd <- panel_card_data(); req(pcd)
      item_label <- if (length(pcd$gene_ids) > 0) "Disease" else "Trait"
      tagList(
        bc_section_panel_evidence_summary(panel_disease_conv(), panel_pathway_conv()),
        bc_section_panel_convergence(panel_disease_conv(), item_label, "Look Up Disease/Trait Convergence"),
        bc_section_panel_convergence(panel_pathway_conv(), "Pathway", "Look Up Pathway Convergence")
      )
    })
    output$bc_panel_tab_source_info <- renderUI({
      div(class = "card", div(class = "card-title", icon("database"), "About Panel Mode"),
          p(class = "submodule-desc", "Panel mode aggregates evidence across multiple biomarkers at once (disease/trait and pathway convergence). For the full per-database evidence dashboard (Gene Study, Gene Ontology, EWAS Catalog, EWAS Atlas, KEGG, Reactome, WikiPathways, Disease/Genetics, Regulatory/Epigenomics, Expression, External Datasets, Literature) on one specific biomarker, switch to a single-identifier search mode and generate its individual Biomarker Card.")
      )
    })

    card_data <- eventReactive(input$bc_generate_btn, {
      cpg <- switch(input$bc_search_mode,
        cpg = trimws(input$bc_cpg_input %||% ""),
        gene = bc_picked_cpg(),
        browse = bc_picked_cpg(),
        panel = bc_picked_cpg(),
        upload = bc_picked_cpg(),
        NULL
      )
      validate(need(!is.null(cpg) && nzchar(cpg),
                    "Select or enter a CpG ID before generating the Biomarker Card - type a CpG ID, or pick one from a gene search / the preloaded results browser / the diagnostic panel / your uploaded list."))
      source_mode <- isolate(input$bc_source)
      array_type <- if (identical(source_mode, "preloaded")) "450K" else (dataset$array_type %||% "450K")
      resolved <- bc_resolve_cpg(cpg, array_type)
      validate(need(isTRUE(resolved$found_in_manifest) || isTRUE(resolved$found_in_champ),
                    sprintf("\"%s\" could not be found in the %s manifest or the ChAMP reference annotation - check the CpG ID.", cpg, array_type)))

      primary_gene <- if (!is.na(resolved$champ_gene) && nzchar(resolved$champ_gene)) resolved$champ_gene
                       else if (!is.na(resolved$gene_names) && nzchar(resolved$gene_names)) trimws(strsplit(resolved$gene_names, ";")[[1]][1])
                       else NA_character_

      island <- bc_island_parse(resolved$island_name)
      gene_struct <- if (!is.na(primary_gene)) bc_gene_structure(primary_gene) else list(ok = FALSE, reason = "No associated gene to resolve.")
      go_terms <- if (isTRUE(gene_struct$ok)) bc_go_terms(gene_struct$entrez) else NULL

      dmp_female <- bc_dmp_lookup(cpg, "female")
      dmp_male <- bc_dmp_lookup(cpg, "male")
      dmr_female <- bc_dmr_lookup(resolved$chr, resolved$pos, "female")
      dmr_male <- bc_dmr_lookup(resolved$chr, resolved$pos, "male")

      live <- bc_dataset_evidence(cpg, dataset)

      single_cpg_diag <- if (isTRUE(live$ok) && isTRUE(live$overall$ok)) {
        bc_single_cpg_roc(live$beta_row, live$group_vec, live$case_label, live$control_label)
      } else list(ok = FALSE, reason = "No case/control split for this CpG (see Dataset tab).")
      single_cpg_cv <- if (isTRUE(live$ok) && isTRUE(live$overall$ok)) {
        bc_single_cpg_cv(live$beta_row, live$group_vec, live$case_label, live$control_label)
      } else list(ok = FALSE, reason = "No case/control split for this CpG (see Dataset tab).")
      panel_membership <- bc_panel_membership_lookup(cpg)
      panel_perprobe <- bc_panel_perprobe_lookup(cpg)

      dataset_snapshot <- list(beta = dataset$beta, source = dataset$source,
                                array_type = array_type, input_scale = dataset$input_scale)

      list(cpg = cpg, source_mode = source_mode, array_type = array_type, resolved = resolved,
           primary_gene = primary_gene, island = island, gene_struct = gene_struct, go_terms = go_terms,
           dmp_female = dmp_female, dmp_male = dmp_male, dmr_female = dmr_female, dmr_male = dmr_male, live = live,
           single_cpg_diag = single_cpg_diag, single_cpg_cv = single_cpg_cv,
           panel_membership = panel_membership, panel_perprobe = panel_perprobe,
           dataset_snapshot = dataset_snapshot)
    }, ignoreInit = TRUE)

    observeEvent(card_data(), {
      d <- card_data()
      if (!is.null(results)) {
        results$biomarkercard <- list(
          cpg = d$cpg, gene = d$primary_gene %||% NA_character_, chr = d$resolved$chr, pos = d$resolved$pos,
          female_fdr = if (!is.null(d$dmp_female)) d$dmp_female$fdr_bacon else NA_real_,
          male_fdr = if (!is.null(d$dmp_male)) d$dmp_male$fdr_bacon else NA_real_
        )
      }
    }, ignoreInit = TRUE)

    ewascatalog_data <- reactiveVal(NULL); ewasatlas_data <- reactiveVal(NULL)
    kegg_data <- reactiveVal(NULL); kegg_ra_data <- reactiveVal(NULL)
    reactome_data <- reactiveVal(NULL); wikipathways_data <- reactiveVal(NULL)
    observeEvent(card_data(), {
      ewascatalog_data(NULL); ewasatlas_data(NULL)
      kegg_data(NULL); kegg_ra_data(NULL); reactome_data(NULL); wikipathways_data(NULL)
    }, ignoreInit = TRUE)

    observeEvent(input$bc_run_ewascatalog, {
      d <- card_data(); req(d)
      ewascatalog_data(bc_ewascatalog_query("cpg", d$cpg))
    }, ignoreInit = TRUE)

    observeEvent(input$bc_run_ewasatlas, {
      d <- card_data(); req(d)
      ewasatlas_data(bc_ewasatlas_probe(d$cpg))
    }, ignoreInit = TRUE)

    observeEvent(input$bc_run_kegg, {
      d <- card_data(); req(d)
      entrez <- d$gene_struct$entrez %||% NA_character_
      kegg_data(bc_kegg_pathways_for_gene(entrez))
      kegg_ra_data(bc_kegg_ra_pathway_check(entrez, d$primary_gene))
    }, ignoreInit = TRUE)

    observeEvent(input$bc_run_reactome, {
      d <- card_data(); req(d)
      reactome_data(if (!is.na(d$primary_gene)) bc_reactome_pathways_for_gene(d$primary_gene) else list(ok = FALSE, pathways = NULL, reason = "No associated gene."))
    }, ignoreInit = TRUE)

    observeEvent(input$bc_run_wikipathways, {
      d <- card_data(); req(d)
      wikipathways_data(bc_wikipathways_pathways_for_gene(d$gene_struct$entrez %||% NA_character_))
    }, ignoreInit = TRUE)

    ewas_combined <- reactive({
      ewc_cpg <- ewascatalog_data(); ewa <- ewasatlas_data()
      ewa_assoc <- NULL
      if (isTRUE(ewa$ok) && !is.null(ewa$data$associationList) && is.data.frame(ewa$data$associationList) && nrow(ewa$data$associationList) > 0) {
        a <- ewa$data$associationList
        ewa_assoc <- data.frame(source = "EWAS Atlas", trait = a$trait, tissue = NA_character_,
                                 effect = a$correlation, p = NA_character_, pmid = as.character(a$pmid), stringsAsFactors = FALSE)
      }
      ewc_assoc <- NULL
      if (isTRUE(ewc_cpg$ok) && nrow(ewc_cpg$df) > 0) {
        dd <- ewc_cpg$df
        ewc_assoc <- data.frame(source = "EWAS Catalog", trait = dd$trait, tissue = dd$tissue,
                                 effect = dd$beta, p = dd$p, pmid = as.character(dd$pmid), stringsAsFactors = FALSE)
      }
      combined <- do.call(rbind, list(ewc_assoc, ewa_assoc))
      if (!is.null(combined)) combined <- combined[!is.na(combined$trait) & nzchar(combined$trait), , drop = FALSE]

      pmids <- if (!is.null(combined)) unique(combined$pmid[!is.na(combined$pmid) & grepl("^[0-9]+$", combined$pmid)]) else character(0)
      pubs <- bc_pubmed_summaries(pmids)

      ra_rows <- if (!is.null(combined)) combined[grepl("rheumatoid|arthritis", combined$trait, ignore.case = TRUE), , drop = FALSE] else NULL

      disease_counts <- NULL
      if (!is.null(combined) && nrow(combined) > 0) {
        cat_counts <- vapply(BC_DISEASE_CATEGORIES, function(pat) sum(grepl(pat, combined$trait, ignore.case = TRUE)), integer(1))
        other_n <- nrow(combined) - sum(cat_counts)
        disease_counts <- data.frame(category = c(names(BC_DISEASE_CATEGORIES), "Other"), n = c(unname(cat_counts), max(other_n, 0)), stringsAsFactors = FALSE)
        disease_counts <- disease_counts[disease_counts$n > 0, , drop = FALSE]
      }

      tissue_counts <- NULL
      if (isTRUE(ewc_cpg$ok) && nrow(ewc_cpg$df) > 0) {
        tt <- ewc_cpg$df$tissue[!is.na(ewc_cpg$df$tissue) & nzchar(ewc_cpg$df$tissue)]
        if (length(tt) > 0) tissue_counts <- as.data.frame(table(tissue = tt), stringsAsFactors = FALSE)
      }

      list(combined_disease = combined, ra_rows = ra_rows, disease_counts = disease_counts,
           tissue_counts = tissue_counts, publications = pubs)
    })

    output$bc_tab_ewascatalog <- renderUI({
      tagList(
        bc_section_ewascatalog(ewascatalog_data()),
        bc_section_disease_evidence(ext_all()), bc_section_ra_evidence(ext_all()),
        bc_section_disease_comparison(ext_all()), bc_section_tissue_evidence(ext_all()),
        bc_section_publications(ext_all()), bc_section_methbank(ext_all())
      )
    })

    output$bc_tab_ewasatlas <- renderUI({
      tagList(
        bc_section_ewasatlas(ewasatlas_data()),
        bc_section_disease_evidence(ext_all()), bc_section_ra_evidence(ext_all()),
        bc_section_disease_comparison(ext_all()), bc_section_tissue_evidence(ext_all()),
        bc_section_publications(ext_all()), bc_section_methbank(ext_all())
      )
    })

    output$bc_tab_kegg <- renderUI({
      tagList(bc_section_kegg(ext_all()), bc_section_kegg_ra_pathway(ext_all()))
    })

    output$bc_tab_reactome <- renderUI({ bc_section_reactome(ext_all()) })

    output$bc_tab_wikipathways <- renderUI({ bc_section_wikipathways(ext_all()) })

    ncbi_gene_data <- reactiveVal(NULL); ensembl_data <- reactiveVal(NULL)
    genetics_data <- reactiveVal(NULL); gwas_catalog_data <- reactiveVal(NULL)
    regulatory_data <- reactiveVal(NULL); encode_data <- reactiveVal(NULL)
    hpa_data <- reactiveVal(NULL)
    geo_data <- reactiveVal(NULL); biostudies_data <- reactiveVal(NULL)
    literature_data <- reactiveVal(NULL); literature_query_used <- reactiveVal(NULL)
    observeEvent(card_data(), {
      ncbi_gene_data(NULL); ensembl_data(NULL); genetics_data(NULL); gwas_catalog_data(NULL)
      regulatory_data(NULL); encode_data(NULL); hpa_data(NULL)
      geo_data(NULL); biostudies_data(NULL); literature_data(NULL); literature_query_used(NULL)
    }, ignoreInit = TRUE)

    observeEvent(input$bc_gene_genome_btn, {
      d <- card_data(); req(d)
      ncbi_gene_data(bc_ncbi_gene_summary(d$gene_struct$entrez %||% NA_character_))
      ensembl_data(bc_ensembl_gene_lookup(d$primary_gene %||% NA_character_))
    }, ignoreInit = TRUE)

    observeEvent(input$bc_disease_genetics_btn, {
      d <- card_data(); req(d)
      ens_id <- if (isTRUE(ensembl_data()$ok)) ensembl_data()$ensembl_id else d$gene_struct$ensembl
      genetics_data(bc_opentargets_evidence_for_gene(ens_id %||% NA_character_))
      gwas_catalog_data(bc_gwas_catalog_by_gene(d$primary_gene %||% NA_character_))
    }, ignoreInit = TRUE)

    observeEvent(input$bc_regulatory_btn, {
      d <- card_data(); req(d)
      regulatory_data(bc_ensembl_regulatory_overlap(d$resolved$chr, d$resolved$pos, d$resolved$pos))
      encode_data(bc_encode_search(d$primary_gene %||% NA_character_))
    }, ignoreInit = TRUE)

    observeEvent(input$bc_expression_btn, {
      d <- card_data(); req(d)
      ens_id <- if (isTRUE(ensembl_data()$ok)) ensembl_data()$ensembl_id else d$gene_struct$ensembl
      hpa_data(bc_hpa_evidence_for_gene(ens_id %||% NA_character_))
    }, ignoreInit = TRUE)

    observeEvent(input$bc_datasets_btn, {
      d <- card_data(); req(d)
      q <- if (!is.na(d$primary_gene)) sprintf("%s methylation", d$primary_gene) else sprintf("%s methylation", d$cpg)
      geo_data(bc_geo_search(q)); biostudies_data(bc_biostudies_search(q))
      if (is.null(encode_data())) encode_data(bc_encode_search(d$primary_gene %||% NA_character_))
    }, ignoreInit = TRUE)

    observeEvent(input$bc_literature_btn, {
      d <- card_data(); req(d)
      identifier <- if (!is.na(d$primary_gene)) d$primary_gene else d$cpg
      q <- bc_literature_query(identifier, input$bc_lit_preset %||% "Gene only")
      literature_query_used(q)
      literature_data(bc_literature_search(q))
    }, ignoreInit = TRUE)

    ext_all <- reactive({
      d <- card_data()
      cmb <- ewas_combined()
      list(
        gene_symbol = if (!is.null(d)) d$primary_gene else NA_character_,
        ncbi_gene = ncbi_gene_data(), ensembl = ensembl_data(),
        ewascatalog_cpg = ewascatalog_data(), ewasatlas = ewasatlas_data(),
        kegg = kegg_data(), kegg_ra = kegg_ra_data(), reactome = reactome_data(), wikipathways = wikipathways_data(),
        combined_disease = cmb$combined_disease, ra_rows = cmb$ra_rows, disease_counts = cmb$disease_counts,
        tissue_counts = cmb$tissue_counts, publications = cmb$publications, methbank_link = bc_methbank_link(),
        genetics = genetics_data(), gwas_catalog = gwas_catalog_data(),
        regulatory = regulatory_data(), encode = encode_data(), hpa = hpa_data(),
        geo = geo_data(), biostudies = biostudies_data(),
        literature = literature_data(), literature_query = literature_query_used()
      )
    })

    output$bc_identifier_resolution_ui <- renderUI({
      d <- card_data(); req(d)
      raw <- if (!is.na(d$primary_gene)) c(d$cpg, d$primary_gene) else d$cpg
      res <- bc_resolve_identifiers(raw, d$array_type)
      if (!isTRUE(res$ok) || is.null(res$df)) return(NULL)
      df <- res$df[, c("input_id", "detected_type", "status_label")]
      colnames(df) <- c("Submitted identifier", "Detected type", "Resolution status")
      div(class = "card", div(class = "card-title", icon("list-check"), "Identifier Resolution"),
          p(class = "submodule-desc", sprintf("%d identifier(s) checked, %d resolved, %d unresolved.", res$n_submitted, res$n_resolved, res$n_unresolved)),
          DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE), class = "stripe hover compact"))
    })

    output$bc_tab_overview <- renderUI({
      d <- card_data(); req(d)
      tagList(
        uiOutput(ns("bc_identifier_resolution_ui")),
        bc_section_summary(d),
        bc_section_evidence_summary(d, ext_all())
      )
    })

    output$bc_tab_gene_genome <- renderUI({
      d <- card_data(); req(d)
      bc_section_gene_genome(d, ncbi_gene_data(), ensembl_data())
    })

    output$bc_tab_disease <- renderUI({
      d <- card_data(); req(d)
      ens_id <- if (isTRUE(ensembl_data()$ok)) ensembl_data()$ensembl_id else d$gene_struct$ensembl
      bc_section_disease_associations(ext_all(), d$primary_gene, ens_id)
    })

    output$bc_tab_go <- renderUI({
      d <- card_data(); req(d)
      bc_section_go(d)
    })

    output$bc_tab_regulatory <- renderUI({
      d <- card_data(); req(d)
      bc_section_regulatory(ext_all(), d)
    })

    output$bc_tab_expression <- renderUI({
      d <- card_data(); req(d)
      ens_id <- if (isTRUE(ensembl_data()$ok)) ensembl_data()$ensembl_id else d$gene_struct$ensembl
      bc_section_expression(ext_all(), d$primary_gene, ens_id)
    })

    output$bc_tab_datasets <- renderUI({
      d <- card_data(); req(d)
      label <- if (!is.na(d$primary_gene)) sprintf("%s methylation", d$primary_gene) else sprintf("%s methylation", d$cpg)
      bc_section_external_datasets(ext_all(), label)
    })

    output$bc_tab_literature <- renderUI({
      d <- card_data(); req(d)
      identifier <- if (!is.na(d$primary_gene)) d$primary_gene else d$cpg
      bc_section_pubmed_literature(ext_all(), identifier)
    })

    output$bc_tab_source_info <- renderUI({
      tagList(bc_section_source_info(ext_all()), bc_section_db_comparison(ext_all()), bc_section_sources(ext_all()))
    })

    output$bc_tab_biomarker_status <- renderUI({
      d <- card_data(); req(d)
      tagList(
        bc_section_summary(d),
        bc_section_interpretation(d, ext_all()),
        bc_section_evidence_summary(d, ext_all())
      )
    })

    output$bc_card_ui <- renderUI({
      if (!has_card()) return(div(class = "empty-note", icon("circle-info"), "Select a biomarker and click \"Generate Biomarker Card\" on the \"Select Biomarker\" tab."))
      if (isTRUE(panel_mode_active())) {
        pcd <- panel_card_data(); req(pcd)
        return(tagList(
          tabsetPanel(id = ns("bc_panel_subtabs"), type = "tabs",
            tabPanel("Panel Overview", br(), uiOutput(ns("bc_panel_tab_overview"))),
            tabPanel("Disease/Trait & Pathway Convergence", br(),
              fluidRow(
                column(6, actionButton(ns("bc_panel_disease_btn"), "Look Up Disease/Trait Convergence", icon = icon("magnifying-glass-location"), class = "btn-primary btn-sm")),
                column(6, actionButton(ns("bc_panel_pathway_btn"), "Look Up Pathway Convergence", icon = icon("magnifying-glass-location"), class = "btn-primary btn-sm"))
              ),
              withSpinner(uiOutput(ns("bc_panel_tab_convergence")), color = "#2563EB", type = 6)),
            tabPanel("Source & API Information", br(), uiOutput(ns("bc_panel_tab_source_info")))
          ),
          div(class = "card", div(class = "card-title", icon("download"), "Download"),
              p(class = "submodule-desc", "A self-contained HTML panel report (identifier resolution + panel overview + whichever convergence lookups you've already run)."),
              downloadButton(ns("bc_panel_download_report"), "Download Panel Report (HTML)", class = "btn-primary btn-sm"))
        ))
      }
      d <- card_data(); req(d)
      tagList(
        tabsetPanel(id = ns("bc_evidence_subtabs"), type = "tabs",
          tabPanel("CpG description", br(),
            bc_section_genomic_context(d),
            div(class = "card",
                div(class = "card-title", icon("magnifying-glass-chart"), "Local Region View"),
                p(class = "submodule-desc", "Source: Illumina manifest annotation (hg19) + TxDb.Hsapiens.UCSC.hg19.knownGene gene structure, ChAMPdata CpG island context."),
                fluidRow(
                  column(4, numericInput(ns("bc_flank_kb"), "Window (kb, +/- around the CpG)", value = 5, min = 1, max = 200, step = 1)),
                  column(4, br(), actionButton(ns("bc_zoom_tight_btn"), "Zoom to Biomarker", icon = icon("magnifying-glass-plus"), class = "btn-sm")),
                  column(4, br(), actionButton(ns("bc_zoom_local_btn"), "View Local Region", icon = icon("expand"), class = "btn-sm"))
                ),
                withSpinner(plotly::plotlyOutput(ns("bc_region_plot"), height = "260px"), color = "#2563EB", type = 6)
            ),
            div(class = "card",
                div(class = "card-title", icon("dna"), "Chromosome View"),
                p(class = "submodule-desc", sprintf("Source: UCSC cytoBandIdeo (hg19), chromosome %s.", d$resolved$chr %||% "?")),
                withSpinner(plotly::plotlyOutput(ns("bc_ideogram_plot"), height = "220px"), color = "#2563EB", type = 6)
            ),
            bc_section_genes(d)
          ),
          tabPanel("Biomarker Status", br(),
            withSpinner(uiOutput(ns("bc_tab_biomarker_status")), color = "#2563EB", type = 6)
          ),
          tabPanel("Dataset", br(), bc_section_dataset_cohort(d$dataset_snapshot, d$live)),
          tabPanel("Differential Methylation (DMP/DMR)", br(),
            bc_section_sex_specific(d, withSpinner(plotly::plotlyOutput(ns("bc_sex_dist_plot"), height = "320px"), color = "#2563EB", type = 6))
          ),
          tabPanel("Single-Gene vs Multi-Gene Signature", br(),
            bc_section_signature_comparison(d, d$single_cpg_diag,
              if (isTRUE(d$single_cpg_diag$ok)) withSpinner(plotOutput(ns("bc_sg_train_roc_plot"), height = "320px"), color = "#2563EB", type = 6) else NULL)
          ),
          tabPanel("Biomarker Performance", br(),
            bc_section_biomarker_performance(d, d$single_cpg_diag, d$single_cpg_cv,
              train_roc_widget = if (isTRUE(d$single_cpg_diag$ok)) withSpinner(plotOutput(ns("bc_sg_train_roc_plot2"), height = "320px"), color = "#2563EB", type = 6) else NULL,
              train_pr_widget = if (isTRUE(d$single_cpg_diag$ok)) withSpinner(plotOutput(ns("bc_sg_train_pr_plot"), height = "320px"), color = "#2563EB", type = 6) else NULL,
              internal_roc_widget = if (isTRUE(d$single_cpg_cv$ok)) withSpinner(plotOutput(ns("bc_sg_internal_roc_plot"), height = "320px"), color = "#2563EB", type = 6) else NULL)
          ),
          tabPanel("External Databases", br(),
            div(class = "card",
                div(class = "card-title", icon("globe"), "Deep Dive: External Databases"),
                p(class = "submodule-desc", "Pick a database, set your own query parameters, and run it - each is independent and links straight to the real record so you can verify it yourself."),
                radioButtons(ns("bc_db_choice"), NULL, inline = TRUE, choices = c(
                  "Gene Study (NCBI Gene, Ensembl)" = "gene_genome",
                  "Gene Ontology (GO.db, QuickGO)" = "go",
                  "EWAS Catalog" = "ewascatalog",
                  "EWAS Atlas" = "ewasatlas",
                  "KEGG" = "kegg",
                  "Reactome" = "reactome",
                  "WikiPathways" = "wikipathways",
                  "Disease / Genetics (Open Targets, GWAS Catalog)" = "disease",
                  "Regulatory / Epigenomics (Ensembl Regulatory Build, ENCODE)" = "regulatory",
                  "Expression (Human Protein Atlas)" = "expression",
                  "GEO, BioStudies/ArrayExpress" = "datasets",
                  "Literature (PubMed)" = "literature"
                )),
                uiOutput(ns("bc_db_controls_ui"))
            ),
            withSpinner(uiOutput(ns("bc_db_result_ui")), color = "#2563EB", type = 6),
            withSpinner(uiOutput(ns("bc_db_comparison_ui")), color = "#2563EB", type = 6)
          ),
          tabPanel("Download", br(),
            div(class = "card", div(class = "card-title", icon("download"), "Download"),
                p(class = "submodule-desc", "The downloaded report includes every section shown across all tabs above, for this CpG."),
                downloadButton(ns("bc_download_report"), "Download Biomarker Report (HTML)", class = "btn-primary btn-sm"))
          )
        )
      )
    })

    output$bc_db_controls_ui <- renderUI({
      choice <- input$bc_db_choice %||% "gene_genome"
      switch(choice,
        gene_genome = actionButton(ns("bc_gene_genome_btn"), "Look Up Gene Study Evidence", icon = icon("magnifying-glass-location"), class = "btn-primary btn-sm"),
        go = div(class = "empty-note", icon("circle-info"), "Local annotation (GO.db/org.Hs.eg.db) - no query needed, result is shown below."),
        ewascatalog = actionButton(ns("bc_run_ewascatalog"), "Run EWAS Catalog Query", icon = icon("play"), class = "btn-primary btn-sm"),
        ewasatlas = actionButton(ns("bc_run_ewasatlas"), "Run EWAS Atlas Query", icon = icon("play"), class = "btn-primary btn-sm"),
        kegg = actionButton(ns("bc_run_kegg"), "Run KEGG Query", icon = icon("play"), class = "btn-primary btn-sm"),
        reactome = actionButton(ns("bc_run_reactome"), "Run Reactome Query", icon = icon("play"), class = "btn-primary btn-sm"),
        wikipathways = actionButton(ns("bc_run_wikipathways"), "Run WikiPathways Query", icon = icon("play"), class = "btn-primary btn-sm"),
        disease = actionButton(ns("bc_disease_genetics_btn"), "Look Up Disease/Genetics Evidence", icon = icon("magnifying-glass-location"), class = "btn-primary btn-sm"),
        regulatory = actionButton(ns("bc_regulatory_btn"), "Look Up Regulatory Evidence", icon = icon("magnifying-glass-location"), class = "btn-primary btn-sm"),
        expression = actionButton(ns("bc_expression_btn"), "Look Up Expression Evidence", icon = icon("magnifying-glass-location"), class = "btn-primary btn-sm"),
        datasets = actionButton(ns("bc_datasets_btn"), "Search External Datasets", icon = icon("magnifying-glass-location"), class = "btn-primary btn-sm"),
        literature = tagList(
          fluidRow(
            column(6, selectInput(ns("bc_lit_preset"), "Query preset", choices = names(BC_LITERATURE_PRESETS))),
            column(6, br(), actionButton(ns("bc_literature_btn"), "Search Literature", icon = icon("magnifying-glass"), class = "btn-primary btn-sm"))
          )
        )
      )
    })

    output$bc_db_result_ui <- renderUI({
      choice <- input$bc_db_choice %||% "gene_genome"
      switch(choice,
        gene_genome = uiOutput(ns("bc_tab_gene_genome")),
        go = uiOutput(ns("bc_tab_go")),
        ewascatalog = uiOutput(ns("bc_tab_ewascatalog")),
        ewasatlas = uiOutput(ns("bc_tab_ewasatlas")),
        kegg = uiOutput(ns("bc_tab_kegg")),
        reactome = uiOutput(ns("bc_tab_reactome")),
        wikipathways = uiOutput(ns("bc_tab_wikipathways")),
        disease = uiOutput(ns("bc_tab_disease")),
        regulatory = uiOutput(ns("bc_tab_regulatory")),
        expression = uiOutput(ns("bc_tab_expression")),
        datasets = uiOutput(ns("bc_tab_datasets")),
        literature = uiOutput(ns("bc_tab_literature"))
      )
    })

    output$bc_db_comparison_ui <- renderUI({
      bc_section_db_comparison(ext_all())
    })

    observeEvent(input$bc_zoom_tight_btn, updateNumericInput(session, "bc_flank_kb", value = 1))
    observeEvent(input$bc_zoom_local_btn, updateNumericInput(session, "bc_flank_kb", value = 10))

    output$bc_region_plot <- plotly::renderPlotly({
      d <- card_data(); req(d)
      flank_bp <- (input$bc_flank_kb %||% 5) * 1000
      p <- bc_plot_region(d$resolved$chr, d$resolved$pos, d$island, d$gene_struct, flank_bp)
      plotly::ggplotly(p) %>% plotly::layout(showlegend = FALSE)
    })

    output$bc_ideogram_plot <- plotly::renderPlotly({
      d <- card_data(); req(d)
      cb <- bc_cytoband_hg19()
      p <- bc_plot_ideogram(d$resolved$chr, d$resolved$pos, cb)
      plotly::ggplotly(p) %>% plotly::layout(showlegend = FALSE)
    })

    output$bc_dist_plot <- plotly::renderPlotly({
      d <- card_data(); req(d)
      validate(need(isTRUE(d$live$ok) && isTRUE(d$live$overall$ok), "No dataset evidence is available to plot for this CpG."))
      df <- data.frame(beta = d$live$beta_row, group = d$live$group_vec, stringsAsFactors = FALSE)
      plotly::ggplotly(bc_plot_methylation_dist(df))
    })

    output$bc_sex_dist_plot <- plotly::renderPlotly({
      d <- card_data(); req(d)
      validate(need(isTRUE(d$live$ok) && !is.null(d$live$sex_vec), "No sex information is available in the loaded sample sheet to compare female vs male methylation."))
      df <- data.frame(beta = d$live$beta_row, group = d$live$group_vec, sex = d$live$sex_vec, stringsAsFactors = FALSE)
      plotly::ggplotly(bc_plot_methylation_dist(df, facet_sex = TRUE))
    })

    output$bc_sg_train_roc_plot <- renderPlot({
      d <- card_data(); req(d); sgd <- d$single_cpg_diag
      validate(need(isTRUE(sgd$ok), sgd$reason %||% "Single-CpG ROC unavailable."))
      bc_plot_roc(sgd$roc_obj, "Single-CpG ROC (training)", sgd$auc, sgd$ci_lo, sgd$ci_hi)
    })
    output$bc_sg_train_roc_plot2 <- renderPlot({
      d <- card_data(); req(d); sgd <- d$single_cpg_diag
      validate(need(isTRUE(sgd$ok), sgd$reason %||% "Single-CpG ROC unavailable."))
      bc_plot_roc(sgd$roc_obj, "Single-CpG ROC (training)", sgd$auc, sgd$ci_lo, sgd$ci_hi)
    })
    output$bc_sg_train_pr_plot <- renderPlot({
      d <- card_data(); req(d); sgd <- d$single_cpg_diag
      validate(need(isTRUE(sgd$ok), sgd$reason %||% "Single-CpG ROC unavailable."))
      pr <- bc_pr_from_roc(sgd$roc_obj)
      validate(need(isTRUE(pr$ok), pr$reason %||% "Precision-recall curve unavailable."))
      bc_plot_pr(pr, "Single-CpG Precision-Recall (training)")
    })
    output$bc_sg_internal_roc_plot <- renderPlot({
      d <- card_data(); req(d); sgcv <- d$single_cpg_cv
      validate(need(isTRUE(sgcv$ok), sgcv$reason %||% "Cross-validated ROC unavailable."))
      bc_plot_roc(sgcv$roc_obj, "Single-CpG ROC (internal CV, pooled out-of-fold)", sgcv$auc)
    })

    output$bc_download_report <- downloadHandler(
      filename = function() sprintf("biomarker_card_%s.html", card_data()$cpg),
      content = function(file) {
        d <- card_data(); req(d)
        cb <- bc_cytoband_hg19()
        body <- bc_build_report_tags(d, cb, ext_all())
        page <- tags$html(
          tags$head(tags$meta(charset = "utf-8"), tags$title(sprintf("Biomarker Card %s", d$cpg)), tags$style(bc_report_css())),
          tags$body(body)
        )
        htmltools::save_html(page, file = file)
      }
    )

    output$bc_panel_download_report <- downloadHandler(
      filename = function() sprintf("biomarker_panel_report_%d_ids.html", panel_card_data()$resolution$n_submitted),
      content = function(file) {
        pcd <- panel_card_data(); req(pcd)
        body <- bc_build_panel_report_tags(pcd, panel_disease_conv(), panel_pathway_conv())
        page <- tags$html(
          tags$head(tags$meta(charset = "utf-8"), tags$title("Biomarker Panel Report"), tags$style(bc_report_css())),
          tags$body(body)
        )
        htmltools::save_html(page, file = file)
      }
    )
  })
}
