## Biomarker Card submodule: single-gene transcriptomic profile combining
## live dataset evidence, saved DGE/candidate/signature/classifier results,
## and opt-in pathway/external-database lookups into one report.

## ---- Gene identity (org.Hs.eg.db - already library()'d in global.R) -------

.tbc_identity_cache <- new.env(parent = emptyenv())

tbc_gene_identity <- function(symbol) {
  if (is.null(symbol) || is.na(symbol) || !nzchar(symbol)) return(list(ok = FALSE, reason = "No gene symbol given."))
  cached <- .tbc_identity_cache[[symbol]]
  if (!is.null(cached)) return(cached)
  res <- tryCatch({
    map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = symbol, keytype = "SYMBOL",
                                                    columns = c("ENTREZID", "ENSEMBL", "GENENAME", "CHR", "MAP", "GENETYPE")))
    map <- map[!is.na(map$ENTREZID), , drop = FALSE]
    if (nrow(map) == 0) return(list(ok = FALSE, reason = sprintf("No NCBI Entrez Gene entry found for symbol \"%s\".", symbol)))
    alias <- tryCatch({
      a <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = symbol, keytype = "SYMBOL", columns = "ALIAS"))
      unique(stats::na.omit(a$ALIAS[a$ALIAS != symbol]))
    }, error = function(e) character(0))
    list(ok = TRUE, symbol = symbol, entrez = map$ENTREZID[1], ensembl = map$ENSEMBL[1], genename = map$GENENAME[1],
         chr = map$CHR[1], map_location = map$MAP[1], genetype = map$GENETYPE[1], aliases = alias)
  }, error = function(e) e)
  if (inherits(res, "error")) res <- list(ok = FALSE, reason = sprintf("Could not resolve gene identity: %s", conditionMessage(res)))
  .tbc_identity_cache[[symbol]] <- res
  res
}

## ---- NCBI Gene summary (basic description) + UniProt protein name --------
## Both are own small, cached, fail-soft clients on the same list(ok, reason)
## contract as every other external lookup in this file - identity fields
## the user can see are "Not available" rather than missing silently.
.tbc_ncbi_summary_cache <- new.env(parent = emptyenv())
tbc_ncbi_gene_summary <- function(entrez) {
  if (is.null(entrez) || is.na(entrez) || !nzchar(entrez)) return(list(ok = FALSE, reason = "No NCBI Entrez Gene ID available."))
  cached <- .tbc_ncbi_summary_cache[[entrez]]
  if (!is.null(cached)) return(cached)
  res <- if (!requireNamespace("httr2", quietly = TRUE)) list(ok = FALSE, reason = "httr2 is not installed in this deployment.") else tryCatch({
    resp <- httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi") %>%
      httr2::req_url_query(db = "gene", id = entrez, retmode = "json") %>%
      httr2::req_url_query(tool = "arthomix-biomarkercard") %>%
      httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
    if (httr2::resp_status(resp) != 200) return(list(ok = FALSE, reason = "NCBI Gene lookup failed (network error or timeout)."))
    rec <- httr2::resp_body_json(resp)$result[[as.character(entrez)]]
    if (is.null(rec)) return(list(ok = FALSE, reason = sprintf("No NCBI Gene summary found for Entrez ID \"%s\".", entrez)))
    list(ok = TRUE, description = rec$summary %||% NA_character_,
         other_designations = rec$otherdesignations %||% NA_character_,
         other_aliases = rec$otheraliases %||% NA_character_)
  }, error = function(e) list(ok = FALSE, reason = sprintf("NCBI Gene lookup failed: %s", conditionMessage(e))))
  .tbc_ncbi_summary_cache[[entrez]] <- res
  res
}

.tbc_uniprot_cache <- new.env(parent = emptyenv())
tbc_uniprot_protein_name <- function(symbol) {
  if (is.null(symbol) || is.na(symbol) || !nzchar(symbol)) return(list(ok = FALSE, reason = "No gene symbol available for UniProt lookup."))
  cached <- .tbc_uniprot_cache[[symbol]]
  if (!is.null(cached)) return(cached)
  res <- if (!requireNamespace("httr2", quietly = TRUE)) list(ok = FALSE, reason = "httr2 is not installed in this deployment.") else tryCatch({
    uniprot_ids <- {
      m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = symbol, keytype = "SYMBOL", columns = "UNIPROT"))
      unique(stats::na.omit(m$UNIPROT))
    }
    if (length(uniprot_ids) == 0) return(list(ok = FALSE, reason = sprintf("No UniProt ID found for gene symbol \"%s\".", symbol)))
    for (uid in uniprot_ids) {
      resp <- httr2::request(sprintf("https://rest.uniprot.org/uniprotkb/%s.json", uid)) %>%
        httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
      if (httr2::resp_status(resp) != 200) next
      body <- httr2::resp_body_json(resp)
      nm <- body$proteinDescription$recommendedName$fullName$value
      if (!is.null(nm)) return(list(ok = TRUE, name = nm, uniprot = uid))
    }
    list(ok = FALSE, reason = sprintf("UniProt returned no protein name for \"%s\" (UniProt ID(s) tried: %s).", symbol, paste(uniprot_ids, collapse = ", ")))
  }, error = function(e) list(ok = FALSE, reason = sprintf("UniProt lookup failed: %s", conditionMessage(e))))
  .tbc_uniprot_cache[[symbol]] <- res
  res
}

## Direct links to authoritative gene resources - pure, no network.
tbc_gene_identity_links <- function(symbol, entrez = NULL, ensembl = NULL) {
  links <- list("NCBI Gene" = if (!is.null(entrez) && !is.na(entrez)) sprintf("https://www.ncbi.nlm.nih.gov/gene/%s", entrez) else NA_character_)
  links[["Ensembl"]] <- if (!is.null(ensembl) && !is.na(ensembl)) sprintf("https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=%s", ensembl) else NA_character_
  links[["GeneCards"]] <- if (!is.null(symbol) && nzchar(symbol)) sprintf("https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s", symbol) else NA_character_
  links[["UniProt"]] <- if (!is.null(symbol) && nzchar(symbol)) sprintf("https://www.uniprot.org/uniprotkb?query=gene:%s+AND+organism_id:9606", symbol) else NA_character_
  links
}

## Top biological-process GO terms for a gene, with GOID kept for QuickGO links.
tbc_go_terms <- function(entrez, n = 8) {
  if (is.null(entrez) || is.na(entrez) || !nzchar(entrez)) return(NULL)
  if (!requireNamespace("GO.db", quietly = TRUE)) return(NULL)
  res <- tryCatch({
    go <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = entrez, keytype = "ENTREZID", columns = c("GOALL", "ONTOLOGYALL")))
    go <- go[!is.na(go$GOALL) & go$ONTOLOGYALL == "BP", , drop = FALSE]
    ids <- unique(go$GOALL)
    if (length(ids) == 0) return(NULL)
    terms <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = ids, keytype = "GOID", columns = "TERM"))
    terms <- unique(stats::na.omit(terms[, c("GOID", "TERM")]))
    utils::head(terms, n)
  }, error = function(e) NULL)
  if (inherits(res, "error")) return(NULL)
  res
}

## Process-wide cache of keggList()'s pathway-name table (own copy, not shared with methylation card).
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

## Maps SYMBOL -> UNIPROT (1:many) and tries each ID until one returns Reactome pathways.
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

## ---- WikiPathways: pathway membership via the same msigdbr-cached ------
## term<->gene table the Multi-Omics Pathways module already loads
## (mp_get_wikipathways_termgene(), R/multiomics/multiomics_pathway_helpers.R) -
## no separate client, no separate license. Same "which pathways contain this
## gene" contract as tbc_kegg_pathways_for_gene / tbc_reactome_pathways_for_gene above.
tbc_wikipathways_pathways_for_gene <- function(entrez) {
  if (is.null(entrez) || is.na(entrez) || !nzchar(entrez)) return(list(ok = FALSE, pathways = NULL, reason = "No NCBI Gene ID available for WikiPathways lookup."))
  t2g <- tryCatch(mp_get_wikipathways_termgene(), error = function(e) NULL)
  if (is.null(t2g)) return(list(ok = FALSE, pathways = NULL, reason = "WikiPathways gene sets (msigdbr) are not available in this deployment."))
  hit <- t2g$TERM2GENE[!is.na(t2g$TERM2GENE$ncbi_gene) & as.character(t2g$TERM2GENE$ncbi_gene) == as.character(entrez), , drop = FALSE]
  if (nrow(hit) == 0) return(list(ok = TRUE, pathways = data.frame(id = character(0), name = character(0)), reason = NULL))
  ids <- unique(hit$gs_exact_source)
  nm <- t2g$TERM2NAME$gs_name[match(ids, t2g$TERM2NAME$gs_exact_source)]
  list(ok = TRUE, pathways = data.frame(id = ids, name = nm, stringsAsFactors = FALSE), reason = NULL)
}

## ---- Open Targets: genetic association + druggability tractability ----
TBC_OT_TRACTABILITY_TIERS <- c("Approved Drug", "Advanced Clinical", "Phase 1 Clinical")
TBC_OT_MODALITY_LABELS <- c(SM = "Small molecule", AB = "Antibody / biologic", PR = "PROTAC-type degrader", OC = "Other clinical modality")

tbc_opentargets_evidence_for_gene <- function(ensembl, top_n_diseases = 8) {
  if (is.null(ensembl) || is.na(ensembl) || !nzchar(ensembl)) return(list(ok = FALSE, reason = "No Ensembl Gene ID available for Open Targets lookup."))
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, reason = "httr2 is not installed in this deployment."))
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
  if (is.null(tgt)) return(list(ok = FALSE, reason = sprintf("No Open Targets entry found for Ensembl ID \"%s\".", ensembl)))

  tract <- tgt$tractability
  modality_summary <- stats::setNames(lapply(names(TBC_OT_MODALITY_LABELS), function(mod) {
    rows <- Filter(function(r) identical(r$modality, mod), tract)
    approved <- any(vapply(rows, function(r) identical(r$label, "Approved Drug") && isTRUE(r$value), logical(1)))
    tier <- Find(function(t) any(vapply(rows, function(r) identical(r$label, t) && isTRUE(r$value), logical(1))), TBC_OT_TRACTABILITY_TIERS)
    feasible <- any(vapply(rows, function(r) isTRUE(r$value), logical(1)))
    if (approved) "Approved drug exists" else if (!is.null(tier)) tier else if (feasible) "Structurally feasible" else "No tractability evidence"
  }), names(TBC_OT_MODALITY_LABELS))

  ## Per-disease datatype-score columns, built dynamically from whatever datatype ids are present.
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

  list(ok = TRUE, n_diseases = tgt$associatedDiseases$count %||% 0, diseases = diseases, tractability = modality_summary)
}

## Human Protein Atlas baseline tissue/blood-lineage RNA expression, keyed by Ensembl ID.
tbc_hpa_evidence_for_gene <- function(ensembl) {
  if (is.null(ensembl) || is.na(ensembl) || !nzchar(ensembl)) return(list(ok = FALSE, reason = "No Ensembl Gene ID available for Human Protein Atlas lookup."))
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, reason = "httr2 is not installed in this deployment."))
  res <- tryCatch({
    resp <- httr2::request(sprintf("https://www.proteinatlas.org/%s.json", ensembl)) %>%
      httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
    if (httr2::resp_status(resp) != 200) NULL else httr2::resp_body_json(resp, simplifyVector = FALSE)
  }, error = function(e) NULL)
  if (is.null(res)) return(list(ok = FALSE, reason = sprintf("No Human Protein Atlas entry found for Ensembl ID \"%s\".", ensembl)))

  named_list_to_df <- function(x, value_label) {
    if (is.null(x) || length(x) == 0) return(NULL)
    data.frame(Tissue = names(x), Value = vapply(x, function(v) as.character(v %||% NA), character(1)), check.names = FALSE, stringsAsFactors = FALSE) %>%
      stats::setNames(c("Tissue / cell type", value_label))
  }
  list(ok = TRUE,
       tissue_specificity = res[["RNA tissue specificity"]] %||% "Not available",
       tissue_top = named_list_to_df(res[["RNA tissue specific nTPM"]], "nTPM"),
       blood_specificity = res[["RNA blood lineage specificity"]] %||% "Not available",
       blood_top = named_list_to_df(res[["RNA blood lineage specific nTPM"]], "nTPM"),
       blood_cluster = res[["Blood expression cluster"]] %||% NA_character_,
       secretome = res[["Secretome location"]] %||% NA_character_,
       protein_class = if (!is.null(res[["Protein class"]])) paste(unlist(res[["Protein class"]]), collapse = ", ") else NA_character_)
}

## STRING PPI neighborhood; a 404 means "no partners found", not a request failure.
tbc_string_ppi_for_gene <- function(symbol, top_n = 10, required_score = 400) {
  if (is.null(symbol) || is.na(symbol) || !nzchar(symbol)) return(list(ok = FALSE, reason = "No gene symbol available for STRING lookup."))
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, reason = "httr2 is not installed in this deployment."))
  out <- tryCatch({
    resp <- httr2::request("https://string-db.org/api/json/interaction_partners") %>%
      httr2::req_url_query(identifiers = symbol, species = 9606, limit = top_n, required_score = required_score) %>%
      httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
    if (httr2::resp_status(resp) == 404) list(status = 404, body = NULL)
    ## check_type = FALSE: STRING serves "text/json", not "application/json".
    else if (httr2::resp_status(resp) == 200) list(status = 200, body = httr2::resp_body_json(resp, check_type = FALSE, simplifyVector = TRUE))
    else list(status = httr2::resp_status(resp), body = NULL)
  }, error = function(e) NULL)
  if (is.null(out)) return(list(ok = FALSE, reason = "STRING lookup failed (network error or timeout)."))
  if (identical(out$status, 404)) return(list(ok = TRUE, partners = NULL, reason = NULL))
  if (!identical(out$status, 200)) return(list(ok = FALSE, reason = sprintf("STRING lookup failed (HTTP %s).", out$status)))
  res <- out$body
  if (!is.data.frame(res) || nrow(res) == 0) return(list(ok = TRUE, partners = NULL, reason = NULL))
  ## STRING's 7 evidence channels: nscore=neighborhood, fscore=fusion, pscore=cooccurrence,
  ## ascore=coexpression, escore=experimental, dscore=database, tscore=textmining.
  chan_cols <- c("score", "nscore", "fscore", "pscore", "ascore", "escore", "dscore", "tscore")
  df <- res[order(-res$score), c("preferredName_B", intersect(chan_cols, colnames(res))), drop = FALSE]
  colnames(df) <- c("Partner gene", "Combined", "Neighborhood", "Fusion", "Cooccurrence", "Coexpression", "Experimental", "Database", "Textmining")[seq_len(ncol(df))]
  num_cols <- setdiff(colnames(df), "Partner gene")
  df[num_cols] <- lapply(df[num_cols], round, digits = 3)
  list(ok = TRUE, partners = utils::head(df, top_n), reason = NULL)
}

## Fetches the actual network diagram PNG STRING renders server-side (not redrawn locally).
tbc_string_network_image <- function(symbol, top_n = 10, required_score = 400) {
  if (is.null(symbol) || is.na(symbol) || !nzchar(symbol)) return(list(ok = FALSE, reason = "No gene symbol available for STRING image."))
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, reason = "httr2 is not installed in this deployment."))
  resp <- tryCatch({
    httr2::request("https://string-db.org/api/image/network") %>%
      httr2::req_url_query(identifiers = symbol, species = 9606, limit = top_n, required_score = required_score, network_flavor = "confidence") %>%
      httr2::req_timeout(20) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
  }, error = function(e) NULL)
  if (is.null(resp) || httr2::resp_status(resp) != 200) return(list(ok = FALSE, reason = "STRING network image unavailable (network error, timeout, or unresolvable gene symbol)."))
  tf <- tempfile(fileext = ".png")
  writeBin(httr2::resp_body_raw(resp), tf)
  list(ok = TRUE, path = tf, reason = NULL)
}

## DGIdb drug-gene interactions, sorted by DGIdb's interactionScore (not a p-value).
tbc_dgidb_drugs_for_gene <- function(symbol, top_n = 12) {
  if (is.null(symbol) || is.na(symbol) || !nzchar(symbol)) return(list(ok = FALSE, reason = "No gene symbol available for DGIdb lookup."))
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, reason = "httr2 is not installed in this deployment."))
  query <- 'query($name:String!){ genes(names:[$name]) { nodes { name interactions {
    drug { name } interactionScore evidenceScore interactionTypes { type } sources { sourceDbName } } } } }'
  res <- tryCatch({
    resp <- httr2::request("https://dgidb.org/api/graphql") %>%
      httr2::req_headers(`Content-Type` = "application/json") %>%
      httr2::req_body_json(list(query = query, variables = list(name = symbol))) %>%
      httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
    if (httr2::resp_status(resp) != 200) NULL else httr2::resp_body_json(resp, simplifyVector = FALSE)
  }, error = function(e) NULL)
  if (is.null(res)) return(list(ok = FALSE, reason = "DGIdb lookup failed (network error or timeout)."))
  nodes <- res$data$genes$nodes
  if (length(nodes) == 0 || length(nodes[[1]]$interactions) == 0) return(list(ok = TRUE, drugs = NULL, reason = NULL))
  inter <- nodes[[1]]$interactions
  df <- do.call(rbind, lapply(inter, function(x) {
    types <- vapply(x$interactionTypes, function(t) t$type, character(1))
    srcs <- vapply(x$sources, function(s) s$sourceDbName %||% NA_character_, character(1))
    data.frame(Drug = x$drug$name %||% NA_character_, Score = x$interactionScore %||% NA_real_,
               `Evidence (# sources)` = x$evidenceScore %||% NA_integer_,
               `Interaction type` = if (length(types) > 0) paste(types, collapse = ", ") else "Unspecified",
               `Source database(s)` = if (length(srcs) > 0) paste(srcs, collapse = ", ") else "Not recorded",
               check.names = FALSE, stringsAsFactors = FALSE)
  }))
  df <- df[order(-df$Score), , drop = FALSE]
  df$Score <- round(df$Score, 3)
  list(ok = TRUE, drugs = utils::head(df, top_n), reason = NULL)
}

## ---- Literature: PubMed via NCBI E-utilities (free, keyless - the same ----
## esearch/esummary endpoints ArthOChat's own pubmed_search() tool already
## calls, see global.R) - kept structured here (title/authors/journal/year/PMID)
## for table display, rather than pubmed_search()'s flattened citation string.
tbc_literature_search <- function(query, max_results = 12) {
  if (is.null(query) || !nzchar(trimws(query))) return(list(ok = FALSE, papers = NULL, reason = "No search query given."))
  if (!requireNamespace("httr2", quietly = TRUE)) return(list(ok = FALSE, papers = NULL, reason = "httr2 is not installed in this deployment."))
  max_results <- min(max(as.integer(max_results %||% 12), 1L), 30L)
  res <- tryCatch({
    esearch <- httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi") %>%
      httr2::req_url_query(db = "pubmed", term = query, retmode = "json", retmax = max_results, sort = "relevance") %>%
      httr2::req_url_query(tool = "arthomix-biomarkercard") %>%
      httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
    if (httr2::resp_status(esearch) != 200) return(NULL)
    ids <- unlist(httr2::resp_body_json(esearch)$esearchresult$idlist)
    if (length(ids) == 0) return(data.frame())
    esummary <- httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi") %>%
      httr2::req_url_query(db = "pubmed", id = paste(ids, collapse = ","), retmode = "json") %>%
      httr2::req_timeout(15) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform()
    if (httr2::resp_status(esummary) != 200) return(NULL)
    recs <- httr2::resp_body_json(esummary)$result
    rows <- lapply(ids, function(pmid) {
      r <- recs[[pmid]]
      if (is.null(r)) return(NULL)
      authors <- vapply(r$authors, function(a) a$name %||% "?", character(1))
      author_str <- if (length(authors) > 1) paste0(authors[1], " et al.") else if (length(authors) == 1) authors[1] else "Unknown author"
      data.frame(Title = r$title %||% "Untitled", Authors = author_str, Journal = r$fulljournalname %||% r$source %||% NA_character_,
                 Year = substr(r$pubdate %||% "", 1, 4), PMID = pmid, stringsAsFactors = FALSE)
    })
    do.call(rbind, Filter(Negate(is.null), rows))
  }, error = function(e) NULL)
  if (is.null(res)) return(list(ok = FALSE, papers = NULL, reason = "PubMed lookup failed (network error, timeout, or malformed response)."))
  list(ok = TRUE, papers = if (is.data.frame(res) && nrow(res) > 0) res else NULL, reason = NULL)
}

## Preset query templates (spec: gene+biomarker, gene+pathway, gene+sex,
## gene+transcriptomics, gene+rheumatoid arthritis; "disease" uses whatever
## case/control context the loaded cohort provides, or free text).
TBC_LITERATURE_PRESETS <- c(
  "Biomarker" = "%s biomarker", "Pathway" = "%s pathway", "Sex differences" = "%s sex",
  "Transcriptomics" = "%s transcriptomics", "Rheumatoid arthritis" = "%s rheumatoid arthritis"
)
tbc_literature_query <- function(gene, preset, custom_disease = NULL) {
  if (identical(preset, "disease")) {
    ctx <- if (!is.null(custom_disease) && nzchar(trimws(custom_disease))) trimws(custom_disease) else "disease"
    return(sprintf("%s %s", gene, ctx))
  }
  tmpl <- TBC_LITERATURE_PRESETS[[preset]] %||% "%s"
  sprintf(tmpl, gene)
}

## Renders a real KEGG pathway diagram (pathview) with this gene colored by its log2FC.
tbc_kegg_diagram_for_gene <- function(kegg_id, entrez, log2fc, pathway_name = NULL) {
  if (is.null(kegg_id) || !nzchar(kegg_id)) return(list(ok = FALSE, path = NULL, error = "No KEGG pathway selected."))
  if (is.null(entrez) || is.na(entrez) || !nzchar(entrez)) return(list(ok = FALSE, path = NULL, error = "No NCBI Gene ID available to color onto the pathway map."))
  val <- suppressWarnings(as.numeric(log2fc))
  if (length(val) == 0 || is.na(val)) val <- 0
  effect <- stats::setNames(val, entrez)
  out <- tryCatch(mp_kegg_pathway_map(kegg_id, effect, out_dir = tempdir()), error = function(e) list(ok = FALSE, path = NULL, error = conditionMessage(e)))
  c(out, list(pathway_id = kegg_id, pathway_name = pathway_name))
}

tbc_reactome_diagram_for_pathway <- function(stable_id, pathway_name = NULL) {
  if (is.null(stable_id) || !nzchar(stable_id)) return(list(ok = FALSE, path = NULL, error = "No Reactome pathway selected."))
  out <- tryCatch(mp_fetch_reactome_diagram_png(stable_id, out_dir = tempdir()), error = function(e) list(ok = FALSE, path = NULL, error = conditionMessage(e)))
  c(out, list(pathway_id = stable_id, pathway_name = pathway_name))
}

## ---- Sample-sheet column detection ----

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

## ---- Expression-scale-aware live evidence ----

## Quick single-gene preview only, not a substitute for the DGE tab's limma/DESeq2 model.
## Raw counts get log2(CPM+1)-transformed (edgeR::cpm()); already-normalised data is used as-is.
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

## Looks up this gene in every saved DGE run (results$dge_runs, written live by mod_dge.R).
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

## ---- Single-gene diagnostic performance (live, computed here) ------------
## Answers "is this gene useful alone as a classifier?" using the same
## statistical primitive (pROC) mod_diagnostic.R already depends on, computed
## directly on the live per-gene values tbc_dataset_evidence() already reads -
## no other module's reactives are touched. Training = full-fit on the whole
## loaded dataset; tbc_single_gene_cv() below is the "Internal Validation"
## counterpart (pooled out-of-fold predictions from k-fold CV).
tbc_single_gene_roc <- function(expr_row, group_vec, case_label, control_label) {
  if (!requireNamespace("pROC", quietly = TRUE)) return(list(ok = FALSE, reason = "pROC is not installed in this deployment."))
  keep <- !is.na(expr_row) & !is.na(group_vec) & group_vec %in% c(case_label, control_label)
  x <- expr_row[keep]; y <- group_vec[keep]
  n_case <- sum(y == case_label); n_control <- sum(y == control_label)
  if (n_case < 3 || n_control < 3) {
    return(list(ok = FALSE, reason = sprintf("Fewer than 3 samples in one of the two groups (case=%d, control=%d) - cannot compute a reliable single-gene ROC curve.", n_case, n_control)))
  }
  roc_obj <- tryCatch(pROC::roc(y, x, levels = c(control_label, case_label), direction = "auto", quiet = TRUE), error = function(e) NULL)
  if (is.null(roc_obj)) return(list(ok = FALSE, reason = "pROC could not fit a ROC curve for this gene (e.g. constant expression)."))
  ci <- tryCatch(as.numeric(pROC::ci.auc(roc_obj, quiet = TRUE)), error = function(e) c(NA_real_, NA_real_, NA_real_))
  best <- tryCatch(pROC::coords(roc_obj, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = FALSE), error = function(e) NULL)
  if (is.null(best) || nrow(best) == 0) return(list(ok = FALSE, reason = "Could not determine an optimal threshold for this gene."))
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

tbc_single_gene_cv <- function(expr_row, group_vec, case_label, control_label, k = 5, seed = 1234) {
  if (!requireNamespace("pROC", quietly = TRUE)) return(list(ok = FALSE, reason = "pROC is not installed in this deployment."))
  keep <- !is.na(expr_row) & !is.na(group_vec) & group_vec %in% c(case_label, control_label)
  x <- expr_row[keep]; y <- group_vec[keep]
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
  if (sum(usable) < 6) return(list(ok = FALSE, reason = "Cross-validation could not produce enough out-of-fold predictions (folds too small, or this gene wasn't separable in any training fold)."))
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

## Precision-recall from an already-fit pROC ROC object (recall = sensitivity,
## precision = PPV) - no new package, PR-AUC via trapezoidal integration.
tbc_pr_from_roc <- function(roc_obj) {
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

tbc_plot_roc <- function(roc_obj, label = "ROC curve", auc = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_) {
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

tbc_plot_pr <- function(pr, label = "Precision-recall curve") {
  if (is.null(pr) || !isTRUE(pr$ok)) return(NULL)
  ggplot(pr$df, aes(x = recall, y = precision)) +
    geom_line(color = ARTHOMIX_COLORS$violet, linewidth = 1) +
    labs(title = label, subtitle = sprintf("PR-AUC %.3f", pr$pr_auc), x = "Recall (sensitivity)", y = "Precision (PPV)") +
    ylim(0, 1) + theme_arthomix()
}

tbc_fmt_num <- function(x, digits = 3) if (is.null(x) || length(x) == 0 || is.na(x)) "Not available" else sprintf(sprintf("%%.%df", digits), x)

tbc_confusion_table <- function(cm) {
  if (is.null(cm)) return(NULL)
  df <- as.data.frame.matrix(cm)
  df <- cbind(Predicted = rownames(df), df)
  DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE, ordering = FALSE), class = "stripe hover compact")
}

## ---- Candidate Gene Identification / ML Feature Selection membership ----

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

## Reads classifier performance verbatim from results$diagnostic for panels containing this gene.
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

## Volcano plot from an already-saved DGE run table, highlighting the selected gene.
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

## Renders a link to the live source page so numbers can be checked against the real database.
tbc_db_provenance <- function(api_domain, live_url, live_label) {
  div(class = "empty-note", style = "display:flex; align-items:center; justify-content:space-between; gap:10px; flex-wrap:wrap;",
      tagList(icon("satellite-dish"), sprintf("Live query to %s - nothing here is precomputed or cached in this app.", api_domain)),
      tags$a(href = live_url, target = "_blank", rel = "noopener", class = "btn btn-sm btn-default",
             icon("arrow-up-right-from-square"), sprintf(" Open %s ", live_label))
  )
}

## Static light-blue info box for the External Databases tab (spec §9) -
## reuses the app's existing (previously unused) .data-source-callout class
## from www/custom.css rather than inventing new CSS/colors.
TBC_EXTERNAL_DB_NAMES <- c(
  "GO", "KEGG", "Reactome", "WikiPathways", "Disease / Genetics — Open Targets",
  "Drug / Target — DGIdb", "Protein / Interaction — STRING", "Expression — Human Protein Atlas", "Literature — PubMed"
)
tbc_external_db_banner <- function() {
  div(class = "data-source-callout",
      icon("circle-info"),
      tagList(tags$b("External databases queried live for the selected gene: "), paste(TBC_EXTERNAL_DB_NAMES, collapse = " | "),
              tags$div(style = "margin-top:4px;", "Each is opt-in below - pick one, click Run, and it queries the real database directly. A result you don't see means the database itself returned nothing, or hasn't been run yet - never a fabricated value."))
  )
}

## ---- Section builders (plain tags, reused on-screen and in the downloadable report) ----

## Full gene-identity card: local org.Hs.eg.db fields always shown; the two
## live lookups (NCBI summary, UniProt protein name) degrade to "Not available"
## independently if either fails, without blocking the rest of the card.
tbc_section_identity <- function(d) {
  gi <- d$gene_identity
  ncbi <- d$ncbi_summary
  uni <- d$uniprot_name
  pairs <- list(
    "Gene symbol" = d$gene,
    "Full name" = if (isTRUE(gi$ok)) gi$genename else NA,
    "NCBI Entrez Gene ID" = if (isTRUE(gi$ok)) gi$entrez else NA,
    "Ensembl Gene ID" = if (isTRUE(gi$ok)) gi$ensembl else NA,
    "Chromosome" = if (isTRUE(gi$ok)) gi$chr else NA,
    "Genomic location (cytogenetic band)" = if (isTRUE(gi$ok)) gi$map_location else NA,
    "Gene type" = if (isTRUE(gi$ok)) gi$genetype else NA,
    "Aliases" = if (isTRUE(gi$ok) && length(gi$aliases) > 0) paste(gi$aliases, collapse = ", ") else NA,
    "Protein name" = if (isTRUE(uni$ok)) uni$name else NA,
    "Basic description" = if (isTRUE(ncbi$ok)) ncbi$description else NA,
    "Present in currently loaded expression matrix" = if (d$in_dataset) "Yes" else "No"
  )
  links <- tbc_gene_identity_links(d$gene, if (isTRUE(gi$ok)) gi$entrez else NULL, if (isTRUE(gi$ok)) gi$ensembl else NULL)
  link_items <- Filter(Negate(is.null), lapply(names(links), function(nm) {
    url <- links[[nm]]
    if (is.na(url)) return(NULL)
    tags$a(href = url, target = "_blank", rel = "noopener", class = "btn btn-sm btn-default", style = "margin-right:6px; margin-bottom:6px;",
           icon("arrow-up-right-from-square"), sprintf(" %s", nm))
  }))
  div(class = "card",
      div(class = "card-title", icon("dna"), "Gene description"),
      if (!isTRUE(gi$ok)) div(class = "empty-note", icon("triangle-exclamation"), gi$reason) else NULL,
      tbc_kv_table(pairs),
      if (!isTRUE(ncbi$ok) && !is.null(ncbi)) div(class = "empty-note", style = "margin-top:8px;", icon("circle-info"), sprintf("Basic description: %s", ncbi$reason %||% "unavailable.")) else NULL,
      if (!isTRUE(uni$ok) && !is.null(uni)) div(class = "empty-note", style = "margin-top:8px;", icon("circle-info"), sprintf("Protein name: %s", uni$reason %||% "unavailable.")) else NULL,
      if (length(link_items) > 0) tagList(tags$div(style = "margin-top:10px;", tags$b("Authoritative gene resources")), div(style = "margin-top:6px;", link_items)) else NULL
  )
}

## Merges what used to be two separate cards (Expression Data Input,
## Transcriptomic Biomarker Discovery) - same dataset/cohort, so one card.
## GEO accession / platform are never stored as separate fields anywhere in
## this app (confirmed: dataset$source is free text, dataset$meta has no
## platform/tissue/species column) - parsed out of dataset$source only when
## the exact "NCBI GEO: <acc> (<platform>, ...)" pattern mod_dataset.R writes
## is actually present, never guessed otherwise.
tbc_parse_geo_source <- function(source_label) {
  if (is.null(source_label) || !nzchar(source_label)) return(list(accession = NA_character_, platform = NA_character_))
  m <- regmatches(source_label, regexec("NCBI GEO:\\s*(\\S+)\\s*\\(([^,)]+)", source_label))[[1]]
  if (length(m) == 3) list(accession = m[2], platform = m[3]) else list(accession = NA_character_, platform = NA_character_)
}

tbc_section_dataset_cohort <- function(dataset, live) {
  n_genes <- tryCatch(nrow(dataset$expr), error = function(e) NA_integer_)
  n_samples <- tryCatch(ncol(dataset$expr), error = function(e) NA_integer_)
  data_type <- if (isTRUE(live$ok) || !is.null(live$data_type)) live$data_type else tryCatch(detect_expr_data_type(dataset$expr), error = function(e) NA_character_)
  type_label <- switch(data_type %||% "unknown",
    counts = "Raw counts (log2-CPM transformed for this preview) - preprocessed/normalised for analysis",
    already_normalised = "Already normalised / log-scale",
    expression = "Linear/log expression, not yet quantile-normalised",
    "Not determined"
  )
  group_label <- if (isTRUE(live$ok) || !is.null(live$group_col)) live$group_col else "Not auto-detected"
  case_label <- if (isTRUE(live$ok)) sprintf("%s vs %s", live$case_label, live$control_label) else "Not available"
  sex_label <- if (isTRUE(live$ok) && !is.null(live$sex_col)) live$sex_col else "Not recorded"
  n_case <- if (isTRUE(live$ok) && isTRUE(live$overall$ok)) live$overall$n_case else NA
  n_ctrl <- if (isTRUE(live$ok) && isTRUE(live$overall$ok)) live$overall$n_control else NA
  sex_breakdown <- if (isTRUE(live$ok) && !is.null(live$sex_vec)) {
    tb <- table(live$sex_vec[!is.na(live$sex_vec)])
    paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = ", ")
  } else NA_character_
  geo <- tbc_parse_geo_source(dataset$source)
  pairs <- list(
    "Dataset / cohort label" = dataset$source,
    "GEO accession" = geo$accession,
    "Platform" = geo$platform,
    "Genes x samples" = if (!is.na(n_genes)) sprintf("%s x %s", format(n_genes, big.mark = ","), n_samples) else NA,
    "Number of samples" = n_samples,
    "Data type / preprocessing status" = type_label,
    "Species" = "Not available (not recorded for this dataset)",
    "Tissue / sample type" = "Not available (not recorded for this dataset)",
    "Analysis group / contrast" = sprintf("%s (grouping column: %s)", case_label, group_label),
    "Cases / Controls" = if (!is.na(n_case)) sprintf("%s / %s", n_case, n_ctrl) else NA,
    "Sex / group breakdown" = sex_breakdown,
    "Dataset designation" = "Training / discovery dataset (this session's loaded dataset). Internal validation = cross-validation folds within it; external validation = a separate cohort, not available in this session's shared results.",
    "Feature identifier" = "Gene symbol (rows of the loaded expression matrix)",
    "Gene mapping" = "org.Hs.eg.db (NCBI Entrez / Ensembl Gene ID resolution) - see Gene description tab"
  )
  div(class = "card",
      div(class = "card-title", icon("table-cells"), "Dataset & Cohort"),
      p(class = "submodule-desc", "Where this biomarker's evidence comes from. Fields with no equivalent recorded anywhere in this app's dataset metadata are shown as \"Not available\", never guessed."),
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
      Gene = d$gene, `Feature/Probe ID` = "Not applicable (gene-symbol level)",
      `Log2 fold-change` = round(ov$log2fc, 3), `Average expression` = "Not available", `P-value` = signif(ov$p_value, 4), `Adjusted p-value` = NA_real_,
      Direction = ov$direction, Rank = "Not available (no saved Differential Expression run loaded)",
      `Biomarker candidate status` = if (!is.na(ov$p_value) && ov$p_value <= 0.05) "Candidate (uncorrected p <= 0.05)" else "Not significant",
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (!is.null(dge_hits) && nrow(dge_hits) > 0) {
    for (i in seq_len(nrow(dge_hits))) {
      r <- dge_hits[i, ]
      rk <- if (identical(r$run, d$selected_run_id) && !is.null(d$selected_run_table)) tbc_de_rank(d$gene, d$selected_run_table) else list(rank = NA_integer_, n = NA_integer_)
      rows[[length(rows) + 1]] <- data.frame(
        Source = sprintf("Saved Differential Expression run (%s, %s)", r$contrast, r$method),
        Gene = d$gene, `Feature/Probe ID` = "Not applicable (gene-symbol level)",
        `Log2 fold-change` = round(r$logFC, 3), `Average expression` = "Not available in the saved run", `P-value` = "Not available in the saved run", `Adjusted p-value` = signif(r$adj.P.Val, 4),
        Direction = r$direction, Rank = if (!is.na(rk$rank)) sprintf("%s of %s", rk$rank, rk$n) else "Not available (open Differential Expression tab)",
        `Biomarker candidate status` = if (identical(r$direction, "Not significant")) "Not significant" else "Candidate transcriptomic biomarker",
        check.names = FALSE, stringsAsFactors = FALSE)
    }
  }
  body <- if (length(rows) == 0) {
    div(class = "empty-note", icon("circle-info"),
        "No differential-expression evidence is available for this gene yet - it is either absent from the loaded expression matrix, the grouping column could not be auto-detected, or no Differential Expression run this session included it. Run the Differential Expression tab, or check the Dataset tab.")
  } else {
    DT::datatable(do.call(rbind, rows), rownames = FALSE, options = list(dom = "ft", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")
  }
  div(class = "card",
      div(class = "card-title", icon("chart-column"), "Differential Expression Table"),
      p(class = "submodule-desc", "Preview: uncorrected p. Saved runs: that run's own adjusted p / log2FC cutoff, rank computed by ordering that run's genes by FDR."),
      body
  )
}

## One status pill (done/pending/neutral), using the shared .pipeline-status-chip CSS.
tbc_status_chip <- function(label, state, value = NULL) {
  cls <- switch(state, done = "status-done", pending = "status-pending", "status-neutral")
  ic <- switch(state, done = "circle-check", pending = "circle-xmark", "circle-minus")
  span(class = paste("pipeline-status-chip", cls), icon(ic),
       if (!is.null(value)) sprintf("%s: %s", label, value) else label)
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
      p(class = "submodule-desc", "Sex-stratified consensus signature membership (ML Feature Selection: LASSO/RF/SVM-RFE agreement)."),
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE), class = "stripe hover compact"),
      if (in_any) tagList(full_lists) else div(class = "empty-note", icon("circle-info"), "This gene is not part of any consensus signature computed so far this session.")
  )
}

## Compact multi-gene panel table: best model per stratum, its stored
## full-fit ("Training") AUC and cross-validated ("Internal Validation") AUC.
## Sensitivity/specificity/threshold/confusion matrix are NOT stored in
## results$diagnostic (mod_diagnostic.R keeps those in its own private
## reactives) - reported honestly as "Not available" rather than guessed.
tbc_multi_gene_perf_table <- function(diag) {
  if (length(diag) == 0) return(NULL)
  do.call(rbind, lapply(names(diag), function(s) {
    x <- diag[[s]]
    aucs <- c(lr = x$lr_cv_auc, enet = x$enet_cv_auc, rf = x$rf_cv_auc, svm = x$svm_cv_auc)
    train_aucs <- c(lr = x$lr_auc, enet = x$enet_auc, rf = x$rf_auc, svm = x$svm_auc)
    best_key <- if (all(is.na(aucs))) NA_character_ else names(aucs)[which.max(aucs)]
    best_label <- if (is.na(best_key)) "Not available" else c(lr = "Logistic Regression", enet = "Elastic Net", rf = "Random Forest", svm = "SVM")[[best_key]]
    data.frame(
      Stratum = tools::toTitleCase(s), `Panel size` = x$panel_size, `Best model` = best_label,
      `Training AUC (full fit)` = if (is.na(best_key)) "Not available" else tbc_fmt_field(train_aucs[[best_key]]),
      `Internal validation AUC (CV)` = if (is.na(best_key)) "Not available" else tbc_fmt_field(aucs[[best_key]]),
      Sensitivity = "Not available (only AUC is stored for saved panels)",
      Specificity = "Not available (only AUC is stored for saved panels)",
      check.names = FALSE, stringsAsFactors = FALSE)
  }))
}

## ---- Single-vs-multi-gene comparison (tab: "Single-Gene vs Multi-Gene ------
## Signature") - answers "is this gene useful alone?" with a real, live,
## honestly-labeled single-gene classifier next to whatever multi-gene panel
## evidence already exists in the shared results store.
tbc_section_signature_comparison <- function(d, sgd, roc_widget = NULL) {
  single_body <- if (isTRUE(sgd$ok)) {
    tbc_kv_table(list(
      AUC = tbc_fmt_num(sgd$auc), "95% CI" = sprintf("%s - %s", tbc_fmt_num(sgd$ci_lo), tbc_fmt_num(sgd$ci_hi)),
      Threshold = tbc_fmt_num(sgd$threshold), Sensitivity = tbc_fmt_num(sgd$sensitivity), Specificity = tbc_fmt_num(sgd$specificity),
      Accuracy = tbc_fmt_num(sgd$accuracy), `Balanced accuracy` = tbc_fmt_num(sgd$balanced_accuracy),
      `n (case / control)` = sprintf("%d / %d", sgd$n_case, sgd$n_control)
    ))
  } else div(class = "empty-note", icon("circle-info"), sgd$reason %||% "Single-gene performance unavailable.")

  multi_table <- tbc_multi_gene_perf_table(d$diagnostic_match)
  multi_body <- if (is.null(multi_table)) {
    div(class = "empty-note", icon("circle-info"), "This gene is not part of any Diagnostic Classifier panel evaluated this session. Run the Diagnostic Classifier tab on a panel that includes it (e.g. its consensus signature below) to populate this section.")
  } else DT::datatable(multi_table, rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")

  comparison <- NULL
  if (isTRUE(sgd$ok) && !is.null(multi_table)) {
    best_multi_cv <- suppressWarnings(max(as.numeric(gsub("[^0-9.]", "", multi_table$`Internal validation AUC (CV)`)), na.rm = TRUE))
    if (is.finite(best_multi_cv)) {
      diff <- best_multi_cv - sgd$auc
      comparison <- div(class = "empty-note", icon("scale-balanced"),
        if (abs(diff) < 0.01) sprintf("The best multi-gene panel's internal-validation AUC (%.3f) is essentially the same as this gene alone (%.3f) - no meaningful signature benefit is shown by the numbers here.", best_multi_cv, sgd$auc)
        else if (diff > 0) sprintf("The best multi-gene panel's internal-validation AUC (%.3f) is %.3f points higher than this gene alone (%.3f) - the signature adds real, measured value here.", best_multi_cv, diff, sgd$auc)
        else sprintf("This gene alone (AUC %.3f) actually scores %.3f points higher than the best multi-gene panel's internal-validation AUC (%.3f) in this session's saved runs.", sgd$auc, -diff, best_multi_cv))
    }
  }

  tagList(
    div(class = "card",
        div(class = "card-title", icon("layer-group"), "Single-Gene Performance"),
        p(class = "submodule-desc", "Live computation: this gene's own expression used as a single classifier on the currently loaded dataset (\"Training\" / full-fit sense - see Biomarker Performance tab for internal-validation numbers)."),
        single_body,
        if (!is.null(roc_widget)) tagList(tags$div(style = "margin-top:10px;"), roc_widget) else NULL
    ),
    div(class = "card",
        div(class = "card-title", icon("layer-group"), "Multi-Gene Signature Performance"),
        p(class = "submodule-desc", "From this session's saved Diagnostic Classifier runs - never recomputed here."),
        multi_body
    ),
    tbc_section_signature(d),
    if (!is.null(comparison)) div(class = "card", div(class = "card-title", icon("scale-balanced"), "Comparison"), comparison) else NULL
  )
}

## ---- Expression / Diagnostic / Validation evidence cards -------------------
## (feed the "Biomarker Status" tab; §4/§12 of the spec)

tbc_de_rank <- function(gene, table) {
  if (is.null(table) || nrow(table) == 0) return(list(rank = NA_integer_, n = NA_integer_))
  ord <- table[order(table$adj.P.Val), , drop = FALSE]
  idx <- which(ord$gene == gene)
  list(rank = if (length(idx) > 0) idx[1] else NA_integer_, n = nrow(table))
}

tbc_section_expression_evidence <- function(d) {
  row_run <- if (!is.null(d$selected_run_table)) d$selected_run_table[d$selected_run_table$gene == d$gene, , drop = FALSE] else NULL
  if (!is.null(row_run) && nrow(row_run) > 0) {
    rk <- tbc_de_rank(d$gene, d$selected_run_table)
    pairs <- list(
      "Log2 fold-change" = round(row_run$logFC[1], 3),
      "Adjusted p-value (FDR)" = signif(row_run$adj.P.Val[1], 4),
      "Raw p-value" = "Not available (only gene/logFC/adjusted-p/direction are kept in a saved run)",
      "Direction" = row_run$direction[1],
      "Effect size (|log2FC|)" = round(abs(row_run$logFC[1]), 3),
      "Rank among tested genes (by FDR)" = if (!is.na(rk$rank)) sprintf("%s of %s", rk$rank, rk$n) else "Not available"
    )
  } else if (isTRUE(d$live$ok) && isTRUE(d$live$overall$ok)) {
    ov <- d$live$overall
    pairs <- list(
      "Log2 fold-change (live preview, uncorrected)" = round(ov$log2fc, 3),
      "Adjusted p-value (FDR)" = "Not computed (live preview is an uncorrected quick check - run the Differential Expression tab for FDR)",
      "Raw p-value" = signif(ov$p_value, 4),
      "Direction" = ov$direction,
      "Effect size (|log2FC|)" = round(abs(ov$log2fc), 3),
      "Rank among tested genes" = "Not available (no saved Differential Expression run loaded)"
    )
  } else {
    return(div(class = "card", div(class = "card-title", icon("chart-column"), "Expression Evidence"),
               div(class = "empty-note", icon("circle-info"), "No differential-expression evidence is available for this gene yet.")))
  }
  div(class = "card", div(class = "card-title", icon("chart-column"), "Expression Evidence"), tbc_kv_table(pairs))
}

tbc_section_diagnostic_evidence <- function(d, sgd) {
  single_body <- if (isTRUE(sgd$ok)) {
    tbc_kv_table(list(AUC = tbc_fmt_num(sgd$auc), Sensitivity = tbc_fmt_num(sgd$sensitivity), Specificity = tbc_fmt_num(sgd$specificity),
                       Accuracy = tbc_fmt_num(sgd$accuracy), Threshold = tbc_fmt_num(sgd$threshold)))
  } else div(class = "empty-note", icon("circle-info"), sgd$reason %||% "Not available.")
  multi_table <- tbc_multi_gene_perf_table(d$diagnostic_match)
  multi_body <- if (is.null(multi_table)) div(class = "empty-note", icon("circle-info"), "Not evaluated in any saved Diagnostic Classifier panel this session.")
                else DT::datatable(multi_table[, c("Stratum", "Panel size", "Best model", "Internal validation AUC (CV)")],
                                   rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")
  div(class = "card", div(class = "card-title", icon("gauge-high"), "Diagnostic Evidence"),
      tags$b("Single-gene (live, this dataset)"), single_body,
      tags$div(style = "margin-top:10px;", tags$b("Multi-gene panel")), multi_body)
}

tbc_section_validation_evidence <- function(sgd, sgcv) {
  rows <- data.frame(
    Stage = c("Training", "Internal Validation", "External Validation"),
    Status = c(if (isTRUE(sgd$ok)) "Available" else "Not available", if (isTRUE(sgcv$ok)) "Available" else "Not available", "Not available"),
    AUC = c(tbc_fmt_num(sgd$auc), tbc_fmt_num(sgcv$auc), "Not available"),
    Note = c("Single-gene, full-fit on the whole currently loaded dataset.",
             if (isTRUE(sgcv$ok)) sprintf("Single-gene, %d-fold cross-validation, pooled out-of-fold predictions.", sgcv$k) else (sgcv$reason %||% "Not available."),
             "Not persisted to this app's shared results store - open the Diagnostic Classifier tab's own External Validation panel to run and view it directly."),
    stringsAsFactors = FALSE
  )
  div(class = "card", div(class = "card-title", icon("shield-halved"), "Validation Evidence"),
      p(class = "submodule-desc", "Training, internal, and external validation are always kept separate below - never merged into one number."),
      DT::datatable(rows, rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact"))
}

## ---- Evidence tier classification (pure, unit-testable; §13) --------------
## Exactly the spec's 3-tier scale, plus an honest 4th "Insufficient evidence"
## floor for a gene with no significant DE at all (never forced into one of
## the 3 named tiers). The tier is always returned alongside the literal
## per-domain checklist that produced it - never a bare, unexplained badge.
tbc_evidence_classification <- function(d, ext, sgd = NULL, sgcv = NULL) {
  statistical_ok <- (isTRUE(d$live$ok) && isTRUE(d$live$overall$ok) && !is.na(d$live$overall$p_value) && d$live$overall$p_value <= 0.05) ||
    (!is.null(d$dge_hits) && any(d$dge_hits$direction != "Not significant"))
  diagnostic_ok <- isTRUE(sgd$ok) || length(d$diagnostic_match %||% list()) > 0
  biological_ok <- !is.null(ext) && (
    (!is.null(ext$go) && nrow(ext$go) > 0) ||
    any(vapply(TBC_EVIDENCE_DBS, function(x) identical(tbc_evidence_status(ext[[x$key]], x$field), "Results found"), logical(1)))
  )
  disease_ok <- isTRUE(ext$genetics$ok) && !is.na(ext$genetics$n_diseases %||% NA) && (ext$genetics$n_diseases %||% 0) > 0
  therapeutic_ok <- isTRUE(ext$drugs$ok) && !is.null(ext$drugs$drugs) && nrow(ext$drugs$drugs) > 0
  validation_training_ok <- diagnostic_ok
  cv_auc_present <- length(d$diagnostic_match %||% list()) > 0 &&
    any(!is.na(unlist(lapply(d$diagnostic_match, function(x) c(x$lr_cv_auc, x$enet_cv_auc, x$rf_cv_auc, x$svm_cv_auc)))))
  validation_internal_ok <- isTRUE(sgcv$ok) || cv_auc_present
  validation_external_ok <- FALSE ## never persisted to the shared results store in this deployment - see plan/spec notes.

  tier <- if (!statistical_ok) "Insufficient evidence"
          else if (validation_external_ok) "Strong candidate"
          else if ((diagnostic_ok || biological_ok) && validation_internal_ok) "Supported candidate"
          else "Candidate biomarker"

  checklist <- list(
    list(label = "Significant differential expression", met = statistical_ok),
    list(label = "Diagnostic evidence (single-gene or panel)", met = diagnostic_ok),
    list(label = "Biological evidence (pathway/GO/disease/drug/PPI/expression/literature)", met = biological_ok),
    list(label = "Disease association evidence", met = disease_ok),
    list(label = "Therapeutic/drug-target evidence", met = therapeutic_ok),
    list(label = "Training evidence", met = validation_training_ok),
    list(label = "Internal validation (cross-validation)", met = validation_internal_ok),
    list(label = "External validation", met = validation_external_ok)
  )
  list(tier = tier, checklist = checklist)
}

TBC_TIER_CLASS <- c("Strong candidate" = "status-done", "Supported candidate" = "status-pending",
                    "Candidate biomarker" = "status-pending", "Insufficient evidence" = "status-neutral")

## ---- "Biomarker Status" tab: the main evidence summary (§4) ---------------
tbc_section_evidence_glance <- function(d, ext, sgd, sgcv) {
  cl <- tbc_evidence_classification(d, ext, sgd, sgcv)
  n_met <- sum(vapply(cl$checklist, function(it) isTRUE(it$met), logical(1)))
  tier_card <- div(class = "card",
      div(class = "card-title", icon("certificate"), "Evidence Status"),
      div(class = "evidence-tier-banner",
          span(class = paste("pipeline-status-chip", TBC_TIER_CLASS[[cl$tier]]), icon("certificate"), sprintf("Evidence tier: %s", cl$tier)),
          span(class = "evidence-tier-count", sprintf("%d of %d evidence domains supported", n_met, length(cl$checklist)))
      ),
      div(class = "evidence-checklist", lapply(cl$checklist, function(it) {
        state <- if (grepl("External validation", it$label) && !it$met) "neutral" else if (it$met) "done" else "pending"
        ic <- switch(state, done = "circle-check", pending = "circle", "circle-minus")
        div(class = paste("evidence-checklist-item", paste0("is-", state)), icon(ic), span(it$label))
      }))
  )
  tagList(
    tier_card,
    tbc_section_expression_evidence(d),
    tbc_section_diagnostic_evidence(d, sgd),
    tbc_section_validation_evidence(sgd, sgcv)
  )
}

## ---- Biomarker Performance tab: Training / Internal / External kept -------
## visibly separate throughout, never overwritten or combined (§8).
tbc_section_biomarker_performance <- function(d, sgd, sgcv, train_roc_widget = NULL, train_pr_widget = NULL, internal_roc_widget = NULL) {
  train_body <- if (isTRUE(sgd$ok)) {
    tagList(
      tbc_kv_table(list(AUC = tbc_fmt_num(sgd$auc), "95% CI" = sprintf("%s - %s", tbc_fmt_num(sgd$ci_lo), tbc_fmt_num(sgd$ci_hi)),
                         Threshold = tbc_fmt_num(sgd$threshold), Sensitivity = tbc_fmt_num(sgd$sensitivity), Specificity = tbc_fmt_num(sgd$specificity),
                         Accuracy = tbc_fmt_num(sgd$accuracy), `Balanced accuracy` = tbc_fmt_num(sgd$balanced_accuracy),
                         PPV = tbc_fmt_num(sgd$ppv), NPV = tbc_fmt_num(sgd$npv), `n (case / control)` = sprintf("%d / %d", sgd$n_case, sgd$n_control))),
      tags$div(style = "margin-top:10px;", tags$b("Confusion matrix (single-gene, at optimal threshold)")), tbc_confusion_table(sgd$confusion),
      if (!is.null(train_roc_widget)) tagList(tags$div(style = "margin-top:10px;", tags$b("ROC curve")), train_roc_widget) else NULL,
      if (!is.null(train_pr_widget)) tagList(tags$div(style = "margin-top:10px;", tags$b("Precision-recall curve")), train_pr_widget) else NULL
    )
  } else div(class = "empty-note", icon("circle-info"), sgd$reason %||% "Not available.")

  internal_body <- if (isTRUE(sgcv$ok)) {
    tagList(
      tbc_kv_table(list(AUC = tbc_fmt_num(sgcv$auc), Sensitivity = tbc_fmt_num(sgcv$sensitivity), Specificity = tbc_fmt_num(sgcv$specificity),
                         Accuracy = tbc_fmt_num(sgcv$accuracy), `Balanced accuracy` = tbc_fmt_num(sgcv$balanced_accuracy),
                         `Folds (k)` = sgcv$k, `n used` = sgcv$n_used)),
      tags$div(style = "margin-top:10px;", tags$b("Confusion matrix (single-gene, pooled out-of-fold predictions)")), tbc_confusion_table(sgcv$confusion),
      if (!is.null(internal_roc_widget)) tagList(tags$div(style = "margin-top:10px;", tags$b("ROC curve (pooled out-of-fold)")), internal_roc_widget) else NULL
    )
  } else div(class = "empty-note", icon("circle-info"), sgcv$reason %||% "Not available.")

  external_body <- div(class = "empty-note", icon("triangle-exclamation"),
    "No data available in this session. The Diagnostic Classifier tab has its own External Validation upload panel, but its results are not persisted to this app's shared results store, so this Biomarker Card cannot display them here - open that tab directly to run and view external validation for a chosen gene panel.")

  multi_table <- tbc_multi_gene_perf_table(d$diagnostic_match)

  div(class = "card",
      div(class = "card-title", icon("vial-circle-check"), "Biomarker Performance"),
      p(class = "submodule-desc", "Single-gene metrics are computed live on the currently loaded dataset; multi-gene panel numbers are read verbatim from the Diagnostic Classifier tab's saved runs."),
      tags$h5("Training (full-fit)"), train_body,
      tags$hr(),
      tags$h5("Internal Validation (cross-validation)"), internal_body,
      tags$hr(),
      tags$h5("External Validation"), external_body,
      tags$hr(),
      tags$b("Calibration"), div(class = "empty-note", icon("circle-info"), "Not applicable - a single raw expression value is not a calibrated probability model."),
      if (!is.null(multi_table)) tagList(tags$hr(), tags$b("Multi-gene panel context (this gene is included in an evaluated panel)"),
        DT::datatable(multi_table, rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")) else NULL
  )
}

## ---- Gene Ontology (biological process) -----------------------------------
tbc_section_go <- function(ext) {
  if (is.null(ext)) return(div(class = "card", div(class = "card-title", icon("circle-nodes"), "Gene Ontology"),
                                div(class = "empty-note", icon("circle-info"), "Not yet looked up - click Run below.")))
  go_body <- if (!is.null(ext$go) && nrow(ext$go) > 0) {
    go_df <- ext$go
    go_df$Link <- sprintf('<a href="https://www.ebi.ac.uk/QuickGO/term/%s" target="_blank" rel="noopener">%s</a>', go_df$GOID, go_df$GOID)
    DT::datatable(go_df[, c("Link", "TERM")], colnames = c("GO ID", "Biological process"), rownames = FALSE, escape = FALSE,
                  options = list(dom = "t", paging = FALSE), class = "stripe hover compact")
  } else div(class = "empty-note", icon("circle-info"), "No Gene Ontology (biological process) terms found, or GO.db is not installed in this deployment.")
  div(class = "card", div(class = "card-title", icon("circle-nodes"), "Gene Ontology (Biological Process)"),
      p(class = "submodule-desc", "Functional annotation via GO.db. Only Biological Process terms are curated for a single gene here (Molecular Function / Cellular Component are shown for gene panels - see the Gene Panel Enrichment section)."),
      go_body)
}

## ---- KEGG ------------------------------------------------------------------
tbc_section_kegg <- function(ext, kegg_map = NULL) {
  if (is.null(ext)) return(div(class = "card", div(class = "card-title", icon("diagram-project"), "KEGG Pathways"),
                                div(class = "empty-note", icon("circle-info"), "Not yet looked up - click Run below.")))
  kegg_body <- if (isTRUE(ext$kegg$ok) && nrow(ext$kegg$pathways) > 0) {
    kdf <- ext$kegg$pathways
    kdf$Link <- sprintf('<a href="https://www.kegg.jp/pathway/%s" target="_blank" rel="noopener">%s</a>', kdf$id, kdf$id)
    DT::datatable(kdf[, c("Link", "name")], colnames = c("KEGG ID", "Pathway"), rownames = FALSE, escape = FALSE,
                  options = list(dom = "t", paging = FALSE), class = "stripe hover compact")
  } else div(class = "empty-note", icon("circle-info"), if (isTRUE(ext$kegg$ok)) "No KEGG pathways found for this gene." else (ext$kegg$reason %||% "KEGG lookup unavailable."))
  kegg_map_body <- if (is.null(kegg_map)) NULL else if (isTRUE(kegg_map$ok) && !is.null(kegg_map$path) && file.exists(kegg_map$path) && requireNamespace("base64enc", quietly = TRUE)) {
    tagList(
      tags$div(style = "margin-top:10px;", tags$b(sprintf("Pathway diagram - %s (real pathview render of the actual KEGG map, this gene colored by its own log2FC)", kegg_map$pathway_name %||% kegg_map$pathway_id %||% ""))),
      tags$img(src = base64enc::dataURI(file = kegg_map$path, mime = "image/png"), style = "max-width:100%; height:auto; border:1px solid #eee; border-radius:6px; background:#fff;")
    )
  } else div(class = "empty-note", icon("circle-info"), kegg_map$error %||% "Pathway diagram unavailable.")
  div(class = "card", div(class = "card-title", icon("diagram-project"), "KEGG Pathways"), kegg_body, kegg_map_body)
}

## ---- Reactome ---------------------------------------------------------------
tbc_section_reactome <- function(ext, reactome_map = NULL) {
  if (is.null(ext)) return(div(class = "card", div(class = "card-title", icon("route"), "Reactome Pathways"),
                                div(class = "empty-note", icon("circle-info"), "Not yet looked up - click Run below.")))
  reactome_body <- if (isTRUE(ext$reactome$ok) && nrow(ext$reactome$pathways) > 0) {
    rdf <- ext$reactome$pathways
    rdf$Link <- sprintf('<a href="https://reactome.org/PathwayBrowser/#/%s" target="_blank" rel="noopener">%s</a>', rdf$stId, rdf$stId)
    DT::datatable(rdf[, c("Link", "displayName")], colnames = c("Reactome ID", "Pathway"), rownames = FALSE, escape = FALSE,
                  options = list(dom = "t", paging = FALSE), class = "stripe hover compact")
  } else div(class = "empty-note", icon("circle-info"), if (isTRUE(ext$reactome$ok)) "No Reactome pathways found for this gene." else (ext$reactome$reason %||% "Reactome lookup unavailable."))
  reactome_map_body <- if (is.null(reactome_map)) NULL else if (isTRUE(reactome_map$ok) && !is.null(reactome_map$path) && file.exists(reactome_map$path) && requireNamespace("base64enc", quietly = TRUE)) {
    tagList(
      tags$div(style = "margin-top:10px;", tags$b(sprintf("Pathway diagram - %s (rendered server-side by Reactome itself)", reactome_map$pathway_name %||% reactome_map$pathway_id %||% ""))),
      tags$img(src = base64enc::dataURI(file = reactome_map$path, mime = "image/png"), style = "max-width:100%; height:auto; border:1px solid #eee; border-radius:6px; background:#fff;")
    )
  } else div(class = "empty-note", icon("circle-info"), reactome_map$error %||% "Pathway diagram unavailable.")
  div(class = "card", div(class = "card-title", icon("route"), "Reactome Pathways"), reactome_body, reactome_map_body)
}

## ---- WikiPathways -----------------------------------------------------------
tbc_section_wikipathways <- function(ext) {
  if (is.null(ext)) return(div(class = "card", div(class = "card-title", icon("map"), "WikiPathways"),
                                div(class = "empty-note", icon("circle-info"), "Not yet looked up - click Run below.")))
  wp <- ext$wikipathways
  body <- if (!isTRUE(wp$ok)) div(class = "empty-note", icon("circle-info"), wp$reason %||% "WikiPathways lookup unavailable.")
          else if (is.null(wp$pathways) || nrow(wp$pathways) == 0) div(class = "empty-note", icon("circle-info"), "No WikiPathways pathways found for this gene.")
          else {
            wdf <- wp$pathways
            wdf$Link <- sprintf('<a href="https://www.wikipathways.org/instance/%s" target="_blank" rel="noopener">%s</a>', wdf$id, wdf$id)
            DT::datatable(wdf[, c("Link", "name")], colnames = c("WikiPathways ID", "Pathway"), rownames = FALSE, escape = FALSE,
                          options = list(dom = "t", paging = FALSE), class = "stripe hover compact")
          }
  div(class = "card", div(class = "card-title", icon("map"), "WikiPathways"),
      p(class = "submodule-desc", "Community-curated pathway diagrams, via the same gene-set source (msigdbr) the Multi-Omics Pathways module uses."),
      body)
}

## ---- Literature (PubMed) -----------------------------------------------------
tbc_section_literature <- function(lit, gene = NULL, query_used = NULL) {
  if (is.null(lit)) return(div(class = "card", div(class = "card-title", icon("book-open"), "Literature (PubMed)"),
                                div(class = "empty-note", icon("circle-info"), "Not yet searched - pick a query below and click Run.")))
  prov <- tbc_db_provenance("eutils.ncbi.nlm.nih.gov", sprintf("https://pubmed.ncbi.nlm.nih.gov/?term=%s", utils::URLencode(query_used %||% gene %||% "", reserved = TRUE)), "on PubMed")
  body <- if (!isTRUE(lit$ok)) div(class = "empty-note", icon("circle-info"), lit$reason %||% "PubMed lookup unavailable.")
          else if (is.null(lit$papers)) div(class = "empty-note", icon("circle-info"), "No matching PubMed records were returned for this query.")
          else {
            df <- lit$papers
            df$Link <- sprintf('<a href="https://pubmed.ncbi.nlm.nih.gov/%s/" target="_blank" rel="noopener">%s</a>', df$PMID, df$PMID)
            DT::datatable(df[, c("Title", "Authors", "Journal", "Year", "Link")], colnames = c("Title", "Authors", "Journal", "Year", "PMID"),
                          rownames = FALSE, escape = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")
          }
  div(class = "card", div(class = "card-title", icon("book-open"), "Literature (PubMed)"),
      prov,
      if (!is.null(query_used)) p(class = "submodule-desc", sprintf('Query: "%s"', query_used)) else NULL,
      body)
}

## Combined view for the downloadable report only (screen UI uses the standalone
## sections above, one per database button, per spec's organized-buttons request).
tbc_section_pathways <- function(ext, kegg_map = NULL, reactome_map = NULL) {
  tagList(tbc_section_go(ext), tbc_section_kegg(ext, kegg_map), tbc_section_reactome(ext, reactome_map), tbc_section_wikipathways(ext))
}

## ---- Genetic association + tractability (Open Targets) -------------------
tbc_section_genetics <- function(gen, gene = NULL, ensembl = NULL) {
  if (is.null(gen)) return(NULL)
  prov <- tbc_db_provenance("api.platform.opentargets.org", sprintf("https://platform.opentargets.org/target/%s", ensembl %||% ""), "on Open Targets")
  if (!isTRUE(gen$ok)) {
    return(div(class = "card", div(class = "card-title", icon("dna"), "Genetic Association & Tractability (Open Targets)"), prov,
               div(class = "empty-note", icon("circle-info"), gen$reason %||% "Open Targets lookup unavailable.")))
  }
  tract_df <- data.frame(Modality = unname(TBC_OT_MODALITY_LABELS[names(gen$tractability)]),
                          Status = unlist(gen$tractability), stringsAsFactors = FALSE)
  disease_body <- if (!is.null(gen$diseases) && nrow(gen$diseases) > 0) {
    DT::datatable(gen$diseases, rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")
  } else div(class = "empty-note", icon("circle-info"), "No disease associations found in Open Targets for this gene.")
  div(class = "card",
      div(class = "card-title", icon("dna"), "Genetic Association & Tractability (Open Targets)"),
      prov,
      p(class = "submodule-desc", sprintf("Aggregates GWAS Catalog + other genetic evidence across %s associated disease(s).", format(gen$n_diseases, big.mark = ","))),
      tags$b("Top disease associations"), disease_body,
      tags$div(style = "margin-top:10px;", tags$b("Druggability by modality")),
      DT::datatable(tract_df, rownames = FALSE, options = list(dom = "t", paging = FALSE), class = "stripe hover compact")
  )
}

## ---- Baseline tissue / blood-lineage expression (Human Protein Atlas) ----
tbc_section_tissue <- function(hpa, gene = NULL, ensembl = NULL) {
  if (is.null(hpa)) return(NULL)
  prov <- tbc_db_provenance("www.proteinatlas.org", sprintf("https://www.proteinatlas.org/%s-%s", ensembl %||% "", gene %||% ""), "on Human Protein Atlas")
  if (!isTRUE(hpa$ok)) {
    return(div(class = "card", div(class = "card-title", icon("layer-group"), "Baseline Tissue Expression (Human Protein Atlas)"), prov,
               div(class = "empty-note", icon("circle-info"), hpa$reason %||% "Human Protein Atlas lookup unavailable.")))
  }
  pairs <- list(
    "Tissue specificity" = hpa$tissue_specificity,
    "Blood lineage specificity" = hpa$blood_specificity,
    "Blood expression cluster" = hpa$blood_cluster,
    "Secretome" = hpa$secretome,
    "Protein class" = hpa$protein_class
  )
  div(class = "card",
      div(class = "card-title", icon("layer-group"), "Baseline Tissue Expression (Human Protein Atlas)"),
      prov,
      p(class = "submodule-desc", "Sanity-check: does this gene's baseline expression pattern plausibly match your loaded cohort's tissue/cell type?"),
      tbc_kv_table(pairs),
      if (!is.null(hpa$tissue_top)) tagList(tags$div(style = "margin-top:10px;", tags$b("Top tissues (nTPM)")),
        DT::datatable(hpa$tissue_top, rownames = FALSE, options = list(dom = "t", paging = FALSE), class = "stripe hover compact")) else NULL,
      if (!is.null(hpa$blood_top)) tagList(tags$div(style = "margin-top:10px;", tags$b("Top blood lineages (nTPM)")),
        DT::datatable(hpa$blood_top, rownames = FALSE, options = list(dom = "t", paging = FALSE), class = "stripe hover compact")) else NULL
  )
}

## ---- PPI neighborhood (STRING) - includes the actual diagram STRING
## itself renders, not just the parsed partner table. -----------------------
tbc_section_string <- function(net, image_path = NULL, gene = NULL) {
  if (is.null(net)) return(NULL)
  prov <- tbc_db_provenance("string-db.org", sprintf("https://string-db.org/cgi/network?identifiers=%s&species=9606", gene %||% ""), "on STRING")
  img_body <- if (!is.null(image_path) && file.exists(image_path) && requireNamespace("base64enc", quietly = TRUE)) {
    tags$img(src = base64enc::dataURI(file = image_path, mime = "image/png"), style = "max-width:100%; height:auto; border:1px solid #eee; border-radius:6px; background:#fff;")
  } else div(class = "empty-note", icon("circle-info"), "Network diagram unavailable.")
  net_body <- if (!isTRUE(net$ok)) div(class = "empty-note", icon("circle-info"), net$reason %||% "STRING lookup unavailable.")
              else if (is.null(net$partners)) div(class = "empty-note", icon("circle-info"), "No STRING interaction partners found for this gene.")
              else DT::datatable(net$partners, rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")
  div(class = "card",
      div(class = "card-title", icon("share-nodes"), "Protein Network (STRING)"),
      prov,
      tags$b("Network diagram - rendered server-side by STRING itself"), img_body,
      tags$div(style = "margin-top:10px;", tags$b("Interaction partners (all evidence channels)")), net_body
  )
}

## ---- Drug-gene interactions (DGIdb) --------------------------------------
tbc_section_dgidb <- function(drg, gene = NULL) {
  if (is.null(drg)) return(NULL)
  ## DGIdb is a client-rendered SPA; /genes/<symbol> looks valid but 404s client-side,
  ## so /results with these query params is the confirmed-working gene page URL.
  prov <- tbc_db_provenance("dgidb.org", sprintf("https://dgidb.org/results?searchType=gene&searchTerms=%s", gene %||% ""), "on DGIdb")
  drg_body <- if (!isTRUE(drg$ok)) div(class = "empty-note", icon("circle-info"), drg$reason %||% "DGIdb lookup unavailable.")
              else if (is.null(drg$drugs)) div(class = "empty-note", icon("circle-info"), "No known drug interactions found in DGIdb for this gene.")
              else DT::datatable(drg$drugs, rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact")
  div(class = "card",
      div(class = "card-title", icon("pills"), "Drug-Gene Interactions (DGIdb)"),
      prov,
      drg_body
  )
}

## ---- Evidence Summary / Database Coverage (spec: honest per-database ------
## status, never a checkmark just because a database exists) ----------------
tbc_evidence_status <- function(res, result_field) {
  if (is.null(res)) return("Not yet run")
  if (!isTRUE(res$ok)) return("Failed")
  val <- res[[result_field]]
  has_rows <- !is.null(val) && (!is.data.frame(val) || nrow(val) > 0)
  if (isTRUE(has_rows)) "Results found" else "No results"
}

TBC_EVIDENCE_DBS <- list(
  list(key = "kegg", label = "KEGG", field = "pathways"),
  list(key = "reactome", label = "Reactome", field = "pathways"),
  list(key = "wikipathways", label = "WikiPathways", field = "pathways"),
  list(key = "genetics", label = "Disease / genetic association (Open Targets)", field = "diseases"),
  list(key = "drugs", label = "Drug / target (DGIdb)", field = "drugs"),
  list(key = "network", label = "Protein interaction (STRING)", field = "partners"),
  list(key = "tissue", label = "Expression (Human Protein Atlas)", field = "tissue_top"),
  list(key = "literature", label = "Literature (PubMed)", field = "papers")
)

## ---- Database Comparison (spec: compare independent resources side by side) --
tbc_section_db_comparison <- function(ext) {
  if (is.null(ext)) return(NULL)
  rows <- list()
  add <- function(resource, evidence_type, results_n, significance, source) {
    rows[[length(rows) + 1]] <<- data.frame(Resource = resource, `Evidence type` = evidence_type, Results = results_n,
                                             Significance = significance, Source = source, check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (!is.null(ext$go)) add("Gene Ontology (BP)", "Functional annotation", nrow(ext$go), "Not applicable (annotation, not enrichment)", "GO.db")
  if (isTRUE(ext$kegg$ok)) add("KEGG", "Curated pathway membership", nrow(ext$kegg$pathways), "Not applicable (membership, not enrichment)", "KEGGREST")
  if (isTRUE(ext$reactome$ok)) add("Reactome", "Curated pathway membership", nrow(ext$reactome$pathways), "Not applicable (membership, not enrichment)", "Reactome ContentService")
  if (isTRUE(ext$wikipathways$ok)) add("WikiPathways", "Curated pathway membership", nrow(ext$wikipathways$pathways), "Not applicable (membership, not enrichment)", "msigdbr")
  if (isTRUE(ext$genetics$ok)) add("Open Targets", "Genetic association", ext$genetics$n_diseases %||% 0L, "Association score (not a p-value)", "Open Targets Platform")
  if (isTRUE(ext$drugs$ok)) add("DGIdb", "Drug-gene interaction", if (!is.null(ext$drugs$drugs)) nrow(ext$drugs$drugs) else 0L, "Interaction score (not a p-value)", "DGIdb")
  if (isTRUE(ext$network$ok)) add("STRING", "Protein-protein interaction", if (!is.null(ext$network$partners)) nrow(ext$network$partners) else 0L, "Combined confidence score", "STRING")
  if (isTRUE(ext$tissue$ok)) add("Human Protein Atlas", "Baseline tissue/blood expression", if (!is.null(ext$tissue$tissue_top)) nrow(ext$tissue$tissue_top) else 0L, "Not applicable", "proteinatlas.org")
  if (isTRUE(ext$literature$ok)) add("PubMed", "Literature", if (!is.null(ext$literature$papers)) nrow(ext$literature$papers) else 0L, "Not applicable", "NCBI E-utilities")
  if (length(rows) == 0) {
    return(div(class = "card", div(class = "card-title", icon("scale-balanced"), "Database Comparison"),
               div(class = "empty-note", icon("circle-info"), "Run at least one database above to populate this comparison.")))
  }
  df <- do.call(rbind, rows)
  div(class = "card", div(class = "card-title", icon("scale-balanced"), "Database Comparison"),
      p(class = "submodule-desc", "Side-by-side view of what each already-run database returned, so you can judge whether independent resources converge on the same biology."),
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact"))
}

tbc_report_css <- function() {
  "body{font-family:-apple-system,Helvetica,Arial,sans-serif; max-width:900px; margin:24px auto; color:#222;}
   .card{border:1px solid #ddd; border-radius:10px; padding:14px 18px; margin-bottom:16px;}
   .card-title{font-weight:700; font-size:15px; margin-bottom:8px;}
   .submodule-desc{color:#666; font-size:12.5px;}
   .empty-note{background:#f6f6f6; border-left:3px solid #999; padding:8px 12px; border-radius:4px; font-size:13px;}
   table{border-collapse:collapse; width:100%;} td,th{border:1px solid #eee; padding:4px 8px; font-size:13px; text-align:left;}
   .pipeline-status-strip{display:flex; flex-wrap:wrap; gap:8px;}
   .pipeline-status-chip{display:inline-flex; align-items:center; gap:6px; font-size:12.5px; font-weight:600; padding:5px 12px; border-radius:999px; border:1px solid transparent;}
   .pipeline-status-chip.status-done{background:#E9F7EF; color:#1E7A3C; border-color:#C3E6CD;}
   .pipeline-status-chip.status-pending{background:#FDF3E3; color:#8A5A00; border-color:#F2DDB0;}
   .pipeline-status-chip.status-neutral{background:#EEF1F4; color:#5A6472; border-color:#DDE2E7;}"
}

## Mirrors the on-screen 9-tab structure exactly (spec §11: "downloads contain
## the actual displayed evidence") - same section-builder functions, plots
## flattened to static images since a downloaded HTML file has no live Shiny outputs.
tbc_build_report_tags <- function(d, dataset, ext = NULL) {
  sgd <- d$single_gene_diag; sgcv <- d$single_gene_cv
  dist_plot <- if (isTRUE(d$live$ok) && isTRUE(d$live$overall$ok))
    tryCatch(tbc_plot_expression_dist(data.frame(expr = d$live$values, group = d$live$group_vec), "Expression (analysis scale)"), error = function(e) NULL) else NULL
  volcano_plot <- if (!is.null(d$selected_run_table))
    tryCatch(tbc_plot_volcano_highlight(d$selected_run_table, d$gene), error = function(e) NULL) else NULL
  train_roc_plot <- if (isTRUE(sgd$ok)) tryCatch(tbc_plot_roc(sgd$roc_obj, "Single-gene ROC (Training)", sgd$auc, sgd$ci_lo, sgd$ci_hi), error = function(e) NULL) else NULL
  train_pr_plot <- if (isTRUE(sgd$ok)) tryCatch(tbc_plot_pr(tbc_pr_from_roc(sgd$roc_obj)), error = function(e) NULL) else NULL
  internal_roc_plot <- if (isTRUE(sgcv$ok)) tryCatch(tbc_plot_roc(sgcv$roc_obj, "Single-gene ROC (Internal Validation, pooled CV)", sgcv$auc), error = function(e) NULL) else NULL

  img_tag <- function(p) {
    uri <- tbc_ggsave_datauri(p)
    if (is.null(uri)) div(class = "empty-note", "Plot unavailable.") else tags$img(src = uri, style = "max-width:100%; height:auto;")
  }
  img_widget <- function(p) if (is.null(p)) NULL else img_tag(p)

  tagList(
    tags$h2(sprintf("Transcriptomic Biomarker Card: %s", d$gene)),
    tags$h3("1. Gene description"), tbc_section_identity(d),
    tags$h3("2. Biomarker Status"), tbc_section_evidence_glance(d, ext, sgd, sgcv),
    tags$h3("3. Dataset"), tbc_section_dataset_cohort(dataset, d$live),
    tags$h3("4. Differential Expression"), tbc_section_differential_expression(d),
    if (!is.null(volcano_plot)) div(class = "card", div(class = "card-title", "Differential Expression (Volcano)"), img_tag(volcano_plot)) else NULL,
    if (!is.null(dist_plot)) div(class = "card", div(class = "card-title", "Expression Distribution (Your Dataset)"), img_tag(dist_plot)) else NULL,
    tags$h3("5. Single-Gene vs Multi-Gene Signature"), tbc_section_signature_comparison(d, sgd, img_widget(train_roc_plot)),
    tags$h3("6. Biomarker Performance"), tbc_section_biomarker_performance(d, sgd, sgcv, img_widget(train_roc_plot), img_widget(train_pr_plot), img_widget(internal_roc_plot)),
    tags$h3("7. External Databases"), tbc_external_db_banner(),
    tbc_section_genetics(ext$genetics, d$gene, d$gene_identity$ensembl),
    tbc_section_tissue(ext$tissue, d$gene, d$gene_identity$ensembl),
    ## Network image not re-fetched for the static report (would require a
    ## fresh network call at download time); the live card shows it.
    tbc_section_string(ext$network, NULL, d$gene),
    tbc_section_dgidb(ext$drugs, d$gene),
    tbc_section_pathways(ext),
    tbc_section_literature(ext$literature, d$gene, ext$literature_query),
    tbc_section_db_comparison(ext),
    tags$p(style = "color:#888; font-size:12px; margin-top:16px;",
           if (is.null(ext)) "External database evidence (Open Targets/HPA/STRING/DGIdb/GO/KEGG/Reactome) not looked up before this report was generated."
           else "Multi-gene panel sensitivity/specificity/full ROC curves: see the Diagnostic Classifier tab's Result panel (not stored in this session's shared results).")
  )
}

## =============================================================================
## Gene Panel mode: identity resolution, panel-level enrichment/network/
## convergence, and the section builders + report that present them.
## =============================================================================

## ---- Identifier resolution (spec: never silently drop an unresolved gene) --
## Reuses the app's shared harmonizer (cx_harmonize_gene_ids(),
## R/crossomics/crossomics_integration_helpers.R) - the same one
## mod_enrichment.R's own gene-list flow already uses.
tbc_panel_identity <- function(genes_raw) {
  genes_raw <- unique(trimws(as.character(genes_raw)))
  genes_raw <- genes_raw[nzchar(genes_raw)]
  if (length(genes_raw) == 0) return(list(ok = FALSE, reason = "No gene identifiers were provided.", df = NULL))
  harm <- cx_harmonize_gene_ids(genes_raw)
  if (!isTRUE(harm$ok)) return(list(ok = FALSE, reason = harm$error %||% "Gene identifier harmonization is unavailable in this deployment.", df = NULL))
  df <- harm$df
  df$resolved <- df$match_type %in% c("exact_symbol", "exact_entrez", "exact_ensembl", "alias_resolved")
  status_labels <- c(exact_symbol = "Resolved (exact gene symbol)", exact_entrez = "Resolved (NCBI Entrez ID)",
                      exact_ensembl = "Resolved (Ensembl Gene ID)", alias_resolved = "Resolved (via known alias)",
                      ambiguous = "Ambiguous - multiple candidate symbols, not guessed", unmatched = "Unresolved - identifier not recognized")
  df$status_label <- unname(status_labels[df$match_type])
  list(ok = TRUE, reason = NULL, df = df, n_submitted = length(genes_raw), n_resolved = sum(df$resolved), n_unresolved = sum(!df$resolved))
}

## ---- Per-gene status rows, reusing the exact same single-gene helpers ------
## this file already uses for the single-gene card (no recomputation logic
## duplicated - only looped across the panel).
tbc_panel_gene_rows <- function(symbols, dataset, results) {
  if (length(symbols) == 0) return(NULL)
  do.call(rbind, lapply(symbols, function(g) {
    in_ds <- tryCatch(g %in% rownames(dataset$expr), error = function(e) FALSE)
    dge_hits <- tbc_dge_matches(g, results)
    cand <- tbc_candidate_status(g, results)
    sig <- tbc_signature_membership(g, results)
    diag <- tbc_diagnostic_lookup(g, results)
    dge_best_fdr <- if (!is.null(dge_hits) && nrow(dge_hits) > 0) min(dge_hits$adj.P.Val, na.rm = TRUE) else NA_real_
    data.frame(
      Gene = g,
      `In loaded dataset` = if (in_ds) "Yes" else "No",
      `Best saved DGE FDR` = if (!is.na(dge_best_fdr)) signif(dge_best_fdr, 3) else NA,
      `Candidate gene` = if (isTRUE(cand$any)) "Yes" else "No",
      `In a consensus signature` = if (length(sig) > 0 && any(vapply(sig, function(x) isTRUE(x$in_signature), logical(1)))) "Yes" else "No",
      `In an evaluated classifier panel` = if (length(diag) > 0) "Yes" else "No",
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }))
}

## ---- Panel-level enrichment: GO(BP/MF/CC)/KEGG/Reactome/WikiPathways -------
## Reuses the Multi-Omics Pathways module's own ORA runners
## (mp_run_ora_go/kegg/reactome/wikipathways(), multiomics_pathway_helpers.R)
## rather than re-implementing clusterProfiler/ReactomePA calls a second time.
tbc_panel_enrichment <- function(genes_entrez, universe_entrez) {
  list(
    go_bp = mp_run_ora_go(genes_entrez, universe_entrez, "BP"),
    go_mf = mp_run_ora_go(genes_entrez, universe_entrez, "MF"),
    go_cc = mp_run_ora_go(genes_entrez, universe_entrez, "CC"),
    kegg = mp_run_ora_kegg(genes_entrez, universe_entrez),
    reactome = mp_run_ora_reactome(genes_entrez, universe_entrez),
    wikipathways = mp_run_ora_wikipathways(genes_entrez, universe_entrez)
  )
}

## ---- Panel network/hubs: reuses mod_enrichment.R's own top-level helpers ---
## (fetch_string_degrees/fetch_string_network_png/build_wgcna_hub_lookup -
## globally callable functions defined outside its moduleServer) unmodified.
tbc_panel_network <- function(genes) {
  degrees <- tryCatch(fetch_string_degrees(genes), error = function(e) list())
  hub_lookup <- tryCatch(build_wgcna_hub_lookup(), error = function(e) NULL)
  png <- tryCatch(fetch_string_network_png(genes), error = function(e) NULL)
  hub_table <- do.call(rbind, lapply(genes, function(g) {
    sd <- degrees[[g]]
    w <- if (!is.null(hub_lookup)) hub_lookup(g) else list(module = NA_character_, disease_module = FALSE, kme = NA_real_, is_hub = FALSE, connectivity = NA_real_)
    data.frame(Gene = g, `STRING partners (this panel)` = if (!is.null(sd)) sd$degree else 0L,
               `Partner list` = if (!is.null(sd) && sd$degree > 0) sd$partners else "",
               `WGCNA module` = w$module %||% NA_character_, `WGCNA hub gene` = isTRUE(w$is_hub), `WGCNA kME` = round(w$kme %||% NA_real_, 3),
               check.names = FALSE, stringsAsFactors = FALSE)
  }))
  if (!is.null(hub_table)) hub_table <- hub_table[order(-hub_table$`STRING partners (this panel)`), , drop = FALSE]
  list(ok = TRUE, hub_table = hub_table, image_path = png)
}

## ---- Panel disease/drug convergence: loops the existing single-gene --------
## Open Targets / DGIdb calls (unmodified) and aggregates locally - no new
## external client. The aggregation step itself is a separate pure function
## (tbc_aggregate_convergence) so it's unit-testable without a live network call.
tbc_aggregate_convergence <- function(long_df, item_col) {
  if (is.null(long_df) || nrow(long_df) == 0) return(NULL)
  ## Gene ~ <item_col>: groups by the item (e.g. Disease/Drug) and collapses
  ## the Gene column - one row per item, listing which genes share it.
  agg <- stats::aggregate(stats::reformulate(item_col, response = "Gene"), long_df, function(x) paste(sort(unique(x)), collapse = ", "))
  agg$`Gene count` <- vapply(strsplit(agg$Gene, ", "), length, integer(1))
  agg <- agg[order(-agg$`Gene count`), , drop = FALSE]
  colnames(agg)[colnames(agg) == "Gene"] <- "Genes"
  agg
}

tbc_panel_disease_convergence <- function(gene_ensembl_map) {
  per_gene <- lapply(names(gene_ensembl_map), function(g) {
    r <- tbc_opentargets_evidence_for_gene(gene_ensembl_map[[g]], top_n_diseases = 10)
    if (!isTRUE(r$ok) || is.null(r$diseases) || nrow(r$diseases) == 0) return(NULL)
    data.frame(Gene = g, Disease = r$diseases$Disease, stringsAsFactors = FALSE)
  })
  per_gene <- Filter(Negate(is.null), per_gene)
  if (length(per_gene) == 0) return(list(ok = FALSE, reason = "No Open Targets disease associations were found for any resolved gene in this panel.", table = NULL))
  list(ok = TRUE, reason = NULL, table = tbc_aggregate_convergence(do.call(rbind, per_gene), "Disease"))
}

tbc_panel_drug_convergence <- function(genes) {
  per_gene <- lapply(genes, function(g) {
    r <- tbc_dgidb_drugs_for_gene(g, top_n = 20)
    if (!isTRUE(r$ok) || is.null(r$drugs) || nrow(r$drugs) == 0) return(NULL)
    data.frame(Gene = g, Drug = r$drugs$Drug, stringsAsFactors = FALSE)
  })
  per_gene <- Filter(Negate(is.null), per_gene)
  if (length(per_gene) == 0) return(list(ok = FALSE, reason = "No DGIdb drug interactions were found for any resolved gene in this panel.", table = NULL))
  list(ok = TRUE, reason = NULL, table = tbc_aggregate_convergence(do.call(rbind, per_gene), "Drug"))
}

## ---- Panel section builders -------------------------------------------------

tbc_section_panel_identity <- function(pid) {
  if (is.null(pid) || !isTRUE(pid$ok)) {
    return(div(class = "card", div(class = "card-title", icon("list-check"), "Identifier Resolution"),
               div(class = "empty-note", icon("triangle-exclamation"), pid$reason %||% "No genes to resolve.")))
  }
  df <- pid$df[, c("input_id", "status_label", "canonical_symbol", "entrez_id", "ensembl_id")]
  colnames(df) <- c("Submitted identifier", "Resolution status", "Canonical symbol", "NCBI Entrez ID", "Ensembl Gene ID")
  div(class = "card",
      div(class = "card-title", icon("list-check"), "Identifier Resolution"),
      p(class = "submodule-desc", sprintf("%d gene(s) submitted, %d resolved, %d unresolved. Every submitted identifier is listed below - none are silently dropped.", pid$n_submitted, pid$n_resolved, pid$n_unresolved)),
      DT::datatable(df, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
  )
}

tbc_section_panel_overview <- function(pid, gene_rows) {
  if (is.null(pid) || !isTRUE(pid$ok)) return(NULL)
  pairs <- list("Genes submitted" = pid$n_submitted, "Genes resolved" = pid$n_resolved, "Genes unresolved" = pid$n_unresolved)
  div(class = "card",
      div(class = "card-title", icon("layer-group"), "Panel Overview"),
      tbc_kv_table(pairs),
      tags$div(style = "margin-top:10px;", tags$b("Per-gene status (this app's own session results)")),
      if (is.null(gene_rows) || nrow(gene_rows) == 0) div(class = "empty-note", icon("circle-info"), "No resolved genes to summarize.")
      else DT::datatable(gene_rows, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
  )
}

tbc_section_panel_enrichment <- function(genrich, universe_label = NULL, n_tested = NULL) {
  if (is.null(genrich)) return(div(class = "card", div(class = "card-title", icon("diagram-project"), "Gene Panel Enrichment"),
                                    div(class = "empty-note", icon("circle-info"), "Not yet run - click Run below (requires at least 3 resolved genes measured in the loaded dataset).")))
  block <- function(res, label) {
    body <- if (!isTRUE(res$ok)) div(class = "empty-note", icon("circle-info"), res$error %||% "No result.")
            else DT::datatable(
              data.frame(ID = res$df$ID, Description = res$df$Description, `Gene ratio` = res$df$GeneRatio, Count = res$df$Count,
                         `P-value` = signif(res$df$pvalue, 3), FDR = signif(res$df$p.adjust, 3), `Matched genes` = res$df$geneID,
                         check.names = FALSE, stringsAsFactors = FALSE)[order(res$df$p.adjust), , drop = FALSE],
              rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE), class = "stripe hover compact")
    tagList(tags$div(style = "margin-top:14px;", tags$b(label)), body)
  }
  n_sig <- sum(vapply(genrich, function(r) isTRUE(r$ok), logical(1)))
  div(class = "card",
      div(class = "card-title", icon("diagram-project"), "Gene Panel Enrichment"),
      p(class = "submodule-desc", sprintf("Over-representation analysis (hypergeometric test, Benjamini-Hochberg FDR) against %s.%s",
                                           universe_label %||% "the currently loaded dataset as background",
                                           if (!is.null(n_tested)) sprintf(" %d gene(s) tested.", n_tested) else "")),
      if (n_sig == 0) div(class = "empty-note", icon("circle-info"), "None of GO/KEGG/Reactome/WikiPathways returned enriched terms for this panel against this background.") else NULL,
      block(genrich$go_bp, "Gene Ontology - Biological Process"),
      block(genrich$go_mf, "Gene Ontology - Molecular Function"),
      block(genrich$go_cc, "Gene Ontology - Cellular Component"),
      block(genrich$kegg, "KEGG"),
      block(genrich$reactome, "Reactome"),
      block(genrich$wikipathways, "WikiPathways")
  )
}

tbc_section_panel_network <- function(net, image_path = NULL) {
  if (is.null(net)) return(div(class = "card", div(class = "card-title", icon("share-nodes"), "Network & Hub Genes"),
                                div(class = "empty-note", icon("circle-info"), "Not yet run - click Run below.")))
  img_body <- if (!is.null(image_path) && file.exists(image_path) && requireNamespace("base64enc", quietly = TRUE)) {
    tags$img(src = base64enc::dataURI(file = image_path, mime = "image/png"), style = "max-width:100%; height:auto; border:1px solid #eee; border-radius:6px; background:#fff;")
  } else div(class = "empty-note", icon("circle-info"), "Network diagram unavailable (STRING lookup failed or too few resolvable genes).")
  hub_body <- if (is.null(net$hub_table) || nrow(net$hub_table) == 0) div(class = "empty-note", icon("circle-info"), "No hub information available.")
              else DT::datatable(net$hub_table, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
  div(class = "card",
      div(class = "card-title", icon("share-nodes"), "Network & Hub Genes"),
      p(class = "submodule-desc", "How the panel's genes interact with each other: live STRING PPI degree within just this list, plus this project's own precomputed WGCNA co-expression hub status."),
      tags$b("STRING network diagram (this panel only)"), img_body,
      tags$div(style = "margin-top:10px;", tags$b("Hub genes (STRING PPI degree + WGCNA hub status)")), hub_body
  )
}

tbc_section_panel_disease <- function(dc) {
  if (is.null(dc)) return(div(class = "card", div(class = "card-title", icon("dna"), "Disease Convergence (Open Targets)"),
                               div(class = "empty-note", icon("circle-info"), "Not yet run - click Run below.")))
  body <- if (!isTRUE(dc$ok)) div(class = "empty-note", icon("circle-info"), dc$reason %||% "Unavailable.")
          else DT::datatable(dc$table, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
  div(class = "card", div(class = "card-title", icon("dna"), "Disease Convergence (Open Targets)"),
      p(class = "submodule-desc", "Diseases associated with more than one panel gene in Open Targets - convergence across genes, not proof of a shared causal mechanism."),
      body)
}

tbc_section_panel_drug <- function(dc) {
  if (is.null(dc)) return(div(class = "card", div(class = "card-title", icon("pills"), "Drug/Target Convergence (DGIdb)"),
                               div(class = "empty-note", icon("circle-info"), "Not yet run - click Run below.")))
  body <- if (!isTRUE(dc$ok)) div(class = "empty-note", icon("circle-info"), dc$reason %||% "Unavailable.")
          else DT::datatable(dc$table, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact")
  div(class = "card", div(class = "card-title", icon("pills"), "Drug/Target Convergence (DGIdb)"),
      p(class = "submodule-desc", "Drugs targeting more than one panel gene in DGIdb."),
      body)
}

tbc_panel_evidence_status <- function(res, field) {
  if (is.null(res)) return("Not yet run")
  if (!isTRUE(res$ok)) return("Failed")
  val <- res[[field]]
  if (!is.null(val) && is.data.frame(val) && nrow(val) > 0) "Results found" else "No results"
}

tbc_section_panel_evidence_summary <- function(genrich, net, dc_disease, dc_drug) {
  rows <- list(
    list(label = "GO Biological Process", status = tbc_panel_evidence_status(genrich$go_bp, "df")),
    list(label = "GO Molecular Function", status = tbc_panel_evidence_status(genrich$go_mf, "df")),
    list(label = "GO Cellular Component", status = tbc_panel_evidence_status(genrich$go_cc, "df")),
    list(label = "KEGG", status = tbc_panel_evidence_status(genrich$kegg, "df")),
    list(label = "Reactome", status = tbc_panel_evidence_status(genrich$reactome, "df")),
    list(label = "WikiPathways", status = tbc_panel_evidence_status(genrich$wikipathways, "df")),
    list(label = "Network & Hub Genes", status = tbc_panel_evidence_status(net, "hub_table")),
    list(label = "Disease Convergence", status = tbc_panel_evidence_status(dc_disease, "table")),
    list(label = "Drug/Target Convergence", status = tbc_panel_evidence_status(dc_drug, "table"))
  )
  status_class <- c("Results found" = "status-done", "No results" = "status-neutral", "Failed" = "status-pending", "Not yet run" = "status-neutral")
  status_icon <- c("Results found" = "circle-check", "No results" = "circle-minus", "Failed" = "triangle-exclamation", "Not yet run" = "circle-info")
  chips <- lapply(rows, function(it) span(class = paste("pipeline-status-chip", status_class[[it$status]]), icon(status_icon[[it$status]]), sprintf("%s: %s", it$label, it$status)))
  div(class = "card",
      div(class = "card-title", icon("clipboard-check"), "Panel Evidence Summary"),
      p(class = "submodule-desc", "Actual retrieval status per panel-level database section."),
      div(class = "pipeline-status-strip", chips))
}

tbc_section_panel_comparison <- function(genrich, net, dc_disease, dc_drug) {
  rows <- list()
  add <- function(resource, evidence_type, results_n, best_fdr, source) {
    rows[[length(rows) + 1]] <<- data.frame(Resource = resource, `Evidence type` = evidence_type, `Terms/rows` = results_n,
                                             `Best FDR` = best_fdr, Source = source, check.names = FALSE, stringsAsFactors = FALSE)
  }
  add_enrich <- function(res, resource, source) {
    if (isTRUE(res$ok) && !is.null(res$df) && nrow(res$df) > 0) add(resource, "Over-representation (ORA)", nrow(res$df), signif(min(res$df$p.adjust, na.rm = TRUE), 3), source)
  }
  if (!is.null(genrich)) {
    add_enrich(genrich$go_bp, "GO Biological Process", "clusterProfiler::enrichGO")
    add_enrich(genrich$go_mf, "GO Molecular Function", "clusterProfiler::enrichGO")
    add_enrich(genrich$go_cc, "GO Cellular Component", "clusterProfiler::enrichGO")
    add_enrich(genrich$kegg, "KEGG", "clusterProfiler::enrichKEGG")
    add_enrich(genrich$reactome, "Reactome", "ReactomePA::enrichPathway")
    add_enrich(genrich$wikipathways, "WikiPathways", "clusterProfiler::enricher (msigdbr)")
  }
  if (!is.null(net) && !is.null(net$hub_table) && nrow(net$hub_table) > 0) add("STRING / WGCNA", "Protein interaction + co-expression hubs", nrow(net$hub_table), NA, "STRING + this project's WGCNA")
  if (!is.null(dc_disease) && isTRUE(dc_disease$ok)) add("Open Targets", "Disease convergence", nrow(dc_disease$table), NA, "Open Targets Platform")
  if (!is.null(dc_drug) && isTRUE(dc_drug$ok)) add("DGIdb", "Drug/target convergence", nrow(dc_drug$table), NA, "DGIdb")
  if (length(rows) == 0) {
    return(div(class = "card", div(class = "card-title", icon("scale-balanced"), "Database Comparison"),
               div(class = "empty-note", icon("circle-info"), "Run the panel enrichment/network/convergence sections above to populate this comparison.")))
  }
  df <- do.call(rbind, rows)
  div(class = "card", div(class = "card-title", icon("scale-balanced"), "Database Comparison"),
      p(class = "submodule-desc", "Which independent resources found evidence for this panel, and how strong."),
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "stripe hover compact"))
}

tbc_build_panel_report_tags <- function(pid, gene_rows, genrich, net, dc_disease, dc_drug, universe_label = NULL, n_tested = NULL) {
  tagList(
    tags$h2(sprintf("Gene Panel Biomarker Report (%d gene(s) submitted)", pid$n_submitted %||% 0)),
    tbc_section_panel_evidence_summary(genrich, net, dc_disease, dc_drug),
    tbc_section_panel_identity(pid),
    tbc_section_panel_overview(pid, gene_rows),
    tbc_section_panel_enrichment(genrich, universe_label, n_tested),
    tbc_section_panel_network(net, NULL),
    tbc_section_panel_disease(dc_disease),
    tbc_section_panel_drug(dc_drug),
    tbc_section_panel_comparison(genrich, net, dc_disease, dc_drug)
  )
}

## ---- Free-text gene list parsing (same token rule as mod_enrichment.R's
## own gene_list textarea: split on comma/newline/tab/space, dedupe, drop blanks).
tbc_split_gene_text <- function(text) {
  toks <- trimws(unlist(strsplit(text %||% "", "[,\n\t ]+")))
  unique(toks[nzchar(toks)])
}

## ---- Config / UI -----------------------------------------------------

mod_biomarkercard_config <- list(
  id = "biomarkercard", group = "Interpretation",
  title = "Biomarker Card",
  description = "Single-gene or gene-panel biomarker intelligence workspace: preloaded or uploaded data, live GO/KEGG/Reactome/WikiPathways/Open Targets/HPA/STRING/DGIdb/PubMed lookups, and a downloadable report.",
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
            radioButtons(ns("bmc_mode"), "Mode", inline = TRUE,
                         choices = c("Single gene" = "single", "Gene panel (multiple genes)" = "panel")),
            div(class = "card-title", icon("list-check"), "Upload Biomarker"),
            radioButtons(ns("bmc_search_mode"), NULL, inline = TRUE,
                         choices = c("Type gene symbol(s)" = "gene",
                                     "Browse Differential Expression results" = "dge",
                                     "Browse candidate genes" = "candidates",
                                     "Browse a feature-selected signature" = "signature",
                                     "Upload a gene list" = "upload")),
            conditionalPanel(condition = sprintf("input['%s'] == 'gene' && input['%s'] != 'panel'", ns("bmc_search_mode"), ns("bmc_mode")),
                              textInput(ns("bmc_gene_input"), "Gene symbol", placeholder = "TNF")),
            conditionalPanel(condition = sprintf("input['%s'] == 'gene' && input['%s'] == 'panel'", ns("bmc_search_mode"), ns("bmc_mode")),
                              textAreaInput(ns("bmc_panel_gene_input"), "Gene symbols / IDs (one per line, or comma/space separated)",
                                            placeholder = "TNF\nIL6\nSTAT3\n... (gene symbols, NCBI Entrez IDs, or Ensembl Gene IDs - any mix)", rows = 5)),
            conditionalPanel(condition = sprintf("input['%s'] == 'dge'", ns("bmc_search_mode")),
                              if (length(dge_runs) == 0) div(class = "empty-note", icon("circle-info"), "No Differential Expression run yet this session - run the Differential Expression tab first, or type a gene symbol directly.")
                              else tagList(
                                selectInput(ns("bmc_dge_run"), "Differential Expression run", choices = stats::setNames(names(dge_runs), vapply(dge_runs, function(r) r$contrast, character(1)))),
                                conditionalPanel(condition = sprintf("input['%s'] == 'panel'", ns("bmc_mode")),
                                                  actionButton(ns("bmc_dge_use_panel"), "Use all listed genes as panel", icon = icon("layer-group"), class = "btn-default btn-sm")),
                                uiOutput(ns("bmc_dge_results_ui"))
                              )),
            conditionalPanel(condition = sprintf("input['%s' ] == 'candidates'", ns("bmc_search_mode")),
                              if (is.null(cand)) div(class = "empty-note", icon("circle-info"), "Candidate Gene Identification has not been run this session - run that tab first, or type a gene symbol directly.")
                              else tagList(
                                radioButtons(ns("bmc_cand_set"), "Candidate set", inline = TRUE,
                                             choices = Filter(function(v) !is.null(cand[[v]]),
                                                              c("Female" = "female", "Male" = "male", "Final (combined)" = "final"))),
                                conditionalPanel(condition = sprintf("input['%s'] == 'panel'", ns("bmc_mode")),
                                                  actionButton(ns("bmc_cand_use_panel"), "Use all listed genes as panel", icon = icon("layer-group"), class = "btn-default btn-sm")),
                                uiOutput(ns("bmc_cand_results_ui"))
                              )),
            conditionalPanel(condition = sprintf("input['%s'] == 'signature'", ns("bmc_search_mode")),
                              if (is.null(fs)) div(class = "empty-note", icon("circle-info"), "ML Feature Selection has not been run this session - run that tab first, or type a gene symbol directly.")
                              else tagList(
                                radioButtons(ns("bmc_sig_sex"), "Stratum", inline = TRUE,
                                             choices = stats::setNames(intersect(c("female", "male", "pooled"), names(fs)), tools::toTitleCase(intersect(c("female", "male", "pooled"), names(fs))))),
                                conditionalPanel(condition = sprintf("input['%s'] == 'panel'", ns("bmc_mode")),
                                                  actionButton(ns("bmc_sig_use_panel"), "Use all listed genes as panel", icon = icon("layer-group"), class = "btn-default btn-sm")),
                                uiOutput(ns("bmc_sig_results_ui"))
                              )),
            conditionalPanel(condition = sprintf("input['%s'] == 'upload'", ns("bmc_search_mode")),
                              p(class = "submodule-desc", "Upload a gene-identifier list (one per line, or the first column of a CSV/TSV - gene symbols, NCBI Entrez IDs, or Ensembl Gene IDs, any mix), or a Diagnostic Classifier RDS export (that tab's own \"Save trained model\" download) - auto-detected by file extension. Every identifier is resolved and reported below; none are silently dropped."),
                              fileInput(ns("bmc_upload_file"), "Biomarker list (.csv, .txt, or .rds)", accept = c(".csv", ".txt", ".rds")),
                              actionButton(ns("bmc_upload_load_btn"), "Load Uploaded List", icon = icon("play"), class = "btn-sm"),
                              conditionalPanel(condition = sprintf("input['%s'] == 'panel'", ns("bmc_mode")),
                                                actionButton(ns("bmc_upload_use_panel"), "Use entire uploaded list as panel", icon = icon("layer-group"), class = "btn-default btn-sm")),
                              uiOutput(ns("bmc_upload_results_ui")))
        ),
        div(class = "empty-note", style = "display:flex; align-items:center; justify-content:space-between; gap:12px; flex-wrap:wrap;",
            uiOutput(ns("bmc_selection_status_ui"), inline = TRUE),
            conditionalPanel(condition = sprintf("input['%s'] != 'panel'", ns("bmc_mode")),
                              actionButton(ns("bmc_generate_btn"), "Generate Biomarker Card", icon = icon("id-card"), class = "btn-primary btn-sm")),
            conditionalPanel(condition = sprintf("input['%s'] == 'panel'", ns("bmc_mode")),
                              actionButton(ns("bmc_generate_panel_btn"), "Generate Gene Panel Report", icon = icon("layer-group"), class = "btn-primary btn-sm")))
      )
    })

    bmc_picked_gene <- reactiveVal(NULL)
    bmc_panel_genes <- reactiveVal(character(0))
    has_card <- reactiveVal(FALSE)
    has_panel_card <- reactiveVal(FALSE)
    observeEvent(input$bmc_search_mode, { bmc_picked_gene(NULL); bmc_panel_genes(character(0)) }, ignoreInit = TRUE)
    observeEvent(input$bmc_mode, { bmc_picked_gene(NULL); bmc_panel_genes(character(0)) }, ignoreInit = TRUE)

    output$bmc_selection_status_ui <- renderUI({
      pmode <- input$bmc_mode %||% "single"
      smode <- input$bmc_search_mode %||% "gene"
      if (identical(pmode, "panel")) {
        genes <- if (identical(smode, "gene")) tbc_split_gene_text(input$bmc_panel_gene_input) else bmc_panel_genes()
        if (length(genes) > 0) {
          tagList(icon("circle-check", style = "color:#0ca30c;"), tags$b(sprintf("%d gene(s) queued for the panel", length(genes))), " - click Generate to build the panel report.")
        } else {
          tagList(icon("circle-info"), "Enter or select genes for the panel above (paste a list, or use \"Use all listed genes as panel\") before clicking Generate.")
        }
      } else {
        gene <- if (identical(smode, "gene")) trimws(input$bmc_gene_input %||% "") else bmc_picked_gene()
        if (!is.null(gene) && nzchar(gene)) {
          tagList(icon("circle-check", style = "color:#0ca30c;"), tags$b(sprintf("Selected: %s", gene)), " - click Generate to build the card.")
        } else {
          tagList(icon("circle-info"), "Select a biomarker above (click a table row, or type a gene symbol) before clicking Generate.")
        }
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
    observeEvent(input$bmc_dge_use_panel, {
      run <- (results$dge_runs %||% list())[[input$bmc_dge_run]]; req(run)
      df <- run$table[order(run$table$adj.P.Val), , drop = FALSE]
      bmc_panel_genes(utils::head(df, 500)$gene)
    }, ignoreInit = TRUE)

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
    observeEvent(input$bmc_cand_use_panel, {
      cand <- results$candidates; req(cand); set <- cand[[input$bmc_cand_set %||% "final"]]; req(set)
      bmc_panel_genes(set$genes)
    }, ignoreInit = TRUE)

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
    observeEvent(input$bmc_sig_use_panel, {
      fs <- results$featureselection; req(fs); s <- fs[[input$bmc_sig_sex %||% ""]]; req(s)
      bmc_panel_genes(s$consensus_genes)
    }, ignoreInit = TRUE)

    ## ---- Upload a gene list, or a Diagnostic Classifier RDS export ($genes/$model_type) ----
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

    ## Identifier resolution for the uploaded list (spec: report resolution
    ## status per identifier, never silently drop one) - reuses tbc_panel_identity(),
    ## the same harmonizer the Gene Panel mode uses.
    upload_resolution <- reactive({
      df <- upload_table(); req(df)
      tbc_panel_identity(df$gene)
    })

    output$bmc_upload_results_ui <- renderUI({
      req(input$bmc_upload_load_btn)
      res <- upload_resolution()
      tagList(p(class = "submodule-desc",
                if (isTRUE(res$ok))
                  sprintf("%d gene(s) loaded from the uploaded file - %d resolved, %d unresolved (see status column below). Click a row to select a single gene, then click \"Generate Biomarker Card\" - or switch Mode to \"Gene panel\" above and use the whole list.", res$n_submitted, res$n_resolved, res$n_unresolved)
                else
                  sprintf("%d gene(s) loaded from the uploaded file (identifier resolution unavailable: %s). Click a row to select it, then click \"Generate Biomarker Card\".", nrow(upload_table()), res$reason %||% "unknown reason")),
              DT::dataTableOutput(ns("bmc_upload_table")))
    })
    output$bmc_upload_table <- DT::renderDataTable({
      df <- upload_table(); req(df)
      res <- upload_resolution()
      show <- if (isTRUE(res$ok)) {
        m <- res$df[match(df$gene, res$df$input_id), c("status_label", "canonical_symbol")]
        data.frame(gene = df$gene, `Resolution status` = m$status_label, `Canonical symbol` = m$canonical_symbol, check.names = FALSE, stringsAsFactors = FALSE)
      } else df
      DT::datatable(show, rownames = FALSE, selection = "single", options = list(pageLength = 10, scrollX = TRUE))
    })
    outputOptions(output, "bmc_upload_table", suspendWhenHidden = FALSE)
    observeEvent(input$bmc_upload_table_rows_selected, {
      df <- upload_table(); req(df)
      idx <- input$bmc_upload_table_rows_selected
      if (length(idx) == 1) bmc_picked_gene(df$gene[idx])
    })
    observeEvent(input$bmc_upload_use_panel, {
      df <- upload_table(); req(df)
      bmc_panel_genes(df$gene)
    }, ignoreInit = TRUE)

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
      selected_run_id <- if (identical(isolate(input$bmc_search_mode), "dge") && !is.null(isolate(input$bmc_dge_run))) {
        isolate(input$bmc_dge_run)
      } else if (!is.null(results$dge_runs) && length(results$dge_runs) > 0) {
        utils::tail(names(results$dge_runs), 1)
      } else NULL
      selected_run_table <- if (!is.null(selected_run_id)) (results$dge_runs %||% list())[[selected_run_id]]$table else NULL

      candidate_status <- tbc_candidate_status(gene, results)
      signature_membership <- tbc_signature_membership(gene, results)
      diagnostic_match <- tbc_diagnostic_lookup(gene, results)

      ## Gene-identity enrichment (§3) - own small, cached, fail-soft lookups.
      ncbi_summary <- if (isTRUE(gene_identity$ok)) tbc_ncbi_gene_summary(gene_identity$entrez) else list(ok = FALSE, reason = "Gene identity unresolved.")
      uniprot_name <- tbc_uniprot_protein_name(gene)

      ## Single-gene diagnostic performance (§7/§8) - live, computed here on
      ## the currently loaded dataset; NULL/ok=FALSE if no clean case/control split exists.
      single_gene_diag <- if (isTRUE(live$ok) && isTRUE(live$overall$ok)) {
        tbc_single_gene_roc(live$values, live$group_vec, live$case_label, live$control_label)
      } else list(ok = FALSE, reason = "No usable case/control split in the currently loaded dataset for this gene (see Dataset tab).")
      single_gene_cv <- if (isTRUE(live$ok) && isTRUE(live$overall$ok)) {
        tbc_single_gene_cv(live$values, live$group_vec, live$case_label, live$control_label)
      } else list(ok = FALSE, reason = "No usable case/control split in the currently loaded dataset for this gene (see Dataset tab).")

      list(gene = gene, in_dataset = in_dataset, gene_identity = gene_identity, live = live,
           dge_hits = dge_hits, selected_run_table = selected_run_table, selected_run_id = selected_run_id,
           candidate_status = candidate_status, signature_membership = signature_membership,
           diagnostic_match = diagnostic_match, ncbi_summary = ncbi_summary, uniprot_name = uniprot_name,
           single_gene_diag = single_gene_diag, single_gene_cv = single_gene_cv)
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

    ## ---- Generate (Gene Panel mode) ----
    observeEvent(input$bmc_generate_panel_btn, {
      has_panel_card(TRUE)
      updateTabsetPanel(session, "bmc_subtabs", selected = "Biomarker Card")
    }, ignoreInit = TRUE)

    panel_card_data <- eventReactive(input$bmc_generate_panel_btn, {
      raw_genes <- switch(input$bmc_search_mode,
        gene = tbc_split_gene_text(input$bmc_panel_gene_input),
        dge = , candidates = , signature = , upload = bmc_panel_genes(),
        character(0)
      )
      validate(need(length(raw_genes) > 0, "Enter or select at least one gene for the panel first - paste a gene list, or use \"Use all listed genes as panel\" on a browse tab."))

      pid <- tbc_panel_identity(raw_genes)
      validate(need(isTRUE(pid$ok), pid$reason %||% "Gene identifier resolution failed."))
      validate(need(pid$n_resolved > 0, sprintf("None of the %d submitted identifier(s) could be resolved to a known human gene - check spelling or identifier type.", pid$n_submitted)))

      resolved_symbols <- unique(pid$df$canonical_symbol[pid$df$resolved])
      gene_rows <- tbc_panel_gene_rows(resolved_symbols, dataset, results)

      ## Universe = genes measured in the currently loaded dataset, same
      ## convention mod_enrichment.R's own panel-enrichment flow already uses.
      universe_symbols <- tryCatch(rownames(dataset$expr), error = function(e) character(0))
      universe_map <- if (length(universe_symbols) > 0) suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = universe_symbols, keytype = "SYMBOL", columns = "ENTREZID")) else NULL
      universe_entrez <- if (!is.null(universe_map)) unique(stats::na.omit(universe_map$ENTREZID)) else character(0)

      entrez_in_universe <- unique(stats::na.omit(pid$df$entrez_id[pid$df$resolved & pid$df$entrez_id %in% universe_entrez]))

      ensembl_map <- stats::setNames(pid$df$ensembl_id, pid$df$canonical_symbol)
      ensembl_map <- ensembl_map[pid$df$resolved & !is.na(ensembl_map) & !duplicated(names(ensembl_map))]

      list(pid = pid, resolved_symbols = resolved_symbols, gene_rows = gene_rows,
           universe_entrez = universe_entrez, entrez_in_universe = entrez_in_universe,
           ensembl_map = as.list(ensembl_map),
           universe_label = sprintf("%s gene(s) measured in the currently loaded dataset", format(length(universe_entrez), big.mark = ",")))
    }, ignoreInit = TRUE)

    panel_enrich_data <- reactiveVal(NULL)
    panel_net_data <- reactiveVal(NULL)
    panel_disease_data <- reactiveVal(NULL)
    panel_drug_data <- reactiveVal(NULL)
    observeEvent(panel_card_data(), {
      panel_enrich_data(NULL); panel_net_data(NULL); panel_disease_data(NULL); panel_drug_data(NULL)
    }, ignoreInit = TRUE)

    observeEvent(input$run_panel_enrichment, {
      pcd <- panel_card_data(); req(pcd)
      panel_enrich_data(tbc_panel_enrichment(pcd$entrez_in_universe, pcd$universe_entrez))
    }, ignoreInit = TRUE)

    observeEvent(input$run_panel_network, {
      pcd <- panel_card_data(); req(pcd)
      panel_net_data(tbc_panel_network(pcd$resolved_symbols))
    }, ignoreInit = TRUE)

    observeEvent(input$run_panel_disease, {
      pcd <- panel_card_data(); req(pcd)
      panel_disease_data(tbc_panel_disease_convergence(pcd$ensembl_map))
    }, ignoreInit = TRUE)

    observeEvent(input$run_panel_drug, {
      pcd <- panel_card_data(); req(pcd)
      panel_drug_data(tbc_panel_drug_convergence(pcd$resolved_symbols))
    }, ignoreInit = TRUE)

    ## Each external database is its own independent, user-triggered live API call.
    ot_data <- reactiveVal(NULL)
    hpa_data <- reactiveVal(NULL)
    string_data <- reactiveVal(NULL)
    string_image_path <- reactiveVal(NULL)
    dgidb_data <- reactiveVal(NULL)
    ## GO/KEGG/Reactome/WikiPathways are each their own independent, opt-in
    ## lookup - one Run button per database (matching Open Targets/DGIdb/
    ## STRING/HPA/PubMed below), not one combined query for all four.
    go_data <- reactiveVal(NULL)
    kegg_data <- reactiveVal(NULL)
    reactome_data <- reactiveVal(NULL)
    wikipathways_data <- reactiveVal(NULL)
    kegg_map_data <- reactiveVal(NULL)
    reactome_map_data <- reactiveVal(NULL)
    literature_data <- reactiveVal(NULL)
    literature_query_used <- reactiveVal(NULL)
    observeEvent(card_data(), {
      ot_data(NULL); hpa_data(NULL); string_data(NULL); string_image_path(NULL); dgidb_data(NULL)
      go_data(NULL); kegg_data(NULL); reactome_data(NULL); wikipathways_data(NULL)
      kegg_map_data(NULL); reactome_map_data(NULL); literature_data(NULL); literature_query_used(NULL)
    }, ignoreInit = TRUE)

    observeEvent(input$run_ot, {
      d <- card_data(); req(d)
      ot_data(tbc_opentargets_evidence_for_gene(d$gene_identity$ensembl %||% NA_character_, top_n_diseases = input$ot_topn %||% 8))
    }, ignoreInit = TRUE)

    observeEvent(input$run_hpa, {
      d <- card_data(); req(d)
      hpa_data(tbc_hpa_evidence_for_gene(d$gene_identity$ensembl %||% NA_character_))
    }, ignoreInit = TRUE)

    observeEvent(input$run_string, {
      d <- card_data(); req(d)
      cutoff <- input$string_cutoff %||% 400; topn <- input$string_topn %||% 10
      string_data(tbc_string_ppi_for_gene(d$gene, top_n = topn, required_score = cutoff))
      string_image_path(NULL)
      img <- tbc_string_network_image(d$gene, top_n = topn, required_score = cutoff)
      if (isTRUE(img$ok)) string_image_path(img$path)
    }, ignoreInit = TRUE)

    observeEvent(input$run_dgidb, {
      d <- card_data(); req(d)
      dgidb_data(tbc_dgidb_drugs_for_gene(d$gene, top_n = input$dgidb_topn %||% 12))
    }, ignoreInit = TRUE)

    observeEvent(input$run_go, {
      d <- card_data(); req(d)
      go_data(tbc_go_terms(d$gene_identity$entrez %||% NA_character_))
    }, ignoreInit = TRUE)

    observeEvent(input$run_kegg, {
      d <- card_data(); req(d)
      kegg_map_data(NULL)
      kegg_data(tbc_kegg_pathways_for_gene(d$gene_identity$entrez %||% NA_character_))
    }, ignoreInit = TRUE)

    observeEvent(input$run_reactome, {
      d <- card_data(); req(d)
      reactome_map_data(NULL)
      reactome_data(tbc_reactome_pathways_for_gene(d$gene))
    }, ignoreInit = TRUE)

    observeEvent(input$run_wikipathways, {
      d <- card_data(); req(d)
      wikipathways_data(tbc_wikipathways_pathways_for_gene(d$gene_identity$entrez %||% NA_character_))
    }, ignoreInit = TRUE)

    ## Diagram rendering is opt-in per pathway; KEGG map is colored by whichever
    ## log2FC evidence is already on the card (live preview, else saved DGE logFC).
    observeEvent(input$render_kegg_map, {
      d <- card_data(); req(d)
      kd <- kegg_data(); req(kd); req(input$kegg_map_pick)
      log2fc <- if (isTRUE(d$live$ok) && isTRUE(d$live$overall$ok)) d$live$overall$log2fc
                else if (!is.null(d$dge_hits) && nrow(d$dge_hits) > 0) d$dge_hits$logFC[1] else NA_real_
      pw_name <- kd$pathways$name[kd$pathways$id == input$kegg_map_pick][1]
      kegg_map_data(tbc_kegg_diagram_for_gene(input$kegg_map_pick, d$gene_identity$entrez %||% NA_character_, log2fc, pw_name))
    }, ignoreInit = TRUE)

    observeEvent(input$render_reactome_map, {
      rd <- reactome_data(); req(rd); req(input$reactome_map_pick)
      pw_name <- rd$pathways$displayName[rd$pathways$stId == input$reactome_map_pick][1]
      reactome_map_data(tbc_reactome_diagram_for_pathway(input$reactome_map_pick, pw_name))
    }, ignoreInit = TRUE)

    observeEvent(input$run_literature, {
      d <- card_data(); req(d)
      q <- tbc_literature_query(d$gene, input$literature_preset %||% "Biomarker", input$literature_disease_text)
      literature_query_used(q)
      literature_data(tbc_literature_search(q, max_results = input$literature_topn %||% 12))
    }, ignoreInit = TRUE)

    ## Combined view of whichever per-database results have been run so far this session.
    ext_data <- reactive({
      list(kegg = kegg_data(), reactome = reactome_data(), wikipathways = wikipathways_data(), go = go_data(),
           genetics = ot_data(), tissue = hpa_data(), network = string_data(), drugs = dgidb_data(),
           literature = literature_data(), literature_query = literature_query_used())
    })

    output$bmc_card_ui <- renderUI({
      mode <- input$bmc_mode %||% "single"

      if (identical(mode, "panel")) {
        if (!has_panel_card()) return(div(class = "empty-note", icon("circle-info"), "Select genes for a panel and click \"Generate Gene Panel Report\" on the \"Select Biomarker\" tab."))
        pcd <- panel_card_data(); req(pcd)
        return(tagList(
          tbc_section_panel_identity(pcd$pid),
          tbc_section_panel_overview(pcd$pid, pcd$gene_rows),
          div(class = "card",
              div(class = "card-title", icon("globe"), "Deep Dive: Panel-Level Databases"),
              p(class = "submodule-desc", "Each section below is computed once, on demand, against a background of genes actually measured in your loaded dataset."),
              radioButtons(ns("bmc_panel_db_choice"), NULL, inline = TRUE, choices = c(
                "GO / KEGG / Reactome / WikiPathways (enrichment)" = "enrichment",
                "Network & Hub Genes" = "network",
                "Disease Convergence (Open Targets)" = "disease",
                "Drug/Target Convergence (DGIdb)" = "drug"
              )),
              uiOutput(ns("bmc_panel_db_controls_ui"))
          ),
          withSpinner(uiOutput(ns("bmc_panel_db_result_ui")), color = "#2563EB", type = 6),
          tbc_section_panel_evidence_summary(panel_enrich_data(), panel_net_data(), panel_disease_data(), panel_drug_data()),
          tbc_section_panel_comparison(panel_enrich_data(), panel_net_data(), panel_disease_data(), panel_drug_data()),
          div(class = "card", div(class = "card-title", icon("download"), "Download"),
              downloadButton(ns("bmc_download_panel_report"), "Download Gene Panel Report (HTML)", class = "btn-primary btn-sm"))
        ))
      }

      if (!has_card()) return(div(class = "empty-note", icon("circle-info"), "Select a biomarker and click \"Generate Biomarker Card\" on the \"Select Biomarker\" tab."))
      d <- card_data(); req(d)
      tabsetPanel(id = ns("bmc_result_tabs"), type = "tabs",
        tabPanel("Gene description", br(), tbc_section_identity(d)),
        tabPanel("Biomarker Status", br(), withSpinner(uiOutput(ns("bmc_evidence_glance_ui")), color = "#2563EB", type = 6)),
        tabPanel("Dataset", br(), tbc_section_dataset_cohort(dataset, d$live)),
        tabPanel("Differential Expression", br(),
          tbc_section_differential_expression(d),
          div(class = "card",
              div(class = "card-title", icon("chart-column"), "Volcano Plot"),
              if (is.null(d$selected_run_table)) div(class = "empty-note", icon("circle-info"), "No Differential Expression run available - run the Differential Expression tab first.")
              else withSpinner(plotOutput(ns("bmc_volcano_plot"), height = "380px"), color = "#2563EB", type = 6)
          ),
          fluidRow(
            column(6, div(class = "card",
                div(class = "card-title", icon("chart-simple"), "Expression Distribution"),
                withSpinner(plotly::plotlyOutput(ns("bmc_dist_plot"), height = "300px"), color = "#2563EB", type = 6)
            )),
            column(6, div(class = "card",
                div(class = "card-title", icon("venus-mars"), "Sex-Specific Expression"),
                withSpinner(plotly::plotlyOutput(ns("bmc_sex_dist_plot"), height = "300px"), color = "#2563EB", type = 6)
            ))
          )
        ),
        tabPanel("Single-Gene vs Multi-Gene Signature", br(),
          tbc_section_signature_comparison(d, d$single_gene_diag,
            if (isTRUE(d$single_gene_diag$ok)) withSpinner(plotOutput(ns("bmc_sg_train_roc_plot"), height = "320px"), color = "#2563EB", type = 6) else NULL)
        ),
        tabPanel("Biomarker Performance", br(),
          tbc_section_biomarker_performance(d, d$single_gene_diag, d$single_gene_cv,
            train_roc_widget = if (isTRUE(d$single_gene_diag$ok)) withSpinner(plotOutput(ns("bmc_sg_train_roc_plot2"), height = "320px"), color = "#2563EB", type = 6) else NULL,
            train_pr_widget = if (isTRUE(d$single_gene_diag$ok)) withSpinner(plotOutput(ns("bmc_sg_train_pr_plot"), height = "320px"), color = "#2563EB", type = 6) else NULL,
            internal_roc_widget = if (isTRUE(d$single_gene_cv$ok)) withSpinner(plotOutput(ns("bmc_sg_internal_roc_plot"), height = "320px"), color = "#2563EB", type = 6) else NULL)
        ),
        tabPanel("External Databases", br(),
          tbc_external_db_banner(),
          div(class = "card",
              div(class = "card-title", icon("globe"), "Deep Dive: External Databases"),
              p(class = "submodule-desc", "Pick a database, set your own query parameters, and run it live - each is independent and links straight to the real record so you can verify it yourself."),
              radioButtons(ns("bmc_db_choice"), NULL, inline = TRUE, choices = c(
                "GO" = "go", "KEGG" = "kegg", "Reactome" = "reactome", "WikiPathways" = "wikipathways",
                "Disease / Genetics (Open Targets)" = "opentargets", "Drug / Target (DGIdb)" = "dgidb",
                "Protein / Interaction (STRING)" = "string", "Expression (Human Protein Atlas)" = "hpa",
                "Literature (PubMed)" = "literature"
              )),
              uiOutput(ns("bmc_db_controls_ui"))
          ),
          withSpinner(uiOutput(ns("bmc_db_result_ui")), color = "#2563EB", type = 6),
          withSpinner(uiOutput(ns("bmc_db_comparison_ui")), color = "#2563EB", type = 6)
        ),
        tabPanel("Download", br(),
          div(class = "card", div(class = "card-title", icon("download"), "Download"),
              p(class = "submodule-desc", "The downloaded report includes every section shown across all tabs above, for this gene."),
              downloadButton(ns("bmc_download_report"), "Download Biomarker Report (HTML)", class = "btn-primary btn-sm"))
        )
      )
    })

    output$bmc_db_controls_ui <- renderUI({
      choice <- input$bmc_db_choice %||% "go"
      switch(choice,
        string = tagList(
          fluidRow(
            column(6, numericInput(ns("string_cutoff"), "Confidence cutoff (required_score, 0-1000)", value = 400, min = 0, max = 1000, step = 50)),
            column(6, numericInput(ns("string_topn"), "Max partners", value = 10, min = 1, max = 50, step = 1))
          ),
          actionButton(ns("run_string"), "Run STRING Query", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        opentargets = tagList(
          numericInput(ns("ot_topn"), "Number of diseases to show", value = 8, min = 1, max = 25, step = 1),
          actionButton(ns("run_ot"), "Run Open Targets Query", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        hpa = tagList(
          p(class = "submodule-desc", "Human Protein Atlas has no adjustable query parameters here - it's a fixed gene-profile lookup. Shown here as an external reference; \"Your Dataset\" expression is on the Differential Expression section above."),
          actionButton(ns("run_hpa"), "Run Human Protein Atlas Query", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        dgidb = tagList(
          numericInput(ns("dgidb_topn"), "Max drug interactions", value = 12, min = 1, max = 50, step = 1),
          actionButton(ns("run_dgidb"), "Run DGIdb Query", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        literature = tagList(
          selectInput(ns("literature_preset"), "Query template",
                      choices = c(stats::setNames(names(TBC_LITERATURE_PRESETS), names(TBC_LITERATURE_PRESETS)), c("Disease (custom text)" = "disease"))),
          conditionalPanel(condition = sprintf("input['%s'] == 'disease'", ns("literature_preset")),
                            textInput(ns("literature_disease_text"), "Disease / condition", placeholder = "e.g. rheumatoid arthritis")),
          numericInput(ns("literature_topn"), "Max results", value = 12, min = 1, max = 30, step = 1),
          actionButton(ns("run_literature"), "Run PubMed Query", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        ## GO/KEGG/Reactome/WikiPathways are each independent, opt-in lookups -
        ## selecting one shows only that database's own Run button, exactly
        ## like Open Targets/DGIdb/STRING/HPA/PubMed above.
        go = tagList(
          p(class = "submodule-desc", "Gene Ontology (Biological Process) terms for this gene."),
          actionButton(ns("run_go"), "Run GO Query", icon = icon("play"), class = "btn-primary btn-sm")
        ),
        kegg = {
          kd <- kegg_data()
          kegg_hits <- if (isTRUE(kd$ok)) kd$pathways else NULL
          tagList(
            p(class = "submodule-desc", "KEGG pathway membership for this gene."),
            actionButton(ns("run_kegg"), "Run KEGG Query", icon = icon("play"), class = "btn-primary btn-sm"),
            if (!is.null(kegg_hits) && nrow(kegg_hits) > 0) tagList(
              tags$div(style = "margin-top:14px;", tags$b("Render a KEGG pathway diagram - real pathview render, gene colored by its own log2FC")),
              fluidRow(
                column(8, selectInput(ns("kegg_map_pick"), NULL, choices = stats::setNames(kegg_hits$id, sprintf("%s - %s", kegg_hits$id, kegg_hits$name)))),
                column(4, actionButton(ns("render_kegg_map"), "Render Diagram", icon = icon("image"), class = "btn-default btn-sm"))
              )
            ) else NULL
          )
        },
        reactome = {
          rd <- reactome_data()
          reactome_hits <- if (isTRUE(rd$ok)) rd$pathways else NULL
          tagList(
            p(class = "submodule-desc", "Reactome pathway membership for this gene."),
            actionButton(ns("run_reactome"), "Run Reactome Query", icon = icon("play"), class = "btn-primary btn-sm"),
            if (!is.null(reactome_hits) && nrow(reactome_hits) > 0) tagList(
              tags$div(style = "margin-top:14px;", tags$b("Render a Reactome pathway diagram - real server-side Reactome image")),
              fluidRow(
                column(8, selectInput(ns("reactome_map_pick"), NULL, choices = stats::setNames(reactome_hits$stId, sprintf("%s - %s", reactome_hits$stId, reactome_hits$displayName)))),
                column(4, actionButton(ns("render_reactome_map"), "Render Diagram", icon = icon("image"), class = "btn-default btn-sm"))
              )
            ) else NULL
          )
        },
        wikipathways = tagList(
          p(class = "submodule-desc", "WikiPathways pathway membership for this gene."),
          actionButton(ns("run_wikipathways"), "Run WikiPathways Query", icon = icon("play"), class = "btn-primary btn-sm")
        )
      )
    })

    output$bmc_db_result_ui <- renderUI({
      d <- card_data(); req(d)
      choice <- input$bmc_db_choice %||% "go"
      not_yet <- div(class = "empty-note", icon("circle-info"), "Not yet run - set your parameters above and click Run.")
      switch(choice,
        string = if (is.null(string_data())) not_yet else tbc_section_string(string_data(), string_image_path(), d$gene),
        opentargets = if (is.null(ot_data())) not_yet else tbc_section_genetics(ot_data(), d$gene, d$gene_identity$ensembl),
        hpa = if (is.null(hpa_data())) not_yet else tbc_section_tissue(hpa_data(), d$gene, d$gene_identity$ensembl),
        dgidb = if (is.null(dgidb_data())) not_yet else tbc_section_dgidb(dgidb_data(), d$gene),
        literature = if (is.null(literature_data())) not_yet else tbc_section_literature(literature_data(), d$gene, literature_query_used()),
        go = if (is.null(go_data())) not_yet else tbc_section_go(list(go = go_data())),
        kegg = if (is.null(kegg_data())) not_yet else tbc_section_kegg(list(kegg = kegg_data()), kegg_map_data()),
        reactome = if (is.null(reactome_data())) not_yet else tbc_section_reactome(list(reactome = reactome_data()), reactome_map_data()),
        wikipathways = if (is.null(wikipathways_data())) not_yet else tbc_section_wikipathways(list(wikipathways = wikipathways_data()))
      )
    })

    ## These two live in their own renderUI (not inline in bmc_card_ui) so
    ## that running an external-database lookup only re-renders this one tab's
    ## content, instead of rebuilding the whole bmc_result_tabs tabsetPanel and
    ## resetting the user back to its first tab (regression found during live
    ## Playwright verification - ext_data() must never be read directly inside
    ## the tabsetPanel-building renderUI).
    output$bmc_evidence_glance_ui <- renderUI({
      d <- card_data(); req(d)
      tbc_section_evidence_glance(d, ext_data(), d$single_gene_diag, d$single_gene_cv)
    })
    output$bmc_db_comparison_ui <- renderUI({
      tbc_section_db_comparison(ext_data())
    })

    output$bmc_panel_db_controls_ui <- renderUI({
      choice <- input$bmc_panel_db_choice %||% "enrichment"
      switch(choice,
        enrichment = tagList(p(class = "submodule-desc", "Over-representation analysis against genes measured in your loaded dataset (requires >=3 resolved genes measured there)."),
                              actionButton(ns("run_panel_enrichment"), "Run Enrichment (GO / KEGG / Reactome / WikiPathways)", icon = icon("play"), class = "btn-primary btn-sm")),
        network = actionButton(ns("run_panel_network"), "Run Network & Hub Analysis", icon = icon("play"), class = "btn-primary btn-sm"),
        disease = actionButton(ns("run_panel_disease"), "Run Disease Convergence (Open Targets)", icon = icon("play"), class = "btn-primary btn-sm"),
        drug = actionButton(ns("run_panel_drug"), "Run Drug/Target Convergence (DGIdb)", icon = icon("play"), class = "btn-primary btn-sm")
      )
    })

    output$bmc_panel_db_result_ui <- renderUI({
      pcd <- panel_card_data(); req(pcd)
      choice <- input$bmc_panel_db_choice %||% "enrichment"
      not_yet <- div(class = "empty-note", icon("circle-info"), "Not yet run - click Run above.")
      switch(choice,
        enrichment = if (is.null(panel_enrich_data())) not_yet else tbc_section_panel_enrichment(panel_enrich_data(), pcd$universe_label, length(pcd$entrez_in_universe)),
        network = if (is.null(panel_net_data())) not_yet else tbc_section_panel_network(panel_net_data(), panel_net_data()$image_path),
        disease = if (is.null(panel_disease_data())) not_yet else tbc_section_panel_disease(panel_disease_data()),
        drug = if (is.null(panel_drug_data())) not_yet else tbc_section_panel_drug(panel_drug_data())
      )
    })

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

    ## Single-gene ROC/PR plots (Training = full-fit, Internal Validation =
    ## pooled out-of-fold CV) - the same plot data is shown on both the
    ## Single-Gene vs Multi-Gene Signature tab and the Biomarker Performance
    ## tab's Training block, via two output ids (bmc_sg_train_roc_plot/plot2).
    output$bmc_sg_train_roc_plot <- renderPlot({
      d <- card_data(); req(d); sgd <- d$single_gene_diag
      validate(need(isTRUE(sgd$ok), sgd$reason %||% "Single-gene ROC unavailable."))
      tbc_plot_roc(sgd$roc_obj, "Single-gene ROC (Training)", sgd$auc, sgd$ci_lo, sgd$ci_hi)
    })
    output$bmc_sg_train_roc_plot2 <- renderPlot({
      d <- card_data(); req(d); sgd <- d$single_gene_diag
      validate(need(isTRUE(sgd$ok), sgd$reason %||% "Single-gene ROC unavailable."))
      tbc_plot_roc(sgd$roc_obj, "Single-gene ROC (Training)", sgd$auc, sgd$ci_lo, sgd$ci_hi)
    })
    output$bmc_sg_train_pr_plot <- renderPlot({
      d <- card_data(); req(d); sgd <- d$single_gene_diag
      validate(need(isTRUE(sgd$ok), sgd$reason %||% "Single-gene PR curve unavailable."))
      pr <- tbc_pr_from_roc(sgd$roc_obj)
      validate(need(isTRUE(pr$ok), pr$reason %||% "Precision-recall curve unavailable."))
      tbc_plot_pr(pr)
    })
    output$bmc_sg_internal_roc_plot <- renderPlot({
      d <- card_data(); req(d); sgcv <- d$single_gene_cv
      validate(need(isTRUE(sgcv$ok), sgcv$reason %||% "Cross-validated ROC unavailable."))
      tbc_plot_roc(sgcv$roc_obj, "Single-gene ROC (Internal Validation, pooled CV)", sgcv$auc)
    })

    output$bmc_download_report <- downloadHandler(
      filename = function() sprintf("transcriptomic_biomarker_card_%s.html", card_data()$gene),
      content = function(file) {
        d <- card_data(); req(d)
        body <- tbc_build_report_tags(d, dataset, ext_data())
        page <- tags$html(
          tags$head(tags$meta(charset = "utf-8"), tags$title(sprintf("Transcriptomic Biomarker Card %s", d$gene)), tags$style(tbc_report_css())),
          tags$body(body)
        )
        htmltools::save_html(page, file = file)
      }
    )

    output$bmc_download_panel_report <- downloadHandler(
      filename = function() sprintf("gene_panel_biomarker_report_%s.html", format(Sys.Date())),
      content = function(file) {
        pcd <- panel_card_data(); req(pcd)
        body <- tbc_build_panel_report_tags(pcd$pid, pcd$gene_rows, panel_enrich_data(), panel_net_data(), panel_disease_data(), panel_drug_data(),
                                             pcd$universe_label, length(pcd$entrez_in_universe))
        page <- tags$html(
          tags$head(tags$meta(charset = "utf-8"), tags$title("Gene Panel Biomarker Report"), tags$style(tbc_report_css())),
          tags$body(body)
        )
        htmltools::save_html(page, file = file)
      }
    )
  })
}
