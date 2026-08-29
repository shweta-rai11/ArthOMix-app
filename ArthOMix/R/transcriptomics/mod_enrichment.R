## Functional Enrichment submodule: live GO/KEGG/Reactome over-representation
## analysis on a user gene list or a bundled biomarker panel, plus MyGene.info
## gene cards and a STRING PPI network for the same genes.

mod_enrichment_config <- list(
  id = "enrichment", group = "Interpretation",
  title = "Functional Enrichment",
  description = "Perform GO, KEGG or Reactome over-representation analysis on a gene list either uploaded or pre-loaded.",
  icon = "diagram-project"
)

## Live external lookups (called once per "Run enrichment" click); wrapped so
## an unreachable service degrades to an empty result rather than failing.

## Batched MyGene.info query for the whole gene list; returns GeneCards-style fields.
fetch_gene_cards <- function(genes) {
  tryCatch({
    resp <- httr::POST(
      "https://mygene.info/v3/query",
      body = list(q = paste(genes, collapse = ","), scopes = "symbol", species = "human",
                  fields = "symbol,name,summary,alias,genomic_pos,entrezgene"),
      encode = "form", httr::timeout(20)
    )
    if (httr::http_error(resp)) return(list())
    hits <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
    if (!is.data.frame(hits) || !nrow(hits)) return(list())
    lapply(seq_len(nrow(hits)), function(i) {
      found <- !isTRUE(hits$notfound[i])
      alias_i <- if ("alias" %in% names(hits)) hits$alias[[i]] else NULL
      list(
        symbol  = hits$query[i],
        found   = found,
        name    = if (found && "name" %in% names(hits)) hits$name[i] else NA_character_,
        summary = if (found && "summary" %in% names(hits)) hits$summary[i] else NA_character_,
        alias   = if (is.null(alias_i) || !length(alias_i)) NA_character_ else paste(unlist(alias_i), collapse = ", "),
        chr     = if (found && "genomic_pos" %in% names(hits)) hits$genomic_pos$chr[i] else NA_character_,
        entrez  = if (found && "entrezgene" %in% names(hits)) hits$entrezgene[i] else NA_character_
      )
    })
  }, error = function(e) list())
}

## Fetches the STRING PPI network as a pre-rendered PNG from string-db.org.
fetch_string_network_png <- function(genes) {
  tryCatch({
    ids <- utils::URLencode(paste(genes, collapse = "\r"), reserved = TRUE)
    url <- sprintf("https://string-db.org/api/image/network?identifiers=%s&species=9606&required_score=400", ids)
    tmp <- tempfile(fileext = ".png")
    resp <- httr::GET(url, httr::write_disk(tmp, overwrite = TRUE), httr::timeout(20))
    if (httr::http_error(resp) || !file.exists(tmp) || file.info(tmp)$size < 500) return(NULL)
    tmp
  }, error = function(e) NULL)
}

## Hub genes: two independent measures - live STRING PPI degree within just
## this gene list, vs. this project's precomputed WGCNA co-expression hub
## status across the full training cohort (see build_wgcna_hub_lookup below).

## STRING PPI degree of each submitted gene against the rest of the list (confidence >= 0.4).
fetch_string_degrees <- function(genes) {
  tryCatch({
    ids <- utils::URLencode(paste(genes, collapse = "\r"), reserved = TRUE)
    url <- sprintf("https://string-db.org/api/tsv/network?identifiers=%s&species=9606&required_score=400", ids)
    resp <- httr::GET(url, httr::timeout(20))
    if (httr::http_error(resp)) return(list())
    txt <- httr::content(resp, "text", encoding = "UTF-8")
    if (!nzchar(trimws(txt))) return(list())
    edges <- utils::read.delim(text = txt, sep = "\t", stringsAsFactors = FALSE)
    if (!nrow(edges) || !all(c("preferredName_A", "preferredName_B") %in% names(edges))) return(list())
    partners <- setNames(vector("list", length(genes)), genes)
    for (g in genes) partners[[g]] <- character(0)
    for (i in seq_len(nrow(edges))) {
      a <- edges$preferredName_A[i]; b <- edges$preferredName_B[i]
      if (a %in% genes && b %in% genes) {
        partners[[a]] <- union(partners[[a]], b)
        partners[[b]] <- union(partners[[b]], a)
      }
    }
    lapply(partners, function(p) list(degree = length(p), partners = paste(sort(p), collapse = ", ")))
  }, error = function(e) list())
}

## Loads the bundled WGCNA module/hub tables once per session, returns a per-gene lookup fn.
build_wgcna_hub_lookup <- function() {
  modules <- read_table_safe("WGCNA_05_gene_module_assignment.csv")
  hubs    <- read_table_safe("WGCNA_07_hub_genes_only.csv")
  function(gene) {
    m <- if (!is.null(modules)) modules[modules$gene == gene, , drop = FALSE] else NULL
    h <- if (!is.null(hubs)) hubs[hubs$gene == gene, , drop = FALSE] else NULL
    list(
      module         = if (!is.null(m) && nrow(m)) m$module[1] else NA_character_,
      disease_module = if (!is.null(m) && nrow(m)) isTRUE(as.character(m$is_disease_module[1]) == "TRUE") else FALSE,
      kme            = if (!is.null(m) && nrow(m)) suppressWarnings(as.numeric(m$kME_own[1])) else NA_real_,
      is_hub         = if (!is.null(h) && nrow(h)) isTRUE(as.character(h$is_hub[1]) == "TRUE") else FALSE,
      connectivity   = if (!is.null(h) && nrow(h)) suppressWarnings(as.numeric(h$connectivity[1])) else NA_real_
    )
  }
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------

mod_enrichment_ui <- function(id) {
  ns <- NS(id)
  wrap_id <- ns("wrap")

  ## gsub() instead of sprintf() - the CSS below has literal "%" widths.
  enrichment_css <- gsub("__WRAP__", wrap_id, fixed = TRUE, "
      #__WRAP__ .gene-chip-row { display: flex; flex-wrap: wrap; gap: 6px; margin: 10px 0; }
      #__WRAP__ .gene-chip {
        background: var(--color-surface-blue); color: var(--color-primary-dark);
        font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 12px; font-weight: 600;
        padding: 4px 9px; border-radius: 20px; line-height: 1.4;
      }
      #__WRAP__ .preset-panel {
        background: #F8FAFC; border: 1px solid var(--color-border); border-radius: var(--radius-card);
        padding: 14px 16px; margin-bottom: 14px;
      }
      #__WRAP__ .preset-panel .btn-group.btn-group-justified { margin-top: 8px; }
      /* 4-way gene-source picker wraps into a 2x2 grid; labels clip in a single row at this sidebar width. */
      #__WRAP__ .source-picker-grid .btn-group-justified {
        display: grid; grid-template-columns: 1fr 1fr; gap: 6px;
      }
      #__WRAP__ .source-picker-grid .btn-group-justified > .btn-group { display: block; width: 100%; }
      #__WRAP__ .source-picker-grid .btn-group-justified .btn {
        width: 100%; white-space: normal; font-size: 12px; line-height: 1.3; padding: 8px 6px;
      }
      #__WRAP__ .gene-card-item { padding: 14px 0; border-bottom: 1px solid var(--color-border); }
      #__WRAP__ .gene-card-item:last-child { border-bottom: none; }
      #__WRAP__ .gene-card-item .card-title a { color: var(--color-ink); text-decoration: none; }
      #__WRAP__ .gene-card-item .card-title a:hover { color: var(--color-primary); }
      #__WRAP__ .gene-card-chip-row { display: flex; flex-wrap: wrap; gap: 6px; margin: 4px 0 8px; }
      #__WRAP__ .gene-card-chip {
        background: #F4F6F9; color: var(--color-ink-secondary); font-size: 11.5px;
        font-weight: 600; padding: 2px 8px; border-radius: 4px;
      }
      #__WRAP__ .gene-card-summary { font-size: 13px; color: var(--color-ink-secondary); line-height: 1.55; margin: 0 0 10px; }
      #__WRAP__ .run-btn-row { margin-top: 14px; }
  ")

  tagList(
    tags$style(HTML(enrichment_css)),

    div(id = wrap_id,
      fluidRow(
        column(
          4,
          box(
            width = NULL, title = tagList(icon("dna"), "Gene list"), status = "primary", solidHeader = FALSE,
            div(class = "filter-section-header", "GENE SOURCE"),
            div(class = "source-picker-grid",
              shinyWidgets::radioGroupButtons(
                ns("gene_source"), NULL,
                choices = c("Your own" = "own", "Same-tissue" = "same_tissue",
                            "Cross-tissue" = "cross_tissue", "Cross-ancestry" = "cross_ancestry"),
                selected = "own", status = "default", size = "sm", justified = TRUE
              )
            ),
            uiOutput(ns("preset_panel_ui")),
            p(class = "card-subtitle", style = "margin-top: 4px;",
              "Paste gene symbols below (one per line or comma separated), or load a panel above."),
            textAreaInput(ns("gene_list"), NULL, rows = 7, placeholder = "TNF\nIL6\nSTAT3\n..."),
            hr(),
            div(class = "filter-section-header", "DATABASE"),
            radioButtons(
              ns("ontology"), NULL,
              choices = c("GO Biological Process" = "GO_BP", "KEGG pathways" = "KEGG", "Reactome pathways" = "REACTOME"),
              selected = "GO_BP"
            ),
            numericInput(ns("qval_cut"), "Q-value cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
            div(class = "run-btn-row",
                actionButton(ns("run_btn"), "Run enrichment", icon = icon("play"), class = "btn-primary btn-sm", width = "100%"))
          )
        ),
        column(
          8,
          uiOutput(ns("result_box_ui"))
        )
      ),
      uiOutput(ns("below_result_ui"))
    )
  )
}

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------

mod_enrichment_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Gene-source panels (same-tissue, cross-tissue, cross-ancestry), each
    ## sex-stratified and preferring this session's live results over the
    ## bundled fallback, same convention as mod_diagnostic.R/mod_crosstissue.R.
    ## Every bundled_* fallback below was only ever computed from the app's own
    ## default merged cohort - only read/offer it when that exact dataset is
    ## still active (dataset$is_bundled_reference), never for an uploaded,
    ## GEO-fetched, or individual raw preloaded dataset. The reactiveValues
    ## read itself can only happen inside a reactive consumer, so loading these
    ## (cheap, static files - harmless regardless of active dataset) stays
    ## unconditional here at module-setup time; the gate is applied instead at
    ## each function's own call site below (same_tissue_panel/etc. are already
    ## only invoked lazily from inside a reactive context) and around the
    ## wgcna_hub_lookup() call further down.
    bundled_synovium <- tryCatch(readRDS(VAL_SYNOVIUM_RDS), error = function(e) NULL)
    bundled_venn <- read_table_safe("FS_venn_membership.csv")
    wgcna_hub_lookup <- build_wgcna_hub_lookup()
    wgcna_hub_lookup_gated <- function(gene) {
      if (isTRUE(dataset$is_bundled_reference)) wgcna_hub_lookup(gene)
      else list(module = NA_character_, disease_module = FALSE, kme = NA_real_, is_hub = FALSE, connectivity = NA_real_)
    }

    same_tissue_panel <- function(sex_label) {
      live <- results$featureselection[[sex_label]]$consensus_genes
      if (!is.null(live) && length(live) >= 2) {
        return(list(genes = live, is_live = TRUE,
          note = sprintf("%d genes - this session's live Feature Selection %s consensus panel.", length(live), sex_label)))
      }
      if (!is.null(bundled_venn)) {
        genes <- unique(bundled_venn$gene[
          bundled_venn$candidate_set == "noMHC" & bundled_venn$sex == sex_label &
          as.character(bundled_venn$in_consensus) == "TRUE"
        ])
        if (length(genes) >= 2) {
          return(list(genes = genes, is_live = FALSE,
            note = sprintf("%d genes - this project's bundled %s diagnostic consensus panel (MHC-excluded, recommended in NESTED_CV_AUTHORITATIVE.csv).", length(genes), sex_label)))
        }
      }
      list(genes = character(0), is_live = FALSE,
           note = sprintf("No live %s panel yet - run Feature Selection on the currently loaded dataset first.", sex_label))
    }

    cross_tissue_panel <- function(sex_label) {
      live <- results$featureselection[[sex_label]]$consensus_genes
      if (!is.null(live) && length(live) >= 2) {
        return(list(genes = live, is_live = TRUE,
          note = sprintf("%d genes - this session's live Feature Selection %s consensus panel.", length(live), sex_label)))
      }
      bundled <- switch(sex_label, female = bundled_synovium$fsig, male = bundled_synovium$msig, NULL)
      bundled <- unique(as.character(bundled))
      if (length(bundled) >= 2) {
        return(list(genes = bundled, is_live = FALSE,
          note = sprintf("%d genes - this project's bundled %s panel, the one its own synovial validation used (MHC-retained).", length(bundled), sex_label)))
      }
      list(genes = character(0), is_live = FALSE,
           note = sprintf("No live %s panel yet - run Feature Selection on the currently loaded dataset first.", sex_label))
    }

    cross_ancestry_panel <- function(sex_label) {
      live <- results$crossancestry[[sex_label]]$biomarker_genes
      if (!is.null(live) && length(live) >= 2) {
        return(list(genes = live, is_live = TRUE,
          note = sprintf("%d genes - this session's live Cross-Ancestry replicated-biomarker list.", length(live))))
      }
      if (isTRUE(dataset$is_bundled_reference)) {
        bundled <- read_table_safe(sprintf("MR35_crossancestry_%s.csv", sex_label))
        if (!is.null(bundled)) {
          genes <- unique(bundled$gene[bundled$ancestry_class == "shared EUR+EAS"])
          if (length(genes) >= 2) {
            return(list(genes = genes, is_live = FALSE,
              note = sprintf("%d genes - this project's bundled %s genes replicated in both European and East-Asian ancestry GWAS.", length(genes), sex_label)))
          }
        }
      }
      list(genes = character(0), is_live = FALSE,
           note = sprintf("No live %s panel yet - run Cross-Ancestry MR on the currently loaded dataset first.", sex_label))
    }

    panel_fns <- list(same_tissue = same_tissue_panel, cross_tissue = cross_tissue_panel, cross_ancestry = cross_ancestry_panel)
    panel_labels <- c(
      same_tissue = "Same-tissue panel (blood train + test)",
      cross_tissue = "Cross-tissue panel (blood-train, synovium-test)",
      cross_ancestry = "Cross-ancestry panel"
    )
    panel_blurbs <- c(
      same_tissue = "Three-method consensus (LASSO ∩ Random Forest ∩ SVM-RFE), MHC-excluded.",
      cross_tissue = "The panel carried from blood training into this project's own synovial (GSE89408) validation, MHC-retained.",
      cross_ancestry = "Genes whose Mendelian-randomisation causal signal replicates in both European and East-Asian ancestry GWAS."
    )

    ## Dynamic panel above the always-present gene_list textarea; avoids the
    ## race of updating an input just removed/re-added to the DOM.
    output$preset_panel_ui <- renderUI({
      src <- input$gene_source %||% "own"
      if (identical(src, "own")) return(NULL)

      sex_label <- input$panel_sex %||% "female"
      panel <- panel_fns[[src]](sex_label)

      div(class = "preset-panel",
        div(class = "card-title", style = "font-size: 13.5px;", panel_labels[[src]]),
        div(class = "card-subtitle", panel_blurbs[[src]]),
        shinyWidgets::radioGroupButtons(
          ns("panel_sex"), NULL, choices = c("Pooled (all)" = "pooled", "Female" = "female", "Male" = "male"),
          selected = sex_label, status = "default", size = "xs", justified = TRUE
        ),
        if (!length(panel$genes)) {
          div(class = "empty-note", style = "margin-top: 10px;", icon("circle-info"), panel$note)
        } else {
          div(class = "gene-chip-row", lapply(panel$genes, function(g) span(class = "gene-chip", g)))
        },
        div(style = "margin-top: 10px;",
            actionButton(ns("load_panel"), "Load these genes into the list below", icon = icon("arrow-down"), class = "btn-sm"))
      )
    })

    observeEvent(input$load_panel, {
      src <- input$gene_source
      req(src %in% names(panel_fns))
      p <- panel_fns[[src]](input$panel_sex %||% "female")
      validate(need(length(p$genes) > 0, "No genes available for this panel/sex."))
      updateTextAreaInput(session, "gene_list", value = paste(p$genes, collapse = "\n"))
    })

    ontology_labels <- c(GO_BP = "GO Biological Process", KEGG = "KEGG pathways", REACTOME = "Reactome pathways")

    result <- eventReactive(input$run_btn, {
      genes_raw <- trimws(unlist(strsplit(input$gene_list, "[,\n\t ]+")))
      genes_raw <- genes_raw[nzchar(genes_raw)]
      validate(need(length(genes_raw) > 0, "Paste gene symbols below (one per line or comma separated), or load a panel above."))

      ## Normalize to uppercase HGNC and dedupe; genes_raw is kept for the submitted-count report.
      genes <- unique(toupper(genes_raw))
      validate(need(length(genes) >= 3, sprintf(
        "Provide at least 3 distinct gene symbols (got %d after removing duplicates).", length(genes)
      )))

      ## Map symbols via this app's shared HGNC harmonizer (cx_harmonize_gene_ids(),
      ## also used by Multi-Omics Pathways) so aliases resolve and unmatched IDs are reported.
      harm <- cx_harmonize_gene_ids(genes)
      validate(need(isTRUE(harm$ok), harm$error %||% "Gene identifier mapping is unavailable in this deployment."))
      mapped_mask <- harm$df$match_type %in% c("exact_symbol", "exact_entrez", "exact_ensembl", "alias_resolved")
      unmapped_genes <- sort(harm$df$input_id[!mapped_mask])
      mapped_df <- harm$df[mapped_mask, , drop = FALSE]
      validate(need(nrow(mapped_df) >= 3, sprintf(
        "Fewer than 3 of the %d submitted gene symbol(s) could be mapped to a valid Entrez gene ID (%d unmapped: %s).",
        length(genes), length(unmapped_genes), paste(utils::head(unmapped_genes, 10), collapse = ", ")
      )))

      ## Background/universe is genes measured in the loaded dataset, not the whole genome.
      universe_symbols <- rownames(dataset$expr)
      map <- suppressMessages(AnnotationDbi::select(
        org.Hs.eg.db, keys = universe_symbols, keytype = "SYMBOL", columns = "ENTREZID"
      ))
      map <- map[!is.na(map$ENTREZID) & !duplicated(map$SYMBOL), ]
      universe_entrez <- unique(map$ENTREZID)
      validate(need(length(universe_entrez) > 0, "The currently loaded dataset has no genes that map to an Entrez ID - cannot build a background/universe."))

      mapped_df$in_universe <- mapped_df$entrez_id %in% universe_entrez
      excluded_not_measured <- sort(unique(mapped_df$input_id[!mapped_df$in_universe]))
      gene_entrez <- unique(mapped_df$entrez_id[mapped_df$in_universe])
      validate(need(length(gene_entrez) >= 3, sprintf(
        "Fewer than 3 mapped genes are measured in the currently loaded dataset (the enrichment background) - %d mapped gene(s) excluded as not measured here: %s.",
        length(excluded_not_measured), paste(utils::head(excluded_not_measured, 10), collapse = ", ")
      )))

      db_label <- ontology_labels[[input$ontology]] %||% input$ontology

      ego <- if (identical(input$ontology, "GO_BP")) {
        tryCatch(clusterProfiler::enrichGO(
          gene = gene_entrez, universe = universe_entrez, OrgDb = org.Hs.eg.db,
          keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
          pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
        ), error = function(e) e)
      } else if (identical(input$ontology, "KEGG")) {
        tryCatch(clusterProfiler::enrichKEGG(
          gene = gene_entrez, universe = universe_entrez, organism = "hsa",
          pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1
        ), error = function(e) e)
      } else if (!requireNamespace("ReactomePA", quietly = TRUE) || !requireNamespace("reactome.db", quietly = TRUE)) {
        simpleError("Reactome pathway enrichment requires the ReactomePA and reactome.db packages, which are not installed in this deployment.")
      } else {
        tryCatch(ReactomePA::enrichPathway(
          gene = gene_entrez, universe = universe_entrez, organism = "human",
          pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
        ), error = function(e) e)
      }
      ## need() forces its message arg eagerly, so conditionMessage(ego) must stay
      ## conditional here - it has no method for a successful enrichResult object.
      validate(need(!inherits(ego, "error"),
        if (inherits(ego, "error")) sprintf("%s enrichment failed: %s", db_label, conditionMessage(ego)) else NULL
      ))

      ## enrichKEGG() has no `readable` arg, so convert its Entrez-keyed geneID column post hoc.
      if (identical(input$ontology, "KEGG")) {
        ego <- tryCatch(clusterProfiler::setReadable(ego, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"), error = function(e) ego)
      }

      df <- as.data.frame(ego)
      validate(need(nrow(df) > 0, sprintf(
        "%s returned no gene sets for these genes against this background (0 tested) - try a larger gene list or a different database.", db_label
      )))
      df <- df[order(df$pvalue, df$p.adjust, df$qvalue), , drop = FALSE]

      string_deg <- fetch_string_degrees(genes)
      hub_table <- do.call(rbind, lapply(genes, function(g) {
        w <- wgcna_hub_lookup_gated(g)
        sd <- string_deg[[g]]
        data.frame(
          gene = g,
          string_degree = if (!is.null(sd)) sd$degree else 0L,
          string_partners = if (!is.null(sd) && sd$degree > 0) sd$partners else "",
          wgcna_module = w$module, wgcna_disease_module = w$disease_module,
          wgcna_kme = round(w$kme, 3), wgcna_hub = w$is_hub, wgcna_connectivity = round(w$connectivity, 2),
          stringsAsFactors = FALSE
        )
      }))
      hub_table <- hub_table[order(-hub_table$string_degree, -hub_table$wgcna_hub, -hub_table$wgcna_kme), ]

      list(
        table = df,
        database_label = db_label,
        n_submitted = length(genes_raw), n_unique = length(genes),
        n_duplicates_removed = length(genes_raw) - length(genes),
        n_mapped = nrow(mapped_df), n_tested = length(gene_entrez),
        unmapped_genes = unmapped_genes, excluded_not_measured = excluded_not_measured,
        universe_label = sprintf("%s gene(s) measured in the currently loaded dataset", format(length(universe_entrez), big.mark = ",")),
        gene_cards = fetch_gene_cards(genes),
        string_png = fetch_string_network_png(genes),
        hub_table = hub_table
      )
    }, ignoreInit = TRUE)

    enrich_has_run <- reactiveVal(FALSE)
    observeEvent(input$run_btn, enrich_has_run(TRUE), ignoreInit = TRUE)

    ## Result boxes stay out of the DOM until "Run enrichment" is clicked once.
    output$result_box_ui <- renderUI({
      if (!enrich_has_run()) {
        return(box(
          width = NULL, title = tagList(icon("chart-column"), "Result"), status = "primary", solidHeader = FALSE,
          div(class = "empty-note", icon("circle-info"),
              "Not run yet. Enter a gene list on the left, then click \"Run enrichment\".")
        ))
      }
      box(
        width = NULL, title = tagList(icon("chart-column"), "Result"), status = "primary", solidHeader = FALSE,
        withSpinner(uiOutput(ns("summary_ui")), color = "#2c6fbb", type = 6),
        withSpinner(plotOutput(ns("bar_plot"), height = 420), color = "#2c6fbb", type = 6)
      )
    })

    output$below_result_ui <- renderUI({
      req(enrich_has_run())
      tagList(
        box(
          width = 12, title = tagList(icon("table"), "Enriched terms"), status = "primary", solidHeader = FALSE,
          div(class = "table-toolbar", downloadButton(ns("download_enrich"), "Download CSV", class = "btn-sm")),
          DT::dataTableOutput(ns("enrich_table"))
        ),
        fluidRow(
          column(
            7,
            box(
              width = NULL, title = tagList(icon("id-card"), "Gene cards"), status = "primary", solidHeader = FALSE,
              p(class = "card-subtitle", "Live lookup from MyGene.info - the same fields GeneCards.org leads with."),
              withSpinner(uiOutput(ns("gene_cards_ui")), color = "#2c6fbb", type = 6)
            )
          ),
          column(
            5,
            box(
              width = NULL, title = tagList(icon("circle-nodes"), "Protein interaction network"), status = "primary", solidHeader = FALSE,
              p(class = "card-subtitle", "Live network from STRING (string-db.org), confidence ≥ 0.4."),
              div(class = "figure-card", withSpinner(imageOutput(ns("string_network"), height = 300), color = "#2c6fbb", type = 6))
            )
          )
        ),
        box(
          width = 12, title = tagList(icon("share-nodes"), "Hub genes"), status = "primary", solidHeader = FALSE,
          p(class = "card-subtitle",
            "Two independent hub definitions, neither a substitute for the other. ",
            strong("STRING degree"), " counts, live, how many ", em("other genes in this list"),
            " each gene interacts with (PPI confidence ≥ 0.4) - it only sees the genes you entered. ",
            strong("WGCNA hub"), " is this project's own precomputed co-expression network across the full ",
            "15,763-gene training cohort: a gene with intramodular connectivity (kME) above the module average ",
            "in the RA disease-associated module - a property of the gene in the whole transcriptome, independent of your list."),
          withSpinner(DT::dataTableOutput(ns("hub_table_ui")), color = "#2c6fbb", type = 6)
        )
      )
    })

    output$summary_ui <- renderUI({
      if (!enrich_has_run()) {
        return(div(class = "empty-note", icon("circle-info"),
                    "Not run yet. Enter a gene list on the left, then click \"Run enrichment\"."))
      }
      res <- result()
      df <- res$table %>% filter(qvalue < input$qval_cut)
      tagList(
        p(strong(res$n_submitted), " gene(s) submitted",
          if (res$n_duplicates_removed > 0) sprintf(" (%d duplicate(s) removed, %d unique)", res$n_duplicates_removed, res$n_unique) else "",
          "."),
        p(strong(res$n_mapped), " of ", res$n_unique, " unique gene symbol(s) mapped to a valid Entrez gene ID; ",
          strong(res$n_tested), " are measured in the currently loaded dataset and used in the test."),
        if (length(res$unmapped_genes) > 0)
          p(class = "empty-note", icon("triangle-exclamation"),
            sprintf(" %d unmapped/invalid gene symbol(s): ", length(res$unmapped_genes)), paste(res$unmapped_genes, collapse = ", ")),
        if (length(res$excluded_not_measured) > 0)
          p(class = "empty-note", icon("circle-info"),
            sprintf(" %d gene(s) mapped but not measured in the current dataset (excluded from the test): ", length(res$excluded_not_measured)), paste(res$excluded_not_measured, collapse = ", ")),
        p(strong("Database: "), res$database_label, " · ", strong("Q-value cutoff: "), input$qval_cut, " · ", strong("Background/universe: "), res$universe_label),
        p(strong(nrow(df)), " terms significant at q < ", input$qval_cut, " (", nrow(res$table), " tested in total).")
      )
    })

    output$bar_plot <- renderPlot({
      req(enrich_has_run())
      res <- result()
      df <- res$table %>% filter(qvalue < input$qval_cut) %>% arrange(qvalue) %>% head(15)
      validate(need(nrow(df) > 0, "No terms pass the current q-value cutoff."))
      ggplot(df, aes(x = reorder(Description, -log10(qvalue)), y = -log10(qvalue))) +
        geom_col(fill = "#2563EB") +
        coord_flip() +
        labs(x = NULL, y = "-log10 q-value") +
        theme_minimal(base_size = 12)
    })

    output$enrich_table <- DT::renderDataTable({
      req(enrich_has_run())
      res <- result()
      DT::datatable(res$table, rownames = FALSE, filter = "top",
                     options = list(pageLength = 15, scrollX = TRUE), class = "stripe hover compact")
    })

    output$download_enrich <- downloadHandler(
      filename = function() "functional_enrichment.csv",
      content = function(file) {
        res <- result()
        write.csv(res$table, file, row.names = FALSE)
      }
    )

    output$gene_cards_ui <- renderUI({
      req(enrich_has_run())
      cards <- result()$gene_cards
      validate(need(length(cards) > 0, "Gene card lookup unavailable (MyGene.info did not respond) or none of the input genes were found."))
      tagList(lapply(cards, function(g) {
        if (!isTRUE(g$found)) {
          return(div(class = "gene-card-item",
                      div(class = "card-title", g$symbol),
                      div(class = "empty-note", icon("circle-info"), "Not found in MyGene.info.")))
        }
        div(
          class = "gene-card-item",
          div(class = "card-title",
              tags$a(g$symbol, href = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", g$symbol), target = "_blank")),
          div(class = "card-subtitle", g$name %||% ""),
          div(class = "gene-card-chip-row",
              span(class = "gene-card-chip", sprintf("chr%s", g$chr %||% "?")),
              span(class = "gene-card-chip", sprintf("Entrez %s", g$entrez %||% "?")),
              span(class = "gene-card-chip", sprintf("aliases: %s", g$alias %||% "none"))
          ),
          tags$p(class = "gene-card-summary", g$summary %||% "No summary available."),
          tags$a("GeneCards", href = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", g$symbol), target = "_blank", class = "btn btn-xs btn-default"),
          " ",
          tags$a("NCBI Gene", href = paste0("https://www.ncbi.nlm.nih.gov/gene/?term=", g$entrez %||% g$symbol), target = "_blank", class = "btn btn-xs btn-default")
        )
      }))
    })

    output$hub_table_ui <- DT::renderDataTable({
      req(enrich_has_run())
      ht <- result()$hub_table
      validate(need(!is.null(ht) && nrow(ht) > 0, "No hub-gene data for this list."))
      disp <- ht
      names(disp) <- c("Gene", "STRING degree (this list)", "Interacts with", "WGCNA module",
                        "RA disease module?", "WGCNA kME", "WGCNA hub?", "WGCNA connectivity")
      DT::datatable(disp, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "stripe hover compact") %>%
        DT::formatStyle("WGCNA hub?", target = "row",
                         backgroundColor = DT::styleEqual(c(TRUE, FALSE), c("#F0FDF4", "white")))
    })

    output$string_network <- renderImage({
      req(enrich_has_run())
      png_path <- result()$string_png
      validate(need(!is.null(png_path), "STRING network unavailable (string-db.org did not respond) or no interactions were found at this confidence threshold."))
      list(src = png_path, contentType = "image/png", width = "100%", alt = "STRING protein interaction network")
    }, deleteFile = FALSE)
  })
}
