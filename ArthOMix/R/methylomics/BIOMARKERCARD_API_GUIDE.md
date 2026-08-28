# Methylomics Biomarker Card — API Evidence Explorer guide

File: [`R/methylomics/mod_methyl_biomarkercard.R`](mod_methyl_biomarkercard.R) (~3000 lines)
Registered as: `MX_MODULES`, [`R/submodules_registry.R:51`](../submodules_registry.R#L51)
Wired up in: [`server.R`](../../server.R) — `mod_methyl_biomarkercard_server("mx_biomarkercard", methyl_dataset, methyl_results)`

---

## 1. What this document covers

This module started as a single-CpG epigenomic profile page (genomic location, CpG
island/gene context, evidence from the loaded dataset). It has been extended into an
**API-driven Methylomic Biomarker Evidence Explorer**: a gene, gene list, CpG, CpG
list, or existing candidate-biomarker table goes in; real external database records
come back, each with a visible source, a live "Open in {Database}" link, and a
retrieval timestamp — never a fabricated or recreated version of what the database
itself already shows.

This document lists every function **added or materially changed** for that
extension. Functions inherited unchanged from the original single-CpG card (Illumina
manifest annotation, ChAMPdata, TxDb/org.Hs.eg.db gene structure, the preloaded
DMP/DMR lookups, the live-dataset Welch t-test, the plots) are not repeated here —
see the file's own top-of-file comment for those.

---

## 2. Core design

- **Every external client returns the same contract**: `list(ok, <data fields>,
  reason, meta)`. `ok=FALSE` + a human-readable `reason` on any failure (network
  error, timeout, no record found, package missing) — never a partial guess.
- **`meta` is a provenance envelope** (`bc_meta()`) attached to every external
  call, old and new: `source`, `query`, `endpoint`, `retrieved_at`, `status`,
  `n_records`. The **Source & API Information** subtab (`bc_section_source_info()`)
  renders one row per query actually made this session, straight from these
  envelopes — this is what makes the card scientifically auditable (spec §22/§34).
- **Every external lookup is opt-in** (its own button, its own `reactiveVal`,
  reset to `NULL` whenever a new biomarker/panel is generated). No external call
  ever fires as a side effect of "Generate Biomarker Card" — this keeps one
  slow/down database from blocking the rest of the card (spec §24/§32), and
  matches the original single-CpG card's own "Look Up External Evidence" pattern.
- **No response is ever recomputed into a synthetic score.** Evidence-summary
  chips report *retrieval status* (Not yet run / Failed / No results / Results
  found), not a combined "biomarker validated" verdict (spec §29/§30).
- **Never fabricate an API.** Databases without a usable public API (DisGeNET,
  GTEx, MethBank, UCSC live tracks) are either linked out to directly, or listed
  under "Considered but not integrated" with the actual reason — see
  `BC_NOT_INTEGRATED_DBS` and `bc_section_sources()`.

---

## 3. Identifier handling

| Function | Purpose | Input → Output |
|---|---|---|
| `bc_detect_identifier_type(x)` | Regex router: is a submitted token a CpG probe, an Ensembl ID, an Entrez ID, or a gene symbol? | `"cg00000029"` → `"cpg"`; `"ENSG00000012048"` → `"ensembl"`; `"672"` → `"entrez"`; `"BRCA1"` → `"gene_symbol"` |
| `bc_split_tokens(text)` | Splits a pasted list on comma/newline/tab/space, dedupes | `"BRCA1, TP53\ncg00000029"` → `c("BRCA1","TP53","cg00000029")` |
| `bc_resolve_identifiers(raw_ids, array_type)` | Resolves a **mixed** list of genes and CpGs in one pass. CpG-like tokens go through the existing `bc_resolve_cpg()` (Illumina manifest/ChAMP); gene-like tokens are batched through the app's shared harmonizer `cx_harmonize_gene_ids()` (one call, not one per gene). Every submitted identifier appears in the result table, resolved or not — spec §6's "never silently discard an unresolved biomarker" | `c("BRCA1","cg00000029","NOTAREAL123")` → a data.frame with one row per input: `input_id, detected_type, resolved, resolved_id, status_label, annotated_gene` |

**Error handling**: empty/blank input → `ok=FALSE` with a reason, never an empty
table pretending nothing was submitted. Unresolvable tokens get
`status_label = "Unresolved - identifier not recognized"` and stay in the table.

**Caching**: none needed — `cx_harmonize_gene_ids()` and `bc_resolve_cpg()` already
have their own process-wide caches (see the file's existing `.bc_*_cache`
environments).

**Used by**: the Biomarker Overview subtab's "Identifier Resolution" card (single
mode) and the Panel Overview subtab's "Identifier Resolution" card (list mode).

---

## 4. Live external-database clients

Every function below follows the same `httr2::request() %>% req_timeout() %>%
req_error(is_error = \(resp) FALSE) %>% req_perform()`, outer-`tryCatch`, `list(ok,
..., reason, meta)` pattern already established in this file's original EWAS
Catalog/EWAS Atlas/KEGG/Reactome/PubMed clients. Each was verified against a real
identifier (BRCA1 / cg00000029) during implementation — no endpoint here was
guessed from documentation alone.

| Function | Database | Endpoint | What it retrieves | Scientific interpretation |
|---|---|---|---|---|
| `bc_ncbi_gene_summary(entrez)` | NCBI Gene | `eutils.ncbi.nlm.nih.gov/esummary.fcgi?db=gene` | Official name, aliases, cytogenetic location, gene summary text | Live confirmation of gene identity/function, distinct from the locally-bundled org.Hs.eg.db mapping shown alongside it |
| `bc_ensembl_gene_lookup(symbol)` | Ensembl (GRCh37/hg19 archive) | `grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}` | Ensembl Gene ID, biotype, genomic span, strand, description | Cross-database identity check against NCBI; the GRCh37 mirror is used deliberately so coordinates match this app's hg19 build (the default `rest.ensembl.org` serves GRCh38 and would silently mismatch) |
| `bc_ensembl_regulatory_overlap(chr, start, end)` | Ensembl Regulatory Build (GRCh37) | `grch37.rest.ensembl.org/overlap/region/human/{region}?feature=regulatory` | Promoter/enhancer features within ±1kb of **the CpG's own coordinate** (not the gene's) | The methylation-specific "CpG → regulatory region" relationship (spec §14/§27) |
| `bc_wikipathways_pathways_for_gene(entrez)` | WikiPathways | reuses `mp_get_wikipathways_termgene()` (msigdbr cache, `R/multiomics/multiomics_pathway_helpers.R`) | Pathway membership | Curated pathway membership, same contract as the existing KEGG/Reactome clients |
| `bc_opentargets_evidence_for_gene(ensembl)` | Open Targets | `POST api.platform.opentargets.org/api/v4/graphql` | Top associated diseases (aggregates GWAS Catalog + other sources via a curated gene→disease join) + druggability tractability by modality | Primary curated disease-association source on this card |
| `bc_hpa_evidence_for_gene(ensembl)` | Human Protein Atlas | `www.proteinatlas.org/{ensembl}.json` | Baseline tissue/blood-lineage RNA expression specificity | Biological-plausibility check: does this gene's baseline expression pattern match the tissue your methylation data comes from? |
| `bc_gwas_catalog_by_gene(symbol)` | GWAS Catalog (NHGRI-EBI) | `www.ebi.ac.uk/gwas/api/search?q={symbol}` | Trait records whose name text-matches the gene symbol, with study/association counts | **Supplementary** trait-name search, explicitly labeled as such — not a curated gene→trait join (the fully denormalized GWAS Catalog REST chain needs 3 sequential calls per SNP; Open Targets above is the primary curated source, see §6 below) |
| `bc_encode_search(query)` | ENCODE | `www.encodeproject.org/search/?type=Experiment` | Regulatory/epigenomic experiments (ChIP-seq, DNase-seq, ATAC-seq, etc.) for a gene | Doubles as Regulatory/Epigenomics evidence and an External Dataset source |
| `bc_geo_search(query)` | GEO (NCBI E-utilities, `db=gds`) | `eutils.ncbi.nlm.nih.gov/esearch.fcgi` → `esummary.fcgi` | Series/dataset accession, title, data type, organism, sample count, date | **Search/metadata only — never downloads a dataset automatically** (spec §16/§17) |
| `bc_biostudies_search(query)` | BioStudies / ArrayExpress (EBI) | `www.ebi.ac.uk/biostudies/api/v1/search` | Accession, title, record type, release date | Independent-dataset discovery, covers ArrayExpress accessions too |
| `bc_literature_search(query)` | PubMed (NCBI E-utilities) | `esearch.fcgi` → `esummary.fcgi`, `db=pubmed` | Title, authors, journal, year, PMID, publication type | Structured (table-ready) search, distinct from the file's pre-existing `bc_pubmed_summaries()` which is keyed by already-known PMIDs, not a free-text query |
| `bc_literature_query(id, preset)` | — (local) | — | Builds the actual query string from a preset (`BC_LITERATURE_PRESETS`: gene only / +methylation / +CpG / +epigenetic biomarker / +rheumatoid arthritis) or free-text disease context | Matches spec §18's requested query combinations |
| `bc_literature_classify(title, abstract)` | — (local) | — | Keyword-heuristic tags: EWAS / DNA methylation / Biomarker study / Disease association / Functional / Mechanistic / Review / Clinical / Validation study | **Always rendered as "Classification (automated)"** in the UI — never presented as a database-provided fact (spec §19) |

### Why GWAS Catalog is a "supplementary" source, not primary

The GWAS Catalog's transactional REST API only exposes a fully denormalized
gene→SNP→association→trait chain (verified directly during implementation): one
call for the gene's SNPs, then one *more* call per SNP for its associations, then
one *more* call per association for its resolved trait name. Doing that for a gene
with dozens of SNPs would mean dozens of sequential live requests from a single
button click — exactly what spec §25 says to avoid. `bc_gwas_catalog_by_gene()`
instead uses GWAS Catalog's own unified trait-search index (`gwas/api/search`, one
call), which returns real, live GWAS Catalog trait/study/association counts, but
matched by trait-name text rather than a curated gene column. This is disclosed
directly in the UI copy next to the table. Open Targets (`bc_opentargets_evidence_for_gene()`)
is the card's primary curated genetic-association source; it aggregates GWAS
Catalog server-side via a proper gene→disease join.

---

## 5. Not integrated (disclosed, not omitted)

Per spec §3/§35, a database without a usable public API is disclosed with its
actual reason rather than silently left out:

- **DisGeNET** — its current API requires registration/licensing not configured in
  this deployment (same decision already made for the sibling transcriptomics
  Biomarker Card).
- **GTEx** — would duplicate the baseline tissue/blood expression already shown via
  Human Protein Atlas above.
- **UCSC live track API** — no live UCSC call is made anywhere in this codebase
  (all existing UCSC usage is the bundled `cytoBandIdeo`/`TxDb.Hsapiens.UCSC.hg19.knownGene`
  data). Instead, `bc_ucsc_browser_link()` builds a direct deep link to the real
  UCSC Genome Browser at the CpG's exact hg19 coordinate — a real external
  visualization, linked rather than redrawn (spec §28).
- **MethBank** — unchanged from the original card: no public API, link-only.

See `BC_NOT_INTEGRATED_DBS` and `bc_section_sources()` for the full, currently-shown list.

---

## 6. Provenance, evidence-status, and comparison machinery

| Function | Purpose |
|---|---|
| `bc_meta(source, query, endpoint, status, n_records)` | Builds the provenance envelope every client attaches as `$meta` |
| `bc_collect_meta_rows(ext_all)` | Walks a merged evidence bundle and extracts every `$meta` into one audit table |
| `bc_section_source_info(ext_all)` | Renders that audit table (the "Source & API Information" subtab) |
| `bc_db_provenance(api_domain, live_url, live_label)` | The "Live query to X — nothing here is precomputed" banner + "Open in X" button shown on every external-DB card |
| `bc_evidence_status(res, field)` | Classifies one database's result into exactly 4 states: Not yet run / Failed / No results / Results found — a checkmark never appears just because a database exists |
| `BC_EVIDENCE_DBS` | Registry of `(key, label, field)` driving the generalized evidence-status chip strip and DB comparison table across all 15 external databases on this card |
| `bc_section_db_evidence_summary(ext_all)` | The chip-strip card ("Evidence Summary (Database Coverage)") |
| `bc_section_db_comparison(ext_all)` | Side-by-side `Resource / Evidence type / Results / Significance / Source` table, explicit about which numbers are p-values vs. association/interaction scores vs. "not applicable" |

---

## 7. Panel mode (multi-gene / multi-CpG list)

| Function | Purpose |
|---|---|
| `bc_aggregate_convergence(long_df, item_col)` | Generic long→wide collapse: groups by the shared item (Disease/Trait/Pathway) and lists every biomarker that hit it, sorted by convergence count |
| `bc_panel_gene_rows(gene_symbols)` / `bc_panel_cpg_rows(cpg_ids, array_type)` | One identity/annotation row per resolved gene or CpG, reusing `bc_gene_structure()`/`bc_resolve_cpg()` unmodified (looped, not reimplemented) |
| `bc_panel_disease_convergence_genes(genes)` | Loops `bc_opentargets_evidence_for_gene()` (capped to 15 genes, spec §25) and aggregates |
| `bc_panel_ewas_trait_convergence_cpgs(cpgs)` | Loops `bc_ewascatalog_query()` (capped to 15 CpGs) and aggregates — the methylation-specific analogue for a CpG list |
| `bc_panel_pathway_convergence_genes(genes)` | Loops the existing KEGG/Reactome clients and aggregates |

**Example**: pasting `BRCA1, BRCA2, TP53` and clicking "Look Up Disease/Trait
Convergence" queries Open Targets for each (capped at 15), then shows which
diseases are shared across more than one of the three genes, with a real count —
never a synthetic combined score.

**Caching**: none beyond the underlying single-item clients' own — panel mode is a
loop over already-fail-soft functions, so one gene/CpG failing never aborts the
whole panel.

---

## 8. Testing

Pure-logic functions (identifier detection/resolution, the provenance envelope,
evidence-status classification, convergence aggregation, literature
classification/query-building, panel row builders) are covered by
[`tests/testthat/test-methyl-biomarkercard-evidence-explorer.R`](../../tests/testthat/test-methyl-biomarkercard-evidence-explorer.R) —
no live network call is made in any test, matching the convention already
established by `tests/testthat/test-biomarkercard-panel.R` for the sibling
transcriptomics card.
