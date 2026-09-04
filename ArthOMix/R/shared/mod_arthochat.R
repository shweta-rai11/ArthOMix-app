## ArthOChat: chat assistant (ellmer + shinychat), living in its own app-wide
## slide-out drawer (ui.R) rather than nested in any one module. Uses a
## hosted Anthropic model when ANTHROPIC_API_KEY is set (arthochat_backend()
## in global.R), falling back to a local Ollama server for offline dev.
## Grounded via four global.R tools: project_methods(),

ARTHOCHAT_MAX_TURNS <- 40L

ARTHOCHAT_MAX_EXECUTIONS <- 5L

## A fixed seed and temperature = 0 make ArthOChat's sampling deterministic
## for a given prompt/context - the verification harness (tests/
## arthochat_verification/README.md) documented the same unanswerable
## question producing two different fabricated numbers on retry when Ollama's
## own (non-zero) default temperature applied. This doesn't fix fabrication
## itself (see arthochat_detect_ungrounded_reference() below for that), but it
## does mean a given session/context no longer has run-to-run answer drift.
ARTHOCHAT_TEMPERATURE <- 0
ARTHOCHAT_SEED <- 20260904L

## Module labels ArthOChat's context builders use as "### <label>" (and, for
## the vertical-level dataset/results headers, "## <label>") section headers
## (see build_tx_context()/build_mx_context()/build_cx_context()/
## build_mo_context() and .format_results_block() in R/modules_index.R).
## Used by arthochat_detect_ungrounded_reference() to recognise when the
## model's response names a specific sub-module, and by
## arthochat_grounded_modules_label() to summarise which ones actually had
## live data this turn.
##
## Computed from the live TX_MODULES/MX_MODULES/CX_MODULES/MULTI_MODULES
## config$title values (R/modules_index.R) rather than hand-maintained as a
## parallel literal list: a hardcoded copy silently drifts out of sync with
## the real titles (mismatched US/UK spelling, renamed titles, punctuation
## differences) and the substring match below then never fires for that
## module - i.e. exactly the bug this guard exists to prevent, just moved
## into the guard's own module list instead of the model's answer. This is a
## function, not a top-level constant, because R/shared/*.R (where this file
## lives) is sourced before R/modules_index.R (0b_load_shared_modules.R runs
## before modules_index.R); it's only ever called at chat-turn time, by which
## point every R/*.R file has been sourced and the *_MODULES lists exist.
.arthochat_known_modules <- function() {
  registries <- list(
    mget("TX_MODULES", envir = .GlobalEnv, ifnotfound = list(NULL))[[1]],
    mget("MX_MODULES", envir = .GlobalEnv, ifnotfound = list(NULL))[[1]],
    mget("CX_MODULES", envir = .GlobalEnv, ifnotfound = list(NULL))[[1]],
    mget("MULTI_MODULES", envir = .GlobalEnv, ifnotfound = list(NULL))[[1]]
  )
  titles <- unlist(lapply(registries, function(reg) {
    if (is.null(reg)) return(character(0))
    vapply(reg, function(m) m$config$title %||% NA_character_, character(1))
  }))
  unique(titles[!is.na(titles) & nzchar(titles)])
}

## Scans a draft assistant response for a mention of a known sub-module that
## the CURRENT context (system prompt, including the "## <module>" sections
## built by build_scoped_assistant_context()) marks as not yet run/loaded in
## this session - i.e. a claim the model could only have produced from its
## own training knowledge, not this session's data, without hedging that it's
## doing so. This is a second, independent safety net alongside the system
## prompt's own "say plainly when unavailable" instruction (which a small
## local model does not reliably follow under pressure - see
## tests/arthochat_verification/README.md) and the context-builder fix in
## R/modules_index.R (which addresses the specific root cause the harness
## found, but not every way a response could still reference an unrun module).
## Splits context_text into "## <header>" sections (the convention every
## context builder in R/modules_index.R uses) and, for each section whose
## header names a known module, classifies it as run/not-run based on
## whether a not-yet-run/loaded marker appears WITHIN THAT SECTION ONLY -
## scoping the check per-section (rather than a flat proximity regex over the
## whole context) so one module's own marker can never be mistaken for a
## different, unrelated module's status just because they're both short
## sections close together in the same context block.
.arthochat_classify_context_modules <- function(context_text, known_modules = .arthochat_known_modules()) {
  empty <- list(not_run = character(0), grounded = character(0))
  if (is.null(context_text) || !nzchar(trimws(context_text %||% ""))) return(empty)
  lines <- strsplit(context_text, "\n", fixed = TRUE)[[1]]
  ## Vertical-level headers use "## ", per-sub-module headers (the ones this
  ## function actually needs to classify) use "### " - matching only "## "
  ## silently never found any per-sub-module header, so this guard never
  ## fired against real context text (see tests/arthochat_verification).
  header_idx <- grep("^#{2,3} ", lines)
  if (!length(header_idx)) return(empty)

  not_run <- character(0); grounded <- character(0)
  for (i in seq_along(header_idx)) {
    start <- header_idx[i]
    end <- if (i < length(header_idx)) header_idx[i + 1] - 1 else length(lines)
    label_raw <- trimws(sub("^#{2,3} ", "", lines[start]))
    hits <- known_modules[vapply(known_modules, function(m) grepl(tolower(m), tolower(label_raw), fixed = TRUE), logical(1))]
    if (!length(hits)) next
    section_text <- paste(lines[start:end], collapse = " ")
    if (grepl("not yet run|not yet loaded", section_text, ignore.case = TRUE)) {
      not_run <- c(not_run, hits[[1]])
    } else {
      grounded <- c(grounded, hits[[1]])
    }
  }
  list(not_run = unique(not_run), grounded = unique(grounded))
}

## Scans a draft assistant response for a mention of a known sub-module that
## the CURRENT context (system prompt, including the "## <module>" sections
## built by build_scoped_assistant_context()) marks as not yet run/loaded in
## this session - i.e. a claim the model could only have produced from its
## own training knowledge, not this session's data, without hedging that it's
## doing so. This is a second, independent safety net alongside the system
## prompt's own "say plainly when unavailable" instruction (which a small
## local model does not reliably follow under pressure - see
## tests/arthochat_verification/README.md) and the context-builder fix in
## R/modules_index.R (which addresses the specific root cause the harness
## found, but not every way a response could still reference an unrun module).
arthochat_detect_ungrounded_reference <- function(response_text, context_text,
                                                   known_modules = .arthochat_known_modules()) {
  empty <- list(flagged = FALSE, modules = character(0))
  if (is.null(response_text) || !nzchar(trimws(response_text %||% ""))) return(empty)
  resp_lower <- tolower(response_text)
  not_run_modules <- .arthochat_classify_context_modules(context_text, known_modules)$not_run
  if (!length(not_run_modules)) return(empty)

  hedge_pattern <- paste(
    "hasn't been run", "has not been run", "not available", "no live result",
    "not yet run", "not yet loaded", "isn't available", "is not available",
    "haven't run", "hasn't run", "no results yet", "not been computed",
    sep = "|"
  )
  hedged <- grepl(hedge_pattern, resp_lower, perl = TRUE)
  ## A response naming a module rarely echoes its full canonical title verbatim
  ## (e.g. it says "WGCNA" or "the DMR analysis", not "WGCNA (Co-Methylation
  ## Network)" or "Differentially Methylated Regions (DMRs)"). Matching only the
  ## full title would let those fabricated answers pass unflagged, so also try
  ## the title with its parenthetical qualifier stripped, and a leading
  ## all-caps acronym if the title has one.
  .mod_variants <- function(mod) {
    core <- trimws(sub("\\s*\\([^)]*\\)\\s*$", "", mod))
    acronym <- if (grepl("^[A-Z]{2,}\\b", mod)) sub("^([A-Z]{2,})\\b.*", "\\1", mod) else NA_character_
    unique(stats::na.omit(c(mod, core, acronym)))
  }
  flagged_modules <- Filter(function(mod) {
    any(vapply(.mod_variants(mod), function(v) grepl(tolower(v), resp_lower, fixed = TRUE), logical(1)))
  }, not_run_modules)
  if (hedged) flagged_modules <- character(0)
  list(flagged = length(flagged_modules) > 0, modules = unique(flagged_modules))
}

## Summarises, for the transparency footer, which known sub-modules the
## CURRENT context actually had live session data for - shown to the user
## only when arthochat_detect_ungrounded_reference() did NOT flag the
## response, so the footer never contradicts a caveat already given.
arthochat_grounded_modules_label <- function(context_text, known_modules = .arthochat_known_modules()) {
  paste(.arthochat_classify_context_modules(context_text, known_modules)$grounded, collapse = ", ")
}

ARTHOCHAT_SYSTEM_PROMPT <- paste(
  "You are ArthOChat, the assistant embedded in the ArthOMix Shiny app",
  "for rheumatoid arthritis multi-omics analysis. Answer anything the user",
  "asks about this project: the currently loaded dataset (the bundled example",
  "cohort, or their own uploaded/merged data - always whatever is actually",
  "loaded right now, described in the context below), this session's analysis",
  "results, how to use or interpret a particular sub-module, or the",
  "underlying biology and methodology behind it. Use the dataset/results",
  "context below and cite specific numbers from it rather than guessing; say",
  "plainly when a sub-module hasn't been run yet instead of inventing",
  "results. You cannot run analyses yourself - if the user needs a result",
  "that isn't in the context, tell them which sub-module to run and what to",
  "set - except Differential Expression, which you CAN run yourself via the",
  "propose_run_dge/execute_confirmed_run tools described below.",
  "",
  "The context below is scoped to whichever module and sub-module the user",
  "currently has open (shown under \"## Current view\") - it refreshes",
  "automatically every time they navigate, so always trust it over anything",
  "said earlier in this conversation about a different view. If a question is",
  "clearly about a different module (e.g. a Methylomics question while",
  "Transcriptomics is open), say the current view doesn't cover that and",
  "either tell them which module to switch to, or call the",
  "other_module_context tool if they explicitly want you to look there",
  "without switching. Never answer using a different module's results than",
  "the ones actually shown to you (in the current context or a tool result) -",
  "if it isn't there, say it's unavailable rather than reusing a stale or",
  "unrelated result.",
  "",
  "Critical distinction, easy to get wrong: project_methods() and",
  "project_methods_methylomics() return the PUBLISHED manuscript's own",
  "write-up and numbers for how that pipeline was originally run - they are",
  "NOT this session's live results, no matter how specific or numeric they",
  "sound. This session's actual results live ONLY in the \"## Computed",
  "analysis results (this session)\" part of the context below (or in an",
  "other_module_context result). A sub-module block that says \"(not yet run",
  "in this session)\" means exactly that - it has not been run in THIS",
  "session, even if the methodology tool or literature describes what",
  "running it typically produces. Never state, imply, or quote a specific",
  "number (a gene count, DEG count, DMP count, p-value, etc.) as \"this",
  "session's\" result unless it came from that Computed-results block or an",
  "other_module_context/other-module Computed-results block - if the",
  "question asks for a live number and the block says not yet run, the",
  "correct answer is that it hasn't been run yet, never a number borrowed",
  "from the methodology write-up or literature.",
  "",
  "You have eight tools:",
  "",
  "- project_methods(module): looks up THIS project's own written methodology",
  "  for a specific transcriptomics sub-module (e.g. \"WGCNA\", \"Mendelian",
  "  randomisation\", \"feature selection\", or a section number like \"2.6\")",
  "  plus its curated reference list - describing how the PUBLISHED pipeline",
  "  was run, not this session's own results. This is the authoritative",
  "  source for \"how does this project do X\" and \"how do I perform or",
  "  interpret module Y\" - use it first whenever the question is about how a",
  "  specific transcriptomics analysis/sub-module works, even if the user",
  "  doesn't name the module explicitly (infer it from what they're asking",
  "  about) - but never for \"what did MY run just produce\", which only the",
  "  Computed-results context below can answer.",
  "- project_methods_methylomics(module): the same idea, but for the",
  "  Methylomics module's own pipeline (e.g. \"DMP\", \"DMR\", \"WGCNA",
  "  methylomics\", \"cell-type deconvolution\", \"feature selection\",",
  "  \"Mendelian randomization\", \"diagnostic classifier\"). Use this instead",
  "  of project_methods whenever the question is clearly about methylation",
  "  data/analysis rather than gene expression.",
  "- pubmed_search(query): a live, broader PubMed search for scientific claims",
  "  project_methods doesn't cover, or when the user wants more/newer external",
  "  literature.",
  "- gwas_catalog_search(query): searches the OpenGWAS catalogue for candidate",
  "  exposure/outcome GWAS datasets matching a trait, tissue or consortium",
  "  name, returning each dataset's OpenGWAS ID, population and sample size.",
  "  Use this whenever the user asks which GWAS/eQTL dataset to use for a",
  "  trait, especially before uploading their own summary statistics on the",
  "  Mendelian Randomization tab. It needs a configured access token; if it",
  "  reports one is missing, relay that plainly and mention the tab's upload",
  "  option as the immediate alternative.",
  "- other_module_context(module): fetches the full context - every",
  "  sub-module, not just the one currently open - for a module OTHER than",
  "  the one shown under \"## Current view\" below (\"transcriptomics\",",
  "  \"methylomics\", \"crossomics\", or \"multiomics\"). Use this only when the",
  "  user explicitly asks about a module they aren't currently viewing; don't",
  "  call it just to pad an answer, and never invent what it would return",
  "  without calling it.",
  "- propose_run_dge(...): proposes a live Differential Expression run",
  "  (contrast_col, ref_group, comp_group, method, plus optional covariate/",
  "  cutoff settings) and returns a plain-language summary of exactly what",
  "  would run - it does NOT run anything yet. Relay that summary to the",
  "  user verbatim (or close to it) and wait for their next message.",
  "- execute_confirmed_run(): actually runs whatever propose_run_dge most",
  "  recently proposed, and reports back real numbers (genes tested,",
  "  significant, up/down, top hits). Only call this if propose_run_dge was",
  "  called earlier in THIS conversation and the user's most recent message",
  "  clearly agrees to proceed - never call it on an assumption, and never",
  "  call it before proposing. If their reply is ambiguous, declines, or",
  "  changes the subject, call cancel_pending_action instead of executing.",
  "- cancel_pending_action(): discards a pending propose_run_dge proposal",
  "  without running it - call this when the user declines or the",
  "  conversation moves on before they confirm.",
  "",
  "For all three: contrast_col/ref_group/comp_group must be values that",
  "actually appear in the dataset context below (e.g. the \"Groups:\" line) -",
  "never invent a column or group name that isn't shown there; if the user",
  "names something not in that context, say so and point them at the Dataset",
  "tab instead of guessing. method must be \"limma\" or \"deseq2\".",
  "",
  "Search before answering, not after. When you cite a paper, give its author,",
  "year, title, journal and PMID (linking to",
  "https://pubmed.ncbi.nlm.nih.gov/PMID/). If a tool returns nothing that",
  "actually supports the claim, say so plainly instead of citing an unrelated",
  "result - never fabricate a paper, PMID, or finding. Plain conversational or",
  "app-navigation questions don't need either tool.",
  "",
  "Keep answers focused: lead with the direct answer, then only the caveats",
  "or citations that actually matter. Two or three references are usually",
  "enough - don't pad the answer with every result a tool returned.",
  sep = "\n"
)

build_arthochat_system_prompt <- function(view, dataset, results,
                                           methyl_dataset, methyl_results,
                                           cross_dataset, cross_results,
                                           multi_dataset, multi_results) {
  ctx <- build_scoped_assistant_context(
    view$module, view$submodule_id,
    dataset, results, methyl_dataset, methyl_results,
    cross_dataset, cross_results, multi_dataset, multi_results
  )
  paste(ARTHOCHAT_SYSTEM_PROMPT, "", sprintf("## Current view: %s", view$view_label), "", ctx, sep = "\n")
}

mod_arthochat_ui <- function(id) {
  ns <- NS(id)
  if (identical(arthochat_backend(), "none")) {
    return(
      div(
        class = "coming-soon",
        icon("robot", class = "coming-soon-icon"),
        h4("ArthOChat isn't reachable"),
        p(sprintf(
          "No AI backend is configured. Either set ANTHROPIC_API_KEY, or run a local Ollama server (install Ollama, run \"ollama pull %s\", make sure it's running at %s), then reload this page.",
          ARTHOMIX_OLLAMA_MODEL, ollama_base_url()
        ))
      )
    )
  }
  tagList(
    p(class = "submodule-desc",
      "Ask about your dataset, results, or the science behind them."),
    shinychat::chat_ui(
      ns("chat"),
      placeholder = "Ask about your dataset, results, or the science behind them...",
      height = "100%", fill = TRUE
    )
  )
}

mod_arthochat_server <- function(id, dataset, results = NULL,
                                  methyl_dataset = NULL, methyl_results = NULL,
                                  cross_dataset = NULL, cross_results = NULL,
                                  multi_dataset = NULL, multi_results = NULL,
                                  current_context = NULL,
                                  run_hooks = new.env(parent = emptyenv())) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    if (identical(arthochat_backend(), "none")) return(invisible(NULL))

    view_r <- if (is.null(current_context)) {
      reactive(list(module = "app", view_label = "ArthOMix", submodule_id = NULL))
    } else {
      current_context
    }

    system_prompt_r <- reactive({
      build_arthochat_system_prompt(
        view_r(), dataset, results,
        methyl_dataset, methyl_results,
        cross_dataset, cross_results,
        multi_dataset, multi_results
      )
    })

    client <- NULL
    get_client <- function() {
      if (is.null(client)) {
        cl <- if (identical(arthochat_backend(), "anthropic")) {
          ellmer::chat_anthropic(
            model = ARTHOCHAT_ANTHROPIC_MODEL,
            system_prompt = system_prompt_r(),
            params = ellmer::params(temperature = ARTHOCHAT_TEMPERATURE)
          )
        } else {
          ellmer::chat_ollama(
            model = ARTHOMIX_OLLAMA_MODEL,
            base_url = ollama_base_url(),
            system_prompt = system_prompt_r(),
            params = ellmer::params(temperature = ARTHOCHAT_TEMPERATURE, seed = ARTHOCHAT_SEED),
            api_args = list(think = FALSE)
          )
        }
        cl$register_tool(ellmer::tool(
          pubmed_search,
          paste(
            "Search PubMed for published literature relevant to a methodology",
            "question (e.g. normalisation, batch correction, a specific",
            "technique or biomarker) and return matching papers with author,",
            "year, title, journal and PMID."
          ),
          arguments = list(
            query = ellmer::type_string("The PubMed search query - keywords, not a full sentence."),
            max_results = ellmer::type_integer("Number of references to return (1-10). Defaults to 5.", required = FALSE)
          )
        ))
        cl$register_tool(ellmer::tool(
          project_methods,
          paste(
            "Look up this project's own methodology write-up and curated",
            "reference list for a specific analysis sub-module - e.g. \"WGCNA\",",
            "\"Mendelian randomisation\", \"feature selection\", \"diagnostic model\",",
            "or a section number like \"2.6\". Use this before pubmed_search",
            "whenever the question is about how a specific sub-module works,",
            "or how to perform or interpret it."
          ),
          arguments = list(
            module = ellmer::type_string("The sub-module name or topic to look up, e.g. \"WGCNA\" or \"2.6\".")
          )
        ))
        cl$register_tool(ellmer::tool(
          project_methods_methylomics,
          paste(
            "Look up the Methylomics module's own methodology write-up and",
            "curated reference list for a specific analysis sub-module - e.g.",
            "\"DMP\", \"DMR\", \"WGCNA methylomics\", \"cell-type deconvolution\",",
            "\"feature selection\", \"Mendelian randomization\", or \"diagnostic",
            "classifier\". Use this (not project_methods) whenever the question",
            "is about methylation data or the Methylomics module specifically."
          ),
          arguments = list(
            module = ellmer::type_string("The Methylomics sub-module name or topic to look up, e.g. \"DMR\" or \"cell-type deconvolution\".")
          )
        ))
        cl$register_tool(ellmer::tool(
          gwas_catalog_search,
          paste(
            "Search the OpenGWAS catalogue for candidate exposure/outcome GWAS",
            "or eQTL datasets matching a trait, tissue, or consortium name.",
            "Returns each match's OpenGWAS ID, population, and sample size -",
            "use before recommending an OpenGWAS ID or before the user",
            "uploads their own summary statistics on the Mendelian",
            "Randomization tab."
          ),
          arguments = list(
            query = ellmer::type_string("Trait, tissue, or consortium keywords, e.g. \"rheumatoid arthritis\" or \"whole blood eQTL\"."),
            max_results = ellmer::type_integer("Number of datasets to return (1-25). Defaults to 10.", required = FALSE)
          )
        ))
        cl$register_tool(ellmer::tool(
          function(module) {
            isolate(build_scoped_assistant_context(
              tolower(trimws(module)), NULL,
              dataset, results, methyl_dataset, methyl_results,
              cross_dataset, cross_results, multi_dataset, multi_results
            ))
          },
          paste(
            "Fetches the full context (every sub-module, not just one) for a",
            "module OTHER than the one currently open. Use only when the user",
            "explicitly asks about a module they aren't currently viewing."
          ),
          arguments = list(
            module = ellmer::type_string(
              "One of: \"transcriptomics\", \"methylomics\", \"crossomics\", \"multiomics\"."
            )
          ),
          name = "other_module_context"
        ))
        cl$register_tool(ellmer::tool(
          function(contrast_col, ref_group, comp_group, method,
                   covariate_col = NULL, covariate_mode = NULL, covariate_level = NULL,
                   padj_cut = NULL, lfc_cut = NULL) {
            params <- list(
              contrast_col = contrast_col, ref_group = ref_group, comp_group = comp_group, method = method,
              covariate_col = covariate_col %||% "(none)", covariate_mode = covariate_mode %||% "filter",
              covariate_level = covariate_level,
              padj_cut = padj_cut %||% 0.05, lfc_cut = lfc_cut %||% 0.1
            )
            summary_txt <- sprintf(
              "Run Differential Expression: %s vs %s on \"%s\"%s, method = %s, adj.P cutoff = %s, |log2FC| cutoff = %s.",
              params$comp_group, params$ref_group, params$contrast_col,
              if (!identical(params$covariate_col, "(none)")) sprintf(", covariate %s (%s)", params$covariate_col, params$covariate_mode) else "",
              params$method, params$padj_cut, params$lfc_cut
            )
            pending_action(list(params = params, summary = summary_txt))
            sprintf("Proposed but NOT run yet. Relay this to the user and wait for their explicit agreement before calling execute_confirmed_run: \"%s\"", summary_txt)
          },
          paste(
            "Proposes a live Differential Expression run with the given",
            "contrast/method (and optional covariate/cutoff settings) and",
            "returns a plain-language summary of exactly what would run -",
            "does NOT run anything. Relay the summary to the user and wait",
            "for their next message before ever calling execute_confirmed_run."
          ),
          arguments = list(
            contrast_col = ellmer::type_string("Metadata column defining the two groups to compare, e.g. \"group\" - must appear in the dataset context."),
            ref_group = ellmer::type_string("Reference (baseline) level of contrast_col."),
            comp_group = ellmer::type_string("Comparison level of contrast_col."),
            method = ellmer::type_string("\"limma\" or \"deseq2\"."),
            covariate_col = ellmer::type_string("Metadata column to filter/adjust for, or omit for none.", required = FALSE),
            covariate_mode = ellmer::type_string("\"filter\" or \"adjust\" - required if covariate_col is set.", required = FALSE),
            covariate_level = ellmer::type_string("Level to filter covariate_col to - required if covariate_mode is \"filter\".", required = FALSE),
            padj_cut = ellmer::type_number("Adjusted p-value cutoff for significance. Defaults to 0.05.", required = FALSE),
            lfc_cut = ellmer::type_number("Absolute log2 fold-change cutoff. Defaults to 0.1.", required = FALSE)
          ),
          name = "propose_run_dge"
        ))
        cl$register_tool(ellmer::tool(
          function() {
            isolate({
              pa <- pending_action()
              if (is.null(pa)) {
                return("Nothing is pending confirmation. Call propose_run_dge first, relay it to the user, and only call this after they clearly agree.")
              }
              if (n_executions() >= ARTHOCHAT_MAX_EXECUTIONS) {
                return("This session's limit for agent-triggered analysis runs has been reached. Tell the user to use the Differential Expression tab directly, or reload the app to reset the limit.")
              }
              run_fn <- run_hooks$transcriptomics$dge
              if (is.null(run_fn)) {
                return("Differential Expression isn't available to run from chat in this deployment.")
              }
              result <- tryCatch(do.call(run_fn, pa$params), error = function(e) conditionMessage(e))
              pending_action(NULL)
              if (is.character(result)) {
                return(paste("Could not run it:", result))
              }
              n_executions(n_executions() + 1L)
              sprintf(
                "Done. %s: %d genes tested, %d significant (%d up, %d down). Top hits: %s. Saved as this session's results - visible in Differential Expression's Result panel and Candidate Gene Identification's contrast picker.",
                result$contrast, result$n_tested, result$n_significant, result$n_up, result$n_down, paste(result$top_hits, collapse = ", ")
              )
            })
          },
          paste(
            "Actually runs whatever propose_run_dge most recently proposed,",
            "and reports back real numbers. Only call this if propose_run_dge",
            "was called earlier in this conversation and the user's most",
            "recent message clearly agrees to proceed - never on an",
            "assumption. Takes no arguments; it runs exactly what was",
            "proposed, not a new guess at the params."
          ),
          arguments = list(),
          name = "execute_confirmed_run"
        ))
        cl$register_tool(ellmer::tool(
          function() {
            pending_action(NULL)
            "Pending run cancelled - nothing was executed."
          },
          paste(
            "Discards a pending propose_run_dge proposal without running it.",
            "Call this when the user declines, or the conversation moves on",
            "before they confirm."
          ),
          arguments = list(),
          name = "cancel_pending_action"
        ))
        client <<- cl
      }
      client
    }

    n_turns <- reactiveVal(0L)
    last_view_key <- NULL

    pending_action <- reactiveVal(NULL)
    n_executions <- reactiveVal(0L)

    observeEvent(input$chat_user_input, {
      if (n_turns() >= ARTHOCHAT_MAX_TURNS) {
        shinychat::chat_append("chat", "You've reached this session's message limit for ArthOChat. Reload the app to reset it.")
        return()
      }
      n_turns(n_turns() + 1L)

      view <- view_r()
      view_key <- paste(view$module, view$submodule_id %||% "")
      if (!is.null(last_view_key) && !identical(last_view_key, view_key)) {
        client <<- NULL
        pending_action(NULL)
        shinychat::chat_clear("chat", session = session)
      }
      last_view_key <<- view_key

      cl <- get_client()
      ctx_text <- system_prompt_r()
      cl$set_system_prompt(ctx_text)

      stream <- cl$stream_async(input$chat_user_input)
      p <- shinychat::chat_append("chat", stream, session = session)

      ## Post-hoc, non-blocking checks once the full response has streamed:
      ## (1) the fabrication guard appends a caveat if the response named a
      ## sub-module this session's context marked as not-yet-run/loaded
      ## without hedging; (2) otherwise, a short transparency footer names
      ## which sub-modules' live data the response could actually draw on.
      promises::then(p, onFulfilled = function(full_text) {
        chk <- arthochat_detect_ungrounded_reference(full_text, ctx_text)
        if (isTRUE(chk$flagged)) {
          shinychat::chat_append(
            "chat",
            sprintf(
              "_Note: I don't have live results for %s in this session - the above may reflect general/training knowledge, not this session's data._",
              paste(chk$modules, collapse = ", ")
            ),
            session = session
          )
        } else {
          grounded <- arthochat_grounded_modules_label(ctx_text)
          if (nzchar(grounded)) {
            shinychat::chat_append("chat", sprintf("_Grounded in: %s_", grounded), session = session)
          }
        }
      })
      ## Error handling: if the backend becomes unreachable mid-session (after
      ## the initial arthochat_backend() check passed) or the stream
      ## otherwise fails, show a clear chat message instead of an
      ## unhandled/silent failure.
      promises::catch(p, function(e) {
        shinychat::chat_append(
          "chat",
          "The AI assistant hit an error while responding - please try again.",
          session = session
        )
      })
    })
  })
}
