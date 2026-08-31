## ArthOChat: chat assistant (ellmer + shinychat) over a local Ollama server,
## living in its own app-wide slide-out drawer (ui.R) rather than nested in
## any one module. Grounded via four global.R tools: project_methods(),
## project_methods_methylomics(), pubmed_search(), gwas_catalog_search() -
## plus one tool defined below, other_module_context(), for an explicit
## cross-module lookup.
##
## Context is module-scoped: `current_context` (passed in from server.R,
## itself a plain reactive() over the existing input$sidebar_tabs/tx_menu/
## mx_menu/cx_menu/mo_menu navigation inputs - no second navigation system)
## tells this module which of the four omics verticals, and which sub-module
## within it, the user is actually looking at right now. system_prompt_r()
## below rebuilds the context block from that plus the matching
## dataset/results reactiveValues via R/submodules_registry.R's
## build_scoped_assistant_context() - a plain reactive(), so Shiny's own
## dependency tracking caches it and only recomputes when the active module/
## sub-module or the data it actually reads changes, not on every chat
## message (see the observeEvent below, which now just reads the cached
## reactive instead of rebuilding the whole context string from scratch).

ARTHOCHAT_MAX_TURNS <- 40L

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
  "set.",
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
  "You have five tools:",
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
  "",
  "Search before answering, not after. When you cite a paper, give its author,",
  "year, title, journal and PMID (linking to",
  "https://pubmed.ncbi.nlm.nih.gov/PMID/). If a tool returns nothing that",
  "actually supports the claim, say so plainly instead of citing an unrelated",
  "result - never fabricate a paper, PMID, or finding. Plain conversational or",
  "app-navigation questions don't need either tool.",
  "",
  "You're running on local hardware with limited generation speed, so keep",
  "answers focused: lead with the direct answer, then only the caveats or",
  "citations that actually matter. Two or three references are usually enough",
  "- don't pad the answer with every result a tool returned.",
  sep = "\n"
)

## Builds the full system prompt for one module-scoped `view` - a list(module=,
## view_label=, submodule_id=) as produced by server.R's `current_context`
## reactive - via R/submodules_registry.R's build_scoped_assistant_context().
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

## Drawer body UI; the drawer's own header (ui.R) already shows the title and close button.
mod_arthochat_ui <- function(id) {
  ns <- NS(id)
  if (!ollama_available()) {
    return(
      div(
        class = "coming-soon",
        icon("robot", class = "coming-soon-icon"),
        h4("Ollama isn't reachable"),
        p(sprintf(
          "Couldn't reach %s at %s. Install Ollama, run \"ollama pull %s\", make sure the Ollama app/server is running, then reload this page.",
          ARTHOMIX_OLLAMA_MODEL, ollama_base_url(), ARTHOMIX_OLLAMA_MODEL
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

## `current_context` is a reactive (server.R's `current_module_context`)
## returning list(module=, view_label=, submodule_id=) for whichever
## sidebar_tabs/tx_menu/mx_menu/cx_menu/mo_menu selection is live right now.
## The methyl_*/cross_*/multi_* reactiveValues are optional so any existing
## call site that only passes dataset/results still works, falling back to
## build_assistant_context()'s whole-app view for every module.
mod_arthochat_server <- function(id, dataset, results = NULL,
                                  methyl_dataset = NULL, methyl_results = NULL,
                                  cross_dataset = NULL, cross_results = NULL,
                                  multi_dataset = NULL, multi_results = NULL,
                                  current_context = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    if (!ollama_available()) return(invisible(NULL))

    ## Falls back to a fixed "whole app" view when no navigation-aware caller
    ## passed current_context (keeps this module usable standalone/in tests).
    view_r <- if (is.null(current_context)) {
      reactive(list(module = "app", view_label = "ArthOMix", submodule_id = NULL))
    } else {
      current_context
    }

    ## Plain reactive(): Shiny caches its value and only recomputes when the
    ## active module/sub-module (view_r()) or the specific dataset/results
    ## fields the matching build_*_context() actually reads change - not on
    ## every chat message, and not because of an unrelated module's results
    ## changing while the user is looking at a different one.
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
        cl <- ellmer::chat_ollama(
          model = ARTHOMIX_OLLAMA_MODEL,
          system_prompt = system_prompt_r(),
          ## qwen3's reasoning mode is ~15x slower with little benefit here (see mod_assistant.R).
          api_args = list(think = FALSE)
        )
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
        ## Explicit escape hatch for a cross-module question - isolate()d
        ## since ellmer may invoke this outside a reactive tick; reads the
        ## same reactiveValues build_arthochat_system_prompt() does, just
        ## unscoped (focus_id = NULL) for whichever module is named.
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
        client <<- cl
      }
      client
    }

    n_turns <- reactiveVal(0L)
    last_module <- NULL

    observeEvent(input$chat_user_input, {
      if (n_turns() >= ARTHOCHAT_MAX_TURNS) {
        shinychat::chat_append("chat", "You've reached this session's message limit for ArthOChat. Reload the app to reset it.")
        return()
      }
      n_turns(n_turns() + 1L)

      view <- view_r()
      ## Crossing into a different top-level module (not just a different
      ## sub-module of the same one) drops the ellmer client so get_client()
      ## rebuilds a fresh conversation grounded only in the new context -
      ## observed live (qwen3, think=FALSE) to sometimes keep answering from
      ## a previous module's turns even after set_system_prompt() below
      ## updates the context text. The visible chat_ui transcript is
      ## untouched (chat_append() below just keeps appending to it), so the
      ## user still sees continuous history - only the model's own backend
      ## turn history resets, per this app's "correct current-module answers
      ## over stale conversation context" priority.
      if (!is.null(last_module) && !identical(last_module, view$module)) {
        client <<- NULL
      }
      last_module <<- view$module

      cl <- get_client()
      cl$set_system_prompt(system_prompt_r())

      stream <- cl$stream_async(input$chat_user_input)
      shinychat::chat_append("chat", stream)
    })
  })
}
