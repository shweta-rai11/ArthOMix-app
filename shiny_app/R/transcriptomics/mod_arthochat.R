## R/mod_arthochat.R
## ArthOChat: the app's one chat assistant (ellmer + shinychat, talking to a
## local Ollama server - see ARTHOMIX_OLLAMA_MODEL / ollama_available() in
## global.R). Nested inside the Overview and Datasets module (mod_overview.R
## renders mod_arthochat_ui()/calls mod_arthochat_server() directly) rather
## than being its own top-level sub-module, so it doesn't have a
## mod_arthochat_config / TX_MODULES entry. Answers anything about the
## project: the currently loaded dataset (the bundled example cohort, or the
## user's own upload from the Dataset tab - `dataset` is the same shared
## reactiveValues either way, so this always reflects whatever's actually
## loaded), this session's analysis results, how to use a particular
## sub-module, or the underlying biology/methodology. Four tools ground it
## in real sources instead of memory: project_methods() (global.R) surfaces
## this project's own written methodology and curated references per
## transcriptomics sub-module; project_methods_methylomics() (global.R) does
## the same against the methylomics pipeline's own per-stage write-ups (both
## the Methylomics and Transcriptomics sidebars open this one shared drawer -
## see arthochat_shortcut_ui() - rather than Methylomics having its own
## separate scoped session, so this tool is what keeps methylomics
## methodology questions grounded here); pubmed_search() (global.R) does a
## live PubMed lookup for anything broader; gwas_catalog_search() (global.R)
## looks up candidate OpenGWAS exposure/outcome datasets for a trait - e.g.
## "which GWAS should I use for RA in the Mendelian Randomization tab's
## upload mode" - degrading to an instructional message if no OPENGWAS_JWT
## token is configured.

ARTHOCHAT_MAX_TURNS <- 40L

ARTHOCHAT_SYSTEM_PROMPT <- paste(
  "You are ArthOChat, the assistant embedded in the ArthOMix Explorer Shiny app",
  "for rheumatoid arthritis transcriptomics analysis. Answer anything the user",
  "asks about this project: the currently loaded dataset (the bundled example",
  "cohort, or their own uploaded expression matrix and metadata - always",
  "whatever is actually loaded right now, described in the context below),",
  "this session's analysis results, how to use or interpret a particular",
  "sub-module, or the underlying biology and methodology behind it. Use the",
  "dataset/results context below and cite specific numbers from it rather",
  "than guessing; say plainly when a sub-module hasn't been run yet instead",
  "of inventing results. You cannot run analyses yourself - if the user needs",
  "a result that isn't in the context, tell them which sub-module to run and",
  "what to set.",
  "",
  "You have four research tools:",
  "",
  "- project_methods(module): looks up THIS project's own written methodology",
  "  for a specific transcriptomics sub-module (e.g. \"WGCNA\", \"Mendelian",
  "  randomisation\", \"feature selection\", or a section number like \"2.6\")",
  "  plus its curated reference list. This is the authoritative source for",
  "  \"how does this project do X\" and \"how do I perform or interpret module",
  "  Y\" - use it first whenever the question is about a specific",
  "  transcriptomics analysis/sub-module, even if the user doesn't name the",
  "  module explicitly (infer it from what they're asking about).",
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

build_arthochat_system_prompt <- function(dataset, results) {
  paste(ARTHOCHAT_SYSTEM_PROMPT, "", build_assistant_context(dataset, results), sep = "\n")
}

## No page-header h2 here - the drawer's own header (see ui.R) already
## shows "ArthOChat" with a close button; this UI is now the drawer's
## body content, not a full standalone page.
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

mod_arthochat_server <- function(id, dataset, results = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    if (!ollama_available()) return(invisible(NULL))

    client <- NULL
    get_client <- function() {
      if (is.null(client)) {
        cl <- ellmer::chat_ollama(
          model = ARTHOMIX_OLLAMA_MODEL,
          system_prompt = build_arthochat_system_prompt(dataset, results),
          ## qwen3's hybrid reasoning mode is ~15x slower for little benefit on
          ## this task (see mod_assistant.R) - keep it off for a responsive chat.
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
        client <<- cl
      }
      client
    }

    n_turns <- reactiveVal(0L)

    observeEvent(input$chat_user_input, {
      if (n_turns() >= ARTHOCHAT_MAX_TURNS) {
        shinychat::chat_append("chat", "You've reached this session's message limit for ArthOChat. Reload the app to reset it.")
        return()
      }
      n_turns(n_turns() + 1L)

      cl <- get_client()
      cl$set_system_prompt(build_arthochat_system_prompt(dataset, results))

      stream <- cl$stream_async(input$chat_user_input)
      shinychat::chat_append("chat", stream)
    })
  })
}
