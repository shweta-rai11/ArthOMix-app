## R/shared/mod_arthochat.R
## ArthOChat: app-wide chat assistant (ellmer + shinychat), in its own slide-out drawer.

ARTHOCHAT_MAX_TURNS <- 40L

ARTHOCHAT_MAX_EXECUTIONS <- 5L

## Fixed seed/temperature=0 for deterministic sampling.
ARTHOCHAT_TEMPERATURE <- 0
ARTHOCHAT_SEED <- 20260904L

## ARTHOCHAT_MAX_TURNS above is a per-session counter: a fresh browser tab or
## reload gets a fresh Shiny session and therefore a fresh counter, so on its
## own it does not bound total spend on the metered LLM backend. This adds a
## session-independent budget shared by every session in the R process, plus
## a per-visitor (IP) budget so a single visitor reloading repeatedly cannot
## alone exhaust the global one. Single-threaded per R process, so no locking.
.arthochat_rate_env <- new.env(parent = emptyenv())
ARTHOCHAT_RATE_WINDOW_SECS <- 3600
ARTHOCHAT_GLOBAL_MAX_TURNS_PER_WINDOW <- 500L
ARTHOCHAT_IP_MAX_TURNS_PER_WINDOW <- 60L

.arthochat_rate_reset_if_stale <- function() {
  now <- Sys.time()
  start <- .arthochat_rate_env[["window_start"]]
  if (is.null(start) || as.numeric(difftime(now, start, units = "secs")) >= ARTHOCHAT_RATE_WINDOW_SECS) {
    .arthochat_rate_env[["window_start"]] <- now
    .arthochat_rate_env[["global_count"]] <- 0L
    .arthochat_rate_env[["ip_counts"]] <- list()
  }
}

## Returns TRUE and records the turn if under both budgets; FALSE if the
## caller should refuse it. `ip` may be NULL/"" (e.g. in tests) and is then
## pooled under "unknown", which still counts against the global budget.
arthochat_rate_allow <- function(ip = NULL) {
  .arthochat_rate_reset_if_stale()
  ip <- if (is.null(ip) || !nzchar(ip)) "unknown" else ip

  global_count <- .arthochat_rate_env[["global_count"]] %||% 0L
  if (global_count >= ARTHOCHAT_GLOBAL_MAX_TURNS_PER_WINDOW) return(FALSE)

  ip_counts <- .arthochat_rate_env[["ip_counts"]] %||% list()
  ip_count <- ip_counts[[ip]] %||% 0L
  if (ip_count >= ARTHOCHAT_IP_MAX_TURNS_PER_WINDOW) return(FALSE)

  .arthochat_rate_env[["global_count"]] <- global_count + 1L
  ip_counts[[ip]] <- ip_count + 1L
  .arthochat_rate_env[["ip_counts"]] <- ip_counts
  TRUE
}

## Best-effort client identifier from the Shiny request, for rate limiting
## only (not authentication). Behind a reverse proxy the first hop's address
## is in X-Forwarded-For; falls back to the direct socket address.
.arthochat_client_ip <- function(session) {
  req <- tryCatch(session$request, error = function(e) NULL)
  if (is.null(req)) return(NULL)
  xff <- req$HTTP_X_FORWARDED_FOR
  if (!is.null(xff) && nzchar(xff)) {
    return(trimws(strsplit(xff, ",", fixed = TRUE)[[1]][1]))
  }
  req$REMOTE_ADDR
}

## Live sub-module titles from *_MODULES, not a hardcoded list, to avoid drift.
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

## Classifies each "##/### <header>" section as run/not-run, section by section.
.arthochat_classify_context_modules <- function(context_text, known_modules = .arthochat_known_modules()) {
  empty <- list(not_run = character(0), grounded = character(0))
  if (is.null(context_text) || !nzchar(trimws(context_text %||% ""))) return(empty)
  lines <- strsplit(context_text, "\n", fixed = TRUE)[[1]]
  ## Matches both "## " (vertical) and "### " (per-sub-module) headers.
  header_idx <- grep("^#{2,3} ", lines)
  if (!length(header_idx)) return(empty)

  not_run <- character(0); grounded <- character(0)
  for (i in seq_along(header_idx)) {
    start <- header_idx[i]
    end <- if (i < length(header_idx)) header_idx[i + 1] - 1 else length(lines)
    label_raw <- trimws(sub("^#{2,3} ", "", lines[start]))
    hits <- known_modules[vapply(known_modules, function(m) grepl(tolower(m), tolower(label_raw), fixed = TRUE), logical(1))]
    if (!length(hits)) next
    ## Prefer the longest (most specific) matching title - e.g. a "ML Feature
    ## Selection" header should resolve to itself, not fall back to the
    ## shorter "Feature Selection" title it also happens to contain.
    best_hit <- hits[[which.max(nchar(hits))]]
    section_text <- paste(lines[start:end], collapse = " ")
    if (grepl("not yet run|not yet loaded", section_text, ignore.case = TRUE)) {
      not_run <- c(not_run, best_hit)
    } else {
      grounded <- c(grounded, best_hit)
    }
  }
  list(not_run = unique(not_run), grounded = unique(grounded))
}

## Flags a response naming a sub-module the current context marks not-yet-run, unless it hedges.
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
  ## Sentence-scoped: a hedge only suppresses the module(s) it's actually
  ## said alongside, not every not-run module named anywhere in the response
  ## (a response can legitimately hedge about one module while asserting a
  ## fabricated result for another in the same turn).
  sentences_lower <- tolower(strsplit(response_text, "(?<=[.!?])\\s+", perl = TRUE)[[1]])
  if (!length(sentences_lower)) sentences_lower <- resp_lower

  .escape_regex <- function(x) gsub("([][{}()+*^$|\\.?])", "\\\\\\1", x, perl = TRUE)
  ## Also match a title's stripped form and leading acronym, not just the full title -
  ## on a word boundary, so a short acronym (e.g. "ML") can't substring-match an
  ## unrelated word (e.g. "html", "normally") and false-flag an unrun module.
  .mod_variants <- function(mod) {
    core <- trimws(sub("\\s*\\([^)]*\\)\\s*$", "", mod))
    acronym <- if (grepl("^[A-Z]{2,}\\b", mod)) sub("^([A-Z]{2,})\\b.*", "\\1", mod) else NA_character_
    unique(stats::na.omit(c(mod, core, acronym)))
  }
  .mentions_variant <- function(variants, text_lower) {
    any(vapply(variants, function(v) grepl(paste0("\\b", .escape_regex(tolower(v)), "\\b"), text_lower, perl = TRUE), logical(1)))
  }

  flagged_modules <- Filter(function(mod) {
    variants <- .mod_variants(mod)
    if (!.mentions_variant(variants, resp_lower)) return(FALSE)
    hedged_here <- any(vapply(sentences_lower, function(s) {
      grepl(hedge_pattern, s, perl = TRUE) && .mentions_variant(variants, s)
    }, logical(1)))
    !hedged_here
  }, not_run_modules)
  list(flagged = length(flagged_modules) > 0, modules = unique(flagged_modules))
}

## Code-level consent check for execute_confirmed_run: whether the user's own
## most recent chat message plausibly affirms proceeding, rather than relying
## solely on the model's own judgment of an ambiguous reply.
arthochat_affirms_pending_run <- function(text) {
  if (is.null(text) || !nzchar(trimws(text %||% ""))) return(FALSE)
  t <- tolower(trimws(text))
  negation_pattern <- paste(
    "\\bno\\b", "\\bnope\\b", "don't", "do not", "not now", "\\bwait\\b",
    "\\bcancel\\b", "\\bstop\\b", "hold off", "not yet", "never ?mind",
    sep = "|"
  )
  if (grepl(negation_pattern, t, perl = TRUE)) return(FALSE)
  affirm_pattern <- paste(
    "\\byes\\b", "\\byeah\\b", "\\byep\\b", "\\byup\\b", "\\bsure\\b",
    "\\bok\\b", "\\bokay\\b", "go ahead", "sounds good", "please do",
    "please run", "please proceed", "\\brun it\\b", "\\bdo it\\b",
    "\\bproceed\\b", "\\bconfirm(ed)?\\b", "let's do (it|this)", "go for it",
    sep = "|"
  )
  grepl(affirm_pattern, t, perl = TRUE)
}

## Sub-modules the current context had live session data for, for the transparency footer.
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
  "This session's actual results live ONLY in the \"## Computed analysis",
  "results (this session)\" part of the context below (or in an",
  "other_module_context result). A sub-module block that says \"(not yet run",
  "in this session)\" means exactly that - it has not been run in THIS",
  "session, even if the literature describes what running it typically",
  "produces. Never state, imply, or quote a specific number (a gene count,",
  "DEG count, DMP count, p-value, etc.) as \"this session's\" result unless it",
  "came from that Computed-results block or an other_module_context/",
  "other-module Computed-results block - if the question asks for a live",
  "number and the block says not yet run, the correct answer is that it",
  "hasn't been run yet, never a number borrowed from the literature.",
  "",
  "You have six tools:",
  "",
  "- pubmed_search(query): a live PubMed search for scientific claims or",
  "  methodology questions, or when the user wants more/newer external",
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
          "No AI backend is configured. Set ANTHROPIC_API_KEY, or run a local Ollama server: install Ollama, run \"ollama pull %s\", and make sure it's running at %s. Then reload this page.",
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
    last_user_message <- NULL
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
              padj_cut = padj_cut %||% 0.05, lfc_cut = lfc_cut %||% 0.5
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
            lfc_cut = ellmer::type_number("Absolute log2 fold-change cutoff. Defaults to 0.5.", required = FALSE)
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
              if (!arthochat_affirms_pending_run(last_user_message)) {
                return("The user's most recent message does not contain a clear affirmative reply (e.g. \"yes\", \"go ahead\", \"run it\"). Do not execute - relay the proposal again and wait for explicit agreement, or call cancel_pending_action if they've declined or changed the subject.")
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
        shinychat::chat_append("chat", "You've reached this session's message limit for ArthOChat.")
        return()
      }
      if (!arthochat_rate_allow(.arthochat_client_ip(session))) {
        shinychat::chat_append("chat", "ArthOChat is at its usage limit for this period. Please try again shortly.")
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
      last_user_message <<- input$chat_user_input

      cl <- get_client()
      ctx_text <- system_prompt_r()
      cl$set_system_prompt(ctx_text)

      stream <- cl$stream_async(input$chat_user_input)
      p <- shinychat::chat_append("chat", stream, session = session)

      ## Post-hoc: append a fabrication caveat, or else a "grounded in" footer.
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
