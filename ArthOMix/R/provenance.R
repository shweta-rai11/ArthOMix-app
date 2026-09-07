## R/provenance.R
## Shared, module-agnostic provenance-manifest helpers.

arthomix_provenance_record <- function(module, checksum_input, params = list(), seed = NULL,
                                        packages = character(0), extra = list()) {
  checksum <- tryCatch(
    digest::digest(checksum_input, algo = "xxhash64"),
    error = function(e) paste0("unavailable: ", conditionMessage(e))
  )

  packages <- unique(as.character(packages))
  package_versions <- stats::setNames(vector("list", length(packages)), packages)
  for (pkg in packages) {
    package_versions[[pkg]] <- if (requireNamespace(pkg, quietly = TRUE)) {
      tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "version unavailable")
    } else {
      "not installed"
    }
  }

  list(
    schema_version = "1.0",
    module = module,
    run_at = Sys.time(),
    checksum = checksum,
    params = params,
    seed = seed,
    software = list(
      r_version = as.character(getRversion()),
      packages = package_versions
    ),
    extra = extra
  )
}

arthomix_provenance_json_safe <- function(x) {
  if (is.list(x)) {
    return(lapply(x, arthomix_provenance_json_safe))
  }
  if (inherits(x, "POSIXt") || inherits(x, "Date")) {
    return(format(x))
  }
  if (is.factor(x)) {
    return(as.character(x))
  }
  if (is.matrix(x) || is.array(x)) {
    return(sprintf("<%s %s, omitted from provenance record>", class(x)[1], paste(dim(x), collapse = "x")))
  }
  if (is.atomic(x) || is.null(x)) {
    return(x)
  }
  tryCatch(as.character(x), error = function(e) sprintf("<%s, not serializable>", paste(class(x), collapse = "/")))
}

arthomix_provenance_download_handler <- function(record_fn, base_name) {
  shiny::downloadHandler(
    filename = function() sprintf("%s_%s.json", base_name, format(Sys.time(), "%Y%m%d_%H%M%S")),
    content = function(file) {
      rec <- tryCatch(
        record_fn(),
        error = function(e) list(error = sprintf("Could not build the analysis record: %s", conditionMessage(e)))
      )
      safe_rec <- arthomix_provenance_json_safe(rec)
      json <- tryCatch(
        jsonlite::toJSON(safe_rec, pretty = TRUE, auto_unbox = TRUE, force = TRUE),
        error = function(e) jsonlite::toJSON(
          list(error = sprintf("Could not serialize the analysis record: %s", conditionMessage(e))),
          pretty = TRUE, auto_unbox = TRUE
        )
      )
      writeLines(json, file)
    }
  )
}

## ---- Session-level analysis-record store ------------------------------------------------
## Each run's record lands here for the header's "Analysis records" panel to list/export.

.arthomix_session_fallback <- new.env(parent = emptyenv())
ARTHOMIX_DIAGNOSTICS_MAX <- 300L

arthomix_session_store <- function(session = shiny::getDefaultReactiveDomain()) {
  if (is.null(session)) return(.arthomix_session_fallback)
  root <- tryCatch(session$rootScope(), error = function(e) session)
  ud <- tryCatch(root$userData, error = function(e) NULL)
  if (is.null(ud) || !is.environment(ud)) return(.arthomix_session_fallback)
  ud
}

arthomix_session_reactive <- function(name, session = shiny::getDefaultReactiveDomain()) {
  store <- arthomix_session_store(session)
  if (is.null(store[[name]])) store[[name]] <- shiny::reactiveVal(list())
  store[[name]]
}

## Appends one record; with dedupe = TRUE, skips it if identical to the module's last record.
arthomix_provenance_push <- function(record, session = shiny::getDefaultReactiveDomain(), dedupe = FALSE) {
  if (is.null(record) || !is.list(record)) return(invisible(NULL))
  rv <- arthomix_session_reactive("arthomix_provenance", session)
  current <- shiny::isolate(rv())
  if (isTRUE(dedupe) && length(current) > 0) {
    same_module <- Filter(function(r) identical(r$module, record$module), current)
    if (length(same_module) > 0) {
      last <- same_module[[length(same_module)]]
      if (identical(last$checksum, record$checksum) && identical(last$params, record$params)) return(invisible(NULL))
    }
  }
  current[[length(current) + 1L]] <- record
  rv(current)
  invisible(record)
}

arthomix_provenance_records <- function(session = shiny::getDefaultReactiveDomain()) {
  rv <- arthomix_session_reactive("arthomix_provenance", session)
  rv()
}

arthomix_provenance_clear <- function(session = shiny::getDefaultReactiveDomain()) {
  rv <- arthomix_session_reactive("arthomix_provenance", session)
  rv(list())
  invisible(NULL)
}

## One-row-per-run summary of the session log, for the Analysis records table.
arthomix_provenance_summary_table <- function(records) {
  if (length(records) == 0) {
    return(data.frame(`#` = integer(0), Module = character(0), `Run at` = character(0), Parameters = character(0),
                      Checksum = character(0), check.names = FALSE, stringsAsFactors = FALSE))
  }
  fmt_param <- function(v) {
    if (is.null(v)) return("NULL")
    if (is.list(v)) return(sprintf("<%d items>", length(v)))
    v <- as.character(v)
    if (length(v) > 4) v <- c(v[1:4], sprintf("... (%d)", length(v)))
    paste(v, collapse = "/")
  }
  rows <- lapply(seq_along(records), function(i) {
    r <- records[[i]]
    params <- r$params %||% list()
    param_str <- if (length(params) == 0) "" else
      paste(sprintf("%s = %s", names(params), vapply(params, fmt_param, character(1))), collapse = "; ")
    data.frame(`#` = i, Module = r$module %||% "?",
               `Run at` = format(r$run_at %||% Sys.time(), "%Y-%m-%d %H:%M:%S"),
               Parameters = param_str, Checksum = r$checksum %||% "",
               check.names = FALSE, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

## ---- Swallowed-error / muffled-warning diagnostics log -------------------------------------
## Drop-in replacements that log the swallowed error/warning instead of losing it silently.
## Both feed the Diagnostics section of the Analysis records panel.

arthomix_diag_push <- function(kind, message, call = NULL, session = shiny::getDefaultReactiveDomain()) {
  entry <- list(kind = kind, message = as.character(message)[1], call = call, at = Sys.time())
  rv <- tryCatch(arthomix_session_reactive("arthomix_diagnostics", session), error = function(e) NULL)
  if (is.null(rv)) return(invisible(entry))
  current <- tryCatch(shiny::isolate(rv()), error = function(e) list())
  current[[length(current) + 1L]] <- entry
  if (length(current) > ARTHOMIX_DIAGNOSTICS_MAX) current <- utils::tail(current, ARTHOMIX_DIAGNOSTICS_MAX)
  tryCatch(rv(current), error = function(e) NULL)
  invisible(entry)
}

arthomix_diag_entries <- function(session = shiny::getDefaultReactiveDomain()) {
  rv <- arthomix_session_reactive("arthomix_diagnostics", session)
  rv()
}

arthomix_diag_clear <- function(session = shiny::getDefaultReactiveDomain()) {
  rv <- arthomix_session_reactive("arthomix_diagnostics", session)
  rv(list())
  invisible(NULL)
}

arthomix_condition_call_label <- function(cond) {
  cl <- tryCatch(conditionCall(cond), error = function(e) NULL)
  if (is.null(cl)) return(NA_character_)
  txt <- tryCatch(paste(deparse(cl, width.cutoff = 120L), collapse = " "), error = function(e) NA_character_)
  if (is.na(txt)) return(NA_character_)
  if (nchar(txt) > 160) txt <- paste0(substr(txt, 1, 157), "...")
  txt
}

## Drop-in for `error = function(e) NULL` that logs the reason (skips shiny.silent.error).
arthomix_null_on_error <- function(e) {
  if (!inherits(e, "shiny.silent.error")) {
    arthomix_diag_push("error", conditionMessage(e), call = arthomix_condition_call_label(e))
  }
  NULL
}

## Drop-in for suppressWarnings(expr): warnings are still muffled, but each one is logged first.
arthomix_quiet <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    arthomix_diag_push("warning", conditionMessage(w), call = arthomix_condition_call_label(w))
    tryCatch(invokeRestart("muffleWarning"), error = function(e) NULL)
  })
}

## Diagnostics table for the panel: newest first, identical messages collapsed with a count.
arthomix_diag_summary_table <- function(entries) {
  if (length(entries) == 0) {
    return(data.frame(Kind = character(0), Count = integer(0), `Last seen` = character(0), Where = character(0),
                      Message = character(0), check.names = FALSE, stringsAsFactors = FALSE))
  }
  df <- data.frame(
    Kind = vapply(entries, function(x) x$kind %||% "?", character(1)),
    Message = vapply(entries, function(x) x$message %||% "", character(1)),
    Where = vapply(entries, function(x) { v <- x$call; if (is.null(v) || is.na(v)) "" else as.character(v) }, character(1)),
    at = as.POSIXct(vapply(entries, function(x) as.numeric(x$at %||% Sys.time()), numeric(1)), origin = "1970-01-01"),
    stringsAsFactors = FALSE
  )
  key <- paste(df$Kind, df$Message, df$Where, sep = "\r")
  agg <- lapply(split(seq_len(nrow(df)), key), function(idx) {
    data.frame(Kind = df$Kind[idx[1]], Count = length(idx),
               `Last seen` = format(max(df$at[idx]), "%H:%M:%S"),
               Where = df$Where[idx[1]], Message = df$Message[idx[1]],
               last_at = max(df$at[idx]), check.names = FALSE, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, agg)
  out <- out[order(out$last_at, decreasing = TRUE), , drop = FALSE]
  out$last_at <- NULL
  rownames(out) <- NULL
  out
}
