## R/provenance.R
## Shared, module-agnostic provenance-manifest helpers.
##
## Nothing in this app records what actually went into an analysis result
## (input identity, chosen parameters, package/R versions, random seeds,
## warnings, when it ran) in one machine-readable place - a scientific
## review flagged this gap. Rather than build a bespoke record shape per
## module, this file gives every module the same two building blocks:
##
##   arthomix_provenance_record()          - build the record (a plain list)
##   arthomix_provenance_download_handler()- wrap it as a downloadHandler
##
## This file lives under R/ like every other module file, so Shiny's own
## app-loading mechanism (loadSupport()'s automatic sourcing of R/*.R, the
## same mechanism every mod_*.R file already relies on) picks it up with no
## extra registration - see global.R's comment on data_paths.R for why THAT
## one file needs an explicit source() instead (it must exist before
## anything else in R/ loads, which loadSupport()'s alphabetical/unordered
## sourcing can't guarantee).
##
## Checksum primitive: digest::digest(x, algo = "xxhash64") - the exact call
## shape global.R already uses (twice - get_or_compute_wgcna_blocks() and
## get_or_compute_meth_wgcna_blocks()) as a cache key over large matrices,
## so it's both already a project dependency and already proven fast enough
## at this data scale. Reused as-is here, not reinvented.
##
## Scope: this file only provides the shared helper and is wired into three
## representative modules (mod_dge.R, mod_methyl_dmp.R,
## mod_methyl_featureselection.R) to prove the pattern end-to-end. Wiring a
## "Download analysis record" button into the app's other ~40 download
## handlers is deliberately out of scope here - see the memory/report notes
## for this task for that follow-up.

## Builds one provenance record as a plain list - JSON- and RDS-serializable
## either way, so a module can either hand it straight to saveRDS() (as
## mod_methyl_featureselection.R's fs_model_export() already does for its
## own larger record, this one nested inside it) or pass a closure returning
## it to arthomix_provenance_download_handler() below for a standalone
## per-run "analysis record" download.
##
##   module          - character(1) identifying the calling module, e.g. "mod_dge".
##   checksum_input  - anything digest::digest() accepts (a list is typical:
##                      matrix dimensions/content, group vectors, choices) -
##                      hashed with the same xxhash64 primitive global.R uses
##                      for its WGCNA cache keys. This is a content fingerprint
##                      of the ANALYSIS INPUT, not a full data export.
##   params          - named list of the parameters actually used (thresholds,
##                      chosen columns/levels, method name, etc).
##   seed            - the actual seed value(s) passed into the computation,
##                      or NULL if the computation is deterministic/no seed
##                      was used. Some libraries (mixOmics/SNFtool) take
##                      their own internal seed= and ignore an external
##                      set.seed() - callers must pass whatever seed value
##                      was ACTUALLY threaded into the relevant call(s), not
##                      assume a global set.seed() controls everything.
##   packages        - character vector of package names this run's methods
##                      depend on (module decides which, e.g. only "DESeq2"
##                      when input$method == "deseq2"). Each is resolved to
##                      an installed version via utils::packageVersion(),
##                      guarded by requireNamespace() the same way global.R
##                      already gates optional packages (e.g.
##                      MULTI_DIABLO_LIVE_AVAILABLE) - an uninstalled/
##                      unresolvable package degrades to a note, it never
##                      errors the whole record.
##   extra           - named list for anything module-specific that doesn't
##                      fit params/seed (e.g. dataset$declared_data_type,
##                      warnings/notes strings).
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

## Recursively coerces a provenance record into something jsonlite::toJSON()
## renders sensibly: Sys.time()/Date values are formatted to strings (raw
## POSIXct/Date encode as opaque epoch numbers otherwise), factors drop to
## character, and any other object jsonlite can't reasonably serialize
## (matrices/arrays, S4/R6 objects, closures, etc) is reduced to a short
## descriptive string rather than being allowed to fail toJSON() for the
## whole record or silently balloon the file with a full data dump - this
## record is a manifest of what went into a run, not a data export.
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

## Wraps a zero-argument record-building function (typically a reactive, or
## a closure over one) as a shiny::downloadHandler() that writes it out as
## pretty-printed JSON - the reusable-factory pattern mod_methyl_mr.R's
## make_plot_dl() already uses for its six near-identical plot download
## handlers, applied here so every module wiring this in writes one line
## (output$download_provenance <- arthomix_provenance_download_handler(...))
## instead of duplicating downloadHandler boilerplate.
##
##   record_fn - zero-arg function returning an arthomix_provenance_record() list.
##   base_name - filename prefix; the full filename is "<base_name>_<timestamp>.json".
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
