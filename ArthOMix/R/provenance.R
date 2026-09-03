## R/provenance.R
## Shared, module-agnostic provenance-manifest helpers.
##

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
