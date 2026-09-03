## R/auth/auth_api.R
## Hand-rolled Supabase GoTrue (Auth) REST client, via httr2 - no shinymanager,
## no supabase R package (none reliably on CRAN). Every function here talks

SUPABASE_URL <- Sys.getenv("SUPABASE_URL", "")
SUPABASE_ANON_KEY <- Sys.getenv("SUPABASE_ANON_KEY", "")

supabase_configured <- function() nzchar(SUPABASE_URL) && nzchar(SUPABASE_ANON_KEY)

.supabase_extract_error <- function(resp_data) {
  msg <- resp_data$msg %||% resp_data$error_description %||% resp_data$error %||% NULL
  if (is.null(msg) || !nzchar(msg)) "Something went wrong. Please try again." else msg
}

.supabase_req <- function(path, body = NULL, method = "POST", access_token = NULL, query = NULL) {
  if (!supabase_configured()) {
    return(list(ok = FALSE, status = NA, data = NULL,
                error = "Authentication isn't configured. Set SUPABASE_URL and SUPABASE_ANON_KEY in .Renviron and restart the app."))
  }
  req <- httr2::request(paste0(SUPABASE_URL, path)) |>
    httr2::req_headers(apikey = SUPABASE_ANON_KEY, `Content-Type` = "application/json") |>
    httr2::req_method(method) |>
    httr2::req_timeout(15) |>
    httr2::req_error(is_error = function(resp) FALSE)
  if (!is.null(access_token)) {
    req <- httr2::req_headers(req, Authorization = paste("Bearer", access_token))
  } else {
    req <- httr2::req_headers(req, Authorization = paste("Bearer", SUPABASE_ANON_KEY))
  }
  if (!is.null(query)) req <- httr2::req_url_query(req, !!!query)
  if (!is.null(body)) req <- httr2::req_body_json(req, body)

  resp <- tryCatch(httr2::req_perform(req), error = function(e) NULL)
  if (is.null(resp)) {
    return(list(ok = FALSE, status = NA, data = NULL,
                error = "Couldn't reach the authentication service. Check your connection and try again."))
  }
  status <- httr2::resp_status(resp)
  data <- tryCatch(httr2::resp_body_json(resp), error = function(e) NULL)
  ok <- status >= 200 && status < 300
  list(ok = ok, status = status, data = data,
       error = if (ok) NULL else .supabase_extract_error(data))
}

supabase_sign_up <- function(email, password) {
  .supabase_req("/auth/v1/signup", body = list(email = email, password = password))
}

supabase_sign_in <- function(email, password) {
  res <- .supabase_req("/auth/v1/token", query = list(grant_type = "password"),
                        body = list(email = email, password = password))
  if (!res$ok) res$error <- "Invalid email or password."
  res
}

supabase_verify_otp <- function(token_hash, type) {
  .supabase_req("/auth/v1/verify", body = list(token_hash = token_hash, type = type))
}

supabase_recover <- function(email) {
  .supabase_req("/auth/v1/recover", body = list(email = email))
}

supabase_resend <- function(email, type = "signup") {
  .supabase_req("/auth/v1/resend", body = list(email = email, type = type))
}

supabase_update_password <- function(access_token, new_password) {
  .supabase_req("/auth/v1/user", method = "PUT", access_token = access_token,
                body = list(password = new_password))
}

supabase_sign_out <- function(access_token) {
  tryCatch(.supabase_req("/auth/v1/logout", access_token = access_token, body = list()),
           error = function(e) NULL)
  invisible(NULL)
}
