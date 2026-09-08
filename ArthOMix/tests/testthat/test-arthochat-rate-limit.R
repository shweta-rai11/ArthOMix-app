## Regression coverage for the ArthOChat rate-limit fix (2026-09-08, security
## audit finding: the 40-turn-per-session cap reset on reload/new tab because
## it lived only in a Shiny-session-scoped reactiveVal, so it did not bound
## total spend on the metered Anthropic backend). These test the pure,
## session-independent helpers directly - no live LLM or Shiny session needed.
source_from_app_root(file.path("R", "shared", "mod_arthochat.R"))

reset_arthochat_rate_state <- function() {
  rm(list = ls(envir = .arthochat_rate_env), envir = .arthochat_rate_env)
}

test_that("arthochat_rate_allow grants turns up to the per-IP budget, then refuses", {
  reset_arthochat_rate_state()
  ip <- "203.0.113.5"
  results <- vapply(seq_len(ARTHOCHAT_IP_MAX_TURNS_PER_WINDOW + 3), function(i) {
    arthochat_rate_allow(ip)
  }, logical(1))

  expect_true(all(results[seq_len(ARTHOCHAT_IP_MAX_TURNS_PER_WINDOW)]))
  expect_false(any(results[(ARTHOCHAT_IP_MAX_TURNS_PER_WINDOW + 1):length(results)]))
})

test_that("the per-IP budget survives a simulated reload (a fresh call with the same IP, no session state)", {
  reset_arthochat_rate_state()
  ip <- "198.51.100.9"
  for (i in seq_len(ARTHOCHAT_IP_MAX_TURNS_PER_WINDOW)) arthochat_rate_allow(ip)

  ## Simulates the exact bypass the audit found: a brand-new Shiny session
  ## (fresh n_turns reactiveVal) from the same visitor. The old code allowed
  ## this unconditionally; the fix must not, because nothing session-scoped
  ## is consulted here at all.
  expect_false(arthochat_rate_allow(ip))
})

test_that("different IPs are budgeted independently until the global cap", {
  reset_arthochat_rate_state()
  expect_true(arthochat_rate_allow("10.0.0.1"))
  expect_true(arthochat_rate_allow("10.0.0.2"))
  expect_true(arthochat_rate_allow("10.0.0.3"))
})

test_that("the global budget refuses once exhausted, even across many distinct IPs", {
  reset_arthochat_rate_state()
  ## Use more IPs than the global cap so no single IP's own budget is what
  ## trips the refusal - only the shared global counter can be responsible.
  n <- ARTHOCHAT_GLOBAL_MAX_TURNS_PER_WINDOW + 5L
  ips <- sprintf("10.1.0.%d", seq_len(n) %% 250)
  results <- vapply(seq_len(n), function(i) arthochat_rate_allow(ips[i]), logical(1))

  expect_equal(sum(results), ARTHOCHAT_GLOBAL_MAX_TURNS_PER_WINDOW)
  expect_false(results[n])
})

test_that("a NULL/empty IP is pooled under a shared bucket rather than erroring or bypassing the cap", {
  reset_arthochat_rate_state()
  expect_true(arthochat_rate_allow(NULL))
  expect_true(arthochat_rate_allow(""))
  for (i in seq_len(ARTHOCHAT_IP_MAX_TURNS_PER_WINDOW - 2)) arthochat_rate_allow(NULL)
  expect_false(arthochat_rate_allow(NULL))
})

test_that(".arthochat_client_ip prefers X-Forwarded-For and falls back to REMOTE_ADDR", {
  fake_session_xff <- list(request = list(
    HTTP_X_FORWARDED_FOR = "203.0.113.7, 10.0.0.1",
    REMOTE_ADDR = "10.0.0.1"
  ))
  expect_equal(.arthochat_client_ip(fake_session_xff), "203.0.113.7")

  fake_session_direct <- list(request = list(
    HTTP_X_FORWARDED_FOR = "",
    REMOTE_ADDR = "192.0.2.55"
  ))
  expect_equal(.arthochat_client_ip(fake_session_direct), "192.0.2.55")

  fake_session_none <- list(request = NULL)
  expect_null(.arthochat_client_ip(fake_session_none))
})

test_that("the rate window resets stale counters instead of accumulating forever", {
  reset_arthochat_rate_state()
  ip <- "172.16.0.1"
  arthochat_rate_allow(ip)
  ## Force the window to look expired without waiting a real hour.
  .arthochat_rate_env[["window_start"]] <- Sys.time() - (ARTHOCHAT_RATE_WINDOW_SECS + 1)

  for (i in seq_len(ARTHOCHAT_IP_MAX_TURNS_PER_WINDOW)) {
    expect_true(arthochat_rate_allow(ip))
  }
  expect_false(arthochat_rate_allow(ip))
})
