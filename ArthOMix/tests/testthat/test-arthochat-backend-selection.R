## Backend precedence for ArthOChat: Anthropic key > Hugging Face token (Qwen3-8B on
## Inference Providers, used by the Hugging Face Space) > local Ollama > none.
## Pure environment-variable logic, so no live LLM or network call is made.
source_from_app_root("global.R")
source_from_app_root(file.path("R", "shared", "mod_arthochat.R"))

with_backend_env <- function(code, anthropic = "", hf = "") {
  old <- Sys.getenv(c("ANTHROPIC_API_KEY", "ARTHOCHAT_HF_TOKEN"), unset = NA)
  on.exit(for (nm in names(old)) {
    if (is.na(old[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, stats::setNames(list(old[[nm]]), nm))
  }, add = TRUE)
  Sys.setenv(ANTHROPIC_API_KEY = anthropic, ARTHOCHAT_HF_TOKEN = hf)
  force(code)
}

## Swap ollama_available() in the global env (this project is not a package, so
## testthat::local_mocked_bindings() has no namespace to patch).
with_ollama <- function(value, code) {
  old <- get("ollama_available", envir = globalenv())
  on.exit(assign("ollama_available", old, envir = globalenv()), add = TRUE)
  assign("ollama_available", function() value, envir = globalenv())
  force(code)
}

test_that("a Hugging Face token selects the huggingface backend when no Anthropic key is set", {
  with_backend_env(hf = "hf_dummy", {
    expect_true(hf_available())
    expect_identical(arthochat_backend(), "huggingface")
  })
})

test_that("an Anthropic key still wins over a Hugging Face token", {
  with_backend_env(anthropic = "sk-dummy", hf = "hf_dummy", {
    expect_identical(arthochat_backend(), "anthropic")
  })
})

test_that("the Hugging Face backend outranks local Ollama (Ollama is not even probed)", {
  with_backend_env(hf = "hf_dummy", {
    probed <- FALSE
    old <- get("ollama_available", envir = globalenv())
    on.exit(assign("ollama_available", old, envir = globalenv()), add = TRUE)
    assign("ollama_available", function() { probed <<- TRUE; TRUE }, envir = globalenv())
    expect_identical(arthochat_backend(), "huggingface")
    expect_false(probed)
  })
})

test_that("with neither key set the backend falls through to ollama or none", {
  with_backend_env({
    with_ollama(FALSE, expect_identical(arthochat_backend(), "none"))
    with_ollama(TRUE,  expect_identical(arthochat_backend(), "ollama"))
  })
})

test_that("the pinned Hugging Face model is Qwen3-8B with a provider suffix", {
  expect_match(ARTHOCHAT_HF_MODEL, "^Qwen/Qwen3-8B(:[a-z0-9-]+)?$")
})

## --- Visitor-facing wording (public Hugging Face Space) -------------------------------

with_space_id <- function(value, code) {
  old <- Sys.getenv("SPACE_ID", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("SPACE_ID") else Sys.setenv(SPACE_ID = old), add = TRUE)
  if (is.na(value)) Sys.unsetenv("SPACE_ID") else Sys.setenv(SPACE_ID = value)
  force(code)
}

test_that("arthochat_on_hosted_space() follows the SPACE_ID variable Hugging Face sets", {
  with_space_id("user/ArthOMix", expect_true(arthochat_on_hosted_space()))
  with_space_id(NA,              expect_false(arthochat_on_hosted_space()))
})

test_that("with no backend, a Space visitor is not told to install Ollama or set env vars", {
  with_backend_env({
    with_ollama(FALSE, {
      html <- function() as.character(mod_arthochat_ui("chat"))
      with_space_id("user/ArthOMix", {
        expect_match(html(), "not enabled on this server", fixed = TRUE)
        expect_no_match(html(), "Ollama", fixed = TRUE)
        expect_no_match(html(), "ANTHROPIC_API_KEY", fixed = TRUE)
      })
      ## Local development keeps the setup instructions.
      with_space_id(NA, {
        expect_match(html(), "ollama pull", fixed = TRUE)
        expect_match(html(), "ARTHOCHAT_HF_TOKEN", fixed = TRUE)
      })
    })
  })
})

test_that("the privacy note appears only for hosted backends, not for local Ollama", {
  with_backend_env(hf = "hf_dummy",        expect_match(as.character(arthochat_privacy_note()), "external AI service", fixed = TRUE))
  with_backend_env(anthropic = "sk-dummy", expect_match(as.character(arthochat_privacy_note()), "patient-identifiable", fixed = TRUE))
  with_backend_env(with_ollama(TRUE,  expect_null(arthochat_privacy_note())))
  with_backend_env(with_ollama(FALSE, expect_null(arthochat_privacy_note())))
})

test_that("the Hugging Face system prompt ends with the Qwen3 /no_think switch and keeps the original prompt", {
  out <- arthochat_hf_system_prompt("You are ArthOChat.\n\n## Current view: Home")
  expect_match(out, "^You are ArthOChat\\.", perl = TRUE)
  expect_match(out, "Current view: Home", fixed = TRUE)
  expect_match(out, "/no_think$", perl = TRUE)
})
