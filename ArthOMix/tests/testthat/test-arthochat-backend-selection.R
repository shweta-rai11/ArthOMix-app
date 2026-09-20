## Backend precedence for ArthOChat: Anthropic key > Hugging Face token (Qwen3-8B on
## Inference Providers, used by the Hugging Face Space) > local Ollama > none.
## Pure environment-variable logic, so no live LLM or network call is made.
source_from_app_root("global.R")

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
