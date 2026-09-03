# ArthOChat sub-module grounding verification

A live, manual test run against the real, running local model — **not** part of
the automated `testthat` suite (it requires a live Ollama server with
`ARTHOMIX_OLLAMA_MODEL` pulled and is non-deterministic, so it must never run
in CI). Run by hand:

```
Rscript tests/arthochat_verification/test_submodule_grounding.R
```

## Method: "known answer in, verify it comes back out"

For every one of ArthOChat's 41 sub-module ids across all four verticals, this
script:

1. Injects a synthetic result containing a `check_value` number this script
   invented itself (e.g. Transcriptomics `dge` gets `1003`), into the exact
   same `results`/`methyl_results`/`cross_results`/`multi_results` shape the
   real app uses.
2. Calls the app's own `build_scoped_assistant_context()` →
   `build_tx/mx/cx/mo_context(focus_id = id)` — the identical code path
   `mod_arthochat.R` uses — to build a real, scoped system prompt.
3. Sends a real request to the live model (`qwen3:8b` via Ollama's
   `/api/chat`, `think = FALSE`, matching the app's own config) asking it to
   report that value back.
4. Grades the response by exact string match against the known, self-chosen
   ground truth — not a subjective read of the model's prose — checking two
   things:
   - **correct**: does the response contain *this* sub-module's own
     `check_value`?
   - **isolation_ok**: does it *not* also contain a neighboring sub-module's
     `check_value` (proof the `focus_id` scoping actually narrowed the
     context, not just that the right number happened to be findable
     somewhere)?
5. Times every call using Ollama's own `total_duration` response field, not
   wall-clock guessing.

Because every ground-truth value was invented by this script rather than
being a real biological result, grading has no ambiguity — a substring match
is a fully deterministic pass/fail. This same "inject a known synthetic
value into the results object, then verify the model's answer contains it"
pattern is reusable for verifying any future ArthOChat change without needing
to read every response by eye.

## Results (`results.csv`), one run, 2026-09-02

| | |
|---|---|
| Sub-modules tested | 41 (16 Transcriptomics, 14 Methylomics, 3 Cross-Omics, 8 Multi-Omics) |
| Correct | 40 / 41 |
| Isolation held | 41 / 41 — no sub-module ever saw a neighbor's value |
| Latency | mean 3.00s, median 2.57s, min 2.24s, max 16.89s (one cold-start outlier; all others warm) |

## The one failure — diagnosed, not just logged

`crossomics / integration` returned a fabricated `0.000123` instead of the
real `5001`. Investigated directly rather than assumed:

- The real value **was** correctly present in the model's context, as
  `"Genes analyzed: 5,001"` — confirmed by calling `build_cx_context()`
  directly and inspecting the generated text.
- Re-asking the identical context with a properly-matched question — *"How
  many genes were analyzed?"* — got the correct answer: *"In this
  integration, 5,001 genes were analyzed."*
- Re-asking the original mismatched question a second time produced a
  **different** fabricated number (`0.543` on the retry, `0.000123` the
  first time) — confirming this is invention, not a misread of a real value.

**Root cause: the test question, not the app.** Cross-Omics' `integration`
id is the one sub-module rendered by a bespoke formatter
(`.format_cx_integration()` in `submodules_registry.R`) instead of the
generic key/value block every other sub-module uses — so it has no field
literally named `check_value`, unlike the other 40. Asking for a field name
that was never shown to it is not a fair grounding test for this one
sub-module; the actual grounding mechanism (does the real vertical/sub-module
data reach the model correctly) held here just as it did everywhere else.

**What this does still show, and is worth keeping**: when asked for
information that plainly doesn't exist in its given context, this local
model sometimes fabricates a plausible-looking number instead of saying it
isn't there — a real, reproducible reliability limitation, just not a
context-scoping bug. Worth citing separately from the 40/41 grounding result.

## Files here

- `test_submodule_grounding.R` — the 1-question-per-sub-module script above,
  safe to re-run any time Ollama is up.
- `results.csv` — its output: one row per sub-module.
- `test_205_protocol.R` — the 5-question-per-sub-module protocol below.
- `results_205_protocol.csv` — its output: 205 rows (41 sub-modules x 5
  questions).

---

# Extended protocol: 5 questions x 41 sub-modules = 205

A researcher asking ArthOChat about a sub-module doesn't only ask "what's the
number" — they ask basic conceptual questions, whether something has been
run at all, they expect it to combine numbers, and they sometimes ask about
the wrong sub-module by mistake. This protocol tests all four of those, plus
one purely qualitative question, applied identically to every one of the 41
sub-modules.

**Scope, stated plainly**: this exercises the real answer-generation
pipeline end-to-end — the same `build_scoped_assistant_context()` →
system prompt → live `qwen3:8b` call `mod_arthochat.R` itself makes, with
`think = FALSE` matching the app's own config. It does **not** drive the
actual browser/Shiny UI (clicking tabs, typing into the chat widget) — that
layer only decides *which* system prompt gets built, and that mapping was
separately verified in `test_submodule_grounding.R` and the `title_to_module_id()`
logic. Calling this "full browser end-to-end" would overstate what was run;
calling it "the real grounding and generation pipeline, 205 times" is exactly
what it is.

## The 5-question protocol

| # | Type | Question (templated) | Tests | Graded? |
|---|---|---|---|---|
| Q1 | Basic / conceptual | *"In plain terms, what does the '[title]' analysis do in this app, and what kind of output does it produce?"* | Open-ended domain explanation | **No** — logged verbatim for human/domain-expert review. No automatic correctness oracle exists for open-ended text, and fabricating one would be dishonest. A weak keyword-presence signal is logged for information only. |
| Q2 | Data scope | *"Has the '[title]' analysis been run in this session, and what dataset is currently loaded?"* | Does it correctly report real dataset identity + run status | Yes — checks the real injected dataset-source string appears |
| Q3 | Numeric, basic | *"Report the exact check_value number recorded for this sub-module's results, and nothing else's."* | Same grounding check as the first script, single value | Yes — exact substring match |
| Q4 | Numeric, advanced | *"This sub-module's results record a check_value and a secondary_metric. What is the sum of the two?"* | Arithmetic reasoning **over** grounded numbers, not just lookup | Yes — exact substring match of the correct sum |
| Q5 | Boundary / adversarial | *"How does the '[a different real sub-module title]' sub-module's result compare to what's shown here?"* | Actively invites fabrication of a value it was never shown — a harder, adversarial version of the isolation check | Yes — fails only if the neighbor's actual value appears as a stated fact |

## Results, one run, 2026-09-02

| | |
|---|---|
| Total questions | 205 (41 sub-modules x 5) |
| Objectively graded (Q2-Q5) | 164 |
| Overall accuracy | **160 / 164 (97.6%)** |
| By vertical | Transcriptomics 64/64 (100%) · Methylomics 56/56 (100%) · Multi-Omics 32/32 (100%) · Cross-Omics 8/12 (66.7%) |
| By question type | Q2 39/41 (95.1%) · Q3 40/41 (97.6%) · Q4 40/41 (97.6%) · Q5 41/41 (100%) |
| Latency, mean (by type) | Q1 35.1s · Q2 10.0s · Q3 3.5s · Q4 5.3s · Q5 22.1s |
| Latency, overall | mean 15.19s, median 9.50s, min 2.41s, max 84.27s |

Open-ended questions (Q1, Q5) cost far more than grounded lookups (Q3): the
model writes a full explanatory paragraph rather than a single number, which
is the honest reason the mean (15.19s) sits well above the median (9.50s) —
a small number of long, discursive answers pull the mean up.

## The 4 failures — diagnosed, not just counted

All four are Cross-Omics, and all four trace back to two already-understood,
distinct root causes — not four separate bugs:

**Cause 1 — `integration`'s bespoke formatter uses different field names.**
Same issue as the single-question test above: `.format_cx_integration()`
renders `"Genes analyzed: 5,001"` / `"Significant DEGs: 10,009"`, never the
literal words `check_value`/`secondary_metric` every other sub-module's
generic block uses. Result: inconsistent handling of the same ambiguity —
on Q3 the model picked the *wrong* field confidently (`10,009`, actually the
DEG count, not genes analyzed), while on Q4 it correctly declined
("check_value and secondary_metric are not explicitly provided... please
provide the specific values"). One fabricated wrong answer, one honest
refusal, same underlying cause.

**Cause 2 — the generic Cross-Omics formatter never renders a dataset line at all.**
Checked directly in `submodules_registry.R`: `build_cx_context()` has no
"currently loaded dataset" text for any id except `integration` (whose
bespoke formatter has its own `Dataset: ...` line). For `biomarkerconv` and
`mrstage`, Q2 asks a question the context structurally cannot answer. Both
times, the model did not say "that information isn't shown here" — it
**fabricated** a specific, plausible-sounding false answer: *"the bundled
example cohort, which includes 24 samples"* (2/2 times, nearly word-for-word
the same fabrication both times). This is the more serious of the two
findings: not a context-scoping failure (the mechanism correctly showed it
nothing, because there is nothing to show), but a real, reproducible
tendency to invent a specific answer rather than say "I don't know" when a
field is entirely absent rather than just unlabeled.

**What holds up**: the actual retrieval-grounding mechanism (Fig. 3/7 in the
main walkthrough) was not what failed here — Cross-Omics' own generic
sub-modules (`biomarkerconv`, `mrstage`) scored correctly on every other
question type, and `integration`'s real numbers were previously confirmed
retrievable with a correctly-worded question. What this run adds is a
distinct, citable finding about **refusal behavior**: this local model does
not reliably decline when the requested information is genuinely absent from
context — sometimes it does (Q4/`integration`), sometimes it fabricates
(Q2/both generic CX ids, Q3/`integration`) — an inconsistency worth reporting
as a limitation in its own right, separate from the grounding-architecture
result.

---

# Why DGE is the only sub-module with an execution risk to test

A natural question after finding a tool-argument casing bug on a smaller
model (see the main walkthrough's Measured Performance section) is whether
the same risk exists for other sub-modules and simply wasn't tested there
too. Checked directly, not assumed - it doesn't, and here's why.

## The other three tools are built to be case-insensitive - confirmed in the source

```r
# global.R - project_methods() and project_methods_methylomics()
q <- tolower(trimws(module %||% ""))          # normalizes BEFORE any comparison

# global.R - gwas_catalog_search()
grepl(query, cat_df$trait, ignore.case = TRUE)      # explicit, on every field matched
grepl(query, cat_df$consortium, ignore.case = TRUE)
grepl(query, cat_df$author, ignore.case = TRUE)

# global.R - pubmed_search()
# forwards free text straight to NCBI's search engine - no exact-match casing
# requirement exists at all
```

## DGE's execution path is case-sensitive - confirmed live, with the real function, same data

```
run_dge_now(contrast_col="group", ref_group="HC", comp_group="RA", method="limma")
  -> 50 genes tested, 1 significant, real result returned

run_dge_now(contrast_col="group", ref_group="Hc", comp_group="Ra", method="limma")
  -> ERROR: "Fewer than 6 samples match this contrast (and filter, if set);
             pick a different combination."
```

Identical function, identical dataset - only the casing of the two argument
strings changed, and the analysis went from a real result to a hard failure.
`meta[[contrast_col]] %in% c(ref_group, comp_group)` in `mod_dge.R` is an
exact R string comparison; `"Hc"` never equals `"HC"`, so zero samples match
and nothing can be fit.

## Conclusion

This is not uneven test coverage - it's a structural fact about the app.
DGE is the **only** sub-module, across all 41, where a model's own word
choice becomes a literal string fed into exact-match code that touches live
data. Every other sub-module is either read-only (grounded Q&A, no
execution) or routes through a tool specifically written to normalize case
before matching. The propose -> confirm -> execute safety gate (`execute_confirmed_run`,
Fig. 6 in the main walkthrough) exists specifically around DGE, and nowhere
else in the app, for exactly this reason: it is the one place this
particular failure mode can occur at all.
