# ArthOChat sub-module grounding tests

These are manual test scripts that check whether ArthOChat (the app's AI
assistant) actually uses the real data it's given, instead of making things
up. They are **not** run in CI — they need a local Ollama server running the
model set in `ARTHOMIX_OLLAMA_MODEL`, and results can vary slightly each run.

## How it works

Each script:

1. Creates a fake result containing a made-up number (e.g. `1003`) and puts
   it into the app's normal results data, exactly like a real analysis would.
2. Builds the same system prompt the app itself builds for that sub-module
   (`build_scoped_assistant_context()`, the same code `mod_arthochat.R` uses).
3. Asks the live model a question and checks whether its answer contains the
   right number — and does **not** contain a neighboring sub-module's number
   (proof the model only sees the one sub-module it should).

Because the "correct answer" is a number we invented ourselves, grading is
just an exact text match — no guesswork or reading the model's prose by eye.

## Files here

- `test_submodule_grounding.R` — asks 1 question per sub-module (41 total).
  Run with:
  ```
  Rscript tests/arthochat_verification/test_submodule_grounding.R
  ```
- `results.csv` — its output, one row per sub-module.
- `test_205_protocol.R` — a bigger version: 5 different kinds of questions
  per sub-module (205 total). See below for what each question type checks.
- `results_205_protocol.csv` — its output, one row per question.

## Results (one run, 2026-09-02)

**Quick test — 1 question x 41 sub-modules:**

| | |
|---|---|
| Correct | 40 / 41 |
| Stayed isolated (never leaked another sub-module's data) | 41 / 41 |
| Answer time | 2.2s–16.9s (one slow first call, rest were fast) |

**Full test — 5 question types x 41 sub-modules = 205 questions:**

| Question type | What it checks |
|---|---|
| Basic / conceptual | Can it explain what the analysis does, in plain terms? (not auto-graded — logged for a human to read) |
| Data scope | Does it correctly say what's loaded and whether the analysis has run? |
| Numeric, basic | Can it report the exact planted number back? |
| Numeric, advanced | Can it add two planted numbers together correctly? |
| Boundary / adversarial | Asked about a *different* sub-module — does it resist making up an answer? |

| | |
|---|---|
| Gradeable answers correct | 160 / 164 (97.6%) |
| Transcriptomics / Methylomics / Multi-Omics | 100% correct |
| Cross-Omics | 8 / 12 (66.7%) — see below |

## What went wrong (Cross-Omics only)

Every wrong answer was in Cross-Omics, and both problems come down to how
that vertical formats its data for the model — not the model misreading real
numbers:

1. **One sub-module (`integration`) uses different field names.** Unlike
   every other sub-module, it has no field literally called `check_value`,
   so asking for that exact name isn't a fair question for it. The model
   sometimes guessed a wrong number, sometimes correctly admitted it didn't
   know.
2. **Two sub-modules (`biomarkerconv`, `mrstage`) never mention which
   dataset is loaded.** Checked directly in `modules_index.R`: their shared
   formatter simply has no "currently loaded dataset" line. When asked what
   dataset was loaded, the model didn't say "that's not shown to me" — both
   times it confidently made up the same specific, wrong answer. This is the
   more important finding: **this model sometimes invents a plausible answer
   instead of saying "I don't know" when the information is genuinely
   missing**, as opposed to just not scoped correctly.

Everything else — the actual mechanism that scopes each sub-module's data —
worked correctly in every test, including every other question type on
these same two sub-modules.

## Why only Differential Expression (DGE) needs this kind of testing

Every other tool the assistant can run (gene search, GWAS search, PubMed
search) ignores letter case in its inputs, so a typo like `hc` vs `HC`
doesn't matter. DGE is the one exception — it compares group names as
**exact, case-sensitive strings**. Tested live, same function, same data:

```
ref_group = "HC", comp_group = "RA"  ->  works, real result returned
ref_group = "Hc", comp_group = "Ra"  ->  fails: "Fewer than 6 samples match"
```

Only the letter casing changed, and it broke. That's why DGE — and only
DGE — goes through an extra "propose → confirm → run" safety step before
executing anything for real.
