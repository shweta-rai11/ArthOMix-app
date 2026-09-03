# Discussion

## 5.1 Overview

ArthOMix was built to answer a practical problem in rheumatoid arthritis (RA)
bioinformatics: the transcriptomic, epigenomic and multi-omic evidence for a
candidate biomarker or pathway is routinely scattered across separate
scripts, separate cohorts, and separate analysts, with no single environment
in which a result can be generated, stratified by sex, validated against an
independent cohort or ancestry, and interrogated again without re-running
code by hand. The platform assembles these steps — transcriptomics,
methylomics, cross-omics integration, and multi-omics integration — into one
browser-based application built around a bundled reference cohort (an
anti-TNF-treated RA cohort merged from GSE93272 and GSE110169, with GSE15573
used for cross-ancestry validation and GSE89408 for cross-tissue
validation), while remaining fully usable on a user's own uploaded data or a
freshly fetched GEO series. The discussion below considers what this design
achieved, what it revealed about the underlying biology and the underlying
software, and where its current boundaries lie — organised around the four
themes the work centred on: single-omics analysis, multi-omics integration,
sex-stratified analysis, and the ArthOChat assistant.

## 5.2 Single-omics analysis

The transcriptomics and methylomics verticals were designed as parallel,
symmetrical pipelines — dataset intake, quality control, differential
analysis (DGE / DMP–DMR), co-expression or co-methylation network analysis
(WGCNA), candidate feature identification, feature selection, diagnostic
modelling, external validation, and biomarker reporting — rather than as
one privileged pipeline with a lesser methylation add-on. This symmetry
matters for a rheumatoid arthritis biomarker panel specifically, because
expression and methylation evidence are only interpretable together when
each side has been through an equivalently rigorous pipeline; an asymmetric
design would have made any later cross-omics concordance claim
uninterpretable (a "significant" expression hit could simply reflect a more
permissive expression pipeline rather than genuine biology).

A structural decision that shaped both verticals was to give every dataset
tab three genuinely isolated intake routes — upload, GEO accession, and a
bundled reference/preloaded dataset — gated by a `source_type` /
`is_bundled_reference` flag rather than a single shared code path with
conditional branches. This was not incidental plumbing: three separate
audit passes (2026-08-28 for transcriptomics, 2026-08-29 for methylomics,
and a full-review pass immediately after) found and fixed several modules
(Enrichment, Nomogram, WGCNA's local caches, GEO column auto-detection)
that silently fell back to the bundled RA gene panel, or failed to reset
stale results, when a user's own uploaded or GEO-fetched data was active.
The generalisable point for a discussion of this kind of platform is that
"three data sources sharing one UI" is a correctness hazard by default —
every module that reads a convenience fallback path has to be positively
proven never to reach it on non-reference data, and this proof degrades
silently under later refactoring unless it is pinned by a regression test.

The methylomics diagnostic/validation split produced one result worth
discussing on its own scientific merits, independent of the software that
produced it. A 12-CpG logistic regression classifier trained on the female
stratum (Control vs. RA) reached a training AUC of 0.754, an internal
held-out test AUC of 0.604, and an *external*-cohort AUC of 0.804. Taken at
face value the external AUC exceeding the internal-test AUC looks like an
error; the platform's own shift-diagnostics view (a joint PCA of training
and external samples) showed instead that the external cohort is genuinely
displaced from the training distribution in feature space, which explains
why sensitivity and specificity at the frozen 0.5 threshold were poor
despite a respectable AUC — the model still ranks cases correctly but its
calibrated decision boundary does not transfer. This is a useful negative
result for the wider thesis: it demonstrates that external validation
without a cohort-shift diagnostic can be actively misleading (a single AUC
number would have been reported as "generalises well"), and it argues that
any biomarker panel produced by this pipeline should be reported with its
shift diagnostic attached, not as a bare AUC.

Several genuine defects were also caught by treating the single-omics
pipelines as targets for systematic, code-independent testing rather than
only manual spot-checking. Two are worth naming because they are the kind
of bug that produces confidently wrong output rather than a crash:
`eda_parse_upload()` throws a raw, undocumented error (rather than its own
`list(ok = FALSE, ...)` contract) on a single-feature-row upload with
multiple sample columns, and a GEO metadata column auto-guesser
(`guess_col()`) checks candidate column names in declaration order rather
than in priority order, so an administrative `status` or generic
`characteristics_ch1` column can be selected ahead of the genuinely
informative `disease state:ch1` column on a real series (confirmed on
GSE93272). Neither failure mode announces itself; both would quietly change
which samples are labelled "case" and "control" for every downstream
analysis. Their discovery, rather than their existence, is the point worth
making in a discussion of methodology: a platform that recomputes live from
arbitrary user-supplied metadata needs its column-inference logic under the
same test discipline as its statistics, because a labelling bug upstream is
invisible in every downstream number.

## 5.3 Multi-omics integration

The multi-omics vertical implements three conceptually different
integration strategies side by side — DIABLO (`mixOmics::block.splsda`,
supervised, data-adaptive keepX selection via `tune.block.splsda`),
Similarity Network Fusion (`SNFtool::SNF`, unsupervised, with eigengap-based
cluster-number estimation rather than a fixed K), and MOFA2 (unsupervised
factor analysis) — plus a head-to-head Compare tab that benchmarks DIABLO
and SNF against single-omics baselines through the same leakage-safe nested
cross-validation used elsewhere in the app. Offering three integration
philosophies rather than one forces an explicit answer to a question that
multi-omics papers often leave implicit: whether the fused representation
should be optimised to predict a known outcome (DIABLO), to discover
unsupervised patient structure without appealing to outcome labels at all
(SNF), or to decompose shared and modality-specific variance without a
supervised target (MOFA2). Making all three available on the same dataset,
with the same held-out evaluation, is what allows the Compare tab's claim —
that a fused model outperforms (or does not outperform) any single omic
layer — to mean something, rather than being a comparison across
differently-tuned, differently-evaluated pipelines.

A deliberate design choice in this vertical is that even its "precomputed"
route is not a static table browser: selecting the bundled reference cohort
re-runs the identical live DIABLO/SNF pipeline against that cohort's own
already-feature-selected matrices (loaded from a saved `.rds` fit), rather
than displaying cached output. This matters for reproducibility framing —
a user inspecting the "preloaded" result is looking at the same code path
that runs on their own uploaded data, with a provenance note distinguishing
the two, rather than a frozen figure that could silently drift out of sync
with the live engine as the codebase evolves.

The cross-omics evidence classifier, `cx_classify_evidence()`, is the
clearest illustration of a reuse principle that shaped the whole multi-omics
vertical: rather than building a second, competing "biomarker tier" scheme
for the Multi-Omics Biomarker Card, the card was made strictly read-only
over the same evidence table (`multi_results$concordance$df`) and the same
six-tier classification (Strong candidate / Moderate candidate /
Expression-only / Methylation-only / Discordant / Insufficient evidence)
that Cross-Omics Integration itself produces. The alternative — an
independently defined "multi-omics-specific" scoring rule — would have
created two answers to "how strong is the evidence for gene X," diverging
silently as either implementation changed. Reuse here is not merely an
engineering economy; it is what keeps "why is this candidate tier X"
traceable to one place in the codebase, which is a precondition for the
result being auditable at all in a thesis or publication context.

The clearest architectural limitation the multi-omics work surfaced is that
the Methylomics vertical and the Multiomics vertical share no reactive
state whatsoever: a DMP/DMR run performed live in the Methylomics tab does
not appear as a methylation block in the Multi-Omics Dataset Workspace.
Methylation data currently enters the multi-omics engine only via manual
re-export/re-upload of a beta- or M-value matrix, or via the bundled
reference cohort's own precomputed fit. This is a real gap rather than a
subtle bug, and it means that, as implemented, a fully live "run DMPs, then
immediately integrate them with live DGE results" workflow is not possible
within a single session — a limitation discussed further in §5.6.

## 5.4 Sex-stratified and sex-differential analysis

Sex is treated throughout ArthOMix as a first-class analytical axis rather
than as an optional covariate bolted onto a single pooled model, and the
platform is careful — more careful than most exploratory pipelines — to
keep three distinct claims from being conflated: a *sex-stratified* result
(the same analysis run separately within each sex, which licenses no claim
about a sex difference on its own), a *sex-differential* result (a formal
statistical test that the effect itself differs by sex), and "sex-specific"
(deliberately avoided as a term in the platform's own documentation,
precisely because it invites the reader to infer a differential finding
from a merely stratified one). This distinction is not just terminological
housekeeping: two lists thresholded separately within strata of unequal
size will differ from one another even when the true underlying effect is
identical in both sexes, purely as an artefact of the smaller stratum's
lower power. A platform that reports separate female and male gene lists
without also offering a genuine interaction test invites exactly this
misreading.

ArthOMix addresses this by implementing both levels explicitly and keeping
them architecturally distinct. Differential expression, feature selection,
diagnostic modelling, the nomogram, and cross-tissue/cross-ancestral
validation all support **sex-stratified** operation — in feature selection
this is not even optional (LASSO, random forest and SVM-RFE are "always
modelled separately, never pooled" per sex, with a pooled candidate list
built afterward as the union of the two strata's panels rather than as its
own independently-selected set). Separately, a dedicated **Sex Interaction
Analysis** submodule exists in both transcriptomics and methylomics,
fitting a single joint model with an explicit diagnosis-by-sex interaction
term (`~ diagnosis * sex`, limma/eBayes) and testing the interaction
coefficient itself — the only construct in the platform that can support an
actual "this gene's disease effect differs by sex" claim, gated on a
minimum of two samples in every diagnosis-by-sex cell and at least twelve
samples overall.

The same discipline carries into multi-omics: a dedicated Sex-Stratified
engine (`multiomics_sexstratified_engine.R`) runs the same nested,
leakage-safe DIABLO/random-forest cross-validation independently for
pooled, female-only, male-only, or both-separately modes (minimum six
samples per stratum), reporting AUROCs with confidence intervals per
stratum rather than a single pooled figure, and the bundled reference
cohort itself is fit as **six separate DIABLO models** — female and male,
crossed with Adalimumab-responder, Etanercept-responder, and a
drug-pooled-responder cell — each fit independently on its own matched
samples and its own selected feature subset, rather than one model sliced
six ways after the fact. Cross-omics integration follows the weaker but
still useful stratified pattern: gene–CpG concordance is computed and
cached independently for the female, male, and all-sample subsets, letting
a user compare the three concordance tables side by side without any
formal claim of divergence being made by the platform itself.

The rationale documented alongside this design is explicitly statistical
rather than biological: the bundled training cohort is markedly
sex-imbalanced (183 whole-blood samples in total — 145 female, comprising
86 RA and 59 control, against only 38 male, comprising 17 RA and 21
control), and the accompanying methodology notes repeatedly flag the small
male stratum as the reason several downstream steps take extra
methodological care — for example, using the (larger, more stable) female
network as the reference for module-preservation testing rather than
attempting a symmetric two-way comparison, and using a permutation-based
null to avoid mistaking the male stratum's smaller size for weaker
preservation. It is also notable, and worth stating plainly in this
discussion, that no disease-biology rationale for stratifying by sex (for
example, RA's well-documented female predominance, or literature on
sex-differential anti-TNF response) is present anywhere in the codebase's
own comments or methodology notes; the platform's justification for
sex-stratification is entirely about protecting against confounding and
under-powered artefacts in an unbalanced cohort, not about a specific
biological hypothesis being tested. A published or defended version of
this work should supply that epidemiological framing from the literature
directly, since the tool itself only argues the statistical case for
*how* to stratify, not the biological case for *why* a sex difference
might be expected in this disease and this drug class.

## 5.5 ArthOChat and the reliability of an embedded LLM assistant

ArthOChat was included on the premise that a working analyst benefits from
being able to ask "what did this panel show" in natural language without
leaving the analysis, but embedding a generative model inside a scientific
tool creates an obvious risk that is easy to gesture at and hard to
quantify: the assistant fabricating a plausible-sounding number that was
never actually computed. Rather than treating this as an unmeasurable risk
to be mitigated by prompt wording alone, the platform builds in two
independent defences and then measures them directly. The first is a
system prompt that requires the model to distinguish this session's
*live, computed* results from *published-methodology* background knowledge
and to say plainly when a given sub-module has not yet been run. The second
is a code-level, post-hoc scanner
(`arthochat_detect_ungrounded_reference()`) that parses the model's own
draft response for a mention of a specific sub-module the current context
marks as not-yet-run, appending an explicit caveat if one is found — a
safety net that does not depend on the underlying model reliably following
its instructions under pressure, since it was added specifically because a
small local model does not.

The value of this design is best seen in what the platform's own
verification harness found, since it separates "the mechanism worked" from
"the mechanism, as implemented, has a specific weak spot." Across a 205-
question protocol spanning five question types and all 41 registered
sub-modules, the assistant answered 160/164 gradeable questions correctly
(97.6%), with transcriptomics, methylomics and multi-omics all at 100% and
cross-omics at 8/12 (66.7%). Every wrong answer traced to the same
proximate cause and it is a genuinely informative one for anyone building
an LLM assistant over structured application state: two cross-omics
sub-modules' context formatters simply never include a "currently loaded
dataset" line, and when asked what dataset was loaded, the model did not
say "that's not shown to me" — in both cases it confidently invented a
specific, wrong answer. This is a cleaner and more useful finding than a
generic "LLMs can hallucinate" caveat, because it localises the failure
mode precisely: the model degrades to confabulation not when it is asked
something outside its scope (the isolation test showed 41/41 sub-modules
never leaked a neighbouring sub-module's data), but when it is asked
something that is *silently absent* from its own context rather than
*explicitly marked absent*. The practical implication carried forward into
the platform's own design is that every context builder must positively
state what is unavailable, rather than relying on the model to infer
absence from omission — omission is exactly the condition under which this
model chose to fabricate.

A second, smaller but methodologically relevant decision was to run the
assistant at temperature zero with a fixed seed, specifically because the
verification harness had itself observed the same unanswerable question
producing two different fabricated numbers on repeated queries under the
model's non-zero default sampling temperature. This does not address
fabrication as a phenomenon — a deterministic model can still confidently
state a wrong answer, and did, in the missing-dataset-line cases above —
but it does mean that a given fabricated or correct answer no longer drifts
from one run to the next for the same session state, which is a
precondition for the kind of exact-match, automated grading the
verification harness itself relies on. Differential Expression is
additionally the only tool exposed to the assistant that is gated behind an
explicit propose-then-confirm step before it executes anything, because it
is the one tool in the platform whose inputs are matched as exact,
case-sensitive strings rather than case-insensitively — a live test found
that `ref_group = "HC"` succeeds while `ref_group = "Hc"` fails silently
with an unrelated "fewer than 6 samples match" message. Treating that one
tool differently from the seven others is a direct, evidence-based response
to a concrete failure mode observed in this application, not a generic
caution applied uniformly.

## 5.6 Software engineering discipline as a contributor to scientific validity

A point that belongs in this discussion, even though it is not itself a
biological result, is that the platform's testing programme (94 files,
941 `test_that()` blocks at the time of writing, spanning unit,
`testServer`, UI, end-to-end and scientific-contract layers) functioned as
a genuine discovery process rather than a confirmation exercise, and it is
worth being explicit about what it found because several of the defects
would otherwise have silently produced wrong scientific conclusions rather
than visible crashes. Beyond the GEO column-priority and single-row-upload
bugs already discussed in §5.2, the multi-omics cohort-harmonization
duplicate-sample-ID guard was found to be dead code (`intersect()` already
returns unique values, so the `duplicated()` check downstream of it can
never fire), meaning a matrix containing the same sample ID five times
would have five-fold-duplicated evidence silently collapsed to its first
row rather than being flagged; a related identifier-harmonization table
never marked a blank sample ID "Invalid" because of an R-specific quirk
(`list[[""]]` always returns `NULL`, even for a genuinely-empty-named
element, so the assignment meant to flag it silently no-ops); and a fully
implemented Results Summary and reproducibility-report submodule was found
to be unreachable in the running application because it had never been
added to the multi-omics module registry. None of these are the kind of
error a user would notice from a plot looking "wrong" — they are exactly
the class of defect that erodes trust in a platform's numbers without ever
producing a visible symptom, and finding them by systematic testing rather
than by a downstream user noticing an implausible result is itself part of
the argument that the platform's outputs can be relied upon.

## 5.7 Limitations and future work

Several boundaries of the current implementation are worth stating plainly
rather than leaving implicit. First, the Methylomics and Multiomics
verticals remain reactively disconnected (§5.3); wiring a live DMP/DMR
result directly into the Multi-Omics Dataset Workspace, rather than
requiring a manual re-export, is the clearest concrete piece of future
work the architecture points to. Second, end-to-end (`shinytest2`) coverage
is written but auth-gated behind a Supabase login screen added concurrently
with this testing work, and skips cleanly rather than running in any
environment without a configured test account — a disclosed limitation of
the test suite's actual current coverage, not a claim that E2E behaviour
has been verified everywhere it is written. Third, the sex-stratification
design, while statistically well-motivated, currently rests on a training
cohort with a genuinely small male stratum (n = 38); the platform's own
minimum-sample gates (for example, six samples per stratum in the
multi-omics sex-stratified engine) are the floor for the pipeline to run at
all, not evidence that estimates at that floor are stable, and any
biomarker claim drawn from the male stratum specifically should be read
with that in mind. Finally, as noted in §5.4, the case for *why* sex
differences might be expected in RA and in anti-TNF response is not
currently documented anywhere in the platform itself and should be
supplied from the epidemiological and pharmacogenomic literature when this
work is written up, rather than left to be inferred from the presence of
the stratification machinery alone.

## 5.8 Summary

Taken together, the single-omics, multi-omics, sex-stratified, and
ArthOChat components of this work support one overall claim: that a
scientifically defensible multi-omics platform is as much a discipline of
architectural isolation, statistical honesty about what a given analysis
mode can and cannot claim, and systematic adversarial testing of its own
outputs, as it is a matter of which integration algorithm is chosen. The
platform's own negative and boundary-case findings — a shift-diagnosed
external validation, a documented distinction between stratified and
differential sex effects, a precisely localised LLM fabrication mode, and
several silent-failure bugs caught before they reached a user — are, in
that sense, as much a part of its contribution as the pipelines that
produced clean results.
