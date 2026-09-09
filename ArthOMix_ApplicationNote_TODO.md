# ArthOMix Application Note — master to-do list

Companion to `ArthOMix_BioinfAdv_ApplicationNote_template.md` (the manuscript skeleton). This file is the execution checklist — work top to bottom within each track, tracks run in parallel. Check items off as you go (`- [x]`).

**2026-09-08 update:** scope expanded from 2 case studies to 4 (one per pipeline — Transcriptomics, Methylomics, Cross-Omics, Multi-Omics), replacing Track B/B2 below with Track B1-B4. This is the right call for the paper's argument (it directly demonstrates the four-pipeline claim in Introduction §2) but it is two more live runs than the original plan — **add ~1-2 days to whatever timeline you were working toward.** Methods §2.1-2.6 are now fully drafted (done, see Track C); the four live runs (Track B1-B4) are the remaining critical path.

---

## TRACK A — Code verification (do first, ~half a day, gates Methods §2.4/2.5 and Availability)

- [ ] Re-run the grep from this session to confirm GSE89253 is still unused: `grep -rn "GSE89253" R/` → expect no hits. If it now has hits, someone wired it in — verify how before citing it.
- [ ] Confirm scorecard fix is still intact: open `R/multiomics/functions/multiomics_helpers.R`, find `multi_qc_scorecard()`, confirm `n_below > 0` still returns `"fail"` (not overridden since 2026-09-08).
- [ ] Confirm cross-omics MR/backfill fix (commit `cba2453`): read current `mod_cross_integration.R` and `mod_cross_mr_stage.R` in full, write one sentence for your own notes on exactly what significance/FDR gate is applied to instruments now — you'll need this sentence verbatim for Methods §2.4.
- [ ] Confirm renv.lock is clean: `cd ArthOMix && Rscript -e 'renv::restore()'` on a clean checkout (or at minimum `renv::status()`) — must complete with no version-resolution errors.
- [ ] Confirm Docker build succeeds from a clean clone (do this in a scratch directory, not your working copy):
  ```
  git clone <repo> /tmp/arthomix-clean-check && cd /tmp/arthomix-clean-check/ArthOMix
  docker build -t arthomix-check .
  docker run --rm -p 3839:3838 arthomix-check
  ```
  Open `http://localhost:3839` and confirm the app actually loads. **Do not write "Docker: builds and runs" in Availability until this succeeds.**
- [ ] Confirm LICENSE file content matches what you intend to state in the manuscript (file exists at repo root, already confirmed present — just read it and note the license name/version for §4/§9 of the template).
- [ ] Check `caret::createDataPartition()` silent-failure bug (methylomics feature selection, `mod_methyl_featureselection.R:~788-796` per the 2026-09-07 audit) — if you plan to describe methylomics feature selection with a train/holdout split in Methods, confirm this doesn't silently produce a 100/0 split on your chosen dataset.
- [ ] Decide whether a live public demo instance exists or will exist by submission. If not, Availability states "run locally" only — do not imply a hosted demo.
- [ ] Grep for `ANTHROPIC_API_KEY` / auth gating one more time if ArthoChat (§2.6) will be described as a headline feature — confirm the session-independent spend cap (not just per-browser-session) is deployed, since the earlier version's own error message told users how to bypass it by reloading.

**Go/no-go:** if any item above surfaces a *new* problem (not already known), stop and fix or scope it out of the paper before drafting Methods for that sub-module. Don't draft prose describing behavior you haven't just verified.

---

## TRACK B1 — Case study 1: Transcriptomics (start Day 1, lower risk)

- [ ] Pick a public GEO dataset. Candidates already verified in memory: `reference_geo_transcriptomics_test_datasets.md` — **re-verify the accession still fetches** before committing (GEO series occasionally get withdrawn/updated).
- [ ] Run live in the app: intake → preprocessing → DE (limma/voom) → WGCNA → feature selection → diagnostic model, whichever subset you'll describe in Results.
- [ ] Record: n samples, n genes tested, FDR/logFC thresholds used, top hits, one representative plot (volcano).
- [ ] Export the plot at publication resolution; save the provenance manifest.
- [ ] Write 2-3 sentences for template §7 Case study 1.

---

## TRACK B2 — Case study 2: Methylomics (start Day 1-2, parallel with B1)

- [ ] Pick a public GEO methylation dataset. Candidates in `reference_geo_methylomics_test_datasets.md` — re-verify accession and array type (450K/EPIC/EPICv2) before committing.
- [ ] Run live: intake → QC (sex check, detection p-value) → normalization → DMP → DMR (if time allows) → feature selection.
- [ ] Record: n samples, array type, normalization method used, FDR threshold, top DMPs/DMRs.
- [ ] Export one representative plot; save the provenance manifest.
- [ ] Write 2-3 sentences for template §7 Case study 2.

---

## TRACK B3 — Case study 3: Cross-Omics (start Day 2-3, needs paired data)

- [ ] Pick a paired dataset (matched expression + methylation, same patients) from `reference_geo_multiomics_test_datasets.md` — re-verify against current GEO state.
- [ ] Run live: Expression-Methylation Integration → Biomarker Convergence → (MR, if time allows).
- [ ] Record: n genes classified per significance category, n genes significant by convergence category, the Bonferroni-within-gene/BH-FDR result described in Methods §2.4 — confirm it actually runs this way on your data.
- [ ] Export one representative plot (quadrant/scatter plot); save the provenance manifest.
- [ ] Write 2-3 sentences for template §7 Case study 3.

---

## TRACK B4 — Case study 4: Multi-Omics (start Day 3, the schedule risk — heaviest compute, most audit history)

- [ ] Use the same paired dataset as B3, or a second verified pair.
- [ ] Run live: Cohort Harmonization → DIABLO/SNF Integration → Biomarker Discovery with sealed hold-out (confirm the 20% default hold-out actually runs, per Methods §2.5's caution note) → Pathways.
- [ ] Record the actual pooled OOF AUROC + DeLong CI from Biomarker Discovery (not `perf()`'s number). **Decision gate:**
  - [ ] If CI excludes 0.5 (genuinely above chance) → use as headline result.
  - [ ] If null/negative → write it up honestly as a negative biological finding on a working, reproducible pipeline. Do not re-run with different parameters searching for a better number — that's p-hacking and it will blow your timeline.
- [ ] Record SNF cluster association test (if used) with its p-value, whatever it is.
- [ ] Export one representative plot (DIABLO sample plot, SNF cluster heatmap, or ROC curve).
- [ ] Write 2-3 sentences for template §7 Case study 4.

---

## TRACK C — Writing (parallel with B1-B4, gated only by Track A per subsection)

- [x] Methods §2.1 Data formats — drafted 2026-09-08 from code.
- [x] Methods §2.2 Transcriptomics — drafted 2026-09-08 from code.
- [x] Methods §2.3 Methylomics — drafted 2026-09-08 from code.
- [x] Methods §2.4 Cross-Omics — drafted 2026-09-08 from code, upgraded with the verified min-p Bonferroni fix.
- [x] Methods §2.5 Multi-Omics — drafted 2026-09-08 from code; still needs a live-run sanity check (Track B4) before final.
- [x] Methods §2.6 ArthoChat — drafted 2026-09-08 from code.
- [ ] Methods §2.5 Multi-Omics — write **last**, after Track B2 numbers exist.
- [ ] Methods §2.6 ArthoChat — state turn/spend limits, confirm session-independent budget wording matches Track A item 9.
- [ ] Introduction paragraph 1 (field context, integration strategies, sex-in-disease motivation) — no dependency on results, can be drafted Day 1.
- [ ] Introduction paragraph 2 (competitive gap vs. named comparator tools, closing pitch sentence) — needs the feature-comparison table drafted alongside it (see Track D).
- [ ] Abstract — write last, after Results numbers exist. Draft to the ~250-word Motivation/Results/Availability structure; have a trimmed 200-word version ready in case the stricter cap is confirmed.
- [ ] Results — both case studies, ~150-250 words total, from Track B/B2 output.
- [ ] Conclusion — ~90-100 words, no new claims.
- [ ] Data Availability statement — GEO accessions used, repo URL, confirm against Track A.
- [ ] Funding statement.
- [ ] Conflict of Interest statement.
- [ ] Acknowledgements.
- [ ] Reference list — budget against the 15-reference cap per template §13's allocation; compile author-year formatted list.

---

## TRACK D — Figures, tables, supplementary (start once Methods module list is locked)

- [ ] Figure 1: inputs → outputs schematic, one numbered box per Methods subsection (2.2-2.6), in the same order as the app's tabs (Transcriptomics → Methylomics → Cross-Omics → Multi-Omics), ArthoChat drawn cross-cutting.
- [ ] Feature-comparison table (ArthOMix vs. named comparator tools) — every checkbox must be true of current code, verify each row against actual behavior, not intent.
- [ ] Supplementary: extended Methods for anything trimmed from main text.
- [ ] Supplementary: all case-study figures not selected as the one main figure.
- [ ] Supplementary: deployment/architecture diagram if relevant.
- [ ] Confirm final figure/table count in main text ≤ 3 combined (per journal limit).

---

## TRACK E — Admin (any time, doesn't block writing)

- [ ] Email Bioinformatics Advances editorial office to resolve the 200-word vs. Motivation/Results/Availability abstract-structure discrepancy; also confirm current page/reference limits are still 4 pages / 15 refs / 3 figs-tables.
- [ ] Get the actual submission template (LaTeX or Word) from the OUP submission portal — do not hand-format from scratch.
- [ ] Confirm author list, affiliations, ORCID iDs, corresponding author.
- [ ] Decide and confirm license terms match what's in `LICENSE` (Track A item 6).
- [ ] Draft cover letter (1 paragraph: what the tool does, why it fits the journal's scope, confirm no prior/concurrent submission elsewhere).
- [ ] If a supervisor/coauthor needs to review — get them a draft by end of Day 8 at the latest; this is the step most likely to consume your buffer.

---

## FINAL PASS — before hitting submit

- [ ] Every `[AUTHOR]` and `[VERIFY]` bracket removed from the manuscript.
- [ ] Every quantitative claim in Results traced to a run you personally did this week (not a bundled table).
- [ ] No reference to GSE89253 as a validation cohort.
- [ ] "Uniquely supports" / comparator claims tied 1:1 to feature-comparison table checkboxes.
- [ ] Docker build re-confirmed on a truly clean clone (not the machine you've been developing on).
- [ ] Word count / figure count / reference count within confirmed journal limits.
- [ ] Repo link, license, and data-availability links all resolve publicly (test in a private/incognito browser window, logged out).
- [ ] Abstract structure matches whatever Track E's editorial-office email confirmed.
- [ ] Final proofread pass for any leftover "strong candidate" / "RA-associated" language sourced from a Biomarker Card screenshot without a checked significance gate (audit finding #5).

---

*Built 2026-09-08. Cross-reference `ArthOMix_BioinfAdv_ApplicationNote_template.md` §0 for the underlying audit findings behind each Track A/Final Pass item.*
