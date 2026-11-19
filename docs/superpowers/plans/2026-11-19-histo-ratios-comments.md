# Histology ratios, comments, slimmer QC Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Several ratio metrics in one `histo_calc_ratio()` call, a free-text `comment` column in curated histology sheets shown in the QC report, and no batch / time outlier plots in the histology QC report.

**Architecture:** Vectorise `histo_calc_ratio()` over named numerators; `histo_read_curated()` peels off a `comment`/`comments` column and attaches it as a `comment` column of the long table (only when present); the QC skeleton gains a Comments section and loses two sections.

**Tech Stack:** R package (dplyr, tidyr, readxl, testthat 3, writexl fixtures).

**Spec:** `docs/superpowers/specs/2026-11-19-histo-ratios-comments-design.md`

## Global Constraints

- CLAUDE.md code style: `|>`, `pkg::fn()` namespacing, `\uXXXX` for non-ASCII in R strings, ASCII in roxygen.
- `devtools::document()` after roxygen changes; `devtools::check()` 0 errors / 0 warnings.
- Comments are not exported to the database; histology only (ELISA untouched).

## Review Focus

- Sheet without a comment column: long table keeps exactly the 5 existing columns (old tests, old preps).
- Ratio rows when `comment` exists: the new rows carry the sample's comment (no `all_of()` failure).
- Header case `Comment` / `COMMENTS`: treated as comment, not a metric.
- Numeric-looking comment (e.g. `5`): stays text, no metric created.
- Old prep `.rds` without `comment`: QC report renders with "No comments".

---

### Task 1: Vectorised `histo_calc_ratio()`

**Files:** Modify `R/histo_calc.R`; Test `tests/testthat/test-histo_calc.R`

**Interfaces:** Produces `histo_calc_ratio(data, numerator, denominator, name = NULL, scale = 100)`. `numerator` is a character scalar (with `name`) or a named character vector (names = new metrics, `name` ignored / must be NULL). `denominator` length 1 or `length(numerator)`. Extra non-key columns of `data` (e.g. `comment`, group cols) are copied from the numerator rows.

- [ ] Test `histo_calc_ratio computes several named ratios`: on the fixture, `numerator = c(pct_cresc = "crescents", pct_nec = "necrosis")`, `denominator = "gloms"` -> M1 `pct_cresc == 22/85*100`, both metrics present, row count = original + 2 x mice.
- [ ] Test `histo_calc_ratio errors on bad vector input`: unnamed vector of length 2 -> error "named"; duplicated names -> error "unique"; `denominator` of length 3 for 2 numerators -> error "denominator"; existing metric name -> "already exists".
- [ ] Test `histo_calc_ratio keeps extra columns`: add `comment` column to fixture data (M1 = "torn") -> ratio row of M1 has `comment == "torn"`.
- [ ] Run tests, see them fail; implement; run `devtools::test(filter = "histo_calc")`, all pass (old tests too).
- [ ] Commit `feat(histo): several ratios in one histo_calc_ratio() call`.

### Task 2: Comment column in `histo_read_curated()`

**Files:** Modify `R/histo_read.R`, `data-raw/make_curated_templates.R` (README line: optional `comment` column, free text, shown in the QC report, not exported), regenerate `inst/templates/histo_curated_template.xlsx`; Test `tests/testthat/test-histo_read.R`

**Interfaces:** Long table gets column `comment` (character, `NA` for empty, after `value`) only when the sheet has a column matching `^comments?$` ignoring case. More than one such column -> error.

- [ ] Test `histo_read_curated reads a comment column as text`: writexl sheet `mouse_ID = c("M1","M2"), gloms = c(80, 70), Comment = c("torn section", NA)` -> metrics only `gloms`; `comment` of M1 `"torn section"`, M2 `NA`.
- [ ] Test numeric comment `comments = 5` -> no metric `comments`, `comment == "5"`.
- [ ] Existing `expect_named` test stays (no comment column -> 5 columns).
- [ ] Implement, run `devtools::test(filter = "histo_read")`; rebuild template via `Rscript data-raw/make_curated_templates.R`; commit `feat(histo): optional comment column in curated sheets`.

### Task 3: QC report -- Comments section, drop batch / time plots

**Files:** Modify `inst/rmarkdown/templates/histo-qc-template/skeleton/skeleton.Rmd`, `inst/scripts/histo_prep_template.R` (vector-form `prepare`, comment on batch plots removed), roxygen of `histo_prep_save()` (`data` may carry `comment`); Test `tests/testthat/test-histo_prep.R`

**Interfaces:** Internal `histo_qc_comments_(data, excluded_keys)` in `R/histo_qc.R` -> tibble `mouse_ID`, `tissue`, `timepoint`, `comment`, `excluded` (distinct per sample, non-empty comments only; zero rows if no `comment` column).

- [ ] Test `histo_qc_comments_` on small tibble: two metric rows of M1 with same comment -> one row; data without `comment` -> 0 rows.
- [ ] Skeleton: delete `### By batch` and `### Over time` chunks; add `## Comments` after "Empty values" using the helper ("No comments." when empty).
- [ ] Extend the existing render test (or add one) so a prep with a comment renders and the HTML contains "Comments" and not "By batch".
- [ ] Run `devtools::test(filter = "histo")`, `devtools::document()`, `devtools::check()`; update CLAUDE.md data-structures line for histo; commit `feat(histo-qc): comments section, drop batch and time plots`.
