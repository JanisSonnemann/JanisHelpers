# QC excluded-sample markers Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Excluded samples are drawn as crosses on every one-point-per-sample QC plot (FACS, ELISA, histology).

**Architecture:** One shared shape scale helper; each point plot maps `shape = .excluded`. FACS reads `scored$excluded`; table-style plots take `excluded_keys` and use `qc_is_excluded_()`. Templates pass the keys and say so in the intro.

**Tech Stack:** R, ggplot2, testthat 3, rmarkdown templates.

**Spec:** `docs/superpowers/specs/2026-11-20-qc-excluded-markers-design.md`

## Global Constraints

- Statistics (boxes, whiskers, `flag`, labels) unchanged.
- Shapes: not excluded = 16, excluded = 4; legend title "excluded"; only the excluded key shown (`breaks = TRUE`), so no legend entry when nothing is excluded.
- `pkg::fn()` namespacing, `|>`, ASCII in R source.
- `devtools::check()`: 0 errors, 0 warnings.

## Review Focus

- No exclusions at all (`excluded_keys = NULL` / empty) -> plots render, no "excluded" legend.
- Batch plots: shape scale must not clash with the `colour = .batch` legend.
- ELISA CV plot excluded points above the limit stay red (colour) and become crosses (shape) -- both aesthetics independent.
- Scored data without an `excluded` column (older callers) -> all circles.
- End-point samples whose dates were resolved: keys and plot rows both use resolved dates (template already resolves both), so matching works.

---

### Task 1: Shape helper + FACS plots

**Files:** Modify `R/facs_qc.R` (`qc_plot_rows_`, `qc_plot_individuals_`, `qc_plot_batch_`, new helper); Test `tests/testthat/test-facs_qc.R`

**Produces:** `qc_excluded_shape_()` -> a `ggplot2::scale_shape_manual()` (values `c("FALSE" = 16, "TRUE" = 4)`, `breaks = "TRUE"`, `labels = "on exclusions sheet"`, `name = "excluded"`, `drop = FALSE`... limits `c("FALSE","TRUE")`); `qc_plot_rows_()` adds `.excluded` (factor/char of `excluded`, missing column -> FALSE).

- [ ] Test `qc plots draw excluded files as crosses`: score with `excluded_files = "<file of m-3>"`; built point layer of individuals and batch plots has shape 4 for that mouse's rows and 16 elsewhere; without exclusions all 16.
- [ ] Run, expect FAIL. Implement. Run `devtools::test(filter = "facs_qc")`, expect PASS. Commit.

### Task 2: Table plots + ELISA CV plot + templates

**Files:** Modify `R/qc_table.R` (`qc_table_plot_(..., excluded_keys = NULL)`), `R/elisa_qc.R` (`elisa_qc_cv_plot_(cv, max_cv = 20, excluded_keys = NULL)`), the three QC `skeleton.Rmd` intros + ELISA/histo calls; Test `tests/testthat/test-elisa_qc.R`

**Consumes:** `qc_excluded_shape_()`, `qc_is_excluded_()`.

- [ ] Test `qc_table_plot_ draws excluded samples as crosses` (synthetic_scores, keys `S6/spleen/end` -> S6 shape 4, others 16; `NULL` -> all 16; batch view too) and `elisa_qc_cv_plot_ draws excluded samples as crosses` (fixture, keys M1/spleen/end).
- [ ] Run, expect FAIL. Implement. Pass `excluded_keys` in templates; intro sentence: "...only flagged here (crosses on the point plots); they are dropped at the final export." Run `devtools::test()`, then `devtools::check()`. Commit.
