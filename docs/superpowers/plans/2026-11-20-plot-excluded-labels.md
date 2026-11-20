# Plot helpers: excluded mice, labels, bubble legend, paired layout -- Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Opt-in excluded-mouse highlighting and mouse-ID labels on the per-mouse plots, a self-explaining bubble legend, and a group-clustered `plot_paired()`.

**Architecture:** Shared internals in `R/plot_box.R` (`plot_split_excluded_()`, `plot_label_layer_()`); each plot keeps its boxes/medians/stats on kept rows and draws one point layer over all rows with shape/colour by `.excluded`.

**Tech Stack:** R, ggplot2, ggrepel (Suggests, optional), testthat 3.

**Spec:** `docs/superpowers/specs/2026-11-20-plot-excluded-labels-design.md`

## Global Constraints

- Base pipe `|>`; every external call `pkg::fn()`; non-ASCII in R strings as `\uXXXX`; roxygen ASCII only.
- New arguments default to current behaviour except: rows with `excluded == TRUE` are dropped (message) by default.
- `ggrepel` only behind `qc_has_ggrepel_()`.
- `devtools::check()`: 0 errors, 0 warnings, no new NOTE.
- Excluded colour `"red3"`, cross = shape 4 (as the QC reports).

## Review Focus

- Data without an `excluded` column and without an `id`/`mouse_id` column, `label = FALSE`: must work exactly as before (no column lookup). Test in Task 1.
- `excluded` with `NA` values: treated as not excluded. Test in Task 1.
- `show_excluded = TRUE` with no excluded rows: no red layer, no legend entry, no error. Test in Task 1.
- `plot_population(scale = "log")` with excluded rows: pseudocount applied to them too. Covered by Task 1 test via `plot_population`.
- Bubble plot where all log2FC are 0 / NA: legend limits must not be `c(0, 0)`. Test in Task 3.

---

### Task 1: excluded + labels in `plot_box()` / `plot_population()`

**Files:** Modify `R/plot_box.R`; Test `tests/testthat/test-plot_box.R`

**Interfaces — Produces:**
- `plot_split_excluded_(data, show_excluded) -> data` with logical `.excluded` column (never NA); drops excluded rows with message when `!show_excluded`; message when `show_excluded` and no `excluded` column.
- `plot_label_layer_(mapping, position, colour = NULL, data = NULL) -> layer` (ggrepel or `geom_text`).
- `plot_box(..., show_excluded = FALSE, label = FALSE, id = mouse_id)`; `plot_population()` passes all three through.

- [ ] Tests: `plot_box drops excluded rows with a message by default` (`expect_message(..., "excluded")`, `nrow(p$data)` = kept); `plot_box draws excluded rows as crosses with red labels` (point layer data has both, layer count +1 text layer, box layer data excludes them, shape scale present); `plot_box labels every point on request`; `plot_box works without id column when nothing is labelled`; `NA in excluded counts as kept`.
- [ ] Implement; run `devtools::test(filter = "plot_box")` -> PASS; commit.

### Task 2: `plot_paired()` clusters, shades, excluded, labels

**Files:** Modify `R/plot_box.R`; Test `tests/testthat/test-plot_box.R`

- Signature: `plot_paired(data, levels, condition, fill, pair, facet, y, scale, colors, show_excluded = FALSE, label = FALSE)`.
- x position = `(group_index - 1) * (n_levels + 0.8) + condition_index`; `scale_x_continuous(breaks = positions, labels = condition)`; group name `geom_text` at `y = Inf`, `vjust = 1.3` over each cluster centre, y expansion `c(0.05, 0.15)`; alpha by condition `seq(0.35, 0.95, length.out = n_levels)` with legend "condition"; lines grey per `pair` (kept), red dashed for excluded.
- [ ] Tests: existing `plot_paired` tests still pass (condition levels, line layer); `plot_paired x labels are the condition names` (`ggplot_build` x labels == rep(levels, n_groups)); `plot_paired shows excluded wells as red crosses`.
- [ ] Implement; tests PASS; commit.

### Task 3: bubble legend

**Files:** Modify `R/plot_bubble.R`; Test `tests/testthat/test-plot_bubble.R`

- Colour limits `c(-m, m)`, `m = max(abs(log2fc), na.rm = TRUE)`, `m <- 1` if not finite or 0; breaks `c(-m, 0, m)`, labels `"<-m> lower than <control>"`, `"0"`, `"+<m> higher than <control>"` (2 significant digits); title `"log2FC vs <control>"`; `guide_colourbar(order = 1)`; size `guide_legend(order = 2, override.aes = list(colour = "grey50"))`.
- [ ] Tests: `plot_bubble_fc legend names higher/lower than control` (colour scale labels match); `plot_bubble_fc builds when all log2FC are 0`.
- [ ] Implement; tests PASS; commit.

### Task 4: `plot_timecourse()` excluded + labels

**Files:** Modify `R/plot_timecourse.R`; Test `tests/testthat/test-plot_timecourse.R`

- `timecourse_prep_()` carries `.excluded` (from `plot_split_excluded_()`) in its summarise `by`. `plot_timecourse(..., show_excluded = FALSE, label = FALSE)`: medians/stars from kept; excluded = red3 dashed line + crosses + one red label per mouse at its last x; `label = TRUE` labels each kept mouse at its last x in its group colour. `plot_timecourse_summary()` drops excluded (default) only.
- [ ] Tests: default drops with message; `show_excluded` adds a red dashed line layer whose data has only excluded mice, medians unchanged; `label = TRUE` gives one label per mouse.
- [ ] Implement; tests PASS; commit.

### Task 5: docs + check

- [ ] `devtools::document()`; CLAUDE.md bullets for `plot_box()`/`plot_paired()`/`plot_timecourse()`/`plot_bubble_fc()`; `devtools::check()` 0/0; commit.
