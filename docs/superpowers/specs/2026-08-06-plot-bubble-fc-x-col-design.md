# Design: `plot_bubble_fc()` — configurable x-axis (`x_col`)

**Date:** 2026-08-06
**Status:** Approved

---

## Goal

`plot_bubble_fc()` (see [[2026-07-10-plot-bubble-fc-design]]) currently puts
one column per non-control group on the x-axis. When there are only 2
groups (control + 1 treatment), that collapses to a single column — wasted
horizontal space. This change lets the caller redirect the x-axis to a
different dimension in that case: additional metrics from the same tissue,
or additional tissues with the same metric, instead of (or in addition to
— generalized to any column) group comparisons.

`group_col`/`control` keep their current meaning and defaults; this is an
additive change to `plot_bubble_fc()`/`calc_bubble_fc_()`, not a redesign.

---

## Public API

### New argument: `x_col`

```r
plot_bubble_fc(data, control,
               group_col = "group",
               population_col = "population",
               value_col = "value",
               x_col = "comparison",
               test = "auto",
               p_adjust_method = "BH",
               summary_fun = mean,
               pseudocount = 0.5)
```

| Value | Behavior |
|---|---|
| `"comparison"` (default) | Unchanged: one x-axis column per non-control group. Works with any number of groups (2+), exactly as today. |
| Any other column name present in `data` | One x-axis column per distinct value of that column (typically `"tissue"` or `"metric"`, but not restricted to those — generic). **Requires `group_col` to have exactly 2 distinct values** (control + one treatment); `stop()` otherwise. |

`calc_bubble_fc_()` gains the same `x_col` parameter, threaded through from
`plot_bubble_fc()`.

---

## Mechanics

Generalizes "the dimension that varies across x" from *group name* to *any
column's distinct values*.

**`x_col == "comparison"` (default):** identical to today —
`comparison_levels <- setdiff(group_levels, control)`; for each
`population x comparison` cell, `control_values`/`group_values` are pulled
from `data` filtered by `group_col`, straight from the (pre-filtered to
one metric + one tissue) `data`.

**`x_col != "comparison"`:**

1. Validate `x_col` exists in `data` and is not `group_col`/
   `population_col`/`value_col` — `stop()` if so (ambiguous role).
2. Validate `group_col` has exactly 2 distinct values in `data` — `stop()`
   otherwise: `"x_col can only be set to a column other than 'comparison'
   when group_col has exactly 2 distinct values (control + one treatment);
   found N. Use x_col = 'comparison' (default) for more groups."`
3. `comparison_group <- setdiff(group_levels, control)` — the single
   treatment group's name, used only in the plot caption (see Plot
   section).
4. `x_values <- unique(data[[x_col]])`, in first-appearance order (same
   convention already used for `population`/`comparison` levels).
5. For each `population x x_value` cell: slice `data` to
   `data[[x_col]] == x_value` **first**, then within that slice pull
   `control_values` (`group_col == control`) and `group_values`
   (`group_col == comparison_group`). Same log2fc/testable/`warning()`
   logic as today, just computed per-slice instead of straight from the
   whole (pre-filtered) `data`.
6. The existing `metric`/`tissue` >1-distinct-value guards stay as today,
   **except** whichever of `metric`/`tissue` equals `x_col` is exempted
   (that's precisely the column now allowed to vary). If `x_col` is
   neither `"metric"` nor `"tissue"` (e.g. some other keyword column),
   both guards still apply unchanged.
7. Test dispatch (`"auto"`/`"wilcox"`/`"t.test"`/`"kruskal"`) works
   identically to today but keyed on `x_value` instead of `comparison`:
   - `p_adjust_method` grouping becomes `dplyr::group_by(x_value)` (was
     `group_by(comparison)`) for the wilcox/t.test paths.
   - The `kruskal`-forced path (and `"auto"` when it would resolve to
     kruskal — moot here since `group_col` has exactly 2 values in this
     mode, so `"auto"` always resolves to `wilcox`) slices by
     `population x x_value` before running `dunn_test()`; still
     degenerates to a 2-group comparison, same as forcing
     `test = "kruskal"` today with `k == 2`.
8. The output tibble's x-axis column is **named after `x_col` itself**
   (`comparison` when default, `tissue`/`metric`/whatever otherwise) —
   `plot_bubble_fc()` builds the plot with
   `ggplot2::aes(x = .data[[x_col]])`, using the `.data` pronoun since the
   column name is held in a variable (per the package's existing tidy-eval
   convention).

---

## Plot & labeling

- `ggplot2::aes(x = .data[[x_col]])` (generalized from the hardcoded
  `comparison`).
- X-axis label: `"Comparison"` when `x_col == "comparison"` (unchanged);
  otherwise `tools::toTitleCase(x_col)` (e.g. `"Tissue"`, `"Metric"`).
  `tools` is already an Import (used elsewhere in the package).
- Caption: unchanged `"Compared to control: {control}"` when
  `x_col == "comparison"`. When `x_col != "comparison"`, becomes
  `"{comparison_group} vs control: {control}"`, naming the single
  treatment group since it's no longer visible as separate x-axis columns.
- Color/size scales, star labels, NA/untestable shape treatment, and theme
  are unchanged.

---

## Error handling (additions to the existing table)

| Condition | Behavior |
|---|---|
| `x_col` not `"comparison"` and not present in `data` | `stop()` |
| `x_col` equals `group_col`, `population_col`, or `value_col` | `stop()` — ambiguous role |
| `x_col` not `"comparison"` and `group_col` has != 2 distinct values in `data` | `stop()` |

All existing error conditions (missing `group_col`/`population_col`/
`value_col`, `control` not found, invalid `test`, `metric`/`tissue`
>1-distinct-value guards) are unchanged, except the `metric`/`tissue`
guard is skipped for whichever one equals `x_col`.

---

## Testing

Extends `tests/testthat/test-plot_bubble.R`. Existing tests for the
`x_col = "comparison"` default path are unaffected — the default
reproduces current behavior exactly, so no existing test needs updating.

**`calc_bubble_fc_()` new cases:**
- `x_col = "tissue"` with 2 groups x 2 tissues x 1 metric: correct
  log2fc/p-value computed independently per tissue slice; output tibble
  has a `tissue` column (not `comparison`) with the two tissue values as
  levels, in first-appearance order.
- `x_col = "metric"`: same shape, keyed on metric instead.
- `x_col = "tissue"` errors when `group_col` has 3 distinct values.
- `x_col = "nonexistent_col"` errors (column not in `data`).
- `x_col` set to `group_col`/`population_col`/`value_col` errors.
- `x_col = "tissue"`: the `tissue` >1-distinct-value guard is lifted, but
  the `metric` >1-distinct-value guard still fires if `data` has >1
  distinct `metric`.
- `p_adjust_method` applied per `x_value` (not globally) when
  `x_col != "comparison"` (mirrors the existing per-`comparison` test).
- `test = "kruskal"` forced with `x_col = "tissue"` still produces a
  p-value per `population x tissue` cell.

**`plot_bubble_fc()` new cases:**
- `x_col = "tissue"`: builds without error (`ggplot2::ggplot_build()`);
  caption contains both `comparison_group` and `control`; x-axis label is
  `"Tissue"`.

---

## Dependencies

No new dependencies. `tools` is already used elsewhere in the package
(e.g. `report_knit_*()`).

---

## Roxygen / CLAUDE.md updates

- `plot_bubble_fc()`'s roxygen gains `@param x_col` documenting the
  default and the "exactly 2 groups" requirement for non-default values.
- CLAUDE.md's `plot_bubble_fc()` data-structure entry gets a short note on
  `x_col` and its group-count requirement.
- "Known check output" gets a bullet for any new `x_col`/`x_value`
  variable-binding notes surfaced by `devtools::check()`, matching the
  existing bullets for other domains (expected: none, since `x_col`
  access goes through `.data[[x_col]]`/`[[x_col]]` subsetting rather than
  bare column names).
