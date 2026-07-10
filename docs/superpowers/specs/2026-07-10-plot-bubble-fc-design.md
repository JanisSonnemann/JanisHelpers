# Design: `plot_bubble_fc()` — group-vs-control fold-change bubble plot

**Date:** 2026-07-10
**Status:** Approved

---

## Goal

A quick visual overview of how many populations shift between one or more
experimental groups and an explicitly named control, in a single glance:
one column per group-vs-control comparison, one row per population, a
bubble at each intersection sized and colored by fold-change, annotated
with a significance-star summary from a named statistical test.

This is a new capability, not a variant of an existing one:
`facs_calc_log2fc()`/`facs_calc_diff()` (see
[[2026-07-03-facs-calc-log2fc-design]]) compare paired restimulation
conditions *within the same mouse*; this compares *independent* groups of
mice against a control, and renders the result as a plot rather than a
tibble.

---

## New domain: `plot_`

No existing domain fits: `facs_` is import/calc, `analysis_` produces
`gt`/`gtsummary` tables (not `ggplot` objects), the rest are unrelated.
`CLAUDE.md`'s domain table gets a new row:

| Domain | Purpose |
|---|---|
| `plot_` | ggplot2 visualization of FACS/analysis data |

File: `R/plot_bubble.R` (`domain_verb` pattern — `bubble` names the chart
type, matching how `facs_calc_log2fc` names the comparison it computes).

---

## Architecture

One exported function, `plot_bubble_fc()`, backed by one unexported helper,
`calc_bubble_fc_()`, that does all the stats and returns a tidy tibble
(testable independently of the plot). This mirrors the existing
`calc_restim_proportions_()` / `facs_calc_log2fc()` split.

```
data (long tibble, pre-filtered to one metric + one tissue)
  -> calc_bubble_fc_()   [stats: fold-change, test dispatch, p-adjust, stars]
  -> plot_bubble_fc()    [ggplot2: bubble geom + star labels + theme + caption]
```

New dependencies: `ggplot2` and `stats` (for `wilcox.test`/`t.test`/
`kruskal.test`/`p.adjust`, all namespaced `stats::`), added to `Imports`.
`rstatix` (already a dependency) supplies `dunn_test()` for the Kruskal
post-hoc.

---

## Public API

### `plot_bubble_fc(data, control, group_col = "group", population_col = "population", value_col = "value", test = "auto", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5)`

| Argument | Type | Default | Notes |
|---|---|---|---|
| `data` | tibble | — | Long tibble already filtered to one `metric` and one `tissue` (e.g. `facs_read_wsp(...)$data` after `dplyr::filter()`). Must contain `group_col`, `population_col`, `value_col`. |
| `control` | chr | — | Value in `group_col` treated as the baseline. Every other distinct value is compared against it and named in the plot caption. |
| `group_col` | chr | `"group"` | Column holding group/treatment labels. |
| `population_col` | chr | `"population"` | Column holding population names (plot rows). |
| `value_col` | chr | `"value"` | Column holding the numeric measurement (already filtered to one metric). |
| `test` | chr | `"auto"` | One of `"auto"`, `"wilcox"`, `"t.test"`, `"kruskal"`. See dispatch table below. |
| `p_adjust_method` | chr | `"BH"` | Passed to `stats::p.adjust()` (wilcox/t.test path) or `rstatix::dunn_test(p.adjust.method = )` (kruskal path). Applied per comparison column, across all populations in that column. |
| `summary_fun` | function | `mean` | Aggregator used on each group's `value_col` vector when computing fold-change (not used for the test itself, which always runs on the raw per-mouse values). |
| `pseudocount` | dbl | `0.5` | Added to both the group and control `summary_fun()` results before dividing, avoiding `log2(0)`/`Inf` when a population is entirely absent in one arm. Same role as `facs_calc_log2fc()`'s `pseudocount`. |

**Returns:** a `ggplot` object (visible return — this is a primary-purpose
plot function, not a data-import side effect, so no `invisible()` per the
`facs_import_wsp()`/`meta_read()` convention).

---

## Test dispatch

Let *k* = number of distinct values in `data[[group_col]]` (including
`control`).

| `test` | Behavior |
|---|---|
| `"auto"` | `k == 2`: pairwise `stats::wilcox.test()` per comparison column. `k > 2`: one `stats::kruskal.test()` per population across *all* groups at once, followed by `rstatix::dunn_test(p.adjust.method = p_adjust_method)`, keeping only the control-vs-group rows. |
| `"wilcox"` | Forced pairwise `stats::wilcox.test(control_values, group_values)` per comparison column, regardless of *k*. |
| `"t.test"` | Forced pairwise `stats::t.test(control_values, group_values)` per comparison column, regardless of *k*. |
| `"kruskal"` | Forced omnibus `stats::kruskal.test()` + `rstatix::dunn_test()` post-hoc even when `k == 2` (degenerates to the 2-group case, statistically equivalent to wilcox but computed via the Kruskal/Dunn path). |

For the wilcox/t.test paths, `p.adjust(p, method = p_adjust_method)` is
applied per comparison column (across all populations tested in that
column) after the raw p-values are collected. For the kruskal path,
`rstatix::dunn_test()`'s own `p.adjust.method` argument does this per
population (across all its pairwise comparisons, including non-control
ones); only the control-vs-group rows are kept afterward, so those already
carry an adjusted p-value from the full comparison set.

---

## Internal helper: `calc_bubble_fc_(data, control, group_col, population_col, value_col, test, p_adjust_method, summary_fun, pseudocount)`

Not exported. Returns a tidy tibble: one row per `population_col` value x
non-control group, columns `population`, `comparison` (the group name),
`log2fc`, `p_value` (adjusted), `stars`.

**Processing steps, in order:**

1. Validate `group_col`/`population_col`/`value_col` exist in `data`;
   `stop()` listing any missing.
2. Validate `control` is present in `data[[group_col]]`; `stop()` if not.
3. If `data` has a `metric` column with >1 distinct value, or a `tissue`
   column with >1 distinct value: `stop()` — caller must pre-filter to one
   metric and one tissue.
4. For each `population_col` value:
   a. `control_values <- value_col` where `group_col == control`.
   b. For each non-control group: `group_values <- value_col` where
      `group_col == that group`.
   c. `log2fc = log2((summary_fun(group_values) + pseudocount) / (summary_fun(control_values) + pseudocount))`.
   d. If either `control_values` or `group_values` has fewer than 2
      observations: `warning()` naming the population/group combo, `p_value`
      and `stars` become `NA` for that cell (point still drawn via
      `log2fc`, but with no star and, in the plot, a distinct "untested"
      visual treatment — see Plot section).
5. Compute raw p-values per the test-dispatch table above.
6. Apply `p_adjust_method` per comparison column (wilcox/t.test) or via
   `rstatix::dunn_test()`'s own argument (kruskal), as described above.
7. Map adjusted p-values to stars: `< 0.001` -> `"***"`, `< 0.01` -> `"**"`,
   `< 0.05` -> `"*"`, otherwise `""` (fixed thresholds, not configurable).
8. Return the tidy tibble, with `population` and `comparison` as factors
   ordered by first appearance in `data` (population: order in
   `data[[population_col]]`; comparison: order of first appearance in
   `data[[group_col]]`, excluding `control`).

---

## Plot

- `ggplot2::ggplot(data, aes(x = comparison, y = population))`.
- `geom_point(aes(color = log2fc, size = abs(log2fc)))`.
- Color: diverging scale via `ggplot2::scale_color_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0)` — blue↔red diverging pair with a near-white neutral midpoint at zero, per the dataviz skill's palette (`references/palette.md`). Light-mode only; this is a static plot, not a themed web surface.
- Size: `ggplot2::scale_size_continuous(range = c(2, 10))`, own legend (not merged with color — two separate legends, simpler and avoids guide-merging edge cases with a diverging color scale).
- Significance stars: `geom_text(aes(label = stars), vjust = -1.2)` (or similar offset so the star sits above the bubble rather than overlapping it).
- Cells with `NA` p-value (untested, <2 observations in one arm): bubble still drawn from `log2fc` (or omitted entirely if `log2fc` is also `NA`), no star, and visually flagged — e.g. a dashed/hollow outline via a secondary `shape` mapping — so "not significant" and "not testable" are never visually confused.
- Theme: `ggplot2::theme_minimal()` base, then hairline gridlines
  (`panel.grid = element_line(color = "#e1e0d9", linewidth = 0.3)`, major
  only), muted axis text (`#898781`), primary-ink axis titles (`#0b0b0b`),
  no chart junk (no vertical gridlines needed since x is categorical).
- Caption naming the control group:
  `ggplot2::labs(caption = glue::glue("Compared to control: {control}"))`.
- Axis labels: x = `group_col`'s name (or a fixed `"Comparison"`), y =
  `population_col`'s name (or a fixed `"Population"`).

---

## Error handling

| Condition | Behavior |
|---|---|
| `group_col`/`population_col`/`value_col` not in `data` | `stop()`, lists missing columns |
| `control` not found in `data[[group_col]]` | `stop()` |
| `data` has a `metric` or `tissue` column with >1 distinct value | `stop()` — caller must pre-filter |
| A population/group combo has <2 observations in either arm | `warning()` naming the combo; `p_value`/`stars` -> `NA` for that cell |
| `test` not one of `"auto"`/`"wilcox"`/`"t.test"`/`"kruskal"` | `stop()` |

---

## Testing

New file `tests/testthat/test-plot_bubble.R`, synthetic tibbles only (no
fixture file needed — unlike `.wsp`-dependent tests, this only needs a
small hand-built long tibble with `group`, `population`, `value`
columns).

**`calc_bubble_fc_()` (tested directly, unlike `calc_restim_proportions_()`,
since it's the natural unit boundary for all the stats logic):**
- Correct `log2fc` sign and magnitude for a simple 2-group case (one group
  vs `control`)
- `pseudocount` avoids `-Inf`/`NaN` when a population's control-group
  values are all `0`
- `k == 2` with `test = "auto"` dispatches to `wilcox.test()`
- `k > 2` with `test = "auto"` dispatches to `kruskal.test()` +
  `rstatix::dunn_test()`, and only control-vs-group rows survive
- `test = "wilcox"`/`"t.test"` forces pairwise tests even when `k > 2`
- `test = "kruskal"` forces the omnibus+Dunn path even when `k == 2`
- `p_adjust_method` changes the resulting `p_value` (compare `"BH"` vs
  `"none"` on a case with >1 population)
- Star thresholds: `0.001`/`0.01`/`0.05` boundaries map correctly, `""`
  for non-significant
- <2 observations in one arm: `warning()` raised, `p_value`/`stars` are
  `NA` for that cell, other cells unaffected
- `population`/`comparison` ordering matches first-appearance order in
  input data, and `control` never appears as a `comparison` value
- Errors: missing `group_col`/`population_col`/`value_col`; `control` not
  present in `group_col`; `metric` or `tissue` column present with >1
  distinct value; invalid `test` value

**`plot_bubble_fc()` (thin wrapper — verify plot construction, not
re-verify stats already covered above):**
- Returns a `ggplot` object
- `caption` in the built plot contains the `control` value
- Smoke test: builds without error for both the `k == 2` and `k > 2`
  cases, and when some cells are `NA` (untested)

---

## Dependencies

Add `ggplot2` and `stats` to `Imports`. `rstatix`, `dplyr`, `tidyr`,
`glue`, `purrr` already present and reused.

---

## Roxygen

`plot_bubble_fc()` gets full `@param`/`@returns`/`@export`. `@examples`
uses `\dontrun{}` (needs realistic multi-column input not easily
constructed inline, consistent with `facs_calc_log2fc()`). `calc_bubble_fc_()`
is unexported — no roxygen block required, but gets a one-line comment
noting it's the stats engine behind `plot_bubble_fc()`.

---

## CLAUDE.md updates

- New `plot_` row in the domain table.
- New "Data structures" entry for `plot_bubble_fc()` describing its input
  shape and output (a `ggplot` object, not a tibble — first non-tibble
  return in the package, worth calling out explicitly).
- `Imports` list updated with `ggplot2`, `stats`.
- "Known check output" gets a new bullet for any `plot_bubble_fc`/
  `calc_bubble_fc_`-specific variable-binding notes surfaced by
  `devtools::check()`, matching the existing bullets for other domains.
