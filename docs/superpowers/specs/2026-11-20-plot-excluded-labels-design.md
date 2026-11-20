# Plot helpers: excluded mice, ID labels, clearer bubble legend, paired layout

Date: 2026-11-20. Size: small change (approved in conversation).

## Intent

Analysis documents need to show *why* a mouse is missing from a comparison
(excluded mice visible but clearly marked), to identify individual mice on
the plot, to read the bubble plot's colours without guessing, and to get a
readable restimulation plot.

## Scope

Plots that draw individual mice: `plot_box()`, `plot_population()`,
`plot_paired()`, `plot_timecourse()`. `plot_timecourse_summary()` and
`plot_bubble_fc()` draw no individual mice (only the default drop of excluded
rows applies to `plot_timecourse_summary()`).

## Design

1. **`show_excluded = FALSE`**. Data with an `excluded` column (from
   `db_load_experiment(keep_excluded = TRUE)` / `db_query_*()`): excluded
   rows are dropped with a message by default. `TRUE` draws them as red
   crosses with a red mouse-ID label; they never enter boxes, medians or
   p-values. `plot_timecourse()`: red dashed line per excluded mouse, one
   label at its last time point. `TRUE` without an `excluded` column = message,
   nothing to show.
2. **`label = FALSE`**. `TRUE` labels every point with the mouse ID (column
   `id` in `plot_box()` / `plot_population()`, default `mouse_id`; `pair` in
   `plot_paired()`; `subject` in `plot_timecourse()`, one label per mouse at
   its last point). Labels share the points' jitter; `ggrepel` when installed
   (Suggests), plain text otherwise.
3. **Bubble legend**. Colour gradient kept, limits symmetric around 0, title
   `log2FC vs <control>`, ends labelled `higher than <control>` (red) /
   `lower than <control>` (blue); size legend drawn in grey and placed after
   the colour bar.
4. **`plot_paired()`**: one cluster per group, conditions side by side within
   it (numeric x positions), lines join each mouse inside its cluster; x labels
   = condition names only; group name printed above its cluster. Fill = group
   colour (option `JanisHelpers.group_colors`), conditions as shades: first
   level (baseline) lightest, last full (alpha ramp, legend "condition").

## Approaches considered

- Paired: condition on the x axis with groups dodged inside (rejected by the
  user in favour of group clusters: a mouse's line stays inside its group).
- Group name *below* the cluster needs nested axes (not in ggplot2) or
  `facet_grid(switch = "x")`, which breaks the current `facet_wrap` population
  layout; builder's call: name above the cluster inside the panel.
- Excluded points as a separate layer would dodge differently from the kept
  points; builder's call: one point layer over all rows with shape/colour by
  `excluded`, boxes/stats from kept rows only.
