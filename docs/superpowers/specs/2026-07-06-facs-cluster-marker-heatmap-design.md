# Design: cluster marker-expression summary + heatmap (`facs_calc_cluster_marker_medians()` / `facs_plot_cluster_heatmap()`)

**Date:** 2026-07-06
**Status:** Approved by user, not yet planned/implemented.

Companion to Stage 2's `facs_cluster_flowsom()` (`docs/superpowers/specs/2026-07-05-facs-cluster-flowsom-design.md`).
Surfaced while planning a real analysis workflow for an anti-MPO vasculitis
experiment (2 groups: Treg-induced/protected vs. control, kidney+spleen+lung,
CD4+ T cell panel) -- see the project memory entry for that experiment's
context. Not a new numbered pipeline stage; a supporting annotation utility
usable after any `facs_cluster_flowsom()` run.

---

## Goal

Nothing in the pipeline currently lets a caller **name** a metacluster --
i.e. answer "what is metacluster 5, biologically?" before interpreting or
reporting a `facs_test_cluster_abundance()` hit. This adds the standard
CyTOF/flow annotation step: per-cluster median marker expression, visualized
as a heatmap.

---

## Decisions made in brainstorming (2026-07-06)

- **Median, not mean.** Matches the diffcyt/Nowicka et al. convention this
  package already follows for cluster-level summaries, and is robust to
  outlier events within a cluster. Not made caller-selectable -- no
  identified need for a mean option, and an unused parameter is scope the
  package doesn't need yet.
- **Scaling lives in the plot function, not the calc function.** The calc
  function returns raw per-cluster medians (one clear, reusable number);
  the plot function decides how to visualize it (`scale = "zscore"` vs.
  `"raw"`). Mirrors the existing calc-computes-raw / plot-decides-view split
  (e.g. `facs_calc_cluster_freq()` returns raw `n`/`fraction`;
  `facs_plot_cluster_abundance()`'s `significant_only` filtering happens at
  plot time, not calc time).
- **Long output format**, not wide -- one row per cluster x marker. Feeds
  `ggplot2::geom_tile()` directly without a pivot step, same reasoning as
  `facs_calc_cluster_freq()`'s long shape feeding `facs_plot_cluster_abundance()`.
- **Diverging red/blue fill scale** (`ggplot2::scale_fill_gradient2()`), not
  viridis. This is a deliberate exception to the package's
  viridis-everywhere convention (`facs_plot_umap()`'s explicit
  colorblind-safe/perceptually-uniform rationale): a diverging scale
  centered at a meaningful midpoint is the field-standard way to read a
  cluster-marker heatmap (which markers are *above* vs. *below* a
  reference, not just relative magnitude). User confirmed this trade-off
  explicitly when reviewing the design.
- **No dendrogram/reordering.** Clusters and markers plot in
  factor/input order -- no new dependency (e.g. `pheatmap`,
  `ComplexHeatmap`), consistent with every other `facs_plot_*()` function
  being a plain `ggplot2` wrapper.
- **No new dependencies.** Both functions only need `dplyr`, `tidyr`,
  `purrr`, `glue`, `ggplot2`, `stats` -- all already `Imports`.

---

## Function 1: `facs_calc_cluster_marker_medians()`

Added to the *existing* `R/facs_cluster.R` (alongside `facs_cluster_flowsom()`
and `facs_calc_cluster_freq()`) -- these three functions form one
cluster-then-summarize family in one file, matching the "one file per verb
family" convention already established (`facs_reduce.R`, `facs_plot.R`).

```r
facs_calc_cluster_marker_medians(
  data,                      # facs_cluster_flowsom() output (per-event tibble)
  markers     = NULL,        # chr vector or NULL; NULL resolves to every dbl-typed column
  cluster_col = "metacluster"
)
```

### Data flow

1. Resolve `markers`: `NULL` (default) uses every `dbl`-typed column in
   `data`, identical convention to `facs_cluster_flowsom()`'s own `markers`
   resolution (including its explicit-column validation below). This
   naturally excludes `cluster_col` (a factor or integer column, never
   `dbl`).
2. `dplyr::group_by(cluster_col) |> dplyr::summarise(across(all_of(markers), stats::median))`
   -- one row per cluster, one column per marker (wide, intermediate).
3. `tidyr::pivot_longer()` the marker columns to long format: one row per
   cluster x marker.
4. Aggregation is over the *entire* input as given -- if `data` pools
   multiple tissues/groups before this call (as recommended for the real
   workflow, so metacluster identity is comparable across tissues), the
   returned medians describe each metacluster's phenotype across the whole
   dataset. No stratification argument (e.g. per-tissue medians) -- not
   needed for the annotation use case, which asks "what is this cluster,
   generally," not "how does this cluster differ by group."

### Validation (mirrors `facs_cluster_flowsom()`'s validation wording)

- `cluster_col` not in `data` -> `stop()`.
- Explicit `markers` containing a name absent from `data` -> `stop()`,
  listing the missing column(s).
- Explicit `markers` containing a non-`double`-typed column -> `stop()`,
  listing the offending column(s).
- Resolved `markers` (explicit or defaulted) is empty -> `stop()`.

### Returns

Long tibble: `{cluster_col name}` (type unchanged from input, e.g. `factor`
for `"metacluster"`), `marker` (chr), `median` (dbl). One row per cluster x
marker.

---

## Function 2: `facs_plot_cluster_heatmap()`

Added to the *existing* `R/facs_plot.R` (alongside `facs_plot_umap()` and
`facs_plot_cluster_abundance()`).

```r
facs_plot_cluster_heatmap(
  marker_medians,                    # facs_calc_cluster_marker_medians() output
  cluster_col = "metacluster",
  scale       = c("zscore", "raw")
)
```

### Data flow

1. Validate `cluster_col`, `marker`, and `median` are columns in
   `marker_medians`.
2. If `scale = "zscore"` (default): z-score the `median` column
   *per `marker`* (`(median - mean(median)) / sd(median)`, grouped by
   `marker`, across clusters) before plotting. If `scale = "raw"`: plot
   `median` unscaled.
3. `ggplot2::geom_tile(aes(x = marker, y = .data[[cluster_col]], fill = value))`
   + `ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)`
   for both `scale` modes. `midpoint = 0` is meaningful in both cases: for
   `"zscore"` it's the across-cluster mean by construction; for `"raw"` it
   aligns with the logicle/biexponential-transformed scale
   `facs_read_fcs_gated()` already returns, where 0 sits at the
   negative/background boundary (see `facs_cluster_flowsom()`'s header
   comment on this transformed scale).
4. `ggplot2::theme_minimal()`, matching the other `facs_plot_*()` functions.

### Validation

- `cluster_col`, `"marker"`, or `"median"` missing from `marker_medians` ->
  `stop()`.
- `scale` validated via `match.arg()`.

### Returns

A `ggplot` object (not printed or saved), matching every other
`facs_plot_*()` function's convention.

---

## Testing plan

Synthetic tibbles only, no fixture dependency -- same precedent as Stages
3-4 and `facs_plot_umap()`'s own tests (built from a synthetic
`metacluster` column, not real `facs_cluster_flowsom()` output).

Planned cases:

1. `facs_calc_cluster_marker_medians()`: toy per-event data with known
   per-cluster values -> exact expected medians (this is plain
   `stats::median()`, safe to pin exact numbers against, unlike Stage 4's
   `diffcyt`-internal p-values).
2. `NULL` markers resolves to every `dbl`-typed column, excluding a `chr`
   keyword column and the (non-`dbl`) `cluster_col`.
3. Validation errors: missing `cluster_col`, explicit `markers` missing
   from data, explicit `markers` non-double, empty resolved `markers`.
4. `facs_plot_cluster_heatmap()`: returns a `ggplot` object for both
   `scale = "zscore"` and `scale = "raw"`; validation errors for each
   missing required column; `match.arg()` rejects an invalid `scale` value.

---

## Open questions for implementation time (not blocking this design)

None identified -- both functions build only on already-verified base/tidyverse/
ggplot2 mechanics (`stats::median()`, `dplyr::summarise()`,
`tidyr::pivot_longer()`, `ggplot2::scale_fill_gradient2()`), unlike Stages
1/2/4's real risk of unverified third-party (`FlowSOM`, `CytoML`, `diffcyt`)
behavior.
