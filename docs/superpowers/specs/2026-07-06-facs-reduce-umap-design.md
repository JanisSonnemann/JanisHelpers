# Design: Stage 3 — UMAP dimensionality reduction

**Date:** 2026-07-06
**Status:** Design approved by user, pending write-up review.
**Roadmap:** `docs/superpowers/specs/2026-07-03-facs-unsupervised-analysis-roadmap.md`
(Stage 3 of 4). Follows Stage 1 (`facs_read_fcs_gated()`, merged) and Stage 2
(`facs_cluster_flowsom()` / `facs_calc_cluster_freq()`, PR #1 open, not yet
merged).

---

## Goal

Compute a 2D UMAP embedding of gated single-cell events (Stage 1's or Stage
2's output) for visualization, and provide a plotting helper that renders
the embedding colored by metacluster, group, or marker expression.

Two new exported functions:
- `facs_reduce_umap()` in a new `R/facs_reduce.R`
- `facs_plot_umap()` in a new `R/facs_plot.R`

New dependencies: `uwot` (UMAP implementation) and `ggplot2` (plotting),
both CRAN, added to `Imports`.

---

## `facs_reduce_umap()`

```r
facs_reduce_umap(data, markers = NULL, max_events = NULL,
                  n_neighbors = 15, min_dist = 0.1, seed = NULL)
```

### Parameters

- **`data`**: tibble shaped like `facs_read_fcs_gated()`'s or
  `facs_cluster_flowsom()`'s output — one row per event, `dbl` marker
  columns, optional `chr` keyword columns, must include `file_name`.
- **`markers`**: character vector of column names to embed on. `NULL`
  (default) uses every `dbl`-typed column in `data` — identical convention
  to `facs_cluster_flowsom()`'s `markers` argument, including its
  validation rules.
- **`max_events`**: integer or `NULL` (default). If set, downsamples the
  combined tibble to (approximately) this many total rows before running
  UMAP, **stratified per `file_name`**: each sample contributes
  `floor(max_events / n_samples)` rows, with the remainder (from
  `max_events %% n_samples`) distributed one extra row each to the first
  samples (by `file_name` sort order) that have enough rows to supply it.
  A sample with fewer rows than its share contributes all of its rows (no
  error, no attempt to redistribute its shortfall to other samples) — so
  the returned row count can be slightly under `max_events` when at least
  one sample is short. If `n_samples > max_events` (share rounds to 0),
  only the first `max_events` samples (by `file_name` sort order) get one
  row each; the rest contribute none. If `max_events >= nrow(data)`,
  downsampling is skipped silently. Rationale:
  UMAP on a naive random pool would let a sample with more input events
  (e.g. higher live-cell yield) dominate the embedding independent of
  biology; stratifying keeps every sample visible.
- **`n_neighbors`**, **`min_dist`**: passed through to `uwot::umap()`
  unchanged (default `15` / `0.1`, `uwot`'s own defaults). Exposed because
  these are the two hyperparameters researchers most often need to retune
  per panel/dataset to get a readable embedding.
- **`seed`**: integer or `NULL`. If set, seeds both the stratified draw (so
  the same subset of rows is chosen) and `uwot::umap()`'s internal RNG (so
  the embedding itself is reproducible).

### Behavior

1. Resolve `markers` (default-all-double or explicit), validating exactly
   like `facs_cluster_flowsom()`: error listing any name absent from
   `data`, error listing any explicitly-supplied non-double column, error
   if resolved `markers` has fewer than 2 columns, error listing any
   marker column containing `NA`.
2. If `max_events` is set and less than `nrow(data)`, perform the
   stratified-per-`file_name` draw described above.
3. Run `uwot::umap()` on the resolved marker matrix with `n_neighbors`,
   `min_dist`, and `seed`.
4. Append `UMAP1`, `UMAP2` (`dbl`) to the (possibly downsampled) tibble and
   return it.

### Returns

The input tibble — or its downsampled subset, if `max_events` triggered
downsampling — with `UMAP1` and `UMAP2` columns appended. Every returned
row has a real embedding (no `NA` coordinates): rows excluded by
downsampling are dropped from the output entirely, not kept with `NA`.

### Errors

Same four marker-validation errors as `facs_cluster_flowsom()` (missing
column, non-double column, fewer than 2 resolved markers, `NA` in a marker
column), plus: error if `max_events` is supplied and is not a positive
integer.

---

## `facs_plot_umap()`

```r
facs_plot_umap(data, color_by = "metacluster", facet_by = NULL)
```

### Parameters

- **`data`**: tibble shaped like `facs_reduce_umap()`'s output — must
  contain `UMAP1` and `UMAP2`.
- **`color_by`**: character; column name in `data` to color points by.
  Default `"metacluster"`. Any column works — a marker (continuous
  expression gradient), `group`, `file_name`, etc.
- **`facet_by`**: character or `NULL` (default). Column name in `data` to
  facet panels by (e.g. `"group"`) via `ggplot2::facet_wrap()`. `NULL`
  produces a single panel.

### Behavior

1. Error if `UMAP1`/`UMAP2` are not both present in `data`.
2. Error if `color_by` is not a column in `data`.
3. Error if `facet_by` is supplied and is not a column in `data`.
4. Build a `ggplot2` scatter: `geom_point()` with small size and partial
   transparency (tuned for event-density overplotting at typical FACS
   event counts), `UMAP1`/`UMAP2` on the axes (labeled `"UMAP1"` /
   `"UMAP2"`), `theme_minimal()` with recessive gridlines.
5. Color scale chosen by `color_by`'s type: `ggplot2::scale_color_viridis_c()`
   if `data[[color_by]]` is `double` (continuous expression gradient),
   `ggplot2::scale_color_viridis_d()` otherwise (discrete — factor/character,
   e.g. `metacluster`, `group`). Both are colorblind-safe and perceptually
   uniform, and ship inside `ggplot2` (via `viridisLite`) — no extra
   dependency beyond `ggplot2` itself.
6. If `facet_by` is set, add `ggplot2::facet_wrap()` on that column (via
   `stats::as.formula(paste0("~", facet_by))`, avoiding an `rlang`
   dependency for tidy-eval faceting).

### Returns

The `ggplot` object (not printed or saved) — caller can further customize
with `+` or save with `ggplot2::ggsave()`.

### Errors

Missing `UMAP1`/`UMAP2`, missing `color_by` column, missing `facet_by`
column (when supplied).

---

## File and dependency changes

- New `R/facs_reduce.R`: `facs_reduce_umap()`.
- New `R/facs_plot.R`: `facs_plot_umap()`.
- `DESCRIPTION` `Imports`: add `uwot`, `ggplot2`.
- Both functions namespaced (`uwot::umap()`, `ggplot2::geom_point()`, etc.)
  per the package's namespacing convention.

---

## Testing

Mirrors Stage 2's test approach: `tests/fixtures/Treg.wsp` +
`tests/fixtures/Treg/*.fcs`, piped `facs_read_fcs_gated()` →
`facs_cluster_flowsom()` → `facs_reduce_umap()` → `facs_plot_umap()`,
skipped with `skip_if_not(dir.exists(fcs_dir), ...)` when the fixture is
absent.

`test-facs_reduce.R`:
- Appends `UMAP1`/`UMAP2` (`dbl`) columns; row count unchanged when
  `max_events` is `NULL`.
- Default `markers` resolves to every `dbl` column; explicit `markers`
  subset works.
- Marker-validation errors: missing column, non-double column, fewer than
  2 resolved columns, `NA` in a marker column.
- `max_events` downsampling: total row count matches `max_events` (when
  less than `nrow(data)`), each `file_name` is represented and no sample
  is entirely dropped (given enough rows to supply its share), skipped
  silently when `max_events >= nrow(data)`.
- `max_events` not a positive integer errors.
- Same `seed` on two calls produces identical `UMAP1`/`UMAP2` (and the
  same downsampled rows, when `max_events` is also set).

`test-facs_plot.R`:
- Returns a `ggplot` object.
- Errors when `UMAP1`/`UMAP2` missing, when `color_by` column missing,
  when `facet_by` column missing.
- Continuous `color_by` (a marker) selects `scale_color_viridis_c`;
  discrete `color_by` (`metacluster`, `group`) selects
  `scale_color_viridis_d` — checked via the built plot's scale class
  (`ggplot2::ggplot_build()`).
- `facet_by` produces the expected number of panels (checked via the
  built plot's layout).

---

## Open items deferred to Stage 4 or later

- Domain prefix for Stage 4 (diffcyt-based differential abundance) is
  still undecided — noted in the roadmap, not affected by this design.
