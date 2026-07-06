# Design: Stage 4 -- differential cluster abundance testing (`facs_test_cluster_abundance()`)

**Date:** 2026-07-06
**Status:** Approved by user, not yet planned/implemented.

Stage 4 of the roadmap at
`docs/superpowers/specs/2026-07-03-facs-unsupervised-analysis-roadmap.md`.
Consumes Stage 2's `facs_calc_cluster_freq()` output directly. See
`docs/superpowers/specs/2026-07-05-facs-cluster-flowsom-design.md` for that
function's exact shape and the `meta_` domain's `meta_annotate()` in
`R/meta_wrangle.R` for the join convention reused here.

---

## Goal

Test whether Stage 2's per-sample cluster/metacluster frequencies differ
between experimental groups, using Bioconductor's `diffcyt` package (the
standard, peer-reviewed differential-abundance framework for CyTOF/flow
cytometry), rather than a hand-rolled GLMM or the package's existing
`rstatix`-based `analysis_` domain.

---

## Decisions made in brainstorming (2026-07-06)

- **Statistical framework:** full `diffcyt` dependency (not a custom
  `lme4`/`edgeR` implementation bypassing it). Trade-off accepted knowingly:
  this adds a real translation layer between our tidy tibbles and
  `diffcyt`'s `SummarizedExperiment`-based container, and a heavier
  Bioconductor dependency chain, in exchange for a peer-reviewed,
  citable methodology and three ready-made test methods instead of one
  we'd own and validate ourselves.
- **Domain prefix:** `facs_` (not `analysis_`) -- Stage 4 is the fourth
  sequential stage of one pipeline that consumes the previous stage's
  output directly, consistent with `facs_cluster_flowsom()`,
  `facs_reduce_umap()`, `facs_calc_cluster_freq()`.
- **Test method:** caller-selectable via `method = c("glmm", "edgeR",
  "voom")`, matching `diffcyt`'s own flexibility rather than locking in
  one method.
- **Design/group specification:** caller passes an optional `meta` tibble
  + explicit fixed/random column names (mirrors `meta_annotate()`'s
  existing convention), not a raw R formula -- consistent with the rest of
  the package's explicit-args-only philosophy (no formula-object API
  surface, no auto-detection).
- **Contrast scope:** `fixed` may have more than 2 levels. An explicit
  `ref_level` argument (default: the column's first factor level) picks
  the reference; every other level is tested against it as its own
  contrast. Chosen over restricting to 2-level factors only, to support
  e.g. 3-arm dose-response designs without forcing callers to pre-filter.
- **Companion visualization:** ship alongside a plotting function in the
  same PR, mirroring Stage 3's `facs_reduce_umap()` + `facs_plot_umap()`
  pairing.

---

## Function 1: `facs_test_cluster_abundance()`

New file `R/facs_test.R`.

```r
facs_test_cluster_abundance(
  freq_data,                 # facs_calc_cluster_freq() output
  meta       = NULL,         # optional tibble to left-join via meta_annotate()
  fixed,                     # chr, name of the fixed-effect column (e.g. "group")
  random     = NULL,         # chr vector or NULL, random-effect column(s) (e.g. "mouse_ID"); GLMM only
  by         = "mouse_ID",   # join key(s) forwarded to meta_annotate(); ignored if meta is NULL
  cluster_col = "metacluster", # column in freq_data identifying the cluster; matches
                                # facs_calc_cluster_freq()'s own cluster_col default
  ref_level  = NULL,         # reference level of `fixed`; default = first factor level
  method     = c("glmm", "edgeR", "voom")
)
```

### Data flow

1. If `meta` is supplied, left-join it onto `freq_data` via the package's
   existing `meta_annotate(freq_data, meta, by)` (reuses its validation
   and unmatched-key warnings rather than reimplementing a join). If
   `meta` is `NULL`, `fixed`/`random` are assumed already present in
   `freq_data` -- which they can be for free today, since keyword columns
   requested in Stage 1's `facs_read_fcs_gated(keywords = ...)` survive
   through Stage 2's clustering and `facs_calc_cluster_freq()`'s
   constant-within-`file_name` passthrough.
2. Recompute each sample's total event count as `sum(n)` grouped by
   `file_name` (needed for the binomial GLMM's
   `cbind(cluster_n, total - cluster_n)` response; `facs_calc_cluster_freq()`
   drops its own `total` column after computing `fraction`, so this is
   recomputed here, not carried through).
3. Pivot `n` wide into a clusters x samples count matrix (`cluster_col`
   values as rows, `file_name` as columns).
4. Build an `experiment_info` tibble: one distinct row per `file_name`,
   holding the `fixed` and `random` columns.
5. Wrap into a `SummarizedExperiment` (assay = count matrix, `rowData` =
   cluster IDs, experiment info attached per `diffcyt`'s expected slot --
   see "Known implementation risk" below) and hand off to
   `diffcyt::createFormula()` (GLMM) or `createDesignMatrix()`
   (edgeR/voom), plus `createContrast()`, then the chosen `testDA_*()`.
6. One `testDA_*()` call per non-reference level of `fixed` (relative to
   `ref_level`), row-bound into one result tibble with a `contrast` column
   identifying each comparison.

### Method / random-effect interaction

`testDA_edgeR()`/`testDA_voom()` do not support random effects in
`diffcyt`'s exposed API -- only `testDA_GLMM()` does. If `method != "glmm"`
and `random` is non-`NULL`, `facs_test_cluster_abundance()` errors rather
than silently dropping the random-effect column:

> "random effects are only supported with method = 'glmm'; pass random =
> NULL, or fold this column into `fixed`, for edgeR/voom."

### Multiple-testing correction

Each `testDA_*()` call returns its own BH-adjusted `p_adj`, computed
across clusters *within that one contrast*. `facs_test_cluster_abundance()`
passes this through as-is (documented explicitly in its `@returns`) rather
than silently recomputing a joint correction across contrasts. A caller
testing a multi-level `fixed` factor who wants correction across all
contrasts combined applies `p.adjust()` to the combined `p_val` column
themselves.

### Validation (fail fast, mirroring Stage 2's FlowSOM guards)

- `fixed`, `random` (if supplied), or `cluster_col` missing from the
  (possibly joined) `freq_data` -> `stop()`.
- `fixed` column has fewer than 2 levels -> `stop()`.
- `ref_level` supplied but not among `fixed`'s levels -> `stop()` listing
  the valid levels.
- Any `NA` in `fixed`/`random` after the `meta` join (i.e. an unmatched
  join key `meta_annotate()` already warned about) -> `stop()`, since
  `diffcyt` would otherwise fail opaquely on `NA` factor levels.
- `method` validated via `match.arg()`.

### Returns

Tibble, one row per cluster x contrast: `{cluster_col name}`, `contrast`,
`p_val`, `p_adj`, plus whatever effect-size column the chosen `method`
provides (documented per-method in the `@returns` tag rather than forced
into one artificial shared name, since GLMM/edgeR/voom do not report
identical statistics). Returned visibly (not `invisible()`) -- the tibble
is the primary output the caller consumes, not a side effect of an
import/write operation.

### Known implementation risk

`diffcyt` is **not currently installed** in this environment (confirmed:
`requireNamespace("diffcyt", quietly = TRUE)` is `FALSE` as of this
writing) and needs `BiocManager::install("diffcyt")` before implementation
or testing can start. Its exact `SummarizedExperiment`/`metadata()` slot
layout for holding `experiment_info`, and its per-method output column
names (`p_val`/`p_adj`/`logFC` naming differs between GLMM, edgeR, and
voom), must be verified against the real installed version at
implementation time, not assumed from documentation memory -- same
category of risk as Stage 1's `flowjo_to_gatingset()` note and Stage 2's
three documented `FlowSOM`/`ConsensusClusterPlus` upstream bugs.

---

## Function 2: `facs_plot_cluster_abundance()`

Added to the *existing* `R/facs_plot.R` (alongside `facs_plot_umap()`) --
Stage 3 established "one file per verb" (`facs_reduce.R` for
`facs_reduce_*`, `facs_plot.R` for `facs_plot_*`), so a second
`facs_plot_` function belongs there, not in a new file.

```r
facs_plot_cluster_abundance(
  freq_data,
  test_result       = NULL,   # optional facs_test_cluster_abundance() output
  group_col,                  # column in freq_data to plot on the x-axis
  cluster_col       = "metacluster",
  significant_only  = FALSE,
  p_adj_threshold   = 0.05
)
```

`ggplot2` boxplot + jitter of `fraction` by `group_col`, faceted by
`cluster_col`. If `test_result` is supplied and `significant_only = TRUE`,
filters facets to clusters with any `p_adj <= p_adj_threshold` in
`test_result`.

---

## Dependencies

New `Imports`: `diffcyt` (Bioconductor) and `SummarizedExperiment`
(Bioconductor -- needed directly since the `SummarizedExperiment` container
is built by hand rather than via `diffcyt::prepareData()`). `lme4`/`edgeR`/
`limma` are pulled in transitively by `diffcyt` but never called directly,
so no direct `Imports` entry for them. `ggplot2` is already a dependency
(Stage 3).

---

## Testing plan

Synthetic data only, no fixture dependency -- following Stage 3's
precedent (its tests build a synthetic `metacluster` column rather than
depending on real `facs_cluster_flowsom()` output). `facs_test_cluster_abundance()`
operates on a plain tibble shape, so there is no reason to route through
the real `.wsp`/`.fcs` fixture chain (`tests/fixtures/Treg.wsp` +
`tests/fixtures/Treg/*.fcs`) at all.

Planned cases:

1. Two-group synthetic data with a deliberate abundance shift in one
   cluster -> that cluster's `p_val` small, others not (sanity-direction
   test, not a hardcoded exact p-value -- `diffcyt`'s internals are not
   something to pin exact numbers against).
2. 3-level `fixed` factor -> correct number of `contrast` rows produced
   (`(nlevels - 1) x n_clusters`), and changing `ref_level` changes which
   comparisons appear.
3. `random` + `method = "glmm"` runs without error; `random` +
   `method = "edgeR"`/`"voom"` errors with the documented message.
4. Validation errors: missing `fixed`/`random`/`cluster_col`, `fixed` with
   fewer than 2 levels, invalid `ref_level`, `NA` in `fixed`/`random`
   after an unmatched `meta` join.
5. Equivalence: passing `meta` + `by` produces the same result as when
   those columns are already present in `freq_data` (the
   keyword-passthrough case) -- checks the two paths agree.
6. `facs_plot_cluster_abundance()`: returns a `ggplot` object;
   `significant_only = TRUE` correctly reduces faceted clusters given a
   synthetic `test_result`.

---

## Open questions for implementation time (not blocking this design)

- Exact `SummarizedExperiment` construction call and `experiment_info`
  attachment point in the installed `diffcyt` version (see "Known
  implementation risk" above) -- verify against real package docs/source
  once installed, do not assume.
- Exact effect-size column name(s) per method (`logFC` for edgeR/voom;
  GLMM's equivalent) -- document precisely once verified against the
  installed version's actual output.
