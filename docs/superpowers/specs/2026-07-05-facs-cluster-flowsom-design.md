# Design: `facs_cluster_flowsom()` + `facs_calc_cluster_freq()` — FlowSOM clustering

**Date:** 2026-07-05
**Status:** Approved

This is Stage 2 of the unsupervised FACS analysis pipeline described in
`docs/superpowers/specs/2026-07-03-facs-unsupervised-analysis-roadmap.md`.
It consumes Stage 1's (`facs_read_fcs_gated()`) per-event tibble and
produces cluster assignments and per-sample cluster frequencies for Stage
3 (UMAP) and Stage 4 (diffcyt differential abundance) to consume.

---

## Goal

Run FlowSOM (self-organizing map + consensus metaclustering) on Stage 1's
raw single-cell events to discover phenotypes, without requiring the
caller to hand-tune a marker list beyond what Stage 1 already selected.
Also compute per-sample cluster frequencies in a shape ready for
differential abundance testing (Stage 4), with zero counts represented
explicitly rather than left as missing rows.

---

## Architecture

New file `R/facs_cluster.R`, containing two exported functions:
`facs_cluster_flowsom()` and `facs_calc_cluster_freq()`.

New dependency: `FlowSOM` (Bioconductor, non-CRAN — same treatment as
`CytoML`/`flowWorkspace`: listed in `Imports`, install requirement
documented in the README, no `remotes::`/`BiocManager::` calls inside
function bodies).

---

## Public API

### `facs_cluster_flowsom()`

```r
facs_cluster_flowsom(
  data,
  markers = NULL,
  grid_xdim = 10,
  grid_ydim = 10,
  n_metaclusters = 10,
  seed = NULL
)
```

| Argument | Type | Default | Notes |
|---|---|---|---|
| `data` | tibble | — | Stage 1's per-event tibble (`file_name` + marker columns + keyword columns), or any tibble with the same shape |
| `markers` | chr vector | `NULL` | Column names to cluster on. Default `NULL` auto-selects every `dbl`-typed column in `data` (markers are `dbl`, keyword columns are `chr` per Stage 1's design — this split requires no extra bookkeeping). Pass an explicit vector to restrict to a subset (e.g. to exclude a viability dye carried along for gating but not wanted for clustering). |
| `grid_xdim` | int | `10` | SOM grid width (number of nodes) |
| `grid_ydim` | int | `10` | SOM grid height |
| `n_metaclusters` | int | `10` | Target consensus metacluster count |
| `seed` | int | `NULL` | If set, seeds SOM training and consensus metaclustering for reproducible assignments |

#### Output

`data` with two columns appended:

| Column | Type | Notes |
|---|---|---|
| `cluster` | int | Raw SOM node assignment (1 to `grid_xdim * grid_ydim`) — kept for diagnostics / re-metaclustering without re-running the SOM step |
| `metacluster` | factor | Consensus metacluster assignment (1 to `n_metaclusters`) — the level downstream stages (UMAP coloring, Stage 4 GLMM grouping) consume. Factor, not integer, since it is categorical by nature. |

Returned visibly (not `invisible()`) — primary deliverable, same convention
as `facs_read_fcs_gated()`.

#### Processing steps

1. Resolve marker columns: use `markers` if supplied (error if any name is
   absent from `data`, or if any supplied name is not a `double`-typed
   column in `data`), else default to every `dbl`-typed column in `data`.
   Error if the resolved `markers` has fewer than 2 columns.
2. Error if any selected marker column contains `NA`, listing the
   affected column name(s).
3. Error if `n_metaclusters < 3`, or if `n_metaclusters` is not strictly
   less than `grid_xdim * grid_ydim` (i.e. reject `>=`, not just `>` --
   can't have as many or more metaclusters than SOM nodes).
4. Run `FlowSOM::FlowSOM()` (or the equivalent lower-level
   `ReadInput`/`BuildSOM`/`metaClustering_consensus` call chain — exact
   signature to be confirmed against the real fixture at implementation
   time; FlowSOM's public API has changed across package versions, same
   category of risk Stage 1 flagged for `flowjo_to_gatingset()`) on the
   selected marker matrix, with `xdim = grid_xdim`, `ydim = grid_ydim`,
   target metacluster count `n_metaclusters`, and `seed`.
5. Append `cluster` (raw SOM node) and `metacluster` (consensus grouping,
   as a factor) columns to `data`, preserving row order.
6. Return the augmented tibble.

---

### `facs_calc_cluster_freq()`

```r
facs_calc_cluster_freq(data, cluster_col = "metacluster")
```

| Argument | Type | Default | Notes |
|---|---|---|---|
| `data` | tibble | — | Output of `facs_cluster_flowsom()` (or any tibble with `file_name` and `cluster_col`) |
| `cluster_col` | chr | `"metacluster"` | Name of the column to aggregate on — pass `"cluster"` to get raw-SOM-node frequencies instead |

#### Output

Long tibble, one row per `file_name` x `cluster_col` value:

| Column | Type | Notes |
|---|---|---|
| `file_name` | chr | Same convention as elsewhere in the package |
| `<cluster_col>` | as input | e.g. `metacluster` |
| `n` | int | Event count for that sample x cluster combination |
| `fraction` | dbl | `n` divided by that sample's total event count |
| `<keyword columns>` | chr | Any column that is constant within `file_name` in the input is carried through automatically — marker columns vary per event and fail this test, so only true keyword columns end up passed through, with no extra configuration needed |

Zero-filled complete grid: every `file_name` present in the input gets a
row for every `cluster_col` value present anywhere in the input, with
`n = 0`/`fraction = 0` where that sample had no events in that
cluster/metacluster. This matters for Stage 4's diffcyt GLMM, which needs
true zero counts rather than missing rows to model abundance correctly.

Returned visibly — primary calculated deliverable, not a side effect.

#### Processing steps

1. Error if `cluster_col` is absent from `data`.
2. Compute each sample's total event count from the ungrouped input
   (before any zero-filling).
3. Group by `file_name` and `cluster_col`, count events (`n`).
4. Complete the grid over all `file_name` x `cluster_col` combinations
   seen in the input (`tidyr::complete()`), filling `n = 0`.
5. Compute `fraction = n / <that sample's total from step 2>`.
6. Identify columns constant within `file_name` in the original input
   (excluding `file_name` and `cluster_col` itself) and left-join them
   back in per sample.
7. Return the resulting long tibble.

---

## Error handling

| Condition | Function | Behavior |
|---|---|---|
| `markers` names a column absent from `data` | `facs_cluster_flowsom()` | Error, listing missing marker(s) |
| Explicit `markers` column is not `double`-typed in `data` | `facs_cluster_flowsom()` | Error, listing offending column(s) |
| Resolved `markers` has fewer than 2 columns | `facs_cluster_flowsom()` | Error (guards a `FlowSOM::BuildSOM()` crash when exactly one marker column is selected) |
| `NA` present in a selected marker column | `facs_cluster_flowsom()` | Error, listing affected column(s) |
| `n_metaclusters` is less than 3 | `facs_cluster_flowsom()` | Error (guards a `ConsensusClusterPlus` crash whenever `n_metaclusters < 3`) |
| `n_metaclusters` is not strictly less than grid node count (`grid_xdim * grid_ydim`) | `facs_cluster_flowsom()` | Error (rejects `>=`, not just `>` -- exact equality also crashes `ConsensusClusterPlus`'s `cutree()`) |
| `cluster_col` absent from `data` | `facs_calc_cluster_freq()` | Error |

---

## Testing

Chain off the existing real fixture (`tests/fixtures/Treg.wsp` +
`tests/fixtures/Treg/*.fcs`, gitignored — see `tests/fixtures/README.md`
for regeneration) by feeding `facs_read_fcs_gated()`'s output into
`facs_cluster_flowsom()`. No mocking `FlowSOM`, same no-mock convention as
`fcexpr`/`CytoML`/`flowWorkspace`. Use a small grid (e.g. 3x3) and a low
`n_metaclusters` in tests to keep `R CMD check` runtime reasonable — a
test-tuning choice, not a design constraint.

Test cases:
- Happy path: default `markers`, correct output shape (`cluster`,
  `metacluster` columns appended, expected value ranges).
- `markers` override to a subset.
- Unmatched marker name — error, message lists it.
- Explicit marker naming a non-`double` column (e.g. a `chr` keyword
  column) — error, message lists it.
- Resolved `markers` of fewer than 2 columns (explicit single marker, or
  auto-selection finding only one `dbl` column) — error.
- `NA` in a marker column — error, message lists the column (likely a
  small synthetic tibble to inject the `NA` cleanly rather than the real
  fixture).
- `n_metaclusters` below 3, or at/above grid node count — error in both
  directions.
- Seed reproducibility: same `seed` → identical `cluster`/`metacluster`
  assignments across two calls.
- `facs_calc_cluster_freq()`: zero-fill grid present for at least one
  sample/metacluster combination with no events; `fraction` sums to 1 per
  sample; keyword columns passed through correctly; `cluster_col`
  override (`"cluster"` vs. `"metacluster"`) produces the expected
  column.

---

## Known implementation risk

`FlowSOM::FlowSOM()`'s exact call signature and return structure (how to
extract per-event SOM node vs. consensus metacluster assignments) differs
across package versions and should be verified against the real fixture
during implementation rather than assumed from documentation alone —
same category of risk Stage 1 flagged for `flowjo_to_gatingset()`.

---

## Open items for later stages

- Domain prefix for the diffcyt-based Stage 4 function (`facs_` vs.
  `analysis_`) — still undecided, to be resolved when Stage 4 is
  brainstormed.
- Exact `FlowSOM` package name/version and its transitive dependencies —
  confirmed and added to `Imports` at implementation time.
- Stage 3 (UMAP) is expected to consume this stage's `metacluster` column
  for embedding coloring — no interface changes anticipated, but not
  confirmed until Stage 3 is brainstormed.
