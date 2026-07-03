# Design: `facs_read_fcs_gated()` — raw single-cell event import with FlowJo gate replay

**Date:** 2026-07-03
**Status:** Approved

This is Stage 1 of the unsupervised FACS analysis pipeline described in
`docs/superpowers/specs/2026-07-03-facs-unsupervised-analysis-roadmap.md`.
It is the foundation every later stage (FlowSOM clustering, UMAP,
diffcyt differential abundance) reads from.

---

## Goal

Every existing `facs_` function operates on FlowJo's already-gated summary
statistics (`facs_read_wsp()`'s long tibble of counts/fractions/stats).
Unsupervised, cell-level analysis needs the raw single-cell events
themselves. `facs_read_fcs_gated()` reads those raw events straight from
`.fcs` files, filtered down to an arbitrary gate already drawn in the
matching `.wsp` workspace, for a caller-specified marker panel — reusing
the existing gating investment (debris/dead cells/doublets already
excluded) instead of starting from ungated events.

---

## Architecture

New file `R/facs_read_fcs.R`, containing one exported function,
`facs_read_fcs_gated()`.

New dependencies: `CytoML` and `flowWorkspace` (Bioconductor, non-CRAN —
same treatment as `fcexpr`: listed in `Imports`, install requirement
documented in the README, no `remotes::`/`BiocManager::` calls inside
function bodies). `CytoML::open_flowjo_xml()` +
`flowWorkspace::flowjo_to_gatingset()` replay the workspace's
compensation, transformation, and gating tree against the raw `.fcs`
files; `flowWorkspace::gh_pop_get_data()` (or equivalent) extracts the
gated event matrix per sample.

---

## Public API

```r
facs_read_fcs_gated(
  wsp_path,
  gate_path,
  markers,
  keywords = NULL,
  fcs_dir = NULL,
  group = "All Samples",
  max_events = NULL,
  seed = NULL
)
```

| Argument | Type | Default | Notes |
|---|---|---|---|
| `wsp_path` | chr | — | Path to the `.wsp` file |
| `gate_path` | chr | — | Full gating path, e.g. `"/Live/Singlets/CD45+"` — same format as `PopulationFullPath` from `facs_read_wsp()`. Applied to every sample in the workspace. |
| `markers` | chr vector | — | Matched per sample against stain label or channel name, `coalesce(stain, channel)` — same priority as `facs_read_wsp()` |
| `keywords` | chr vector | `NULL` | FlowJo keyword names appended as columns, same convention as `facs_read_wsp()` |
| `fcs_dir` | chr | `NULL` | Folder to search for this workspace's `.fcs` files. Default `NULL` auto-derives it as the subfolder named after the `.wsp` file (sans extension), sitting next to it — matches the `/FACS/<panel>.wsp` + `/FACS/<panel>/*.fcs` layout. Pass an explicit path to override (e.g. if files were reorganized after acquisition). |
| `group` | chr | `"All Samples"` | FlowJo sample group to load, passed through to `flowWorkspace::flowjo_to_gatingset()` |
| `max_events` | int | `NULL` | If set, randomly downsample each sample to at most this many events |
| `seed` | int | `NULL` | If set, seeds the random draw for reproducible downsampling |

### Output

Wide tibble, one row per event:

| Column | Type | Notes |
|---|---|---|
| `file_name` | chr | FCS filename, same convention as `facs_read_wsp()` |
| `<marker>` | dbl | One column per requested marker, on FlowJo's transformed scale (e.g. logicle/biexponential) — the scale the gates were drawn on, and what FlowSOM/UMAP expect |
| `<keyword>` | chr | One column per requested keyword, appended on the right, same convention as `facs_read_wsp()` |

Returned visibly (not `invisible()`) — this function's output is the
primary deliverable, not a side effect of an import+message step.

---

## Processing steps

1. Resolve `fcs_dir` (auto-derived default or override).
2. `CytoML::open_flowjo_xml(wsp_path)`, then
   `flowWorkspace::flowjo_to_gatingset(ws, name = group, path = fcs_dir)`
   to build a `GatingSet` with compensation, transformation, and gating
   tree replayed.
3. For each sample (`GatingHierarchy`) in the `GatingSet`:
   a. Check `gate_path` exists in that sample's gating tree — if not,
      record it as skipped (see error handling).
   b. Extract the gated event matrix at `gate_path` (already compensated
      + transformed).
   c. Match `markers` against that sample's panel (stain label / channel
      name); subset to the matched columns.
   d. If `max_events` is set and the sample exceeds it, subsample down to
      `max_events` (using `seed` if provided).
   e. Convert to a tibble, add `file_name`, append keyword columns (warn +
      `NA`-fill for any missing, matching `facs_read_wsp()`).
4. `dplyr::bind_rows()` every sample's tibble and return the combined
   result.

---

## Error handling

| Condition | Behavior |
|---|---|
| `gate_path` missing for a sample | Warn (listing affected `file_name`s), skip that sample, continue processing the rest |
| A requested marker doesn't match any channel/stain in a sample's panel | Error immediately, listing the unmatched marker(s) and affected sample |
| A requested keyword is missing for every sample | Warn (matching `facs_read_wsp()`'s whole-workspace check), fill `NA`. A keyword missing for only some samples is filled `NA` for those without a warning — legitimate per-sample variation, not an error condition. |
| An `.fcs` file referenced by the workspace can't be found under `fcs_dir` | Error — indicates a broken file layout, not legitimate per-sample variation, so it does not get warn-and-skip treatment |

---

## Testing

New fixture needed in `tests/fixtures/`, following the exact convention
already established for `minimal.wsp` (`tests/fixtures/README.md`): a
small real pilot experiment, not a synthetic/mocked one (no mocking
`CytoML`/`flowWorkspace`, same rule as the existing no-mock-`fcexpr`
rule). Commit it to git only if it contains no patient/animal-identifiable
metadata; otherwise add it to `.gitignore` (as `minimal.wsp` and
`meta_minimal.xlsx` already are) and document regeneration steps in the
README instead.

Unlike `facs_read_wsp()`'s fixture (which only needs the `.wsp`, since it
reads exported stats from the XML), this stage needs real `.fcs` files
too, placed in a subfolder named after the workspace — e.g.
`tests/fixtures/minimal_fcs/<subfolder-matching-wsp-name>/*.fcs` —
mirroring the real `/FACS/<panel>.wsp` + `/FACS/<panel>/*.fcs` layout so
the `fcs_dir` auto-derivation logic has something real to resolve
against.

Fixture requirements (finalized at plan time): at least 2 `.fcs` files,
a gating hierarchy at least 2 levels deep so a non-trivial `gate_path`
can be exercised, a panel with at least one channel carrying a stain
label (`$PnS`) so label-vs-channel matching can be tested, and at least
one custom keyword.

Test cases to cover: happy path (valid gate + markers → correct
row/column shape and transformed-scale values), missing gate on one
sample (warn + skip, others still returned), unmatched marker (error,
message lists the marker), missing keyword (warn + `NA`-fill),
`max_events` downsampling (row count capped, reproducible with `seed`),
`fcs_dir` auto-derivation vs. explicit override, and non-default `group`.

---

## Known implementation risk

`flowjo_to_gatingset()`'s exact behavior around sample groups, and
whether `path`-based recursive `.fcs` lookup behaves as documented when
pointed at a per-panel subfolder, should be verified against the real
fixture during implementation rather than assumed from documentation
alone.

---

## Open items for later stages

- Domain prefix for the diffcyt-based Stage 4 function (`facs_` vs.
  `analysis_`) — decided when Stage 4 is brainstormed.
- Exact package names/versions for `FlowSOM`, `uwot`/`umap`, `diffcyt`
  and their transitive dependencies — confirmed and added to `Imports`
  at each stage's implementation time.
