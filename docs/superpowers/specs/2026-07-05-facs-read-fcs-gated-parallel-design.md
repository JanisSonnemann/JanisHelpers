# Design: Parallelize `facs_read_fcs_gated()` across samples

**Date:** 2026-07-05
**Status:** Approved

---

## Goal

`facs_read_fcs_gated()` (added in
`docs/superpowers/specs/2026-07-03-facs-read-fcs-gated-design.md`) takes
~3 minutes to import a ~5GB raw FACS dataset with no down-sampling, on a
12-core machine, single-threaded. The user wants this as fast as
possible with no specific memory constraint (nothing crashes today).
This design parallelizes the per-sample work across CPU cores.

---

## Evidence gathered before designing

Profiled against `tests/fixtures/Treg.wsp` + `tests/fixtures/Treg/*.fcs`
(6 samples):

- `CytoML::flowjo_to_gatingset()` (compensation + transformation +
  gating tree replay) accounts for ~75% of total runtime and scales
  **linearly per sample** (2 samples: 0.100s elapsed; 6 samples: 0.265s
  elapsed) — confirming it is effectively an independent per-file loop
  internally, not batched/shared work.
- The post-build per-sample extraction loop (`gh_pop_get_data()` +
  `realize_view()` + `exprs()`) is only ~25% of total runtime.
- Building a single-sample `GatingSet` via
  `CytoML::flowjo_to_gatingset(ws, name = group, path = fcs_dir, subset = sample_name)`
  produces a **byte-identical** gated event matrix (verified with
  `all.equal()`) to building that same sample inside the full
  multi-sample group — gating is purely per-sample, no cross-sample
  state, so per-sample GatingSet construction is safe to parallelize.
- `subset` accepts the sample's clean file name (the same string used
  as `pData(gh)$name` / the `file_name` column), not just a numeric
  index.
- `CytoML::fj_ws_get_sample_groups(ws)` (groupName -> sampleID) joined
  with `CytoML::fj_ws_get_samples(ws)` (sampleID -> name) resolves the
  list of sample names belonging to a `group` directly from the
  workspace XML metadata — **no GatingSet build required** just to
  enumerate samples, so this can happen once in the parent process
  before dispatching parallel workers.
- `parallel::mclapply(..., mc.set.seed = TRUE)` (the default) gives
  independent, non-identical random draws across forked workers —
  confirmed empirically.

Because the dominant cost (~75%) is linear-per-sample and independent
across samples, parallelizing the *entire* per-sample pipeline (not
just the smaller post-build extraction step) is the approach with real
payoff.

---

## Architecture

Modify `R/facs_read_fcs.R`. No new files, no new exported functions —
`facs_read_fcs_gated()` keeps its existing signature plus one new
argument.

**New dependency: none.** `parallel` is a base R package (ships with
every R install) — no `Imports` change needed.

### Changed public API

```r
facs_read_fcs_gated(
  wsp_path,
  gate_path,
  markers,
  keywords = NULL,
  fcs_dir = NULL,
  group = "All Samples",
  max_events = NULL,
  seed = NULL,
  workers = 1L        # new
)
```

| Argument | Type | Default | Notes |
|---|---|---|---|
| `workers` | int | `1L` | Number of samples to process in parallel via `parallel::mclapply()`. `1` (default) preserves today's exact sequential behavior — existing tests/output are unaffected. `> 1` forks one worker per in-flight sample, each building its own single-sample `GatingSet`. Fork-based (`parallel::mclapply()`), so on Windows this silently runs sequentially regardless of the value passed — documented, no platform-detection code needed. |

Output shape, `keywords` handling, and all other arguments are
unchanged from the existing design.

### Processing steps (revised)

1. Resolve `fcs_dir` (unchanged).
2. `ws <- CytoML::open_flowjo_xml(wsp_path)` once, in the parent process
   (cheap: ~0.04s observed).
3. Resolve the sample name list for `group` via
   `CytoML::fj_ws_get_sample_groups(ws)` joined to
   `CytoML::fj_ws_get_samples(ws)` on `sampleID` — no GatingSet build.
4. For each sample name, dispatch a self-contained worker
   (`read_one_sample_()`, refactored) that:
   a. Builds its own single-sample `GatingSet`:
      `CytoML::flowjo_to_gatingset(ws, name = group, path = fcs_dir, subset = sample_name)`.
   b. Looks up `gate_path` in that sample's gating tree (unchanged
      warn-and-skip-on-missing behavior).
   c. Extracts the gated matrix, matches `markers` (unchanged).
   d. If `max_events` is set, downsamples using a **per-sample derived
      seed** (see below), independent of execution order.
   e. Converts to a tibble with `file_name` (unchanged).
   - When `workers == 1`: dispatched via `lapply()` (today's behavior,
     same as the current `purrr::map()`, just avoiding fork overhead
     entirely).
   - When `workers > 1`: dispatched via
     `parallel::mclapply(sample_names, ..., mc.cores = workers, mc.set.seed = TRUE)`.
5. After dispatch, scan results for `inherits(x, "try-error")` (only
   possible when `workers > 1`, since `mclapply` captures per-worker
   errors instead of propagating them immediately). If any exist,
   `stop()` re-raising the first captured error's message — preserves
   today's "unmatched marker halts the whole read" semantics, just
   delayed until every in-flight worker finishes rather than instant.
6. `dplyr::bind_rows()` every sample's tibble (unchanged).
7. Keyword lookup/join (unchanged — already happens once, top-level,
   reading the `.wsp` XML directly).

---

## Reproducible downsampling under parallelism

Today: `set.seed(seed)` once before the sequential loop; each sample's
`sample.int()` call consumes the next slice of one shared RNG stream.
This is already order-dependent (inserting a sample earlier in the list
shifts every later sample's draw) even in the sequential version.

New: derive a per-sample seed as `seed + match(sample_name, sample_names)`
(i.e. `seed` offset by that sample's position in the resolved sample
list) and call `set.seed()` with it *inside* the worker, before
`sample.int()`. This makes each sample's downsampling draw depend only
on its own identity/position, not on fork scheduling or how many other
samples ran before it — strictly more robust than the current
behavior, not merely a parallel-safety shim.

**Caveat to document in `@param seed`**: reproducible for a fixed
`sample_names` ordering (which is stable for a given `wsp_path` +
`group`), but the exact per-sample seed is an internal implementation
detail, not a public contract — do not rely on the specific values
downsampled events, only on repeatability across repeated calls with
the same inputs.

---

## Error handling (updated)

| Condition | Behavior |
|---|---|
| `gate_path` missing for a sample | Unchanged: warn (listing affected `file_name`s), skip that sample |
| A requested marker doesn't match any channel/stain in a sample's panel | Unchanged in effect (error, listing unmatched marker + sample) but now surfaces after all in-flight parallel workers complete rather than the instant one sample fails |
| A requested keyword is missing for every sample | Unchanged |
| An `.fcs` file referenced by the workspace can't be found under `fcs_dir` | Unchanged |
| `workers > 1` on Windows | No error — `parallel::mclapply()` runs sequentially internally; document as a known limitation, not a runtime warning (avoids noisy output on every call) |

---

## Testing

- All existing tests keep the default `workers = 1L` — output must
  remain byte-identical to the current implementation (this is a
  regression bar, not just a nice-to-have, since the internal
  restructuring from "one shared multi-sample GatingSet" to "one
  GatingSet per sample" must not change results).
- New test: run with `workers = 2` against `tests/fixtures/Treg.wsp` and
  assert the result matches the `workers = 1` run (compare as sets of
  rows — `dplyr::arrange()` both before comparing — since parallel
  scheduling may return samples in a different order than input order).
- New test: unmatched-marker error still surfaces (as an error, not
  silently swallowed) when `workers = 2`.
- No new test needed for the Windows sequential-fallback behavior — not
  practically testable in this project's CI/dev environment (macOS).

---

## Explicitly out of scope

Ruled out during brainstorming because the user confirmed there is no
current memory pressure (nothing crashes) and no specific target scale
beyond "as fast as possible":

- Disk-backed / on-demand (`h5`) cytoframe backends for the `GatingSet`.
- Streaming or incrementally writing per-sample output to disk instead
  of collecting one in-memory tibble.
- Progress bar / progress reporting during the parallel dispatch.
- Auto-detecting a default `workers` count (e.g. `parallel::detectCores() - 1`)
  — kept explicit and conservative (`1L` default) since this is a
  shared personal package, not a single-purpose script.

If dataset sizes grow enough to reintroduce memory pressure, revisit
disk-backed backends and/or streaming output as a separate design.
