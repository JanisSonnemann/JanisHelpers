# Design: multi-sheet meta import + organ-weight/facs-volume wrangling for `facs_calc_count_per_g()`

**Date:** 2026-07-02
**Status:** Approved by user, ready for implementation plan

## Problem

`meta_minimal.xlsx` now has three sheets: `meta` (per-mouse subject metadata),
`organ_weights` (one row per `mouse_ID`, wide across tissue — e.g.
`kidney_total_weight`, `kidney_facs_weight`, `lung_total_weight`, ...), and
`facs_volumes` (one row per `mouse_ID` x `tissue`, wide across staining panel
— e.g. `overview_vol`, `overview_resuspended_vol`, `overview_measured_vol`,
`Treg_vol`, ...).

`facs_calc_count_per_g()` needs a `meta` tibble with exactly one row per
`mouse_ID` (today) containing `vol_total`, `vol_stained`, `vol_resuspended`,
`vol_measured`, `organ_piece_weight`-equivalent columns. Its current
single-`tissue` argument requires the caller to pre-filter `meta` to one
tissue per call, which is fragile (silent row fan-out on the join if
forgotten) and doesn't fit `organ_weights`/`facs_volumes`'s natural shape.

This design changes `facs_calc_count_per_g()` to process every tissue
present in `data` in a single call (joining `meta` on `mouse_ID` + `tissue`),
and adds the import/wrangling functions needed to turn the raw
`organ_weights`/`facs_volumes` sheets into that shape.

## Changes

### 1. `meta_read(path)` — reads all sheets

- Drops the `sheet` parameter.
- Reads every sheet via `readxl::excel_sheets(path)`.
- Applies the existing per-sheet cleaning independently to each sheet:
  `janitor::remove_empty()`, `janitor::clean_names()`, `mouse_id` ->
  `mouse_ID` rename (only if present), character-column trimming, date/time
  column coercion to `Date`, `group` -> factor coercion (only if present).
- Returns (invisibly) a named list of tibbles, one per sheet, named after
  the sheet (e.g. `list(meta = ..., organ_weights = ..., facs_volumes = ...)`).

### 2. `meta_clean(meta_list)` — new function

- New exported function, domain `meta_`, added to `R/meta_wrangle.R`.
- Input: the named list returned by `meta_read()`. Errors clearly if
  `organ_weights` or `facs_volumes` elements are missing.
- Behavior:
  1. Pivots `meta_list$organ_weights` from wide (tissue-prefixed columns)
     to long, via an unexported helper. Column-name pattern:
     `^(.*)_(total_weight|facs_weight)$` (tissue-prefix capture group +
     hardcoded measure-suffix capture group, via
     `tidyr::pivot_longer(cols = !mouse_ID, names_to = c("tissue", ".value"), names_pattern = ...)`).
     Suffixes are hardcoded (`total_weight`, `facs_weight`) rather than
     generic, matching only the current `organ_weights` sheet shape —
     acceptable since this function is purpose-built for that sheet.
     Result: `mouse_ID`, `tissue`, `total_weight`, `facs_weight`.
  2. Leaves `meta_list$facs_volumes` as-is (already long by `tissue`,
     already cleaned by `meta_read()`).
  3. Joins the two via `meta_annotate(organ_weights_long, facs_volumes, by = c("mouse_ID", "tissue"))`,
     reusing `meta_annotate()`'s collision-check and unmatched-combo warning
     (see change 3 below).
- Returns the combined tibble: one row per `mouse_ID` x `tissue`, with
  `total_weight`, `facs_weight`, `total_vol`, `overview_vol`,
  `overview_resuspended_vol`, `overview_measured_vol`, `Treg_vol`,
  `Treg_resuspended_vol`, `Treg_measured_vol`. Ready to pass as `meta` to
  `facs_calc_count_per_g()` (selecting the relevant panel's volume columns
  by name at the call site).

### 3. `meta_annotate()` — `by` accepts a character vector

- `by` parameter generalized from a single string to a character vector,
  default unchanged (`"mouse_ID"`).
- Missing-column check: loop over each element of `by`, checking presence
  in both `data` and `meta` (error naming which side is missing which
  column).
- Collision check unchanged — `setdiff(intersect(names(data), names(meta)), by)`
  already works correctly with a vector `by`.
- Unmatched-value check: generalized from comparing single-column unique
  values to comparing distinct combinations of the `by` columns (e.g. via
  anti-join of distinct `by`-column tuples), warning with the unmatched
  combinations listed (e.g. `"mouse_ID=m1, tissue=lung"`).

### 4. `facs_calc_count_per_g()` — remove `tissue` argument

- `tissue` parameter removed from the signature entirely.
- The join between `data` and resolved `meta` (`m`) changes from
  `by = "mouse_ID"` to `by = c("mouse_ID", "tissue")`. `meta` must now
  carry a `tissue` column.
- No more single-tissue filter on `data` — every `mouse_ID` x `tissue`
  combination present in `data` (after filtering `metric == "count"`,
  excluding `bead_pop`) is processed in one call.
- Unmatched-mouse warning generalizes to unmatched `(mouse_ID, tissue)`
  pairs (`data` combos with no match in `meta`) — same warn + NA-fill
  behavior as today, keyed on the pair instead of `mouse_ID` alone.
- All other behavior (bead vs. HTS formula, `method_col` resolution order,
  `resolve_var_()` column-or-constant handling) is unchanged.

## End-to-end example

```r
meta_list     <- meta_read("meta.xlsx")
meta_combined <- meta_clean(meta_list)

facs_data <- facs_read_wsp("experiment.wsp", keywords = c("mouse_ID", "tissue"))$data

result <- facs_calc_count_per_g(
  facs_data, meta_combined,
  vol_total = "total_vol", vol_stained = "overview_vol",
  vol_resuspended = "overview_resuspended_vol", vol_measured = "overview_measured_vol",
  organ_piece_weight = "facs_weight"
)
```

One call computes `count_per_g` across every tissue present in `facs_data`
for the "overview" panel. A second call with `vol_stained = "Treg_vol"`,
`vol_resuspended = "Treg_resuspended_vol"`, `vol_measured = "Treg_measured_vol"`
covers the "Treg" panel.

Separately, the `meta` sheet (subject metadata: sex, group, dates, ...)
still joins onto FACS data by `mouse_ID` alone, e.g.
`meta_annotate(facs_data, meta_list$meta)` — unrelated to this workflow.

## Breaking changes / test impact

- `meta_read()`'s return type changes from a single tibble to a named list
  of tibbles. All 9 existing tests in `tests/testthat/test-meta_wrangle.R`
  that call `meta_read()` need rewriting for the new return shape.
- `meta_annotate()`'s unmatched-warning message format changes for the
  default single-column case only if the combination-based implementation
  changes wording; existing single-`by` tests must still pass unchanged in
  behavior (warn + NA-fill on unmatched values), only the internal
  implementation generalizes.
- `facs_calc_count_per_g()` signature changes (`tissue` removed). All
  `facs_calc_count_per_g()`-related tests in `tests/testthat/test-facs_calc.R`
  (~13 tests, including the end-to-end fixture test) need their `meta`
  fixtures updated to include a `tissue` column and their calls updated to
  drop `tissue = ...`. At least one new test should cover processing
  multiple tissues in a single call.
- `CLAUDE.md`'s "Data structures" section needs updating: `meta_read()`
  return shape, new `meta_clean()` entry, `meta_annotate()`'s `by` now a
  vector, `facs_calc_count_per_g()`'s signature (no `tissue` arg, `meta`
  now requires `tissue` column).
- No new package dependencies — `tidyr`, `janitor`, `readxl` already in
  `Imports`.

## Out of scope

- No change to `facs_calc_pct_of()`.
- No change to the bead-vs-HTS formula math itself.
- No generic/reusable "wide-by-tissue" pivot utility — the pivot helper is
  purpose-built for `organ_weights`'s exact column shape (hardcoded
  `total_weight`/`facs_weight` suffixes), per YAGNI.
- No handling of tissue names containing underscores beyond what the
  hardcoded-suffix regex already tolerates (tissue prefix is captured via
  `(.*)`, not `[^_]+`, so multi-word tissue names are safe).
