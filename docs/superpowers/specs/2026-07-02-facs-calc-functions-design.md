# Design: `facs_calc_pct_of()` + `facs_calc_count_per_g()` — derived FACS metrics

**Date:** 2026-07-02
**Status:** Approved

---

## Goal

`facs_read_wsp()` produces raw counts and immediate-parent fractions per population. Two further derived metrics are needed for downstream analysis:

1. **Percentage of an arbitrary ancestor population** (not just the immediate parent) — e.g. every population's fraction of `"CD45+"`, regardless of how many gating levels separate them.
2. **Absolute cell counts per gram of tissue** — converting a raw event count into a physiologically meaningful quantity, using per-mouse organ weights and staining volumes from metadata, with two possible counting methods (bead-based or HTS/volumetric).

These replace the draft in `new_facs_functions.R` (`normalize_to()`, `calc_absolute_counts()` — the latter defined twice, the first definition dead code referencing an undefined `count` variable and a global `meta`). Both drafts also targeted the deprecated `facs_import_wsp()` CamelCase/`"Count"` schema rather than the current `facs_read_wsp()` snake_case/lowercase-metric schema.

---

## Architecture

New file `R/facs_calc.R`, replacing `new_facs_functions.R`, containing two exported functions and one shared unexported helper (`resolve_var_()`, for the column-or-constant argument pattern).

No new dependencies — uses `dplyr`, `tidyr`, `glue`, already in `Imports`.

---

## Public API

### `facs_calc_pct_of(data, ref_pop)`

| Argument | Type | Default | Notes |
|---|---|---|---|
| `data` | tibble | — | `facs_read_wsp(...)$data`-shaped: `file_name`, `population_full_path`, `population`, `metric`, `value` |
| `ref_pop` | chr | — | Leaf population name (matches `population`, not `population_full_path`) to use as the denominator |

Returns `data` with additional rows appended, invisibly.

**Processing steps, in order:**

1. `ref_counts <- data |> dplyr::filter(population == ref_pop, metric == "count")` per `file_name`.
2. Count matches per `file_name`. If any `file_name` has **more than one** match: `stop()` listing the offending `file_name`s and their `population_full_path`s — ambiguous leaf name, caller must resolve upstream (e.g. rename the gate in FlowJo).
3. If any `file_name` has **zero** matches: `warning()` listing those `file_name`s; their new rows get `value = NA_real_`.
4. For every row with `metric == "count"` and `population != ref_pop`: compute `value / ref_count`, as a plain 0–1 fraction (consistent with `fraction_of_parent`, not ×100 despite the `pct` name).
5. New rows get `metric = paste0("pct_of_", ref_pop)`.
6. `dplyr::bind_rows(data, new_rows) |> dplyr::arrange(file_name, population_full_path)`.

### `facs_calc_count_per_g(data, meta, tissue, vol_total, vol_stained, vol_resuspended, vol_measured, organ_piece_weight, method_col = NULL, bead_pop = "beads", bead_concentration = 10400)`

| Argument | Type | Default | Notes |
|---|---|---|---|
| `data` | tibble | — | `facs_read_wsp(...)$data`-shaped, must include a `tissue` column (e.g. joined via `keywords = "tissue"`) |
| `meta` | tibble | — | Per-mouse metadata, keyed by `mouse_ID`, e.g. from `meta_read()` merged with organ-weight/volume columns |
| `tissue` | chr | — | Value to filter `data$tissue` on (e.g. `"kidney"`) — makes the function organ-agnostic |
| `vol_total` | chr or dbl | — | Column name in `meta`, or a single numeric constant: total organ digest volume |
| `vol_stained` | chr or dbl | — | Column name or constant: volume of digest taken for staining |
| `vol_resuspended` | chr or dbl | — | Column name or constant: volume the stained pellet was resuspended in |
| `vol_measured` | chr or dbl | — | Column name or constant: volume actually run/measured on the cytometer |
| `organ_piece_weight` | chr or dbl | — | Column name or constant: weight (mg) of the organ piece processed |
| `method_col` | chr | `NULL` | Column name in `meta` with values `"beads"` or `"hts"` per mouse. `NULL` (default) → every sample uses the HTS formula |
| `bead_pop` | chr | `"beads"` | Population name (leaf) used to look up bead counts |
| `bead_concentration` | dbl | `10400` | Reference bead concentration (beads/µL) used in the bead formula |

Returns `data` with additional rows appended, invisibly.

**Method resolution per mouse (in order):**

1. If `method_col` is `NULL`: every sample uses **HTS**.
2. If `method_col` is given: look up `meta[[method_col]]` per `mouse_ID`.
   - `"hts"` → HTS formula.
   - `"beads"` → bead formula.
   - `NA` → defaults to **HTS** (documented fallback, not a warning — HTS is the standard method).
   - Any other value → `stop()` listing the offending `mouse_ID`s and values.

**Processing steps, in order:**

1. `resolve_var_(meta, x)` helper: if `x` is a length-1 character naming a column in `meta`, return that column; if `x` is a length-1 numeric, recycle it to `nrow(meta)`; else `stop()`.
2. Build `m <- meta` with resolved `vol_total`, `vol_stained`, `vol_resuspended`, `vol_measured`, `organ_piece_weight`, and (if `method_col` supplied) the resolved method per mouse, defaulting `NA` to `"hts"`.
3. Extract bead counts from the **unfiltered** `data`: `dplyr::filter(population == bead_pop, metric == "count")`, keyed by `file_name` (bead counts live in the same per-tissue FCS file as the sample, so `file_name` already disambiguates tissue).
4. Filter `data` to `metric == "count" & tissue == tissue`, left-join bead counts (by `file_name`) and `m` (by `mouse_ID`).
5. If any sample resolves to the bead method but has no matching bead count: `warning()` listing those `file_name`s; `value = NA_real_` for those rows.
6. Compute, row-wise based on resolved method:
   - HTS: `((value / (vol_measured / vol_resuspended)) / (vol_stained / vol_total)) / (organ_piece_weight / 1000)`
   - Beads: `((value / (bead_count / bead_concentration)) / (vol_stained / vol_total)) / (organ_piece_weight / 1000)`
7. New rows get `metric = "count_per_g"`.
8. `dplyr::bind_rows(data, new_rows) |> dplyr::arrange(file_name, population_full_path)`.

---

## Error handling

| Condition | Behaviour |
|---|---|
| `facs_calc_pct_of()`: `ref_pop` matches >1 row for a `file_name` | `stop()`, lists `file_name`s and full paths |
| `facs_calc_pct_of()`: `ref_pop` matches 0 rows for a `file_name` | `warning()`, lists `file_name`s; result rows `NA` |
| `facs_calc_count_per_g()`: `method_col` value outside `{"beads","hts",NA}` | `stop()`, lists `mouse_ID`s and bad values |
| `facs_calc_count_per_g()`: `method_col` value is `NA` | Defaults to `"hts"`, no warning |
| `facs_calc_count_per_g()`: bead method resolved but no bead count found | `warning()`, lists `file_name`s; result rows `NA` |
| `resolve_var_()`: argument is neither a valid column name nor a single numeric | `stop()` |

---

## Testing

New file `tests/testthat/test-facs_calc.R`, following `test-meta_wrangle.R`'s pattern: synthetic tibbles for unit-level logic, plus one end-to-end test against the real fixtures.

**`facs_calc_pct_of()`:**
- Correctly computes `value / ref_count` as a 0–1 fraction for a simple two-population case
- Errors when `ref_pop` matches more than one row for a `file_name`
- Warns and fills `NA` when `ref_pop` has no match for a `file_name`
- Excludes `ref_pop`'s own row from the output rows added

**`facs_calc_count_per_g()`:**
- HTS formula produces the expected value for a simple case (`method_col = NULL`)
- Bead formula produces the expected value when `method_col` selects `"beads"`
- `NA` in `method_col` falls back to HTS with no warning
- Errors on an invalid `method_col` value
- Warns and fills `NA` when bead method resolved but no bead count present
- `tissue` argument correctly restricts which rows are processed
- `vol_total`/`vol_stained`/etc. accept both a column name and a numeric constant

**End-to-end:**
- `facs_read_wsp(minimal.wsp, keywords = c("mouse_ID", "tissue"))$data` piped through both functions using `meta_minimal.xlsx` (extended with the organ-weight/volume columns needed, or a synthetic `meta` merged with real `mouse_ID`s from the fixture).

---

## Dependencies

No changes — `dplyr`, `tidyr`, `glue` already in `Imports`.

---

## Roxygen

Both functions get full `@param`/`@returns`/`@export` per CLAUDE.md convention. `@examples` use `\dontrun{}` since both require realistic multi-column input not easily constructed inline.

---

## Naming / conventions housekeeping

- `new_facs_functions.R` is deleted; replaced by `R/facs_calc.R` (matches `domain_verb.R` file convention).
- All calls namespaced (`dplyr::`, `tidyr::`, `glue::`) per CLAUDE.md.
- CLAUDE.md's `facs_read_wsp()` data-structure table is unaffected (no changes to that function); no CLAUDE.md updates needed for this work beyond it already being current.
