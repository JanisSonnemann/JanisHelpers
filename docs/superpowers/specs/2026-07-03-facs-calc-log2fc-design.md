# Design: `facs_calc_log2fc()` + `facs_calc_diff()` — T-cell restimulation normalization

**Date:** 2026-07-03
**Status:** Approved

---

## Goal

Restimulation panels (e.g. T-cell intracellular cytokine staining) run each mouse/tissue through several conditions — one or more stimuli (e.g. `"MPO"`, `"PMA-Iono"`) plus an unstimulated control (`"unstim"`). Raw or ref-pop-normalized counts vary mouse-to-mouse independent of the stimulus, so the standard readout compares each stimulated condition against that same mouse's unstimulated control. Two comparison forms are needed: a log2 fold-change (`facs_calc_log2fc()`) and a plain difference (`facs_calc_diff()`).

This replaces a manually-written draft (`normalized_log2FC()`) that: hardcoded exactly two stim conditions (`MPO`, `PMA-Iono`) by name, hardcoded `group` as the only passthrough column, returned a wide tibble, and had a pseudocount argument whose default (`0`) contradicted its own comment ("add pseudocount of 0.5 to avoid zeros"). Both new functions generalize all four issues; `facs_calc_diff()` is new, not present in the draft.

---

## Architecture

Two new exported functions, `facs_calc_log2fc()` and `facs_calc_diff()`, added to the existing `R/facs_calc.R` (alongside `facs_calc_pct_of()` and `facs_calc_count_per_g()`). They share all logic through the ref-pop-normalized-proportion step — computing `ref_counts`, the ambiguity/missing-`ref_pop` checks, the proportion itself, passthrough-column auto-detection, the wide pivot, and the missing-`ref_level` check — via one new unexported helper, `calc_restim_proportions_()`. Each public function then applies its own final comparison (`log2(x / ref)` vs. `x - ref`) and its own `pivot_longer()` (with its own value-column name: `log2fc` vs. `diff`). This mirrors the existing `resolve_var_()` pattern shared by `facs_calc_count_per_g()`.

No new dependencies — uses `dplyr`, `tidyr`, `glue`, `purrr`, already in `Imports`.

Unlike `facs_calc_pct_of()`/`facs_calc_count_per_g()`, neither function appends rows to `data` and returns it. Their output grain (one row per `mouse_ID` x `tissue` x `population` x non-reference restimulation level) collapses the restimulation dimension via comparison, so it doesn't fit the `metric`/`value` long-row-per-`file_name` shape of the input. Each returns a new, differently-shaped tibble instead.

---

## Internal helper

### `calc_restim_proportions_(data, ref_pop, restim_col, ref_level, pseudocount)`

Not exported. Returns a wide tibble: one row per `mouse_ID` x `tissue` x `population`, columns `mouse_ID`, `tissue`, `population`, `<passthrough cols>`, plus one column per distinct value of `restim_col` (holding that condition's `ref_pop`-normalized proportion). `ref_level`'s column is guaranteed present (filled `NA` if the level never occurs in `data`).

**Processing steps, in order:**

1. If `restim_col` is not a column in `data`: `stop()`.
2. `ref_counts <- data |> dplyr::filter(population == ref_pop, metric == "count") |> dplyr::select(mouse_ID, tissue, {restim_col}, ref_count = value)`, keyed by `mouse_ID`, `tissue`, `restim_col`.
3. Ambiguity check: if `ref_pop` matches more than one row for any `mouse_ID`/`tissue`/restim-level combo: `stop()` listing the offending combos.
4. Missing check: if `ref_pop` has zero matches for a combo present in `data`: `warning()` listing the combos; those rows' `proportion` becomes `NA`.
5. `proportions <- data |> dplyr::filter(metric == "count", population != ref_pop) |> dplyr::left_join(ref_counts, by = c("mouse_ID", "tissue", restim_col)) |> dplyr::mutate(proportion = (value + pseudocount) / ref_count)`.
6. Passthrough-column auto-detection: for every column in `data` other than `mouse_ID`, `tissue`, `restim_col`, `population`, `population_full_path`, `metric`, `value`, `file_name`, keep it if it has exactly one distinct value per `mouse_ID` x `tissue` group (e.g. `group`); drop it otherwise.
7. `tidyr::pivot_wider(id_cols = c(mouse_ID, tissue, population, <passthrough cols>), names_from = {restim_col}, values_from = proportion)`.
8. Reference-level check: if the `ref_level` column is missing entirely from the pivoted result (add it, filled `NA`), or is `NA` for some `mouse_ID`/`tissue` group: `warning()` listing the affected combos.
9. Return the wide tibble.

---

## Public API

### `facs_calc_log2fc(data, ref_pop, restim_col = "restimulation", ref_level = "unstim", pseudocount = 0.5)`

| Argument | Type | Default | Notes |
|---|---|---|---|
| `data` | tibble | — | `facs_read_wsp(...)$data`-shaped: must contain `mouse_ID`, `tissue`, `population`, `metric`, `value`, and a `restim_col` column |
| `ref_pop` | chr | — | Leaf population name (matches `population`) used as the denominator for a proportion, e.g. `"CD3+"` |
| `restim_col` | chr | `"restimulation"` | Column in `data` holding the stimulation condition label |
| `ref_level` | chr | `"unstim"` | Value in `restim_col` treated as the baseline; every other distinct value is compared against it |
| `pseudocount` | dbl | `0.5` | Added to every `value` before computing the `ref_pop` proportion, avoiding `log2(0)` when a condition has zero positive events |

**Processing:** calls `calc_restim_proportions_(data, ref_pop, restim_col, ref_level, pseudocount)`, then for every column except `ref_level` and the id columns computes `log2(<col> / <ref_level column>)`, then `tidyr::pivot_longer()` into a tidy `{restim_col}` / `log2fc` pair (excluding `ref_level` itself).

Returns a new tibble: `mouse_ID`, `tissue`, `population`, `<passthrough cols>`, `{restim_col}`, `log2fc`.

### `facs_calc_diff(data, ref_pop, restim_col = "restimulation", ref_level = "unstim", pseudocount = 0)`

Identical signature and shared processing to `facs_calc_log2fc()`, except:
- `pseudocount` defaults to `0` (subtraction has no divide-by-zero/log-of-zero failure mode, so no pseudocount is needed by default; still available for callers who want proportions comparable to a `facs_calc_log2fc()` call using a nonzero `pseudocount`).
- Final comparison is `<col> - <ref_level column>` instead of `log2(<col> / <ref_level column>)`.
- Output value column is named `diff` instead of `log2fc`.

Returns a new tibble: `mouse_ID`, `tissue`, `population`, `<passthrough cols>`, `{restim_col}`, `diff`.

---

## Error handling

Both functions share identical error/warning behavior, since it all happens inside `calc_restim_proportions_()`:

| Condition | Behaviour |
|---|---|
| `restim_col` not found in `data` | `stop()` |
| `ref_pop` matches >1 row for a `mouse_ID`/`tissue`/restim-level combo | `stop()`, lists combos |
| `ref_pop` has 0 matches for a `mouse_ID`/`tissue`/restim-level combo present in `data` | `warning()`, lists combos; `NA` fill |
| `ref_level` missing (no row, or `NA` proportion) for a `mouse_ID`/`tissue` group | `warning()`, lists combos; `NA` fill |

---

## Testing

Tests appended to `tests/testthat/test-facs_calc.R`, following the existing `facs_calc_pct_of()`/`facs_calc_count_per_g()` pattern: synthetic tibbles for unit-level logic. `calc_restim_proportions_()` is not tested directly — it's exercised indirectly through both public functions, matching the existing `resolve_var_()` precedent (also untested directly).

**`facs_calc_log2fc()`:**
- Correctly computes `log2((value + pseudocount) / ref_count)` fold-change for a simple two-condition case (one stim vs `"unstim"`)
- Handles more than one non-reference `restim_col` level in a single call (e.g. `"MPO"` and `"PMA-Iono"` both present), producing one row per level
- Excludes `ref_level` from the output rows (no `"unstim"` row in `restim_col`)
- Excludes `ref_pop`'s own row from consideration (not present as a `population` in the output)
- `pseudocount` default (`0.5`) avoids `-Inf`/`NaN` when a condition's `value` is `0`
- Auto-detects and carries through a constant passthrough column (e.g. `group`)
- Drops a passthrough-candidate column that varies within a `mouse_ID`/`tissue` group (does not error, just excluded)
- Errors when `ref_pop` matches more than one row for a `mouse_ID`/`tissue`/restim-level combo
- Warns and fills `NA` when `ref_pop` has no match for a combo
- Warns and fills `NA` when `ref_level` (`"unstim"`) is missing for a `mouse_ID`/`tissue` group
- Errors when `restim_col` is not a column in `data`
- `restim_col`/`ref_level` arguments work with non-default names/values (not just `"restimulation"`/`"unstim"`)

**`facs_calc_diff()`:**
- Correctly computes `(value + pseudocount) / ref_count` difference (`proportion_stim - proportion_ref`) for a simple two-condition case
- `pseudocount` defaults to `0` (not `0.5`) — verify the unadjusted proportion is used unless the caller overrides it
- Excludes `ref_level`/`ref_pop` from the output, same as `facs_calc_log2fc()`
- Errors when `ref_pop` matches more than one row for a combo (same shared-helper behavior as `facs_calc_log2fc()` — one test suffices to confirm the helper is wired up; full branch coverage already lives on the `facs_calc_log2fc()` tests)
- Warns and fills `NA` when `ref_level` is missing for a group

No end-to-end fixture test for either function — `tests/fixtures/minimal.wsp` doesn't carry `mouse_ID`/`tissue`/`restimulation` keywords (consistent with `facs_calc_pct_of()`, which is also unit-tested only).

---

## Dependencies

No changes — `dplyr`, `tidyr`, `glue`, `purrr` already in `Imports`.

---

## Roxygen

Both exported functions get full `@param`/`@returns`/`@export` per CLAUDE.md convention. `@examples` uses `\dontrun{}` since both require realistic multi-column input not easily constructed inline. `calc_restim_proportions_()` is unexported — no roxygen block required (matches `resolve_var_()` precedent), but gets a one-line `#` comment above it noting it's shared by both public functions.

---

## CLAUDE.md updates

None needed for the "Data structures" section: `facs_calc_pct_of()` and `facs_calc_count_per_g()` — the two existing `facs_calc_*()` functions — have no entries there either (only roxygen `@returns` documents their shape), so `facs_calc_log2fc()`/`facs_calc_diff()` follow the same precedent. The "Known check output" section does get a new bullet for any `facs_calc_log2fc`/`facs_calc_diff`/`calc_restim_proportions_`-specific variable-binding notes surfaced by `devtools::check()`, matching the existing bullets for `facs_calc_pct_of`/`facs_calc_count_per_g`.
