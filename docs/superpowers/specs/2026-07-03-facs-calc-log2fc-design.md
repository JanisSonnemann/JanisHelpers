# Design: `facs_calc_log2fc()` — T-cell restimulation normalization

**Date:** 2026-07-03
**Status:** Approved

---

## Goal

Restimulation panels (e.g. T-cell intracellular cytokine staining) run each mouse/tissue through several conditions — one or more stimuli (e.g. `"MPO"`, `"PMA-Iono"`) plus an unstimulated control (`"unstim"`). Raw or ref-pop-normalized counts vary mouse-to-mouse independent of the stimulus, so the standard readout is a per-mouse fold-change of each stimulated condition over that same mouse's unstimulated control, log2-transformed.

This replaces a manually-written draft (`normalized_log2FC()`) that: hardcoded exactly two stim conditions (`MPO`, `PMA-Iono`) by name, hardcoded `group` as the only passthrough column, returned a wide tibble, and had a pseudocount argument whose default (`0`) contradicted its own comment ("add pseudocount of 0.5 to avoid zeros"). `facs_calc_log2fc()` generalizes all four.

---

## Architecture

New exported function `facs_calc_log2fc()` added to the existing `R/facs_calc.R` (alongside `facs_calc_pct_of()` and `facs_calc_count_per_g()`). No new unexported helpers — the passthrough-column auto-detection is a few lines, not worth extracting.

No new dependencies — uses `dplyr`, `tidyr`, `glue`, `purrr`, already in `Imports`.

Unlike `facs_calc_pct_of()`/`facs_calc_count_per_g()`, this function does **not** append rows to `data` and return it. Its output grain (one row per `mouse_ID` x `tissue` x `population` x non-reference `restimulation` level) collapses the restimulation dimension via comparison, so it doesn't fit the `metric`/`value` long-row-per-file_name shape of the input. It returns a new, differently-shaped tibble instead.

---

## Public API

### `facs_calc_log2fc(data, ref_pop, restim_col = "restimulation", ref_level = "unstim", pseudocount = 0.5)`

| Argument | Type | Default | Notes |
|---|---|---|---|
| `data` | tibble | — | `facs_read_wsp(...)$data`-shaped: must contain `mouse_ID`, `tissue`, `population`, `metric`, `value`, and a `restim_col` column (e.g. joined via `facs_read_wsp(keywords = c("mouse_ID", "tissue", "restimulation"))`) |
| `ref_pop` | chr | — | Leaf population name (matches `population`) used as the denominator for a proportion, e.g. `"CD3+"` |
| `restim_col` | chr | `"restimulation"` | Column in `data` holding the stimulation condition label |
| `ref_level` | chr | `"unstim"` | Value in `restim_col` treated as the baseline; every other distinct value is compared against it |
| `pseudocount` | dbl | `0.5` | Added to every `value` before computing the `ref_pop` proportion, avoiding `log2(0)` when a condition has zero positive events |

Returns a new tibble (not `data` with rows appended).

**Processing steps, in order:**

1. If `restim_col` is not a column in `data`: `stop()`.
2. `ref_counts <- data |> dplyr::filter(population == ref_pop, metric == "count") |> dplyr::select(mouse_ID, tissue, {restim_col}, ref_count = value)`, keyed by `mouse_ID`, `tissue`, `restim_col`.
3. Ambiguity check: if `ref_pop` matches more than one row for any `mouse_ID`/`tissue`/restim-level combo: `stop()` listing the offending combos (mirrors `facs_calc_pct_of()`'s per-`file_name` check, but keyed on the three grouping columns since one `ref_pop` count is expected per condition, not per file).
4. Missing check: if `ref_pop` has zero matches for a combo present in `data`: `warning()` listing the combos; those rows' `proportion` (and downstream `log2fc`) become `NA`.
5. `proportions <- data |> dplyr::filter(metric == "count", population != ref_pop) |> dplyr::left_join(ref_counts, by = c("mouse_ID", "tissue", restim_col)) |> dplyr::mutate(proportion = (value + pseudocount) / ref_count)`.
6. Passthrough-column auto-detection: for every column in `data` other than `mouse_ID`, `tissue`, `restim_col`, `population`, `population_full_path`, `metric`, `value`, `file_name`, keep it if it has exactly one distinct value per `mouse_ID` x `tissue` group (e.g. `group`); drop it otherwise (varies within the group, can't be carried through a pivot keyed on `mouse_ID`/`tissue`/`population`).
7. `tidyr::pivot_wider(id_cols = c(mouse_ID, tissue, population, <passthrough cols>), names_from = {restim_col}, values_from = proportion)`.
8. Reference-level check: if the `ref_level` column is missing entirely from the pivoted result, or is `NA` for some `mouse_ID`/`tissue` group: `warning()` listing the affected combos; `log2fc` becomes `NA` for those rows (via ordinary `NA` propagation through `log2()`).
9. For every pivoted column except `ref_level` and the id columns: compute `log2(<col> / <ref_level column>)`.
10. `tidyr::pivot_longer()` those log2FC columns into a tidy `{restim_col}` / `log2fc` pair (column name reused, e.g. `restimulation`), excluding `ref_level` itself (no self-comparison row).
11. Return the resulting tibble: `mouse_ID`, `tissue`, `population`, `<passthrough cols>`, `{restim_col}`, `log2fc`.

---

## Error handling

| Condition | Behaviour |
|---|---|
| `restim_col` not found in `data` | `stop()` |
| `ref_pop` matches >1 row for a `mouse_ID`/`tissue`/restim-level combo | `stop()`, lists combos |
| `ref_pop` has 0 matches for a `mouse_ID`/`tissue`/restim-level combo present in `data` | `warning()`, lists combos; `NA` fill |
| `ref_level` missing (no row, or `NA` proportion) for a `mouse_ID`/`tissue` group | `warning()`, lists combos; `NA` fill |

---

## Testing

Tests appended to `tests/testthat/test-facs_calc.R`, following the existing `facs_calc_pct_of()`/`facs_calc_count_per_g()` pattern: synthetic tibbles for unit-level logic.

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

No end-to-end fixture test — `tests/fixtures/minimal.wsp` doesn't carry `mouse_ID`/`tissue`/`restimulation` keywords (consistent with `facs_calc_pct_of()`, which is also unit-tested only).

---

## Dependencies

No changes — `dplyr`, `tidyr`, `glue`, `purrr` already in `Imports`.

---

## Roxygen

Full `@param`/`@returns`/`@export` per CLAUDE.md convention. `@examples` uses `\dontrun{}` since it requires realistic multi-column input not easily constructed inline.

---

## CLAUDE.md updates

None needed for the "Data structures" section: `facs_calc_pct_of()` and `facs_calc_count_per_g()` — the two existing `facs_calc_*()` functions — have no entries there either (only roxygen `@returns` documents their shape), so `facs_calc_log2fc()` follows the same precedent. The "Known check output" section does get a new bullet for any `facs_calc_log2fc`-specific variable-binding notes surfaced by `devtools::check()`, matching the existing bullets for `facs_calc_pct_of`/`facs_calc_count_per_g`.
