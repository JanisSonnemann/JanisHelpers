# Design: `meta_read()` + `meta_annotate()` — experiment metadata import and annotation

**Date:** 2026-07-02
**Status:** Approved

---

## Goal

Researchers keep per-mouse experiment metadata (ID, cage, line, sex, group, dates, dropout reason, etc.) in an Excel sheet, separate from FlowJo workspaces. This design adds two composable functions to read that metadata into a clean tibble and to join it onto any experimental data tibble (typically FACS data from `facs_read_wsp()`) by a shared identifier column, so downstream analysis functions (`analysis_*`) can filter/group by mouse-level attributes.

This introduces a new domain, `meta_`, since the reader and joiner are not FACS-specific — they apply to any tabular experimental data keyed by a subject identifier.

---

## Architecture

New file `R/meta_wrangle.R` containing two exported functions, no shared unexported helpers (each is simple enough standalone).

`janitor` added to `Imports` in `DESCRIPTION`, for `clean_names()` and `remove_empty()`.

CLAUDE.md's domain table gets a new row:

| Domain | Purpose |
|---|---|
| `meta_` | Experiment/subject metadata import and annotation |

---

## Public API

### `meta_read(path, sheet = 1)`

| Argument | Type | Default | Notes |
|---|---|---|---|
| `path` | chr | — | Path to `.xlsx` metadata file |
| `sheet` | chr or int | `1` | Sheet name or index, passed to `readxl::read_excel()` |

Returns a cleaned tibble, invisibly.

**Processing steps, in order:**

1. `readxl::read_excel(path, sheet = sheet)`
2. `janitor::remove_empty(c("rows", "cols"))` — drops fully-blank rows/columns left over from manual Excel editing
3. `janitor::clean_names()` — standardizes all column names to `snake_case`
4. Rename `mouse_id` back to `mouse_ID` — this is the literal keyword name FlowJo stores and that `facs_read_wsp(keywords = "mouse_ID")` produces; `meta_annotate()` joins on it verbatim, so it must not be mangled by step 3. If no `mouse_id` column exists after cleaning, skip this step (the sheet may not have that column).
5. Trim whitespace on every character column: `dplyr::mutate(dplyr::across(dplyr::where(is.character), stringr::str_trim))`
6. Coerce every `POSIXct`/`dttm` column to `Date`: `dplyr::mutate(dplyr::across(dplyr::where(lubridate::is.POSIXct), as.Date))` — implemented with base `as.Date()`; column detection via `inherits(x, "POSIXct")` to avoid adding a `lubridate` dependency for a one-line check
7. Coerce `group` column to a factor via `factor()`, if a `group` column is present after cleaning (levels = sorted unique values; no fixed level order assumed since group labels vary by experiment)

### `meta_annotate(data, meta, by = "mouse_ID")`

| Argument | Type | Default | Notes |
|---|---|---|---|
| `data` | tibble | — | Any tibble to annotate, e.g. `facs_read_wsp(...)$data` |
| `meta` | tibble | — | Metadata tibble, typically from `meta_read()` |
| `by` | chr | `"mouse_ID"` | Join column name; must exist in both `data` and `meta` |

Returns the joined tibble, invisibly.

**Processing steps, in order:**

1. If `by` is not a column in `data`, or not a column in `meta`: `stop()` with a message naming which side is missing it.
2. Compute `intersect(names(data), names(meta))` excluding `by`. If non-empty: `stop()` listing the colliding column names and asking the caller to rename/drop them before joining.
3. `dplyr::left_join(data, meta, by = by)`.
4. Compute `setdiff(unique(data[[by]]), unique(meta[[by]]))`. If non-empty: `warning()` listing these unmatched values (their rows keep `NA` for all meta columns).
5. Return the joined tibble, invisibly.

---

## Error handling

| Condition | Behaviour |
|---|---|
| `by` column missing from `data` or `meta` | `stop()`, names which side |
| Non-`by` column name collision between `data` and `meta` | `stop()`, lists colliding names |
| Value of `by` present in `data` but absent from `meta` | `warning()`, lists unmatched values; row kept with `NA` meta columns |
| Value of `by` present in `meta` but absent from `data` | Silently dropped (meta may cover more subjects than were FACS-processed) |
| `meta_read()` sheet has no `mouse_id`-cleaned column | Step 4 skipped, no error — not every metadata sheet is mouse-keyed |
| `meta_read()` sheet has no `group` column | Step 7 skipped, no error |

---

## Testing

New file `tests/testthat/test-meta_wrangle.R`. Uses the real fixture `tests/fixtures/meta_minimal.xlsx` (no mocking, consistent with project convention for file-reading functions).

**`meta_read()`:**
- Returns a tibble with expected columns, `mouse_ID` present and exactly named (not `mouse_id`)
- `dob`, `start_date`, `bmt_date`, `death_date` (note: `DOB` is cleaned to `dob` like any other column — only `mouse_ID` is exempted from `clean_names()`) are class `Date`
- `group` is a factor
- No fully-empty rows/columns remain
- Character columns have no leading/trailing whitespace

**`meta_annotate()`** — unit tests with small inline tibbles:
- Successful join adds meta columns to data
- Warning issued listing `by` values in `data` missing from `meta`
- Error when `by` column missing from either side
- Error when a non-`by` column collides between `data` and `meta`

**End-to-end:**
- `facs_read_wsp(testthat::test_path("fixtures", "minimal.wsp"), keywords = "mouse_ID")$data |> meta_annotate(meta_read(testthat::test_path("fixtures", "meta_minimal.xlsx")))` — both fixture FCS files carry keyword `mouse_ID = "26-1-17"`, which exists in `meta_minimal.xlsx`, proving the real pipeline joins correctly with no warnings.

---

## Dependencies

| Package | Change |
|---|---|
| `janitor` | Add to `Imports` |
| `readxl` | Add to `Imports` (not currently listed, needed for `meta_read()`) |
| `dplyr`, `stringr`, `tibble` | No change |

---

## Roxygen

Both functions get full `@param`/`@returns`/`@export` per CLAUDE.md convention. `@examples` use `\dontrun{}` since both require external files.
