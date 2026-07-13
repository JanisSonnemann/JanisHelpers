# `elisa_read_results()` design

Date: 2026-07-13

## Problem

Multiplex bead ELISA (Luminex/Bio-Plex) runs produce one `Results_<cytokine>.xlsx`
workbook per cytokine, each with `Curve`, `Standards`, and `Unknowns` sheets (the
latter two prefixed with the cytokine name for some exports, e.g. `TNFa_Unknowns`,
but not for others, e.g. plain `Unknowns`). The `Unknowns` sheet holds two
side-by-side tables: a per-replicate table (`Sample`, `Rep`, `NetMFI`, `Backcalc
(pg/ml)`, `Result Status`) and a per-sample summary table (`Sample`, `n`, `Mean`,
`SD`, `SEM`, `CV %`, `Final Result`). There is currently no package function to
extract the back-calculated cytokine concentrations for downstream analysis.

## Domain

This is multiplex bead ELISA import, distinct from flow cytometry (`facs_`),
statistical summarization (`analysis_`), or metadata (`meta_`). Introduces a new
domain prefix: **`elisa_`**. `CLAUDE.md`'s domain table is updated to include it
before implementation, per the package's naming-convention rule.

## Function

`elisa_read_results(path, cytokine = NULL)` in `R/elisa_read.R`.

- **`path`**: path to a `Results_<cytokine>.xlsx` workbook.
- **`cytokine`**: optional cytokine label. If `NULL` (default), derived by
  stripping a `Results_` prefix and `.xlsx` suffix from `basename(path)`
  (e.g. `Results_IL-17A.xlsx` -> `"IL-17A"`). If supplied, used verbatim.

### Behavior

1. **Locate the samples sheet.** `grepl("Unknowns$", readxl::excel_sheets(path))`
   — matches both `"Unknowns"` and `"TNFa_Unknowns"` without assuming a fixed
   prefix. Error if zero or more than one sheet matches.
2. **Read and select columns.** Read the matched sheet with
   `readxl::read_excel()`. Keep the first `Sample` column (the sheet has a
   second `Sample` column further right, for the summary block — ignored),
   `Rep`, and the Backcalc/Result-Status columns. The latter two are located
   by regex on `names()` (`"^Backcalc"`, `"^Result.*Status"`) rather than by
   literal header string, since the real headers embed `\r\n`
   (`"Backcalc\r\n(pg/ml)"`, `"Result\r\nStatus"`) and matching on that
   literal text would be fragile.
3. **Drop padding rows.** Rows where `Sample` is `NA` are dropped — the sheet
   is padded with blank rows past the real data (observed: 1001 rows in
   `Curve`, 56 real vs many blank in `Unknowns`).
4. **Parse the unit.** Extract the parenthesized unit from the Backcalc
   column's header text via `stringr::str_extract()` (e.g. `"(pg/ml)"` ->
   `"pg/ml"`). Warn and set `NA` if no match is found.
5. **Coerce values.** `Backcalc` -> numeric via
   `suppressWarnings(as.numeric())`; non-numeric flags (`"OOR<"`, `"OOR>"`)
   become `NA`. `Result Status` is trimmed of embedded `\r\n` and kept as-is
   (`"OK"`, `"OOR<"`, `"OOR>"`, `"<LLOQ"`, etc.) — this is the existing status
   flag already present in the sheet, not re-derived from the Backcalc text.

### Output

Long tibble, one row per `Sample x Rep`:

| Column | Type | Notes |
|---|---|---|
| `cytokine` | chr | from filename or `cytokine` argument |
| `sample_id` | chr | raw `Sample` value, e.g. `"25-7-1"` — not assumed to equal `mouse_ID` |
| `replicate` | int | from `Rep` |
| `value` | dbl | back-calculated concentration; `NA` if out of standard-curve range |
| `unit` | chr | parsed from the Backcalc column header, e.g. `"pg/ml"` |
| `result_status` | chr | `"OK"` / `"OOR<"` / `"OOR>"` / `"<LLOQ"` etc., trimmed |

### Explicitly out of scope

- The `Curve` and `Standards` sheets (standard-curve fitting/QC data).
- The per-sample summary table (`Mean`/`SD`/`SEM`/`CV %`/`Final Result`) —
  downstream summarization across replicates is left to the package's
  existing `analysis_` functions, so the sheet's own precomputed mean isn't
  duplicated.
- Combining multiple cytokine files into one tibble — `elisa_read_results()`
  reads one file per call (mirrors `facs_import_wsp()`); callers combine
  multiple cytokines with `dplyr::bind_rows()`.
- Joining to `mouse_ID`/experiment metadata — `sample_id` is left as the raw
  value from the sheet; mapping to `mouse_ID` (if 1:1) is a caller concern,
  not assumed by this function.

### Error handling

- Zero or multiple sheets matching `"Unknowns$"` -> error.
- `Backcalc`/`Result Status` columns not found by regex -> error.

### Dependencies

None new. `readxl` and `stringr` are already in `Imports`.

## Testing

`tests/testthat/test-elisa_read.R`, smoke tests against the real
`tests/fixtures/Results_IL-17A.xlsx` and `tests/fixtures/Results_TNFa.xlsx`
fixtures (no mocking of `readxl`, consistent with the package's fixture-based
testing philosophy for `fcexpr`/`CytoML`):

- Required columns (`cytokine`, `sample_id`, `replicate`, `value`, `unit`,
  `result_status`) are present for both files.
- `cytokine` is correctly derived from filename for both the unprefixed
  (`Results_IL-17A.xlsx`) and prefixed-sheet-name (`Results_TNFa.xlsx`) cases.
- Rows with `result_status` of `"OOR<"`/`"OOR>"` have `value == NA`.
- `cytokine` argument override works when supplied explicitly.
