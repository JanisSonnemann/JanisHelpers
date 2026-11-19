# ELISA: reshape wide-by-date sheets before prep -- design

Date: 2026-11-19. Size: small change.

## Intent

Longitudinal ELISA readouts (e.g. anti-MPO titer, `anti-mpo-titer_26-1_wide.xlsx`)
arrive as one row per mouse with one column per sampling date plus an end-point
column (`final`). The prep -> QC -> export pipeline already handles repeated dated
measurements, but only in the curated long layout (sheet `values`: `mouse_ID`,
`date`, one column per analyte; empty `date` = end-point). Add an **optional**
step that converts the wide file into a curated file; the pipeline itself is not
changed.

## Goals / non-goals

- Goal: one optional line in the prep script turns a wide file into a curated
  file that `elisa_prep_run()` reads unchanged.
- Non-goal: changing `elisa_read_curated()`, `elisa_prep_*()`, the QC report or
  the database export.
- Anti-goal: the converter breaking the stale check -- rerunning the script with
  `final = TRUE` must not make `elisa_prep_load()` refuse the saved prep.

## Approaches

- **A. Converter writes a curated `.xlsx` next to the raw file (picked).** Pipeline
  untouched; the long file can be opened and checked; clear import unit name.
  Costs a second file and moving `writexl` to Imports.
- B. Reshape in memory and pass a table to `elisa_prep_run()`. Rejected: prep
  save/load, the stale check and the QC report all work from a source file.
- C. Teach `elisa_read_curated()` to detect the wide layout. Rejected: changes
  the core reader and gives it two layouts to guess between.

## Design

New file `R/elisa_reshape.R`, exported function:

```r
elisa_reshape_wide(path, analyte = NULL, tissue = NULL, id_col = "Lab ID", out = NULL,
                   end_labels = c("final", "end", "endpoint"),
                   unit = NA, lloq = NA, uloq = NA)
```

- Reads the first sheet of `path`. `id_col` -> `mouse_ID` (error if missing).
- Optional `analyte` and `tissue` columns (header matched case-insensitively).
  With an `analyte` column, rows are mouse x analyte (x tissue): each analyte
  becomes its own column in the long file, so one wide file can hold several
  analytes. Without it, the `analyte` argument names the single analyte (error
  if neither is given). The column wins; the argument only fills empty cells.
- A `tissue` column is carried into the long file's `tissue` column; the
  `tissue` argument fills a missing column or empty cells. With neither, no
  `tissue` column is written and tissue comes from `elisa_prep_run(tissue = )`
  as today.
- Repeated mouse x tissue x analyte rows are replicates: numbered in row order
  and written as repeated mouse x tissue x date rows, which the curated reader
  already reads as replicates.
- Timepoint columns: headers that are dates (Excel serial numbers as readxl
  returns them, or `YYYY-MM-DD` text) -> that date; headers matching
  `end_labels` (case-insensitive) -> empty `date` (end-point). Every other
  column (`ID`, `group`, ...) is dropped and listed in a message.
- Errors if no timepoint column is found, or a date appears twice.
- Cells are passed through as text (numbers, `OOR<`, `OOR>`, empty), so the
  curated reader keeps doing all validation.
- Output: `out` defaults to `path` with a trailing `_wide` replaced by `_long`
  (else `_long` appended), same folder. Sheet `values` = `mouse_ID`,
  optional `tissue`, `date`, one column per analyte; rows ordered by tissue,
  date (end-point last), then input row order.
  If any of `unit`/`lloq`/`uloq` is given, a `range` sheet is written too;
  each is a single value (all analytes) or a vector named by analyte.
- Only rewrites when `out` is missing or older than `path`; otherwise messages
  "up to date" and does nothing. This keeps the prep stale check intact.
- Returns `out` invisibly, so the script can do
  `path <- elisa_reshape_wide(raw_path, analyte = "anti-MPO-titer")`.

Other changes:

- `DESCRIPTION`: `writexl` moves from Suggests to Imports; CLAUDE.md
  dependency lists updated.
- `inst/scripts/elisa_prep_template.R`: commented optional line in the assay
  block showing the reshape call.
- CLAUDE.md data-structures map: one line for `elisa_reshape_wide()`.
- Tests `tests/testthat/test-elisa_reshape.R`: wide fixture built with writexl
  in a temp dir -> reshape -> `elisa_read_curated()` round trip (dates,
  end-point NA, dropped columns, NA cells); multi-analyte + tissue columns
  -> one column per analyte and per-row tissue; replicate rows; no-rewrite
  when up to date; errors
  for missing `id_col` and no timepoint columns.

## Builder's calls (mine, correctable)

- Function name `elisa_reshape_wide()`; with no `analyte` column the
  `analyte` argument is required (no guessing from the file name).
- Optional `unit`/`lloq`/`uloq` args, so range info isn't lost by hand-editing a
  generated file that would be overwritten.
- `analyte` / `tissue` columns are optional and matched in any case; the
  arguments are fallbacks for missing columns or empty cells.
