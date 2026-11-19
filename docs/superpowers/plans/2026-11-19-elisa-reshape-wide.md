# ELISA wide-to-long reshape Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `elisa_reshape_wide()` turns a one-row-per-mouse, one-column-per-date ELISA sheet into a curated long `.xlsx` that `elisa_prep_run()` reads unchanged.

**Architecture:** A pure internal function reshapes the read sheet into the curated `values` table (+ optional `range` table); the exported function wraps it with file reading, output naming, write-only-when-changed, and returns the output path. No existing pipeline function changes.

**Tech Stack:** R, readxl, writexl (moves to Imports), dplyr/tidyr/purrr/stringr/glue, testthat 3.

**Spec:** `docs/superpowers/specs/2026-11-19-elisa-reshape-wide-design.md`

## Global Constraints

- Follow CLAUDE.md: base pipe `|>`, every external call `pkg::fn()`, non-ASCII escaped in R strings, plain ASCII in roxygen, `@param`/`@returns`/`@export` on the exported function, `\dontrun{}` example (needs a file).
- File `R/elisa_reshape.R`, tests `tests/testthat/test-elisa_reshape.R`.
- Signature exactly: `elisa_reshape_wide(path, analyte = NULL, tissue = NULL, id_col = "Lab ID", out = NULL, end_labels = c("final", "end", "endpoint"), unit = NA, lloq = NA, uloq = NA)`.
- Output sheet `values` columns in order: `mouse_ID`, `tissue` (only when a tissue column or argument exists), `date`, one column per analyte (first-appearance order). End-point rows have empty `date`.
- `devtools::check()`: 0 errors, 0 warnings.

## Review Focus

1. Changing an argument (e.g. adding `lloq`) without touching the wide file must take effect -> write decision is content-based: write when `out` is missing or its content differs from the new tables; otherwise message "up to date" and leave the file (and its mtime) alone. (Spec said "older than path"; content comparison covers that case and the argument-change case, and still never touches an unchanged file. Update the spec line in Task 2.)
2. Numeric mouse IDs or ID columns read by Excel as numbers must stay `"21560"`, not `"21560.0"` -> read with `col_types = "text"`.
3. Two headers mapping to the same timepoint (two `final` columns, a serial and an ISO header for the same day) -> error naming both headers.
4. Numbers must stay numbers in the written file (no "number stored as text") -> an analyte column whose non-empty cells all parse as numeric is written numeric; otherwise text (keeps `OOR<`/`OOR>`).
5. Rows with an empty ID cell (padding rows) are dropped silently, like the curated reader.

---

### Task 1: Pure reshape `elisa_reshape_wide_table_()`

**Files:**
- Create: `R/elisa_reshape.R`
- Test: `tests/testthat/test-elisa_reshape.R`

**Interfaces:**
- Produces: `elisa_reshape_wide_table_(raw, analyte = NULL, tissue = NULL, id_col = "Lab ID", end_labels = c("final", "end", "endpoint"), unit = NA, lloq = NA, uloq = NA)` -> `list(values = <tibble>, range = <tibble or NULL>, dropped = <character>)`. `raw` is a tibble of character columns (as from `readxl::read_excel(col_types = "text")`). `values$date` is class `Date` (NA = end-point). `range` (only when any of `unit`/`lloq`/`uloq` is non-NA) has first column `...1` with rows `"LLOQ"`, `"ULOQ"`, `"unit"` and one column per analyte.
- Produces: `parse_wide_date_header_(x)` -> `Date` vector; Excel serial text (`"46147"`) or `YYYY-MM-DD` -> date, anything else -> NA.

- [ ] **Step 1: Write failing tests** (build `raw` with `tibble::tibble(..., .name_repair = "minimal")`, all character)

```r
test_that("single analyte: date headers become rows, final is end-point, extras dropped", {
  raw <- tibble::tibble(`Lab ID` = c("26-1-1", "26-1-2"), ID = c("21560", "21561"),
    group = c("irrel", "MPO409"), `46147` = c("2.6658", "2.6729"),
    `2026-05-19` = c("2.6712", NA), final = c("2.5228", "OOR<"))
  r <- elisa_reshape_wide_table_(raw, analyte = "anti-MPO-titer")
  expect_named(r$values, c("mouse_ID", "date", "anti-MPO-titer"))
  expect_equal(r$values$date, as.Date(c("2026-05-05", "2026-05-05", "2026-05-19", "2026-05-19", NA, NA)))
  expect_equal(r$values$mouse_ID, rep(c("26-1-1", "26-1-2"), 3))
  expect_equal(r$values$`anti-MPO-titer`, c("2.6658", "2.6729", "2.6712", NA, "2.5228", "OOR<"))
  expect_setequal(r$dropped, c("ID", "group"))
  expect_null(r$range)
})
```

Plus, each its own `test_that`:
- `analyte`/`tissue` columns (header `Analyte`, `TISSUE` -- any case): two analytes x two mice x tissue `serum` -> `values` named `mouse_ID, tissue, date, TNFa, IFNg`, one row per mouse x date, tissue column `"serum"`.
- `tissue` argument fills empty tissue cells and a missing column; with neither, no `tissue` column.
- `analyte` argument fills empty `analyte` cells; neither column nor argument -> error matching `"analyte"`.
- Replicates: mouse `26-1-1` twice for analyte `TNFa` -> two `values` rows for that mouse x date (row order kept).
- Numeric column: all-numeric analyte column is `numeric`; column with `"OOR<"` is `character`.
- `range`: `lloq = c(TNFa = 7.8, IFNg = 5)`, `unit = "pg/ml"` -> `range[[1]] == c("LLOQ", "ULOQ", "unit")`, `range$TNFa == c("7.8", NA, "pg/ml")`.
- Errors: missing `id_col` (message names it); no timepoint column; duplicate timepoint (`46147` and `2026-05-05`; two `final`-like headers `final` + `End`) with both headers in the message.
- Rows with empty `Lab ID` dropped.

- [ ] **Step 2: Run, expect FAIL** -- `devtools::test(filter = "elisa_reshape")` -> "could not find function".

- [ ] **Step 3: Implement** both functions in `R/elisa_reshape.R`. Approach: locate `id_col` exactly, `analyte`/`tissue` via `tolower(names) == "analyte"/"tissue"`; classify remaining headers with `parse_wide_date_header_()` and `tolower(trimws(h)) %in% tolower(end_labels)`; others -> `dropped`. Number replicates with `dplyr::row_number()` by `mouse_ID x tissue x analyte`; `tidyr::pivot_longer()` the timepoint columns to `date`/`value`, then `tidyr::pivot_wider(names_from = analyte, values_from = value)` keyed by `mouse_ID, tissue, date, replicate`; drop `replicate`; arrange by tissue, date (`NA` last), replicate, then first-appearance mouse order. Convert all-numeric analyte columns with `as.numeric()`. `unit`/`lloq`/`uloq`: length-1 unnamed -> recycled to all analytes; named -> matched by analyte (unknown names -> error).

- [ ] **Step 4: Run, expect PASS** -- `devtools::test(filter = "elisa_reshape")`.

- [ ] **Step 5: Commit** -- `git add R/elisa_reshape.R tests/testthat/test-elisa_reshape.R && git commit -m "feat(elisa): pure wide-by-date reshape helper"`

### Task 2: Exported `elisa_reshape_wide()` + wiring

**Files:**
- Modify: `R/elisa_reshape.R`, `DESCRIPTION` (writexl Suggests -> Imports, between `uwot` and `xml2`), `NAMESPACE`/`man/` via `devtools::document()`, `inst/scripts/elisa_prep_template.R`, `CLAUDE.md`, spec (Review Focus 1 line)
- Test: `tests/testthat/test-elisa_reshape.R`

**Interfaces:**
- Consumes: `elisa_reshape_wide_table_()` from Task 1.
- Produces: `elisa_reshape_wide(...)` (signature in Global Constraints) -> invisible character path of `out`.

- [ ] **Step 1: Write failing tests** (fixture: `writexl::write_xlsx(list(Tabelle1 = raw), tempfile("titer_wide", fileext = ".xlsx"))`, `skip_if_not_installed("writexl")`)

- `"round trip through elisa_read_curated"`: `p <- elisa_reshape_wide(f, analyte = "anti-MPO-titer")`; `expect_equal(basename(p), sub("_wide", "_long", basename(f)))`; `d <- suppressMessages(elisa_read_curated(p, tissue = "serum"))`; expect 6 rows, `sum(is.na(d$collected_at)) == 2`, `d$result_status` contains `"OOR<"`, `unique(d$analyte) == "anti-MPO-titer"`.
- `"out without _wide gets _long appended"` and explicit `out` honoured.
- `"unchanged content is not rewritten"`: call twice; second call messages `"up to date"` and `file.mtime(p)` unchanged (set the first mtime back with `Sys.setFileTime(p, Sys.time() - 60)` before the second call).
- `"changed arguments rewrite"`: second call with `lloq = 0.1` rewrites; `readxl::excel_sheets(p)` contains `"range"`.
- `"dropped columns are messaged"`: `expect_message(..., "ID, group")`.

- [ ] **Step 2: Run, expect FAIL** -- `devtools::test(filter = "elisa_reshape")`.

- [ ] **Step 3: Implement `elisa_reshape_wide()`**: read first sheet with `readxl::read_excel(path, sheet = 1, col_types = "text", .name_repair = "minimal")`; call the helper; default `out` = `sub("_wide$", "_long", stem)` or `paste0(stem, "_long")` + `.xlsx` in `dirname(path)`; compare with existing `out` by reading its sheets back as text (`values`, and `range` if present) and comparing with the new tables formatted as text (dates as `YYYY-MM-DD`); write with `writexl::write_xlsx(list(values = ..., range = ...))` (omit `range` when NULL) only if missing/different; message dropped columns, rows/analytes written or "up to date". Roxygen per CLAUDE.md (plain ASCII), `@returns` the output path invisibly.

- [ ] **Step 4: Wire up** -- move `writexl` to Imports in `DESCRIPTION`; in CLAUDE.md add `writexl` to the "Current Imports" list, remove it from the Suggests sentence, add one data-structures bullet: "`elisa_reshape_wide()` -- optional step before ELISA prep: wide sheet (one row per mouse [x analyte x tissue], one column per date, `final` = end-point) -> curated `<stem>_long.xlsx` next to it; rewrites only when content changed (keeps the prep stale check); returns the path."; in `inst/scripts/elisa_prep_template.R` after the `path <-` line add:

```r
# Wide sheet (one column per date, `final` = end-point)? Reshape it first:
# path <- elisa_reshape_wide(file.path(exp_path, "raw", "ELISA", "anti-mpo-titer_26-1_wide.xlsx"),
#                            analyte = "anti-MPO-titer")
```

and update spec bullet "Only rewrites when `out` is missing or older than `path`" to the content-based rule.

- [ ] **Step 5: Run** `devtools::document()` then `devtools::test(filter = "elisa")` -> all PASS.

- [ ] **Step 6: Real-file smoke test** (not committed): run `elisa_reshape_wide()` on the 26-1 wide file with `out` in `$CLAUDE_JOB_DIR/tmp`, then `elisa_read_curated(out, tissue = "serum")` -> 17 mice, 68 rows, 17 end-point rows; diff values against `anti-mpo-titer_26-1_long.xlsx`.

- [ ] **Step 7: `devtools::check()`** -> 0 errors, 0 warnings, only NOTEs listed in CLAUDE.md.

- [ ] **Step 8: Commit** -- `git add -A R tests man NAMESPACE DESCRIPTION CLAUDE.md inst/scripts docs && git commit -m "feat(elisa): elisa_reshape_wide() optional step before prep"`
