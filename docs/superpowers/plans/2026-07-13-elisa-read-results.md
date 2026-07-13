# `elisa_read_results()` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `elisa_read_results()`, a new exported function that extracts per-replicate back-calculated cytokine concentrations from a multiplex bead ELISA (Luminex/Bio-Plex) `Results_<cytokine>.xlsx` export, for downstream analysis.

**Architecture:** One exported function `elisa_read_results()` in a new `R/elisa_read.R` file. It locates the sample-level sheet (name ends in `"Unknowns"`), reads it with `readxl::read_excel()`, and reshapes the per-replicate columns into a tidy long tibble. Introduces a new `elisa_` domain.

**Tech Stack:** R, `readxl` and `stringr` (both already dependencies — no new `Imports`), `dplyr`/`glue` (already dependencies), `testthat` 3rd edition.

## Global Constraints

- Pipe: use `|>` (base pipe) everywhere; no `%>%`.
- Namespacing: every external call inside the function body is `pkg::fn()` — no bare `readxl::read_excel()`/`dplyr::filter()`/etc. calls without a namespace. Bare column names inside `dplyr` verbs are fine (tidy eval), not a namespacing violation.
- Every exported function needs `@param`, `@returns`, `@export`; `@examples` uses `\dontrun{}` since it needs a real `.xlsx` file.
- `elisa_read_results()` returns its tibble via `invisible()` — it's an import-plus-message function, matching `meta_read()`'s convention (see `CLAUDE.md`'s `invisible()` rule).
- No mocking of `readxl` — test against the real fixtures in `tests/fixtures/`, matching the package's existing fixture-based testing philosophy for `fcexpr`/`CytoML`.
- Target after all tasks: `devtools::check()` → 0 errors, 0 warnings.
- Run `devtools::load_all()` before interactively testing; run the specific `testthat::test_file()` after each implementation step.

---

### Task 1: Declare the `elisa_` domain in `CLAUDE.md`

**Files:**
- Modify: `CLAUDE.md`

**Interfaces:**
- Produces: `elisa_` domain declared in `CLAUDE.md`, ready for `elisa_read_results()` to follow the naming convention required by the package's own rules.

- [ ] **Step 1: Add the `elisa_` domain to `CLAUDE.md`'s domain list**

Edit `CLAUDE.md`, in the "File and function naming" section, replacing:

```
- `facs_` — FlowJo / flow cytometry
- `analysis_` — statistical summaries and tests
- `plot_` — ggplot2 visualization of FACS/analysis data
- `report_` — RMarkdown rendering
- `meta_` — experiment/subject metadata import and annotation
- `wrangle_` — general data wrangling [stub — no functions yet]
- `db_` — database access [stub — no functions yet]
```

with:

```
- `facs_` — FlowJo / flow cytometry
- `analysis_` — statistical summaries and tests
- `plot_` — ggplot2 visualization of FACS/analysis data
- `report_` — RMarkdown rendering
- `meta_` — experiment/subject metadata import and annotation
- `elisa_` — multiplex bead ELISA (Luminex/Bio-Plex) import
- `wrangle_` — general data wrangling [stub — no functions yet]
- `db_` — database access [stub — no functions yet]
```

- [ ] **Step 2: Verify the package still loads**

Run: `Rscript -e "devtools::load_all()"`
Expected: loads with no errors (doc-only change so far).

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: declare elisa_ domain for multiplex ELISA import"
```

---

### Task 2: `elisa_read_results()` — read and reshape the Unknowns sheet

**Files:**
- Create: `R/elisa_read.R`
- Create: `tests/testthat/test-elisa_read.R`

**Interfaces:**
- Produces: exported `elisa_read_results(path, cytokine = NULL)`, returning (invisibly) a tibble with columns `cytokine` (chr), `sample_id` (chr), `replicate` (int), `value` (dbl), `unit` (chr), `result_status` (chr).

Verified against the real fixtures during design (both files have an `Unknowns`-suffixed sheet with 56 real data rows, no `NA`-`Sample` padding rows):
- `tests/fixtures/Results_IL-17A.xlsx`: sheet literally named `"Unknowns"`; `Backcalc` column is `chr` (contains `"OOR<"` text for 20 of 56 rows, across 10 distinct `sample_id`s: `25-7-5`, `25-7-11`, `25-7-12`, `25-7-13`, `25-7-14`, `25-7-20`, `25-7-21`, `25-7-24`, `25-7-25`, `25-7-28`); `result_status` for those rows is `"OOR<"`. Unit parses to `"pg/ml"`.
- `tests/fixtures/Results_TNFa.xlsx`: sheet named `"TNFa_Unknowns"`; `Backcalc` column is fully numeric (no `"OOR<"`/`"OOR>"` text); 6 rows have `result_status == "<LLOQ"` but still carry a non-`NA` numeric `value` (e.g. `sample_id == "25-7-3"`, `replicate == 1`, `value == 1.303633`) — `<LLOQ` is a below-quantification-limit flag, distinct from out-of-range (`OOR`). Unit parses to `"pg/ml"`.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-elisa_read.R`:

```r
library(testthat)
library(JanisHelpers)

il17a_path <- testthat::test_path("../fixtures/Results_IL-17A.xlsx")
tnfa_path  <- testthat::test_path("../fixtures/Results_TNFa.xlsx")
il17a_skip_msg <- "Results_IL-17A.xlsx fixture not available"
tnfa_skip_msg  <- "Results_TNFa.xlsx fixture not available"

test_that("elisa_read_results() returns the expected columns and types", {
  skip_if_not(file.exists(il17a_path), il17a_skip_msg)
  dat <- elisa_read_results(il17a_path)
  expect_s3_class(dat, "tbl_df")
  expect_equal(
    names(dat),
    c("cytokine", "sample_id", "replicate", "value", "unit", "result_status")
  )
  expect_type(dat$cytokine, "character")
  expect_type(dat$sample_id, "character")
  expect_type(dat$replicate, "integer")
  expect_type(dat$value, "double")
  expect_type(dat$unit, "character")
  expect_type(dat$result_status, "character")
  expect_equal(nrow(dat), 56L)
})

test_that("elisa_read_results() derives cytokine from filename (no cytokine-prefixed sheet)", {
  skip_if_not(file.exists(il17a_path), il17a_skip_msg)
  dat <- elisa_read_results(il17a_path)
  expect_equal(unique(dat$cytokine), "IL-17A")
})

test_that("elisa_read_results() derives cytokine from filename (cytokine-prefixed sheet)", {
  skip_if_not(file.exists(tnfa_path), tnfa_skip_msg)
  dat <- elisa_read_results(tnfa_path)
  expect_equal(unique(dat$cytokine), "TNFa")
  expect_equal(nrow(dat), 56L)
})

test_that("elisa_read_results() parses the unit from the Backcalc column header", {
  skip_if_not(file.exists(il17a_path), il17a_skip_msg)
  dat <- elisa_read_results(il17a_path)
  expect_equal(unique(dat$unit), "pg/ml")
})

test_that("elisa_read_results() sets value to NA for OOR rows and keeps result_status", {
  skip_if_not(file.exists(il17a_path), il17a_skip_msg)
  dat <- elisa_read_results(il17a_path)
  oor <- dplyr::filter(dat, result_status == "OOR<")
  expect_equal(nrow(oor), 20L)
  expect_true(all(is.na(oor$value)))
  expect_setequal(
    unique(oor$sample_id),
    c("25-7-5", "25-7-11", "25-7-12", "25-7-13", "25-7-14",
      "25-7-20", "25-7-21", "25-7-24", "25-7-25", "25-7-28")
  )
})

test_that("elisa_read_results() keeps a non-NA value for <LLOQ rows (below quantification, not out of range)", {
  skip_if_not(file.exists(tnfa_path), tnfa_skip_msg)
  dat <- elisa_read_results(tnfa_path)
  lloq <- dplyr::filter(dat, result_status == "<LLOQ")
  expect_equal(nrow(lloq), 6L)
  expect_true(all(!is.na(lloq$value)))
  row <- dplyr::filter(lloq, sample_id == "25-7-3", replicate == 1L)
  expect_equal(row$value, 1.303633, tolerance = 1e-6)
})

test_that("elisa_read_results() uses the cytokine argument when supplied, overriding the filename", {
  skip_if_not(file.exists(tnfa_path), tnfa_skip_msg)
  dat <- elisa_read_results(tnfa_path, cytokine = "Custom")
  expect_equal(unique(dat$cytokine), "Custom")
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-elisa_read.R')"`
Expected: FAIL with `could not find function "elisa_read_results"`.

- [ ] **Step 3: Implement `elisa_read_results()`**

Create `R/elisa_read.R`:

```r
#' Read back-calculated cytokine concentrations from a multiplex ELISA results file
#'
#' @description
#' Reads the sample-level ("Unknowns") sheet of a multiplex bead ELISA
#' (Luminex/Bio-Plex) \code{Results_<cytokine>.xlsx} export and extracts the
#' per-replicate back-calculated concentration. The samples sheet is located
#' by matching a sheet name ending in \code{"Unknowns"}, since some exports
#' name it plainly \code{"Unknowns"} while others prefix it with the
#' cytokine (e.g. \code{"TNFa_Unknowns"}). Rows with no \code{Sample} value
#' are dropped (a defensive safeguard against blank padding rows, a pattern
#' seen in these workbooks' other sheets). Out-of-range replicates
#' (\code{"OOR<"}, \code{"OOR>"} in the back-calculated column) become
#' \code{NA} in \code{value}, with the original flag preserved in
#' \code{result_status}.
#'
#' @param path path to a \code{Results_<cytokine>.xlsx} workbook
#' @param cytokine character; cytokine label for the \code{cytokine} output
#'   column. Default \code{NULL} derives it from \code{basename(path)} by
#'   stripping a leading \code{Results_} and trailing \code{.xlsx}.
#'
#' @returns tibble, one row per sample x replicate, with columns
#'   \code{cytokine}, \code{sample_id}, \code{replicate}, \code{value},
#'   \code{unit}, \code{result_status}; returned invisibly -- assign the
#'   result explicitly
#' @export
#'
#' @examples
#' \dontrun{
#'   elisa_read_results("Results_IL-17A.xlsx")
#' }
elisa_read_results <- function(path, cytokine = NULL) {
  if (is.null(cytokine)) {
    cytokine <- basename(path) |>
      stringr::str_remove("^Results_") |>
      stringr::str_remove("\\.xlsx$")
  }

  sheet_names <- readxl::excel_sheets(path)
  unknowns_sheet <- sheet_names[grepl("Unknowns$", sheet_names)]

  if (length(unknowns_sheet) != 1L) {
    stop(glue::glue(
      "Expected exactly one sheet ending in 'Unknowns' in '{path}', found ",
      "{length(unknowns_sheet)}: {paste(sheet_names, collapse = ', ')}."
    ))
  }

  raw <- readxl::read_excel(path, sheet = unknowns_sheet)

  backcalc_col <- names(raw)[grepl("^Backcalc", names(raw))]
  status_col   <- names(raw)[grepl("^Result.*Status", names(raw))]

  if (length(backcalc_col) != 1L) {
    stop(glue::glue(
      "Expected exactly one 'Backcalc' column in sheet '{unknowns_sheet}' ",
      "of '{path}', found {length(backcalc_col)}."
    ))
  }
  if (length(status_col) != 1L) {
    stop(glue::glue(
      "Expected exactly one 'Result Status' column in sheet '{unknowns_sheet}' ",
      "of '{path}', found {length(status_col)}."
    ))
  }

  unit <- stringr::str_extract(backcalc_col, "(?<=\\()[^)]+(?=\\))")
  if (is.na(unit)) {
    warning(glue::glue(
      "Could not parse a unit from Backcalc column header '{backcalc_col}'; ",
      "`unit` will be NA."
    ))
  }

  result <- raw |>
    dplyr::select(
      sample_id     = 1,
      replicate     = Rep,
      value         = dplyr::all_of(backcalc_col),
      result_status = dplyr::all_of(status_col)
    ) |>
    dplyr::filter(!is.na(sample_id)) |>
    dplyr::mutate(
      cytokine      = cytokine,
      replicate     = as.integer(replicate),
      value         = suppressWarnings(as.numeric(value)),
      unit          = unit,
      result_status = stringr::str_trim(result_status)
    ) |>
    dplyr::relocate(cytokine, sample_id, replicate, value, unit, result_status)

  message(glue::glue(
    "\nExtraction Summary",
    "\n----------------------------------------------",
    "\nCytokine:          {cytokine}",
    "\nNumber of samples: {dplyr::n_distinct(result$sample_id)}",
    "\nNumber of rows:    {nrow(result)}",
    "\n----------------------------------------------\n"
  ))

  invisible(result)
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-elisa_read.R')"`
Expected: PASS, 7 tests.

- [ ] **Step 5: Commit**

```bash
git add R/elisa_read.R tests/testthat/test-elisa_read.R
git commit -m "feat: add elisa_read_results() for multiplex ELISA back-calculated values"
```

---

### Task 3: Docs and verification

**Files:**
- Modify: `CLAUDE.md`
- Create/Modify (generated): `NAMESPACE`, `man/elisa_read_results.Rd`

**Interfaces:**
- Consumes: the finished `elisa_read_results()` from Task 2.
- Produces: exported symbol registered in `NAMESPACE`, `man/elisa_read_results.Rd` generated, `CLAUDE.md` documents the new data structure, `devtools::check()` passes with 0 errors/0 warnings.

- [ ] **Step 1: Add a Data structures entry to CLAUDE.md**

In `CLAUDE.md`'s "Data structures" section, after the `plot_bubble_fc()` entry and before the closing `---` of that section, add:

```
### `elisa_read_results()` — multiplex ELISA back-calculated concentrations
- **Input**: `path` = path to a `Results_<cytokine>.xlsx` multiplex bead
  ELISA (Luminex/Bio-Plex) export. `cytokine` = optional label override,
  default `NULL` derives it from the filename (`Results_<cytokine>.xlsx`).
- **Output**: long tibble, one row per `sample_id x replicate`, read from
  the sheet whose name ends in `"Unknowns"` (naming varies: plain
  `"Unknowns"` or `"<cytokine>_Unknowns"`).

| Column | Type | Notes |
|---|---|---|
| `cytokine` | chr | from filename or `cytokine` argument |
| `sample_id` | chr | raw `Sample` value, e.g. `"25-7-1"` -- not assumed to equal `mouse_ID` |
| `replicate` | int | from `Rep` |
| `value` | dbl | back-calculated concentration; `NA` if out of standard-curve range (`"OOR<"`/`"OOR>"`) |
| `unit` | chr | parsed from the Backcalc column header, e.g. `"pg/ml"` |
| `result_status` | chr | `"OK"` / `"OOR<"` / `"OOR>"` / `"<LLOQ"` etc. |

- Only the `Unknowns` sheet's per-replicate table is read; the `Curve` and
  `Standards` sheets, and the sheet's own per-sample summary columns
  (`Mean`/`SD`/`SEM`/`CV %`/`Final Result`), are out of scope -- use the
  package's `analysis_` functions to summarize across replicates instead.
```

- [ ] **Step 2: Regenerate documentation**

Run: `Rscript -e "devtools::document()"`
Expected: creates `man/elisa_read_results.Rd`, adds `export(elisa_read_results)` to `NAMESPACE`, no roxygen errors.

- [ ] **Step 3: Run the full test suite**

Run: `Rscript -e "devtools::test()"`
Expected: all tests pass, including the pre-existing suite (no regressions).

- [ ] **Step 4: Run `devtools::check()` and reconcile CLAUDE.md's "Known check output" section**

Run: `Rscript -e "devtools::check()"`
Expected: 0 errors, 0 warnings. Read the NOTEs output. If any new "R code for possible problems" NOTEs reference `elisa_read_results` (expected: bare column names like `sample_id`, `Rep`, `cytokine`, `replicate`, `value`, `unit`, `result_status` inside `dplyr` verbs -- the same tidy-eval false-positive pattern as every other domain), add one bullet for them to CLAUDE.md's "Known check output" section, following the existing bullet style (e.g. the `facs_read_fcs_gated` bullet).

- [ ] **Step 5: Commit**

```bash
git add CLAUDE.md NAMESPACE man/elisa_read_results.Rd
git commit -m "docs: document elisa_read_results() and reconcile check output notes"
```
