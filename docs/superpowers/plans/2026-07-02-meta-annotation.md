# meta_read() + meta_annotate() Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a new `meta_` domain with two functions — `meta_read()` to import and clean a per-subject experiment metadata spreadsheet, and `meta_annotate()` to left-join that metadata onto any experimental data tibble (e.g. FACS data from `facs_read_wsp()`) by a shared identifier column.

**Architecture:** One new file `R/meta_wrangle.R` with two exported functions and one small unexported helper. Tests in `tests/testthat/test-meta_wrangle.R`, using the real fixtures `tests/fixtures/meta_minimal.xlsx` and `tests/fixtures/minimal.wsp` (no mocking). `readxl` and `janitor` added to `Imports`. CLAUDE.md's domain table gets a new `meta_` row.

**Tech Stack:** R, `readxl`, `janitor`, `dplyr`, `stringr`, `glue`, `testthat` (edition 3).

## Global Constraints

- Every external function call inside an exported function body must be namespaced `pkg::fn()` — no bare `map()`, `filter()`, etc. (CLAUDE.md).
- Use the base pipe `|>`, not `%>%`.
- `dplyr::select(!col)` for deselection, never `-col`.
- Every exported function needs `@param`, `@returns`, and `@export` roxygen tags; `@returns` must not be empty.
- `@examples` must use `\dontrun{}` since both functions require external files.
- Functions that return tibbles primarily for later use but follow this package's import-function convention return `invisible(result)`.
- Run `devtools::document()` after any roxygen change, before committing.
- Target `devtools::check()`: 0 errors, 0 warnings (existing 2 notes are pre-existing and unrelated).
- The join key `mouse_ID` must be preserved with exact casing — it is the literal FlowJo keyword name and must match what `facs_read_wsp(keywords = "mouse_ID")` produces.
- Test fixtures `tests/fixtures/meta_minimal.xlsx` and `tests/fixtures/minimal.wsp` already exist locally but are untracked in git (consistent with existing project convention of not committing fixtures that may contain identifiable lab data — see `tests/fixtures/README.md`). Do not add them to git as part of this plan. All tests must `skip_if_not(file.exists(...))` so the suite passes even without the fixtures present.

---

### Task 1: `meta_read()`

**Files:**
- Create: `R/meta_wrangle.R`
- Modify: `DESCRIPTION` (add `readxl`, `janitor` to `Imports`)
- Modify: `CLAUDE.md` (add `meta_` domain row)
- Test: `tests/testthat/test-meta_wrangle.R`

**Interfaces:**
- Produces: `meta_read(path, sheet = 1)` — returns a cleaned tibble invisibly. Columns are `snake_case` except `mouse_ID` (exact casing preserved if present after cleaning). Any originally-`POSIXct` column is `Date`. `group` column (if present) is a `factor`.

- [ ] **Step 1: Add `readxl` and `janitor` to `DESCRIPTION` Imports**

Open `DESCRIPTION` and update the `Imports:` block to insert `readxl` and `janitor` in alphabetical order:

```
Imports:
    dplyr,
    fcexpr,
    glue,
    gt,
    gtsummary,
    janitor,
    purrr,
    readxl,
    rmarkdown,
    rstatix,
    stringr,
    tibble,
    tidyr,
    xml2,
    xfun
```

- [ ] **Step 2: Add the `meta_` domain row to `CLAUDE.md`**

In `CLAUDE.md`, find this block (around line 35-40):

```
**Domains** (use exactly these prefixes):
- `facs_` — FlowJo / flow cytometry
- `analysis_` — statistical summaries and tests
- `report_` — RMarkdown rendering
- `wrangle_` — general data wrangling [stub — no functions yet]
- `db_` — database access [stub — no functions yet]
```

Replace it with:

```
**Domains** (use exactly these prefixes):
- `facs_` — FlowJo / flow cytometry
- `analysis_` — statistical summaries and tests
- `report_` — RMarkdown rendering
- `meta_` — experiment/subject metadata import and annotation
- `wrangle_` — general data wrangling [stub — no functions yet]
- `db_` — database access [stub — no functions yet]
```

- [ ] **Step 3: Write the failing test for `meta_read()`**

Create `tests/testthat/test-meta_wrangle.R`:

```r
library(testthat)
library(JanisHelpers)

meta_path <- testthat::test_path("../fixtures/meta_minimal.xlsx")
wsp_path  <- testthat::test_path("../fixtures/minimal.wsp")
meta_skip_msg <- "meta_minimal.xlsx fixture not available"
wsp_skip_msg  <- "minimal.wsp fixture not available"

test_that("meta_read() returns a tibble with mouse_ID preserved verbatim", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta <- meta_read(meta_path)
  expect_s3_class(meta, "tbl_df")
  expect_true("mouse_ID" %in% names(meta))
  expect_false("mouse_id" %in% names(meta))
})

test_that("meta_read() coerces date columns to Date", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta <- meta_read(meta_path)
  date_cols <- intersect(names(meta), c("dob", "start_date", "bmt_date", "death_date"))
  expect_true(length(date_cols) > 0L)
  for (col in date_cols) {
    expect_s3_class(meta[[col]], "Date")
  }
})

test_that("meta_read() coerces group to a factor", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta <- meta_read(meta_path)
  expect_true("group" %in% names(meta))
  expect_s3_class(meta$group, "factor")
})

test_that("meta_read() trims whitespace and drops empty rows/cols", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta <- meta_read(meta_path)
  chr_cols <- names(meta)[vapply(meta, is.character, logical(1))]
  for (col in chr_cols) {
    vals <- stats::na.omit(meta[[col]])
    expect_equal(vals, stringr::str_trim(vals))
  }
  expect_true(all(!vapply(meta, function(x) all(is.na(x)), logical(1))))
})
```

- [ ] **Step 4: Run the test to verify it fails**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-meta_wrangle.R')"`
Expected: FAIL with `could not find function "meta_read"`

- [ ] **Step 5: Implement `meta_read()`**

Create `R/meta_wrangle.R`:

```r
is_posixct_ <- function(x) inherits(x, "POSIXct")

#' Read and clean an experiment metadata spreadsheet
#'
#' @description
#' Reads an Excel sheet of per-subject experiment metadata (e.g. mouse ID,
#' cage, sex, group, dates) and applies standard cleaning: blank rows and
#' columns removed, column names standardized to snake_case (except
#' \code{mouse_ID}, preserved verbatim so it stays joinable against FlowJo
#' keyword data such as \code{facs_read_wsp(keywords = "mouse_ID")}),
#' character columns trimmed of whitespace, date/time columns coerced to
#' \code{Date}, and a \code{group} column (if present) coerced to a factor.
#'
#' @param path path to \code{.xlsx} metadata file
#' @param sheet sheet name or index to read, default = 1
#'
#' @returns cleaned tibble, one row per subject, returned invisibly --
#'   assign the result explicitly
#' @export
#'
#' @examples
#' \dontrun{
#'   meta <- meta_read("meta.xlsx")
#' }
meta_read <- function(path, sheet = 1) {
  dat <- readxl::read_excel(path, sheet = sheet)
  dat <- janitor::remove_empty(dat, which = c("rows", "cols"))
  dat <- janitor::clean_names(dat)

  if ("mouse_id" %in% names(dat)) {
    dat <- dplyr::rename(dat, mouse_ID = mouse_id)
  }

  dat <- dplyr::mutate(
    dat,
    dplyr::across(dplyr::where(is.character), stringr::str_trim)
  )

  dat <- dplyr::mutate(
    dat,
    dplyr::across(dplyr::where(is_posixct_), as.Date)
  )

  if ("group" %in% names(dat)) {
    dat <- dplyr::mutate(dat, group = factor(group))
  }

  invisible(dat)
}
```

- [ ] **Step 6: Run `devtools::document()` then the test to verify it passes**

Run: `Rscript -e "devtools::document()"`
Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-meta_wrangle.R')"`
Expected: all 4 tests PASS (the `group`/date/trim tests depend on `meta_minimal.xlsx`'s actual columns — if a test fails because a column genuinely isn't present in the fixture, inspect the fixture with `readxl::read_excel()` and adjust the test's column expectations to match, don't weaken the implementation)

- [ ] **Step 7: Commit**

```bash
git add DESCRIPTION CLAUDE.md R/meta_wrangle.R tests/testthat/test-meta_wrangle.R man/meta_read.Rd NAMESPACE
git commit -m "feat: add meta_read() for cleaned experiment metadata import"
```

---

### Task 2: `meta_annotate()`

**Files:**
- Modify: `R/meta_wrangle.R` (append function)
- Modify: `tests/testthat/test-meta_wrangle.R` (append tests)

**Interfaces:**
- Consumes: nothing from Task 1 directly (works on any tibble), but is tested against `meta_read()`'s output.
- Produces: `meta_annotate(data, meta, by = "mouse_ID")` — returns `data` left-joined with `meta` on `by`, invisibly. Errors if `by` missing from either side or if non-`by` column names collide. Warns (does not error) if any `by` value in `data` is unmatched in `meta`.

- [ ] **Step 1: Write the failing tests for `meta_annotate()`**

Append to `tests/testthat/test-meta_wrangle.R`:

```r
test_that("meta_annotate() left-joins data and meta on the by column", {
  data <- tibble::tibble(mouse_ID = c("A", "B"), value = c(1, 2))
  meta <- tibble::tibble(mouse_ID = c("A", "B"), sex = c("m", "f"))
  result <- meta_annotate(data, meta)
  expect_true("sex" %in% names(result))
  expect_equal(result$sex, c("m", "f"))
})

test_that("meta_annotate() warns and keeps NA for unmatched by values", {
  data <- tibble::tibble(mouse_ID = c("A", "B"), value = c(1, 2))
  meta <- tibble::tibble(mouse_ID = "A", sex = "m")
  expect_warning(
    result <- meta_annotate(data, meta),
    regexp = "B",
    fixed = TRUE
  )
  expect_true(is.na(result$sex[result$mouse_ID == "B"]))
})

test_that("meta_annotate() errors when by column missing from data", {
  data <- tibble::tibble(id = "A", value = 1)
  meta <- tibble::tibble(mouse_ID = "A", sex = "m")
  expect_error(meta_annotate(data, meta), regexp = "data", ignore.case = TRUE)
})

test_that("meta_annotate() errors when by column missing from meta", {
  data <- tibble::tibble(mouse_ID = "A", value = 1)
  meta <- tibble::tibble(id = "A", sex = "m")
  expect_error(meta_annotate(data, meta), regexp = "meta", ignore.case = TRUE)
})

test_that("meta_annotate() errors on colliding non-by column names", {
  data <- tibble::tibble(mouse_ID = "A", sex = "unknown")
  meta <- tibble::tibble(mouse_ID = "A", sex = "m")
  expect_error(meta_annotate(data, meta), regexp = "sex", fixed = TRUE)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-meta_wrangle.R')"`
Expected: FAIL with `could not find function "meta_annotate"` (5 new failures)

- [ ] **Step 3: Implement `meta_annotate()`**

Append to `R/meta_wrangle.R`:

```r
#' Annotate experimental data with subject metadata
#'
#' @description
#' Left-joins a metadata tibble (typically from \code{meta_read()}) onto any
#' experimental data tibble (e.g. \code{facs_read_wsp(...)$data}) by a shared
#' identifier column. Errors if the join column is missing from either side,
#' or if any non-join column names collide between the two tibbles. Warns if
#' any \code{by} value present in \code{data} has no match in \code{meta}
#' (those rows keep \code{NA} for all meta columns).
#'
#' @param data tibble to annotate, e.g. \code{facs_read_wsp(...)$data}
#' @param meta metadata tibble, e.g. from \code{meta_read()}
#' @param by name of the shared identifier column, default = "mouse_ID"
#'
#' @returns \code{data} left-joined with \code{meta}, returned invisibly --
#'   assign the result explicitly
#' @export
#'
#' @examples
#' \dontrun{
#'   dat <- meta_annotate(facs_read_wsp("experiment.wsp")$data, meta_read("meta.xlsx"))
#' }
meta_annotate <- function(data, meta, by = "mouse_ID") {
  if (!by %in% names(data)) {
    stop(glue::glue("Join column '{by}' not found in `data`."))
  }
  if (!by %in% names(meta)) {
    stop(glue::glue("Join column '{by}' not found in `meta`."))
  }

  colliding <- setdiff(intersect(names(data), names(meta)), by)
  if (length(colliding) > 0L) {
    stop(glue::glue(
      "Column name(s) present in both `data` and `meta`: ",
      "{paste(colliding, collapse = ', ')}. ",
      "Rename or drop them before calling meta_annotate()."
    ))
  }

  joined <- dplyr::left_join(data, meta, by = by)

  unmatched <- setdiff(unique(data[[by]]), unique(meta[[by]]))
  if (length(unmatched) > 0L) {
    warning(glue::glue(
      "The following '{by}' values in `data` have no match in `meta`: ",
      "{paste(unmatched, collapse = ', ')}"
    ))
  }

  invisible(joined)
}
```

- [ ] **Step 4: Run `devtools::document()` then the tests to verify they pass**

Run: `Rscript -e "devtools::document()"`
Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-meta_wrangle.R')"`
Expected: all tests PASS (9 total: 4 from Task 1, 5 from this task)

- [ ] **Step 5: Commit**

```bash
git add R/meta_wrangle.R tests/testthat/test-meta_wrangle.R man/meta_annotate.Rd NAMESPACE
git commit -m "feat: add meta_annotate() to join subject metadata onto experimental data"
```

---

### Task 3: End-to-end integration test and final check

**Files:**
- Modify: `tests/testthat/test-meta_wrangle.R` (append one test)

**Interfaces:**
- Consumes: `meta_read()` (Task 1), `meta_annotate()` (Task 2), `facs_read_wsp()` (existing, `R/facs_read.R`)

- [ ] **Step 1: Write the failing end-to-end test**

Append to `tests/testthat/test-meta_wrangle.R`:

```r
test_that("facs_read_wsp() data can be annotated end-to-end with meta_read()", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  skip_if_not(file.exists(wsp_path), wsp_skip_msg)

  facs_data <- suppressMessages(
    facs_read_wsp(wsp_path, keywords = "mouse_ID")
  )$data
  meta <- meta_read(meta_path)

  result <- expect_no_warning(meta_annotate(facs_data, meta))

  expect_true(all(c("sex", "group", "cage") %in% names(result)))
  expect_false(any(is.na(result$sex)))
})
```

- [ ] **Step 2: Run the test to verify it fails or passes for the right reason**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-meta_wrangle.R')"`
Expected: PASS if fixtures already line up as designed (both fixture FCS files carry `mouse_ID = "26-1-17"`, present in `meta_minimal.xlsx`). If it fails, the failure message will show which column/value mismatched — inspect actual fixture contents with `readxl::read_excel(meta_path)` and `xml2::read_xml(wsp_path)` rather than adjusting `meta_annotate()`'s behavior to fit.

- [ ] **Step 3: Run the full test suite**

Run: `Rscript -e "devtools::test()"`
Expected: 0 failures across all test files (not just `test-meta_wrangle.R`)

- [ ] **Step 4: Run `devtools::check()`**

Run: `Rscript -e "devtools::check()"`
Expected: 0 errors, 0 warnings. New notes are acceptable only if they are the same class of pre-existing false positives documented in CLAUDE.md's "Known check output" section (tidy-eval / bare column name NOTEs for `meta_read`/`meta_annotate`'s dplyr verbs, e.g. `mouse_ID`, `group`, `sex`). If a genuinely new warning or error appears, fix it before proceeding.

- [ ] **Step 5: Update CLAUDE.md's "Known check output" section if new NOTEs appeared**

If `devtools::check()` produced new NOTEs for bare column names in `meta_wrangle.R` (e.g. `mouse_ID`, `group`, `sex`, `mouse_id`), add a line to the bulleted list under "Known check output" in `CLAUDE.md`, following the existing pattern for `facs_read_wsp`'s variable-binding notes:

```
  - All `meta_read`/`meta_annotate` variable-binding notes (`mouse_id`, `mouse_ID`, `group`) — bare column names inside dplyr verbs; valid tidy eval, flagged as a static-analysis false positive.
```

- [ ] **Step 6: Add the Data Structures section for the new domain to CLAUDE.md**

In `CLAUDE.md`, after the `### report_knit_*()` section (before `## Dependency philosophy`), add:

```
### `meta_read()` / `meta_annotate()` — experiment metadata
- **`meta_read(path, sheet = 1)`**: reads an Excel sheet into a cleaned tibble --
  blank rows/columns dropped, column names snake_cased (except `mouse_ID`,
  preserved verbatim), character columns trimmed, date columns coerced to
  `Date`, `group` coerced to a factor if present.
- **`meta_annotate(data, meta, by = "mouse_ID")`**: left-joins `meta` onto
  `data` by `by`. Errors if `by` is missing from either side or if non-`by`
  column names collide. Warns (row kept, `NA` filled) if a `by` value in
  `data` has no match in `meta`.
```

- [ ] **Step 7: Commit**

```bash
git add tests/testthat/test-meta_wrangle.R CLAUDE.md
git commit -m "test: add end-to-end meta_annotate/facs_read_wsp integration test"
```
