# Multi-sheet meta import + organ-weight/facs-volume wrangling Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let `meta_read()` import every sheet of `meta_minimal.xlsx` (`meta`, `organ_weights`, `facs_volumes`) at once, add a `meta_clean()` function that reshapes and joins `organ_weights`/`facs_volumes` into one tibble, and change `facs_calc_count_per_g()` to process every tissue present in `data` in a single call instead of requiring a pre-filtered, single-tissue `meta`.

**Architecture:** Two files change. `R/meta_wrangle.R`: `meta_annotate()`'s `by` argument generalizes to a character vector (reused later by `meta_clean()`); `meta_read()` drops its `sheet` parameter and returns a named list of cleaned tibbles (one per sheet); a new `meta_clean()` function pivots the wide-by-tissue `organ_weights` sheet to long and joins it with `facs_volumes` on `mouse_ID` + `tissue`. `R/facs_calc.R`: `facs_calc_count_per_g()` drops its `tissue` parameter and joins `meta` on `mouse_ID` + `tissue` instead of `mouse_ID` alone, processing every tissue present in `data` in one call.

**Tech Stack:** R, dplyr, tidyr, purrr, glue, janitor, readxl, testthat 3e — all already in `DESCRIPTION` `Imports`/`Suggests`, no dependency changes.

**Reference spec:** `docs/superpowers/specs/2026-07-02-meta-organ-facs-import-design.md`

## Global Constraints

- Pipe: `|>` (base pipe) everywhere.
- Namespacing: every external function call inside an exported function body must be `pkg::fn()` — no bare `dplyr`/`tidyr`/`purrr`/`glue` calls.
- Deselection uses `dplyr::select(!col)`, never `-col`.
- `meta_read()`, `meta_annotate()`, `meta_clean()` all return their result via `invisible()` (existing convention for this domain).
- Every exported function needs `@param`, `@returns`, `@export`; `@examples` use `\dontrun{}` since all three require realistic multi-column/file input.
- After any roxygen change, run `devtools::document()` before committing.
- Run `devtools::check()` before the final commit — target 0 errors, 0 warnings.
- Non-ASCII characters in R source must be escaped as `\uXXXX`.
- `janitor::clean_names()` lowercases mixed-case column names — e.g. the real fixture's `Treg_vol` column becomes `treg_vol` after `meta_read()`. Any test or example referencing that column must use the lowercase form.

---

## Task 1: `meta_annotate()` — `by` accepts a character vector

**Files:**
- Modify: `R/meta_wrangle.R:51-104` (the `meta_annotate()` roxygen block + function)
- Test: `tests/testthat/test-meta_wrangle.R`

**Interfaces:**
- Produces: `meta_annotate(data, meta, by = "mouse_ID")` — `by` now accepts a character vector of one or more column names (default unchanged). Used by Task 3's `meta_clean()` as `meta_annotate(organ_weights_long, facs_volumes, by = c("mouse_ID", "tissue"))`.

- [ ] **Step 1: Write the failing tests**

Insert the following two tests into `tests/testthat/test-meta_wrangle.R` immediately after the `test_that("meta_annotate() silently drops by values present only in meta", ...)` block (i.e. right before the `test_that("meta_read() skips the mouse_id rename when no such column exists", ...)` block):

```r
test_that("meta_annotate() joins on multiple by columns", {
  data <- tibble::tibble(mouse_ID = c("A", "A"), tissue = c("kidney", "lung"), value = c(1, 2))
  meta <- tibble::tibble(mouse_ID = c("A", "A"), tissue = c("kidney", "lung"), weight = c(100, 50))
  result <- meta_annotate(data, meta, by = c("mouse_ID", "tissue"))
  expect_equal(result$weight, c(100, 50))
})

test_that("meta_annotate() warns listing unmatched mouse_ID/tissue combinations", {
  data <- tibble::tibble(mouse_ID = c("A", "A"), tissue = c("kidney", "lung"), value = c(1, 2))
  meta <- tibble::tibble(mouse_ID = "A", tissue = "kidney", weight = 100)
  expect_warning(
    result <- meta_annotate(data, meta, by = c("mouse_ID", "tissue")),
    regexp = "tissue=lung",
    fixed = TRUE
  )
  expect_true(is.na(result$weight[result$tissue == "lung"]))
})
```

- [ ] **Step 2: Run tests to verify the new tests fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-meta_wrangle.R")'`
Expected: the two new tests FAIL (multi-column `by` isn't supported yet — `left_join` errors or the unmatched-combination logic doesn't run correctly for a vector `by`); all pre-existing tests in the file still PASS.

- [ ] **Step 3: Replace `meta_annotate()`**

Replace lines 51-104 of `R/meta_wrangle.R` (the `meta_annotate()` roxygen block and function) with:

```r
#' Annotate experimental data with subject metadata
#'
#' @description
#' Left-joins a metadata tibble (typically from \code{meta_read()}) onto any
#' experimental data tibble (e.g. \code{facs_read_wsp(...)$data}) by one or
#' more shared identifier columns. Errors if any \code{by} column is missing
#' from either side, or if any non-\code{by} column names collide between
#' the two tibbles. Warns if any \code{by} combination present in
#' \code{data} has no match in \code{meta} (those rows keep \code{NA} for
#' all meta columns). \code{group} is a particularly likely candidate for
#' the column-collision error, since \code{facs_read_wsp(..., keywords =
#' "group")} and \code{meta_read()} both commonly produce a \code{group}
#' column.
#'
#' @param data tibble to annotate, e.g. \code{facs_read_wsp(...)$data}
#' @param meta metadata tibble, e.g. one element of \code{meta_read()}'s
#'   result
#' @param by character vector of shared identifier column name(s), default
#'   = \code{"mouse_ID"}. Pass e.g. \code{c("mouse_ID", "tissue")} to join
#'   on multiple columns.
#'
#' @returns \code{data} left-joined with \code{meta}, returned invisibly --
#'   assign the result explicitly
#' @export
#'
#' @examples
#' \dontrun{
#'   meta_list <- meta_read("meta.xlsx")
#'   dat <- meta_annotate(facs_read_wsp("experiment.wsp")$data, meta_list$meta)
#' }
meta_annotate <- function(data, meta, by = "mouse_ID") {
  missing_in_data <- setdiff(by, names(data))
  if (length(missing_in_data) > 0L) {
    stop(glue::glue(
      "Join column(s) not found in `data`: {paste(missing_in_data, collapse = ', ')}."
    ))
  }
  missing_in_meta <- setdiff(by, names(meta))
  if (length(missing_in_meta) > 0L) {
    stop(glue::glue(
      "Join column(s) not found in `meta`: {paste(missing_in_meta, collapse = ', ')}."
    ))
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

  data_keys <- dplyr::distinct(data, dplyr::across(dplyr::all_of(by)))
  meta_keys <- dplyr::distinct(meta, dplyr::across(dplyr::all_of(by)))
  unmatched <- dplyr::anti_join(data_keys, meta_keys, by = by)

  if (nrow(unmatched) > 0L) {
    unmatched_desc <- purrr::pmap_chr(
      unmatched,
      function(...) {
        row <- list(...)
        paste(paste0(names(row), "=", row), collapse = ", ")
      }
    )
    warning(glue::glue(
      "The following '{paste(by, collapse = ', ')}' combination(s) in `data` have no match in `meta`: ",
      "{paste(unmatched_desc, collapse = '; ')}"
    ))
  }

  invisible(joined)
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-meta_wrangle.R")'`
Expected: PASS, all tests in the file (including the 2 new ones), 0 failures.

- [ ] **Step 5: Document and commit**

```bash
Rscript -e 'devtools::document()'
git add R/meta_wrangle.R tests/testthat/test-meta_wrangle.R man/meta_annotate.Rd
git commit -m "feat: generalize meta_annotate() to accept a multi-column by"
```

---

## Task 2: `meta_read()` — reads every sheet into a named list

**Files:**
- Modify: `R/meta_wrangle.R:1-49` (`is_posixct_()` + `meta_read()` roxygen block + function)
- Test: `tests/testthat/test-meta_wrangle.R`

**Interfaces:**
- Consumes: nothing from Task 1.
- Produces: `meta_read(path)` — exported, `sheet` parameter removed. Returns (invisibly) a named list of tibbles, one per sheet, named after the sheet (e.g. `list(meta = ..., organ_weights = ..., facs_volumes = ...)`). Each tibble is cleaned via a new unexported helper `meta_clean_sheet_(dat)`. Task 3's `meta_clean()` consumes this list shape directly (`meta_list$organ_weights`, `meta_list$facs_volumes`).

- [ ] **Step 1: Rewrite the failing tests**

Replace the entire contents of `tests/testthat/test-meta_wrangle.R` with the following (this replaces all `meta_read()`-specific tests for the new list-of-sheets return; the Task 1 `meta_annotate()` tests are carried over unchanged; the two mocked-binding tests now also mock `excel_sheets`):

```r
library(testthat)
library(JanisHelpers)

meta_path <- testthat::test_path("../fixtures/meta_minimal.xlsx")
wsp_path  <- testthat::test_path("../fixtures/minimal.wsp")
meta_skip_msg <- "meta_minimal.xlsx fixture not available"
wsp_skip_msg  <- "minimal.wsp fixture not available"

test_that("meta_read() returns a named list with one tibble per sheet", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta_list <- meta_read(meta_path)
  expect_type(meta_list, "list")
  expect_setequal(names(meta_list), c("meta", "organ_weights", "facs_volumes"))
  expect_s3_class(meta_list$meta, "tbl_df")
  expect_s3_class(meta_list$organ_weights, "tbl_df")
  expect_s3_class(meta_list$facs_volumes, "tbl_df")
})

test_that("meta_read() preserves mouse_ID verbatim in every sheet", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta_list <- meta_read(meta_path)
  for (sheet in meta_list) {
    expect_true("mouse_ID" %in% names(sheet))
    expect_false("mouse_id" %in% names(sheet))
  }
})

test_that("meta_read() coerces date columns to Date in the meta sheet", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta <- meta_read(meta_path)$meta
  date_cols <- intersect(names(meta), c("dob", "start_date", "bmt_date", "death_date"))
  expect_true(length(date_cols) > 0L)
  for (col in date_cols) {
    expect_s3_class(meta[[col]], "Date")
  }
})

test_that("meta_read() coerces group to a factor in the meta sheet", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta <- meta_read(meta_path)$meta
  expect_true("group" %in% names(meta))
  expect_s3_class(meta$group, "factor")
})

test_that("meta_read() trims whitespace and drops empty rows/cols", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta <- meta_read(meta_path)$meta
  chr_cols <- names(meta)[vapply(meta, is.character, logical(1))]
  for (col in chr_cols) {
    vals <- stats::na.omit(meta[[col]])
    # ignore_attr = TRUE: na.omit() attaches an `na.action` attribute that
    # stringr::str_trim() does not propagate; only the values matter here.
    expect_equal(vals, stringr::str_trim(vals), ignore_attr = TRUE)
  }
  expect_true(all(!vapply(meta, function(x) all(is.na(x)), logical(1))))
})

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

test_that("facs_read_wsp() data can be annotated end-to-end with meta_read()", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  skip_if_not(file.exists(wsp_path), wsp_skip_msg)

  facs_data <- suppressMessages(
    facs_read_wsp(wsp_path, keywords = "mouse_ID")
  )$data
  meta <- meta_read(meta_path)$meta

  result <- expect_no_warning(meta_annotate(facs_data, meta))

  expect_true(all(c("sex", "group", "cage") %in% names(result)))
  expect_false(any(is.na(result$sex)))
})

test_that("meta_annotate() silently drops by values present only in meta", {
  data <- tibble::tibble(mouse_ID = "A", value = 1)
  meta <- tibble::tibble(mouse_ID = c("A", "B"), sex = c("m", "f"))
  result <- expect_no_warning(meta_annotate(data, meta))
  expect_equal(nrow(result), 1L)
  expect_equal(result$mouse_ID, "A")
  expect_false("B" %in% result$mouse_ID)
})

test_that("meta_annotate() joins on multiple by columns", {
  data <- tibble::tibble(mouse_ID = c("A", "A"), tissue = c("kidney", "lung"), value = c(1, 2))
  meta <- tibble::tibble(mouse_ID = c("A", "A"), tissue = c("kidney", "lung"), weight = c(100, 50))
  result <- meta_annotate(data, meta, by = c("mouse_ID", "tissue"))
  expect_equal(result$weight, c(100, 50))
})

test_that("meta_annotate() warns listing unmatched mouse_ID/tissue combinations", {
  data <- tibble::tibble(mouse_ID = c("A", "A"), tissue = c("kidney", "lung"), value = c(1, 2))
  meta <- tibble::tibble(mouse_ID = "A", tissue = "kidney", weight = 100)
  expect_warning(
    result <- meta_annotate(data, meta, by = c("mouse_ID", "tissue")),
    regexp = "tissue=lung",
    fixed = TRUE
  )
  expect_true(is.na(result$weight[result$tissue == "lung"]))
})

test_that("meta_read() skips the mouse_id rename when no such column exists", {
  testthat::local_mocked_bindings(
    excel_sheets = function(...) "sheet1",
    read_excel   = function(...) tibble::tibble(subject = "A", sex = "m"),
    .package = "readxl"
  )
  meta_list <- meta_read("fake/path.xlsx")
  expect_false("mouse_ID" %in% names(meta_list$sheet1))
  expect_false("mouse_id" %in% names(meta_list$sheet1))
  expect_true("subject" %in% names(meta_list$sheet1))
})

test_that("meta_read() skips group factor coercion when no group column exists", {
  testthat::local_mocked_bindings(
    excel_sheets = function(...) "sheet1",
    read_excel   = function(...) tibble::tibble(mouse_id = "A", sex = "m"),
    .package = "readxl"
  )
  meta_list <- meta_read("fake/path.xlsx")
  expect_false("group" %in% names(meta_list$sheet1))
  expect_true("mouse_ID" %in% names(meta_list$sheet1))
})
```

- [ ] **Step 2: Run tests to verify the rewritten tests fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-meta_wrangle.R")'`
Expected: FAIL — `meta_read()` still takes a `sheet` argument and returns a single tibble, not a named list (e.g. `meta_read(meta_path)$meta` will error or return the wrong shape; the mocked tests fail because `excel_sheets` isn't called by the current implementation).

- [ ] **Step 3: Replace `is_posixct_()` and `meta_read()`**

Replace lines 1-49 of `R/meta_wrangle.R` (the `is_posixct_()` helper and the `meta_read()` roxygen block + function) with:

```r
is_posixct_ <- function(x) inherits(x, "POSIXct")

meta_clean_sheet_ <- function(dat) {
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

  dat
}

#' Read and clean every sheet of an experiment metadata spreadsheet
#'
#' @description
#' Reads every sheet of an Excel workbook of experiment metadata and applies
#' standard cleaning to each: blank rows and columns removed, column names
#' standardized to snake_case (except \code{mouse_ID}, preserved verbatim so
#' it stays joinable against FlowJo keyword data such as
#' \code{facs_read_wsp(keywords = "mouse_ID")}), character columns trimmed
#' of whitespace, date/time columns coerced to \code{Date}, and a
#' \code{group} column (if present) coerced to a factor.
#'
#' @param path path to \code{.xlsx} metadata file
#'
#' @returns named list of cleaned tibbles, one per sheet, named after the
#'   sheet (e.g. \code{list(meta = ..., organ_weights = ..., facs_volumes =
#'   ...)}); returned invisibly -- assign the result explicitly
#' @export
#'
#' @examples
#' \dontrun{
#'   meta_list <- meta_read("meta.xlsx")
#'   meta_list$meta
#' }
meta_read <- function(path) {
  sheet_names <- readxl::excel_sheets(path)

  result <- sheet_names |>
    purrr::map(~ meta_clean_sheet_(readxl::read_excel(path, sheet = .x))) |>
    purrr::set_names(sheet_names)

  invisible(result)
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-meta_wrangle.R")'`
Expected: PASS, all tests in the file, 0 failures.

- [ ] **Step 5: Document and commit**

```bash
Rscript -e 'devtools::document()'
git add R/meta_wrangle.R tests/testthat/test-meta_wrangle.R man/meta_read.Rd
git commit -m "feat: meta_read() imports every sheet into a named list"
```

---

## Task 3: `meta_clean()` — pivot `organ_weights` long and join with `facs_volumes`

**Files:**
- Modify: `R/meta_wrangle.R` (append at end of file)
- Test: `tests/testthat/test-meta_wrangle.R` (append at end of file)

**Interfaces:**
- Consumes: `meta_read()`'s list-of-sheets shape (Task 2), `meta_annotate()`'s vector `by` (Task 1).
- Produces: `meta_clean(meta_list)` — exported. `pivot_organ_weights_long_(organ_weights)` — unexported helper. Output shape: one row per `mouse_ID` x `tissue`, columns `total_weight`, `facs_weight` (from `organ_weights`) plus every `facs_volumes` column (`total_vol`, `overview_vol`, `overview_resuspended_vol`, `overview_measured_vol`, `treg_vol`, `treg_resuspended_vol`, `treg_measured_vol` — note lowercase `treg_*`, see Global Constraints). This is the exact shape Task 4's `facs_calc_count_per_g()` expects for its `meta` argument.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-meta_wrangle.R`:

```r
test_that("meta_clean() pivots organ_weights long and joins with facs_volumes", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta_list <- meta_read(meta_path)
  result <- meta_clean(meta_list)

  expect_equal(nrow(result), 9L)
  expect_true(all(c("mouse_ID", "tissue", "total_weight", "facs_weight",
                     "total_vol", "overview_vol", "treg_vol") %in% names(result)))
  expect_setequal(unique(result$tissue), c("kidney", "lung", "spleen"))

  kidney_row <- dplyr::filter(result, mouse_ID == "26-1-1", tissue == "kidney")
  expect_equal(kidney_row$total_weight, 300)
  expect_equal(kidney_row$facs_weight, 200)
})

test_that("meta_clean() errors when organ_weights or facs_volumes is missing", {
  expect_error(
    meta_clean(list(meta = tibble::tibble(mouse_ID = "A"))),
    regexp = "organ_weights",
    fixed = TRUE
  )
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-meta_wrangle.R")'`
Expected: FAIL with `could not find function "meta_clean"`.

- [ ] **Step 3: Append `pivot_organ_weights_long_()` and `meta_clean()`**

Append to `R/meta_wrangle.R`:

```r
pivot_organ_weights_long_ <- function(organ_weights) {
  tidyr::pivot_longer(
    organ_weights,
    cols = !mouse_ID,
    names_to = c("tissue", ".value"),
    names_pattern = "^(.*)_(total_weight|facs_weight)$"
  )
}

#' Combine cleaned organ-weight and FACS-volume sheets into one meta tibble
#'
#' @description
#' Takes the named list returned by \code{meta_read()} and prepares the
#' \code{meta} argument expected by \code{facs_calc_count_per_g()}: pivots
#' the \code{organ_weights} sheet from wide-by-tissue (e.g.
#' \code{kidney_total_weight}, \code{kidney_facs_weight}) to long
#' (\code{mouse_ID}, \code{tissue}, \code{total_weight}, \code{facs_weight}),
#' then left-joins it with the \code{facs_volumes} sheet (already long by
#' \code{tissue}) on \code{mouse_ID} and \code{tissue} via
#' \code{meta_annotate()}.
#'
#' @param meta_list named list as returned by \code{meta_read()}; must
#'   contain \code{organ_weights} and \code{facs_volumes} elements.
#'
#' @returns tibble, one row per \code{mouse_ID} x \code{tissue}, combining
#'   \code{total_weight}/\code{facs_weight} (from \code{organ_weights}) with
#'   the staining-volume columns (from \code{facs_volumes}); returned
#'   invisibly -- assign the result explicitly. Errors if
#'   \code{organ_weights} or \code{facs_volumes} is missing from
#'   \code{meta_list}.
#' @export
#'
#' @examples
#' \dontrun{
#'   meta_combined <- meta_read("meta.xlsx") |> meta_clean()
#' }
meta_clean <- function(meta_list) {
  required <- c("organ_weights", "facs_volumes")
  missing_sheets <- setdiff(required, names(meta_list))
  if (length(missing_sheets) > 0L) {
    stop(glue::glue(
      "`meta_list` is missing required sheet(s): ",
      "{paste(missing_sheets, collapse = ', ')}. Did you read the metadata ",
      "file with meta_read()?"
    ))
  }

  organ_weights_long <- pivot_organ_weights_long_(meta_list$organ_weights)

  meta_annotate(organ_weights_long, meta_list$facs_volumes, by = c("mouse_ID", "tissue"))
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-meta_wrangle.R")'`
Expected: PASS, all tests in the file, 0 failures.

- [ ] **Step 5: Document and commit**

```bash
Rscript -e 'devtools::document()'
git add R/meta_wrangle.R tests/testthat/test-meta_wrangle.R man/meta_clean.Rd
git commit -m "feat: add meta_clean() to pivot organ_weights and join with facs_volumes"
```

---

## Task 4: `facs_calc_count_per_g()` — remove `tissue`, join on `mouse_ID` + `tissue`

**Files:**
- Modify: `R/facs_calc.R:78-213` (the `facs_calc_count_per_g()` roxygen block + function)
- Test: `tests/testthat/test-facs_calc.R`

**Interfaces:**
- Consumes: nothing from Tasks 1-3 directly in its own tests (test fixtures stay inline tibbles, as before); the real fixture end-to-end test at the bottom does not use `meta_clean()` (the fixture's `organ_weights`/`facs_volumes` mouse_IDs — `26-1-1/2/3` — don't match the `minimal.wsp` fixture's mouse_ID `26-1-17`, so it continues to build its `meta` tibble inline, same as today).
- Produces: `facs_calc_count_per_g(data, meta, vol_total, vol_stained, vol_resuspended, vol_measured, organ_piece_weight, method_col = NULL, bead_pop = "beads", bead_concentration = 10400)` — `tissue` parameter removed; `meta` must now include a `tissue` column.

- [ ] **Step 1: Rewrite the failing tests**

Replace the entire contents of `tests/testthat/test-facs_calc.R` with:

```r
library(testthat)
library(JanisHelpers)

wsp_path  <- testthat::test_path("../fixtures/minimal.wsp")
wsp_skip_msg <- "minimal.wsp fixture not available"

test_that("facs_calc_pct_of() computes value / ref_count as a 0-1 fraction", {
  data <- tibble::tibble(
    file_name             = c("f1", "f1", "f1"),
    population_full_path  = c("CD45+", "CD45+/CD3+", "CD45+/CD3+/CD4+"),
    population             = c("CD45+", "CD3+", "CD4+"),
    metric                 = c("count", "count", "count"),
    value                  = c(1000, 400, 100)
  )

  result <- facs_calc_pct_of(data, ref_pop = "CD45+")
  new_rows <- dplyr::filter(result, metric == "pct_of_CD45+") |>
    dplyr::arrange(population)

  expect_equal(new_rows$population, c("CD3+", "CD4+"))
  expect_equal(new_rows$value, c(0.4, 0.1))
  expect_false("CD45+" %in% new_rows$population)
})

test_that("facs_calc_pct_of() errors when ref_pop matches more than one row per file_name", {
  data <- tibble::tibble(
    file_name             = c("f1", "f1"),
    population_full_path  = c("A/Live", "B/Live"),
    population             = c("Live", "Live"),
    metric                 = c("count", "count"),
    value                  = c(100, 200)
  )

  expect_error(facs_calc_pct_of(data, ref_pop = "Live"), regexp = "f1", fixed = TRUE)
})

test_that("facs_calc_pct_of() warns and fills NA when ref_pop has no match for a file_name", {
  data <- tibble::tibble(
    file_name             = c("f1", "f1", "f2", "f2"),
    population_full_path  = c("CD45+", "CD45+/CD3+", "CD3+", "CD3+/CD4+"),
    population             = c("CD45+", "CD3+", "CD3+", "CD4+"),
    metric                 = c("count", "count", "count", "count"),
    value                  = c(1000, 300, 50, 20)
  )

  expect_warning(
    result <- facs_calc_pct_of(data, ref_pop = "CD45+"),
    regexp = "f2",
    fixed = TRUE
  )

  new_rows <- dplyr::filter(result, metric == "pct_of_CD45+")
  f1_row <- dplyr::filter(new_rows, file_name == "f1", population == "CD3+")
  f2_rows <- dplyr::filter(new_rows, file_name == "f2")

  expect_equal(f1_row$value, 0.3)
  expect_true(all(is.na(f2_rows$value)))
})

# ── facs_calc_count_per_g ────────────────────────────────────────────────────

count_per_g_data <- tibble::tibble(
  file_name             = "s1",
  population_full_path  = "CD45+",
  population             = "CD45+",
  metric                 = "count",
  value                  = 1000,
  mouse_ID               = "m1",
  tissue                  = "kidney"
)

count_per_g_meta <- tibble::tibble(
  mouse_ID           = "m1",
  tissue              = "kidney",
  vol_total          = 1000,
  vol_stained        = 100,
  vol_resuspended    = 500,
  vol_measured       = 50,
  organ_piece_weight = 200
)

test_that("facs_calc_count_per_g() HTS formula matches manual calculation (columns)", {
  result <- facs_calc_count_per_g(
    count_per_g_data, count_per_g_meta,
    vol_total = "vol_total", vol_stained = "vol_stained",
    vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
    organ_piece_weight = "organ_piece_weight"
  )
  new_row <- dplyr::filter(result, metric == "count_per_g")
  expect_equal(new_row$value, 500000)
})

test_that("facs_calc_count_per_g() accepts numeric constants in place of column names", {
  meta_no_total <- dplyr::select(count_per_g_meta, !vol_total)

  result <- facs_calc_count_per_g(
    count_per_g_data, meta_no_total,
    vol_total = 1000, vol_stained = "vol_stained",
    vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
    organ_piece_weight = "organ_piece_weight"
  )
  new_row <- dplyr::filter(result, metric == "count_per_g")
  expect_equal(new_row$value, 500000)
})

test_that("facs_calc_count_per_g() computes count_per_g for every tissue present in data", {
  data_multi <- dplyr::bind_rows(
    count_per_g_data,
    tibble::tibble(
      file_name = "s2", population_full_path = "CD45+", population = "CD45+",
      metric = "count", value = 2000, mouse_ID = "m1", tissue = "spleen"
    )
  )
  meta_multi <- dplyr::bind_rows(
    count_per_g_meta,
    tibble::tibble(
      mouse_ID = "m1", tissue = "spleen",
      vol_total = 1000, vol_stained = 100, vol_resuspended = 500,
      vol_measured = 50, organ_piece_weight = 200
    )
  )

  result <- facs_calc_count_per_g(
    data_multi, meta_multi,
    vol_total = "vol_total", vol_stained = "vol_stained",
    vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
    organ_piece_weight = "organ_piece_weight"
  )
  new_rows <- dplyr::filter(result, metric == "count_per_g") |> dplyr::arrange(file_name)
  expect_equal(nrow(new_rows), 2L)
  expect_equal(new_rows$file_name, c("s1", "s2"))
  expect_equal(new_rows$value, c(500000, 1000000))
})

test_that("facs_calc_count_per_g() warns and fills NA when a mouse_ID/tissue combination has no match in meta", {
  meta_wrong_mouse <- dplyr::mutate(count_per_g_meta, mouse_ID = "m2")

  expect_warning(
    result <- facs_calc_count_per_g(
      count_per_g_data, meta_wrong_mouse,
      vol_total = "vol_total", vol_stained = "vol_stained",
      vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
      organ_piece_weight = "organ_piece_weight"
    ),
    regexp = "m1",
    fixed = TRUE
  )
  new_row <- dplyr::filter(result, metric == "count_per_g")
  expect_true(is.na(new_row$value))
})

test_that("facs_calc_count_per_g() errors when a volume argument is neither a column nor a numeric constant", {
  expect_error(
    facs_calc_count_per_g(
      count_per_g_data, count_per_g_meta,
      vol_total = TRUE, vol_stained = "vol_stained",
      vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
      organ_piece_weight = "organ_piece_weight"
    ),
    regexp = "vol_total",
    fixed = TRUE
  )
})

# ── bead path + method_col resolution ────────────────────────────────────────

bead_row <- tibble::tibble(
  file_name = "s1", population_full_path = "beads", population = "beads",
  metric = "count", value = 5200, mouse_ID = "m1", tissue = "kidney"
)

test_that("facs_calc_count_per_g() bead formula matches manual calculation (method_col in meta)", {
  data_beads <- dplyr::bind_rows(count_per_g_data, bead_row)
  meta_beads <- dplyr::mutate(count_per_g_meta, method = "beads")

  result <- facs_calc_count_per_g(
    data_beads, meta_beads,
    vol_total = "vol_total", vol_stained = "vol_stained",
    vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
    organ_piece_weight = "organ_piece_weight",
    method_col = "method"
  )
  new_row <- dplyr::filter(result, metric == "count_per_g", population == "CD45+")
  expect_equal(new_row$value, 100000)
})

test_that("facs_calc_count_per_g() resolves method_col from data (per-sample) over meta", {
  data_beads_kw <- dplyr::bind_rows(count_per_g_data, bead_row) |>
    dplyr::mutate(method = "beads")
  meta_conflicting <- dplyr::mutate(count_per_g_meta, method = "hts")

  result <- facs_calc_count_per_g(
    data_beads_kw, meta_conflicting,
    vol_total = "vol_total", vol_stained = "vol_stained",
    vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
    organ_piece_weight = "organ_piece_weight",
    method_col = "method"
  )
  new_row <- dplyr::filter(result, metric == "count_per_g", population == "CD45+")
  expect_equal(new_row$value, 100000)
})

test_that("facs_calc_count_per_g() defaults NA method_col values to hts with no warning", {
  meta_na_method <- dplyr::mutate(count_per_g_meta, method = NA_character_)

  expect_no_warning(
    result <- facs_calc_count_per_g(
      count_per_g_data, meta_na_method,
      vol_total = "vol_total", vol_stained = "vol_stained",
      vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
      organ_piece_weight = "organ_piece_weight",
      method_col = "method"
    )
  )
  new_row <- dplyr::filter(result, metric == "count_per_g")
  expect_equal(new_row$value, 500000)
})

test_that("facs_calc_count_per_g() errors on an invalid method_col value", {
  meta_bad_method <- dplyr::mutate(count_per_g_meta, method = "unknown")

  expect_error(
    facs_calc_count_per_g(
      count_per_g_data, meta_bad_method,
      vol_total = "vol_total", vol_stained = "vol_stained",
      vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
      organ_piece_weight = "organ_piece_weight",
      method_col = "method"
    ),
    regexp = "unknown",
    fixed = TRUE
  )
})

test_that("facs_calc_count_per_g() errors when method_col is not found in data or meta", {
  expect_error(
    facs_calc_count_per_g(
      count_per_g_data, count_per_g_meta,
      vol_total = "vol_total", vol_stained = "vol_stained",
      vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
      organ_piece_weight = "organ_piece_weight",
      method_col = "nonexistent"
    ),
    regexp = "nonexistent",
    fixed = TRUE
  )
})

test_that("facs_calc_count_per_g() warns and fills NA when bead method resolved but no bead count found", {
  meta_beads_missing <- dplyr::mutate(count_per_g_meta, method = "beads")

  expect_warning(
    result <- facs_calc_count_per_g(
      count_per_g_data, meta_beads_missing,
      vol_total = "vol_total", vol_stained = "vol_stained",
      vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
      organ_piece_weight = "organ_piece_weight",
      method_col = "method"
    ),
    regexp = "s1",
    fixed = TRUE
  )
  new_row <- dplyr::filter(result, metric == "count_per_g")
  expect_true(is.na(new_row$value))
})

# ── end-to-end fixture test ─────────────────────────────────────────────────

test_that("facs_calc_pct_of() and facs_calc_count_per_g() work end-to-end on real fixture data", {
  skip_if_not(file.exists(wsp_path), wsp_skip_msg)

  facs_data <- suppressMessages(
    facs_read_wsp(wsp_path, keywords = c("mouse_ID", "tissue"))
  )$data

  pct_result <- facs_calc_pct_of(facs_data, ref_pop = "Singlets")
  pct_rows <- dplyr::filter(pct_result, metric == "pct_of_Singlets")
  expect_true(nrow(pct_rows) > 0L)
  expect_true(all(pct_rows$value >= 0 & pct_rows$value <= 1, na.rm = TRUE))

  # minimal.wsp has one mouse (26-1-17) with two tissue values: "kidney" and
  # "percoll-kidney" (one FCS file each) -- covering both in `meta` exercises
  # facs_calc_count_per_g() processing multiple tissues in a single call.
  meta <- tibble::tibble(
    mouse_ID           = "26-1-17",
    tissue              = c("kidney", "percoll-kidney"),
    vol_total          = 1000,
    vol_stained        = 100,
    vol_resuspended    = 500,
    vol_measured       = 50,
    organ_piece_weight = 200
  )

  count_result <- expect_no_warning(facs_calc_count_per_g(
    facs_data, meta,
    vol_total = "vol_total", vol_stained = "vol_stained",
    vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
    organ_piece_weight = "organ_piece_weight"
  ))
  count_rows <- dplyr::filter(count_result, metric == "count_per_g")
  expect_equal(dplyr::n_distinct(count_rows$file_name), 2L)
  expect_setequal(
    count_rows$file_name,
    c("26-1-17_whole_kidney_E05.fcs", "26-1-17_percoll-kidney_E06.fcs")
  )
})
```

- [ ] **Step 2: Run tests to verify the rewritten tests fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_calc.R")'`
Expected: FAIL — `facs_calc_count_per_g()` still requires a `tissue` argument (calls omitting it error with a missing-argument message), so most `facs_calc_count_per_g()` tests fail; the 3 `facs_calc_pct_of()` tests still PASS (unaffected).

- [ ] **Step 3: Replace `facs_calc_count_per_g()`**

Replace lines 78-213 of `R/facs_calc.R` (the `facs_calc_count_per_g()` roxygen block and function; leave `facs_calc_pct_of()` and `resolve_var_()` above it untouched) with:

```r
#' Compute absolute cell counts per gram of tissue
#'
#' @description
#' Converts raw event counts into cells per gram of processed tissue, using
#' per-mouse, per-tissue organ weights and staining volumes from
#' \code{meta}. Every \code{mouse_ID} x \code{tissue} combination present in
#' \code{data} is processed in one call (joined against \code{meta} on
#' \code{mouse_ID} and \code{tissue}). Supports two counting methods:
#' HTS/volumetric (the default for every sample unless \code{method_col}
#' says otherwise) and bead-based (using a reference bead population's
#' count and a known bead concentration).
#'
#' @param data tibble shaped like \code{facs_read_wsp(...)$data}, must
#'   include \code{mouse_ID} and \code{tissue} columns (e.g. joined via
#'   \code{facs_read_wsp(keywords = c("mouse_ID", "tissue"))}).
#' @param meta tibble of per-mouse, per-tissue metadata with \code{mouse_ID}
#'   and \code{tissue} columns, e.g. from \code{meta_read()}'s sheets
#'   combined via \code{meta_clean()}.
#' @param vol_total column name in \code{meta} (character) or a single
#'   numeric constant: total organ digest volume.
#' @param vol_stained column name or constant: volume of digest taken for
#'   staining.
#' @param vol_resuspended column name or constant: volume the stained pellet
#'   was resuspended in.
#' @param vol_measured column name or constant: volume actually run/measured
#'   on the cytometer.
#' @param organ_piece_weight column name or constant: weight (mg) of the
#'   organ piece processed.
#' @param method_col character; column name found in \code{data} (checked
#'   first, per-sample) or \code{meta} (checked second, per-mouse) with
#'   values \code{"beads"} or \code{"hts"}. \code{NA} defaults to
#'   \code{"hts"}. \code{NULL} (default) applies \code{"hts"} to every
#'   sample.
#' @param bead_pop character; leaf population name used to look up bead
#'   counts, default \code{"beads"}.
#' @param bead_concentration numeric; reference bead concentration
#'   (beads per microliter), default \code{10400}.
#'
#' @returns \code{data} with additional rows appended:
#'   \code{metric = "count_per_g"}. Errors if \code{method_col} is not found
#'   in \code{data} or \code{meta}, or contains a value outside
#'   \code{{"beads", "hts", NA}}. Warns and fills \code{NA} if a
#'   \code{mouse_ID}/\code{tissue} combination in \code{data} has no match
#'   in \code{meta}, or if the bead method is resolved for a sample with no
#'   matching bead count.
#' @export
#'
#' @examples
#' \dontrun{
#'   dat <- facs_read_wsp("experiment.wsp", keywords = c("mouse_ID", "tissue"))$data
#'   meta_combined <- meta_read("meta.xlsx") |> meta_clean()
#'   facs_calc_count_per_g(
#'     dat, meta_combined,
#'     vol_total = "total_vol", vol_stained = "overview_vol",
#'     vol_resuspended = "overview_resuspended_vol", vol_measured = "overview_measured_vol",
#'     organ_piece_weight = "facs_weight"
#'   )
#' }
facs_calc_count_per_g <- function(
    data,
    meta,
    vol_total,
    vol_stained,
    vol_resuspended,
    vol_measured,
    organ_piece_weight,
    method_col = NULL,
    bead_pop = "beads",
    bead_concentration = 10400
) {
  m <- meta |>
    dplyr::mutate(
      vol_total          = resolve_var_(meta, .env$vol_total, "vol_total"),
      vol_stained        = resolve_var_(meta, .env$vol_stained, "vol_stained"),
      vol_resuspended    = resolve_var_(meta, .env$vol_resuspended, "vol_resuspended"),
      vol_measured       = resolve_var_(meta, .env$vol_measured, "vol_measured"),
      organ_piece_weight = resolve_var_(meta, .env$organ_piece_weight, "organ_piece_weight")
    ) |>
    dplyr::select(mouse_ID, tissue, vol_total, vol_stained, vol_resuspended, vol_measured, organ_piece_weight)

  if (!is.null(method_col) && !method_col %in% names(data) && method_col %in% names(meta)) {
    m[[method_col]] <- meta[[method_col]]
  }
  if (!is.null(method_col) && !method_col %in% names(data) && !method_col %in% names(meta)) {
    stop(glue::glue("`method_col` ('{method_col}') not found in `data` or `meta`."))
  }

  beads <- data |>
    dplyr::filter(population == bead_pop, metric == "count") |>
    dplyr::select(file_name, bead_count = value)

  filtered <- data |>
    dplyr::filter(metric == "count", population != bead_pop) |>
    dplyr::left_join(beads, by = "file_name") |>
    dplyr::left_join(m, by = c("mouse_ID", "tissue"))

  unmatched_combos <- dplyr::anti_join(
    dplyr::distinct(filtered, mouse_ID, tissue),
    dplyr::distinct(m, mouse_ID, tissue),
    by = c("mouse_ID", "tissue")
  )
  if (nrow(unmatched_combos) > 0L) {
    unmatched_desc <- purrr::pmap_chr(
      unmatched_combos,
      function(mouse_ID, tissue) glue::glue("mouse_ID={mouse_ID}, tissue={tissue}")
    )
    warning(glue::glue(
      "The following mouse_ID/tissue combination(s) in `data` have no match in `meta`: ",
      "{paste(unmatched_desc, collapse = '; ')}"
    ))
  }

  filtered$method <- if (is.null(method_col)) {
    "hts"
  } else {
    dplyr::coalesce(filtered[[method_col]], "hts")
  }

  bad_methods <- unique(filtered$method[!filtered$method %in% c("beads", "hts")])
  if (length(bad_methods) > 0L) {
    stop(glue::glue(
      "`method_col` contains value(s) other than 'beads'/'hts': ",
      "{paste(bad_methods, collapse = ', ')}"
    ))
  }

  missing_bead <- filtered$method == "beads" & is.na(filtered$bead_count)
  if (any(missing_bead)) {
    warning(glue::glue(
      "Bead method resolved but no bead count found for file_name(s): ",
      "{paste(unique(filtered$file_name[missing_bead]), collapse = ', ')}. Result filled with NA."
    ))
  }

  new_rows <- filtered |>
    dplyr::mutate(
      metric = "count_per_g",
      value  = dplyr::if_else(
        method == "beads",
        ((value / (bead_count / bead_concentration)) / (vol_stained / vol_total)) / (organ_piece_weight / 1000),
        ((value / (vol_measured / vol_resuspended)) / (vol_stained / vol_total)) / (organ_piece_weight / 1000)
      )
    ) |>
    dplyr::select(dplyr::all_of(names(data)))

  dplyr::bind_rows(data, new_rows) |>
    dplyr::arrange(file_name, population_full_path)
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_calc.R")'`
Expected: PASS, all tests in the file, 0 failures.

- [ ] **Step 5: Document and commit**

```bash
Rscript -e 'devtools::document()'
git add R/facs_calc.R tests/testthat/test-facs_calc.R man/facs_calc_count_per_g.Rd
git commit -m "feat: facs_calc_count_per_g() processes every tissue in one call"
```

---

## Task 5: CLAUDE.md housekeeping, full check, final commit

**Files:**
- Modify: `CLAUDE.md`

**Interfaces:**
- Consumes: all of Tasks 1-4 (final verification pass across the whole package).

- [ ] **Step 1: Update the `meta_read()`/`meta_annotate()` data-structures section**

In `CLAUDE.md`, find this block:

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

Replace it with:

```
### `meta_read()` / `meta_clean()` / `meta_annotate()` — experiment metadata
- **`meta_read(path)`**: reads every sheet of an Excel workbook into a named
  list of cleaned tibbles (one per sheet, named after the sheet) -- blank
  rows/columns dropped, column names snake_cased (except `mouse_ID`,
  preserved verbatim), character columns trimmed, date columns coerced to
  `Date`, `group` coerced to a factor if present (each rule applied
  independently per sheet, only where applicable columns exist).
- **`meta_clean(meta_list)`**: takes the list from `meta_read()`, pivots its
  `organ_weights` element from wide-by-tissue (e.g. `kidney_total_weight`,
  `kidney_facs_weight`) to long (`mouse_ID`, `tissue`, `total_weight`,
  `facs_weight`), then joins it with `facs_volumes` (already long by
  `tissue`) via `meta_annotate(by = c("mouse_ID", "tissue"))`. Returns one
  tibble ready to pass as `meta` to `facs_calc_count_per_g()`. Errors if
  `organ_weights` or `facs_volumes` is missing from `meta_list`.
- **`meta_annotate(data, meta, by = "mouse_ID")`**: left-joins `meta` onto
  `data` by `by`, a character vector of one or more shared column names
  (e.g. `c("mouse_ID", "tissue")`). Errors if any `by` column is missing
  from either side or if non-`by` column names collide. Warns (rows kept,
  `NA` filled) if a `by` combination in `data` has no match in `meta`.
```

- [ ] **Step 2: Run the full test suite**

Run: `Rscript -e 'devtools::test()'`
Expected: all tests pass, including every test in `test-meta_wrangle.R` and `test-facs_calc.R`.

- [ ] **Step 3: Run `devtools::check()` and compare NOTEs against CLAUDE.md**

Run: `Rscript -e 'devtools::check()'`
Expected: 0 errors, 0 warnings. Compare the "no visible binding for global variable" NOTEs against `CLAUDE.md`'s existing "Known check output" bullet for `facs_calc_pct_of`/`facs_calc_count_per_g`. Since `tissue` is no longer a formal argument of `facs_calc_count_per_g()`, it should now appear as a flagged bare-column-name NOTE (it wasn't flagged before, when it was a parameter). If any other new variable names appear in the NOTE output that aren't covered by the update in Step 4 below, add them to the same bullet using the same phrasing style.

- [ ] **Step 4: Update the "Known check output" section**

In `CLAUDE.md`, find this line:

```
  - All `facs_calc_pct_of`/`facs_calc_count_per_g` variable-binding notes (`population`, `metric`, `value`, `file_name`, `population_full_path`, `ref_count`, `mouse_ID`, `method`, `bead_count`, `.env`) — bare column names inside dplyr verbs; valid tidy eval, flagged as a static-analysis false positive. Note: the function's own formal arguments (`tissue`, `vol_total`, `vol_stained`, `vol_resuspended`, `vol_measured`, `organ_piece_weight`, `method_col`, `bead_pop`, `bead_concentration`) are also referenced bare inside `dplyr` verbs but are not flagged, since codetools resolves them against the matching parameter name in scope.
```

Replace it with:

```
  - All `facs_calc_pct_of`/`facs_calc_count_per_g` variable-binding notes (`population`, `metric`, `value`, `file_name`, `population_full_path`, `ref_count`, `mouse_ID`, `tissue`, `method`, `bead_count`, `.env`) — bare column names inside dplyr verbs; valid tidy eval, flagged as a static-analysis false positive. `tissue` is a bare column reference (not a formal argument, unlike in the old single-tissue signature) — `facs_calc_count_per_g()` now joins `meta` on `mouse_ID` + `tissue` and processes every tissue present in `data` in one call. Note: the function's own formal arguments (`vol_total`, `vol_stained`, `vol_resuspended`, `vol_measured`, `organ_piece_weight`, `method_col`, `bead_pop`, `bead_concentration`) are also referenced bare inside `dplyr` verbs but are not flagged, since codetools resolves them against the matching parameter name in scope.
  - All `meta_clean`/`pivot_organ_weights_long_` variable-binding notes (`mouse_ID`) — bare column name inside a `tidyr::pivot_longer()` deselection (`cols = !mouse_ID`); valid tidy eval, flagged as a static-analysis false positive, same pattern as the existing `meta_read`/`meta_annotate` note above.
```

If Step 3's actual `devtools::check()` output differs from this prediction (e.g. a NOTE variable name isn't listed here, or one listed here doesn't actually appear), adjust the bullet text to match the real output rather than leaving it inaccurate.

- [ ] **Step 5: Final commit**

```bash
git add CLAUDE.md
git commit -m "docs: update CLAUDE.md for meta_read()/meta_clean()/facs_calc_count_per_g() changes"
```

---

## Self-Review Notes

- **Spec coverage:** Section 1 (`meta_clean()` + `pivot_organ_weights_long_()`) — Task 3. Section 2 (`meta_read()` returns all sheets) — Task 2. Section 3 (`meta_annotate()` vector `by`) — Task 1. Section 4 (`facs_calc_count_per_g()` drops `tissue`, joins on `mouse_ID`+`tissue`) — Task 4. End-to-end example from the spec — exercised by Task 4's rewritten end-to-end test (via the real `minimal.wsp` fixture's two tissues) and Task 3's `meta_clean()` fixture test (via the real `meta_minimal.xlsx` fixture). CLAUDE.md housekeeping — Task 5. All spec sections covered.
- **Type consistency:** `meta_read()`'s new list-of-tibbles return type, `meta_clean()`'s `meta_list` parameter, `meta_annotate()`'s vector `by`, and `facs_calc_count_per_g()`'s `meta`-with-`tissue`-column requirement are used identically across Tasks 2-5 and in every roxygen example.
- **No placeholders:** every step has complete, runnable code or an exact command with expected output. Task 5's CLAUDE.md NOTE prediction includes an explicit fallback instruction (adjust to real `devtools::check()` output) rather than asserting a guess as fact.
- **Fixture-shape verification performed during planning:** confirmed via `Rscript` that `janitor::clean_names()` lowercases `Treg_vol` → `treg_vol` (documented in Global Constraints and used correctly in Task 3's test), and that `minimal.wsp` contains one mouse (`26-1-17`) with two tissue values (`kidney`, `percoll-kidney`), each a separate FCS file — used to build a genuine multi-tissue end-to-end test in Task 4 instead of a synthetic one.
