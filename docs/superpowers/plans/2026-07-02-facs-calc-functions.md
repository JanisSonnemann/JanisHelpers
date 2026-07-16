# facs_calc_pct_of() + facs_calc_count_per_g() Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the abandoned `R/new_facs_functions.R` draft with two properly namespaced, tested, roxygen-documented functions in `R/facs_calc.R`: `facs_calc_pct_of()` (percentage of an arbitrary ancestor population) and `facs_calc_count_per_g()` (absolute cell counts per gram of tissue, HTS or bead-based).

**Architecture:** Single new file `R/facs_calc.R` with one unexported helper (`resolve_var_()`) and two exported functions, both operating on `facs_read_wsp()`-shaped tibbles (snake_case columns, lowercase `metric` values). `R/new_facs_functions.R` is deleted.

**Tech Stack:** R, dplyr, tidyr, glue, testthat 3e — all already in `DESCRIPTION` `Imports`/`Suggests`, no dependency changes.

## Global Constraints

- Pipe: `|>` (base pipe) everywhere.
- Namespacing: every external function call inside an exported function body must be `pkg::fn()` — no bare `dplyr`/`tidyr`/`glue` calls.
- Deselection uses `dplyr::select(!col)`, never `-col`.
- Functions returning data frames primarily for use (not side-effect) — both these do — return visibly (do not wrap in `invisible()`).
- File/function naming: `domain_verb.R` / `domain_verb_qualifier()`, domain prefix `facs_`.
- Every exported function needs `@param`, `@returns`, `@export`; `@examples` use `\dontrun{}` since both require realistic multi-column input.
- After any roxygen change, run `devtools::document()` before committing.
- Run `devtools::check()` before the final commit — target 0 errors, 0 warnings.
- Non-ASCII characters in R source must be escaped as `\uXXXX`.

---

## Task 1: `facs_calc_pct_of()` + delete the dead draft

**Files:**
- Create: `R/facs_calc.R`
- Delete: `R/new_facs_functions.R`
- Test: `tests/testthat/test-facs_calc.R`

**Interfaces:**
- Produces: `facs_calc_pct_of(data, ref_pop)` — exported. `data` is `facs_read_wsp(...)$data`-shaped (`file_name`, `population_full_path`, `population`, `metric`, `value`, plus any keyword columns). `ref_pop` is a single character string matching a `population` leaf name. Returns `data` with additional rows appended (visibly).

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-facs_calc.R`:

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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_calc.R")'`
Expected: FAIL with `could not find function "facs_calc_pct_of"`.

- [ ] **Step 3: Delete the dead draft**

```bash
git rm R/new_facs_functions.R
```

- [ ] **Step 4: Create `R/facs_calc.R` with `facs_calc_pct_of()`**

```r
#' Compute a population's percentage of an arbitrary reference population
#'
#' @description
#' Computes each population's count as a fraction of a named ancestor
#' population's count, regardless of how many gating levels separate them
#' (unlike \code{fraction_of_parent} from \code{facs_read_wsp()}, which is
#' always relative to the immediate parent gate).
#'
#' @param data tibble shaped like \code{facs_read_wsp(...)$data}: must
#'   contain \code{file_name}, \code{population_full_path}, \code{population},
#'   \code{metric}, \code{value}.
#' @param ref_pop character; leaf population name (matches \code{population},
#'   not \code{population_full_path}) to use as the denominator.
#'
#' @returns \code{data} with additional rows appended: one row per
#'   \code{file_name x population} (excluding \code{ref_pop} itself), with
#'   \code{metric = paste0("pct_of_", ref_pop)} and \code{value} as a 0-1
#'   fraction of \code{ref_pop}'s count in that file. Errors if \code{ref_pop}
#'   matches more than one population per file; warns and fills \code{NA} if
#'   \code{ref_pop} has no match for a file.
#' @export
#'
#' @examples
#' \dontrun{
#'   dat <- facs_read_wsp("experiment.wsp")$data
#'   facs_calc_pct_of(dat, ref_pop = "CD45+")
#' }
facs_calc_pct_of <- function(data, ref_pop) {
  ref_counts <- data |>
    dplyr::filter(population == ref_pop, metric == "count") |>
    dplyr::select(file_name, ref_count = value)

  dup_files <- unique(ref_counts$file_name[duplicated(ref_counts$file_name)])
  if (length(dup_files) > 0L) {
    stop(glue::glue(
      "ref_pop '{ref_pop}' matches more than one population for file_name(s): ",
      "{paste(dup_files, collapse = ', ')}. Leaf population name is ambiguous ",
      "— rename the gate(s) in FlowJo or disambiguate before calling facs_calc_pct_of()."
    ))
  }

  missing_files <- setdiff(unique(data$file_name), ref_counts$file_name)
  if (length(missing_files) > 0L) {
    warning(glue::glue(
      "ref_pop '{ref_pop}' not found for file_name(s): ",
      "{paste(missing_files, collapse = ', ')}. Result filled with NA."
    ))
  }

  new_rows <- data |>
    dplyr::filter(metric == "count", population != ref_pop) |>
    dplyr::left_join(ref_counts, by = "file_name") |>
    dplyr::mutate(
      metric = paste0("pct_of_", ref_pop),
      value  = value / ref_count
    ) |>
    dplyr::select(dplyr::all_of(names(data)))

  dplyr::bind_rows(data, new_rows) |>
    dplyr::arrange(file_name, population_full_path)
}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_calc.R")'`
Expected: PASS, 3 tests, 0 failures.

- [ ] **Step 6: Document and commit**

```bash
Rscript -e 'devtools::document()'
git add R/facs_calc.R tests/testthat/test-facs_calc.R man/facs_calc_pct_of.Rd NAMESPACE
git commit -m "feat: add facs_calc_pct_of(); remove abandoned new_facs_functions.R draft"
```

---

## Task 2: `facs_calc_count_per_g()` — HTS path

**Files:**
- Modify: `R/facs_calc.R`
- Test: `tests/testthat/test-facs_calc.R`

**Interfaces:**
- Consumes: nothing from Task 1's `facs_calc_pct_of()` — independent function in the same file.
- Produces: `resolve_var_(df, var, arg_name)` (unexported) — returns `df[[var]]` if `var` is a length-1 character naming a column of `df`, `rep(var, nrow(df))` if `var` is a length-1 numeric, else `stop()`. `facs_calc_count_per_g(data, meta, tissue, vol_total, vol_stained, vol_resuspended, vol_measured, organ_piece_weight, method_col = NULL, bead_pop = "beads", bead_concentration = 10400)` — exported (bead-method behavior added in Task 3; this task implements the HTS-only path with `method_col` accepted but not yet wired to bead logic beyond always resolving to `"hts"` when `NULL`).

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-facs_calc.R`:

```r
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
  vol_total          = 1000,
  vol_stained        = 100,
  vol_resuspended    = 500,
  vol_measured       = 50,
  organ_piece_weight = 200
)

test_that("facs_calc_count_per_g() HTS formula matches manual calculation (columns)", {
  result <- facs_calc_count_per_g(
    count_per_g_data, count_per_g_meta, tissue = "kidney",
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
    count_per_g_data, meta_no_total, tissue = "kidney",
    vol_total = 1000, vol_stained = "vol_stained",
    vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
    organ_piece_weight = "organ_piece_weight"
  )
  new_row <- dplyr::filter(result, metric == "count_per_g")
  expect_equal(new_row$value, 500000)
})

test_that("facs_calc_count_per_g() only processes rows matching the tissue argument", {
  data_multi <- dplyr::bind_rows(
    count_per_g_data,
    tibble::tibble(
      file_name = "s2", population_full_path = "CD45+", population = "CD45+",
      metric = "count", value = 2000, mouse_ID = "m1", tissue = "spleen"
    )
  )

  result <- facs_calc_count_per_g(
    data_multi, count_per_g_meta, tissue = "kidney",
    vol_total = "vol_total", vol_stained = "vol_stained",
    vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
    organ_piece_weight = "organ_piece_weight"
  )
  new_rows <- dplyr::filter(result, metric == "count_per_g")
  expect_equal(nrow(new_rows), 1L)
  expect_equal(new_rows$file_name, "s1")
})

test_that("facs_calc_count_per_g() errors when a volume argument is neither a column nor a numeric constant", {
  expect_error(
    facs_calc_count_per_g(
      count_per_g_data, count_per_g_meta, tissue = "kidney",
      vol_total = TRUE, vol_stained = "vol_stained",
      vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
      organ_piece_weight = "organ_piece_weight"
    ),
    regexp = "vol_total",
    fixed = TRUE
  )
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_calc.R")'`
Expected: FAIL with `could not find function "facs_calc_count_per_g"`.

- [ ] **Step 3: Implement `resolve_var_()` and `facs_calc_count_per_g()`**

Append to `R/facs_calc.R`:

```r
resolve_var_ <- function(df, var, arg_name) {
  if (is.character(var) && length(var) == 1L) {
    if (!var %in% names(df)) {
      stop(glue::glue("Column '{var}' (from `{arg_name}`) not found in `meta`."))
    }
    df[[var]]
  } else if (is.numeric(var) && length(var) == 1L) {
    rep(var, nrow(df))
  } else {
    stop(glue::glue(
      "`{arg_name}` must be a single column name (character) or a single numeric value."
    ))
  }
}

#' Compute absolute cell counts per gram of tissue
#'
#' @description
#' Converts raw event counts into cells per gram of processed tissue, using
#' per-mouse organ weights and staining volumes from \code{meta}. Supports
#' two counting methods: HTS/volumetric (the default for every sample unless
#' \code{method_col} says otherwise) and bead-based (using a reference bead
#' population's count and a known bead concentration).
#'
#' @param data tibble shaped like \code{facs_read_wsp(...)$data}, must
#'   include \code{mouse_ID} and \code{tissue} columns (e.g. joined via
#'   \code{facs_read_wsp(keywords = c("mouse_ID", "tissue"))}).
#' @param meta tibble of per-mouse metadata keyed by \code{mouse_ID}, e.g.
#'   from \code{meta_read()} merged with organ-weight/volume columns.
#' @param tissue character; value to filter \code{data$tissue} on (e.g.
#'   \code{"kidney"}).
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
#'   (beads/µL), default \code{10400}.
#'
#' @returns \code{data} with additional rows appended:
#'   \code{metric = "count_per_g"}. Errors if \code{method_col} is not found
#'   in \code{data} or \code{meta}, or contains a value outside
#'   \code{{"beads", "hts", NA}}. Warns and fills \code{NA} if the bead
#'   method is resolved for a sample with no matching bead count.
#' @export
#'
#' @examples
#' \dontrun{
#'   dat <- facs_read_wsp("experiment.wsp", keywords = c("mouse_ID", "tissue"))$data
#'   facs_calc_count_per_g(
#'     dat, meta_read("meta.xlsx"), tissue = "kidney",
#'     vol_total = "kidney_vol_total", vol_stained = "kidney_vol_stained",
#'     vol_resuspended = "kidney_vol_resuspended", vol_measured = "kidney_vol_measured",
#'     organ_piece_weight = "kidney_piece_weight"
#'   )
#' }
facs_calc_count_per_g <- function(
    data,
    meta,
    tissue,
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
      vol_total          = resolve_var_(meta, vol_total, "vol_total"),
      vol_stained        = resolve_var_(meta, vol_stained, "vol_stained"),
      vol_resuspended    = resolve_var_(meta, vol_resuspended, "vol_resuspended"),
      vol_measured       = resolve_var_(meta, vol_measured, "vol_measured"),
      organ_piece_weight = resolve_var_(meta, organ_piece_weight, "organ_piece_weight")
    ) |>
    dplyr::select(mouse_ID, vol_total, vol_stained, vol_resuspended, vol_measured, organ_piece_weight)

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
    dplyr::filter(metric == "count", tissue == .env$tissue) |>
    dplyr::left_join(beads, by = "file_name") |>
    dplyr::left_join(m, by = "mouse_ID")

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
Expected: PASS, 7 tests total, 0 failures.

- [ ] **Step 5: Document and commit**

```bash
Rscript -e 'devtools::document()'
git add R/facs_calc.R tests/testthat/test-facs_calc.R man/facs_calc_count_per_g.Rd NAMESPACE
git commit -m "feat: add facs_calc_count_per_g() HTS path"
```

---

## Task 3: `facs_calc_count_per_g()` — bead path + `method_col` resolution

**Files:**
- Modify: none (implementation already complete from Task 2 — this task is test-only, verifying the bead branch, `method_col` source resolution, and error/warning paths already implemented)
- Test: `tests/testthat/test-facs_calc.R`

**Interfaces:**
- Consumes: `facs_calc_count_per_g()` from Task 2 (already implements the full bead/`method_col` logic — Task 2's implementation was written ahead to keep the function whole and testable in one piece; this task's job is proving every branch behaves as specified).

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-facs_calc.R`:

```r
bead_row <- tibble::tibble(
  file_name = "s1", population_full_path = "beads", population = "beads",
  metric = "count", value = 5200, mouse_ID = "m1", tissue = "kidney"
)

test_that("facs_calc_count_per_g() bead formula matches manual calculation (method_col in meta)", {
  data_beads <- dplyr::bind_rows(count_per_g_data, bead_row)
  meta_beads <- dplyr::mutate(count_per_g_meta, method = "beads")

  result <- facs_calc_count_per_g(
    data_beads, meta_beads, tissue = "kidney",
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

  result <- facs_calc_count_per_g(
    data_beads_kw, count_per_g_meta, tissue = "kidney",
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
      count_per_g_data, meta_na_method, tissue = "kidney",
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
      count_per_g_data, meta_bad_method, tissue = "kidney",
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
      count_per_g_data, count_per_g_meta, tissue = "kidney",
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
      count_per_g_data, meta_beads_missing, tissue = "kidney",
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
```

- [ ] **Step 2: Run tests to verify they fail or pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_calc.R")'`
Expected: PASS, 13 tests total, 0 failures — Task 2's implementation already covers these branches. If any of these 6 new tests fail, fix `facs_calc_count_per_g()` in `R/facs_calc.R` (do not weaken the test) before proceeding.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-facs_calc.R
git commit -m "test: cover bead formula and method_col resolution in facs_calc_count_per_g()"
```

---

## Task 4: End-to-end fixture test, CLAUDE.md housekeeping, final check

**Files:**
- Modify: `tests/testthat/test-facs_calc.R`, `CLAUDE.md`

**Interfaces:**
- Consumes: `facs_calc_pct_of()`, `facs_calc_count_per_g()` (Tasks 1–3), `facs_read_wsp()` (existing, `R/facs_read.R`).

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-facs_calc.R`:

```r
test_that("facs_calc_pct_of() and facs_calc_count_per_g() work end-to-end on real fixture data", {
  skip_if_not(file.exists(wsp_path), wsp_skip_msg)

  facs_data <- suppressMessages(
    facs_read_wsp(wsp_path, keywords = c("mouse_ID", "tissue"))
  )$data

  pct_result <- facs_calc_pct_of(facs_data, ref_pop = "Singlets")
  pct_rows <- dplyr::filter(pct_result, metric == "pct_of_Singlets")
  expect_true(nrow(pct_rows) > 0L)
  expect_true(all(pct_rows$value >= 0 & pct_rows$value <= 1, na.rm = TRUE))

  meta <- tibble::tibble(
    mouse_ID           = "26-1-17",
    vol_total          = 1000,
    vol_stained        = 100,
    vol_resuspended    = 500,
    vol_measured       = 50,
    organ_piece_weight = 200
  )

  count_result <- facs_calc_count_per_g(
    facs_data, meta, tissue = "kidney",
    vol_total = "vol_total", vol_stained = "vol_stained",
    vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
    organ_piece_weight = "organ_piece_weight"
  )
  count_rows <- dplyr::filter(count_result, metric == "count_per_g")
  expect_equal(dplyr::n_distinct(count_rows$file_name), 1L)
  expect_equal(unique(count_rows$file_name), "26-1-17_whole_kidney_E05.fcs")
})
```

- [ ] **Step 2: Run test to verify it fails or passes**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_calc.R")'`
Expected: PASS, 14 tests total, 0 failures (this test exercises only already-implemented code against real fixture data — it should pass immediately; if `mouse_ID`/`tissue` values in `tests/fixtures/minimal.wsp` don't match `"26-1-17"`/`"kidney"`, adjust the test's `meta$mouse_ID` and `tissue =` argument to match the actual fixture keyword values, verified via `facs_read_wsp(wsp_path, keywords = c("mouse_ID", "tissue"))$data`).

- [ ] **Step 3: Run `devtools::check()` and record the actual NOTEs**

Run: `Rscript -e 'devtools::check()'`
Expected: 0 errors, 0 warnings. Note any new "no visible binding for global variable" NOTEs for `R/facs_calc.R` — these are tidy-eval false positives (bare column names inside `dplyr` verbs: `population`, `metric`, `value`, `file_name`, `ref_count`, `mouse_ID`, `tissue`, `vol_total`, `vol_stained`, `vol_resuspended`, `vol_measured`, `organ_piece_weight`, `bead_count`, `method`), matching the existing pattern already documented for `facs_read_wsp`/`meta_annotate` in CLAUDE.md.

- [ ] **Step 4: Update CLAUDE.md's "Known check output" section**

In `CLAUDE.md`, find this line (near the end of the "Known check output" section):

```
  - All `meta_read`/`meta_annotate` variable-binding notes (`mouse_id`, `mouse_ID`, `group`) — bare column names inside dplyr verbs; valid tidy eval, flagged as a static-analysis false positive.
```

Add immediately after it:

```
  - All `facs_calc_pct_of`/`facs_calc_count_per_g` variable-binding notes (`population`, `metric`, `value`, `file_name`, `ref_count`, `mouse_ID`, `tissue`, `vol_total`, `vol_stained`, `vol_resuspended`, `vol_measured`, `organ_piece_weight`, `bead_count`, `method`) — bare column names inside dplyr verbs; valid tidy eval, flagged as a static-analysis false positive.
```

- [ ] **Step 5: Final full test suite + check run**

Run: `Rscript -e 'devtools::test()'`
Expected: All tests pass, including the 14 in `test-facs_calc.R` and all pre-existing tests in the other `test-*.R` files.

Run: `Rscript -e 'devtools::check()'`
Expected: 0 errors, 0 warnings (NOTEs limited to the pre-existing ones already listed in CLAUDE.md, plus the new `facs_calc_*` variable-binding NOTE now documented in Step 4).

- [ ] **Step 6: Commit**

```bash
git add tests/testthat/test-facs_calc.R CLAUDE.md
git commit -m "test: add end-to-end facs_calc test against real fixture; docs: note facs_calc variable-binding NOTEs"
```

---

## Self-Review Notes

- **Spec coverage:** `facs_calc_pct_of()` (ambiguity error, missing-match warning, 0-1 fraction) — Task 1. `facs_calc_count_per_g()` HTS path, `resolve_var_()` column-or-constant, `tissue` filtering — Task 2. Bead formula, `method_col` resolution order (data then meta), `NA`-defaults-to-hts, invalid-value error, not-found error, missing-bead-count warning — Task 3. End-to-end fixture test, CLAUDE.md housekeeping — Task 4. All spec sections covered.
- **Type consistency:** `facs_calc_count_per_g()`'s signature, `resolve_var_()`'s signature, and every call site across Tasks 2–4 use identical argument names and order throughout.
- **No placeholders:** every step has complete, runnable code or an exact command with expected output.
