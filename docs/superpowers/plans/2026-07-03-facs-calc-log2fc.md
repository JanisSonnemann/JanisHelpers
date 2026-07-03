# facs_calc_log2fc() Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `facs_calc_log2fc()` to `R/facs_calc.R` — a T-cell restimulation normalizer that computes, per mouse/tissue/population, the log2 fold-change of each stimulated condition's `ref_pop`-normalized proportion against an unstimulated reference condition. Replaces a manually-written draft that hardcoded stim condition names, hardcoded a `group` passthrough column, and had a pseudocount default (`0`) that contradicted its own comment.

**Architecture:** One new exported function appended to the existing `R/facs_calc.R` (alongside `facs_calc_pct_of()`/`facs_calc_count_per_g()`), no new unexported helpers. Unlike its siblings, it does **not** append rows to `data` and return it — its output grain (one row per `mouse_ID` x `tissue` x `population` x non-reference restimulation level) collapses the restimulation dimension via comparison, so it returns a new, differently-shaped tibble.

**Tech Stack:** R, dplyr, tidyr, glue, purrr, testthat 3e — all already in `DESCRIPTION` `Imports`/`Suggests`, no dependency changes.

## Global Constraints

- Pipe: `|>` (base pipe) everywhere.
- Namespacing: every external function call inside the exported function body must be `pkg::fn()` — no bare `dplyr`/`tidyr`/`glue`/`purrr` calls. Avoid `\(x)` backslash-lambda syntax (not used elsewhere in this codebase) — use `function(x) ...` instead.
- Deselection uses `dplyr::select(!col)`, never `-col`.
- `facs_calc_log2fc()` returns a **new** tibble (not `data` with rows appended, unlike `facs_calc_pct_of()`/`facs_calc_count_per_g()`) — do not wrap in `dplyr::bind_rows(data, ...)`.
- The function returns visibly (primary output, not a side effect) — do not wrap in `invisible()`.
- File/function naming: `domain_verb.R` / `domain_verb_qualifier()`, domain prefix `facs_`.
- The exported function needs `@param`, `@returns`, `@export`; `@examples` uses `\dontrun{}` since it requires realistic multi-column input.
- After any roxygen change, run `devtools::document()` before committing.
- Run `devtools::check()` before the final commit — target 0 errors, 0 warnings.
- Non-ASCII characters in R source must be escaped as `\uXXXX`.

---

## Task 1: `facs_calc_log2fc()` — core implementation + happy-path/error tests

**Files:**
- Modify: `R/facs_calc.R`
- Test: `tests/testthat/test-facs_calc.R`

**Interfaces:**
- Produces: `facs_calc_log2fc(data, ref_pop, restim_col = "restimulation", ref_level = "unstim", pseudocount = 0.5)` — exported. `data` is `facs_read_wsp(...)$data`-shaped, must include `mouse_ID`, `tissue`, `population`, `metric`, `value`, and a `restim_col` column. Returns a new tibble: `mouse_ID`, `tissue`, `population`, `<passthrough cols>`, `restim_col` (name reused), `log2fc`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-facs_calc.R`:

```r
# ── facs_calc_log2fc ─────────────────────────────────────────────────────────

log2fc_data <- tibble::tibble(
  file_name             = c("f1", "f1", "f2", "f2"),
  population_full_path  = c("CD3+", "CD3+/CD4+", "CD3+", "CD3+/CD4+"),
  population             = c("CD3+", "CD4+", "CD3+", "CD4+"),
  metric                 = "count",
  value                  = c(1000, 100, 1000, 400),
  mouse_ID               = "m1",
  tissue                  = "spleen",
  restimulation           = c("unstim", "unstim", "MPO", "MPO")
)

test_that("facs_calc_log2fc() computes log2 fold-change of a stim condition vs unstim", {
  result <- facs_calc_log2fc(log2fc_data, ref_pop = "CD3+")

  expect_equal(nrow(result), 1L)
  expect_equal(result$restimulation, "MPO")
  expect_equal(result$population, "CD4+")

  unstim_prop <- (100 + 0.5) / 1000
  mpo_prop    <- (400 + 0.5) / 1000
  expect_equal(result$log2fc, log2(mpo_prop / unstim_prop))
})

test_that("facs_calc_log2fc() excludes ref_pop and ref_level from the output", {
  result <- facs_calc_log2fc(log2fc_data, ref_pop = "CD3+")

  expect_false("CD3+" %in% result$population)
  expect_false("unstim" %in% result$restimulation)
})

test_that("facs_calc_log2fc() pseudocount avoids -Inf when a condition has zero events", {
  data_zero <- dplyr::mutate(
    log2fc_data,
    value = dplyr::if_else(restimulation == "MPO" & population == "CD4+", 0, value)
  )

  result <- facs_calc_log2fc(data_zero, ref_pop = "CD3+")
  expect_true(is.finite(result$log2fc))
})

test_that("facs_calc_log2fc() errors when ref_pop matches more than one row per mouse_ID/tissue/restim-level combo", {
  data_dup <- dplyr::bind_rows(
    log2fc_data,
    tibble::tibble(
      file_name = "f1", population_full_path = "OtherPath/CD3+", population = "CD3+",
      metric = "count", value = 500, mouse_ID = "m1", tissue = "spleen", restimulation = "unstim"
    )
  )

  expect_error(facs_calc_log2fc(data_dup, ref_pop = "CD3+"), regexp = "m1", fixed = TRUE)
})

test_that("facs_calc_log2fc() warns and fills NA when ref_pop has no match for a combo", {
  data_missing_refpop <- dplyr::filter(log2fc_data, !(restimulation == "MPO" & population == "CD3+"))

  expect_warning(
    result <- facs_calc_log2fc(data_missing_refpop, ref_pop = "CD3+"),
    regexp = "m1",
    fixed = TRUE
  )
  expect_true(is.na(result$log2fc))
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_calc.R")'`
Expected: FAIL with `could not find function "facs_calc_log2fc"`.

- [ ] **Step 3: Implement `facs_calc_log2fc()`**

Append to `R/facs_calc.R`:

```r
#' Compute per-mouse log2 fold-change of restimulated vs unstimulated conditions
#'
#' @description
#' Normalizes each population's count to a reference population's count (an
#' event-count proportion, with a pseudocount to avoid dividing into zero),
#' then compares every non-reference \code{restim_col} level against
#' \code{ref_level} within the same \code{mouse_ID} x \code{tissue} group,
#' returning a log2 fold-change.
#'
#' @param data tibble shaped like \code{facs_read_wsp(...)$data}: must
#'   contain \code{mouse_ID}, \code{tissue}, \code{population}, \code{metric},
#'   \code{value}, and a \code{restim_col} column (e.g. joined via
#'   \code{facs_read_wsp(keywords = c("mouse_ID", "tissue", "restimulation"))}).
#' @param ref_pop character; leaf population name (matches \code{population})
#'   used as the denominator for a proportion, e.g. \code{"CD3+"}.
#' @param restim_col character; column in \code{data} holding the
#'   stimulation condition label. Default \code{"restimulation"}.
#' @param ref_level character; value in \code{restim_col} treated as the
#'   baseline; every other distinct value is compared against it. Default
#'   \code{"unstim"}.
#' @param pseudocount numeric; added to every \code{value} before computing
#'   the \code{ref_pop} proportion, avoiding \code{log2(0)}. Default
#'   \code{0.5}.
#'
#' @returns A tibble with one row per \code{mouse_ID} x \code{tissue} x
#'   \code{population} x non-reference \code{restim_col} level: columns
#'   \code{mouse_ID}, \code{tissue}, \code{population}, any passthrough
#'   columns constant within \code{mouse_ID} x \code{tissue} (e.g.
#'   \code{group}), \code{restim_col} (name reused), and \code{log2fc}.
#'   Errors if \code{restim_col} is not a column in \code{data}, or if
#'   \code{ref_pop} matches more than one row per \code{mouse_ID}/
#'   \code{tissue}/restim-level combo. Warns and fills \code{NA} if
#'   \code{ref_pop} has no match for a combo, or if \code{ref_level} is
#'   missing for a \code{mouse_ID}/\code{tissue} group.
#' @export
#'
#' @examples
#' \dontrun{
#'   dat <- facs_read_wsp(
#'     "experiment.wsp",
#'     keywords = c("mouse_ID", "tissue", "restimulation")
#'   )$data
#'   facs_calc_log2fc(dat, ref_pop = "CD3+")
#' }
facs_calc_log2fc <- function(
    data,
    ref_pop,
    restim_col = "restimulation",
    ref_level = "unstim",
    pseudocount = 0.5
) {
  if (!restim_col %in% names(data)) {
    stop(glue::glue("`restim_col` ('{restim_col}') not found in `data`."))
  }

  group_cols <- c("mouse_ID", "tissue", restim_col)

  ref_counts <- data |>
    dplyr::filter(population == ref_pop, metric == "count") |>
    dplyr::select(dplyr::all_of(group_cols), ref_count = value)

  dup_combos <- ref_counts |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    dplyr::filter(n > 1L)
  if (nrow(dup_combos) > 0L) {
    dup_desc <- purrr::pmap_chr(dup_combos[group_cols], function(...) paste(..., sep = "/"))
    stop(glue::glue(
      "ref_pop '{ref_pop}' matches more than one population for mouse_ID/tissue/{restim_col} combo(s): ",
      "{paste(dup_desc, collapse = ', ')}. Leaf population name is ambiguous."
    ))
  }

  proportions <- data |>
    dplyr::filter(metric == "count", population != ref_pop) |>
    dplyr::left_join(ref_counts, by = group_cols) |>
    dplyr::mutate(proportion = (value + pseudocount) / ref_count)

  missing_ref_pop <- proportions |>
    dplyr::filter(is.na(ref_count)) |>
    dplyr::distinct(dplyr::across(dplyr::all_of(group_cols)))
  if (nrow(missing_ref_pop) > 0L) {
    missing_desc <- purrr::pmap_chr(missing_ref_pop[group_cols], function(...) paste(..., sep = "/"))
    warning(glue::glue(
      "ref_pop '{ref_pop}' not found for mouse_ID/tissue/{restim_col} combo(s): ",
      "{paste(missing_desc, collapse = ', ')}. Result filled with NA."
    ))
  }

  id_candidates <- setdiff(
    names(data),
    c("mouse_ID", "tissue", restim_col, "population", "population_full_path", "metric", "value", "file_name")
  )
  is_constant_ <- function(col) {
    n_distinct_per_group <- proportions |>
      dplyr::group_by(mouse_ID, tissue) |>
      dplyr::summarise(n_distinct = dplyr::n_distinct(.data[[col]]), .groups = "drop") |>
      dplyr::pull(n_distinct)
    all(n_distinct_per_group == 1L)
  }
  passthrough_cols <- id_candidates[purrr::map_lgl(id_candidates, is_constant_)]
  id_cols <- c("mouse_ID", "tissue", "population", passthrough_cols)

  wide <- proportions |>
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(id_cols),
      names_from = dplyr::all_of(restim_col),
      values_from = proportion
    )

  if (!ref_level %in% names(wide)) {
    wide[[ref_level]] <- NA_real_
  }

  missing_ref_level <- dplyr::filter(wide, is.na(.data[[ref_level]]))
  if (nrow(missing_ref_level) > 0L) {
    missing_desc <- purrr::pmap_chr(
      missing_ref_level[c("mouse_ID", "tissue")],
      function(...) paste(..., sep = "/")
    )
    warning(glue::glue(
      "ref_level '{ref_level}' missing for mouse_ID/tissue combo(s): ",
      "{paste(unique(missing_desc), collapse = ', ')}. Result filled with NA."
    ))
  }

  stim_cols <- setdiff(names(wide), c(id_cols, ref_level))

  wide |>
    dplyr::mutate(dplyr::across(
      dplyr::all_of(stim_cols),
      function(x) log2(x / .data[[ref_level]])
    )) |>
    dplyr::select(!dplyr::all_of(ref_level)) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(stim_cols),
      names_to = restim_col,
      values_to = "log2fc"
    )
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_calc.R")'`
Expected: PASS, 20 tests total, 0 failures.

- [ ] **Step 5: Document and commit**

```bash
Rscript -e 'devtools::document()'
git add R/facs_calc.R tests/testthat/test-facs_calc.R man/facs_calc_log2fc.Rd NAMESPACE
git commit -m "feat: add facs_calc_log2fc() for T-cell restimulation normalization"
```

---

## Task 2: Passthrough columns, multi-condition, missing ref_level, custom column names

**Files:**
- Modify: none (implementation already complete from Task 1 — this task is test-only, verifying branches not yet exercised: multiple simultaneous non-reference conditions, passthrough-column auto-detection, missing `ref_level`, `restim_col` validation, and non-default argument names)
- Test: `tests/testthat/test-facs_calc.R`

**Interfaces:**
- Consumes: `facs_calc_log2fc()` from Task 1 (already implements the full logic — this task's job is proving every remaining branch behaves as specified).

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-facs_calc.R`:

```r
test_that("facs_calc_log2fc() warns and fills NA when ref_level is missing for a mouse_ID/tissue group", {
  data_no_unstim <- dplyr::filter(log2fc_data, restimulation != "unstim")

  expect_warning(
    result <- facs_calc_log2fc(data_no_unstim, ref_pop = "CD3+"),
    regexp = "m1",
    fixed = TRUE
  )
  expect_true(is.na(result$log2fc))
})

test_that("facs_calc_log2fc() errors when restim_col is not a column in data", {
  data_no_restim <- dplyr::select(log2fc_data, !restimulation)

  expect_error(facs_calc_log2fc(data_no_restim, ref_pop = "CD3+"), regexp = "restim_col", fixed = TRUE)
})

test_that("facs_calc_log2fc() handles more than one non-reference restimulation level", {
  data_multi_stim <- dplyr::bind_rows(
    log2fc_data,
    tibble::tibble(
      file_name = "f3", population_full_path = c("CD3+", "CD3+/CD4+"), population = c("CD3+", "CD4+"),
      metric = "count", value = c(1000, 250), mouse_ID = "m1", tissue = "spleen", restimulation = "PMA-Iono"
    )
  )

  result <- facs_calc_log2fc(data_multi_stim, ref_pop = "CD3+")

  expect_setequal(result$restimulation, c("MPO", "PMA-Iono"))
  expect_equal(nrow(result), 2L)
})

test_that("facs_calc_log2fc() carries through a passthrough column constant within mouse_ID/tissue", {
  data_group <- dplyr::mutate(log2fc_data, group = "WT")

  result <- facs_calc_log2fc(data_group, ref_pop = "CD3+")
  expect_true("group" %in% names(result))
  expect_equal(unique(result$group), "WT")
})

test_that("facs_calc_log2fc() drops a passthrough-candidate column that varies within mouse_ID/tissue", {
  data_varying <- dplyr::mutate(
    log2fc_data,
    batch = dplyr::if_else(restimulation == "unstim", "batch1", "batch2")
  )

  result <- facs_calc_log2fc(data_varying, ref_pop = "CD3+")
  expect_false("batch" %in% names(result))
})

test_that("facs_calc_log2fc() works with non-default restim_col/ref_level names", {
  data_custom <- log2fc_data |>
    dplyr::rename(condition = restimulation) |>
    dplyr::mutate(condition = dplyr::if_else(condition == "unstim", "baseline", "MPO"))

  result <- facs_calc_log2fc(data_custom, ref_pop = "CD3+", restim_col = "condition", ref_level = "baseline")
  expect_equal(result$condition, "MPO")
})
```

- [ ] **Step 2: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-facs_calc.R")'`
Expected: PASS, 26 tests total, 0 failures — Task 1's implementation already covers these branches. If any of these 6 new tests fail, fix `facs_calc_log2fc()` in `R/facs_calc.R` (do not weaken the test) before proceeding.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-facs_calc.R
git commit -m "test: cover multi-condition, passthrough columns, and custom column names in facs_calc_log2fc()"
```

---

## Task 3: Final check, CLAUDE.md housekeeping

**Files:**
- Modify: `CLAUDE.md`

**Interfaces:**
- Consumes: `facs_calc_log2fc()` (Tasks 1–2).

- [ ] **Step 1: Run full test suite**

Run: `Rscript -e 'devtools::test()'`
Expected: All tests pass, including the 26 in `test-facs_calc.R` and all pre-existing tests in the other `test-*.R` files.

- [ ] **Step 2: Run `devtools::check()` and record the actual NOTEs**

Run: `Rscript -e 'devtools::check()'`
Expected: 0 errors, 0 warnings. Note any new "no visible binding for global variable" NOTEs for `R/facs_calc.R` attributable to `facs_calc_log2fc()` — these are tidy-eval false positives (bare column names inside `dplyr` verbs, e.g. `population`, `metric`, `value`, `mouse_ID`, `tissue`, `ref_count`, `n`, plus the `.data` pronoun), matching the existing pattern already documented for `facs_calc_pct_of`/`facs_calc_count_per_g` in CLAUDE.md. Write down the exact variable names reported — Step 3 needs them.

- [ ] **Step 3: Update CLAUDE.md's "Known check output" section**

In `CLAUDE.md`, find this line (near the end of the "Known check output" section):

```
  - `pivot_organ_weights_long_` (the internal `meta_clean()` helper that pivots `organ_weights` to long format) variable-binding note (`mouse_ID`) — bare column name inside a `tidyr::pivot_longer()` deselection (`cols = !mouse_ID`); valid tidy eval, flagged as a static-analysis false positive, same pattern as the existing `meta_read`/`meta_annotate` note above. (`meta_clean_sheet_`, the internal `meta_read()` per-sheet helper, also produces `mouse_id`/`group` notes -- already covered by the `meta_read`/`meta_annotate` bullet above.)
```

Add immediately after it, using the exact variable names recorded in Step 2 (adjust the list below if `devtools::check()` reported different names):

```
  - All `facs_calc_log2fc` variable-binding notes (`population`, `metric`, `value`, `mouse_ID`, `tissue`, `ref_count`, `n`) — bare column names inside dplyr verbs, same false-positive pattern as `facs_calc_pct_of`/`facs_calc_count_per_g` above. The internal `is_constant_` closure's `mouse_ID`/`tissue` references (bare inside `dplyr::group_by()`) are also flagged for the same reason.
```

- [ ] **Step 4: Final full test suite + check run**

Run: `Rscript -e 'devtools::test()'`
Expected: All tests pass.

Run: `Rscript -e 'devtools::check()'`
Expected: 0 errors, 0 warnings (NOTEs limited to the pre-existing ones already listed in CLAUDE.md, plus the new `facs_calc_log2fc` variable-binding NOTE now documented in Step 3).

- [ ] **Step 5: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: note facs_calc_log2fc variable-binding NOTEs in CLAUDE.md"
```

---

## Self-Review Notes

- **Spec coverage:** Core log2FC computation, pseudocount, ref_pop/ref_level exclusion, ambiguity error, missing-ref_pop warning — Task 1. Missing-ref_level warning, restim_col-not-found error, multi-condition handling, passthrough-column auto-detection (both keep and drop cases), non-default `restim_col`/`ref_level` names — Task 2. `devtools::check()` verification and CLAUDE.md housekeeping — Task 3. All design-doc sections (Public API, Error handling, Testing) covered. Per the design doc, no end-to-end fixture test is included — `tests/fixtures/minimal.wsp` doesn't carry `mouse_ID`/`tissue`/`restimulation` keywords, consistent with `facs_calc_pct_of()` also being unit-tested only for that reason (`facs_calc_count_per_g()` is the one exception with fixture keywords already present).
- **Type consistency:** `facs_calc_log2fc()`'s signature (`data, ref_pop, restim_col, ref_level, pseudocount`) is identical across Task 1's implementation and every call site in Tasks 1–2's tests.
- **No placeholders:** every step has complete, runnable code or an exact command with expected output.
