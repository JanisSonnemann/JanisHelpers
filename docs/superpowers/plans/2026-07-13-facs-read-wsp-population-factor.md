# facs_read_wsp() population Gating-Order Factor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `facs_read_wsp()`'s `data$population` column a plain `factor` (not `ordered`), with level order matching the depth-first traversal order of `population_full_path` (i.e. the FlowJo gating hierarchy), instead of a plain character column that sorts alphabetically downstream.

**Architecture:** `walk_pops_()` already emits `pops` tibble rows in depth-first pre-order relative to `population_full_path` (a population's own `count`/`fraction_of_parent` rows precede its descendants' rows). Add one small internal helper, `population_gating_order_()`, that derives the leaf-name level order from that existing row order, and apply it as a `factor()` conversion inside `facs_read_wsp()` right after `data <- pops`. No changes needed to `walk_pops_()` or `parse_populations_()`.

**Tech Stack:** R, `xml2`, `dplyr`, `tibble`, `testthat`.

## Global Constraints

- Namespace every external call as `pkg::fn()` inside exported function bodies (no bare `dplyr`/`purrr`/etc. calls) — per `CLAUDE.md` code style.
- Non-ASCII characters in R source must be escaped as `\uXXXX`.
- `devtools::document()` must be run after any roxygen `@param`/`@returns` change.
- `devtools::check()` target: 0 errors, 0 warnings (pre-existing NOTEs listed in `CLAUDE.md` are expected and not regressions).
- `population` must be a plain `factor` (level order fixed), **not** an `ordered` factor — avoids changing `rstatix`/`gtsummary` statistical treatment in `analysis_stats.R`.
- Design spec: `docs/superpowers/specs/2026-07-13-facs-read-wsp-population-factor-design.md`.

---

### Task 1: `population_gating_order_()` helper + factor conversion in `facs_read_wsp()`

**Files:**
- Modify: `R/facs_read.R` (insert helper immediately after `parse_populations_()` ends, before the `# Exported function` section marker — around line 190)
- Modify: `R/facs_read.R:295` (apply factor conversion right after `data <- pops`)
- Test: `tests/testthat/test-facs_read.R`

**Interfaces:**
- Consumes: `pops` tibble produced by `parse_populations_(doc, sample_ids)` (existing, unchanged) — columns `file_name` (chr), `population_full_path` (chr), `population` (chr), `metric` (chr), `value` (dbl).
- Produces: `population_gating_order_(pops)` — takes the `pops` tibble, returns a character vector of unique leaf population names in first-encountered depth-first traversal order. Consumed only within `facs_read_wsp()` in this same task.
- After this task, `facs_read_wsp(...)$data$population` is `factor` (not `ordered`), with `levels()` equal to `population_gating_order_(pops)`.

- [ ] **Step 1: Write the failing tests**

Add to `tests/testthat/test-facs_read.R`, after the existing test block ending at line 148 (`"populations gated downstream of a boolean gate are exported"`):

```r
test_that("population is a factor ordered by depth-first gating hierarchy traversal", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))

  expect_true(is.factor(result$data$population))
  expect_false(is.ordered(result$data$population))

  expected_levels <- unique(basename(unique(result$data$population_full_path)))
  expect_equal(levels(result$data$population), expected_levels)
})

test_that("ancestor populations sort before their descendants in the population factor", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))

  paths <- unique(result$data$population_full_path)
  level_index <- function(full_path) {
    match(basename(full_path), levels(result$data$population))
  }

  for (child_path in paths) {
    parent_path <- dirname(child_path)
    if (parent_path %in% paths && parent_path != child_path) {
      expect_lt(level_index(parent_path), level_index(child_path))
    }
  }
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-facs_read.R", reporter = "summary")'`

Expected: FAIL on both new tests — `is.factor(result$data$population)` is `FALSE` (currently character), since `population` is not yet converted to a factor.

- [ ] **Step 3: Add the `population_gating_order_()` helper**

In `R/facs_read.R`, insert immediately before the `# Exported function` section marker (currently lines 191-193, right after `parse_populations_()` ends):

```r
population_gating_order_ <- function(pops) {
  full_path_order <- unique(pops$population_full_path)
  unique(basename(full_path_order))
}

```

- [ ] **Step 4: Apply the factor conversion in `facs_read_wsp()`**

In `R/facs_read.R`, at line 295 (`data <- pops`), change:

```r
  # Build data: population rows + optional keyword join
  data <- pops
```

to:

```r
  # Build data: population rows + optional keyword join
  data <- pops
  data$population <- factor(data$population, levels = population_gating_order_(pops))
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-facs_read.R", reporter = "summary")'`

Expected: PASS, all tests in the file green (22 existing + 2 new = 24 tests, 0 failures).

- [ ] **Step 6: Run the full test suite to check for regressions in downstream consumers**

Run: `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_dir("tests/testthat", reporter = "summary")'`

Expected: no new failures in `test-facs_calc.R` or `test-analysis_stats.R` (both consume `population` via `==`/`!=`/`paste()`/`group_by()`, all factor-compatible — see design spec's "Downstream impact" section). Any fixture-skip messages (e.g. `Treg fixture not available`) are expected and unrelated.

- [ ] **Step 7: Commit**

```bash
git add R/facs_read.R tests/testthat/test-facs_read.R
git commit -m "feat: order facs_read_wsp() population as a gating-hierarchy factor"
```

---

### Task 2: Documentation updates

**Files:**
- Modify: `R/facs_read.R:207-221` (roxygen `@returns` for `facs_read_wsp()`)
- Modify: `man/facs_read_wsp.Rd` (regenerated, not hand-edited)
- Modify: `CLAUDE.md:73-83` (Data structures table + bullets)

**Interfaces:**
- Consumes: nothing new — this task only updates documentation to describe Task 1's already-implemented behavior.
- Produces: nothing consumed by later tasks (this is the final task in the plan).

- [ ] **Step 1: Update the `facs_read_wsp()` roxygen `@returns`**

In `R/facs_read.R`, the `data` slot's `@returns` description (currently lines 209-213):

```r
#'     \item{\code{data}}{Long-format tibble, one row per
#'       \code{file_name x population_full_path x metric}.
#'       Columns: \code{file_name}, \code{population_full_path},
#'       \code{population}, \code{metric}, \code{value}, plus any
#'       requested \code{keywords}.}
```

replace with:

```r
#'     \item{\code{data}}{Long-format tibble, one row per
#'       \code{file_name x population_full_path x metric}.
#'       Columns: \code{file_name}, \code{population_full_path} (chr),
#'       \code{population} (factor; leaf gate name, levels ordered by
#'       first-encountered depth-first traversal of the gating hierarchy,
#'       matching \code{population_full_path} order -- if the same leaf
#'       name occurs under more than one parent it collapses to a single
#'       level positioned at its first occurrence; use
#'       \code{population_full_path} to disambiguate), \code{metric},
#'       \code{value}, plus any requested \code{keywords}.}
```

- [ ] **Step 2: Regenerate roxygen docs**

Run: `Rscript -e 'devtools::document()'`

Expected: `man/facs_read_wsp.Rd` updated to reflect the new `@returns` text; no roxygen warnings/errors.

- [ ] **Step 3: Update `CLAUDE.md`'s Data structures table**

In `CLAUDE.md`, in the `facs_import_wsp()` table (currently lines 73-80), change:

```
| `Population` | chr | Leaf gate name (`basename(PopulationFullPath)`) |
```

to:

```
| `Population` | factor | Leaf gate name (`basename(PopulationFullPath)`); levels ordered by first-encountered depth-first gating-hierarchy traversal (matching `PopulationFullPath` order), not alphabetically |
```

- [ ] **Step 4: Run the full test suite and `devtools::check()`**

Run: `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_dir("tests/testthat", reporter = "summary")'`

Expected: all tests pass (same as Task 1 Step 6).

Run: `Rscript -e 'devtools::check()'`

Expected: 0 errors, 0 warnings; only the pre-existing NOTEs already documented in `CLAUDE.md`'s "Known check output" section.

- [ ] **Step 5: Commit**

```bash
git add R/facs_read.R man/facs_read_wsp.Rd CLAUDE.md
git commit -m "docs: document facs_read_wsp() population gating-order factor"
```
