# plot_bubble_fc() Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `plot_bubble_fc()`, a new exported function producing a bubble plot that compares one or more groups against a named control group across FACS populations, colored/sized by log2 fold-change and annotated with significance stars.

**Architecture:** One exported function `plot_bubble_fc()` in a new `R/plot_bubble.R` file, backed by an unexported stats helper `calc_bubble_fc_()` (plus its own helpers `calc_bubble_p_values_()` and `adjust_p_()`). `calc_bubble_fc_()` returns a tidy tibble (population, comparison, log2fc, p_value, stars); `plot_bubble_fc()` turns that into a `ggplot2` bubble plot. This introduces a new `plot_` domain.

**Tech Stack:** R, `ggplot2` (new dependency), `stats` (new dependency, namespaced), `rstatix::dunn_test()` (already a dependency), `dplyr`/`tidyr`/`purrr`/`tibble`/`glue` (already dependencies), `testthat` 3rd edition.

## Global Constraints

- Pipe: use `|>` (base pipe) everywhere; no `%>%`.
- Namespacing: every external call inside an exported/internal function body is `pkg::fn()` — no bare `map()`, `filter()`, etc.
- Bare column names inside `dplyr`/`tidyr` verbs are fine (tidy eval), not a namespacing violation.
- Non-ASCII characters in R source must be escaped as `\uXXXX`.
- Every exported function needs `@param`, `@returns`, `@export`; `@examples` uses `\dontrun{}` when it needs a realistic multi-column input not easily constructed inline.
- `plot_bubble_fc()` returns a `ggplot` object **visibly** (no `invisible()` — that convention is for data-import side-effect functions, not primary-purpose plot output).
- Target after all tasks: `devtools::check()` → 0 errors, 0 warnings.
- Run `devtools::load_all()` before interactively testing; run the specific `testthat::test_file()` after each implementation step.

---

### Task 1: Package scaffolding — DESCRIPTION + CLAUDE.md domain

**Files:**
- Modify: `DESCRIPTION`
- Modify: `CLAUDE.md`

**Interfaces:**
- Produces: `ggplot2` and `stats` available as namespaced dependencies for all later tasks.

- [ ] **Step 1: Add `ggplot2` and `stats` to `DESCRIPTION`'s `Imports`**

Edit `DESCRIPTION`, replacing the `Imports:` block:

```
Imports:
    CytoML,
    dplyr,
    fcexpr,
    flowCore,
    flowWorkspace,
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
    tools,
    xml2,
    xfun
```

with:

```
Imports:
    CytoML,
    dplyr,
    fcexpr,
    flowCore,
    flowWorkspace,
    ggplot2,
    glue,
    gt,
    gtsummary,
    janitor,
    purrr,
    readxl,
    rmarkdown,
    rstatix,
    stats,
    stringr,
    tibble,
    tidyr,
    tools,
    xml2,
    xfun
```

- [ ] **Step 2: Add the `plot_` domain to `CLAUDE.md`'s domain list**

Edit `CLAUDE.md`, in the "File and function naming" section, replacing:

```
- `facs_` — FlowJo / flow cytometry
- `analysis_` — statistical summaries and tests
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
- `wrangle_` — general data wrangling [stub — no functions yet]
- `db_` — database access [stub — no functions yet]
```

- [ ] **Step 3: Update the "Current Imports" line in CLAUDE.md's Dependency philosophy section**

Replace:

```
Current `Imports`: `fcexpr`, `dplyr`, `tidyr`, `stringr`, `tibble`, `purrr`, `glue`, `rmarkdown`, `xfun`, `gt`, `gtsummary`, `rstatix`.
```

with:

```
Current `Imports`: `fcexpr`, `dplyr`, `tidyr`, `stringr`, `tibble`, `purrr`, `glue`, `rmarkdown`, `xfun`, `gt`, `gtsummary`, `rstatix`, `ggplot2`, `stats`.
```

- [ ] **Step 4: Verify the package still loads**

Run: `Rscript -e "devtools::load_all()"`
Expected: loads with no errors (no code changes yet, this just confirms `DESCRIPTION` is still well-formed).

- [ ] **Step 5: Commit**

```bash
git add DESCRIPTION CLAUDE.md
git commit -m "chore: add plot_ domain and ggplot2/stats dependencies"
```

---

### Task 2: `calc_bubble_fc_()` — input validation

**Files:**
- Create: `R/plot_bubble.R`
- Create: `tests/testthat/test-plot_bubble.R`

**Interfaces:**
- Produces: `calc_bubble_fc_(data, control, group_col, population_col, value_col, test, p_adjust_method, summary_fun, pseudocount)` — unexported, currently returns `NULL` after validation passes. Later tasks fill in the real return value.

- [ ] **Step 1: Write the failing validation tests**

Create `tests/testthat/test-plot_bubble.R`:

```r
test_that("calc_bubble_fc_() errors when a required column is missing", {
  data <- tibble::tibble(group = c("control", "treated"), value = c(1, 2))
  expect_error(
    calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value",
      test = "auto", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5
    ),
    "population"
  )
})

test_that("calc_bubble_fc_() errors when control is not present in group_col", {
  data <- tibble::tibble(
    group = c("a", "b"), population = c("CD4", "CD4"), value = c(1, 2)
  )
  expect_error(
    calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value",
      test = "auto", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5
    ),
    "control"
  )
})

test_that("calc_bubble_fc_() errors when metric has more than one distinct value", {
  data <- tibble::tibble(
    group = c("control", "treated"), population = c("CD4", "CD4"),
    value = c(1, 2), metric = c("Count", "FractionOfParent")
  )
  expect_error(
    calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value",
      test = "auto", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5
    ),
    "metric"
  )
})

test_that("calc_bubble_fc_() errors when tissue has more than one distinct value", {
  data <- tibble::tibble(
    group = c("control", "treated"), population = c("CD4", "CD4"),
    value = c(1, 2), tissue = c("Spleen", "Blood")
  )
  expect_error(
    calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value",
      test = "auto", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5
    ),
    "tissue"
  )
})

test_that("calc_bubble_fc_() errors on an invalid test argument", {
  data <- tibble::tibble(
    group = c("control", "treated"), population = c("CD4", "CD4"), value = c(1, 2)
  )
  expect_error(
    calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value",
      test = "not-a-test", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5
    ),
    "test"
  )
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-plot_bubble.R')"`
Expected: FAIL with `could not find function "calc_bubble_fc_"`.

- [ ] **Step 3: Implement validation in `R/plot_bubble.R`**

Create `R/plot_bubble.R`:

```r
# Stats engine behind plot_bubble_fc(): computes per-population,
# per-comparison log2 fold-change and test statistics.
calc_bubble_fc_ <- function(data, control, group_col, population_col,
                             value_col, test, p_adjust_method,
                             summary_fun, pseudocount) {

  missing_cols <- setdiff(c(group_col, population_col, value_col), names(data))
  if (length(missing_cols) > 0L) {
    stop(glue::glue("`data` is missing column(s): {paste(missing_cols, collapse = ', ')}."))
  }
  if (!control %in% data[[group_col]]) {
    stop(glue::glue("`control` ('{control}') not found in `data${group_col}`."))
  }
  if ("metric" %in% names(data) && dplyr::n_distinct(data[["metric"]]) > 1L) {
    stop("`data` contains more than one distinct `metric`; pre-filter to a single metric.")
  }
  if ("tissue" %in% names(data) && dplyr::n_distinct(data[["tissue"]]) > 1L) {
    stop("`data` contains more than one distinct `tissue`; pre-filter to a single tissue.")
  }
  if (!test %in% c("auto", "wilcox", "t.test", "kruskal")) {
    stop('`test` must be one of "auto", "wilcox", "t.test", "kruskal".')
  }

  NULL
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-plot_bubble.R')"`
Expected: PASS, 5 tests.

- [ ] **Step 5: Commit**

```bash
git add R/plot_bubble.R tests/testthat/test-plot_bubble.R
git commit -m "feat: validate inputs in calc_bubble_fc_()"
```

---

### Task 3: `calc_bubble_fc_()` — log2fc, ordering, untestable-cell warning

**Files:**
- Modify: `R/plot_bubble.R`
- Modify: `tests/testthat/test-plot_bubble.R`

**Interfaces:**
- Consumes: the validation block from Task 2 (unchanged).
- Produces: `calc_bubble_fc_()` now returns a tibble with columns `population` (factor, first-appearance order), `comparison` (factor, first-appearance order, excludes `control`), `log2fc` (dbl), `p_value` (dbl, currently always `NA_real_` — filled for real in Task 4), `stars` (chr, currently always `""` — filled for real in Task 4). Emits a `warning()` per population/comparison combo with <2 observations in either arm.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-plot_bubble.R`:

```r
make_data_ <- function(...) {
  groups <- list(...)
  purrr::imap_dfr(groups, function(values, group_pop) {
    parts <- strsplit(group_pop, "\\.")[[1]]
    tibble::tibble(group = parts[1], population = parts[2], value = values)
  })
}

test_that("calc_bubble_fc_() computes log2 fold-change per population", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    treated.CD4 = c(19.6, 19.9, 20.1, 20.4),
    control.CD8 = c(4.8, 4.9, 5.1, 5.2),
    treated.CD8 = c(4.85, 4.95, 5.05, 5.15)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "wilcox", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0
  )
  expect_equal(result$log2fc[result$population == "CD4"], 1, tolerance = 1e-8)
  expect_equal(result$log2fc[result$population == "CD8"], 0, tolerance = 1e-8)
})

test_that("calc_bubble_fc_() orders population/comparison by first appearance and excludes control", {
  data <- make_data_(
    control.CD8 = c(4.8, 4.9, 5.1, 5.2),
    treated.CD8 = c(4.85, 4.95, 5.05, 5.15),
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    treated.CD4 = c(19.6, 19.9, 20.1, 20.4)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "wilcox", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0
  )
  expect_equal(levels(result$population), c("CD8", "CD4"))
  expect_true(is.factor(result$comparison))
  expect_false("control" %in% as.character(result$comparison))
})

test_that("calc_bubble_fc_() warns and computes a finite log2fc for <2 observations in one arm", {
  data <- tibble::tibble(
    group = c("control", "treated", "treated", "treated"),
    population = c("CD3", "CD3", "CD3", "CD3"),
    value = c(10, 19.6, 19.9, 20.1)
  )
  expect_warning(
    result <- calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value",
      test = "wilcox", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5
    ),
    "Fewer than 2 observations"
  )
  expect_true(is.finite(result$log2fc))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-plot_bubble.R')"`
Expected: FAIL — `calc_bubble_fc_()` currently returns `NULL`, so `result$log2fc` etc. error or return `NULL`.

- [ ] **Step 3: Implement the fc_table shape**

In `R/plot_bubble.R`, replace the final `NULL` in `calc_bubble_fc_()` with:

```r
  population_levels <- unique(data[[population_col]])
  group_levels <- unique(data[[group_col]])
  comparison_levels <- setdiff(group_levels, control)

  fc_table <- tidyr::expand_grid(
    population = population_levels,
    comparison = comparison_levels
  ) |>
    purrr::pmap_dfr(function(population, comparison) {
      control_values <- data[[value_col]][
        data[[population_col]] == population & data[[group_col]] == control
      ]
      group_values <- data[[value_col]][
        data[[population_col]] == population & data[[group_col]] == comparison
      ]
      log2fc <- log2(
        (summary_fun(group_values) + pseudocount) / (summary_fun(control_values) + pseudocount)
      )
      testable <- length(control_values) >= 2L && length(group_values) >= 2L
      if (!testable) {
        warning(glue::glue(
          "Fewer than 2 observations for population '{population}', ",
          "comparison '{comparison}' vs control '{control}'; p-value set to NA."
        ))
      }
      tibble::tibble(
        population = population, comparison = comparison,
        log2fc = log2fc, testable = testable
      )
    })

  fc_table |>
    dplyr::mutate(p_value = NA_real_, stars = "") |>
    dplyr::select(!testable) |>
    dplyr::mutate(
      population = factor(population, levels = population_levels),
      comparison = factor(comparison, levels = comparison_levels)
    )
```

(The complete `calc_bubble_fc_()` function now ends with this block instead of `NULL`.)

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-plot_bubble.R')"`
Expected: PASS, 8 tests.

- [ ] **Step 5: Commit**

```bash
git add R/plot_bubble.R tests/testthat/test-plot_bubble.R
git commit -m "feat: compute log2fc and ordering in calc_bubble_fc_()"
```

---

### Task 4: `calc_bubble_fc_()` — wilcox/t.test dispatch, p-adjustment, stars

**Files:**
- Modify: `R/plot_bubble.R`
- Modify: `tests/testthat/test-plot_bubble.R`

**Interfaces:**
- Consumes: `fc_table` (with `testable` column) from Task 3's construction step, before it's dropped.
- Produces: two new unexported helpers — `adjust_p_(p, method)` (NA-aware wrapper around `stats::p.adjust()`) and `calc_bubble_p_values_(data, control, group_col, population_col, value_col, resolved_test, p_adjust_method, fc_table)` (currently handles `resolved_test %in% c("wilcox", "t.test")`; `"kruskal"` added in Task 5). `calc_bubble_fc_()` now computes real `p_value`/`stars`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-plot_bubble.R`:

```r
test_that("calc_bubble_fc_() computes a wilcox p-value for a 2-group comparison", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    treated.CD4 = c(19.6, 19.9, 20.1, 20.4)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "wilcox", p_adjust_method = "none", summary_fun = mean, pseudocount = 0
  )
  expected_p <- stats::wilcox.test(
    data$value[data$group == "control"], data$value[data$group == "treated"]
  )$p.value
  expect_equal(result$p_value, expected_p, tolerance = 1e-8)
})

test_that("calc_bubble_fc_() uses stats::t.test when test = \"t.test\"", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    treated.CD4 = c(19.6, 19.9, 20.1, 20.4)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "t.test", p_adjust_method = "none", summary_fun = mean, pseudocount = 0
  )
  expected_p <- stats::t.test(
    data$value[data$group == "control"], data$value[data$group == "treated"]
  )$p.value
  expect_equal(result$p_value, expected_p, tolerance = 1e-8)
})

test_that("calc_bubble_fc_() adjusts p-values per comparison column", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    treated.CD4 = c(19.6, 19.9, 20.1, 20.4),
    control.CD8 = c(4.8, 4.9, 5.1, 5.2),
    treated.CD8 = c(4.85, 4.95, 5.05, 5.15)
  )
  raw <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "wilcox", p_adjust_method = "none", summary_fun = mean, pseudocount = 0
  )
  adjusted <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "wilcox", p_adjust_method = "bonferroni", summary_fun = mean, pseudocount = 0
  )
  expect_equal(
    sort(adjusted$p_value),
    sort(stats::p.adjust(raw$p_value, method = "bonferroni")),
    tolerance = 1e-8
  )
})

test_that("calc_bubble_fc_() assigns stars based on p-value thresholds", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    treated.CD4 = c(119.6, 119.9, 120.1, 120.4),
    control.CD8 = c(10, 10, 10, 10.001),
    treated.CD8 = c(10, 10.001, 10, 10)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "t.test", p_adjust_method = "none", summary_fun = mean, pseudocount = 0
  )
  expect_equal(result$stars[result$population == "CD4"], "***")
  expect_equal(result$stars[result$population == "CD8"], "")
})

test_that("calc_bubble_fc_() leaves p_value/stars NA/blank for <2-observation cells", {
  data <- tibble::tibble(
    group = c("control", "treated", "treated", "treated"),
    population = c("CD3", "CD3", "CD3", "CD3"),
    value = c(10, 19.6, 19.9, 20.1)
  )
  result <- suppressWarnings(calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "wilcox", p_adjust_method = "none", summary_fun = mean, pseudocount = 0.5
  ))
  expect_true(is.na(result$p_value))
  expect_equal(result$stars, "")
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-plot_bubble.R')"`
Expected: FAIL — `result$p_value` is always `NA_real_` and `result$stars` always `""` from Task 3's placeholder.

- [ ] **Step 3: Implement dispatch, p-adjustment, and stars**

In `R/plot_bubble.R`, replace the final return block of `calc_bubble_fc_()` (the `fc_table |> dplyr::mutate(p_value = NA_real_, stars = "") |> ...` block from Task 3) with:

```r
  resolved_test <- test
  if (test == "auto") {
    resolved_test <- if (length(group_levels) > 2L) "kruskal" else "wilcox"
  }

  p_values <- calc_bubble_p_values_(
    data, control, group_col, population_col, value_col,
    resolved_test, p_adjust_method, fc_table
  )

  fc_table |>
    dplyr::select(!testable) |>
    dplyr::left_join(p_values, by = c("population", "comparison")) |>
    dplyr::mutate(
      stars = dplyr::case_when(
        is.na(p_value) ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE ~ ""
      ),
      population = factor(population, levels = population_levels),
      comparison = factor(comparison, levels = comparison_levels)
    )
```

Note `resolved_test` must be computed before `fc_table` is built (move the `resolved_test <- test; if (test == "auto") ...` block up, directly after the `comparison_levels <- setdiff(group_levels, control)` line from Task 3, since it doesn't depend on `fc_table`).

Then add the two new helpers at the bottom of `R/plot_bubble.R`:

```r
adjust_p_ <- function(p, method) {
  out <- rep(NA_real_, length(p))
  non_na <- !is.na(p)
  out[non_na] <- stats::p.adjust(p[non_na], method = method)
  out
}

calc_bubble_p_values_ <- function(data, control, group_col, population_col,
                                   value_col, resolved_test, p_adjust_method,
                                   fc_table) {
  if (resolved_test %in% c("wilcox", "t.test")) {
    raw <- fc_table |>
      purrr::pmap_dfr(function(population, comparison, log2fc, testable) {
        p_value <- NA_real_
        if (testable) {
          control_values <- data[[value_col]][
            data[[population_col]] == population & data[[group_col]] == control
          ]
          group_values <- data[[value_col]][
            data[[population_col]] == population & data[[group_col]] == comparison
          ]
          p_value <- if (resolved_test == "wilcox") {
            stats::wilcox.test(control_values, group_values)$p.value
          } else {
            stats::t.test(control_values, group_values)$p.value
          }
        }
        tibble::tibble(population = population, comparison = comparison, p_value = p_value)
      })

    raw |>
      dplyr::group_by(comparison) |>
      dplyr::mutate(p_value = adjust_p_(p_value, p_adjust_method)) |>
      dplyr::ungroup()
  }
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-plot_bubble.R')"`
Expected: PASS, 13 tests.

- [ ] **Step 5: Commit**

```bash
git add R/plot_bubble.R tests/testthat/test-plot_bubble.R
git commit -m "feat: wilcox/t.test dispatch, p-adjustment, and stars in calc_bubble_fc_()"
```

---

### Task 5: `calc_bubble_fc_()` — kruskal + Dunn post-hoc dispatch

**Files:**
- Modify: `R/plot_bubble.R`
- Modify: `tests/testthat/test-plot_bubble.R`

**Interfaces:**
- Consumes: `calc_bubble_p_values_()` from Task 4 (extended with a new branch).
- Produces: `resolved_test == "kruskal"` now returns real Dunn-test-derived p-values (control-vs-group rows only), whether reached via `test = "auto"` with >2 groups or `test = "kruskal"` explicitly.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-plot_bubble.R`:

```r
test_that("calc_bubble_fc_() dispatches to kruskal+dunn when test = \"auto\" and there are >2 groups", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    low.CD4 = c(14.6, 14.9, 15.1, 15.4),
    high.CD4 = c(19.6, 19.9, 20.1, 20.4)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "auto", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0
  )
  expect_equal(sort(as.character(result$comparison)), c("high", "low"))
  expect_true(all(!is.na(result$p_value)))
  expect_equal(result$log2fc[result$comparison == "low"], log2(1.5), tolerance = 1e-8)
  expect_equal(result$log2fc[result$comparison == "high"], log2(2), tolerance = 1e-8)
})

test_that("calc_bubble_fc_() forces kruskal+dunn even with only 2 groups when test = \"kruskal\"", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    treated.CD4 = c(19.6, 19.9, 20.1, 20.4)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "kruskal", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0
  )
  expect_false(is.na(result$p_value))
})

test_that("calc_bubble_fc_() forces pairwise wilcox even with >2 groups when test = \"wilcox\"", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    low.CD4 = c(14.6, 14.9, 15.1, 15.4),
    high.CD4 = c(19.6, 19.9, 20.1, 20.4)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "wilcox", p_adjust_method = "none", summary_fun = mean, pseudocount = 0
  )
  expected_low_p <- stats::wilcox.test(
    data$value[data$group == "control"], data$value[data$group == "low"]
  )$p.value
  expect_equal(result$p_value[result$comparison == "low"], expected_low_p, tolerance = 1e-8)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-plot_bubble.R')"`
Expected: FAIL — the first two tests get `p_value` all `NA` (no `"kruskal"` branch yet in `calc_bubble_p_values_()`, which currently returns `NULL` for unmatched `resolved_test`, causing a `left_join()` error or all-NA `p_value`). The third test should already pass (pairwise wilcox unaffected).

- [ ] **Step 3: Add the kruskal branch**

In `R/plot_bubble.R`, extend `calc_bubble_p_values_()` by adding an `else if` branch after the closing `}` of the `if (resolved_test %in% c("wilcox", "t.test")) { ... }` block:

```r
  } else if (resolved_test == "kruskal") {
    population_levels <- unique(fc_table$population)
    purrr::map_dfr(population_levels, function(pop) {
      pop_data <- data[data[[population_col]] == pop, ]
      group_counts <- table(pop_data[[group_col]])
      testable_groups <- names(group_counts)[group_counts >= 2L]
      if (!(control %in% testable_groups) || length(testable_groups) < 2L) {
        return(tibble::tibble(
          population = character(), comparison = character(), p_value = double()
        ))
      }
      pop_data_testable <- pop_data[pop_data[[group_col]] %in% testable_groups, ]
      dunn <- rstatix::dunn_test(
        pop_data_testable,
        stats::as.formula(paste(value_col, "~", group_col)),
        p.adjust.method = p_adjust_method
      )
      dunn |>
        dplyr::filter(group1 == control | group2 == control) |>
        dplyr::mutate(
          comparison = ifelse(group1 == control, group2, group1),
          population = pop
        ) |>
        dplyr::select(population, comparison, p_value = p.adj)
    })
  }
```

(This becomes the full body of `calc_bubble_p_values_()`: the existing `if (resolved_test %in% c("wilcox", "t.test")) { ... }` block from Task 4, immediately followed by this `else if (resolved_test == "kruskal") { ... }` block.)

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-plot_bubble.R')"`
Expected: PASS, 16 tests.

- [ ] **Step 5: Commit**

```bash
git add R/plot_bubble.R tests/testthat/test-plot_bubble.R
git commit -m "feat: kruskal + Dunn post-hoc dispatch in calc_bubble_fc_()"
```

---

### Task 6: `plot_bubble_fc()` — the exported bubble plot

**Files:**
- Modify: `R/plot_bubble.R`
- Modify: `tests/testthat/test-plot_bubble.R`

**Interfaces:**
- Consumes: `calc_bubble_fc_(data, control, group_col, population_col, value_col, test, p_adjust_method, summary_fun, pseudocount)` from Tasks 2–5 (final form).
- Produces: exported `plot_bubble_fc(data, control, group_col = "group", population_col = "population", value_col = "value", test = "auto", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5)` returning a `ggplot` object.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-plot_bubble.R`:

```r
test_that("plot_bubble_fc() returns a ggplot object", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    treated.CD4 = c(19.6, 19.9, 20.1, 20.4)
  )
  p <- plot_bubble_fc(data, control = "control")
  expect_s3_class(p, "ggplot")
})

test_that("plot_bubble_fc() captions the control group", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    treated.CD4 = c(19.6, 19.9, 20.1, 20.4)
  )
  p <- plot_bubble_fc(data, control = "control")
  expect_match(p$labels$caption, "control")
})

test_that("plot_bubble_fc() builds for a 3-group (kruskal) comparison without error", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    low.CD4 = c(14.6, 14.9, 15.1, 15.4),
    high.CD4 = c(19.6, 19.9, 20.1, 20.4)
  )
  p <- plot_bubble_fc(data, control = "control")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot_bubble_fc() builds when some cells are untestable (NA p-value)", {
  data <- tibble::tibble(
    group = c("control", "treated", "treated", "treated"),
    population = c("CD3", "CD3", "CD3", "CD3"),
    value = c(10, 19.6, 19.9, 20.1)
  )
  p <- suppressWarnings(plot_bubble_fc(data, control = "control"))
  expect_no_error(ggplot2::ggplot_build(p))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-plot_bubble.R')"`
Expected: FAIL with `could not find function "plot_bubble_fc"`.

- [ ] **Step 3: Implement `plot_bubble_fc()`**

Add to the top of `R/plot_bubble.R` (before `calc_bubble_fc_()`):

```r
#' Bubble plot of group-vs-control fold-change across FACS populations
#'
#' @description
#' Compares one or more groups against an explicitly named control group,
#' one column per group-vs-control comparison and one row per population.
#' Each point is colored and sized by the log2 fold-change of the group's
#' mean value versus the control's, and annotated with a significance-star
#' summary from `test`. `data` must already be filtered to a single
#' `metric` and a single `tissue` (if those columns are present).
#'
#' @param data tibble already filtered to one metric and one tissue; must
#'   contain `group_col`, `population_col`, `value_col`.
#' @param control character; value in `group_col` treated as the
#'   baseline. Every other distinct value is compared against it.
#' @param group_col character; column in `data` holding group labels.
#'   Default `"group"`.
#' @param population_col character; column in `data` holding
#'   population names. Default `"population"`.
#' @param value_col character; column in `data` holding the numeric
#'   measurement. Default `"value"`.
#' @param test character; one of `"auto"` (default), `"wilcox"`,
#'   `"t.test"`, `"kruskal"`. `"auto"` uses pairwise
#'   `stats::wilcox.test()` when there are 2 groups (incl. control) and
#'   `stats::kruskal.test()` + `rstatix::dunn_test()` post-hoc when there
#'   are more than 2.
#' @param p_adjust_method character; passed to `stats::p.adjust()`
#'   (wilcox/t.test) or `rstatix::dunn_test()`'s own
#'   `p.adjust.method` (kruskal). Default `"BH"`.
#' @param summary_fun function; aggregator applied to each group's values
#'   when computing fold-change. Default `mean`.
#' @param pseudocount numeric; added to both summary values before
#'   dividing, avoiding `log2(0)`/`Inf`. Default `0.5`.
#'
#' @returns A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#'   dat <- facs_read_wsp("experiment.wsp", keywords = c("group"))$data |>
#'     dplyr::filter(metric == "FractionOfParent", tissue == "Spleen")
#'   plot_bubble_fc(dat, control = "control")
#' }
plot_bubble_fc <- function(data,
                            control,
                            group_col = "group",
                            population_col = "population",
                            value_col = "value",
                            test = "auto",
                            p_adjust_method = "BH",
                            summary_fun = mean,
                            pseudocount = 0.5) {

  fc_data <- calc_bubble_fc_(
    data, control, group_col, population_col, value_col,
    test, p_adjust_method, summary_fun, pseudocount
  )

  ggplot2::ggplot(fc_data, ggplot2::aes(x = comparison, y = population)) +
    ggplot2::geom_point(ggplot2::aes(
      color = log2fc, size = abs(log2fc), shape = is.na(p_value)
    )) +
    ggplot2::geom_text(ggplot2::aes(label = stars), vjust = -1.2, size = 3.5) +
    ggplot2::scale_color_gradient2(
      low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0,
      name = "log2FC"
    ) +
    ggplot2::scale_size_continuous(range = c(2, 10), name = "|log2FC|") +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
    ggplot2::labs(
      x = "Comparison",
      y = "Population",
      caption = glue::glue("Compared to control: {control}")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "#e1e0d9", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "#898781"),
      axis.title = ggplot2::element_text(color = "#0b0b0b"),
      plot.caption = ggplot2::element_text(color = "#52514e")
    )
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-plot_bubble.R')"`
Expected: PASS, 20 tests.

- [ ] **Step 5: Commit**

```bash
git add R/plot_bubble.R tests/testthat/test-plot_bubble.R
git commit -m "feat: add plot_bubble_fc() bubble plot"
```

---

### Task 7: Docs and verification

**Files:**
- Modify: `CLAUDE.md`
- Create/Modify (generated): `NAMESPACE`, `man/plot_bubble_fc.Rd`

**Interfaces:**
- Consumes: the finished `plot_bubble_fc()` from Task 6.
- Produces: exported symbol registered in `NAMESPACE`, `man/plot_bubble_fc.Rd` generated, `CLAUDE.md` documents the new data structure, `devtools::check()` passes with 0 errors/0 warnings.

- [ ] **Step 1: Add a Data structures entry to CLAUDE.md**

In `CLAUDE.md`'s "Data structures" section, after the `meta_annotate()` entry and before the closing `---` of that section, add:

```
### `plot_bubble_fc()` — group-vs-control fold-change bubble plot
- **Input**: long tibble already filtered to one `metric` and one `tissue`
  (if those columns are present); must contain a group column (default
  `group`), a population column (default `population`), and a numeric
  value column (default `value`).
- **Output**: a `ggplot` object (not a tibble) — a bubble plot with one
  column per non-control group and one row per population. Point color =
  signed log2 fold-change (group mean vs `control` mean, diverging
  blue/red), point size = `abs(log2fc)`, overlaid with significance stars
  (`*`/`**`/`***`) from `test` (`"auto"`: `wilcox.test` for 2 groups,
  `kruskal.test` + `rstatix::dunn_test()` post-hoc for >2). Caption names
  the `control` group.
```

- [ ] **Step 2: Regenerate documentation**

Run: `Rscript -e "devtools::document()"`
Expected: creates `man/plot_bubble_fc.Rd`, adds `export(plot_bubble_fc)` to `NAMESPACE`, no roxygen errors.

- [ ] **Step 3: Run the full test suite**

Run: `Rscript -e "devtools::test()"`
Expected: all tests pass, including the pre-existing suite (no regressions).

- [ ] **Step 4: Run `devtools::check()` and reconcile CLAUDE.md's "Known check output" section**

Run: `Rscript -e "devtools::check()"`
Expected: 0 errors, 0 warnings. Read the NOTEs output. If any new "R code for possible problems" NOTEs reference `plot_bubble_fc`, `calc_bubble_fc_`, `calc_bubble_p_values_`, or `adjust_p_` (expected: bare column names like `population`, `comparison`, `log2fc`, `p_value`, `testable`, `group1`, `group2`, `p.adj` inside `dplyr`/`tidyr` verbs — the same tidy-eval false-positive pattern as every other domain), add one bullet for them to CLAUDE.md's "Known check output" section, following the existing bullet style (e.g. the `facs_read_fcs_gated` bullet).

- [ ] **Step 5: Commit**

```bash
git add CLAUDE.md NAMESPACE man/plot_bubble_fc.Rd
git commit -m "docs: document plot_bubble_fc() and reconcile check output notes"
```
