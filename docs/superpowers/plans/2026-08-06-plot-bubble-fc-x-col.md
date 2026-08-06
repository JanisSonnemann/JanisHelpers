# `plot_bubble_fc()` x_col Argument Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an `x_col` argument to `plot_bubble_fc()`/`calc_bubble_fc_()` that lets the caller redirect the plot's x-axis to any column (typically `"tissue"` or `"metric"`) instead of the default group-vs-control `"comparison"` axis, when `group_col` has exactly 2 distinct values.

**Architecture:** Generalize the existing "one column per non-control group" logic in `R/plot_bubble.R` to "one column per distinct value of `x_col`". Internally, both `calc_bubble_fc_()` and `calc_bubble_p_values_()` switch from an `x_col == "comparison"` branch (unchanged existing behavior: slice-free, iterate `group_col` values) to an `x_col != "comparison"` branch (slice `data` by each `x_col` value first, then compare the single non-control group against control within that slice). The output tibble's x-axis column is renamed to whatever `x_col` is, so `plot_bubble_fc()` can build `aes(x = .data[[x_col]])` generically.

**Tech Stack:** R, `dplyr`, `tidyr`, `purrr`, `tibble`, `ggplot2`, `rstatix`, `stats`, `glue`, `tools` — all already Imports of this package, no new dependencies.

## Global Constraints

- Follow `CLAUDE.md`: base pipe `|>`, `pkg::fn()` namespacing for every external call inside exported/internal function bodies, `dplyr::select(!col)` not `-col`, non-ASCII escaped as `\uXXXX` in real R string literals only (not in `#'` roxygen prose — use plain ASCII there, per [[feedback_roxygen_non_ascii_escape]]).
- `devtools::check()` must stay at 0 errors, 0 warnings after every task that changes `R/plot_bubble.R`.
- Every exported function keeps `@param`/`@returns`/`@export`; `@examples` uses `\dontrun{}` (unchanged — `plot_bubble_fc()` already does).
- Do not modify the `x_col == "comparison"` (default) code path's *output* in any way — every existing test in `tests/testthat/test-plot_bubble.R` must keep passing unmodified.
- Spec: `docs/superpowers/specs/2026-08-06-plot-bubble-fc-x-col-design.md`.

---

### Task 1: `x_col` argument wiring + validation guards

**Files:**
- Modify: `R/plot_bubble.R` (signature of `plot_bubble_fc()` and `calc_bubble_fc_()`, plus `calc_bubble_fc_()`'s validation block)
- Test: `tests/testthat/test-plot_bubble.R`

**Interfaces:**
- Consumes: nothing new — this task only adds a parameter and guard clauses ahead of the existing (untouched) computation.
- Produces: `calc_bubble_fc_(data, control, group_col, population_col, value_col, x_col, test, p_adjust_method, summary_fun, pseudocount)` — note `x_col` is the new 6th positional parameter, inserted after `value_col` and before `test`. `plot_bubble_fc()` gains the same `x_col = "comparison"` parameter and threads it through positionally. `group_levels <- unique(data[[group_col]])` is now computed early in `calc_bubble_fc_()` (previously computed later, at what is currently line 108) — Task 2 will remove that now-redundant later computation.

- [ ] **Step 1: Write the failing tests**

Add to `tests/testthat/test-plot_bubble.R` (after the existing `test_that("calc_bubble_fc_() errors on an invalid test argument", ...)` block, before `make_data_ <- function(...) {...}`):

```r
test_that("calc_bubble_fc_() errors when x_col is set and group_col has more than 2 distinct values", {
  data <- tibble::tibble(
    group = c("control", "control", "low", "low", "high", "high"),
    population = "CD4",
    tissue = "Spleen",
    value = c(9.7, 9.9, 14.6, 14.9, 19.6, 19.9)
  )
  expect_error(
    calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value", x_col = "tissue",
      test = "auto", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5
    ),
    "exactly 2 distinct values"
  )
})

test_that("calc_bubble_fc_() errors when x_col is not present in data", {
  data <- tibble::tibble(
    group = c("control", "control", "treated", "treated"),
    population = "CD4",
    value = c(9.7, 9.9, 19.6, 19.9)
  )
  expect_error(
    calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value", x_col = "tissue",
      test = "auto", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5
    ),
    "not found in `data`"
  )
})

test_that("calc_bubble_fc_() errors when x_col equals group_col, population_col, or value_col", {
  data <- tibble::tibble(
    group = c("control", "control", "treated", "treated"),
    population = "CD4",
    value = c(9.7, 9.9, 19.6, 19.9)
  )
  expect_error(
    calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value", x_col = "group",
      test = "auto", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5
    ),
    "cannot be the same as"
  )
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-plot_bubble.R")'`
Expected: FAIL — `calc_bubble_fc_()` doesn't accept an `x_col` argument yet (`unused argument (x_col = "tissue")`).

- [ ] **Step 3: Implement `x_col` wiring and guards**

In `R/plot_bubble.R`, change `plot_bubble_fc()`'s signature from:

```r
plot_bubble_fc <- function(data,
                            control,
                            group_col = "group",
                            population_col = "population",
                            value_col = "value",
                            test = "auto",
                            p_adjust_method = "BH",
                            summary_fun = mean,
                            pseudocount = 0.5) {
```

to:

```r
plot_bubble_fc <- function(data,
                            control,
                            group_col = "group",
                            population_col = "population",
                            value_col = "value",
                            x_col = "comparison",
                            test = "auto",
                            p_adjust_method = "BH",
                            summary_fun = mean,
                            pseudocount = 0.5) {
```

Change the `calc_bubble_fc_()` call inside `plot_bubble_fc()` from:

```r
  fc_data <- calc_bubble_fc_(
    data, control, group_col, population_col, value_col,
    test, p_adjust_method, summary_fun, pseudocount
  )
```

to:

```r
  fc_data <- calc_bubble_fc_(
    data, control, group_col, population_col, value_col, x_col,
    test, p_adjust_method, summary_fun, pseudocount
  )
```

Change `calc_bubble_fc_()`'s signature from:

```r
calc_bubble_fc_ <- function(data, control, group_col, population_col,
                             value_col, test, p_adjust_method,
                             summary_fun, pseudocount) {
```

to:

```r
calc_bubble_fc_ <- function(data, control, group_col, population_col,
                             value_col, x_col, test, p_adjust_method,
                             summary_fun, pseudocount) {
```

Change `calc_bubble_fc_()`'s validation block from:

```r
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
```

to:

```r
  missing_cols <- setdiff(c(group_col, population_col, value_col), names(data))
  if (length(missing_cols) > 0L) {
    stop(glue::glue("`data` is missing column(s): {paste(missing_cols, collapse = ', ')}."))
  }
  if (!control %in% data[[group_col]]) {
    stop(glue::glue("`control` ('{control}') not found in `data${group_col}`."))
  }

  group_levels <- unique(data[[group_col]])

  if (x_col != "comparison") {
    if (x_col %in% c(group_col, population_col, value_col)) {
      stop(glue::glue(
        "`x_col` ('{x_col}') cannot be the same as `group_col`, ",
        "`population_col`, or `value_col`."
      ))
    }
    if (!x_col %in% names(data)) {
      stop(glue::glue("`x_col` ('{x_col}') not found in `data`."))
    }
    if (length(group_levels) != 2L) {
      stop(glue::glue(
        "`x_col` can only be set to a column other than \"comparison\" when ",
        "`group_col` has exactly 2 distinct values (control + one treatment); ",
        "found {length(group_levels)}. Use x_col = \"comparison\" (default) ",
        "for more groups."
      ))
    }
  }
  if ("metric" %in% names(data) && x_col != "metric" &&
      dplyr::n_distinct(data[["metric"]]) > 1L) {
    stop("`data` contains more than one distinct `metric`; pre-filter to a single metric.")
  }
  if ("tissue" %in% names(data) && x_col != "tissue" &&
      dplyr::n_distinct(data[["tissue"]]) > 1L) {
    stop("`data` contains more than one distinct `tissue`; pre-filter to a single tissue.")
  }
  if (!test %in% c("auto", "wilcox", "t.test", "kruskal")) {
    stop('`test` must be one of "auto", "wilcox", "t.test", "kruskal".')
  }
```

Leave the rest of `calc_bubble_fc_()` and all of `calc_bubble_p_values_()` untouched for now — the later `group_levels <- unique(data[[group_col]])` line further down the function (in the `population_levels`/`comparison_levels` block) becomes a harmless redundant recomputation until Task 2 removes it.

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-plot_bubble.R")'`
Expected: PASS (all tests, including the 3 new ones and every pre-existing test).

- [ ] **Step 5: Commit**

```bash
git add R/plot_bubble.R tests/testthat/test-plot_bubble.R
git commit -m "feat: add x_col argument wiring and validation to plot_bubble_fc()"
```

---

### Task 2: Generalize `fc_table` and the wilcox/t.test p-value path to `x_col`

**Files:**
- Modify: `R/plot_bubble.R` (`calc_bubble_fc_()`'s `population_levels`/`fc_table`/final-pipe block; `calc_bubble_p_values_()`'s signature and wilcox/t.test branch, plus a column-name-only rename in the kruskal branch)
- Test: `tests/testthat/test-plot_bubble.R`

**Interfaces:**
- Consumes: `x_col`/`group_levels` from Task 1.
- Produces: `calc_bubble_p_values_(data, control, group_col, population_col, value_col, x_col, comparison_group, resolved_test, p_adjust_method, fc_table)` — note the 2 new parameters `x_col`, `comparison_group` inserted before `resolved_test`. `calc_bubble_fc_()`'s returned tibble now has a column named after the literal value of `x_col` (e.g. `comparison` by default, or `tissue`/`metric`/etc.) instead of a hardcoded `comparison` column — later tasks and callers must read that column dynamically via `.data[[x_col]]`, not assume it's literally called `comparison`.

- [ ] **Step 1: Write the failing tests**

Add to `tests/testthat/test-plot_bubble.R`, after the existing `make_data_()`-based tests (a good spot: right after the `"calc_bubble_fc_() adjusts p-values per comparison column"` test block):

```r
test_that("calc_bubble_fc_() supports x_col = \"tissue\" for a 2-group comparison", {
  data <- tibble::tibble(
    group = c("control", "control", "treated", "treated",
              "control", "control", "treated", "treated"),
    population = "CD4",
    tissue = c("Spleen", "Spleen", "Spleen", "Spleen",
               "Blood", "Blood", "Blood", "Blood"),
    value = c(9.7, 10.3, 19.6, 20.4, 4.8, 5.2, 4.85, 5.15)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value", x_col = "tissue",
    test = "wilcox", p_adjust_method = "none", summary_fun = mean, pseudocount = 0
  )
  expect_true("tissue" %in% names(result))
  expect_false("comparison" %in% names(result))
  expect_equal(levels(result$tissue), c("Spleen", "Blood"))
  expect_equal(result$log2fc[result$tissue == "Spleen"], 1, tolerance = 1e-8)
  expect_equal(result$log2fc[result$tissue == "Blood"], 0, tolerance = 1e-8)
})

test_that("calc_bubble_fc_() supports x_col = \"metric\" for a 2-group comparison", {
  data <- tibble::tibble(
    group = c("control", "control", "treated", "treated",
              "control", "control", "treated", "treated"),
    population = "CD4",
    metric = c("Count", "Count", "Count", "Count",
               "FractionOfParent", "FractionOfParent", "FractionOfParent", "FractionOfParent"),
    value = c(9.7, 10.3, 19.6, 20.4, 4.8, 5.2, 4.85, 5.15)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value", x_col = "metric",
    test = "wilcox", p_adjust_method = "none", summary_fun = mean, pseudocount = 0
  )
  expect_true("metric" %in% names(result))
  expect_equal(levels(result$metric), c("Count", "FractionOfParent"))
  expect_equal(result$log2fc[result$metric == "Count"], 1, tolerance = 1e-8)
  expect_equal(result$log2fc[result$metric == "FractionOfParent"], 0, tolerance = 1e-8)
})

test_that("calc_bubble_fc_() lifts the tissue guard but keeps the metric guard when x_col = \"tissue\"", {
  data <- tibble::tibble(
    group = rep(c("control", "control", "treated", "treated"), 2),
    population = "CD4",
    tissue = rep(c("Spleen", "Blood"), each = 4),
    metric = c("Count", "Count", "Count", "Count",
               "FractionOfParent", "FractionOfParent", "FractionOfParent", "FractionOfParent"),
    value = c(9.7, 10.3, 19.6, 20.4, 4.8, 5.2, 4.85, 5.15)
  )
  expect_error(
    calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value", x_col = "tissue",
      test = "wilcox", p_adjust_method = "none", summary_fun = mean, pseudocount = 0
    ),
    "metric"
  )
})

test_that("calc_bubble_fc_() applies p_adjust_method per x_value when x_col is set", {
  data <- tibble::tibble(
    group = rep(c("control", "control", "treated", "treated"), 2),
    population = rep(c("CD4", "CD8"), each = 4),
    tissue = "Spleen",
    value = c(9.7, 10.3, 19.6, 20.4, 4.8, 5.2, 4.85, 5.15)
  )
  raw <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value", x_col = "tissue",
    test = "wilcox", p_adjust_method = "none", summary_fun = mean, pseudocount = 0
  )
  adjusted <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value", x_col = "tissue",
    test = "wilcox", p_adjust_method = "bonferroni", summary_fun = mean, pseudocount = 0
  )
  expect_equal(
    sort(adjusted$p_value),
    sort(stats::p.adjust(raw$p_value, method = "bonferroni")),
    tolerance = 1e-8
  )
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-plot_bubble.R")'`
Expected: FAIL — `result$tissue`/`result$metric` don't exist yet (output is still a `comparison` column regardless of `x_col`), or errors about unused `x_col` semantics.

- [ ] **Step 3: Implement the generalization**

In `R/plot_bubble.R`, change the block in `calc_bubble_fc_()` starting at `population_levels <- unique(data[[population_col]])` through the end of the function, from:

```r
  population_levels <- unique(data[[population_col]])
  group_levels <- unique(data[[group_col]])
  comparison_levels <- setdiff(group_levels, control)

  resolved_test <- test
  if (test == "auto") {
    resolved_test <- if (length(group_levels) > 2L) "kruskal" else "wilcox"
  }

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
}
```

to:

```r
  population_levels <- unique(data[[population_col]])

  if (x_col == "comparison") {
    x_values <- setdiff(group_levels, control)
    comparison_group <- NA_character_
  } else {
    x_values <- unique(data[[x_col]])
    comparison_group <- setdiff(group_levels, control)
  }

  resolved_test <- test
  if (test == "auto") {
    resolved_test <- if (length(group_levels) > 2L) "kruskal" else "wilcox"
  }

  fc_table <- tidyr::expand_grid(
    population = population_levels,
    x_value = x_values
  ) |>
    purrr::pmap_dfr(function(population, x_value) {
      if (x_col == "comparison") {
        slice <- data
        comparison <- x_value
      } else {
        slice <- data[data[[x_col]] == x_value, ]
        comparison <- comparison_group
      }
      control_values <- slice[[value_col]][
        slice[[population_col]] == population & slice[[group_col]] == control
      ]
      group_values <- slice[[value_col]][
        slice[[population_col]] == population & slice[[group_col]] == comparison
      ]
      log2fc <- log2(
        (summary_fun(group_values) + pseudocount) / (summary_fun(control_values) + pseudocount)
      )
      testable <- length(control_values) >= 2L && length(group_values) >= 2L
      if (!testable) {
        warning(glue::glue(
          "Fewer than 2 observations for population '{population}', ",
          "{x_col} '{x_value}' vs control '{control}'; p-value set to NA."
        ))
      }
      tibble::tibble(
        population = population, x_value = x_value,
        log2fc = log2fc, testable = testable
      )
    })

  p_values <- calc_bubble_p_values_(
    data, control, group_col, population_col, value_col, x_col, comparison_group,
    resolved_test, p_adjust_method, fc_table
  )

  fc_table |>
    dplyr::select(!testable) |>
    dplyr::left_join(p_values, by = c("population", "x_value")) |>
    dplyr::mutate(
      stars = dplyr::case_when(
        is.na(p_value) ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE ~ ""
      ),
      population = factor(population, levels = population_levels),
      x_value = factor(x_value, levels = x_values)
    ) |>
    dplyr::rename(!!x_col := x_value)
}
```

Then update `calc_bubble_p_values_()`. Change its signature and the `wilcox`/`t.test` branch from:

```r
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
  } else if (resolved_test == "kruskal") {
```

to:

```r
calc_bubble_p_values_ <- function(data, control, group_col, population_col,
                                   value_col, x_col, comparison_group,
                                   resolved_test, p_adjust_method, fc_table) {
  if (resolved_test %in% c("wilcox", "t.test")) {
    raw <- fc_table |>
      purrr::pmap_dfr(function(population, x_value, log2fc, testable) {
        p_value <- NA_real_
        if (testable) {
          if (x_col == "comparison") {
            slice <- data
            comparison <- x_value
          } else {
            slice <- data[data[[x_col]] == x_value, ]
            comparison <- comparison_group
          }
          control_values <- slice[[value_col]][
            slice[[population_col]] == population & slice[[group_col]] == control
          ]
          group_values <- slice[[value_col]][
            slice[[population_col]] == population & slice[[group_col]] == comparison
          ]
          p_value <- if (resolved_test == "wilcox") {
            stats::wilcox.test(control_values, group_values)$p.value
          } else {
            stats::t.test(control_values, group_values)$p.value
          }
        }
        tibble::tibble(population = population, x_value = x_value, p_value = p_value)
      })

    raw |>
      dplyr::group_by(x_value) |>
      dplyr::mutate(p_value = adjust_p_(p_value, p_adjust_method)) |>
      dplyr::ungroup()
  } else if (resolved_test == "kruskal") {
```

Finally, rename the two `comparison` occurrences in the (still `x_col`-unaware — Task 3 handles that) kruskal branch to `x_value`, so its output columns match the new `by = c("population", "x_value")` join key. Change:

```r
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
```

to:

```r
        return(tibble::tibble(
          population = character(), x_value = character(), p_value = double()
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
          x_value = ifelse(group1 == control, group2, group1),
          population = pop
        ) |>
        dplyr::select(population, x_value, p_value = p.adj)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-plot_bubble.R")'`
Expected: PASS — all new tests plus every pre-existing test (the default `x_col = "comparison"` path is untouched in output shape, just renamed internally and then renamed back via `dplyr::rename(!!x_col := x_value)`).

- [ ] **Step 5: Commit**

```bash
git add R/plot_bubble.R tests/testthat/test-plot_bubble.R
git commit -m "feat: generalize calc_bubble_fc_() fold-change and wilcox/t.test paths to x_col"
```

---

### Task 3: Generalize the kruskal-forced path to `x_col`

**Files:**
- Modify: `R/plot_bubble.R` (`calc_bubble_p_values_()`'s kruskal branch)
- Test: `tests/testthat/test-plot_bubble.R`

**Interfaces:**
- Consumes: `x_col`, `comparison_group` params from Task 2.
- Produces: no new interface — same `calc_bubble_p_values_()` signature as Task 2, now fully `x_col`-aware in every branch.

- [ ] **Step 1: Write the failing test**

Add to `tests/testthat/test-plot_bubble.R`, after the `"calc_bubble_fc_() applies p_adjust_method per x_value when x_col is set"` test from Task 2:

```r
test_that("calc_bubble_fc_() forces kruskal+dunn per x_value when x_col is set", {
  data <- tibble::tibble(
    group = rep(c("control", "control", "treated", "treated"), 2),
    population = "CD4",
    tissue = rep(c("Spleen", "Blood"), each = 4),
    value = c(9.7, 9.9, 19.6, 19.9, 4.8, 4.9, 4.85, 4.95)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value", x_col = "tissue",
    test = "kruskal", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0
  )
  expect_equal(levels(result$tissue), c("Spleen", "Blood"))
  expect_true(all(!is.na(result$p_value)))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-plot_bubble.R")'`
Expected: FAIL — the kruskal branch ignores `x_col`/`tissue` slicing and runs `dunn_test()` on the full (both-tissue) data per population, producing wrong/`NA` p-values or a row-count mismatch against `fc_table`.

- [ ] **Step 3: Implement the kruskal-branch generalization**

In `R/plot_bubble.R`, change the `else if (resolved_test == "kruskal")` branch of `calc_bubble_p_values_()` from:

```r
  } else if (resolved_test == "kruskal") {
    population_levels <- unique(fc_table$population)
    purrr::map_dfr(population_levels, function(pop) {
      pop_data <- data[data[[population_col]] == pop, ]
      group_counts <- table(pop_data[[group_col]])
      testable_groups <- names(group_counts)[group_counts >= 2L]
      if (!(control %in% testable_groups) || length(testable_groups) < 2L) {
        return(tibble::tibble(
          population = character(), x_value = character(), p_value = double()
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
          x_value = ifelse(group1 == control, group2, group1),
          population = pop
        ) |>
        dplyr::select(population, x_value, p_value = p.adj)
    })
  }
}
```

to:

```r
  } else if (resolved_test == "kruskal") {
    population_levels <- unique(fc_table$population)

    if (x_col == "comparison") {
      purrr::map_dfr(population_levels, function(pop) {
        pop_data <- data[data[[population_col]] == pop, ]
        group_counts <- table(pop_data[[group_col]])
        testable_groups <- names(group_counts)[group_counts >= 2L]
        if (!(control %in% testable_groups) || length(testable_groups) < 2L) {
          return(tibble::tibble(
            population = character(), x_value = character(), p_value = double()
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
            x_value = ifelse(group1 == control, group2, group1),
            population = pop
          ) |>
          dplyr::select(population, x_value, p_value = p.adj)
      })
    } else {
      x_values <- unique(fc_table$x_value)
      purrr::map_dfr(x_values, function(x_value) {
        slice <- data[data[[x_col]] == x_value, ]
        purrr::map_dfr(population_levels, function(pop) {
          pop_data <- slice[slice[[population_col]] == pop, ]
          group_counts <- table(pop_data[[group_col]])
          testable_groups <- names(group_counts)[group_counts >= 2L]
          if (!(control %in% testable_groups) || length(testable_groups) < 2L) {
            return(tibble::tibble(
              population = character(), x_value = character(), p_value = double()
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
            dplyr::mutate(population = pop, x_value = x_value) |>
            dplyr::select(population, x_value, p_value = p.adj)
        })
      })
    }
  }
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-plot_bubble.R")'`
Expected: PASS — all tests including the new kruskal/`x_col` test.

- [ ] **Step 5: Commit**

```bash
git add R/plot_bubble.R tests/testthat/test-plot_bubble.R
git commit -m "feat: generalize calc_bubble_p_values_() kruskal path to x_col"
```

---

### Task 4: Generalize `plot_bubble_fc()`'s rendering (aes, axis label, caption)

**Files:**
- Modify: `R/plot_bubble.R` (`plot_bubble_fc()`'s ggplot-building body)
- Test: `tests/testthat/test-plot_bubble.R`

**Interfaces:**
- Consumes: `fc_data` (from `calc_bubble_fc_()`, Task 2/3) with its x-axis column named after `x_col`; `x_col`, `data`, `group_col`, `control` (already in scope in `plot_bubble_fc()`).
- Produces: a `ggplot` object whose x aesthetic, x-axis label, and caption reflect `x_col`. No new exported interface.

- [ ] **Step 1: Write the failing test**

Add to `tests/testthat/test-plot_bubble.R`, after the existing `"plot_bubble_fc() builds when some cells are untestable (NA p-value)"` test at the end of the file:

```r
test_that("plot_bubble_fc() builds with x_col = \"tissue\" and captions the single comparison group", {
  data <- tibble::tibble(
    group = c("control", "control", "treated", "treated",
              "control", "control", "treated", "treated"),
    population = "CD4",
    tissue = c("Spleen", "Spleen", "Spleen", "Spleen",
               "Blood", "Blood", "Blood", "Blood"),
    value = c(9.7, 10.3, 19.6, 20.4, 4.8, 5.2, 4.85, 5.15)
  )
  p <- plot_bubble_fc(data, control = "control", x_col = "tissue")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_match(p$labels$caption, "treated")
  expect_match(p$labels$caption, "control")
  expect_equal(p$labels$x, "Tissue")
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-plot_bubble.R")'`
Expected: FAIL — `ggplot2::aes(x = comparison, ...)` looks up a literal `comparison` column, which no longer exists in `fc_data` when `x_col = "tissue"` (it's now named `tissue`); build errors with "object 'comparison' not found" (or similar), and `p$labels$x`/`caption` don't reflect `x_col`.

- [ ] **Step 3: Implement the rendering generalization**

In `R/plot_bubble.R`, change the plot-building body of `plot_bubble_fc()` from:

```r
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
```

to:

```r
  x_label <- if (x_col == "comparison") "Comparison" else tools::toTitleCase(x_col)
  caption <- if (x_col == "comparison") {
    glue::glue("Compared to control: {control}")
  } else {
    comparison_group <- setdiff(unique(data[[group_col]]), control)
    glue::glue("{comparison_group} vs control: {control}")
  }

  ggplot2::ggplot(fc_data, ggplot2::aes(x = .data[[x_col]], y = population)) +
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
      x = x_label,
      y = "Population",
      caption = caption
    ) +
    ggplot2::theme_minimal() +
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-plot_bubble.R")'`
Expected: PASS — every test in the file, including all pre-existing `x_col = "comparison"` (default) plot tests, which are unaffected since `x_label`/`caption` reduce to their original values in that branch.

- [ ] **Step 5: Run the full test suite and `devtools::check()`**

Run: `Rscript -e 'devtools::load_all(); devtools::test()'`
Expected: 0 failures across the whole package (not just `test-plot_bubble.R`).

Run: `Rscript -e 'devtools::document(); devtools::check()'`
Expected: 0 errors, 0 warnings. Note any new NOTE text for `x_col`/`x_value`/`comparison_group` — Task 5 documents it in `CLAUDE.md`.

- [ ] **Step 6: Commit**

```bash
git add R/plot_bubble.R tests/testthat/test-plot_bubble.R
git commit -m "feat: generalize plot_bubble_fc() aes/label/caption to x_col"
```

---

### Task 5: Documentation — roxygen and `CLAUDE.md`

**Files:**
- Modify: `R/plot_bubble.R` (roxygen block above `plot_bubble_fc()`)
- Modify: `man/plot_bubble_fc.Rd` (regenerated, not hand-edited)
- Modify: `CLAUDE.md`

**Interfaces:**
- Consumes: final behavior of `plot_bubble_fc()`/`calc_bubble_fc_()` from Tasks 1-4.
- Produces: nothing consumed by other tasks — this is the terminal documentation task.

- [ ] **Step 1: Update the roxygen block**

In `R/plot_bubble.R`, change the `@description` paragraph and add an `@param x_col` entry. Change:

```r
#' @description
#' Compares one or more groups against an explicitly named control group,
#' one column per group-vs-control comparison and one row per population.
#' Each point is colored and sized by the log2 fold-change of the group's
#' mean value versus the control's, and annotated with a significance-star
#' summary from `test`. `data` must already be filtered to a single
#' `metric` and a single `tissue` (if those columns are present).
```

to:

```r
#' @description
#' Compares one or more groups against an explicitly named control group,
#' one column per group-vs-control comparison and one row per population.
#' Each point is colored and sized by the log2 fold-change of the group's
#' mean value versus the control's, and annotated with a significance-star
#' summary from `test`. `data` must already be filtered to a single
#' `metric` and a single `tissue` (if those columns are present), unless
#' `x_col` is set to `"metric"`/`"tissue"` (see `x_col`).
```

Change:

```r
#' @param value_col character; column in `data` holding the numeric
#'   measurement. Default `"value"`.
#' @param test character; one of `"auto"` (default), `"wilcox"`,
```

to:

```r
#' @param value_col character; column in `data` holding the numeric
#'   measurement. Default `"value"`.
#' @param x_col character; which dimension varies across the plot's
#'   x-axis. Default `"comparison"`: one column per non-control group, as
#'   today, and works with any number of groups. Any other column name
#'   present in `data` (typically `"tissue"` or `"metric"`, but not
#'   restricted to those): one column per distinct value of that column
#'   instead, holding the single non-control group fixed. Only valid when
#'   `group_col` has exactly 2 distinct values (control + one treatment).
#' @param test character; one of `"auto"` (default), `"wilcox"`,
```

- [ ] **Step 2: Regenerate `man/plot_bubble_fc.Rd`**

Run: `Rscript -e 'devtools::document()'`
Expected: `man/plot_bubble_fc.Rd` is updated to include the new `\item{x_col}{...}` entry; no roxygen warnings.

- [ ] **Step 3: Update `CLAUDE.md`**

In `CLAUDE.md`, under `### plot_bubble_fc() — group-vs-control fold-change bubble plot`, add a sentence about `x_col` to the existing description paragraph — change:

```markdown
- **Output**: a `ggplot` object (not a tibble) — a bubble plot with one
  column per non-control group and one row per population. Point color =
```

to:

```markdown
- **Output**: a `ggplot` object (not a tibble) — a bubble plot with one
  column per non-control group (or, if `x_col` is set to a column other
  than `"comparison"` and `group_col` has exactly 2 distinct values, one
  column per distinct value of that column instead — e.g. `x_col =
  "tissue"` to compare one treatment vs control across multiple tissues
  in a single plot) and one row per population. Point color =
```

Then run `devtools::check()` (Task 4, Step 5 already ran it once) and compare its NOTE output for `R/plot_bubble.R` against the existing "Known check output" bullet for `plot_bubble_fc`/`calc_bubble_fc_`/`calc_bubble_p_values_` in `CLAUDE.md`. The variable names that bullet lists (`testable`, `population`, `comparison`, `log2fc`, `p_value`, `group1`, `group2`, `p.adj`, `stars`) should now read `x_value` in place of `comparison` (since that's the renamed bare-column reference used throughout the rewritten code) — update the bullet text to match the actual `R CMD check` output verbatim. `x_col` and `comparison_group` should **not** appear in the NOTE list (they're formal arguments/plain local variables resolved by codetools, not bare dplyr-verb column references) — if they unexpectedly do appear, investigate before editing `CLAUDE.md` (that would indicate a bug, e.g. a bare column reference that should have used `.data[[x_col]]`).

- [ ] **Step 4: Run `devtools::check()` one final time**

Run: `Rscript -e 'devtools::check()'`
Expected: 0 errors, 0 warnings, and the NOTE list matches what `CLAUDE.md` now documents.

- [ ] **Step 5: Commit**

```bash
git add R/plot_bubble.R man/plot_bubble_fc.Rd CLAUDE.md
git commit -m "docs: document x_col argument for plot_bubble_fc()"
```

---

## Self-Review Notes

- **Spec coverage:** API shape (Task 1), mechanics/slicing generalization for wilcox/t.test (Task 2), kruskal path (Task 3), plot/labeling (Task 4), error handling (Tasks 1-2, folded into their guard clauses rather than a separate task since each guard belongs with the code path it protects), testing (test blocks embedded in every task), roxygen/CLAUDE.md (Task 5). All spec sections have a corresponding task.
- **Type consistency:** `calc_bubble_fc_()`'s signature (`data, control, group_col, population_col, value_col, x_col, test, p_adjust_method, summary_fun, pseudocount`) is introduced in Task 1 and never changes shape afterward. `calc_bubble_p_values_()`'s signature (`data, control, group_col, population_col, value_col, x_col, comparison_group, resolved_test, p_adjust_method, fc_table`) is introduced whole in Task 2 (Task 3 only edits its kruskal branch's body, not its signature). The output tibble's x-axis column name (`x_col`'s value) is consistent from Task 2 (where it's produced) through Task 4 (where `plot_bubble_fc()` consumes it via `.data[[x_col]]`).
