# Population order, subjects query, plot helpers -- Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Store each FACS panel's population display order in the DB, add a once-per-doc subjects query, make the restim calc functions run on DB output, and move the analysis doc's generic plot helpers into the package.

**Architecture:** One additive column (`facs_measurements.population_order`) written by the existing export path and returned by `db_query_facs()`; `facs_order_populations()` rebuilds the factor after `collect()`. `db_query_subjects()` is a new lazy query; meta is joined per dataset with the existing `meta_annotate()`. Plot helpers are new `plot_` functions in `R/plot_box.R`.

**Tech Stack:** R package (roxygen2, testthat 3e), DuckDB via DBI/dbplyr, ggplot2, ggpubr.

**Spec:** `docs/superpowers/specs/2026-11-19-population-order-plot-helpers-design.md`

## Global Constraints

- Follow `CLAUDE.md`: `|>` pipe; every external call namespaced `pkg::fn()`; `dplyr::select(!col)` not `-col`; non-ASCII only as `\uXXXX` in R strings, plain ASCII in roxygen.
- Every exported function: `@param`, `@returns` (non-empty), `@export`; `\dontrun{}` for examples needing a DB or files.
- Schema change is additive: an existing DB must open with `db_connect()` without recreating; `db_old_schema_()` must not reject a DB lacking `population_order`.
- `subjects.treatment_group` is NOT renamed in the DB; only `db_query_subjects()` output calls it `group`.
- No back-fill of `population_order` outside `db_export_facs()`.
- `ggpubr` and `scales` go into `Imports` (DESCRIPTION, alphabetical position, and the CLAUDE.md Imports list). `scales` is a ggplot2 dependency already installed; it must be declared because `plot_box()` calls `scales::percent` (my addition -- the spec only names ggpubr).
- Run R with the library that has duckdb/devtools/ggpubr (the worktree's renv lib is empty):
  `export R_LIBS=~/Library/Caches/org.R-project.R/R/renv/library/facs-pipeline-workflow-0ff3068e/macos/R-4.6/aarch64-apple-darwin23`
  then `Rscript --vanilla -e 'devtools::test(filter = "<file>")'`.
- Final `devtools::check()`: 0 errors, 0 warnings; only the NOTEs listed in CLAUDE.md.

## Review Focus

1. A prep whose `population` is character (no `facs_select_populations()`) must export with `population_order` `NULL`, not error -- Task 2 test.
2. Re-exporting a root that was stored before the column existed (old rows `NULL`) must leave only non-`NULL` orders for that root -- Task 2 test.
3. Filtering `db_query_subjects()` lazily on `group` (SQL reserved word) must work before `collect()` -- Task 3 test.
4. `facs_calc_log2fc()` on DB output carries columns that vary per population (`population_order`, `population_full_path`) -- they must be dropped, not break the pivot -- Task 5 test.
5. `facs_order_populations()` on data missing `population_order` / `wsp_root` (e.g. not from `db_query_facs()`) must error naming the column -- Task 4 test.

---

### Task 1: Schema -- `population_order` column

**Files:**
- Modify: `R/db_connect.R` (`create_schema_()`, `facs_measurements` block ~L152-168; roxygen of `db_connect()`)
- Test: `tests/testthat/test-db_connect.R`

**Interfaces:**
- Produces: `facs_measurements.population_order INTEGER` (nullable) in every DB opened by `db_connect()`.

- [ ] **Step 1: Write the failing tests**

```r
test_that("facs_measurements has a nullable population_order column", {
  con <- local_test_db()
  expect_true("population_order" %in% DBI::dbListFields(con, "facs_measurements"))
})

test_that("db_connect adds population_order to an existing database, keeping rows", {
  path <- tempfile(fileext = ".duckdb")
  con <- db_connect(path)
  DBI::dbExecute(con, "ALTER TABLE facs_measurements DROP COLUMN population_order")
  DBI::dbExecute(con, "INSERT INTO facs_measurements
    (facs_stain_id, population_full_path, population, metric, value)
    VALUES (1, 'root/CD4+', 'CD4+', 'count', 10)")
  DBI::dbDisconnect(con, shutdown = TRUE)

  con <- db_connect(path)
  withr::defer(DBI::dbDisconnect(con, shutdown = TRUE))
  rows <- DBI::dbGetQuery(con, "SELECT population, population_order FROM facs_measurements")
  expect_equal(rows$population, "CD4+")
  expect_true(is.na(rows$population_order))
})
```

- [ ] **Step 2: Run to verify they fail**

Run: `Rscript --vanilla -e 'devtools::test(filter = "db_connect")'`
Expected: both new tests FAIL (column absent).

- [ ] **Step 3: Implement**

Add `population_order INTEGER,` after `value` in the `CREATE TABLE`, and after it
`DBI::dbExecute(con, "ALTER TABLE facs_measurements ADD COLUMN IF NOT EXISTS population_order INTEGER")`.
Document the column in `db_connect()`'s roxygen (one sentence: added on open if missing).

- [ ] **Step 4: Run to verify they pass**

Run: `Rscript --vanilla -e 'devtools::test(filter = "db_connect")'`
Expected: all PASS (incl. the existing "a database with data can be snapshotted" test).

- [ ] **Step 5: Commit**

```bash
git add R/db_connect.R tests/testthat/test-db_connect.R man/
git commit -m "feat(db): additive facs_measurements.population_order column"
```

---

### Task 2: Export writes the order, `db_query_facs()` returns it

**Files:**
- Modify: `R/db_write_measurements.R` (`insert_facs_stain_()` ~L304-346; `db_write_facs()` ~L285-290)
- Modify: `R/db_export.R` (only if the per-stain data passed to `insert_facs_stain_()` loses the factor; today it is `dplyr::filter(data, ...)`, which keeps it)
- Modify: `R/db_query.R` (`db_query_facs()` select list + roxygen)
- Test: `tests/testthat/test-db_export.R`, `tests/testthat/test-db_query.R`

**Interfaces:**
- Consumes: Task 1 column.
- Produces: `insert_facs_stain_(con, sample_id, assay_id, data, ...)` -- `data$population` may be factor or character; when factor, `population_order = as.integer(data$population)`, else `NA_integer_`. `db_query_facs()` returns column `population_order` (integer, after `population`).

- [ ] **Step 1: Write the failing tests**

In `test-db_export.R` (fixture `fake_prep_data()` has levels `live_CD45+`, `CD4+`):

```r
test_that("export stores each population's display position", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep <- suppressMessages(local_loaded_prep(d, two_mice))
  suppressMessages(db_export_facs(con, prep, data_dir = d$data_dir))
  rows <- DBI::dbGetQuery(con,
    "SELECT DISTINCT population, population_order FROM facs_measurements ORDER BY 2")
  expect_equal(rows$population, c("live_CD45+", "CD4+"))
  expect_equal(rows$population_order, c(1L, 2L))
})

test_that("character population exports with NULL population_order", {
  # build the prep with population as character (as.character() before saving)
  # via the same helpers; expect all(is.na(population_order)) and 8 rows.
})

test_that("re-exporting a root stored without order fills population_order", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep <- suppressMessages(local_loaded_prep(d, two_mice))
  suppressMessages(db_export_facs(con, prep, data_dir = d$data_dir))
  DBI::dbExecute(con, "UPDATE facs_measurements SET population_order = NULL")
  suppressMessages(db_export_facs(con, prep, data_dir = d$data_dir, overwrite = TRUE))
  expect_equal(
    DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM facs_measurements WHERE population_order IS NULL")$n,
    0
  )
})
```

The character-population test: read how `local_loaded_prep()` (helper-prep.R) builds the prep and pass a data tibble whose `population` is `as.character()`; assert `all(is.na(...))` and `n_rows(con, "facs_measurements") == 8`.

In `test-db_query.R`, extend the first test: `expect_true(is.na(result$population_order))` (written via `db_write_facs()` with character `Population`), and add:

```r
test_that("db_write_facs stores population_order when Population is a factor", {
  # same call as the first test but Population = factor("CD4+", levels = c("CD3+", "CD4+"))
  # expect db_query_facs(con) |> dplyr::collect() |> dplyr::pull(population_order) == 2L
})
```

- [ ] **Step 2: Run to verify they fail**

Run: `Rscript --vanilla -e 'devtools::test(filter = "db_export|db_query")'`
Expected: new tests FAIL (column missing from insert / select).

- [ ] **Step 3: Implement**

`insert_facs_stain_()`: add `population_order` to `to_write` (computed from `data$population` before `as.character()`), and to the `INSERT ... SELECT` column lists. `db_write_facs()`: pass `population = data$Population` unchanged (factor kept) into the tibble so the same rule applies. `db_query_facs()`: add `population_order` after `population` in `select()`; mention it in roxygen ("position in the panel's display order; `NA` for imports before this column existed -- see [facs_order_populations()]").

- [ ] **Step 4: Run to verify they pass**

Run: `Rscript --vanilla -e 'devtools::test(filter = "db_export|db_query|db_write_measurements")'`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add R/db_write_measurements.R R/db_export.R R/db_query.R tests/testthat/test-db_export.R tests/testthat/test-db_query.R man/
git commit -m "feat(db): export and query FACS population display order"
```

---

### Task 3: `db_query_subjects()`

**Files:**
- Modify: `R/db_query.R`
- Test: `tests/testthat/test-db_query.R`

**Interfaces:**
- Produces: `db_query_subjects(con)` -> lazy `tbl` with columns, in order: `experiment_code`, `mouse_id`, `group`, `mouse_treatment`, `sex`, `cage`, `mouse_strain`, `generation`, `dob`, `start_date`, `bmt_date`, `end_date`. `group` = `subjects.treatment_group`.

- [ ] **Step 1: Write the failing test**

```r
test_that("db_query_subjects returns one row per mouse with group", {
  con <- local_test_db()
  db_write_experiment(con, experiment_code = "25-7")
  db_write_subjects(con, tibble::tibble(
    mouse_id = c("25-7-1", "25-7-2"), experiment_code = "25-7",
    treatment_group = c("irrel", "MPO409"), end_date = as.Date("2026-06-10")
  ))

  lazy <- db_query_subjects(con)
  expect_s3_class(lazy, "tbl_lazy")
  irrel <- lazy |> dplyr::filter(group == "irrel") |> dplyr::collect()
  expect_equal(irrel$mouse_id, "25-7-1")

  all <- dplyr::collect(lazy)
  expect_equal(names(all), c("experiment_code", "mouse_id", "group",
    "mouse_treatment", "sex", "cage", "mouse_strain", "generation", "dob",
    "start_date", "bmt_date", "end_date"))
  expect_false("treatment_group" %in% names(all))
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `Rscript --vanilla -e 'devtools::test(filter = "db_query")'`
Expected: FAIL, `could not find function "db_query_subjects"`.

- [ ] **Step 3: Implement `db_query_subjects(con)` in `R/db_query.R`**

`subjects` left-join `experiments` by `experiment_id`, `dplyr::select(..., group = treatment_group, ...)`. Roxygen `@details` describes the workflow from spec decision 3: connect read-only with `DBI::dbConnect(duckdb::duckdb(), path, read_only = TRUE)`, query subjects once, filter to the experiment, `collect()`, set group levels once, then join needed columns per dataset with `meta_annotate(data, subjects_cols, by = c("experiment_code", "mouse_id"))`. Put that as the `\dontrun{}` example.

- [ ] **Step 4: Run to verify it passes**

Run: `Rscript --vanilla -e 'devtools::test(filter = "db_query")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/db_query.R tests/testthat/test-db_query.R man/ NAMESPACE
git commit -m "feat(db): db_query_subjects() -- subject meta once, group column"
```

---

### Task 4: `facs_order_populations()`

**Files:**
- Modify: `R/facs_calc.R` (new exported function at the end)
- Test: `tests/testthat/test-facs_calc.R`

**Interfaces:**
- Consumes: columns `population` (chr), `population_order` (int), `wsp_root` (chr) from collected `db_query_facs()` output (Task 2).
- Produces: `facs_order_populations(data) -> data` with `population` as factor; nothing else changes.

Algorithm (decided by the spec, not by the tests alone): sort distinct `wsp_root`, `population`, `population_order` by `wsp_root` (alphabetical) then `population_order`; ordered rows first, `NA`-order rows last sorted by `population`; levels = `unique(population)` of that sequence (first occurrence wins). If any `NA` order: `message()` naming the affected `wsp_root`s.

- [ ] **Step 1: Write the failing tests**

```r
db_shaped <- function(root, pops, order) {
  tibble::tibble(wsp_root = root, population = pops, population_order = order,
                 metric = "count", value = 1)
}

test_that("facs_order_populations orders one panel by population_order", {
  d <- db_shaped("ov_26-1", c("CD8+", "live", "CD4+"), c(3L, 1L, 2L))
  out <- facs_order_populations(d)
  expect_equal(levels(out$population), c("live", "CD4+", "CD8+"))
  expect_equal(as.character(out$population), d$population)
})

test_that("a label shared by two panels keeps the first panel's position", {
  d <- dplyr::bind_rows(
    db_shaped("a_panel", c("live", "PMN", "CD4+"), 1:3),
    db_shaped("b_panel", c("PMN", "Tc"), 1:2)
  )
  expect_no_error(out <- facs_order_populations(d))
  expect_equal(levels(out$population), c("live", "PMN", "CD4+", "Tc"))
})

test_that("NA orders go last alphabetically with a message", {
  d <- dplyr::bind_rows(
    db_shaped("a_panel", c("live", "CD4+"), 1:2),
    db_shaped("old_panel", c("Zeta", "Beta"), NA_integer_)
  )
  expect_message(out <- facs_order_populations(d), "old_panel")
  expect_equal(levels(out$population), c("live", "CD4+", "Beta", "Zeta"))
})

test_that("facs_order_populations errors naming a missing column", {
  d <- db_shaped("a", "live", 1L) |> dplyr::select(!wsp_root)
  expect_error(facs_order_populations(d), "wsp_root")
})
```

- [ ] **Step 2: Run to verify they fail**

Run: `Rscript --vanilla -e 'devtools::test(filter = "facs_calc")'`
Expected: FAIL, function not found.

- [ ] **Step 3: Implement `facs_order_populations(data)` in `R/facs_calc.R`** (algorithm above; roxygen notes it is for collected `db_query_facs()` output and explains the multi-panel rule).

- [ ] **Step 4: Run to verify they pass**

Run: `Rscript --vanilla -e 'devtools::test(filter = "facs_calc")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/facs_calc.R tests/testthat/test-facs_calc.R man/ NAMESPACE
git commit -m "feat(facs): facs_order_populations() rebuilds display order after a DB query"
```

---

### Task 5: `mouse_col` for `facs_calc_log2fc()` / `facs_calc_diff()`

**Files:**
- Modify: `R/facs_calc.R` (`calc_restim_proportions_()` ~L264-340, `facs_calc_log2fc()`, `facs_calc_diff()` + roxygen)
- Test: `tests/testthat/test-facs_calc.R`

**Interfaces:**
- Produces: `facs_calc_log2fc(data, ref_pop, restim_col = "restimulation", ref_level = "unstim", pseudocount = 0.5, mouse_col = "mouse_ID")`; `facs_calc_diff(data, ref_pop, restim_col = "restimulation", ref_level = "unstim", pseudocount = 0, mouse_col = "mouse_ID")`. Output keeps the mouse column under its input name.

- [ ] **Step 1: Write the failing tests**

```r
db_log2fc_data <- log2fc_data |>
  dplyr::rename(mouse_id = mouse_ID, condition = restimulation) |>
  dplyr::select(!file_name) |>
  dplyr::mutate(
    experiment_code = "26-1", collected_at = as.Date("2026-06-16"),
    assay_name = "mpo-restim", wsp_root = "MPO-Restim_26-1",
    condition_keyword = "restimulation", source_file = "x.wsp",
    imported_at = as.POSIXct("2026-11-19 08:00:00"),
    population_order = ifelse(population == "CD3+", 1L, 2L),
    group = "irrel"
  )

test_that("facs_calc_log2fc runs on DB-shaped data with mouse_col/restim_col", {
  ref <- facs_calc_log2fc(log2fc_data, ref_pop = "CD3+")
  out <- facs_calc_log2fc(db_log2fc_data, ref_pop = "CD3+",
                          mouse_col = "mouse_id", restim_col = "condition")
  expect_equal(out$log2fc, ref$log2fc)
  expect_equal(out$mouse_id, "m1")
  expect_equal(out$group, "irrel")
  expect_false("population_order" %in% names(out))
})

test_that("facs_calc_diff runs on DB-shaped data with mouse_col/restim_col", {
  ref <- facs_calc_diff(log2fc_data, ref_pop = "CD3+")
  out <- facs_calc_diff(db_log2fc_data, ref_pop = "CD3+",
                        mouse_col = "mouse_id", restim_col = "condition")
  expect_equal(out$diff, ref$diff)
})
```

- [ ] **Step 2: Run to verify they fail**

Run: `Rscript --vanilla -e 'devtools::test(filter = "facs_calc")'`
Expected: FAIL, unused argument `mouse_col`.

- [ ] **Step 3: Implement**

Thread `mouse_col` into `calc_restim_proportions_()`; replace every hard-coded `"mouse_ID"` / bare `mouse_ID` there (group cols, passthrough exclusion list, `is_constant_()` grouping via `dplyr::across(dplyr::all_of(c(mouse_col, "tissue")))`, `id_cols`) and in the error/warning text (`{mouse_col}/tissue/{restim_col}`). Check `mouse_col %in% names(data)` with an error like the existing `restim_col` check. Existing tests must still pass (they match on `"m1"` only).

- [ ] **Step 4: Run to verify they pass**

Run: `Rscript --vanilla -e 'devtools::test(filter = "facs_calc")'`
Expected: all PASS, including the old log2fc/diff tests.

- [ ] **Step 5: Commit**

```bash
git add R/facs_calc.R tests/testthat/test-facs_calc.R man/
git commit -m "feat(facs): mouse_col lets facs_calc_log2fc()/diff() run on DB output"
```

---

### Task 6: `plot_box()` and `plot_population()` (+ ggpubr)

**Files:**
- Create: `R/plot_box.R`
- Modify: `DESCRIPTION` (add `ggpubr` between `ggplot2` and `glue`, `scales` between `rstatix` and `stats`)
- Test: `tests/testthat/test-plot_box.R`

**Interfaces:**
- Produces:
  - `plot_box(data, x, y, fill, scale = c("linear", "percent", "log"), colors = NULL) -> ggplot`
  - `plot_population(data, populations, metric, tissues = NULL, x = tissue, fill = group, facet = population, scale = c("linear", "percent", "log"), p_val = FALSE, colors = NULL, pseudocount = 1) -> ggplot`
  - `x`, `y`, `fill`, `facet` are bare column names (`{{ }}`); `metric` may be a vector (MFI panels).

Behaviour from spec decision 6: `plot_box` = `geom_boxplot(position_dodge(0.75), outlier.shape = NA)` + `geom_point(position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), alpha = 0.7)` + `scale_x_discrete(guide = guide_axis(n.dodge = 2))` + `theme_bw()`; percent -> `scale_y_continuous(labels = scales::percent)`; log -> `scale_y_log10(guide = "axis_logticks", minor_breaks = 0)`; `colors` non-NULL -> `scale_fill_manual(values = colors)`. `plot_population` filters, adds `pseudocount` to `value` for log, calls `plot_box(y = value)`, `facet_wrap(ggplot2::vars({{ facet }}), scales = "free_y")`, `p_val` -> `ggpubr::stat_compare_means(label = "p.format")`.

- [ ] **Step 1: Write the failing tests**

```r
pop_data <- tibble::tibble(
  mouse_id = rep(paste0("m", 1:4), each = 2),
  group = rep(c("irrel", "MPO409"), each = 4),
  tissue = rep(c("kidney", "lung"), 4),
  population = "CD4+", metric = "count_per_g",
  value = c(0, 10, 5, 20, 50, 60, 70, 80)
) |> dplyr::bind_rows(tibble::tibble(
  mouse_id = "m1", group = "irrel", tissue = "kidney",
  population = "CD8+", metric = "count_per_g", value = 3
))

test_that("plot_box returns a ggplot with the requested scales", {
  p <- plot_box(pop_data, x = tissue, y = value, fill = group, scale = "log",
                colors = c(irrel = "#F8766D", MPO409 = "#00BFC4"))
  expect_s3_class(p, "ggplot")
  expect_true(any(vapply(p$scales$scales, function(s) "fill" %in% s$aesthetics, logical(1))))
  expect_error(plot_box(pop_data, tissue, value, group, scale = "sqrt"))
})

test_that("plot_population filters and offsets values for log scale", {
  p <- plot_population(pop_data, populations = "CD4+", metric = "count_per_g",
                       tissues = "kidney", x = group, scale = "log")
  expect_equal(nrow(p$data), 4)
  expect_equal(sort(p$data$value), c(1, 6, 51, 71))
})

test_that("plot_population adds a p-value layer on request", {
  p <- plot_population(pop_data, "CD4+", "count_per_g", x = group, p_val = TRUE)
  expect_true(any(vapply(p$layers, function(l) inherits(l$stat, "StatCompareMeans"), logical(1))))
})
```

- [ ] **Step 2: Run to verify they fail**

Run: `Rscript --vanilla -e 'devtools::test(filter = "plot_box")'`
Expected: FAIL, functions not found.

- [ ] **Step 3: Implement both functions in `R/plot_box.R`, add `ggpubr` and `scales` to DESCRIPTION, `devtools::document()`.**

- [ ] **Step 4: Run to verify they pass**

Run: `Rscript --vanilla -e 'devtools::test(filter = "plot_box")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/plot_box.R tests/testthat/test-plot_box.R DESCRIPTION NAMESPACE man/
git commit -m "feat(plot): plot_box() and plot_population() from the analysis doc"
```

---

### Task 7: `plot_paired()`

**Files:**
- Modify: `R/plot_box.R`
- Test: `tests/testthat/test-plot_box.R`

**Interfaces:**
- Produces: `plot_paired(data, levels, condition = condition, fill = group, pair = mouse_id, facet = population, y = value, scale = c("percent", "linear", "log"), colors = NULL) -> ggplot`. `levels` (character, required, no default) = conditions to keep, in x order.

Behaviour from spec decision 6: keep rows with `condition %in% levels`, set `condition` factor levels to `levels`; `aes(x = interaction(condition, fill), y, fill)`; `geom_boxplot(outlier.shape = NA, alpha = 0.6)`; `geom_line(aes(group = pair), colour = "grey40", alpha = 0.4, linewidth = 0.3)`; `geom_point(shape = 21, colour = "black", alpha = 0.8, size = 2)`; `scale_x_discrete(guide = guide_axis(n.dodge = 2))`; `facet_wrap(vars(facet), scales = "free_y", nrow = 2)`; y scale as in `plot_box`; `colors` -> `scale_fill_manual`; `theme_bw()`. Error if a value of `levels` is absent from the data's condition column (message naming it).

- [ ] **Step 1: Write the failing tests**

```r
paired_data <- tibble::tibble(
  mouse_id = rep(c("m1", "m2"), each = 3),
  group = rep(c("irrel", "MPO409"), each = 3),
  condition = rep(c("unstim", "MPO", "PMA-Iono"), 2),
  population = "CD4+_CD44+_IFNg+", value = c(0.01, 0.03, 0.2, 0.02, 0.05, 0.3)
)

test_that("plot_paired keeps only the requested conditions, in order", {
  p <- plot_paired(paired_data, levels = c("unstim", "MPO"))
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$condition), c("unstim", "MPO"))
  expect_false("PMA-Iono" %in% p$data$condition)
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomLine"), logical(1))))
})

test_that("plot_paired errors on a level not in the data", {
  expect_error(plot_paired(paired_data, levels = c("unstim", "OVA")), "OVA")
})
```

- [ ] **Step 2: Run to verify they fail**

Run: `Rscript --vanilla -e 'devtools::test(filter = "plot_box")'`
Expected: new tests FAIL, function not found.

- [ ] **Step 3: Implement `plot_paired()` in `R/plot_box.R`; `devtools::document()`.**

- [ ] **Step 4: Run to verify they pass**

Run: `Rscript --vanilla -e 'devtools::test(filter = "plot_box")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/plot_box.R tests/testthat/test-plot_box.R NAMESPACE man/
git commit -m "feat(plot): plot_paired() for stim-vs-unstim per mouse"
```

---

### Task 8: Docs and full check

**Files:**
- Modify: `CLAUDE.md` (Data structures: `db_query_*()` bullet -- `population_order`, `db_query_subjects()` returning `group`, `facs_order_populations()`; `facs_calc_*()` bullet -- `mouse_col`; new bullet for `plot_box()` / `plot_population()` / `plot_paired()`; Imports list gets `ggpubr` and `scales`)
- Modify: `docs/superpowers/specs/2026-07-17-facs-db-schema-design.md` (`facs_measurements` gets `population_order INTEGER -- display position within the panel, NULL before 2026-11-19 exports`)

- [ ] **Step 1: Update both docs** (plain ASCII, matching the existing bullet style).

- [ ] **Step 2: Full test suite**

Run: `Rscript --vanilla -e 'devtools::test()'`
Expected: 0 failures.

- [ ] **Step 3: R CMD check**

Run: `Rscript --vanilla -e 'devtools::check()'`
Expected: 0 errors, 0 warnings; NOTEs only from the CLAUDE.md "Known check output" list (the `.git` hidden-file NOTE is expected in a worktree). Fix any new NOTE.

- [ ] **Step 4: Commit**

```bash
git add CLAUDE.md docs/superpowers/specs/2026-07-17-facs-db-schema-design.md
git commit -m "docs: population order, db_query_subjects, plot helpers"
```

---

## After merge (user, outside the repo)

Rerun `overview-prep_26-1.R`, `treg_prep_26-1.R`, `mpo-restim-facs_prep_26-1.R` with `final <- TRUE`, confirm each overwrite. Check: per-root row counts and values equal the previous imports; `SELECT COUNT(*) FROM facs_measurements WHERE population_order IS NULL` = 0.
