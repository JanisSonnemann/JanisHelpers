# Cohort tables, week labels and timecourse plots Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `meta_label_weeks()`, `meta_summary_table()`, `meta_mouse_table()`, `plot_timecourse()` and `plot_timecourse_summary()` to JanisHelpers and wire them into the analysis docs.

**Architecture:** Three new source files: `R/meta_weeks.R` (week labels), `R/meta_tables.R` (cohort normaliser + two `gt` tables), `R/plot_timecourse.R` (shared prep internals + two plots). All work on `db_load_experiment()` output plus the meta workbook's `meta` sheet; no DB connection.

**Tech Stack:** R, dplyr/tidyr, ggplot2, gt, rstatix, testthat 3.

**Spec:** `docs/superpowers/specs/2026-11-20-cohort-tables-weeks-timecourse-design.md`

## Global Constraints

- CLAUDE.md rules: `|>` only; every external call namespaced `pkg::fn()`; `.data$col` only for column names held in variables; `dplyr::select(!col)`; non-ASCII in R strings as `\uXXXX`, ASCII only in roxygen.
- Every exported function: `@param`, `@returns` (non-empty), `@export`; examples runnable without files or omitted.
- No new Imports (all needed packages -- `dplyr`, `tidyr`, `ggplot2`, `gt`, `rstatix`, `scales`, `stats`, `tibble`, `stringr` -- are already there). `ggrepel` stays Suggests.
- `devtools::document()` after roxygen changes; `devtools::check()` 0 errors / 0 warnings before the final commit.
- Test command (sandbox libraries):
  `L1=$HOME/Library/Caches/org.R-project.R/R/renv/library/facs-pipeline-workflow-0ff3068e/macos/R-4.6/aarch64-apple-darwin23; L2=$HOME/Library/Caches/org.R-project.R/R/renv/library/JanisHelpers-b26afe2b/macos/R-4.6/aarch64-apple-darwin23; R_LIBS=$L1:$L2:$TMPDIR/rlib Rscript --vanilla -e 'devtools::test(filter="<name>")'`
- Week rule: `week_num = round(days_post_bmt / 7)`; end-point samples (`collected_at == end_date`) of one experiment in the same ISO calendar week get `round(median(days) / 7)`.
- Excluded = mouse-wide (`tissue`, `panel`, `condition`, `file_name`, `collected_at` all NA) in every import unit -- reuse `fully_excluded_mice_(exclusions, units)` from `R/db_load.R`, with `units` = distinct `import_unit` / `wsp_root` values over the data tibbles of `dat`.
- Analysed = mouse with end-point data in >= 1 import unit and not excluded; early-ending mice with end-point data count.
- Summary plot: median + IQR (`stats::quantile(type = 7)`); groups with n < `min_n` drawn hollow, no error bar, bridged by the median line.

## Review Focus

1. Workbook from `meta_read()` has cleaned names (`id`, `dob`, `death_date`, `mouse_ID`) and POSIXct dates -- normaliser must accept these and DB names (`mouse_id`, `end_date`, Date). Test in Task 2.
2. Mouse ids sort naturally (`26-1-2` before `26-1-10`). Test in Task 3.
3. ELISA `OOR<` rows with `value = NA` -- dropped with a message, never plotted at 0. Test in Task 4.
4. Several replicates per mouse x week (restim ELISA has 2) -- averaged before medians, so n counts mice not wells. Test in Task 5.
5. `data` with mice not in `meta`, or mice lacking `bmt_date` -> clear error / NA-with-warning, not silent wrong weeks. Test in Task 1.

---

### Task 1: `meta_label_weeks()`

**Files:**
- Create: `R/meta_weeks.R`
- Test: `tests/testthat/test-meta_weeks.R`

**Interfaces:**
- Produces: `meta_label_weeks(data, meta, by = c("experiment_code", "mouse_id"), ref_col = "bmt_date", date_col = "collected_at")` -> `data` + `days_post_bmt` (integer), `week_num` (integer), `week` (ordered factor `"week <n>"`, levels numeric order, only occurring levels).

- [ ] **Step 1: Write the failing tests** (synthetic 26-1 pattern: `meta` with `experiment_code = "26-1"`, `mouse_id` m1..m3, `bmt_date = 2026-05-12`, `end_date` m1 `2026-06-10`, m2 `2026-06-15`, m3 `2026-06-16`; `data` bleeds on 05-05, 05-19, 06-02 for all, plus each mouse's end date)

```r
test_that("26-1 dates give weeks -1, 1, 3, 4, 5", {
  out <- meta_label_weeks(data, meta)
  expect_equal(out$week_num[out$collected_at == as.Date("2026-05-05")], rep(-1L, 3))
  expect_equal(unique(out$week_num[out$collected_at == as.Date("2026-06-10")]), 4L)
  expect_equal(unique(out$week_num[out$collected_at %in% as.Date(c("2026-06-15", "2026-06-16"))]), 5L)
  expect_equal(levels(out$week), c("week -1", "week 1", "week 3", "week 4", "week 5"))
  expect_true(is.ordered(out$week))
  expect_type(out$days_post_bmt, "integer")
})
test_that("end-points in one calendar week straddling a rounding boundary pool", {
  # bmt 2026-05-12 (Tue); ends d31 = 2026-06-12 (Fri) and d32 = 2026-06-13 (Sat)
  # round(31/7)=4, round(32/7)=5 -> pooled round(31.5/7) = 5 for both
  expect_equal(unique(out$week_num[is_end]), 5L)
})
test_that("serial samples are not pooled", { ... d31 and d32 as non-end-point bleeds -> 4 and 5 })
test_that("mice missing from meta error, naming them", {
  expect_error(meta_label_weeks(data_with_m9, meta), "m9")
})
test_that("missing bmt_date gives NA week with a warning naming the mouse", {
  expect_warning(out <- meta_label_weeks(data, meta_no_bmt_m1), "m1")
  expect_true(all(is.na(out$week_num[out$mouse_id == "m1"])))
})
test_that("existing output columns error", {
  expect_error(meta_label_weeks(dplyr::mutate(data, week = 1), meta), "week")
})
```

- [ ] **Step 2: Run** `devtools::test(filter="meta_weeks")` -- expect FAIL (`could not find function "meta_label_weeks"`).

- [ ] **Step 3: Implement** in `R/meta_weeks.R`. Coerce `collected_at`, `ref_col`, `end_date` with `as.Date()` (POSIXct-safe). ISO calendar week key: `format(date, "%G-%V")`. Pool only rows where `collected_at == end_date`, grouped by `experiment_code` (if present in `by`) x ISO week.

- [ ] **Step 4: Run tests** -- expect PASS.

- [ ] **Step 5: `devtools::document()`, commit** `feat(meta): meta_label_weeks() labels samples by week after BMT`.

---

### Task 2: cohort normaliser + `meta_summary_table()`

**Files:**
- Create: `R/meta_tables.R`
- Test: `tests/testthat/test-meta_tables.R`

**Interfaces:**
- Consumes: `fully_excluded_mice_()` (`R/db_load.R`).
- Produces (internal, used by Task 3): `cohort_normalise_(cohort)` -> tibble `mouse_id` (chr, NA allowed), `ear_tag` (chr), `group` (chr/factor, NA allowed), `sex`, `cage`, `strain`, `generation`, `dob`, `bmt_date`, `end_date` (Date), `early` (logical: `!is.na(early)` for workbook, `FALSE` if absent), `reason` (chr). Accepts workbook names (`mouse_ID`, `id`, `line`, `death_date`, `dob`) or DB names (`mouse_id`, `mouse_strain`, `end_date`).
- Produces (internal): `analysis_status_(cohort, dat)` -> cohort + `in_db`, `excluded` (lgl), `exclusion_reason`, `has_endpoint` (lgl), `analysed` (lgl).
- Produces: `meta_summary_table(workbook, dat = NULL, group_col = "group")` -> `gt_tbl`. `workbook = NULL` uses `dat$meta` with a message.

- [ ] **Step 1: Write failing tests** with a synthetic workbook (cleaned names, POSIXct dates) mirroring 26-1 in miniature: 4 DB mice (a: ok MPO409, b: ok irrel, c: early irrel after BMT with end-point titer, d: excluded mouse-wide MPO409) + 3 non-DB rows (`mouse_ID` NA: one pre-BMT death `bmt_date` NA, one removed at BMT with `group` NA and `death_date == bmt_date`, one died after BMT without data). `dat` = list(meta (a,b,c), exclusions (d, all scopes NA, import_unit "titer_x"), titer (a,b,c,d rows incl. end-points, `import_unit = "titer_x"`)).

```r
test_that("cohort normaliser accepts workbook and DB names", {
  expect_named(cohort_normalise_(wb), cohort_cols, ignore.order = TRUE)
  expect_equal(cohort_normalise_(dat$meta)$mouse_id, dat$meta$mouse_id)
  expect_s3_class(cohort_normalise_(wb)$end_date, "Date")
})
test_that("analysis status: excluded mouse-wide, early with end-point counts", {
  st <- analysis_status_(cohort_normalise_(wb), dat)
  expect_equal(st$analysed[match(c("a","b","c","d"), st$mouse_id)], c(TRUE, TRUE, TRUE, FALSE))
  expect_true(st$excluded[st$mouse_id %in% "d"])
})
test_that("partial exclusion does not exclude the mouse", { ... exclusion with tissue = "kidney" -> excluded FALSE })
test_that("summary table counts the whole cohort", {
  tb <- meta_summary_table(wb, dat)
  expect_s3_class(tb, "gt_tbl")
  d <- tb[["_data"]]
  expect_equal(row_value(d, "Enrolled"), c(MPO409 = "3", irrel = "3", `not assigned` = "1"))
  expect_equal(row_total(d, "Analysed"), 3)        # a, b, c
  expect_equal(row_total(d, "Dropped before / at BMT"), 2)
  expect_equal(row_total(d, "Ended early after BMT"), 2)   # c + non-DB post-BMT death
  expect_equal(row_total(d, "Excluded from analysis"), 1)
})
test_that("without dat the Analysed and Excluded rows are omitted", { ... })
test_that("workbook = NULL falls back to dat$meta with a message", {
  expect_message(meta_summary_table(NULL, dat), "not shown")
})
```
(`row_value` / `row_total` are small test helpers reading the `characteristic` row of `_data`; columns of `_data` are `characteristic` + one per group.)

- [ ] **Step 2: Run** -- expect FAIL.
- [ ] **Step 3: Implement** the two internals and `meta_summary_table()` in `R/meta_tables.R`. Build a character tibble (`characteristic`, one column per group level, `not assigned` only if any group NA), then `gt::gt()`; reasons of dropped / early mice as `gt::tab_footnote()` on their rows. Age at BMT in weeks, `median (min-max)`, 1 decimal. Row labels exactly as in the test.
- [ ] **Step 4: Run tests** -- PASS.
- [ ] **Step 5: document, commit** `feat(meta): meta_summary_table() cohort summary per group`.

---

### Task 3: `meta_mouse_table()`

**Files:**
- Modify: `R/meta_tables.R`
- Test: `tests/testthat/test-meta_tables.R`

**Interfaces:**
- Consumes: `cohort_normalise_()`, `analysis_status_()` (Task 2).
- Produces: `meta_mouse_table(dat, workbook = NULL)` -> `gt_tbl`; `_data` columns `mouse_id`, `ear_tag`, `group`, `sex`, `cage`, `age_at_bmt_wk`, `end_day`, one column per data unit of `dat` (its name), `status`.

- [ ] **Step 1: Write failing tests** (reuse Task 2 fixtures; add `dat$overview` with tissues kidney/spleen/lung for a and b, kidney only for c; mouse ids `26-1-2`, `26-1-10` in one extra fixture for sorting)

```r
test_that("mouse grid shows tissues, date counts and status", {
  d <- meta_mouse_table(dat, wb)[["_data"]]
  expect_equal(d$overview[d$mouse_id %in% "a"], "K S L")
  expect_equal(d$overview[d$mouse_id %in% "c"], "K")
  expect_equal(d$titer[d$mouse_id %in% "a"], "4/4")      # 3 bleeds + end-point
  expect_equal(d$titer[d$mouse_id %in% "c"], "4/4")      # early end-point still counts
  expect_match(d$status[d$mouse_id %in% "c"], "^early: .* \\(d [0-9]+\\)$")
  expect_match(d$status[d$mouse_id %in% "d"], "^excluded: ")
  expect_equal(sum(grepl("^not in DB", d$status)), 3)
  expect_equal(d$status[d$mouse_id %in% "a"], "ok")
})
test_that("partial exclusion shows as a note", {
  expect_match(status_of_partial, "^ok; excl\\. kidney")
})
test_that("mouse ids sort naturally within group", {
  expect_equal(d2$mouse_id, c("26-1-2", "26-1-10"))
})
test_that("blank cell for a unit without data for the mouse", { expect_equal(d$overview[d$mouse_id %in% "d"], "") })
```

- [ ] **Step 2: Run** -- FAIL.
- [ ] **Step 3: Implement.** Data units = elements of `dat` other than `meta`, `exclusions`, `workbook` that are data frames with `mouse_id`. Cell rule: unit has > 1 distinct tissue -> first letters of tissues present (upper case, sorted by the unit's tissue order of appearance... alphabetical is fine); else unit has serial samples (any `collected_at != end_date`) -> `present/expected`, expected = distinct serial dates + 1; else `"✓"`. Natural sort: order by group then `as.integer(stringr::str_extract(mouse_id, "[0-9]+$"))`, then ear tag for non-DB rows. Non-`ok` rows: `gt::tab_style(gt::cell_fill("#FFF3CD"), gt::cells_body(rows = status != "ok"))`.
- [ ] **Step 4: Run tests** -- PASS.
- [ ] **Step 5: document, commit** `feat(meta): meta_mouse_table() per-mouse assay accounting grid`.

---

### Task 4: `plot_timecourse()` and shared internals

**Files:**
- Create: `R/plot_timecourse.R`
- Modify: `R/plot_box.R` (generalise the fill-scale helper to also set colour)
- Test: `tests/testthat/test-plot_timecourse.R`

**Interfaces:**
- Consumes: `add_y_scale_()`, `add_fill_scale_()` (`R/plot_box.R`).
- Produces: `plot_timecourse(data, x = week_num, y = value, fill = group, subject = mouse_id, facet = NULL, scale = c("linear", "log"), colors = getOption("JanisHelpers.group_colors"), p_val = FALSE)` -> `ggplot`.
- Produces (internal, used by Task 5): `timecourse_prep_(data, x, y, fill, subject, facet)` (quosures) -> tibble `.x`, `.y`, `.fill`, `.subject`, `.facet` (or absent), one row per subject x x (x facet), `.y` = mean of replicates, NA `y` dropped with message `"Dropped <n> row(s) with missing values"`; `timecourse_axis_(p, x_values)`; `timecourse_lloq_(data, facet)` -> layer or `NULL`; `timecourse_stars_(prepped, y_pos)` -> tibble `.x`, `.facet`, `p_adj`, `label`.
- Modify: `add_fill_scale_(p, colors, fill_values, aesthetics = "fill")` -- passes `aesthetics` to `ggplot2::scale_fill_manual()` / use `ggplot2::scale_discrete_manual(aesthetics, values = colors)`; existing callers unchanged.

- [ ] **Step 1: Write failing tests** (fixture: 2 groups x 4 mice, weeks -1, 1, 3, 5; restim-like variant with 2 replicates and one `NA` value; `lloq` column variant)

```r
test_that("plot_timecourse returns a ggplot with mouse and median layers", {
  p <- plot_timecourse(tc, colors = c(MPO409 = "#00BFC4", irrel = "#F8766D"))
  expect_s3_class(p, "ggplot")
  expect_true(has_scale_(p, "colour"))
})
test_that("missing values are dropped with a message", {
  expect_message(plot_timecourse(tc_na), "Dropped 1 row")
})
test_that("replicates are averaged per mouse and week", {
  pr <- timecourse_prep_(tc_reps, ...)   # via rlang::quo
  expect_equal(nrow(pr), nrow(dplyr::distinct(tc_reps, mouse_id, week_num)))
})
test_that("lloq line only when lloq present", {
  expect_true(any(vapply(plot_timecourse(tc_lloq)$layers, \(l) inherits(l$geom, "GeomHline"), logical(1))))
  expect_false(any(vapply(plot_timecourse(tc)$layers, \(l) inherits(l$geom, "GeomHline"), logical(1))))
})
test_that("numeric x gets week labels at the weeks present", {
  b <- ggplot2::ggplot_build(plot_timecourse(tc))$layout$panel_params[[1]]$x$get_labels()
  expect_equal(b, c("week -1", "week 1", "week 3", "week 5"))
})
test_that("log scale with non-positive values errors", {
  expect_error(plot_timecourse(dplyr::mutate(tc, value = value - 100), scale = "log"), "positive")
})
test_that("p_val adds star labels", {
  expect_true(any(vapply(plot_timecourse(tc, p_val = TRUE)$layers, \(l) inherits(l$geom, "GeomText"), logical(1))))
})
test_that("facet splits by analyte", { expect_s3_class(plot_timecourse(tc2, facet = analyte)$facet, "FacetWrap") })
test_that("plot_box colours are unchanged", { ... existing has_scale_ fill test still passes })
```

- [ ] **Step 2: Run** -- FAIL.
- [ ] **Step 3: Implement.** Layers: `geom_line(aes(group = .subject, colour = .fill), alpha = 0.4, linewidth = 0.4)` + small points; group medians computed in the function (median of `.y` per `.fill` x `.x` x facet), drawn `linewidth = 1.2` + points size 3. Numeric `.x`: `scale_x_continuous(breaks = sort(unique(.x)), labels = paste("week", breaks))`; factor `.x`: discrete. Stars: per `.x` (x facet), 2 levels -> `rstatix::wilcox_test`, > 2 -> `rstatix::kruskal_test`; skip weeks with a group n < 2; `stats::p.adjust(method = "BH")` across weeks within facet; label `***`/`**`/`*`, `ns` dropped; y = 1.05 x max of the week (log: x 1.3). Log check: `any(.y <= 0)` -> `stop("scale = \"log\" needs positive values")`.
- [ ] **Step 4: Run tests** (`plot_timecourse` and `plot_box`) -- PASS.
- [ ] **Step 5: document, commit** `feat(plot): plot_timecourse() per-mouse timecourse with group medians`.

---

### Task 5: `plot_timecourse_summary()`

**Files:**
- Modify: `R/plot_timecourse.R`
- Test: `tests/testthat/test-plot_timecourse.R`

**Interfaces:**
- Consumes: Task 4 internals.
- Produces: `plot_timecourse_summary(data, x = week_num, y = value, fill = group, subject = mouse_id, facet = NULL, scale = c("linear", "log"), colors = getOption("JanisHelpers.group_colors"), p_val = FALSE, min_n = 2)` -> `ggplot`; internal `timecourse_summary_(prepped)` -> tibble `.fill`, `.x`, (`.facet`), `n`, `median`, `q25`, `q75`, `low_n` (`n < min_n`).

- [ ] **Step 1: Write failing tests** (fixture: irrel has one mouse at week 4, others n = 4; replicates variant)

```r
test_that("summary computes median and IQR over mice, not wells", {
  s <- timecourse_summary_(timecourse_prep_(tc_reps, ...), min_n = 2)
  expect_equal(s$n[s$.fill == "irrel" & s$.x == 1], 4L)
  expect_equal(s$q25[1], unname(stats::quantile(expected_vals, 0.25, type = 7)))
})
test_that("low-n points are hollow, without error bar, and bridged", {
  p <- plot_timecourse_summary(tc_low)
  ld <- ggplot2::layer_data(p, i = which_errorbar(p))
  expect_false(any(ld$x == 4 & ld$colour == irrel_colour))  # no bar at week 4 irrel
  line <- ggplot2::layer_data(p, i = which_line(p))
  expect_false(4 %in% line$x[line$group == irrel_group])     # bridged: 3 -> 5
  expect_true(all(c(3, 5) %in% line$x[line$group == irrel_group]))
})
test_that("n labels are printed per group and week", {
  expect_true(any(vapply(p$layers, \(l) inherits(l$geom, "GeomText"), logical(1))))
})
test_that("default caption names the measure", {
  expect_equal(p$labels$caption, "median and IQR; n per group below")
})
```
(`which_errorbar` / `which_line` find the layer index by geom class.)

- [ ] **Step 2: Run** -- FAIL.
- [ ] **Step 3: Implement.** `position_dodge(width = 0.4)` for points, error bars (`width = 0.3`), and lines (line drawn from rows with `!low_n` only, so it bridges). Low-n points `shape = 21, fill = "white"`. n labels: `geom_text(aes(y = -Inf, label = n), vjust = -0.5, position = dodge, size = 3)`. Stars reuse `timecourse_stars_()` on prepped data.
- [ ] **Step 4: Run tests** -- PASS.
- [ ] **Step 5: document, commit** `feat(plot): plot_timecourse_summary() median and IQR timecourse`.

---

### Task 6: Template, CLAUDE.md, check

**Files:**
- Modify: `inst/rmarkdown/templates/analysis-db/skeleton/skeleton.Rmd`
- Modify: `CLAUDE.md` (Data structures map)
- Test: `tests/testthat/test-report_templates.R` (existing parse test must still pass)

- [ ] **Step 1: Template.** Setup: commented `# dat$workbook <- meta_read(here::here("data", "<exp>", "meta", "meta_<exp>.xlsx"))$meta` inside `local()` before returning `dat`. Overview section: `meta_summary_table(dat$workbook, dat)` and `meta_mouse_table(dat, dat$workbook)` (both work with `dat$workbook` NULL). ELISA section: commented example `dat$<elisa> <- meta_label_weeks(dat$<elisa>, dat$meta)`, `plot_timecourse(dat$<elisa>)`, `plot_timecourse_summary(dat$<elisa>)`.
- [ ] **Step 2: CLAUDE.md** -- add two bullets to "Data structures": meta tables + `meta_label_weeks()` (workbook as cohort source, mouse-wide exclusion rule, analysed definition, week rule); timecourse plots (`week_num` numeric x, replicate mean, median+IQR, `min_n` bridge). Update the `analysis-db` template row.
- [ ] **Step 3: Run** `devtools::test(filter="report_templates")` -- PASS; then `devtools::check()` -- 0 errors, 0 warnings, only the known NOTEs.
- [ ] **Step 4: Commit** `feat(report): analysis-db template uses cohort tables and timecourse plots`.

---

### Task 7: Wire into `analysis_26-1_db.Rmd` (OneDrive, outside repo)

**Files:**
- Modify: `~/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/Shared Documents/MPO_mRNA/data/26-1_BMT-mRNA-V4/analysis/analysis_26-1_db.Rmd` (needs a sandbox-override write; not in git)

- [ ] **Step 1:** Setup: replace the `dat$dropouts` block with `dat$workbook <- meta_read(...)$meta`; add `dat$titer <- meta_label_weeks(dat$titer, dat$meta)`.
- [ ] **Step 2:** New `# Cohort` section before `# Changes vs. the wsp version`: `meta_summary_table(dat$workbook, dat)`. QC: replace chunks `qc-overview`, `qc-dropouts`, `qc-missing-samples`, `qc-exclusions` with one `qc-mice` chunk `meta_mouse_table(dat, dat$workbook)`.
- [ ] **Step 3:** Titer section: replace `mpo-elisa` chunk with `plot_timecourse(dat$titer) + labs(title = "anti-MPO titer", y = "AU")` and `plot_timecourse_summary(dat$titer) + labs(...)`.
- [ ] **Step 4:** Install the branch into a scratch lib and knit the doc to a temp output dir; confirm it renders, Enrolled = 23 (MPO409 / irrel / not assigned), Analysed = 16, titer weeks -1/1/3/4/5.
