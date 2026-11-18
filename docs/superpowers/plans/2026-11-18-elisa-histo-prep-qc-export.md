# ELISA + histology prep / QC / export Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Curated ELISA and histology sheets go through a FACS-style prep -> QC report (with outlier plots) -> `final` export loop into the DuckDB lab database, whose samples now carry a collection date.

**Architecture:** Samples become `mouse x tissue x collected_at` (empty date = the subject's `end_date`). ELISA and histology each get their own reader, prep save/load/run, QC template, knit wrapper, export and script template (option 2B). Identical low-level pieces (exclusion matching, replace-import transaction, whisker-outlier rule, QC render) are shared internal helpers.

**Tech Stack:** R package (roxygen2, testthat 3e), DuckDB via DBI/dbplyr, readxl, ggplot2, rmarkdown.

**Spec:** `docs/superpowers/specs/2026-11-18-elisa-histo-prep-qc-export-design.md`

## Global Constraints

- Follow CLAUDE.md: `|>` pipe; every external call `pkg::fn()`; `dplyr::select(!col)`; non-ASCII in R strings as `\uXXXX`, plain ASCII in roxygen.
- Function names: `elisa_`, `histo_`, `db_`, `report_` domains only. `histo_` must be declared in CLAUDE.md before the first `histo_` function (Task 7).
- Every exported function: `@param`, `@returns`, `@export`; `@examples` in `\dontrun{}` when it needs a file or a DB.
- No new `Imports`. `ggrepel` stays optional behind `requireNamespace()`.
- `db_` tests use `local_test_db()` from `tests/testthat/helper-db.R`.
- Run `devtools::document()` after roxygen changes; `devtools::test()` per task; `devtools::check()` at the end: 0 errors, 0 warnings, only the NOTEs listed in CLAUDE.md.
- Curated sheet `values`; ELISA optional sheet `range` with row labels `LLOQ`, `ULOQ`, optional `unit` in column 1.
- `result_status` values exactly `"ok"`, `"OOR<"`, `"OOR>"`. Replicate CV threshold default `max_cv = 20` (percent).
- Import unit (root) = `tools::file_path_sans_ext(basename(path))`. Prep at `<dirname(path)>/processed/<root>_prep.rds`, report at `<dirname(path)>/qc/<root>_qc.html`.
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

1. The exclusions sheet names an ELISA or histology assay in `panel` -> `db_export_facs()` must not reject it as an unknown panel (its check must accept assays of every domain). Test in Task 3.
2. A subject without `end_date` (no `death_date` in meta) gets an end-point sample -> error naming every such mouse, nothing written. Test in Task 1.
3. A sheet has a row with an empty date and a row whose date equals that mouse's `end_date` -> after resolution they are one sample: the export errors listing the clash and writes nothing. Test in Task 5.
4. An analyte cell holds text other than a number, `OOR<` or `OOR>` (e.g. `n.d.`) -> the reader errors naming analyte, mouse and value, rather than storing `NA` silently. Test in Task 4.
5. The `date` column holds Excel dates, ISO text or a mix -> all become `Date`; an unparseable entry errors naming the mouse and value. Test in Task 4 (shared parser used by Task 7).

---

### Task 1: Samples carry a collection date

**Files:**
- Modify: `R/db_connect.R` (samples table, old-schema detection)
- Modify: `R/db_write_dimensions.R:102-130` (`db_write_samples()`)
- Modify: `R/db_write_measurements.R:175-223` (`lookup_sample_id_()`, `lookup_sample_ids_()`)
- Modify: `R/db_export.R:154-161` (FACS sample resolution), `R/db_query.R` (`db_query_facs()` adds `collected_at`)
- Modify: `tests/testthat/helper-db.R`, `tests/testthat/test-db_write_dimensions.R`, `test-db_write_measurements.R`, `test-db_connect.R`, `test-db_query.R`

**Interfaces:**
- Produces: `samples.collected_at DATE NOT NULL`, `UNIQUE(subject_id, tissue, collected_at)`.
- Produces: `db_write_samples(con, data)` -- `data` has `mouse_id`, `tissue`, optional `collected_at`; `NA`/absent `collected_at` = subject's `end_date`. Returns invisibly the rows inserted.
- Produces: `lookup_sample_id_(con, mouse_id, tissue, collected_at = NA)` and `lookup_sample_ids_(con, mouse_id, tissue, collected_at = NA)` (vectorised; `NA` = `end_date`; error lists every unresolved `mouse x tissue x date`).
- Produces: internal `db_old_schema_(con)` -> character vector of reasons (empty = current).

- [ ] **Step 1: Update fixtures and write failing tests**

`helper-db.R`: `seed_minimal_fixture()` registers `25-7-1` with `end_date = as.Date("2026-06-10")`. Give every subject other tests create an `end_date` too.

In `test-db_write_dimensions.R`:
```r
test_that("db_write_samples defaults collected_at to the subject's end_date", {
  con <- local_test_db(); seed_minimal_fixture(con)
  got <- DBI::dbGetQuery(con, "SELECT collected_at FROM samples")$collected_at
  expect_equal(got, as.Date("2026-06-10"))
})
test_that("db_write_samples keeps one sample per mouse x tissue x date", {
  con <- local_test_db(); seed_minimal_fixture(con)
  n <- db_write_samples(con, tibble::tibble(
    mouse_id = "25-7-1", tissue = "serum",
    collected_at = as.Date(c("2026-05-06", "2026-05-20"))))
  expect_equal(n, 2)
  expect_equal(db_write_samples(con, tibble::tibble(mouse_id = "25-7-1", tissue = "spleen")), 0)
})
test_that("db_write_samples names every subject without end_date and writes nothing", {
  con <- local_test_db()
  db_write_experiment(con, "25-7")
  db_write_subjects(con, tibble::tibble(mouse_id = c("25-7-8", "25-7-9"), experiment_code = "25-7"))
  err <- expect_error(db_write_samples(con, tibble::tibble(mouse_id = c("25-7-8", "25-7-9"), tissue = "kidney")))
  expect_match(conditionMessage(err), "25-7-8.*25-7-9")
  expect_match(conditionMessage(err), "end_date")
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM samples")$n, 0)
})
```
In `test-db_write_measurements.R`: `lookup_sample_ids_()` with `collected_at = c(NA, as.Date("2026-05-06"))` resolves both (end + dated serum), and an unregistered date errors with `"25-7-1 x serum x 2026-05-07"` in the message.
In `test-db_connect.R`: a DuckDB file whose `samples` was created with nullable `collected_at` and `UNIQUE(subject_id, tissue)` makes `db_connect()` error with `"older schema"` and `"Recreate"`.
In `test-db_query.R`: `db_query_facs()` has a `collected_at` column.

- [ ] **Step 2: Run tests, verify they fail**

Run: `Rscript -e 'devtools::test(filter = "db_")'`
Expected: FAIL on the new expectations.

- [ ] **Step 3: Implement**

- `samples`: `collected_at DATE NOT NULL`, `UNIQUE(subject_id, tissue, collected_at)`.
- `db_old_schema_()`: reasons from `information_schema.columns` -- `facs_stains` missing `wsp_root`/`condition` (existing check, moved here), `samples.collected_at` nullable, `elisa_measurements` without `analyte` (used by Task 2). `db_connect()` disconnects and errors `"'<path>' was created by an older schema (<reasons>). Recreate the database file and re-import."` when non-empty.
- `db_write_samples()`: fill `collected_at` via `COALESCE(t.collected_at, sub.end_date)`; before inserting, error listing mice where both are NULL; `ON CONFLICT (subject_id, tissue, collected_at) DO NOTHING`.
- Lookups join on `sm.collected_at = COALESCE(k.collected_at, sub.end_date)`.
- `db_export_facs()`: unchanged calls (no date = end date) -- verify `test-db_export.R` still passes.

- [ ] **Step 4: Run tests, verify they pass**

Run: `Rscript -e 'devtools::test(filter = "db_")'`
Expected: PASS (all `db_` files, incl. `test-db_export.R`).

- [ ] **Step 5: Commit**

```bash
git add R/db_connect.R R/db_write_dimensions.R R/db_write_measurements.R R/db_export.R R/db_query.R tests/testthat man
git commit -m "feat(db): samples keyed by mouse x tissue x collection date"
```

---

### Task 2: New ELISA / histology tables, queries; remove append-only writers

**Files:**
- Modify: `R/db_connect.R` (`elisa_measurements`, `histo_measurements`, `import_exclusions`)
- Modify: `R/db_write_measurements.R` (delete `db_write_elisa()`, `db_write_histo()` and their roxygen)
- Modify: `R/db_query.R` (`db_query_elisa()`, `db_query_histo()`)
- Modify: `R/db_import_log.R:19` (comment: callers are the three `db_export_*()`)
- Modify: `tests/testthat/test-db_write_measurements.R` (drop `db_write_elisa`/`histo` tests), `test-db_query.R`, `test-db_connect.R`

**Interfaces:**
- Produces tables:
  - `elisa_measurements(elisa_measurement_id, sample_id NOT NULL REFERENCES samples, assay_id NOT NULL REFERENCES assays, import_unit VARCHAR NOT NULL, analyte VARCHAR NOT NULL, replicate INTEGER NOT NULL, value DOUBLE, unit VARCHAR, result_status VARCHAR NOT NULL, lloq DOUBLE, uloq DOUBLE, source_file VARCHAR, imported_at TIMESTAMP DEFAULT current_timestamp, UNIQUE(sample_id, assay_id, analyte, replicate))`
  - `histo_measurements(histo_measurement_id, sample_id NOT NULL, assay_id NOT NULL, import_unit VARCHAR NOT NULL, metric VARCHAR NOT NULL, value DOUBLE, source_file VARCHAR, imported_at ..., UNIQUE(sample_id, assay_id, metric))`
  - `import_exclusions` gains `collected_at DATE`.
- Produces: `db_query_elisa(con)` columns `experiment_code, mouse_id, tissue, collected_at, assay_name, import_unit, analyte, replicate, value, unit, result_status, lloq, uloq, source_file, imported_at`; `db_query_histo(con)` columns `experiment_code, mouse_id, tissue, collected_at, assay_name, import_unit, metric, value, source_file, imported_at`.

- [ ] **Step 1: Write failing tests**

`test-db_query.R`: insert one ELISA and one histo row with plain SQL `INSERT` (sample from `seed_minimal_fixture()`), then
```r
expect_named(dplyr::collect(db_query_elisa(con)), c("experiment_code", "mouse_id", "tissue", "collected_at", "assay_name", "import_unit", "analyte", "replicate", "value", "unit", "result_status", "lloq", "uloq", "source_file", "imported_at"))
expect_named(dplyr::collect(db_query_histo(con)), c("experiment_code", "mouse_id", "tissue", "collected_at", "assay_name", "import_unit", "metric", "value", "source_file", "imported_at"))
```
`test-db_connect.R`: `"import_exclusions"` fields include `"collected_at"`; `exists("db_write_elisa", asNamespace("JanisHelpers"))` is `FALSE`.

- [ ] **Step 2: Run, verify failure** -- `Rscript -e 'devtools::test(filter = "db_")'` -> FAIL.

- [ ] **Step 3: Implement** the schema, queries and removals above; `devtools::document()` (removes `man/db_write_elisa.Rd`, `man/db_write_histo.Rd`, NAMESPACE entries).

- [ ] **Step 4: Run, verify pass** -- same command -> PASS.

- [ ] **Step 5: Commit**

```bash
git add R tests/testthat man NAMESPACE
git commit -m "feat(db): replaceable ELISA/histo tables with dates; drop append-only writers"
```

---

### Task 3: Shared internals -- exclusions, replace-import, outliers, render

**Files:**
- Create: `R/db_export_table.R` (replace-import + overwrite prompt + panel check)
- Create: `R/db_exclude.R` (internal sample-level exclusion matching)
- Modify: `R/db_export.R` (use `check_exclusion_panels_()`; write `collected_at = NA` to `import_exclusions`)
- Modify: `R/facs_qc.R:80-138` (extract the whisker rule)
- Modify: `R/report_knit.R:122-145` (extract render helper)
- Test: `tests/testthat/test-db_export_table.R`, `test-db_exclude.R`, `test-facs_qc.R`, `test-db_export.R`

**Interfaces:**
- Produces: `match_sample_exclusions_(data, exclusions, assay)` -- `data` has `mouse_ID`, `tissue`, `collected_at` (Date, `NA` = end). Returns `list(applicable, excluded_samples, unmatched)`: `applicable` columns `mouse_ID, tissue, panel, collected_at (Date), reason, date (Date), n_samples_matched (int)`; `excluded_samples` = distinct `mouse_ID, tissue, collected_at, reason`. Rows apply when `panel` is empty or equals `assay`; non-empty `tissue`/`collected_at` must equal the sample's (an exclusion `collected_at` never matches an `NA` sample date). Errors: missing `mouse_ID`/`reason` column or value (same wording as `facs_match_exclusions()`); an applicable row with non-empty `condition` or `file_name` -> `"condition / file_name don't apply to <assay>"`.
- Produces: `check_exclusion_panels_(con, exclusions, assay)` -- errors naming panels that are neither `assay` nor any registered assay of any domain. Replaces the FACS-only check in `db_export_facs()`.
- Produces: `confirm_replace_(con, table, root, overwrite)` -> `TRUE` to proceed / `FALSE` declined (asks via `utils::askYesNo()` when the table has rows with `import_unit = root`; same non-interactive error as FACS).
- Produces: `replace_import_(con, domain, table, rows, root, source_file, qc_report, applicable)` -> `import_id`. One `DBI::dbWithTransaction()`: `DELETE FROM <table> WHERE import_unit = root`; insert `rows` (columns = table columns); `write_import_log_()` with `qc_passed = TRUE`, `source_version = NA`, `supersedes_import_id` = latest log id of `domain` + `root`; insert `applicable` into `import_exclusions` (`file_name`, `condition` NULL, `n_files_matched` = `n_samples_matched`).
- Produces: `qc_whisker_side_(t)` -> character vector `"low"` / `"high"` / `NA` (1.5 x IQR rule on `t`, quantile type 7, `na.rm`). `qc_score_()` uses it per group.
- Produces: `render_qc_(template, prep_path, exclusions)` -> html path in `<dirname(dirname(prep_path))>/qc/`, `sub("_prep\\.rds$", "_qc.html", ...)`. `report_knit_facs_qc()` calls it with `"facs-qc-template"`.

- [ ] **Step 1: Write failing tests**

- `test-db_exclude.R`: a row with only `mouse_ID` excludes all that mouse's samples; a row with `collected_at = "2026-05-06"` excludes only that date; a row for panel `"overview"` is ignored for assay `"titer"`; a row matching nothing lands in `unmatched`; a row with `file_name` set errors with `"don't apply"`.
- `test-db_export.R`: exclusions with `panel = "titer"` where `db_write_assay(con, "titer", "elisa")` is registered -> `db_export_facs()` succeeds (Review Focus 1); `panel = "Treg"` unregistered still errors `"unknown panel"`.
- `test-db_export_table.R` (using a registered sample and a temp `histo_measurements` row set): two `replace_import_()` calls for root `"histo_25-7"` with different values leave only the second call's rows; the second log row's `supersedes_import_id` equals the first id; a `rows` tibble with a duplicate key errors and leaves both table and `import_log` unchanged; `confirm_replace_(con, ..., overwrite = FALSE)` returns `FALSE` when the root exists.
- `test-facs_qc.R`: `qc_whisker_side_(c(1, 2, 2, 3, 2, 50))` gives `NA` x5 then `"high"`; existing `qc_score_` tests still pass.

- [ ] **Step 2: Run, verify failure** -- `Rscript -e 'devtools::test(filter = "db_exclude|db_export|facs_qc")'` -> FAIL.

- [ ] **Step 3: Implement** the six helpers with the signatures above. In `db_export_facs()`, swap the inline panel check for `check_exclusion_panels_()` and add `collected_at` (NA) to the `import_exclusions` insert.

- [ ] **Step 4: Run, verify pass** -- same command plus `filter = "report_knit"` -> PASS.

- [ ] **Step 5: Commit**

```bash
git add R tests/testthat
git commit -m "refactor(db): shared exclusion, replace-import and QC render helpers"
```

---

### Task 4: `elisa_read_curated()`

**Files:**
- Modify: `R/elisa_read.R` (add function + internal `parse_curated_dates_()`)
- Create: `tests/fixtures/make_curated_fixtures.R` (generator, uses `writexl` -- install locally, not a package dependency), `tests/fixtures/elisa_curated.xlsx`, `tests/fixtures/histo_curated.xlsx`; note both in `tests/fixtures/README.md`
- Test: `tests/testthat/test-elisa_read.R`

**Interfaces:**
- Produces: `elisa_read_curated(path, tissue = NA_character_)` -> tibble, one row per sample x analyte x replicate: `mouse_ID` (chr), `tissue` (chr), `collected_at` (Date, `NA` = end), `analyte` (chr, column header verbatim), `replicate` (int, row order within `mouse_ID x tissue x collected_at`), `value` (dbl, `NA` for OOR), `result_status` (`"ok"`/`"OOR<"`/`"OOR>"`), `unit` (chr), `lloq`, `uloq` (dbl). Returned invisibly with a summary message (like `elisa_read_results()`).
- Produces: internal `parse_curated_dates_(x, mouse_ID)` -> Date; accepts POSIXct, Date, Excel serial numbers, ISO text; empty -> `NA`; errors naming mouse + value otherwise. Reused by Task 7.

Fixture `elisa_curated.xlsx`: sheet `values` columns `mouse_ID, date, TNFa, IL-2`; rows: `M1` x2 (empty date, values 4.2/3.9 and `OOR<`/10.1), `M2` x2 (empty date), `M1` with date `2026-05-06` x1; sheet `range` rows `LLOQ 7.81 15.6`, `ULOQ 31.3 250`, `unit pg/ml pg/ml`.

- [ ] **Step 1: Write failing tests**

```r
test_that("elisa_read_curated returns long rows with replicates and range", {
  d <- elisa_read_curated(test_path("../fixtures/elisa_curated.xlsx"), tissue = "spleen")
  expect_named(d, c("mouse_ID", "tissue", "collected_at", "analyte", "replicate", "value", "result_status", "unit", "lloq", "uloq"))
  m1 <- dplyr::filter(d, mouse_ID == "M1", is.na(collected_at), analyte == "IL-2")
  expect_equal(m1$replicate, 1:2)
  expect_equal(m1$result_status, c("OOR<", "ok"))
  expect_true(is.na(m1$value[1]))
  expect_equal(unique(d$lloq[d$analyte == "TNFa"]), 7.81)
  expect_equal(unique(d$unit), "pg/ml")
  expect_equal(unique(d$tissue), "spleen")
  expect_equal(sum(!is.na(d$collected_at)), 2)
})
```
Plus: a `tissue` column in the sheet overrides the argument only where non-empty; no `tissue` column and `tissue = NA` errors `"tissue"`; no `range` sheet -> `lloq`/`uloq`/`unit` `NA`; analyte cell `"n.d."` errors with `"n.d."`, the analyte and the mouse (Review Focus 4; build the bad sheet in a temp file with `writexl` behind `skip_if_not_installed("writexl")`); `parse_curated_dates_()` on `list(as.POSIXct("2026-05-06", tz = "UTC"), "2026-05-20", 46147, NA)` gives the three dates + `NA`, and on `"15/16.6.2026"` errors naming the mouse (Review Focus 5).

- [ ] **Step 2: Run, verify failure** -- `Rscript -e 'devtools::test(filter = "elisa_read")'` -> FAIL.

- [ ] **Step 3: Implement** with `readxl::read_excel(..., col_types = "list")` for `values` so mixed number/`OOR` cells survive; empty `mouse_ID` rows dropped.

- [ ] **Step 4: Run, verify pass** -- same command -> PASS.

- [ ] **Step 5: Commit**

```bash
git add R/elisa_read.R tests/fixtures tests/testthat/test-elisa_read.R man NAMESPACE
git commit -m "feat(elisa): read curated ELISA sheets"
```

---

### Task 5: ELISA prep save / load and `db_export_elisa()`

**Files:**
- Create: `R/elisa_prep.R`, `R/db_export_elisa.R`
- Test: `tests/testthat/test-elisa_prep.R`, `tests/testthat/test-db_export_elisa.R`

**Interfaces:**
- Consumes: Task 3 helpers; Task 4 data shape.
- Produces: `elisa_prep_save(data, path, assay, group_col = NULL, batch_col = NULL, expected_mice = NULL, max_cv = 20, log_scale = FALSE)` -> invisibly the `.rds` path. Saves `list(data, info)`; `info` = `source_path`, `root`, `assay`, `created_at`, `package_version`, `qc` (list of the QC arguments; `expected_mice` default = mice in `data`). Errors on missing Task 4 columns.
- Produces: `elisa_prep_load(path)` -> prep plus `paths` (`prep`, `qc`); refuses (message says rerun with `final <- FALSE`) if either file is missing, the report is older than the prep, or the `.xlsx` is newer than the prep.
- Produces: `db_export_elisa(con, prep, exclusions = NULL, data_dir = "data", overwrite = NULL)` -> invisibly `import_id` or `NULL` if declined. Order: require `prep$paths`; `check_exclusion_panels_()`; `match_sample_exclusions_()` (unmatched -> error, nothing written); drop excluded samples; error if nothing left; resolve `NA` dates to `end_date` (error listing mice without one); error if any `mouse_ID x tissue x collected_at x analyte x replicate` occurs twice after resolution; relative `source_file`/`qc_report` via `relative_to_data_dir_()`; `confirm_replace_()`; `db_write_samples()`; `db_write_assay(con, assay, "elisa")`; build rows with `lookup_sample_ids_()`; `replace_import_(con, "elisa", "elisa_measurements", ...)`; summary message.

- [ ] **Step 1: Write failing tests**

- `test-elisa_prep.R`: save writes `<dir>/processed/elisa_curated_prep.rds` (copy fixture to a temp dir); load errors `"final <- FALSE"` with no report; load errors when the `.xlsx` is touched after saving (`Sys.setFileTime()`).
- `test-db_export_elisa.R` (subjects with `end_date`, a fake `qc/` html written beside the prep, `data_dir` = temp root): rows written = non-excluded rows; `db_query_elisa()` shows `collected_at` = `end_date` for empty dates and `2026-05-06` for the dated row; a second export after changing one value with `overwrite = TRUE` leaves the new value only and the log row supersedes the first; `overwrite = FALSE` returns `NULL` and changes nothing; exclusion of `M2` drops M2 and writes one `import_exclusions` row with `n_files_matched = 2`; an unmatched exclusion errors and writes nothing; a sheet row dated the mouse's `end_date` plus an empty-date row for the same mouse/tissue errors with `"twice"` and writes nothing (Review Focus 3).

- [ ] **Step 2: Run, verify failure** -- `Rscript -e 'devtools::test(filter = "elisa_prep|db_export_elisa")'` -> FAIL.

- [ ] **Step 3: Implement** as specified.

- [ ] **Step 4: Run, verify pass** -- same -> PASS.

- [ ] **Step 5: Commit**

```bash
git add R tests/testthat man NAMESPACE
git commit -m "feat(elisa): prep save/load and replaceable database export"
```

---

### Task 6: ELISA QC report, `elisa_prep_run()`, script template

**Files:**
- Create: `R/elisa_qc.R` (internal checks + plots), `inst/rmarkdown/templates/elisa-qc-template/{template.yaml,skeleton/skeleton.Rmd}`, `inst/scripts/elisa_prep_template.R`
- Modify: `R/report_knit.R` (`report_knit_elisa_qc()`), `R/elisa_prep.R` (`elisa_prep_run()`)
- Create: `R/qc_table.R` (internal, shared ELISA/histo: completeness, empty grid, outlier scoring + plots)
- Test: `tests/testthat/test-elisa_qc.R`, `test-report_knit.R`, `test-elisa_prep.R`

**Interfaces:**
- Produces (internal, `R/qc_table.R`): `qc_timepoint_(collected_at)` -> chr (`format(date)` or `"end"`); `qc_table_completeness_(data, expected_mice)` -> tibble `tissue, timepoint, missing` (comma-joined mouse IDs) plus attribute `unknown` (mice not in `expected_mice`); `qc_table_empty_(data, var)` -> counts of `NA` `value` per `var x tissue x timepoint`; `qc_table_score_(data, var, group_col, log_scale)` -> `data` plus `side`, `flag` using `qc_whisker_side_()` within `group_col x tissue x timepoint x var` (no group -> no flags); `qc_table_plot_(scored, var, group_col, log_scale, batch_col = NULL)` -> ggplot, boxplot + jitter, facets `var ~ timepoint`, flagged points labelled with `mouse_ID` (`ggrepel` if available), points coloured by `batch_col` when given (the report's batch view); `qc_table_time_plot_(data, var, group_col)` -> per-mouse lines over `collected_at`, `NULL` if fewer than 2 dates.
- Produces (internal, `R/elisa_qc.R`): `elisa_qc_range_(data)` -> per analyte `n`, `n_oor_low`, `n_oor_high`, `n_below_lloq`, `n_above_uloq`; `elisa_qc_replicates_(data, max_cv)` -> samples x analyte with >= 2 numeric replicates and CV (sd / mean x 100) > `max_cv`, with `cv` and `values` (comma-joined); `elisa_qc_means_(data)` -> replicate means per sample x analyte (input to scoring/plots).
- Produces: `report_knit_elisa_qc(prep_path, exclusions = NULL)` -> html path (via `render_qc_("elisa-qc-template", ...)`).
- Produces: `elisa_prep_run(final, path, assay, tissue = NA_character_, prepare = identity, mice = NULL, exclusions = NULL, con = NULL, group_col = NULL, batch_col = NULL, ..., data_dir = "data", snapshot_root = NULL, overwrite = NULL)` -- mirrors `facs_prep_run()`: `final = FALSE` reads (`elisa_read_curated()`), `prepare()`, joins `group_col`/`batch_col` from `mice` with `meta_annotate(by = "mouse_ID")`, `elisa_prep_save(..., ...)`, renders, returns report path; `final = TRUE` requires `con`, `elisa_prep_load()`, `db_export_elisa()`, snapshot if `snapshot_root` and an id came back.
- Report sections, in order: Overview (samples per tissue x timepoint; file + prep times), Completeness, Exclusions (flagged samples; unmatched rows as warnings), Empty values, Range (OOR / LLOQ / ULOQ table), Replicate agreement, Outliers (`qc_table_plot_()` on replicate means; time plot when present), Flagged samples (ranked by analyte count).

- [ ] **Step 1: Write failing tests**

- `test-elisa_qc.R` on the fixture data: `elisa_qc_range_()` gives TNFa `n_below_lloq` = number of TNFa values < 7.81 and IL-2 `n_oor_low = 1`; `elisa_qc_replicates_()` with `max_cv = 5` returns M1/TNFa (4.2 vs 3.9 -> CV 5.2) and with `max_cv = 20` returns none; `qc_table_completeness_()` with `expected_mice = c("M1", "M2", "M3")` lists `M3` missing for `end` and `M2, M3` for `2026-05-06`; `qc_table_score_()` flags the one extreme point in a synthetic group of `c(1, 2, 2, 3, 2, 50)`.
- `test-report_knit.R`: `report_knit_elisa_qc()` on a saved fixture prep with `group_col` and one exclusion renders `qc/elisa_curated_qc.html` whose text contains `"Replicate agreement"`, `"Outliers"` and the excluded mouse ID.
- `test-elisa_prep.R`: `elisa_prep_run(FALSE, ...)` returns an existing html path; `elisa_prep_run(TRUE, ..., con = NULL)` errors `"con"`.

- [ ] **Step 2: Run, verify failure** -- `Rscript -e 'devtools::test(filter = "elisa|report_knit")'` -> FAIL.

- [ ] **Step 3: Implement** internals, template (same params `prep`, `exclusions` and figure settings as `facs-qc-template`), wrapper, `elisa_prep_run()`, and `inst/scripts/elisa_prep_template.R` shaped like `facs_prep_template.R`: `final`, experiment block (meta, mice, db_path), assay block (`path`, `assay`, `tissue`), `prepare`, run block with `db_write_meta()` on `final`, `expected_mice`, `max_cv`, `snapshot_root`. Add `elisa-qc-template` to CLAUDE.md's template table.

- [ ] **Step 4: Run, verify pass** -- same -> PASS.

- [ ] **Step 5: Commit**

```bash
git add R inst tests/testthat man NAMESPACE CLAUDE.md
git commit -m "feat(elisa): QC report with outlier plots and elisa_prep_run()"
```

---

### Task 7: `histo_` domain -- reader and ratio helper

**Files:**
- Modify: `CLAUDE.md` (add `histo_` -- histology scoring import, QC and export -- to the domain list, first)
- Create: `R/histo_read.R`, `R/histo_calc.R`
- Test: `tests/testthat/test-histo_read.R`, `tests/testthat/test-histo_calc.R`

**Interfaces:**
- Consumes: `parse_curated_dates_()` (Task 4); fixture `histo_curated.xlsx` (Task 4: `values` columns `mouse_ID, ID, gloms, crescents, necrosis`; `M1 85 22 15`, `M2 68 3 1`, `M3` all empty, `M4 70 80 2` (crescents > gloms)).
- Produces: `histo_read_curated(path, tissue = NA_character_)` -> tibble `mouse_ID, tissue, collected_at, metric, value`; non-numeric metric cells error like Task 4; ignored columns: any column named `ID` or ending in `_ID` other than `mouse_ID` (*plan's call, from the spec's "ID slide label"*); an all-empty row yields `NA` values (kept).
- Produces: `histo_calc_ratio(data, numerator, denominator, name, scale = 100)` -> `data` with extra rows `metric = name`, `value = numerator / denominator * scale` per `mouse_ID x tissue x collected_at`; `NA` when either is `NA` or the denominator is 0; errors if `name` already exists or a metric is absent.

- [ ] **Step 1: Write failing tests**

```r
test_that("histo_read_curated returns long metrics and keeps empty mice", {
  d <- histo_read_curated(test_path("../fixtures/histo_curated.xlsx"), tissue = "kidney")
  expect_named(d, c("mouse_ID", "tissue", "collected_at", "metric", "value"))
  expect_setequal(unique(d$metric), c("gloms", "crescents", "necrosis"))
  expect_true(all(is.na(d$value[d$mouse_ID == "M3"])))
})
test_that("histo_calc_ratio adds a percentage metric", {
  d <- histo_read_curated(test_path("../fixtures/histo_curated.xlsx"), tissue = "kidney") |>
    histo_calc_ratio("crescents", "gloms", "pct_crescents")
  expect_equal(d$value[d$mouse_ID == "M1" & d$metric == "pct_crescents"], 22 / 85 * 100)
  expect_true(is.na(d$value[d$mouse_ID == "M3" & d$metric == "pct_crescents"]))
})
```

- [ ] **Step 2: Run, verify failure** -- `Rscript -e 'devtools::test(filter = "histo")'` -> FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run, verify pass** -- same -> PASS.
- [ ] **Step 5: Commit**

```bash
git add CLAUDE.md R/histo_read.R R/histo_calc.R tests/testthat man NAMESPACE
git commit -m "feat(histo): histo_ domain with curated reader and ratio helper"
```

---

### Task 8: Histology prep save / load and `db_export_histo()`

**Files:**
- Create: `R/histo_prep.R`, `R/db_export_histo.R`
- Test: `tests/testthat/test-histo_prep.R`, `tests/testthat/test-db_export_histo.R`

**Interfaces:**
- Produces: `histo_prep_save(data, path, assay, group_col = NULL, batch_col = NULL, expected_mice = NULL, rules = character())` and `histo_prep_load(path)` -- same contract as Task 5 (info has `rules`). Save warns (doesn't error) when one `mouse_ID x tissue x collected_at x metric` occurs twice, listing them.
- Produces: `db_export_histo(con, prep, exclusions = NULL, data_dir = "data", overwrite = NULL)` -- same order as `db_export_elisa()`, duplicate key `mouse_ID x tissue x collected_at x metric`, table `histo_measurements`, domain `"histo"`. Rows with `NA` value are written (an empty score is information).

- [ ] **Step 1: Write failing tests** mirroring Task 5's list for histology: paths, staleness, write + replace + supersede, decline, exclusion row recorded, unmatched exclusion error, duplicate after date resolution error, and `pct_crescents` rows from `prepare` stored.
- [ ] **Step 2: Run, verify failure** -- `Rscript -e 'devtools::test(filter = "histo_prep|db_export_histo")'` -> FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run, verify pass** -- same -> PASS.
- [ ] **Step 5: Commit**

```bash
git add R tests/testthat man NAMESPACE
git commit -m "feat(histo): prep save/load and replaceable database export"
```

---

### Task 9: Histology QC report, `histo_prep_run()`, script template

**Files:**
- Create: `R/histo_qc.R`, `inst/rmarkdown/templates/histo-qc-template/{template.yaml,skeleton/skeleton.Rmd}`, `inst/scripts/histo_prep_template.R`
- Modify: `R/report_knit.R` (`report_knit_histo_qc()`), `R/histo_prep.R` (`histo_prep_run()`), `CLAUDE.md` (template table)
- Test: `tests/testthat/test-histo_qc.R`, `test-report_knit.R`, `test-histo_prep.R`

**Interfaces:**
- Consumes: `R/qc_table.R` helpers (Task 6).
- Produces: internal `histo_qc_rules_(data, rules)` -> tibble `rule, mouse_ID, tissue, timepoint, values` for rows where the rule is `FALSE` (`NA` = not a violation). Pivot wide by metric, evaluate each rule with `eval(parse(text = rule), wide, baseenv())`; a rule naming an absent metric errors naming the rule.
- Produces: `report_knit_histo_qc(prep_path, exclusions = NULL)`; `histo_prep_run(final, path, assay, tissue = NA_character_, prepare = identity, mice = NULL, exclusions = NULL, con = NULL, group_col = NULL, batch_col = NULL, ..., data_dir = "data", snapshot_root = NULL, overwrite = NULL)` mirroring `elisa_prep_run()`.
- Report sections: Overview, Completeness, Exclusions, Empty values, Sanity rules, Outliers (per metric), Flagged samples.

- [ ] **Step 1: Write failing tests** -- `histo_qc_rules_(d, "crescents <= gloms")` returns exactly `M4`; `"foo <= gloms"` errors `"foo <= gloms"`; report renders `qc/histo_curated_qc.html` containing `"Sanity rules"`, `"M4"`, `"Outliers"`; `histo_prep_run(FALSE, ...)` returns the html path.
- [ ] **Step 2: Run, verify failure** -- `Rscript -e 'devtools::test(filter = "histo|report_knit")'` -> FAIL.
- [ ] **Step 3: Implement**, including `inst/scripts/histo_prep_template.R` (like the ELISA one, with `prepare <- function(d) histo_calc_ratio(d, "crescents", "gloms", "pct_crescents")` and `rules = c("crescents <= gloms", "necrosis <= gloms")`).
- [ ] **Step 4: Run, verify pass** -- same -> PASS.
- [ ] **Step 5: Commit**

```bash
git add R inst tests/testthat man NAMESPACE CLAUDE.md
git commit -m "feat(histo): QC report with sanity rules and histo_prep_run()"
```

---

### Task 10: Docs and full check

**Files:**
- Modify: `CLAUDE.md` (Data structures: samples keyed by date; `elisa_read_curated()`, ELISA/histo prep loop, `db_export_elisa/histo()`, removal of `db_write_elisa/histo()`; `db_` bullet), `docs/superpowers/specs/2026-07-17-facs-db-schema-design.md` (samples, ELISA/histo tables, `import_exclusions.collected_at`), `docs/superpowers/specs/2026-11-18-elisa-histo-prep-qc-export-design.md` (Status: approved 2026-11-18)

- [ ] **Step 1: Update the docs above.**
- [ ] **Step 2: Run** `Rscript -e 'devtools::document(); devtools::check()'` -- Expected: `0 errors | 0 warnings`, NOTEs only from CLAUDE.md's known list.
- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md docs man NAMESPACE
git commit -m "docs: ELISA/histology pipeline and dated samples"
```

### After merge (with the user, outside the package)

Snapshot the live DB, recreate `~/databases/MPO_mRNA/MPO_mRNA.duckdb`, rerun the 26-1 FACS panel scripts with `final <- TRUE`, reshape `anti-mpo-titer_26-1.xlsx` into the curated layout (keep the original), add ELISA/histo prep scripts under the experiment's `src/`.
