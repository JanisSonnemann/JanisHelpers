# DB schema: wsp_root + import_log Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Drop FACS calculation inputs from the DuckDB schema, record each stain's workspace root, and log every import in a new `import_log` table.

**Architecture:** All DDL lives in `create_schema_()` (`R/db_connect.R`); no migration (no real database exists). Writers in `R/db_write_measurements.R` validate, then write data + log row in one `DBI::dbWithTransaction()`. A shared internal `write_import_log_()` is reused later by the Stage 3 export.

**Tech Stack:** R, DBI, duckdb, dbplyr, testthat 3 (`withr` for cleanup).

**Spec:** `docs/superpowers/specs/2026-11-12-db-schema-import-log-design.md`

## Global Constraints

- Work on branch `worktree-facs-pipeline-workflow`, in the worktree root. Never `cd` to the main checkout.
- CLAUDE.md rules apply: base pipe `|>`, every external call namespaced `pkg::fn()`, roxygen prose plain ASCII (no `\uXXXX`), `@param`/`@returns`/`@export` on exported functions.
- Run tests with `Rscript -e 'testthat::test_local(filter = "<name>", reporter = "summary", stop_on_failure = FALSE)'`; full suite without `filter`. Run `Rscript -e 'devtools::document()'` after roxygen changes.
- Test DBs come from `local_test_db()` / `seed_minimal_fixture()` in `tests/testthat/helper-db.R` (fixture: experiment `25-7`, mouse `25-7-1`, sample spleen, assays `overview` (facs) and `anti-MPO` (elisa)).
- Target `devtools::check()`: 0 errors, 0 warnings.
- `import_log.package_version` = `as.character(utils::packageVersion("JanisHelpers"))`.

## Review Focus

1. **Regex characters in a root** (`cd4+_26-2`): the filename check must match the root literally -- `cd44_26-2_20261112.wsp` must be rejected for root `cd4+_26-2`. Test in Task 2.
2. **Second sample of the same workspace**: a per-sample `db_write_facs()` call for another mouse/tissue with the same root and same `source_file` must succeed. Test in Task 2.
3. **New version arriving on a different sample**: root `overview_25-7` holds `_20261112` data for spleen; writing `_20261120` for a *new* sample must also be rejected (version check is root-wide, not per stain -- stricter reading of the spec's "one workspace version per root"). Test in Task 2.
4. **ELISA plate spanning two experiments**: log `experiment_code` = distinct codes sorted, joined `", "`. Test in Task 3.
5. **Opening a database file built by the old schema**: `db_connect()` errors with a "recreate" message instead of silently keeping old tables. Test in Task 1.

---

### Task 1: New schema + adapt everything that touches dropped columns

**Files:**
- Modify: `R/db_connect.R` (`create_schema_()`, `db_connect()`)
- Modify: `R/db_write_measurements.R` (`db_write_facs()`)
- Modify: `R/db_write_dimensions.R` (`db_write_samples()`)
- Modify: `R/db_query.R` (`db_query_facs()`)
- Test: `tests/testthat/test-db_connect.R`, `test-db_write_measurements.R`, `test-db_write_dimensions.R`, `test-db_query.R`

**Interfaces:**
- Produces: tables `facs_stains(facs_stain_id, sample_id, assay_id, wsp_root VARCHAR NOT NULL, stain_date, UNIQUE(sample_id, assay_id))`; `facs_measurements` without FK on `facs_stain_id`; `import_log` exactly as in the spec's DDL (`seq_import_log`); no `sample_digestions`.
- Produces: `db_write_facs(con, data, mouse_id, tissue, assay_name, wsp_root, source_file, stain_date = NA)` -- returns invisibly the number of measurement rows inserted. `wsp_root` validated as a single non-empty string.
- Produces: internal `check_string_arg_(x, arg)` -- errors unless `x` is a single non-`NA`, non-empty character; message contains the backticked arg name. Replaces the inline `source_file` check from step 3 (keep its message wording for `source_file` via an optional `hint` argument or by passing the example in).

- [ ] **Step 1: Write the failing tests**

`test-db_connect.R`: replace `"sample_digestions"` in the table list with `"import_log"`, and add:

```r
test_that("schema has wsp_root, no calc-input columns, no sample_digestions", {
  con <- local_test_db()
  cols <- DBI::dbListFields(con, "facs_stains")
  expect_true("wsp_root" %in% cols)
  expect_false(any(c("vol_stained", "count_method", "vol_resuspended", "vol_measured",
    "bead_volume_added", "bead_concentration", "bead_population_path") %in% cols))
  expect_false("sample_digestions" %in% DBI::dbListTables(con))
  expect_setequal(DBI::dbListFields(con, "import_log"), c(
    "import_id", "domain_id", "experiment_code", "import_unit", "source_file",
    "source_version", "qc_report", "qc_passed", "rows_written",
    "package_version", "supersedes_import_id", "imported_at"))
})

test_that("db_connect errors on a file created by the old schema", {
  path <- tempfile(fileext = ".duckdb")
  old <- DBI::dbConnect(duckdb::duckdb(), path)
  DBI::dbExecute(old, "CREATE TABLE facs_stains (facs_stain_id INTEGER, count_method VARCHAR)")
  DBI::dbDisconnect(old, shutdown = TRUE)
  expect_error(db_connect(path), "older schema")
})

test_that("a root's measurements and stains can be replaced in one transaction", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  write_overview_ <- function(src) db_write_facs(
    con, tibble::tibble(PopulationFullPath = "root/CD4+", Population = "CD4+",
      metric = "Count", value = 1),
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    wsp_root = "overview_25-7", source_file = src)
  write_overview_("25-7/facs/overview_25-7_20261112.wsp")
  expect_no_error(DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "DELETE FROM facs_measurements WHERE facs_stain_id IN
      (SELECT facs_stain_id FROM facs_stains WHERE wsp_root = 'overview_25-7')")
    DBI::dbExecute(con, "DELETE FROM facs_stains WHERE wsp_root = 'overview_25-7'")
    write_overview_("25-7/facs/overview_25-7_20261120.wsp")
  }))
  expect_equal(DBI::dbGetQuery(con, "SELECT source_file FROM facs_measurements")$source_file,
    "25-7/facs/overview_25-7_20261120.wsp")
})
```

`test-db_write_measurements.R`: in every existing `db_write_facs()` call, remove `count_method`, `vol_*`, `bead_*` arguments, add `wsp_root = "overview_25-7"` (or `"treg_25-7"` for the `treg` call) and change `source_file` to `"25-7/facs/<root>_20261112.wsp"`. Delete the bead/mixed-fields test (the CHECK no longer exists). Rename the volumetric test to `"db_write_facs creates a stain with wsp_root and inserts measurements"` and replace `expect_equal(stains$count_method, "volumetric")` with `expect_equal(stains$wsp_root, "overview_25-7")`. Add `wsp_root` validation, mirroring the `source_file` test:

```r
test_that("db_write_facs requires a single non-empty wsp_root", {
  # same setup/args as the source_file test, with source_file set
  expect_error(do.call(db_write_facs, args), "wsp_root")
  expect_error(do.call(db_write_facs, c(args, wsp_root = NA_character_)), "wsp_root")
  expect_error(do.call(db_write_facs, c(args, wsp_root = "")), "wsp_root")
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM facs_stains")$n, 0)
})
```

`test-db_write_dimensions.R`: replace `"db_write_samples inserts samples and digestions where weights are given"` with `"db_write_samples inserts samples and ignores weight columns"` -- pass `total_weight`/`facs_weight`/`vol_total`, expect the `samples` rows and no error.

`test-db_query.R`: update the `db_write_facs()` call as above; add `expect_equal(result$wsp_root, "overview_25-7")` and `expect_false("count_method" %in% names(result))`.

- [ ] **Step 2: Run to verify they fail**

Run: `Rscript -e 'testthat::test_local(filter = "db_", reporter = "summary", stop_on_failure = FALSE)'`
Expected: failures on missing `import_log`/`wsp_root`, unused-argument errors for `wsp_root`.

- [ ] **Step 3: Implement**

- `create_schema_()`: DDL exactly as in the spec (facs_stains, facs_measurements without the FK plus the `-- no FK` comment, import_log with sequence; remove `sample_digestions`).
- `db_connect()`: **before** `create_schema_()`, if table `facs_stains` exists and `"wsp_root"` is not in `DBI::dbListFields(con, "facs_stains")`, disconnect and `stop()` with a message containing `"older schema"` and telling the user to recreate the file. (Checking first keeps `create_schema_()` from failing with an unrelated error on an old file.)
- `db_write_facs()`: new signature; stain INSERT writes `(sample_id, assay_id, wsp_root, stain_date)`; roxygen updated (drop removed params, document `wsp_root`).
- `db_write_samples()`: drop the three weight columns from `fill_missing_cols_()` and the `sample_digestions` INSERT; roxygen says weight columns are ignored.
- `db_query_facs()`: select `wsp_root` instead of the seven dropped columns.

- [ ] **Step 4: Run tests, then the full suite**

Run the Step 2 command, then the full suite. Expected: 0 failures.

- [ ] **Step 5: Document and commit**

```bash
Rscript -e 'devtools::document()'
git add R/ man/ tests/testthat/
git commit -m "feat(db): drop FACS calc inputs, add wsp_root and import_log tables"
```

---

### Task 2: `db_write_facs()` root and version guards

**Files:**
- Modify: `R/db_write_measurements.R` (`db_write_facs()`)
- Test: `tests/testthat/test-db_write_measurements.R`

**Interfaces:**
- Consumes: Task 1's `db_write_facs()` signature and `facs_stains.wsp_root`.
- Produces: no new interface; `db_write_facs()` errors (writing nothing) on the four conditions below.

Guards, all checked before the transaction writes anything:
1. `basename(source_file)` must equal `paste0(wsp_root, "_", <8 digits>, ".wsp")` -- match the root **literally** (e.g. `startsWith()` + a check that the remainder matches `^_[0-9]{8}\.wsp$`, or escape the root). Error mentions both `wsp_root` and the filename.
2. A stain exists for the same experiment (sample -> subject -> experiment) and `assay_id` with a different `wsp_root` -> error naming both roots.
3. A stain exists with the same `wsp_root` but a different experiment or `assay_id` -> error.
4. Any measurement under any stain with this `wsp_root` has a `source_file` different from this call's -> error saying replacing a root is the Stage 3 export's job (root-wide, see Review Focus 3).

- [ ] **Step 1: Write the failing tests**

Add a second mouse to the setup where needed: `db_write_subjects(con, tibble::tibble(mouse_id = "25-7-2", experiment_code = "25-7"))`, `db_write_samples(con, tibble::tibble(mouse_id = "25-7-2", tissue = "spleen"))`; for a second experiment register `25-8` / mouse `25-8-1` / spleen. Each test asserts the error *and* that `facs_stains`/`facs_measurements` row counts are unchanged.

- `"db_write_facs rejects a source_file that doesn't match wsp_root"`: root `overview_25-7`, file `"25-7/facs/treg_25-7_20261112.wsp"` -> error matching `"wsp_root"`; file `"25-7/facs/overview_25-7.wsp"` (no date) -> error.
- `"db_write_facs matches the root literally"`: register facs assay `cd4+`; root `"cd4+_25-7"` with file `"25-7/facs/cd4+_25-7_20261112.wsp"` succeeds; with file `"25-7/facs/cd44_25-7_20261112.wsp"` errors.
- `"db_write_facs rejects a second root for the same experiment x panel"`: write mouse 25-7-1 under `overview_25-7`, then mouse 25-7-2 / `overview` under `overview2_25-7` (file `overview2_25-7_20261112.wsp`) -> error matching `"overview_25-7"`.
- `"db_write_facs rejects one root across two panels or experiments"`: after writing `overview_25-7`, writing mouse 25-7-2 with `assay_name = "treg"` (registered) and `wsp_root = "overview_25-7"` errors; writing mouse 25-8-1 with `overview` and `overview_25-7` errors.
- `"db_write_facs rejects a second workspace version of a root"`: write 25-7-1 from `_20261112`; then 25-7-1 from `_20261120` errors matching `"Stage 3|export"`; then 25-7-2 (new sample) from `_20261120` also errors.
- `"db_write_facs accepts further samples of the same workspace"`: 25-7-1 and 25-7-2 from the same root and same file both succeed; `facs_stains` has 2 rows.
- `"every facs measurement has a matching stain"`: after the previous writes, `SELECT COUNT(*) FROM facs_measurements m LEFT JOIN facs_stains s USING (facs_stain_id) WHERE s.facs_stain_id IS NULL` is 0.

- [ ] **Step 2: Run to verify they fail**

Run: `Rscript -e 'testthat::test_local(filter = "db_write_measurements", reporter = "summary", stop_on_failure = FALSE)'`
Expected: the guard tests fail (no error thrown); the "accepts further samples" and "matching stain" tests pass.

- [ ] **Step 3: Implement the four guards in `db_write_facs()`**

One SQL query per guard against `facs_stains` joined to `samples`/`subjects` (for experiment) and `facs_measurements` (for guard 4), after resolving `sample_id`/`assay_id` and before `DBI::dbWithTransaction()`. Update roxygen `@details` to list the guards.

- [ ] **Step 4: Run tests, then the full suite** -- Expected: 0 failures.

- [ ] **Step 5: Document and commit**

```bash
Rscript -e 'devtools::document()'
git add R/db_write_measurements.R man/db_write_facs.Rd tests/testthat/test-db_write_measurements.R
git commit -m "feat(db): db_write_facs guards one root and one version per experiment x panel"
```

---

### Task 3: `write_import_log_()` + ELISA logging

**Files:**
- Create: `R/db_import_log.R`
- Modify: `R/db_write_measurements.R` (`db_write_elisa()`), `DESCRIPTION` (add `utils` to Imports, alphabetical: after `tools`)
- Test: `tests/testthat/test-db_import_log.R`, `tests/testthat/test-db_write_measurements.R`, `tests/testthat/test-db_query.R`

**Interfaces:**
- Consumes: Task 1's `import_log` table, `check_string_arg_()`; existing `lookup_sample_ids_()`.
- Produces: `write_import_log_(con, domain, experiment_code, import_unit, source_file, rows_written, source_version = NA, qc_report = NA_character_, qc_passed = NA, supersedes_import_id = NA_integer_)` -> integer `import_id`. Resolves `domain` (`"facs"`/`"elisa"`/`"histo"`) to `domain_id`; fills `package_version`. Opens no transaction of its own.
- Produces: `db_write_elisa(con, data, assay_name, source_file)` -- `source_file` now required; writes data + one log row in one transaction; returns invisibly rows inserted.

- [ ] **Step 1: Write the failing tests**

`test-db_import_log.R`:
- `"write_import_log_ writes one row with package version and returns its id"`: call with `domain = "facs"`, `experiment_code = "25-7"`, `import_unit = "overview_25-7"`, `source_file = "25-7/facs/overview_25-7_20261112.wsp"`, `rows_written = 5L`, `source_version = as.Date("2026-11-12")`, `qc_report = "25-7/facs/qc/2026-11-12_qc.html"`, `qc_passed = TRUE`. Expect returned id equals the row's `import_id`; `package_version == as.character(utils::packageVersion("JanisHelpers"))`; domain joins to `"facs"`; `imported_at` not `NA`.
- `"write_import_log_ links supersedes_import_id"`: second call with `supersedes_import_id = <first id>` stores it.

`test-db_write_measurements.R` (ELISA block; `seed_two_mouse_plate_()` / `elisa_plate_()` already exist): pass `source_file = "25-7/elisa/Results_MPO.xlsx"` in every call, then add:
- `"db_write_elisa writes one import_log row per call"`: after one call, log has 1 row: domain `elisa`, `experiment_code == "25-7"`, `import_unit == "Results_MPO.xlsx"`, `rows_written == 4`, `qc_passed`/`qc_report`/`source_version`/`supersedes_import_id` all `NA`. A second identical call adds a row with `rows_written == 0`.
- `"db_write_elisa logs every experiment on a mixed plate"`: register experiment `25-8`, mouse `25-8-1`, sample serum; plate rows for `25-7-1` and `25-8-1` -> `experiment_code == "25-7, 25-8"`.
- `"db_write_elisa writes no log row when the data insert fails"`: plate with an unresolved mouse -> error, `import_log` empty.
- `"db_write_elisa requires source_file"`: omitted / `NA` / `""` -> error matching `"source_file"`.

`test-db_query.R`: add `source_file` to the `db_write_elisa()` call.

- [ ] **Step 2: Run to verify they fail**

Run: `Rscript -e 'testthat::test_local(filter = "db_", reporter = "summary", stop_on_failure = FALSE)'`
Expected: `write_import_log_` not found; no log rows written.

- [ ] **Step 3: Implement**

- `R/db_import_log.R`: `write_import_log_()` as in Interfaces (no roxygen export; a plain `#` comment block).
- `db_write_elisa()`: validate `source_file` with `check_string_arg_()`; compute `experiment_code` from the resolved sample ids (distinct, sorted, `paste(collapse = ", ")`); wrap the INSERT and `write_import_log_(con, "elisa", experiment_code, basename(source_file), source_file, n)` in `DBI::dbWithTransaction()`. Roxygen: `source_file` required, mentions the log row.
- `DESCRIPTION`: `utils` in Imports.

- [ ] **Step 4: Run tests, then the full suite** -- Expected: 0 failures.

- [ ] **Step 5: Document and commit**

```bash
Rscript -e 'devtools::document()'
git add R/db_import_log.R R/db_write_measurements.R DESCRIPTION man/ tests/testthat/
git commit -m "feat(db): import_log helper; db_write_elisa logs every import"
```

---

### Task 4: `db_write_histo()` whole-table writer with logging

**Files:**
- Modify: `R/db_write_measurements.R` (`db_write_histo()`)
- Test: `tests/testthat/test-db_write_measurements.R`, `tests/testthat/test-db_query.R`

**Interfaces:**
- Consumes: `lookup_sample_ids_()`, `write_import_log_()`, `check_string_arg_()`.
- Produces: `db_write_histo(con, data, assay_name = NA_character_, source_file)` -- `data` has `metric`, `value`, `mouse_id`, `tissue`; returns invisibly rows inserted; one `import_log` row per call (domain `histo`, `import_unit = basename(source_file)`), same transaction.

- [ ] **Step 1: Write the failing tests**

Rewrite the existing histo tests to the new signature (columns `mouse_id`, `tissue` in `data`; `source_file = "25-7/histo/scores.xlsx"`), keeping their intent (no-assay insert, unknown `assay_name` errors). Add:
- `"db_write_histo writes a whole table across mice"`: with `seed_two_mouse_plate_()`, rows for 25-7-1/serum and 25-7-2/serum -> each row lands on the right sample (join back through `samples`/`subjects`).
- `"db_write_histo lists every unresolved mouse_id x tissue and writes nothing"`: rows with `25-7-9`/serum and `NA`/kidney -> message contains both `'25-7-9 x serum'` and `'NA x kidney'`; `histo_measurements` and `import_log` empty.
- `"db_write_histo requires mouse_id and tissue columns"`: missing `tissue` errors matching `"tissue"`; `mouse_ID` instead of `mouse_id` errors matching `"mouse_ID"`.
- `"db_write_histo writes one import_log row per call"`: domain `histo`, `import_unit == "scores.xlsx"`, `qc_passed` `NA`.

Update the `db_write_histo()` call in `test-db_query.R`.

- [ ] **Step 2: Run to verify they fail** -- filter `"db_"`; expected: unused/missing-argument errors.

- [ ] **Step 3: Implement**

Reuse `db_write_elisa()`'s column check (extract it into an internal `check_mouse_tissue_cols_(data, fn_name)` used by both) and its transaction + log pattern. Roxygen updated.

- [ ] **Step 4: Run tests, then the full suite** -- Expected: 0 failures.

- [ ] **Step 5: Document and commit**

```bash
Rscript -e 'devtools::document()'
git add R/db_write_measurements.R man/db_write_histo.Rd tests/testthat/
git commit -m "feat(db): db_write_histo writes a whole table and logs the import"
```

---

### Task 5: Docs, check, tick the plan

**Files:**
- Modify: `docs/superpowers/specs/2026-07-17-facs-db-schema-design.md`, `CLAUDE.md`, `docs/superpowers/plans/2026-11-12-db-setup.md`

- [ ] **Step 1: Update the schema spec** -- replace the `facs_stains`, `facs_measurements` and `sample_digestions` DDL with the new versions, add `import_log`, and add a dated note at the top: "2026-11-12: calc inputs removed, `wsp_root` + `import_log` added, `facs_measurements` FK dropped -- see `2026-11-12-db-schema-import-log-design.md`."
- [ ] **Step 2: Update CLAUDE.md** -- `db_` data-structure line: add "every import writes an `import_log` row (FACS: from the Stage 3 export); `db_write_elisa()`/`db_write_histo()` take whole tables with `mouse_id`/`tissue` columns". Add `utils` to the `Current Imports` line in DESCRIPTION order.
- [ ] **Step 3: Tick step 4b** in `2026-11-12-db-setup.md`.
- [ ] **Step 4: Full check**

Run: `Rscript -e 'devtools::check(error_on = "never", args = "--no-manual")'`
Expected: `0 errors | 0 warnings`; NOTEs only of the kinds listed in CLAUDE.md "Known check output".

- [ ] **Step 5: Commit**

```bash
git add CLAUDE.md docs/superpowers/
git commit -m "docs: schema spec, CLAUDE.md and db-setup plan for step 4b"
```
