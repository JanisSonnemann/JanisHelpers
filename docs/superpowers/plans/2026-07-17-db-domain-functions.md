# db_ Domain Functions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the `db_` domain -- R functions to connect to, write to, and
query the DuckDB lab-measurements schema specified in
`docs/superpowers/specs/2026-07-17-facs-db-schema-design.md`.

**Architecture:** One connect/schema-init function, one file of dimension-table
writers (experiments/subjects/samples/assays), one file of fact-table writers
(facs/elisa/histo measurements, each resolving natural keys to surrogate IDs
internally), and one file of query helpers that join fact tables back to
natural keys (`mouse_id`/`tissue`/`experiment_code`/`assay_name`) so callers
never see surrogate integer IDs. Every write is idempotent via the schema's
`UNIQUE` constraints and `INSERT ... ON CONFLICT DO NOTHING`.

**Tech Stack:** R, DBI, duckdb, dplyr, tibble, testthat (edition 3).

## Global Constraints

- File/function naming: `db_verb()` in files named `db_verb.R` (existing
  package convention, see `CLAUDE.md`).
- Every exported function: `@param`, `@returns`, `@export` roxygen tags; run
  `devtools::document()` after any roxygen change.
- Every external package call inside an exported function body is namespaced
  `pkg::fn()` -- no bare calls.
- Base pipe `|>`, not `%>%`.
- Functions that write data primarily for a side effect return
  `invisible(...)`.
- `devtools::check()` before every commit, target 0 errors / 0 warnings (new
  NOTEs for bare column names inside `dplyr` verbs are expected and match the
  existing pattern documented in `CLAUDE.md`'s "Known check output" section --
  Task 10 folds the new NOTEs into that section).
- Requires the `duckdb` R package version with `INSERT ... ON CONFLICT`
  support (current CRAN releases; this has been present for several years).
- Migrating existing `data/MPOmRNA.db` (SQLite) data into this schema is
  explicitly out of scope for this plan.

---

## File Overview

| File | Responsibility |
|---|---|
| `R/db_connect.R` | `db_connect()` + internal schema creation (all `CREATE TABLE`/`CREATE SEQUENCE`/domain seeding) |
| `R/db_write_dimensions.R` | `db_write_experiment()`, `db_write_subjects()`, `db_write_samples()`, `db_write_assay()` |
| `R/db_write_measurements.R` | `db_write_facs()`, `db_write_elisa()`, `db_write_histo()` + shared internal lookup helpers |
| `R/db_query.R` | `db_query_facs()`, `db_query_elisa()`, `db_query_histo()` |
| `tests/testthat/helper-db.R` | Shared test fixtures (`local_test_db()`, `seed_minimal_fixture()`) |

---

### Task 1: Database connection and schema creation

**Files:**
- Modify: `DESCRIPTION` (add `DBI`, `duckdb` to `Imports`)
- Create: `R/db_connect.R`
- Create: `tests/testthat/helper-db.R`
- Test: `tests/testthat/test-db_connect.R`

**Interfaces:**
- Produces: `db_connect(path)` returns a `DBI` connection (class
  `duckdb_connection`) with the full schema created and the three known
  domains (`"facs"`, `"histo"`, `"elisa"`) seeded into `domains`. Safe to call
  repeatedly on the same path.
- Produces (test helper): `local_test_db(env = parent.frame())` returns a
  connection to a fresh temp-file database, registering cleanup via
  `withr::defer()` in the caller's test environment.

- [ ] **Step 1: Add `DBI` and `duckdb` to `DESCRIPTION` Imports**

Edit the `Imports:` block in `DESCRIPTION` to read (alphabetical, case-insensitive):

```
Imports:
    CytoML,
    DBI,
    diffcyt,
    dplyr,
    duckdb,
    fcexpr,
    flowCore,
    FlowSOM,
    flowWorkspace,
    ggplot2,
    glue,
    gt,
    gtsummary,
    janitor,
    parallel,
    purrr,
    readxl,
    rmarkdown,
    rstatix,
    stats,
    stringr,
    SummarizedExperiment,
    tibble,
    tidyr,
    tools,
    uwot,
    xml2,
    xfun
```

- [ ] **Step 2: Write the test helper**

Create `tests/testthat/helper-db.R`:

```r
local_test_db <- function(env = parent.frame()) {
  path <- tempfile(fileext = ".duckdb")
  con <- db_connect(path)
  withr::defer(DBI::dbDisconnect(con, shutdown = TRUE), envir = env)
  con
}
```

- [ ] **Step 3: Write the failing test**

Create `tests/testthat/test-db_connect.R`:

```r
test_that("db_connect creates all tables and seeds domains", {
  con <- local_test_db()

  expect_true(all(c(
    "domains", "experiments", "subjects", "samples", "sample_digestions",
    "assays", "facs_stains", "facs_measurements", "elisa_measurements",
    "histo_measurements"
  ) %in% DBI::dbListTables(con)))

  domain_names <- DBI::dbGetQuery(
    con, "SELECT domain_name FROM domains ORDER BY domain_name"
  )$domain_name
  expect_equal(domain_names, c("elisa", "facs", "histo"))
})

test_that("db_connect is idempotent", {
  path <- tempfile(fileext = ".duckdb")
  con1 <- db_connect(path)
  DBI::dbDisconnect(con1, shutdown = TRUE)

  con2 <- db_connect(path)
  withr::defer(DBI::dbDisconnect(con2, shutdown = TRUE))

  domain_count <- DBI::dbGetQuery(con2, "SELECT COUNT(*) AS n FROM domains")$n
  expect_equal(domain_count, 3)
})
```

- [ ] **Step 4: Run test to verify it fails**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-db_connect.R')"`
Expected: FAIL -- `could not find function "db_connect"`

- [ ] **Step 5: Implement `db_connect()`**

Create `R/db_connect.R`:

```r
#' Connect to the lab measurements database
#'
#' Opens (creating if necessary) a DuckDB database file and ensures the full
#' dimension/fact schema exists, creating any missing tables, sequences, and
#' seed rows. Safe to call repeatedly -- existing tables and data are left
#' untouched.
#'
#' @param path Path to the DuckDB database file. Created if it doesn't exist.
#' @returns A `DBI` connection (class `duckdb_connection`).
#' @export
db_connect <- function(path) {
  con <- DBI::dbConnect(duckdb::duckdb(), path)
  create_schema_(con)
  con
}

create_schema_ <- function(con) {
  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_domains START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS domains (
      domain_id   INTEGER PRIMARY KEY DEFAULT nextval('seq_domains'),
      domain_name VARCHAR NOT NULL UNIQUE
    )
  ")
  DBI::dbExecute(con, "
    INSERT INTO domains (domain_name)
    SELECT * FROM (VALUES ('facs'), ('histo'), ('elisa')) AS v(domain_name)
    ON CONFLICT (domain_name) DO NOTHING
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_experiments START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS experiments (
      experiment_id   INTEGER PRIMARY KEY DEFAULT nextval('seq_experiments'),
      experiment_code VARCHAR NOT NULL UNIQUE,
      experiment_name VARCHAR,
      project         VARCHAR,
      description     VARCHAR,
      created_at      TIMESTAMP DEFAULT current_timestamp
    )
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_subjects START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS subjects (
      subject_id      INTEGER PRIMARY KEY DEFAULT nextval('seq_subjects'),
      mouse_id        VARCHAR NOT NULL UNIQUE,
      experiment_id   INTEGER NOT NULL REFERENCES experiments(experiment_id),
      mdc_id          VARCHAR,
      cage            VARCHAR,
      mouse_strain    VARCHAR,
      generation      VARCHAR,
      sex             VARCHAR,
      mouse_treatment VARCHAR,
      treatment_group VARCHAR,
      dob             DATE,
      start_date      DATE,
      bmt_date        DATE,
      end_date        DATE,
      created_at      TIMESTAMP DEFAULT current_timestamp
    )
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_samples START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS samples (
      sample_id    INTEGER PRIMARY KEY DEFAULT nextval('seq_samples'),
      subject_id   INTEGER NOT NULL REFERENCES subjects(subject_id),
      tissue       VARCHAR NOT NULL,
      collected_at DATE,
      created_at   TIMESTAMP DEFAULT current_timestamp,
      UNIQUE(subject_id, tissue)
    )
  ")

  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS sample_digestions (
      sample_id    INTEGER PRIMARY KEY REFERENCES samples(sample_id),
      total_weight DOUBLE,
      facs_weight  DOUBLE,
      vol_total    DOUBLE
    )
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_assays START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS assays (
      assay_id    INTEGER PRIMARY KEY DEFAULT nextval('seq_assays'),
      assay_name  VARCHAR NOT NULL,
      domain_id   INTEGER NOT NULL REFERENCES domains(domain_id),
      description VARCHAR,
      UNIQUE(assay_name, domain_id)
    )
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_facs_stains START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS facs_stains (
      facs_stain_id        INTEGER PRIMARY KEY DEFAULT nextval('seq_facs_stains'),
      sample_id            INTEGER NOT NULL REFERENCES samples(sample_id),
      assay_id             INTEGER NOT NULL REFERENCES assays(assay_id),
      vol_stained          DOUBLE,
      count_method         VARCHAR NOT NULL CHECK (count_method IN ('volumetric', 'bead')),
      vol_resuspended      DOUBLE,
      vol_measured         DOUBLE,
      bead_volume_added    DOUBLE,
      bead_concentration   DOUBLE,
      bead_population_path VARCHAR,
      stain_date           DATE,
      UNIQUE(sample_id, assay_id),
      CHECK (
        (count_method = 'volumetric'
          AND vol_resuspended IS NOT NULL AND vol_measured IS NOT NULL
          AND bead_volume_added IS NULL AND bead_concentration IS NULL)
        OR
        (count_method = 'bead'
          AND bead_volume_added IS NOT NULL AND bead_concentration IS NOT NULL
          AND vol_resuspended IS NULL AND vol_measured IS NULL)
      )
    )
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_facs_measurements START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS facs_measurements (
      facs_measurement_id  BIGINT PRIMARY KEY DEFAULT nextval('seq_facs_measurements'),
      facs_stain_id        INTEGER NOT NULL REFERENCES facs_stains(facs_stain_id),
      population_full_path VARCHAR NOT NULL,
      population           VARCHAR NOT NULL,
      metric               VARCHAR NOT NULL,
      value                DOUBLE,
      source_file          VARCHAR,
      imported_at          TIMESTAMP DEFAULT current_timestamp,
      UNIQUE(facs_stain_id, population_full_path, metric)
    )
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_elisa_measurements START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS elisa_measurements (
      elisa_measurement_id INTEGER PRIMARY KEY DEFAULT nextval('seq_elisa_measurements'),
      sample_id            INTEGER NOT NULL REFERENCES samples(sample_id),
      assay_id             INTEGER NOT NULL REFERENCES assays(assay_id),
      cytokine             VARCHAR NOT NULL,
      sample_id_raw        VARCHAR,
      replicate            INTEGER,
      value                DOUBLE,
      unit                 VARCHAR,
      result_status        VARCHAR,
      source_file          VARCHAR,
      imported_at          TIMESTAMP DEFAULT current_timestamp,
      UNIQUE(sample_id, assay_id, cytokine, sample_id_raw, replicate)
    )
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_histo_measurements START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS histo_measurements (
      histo_measurement_id INTEGER PRIMARY KEY DEFAULT nextval('seq_histo_measurements'),
      sample_id            INTEGER NOT NULL REFERENCES samples(sample_id),
      assay_id             INTEGER REFERENCES assays(assay_id),
      metric               VARCHAR NOT NULL,
      value                DOUBLE,
      source_file          VARCHAR,
      imported_at          TIMESTAMP DEFAULT current_timestamp,
      UNIQUE(sample_id, assay_id, metric)
    )
  ")

  invisible(con)
}
```

- [ ] **Step 6: Run `devtools::document()`**

Run: `Rscript -e "devtools::document()"`
Expected: `man/db_connect.Rd` created, no warnings.

- [ ] **Step 7: Run test to verify it passes**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-db_connect.R')"`
Expected: PASS, 0 failures.

- [ ] **Step 8: Commit**

```bash
git add DESCRIPTION R/db_connect.R tests/testthat/helper-db.R tests/testthat/test-db_connect.R man/db_connect.Rd
git commit -m "feat: add db_connect() and DuckDB schema creation"
```

---

### Task 2: `db_write_experiment()`

**Files:**
- Create: `R/db_write_dimensions.R`
- Test: `tests/testthat/test-db_write_dimensions.R`

**Interfaces:**
- Consumes: `db_connect(path)` from Task 1; `local_test_db()` from Task 1.
- Produces: `db_write_experiment(con, experiment_code, experiment_name = NA_character_, project = NA_character_, description = NA_character_)` returns invisibly the number of rows inserted (`0` or `1`).

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-db_write_dimensions.R`:

```r
test_that("db_write_experiment inserts a new experiment", {
  con <- local_test_db()

  n <- db_write_experiment(con, experiment_code = "25-7", experiment_name = "Test experiment")
  expect_equal(n, 1)

  row <- DBI::dbGetQuery(con, "SELECT * FROM experiments WHERE experiment_code = '25-7'")
  expect_equal(nrow(row), 1)
  expect_equal(row$experiment_name, "Test experiment")
})

test_that("db_write_experiment is idempotent", {
  con <- local_test_db()

  db_write_experiment(con, experiment_code = "25-7")
  n <- db_write_experiment(con, experiment_code = "25-7", experiment_name = "Different name")
  expect_equal(n, 0)

  row <- DBI::dbGetQuery(con, "SELECT * FROM experiments WHERE experiment_code = '25-7'")
  expect_equal(nrow(row), 1)
  expect_true(is.na(row$experiment_name))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-db_write_dimensions.R')"`
Expected: FAIL -- `could not find function "db_write_experiment"`

- [ ] **Step 3: Implement `db_write_experiment()`**

Create `R/db_write_dimensions.R`:

```r
#' Register an experiment
#'
#' Inserts a new experiment row identified by `experiment_code`. If an
#' experiment with that code already exists, the call is a no-op.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param experiment_code Character scalar, the experiment's natural key
#'   (e.g. `"25-7"`).
#' @param experiment_name Character scalar, optional.
#' @param project Character scalar, optional.
#' @param description Character scalar, optional.
#' @returns Invisibly, the number of rows inserted (`0` or `1`).
#' @export
db_write_experiment <- function(con, experiment_code, experiment_name = NA_character_,
                                 project = NA_character_, description = NA_character_) {
  n <- DBI::dbExecute(
    con,
    "INSERT INTO experiments (experiment_code, experiment_name, project, description)
     VALUES (?, ?, ?, ?)
     ON CONFLICT (experiment_code) DO NOTHING",
    params = list(experiment_code, experiment_name, project, description)
  )
  invisible(n)
}
```

- [ ] **Step 4: Run `devtools::document()`**

Run: `Rscript -e "devtools::document()"`
Expected: `man/db_write_experiment.Rd` created, no warnings.

- [ ] **Step 5: Run test to verify it passes**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-db_write_dimensions.R')"`
Expected: PASS, 0 failures.

- [ ] **Step 6: Commit**

```bash
git add R/db_write_dimensions.R tests/testthat/test-db_write_dimensions.R man/db_write_experiment.Rd
git commit -m "feat: add db_write_experiment()"
```

---

### Task 3: `db_write_subjects()`

**Files:**
- Modify: `R/db_write_dimensions.R`
- Modify: `tests/testthat/test-db_write_dimensions.R`

**Interfaces:**
- Consumes: `db_write_experiment()` from Task 2.
- Produces: `db_write_subjects(con, data)` where `data` has columns
  `mouse_id`, `experiment_code`, and optionally `mdc_id`, `cage`,
  `mouse_strain`, `generation`, `sex`, `mouse_treatment`, `treatment_group`,
  `dob`, `start_date`, `bmt_date`, `end_date`. Returns invisibly the number
  of rows inserted. Also produces the internal helper `fill_missing_cols_(data, cols)`.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-db_write_dimensions.R`:

```r
test_that("db_write_subjects inserts subjects linked to an experiment", {
  con <- local_test_db()
  db_write_experiment(con, experiment_code = "25-7")

  n <- db_write_subjects(con, tibble::tibble(
    mouse_id = c("25-7-1", "25-7-2"),
    experiment_code = "25-7",
    sex = c("f", "m")
  ))
  expect_equal(n, 2)

  rows <- DBI::dbGetQuery(con, "
    SELECT s.mouse_id, s.sex, e.experiment_code
    FROM subjects s JOIN experiments e ON e.experiment_id = s.experiment_id
    ORDER BY s.mouse_id
  ")
  expect_equal(rows$mouse_id, c("25-7-1", "25-7-2"))
  expect_equal(rows$sex, c("f", "m"))
  expect_equal(rows$experiment_code, c("25-7", "25-7"))
})

test_that("db_write_subjects is idempotent and fills missing optional columns", {
  con <- local_test_db()
  db_write_experiment(con, experiment_code = "25-7")

  db_write_subjects(con, tibble::tibble(mouse_id = "25-7-1", experiment_code = "25-7"))
  n <- db_write_subjects(con, tibble::tibble(mouse_id = "25-7-1", experiment_code = "25-7"))
  expect_equal(n, 0)

  row <- DBI::dbGetQuery(con, "SELECT * FROM subjects WHERE mouse_id = '25-7-1'")
  expect_equal(nrow(row), 1)
  expect_true(is.na(row$sex))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-db_write_dimensions.R')"`
Expected: FAIL -- `could not find function "db_write_subjects"`

- [ ] **Step 3: Implement `db_write_subjects()`**

Append to `R/db_write_dimensions.R`:

```r
fill_missing_cols_ <- function(data, cols) {
  missing <- setdiff(cols, names(data))
  data[missing] <- NA
  data
}

#' Register subjects (mice)
#'
#' Bulk-inserts subject rows. `data` must already be shaped to match the
#' `subjects` table: one row per mouse, with `mouse_id` and an
#' `experiment_code` identifying an experiment already registered via
#' [db_write_experiment()]. Rows whose `mouse_id` already exists are skipped.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param data A tibble with columns `mouse_id`, `experiment_code`, and
#'   optionally `mdc_id`, `cage`, `mouse_strain`, `generation`, `sex`,
#'   `mouse_treatment`, `treatment_group`, `dob`, `start_date`, `bmt_date`,
#'   `end_date`. Missing optional columns are filled with `NA`.
#' @returns Invisibly, the number of rows inserted.
#' @export
db_write_subjects <- function(con, data) {
  data <- fill_missing_cols_(data, c(
    "mouse_id", "experiment_code", "mdc_id", "cage", "mouse_strain",
    "generation", "sex", "mouse_treatment", "treatment_group",
    "dob", "start_date", "bmt_date", "end_date"
  ))
  duckdb::duckdb_register(con, "tmp_subjects", data)
  on.exit(duckdb::duckdb_unregister(con, "tmp_subjects"), add = TRUE)

  n <- DBI::dbExecute(con, "
    INSERT INTO subjects (
      mouse_id, experiment_id, mdc_id, cage, mouse_strain, generation,
      sex, mouse_treatment, treatment_group, dob, start_date, bmt_date, end_date
    )
    SELECT
      s.mouse_id, e.experiment_id, s.mdc_id, s.cage, s.mouse_strain, s.generation,
      s.sex, s.mouse_treatment, s.treatment_group,
      CAST(s.dob AS DATE), CAST(s.start_date AS DATE),
      CAST(s.bmt_date AS DATE), CAST(s.end_date AS DATE)
    FROM tmp_subjects s
    JOIN experiments e ON e.experiment_code = s.experiment_code
    ON CONFLICT (mouse_id) DO NOTHING
  ")
  invisible(n)
}
```

- [ ] **Step 4: Run `devtools::document()`**

Run: `Rscript -e "devtools::document()"`
Expected: `man/db_write_subjects.Rd` created, no warnings.

- [ ] **Step 5: Run test to verify it passes**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-db_write_dimensions.R')"`
Expected: PASS, 0 failures.

- [ ] **Step 6: Commit**

```bash
git add R/db_write_dimensions.R tests/testthat/test-db_write_dimensions.R man/db_write_subjects.Rd
git commit -m "feat: add db_write_subjects()"
```

---

### Task 4: `db_write_samples()`

**Files:**
- Modify: `R/db_write_dimensions.R`
- Modify: `tests/testthat/test-db_write_dimensions.R`

**Interfaces:**
- Consumes: `db_write_subjects()` from Task 3; `fill_missing_cols_()` from Task 3.
- Produces: `db_write_samples(con, data)` where `data` has columns
  `mouse_id`, `tissue`, and optionally `collected_at`, `total_weight`,
  `facs_weight`, `vol_total`. Returns invisibly the number of rows inserted
  into `samples`. Writes to `sample_digestions` only for rows where at least
  one of `total_weight`/`facs_weight`/`vol_total` is non-`NA`.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-db_write_dimensions.R`:

```r
test_that("db_write_samples inserts samples and digestions where weights are given", {
  con <- local_test_db()
  db_write_experiment(con, experiment_code = "25-7")
  db_write_subjects(con, tibble::tibble(mouse_id = "25-7-1", experiment_code = "25-7"))

  n <- db_write_samples(con, tibble::tibble(
    mouse_id = c("25-7-1", "25-7-1"),
    tissue = c("spleen", "serum"),
    total_weight = c(0.12, NA),
    facs_weight = c(0.05, NA),
    vol_total = c(2.5, NA)
  ))
  expect_equal(n, 2)

  digestions <- DBI::dbGetQuery(con, "
    SELECT sm.tissue, d.total_weight
    FROM sample_digestions d JOIN samples sm ON sm.sample_id = d.sample_id
  ")
  expect_equal(digestions$tissue, "spleen")
  expect_equal(digestions$total_weight, 0.12)
})

test_that("db_write_samples is idempotent", {
  con <- local_test_db()
  db_write_experiment(con, experiment_code = "25-7")
  db_write_subjects(con, tibble::tibble(mouse_id = "25-7-1", experiment_code = "25-7"))

  db_write_samples(con, tibble::tibble(mouse_id = "25-7-1", tissue = "spleen"))
  n <- db_write_samples(con, tibble::tibble(mouse_id = "25-7-1", tissue = "spleen"))
  expect_equal(n, 0)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-db_write_dimensions.R')"`
Expected: FAIL -- `could not find function "db_write_samples"`

- [ ] **Step 3: Implement `db_write_samples()`**

Append to `R/db_write_dimensions.R`:

```r
#' Register samples (subject x tissue harvests)
#'
#' Bulk-inserts sample rows, and -- for tissue that went through FACS
#' dissociation -- the matching `sample_digestions` row. `data` must already
#' be shaped with one row per subject x tissue, referencing subjects already
#' registered via [db_write_subjects()].
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param data A tibble with columns `mouse_id`, `tissue`, and optionally
#'   `collected_at`, `total_weight`, `facs_weight`, `vol_total`. The three
#'   weight/volume columns are only written to `sample_digestions` for rows
#'   where at least one of them is non-`NA`.
#' @returns Invisibly, the number of rows inserted into `samples`.
#' @export
db_write_samples <- function(con, data) {
  data <- fill_missing_cols_(data, c(
    "mouse_id", "tissue", "collected_at", "total_weight", "facs_weight", "vol_total"
  ))
  duckdb::duckdb_register(con, "tmp_samples", data)
  on.exit(duckdb::duckdb_unregister(con, "tmp_samples"), add = TRUE)

  n <- DBI::dbExecute(con, "
    INSERT INTO samples (subject_id, tissue, collected_at)
    SELECT sub.subject_id, t.tissue, CAST(t.collected_at AS DATE)
    FROM tmp_samples t
    JOIN subjects sub ON sub.mouse_id = t.mouse_id
    ON CONFLICT (subject_id, tissue) DO NOTHING
  ")

  DBI::dbExecute(con, "
    INSERT INTO sample_digestions (sample_id, total_weight, facs_weight, vol_total)
    SELECT sm.sample_id, t.total_weight, t.facs_weight, t.vol_total
    FROM tmp_samples t
    JOIN subjects sub ON sub.mouse_id = t.mouse_id
    JOIN samples sm ON sm.subject_id = sub.subject_id AND sm.tissue = t.tissue
    WHERE t.total_weight IS NOT NULL OR t.facs_weight IS NOT NULL OR t.vol_total IS NOT NULL
    ON CONFLICT (sample_id) DO NOTHING
  ")

  invisible(n)
}
```

- [ ] **Step 4: Run `devtools::document()`**

Run: `Rscript -e "devtools::document()"`
Expected: `man/db_write_samples.Rd` created, no warnings.

- [ ] **Step 5: Run test to verify it passes**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-db_write_dimensions.R')"`
Expected: PASS, 0 failures.

- [ ] **Step 6: Commit**

```bash
git add R/db_write_dimensions.R tests/testthat/test-db_write_dimensions.R man/db_write_samples.Rd
git commit -m "feat: add db_write_samples()"
```

---

### Task 5: `db_write_assay()`

**Files:**
- Modify: `R/db_write_dimensions.R`
- Modify: `tests/testthat/test-db_write_dimensions.R`

**Interfaces:**
- Consumes: domains seeded by `db_connect()` (Task 1: `"facs"`, `"histo"`, `"elisa"`).
- Produces: `db_write_assay(con, assay_name, domain, description = NA_character_)`
  returns invisibly the number of rows inserted. Errors clearly if `domain`
  doesn't match a seeded row in `domains`.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-db_write_dimensions.R`:

```r
test_that("db_write_assay inserts an assay scoped to a domain", {
  con <- local_test_db()

  n <- db_write_assay(con, assay_name = "overview", domain = "facs")
  expect_equal(n, 1)

  row <- DBI::dbGetQuery(con, "
    SELECT a.assay_name, d.domain_name
    FROM assays a JOIN domains d ON d.domain_id = a.domain_id
    WHERE a.assay_name = 'overview'
  ")
  expect_equal(row$domain_name, "facs")
})

test_that("db_write_assay is idempotent", {
  con <- local_test_db()
  db_write_assay(con, assay_name = "overview", domain = "facs")
  n <- db_write_assay(con, assay_name = "overview", domain = "facs")
  expect_equal(n, 0)
})

test_that("db_write_assay errors clearly on an unknown domain", {
  con <- local_test_db()
  expect_error(
    db_write_assay(con, assay_name = "overview", domain = "not-a-domain"),
    "Unknown domain"
  )
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-db_write_dimensions.R')"`
Expected: FAIL -- `could not find function "db_write_assay"`

- [ ] **Step 3: Implement `db_write_assay()`**

Append to `R/db_write_dimensions.R`:

```r
#' Register an assay/panel
#'
#' Inserts a new assay identified by `assay_name` within `domain`. If an
#' assay with that name already exists in that domain, the call is a no-op.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param assay_name Character scalar, e.g. `"overview"` or `"anti-MPO"`.
#' @param domain Character scalar, one of `"facs"`, `"histo"`, `"elisa"` (or
#'   any additional domain later added to the `domains` table).
#' @param description Character scalar, optional.
#' @returns Invisibly, the number of rows inserted (`0` or `1`).
#' @export
db_write_assay <- function(con, assay_name, domain, description = NA_character_) {
  domain_row <- DBI::dbGetQuery(
    con, "SELECT domain_id FROM domains WHERE domain_name = ?",
    params = list(domain)
  )
  if (nrow(domain_row) == 0) {
    stop(
      "Unknown domain '", domain, "' -- expected one already present in the domains table.",
      call. = FALSE
    )
  }

  n <- DBI::dbExecute(
    con,
    "INSERT INTO assays (assay_name, domain_id, description)
     VALUES (?, ?, ?)
     ON CONFLICT (assay_name, domain_id) DO NOTHING",
    params = list(assay_name, domain_row$domain_id[[1]], description)
  )
  invisible(n)
}
```

- [ ] **Step 4: Run `devtools::document()`**

Run: `Rscript -e "devtools::document()"`
Expected: `man/db_write_assay.Rd` created, no warnings.

- [ ] **Step 5: Run test to verify it passes**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-db_write_dimensions.R')"`
Expected: PASS, 0 failures.

- [ ] **Step 6: Commit**

```bash
git add R/db_write_dimensions.R tests/testthat/test-db_write_dimensions.R man/db_write_assay.Rd
git commit -m "feat: add db_write_assay()"
```

---

### Task 6: `db_write_facs()`

**Files:**
- Create: `R/db_write_measurements.R`
- Create: `tests/testthat/test-db_write_measurements.R`
- Modify: `tests/testthat/helper-db.R` (add `seed_minimal_fixture()`)

**Interfaces:**
- Consumes: `db_write_experiment()`, `db_write_subjects()`, `db_write_samples()`, `db_write_assay()` from Tasks 2-5.
- Produces: internal helpers `lookup_sample_id_(con, mouse_id, tissue)` and
  `lookup_assay_id_(con, assay_name, domain)`, both returning a scalar integer
  ID or erroring with a clear message if no match exists.
- Produces: `db_write_facs(con, data, mouse_id, tissue, assay_name, count_method, vol_stained, vol_resuspended = NA_real_, vol_measured = NA_real_, bead_volume_added = NA_real_, bead_concentration = NA_real_, bead_population_path = NA_character_, stain_date = NA, source_file = NULL)`.
  `data` has columns `PopulationFullPath`, `Population`, `metric`, `value`
  (the shape `facs_import_wsp()` produces, filtered to one `FileName`).
  Returns invisibly the number of measurement rows inserted.

- [ ] **Step 1: Add the shared fixture helper**

Append to `tests/testthat/helper-db.R`:

```r
seed_minimal_fixture <- function(con) {
  db_write_experiment(con, experiment_code = "25-7")
  db_write_subjects(con, tibble::tibble(mouse_id = "25-7-1", experiment_code = "25-7"))
  db_write_samples(con, tibble::tibble(mouse_id = "25-7-1", tissue = "spleen"))
  db_write_assay(con, assay_name = "overview", domain = "facs")
  invisible(con)
}
```

- [ ] **Step 2: Write the failing test**

Create `tests/testthat/test-db_write_measurements.R`:

```r
test_that("db_write_facs creates a stain and inserts measurements (volumetric)", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  facs_data <- tibble::tibble(
    PopulationFullPath = c("root/Lymphocytes/CD4+", "root/Lymphocytes/CD8+"),
    Population = c("CD4+", "CD8+"),
    metric = c("Count", "Count"),
    value = c(1000, 800)
  )

  n <- db_write_facs(
    con, facs_data,
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    count_method = "volumetric", vol_stained = 0.1,
    vol_resuspended = 0.5, vol_measured = 0.05
  )
  expect_equal(n, 2)

  stains <- DBI::dbGetQuery(con, "SELECT * FROM facs_stains")
  expect_equal(nrow(stains), 1)
  expect_equal(stains$count_method, "volumetric")

  rows <- DBI::dbGetQuery(con, "SELECT population, value FROM facs_measurements ORDER BY population")
  expect_equal(rows$population, c("CD4+", "CD8+"))
  expect_equal(rows$value, c(1000, 800))
})

test_that("db_write_facs supports bead counting and rejects mixed fields", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  db_write_assay(con, assay_name = "treg", domain = "facs")

  facs_data <- tibble::tibble(
    PopulationFullPath = "root/Lymphocytes/CD4+",
    Population = "CD4+",
    metric = "Count",
    value = 1000
  )

  n <- db_write_facs(
    con, facs_data,
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    count_method = "bead", vol_stained = 0.1,
    bead_volume_added = 0.02, bead_concentration = 5000
  )
  expect_equal(n, 1)

  # Uses a different assay ("treg") than the call above, so this is a genuinely
  # new (sample_id, assay_id) pair -- the facs_stains INSERT is actually
  # attempted (not skipped by ON CONFLICT DO NOTHING), and the CHECK
  # constraint rejecting a mixed volumetric/bead field set fires as intended.
  expect_error(
    db_write_facs(
      con, facs_data,
      mouse_id = "25-7-1", tissue = "spleen", assay_name = "treg",
      count_method = "bead", vol_stained = 0.1,
      bead_volume_added = 0.02, bead_concentration = 5000,
      vol_measured = 0.05
    )
  )
})

test_that("db_write_facs errors clearly when the sample doesn't exist", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  facs_data <- tibble::tibble(
    PopulationFullPath = "root/Lymphocytes/CD4+", Population = "CD4+",
    metric = "Count", value = 1000
  )

  expect_error(
    db_write_facs(
      con, facs_data,
      mouse_id = "25-7-1", tissue = "kidney", assay_name = "overview",
      count_method = "volumetric", vol_stained = 0.1,
      vol_resuspended = 0.5, vol_measured = 0.05
    ),
    "No sample found"
  )
})

test_that("db_write_facs is idempotent", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  facs_data <- tibble::tibble(
    PopulationFullPath = "root/Lymphocytes/CD4+", Population = "CD4+",
    metric = "Count", value = 1000
  )
  args <- list(
    con = con, data = facs_data,
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    count_method = "volumetric", vol_stained = 0.1,
    vol_resuspended = 0.5, vol_measured = 0.05
  )
  do.call(db_write_facs, args)
  n <- do.call(db_write_facs, args)
  expect_equal(n, 0)
})
```

- [ ] **Step 3: Run test to verify it fails**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-db_write_measurements.R')"`
Expected: FAIL -- `could not find function "db_write_facs"`

- [ ] **Step 4: Implement shared lookup helpers and `db_write_facs()`**

Create `R/db_write_measurements.R`:

```r
lookup_sample_id_ <- function(con, mouse_id, tissue) {
  result <- DBI::dbGetQuery(
    con,
    "SELECT sm.sample_id
     FROM samples sm
     JOIN subjects sub ON sub.subject_id = sm.subject_id
     WHERE sub.mouse_id = ? AND sm.tissue = ?",
    params = list(mouse_id, tissue)
  )
  if (nrow(result) == 0) {
    stop(
      "No sample found for mouse_id '", mouse_id, "', tissue '", tissue,
      "' -- register the subject (db_write_subjects()) and sample (db_write_samples()) first.",
      call. = FALSE
    )
  }
  result$sample_id[[1]]
}

lookup_assay_id_ <- function(con, assay_name, domain) {
  result <- DBI::dbGetQuery(
    con,
    "SELECT a.assay_id
     FROM assays a
     JOIN domains d ON d.domain_id = a.domain_id
     WHERE a.assay_name = ? AND d.domain_name = ?",
    params = list(assay_name, domain)
  )
  if (nrow(result) == 0) {
    stop(
      "No assay found for assay_name '", assay_name, "' in domain '", domain,
      "' -- register it first via db_write_assay().",
      call. = FALSE
    )
  }
  result$assay_id[[1]]
}

#' Import FACS measurements for one stain
#'
#' Writes one staining event's tidy long measurements (as produced by
#' `facs_import_wsp()`, filtered to a single `FileName`) into
#' `facs_measurements`, creating the matching `facs_stains` row first if it
#' doesn't already exist. Which of `vol_resuspended`/`vol_measured` versus
#' `bead_volume_added`/`bead_concentration`/`bead_population_path` are
#' required depends on `count_method`; the database's `CHECK` constraint on
#' `facs_stains` rejects a call that supplies the wrong subset.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param data A tibble with columns `PopulationFullPath`, `Population`,
#'   `metric`, `value` -- the shape `facs_import_wsp()` produces, filtered to
#'   one `FileName`.
#' @param mouse_id Character scalar identifying the subject (already
#'   registered via [db_write_subjects()]).
#' @param tissue Character scalar identifying the sample (already registered
#'   via [db_write_samples()]).
#' @param assay_name Character scalar identifying the panel (already
#'   registered via [db_write_assay()] with `domain = "facs"`).
#' @param count_method Either `"volumetric"` or `"bead"`.
#' @param vol_stained Numeric scalar, volume of suspension stained.
#' @param vol_resuspended Numeric scalar, required when
#'   `count_method = "volumetric"`.
#' @param vol_measured Numeric scalar, required when
#'   `count_method = "volumetric"`.
#' @param bead_volume_added Numeric scalar, required when
#'   `count_method = "bead"`.
#' @param bead_concentration Numeric scalar, required when
#'   `count_method = "bead"`.
#' @param bead_population_path Character scalar, required when
#'   `count_method = "bead"`.
#' @param stain_date Date scalar, optional.
#' @param source_file Character scalar, optional -- defaults to
#'   `data$FileName[1]` when `data` has a `FileName` column.
#' @returns Invisibly, the number of measurement rows inserted.
#' @export
db_write_facs <- function(con, data, mouse_id, tissue, assay_name, count_method,
                           vol_stained,
                           vol_resuspended = NA_real_, vol_measured = NA_real_,
                           bead_volume_added = NA_real_, bead_concentration = NA_real_,
                           bead_population_path = NA_character_,
                           stain_date = NA, source_file = NULL) {
  count_method <- match.arg(count_method, c("volumetric", "bead"))

  sample_id <- lookup_sample_id_(con, mouse_id, tissue)
  assay_id <- lookup_assay_id_(con, assay_name, domain = "facs")

  DBI::dbExecute(
    con,
    "INSERT INTO facs_stains (
       sample_id, assay_id, vol_stained, count_method,
       vol_resuspended, vol_measured,
       bead_volume_added, bead_concentration, bead_population_path, stain_date
     ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
     ON CONFLICT (sample_id, assay_id) DO NOTHING",
    params = list(
      sample_id, assay_id, vol_stained, count_method,
      vol_resuspended, vol_measured,
      bead_volume_added, bead_concentration, bead_population_path, stain_date
    )
  )

  facs_stain_id <- DBI::dbGetQuery(
    con,
    "SELECT facs_stain_id FROM facs_stains WHERE sample_id = ? AND assay_id = ?",
    params = list(sample_id, assay_id)
  )$facs_stain_id[[1]]

  if (is.null(source_file)) {
    source_file <- if ("FileName" %in% names(data)) as.character(data$FileName[[1]]) else NA_character_
  }

  to_write <- tibble::tibble(
    facs_stain_id = facs_stain_id,
    population_full_path = as.character(data$PopulationFullPath),
    population = as.character(data$Population),
    metric = as.character(data$metric),
    value = as.double(data$value),
    source_file = source_file
  )
  duckdb::duckdb_register(con, "tmp_facs_measurements", to_write)
  on.exit(duckdb::duckdb_unregister(con, "tmp_facs_measurements"), add = TRUE)

  n <- DBI::dbExecute(con, "
    INSERT INTO facs_measurements (
      facs_stain_id, population_full_path, population, metric, value, source_file
    )
    SELECT facs_stain_id, population_full_path, population, metric, value, source_file
    FROM tmp_facs_measurements
    ON CONFLICT (facs_stain_id, population_full_path, metric) DO NOTHING
  ")

  invisible(n)
}
```

- [ ] **Step 5: Run `devtools::document()`**

Run: `Rscript -e "devtools::document()"`
Expected: `man/db_write_facs.Rd` created, no warnings.

- [ ] **Step 6: Run test to verify it passes**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-db_write_measurements.R')"`
Expected: PASS, 0 failures.

- [ ] **Step 7: Commit**

```bash
git add R/db_write_measurements.R tests/testthat/test-db_write_measurements.R tests/testthat/helper-db.R man/db_write_facs.Rd
git commit -m "feat: add db_write_facs()"
```

---

### Task 7: `db_write_elisa()`

**Files:**
- Modify: `R/db_write_measurements.R`
- Modify: `tests/testthat/test-db_write_measurements.R`
- Modify: `tests/testthat/helper-db.R` (seed an `"elisa"` assay too)

**Interfaces:**
- Consumes: `lookup_sample_id_()`, `lookup_assay_id_()` from Task 6.
- Produces: `db_write_elisa(con, data, mouse_id, tissue, assay_name, source_file = NA_character_)`
  where `data` has columns `cytokine`, `sample_id`, `replicate`, `value`,
  `unit`, `result_status` (the shape `elisa_read_results()` produces).
  Returns invisibly the number of rows inserted.

- [ ] **Step 1: Extend the shared fixture helper**

Update `seed_minimal_fixture()` in `tests/testthat/helper-db.R` to also seed
an ELISA assay:

```r
seed_minimal_fixture <- function(con) {
  db_write_experiment(con, experiment_code = "25-7")
  db_write_subjects(con, tibble::tibble(mouse_id = "25-7-1", experiment_code = "25-7"))
  db_write_samples(con, tibble::tibble(mouse_id = "25-7-1", tissue = "spleen"))
  db_write_assay(con, assay_name = "overview", domain = "facs")
  db_write_assay(con, assay_name = "anti-MPO", domain = "elisa")
  invisible(con)
}
```

- [ ] **Step 2: Write the failing test**

Append to `tests/testthat/test-db_write_measurements.R`:

```r
test_that("db_write_elisa inserts measurements", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  elisa_data <- tibble::tibble(
    cytokine = "MPO",
    sample_id = c("25-7-1", "25-7-1"),
    replicate = c(1L, 2L),
    value = c(12.3, 12.7),
    unit = "pg/ml",
    result_status = "OK"
  )

  n <- db_write_elisa(con, elisa_data, mouse_id = "25-7-1", tissue = "spleen", assay_name = "anti-MPO")
  expect_equal(n, 2)

  rows <- DBI::dbGetQuery(con, "SELECT replicate, value FROM elisa_measurements ORDER BY replicate")
  expect_equal(rows$replicate, c(1L, 2L))
  expect_equal(rows$value, c(12.3, 12.7))
})

test_that("db_write_elisa is idempotent", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  elisa_data <- tibble::tibble(
    cytokine = "MPO", sample_id = "25-7-1", replicate = 1L,
    value = 12.3, unit = "pg/ml", result_status = "OK"
  )
  args <- list(con = con, data = elisa_data, mouse_id = "25-7-1", tissue = "spleen", assay_name = "anti-MPO")
  do.call(db_write_elisa, args)
  n <- do.call(db_write_elisa, args)
  expect_equal(n, 0)
})
```

- [ ] **Step 3: Run test to verify it fails**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-db_write_measurements.R')"`
Expected: FAIL -- `could not find function "db_write_elisa"`

- [ ] **Step 4: Implement `db_write_elisa()`**

Append to `R/db_write_measurements.R`:

```r
#' Import ELISA measurements
#'
#' Writes a tidy tibble of ELISA measurements (as produced by
#' `elisa_read_results()`) into `elisa_measurements`.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param data A tibble with columns `cytokine`, `sample_id` (raw plate
#'   label), `replicate`, `value`, `unit`, `result_status` -- the shape
#'   `elisa_read_results()` produces.
#' @param mouse_id Character scalar identifying the subject.
#' @param tissue Character scalar identifying the sample.
#' @param assay_name Character scalar identifying the assay/panel (already
#'   registered via [db_write_assay()] with `domain = "elisa"`).
#' @param source_file Character scalar, optional.
#' @returns Invisibly, the number of rows inserted.
#' @export
db_write_elisa <- function(con, data, mouse_id, tissue, assay_name, source_file = NA_character_) {
  sample_id <- lookup_sample_id_(con, mouse_id, tissue)
  assay_id <- lookup_assay_id_(con, assay_name, domain = "elisa")

  to_write <- tibble::tibble(
    sample_id = sample_id,
    assay_id = assay_id,
    cytokine = as.character(data$cytokine),
    sample_id_raw = as.character(data$sample_id),
    replicate = as.integer(data$replicate),
    value = as.double(data$value),
    unit = as.character(data$unit),
    result_status = as.character(data$result_status),
    source_file = source_file
  )
  duckdb::duckdb_register(con, "tmp_elisa_measurements", to_write)
  on.exit(duckdb::duckdb_unregister(con, "tmp_elisa_measurements"), add = TRUE)

  n <- DBI::dbExecute(con, "
    INSERT INTO elisa_measurements (
      sample_id, assay_id, cytokine, sample_id_raw, replicate, value, unit, result_status, source_file
    )
    SELECT sample_id, assay_id, cytokine, sample_id_raw, replicate, value, unit, result_status, source_file
    FROM tmp_elisa_measurements
    ON CONFLICT (sample_id, assay_id, cytokine, sample_id_raw, replicate) DO NOTHING
  ")
  invisible(n)
}
```

- [ ] **Step 5: Run `devtools::document()`**

Run: `Rscript -e "devtools::document()"`
Expected: `man/db_write_elisa.Rd` created, no warnings.

- [ ] **Step 6: Run test to verify it passes**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-db_write_measurements.R')"`
Expected: PASS, 0 failures.

- [ ] **Step 7: Commit**

```bash
git add R/db_write_measurements.R tests/testthat/test-db_write_measurements.R tests/testthat/helper-db.R man/db_write_elisa.Rd
git commit -m "feat: add db_write_elisa()"
```

---

### Task 8: `db_write_histo()`

**Files:**
- Modify: `R/db_write_measurements.R`
- Modify: `tests/testthat/test-db_write_measurements.R`

**Interfaces:**
- Consumes: `lookup_sample_id_()`, `lookup_assay_id_()` from Task 6.
- Produces: `db_write_histo(con, data, mouse_id, tissue, assay_name = NA_character_, source_file = NA_character_)`
  where `data` has columns `metric`, `value`. Returns invisibly the number
  of rows inserted.

Note: because `assay_id` is nullable and SQL treats `NULL` as never equal to
another `NULL`, rows written with `assay_name = NA` are **not** deduplicated
against each other by the `UNIQUE(sample_id, assay_id, metric)` constraint --
this only matters once the histo domain tracks more than one un-named stain
type per sample/metric. This is documented on the function; no workaround is
implemented here (YAGNI until it's an actual problem).

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-db_write_measurements.R`:

```r
test_that("db_write_histo inserts measurements without an assay", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  histo_data <- tibble::tibble(metric = c("damage_score", "fibrosis_pct"), value = c(2, 15.5))

  n <- db_write_histo(con, histo_data, mouse_id = "25-7-1", tissue = "spleen")
  expect_equal(n, 2)

  rows <- DBI::dbGetQuery(con, "SELECT metric, value, assay_id FROM histo_measurements ORDER BY metric")
  expect_equal(rows$metric, c("damage_score", "fibrosis_pct"))
  expect_true(all(is.na(rows$assay_id)))
})

test_that("db_write_histo errors clearly on an unknown assay_name", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  histo_data <- tibble::tibble(metric = "damage_score", value = 2)
  expect_error(
    db_write_histo(con, histo_data, mouse_id = "25-7-1", tissue = "spleen", assay_name = "not-registered"),
    "No assay found"
  )
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-db_write_measurements.R')"`
Expected: FAIL -- `could not find function "db_write_histo"`

- [ ] **Step 3: Implement `db_write_histo()`**

Append to `R/db_write_measurements.R`:

```r
#' Import histology measurements
#'
#' Writes a tidy tibble of histology measurements into `histo_measurements`.
#'
#' Note: `assay_name` is optional (histology has traditionally used a single
#' un-named stain/scoring protocol). Rows written with `assay_name = NA` are
#' not deduplicated against each other -- SQL's `NULL <> NULL` means the
#' `UNIQUE(sample_id, assay_id, metric)` constraint doesn't catch a repeated
#' import when `assay_id` is `NULL`. Pass an `assay_name` (registered via
#' [db_write_assay()] with `domain = "histo"`) to get full dedup protection.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param data A tibble with columns `metric`, `value`.
#' @param mouse_id Character scalar identifying the subject.
#' @param tissue Character scalar identifying the sample.
#' @param assay_name Character scalar, optional -- identifies the
#'   stain/scoring protocol.
#' @param source_file Character scalar, optional.
#' @returns Invisibly, the number of rows inserted.
#' @export
db_write_histo <- function(con, data, mouse_id, tissue, assay_name = NA_character_,
                            source_file = NA_character_) {
  sample_id <- lookup_sample_id_(con, mouse_id, tissue)
  assay_id <- if (!is.na(assay_name)) lookup_assay_id_(con, assay_name, domain = "histo") else NA_integer_

  to_write <- tibble::tibble(
    sample_id = sample_id,
    assay_id = assay_id,
    metric = as.character(data$metric),
    value = as.double(data$value),
    source_file = source_file
  )
  duckdb::duckdb_register(con, "tmp_histo_measurements", to_write)
  on.exit(duckdb::duckdb_unregister(con, "tmp_histo_measurements"), add = TRUE)

  n <- DBI::dbExecute(con, "
    INSERT INTO histo_measurements (sample_id, assay_id, metric, value, source_file)
    SELECT sample_id, assay_id, metric, value, source_file
    FROM tmp_histo_measurements
    ON CONFLICT (sample_id, assay_id, metric) DO NOTHING
  ")
  invisible(n)
}
```

- [ ] **Step 4: Run `devtools::document()`**

Run: `Rscript -e "devtools::document()"`
Expected: `man/db_write_histo.Rd` created, no warnings.

- [ ] **Step 5: Run test to verify it passes**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-db_write_measurements.R')"`
Expected: PASS, 0 failures.

- [ ] **Step 6: Commit**

```bash
git add R/db_write_measurements.R tests/testthat/test-db_write_measurements.R man/db_write_histo.Rd
git commit -m "feat: add db_write_histo()"
```

---

### Task 9: Query helpers

**Files:**
- Create: `R/db_query.R`
- Create: `tests/testthat/test-db_query.R`

**Interfaces:**
- Consumes: all write functions from Tasks 2-8 (used to seed test data).
- Produces: `db_query_facs(con)`, `db_query_elisa(con)`, `db_query_histo(con)`,
  each returning a lazy `dbplyr` `tbl` exposing natural keys
  (`experiment_code`, `mouse_id`, `tissue`, `assay_name`) instead of
  surrogate IDs. Callers `dplyr::collect()` as needed.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-db_query.R`:

```r
test_that("db_query_facs returns a flat tibble keyed by natural keys", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  db_write_facs(
    con,
    tibble::tibble(
      PopulationFullPath = "root/Lymphocytes/CD4+", Population = "CD4+",
      metric = "Count", value = 1000
    ),
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    count_method = "volumetric", vol_stained = 0.1,
    vol_resuspended = 0.5, vol_measured = 0.05
  )

  result <- db_query_facs(con) |> dplyr::collect()

  expect_equal(result$mouse_id, "25-7-1")
  expect_equal(result$experiment_code, "25-7")
  expect_equal(result$tissue, "spleen")
  expect_equal(result$assay_name, "overview")
  expect_equal(result$population, "CD4+")
  expect_equal(result$value, 1000)
})

test_that("db_query_elisa returns a flat tibble keyed by natural keys", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  db_write_elisa(
    con,
    tibble::tibble(
      cytokine = "MPO", sample_id = "25-7-1", replicate = 1L,
      value = 12.3, unit = "pg/ml", result_status = "OK"
    ),
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "anti-MPO"
  )

  result <- db_query_elisa(con) |> dplyr::collect()

  expect_equal(result$mouse_id, "25-7-1")
  expect_equal(result$assay_name, "anti-MPO")
  expect_equal(result$cytokine, "MPO")
  expect_equal(result$value, 12.3)
})

test_that("db_query_histo returns a flat tibble and tolerates a NULL assay_id", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  db_write_histo(
    con, tibble::tibble(metric = "damage_score", value = 2),
    mouse_id = "25-7-1", tissue = "spleen"
  )

  result <- db_query_histo(con) |> dplyr::collect()

  expect_equal(result$mouse_id, "25-7-1")
  expect_equal(result$metric, "damage_score")
  expect_true(is.na(result$assay_name))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-db_query.R')"`
Expected: FAIL -- `could not find function "db_query_facs"`

- [ ] **Step 3: Implement the query helpers**

Create `R/db_query.R`:

```r
#' Query FACS measurements
#'
#' Returns a lazy `dbplyr` `tbl` joining `facs_measurements` back through
#' `facs_stains`/`samples`/`subjects`/`experiments`/`assays`, exposing
#' natural keys (`mouse_id`, `tissue`, `experiment_code`, `assay_name`)
#' instead of surrogate IDs. Call [dplyr::collect()] to materialize.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @returns A lazy `tbl` -- filter/collect as needed.
#' @export
db_query_facs <- function(con) {
  dplyr::tbl(con, "facs_measurements") |>
    dplyr::left_join(dplyr::tbl(con, "facs_stains"), by = "facs_stain_id") |>
    dplyr::left_join(dplyr::tbl(con, "samples"), by = "sample_id") |>
    dplyr::left_join(dplyr::tbl(con, "subjects"), by = "subject_id") |>
    dplyr::left_join(dplyr::tbl(con, "experiments"), by = "experiment_id") |>
    dplyr::left_join(dplyr::tbl(con, "assays"), by = "assay_id") |>
    dplyr::select(
      experiment_code, mouse_id, tissue, assay_name,
      population_full_path, population, metric, value,
      count_method, vol_stained, vol_resuspended, vol_measured,
      bead_volume_added, bead_concentration, bead_population_path,
      source_file, imported_at
    )
}

#' Query ELISA measurements
#'
#' Returns a lazy `dbplyr` `tbl` joining `elisa_measurements` back through
#' `samples`/`subjects`/`experiments`/`assays`, exposing natural keys instead
#' of surrogate IDs. Call [dplyr::collect()] to materialize.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @returns A lazy `tbl` -- filter/collect as needed.
#' @export
db_query_elisa <- function(con) {
  dplyr::tbl(con, "elisa_measurements") |>
    dplyr::left_join(dplyr::tbl(con, "samples"), by = "sample_id") |>
    dplyr::left_join(dplyr::tbl(con, "subjects"), by = "subject_id") |>
    dplyr::left_join(dplyr::tbl(con, "experiments"), by = "experiment_id") |>
    dplyr::left_join(dplyr::tbl(con, "assays"), by = "assay_id") |>
    dplyr::select(
      experiment_code, mouse_id, tissue, assay_name,
      cytokine, sample_id_raw, replicate, value, unit, result_status,
      source_file, imported_at
    )
}

#' Query histology measurements
#'
#' Returns a lazy `dbplyr` `tbl` joining `histo_measurements` back through
#' `samples`/`subjects`/`experiments`/`assays`, exposing natural keys instead
#' of surrogate IDs. `assay_name` is `NA` for rows written without an
#' `assay_name`. Call [dplyr::collect()] to materialize.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @returns A lazy `tbl` -- filter/collect as needed.
#' @export
db_query_histo <- function(con) {
  dplyr::tbl(con, "histo_measurements") |>
    dplyr::left_join(dplyr::tbl(con, "samples"), by = "sample_id") |>
    dplyr::left_join(dplyr::tbl(con, "subjects"), by = "subject_id") |>
    dplyr::left_join(dplyr::tbl(con, "experiments"), by = "experiment_id") |>
    dplyr::left_join(dplyr::tbl(con, "assays"), by = "assay_id") |>
    dplyr::select(
      experiment_code, mouse_id, tissue, assay_name,
      metric, value, source_file, imported_at
    )
}
```

- [ ] **Step 4: Run `devtools::document()`**

Run: `Rscript -e "devtools::document()"`
Expected: `man/db_query_facs.Rd`, `man/db_query_elisa.Rd`, `man/db_query_histo.Rd` created, no warnings.

- [ ] **Step 5: Run test to verify it passes**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-db_query.R')"`
Expected: PASS, 0 failures.

- [ ] **Step 6: Commit**

```bash
git add R/db_query.R tests/testthat/test-db_query.R man/db_query_facs.Rd man/db_query_elisa.Rd man/db_query_histo.Rd
git commit -m "feat: add db_query_facs/elisa/histo()"
```

---

### Task 10: Full check and documentation update

**Files:**
- Modify: `CLAUDE.md`

**Interfaces:**
- Consumes: nothing new -- this task verifies and documents the completed feature.

- [ ] **Step 1: Run the full test suite**

Run: `Rscript -e "devtools::load_all(); devtools::test()"`
Expected: PASS, 0 failures across all test files (existing domains plus the four new `test-db_*.R` files).

- [ ] **Step 2: Run full `R CMD check`**

Run: `Rscript -e "devtools::check()"`
Expected: 0 errors, 0 warnings. Record any new NOTEs produced by the `db_`
files (expected: bare column-name NOTEs inside `dplyr` verbs in
`R/db_query.R`, following the same false-positive pattern as every other
domain).

- [ ] **Step 3: Update `CLAUDE.md`'s domain table**

In the "File and function naming" section's domain list, change:

```
- `db_` — database access [stub — no functions yet]
```

to:

```
- `db_` — database access (DuckDB: connect/write/query the FACS/ELISA/histology schema)
```

- [ ] **Step 4: Add a `db_` data-structure entry**

In the "Data structures" section of `CLAUDE.md`, add a new subsection
(matching the style of the existing `facs_import_wsp()`/`meta_read()`
subsections) describing the schema at a level of detail useful to a future
session -- reference `docs/superpowers/specs/2026-07-17-facs-db-schema-design.md`
for the full column-level DDL rather than duplicating it:

```markdown
### `db_connect()` / `db_write_*()` / `db_query_*()` — lab measurements database
- **Schema**: see `docs/superpowers/specs/2026-07-17-facs-db-schema-design.md`
  for full DDL. Two layers: dimensions (`domains`, `experiments`, `subjects`,
  `samples`, `sample_digestions`, `assays`) and per-domain fact tables
  (`facs_stains` + `facs_measurements`, `elisa_measurements`,
  `histo_measurements`), on DuckDB.
- **`db_connect(path)`**: opens/creates the database file and idempotently
  ensures the full schema exists, seeding `domains` with `"facs"`/`"histo"`/`"elisa"`.
- **`db_write_experiment()`/`db_write_subjects()`/`db_write_samples()`/`db_write_assay()`**:
  register dimension rows. All idempotent (`ON CONFLICT DO NOTHING` against
  the schema's `UNIQUE` constraints) -- safe to re-run.
- **`db_write_facs()`/`db_write_elisa()`/`db_write_histo()`**: import a
  domain's tidy long tibble (matching that domain's existing read function's
  output shape) against an already-registered `mouse_id`/`tissue`/`assay_name`.
  All idempotent. `db_write_facs()` additionally creates the `facs_stains` row
  (one per sample x panel) on first write, branching on `count_method`
  (`"volumetric"` or `"bead"`) -- the database's `CHECK` constraint enforces
  that only the fields matching the chosen method are populated.
- **`db_query_facs()`/`db_query_elisa()`/`db_query_histo()`**: return a lazy
  `dbplyr` `tbl` joining a fact table back to natural keys (`experiment_code`,
  `mouse_id`, `tissue`, `assay_name`) -- callers never see surrogate IDs.
  `dplyr::collect()` to materialize.
- Migrating existing `data/MPOmRNA.db` (SQLite) data into this schema is a
  separate, not-yet-started task.
```

- [ ] **Step 5: Update the "Known check output" section**

Append a bullet documenting the NOTEs actually observed in Step 2 (fill in
the exact variable names once Step 2 has run -- these are expected to be the
bare column names selected in `R/db_query.R`'s `dplyr::select()` calls,
e.g. `experiment_code`, `mouse_id`, `tissue`, `assay_name`, `population`,
`metric`, `value`, `cytokine`, `sample_id_raw`, `replicate`, `unit`,
`result_status`, `source_file`, `imported_at`, following the same
false-positive pattern as every other domain in this section).

- [ ] **Step 6: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: document the db_ domain in CLAUDE.md"
```
