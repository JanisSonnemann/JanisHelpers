test_that("db_connect creates all tables and seeds domains", {
  con <- local_test_db()

  expect_true(all(c(
    "domains", "experiments", "subjects", "samples", "import_log",
    "assays", "facs_stains", "facs_measurements", "elisa_measurements",
    "histo_measurements", "import_exclusions"
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

test_that("schema has wsp_root, no calc-input columns, no sample_digestions", {
  con <- local_test_db()
  cols <- DBI::dbListFields(con, "facs_stains")
  expect_true(all(c("wsp_root", "condition", "condition_keyword") %in% cols))
  expect_false(any(c(
    "vol_stained", "count_method", "vol_resuspended", "vol_measured",
    "bead_volume_added", "bead_concentration", "bead_population_path"
  ) %in% cols))
  expect_false("sample_digestions" %in% DBI::dbListTables(con))
  expect_setequal(DBI::dbListFields(con, "import_log"), c(
    "import_id", "domain_id", "experiment_code", "import_unit", "source_file",
    "source_version", "qc_report", "qc_passed", "rows_written",
    "package_version", "supersedes_import_id", "imported_at"
  ))
})

test_that("db_connect errors on a file created by the old schema", {
  path <- tempfile(fileext = ".duckdb")
  old <- DBI::dbConnect(duckdb::duckdb(), path)
  DBI::dbExecute(old, "CREATE TABLE facs_stains (facs_stain_id INTEGER, count_method VARCHAR)")
  DBI::dbDisconnect(old, shutdown = TRUE)
  expect_error(db_connect(path), "older schema")
})

test_that("db_connect errors on a file whose facs_stains has no condition column", {
  path <- tempfile(fileext = ".duckdb")
  old <- DBI::dbConnect(duckdb::duckdb(), path)
  DBI::dbExecute(old, "CREATE TABLE facs_stains (facs_stain_id INTEGER, wsp_root VARCHAR)")
  DBI::dbDisconnect(old, shutdown = TRUE)
  expect_error(db_connect(path), "older schema")
})

test_that("a database with data can be snapshotted with EXPORT DATABASE", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  import_id <- JanisHelpers:::write_import_log_(
    con, "facs", "25-7", "overview_25-7",
    "25-7/facs/overview_25-7_20261112.wsp", rows_written = 0L
  )
  DBI::dbExecute(con, "
    INSERT INTO import_exclusions (import_id, mouse_id, reason, n_files_matched)
    VALUES (?, '25-7-1', 'clots', 1)", params = list(import_id))
  out <- tempfile("snapshot")
  withr::defer(unlink(out, recursive = TRUE))

  DBI::dbExecute(con, sprintf("EXPORT DATABASE '%s' (FORMAT PARQUET)", out))

  expect_true(file.exists(file.path(out, "schema.sql")))
  expect_true(file.exists(file.path(out, "import_log.parquet")))
  expect_true(file.exists(file.path(out, "import_exclusions.parquet")))
})

test_that("a root's measurements and stains can be replaced in one transaction", {
  # Mirrors what the Stage 3 export will do in ONE transaction (DuckDB does
  # not nest transactions, so it can't wrap db_write_facs() calls): delete
  # the root's measurements, then its stains, then insert the new version.
  con <- local_test_db()
  seed_minimal_fixture(con)
  db_write_facs(
    con,
    tibble::tibble(
      PopulationFullPath = "root/CD4+", Population = "CD4+",
      metric = "Count", value = 1
    ),
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    wsp_root = "overview_25-7", source_file = "25-7/facs/overview_25-7_20261112.wsp"
  )
  expect_no_error(DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "DELETE FROM facs_measurements WHERE facs_stain_id IN
      (SELECT facs_stain_id FROM facs_stains WHERE wsp_root = 'overview_25-7')")
    DBI::dbExecute(con, "DELETE FROM facs_stains WHERE wsp_root = 'overview_25-7'")
    DBI::dbExecute(con, "
      INSERT INTO facs_stains (sample_id, assay_id, wsp_root)
      SELECT sm.sample_id, a.assay_id, 'overview_25-7'
      FROM samples sm, assays a WHERE a.assay_name = 'overview'")
    DBI::dbExecute(con, "
      INSERT INTO facs_measurements
        (facs_stain_id, population_full_path, population, metric, value, source_file)
      SELECT facs_stain_id, 'root/CD4+', 'CD4+', 'Count', 2,
        '25-7/facs/overview_25-7_20261120.wsp'
      FROM facs_stains WHERE wsp_root = 'overview_25-7'")
  }))
  expect_equal(
    DBI::dbGetQuery(con, "SELECT source_file FROM facs_measurements")$source_file,
    "25-7/facs/overview_25-7_20261120.wsp"
  )
})

test_that("db_connect errors on a samples table with nullable collected_at", {
  path <- tempfile(fileext = ".duckdb")
  old <- DBI::dbConnect(duckdb::duckdb(), path)
  DBI::dbExecute(old, "CREATE TABLE samples (sample_id INTEGER, subject_id INTEGER,
    tissue VARCHAR, collected_at DATE, UNIQUE(subject_id, tissue))")
  DBI::dbDisconnect(old, shutdown = TRUE)
  err <- expect_error(db_connect(path), "older schema")
  expect_match(conditionMessage(err), "Recreate")
})

test_that("import_exclusions has collected_at and the append-only writers are gone", {
  con <- local_test_db()
  expect_true("collected_at" %in% DBI::dbListFields(con, "import_exclusions"))
  expect_false(exists("db_write_elisa", asNamespace("JanisHelpers")))
  expect_false(exists("db_write_histo", asNamespace("JanisHelpers")))
})

test_that("db_connect errors on an elisa_measurements table without analyte", {
  path <- tempfile(fileext = ".duckdb")
  old <- DBI::dbConnect(duckdb::duckdb(), path)
  DBI::dbExecute(old, "CREATE TABLE elisa_measurements (elisa_measurement_id INTEGER, cytokine VARCHAR)")
  DBI::dbDisconnect(old, shutdown = TRUE)
  err <- expect_error(db_connect(path), "older schema")
  expect_match(conditionMessage(err), "analyte")
})

test_that("facs_measurements has a nullable population_order column", {
  con <- local_test_db()
  expect_true("population_order" %in% DBI::dbListFields(con, "facs_measurements"))
})

test_that("db_connect adds population_order to an existing database, keeping rows", {
  path <- tempfile(fileext = ".duckdb")
  con <- db_connect(path)
  if ("population_order" %in% DBI::dbListFields(con, "facs_measurements")) {
    DBI::dbExecute(con, "ALTER TABLE facs_measurements DROP COLUMN population_order")
  }
  DBI::dbExecute(con, "INSERT INTO facs_measurements
    (facs_stain_id, population_full_path, population, metric, value)
    VALUES (1, 'root/CD4+', 'CD4+', 'count', 10)")
  DBI::dbDisconnect(con, shutdown = TRUE)

  con <- db_connect(path)
  withr::defer(DBI::dbDisconnect(con, shutdown = TRUE))
  rows <- DBI::dbGetQuery(con, "SELECT * FROM facs_measurements")
  expect_equal(rows$population, "CD4+")
  expect_true("population_order" %in% names(rows))
  expect_true(is.na(rows$population_order))
})

test_that("new schema has exclusion tables and subject flags", {
  con <- local_test_db()
  expect_true(DBI::dbExistsTable(con, "exclusion_sets"))
  expect_true(DBI::dbExistsTable(con, "exclusion_rules"))
  expect_true(all(
    c("early", "early_reason", "id_source") %in% DBI::dbListFields(con, "subjects")
  ))
  expect_true("file_name" %in% DBI::dbListFields(con, "facs_stains"))
})

test_that("files from the previous schema are migrated", {
  path <- tempfile(fileext = ".duckdb")
  # DuckDB refuses DROP COLUMN on tables other tables depend on, so build the
  # previous-schema tables by hand (pre-change DDL, no dependents).
  old <- DBI::dbConnect(duckdb::duckdb(), path)
  DBI::dbExecute(old, "CREATE SEQUENCE seq_experiments START 1")
  DBI::dbExecute(old, "CREATE TABLE experiments (
    experiment_id INTEGER PRIMARY KEY DEFAULT nextval('seq_experiments'),
    experiment_code VARCHAR NOT NULL UNIQUE, experiment_name VARCHAR,
    project VARCHAR, description VARCHAR,
    created_at TIMESTAMP DEFAULT current_timestamp)")
  DBI::dbExecute(old, "CREATE SEQUENCE seq_subjects START 1")
  DBI::dbExecute(old, "CREATE TABLE subjects (
    subject_id INTEGER PRIMARY KEY DEFAULT nextval('seq_subjects'),
    mouse_id VARCHAR NOT NULL UNIQUE,
    experiment_id INTEGER NOT NULL REFERENCES experiments(experiment_id),
    mdc_id VARCHAR, cage VARCHAR, mouse_strain VARCHAR, generation VARCHAR,
    sex VARCHAR, mouse_treatment VARCHAR, treatment_group VARCHAR,
    dob DATE, start_date DATE, bmt_date DATE, end_date DATE,
    created_at TIMESTAMP DEFAULT current_timestamp)")
  DBI::dbExecute(old, "CREATE TABLE facs_stains (
    facs_stain_id INTEGER PRIMARY KEY, sample_id INTEGER NOT NULL,
    assay_id INTEGER NOT NULL, wsp_root VARCHAR NOT NULL,
    condition VARCHAR NOT NULL DEFAULT 'none', condition_keyword VARCHAR,
    stain_date DATE)")
  DBI::dbExecute(old, "INSERT INTO experiments (experiment_code) VALUES ('25-7')")
  DBI::dbExecute(old, "INSERT INTO subjects (mouse_id, experiment_id)
    SELECT '25-7-1', experiment_id FROM experiments")
  DBI::dbDisconnect(old, shutdown = TRUE)

  con <- db_connect(path)
  withr::defer(DBI::dbDisconnect(con, shutdown = TRUE))
  expect_true(all(
    c("early", "early_reason", "id_source") %in% DBI::dbListFields(con, "subjects")
  ))
  expect_true("file_name" %in% DBI::dbListFields(con, "facs_stains"))
  expect_true(DBI::dbExistsTable(con, "exclusion_sets"))
  expect_true(DBI::dbExistsTable(con, "exclusion_rules"))
  s <- DBI::dbGetQuery(con, "SELECT early, id_source FROM subjects")
  expect_false(s$early)
  expect_equal(s$id_source, "mouse_ID")
})
