test_that("db_connect creates all tables and seeds domains", {
  con <- local_test_db()

  expect_true(all(c(
    "domains", "experiments", "subjects", "samples", "import_log",
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

test_that("schema has wsp_root, no calc-input columns, no sample_digestions", {
  con <- local_test_db()
  cols <- DBI::dbListFields(con, "facs_stains")
  expect_true("wsp_root" %in% cols)
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

test_that("a database with data can be snapshotted with EXPORT DATABASE", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  out <- tempfile("snapshot")
  withr::defer(unlink(out, recursive = TRUE))

  DBI::dbExecute(con, sprintf("EXPORT DATABASE '%s' (FORMAT PARQUET)", out))

  expect_true(file.exists(file.path(out, "schema.sql")))
  expect_true(file.exists(file.path(out, "import_log.parquet")))
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
