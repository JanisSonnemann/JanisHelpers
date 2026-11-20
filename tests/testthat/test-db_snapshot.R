test_that("db_snapshot exports a restorable dated snapshot and keeps con open", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  root <- withr::local_tempdir()

  out <- suppressMessages(db_snapshot(con, root))

  expect_equal(
    normalizePath(out),
    normalizePath(file.path(root, format(Sys.Date(), "%Y-%m-%d")))
  )
  expect_true(DBI::dbIsValid(con))

  restored <- DBI::dbConnect(duckdb::duckdb(), tempfile(fileext = ".duckdb"))
  withr::defer(DBI::dbDisconnect(restored, shutdown = TRUE))
  DBI::dbExecute(restored, sprintf("IMPORT DATABASE '%s'", out))
  expect_equal(
    DBI::dbGetQuery(restored, "SELECT experiment_code FROM experiments")$experiment_code,
    "25-7"
  )
})

test_that("db_snapshot replaces a same-day snapshot", {
  con <- local_test_db()
  root <- withr::local_tempdir()
  stale <- file.path(root, format(Sys.Date(), "%Y-%m-%d"), "stale.txt")
  dir.create(dirname(stale), recursive = TRUE)
  file.create(stale)

  suppressMessages(db_snapshot(con, root))

  expect_false(file.exists(stale))
})

test_that("snapshot works with exclusion tables", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  DBI::dbExecute(con, "INSERT INTO exclusion_sets
    (experiment_id, package_version)
    SELECT experiment_id, '0.0.0' FROM experiments")
  DBI::dbExecute(con, "INSERT INTO exclusion_rules
    (exclusion_set_id, mouse_id, reason)
    SELECT exclusion_set_id, '25-7-1', 'test' FROM exclusion_sets")
  root <- withr::local_tempdir()

  out <- suppressMessages(db_snapshot(con, root))

  expect_true(any(grepl("exclusion_rules", list.files(out))))
})

test_that("snapshot with flag views restores and reopens", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  db_write_facs(
    con,
    tibble::tibble(PopulationFullPath = "root/CD4+", Population = "CD4+",
                   metric = "count", value = 1),
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    wsp_root = "overview_25-7", source_file = "25-7/overview_25-7_20261101.wsp"
  )
  suppressMessages(db_write_exclusions(
    con, tibble::tibble(mouse_ID = "25-7-1", reason = "test"), "25-7"))
  root <- withr::local_tempdir()

  out <- suppressMessages(db_snapshot(con, root))

  path <- tempfile(fileext = ".duckdb")
  restored <- DBI::dbConnect(duckdb::duckdb(), path)
  DBI::dbExecute(restored, sprintf("IMPORT DATABASE '%s'", out))
  expect_equal(
    DBI::dbGetQuery(restored, "SELECT exclusion_reason FROM facs_exclusion_flags")$exclusion_reason,
    "test")
  DBI::dbDisconnect(restored, shutdown = TRUE)
  reopened <- db_connect(path)
  withr::defer(DBI::dbDisconnect(reopened, shutdown = TRUE))
  expect_true(dplyr::collect(db_query_facs(reopened))$excluded)
})
