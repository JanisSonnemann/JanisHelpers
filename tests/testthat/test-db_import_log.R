test_that("write_import_log_ writes one row with package version and returns its id", {
  con <- local_test_db()

  id <- write_import_log_(
    con,
    domain = "facs",
    experiment_code = "25-7",
    import_unit = "overview_25-7",
    source_file = "25-7/facs/overview_25-7_20261112.wsp",
    rows_written = 5L,
    source_version = as.Date("2026-11-12"),
    qc_report = "25-7/facs/qc/2026-11-12_qc.html",
    qc_passed = TRUE
  )

  log <- DBI::dbGetQuery(con, "
    SELECT l.*, d.domain_name FROM import_log l
    JOIN domains d ON d.domain_id = l.domain_id
  ")
  expect_equal(nrow(log), 1)
  expect_equal(id, log$import_id)
  expect_equal(log$domain_name, "facs")
  expect_equal(log$experiment_code, "25-7")
  expect_equal(log$import_unit, "overview_25-7")
  expect_equal(log$source_version, as.Date("2026-11-12"))
  expect_equal(log$qc_report, "25-7/facs/qc/2026-11-12_qc.html")
  expect_true(log$qc_passed)
  expect_equal(log$rows_written, 5L)
  expect_equal(log$package_version, as.character(utils::packageVersion("JanisHelpers")))
  expect_false(is.na(log$imported_at))
  expect_true(is.na(log$supersedes_import_id))
})

test_that("write_import_log_ links supersedes_import_id", {
  con <- local_test_db()
  args <- list(
    con = con, domain = "facs", experiment_code = "25-7",
    import_unit = "overview_25-7", rows_written = 1L
  )

  first <- do.call(write_import_log_, c(args, source_file = "25-7/facs/overview_25-7_20261112.wsp"))
  second <- do.call(write_import_log_, c(args,
    source_file = "25-7/facs/overview_25-7_20261120.wsp",
    supersedes_import_id = first
  ))

  linked <- DBI::dbGetQuery(
    con, "SELECT supersedes_import_id FROM import_log WHERE import_id = ?",
    params = list(second)
  )$supersedes_import_id
  expect_equal(linked, first)
})
