two_mice <- tibble::tibble(
  file_name = c("a.fcs", "b.fcs"),
  mouse_ID = c("25-7-1", "25-7-2"),
  tissue = "spleen"
)

seed_export_fixture <- function(con) {
  db_write_experiment(con, experiment_code = "25-7")
  db_write_subjects(con, tibble::tibble(
    mouse_id = c("25-7-1", "25-7-2"), experiment_code = "25-7"
  ))
  invisible(con)
}

n_rows <- function(con, table) {
  DBI::dbGetQuery(con, paste("SELECT COUNT(*) AS n FROM", table))$n
}

# A second, later workspace version of the same root with its own prep.
second_version <- function(d, value = 200) {
  wsp2 <- file.path(d$facs_dir, "overview_25-7_20261105.wsp")
  file.create(wsp2)
  Sys.setFileTime(wsp2, Sys.time() - 1800)
  suppressMessages(local_loaded_prep(d, two_mice, value = value, wsp_path = wsp2))
}

test_that("exports all files, logs the import with relative paths", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep <- suppressMessages(local_loaded_prep(d, two_mice))

  id <- suppressMessages(db_export_facs(con, prep, data_dir = d$data_dir))

  log <- DBI::dbGetQuery(con, "SELECT * FROM import_log WHERE import_id = ?",
                         params = list(id))
  expect_equal(log$import_unit, "overview_25-7")
  expect_equal(log$source_file, "25-7_exp/raw/FACS/overview_25-7_20261101.wsp")
  expect_equal(
    log$qc_report,
    "25-7_exp/raw/FACS/overview/qc/overview_25-7_20261101_qc.html"
  )
  expect_true(log$qc_passed)
  expect_equal(as.Date(log$source_version), as.Date("2026-11-01"))
  expect_equal(log$experiment_code, "25-7")
  expect_equal(n_rows(con, "facs_stains"), 2)
  expect_equal(log$rows_written, n_rows(con, "facs_measurements"))
  expect_equal(log$rows_written, 8)
})

test_that("applies exclusions and records them with the import", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep <- suppressMessages(local_loaded_prep(d, two_mice))
  ex <- tibble::tibble(mouse_ID = "25-7-2", reason = "clots",
                       date = as.Date("2026-11-10"))

  id <- suppressMessages(db_export_facs(con, prep, ex, data_dir = d$data_dir))

  expect_equal(n_rows(con, "facs_stains"), 1)
  mice <- dplyr::collect(db_query_facs(con))$mouse_id
  expect_equal(unique(mice), "25-7-1")
  exc <- DBI::dbGetQuery(con, "SELECT * FROM import_exclusions")
  expect_equal(nrow(exc), 1)
  expect_equal(exc$import_id, id)
  expect_equal(exc$mouse_id, "25-7-2")
  expect_equal(exc$reason, "clots")
  expect_equal(exc$n_files_matched, 1)
  expect_equal(as.Date(exc$decided_on), as.Date("2026-11-10"))
  expect_true(is.na(exc$tissue) && is.na(exc$panel) && is.na(exc$file_name))
})

test_that("unmatched exclusion row stops the export, nothing written", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep <- suppressMessages(local_loaded_prep(d, two_mice))
  ex <- tibble::tibble(mouse_ID = "25-7-1", tissue = "Spleen", reason = "typo")

  expect_error(db_export_facs(con, prep, ex, data_dir = d$data_dir),
               "match no file")
  expect_equal(n_rows(con, "facs_stains"), 0)
  expect_equal(n_rows(con, "import_log"), 0)
})

test_that("all files excluded stops the export", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep <- suppressMessages(local_loaded_prep(d, two_mice))
  ex <- tibble::tibble(mouse_ID = c("25-7-1", "25-7-2"), reason = "x")

  expect_error(db_export_facs(con, prep, ex, data_dir = d$data_dir),
               "nothing left")
  expect_equal(n_rows(con, "import_log"), 0)
})

test_that("keyword clash stops the export", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  clash <- dplyr::mutate(two_mice, mouse_ID = "25-7-1")
  prep <- suppressWarnings(suppressMessages(local_loaded_prep(d, clash)))

  expect_error(db_export_facs(con, prep, data_dir = d$data_dir), "share")
  expect_equal(n_rows(con, "facs_stains"), 0)
})

test_that("a prep not loaded with facs_prep_load() is refused", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep <- suppressMessages(local_loaded_prep(d, two_mice))
  prep$paths <- NULL
  expect_error(db_export_facs(con, prep, data_dir = d$data_dir), "facs_prep_load")
})

test_that("a workspace outside data_dir is refused", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep <- suppressMessages(local_loaded_prep(d, two_mice))
  expect_error(
    db_export_facs(con, prep, data_dir = withr::local_tempdir()),
    "data_dir"
  )
  expect_equal(n_rows(con, "import_log"), 0)
})

test_that("existing root: overwrite = FALSE aborts without changes", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep1 <- suppressMessages(local_loaded_prep(d, two_mice))
  suppressMessages(db_export_facs(con, prep1, data_dir = d$data_dir))
  prep2 <- second_version(d)

  expect_message(
    res <- db_export_facs(con, prep2, data_dir = d$data_dir, overwrite = FALSE),
    "nothing changed"
  )
  expect_null(res)
  expect_equal(unique(DBI::dbGetQuery(con, "SELECT value FROM facs_measurements")$value), 100)
  expect_equal(n_rows(con, "import_log"), 1)
})

test_that("existing root: overwrite = TRUE replaces data and links the log", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep1 <- suppressMessages(local_loaded_prep(d, two_mice))
  ex <- tibble::tibble(mouse_ID = "25-7-2", reason = "clots")
  id1 <- suppressMessages(db_export_facs(con, prep1, ex, data_dir = d$data_dir))
  prep2 <- second_version(d)

  id2 <- suppressMessages(
    db_export_facs(con, prep2, data_dir = d$data_dir, overwrite = TRUE)
  )

  m <- DBI::dbGetQuery(con, "SELECT DISTINCT value, source_file FROM facs_measurements")
  expect_equal(m$value, 200)
  expect_equal(m$source_file, "25-7_exp/raw/FACS/overview_25-7_20261105.wsp")
  expect_equal(n_rows(con, "facs_stains"), 2)
  log <- DBI::dbGetQuery(con, "SELECT * FROM import_log ORDER BY import_id")
  expect_equal(nrow(log), 2)
  expect_equal(log$supersedes_import_id[[2]], id1)
  exc <- DBI::dbGetQuery(con, "SELECT import_id FROM import_exclusions")
  expect_equal(exc$import_id, id1)
})

test_that("re-exporting the same version overwrites it in place", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep <- suppressMessages(local_loaded_prep(d, two_mice))
  suppressMessages(db_export_facs(con, prep, data_dir = d$data_dir))
  expect_no_error(suppressMessages(
    db_export_facs(con, prep, data_dir = d$data_dir, overwrite = TRUE)
  ))
  expect_equal(n_rows(con, "facs_measurements"), 8)
})

test_that("warns when a newer dated wsp of the root exists on disk", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep <- suppressMessages(local_loaded_prep(d, two_mice))
  file.create(file.path(d$facs_dir, "overview_25-7_20261201.wsp"))
  expect_warning(
    suppressMessages(db_export_facs(con, prep, data_dir = d$data_dir)),
    "newer"
  )
})

test_that("overwrite = NULL outside an interactive session refuses to overwrite", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep1 <- suppressMessages(local_loaded_prep(d, two_mice))
  suppressMessages(db_export_facs(con, prep1, data_dir = d$data_dir))
  prep2 <- second_version(d)
  expect_error(db_export_facs(con, prep2, data_dir = d$data_dir), "overwrite")
  expect_equal(n_rows(con, "import_log"), 1)
})

test_that("an exclusion row naming an unknown panel stops the export", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep <- suppressMessages(local_loaded_prep(d, two_mice))
  ex <- tibble::tibble(mouse_ID = "25-7-2", panel = "Overview", reason = "case typo")
  expect_error(db_export_facs(con, prep, ex, data_dir = d$data_dir), "Overview")
  expect_equal(n_rows(con, "import_log"), 0)
  # a registered other FACS panel is fine (row ignored here)
  db_write_assay(con, "treg", domain = "facs")
  ex$panel <- "treg"
  expect_no_error(suppressMessages(
    db_export_facs(con, prep, ex, data_dir = d$data_dir)
  ))
})
