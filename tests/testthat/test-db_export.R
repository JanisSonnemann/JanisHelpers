two_mice <- tibble::tibble(
  file_name = c("a.fcs", "b.fcs"),
  mouse_ID = c("25-7-1", "25-7-2"),
  tissue = "spleen"
)

seed_export_fixture <- function(con) {
  db_write_experiment(con, experiment_code = "25-7")
  db_write_subjects(con, tibble::tibble(
    mouse_id = c("25-7-1", "25-7-2"), experiment_code = "25-7",
    end_date = as.Date("2026-06-10")
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

test_that("db_export_facs ignores extra meta columns in prep data", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  with_meta <- dplyr::mutate(two_mice, grp = "A", day = as.Date("2026-06-15"))
  prep <- suppressMessages(local_loaded_prep(d, with_meta))
  suppressMessages(db_export_facs(con, prep, data_dir = d$data_dir))
  expect_equal(n_rows(con, "facs_measurements"), nrow(prep$data))
})

test_that("exclusion rules do not drop files and write no import_exclusions", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep <- suppressMessages(local_loaded_prep(d, two_mice))
  ex <- tibble::tibble(mouse_ID = "25-7-2", reason = "clots",
                       date = as.Date("2026-11-10"))

  expect_message(
    db_export_facs(con, prep, ex, data_dir = d$data_dir),
    "1 file\\(s\\) match an exclusion rule \\(stored, flagged on read\\)"
  )

  expect_equal(n_rows(con, "facs_stains"), 2)
  mice <- sort(unique(dplyr::collect(db_query_facs(con))$mouse_id))
  expect_equal(mice, c("25-7-1", "25-7-2"))
  expect_equal(n_rows(con, "import_exclusions"), 0)
})

test_that("db_export_facs stores file_name per stain", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep <- suppressMessages(local_loaded_prep(d, two_mice))
  suppressMessages(db_export_facs(con, prep, data_dir = d$data_dir))
  got <- DBI::dbGetQuery(con, "SELECT file_name FROM facs_stains ORDER BY file_name")$file_name
  expect_equal(got, sort(unique(prep$data$file_name)))
  expect_false(anyNA(got))
})

test_that("a dated exclusion row (ELISA bleed) is only checked, not applied", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep <- suppressMessages(local_loaded_prep(d, two_mice))
  ex <- tibble::tibble(
    mouse_ID = c("25-7-1", "25-7-2"),
    collected_at = as.Date(c("2026-05-01", NA)),
    reason = c("hemolysed bleed", "clots")
  )
  suppressMessages(db_export_facs(con, prep, ex, data_dir = d$data_dir))
  expect_equal(n_rows(con, "facs_stains"), 2)
  # a dated row naming the FACS panel stops the export
  ex2 <- tibble::tibble(mouse_ID = "25-7-1", panel = "overview",
                        collected_at = as.Date("2026-05-01"), reason = "x")
  expect_error(db_export_facs(con, prep, ex2, data_dir = d$data_dir,
                              overwrite = TRUE), "collected_at")
})

test_that("unmatched exclusion row gives a message, export succeeds", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep <- suppressMessages(local_loaded_prep(d, two_mice))
  ex <- tibble::tibble(mouse_ID = "25-7-1", tissue = "Spleen", reason = "typo")

  expect_message(db_export_facs(con, prep, ex, data_dir = d$data_dir),
                 "match no file.*Spleen")
  expect_equal(n_rows(con, "facs_stains"), 2)
  expect_equal(n_rows(con, "import_log"), 1)
})

test_that("all files matching a rule still export", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep <- suppressMessages(local_loaded_prep(d, two_mice))
  ex <- tibble::tibble(mouse_ID = c("25-7-1", "25-7-2"), reason = "x")

  suppressMessages(db_export_facs(con, prep, ex, data_dir = d$data_dir))
  expect_equal(n_rows(con, "facs_stains"), 2)
  expect_equal(n_rows(con, "import_log"), 1)
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
  expect_equal(n_rows(con, "import_exclusions"), 0)
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
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep <- suppressMessages(local_loaded_prep(d, two_mice))
  prep$data$population <- as.character(prep$data$population)
  suppressMessages(db_export_facs(con, prep, data_dir = d$data_dir))
  rows <- DBI::dbGetQuery(con, "SELECT population_order FROM facs_measurements")
  expect_equal(nrow(rows), 8)
  expect_true(all(is.na(rows$population_order)))
})

test_that("re-exporting a root stored without order fills population_order", {
  con <- local_test_db(); seed_export_fixture(con)
  d <- local_prep_dirs()
  prep <- suppressMessages(local_loaded_prep(d, two_mice))
  suppressMessages(db_export_facs(con, prep, data_dir = d$data_dir))
  DBI::dbExecute(con, "UPDATE facs_measurements SET population_order = NULL")
  suppressMessages(db_export_facs(con, prep, data_dir = d$data_dir, overwrite = TRUE))
  expect_equal(
    DBI::dbGetQuery(con,
      "SELECT COUNT(*) AS n FROM facs_measurements WHERE population_order IS NULL")$n,
    0
  )
})
