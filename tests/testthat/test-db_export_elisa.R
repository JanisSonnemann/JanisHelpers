seed_elisa_db <- function(con) {
  db_write_experiment(con, experiment_code = "26-1")
  db_write_subjects(con, tibble::tibble(
    mouse_id = c("M1", "M2"), experiment_code = "26-1",
    end_date = as.Date("2026-06-10")
  ))
  invisible(con)
}

# Writes the data (default: the fixture) as a curated file's prep result plus a
# fake QC report inside a temp data/ root; returns the loaded prep + data_dir.
elisa_prep_in_root <- function(data = NULL, env = parent.frame()) {
  data_dir <- withr::local_tempdir(.local_envir = env)
  folder <- file.path(data_dir, "26-1_exp", "raw", "ELISA")
  dir.create(folder, recursive = TRUE)
  path <- file.path(folder, "elisa_curated.xlsx")
  file.copy(test_path("..", "fixtures", "elisa_curated.xlsx"), path)
  Sys.setFileTime(path, Sys.time() - 3600)
  if (is.null(data)) {
    data <- suppressMessages(elisa_read_curated(path, tissue = "spleen"))
  }
  suppressMessages(elisa_prep_save(data, path, "restim"))
  paths <- elisa_prep_paths_(path)
  dir.create(dirname(paths$qc))
  writeLines("<html></html>", paths$qc)
  Sys.setFileTime(paths$qc, Sys.time() + 10)
  list(prep = elisa_prep_load(path), data_dir = data_dir)
}

fixture_data <- function() {
  suppressMessages(elisa_read_curated(
    test_path("..", "fixtures", "elisa_curated.xlsx"), tissue = "spleen"))
}

export_quiet <- function(...) suppressMessages(db_export_elisa(...))

test_that("exports every row and resolves empty dates to end_date", {
  con <- local_test_db(); seed_elisa_db(con)
  p <- elisa_prep_in_root()
  id <- export_quiet(con, p$prep, data_dir = p$data_dir)
  expect_equal(count_rows_(con, "elisa_measurements"), nrow(p$prep$data))
  res <- dplyr::collect(db_query_elisa(con))
  expect_equal(sort(unique(as.character(res$collected_at))),
               c("2026-05-06", "2026-06-10"))
  expect_true(all(res$import_unit == "elisa_curated"))
  log <- DBI::dbGetQuery(con, "SELECT * FROM import_log WHERE import_id = ?",
                         params = list(id))
  expect_equal(log$source_file, "26-1_exp/raw/ELISA/elisa_curated.xlsx")
  expect_equal(log$qc_report, "26-1_exp/raw/ELISA/qc/elisa_curated_qc.html")
  expect_true(log$qc_passed)
})

test_that("overwrite replaces the root; FALSE changes nothing", {
  con <- local_test_db(); seed_elisa_db(con)
  p <- elisa_prep_in_root()
  id1 <- export_quiet(con, p$prep, data_dir = p$data_dir)
  n <- count_rows_(con, "elisa_measurements")

  changed <- fixture_data()
  changed$value[which(!is.na(changed$value))[1]] <- 12345
  p2 <- elisa_prep_in_root(changed)
  expect_null(export_quiet(con, p2$prep, data_dir = p2$data_dir, overwrite = FALSE))
  expect_equal(DBI::dbGetQuery(
    con, "SELECT COUNT(*) AS n FROM elisa_measurements WHERE value = 12345")$n, 0)

  id2 <- export_quiet(con, p2$prep, data_dir = p2$data_dir, overwrite = TRUE)
  expect_equal(count_rows_(con, "elisa_measurements"), n)
  expect_equal(DBI::dbGetQuery(
    con, "SELECT COUNT(*) AS n FROM elisa_measurements WHERE value = 12345")$n, 1)
  sup <- DBI::dbGetQuery(con, "SELECT supersedes_import_id AS s FROM import_log WHERE import_id = ?",
                         params = list(id2))$s
  expect_equal(sup, id1)
})

test_that("exclusion rule does not drop the sample and logs no exclusion", {
  con <- local_test_db(); seed_elisa_db(con)
  p <- elisa_prep_in_root()
  ex <- tibble::tibble(mouse_ID = "M2", reason = "bad plate")
  expect_message(
    id <- db_export_elisa(con, p$prep, ex, data_dir = p$data_dir),
    "1 sample\\(s\\) match an exclusion rule \\(stored, flagged on read\\)"
  )
  expect_equal(DBI::dbGetQuery(con, "
    SELECT COUNT(*) AS n FROM elisa_measurements e
    JOIN samples s ON s.sample_id = e.sample_id
    JOIN subjects u ON u.subject_id = s.subject_id WHERE u.mouse_id = 'M2'")$n > 0, TRUE)
  expect_equal(count_rows_(con, "elisa_measurements"), nrow(p$prep$data))
  expect_equal(count_rows_(con, "import_exclusions"), 0)
})

test_that("an unmatched exclusion gives a message and the export succeeds", {
  con <- local_test_db(); seed_elisa_db(con)
  p <- elisa_prep_in_root()
  ex <- tibble::tibble(mouse_ID = "M9", reason = "typo")
  expect_message(export_quiet_msg <- db_export_elisa(
    con, p$prep, ex, data_dir = p$data_dir), "match no sample.*M9")
  expect_equal(count_rows_(con, "elisa_measurements"), nrow(p$prep$data))
  expect_equal(count_rows_(con, "import_log"), 1)
})

test_that("a rule dated on the end date matches an end-point sample", {
  con <- local_test_db(); seed_elisa_db(con)
  p <- elisa_prep_in_root()
  ex <- tibble::tibble(mouse_ID = "M2", collected_at = as.Date("2026-06-10"),
                       reason = "end bleed")
  msgs <- testthat::capture_messages(
    db_export_elisa(con, p$prep, ex, data_dir = p$data_dir))
  expect_false(any(grepl("match no sample", msgs)))
  expect_true(any(grepl("1 sample\\(s\\) match an exclusion rule", msgs)))
  # an unmatched dated row carries the new hint
  ex$mouse_ID <- "M9"
  expect_message(
    db_export_elisa(con, p$prep, ex, data_dir = p$data_dir, overwrite = TRUE),
    "stored with the mouse's end date -- a rule dated on that day matches them"
  )
})

test_that("a dated row equal to end_date plus an empty-date row errors", {
  con <- local_test_db(); seed_elisa_db(con)
  d <- fixture_data()
  dup <- dplyr::filter(d, mouse_ID == "M2")
  dup$collected_at <- as.Date("2026-06-10")
  p <- elisa_prep_in_root(dplyr::bind_rows(d, dup))
  expect_error(export_quiet(con, p$prep, data_dir = p$data_dir), "twice")
  expect_equal(count_rows_(con, "elisa_measurements"), 0)
  expect_equal(count_rows_(con, "samples"), 0)
})

test_that("mice without end_date error before anything is written", {
  con <- local_test_db()
  db_write_experiment(con, experiment_code = "26-1")
  db_write_subjects(con, tibble::tibble(
    mouse_id = c("M1", "M2"), experiment_code = "26-1",
    end_date = as.Date(c("2026-06-10", NA))))
  p <- elisa_prep_in_root()
  expect_error(export_quiet(con, p$prep, data_dir = p$data_dir), "'M2'")
  expect_equal(count_rows_(con, "samples"), 0)
})

test_that("prep must come from elisa_prep_load()", {
  con <- local_test_db()
  expect_error(db_export_elisa(con, list(info = list())), "elisa_prep_load")
})
