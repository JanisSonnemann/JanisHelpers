histo_rows <- function(con, value, extra = NULL) {
  sid <- lookup_sample_id_(con, "25-7-1", "spleen")
  aid <- lookup_assay_id_(con, "hist", "histo")
  tibble::tibble(
    sample_id = as.integer(sid), assay_id = as.integer(aid),
    import_unit = "histo_25-7", metric = c("area", "count"), value = value,
    source_file = "histo/histo_25-7.xlsx"
  )
}

setup_histo <- function(con) {
  seed_minimal_fixture(con)
  db_write_assay(con, "hist", domain = "histo")
}

do_replace <- function(con, rows, applicable = NULL) {
  replace_import_(con, "histo", "histo_measurements", rows, "histo_25-7",
                  "histo/histo_25-7.xlsx", "histo/qc/histo_25-7_qc.html",
                  applicable)
}

test_that("replace_import_() keeps only the latest rows and chains the log", {
  con <- local_test_db(); setup_histo(con)
  id1 <- do_replace(con, histo_rows(con, c(1, 2)))
  id2 <- do_replace(con, histo_rows(con, c(3, 4)))
  got <- DBI::dbGetQuery(con, "SELECT metric, value FROM histo_measurements ORDER BY metric")
  expect_equal(got$value, c(3, 4))
  log <- DBI::dbGetQuery(con, "SELECT * FROM import_log ORDER BY import_id")
  expect_equal(log$import_id, c(id1, id2))
  expect_true(is.na(log$supersedes_import_id[1]))
  expect_equal(log$supersedes_import_id[2], id1)
  expect_true(all(log$qc_passed))
  expect_equal(log$rows_written, c(2L, 2L))
})

test_that("replace_import_() stores applicable exclusions", {
  con <- local_test_db(); setup_histo(con)
  app <- tibble::tibble(
    mouse_ID = "25-7-9", tissue = NA_character_, panel = "hist",
    collected_at = as.Date(NA), reason = "bad", date = as.Date("2026-01-02"),
    n_samples_matched = 0L
  )
  do_replace(con, histo_rows(con, c(1, 2)), app)
  got <- DBI::dbGetQuery(con, "SELECT * FROM import_exclusions")
  expect_equal(nrow(got), 1)
  expect_equal(got$n_files_matched, 0L)
  expect_true(is.na(got$file_name))
})

test_that("a failing replace leaves table and log unchanged", {
  con <- local_test_db(); setup_histo(con)
  do_replace(con, histo_rows(con, c(1, 2)))
  bad <- dplyr::bind_rows(histo_rows(con, c(5, 6)), histo_rows(con, c(5, 6)))
  expect_error(do_replace(con, bad))
  got <- DBI::dbGetQuery(con, "SELECT metric, value FROM histo_measurements ORDER BY metric")
  expect_equal(got$value, c(1, 2))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM import_log")$n, 1)
})

test_that("replace_import_() rejects an unknown table", {
  con <- local_test_db(); setup_histo(con)
  expect_error(
    replace_import_(con, "histo", "samples; DROP", tibble::tibble(), "r", "s", "q", NULL),
    "table"
  )
})

test_that("confirm_replace_() proceeds on a new root, honours overwrite", {
  con <- local_test_db(); setup_histo(con)
  expect_true(confirm_replace_(con, "histo_measurements", "histo_25-7", NULL))
  do_replace(con, histo_rows(con, c(1, 2)))
  expect_false(confirm_replace_(con, "histo_measurements", "histo_25-7", FALSE))
  expect_true(confirm_replace_(con, "histo_measurements", "histo_25-7", TRUE))
  expect_error(confirm_replace_(con, "histo_measurements", "histo_25-7", NULL), "overwrite")
})

test_that("check_exclusion_panels_() accepts registered assays of any domain", {
  con <- local_test_db(); setup_histo(con)
  ex <- tibble::tibble(mouse_ID = "a", panel = "hist", reason = "x")
  expect_no_error(check_exclusion_panels_(con, ex, "titer"))
  ex$panel <- "Treg"
  expect_warning(check_exclusion_panels_(con, ex, "titer"), "'Treg'.*Known panels / assays: titer")
  expect_no_error(check_exclusion_panels_(con, NULL, "titer"))
})
