seed_histo_db <- function(con, end_dates = rep(as.Date("2026-06-10"), 4)) {
  db_write_experiment(con, experiment_code = "26-1")
  db_write_subjects(con, tibble::tibble(
    mouse_id = paste0("M", 1:4), experiment_code = "26-1", end_date = end_dates
  ))
  invisible(con)
}

histo_fixture_data <- function() {
  suppressMessages(histo_read_curated(
    test_path("..", "fixtures", "histo_curated.xlsx"), tissue = "kidney")) |>
    histo_calc_ratio("crescents", "gloms", "pct_crescents")
}

histo_prep_in_root <- function(data = NULL, env = parent.frame()) {
  data_dir <- withr::local_tempdir(.local_envir = env)
  folder <- file.path(data_dir, "26-1_exp", "raw", "histo")
  dir.create(folder, recursive = TRUE)
  path <- file.path(folder, "histo_curated.xlsx")
  file.copy(test_path("..", "fixtures", "histo_curated.xlsx"), path)
  Sys.setFileTime(path, Sys.time() - 3600)
  if (is.null(data)) data <- histo_fixture_data()
  suppressMessages(histo_prep_save(data, path, "PAS"))
  paths <- histo_prep_paths_(path)
  dir.create(dirname(paths$qc))
  writeLines("<html></html>", paths$qc)
  Sys.setFileTime(paths$qc, Sys.time() + 10)
  list(prep = histo_prep_load(path), data_dir = data_dir)
}

histo_export_quiet <- function(...) suppressMessages(db_export_histo(...))

test_that("histo export writes every row incl. NA values and pct_crescents", {
  con <- local_test_db(); seed_histo_db(con)
  p <- histo_prep_in_root()
  id <- histo_export_quiet(con, p$prep, data_dir = p$data_dir)
  expect_equal(count_rows_(con, "histo_measurements"), nrow(p$prep$data))
  res <- dplyr::collect(db_query_histo(con))
  expect_equal(sort(unique(as.character(res$collected_at))), "2026-06-10")
  expect_true(all(res$import_unit == "histo_curated"))
  expect_true("pct_crescents" %in% res$metric)
  m3 <- res[res$mouse_id == "M3", ]
  expect_equal(nrow(m3), 4)
  expect_true(all(is.na(m3$value)))
  log <- DBI::dbGetQuery(con, "SELECT * FROM import_log WHERE import_id = ?",
                         params = list(id))
  expect_equal(log$source_file, "26-1_exp/raw/histo/histo_curated.xlsx")
  expect_equal(log$qc_report, "26-1_exp/raw/histo/qc/histo_curated_qc.html")
  expect_true(log$qc_passed)
})

test_that("histo overwrite replaces the root; FALSE changes nothing", {
  con <- local_test_db(); seed_histo_db(con)
  p <- histo_prep_in_root()
  id1 <- histo_export_quiet(con, p$prep, data_dir = p$data_dir)
  n <- count_rows_(con, "histo_measurements")

  changed <- histo_fixture_data()
  changed$value[which(!is.na(changed$value))[1]] <- 12345
  p2 <- histo_prep_in_root(changed)
  expect_null(histo_export_quiet(con, p2$prep, data_dir = p2$data_dir,
                                 overwrite = FALSE))
  expect_equal(DBI::dbGetQuery(
    con, "SELECT COUNT(*) AS n FROM histo_measurements WHERE value = 12345")$n, 0)

  id2 <- histo_export_quiet(con, p2$prep, data_dir = p2$data_dir, overwrite = TRUE)
  expect_equal(count_rows_(con, "histo_measurements"), n)
  expect_equal(DBI::dbGetQuery(
    con, "SELECT COUNT(*) AS n FROM histo_measurements WHERE value = 12345")$n, 1)
  sup <- DBI::dbGetQuery(
    con, "SELECT supersedes_import_id AS s FROM import_log WHERE import_id = ?",
    params = list(id2))$s
  expect_equal(sup, id1)
})

test_that("histo exclusion drops the sample and logs one row", {
  con <- local_test_db(); seed_histo_db(con)
  p <- histo_prep_in_root()
  ex <- tibble::tibble(mouse_ID = "M2", reason = "bad slide")
  id <- histo_export_quiet(con, p$prep, ex, data_dir = p$data_dir)
  expect_equal(DBI::dbGetQuery(con, "
    SELECT COUNT(*) AS n FROM histo_measurements e
    JOIN samples s ON s.sample_id = e.sample_id
    JOIN subjects u ON u.subject_id = s.subject_id WHERE u.mouse_id = 'M2'")$n, 0)
  log <- DBI::dbGetQuery(con, "SELECT * FROM import_exclusions WHERE import_id = ?",
                         params = list(id))
  expect_equal(nrow(log), 1)
  expect_equal(log$mouse_id, "M2")
  expect_equal(log$n_files_matched, 1)
})

test_that("histo unmatched exclusion errors and writes nothing", {
  con <- local_test_db(); seed_histo_db(con)
  p <- histo_prep_in_root()
  ex <- tibble::tibble(mouse_ID = "M9", reason = "typo")
  expect_error(histo_export_quiet(con, p$prep, ex, data_dir = p$data_dir),
               "match no sample")
  expect_equal(count_rows_(con, "histo_measurements"), 0)
  expect_equal(count_rows_(con, "import_log"), 0)
})

test_that("histo dated row equal to end_date plus an empty-date row errors", {
  con <- local_test_db(); seed_histo_db(con)
  d <- histo_fixture_data()
  dup <- dplyr::filter(d, mouse_ID == "M2")
  dup$collected_at <- as.Date("2026-06-10")
  p <- histo_prep_in_root(dplyr::bind_rows(d, dup))
  expect_error(histo_export_quiet(con, p$prep, data_dir = p$data_dir), "twice")
  expect_equal(count_rows_(con, "histo_measurements"), 0)
  expect_equal(count_rows_(con, "samples"), 0)
})

test_that("histo mice without end_date error before anything is written", {
  con <- local_test_db()
  seed_histo_db(con, as.Date(c("2026-06-10", NA, "2026-06-10", "2026-06-10")))
  p <- histo_prep_in_root()
  expect_error(histo_export_quiet(con, p$prep, data_dir = p$data_dir), "'M2'")
  expect_equal(count_rows_(con, "samples"), 0)
})

test_that("prep must come from histo_prep_load()", {
  con <- local_test_db()
  expect_error(db_export_histo(con, list(info = list())), "histo_prep_load")
})
