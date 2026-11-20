local_elisa_dir <- function(env = parent.frame()) {
  d <- withr::local_tempdir(.local_envir = env)
  path <- file.path(d, "elisa_curated.xlsx")
  file.copy(test_path("..", "fixtures", "elisa_curated.xlsx"), path)
  path
}

read_elisa_fixture <- function(path) {
  suppressMessages(elisa_read_curated(path, tissue = "spleen"))
}

test_that("elisa_prep_save writes list(data, info) to processed/", {
  path <- local_elisa_dir()
  out <- suppressMessages(elisa_prep_save(read_elisa_fixture(path), path, "restim"))
  expect_equal(out, file.path(dirname(path), "processed", "elisa_curated_prep.rds"))
  prep <- readRDS(out)
  expect_named(prep, c("data", "info"))
  expect_equal(prep$info$root, "elisa_curated")
  expect_equal(prep$info$assay, "restim")
  expect_equal(prep$info$qc$expected_mice, c("M1", "M2"))
  expect_equal(prep$info$qc$max_cv, 20)
})

test_that("elisa_prep_save errors on missing columns and bad QC columns", {
  path <- local_elisa_dir()
  d <- read_elisa_fixture(path)
  expect_error(elisa_prep_save(dplyr::select(d, !lloq), path, "a"), "lloq")
  expect_error(elisa_prep_save(d, path, "a", group_col = "nope"), "group_col")
})

test_that("elisa_prep_load asks for a rerun without a report", {
  path <- local_elisa_dir()
  suppressMessages(elisa_prep_save(read_elisa_fixture(path), path, "restim"))
  expect_error(elisa_prep_load(path), "final <- FALSE")
})

test_that("elisa_prep_load returns paths and refuses a touched xlsx", {
  path <- local_elisa_dir()
  Sys.setFileTime(path, Sys.time() - 3600)
  suppressMessages(elisa_prep_save(read_elisa_fixture(path), path, "restim"))
  paths <- elisa_prep_paths_(path)
  dir.create(dirname(paths$qc))
  writeLines("<html></html>", paths$qc)
  Sys.setFileTime(paths$qc, Sys.time() + 10)
  prep <- elisa_prep_load(path)
  expect_equal(prep$paths$qc, paths$qc)
  expect_s3_class(prep$data, "tbl_df")

  Sys.setFileTime(path, Sys.time() + 60)
  expect_error(elisa_prep_load(path), "newer than its prep result")
})

test_that("elisa_prep_load refuses a report older than the prep", {
  path <- local_elisa_dir()
  Sys.setFileTime(path, Sys.time() - 3600)
  suppressMessages(elisa_prep_save(read_elisa_fixture(path), path, "restim"))
  paths <- elisa_prep_paths_(path)
  dir.create(dirname(paths$qc))
  writeLines("x", paths$qc)
  Sys.setFileTime(paths$qc, Sys.time() - 1800)
  expect_error(elisa_prep_load(path), "older than the prep")
})

# ---- elisa_prep_run ----------------------------------------------------------

test_that("elisa_prep_run(final = FALSE) saves and renders the report", {
  skip_if_not(rmarkdown::pandoc_available())
  path <- local_elisa_dir()
  mice <- tibble::tibble(mouse_ID = c("M1", "M2"), group = c("a", "b"))
  out <- suppressMessages(elisa_prep_run(
    FALSE, path, "restim", tissue = "spleen", mice = mice, group_col = "group",
    max_cv = 10
  ))
  expect_true(file.exists(out))
  prep <- readRDS(elisa_prep_paths_(path)$prep)
  expect_equal(prep$info$qc$group_col, "group")
  expect_equal(prep$info$qc$max_cv, 10)
  expect_true("group" %in% names(prep$data))
})

test_that("elisa_prep_save stores end_dates; elisa_prep_run fills them from death_date", {
  path <- local_elisa_dir()
  out <- suppressMessages(elisa_prep_save(
    read_elisa_fixture(path), path, "restim",
    end_dates = c(M1 = as.Date("2026-06-10"))))
  expect_equal(readRDS(out)$info$end_dates, c(M1 = as.Date("2026-06-10")))
  out <- suppressMessages(elisa_prep_save(read_elisa_fixture(path), path, "restim"))
  expect_null(readRDS(out)$info$end_dates)

  skip_if_not(rmarkdown::pandoc_available())
  mice <- tibble::tibble(mouse_ID = c("M1", "M2"), group = c("a", "b"),
                         death_date = as.Date(c("2026-06-10", "2026-06-11")))
  suppressMessages(elisa_prep_run(FALSE, path, "restim", tissue = "spleen",
                                  mice = mice, group_col = "group"))
  prep <- readRDS(elisa_prep_paths_(path)$prep)
  expect_equal(prep$info$end_dates,
               c(M1 = as.Date("2026-06-10"), M2 = as.Date("2026-06-11")))
  expect_false("death_date" %in% names(prep$data))
  # explicit end_dates in ... win over mice$death_date
  suppressMessages(elisa_prep_run(FALSE, path, "restim", tissue = "spleen",
                                  mice = mice,
                                  end_dates = c(M2 = as.Date("2026-07-01"))))
  prep <- readRDS(elisa_prep_paths_(path)$prep)
  expect_equal(prep$info$end_dates, c(M2 = as.Date("2026-07-01")))
})

test_that("elisa_prep_run(final = TRUE) needs a connection", {
  path <- local_elisa_dir()
  expect_error(elisa_prep_run(TRUE, path, "restim"), "con")
  expect_error(elisa_prep_run("yes", path, "restim"), "final")
})

# ---- final run records the exclusion rule set ---------------------------------

elisa_final_setup <- function(env = parent.frame()) {
  con <- local_test_db(env)
  db_write_experiment(con, experiment_code = "26-1")
  db_write_subjects(con, tibble::tibble(
    mouse_id = c("M1", "M2"), experiment_code = "26-1",
    end_date = as.Date("2026-06-10")
  ))
  path <- local_elisa_dir(env)
  Sys.setFileTime(path, Sys.time() - 3600)
  suppressMessages(elisa_prep_save(read_elisa_fixture(path), path, "restim"))
  paths <- elisa_prep_paths_(path)
  dir.create(dirname(paths$qc))
  writeLines("<html></html>", paths$qc)
  Sys.setFileTime(paths$qc, Sys.time() + 10)
  list(con = con, path = path, data_dir = dirname(path))
}

test_that("elisa_prep_run(final = TRUE) writes the exclusion set", {
  s <- elisa_final_setup()
  ex <- tibble::tibble(mouse_ID = "M2", reason = "bad plate")
  suppressMessages(elisa_prep_run(TRUE, s$path, "restim", exclusions = ex,
                                  con = s$con, data_dir = s$data_dir))
  expect_equal(count_rows_(s$con, "exclusion_sets"), 1)
  expect_equal(count_rows_(s$con, "exclusion_rules"), 1)
  expect_equal(
    DBI::dbGetQuery(s$con, "SELECT reason FROM exclusion_rules")$reason,
    "bad plate")
})

test_that("elisa_prep_run(final = TRUE) without exclusions keeps the rule set", {
  s <- elisa_final_setup()
  suppressMessages(db_write_exclusions(
    s$con, tibble::tibble(mouse_ID = "M2", reason = "bad plate"), "26-1"))
  msgs <- capture_messages(
    elisa_prep_run(TRUE, s$path, "restim", con = s$con, data_dir = s$data_dir))
  expect_true(any(grepl("exclusions not given -- rule set of 26-1 left unchanged.",
                        msgs, fixed = TRUE)))
  expect_equal(count_rows_(s$con, "exclusion_sets"), 1)
  expect_equal(count_rows_(s$con, "exclusion_rules"), 1)
  expect_equal(count_rows_(s$con, "elisa_measurements") > 0, TRUE)
})

test_that("elisa_prep_run meta_path fills source_file", {
  s <- elisa_final_setup()
  meta_file <- file.path(s$data_dir, "meta.xlsx")
  file.create(meta_file)
  ex <- tibble::tibble(mouse_ID = "M2", reason = "bad plate")
  suppressMessages(elisa_prep_run(TRUE, s$path, "restim", exclusions = ex,
                                  con = s$con, meta_path = meta_file,
                                  data_dir = s$data_dir))
  expect_equal(
    DBI::dbGetQuery(s$con, "SELECT source_file FROM exclusion_sets")$source_file,
    "meta.xlsx")
})
