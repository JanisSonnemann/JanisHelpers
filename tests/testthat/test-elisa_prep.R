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

test_that("elisa_prep_run(final = TRUE) needs a connection", {
  path <- local_elisa_dir()
  expect_error(elisa_prep_run(TRUE, path, "restim"), "con")
  expect_error(elisa_prep_run("yes", path, "restim"), "final")
})
