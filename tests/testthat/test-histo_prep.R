local_histo_dir <- function(env = parent.frame()) {
  d <- withr::local_tempdir(.local_envir = env)
  path <- file.path(d, "histo_curated.xlsx")
  file.copy(test_path("..", "fixtures", "histo_curated.xlsx"), path)
  path
}

read_histo_fixture <- function(path) {
  suppressMessages(histo_read_curated(path, tissue = "kidney"))
}

test_that("histo_prep_paths_ derives root, prep and report paths", {
  p <- histo_prep_paths_("/x/y/kidney_26-1.xlsx")
  expect_equal(p$root, "kidney_26-1")
  expect_equal(p$prep, "/x/y/processed/kidney_26-1_prep.rds")
  expect_equal(p$qc, "/x/y/qc/kidney_26-1_qc.html")
})

test_that("histo_prep_save writes list(data, info) to processed/", {
  path <- local_histo_dir()
  out <- suppressMessages(histo_prep_save(read_histo_fixture(path), path, "PAS",
                                          rules = "crescents <= gloms"))
  expect_equal(out, file.path(dirname(path), "processed", "histo_curated_prep.rds"))
  prep <- readRDS(out)
  expect_named(prep, c("data", "info"))
  expect_equal(prep$info$root, "histo_curated")
  expect_equal(prep$info$assay, "PAS")
  expect_equal(prep$info$qc$expected_mice, c("M1", "M2", "M3", "M4"))
  expect_equal(prep$info$qc$rules, "crescents <= gloms")
  expect_equal(
    readRDS(suppressMessages(histo_prep_save(read_histo_fixture(path), path, "PAS")))$info$qc$rules,
    character())
})

test_that("histo_prep_save errors on bad input and warns on duplicates", {
  path <- local_histo_dir()
  d <- read_histo_fixture(path)
  expect_error(histo_prep_save(dplyr::select(d, !metric), path, "a"), "metric")
  expect_error(histo_prep_save(d, path, "a", group_col = "nope"), "group_col")
  expect_error(histo_prep_save(d, path, "a", rules = 1), "rules")
  expect_warning(
    suppressMessages(histo_prep_save(dplyr::bind_rows(d, d[1, ]), path, "a")),
    "occur twice")
})

test_that("histo_prep_load asks for a rerun without a report", {
  path <- local_histo_dir()
  suppressMessages(histo_prep_save(read_histo_fixture(path), path, "PAS"))
  expect_error(histo_prep_load(path), "final <- FALSE")
})

test_that("histo_prep_load returns paths and refuses a touched xlsx", {
  path <- local_histo_dir()
  Sys.setFileTime(path, Sys.time() - 3600)
  suppressMessages(histo_prep_save(read_histo_fixture(path), path, "PAS"))
  paths <- histo_prep_paths_(path)
  dir.create(dirname(paths$qc))
  writeLines("<html></html>", paths$qc)
  Sys.setFileTime(paths$qc, Sys.time() + 10)
  prep <- histo_prep_load(path)
  expect_equal(prep$paths$qc, paths$qc)
  expect_s3_class(prep$data, "tbl_df")

  Sys.setFileTime(path, Sys.time() + 60)
  expect_error(histo_prep_load(path), "newer than its prep result")
})

test_that("histo_prep_load refuses a report older than the prep", {
  path <- local_histo_dir()
  Sys.setFileTime(path, Sys.time() - 3600)
  suppressMessages(histo_prep_save(read_histo_fixture(path), path, "PAS"))
  paths <- histo_prep_paths_(path)
  dir.create(dirname(paths$qc))
  writeLines("x", paths$qc)
  Sys.setFileTime(paths$qc, Sys.time() - 1800)
  expect_error(histo_prep_load(path), "older than the prep")
})

# ---- histo_prep_run ----------------------------------------------------------

test_that("histo_prep_run(final = FALSE) saves and renders the report", {
  skip_if_not(rmarkdown::pandoc_available())
  path <- local_histo_dir()
  mice <- tibble::tibble(mouse_ID = c("M1", "M2", "M3", "M4"),
                         group = c("a", "b", "a", "b"))
  out <- suppressMessages(histo_prep_run(
    FALSE, path, "PAS", tissue = "kidney", mice = mice, group_col = "group",
    prepare = function(d) histo_calc_ratio(d, "crescents", "gloms", "pct_crescents"),
    rules = "crescents <= gloms"
  ))
  expect_true(file.exists(out))
  prep <- readRDS(histo_prep_paths_(path)$prep)
  expect_equal(prep$info$qc$group_col, "group")
  expect_equal(prep$info$qc$rules, "crescents <= gloms")
  expect_true("pct_crescents" %in% prep$data$metric)
})

test_that("histo_prep_run(final = TRUE) needs a connection", {
  path <- local_histo_dir()
  expect_error(histo_prep_run(TRUE, path, "PAS"), "con")
  expect_error(histo_prep_run("yes", path, "PAS"), "final")
})
