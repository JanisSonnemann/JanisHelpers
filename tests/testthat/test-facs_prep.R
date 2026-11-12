two_files <- tibble::tibble(
  file_name = c("a.fcs", "b.fcs"),
  mouse_ID = c("25-7-1", "25-7-2"),
  tissue = "spleen"
)

test_that("facs_prep_save writes rds with data + info and returns its path", {
  d <- local_prep_dirs()
  path <- facs_prep_save(fake_prep_data(two_files), d$wsp_path, "overview", d$panel_dir)
  expect_equal(
    path,
    file.path(d$panel_dir, "processed", "overview_25-7_20261101_prep.rds")
  )
  prep <- readRDS(path)
  expect_equal(prep$info$wsp_root, "overview_25-7")
  expect_equal(prep$info$wsp_version, as.Date("2026-11-01"))
  expect_equal(prep$info$assay, "overview")
  expect_equal(prep$info$wsp_path, d$wsp_path)
  expect_s3_class(prep$info$created_at, "POSIXct")
  expect_type(prep$info$package_version, "character")
  expect_equal(unique(prep$data$condition), "none")
  expect_setequal(unique(prep$data$file_name), c("a.fcs", "b.fcs"))
})

test_that("facs_prep_save drops files without mouse_ID or tissue", {
  d <- local_prep_dirs()
  files <- tibble::tibble(
    file_name = c("a.fcs", "comp.fcs"),
    mouse_ID = c("25-7-1", NA), tissue = c("spleen", NA)
  )
  path <- facs_prep_save(fake_prep_data(files), d$wsp_path, "overview", d$panel_dir)
  expect_equal(unique(readRDS(path)$data$file_name), "a.fcs")
})

test_that("facs_prep_save takes condition from the keyword column", {
  d <- local_prep_dirs()
  files <- tibble::tibble(
    file_name = c("a.fcs", "b.fcs"), mouse_ID = "25-7-1", tissue = "spleen",
    restimulation = c("unstim", "MPO")
  )
  path <- facs_prep_save(fake_prep_data(files), d$wsp_path, "overview",
                         d$panel_dir, condition_keyword = "restimulation")
  data <- readRDS(path)$data
  expect_equal(data$condition, data$restimulation)
  expect_equal(readRDS(path)$info$condition_keyword, "restimulation")

  files$restimulation[2] <- NA
  expect_error(
    facs_prep_save(fake_prep_data(files), d$wsp_path, "overview",
                   d$panel_dir, condition_keyword = "restimulation"),
    "restimulation"
  )
})

test_that("facs_prep_save errors on missing columns", {
  d <- local_prep_dirs()
  data <- dplyr::select(fake_prep_data(two_files), !tissue)
  expect_error(
    facs_prep_save(data, d$wsp_path, "overview", d$panel_dir),
    "tissue"
  )
})

test_that("facs_prep_save rejects a wsp name without _YYYYMMDD", {
  d <- local_prep_dirs()
  bad <- file.path(d$facs_dir, "overview_25-7.wsp")
  file.create(bad)
  expect_error(
    facs_prep_save(fake_prep_data(two_files), bad, "overview", d$panel_dir),
    "_YYYYMMDD"
  )
})

test_that("facs_prep_save warns on files sharing mouse x tissue x condition", {
  d <- local_prep_dirs()
  files <- tibble::tibble(file_name = c("a.fcs", "b.fcs"),
                          mouse_ID = "25-7-1", tissue = "spleen")
  expect_warning(
    path <- facs_prep_save(fake_prep_data(files), d$wsp_path, "overview", d$panel_dir),
    "share"
  )
  expect_true(file.exists(path))
})

test_that("facs_check_file_clashes lists every clashing file", {
  data <- tibble::tibble(
    file_name = c("a.fcs", "b.fcs", "c.fcs"),
    mouse_ID = c("25-7-1", "25-7-1", "25-7-2"),
    tissue = "spleen", condition = "none"
  )
  clashes <- facs_check_file_clashes(data)
  expect_named(clashes, c("file_name", "mouse_ID", "tissue", "condition"))
  expect_setequal(clashes$file_name, c("a.fcs", "b.fcs"))
  expect_equal(nrow(facs_check_file_clashes(data[c(1, 3), ])), 0)
})

test_that("facs_prep_load returns data, info and paths", {
  d <- local_prep_dirs()
  prep <- local_loaded_prep(d, two_files)
  expect_named(prep, c("data", "info", "paths"))
  expect_equal(
    prep$paths$qc,
    file.path(d$panel_dir, "qc", "overview_25-7_20261101_qc.html")
  )
  expect_equal(
    prep$paths$prep,
    file.path(d$panel_dir, "processed", "overview_25-7_20261101_prep.rds")
  )
})

test_that("facs_prep_load stops when the rds or report is missing", {
  d <- local_prep_dirs()
  expect_error(facs_prep_load(d$wsp_path, d$panel_dir), "final <- FALSE")
  facs_prep_save(fake_prep_data(two_files), d$wsp_path, "overview", d$panel_dir)
  expect_error(facs_prep_load(d$wsp_path, d$panel_dir), "final <- FALSE")
})

test_that("facs_prep_load stops when the wsp is newer than the rds", {
  d <- local_prep_dirs()
  local_loaded_prep(d, two_files)
  Sys.setFileTime(d$wsp_path, Sys.time() + 60)
  expect_error(facs_prep_load(d$wsp_path, d$panel_dir), "newer")
})

test_that("the prep script template ships and parses", {
  path <- system.file("scripts", "facs_prep_template.R", package = "JanisHelpers")
  expect_true(nzchar(path))
  expect_no_error(parse(path))
})
