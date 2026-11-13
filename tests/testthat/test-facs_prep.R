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

# -- QC settings ---------------------------------------------------------------

meta_files <- two_files |>
  dplyr::mutate(grp = c("A", "B"), day = as.Date(c("2026-06-15", "2026-06-16")))

test_that("facs_prep_save stores resolved qc settings in info", {
  d <- local_prep_dirs()
  path <- suppressMessages(facs_prep_save(
    fake_prep_data(meta_files), d$wsp_path, "overview", d$panel_dir,
    group_col = "grp", batch_col = "day",
    expected_mice = c("25-7-1", "25-7-2")
  ))
  qc <- readRDS(path)$info$qc
  expect_equal(qc$group_col, "grp")
  expect_equal(qc$batch_col, "day")
  expect_equal(qc$expected_mice, c("25-7-1", "25-7-2"))
  expect_equal(qc$qc_populations, c("live_CD45+", "CD4+"))
  expect_equal(qc$expected_tissues, "spleen")
  expect_equal(qc$qc_mfi_populations, character())
  expect_equal(qc$qc_metrics, character())
})

test_that("facs_prep_save errors on unknown qc columns or populations", {
  d <- local_prep_dirs()
  data <- fake_prep_data(meta_files)
  save <- function(...) {
    suppressMessages(facs_prep_save(data, d$wsp_path, "overview", d$panel_dir, ...))
  }
  expect_error(save(group_col = "nope"), "nope")
  expect_error(save(batch_col = "nope_day"), "nope_day")
  expect_error(save(qc_populations = "B cells"), "B cells")
  expect_error(save(qc_mfi_populations = "Treg"), "Treg")
})

test_that("facs_prep_save keeps extra meta columns in data", {
  d <- local_prep_dirs()
  path <- suppressMessages(facs_prep_save(
    fake_prep_data(meta_files), d$wsp_path, "overview", d$panel_dir,
    group_col = "grp"
  ))
  expect_true(all(c("grp", "day") %in% names(readRDS(path)$data)))
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

test_that("facs_prep_load stops when the report is older than the rds", {
  # e.g. prep saved, then rendering the report failed: the old report stays
  d <- local_prep_dirs()
  prep <- local_loaded_prep(d, two_files)
  Sys.setFileTime(prep$paths$qc, Sys.time() - 600)
  expect_error(facs_prep_load(d$wsp_path, d$panel_dir), "report")
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

# ---- facs_prep_run ----------------------------------------------------------

# Temp layout with the real minimal.wsp fixture as a dated workspace (two
# files of mouse 26-1-17, tissues kidney and percoll-kidney).
local_run_dirs <- function(env = parent.frame()) {
  d <- local_prep_dirs(env)
  file.copy(testthat::test_path("../fixtures/minimal.wsp"), d$wsp_path,
            overwrite = TRUE)
  Sys.setFileTime(d$wsp_path, Sys.time() - 3600)
  d
}

test_that("facs_prep_run(final = FALSE) preps, saves with qc settings and renders", {
  skip_if_not(rmarkdown::pandoc_available())
  d <- local_run_dirs()
  mice <- tibble::tibble(mouse_ID = "26-1-17", group = "A", day = "d1",
                         unused = "x")
  keep_cd45 <- function(data) dplyr::filter(data, population == "CD45+")

  out <- suppressMessages(facs_prep_run(
    final = FALSE, wsp_path = d$wsp_path, panel_dir = d$panel_dir,
    assay = "overview", prepare = keep_cd45, mice = mice,
    group_col = "group", batch_col = "day", expected_mice = "26-1-17"
  ))

  expect_equal(
    normalizePath(out),
    normalizePath(file.path(d$panel_dir, "qc", "overview_25-7_20261101_qc.html"))
  )
  prep <- readRDS(file.path(d$panel_dir, "processed",
                            "overview_25-7_20261101_prep.rds"))
  expect_equal(unique(as.character(prep$data$population)), "CD45+")
  expect_true(all(c("group", "day") %in% names(prep$data)))
  expect_false("unused" %in% names(prep$data))
  expect_equal(prep$info$qc$group_col, "group")
  expect_equal(prep$info$qc$batch_col, "day")
  expect_equal(prep$info$qc$expected_mice, "26-1-17")
})

test_that("facs_prep_run(final = FALSE) reads the condition keyword", {
  d <- local_run_dirs()
  # Absent from the workspace -> read as NA -> facs_prep_save() refuses.
  expect_error(
    suppressWarnings(suppressMessages(facs_prep_run(
      final = FALSE, wsp_path = d$wsp_path, panel_dir = d$panel_dir,
      assay = "overview", condition_keyword = "KW_DOES_NOT_EXIST_XYZ"
    ))),
    "KW_DOES_NOT_EXIST_XYZ"
  )
})

test_that("facs_prep_run(final = TRUE) needs a connection", {
  d <- local_prep_dirs()
  expect_error(
    facs_prep_run(final = TRUE, wsp_path = d$wsp_path,
                  panel_dir = d$panel_dir, assay = "overview"),
    "con"
  )
})

test_that("facs_prep_run rejects a non-logical final", {
  d <- local_prep_dirs()
  expect_error(
    facs_prep_run(final = "yes", wsp_path = d$wsp_path,
                  panel_dir = d$panel_dir, assay = "overview"),
    "final"
  )
})

test_that("facs_prep_run(final = TRUE) exports the saved prep and snapshots", {
  con <- local_test_db()
  db_write_experiment(con, experiment_code = "25-7")
  db_write_subjects(con, tibble::tibble(
    mouse_id = c("25-7-1", "25-7-2"), experiment_code = "25-7"
  ))
  d <- local_prep_dirs()
  suppressMessages(local_loaded_prep(d, two_files))
  snap_root <- file.path(d$data_dir, "z_snapshots")
  ex <- tibble::tibble(mouse_ID = "25-7-2", reason = "clots")

  id <- suppressMessages(facs_prep_run(
    final = TRUE, wsp_path = d$wsp_path, panel_dir = d$panel_dir,
    assay = "overview", exclusions = ex, con = con,
    data_dir = d$data_dir, snapshot_root = snap_root
  ))

  expect_type(id, "integer")
  expect_equal(
    DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM facs_stains")$n, 1
  )
  expect_true(DBI::dbIsValid(con))
  snap <- file.path(snap_root, format(Sys.Date(), "%Y-%m-%d"))
  expect_true(file.exists(file.path(snap, "schema.sql")))
})

test_that("facs_prep_run(final = TRUE) skips the snapshot without snapshot_root", {
  con <- local_test_db()
  db_write_experiment(con, experiment_code = "25-7")
  db_write_subjects(con, tibble::tibble(
    mouse_id = c("25-7-1", "25-7-2"), experiment_code = "25-7"
  ))
  d <- local_prep_dirs()
  suppressMessages(local_loaded_prep(d, two_files))

  suppressMessages(facs_prep_run(
    final = TRUE, wsp_path = d$wsp_path, panel_dir = d$panel_dir,
    assay = "overview", con = con, data_dir = d$data_dir
  ))

  expect_false(dir.exists(file.path(d$data_dir, "z_snapshots")))
})
