files <- tibble::tibble(
  file_name = c("m1_spl.fcs", "m1_lung.fcs", "m2_spl.fcs"),
  mouse_ID = c("25-7-1", "25-7-1", "25-7-2"),
  tissue = c("spleen", "lung", "spleen"),
  condition = "none"
)

test_that("whole-mouse row excludes all its files", {
  ex <- tibble::tibble(mouse_ID = "25-7-1", reason = "died early")
  m <- facs_match_exclusions(files, ex, panel = "overview")
  expect_setequal(m$excluded_files$file_name, c("m1_spl.fcs", "m1_lung.fcs"))
  expect_equal(m$applicable$n_files_matched, 2L)
  expect_equal(nrow(m$unmatched), 0)
  expect_named(
    m$applicable,
    c("mouse_ID", "tissue", "panel", "condition", "file_name", "reason",
      "date", "n_files_matched")
  )
})

test_that("tissue and file_name narrow the match", {
  ex <- tibble::tibble(
    mouse_ID = c("25-7-1", "25-7-2"),
    tissue = c("lung", NA),
    file_name = c(NA, "m2_spl.fcs"),
    reason = c("bad stain", "clog")
  )
  m <- facs_match_exclusions(files, ex, panel = "overview")
  expect_setequal(m$excluded_files$file_name, c("m1_lung.fcs", "m2_spl.fcs"))
  expect_equal(m$applicable$n_files_matched, c(1L, 1L))
})

test_that("condition narrows the match", {
  f <- dplyr::mutate(files, condition = c("MPO", "unstim", "MPO"))
  ex <- tibble::tibble(mouse_ID = "25-7-1", condition = "MPO", reason = "x")
  m <- facs_match_exclusions(f, ex, panel = "restim")
  expect_equal(m$excluded_files$file_name, "m1_spl.fcs")
})

test_that("rows for another panel are ignored, empty panel applies", {
  ex <- tibble::tibble(
    mouse_ID = c("25-7-2", "25-7-1"),
    panel = c("treg", NA),
    reason = c("treg only", "all panels")
  )
  m <- facs_match_exclusions(files, ex, panel = "overview")
  expect_equal(nrow(m$applicable), 1)
  expect_equal(m$applicable$reason, "all panels")
  expect_setequal(m$excluded_files$file_name, c("m1_spl.fcs", "m1_lung.fcs"))
})

test_that("several rows hitting one file join their reasons", {
  ex <- tibble::tibble(
    mouse_ID = "25-7-2", tissue = c("spleen", NA), reason = c("clog", "sick")
  )
  m <- facs_match_exclusions(files, ex, panel = "overview")
  expect_equal(m$excluded_files$file_name, "m2_spl.fcs")
  expect_equal(m$excluded_files$reason, "clog; sick")
})

test_that("unmatched applicable rows are returned, not an error", {
  ex <- tibble::tibble(mouse_ID = "25-7-1", tissue = "Spleen", reason = "typo")
  m <- facs_match_exclusions(files, ex, panel = "overview")
  expect_equal(nrow(m$unmatched), 1)
  expect_equal(m$unmatched$tissue, "Spleen")
  expect_equal(nrow(m$excluded_files), 0)
})

test_that("logical NA columns, Date column and empty sheets work", {
  ex <- tibble::tibble(
    mouse_ID = "25-7-2", tissue = NA, panel = NA, reason = "clog",
    date = as.Date("2026-11-01")
  )
  m <- facs_match_exclusions(files, ex, panel = "overview")
  expect_equal(m$excluded_files$file_name, "m2_spl.fcs")
  expect_equal(m$applicable$date, as.Date("2026-11-01"))
  expect_type(m$applicable$tissue, "character")

  for (empty in list(NULL, tibble::tibble())) {
    m <- facs_match_exclusions(files, empty, panel = "overview")
    expect_equal(nrow(m$applicable), 0)
    expect_equal(nrow(m$excluded_files), 0)
    expect_equal(nrow(m$unmatched), 0)
    expect_named(m$excluded_files, c("file_name", "reason"))
  }
})

test_that("missing mouse_ID/reason column or empty reason errors", {
  expect_error(
    facs_match_exclusions(files, tibble::tibble(mouse_ID = "25-7-1"), "overview"),
    "reason"
  )
  expect_error(
    facs_match_exclusions(files, tibble::tibble(reason = "x"), "overview"),
    "mouse_ID"
  )
  expect_error(
    facs_match_exclusions(
      files, tibble::tibble(mouse_ID = "25-7-1", reason = NA), "overview"
    ),
    "reason"
  )
})
