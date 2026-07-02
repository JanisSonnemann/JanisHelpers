library(testthat)
library(JanisHelpers)

meta_path <- testthat::test_path("../fixtures/meta_minimal.xlsx")
wsp_path  <- testthat::test_path("../fixtures/minimal.wsp")
meta_skip_msg <- "meta_minimal.xlsx fixture not available"
wsp_skip_msg  <- "minimal.wsp fixture not available"

test_that("meta_read() returns a tibble with mouse_ID preserved verbatim", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta <- meta_read(meta_path)
  expect_s3_class(meta, "tbl_df")
  expect_true("mouse_ID" %in% names(meta))
  expect_false("mouse_id" %in% names(meta))
})

test_that("meta_read() coerces date columns to Date", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta <- meta_read(meta_path)
  date_cols <- intersect(names(meta), c("dob", "start_date", "bmt_date", "death_date"))
  expect_true(length(date_cols) > 0L)
  for (col in date_cols) {
    expect_s3_class(meta[[col]], "Date")
  }
})

test_that("meta_read() coerces group to a factor", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta <- meta_read(meta_path)
  expect_true("group" %in% names(meta))
  expect_s3_class(meta$group, "factor")
})

test_that("meta_read() trims whitespace and drops empty rows/cols", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta <- meta_read(meta_path)
  chr_cols <- names(meta)[vapply(meta, is.character, logical(1))]
  for (col in chr_cols) {
    vals <- stats::na.omit(meta[[col]])
    # ignore_attr = TRUE: na.omit() attaches an `na.action` attribute that
    # stringr::str_trim() does not propagate; only the values matter here.
    expect_equal(vals, stringr::str_trim(vals), ignore_attr = TRUE)
  }
  expect_true(all(!vapply(meta, function(x) all(is.na(x)), logical(1))))
})
