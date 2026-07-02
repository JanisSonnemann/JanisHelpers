library(testthat)
library(JanisHelpers)

meta_path <- testthat::test_path("../fixtures/meta_minimal.xlsx")
wsp_path  <- testthat::test_path("../fixtures/minimal.wsp")
meta_skip_msg <- "meta_minimal.xlsx fixture not available"
wsp_skip_msg  <- "minimal.wsp fixture not available"

test_that("meta_read() returns a named list with one tibble per sheet", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta_list <- meta_read(meta_path)
  expect_type(meta_list, "list")
  expect_setequal(names(meta_list), c("meta", "organ_weights", "facs_volumes"))
  expect_s3_class(meta_list$meta, "tbl_df")
  expect_s3_class(meta_list$organ_weights, "tbl_df")
  expect_s3_class(meta_list$facs_volumes, "tbl_df")
})

test_that("meta_read() preserves mouse_ID verbatim in every sheet", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta_list <- meta_read(meta_path)
  for (sheet in meta_list) {
    expect_true("mouse_ID" %in% names(sheet))
    expect_false("mouse_id" %in% names(sheet))
  }
})

test_that("meta_read() coerces date columns to Date in the meta sheet", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta <- meta_read(meta_path)$meta
  date_cols <- intersect(names(meta), c("dob", "start_date", "bmt_date", "death_date"))
  expect_true(length(date_cols) > 0L)
  for (col in date_cols) {
    expect_s3_class(meta[[col]], "Date")
  }
})

test_that("meta_read() coerces group to a factor in the meta sheet", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta <- meta_read(meta_path)$meta
  expect_true("group" %in% names(meta))
  expect_s3_class(meta$group, "factor")
})

test_that("meta_read() trims whitespace and drops empty rows/cols", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta <- meta_read(meta_path)$meta
  chr_cols <- names(meta)[vapply(meta, is.character, logical(1))]
  for (col in chr_cols) {
    vals <- stats::na.omit(meta[[col]])
    # ignore_attr = TRUE: na.omit() attaches an `na.action` attribute that
    # stringr::str_trim() does not propagate; only the values matter here.
    expect_equal(vals, stringr::str_trim(vals), ignore_attr = TRUE)
  }
  expect_true(all(!vapply(meta, function(x) all(is.na(x)), logical(1))))
})

test_that("meta_annotate() left-joins data and meta on the by column", {
  data <- tibble::tibble(mouse_ID = c("A", "B"), value = c(1, 2))
  meta <- tibble::tibble(mouse_ID = c("A", "B"), sex = c("m", "f"))
  result <- meta_annotate(data, meta)
  expect_true("sex" %in% names(result))
  expect_equal(result$sex, c("m", "f"))
})

test_that("meta_annotate() warns and keeps NA for unmatched by values", {
  data <- tibble::tibble(mouse_ID = c("A", "B"), value = c(1, 2))
  meta <- tibble::tibble(mouse_ID = "A", sex = "m")
  expect_warning(
    result <- meta_annotate(data, meta),
    regexp = "B",
    fixed = TRUE
  )
  expect_true(is.na(result$sex[result$mouse_ID == "B"]))
})

test_that("meta_annotate() errors when by column missing from data", {
  data <- tibble::tibble(id = "A", value = 1)
  meta <- tibble::tibble(mouse_ID = "A", sex = "m")
  expect_error(meta_annotate(data, meta), regexp = "data", ignore.case = TRUE)
})

test_that("meta_annotate() errors when by column missing from meta", {
  data <- tibble::tibble(mouse_ID = "A", value = 1)
  meta <- tibble::tibble(id = "A", sex = "m")
  expect_error(meta_annotate(data, meta), regexp = "meta", ignore.case = TRUE)
})

test_that("meta_annotate() errors on colliding non-by column names", {
  data <- tibble::tibble(mouse_ID = "A", sex = "unknown")
  meta <- tibble::tibble(mouse_ID = "A", sex = "m")
  expect_error(meta_annotate(data, meta), regexp = "sex", fixed = TRUE)
})

test_that("facs_read_wsp() data can be annotated end-to-end with meta_read()", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  skip_if_not(file.exists(wsp_path), wsp_skip_msg)

  facs_data <- suppressMessages(
    facs_read_wsp(wsp_path, keywords = "mouse_ID")
  )$data
  meta <- meta_read(meta_path)$meta

  result <- expect_no_warning(meta_annotate(facs_data, meta))

  expect_true(all(c("sex", "group", "cage") %in% names(result)))
  expect_false(any(is.na(result$sex)))
})

test_that("meta_annotate() silently drops by values present only in meta", {
  data <- tibble::tibble(mouse_ID = "A", value = 1)
  meta <- tibble::tibble(mouse_ID = c("A", "B"), sex = c("m", "f"))
  result <- expect_no_warning(meta_annotate(data, meta))
  expect_equal(nrow(result), 1L)
  expect_equal(result$mouse_ID, "A")
  expect_false("B" %in% result$mouse_ID)
})

test_that("meta_annotate() joins on multiple by columns", {
  data <- tibble::tibble(mouse_ID = c("A", "A"), tissue = c("kidney", "lung"), value = c(1, 2))
  meta <- tibble::tibble(mouse_ID = c("A", "A"), tissue = c("kidney", "lung"), weight = c(100, 50))
  result <- meta_annotate(data, meta, by = c("mouse_ID", "tissue"))
  expect_equal(result$weight, c(100, 50))
})

test_that("meta_annotate() warns listing unmatched mouse_ID/tissue combinations", {
  data <- tibble::tibble(mouse_ID = c("A", "A"), tissue = c("kidney", "lung"), value = c(1, 2))
  meta <- tibble::tibble(mouse_ID = "A", tissue = "kidney", weight = 100)
  expect_warning(
    result <- meta_annotate(data, meta, by = c("mouse_ID", "tissue")),
    regexp = "tissue=lung",
    fixed = TRUE
  )
  expect_true(is.na(result$weight[result$tissue == "lung"]))
})

test_that("meta_read() skips the mouse_id rename when no such column exists", {
  testthat::local_mocked_bindings(
    excel_sheets = function(...) "sheet1",
    read_excel   = function(...) tibble::tibble(subject = "A", sex = "m"),
    .package = "readxl"
  )
  meta_list <- meta_read("fake/path.xlsx")
  expect_false("mouse_ID" %in% names(meta_list$sheet1))
  expect_false("mouse_id" %in% names(meta_list$sheet1))
  expect_true("subject" %in% names(meta_list$sheet1))
})

test_that("meta_read() skips group factor coercion when no group column exists", {
  testthat::local_mocked_bindings(
    excel_sheets = function(...) "sheet1",
    read_excel   = function(...) tibble::tibble(mouse_id = "A", sex = "m"),
    .package = "readxl"
  )
  meta_list <- meta_read("fake/path.xlsx")
  expect_false("group" %in% names(meta_list$sheet1))
  expect_true("mouse_ID" %in% names(meta_list$sheet1))
})

test_that("meta_clean() pivots organ_weights long and joins with facs_volumes", {
  skip_if_not(file.exists(meta_path), meta_skip_msg)
  meta_list <- meta_read(meta_path)
  result <- meta_clean(meta_list)

  expect_equal(nrow(result), 9L)
  expect_true(all(c("mouse_ID", "tissue", "total_weight", "facs_weight",
                     "total_vol", "overview_vol", "treg_vol") %in% names(result)))
  expect_setequal(unique(result$tissue), c("kidney", "lung", "spleen"))

  kidney_row <- dplyr::filter(result, mouse_ID == "26-1-1", tissue == "kidney")
  expect_equal(kidney_row$total_weight, 300)
  expect_equal(kidney_row$facs_weight, 200)
})

test_that("meta_clean() errors when organ_weights or facs_volumes is missing", {
  expect_error(
    meta_clean(list(meta = tibble::tibble(mouse_ID = "A"))),
    regexp = "organ_weights",
    fixed = TRUE
  )
})
