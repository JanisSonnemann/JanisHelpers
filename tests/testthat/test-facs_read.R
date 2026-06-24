library(testthat)
library(JanisHelpers)

wsp_path <- testthat::test_path("../fixtures/minimal.wsp")
skip_msg  <- "WSP fixture not available"

test_that("facs_read_wsp() returns a named list with slots data, meta, panel", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))
  expect_type(result, "list")
  expect_named(result, c("data", "meta", "panel"))
})

test_that("data slot is a tibble with required columns", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))
  expect_s3_class(result$data, "tbl_df")
  expect_true(all(
    c("file_name", "population_full_path", "population", "metric", "value")
    %in% names(result$data)
  ))
})

test_that("data metric contains count and fraction_of_parent", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))
  expect_true("count" %in% result$data$metric)
  expect_true("fraction_of_parent" %in% result$data$metric)
})

test_that("data value column is numeric", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))
  expect_type(result$data$value, "double")
})

test_that("meta slot is a wide tibble with one row per file", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))
  expect_s3_class(result$meta, "tbl_df")
  n_files <- dplyr::n_distinct(result$data$file_name)
  expect_equal(nrow(result$meta), n_files)
  expect_true("file_name" %in% names(result$meta))
})

test_that("panel slot is a wide tibble with one row per file", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))
  expect_s3_class(result$panel, "tbl_df")
  n_files <- dplyr::n_distinct(result$data$file_name)
  expect_equal(nrow(result$panel), n_files)
  expect_true("file_name" %in% names(result$panel))
})

test_that("requested keywords are joined into data", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  doc <- xml2::read_xml(wsp_path)
  kw_nodes <- xml2::xml_find_all(doc, ".//SampleList/Sample[1]//Keywords/Keyword")
  all_kws  <- xml2::xml_attr(kw_nodes, "name")
  user_kw  <- all_kws[!grepl("^\\$", all_kws)][1]
  skip_if(is.na(user_kw), "No user-level keywords in fixture")

  result <- suppressMessages(facs_read_wsp(wsp_path, keywords = user_kw))
  expect_true(user_kw %in% names(result$data))
})

test_that("missing keyword warns and fills NA", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  expect_warning(
    result <- suppressMessages(
      facs_read_wsp(wsp_path, keywords = "KW_DOES_NOT_EXIST_XYZ")
    ),
    regexp = "not found",
    ignore.case = TRUE
  )
  expect_true("KW_DOES_NOT_EXIST_XYZ" %in% names(result$data))
  expect_true(all(is.na(result$data$KW_DOES_NOT_EXIST_XYZ)))
})

test_that("group filtering restricts to that group's samples", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  doc    <- xml2::read_xml(wsp_path)
  groups <- xml2::xml_attr(xml2::xml_find_all(doc, ".//Groups/GroupNode"), "name")
  skip_if(length(groups) == 0L, "No groups in fixture")

  result_all   <- suppressMessages(facs_read_wsp(wsp_path))
  result_group <- suppressMessages(facs_read_wsp(wsp_path, group = groups[[1]]))
  expect_lte(
    dplyr::n_distinct(result_group$data$file_name),
    dplyr::n_distinct(result_all$data$file_name)
  )
})

test_that("invalid group name stops with informative error", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  expect_error(
    facs_read_wsp(wsp_path, group = "GROUP_DOES_NOT_EXIST_XYZ"),
    regexp = "not found",
    ignore.case = TRUE
  )
})
