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

  result_all <- suppressMessages(facs_read_wsp(wsp_path))
  n_all      <- dplyr::n_distinct(result_all$data$file_name)

  # Find a group that has strictly fewer samples than the total
  strict_group <- NULL
  for (g in groups) {
    res_g <- tryCatch(
      suppressMessages(facs_read_wsp(wsp_path, group = g)),
      error = function(e) NULL
    )
    if (!is.null(res_g) && dplyr::n_distinct(res_g$data$file_name) < n_all) {
      strict_group <- g
      break
    }
  }
  skip_if(is.null(strict_group), "No group with fewer samples than total in fixture")

  result_group <- suppressMessages(facs_read_wsp(wsp_path, group = strict_group))
  expect_lt(dplyr::n_distinct(result_group$data$file_name), n_all)
})

test_that("invalid group name stops with informative error", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  expect_error(
    facs_read_wsp(wsp_path, group = "GROUP_DOES_NOT_EXIST_XYZ"),
    regexp = "not found",
    ignore.case = TRUE
  )
})

test_that("boolean gate populations (NotNode/AndNode/OrNode) are exported", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))

  boolean_pop <- result$data |>
    dplyr::filter(
      file_name == "26-1-17_percoll-kidney_E06.fcs",
      population_full_path == "Singlets/non-debris-"
    )

  expect_true(nrow(boolean_pop) > 0)
  expect_equal(
    boolean_pop$value[boolean_pop$metric == "count"],
    378685
  )
})

test_that("populations gated downstream of a boolean gate are exported", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))

  child_pop <- result$data |>
    dplyr::filter(
      file_name == "26-1-17_percoll-kidney_E06.fcs",
      population_full_path == "Singlets/non-debris-/CD45+"
    )

  expect_true(nrow(child_pop) > 0)
  expect_equal(
    child_pop$value[child_pop$metric == "count"],
    984
  )
})

test_that("population is a factor ordered by depth-first gating hierarchy traversal", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))

  expect_true(is.factor(result$data$population))
  expect_false(is.ordered(result$data$population))

  expected_levels <- unique(basename(unique(result$data$population_full_path)))
  expect_equal(levels(result$data$population), expected_levels)
})

test_that("ancestor populations sort before their descendants in the population factor", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))

  paths <- unique(result$data$population_full_path)
  level_index <- function(full_path) {
    match(basename(full_path), levels(result$data$population))
  }

  for (child_path in paths) {
    parent_path <- dirname(child_path)
    if (parent_path %in% paths && parent_path != child_path) {
      expect_lt(level_index(parent_path), level_index(child_path))
    }
  }
})
