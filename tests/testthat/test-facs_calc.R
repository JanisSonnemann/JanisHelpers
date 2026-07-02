library(testthat)
library(JanisHelpers)

wsp_path  <- testthat::test_path("../fixtures/minimal.wsp")
wsp_skip_msg <- "minimal.wsp fixture not available"

test_that("facs_calc_pct_of() computes value / ref_count as a 0-1 fraction", {
  data <- tibble::tibble(
    file_name             = c("f1", "f1", "f1"),
    population_full_path  = c("CD45+", "CD45+/CD3+", "CD45+/CD3+/CD4+"),
    population             = c("CD45+", "CD3+", "CD4+"),
    metric                 = c("count", "count", "count"),
    value                  = c(1000, 400, 100)
  )

  result <- facs_calc_pct_of(data, ref_pop = "CD45+")
  new_rows <- dplyr::filter(result, metric == "pct_of_CD45+") |>
    dplyr::arrange(population)

  expect_equal(new_rows$population, c("CD3+", "CD4+"))
  expect_equal(new_rows$value, c(0.4, 0.1))
  expect_false("CD45+" %in% new_rows$population)
})

test_that("facs_calc_pct_of() errors when ref_pop matches more than one row per file_name", {
  data <- tibble::tibble(
    file_name             = c("f1", "f1"),
    population_full_path  = c("A/Live", "B/Live"),
    population             = c("Live", "Live"),
    metric                 = c("count", "count"),
    value                  = c(100, 200)
  )

  expect_error(facs_calc_pct_of(data, ref_pop = "Live"), regexp = "f1", fixed = TRUE)
})

test_that("facs_calc_pct_of() warns and fills NA when ref_pop has no match for a file_name", {
  data <- tibble::tibble(
    file_name             = c("f1", "f1", "f2", "f2"),
    population_full_path  = c("CD45+", "CD45+/CD3+", "CD3+", "CD3+/CD4+"),
    population             = c("CD45+", "CD3+", "CD3+", "CD4+"),
    metric                 = c("count", "count", "count", "count"),
    value                  = c(1000, 300, 50, 20)
  )

  expect_warning(
    result <- facs_calc_pct_of(data, ref_pop = "CD45+"),
    regexp = "f2",
    fixed = TRUE
  )

  new_rows <- dplyr::filter(result, metric == "pct_of_CD45+")
  f1_row <- dplyr::filter(new_rows, file_name == "f1", population == "CD3+")
  f2_rows <- dplyr::filter(new_rows, file_name == "f2")

  expect_equal(f1_row$value, 0.3)
  expect_true(all(is.na(f2_rows$value)))
})
