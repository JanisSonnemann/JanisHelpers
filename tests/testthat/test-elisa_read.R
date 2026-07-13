library(testthat)
library(JanisHelpers)

il17a_path <- testthat::test_path("../fixtures/Results_IL-17A.xlsx")
tnfa_path  <- testthat::test_path("../fixtures/Results_TNFa.xlsx")
il17a_skip_msg <- "Results_IL-17A.xlsx fixture not available"
tnfa_skip_msg  <- "Results_TNFa.xlsx fixture not available"

test_that("elisa_read_results() returns the expected columns and types", {
  skip_if_not(file.exists(il17a_path), il17a_skip_msg)
  dat <- elisa_read_results(il17a_path)
  expect_s3_class(dat, "tbl_df")
  expect_equal(
    names(dat),
    c("cytokine", "sample_id", "replicate", "value", "unit", "result_status")
  )
  expect_type(dat$cytokine, "character")
  expect_type(dat$sample_id, "character")
  expect_type(dat$replicate, "integer")
  expect_type(dat$value, "double")
  expect_type(dat$unit, "character")
  expect_type(dat$result_status, "character")
  expect_equal(nrow(dat), 56L)
})

test_that("elisa_read_results() derives cytokine from filename (no cytokine-prefixed sheet)", {
  skip_if_not(file.exists(il17a_path), il17a_skip_msg)
  dat <- elisa_read_results(il17a_path)
  expect_equal(unique(dat$cytokine), "IL-17A")
})

test_that("elisa_read_results() derives cytokine from filename (cytokine-prefixed sheet)", {
  skip_if_not(file.exists(tnfa_path), tnfa_skip_msg)
  dat <- elisa_read_results(tnfa_path)
  expect_equal(unique(dat$cytokine), "TNFa")
  expect_equal(nrow(dat), 56L)
})

test_that("elisa_read_results() parses the unit from the Backcalc column header", {
  skip_if_not(file.exists(il17a_path), il17a_skip_msg)
  dat <- elisa_read_results(il17a_path)
  expect_equal(unique(dat$unit), "pg/ml")
})

test_that("elisa_read_results() sets value to NA for OOR rows and keeps result_status", {
  skip_if_not(file.exists(il17a_path), il17a_skip_msg)
  dat <- elisa_read_results(il17a_path)
  oor <- dplyr::filter(dat, result_status == "OOR<")
  expect_equal(nrow(oor), 20L)
  expect_true(all(is.na(oor$value)))
  expect_setequal(
    unique(oor$sample_id),
    c("25-7-5", "25-7-11", "25-7-12", "25-7-13", "25-7-14",
      "25-7-20", "25-7-21", "25-7-24", "25-7-25", "25-7-28")
  )
})

test_that("elisa_read_results() keeps a non-NA value for <LLOQ rows (below quantification, not out of range)", {
  skip_if_not(file.exists(tnfa_path), tnfa_skip_msg)
  dat <- elisa_read_results(tnfa_path)
  lloq <- dplyr::filter(dat, result_status == "<LLOQ")
  expect_equal(nrow(lloq), 11L)
  expect_true(all(!is.na(lloq$value)))
  row <- dplyr::filter(lloq, sample_id == "25-7-3", replicate == 1L)
  expect_equal(row$value, 1.303633, tolerance = 1e-6)
})

test_that("elisa_read_results() uses the cytokine argument when supplied, overriding the filename", {
  skip_if_not(file.exists(tnfa_path), tnfa_skip_msg)
  dat <- elisa_read_results(tnfa_path, cytokine = "Custom")
  expect_equal(unique(dat$cytokine), "Custom")
})
