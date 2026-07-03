library(testthat)
library(JanisHelpers)

wsp_path <- testthat::test_path("../fixtures/Treg.wsp")
fcs_dir  <- testthat::test_path("../fixtures/Treg")
skip_msg <- "Treg fixture not available"
cd45_gate <- "Singlets/Lymphocytes/live/CD45+"

test_that("facs_read_fcs_gated() returns a wide tibble with one row per event", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  result <- facs_read_fcs_gated(
    wsp_path  = wsp_path,
    gate_path = cd45_gate,
    markers   = c("CD4", "CD45")
  )

  expect_s3_class(result, "tbl_df")
  expect_true(all(c("file_name", "CD4", "CD45") %in% names(result)))
  expect_type(result$CD4, "double")
  expect_type(result$CD45, "double")
  expect_equal(dplyr::n_distinct(result$file_name), 6L)
  expect_gt(nrow(result), 1000L)
  # Transformed (biexponential/logicle) scale, not raw fluorescence intensity
  expect_true(max(abs(result$CD4)) < 10000)
})

test_that("facs_read_fcs_gated() attaches requested keyword columns", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  result <- facs_read_fcs_gated(
    wsp_path  = wsp_path,
    gate_path = cd45_gate,
    markers   = c("CD4", "CD45"),
    keywords  = c("mouse_ID", "tissue")
  )

  expect_true(all(c("mouse_ID", "tissue") %in% names(result)))
  expect_false(any(is.na(result$mouse_ID)))
  expect_false(any(is.na(result$tissue)))
  expect_setequal(unique(result$tissue), c("kidney", "lung", "spleen"))
})

test_that("facs_read_fcs_gated() warns and skips every sample when gate_path matches none", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  expect_warning(
    result <- facs_read_fcs_gated(
      wsp_path  = wsp_path,
      gate_path = "Nonexistent/Path",
      markers   = c("CD4")
    ),
    "not found"
  )
  expect_equal(nrow(result), 0L)
})

test_that("facs_read_fcs_gated() errors immediately when a marker cannot be matched", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  expect_error(
    facs_read_fcs_gated(
      wsp_path  = wsp_path,
      gate_path = cd45_gate,
      markers   = c("CD4", "NotARealMarker")
    ),
    "NotARealMarker"
  )
})

test_that("facs_read_fcs_gated() warns and fills NA when a keyword is missing for every sample", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  expect_warning(
    result <- facs_read_fcs_gated(
      wsp_path  = wsp_path,
      gate_path = cd45_gate,
      markers   = c("CD4"),
      keywords  = c("mouse_ID", "not_a_real_keyword")
    ),
    "not_a_real_keyword"
  )
  expect_true(all(is.na(result$not_a_real_keyword)))
  expect_false(any(is.na(result$mouse_ID)))
})
