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

  warns <- testthat::capture_warnings(
    result <- facs_read_fcs_gated(
      wsp_path  = wsp_path,
      gate_path = "Nonexistent/Path",
      markers   = c("CD4")
    )
  )
  expect_true(all(grepl("not found", warns)))
  expect_length(warns, 6L)
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

test_that("facs_read_fcs_gated() matches a marker by raw channel name when it has no stain label", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  result <- facs_read_fcs_gated(
    wsp_path  = wsp_path,
    gate_path = cd45_gate,
    markers   = c("CD4", "FSC-A")
  )

  expect_true("FSC-A" %in% names(result))
  expect_type(result[["FSC-A"]], "double")
})

test_that("facs_read_fcs_gated() downsamples to max_events per sample, reproducibly with seed", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  result_a <- facs_read_fcs_gated(
    wsp_path   = wsp_path,
    gate_path  = cd45_gate,
    markers    = c("CD4", "CD45"),
    max_events = 500,
    seed       = 42
  )
  result_b <- facs_read_fcs_gated(
    wsp_path   = wsp_path,
    gate_path  = cd45_gate,
    markers    = c("CD4", "CD45"),
    max_events = 500,
    seed       = 42
  )

  expect_equal(nrow(result_a), 6L * 500L)
  expect_identical(result_a, result_b)
})

test_that("facs_read_fcs_gated() auto-derives fcs_dir from the wsp filename", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  result_auto <- facs_read_fcs_gated(
    wsp_path  = wsp_path,
    gate_path = cd45_gate,
    markers   = c("CD4")
  )
  result_explicit <- facs_read_fcs_gated(
    wsp_path  = wsp_path,
    gate_path = cd45_gate,
    markers   = c("CD4"),
    fcs_dir   = fcs_dir
  )

  expect_identical(result_auto, result_explicit)
})

test_that("facs_read_fcs_gated() loads only the samples in a non-default FlowJo group", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  result <- facs_read_fcs_gated(
    wsp_path  = wsp_path,
    gate_path = cd45_gate,
    markers   = c("CD4"),
    group     = "kidney"
  )

  expect_equal(dplyr::n_distinct(result$file_name), 2L)
  expect_true(all(grepl("kidney", result$file_name)))
})

test_that("facs_read_fcs_gated() reads keywords from the .wsp XML, including ones absent from the raw .fcs file", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  doc <- xml2::read_xml(wsp_path)
  kw_nodes <- xml2::xml_find_all(doc, ".//SampleList/Sample//Keywords")
  for (kn in kw_nodes) {
    xml2::xml_add_child(kn, "Keyword", name = "workspace_only_key", value = "hello")
  }
  tmp_wsp <- tempfile(fileext = ".wsp")
  xml2::write_xml(doc, tmp_wsp)
  on.exit(unlink(tmp_wsp), add = TRUE)

  result <- facs_read_fcs_gated(
    wsp_path  = tmp_wsp,
    gate_path = cd45_gate,
    markers   = c("CD4"),
    keywords  = c("workspace_only_key"),
    fcs_dir   = fcs_dir
  )

  expect_false(any(is.na(result$workspace_only_key)))
  expect_true(all(result$workspace_only_key == "hello"))
})
