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

# ── facs_calc_count_per_g ────────────────────────────────────────────────────

count_per_g_data <- tibble::tibble(
  file_name             = "s1",
  population_full_path  = "CD45+",
  population             = "CD45+",
  metric                 = "count",
  value                  = 1000,
  mouse_ID               = "m1",
  tissue                  = "kidney"
)

count_per_g_meta <- tibble::tibble(
  mouse_ID           = "m1",
  vol_total          = 1000,
  vol_stained        = 100,
  vol_resuspended    = 500,
  vol_measured       = 50,
  organ_piece_weight = 200
)

test_that("facs_calc_count_per_g() HTS formula matches manual calculation (columns)", {
  result <- facs_calc_count_per_g(
    count_per_g_data, count_per_g_meta, tissue = "kidney",
    vol_total = "vol_total", vol_stained = "vol_stained",
    vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
    organ_piece_weight = "organ_piece_weight"
  )
  new_row <- dplyr::filter(result, metric == "count_per_g")
  expect_equal(new_row$value, 500000)
})

test_that("facs_calc_count_per_g() accepts numeric constants in place of column names", {
  meta_no_total <- dplyr::select(count_per_g_meta, !vol_total)

  result <- facs_calc_count_per_g(
    count_per_g_data, meta_no_total, tissue = "kidney",
    vol_total = 1000, vol_stained = "vol_stained",
    vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
    organ_piece_weight = "organ_piece_weight"
  )
  new_row <- dplyr::filter(result, metric == "count_per_g")
  expect_equal(new_row$value, 500000)
})

test_that("facs_calc_count_per_g() only processes rows matching the tissue argument", {
  data_multi <- dplyr::bind_rows(
    count_per_g_data,
    tibble::tibble(
      file_name = "s2", population_full_path = "CD45+", population = "CD45+",
      metric = "count", value = 2000, mouse_ID = "m1", tissue = "spleen"
    )
  )

  result <- facs_calc_count_per_g(
    data_multi, count_per_g_meta, tissue = "kidney",
    vol_total = "vol_total", vol_stained = "vol_stained",
    vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
    organ_piece_weight = "organ_piece_weight"
  )
  new_rows <- dplyr::filter(result, metric == "count_per_g")
  expect_equal(nrow(new_rows), 1L)
  expect_equal(new_rows$file_name, "s1")
})

test_that("facs_calc_count_per_g() errors when a volume argument is neither a column nor a numeric constant", {
  expect_error(
    facs_calc_count_per_g(
      count_per_g_data, count_per_g_meta, tissue = "kidney",
      vol_total = TRUE, vol_stained = "vol_stained",
      vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
      organ_piece_weight = "organ_piece_weight"
    ),
    regexp = "vol_total",
    fixed = TRUE
  )
})

# ── bead path + method_col resolution ────────────────────────────────────────

bead_row <- tibble::tibble(
  file_name = "s1", population_full_path = "beads", population = "beads",
  metric = "count", value = 5200, mouse_ID = "m1", tissue = "kidney"
)

test_that("facs_calc_count_per_g() bead formula matches manual calculation (method_col in meta)", {
  data_beads <- dplyr::bind_rows(count_per_g_data, bead_row)
  meta_beads <- dplyr::mutate(count_per_g_meta, method = "beads")

  result <- facs_calc_count_per_g(
    data_beads, meta_beads, tissue = "kidney",
    vol_total = "vol_total", vol_stained = "vol_stained",
    vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
    organ_piece_weight = "organ_piece_weight",
    method_col = "method"
  )
  new_row <- dplyr::filter(result, metric == "count_per_g", population == "CD45+")
  expect_equal(new_row$value, 100000)
})

test_that("facs_calc_count_per_g() resolves method_col from data (per-sample) over meta", {
  data_beads_kw <- dplyr::bind_rows(count_per_g_data, bead_row) |>
    dplyr::mutate(method = "beads")

  result <- facs_calc_count_per_g(
    data_beads_kw, count_per_g_meta, tissue = "kidney",
    vol_total = "vol_total", vol_stained = "vol_stained",
    vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
    organ_piece_weight = "organ_piece_weight",
    method_col = "method"
  )
  new_row <- dplyr::filter(result, metric == "count_per_g", population == "CD45+")
  expect_equal(new_row$value, 100000)
})

test_that("facs_calc_count_per_g() defaults NA method_col values to hts with no warning", {
  meta_na_method <- dplyr::mutate(count_per_g_meta, method = NA_character_)

  expect_no_warning(
    result <- facs_calc_count_per_g(
      count_per_g_data, meta_na_method, tissue = "kidney",
      vol_total = "vol_total", vol_stained = "vol_stained",
      vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
      organ_piece_weight = "organ_piece_weight",
      method_col = "method"
    )
  )
  new_row <- dplyr::filter(result, metric == "count_per_g")
  expect_equal(new_row$value, 500000)
})

test_that("facs_calc_count_per_g() errors on an invalid method_col value", {
  meta_bad_method <- dplyr::mutate(count_per_g_meta, method = "unknown")

  expect_error(
    facs_calc_count_per_g(
      count_per_g_data, meta_bad_method, tissue = "kidney",
      vol_total = "vol_total", vol_stained = "vol_stained",
      vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
      organ_piece_weight = "organ_piece_weight",
      method_col = "method"
    ),
    regexp = "unknown",
    fixed = TRUE
  )
})

test_that("facs_calc_count_per_g() errors when method_col is not found in data or meta", {
  expect_error(
    facs_calc_count_per_g(
      count_per_g_data, count_per_g_meta, tissue = "kidney",
      vol_total = "vol_total", vol_stained = "vol_stained",
      vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
      organ_piece_weight = "organ_piece_weight",
      method_col = "nonexistent"
    ),
    regexp = "nonexistent",
    fixed = TRUE
  )
})

test_that("facs_calc_count_per_g() warns and fills NA when bead method resolved but no bead count found", {
  meta_beads_missing <- dplyr::mutate(count_per_g_meta, method = "beads")

  expect_warning(
    result <- facs_calc_count_per_g(
      count_per_g_data, meta_beads_missing, tissue = "kidney",
      vol_total = "vol_total", vol_stained = "vol_stained",
      vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
      organ_piece_weight = "organ_piece_weight",
      method_col = "method"
    ),
    regexp = "s1",
    fixed = TRUE
  )
  new_row <- dplyr::filter(result, metric == "count_per_g")
  expect_true(is.na(new_row$value))
})
