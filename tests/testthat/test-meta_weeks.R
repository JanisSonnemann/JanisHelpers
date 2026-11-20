wk_meta <- tibble::tibble(
  experiment_code = "26-1",
  mouse_id = c("m1", "m2", "m3"),
  bmt_date = as.Date("2026-05-12"),
  end_date = as.Date(c("2026-06-10", "2026-06-15", "2026-06-16"))
)

wk_data <- tidyr::expand_grid(
  experiment_code = "26-1",
  mouse_id = c("m1", "m2", "m3"),
  collected_at = as.Date(c("2026-05-05", "2026-05-19", "2026-06-02"))
) |>
  dplyr::bind_rows(dplyr::select(wk_meta, experiment_code, mouse_id, collected_at = end_date)) |>
  dplyr::mutate(value = seq_len(dplyr::n()))

test_that("26-1 dates give weeks -1, 1, 3, 4, 5", {
  out <- meta_label_weeks(wk_data, wk_meta)
  expect_equal(out$week_num[out$collected_at == as.Date("2026-05-05")], rep(-1L, 3))
  expect_equal(unique(out$week_num[out$collected_at == as.Date("2026-06-10")]), 4L)
  expect_equal(
    unique(out$week_num[out$collected_at %in% as.Date(c("2026-06-15", "2026-06-16"))]),
    5L
  )
  expect_equal(levels(out$week), c("week -1", "week 1", "week 3", "week 4", "week 5"))
  expect_true(is.ordered(out$week))
  expect_type(out$days_post_bmt, "integer")
  expect_equal(out$days_post_bmt[out$collected_at == as.Date("2026-06-16")], 35L)
})

test_that("end-points in one calendar week straddling a rounding boundary pool", {
  # BMT 2026-05-12; ends day 31 (Fri 06-12) and day 32 (Sat 06-13), same ISO week:
  # round(31/7) = 4, round(32/7) = 5 -> pooled round(31.5/7) = 5 for both
  meta <- dplyr::mutate(wk_meta[1:2, ], end_date = as.Date(c("2026-06-12", "2026-06-13")))
  data <- dplyr::transmute(meta, experiment_code, mouse_id, collected_at = end_date)
  out <- meta_label_weeks(data, meta)
  expect_equal(out$week_num, c(5L, 5L))
})

test_that("serial samples are not pooled", {
  data <- tibble::tibble(
    experiment_code = "26-1", mouse_id = c("m2", "m3"),
    collected_at = as.Date(c("2026-06-12", "2026-06-13"))
  )
  out <- meta_label_weeks(data, wk_meta)
  expect_equal(out$week_num, c(4L, 5L))
})

test_that("mice missing from meta error, naming them", {
  data <- dplyr::bind_rows(
    wk_data,
    tibble::tibble(experiment_code = "26-1", mouse_id = "m9",
                   collected_at = as.Date("2026-05-19"))
  )
  expect_error(meta_label_weeks(data, wk_meta), "m9")
})

test_that("missing bmt_date gives NA week with a warning naming the mouse", {
  meta <- dplyr::mutate(wk_meta, bmt_date = dplyr::if_else(mouse_id == "m1", as.Date(NA), bmt_date))
  expect_warning(out <- meta_label_weeks(wk_data, meta), "m1")
  expect_true(all(is.na(out$week_num[out$mouse_id == "m1"])))
  expect_false(anyNA(out$week_num[out$mouse_id != "m1"]))
})

test_that("existing output columns error", {
  expect_error(meta_label_weeks(dplyr::mutate(wk_data, week = 1), wk_meta),
               "already has column.*week")
})

test_that("POSIXct dates are accepted", {
  meta <- dplyr::mutate(wk_meta, bmt_date = as.POSIXct(bmt_date), end_date = as.POSIXct(end_date))
  out <- meta_label_weeks(wk_data, meta)
  expect_equal(sort(unique(out$week_num)), c(-1L, 1L, 3L, 4L, 5L))
})

test_that("missing end_date warns, naming the mouse", {
  meta <- dplyr::mutate(wk_meta, end_date = dplyr::if_else(mouse_id == "m2", as.Date(NA), end_date))
  expect_warning(meta_label_weeks(wk_data, meta), "end_date.*m2")
})
