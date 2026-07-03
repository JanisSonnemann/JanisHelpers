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
  tissue              = "kidney",
  vol_total          = 1000,
  vol_stained        = 100,
  vol_resuspended    = 500,
  vol_measured       = 50,
  organ_piece_weight = 200
)

test_that("facs_calc_count_per_g() HTS formula matches manual calculation (columns)", {
  result <- facs_calc_count_per_g(
    count_per_g_data, count_per_g_meta,
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
    count_per_g_data, meta_no_total,
    vol_total = 1000, vol_stained = "vol_stained",
    vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
    organ_piece_weight = "organ_piece_weight"
  )
  new_row <- dplyr::filter(result, metric == "count_per_g")
  expect_equal(new_row$value, 500000)
})

test_that("facs_calc_count_per_g() computes count_per_g for every tissue present in data", {
  data_multi <- dplyr::bind_rows(
    count_per_g_data,
    tibble::tibble(
      file_name = "s2", population_full_path = "CD45+", population = "CD45+",
      metric = "count", value = 2000, mouse_ID = "m1", tissue = "spleen"
    )
  )
  meta_multi <- dplyr::bind_rows(
    count_per_g_meta,
    tibble::tibble(
      mouse_ID = "m1", tissue = "spleen",
      vol_total = 1000, vol_stained = 100, vol_resuspended = 500,
      vol_measured = 50, organ_piece_weight = 200
    )
  )

  result <- facs_calc_count_per_g(
    data_multi, meta_multi,
    vol_total = "vol_total", vol_stained = "vol_stained",
    vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
    organ_piece_weight = "organ_piece_weight"
  )
  new_rows <- dplyr::filter(result, metric == "count_per_g") |> dplyr::arrange(file_name)
  expect_equal(nrow(new_rows), 2L)
  expect_equal(new_rows$file_name, c("s1", "s2"))
  expect_equal(new_rows$value, c(500000, 1000000))
})

test_that("facs_calc_count_per_g() warns and fills NA when a mouse_ID/tissue combination has no match in meta", {
  meta_wrong_mouse <- dplyr::mutate(count_per_g_meta, mouse_ID = "m2")

  expect_warning(
    result <- facs_calc_count_per_g(
      count_per_g_data, meta_wrong_mouse,
      vol_total = "vol_total", vol_stained = "vol_stained",
      vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
      organ_piece_weight = "organ_piece_weight"
    ),
    regexp = "m1",
    fixed = TRUE
  )
  new_row <- dplyr::filter(result, metric == "count_per_g")
  expect_true(is.na(new_row$value))
})

test_that("facs_calc_count_per_g() errors when a volume argument is neither a column nor a numeric constant", {
  expect_error(
    facs_calc_count_per_g(
      count_per_g_data, count_per_g_meta,
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
    data_beads, meta_beads,
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
  meta_conflicting <- dplyr::mutate(count_per_g_meta, method = "hts")

  result <- facs_calc_count_per_g(
    data_beads_kw, meta_conflicting,
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
      count_per_g_data, meta_na_method,
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
      count_per_g_data, meta_bad_method,
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
      count_per_g_data, count_per_g_meta,
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
      count_per_g_data, meta_beads_missing,
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

# ── end-to-end fixture test ─────────────────────────────────────────────────

test_that("facs_calc_pct_of() and facs_calc_count_per_g() work end-to-end on real fixture data", {
  skip_if_not(file.exists(wsp_path), wsp_skip_msg)

  facs_data <- suppressMessages(
    facs_read_wsp(wsp_path, keywords = c("mouse_ID", "tissue"))
  )$data

  pct_result <- facs_calc_pct_of(facs_data, ref_pop = "Singlets")
  pct_rows <- dplyr::filter(pct_result, metric == "pct_of_Singlets")
  expect_true(nrow(pct_rows) > 0L)
  expect_true(all(pct_rows$value >= 0 & pct_rows$value <= 1, na.rm = TRUE))

  # minimal.wsp has one mouse (26-1-17) with two tissue values: "kidney" and
  # "percoll-kidney" (one FCS file each) -- covering both in `meta` exercises
  # facs_calc_count_per_g() processing multiple tissues in a single call.
  meta <- tibble::tibble(
    mouse_ID           = "26-1-17",
    tissue              = c("kidney", "percoll-kidney"),
    vol_total          = 1000,
    vol_stained        = 100,
    vol_resuspended    = 500,
    vol_measured       = 50,
    organ_piece_weight = 200
  )

  count_result <- expect_no_warning(facs_calc_count_per_g(
    facs_data, meta,
    vol_total = "vol_total", vol_stained = "vol_stained",
    vol_resuspended = "vol_resuspended", vol_measured = "vol_measured",
    organ_piece_weight = "organ_piece_weight"
  ))
  count_rows <- dplyr::filter(count_result, metric == "count_per_g")
  expect_equal(dplyr::n_distinct(count_rows$file_name), 2L)
  expect_setequal(
    count_rows$file_name,
    c("26-1-17_whole_kidney_E05.fcs", "26-1-17_percoll-kidney_E06.fcs")
  )
})

# ── facs_calc_log2fc / facs_calc_diff (shared calc_restim_proportions_) ─────

log2fc_data <- tibble::tibble(
  file_name             = c("f1", "f1", "f2", "f2"),
  population_full_path  = c("CD3+", "CD3+/CD4+", "CD3+", "CD3+/CD4+"),
  population             = c("CD3+", "CD4+", "CD3+", "CD4+"),
  metric                 = "count",
  value                  = c(1000, 100, 1000, 400),
  mouse_ID               = "m1",
  tissue                  = "spleen",
  restimulation           = c("unstim", "unstim", "MPO", "MPO")
)

test_that("facs_calc_log2fc() computes log2 fold-change of a stim condition vs unstim", {
  result <- facs_calc_log2fc(log2fc_data, ref_pop = "CD3+")

  expect_equal(nrow(result), 1L)
  expect_equal(result$restimulation, "MPO")
  expect_equal(result$population, "CD4+")

  unstim_prop <- (100 + 0.5) / 1000
  mpo_prop    <- (400 + 0.5) / 1000
  expect_equal(result$log2fc, log2(mpo_prop / unstim_prop))
})

test_that("facs_calc_log2fc() excludes ref_pop and ref_level from the output", {
  result <- facs_calc_log2fc(log2fc_data, ref_pop = "CD3+")

  expect_false("CD3+" %in% result$population)
  expect_false("unstim" %in% result$restimulation)
})

test_that("facs_calc_log2fc() pseudocount avoids -Inf when a condition has zero events", {
  data_zero <- dplyr::mutate(
    log2fc_data,
    value = dplyr::if_else(restimulation == "MPO" & population == "CD4+", 0, value)
  )

  result <- facs_calc_log2fc(data_zero, ref_pop = "CD3+")
  expect_true(is.finite(result$log2fc))
})

test_that("facs_calc_log2fc() errors when ref_pop matches more than one row per mouse_ID/tissue/restim-level combo", {
  data_dup <- dplyr::bind_rows(
    log2fc_data,
    tibble::tibble(
      file_name = "f1", population_full_path = "OtherPath/CD3+", population = "CD3+",
      metric = "count", value = 500, mouse_ID = "m1", tissue = "spleen", restimulation = "unstim"
    )
  )

  expect_error(facs_calc_log2fc(data_dup, ref_pop = "CD3+"), regexp = "m1", fixed = TRUE)
})

test_that("facs_calc_log2fc() warns and fills NA when ref_pop has no match for a combo", {
  data_missing_refpop <- dplyr::filter(log2fc_data, !(restimulation == "MPO" & population == "CD3+"))

  expect_warning(
    result <- facs_calc_log2fc(data_missing_refpop, ref_pop = "CD3+"),
    regexp = "m1",
    fixed = TRUE
  )
  expect_true(is.na(result$log2fc))
})

test_that("facs_calc_log2fc() warns and fills NA when ref_level is missing for a mouse_ID/tissue group", {
  data_no_unstim <- dplyr::filter(log2fc_data, restimulation != "unstim")

  expect_warning(
    result <- facs_calc_log2fc(data_no_unstim, ref_pop = "CD3+"),
    regexp = "m1",
    fixed = TRUE
  )
  expect_true(is.na(result$log2fc))
})

test_that("facs_calc_log2fc() errors when restim_col is not a column in data", {
  data_no_restim <- dplyr::select(log2fc_data, !restimulation)

  expect_error(facs_calc_log2fc(data_no_restim, ref_pop = "CD3+"), regexp = "restim_col", fixed = TRUE)
})

test_that("facs_calc_log2fc() handles more than one non-reference restimulation level", {
  data_multi_stim <- dplyr::bind_rows(
    log2fc_data,
    tibble::tibble(
      file_name = "f3", population_full_path = c("CD3+", "CD3+/CD4+"), population = c("CD3+", "CD4+"),
      metric = "count", value = c(1000, 250), mouse_ID = "m1", tissue = "spleen", restimulation = "PMA-Iono"
    )
  )

  result <- facs_calc_log2fc(data_multi_stim, ref_pop = "CD3+")

  expect_setequal(result$restimulation, c("MPO", "PMA-Iono"))
  expect_equal(nrow(result), 2L)
})

test_that("facs_calc_log2fc() carries through a passthrough column constant within mouse_ID/tissue", {
  data_group <- dplyr::mutate(log2fc_data, group = "WT")

  result <- facs_calc_log2fc(data_group, ref_pop = "CD3+")
  expect_true("group" %in% names(result))
  expect_equal(unique(result$group), "WT")
})

test_that("facs_calc_log2fc() drops a passthrough-candidate column that varies within mouse_ID/tissue", {
  data_varying <- dplyr::mutate(
    log2fc_data,
    batch = dplyr::if_else(restimulation == "unstim", "batch1", "batch2")
  )

  result <- facs_calc_log2fc(data_varying, ref_pop = "CD3+")
  expect_false("batch" %in% names(result))
})

test_that("facs_calc_log2fc() works with non-default restim_col/ref_level names", {
  data_custom <- log2fc_data |>
    dplyr::rename(condition = restimulation) |>
    dplyr::mutate(condition = dplyr::if_else(condition == "unstim", "baseline", "MPO"))

  result <- facs_calc_log2fc(data_custom, ref_pop = "CD3+", restim_col = "condition", ref_level = "baseline")
  expect_equal(result$condition, "MPO")
})

# ── facs_calc_diff ───────────────────────────────────────────

test_that("facs_calc_diff() computes proportion_stim - proportion_unstim with pseudocount = 0 by default", {
  result <- facs_calc_diff(log2fc_data, ref_pop = "CD3+")

  expect_equal(nrow(result), 1L)
  expect_equal(result$restimulation, "MPO")
  expect_equal(result$population, "CD4+")

  unstim_prop <- 100 / 1000
  mpo_prop    <- 400 / 1000
  expect_equal(result$diff, mpo_prop - unstim_prop)
})

test_that("facs_calc_diff() accepts a nonzero pseudocount", {
  result <- facs_calc_diff(log2fc_data, ref_pop = "CD3+", pseudocount = 0.5)

  unstim_prop <- (100 + 0.5) / 1000
  mpo_prop    <- (400 + 0.5) / 1000
  expect_equal(result$diff, mpo_prop - unstim_prop)
})

test_that("facs_calc_diff() excludes ref_pop and ref_level from the output", {
  result <- facs_calc_diff(log2fc_data, ref_pop = "CD3+")

  expect_false("CD3+" %in% result$population)
  expect_false("unstim" %in% result$restimulation)
})

test_that("facs_calc_diff() errors when ref_pop matches more than one row per mouse_ID/tissue/restim-level combo", {
  data_dup <- dplyr::bind_rows(
    log2fc_data,
    tibble::tibble(
      file_name = "f1", population_full_path = "OtherPath/CD3+", population = "CD3+",
      metric = "count", value = 500, mouse_ID = "m1", tissue = "spleen", restimulation = "unstim"
    )
  )

  expect_error(facs_calc_diff(data_dup, ref_pop = "CD3+"), regexp = "m1", fixed = TRUE)
})

test_that("facs_calc_diff() warns and fills NA when ref_level is missing for a mouse_ID/tissue group", {
  data_no_unstim <- dplyr::filter(log2fc_data, restimulation != "unstim")

  expect_warning(
    result <- facs_calc_diff(data_no_unstim, ref_pop = "CD3+"),
    regexp = "m1",
    fixed = TRUE
  )
  expect_true(is.na(result$diff))
})
