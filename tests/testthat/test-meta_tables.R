# Miniature 26-1 cohort. Workbook (meta_read() names): a, b ok; c ended early
# after BMT but has end-point data; d excluded mouse-wide; e died before BMT;
# f removed at BMT without group; g died after BMT without data.
tb_wb <- tibble::tibble(
  mouse_ID = c("x-1", "x-2", "x-3", "x-4", NA, NA, NA),
  id = paste0("T", 1:7),
  cage = c("c1", "c1", "c2", "c2", "c3", "c3", "c3"),
  line = "B6-Mpo", generation = "F3",
  sex = c("f", "m", "f", "f", "f", "m", "f"),
  group = factor(c("MPO409", "irrel", "irrel", "MPO409", "irrel", NA, "MPO409")),
  dob = as.Date("2026-01-01"),
  start_date = as.Date("2026-03-30"),
  bmt_date = as.Date(c(rep("2026-05-12", 4), NA, "2026-05-12", "2026-05-12")),
  death_date = as.Date(c("2026-06-15", "2026-06-15", "2026-06-10", "2026-06-16",
                         "2026-05-05", "2026-05-12", "2026-05-22")),
  early = c(NA, NA, "yes", NA, "yes", "yes", "yes"),
  reason = c(NA, NA, "pain", NA, "died at bleed", "bad titer", "found dead")
)

tb_meta <- tibble::tibble(
  experiment_code = "x", mouse_id = c("x-1", "x-2", "x-3"),
  group = factor(c("MPO409", "irrel", "irrel"), levels = c("MPO409", "irrel")),
  mouse_treatment = NA_character_, sex = c("f", "m", "f"), cage = c("c1", "c1", "c2"),
  mouse_strain = "B6-Mpo", generation = "F3", dob = as.Date("2026-01-01"),
  start_date = as.Date("2026-03-30"), bmt_date = as.Date("2026-05-12"),
  end_date = as.Date(c("2026-06-15", "2026-06-15", "2026-06-10"))
)

tb_excl <- function(tissue = NA_character_) {
  tibble::tibble(
    import_unit = "titer_x", mouse_id = "x-4", tissue = tissue, panel = NA_character_,
    condition = NA_character_, file_name = NA_character_, collected_at = as.Date(NA),
    reason = "no response", decided_on = as.Date("2026-11-17")
  )
}

tb_titer <- tidyr::expand_grid(
  mouse_id = c("x-1", "x-2", "x-3"),
  collected_at = as.Date(c("2026-05-05", "2026-05-19", "2026-06-02"))
) |>
  dplyr::bind_rows(dplyr::select(tb_meta, mouse_id, collected_at = end_date)) |>
  dplyr::mutate(experiment_code = "x", tissue = "serum", import_unit = "titer_x",
                analyte = "anti-MPO", value = 1)

tb_dat <- list(meta = tb_meta, exclusions = tb_excl(), titer = tb_titer)

tb_row <- function(tbl, label) {
  d <- tbl[["_data"]]
  unlist(d[startsWith(d$characteristic, label), setdiff(names(d), "characteristic")])
}
tb_total <- function(tbl, label) as.integer(sub("^([0-9]+).*", "\\1", tb_row(tbl, label)[["total"]]))

test_that("cohort normaliser accepts workbook and DB names", {
  cols <- c("mouse_id", "ear_tag", "group", "sex", "cage", "strain", "generation",
            "dob", "bmt_date", "end_date", "early", "reason")
  wb <- cohort_normalise_(tb_wb)
  expect_named(wb, cols, ignore.order = TRUE)
  expect_equal(wb$ear_tag, paste0("T", 1:7))
  expect_equal(wb$early, c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE))
  db <- cohort_normalise_(tb_meta)
  expect_named(db, cols, ignore.order = TRUE)
  expect_equal(db$mouse_id, tb_meta$mouse_id)
  expect_equal(db$strain, rep("B6-Mpo", 3))
  expect_false(any(db$early))
  expect_s3_class(cohort_normalise_(dplyr::mutate(tb_wb, death_date = as.POSIXct(death_date)))$end_date, "Date")
})

test_that("analysis status: excluded mouse-wide, early with end-point counts", {
  st <- analysis_status_(cohort_normalise_(tb_wb), tb_dat)
  i <- match(c("x-1", "x-2", "x-3", "x-4"), st$mouse_id)
  expect_equal(st$analysed[i], c(TRUE, TRUE, TRUE, FALSE))
  expect_equal(st$excluded[i], c(FALSE, FALSE, FALSE, TRUE))
  expect_equal(st$exclusion_reason[i[4]], "no response")
  expect_equal(st$in_db, c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_false(any(st$analysed[is.na(st$mouse_id)]))
})

test_that("partial exclusion does not exclude the mouse", {
  dat <- tb_dat
  dat$exclusions <- tb_excl(tissue = "kidney")
  st <- analysis_status_(cohort_normalise_(tb_wb), dat)
  expect_false(st$excluded[st$mouse_id %in% "x-4"])
})

test_that("summary table counts the whole cohort", {
  tb <- meta_summary_table(tb_wb, tb_dat)
  expect_s3_class(tb, "gt_tbl")
  expect_equal(names(tb[["_data"]]), c("characteristic", "MPO409", "irrel", "not assigned", "total"))
  expect_equal(unname(tb_row(tb, "Enrolled")), c("3", "3", "1", "7"))
  expect_equal(unname(tb_row(tb, "Sex")), c("3 / 0", "2 / 1", "0 / 1", "5 / 2"))
  expect_equal(tb_total(tb, "Analysed"), 3)
  expect_equal(tb_total(tb, "Dropped before / at BMT"), 2)
  expect_equal(tb_total(tb, "Ended early after BMT"), 2)
  expect_equal(tb_total(tb, "Excluded from analysis"), 1)
  expect_match(tb_row(tb, "Excluded from analysis")[["MPO409"]], "x-4")
  expect_equal(tb_row(tb, "Age at BMT")[["total"]], "18.7 (18.7-18.7)")
})

test_that("without dat the Analysed and Excluded rows are omitted", {
  tb <- meta_summary_table(tb_wb)
  labels <- tb[["_data"]]$characteristic
  expect_false(any(startsWith(labels, "Analysed")))
  expect_false(any(startsWith(labels, "Excluded")))
  expect_true(any(startsWith(labels, "Enrolled")))
})

test_that("workbook = NULL falls back to dat$meta with a message", {
  expect_message(tb <- meta_summary_table(NULL, tb_dat), "not shown")
  expect_equal(tb_total(tb, "Enrolled"), 3)
})

test_that("no 'not assigned' column when every mouse has a group", {
  tb <- meta_summary_table(tb_wb[!is.na(tb_wb$group), ], tb_dat)
  expect_false("not assigned" %in% names(tb[["_data"]]))
})
