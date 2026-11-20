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
    exclusion_set_id = 1L, mouse_id = "x-4", tissue = tissue, panel = NA_character_,
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
            "dob", "bmt_date", "end_date", "early", "reason", "excluded_mouse_wide",
            "exclusion_reason")
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

tb_overview <- tidyr::expand_grid(
  mouse_id = c("x-1", "x-2"), tissue = c("kidney", "spleen", "lung")
) |>
  dplyr::bind_rows(tibble::tibble(mouse_id = "x-3", tissue = "kidney")) |>
  dplyr::left_join(dplyr::select(tb_meta, mouse_id, collected_at = end_date), by = "mouse_id") |>
  dplyr::mutate(experiment_code = "x", wsp_root = "ov_x", population = "CD4+",
                metric = "count", value = 1)
tb_dat3 <- c(tb_dat, list(overview = tb_overview))

test_that("mouse grid shows tissues, date counts and status", {
  tb <- meta_mouse_table(tb_dat3, tb_wb)
  expect_s3_class(tb, "gt_tbl")
  d <- tb[["_data"]]
  expect_equal(names(d), c("mouse_id", "ear_tag", "group", "sex", "cage", "age_at_bmt_wk",
                           "end_day", "titer", "overview", "status"))
  expect_equal(nrow(d), 7)
  expect_equal(d$overview[d$mouse_id %in% "x-1"], "3/3")
  expect_equal(d$overview[d$mouse_id %in% "x-3"], "1/3 (no spleen, lung)")
  expect_equal(d$titer[d$mouse_id %in% "x-1"], "4/4")
  expect_equal(d$titer[d$mouse_id %in% "x-3"], "4/4")
  expect_equal(d$end_day[d$mouse_id %in% "x-3"], 29)
  expect_equal(d$status[d$mouse_id %in% "x-3"], "early: pain (d 29)")
  expect_equal(d$status[d$mouse_id %in% "x-4"], "excluded: no response")
  expect_equal(sum(grepl("^not in DB", d$status)), 3)
  expect_equal(d$status[d$ear_tag %in% "T7"], "not in DB: found dead")
  expect_equal(d$status[d$mouse_id %in% "x-1"], "ok")
})

test_that("blank cell for a unit without data for the mouse", {
  d <- meta_mouse_table(tb_dat3, tb_wb)[["_data"]]
  expect_equal(d$overview[d$mouse_id %in% "x-4"], "")
  expect_equal(d$titer[d$ear_tag %in% "T5"], "")
})

test_that("partial exclusion shows as a note", {
  dat <- tb_dat3
  dat$exclusions <- dplyr::mutate(tb_excl(tissue = "kidney"), mouse_id = "x-1")
  d <- meta_mouse_table(dat, tb_wb)[["_data"]]
  expect_equal(d$status[d$mouse_id %in% "x-1"], "ok; excl. kidney (no response)")
})

test_that("mice sort by group, then naturally by mouse id", {
  wb <- dplyr::mutate(tb_wb[1:2, ], mouse_ID = c("26-1-10", "26-1-2"), group = factor("irrel"))
  dat <- list(meta = dplyr::mutate(tb_meta[1:2, ], mouse_id = c("26-1-10", "26-1-2")),
              exclusions = tb_excl()[0, ])
  d <- meta_mouse_table(dat, wb)[["_data"]]
  expect_equal(d$mouse_id, c("26-1-2", "26-1-10"))
})

test_that("group order follows dat$meta levels", {
  d <- meta_mouse_table(tb_dat3, tb_wb)[["_data"]]
  expect_equal(unique(d$group), c("MPO409", "irrel", NA))
})

test_that("sex counts ignore case and report other codings", {
  wb <- dplyr::mutate(tb_wb, sex = c("F", "m", "f", "F", "w", "M", "?"))
  tb <- meta_summary_table(wb, tb_dat)
  expect_equal(tb_row(tb, "Sex")[["total"]], "4 / 2 / 1 other")
  expect_equal(tb_row(meta_summary_table(tb_wb, tb_dat), "Sex")[["total"]], "5 / 2")
})

test_that("early mice without death_date warn and never print NA", {
  wb <- dplyr::mutate(tb_wb, death_date = dplyr::if_else(id == "T7", as.Date(NA), death_date))
  expect_warning(tb <- meta_summary_table(wb, tb_dat), "T7")
  d <- tb[["_data"]]
  expect_false(any(grepl("NA", unlist(d[setdiff(names(d), "characteristic")]))))
})

test_that("a mouse-wide rule excludes the mouse whatever units exist", {
  d <- meta_mouse_table(c(tb_dat, list(overview = tb_overview)), tb_wb)[["_data"]]
  expect_equal(d$status[d$mouse_id %in% "x-4"], "excluded: no response")
})

test_that("grid end dates come from the DB when the workbook lacks them", {
  wb <- dplyr::mutate(tb_wb, death_date = dplyr::if_else(id == "T1", as.Date(NA), death_date))
  d <- meta_mouse_table(tb_dat3, wb)[["_data"]]
  expect_equal(d$titer[!is.na(d$mouse_id) & d$mouse_id %in% c("x-1", "x-2", "x-3")],
               rep("4/4", 3))
})

test_that("derived tibbles without import unit are not grid columns", {
  dat <- c(tb_dat3, list(titer_fc = dplyr::select(tb_titer, mouse_id, collected_at, value)))
  d <- meta_mouse_table(dat, tb_wb)[["_data"]]
  expect_false("titer_fc" %in% names(d))
})

# DB-sourced cohort: x-1..x-3 ok/early, x-4 excluded mouse-wide, D-1 a facility-ID
# dropout before BMT.
tb_cohort <- dplyr::bind_rows(
  dplyr::mutate(tb_meta, early = c(FALSE, FALSE, TRUE), early_reason = c(NA, NA, "pain")),
  tibble::tibble(experiment_code = "x", mouse_id = c("x-4", "D-1"),
                 group = factor(c("MPO409", "irrel"), levels = c("MPO409", "irrel")),
                 sex = "f", bmt_date = as.Date("2026-05-12"),
                 end_date = as.Date(c("2026-06-16", "2026-05-05")),
                 early = c(FALSE, TRUE), early_reason = c(NA, "died at bleed"))
) |>
  dplyr::mutate(id_source = dplyr::if_else(mouse_id == "D-1", "facility", "mouse_ID"),
                excluded_mouse_wide = mouse_id == "x-4",
                exclusion_reason = dplyr::if_else(mouse_id == "x-4", "no response", NA_character_))
tb_dat_db <- c(list(meta = tb_meta, cohort = tb_cohort), tb_dat3[c("exclusions", "titer")])

test_that("meta_summary_table uses dat$cohort when workbook is NULL", {
  expect_no_message(tb <- meta_summary_table(NULL, tb_dat_db))
  expect_equal(tb_total(tb, "Enrolled"), 5)
  expect_equal(tb_total(tb, "Dropped before"), 1)
  expect_equal(tb_total(tb, "Excluded"), 1)
})

test_that("facility-ID subjects get no mouse_id and the facility ID as ear tag", {
  co <- cohort_normalise_(tb_cohort)
  expect_true(is.na(co$mouse_id[5]))
  expect_equal(co$ear_tag[5], "D-1")
  expect_equal(co$mouse_id[1:4], tb_cohort$mouse_id[1:4])
  # a stored mdc_id stays the ear tag
  co2 <- cohort_normalise_(dplyr::mutate(tb_cohort, mdc_id = paste0("F", 1:5)))
  expect_equal(co2$ear_tag, paste0("F", 1:5))
  expect_true(is.na(co2$mouse_id[5]))
})

test_that("meta_mouse_table shows DB dropouts", {
  expect_no_message(d <- meta_mouse_table(tb_dat_db)[["_data"]])
  expect_match(d$status[d$ear_tag %in% "D-1"], "^early: died at bleed")
  expect_equal(d$status[d$mouse_id %in% "x-4"], "excluded: no response")
  expect_false(any(grepl("not in DB", d$status)))
})

test_that("cohort is not a grid column", {
  expect_false("cohort" %in% names(data_units_(tb_dat_db)))
})

test_that("cohort tables work on a real db_load_experiment() list", {
  con <- local_test_db()
  db_write_experiment(con, experiment_code = "26-1")
  db_write_subjects(con, tibble::tibble(
    mouse_id = c("26-1-1", "26-1-2", "D-1"), experiment_code = "26-1",
    treatment_group = c("a", "a", "b"), bmt_date = as.Date("2026-05-12"),
    end_date = as.Date(c("2026-06-15", "2026-06-15", "2026-05-05"))
  ))
  DBI::dbExecute(con, "UPDATE subjects SET id_source = 'facility', early = TRUE,
                       early_reason = 'died at bleed' WHERE mouse_id = 'D-1'")
  db_write_samples(con, tibble::tibble(mouse_id = c("26-1-1", "26-1-2"), tissue = "serum"))
  db_write_assay(con, assay_name = "anti-MPO", domain = "elisa")
  DBI::dbExecute(con, "INSERT INTO elisa_measurements
    (sample_id, assay_id, import_unit, analyte, replicate, value, unit, result_status)
    SELECT sm.sample_id, a.assay_id, 'titer_26-1', 'anti-MPO', 1, 5, 'U', 'ok'
    FROM samples sm, assays a WHERE a.assay_name = 'anti-MPO'")
  suppressMessages(db_write_exclusions(con, tibble::tibble(
    mouse_ID = "26-1-2", reason = "bad", decided_on = as.Date("2026-11-17")), "26-1"))
  dat <- suppressMessages(db_load_experiment(con, "26-1", group_levels = c("a", "b")))
  expect_no_message(tb <- meta_summary_table(NULL, dat))
  expect_equal(tb_total(tb, "Enrolled"), 3)
  expect_equal(tb_total(tb, "Dropped before"), 1)
  expect_equal(tb_total(tb, "Excluded"), 1)
  expect_equal(tb_total(tb, "Analysed"), 1)
  d <- meta_mouse_table(dat)[["_data"]]
  expect_equal(d$status[d$mouse_id %in% "26-1-2"], "excluded: bad")
  expect_match(d$status[d$ear_tag %in% "D-1"], "^early: died at bleed")
  expect_true(is.na(d$mouse_id[d$ear_tag %in% "D-1"]))
  expect_false(any(grepl("not in DB", d$status)))
})
