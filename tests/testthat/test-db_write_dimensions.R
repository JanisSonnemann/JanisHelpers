test_that("db_write_experiment inserts a new experiment", {
  con <- local_test_db()

  n <- db_write_experiment(con, experiment_code = "25-7", experiment_name = "Test experiment")
  expect_equal(n, 1)

  row <- DBI::dbGetQuery(con, "SELECT * FROM experiments WHERE experiment_code = '25-7'")
  expect_equal(nrow(row), 1)
  expect_equal(row$experiment_name, "Test experiment")
})

test_that("db_write_experiment is idempotent", {
  con <- local_test_db()

  db_write_experiment(con, experiment_code = "25-7")
  n <- db_write_experiment(con, experiment_code = "25-7", experiment_name = "Different name")
  expect_equal(n, 0)

  row <- DBI::dbGetQuery(con, "SELECT * FROM experiments WHERE experiment_code = '25-7'")
  expect_equal(nrow(row), 1)
  expect_true(is.na(row$experiment_name))
})

test_that("db_write_subjects inserts subjects linked to an experiment", {
  con <- local_test_db()
  db_write_experiment(con, experiment_code = "25-7")

  n <- db_write_subjects(con, tibble::tibble(
    mouse_id = c("25-7-1", "25-7-2"),
    experiment_code = "25-7",
    sex = c("f", "m")
  ))
  expect_equal(n, 2)

  rows <- DBI::dbGetQuery(con, "
    SELECT s.mouse_id, s.sex, e.experiment_code
    FROM subjects s JOIN experiments e ON e.experiment_id = s.experiment_id
    ORDER BY s.mouse_id
  ")
  expect_equal(rows$mouse_id, c("25-7-1", "25-7-2"))
  expect_equal(rows$sex, c("f", "m"))
  expect_equal(rows$experiment_code, c("25-7", "25-7"))
})

test_that("db_write_subjects is idempotent and fills missing optional columns", {
  con <- local_test_db()
  db_write_experiment(con, experiment_code = "25-7")

  db_write_subjects(con, tibble::tibble(mouse_id = "25-7-1", experiment_code = "25-7", end_date = as.Date("2026-06-10")))
  n <- db_write_subjects(con, tibble::tibble(mouse_id = "25-7-1", experiment_code = "25-7", end_date = as.Date("2026-06-10")))
  expect_equal(n, 0)

  row <- DBI::dbGetQuery(con, "SELECT * FROM subjects WHERE mouse_id = '25-7-1'")
  expect_equal(nrow(row), 1)
  expect_true(is.na(row$sex))
})

test_that("db_write_subjects errors clearly on an unknown experiment_code", {
  con <- local_test_db()

  expect_error(
    db_write_subjects(con, tibble::tibble(mouse_id = "25-7-1", experiment_code = "not-registered", end_date = as.Date("2026-06-10"))),
    "Unknown experiment_code"
  )

  row <- DBI::dbGetQuery(con, "SELECT * FROM subjects WHERE mouse_id = '25-7-1'")
  expect_equal(nrow(row), 0)
})

test_that("db_write_samples inserts samples and ignores weight columns", {
  con <- local_test_db()
  db_write_experiment(con, experiment_code = "25-7")
  db_write_subjects(con, tibble::tibble(mouse_id = "25-7-1", experiment_code = "25-7", end_date = as.Date("2026-06-10")))

  n <- db_write_samples(con, tibble::tibble(
    mouse_id = c("25-7-1", "25-7-1"),
    tissue = c("spleen", "serum"),
    total_weight = c(0.12, NA),
    facs_weight = c(0.05, NA),
    vol_total = c(2.5, NA)
  ))
  expect_equal(n, 2)

  samples <- DBI::dbGetQuery(con, "SELECT tissue FROM samples ORDER BY tissue")
  expect_equal(samples$tissue, c("serum", "spleen"))
  expect_false(any(c("total_weight", "facs_weight", "vol_total") %in%
    DBI::dbListFields(con, "samples")))
})

test_that("db_write_samples is idempotent", {
  con <- local_test_db()
  db_write_experiment(con, experiment_code = "25-7")
  db_write_subjects(con, tibble::tibble(mouse_id = "25-7-1", experiment_code = "25-7", end_date = as.Date("2026-06-10")))

  db_write_samples(con, tibble::tibble(mouse_id = "25-7-1", tissue = "spleen"))
  n <- db_write_samples(con, tibble::tibble(mouse_id = "25-7-1", tissue = "spleen"))
  expect_equal(n, 0)
})

test_that("db_write_samples errors clearly on an unknown mouse_id", {
  con <- local_test_db()

  expect_error(
    db_write_samples(con, tibble::tibble(mouse_id = "not-registered", tissue = "spleen")),
    "Unknown mouse_id"
  )

  row <- DBI::dbGetQuery(con, "SELECT * FROM samples")
  expect_equal(nrow(row), 0)
})

test_that("db_write_assay inserts an assay scoped to a domain", {
  con <- local_test_db()

  n <- db_write_assay(con, assay_name = "overview", domain = "facs")
  expect_equal(n, 1)

  row <- DBI::dbGetQuery(con, "
    SELECT a.assay_name, d.domain_name
    FROM assays a JOIN domains d ON d.domain_id = a.domain_id
    WHERE a.assay_name = 'overview'
  ")
  expect_equal(row$domain_name, "facs")
})

test_that("db_write_assay is idempotent", {
  con <- local_test_db()
  db_write_assay(con, assay_name = "overview", domain = "facs")
  n <- db_write_assay(con, assay_name = "overview", domain = "facs")
  expect_equal(n, 0)
})

test_that("db_write_assay errors clearly on an unknown domain", {
  con <- local_test_db()
  expect_error(
    db_write_assay(con, assay_name = "overview", domain = "not-a-domain"),
    "Unknown domain"
  )
})

# ---- db_write_meta ------------------------------------------------------------

test_that("db_write_meta registers experiment + subjects from the meta sheet", {
  con <- local_test_db()
  meta <- meta_read(testthat::test_path("../fixtures/meta_minimal.xlsx"))

  out <- suppressMessages(
    db_write_meta(con, meta, "26-1", "BMT-mRNA-V4", project = "MPO_mRNA")
  )

  expect_identical(out, con)
  exp <- DBI::dbGetQuery(con, "SELECT * FROM experiments")
  expect_equal(exp$experiment_code, "26-1")
  expect_equal(exp$experiment_name, "BMT-mRNA-V4")
  expect_equal(exp$project, "MPO_mRNA")
  s <- DBI::dbGetQuery(con, "SELECT * FROM subjects WHERE mouse_id = '26-1-1'")
  expect_equal(s$mdc_id, "LFS0-21559")
  expect_equal(s$cage, "ELFS0-00286")
  expect_equal(s$mouse_strain, "B6-Mpo")
  expect_equal(s$generation, "F3")
  expect_equal(s$sex, "m")
  expect_equal(s$treatment_group, "irrel")
  expect_equal(s$dob, as.Date("2025-11-26"))
  expect_equal(s$start_date, as.Date("2026-03-30"))
  expect_equal(s$bmt_date, as.Date("2026-05-01"))
  expect_equal(s$end_date, as.Date("2026-06-15"))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM subjects")$n, 4)
})

test_that("db_write_meta also takes the meta sheet itself and skips fully blank rows", {
  con <- local_test_db()
  sheet <- meta_read(testthat::test_path("../fixtures/meta_minimal.xlsx"))$meta
  sheet[2, ] <- NA  # a fully blank row

  suppressMessages(db_write_meta(con, sheet, "26-1"))

  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM subjects")$n, 3)
})

test_that("db_write_meta is idempotent", {
  con <- local_test_db()
  meta <- meta_read(testthat::test_path("../fixtures/meta_minimal.xlsx"))
  suppressMessages(db_write_meta(con, meta, "26-1"))
  suppressMessages(db_write_meta(con, meta, "26-1"))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM subjects")$n, 4)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM experiments")$n, 1)
})

test_that("db_write_meta lists every missing template column", {
  con <- local_test_db()
  sheet <- meta_read(testthat::test_path("../fixtures/meta_minimal.xlsx"))$meta
  sheet <- dplyr::select(sheet, !c(line, bmt_date))
  expect_error(db_write_meta(con, sheet, "26-1"), "line, bmt_date")
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM experiments")$n, 0)
})

meta_sheet_ <- function(...) {
  tibble::tibble(
    mouse_ID = c("m1", NA), id = c("F1", "F2"), cage = "c1", line = "L",
    generation = "F1", sex = "m", group = "g", dob = as.Date("2025-01-01"),
    start_date = as.Date("2025-02-01"), bmt_date = as.Date("2025-03-01"),
    death_date = as.Date(c("2025-06-01", NA)),
    early = c(NA, "yes"), reason = c(NA, "found dead")
  ) |> dplyr::mutate(...)
}

test_that("db_write_meta registers dropouts without mouse_ID by facility ID", {
  con <- local_test_db()
  suppressMessages(db_write_meta(con, meta_sheet_(), "26-1"))
  s <- DBI::dbGetQuery(con, "SELECT * FROM subjects ORDER BY subject_id")
  expect_equal(nrow(s), 2)
  expect_equal(s$mouse_id, c("m1", "F2"))
  expect_equal(s$id_source, c("mouse_ID", "facility"))
  expect_equal(s$early, c(FALSE, TRUE))
  expect_equal(s$early_reason, c(NA, "found dead"))
})

test_that("db_write_meta errors on a row with neither mouse_ID nor id", {
  con <- local_test_db()
  sheet <- meta_sheet_()
  sheet$id[2] <- NA
  expect_error(db_write_meta(con, sheet, "26-1"), "neither mouse_ID nor id")
  expect_equal(count_rows_(con, "subjects"), 0)
})

test_that("db_write_meta errors on unknown early values and accepts known ones", {
  con <- local_test_db()
  expect_error(db_write_meta(con, meta_sheet_(early = c(NA, "maybe")), "26-1"), "maybe")
  yes <- c("yes", "ja", "y", "x", "true", "TRUE", " Yes ")
  no <- c(NA, "", "no", "nein", "n", "false")
  expect_equal(parse_early_(c(yes, no)), rep(c(TRUE, FALSE), c(length(yes), length(no))))
  expect_equal(parse_early_(c(TRUE, FALSE, NA)), c(TRUE, FALSE, FALSE))
})

test_that("re-running db_write_meta updates early but refuses end_date changes", {
  con <- local_test_db()
  suppressMessages(db_write_meta(con, meta_sheet_(), "26-1"))
  suppressMessages(db_write_meta(
    con, meta_sheet_(early = c("yes", "yes"), reason = c("sick", "found dead")), "26-1"))
  s <- DBI::dbGetQuery(con, "SELECT * FROM subjects ORDER BY subject_id")
  expect_equal(s$early, c(TRUE, TRUE))
  expect_equal(s$early_reason, c("sick", "found dead"))

  err <- expect_error(
    db_write_meta(con, meta_sheet_(death_date = as.Date(c("2025-07-01", NA)),
                                   early = c("no", "no")), "26-1"),
    "end_date differs")
  expect_match(conditionMessage(err), "m1")
  expect_match(conditionMessage(err), "2025-06-01")
  expect_match(conditionMessage(err), "2025-07-01")
  s2 <- DBI::dbGetQuery(con, "SELECT * FROM subjects ORDER BY subject_id")
  expect_equal(s2$early, c(TRUE, TRUE))
  expect_equal(s2$end_date[1], as.Date("2025-06-01"))
})

test_that("db_write_meta refuses a dropout that later gets a mouse_ID", {
  con <- local_test_db()
  suppressMessages(db_write_meta(con, meta_sheet_(), "26-1"))
  err <- expect_error(
    db_write_meta(con, meta_sheet_(mouse_ID = c("m1", "m9")), "26-1"),
    "already registered")
  expect_match(conditionMessage(err), "'F2'")
  expect_match(conditionMessage(err), "m9")
  expect_equal(count_rows_(con, "subjects"), 2)
})

test_that("db_write_subjects still returns only newly inserted rows", {
  con <- local_test_db()
  db_write_experiment(con, "25-7")
  d <- tibble::tibble(mouse_id = "a", experiment_code = "25-7")
  expect_equal(db_write_subjects(con, d), 1)
  expect_equal(db_write_subjects(con, dplyr::mutate(d, early = TRUE)), 0)
})

test_that("db_write_samples defaults collected_at to the subject's end_date", {
  con <- local_test_db(); seed_minimal_fixture(con)
  got <- DBI::dbGetQuery(con, "SELECT collected_at FROM samples")$collected_at
  expect_equal(got, as.Date("2026-06-10"))
})

test_that("db_write_samples keeps one sample per mouse x tissue x date", {
  con <- local_test_db(); seed_minimal_fixture(con)
  n <- db_write_samples(con, tibble::tibble(
    mouse_id = "25-7-1", tissue = "serum",
    collected_at = as.Date(c("2026-05-06", "2026-05-20"))))
  expect_equal(n, 2)
  expect_equal(db_write_samples(con, tibble::tibble(mouse_id = "25-7-1", tissue = "spleen")), 0)
})

test_that("db_write_samples names every subject without end_date and writes nothing", {
  con <- local_test_db()
  db_write_experiment(con, "25-7")
  db_write_subjects(con, tibble::tibble(mouse_id = c("25-7-8", "25-7-9"), experiment_code = "25-7"))
  err <- expect_error(db_write_samples(con, tibble::tibble(mouse_id = c("25-7-8", "25-7-9"), tissue = "kidney")))
  expect_match(conditionMessage(err), "25-7-8.*25-7-9")
  expect_match(conditionMessage(err), "end_date")
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM samples")$n, 0)
})
