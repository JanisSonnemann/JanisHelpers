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

test_that("db_write_meta also takes the meta sheet itself and skips rows without mouse_ID", {
  con <- local_test_db()
  sheet <- meta_read(testthat::test_path("../fixtures/meta_minimal.xlsx"))$meta
  sheet$mouse_ID[2] <- NA

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
