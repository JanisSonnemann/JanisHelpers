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

  db_write_subjects(con, tibble::tibble(mouse_id = "25-7-1", experiment_code = "25-7"))
  n <- db_write_subjects(con, tibble::tibble(mouse_id = "25-7-1", experiment_code = "25-7"))
  expect_equal(n, 0)

  row <- DBI::dbGetQuery(con, "SELECT * FROM subjects WHERE mouse_id = '25-7-1'")
  expect_equal(nrow(row), 1)
  expect_true(is.na(row$sex))
})
