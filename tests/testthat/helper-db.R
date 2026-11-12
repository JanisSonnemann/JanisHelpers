local_test_db <- function(env = parent.frame()) {
  path <- tempfile(fileext = ".duckdb")
  con <- db_connect(path)
  withr::defer(DBI::dbDisconnect(con, shutdown = TRUE), envir = env)
  con
}

seed_minimal_fixture <- function(con) {
  db_write_experiment(con, experiment_code = "25-7")
  db_write_subjects(con, tibble::tibble(mouse_id = "25-7-1", experiment_code = "25-7"))
  db_write_samples(con, tibble::tibble(mouse_id = "25-7-1", tissue = "spleen"))
  db_write_assay(con, assay_name = "overview", domain = "facs")
  db_write_assay(con, assay_name = "anti-MPO", domain = "elisa")
  invisible(con)
}
