local_test_db <- function(env = parent.frame()) {
  path <- tempfile(fileext = ".duckdb")
  con <- db_connect(path)
  withr::defer(DBI::dbDisconnect(con, shutdown = TRUE), envir = env)
  con
}
