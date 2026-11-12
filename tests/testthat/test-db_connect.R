test_that("db_connect creates all tables and seeds domains", {
  con <- local_test_db()

  expect_true(all(c(
    "domains", "experiments", "subjects", "samples", "sample_digestions",
    "assays", "facs_stains", "facs_measurements", "elisa_measurements",
    "histo_measurements"
  ) %in% DBI::dbListTables(con)))

  domain_names <- DBI::dbGetQuery(
    con, "SELECT domain_name FROM domains ORDER BY domain_name"
  )$domain_name
  expect_equal(domain_names, c("elisa", "facs", "histo"))
})

test_that("db_connect is idempotent", {
  path <- tempfile(fileext = ".duckdb")
  con1 <- db_connect(path)
  DBI::dbDisconnect(con1, shutdown = TRUE)

  con2 <- db_connect(path)
  withr::defer(DBI::dbDisconnect(con2, shutdown = TRUE))

  domain_count <- DBI::dbGetQuery(con2, "SELECT COUNT(*) AS n FROM domains")$n
  expect_equal(domain_count, 3)
})
