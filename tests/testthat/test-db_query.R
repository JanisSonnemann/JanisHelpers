test_that("db_query_facs returns a flat tibble keyed by natural keys", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  db_write_facs(
    con,
    tibble::tibble(
      PopulationFullPath = "root/Lymphocytes/CD4+", Population = "CD4+",
      metric = "Count", value = 1000
    ),
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    wsp_root = "overview_25-7",
    source_file = "25-7/facs/overview_25-7_20261112.wsp"
  )

  result <- db_query_facs(con) |> dplyr::collect()

  expect_equal(result$wsp_root, "overview_25-7")
  expect_equal(result$condition, "none")
  expect_true(is.na(result$condition_keyword))
  expect_false("count_method" %in% names(result))
  expect_equal(result$mouse_id, "25-7-1")
  expect_equal(result$experiment_code, "25-7")
  expect_equal(result$tissue, "spleen")
  expect_equal(result$collected_at, as.Date("2026-06-10"))
  expect_equal(result$assay_name, "overview")
  expect_equal(result$population, "CD4+")
  expect_equal(result$value, 1000)
})

test_that("db_query_elisa and db_query_histo return flat tibbles keyed by natural keys", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  db_write_assay(con, "anti-MPO", "elisa")
  db_write_assay(con, "scores", "histo")
  DBI::dbExecute(con, "INSERT INTO elisa_measurements
    (sample_id, assay_id, import_unit, analyte, replicate, value, unit, result_status)
    SELECT sm.sample_id, a.assay_id, 'Results_MPO', 'MPO', 1, 12.3, 'pg/ml', 'ok'
    FROM samples sm, assays a WHERE a.assay_name = 'anti-MPO' LIMIT 1")
  DBI::dbExecute(con, "INSERT INTO histo_measurements
    (sample_id, assay_id, import_unit, metric, value)
    SELECT sm.sample_id, a.assay_id, 'scores', 'damage_score', 2
    FROM samples sm, assays a WHERE a.assay_name = 'scores' LIMIT 1")

  elisa <- dplyr::collect(db_query_elisa(con))
  expect_named(elisa, c("experiment_code", "mouse_id", "tissue", "collected_at",
    "assay_name", "import_unit", "analyte", "replicate", "value", "unit",
    "result_status", "lloq", "uloq", "source_file", "imported_at"))
  expect_equal(elisa$analyte, "MPO")
  expect_equal(elisa$assay_name, "anti-MPO")
  histo <- dplyr::collect(db_query_histo(con))
  expect_named(histo, c("experiment_code", "mouse_id", "tissue", "collected_at",
    "assay_name", "import_unit", "metric", "value", "source_file", "imported_at"))
  expect_equal(histo$metric, "damage_score")
})

test_that("db_query_facs returns population_order (NA for character labels)", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  db_write_facs(
    con,
    tibble::tibble(PopulationFullPath = "root/CD4+", Population = "CD4+",
                   metric = "count", value = 1),
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    wsp_root = "overview_25-7", source_file = "25-7/overview_25-7_20261112.wsp"
  )
  result <- dplyr::collect(db_query_facs(con))
  expect_true("population_order" %in% names(result))
  expect_true(is.na(result$population_order))
})

test_that("db_write_facs stores population_order when Population is a factor", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  db_write_facs(
    con,
    tibble::tibble(PopulationFullPath = "root/CD3+/CD4+",
                   Population = factor("CD4+", levels = c("CD3+", "CD4+")),
                   metric = "count", value = 1),
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    wsp_root = "overview_25-7", source_file = "25-7/overview_25-7_20261112.wsp"
  )
  result <- dplyr::collect(db_query_facs(con))
  expect_equal(result$population_order, 2L)
  expect_equal(result$population, "CD4+")
})

test_that("db_query_subjects returns one row per mouse with group", {
  con <- local_test_db()
  db_write_experiment(con, experiment_code = "25-7")
  db_write_subjects(con, tibble::tibble(
    mouse_id = c("25-7-1", "25-7-2"), experiment_code = "25-7",
    treatment_group = c("irrel", "MPO409"), end_date = as.Date("2026-06-10")
  ))

  lazy <- db_query_subjects(con)
  expect_s3_class(lazy, "tbl_lazy")
  irrel <- lazy |> dplyr::filter(group == "irrel") |> dplyr::collect()
  expect_equal(irrel$mouse_id, "25-7-1")

  all <- dplyr::collect(lazy)
  expect_equal(names(all), c("experiment_code", "mouse_id", "group",
    "mouse_treatment", "sex", "cage", "mouse_strain", "generation", "dob",
    "start_date", "bmt_date", "end_date"))
  expect_equal(nrow(all), 2)
  expect_false("treatment_group" %in% names(all))
})

test_that("db_query_facs falls back when population_order is missing (read-only)", {
  path <- tempfile(fileext = ".duckdb")
  con <- db_connect(path)
  seed_minimal_fixture(con)
  db_write_facs(
    con,
    tibble::tibble(PopulationFullPath = "root/CD4+", Population = "CD4+",
                   metric = "count", value = 1),
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    wsp_root = "overview_25-7", source_file = "25-7/overview_25-7_20261112.wsp"
  )
  DBI::dbExecute(con, "ALTER TABLE facs_measurements DROP COLUMN population_order")
  DBI::dbDisconnect(con, shutdown = TRUE)

  ro <- DBI::dbConnect(duckdb::duckdb(), path, read_only = TRUE)
  withr::defer(DBI::dbDisconnect(ro, shutdown = TRUE))
  res <- dplyr::collect(db_query_facs(ro))
  expect_true("population_order" %in% names(res))
  expect_true(all(is.na(res$population_order)))
  expect_type(res$population_order, "integer")
  expect_equal(res$population, "CD4+")
})
