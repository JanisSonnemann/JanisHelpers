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
    count_method = "volumetric", vol_stained = 0.1,
    vol_resuspended = 0.5, vol_measured = 0.05
  )

  result <- db_query_facs(con) |> dplyr::collect()

  expect_equal(result$mouse_id, "25-7-1")
  expect_equal(result$experiment_code, "25-7")
  expect_equal(result$tissue, "spleen")
  expect_equal(result$assay_name, "overview")
  expect_equal(result$population, "CD4+")
  expect_equal(result$value, 1000)
})

test_that("db_query_elisa returns a flat tibble keyed by natural keys", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  db_write_elisa(
    con,
    tibble::tibble(
      cytokine = "MPO", sample_id = "25-7-1", replicate = 1L,
      value = 12.3, unit = "pg/ml", result_status = "OK"
    ),
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "anti-MPO"
  )

  result <- db_query_elisa(con) |> dplyr::collect()

  expect_equal(result$mouse_id, "25-7-1")
  expect_equal(result$assay_name, "anti-MPO")
  expect_equal(result$cytokine, "MPO")
  expect_equal(result$value, 12.3)
})

test_that("db_query_histo returns a flat tibble and tolerates a NULL assay_id", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  db_write_histo(
    con, tibble::tibble(metric = "damage_score", value = 2),
    mouse_id = "25-7-1", tissue = "spleen"
  )

  result <- db_query_histo(con) |> dplyr::collect()

  expect_equal(result$mouse_id, "25-7-1")
  expect_equal(result$metric, "damage_score")
  expect_true(is.na(result$assay_name))
})
