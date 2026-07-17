test_that("db_write_facs creates a stain and inserts measurements (volumetric)", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  facs_data <- tibble::tibble(
    PopulationFullPath = c("root/Lymphocytes/CD4+", "root/Lymphocytes/CD8+"),
    Population = c("CD4+", "CD8+"),
    metric = c("Count", "Count"),
    value = c(1000, 800)
  )

  n <- db_write_facs(
    con, facs_data,
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    count_method = "volumetric", vol_stained = 0.1,
    vol_resuspended = 0.5, vol_measured = 0.05
  )
  expect_equal(n, 2)

  stains <- DBI::dbGetQuery(con, "SELECT * FROM facs_stains")
  expect_equal(nrow(stains), 1)
  expect_equal(stains$count_method, "volumetric")

  rows <- DBI::dbGetQuery(con, "SELECT population, value FROM facs_measurements ORDER BY population")
  expect_equal(rows$population, c("CD4+", "CD8+"))
  expect_equal(rows$value, c(1000, 800))
})

test_that("db_write_facs supports bead counting and rejects mixed fields", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  db_write_assay(con, assay_name = "treg", domain = "facs")

  facs_data <- tibble::tibble(
    PopulationFullPath = "root/Lymphocytes/CD4+",
    Population = "CD4+",
    metric = "Count",
    value = 1000
  )

  n <- db_write_facs(
    con, facs_data,
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    count_method = "bead", vol_stained = 0.1,
    bead_volume_added = 0.02, bead_concentration = 5000
  )
  expect_equal(n, 1)

  # Uses a different assay ("treg") than the call above, so this is a genuinely
  # new (sample_id, assay_id) pair -- the facs_stains INSERT is actually
  # attempted (not skipped by ON CONFLICT DO NOTHING), and the CHECK
  # constraint rejecting a mixed volumetric/bead field set fires as intended.
  expect_error(
    db_write_facs(
      con, facs_data,
      mouse_id = "25-7-1", tissue = "spleen", assay_name = "treg",
      count_method = "bead", vol_stained = 0.1,
      bead_volume_added = 0.02, bead_concentration = 5000,
      vol_measured = 0.05
    )
  )
})

test_that("db_write_facs errors clearly when the sample doesn't exist", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  facs_data <- tibble::tibble(
    PopulationFullPath = "root/Lymphocytes/CD4+", Population = "CD4+",
    metric = "Count", value = 1000
  )

  expect_error(
    db_write_facs(
      con, facs_data,
      mouse_id = "25-7-1", tissue = "kidney", assay_name = "overview",
      count_method = "volumetric", vol_stained = 0.1,
      vol_resuspended = 0.5, vol_measured = 0.05
    ),
    "No sample found"
  )
})

test_that("db_write_facs is idempotent", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  facs_data <- tibble::tibble(
    PopulationFullPath = "root/Lymphocytes/CD4+", Population = "CD4+",
    metric = "Count", value = 1000
  )
  args <- list(
    con = con, data = facs_data,
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    count_method = "volumetric", vol_stained = 0.1,
    vol_resuspended = 0.5, vol_measured = 0.05
  )
  do.call(db_write_facs, args)
  n <- do.call(db_write_facs, args)
  expect_equal(n, 0)
})

test_that("db_write_elisa inserts measurements", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  elisa_data <- tibble::tibble(
    cytokine = "MPO",
    sample_id = c("25-7-1", "25-7-1"),
    replicate = c(1L, 2L),
    value = c(12.3, 12.7),
    unit = "pg/ml",
    result_status = "OK"
  )

  n <- db_write_elisa(con, elisa_data, mouse_id = "25-7-1", tissue = "spleen", assay_name = "anti-MPO")
  expect_equal(n, 2)

  rows <- DBI::dbGetQuery(con, "SELECT replicate, value FROM elisa_measurements ORDER BY replicate")
  expect_equal(rows$replicate, c(1L, 2L))
  expect_equal(rows$value, c(12.3, 12.7))
})

test_that("db_write_elisa is idempotent", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  elisa_data <- tibble::tibble(
    cytokine = "MPO", sample_id = "25-7-1", replicate = 1L,
    value = 12.3, unit = "pg/ml", result_status = "OK"
  )
  args <- list(con = con, data = elisa_data, mouse_id = "25-7-1", tissue = "spleen", assay_name = "anti-MPO")
  do.call(db_write_elisa, args)
  n <- do.call(db_write_elisa, args)
  expect_equal(n, 0)
})
