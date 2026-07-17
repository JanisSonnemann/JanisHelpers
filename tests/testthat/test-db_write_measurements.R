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

test_that("db_write_facs rolls back the facs_stains row when the measurements insert fails", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  bad_facs_data <- tibble::tibble(
    PopulationFullPath = c("root/Lymphocytes/CD4+", NA_character_),
    Population = c("CD4+", "CD8+"),
    metric = c("Count", NA_character_),
    value = c(1000, 800)
  )

  expect_error(
    db_write_facs(
      con, bad_facs_data,
      mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
      count_method = "volumetric", vol_stained = 0.1,
      vol_resuspended = 0.5, vol_measured = 0.05
    )
  )

  stains <- DBI::dbGetQuery(
    con,
    "SELECT fs.* FROM facs_stains fs
     JOIN samples sm ON sm.sample_id = fs.sample_id
     JOIN subjects sub ON sub.subject_id = sm.subject_id
     JOIN assays a ON a.assay_id = fs.assay_id
     WHERE sub.mouse_id = '25-7-1' AND sm.tissue = 'spleen' AND a.assay_name = 'overview'"
  )
  expect_equal(nrow(stains), 0)

  measurements <- DBI::dbGetQuery(con, "SELECT * FROM facs_measurements")
  expect_equal(nrow(measurements), 0)
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

test_that("db_write_elisa errors on NA sample_id or replicate and inserts nothing", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  elisa_data_na_sample <- tibble::tibble(
    cytokine = "MPO", sample_id = NA_character_, replicate = 1L,
    value = 12.3, unit = "pg/ml", result_status = "OK"
  )
  expect_error(
    db_write_elisa(con, elisa_data_na_sample, mouse_id = "25-7-1", tissue = "spleen", assay_name = "anti-MPO"),
    "requires non-NA sample_id and replicate"
  )

  elisa_data_na_replicate <- tibble::tibble(
    cytokine = "MPO", sample_id = "25-7-1", replicate = NA_integer_,
    value = 12.3, unit = "pg/ml", result_status = "OK"
  )
  expect_error(
    db_write_elisa(con, elisa_data_na_replicate, mouse_id = "25-7-1", tissue = "spleen", assay_name = "anti-MPO"),
    "requires non-NA sample_id and replicate"
  )

  rows <- DBI::dbGetQuery(con, "SELECT * FROM elisa_measurements")
  expect_equal(nrow(rows), 0)
})

test_that("db_write_histo inserts measurements without an assay", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  histo_data <- tibble::tibble(metric = c("damage_score", "fibrosis_pct"), value = c(2, 15.5))

  n <- db_write_histo(con, histo_data, mouse_id = "25-7-1", tissue = "spleen")
  expect_equal(n, 2)

  rows <- DBI::dbGetQuery(con, "SELECT metric, value, assay_id FROM histo_measurements ORDER BY metric")
  expect_equal(rows$metric, c("damage_score", "fibrosis_pct"))
  expect_true(all(is.na(rows$assay_id)))
})

test_that("db_write_histo errors clearly on an unknown assay_name", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  histo_data <- tibble::tibble(metric = "damage_score", value = 2)
  expect_error(
    db_write_histo(con, histo_data, mouse_id = "25-7-1", tissue = "spleen", assay_name = "not-registered"),
    "No assay found"
  )
})
