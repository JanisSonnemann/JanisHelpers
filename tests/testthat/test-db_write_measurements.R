overview_wsp_ <- "25-7/facs/overview_25-7_20261112.wsp"

cd4_data_ <- function() {
  tibble::tibble(
    PopulationFullPath = "root/Lymphocytes/CD4+", Population = "CD4+",
    metric = "Count", value = 1000
  )
}

test_that("db_write_facs creates a stain with wsp_root and inserts measurements", {
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
    wsp_root = "overview_25-7", source_file = overview_wsp_
  )
  expect_equal(n, 2)

  stains <- DBI::dbGetQuery(con, "SELECT * FROM facs_stains")
  expect_equal(nrow(stains), 1)
  expect_equal(stains$wsp_root, "overview_25-7")

  rows <- DBI::dbGetQuery(con, "SELECT population, value FROM facs_measurements ORDER BY population")
  expect_equal(rows$population, c("CD4+", "CD8+"))
  expect_equal(rows$value, c(1000, 800))
})

test_that("db_write_facs errors clearly when the sample doesn't exist", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  expect_error(
    db_write_facs(
      con, cd4_data_(),
      mouse_id = "25-7-1", tissue = "kidney", assay_name = "overview",
      wsp_root = "overview_25-7", source_file = overview_wsp_
    ),
    "No sample found"
  )
})

test_that("db_write_facs is idempotent", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  args <- list(
    con = con, data = cd4_data_(),
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    wsp_root = "overview_25-7", source_file = overview_wsp_
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
      wsp_root = "overview_25-7", source_file = overview_wsp_
    )
  )

  expect_equal(nrow(DBI::dbGetQuery(con, "SELECT * FROM facs_stains")), 0)
  expect_equal(nrow(DBI::dbGetQuery(con, "SELECT * FROM facs_measurements")), 0)
})

test_that("db_write_facs records the .wsp path as source_file, not the FCS filename", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  facs_data <- dplyr::mutate(cd4_data_(), FileName = "spleen_001.fcs")
  db_write_facs(
    con, facs_data,
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    wsp_root = "overview_25-7", source_file = overview_wsp_
  )

  rows <- DBI::dbGetQuery(con, "SELECT source_file FROM facs_measurements")
  expect_equal(rows$source_file, overview_wsp_)
})

test_that("db_write_facs requires a single non-empty source_file and writes nothing otherwise", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  args <- list(
    con = con, data = cd4_data_(),
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    wsp_root = "overview_25-7"
  )

  expect_error(do.call(db_write_facs, args), "source_file")
  expect_error(do.call(db_write_facs, c(args, source_file = NA_character_)), "source_file")
  expect_error(do.call(db_write_facs, c(args, source_file = "")), "source_file")
  expect_error(
    do.call(db_write_facs, c(args, list(source_file = c("a.wsp", "b.wsp")))),
    "source_file"
  )
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM facs_stains")$n, 0)
})

test_that("db_write_facs requires a single non-empty wsp_root", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  args <- list(
    con = con, data = cd4_data_(),
    mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    source_file = overview_wsp_
  )

  expect_error(do.call(db_write_facs, args), "wsp_root")
  expect_error(do.call(db_write_facs, c(args, wsp_root = NA_character_)), "wsp_root")
  expect_error(do.call(db_write_facs, c(args, wsp_root = "")), "wsp_root")
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM facs_stains")$n, 0)
})

seed_guard_fixture_ <- function(con) {
  seed_minimal_fixture(con)
  db_write_subjects(con, tibble::tibble(mouse_id = "25-7-2", experiment_code = "25-7"))
  db_write_samples(con, tibble::tibble(mouse_id = "25-7-2", tissue = "spleen"))
  db_write_experiment(con, experiment_code = "25-8")
  db_write_subjects(con, tibble::tibble(mouse_id = "25-8-1", experiment_code = "25-8"))
  db_write_samples(con, tibble::tibble(mouse_id = "25-8-1", tissue = "spleen"))
  db_write_assay(con, assay_name = "treg", domain = "facs")
  invisible(con)
}

write_cd4_ <- function(con, mouse_id, wsp_root, source_file, assay_name = "overview") {
  db_write_facs(
    con, cd4_data_(),
    mouse_id = mouse_id, tissue = "spleen", assay_name = assay_name,
    wsp_root = wsp_root, source_file = source_file
  )
}

facs_row_counts_ <- function(con) {
  c(
    stains = DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM facs_stains")$n,
    measurements = DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM facs_measurements")$n
  )
}

test_that("db_write_facs rejects a source_file that doesn't match wsp_root", {
  con <- local_test_db()
  seed_guard_fixture_(con)

  expect_error(
    write_cd4_(con, "25-7-1", "overview_25-7", "25-7/facs/treg_25-7_20261112.wsp"),
    "wsp_root"
  )
  expect_error(
    write_cd4_(con, "25-7-1", "overview_25-7", "25-7/facs/overview_25-7.wsp"),
    "wsp_root"
  )
  expect_equal(facs_row_counts_(con), c(stains = 0, measurements = 0))
})

test_that("db_write_facs matches the root literally", {
  con <- local_test_db()
  seed_guard_fixture_(con)
  db_write_assay(con, assay_name = "cd4+", domain = "facs")

  expect_error(
    write_cd4_(con, "25-7-1", "cd4+_25-7", "25-7/facs/cd44_25-7_20261112.wsp", "cd4+"),
    "wsp_root"
  )
  expect_equal(facs_row_counts_(con), c(stains = 0, measurements = 0))
  expect_no_error(
    write_cd4_(con, "25-7-1", "cd4+_25-7", "25-7/facs/cd4+_25-7_20261112.wsp", "cd4+")
  )
})

test_that("db_write_facs rejects a second root for the same experiment x panel", {
  con <- local_test_db()
  seed_guard_fixture_(con)
  write_cd4_(con, "25-7-1", "overview_25-7", overview_wsp_)
  before <- facs_row_counts_(con)

  expect_error(
    write_cd4_(con, "25-7-2", "overview2_25-7", "25-7/facs/overview2_25-7_20261112.wsp"),
    "overview_25-7"
  )
  expect_equal(facs_row_counts_(con), before)
})

test_that("db_write_facs rejects one root across two panels or experiments", {
  con <- local_test_db()
  seed_guard_fixture_(con)
  write_cd4_(con, "25-7-1", "overview_25-7", overview_wsp_)
  before <- facs_row_counts_(con)

  expect_error(
    write_cd4_(con, "25-7-2", "overview_25-7", overview_wsp_, assay_name = "treg"),
    "overview_25-7"
  )
  expect_error(
    write_cd4_(con, "25-8-1", "overview_25-7", overview_wsp_),
    "overview_25-7"
  )
  expect_equal(facs_row_counts_(con), before)
})

test_that("db_write_facs rejects a second workspace version of a root", {
  con <- local_test_db()
  seed_guard_fixture_(con)
  write_cd4_(con, "25-7-1", "overview_25-7", overview_wsp_)
  before <- facs_row_counts_(con)
  newer <- "25-7/facs/overview_25-7_20261120.wsp"

  expect_error(write_cd4_(con, "25-7-1", "overview_25-7", newer), "Stage 3|export")
  expect_error(write_cd4_(con, "25-7-2", "overview_25-7", newer), "Stage 3|export")
  expect_equal(facs_row_counts_(con), before)
})

test_that("db_write_facs accepts further samples of the same workspace", {
  con <- local_test_db()
  seed_guard_fixture_(con)

  write_cd4_(con, "25-7-1", "overview_25-7", overview_wsp_)
  expect_no_error(write_cd4_(con, "25-7-2", "overview_25-7", overview_wsp_))
  expect_equal(facs_row_counts_(con), c(stains = 2, measurements = 2))
})

test_that("every facs measurement has a matching stain", {
  con <- local_test_db()
  seed_guard_fixture_(con)
  write_cd4_(con, "25-7-1", "overview_25-7", overview_wsp_)
  write_cd4_(con, "25-7-2", "overview_25-7", overview_wsp_)
  write_cd4_(con, "25-7-1", "treg_25-7", "25-7/facs/treg_25-7_20261112.wsp", "treg")

  orphans <- DBI::dbGetQuery(con, "
    SELECT COUNT(*) AS n FROM facs_measurements m
    LEFT JOIN facs_stains s USING (facs_stain_id)
    WHERE s.facs_stain_id IS NULL
  ")$n
  expect_equal(orphans, 0)
})

seed_two_mouse_plate_ <- function(con) {
  seed_minimal_fixture(con)
  db_write_subjects(con, tibble::tibble(mouse_id = "25-7-2", experiment_code = "25-7"))
  db_write_samples(con, tibble::tibble(
    mouse_id = c("25-7-1", "25-7-2"), tissue = c("serum", "serum")
  ))
  invisible(con)
}

elisa_plate_ <- function(...) {
  plate <- tibble::tibble(
    cytokine = "MPO",
    sample_id = c("P1", "P1", "P2", "P3"),
    replicate = c(1L, 2L, 1L, 1L),
    value = c(12.3, 12.7, 30.1, 4.2),
    unit = "pg/ml",
    result_status = "OK",
    mouse_id = c("25-7-1", "25-7-1", "25-7-2", "25-7-1"),
    tissue = c("serum", "serum", "serum", "spleen")
  )
  overrides <- list(...)
  for (nm in names(overrides)) plate[[nm]] <- overrides[[nm]]
  plate
}

test_that("db_write_elisa writes a whole plate, resolving each row's sample", {
  con <- local_test_db()
  seed_two_mouse_plate_(con)

  n <- db_write_elisa(
    con, elisa_plate_(), assay_name = "anti-MPO",
    source_file = "25-7/elisa/Results_MPO.xlsx"
  )
  expect_equal(n, 4)

  rows <- DBI::dbGetQuery(con, "
    SELECT sub.mouse_id, sm.tissue, e.sample_id_raw, e.replicate, e.value, e.source_file
    FROM elisa_measurements e
    JOIN samples sm ON sm.sample_id = e.sample_id
    JOIN subjects sub ON sub.subject_id = sm.subject_id
    ORDER BY e.sample_id_raw, e.replicate
  ")
  expect_equal(rows$mouse_id, c("25-7-1", "25-7-1", "25-7-2", "25-7-1"))
  expect_equal(rows$tissue, c("serum", "serum", "serum", "spleen"))
  expect_equal(rows$sample_id_raw, c("P1", "P1", "P2", "P3"))
  expect_equal(rows$value, c(12.3, 12.7, 30.1, 4.2))
  expect_equal(unique(rows$source_file), "25-7/elisa/Results_MPO.xlsx")
})

test_that("db_write_elisa is idempotent", {
  con <- local_test_db()
  seed_two_mouse_plate_(con)

  db_write_elisa(con, elisa_plate_(), assay_name = "anti-MPO")
  n <- db_write_elisa(con, elisa_plate_(), assay_name = "anti-MPO")
  expect_equal(n, 0)
})

test_that("db_write_elisa lists every unresolved mouse_id x tissue and writes nothing", {
  con <- local_test_db()
  seed_two_mouse_plate_(con)

  plate <- elisa_plate_(
    mouse_id = c("25-7-1", "25-7-1", "25-7-9", NA),
    tissue = c("serum", "serum", "serum", "kidney")
  )
  err <- expect_error(db_write_elisa(con, plate, assay_name = "anti-MPO"))
  expect_match(conditionMessage(err), "25-7-9 x serum", fixed = TRUE)
  expect_match(conditionMessage(err), "NA x kidney", fixed = TRUE)
  expect_no_match(conditionMessage(err), "25-7-1 x serum", fixed = TRUE)

  expect_equal(nrow(DBI::dbGetQuery(con, "SELECT * FROM elisa_measurements")), 0)
})

test_that("db_write_elisa requires mouse_id and tissue columns, hinting at mouse_ID", {
  con <- local_test_db()
  seed_two_mouse_plate_(con)

  expect_error(
    db_write_elisa(con, dplyr::select(elisa_plate_(), !tissue), assay_name = "anti-MPO"),
    "tissue"
  )
  expect_error(
    db_write_elisa(con, dplyr::rename(elisa_plate_(), mouse_ID = mouse_id), assay_name = "anti-MPO"),
    "mouse_ID"
  )
})

test_that("db_write_elisa errors on NA sample_id or replicate and inserts nothing", {
  con <- local_test_db()
  seed_two_mouse_plate_(con)

  expect_error(
    db_write_elisa(con, elisa_plate_(sample_id = c("P1", NA, "P2", "P3")), assay_name = "anti-MPO"),
    "requires non-NA sample_id and replicate"
  )
  expect_error(
    db_write_elisa(con, elisa_plate_(replicate = c(1L, NA, 1L, 1L)), assay_name = "anti-MPO"),
    "requires non-NA sample_id and replicate"
  )

  expect_equal(nrow(DBI::dbGetQuery(con, "SELECT * FROM elisa_measurements")), 0)
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
