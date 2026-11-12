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

test_that("db_write_facs rejects wrong-shaped or empty data and writes nothing", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  args <- list(
    con = con, mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
    wsp_root = "overview_25-7", source_file = overview_wsp_
  )

  snake <- tibble::tibble(
    population_full_path = "root/CD4+", population = "CD4+", metric = "Count", value = 1
  )
  expect_error(do.call(db_write_facs, c(args, list(data = snake))), "PopulationFullPath")
  expect_error(do.call(db_write_facs, c(args, list(data = cd4_data_()[0, ]))), "no rows")
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM facs_stains")$n, 0)
})

test_that("db_write_facs rejects data from more than one FCS file", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  two_files <- dplyr::bind_rows(
    dplyr::mutate(cd4_data_(), FileName = "spleen_001.fcs"),
    dplyr::mutate(cd4_data_(), FileName = "spleen_002.fcs", value = 2000)
  )
  expect_error(
    db_write_facs(
      con, two_files, mouse_id = "25-7-1", tissue = "spleen", assay_name = "overview",
      wsp_root = "overview_25-7", source_file = overview_wsp_
    ),
    "FileName"
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

restim_wsp_ <- "25-7/facs/restim_25-7_20261112.wsp"

write_restim_ <- function(con, condition, condition_keyword = "restimulation",
                          mouse_id = "25-7-1", wsp_root = "restim_25-7",
                          source_file = restim_wsp_) {
  db_write_facs(
    con, cd4_data_(),
    mouse_id = mouse_id, tissue = "spleen", assay_name = "restim",
    wsp_root = wsp_root, source_file = source_file,
    condition = condition, condition_keyword = condition_keyword
  )
}

test_that("db_write_facs defaults to condition 'none' with no keyword", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  write_cd4_(con, "25-7-1", "overview_25-7", overview_wsp_)

  stain <- DBI::dbGetQuery(con, "SELECT condition, condition_keyword FROM facs_stains")
  expect_equal(stain$condition, "none")
  expect_true(is.na(stain$condition_keyword))
})

test_that("db_write_facs stores one stain per condition of the same sample x panel", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  db_write_assay(con, assay_name = "restim", domain = "facs")

  for (cond in c("unstim", "MPO409", "PMA-Iono")) {
    expect_equal(write_restim_(con, cond), 1)
  }

  stains <- DBI::dbGetQuery(con, "
    SELECT condition, condition_keyword, wsp_root FROM facs_stains ORDER BY condition
  ")
  expect_equal(stains$condition, c("MPO409", "PMA-Iono", "unstim"))
  expect_equal(unique(stains$condition_keyword), "restimulation")
  expect_equal(unique(stains$wsp_root), "restim_25-7")
  expect_equal(facs_row_counts_(con), c(stains = 3, measurements = 3))
  expect_equal(write_restim_(con, "MPO409"), 0)
})

test_that("db_write_facs requires condition and condition_keyword together", {
  con <- local_test_db()
  seed_minimal_fixture(con)
  db_write_assay(con, assay_name = "restim", domain = "facs")

  expect_error(write_restim_(con, "MPO409", condition_keyword = NA_character_), "condition_keyword")
  expect_error(write_restim_(con, "none", condition_keyword = "restimulation"), "condition")
  expect_error(write_restim_(con, c("a", "b")), "condition")
  expect_error(write_restim_(con, NA_character_), "condition")
  expect_equal(facs_row_counts_(con), c(stains = 0, measurements = 0))
})

test_that("db_write_facs rejects a root that mixes condition keywords", {
  con <- local_test_db()
  seed_guard_fixture_(con)
  db_write_assay(con, assay_name = "restim", domain = "facs")
  write_restim_(con, "MPO409")
  before <- facs_row_counts_(con)

  expect_error(write_restim_(con, "24h", condition_keyword = "timepoint"), "keyword")
  expect_error(
    write_restim_(con, "none", condition_keyword = NA_character_, mouse_id = "25-7-2"),
    "keyword"
  )
  expect_equal(facs_row_counts_(con), before)
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

elisa_src_ <- "25-7/elisa/Results_MPO.xlsx"

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

  db_write_elisa(con, elisa_plate_(), assay_name = "anti-MPO", source_file = elisa_src_)
  n <- db_write_elisa(con, elisa_plate_(), assay_name = "anti-MPO", source_file = elisa_src_)
  expect_equal(n, 0)
})

test_that("db_write_elisa lists every unresolved mouse_id x tissue and writes nothing", {
  con <- local_test_db()
  seed_two_mouse_plate_(con)

  plate <- elisa_plate_(
    mouse_id = c("25-7-1", "25-7-1", "25-7-9", NA),
    tissue = c("serum", "serum", "serum", "kidney")
  )
  err <- expect_error(db_write_elisa(con, plate, assay_name = "anti-MPO", source_file = elisa_src_))
  expect_match(conditionMessage(err), "25-7-9 x serum", fixed = TRUE)
  expect_match(conditionMessage(err), "NA x kidney", fixed = TRUE)
  expect_no_match(conditionMessage(err), "25-7-1 x serum", fixed = TRUE)

  expect_equal(nrow(DBI::dbGetQuery(con, "SELECT * FROM elisa_measurements")), 0)
})

test_that("db_write_elisa requires mouse_id and tissue columns, hinting at mouse_ID", {
  con <- local_test_db()
  seed_two_mouse_plate_(con)

  expect_error(
    db_write_elisa(con, dplyr::select(elisa_plate_(), !tissue), assay_name = "anti-MPO", source_file = elisa_src_),
    "tissue"
  )
  expect_error(
    db_write_elisa(con, dplyr::rename(elisa_plate_(), mouse_ID = mouse_id), assay_name = "anti-MPO", source_file = elisa_src_),
    "mouse_ID"
  )
})

test_that("db_write_elisa errors on NA sample_id or replicate and inserts nothing", {
  con <- local_test_db()
  seed_two_mouse_plate_(con)

  expect_error(
    db_write_elisa(con, elisa_plate_(sample_id = c("P1", NA, "P2", "P3")), assay_name = "anti-MPO", source_file = elisa_src_),
    "requires non-NA sample_id and replicate"
  )
  expect_error(
    db_write_elisa(con, elisa_plate_(replicate = c(1L, NA, 1L, 1L)), assay_name = "anti-MPO", source_file = elisa_src_),
    "requires non-NA sample_id and replicate"
  )

  expect_equal(nrow(DBI::dbGetQuery(con, "SELECT * FROM elisa_measurements")), 0)
})

read_log_ <- function(con) {
  DBI::dbGetQuery(con, "
    SELECT l.*, d.domain_name FROM import_log l
    JOIN domains d ON d.domain_id = l.domain_id
    ORDER BY l.import_id
  ")
}

test_that("db_write_elisa writes one import_log row per call", {
  con <- local_test_db()
  seed_two_mouse_plate_(con)

  db_write_elisa(con, elisa_plate_(), assay_name = "anti-MPO", source_file = elisa_src_)
  log <- read_log_(con)
  expect_equal(nrow(log), 1)
  expect_equal(log$domain_name, "elisa")
  expect_equal(log$experiment_code, "25-7")
  expect_equal(log$import_unit, "Results_MPO.xlsx")
  expect_equal(log$source_file, elisa_src_)
  expect_equal(log$rows_written, 4L)
  expect_true(is.na(log$qc_passed))
  expect_true(is.na(log$qc_report))
  expect_true(is.na(log$source_version))
  expect_true(is.na(log$supersedes_import_id))

  db_write_elisa(con, elisa_plate_(), assay_name = "anti-MPO", source_file = elisa_src_)
  log <- read_log_(con)
  expect_equal(nrow(log), 2)
  expect_equal(log$rows_written[[2]], 0L)
})

test_that("db_write_elisa logs every experiment on a mixed plate", {
  con <- local_test_db()
  seed_two_mouse_plate_(con)
  db_write_experiment(con, experiment_code = "25-8")
  db_write_subjects(con, tibble::tibble(mouse_id = "25-8-1", experiment_code = "25-8"))
  db_write_samples(con, tibble::tibble(mouse_id = "25-8-1", tissue = "serum"))

  plate <- elisa_plate_(
    mouse_id = c("25-8-1", "25-8-1", "25-7-2", "25-7-1"),
    tissue = c("serum", "serum", "serum", "serum")
  )
  db_write_elisa(con, plate, assay_name = "anti-MPO", source_file = elisa_src_)
  expect_equal(read_log_(con)$experiment_code, "25-7, 25-8")
})

test_that("db_write_elisa writes no log row when the data insert fails", {
  con <- local_test_db()
  seed_two_mouse_plate_(con)

  plate <- elisa_plate_(mouse_id = c("25-7-1", "25-7-1", "25-7-9", "25-7-1"))
  expect_error(db_write_elisa(con, plate, assay_name = "anti-MPO", source_file = elisa_src_))
  expect_equal(nrow(read_log_(con)), 0)

  # fails inside the transaction: cytokine is NOT NULL
  plate <- elisa_plate_(cytokine = c("MPO", NA, "MPO", "MPO"))
  expect_error(db_write_elisa(con, plate, assay_name = "anti-MPO", source_file = elisa_src_))
  expect_equal(nrow(read_log_(con)), 0)
  expect_equal(nrow(DBI::dbGetQuery(con, "SELECT * FROM elisa_measurements")), 0)
})

test_that("db_write_elisa requires source_file", {
  con <- local_test_db()
  seed_two_mouse_plate_(con)

  expect_error(db_write_elisa(con, elisa_plate_(), assay_name = "anti-MPO"), "source_file")
  expect_error(
    db_write_elisa(con, elisa_plate_(), assay_name = "anti-MPO", source_file = NA_character_),
    "source_file"
  )
  expect_error(
    db_write_elisa(con, elisa_plate_(), assay_name = "anti-MPO", source_file = ""),
    "source_file"
  )
  expect_equal(nrow(DBI::dbGetQuery(con, "SELECT * FROM elisa_measurements")), 0)
})

histo_src_ <- "25-7/histo/scores.xlsx"

test_that("db_write_histo inserts measurements without an assay", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  histo_data <- tibble::tibble(
    metric = c("damage_score", "fibrosis_pct"), value = c(2, 15.5),
    mouse_id = "25-7-1", tissue = "spleen"
  )

  n <- db_write_histo(con, histo_data, source_file = histo_src_)
  expect_equal(n, 2)

  rows <- DBI::dbGetQuery(con, "SELECT metric, value, assay_id FROM histo_measurements ORDER BY metric")
  expect_equal(rows$metric, c("damage_score", "fibrosis_pct"))
  expect_true(all(is.na(rows$assay_id)))
})

test_that("db_write_histo errors clearly on an unknown assay_name", {
  con <- local_test_db()
  seed_minimal_fixture(con)

  histo_data <- tibble::tibble(
    metric = "damage_score", value = 2, mouse_id = "25-7-1", tissue = "spleen"
  )
  expect_error(
    db_write_histo(con, histo_data, assay_name = "not-registered", source_file = histo_src_),
    "No assay found"
  )
})

test_that("db_write_histo writes a whole table across mice", {
  con <- local_test_db()
  seed_two_mouse_plate_(con)

  histo_data <- tibble::tibble(
    metric = "damage_score", value = c(1, 3),
    mouse_id = c("25-7-1", "25-7-2"), tissue = "serum"
  )
  db_write_histo(con, histo_data, source_file = histo_src_)

  rows <- DBI::dbGetQuery(con, "
    SELECT sub.mouse_id, sm.tissue, h.value FROM histo_measurements h
    JOIN samples sm ON sm.sample_id = h.sample_id
    JOIN subjects sub ON sub.subject_id = sm.subject_id
    ORDER BY sub.mouse_id
  ")
  expect_equal(rows$mouse_id, c("25-7-1", "25-7-2"))
  expect_equal(rows$tissue, c("serum", "serum"))
  expect_equal(rows$value, c(1, 3))
})

test_that("db_write_histo lists every unresolved mouse_id x tissue and writes nothing", {
  con <- local_test_db()
  seed_two_mouse_plate_(con)

  histo_data <- tibble::tibble(
    metric = "damage_score", value = c(1, 2, 3),
    mouse_id = c("25-7-1", "25-7-9", NA), tissue = c("serum", "serum", "kidney")
  )
  err <- expect_error(db_write_histo(con, histo_data, source_file = histo_src_))
  expect_match(conditionMessage(err), "'25-7-9 x serum'", fixed = TRUE)
  expect_match(conditionMessage(err), "'NA x kidney'", fixed = TRUE)
  expect_equal(nrow(DBI::dbGetQuery(con, "SELECT * FROM histo_measurements")), 0)
  expect_equal(nrow(read_log_(con)), 0)
})

test_that("db_write_histo requires mouse_id and tissue columns", {
  con <- local_test_db()
  seed_two_mouse_plate_(con)

  histo_data <- tibble::tibble(
    metric = "damage_score", value = 1, mouse_id = "25-7-1", tissue = "serum"
  )
  expect_error(
    db_write_histo(con, dplyr::select(histo_data, !tissue), source_file = histo_src_),
    "tissue"
  )
  expect_error(
    db_write_histo(con, dplyr::rename(histo_data, mouse_ID = mouse_id), source_file = histo_src_),
    "mouse_ID"
  )
})

test_that("db_write_histo writes one import_log row per call", {
  con <- local_test_db()
  seed_two_mouse_plate_(con)

  histo_data <- tibble::tibble(
    metric = "damage_score", value = c(1, 3),
    mouse_id = c("25-7-1", "25-7-2"), tissue = "serum"
  )
  db_write_histo(con, histo_data, source_file = histo_src_)

  log <- read_log_(con)
  expect_equal(nrow(log), 1)
  expect_equal(log$domain_name, "histo")
  expect_equal(log$experiment_code, "25-7")
  expect_equal(log$import_unit, "scores.xlsx")
  expect_equal(log$rows_written, 2L)
  expect_true(is.na(log$qc_passed))
})

test_that("db_write_histo requires source_file", {
  con <- local_test_db()
  seed_two_mouse_plate_(con)

  histo_data <- tibble::tibble(
    metric = "damage_score", value = 1, mouse_id = "25-7-1", tissue = "serum"
  )
  expect_error(db_write_histo(con, histo_data), "source_file")
  expect_error(db_write_histo(con, histo_data, source_file = ""), "source_file")
})
