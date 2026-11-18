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
  db_write_subjects(con, tibble::tibble(mouse_id = "25-7-2", experiment_code = "25-7", end_date = as.Date("2026-06-10")))
  db_write_samples(con, tibble::tibble(mouse_id = "25-7-2", tissue = "spleen"))
  db_write_experiment(con, experiment_code = "25-8")
  db_write_subjects(con, tibble::tibble(mouse_id = "25-8-1", experiment_code = "25-8", end_date = as.Date("2026-06-10")))
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

test_that("lookup_sample_ids_ resolves NA to end_date and dated samples", {
  con <- local_test_db(); seed_minimal_fixture(con)
  db_write_samples(con, tibble::tibble(
    mouse_id = "25-7-1", tissue = "serum", collected_at = as.Date("2026-05-06")))
  db_write_samples(con, tibble::tibble(mouse_id = "25-7-1", tissue = "serum"))
  ids <- lookup_sample_ids_(con, c("25-7-1", "25-7-1"), c("serum", "serum"),
                            collected_at = as.Date(c(NA, "2026-05-06")))
  expect_length(unique(ids), 2)
  err <- expect_error(lookup_sample_ids_(con, "25-7-1", "serum",
                                         collected_at = as.Date("2026-05-07")))
  expect_match(conditionMessage(err), "25-7-1 x serum x 2026-05-07")
})

test_that("lookup_sample_ids_ keeps Date input as Date", {
  con <- local_test_db(); seed_minimal_fixture(con)
  db_write_samples(con, tibble::tibble(
    mouse_id = "25-7-1", tissue = "serum", collected_at = as.Date("2026-05-06")))
  expect_length(lookup_sample_ids_(con, c("25-7-1", "25-7-1"), c("serum", "serum"),
                                   collected_at = as.Date("2026-05-06")), 2)
})
