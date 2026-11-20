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
    "result_status", "lloq", "uloq", "source_file", "imported_at", "excluded",
    "exclusion_reason"))
  expect_equal(elisa$analyte, "MPO")
  expect_equal(elisa$assay_name, "anti-MPO")
  histo <- dplyr::collect(db_query_histo(con))
  expect_named(histo, c("experiment_code", "mouse_id", "tissue", "collected_at",
    "assay_name", "import_unit", "metric", "value", "source_file", "imported_at",
    "excluded", "exclusion_reason"))
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
    "start_date", "bmt_date", "end_date", "mdc_id", "early", "early_reason",
    "id_source", "excluded_mouse_wide", "exclusion_reason"))
  expect_equal(nrow(all), 2)
  expect_false("treatment_group" %in% names(all))
})

test_that("db_query_subjects returns dropout columns", {
  con <- local_test_db()
  db_write_experiment(con, "25-7")
  db_write_subjects(con, tibble::tibble(
    mouse_id = "F2", experiment_code = "25-7", mdc_id = "F2", early = TRUE,
    early_reason = "found dead", id_source = "facility"))
  got <- dplyr::collect(db_query_subjects(con))
  expect_true(all(c("mdc_id", "early", "early_reason", "id_source") %in% names(got)))
  expect_equal(got$early, TRUE)
  expect_equal(got$early_reason, "found dead")
  expect_equal(got$id_source, "facility")
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

# ---- Exclusion flags on read -------------------------------------------------

# One experiment "26-1" (mice M1-M3, end date 2026-06-10) with a FACS panel
# `overview` (two wells per spleen via keyword `restim`, plus kidney) and the
# ELISA fixture (assay `restim`: M1 serial + end point, M2 end point), with an
# M3 end point copied from M2. Returns the FACS and ELISA preps.
flag_files <- tibble::tribble(
  ~file_name,     ~mouse_ID, ~tissue,   ~restim,
  "M1_spl_0.fcs", "M1",      "spleen",  "none",
  "M1_spl_s.fcs", "M1",      "spleen",  "stim",
  "M1_kid_0.fcs", "M1",      "kidney",  "none",
  "M2_spl_0.fcs", "M2",      "spleen",  "none",
  "M2_spl_s.fcs", "M2",      "spleen",  "stim",
  "M2_kid_0.fcs", "M2",      "kidney",  "none",
  "M3_spl_0.fcs", "M3",      "spleen",  "none",
  "M3_spl_s.fcs", "M3",      "spleen",  "stim",
  "M3_kid_0.fcs", "M3",      "kidney",  "none"
)

seed_flag_fixture <- function(con, env = parent.frame()) {
  db_write_experiment(con, experiment_code = "26-1")
  db_write_subjects(con, tibble::tibble(
    mouse_id = c("M1", "M2", "M3"), experiment_code = "26-1",
    end_date = as.Date("2026-06-10")
  ))
  d <- local_prep_dirs(env)
  facs <- suppressMessages(
    local_loaded_prep(d, flag_files, condition_keyword = "restim")
  )
  suppressMessages(db_export_facs(con, facs, data_dir = d$data_dir))

  folder <- file.path(d$data_dir, "26-1_exp", "raw", "ELISA")
  dir.create(folder, recursive = TRUE)
  path <- file.path(folder, "elisa_curated.xlsx")
  file.copy(test_path("..", "fixtures", "elisa_curated.xlsx"), path)
  Sys.setFileTime(path, Sys.time() - 3600)
  data <- suppressMessages(elisa_read_curated(path, tissue = "spleen"))
  data <- dplyr::bind_rows(
    data, dplyr::mutate(dplyr::filter(data, mouse_ID == "M2"), mouse_ID = "M3")
  )
  suppressMessages(elisa_prep_save(data, path, "restim"))
  paths <- elisa_prep_paths_(path)
  dir.create(dirname(paths$qc))
  writeLines("<html></html>", paths$qc)
  Sys.setFileTime(paths$qc, Sys.time() + 10)
  elisa <- elisa_prep_load(path)
  suppressMessages(db_export_elisa(con, elisa, data_dir = d$data_dir))
  list(facs = facs, elisa = elisa)
}

flag_sheet <- tibble::tribble(
  ~mouse_ID, ~tissue,  ~panel,     ~condition, ~file_name,     ~collected_at,         ~reason,
  "M2",      NA,       NA,         NA,         NA,             NA,                    "r_mouse",
  "M1",      "kidney", NA,         NA,         NA,             NA,                    "r_tissue",
  "M3",      NA,       "overview", "stim",     NA,             NA,                    "r_cond",
  "M3",      NA,       "overview", NA,         "M3_kid_0.fcs", NA,                    "r_file",
  "M1",      NA,       NA,         NA,         NA,             as.Date("2026-05-06"), "r_serial",
  "M1",      NA,       "restim",   NA,         NA,             as.Date("2026-06-10"), "r_end"
)

test_that("exclusion flags match the R matchers", {
  con <- local_test_db()
  preps <- seed_flag_fixture(con)
  suppressMessages(db_write_exclusions(con, flag_sheet, "26-1"))

  facs <- dplyr::collect(db_query_facs(con))
  expect_type(facs$excluded, "logical")
  stains <- DBI::dbGetQuery(con, "
    SELECT fs.file_name, su.mouse_id, sa.tissue, fs.condition
    FROM facs_stains fs JOIN samples sa USING (sample_id)
    JOIN subjects su USING (subject_id)")
  flagged_keys <- dplyr::distinct(dplyr::filter(facs, excluded),
                                  mouse_id, tissue, condition)
  flagged_files <- sort(dplyr::semi_join(
    stains, flagged_keys, by = c("mouse_id", "tissue", "condition"))$file_name)
  r_files <- facs_match_exclusions(preps$facs$data, flag_sheet,
                                   "overview")$excluded_files$file_name
  expect_equal(flagged_files, sort(r_files))
  expect_equal(flagged_files, sort(c("M1_kid_0.fcs", "M2_spl_0.fcs",
    "M2_spl_s.fcs", "M2_kid_0.fcs", "M3_spl_s.fcs", "M3_kid_0.fcs")))
  # a stain is flagged in all its measurements or none
  per_stain <- dplyr::summarise(facs, n = dplyr::n_distinct(excluded),
                                .by = c(mouse_id, tissue, condition))
  expect_true(all(per_stain$n == 1))

  elisa <- dplyr::collect(db_query_elisa(con))
  flagged <- dplyr::distinct(dplyr::filter(elisa, excluded),
                             mouse_ID = mouse_id, tissue, collected_at) |>
    dplyr::arrange(mouse_ID, collected_at)
  resolved <- preps$elisa$data
  empty <- is.na(resolved$collected_at)
  resolved$collected_at[empty] <- as.Date("2026-06-10")
  r_samples <- match_sample_exclusions_(resolved, flag_sheet, "restim")$excluded_samples |>
    dplyr::distinct(mouse_ID, tissue, collected_at) |>
    dplyr::arrange(mouse_ID, collected_at)
  expect_equal(flagged, r_samples)
  # the end-date rule (r_end) flags M1's end-point sample
  expect_equal(
    unique(elisa$exclusion_reason[elisa$mouse_id == "M1" &
                                  elisa$collected_at == as.Date("2026-06-10")]),
    "r_end")
  expect_false(any(elisa$excluded[elisa$mouse_id == "M3"]))
  expect_false(anyNA(elisa$excluded))
})

test_that("multiple rules give one row per measurement", {
  con <- local_test_db()
  seed_flag_fixture(con)
  n_before <- nrow(dplyr::collect(db_query_facs(con)))
  sheet <- tibble::tibble(mouse_ID = "M1", tissue = c("spleen", NA),
                          reason = c("r_tissue", "r_mouse"))
  suppressMessages(db_write_exclusions(con, sheet, "26-1"))
  facs <- dplyr::collect(db_query_facs(con))
  expect_equal(nrow(facs), n_before)
  m1 <- dplyr::filter(facs, mouse_id == "M1")
  expect_true(all(m1$excluded))
  expect_equal(unique(m1$exclusion_reason[m1$tissue == "spleen"]), "r_mouse; r_tissue")
  expect_equal(unique(m1$exclusion_reason[m1$tissue == "kidney"]), "r_mouse")
  expect_true(all(is.na(facs$exclusion_reason[!facs$excluded])))
})

test_that("no set means nothing excluded", {
  con <- local_test_db()
  seed_flag_fixture(con)
  db_write_assay(con, "scores", "histo")
  DBI::dbExecute(con, "INSERT INTO histo_measurements
    (sample_id, assay_id, import_unit, metric, value)
    SELECT sm.sample_id, a.assay_id, 'scores', 'damage_score', 2
    FROM samples sm, assays a WHERE a.assay_name = 'scores' LIMIT 1")
  check_none <- function() {
    for (q in list(db_query_facs, db_query_elisa, db_query_histo)) {
      got <- dplyr::collect(q(con))
      expect_true(nrow(got) > 0)
      expect_identical(unique(got$excluded), FALSE)
      expect_true(all(is.na(got$exclusion_reason)))
    }
    subj <- dplyr::collect(db_query_subjects(con))
    expect_identical(unique(subj$excluded_mouse_wide), FALSE)
  }
  check_none()
  suppressMessages(db_write_exclusions(con, NULL, "26-1"))
  expect_equal(count_rows_(con, "exclusion_sets"), 1)
  check_none()
})

test_that("histology rows are flagged; condition / file_name rules never match", {
  con <- local_test_db()
  seed_flag_fixture(con)
  db_write_assay(con, "scores", "histo")
  DBI::dbExecute(con, "INSERT INTO histo_measurements
    (sample_id, assay_id, import_unit, metric, value)
    SELECT sm.sample_id, a.assay_id, 'scores', 'damage_score', 2
    FROM samples sm JOIN subjects su USING (subject_id), assays a
    WHERE a.assay_name = 'scores' AND sm.tissue = 'kidney'")
  sheet <- tibble::tibble(
    mouse_ID = c("M1", "M2", "M3"), tissue = "kidney",
    condition = c(NA, "none", NA), file_name = c(NA, NA, "M3_kid_0.fcs"),
    reason = c("r1", "r2", "r3"))
  suppressMessages(db_write_exclusions(con, sheet, "26-1"))
  histo <- dplyr::collect(db_query_histo(con)) |> dplyr::arrange(mouse_id)
  expect_equal(histo$excluded, c(TRUE, FALSE, FALSE))
  expect_equal(histo$exclusion_reason, c("r1", NA, NA))
})

test_that("rules do not cross experiments", {
  con <- local_test_db()
  db_write_experiment(con, "E1")
  db_write_experiment(con, "E2")
  db_write_subjects(con, tibble::tibble(
    mouse_id = c("a1", "b1"), experiment_code = c("E1", "E2"),
    end_date = as.Date("2026-06-10")))
  db_write_samples(con, tibble::tibble(mouse_id = c("a1", "b1"), tissue = "spleen"))
  db_write_assay(con, "overview", "facs")
  for (m in c("a1", "b1")) {
    db_write_facs(
      con,
      tibble::tibble(PopulationFullPath = "root/CD4+", Population = "CD4+",
                     metric = "count", value = 1),
      mouse_id = m, tissue = "spleen", assay_name = "overview",
      wsp_root = paste0("overview_", m), source_file = paste0("overview_", m, "_20261101.wsp")
    )
  }
  suppressMessages(db_write_exclusions(
    con, tibble::tibble(mouse_ID = "a1", tissue = "spleen", panel = "overview",
                        reason = "e1"), "E1"))
  suppressMessages(db_write_exclusions(
    con, tibble::tibble(mouse_ID = "b1", tissue = "kidney", reason = "e2"), "E2"))
  # A rule naming E1's mouse inside E2's set must not reach E1's rows.
  DBI::dbExecute(con, "INSERT INTO exclusion_rules (exclusion_set_id, mouse_id, reason)
    SELECT MAX(s.exclusion_set_id), 'a1', 'leak' FROM exclusion_sets s
    JOIN experiments e USING (experiment_id) WHERE e.experiment_code = 'E2'")
  facs <- dplyr::collect(db_query_facs(con)) |> dplyr::arrange(mouse_id)
  expect_equal(facs$excluded, c(TRUE, FALSE))
  expect_equal(facs$exclusion_reason, c("e1", NA))
  subj <- dplyr::collect(db_query_subjects(con)) |> dplyr::arrange(mouse_id)
  expect_equal(subj$excluded_mouse_wide, c(FALSE, FALSE))
})

test_that("old stains without file_name are not flagged by file rules", {
  con <- local_test_db()
  seed_flag_fixture(con)
  DBI::dbExecute(con, "UPDATE facs_stains SET file_name = NULL")
  sheet <- tibble::tibble(mouse_ID = "M3", file_name = "M3_kid_0.fcs",
                          reason = "r_file")
  expect_message(db_write_exclusions(con, sheet, "26-1"),
                 "re-export[^\n]*\n.*M3_kid_0.fcs")
  facs <- dplyr::collect(db_query_facs(con))
  expect_false(any(facs$excluded))
})

test_that("changing the sheet changes the flags without re-export", {
  con <- local_test_db()
  seed_flag_fixture(con)
  suppressMessages(db_write_exclusions(
    con, tibble::tibble(mouse_ID = "M1", reason = "r"), "26-1"))
  expect_true(any(dplyr::collect(db_query_facs(con))$excluded))
  suppressMessages(db_write_exclusions(
    con, tibble::tibble(mouse_ID = "M2", tissue = "liver", reason = "r"), "26-1"))
  expect_false(any(dplyr::collect(db_query_facs(con))$excluded))
  expect_false(any(dplyr::collect(db_query_elisa(con))$excluded))
  expect_equal(count_rows_(con, "exclusion_sets"), 2)
})

test_that("db_query_subjects flags mouse-wide exclusions", {
  con <- local_test_db()
  seed_flag_fixture(con)
  suppressMessages(db_write_exclusions(con, flag_sheet, "26-1"))
  subj <- dplyr::collect(db_query_subjects(con)) |> dplyr::arrange(mouse_id)
  expect_equal(subj$excluded_mouse_wide, c(FALSE, TRUE, FALSE))
  expect_equal(subj$exclusion_reason, c(NA, "r_mouse", NA))
  expect_equal(tail(names(subj), 2), c("excluded_mouse_wide", "exclusion_reason"))
})

test_that("db_query_* fall back without flag views (old file, read-only)", {
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
  for (v in c("facs", "elisa", "histo")) {
    DBI::dbExecute(con, paste0("DROP VIEW IF EXISTS ", v, "_exclusion_flags"))
  }
  DBI::dbExecute(con, "DROP TABLE exclusion_rules")
  DBI::dbExecute(con, "DROP TABLE exclusion_sets")
  DBI::dbDisconnect(con, shutdown = TRUE)

  ro <- DBI::dbConnect(duckdb::duckdb(), path, read_only = TRUE)
  withr::defer(DBI::dbDisconnect(ro, shutdown = TRUE))
  facs <- dplyr::collect(db_query_facs(ro))
  expect_identical(facs$excluded, FALSE)
  expect_true(is.na(facs$exclusion_reason))
  expect_identical(dplyr::collect(db_query_elisa(ro))$excluded, logical())
  subj <- dplyr::collect(db_query_subjects(ro))
  expect_identical(subj$excluded_mouse_wide, FALSE)
})
