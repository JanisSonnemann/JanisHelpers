# Fixture: experiment 26-1 with three mice (MPO409, irrel, NA group), two FACS
# panels sharing the label CD4+ at different positions, one ELISA and one
# histology import unit; experiment 26-10 with one FACS panel only.
seed_load_fixture <- function(con) {
  db_write_experiment(con, experiment_code = "26-1")
  db_write_experiment(con, experiment_code = "26-10")
  db_write_subjects(con, tibble::tibble(
    mouse_id = c("26-1-1", "26-1-2", "26-1-3"), experiment_code = "26-1",
    treatment_group = c("MPO409", "irrel", NA), sex = c("f", "m", "f"),
    end_date = as.Date("2026-06-15")
  ))
  db_write_subjects(con, tibble::tibble(
    mouse_id = "26-10-1", experiment_code = "26-10",
    treatment_group = "irrel", end_date = as.Date("2026-07-01")
  ))
  db_write_samples(con, tibble::tibble(
    mouse_id = c("26-1-1", "26-1-2", "26-1-3", "26-10-1"), tissue = "spleen"
  ))
  db_write_samples(con, tibble::tibble(mouse_id = "26-1-1", tissue = "serum"))
  db_write_assay(con, assay_name = "overview", domain = "facs")
  db_write_assay(con, assay_name = "treg", domain = "facs")
  db_write_assay(con, assay_name = "anti-MPO", domain = "elisa")
  db_write_assay(con, assay_name = "scores", domain = "histo")

  facs_rows <- function(levels) {
    tibble::tibble(
      PopulationFullPath = paste0("root/", levels),
      Population = factor(levels, levels = levels),
      metric = "count", value = seq_along(levels)
    )
  }
  for (m in c("26-1-1", "26-1-2", "26-1-3")) {
    db_write_facs(con, facs_rows(c("live", "CD4+", "CD8+")),
      mouse_id = m, tissue = "spleen", assay_name = "overview",
      wsp_root = "overview-26-1", source_file = "26-1/overview-26-1_20260806.wsp")
    db_write_facs(con, facs_rows(c("CD4+", "FoxP3+", "live")),
      mouse_id = m, tissue = "spleen", assay_name = "treg",
      wsp_root = "Treg_26-1", source_file = "26-1/Treg_26-1_20260828.wsp")
  }
  db_write_facs(con, facs_rows(c("live", "CD4+")),
    mouse_id = "26-10-1", tissue = "spleen", assay_name = "overview",
    wsp_root = "overview-26-10", source_file = "26-10/overview-26-10_20260901.wsp")

  DBI::dbExecute(con, "INSERT INTO elisa_measurements
    (sample_id, assay_id, import_unit, analyte, replicate, value, unit, result_status)
    SELECT sm.sample_id, a.assay_id, 'anti-mpo-titer_26-1_long', 'anti-MPO', 1, 5, 'U', 'ok'
    FROM samples sm JOIN subjects s USING (subject_id), assays a
    WHERE a.assay_name = 'anti-MPO' AND s.mouse_id = '26-1-1' AND sm.tissue = 'serum'")
  DBI::dbExecute(con, "INSERT INTO histo_measurements
    (sample_id, assay_id, import_unit, metric, value)
    SELECT sm.sample_id, a.assay_id, 'histo_26-1', 'pct_crescents', 10
    FROM samples sm JOIN subjects s USING (subject_id), assays a
    WHERE a.assay_name = 'scores' AND s.experiment_id =
      (SELECT experiment_id FROM experiments WHERE experiment_code = '26-1')
      AND sm.tissue = 'spleen'")
  invisible(con)
}

local_load_db <- function(env = parent.frame()) {
  con <- local_test_db(env)
  seed_load_fixture(con)
  con
}

load_quiet_ <- function(...) suppressMessages(db_load_experiment(...))

test_that("returns meta + exclusions + one tibble per import unit in domain order", {
  con <- local_load_db()
  dat <- load_quiet_(con, "26-1")
  expect_equal(names(dat), c("meta", "cohort", "exclusions", "overview", "Treg", "anti-mpo-titer", "histo"))
  expect_equal(nrow(dat$meta), 3)
  for (el in dat) expect_s3_class(el, "tbl_df")
  expect_true(all(dat$overview$experiment_code == "26-1"))
  expect_equal(nrow(dat$overview), 9)
  expect_equal(nrow(dat$histo), 3)
})

test_that("automatic names strip the code literally and only at the end", {
  expect_equal(
    JanisHelpers:::unit_names_(c("overview-26-1", "Treg_26-1", "anti-mpo-titer_26-1_long",
                                 "MPO409-spleen-restim-48h_26-1"), "26-1"),
    c("overview", "Treg", "anti-mpo-titer", "MPO409-spleen-restim-48h")
  )
  expect_equal(JanisHelpers:::unit_names_("overview-26-10", "26-1"), "overview-26-10")
  expect_equal(JanisHelpers:::unit_names_("overview-2621", "26-1"), "overview-2621")
  con <- local_load_db()
  expect_equal(names(load_quiet_(con, "26-10")), c("meta", "cohort", "exclusions", "overview"))
})

test_that("names overrides selected units; unknown unit errors", {
  con <- local_load_db()
  dat <- load_quiet_(con, "26-1", names = c(titer = "anti-mpo-titer_26-1_long"))
  expect_equal(names(dat), c("meta", "cohort", "exclusions", "overview", "Treg", "titer", "histo"))
  expect_error(load_quiet_(con, "26-1", names = c(x = "nope")), "nope.*histo_26-1")
})

test_that("colliding, reserved or empty names error", {
  con <- local_load_db()
  expect_error(load_quiet_(con, "26-1", names = c(overview = "Treg_26-1")), "overview")
  expect_error(load_quiet_(con, "26-1", names = c(meta = "histo_26-1")), "meta")
  expect_error(load_quiet_(con, "26-1", names = c(exclusions = "histo_26-1")), "exclusions")
  expect_error(load_quiet_(con, "26-1", names = c(cohort = "histo_26-1")), "cohort")
  expect_error(load_quiet_(con, "26-1", names = stats::setNames("histo_26-1", "")), "name")
})

test_that("unknown experiment errors listing available codes", {
  con <- local_load_db()
  expect_error(db_load_experiment(con, "99-9"), "26-1")
})

test_that("group_levels sets the factor in meta and joined tibbles; NA group allowed", {
  con <- local_load_db()
  dat <- load_quiet_(con, "26-1", group_levels = c("MPO409", "irrel"))
  expect_equal(levels(dat$meta$group), c("MPO409", "irrel"))
  expect_equal(levels(dat$overview$group), c("MPO409", "irrel"))
  expect_true(is.na(dat$meta$group[dat$meta$mouse_id == "26-1-3"]))
  expect_error(load_quiet_(con, "26-1", group_levels = "MPO409"), "irrel")
})

test_that("subject_cols joins requested columns only", {
  con <- local_load_db()
  dat <- load_quiet_(con, "26-1")
  expect_true("group" %in% names(dat$histo))
  expect_false("sex" %in% names(dat$histo))
  dat2 <- load_quiet_(con, "26-1", subject_cols = c("group", "sex"))
  expect_true(all(c("group", "sex") %in% names(dat2$overview)))
  expect_error(load_quiet_(con, "26-1", subject_cols = "nope"), "nope")
})

test_that("FACS population factor is per panel", {
  con <- local_load_db()
  dat <- load_quiet_(con, "26-1")
  expect_equal(levels(dat$overview$population), c("live", "CD4+", "CD8+"))
  expect_equal(levels(dat$Treg$population), c("CD4+", "FoxP3+", "live"))
})

test_that("works on a read-only connection", {
  path <- tempfile(fileext = ".duckdb")
  con <- db_connect(path)
  seed_load_fixture(con)
  DBI::dbDisconnect(con, shutdown = TRUE)
  ro <- DBI::dbConnect(duckdb::duckdb(), path, read_only = TRUE)
  withr::defer(DBI::dbDisconnect(ro, shutdown = TRUE))
  dat <- load_quiet_(ro, "26-1")
  expect_equal(nrow(dat$meta), 3)
})

# Rule set on top of seed_load_fixture(): 26-1-1 excluded mouse-wide, 26-1-2
# only for the kidney tissue (matches no stored data, so nothing is flagged).
seed_exclusions_ <- function(con, rules = NULL) {
  if (is.null(rules)) {
    rules <- tibble::tibble(
      mouse_ID = c("26-1-1", "26-1-2"), tissue = c(NA, "kidney"),
      reason = "bad", decided_on = as.Date("2026-11-17")
    )
  }
  suppressMessages(db_write_exclusions(con, rules, "26-1"))
  invisible(con)
}

test_that("exclusions holds the current rule set", {
  con <- local_load_db(); seed_exclusions_(con)
  dat <- load_quiet_(con, "26-1")
  expect_named(dat$exclusions, c("exclusion_set_id", "mouse_id", "tissue", "panel",
    "condition", "file_name", "collected_at", "reason", "decided_on"))
  expect_equal(nrow(dat$exclusions), 2)
})

test_that("default load drops flagged rows and facility-ID dropouts", {
  con <- local_load_db(); seed_exclusions_(con)
  DBI::dbExecute(con, "UPDATE subjects SET id_source = 'facility' WHERE mouse_id = '26-1-3'")
  expect_message(dat <- db_load_experiment(con, "26-1"), "26-1-1")
  expect_equal(sort(dat$meta$mouse_id), "26-1-2")
  expect_false(any(c("id_source", "excluded_mouse_wide", "exclusion_reason") %in% names(dat$meta)))
  for (u in c("overview", "Treg", "histo")) {
    expect_false("26-1-1" %in% dat[[u]]$mouse_id)
    expect_false(any(c("excluded", "exclusion_reason") %in% names(dat[[u]])))
  }
  # only 26-1-2 remains (26-1-1 excluded, 26-1-3 facility ID)
  expect_equal(nrow(dat$overview), 3)
  expect_equal(nrow(dat$histo), 1)
  expect_equal(nrow(dat[["anti-mpo-titer"]]), 0)
})

test_that("keep_excluded keeps flagged rows with flags", {
  con <- local_load_db(); seed_exclusions_(con)
  dat <- load_quiet_(con, "26-1", keep_excluded = TRUE)
  expect_equal(nrow(dat$meta), 3)
  expect_true(all(c("id_source", "excluded_mouse_wide") %in% names(dat$meta)))
  expect_equal(nrow(dat$overview), 9)
  expect_true(all(c("excluded", "exclusion_reason") %in% names(dat$overview)))
  expect_equal(sum(dat$overview$excluded), 3)
  expect_true(all(dat$overview$mouse_id[dat$overview$excluded] == "26-1-1"))
  expect_true(dat[["anti-mpo-titer"]]$excluded)
})

test_that("cohort holds every subject including dropouts", {
  con <- local_load_db(); seed_exclusions_(con)
  DBI::dbExecute(con, "UPDATE subjects SET id_source = 'facility' WHERE mouse_id = '26-1-3'")
  dat <- load_quiet_(con, "26-1")
  expect_equal(nrow(dat$cohort), 3)
  expect_true(all(c("id_source", "excluded_mouse_wide") %in% names(dat$cohort)))
  expect_equal(nrow(dat$meta), 1)
})

test_that("no rule set gives a message; an empty set counts as recorded", {
  con <- local_load_db()
  expect_message(dat <- db_load_experiment(con, "26-1"),
                 "No exclusions recorded for experiment")
  expect_equal(nrow(dat$exclusions), 0)
  expect_equal(nrow(dat$meta), 3)
  expect_equal(nrow(dat$overview), 9)
  seed_exclusions_(con, tibble::tibble(mouse_ID = character(), reason = character()))
  expect_no_message(db_load_experiment(con, "26-1"))
})
