samples <- tibble::tibble(
  mouse_ID = c("m1", "m1", "m2", "m3"),
  tissue = c("spleen", "spleen", "spleen", "spleen"),
  collected_at = as.Date(c("2026-05-06", "2026-06-10", "2026-06-10", NA))
)

test_that("a mouse_ID-only row excludes all that mouse's samples", {
  ex <- tibble::tibble(mouse_ID = "m1", reason = "bad")
  res <- match_sample_exclusions_(samples, ex, "titer")
  expect_equal(nrow(res$excluded_samples), 2)
  expect_equal(res$applicable$n_samples_matched, 2L)
  expect_s3_class(res$applicable$collected_at, "Date")
  expect_equal(nrow(res$unmatched), 0)
})

test_that("a collected_at row excludes only that date", {
  ex <- tibble::tibble(mouse_ID = "m1", collected_at = "2026-05-06", reason = "bad")
  res <- match_sample_exclusions_(samples, ex, "titer")
  expect_equal(res$excluded_samples$collected_at, as.Date("2026-05-06"))
})

test_that("an exclusion date never matches an NA sample date", {
  ex <- tibble::tibble(mouse_ID = "m3", collected_at = "2026-06-10", reason = "bad")
  res <- match_sample_exclusions_(samples, ex, "titer")
  expect_equal(nrow(res$unmatched), 1)
})

test_that("rows for other panels are ignored; matching nothing is unmatched", {
  ex <- tibble::tibble(
    mouse_ID = c("m1", "zz"), panel = c("overview", NA), reason = "bad"
  )
  res <- match_sample_exclusions_(samples, ex, "titer")
  expect_equal(nrow(res$applicable), 1)
  expect_equal(res$unmatched$mouse_ID, "zz")
})

test_that("file_name / condition rows error, missing fields error", {
  ex <- tibble::tibble(mouse_ID = "m1", file_name = "a.fcs", reason = "bad")
  expect_error(match_sample_exclusions_(samples, ex, "titer"), "don't apply")
  ex2 <- tibble::tibble(mouse_ID = "m1", panel = "overview", file_name = "a", reason = "x")
  expect_no_error(match_sample_exclusions_(samples, ex2, "titer"))
  expect_error(
    match_sample_exclusions_(samples, tibble::tibble(mouse_ID = "m1"), "titer"),
    "needs column"
  )
  expect_error(
    match_sample_exclusions_(samples, tibble::tibble(mouse_ID = "m1", reason = ""), "titer"),
    "needs a mouse_ID and a reason"
  )
})

test_that("NULL exclusions give empty results", {
  res <- match_sample_exclusions_(samples, NULL, "titer")
  expect_equal(nrow(res$applicable), 0)
  expect_equal(nrow(res$excluded_samples), 0)
})

# ---- db_write_exclusions() ---------------------------------------------------

seed_exclusion_fixture <- function(con) {
  db_write_experiment(con, experiment_code = "E1")
  db_write_subjects(con, tibble::tibble(
    mouse_id = c("m1", "m2"), experiment_code = "E1",
    end_date = as.Date("2026-06-10")
  ))
  db_write_assay(con, assay_name = "overview", domain = "facs")
  invisible(con)
}
ex_sheet <- tibble::tibble(mouse_ID = "m1", reason = "r")

test_that("db_write_exclusions writes a first set and returns its id", {
  con <- local_test_db()
  seed_exclusion_fixture(con)
  id <- suppressMessages(db_write_exclusions(con, ex_sheet, "E1"))
  expect_identical(id, 1L)
  rules <- current_exclusion_rules_(con, "E1")
  expect_equal(nrow(rules), 1)
  expect_true(all(is.na(rules[c("tissue", "panel", "condition", "file_name", "collected_at")])))
})

test_that("an unchanged sheet writes no new set", {
  con <- local_test_db()
  seed_exclusion_fixture(con)
  suppressMessages(db_write_exclusions(con, ex_sheet, "E1"))
  expect_message(id <- db_write_exclusions(con, ex_sheet, "E1"), "unchanged")
  expect_identical(id, 1L)
  expect_equal(count_rows_(con, "exclusion_sets"), 1)
})

test_that("blank strings equal NA", {
  con <- local_test_db()
  seed_exclusion_fixture(con)
  suppressMessages(db_write_exclusions(con, ex_sheet, "E1"))
  expect_message(db_write_exclusions(con, dplyr::mutate(ex_sheet, tissue = ""), "E1"),
                 "unchanged")
})

test_that("an empty sheet writes an empty set", {
  con <- local_test_db()
  seed_exclusion_fixture(con)
  suppressMessages(db_write_exclusions(con, ex_sheet, "E1"))
  expect_identical(db_write_exclusions(con, ex_sheet[0, ], "E1"), 2L)
  expect_equal(nrow(current_exclusion_rules_(con, "E1")), 0)
  expect_message(db_write_exclusions(con, NULL, "E1"), "unchanged")
  expect_equal(count_rows_(con, "exclusion_sets"), 2)
  suppressMessages(db_write_exclusions(con, ex_sheet, "E1"))
  expect_identical(db_write_exclusions(con, NULL, "E1"), 4L)
  expect_equal(nrow(current_exclusion_rules_(con, "E1")), 0)
})

test_that("db_write_exclusions validates and writes nothing on error", {
  con <- local_test_db()
  seed_exclusion_fixture(con)
  expect_error(db_write_exclusions(con, tibble::tibble(mouse_ID = "m1"), "E1"), "reason")
  expect_error(
    db_write_exclusions(con, tibble::tibble(mouse_ID = c("m1", NA), reason = c("r", "r")), "E1"),
    "mouse_ID")
  expect_error(
    db_write_exclusions(con, tibble::tibble(mouse_ID = "m1", reason = ""), "E1"),
    "reason")
  expect_error(
    db_write_exclusions(con, tibble::tibble(mouse_ID = "zz", reason = "r"), "E1"), "zz")
  expect_error(db_write_exclusions(con, ex_sheet, "nope"), "nope")
  expect_equal(count_rows_(con, "exclusion_sets"), 0)
  expect_equal(count_rows_(con, "exclusion_rules"), 0)
})

test_that("an unknown panel only messages", {
  con <- local_test_db()
  seed_exclusion_fixture(con)
  expect_message(
    db_write_exclusions(con, dplyr::mutate(ex_sheet, panel = "elisaX"), "E1"),
    "not registered yet")
  expect_equal(count_rows_(con, "exclusion_sets"), 1)
})

test_that("rules for mice without samples are reported", {
  con <- local_test_db()
  seed_exclusion_fixture(con)
  expect_message(db_write_exclusions(con, ex_sheet, "E1"), "match no stored row.*m1|m1.*match no stored row")
})

seed_rule_match_fixture <- function(con) {
  seed_exclusion_fixture(con)
  db_write_samples(con, tibble::tibble(mouse_id = "m1", tissue = "spleen"))
  db_write_facs(
    con,
    tibble::tibble(PopulationFullPath = "root/CD4+", Population = "CD4+",
                   metric = "count", value = 1),
    mouse_id = "m1", tissue = "spleen", assay_name = "overview",
    wsp_root = "overview_E1", source_file = "E1/overview_E1_20261101.wsp"
  )
  DBI::dbExecute(con, "UPDATE facs_stains SET file_name = 'm1.fcs'")
  invisible(con)
}

test_that("rules_without_match_ uses the full match, not just the mouse", {
  con <- local_test_db()
  seed_rule_match_fixture(con)
  sheet <- tibble::tibble(
    mouse_ID = c("m1", "m1", "m1", "m1"),
    tissue = c("spleen", "kidney", NA, NA),
    file_name = c(NA, NA, "m1.fcs", "other.fcs"),
    reason = c("ok", "wrong tissue", "ok file", "wrong file"))
  suppressMessages(db_write_exclusions(con, sheet, "E1"))
  res <- rules_without_match_(con, "E1")
  expect_setequal(res$reason, c("wrong tissue", "wrong file"))
  expect_false(any(res$needs_reexport))
  expect_true(all(c("mouse_id", "tissue", "panel", "condition", "file_name",
                    "collected_at", "reason", "needs_reexport") %in% names(res)))
})

test_that("file rules on stains without file_name need a re-export", {
  con <- local_test_db()
  seed_rule_match_fixture(con)
  DBI::dbExecute(con, "UPDATE facs_stains SET file_name = NULL")
  sheet <- tibble::tibble(mouse_ID = c("m1", "m1"), file_name = c("m1.fcs", NA),
                          tissue = c(NA, "kidney"), reason = c("file", "kid"))
  expect_message(db_write_exclusions(con, sheet, "E1"), "re-export")
  res <- rules_without_match_(con, "E1") |> dplyr::arrange(reason)
  expect_equal(res$reason, c("file", "kid"))
  expect_equal(res$needs_reexport, c(TRUE, FALSE))
})

test_that("rules_without_match_ is empty without rules or set", {
  con <- local_test_db()
  seed_rule_match_fixture(con)
  expect_equal(nrow(rules_without_match_(con, "E1")), 0)
  suppressMessages(db_write_exclusions(con, NULL, "E1"))
  expect_equal(nrow(rules_without_match_(con, "E1")), 0)
})
