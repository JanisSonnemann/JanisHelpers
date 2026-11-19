test_that("histo_qc_rules_ returns the rows violating a rule", {
  d <- suppressMessages(histo_read_curated(
    test_path("..", "fixtures", "histo_curated.xlsx"), tissue = "kidney"))
  out <- histo_qc_rules_(d, "crescents <= gloms")
  expect_equal(out$mouse_ID, "M4")
  expect_equal(out$rule, "crescents <= gloms")
  expect_named(out, c("rule", "mouse_ID", "tissue", "timepoint", "values"))
  expect_equal(nrow(histo_qc_rules_(d, "necrosis <= gloms")), 0)
  expect_equal(nrow(histo_qc_rules_(d, character())), 0)
  expect_error(histo_qc_rules_(d, "foo <= gloms"), "foo <= gloms", fixed = TRUE)
})

test_that("histo_qc_comments_ lists one row per commented sample", {
  d <- tibble::tibble(
    mouse_ID = c("M1", "M1", "M2", "M2"), tissue = "kidney",
    collected_at = as.Date(NA), metric = c("gloms", "crescents", "gloms", "crescents"),
    value = 1, comment = c("torn", "torn", NA, NA)
  )
  ex <- tibble::tibble(mouse_ID = "M1", tissue = "kidney", timepoint = "end")
  out <- histo_qc_comments_(d, ex)
  expect_named(out, c("mouse_ID", "tissue", "timepoint", "comment", "excluded"))
  expect_equal(out$mouse_ID, "M1")
  expect_equal(out$comment, "torn")
  expect_equal(out$excluded, "excluded")
  expect_equal(nrow(histo_qc_comments_(dplyr::select(d, !comment), ex)), 0)
})

test_that("rule violations are marked excluded per mouse x tissue x timepoint", {
  d <- suppressMessages(histo_read_curated(
    test_path("..", "fixtures", "histo_curated.xlsx"), tissue = "kidney"))
  ends <- c(M4 = as.Date("2026-06-10"))
  viol <- histo_qc_rules_(qc_resolve_end_dates_(d, ends), "crescents <= gloms")
  expect_equal(viol$timepoint, "2026-06-10")
  ex <- tibble::tibble(mouse_ID = "M4", tissue = "kidney",
                       collected_at = as.Date(NA), reason = "x")
  expect_true(qc_is_excluded_(viol, qc_excluded_keys_(ex, ends)))
  ex$tissue <- "spleen"
  expect_false(qc_is_excluded_(viol, qc_excluded_keys_(ex, ends)))
  ex$tissue <- "kidney"
  ex$collected_at <- as.Date("2026-05-01")
  expect_false(qc_is_excluded_(viol, qc_excluded_keys_(ex, ends)))
})
