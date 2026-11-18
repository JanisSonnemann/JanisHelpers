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
