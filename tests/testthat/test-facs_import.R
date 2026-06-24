test_that("facs_import_wsp() emits a deprecation warning", {
  skip_if_not(
    file.exists(testthat::test_path("../fixtures/minimal.wsp")),
    "WSP fixture not available"
  )
  expect_warning(
    suppressMessages(
      facs_import_wsp(testthat::test_path("../fixtures/minimal.wsp"))
    ),
    regexp = "deprecated",
    ignore.case = TRUE
  )
})
