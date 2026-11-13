qc_files <- tibble::tibble(
  file_name = c("a.fcs", "b.fcs"),
  mouse_ID = c("25-7-1", "25-7-2"),
  tissue = "spleen"
)

# -- settings -------------------------------------------------------------------

test_that("qc_settings_ fills defaults for a prep saved without qc settings", {
  prep <- list(data = fake_prep_data(qc_files), info = list())
  s <- qc_settings_(prep)
  expect_null(s$group_col)
  expect_null(s$batch_col)
  expect_null(s$expected_mice)
  expect_equal(s$expected_tissues, "spleen")
  expect_equal(s$qc_populations, c("live_CD45+", "CD4+"))
  expect_equal(s$qc_mfi_populations, character())
  expect_equal(s$qc_metrics, character())
})

test_that("qc_settings_ returns stored settings when present", {
  stored <- list(group_col = "grp", batch_col = NULL, expected_mice = NULL,
                 expected_tissues = "spleen", qc_populations = "CD4+",
                 qc_mfi_populations = character(), qc_metrics = character())
  prep <- list(data = fake_prep_data(qc_files), info = list(qc = stored))
  expect_equal(qc_settings_(prep), stored)
})
