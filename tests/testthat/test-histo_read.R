test_that("histo_read_curated returns long metrics and keeps empty mice", {
  d <- suppressMessages(histo_read_curated(test_path("../fixtures/histo_curated.xlsx"), tissue = "kidney"))
  expect_named(d, c("mouse_ID", "tissue", "collected_at", "metric", "value"))
  expect_setequal(unique(d$metric), c("gloms", "crescents", "necrosis"))
  expect_true(all(is.na(d$value[d$mouse_ID == "M3"])))
  expect_equal(d$value[d$mouse_ID == "M1" & d$metric == "gloms"], 85)
})

test_that("histo_read_curated errors without a tissue", {
  expect_error(suppressMessages(histo_read_curated(test_path("../fixtures/histo_curated.xlsx"))),
               "No tissue")
})

test_that("histo_read_curated errors on non-numeric cells", {
  skip_if_not_installed("writexl")
  f <- withr::local_tempfile(fileext = ".xlsx")
  writexl::write_xlsx(list(values = data.frame(mouse_ID = "M1", gloms = "abc")), f)
  expect_error(suppressMessages(histo_read_curated(f, tissue = "kidney")),
               "gloms.*M1|abc")
})
