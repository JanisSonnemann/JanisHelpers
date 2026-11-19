test_that("histo_read_curated returns long metrics and keeps empty mice", {
  d <- suppressMessages(histo_read_curated(test_path("../fixtures/histo_curated.xlsx"), tissue = "kidney"))
  expect_named(d, c("mouse_ID", "tissue", "collected_at", "metric", "value"))
  expect_setequal(unique(d$metric), c("gloms", "crescents", "necrosis"))
  expect_true(all(is.na(d$value[d$mouse_ID == "M3"])))
  expect_equal(d$value[d$mouse_ID == "M1" & d$metric == "gloms"], 85)
})

test_that("histo_read_curated reads a comment column as text", {
  skip_if_not_installed("writexl")
  f <- withr::local_tempfile(fileext = ".xlsx")
  writexl::write_xlsx(list(values = data.frame(
    mouse_ID = c("M1", "M2"), gloms = c(80, 70), Comment = c("torn section", NA)
  )), f)
  d <- suppressMessages(histo_read_curated(f, tissue = "kidney"))
  expect_named(d, c("mouse_ID", "tissue", "collected_at", "metric", "value", "comment"))
  expect_equal(unique(d$metric), "gloms")
  expect_equal(d$comment[d$mouse_ID == "M1"], "torn section")
  expect_true(is.na(d$comment[d$mouse_ID == "M2"]))

  writexl::write_xlsx(list(values = data.frame(mouse_ID = "M1", gloms = 80, comments = 5)), f)
  d <- suppressMessages(histo_read_curated(f, tissue = "kidney"))
  expect_equal(unique(d$metric), "gloms")
  expect_equal(d$comment, "5")

  writexl::write_xlsx(list(values = data.frame(mouse_ID = "M1", gloms = 80,
                                               comment = "a", COMMENTS = "b")), f)
  expect_error(suppressMessages(histo_read_curated(f, tissue = "kidney")),
               "comment column")
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

test_that("histo_read_curated errors on empty or duplicated metric headers", {
  skip_if_not_installed("writexl")
  f <- withr::local_tempfile(fileext = ".xlsx")
  v <- data.frame(mouse_ID = "M1", gloms = 1, gloms = 2, check.names = FALSE)
  writexl::write_xlsx(list(values = v), f)
  expect_error(suppressMessages(histo_read_curated(f, tissue = "kidney")),
               "duplicated metric column\\(s\\): gloms")
  v <- data.frame(mouse_ID = "M1", gloms = 1, x = 2)
  names(v)[3] <- ""
  writexl::write_xlsx(list(values = v), f)
  expect_error(suppressMessages(histo_read_curated(f, tissue = "kidney")),
               "without a header")
})

test_that("the shipped histology template reads once filled in", {
  skip_if_not_installed("writexl")
  tpl <- system.file("templates", "histo_curated_template.xlsx", package = "JanisHelpers")
  expect_true(nzchar(tpl))
  expect_equal(readxl::excel_sheets(tpl), c("values", "README"))

  values <- readxl::read_excel(tpl, "values")
  expect_equal(names(values)[1:4], c("mouse_ID", "ID", "tissue", "date"))
  filled <- tibble::tibble(
    mouse_ID = c("M1", "M2"), ID = c("S1", "S2"), tissue = NA_character_,
    date = as.Date(NA), gloms = c(80, NA), crescents = c(8, NA),
    necrosis = c(2, NA)
  )
  path <- tempfile(fileext = ".xlsx")
  writexl::write_xlsx(list(values = filled,
                           README = readxl::read_excel(tpl, "README")), path)

  d <- suppressMessages(histo_read_curated(path, tissue = "kidney"))
  expect_setequal(unique(d$metric), c("gloms", "crescents", "necrosis"))
  expect_equal(d$value[d$mouse_ID == "M1" & d$metric == "gloms"], 80)
  expect_true(all(is.na(d$value[d$mouse_ID == "M2"])))
})
