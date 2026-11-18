library(testthat)
library(JanisHelpers)

il17a_path <- testthat::test_path("../fixtures/Results_IL-17A.xlsx")
tnfa_path  <- testthat::test_path("../fixtures/Results_TNFa.xlsx")
il17a_skip_msg <- "Results_IL-17A.xlsx fixture not available"
tnfa_skip_msg  <- "Results_TNFa.xlsx fixture not available"

test_that("elisa_read_results() returns the expected columns and types", {
  skip_if_not(file.exists(il17a_path), il17a_skip_msg)
  dat <- elisa_read_results(il17a_path)
  expect_s3_class(dat, "tbl_df")
  expect_equal(
    names(dat),
    c("cytokine", "sample_id", "replicate", "value", "unit", "result_status")
  )
  expect_type(dat$cytokine, "character")
  expect_type(dat$sample_id, "character")
  expect_type(dat$replicate, "integer")
  expect_type(dat$value, "double")
  expect_type(dat$unit, "character")
  expect_type(dat$result_status, "character")
  expect_equal(nrow(dat), 56L)
})

test_that("elisa_read_results() derives cytokine from filename (no cytokine-prefixed sheet)", {
  skip_if_not(file.exists(il17a_path), il17a_skip_msg)
  dat <- elisa_read_results(il17a_path)
  expect_equal(unique(dat$cytokine), "IL-17A")
})

test_that("elisa_read_results() derives cytokine from filename (cytokine-prefixed sheet)", {
  skip_if_not(file.exists(tnfa_path), tnfa_skip_msg)
  dat <- elisa_read_results(tnfa_path)
  expect_equal(unique(dat$cytokine), "TNFa")
  expect_equal(nrow(dat), 56L)
})

test_that("elisa_read_results() parses the unit from the Backcalc column header", {
  skip_if_not(file.exists(il17a_path), il17a_skip_msg)
  dat <- elisa_read_results(il17a_path)
  expect_equal(unique(dat$unit), "pg/ml")
})

test_that("elisa_read_results() sets value to NA for OOR rows and keeps result_status", {
  skip_if_not(file.exists(il17a_path), il17a_skip_msg)
  dat <- elisa_read_results(il17a_path)
  oor <- dplyr::filter(dat, result_status == "OOR<")
  expect_equal(nrow(oor), 20L)
  expect_true(all(is.na(oor$value)))
  expect_setequal(
    unique(oor$sample_id),
    c("25-7-5", "25-7-11", "25-7-12", "25-7-13", "25-7-14",
      "25-7-20", "25-7-21", "25-7-24", "25-7-25", "25-7-28")
  )
})

test_that("elisa_read_results() keeps a non-NA value for <LLOQ rows (below quantification, not out of range)", {
  skip_if_not(file.exists(tnfa_path), tnfa_skip_msg)
  dat <- elisa_read_results(tnfa_path)
  lloq <- dplyr::filter(dat, result_status == "<LLOQ")
  expect_equal(nrow(lloq), 11L)
  expect_true(all(!is.na(lloq$value)))
  row <- dplyr::filter(lloq, sample_id == "25-7-3", replicate == 1L)
  expect_equal(row$value, 1.303633, tolerance = 1e-6)
})

test_that("elisa_read_results() uses the cytokine argument when supplied, overriding the filename", {
  skip_if_not(file.exists(tnfa_path), tnfa_skip_msg)
  dat <- elisa_read_results(tnfa_path, cytokine = "Custom")
  expect_equal(unique(dat$cytokine), "Custom")
})

# ---- elisa_read_curated() ---------------------------------------------------

curated_path <- testthat::test_path("../fixtures/elisa_curated.xlsx")

write_curated_ <- function(values, range = NULL) {
  testthat::skip_if_not_installed("writexl")
  f <- tempfile(fileext = ".xlsx")
  sheets <- list(values = values)
  if (!is.null(range)) sheets$range <- range
  writexl::write_xlsx(sheets, f)
  f
}

test_that("elisa_read_curated returns long rows with replicates and range", {
  d <- suppressMessages(elisa_read_curated(curated_path, tissue = "spleen"))
  expect_named(d, c("mouse_ID", "tissue", "collected_at", "analyte", "replicate", "value", "result_status", "unit", "lloq", "uloq"))
  m1 <- dplyr::filter(d, mouse_ID == "M1", is.na(collected_at), analyte == "IL-2")
  expect_equal(m1$replicate, 1:2)
  expect_equal(m1$result_status, c("OOR<", "ok"))
  expect_true(is.na(m1$value[1]))
  expect_equal(unique(d$lloq[d$analyte == "TNFa"]), 7.81)
  expect_equal(unique(d$unit), "pg/ml")
  expect_equal(unique(d$tissue), "spleen")
  expect_equal(sum(!is.na(d$collected_at)), 2)
  expect_equal(d$replicate[d$analyte == "TNFa"], c(1L, 2L, 1L, 2L, 1L))
})

test_that("elisa_read_curated: tissue column overrides the argument where non-empty", {
  f <- write_curated_(data.frame(
    mouse_ID = c("M1", "M2"), tissue = c("serum", NA), A = c(1, 2)
  ))
  d <- suppressMessages(elisa_read_curated(f, tissue = "spleen"))
  expect_equal(d$tissue, c("serum", "spleen"))
})

test_that("elisa_read_curated errors without any tissue", {
  f <- write_curated_(data.frame(mouse_ID = "M1", A = 1))
  expect_error(suppressMessages(elisa_read_curated(f)), "tissue")
})

test_that("elisa_read_curated without range sheet gives NA lloq/uloq/unit", {
  f <- write_curated_(data.frame(mouse_ID = "M1", A = 1))
  d <- suppressMessages(elisa_read_curated(f, tissue = "serum"))
  expect_true(all(is.na(d$lloq)))
  expect_true(all(is.na(d$uloq)))
  expect_true(all(is.na(d$unit)))
})

test_that("elisa_read_curated errors on unreadable analyte cells", {
  f <- write_curated_(data.frame(mouse_ID = c("M1", "M7"), TNFa = c("1.5", "n.d.")))
  expect_error(suppressMessages(elisa_read_curated(f, tissue = "serum")), "n\\.d\\.")
  expect_error(suppressMessages(elisa_read_curated(f, tissue = "serum")), "TNFa")
  expect_error(suppressMessages(elisa_read_curated(f, tissue = "serum")), "M7")
})

test_that("parse_curated_dates_ handles the supported types", {
  x <- list(as.POSIXct("2026-05-06", tz = "UTC"), "2026-05-20", 46147, NA)
  r <- parse_curated_dates_(x, paste0("M", 1:4))
  expect_s3_class(r, "Date")
  expect_equal(r, as.Date(c("2026-05-06", "2026-05-20", "2026-05-05", NA)))
  expect_error(parse_curated_dates_(list("15/16.6.2026"), "M9"), "M9")
})
