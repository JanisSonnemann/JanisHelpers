wide_raw_ <- function() {
  tibble::tibble(
    `Lab ID` = c("26-1-1", "26-1-2"), ID = c("21560", "21561"),
    group = c("irrel", "MPO409"), `46147` = c("2.6658", "2.6729"),
    `2026-05-19` = c("2.6712", NA), final = c("2.5228", "OOR<")
  )
}

test_that("single analyte: date headers become rows, final is end-point, extras dropped", {
  r <- elisa_reshape_wide_table_(wide_raw_(), analyte = "anti-MPO-titer")
  expect_named(r$values, c("mouse_ID", "date", "anti-MPO-titer"))
  expect_equal(r$values$date, as.Date(c("2026-05-05", "2026-05-05", "2026-05-19",
                                        NA, NA)))
  expect_equal(r$values$mouse_ID, c("26-1-1", "26-1-2", "26-1-1", "26-1-1", "26-1-2"))
  expect_equal(r$values$`anti-MPO-titer`,
               c("2.6658", "2.6729", "2.6712", "2.5228", "OOR<"))
  expect_setequal(r$dropped, c("ID", "group"))
  expect_null(r$range)
})

test_that("analyte and tissue columns (any case) give one column per analyte", {
  raw <- tibble::tibble(
    `Lab ID` = c("m1", "m1", "m2", "m2"), Analyte = c("TNFa", "IFNg", "TNFa", "IFNg"),
    TISSUE = "serum", `46147` = c("1", "2", "3", "4"), final = c("5", "6", "7", "8")
  )
  r <- elisa_reshape_wide_table_(raw)
  expect_named(r$values, c("mouse_ID", "tissue", "date", "TNFa", "IFNg"))
  expect_equal(nrow(r$values), 4)
  expect_equal(unique(r$values$tissue), "serum")
  expect_equal(r$values$TNFa, c(1, 3, 5, 7))
  expect_equal(r$values$IFNg, c(2, 4, 6, 8))
})

test_that("tissue argument fills empty cells and a missing column; neither = no column", {
  raw <- tibble::tibble(`Lab ID` = c("m1", "m2"), tissue = c("serum", NA),
                        `46147` = c("1", "2"))
  r <- elisa_reshape_wide_table_(raw, analyte = "a", tissue = "plasma")
  expect_equal(r$values$tissue, c("serum", "plasma"))
  raw2 <- dplyr::select(raw, !tissue)
  expect_equal(elisa_reshape_wide_table_(raw2, analyte = "a", tissue = "plasma")$values$tissue,
               c("plasma", "plasma"))
  expect_false("tissue" %in% names(elisa_reshape_wide_table_(raw2, analyte = "a")$values))
})

test_that("analyte argument fills empty analyte cells; neither errors", {
  raw <- tibble::tibble(`Lab ID` = c("m1", "m2"), analyte = c("TNFa", NA),
                        `46147` = c("1", "2"))
  r <- elisa_reshape_wide_table_(raw, analyte = "TNFa")
  expect_named(r$values, c("mouse_ID", "date", "TNFa"))
  expect_equal(r$values$TNFa, c(1, 2))
  expect_error(elisa_reshape_wide_table_(raw), "analyte")
  expect_error(elisa_reshape_wide_table_(dplyr::select(raw, !analyte)), "analyte")
})

test_that("repeated mouse x analyte rows are replicates in row order", {
  raw <- tibble::tibble(`Lab ID` = c("m1", "m1", "m2"), `46147` = c("1", "2", "3"))
  r <- elisa_reshape_wide_table_(raw, analyte = "TNFa")
  expect_equal(r$values$mouse_ID, c("m1", "m1", "m2"))
  expect_equal(r$values$TNFa, c(1, 2, 3))
})

test_that("all-numeric analyte columns become numeric, others stay text", {
  raw <- tibble::tibble(`Lab ID` = c("m1", "m1", "m2", "m2"),
                        analyte = c("A", "B", "A", "B"), `46147` = c("1", "OOR<", "2", "3"))
  r <- elisa_reshape_wide_table_(raw)
  expect_type(r$values$A, "double")
  expect_type(r$values$B, "character")
})

test_that("unit / lloq / uloq build a range table (scalar or named by analyte)", {
  raw <- tibble::tibble(`Lab ID` = c("m1", "m1"), analyte = c("TNFa", "IFNg"),
                        `46147` = c("1", "2"))
  r <- elisa_reshape_wide_table_(raw, lloq = c(TNFa = 7.8, IFNg = 5), unit = "pg/ml")
  expect_equal(r$range[[1]], c("LLOQ", "ULOQ", "unit"))
  expect_equal(r$range$TNFa, c("7.8", NA, "pg/ml"))
  expect_equal(r$range$IFNg, c("5", NA, "pg/ml"))
  expect_error(elisa_reshape_wide_table_(raw, lloq = c(IL6 = 1)), "IL6")
})

test_that("errors: missing id column, no timepoints, duplicate timepoints", {
  expect_error(elisa_reshape_wide_table_(wide_raw_(), analyte = "a", id_col = "mouse"),
               "mouse")
  raw <- tibble::tibble(`Lab ID` = "m1", group = "x")
  expect_error(elisa_reshape_wide_table_(raw, analyte = "a"), "timepoint")
  dup <- tibble::tibble(`Lab ID` = "m1", `46147` = "1", `2026-05-05` = "2")
  expect_error(elisa_reshape_wide_table_(dup, analyte = "a"), "46147.*2026-05-05")
  dup2 <- tibble::tibble(`Lab ID` = "m1", final = "1", End = "2")
  expect_error(elisa_reshape_wide_table_(dup2, analyte = "a"), "final.*End")
})

test_that("default id_col accepts mouse_ID or Lab ID (first found)", {
  raw <- tibble::tibble(mouse_ID = c("m1", "m2"), `46147` = c("1", "2"))
  r <- elisa_reshape_wide_table_(raw, analyte = "a")
  expect_equal(r$values$mouse_ID, c("m1", "m2"))
  expect_error(elisa_reshape_wide_table_(dplyr::rename(raw, id = mouse_ID), analyte = "a"),
               "mouse_ID, Lab ID")
})

test_that("rows without an ID are dropped (silently when empty)", {
  raw <- tibble::tibble(`Lab ID` = c("m1", NA, " "), `46147` = c("1", "2", "3"))
  expect_warning(r <- elisa_reshape_wide_table_(raw, analyte = "a"), "2 row")
  expect_equal(r$values$mouse_ID, "m1")
  pad <- tibble::tibble(`Lab ID` = c("m1", NA), `46147` = c("1", NA))
  expect_no_warning(elisa_reshape_wide_table_(pad, analyte = "a"))
})

test_that("empty cells give no row (no fake samples)", {
  r <- elisa_reshape_wide_table_(wide_raw_(), analyte = "a")
  expect_equal(nrow(r$values), 5)
  expect_false(any(r$values$mouse_ID == "26-1-2" & r$values$date %in% as.Date("2026-05-19")))
  raw <- tibble::tibble(`Lab ID` = c("m1", "m1"), analyte = c("A", "B"),
                        `46147` = c("1", NA), final = c(NA, NA))
  r2 <- elisa_reshape_wide_table_(raw)
  expect_equal(nrow(r2$values), 1)
  expect_equal(r2$values$A, 1)
})

test_that("date-time serial and dd.mm.yyyy headers are timepoints", {
  expect_equal(parse_wide_date_header_(c("46161.4375", "19.05.2026")),
               as.Date(c("2026-05-19", "2026-05-19")))
})

test_that("rows without an ID but with values warn", {
  raw <- tibble::tibble(`Lab ID` = c("m1", NA, NA), `46147` = c("1", "2", NA))
  expect_warning(r <- elisa_reshape_wide_table_(raw, analyte = "a"), "1 row")
  expect_equal(r$values$mouse_ID, "m1")
})

test_that("parse_wide_date_header_ reads serials and ISO text only", {
  expect_equal(parse_wide_date_header_(c("46147", "2026-05-19", "final", "ID", NA)),
               as.Date(c("2026-05-05", "2026-05-19", NA, NA, NA)))
})

wide_file_ <- function(name = "titer_wide", env = parent.frame()) {
  testthat::skip_if_not_installed("writexl")
  f <- file.path(withr::local_tempdir(.local_envir = env), paste0(name, ".xlsx"))
  writexl::write_xlsx(list(Tabelle1 = wide_raw_()), f)
  f
}

test_that("round trip through elisa_read_curated", {
  f <- wide_file_()
  p <- suppressMessages(elisa_reshape_wide(f, analyte = "anti-MPO-titer"))
  expect_equal(basename(p), "titer_long.xlsx")
  d <- suppressMessages(elisa_read_curated(p, tissue = "serum"))
  expect_equal(nrow(d), 5)
  expect_equal(sum(is.na(d$collected_at)), 2)
  expect_true("OOR<" %in% d$result_status)
  expect_equal(unique(d$analyte), "anti-MPO-titer")
  expect_equal(sort(unique(d$collected_at)), as.Date(c("2026-05-05", "2026-05-19")))
})

test_that("out without _wide gets _long appended; explicit out honoured", {
  f <- wide_file_("titer")
  p <- suppressMessages(elisa_reshape_wide(f, analyte = "a"))
  expect_equal(basename(p), "titer_long.xlsx")
  o <- file.path(dirname(f), "custom.xlsx")
  expect_equal(suppressMessages(elisa_reshape_wide(f, analyte = "a", out = o)), o)
  expect_true(file.exists(o))
})

test_that("unchanged content is not rewritten", {
  f <- wide_file_()
  p <- suppressMessages(elisa_reshape_wide(f, analyte = "a"))
  Sys.setFileTime(p, Sys.time() - 60)
  before <- file.mtime(p)
  expect_message(elisa_reshape_wide(f, analyte = "a"), "up to date")
  expect_equal(file.mtime(p), before)
})

test_that("changed arguments rewrite the file", {
  f <- wide_file_()
  p <- suppressMessages(elisa_reshape_wide(f, analyte = "a"))
  expect_false("range" %in% readxl::excel_sheets(p))
  suppressMessages(elisa_reshape_wide(f, analyte = "a", lloq = 0.1))
  expect_true("range" %in% readxl::excel_sheets(p))
  d <- suppressMessages(elisa_read_curated(p, tissue = "serum"))
  expect_equal(unique(d$lloq), 0.1)
})

test_that("dropped columns are messaged", {
  f <- wide_file_()
  expect_message(elisa_reshape_wide(f, analyte = "a"), "ID, group")
})
