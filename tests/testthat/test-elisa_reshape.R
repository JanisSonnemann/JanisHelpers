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
                                        "2026-05-19", NA, NA)))
  expect_equal(r$values$mouse_ID, rep(c("26-1-1", "26-1-2"), 3))
  expect_equal(r$values$`anti-MPO-titer`,
               c("2.6658", "2.6729", "2.6712", NA, "2.5228", "OOR<"))
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

test_that("rows without an ID are dropped", {
  raw <- tibble::tibble(`Lab ID` = c("m1", NA, " "), `46147` = c("1", "2", "3"))
  r <- elisa_reshape_wide_table_(raw, analyte = "a")
  expect_equal(r$values$mouse_ID, "m1")
})

test_that("parse_wide_date_header_ reads serials and ISO text only", {
  expect_equal(parse_wide_date_header_(c("46147", "2026-05-19", "final", "ID", NA)),
               as.Date(c("2026-05-05", "2026-05-19", NA, NA, NA)))
})
