fixture_data <- function() {
  suppressMessages(elisa_read_curated(
    test_path("..", "fixtures", "elisa_curated.xlsx"), tissue = "spleen"
  ))
}

test_that("qc_timepoint_ formats dates and NA as end", {
  expect_equal(qc_timepoint_(as.Date(c("2026-05-06", NA))), c("2026-05-06", "end"))
})

test_that("elisa_qc_range_ counts out-of-range values per analyte", {
  d <- fixture_data()
  r <- elisa_qc_range_(d)
  tnfa <- d[d$analyte == "TNFa", ]
  expect_equal(r$n_below_lloq[r$analyte == "TNFa"],
               sum(tnfa$value < tnfa$lloq, na.rm = TRUE))
  expect_equal(r$n_oor_low[r$analyte == "IL-2"], 1L)
  expect_named(r, c("analyte", "n", "n_oor_low", "n_oor_high",
                    "n_below_lloq", "n_above_uloq"))
})

test_that("elisa_qc_replicates_ lists samples above the CV threshold", {
  d <- fixture_data()
  hit <- elisa_qc_replicates_(d, max_cv = 5)
  expect_equal(hit$mouse_ID, "M1")
  expect_equal(hit$analyte, "TNFa")
  expect_equal(round(hit$cv, 1), 5.2)
  expect_match(hit$values, "4.2")
  expect_equal(nrow(elisa_qc_replicates_(d, max_cv = 20)), 0L)
})

test_that("elisa_qc_means_ averages replicates per sample x analyte", {
  m <- elisa_qc_means_(fixture_data())
  v <- m$value[m$mouse_ID == "M1" & m$analyte == "TNFa" & is.na(m$collected_at)]
  expect_equal(v, mean(c(4.2, 3.9)))
  expect_false("replicate" %in% names(m))
})

test_that("qc_table_completeness_ lists missing mice per tissue x timepoint", {
  d <- fixture_data()
  comp <- qc_table_completeness_(d, c("M1", "M2", "M3"))
  expect_equal(comp$missing[comp$timepoint == "end"], "M3")
  expect_equal(comp$missing[comp$timepoint == "2026-05-06"], "M2, M3")
  expect_length(attr(comp, "unknown"), 0L)
  unk <- qc_table_completeness_(d, "M1")
  expect_equal(attr(unk, "unknown"), "M2")
})

test_that("qc_table_empty_ counts NA values", {
  d <- fixture_data()
  e <- qc_table_empty_(d, "analyte")
  expect_equal(sum(e$n_empty), sum(is.na(d$value)))
  expect_equal(sum(e$n_values), nrow(d))
})

synthetic_scores <- function() {
  tibble::tibble(
    mouse_ID = paste0("S", 1:6), tissue = "spleen",
    collected_at = as.Date(NA), analyte = "A", grp = "g",
    value = c(1, 2, 2, 3, 2, 50)
  )
}

test_that("qc_table_score_ flags the extreme point only", {
  s <- qc_table_score_(synthetic_scores(), "analyte", "grp", FALSE)
  expect_equal(s$mouse_ID[s$flag], "S6")
  expect_equal(s$side[s$flag], "high")
  none <- qc_table_score_(synthetic_scores(), "analyte", NULL, FALSE)
  expect_false(any(none$flag))
})

test_that("qc_table_score_ scores on the log axis when asked", {
  s <- qc_table_score_(synthetic_scores(), "analyte", "grp", TRUE)
  expect_equal(s$t, log10(synthetic_scores()$value))
})

test_that("qc_table_plot_ and flag table work", {
  s <- qc_table_score_(synthetic_scores(), "analyte", "grp", FALSE)
  s$batch <- rep(c("b1", "b2"), 3)
  expect_s3_class(qc_table_plot_(s, "analyte", "grp", FALSE), "ggplot")
  expect_s3_class(qc_table_plot_(s, "analyte", "grp", TRUE, batch_col = "batch"),
                  "ggplot")
  ft <- qc_table_flag_table_(s, "analyte")
  expect_equal(ft$mouse_ID, "S6")
  expect_equal(ft$n_flagged, 1L)
})

test_that("qc_table_time_plot_ needs two dates", {
  d <- synthetic_scores()
  expect_null(qc_table_time_plot_(d, "analyte", "grp"))
  d$collected_at <- as.Date("2026-05-01") + c(0, 0, 0, 7, 7, 7)
  expect_s3_class(qc_table_time_plot_(d, "analyte", "grp"), "ggplot")
})
