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

# ---- end dates / duplicates / exclusion keys ---------------------------------

test_that("qc_end_dates_ normalises named vectors and data frames", {
  v <- qc_end_dates_(c(M1 = "2026-06-10", M2 = NA))
  expect_equal(v, c(M1 = as.Date("2026-06-10")))
  df <- tibble::tibble(mouse_ID = c("M1", "M2", NA),
                       death_date = as.POSIXct(c("2026-06-10", "2026-06-11", "2026-01-01"),
                                               tz = "UTC"))
  expect_equal(qc_end_dates_(df),
               c(M1 = as.Date("2026-06-10"), M2 = as.Date("2026-06-11")))
  expect_null(qc_end_dates_(NULL))
  expect_error(qc_end_dates_(as.Date("2026-06-10")), "named")
  expect_error(qc_end_dates_(tibble::tibble(mouse_ID = "M1", a = 1, b = 2)),
               "one date column")
  expect_error(qc_end_dates_(c(M1 = "2026-06-10", M1 = "2026-06-11")),
               "several dates.*M1")
})

test_that("resolved end dates put the end bleed on the time plot", {
  d <- fixture_data()
  ends <- c(M1 = as.Date("2026-06-10"), M2 = as.Date("2026-06-11"))
  r <- qc_resolve_end_dates_(d, ends)
  expect_false(anyNA(r$collected_at))
  expect_equal(r$collected_at[d$mouse_ID == "M2"],
               rep(as.Date("2026-06-11"), sum(d$mouse_ID == "M2")))
  # unchanged without end dates (older preps); unknown mice stay NA
  expect_identical(qc_resolve_end_dates_(d, NULL)$collected_at, d$collected_at)
  expect_equal(sum(is.na(qc_resolve_end_dates_(d, ends["M1"])$collected_at)),
               sum(d$mouse_ID == "M2"))

  expect_null(qc_table_time_plot_(elisa_qc_means_(d), "analyte"))
  p <- qc_table_time_plot_(elisa_qc_means_(r), "analyte")
  expect_s3_class(p, "ggplot")
  m1 <- p$data[p$data$mouse_ID == "M1" & p$data$analyte == "TNFa", ]
  expect_setequal(m1$collected_at, as.Date(c("2026-05-06", "2026-06-10")))
  expect_true(as.Date("2026-06-11") %in% p$data$collected_at)
})

test_that("qc_duplicate_keys_ finds an empty-date row dated at the end date", {
  d <- fixture_data()
  key <- c("mouse_ID", "tissue", "collected_at", "analyte", "replicate")
  expect_equal(nrow(qc_duplicate_keys_(d, key)), 0)
  r <- qc_resolve_end_dates_(d, c(M1 = as.Date("2026-05-06")))
  dups <- qc_duplicate_keys_(r, key)
  expect_equal(nrow(dups), 2)   # M1 x 2026-05-06 x 2 analytes, replicate 1
  expect_equal(unique(dups$mouse_ID), "M1")
})

test_that("flag table and exclusion keys match on mouse x tissue x timepoint", {
  s <- qc_table_score_(synthetic_scores(), "analyte", "grp", FALSE)
  s$tissue <- "serum"
  s$collected_at <- as.Date(NA)
  s$timepoint <- "end"
  ex_samples <- tibble::tibble(mouse_ID = "S6", tissue = "serum",
                               collected_at = as.Date(NA), reason = "x")
  keys <- qc_excluded_keys_(ex_samples, c(S6 = as.Date("2026-06-10")))
  expect_equal(keys$timepoint, "2026-06-10")
  expect_equal(qc_table_flag_table_(s, "analyte", keys)$excluded, "")
  keys <- qc_excluded_keys_(ex_samples, NULL)
  expect_equal(keys$timepoint, "end")
  expect_equal(qc_table_flag_table_(s, "analyte", keys)$excluded, "excluded")
  # another tissue of the same mouse is not this sample
  keys$tissue <- "spleen"
  expect_equal(qc_table_flag_table_(s, "analyte", keys)$excluded, "")
  expect_false("excluded" %in% names(qc_table_flag_table_(s, "analyte")))
})
