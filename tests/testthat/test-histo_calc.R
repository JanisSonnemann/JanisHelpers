test_that("histo_calc_ratio adds a percentage metric", {
  d <- suppressMessages(histo_read_curated(test_path("../fixtures/histo_curated.xlsx"), tissue = "kidney")) |>
    histo_calc_ratio("crescents", "gloms", "pct_crescents")
  expect_equal(d$value[d$mouse_ID == "M1" & d$metric == "pct_crescents"], 22 / 85 * 100)
  expect_true(is.na(d$value[d$mouse_ID == "M3" & d$metric == "pct_crescents"]))
})

test_that("histo_calc_ratio guards names and zero denominators", {
  d <- suppressMessages(histo_read_curated(test_path("../fixtures/histo_curated.xlsx"), tissue = "kidney"))
  expect_error(histo_calc_ratio(d, "crescents", "gloms", "gloms"), "already exists")
  expect_error(histo_calc_ratio(d, "nope", "gloms", "x"), "not in")
  d0 <- d
  d0$value[d0$mouse_ID == "M1" & d0$metric == "gloms"] <- 0
  r <- histo_calc_ratio(d0, "crescents", "gloms", "p")
  expect_true(is.na(r$value[r$mouse_ID == "M1" & r$metric == "p"]))
})

test_that("histo_calc_ratio computes several named ratios", {
  d <- suppressMessages(histo_read_curated(test_path("../fixtures/histo_curated.xlsx"), tissue = "kidney"))
  r <- histo_calc_ratio(d, c(pct_cresc = "crescents", pct_nec = "necrosis"), "gloms")
  expect_equal(r$value[r$mouse_ID == "M1" & r$metric == "pct_cresc"], 22 / 85 * 100)
  expect_true(all(c("pct_cresc", "pct_nec") %in% r$metric))
  expect_equal(nrow(r), nrow(d) + 2 * dplyr::n_distinct(d$mouse_ID))
})

test_that("histo_calc_ratio errors on bad vector input", {
  d <- suppressMessages(histo_read_curated(test_path("../fixtures/histo_curated.xlsx"), tissue = "kidney"))
  expect_error(histo_calc_ratio(d, c("crescents", "necrosis"), "gloms"), "named")
  expect_error(histo_calc_ratio(d, c(a = "crescents", a = "necrosis"), "gloms"), "unique")
  expect_error(histo_calc_ratio(d, c(a = "crescents", b = "necrosis"),
                                c("gloms", "gloms", "gloms")), "denominator")
  expect_error(histo_calc_ratio(d, c(gloms = "crescents"), "gloms"), "already exists")
})

test_that("histo_calc_ratio keeps extra columns", {
  d <- suppressMessages(histo_read_curated(test_path("../fixtures/histo_curated.xlsx"), tissue = "kidney"))
  d$comment <- ifelse(d$mouse_ID == "M1", "torn", NA_character_)
  r <- histo_calc_ratio(d, c(p = "crescents"), "gloms")
  expect_equal(r$comment[r$mouse_ID == "M1" & r$metric == "p"], "torn")
})
