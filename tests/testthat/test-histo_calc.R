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
