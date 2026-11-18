samples <- tibble::tibble(
  mouse_ID = c("m1", "m1", "m2", "m3"),
  tissue = c("spleen", "spleen", "spleen", "spleen"),
  collected_at = as.Date(c("2026-05-06", "2026-06-10", "2026-06-10", NA))
)

test_that("a mouse_ID-only row excludes all that mouse's samples", {
  ex <- tibble::tibble(mouse_ID = "m1", reason = "bad")
  res <- match_sample_exclusions_(samples, ex, "titer")
  expect_equal(nrow(res$excluded_samples), 2)
  expect_equal(res$applicable$n_samples_matched, 2L)
  expect_s3_class(res$applicable$collected_at, "Date")
  expect_equal(nrow(res$unmatched), 0)
})

test_that("a collected_at row excludes only that date", {
  ex <- tibble::tibble(mouse_ID = "m1", collected_at = "2026-05-06", reason = "bad")
  res <- match_sample_exclusions_(samples, ex, "titer")
  expect_equal(res$excluded_samples$collected_at, as.Date("2026-05-06"))
})

test_that("an exclusion date never matches an NA sample date", {
  ex <- tibble::tibble(mouse_ID = "m3", collected_at = "2026-06-10", reason = "bad")
  res <- match_sample_exclusions_(samples, ex, "titer")
  expect_equal(nrow(res$unmatched), 1)
})

test_that("rows for other panels are ignored; matching nothing is unmatched", {
  ex <- tibble::tibble(
    mouse_ID = c("m1", "zz"), panel = c("overview", NA), reason = "bad"
  )
  res <- match_sample_exclusions_(samples, ex, "titer")
  expect_equal(nrow(res$applicable), 1)
  expect_equal(res$unmatched$mouse_ID, "zz")
})

test_that("file_name / condition rows error, missing fields error", {
  ex <- tibble::tibble(mouse_ID = "m1", file_name = "a.fcs", reason = "bad")
  expect_error(match_sample_exclusions_(samples, ex, "titer"), "don't apply")
  ex2 <- tibble::tibble(mouse_ID = "m1", panel = "overview", file_name = "a", reason = "x")
  expect_no_error(match_sample_exclusions_(samples, ex2, "titer"))
  expect_error(
    match_sample_exclusions_(samples, tibble::tibble(mouse_ID = "m1"), "titer"),
    "needs column"
  )
  expect_error(
    match_sample_exclusions_(samples, tibble::tibble(mouse_ID = "m1", reason = ""), "titer"),
    "needs a mouse_ID and a reason"
  )
})

test_that("NULL exclusions give empty results", {
  res <- match_sample_exclusions_(samples, NULL, "titer")
  expect_equal(nrow(res$applicable), 0)
  expect_equal(nrow(res$excluded_samples), 0)
})
