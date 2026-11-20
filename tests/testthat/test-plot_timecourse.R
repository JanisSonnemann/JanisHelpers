tc <- tidyr::expand_grid(
  group = factor(c("MPO409", "irrel"), levels = c("MPO409", "irrel")),
  m = 1:4,
  week_num = c(-1L, 1L, 3L, 5L)
) |>
  dplyr::mutate(
    mouse_id = paste0(group, "-", m),
    value = dplyr::if_else(group == "MPO409", 10, 1) * (week_num + 2) + m
  ) |>
  dplyr::select(!m)

tc_cols <- c(MPO409 = "#00BFC4", irrel = "#F8766D")

tc_has_scale <- function(p, aes) {
  any(vapply(p$scales$scales, function(s) aes %in% s$aesthetics, logical(1)))
}
tc_has_geom <- function(p, geom) {
  any(vapply(p$layers, function(l) inherits(l$geom, geom), logical(1)))
}

test_that("plot_timecourse returns a ggplot with group colours on lines", {
  p <- plot_timecourse(tc, colors = tc_cols)
  expect_s3_class(p, "ggplot")
  expect_true(tc_has_scale(p, "colour"))
  expect_true(tc_has_geom(p, "GeomLine"))
  expect_false(tc_has_scale(plot_timecourse(tc, colors = NULL), "colour"))
})

test_that("missing values are dropped with a message", {
  tc_na <- dplyr::mutate(tc, value = dplyr::if_else(dplyr::row_number() == 1, NA_real_, value))
  expect_message(plot_timecourse(tc_na), "Dropped 1 row")
})

test_that("replicates are averaged per mouse and week", {
  tc_reps <- dplyr::bind_rows(tc, dplyr::mutate(tc, value = value + 2))
  pr <- timecourse_prep_(tc_reps, week_num, value, group, mouse_id, NULL)
  expect_equal(nrow(pr), nrow(tc))
  expect_equal(sort(pr$.y), sort(tc$value + 1))
})

test_that("lloq line only when lloq present", {
  expect_true(tc_has_geom(plot_timecourse(dplyr::mutate(tc, lloq = 2)), "GeomHline"))
  expect_false(tc_has_geom(plot_timecourse(tc), "GeomHline"))
  expect_false(tc_has_geom(plot_timecourse(dplyr::mutate(tc, lloq = NA_real_)), "GeomHline"))
})

test_that("numeric x gets week labels at the weeks present", {
  b <- ggplot2::ggplot_build(plot_timecourse(tc))
  expect_equal(b$layout$panel_params[[1]]$x$get_labels(),
               c("week -1", "week 1", "week 3", "week 5"))
})

test_that("log scale with non-positive values errors", {
  expect_error(plot_timecourse(dplyr::mutate(tc, value = value - 100), scale = "log"), "positive")
  expect_s3_class(plot_timecourse(tc, scale = "log"), "ggplot")
})

test_that("p_val adds star labels", {
  expect_true(tc_has_geom(plot_timecourse(tc, p_val = TRUE), "GeomText"))
  expect_false(tc_has_geom(plot_timecourse(tc), "GeomText"))
})

test_that("facet splits by analyte", {
  tc2 <- dplyr::bind_rows(dplyr::mutate(tc, analyte = "a"), dplyr::mutate(tc, analyte = "b"))
  expect_s3_class(plot_timecourse(tc2, facet = analyte)$facet, "FacetWrap")
})

test_that("factor x is drawn as categories", {
  tcf <- dplyr::mutate(tc, week = factor(paste("week", week_num),
                                         levels = paste("week", c(-1, 1, 3, 5))))
  expect_s3_class(plot_timecourse(tcf, x = week), "ggplot")
  b <- ggplot2::ggplot_build(plot_timecourse(tcf, x = week))
  expect_equal(b$layout$panel_params[[1]]$x$get_labels(), levels(tcf$week))
})
