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

# one irrel mouse ends early at week 4 (the 26-1-17 case)
tc_low <- dplyr::bind_rows(
  tc,
  tibble::tibble(group = factor("irrel", levels = c("MPO409", "irrel")),
                 week_num = 4L, mouse_id = "irrel-9", value = 7)
)
tc_layer <- function(p, geom) {
  p$layers[[which(vapply(p$layers, function(l) inherits(l$geom, geom), logical(1)))[1]]]
}

test_that("summary computes median and IQR over mice, not wells", {
  tc_reps <- dplyr::bind_rows(tc, dplyr::mutate(tc, value = value + 2))
  s <- timecourse_summary_(timecourse_prep_(tc_reps, week_num, value, group, mouse_id, NULL),
                           min_n = 2)
  one <- s[s$.fill == "irrel" & s$.x == 1, ]
  expect_equal(one$n, 4L)
  vals <- tc$value[tc$group == "irrel" & tc$week_num == 1] + 1
  expect_equal(one$median, stats::median(vals))
  expect_equal(one$q25, unname(stats::quantile(vals, 0.25, type = 7)))
  expect_equal(one$q75, unname(stats::quantile(vals, 0.75, type = 7)))
  expect_false(any(s$low_n))
})

test_that("low-n points are hollow, without error bar, and bridged by the line", {
  p <- plot_timecourse_summary(tc_low, colors = tc_cols)
  expect_s3_class(p, "ggplot")
  bars <- tc_layer(p, "GeomErrorbar")$data
  expect_false(any(bars$.x == 4 & bars$.fill == "irrel"))
  line <- tc_layer(p, "GeomLine")$data
  irrel_x <- line$.x[line$.fill == "irrel"]
  expect_false(4 %in% irrel_x)
  expect_true(all(c(3, 5) %in% irrel_x))
  hollow <- Filter(function(l) inherits(l$geom, "GeomPoint") && identical(l$aes_params$shape, 21),
                   p$layers)
  expect_length(hollow, 1)
  expect_equal(hollow[[1]]$data$.x, 4)
})

test_that("n labels are printed per group and week", {
  p <- plot_timecourse_summary(tc_low)
  txt <- tc_layer(p, "GeomText")$data
  expect_equal(nrow(txt), 9)
  expect_equal(txt$n[txt$.x == 4], 1L)
})

test_that("default caption names the measure", {
  expect_equal(plot_timecourse_summary(tc)$labels$caption, "median and IQR; n per group below")
})

test_that("summary p_val adds stars without dropping the n labels", {
  p <- plot_timecourse_summary(tc, p_val = TRUE)
  n_text <- sum(vapply(p$layers, function(l) inherits(l$geom, "GeomText"), logical(1)))
  expect_equal(n_text, 2)
})

test_that("full plot median line bridges time points with fewer than min_n mice", {
  p <- plot_timecourse(tc_low)
  thick <- Filter(function(l) inherits(l$geom, "GeomLine") && identical(l$aes_params$linewidth, 1.2),
                  p$layers)[[1]]$data
  expect_false(4 %in% thick$.x[thick$.fill == "irrel"])
  expect_true(all(c(3, 5) %in% thick$.x[thick$.fill == "irrel"]))
})

test_that("hollow low-n points stay out of the legend", {
  p <- plot_timecourse_summary(tc_low)
  hollow <- Filter(function(l) identical(l$aes_params$shape, 21), p$layers)[[1]]
  expect_false(hollow$show.legend)
})

test_that("summary keeps groups side by side when one group is below min_n", {
  # irrel has one mouse at week 3, MPO409 four
  d <- dplyr::filter(tc, !(group == "irrel" & week_num == 3 & mouse_id != "irrel-1"))
  b <- ggplot2::ggplot_build(plot_timecourse_summary(d))
  pts <- dplyr::bind_rows(lapply(seq_along(b$data), function(i) {
    if (inherits(b$plot$layers[[i]]$geom, "GeomPoint")) b$data[[i]][, c("x", "y", "colour")]
  }))
  at3 <- pts[abs(pts$x - 3) < 0.5, ]
  expect_equal(length(unique(round(at3$x, 6))), 2)
  txt <- b$data[[which(vapply(b$plot$layers, function(l) inherits(l$geom, "GeomText"), logical(1)))]]
  expect_setequal(round(txt$x[abs(txt$x - 3) < 0.5], 6), round(at3$x, 6))
})

test_that("several analytes per mouse and time point error instead of averaging", {
  tc2 <- dplyr::bind_rows(dplyr::mutate(tc, analyte = "IFNg"), dplyr::mutate(tc, analyte = "IL-17"))
  expect_error(plot_timecourse(tc2), "several analyte")
  expect_error(plot_timecourse_summary(tc2), "several analyte")
  expect_s3_class(plot_timecourse(tc2, facet = analyte), "ggplot")
})
