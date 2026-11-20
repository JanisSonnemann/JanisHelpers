pop_data <- tibble::tibble(
  mouse_id = rep(paste0("m", 1:4), each = 2),
  group = rep(c("irrel", "MPO409"), each = 4),
  tissue = rep(c("kidney", "lung"), 4),
  population = "CD4+", metric = "count_per_g",
  value = c(0, 10, 5, 20, 50, 60, 70, 80)
) |> dplyr::bind_rows(tibble::tibble(
  mouse_id = "m1", group = "irrel", tissue = "kidney",
  population = "CD8+", metric = "count_per_g", value = 3
))

has_scale_ <- function(p, aes) {
  any(vapply(p$scales$scales, function(s) aes %in% s$aesthetics, logical(1)))
}

test_that("plot_box returns a ggplot with the requested scales", {
  p <- plot_box(pop_data, x = tissue, y = value, fill = group, scale = "log",
                colors = c(irrel = "#F8766D", MPO409 = "#00BFC4"))
  expect_s3_class(p, "ggplot")
  expect_true(has_scale_(p, "fill"))
  expect_true(has_scale_(p, "y"))
  expect_false(has_scale_(plot_box(pop_data, tissue, value, group), "fill"))
  expect_error(plot_box(pop_data, tissue, value, group, scale = "sqrt"))
})

test_that("plot_population filters and offsets values for log scale", {
  p <- plot_population(pop_data, populations = "CD4+", metric = "count_per_g",
                       tissues = "kidney", x = group, scale = "log")
  expect_equal(nrow(p$data), 4)
  expect_equal(sort(p$data$value), c(1, 6, 51, 71))
})

test_that("plot_population adds a p-value layer on request", {
  p <- plot_population(pop_data, "CD4+", "count_per_g", x = group, p_val = TRUE)
  expect_true(any(vapply(p$layers, function(l) inherits(l$stat, "StatCompareMeans"), logical(1))))
  p0 <- plot_population(pop_data, "CD4+", "count_per_g", x = group)
  expect_false(any(vapply(p0$layers, function(l) inherits(l$stat, "StatCompareMeans"), logical(1))))
})

paired_data <- tibble::tibble(
  mouse_id = rep(c("m1", "m2"), each = 3),
  group = rep(c("irrel", "MPO409"), each = 3),
  condition = rep(c("unstim", "MPO", "PMA-Iono"), 2),
  population = "CD4+_CD44+_IFNg+", value = c(0.01, 0.03, 0.2, 0.02, 0.05, 0.3)
)

test_that("plot_paired keeps only the requested conditions, in order", {
  p <- plot_paired(paired_data, levels = c("unstim", "MPO"))
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$condition), c("unstim", "MPO"))
  expect_false("PMA-Iono" %in% p$data$condition)
  expect_equal(nrow(p$data), 4)
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomLine"), logical(1))))
})

test_that("every plot helper renders without ggpubr attached", {
  expect_false("package:ggpubr" %in% search())
  expect_no_error(ggplot2::ggplot_build(
    plot_box(pop_data, tissue, value, group, scale = "log")
  ))
  p <- plot_population(pop_data, "CD4+", "count_per_g", x = group, p_val = TRUE)
  built <- expect_no_error(ggplot2::ggplot_build(p))
  p_layer <- which(vapply(p$layers, function(l) inherits(l$stat, "StatCompareMeans"), logical(1)))
  expect_match(built$data[[p_layer]]$label, "^p = ")
  expect_no_error(ggplot2::ggplot_build(
    plot_paired(paired_data, levels = c("unstim", "MPO"))
  ))
})

test_that("plot_paired errors on a level not in the data", {
  expect_error(plot_paired(paired_data, levels = c("unstim", "OVA")), "OVA")
})

paired_colour_data_ <- tibble::tibble(
  mouse_id = rep(c("m1", "m2"), each = 2), group = rep(c("a", "b"), each = 2),
  condition = rep(c("unstim", "MPO"), 2), population = "CD4+", value = c(0.1, 0.2, 0.3, 0.4)
)
colour_data_ <- dplyr::mutate(pop_data, group = ifelse(group == "irrel", "a", "b"))
box_fills_ <- function(p) unique(ggplot2::ggplot_build(p)$data[[1]]$fill)
three_plots_ <- function(...) list(
  plot_box(colour_data_, x = tissue, y = value, fill = group, ...),
  plot_population(colour_data_, populations = "CD4+", metric = "count_per_g", ...),
  plot_paired(paired_colour_data_, levels = c("unstim", "MPO"), ...)
)

test_that("plot helpers use the group colour option when colors is missing", {
  withr::local_options(JanisHelpers.group_colors = c(a = "#111111", b = "#222222"))
  for (p in three_plots_()) expect_setequal(box_fills_(p), c("#111111", "#222222"))
})

test_that("explicit colors beat the option", {
  withr::local_options(JanisHelpers.group_colors = c(a = "#111111", b = "#222222"))
  for (p in three_plots_(colors = c(a = "#333333", b = "#444444"))) {
    expect_setequal(box_fills_(p), c("#333333", "#444444"))
  }
})

test_that("group colour option is skipped when fill is not the group", {
  withr::local_options(JanisHelpers.group_colors = c(a = "#111111", b = "#222222"))
  by_tissue <- box_fills_(plot_box(colour_data_, x = group, y = value, fill = tissue))
  by_cond <- box_fills_(plot_paired(paired_colour_data_, levels = c("unstim", "MPO"),
                                    fill = condition))
  for (f in list(by_tissue, by_cond)) {
    expect_length(f, 2)
    expect_false(any(f %in% c("grey50", "#7F7F7F")))
  }
})

test_that("explicit NULL gives ggplot default fills", {
  withr::local_options(JanisHelpers.group_colors = c(a = "#111111", b = "#222222"))
  for (p in three_plots_(colors = NULL)) {
    expect_length(intersect(box_fills_(p), c("#111111", "#222222")), 0)
  }
})

# -- excluded mice and labels -------------------------------------------------

excl_data_ <- dplyr::mutate(
  dplyr::filter(pop_data, population == "CD4+"),
  excluded = c(FALSE, FALSE, TRUE, NA, FALSE, FALSE, FALSE, TRUE)
)
# GeomText matches ggrepel's GeomTextRepel too (labels with or without ggrepel)
layer_idx_ <- function(p, geom) {
  if (geom == "GeomText") geom <- c("GeomText", "GeomTextRepel")
  which(vapply(p$layers, function(l) inherits(l$geom, geom), logical(1)))
}

test_that("plot_box drops excluded rows with a message by default", {
  expect_message(p <- plot_box(excl_data_, tissue, value, group), "Dropped 2 excluded")
  expect_equal(nrow(p$data), 6)
  expect_length(layer_idx_(p, "GeomText"), 0)
})

test_that("NA in excluded counts as kept", {
  p <- suppressMessages(plot_box(excl_data_, tissue, value, group))
  expect_true(any(is.na(p$data$excluded)))
})

test_that("plot_box draws excluded rows as red crosses with labels, outside the boxes", {
  p <- plot_box(excl_data_, tissue, value, group, show_excluded = TRUE)
  expect_equal(nrow(p$data), 6)                       # boxes: kept rows only
  pt_i <- layer_idx_(p, "GeomPoint")
  lab_i <- layer_idx_(p, "GeomText")
  expect_equal(nrow(p$layers[[pt_i]]$data), 8)
  b <- ggplot2::ggplot_build(p)
  pts <- b$data[[pt_i]]
  expect_equal(sum(pts$colour == "red3"), 2)
  expect_true(all(pts$shape[pts$colour == "red3"] == 4))
  labs <- b$data[[lab_i]]
  expect_setequal(labs$label[labs$label != ""], c("m2", "m4"))
  expect_true(all(labs$colour[labs$label != ""] == "red3"))
  expect_equal(labs$x, pts$x)                         # labels on their points
})

test_that("show_excluded without an excluded column adds nothing", {
  expect_message(p <- plot_box(pop_data, tissue, value, group, show_excluded = TRUE),
                 "no `excluded` column")
  expect_length(layer_idx_(p, "GeomText"), 0)
  expect_false(has_scale_(p, "shape"))
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot_box labels every point on request", {
  p <- plot_box(pop_data, tissue, value, group, label = TRUE)
  b <- ggplot2::ggplot_build(p)
  expect_setequal(b$data[[layer_idx_(p, "GeomText")]]$label, pop_data$mouse_id)
})

test_that("plot_box needs no id column when nothing is labelled", {
  d <- dplyr::select(pop_data, !mouse_id)
  expect_no_error(ggplot2::ggplot_build(plot_box(d, tissue, value, group)))
})

test_that("plot_population passes show_excluded and label through", {
  p <- plot_population(excl_data_, "CD4+", "count_per_g", x = group, scale = "log",
                       show_excluded = TRUE, label = TRUE)
  pts <- p$layers[[layer_idx_(p, "GeomPoint")]]$data
  expect_equal(sort(pts$value), sort(excl_data_$value + 1))
  expect_length(layer_idx_(p, "GeomText"), 1)
})

test_that("plot_paired x labels are the condition names, one cluster per group", {
  p <- plot_paired(paired_data, levels = c("unstim", "MPO"))
  b <- ggplot2::ggplot_build(p)
  expect_equal(b$layout$panel_params[[1]]$x$get_labels(), rep(c("unstim", "MPO"), 2))
  expect_equal(sort(unique(p$data$.x)), c(1, 2, 3.8, 4.8))
  grp <- p$layers[[layer_idx_(p, "GeomText")]]$data
  expect_setequal(grp$.grp, c("irrel", "MPO409"))
})

test_that("plot_paired shades conditions and shows excluded wells as red crosses", {
  d <- dplyr::mutate(paired_data, excluded = mouse_id == "m2" & condition == "MPO")
  p <- plot_paired(d, levels = c("unstim", "MPO"), show_excluded = TRUE, label = TRUE)
  expect_true(has_scale_(p, "alpha"))
  expect_equal(nrow(p$data), 3)
  b <- ggplot2::ggplot_build(p)
  expect_equal(sum(b$data[[layer_idx_(p, "GeomPoint")]]$colour == "red3"), 1)
})
