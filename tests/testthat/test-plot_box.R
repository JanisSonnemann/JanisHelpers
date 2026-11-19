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
