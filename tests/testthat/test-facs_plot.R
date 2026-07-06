library(testthat)
library(JanisHelpers)

umap_input <- function() {
  tibble::tibble(
    file_name   = rep(c("a.fcs", "b.fcs"), each = 5),
    UMAP1       = rnorm(10),
    UMAP2       = rnorm(10),
    CD4         = rnorm(10),
    metacluster = factor(rep(1:2, 5)),
    group       = rep(c("control", "treated"), 5)
  )
}

test_that("facs_plot_umap() returns a ggplot object", {
  p <- facs_plot_umap(umap_input())
  expect_s3_class(p, "ggplot")
})

test_that("facs_plot_umap() errors when UMAP1/UMAP2 are missing", {
  dat <- umap_input() |> dplyr::select(!UMAP1)
  expect_error(facs_plot_umap(dat), "UMAP1")
})

test_that("facs_plot_umap() errors when color_by is not a column", {
  expect_error(facs_plot_umap(umap_input(), color_by = "NotAColumn"), "NotAColumn")
})

test_that("facs_plot_umap() errors when facet_by is not a column", {
  expect_error(facs_plot_umap(umap_input(), facet_by = "NotAColumn"), "NotAColumn")
})

test_that("facs_plot_umap() uses a continuous viridis scale for a double color_by", {
  p <- facs_plot_umap(umap_input(), color_by = "CD4")
  expect_s3_class(p$scales$get_scales("colour"), "ScaleContinuous")
})

test_that("facs_plot_umap() uses a discrete viridis scale for the default metacluster color_by", {
  p <- facs_plot_umap(umap_input())
  expect_s3_class(p$scales$get_scales("colour"), "ScaleDiscrete")
})

test_that("facs_plot_umap() facets by facet_by", {
  p <- facs_plot_umap(umap_input(), facet_by = "group")
  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$layout$layout), 2L)
})
