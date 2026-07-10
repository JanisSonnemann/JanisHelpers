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

abundance_input <- function() {
  tibble::tibble(
    file_name   = rep(c("a.fcs", "b.fcs", "c.fcs", "d.fcs"), each = 3),
    metacluster = factor(rep(1:3, 4)),
    fraction    = runif(12),
    group       = rep(c("control", "treated"), each = 6)
  )
}

test_that("facs_plot_cluster_abundance() returns a ggplot object", {
  p <- facs_plot_cluster_abundance(abundance_input(), group_col = "group")
  expect_s3_class(p, "ggplot")
})

test_that("facs_plot_cluster_abundance() errors when group_col is not a column", {
  expect_error(
    facs_plot_cluster_abundance(abundance_input(), group_col = "NotAColumn"),
    "NotAColumn"
  )
})

test_that("facs_plot_cluster_abundance() errors when cluster_col is not a column", {
  expect_error(
    facs_plot_cluster_abundance(abundance_input(), group_col = "group", cluster_col = "NotAColumn"),
    "NotAColumn"
  )
})

test_that("facs_plot_cluster_abundance() facets by cluster_col", {
  p <- facs_plot_cluster_abundance(abundance_input(), group_col = "group")
  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$layout$layout), 3L)
})

test_that("facs_plot_cluster_abundance() filters to significant clusters when significant_only = TRUE", {
  test_result <- tibble::tibble(
    metacluster = factor(1:3),
    contrast    = "treated_vs_control",
    p_val       = c(0.001, 0.2, 0.5),
    p_adj       = c(0.003, 0.3, 0.5)
  )

  p <- facs_plot_cluster_abundance(
    abundance_input(),
    test_result      = test_result,
    group_col        = "group",
    significant_only = TRUE,
    p_adj_threshold  = 0.05
  )
  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$layout$layout), 1L)
})

heatmap_input <- function() {
  tibble::tibble(
    metacluster = factor(c(1, 1, 2, 2)),
    marker      = c("CD44", "Foxp3", "CD44", "Foxp3"),
    median      = c(2, 5, 20, 1)
  )
}

test_that("facs_plot_cluster_heatmap() returns a ggplot object for scale = 'zscore'", {
  p <- facs_plot_cluster_heatmap(heatmap_input())
  expect_s3_class(p, "ggplot")
})

test_that("facs_plot_cluster_heatmap() returns a ggplot object for scale = 'raw'", {
  p <- facs_plot_cluster_heatmap(heatmap_input(), scale = "raw")
  expect_s3_class(p, "ggplot")
})

test_that("facs_plot_cluster_heatmap() z-scores per marker, reaching pure blue/red at the extremes", {
  # heatmap_input()'s per-marker values are symmetric around each marker's
  # own mean by construction (CD44: 2/20, mean 11; Foxp3: 5/1, mean 3), so
  # the z-scored fill data is exactly symmetric around 0 -- scale_fill_
  # gradient2()'s default midpoint = 0 rescaling then reaches pure "low"/
  # "high" colors at the min/max, confirmed directly (see the plan's
  # "Verified mechanics" section).
  p <- facs_plot_cluster_heatmap(heatmap_input(), scale = "zscore")
  built <- ggplot2::ggplot_build(p)
  expect_setequal(built$data[[1]]$fill, c("#0000FF", "#FF0000"))
})

test_that("facs_plot_cluster_heatmap() errors when cluster_col is not found", {
  expect_error(
    facs_plot_cluster_heatmap(heatmap_input(), cluster_col = "NotAColumn"),
    "NotAColumn"
  )
})

test_that("facs_plot_cluster_heatmap() errors when marker column is missing", {
  dat <- dplyr::select(heatmap_input(), !marker)
  expect_error(facs_plot_cluster_heatmap(dat), "marker")
})

test_that("facs_plot_cluster_heatmap() errors when median column is missing", {
  dat <- dplyr::select(heatmap_input(), !median)
  expect_error(facs_plot_cluster_heatmap(dat), "median")
})

test_that("facs_plot_cluster_heatmap() rejects an invalid scale value", {
  expect_error(facs_plot_cluster_heatmap(heatmap_input(), scale = "NotAScale"))
})
