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
