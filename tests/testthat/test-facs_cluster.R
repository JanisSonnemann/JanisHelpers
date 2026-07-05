library(testthat)
library(JanisHelpers)

wsp_path <- testthat::test_path("../fixtures/Treg.wsp")
fcs_dir  <- testthat::test_path("../fixtures/Treg")
skip_msg <- "Treg fixture not available"
cd45_gate <- "Singlets/Lymphocytes/live/CD45+"

# 3 markers (not 2) so the "explicit markers override" test below can
# restrict to a genuine subset -- see the n_metaclusters comment further
# down for why grid/n_metaclusters values avoid n_metaclusters = 2.
cluster_input <- function(markers = c("CD4", "CD45", "TCRb"), max_events = 200, seed = 1) {
  facs_read_fcs_gated(
    wsp_path   = wsp_path,
    gate_path  = cd45_gate,
    markers    = markers,
    keywords   = c("mouse_ID", "tissue"),
    max_events = max_events,
    seed       = seed
  )
}

test_that("facs_cluster_flowsom() appends cluster and metacluster columns", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- cluster_input()
  result <- facs_cluster_flowsom(
    dat,
    grid_xdim = 3, grid_ydim = 3, n_metaclusters = 3, seed = 1
  )

  expect_true(all(c("cluster", "metacluster") %in% names(result)))
  expect_type(result$cluster, "integer")
  expect_s3_class(result$metacluster, "factor")
  expect_true(all(result$cluster >= 1L & result$cluster <= 9L))
  expect_equal(nrow(result), nrow(dat))
})

test_that("facs_cluster_flowsom() defaults to clustering on every double column", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  # n_metaclusters = 3, not 2: the installed ConsensusClusterPlus (1.72.0,
  # current Bioc 3.21 release, a FlowSOM::FlowSOM() dependency) throws
  # "Argument der Laenge 0" / "argument of length 0" from its internal
  # clusterTrackingPlot() whenever maxK = 2 -- its consensus-color matrix
  # only ever accumulates one row in that case, which R's default `[`
  # silently drops from a matrix to a bare vector, and nrow() on that
  # vector is NULL. Reproduces independent of FlowSOM/our code with e.g.
  # ConsensusClusterPlus::ConsensusClusterPlus(matrix(rnorm(40), 10), maxK = 2, ...).
  # n_metaclusters = 3 (grid still 2x2 = 4 nodes) avoids the bug.
  dat <- cluster_input()
  result <- facs_cluster_flowsom(
    dat,
    grid_xdim = 2, grid_ydim = 2, n_metaclusters = 3, seed = 1
  )
  expect_true(all(c("CD4", "CD45", "TCRb", "cluster", "metacluster") %in% names(result)))
})

test_that("facs_cluster_flowsom() clusters on an explicit markers override", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  # markers = c("CD4", "CD45") (a 2-out-of-3 subset), not a single marker:
  # FlowSOM::BuildSOM() does `fsom$data[, colsToUse]` with no `drop = FALSE`
  # (see FlowSOM 2.16.0 source), so a single-column `colsToUse` collapses
  # the matrix to a bare vector and the subsequent `seq_len(ncol(data))`
  # inside SOM() fails on that vector's NULL ncol(). This is independent of
  # n_metaclusters -- reproduces with any grid/n_metaclusters as long as
  # exactly one marker is selected. A 2-marker override still exercises the
  # "restrict to a subset of `data`'s columns" behavior this test targets.
  # n_metaclusters = 3 for the same reason as the test above.
  dat <- cluster_input()
  result <- facs_cluster_flowsom(
    dat, markers = c("CD4", "CD45"),
    grid_xdim = 2, grid_ydim = 2, n_metaclusters = 3, seed = 1
  )
  expect_true(all(c("cluster", "metacluster") %in% names(result)))
})

test_that("facs_cluster_flowsom() errors when a marker override is not found in data", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- cluster_input()
  expect_error(
    facs_cluster_flowsom(dat, markers = c("CD4", "NotAColumn")),
    "NotAColumn"
  )
})

test_that("facs_cluster_flowsom() errors when a marker column contains NA", {
  dat <- tibble::tibble(
    file_name = c("a.fcs", "a.fcs", "b.fcs"),
    CD4       = c(1.1, NA, 2.2),
    CD45      = c(3.3, 4.4, 5.5)
  )
  expect_error(
    facs_cluster_flowsom(dat, grid_xdim = 2, grid_ydim = 2, n_metaclusters = 2),
    "CD4"
  )
})

test_that("facs_cluster_flowsom() errors when n_metaclusters exceeds the grid's node count", {
  dat <- tibble::tibble(
    file_name = c("a.fcs", "b.fcs"),
    CD4       = c(1.1, 2.2),
    CD45      = c(3.3, 4.4)
  )
  expect_error(
    facs_cluster_flowsom(dat, grid_xdim = 1, grid_ydim = 1, n_metaclusters = 5),
    "n_metaclusters"
  )
})

test_that("facs_cluster_flowsom() is reproducible with the same seed", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- cluster_input()
  result_a <- facs_cluster_flowsom(dat, grid_xdim = 3, grid_ydim = 3, n_metaclusters = 3, seed = 42)
  result_b <- facs_cluster_flowsom(dat, grid_xdim = 3, grid_ydim = 3, n_metaclusters = 3, seed = 42)

  expect_identical(result_a$cluster, result_b$cluster)
  expect_identical(result_a$metacluster, result_b$metacluster)
})
