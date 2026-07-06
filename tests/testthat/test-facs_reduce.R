library(testthat)
library(JanisHelpers)

wsp_path <- testthat::test_path("../fixtures/Treg.wsp")
fcs_dir  <- testthat::test_path("../fixtures/Treg")
skip_msg <- "Treg fixture not available"
cd45_gate <- "Singlets/Lymphocytes/live/CD45+"

reduce_input <- function(markers = c("CD4", "CD45", "TCRb"), max_events = 50, seed = 1) {
  facs_read_fcs_gated(
    wsp_path   = wsp_path,
    gate_path  = cd45_gate,
    markers    = markers,
    keywords   = c("mouse_ID", "tissue"),
    max_events = max_events,
    seed       = seed
  )
}

test_that("facs_reduce_umap() appends UMAP1 and UMAP2 columns", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- reduce_input()
  result <- facs_reduce_umap(dat, n_neighbors = 5, seed = 1)

  expect_true(all(c("UMAP1", "UMAP2") %in% names(result)))
  expect_type(result$UMAP1, "double")
  expect_type(result$UMAP2, "double")
  expect_equal(nrow(result), nrow(dat))
})

test_that("facs_reduce_umap() defaults to embedding on every double column", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- reduce_input()
  result <- facs_reduce_umap(dat, n_neighbors = 5, seed = 1)
  expect_true(all(c("CD4", "CD45", "TCRb", "UMAP1", "UMAP2") %in% names(result)))
})

test_that("facs_reduce_umap() embeds on an explicit markers override", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- reduce_input()
  result <- facs_reduce_umap(dat, markers = c("CD4", "CD45"), n_neighbors = 5, seed = 1)
  expect_true(all(c("UMAP1", "UMAP2") %in% names(result)))
})

test_that("facs_reduce_umap() errors when a marker override is not found in data", {
  dat <- tibble::tibble(
    file_name = c("a.fcs", "b.fcs"),
    CD4       = c(1.1, 2.2),
    CD45      = c(3.3, 4.4)
  )
  expect_error(
    facs_reduce_umap(dat, markers = c("CD4", "NotAColumn")),
    "NotAColumn"
  )
})

test_that("facs_reduce_umap() errors when an explicit marker is not double-typed", {
  dat <- tibble::tibble(
    file_name = c("a.fcs", "b.fcs"),
    CD4       = c(1.1, 2.2),
    mouse_ID  = c("M1", "M2")
  )
  expect_error(
    facs_reduce_umap(dat, markers = c("CD4", "mouse_ID")),
    "mouse_ID"
  )
})

test_that("facs_reduce_umap() errors when fewer than 2 markers resolve", {
  dat <- tibble::tibble(
    file_name = c("a.fcs", "b.fcs"),
    CD4       = c(1.1, 2.2)
  )
  expect_error(facs_reduce_umap(dat), "at least 2")
})

test_that("facs_reduce_umap() errors when a marker column contains NA", {
  dat <- tibble::tibble(
    file_name = c("a.fcs", "a.fcs", "b.fcs"),
    CD4       = c(1.1, NA, 2.2),
    CD45      = c(3.3, 4.4, 5.5)
  )
  expect_error(facs_reduce_umap(dat, n_neighbors = 2), "CD4")
})

test_that("facs_reduce_umap() errors when max_events is not a positive integer", {
  dat <- tibble::tibble(
    file_name = c("a.fcs", "b.fcs"),
    CD4       = c(1.1, 2.2),
    CD45      = c(3.3, 4.4)
  )
  expect_error(facs_reduce_umap(dat, max_events = -1), "positive integer")
  expect_error(facs_reduce_umap(dat, max_events = 1.5), "positive integer")
})

test_that("facs_reduce_umap() stratifies downsampling per file_name", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- reduce_input(max_events = 50)
  n_samples <- length(unique(dat$file_name))
  result <- facs_reduce_umap(dat, max_events = 12, n_neighbors = 5, seed = 1)

  expect_equal(nrow(result), 12L)
  counts <- table(result$file_name)
  expect_equal(length(counts), n_samples)
  expect_true(all(counts == 2L))
})

test_that("facs_reduce_umap() gives short samples all their rows without redistributing the shortfall", {
  dat <- tibble::tibble(
    file_name = c(rep("a.fcs", 2), rep("b.fcs", 10), rep("c.fcs", 10)),
    CD4       = rnorm(22),
    CD45      = rnorm(22)
  )
  result <- facs_reduce_umap(dat, max_events = 15, n_neighbors = 3, seed = 1)

  counts <- table(result$file_name)
  expect_equal(unname(counts[["a.fcs"]]), 2L)
  expect_lt(nrow(result), 15L)
})

test_that("facs_reduce_umap() drops samples entirely once max_events is smaller than n_samples", {
  file_names <- sprintf("s%02d.fcs", 1:20)
  dat <- tibble::tibble(
    file_name = file_names,
    CD4       = rnorm(20),
    CD45      = rnorm(20)
  )
  result <- facs_reduce_umap(dat, max_events = 16, n_neighbors = 5, seed = 1)

  expect_equal(nrow(result), 16L)
  expect_setequal(unique(result$file_name), file_names[1:16])
})

test_that("facs_reduce_umap() skips downsampling silently when max_events exceeds nrow", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- reduce_input(max_events = 20)
  result <- facs_reduce_umap(dat, max_events = 10000L, n_neighbors = 5, seed = 1)
  expect_equal(nrow(result), nrow(dat))
})

test_that("facs_reduce_umap() is reproducible with the same seed", {
  skip_if_not(dir.exists(fcs_dir), skip_msg)

  dat <- reduce_input(max_events = 50)
  r1 <- facs_reduce_umap(dat, max_events = 12, n_neighbors = 5, seed = 42)
  r2 <- facs_reduce_umap(dat, max_events = 12, n_neighbors = 5, seed = 42)

  expect_equal(r1$file_name, r2$file_name)
  expect_equal(r1$UMAP1, r2$UMAP1)
  expect_equal(r1$UMAP2, r2$UMAP2)
})

test_that("facs_reduce_umap() is reproducible with the same seed (synthetic data)", {
  dat <- tibble::tibble(
    file_name = rep(c("a.fcs", "b.fcs", "c.fcs"), each = 8),
    CD4       = rnorm(24),
    CD45      = rnorm(24),
    TCRb      = rnorm(24)
  )
  r1 <- facs_reduce_umap(dat, n_neighbors = 5, seed = 42)
  r2 <- facs_reduce_umap(dat, n_neighbors = 5, seed = 42)

  expect_equal(r1$UMAP1, r2$UMAP1)
  expect_equal(r1$UMAP2, r2$UMAP2)
})
