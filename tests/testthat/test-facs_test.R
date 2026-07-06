library(testthat)
library(JanisHelpers)

# Synthetic tibble shaped like facs_calc_cluster_freq()'s output: file_name x
# metacluster x n/fraction, plus mouse_ID/group passthrough columns. Cluster
# "1" gets a deliberate abundance shift for the non-reference group(s), so
# tests can assert it has the smallest p_val without pinning exact values.
freq_input <- function(groups = c("control", "control", "control", "control",
                                   "treated", "treated", "treated", "treated")) {
  mice  <- paste0("mouse", seq_along(groups))
  files <- paste0(mice, ".fcs")
  clusters <- factor(1:4)

  base <- tidyr::expand_grid(file_name = files, metacluster = clusters)
  set.seed(7)
  base$n <- stats::rpois(nrow(base), lambda = 50)

  shifted_files <- files[groups != groups[1]]
  base$n[base$metacluster == "1" & base$file_name %in% shifted_files] <-
    base$n[base$metacluster == "1" & base$file_name %in% shifted_files] + 40L

  lookup <- tibble::tibble(file_name = files, mouse_ID = mice, group = factor(groups))
  base <- dplyr::left_join(base, lookup, by = "file_name")

  totals <- base |> dplyr::group_by(file_name) |> dplyr::summarise(total = sum(n), .groups = "drop")
  base |>
    dplyr::left_join(totals, by = "file_name") |>
    dplyr::mutate(fraction = n / total) |>
    dplyr::select(!total)
}

test_that("facs_test_cluster_abundance() flags the shifted cluster with method='glmm'", {
  result <- facs_test_cluster_abundance(freq_input(), fixed = "group", method = "glmm")

  expect_true(all(c("metacluster", "contrast", "p_val", "p_adj") %in% names(result)))
  expect_equal(nrow(result), 4L)
  expect_equal(unique(result$contrast), "treated_vs_control")

  shifted_p <- result$p_val[result$metacluster == "1"]
  expect_true(shifted_p == min(result$p_val))
})

test_that("facs_test_cluster_abundance() supports random effects with method='glmm'", {
  result <- facs_test_cluster_abundance(
    freq_input(), fixed = "group", random = "mouse_ID", method = "glmm"
  )
  expect_true(all(c("metacluster", "contrast", "p_val", "p_adj") %in% names(result)))
  expect_equal(nrow(result), 4L)
})

test_that("facs_test_cluster_abundance() flags the shifted cluster with method='edgeR'", {
  result <- facs_test_cluster_abundance(freq_input(), fixed = "group", method = "edgeR")

  expect_true(all(c("metacluster", "contrast", "logFC", "logCPM", "p_val", "p_adj") %in% names(result)))
  shifted_p <- result$p_val[result$metacluster == "1"]
  expect_true(shifted_p == min(result$p_val))
})

test_that("facs_test_cluster_abundance() flags the shifted cluster with method='voom'", {
  result <- facs_test_cluster_abundance(freq_input(), fixed = "group", method = "voom")

  expect_true(all(c("metacluster", "contrast", "logFC", "t", "p_val", "p_adj") %in% names(result)))
  shifted_p <- result$p_val[result$metacluster == "1"]
  expect_true(shifted_p == min(result$p_val))
})

test_that("facs_test_cluster_abundance() tests one contrast per non-reference level", {
  groups <- rep(c("low", "mid", "high"), length.out = 8)
  result <- facs_test_cluster_abundance(
    freq_input(groups), fixed = "group", ref_level = "mid", method = "glmm"
  )

  expect_equal(nrow(result), 8L)
  expect_setequal(unique(result$contrast), c("low_vs_mid", "high_vs_mid"))
})

test_that("facs_test_cluster_abundance() defaults ref_level to the first factor level", {
  result <- facs_test_cluster_abundance(freq_input(), fixed = "group", method = "glmm")
  expect_equal(unique(result$contrast), "treated_vs_control")
})
