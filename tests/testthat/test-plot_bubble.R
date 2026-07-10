test_that("calc_bubble_fc_() errors when a required column is missing", {
  data <- tibble::tibble(group = c("control", "treated"), value = c(1, 2))
  expect_error(
    calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value",
      test = "auto", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5
    ),
    "population"
  )
})

test_that("calc_bubble_fc_() errors when control is not present in group_col", {
  data <- tibble::tibble(
    group = c("a", "b"), population = c("CD4", "CD4"), value = c(1, 2)
  )
  expect_error(
    calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value",
      test = "auto", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5
    ),
    "control"
  )
})

test_that("calc_bubble_fc_() errors when metric has more than one distinct value", {
  data <- tibble::tibble(
    group = c("control", "treated"), population = c("CD4", "CD4"),
    value = c(1, 2), metric = c("Count", "FractionOfParent")
  )
  expect_error(
    calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value",
      test = "auto", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5
    ),
    "metric"
  )
})

test_that("calc_bubble_fc_() errors when tissue has more than one distinct value", {
  data <- tibble::tibble(
    group = c("control", "treated"), population = c("CD4", "CD4"),
    value = c(1, 2), tissue = c("Spleen", "Blood")
  )
  expect_error(
    calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value",
      test = "auto", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5
    ),
    "tissue"
  )
})

test_that("calc_bubble_fc_() errors on an invalid test argument", {
  data <- tibble::tibble(
    group = c("control", "treated"), population = c("CD4", "CD4"), value = c(1, 2)
  )
  expect_error(
    calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value",
      test = "not-a-test", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5
    ),
    "test"
  )
})

make_data_ <- function(...) {
  groups <- list(...)
  purrr::imap_dfr(groups, function(values, group_pop) {
    parts <- strsplit(group_pop, "\\.")[[1]]
    tibble::tibble(group = parts[1], population = parts[2], value = values)
  })
}

test_that("calc_bubble_fc_() computes log2 fold-change per population", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    treated.CD4 = c(19.6, 19.9, 20.1, 20.4),
    control.CD8 = c(4.8, 4.9, 5.1, 5.2),
    treated.CD8 = c(4.85, 4.95, 5.05, 5.15)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "wilcox", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0
  )
  expect_equal(result$log2fc[result$population == "CD4"], 1, tolerance = 1e-8)
  expect_equal(result$log2fc[result$population == "CD8"], 0, tolerance = 1e-8)
})

test_that("calc_bubble_fc_() orders population/comparison by first appearance and excludes control", {
  data <- make_data_(
    control.CD8 = c(4.8, 4.9, 5.1, 5.2),
    treated.CD8 = c(4.85, 4.95, 5.05, 5.15),
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    treated.CD4 = c(19.6, 19.9, 20.1, 20.4)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "wilcox", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0
  )
  expect_equal(levels(result$population), c("CD8", "CD4"))
  expect_true(is.factor(result$comparison))
  expect_false("control" %in% as.character(result$comparison))
})

test_that("calc_bubble_fc_() warns and computes a finite log2fc for <2 observations in one arm", {
  data <- tibble::tibble(
    group = c("control", "treated", "treated", "treated"),
    population = c("CD3", "CD3", "CD3", "CD3"),
    value = c(10, 19.6, 19.9, 20.1)
  )
  expect_warning(
    result <- calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value",
      test = "wilcox", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5
    ),
    "Fewer than 2 observations"
  )
  expect_true(is.finite(result$log2fc))
})

test_that("calc_bubble_fc_() computes a wilcox p-value for a 2-group comparison", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    treated.CD4 = c(19.6, 19.9, 20.1, 20.4)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "wilcox", p_adjust_method = "none", summary_fun = mean, pseudocount = 0
  )
  expected_p <- stats::wilcox.test(
    data$value[data$group == "control"], data$value[data$group == "treated"]
  )$p.value
  expect_equal(result$p_value, expected_p, tolerance = 1e-8)
})

test_that("calc_bubble_fc_() uses stats::t.test when test = \"t.test\"", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    treated.CD4 = c(19.6, 19.9, 20.1, 20.4)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "t.test", p_adjust_method = "none", summary_fun = mean, pseudocount = 0
  )
  expected_p <- stats::t.test(
    data$value[data$group == "control"], data$value[data$group == "treated"]
  )$p.value
  expect_equal(result$p_value, expected_p, tolerance = 1e-8)
})

test_that("calc_bubble_fc_() adjusts p-values per comparison column", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    treated.CD4 = c(19.6, 19.9, 20.1, 20.4),
    control.CD8 = c(4.8, 4.9, 5.1, 5.2),
    treated.CD8 = c(4.85, 4.95, 5.05, 5.15)
  )
  raw <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "wilcox", p_adjust_method = "none", summary_fun = mean, pseudocount = 0
  )
  adjusted <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "wilcox", p_adjust_method = "bonferroni", summary_fun = mean, pseudocount = 0
  )
  expect_equal(
    sort(adjusted$p_value),
    sort(stats::p.adjust(raw$p_value, method = "bonferroni")),
    tolerance = 1e-8
  )
})

test_that("calc_bubble_fc_() assigns stars based on p-value thresholds", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    treated.CD4 = c(119.6, 119.9, 120.1, 120.4),
    control.CD8 = c(10, 10, 10, 10.001),
    treated.CD8 = c(10, 10.001, 10, 10)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "t.test", p_adjust_method = "none", summary_fun = mean, pseudocount = 0
  )
  expect_equal(result$stars[result$population == "CD4"], "***")
  expect_equal(result$stars[result$population == "CD8"], "")
})

test_that("calc_bubble_fc_() leaves p_value/stars NA/blank for <2-observation cells", {
  data <- tibble::tibble(
    group = c("control", "treated", "treated", "treated"),
    population = c("CD3", "CD3", "CD3", "CD3"),
    value = c(10, 19.6, 19.9, 20.1)
  )
  result <- suppressWarnings(calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "wilcox", p_adjust_method = "none", summary_fun = mean, pseudocount = 0.5
  ))
  expect_true(is.na(result$p_value))
  expect_equal(result$stars, "")
})

test_that("calc_bubble_fc_() dispatches to kruskal+dunn when test = \"auto\" and there are >2 groups", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    low.CD4 = c(14.6, 14.9, 15.1, 15.4),
    high.CD4 = c(19.6, 19.9, 20.1, 20.4)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "auto", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0
  )
  expect_equal(sort(as.character(result$comparison)), c("high", "low"))
  expect_true(all(!is.na(result$p_value)))
  expect_equal(result$log2fc[result$comparison == "low"], log2(1.5), tolerance = 1e-8)
  expect_equal(result$log2fc[result$comparison == "high"], log2(2), tolerance = 1e-8)
})

test_that("calc_bubble_fc_() forces kruskal+dunn even with only 2 groups when test = \"kruskal\"", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    treated.CD4 = c(19.6, 19.9, 20.1, 20.4)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "kruskal", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0
  )
  expect_false(is.na(result$p_value))
})

test_that("calc_bubble_fc_() forces pairwise wilcox even with >2 groups when test = \"wilcox\"", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    low.CD4 = c(14.6, 14.9, 15.1, 15.4),
    high.CD4 = c(19.6, 19.9, 20.1, 20.4)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value",
    test = "wilcox", p_adjust_method = "none", summary_fun = mean, pseudocount = 0
  )
  expected_low_p <- stats::wilcox.test(
    data$value[data$group == "control"], data$value[data$group == "low"]
  )$p.value
  expect_equal(result$p_value[result$comparison == "low"], expected_low_p, tolerance = 1e-8)
})
