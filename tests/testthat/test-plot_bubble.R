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

test_that("calc_bubble_fc_() errors when x_col is set and group_col has more than 2 distinct values", {
  data <- tibble::tibble(
    group = c("control", "control", "low", "low", "high", "high"),
    population = "CD4",
    tissue = "Spleen",
    value = c(9.7, 9.9, 14.6, 14.9, 19.6, 19.9)
  )
  expect_error(
    calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value", x_col = "tissue",
      test = "auto", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5
    ),
    "exactly 2 distinct values"
  )
})

test_that("calc_bubble_fc_() errors when x_col is not present in data", {
  data <- tibble::tibble(
    group = c("control", "control", "treated", "treated"),
    population = "CD4",
    value = c(9.7, 9.9, 19.6, 19.9)
  )
  expect_error(
    calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value", x_col = "tissue",
      test = "auto", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5
    ),
    "not found in `data`"
  )
})

test_that("calc_bubble_fc_() errors when x_col equals group_col, population_col, or value_col", {
  data <- tibble::tibble(
    group = c("control", "control", "treated", "treated"),
    population = "CD4",
    value = c(9.7, 9.9, 19.6, 19.9)
  )
  expect_error(
    calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value", x_col = "group",
      test = "auto", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0.5
    ),
    "cannot be the same as"
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

test_that("calc_bubble_fc_() supports x_col = \"tissue\" for a 2-group comparison", {
  data <- tibble::tibble(
    group = c("control", "control", "treated", "treated",
              "control", "control", "treated", "treated"),
    population = "CD4",
    tissue = c("Spleen", "Spleen", "Spleen", "Spleen",
               "Blood", "Blood", "Blood", "Blood"),
    value = c(9.7, 10.3, 19.6, 20.4, 4.8, 5.2, 4.85, 5.15)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value", x_col = "tissue",
    test = "wilcox", p_adjust_method = "none", summary_fun = mean, pseudocount = 0
  )
  expect_true("tissue" %in% names(result))
  expect_false("comparison" %in% names(result))
  expect_equal(levels(result$tissue), c("Spleen", "Blood"))
  expect_equal(result$log2fc[result$tissue == "Spleen"], 1, tolerance = 1e-8)
  expect_equal(result$log2fc[result$tissue == "Blood"], 0, tolerance = 1e-8)
})

test_that("calc_bubble_fc_() supports x_col = \"metric\" for a 2-group comparison", {
  data <- tibble::tibble(
    group = c("control", "control", "treated", "treated",
              "control", "control", "treated", "treated"),
    population = "CD4",
    metric = c("Count", "Count", "Count", "Count",
               "FractionOfParent", "FractionOfParent", "FractionOfParent", "FractionOfParent"),
    value = c(9.7, 10.3, 19.6, 20.4, 4.8, 5.2, 4.85, 5.15)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value", x_col = "metric",
    test = "wilcox", p_adjust_method = "none", summary_fun = mean, pseudocount = 0
  )
  expect_true("metric" %in% names(result))
  expect_equal(levels(result$metric), c("Count", "FractionOfParent"))
  expect_equal(result$log2fc[result$metric == "Count"], 1, tolerance = 1e-8)
  expect_equal(result$log2fc[result$metric == "FractionOfParent"], 0, tolerance = 1e-8)
})

test_that("calc_bubble_fc_() lifts the tissue guard but keeps the metric guard when x_col = \"tissue\"", {
  data <- tibble::tibble(
    group = rep(c("control", "control", "treated", "treated"), 2),
    population = "CD4",
    tissue = rep(c("Spleen", "Blood"), each = 4),
    metric = c("Count", "Count", "Count", "Count",
               "FractionOfParent", "FractionOfParent", "FractionOfParent", "FractionOfParent"),
    value = c(9.7, 10.3, 19.6, 20.4, 4.8, 5.2, 4.85, 5.15)
  )
  expect_error(
    calc_bubble_fc_(
      data, control = "control", group_col = "group",
      population_col = "population", value_col = "value", x_col = "tissue",
      test = "wilcox", p_adjust_method = "none", summary_fun = mean, pseudocount = 0
    ),
    "metric"
  )
})

test_that("calc_bubble_fc_() applies p_adjust_method per x_value when x_col is set", {
  data <- tibble::tibble(
    group = rep(c("control", "control", "treated", "treated"), 2),
    population = rep(c("CD4", "CD8"), each = 4),
    tissue = "Spleen",
    value = c(9.7, 10.3, 19.6, 20.4, 4.8, 5.2, 4.85, 5.15)
  )
  raw <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value", x_col = "tissue",
    test = "wilcox", p_adjust_method = "none", summary_fun = mean, pseudocount = 0
  )
  adjusted <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value", x_col = "tissue",
    test = "wilcox", p_adjust_method = "bonferroni", summary_fun = mean, pseudocount = 0
  )
  expect_equal(
    sort(adjusted$p_value),
    sort(stats::p.adjust(raw$p_value, method = "bonferroni")),
    tolerance = 1e-8
  )
})

test_that("calc_bubble_fc_() forces kruskal+dunn per x_value when x_col is set", {
  data <- tibble::tibble(
    group = rep(c("control", "control", "treated", "treated"), 2),
    population = "CD4",
    tissue = rep(c("Spleen", "Blood"), each = 4),
    value = c(9.7, 9.9, 19.6, 19.9, 4.8, 4.9, 4.85, 4.95)
  )
  result <- calc_bubble_fc_(
    data, control = "control", group_col = "group",
    population_col = "population", value_col = "value", x_col = "tissue",
    test = "kruskal", p_adjust_method = "BH", summary_fun = mean, pseudocount = 0
  )
  expect_equal(levels(result$tissue), c("Spleen", "Blood"))
  expect_true(all(!is.na(result$p_value)))
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

  expected_dunn <- rstatix::dunn_test(
    data, value ~ group, p.adjust.method = "BH"
  ) |>
    dplyr::filter(group1 == "control" | group2 == "control")
  expected_low_p <- expected_dunn$p.adj[
    expected_dunn$group1 == "low" | expected_dunn$group2 == "low"
  ]
  expected_high_p <- expected_dunn$p.adj[
    expected_dunn$group1 == "high" | expected_dunn$group2 == "high"
  ]
  expect_equal(result$p_value[result$comparison == "low"], expected_low_p, tolerance = 1e-8)
  expect_equal(result$p_value[result$comparison == "high"], expected_high_p, tolerance = 1e-8)
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

test_that("plot_bubble_fc() returns a ggplot object", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    treated.CD4 = c(19.6, 19.9, 20.1, 20.4)
  )
  p <- plot_bubble_fc(data, control = "control")
  expect_s3_class(p, "ggplot")
})

test_that("plot_bubble_fc() captions the control group", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    treated.CD4 = c(19.6, 19.9, 20.1, 20.4)
  )
  p <- plot_bubble_fc(data, control = "control")
  expect_match(p$labels$caption, "control")
})

test_that("plot_bubble_fc() builds for a 3-group (kruskal) comparison without error", {
  data <- make_data_(
    control.CD4 = c(9.7, 9.9, 10.1, 10.3),
    low.CD4 = c(14.6, 14.9, 15.1, 15.4),
    high.CD4 = c(19.6, 19.9, 20.1, 20.4)
  )
  p <- plot_bubble_fc(data, control = "control")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot_bubble_fc() builds when some cells are untestable (NA p-value)", {
  data <- tibble::tibble(
    group = c("control", "treated", "treated", "treated"),
    population = c("CD3", "CD3", "CD3", "CD3"),
    value = c(10, 19.6, 19.9, 20.1)
  )
  p <- suppressWarnings(plot_bubble_fc(data, control = "control"))
  expect_no_error(ggplot2::ggplot_build(p))
})
