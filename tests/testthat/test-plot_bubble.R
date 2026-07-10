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
