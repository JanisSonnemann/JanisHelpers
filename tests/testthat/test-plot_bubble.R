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
