# Stats engine behind plot_bubble_fc(): computes per-population,
# per-comparison log2 fold-change and test statistics.
calc_bubble_fc_ <- function(data, control, group_col, population_col,
                             value_col, test, p_adjust_method,
                             summary_fun, pseudocount) {

  missing_cols <- setdiff(c(group_col, population_col, value_col), names(data))
  if (length(missing_cols) > 0L) {
    stop(glue::glue("`data` is missing column(s): {paste(missing_cols, collapse = ', ')}."))
  }
  if (!control %in% data[[group_col]]) {
    stop(glue::glue("`control` ('{control}') not found in `data${group_col}`."))
  }
  if ("metric" %in% names(data) && dplyr::n_distinct(data[["metric"]]) > 1L) {
    stop("`data` contains more than one distinct `metric`; pre-filter to a single metric.")
  }
  if ("tissue" %in% names(data) && dplyr::n_distinct(data[["tissue"]]) > 1L) {
    stop("`data` contains more than one distinct `tissue`; pre-filter to a single tissue.")
  }
  if (!test %in% c("auto", "wilcox", "t.test", "kruskal")) {
    stop('`test` must be one of "auto", "wilcox", "t.test", "kruskal".')
  }

  population_levels <- unique(data[[population_col]])
  group_levels <- unique(data[[group_col]])
  comparison_levels <- setdiff(group_levels, control)

  fc_table <- tidyr::expand_grid(
    population = population_levels,
    comparison = comparison_levels
  ) |>
    purrr::pmap_dfr(function(population, comparison) {
      control_values <- data[[value_col]][
        data[[population_col]] == population & data[[group_col]] == control
      ]
      group_values <- data[[value_col]][
        data[[population_col]] == population & data[[group_col]] == comparison
      ]
      log2fc <- log2(
        (summary_fun(group_values) + pseudocount) / (summary_fun(control_values) + pseudocount)
      )
      testable <- length(control_values) >= 2L && length(group_values) >= 2L
      if (!testable) {
        warning(glue::glue(
          "Fewer than 2 observations for population '{population}', ",
          "comparison '{comparison}' vs control '{control}'; p-value set to NA."
        ))
      }
      tibble::tibble(
        population = population, comparison = comparison,
        log2fc = log2fc, testable = testable
      )
    })

  fc_table |>
    dplyr::mutate(p_value = NA_real_, stars = "") |>
    dplyr::select(!testable) |>
    dplyr::mutate(
      population = factor(population, levels = population_levels),
      comparison = factor(comparison, levels = comparison_levels)
    )
}
