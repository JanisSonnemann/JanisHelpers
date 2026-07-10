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

  resolved_test <- test
  if (test == "auto") {
    resolved_test <- if (length(group_levels) > 2L) "kruskal" else "wilcox"
  }

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

  p_values <- calc_bubble_p_values_(
    data, control, group_col, population_col, value_col,
    resolved_test, p_adjust_method, fc_table
  )

  fc_table |>
    dplyr::select(!testable) |>
    dplyr::left_join(p_values, by = c("population", "comparison")) |>
    dplyr::mutate(
      stars = dplyr::case_when(
        is.na(p_value) ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE ~ ""
      ),
      population = factor(population, levels = population_levels),
      comparison = factor(comparison, levels = comparison_levels)
    )
}

adjust_p_ <- function(p, method) {
  out <- rep(NA_real_, length(p))
  non_na <- !is.na(p)
  out[non_na] <- stats::p.adjust(p[non_na], method = method)
  out
}

calc_bubble_p_values_ <- function(data, control, group_col, population_col,
                                   value_col, resolved_test, p_adjust_method,
                                   fc_table) {
  if (resolved_test %in% c("wilcox", "t.test")) {
    raw <- fc_table |>
      purrr::pmap_dfr(function(population, comparison, log2fc, testable) {
        p_value <- NA_real_
        if (testable) {
          control_values <- data[[value_col]][
            data[[population_col]] == population & data[[group_col]] == control
          ]
          group_values <- data[[value_col]][
            data[[population_col]] == population & data[[group_col]] == comparison
          ]
          p_value <- if (resolved_test == "wilcox") {
            stats::wilcox.test(control_values, group_values)$p.value
          } else {
            stats::t.test(control_values, group_values)$p.value
          }
        }
        tibble::tibble(population = population, comparison = comparison, p_value = p_value)
      })

    raw |>
      dplyr::group_by(comparison) |>
      dplyr::mutate(p_value = adjust_p_(p_value, p_adjust_method)) |>
      dplyr::ungroup()
  } else if (resolved_test == "kruskal") {
    population_levels <- unique(fc_table$population)
    purrr::map_dfr(population_levels, function(pop) {
      pop_data <- data[data[[population_col]] == pop, ]
      group_counts <- table(pop_data[[group_col]])
      testable_groups <- names(group_counts)[group_counts >= 2L]
      if (!(control %in% testable_groups) || length(testable_groups) < 2L) {
        return(tibble::tibble(
          population = character(), comparison = character(), p_value = double()
        ))
      }
      pop_data_testable <- pop_data[pop_data[[group_col]] %in% testable_groups, ]
      dunn <- rstatix::dunn_test(
        pop_data_testable,
        stats::as.formula(paste(value_col, "~", group_col)),
        p.adjust.method = p_adjust_method
      )
      dunn |>
        dplyr::filter(group1 == control | group2 == control) |>
        dplyr::mutate(
          comparison = ifelse(group1 == control, group2, group1),
          population = pop
        ) |>
        dplyr::select(population, comparison, p_value = p.adj)
    })
  }
}
