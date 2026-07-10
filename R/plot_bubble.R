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

  NULL
}
