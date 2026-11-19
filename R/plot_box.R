#' Boxplot with jittered points per group
#'
#' Dodged boxplots (no outlier points) with the individual values jittered on
#' top -- the basic group-comparison plot of the analysis documents. Returns a
#' `ggplot` to extend with `+ ggplot2::labs()`, `+ ggplot2::facet_wrap()` etc.
#'
#' @param data A data frame.
#' @param x,y,fill Bare column names for the x axis, the y values and the box
#'   fill (e.g. treatment group).
#' @param scale `"linear"` (default), `"percent"` (y labels as percent, for
#'   0-1 fractions) or `"log"` (log10 y axis with log ticks; values must be
#'   positive -- no offset is added here, see [plot_population()]).
#' @param colors Optional named character vector of fill colours
#'   (`c(MPO409 = "#00BFC4", irrel = "#F8766D")`). Defaults to the option
#'   `JanisHelpers.group_colors`, so a document can set its group colours
#'   once with `options(JanisHelpers.group_colors = c(...))`. `NULL` (option
#'   unset, or passed explicitly) = ggplot2 default fills.
#' @returns A `ggplot` object.
#' @examples
#' d <- data.frame(group = rep(c("a", "b"), each = 4), tissue = "kidney",
#'                 value = c(1, 2, 3, 4, 3, 4, 5, 6))
#' plot_box(d, x = tissue, y = value, fill = group)
#' @export
plot_box <- function(data, x, y, fill, scale = c("linear", "percent", "log"),
                     colors = getOption("JanisHelpers.group_colors")) {
  scale <- match.arg(scale)
  p <- ggplot2::ggplot(data, ggplot2::aes(x = {{ x }}, y = {{ y }}, fill = {{ fill }})) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.75),
                          outlier.shape = NA) +
    ggplot2::geom_point(
      position = ggplot2::position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
      alpha = 0.7
    ) +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2)) +
    ggplot2::theme_bw()
  p <- add_y_scale_(p, scale)
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }
  p
}

#' Plot selected populations of long FACS data by group
#'
#' Filters long FACS data (one row per sample x population x metric, e.g.
#' collected [db_query_facs()] output joined with the subjects' `group`) to
#' some populations, metric(s) and tissues and draws [plot_box()], one facet
#' per population.
#'
#' @param data Long tibble with columns `population`, `metric`, `value`,
#'   `tissue` (when `tissues` is given) and the columns named in `x`, `fill`,
#'   `facet`.
#' @param populations Character vector of populations to keep.
#' @param metric Character vector of metrics to keep (several for an MFI
#'   panel faceted by `metric`).
#' @param tissues Optional character vector of tissues to keep; `NULL` = all.
#' @param x,fill,facet Bare column names: x axis (default `tissue`), box fill
#'   (default `group`, as returned by [db_query_subjects()]) and facet
#'   variable (default `population`).
#' @param scale `"linear"`, `"percent"` or `"log"`, see [plot_box()]. For
#'   `"log"`, `pseudocount` is added to every value first so zeros stay on the
#'   plot.
#' @param p_val `TRUE` adds `ggpubr::stat_compare_means()` p-values (Wilcoxon
#'   for two groups, Kruskal-Wallis for more) per facet.
#' @param colors Optional named fill colours; defaults to the option
#'   `JanisHelpers.group_colors`, see [plot_box()].
#' @param pseudocount Number added to `value` for `scale = "log"`. Default 1.
#' @returns A `ggplot` object.
#' @examples
#' \dontrun{
#' plot_population(ov, populations = c("CD4+", "CD8+"), metric = "count_per_g",
#'                 tissues = "kidney", x = group, scale = "log", p_val = TRUE)
#' }
#' @export
plot_population <- function(data, populations, metric, tissues = NULL,
                            x = tissue, fill = group, facet = population,
                            scale = c("linear", "percent", "log"),
                            p_val = FALSE,
                            colors = getOption("JanisHelpers.group_colors"),
                            pseudocount = 1) {
  scale <- match.arg(scale)
  if (!is.null(tissues)) {
    data <- dplyr::filter(data, .data$tissue %in% .env$tissues)
  }
  data <- dplyr::filter(
    data,
    .data$metric %in% .env$metric,
    .data$population %in% .env$populations
  )
  if (scale == "log") {
    data <- dplyr::mutate(data, value = .data$value + .env$pseudocount)
  }
  p <- plot_box(data, x = {{ x }}, y = value, fill = {{ fill }},
                scale = scale, colors = colors) +
    ggplot2::facet_wrap(ggplot2::vars({{ facet }}), scales = "free_y")
  if (isTRUE(p_val)) {
    # label built here, not via label = "p.format": that expands to ggpubr's
    # unexported create_p_label(), which only resolves when ggpubr is attached
    p <- p + ggpubr::stat_compare_means(
      ggplot2::aes(label = ggplot2::after_stat(paste0("p = ", p.format)))
    )
  }
  p
}

#' Paired plot of conditions per subject (e.g. stim vs unstim)
#'
#' Boxplots per condition x group with one grey line per subject joining its
#' values across conditions -- for wells of the same sample that differ in one
#' variable (restimulation). Returns a `ggplot` to extend with `+ labs()` etc.
#'
#' @param data Long tibble, already filtered to the metric (and tissue) to
#'   plot, with the columns named below.
#' @param levels Character vector of conditions to keep, in x-axis order
#'   (e.g. `c("unstim", "MPO")`). Errors if one is absent from the data.
#' @param condition,fill,pair,facet,y Bare column names: condition (default
#'   `condition`, as in [db_query_facs()]), box fill / group (default
#'   `group`), subject joined by lines (default `mouse_id`), facet variable
#'   (default `population`) and values (default `value`).
#' @param scale `"percent"` (default, for 0-1 fractions), `"linear"` or
#'   `"log"`, see [plot_box()].
#' @param colors Optional named fill colours; defaults to the option
#'   `JanisHelpers.group_colors`, see [plot_box()].
#' @returns A `ggplot` object.
#' @examples
#' \dontrun{
#' restim |>
#'   dplyr::filter(metric == "fraction_of_parent", tissue == "kidney") |>
#'   plot_paired(levels = c("unstim", "MPO")) +
#'   ggplot2::labs(y = "% of CD4+")
#' }
#' @export
plot_paired <- function(data, levels, condition = condition, fill = group,
                        pair = mouse_id, facet = population, y = value,
                        scale = c("percent", "linear", "log"),
                        colors = getOption("JanisHelpers.group_colors")) {
  scale <- match.arg(scale)
  cond_col <- names(dplyr::select(data, {{ condition }}))
  absent <- setdiff(levels, unique(as.character(data[[cond_col]])))
  if (length(absent) > 0) {
    stop(
      "Condition level(s) not in `", cond_col, "`: ",
      paste0("'", absent, "'", collapse = ", "), ".",
      call. = FALSE
    )
  }
  data <- data[as.character(data[[cond_col]]) %in% levels, , drop = FALSE]
  data[[cond_col]] <- factor(as.character(data[[cond_col]]), levels = levels)

  p <- ggplot2::ggplot(data, ggplot2::aes(
    x = interaction({{ condition }}, {{ fill }}), y = {{ y }}, fill = {{ fill }}
  )) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    ggplot2::geom_line(ggplot2::aes(group = {{ pair }}),
                       colour = "grey40", alpha = 0.4, linewidth = 0.3) +
    ggplot2::geom_point(shape = 21, colour = "black", alpha = 0.8, size = 2) +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2)) +
    ggplot2::facet_wrap(ggplot2::vars({{ facet }}), scales = "free_y", nrow = 2) +
    ggplot2::theme_bw()
  p <- add_y_scale_(p, scale)
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }
  p
}

# Internal: y scale shared by plot_box() and plot_paired().
add_y_scale_ <- function(p, scale) {
  switch(scale,
    percent = p + ggplot2::scale_y_continuous(labels = scales::percent),
    log = p + ggplot2::scale_y_log10(guide = "axis_logticks", minor_breaks = 0),
    linear = p
  )
}
