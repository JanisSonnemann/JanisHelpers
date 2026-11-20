#' Boxplot with jittered points per group
#'
#' Dodged boxplots (no outlier points) with the individual values jittered on
#' top -- the basic group-comparison plot of the analysis documents. Returns a
#' `ggplot` to extend with `+ ggplot2::labs()`, `+ ggplot2::facet_wrap()` etc.
#'
#' Rows flagged in an `excluded` column (from
#' `db_load_experiment(keep_excluded = TRUE)` or `db_query_*()`) are dropped
#' with a message, or with `show_excluded = TRUE` drawn as red crosses with a
#' red mouse-ID label. Excluded rows never enter the boxes (or p-values).
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
#'   once with `options(JanisHelpers.group_colors = c(...))`; the option is
#'   ignored when no fill value is one of its names (fill by condition, tissue,
#'   ...). `NULL` (option unset, or passed explicitly) = ggplot2 default fills.
#' @param show_excluded `FALSE` (default) drops rows with `excluded == TRUE`
#'   (message); `TRUE` draws them as red crosses labelled with their mouse ID.
#' @param label `TRUE` labels every point with its mouse ID (uses `ggrepel`
#'   when installed).
#' @param id Bare column name of the mouse ID used for labels (default
#'   `mouse_id`); only read when something is labelled.
#' @returns A `ggplot` object.
#' @examples
#' d <- data.frame(group = rep(c("a", "b"), each = 4), tissue = "kidney",
#'                 mouse_id = paste0("m", 1:8),
#'                 value = c(1, 2, 3, 4, 3, 4, 5, 6))
#' plot_box(d, x = tissue, y = value, fill = group, label = TRUE)
#' @export
plot_box <- function(data, x, y, fill, scale = c("linear", "percent", "log"),
                     colors = getOption("JanisHelpers.group_colors"),
                     show_excluded = FALSE, label = FALSE, id = mouse_id) {
  scale <- match.arg(scale)
  data <- plot_split_excluded_(data, show_excluded)
  any_excl <- any(data$.excluded)
  if (isTRUE(label) || any_excl) {
    data <- plot_add_labels_(data, dplyr::pull(data, {{ id }}), label)
  }
  kept <- data[!data$.excluded, , drop = FALSE]
  pos <- ggplot2::position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75,
                                       seed = 1)
  # group pinned to the fill: the .excluded colour / shape must not change the
  # dodge, and points and labels must share it (same rows, same jitter seed)
  pt_map <- if (any_excl) {
    ggplot2::aes(group = {{ fill }}, shape = .data$.excluded, colour = .data$.excluded)
  } else {
    ggplot2::aes(group = {{ fill }})
  }
  lab_map <- if (any_excl) {
    ggplot2::aes(group = {{ fill }}, label = .data$.lab, colour = .data$.excluded)
  } else {
    ggplot2::aes(group = {{ fill }}, label = .data$.lab)
  }
  p <- ggplot2::ggplot(kept, ggplot2::aes(x = {{ x }}, y = {{ y }}, fill = {{ fill }})) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.75),
                          outlier.shape = NA) +
    ggplot2::geom_point(data = data, mapping = pt_map, position = pos, alpha = 0.7)
  if (any_excl) p <- p + plot_excluded_scales_(c(`FALSE` = 16, `TRUE` = 4))
  p <- p + plot_label_layer_(data, lab_map, pos) +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2)) +
    ggplot2::theme_bw()
  p <- add_y_scale_(p, scale)
  add_fill_scale_(p, colors, dplyr::transmute(kept, .fill = {{ fill }})$.fill)
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
#' @inheritParams plot_box
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
                            pseudocount = 1, show_excluded = FALSE,
                            label = FALSE, id = mouse_id) {
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
                scale = scale, colors = colors, show_excluded = show_excluded,
                label = label, id = {{ id }}) +
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
#' One cluster of boxes per group, the conditions side by side inside it, and
#' one grey line per subject joining its values across conditions -- for wells
#' of the same sample that differ in one variable (restimulation). The x axis
#' shows the condition names (in two staggered rows); the group is shown by the
#' box colour (legend), conditions are shades of it (first level, the
#' baseline, lightest). Returns a `ggplot` to extend with `+ labs()` etc.
#'
#' @param data Long tibble, already filtered to the metric (and tissue) to
#'   plot, with the columns named below.
#' @param levels Character vector of conditions to keep, in x-axis order
#'   (e.g. `c("unstim", "MPO")`). Errors if one is absent from the data.
#' @param condition,fill,pair,facet,y Bare column names: condition (default
#'   `condition`, as in [db_query_facs()]), box fill / group (default
#'   `group`), subject joined by lines and used for labels (default
#'   `mouse_id`), facet variable (default `population`) and values (default
#'   `value`).
#' @param scale `"percent"` (default, for 0-1 fractions), `"linear"` or
#'   `"log"`, see [plot_box()].
#' @param colors Optional named fill colours; defaults to the option
#'   `JanisHelpers.group_colors`, see [plot_box()].
#' @param show_excluded `FALSE` (default) drops rows with `excluded == TRUE`
#'   (message); `TRUE` draws them as red crosses labelled with the subject
#'   (joined by a red dashed line), outside the boxes.
#' @param label `TRUE` labels every point with its subject (`pair`).
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
                        colors = getOption("JanisHelpers.group_colors"),
                        show_excluded = FALSE, label = FALSE) {
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
  data <- plot_split_excluded_(data, show_excluded)
  data <- dplyr::mutate(data, .cond = {{ condition }}, .grp = {{ fill }},
                        .pair = {{ pair }}, .y = {{ y }})

  # numeric x: one cluster per group, the conditions side by side inside it
  groups <- levels(droplevels(as.factor(data$.grp)))
  n_cond <- length(levels)
  step <- n_cond + 0.8
  data$.x <- (match(as.character(data$.grp), groups) - 1) * step +
    as.integer(data$.cond)
  any_excl <- any(data$.excluded)
  data <- plot_add_labels_(data, data$.pair, label)
  kept <- data[!data$.excluded, , drop = FALSE]
  excl <- data[data$.excluded, , drop = FALSE]
  excl_lines <- excl[excl$.pair %in% excl$.pair[duplicated(excl$.pair)], , drop = FALSE]
  breaks <- as.vector(outer(seq_len(n_cond), (seq_along(groups) - 1) * step, `+`))

  p <- ggplot2::ggplot(kept, ggplot2::aes(x = .data$.x, y = .data$.y, fill = .data$.grp)) +
    ggplot2::geom_boxplot(ggplot2::aes(group = .data$.x, alpha = .data$.cond),
                          outlier.shape = NA, width = 0.7) +
    ggplot2::geom_line(ggplot2::aes(group = .data$.pair),
                       colour = "grey40", alpha = 0.4, linewidth = 0.3)
  if (nrow(excl_lines) > 0) {
    p <- p + ggplot2::geom_line(data = excl_lines, ggplot2::aes(group = .data$.pair),
                                colour = "red3", linetype = "dashed", linewidth = 0.3)
  }
  p <- if (any_excl) {
    p + ggplot2::geom_point(
      data = data, ggplot2::aes(shape = .data$.excluded, colour = .data$.excluded),
      alpha = 0.8, size = 2
    ) + plot_excluded_scales_(c(`FALSE` = 21, `TRUE` = 4))
  } else {
    p + ggplot2::geom_point(data = data, shape = 21, colour = "black", alpha = 0.8, size = 2)
  }
  p <- p +
    plot_label_layer_(
      data,
      if (any_excl) {
        ggplot2::aes(label = .data$.lab, colour = .data$.excluded)
      } else {
        ggplot2::aes(label = .data$.lab)
      },
      "identity"
    ) +
    ggplot2::scale_x_continuous(breaks = breaks, labels = rep(levels, length(groups)),
                                minor_breaks = NULL,
                                guide = ggplot2::guide_axis(n.dodge = 2)) +
    ggplot2::scale_alpha_manual(
      values = stats::setNames(seq(0.3, 0.95, length.out = n_cond), levels),
      name = cond_col,
      guide = ggplot2::guide_legend(override.aes = list(fill = "grey30"))
    ) +
    ggplot2::facet_wrap(ggplot2::vars({{ facet }}), scales = "free_y", nrow = 2) +
    ggplot2::labs(x = NULL, y = names(dplyr::select(data, {{ y }})),
                  fill = names(dplyr::select(data, {{ fill }}))) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  p <- add_y_scale_(p, scale)
  add_fill_scale_(p, colors, kept$.grp)
}

# Internal: manual fill scale. The document-wide group colours
# (option JanisHelpers.group_colors) are skipped when none of the fill values
# is named in them -- the plot is filled by something other than the group
# (condition, tissue, ...), and a manual scale would grey out every box.
# `aesthetics` = c("colour", "fill") colours lines too (plot_timecourse()).
add_fill_scale_ <- function(p, colors, fill_values, aesthetics = "fill") {
  if (is.null(colors)) return(p)
  from_option <- identical(colors, getOption("JanisHelpers.group_colors"))
  if (from_option && !any(as.character(unique(fill_values)) %in% names(colors))) {
    return(p)
  }
  p + ggplot2::scale_discrete_manual(aesthetics = aesthetics, values = colors)
}

# Internal: y scale shared by plot_box() and plot_paired().
add_y_scale_ <- function(p, scale) {
  switch(scale,
    percent = p + ggplot2::scale_y_continuous(labels = scales::percent),
    log = p + ggplot2::scale_y_log10(guide = "axis_logticks", minor_breaks = 0),
    linear = p
  )
}

# Internal: logical `.excluded` (never NA) from an `excluded` column. Without
# show_excluded, excluded rows are dropped with a message.
plot_split_excluded_ <- function(data, show_excluded) {
  has_col <- "excluded" %in% names(data)
  data$.excluded <- if (has_col) data$excluded %in% TRUE else rep(FALSE, nrow(data))
  if (isTRUE(show_excluded) && !has_col) {
    message("show_excluded = TRUE, but `data` has no `excluded` column ",
            "(load with db_load_experiment(keep_excluded = TRUE)).")
  }
  n_excl <- sum(data$.excluded)
  if (!isTRUE(show_excluded) && n_excl > 0) {
    message("Dropped ", n_excl, " excluded row(s) ",
            "(show_excluded = TRUE draws them as red crosses).")
    data <- data[!data$.excluded, , drop = FALSE]
  }
  data
}

# Internal: label column `.lab` -- the ID for excluded rows always, for kept
# rows when `label`; "" = no label. Every row is kept, so the label layer has
# the same rows (and jitter) as the point layer.
plot_add_labels_ <- function(data, ids, label) {
  ids <- as.character(ids)
  data$.lab <- dplyr::if_else(data$.excluded | isTRUE(label), ids, "")
  data
}

# Internal: label layer over all rows of `data` (like the point layer, same
# `position`); NULL when nothing is labelled. With excluded rows `mapping`
# maps colour to `.excluded` (black / red3 via plot_excluded_scales_()),
# otherwise labels are black.
plot_label_layer_ <- function(data, mapping, position) {
  if (!".lab" %in% names(data) || all(data$.lab == "")) return(NULL)
  colour <- if (any(data$.excluded)) NULL else "black"
  args <- list(data = data, mapping = mapping, position = position, size = 3,
               show.legend = FALSE, colour = colour)
  args <- Filter(Negate(is.null), args)
  if (qc_has_ggrepel_()) {
    do.call(ggrepel::geom_text_repel,
            c(args, list(min.segment.length = 0, max.overlaps = Inf)))
  } else {
    do.call(ggplot2::geom_text, c(args, list(hjust = -0.2)))
  }
}

# Internal: shape + colour scales for `.excluded` (red cross); one merged
# legend key "excluded".
plot_excluded_scales_ <- function(shapes) {
  list(
    ggplot2::scale_shape_manual(values = shapes, breaks = "TRUE",
                                labels = "excluded", name = NULL),
    ggplot2::scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "red3"),
                                 breaks = "TRUE", labels = "excluded", name = NULL)
  )
}
