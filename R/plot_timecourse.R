#' Timecourse plot of each mouse with group medians
#'
#' One thin line per mouse (coloured by group) plus a thick line through the
#' group medians, over study weeks -- e.g. anti-MPO titer per bleed. Input is
#' long data with one row per sample (ELISA from [db_load_experiment()] after
#' [meta_label_weeks()], with `group` joined). Rows with a missing value
#' (e.g. `OOR<`, stored as `NA`) are dropped with a message, never
#' substituted; replicates are averaged per mouse and week first. A `lloq`
#' column (as in DB ELISA data) draws a dashed line when it holds a value.
#' Returns a `ggplot` to extend with `+ ggplot2::labs()` etc.
#'
#' @param data Long data frame.
#' @param x,y,fill,subject Bare column names: time (default `week_num`;
#'   numeric = real time, factor = evenly spaced categories), value (default
#'   `value`), group colour (default `group`) and mouse (default `mouse_id`).
#' @param facet Optional bare column name to facet by (e.g. `analyte`), free y
#'   scales.
#' @param scale `"linear"` (default) or `"log"` (log10 y; values must be
#'   positive).
#' @param colors Optional named vector of group colours; defaults to the
#'   option `JanisHelpers.group_colors`, see [plot_box()]. `NULL` = ggplot2
#'   defaults.
#' @param p_val If `TRUE`, compares the groups at each time point (Wilcoxon
#'   rank-sum for two groups, Kruskal-Wallis for more; time points where a
#'   group has fewer than 2 mice are skipped), adjusts the p-values across
#'   time points (Benjamini-Hochberg, per facet) and prints stars (`*` < 0.05,
#'   `**` < 0.01, `***` < 0.001) above significant time points.
#' @param min_n Minimum number of mice for a time point to be part of the
#'   group median line (default 2); the line bridges time points with fewer
#'   mice (e.g. one early death), whose mice are still drawn.
#' @returns A `ggplot` object.
#' @examples
#' d <- expand.grid(group = c("a", "b"), m = 1:3, week_num = c(-1, 1, 3))
#' d$mouse_id <- paste0(d$group, d$m)
#' d$value <- d$week_num + d$m + ifelse(d$group == "a", 5, 0)
#' plot_timecourse(d)
#' @seealso [plot_timecourse_summary()] for medians and IQR only.
#' @export
plot_timecourse <- function(data, x = week_num, y = value, fill = group, subject = mouse_id,
                            facet = NULL, scale = c("linear", "log"),
                            colors = getOption("JanisHelpers.group_colors"), p_val = FALSE,
                            min_n = 2) {
  scale <- match.arg(scale)
  facet_name <- facet_name_(substitute(facet))
  x_name <- deparse(substitute(x))
  d <- timecourse_prep_(data, {{ x }}, {{ y }}, {{ fill }}, {{ subject }}, facet_name)
  check_log_values_(d, scale)
  med <- dplyr::summarise(d, .y = stats::median(.data$.y), n = dplyr::n_distinct(.data$.subject),
                          .by = dplyr::any_of(c(".fill", ".x", ".facet"))) |>
    dplyr::filter(.data$n >= min_n)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$.x, y = .data$.y, colour = .data$.fill)) +
    ggplot2::geom_line(ggplot2::aes(group = .data$.subject), alpha = 0.4, linewidth = 0.4) +
    ggplot2::geom_point(alpha = 0.4, size = 1.2) +
    ggplot2::geom_line(data = med, ggplot2::aes(group = .data$.fill), linewidth = 1.2) +
    ggplot2::geom_point(data = med, size = 3) +
    ggplot2::labs(x = NULL, y = deparse(substitute(y)), colour = deparse(substitute(fill)))
  timecourse_finish_(p, d, data, facet_name, x_name, scale, colors, p_val)
}

#' Timecourse plot of group medians and IQR
#'
#' The group summary of [plot_timecourse()]: per group and time point the
#' median of the mice with an error bar from the 25th to the 75th percentile
#' (IQR), medians joined by a line, groups dodged side by side, and the number
#' of mice per group printed at the bottom of each time point. No individual
#' mice are drawn. Median and IQR suit small, skewed groups and match the
#' rank-based tests of `p_val`. A group with fewer than `min_n` mice at a time
#' point (e.g. one early death) is drawn as a hollow point without error bar
#' and left out of the median line, which bridges it. Same data handling as
#' [plot_timecourse()] (missing values dropped, replicates averaged per mouse).
#'
#' @inheritParams plot_timecourse
#' @param min_n Minimum number of mice for a time point to get an error bar and
#'   be part of the median line (default 2); fewer = hollow point.
#' @returns A `ggplot` object (caption `"median and IQR; n per group below"`,
#'   override with `+ ggplot2::labs(caption = ...)`).
#' @examples
#' d <- expand.grid(group = c("a", "b"), m = 1:3, week_num = c(-1, 1, 3))
#' d$mouse_id <- paste0(d$group, d$m)
#' d$value <- d$week_num + d$m + ifelse(d$group == "a", 5, 0)
#' plot_timecourse_summary(d)
#' @export
plot_timecourse_summary <- function(data, x = week_num, y = value, fill = group,
                                    subject = mouse_id, facet = NULL,
                                    scale = c("linear", "log"),
                                    colors = getOption("JanisHelpers.group_colors"),
                                    p_val = FALSE, min_n = 2) {
  scale <- match.arg(scale)
  facet_name <- facet_name_(substitute(facet))
  x_name <- deparse(substitute(x))
  d <- timecourse_prep_(data, {{ x }}, {{ y }}, {{ fill }}, {{ subject }}, facet_name)
  check_log_values_(d, scale)
  # categories -> positions 1..n, so every layer can share one numeric offset
  x_labels <- NULL
  if (!is.numeric(d$.x)) {
    x_labels <- levels(droplevels(as.factor(d$.x)))
    d$.x <- match(as.character(d$.x), x_labels)
  }
  s <- timecourse_summary_(d, min_n)
  # one sideways offset per group for all layers (position_dodge works per
  # layer and would stack groups split across the ok / low-n layers)
  groups <- levels(droplevels(as.factor(s$.fill)))
  k <- length(groups)
  width <- 0.4
  s$.xd <- s$.x + (match(as.character(s$.fill), groups) - (k + 1) / 2) * width / k
  ok <- s[!s$low_n, ]

  p <- ggplot2::ggplot(s, ggplot2::aes(x = .data$.xd, y = .data$median, colour = .data$.fill,
                                       group = .data$.fill)) +
    ggplot2::geom_errorbar(data = ok, ggplot2::aes(ymin = .data$q25, ymax = .data$q75),
                           width = 0.8 * width / k) +
    ggplot2::geom_line(data = ok, linewidth = 0.8) +
    ggplot2::geom_point(data = ok, size = 2.5)
  if (any(s$low_n)) {
    p <- p + ggplot2::geom_point(data = s[s$low_n, ], size = 2.5,
                                 shape = 21, fill = "white", show.legend = FALSE)
  }
  p <- p +
    ggplot2::geom_text(data = s, ggplot2::aes(y = -Inf, label = .data$n), vjust = -0.5,
                       size = 3, show.legend = FALSE) +
    ggplot2::labs(x = NULL, y = deparse(substitute(y)), colour = deparse(substitute(fill)),
                  caption = "median and IQR; n per group below")
  timecourse_finish_(p, d, data, facet_name, x_name, scale, colors, p_val, x_labels)
}

# Internal: per group x time point (x facet): n mice, median, quartiles,
# low_n flag.
timecourse_summary_ <- function(d, min_n) {
  dplyr::summarise(
    d,
    n = dplyr::n_distinct(.data$.subject),
    median = stats::median(.data$.y),
    q25 = unname(stats::quantile(.data$.y, 0.25, type = 7)),
    q75 = unname(stats::quantile(.data$.y, 0.75, type = 7)),
    .by = dplyr::any_of(c(".fill", ".x", ".facet"))
  ) |>
    dplyr::mutate(low_n = .data$n < min_n)
}

# Internal: name of the facet column from substitute(facet), or NULL.
facet_name_ <- function(expr) {
  if (is.null(expr)) return(NULL)
  if (is.character(expr)) return(expr)
  deparse(expr)
}

# Internal: standard columns .x, .y, .fill, .subject (.facet), missing values
# dropped with a message, replicates averaged per subject x time (x facet).
timecourse_prep_ <- function(data, x, y, fill, subject, facet_name = NULL) {
  d <- dplyr::mutate(data, .x = {{ x }}, .y = {{ y }}, .fill = {{ fill }},
                     .subject = {{ subject }})
  if (!is.null(facet_name)) d$.facet <- d[[facet_name]]
  # averaging is for replicates only: refuse to mix analytes / tissues / ...
  cell <- intersect(c(".subject", ".x", ".facet"), names(d))
  for (col in setdiff(intersect(c("analyte", "tissue", "condition", "population", "metric"),
                                names(d)), facet_name)) {
    n_val <- dplyr::summarise(d, n = dplyr::n_distinct(.data[[col]]),
                              .by = dplyr::all_of(cell))$n
    if (any(n_val > 1)) {
      stop("several ", col, " values per mouse and time point -- filter `data` ",
           "or use `facet = ", col, "`.", call. = FALSE)
    }
  }
  n_na <- sum(is.na(d$.y) | is.na(d$.x))
  if (n_na > 0) {
    message("Dropped ", n_na, " row(s) with missing values (value or time point, ",
            "e.g. below LLOQ).")
  }
  d |>
    dplyr::filter(!is.na(.data$.y), !is.na(.data$.x)) |>
    dplyr::summarise(.y = mean(.data$.y),
                     .by = dplyr::any_of(c(".subject", ".x", ".fill", ".facet")))
}

# Internal: log scale needs positive values.
check_log_values_ <- function(d, scale) {
  if (scale == "log" && any(d$.y <= 0)) {
    stop("scale = \"log\" needs positive values; ", sum(d$.y <= 0),
         " value(s) are <= 0.", call. = FALSE)
  }
}

# Internal: axis, LLOQ line, y scale, colours, facets, stars, theme.
timecourse_finish_ <- function(p, d, data, facet_name, x_name, scale, colors, p_val,
                               x_labels = NULL) {
  if (is.numeric(d$.x)) {
    br <- sort(unique(d$.x))
    labs <- if (!is.null(x_labels)) x_labels[br] else
      if (grepl("week", x_name)) paste("week", br) else format(br)
    p <- p + ggplot2::scale_x_continuous(breaks = br, labels = labs, minor_breaks = NULL)
  }
  lloq <- timecourse_lloq_(data, facet_name)
  if (!is.null(lloq)) {
    p <- p + ggplot2::geom_hline(data = lloq, ggplot2::aes(yintercept = .data$lloq),
                                 linetype = "dashed", colour = "grey40")
  }
  if (p_val) {
    stars <- timecourse_stars_(d, scale)
    if (nrow(stars) > 0) {
      p <- p + ggplot2::geom_text(
        data = stars, ggplot2::aes(x = .data$.x, y = .data$.y_pos, label = .data$label),
        inherit.aes = FALSE, size = 5
      )
    }
  }
  if (!is.null(facet_name)) p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$.facet), scales = "free_y")
  p <- add_y_scale_(p + ggplot2::theme_bw(), scale)
  add_fill_scale_(p, colors, d$.fill, aesthetics = c("colour", "fill"))
}

# Internal: one LLOQ per facet (max if several) from a `lloq` column, or NULL.
timecourse_lloq_ <- function(data, facet_name) {
  if (!"lloq" %in% names(data) || all(is.na(data$lloq))) return(NULL)
  d <- dplyr::filter(data, !is.na(.data$lloq))
  if (!is.null(facet_name)) d$.facet <- d[[facet_name]]
  dplyr::summarise(d, lloq = max(.data$lloq), .by = dplyr::any_of(".facet"))
}

# Internal: per time point (x facet) group test, BH-adjusted per facet;
# returns significant points with star labels and a y position.
timecourse_stars_ <- function(d, scale) {
  by <- intersect(c(".facet", ".x"), names(d))
  tests <- d |>
    dplyr::summarise(
      p = timecourse_test_(.data$.y, .data$.fill),
      .y_pos = max(.data$.y) * if (scale == "log") 1.3 else 1.05,
      .by = dplyr::all_of(by)
    ) |>
    dplyr::filter(!is.na(.data$p)) |>
    dplyr::mutate(p_adj = stats::p.adjust(.data$p, method = "BH"),
                  .by = dplyr::any_of(".facet")) |>
    dplyr::mutate(label = dplyr::case_when(
      .data$p_adj < 0.001 ~ "***", .data$p_adj < 0.01 ~ "**",
      .data$p_adj < 0.05 ~ "*", TRUE ~ NA_character_
    ))
  dplyr::filter(tests, !is.na(.data$label))
}

# Internal: p-value of one time point, NA if fewer than 2 groups or a group
# with fewer than 2 values.
timecourse_test_ <- function(y, g) {
  g <- droplevels(as.factor(g))
  n <- table(g)
  if (length(n) < 2 || any(n < 2)) return(NA_real_)
  if (length(n) == 2) stats::wilcox.test(y ~ g, exact = FALSE)$p.value
  else stats::kruskal.test(y ~ g)$p.value
}
