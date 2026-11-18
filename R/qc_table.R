# Internal QC helpers shared by the ELISA and histology reports. `data` is a
# long table with `mouse_ID`, `tissue`, `collected_at` (Date; NA = end),
# `value` and a variable column `var` ("analyte" for ELISA, "metric" for
# histology).

# Timepoint label: the date, or "end" for the end-point sample (NA date).
qc_timepoint_ <- function(collected_at) {
  ifelse(is.na(collected_at), "end", format(as.Date(collected_at)))
}

# Timepoint labels in display order: dates ascending, then "end".
qc_timepoint_levels_ <- function(timepoint) {
  u <- unique(timepoint)
  c(sort(setdiff(u, "end")), intersect("end", u))
}

# Missing expected mice per tissue x timepoint (comma-joined `missing`);
# attribute `unknown` = mice in `data` that are not expected.
qc_table_completeness_ <- function(data, expected_mice) {
  data <- dplyr::mutate(data, timepoint = qc_timepoint_(collected_at),
                        mouse_ID = as.character(mouse_ID))
  expected <- as.character(expected_mice)
  out <- data |>
    dplyr::summarise(
      missing = paste(setdiff(expected, mouse_ID), collapse = ", "),
      .by = c(tissue, timepoint)
    ) |>
    dplyr::arrange(tissue, match(timepoint, qc_timepoint_levels_(timepoint)))
  attr(out, "unknown") <- sort(setdiff(unique(data$mouse_ID), expected))
  out
}

# Values and empty (NA) values per `var` x tissue x timepoint.
qc_table_empty_ <- function(data, var) {
  data |>
    dplyr::mutate(timepoint = qc_timepoint_(collected_at)) |>
    dplyr::summarise(
      n_values = dplyr::n(), n_empty = sum(is.na(value)),
      .by = dplyr::all_of(c(var, "tissue", "timepoint"))
    ) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(c(var, "tissue"))),
                   match(timepoint, qc_timepoint_levels_(timepoint)))
}

# `data` plus `timepoint`, `t` (value on the plot axis: log10 if `log_scale`,
# NA for values <= 0 then), `side` and `flag` (outside the boxplot whiskers
# within group x tissue x timepoint x var; no group column -> no flags).
qc_table_score_ <- function(data, var, group_col = NULL, log_scale = FALSE) {
  d <- data |>
    dplyr::mutate(
      timepoint = qc_timepoint_(collected_at),
      t = if (isTRUE(log_scale)) {
        dplyr::if_else(!is.na(value) & value > 0, log10(value), NA_real_)
      } else {
        as.numeric(value)
      }
    )
  if (is.null(group_col)) {
    return(dplyr::mutate(d, side = NA_character_, flag = FALSE))
  }
  d |>
    dplyr::group_by(dplyr::across(dplyr::all_of(
      c(group_col, "tissue", "timepoint", var)
    ))) |>
    dplyr::mutate(side = qc_whisker_side_(t), flag = !is.na(side)) |>
    dplyr::ungroup()
}

# Boxplot + jitter of scored values; facets `var` (rows) x timepoint (columns,
# plus tissue if several); flagged points labelled with `mouse_ID`. x = group,
# or `batch_col` (points coloured by it) for the batch view.
qc_table_plot_ <- function(scored, var, group_col = NULL, log_scale = FALSE,
                           batch_col = NULL) {
  d <- scored[!is.na(scored$value) & (!isTRUE(log_scale) | scored$value > 0), ]
  if (nrow(d) == 0) {
    return(NULL)
  }
  d$.row <- as.character(d[[var]])
  d$.x <- if (!is.null(batch_col)) {
    qc_batch_label_(d[[batch_col]])
  } else if (is.null(group_col)) {
    "all"
  } else {
    qc_batch_label_(d[[group_col]])
  }
  d$.label <- dplyr::if_else(d$flag, as.character(d$mouse_ID), "")
  d$timepoint <- factor(d$timepoint, levels = qc_timepoint_levels_(d$timepoint))
  pos <- ggplot2::position_jitter(width = 0.15, height = 0, seed = 1)
  label_layer <- if (qc_has_ggrepel_()) {
    ggrepel::geom_text_repel(
      ggplot2::aes(label = .label), position = pos, size = 3, colour = "red3",
      min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE
    )
  } else {
    ggplot2::geom_text(
      ggplot2::aes(label = .label), position = pos, size = 3, colour = "red3",
      hjust = -0.2, show.legend = FALSE
    )
  }
  cols <- if (length(unique(d$tissue)) > 1) {
    ggplot2::vars(tissue, timepoint)
  } else {
    ggplot2::vars(timepoint)
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .x, y = value, fill = .x)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5)
  p <- if (is.null(batch_col)) {
    p + ggplot2::geom_point(position = pos, alpha = 0.7, size = 1.5)
  } else {
    p + ggplot2::geom_point(ggplot2::aes(colour = .x), position = pos, size = 1.5)
  }
  p <- p +
    label_layer +
    ggplot2::facet_grid(rows = ggplot2::vars(.row), cols = cols,
                        scales = "free_y") +
    ggplot2::labs(x = NULL, y = NULL,
                  fill = if (is.null(batch_col)) group_col else batch_col,
                  colour = batch_col) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  if (isTRUE(log_scale)) p <- p + ggplot2::scale_y_log10()
  p
}

# Per-mouse lines over `collected_at` (dated samples only), facets `var`;
# NULL if fewer than two distinct dates.
qc_table_time_plot_ <- function(data, var, group_col = NULL) {
  d <- data[!is.na(data$collected_at) & !is.na(data$value), ]
  if (length(unique(d$collected_at)) < 2) {
    return(NULL)
  }
  d$.row <- as.character(d[[var]])
  d$.group <- if (is.null(group_col)) "all" else qc_batch_label_(d[[group_col]])
  d$.line <- paste(d$mouse_ID, d$tissue)
  ggplot2::ggplot(d, ggplot2::aes(x = collected_at, y = value,
                                  group = .line, colour = .group)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::facet_grid(rows = ggplot2::vars(.row), scales = "free_y") +
    ggplot2::labs(x = NULL, y = NULL, colour = group_col) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))
}

# Flagged samples ranked by the number of distinct `var` values flagged:
# mouse_ID, n_flagged, <var>s, plus the flagged tissue x timepoint list.
qc_table_flag_table_ <- function(scored, var) {
  flagged <- scored[scored$flag, ]
  if (nrow(flagged) == 0) {
    return(tibble::tibble(mouse_ID = character(), n_flagged = integer(),
                          flagged = character(), where = character()))
  }
  flagged |>
    dplyr::summarise(
      n_flagged = dplyr::n_distinct(.data[[var]]),
      flagged = paste(sort(unique(.data[[var]])), collapse = ", "),
      where = paste(sort(unique(paste(tissue, timepoint))), collapse = ", "),
      .by = mouse_ID
    ) |>
    dplyr::arrange(dplyr::desc(n_flagged), mouse_ID)
}
