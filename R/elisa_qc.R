# Internal ELISA QC checks on the long table from elisa_read_curated().

# Per analyte: values (n), below / above the quantification range
# (OOR< / OOR> status) and numeric values under the LLOQ / over the ULOQ.
elisa_qc_range_ <- function(data) {
  data |>
    dplyr::summarise(
      n = dplyr::n(),
      n_oor_low = sum(result_status == "OOR<", na.rm = TRUE),
      n_oor_high = sum(result_status == "OOR>", na.rm = TRUE),
      n_below_lloq = sum(value < lloq, na.rm = TRUE),
      n_above_uloq = sum(value > uloq, na.rm = TRUE),
      .by = analyte
    ) |>
    dplyr::mutate(dplyr::across(!analyte, as.integer)) |>
    dplyr::arrange(analyte)
}

# `data` plus `timepoint` and `range`, a factor classifying every replicate
# value: OOR< / OOR> (instrument status), < LLOQ / > ULOQ (numeric value
# outside the quantification limits), in range, or empty (no value, status ok).
elisa_qc_range_class_ <- function(data) {
  lv <- c("OOR<", "< LLOQ", "in range", "> ULOQ", "OOR>", "empty")
  data |>
    dplyr::mutate(
      timepoint = qc_timepoint_(collected_at),
      range = dplyr::case_when(
        result_status %in% c("OOR<", "OOR>") ~ result_status,
        is.na(value) ~ "empty",
        !is.na(lloq) & value < lloq ~ "< LLOQ",
        !is.na(uloq) & value > uloq ~ "> ULOQ",
        .default = "in range"
      ),
      range = factor(range, levels = lv)
    )
}

# Horizontal stacked bars of counts `n` per analyte (y) filled by `fill_col`;
# facets timepoint (plus tissue if several). `colours` = named fill values.
elisa_qc_bar_plot_ <- function(counts, fill_col, colours) {
  counts$timepoint <- factor(counts$timepoint,
                             levels = qc_timepoint_levels_(counts$timepoint))
  counts$analyte <- factor(counts$analyte,
                           levels = rev(sort(unique(counts$analyte))))
  cols <- if (length(unique(counts$tissue)) > 1) {
    ggplot2::vars(tissue, timepoint)
  } else {
    ggplot2::vars(timepoint)
  }
  ggplot2::ggplot(counts, ggplot2::aes(x = n, y = analyte,
                                       fill = .data[[fill_col]])) +
    ggplot2::geom_col(position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::facet_grid(cols = cols) +
    ggplot2::scale_fill_manual(values = colours, drop = TRUE) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::labs(x = "replicate values", y = NULL, fill = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}

# Figure height (inches) for a bar plot with `n_analytes` bars per panel.
elisa_qc_bar_height_ <- function(n_analytes) {
  max(2.5, 0.3 * n_analytes + 1.5)
}

# Replicate values per analyte x tissue x timepoint: with a value, empty
# because reported out of range, or missing (empty, status ok).
elisa_qc_empty_plot_ <- function(data) {
  lv <- c("value", "empty: OOR", "empty: missing")
  counts <- data |>
    dplyr::mutate(
      timepoint = qc_timepoint_(collected_at),
      status = dplyr::case_when(
        !is.na(value) ~ "value",
        result_status %in% c("OOR<", "OOR>") ~ "empty: OOR",
        .default = "empty: missing"
      ),
      status = factor(status, levels = lv)
    ) |>
    dplyr::count(analyte, tissue, timepoint, status)
  elisa_qc_bar_plot_(counts, "status", c(
    "value" = "grey75", "empty: OOR" = "darkorange", "empty: missing" = "red3"
  ))
}

# Replicate values per range class (elisa_qc_range_class_()) per analyte x
# tissue x timepoint.
elisa_qc_range_plot_ <- function(data) {
  counts <- elisa_qc_range_class_(data) |>
    dplyr::count(analyte, tissue, timepoint, range)
  elisa_qc_bar_plot_(counts, "range", c(
    "OOR<" = "#08519c", "< LLOQ" = "#6baed6", "in range" = "grey75",
    "> ULOQ" = "#fc9272", "OOR>" = "#cb181d", "empty" = "grey92"
  ))
}

# CV (sd / mean x 100, percent) of every sample x analyte with >= 2 numeric
# replicates: mouse_ID, tissue, collected_at, analyte, n, cv, values.
elisa_qc_cv_ <- function(data) {
  data |>
    dplyr::filter(!is.na(value)) |>
    dplyr::summarise(
      n = dplyr::n(),
      cv = stats::sd(value) / mean(value) * 100,
      values = paste(format(value, trim = TRUE, digits = 4), collapse = ", "),
      .by = c(mouse_ID, tissue, collected_at, analyte)
    ) |>
    dplyr::filter(n >= 2, !is.na(cv))
}

# Samples x analyte whose replicate CV exceeds `max_cv`, highest first.
elisa_qc_replicates_ <- function(data, max_cv = 20) {
  elisa_qc_cv_(data) |>
    dplyr::filter(cv > max_cv) |>
    dplyr::select(!n) |>
    dplyr::arrange(dplyr::desc(cv))
}

# Jittered CVs (elisa_qc_cv_()) per analyte, facets timepoint (plus tissue if
# several); dashed line at `max_cv`, points above it red and labelled with
# mouse_ID. NULL without any CV.
elisa_qc_cv_plot_ <- function(cv, max_cv = 20) {
  if (nrow(cv) == 0) {
    return(NULL)
  }
  d <- cv
  d$timepoint <- qc_timepoint_(d$collected_at)
  d$timepoint <- factor(d$timepoint, levels = qc_timepoint_levels_(d$timepoint))
  d$.high <- d$cv > max_cv
  d$.label <- dplyr::if_else(d$.high, as.character(d$mouse_ID), "")
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
  ggplot2::ggplot(d, ggplot2::aes(x = analyte, y = cv)) +
    ggplot2::geom_hline(yintercept = max_cv, linetype = "dashed",
                        colour = "red3") +
    ggplot2::geom_point(ggplot2::aes(colour = .high), position = pos,
                        alpha = 0.7, size = 1.5, show.legend = FALSE) +
    label_layer +
    ggplot2::scale_colour_manual(values = c("FALSE" = "grey20",
                                            "TRUE" = "red3")) +
    ggplot2::facet_grid(cols = cols) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(x = NULL, y = "replicate CV (%)") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

# Adds a dashed LLOQ line per analyte facet (per tissue too when the plot
# splits by tissue) to the outlier plot from qc_table_plot_(). `scored` must
# carry `lloq`; missing or non-positive LLOQs (log axis) draw nothing.
elisa_qc_add_lloq_ <- function(p, scored, log_scale = FALSE) {
  if (is.null(p) || !"lloq" %in% names(scored)) {
    return(p)
  }
  by <- if (length(unique(p$data$tissue)) > 1) c("analyte", "tissue") else "analyte"
  lines <- scored |>
    dplyr::filter(!is.na(lloq), !isTRUE(log_scale) | lloq > 0) |>
    dplyr::distinct(dplyr::across(dplyr::all_of(c(by, "lloq"))))
  if (nrow(lines) == 0) {
    return(p)
  }
  lines$.row <- as.character(lines$analyte)
  p + ggplot2::geom_hline(
    data = lines, ggplot2::aes(yintercept = lloq),
    linetype = "dashed", colour = "grey30", inherit.aes = FALSE
  )
}

# Replicate means per sample x analyte; keeps every other column that is
# constant per sample (unit, LLOQ / ULOQ, joined group / batch columns).
elisa_qc_means_ <- function(data) {
  keys <- setdiff(names(data), c("replicate", "value", "result_status"))
  data |>
    dplyr::summarise(
      value = if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE),
      .by = dplyr::all_of(keys)
    ) |>
    dplyr::relocate(value, .after = analyte)
}
