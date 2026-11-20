# Internal helpers for the QC checks in the facs-qc-template report
# (spec: docs/superpowers/specs/2026-11-13-facs-qc-checks-design.md).
# Called from the template as JanisHelpers:::fn().

qc_setting_names_ <- c(
  "group_col", "batch_col", "expected_mice", "expected_tissues",
  "qc_populations", "qc_mfi_populations", "qc_metrics", "min_events"
)

# Default for `min_events` (also filled into older stored settings).
qc_min_events_default_ <- 50

qc_resolve_min_events_ <- function(x) {
  if (is.null(x)) {
    return(qc_min_events_default_)
  }
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < 0) {
    stop("`min_events` must be a single non-negative number.", call. = FALSE)
  }
  x
}

# Validate a partial QC settings list against `data` and fill the defaults.
# Errors name every unknown column / population / metric.
qc_resolve_ <- function(data, qc = list()) {
  unknown <- setdiff(names(qc), qc_setting_names_)
  if (length(unknown) > 0) {
    stop("Unknown QC setting(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  for (col in c("group_col", "batch_col")) {
    v <- qc[[col]]
    if (!is.null(v) && (length(v) != 1 || !v %in% names(data))) {
      stop("`", col, "` '", paste(v, collapse = ", "),
           "' is not a column of `data`.", call. = FALSE)
    }
  }
  all_pops <- if (is.factor(data$population)) {
    levels(droplevels(data$population))
  } else {
    unique(as.character(data$population))
  }
  all_mfi <- unique(data$metric[startsWith(data$metric, "median_")])
  check_in_ <- function(value, allowed, what) {
    bad <- setdiff(value, allowed)
    if (length(bad) > 0) {
      stop("`", what, "` not found in `data`: ", paste(bad, collapse = ", "),
           call. = FALSE)
    }
  }
  check_in_(qc$qc_populations, all_pops, "qc_populations")
  check_in_(qc$qc_mfi_populations, all_pops, "qc_mfi_populations")
  check_in_(qc$qc_metrics, all_mfi, "qc_metrics")

  list(
    group_col = qc$group_col,
    batch_col = qc$batch_col,
    expected_mice = if (is.null(qc$expected_mice)) NULL else as.character(qc$expected_mice),
    expected_tissues = if (is.null(qc$expected_tissues)) {
      sort(unique(as.character(data$tissue)))
    } else {
      as.character(qc$expected_tissues)
    },
    qc_populations = if (is.null(qc$qc_populations)) all_pops else as.character(qc$qc_populations),
    qc_mfi_populations = if (is.null(qc$qc_mfi_populations)) character() else as.character(qc$qc_mfi_populations),
    qc_metrics = if (is.null(qc$qc_metrics)) all_mfi else as.character(qc$qc_metrics),
    min_events = qc_resolve_min_events_(qc$min_events)
  )
}

# Tukey whisker rule on one group's plot values `t`: "low" / "high" for points
# beyond 1.5 x IQR from the quartiles (quantile type 7, NAs ignored), else NA.
qc_whisker_side_ <- function(t) {
  q <- stats::quantile(t, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
  iqr <- q[[2]] - q[[1]]
  dplyr::case_when(
    t < q[[1]] - 1.5 * iqr ~ "low",
    t > q[[2]] + 1.5 * iqr ~ "high",
    .default = NA_character_
  )
}

# Boxplot-whisker outliers for the three plot checks. Returns the rows of
# `data` that belong to a check, plus `check`, `n_events`, `t` (the value on
# the plot's axis: fraction as is, log10(count + 1), log10(MFI); NA if not
# plotted), `excluded`, `box_median` (on the value scale), `side`
# ("low" / "high" / NA) and `flag`.
# Same rule and same points as the boxplot (Tukey: beyond 1.5 x IQR from the
# quartiles, per group x tissue x condition x population x metric), so a
# value is flagged exactly when its point lies outside the whiskers. Points
# from populations under `settings$min_events` events still feed the box but
# are never flagged.
qc_score_ <- function(data, settings, excluded_files = character()) {
  min_events <- qc_resolve_min_events_(settings$min_events)
  pops <- settings$qc_populations
  counts <- data |>
    dplyr::filter(metric == "count") |>
    dplyr::distinct(file_name, population_full_path, n_events = value)

  rows <- data |>
    dplyr::mutate(check = dplyr::case_when(
      metric == "fraction_of_parent" & population %in% pops ~ "hierarchy",
      metric == "count" & population %in% pops ~ "counts",
      metric %in% settings$qc_metrics &
        population %in% settings$qc_mfi_populations ~ "mfi",
      .default = NA_character_
    )) |>
    dplyr::filter(!is.na(check)) |>
    dplyr::left_join(counts, by = c("file_name", "population_full_path")) |>
    dplyr::mutate(
      excluded = file_name %in% excluded_files,
      t = dplyr::case_when(
        check == "hierarchy" ~ value,
        check == "counts" ~ log10(pmax(value, 0) + 1),
        check == "mfi" & value > 0 ~ log10(pmax(value, 1e-12)),
        .default = NA_real_
      )
    )

  group_col <- settings$group_col
  if (is.null(group_col)) {
    return(dplyr::mutate(rows, box_median = NA_real_, side = NA_character_,
                         flag = FALSE))
  }

  rows |>
    dplyr::group_by(dplyr::across(dplyr::all_of(
      c(group_col, "tissue", "condition", "population_full_path", "metric")
    ))) |>
    dplyr::mutate(
      side = qc_whisker_side_(t),
      box_median = stats::median(value[!is.na(t)])
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      side = dplyr::if_else(!is.na(n_events) & n_events >= min_events,
                            side, NA_character_),
      flag = !is.na(side)
    )
}

# Display label per row: the leaf population name, or "parent/leaf" when the
# leaf name occurs under more than one parent (facs_read_wsp() keeps such
# gates as one factor level). Ordered by population level, then path.
qc_pop_label_ <- function(population, full_path) {
  leaf <- as.character(population)
  paths_per_leaf <- tapply(full_path, leaf, function(x) length(unique(x)))
  dup <- leaf %in% names(paths_per_leaf)[paths_per_leaf > 1]
  parent <- basename(dirname(full_path))
  label <- ifelse(dup, paste0(parent, "/", leaf), leaf)
  lev_order <- if (is.factor(population)) as.integer(population) else match(leaf, unique(leaf))
  factor(label, levels = unique(label[order(lev_order, full_path)]))
}

# Batch values as display labels: dates (Date / POSIXct) as "YYYY-MM-DD".
qc_batch_label_ <- function(x) {
  if (inherits(x, c("Date", "POSIXt"))) format(x, "%Y-%m-%d") else as.character(x)
}

# Expected mouse x tissue x condition vs. present. NULL without expected mice.
qc_completeness_ <- function(data, settings) {
  if (is.null(settings$expected_mice)) {
    return(NULL)
  }
  expected <- tidyr::expand_grid(
    mouse_ID = settings$expected_mice,
    tissue = settings$expected_tissues,
    condition = sort(unique(data$condition))
  )
  present <- dplyr::distinct(data, mouse_ID, tissue, condition)
  keys <- c("mouse_ID", "tissue", "condition")
  list(
    n_expected = nrow(expected),
    n_present = nrow(dplyr::semi_join(expected, present, by = keys)),
    missing = expected |>
      dplyr::anti_join(present, by = keys) |>
      dplyr::arrange(mouse_ID, tissue, condition)
  )
}

# One row per flagged value; mice ranked by the number of distinct
# tissue x condition x check combinations flagged.
qc_flag_table_ <- function(scored, settings) {
  label_of <- function(col) {
    if (is.null(col)) NA_character_ else qc_batch_label_(scored[[col]])
  }
  scored |>
    dplyr::mutate(group = label_of(settings$group_col),
                  batch = label_of(settings$batch_col),
                  population = qc_pop_label_(population, population_full_path)) |>
    dplyr::filter(flag) |>
    dplyr::mutate(
      score = dplyr::n_distinct(tissue, condition, check),
      .by = mouse_ID
    ) |>
    dplyr::arrange(dplyr::desc(score), mouse_ID, tissue, condition, check,
                   population) |>
    dplyr::transmute(
      mouse_ID, score, tissue, condition, group, batch,
      population = as.character(population), metric, value, box_median, side,
      excluded
    )
}

# Figure height (inches) for a plot with `n_rows` facet rows.
qc_fig_height_ <- function(n_rows) {
  max(3, 1.6 * n_rows + 1)
}

qc_has_ggrepel_ <- function() {
  requireNamespace("ggrepel", quietly = TRUE)
}

# Point shape by exclusion: filled circle, or a cross for samples on the
# exclusions sheet. Only the cross gets a legend key, and only when the plot
# has an excluded point (`.excluded` = the plot's column).
qc_excluded_shape_ <- function(.excluded) {
  shown <- any(.excluded == "TRUE")
  ggplot2::scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 4), limits = c("FALSE", "TRUE"),
    breaks = if (shown) "TRUE" else NULL,
    labels = if (shown) "on exclusions sheet" else NULL,
    name = "excluded"
  )
}

# `.excluded` ("TRUE"/"FALSE") for the shape aesthetic from a logical vector
# (NULL, e.g. no `excluded` column, = none excluded).
qc_excluded_col_ <- function(excluded, n) {
  if (is.null(excluded)) excluded <- rep(FALSE, n)
  as.character(excluded %in% TRUE)
}

# Rows of one check (and condition / MFI population) prepared for plotting:
# `.row` = facet row (population, or marker for MFI), `.group` = group label,
# `.y` = plotted value (counts + 1 and positive MFIs only, for the log axis).
qc_plot_rows_ <- function(scored, settings, check, condition, population) {
  d <- scored[scored$check == check, ]
  if (!is.null(condition)) d <- d[d$condition == condition, ]
  if (check == "mfi") {
    if (is.null(population)) population <- settings$qc_mfi_populations[1]
    d <- d[!is.na(d$value) & d$value > 0 &
             as.character(d$population) %in% population, ]
  }
  if (nrow(d) == 0) {
    return(NULL)
  }
  group_col <- settings$group_col
  pop_label <- qc_pop_label_(d$population, d$population_full_path)
  d$.row <- if (check != "mfi") {
    pop_label
  } else if (length(unique(d$population_full_path)) > 1) {
    paste0(pop_label, ": ", sub("^median_", "", d$metric))
  } else {
    sub("^median_", "", d$metric)
  }
  d$.group <- if (is.null(group_col)) "all" else qc_batch_label_(d[[group_col]])
  d$.y <- if (check == "counts") d$value + 1 else d$value
  d$.excluded <- qc_excluded_col_(
    if ("excluded" %in% names(d)) d$excluded, nrow(d)
  )
  d
}

qc_y_scale_ <- function(check) {
  if (check == "hierarchy") {
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"))
  } else {
    ggplot2::scale_y_log10()
  }
}

# Individuals plot: x = group, boxplot + jittered points, mouse_ID labels on
# flagged values; facets population (MFI: marker) x tissue. NULL if no rows.
qc_plot_individuals_ <- function(scored, settings, check, condition = NULL,
                                 population = NULL) {
  d <- qc_plot_rows_(scored, settings, check, condition, population)
  if (is.null(d)) {
    return(NULL)
  }
  # Label every row ("" if not flagged) so the labels share the points' jitter.
  d$.label <- dplyr::if_else(d$flag, d$mouse_ID, "")
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
  ggplot2::ggplot(d, ggplot2::aes(x = .group, y = .y, fill = .group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    ggplot2::geom_point(ggplot2::aes(shape = .excluded), position = pos,
                        alpha = 0.7, size = 1.5) +
    label_layer +
    qc_excluded_shape_(d$.excluded) +
    qc_y_scale_(check) +
    ggplot2::facet_grid(rows = ggplot2::vars(.row), cols = ggplot2::vars(tissue),
                        scales = "free_y") +
    ggplot2::labs(x = NULL, y = NULL, fill = settings$group_col) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))
}

# Batch plot (layout C) for one tissue: x = batch, coloured by batch; facets
# population (MFI: marker) x group. NULL without batch column or rows.
qc_plot_batch_ <- function(scored, settings, check, tissue, condition = NULL,
                           population = NULL) {
  batch_col <- settings$batch_col
  if (is.null(batch_col)) {
    return(NULL)
  }
  d <- qc_plot_rows_(scored, settings, check, condition, population)
  if (is.null(d)) {
    return(NULL)
  }
  d <- d[d$tissue == tissue, ]
  if (nrow(d) == 0) {
    return(NULL)
  }
  d$.batch <- factor(qc_batch_label_(d[[batch_col]]))
  pos <- ggplot2::position_jitter(width = 0.15, height = 0, seed = 1)
  ggplot2::ggplot(d, ggplot2::aes(x = .batch, y = .y, fill = .batch)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    ggplot2::geom_point(ggplot2::aes(colour = .batch, shape = .excluded),
                        position = pos, size = 1.5) +
    qc_excluded_shape_(d$.excluded) +
    qc_y_scale_(check) +
    ggplot2::facet_grid(rows = ggplot2::vars(.row), cols = ggplot2::vars(.group),
                        scales = "free_y") +
    ggplot2::labs(x = NULL, y = NULL, fill = batch_col, colour = batch_col) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

# QC settings of a prep result: the stored `info$qc` (`min_events` filled in
# for preps saved before it was a setting), or defaults for a prep saved
# before QC settings existed (no group / batch column, no expected mice).
qc_settings_ <- function(prep) {
  qc <- prep$info$qc
  if (!is.null(qc)) {
    if (is.null(qc$min_events)) qc$min_events <- qc_min_events_default_
    return(qc)
  }
  qc_resolve_(prep$data)
}

# -- empty values / missing stats / keywords -------------------------------------

# Values present per metric x tissue x condition: one row per metric x every
# tissue x condition in `data` (a metric missing there is "absent"). Columns:
# metric, tissue, condition, n_values, n_empty (NA values), n_no_events (of
# those, in populations with 0 events -- expected, e.g. an MFI of no cells),
# n_files_empty (files with unexpected empties), status ("ok", "partly
# empty", "all empty", "absent"; only unexpected empties count).
qc_empty_values_ <- function(data) {
  counts <- data |>
    dplyr::filter(metric == "count") |>
    dplyr::distinct(file_name, population_full_path, n_events = value)
  present <- data |>
    dplyr::left_join(counts, by = c("file_name", "population_full_path")) |>
    dplyr::mutate(.empty = is.na(value),
                  .no_events = .empty & dplyr::coalesce(n_events == 0, FALSE)) |>
    dplyr::summarise(
      n_values = dplyr::n(),
      n_empty = sum(.empty),
      n_no_events = sum(.no_events),
      n_files_empty = dplyr::n_distinct(file_name[.empty & !.no_events]),
      .by = c(metric, tissue, condition)
    )
  tidyr::expand_grid(
    metric = unique(data$metric),
    dplyr::distinct(data, tissue, condition)
  ) |>
    dplyr::left_join(present, by = c("metric", "tissue", "condition")) |>
    dplyr::mutate(
      status = dplyr::case_when(
        is.na(n_values) ~ "absent",
        n_empty == n_no_events ~ "ok",
        n_empty == n_values ~ "all empty",
        .default = "partly empty"
      ),
      dplyr::across(c(n_values, n_empty, n_no_events, n_files_empty),
                    ~ dplyr::coalesce(.x, 0L))
    ) |>
    dplyr::arrange(metric, tissue, condition)
}

# Wide display of qc_empty_values_(): one row per metric, one column per
# tissue (x condition, if there are several); cells "ok", "absent",
# "all empty" or "<n>/<N> empty", with " (<k> at 0 events)" appended when
# some empties are populations without events.
qc_empty_grid_ <- function(empty) {
  multi_cond <- length(unique(empty$condition)) > 1
  empty |>
    dplyr::mutate(
      column = if (multi_cond) paste0(tissue, " / ", condition) else tissue,
      cell = dplyr::if_else(status == "partly empty",
                            paste0(n_empty - n_no_events, "/", n_values, " empty"),
                            status),
      cell = dplyr::if_else(n_no_events > 0,
                            paste0(cell, " (", n_no_events, " at 0 events)"), cell)
    ) |>
    dplyr::select(metric, column, cell) |>
    tidyr::pivot_wider(names_from = column, values_from = cell)
}

# Population x metric rows that only some files of a tissue x condition have
# (a gate or statistic missing for some samples in FlowJo). One row per
# tissue x condition x gate x set of files, metrics collapsed: tissue,
# condition, population (leaf name, "parent/leaf" if ambiguous), metrics
# (comma-separated), n_files ("<present>/<total>"), files ("only in: ..."
# when present in at most half the files, else "missing in: ..."; mouse IDs).
qc_missing_stats_ <- function(data) {
  keys <- c("tissue", "condition")
  files <- dplyr::distinct(data, file_name, mouse_ID, tissue, condition)
  n_files <- dplyr::count(files, tissue, condition, name = "n_total")
  out <- data |>
    dplyr::distinct(tissue, condition, population, population_full_path,
                    metric, file_name) |>
    dplyr::summarise(present = list(file_name),
                     .by = c(tissue, condition, population,
                             population_full_path, metric)) |>
    dplyr::left_join(n_files, by = keys) |>
    dplyr::mutate(n_present = lengths(present)) |>
    dplyr::filter(n_present < n_total)
  if (nrow(out) == 0) {
    return(tibble::tibble(tissue = character(), condition = character(),
                          population = character(), metrics = character(),
                          n_files = character(), files = character()))
  }
  out$population <- as.character(
    qc_pop_label_(out$population, out$population_full_path)
  )
  out$files <- purrr::pmap_chr(
    list(out$present, out$tissue, out$condition, out$n_present, out$n_total),
    function(present, tis, cond, n_present, n_total) {
      f <- files[files$tissue == tis & files$condition == cond, ]
      only <- n_present <= n_total / 2
      ids <- sort(f$mouse_ID[(f$file_name %in% present) == only])
      paste0(if (only) "only in: " else "missing in: ", paste(ids, collapse = ", "))
    }
  )
  out |>
    dplyr::mutate(n_files = paste0(n_present, "/", n_total)) |>
    dplyr::summarise(metrics = paste(metric, collapse = ", "),
                     .by = c(tissue, condition, population,
                             population_full_path, n_files, files)) |>
    dplyr::arrange(tissue, condition, population_full_path) |>
    dplyr::select(tissue, condition, population, metrics, n_files, files)
}

# Files with an empty group / batch value (columns joined from the meta
# workbook). Columns: file_name, mouse_ID, tissue, condition, empty (the
# empty column's name). Zero rows if none, or if neither column is set.
qc_empty_keywords_ <- function(data, settings) {
  cols <- intersect(c(settings$group_col, settings$batch_col), names(data))
  if (length(cols) == 0) {
    return(tibble::tibble(file_name = character(), mouse_ID = character(),
                          tissue = character(), condition = character(),
                          empty = character()))
  }
  data |>
    dplyr::distinct(file_name, mouse_ID, tissue, condition,
                    dplyr::across(dplyr::all_of(cols))) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(cols), as.character)) |>
    tidyr::pivot_longer(dplyr::all_of(cols), names_to = "empty") |>
    dplyr::filter(is.na(value)) |>
    dplyr::select(file_name, mouse_ID, tissue, condition, empty) |>
    dplyr::arrange(mouse_ID, tissue, condition, empty)
}

# Files without a mouse_ID or tissue keyword -- dropped by facs_prep_save()
# (compensation / FMO files are expected here). Columns: file_name,
# mouse_ID, tissue.
qc_dropped_files_ <- function(data) {
  data |>
    dplyr::filter(is.na(mouse_ID) | is.na(tissue)) |>
    dplyr::distinct(file_name, mouse_ID, tissue) |>
    dplyr::arrange(file_name)
}
