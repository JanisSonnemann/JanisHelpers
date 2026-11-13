# Internal helpers for the QC checks in the facs-qc-template report
# (spec: docs/superpowers/specs/2026-11-13-facs-qc-checks-design.md).
# Called from the template as JanisHelpers:::fn().

qc_setting_names_ <- c(
  "group_col", "batch_col", "expected_mice", "expected_tissues",
  "qc_populations", "qc_mfi_populations", "qc_metrics"
)

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
    qc_metrics = if (is.null(qc$qc_metrics)) all_mfi else as.character(qc$qc_metrics)
  )
}

# Robust-z outlier scoring for the three plot checks. Returns the rows of
# `data` that belong to a check, plus `check`, `n_events`, `t` (transformed
# value, NA if unscorable), `excluded`, `z` (NA when not scored) and `flag`.
# Scored within group x tissue x condition x population x metric; excluded
# files and unscorable values don't feed the box median / MAD.
qc_score_ <- function(data, settings, excluded_files = character(),
                      min_n = 4, min_events = 50, z_cut = 3.5,
                      mad_floor = c(logit = 0.2, log10 = 0.1)) {
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
        check == "hierarchy" & value > 0 ~
          stats::qlogis((n_events + 0.5) / (n_events / value + 1)),
        check == "counts" ~ log10(value + 1),
        check == "mfi" & value > 0 ~ log10(value),
        .default = NA_real_
      ),
      t = dplyr::if_else(is.na(n_events) | n_events < min_events, NA_real_, t),
      floor = dplyr::if_else(check == "hierarchy", mad_floor[["logit"]],
                             mad_floor[["log10"]])
    )

  group_col <- settings$group_col
  if (is.null(group_col)) {
    return(rows |>
      dplyr::mutate(z = NA_real_, flag = FALSE) |>
      dplyr::select(!floor))
  }

  rows |>
    dplyr::mutate(t = dplyr::if_else(is.na(.data[[group_col]]), NA_real_, t)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(
      c(group_col, "tissue", "condition", "population_full_path", "metric")
    ))) |>
    dplyr::mutate(
      .ref = !is.na(t) & !excluded,
      .n = sum(.ref),
      .med = stats::median(t[.ref]),
      .mad = 1.4826 * stats::median(abs(t[.ref] - .med)),
      z = dplyr::if_else(.n >= min_n & !is.na(t),
                         (t - .med) / pmax(.mad, floor), NA_real_)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(flag = !is.na(z) & abs(z) > z_cut) |>
    dplyr::select(!c(floor, .ref, .n, .med, .mad))
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
                  batch = label_of(settings$batch_col)) |>
    dplyr::filter(flag) |>
    dplyr::mutate(
      score = dplyr::n_distinct(tissue, condition, check),
      .by = mouse_ID
    ) |>
    dplyr::arrange(dplyr::desc(score), mouse_ID, dplyr::desc(abs(z))) |>
    dplyr::transmute(
      mouse_ID, score, tissue, condition, group, batch,
      population = as.character(population), metric, value, z, excluded
    )
}

# QC settings of a prep result: the stored `info$qc`, or defaults for a prep
# saved before QC settings existed (no group / batch column, no expected mice).
qc_settings_ <- function(prep) {
  if (!is.null(prep$info$qc)) {
    return(prep$info$qc)
  }
  qc_resolve_(prep$data)
}
