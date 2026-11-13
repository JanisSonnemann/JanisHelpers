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

# QC settings of a prep result: the stored `info$qc`, or defaults for a prep
# saved before QC settings existed (no group / batch column, no expected mice).
qc_settings_ <- function(prep) {
  if (!is.null(prep$info$qc)) {
    return(prep$info$qc)
  }
  qc_resolve_(prep$data)
}
