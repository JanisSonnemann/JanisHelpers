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

# Samples x analyte with >= 2 numeric replicates whose CV (sd / mean x 100,
# percent) exceeds `max_cv`.
elisa_qc_replicates_ <- function(data, max_cv = 20) {
  data |>
    dplyr::filter(!is.na(value)) |>
    dplyr::summarise(
      n = dplyr::n(),
      cv = stats::sd(value) / mean(value) * 100,
      values = paste(format(value, trim = TRUE, digits = 4), collapse = ", "),
      .by = c(mouse_ID, tissue, collected_at, analyte)
    ) |>
    dplyr::filter(n >= 2, !is.na(cv), cv > max_cv) |>
    dplyr::select(!n) |>
    dplyr::arrange(dplyr::desc(cv))
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
