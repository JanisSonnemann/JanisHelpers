#' Load one experiment from the lab database as a named list
#'
#' Collects everything stored for one experiment into a single list, so an
#' analysis document keeps one object in its workspace: subject meta, the
#' cohort, the current exclusion rules, and one tibble per import unit (FACS panel =
#' `wsp_root`, ELISA / histology file = `import_unit`). Each measurement tibble
#' gets the requested subject columns joined with [meta_annotate()]; FACS
#' tibbles get their panel's population display order via
#' [facs_order_populations()]. Nothing is written to the database, so a
#' read-only connection is the intended use.
#'
#' Element names: the import unit without a trailing `_long` and without the
#' experiment code (plus its `-`/`_` separator) at the end, e.g.
#' `overview-26-1` -> `overview`, `anti-mpo-titer_26-1_long` ->
#' `anti-mpo-titer`. Override any of them with `names`.
#'
#' @param con A `DBI` connection, e.g.
#'   `DBI::dbConnect(duckdb::duckdb(), path, read_only = TRUE)`.
#' @param experiment Experiment code, e.g. `"26-1"`.
#' @param names Optional named character vector `c(new_name = "import_unit")`
#'   renaming those import units' list elements.
#' @param group_levels Optional character vector: `group` becomes a factor
#'   with these levels in `meta` and every joined tibble. A stored group not
#'   listed is an error; `NA` groups stay `NA`.
#' @param subject_cols Columns of `meta` joined onto each measurement tibble
#'   (by `experiment_code` + `mouse_id`). Default `"group"`.
#' @param keep_excluded If `FALSE` (default), rows flagged `excluded` are
#'   dropped from every import unit, mice excluded mouse-wide and mice with a
#'   facility ID (`id_source == "facility"`) are dropped from `meta` (with a
#'   message for the excluded ones), and the flag columns are removed. `TRUE`
#'   keeps all rows and subjects, with `excluded` / `exclusion_reason` in the
#'   units and `id_source` / `excluded_mouse_wide` / `exclusion_reason` in
#'   `meta`.
#' @returns A named list: `meta` (analysed mice, columns as
#'   [db_query_subjects()]), `cohort` (every subject of the experiment, with
#'   all [db_query_subjects()] columns), `exclusions` (the current exclusion
#'   rules: `exclusion_set_id`, `mouse_id`, `tissue`, `panel`, `condition`,
#'   `file_name`, `collected_at`, `reason`, `decided_on`; a message says when
#'   none were recorded), then one collected tibble per import unit -- FACS
#'   panels, then ELISA, then histology, each alphabetical by import unit.
#' @examples
#' \dontrun{
#' dat <- local({
#'   con <- DBI::dbConnect(duckdb::duckdb(), "~/databases/lab.duckdb",
#'                         read_only = TRUE)
#'   on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
#'   db_load_experiment(con, "26-1", group_levels = c("MPO409", "irrel"),
#'                      names = c(titer = "anti-mpo-titer_26-1_long"))
#' })
#' }
#' @export
db_load_experiment <- function(con, experiment, names = NULL, group_levels = NULL,
                               subject_cols = "group", keep_excluded = FALSE) {
  codes <- DBI::dbGetQuery(con, "SELECT experiment_code FROM experiments")$experiment_code
  if (length(experiment) != 1 || !experiment %in% codes) {
    stop(
      "Experiment '", paste(experiment, collapse = ", "), "' not in the database. ",
      "Available: ", paste(sort(codes), collapse = ", "), ".",
      call. = FALSE
    )
  }

  cohort <- db_query_subjects(con) |>
    dplyr::filter(.data$experiment_code == !!experiment) |>
    dplyr::collect()
  if (!is.null(group_levels)) {
    unknown <- setdiff(stats::na.omit(unique(cohort$group)), group_levels)
    if (length(unknown) > 0) {
      stop(
        "Group(s) not in `group_levels`: ", paste(unknown, collapse = ", "), ".",
        call. = FALSE
      )
    }
    cohort <- dplyr::mutate(cohort, group = factor(.data$group, levels = group_levels))
  }
  subject_cols <- setdiff(subject_cols, c("experiment_code", "mouse_id"))
  missing_cols <- setdiff(subject_cols, base::names(cohort))
  if (length(missing_cols) > 0) {
    stop(
      "`subject_cols` not in the subject meta: ", paste(missing_cols, collapse = ", "),
      ". Available: ", paste(base::names(cohort), collapse = ", "), ".",
      call. = FALSE
    )
  }

  domains <- list(
    facs = db_query_facs(con) |>
      dplyr::filter(.data$experiment_code == !!experiment) |>
      dplyr::collect() |>
      dplyr::mutate(unit_ = .data$wsp_root),
    elisa = db_query_elisa(con) |>
      dplyr::filter(.data$experiment_code == !!experiment) |>
      dplyr::collect() |>
      dplyr::mutate(unit_ = .data$import_unit),
    histo = db_query_histo(con) |>
      dplyr::filter(.data$experiment_code == !!experiment) |>
      dplyr::collect() |>
      dplyr::mutate(unit_ = .data$import_unit)
  )
  # alphabetical regardless of case or locale ("overview" before "Treg")
  sort_ci_ <- function(x) x[order(tolower(x), x, method = "radix")]
  units <- unlist(lapply(domains, function(d) sort_ci_(unique(d$unit_))), use.names = FALSE)
  unit_domain <- rep(base::names(domains), vapply(domains, function(d) length(unique(d$unit_)), integer(1)))

  el_names <- resolve_unit_names_(units, experiment, names)

  exclusions <- current_exclusion_rules_(con, experiment)
  if (!DBI::dbGetQuery(con, "
    SELECT COUNT(*) AS n FROM exclusion_sets s JOIN experiments e USING (experiment_id)
    WHERE e.experiment_code = ?", params = list(experiment))$n > 0) {
    message("No exclusions recorded for experiment ", experiment, ".")
  }

  meta <- cohort
  if (!keep_excluded) {
    meta <- dplyr::filter(meta, .data$id_source != "facility")
    dropped <- sort(cohort$mouse_id[cohort$excluded_mouse_wide])
    if (length(dropped) > 0) {
      message(
        "Dropping mouse(s) excluded from every import unit: ",
        paste(dropped, collapse = ", "), " (keep_excluded = TRUE keeps them)."
      )
    }
    meta <- meta |>
      dplyr::filter(!.data$excluded_mouse_wide) |>
      dplyr::select(!c("id_source", "excluded_mouse_wide", "exclusion_reason"))
  }

  join_cols <- dplyr::select(meta, "experiment_code", "mouse_id", dplyr::all_of(subject_cols))
  tibbles <- lapply(seq_along(units), function(i) {
    d <- domains[[unit_domain[i]]]
    d <- dplyr::select(d[d$unit_ == units[i], ], !"unit_")
    if (!keep_excluded) {
      d <- d |>
        dplyr::filter(!.data$excluded, .data$mouse_id %in% meta$mouse_id) |>
        dplyr::select(!c("excluded", "exclusion_reason"))
    }
    d <- meta_annotate(d, join_cols, by = c("experiment_code", "mouse_id"))
    if (unit_domain[i] == "facs") d <- facs_order_populations(d)
    d
  })
  base::names(tibbles) <- el_names

  c(list(meta = meta, cohort = cohort, exclusions = exclusions), tibbles)
}

# Internal: list element names for import units -- automatic names, then the
# user's `names` overrides; errors on unknown units and bad names.
resolve_unit_names_ <- function(units, experiment, overrides) {
  out <- unit_names_(units, experiment)
  if (!is.null(overrides)) {
    unknown <- setdiff(overrides, units)
    if (length(unknown) > 0) {
      stop(
        "Import unit(s) in `names` not in experiment ", experiment, ": ",
        paste(unknown, collapse = ", "), ". Available: ",
        paste(units, collapse = ", "), ".",
        call. = FALSE
      )
    }
    new <- base::names(overrides)
    if (is.null(new)) new <- rep("", length(overrides))
    out[match(overrides, units)] <- new
  }
  bad <- unique(c(
    out[is.na(out) | out == ""],
    out[duplicated(out)],
    intersect(out, c("meta", "cohort", "exclusions"))
  ))
  if (length(bad) > 0 || any(is.na(out) | out == "")) {
    stop(
      "Element name(s) empty, duplicated or reserved (meta, cohort, exclusions): ",
      paste(ifelse(is.na(bad) | bad == "", "<empty name>", bad), collapse = ", "),
      ". Set them with `names = c(new_name = \"import_unit\")`.",
      call. = FALSE
    )
  }
  out
}

# Internal: strip a trailing "_long", then a trailing "-<experiment>" or
# "_<experiment>" (fixed string, no regex -- codes contain "-").
unit_names_ <- function(units, experiment) {
  strip_end_ <- function(x, suffix) {
    hit <- endsWith(x, suffix)
    x[hit] <- substr(x[hit], 1, nchar(x[hit]) - nchar(suffix))
    x
  }
  out <- strip_end_(units, "_long")
  for (sep in c("-", "_")) out <- strip_end_(out, paste0(sep, experiment))
  out
}
