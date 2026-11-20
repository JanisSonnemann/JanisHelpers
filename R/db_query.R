#' @importFrom dbplyr sql
NULL

#' Query FACS measurements
#'
#' Returns a lazy `dbplyr` `tbl` joining `facs_measurements` back through
#' `facs_stains`/`samples`/`subjects`/`experiments`/`assays`, exposing
#' natural keys (`mouse_id`, `tissue`, `experiment_code`, `assay_name`)
#' instead of surrogate IDs. Call [dplyr::collect()] to materialize.
#'
#' `population_order` is the population's position in its panel's display
#' order (the prep's [facs_select_populations()] order); `NA` for imports
#' made before the column existed. Positions are only comparable within one
#' `wsp_root`. After collecting, [facs_order_populations()] turns
#' `population` back into an ordered factor. On a database file created
#' before the column existed and opened read-only (so [db_connect()] could not
#' add it), `population_order` is returned as all `NA`.
#'
#' Every stored row is returned. `excluded` (logical, `FALSE` when no rule
#' matches) and `exclusion_reason` (the matching rules' reasons, sorted and
#' joined with `"; "`) come from the experiment's current exclusion rule set
#' (see [db_write_exclusions()]); filter with `dplyr::filter(!excluded)`.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @returns A lazy `tbl` -- filter/collect as needed. Ends with `excluded` and
#'   `exclusion_reason`.
#' @export
db_query_facs <- function(con) {
  measurements <- dplyr::tbl(con, "facs_measurements")
  if (!"population_order" %in% DBI::dbListFields(con, "facs_measurements")) {
    measurements <- measurements |>
      dplyr::mutate(population_order = dbplyr::sql("CAST(NULL AS INTEGER)"))
  }
  measurements |>
    dplyr::left_join(dplyr::tbl(con, "facs_stains"), by = "facs_stain_id") |>
    dplyr::left_join(dplyr::tbl(con, "samples"), by = "sample_id") |>
    dplyr::left_join(dplyr::tbl(con, "subjects"), by = "subject_id") |>
    dplyr::left_join(dplyr::tbl(con, "experiments"), by = "experiment_id") |>
    dplyr::left_join(dplyr::tbl(con, "assays"), by = "assay_id") |>
    add_exclusion_flags_(con, "facs") |>
    dplyr::select(
      experiment_code, mouse_id, tissue, collected_at, assay_name,
      population_full_path, population, population_order, metric, value,
      wsp_root, condition, condition_keyword, source_file, imported_at,
      excluded, exclusion_reason
    )
}

#' Query ELISA measurements
#'
#' Returns a lazy `dbplyr` `tbl` joining `elisa_measurements` back through
#' `samples`/`subjects`/`experiments`/`assays`, exposing natural keys instead
#' of surrogate IDs. Call [dplyr::collect()] to materialize.
#'
#' Every stored row is returned, flagged by the experiment's current exclusion
#' rule set as in [db_query_facs()]. Rules naming a `condition` or `file_name`
#' never match here.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @returns A lazy `tbl` -- filter/collect as needed. Ends with `excluded` and
#'   `exclusion_reason`.
#' @export
db_query_elisa <- function(con) {
  dplyr::tbl(con, "elisa_measurements") |>
    dplyr::left_join(dplyr::tbl(con, "samples"), by = "sample_id") |>
    dplyr::left_join(dplyr::tbl(con, "subjects"), by = "subject_id") |>
    dplyr::left_join(dplyr::tbl(con, "experiments"), by = "experiment_id") |>
    dplyr::left_join(dplyr::tbl(con, "assays"), by = "assay_id") |>
    add_exclusion_flags_(con, "elisa") |>
    dplyr::select(
      experiment_code, mouse_id, tissue, collected_at, assay_name,
      import_unit, analyte, replicate, value, unit, result_status,
      lloq, uloq, source_file, imported_at, excluded, exclusion_reason
    )
}

#' Query histology measurements
#'
#' Returns a lazy `dbplyr` `tbl` joining `histo_measurements` back through
#' `samples`/`subjects`/`experiments`/`assays`, exposing natural keys instead
#' of surrogate IDs. Call [dplyr::collect()] to materialize.
#'
#' Every stored row is returned, flagged by the experiment's current exclusion
#' rule set as in [db_query_facs()]. Rules naming a `condition` or `file_name`
#' never match here.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @returns A lazy `tbl` -- filter/collect as needed. Ends with `excluded` and
#'   `exclusion_reason`.
#' @export
db_query_histo <- function(con) {
  dplyr::tbl(con, "histo_measurements") |>
    dplyr::left_join(dplyr::tbl(con, "samples"), by = "sample_id") |>
    dplyr::left_join(dplyr::tbl(con, "subjects"), by = "subject_id") |>
    dplyr::left_join(dplyr::tbl(con, "experiments"), by = "experiment_id") |>
    dplyr::left_join(dplyr::tbl(con, "assays"), by = "assay_id") |>
    add_exclusion_flags_(con, "histo") |>
    dplyr::select(
      experiment_code, mouse_id, tissue, collected_at, assay_name,
      import_unit, metric, value, source_file, imported_at, excluded,
      exclusion_reason
    )
}

# Internal: left-join a lazy measurement table (still carrying
# `<domain>_measurement_id`) to its exclusion-flag view: `excluded` (FALSE
# when no current rule matches) and `exclusion_reason` (NA then). A file
# opened read-only before the views existed gets FALSE / NA throughout.
add_exclusion_flags_ <- function(tbl, con, domain) {
  view <- paste0(domain, "_exclusion_flags")
  if (!DBI::dbExistsTable(con, view)) {
    return(dplyr::mutate(
      tbl,
      excluded = FALSE,
      exclusion_reason = dbplyr::sql("CAST(NULL AS VARCHAR)")
    ))
  }
  tbl |>
    dplyr::left_join(dplyr::tbl(con, view),
                     by = paste0(domain, "_measurement_id")) |>
    dplyr::mutate(excluded = dplyr::coalesce(excluded, FALSE))
}

#' Query subject (mouse) meta data
#'
#' Returns a lazy `dbplyr` `tbl` with one row per mouse: `subjects` joined to
#' `experiments`, natural keys only. The stored `treatment_group` column is
#' returned as `group` (the database keeps `treatment_group` because `group`
#' is an SQL reserved word; dbplyr quotes it, so `dplyr::filter(group == ...)`
#' works on the lazy table). Call [dplyr::collect()] to materialize.
#'
#' On a database file created before the dropout columns or the exclusion
#' tables existed and opened read-only (so [db_connect()] could not add
#' them), `mdc_id` / `early_reason` are `NA`, `early` is `FALSE`, `id_source`
#' is `"mouse_ID"` and no mouse is excluded.
#'
#' Intended use in an analysis document: connect read-only, query the subjects
#' once, filter to the experiment, collect, and set the group factor levels
#' once. Then join only the columns a dataset needs with the checked join
#' [meta_annotate()], by `experiment_code` and `mouse_id` (mouse IDs are only
#' unique within an experiment). The measurement queries ([db_query_facs()],
#' [db_query_elisa()], [db_query_histo()]) carry no subject meta.
#'
#' @param con A `DBI` connection -- [db_connect()], or for analysis a
#'   read-only `DBI::dbConnect(duckdb::duckdb(), path, read_only = TRUE)`.
#' @returns A lazy `tbl` with columns `experiment_code`, `mouse_id`, `group`,
#'   `mouse_treatment`, `sex`, `cage`, `mouse_strain`, `generation`, `dob`,
#'   `start_date`, `bmt_date`, `end_date`, `mdc_id` (facility ID), `early`
#'   (logical, left the experiment early), `early_reason`, `id_source`
#'   (`"mouse_ID"` or `"facility"`: where `mouse_id` came from),
#'   `excluded_mouse_wide` (logical: the experiment's current exclusion rule
#'   set holds a rule for the mouse with every scope column empty) and
#'   `exclusion_reason` (the reasons of those mouse-wide rules, `NA` if none).
#' @examples
#' \dontrun{
#' con <- DBI::dbConnect(duckdb::duckdb(), "~/databases/lab.duckdb",
#'                       read_only = TRUE)
#' subjects <- db_query_subjects(con) |>
#'   dplyr::filter(experiment_code == "26-1") |>
#'   dplyr::collect() |>
#'   dplyr::mutate(group = factor(group, levels = c("MPO409", "irrel")))
#' facs <- db_query_facs(con) |>
#'   dplyr::filter(experiment_code == "26-1") |>
#'   dplyr::collect() |>
#'   meta_annotate(dplyr::select(subjects, experiment_code, mouse_id, group),
#'                 by = c("experiment_code", "mouse_id")) |>
#'   facs_order_populations()
#' DBI::dbDisconnect(con, shutdown = TRUE)
#' }
#' @export
db_query_subjects <- function(con) {
  subjects <- dplyr::tbl(con, "subjects")
  # a file created before these columns existed and opened read-only (so
  # db_connect() could not add them): typed NULLs, then the defaults below
  added <- c(mdc_id = "VARCHAR", early = "BOOLEAN", early_reason = "VARCHAR",
             id_source = "VARCHAR")
  missing <- setdiff(base::names(added), DBI::dbListFields(con, "subjects"))
  if (length(missing) > 0) {
    casts <- lapply(missing, function(col) {
      dbplyr::sql(paste0("CAST(NULL AS ", added[[col]], ")"))
    })
    subjects <- dplyr::mutate(subjects, !!!stats::setNames(casts, missing))
  }
  if (DBI::dbExistsTable(con, "exclusion_rules") &&
      DBI::dbExistsTable(con, "exclusion_sets")) {
    # mouse-wide = a rule of the experiment's current set with no scope
    mouse_wide <- dplyr::tbl(con, dbplyr::sql("
      SELECT cur.experiment_id, r.mouse_id, TRUE AS excluded_mouse_wide,
             string_agg(DISTINCT r.reason, '; ' ORDER BY r.reason) AS exclusion_reason
      FROM (SELECT experiment_id, MAX(exclusion_set_id) AS exclusion_set_id
            FROM exclusion_sets GROUP BY experiment_id) cur
      JOIN exclusion_rules r ON r.exclusion_set_id = cur.exclusion_set_id
      WHERE r.tissue IS NULL AND r.panel IS NULL AND r.condition IS NULL
        AND r.file_name IS NULL AND r.collected_at IS NULL
      GROUP BY cur.experiment_id, r.mouse_id
    "))
    subjects <- subjects |>
      dplyr::left_join(mouse_wide, by = c("experiment_id", "mouse_id"))
  } else {
    subjects <- dplyr::mutate(
      subjects,
      excluded_mouse_wide = FALSE,
      exclusion_reason = dbplyr::sql("CAST(NULL AS VARCHAR)")
    )
  }
  subjects |>
    dplyr::left_join(dplyr::tbl(con, "experiments"), by = "experiment_id") |>
    dplyr::select(
      experiment_code, mouse_id, group = treatment_group, mouse_treatment,
      sex, cage, mouse_strain, generation, dob, start_date, bmt_date, end_date,
      mdc_id, early, early_reason, id_source, excluded_mouse_wide,
      exclusion_reason
    ) |>
    dplyr::mutate(
      early = dplyr::coalesce(early, FALSE),
      id_source = dplyr::coalesce(id_source, "mouse_ID"),
      excluded_mouse_wide = dplyr::coalesce(excluded_mouse_wide, FALSE)
    )
}
