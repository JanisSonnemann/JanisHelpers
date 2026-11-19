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
#' @param con A `DBI` connection from [db_connect()].
#' @returns A lazy `tbl` -- filter/collect as needed.
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
    dplyr::select(
      experiment_code, mouse_id, tissue, collected_at, assay_name,
      population_full_path, population, population_order, metric, value,
      wsp_root, condition, condition_keyword, source_file, imported_at
    )
}

#' Query ELISA measurements
#'
#' Returns a lazy `dbplyr` `tbl` joining `elisa_measurements` back through
#' `samples`/`subjects`/`experiments`/`assays`, exposing natural keys instead
#' of surrogate IDs. Call [dplyr::collect()] to materialize.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @returns A lazy `tbl` -- filter/collect as needed.
#' @export
db_query_elisa <- function(con) {
  dplyr::tbl(con, "elisa_measurements") |>
    dplyr::left_join(dplyr::tbl(con, "samples"), by = "sample_id") |>
    dplyr::left_join(dplyr::tbl(con, "subjects"), by = "subject_id") |>
    dplyr::left_join(dplyr::tbl(con, "experiments"), by = "experiment_id") |>
    dplyr::left_join(dplyr::tbl(con, "assays"), by = "assay_id") |>
    dplyr::select(
      experiment_code, mouse_id, tissue, collected_at, assay_name,
      import_unit, analyte, replicate, value, unit, result_status,
      lloq, uloq, source_file, imported_at
    )
}

#' Query histology measurements
#'
#' Returns a lazy `dbplyr` `tbl` joining `histo_measurements` back through
#' `samples`/`subjects`/`experiments`/`assays`, exposing natural keys instead
#' of surrogate IDs. Call [dplyr::collect()] to materialize.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @returns A lazy `tbl` -- filter/collect as needed.
#' @export
db_query_histo <- function(con) {
  dplyr::tbl(con, "histo_measurements") |>
    dplyr::left_join(dplyr::tbl(con, "samples"), by = "sample_id") |>
    dplyr::left_join(dplyr::tbl(con, "subjects"), by = "subject_id") |>
    dplyr::left_join(dplyr::tbl(con, "experiments"), by = "experiment_id") |>
    dplyr::left_join(dplyr::tbl(con, "assays"), by = "assay_id") |>
    dplyr::select(
      experiment_code, mouse_id, tissue, collected_at, assay_name,
      import_unit, metric, value, source_file, imported_at
    )
}

#' Query subject (mouse) meta data
#'
#' Returns a lazy `dbplyr` `tbl` with one row per mouse: `subjects` joined to
#' `experiments`, natural keys only. The stored `treatment_group` column is
#' returned as `group` (the database keeps `treatment_group` because `group`
#' is an SQL reserved word; dbplyr quotes it, so `dplyr::filter(group == ...)`
#' works on the lazy table). Call [dplyr::collect()] to materialize.
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
#'   `start_date`, `bmt_date`, `end_date`.
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
  dplyr::tbl(con, "subjects") |>
    dplyr::left_join(dplyr::tbl(con, "experiments"), by = "experiment_id") |>
    dplyr::select(
      experiment_code, mouse_id, group = treatment_group, mouse_treatment,
      sex, cage, mouse_strain, generation, dob, start_date, bmt_date, end_date
    )
}
