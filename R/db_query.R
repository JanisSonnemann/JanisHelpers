#' Query FACS measurements
#'
#' Returns a lazy `dbplyr` `tbl` joining `facs_measurements` back through
#' `facs_stains`/`samples`/`subjects`/`experiments`/`assays`, exposing
#' natural keys (`mouse_id`, `tissue`, `experiment_code`, `assay_name`)
#' instead of surrogate IDs. Call [dplyr::collect()] to materialize.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @returns A lazy `tbl` -- filter/collect as needed.
#' @export
db_query_facs <- function(con) {
  dplyr::tbl(con, "facs_measurements") |>
    dplyr::left_join(dplyr::tbl(con, "facs_stains"), by = "facs_stain_id") |>
    dplyr::left_join(dplyr::tbl(con, "samples"), by = "sample_id") |>
    dplyr::left_join(dplyr::tbl(con, "subjects"), by = "subject_id") |>
    dplyr::left_join(dplyr::tbl(con, "experiments"), by = "experiment_id") |>
    dplyr::left_join(dplyr::tbl(con, "assays"), by = "assay_id") |>
    dplyr::select(
      experiment_code, mouse_id, tissue, assay_name,
      population_full_path, population, metric, value,
      count_method, vol_stained, vol_resuspended, vol_measured,
      bead_volume_added, bead_concentration, bead_population_path,
      source_file, imported_at
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
      experiment_code, mouse_id, tissue, assay_name,
      cytokine, sample_id_raw, replicate, value, unit, result_status,
      source_file, imported_at
    )
}

#' Query histology measurements
#'
#' Returns a lazy `dbplyr` `tbl` joining `histo_measurements` back through
#' `samples`/`subjects`/`experiments`/`assays`, exposing natural keys instead
#' of surrogate IDs. `assay_name` is `NA` for rows written without an
#' `assay_name`. Call [dplyr::collect()] to materialize.
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
      experiment_code, mouse_id, tissue, assay_name,
      metric, value, source_file, imported_at
    )
}
