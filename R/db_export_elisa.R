#' Export a QC-approved ELISA prep result to the database
#'
#' The final step of the ELISA prep / QC loop, run from the prep script with
#' `final <- TRUE`. Takes the prep result exactly as the last QC report showed
#' it ([elisa_prep_load()]), drops the samples hit by the exclusions sheet and
#' replaces the import unit (the curated file's name) in `elisa_measurements`.
#' Calling it is the user's "QC passed".
#'
#' Checks first and writes nothing if: an exclusion row for this assay matches
#' no sample; no sample is left; a sample without date belongs to a mouse with
#' no `end_date`; one mouse x tissue x date x analyte x replicate occurs twice
#' (after empty dates are resolved to the subject's `end_date`); the source
#' file or QC report lies outside `data_dir`. Experiment and subjects must
#' already be registered; samples and the assay are registered here.
#'
#' If the import unit is already in the database, asks before overwriting
#' ("no" changes nothing). On "yes", one transaction deletes its rows, inserts
#' the new ones, writes an `import_log` row (linked to the previous export via
#' `supersedes_import_id`) and stores the applied exclusion rows in
#' `import_exclusions`.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param prep A prep result from [elisa_prep_load()].
#' @param exclusions Optional tibble -- the meta workbook's `exclusions` sheet.
#'   `NULL` = nothing excluded.
#' @param data_dir The project's `data/` folder; `source_file` and
#'   `qc_report` are stored relative to it.
#' @param overwrite `NULL` (default) asks interactively when the import unit
#'   exists (errors in a non-interactive session); `TRUE` / `FALSE` answers
#'   without asking.
#' @returns Invisibly, the new `import_id`, or `NULL` if the overwrite was
#'   declined.
#' @examples
#' \dontrun{
#' prep <- elisa_prep_load(path)
#' db_export_elisa(con, prep, meta_read(meta_path)$exclusions)
#' }
#' @export
db_export_elisa <- function(con, prep, exclusions = NULL, data_dir = "data",
                            overwrite = NULL) {
  if (is.null(prep$paths)) {
    stop("`prep` must come from elisa_prep_load().", call. = FALSE)
  }
  info <- prep$info
  root <- info$root

  prepared <- prepare_sample_export_(
    con, prep$data, exclusions, info$assay,
    key = c("mouse_ID", "tissue", "collected_at", "analyte", "replicate"),
    key_label = "mouse_ID x tissue x date x analyte x replicate"
  )
  data <- prepared$data
  excluded <- prepared$excluded
  matched <- prepared$matched

  source_file <- relative_to_data_dir_(info$source_path, data_dir)
  qc_report <- relative_to_data_dir_(prep$paths$qc, data_dir)

  if (!confirm_replace_(con, "elisa_measurements", root, overwrite)) {
    message("Not exported -- nothing changed.")
    return(invisible(NULL))
  }

  db_write_samples(con, dplyr::distinct(
    data, mouse_id = mouse_ID, tissue, collected_at))
  db_write_assay(con, info$assay, domain = "elisa")
  assay_id <- lookup_assay_id_(con, info$assay, domain = "elisa")
  rows <- tibble::tibble(
    sample_id = as.integer(lookup_sample_ids_(
      con, data$mouse_ID, data$tissue, data$collected_at)),
    assay_id = as.integer(assay_id),
    import_unit = root,
    analyte = as.character(data$analyte),
    replicate = as.integer(data$replicate),
    value = as.numeric(data$value),
    unit = as.character(data$unit),
    result_status = as.character(data$result_status),
    lloq = as.numeric(data$lloq),
    uloq = as.numeric(data$uloq),
    source_file = source_file
  )

  import_id <- replace_import_(
    con, "elisa", "elisa_measurements", rows, root, source_file, qc_report,
    matched$applicable
  )
  message(
    "Exported '", root, "': ", nrow(dplyr::distinct(data, mouse_ID, tissue, collected_at)),
    " samples, ", n_rows_written_(con, import_id), " measurement rows, ",
    nrow(excluded), " sample(s) excluded."
  )
  invisible(import_id)
}
