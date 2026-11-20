#' Export a QC-approved histology prep result to the database
#'
#' The final step of the histology prep / QC loop, run from the prep script with
#' `final <- TRUE`. Takes the prep result exactly as the last QC report showed
#' it ([histo_prep_load()]), writes every sample (exclusions are flagged on read, see
#' [db_write_exclusions()]) and replaces the import unit (the curated file's
#' name) in `histo_measurements`. Calling it is the user's "QC passed". Rows with an
#' `NA` value (an empty score) are written.
#'
#' Checks first and writes nothing if: a sample without date belongs to a
#' mouse with no `end_date`; one mouse x tissue x date x metric occurs
#' twice (after empty dates are resolved to the subject's `end_date`); the
#' source file or QC report lies outside `data_dir`. Experiment and subjects
#' must already be registered; samples and the assay are registered here.
#'
#' If the import unit is already in the database, asks before overwriting
#' ("no" changes nothing). On "yes", one transaction deletes its rows, inserts
#' the new ones and writes an `import_log` row (linked to the previous export
#' via `supersedes_import_id`). Exclusion rows matching no sample only give a
#' message.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param prep A prep result from [histo_prep_load()].
#' @param exclusions Optional tibble -- the meta workbook's `exclusions` sheet.
#'   `NULL` = nothing reported. Only used for the export message about
#'   this unit's own rows -- rules are stored by [db_write_exclusions()]. A
#'   dated row matches the stored date; end-point samples carry the mouse's
#'   end date.
#' @param data_dir The project's `data/` folder; `source_file` and
#'   `qc_report` are stored relative to it.
#' @param overwrite `NULL` (default) asks interactively when the import unit
#'   exists (errors in a non-interactive session); `TRUE` / `FALSE` answers
#'   without asking.
#' @returns Invisibly, the new `import_id`, or `NULL` if the overwrite was
#'   declined.
#' @examples
#' \dontrun{
#' prep <- histo_prep_load(path)
#' db_export_histo(con, prep, meta_read(meta_path)$exclusions)
#' }
#' @export
db_export_histo <- function(con, prep, exclusions = NULL, data_dir = "data",
                            overwrite = NULL) {
  if (is.null(prep$paths)) {
    stop("`prep` must come from histo_prep_load().", call. = FALSE)
  }
  info <- prep$info
  root <- info$root

  prepared <- prepare_sample_export_(
    con, prep$data, exclusions, info$assay,
    key = c("mouse_ID", "tissue", "collected_at", "metric"),
    key_label = "mouse_ID x tissue x date x metric"
  )
  data <- prepared$data
  matched <- prepared$matched

  source_file <- relative_to_data_dir_(info$source_path, data_dir)
  qc_report <- relative_to_data_dir_(prep$paths$qc, data_dir)

  if (!confirm_replace_(con, "histo_measurements", root, overwrite)) {
    message("Not exported -- nothing changed.")
    return(invisible(NULL))
  }

  db_write_samples(con, dplyr::distinct(
    data, mouse_id = mouse_ID, tissue, collected_at))
  db_write_assay(con, info$assay, domain = "histo")
  assay_id <- lookup_assay_id_(con, info$assay, domain = "histo")
  rows <- tibble::tibble(
    sample_id = as.integer(lookup_sample_ids_(
      con, data$mouse_ID, data$tissue, data$collected_at)),
    assay_id = as.integer(assay_id),
    import_unit = root,
    metric = as.character(data$metric),
    value = as.numeric(data$value),
    source_file = source_file
  )

  import_id <- replace_import_(
    con, "histo", "histo_measurements", rows, root, source_file, qc_report
  )
  message(
    "Exported '", root, "': ",
    nrow(dplyr::distinct(data, mouse_ID, tissue, collected_at)),
    " samples, ", n_rows_written_(con, import_id), " measurement rows, ",
    nrow(dplyr::distinct(matched$excluded_samples, mouse_ID, tissue, collected_at)),
    " sample(s) match an exclusion rule (stored, flagged on read)."
  )
  invisible(import_id)
}
