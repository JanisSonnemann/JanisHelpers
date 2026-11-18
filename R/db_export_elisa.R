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

  check_exclusion_panels_(con, exclusions, info$assay)
  matched <- match_sample_exclusions_(prep$data, exclusions, info$assay)
  if (nrow(matched$unmatched) > 0) {
    stop(
      "These exclusion rows match no sample of assay '", info$assay, "' -- ",
      "fix them in the meta workbook. Nothing was written.\n",
      format_rows_(dplyr::select(matched$unmatched, !n_samples_matched)),
      call. = FALSE
    )
  }
  excluded <- dplyr::distinct(matched$excluded_samples, mouse_ID, tissue, collected_at)
  data <- dplyr::anti_join(prep$data, excluded,
                           by = c("mouse_ID", "tissue", "collected_at"))
  if (nrow(data) == 0) {
    stop("Every sample is excluded -- nothing left to export.", call. = FALSE)
  }

  # Resolve empty dates to the subject's end date (what db_write_samples() does).
  ends <- DBI::dbGetQuery(con, "SELECT mouse_id, end_date FROM subjects")
  end_date <- ends$end_date[match(data$mouse_ID, ends$mouse_id)]
  data$collected_at <- as.Date(ifelse(is.na(data$collected_at),
                                      end_date, data$collected_at))
  if (anyNA(data$collected_at)) {
    stop(
      "No date and no subject end_date (or unknown subject) for mouse_ID(s): ",
      paste(sprintf("'%s'", unique(data$mouse_ID[is.na(data$collected_at)])),
            collapse = ", "),
      " -- register the subject with an end_date or give the sheet a date. ",
      "Nothing was written.",
      call. = FALSE
    )
  }

  key <- c("mouse_ID", "tissue", "collected_at", "analyte", "replicate")
  dups <- dplyr::count(data, dplyr::pick(dplyr::all_of(key))) |>
    dplyr::filter(n > 1)
  if (nrow(dups) > 0) {
    stop(
      "These mouse_ID x tissue x date x analyte x replicate occur twice ",
      "(empty dates count as the subject's end_date) -- fix the curated ",
      "file, rerun the prep script, then export. Nothing was written.\n",
      format_rows_(dups),
      call. = FALSE
    )
  }

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
