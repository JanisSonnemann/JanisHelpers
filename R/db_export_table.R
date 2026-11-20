# Internal helpers shared by db_export_elisa() and db_export_histo(): an
# import unit (root) is replaced wholesale, in one transaction.

.replace_tables <- c("elisa_measurements", "histo_measurements")

check_replace_table_ <- function(table) {
  if (!is.character(table) || length(table) != 1 || !table %in% .replace_tables) {
    stop("`table` must be one of: ", paste(.replace_tables, collapse = ", "), ".",
         call. = FALSE)
  }
  invisible(table)
}

# Internal: hint appended to unmatched-exclusion messages when a row is dated
# (rules match resolved dates, as the database does).
exclusion_date_hint_ <- function(unmatched) {
  if (!any(!is.na(unmatched$collected_at))) {
    return("")
  }
  paste0("\nEnd-point samples are stored with the mouse's end date -- a rule ",
         "dated on that day matches them.")
}

# Internal: TRUE to proceed, FALSE if the user declined. Asks only when
# `table` already holds rows of import unit `root`.
confirm_replace_ <- function(con, table, root, overwrite) {
  check_replace_table_(table)
  existing <- DBI::dbGetQuery(con, paste0(
    "SELECT DISTINCT source_file FROM ", table, " WHERE import_unit = ?"),
    params = list(root))
  if (nrow(existing) == 0) {
    return(TRUE)
  }
  if (is.null(overwrite)) {
    if (!interactive()) {
      stop(
        "Import unit '", root, "' is already in the database and this session ",
        "can't ask -- pass overwrite = TRUE or FALSE.",
        call. = FALSE
      )
    }
    overwrite <- utils::askYesNo(paste0(
      "Import unit '", root, "' is already in the database (from ",
      paste(stats::na.omit(existing$source_file), collapse = ", "),
      "). Overwrite?"
    ), default = FALSE)
  }
  isTRUE(overwrite)
}

# Internal: replace import unit `root` in `table` and log it, atomically.
# `rows` has the table's columns (minus the id and imported_at). Returns the
# new import_id.
replace_import_ <- function(con, domain, table, rows, root, source_file,
                            qc_report) {
  check_replace_table_(table)
  previous_id <- DBI::dbGetQuery(con, "
    SELECT MAX(l.import_id) AS id FROM import_log l
    JOIN domains d ON d.domain_id = l.domain_id
    WHERE d.domain_name = ? AND l.import_unit = ?",
    params = list(domain, root))$id

  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, paste0("DELETE FROM ", table, " WHERE import_unit = ?"),
                   params = list(root))
    DBI::dbAppendTable(con, table, as.data.frame(rows))

    id <- write_import_log_(
      con, domain, experiment_codes_(con, rows$sample_id), root, source_file,
      rows_written = nrow(rows), source_version = NA, qc_report = qc_report,
      qc_passed = TRUE,
      supersedes_import_id = if (is.na(previous_id)) NA_integer_ else previous_id
    )

    id
  })
}

# Internal: the validation sequence shared by db_export_elisa() and
# db_export_histo(). Resolves empty dates to the subject's end_date (what
# db_write_samples() does), matches the exclusion rules against the resolved
# dates (for messages only -- every row is written, flagged on read), and
# checks `key` is unique. Stops before anything is written.
# Returns list(data, matched).
prepare_sample_export_ <- function(con, data, exclusions, assay, key, key_label) {
  # An all-empty date column may arrive as logical NA; keep it a Date so the
  # end-date fill below stays a Date.
  data$collected_at <- as.Date(data$collected_at)

  ends <- DBI::dbGetQuery(con, "SELECT mouse_id, end_date FROM subjects")
  end_date <- ends$end_date[match(data$mouse_ID, ends$mouse_id)]
  empty <- is.na(data$collected_at)
  data$collected_at[empty] <- end_date[empty]
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

  matched <- match_sample_exclusions_(data, exclusions, assay)
  if (nrow(matched$unmatched) > 0) {
    message(
      "These exclusion rows match no sample of assay '", assay, "' -- check ",
      "them in the meta workbook (the export goes ahead).\n",
      format_rows_(dplyr::select(matched$unmatched, !n_samples_matched)),
      exclusion_date_hint_(matched$unmatched)
    )
  }

  dups <- dplyr::count(data, dplyr::pick(dplyr::all_of(key))) |>
    dplyr::filter(n > 1)
  if (nrow(dups) > 0) {
    stop(
      "These ", key_label, " occur twice ",
      "(empty dates count as the subject's end_date) -- fix the curated ",
      "file, rerun the prep script, then export. Nothing was written.\n",
      format_rows_(dups),
      call. = FALSE
    )
  }
  list(data = data, matched = matched)
}
