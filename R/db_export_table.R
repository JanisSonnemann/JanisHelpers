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

# Internal: exclusion rows naming a panel that is neither `assay` nor any
# registered assay (any domain) raise a warning -- a mistyped panel would be
# silently ignored. Not an error: on a fresh database the other assays are
# not registered yet, and their rows must not block this export.
check_exclusion_panels_ <- function(con, exclusions, assay) {
  if (is.null(exclusions) || !"panel" %in% names(exclusions)) {
    return(invisible(TRUE))
  }
  known <- unique(c(assay, DBI::dbGetQuery(con, "SELECT assay_name FROM assays")$assay_name))
  panels <- unique(stats::na.omit(as.character(exclusions$panel)))
  unknown <- setdiff(panels[panels != ""], known)
  if (length(unknown) > 0) {
    warning(
      "The exclusions sheet names panel(s) not registered in the database: ",
      paste0("'", unknown, "'", collapse = ", "), ". Known panels / assays: ",
      paste(known, collapse = ", "), ". Their rows are ignored here -- check ",
      "for typos if that assay is already exported.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Internal: hint appended to unmatched-exclusion messages when a row is dated
# (exclusions match the sheet's dates, where end-point samples are empty).
exclusion_date_hint_ <- function(unmatched) {
  if (!any(!is.na(unmatched$collected_at))) {
    return("")
  }
  paste0("\nEnd-point samples have an empty date in the sheet -- leave ",
         "collected_at empty to target them.")
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
# `rows` has the table's columns (minus the id and imported_at); `applicable`
# is the applied exclusion rows from match_sample_exclusions_() (or NULL).
# Returns the new import_id.
replace_import_ <- function(con, domain, table, rows, root, source_file,
                            qc_report, applicable) {
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

    if (!is.null(applicable) && nrow(applicable) > 0) {
      ex <- data.frame(
        import_id = rep(id, nrow(applicable)),
        mouse_id = applicable$mouse_ID,
        tissue = applicable$tissue,
        panel = applicable$panel,
        collected_at = applicable$collected_at,
        reason = applicable$reason,
        decided_on = applicable$date,
        n_files_matched = applicable$n_samples_matched
      )
      DBI::dbAppendTable(con, "import_exclusions", ex)
    }
    id
  })
}

# Internal: the validation sequence shared by db_export_elisa() and
# db_export_histo(). Matches exclusions, drops the excluded samples (anti-join
# on mouse_ID x tissue x collected_at), resolves empty dates to the subject's
# end_date and checks `key` is unique. Stops before anything is written.
# Returns list(data, excluded, matched).
prepare_sample_export_ <- function(con, data, exclusions, assay, key, key_label) {
  # An all-empty date column may arrive as logical NA; keep it a Date so the
  # end-date fill below stays a Date.
  data$collected_at <- as.Date(data$collected_at)
  check_exclusion_panels_(con, exclusions, assay)
  matched <- match_sample_exclusions_(data, exclusions, assay)
  if (nrow(matched$unmatched) > 0) {
    stop(
      "These exclusion rows match no sample of assay '", assay, "' -- ",
      "fix them in the meta workbook. Nothing was written.\n",
      format_rows_(dplyr::select(matched$unmatched, !n_samples_matched)),
      exclusion_date_hint_(matched$unmatched),
      call. = FALSE
    )
  }
  excluded <- dplyr::distinct(matched$excluded_samples, mouse_ID, tissue, collected_at)
  data <- dplyr::anti_join(data, excluded,
                           by = c("mouse_ID", "tissue", "collected_at"))
  if (nrow(data) == 0) {
    stop("Every sample is excluded -- nothing left to export.", call. = FALSE)
  }

  # Resolve empty dates to the subject's end date (what db_write_samples() does).
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
  list(data = data, excluded = excluded, matched = matched)
}
