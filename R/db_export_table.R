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
# registered assay (any domain) stop the export -- a mistyped panel would be
# silently ignored.
check_exclusion_panels_ <- function(con, exclusions, assay) {
  if (is.null(exclusions) || !"panel" %in% names(exclusions)) {
    return(invisible(TRUE))
  }
  known <- unique(c(assay, DBI::dbGetQuery(con, "SELECT assay_name FROM assays")$assay_name))
  panels <- unique(stats::na.omit(as.character(exclusions$panel)))
  unknown <- setdiff(panels[panels != ""], known)
  if (length(unknown) > 0) {
    stop(
      "The exclusions sheet names unknown panel(s): ",
      paste0("'", unknown, "'", collapse = ", "), ". Known panels / assays: ",
      paste(known, collapse = ", "), ". Nothing was written.",
      call. = FALSE
    )
  }
  invisible(TRUE)
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
