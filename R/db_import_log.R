# Internal: the distinct experiment code(s) of the given sample ids, sorted
# and joined with ", " (an ELISA plate can span experiments).
experiment_codes_ <- function(con, sample_ids) {
  keys <- tibble::tibble(sample_id = unique(as.integer(sample_ids)))
  duckdb::duckdb_register(con, "tmp_log_sample_ids", keys)
  on.exit(duckdb::duckdb_unregister(con, "tmp_log_sample_ids"), add = TRUE)
  codes <- DBI::dbGetQuery(con, "
    SELECT DISTINCT e.experiment_code FROM tmp_log_sample_ids k
    JOIN samples sm ON sm.sample_id = k.sample_id
    JOIN subjects sub ON sub.subject_id = sm.subject_id
    JOIN experiments e ON e.experiment_id = sub.experiment_id
    ORDER BY e.experiment_code
  ")$experiment_code
  paste(codes, collapse = ", ")
}

# Internal: write one row to import_log and return its import_id.
#
# Called by every importer (db_write_elisa(), db_write_histo(), and later the
# Stage 3 FACS export) INSIDE the caller's DBI::dbWithTransaction(), so the
# data and its log row land together or not at all. Opens no transaction of
# its own (DuckDB does not nest transactions). Log rows are never updated or
# deleted by package code.
#
# domain: "facs", "elisa" or "histo" (resolved to domains.domain_id).
# See docs/superpowers/specs/2026-11-12-db-schema-import-log-design.md for
# what each column holds per domain.
write_import_log_ <- function(con, domain, experiment_code, import_unit,
                              source_file, rows_written,
                              source_version = NA, qc_report = NA_character_,
                              qc_passed = NA, supersedes_import_id = NA_integer_) {
  domain_id <- DBI::dbGetQuery(
    con, "SELECT domain_id FROM domains WHERE domain_name = ?",
    params = list(domain)
  )$domain_id
  if (length(domain_id) != 1) {
    stop("Unknown domain '", domain, "'.", call. = FALSE)
  }
  # Stands in for a self-referencing foreign key, which DuckDB's
  # EXPORT DATABASE can't handle.
  if (!is.na(supersedes_import_id)) {
    found <- DBI::dbGetQuery(
      con, "SELECT COUNT(*) AS n FROM import_log WHERE import_id = ?",
      params = list(as.integer(supersedes_import_id))
    )$n
    if (found == 0) {
      stop(
        "supersedes_import_id ", supersedes_import_id,
        " doesn't match any import_log row.",
        call. = FALSE
      )
    }
  }

  DBI::dbGetQuery(con, "
    INSERT INTO import_log (
      domain_id, experiment_code, import_unit, source_file, source_version,
      qc_report, qc_passed, rows_written, package_version, supersedes_import_id
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    RETURNING import_id
  ", params = list(
    domain_id,
    as.character(experiment_code),
    as.character(import_unit),
    as.character(source_file),
    as.Date(source_version),
    as.character(qc_report),
    as.logical(qc_passed),
    as.integer(rows_written),
    as.character(utils::packageVersion("JanisHelpers")),
    as.integer(supersedes_import_id)
  ))$import_id[[1]]
}
