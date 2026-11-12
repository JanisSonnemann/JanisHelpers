# Internal: `path` relative to `data_dir`, erroring if it lies outside.
relative_to_data_dir_ <- function(path, data_dir) {
  root <- normalizePath(data_dir, mustWork = TRUE)
  full <- normalizePath(path, mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  if (!startsWith(full, prefix)) {
    stop("'", path, "' is not inside data_dir '", data_dir, "'.", call. = FALSE)
  }
  substring(full, nchar(prefix) + 1)
}

# Internal: format a tibble for an error/warning message.
format_rows_ <- function(x) {
  paste(utils::capture.output(print(x, n = Inf, width = Inf)), collapse = "\n")
}

#' Export a QC-approved FACS prep result to the database
#'
#' The final step of the prep / QC loop, run from the panel's prep script with
#' `final <- TRUE`. Takes the prep result exactly as the last QC report showed
#' it ([facs_prep_load()]), drops the files hit by the exclusions sheet, and
#' writes the rest -- raw and calculated metrics -- as one workspace root.
#' Calling it is the user's "QC passed": the import is logged with
#' `qc_passed = TRUE`.
#'
#' Checks first and writes nothing if: several files share one mouse x tissue
#' x condition; an exclusion row for this panel matches no file; no file is
#' left after exclusions; the workspace or QC report lies outside `data_dir`.
#' Warns if a newer dated version of the workspace exists beside it.
#'
#' If the root is already in the database, asks before overwriting ("no"
#' changes nothing). On "yes", one transaction deletes the root's stains and
#' measurements, inserts the new ones, writes an `import_log` row (linked to
#' the previous export of the root via `supersedes_import_id`) and stores the
#' applied exclusion rows in `import_exclusions`. Experiment and subjects must
#' already be registered; samples and the assay are registered here.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param prep A prep result from [facs_prep_load()].
#' @param exclusions Optional tibble -- the meta workbook's `exclusions` sheet
#'   (see [facs_match_exclusions()]). `NULL` = nothing excluded.
#' @param data_dir The project's `data/` folder; `source_file` and
#'   `qc_report` are stored relative to it.
#' @param overwrite `NULL` (default) asks interactively when the root already
#'   exists; `TRUE` / `FALSE` answers without asking.
#' @returns Invisibly, the new `import_id`, or `NULL` if the overwrite was
#'   declined.
#' @examples
#' \dontrun{
#' prep <- facs_prep_load(wsp_path, panel_dir)
#' db_export_facs(con, prep, meta_read(meta_path)$exclusions)
#' }
#' @export
db_export_facs <- function(con, prep, exclusions = NULL, data_dir = "data",
                           overwrite = NULL) {
  if (is.null(prep$paths)) {
    stop("`prep` must come from facs_prep_load().", call. = FALSE)
  }
  info <- prep$info
  root <- info$wsp_root

  clashes <- facs_check_file_clashes(prep$data)
  if (nrow(clashes) > 0) {
    stop(
      "Several FCS files share one mouse_ID x tissue x condition -- fix the ",
      "keywords in FlowJo, rerun the prep script, then export. Nothing was ",
      "written.\n", format_rows_(clashes),
      call. = FALSE
    )
  }

  matched <- facs_match_exclusions(prep$data, exclusions, info$assay)
  if (nrow(matched$unmatched) > 0) {
    stop(
      "These exclusion rows match no file of panel '", info$assay, "' -- fix ",
      "them in the meta workbook. Nothing was written.\n",
      format_rows_(dplyr::select(matched$unmatched, !n_files_matched)),
      call. = FALSE
    )
  }
  data <- dplyr::filter(prep$data, !file_name %in% matched$excluded_files$file_name)
  if (nrow(data) == 0) {
    stop("Every file is excluded -- nothing left to export.", call. = FALSE)
  }

  versions <- list.files(dirname(info$wsp_path), pattern = "\\.wsp$")
  versions <- versions[startsWith(versions, root) &
    grepl("^_[0-9]{8}\\.wsp$", substring(versions, nchar(root) + 1))]
  version_dates <- as.Date(sub(".*_([0-9]{8})\\.wsp$", "\\1", versions), "%Y%m%d")
  newer <- versions[version_dates > info$wsp_version]
  if (length(newer) > 0) {
    warning(
      "A newer version of root '", root, "' exists: ",
      paste(newer, collapse = ", "), ". Exporting ",
      basename(info$wsp_path), " anyway.",
      call. = FALSE
    )
  }

  source_file <- relative_to_data_dir_(info$wsp_path, data_dir)
  qc_report <- relative_to_data_dir_(prep$paths$qc, data_dir)

  existing <- DBI::dbGetQuery(con, "
    SELECT DISTINCT m.source_file FROM facs_stains fs
    LEFT JOIN facs_measurements m ON m.facs_stain_id = fs.facs_stain_id
    WHERE fs.wsp_root = ?", params = list(root))
  if (nrow(existing) > 0) {
    if (is.null(overwrite)) {
      overwrite <- utils::askYesNo(paste0(
        "Root '", root, "' is already in the database (from ",
        paste(stats::na.omit(existing$source_file), collapse = ", "),
        "). Overwrite?"
      ))
    }
    if (!isTRUE(overwrite)) {
      message("Not exported -- nothing changed.")
      return(invisible(NULL))
    }
  }
  previous_id <- DBI::dbGetQuery(con, "
    SELECT MAX(l.import_id) AS id FROM import_log l
    JOIN domains d ON d.domain_id = l.domain_id
    WHERE d.domain_name = 'facs' AND l.import_unit = ?",
    params = list(root))$id

  db_write_samples(con, dplyr::distinct(data, mouse_id = mouse_ID, tissue))
  db_write_assay(con, info$assay, domain = "facs")
  assay_id <- lookup_assay_id_(con, info$assay, domain = "facs")
  stains <- dplyr::distinct(data, file_name, mouse_ID, tissue, condition)
  stains$sample_id <- purrr::map2_int(
    stains$mouse_ID, stains$tissue,
    function(m, t) as.integer(lookup_sample_id_(con, m, t))
  )
  kw <- if (is.null(info$condition_keyword)) NA_character_ else info$condition_keyword

  import_id <- DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "
      DELETE FROM facs_measurements WHERE facs_stain_id IN
        (SELECT facs_stain_id FROM facs_stains WHERE wsp_root = ?)",
      params = list(root))
    DBI::dbExecute(con, "DELETE FROM facs_stains WHERE wsp_root = ?",
                   params = list(root))

    n <- 0
    for (i in seq_len(nrow(stains))) {
      s <- stains[i, ]
      check_facs_root_(con, s$sample_id, assay_id, root, source_file)
      check_facs_root_keyword_(con, root, kw)
      n <- n + insert_facs_stain_(
        con, s$sample_id, assay_id,
        dplyr::filter(data, file_name == s$file_name),
        root, source_file, s$condition, kw
      )
    }

    id <- write_import_log_(
      con, "facs", experiment_codes_(con, stains$sample_id), root,
      source_file, rows_written = n, source_version = info$wsp_version,
      qc_report = qc_report, qc_passed = TRUE,
      supersedes_import_id = if (is.na(previous_id)) NA_integer_ else previous_id
    )

    ex <- matched$applicable
    if (nrow(ex) > 0) {
      DBI::dbExecute(con, "
        INSERT INTO import_exclusions (
          import_id, mouse_id, tissue, panel, condition, file_name, reason,
          decided_on, n_files_matched
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        params = list(
          rep(id, nrow(ex)), ex$mouse_ID, ex$tissue, ex$panel, ex$condition,
          ex$file_name, ex$reason, ex$date, ex$n_files_matched
        ))
    }
    id
  })

  message(
    "Exported root '", root, "': ", nrow(stains), " stains, ", n_rows_written_(con, import_id),
    " measurement rows, ", nrow(matched$excluded_files), " file(s) excluded."
  )
  invisible(import_id)
}

# Internal: rows_written of one import_log row (for the summary message).
n_rows_written_ <- function(con, import_id) {
  DBI::dbGetQuery(con, "SELECT rows_written FROM import_log WHERE import_id = ?",
                  params = list(import_id))$rows_written
}
