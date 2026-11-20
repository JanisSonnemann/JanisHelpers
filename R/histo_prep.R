# Internal: paths of the prep result and QC report of a curated histology file.
histo_prep_paths_ <- function(path) {
  root <- tools::file_path_sans_ext(basename(path))
  list(
    root = root,
    prep = file.path(dirname(path), "processed", paste0(root, "_prep.rds")),
    qc = file.path(dirname(path), "qc", paste0(root, "_qc.html"))
  )
}

.histo_prep_cols <- c("mouse_ID", "tissue", "collected_at", "metric", "value")

#' Save a histology prep result for QC and export
#'
#' First step of the histology prep / QC loop. Stores the curated data (all
#' samples; exclusions are only flagged in the QC report) with the QC settings
#' at `<folder of path>/processed/<file name>_prep.rds`. The file name without
#' extension is the import unit ("root") stored in the database. Warns when
#' one mouse x tissue x date x metric occurs twice (after resolving empty
#' dates with `end_dates`).
#'
#' @param data Tibble from [histo_read_curated()] (optionally with ratio
#'   metrics from [histo_calc_ratio()] and group / batch columns joined on) --
#'   needs `mouse_ID`, `tissue`, `collected_at`, `metric`, `value`. An
#'   optional `comment` column is listed in the QC report and not exported.
#' @param path Path of the curated `.xlsx` the data came from.
#' @param assay Assay name the data is exported under.
#' @param group_col Optional column of `data` used to group the outlier plots.
#' @param batch_col Optional batch column of `data`, shown in the QC report
#'   overview.
#' @param expected_mice Character vector of mice expected in the data.
#'   Default: the mice in `data`.
#' @param rules Character vector of plausibility rules for the QC report.
#'   Default `character()` (none).
#' @param end_dates Optional end (death) dates used by the QC report to place
#'   end-point samples (empty `collected_at`) on the time axis: a Date vector
#'   named by `mouse_ID`, or a data frame with `mouse_ID` and one date column.
#'   Default `NULL` (end-point samples are shown as "end"). Stored as
#'   `info$end_dates`; the export resolves dates from the database instead.
#' @returns Invisibly, the path of the written `.rds` (which holds a list with
#'   `data` and `info`: `source_path`, `root`, `assay`, `created_at`,
#'   `package_version`, `end_dates`, `qc`).
#' @examples
#' \dontrun{
#' histo_read_curated("data/26-1_exp/raw/histo/kidney_26-1.xlsx", "kidney") |>
#'   histo_prep_save("data/26-1_exp/raw/histo/kidney_26-1.xlsx", "PAS")
#' }
#' @export
histo_prep_save <- function(data, path, assay, group_col = NULL,
                            batch_col = NULL, expected_mice = NULL,
                            rules = character(), end_dates = NULL) {
  missing_cols <- setdiff(.histo_prep_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("`data` is missing column(s): ", paste(missing_cols, collapse = ", "),
         ".", call. = FALSE)
  }
  for (col in c("group_col", "batch_col")) {
    v <- get(col)
    if (!is.null(v) && (length(v) != 1 || !v %in% names(data))) {
      stop("`", col, "` '", paste(v, collapse = ", "),
           "' is not a column of `data`.", call. = FALSE)
    }
  }
  if (!is.character(rules)) {
    stop("`rules` must be a character vector.", call. = FALSE)
  }
  end_dates <- qc_end_dates_(end_dates)
  dups <- qc_duplicate_keys_(
    qc_resolve_end_dates_(data, end_dates),
    c("mouse_ID", "tissue", "collected_at", "metric")
  )
  if (nrow(dups) > 0) {
    warning(
      "These mouse_ID x tissue x date x metric occur twice (empty dates count ",
      "as the end date, where known) -- the export will refuse them:\n",
      format_rows_(dups), call. = FALSE
    )
  }
  paths <- histo_prep_paths_(path)
  prep <- list(
    data = data,
    info = list(
      source_path = path,
      root = paths$root,
      assay = assay,
      created_at = Sys.time(),
      package_version = as.character(utils::packageVersion("JanisHelpers")),
      end_dates = end_dates,
      qc = list(
        group_col = group_col,
        batch_col = batch_col,
        expected_mice = if (is.null(expected_mice)) {
          sort(unique(as.character(data$mouse_ID)))
        } else {
          as.character(expected_mice)
        },
        rules = rules
      )
    )
  )
  dir.create(dirname(paths$prep), recursive = TRUE, showWarnings = FALSE)
  saveRDS(prep, paths$prep)
  message("Prep result saved to ", paths$prep)
  invisible(paths$prep)
}

#' Load the saved histology prep result for the final export
#'
#' Reads the `.rds` written by [histo_prep_save()] and checks it is current:
#' the QC report must exist and be newer than the prep result, and the curated
#' `.xlsx` must not have been changed since.
#'
#' @param path Path of the curated `.xlsx`, as passed to [histo_prep_save()].
#' @returns The prep list (`data`, `info`) plus `paths` (`prep`, `qc`).
#' @examples
#' \dontrun{
#' prep <- histo_prep_load("data/26-1_exp/raw/histo/kidney_26-1.xlsx")
#' }
#' @export
histo_prep_load <- function(path) {
  paths <- histo_prep_paths_(path)
  for (p in c(paths$prep, paths$qc)) {
    if (!file.exists(p)) {
      stop("'", p, "' not found. Run the prep script with final <- FALSE first.",
           call. = FALSE)
    }
  }
  if (file.mtime(paths$qc) < file.mtime(paths$prep)) {
    stop(
      "The QC report is older than the prep result (did rendering fail?) -- ",
      "rerun the prep script with final <- FALSE and check the report.",
      call. = FALSE
    )
  }
  if (file.exists(path) && file.mtime(path) > file.mtime(paths$prep)) {
    stop(
      "'", basename(path), "' is newer than its prep result -- rerun the ",
      "prep script with final <- FALSE and check the QC report.",
      call. = FALSE
    )
  }
  prep <- readRDS(paths$prep)
  prep$paths <- paths[c("prep", "qc")]
  prep
}

#' Run the histology prep / QC loop for one curated file
#'
#' One call for the whole loop, driven by a `final` switch (see
#' `inst/scripts/histo_prep_template.R`). `final = FALSE` reads the curated
#' sheet ([histo_read_curated()]), applies `prepare()`, joins the group /
#' batch columns from `mice`, saves the prep result ([histo_prep_save()]) and
#' renders the QC report ([report_knit_histo_qc()]). `final = TRUE` loads the
#' saved prep ([histo_prep_load()]), exports it ([db_export_histo()]) and
#' optionally snapshots the database.
#'
#' @param final `FALSE` to prep and render the QC report, `TRUE` to export.
#' @param path Path of the curated `.xlsx`.
#' @param assay Assay name the data is exported under.
#' @param tissue Tissue of the file, passed to [histo_read_curated()].
#' @param prepare Function applied to the long data tibble before saving.
#' @param mice Optional tibble with `mouse_ID` and the `group_col` /
#'   `batch_col` columns, joined into the data with [meta_annotate()]. If it
#'   has a `death_date` column, that becomes the prep's `end_dates` (unless
#'   `end_dates` is passed in `...`).
#' @param exclusions Optional tibble -- the meta workbook's `exclusions`
#'   sheet. With `final = TRUE` it is stored as the experiment's
#'   rule set ([db_write_exclusions()]); `NULL` leaves the stored rule set
#'   unchanged (with a message), a zero-row sheet clears it.
#' @param con A `DBI` connection from [db_connect()], required when
#'   `final = TRUE`; left open.
#' @param group_col,batch_col Character scalars or `NULL`, see
#'   [histo_prep_save()].
#' @param ... Further QC settings passed to [histo_prep_save()]
#'   (`expected_mice`, `rules`, `end_dates`). Used only when
#'   `final = FALSE`.
#' @param meta_path Path of the meta workbook, or `NULL`. With `final = TRUE`
#'   the `exclusions` sheet is also written to the database as the
#'   experiment's rule set ([db_write_exclusions()]), with this path
#'   (relative to `data_dir`) as its source.
#' @param data_dir The project's `data/` folder, see [db_export_histo()].
#' @param snapshot_root Folder for dated snapshots after an export. `NULL`
#'   (default) = none.
#' @param overwrite Passed to [db_export_histo()].
#' @returns Invisibly: with `final = FALSE` the QC report path; with
#'   `final = TRUE` the new `import_id`, or `NULL` if the overwrite was
#'   declined.
#' @examples
#' \dontrun{
#' histo_prep_run(
#'   FALSE, "data/26-1_exp/raw/histo/kidney_26-1.xlsx", "PAS",
#'   tissue = "kidney", mice = meta$meta, group_col = "group"
#' )
#' }
#' @export
histo_prep_run <- function(final, path, assay, tissue = NA_character_,
                           prepare = identity, mice = NULL, exclusions = NULL,
                           con = NULL, group_col = NULL, batch_col = NULL,
                           ..., meta_path = NULL, data_dir = "data", snapshot_root = NULL,
                           overwrite = NULL) {
  if (!is.logical(final) || length(final) != 1L || is.na(final)) {
    stop("`final` must be TRUE or FALSE.", call. = FALSE)
  }
  if (final) {
    if (is.null(con)) {
      stop("`con` is required for the final export (final = TRUE).",
           call. = FALSE)
    }
    prep <- histo_prep_load(path)
    write_prep_exclusions_(con, prep, exclusions, meta_path, data_dir)
    id <- db_export_histo(con, prep, exclusions, data_dir = data_dir,
                          overwrite = overwrite)
    if (!is.null(id) && !is.null(snapshot_root)) db_snapshot(con, snapshot_root)
    return(invisible(id))
  }

  if (!file.exists(path)) {
    stop("Curated file not found: ", path, call. = FALSE)
  }
  data <- histo_read_curated(path, tissue) |>
    prepare()
  cols <- c(group_col, batch_col)
  if (!is.null(mice) && length(cols) > 0L) {
    data <- meta_annotate(
      data,
      mice |>
        dplyr::filter(!is.na(mouse_ID)) |>
        dplyr::select(mouse_ID, dplyr::all_of(unique(cols))),
      by = "mouse_ID"
    )
  }
  args <- list(data, path, assay, group_col = group_col,
               batch_col = batch_col, ...)
  if (is.null(args$end_dates) && !is.null(mice) && "death_date" %in% names(mice)) {
    args$end_dates <- mice |>
      dplyr::filter(!is.na(mouse_ID)) |>
      dplyr::select(mouse_ID, death_date)
  }
  prep_path <- do.call(histo_prep_save, args)
  report <- report_knit_histo_qc(prep_path, exclusions)
  message("QC report: ", report)
  invisible(report)
}
