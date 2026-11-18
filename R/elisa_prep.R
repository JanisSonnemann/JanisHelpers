# Internal: paths of the prep result and QC report of a curated ELISA file.
elisa_prep_paths_ <- function(path) {
  root <- tools::file_path_sans_ext(basename(path))
  list(
    root = root,
    prep = file.path(dirname(path), "processed", paste0(root, "_prep.rds")),
    qc = file.path(dirname(path), "qc", paste0(root, "_qc.html"))
  )
}

.elisa_prep_cols <- c("mouse_ID", "tissue", "collected_at", "analyte",
                      "replicate", "value", "result_status", "unit", "lloq",
                      "uloq")

#' Save an ELISA prep result for QC and export
#'
#' First step of the ELISA prep / QC loop. Stores the curated data (all
#' samples; exclusions are only flagged in the QC report) with the QC settings
#' at `<folder of path>/processed/<file name>_prep.rds`. The file name without
#' extension is the import unit ("root") stored in the database.
#'
#' @param data Tibble from [elisa_read_curated()] (optionally with group /
#'   batch columns joined on) -- needs `mouse_ID`, `tissue`, `collected_at`,
#'   `analyte`, `replicate`, `value`, `result_status`, `unit`, `lloq`, `uloq`.
#' @param path Path of the curated `.xlsx` the data came from.
#' @param assay Assay name the data is exported under.
#' @param group_col Optional column of `data` used to group the outlier plots.
#' @param batch_col Optional column of `data` for batch plots.
#' @param expected_mice Character vector of mice expected in the data.
#'   Default: the mice in `data`.
#' @param max_cv Replicate CV threshold in percent. Default 20.
#' @param log_scale Logical; log-scale the outlier plots. Default `FALSE`.
#' @returns Invisibly, the path of the written `.rds`, a list with `data` and
#'   `info` (`source_path`, `root`, `assay`, `created_at`, `package_version`,
#'   `qc`).
#' @examples
#' \dontrun{
#' elisa_read_curated("data/26-1_exp/raw/ELISA/spleen_26-1.xlsx", "spleen") |>
#'   elisa_prep_save("data/26-1_exp/raw/ELISA/spleen_26-1.xlsx", "restim")
#' }
#' @export
elisa_prep_save <- function(data, path, assay, group_col = NULL,
                            batch_col = NULL, expected_mice = NULL,
                            max_cv = 20, log_scale = FALSE) {
  missing_cols <- setdiff(.elisa_prep_cols, names(data))
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
  if (!is.numeric(max_cv) || length(max_cv) != 1 || is.na(max_cv) || max_cv < 0) {
    stop("`max_cv` must be a single non-negative number.", call. = FALSE)
  }
  paths <- elisa_prep_paths_(path)
  prep <- list(
    data = data,
    info = list(
      source_path = path,
      root = paths$root,
      assay = assay,
      created_at = Sys.time(),
      package_version = as.character(utils::packageVersion("JanisHelpers")),
      qc = list(
        group_col = group_col,
        batch_col = batch_col,
        expected_mice = if (is.null(expected_mice)) {
          sort(unique(as.character(data$mouse_ID)))
        } else {
          as.character(expected_mice)
        },
        max_cv = max_cv,
        log_scale = isTRUE(log_scale)
      )
    )
  )
  dir.create(dirname(paths$prep), recursive = TRUE, showWarnings = FALSE)
  saveRDS(prep, paths$prep)
  message("Prep result saved to ", paths$prep)
  invisible(paths$prep)
}

#' Load the saved ELISA prep result for the final export
#'
#' Reads the `.rds` written by [elisa_prep_save()] and checks it is current:
#' the QC report must exist and be newer than the prep result, and the curated
#' `.xlsx` must not have been changed since.
#'
#' @param path Path of the curated `.xlsx`, as passed to [elisa_prep_save()].
#' @returns The prep list (`data`, `info`) plus `paths` (`prep`, `qc`).
#' @examples
#' \dontrun{
#' prep <- elisa_prep_load("data/26-1_exp/raw/ELISA/spleen_26-1.xlsx")
#' }
#' @export
elisa_prep_load <- function(path) {
  paths <- elisa_prep_paths_(path)
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

#' Run the ELISA prep / QC loop for one curated file
#'
#' One call for the whole loop, driven by a `final` switch (see
#' `inst/scripts/elisa_prep_template.R`). `final = FALSE` reads the curated
#' sheet ([elisa_read_curated()]), applies `prepare()`, joins the group /
#' batch columns from `mice`, saves the prep result ([elisa_prep_save()]) and
#' renders the QC report ([report_knit_elisa_qc()]). `final = TRUE` loads the
#' saved prep ([elisa_prep_load()]), exports it ([db_export_elisa()]) and
#' optionally snapshots the database.
#'
#' @param final `FALSE` to prep and render the QC report, `TRUE` to export.
#' @param path Path of the curated `.xlsx`.
#' @param assay Assay name the data is exported under.
#' @param tissue Tissue of the file, passed to [elisa_read_curated()].
#' @param prepare Function applied to the long data tibble before saving.
#' @param mice Optional tibble with `mouse_ID` and the `group_col` /
#'   `batch_col` columns, joined into the data with [meta_annotate()].
#' @param exclusions Optional tibble -- the meta workbook's `exclusions`
#'   sheet. `NULL` = nothing excluded.
#' @param con A `DBI` connection from [db_connect()], required when
#'   `final = TRUE`; left open.
#' @param group_col,batch_col Character scalars or `NULL`, see
#'   [elisa_prep_save()].
#' @param ... Further QC settings passed to [elisa_prep_save()]
#'   (`expected_mice`, `max_cv`, `log_scale`). Used only when
#'   `final = FALSE`.
#' @param data_dir The project's `data/` folder, see [db_export_elisa()].
#' @param snapshot_root Folder for dated snapshots after an export. `NULL`
#'   (default) = none.
#' @param overwrite Passed to [db_export_elisa()].
#' @returns Invisibly: with `final = FALSE` the QC report path; with
#'   `final = TRUE` the new `import_id`, or `NULL` if the overwrite was
#'   declined.
#' @examples
#' \dontrun{
#' elisa_prep_run(
#'   FALSE, "data/26-1_exp/raw/ELISA/spleen_26-1.xlsx", "restim",
#'   tissue = "spleen", mice = meta$meta, group_col = "group"
#' )
#' }
#' @export
elisa_prep_run <- function(final, path, assay, tissue = NA_character_,
                           prepare = identity, mice = NULL, exclusions = NULL,
                           con = NULL, group_col = NULL, batch_col = NULL,
                           ..., data_dir = "data", snapshot_root = NULL,
                           overwrite = NULL) {
  if (!is.logical(final) || length(final) != 1L || is.na(final)) {
    stop("`final` must be TRUE or FALSE.", call. = FALSE)
  }
  if (final) {
    if (is.null(con)) {
      stop("`con` is required for the final export (final = TRUE).",
           call. = FALSE)
    }
    prep <- elisa_prep_load(path)
    id <- db_export_elisa(con, prep, exclusions, data_dir = data_dir,
                          overwrite = overwrite)
    if (!is.null(id) && !is.null(snapshot_root)) db_snapshot(con, snapshot_root)
    return(invisible(id))
  }

  if (!file.exists(path)) {
    stop("Curated file not found: ", path, call. = FALSE)
  }
  data <- elisa_read_curated(path, tissue) |>
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
  prep_path <- elisa_prep_save(data, path, assay, group_col = group_col,
                               batch_col = batch_col, ...)
  report <- report_knit_elisa_qc(prep_path, exclusions)
  message("QC report: ", report)
  invisible(report)
}
