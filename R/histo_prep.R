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
#' one mouse x tissue x date x metric occurs twice.
#'
#' @param data Tibble from [histo_read_curated()] (optionally with ratio
#'   metrics from [histo_calc_ratio()] and group / batch columns joined on) --
#'   needs `mouse_ID`, `tissue`, `collected_at`, `metric`, `value`.
#' @param path Path of the curated `.xlsx` the data came from.
#' @param assay Assay name the data is exported under.
#' @param group_col Optional column of `data` used to group the outlier plots.
#' @param batch_col Optional column of `data` for batch plots.
#' @param expected_mice Character vector of mice expected in the data.
#'   Default: the mice in `data`.
#' @param rules Character vector of plausibility rules for the QC report.
#'   Default `character()` (none).
#' @returns Invisibly, the path of the written `.rds`, a list with `data` and
#'   `info` (`source_path`, `root`, `assay`, `created_at`, `package_version`,
#'   `qc`).
#' @examples
#' \dontrun{
#' histo_read_curated("data/26-1_exp/raw/histo/kidney_26-1.xlsx", "kidney") |>
#'   histo_prep_save("data/26-1_exp/raw/histo/kidney_26-1.xlsx", "PAS")
#' }
#' @export
histo_prep_save <- function(data, path, assay, group_col = NULL,
                            batch_col = NULL, expected_mice = NULL,
                            rules = character()) {
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
  dups <- dplyr::count(
    data, dplyr::pick(dplyr::all_of(c("mouse_ID", "tissue", "collected_at", "metric")))
  ) |>
    dplyr::filter(n > 1)
  if (nrow(dups) > 0) {
    warning(
      "These mouse_ID x tissue x date x metric occur twice -- the export will ",
      "refuse them:\n", format_rows_(dups), call. = FALSE
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
