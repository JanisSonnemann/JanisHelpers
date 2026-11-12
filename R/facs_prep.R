# Internal: standard locations of a workspace version's prep result and QC
# report inside the panel folder (the folder holding the panel's raw .fcs).
facs_prep_paths_ <- function(wsp_path, panel_dir) {
  base <- tools::file_path_sans_ext(basename(wsp_path))
  list(
    prep = file.path(panel_dir, "processed", paste0(base, "_prep.rds")),
    qc = file.path(panel_dir, "qc", paste0(base, "_qc.html"))
  )
}

# Internal: split "<root>_YYYYMMDD.wsp" into root + version date.
parse_wsp_name_ <- function(wsp_path) {
  file <- basename(wsp_path)
  m <- regmatches(file, regexec("^(.+)_([0-9]{8})\\.wsp$", file))[[1]]
  version <- if (length(m) == 3) as.Date(m[[3]], "%Y%m%d") else NA
  if (length(m) != 3 || is.na(version)) {
    stop(
      "Workspace '", file, "' isn't named '<root>_YYYYMMDD.wsp'.",
      call. = FALSE
    )
  }
  list(root = m[[2]], version = version)
}

#' Find FCS files that share one mouse x tissue x condition
#'
#' Each mouse x tissue x condition should come from exactly one FCS file. Two
#' files with the same combination usually mean a mistyped FlowJo keyword; if
#' they reached the database they would be merged into one stain.
#'
#' @param data Long FACS tibble with columns `file_name`, `mouse_ID`,
#'   `tissue`, `condition` (e.g. the `data` of a prep result from
#'   [facs_prep_save()]).
#' @returns A tibble with columns `file_name`, `mouse_ID`, `tissue`,
#'   `condition`: one row per file whose combination is shared with another
#'   file. Zero rows if there are no clashes.
#' @export
facs_check_file_clashes <- function(data) {
  data |>
    dplyr::distinct(file_name, mouse_ID, tissue, condition) |>
    dplyr::filter(dplyr::n() > 1, .by = c(mouse_ID, tissue, condition)) |>
    dplyr::arrange(mouse_ID, tissue, condition, file_name)
}

#' Save a prep result for one workspace version
#'
#' Step one of the prep / QC loop: after reading a workspace with
#' [facs_read_wsp()], selecting populations and adding calculated metrics,
#' save the result so the QC report ([report_knit_facs_qc()]) and later the
#' final export ([db_export_facs()]) use exactly this data. All samples are
#' kept -- exclusions are applied only at export. Files without a `mouse_ID`
#' or `tissue` keyword (compensation, FMO) are dropped.
#'
#' The file is written to `<panel_dir>/processed/<wsp name>_prep.rds` (one per
#' workspace version, overwritten on rerun). Warns -- but still saves -- if
#' several files share one mouse x tissue x condition (see
#' [facs_check_file_clashes()]).
#'
#' @param data Long tibble with columns `file_name`, `population_full_path`,
#'   `population`, `metric`, `value`, `mouse_ID`, `tissue` (and the
#'   `condition_keyword` column, if given) -- the `data` element of
#'   [facs_read_wsp()] after selection and `facs_calc_*()`.
#' @param wsp_path Path to the workspace, named `<root>_YYYYMMDD.wsp`.
#' @param assay Character scalar, the panel's assay name in the database
#'   (e.g. `"overview"`).
#' @param panel_dir Path to the panel's folder of raw `.fcs` files (e.g.
#'   `"data/26-1_exp/raw/FACS/overview"`); `processed/` is created inside it.
#' @param condition_keyword Character scalar -- for panels with several wells
#'   per mouse x tissue that differ in one variable, the keyword column holding
#'   it (e.g. `"restimulation"`). `NA` (default) = one well, condition
#'   `"none"`.
#' @returns Invisibly, the path of the written `.rds`. It holds a list with
#'   `data` (the tibble plus a `condition` column) and `info` (`wsp_path`,
#'   `wsp_root`, `wsp_version`, `assay`, `condition_keyword`, `panel_dir`,
#'   `created_at`, `package_version`).
#' @examples
#' \dontrun{
#' wsp <- "data/26-1_exp/raw/FACS/overview-26-1_20260806.wsp"
#' facs_read_wsp(wsp, keywords = c("mouse_ID", "tissue"))$data |>
#'   facs_calc_pct_of(ref_pop = "live_CD45+") |>
#'   facs_prep_save(wsp, "overview", "data/26-1_exp/raw/FACS/overview")
#' }
#' @export
facs_prep_save <- function(data, wsp_path, assay, panel_dir,
                           condition_keyword = NA_character_) {
  kw <- condition_keyword
  needed <- c("file_name", "population_full_path", "population", "metric",
              "value", "mouse_ID", "tissue", if (!is.na(kw)) kw)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "`data` is missing column(s): ", paste(missing_cols, collapse = ", "), ".",
      call. = FALSE
    )
  }
  wsp <- parse_wsp_name_(wsp_path)

  data <- dplyr::filter(data, !is.na(mouse_ID), !is.na(tissue))
  data$condition <- if (is.na(kw)) "none" else as.character(data[[kw]])
  if (anyNA(data$condition)) {
    stop("Some files have no '", kw, "' keyword value.", call. = FALSE)
  }

  clashes <- facs_check_file_clashes(data)
  if (nrow(clashes) > 0) {
    warning(
      "Several FCS files share one mouse_ID x tissue x condition -- fix the ",
      "keywords in FlowJo before the final export:\n",
      paste(utils::capture.output(print(clashes, n = Inf)), collapse = "\n"),
      call. = FALSE
    )
  }

  prep <- list(
    data = data,
    info = list(
      wsp_path = wsp_path,
      wsp_root = wsp$root,
      wsp_version = wsp$version,
      assay = assay,
      condition_keyword = kw,
      panel_dir = panel_dir,
      created_at = Sys.time(),
      package_version = as.character(utils::packageVersion("JanisHelpers"))
    )
  )
  path <- facs_prep_paths_(wsp_path, panel_dir)$prep
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(prep, path)
  message("Prep result saved to ", path)
  invisible(path)
}

#' Load the saved prep result for the final export
#'
#' Reads the `.rds` written by [facs_prep_save()] for one workspace version
#' and checks that it is current: the QC report must exist and the workspace
#' must not have been saved since the prep result was made (otherwise the
#' report shows outdated gates).
#'
#' @param wsp_path Path to the workspace, as passed to [facs_prep_save()].
#' @param panel_dir Path to the panel folder, as passed to [facs_prep_save()].
#' @returns The prep list (`data`, `info`) plus `paths` (`prep`, `qc`).
#' @examples
#' \dontrun{
#' prep <- facs_prep_load(
#'   "data/26-1_exp/raw/FACS/overview-26-1_20260806.wsp",
#'   "data/26-1_exp/raw/FACS/overview"
#' )
#' }
#' @export
facs_prep_load <- function(wsp_path, panel_dir) {
  paths <- facs_prep_paths_(wsp_path, panel_dir)
  for (p in c(paths$prep, paths$qc)) {
    if (!file.exists(p)) {
      stop(
        "'", p, "' not found. Run the prep script with final <- FALSE first.",
        call. = FALSE
      )
    }
  }
  if (file.mtime(wsp_path) > file.mtime(paths$prep)) {
    stop(
      "'", basename(wsp_path), "' is newer than its prep result -- rerun the ",
      "prep script with final <- FALSE and check the QC report.",
      call. = FALSE
    )
  }
  prep <- readRDS(paths$prep)
  prep$paths <- paths
  prep
}
