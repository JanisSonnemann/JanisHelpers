is_posixct_ <- function(x) inherits(x, "POSIXct")

meta_clean_sheet_ <- function(dat) {
  dat <- janitor::remove_empty(dat, which = c("rows", "cols"))
  dat <- janitor::clean_names(dat)

  if ("mouse_id" %in% names(dat)) {
    dat <- dplyr::rename(dat, mouse_ID = mouse_id)
  }

  dat <- dplyr::mutate(
    dat,
    dplyr::across(dplyr::where(is.character), stringr::str_trim)
  )

  dat <- dplyr::mutate(
    dat,
    dplyr::across(dplyr::where(is_posixct_), as.Date)
  )

  if ("group" %in% names(dat)) {
    dat <- dplyr::mutate(dat, group = factor(group))
  }

  dat
}

#' Read and clean every sheet of an experiment metadata spreadsheet
#'
#' @description
#' Reads every sheet of an Excel workbook of experiment metadata and applies
#' standard cleaning to each: blank rows and columns removed, column names
#' standardized to snake_case (except \code{mouse_ID}, preserved verbatim so
#' it stays joinable against FlowJo keyword data such as
#' \code{facs_read_wsp(keywords = "mouse_ID")}), character columns trimmed
#' of whitespace, date/time columns coerced to \code{Date}, and a
#' \code{group} column (if present) coerced to a factor.
#'
#' @param path path to \code{.xlsx} metadata file
#'
#' @returns named list of cleaned tibbles, one per sheet, named after the
#'   sheet (e.g. \code{list(meta = ..., organ_weights = ..., facs_volumes =
#'   ...)}); returned invisibly -- assign the result explicitly
#' @export
#'
#' @examples
#' \dontrun{
#'   meta_list <- meta_read("meta.xlsx")
#'   meta_list$meta
#' }
meta_read <- function(path) {
  sheet_names <- readxl::excel_sheets(path)

  result <- sheet_names |>
    purrr::map(~ meta_clean_sheet_(readxl::read_excel(path, sheet = .x))) |>
    purrr::set_names(sheet_names)

  invisible(result)
}

#' Annotate experimental data with subject metadata
#'
#' @description
#' Left-joins a metadata tibble (typically from \code{meta_read()}) onto any
#' experimental data tibble (e.g. \code{facs_read_wsp(...)$data}) by one or
#' more shared identifier columns. Errors if any \code{by} column is missing
#' from either side, or if any non-\code{by} column names collide between
#' the two tibbles. Warns if any \code{by} combination present in
#' \code{data} has no match in \code{meta} (those rows keep \code{NA} for
#' all meta columns). \code{group} is a particularly likely candidate for
#' the column-collision error, since \code{facs_read_wsp(..., keywords =
#' "group")} and \code{meta_read()} both commonly produce a \code{group}
#' column.
#'
#' @param data tibble to annotate, e.g. \code{facs_read_wsp(...)$data}
#' @param meta metadata tibble, e.g. one element of \code{meta_read()}'s
#'   result
#' @param by character vector of shared identifier column name(s), default
#'   = \code{"mouse_ID"}. Pass e.g. \code{c("mouse_ID", "tissue")} to join
#'   on multiple columns.
#'
#' @returns \code{data} left-joined with \code{meta}, returned invisibly --
#'   assign the result explicitly
#' @export
#'
#' @examples
#' \dontrun{
#'   meta_list <- meta_read("meta.xlsx")
#'   dat <- meta_annotate(facs_read_wsp("experiment.wsp")$data, meta_list$meta)
#' }
meta_annotate <- function(data, meta, by = "mouse_ID") {
  missing_in_data <- setdiff(by, names(data))
  if (length(missing_in_data) > 0L) {
    stop(glue::glue(
      "Join column(s) not found in `data`: {paste(missing_in_data, collapse = ', ')}."
    ))
  }
  missing_in_meta <- setdiff(by, names(meta))
  if (length(missing_in_meta) > 0L) {
    stop(glue::glue(
      "Join column(s) not found in `meta`: {paste(missing_in_meta, collapse = ', ')}."
    ))
  }

  colliding <- setdiff(intersect(names(data), names(meta)), by)
  if (length(colliding) > 0L) {
    stop(glue::glue(
      "Column name(s) present in both `data` and `meta`: ",
      "{paste(colliding, collapse = ', ')}. ",
      "Rename or drop them before calling meta_annotate()."
    ))
  }

  joined <- dplyr::left_join(data, meta, by = by)

  data_keys <- dplyr::distinct(data, dplyr::across(dplyr::all_of(by)))
  meta_keys <- dplyr::distinct(meta, dplyr::across(dplyr::all_of(by)))
  unmatched <- dplyr::anti_join(data_keys, meta_keys, by = by)

  if (nrow(unmatched) > 0L) {
    unmatched_desc <- purrr::pmap_chr(
      unmatched,
      function(...) {
        row <- list(...)
        paste(paste0(names(row), "=", row), collapse = ", ")
      }
    )
    warning(glue::glue(
      "The following '{paste(by, collapse = ', ')}' combination(s) in `data` have no match in `meta`: ",
      "{paste(unmatched_desc, collapse = '; ')}"
    ))
  }

  invisible(joined)
}
