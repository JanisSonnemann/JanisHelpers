is_posixct_ <- function(x) inherits(x, "POSIXct")

#' Read and clean an experiment metadata spreadsheet
#'
#' @description
#' Reads an Excel sheet of per-subject experiment metadata (e.g. mouse ID,
#' cage, sex, group, dates) and applies standard cleaning: blank rows and
#' columns removed, column names standardized to snake_case (except
#' \code{mouse_ID}, preserved verbatim so it stays joinable against FlowJo
#' keyword data such as \code{facs_read_wsp(keywords = "mouse_ID")}),
#' character columns trimmed of whitespace, date/time columns coerced to
#' \code{Date}, and a \code{group} column (if present) coerced to a factor.
#'
#' @param path path to \code{.xlsx} metadata file
#' @param sheet sheet name or index to read, default = 1
#'
#' @returns cleaned tibble, one row per subject, returned invisibly --
#'   assign the result explicitly
#' @export
#'
#' @examples
#' \dontrun{
#'   meta <- meta_read("meta.xlsx")
#' }
meta_read <- function(path, sheet = 1) {
  dat <- readxl::read_excel(path, sheet = sheet)
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

  invisible(dat)
}

#' Annotate experimental data with subject metadata
#'
#' @description
#' Left-joins a metadata tibble (typically from \code{meta_read()}) onto any
#' experimental data tibble (e.g. \code{facs_read_wsp(...)$data}) by a shared
#' identifier column. Errors if the join column is missing from either side,
#' or if any non-join column names collide between the two tibbles. Warns if
#' any \code{by} value present in \code{data} has no match in \code{meta}
#' (those rows keep \code{NA} for all meta columns). \code{group} is a
#' particularly likely candidate for this column-collision error, since
#' \code{facs_read_wsp(..., keywords = "group")} and \code{meta_read()} both
#' commonly produce a \code{group} column.
#'
#' @param data tibble to annotate, e.g. \code{facs_read_wsp(...)$data}
#' @param meta metadata tibble, e.g. from \code{meta_read()}
#' @param by name of the shared identifier column, default = "mouse_ID"
#'
#' @returns \code{data} left-joined with \code{meta}, returned invisibly --
#'   assign the result explicitly
#' @export
#'
#' @examples
#' \dontrun{
#'   dat <- meta_annotate(facs_read_wsp("experiment.wsp")$data, meta_read("meta.xlsx"))
#' }
meta_annotate <- function(data, meta, by = "mouse_ID") {
  if (!by %in% names(data)) {
    stop(glue::glue("Join column '{by}' not found in `data`."))
  }
  if (!by %in% names(meta)) {
    stop(glue::glue("Join column '{by}' not found in `meta`."))
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

  unmatched <- setdiff(unique(data[[by]]), unique(meta[[by]]))
  if (length(unmatched) > 0L) {
    warning(glue::glue(
      "The following '{by}' values in `data` have no match in `meta`: ",
      "{paste(unmatched, collapse = ', ')}"
    ))
  }

  invisible(joined)
}
