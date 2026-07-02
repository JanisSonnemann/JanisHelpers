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
