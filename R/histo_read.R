#' Read a curated histology sheet
#'
#' @description
#' Reads the curated histology layout: sheet \code{values} with
#' \code{mouse_ID}, optional \code{tissue}, optional \code{date}, and one
#' column per scored metric (numeric). Columns named \code{ID} or ending in
#' \code{_ID} (other than \code{mouse_ID}) are slide labels and ignored. Rows
#' without \code{mouse_ID} are dropped; mice with empty cells are kept with
#' \code{NA} values. Any non-numeric text in a metric cell is an error.
#'
#' An optional free-text column named \code{comment} or \code{comments} (any
#' case) is not a metric: it is returned as a \code{comment} column (repeated
#' on every metric row of the sample, \code{NA} when empty), shown in the QC
#' report and not exported to the database.
#'
#' A blank template with a README sheet ships with the package:
#' \code{system.file("templates", "histo_curated_template.xlsx",
#' package = "JanisHelpers")}.
#'
#' @param path path to the curated \code{.xlsx} file
#' @param tissue character; tissue used where the sheet has no \code{tissue}
#'   column or the cell is empty. Default \code{NA_character_} (error if the
#'   sheet gives no tissue either).
#'
#' @returns tibble, one row per mouse x metric, with columns \code{mouse_ID},
#'   \code{tissue}, \code{collected_at} (Date, \code{NA} = end-point),
#'   \code{metric}, \code{value} (+ \code{comment} when the sheet has a
#'   comment column); returned invisibly -- assign the result
#'   explicitly
#' @export
#'
#' @examples
#' \dontrun{
#'   histo_read_curated("kidney_26-1.xlsx", tissue = "kidney")
#' }
histo_read_curated <- function(path, tissue = NA_character_) {
  cv <- read_curated_values_(path, tissue)
  raw <- cv$raw
  mouse <- cv$mouse

  nm <- names(raw)
  id_col <- nm == "ID" | (grepl("_ID$", nm) & nm != "mouse_ID")
  comment_col <- !is.na(nm) & grepl("^comments?$", nm, ignore.case = TRUE)
  if (sum(comment_col) > 1L) {
    stop(glue::glue(
      "Sheet 'values' of '{path}' has more than one comment column: ",
      "{paste(nm[comment_col], collapse = ', ')} -- keep one."
    ), call. = FALSE)
  }
  metrics <- nm[!id_col & !comment_col & !nm %in% c("mouse_ID", "tissue", "date")]
  check_curated_headers_(metrics, path, "metric")
  if (length(metrics) == 0L) {
    stop(glue::glue("Sheet 'values' of '{path}' has no metric columns."))
  }

  base <- tibble::tibble(mouse_ID = mouse, tissue = cv$tissue,
                         collected_at = cv$collected_at)
  has_comment <- any(comment_col)
  if (has_comment) {
    cm <- stringr::str_trim(curated_cell_(raw[[which(comment_col)]]))
    cm[!is.na(cm) & !nzchar(cm)] <- NA_character_
  }
  long <- purrr::map(metrics, function(m) {
    txt <- stringr::str_trim(curated_cell_(raw[[m]]))
    num <- suppressWarnings(as.numeric(txt))
    bad <- !is.na(txt) & nzchar(txt) & is.na(num)
    if (any(bad)) {
      stop(glue::glue(
        "Cannot read value '{txt[bad][1]}' for metric '{m}', mouse '{mouse[bad][1]}': ",
        "expected a number."
      ), call. = FALSE)
    }
    out <- dplyr::mutate(base, metric = m, value = num)
    if (has_comment) out$comment <- cm
    out
  }) |>
    purrr::list_rbind()

  message(glue::glue(
    "\nExtraction Summary",
    "\n----------------------------------------------",
    "\nNumber of mice:    {dplyr::n_distinct(long$mouse_ID)}",
    "\nMetrics:           {paste(metrics, collapse = ', ')}",
    "\nNumber of rows:    {nrow(long)}",
    "\n----------------------------------------------\n"
  ))

  invisible(long)
}
