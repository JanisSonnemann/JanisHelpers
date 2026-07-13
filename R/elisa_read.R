#' Read back-calculated cytokine concentrations from a multiplex ELISA results file
#'
#' @description
#' Reads the sample-level ("Unknowns") sheet of a multiplex bead ELISA
#' (Luminex/Bio-Plex) \code{Results_<cytokine>.xlsx} export and extracts the
#' per-replicate back-calculated concentration. The samples sheet is located
#' by matching a sheet name ending in \code{"Unknowns"}, since some exports
#' name it plainly \code{"Unknowns"} while others prefix it with the
#' cytokine (e.g. \code{"TNFa_Unknowns"}). Rows with no \code{Sample} value
#' are dropped (a defensive safeguard against blank padding rows, a pattern
#' seen in these workbooks' other sheets). Out-of-range replicates
#' (\code{"OOR<"}, \code{"OOR>"} in the back-calculated column) become
#' \code{NA} in \code{value}, with the original flag preserved in
#' \code{result_status}.
#'
#' @param path path to a \code{Results_<cytokine>.xlsx} workbook
#' @param cytokine character; cytokine label for the \code{cytokine} output
#'   column. Default \code{NULL} derives it from \code{basename(path)} by
#'   stripping a leading \code{Results_} and trailing \code{.xlsx}.
#'
#' @returns tibble, one row per sample x replicate, with columns
#'   \code{cytokine}, \code{sample_id}, \code{replicate}, \code{value},
#'   \code{unit}, \code{result_status}; returned invisibly -- assign the
#'   result explicitly
#' @export
#'
#' @examples
#' \dontrun{
#'   elisa_read_results("Results_IL-17A.xlsx")
#' }
elisa_read_results <- function(path, cytokine = NULL) {
  if (is.null(cytokine)) {
    cytokine <- basename(path) |>
      stringr::str_remove("^Results_") |>
      stringr::str_remove("\\.xlsx$")
  }

  sheet_names <- readxl::excel_sheets(path)
  unknowns_sheet <- sheet_names[grepl("Unknowns$", sheet_names)]

  if (length(unknowns_sheet) != 1L) {
    stop(glue::glue(
      "Expected exactly one sheet ending in 'Unknowns' in '{path}', found ",
      "{length(unknowns_sheet)}: {paste(sheet_names, collapse = ', ')}."
    ))
  }

  raw <- readxl::read_excel(path, sheet = unknowns_sheet, .name_repair = "unique_quiet")

  backcalc_col <- names(raw)[grepl("^Backcalc", names(raw))]
  status_col   <- names(raw)[grepl("^Result.*Status", names(raw))]

  if (length(backcalc_col) != 1L) {
    stop(glue::glue(
      "Expected exactly one 'Backcalc' column in sheet '{unknowns_sheet}' ",
      "of '{path}', found {length(backcalc_col)}."
    ))
  }
  if (length(status_col) != 1L) {
    stop(glue::glue(
      "Expected exactly one 'Result Status' column in sheet '{unknowns_sheet}' ",
      "of '{path}', found {length(status_col)}."
    ))
  }

  unit <- stringr::str_extract(backcalc_col, "(?<=\\()[^)]+(?=\\))")
  if (is.na(unit)) {
    warning(glue::glue(
      "Could not parse a unit from Backcalc column header '{backcalc_col}'; ",
      "`unit` will be NA."
    ))
  }

  result <- raw |>
    dplyr::select(
      sample_id     = 1,
      replicate     = Rep,
      value         = dplyr::all_of(backcalc_col),
      result_status = dplyr::all_of(status_col)
    ) |>
    dplyr::filter(!is.na(sample_id)) |>
    dplyr::mutate(
      cytokine      = cytokine,
      replicate     = as.integer(replicate),
      value         = suppressWarnings(as.numeric(value)),
      unit          = unit,
      result_status = stringr::str_trim(result_status)
    ) |>
    dplyr::relocate(cytokine, sample_id, replicate, value, unit, result_status)

  message(glue::glue(
    "\nExtraction Summary",
    "\n----------------------------------------------",
    "\nCytokine:          {cytokine}",
    "\nNumber of samples: {dplyr::n_distinct(result$sample_id)}",
    "\nNumber of rows:    {nrow(result)}",
    "\n----------------------------------------------\n"
  ))

  invisible(result)
}
