#' Match exclusion-sheet rows to the files of a prep result
#'
#' Resolves the `exclusions` sheet of an experiment's meta workbook (as
#' returned by [meta_read()]) against the files of one panel's prep data. An
#' exclusion can name a whole mouse or anything narrower, down to one `.fcs`
#' file: every non-empty column of a row must equal the file's value (exact
#' match, no case folding). Rows whose `panel` is empty apply to every panel;
#' rows for other panels are ignored.
#'
#' FACS samples are end-point only, so a row with a filled `collected_at`
#' (meant for one dated ELISA / histology sample) and an empty `panel` does
#' not apply to FACS; a dated row naming this FACS panel is an error.
#'
#' Used by the QC report (to flag samples) and by [db_export_facs()] (to drop
#' them). This function never errors on rows that match nothing -- it returns
#' them in `unmatched`, and the caller decides.
#'
#' Optional columns may be absent ([meta_read()] drops columns that are empty
#' in every row) or logical `NA`; both count as empty.
#'
#' @param data Long FACS tibble with columns `file_name`, `mouse_ID`,
#'   `tissue`, `condition` (the `data` of a prep result).
#' @param exclusions Tibble with required columns `mouse_ID` and `reason`,
#'   optional `tissue`, `panel`, `condition`, `file_name`, `collected_at`,
#'   `date`. `NULL` or zero rows = no exclusions.
#' @param panel Character scalar, the assay name of the panel being checked.
#' @returns A named list:
#'   \describe{
#'     \item{`applicable`}{the rows that apply to `panel`, with columns
#'       `mouse_ID`, `tissue`, `panel`, `condition`, `file_name`, `reason`
#'       (character), `collected_at` (Date, always empty here), `date` (Date)
#'       and `n_files_matched` (integer).}
#'     \item{`excluded_files`}{distinct `file_name` + `reason` (reasons of
#'       several rows hitting one file joined with `"; "`).}
#'     \item{`unmatched`}{the `applicable` rows that match no file.}
#'   }
#' @examples
#' \dontrun{
#' prep <- facs_prep_load(wsp, panel_dir)
#' ex <- meta_read("data/26-1_exp/meta/meta_26-1.xlsx")$exclusions
#' facs_match_exclusions(prep$data, ex, panel = "overview")
#' }
#' @export
facs_match_exclusions <- function(data, exclusions, panel) {
  match_cols <- c("mouse_ID", "tissue", "panel", "condition", "file_name")
  empty_applicable <- tibble::tibble(
    mouse_ID = character(), tissue = character(), panel = character(),
    condition = character(), file_name = character(), reason = character(),
    collected_at = as.Date(character()), date = as.Date(character()),
    n_files_matched = integer()
  )

  if (is.null(exclusions) || nrow(exclusions) == 0) {
    return(list(
      applicable = empty_applicable,
      excluded_files = tibble::tibble(file_name = character(), reason = character()),
      unmatched = empty_applicable
    ))
  }

  missing_cols <- setdiff(c("mouse_ID", "reason"), names(exclusions))
  if (length(missing_cols) > 0) {
    stop(
      "The exclusions sheet needs column(s): ",
      paste(missing_cols, collapse = ", "), ".",
      call. = FALSE
    )
  }
  for (col in c(match_cols, "collected_at", "date")) {
    if (!col %in% names(exclusions)) exclusions[[col]] <- NA
  }
  ex <- exclusions |>
    dplyr::transmute(
      dplyr::across(dplyr::all_of(c(match_cols, "reason")), as.character),
      collected_at = as.Date(collected_at),
      date = as.Date(date)
    ) |>
    dplyr::filter(is.na(panel) | panel == "" | panel == !!panel)

  # FACS samples are end-point only: a dated row targets an ELISA / histology
  # sample. Naming this panel with a date is a mistake; without a panel the
  # row simply doesn't apply here.
  dated <- !is.na(ex$collected_at)
  named <- !is.na(ex$panel) & ex$panel == panel
  if (any(dated & named)) {
    stop(
      "collected_at doesn't apply to FACS panel '", panel, "'; leave it empty ",
      "(row(s) for: ", paste(ex$mouse_ID[dated & named], collapse = ", "), ").",
      call. = FALSE
    )
  }
  ex <- ex[!dated, , drop = FALSE]

  bad <- is.na(ex$mouse_ID) | ex$mouse_ID == "" | is.na(ex$reason) | ex$reason == ""
  if (any(bad)) {
    stop(
      "Every exclusion row needs a mouse_ID and a reason; missing in row(s) ",
      "for: ", paste(ifelse(is.na(ex$mouse_ID[bad]), "<no mouse_ID>",
                            ex$mouse_ID[bad]), collapse = ", "), ".",
      call. = FALSE
    )
  }

  files <- dplyr::distinct(data, file_name, mouse_ID, tissue, condition) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  file_cols <- c("mouse_ID", "tissue", "condition", "file_name")

  hits <- purrr::map(seq_len(nrow(ex)), function(i) {
    keep <- rep(TRUE, nrow(files))
    for (col in file_cols) {
      wanted <- ex[[col]][[i]]
      if (!is.na(wanted)) keep <- keep & !is.na(files[[col]]) & files[[col]] == wanted
    }
    files$file_name[keep]
  })

  applicable <- dplyr::mutate(ex, n_files_matched = lengths(hits))
  excluded_files <- tibble::tibble(
    file_name = as.character(unlist(hits)),
    reason = rep(ex$reason, lengths(hits))
  ) |>
    dplyr::summarise(reason = paste(unique(reason), collapse = "; "), .by = file_name)

  list(
    applicable = applicable,
    excluded_files = excluded_files,
    unmatched = dplyr::filter(applicable, n_files_matched == 0L)
  )
}
