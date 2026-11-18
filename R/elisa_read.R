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
      sample_id     = as.character(sample_id),
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

#' Read a curated ELISA sheet
#'
#' @description
#' Reads the curated ELISA layout: sheet \code{values} with \code{mouse_ID},
#' optional \code{tissue}, optional \code{date}, and one column per analyte
#' (numeric, or \code{OOR<} / \code{OOR>}). Repeated
#' \code{mouse_ID x tissue x date} rows are technical replicates, numbered in
#' row order. The optional sheet \code{range} holds the row labels
#' \code{LLOQ}, \code{ULOQ} and optionally \code{unit} in its first column and
#' one column per analyte. Rows without \code{mouse_ID} are dropped. Empty
#' analyte cells become \code{NA} with status \code{"ok"}; any other
#' non-numeric text is an error.
#'
#' @param path path to the curated \code{.xlsx} file
#' @param tissue character; tissue used where the sheet has no \code{tissue}
#'   column or the cell is empty. Default \code{NA_character_} (error if the
#'   sheet gives no tissue either).
#'
#' @returns tibble, one row per sample x analyte x replicate, with columns
#'   \code{mouse_ID}, \code{tissue}, \code{collected_at} (Date, \code{NA} =
#'   end-point), \code{analyte}, \code{replicate}, \code{value} (\code{NA} for
#'   out-of-range), \code{result_status} (\code{"ok"}, \code{"OOR<"},
#'   \code{"OOR>"}), \code{unit}, \code{lloq}, \code{uloq}; returned
#'   invisibly -- assign the result explicitly
#' @export
#'
#' @examples
#' \dontrun{
#'   elisa_read_curated("spleen-restim-48h_26-1.xlsx", tissue = "spleen")
#' }
elisa_read_curated <- function(path, tissue = NA_character_) {
  sheets <- readxl::excel_sheets(path)
  if (!"values" %in% sheets) {
    stop(glue::glue(
      "No 'values' sheet in '{path}'; found: {paste(sheets, collapse = ', ')}."
    ))
  }
  raw <- readxl::read_excel(path, sheet = "values", col_types = "list",
                            .name_repair = "minimal")
  cell <- function(col) {
    vapply(col, function(z) if (length(z) == 0L || is.na(z[[1]])) NA_character_ else as.character(z[[1]]), character(1))
  }
  if (!"mouse_ID" %in% names(raw)) {
    stop(glue::glue("Sheet 'values' of '{path}' has no 'mouse_ID' column."))
  }
  mouse <- stringr::str_trim(cell(raw$mouse_ID))
  keep <- !is.na(mouse) & nzchar(mouse)
  raw <- raw[keep, , drop = FALSE]
  mouse <- mouse[keep]

  tissue_chr <- rep(as.character(tissue), length(mouse))
  if ("tissue" %in% names(raw)) {
    tc <- stringr::str_trim(cell(raw$tissue))
    ok <- !is.na(tc) & nzchar(tc)
    tissue_chr[ok] <- tc[ok]
  }
  if (any(is.na(tissue_chr))) {
    stop(glue::glue(
      "No tissue for mouse(s) {paste(unique(mouse[is.na(tissue_chr)]), collapse = ', ')}: ",
      "give a 'tissue' column or the `tissue` argument."
    ))
  }

  collected_at <- if ("date" %in% names(raw)) {
    parse_curated_dates_(raw$date, mouse)
  } else {
    as.Date(rep(NA_real_, length(mouse)))
  }

  analytes <- setdiff(names(raw), c("mouse_ID", "tissue", "date"))
  if (length(analytes) == 0L) {
    stop(glue::glue("Sheet 'values' of '{path}' has no analyte columns."))
  }

  # range sheet
  rng <- tibble::tibble(analyte = analytes, lloq = NA_real_, uloq = NA_real_,
                        unit = NA_character_)
  if ("range" %in% sheets) {
    r <- readxl::read_excel(path, sheet = "range", col_types = "text",
                            .name_repair = "minimal")
    lab <- toupper(stringr::str_trim(r[[1]]))
    pick <- function(l) {
      i <- which(lab == l)
      if (length(i) == 0L) return(rep(NA_character_, length(analytes)))
      vapply(analytes, function(a) if (a %in% names(r)) r[[a]][i[1]] else NA_character_, character(1))
    }
    rng$lloq <- suppressWarnings(as.numeric(pick("LLOQ")))
    rng$uloq <- suppressWarnings(as.numeric(pick("ULOQ")))
    rng$unit <- unname(pick("UNIT"))
  }

  base <- tibble::tibble(mouse_ID = mouse, tissue = tissue_chr,
                         collected_at = collected_at)
  base <- base |>
    dplyr::group_by(.data$mouse_ID, .data$tissue, .data$collected_at) |>
    dplyr::mutate(replicate = dplyr::row_number()) |>
    dplyr::ungroup()

  long <- purrr::map(analytes, function(a) {
    txt <- stringr::str_trim(cell(raw[[a]]))
    status <- rep("ok", length(txt))
    status[txt %in% c("OOR<", "OOR>")] <- txt[txt %in% c("OOR<", "OOR>")]
    num <- suppressWarnings(as.numeric(txt))
    bad <- !is.na(txt) & nzchar(txt) & status == "ok" & is.na(num)
    if (any(bad)) {
      stop(glue::glue(
        "Cannot read value '{txt[bad][1]}' for analyte '{a}', mouse '{mouse[bad][1]}': ",
        "expected a number, 'OOR<' or 'OOR>'."
      ))
    }
    num[status != "ok"] <- NA_real_
    dplyr::mutate(base, analyte = a, value = num, result_status = status)
  }) |>
    purrr::list_rbind() |>
    dplyr::left_join(rng, by = "analyte") |>
    dplyr::select(
      "mouse_ID", "tissue", "collected_at", "analyte", "replicate", "value",
      "result_status", "unit", "lloq", "uloq"
    )

  message(glue::glue(
    "\nExtraction Summary",
    "\n----------------------------------------------",
    "\nNumber of mice:    {dplyr::n_distinct(long$mouse_ID)}",
    "\nAnalytes:          {paste(analytes, collapse = ', ')}",
    "\nNumber of rows:    {nrow(long)}",
    "\nOut of range:      {sum(long$result_status != 'ok')}",
    "\n----------------------------------------------\n"
  ))

  invisible(long)
}

# Convert the cells of a curated `date` column to Date. Accepts POSIXct, Date,
# Excel serial numbers and ISO text; empty -> NA; anything else is an error
# naming the mouse and the value. `x` is a vector or a list of cells.
parse_curated_dates_ <- function(x, mouse_ID) {
  cells <- if (is.list(x)) x else as.list(x)
  out <- vapply(seq_along(cells), function(i) {
    z <- cells[[i]]
    if (length(z) == 0L || is.na(z[[1]])) return(NA_real_)
    z <- z[[1]]
    if (inherits(z, "POSIXct")) return(as.numeric(as.Date(z, tz = "UTC")))
    if (inherits(z, "Date")) return(as.numeric(z))
    if (is.numeric(z)) return(as.numeric(as.Date(z, origin = "1899-12-30")))
    txt <- stringr::str_trim(as.character(z))
    if (!nzchar(txt)) return(NA_real_)
    d <- if (grepl("^\\d{4}-\\d{2}-\\d{2}$", txt)) {
      as.Date(txt, format = "%Y-%m-%d")
    } else {
      as.Date(NA)
    }
    if (is.na(d)) {
      stop(glue::glue(
        "Cannot read date '{txt}' for mouse '{mouse_ID[i]}': use a real Excel date or YYYY-MM-DD."
      ), call. = FALSE)
    }
    as.numeric(d)
  }, numeric(1))
  as.Date(out, origin = "1970-01-01")
}
