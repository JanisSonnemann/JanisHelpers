# Reshape a wide-by-date ELISA sheet (one row per mouse [x analyte x tissue],
# one column per sampling date, end-point column e.g. `final`) into the
# curated long layout read by elisa_read_curated().

#' Reshape a wide-by-date ELISA sheet into a curated long file
#'
#' @description
#' Optional step before \code{elisa_prep_run()}: converts a sheet with one row
#' per mouse (optionally per mouse x analyte x tissue) and one column per
#' sampling date into the curated layout read by \code{elisa_read_curated()}
#' (sheet \code{values}: \code{mouse_ID}, optional \code{tissue},
#' \code{date}, one column per analyte; empty \code{date} = end-point).
#'
#' Columns whose header is a date (a real Excel date or \code{YYYY-MM-DD}
#' text) are timepoints; columns named like \code{end_labels} (any case) are
#' the end-point. Optional \code{analyte} / \code{tissue} columns (any case)
#' give each row's analyte and tissue; the arguments of the same name fill a
#' missing column or empty cells. All other columns are dropped (listed in a
#' message). Date headers may be Excel dates (a time part is ignored),
#' \code{YYYY-MM-DD} or \code{dd.mm.yyyy} text. Repeated mouse x tissue x
#' analyte rows are replicates. Empty cells give no row (no sample is
#' created for an unsampled date). Rows without an ID are dropped, with a
#' warning if they hold values. Cell values are passed through unchanged, so
#' \code{elisa_read_curated()} still validates them.
#'
#' The output is only written when it is missing or its content would
#' change, so rerunning a prep script does not make \code{elisa_prep_load()}
#' see the file as newer than the saved prep.
#'
#' @param path path to the wide \code{.xlsx} file (first sheet is read)
#' @param analyte character; analyte name when the sheet has no
#'   \code{analyte} column (or for its empty cells). Default \code{NULL}.
#' @param tissue character; tissue when the sheet has no \code{tissue}
#'   column (or for its empty cells). Default \code{NULL}: no \code{tissue}
#'   column is written unless the sheet has one -- give the tissue to
#'   \code{elisa_prep_run()} instead.
#' @param id_col character; header(s) of the mouse ID column -- the first
#'   one present is used. Default \code{c("mouse_ID", "Lab ID")}.
#' @param out output path. Default \code{NULL}: \code{path} with a trailing
#'   \code{_wide} replaced by \code{_long} (else \code{_long} appended), same
#'   folder.
#' @param end_labels character; headers that mark the end-point column.
#' @param unit,lloq,uloq optional range information for a \code{range} sheet:
#'   one value for all analytes or a vector named by analyte. Default
#'   \code{NA} (no \code{range} sheet).
#'
#' @returns the output path, invisibly
#' @export
#'
#' @examples
#' \dontrun{
#'   path <- elisa_reshape_wide("anti-mpo-titer_26-1_wide.xlsx",
#'                              analyte = "anti-MPO-titer")
#' }
elisa_reshape_wide <- function(path, analyte = NULL, tissue = NULL,
                               id_col = c("mouse_ID", "Lab ID"), out = NULL,
                               end_labels = c("final", "end", "endpoint"),
                               unit = NA, lloq = NA, uloq = NA) {
  if (is.null(out)) {
    stem <- tools::file_path_sans_ext(basename(path))
    stem <- if (grepl("_wide$", stem)) sub("_wide$", "_long", stem) else paste0(stem, "_long")
    out <- file.path(dirname(path), paste0(stem, ".xlsx"))
  }
  raw <- readxl::read_excel(path, sheet = 1, col_types = "text",
                            .name_repair = "minimal")
  r <- elisa_reshape_wide_table_(raw, analyte = analyte, tissue = tissue,
                                 id_col = id_col, end_labels = end_labels,
                                 unit = unit, lloq = lloq, uloq = uloq)
  if (length(r$dropped) > 0L) {
    message(glue::glue("Dropped columns: {paste(r$dropped, collapse = ', ')}"))
  }
  if (r$n_empty > 0L) {
    message(glue::glue("Left out {r$n_empty} mouse x date row(s) without any value."))
  }
  sheets <- list(values = r$values)
  if (!is.null(r$range)) sheets$range <- r$range

  tmp <- tempfile(fileext = ".xlsx")
  on.exit(unlink(tmp), add = TRUE)
  writexl::write_xlsx(sheets, tmp)
  if (file.exists(out) && identical(read_xlsx_text_(tmp), read_xlsx_text_(out))) {
    message(glue::glue("{basename(out)} is up to date."))
    return(invisible(out))
  }
  if (!file.copy(tmp, out, overwrite = TRUE)) {
    stop(glue::glue("Could not write '{out}'."), call. = FALSE)
  }
  analytes <- setdiff(names(r$values), c("mouse_ID", "tissue", "date"))
  message(glue::glue(
    "Wrote {basename(out)}: {nrow(r$values)} rows, ",
    "{dplyr::n_distinct(r$values$mouse_ID)} mice, analytes: {paste(analytes, collapse = ', ')}."
  ))
  invisible(out)
}

# All sheets of an .xlsx as text, for comparing written content.
read_xlsx_text_ <- function(path) {
  s <- readxl::excel_sheets(path)
  stats::setNames(lapply(s, function(x) {
    as.list(readxl::read_excel(path, sheet = x, col_types = "text",
                               .name_repair = "minimal"))
  }), s)
}

# Header text -> Date: Excel serial numbers (as readxl returns date headers,
# five digits; a time part is dropped), YYYY-MM-DD or dd.mm.yyyy text;
# anything else -> NA.
parse_wide_date_header_ <- function(x) {
  x <- stringr::str_trim(as.character(x))
  out <- rep(NA_real_, length(x))
  serial <- !is.na(x) & grepl("^\\d{5}(\\.\\d+)?$", x)
  out[serial] <- as.numeric(as.Date(floor(as.numeric(x[serial])), origin = "1899-12-30"))
  iso <- !is.na(x) & grepl("^\\d{4}-\\d{2}-\\d{2}$", x)
  out[iso] <- as.numeric(as.Date(x[iso], format = "%Y-%m-%d"))
  de <- !is.na(x) & grepl("^\\d{1,2}\\.\\d{1,2}\\.\\d{4}$", x)
  out[de] <- as.numeric(as.Date(x[de], format = "%d.%m.%Y"))
  as.Date(out, origin = "1970-01-01")
}

# Expand `unit` / `lloq` / `uloq` (NA, one value for all analytes, or a vector
# named by analyte) to one character value per analyte.
range_values_ <- function(x, analytes, arg) {
  if (is.null(names(x))) {
    if (length(x) != 1L) {
      stop(glue::glue("`{arg}` must be one value or a vector named by analyte."),
           call. = FALSE)
    }
    return(rep(as.character(x), length(analytes)))
  }
  unknown <- setdiff(names(x), analytes)
  if (length(unknown) > 0L) {
    stop(glue::glue(
      "`{arg}` names analyte(s) not in the sheet: {paste(unknown, collapse = ', ')}."
    ), call. = FALSE)
  }
  unname(as.character(x[analytes]))
}

# `raw`: tibble of character columns (readxl col_types = "text", minimal
# names). Returns list(values, range, dropped) -- `values` in the curated
# layout (`date` = Date, NA = end-point), `range` NULL unless any of
# unit/lloq/uloq is given, `dropped` = ignored column names, `n_empty` =
# mouse x date rows left out because every analyte cell was empty.
elisa_reshape_wide_table_ <- function(raw, analyte = NULL, tissue = NULL,
                                      id_col = c("mouse_ID", "Lab ID"),
                                      end_labels = c("final", "end", "endpoint"),
                                      unit = NA, lloq = NA, uloq = NA) {
  nms <- names(raw)
  found <- id_col[id_col %in% nms]
  if (length(found) == 0L) {
    stop(glue::glue(
      "No ID column ({paste(id_col, collapse = ', ')}) in the sheet; found: ",
      "{paste(nms, collapse = ', ')}. Set `id_col`."
    ), call. = FALSE)
  }
  id_col <- found[1]
  low <- tolower(stringr::str_trim(nms))
  analyte_col <- nms[low == "analyte" & nms != id_col][1]
  tissue_col <- nms[low == "tissue" & nms != id_col][1]

  rest <- setdiff(nms, c(id_col, analyte_col, tissue_col))
  dates <- parse_wide_date_header_(rest)
  is_end <- tolower(stringr::str_trim(rest)) %in% tolower(end_labels)
  tp_cols <- rest[!is.na(dates) | is_end]
  dropped <- rest[is.na(dates) & !is_end]
  if (length(tp_cols) == 0L) {
    stop(glue::glue(
      "No timepoint columns found: headers must be dates or one of ",
      "{paste(end_labels, collapse = ', ')}."
    ), call. = FALSE)
  }
  tp_key <- ifelse(is_end, "end-point", format(dates))[!is.na(dates) | is_end]
  dup <- unique(tp_key[duplicated(tp_key)])
  if (length(dup) > 0L) {
    clash <- vapply(dup, function(k) paste(tp_cols[tp_key == k], collapse = " / "),
                    character(1))
    stop(glue::glue(
      "Several columns give the same timepoint: {paste(clash, collapse = '; ')}."
    ), call. = FALSE)
  }
  tp_date <- stats::setNames(as.Date(ifelse(tp_key == "end-point", NA, tp_key)),
                             tp_cols)

  mouse <- stringr::str_trim(raw[[id_col]])
  keep <- !is.na(mouse) & nzchar(mouse)
  lost <- sum(!keep & rowSums(!is_empty_cell_(raw[tp_cols])) > 0)
  if (lost > 0L) {
    warning(glue::glue(
      "Dropped {lost} row(s) with values but no ID in '{id_col}' ",
      "(merged ID cells?) -- fill in the ID on every row."
    ), call. = FALSE)
  }
  raw <- raw[keep, , drop = FALSE]
  mouse <- mouse[keep]

  fill <- function(col, arg) {
    v <- if (is.na(col)) rep(NA_character_, length(mouse)) else stringr::str_trim(raw[[col]])
    v[is.na(v) | !nzchar(v)] <- if (is.null(arg)) NA_character_ else arg
    v
  }
  an <- fill(analyte_col, analyte)
  if (any(is.na(an))) {
    stop("No analyte for some rows: give an 'analyte' column or the `analyte` argument.",
         call. = FALSE)
  }
  has_tissue <- !is.na(tissue_col) || !is.null(tissue)
  ti <- fill(tissue_col, tissue)

  analytes <- unique(an)
  mice <- unique(mouse)
  tissues <- unique(ti)
  base <- tibble::tibble(mouse_ID = mouse, tissue = ti, analyte = an,
                         dplyr::select(raw, dplyr::all_of(tp_cols))) |>
    dplyr::mutate(replicate = dplyr::row_number(),
                  .by = c("mouse_ID", "tissue", "analyte"))

  values <- base |>
    tidyr::pivot_longer(dplyr::all_of(tp_cols), names_to = "tp", values_to = "value") |>
    dplyr::mutate(date = unname(tp_date[.data$tp])) |>
    dplyr::select(!"tp") |>
    tidyr::pivot_wider(names_from = "analyte", values_from = "value") |>
    dplyr::arrange(match(.data$tissue, tissues), is.na(.data$date), .data$date,
                   match(.data$mouse_ID, mice), .data$replicate) |>
    dplyr::select("mouse_ID", "tissue", "date", dplyr::all_of(analytes))
  empty <- rowSums(!is_empty_cell_(values[analytes])) == 0
  values <- values[!empty, , drop = FALSE]
  if (!has_tissue) values <- dplyr::select(values, !"tissue")

  for (a in analytes) {
    txt <- values[[a]]
    num <- suppressWarnings(as.numeric(txt))
    if (all(is.na(txt) | !is.na(num))) values[[a]] <- num
  }

  range <- NULL
  if (!all(is.na(c(unit, lloq, uloq)))) {
    range <- tibble::as_tibble(
      c(list(c("LLOQ", "ULOQ", "unit")),
        stats::setNames(lapply(seq_along(analytes), function(i) {
          c(range_values_(lloq, analytes, "lloq")[i],
            range_values_(uloq, analytes, "uloq")[i],
            range_values_(unit, analytes, "unit")[i])
        }), analytes)),
      .name_repair = "minimal"
    )
    names(range)[1] <- ""
    range[-1] <- lapply(range[-1], function(x) { x[x %in% "NA"] <- NA; x })
  }

  list(values = values, range = range, dropped = dropped, n_empty = sum(empty))
}

# TRUE for NA or blank cells of a data frame of character columns (as matrix).
is_empty_cell_ <- function(df) {
  m <- as.matrix(df)
  is.na(m) | !nzchar(stringr::str_trim(m))
}
