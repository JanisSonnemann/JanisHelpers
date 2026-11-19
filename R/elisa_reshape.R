# Reshape a wide-by-date ELISA sheet (one row per mouse [x analyte x tissue],
# one column per sampling date, end-point column e.g. `final`) into the
# curated long layout read by elisa_read_curated().

# Header text -> Date: Excel serial numbers (as readxl returns date headers,
# five digits) or YYYY-MM-DD text; anything else -> NA.
parse_wide_date_header_ <- function(x) {
  x <- stringr::str_trim(as.character(x))
  out <- rep(NA_real_, length(x))
  serial <- !is.na(x) & grepl("^\\d{5}(\\.0+)?$", x)
  out[serial] <- as.numeric(as.Date(as.numeric(x[serial]), origin = "1899-12-30"))
  iso <- !is.na(x) & grepl("^\\d{4}-\\d{2}-\\d{2}$", x)
  out[iso] <- as.numeric(as.Date(x[iso], format = "%Y-%m-%d"))
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
# unit/lloq/uloq is given, `dropped` = ignored column names.
elisa_reshape_wide_table_ <- function(raw, analyte = NULL, tissue = NULL,
                                      id_col = "Lab ID",
                                      end_labels = c("final", "end", "endpoint"),
                                      unit = NA, lloq = NA, uloq = NA) {
  nms <- names(raw)
  if (!id_col %in% nms) {
    stop(glue::glue(
      "No ID column '{id_col}' in the sheet; found: {paste(nms, collapse = ', ')}. ",
      "Set `id_col`."
    ), call. = FALSE)
  }
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

  list(values = values, range = range, dropped = dropped)
}
