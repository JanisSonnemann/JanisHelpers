#' Import FlowJo workspace data in long format
#'
#' @description
#' **Deprecated.** Use \code{facs_read_wsp} instead.
#'
#' @param path path to .wsp file
#' @param group group to extract from workspace, default = NULL (all groups)
#' @param r_stats logical, whether to extract statistics such as MFI, default = FALSE
#' @param keywords character vector of FCS keywords to attach to each row, default = NULL
#'
#' @returns tibble in long format with one row per file x population x metric,
#'   returned invisibly -- assign the result explicitly
#' @export
#'
#' @examples
#' \dontrun{
#'   dat <- facs_import_wsp(
#'     path     = "experiment.wsp",
#'     group    = "Spleen",
#'     r_stats  = TRUE,
#'     keywords = c("mouse_ID", "group")
#'   )
#' }
facs_import_wsp <- function(path, group = NULL, r_stats = FALSE, keywords = NULL) {
  .Deprecated("facs_read_wsp")

  # Import raw workspace; fall back to legacy importer if FCS files were renamed after export
  ps_raw <- tryCatch(
    fcexpr::wsx_get_popstats(ws = path, return_stats = r_stats, groups = group),
    error = function(e) {
      if (grepl("object 'out' not found", conditionMessage(e), fixed = TRUE)) {
        message("wsx_get_popstats failed (likely renamed FCS files) \u2014 retrying with legacy importer.")
        fcexpr::wsx_get_popstats_legacy(ws = path, return_stats = r_stats, groups = group)
      } else {
        stop(e)
      }
    }
  )

  # Extract counts and pivot to long
  counts <- ps_raw[["counts"]] |>
    dplyr::select(FileName, PopulationFullPath, Population, Count, FractionOfParent) |>
    tidyr::pivot_longer(
      cols      = c(Count, FractionOfParent),
      names_to  = "metric",
      values_to = "value"
    )

  # Only parse keywords from workspace when they are actually needed
  if (!is.null(keywords) || r_stats) {
    keys <- fcexpr::wsx_get_keywords(ws = path, return = "data.frame") |>
      tibble::enframe() |>
      dplyr::rename(FileName = name) |>
      tidyr::unnest(c(FileName, value)) |>
      dplyr::rename(key = name)
  }

  # Build stats if requested
  if (r_stats) {

    # Extract channel-to-stain lookup from $PnN / $PnS keywords
    stains <- keys |>
      dplyr::filter(grepl("\\$P[0-9]+[NS]", key)) |>
      dplyr::mutate(
        type   = dplyr::case_when(
          grepl("\\$P[0-9]+N", key) ~ "channel",
          grepl("\\$P[0-9]+S", key) ~ "stain"
        ),
        number = stringr::str_extract(key, "(?<=P)[0-9]+")
      ) |>
      dplyr::select(!key) |>
      tidyr::pivot_wider(names_from = type, values_from = value) |>
      dplyr::select(!number)

    stats <- ps_raw[["stats"]] |>
      dplyr::left_join(stains, by = c("FileName", "channel")) |>
      dplyr::mutate(
        stain      = dplyr::na_if(stain, ""),
        label      = dplyr::coalesce(stain, channel),
        metric     = dplyr::if_else(
          !is.na(statistic) & nzchar(statistic),
          paste0(statistic, "_", label),
          NA_character_
        ),
        Population = basename(PopulationFullPath)
      ) |>
      dplyr::select(FileName, PopulationFullPath, Population, metric, value)

    df <- dplyr::bind_rows(counts, stats)

  } else {
    df <- counts
  }

  # Filter, widen, and validate requested keywords
  if (!is.null(keywords) && length(keywords) > 0) {

    keys_clean <- keys |>
      dplyr::filter(key %in% keywords) |>
      tidyr::pivot_wider(names_from = key, values_from = value)

    missing_cols <- setdiff(keywords, names(keys_clean))
    if (length(missing_cols) > 0) {
      warning(
        "The following requested keywords were not found in the workspace and were filled with NA: ",
        paste(missing_cols, collapse = ", ")
      )
      keys_clean[missing_cols] <- NA_character_
    }

    df <- dplyr::left_join(df, keys_clean, by = "FileName")
  }

  # Tidy column order \u2014 core flow columns first, keyword columns follow
  final <- df |>
    dplyr::relocate(FileName, PopulationFullPath, Population, metric, value)

  # Print extraction summary
  n_files <- dplyr::n_distinct(final$FileName)

  keyword_summary <- if (!is.null(keywords) && length(keywords) > 0) {
    purrr::map_chr(keywords, function(k) paste0(k, ": ", dplyr::n_distinct(final[[k]]))) |>
      paste(collapse = "\n")
  } else {
    "none requested"
  }

  message(glue::glue(
    "\nExtraction Summary",
    "\n----------------------------------------------",
    "\nExtracted groups:          {if (is.null(group)) 'all' else group}",
    "\nNumber of samples:         {n_files}",
    "\nKeyword distinct values:\n{keyword_summary}",
    "\n----------------------------------------------\n"
  ))

  invisible(final)
}
