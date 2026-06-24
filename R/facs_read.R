# FlowJo WSP XML structure — verified against tests/fixtures/minimal.wsp
#
# Root sample tree node name (child of <Sample>):          SampleNode
# Population container element:                            Subpopulations
# Individual population element:                           Population
# Population count attribute:                              count
# Population name attribute:                               name
# Keyword container element (child of <Sample>):           Keywords (capital K)
# Keyword child element:                                   Keyword
# Keyword name/value attributes:                           name, value
# Statistic element name:                                  Statistic
# Statistic type attribute (Median/Mean/etc.):             name
# Statistic channel attribute:                             id  (NOT "channel")
# Statistic value attribute:                               value
# NOTE: Statistic nodes are siblings of Population inside Subpopulations, NOT children of Population
# SampleRef ID attribute:                                  sampleID
# DataSet sample ID attribute:                             sampleID
# Namespace prefix required for XPath:                     none (xml_ns_strip not needed)

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

filter_samples_ <- function(samples, ids) {
  if (is.null(ids)) return(samples)
  found_ids <- purrr::map_chr(samples, function(s) {
    xml2::xml_attr(xml2::xml_find_first(s, "DataSet"), "sampleID")
  })
  samples[found_ids %in% ids]
}

parse_keywords_ <- function(doc, sample_ids = NULL) {
  PANEL_PAT <- "^\\$P[0-9]+[NS]$"
  samples   <- filter_samples_(
    xml2::xml_find_all(doc, ".//SampleList/Sample"),
    sample_ids
  )

  purrr::map(samples, function(sample) {
    file_name <- basename(
      xml2::xml_attr(xml2::xml_find_first(sample, "DataSet"), "uri")
    )
    kw_nodes <- xml2::xml_find_all(sample, ".//Keywords/Keyword")
    tibble::tibble(
      file_name = file_name,
      key       = xml2::xml_attr(kw_nodes, "name"),
      value     = xml2::xml_attr(kw_nodes, "value")
    ) |>
      dplyr::filter(!grepl(PANEL_PAT, key))
  }) |>
    dplyr::bind_rows()
}

parse_panel_ <- function(doc, sample_ids = NULL) {
  samples <- filter_samples_(
    xml2::xml_find_all(doc, ".//SampleList/Sample"),
    sample_ids
  )
  purrr::map(samples, function(sample) {
    file_name <- basename(
      xml2::xml_attr(xml2::xml_find_first(sample, "DataSet"), "uri")
    )
    kw_nodes <- xml2::xml_find_all(sample, ".//Keywords/Keyword")
    kw_df <- tibble::tibble(
      key   = xml2::xml_attr(kw_nodes, "name"),
      value = xml2::xml_attr(kw_nodes, "value")
    )
    panel_n <- kw_df |>
      dplyr::filter(grepl("^\\$P[0-9]+N$", key)) |>
      dplyr::mutate(number = stringr::str_extract(key, "(?<=\\$P)[0-9]+")) |>
      dplyr::select(number, channel = value)
    panel_s <- kw_df |>
      dplyr::filter(grepl("^\\$P[0-9]+S$", key)) |>
      dplyr::mutate(number = stringr::str_extract(key, "(?<=\\$P)[0-9]+")) |>
      dplyr::select(number, stain = value)
    panel_long <- dplyr::left_join(panel_n, panel_s, by = "number") |>
      dplyr::mutate(stain = dplyr::na_if(stain, "")) |>
      dplyr::select(!number) |>
      tidyr::pivot_wider(names_from = channel, values_from = stain)
    tibble::add_column(panel_long, file_name = file_name, .before = 1L)
  }) |>
    dplyr::bind_rows()
}

walk_pops_ <- function(node, file_name, path, parent_count, stain_lookup) {
  stop("not yet implemented")
}

parse_populations_ <- function(doc, sample_ids = NULL) {
  stop("not yet implemented")
}

# ---------------------------------------------------------------------------
# Exported function
# ---------------------------------------------------------------------------

#' Read a FlowJo workspace into a structured list
#'
#' Parses a FlowJo \code{.wsp} file directly via \code{xml2} with no dependency
#' on \code{fcexpr}. All population statistics are always extracted.
#'
#' @param path path to \code{.wsp} file
#' @param group character; group name to extract. \code{NULL} (default) extracts
#'   all groups.
#' @param keywords character vector of FCS keyword names to join into
#'   \code{data}. Keywords absent from the workspace are filled
#'   \code{NA_character_} with a warning.
#'
#' @returns A named list with three elements:
#'   \describe{
#'     \item{\code{data}}{Long-format tibble, one row per
#'       \code{file_name x population_full_path x metric}.
#'       Columns: \code{file_name}, \code{population_full_path},
#'       \code{population}, \code{metric}, \code{value}, plus any
#'       requested \code{keywords}.}
#'     \item{\code{meta}}{Wide-format tibble, one row per file.
#'       Columns: \code{file_name}, \code{DATE}, \code{BTIM},
#'       \code{ETIM}, \code{CYT}, \code{INST}, \code{OP}, \code{TOT}.
#'       Missing system keywords filled \code{NA_character_}.}
#'     \item{\code{panel}}{Wide-format tibble, one row per file.
#'       One column per cytometer parameter (channel name); value is
#'       the stain label, or \code{NA} for unlabelled channels.}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#'   res <- facs_read_wsp(
#'     path     = "experiment.wsp",
#'     group    = "Spleen",
#'     keywords = c("mouse_ID", "group")
#'   )
#'   res$data
#'   res$meta
#'   res$panel
#' }
facs_read_wsp <- function(path, group = NULL, keywords = NULL) {
  stop("not yet implemented")
}
