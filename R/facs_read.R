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
  subpops <- xml2::xml_find_first(node, "Subpopulations")
  if (inherits(subpops, "xml_missing")) return(tibble::tibble())

  skip_stats <- c("count", "freq. of parent", "freq of parent", "frequency of parent")

  purrr::map(xml2::xml_children(subpops), function(child) {
    nm <- xml2::xml_name(child)

    if (nm == "Population") {
      pop_name  <- xml2::xml_attr(child, "name")
      pop_count <- suppressWarnings(as.numeric(xml2::xml_attr(child, "count")))
      pop_path  <- if (nzchar(path)) paste0(path, "/", pop_name) else pop_name

      fop <- if (!is.na(parent_count) && parent_count > 0L) {
        pop_count / parent_count
      } else {
        NA_real_
      }

      base_rows <- tibble::tibble(
        file_name            = file_name,
        population_full_path = pop_path,
        population           = pop_name,
        metric               = c("count", "fraction_of_parent"),
        value                = c(pop_count, fop)
      )

      dplyr::bind_rows(
        base_rows,
        walk_pops_(child, file_name, pop_path, pop_count, stain_lookup)
      )

    } else if (nm == "Statistic") {
      stat_type    <- xml2::xml_attr(child, "name")
      stat_channel <- xml2::xml_attr(child, "id")
      stat_value   <- suppressWarnings(as.numeric(xml2::xml_attr(child, "value")))

      if (is.na(stat_type) || tolower(stat_type) %in% skip_stats) return(tibble::tibble())
      if (is.na(stat_channel) || !nzchar(stat_channel))           return(tibble::tibble())

      matched <- stain_lookup$label[stain_lookup$channel == stat_channel]
      label   <- if (length(matched) > 0L && !is.na(matched[[1L]])) matched[[1L]] else stat_channel

      tibble::tibble(
        file_name            = file_name,
        population_full_path = path,
        population           = basename(path),
        metric               = paste0(tolower(stat_type), "_", label),
        value                = stat_value
      )

    } else {
      tibble::tibble()
    }
  }) |>
    dplyr::bind_rows()
}

parse_populations_ <- function(doc, sample_ids = NULL) {
  samples <- filter_samples_(
    xml2::xml_find_all(doc, ".//SampleList/Sample"),
    sample_ids
  )

  purrr::map(samples, function(sample) {
    file_name <- basename(
      xml2::xml_attr(xml2::xml_find_first(sample, "DataSet"), "uri")
    )

    # Build stain lookup: tibble(channel, label) for this sample
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
    stain_lookup <- dplyr::left_join(panel_n, panel_s, by = "number") |>
      dplyr::mutate(
        stain = dplyr::na_if(stain, ""),
        label = dplyr::coalesce(stain, channel)
      ) |>
      dplyr::select(channel, label)

    # Root tree node
    root_node  <- xml2::xml_find_first(sample, "SampleNode")
    root_count <- suppressWarnings(as.numeric(xml2::xml_attr(root_node, "count")))

    walk_pops_(root_node, file_name, "", root_count, stain_lookup)
  }) |>
    dplyr::bind_rows()
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
