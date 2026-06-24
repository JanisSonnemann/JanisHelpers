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
  stop("not yet implemented")
}

parse_keywords_ <- function(doc, sample_ids = NULL) {
  stop("not yet implemented")
}

parse_panel_ <- function(doc, sample_ids = NULL) {
  stop("not yet implemented")
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
