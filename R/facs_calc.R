#' Compute a population's percentage of an arbitrary reference population
#'
#' @description
#' Computes each population's count as a fraction of a named ancestor
#' population's count, regardless of how many gating levels separate them
#' (unlike \code{fraction_of_parent} from \code{facs_read_wsp()}, which is
#' always relative to the immediate parent gate).
#'
#' @param data tibble shaped like \code{facs_read_wsp(...)$data}: must
#'   contain \code{file_name}, \code{population_full_path}, \code{population},
#'   \code{metric}, \code{value}.
#' @param ref_pop character; leaf population name (matches \code{population},
#'   not \code{population_full_path}) to use as the denominator.
#'
#' @returns \code{data} with additional rows appended: one row per
#'   \code{file_name x population} (excluding \code{ref_pop} itself), with
#'   \code{metric = paste0("pct_of_", ref_pop)} and \code{value} as a 0-1
#'   fraction of \code{ref_pop}'s count in that file. Errors if \code{ref_pop}
#'   matches more than one population per file; warns and fills \code{NA} if
#'   \code{ref_pop} has no match for a file.
#' @export
#'
#' @examples
#' \dontrun{
#'   dat <- facs_read_wsp("experiment.wsp")$data
#'   facs_calc_pct_of(dat, ref_pop = "CD45+")
#' }
facs_calc_pct_of <- function(data, ref_pop) {
  ref_counts <- data |>
    dplyr::filter(population == ref_pop, metric == "count") |>
    dplyr::select(file_name, ref_count = value)

  dup_files <- unique(ref_counts$file_name[duplicated(ref_counts$file_name)])
  if (length(dup_files) > 0L) {
    stop(glue::glue(
      "ref_pop '{ref_pop}' matches more than one population for file_name(s): ",
      "{paste(dup_files, collapse = ', ')}. Leaf population name is ambiguous ",
      "\u2014 rename the gate(s) in FlowJo or disambiguate before calling facs_calc_pct_of()."
    ))
  }

  missing_files <- setdiff(unique(data$file_name), ref_counts$file_name)
  if (length(missing_files) > 0L) {
    warning(glue::glue(
      "ref_pop '{ref_pop}' not found for file_name(s): ",
      "{paste(missing_files, collapse = ', ')}. Result filled with NA."
    ))
  }

  new_rows <- data |>
    dplyr::filter(metric == "count", population != ref_pop) |>
    dplyr::left_join(ref_counts, by = "file_name") |>
    dplyr::mutate(
      metric = paste0("pct_of_", ref_pop),
      value  = value / ref_count
    ) |>
    dplyr::select(dplyr::all_of(names(data)))

  dplyr::bind_rows(data, new_rows) |>
    dplyr::arrange(file_name, population_full_path)
}
