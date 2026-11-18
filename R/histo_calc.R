#' Add a ratio metric to long histology data
#'
#' @description
#' Computes \code{numerator / denominator * scale} per
#' \code{mouse_ID x tissue x collected_at} and appends it as new rows with
#' \code{metric = name}. The result is \code{NA} when either value is
#' \code{NA} or the denominator is 0.
#'
#' @param data long tibble from \code{histo_read_curated()}
#' @param numerator,denominator character; metric names in \code{data}
#' @param name character; metric name of the new rows (must not exist yet)
#' @param scale numeric; multiplier, default 100 (percent)
#'
#' @returns \code{data} with the extra rows appended
#' @export
#'
#' @examples
#' \dontrun{
#'   histo_read_curated("kidney_26-1.xlsx", tissue = "kidney") |>
#'     histo_calc_ratio("crescents", "gloms", "pct_crescents")
#' }
histo_calc_ratio <- function(data, numerator, denominator, name, scale = 100) {
  if (name %in% data$metric) {
    stop(glue::glue("Metric '{name}' already exists in `data`."), call. = FALSE)
  }
  absent <- setdiff(c(numerator, denominator), data$metric)
  if (length(absent) > 0L) {
    stop(glue::glue("Metric(s) not in `data`: {paste(absent, collapse = ', ')}."),
         call. = FALSE)
  }
  pick <- function(m, col) {
    data |>
      dplyr::filter(.data$metric == m) |>
      dplyr::select("mouse_ID", "tissue", "collected_at", !!col := "value")
  }
  new <- pick(numerator, "num") |>
    dplyr::left_join(pick(denominator, "den"),
                     by = c("mouse_ID", "tissue", "collected_at")) |>
    dplyr::mutate(
      metric = name,
      value = dplyr::if_else(is.na(.data$den) | .data$den == 0, NA_real_,
                             .data$num / .data$den * scale)
    ) |>
    dplyr::select(dplyr::all_of(names(data)))
  dplyr::bind_rows(data, new)
}
