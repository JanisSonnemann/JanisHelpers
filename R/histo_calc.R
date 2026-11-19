#' Add ratio metrics to long histology data
#'
#' @description
#' Computes \code{numerator / denominator * scale} per
#' \code{mouse_ID x tissue x collected_at} and appends it as new rows. The
#' result is \code{NA} when either value is \code{NA} or the denominator is 0.
#' Several ratios can be added in one call by passing a named
#' \code{numerator} vector; the names become the new metric names. Other
#' columns of \code{data} (e.g. \code{comment}, group columns) are copied
#' from the numerator rows.
#'
#' @param data long tibble from \code{histo_read_curated()}
#' @param numerator character; metric name(s) in \code{data}. Either one
#'   metric (with \code{name}) or a named vector, one element per new ratio
#'   (names = new metric names; \code{name} must then be \code{NULL}).
#' @param denominator character; metric name(s) in \code{data}, length 1
#'   (used for every ratio) or one per numerator
#' @param name character; metric name of the new rows when \code{numerator}
#'   is a single unnamed metric. Default \code{NULL}.
#' @param scale numeric; multiplier, default 100 (percent)
#'
#' @returns \code{data} with the extra rows appended
#' @export
#'
#' @examples
#' \dontrun{
#'   histo_read_curated("kidney_26-1.xlsx", tissue = "kidney") |>
#'     histo_calc_ratio(
#'       c(pct_cresc = "crescents", pct_necrosis = "necrosis",
#'         pct_abnormal = "abnormal_gloms"),
#'       denominator = "gloms"
#'     )
#' }
histo_calc_ratio <- function(data, numerator, denominator, name = NULL,
                             scale = 100) {
  if (!is.null(name)) {
    if (length(numerator) != 1L || length(name) != 1L) {
      stop("`name` is for a single ratio -- for several, pass a named ",
           "`numerator` vector instead.", call. = FALSE)
    }
    numerator <- stats::setNames(numerator, name)
  }
  new_names <- names(numerator)
  if (is.null(new_names) || any(is.na(new_names) | !nzchar(new_names))) {
    stop("`numerator` must be named (names = new metric names), or give ",
         "`name` for a single ratio.", call. = FALSE)
  }
  if (anyDuplicated(new_names)) {
    stop("New metric names must be unique: ",
         paste(unique(new_names[duplicated(new_names)]), collapse = ", "), ".",
         call. = FALSE)
  }
  if (!length(denominator) %in% c(1L, length(numerator))) {
    stop("`denominator` must have length 1 or one per numerator (",
         length(numerator), ").", call. = FALSE)
  }
  denominator <- rep_len(denominator, length(numerator))
  taken <- intersect(new_names, data$metric)
  if (length(taken) > 0L) {
    stop(glue::glue("Metric(s) already exists in `data`: {paste(taken, collapse = ', ')}."),
         call. = FALSE)
  }
  absent <- setdiff(c(numerator, denominator), data$metric)
  if (length(absent) > 0L) {
    stop(glue::glue("Metric(s) not in `data`: {paste(absent, collapse = ', ')}."),
         call. = FALSE)
  }

  keys <- c("mouse_ID", "tissue", "collected_at")
  new <- purrr::map(seq_along(numerator), function(i) {
    num <- data |>
      dplyr::filter(.data$metric == numerator[[i]]) |>
      dplyr::rename(num = "value")
    den <- data |>
      dplyr::filter(.data$metric == denominator[[i]]) |>
      dplyr::select(dplyr::all_of(keys), den = "value")
    num |>
      dplyr::left_join(den, by = keys) |>
      dplyr::mutate(
        metric = new_names[[i]],
        value = dplyr::if_else(is.na(.data$den) | .data$den == 0, NA_real_,
                               .data$num / .data$den * scale)
      ) |>
      dplyr::select(dplyr::all_of(names(data)))
  }) |>
    purrr::list_rbind()
  dplyr::bind_rows(data, new)
}
