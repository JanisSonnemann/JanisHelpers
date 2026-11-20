#' Label samples by week after BMT
#'
#' Adds the study week of each sample, counted from the mouse's bone marrow
#' transplantation (or another reference date). Every sample gets the nearest
#' whole week (days / 7 rounded, halves up; negative before BMT). End-point samples
#' (taken on the mouse's `end_date`) of one experiment harvested in the same
#' calendar week (Mon-Sun) share one week -- the nearest week of their median
#' day -- so one harvest is never split across two labels. An end-point in a
#' different calendar week (an early death) keeps its own week.
#'
#' @param data Long data with the `by` columns and `date_col`, e.g. an ELISA
#'   tibble from [db_load_experiment()].
#' @param meta One row per mouse with the `by` columns, `ref_col` and
#'   `end_date` (e.g. `dat$meta` from [db_load_experiment()]).
#' @param by Key columns joining `data` to `meta`.
#' @param ref_col Column of `meta` with the reference date (default
#'   `"bmt_date"`).
#' @param date_col Column of `data` with the sampling date (default
#'   `"collected_at"`).
#' @returns `data` with three added columns: `days_post_bmt` (integer days
#'   since the reference date), `week_num` (integer week) and `week` (ordered
#'   factor `"week <n>"`, levels in numeric order). Mice without a reference
#'   date get `NA` with a warning.
#' @examples
#' meta <- data.frame(mouse_id = c("m1", "m2"), bmt_date = as.Date("2026-05-12"),
#'                    end_date = as.Date(c("2026-06-15", "2026-06-16")))
#' data <- data.frame(mouse_id = c("m1", "m1", "m2"),
#'                    collected_at = as.Date(c("2026-05-19", "2026-06-15", "2026-06-16")))
#' meta_label_weeks(data, meta, by = "mouse_id")
#' @export
meta_label_weeks <- function(data, meta, by = c("experiment_code", "mouse_id"),
                             ref_col = "bmt_date", date_col = "collected_at") {
  out_cols <- c("days_post_bmt", "week_num", "week")
  clash <- intersect(out_cols, names(data))
  if (length(clash) > 0) {
    stop("`data` already has column(s) ", paste(clash, collapse = ", "), ".", call. = FALSE)
  }
  ref <- meta |>
    dplyr::select(dplyr::all_of(c(by, ref_col, "end_date"))) |>
    dplyr::distinct() |>
    dplyr::mutate(.ref_ = as.Date(.data[[ref_col]]), .end_ = as.Date(end_date)) |>
    dplyr::select(dplyr::all_of(by), ".ref_", ".end_")

  missing <- dplyr::anti_join(dplyr::distinct(data, dplyr::across(dplyr::all_of(by))),
                              ref, by = by)
  if (nrow(missing) > 0) {
    stop("Mice in `data` not found in `meta`: ",
         paste(do.call(paste, c(missing, sep = "/")), collapse = ", "), call. = FALSE)
  }

  d <- dplyr::left_join(data, ref, by = by) |>
    dplyr::mutate(
      .date_ = as.Date(.data[[date_col]]),
      days_post_bmt = as.integer(.date_ - .ref_),
      .end_point_ = !is.na(.end_) & .date_ == .end_,
      week_num = round_half_up_(days_post_bmt / 7)
    )

  no_ref <- unique(d[[by[length(by)]]][is.na(d$.ref_)])
  if (length(no_ref) > 0) {
    warning("No ", ref_col, " for mouse(s) ", paste(no_ref, collapse = ", "),
            ": week is NA.", call. = FALSE)
  }

  pool_by <- intersect("experiment_code", by)
  d <- d |>
    dplyr::mutate(.iso_ = format(.date_, "%G-%V")) |>
    dplyr::mutate(
      .pooled_ = round_half_up_(stats::median(days_post_bmt, na.rm = TRUE) / 7),
      .by = dplyr::all_of(c(pool_by, ".iso_", ".end_point_"))
    ) |>
    dplyr::mutate(week_num = dplyr::if_else(.end_point_ & !is.na(days_post_bmt),
                                            .pooled_, week_num))

  lvls <- sort(unique(stats::na.omit(d$week_num)))
  d |>
    dplyr::mutate(week = factor(dplyr::if_else(is.na(week_num), NA_character_,
                                               paste("week", week_num)),
                                levels = paste("week", lvls), ordered = TRUE)) |>
    dplyr::select(!dplyr::any_of(c(".ref_", ".end_", ".date_", ".end_point_", ".iso_", ".pooled_")))
}

# Internal: round to integer with halves up (base round() rounds halves to
# even, which would put 2.5 and 3.5 weeks on different sides).
round_half_up_ <- function(x) as.integer(floor(x + 0.5))
