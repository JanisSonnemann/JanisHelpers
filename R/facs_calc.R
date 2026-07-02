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

resolve_var_ <- function(df, var, arg_name) {
  if (is.character(var) && length(var) == 1L) {
    if (!var %in% names(df)) {
      stop(glue::glue("Column '{var}' (from `{arg_name}`) not found in `meta`."))
    }
    df[[var]]
  } else if (is.numeric(var) && length(var) == 1L) {
    rep(var, nrow(df))
  } else {
    stop(glue::glue(
      "`{arg_name}` must be a single column name (character) or a single numeric value."
    ))
  }
}

#' Compute absolute cell counts per gram of tissue
#'
#' @description
#' Converts raw event counts into cells per gram of processed tissue, using
#' per-mouse, per-tissue organ weights and staining volumes from
#' \code{meta}. Every \code{mouse_ID} x \code{tissue} combination present in
#' \code{data} is processed in one call (joined against \code{meta} on
#' \code{mouse_ID} and \code{tissue}). Supports two counting methods:
#' HTS/volumetric (the default for every sample unless \code{method_col}
#' says otherwise) and bead-based (using a reference bead population's
#' count and a known bead concentration).
#'
#' @param data tibble shaped like \code{facs_read_wsp(...)$data}, must
#'   include \code{mouse_ID} and \code{tissue} columns (e.g. joined via
#'   \code{facs_read_wsp(keywords = c("mouse_ID", "tissue"))}).
#' @param meta tibble of per-mouse, per-tissue metadata with \code{mouse_ID}
#'   and \code{tissue} columns, e.g. from \code{meta_read()}'s sheets
#'   combined via \code{meta_clean()}.
#' @param vol_total column name in \code{meta} (character) or a single
#'   numeric constant: total organ digest volume.
#' @param vol_stained column name or constant: volume of digest taken for
#'   staining.
#' @param vol_resuspended column name or constant: volume the stained pellet
#'   was resuspended in.
#' @param vol_measured column name or constant: volume actually run/measured
#'   on the cytometer.
#' @param organ_piece_weight column name or constant: weight (mg) of the
#'   organ piece processed.
#' @param method_col character; column name found in \code{data} (checked
#'   first, per-sample) or \code{meta} (checked second, per-mouse) with
#'   values \code{"beads"} or \code{"hts"}. \code{NA} defaults to
#'   \code{"hts"}. \code{NULL} (default) applies \code{"hts"} to every
#'   sample.
#' @param bead_pop character; leaf population name used to look up bead
#'   counts, default \code{"beads"}.
#' @param bead_concentration numeric; reference bead concentration
#'   (beads per microliter), default \code{10400}.
#'
#' @returns \code{data} with additional rows appended:
#'   \code{metric = "count_per_g"}. Errors if \code{method_col} is not found
#'   in \code{data} or \code{meta}, or contains a value outside
#'   \code{{"beads", "hts", NA}}. Warns and fills \code{NA} if a
#'   \code{mouse_ID}/\code{tissue} combination in \code{data} has no match
#'   in \code{meta}, or if the bead method is resolved for a sample with no
#'   matching bead count.
#' @export
#'
#' @examples
#' \dontrun{
#'   dat <- facs_read_wsp("experiment.wsp", keywords = c("mouse_ID", "tissue"))$data
#'   meta_combined <- meta_read("meta.xlsx") |> meta_clean()
#'   facs_calc_count_per_g(
#'     dat, meta_combined,
#'     vol_total = "total_vol", vol_stained = "overview_vol",
#'     vol_resuspended = "overview_resuspended_vol", vol_measured = "overview_measured_vol",
#'     organ_piece_weight = "facs_weight"
#'   )
#' }
facs_calc_count_per_g <- function(
    data,
    meta,
    vol_total,
    vol_stained,
    vol_resuspended,
    vol_measured,
    organ_piece_weight,
    method_col = NULL,
    bead_pop = "beads",
    bead_concentration = 10400
) {
  m <- meta |>
    dplyr::mutate(
      vol_total          = resolve_var_(meta, .env$vol_total, "vol_total"),
      vol_stained        = resolve_var_(meta, .env$vol_stained, "vol_stained"),
      vol_resuspended    = resolve_var_(meta, .env$vol_resuspended, "vol_resuspended"),
      vol_measured       = resolve_var_(meta, .env$vol_measured, "vol_measured"),
      organ_piece_weight = resolve_var_(meta, .env$organ_piece_weight, "organ_piece_weight")
    ) |>
    dplyr::select(mouse_ID, tissue, vol_total, vol_stained, vol_resuspended, vol_measured, organ_piece_weight)

  if (!is.null(method_col) && !method_col %in% names(data) && method_col %in% names(meta)) {
    m[[method_col]] <- meta[[method_col]]
  }
  if (!is.null(method_col) && !method_col %in% names(data) && !method_col %in% names(meta)) {
    stop(glue::glue("`method_col` ('{method_col}') not found in `data` or `meta`."))
  }

  beads <- data |>
    dplyr::filter(population == bead_pop, metric == "count") |>
    dplyr::select(file_name, bead_count = value)

  filtered <- data |>
    dplyr::filter(metric == "count", population != bead_pop) |>
    dplyr::left_join(beads, by = "file_name") |>
    dplyr::left_join(m, by = c("mouse_ID", "tissue"))

  unmatched_combos <- dplyr::anti_join(
    dplyr::distinct(filtered, mouse_ID, tissue),
    dplyr::distinct(m, mouse_ID, tissue),
    by = c("mouse_ID", "tissue")
  )
  if (nrow(unmatched_combos) > 0L) {
    unmatched_desc <- purrr::pmap_chr(
      unmatched_combos,
      function(mouse_ID, tissue) glue::glue("mouse_ID={mouse_ID}, tissue={tissue}")
    )
    warning(glue::glue(
      "The following mouse_ID/tissue combination(s) in `data` have no match in `meta`: ",
      "{paste(unmatched_desc, collapse = '; ')}"
    ))
  }

  filtered$method <- if (is.null(method_col)) {
    "hts"
  } else {
    dplyr::coalesce(filtered[[method_col]], "hts")
  }

  bad_methods <- unique(filtered$method[!filtered$method %in% c("beads", "hts")])
  if (length(bad_methods) > 0L) {
    stop(glue::glue(
      "`method_col` contains value(s) other than 'beads'/'hts': ",
      "{paste(bad_methods, collapse = ', ')}"
    ))
  }

  missing_bead <- filtered$method == "beads" & is.na(filtered$bead_count)
  if (any(missing_bead)) {
    warning(glue::glue(
      "Bead method resolved but no bead count found for file_name(s): ",
      "{paste(unique(filtered$file_name[missing_bead]), collapse = ', ')}. Result filled with NA."
    ))
  }

  new_rows <- filtered |>
    dplyr::mutate(
      metric = "count_per_g",
      value  = dplyr::if_else(
        method == "beads",
        ((value / (bead_count / bead_concentration)) / (vol_stained / vol_total)) / (organ_piece_weight / 1000),
        ((value / (vol_measured / vol_resuspended)) / (vol_stained / vol_total)) / (organ_piece_weight / 1000)
      )
    ) |>
    dplyr::select(dplyr::all_of(names(data)))

  dplyr::bind_rows(data, new_rows) |>
    dplyr::arrange(file_name, population_full_path)
}
