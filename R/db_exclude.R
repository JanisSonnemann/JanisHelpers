# Internal: match exclusion-sheet rows to the samples of an ELISA / histology
# import. `data` has `mouse_ID`, `tissue`, `collected_at` (Date; NA = the
# subject's end date). Rows apply when `panel` is empty or equals `assay`;
# non-empty `tissue` / `collected_at` must equal the sample's (an exclusion
# date never matches an NA sample date). `condition` / `file_name` are FACS
# concepts and error on applicable rows.
# Returns list(applicable, excluded_samples, unmatched).
match_sample_exclusions_ <- function(data, exclusions, assay) {
  empty_applicable <- tibble::tibble(
    mouse_ID = character(), tissue = character(), panel = character(),
    collected_at = as.Date(character()), reason = character(),
    date = as.Date(character()), n_samples_matched = integer()
  )
  empty_samples <- tibble::tibble(
    mouse_ID = character(), tissue = character(),
    collected_at = as.Date(character()), reason = character()
  )
  if (is.null(exclusions) || nrow(exclusions) == 0) {
    return(list(applicable = empty_applicable, excluded_samples = empty_samples,
                unmatched = empty_applicable))
  }

  missing_cols <- setdiff(c("mouse_ID", "reason"), names(exclusions))
  if (length(missing_cols) > 0) {
    stop(
      "The exclusions sheet needs column(s): ",
      paste(missing_cols, collapse = ", "), ".",
      call. = FALSE
    )
  }
  for (col in c("tissue", "panel", "condition", "file_name", "collected_at", "date")) {
    if (!col %in% names(exclusions)) exclusions[[col]] <- NA
  }
  ex <- exclusions |>
    dplyr::transmute(
      dplyr::across(
        dplyr::all_of(c("mouse_ID", "tissue", "panel", "condition", "file_name", "reason")),
        as.character
      ),
      collected_at = as.Date(collected_at),
      date = as.Date(date)
    ) |>
    dplyr::filter(is.na(panel) | panel == "" | panel == !!assay)

  bad <- is.na(ex$mouse_ID) | ex$mouse_ID == "" | is.na(ex$reason) | ex$reason == ""
  if (any(bad)) {
    stop(
      "Every exclusion row needs a mouse_ID and a reason; missing in row(s) ",
      "for: ", paste(ifelse(is.na(ex$mouse_ID[bad]), "<no mouse_ID>",
                            ex$mouse_ID[bad]), collapse = ", "), ".",
      call. = FALSE
    )
  }
  filled <- function(x) !is.na(x) & x != ""
  not_applicable <- filled(ex$condition) | filled(ex$file_name)
  if (any(not_applicable)) {
    stop(
      "condition / file_name don't apply to ", assay, " -- clear them in the ",
      "exclusions sheet (row(s) for: ",
      paste(ex$mouse_ID[not_applicable], collapse = ", "), ").",
      call. = FALSE
    )
  }
  ex <- dplyr::select(ex, !c(condition, file_name))

  samples <- dplyr::distinct(data, mouse_ID, tissue, collected_at) |>
    dplyr::mutate(
      dplyr::across(c(mouse_ID, tissue), as.character),
      collected_at = as.Date(collected_at)
    )

  hits <- purrr::map(seq_len(nrow(ex)), function(i) {
    keep <- !is.na(samples$mouse_ID) & samples$mouse_ID == ex$mouse_ID[[i]]
    if (filled(ex$tissue[[i]])) {
      keep <- keep & !is.na(samples$tissue) & samples$tissue == ex$tissue[[i]]
    }
    if (!is.na(ex$collected_at[[i]])) {
      keep <- keep & !is.na(samples$collected_at) &
        samples$collected_at == ex$collected_at[[i]]
    }
    which(keep)
  })

  applicable <- dplyr::mutate(ex, n_samples_matched = lengths(hits))
  excluded <- purrr::map(seq_along(hits), function(i) {
    dplyr::mutate(samples[hits[[i]], ], reason = ex$reason[[i]])
  }) |>
    dplyr::bind_rows(empty_samples) |>
    dplyr::distinct()

  list(
    applicable = applicable,
    excluded_samples = excluded,
    unmatched = dplyr::filter(applicable, n_samples_matched == 0L)
  )
}
