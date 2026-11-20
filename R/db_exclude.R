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

# ---- Versioned exclusion rule sets -------------------------------------------

#' Store an experiment's exclusion rules as a new versioned set
#'
#' Exclusions live in the database as one rule set per experiment. The
#' current set is the most recent one. A change of mind is editing the
#' exclusions sheet of the meta workbook and calling this function again; no
#' prep, QC report or export is re-run. Every set is kept in `exclusion_sets`
#' (with `written_at`, `package_version`, `source_file` and `note`), so the
#' history of decisions stays visible.
#'
#' A new set is only written when the rules differ from the current one
#' (blank strings count as `NA`, row order is ignored); otherwise a message
#' says "unchanged". An empty or `NULL` sheet writes an empty set, which
#' re-includes everything. Rules naming a panel that is not registered yet, or
#' matching no stored row, are reported as messages, not errors, because the
#' data of that unit may not be imported yet.
#'
#' @param con A DuckDB connection from [db_connect()].
#' @param exclusions The exclusions sheet (from [meta_read()]) with columns
#'   `mouse_ID` and `reason`, and optional `tissue`, `panel`, `condition`,
#'   `file_name`, `collected_at`, `date` (the decision date). Empty or absent
#'   scope columns mean "all". `NULL` or zero rows give an empty set.
#' @param experiment_code Registered experiment the rules belong to. Every
#'   `mouse_ID` must be a registered subject of it.
#' @param source_file Optional path of the workbook, relative to the
#'   project's `data/` folder.
#' @param note Optional free text stored with the set.
#'
#' @returns The `exclusion_set_id` (integer), invisibly -- the new one, or the
#'   current one when nothing changed.
#' @export
db_write_exclusions <- function(con, exclusions, experiment_code,
                                source_file = NA_character_,
                                note = NA_character_) {
  exp_id <- DBI::dbGetQuery(
    con, "SELECT experiment_id FROM experiments WHERE experiment_code = ?",
    params = list(experiment_code)
  )$experiment_id
  if (length(exp_id) == 0) {
    stop("Unknown experiment_code: ", experiment_code, ". Register it first.",
         call. = FALSE)
  }

  rules <- normalise_exclusions_(exclusions)

  subjects <- DBI::dbGetQuery(
    con, "SELECT mouse_id FROM subjects WHERE experiment_id = ?",
    params = list(exp_id)
  )$mouse_id
  unknown <- setdiff(rules$mouse_id, subjects)
  if (length(unknown) > 0) {
    stop(
      "mouse_ID not registered as a subject of ", experiment_code, ": ",
      paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }

  current_id <- current_exclusion_set_id_(con, exp_id)
  if (!is.na(current_id)) {
    current <- current_exclusion_rules_(con, experiment_code) |>
      dplyr::select(!exclusion_set_id)
    if (identical(rules, current)) {
      message("Exclusion rules of ", experiment_code, " unchanged (set ",
              current_id, ").")
      return(invisible(current_id))
    }
  }

  known_panels <- DBI::dbGetQuery(con, "SELECT DISTINCT assay_name FROM assays")$assay_name
  unknown_panels <- setdiff(stats::na.omit(rules$panel), known_panels)
  if (length(unknown_panels) > 0) {
    message(
      "Panel(s) not registered yet: ", paste(unknown_panels, collapse = ", "),
      ". The rules are stored and apply once that panel is exported."
    )
  }

  unmatched <- rules_without_match_(con, rules)
  if (length(unmatched) > 0) {
    message(
      "Rule(s) for mouse/mice ", paste(unmatched, collapse = ", "),
      " match no stored row (yet)."
    )
  }

  version <- as.character(utils::packageVersion("JanisHelpers"))
  new_id <- NA_integer_
  DBI::dbWithTransaction(con, {
    new_id <- DBI::dbGetQuery(con, "
      INSERT INTO exclusion_sets (experiment_id, source_file, package_version, note)
      VALUES (?, ?, ?, ?) RETURNING exclusion_set_id
    ", params = list(exp_id, as.character(source_file), version,
                     as.character(note)))$exclusion_set_id
    if (nrow(rules) > 0) {
      tmp <- dplyr::mutate(rules, exclusion_set_id = new_id)
      duckdb::duckdb_register(con, "tmp_exclusion_rules", tmp)
      on.exit(duckdb::duckdb_unregister(con, "tmp_exclusion_rules"), add = TRUE)
      DBI::dbExecute(con, "
        INSERT INTO exclusion_rules (exclusion_set_id, mouse_id, tissue, panel,
          condition, file_name, collected_at, reason, decided_on)
        SELECT exclusion_set_id, mouse_id, tissue, panel, condition, file_name,
          collected_at, reason, decided_on
        FROM tmp_exclusion_rules
      ")
    }
  })
  invisible(as.integer(new_id))
}

# Internal: id of an experiment's current (highest) set, NA if none.
current_exclusion_set_id_ <- function(con, experiment_id) {
  id <- DBI::dbGetQuery(
    con, "SELECT MAX(exclusion_set_id) AS id FROM exclusion_sets WHERE experiment_id = ?",
    params = list(experiment_id)
  )$id
  as.integer(id)
}

# Internal: the rules of an experiment's current set (0 rows when no set).
current_exclusion_rules_ <- function(con, experiment_code) {
  out <- DBI::dbGetQuery(con, "
    SELECT r.exclusion_set_id, r.mouse_id, r.tissue, r.panel, r.condition,
           r.file_name, r.collected_at, r.reason, r.decided_on
    FROM exclusion_rules r
    JOIN exclusion_sets s USING (exclusion_set_id)
    JOIN experiments e USING (experiment_id)
    WHERE e.experiment_code = ?
      AND s.exclusion_set_id = (
        SELECT MAX(s2.exclusion_set_id) FROM exclusion_sets s2
        WHERE s2.experiment_id = s.experiment_id)
  ", params = list(experiment_code))
  out <- tibble::as_tibble(out)
  out$collected_at <- as.Date(out$collected_at)
  out$decided_on <- as.Date(out$decided_on)
  sort_exclusions_(out)
}

sort_exclusions_ <- function(x) {
  scope <- c("tissue", "panel", "condition", "file_name", "collected_at")
  x <- tibble::as_tibble(x)
  dplyr::arrange(x, dplyr::across(dplyr::all_of(c("mouse_id", scope, "reason"))))
}

# Internal: exclusions sheet -> rule columns (validated, sorted). `""` -> NA;
# absent scope columns are all-NA; extra columns are ignored.
normalise_exclusions_ <- function(exclusions) {
  empty <- tibble::tibble(
    mouse_id = character(), tissue = character(), panel = character(),
    condition = character(), file_name = character(),
    collected_at = as.Date(character()), reason = character(),
    decided_on = as.Date(character())
  )
  if (is.null(exclusions)) return(empty)
  missing_cols <- setdiff(c("mouse_ID", "reason"), names(exclusions))
  if (length(missing_cols) > 0) {
    stop("The exclusions sheet needs column(s): ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  if (nrow(exclusions) == 0) return(empty)

  blank_na <- function(x) {
    x <- as.character(x)
    x[!is.na(x) & trimws(x) == ""] <- NA_character_
    x
  }
  get_col <- function(name) {
    if (name %in% names(exclusions)) exclusions[[name]] else rep(NA, nrow(exclusions))
  }
  out <- tibble::tibble(
    mouse_id = blank_na(exclusions$mouse_ID),
    tissue = blank_na(get_col("tissue")),
    panel = blank_na(get_col("panel")),
    condition = blank_na(get_col("condition")),
    file_name = blank_na(get_col("file_name")),
    collected_at = as.Date(get_col("collected_at")),
    reason = blank_na(exclusions$reason),
    decided_on = as.Date(get_col("date"))
  )
  bad <- is.na(out$mouse_id) | is.na(out$reason)
  if (any(bad)) {
    stop(
      "Every exclusion row needs a mouse_ID and a reason; missing in row(s) ",
      "for: ", paste(ifelse(is.na(out$mouse_id[bad]), "<no mouse_ID>",
                            out$mouse_id[bad]), collapse = ", "), ".",
      call. = FALSE
    )
  }
  sort_exclusions_(out)
}

# Internal: mouse ids of `rules` that match no stored row. Simple version: the
# mouse has no row in `samples`. (Task 5 swaps this for a view-based match.)
rules_without_match_ <- function(con, rules) {
  if (nrow(rules) == 0) return(character())
  have <- DBI::dbGetQuery(con, "
    SELECT DISTINCT s.mouse_id FROM samples sm
    JOIN subjects s USING (subject_id)
  ")$mouse_id
  setdiff(unique(rules$mouse_id), have)
}
