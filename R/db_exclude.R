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
#' says "unchanged". Writing a set always gives a message with the new set id
#' and its number of rules (and the previous set's). An empty or `NULL` sheet
#' writes an empty set, which re-includes everything; replacing a non-empty
#' set by an empty one also gives a warning. (The prep runners, e.g.
#' [facs_prep_run()], skip this call when their `exclusions` is `NULL`, so a
#' final run without the sheet leaves the rule set unchanged.) Rules naming a panel that is not registered yet, or
#' matching no stored row, are reported as messages, not errors, because the
#' data of that unit may not be imported yet. `file_name` rules can't match
#' FACS stains exported before file names were stored; the message then asks
#' for a re-export of that panel.
#'
#' The rules are applied on read: [db_query_facs()], [db_query_elisa()] and
#' [db_query_histo()] flag matching rows with `excluded` / `exclusion_reason`,
#' and [db_query_subjects()] flags mouse-wide rules.
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
  n_before <- 0L
  if (!is.na(current_id)) {
    current <- current_exclusion_rules_(con, experiment_code) |>
      dplyr::select(!exclusion_set_id)
    if (identical(rules, current)) {
      message("Exclusion rules of ", experiment_code, " unchanged (set ",
              current_id, ").")
      return(invisible(current_id))
    }
    n_before <- nrow(current)
  }

  known_panels <- DBI::dbGetQuery(con, "SELECT DISTINCT assay_name FROM assays")$assay_name
  unknown_panels <- setdiff(stats::na.omit(rules$panel), known_panels)
  if (length(unknown_panels) > 0) {
    message(
      "Panel(s) not registered yet: ", paste(unknown_panels, collapse = ", "),
      ". The rules are stored and apply once that panel is exported."
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
  message("Exclusion rules of ", experiment_code, ": set ", new_id, " written, ",
          nrow(rules), " rule(s) (was ", n_before, ").")
  if (nrow(rules) == 0 && n_before > 0) {
    warning("Exclusion rules of ", experiment_code, " replaced by an empty set ",
            "-- every previously excluded row is now included.", call. = FALSE)
  }

  unmatched <- rules_without_match_(con, experiment_code)
  show <- function(x) format_rows_(dplyr::select(x, !needs_reexport))
  stale <- dplyr::filter(unmatched, needs_reexport)
  unmatched <- dplyr::filter(unmatched, !needs_reexport)
  if (nrow(unmatched) > 0) {
    message("These rules match no stored row (yet):\n", show(unmatched))
  }
  if (nrow(stale) > 0) {
    message(
      "These file_name rules can't match yet: the matching stains were ",
      "exported before file names were stored -- re-export the panel ",
      "(final = TRUE).\n", show(stale)
    )
  }
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

# Internal: the rules of the experiment's current set that match no stored
# FACS / ELISA / histology row (same match as the *_exclusion_flags views).
# `needs_reexport` is TRUE for file_name rules that would match a FACS stain
# stored without a file name (exported before file names were kept).
# Returns a tibble with the rule columns plus `needs_reexport` (0 rows when
# there is no set or every rule matches).
rules_without_match_ <- function(con, experiment_code) {
  matched_sql <- function(...) {
    paste("SELECT DISTINCT r.exclusion_rule_id", exclusion_match_from_(...))
  }
  sql <- paste0("
    WITH matched AS (",
      matched_sql("facs"), " UNION ", matched_sql("elisa"), " UNION ",
      matched_sql("histo"), "),
    unnamed AS (",
      matched_sql("facs", "r.file_name IS NOT NULL AND fs.file_name IS NULL"), ")
    SELECT r.mouse_id, r.tissue, r.panel, r.condition, r.file_name,
           r.collected_at, r.reason,
           r.exclusion_rule_id IN (SELECT exclusion_rule_id FROM unnamed)
             AS needs_reexport
    FROM exclusion_rules r
    JOIN exclusion_sets s USING (exclusion_set_id)
    JOIN experiments e USING (experiment_id)
    WHERE e.experiment_code = ?
      AND s.exclusion_set_id = (
        SELECT MAX(s2.exclusion_set_id) FROM exclusion_sets s2
        WHERE s2.experiment_id = s.experiment_id)
      AND r.exclusion_rule_id NOT IN (SELECT exclusion_rule_id FROM matched)
    ORDER BY r.mouse_id, r.tissue, r.panel, r.condition, r.file_name,
             r.collected_at, r.reason
  ")
  out <- tibble::as_tibble(DBI::dbGetQuery(con, sql, params = list(experiment_code)))
  out$collected_at <- as.Date(out$collected_at)
  out
}

# Internal: the single experiment code the prep's mice belong to. Mice not
# registered at all are left to the export's own error.
experiment_code_of_prep_ <- function(con, prep) {
  mice <- unique(stats::na.omit(as.character(prep$data$mouse_ID)))
  found <- if (length(mice) == 0) {
    character()
  } else {
    DBI::dbGetQuery(
      con,
      paste0("SELECT DISTINCT e.experiment_code FROM subjects s ",
             "JOIN experiments e ON e.experiment_id = s.experiment_id ",
             "WHERE s.mouse_id IN (",
             paste(rep("?", length(mice)), collapse = ", "), ")"),
      params = as.list(mice)
    )$experiment_code
  }
  if (length(found) == 0) {
    stop("None of the prep's mice is registered in the database; register ",
         "the experiment first (db_write_meta()).", call. = FALSE)
  }
  if (length(found) > 1) {
    stop("The prep's mice span several experiments: ",
         paste(sort(found), collapse = ", "), ".", call. = FALSE)
  }
  found
}

# Internal: record the exclusions sheet as the experiment's rule set. A NULL
# sheet (not given to the runner) leaves the stored set alone -- only an
# explicit (possibly zero-row) sheet replaces it. Returns the set id, or NULL
# (invisibly) when nothing was written.
write_prep_exclusions_ <- function(con, prep, exclusions, meta_path, data_dir) {
  experiment_code <- experiment_code_of_prep_(con, prep)
  if (is.null(exclusions)) {
    message("exclusions not given -- rule set of ", experiment_code,
            " left unchanged.")
    return(invisible(NULL))
  }
  source_file <- if (is.null(meta_path)) {
    NA_character_
  } else {
    relative_to_data_dir_(meta_path, data_dir)
  }
  db_write_exclusions(con, exclusions, experiment_code,
                      source_file = source_file)
}
