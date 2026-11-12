# Errors unless `x` is a single non-NA, non-empty string. `NULL` stands for a
# missing argument. `what` describes the expected value in the message.
check_string_arg_ <- function(x, arg, what) {
  if (is.null(x) || !is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x)) {
    stop("`", arg, "` must be a single non-empty string: ", what, ".", call. = FALSE)
  }
  invisible(x)
}

# Guards for db_write_facs(): the .wsp filename belongs to the root, one root
# per experiment x panel (and vice versa), one workspace version per root.
# Errors before anything is written.
check_facs_root_ <- function(con, sample_id, assay_id, wsp_root, source_file) {
  file <- basename(source_file)
  suffix <- substring(file, nchar(wsp_root) + 1)
  if (!startsWith(file, wsp_root) || !grepl("^_[0-9]{8}\\.wsp$", suffix)) {
    stop(
      "`source_file` '", file, "' doesn't belong to `wsp_root` '", wsp_root,
      "' -- expected '", wsp_root, "_YYYYMMDD.wsp'.",
      call. = FALSE
    )
  }

  experiment_id <- DBI::dbGetQuery(con, "
    SELECT sub.experiment_id FROM samples sm
    JOIN subjects sub ON sub.subject_id = sm.subject_id
    WHERE sm.sample_id = ?", params = list(sample_id))$experiment_id[[1]]

  other_roots <- DBI::dbGetQuery(con, "
    SELECT DISTINCT fs.wsp_root FROM facs_stains fs
    JOIN samples sm ON sm.sample_id = fs.sample_id
    JOIN subjects sub ON sub.subject_id = sm.subject_id
    WHERE sub.experiment_id = ? AND fs.assay_id = ? AND fs.wsp_root <> ?",
    params = list(experiment_id, assay_id, wsp_root))$wsp_root
  if (length(other_roots) > 0) {
    stop(
      "This experiment x panel is already stored under root '", other_roots[[1]],
      "', not '", wsp_root, "'. Re-export root '", other_roots[[1]],
      "' to replace it, or fix the workspace name.",
      call. = FALSE
    )
  }

  elsewhere <- DBI::dbGetQuery(con, "
    SELECT DISTINCT e.experiment_code, a.assay_name FROM facs_stains fs
    JOIN samples sm ON sm.sample_id = fs.sample_id
    JOIN subjects sub ON sub.subject_id = sm.subject_id
    JOIN experiments e ON e.experiment_id = sub.experiment_id
    JOIN assays a ON a.assay_id = fs.assay_id
    WHERE fs.wsp_root = ? AND (sub.experiment_id <> ? OR fs.assay_id <> ?)",
    params = list(wsp_root, experiment_id, assay_id))
  if (nrow(elsewhere) > 0) {
    stop(
      "Root '", wsp_root, "' is already used for experiment '",
      elsewhere$experiment_code[[1]], "' x panel '", elsewhere$assay_name[[1]],
      "'. One root covers exactly one experiment x panel.",
      call. = FALSE
    )
  }

  other_versions <- DBI::dbGetQuery(con, "
    SELECT DISTINCT m.source_file FROM facs_measurements m
    JOIN facs_stains fs ON fs.facs_stain_id = m.facs_stain_id
    WHERE fs.wsp_root = ? AND m.source_file <> ?",
    params = list(wsp_root, source_file))$source_file
  if (length(other_versions) > 0) {
    stop(
      "Root '", wsp_root, "' already holds data from '", other_versions[[1]],
      "'. Replacing a root with a newer workspace version is the Stage 3 ",
      "export's job, not db_write_facs().",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

lookup_sample_id_ <- function(con, mouse_id, tissue) {
  result <- DBI::dbGetQuery(
    con,
    "SELECT sm.sample_id
     FROM samples sm
     JOIN subjects sub ON sub.subject_id = sm.subject_id
     WHERE sub.mouse_id = ? AND sm.tissue = ?",
    params = list(mouse_id, tissue)
  )
  if (nrow(result) == 0) {
    stop(
      "No sample found for mouse_id '", mouse_id, "', tissue '", tissue,
      "' -- register the subject (db_write_subjects()) and sample (db_write_samples()) first.",
      call. = FALSE
    )
  }
  result$sample_id[[1]]
}

# Vectorised lookup_sample_id_(): one sample_id per (mouse_id, tissue) pair,
# in input order. Errors naming every unresolved pair (NA counts as unresolved).
lookup_sample_ids_ <- function(con, mouse_id, tissue) {
  keys <- tibble::tibble(
    row_ = seq_along(mouse_id),
    mouse_id = as.character(mouse_id),
    tissue = as.character(tissue)
  )
  duckdb::duckdb_register(con, "tmp_sample_keys", keys)
  on.exit(duckdb::duckdb_unregister(con, "tmp_sample_keys"), add = TRUE)

  resolved <- DBI::dbGetQuery(con, "
    SELECT k.row_, k.mouse_id, k.tissue, sm.sample_id
    FROM tmp_sample_keys k
    LEFT JOIN subjects sub ON sub.mouse_id = k.mouse_id
    LEFT JOIN samples sm ON sm.subject_id = sub.subject_id AND sm.tissue = k.tissue
    ORDER BY k.row_
  ")

  unresolved <- unique(resolved[is.na(resolved$sample_id), c("mouse_id", "tissue")])
  if (nrow(unresolved) > 0) {
    stop(
      "No sample found for mouse_id x tissue: ",
      paste(sprintf("'%s x %s'", unresolved$mouse_id, unresolved$tissue), collapse = ", "),
      " -- register the subjects (db_write_subjects()) and samples (db_write_samples()) first.",
      call. = FALSE
    )
  }
  resolved$sample_id
}

lookup_assay_id_ <- function(con, assay_name, domain) {
  result <- DBI::dbGetQuery(
    con,
    "SELECT a.assay_id
     FROM assays a
     JOIN domains d ON d.domain_id = a.domain_id
     WHERE a.assay_name = ? AND d.domain_name = ?",
    params = list(assay_name, domain)
  )
  if (nrow(result) == 0) {
    stop(
      "No assay found for assay_name '", assay_name, "' in domain '", domain,
      "' -- register it first via db_write_assay().",
      call. = FALSE
    )
  }
  result$assay_id[[1]]
}

#' Import FACS measurements for one stain
#'
#' Writes one staining event's tidy long measurements (as produced by
#' `facs_import_wsp()`, filtered to a single `FileName`) into
#' `facs_measurements`, creating the matching `facs_stains` row (tagged with
#' `wsp_root`) first if it doesn't already exist. Stain volumes, bead settings
#' and organ weights are not stored: they are calculation inputs, used by the
#' export script before the data reaches the database.
#'
#' Errors, writing nothing, if:
#' - `basename(source_file)` isn't `<wsp_root>_YYYYMMDD.wsp`;
#' - this experiment x panel is already stored under a different root;
#' - `wsp_root` is already used for a different experiment or panel;
#' - the root already holds data from a different `source_file` (a newer
#'   workspace version replaces a root via the Stage 3 export, not here).
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param data A tibble with columns `PopulationFullPath`, `Population`,
#'   `metric`, `value` -- the shape `facs_import_wsp()` produces, filtered to
#'   one `FileName`.
#' @param mouse_id Character scalar identifying the subject (already
#'   registered via [db_write_subjects()]).
#' @param tissue Character scalar identifying the sample (already registered
#'   via [db_write_samples()]).
#' @param assay_name Character scalar identifying the panel (already
#'   registered via [db_write_assay()] with `domain = "facs"`).
#' @param wsp_root Character scalar, required -- the workspace root (the
#'   `.wsp` filename without its `_YYYYMMDD` version suffix, e.g.
#'   `"treg_25-7"`). One root = one experiment x one panel.
#' @param source_file Character scalar, required -- the FlowJo workspace the
#'   measurements came from, as a path relative to the lab database's `raw/`
#'   folder (e.g. `"25-7/facs/treg_25-7_20261112.wsp"`). Recorded on every
#'   measurement row as provenance. The FCS filename is deliberately not used:
#'   it doesn't identify which workspace version (gates, compensation)
#'   produced the numbers.
#' @param stain_date Date scalar, optional.
#' @returns Invisibly, the number of measurement rows inserted.
#' @export
db_write_facs <- function(con, data, mouse_id, tissue, assay_name, wsp_root,
                           source_file, stain_date = NA) {
  check_string_arg_(
    if (missing(wsp_root)) NULL else wsp_root, "wsp_root",
    "the .wsp filename without its _YYYYMMDD suffix (e.g. \"treg_25-7\")"
  )
  check_string_arg_(
    if (missing(source_file)) NULL else source_file, "source_file",
    "the .wsp path relative to raw/ (e.g. \"25-7/facs/treg_25-7_20261112.wsp\")"
  )

  sample_id <- lookup_sample_id_(con, mouse_id, tissue)
  assay_id <- lookup_assay_id_(con, assay_name, domain = "facs")
  check_facs_root_(con, sample_id, assay_id, wsp_root, source_file)

  to_write <- tibble::tibble(
    population_full_path = as.character(data$PopulationFullPath),
    population = as.character(data$Population),
    metric = as.character(data$metric),
    value = as.double(data$value),
    source_file = source_file
  )
  duckdb::duckdb_register(con, "tmp_facs_measurements", to_write)
  on.exit(duckdb::duckdb_unregister(con, "tmp_facs_measurements"), add = TRUE)

  n <- DBI::dbWithTransaction(con, {
    DBI::dbExecute(
      con,
      "INSERT INTO facs_stains (sample_id, assay_id, wsp_root, stain_date)
       VALUES (?, ?, ?, ?)
       ON CONFLICT (sample_id, assay_id) DO NOTHING",
      params = list(sample_id, assay_id, wsp_root, stain_date)
    )

    facs_stain_id <- DBI::dbGetQuery(
      con,
      "SELECT facs_stain_id FROM facs_stains WHERE sample_id = ? AND assay_id = ?",
      params = list(sample_id, assay_id)
    )$facs_stain_id[[1]]

    DBI::dbExecute(con, "
      INSERT INTO facs_measurements (
        facs_stain_id, population_full_path, population, metric, value, source_file
      )
      SELECT ?, population_full_path, population, metric, value, source_file
      FROM tmp_facs_measurements
      ON CONFLICT (facs_stain_id, population_full_path, metric) DO NOTHING
    ", params = list(facs_stain_id))
  })

  invisible(n)
}

#' Import ELISA measurements
#'
#' Writes a tidy tibble of ELISA measurements (as produced by
#' `elisa_read_results()`) into `elisa_measurements`.
#'
#' `elisa_measurements`'s uniqueness constraint (which backs duplicate-import
#' detection) is `UNIQUE(sample_id, assay_id, cytokine, sample_id_raw,
#' replicate)`; because SQL treats `NULL <> NULL`, an `NA` in `data$sample_id`
#' (written as `sample_id_raw`) or `data$replicate` would let a repeat import
#' insert duplicate rows undetected, so this function errors before writing
#' anything if either column contains `NA`.
#'
#' A whole plate is written in one call: every row carries its own
#' `mouse_id` and `tissue`, and each row's sample is resolved against the
#' registered `subjects`/`samples`. If any `mouse_id` x `tissue` combination
#' isn't registered (including an `NA` `mouse_id`), the function errors
#' naming every unresolved combination and writes nothing. Drop standards,
#' blanks and controls from `data` before calling.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param data A tibble with columns `cytokine`, `sample_id` (raw plate
#'   label), `replicate`, `value`, `unit`, `result_status` -- the shape
#'   `elisa_read_results()` produces -- plus `mouse_id` and `tissue`, added
#'   beforehand by joining a plate-label -> mouse table (e.g. with
#'   [meta_annotate()]). `sample_id` and `replicate` must not contain `NA`.
#' @param assay_name Character scalar identifying the assay/panel (already
#'   registered via [db_write_assay()] with `domain = "elisa"`).
#' @param source_file Character scalar, optional.
#' @returns Invisibly, the number of rows inserted.
#' @export
db_write_elisa <- function(con, data, assay_name, source_file = NA_character_) {
  missing_cols <- setdiff(c("mouse_id", "tissue"), names(data))
  if (length(missing_cols) > 0) {
    hint <- if ("mouse_id" %in% missing_cols && "mouse_ID" %in% names(data)) {
      " `data` has `mouse_ID` -- rename it to `mouse_id` (the database's spelling)."
    } else {
      ""
    }
    stop(
      "db_write_elisa() needs column(s) ",
      paste(sprintf("`%s`", missing_cols), collapse = ", "),
      " in `data`, one value per row.", hint,
      call. = FALSE
    )
  }
  if (anyNA(data$sample_id) || anyNA(data$replicate)) {
    stop(
      "db_write_elisa() requires non-NA sample_id and replicate for every row ",
      "-- these back the table's uniqueness constraint; NA values would break ",
      "duplicate-import detection.",
      call. = FALSE
    )
  }

  assay_id <- lookup_assay_id_(con, assay_name, domain = "elisa")
  sample_ids <- lookup_sample_ids_(con, data$mouse_id, data$tissue)

  to_write <- tibble::tibble(
    sample_id = sample_ids,
    assay_id = assay_id,
    cytokine = as.character(data$cytokine),
    sample_id_raw = as.character(data$sample_id),
    replicate = as.integer(data$replicate),
    value = as.double(data$value),
    unit = as.character(data$unit),
    result_status = as.character(data$result_status),
    source_file = source_file
  )
  duckdb::duckdb_register(con, "tmp_elisa_measurements", to_write)
  on.exit(duckdb::duckdb_unregister(con, "tmp_elisa_measurements"), add = TRUE)

  n <- DBI::dbExecute(con, "
    INSERT INTO elisa_measurements (
      sample_id, assay_id, cytokine, sample_id_raw, replicate, value, unit, result_status, source_file
    )
    SELECT sample_id, assay_id, cytokine, sample_id_raw, replicate, value, unit, result_status, source_file
    FROM tmp_elisa_measurements
    ON CONFLICT (sample_id, assay_id, cytokine, sample_id_raw, replicate) DO NOTHING
  ")
  invisible(n)
}

#' Import histology measurements
#'
#' Writes a tidy tibble of histology measurements into `histo_measurements`.
#'
#' Note: `assay_name` is optional (histology has traditionally used a single
#' un-named stain/scoring protocol). Rows written with `assay_name = NA` are
#' not deduplicated against each other -- SQL's `NULL <> NULL` means the
#' `UNIQUE(sample_id, assay_id, metric)` constraint doesn't catch a repeated
#' import when `assay_id` is `NULL`. Pass an `assay_name` (registered via
#' [db_write_assay()] with `domain = "histo"`) to get full dedup protection.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param data A tibble with columns `metric`, `value`.
#' @param mouse_id Character scalar identifying the subject.
#' @param tissue Character scalar identifying the sample.
#' @param assay_name Character scalar, optional -- identifies the
#'   stain/scoring protocol.
#' @param source_file Character scalar, optional.
#' @returns Invisibly, the number of rows inserted.
#' @export
db_write_histo <- function(con, data, mouse_id, tissue, assay_name = NA_character_,
                            source_file = NA_character_) {
  sample_id <- lookup_sample_id_(con, mouse_id, tissue)
  assay_id <- if (!is.na(assay_name)) lookup_assay_id_(con, assay_name, domain = "histo") else NA_integer_

  to_write <- tibble::tibble(
    sample_id = sample_id,
    assay_id = assay_id,
    metric = as.character(data$metric),
    value = as.double(data$value),
    source_file = source_file
  )
  duckdb::duckdb_register(con, "tmp_histo_measurements", to_write)
  on.exit(duckdb::duckdb_unregister(con, "tmp_histo_measurements"), add = TRUE)

  n <- DBI::dbExecute(con, "
    INSERT INTO histo_measurements (sample_id, assay_id, metric, value, source_file)
    SELECT sample_id, assay_id, metric, value, source_file
    FROM tmp_histo_measurements
    ON CONFLICT (sample_id, assay_id, metric) DO NOTHING
  ")
  invisible(n)
}
