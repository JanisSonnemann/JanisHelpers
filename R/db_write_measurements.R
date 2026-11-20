# Errors unless `x` is a single non-NA, non-empty string. `NULL` stands for a
# missing argument. `what` describes the expected value in the message.
check_string_arg_ <- function(x, arg, what) {
  if (is.null(x) || !is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x)) {
    stop("`", arg, "` must be a single non-empty string: ", what, ".", call. = FALSE)
  }
  invisible(x)
}

# Errors unless `data` has the facs_import_wsp() columns, at least one row,
# and (if it has FileName) rows from a single FCS file. Without this, a
# wrong-shaped or empty tibble would create a stain with no measurements,
# which the version guard can't see.
check_facs_data_ <- function(data) {
  needed <- c("PopulationFullPath", "Population", "metric", "value")
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "db_write_facs() needs column(s) ",
      paste(sprintf("`%s`", missing_cols), collapse = ", "),
      " in `data` (the facs_import_wsp() shape: ",
      paste(needed, collapse = ", "), ").",
      call. = FALSE
    )
  }
  if (nrow(data) == 0) {
    stop("db_write_facs() got `data` with no rows -- nothing to write.", call. = FALSE)
  }
  if ("FileName" %in% names(data) && length(unique(data$FileName)) > 1) {
    stop(
      "db_write_facs() writes one stain (one FCS file) per call, but `data` has ",
      length(unique(data$FileName)), " distinct FileName values -- filter to one first.",
      call. = FALSE
    )
  }
  invisible(TRUE)
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

# condition and condition_keyword come as a pair: a real condition needs the
# keyword it came from; the "none" placeholder must not have one.
check_facs_condition_ <- function(condition, condition_keyword) {
  check_string_arg_(
    condition, "condition",
    "the well's value of the varying variable (e.g. \"MPO409\"), or \"none\""
  )
  if (length(condition_keyword) != 1) {
    stop("`condition_keyword` must be a single string or NA.", call. = FALSE)
  }
  has_keyword <- !is.na(condition_keyword) && nzchar(condition_keyword)
  if (condition == "none" && has_keyword) {
    stop(
      "`condition_keyword` '", condition_keyword, "' was given but `condition` ",
      "is \"none\" -- pass the well's value as `condition`, or drop the keyword.",
      call. = FALSE
    )
  }
  if (condition != "none" && !has_keyword) {
    stop(
      "`condition` '", condition, "' needs `condition_keyword` -- the FlowJo ",
      "keyword it came from (e.g. \"restimulation\").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# One root varies at most one variable: all its stains share one
# condition_keyword (or all have none).
check_facs_root_keyword_ <- function(con, wsp_root, condition_keyword) {
  existing <- DBI::dbGetQuery(con, "
    SELECT DISTINCT COALESCE(condition_keyword, '') AS keyword
    FROM facs_stains WHERE wsp_root = ?",
    params = list(wsp_root))$keyword
  given <- if (is.na(condition_keyword)) "" else condition_keyword
  other <- setdiff(existing, given)
  if (length(other) > 0) {
    show <- function(k) if (nzchar(k)) paste0("'", k, "'") else "no keyword"
    stop(
      "Root '", wsp_root, "' already holds stains with condition_keyword ",
      show(other[[1]]), ", not ", show(given), ". One root varies at most one ",
      "variable -- check the workspace's keywords.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

lookup_sample_id_ <- function(con, mouse_id, tissue, collected_at = NA) {
  lookup_sample_ids_(con, mouse_id, tissue, collected_at)[[1]]
}

# Vectorised: one sample_id per (mouse_id, tissue, collected_at) in input
# order; NA collected_at = the subject's end_date. Errors naming every
# unresolved mouse x tissue x date (NA counts as unresolved).
lookup_sample_ids_ <- function(con, mouse_id, tissue, collected_at = NA) {
  keys <- tibble::tibble(
    row_ = seq_along(mouse_id),
    mouse_id = as.character(mouse_id),
    tissue = as.character(tissue),
    collected_at = as.Date(rep(collected_at, length.out = length(mouse_id)))
  )
  duckdb::duckdb_register(con, "tmp_sample_keys", keys)
  on.exit(duckdb::duckdb_unregister(con, "tmp_sample_keys"), add = TRUE)

  resolved <- DBI::dbGetQuery(con, "
    SELECT k.row_, k.mouse_id, k.tissue, sm.sample_id,
      CAST(COALESCE(k.collected_at, sub.end_date) AS VARCHAR) AS date_
    FROM tmp_sample_keys k
    LEFT JOIN subjects sub ON sub.mouse_id = k.mouse_id
    LEFT JOIN samples sm ON sm.subject_id = sub.subject_id
      AND sm.tissue = k.tissue
      AND sm.collected_at = COALESCE(k.collected_at, sub.end_date)
    ORDER BY k.row_
  ")

  bad <- resolved[is.na(resolved$sample_id), ]
  bad <- unique(bad[c("mouse_id", "tissue", "date_")])
  if (nrow(bad) > 0) {
    stop(
      "No sample found for mouse_id x tissue x date: ",
      paste(sprintf("'%s x %s x %s'", bad$mouse_id, bad$tissue, bad$date_), collapse = ", "),
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
#'   workspace version replaces a root via [db_export_facs()], not here);
#' - the root already holds stains with a different `condition_keyword`
#'   (one root varies at most one variable).
#'
#' Wells of the same sample x panel that differ in one variable (e.g. a
#' restimulation: `unstim` / `MPO409` / `PMA-Iono`) are stored as separate
#' stains: pass the well's value as `condition` and the FlowJo keyword it
#' came from as `condition_keyword`. They share one `wsp_root`.
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
#'   measurements came from, as a path relative to the project's `data/`
#'   folder (e.g. `"25-7_exp/raw/FACS/treg_25-7_20261112.wsp"`). Recorded on every
#'   measurement row as provenance. The FCS filename is deliberately not used:
#'   it doesn't identify which workspace version (gates, compensation)
#'   produced the numbers.
#' @param condition Character scalar -- the value of the variable that differs
#'   between wells of this sample x panel (e.g. `"MPO409"`). Default `"none"`
#'   for panels with one well per sample.
#' @param condition_keyword Character scalar -- the FlowJo keyword `condition`
#'   came from (e.g. `"restimulation"`). Required when `condition` isn't
#'   `"none"`, must be `NA` when it is.
#' @param stain_date Date scalar, optional.
#' @returns Invisibly, the number of measurement rows inserted.
#' @export
db_write_facs <- function(con, data, mouse_id, tissue, assay_name, wsp_root,
                           source_file, condition = "none",
                           condition_keyword = NA_character_, stain_date = NA) {
  check_string_arg_(
    if (missing(wsp_root)) NULL else wsp_root, "wsp_root",
    "the .wsp filename without its _YYYYMMDD suffix (e.g. \"treg_25-7\")"
  )
  check_string_arg_(
    if (missing(source_file)) NULL else source_file, "source_file",
    "the .wsp path relative to data/ (e.g. \"25-7_exp/raw/FACS/treg_25-7_20261112.wsp\")"
  )
  check_facs_condition_(condition, condition_keyword)
  check_facs_data_(data)

  sample_id <- lookup_sample_id_(con, mouse_id, tissue)
  assay_id <- lookup_assay_id_(con, assay_name, domain = "facs")
  check_facs_root_(con, sample_id, assay_id, wsp_root, source_file)
  check_facs_root_keyword_(con, wsp_root, condition_keyword)

  measurements <- tibble::tibble(
    population_full_path = data$PopulationFullPath,
    population = data$Population,
    metric = data$metric,
    value = data$value
  )
  n <- DBI::dbWithTransaction(con, {
    insert_facs_stain_(
      con, sample_id, assay_id, measurements, wsp_root, source_file,
      condition, condition_keyword, stain_date
    )
  })

  invisible(n)
}

# Internal: insert one stain (if new) and its measurements. Opens NO
# transaction -- callers wrap it (db_write_facs() per stain, db_export_facs()
# around a whole root), since DuckDB doesn't nest transactions. `data` has
# snake-case columns population_full_path, population, metric, value. A
# factor `population` also stores each label's position in the display order
# (`population_order`); a character one stores NULL.
# Returns the number of measurement rows inserted.
insert_facs_stain_ <- function(con, sample_id, assay_id, data, wsp_root,
                               source_file, condition, condition_keyword,
                               stain_date = NA, file_name = NA_character_) {
  population_order <- if (is.factor(data$population)) {
    as.integer(data$population)
  } else {
    rep(NA_integer_, length(data$population))
  }
  to_write <- tibble::tibble(
    population_full_path = as.character(data$population_full_path),
    population = as.character(data$population),
    population_order = population_order,
    metric = as.character(data$metric),
    value = as.double(data$value),
    source_file = source_file
  )
  duckdb::duckdb_register(con, "tmp_facs_measurements", to_write)
  on.exit(duckdb::duckdb_unregister(con, "tmp_facs_measurements"), add = TRUE)

  DBI::dbExecute(
    con,
    "INSERT INTO facs_stains (
       sample_id, assay_id, wsp_root, condition, condition_keyword, stain_date,
       file_name
     )
     VALUES (?, ?, ?, ?, ?, ?, ?)
     ON CONFLICT (sample_id, assay_id, condition) DO NOTHING",
    params = list(
      sample_id, assay_id, wsp_root, condition,
      as.character(condition_keyword), stain_date, as.character(file_name)
    )
  )

  facs_stain_id <- DBI::dbGetQuery(
    con,
    "SELECT facs_stain_id FROM facs_stains
     WHERE sample_id = ? AND assay_id = ? AND condition = ?",
    params = list(sample_id, assay_id, condition)
  )$facs_stain_id[[1]]

  DBI::dbExecute(con, "
    INSERT INTO facs_measurements (
      facs_stain_id, population_full_path, population, population_order,
      metric, value, source_file
    )
    SELECT ?, population_full_path, population, population_order,
      metric, value, source_file
    FROM tmp_facs_measurements
    ON CONFLICT (facs_stain_id, population_full_path, metric) DO NOTHING
  ", params = list(facs_stain_id))
}
