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
#' `facs_measurements`, creating the matching `facs_stains` row first if it
#' doesn't already exist. Which of `vol_resuspended`/`vol_measured` versus
#' `bead_volume_added`/`bead_concentration`/`bead_population_path` are
#' required depends on `count_method`; the database's `CHECK` constraint on
#' `facs_stains` rejects a call that supplies the wrong subset.
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
#' @param count_method Either `"volumetric"` or `"bead"`.
#' @param vol_stained Numeric scalar, volume of suspension stained.
#' @param vol_resuspended Numeric scalar, required when
#'   `count_method = "volumetric"`.
#' @param vol_measured Numeric scalar, required when
#'   `count_method = "volumetric"`.
#' @param bead_volume_added Numeric scalar, required when
#'   `count_method = "bead"`.
#' @param bead_concentration Numeric scalar, required when
#'   `count_method = "bead"`.
#' @param bead_population_path Character scalar, required when
#'   `count_method = "bead"`.
#' @param stain_date Date scalar, optional.
#' @param source_file Character scalar, optional -- defaults to
#'   `data$FileName[1]` when `data` has a `FileName` column.
#' @returns Invisibly, the number of measurement rows inserted.
#' @export
db_write_facs <- function(con, data, mouse_id, tissue, assay_name, count_method,
                           vol_stained,
                           vol_resuspended = NA_real_, vol_measured = NA_real_,
                           bead_volume_added = NA_real_, bead_concentration = NA_real_,
                           bead_population_path = NA_character_,
                           stain_date = NA, source_file = NULL) {
  count_method <- match.arg(count_method, c("volumetric", "bead"))

  sample_id <- lookup_sample_id_(con, mouse_id, tissue)
  assay_id <- lookup_assay_id_(con, assay_name, domain = "facs")

  DBI::dbExecute(
    con,
    "INSERT INTO facs_stains (
       sample_id, assay_id, vol_stained, count_method,
       vol_resuspended, vol_measured,
       bead_volume_added, bead_concentration, bead_population_path, stain_date
     ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
     ON CONFLICT (sample_id, assay_id) DO NOTHING",
    params = list(
      sample_id, assay_id, vol_stained, count_method,
      vol_resuspended, vol_measured,
      bead_volume_added, bead_concentration, bead_population_path, stain_date
    )
  )

  facs_stain_id <- DBI::dbGetQuery(
    con,
    "SELECT facs_stain_id FROM facs_stains WHERE sample_id = ? AND assay_id = ?",
    params = list(sample_id, assay_id)
  )$facs_stain_id[[1]]

  if (is.null(source_file)) {
    source_file <- if ("FileName" %in% names(data)) as.character(data$FileName[[1]]) else NA_character_
  }

  to_write <- tibble::tibble(
    facs_stain_id = facs_stain_id,
    population_full_path = as.character(data$PopulationFullPath),
    population = as.character(data$Population),
    metric = as.character(data$metric),
    value = as.double(data$value),
    source_file = source_file
  )
  duckdb::duckdb_register(con, "tmp_facs_measurements", to_write)
  on.exit(duckdb::duckdb_unregister(con, "tmp_facs_measurements"), add = TRUE)

  n <- DBI::dbExecute(con, "
    INSERT INTO facs_measurements (
      facs_stain_id, population_full_path, population, metric, value, source_file
    )
    SELECT facs_stain_id, population_full_path, population, metric, value, source_file
    FROM tmp_facs_measurements
    ON CONFLICT (facs_stain_id, population_full_path, metric) DO NOTHING
  ")

  invisible(n)
}

#' Import ELISA measurements
#'
#' Writes a tidy tibble of ELISA measurements (as produced by
#' `elisa_read_results()`) into `elisa_measurements`.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param data A tibble with columns `cytokine`, `sample_id` (raw plate
#'   label), `replicate`, `value`, `unit`, `result_status` -- the shape
#'   `elisa_read_results()` produces.
#' @param mouse_id Character scalar identifying the subject.
#' @param tissue Character scalar identifying the sample.
#' @param assay_name Character scalar identifying the assay/panel (already
#'   registered via [db_write_assay()] with `domain = "elisa"`).
#' @param source_file Character scalar, optional.
#' @returns Invisibly, the number of rows inserted.
#' @export
db_write_elisa <- function(con, data, mouse_id, tissue, assay_name, source_file = NA_character_) {
  sample_id <- lookup_sample_id_(con, mouse_id, tissue)
  assay_id <- lookup_assay_id_(con, assay_name, domain = "elisa")

  to_write <- tibble::tibble(
    sample_id = sample_id,
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
