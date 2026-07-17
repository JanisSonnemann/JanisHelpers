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
