#' Register an experiment
#'
#' Inserts a new experiment row identified by `experiment_code`. If an
#' experiment with that code already exists, the call is a no-op.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param experiment_code Character scalar, the experiment's natural key
#'   (e.g. `"25-7"`).
#' @param experiment_name Character scalar, optional.
#' @param project Character scalar, optional.
#' @param description Character scalar, optional.
#' @returns Invisibly, the number of rows inserted (`0` or `1`).
#' @export
db_write_experiment <- function(con, experiment_code, experiment_name = NA_character_,
                                 project = NA_character_, description = NA_character_) {
  n <- DBI::dbExecute(
    con,
    "INSERT INTO experiments (experiment_code, experiment_name, project, description)
     VALUES (?, ?, ?, ?)
     ON CONFLICT (experiment_code) DO NOTHING",
    params = list(experiment_code, experiment_name, project, description)
  )
  invisible(n)
}

fill_missing_cols_ <- function(data, cols) {
  missing <- setdiff(cols, names(data))
  data[missing] <- NA
  data
}

#' Register subjects (mice)
#'
#' Bulk-inserts subject rows. `data` must already be shaped to match the
#' `subjects` table: one row per mouse, with `mouse_id` and an
#' `experiment_code` identifying an experiment already registered via
#' [db_write_experiment()]. Rows whose `mouse_id` already exists are skipped.
#' Errors if any row's `experiment_code` isn't registered yet -- naming
#' every unresolved code -- rather than silently dropping those rows.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param data A tibble with columns `mouse_id`, `experiment_code`, and
#'   optionally `mdc_id`, `cage`, `mouse_strain`, `generation`, `sex`,
#'   `mouse_treatment`, `treatment_group`, `dob`, `start_date`, `bmt_date`,
#'   `end_date`. Missing optional columns are filled with `NA`.
#' @returns Invisibly, the number of rows inserted.
#' @export
db_write_subjects <- function(con, data) {
  data <- fill_missing_cols_(data, c(
    "mouse_id", "experiment_code", "mdc_id", "cage", "mouse_strain",
    "generation", "sex", "mouse_treatment", "treatment_group",
    "dob", "start_date", "bmt_date", "end_date"
  ))
  duckdb::duckdb_register(con, "tmp_subjects", data)
  on.exit(duckdb::duckdb_unregister(con, "tmp_subjects"), add = TRUE)

  unknown_codes <- DBI::dbGetQuery(con, "
    SELECT DISTINCT s.experiment_code
    FROM tmp_subjects s
    LEFT JOIN experiments e ON e.experiment_code = s.experiment_code
    WHERE e.experiment_id IS NULL
  ")$experiment_code
  if (length(unknown_codes) > 0) {
    stop(
      "Unknown experiment_code(s): ", paste(sprintf("'%s'", unknown_codes), collapse = ", "),
      " -- register them first via db_write_experiment().",
      call. = FALSE
    )
  }

  n <- DBI::dbExecute(con, "
    INSERT INTO subjects (
      mouse_id, experiment_id, mdc_id, cage, mouse_strain, generation,
      sex, mouse_treatment, treatment_group, dob, start_date, bmt_date, end_date
    )
    SELECT
      s.mouse_id, e.experiment_id, s.mdc_id, s.cage, s.mouse_strain, s.generation,
      s.sex, s.mouse_treatment, s.treatment_group,
      CAST(s.dob AS DATE), CAST(s.start_date AS DATE),
      CAST(s.bmt_date AS DATE), CAST(s.end_date AS DATE)
    FROM tmp_subjects s
    JOIN experiments e ON e.experiment_code = s.experiment_code
    ON CONFLICT (mouse_id) DO NOTHING
  ")
  invisible(n)
}
