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
