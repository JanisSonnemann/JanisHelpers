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
#' [db_write_experiment()]. Rows whose `mouse_id` already exists are not
#' inserted again; their `early` / `early_reason` are updated when `data` has
#' those columns and a stored `NA` `end_date` is filled from `data`; all other
#' columns stay untouched. Errors (writing nothing)
#' if a stored `end_date` differs from a non-`NA` incoming one -- end-point
#' samples are keyed by it. Also errors if an incoming `mdc_id` is already
#' stored for a different `mouse_id` (one facility ID = one subject; a dropout
#' cannot be re-keyed). Errors if any row's `experiment_code` isn't registered yet -- naming
#' every unresolved code -- rather than silently dropping those rows.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param data A tibble with columns `mouse_id`, `experiment_code`, and
#'   optionally `mdc_id`, `cage`, `mouse_strain`, `generation`, `sex`,
#'   `mouse_treatment`, `treatment_group`, `dob`, `start_date`, `bmt_date`,
#'   `end_date`, `early` (logical, default `FALSE`), `early_reason` and
#'   `id_source` (`"mouse_ID"` (default) or `"facility"`). Missing optional
#'   columns are filled with `NA`.
#' @returns Invisibly, the number of rows newly inserted (updates of
#'   `early` / `early_reason` / `end_date` do not count).
#' @export
db_write_subjects <- function(con, data) {
  update_early <- any(c("early", "early_reason") %in% names(data))
  data <- fill_missing_cols_(data, c(
    "early", "early_reason", "id_source",
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

  changed <- DBI::dbGetQuery(con, "
    SELECT s.mouse_id, sub.end_date AS stored, CAST(s.end_date AS DATE) AS new
    FROM tmp_subjects s
    JOIN subjects sub ON sub.mouse_id = s.mouse_id
    WHERE sub.end_date IS NOT NULL AND s.end_date IS NOT NULL
      AND sub.end_date <> CAST(s.end_date AS DATE)
  ")
  if (nrow(changed) > 0) {
    stop(
      "end_date differs from the stored value (end-point samples are keyed by it); nothing written:\n",
      paste0("  ", changed$mouse_id, ": stored ", changed$stored,
             ", new ", changed$new, collapse = "\n"),
      call. = FALSE
    )
  }

  rekeyed <- DBI::dbGetQuery(con, "
    SELECT DISTINCT sub.mouse_id AS stored, s.mouse_id AS new, s.mdc_id
    FROM tmp_subjects s
    JOIN subjects sub ON sub.mdc_id = s.mdc_id
    WHERE s.mdc_id IS NOT NULL AND sub.mouse_id <> s.mouse_id
  ")
  if (nrow(rekeyed) > 0) {
    stop(
      "Facility ID already registered under a different mouse_id; nothing written:\n",
      paste0("  Facility ID '", rekeyed$mdc_id, "' is already registered as '",
             rekeyed$stored, "' -- it now has mouse_ID '", rekeyed$new, "'",
             collapse = "\n"),
      "\nRe-keying a subject is not supported; correct the sheet or remove the old subject.",
      call. = FALSE
    )
  }

  n_new <- DBI::dbGetQuery(con, "
    SELECT COUNT(DISTINCT s.mouse_id) AS n
    FROM tmp_subjects s
    WHERE s.mouse_id NOT IN (SELECT mouse_id FROM subjects)
  ")$n

  # a stored NA end_date is filled from the incoming one (a changed non-NA
  # one was refused above); early / early_reason only when `data` has them
  on_conflict <- paste0(
    "DO UPDATE SET end_date = COALESCE(subjects.end_date, excluded.end_date)",
    if (update_early) ", early = excluded.early, early_reason = excluded.early_reason"
  )
  DBI::dbExecute(con, paste("
    INSERT INTO subjects (
      mouse_id, experiment_id, mdc_id, cage, mouse_strain, generation,
      sex, mouse_treatment, treatment_group, dob, start_date, bmt_date, end_date,
      early, early_reason, id_source
    )
    SELECT
      s.mouse_id, e.experiment_id, s.mdc_id, s.cage, s.mouse_strain, s.generation,
      s.sex, s.mouse_treatment, s.treatment_group,
      CAST(s.dob AS DATE), CAST(s.start_date AS DATE),
      CAST(s.bmt_date AS DATE), CAST(s.end_date AS DATE),
      COALESCE(CAST(s.early AS BOOLEAN), FALSE), CAST(s.early_reason AS VARCHAR),
      COALESCE(CAST(s.id_source AS VARCHAR), 'mouse_ID')
    FROM tmp_subjects s
    JOIN experiments e ON e.experiment_code = s.experiment_code
    ON CONFLICT (mouse_id)", on_conflict))
  invisible(as.numeric(n_new))
}

#' Register samples (subject x tissue x collection date)
#'
#' Bulk-inserts sample rows. `data` must already be shaped with one row per
#' subject x tissue x collection date, referencing subjects already registered via
#' [db_write_subjects()]. Errors if any row's `mouse_id` isn't registered yet
#' -- naming every unresolved `mouse_id` -- rather than silently dropping
#' those rows. Organ weights and volumes are not stored (they are calculation
#' inputs for the FACS export script); such columns in `data` are ignored.
#' A missing (`NA` or absent) `collected_at` means the subject's `end_date`;
#' errors naming every mouse with neither.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param data A tibble with columns `mouse_id`, `tissue`, and optionally
#'   `collected_at`. Other columns are ignored.
#' @returns Invisibly, the number of rows inserted into `samples`.
#' @export
db_write_samples <- function(con, data) {
  data <- fill_missing_cols_(data, c("mouse_id", "tissue", "collected_at"))
  data <- data[c("mouse_id", "tissue", "collected_at")]
  duckdb::duckdb_register(con, "tmp_samples", data)
  on.exit(duckdb::duckdb_unregister(con, "tmp_samples"), add = TRUE)

  unknown_mice <- DBI::dbGetQuery(con, "
    SELECT DISTINCT t.mouse_id
    FROM tmp_samples t
    LEFT JOIN subjects sub ON sub.mouse_id = t.mouse_id
    WHERE sub.subject_id IS NULL
  ")$mouse_id
  if (length(unknown_mice) > 0) {
    stop(
      "Unknown mouse_id(s): ", paste(sprintf("'%s'", unknown_mice), collapse = ", "),
      " -- register them first via db_write_subjects().",
      call. = FALSE
    )
  }

  no_date <- DBI::dbGetQuery(con, "
    SELECT DISTINCT t.mouse_id
    FROM tmp_samples t
    JOIN subjects sub ON sub.mouse_id = t.mouse_id
    WHERE CAST(t.collected_at AS DATE) IS NULL AND sub.end_date IS NULL
  ")$mouse_id
  if (length(no_date) > 0) {
    stop(
      "No collected_at and no subject end_date for mouse_id(s): ",
      paste(sprintf("'%s'", no_date), collapse = ", "),
      " -- give the subject an end_date or pass collected_at.",
      call. = FALSE
    )
  }

  n <- DBI::dbExecute(con, "
    INSERT INTO samples (subject_id, tissue, collected_at)
    SELECT DISTINCT sub.subject_id, t.tissue,
      COALESCE(CAST(t.collected_at AS DATE), sub.end_date)
    FROM tmp_samples t
    JOIN subjects sub ON sub.mouse_id = t.mouse_id
    ON CONFLICT (subject_id, tissue, collected_at) DO NOTHING
  ")

  invisible(n)
}

#' Register an assay/panel
#'
#' Inserts a new assay identified by `assay_name` within `domain`. If an
#' assay with that name already exists in that domain, the call is a no-op.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param assay_name Character scalar, e.g. `"overview"` or `"anti-MPO"`.
#' @param domain Character scalar, one of `"facs"`, `"histo"`, `"elisa"` (or
#'   any additional domain later added to the `domains` table).
#' @param description Character scalar, optional.
#' @returns Invisibly, the number of rows inserted (`0` or `1`).
#' @export
db_write_assay <- function(con, assay_name, domain, description = NA_character_) {
  domain_row <- DBI::dbGetQuery(
    con, "SELECT domain_id FROM domains WHERE domain_name = ?",
    params = list(domain)
  )
  if (nrow(domain_row) == 0) {
    stop(
      "Unknown domain '", domain, "' -- expected one already present in the domains table.",
      call. = FALSE
    )
  }

  n <- DBI::dbExecute(
    con,
    "INSERT INTO assays (assay_name, domain_id, description)
     VALUES (?, ?, ?)
     ON CONFLICT (assay_name, domain_id) DO NOTHING",
    params = list(assay_name, domain_row$domain_id[[1]], description)
  )
  invisible(n)
}

#' Register an experiment and its mice from the meta workbook
#'
#' Registers the experiment ([db_write_experiment()]) and its subjects
#' ([db_write_subjects()]) from the `meta` sheet of a meta workbook built from
#' the standard template (`meta_template.xlsm`), as read by [meta_read()].
#' Every row with a `mouse_ID` or a facility `id` is registered (the key
#' `mouse_id` is the `mouse_ID`, else the `id`; `id_source` records which); a row
#' with neither is an error. Re-running updates `early` / `early_reason` of
#' stored mice but refuses `end_date` changes (see [db_write_subjects()]).
#' Idempotent. For workbooks that don't
#' follow the template, shape the subjects yourself and call
#' [db_write_subjects()].
#'
#' Column mapping (sheet -> `subjects`): `mouse_ID` -> `mouse_id`, `id` ->
#' `mdc_id`, `cage`, `line` -> `mouse_strain`, `generation`, `sex`, `group`
#' -> `treatment_group`, `dob`, `start_date`, `bmt_date`, `death_date` ->
#' `end_date`. Optional columns: `early` (`yes`/`ja`/`y`/`x`/`true` -> `TRUE`;
#' empty/`no`/`nein`/`n`/`false` -> `FALSE`; anything else is an error) and
#' `reason` -> `early_reason`.
#'
#' @param con A `DBI` connection from [db_connect()].
#' @param meta The list returned by [meta_read()] (its `meta` element is
#'   used), or the `meta` sheet itself.
#' @param experiment_code Character scalar, e.g. `"26-1"`.
#' @param experiment_name,project Character scalars, optional.
#' @returns Invisibly, `con` -- so it can follow [db_connect()] in a pipe.
#' @examples
#' \dontrun{
#' meta <- meta_read("data/26-1_exp/meta/meta_26-1.xlsx")
#' con <- db_connect("~/databases/lab.duckdb") |>
#'   db_write_meta(meta, "26-1", "BMT-mRNA-V4", project = "MPO_mRNA")
#' }
#' @export
db_write_meta <- function(con, meta, experiment_code,
                          experiment_name = NA_character_,
                          project = NA_character_) {
  sheet <- if (is.data.frame(meta)) meta else meta$meta
  if (!is.data.frame(sheet)) {
    stop("`meta` must be meta_read() output or its `meta` sheet.", call. = FALSE)
  }
  needed <- c("mouse_ID", "id", "cage", "line", "generation", "sex", "group",
              "dob", "start_date", "bmt_date", "death_date")
  missing_cols <- setdiff(needed, names(sheet))
  if (length(missing_cols) > 0) {
    stop(
      "The meta sheet is missing template column(s): ",
      paste(missing_cols, collapse = ", "),
      " -- fix the workbook or use db_write_subjects() directly.",
      call. = FALSE
    )
  }

  no_key <- is.na(sheet$mouse_ID) & is.na(sheet$id)
  # fully blank rows are skipped; a row with other data but no key is an error
  other <- sheet[setdiff(names(sheet), c("mouse_ID", "id"))]
  stray <- which(no_key & rowSums(!is.na(other)) > 0)
  if (length(stray) > 0) {
    stop(
      "Row(s) ", paste(stray, collapse = ", "),
      " of the meta sheet have neither mouse_ID nor id.",
      call. = FALSE
    )
  }
  sheet$early <- if ("early" %in% names(sheet)) parse_early_(sheet$early) else FALSE
  sheet$reason <- if ("reason" %in% names(sheet)) as.character(sheet$reason) else NA_character_

  subjects <- sheet |>
    dplyr::filter(!no_key) |>
    dplyr::transmute(
      mouse_id        = dplyr::coalesce(as.character(mouse_ID), as.character(id)),
      id_source       = ifelse(is.na(mouse_ID), "facility", "mouse_ID"),
      early           = early,
      early_reason    = reason,
      experiment_code = experiment_code,
      mdc_id          = as.character(id),
      cage            = as.character(cage),
      mouse_strain    = as.character(line),
      generation      = as.character(generation),
      sex             = as.character(sex),
      treatment_group = as.character(group),
      dob             = as.Date(dob),
      start_date      = as.Date(start_date),
      bmt_date        = as.Date(bmt_date),
      end_date        = as.Date(death_date)
    )

  db_write_experiment(con, experiment_code, experiment_name, project = project)
  n <- db_write_subjects(con, subjects)
  message("Experiment ", experiment_code, ": ", n, " new of ",
          nrow(subjects), " mice registered.")
  invisible(con)
}

# Parse the meta sheet's `early` column; errors on values it cannot interpret.
parse_early_ <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  v <- tolower(trimws(as.character(x)))
  yes <- c("yes", "ja", "y", "x", "true")
  no <- c("", "no", "nein", "n", "false")
  bad <- unique(x[!is.na(v) & !(v %in% c(yes, no))])
  if (length(bad) > 0) {
    stop(
      "Unknown `early` value(s): ", paste(sprintf("'%s'", bad), collapse = ", "),
      " -- use yes/no (or leave empty).",
      call. = FALSE
    )
  }
  v %in% yes
}
