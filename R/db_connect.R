#' Connect to the lab measurements database
#'
#' Opens (creating if necessary) a DuckDB database file and ensures the full
#' dimension/fact schema exists, creating any missing tables, sequences, and
#' seed rows. Safe to call repeatedly -- existing tables and data are left
#' untouched. This includes `import_log` and `import_exclusions` (history
#' only, no longer written) and the versioned exclusion rules
#' (`exclusion_sets`, `exclusion_rules`).
#' Additive columns introduced later (`facs_measurements.population_order`,
#' `subjects.early` / `early_reason` / `id_source`, `facs_stains.file_name`)
#' are added to an existing file on open; earlier rows stay `NULL`
#' (`early` = `FALSE`, `id_source` = `"mouse_ID"`).
#' The views `facs_exclusion_flags`, `elisa_exclusion_flags` and
#' `histo_exclusion_flags` (one row per measurement matched by the current
#' exclusion rule set of its experiment) are (re)created on every open.
#'
#' Errors if the file was created by an older schema (one whose
#' `facs_stains` table has no `wsp_root` or `condition` column, or whose
#' `samples.collected_at` is nullable, or whose `elisa_measurements` has no
#' `analyte` column) -- such a file must be recreated, since
#' `CREATE TABLE IF NOT EXISTS` would silently keep the old tables.
#'
#' @param path Path to the DuckDB database file. Created if it doesn't exist.
#' @returns A `DBI` connection (class `duckdb_connection`).
#' @export
db_connect <- function(path) {
  con <- DBI::dbConnect(duckdb::duckdb(), path)
  reasons <- db_old_schema_(con)
  if (length(reasons) > 0) {
    DBI::dbDisconnect(con, shutdown = TRUE)
    stop(
      "'", path, "' was created by an older schema (",
      paste(reasons, collapse = "; "), "). Recreate the database file and ",
      "re-import.",
      call. = FALSE
    )
  }
  create_schema_(con)
  con
}

# Reasons why an existing database file predates the current schema (empty =
# current or new). Checked before create_schema_(), which would silently keep
# old tables.
db_old_schema_ <- function(con) {
  reasons <- character()
  if (DBI::dbExistsTable(con, "facs_stains")) {
    missing_cols <- setdiff(
      c("wsp_root", "condition"), DBI::dbListFields(con, "facs_stains")
    )
    if (length(missing_cols) > 0) {
      reasons <- c(reasons, paste0(
        "facs_stains has no ", paste(missing_cols, collapse = "/"), " column"
      ))
    }
  }
  if (DBI::dbExistsTable(con, "samples")) {
    nullable <- DBI::dbGetQuery(con, "
      SELECT is_nullable FROM information_schema.columns
      WHERE table_name = 'samples' AND column_name = 'collected_at'
    ")$is_nullable
    if (identical(nullable, "YES")) {
      reasons <- c(reasons, "samples.collected_at is nullable")
    }
  }
  if (DBI::dbExistsTable(con, "elisa_measurements") &&
      !"analyte" %in% DBI::dbListFields(con, "elisa_measurements")) {
    reasons <- c(reasons, "elisa_measurements has no analyte column")
  }
  reasons
}

create_schema_ <- function(con) {
  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_domains START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS domains (
      domain_id   INTEGER PRIMARY KEY DEFAULT nextval('seq_domains'),
      domain_name VARCHAR NOT NULL UNIQUE
    )
  ")
  DBI::dbExecute(con, "
    INSERT INTO domains (domain_name)
    SELECT * FROM (VALUES ('facs'), ('histo'), ('elisa')) AS v(domain_name)
    ON CONFLICT (domain_name) DO NOTHING
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_experiments START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS experiments (
      experiment_id   INTEGER PRIMARY KEY DEFAULT nextval('seq_experiments'),
      experiment_code VARCHAR NOT NULL UNIQUE,
      experiment_name VARCHAR,
      project         VARCHAR,
      description     VARCHAR,
      created_at      TIMESTAMP DEFAULT current_timestamp
    )
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_subjects START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS subjects (
      subject_id      INTEGER PRIMARY KEY DEFAULT nextval('seq_subjects'),
      mouse_id        VARCHAR NOT NULL UNIQUE,
      experiment_id   INTEGER NOT NULL REFERENCES experiments(experiment_id),
      mdc_id          VARCHAR,
      cage            VARCHAR,
      mouse_strain    VARCHAR,
      generation      VARCHAR,
      sex             VARCHAR,
      mouse_treatment VARCHAR,
      treatment_group VARCHAR,
      dob             DATE,
      start_date      DATE,
      bmt_date        DATE,
      end_date        DATE,
      created_at      TIMESTAMP DEFAULT current_timestamp,
      -- dropout info: early = left the experiment before the end point;
      -- id_source = 'mouse_ID', or 'facility' when mouse_id was taken from ID
      early           BOOLEAN NOT NULL DEFAULT FALSE,
      early_reason    VARCHAR,
      id_source       VARCHAR NOT NULL DEFAULT 'mouse_ID'
    )
  ")
  # Additive migration for files created before the dropout columns existed.
  # (DuckDB cannot add a NOT NULL column; the DEFAULT fills existing rows.)
  DBI::dbExecute(con, "ALTER TABLE subjects ADD COLUMN IF NOT EXISTS early BOOLEAN DEFAULT FALSE")
  DBI::dbExecute(con, "ALTER TABLE subjects ADD COLUMN IF NOT EXISTS early_reason VARCHAR")
  DBI::dbExecute(con, "ALTER TABLE subjects ADD COLUMN IF NOT EXISTS id_source VARCHAR DEFAULT 'mouse_ID'")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_samples START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS samples (
      sample_id    INTEGER PRIMARY KEY DEFAULT nextval('seq_samples'),
      subject_id   INTEGER NOT NULL REFERENCES subjects(subject_id),
      tissue       VARCHAR NOT NULL,
      collected_at DATE NOT NULL,
      created_at   TIMESTAMP DEFAULT current_timestamp,
      UNIQUE(subject_id, tissue, collected_at)
    )
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_assays START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS assays (
      assay_id    INTEGER PRIMARY KEY DEFAULT nextval('seq_assays'),
      assay_name  VARCHAR NOT NULL,
      domain_id   INTEGER NOT NULL REFERENCES domains(domain_id),
      description VARCHAR,
      UNIQUE(assay_name, domain_id)
    )
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_facs_stains START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS facs_stains (
      facs_stain_id INTEGER PRIMARY KEY DEFAULT nextval('seq_facs_stains'),
      sample_id     INTEGER NOT NULL REFERENCES samples(sample_id),
      assay_id      INTEGER NOT NULL REFERENCES assays(assay_id),
      wsp_root      VARCHAR NOT NULL,
      -- the value of the one variable that differs between wells of the same
      -- sample x panel (e.g. 'MPO409'), and the FlowJo keyword it came from
      -- (e.g. 'restimulation'). 'none' / NULL for panels with no such variable.
      condition         VARCHAR NOT NULL DEFAULT 'none',
      condition_keyword VARCHAR,
      stain_date    DATE,
      -- FlowJo file name, matched by file-level exclusion rules
      file_name     VARCHAR,
      UNIQUE(sample_id, assay_id, condition)
    )
  ")
  DBI::dbExecute(con, "ALTER TABLE facs_stains ADD COLUMN IF NOT EXISTS file_name VARCHAR")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_facs_measurements START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS facs_measurements (
      facs_measurement_id  BIGINT PRIMARY KEY DEFAULT nextval('seq_facs_measurements'),
      -- no FK to facs_stains: DuckDB can't delete a parent row in the same
      -- transaction as its children, which replacing a root requires.
      -- See docs/superpowers/specs/2026-11-12-db-schema-import-log-design.md
      facs_stain_id        INTEGER NOT NULL,
      population_full_path VARCHAR NOT NULL,
      population           VARCHAR NOT NULL,
      metric               VARCHAR NOT NULL,
      value                DOUBLE,
      source_file          VARCHAR,
      imported_at          TIMESTAMP DEFAULT current_timestamp,
      -- position of `population` in the panel's display order (prep factor
      -- levels); last column so new and migrated files have the same layout
      population_order     INTEGER,
      UNIQUE(facs_stain_id, population_full_path, metric)
    )
  ")
  # Additive migration for files created before population_order existed;
  # their rows stay NULL until the panel is re-exported.
  DBI::dbExecute(con, "
    ALTER TABLE facs_measurements ADD COLUMN IF NOT EXISTS population_order INTEGER
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_elisa_measurements START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS elisa_measurements (
      elisa_measurement_id INTEGER PRIMARY KEY DEFAULT nextval('seq_elisa_measurements'),
      sample_id            INTEGER NOT NULL REFERENCES samples(sample_id),
      assay_id             INTEGER NOT NULL REFERENCES assays(assay_id),
      import_unit          VARCHAR NOT NULL,
      analyte              VARCHAR NOT NULL,
      replicate            INTEGER NOT NULL,
      value                DOUBLE,
      unit                 VARCHAR,
      result_status        VARCHAR NOT NULL,
      lloq                 DOUBLE,
      uloq                 DOUBLE,
      source_file          VARCHAR,
      imported_at          TIMESTAMP DEFAULT current_timestamp,
      UNIQUE(sample_id, assay_id, analyte, replicate)
    )
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_histo_measurements START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS histo_measurements (
      histo_measurement_id INTEGER PRIMARY KEY DEFAULT nextval('seq_histo_measurements'),
      sample_id            INTEGER NOT NULL REFERENCES samples(sample_id),
      assay_id             INTEGER NOT NULL REFERENCES assays(assay_id),
      import_unit          VARCHAR NOT NULL,
      metric               VARCHAR NOT NULL,
      value                DOUBLE,
      source_file          VARCHAR,
      imported_at          TIMESTAMP DEFAULT current_timestamp,
      UNIQUE(sample_id, assay_id, metric)
    )
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_import_log START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS import_log (
      import_id            INTEGER PRIMARY KEY DEFAULT nextval('seq_import_log'),
      domain_id            INTEGER NOT NULL REFERENCES domains(domain_id),
      experiment_code      VARCHAR,
      import_unit          VARCHAR NOT NULL,
      source_file          VARCHAR NOT NULL,
      source_version       DATE,
      qc_report            VARCHAR,
      qc_passed            BOOLEAN,
      rows_written         INTEGER NOT NULL,
      package_version      VARCHAR NOT NULL,
      -- no REFERENCES import_log(import_id): DuckDB's EXPORT DATABASE (the
      -- snapshot) rejects self-referencing foreign keys. write_import_log_()
      -- checks the referenced row exists instead.
      supersedes_import_id INTEGER,
      imported_at          TIMESTAMP NOT NULL DEFAULT current_timestamp
    )
  ")

  # Exclusion-sheet rows applied by one FACS export (db_export_facs()). Kept
  # when the root is overwritten later, so every import_log row keeps its own
  # record of what was left out and why. The FK points at another table, which
  # EXPORT DATABASE handles fine.
  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_import_exclusions START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS import_exclusions (
      import_exclusion_id INTEGER PRIMARY KEY DEFAULT nextval('seq_import_exclusions'),
      import_id       INTEGER NOT NULL REFERENCES import_log(import_id),
      mouse_id        VARCHAR NOT NULL,
      tissue          VARCHAR,
      panel           VARCHAR,
      condition       VARCHAR,
      file_name       VARCHAR,
      collected_at    DATE,
      reason          VARCHAR NOT NULL,
      decided_on      DATE,
      n_files_matched INTEGER NOT NULL
    )
  ")

  # NOTE: import_exclusions above is history only -- it is no longer written.
  # The live rules are exclusion_sets / exclusion_rules (spec:
  # docs/superpowers/specs/2026-11-20-exclusions-in-db-design.md). No
  # self-referencing FKs (EXPORT DATABASE rejects them).
  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_exclusion_sets START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS exclusion_sets (
      exclusion_set_id INTEGER PRIMARY KEY DEFAULT nextval('seq_exclusion_sets'),
      experiment_id    INTEGER NOT NULL REFERENCES experiments(experiment_id),
      source_file      VARCHAR,
      written_at       TIMESTAMP NOT NULL DEFAULT current_timestamp,
      package_version  VARCHAR NOT NULL,
      note             VARCHAR
    )
  ")
  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_exclusion_rules START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS exclusion_rules (
      exclusion_rule_id INTEGER PRIMARY KEY DEFAULT nextval('seq_exclusion_rules'),
      exclusion_set_id  INTEGER NOT NULL REFERENCES exclusion_sets(exclusion_set_id),
      mouse_id          VARCHAR NOT NULL,
      tissue            VARCHAR,
      panel             VARCHAR,
      condition         VARCHAR,
      file_name         VARCHAR,
      collected_at      DATE,
      reason            VARCHAR NOT NULL,
      decided_on        DATE
    )
  ")

  # Exclusion flags, computed on read from each experiment's current rule set
  # (spec "Matching on read"). One row per flagged measurement; rows without a
  # match are absent (the db_query_*() left joins turn that into FALSE).
  # CREATE OR REPLACE also migrates older files on open.
  for (domain in c("facs", "elisa", "histo")) {
    id <- paste0(domain, "_measurement_id")
    DBI::dbExecute(con, paste0(
      "CREATE OR REPLACE VIEW ", domain, "_exclusion_flags AS ",
      "SELECT m.", id, ", TRUE AS excluded, ",
      "string_agg(DISTINCT r.reason, '; ' ORDER BY r.reason) AS exclusion_reason ",
      exclusion_match_from_(domain),
      " GROUP BY m.", id
    ))
  }

  invisible(con)
}

# Internal: the FROM / JOIN clause matching each measurement of `domain`
# ("facs", "elisa", "histo"; aliases m = measurement, sa = sample, su =
# subject, a = assay, fs = FACS stain) against the rules `r` of its
# experiment's current set. A rule matches when mouse_id is equal and every
# scope column is NULL or equal (`panel` vs `assay_name`). ELISA / histology
# have no condition / file_name, so rules naming one never match there.
# `file_match` replaces the FACS file_name condition (used to find file rules
# whose stains were stored without a file name).
exclusion_match_from_ <- function(domain,
                                  file_match = "(r.file_name IS NULL OR r.file_name = fs.file_name)") {
  if (domain == "facs") {
    joins <- "
      FROM facs_measurements m
      JOIN facs_stains fs ON fs.facs_stain_id = m.facs_stain_id
      JOIN samples sa     ON sa.sample_id = fs.sample_id
      JOIN assays a       ON a.assay_id = fs.assay_id"
    extra <- paste0("
       AND (r.condition IS NULL OR r.condition = fs.condition)
       AND ", file_match)
  } else {
    joins <- paste0("
      FROM ", domain, "_measurements m
      JOIN samples sa     ON sa.sample_id = m.sample_id
      JOIN assays a       ON a.assay_id = m.assay_id")
    extra <- "
       AND r.condition IS NULL AND r.file_name IS NULL"
  }
  paste0(joins, "
      JOIN subjects su    ON su.subject_id = sa.subject_id
      JOIN (SELECT experiment_id, MAX(exclusion_set_id) AS exclusion_set_id
            FROM exclusion_sets GROUP BY experiment_id) cur
        ON cur.experiment_id = su.experiment_id
      JOIN exclusion_rules r ON r.exclusion_set_id = cur.exclusion_set_id
       AND r.mouse_id = su.mouse_id
       AND (r.tissue IS NULL OR r.tissue = sa.tissue)
       AND (r.panel IS NULL OR r.panel = a.assay_name)
       AND (r.collected_at IS NULL OR r.collected_at = sa.collected_at)",
    extra)
}
