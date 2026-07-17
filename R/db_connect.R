#' Connect to the lab measurements database
#'
#' Opens (creating if necessary) a DuckDB database file and ensures the full
#' dimension/fact schema exists, creating any missing tables, sequences, and
#' seed rows. Safe to call repeatedly — existing tables and data are left
#' untouched.
#'
#' @param path Path to the DuckDB database file. Created if it doesn't exist.
#' @returns A `DBI` connection (class `duckdb_connection`).
#' @export
db_connect <- function(path) {
  con <- DBI::dbConnect(duckdb::duckdb(), path)
  create_schema_(con)
  con
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
      created_at      TIMESTAMP DEFAULT current_timestamp
    )
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_samples START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS samples (
      sample_id    INTEGER PRIMARY KEY DEFAULT nextval('seq_samples'),
      subject_id   INTEGER NOT NULL REFERENCES subjects(subject_id),
      tissue       VARCHAR NOT NULL,
      collected_at DATE,
      created_at   TIMESTAMP DEFAULT current_timestamp,
      UNIQUE(subject_id, tissue)
    )
  ")

  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS sample_digestions (
      sample_id    INTEGER PRIMARY KEY REFERENCES samples(sample_id),
      total_weight DOUBLE,
      facs_weight  DOUBLE,
      vol_total    DOUBLE
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
      facs_stain_id        INTEGER PRIMARY KEY DEFAULT nextval('seq_facs_stains'),
      sample_id            INTEGER NOT NULL REFERENCES samples(sample_id),
      assay_id             INTEGER NOT NULL REFERENCES assays(assay_id),
      vol_stained          DOUBLE,
      count_method         VARCHAR NOT NULL CHECK (count_method IN ('volumetric', 'bead')),
      vol_resuspended      DOUBLE,
      vol_measured         DOUBLE,
      bead_volume_added    DOUBLE,
      bead_concentration   DOUBLE,
      bead_population_path VARCHAR,
      stain_date           DATE,
      UNIQUE(sample_id, assay_id),
      CHECK (
        (count_method = 'volumetric'
          AND vol_resuspended IS NOT NULL AND vol_measured IS NOT NULL
          AND bead_volume_added IS NULL AND bead_concentration IS NULL)
        OR
        (count_method = 'bead'
          AND bead_volume_added IS NOT NULL AND bead_concentration IS NOT NULL
          AND vol_resuspended IS NULL AND vol_measured IS NULL)
      )
    )
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_facs_measurements START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS facs_measurements (
      facs_measurement_id  BIGINT PRIMARY KEY DEFAULT nextval('seq_facs_measurements'),
      facs_stain_id        INTEGER NOT NULL REFERENCES facs_stains(facs_stain_id),
      population_full_path VARCHAR NOT NULL,
      population           VARCHAR NOT NULL,
      metric               VARCHAR NOT NULL,
      value                DOUBLE,
      source_file          VARCHAR,
      imported_at          TIMESTAMP DEFAULT current_timestamp,
      UNIQUE(facs_stain_id, population_full_path, metric)
    )
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_elisa_measurements START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS elisa_measurements (
      elisa_measurement_id INTEGER PRIMARY KEY DEFAULT nextval('seq_elisa_measurements'),
      sample_id            INTEGER NOT NULL REFERENCES samples(sample_id),
      assay_id             INTEGER NOT NULL REFERENCES assays(assay_id),
      cytokine             VARCHAR NOT NULL,
      sample_id_raw        VARCHAR,
      replicate            INTEGER,
      value                DOUBLE,
      unit                 VARCHAR,
      result_status        VARCHAR,
      source_file          VARCHAR,
      imported_at          TIMESTAMP DEFAULT current_timestamp,
      UNIQUE(sample_id, assay_id, cytokine, sample_id_raw, replicate)
    )
  ")

  DBI::dbExecute(con, "CREATE SEQUENCE IF NOT EXISTS seq_histo_measurements START 1")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS histo_measurements (
      histo_measurement_id INTEGER PRIMARY KEY DEFAULT nextval('seq_histo_measurements'),
      sample_id            INTEGER NOT NULL REFERENCES samples(sample_id),
      assay_id             INTEGER REFERENCES assays(assay_id),
      metric               VARCHAR NOT NULL,
      value                DOUBLE,
      source_file          VARCHAR,
      imported_at          TIMESTAMP DEFAULT current_timestamp,
      UNIQUE(sample_id, assay_id, metric)
    )
  ")

  invisible(con)
}
