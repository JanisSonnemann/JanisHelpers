# Lab database schema redesign

Status: draft, pending review
Related: `docs/superpowers/specs/2026-07-17-facs-db-schema-diagram.html` (visual companion — entity cards, relationships, workflow diagram)

## Context

`R/DB-setup.R` and `R/export-to-sql.R` are personal, untracked working scripts that
set up and populate a local SQLite database (`data/MPOmRNA.db`) used to store FACS,
ELISA, and histology data across experiments for cross-experiment querying via
`dplyr`/`dbplyr`. Reviewing the current schema surfaced several structural problems:

- **Foreign keys are declared but never enforced.** SQLite requires
  `PRAGMA foreign_keys = ON` per connection; it's never set, so every
  `FOREIGN KEY` clause in the existing schema is decorative.
- **No `UNIQUE` constraints anywhere.** Nothing stops a re-run of an import
  script from appending duplicate rows — the scripts rely entirely on a manual
  "check via `dplyr::tbl() |> filter() |> collect()`, then append" pattern.
- **`facs.metric` is typed `REAL`** but holds category labels (`"Count"`,
  `"FractionOfParent"`, `"Median_CD4"`, ...) per `facs_import_wsp()`'s documented
  output shape — a modeling error papered over by SQLite's flexible typing.
- **No provenance.** A row can't be traced back to which `.wsp`/`.xlsx` file it
  came from.
- **No controlled vocabulary** for assay/tissue names — a typo (`"MPO-BCells"`)
  required a manual `DELETE` to fix, and a `rename_metric()` helper function has
  a copy-paste bug (checks one table, updates another).
- **`elisa` table shape doesn't match `elisa_read_results()`'s documented
  output** (raw per-replicate rows with `sample_id`/`replicate`/`unit`/
  `result_status`) — the table instead has `duration`/`timepoint` columns.

This document specifies a redesigned schema to replace it.

## Decision: DuckDB instead of SQLite or Postgres

The actual workload here is analytical (OLAP), not transactional: long-format
tidy tibbles, heavy `group_by()`/`summarize()` across experiments, joins on
`mouse_ID`/`tissue`/`experiment_ID`, single researcher, local machine, no
concurrent multi-user access. SQLite is a row-store built for many small
transactions, which fights this workload; a client-server RDBMS (Postgres)
solves a concurrency problem that doesn't exist here at the cost of running a
server. DuckDB is an embedded, single-file, columnar analytical database:

- Same `DBI`/`dplyr::tbl()` interface already in use — swapping the backend
  from `RSQLite::SQLite()` to `duckdb::duckdb()` doesn't require rewriting
  query code.
- Natively enforces `PRIMARY KEY`/`UNIQUE`/`FOREIGN KEY`/`CHECK`/`NOT NULL`
  constraints at insert time — no pragma required, unlike SQLite.
- Real `DATE`/`TIMESTAMP` types (the current schema stores dates as `TEXT` via
  `as.character()`).
- Faster on aggregation-heavy queries, which is most of what `analysis_*()`
  functions do.

Caveat: DuckDB's `FOREIGN KEY` constraints must be declared at `CREATE TABLE`
time (no `ALTER TABLE ADD FOREIGN KEY` later), and it's a newer engine than
SQLite/Postgres — acceptable for a local, single-writer analytical database.

## Schema

Two layers: **dimensions** (who/what/where of an experiment — stable, shared by
everything downstream) and **facts** (the measurements — one table per assay
domain, each shaped to match what that domain's existing import function
already produces as a tidy long tibble).

### Dimensions

```sql
CREATE SEQUENCE seq_domains START 1;
CREATE TABLE domains (
  domain_id   INTEGER PRIMARY KEY DEFAULT nextval('seq_domains'),
  domain_name VARCHAR NOT NULL UNIQUE
);
-- seed: ('facs'), ('histo'), ('elisa'); future modalities (e.g. 'scrnaseq')
-- are added here as a plain INSERT, not a schema migration.

CREATE SEQUENCE seq_experiments START 1;
CREATE TABLE experiments (
  experiment_id   INTEGER PRIMARY KEY DEFAULT nextval('seq_experiments'),
  experiment_code VARCHAR NOT NULL UNIQUE,   -- e.g. "25-7"
  experiment_name VARCHAR,
  project         VARCHAR,
  description     VARCHAR,
  created_at      TIMESTAMP DEFAULT current_timestamp
);

CREATE SEQUENCE seq_subjects START 1;
CREATE TABLE subjects (
  subject_id      INTEGER PRIMARY KEY DEFAULT nextval('seq_subjects'),
  mouse_id        VARCHAR NOT NULL UNIQUE,   -- e.g. "25-7-1"
  experiment_id   INTEGER NOT NULL REFERENCES experiments(experiment_id),
  mdc_id          VARCHAR,
  cage            VARCHAR,
  mouse_strain    VARCHAR,
  generation      VARCHAR,
  sex             VARCHAR,
  mouse_treatment VARCHAR,
  treatment_group VARCHAR,                  -- renamed from `group` (reserved word)
  dob             DATE,
  start_date      DATE,
  bmt_date        DATE,
  end_date        DATE,
  created_at      TIMESTAMP DEFAULT current_timestamp
);

CREATE SEQUENCE seq_samples START 1;
CREATE TABLE samples (
  sample_id    INTEGER PRIMARY KEY DEFAULT nextval('seq_samples'),
  subject_id   INTEGER NOT NULL REFERENCES subjects(subject_id),
  tissue       VARCHAR NOT NULL,   -- any specimen type: organ or fluid
                                    -- ("kidney", "serum", "urine", ...), free text
  collected_at DATE,
  created_at   TIMESTAMP DEFAULT current_timestamp,
  UNIQUE(subject_id, tissue)
);

-- 1:1 extension of `samples`, FACS-only: the master single-cell suspension
-- made from a digested tissue piece, before it's split into stained aliquots.
-- Fluid samples (serum, urine) have no row here.
CREATE TABLE sample_digestions (
  sample_id    INTEGER PRIMARY KEY REFERENCES samples(sample_id),
  total_weight DOUBLE,
  facs_weight  DOUBLE,
  vol_total    DOUBLE
);

CREATE SEQUENCE seq_assays START 1;
CREATE TABLE assays (
  assay_id    INTEGER PRIMARY KEY DEFAULT nextval('seq_assays'),
  assay_name  VARCHAR NOT NULL,
  domain_id   INTEGER NOT NULL REFERENCES domains(domain_id),
  description VARCHAR,
  UNIQUE(assay_name, domain_id)
);
```

### Facts

```sql
-- One row per (sample x panel) staining event -- i.e. one physical tube.
-- A single digested sample is routinely split into several differently
-- stained aliquots (e.g. 100ul of spleen suspension for the overview panel,
-- 200ul for the Treg panel), so vol_stained/vol_resuspended/vol_measured
-- belong here, one level below sample_digestions, not on the sample itself.
CREATE SEQUENCE seq_facs_stains START 1;
CREATE TABLE facs_stains (
  facs_stain_id        INTEGER PRIMARY KEY DEFAULT nextval('seq_facs_stains'),
  sample_id            INTEGER NOT NULL REFERENCES samples(sample_id),
  assay_id             INTEGER NOT NULL REFERENCES assays(assay_id),
  vol_stained          DOUBLE,             -- aliquot drawn from vol_total
  count_method         VARCHAR NOT NULL CHECK (count_method IN ('volumetric', 'bead')),
  vol_resuspended      DOUBLE,             -- volumetric method only
  vol_measured         DOUBLE,             -- volumetric method only: volume acquired
  bead_volume_added    DOUBLE,             -- bead method only: volume of bead reagent added
  bead_concentration   DOUBLE,             -- bead method only: beads per ul of reagent
  bead_population_path VARCHAR,           -- bead method only: which gated population is the bead cluster
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
);

CREATE SEQUENCE seq_facs_measurements START 1;
CREATE TABLE facs_measurements (
  facs_measurement_id  BIGINT PRIMARY KEY DEFAULT nextval('seq_facs_measurements'),
  facs_stain_id        INTEGER NOT NULL REFERENCES facs_stains(facs_stain_id),
  population_full_path VARCHAR NOT NULL,
  population           VARCHAR NOT NULL,
  metric               VARCHAR NOT NULL,   -- "Count" / "FractionOfParent" / "Median_CD4" / ...
  value                DOUBLE,
  source_file          VARCHAR,            -- provenance: which .wsp
  imported_at          TIMESTAMP DEFAULT current_timestamp,
  UNIQUE(facs_stain_id, population_full_path, metric)
);

CREATE SEQUENCE seq_elisa_measurements START 1;
CREATE TABLE elisa_measurements (
  elisa_measurement_id INTEGER PRIMARY KEY DEFAULT nextval('seq_elisa_measurements'),
  sample_id            INTEGER NOT NULL REFERENCES samples(sample_id),
  assay_id             INTEGER NOT NULL REFERENCES assays(assay_id),
  cytokine             VARCHAR NOT NULL,
  sample_id_raw        VARCHAR,            -- raw plate label, e.g. "25-7-1" -- not
                                            -- assumed to equal mouse_id, per
                                            -- elisa_read_results()'s own doc note
  replicate            INTEGER,
  value                DOUBLE,
  unit                 VARCHAR,
  result_status        VARCHAR,            -- "OK" / "OOR<" / "OOR>" / "<LLOQ" / ...
  source_file          VARCHAR,            -- provenance: which Results_<cytokine>.xlsx
  imported_at          TIMESTAMP DEFAULT current_timestamp,
  UNIQUE(sample_id, assay_id, cytokine, sample_id_raw, replicate)
);

CREATE SEQUENCE seq_histo_measurements START 1;
CREATE TABLE histo_measurements (
  histo_measurement_id INTEGER PRIMARY KEY DEFAULT nextval('seq_histo_measurements'),
  sample_id            INTEGER NOT NULL REFERENCES samples(sample_id),
  assay_id             INTEGER REFERENCES assays(assay_id),   -- nullable: only one
                                                               -- histo assay type in use today
  metric               VARCHAR NOT NULL,
  value                DOUBLE,
  source_file          VARCHAR,
  imported_at          TIMESTAMP DEFAULT current_timestamp,
  UNIQUE(sample_id, assay_id, metric)
);
```

Surrogate-key generation is shown above via `CREATE SEQUENCE` +
`DEFAULT nextval(...)` for portability across DuckDB versions; the
implementation plan may instead use `GENERATED ALWAYS AS IDENTITY` if the
targeted DuckDB version supports it. Either way the generation mechanism is an
implementation detail — the columns and constraints are what this document
specifies.

## Design decisions and rationale

| Decision | Why |
|---|---|
| Surrogate integer PKs, natural keys stay `UNIQUE` | Fast joins; `mouse_id`/`experiment_code` remain what any R-facing `db_` function exposes publicly, surrogate IDs stay internal. |
| `UNIQUE` constraint on every fact table's natural key | Re-running an import script is idempotent (DuckDB rejects the collision) instead of relying on a manual check-then-append. |
| `facs_stains` split from `sample_digestions` | One digested tissue sample is routinely split into several differently-stained aliquots, each with its own `vol_stained`/`vol_resuspended`/`vol_measured` (or bead equivalents) -- these can't live at the sample level. |
| `count_method` branch with a cross-field `CHECK` | Some panels are counted volumetrically (HTS: acquired-volume fraction), others via counting beads (bead-to-cell ratio x reagent concentration) -- structurally different fields, and the `CHECK` catches a half-filled-in row at insert time. |
| `domains` lookup table instead of a hard-coded `CHECK` list on `assays.domain` | DuckDB can't cheaply `ALTER` a `CHECK` constraint later. Adding a future modality (e.g. scRNA-seq) becomes a plain `INSERT`, not a schema migration. |
| `source_file` + `imported_at` on every fact table | Every measurement traces back to the raw export it came from. |
| `samples.tissue` is free `VARCHAR`, not an enum | Specimen types span organs and fluids (serum, urine) -- forcing a rigid vocabulary here would fight real lab variation; typo protection instead lives on the `assays` lookup, where it matters more. |

## Workflow

Every domain follows the same pipeline (see the diagram doc for the visual):

1. A raw export (`.wsp`, `Results_<cytokine>.xlsx`, meta workbook) is parsed by
   an existing package function (`facs_import_wsp()`, `elisa_read_results()`,
   `meta_read()`/`meta_clean()`) into a tidy long tibble, in-session.
2. The tibble is written into the matching fact table, behind that table's
   `UNIQUE` constraint.
3. Cross-experiment analysis reads back out through the same
   `dplyr::tbl(con, ...) |> filter() |> collect()` interface already in use,
   feeding into `analysis_*()`/`plot_*()`/`report_knit_*()`.

Onboarding a new experiment follows the dimension chain top-down: register the
experiment, register its subjects, register their samples (and
`sample_digestions` for FACS-processed tissue), then import each domain's
measurements against those samples.

## Out of scope

- The R-facing `db_` functions that will read/write this schema (connection
  management, per-domain import/write functions, query helpers that hide the
  surrogate keys behind `mouse_id`/`tissue`/`experiment_code`). Separate design.
- Migrating the existing `MPOmRNA.db` (SQLite) data into this schema. Separate
  task, follows once the schema and the `db_` functions both exist.
- Extending to additional modalities (e.g. scRNA-seq) beyond confirming the
  schema's `domains`-lookup extension point supports them structurally.
