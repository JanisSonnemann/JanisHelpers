# DB schema changes for the FACS export: drop calc inputs, add `wsp_root` and `import_log`

Status: draft, pending review
Parent: `2026-11-12-facs-pipeline-workflow-design.md` (Stage 0, db-setup step 4b)
Branch: `worktree-facs-pipeline-workflow`

## Why

The FACS pipeline design moved all calculations (counts per g, % of
reference populations, ratios) into a per-panel export script that runs
*before* the database. That means:

- Stain volumes, bead settings and organ weights are calculation **inputs**.
  They live in the meta workbook and are never stored in the database.
- The database needs to know which **workspace root** each piece of FACS
  data came from, so a root can be found and replaced as one unit.
- Every import needs a permanent **log row**, kept even after the data it
  describes is overwritten.

These are schema changes, so they have to land **before the first real
import** (db-setup step 6). No database with real data exists yet, so there
is no migration: the `CREATE TABLE` statements are just changed.

## Goals

- `facs_stains` and `db_write_facs()` carry no volume/bead/count-method
  data. `sample_digestions` is gone.
- Every `facs_stains` row records its `wsp_root`.
- The database can never hold two roots for the same experiment x panel,
  and one root never spans two experiments or two panels.
- One `import_log` table records every import of every data type (FACS,
  ELISA, histology).
- ELISA and histology imports write their log row now. FACS logging comes
  with the Stage 3 export function.

## Non-goals

- The Stage 3 export function itself (`db_export_facs()` or similar), the
  overwrite prompt, and the delete-and-replace transaction. 4b only builds
  what they need.
- Overwrite / replace for ELISA or histology. Re-importing still skips
  existing rows (`ON CONFLICT DO NOTHING`), as today.
- A QC gate for ELISA or histology (decided by the user: not needed for now).
- Migrating any existing database file.
- A query function for the log. `DBI::dbReadTable(con, "import_log")` is
  enough for now.
- Switching `db_write_facs()` from the old CamelCase input
  (`PopulationFullPath`, from the deprecated `facs_import_wsp()`) to the
  snake_case `facs_read_wsp()` shape. The Stage 3 design decides the
  export's input shape. Noted here so nobody is surprised.

## Anti-goals

- Two versions of the same experiment x panel in the database.
- A log row that can be lost or that disagrees with the data (data written,
  log not, or the reverse).
- More ceremony at import time than the pipeline design already accepted.

## Facts this design rests on (from the user)

- **One workspace per panel.** Within an experiment, each `.wsp` holds
  exactly one panel. So a root (`treg_26-2`) = one experiment x one panel,
  covering every mouse and tissue stained with it.
- Workspace filenames are `<root>_<YYYYMMDD>.wsp` (pipeline design,
  decided).
- Every data type gets logged.
- ELISA and histology have no QC check for now.

## Schema changes (in `create_schema_()`, `R/db_connect.R`)

### `facs_stains`

Drop: `vol_stained`, `count_method`, `vol_resuspended`, `vol_measured`,
`bead_volume_added`, `bead_concentration`, `bead_population_path`, and the
table's multi-column `CHECK`.

Add: `wsp_root VARCHAR NOT NULL`.

Keep: `facs_stain_id`, `sample_id`, `assay_id`, `stain_date`,
`UNIQUE(sample_id, assay_id)`.

```sql
CREATE TABLE IF NOT EXISTS facs_stains (
  facs_stain_id INTEGER PRIMARY KEY DEFAULT nextval('seq_facs_stains'),
  sample_id     INTEGER NOT NULL REFERENCES samples(sample_id),
  assay_id      INTEGER NOT NULL REFERENCES assays(assay_id),
  wsp_root      VARCHAR NOT NULL,
  stain_date    DATE,
  UNIQUE(sample_id, assay_id)
)
```

*(My call)* The root string is stored instead of being derived from
experiment + panel: it is the name used in filenames, the export call and
the log, so a direct lookup is clearer than a join.

### `facs_measurements`: drop the foreign key to `facs_stains` *(my call)*

`facs_stain_id` stays (`INTEGER NOT NULL`) but loses
`REFERENCES facs_stains(facs_stain_id)`. Reason: DuckDB refuses to delete a
parent row in the same transaction that deleted its children (verified:
`Constraint Error: Violates foreign key constraint because key ... is still
referenced`). Stage 3 must delete a root's measurements **and** stains and
insert the new ones in one all-or-nothing transaction, which is impossible
with the foreign key in place. Integrity is kept by code instead: the only
writer (`db_write_facs()`, later the Stage 3 export) takes `facs_stain_id`
from the stain row it just found or created. A test asserts that every
`facs_measurements.facs_stain_id` has a matching stain after each write
path. Considered and rejected: keeping stain rows across an overwrite
(leaves orphan stains for samples dropped from the new workspace) and
deleting outside the transaction (breaks the all-or-nothing anti-goal).

```sql
CREATE TABLE IF NOT EXISTS facs_measurements (
  facs_measurement_id  BIGINT PRIMARY KEY DEFAULT nextval('seq_facs_measurements'),
  facs_stain_id        INTEGER NOT NULL,  -- no FK: see design 2026-11-12
  population_full_path VARCHAR NOT NULL,
  population           VARCHAR NOT NULL,
  metric               VARCHAR NOT NULL,
  value                DOUBLE,
  source_file          VARCHAR,
  imported_at          TIMESTAMP DEFAULT current_timestamp,
  UNIQUE(facs_stain_id, population_full_path, metric)
)
```

### `sample_digestions`

Removed from the schema entirely.

### `import_log` (new)

```sql
CREATE SEQUENCE IF NOT EXISTS seq_import_log START 1;
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
  supersedes_import_id INTEGER REFERENCES import_log(import_id),
  imported_at          TIMESTAMP NOT NULL DEFAULT current_timestamp
)
```

| Column | FACS (Stage 3) | ELISA / histology (4b) |
|---|---|---|
| `domain_id` | `facs` | `elisa` / `histo` |
| `experiment_code` | the root's experiment | distinct experiment code(s) of the rows written, joined with `", "` if a plate spans experiments *(my call)* |
| `import_unit` | workspace root, e.g. `treg_26-2` | `basename(source_file)` *(my call: ELISA/histo have no root concept)* |
| `source_file` | `.wsp` path relative to `raw/` | file path relative to `raw/` |
| `source_version` | date parsed from `<root>_<YYYYMMDD>.wsp` | `NULL` |
| `qc_report` | QC report path | `NULL` |
| `qc_passed` | `TRUE` (export refuses otherwise) | `NULL` |
| `rows_written` | rows inserted | rows inserted (may be 0 on a re-run) |
| `package_version` | `as.character(utils::packageVersion("JanisHelpers"))` | same |
| `supersedes_import_id` | previous log row of the same domain + `import_unit`, when overwriting | `NULL` |

`domain_id` references the existing `domains` table, so it matches the rest
of the schema (assays do the same).

Log rows are never updated or deleted by package code. *(My call: a
re-import that writes 0 rows still gets a log row -- "log every import"
taken literally, and it shows the import was attempted.)*

### Old-schema guard *(my call)*

`db_connect()` uses `CREATE TABLE IF NOT EXISTS`, so opening a database file
created by the old schema would silently keep the old tables. After schema
creation, `db_connect()` checks that `facs_stains` has a `wsp_root` column
and that `import_log` exists; if not, it errors saying the file was created
by an older schema and must be recreated. No real database exists yet, so
this only protects test or scratch files.

## Function changes

### `db_write_facs()`

New signature:

```r
db_write_facs(con, data, mouse_id, tissue, assay_name, wsp_root,
              source_file, stain_date = NA)
```

- Removed: `count_method`, `vol_stained`, `vol_resuspended`,
  `vol_measured`, `bead_volume_added`, `bead_concentration`,
  `bead_population_path`.
- `wsp_root`: required, single non-empty string, validated like
  `source_file` (done in step 3).
- **Root / filename consistency** *(my call)*: `basename(source_file)` must
  match `^<wsp_root>_[0-9]{8}\.wsp$`, otherwise error. It's cheap and catches
  passing the wrong root for a file. The user already decided this naming
  convention, so it adds no new ceremony.
- **One root per experiment x panel, one experiment x panel per root.**
  Before inserting the stain, error (writing nothing) if either:
  - a `facs_stains` row exists for the same experiment (via sample ->
    subject -> experiment) and `assay_id` with a **different** `wsp_root`
    -- the message names both roots and says to export the existing root
    instead (Stage 3 replaces it) or fix the name; or
  - a `facs_stains` row exists with the **same** `wsp_root` but a
    different experiment or `assay_id`.
- **One workspace version per root** *(my call)*: error (writing nothing)
  if the stain already has measurements whose `source_file` differs from
  this call's -- e.g. `treg_26-2_20261120.wsp` arriving on top of
  `treg_26-2_20261112.wsp`. Otherwise `ON CONFLICT DO NOTHING` would keep the
  old values and add only new rows, mixing two versions under one root. The
  message says replacing a root is the Stage 3 export's job. No new column:
  `facs_measurements.source_file` already records the version.
- Still idempotent: the same sample x panel x root x `source_file` writes no
  new stain and `ON CONFLICT DO NOTHING` on measurements.
- Does **not** write a log row (it is called per sample; the Stage 3 export
  writes one row per workspace).

### `db_write_elisa()` / `db_write_histo()`

- **`db_write_histo()` becomes a whole-table writer, like `db_write_elisa()`**
  *(my call)*: new signature `db_write_histo(con, data, assay_name = NA,
  source_file)`, with `mouse_id` / `tissue` as columns in `data`, resolved
  per row via `lookup_sample_ids_()` (errors listing every unresolved
  pair, writes nothing). This keeps "one log row per call" equal to "one log
  row per imported file" for both domains. Breaking change; no real data
  exists yet. `assay_name` stays optional as today.
- `source_file` becomes **required** (single non-empty string, path
  relative to `raw/`). *(My call: the log row needs it, and FACS already
  requires it.)*
- Each call writes its data **and** one `import_log` row inside a single
  `DBI::dbWithTransaction()`. If either fails, neither is written.
- No `qc_passed` argument.

### Internal helper `write_import_log_()`

In `R/db_write_measurements.R` (or a new `R/db_import_log.R` -- builder's
choice), unexported. Takes `con`, `domain`, `experiment_code`,
`import_unit`, `source_file`, `rows_written`, and optional
`source_version`, `qc_report`, `qc_passed`, `supersedes_import_id`. Fills
`package_version` itself. Returns the new `import_id`. Must be called
inside the caller's transaction (it does not open its own). Stage 3 reuses
it.

### `db_write_samples()`

Drop the `total_weight`, `facs_weight`, `vol_total` columns and the
`sample_digestions` insert. Extra columns in `data` are ignored, as today.

### `db_query_facs()`

Drop the removed `facs_stains` columns from the output. Add `wsp_root`.

## Docs and dependencies

- `utils` added to `DESCRIPTION` Imports and to CLAUDE.md's Imports line
  (for `utils::packageVersion()`).
- Update `docs/superpowers/specs/2026-07-17-facs-db-schema-design.md` to the
  new DDL: remove `sample_digestions` and the volume columns, add `wsp_root`
  and `import_log`, with a dated note pointing here.
- CLAUDE.md `db_` line: mention `import_log`.
- Roxygen + `devtools::document()` for every changed function.
- db-setup plan: tick step 4b.

## Tests (TDD, `tests/testthat/`)

- Schema: `facs_stains` has `wsp_root` and none of the dropped columns;
  `sample_digestions` doesn't exist; `import_log` exists with the columns
  above.
- Old-schema guard: a file whose `facs_stains` lacks `wsp_root` makes
  `db_connect()` error.
- `db_write_facs()`: writes `wsp_root`; rejects a missing/`NA`/empty
  `wsp_root`; rejects a `source_file` whose basename doesn't match the root;
  rejects a second root for the same experiment x panel; rejects one root
  across two panels or two experiments; rejects a second dated version of
  the same root; all rejections leave the tables unchanged; still
  idempotent; every measurement's `facs_stain_id` has a matching stain.
- Overwrite feasibility (guards the no-FK decision for Stage 3): deleting a
  root's measurements and stains and re-inserting them in one
  `DBI::dbWithTransaction()` succeeds.
- `db_write_histo()`: whole table across several mice, unresolved pairs
  listed, as for ELISA.
- `db_write_elisa()` / `db_write_histo()`: one log row per call with the
  values in the table above (incl. `qc_passed` `NULL`, `package_version`
  filled); a re-run writes a second log row with `rows_written = 0`; a
  failing data insert leaves no log row; missing `source_file` errors.
- `db_write_samples()`: works without weight columns; ignores them if given.
- `db_query_facs()`: returns `wsp_root`, not the dropped columns.
- Existing tests updated to the new signatures.
- `devtools::check()`: 0 errors, 0 warnings.

## Left to the builder

- File placement of the helper, exact error wording, test fixture helpers.
- How the old-schema check is queried (`information_schema` or
  `PRAGMA table_info`).
