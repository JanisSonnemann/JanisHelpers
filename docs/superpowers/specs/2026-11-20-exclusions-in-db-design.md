# Full record in the DB: dropouts and excluded data stored, exclusions flagged live

Date: 2026-11-20
Status: draft, awaiting approval

## Intent

Today the lab DuckDB leaves out:

- **Dropouts without a `mouse_ID`.** `db_write_meta()` filters them out because
  `subjects.mouse_id` is `NOT NULL UNIQUE`. In 26-1 this is 6 of 7 early mice.
- **The meta sheet's `early` / `reason` columns.** They are not stored.
- **Excluded data.** `db_export_facs()` / `db_export_elisa()` / `db_export_histo()`
  drop the files and samples matched by the workbook's `exclusions` sheet before
  writing anything (`R/db_export.R:88`, `R/db_export_table.R:137`). Only the
  rule survives, in `import_exclusions`.

The user wants the DB to be a **full record**, because leaving data out as if it
never happened is wrong. They also want to:

1. analyse dropouts for patterns (group, cage, sex, day, reason);
2. re-check analyses with excluded mice or samples put back in;
3. change their mind about an exclusion, ideally without re-importing data.

Analysis still drops excluded data by default.

## Goals

- Every enrolled mouse is in `subjects`, including dropouts without a
  `mouse_ID`, with `early` and the dropout reason stored.
- Every measurement is written. Exports no longer drop excluded rows.
- Exclusion rules live in their own table, per experiment, separate from the
  data. Whether a row is excluded is computed when it is read.
- Changing an exclusion means editing the workbook sheet and making one call.
  No prep, QC report or export is re-run.
- Every change to the rules is kept. The DB shows which rules applied when, and
  why.
- `db_query_*()` return `excluded` / `exclusion_reason` columns.
  `db_load_experiment()` drops excluded rows by default and keeps them, flagged,
  on request.

## Non-goals

- No plot helper highlights excluded points. With `keep_excluded = TRUE` the
  `excluded` column is there to map or filter by hand. This can be added later.
- No change to the QC reports. They already show all samples and flag
  exclusions.
- No storing of volumes / weights (unchanged policy).

## Anti-goals (would count as failure even if it technically works)

- An excluded value reaches an analysis by default. The default load must give
  the same data as today's DB for the same exclusions sheet.
- Two places say different things about whether a row is excluded. There is one
  source of truth: the current rule set.
- Changing an exclusion requires a re-export.
- An exclusion change leaves no trace.

## Approaches considered

- **A. Flag written into each row at export.** `excluded` / `exclusion_reason`
  columns on stains and ELISA / histology rows. Queries are the simplest, but
  every change of mind needs a fresh prep, a QC report and a re-export.
  Rejected: it fails goal 3.
- **B. Rules table, flag computed on read (chosen by the user).** One current
  rule set per experiment, versioned. `db_query_*()` match rows against it. A
  change is one `db_write_exclusions()` call. Cost: the matching moves into SQL,
  and FACS must store `file_name` per stain to match file-level rules.
- **C. Rules table plus materialised flags, refreshed on every rule change.**
  Fast queries, but two copies that can disagree, and the DB is small.
  Rejected: needs a sync step for no real gain.

## Design

### Subjects: everyone, with dropout info

- `db_write_meta()` registers every row of the `meta` sheet that has a
  `mouse_ID` **or** an `ID` (facility ID, stored as `mdc_id`). A row with
  neither is an error that names the row.
- **Key for ID-less mice (builder's call, flagged for review):** `mouse_id` =
  `mouse_ID` if present, otherwise the facility `ID` (e.g. `LFS0-21559`). It is
  unique and traceable, and keeps `mouse_id NOT NULL UNIQUE` and every existing
  join. In 26-1 every dropout has an `ID`, and `ID` values are unique.
  - When a dropout later gets a `mouse_ID` (does not happen in practice), the
    re-registration should be refused with a clear message rather than creating
    a second subject.
- New `subjects` columns:
  - `early BOOLEAN`: sheet `early` is `"yes"` or `TRUE`, otherwise `FALSE`.
  - `early_reason VARCHAR`: the sheet's `reason` column. Its name avoids a clash
    with `exclusion_reason`.
  - `id_source VARCHAR NOT NULL DEFAULT 'mouse_ID'`: `'mouse_ID'` or
    `'facility'` (key taken from `ID` because `mouse_ID` was empty). Migrated
    rows get the default, which is correct: only mice with a `mouse_ID` were
    stored before.
- `db_write_subjects()` stays idempotent. Re-running `db_write_meta()` updates
  `early` / `early_reason` on existing subjects (today it only adds new ones).
  **`end_date` changes are refused (builder's call):** if a stored subject's
  `end_date` differs from the sheet, `db_write_meta()` errors and lists the
  mice with both dates. End-point samples store `end_date` as their
  `collected_at`, so changing it silently would break the link. A deliberate
  correction is a separate, later function. Today such changes are silently
  ignored.
- `db_query_subjects()` returns `early`, `early_reason` and
  `excluded_mouse_wide` / `exclusion_reason` (from the current rule set, see
  below).

### Measurements: everything written

- `db_export_facs()`, `db_export_elisa()` and `db_export_histo()` write all rows
  of the prep. The "every file / sample is excluded" errors go away.
- `facs_stains` gets a `file_name VARCHAR` column (the FCS file name from the
  prep), so file-level rules can match.
- The `exclusions` argument of the exports is still validated (unknown panels,
  rows matching nothing in this unit -> message, see below). It no longer
  removes data.
- `import_log` and the existing per-import `import_exclusions` rows: the exports
  stop writing `import_exclusions` (the rule set replaces it). The table and its
  old rows are kept as history and not read by the loader any more.

### Exclusion rules: one versioned set per experiment

New tables (no self-referencing FKs, so `EXPORT DATABASE` snapshots keep working):

```
exclusion_sets
  exclusion_set_id  INTEGER PK (sequence)
  experiment_id     INTEGER NOT NULL REFERENCES experiments
  source_file       VARCHAR          -- workbook path relative to data/
  written_at        TIMESTAMP NOT NULL DEFAULT current_timestamp
  package_version   VARCHAR NOT NULL
  note              VARCHAR          -- optional free text ("re-included 26-1-17")

exclusion_rules
  exclusion_rule_id INTEGER PK (sequence)
  exclusion_set_id  INTEGER NOT NULL REFERENCES exclusion_sets
  mouse_id          VARCHAR NOT NULL
  tissue            VARCHAR          -- NULL = all
  panel             VARCHAR          -- assay_name; NULL = all units
  condition         VARCHAR          -- NULL = all
  file_name         VARCHAR          -- FACS only; NULL = all
  collected_at      DATE             -- NULL = all dates
  reason            VARCHAR NOT NULL
  decided_on        DATE
```

- The **current** set of an experiment is the one with the highest
  `exclusion_set_id`. An empty sheet writes an empty set, which re-includes
  everything and is itself recorded.
- **`db_write_exclusions(con, exclusions, experiment_code, source_file = NA,
  note = NA)`** (new, exported, `db_` domain, `R/db_exclude.R`):
  - validates the sheet like today: `mouse_ID` + `reason` required, and every
    `mouse_ID` must be a registered subject of the experiment (error);
  - a `panel` that isn't registered yet gets a **message, not an error**
    (builder's call). Assays are registered only at export, and every unit's
    `final = TRUE` run writes the whole sheet, so rows for units not yet
    exported are normal. This replaces `check_exclusion_panels_()`'s warning
    in the exports. Each export then reports rows naming its own panel that
    match none of its rows, so typos still surface;
  - writes a new set only when the rules differ from the current set (otherwise
    a message: "unchanged"), in one transaction;
  - reports rules that currently match no stored row as a **message, not an
    error**, because data for that unit may not be imported yet;
  - returns the new `exclusion_set_id` invisibly.
- `facs_prep_run()` / `elisa_prep_run()` / `histo_prep_run()` with
  `final = TRUE` call `db_write_exclusions()` with the same `exclusions` sheet
  after `db_write_meta()`, so the normal workflow stays one call. Changing your
  mind later = edit sheet + `db_write_exclusions()` on its own.
  - `experiment_code` is inferred from the prep's mice, which must all belong
    to one registered experiment (error otherwise). The templates and existing
    panel scripts don't change (builder's call).
  - New optional argument `meta_path` on the three `*_prep_run()` functions
    fills the set's `source_file`. The templates pass the workbook path they
    already read; existing scripts get `NA`.

### Matching on read

- A DuckDB **view per domain**, created in `create_schema_()`, matches each
  measurement's stain / sample against the current rule set of its experiment:
  - `mouse_id` equal;
  - each of `tissue`, `panel`, `condition`, `file_name`, `collected_at` either
    `NULL` in the rule or equal to the row's value;
  - `panel` is compared with `assay_name`;
  - `condition` and `file_name` only exist for FACS, so ELISA / histology rules
    with a non-`NULL` `condition` / `file_name` never match there (same as
    today's matcher).
- Views: `facs_exclusion_flags(facs_stain_id, excluded, exclusion_reason)`,
  `elisa_exclusion_flags(elisa_measurement_id, ...)`,
  `histo_exclusion_flags(histo_measurement_id, ...)`. Several matching rules ->
  reasons joined with `"; "`, distinct.
- `db_query_facs()` / `db_query_elisa()` / `db_query_histo()` left-join their
  view and return `excluded` (`FALSE` when nothing matches) and
  `exclusion_reason`. They stay lazy.
- **Equivalence requirement:** the view must flag the rows that
  `facs_match_exclusions()` / `match_sample_exclusions_()` would have dropped
  for the same sheet. A test compares both on the fixtures.
  - **One intended difference (builder's call):** dated rules match the
    *stored* `collected_at`. End-point samples store the subject's `end_date`,
    so a rule dated on the death day now also hits that mouse's end-point
    sample. Today the R matcher compares against the sheet's empty date and
    errors with "matches no sample" (`exclusion_date_hint_()`). "Exclude what
    was collected that day" is the natural reading, and the DB can't tell the
    two apart anyway.
  - The equivalence test covers this case explicitly.
  - `match_sample_exclusions_()` keeps running before export, for the QC report
    and messages. Its date hint changes to say that end-point samples match
    their death date in the DB.
- `excluded_mouse_wide` (for subjects) = the current set holds a rule for the
  mouse with every scope column `NULL`. Under one set per experiment this
  replaces today's "mouse-wide in every import unit"
  (`fully_excluded_mice_()`).

### Analysis

- `db_load_experiment(..., keep_excluded = FALSE)`:
  - every data tibble drops rows with `excluded = TRUE`, and drops the
    `excluded` / `exclusion_reason` columns;
  - `meta` drops subjects with `id_source = 'facility'` (dropouts without a
    `mouse_ID`) and mice with `excluded_mouse_wide`. That is exactly today's
    `dat$meta`: every mouse with a `mouse_ID`, minus mouse-wide exclusions.
    Dropouts that have a `mouse_ID` stay, as today;
  - `exclusions` = the current rule set (same columns as today, minus
    `import_unit`, plus `exclusion_set_id`).
- `keep_excluded = TRUE`: every tibble keeps all rows with `excluded` /
  `exclusion_reason`, and `meta` holds every subject with `early`,
  `early_reason`, `excluded_mouse_wide`.
- **New element `dat$cohort`** (always, regardless of `keep_excluded`): every
  subject of the experiment with all subject columns plus `early`,
  `early_reason`, `excluded_mouse_wide`, `exclusion_reason`. This is the input
  for dropout pattern analysis and the cohort tables. `cohort` becomes a
  reserved name next to `meta` / `exclusions`.
- `meta_summary_table()` / `meta_mouse_table()`: `workbook = NULL` uses
  `dat$cohort` without the "mice missing" message. The workbook keeps working.
  "Excluded" = `excluded_mouse_wide`, and the `not in DB` status disappears for
  DB input.
- `analysis-db` template: the overview chunk uses `dat$cohort`, and
  `dat$workbook` becomes optional.

### Existing database files

- Additive migration in `create_schema_()`, like `population_order`:
  - `ALTER TABLE ... ADD COLUMN IF NOT EXISTS` for `subjects.early`,
    `subjects.early_reason`, `subjects.id_source`, `facs_stains.file_name`;
  - `CREATE TABLE IF NOT EXISTS` for the two rule tables;
  - `CREATE OR REPLACE VIEW` for the flag views.
  - No `db_old_schema_()` refusal.
- Until the user re-runs each panel / assay script with `final = TRUE`:
  - old units lack the excluded data;
  - old FACS stains have `file_name = NULL`, so file-level rules can't match
    them, and `db_write_exclusions()` says so in its message;
  - an experiment with no rule set has nothing excluded, and `db_load_experiment()`
    messages that no exclusions are recorded.
- Migration path for the user: re-run `db_write_meta()` (dropouts + early), then
  each unit's script with `final = TRUE` (full data + first rule set).

### Testing

- `test-db_exclude.R`:
  - set versioning (unchanged sheet writes nothing; an empty sheet writes an
    empty set);
  - validation errors;
  - unmatched-rule message.
- Equivalence of view flags vs the R matchers on `tests/fixtures` (FACS
  file-level, tissue-level, condition-level; ELISA dated vs end-point).
- `db_write_meta()` with an ID-less dropout; refusal of a row with neither ID.
- Exports write excluded rows; the "everything excluded" path no longer errors.
- `db_load_experiment()` default output equals today's for the same sheet;
  `keep_excluded = TRUE` adds rows and flags; `dat$cohort` holds dropouts.
- Migration: open a file built by the old schema, add a rule set, query.
- Snapshot (`db_snapshot()`) still works with the new tables and views.

### Docs

- `CLAUDE.md` data-structures section and the schema spec get a note pointing
  to this spec.
- Roxygen for every touched exported function, then `devtools::document()`.
- `devtools::check()`: 0 errors, 0 warnings.

## Left to the builder

- SQL shape of the views (one view per domain, or a single macro / CTE reused
  by all three).
- How `early` strings from the sheet map to `TRUE` (`"yes"`, `"ja"`, `TRUE`,
  `"x"`?). Be lenient, and error on unknown non-empty values.
- Where the unmatched-rule check lives (in `db_write_exclusions()` and,
  narrowed to the unit, in each export).
