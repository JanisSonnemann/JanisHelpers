# Full Record in the DB: Exclusions Flagged Live -- Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Store every mouse (including dropouts without a `mouse_ID`) and every
measurement in the lab DuckDB. Exclusions become a versioned per-experiment
rule set that `db_query_*()` match on read, and `db_load_experiment()` still
drops excluded data by default.

**Architecture:**
- New `exclusion_sets` / `exclusion_rules` tables, written by
  `db_write_exclusions()`.
- Three SQL views (`facs_exclusion_flags`, `elisa_exclusion_flags`,
  `histo_exclusion_flags`) compute `excluded` / `exclusion_reason` per
  measurement from the experiment's latest set.
- The exports stop dropping rows. `subjects` gains `early`, `early_reason` and
  `id_source`.
- Existing files are migrated additively in `create_schema_()`.

**Tech Stack:** R package (testthat 3e), DBI + duckdb, dplyr/dbplyr, roxygen2.

**Spec:** `docs/superpowers/specs/2026-11-20-exclusions-in-db-design.md` (read it first).

## Global Constraints

- Follow `CLAUDE.md`:
  - base pipe `|>`;
  - every external call `pkg::fn()`;
  - `dplyr::select(!col)`;
  - non-ASCII escaped as `\uXXXX` in R strings and plain ASCII in roxygen;
  - every exported function has `@param` / `@returns` / `@export`.
- No self-referencing foreign keys (DuckDB `EXPORT DATABASE` / `db_snapshot()`
  rejects them).
- Migration is additive only:
  - `ALTER TABLE ... ADD COLUMN IF NOT EXISTS`;
  - `CREATE TABLE IF NOT EXISTS`;
  - `CREATE OR REPLACE VIEW`;
  - no new `db_old_schema_()` reason.
- Default `db_load_experiment()` output must equal today's for the same
  exclusions sheet (anti-goal in the spec).
- Run tests with the environment below; `devtools::check()` at the end:
  0 errors, 0 warnings, only the NOTEs listed in `CLAUDE.md`.

Test environment (once per shell; `writexl` is only in the renv cache):

```bash
C=$HOME/Library/Caches/org.R-project.R/R/renv/cache/v5/macos/R-4.6/aarch64-apple-darwin23
S=$TMPDIR/rlib; mkdir -p $S; ln -sf "$(ls -d $C/writexl/*/*/writexl | tail -1)" $S/writexl
L=$HOME/Library/Caches/org.R-project.R/R/renv/library
export R_LIBS="$S:$L/facs-pipeline-workflow-0ff3068e/macos/R-4.6/aarch64-apple-darwin23:$L/JanisHelpers-b26afe2b/macos/R-4.6/aarch64-apple-darwin23"
# then: Rscript --vanilla -e 'devtools::test(filter = "<file stem>")'
```

## Review Focus

1. **A mouse matched by several rules** (e.g. a mouse-wide rule plus a tissue
   rule) must still give exactly one row per measurement, with reasons joined
   by `"; "`. A naive join duplicates rows. Test in Task 5.
2. **Excel empty cells arrive as `""` or `NA`.** Both mean "all" for scope
   columns, and both must give identical stored rules and flags. They must not
   create a "changed" set. Test in Task 3.
3. **An experiment with no rule set yet, or with an empty set,** gives
   `excluded = FALSE`, never `NA`, so `dplyr::filter(!excluded)` keeps
   everything. Test in Task 5.
4. **Rules of experiment A never flag rows of experiment B**, even with
   identical tissue / panel names. Test in Task 5.
5. **Old FACS stains with `file_name = NULL` and a file-level rule:** not
   flagged, and `db_write_exclusions()` names the rule in its "matches no
   stored row" message. Test in Task 3 (message) and Task 5 (not flagged).

---

### Task 1: Schema -- new columns, rule tables, migration

**Files:**
- Modify: `R/db_connect.R` (`create_schema_()`, roxygen of `db_connect()`)
- Test: `tests/testthat/test-db_connect.R`, `tests/testthat/test-db_snapshot.R`

**Interfaces:**
- Produces:
  - `subjects.early BOOLEAN NOT NULL DEFAULT FALSE`
  - `subjects.early_reason VARCHAR`
  - `subjects.id_source VARCHAR NOT NULL DEFAULT 'mouse_ID'`
  - `facs_stains.file_name VARCHAR`
  - tables `exclusion_sets` and `exclusion_rules`, columns exactly as in the
    spec section "Exclusion rules", with sequences `seq_exclusion_sets` and
    `seq_exclusion_rules`.
  - Views come in Task 5.

- [ ] **Step 1: Write the failing tests** in `test-db_connect.R`:
  - `test_that("new schema has exclusion tables and subject flags")`:
    `local_test_db()`; expect `DBI::dbExistsTable` for both tables; expect
    `c("early", "early_reason", "id_source")` in
    `DBI::dbListFields(con, "subjects")`; expect `"file_name"` in the
    `facs_stains` fields.
  - `test_that("files from the previous schema are migrated")`:
    1. build a file with the current `create_schema_()`;
    2. `ALTER TABLE ... DROP COLUMN` the four new columns and `DROP TABLE` the
       two new tables;
    3. insert one subject;
    4. disconnect, reopen with `db_connect()`;
    5. expect the columns back, the subject's `early == FALSE` and
       `id_source == "mouse_ID"`.
- [ ] **Step 2: Run** `devtools::test(filter = "db_connect")`. Expect FAIL
  (tables / fields missing).
- [ ] **Step 3: Implement.** In `create_schema_()`:
  - add the columns to the `CREATE TABLE` statements;
  - add the `ALTER TABLE ... ADD COLUMN IF NOT EXISTS` lines next to the
    existing `population_order` migration;
  - create the two tables after `import_exclusions`;
  - add a comment there saying `import_exclusions` is history only (no longer
    written, see the spec).
- [ ] **Step 4: Snapshot test.** In `test-db_snapshot.R` add
  `test_that("snapshot works with exclusion tables")`: insert one set and one
  rule by SQL, run `db_snapshot()`, expect no error and that
  `exclusion_rules.csv` (or `.parquet`, whatever the existing snapshot test
  checks) exists.
- [ ] **Step 5: Run** `devtools::test(filter = "db_connect|db_snapshot")`.
  Expect PASS.
- [ ] **Step 6: Commit**: `feat(db): exclusion rule tables and subject dropout columns`.

### Task 2: Subjects -- every mouse, early flags, end-date guard

**Files:**
- Modify: `R/db_write_dimensions.R` (`db_write_subjects()`, `db_write_meta()`)
- Modify: `R/db_query.R` (`db_query_subjects()`)
- Test: `tests/testthat/test-db_write_dimensions.R`, `tests/testthat/test-db_query.R`

**Interfaces:**
- Consumes: Task 1 columns.
- Produces:
  - `db_write_subjects(con, data)` accepts optional `early` (logical),
    `early_reason`, `id_source`.
  - `db_query_subjects(con)` additionally returns `mdc_id`, `early`,
    `early_reason`, `id_source`. Task 5 adds `excluded_mouse_wide` /
    `exclusion_reason`.

- [ ] **Step 1: Write the failing tests**:
  - `db_write_meta registers dropouts without mouse_ID by facility ID`: a
    `meta` sheet tibble with the template columns (`mouse_ID`, `id`, `cage`,
    `line`, `generation`, `sex`, `group`, `dob`, `start_date`, `bmt_date`,
    `death_date`, `early`, `reason`) and rows
    `("m1", "F1", early NA)` and `(NA, "F2", early "yes", reason "found dead")`.
    Expect:
    - two subjects;
    - `mouse_id` `c("m1", "F2")`;
    - `id_source` `c("mouse_ID", "facility")`;
    - `early` `c(FALSE, TRUE)`;
    - `early_reason` `c(NA, "found dead")`.
  - `db_write_meta errors on a row with neither mouse_ID nor id`: expect
    `regexp = "neither mouse_ID nor id"`.
  - `db_write_meta errors on unknown early values`: early `"maybe"` gives an
    error naming `"maybe"`.
    - Accepted values (case-insensitive, trimmed): `yes`, `ja`, `y`, `x`,
      `true`, `TRUE` -> TRUE.
    - `NA`, `""`, `no`, `nein`, `n`, `false` -> FALSE.
  - `re-running db_write_meta updates early but refuses end_date changes`:
    1. write once;
    2. change `early` / `reason` of m1, write again, expect the stored values
       updated;
    3. then change m1's `death_date`, expect an error matching
       `"end_date differs"` that names `m1` and both dates, and that nothing
       was written.
  - In `test-db_query.R`: `db_query_subjects returns dropout columns`. Expect
    the columns `mdc_id`, `early`, `early_reason`, `id_source`.
- [ ] **Step 2: Run** `devtools::test(filter = "db_write_dimensions|db_query")`.
  Expect FAIL.
- [ ] **Step 3: Implement.**
  - `db_write_meta()`:
    - filter `!is.na(mouse_ID) | !is.na(id)`;
    - set `mouse_id = dplyr::coalesce(mouse_ID, id)` and
      `id_source = ifelse(is.na(mouse_ID), "facility", "mouse_ID")`;
    - parse `early` with an internal `parse_early_(x)` that errors on unknown
      values;
    - map `early_reason = reason`;
    - `early` / `reason` stay optional columns: absent means FALSE / NA.
  - `db_write_subjects()`:
    1. **before** inserting, compare `end_date` with stored subjects; error
       (nothing written) if any non-NA stored `end_date` differs from a
       non-NA incoming one;
    2. insert with `ON CONFLICT (mouse_id) DO UPDATE SET early = excluded.early,
       early_reason = excluded.early_reason`. DuckDB's `excluded.` pseudo-table
       is SQL, not our concept. Leave other columns untouched, as today.
  - Default `early` is FALSE and `id_source` is `'mouse_ID'` when the columns
    are absent from `data`.
  - Update the roxygen of both functions and `db_query_subjects()`.
- [ ] **Step 4: Run** the same filter. Expect PASS. Also run
  `devtools::test(filter = "db_")` to confirm no regressions.
- [ ] **Step 5: Commit**: `feat(db): register dropouts by facility ID, store early/reason`.

### Task 3: `db_write_exclusions()` -- versioned rule sets

**Files:**
- Modify: `R/db_exclude.R` (add the exported function plus internals below the
  existing `match_sample_exclusions_()`)
- Test: `tests/testthat/test-db_exclude.R`

**Interfaces:**
- Consumes: Task 1 tables, Task 2 subjects.
- Produces:
  - `db_write_exclusions(con, exclusions, experiment_code, source_file = NA_character_, note = NA_character_)`
    returns the new `exclusion_set_id` (integer) invisibly, or the current one
    when unchanged.
  - Internal `current_exclusion_rules_(con, experiment_code)` returns a tibble
    `exclusion_set_id, mouse_id, tissue, panel, condition, file_name,
    collected_at, reason, decided_on`. It has 0 rows when there is no set.
  - Internal `normalise_exclusions_(exclusions)` maps sheet names to rule
    columns: `mouse_ID -> mouse_id`, `date -> decided_on`. `""` -> `NA`,
    character scope columns, `Date` for `collected_at` / `decided_on`. It
    sorts rows (by `mouse_id`, then the scope columns) so comparison is
    order-independent.

- [ ] **Step 1: Write the failing tests** (seed one experiment `"E1"`, subjects
  `m1`, `m2`, assay `overview` / facs):
  - `writes a first set and returns its id`: a sheet with one row
    `mouse_ID = "m1", reason = "r"`. Expect id `1L`, one rule, all scope
    columns `NA`.
  - `unchanged sheet writes no new set`: write twice. Expect the message
    `"unchanged"` and `count_rows_(con, "exclusion_sets") == 1`.
  - `blank strings equal NA` (Review Focus 2): rewrite the same sheet with
    `tissue = ""`. Expect "unchanged".
  - `empty sheet writes an empty set`: `exclusions = sheet[0, ]` and
    `exclusions = NULL` both write set 2 with 0 rules, and
    `current_exclusion_rules_()` returns 0 rows.
  - `validation errors`:
    - a missing `reason` column errors;
    - a row without `mouse_ID` or `reason` errors;
    - a `mouse_ID` not a subject of `"E1"` errors, naming it;
    - an unknown `experiment_code` errors.
    - Nothing is written in any case.
  - `unknown panel only messages`: `panel = "elisaX"` gives a message matching
    `"not registered yet"`, and the set is written.
  - `rules matching no stored row are reported` (Review Focus 5):
    - with no measurements stored, a message matching `"match no stored row"`
      that names `m1`;
    - a rule with non-NA `file_name` when every matching stain has
      `file_name` NULL also says `"re-export"`.
    - The matching itself uses the Task 5 views, so this test lands in Task 5
      Step 1. Here assert only the no-measurements case: any rule for a mouse
      with no samples is reported.
- [ ] **Step 2: Run** `devtools::test(filter = "db_exclude")`. Expect FAIL.
- [ ] **Step 3: Implement** in one `DBI::dbWithTransaction()`:
  1. validate;
  2. compare `normalise_exclusions_(exclusions)` with
     `current_exclusion_rules_()` (minus the id column, via `identical()` after
     the same sort);
  3. insert the set (`package_version = as.character(utils::packageVersion("JanisHelpers"))`);
  4. insert the rules.
  - Validation and messages happen before the transaction.
  - Roxygen: explain that a change of mind is editing the sheet then calling
    this function, and that history is kept in `exclusion_sets`.
- [ ] **Step 4: Run** the same filter. Expect PASS.
- [ ] **Step 5: Commit**: `feat(db): db_write_exclusions() versioned exclusion rule sets`.

### Task 4: Exports write every row

**Files:**
- Modify: `R/db_export.R` (`db_export_facs()`), `R/db_export_table.R`
  (`replace_import_()`, `prepare_sample_export_()`,
  `exclusion_date_hint_()`, delete `check_exclusion_panels_()`),
  `R/db_write_measurements.R` (`insert_facs_stain_()` gains `file_name`),
  `R/db_export_elisa.R`, `R/db_export_histo.R` (roxygen + messages)
- Test: `tests/testthat/test-db_export.R`, `test-db_export_elisa.R`,
  `test-db_export_histo.R`, `test-db_export_table.R`

**Interfaces:**
- Consumes: Task 1 `facs_stains.file_name`.
- Produces:
  - Exports keep their signatures. `exclusions` is now only used for messages
    about this unit's own rows.
  - Every prep row is written; `facs_stains.file_name` is filled.
  - `import_exclusions` is no longer written.
  - Export messages report `"<n> file(s)/sample(s) match an exclusion rule
    (stored, flagged on read)"`.

- [ ] **Step 1: Update and add tests**:
  - Existing tests that expect excluded rows to be absent, or
    `import_exclusions` rows to be present, now expect the rows **present** and
    `count_rows_(con, "import_exclusions") == 0`.
  - The "every file / sample is excluded" error tests now expect a successful
    export.
  - The "matches no file / sample" tests now expect a **message** (not an
    error) containing the unmatched row, and the export to succeed.
  - New: `db_export_facs stores file_name per stain`. After export, every
    `facs_stains.file_name` equals the prep's file names.
  - The date hint: a rule dated on a mouse's end date no longer appears as
    unmatched for an end-point ELISA sample, because the hint and the matching
    message compare against resolved dates, as the DB will.
    - Concretely: in `prepare_sample_export_()` resolve empty dates
      **before** calling `match_sample_exclusions_()`.
    - New text of `exclusion_date_hint_()`: `"End-point samples are stored
      with the mouse's end date -- a rule dated on that day matches them."`
- [ ] **Step 2: Run** `devtools::test(filter = "db_export")`. Expect FAIL.
- [ ] **Step 3: Implement.**
  - `db_export_facs()`:
    - remove the filter and the empty-data stop;
    - `stop` for unmatched becomes `message`;
    - pass `s$file_name` to `insert_facs_stain_()`;
    - delete the `import_exclusions` insert.
  - `prepare_sample_export_()`:
    - resolve dates first, then match;
    - unmatched becomes a message;
    - no anti-join;
    - return `list(data, matched)`.
  - `replace_import_()`: drop the `applicable` parameter and the insert.
  - Remove `check_exclusion_panels_()` and its calls (its job moves to
    `db_write_exclusions()`) and its tests in `test-db_export_table.R`.
  - Update roxygen ("drops the files hit by the exclusions sheet" ->
    "writes every file; exclusions are flagged on read, see
    db_write_exclusions()").
- [ ] **Step 4: Run** `devtools::test(filter = "db_export|db_write_measurements")`.
  Expect PASS.
- [ ] **Step 5: Commit**: `feat(db): exports write excluded data; store FACS file_name`.

### Task 5: Flag views and `excluded` columns on read

**Files:**
- Modify: `R/db_connect.R` (`create_schema_()` creates the views; it runs on
  every open, so `CREATE OR REPLACE` also migrates)
- Modify: `R/db_query.R` (all four `db_query_*()`)
- Test: `tests/testthat/test-db_query.R`, `tests/testthat/test-db_exclude.R`

**Interfaces:**
- Consumes: Tasks 1-4.
- Produces:
  - Views `facs_exclusion_flags(facs_measurement_id, excluded, exclusion_reason)`,
    `elisa_exclusion_flags(elisa_measurement_id, ...)`,
    `histo_exclusion_flags(histo_measurement_id, ...)`. Each lists only
    measurements with at least one matching rule; `excluded` is always TRUE
    there.
  - `db_query_facs/elisa/histo()` add `excluded` (logical,
    `COALESCE(..., FALSE)`) and `exclusion_reason`.
  - `db_query_subjects()` adds `excluded_mouse_wide` (logical, FALSE when no
    rule) and `exclusion_reason` (the reasons of mouse-wide rules).
  - Internal `rules_without_match_(con, experiment_code)` returns the current
    rules matching no stored row, with a `needs_reexport` flag for file rules
    when matching stains have `file_name` NULL. Used by
    `db_write_exclusions()`.

The match condition is shared by all three views. FACS is shown here; for
ELISA / histology `condition` / `file_name` must be NULL in the rule:

```sql
-- current set per experiment
WITH cur AS (
  SELECT experiment_id, MAX(exclusion_set_id) AS exclusion_set_id
  FROM exclusion_sets GROUP BY experiment_id)
SELECT m.facs_measurement_id, TRUE AS excluded,
       string_agg(DISTINCT r.reason, '; ' ORDER BY r.reason) AS exclusion_reason
FROM facs_measurements m
JOIN facs_stains fs ON fs.facs_stain_id = m.facs_stain_id
JOIN samples sa     ON sa.sample_id = fs.sample_id
JOIN subjects su    ON su.subject_id = sa.subject_id
JOIN assays a       ON a.assay_id = fs.assay_id
JOIN cur            ON cur.experiment_id = su.experiment_id
JOIN exclusion_rules r ON r.exclusion_set_id = cur.exclusion_set_id
 AND r.mouse_id = su.mouse_id
 AND (r.tissue IS NULL OR r.tissue = sa.tissue)
 AND (r.panel IS NULL OR r.panel = a.assay_name)
 AND (r.condition IS NULL OR r.condition = fs.condition)
 AND (r.file_name IS NULL OR r.file_name = fs.file_name)
 AND (r.collected_at IS NULL OR r.collected_at = sa.collected_at)
GROUP BY m.facs_measurement_id
```

FACS rules with a non-NULL `collected_at` and a named FACS panel are rejected
by `facs_match_exclusions()` today. Keep that check in the FACS export message
path. The view needs no special case.

- [ ] **Step 1: Write the failing tests**:
  - **Equivalence** (`test-db_query.R`, `exclusion flags match the R matchers`):
    1. export the FACS fixture prep used by `test-db_export.R` and an ELISA
       prep with serial + end-point samples;
    2. write a sheet with:
       - a mouse-wide rule;
       - a tissue rule;
       - a FACS `condition` rule;
       - a FACS `file_name` rule;
       - an ELISA serial-date rule;
       - an ELISA rule dated on a mouse's end date;
    3. expect the set of flagged `file_name`s from
       `db_query_facs() |> dplyr::filter(excluded)` to equal
       `facs_match_exclusions(prep$data, sheet, assay)$excluded_files$file_name`;
    4. expect the flagged ELISA `mouse_id x tissue x collected_at` to equal
       `match_sample_exclusions_()` on the date-resolved data. The end-date
       rule is the documented intended difference: it **is** flagged.
  - `multiple rules give one row per measurement` (Review Focus 1): the
    mouse-wide and tissue rule on the same mouse. Expect
    `nrow(db_query_facs() |> collect())` unchanged vs no rules and
    `exclusion_reason == "r_mouse; r_tissue"` (sorted).
  - `no set means nothing excluded` (Review Focus 3): before any
    `db_write_exclusions()`, `excluded` is all FALSE (not NA). The same after
    writing an empty set.
  - `rules do not cross experiments` (Review Focus 4): two experiments with the
    same tissue and panel; a rule in E1 flags no E2 row.
  - `old stains without file_name are not flagged by file rules` (Review Focus
    5): set `file_name = NULL` by SQL. The file rule flags nothing, and
    `db_write_exclusions()` messages `"re-export"` naming the rule.
  - `changing the sheet changes the flags without re-export`: write set 1 with
    m1 excluded, then set 2 without. `excluded` all FALSE; `exclusion_sets` has
    2 rows.
  - `db_query_subjects flags mouse-wide exclusions`: `excluded_mouse_wide` TRUE
    only for the mouse-wide rule's mouse.
- [ ] **Step 2: Run** `devtools::test(filter = "db_query|db_exclude")`. Expect
  FAIL.
- [ ] **Step 3: Implement.**
  - Create the three views in `create_schema_()` after the tables.
  - In `db_query_*()`, `dplyr::left_join(dplyr::tbl(con, "<domain>_exclusion_flags"), by = "<id>")`
    then `dplyr::mutate(excluded = dplyr::coalesce(excluded, FALSE))`.
    - The FACS query must now carry `facs_measurement_id` through the join
      before the final `select()`.
    - Add `excluded` and `exclusion_reason` at the end of each `select()`.
  - `db_query_subjects()`: left-join an aggregate of the current set's
    mouse-wide rules (all scope columns NULL) per `mouse_id` and experiment.
  - Wire `rules_without_match_()` into `db_write_exclusions()` (after the
    write), replacing Task 3's simple no-samples check.
- [ ] **Step 4: Run** `devtools::test(filter = "db_")`. Expect PASS.
- [ ] **Step 5: Commit**: `feat(db): flag excluded rows on read via exclusion views`.

### Task 6: Prep runners write the rule set

**Files:**
- Modify: `R/facs_prep.R`, `R/elisa_prep.R`, `R/histo_prep.R` (`*_prep_run()`)
- Modify: `inst/scripts/facs_prep_template.R`, `elisa_prep_template.R`,
  `histo_prep_template.R`
- Test: `tests/testthat/test-facs_prep.R`, `test-elisa_prep.R`, `test-histo_prep.R`

**Interfaces:**
- Consumes: `db_write_exclusions()` (Task 3).
- Produces:
  - New last-but-one argument `meta_path = NULL` on all three `*_prep_run()`
    (before `data_dir`).
  - With `final = TRUE`, before the export, they call
    `db_write_exclusions(con, exclusions, experiment_code_of_prep_(con, prep), source_file = <meta_path relative to data_dir, or NA>)`.
  - Internal `experiment_code_of_prep_(con, prep)` (in `R/db_exclude.R`)
    returns the single experiment code of the prep's `mouse_ID`s. It errors
    when they span several experiments or none is registered.

- [ ] **Step 1: Write the failing tests**:
  - `final run writes the exclusion set` (one per domain, reusing each file's
    existing final-run test setup): after `*_prep_run(final = TRUE, ...,
    exclusions = sheet)`, expect `count_rows_(con, "exclusion_sets") == 1` and
    the rule stored.
  - `meta_path fills source_file`: expect the relative path stored.
  - `prep spanning two experiments errors` (FACS only).
- [ ] **Step 2: Run** `devtools::test(filter = "prep")`. Expect FAIL.
- [ ] **Step 3: Implement.** Templates: add `meta_path = meta_path,` (the
  variable each template already uses for `meta_read()`; introduce it if the
  template inlines the path), with a comment that `final = TRUE` also records
  the exclusions sheet. Update roxygen.
- [ ] **Step 4: Run** `devtools::test(filter = "prep")`. Expect PASS.
- [ ] **Step 5: Commit**: `feat(prep): final runs record the exclusions sheet as a rule set`.

### Task 7: `db_load_experiment()` -- filtering, `keep_excluded`, `cohort`

**Files:**
- Modify: `R/db_load.R`
- Test: `tests/testthat/test-db_load.R`

**Interfaces:**
- Consumes: Task 5 columns.
- Produces:
  - Return list order: `meta`, `cohort`, `exclusions`, then the units.
  - `cohort` = every subject of the experiment (`db_query_subjects()` rows,
    group levels applied).
  - `meta` = `cohort` minus `id_source == "facility"` minus
    `excluded_mouse_wide` (unless `keep_excluded`, which keeps all subjects).
    `meta` drops the `id_source` / `excluded_mouse_wide` columns only when
    `keep_excluded = FALSE`.
  - `exclusions` = `current_exclusion_rules_()`.
  - Units: `keep_excluded = FALSE` drops rows with `excluded` and drops the
    `excluded` / `exclusion_reason` columns; TRUE keeps both.
  - `"cohort"` added to the reserved names in `resolve_unit_names_()`.
  - `load_exclusions_()` and `fully_excluded_mice_()` are deleted; Task 8
    removes the other caller.

- [ ] **Step 1: Update and add tests**:
  - Existing default-load tests must pass unchanged in data content: units have
    the same rows and columns. `meta` has the same rows; the new subject
    columns from Task 2 (`mdc_id`, `early`, `early_reason`) are allowed extras.
    Update only the list names and the `exclusions` columns (now `exclusion_set_id, mouse_id, tissue, panel,
    condition, file_name, collected_at, reason, decided_on`).
  - New: `default load drops flagged rows and facility-ID dropouts`.
  - New: `keep_excluded keeps flagged rows with flags`.
  - New: `cohort holds every subject including dropouts`.
  - New: `no rule set gives a message`, matching
    `"No exclusions recorded for experiment"`.
  - New: `"cohort" is a reserved name`.
- [ ] **Step 2: Run** `devtools::test(filter = "db_load")`. Expect FAIL.
- [ ] **Step 3: Implement** and update the roxygen (`@returns`, the
  `keep_excluded` param).
- [ ] **Step 4: Run** `devtools::test(filter = "db_load")`. Expect PASS.
- [ ] **Step 5: Commit**: `feat(db): db_load_experiment() filters flagged rows, adds cohort`.

### Task 8: Cohort tables and analysis template use `dat$cohort`

**Files:**
- Modify: `R/meta_tables.R` (`cohort_input_()`, `analysis_status_()`,
  `data_units_()`, roxygen of `meta_summary_table()` / `meta_mouse_table()`)
- Modify: `inst/rmarkdown/templates/analysis-db/skeleton/skeleton.Rmd`
- Test: `tests/testthat/test-meta_tables.R`

**Interfaces:**
- Consumes: Task 7 `dat$cohort`, `dat$exclusions` (new columns, no
  `import_unit`).
- Produces:
  - `workbook = NULL` uses `dat$cohort` silently; it falls back to `dat$meta`
    with today's message only when `dat$cohort` is absent (older lists).
  - `cohort_normalise_()` maps `early_reason` to `reason` (pick
    `"reason", "early_reason"`).
  - Excluded = `excluded_mouse_wide` when present in the cohort, else
    mouse-wide rules in `dat$exclusions` (all scope columns NA).
  - `"cohort"` is added to `data_units_()`'s skip list.
  - `not in DB` stays for workbook input.

- [ ] **Step 1: Update and add tests**:
  - Existing workbook-based tests keep passing. Adapt fixtures building
    `dat$exclusions` with `import_unit` to the new columns.
  - New: `meta_summary_table uses dat$cohort when workbook is NULL`: a `dat`
    with a facility-ID dropout in `cohort`. Expect it counted in "Dropped
    before / at BMT" and no message.
  - New: `meta_mouse_table shows DB dropouts`: the dropout row has status
    `early: <reason> ...`, not `not in DB`.
- [ ] **Step 2: Run** `devtools::test(filter = "meta_tables")`. Expect FAIL.
- [ ] **Step 3: Implement.** Template: replace the `dat$workbook` calls with
  `meta_summary_table(NULL, dat)` / `meta_mouse_table(dat)`. Keep the
  commented `dat$workbook <- meta_read(...)` line as an option.
- [ ] **Step 4: Run** `devtools::test(filter = "meta_tables|report")`. Expect
  PASS (the template parse test lives under report).
- [ ] **Step 5: Commit**: `feat(meta): cohort tables read dropouts and exclusions from the DB`.

### Task 9: Docs, NAMESPACE, full check

**Files:**
- Modify: `CLAUDE.md` (data-structures bullets for `db_load_experiment()`,
  `db_connect()/db_write_*()/db_query_*()`, `meta_summary_table()`, ELISA /
  FACS export lines that say "drops excluded rows")
- Modify: `docs/superpowers/specs/2026-07-17-facs-db-schema-design.md` (one
  line pointing to the new spec)
- Generated: `NAMESPACE`, `man/*.Rd`

- [ ] **Step 1:** In `CLAUDE.md`, replace every claim that exports drop
  excluded data or that the DB lacks dropouts with the new behaviour. One line
  each, pointing to `docs/superpowers/specs/2026-11-20-exclusions-in-db-design.md`.
  Add `db_write_exclusions()`.
- [ ] **Step 2:** Run `Rscript --vanilla -e 'devtools::document()'`. Expect
  `db_write_exclusions` exported in `NAMESPACE` and no roxygen warnings.
- [ ] **Step 3:** Run `Rscript --vanilla -e 'devtools::test()'`. Expect 0
  failures.
- [ ] **Step 4:** Run `Rscript --vanilla -e 'devtools::check()'`. Expect 0
  errors, 0 warnings, and only the NOTEs listed in `CLAUDE.md` "Known check
  output".
- [ ] **Step 5: Commit**: `docs: exclusions-in-db in CLAUDE.md and schema spec`.

## After merge (user, not part of the build)

For each experiment, in this order:

1. re-run `db_write_meta()`, which registers dropouts and early flags;
2. re-run each panel / assay script with `final = TRUE`, which writes the full
   data and the first rule set.

Until then, old units lack the excluded rows, and file-level rules can't match
old FACS stains.
