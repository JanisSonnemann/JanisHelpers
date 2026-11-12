# FACS pipeline workflow: workspace -> QC -> database -> analysis

Status: draft, pending review
Type: big-picture working plan. Each stage below becomes its own design/plan
when it is picked up; this document fixes the order, the hand-offs between
stages, and the decisions already made.

## Why

Three recurring problems:

1. **Untraceable workspace changes.** Compensation and gates get changed in
   FlowJo, and later nobody can reconstruct why populations look different.
2. **No standardized QC.** Data moves from FlowJo to analysis without a fixed
   check.
3. **Tedious retrieval.** Every cross-experiment comparison means running
   `facs_read_wsp()` per workspace, then filtering populations and reapplying
   calculations/corrections by hand.

Target: one **single source of truth** -- a database holding only finalized,
QC-approved data for every experiment. Analysis starts with a query, not with
re-reading workspaces.

## Goals

- Every major workspace change is a dated, saved `.wsp` version with a
  changelog entry explaining what changed and why.
- Every workspace version that goes on to the database has a QC report,
  rendered from exactly the data that gets exported.
- Exporting to the database is one deliberate, logged step (the `final`
  switch in the panel's prep script), which is the user's "QC passed".
- Excluded samples are dropped only at export, and the database records what
  was excluded and why.
- Analysis across experiments = connect, select experiments + populations,
  get one tidy tibble.

## Non-goals

- Automating the QC decision. The user decides whether QC passed.
- Keeping multiple versions of the same experiment's data in the database.
- Automated version comparison as part of the pipeline (it's manual, on
  demand -- see Stages 2 + 3).
- Analysis-side audit trail (logging calls to `facs_cluster_flowsom()` etc.,
  idea on branch `worktree-facs-roadmap-audit-trail`). Separate effort.

## Anti-goals (would make this a failure even if it works)

- The database contains stale, superseded, or non-QC'd data -- "only final and
  clean data, to avoid confusion".
- Overwriting data loses the record of *what* was replaced and *when*.
- The workflow adds so much ceremony that it's skipped in practice.

## The workflow

```
FlowJo (.fcs analysed)
   |  major change? -> save dated .wsp + changelog entry        [Stage 1]
   v
Prep script, final <- FALSE  (repeat)                            [Stages 2+3]
   .wsp -> select populations -> facs_calc_*()   all samples
   -> save prep result (.rds) -> render QC report (flags exclusions)
   review -> adjust gates (new .wsp version) / edit exclusions sheet -> rerun
   v
Prep script, final <- TRUE
   load saved .rds -> apply exclusions sheet -> DB
   root exists? no -> add; yes -> ask user to overwrite
   import log row + exclusion rows, then dated DB snapshot
   v
Analysis: db_query_facs(experiments, populations) -> analysis_/plot_  [Stage 4]
```

### Stage 0 -- Finish the database setup (prerequisite)

The existing `docs/superpowers/plans/2026-11-12-db-setup.md` (steps 2-7) is
absorbed as the first stage of this plan, not run in parallel *(my call)*:
trim CLAUDE.md, fix `db_write_facs()` provenance so `source_file` records the
`.wsp` (step 3 -- directly needed by Stages 2 + 3), fix `db_write_elisa()`, merge
branch `worktree-db-domain-functions`, create the OneDrive folders, first
import. Decisions recorded there (DuckDB at `~/lab-db/lab.duckdb`, raw files
under OneDrive `Forschung/LabDB/raw/<experiment_code>/facs/`, snapshots under
`snapshots/YYYY-MM-DD/`) stand.

### Stage 1 -- Workspace versioning + changelog

Decided:
- Work continues in FlowJo as today. A **major change** = compensation
  changed for one or more channels, or one or more gates adjusted across
  multiple samples. On a major change the user saves a new dated copy of the
  `.wsp`; all versions are kept.
- The changelog is a per-experiment `.md` file, written by the user. It stays
  the human record of *why*.

Proposed *(my suggestion, user said "sounds good")*: a function that
compares two `.wsp` versions directly from the XML and drafts the changelog
entry -- which spillover values changed (channel pair, samples affected) and
which gates changed (gate path, samples affected). The user adds the reason.

- **Filename convention (decided):** `<root>_<YYYYMMDD>.wsp`, e.g.
  `treg_26-2_20261112.wsp`. The *root* (`treg_26-2`) names the panel or
  experiment and is the unit of identity throughout the pipeline (QC, export,
  overwrite, log). The date suffix is the version.
- Every workspace carries `mouse_ID`, `tissue` and `panel` as FlowJo keywords
  (confirmed by user).

Left to the builder:
- Where versions live (suggest `raw/<experiment_code>/facs/` beside the
  `.fcs` folder) *(my call, to confirm in the Stage 1 design)*.
- Changelog filename/format (suggest `CHANGELOG.md` per experiment, one
  dated heading per version, an auto-drafted "what changed" block plus a
  free-text "why").
- Function name within the `facs_` domain (e.g. `facs_diff_wsp()`), how gate
  geometry differences are summarized, and tolerance for float noise.

### Stages 2 + 3 -- Prep / QC loop and final export

Designed in detail in
`docs/superpowers/specs/2026-11-12-facs-prep-qc-export-design.md`, which
replaced the earlier separate "QC report" and "export script" stages. In
short:

- One per-panel **prep script**, run interactively, handles the whole loop.
  It reads the `.wsp`, selects populations and calculates derived values.
  It saves the result as `.rds` and renders a QC report from it. You rerun
  it after every gate change.
- Prep and report show **all** samples. Exclusions live in an `exclusions`
  sheet of the meta workbook (from a whole mouse down to one `.fcs` file,
  each with a reason). The report only flags them.
- `final <- TRUE` in the same script exports the saved `.rds` (no re-read),
  applies the exclusions and stores them in a new `import_exclusions` table
  linked to the import's log row.
- The QC checks themselves are supplied later by the user. A version
  comparison (old vs new populations) stays optional and manual.

Unchanged from the original Stage 3 decisions (built in
`2026-11-12-db-schema-import-log-design.md` or carried into the design above):

- **Unit = workspace root.** The database never holds two versions of one
  root. If the root already exists, the export asks before overwriting
  (`utils::askYesNo()`); "no" aborts with no changes. The overwrite is one
  transaction: delete the old stains and measurements, insert the new ones
  and write the log row.
- **Calculations happen before export, results are stored.** Volumes and
  organ weights are calculation inputs and never stored. The database gets raw
  counts, frequencies and calculated values as extra `metric` rows. Which
  calculations run is defined per panel in the script.
- FCS file -> mouse/tissue/panel mapping comes from the workspace keywords
  `mouse_ID`, `tissue`, `panel`. Experiment and subjects must already be
  registered.
- **Import log** (`import_log`): one row per export, kept after overwrite,
  with root, `.wsp` file + version date, QC report path, `qc_passed`,
  timestamp, package version, rows written and the superseded log id.
- The export warns if a newer dated `.wsp` of the same root exists on disk.
  There's no hard check that the report matches the `.wsp` beyond the prep
  staleness check.
- A volume or weight correction goes through the same loop: rerun prep, then
  export with `final <- TRUE`.
- A dated `EXPORT DATABASE` snapshot follows each export.

### Stage 4 -- Analysis from the database

Decided:
- Analysis starts from the database: select experiments and populations,
  run the analysis. No more per-workspace reads.

Because calculated values are stored at export, Stage 4 needs no
calculation step: queries return raw and calculated metrics ready for the
`analysis_`/`plot_` functions. Cross-experiment comparability of population
names still depends on the parked harmonization item.

Left to the builder: extend `db_query_facs()` (on the `db_` branch) with
experiment/population/metric selection, subject metadata (group, tissue)
joined in, and one consistent column naming (CamelCase from
`facs_import_wsp()` vs snake_case elsewhere).

## Parked (decide later)

- **Cleaning before export.** Decided: derived values are calculated in the
  export script and stored. Still open: population-name harmonization across
  experiments (same population, different gate names) and any other
  corrections, and where they happen (in the export script vs. a lookup
  table in the database). Needed for Stage 4 cross-experiment queries.
- **Which calculations per panel** (which reference populations, which
  ratios) -- user defines them in each panel's export script.
- The QC checks themselves (user supplies).
- Migrating the old `MPOmRNA.db` SQLite data (db-setup step 7).

## Order

Stage 0 -> 1 -> 2+3 -> 4. Stage 1 is independent of the database and could
run before or alongside Stage 0; Stages 2+3 need Stage 0 merged. Each stage
gets its own design + planwhen picked up.
