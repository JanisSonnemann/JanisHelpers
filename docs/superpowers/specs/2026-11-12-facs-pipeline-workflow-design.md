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
- Every workspace version that goes on to the database has a dated QC report.
- Exporting to the database is one deliberate, logged step, gated by the
  user's explicit "QC passed".
- Analysis across experiments = connect, select experiments + populations,
  get one tidy tibble.

## Non-goals

- Automating the QC decision. The user decides whether QC passed.
- Keeping multiple versions of the same experiment's data in the database.
- Automated version comparison as part of the pipeline (it's manual, on
  demand -- see Stage 2).
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
QC report per .wsp version (facs-qc-template.rmd)               [Stage 2]
   |  (optional, manual) version comparison old vs new
   |  user decides: QC passed
   v
Export script: .wsp + meta -> raw + calculated metrics -> DB    [Stage 3]
   root exists? no -> add; yes -> ask user to overwrite
   writes import log row (qc_passed, report, .wsp version, time)
   then dated DB snapshot
   v
Analysis: db_query_facs(experiments, populations) -> analysis_/plot_  [Stage 4]
```

### Stage 0 -- Finish the database setup (prerequisite)

The existing `docs/superpowers/plans/2026-11-12-db-setup.md` (steps 2-7) is
absorbed as the first stage of this plan, not run in parallel *(my call)*:
trim CLAUDE.md, fix `db_write_facs()` provenance so `source_file` records the
`.wsp` (step 3 -- directly needed by Stage 3), fix `db_write_elisa()`, merge
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

### Stage 2 -- QC report per workspace version

Decided:
- One dated QC report per workspace version, from a new
  `facs-qc-template.rmd`.
- The actual QC checks/functions are supplied later by the user. This stage
  delivers the **skeleton**: template + parameters + rendering into the
  right place.
- After a major change, QC is redone from scratch on the new version.
- A version comparison (old vs new populations) is **optional and manual**,
  produced only when the user wants it, as a separate file. Not part of the
  pipeline.

Left to the builder:
- Template under `inst/rmarkdown/templates/facs-qc-template/`, parameterized
  by `.wsp` path (Rmd `params`), rendered via the existing
  `report_knit_dated()`/`report_knit_exp()` helpers or a thin `report_`
  wrapper.
- Output location -- suggest `raw/<experiment_code>/facs/qc/<YYYY-MM-DD>_qc.html`
  next to the workspace so Stage 3 can reference it *(my call)*.
- Placeholder sections for the future checks.

### Stage 3 -- Export to database (logged, overwriting)

Decided:
- The user runs the export manually once satisfied with QC.
- The user records "QC passed" manually at export time.
- **Unit = workspace root.** Exporting `treg_26-2_20261120.wsp` concerns
  exactly the data of root `treg_26-2` -- nothing else. Other roots of the
  same experiment are untouched. The database never holds two versions of
  one root.
- **Existing-entry check.** The export first checks whether data for the
  root already exists in the database. If not, it adds it. If it does, it
  **asks the user whether to overwrite**; only on yes is the old data
  replaced.
- **Calculations happen in the export, results are stored.** Stain volumes
  and organ weights are **not** stored in the database. They are read in the
  export script only, to calculate absolute and derived values. The database
  receives all raw counts, frequencies (FlowJo statistics) **and** the
  calculated values (counts per g, % of CD45+ / other reference populations,
  ratios, ...).
- FCS file -> mouse/tissue/panel mapping comes from the workspace keywords
  `mouse_ID`, `tissue`, `panel` (panel = assay).
- Logging with timestamp etc. is required.

Consequence for the schema (on the unmerged `db_` branch): the volume /
count-method / bead columns on `facs_stains`, its cross-field `CHECK`, and
the `sample_digestions` table are no longer needed for FACS, and
`db_write_facs()` loses those arguments. Calculated values are stored as
additional `metric` rows in `facs_measurements`.

Proposed *(my calls)*:
- Export = a per-panel **plain R script** (decided: run by the user
  interactively in R, never knitted), copied from a template and adapted. It
  reads the `.wsp` with
  `facs_read_wsp()`, reads the meta workbook with `meta_read()`/`meta_clean()`,
  calculates derived values with the `facs_calc_*()` functions, then hands one
  long tibble (raw + calculated metrics) to a package export function (e.g.
  `db_export_facs()` -- name to settle in its design). Which calculations
  are done is defined per panel in the script, not in the package.
- The package export function takes `qc_passed = TRUE` (refuses otherwise),
  the QC report path and the `.wsp` path. If the root already exists it
  prompts (e.g. `utils::askYesNo()`); "no" aborts without changes.
- Naming of calculated metrics (e.g. `"CountPerG"`, `"PctOf_CD45+"`,
  `"Ratio_CD4_CD8"`) follows one documented convention, settled in the Stage 3
  design.
- Experiment, subjects and samples must already be registered; the export
  errors listing every unresolved `mouse_ID` x `tissue` / panel (same pattern
  as `db_write_samples()`), it does not create them.
- The root has to be stored with the data so it can be found/deleted by root
  (e.g. a `wsp_root` column on `facs_stains`) -- schema addition.
- Overwrite happens in **one transaction**: delete the root's existing
  `facs_measurements` and `facs_stains` rows, insert the new ones, write the
  log row -- all or nothing.
- A correction of a volume or weight goes through the same re-export of the
  current `.wsp` (recalculates, prompts to overwrite) -- no separate path.
- All schema changes (drop volume columns / `sample_digestions`, add
  `wsp_root`, add `import_log`) are made on the `db_` branch **before the
  first real import** (db-setup step 6), so no data migration is needed.
- No hard check that the QC report matches the `.wsp` (avoid ceremony), but
  the export **warns** if a newer dated version of the same root exists on
  disk, and the log records both paths.
- **Import log table** (new, e.g. `import_log`), one row per export, kept
  even after the data it describes is overwritten: workspace root,
  experiment, domain, `.wsp` file + version date, QC report path,
  `qc_passed`, exported_at, JanisHelpers package version, rows written, and
  the id of the log row it superseded (previous export of the same root).
  Requires a schema addition to the DuckDB design.
- After a successful export, take the dated `EXPORT DATABASE` snapshot already
  described in the db-setup plan.

Left to the builder: exact log columns, whether the log also covers ELISA/
histo imports (suggest yes, `domain` column), how FCS files map to stains.

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

Stage 0 -> 1 -> 2 -> 3 -> 4. Stages 1 and 2 are independent of the database
and could run before or alongside Stage 0; Stage 3 needs Stage 0 merged and
the Stage 2 report location fixed. Each stage gets its own design + plan
when picked up.
