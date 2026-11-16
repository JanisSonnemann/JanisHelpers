# FACS pipeline -- overview & progress

Full design: `docs/superpowers/specs/2026-11-12-facs-pipeline-workflow-design.md`
(restructured 2026-11-13: stages reordered, descriptions cut to goals;
status updated 2026-11-16)

**Goal:** one database with only final, QC-approved FACS data for every
experiment. Analysis starts with a query, not with re-reading workspaces.

```
FlowJo -> dated .wsp -> prep script, final <- FALSE
            (read, calc, save .rds, QC report; repeat until clean)
       -> final <- TRUE: reload .rds, apply exclusions, export -> DB
            -> snapshot -> analysis (db_query_*)
```

## Key rules

- Workspace name = `<root>_<YYYYMMDD>.wsp` (root = experiment x panel,
  date = version). DB holds one version per root; `import_log` keeps history.
- One self-contained prep script per panel, showing **all** samples;
  exclusions live in the meta workbook's `exclusions` sheet and are applied
  only at export.
- Volumes / organ weights are calculation inputs, not stored. DB gets raw
  counts, frequencies and calculated values.
- Wells differing in one variable (e.g. restimulation) are separate stains
  via `condition` + `condition_keyword`; ratios between wells are computed in
  analysis, not stored.

## The pipeline as built (2026-11-16)

One prep script per panel (from `inst/scripts/facs_prep_template.R`),
one switch:

**`final <- FALSE`** -- prep / QC loop, repeat until clean
1. Read the workspace and prepare the data (population selection,
   calculations, metadata).
2. Write the prep `.rds`.
3. Create the QC report.

**`final <- TRUE`** -- export
1. Register experiment + mice in the DB.
2. Reload the `.rds` and apply exclusions.
3. Write to the DB (replaces the stored version, logs the import).
4. Snapshot the DB.

## Stages

### 1. Prep / QC loop + final export -- done

**Goal:** turn a FlowJo workspace into QC-approved data in the DB with one
switch (`final <- TRUE`).

- [x] Prep script template, `.rds` save/load, exclusion matching,
      `db_export_facs()` (overwrite prompt, one transaction, logging).
      Merged at `e15ed9b`. Design / plan: `2026-11-12-facs-prep-qc-export*`.
- [x] QC checks in the report (merged to master 2026-11-13, tuned on 26-1)
      -- completeness, gating hierarchy, raw counts, MFI, batch plots,
      ranked flag table. Design / plan: `2026-11-13-facs-qc-checks*`.
- [x] `facs_prep_run()` (one call per script) and `db_snapshot()`
      (`849e8a2`); `db_write_meta()` registers experiment + mice from the
      standard `meta` sheet (`02e45b0`).
- [x] 2026-11-16: outlier labels switched from robust z-score to the
      boxplot-whisker rule, `min_events` became a stored setting (`88db396`);
      report split into "Fraction of parent", "Raw counts" and "Batch
      effects" sections (`aeb18c9`).

### 2. First data in the DB -- next

**Goal:** a real experiment (26-1) in the DB through the full pipeline, plus
the old data brought over. (Tracked in detail in `2026-11-12-db-setup.md`
step 6.)

- [x] `facs_prep_26-1_{overview,treg,mpo-restim}.R` written in the MPO_mRNA
      project; dry run of all three panels passed (294 stains, 22,675 rows,
      re-run inserts 0). Treg `mouse_ID` keyword fixed.
- [ ] **(user)** Fill in the Treg `.wsp` name, Treg / MPO-Restim population
      selections and QC settings (`group_col`, `batch_col`,
      `expected_mice`, `qc_mfi_populations`, ...); rerun with
      `final <- FALSE` and review the QC reports.
- [ ] **(user)** Fill the `exclusions` sheet, then run each script with
      `final <- TRUE`; check with `db_query_facs()`.
- [ ] ELISA for 26-1 (hand-made wide sheet, not the `elisa_read_results()`
      format; needs plate sample -> mouse map).
- [ ] Migrate old `MPOmRNA.db` (SQLite) into the new schema
      (`2026-11-12-db-setup.md` step 7).

### 3. Workspace versioning

**Goal:** every gate or compensation change is traceable -- which `.wsp`
version is in the DB and what changed between versions.

- [x] Naming convention `<root>_<YYYYMMDD>.wsp`; workspaces live in the
      experiment's `raw/FACS/`, one folder per panel for the `.fcs`.
- [ ] Per-experiment changelog `.md`.
- [ ] Decide how much to automate: changelog template only, or a `facs_`
      function that diffs two `.wsp` versions and drafts the entry.

### 4. Analysis

**Goal:** analysis scripts start from a DB query and feed `analysis_` /
`plot_` directly, also across experiments.

- [ ] `db_query_facs()` with experiment / population / metric selection and
      metadata joined (currently returns the whole lazy table).
- [ ] Population-name harmonization across experiments (prep script vs.
      lookup table in the DB) -- needed for cross-experiment comparisons.

Each stage gets its own design + plan when started.

## Open questions

- Calculations per panel (reference populations, ratios) -- defined by you in
  each prep script.
- Versioning (stage 3) vs. `db_query_facs()` (stage 4) first, once 26-1 is
  in -- the query function pays off as soon as there is data to analyse.
