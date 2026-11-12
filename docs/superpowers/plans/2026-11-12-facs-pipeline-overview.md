# FACS pipeline -- overview & progress

Full design: `docs/superpowers/specs/2026-11-12-facs-pipeline-workflow-design.md`
Branch / worktree: `worktree-facs-pipeline-workflow`

**Goal:** one database with only final, QC-approved FACS data for every
experiment. Analysis starts with a query, not with re-reading workspaces.

```
FlowJo -> dated .wsp + changelog -> prep script (calc, QC report; repeat)
       -> final <- TRUE: apply exclusions, export (log, overwrite) -> DB -> analysis
```

## Key rules

- Workspace name = `<root>_<YYYYMMDD>.wsp`, e.g. `treg_26-2_20261112.wsp`.
  The root (`treg_26-2`) = panel/experiment; the date = version.
- Major change (compensation, gates across samples) -> save a new dated
  `.wsp` + write the changelog entry.
- One prep script per panel: read `.wsp`, select populations, calculate,
  save `.rds`, render the QC report. Rerun after each gate change.
- Prep and report show **all** samples. Exclusions go in the meta workbook's
  `exclusions` sheet (whole mouse down to one `.fcs`, with a reason); the
  report only flags them.
- `final <- TRUE` = QC passed: exports the saved `.rds`, applies exclusions and
  stores them in the DB. Root not in DB -> added; root already in DB ->
  **asks you whether to overwrite**. Writes an import log row, then takes a
  DB snapshot.
- Stain volumes / organ weights stay **out** of the DB: the prep script uses
  them to calculate values. DB gets raw counts, frequencies **and** calculated
  values (per g, % of CD45+, ratios, ...).
- DB holds one version per root, never more. The log keeps the history.
- Workspaces carry `mouse_ID`, `tissue`, `panel` keywords -> used for mapping.

## Stages

- [ ] **0. Finish DB setup** -- steps 2-7 of `2026-11-12-db-setup.md`
      (trim CLAUDE.md, `.wsp` provenance fix, merge `db_` branch, first import).
      Before the first import: schema changes from Stage 3 (drop volume
      columns + `sample_digestions`, add `wsp_root`, add `import_log`).
- [ ] **1. Workspace versioning** -- filename convention, per-experiment
      changelog `.md`, `facs_` function that diffs two `.wsp` versions
      (compensation + gates) and drafts the changelog entry
- [ ] **2+3. Prep / QC loop + final export** -- design:
      `2026-11-12-facs-prep-qc-export-design.md`, plan:
      `2026-11-12-facs-prep-qc-export.md`. Exclusion matching, prep
      `.rds` save/load, `facs-qc-template` + `report_` wrapper, `db_` export
      (overwrite prompt, one transaction, `import_log`, new
      `import_exclusions` table), prep script template; rework the 26-1
      scripts into `facs_prep_26-1_<panel>.R`. QC checks added later.
- [ ] **4. Analysis queries** -- `db_query_facs()` with experiment /
      population / metric selection, metadata joined, ready for `analysis_`/`plot_`

Order: 0 -> 1 -> 2+3 -> 4 (1 doesn't need the DB and can run earlier).
Each stage gets its own design + plan when started.

## Open issues

- [ ] **Population-name harmonization** (parked): same population, different
      gate names across experiments -- fix in prep script or via a lookup
      table in the DB? Needed for cross-experiment comparisons in Stage 4.
- [ ] **Calculations per panel** -- which reference populations / ratios;
      you define them in each panel's prep script.
- [ ] **QC checks** -- user to supply functions/analyses for the template.
- [ ] **Version comparison** -- manual, on demand only; format undecided.
- [ ] **Old `MPOmRNA.db` migration** (db-setup step 7).
- [ ] Where `.wsp` versions live (suggested: `raw/<experiment_code>/facs/`).

## Calls made without you (check when convenient)

- DB setup plan absorbed as Stage 0.
- Prep result and QC report go into new `processed/` and `qc/` folders in the
  panel's raw `.fcs` folder, e.g. `raw/FACS/overview/{processed,qc}/` (your
  call; one per `.wsp` version, overwritten on rerun).
- `final <- TRUE` exports the saved `.rds` without re-reading; refuses if it's
  older than the `.wsp`. No separate QC prompt.
- Exclusion rows need a reason; an applicable row matching no file stops the
  export (only a warning in the report).
- Export requires subjects/samples already registered (errors, doesn't create).
- Prep script = a copyable per-panel R template (run interactively) + package
  helpers; the script decides which calculations run.
- Calculated values stored as extra `metric` rows with one naming convention.
- Volume/weight fixes go through the same re-export.
- Schema changes done before the first real import (no migration).
- No hard QC-report <-> `.wsp` match check; only a newer-version warning.
