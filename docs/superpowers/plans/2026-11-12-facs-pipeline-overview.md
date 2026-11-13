# FACS pipeline -- overview & progress

Full design: `docs/superpowers/specs/2026-11-12-facs-pipeline-workflow-design.md`
(restructured 2026-11-13: stages reordered, descriptions cut to goals)

**Goal:** one database with only final, QC-approved FACS data for every
experiment. Analysis starts with a query, not with re-reading workspaces.

```
FlowJo -> dated .wsp + changelog -> prep script (calc, QC report; repeat)
       -> final <- TRUE: apply exclusions, export -> DB -> analysis
```

## Key rules

- Workspace name = `<root>_<YYYYMMDD>.wsp` (root = experiment x panel,
  date = version). DB holds one version per root; `import_log` keeps history.
- One prep script per panel shows **all** samples; exclusions live in the meta
  workbook's `exclusions` sheet and are applied only at export.
- Volumes / organ weights are calculation inputs, not stored. DB gets raw
  counts, frequencies and calculated values.

## Stages

### 1. Prep / QC loop + final export -- built, QC checks missing

**Goal:** turn a FlowJo workspace into QC-approved data in the DB with one
switch (`final <- TRUE`).

- [x] Prep script template, `.rds` save/load, exclusion matching,
      `db_export_facs()` (overwrite prompt, one transaction, logging).
      Merged at `e15ed9b`. Design / plan: `2026-11-12-facs-prep-qc-export*`.
- [ ] **QC checks in the report** -- the template is a skeleton; you supply
      the checks / plots that decide "QC passed".

### 2. Workspace versioning

**Goal:** every gate or compensation change is traceable -- which `.wsp`
version is in the DB and what changed between versions.

- [ ] Fix the naming convention and where `.wsp` versions live
      (suggested: `raw/<experiment_code>/facs/`).
- [ ] Per-experiment changelog `.md`.
- [ ] Decide how much to automate: changelog template only, or a `facs_`
      function that diffs two `.wsp` versions and drafts the entry.

### 3. First data in the DB

**Goal:** a real experiment (26-1) in the DB through the full pipeline, plus
the old data brought over. (Schema, writers and snapshots already built --
`2026-11-12-db-setup.md` steps 1-5.)

- [ ] 26-1 import: Treg `.wsp` name, Treg / MPO-Restim population selections,
      then each `facs_prep_26-1_<panel>.R` with `final <- TRUE`.
- [ ] ELISA for 26-1 (hand-made wide sheet; needs plate sample -> mouse map).
- [ ] Migrate old `MPOmRNA.db` (SQLite) into the new schema.

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
