# FACS prep -> QC loop -> final export

Status: draft, pending review
Detailed design for Stages 2 + 3 of
`2026-11-12-facs-pipeline-workflow-design.md` (the big-picture plan).

## Why

Today `db_import_26-1_<panel>.R` reads a workspace and writes it straight to
the database. In practice, reaching the final data takes several rounds: read
the workspace, clean it, calculate, look at a QC report, adjust gates in
FlowJo, and repeat. The database should get the data only once, at the end,
and it should be exactly the data the last QC report showed. Excluded samples
get dropped at that point, and the database keeps a record of what was left
out and why.

## Goals

- One per-panel **prep script** handles the whole loop. It reads the `.wsp`,
  selects populations, calculates derived metrics, saves the result and
  renders the QC report. You rerun it after every gate change.
- The prep result and the QC report always show **all** samples. Exclusions
  are only *flagged* in the report.
- A **`final` switch** in the same script runs the database export. The export
  applies the exclusions, writes the data and logs the import together with
  the exclusions that applied.
- The export writes exactly the data from the saved prep result that the last
  report was rendered from. It never re-reads the workspace.

## Non-goals

- Automated QC checks. The report is a skeleton, and you supply the checks
  later (unchanged from the parent spec).
- Removing samples during prep or QC.
- Any change to how ELISA or histology data is imported.

## Anti-goals

- A sample disappearing from the QC view because it's on the exclusion list.
- An exclusion row that silently matches nothing (a typo or a renamed file)
  while the import still goes ahead.
- The database receiving data that differs from what the last report showed.

## The loop

```
prep script (final <- FALSE)
  read .wsp -> select populations -> facs_calc_*()        all samples
  -> save prep result (.rds)
  -> render QC report from that .rds   (flags rows on the exclusions sheet)
        |
  review report -> adjust gates in FlowJo (new dated .wsp if major)
                -> edit exclusions sheet
                -> rerun                                   (repeat)
        |
prep script (final <- TRUE)
  load the saved .rds + report for this .wsp version   (no prep, no render)
  -> apply exclusions sheet -> db export (overwrite prompt, one transaction)
  -> import_log row + exclusion rows -> snapshot
```

## Decisions (yours)

- Prep, calculations and the QC report all live in one per-panel R script.
  You run it interactively and rerun it for each iteration.
- The prep result and the QC report contain every sample.
- Exclusions are recorded in the experiment's meta workbook and applied only
  at the final export. The export also stores them in the database.
- An exclusion can be a whole mouse, or anything narrower down to a single
  `.fcs` file.
- The final import is triggered by a `final` switch in the prep script, not
  by a separate script.

## Exclusions sheet

This is a sheet named `exclusions` in the experiment's meta workbook (e.g.
`meta_26-1.xlsx`), read by `meta_read()`. There's one row per exclusion:

| column | required | meaning |
|---|---|---|
| `mouse_ID` | yes | mouse |
| `tissue` | no | empty = all tissues |
| `panel` | no | assay name as used in the prep script (e.g. `overview`); empty = all panels |
| `condition` | no | well condition (e.g. `MPO`); empty = all conditions |
| `file_name` | no | one `.fcs` file; empty = all files |
| `reason` | yes | why *(my call: an exclusion without a reason defeats the point)* |
| `date` | no | when decided |

Matching rules *(my calls)*:
- A row applies to the panel being exported if `panel` is empty or equals that
  panel's assay name. Rows for other panels are ignored.
- An applicable row excludes every file whose values equal all of the row's
  non-empty columns. The comparison is exact, with no case folding.
- If an applicable row matches **no** file in the prep result, the export
  stops and lists every such row. Nothing is written.
- Missing required columns or an empty `reason` also stop the export.
- If the sheet is absent, nothing is excluded. If the sheet has no rows, the
  export goes ahead with a message.
- The QC report runs the same matching to flag rows. There, a non-matching
  row is only a warning shown in the report, not an error.

## The prep result

An `.rds` file holding a named list:
- `data`: the long tibble from `facs_read_wsp()` after population selection
  and calculations. It holds raw plus calculated metrics for every sample, one
  row per `file_name` x `population_full_path` x `metric`, with `mouse_ID`,
  `tissue` and (if used) the condition keyword as columns.
  Compensation/FMO files without `mouse_ID` or `tissue` are dropped.
- `info`: the `.wsp` path, `wsp_root`, version date, assay name, condition
  keyword, created-at time and JanisHelpers version.

It is saved as `<exp>/processed/facs/<wsp basename>_prep.rds`, with one file
per `.wsp` version that gets overwritten on rerun *(my call)*. The QC report
goes to `<exp>/qc/facs/<wsp basename>_qc.html`, which is also one per version
and overwritten on rerun *(my call)*. The report embeds the prep file's
created-at time. Separate folders keep `raw/` limited to raw files.

Prep also runs the **one-file-per-mouse x tissue x condition check** that is in
`db_import_26-1_common.R` today. In prep it is a warning, shown in the report,
so you can fix the keyword in FlowJo during the loop. At export it is an error.

## `final <- TRUE`

- Skips prep and rendering. It loads the prep `.rds` and report path for the
  `.wsp` version configured in the script. If either is missing, or the `.rds`
  is older than the `.wsp` file (gates changed since the last prep run), it
  stops and tells you to run with `final <- FALSE` first *(my call)*.
- Setting `final <- TRUE` is your "QC passed". The export records
  `qc_passed = TRUE`. There's no extra QC prompt.
- The export follows the parent spec's Stage 3: register samples and assay,
  then check whether the `wsp_root` is already in the database. If it is,
  `utils::askYesNo()` asks before overwriting, and "no" aborts with no
  changes. Inside one transaction it deletes the root's old stains and
  measurements, inserts the new ones, writes the `import_log` row (with
  `supersedes_import_id` on an overwrite) and writes the exclusion rows.
- It warns if a newer dated `.wsp` of the same root exists on disk.
- A dated snapshot follows, as today.

## Package pieces

New functions follow the CLAUDE.md naming rules. Final names get settled in
the plan.

1. **Exclusion matching** (`facs_` domain, e.g. `facs_match_exclusions(data,
   exclusions, panel)`). Returns the applicable exclusion rows with the files
   each one matched, plus the unmatched rows. The QC report and the export
   both use it.
2. **Prep save/load helpers** (e.g. `facs_prep_save()` / `facs_prep_load()`).
   They write and read the `.rds` with its `info`, choose the standard paths
   and run the staleness check.
3. **QC report template** `inst/rmarkdown/templates/facs-qc-template/`. It is
   parameterized by the prep `.rds` path and the meta workbook path. It
   contains an overview of samples and files, the exclusion flags, the
   keyword-clash warning and placeholder sections for your checks. It is
   rendered by a thin `report_` wrapper (e.g. `report_knit_facs_qc()`) to the
   path above, and that wrapper returns the path.
4. **DB export** (`db_` domain, e.g. `db_export_facs(con, prep, exclusions,
   qc_report)`). This is the Stage 3 export with exclusions added. It replaces
   the per-file `db_import_facs_wsp()` loop in `db_import_26-1_common.R`,
   reuses `db_write_facs()` internals and `write_import_log_()`, and returns
   the `import_id` invisibly.
5. **Schema: new `import_exclusions` table.** Columns: `import_id` (FK to
   `import_log`), `mouse_id`, `tissue`, `panel`, `condition`, `file_name`,
   `reason`, `decided_on`, `n_files_matched`. Rows are kept when the root is
   overwritten later, so each export's exclusions stay with its log row. The
   FK points at another table, so `EXPORT DATABASE` is fine. Created
   idempotently by `db_connect()`, with an update to the schema spec and to
   CLAUDE.md's data-structures section.
6. **Script template** `inst/scripts/facs_prep_template.R`: config (wsp,
   assay, condition keyword, `final`), population selection, calculations,
   then the prep/report or export branch. A copy of this is what each panel's
   script starts from.

## Project scripts (outside the package)

Rework `MPO_mRNA/data/26-1_BMT-mRNA-V4/src/db_import_26-1_{overview,treg,
mpo-restim}.R` into this shape, starting with overview, which is already
partly done. Rename them `facs_prep_26-1_<panel>.R` *(my call)*. Keep
`db_import_26-1_common.R` for config, `db_open_26_1()` and the snapshot, and
drop its import loop. Fix what the current overview draft still has open:
`ov_path` is undefined (should be built from `wsp`), `mutate` is bare,
`exclude` and `meta$calculations` are undefined, and the `exclude` filter has
to go, since exclusions now come from the sheet. Add the `exclusions` sheet to
`meta_26-1.xlsx` (header only; your entries).

## Left to the builder

- Exact function names and argument order within the domains above.
- QC report layout beyond the required sections.
- Message wording and how unmatched or clashing rows are printed.
