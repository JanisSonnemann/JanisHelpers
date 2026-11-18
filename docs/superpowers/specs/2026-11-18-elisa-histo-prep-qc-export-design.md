# ELISA + histology prep -> QC loop -> database export (with timepoints)

Status: draft, pending review
Parallel to `2026-11-12-facs-prep-qc-export-design.md` (the FACS loop). Changes
the sample model of `2026-07-17-facs-db-schema-design.md` and the ELISA /
histology parts of `2026-11-12-db-schema-import-log-design.md`.

## Why

Today every analysis reads the ELISA and histology Excel files again by hand,
and nothing checks them. The database should hold these results once, curated
and checked, so analysis starts with a query -- the same single source of truth
the FACS pipeline already gives. Anti-MPO titers are measured on several dates
per mouse, and the database has no notion of a timepoint yet.

## Goals

- One per-assay **prep script** per curated ELISA or histology file, with the
  same `final` switch as the FACS panel scripts:
  - `final <- FALSE`: read the curated sheet, run optional derivations, save a
    prep result and render a QC report (all samples, exclusions only flagged).
  - `final <- TRUE`: load the saved prep result, apply the meta workbook's
    `exclusions` sheet, replace any earlier import of the same file in one
    transaction, log it (`qc_passed = TRUE`), record the exclusions, snapshot.
- A **simple QC report with outlier plots** for each domain.
- **Timepoints in the database**: a sample is one mouse x tissue x collection
  date. Serum titers on five dates are five samples.
- Analysis gets ELISA / histology values (with `collected_at`) from
  `db_query_elisa()` / `db_query_histo()`.

## Non-goals

- Parsing raw instrument exports or arbitrary hand-made layouts. The pipeline
  reads one **curated layout per domain** (below); you bring each file into it.
  `elisa_read_results()` (Bio-Plex `Results_*.xlsx`) stays as it is, outside
  this pipeline.
- Automating the QC decision -- `final <- TRUE` is your "QC passed", as for FACS.
- Timepoints for FACS. FACS stays end-point only (its samples get the mouse's
  end date); a dated FACS panel is a later change.
- Storing substituted values for out-of-range results (e.g. LLOQ / 2). The
  database stores `NA` + the status; substitution is an analysis choice.

## Anti-goals

- Two imports of the same file both in the database (today's append-only
  writers make that possible as soon as a value is corrected).
- A sample silently disappearing from the QC view because it is excluded.
- An exclusion row that matches nothing while the import goes ahead.
- The database receiving data that differs from what the last report showed.

## Decisions (yours)

- Import only curated data that passed a simple QC check; analysis reads the
  database.
- The QC report includes **outlier plots**.
- **Timepoints (option 1A):** `samples` becomes unique on
  `subject_id x tissue x collected_at`, `collected_at` required. End-point
  samples (all FACS, histology, restim ELISA) get the subject's `end_date`
  (filled from `death_date` by `db_write_meta()`).
- **Separate pipelines per domain (option 2B):** ELISA and histology each get
  their own exported prep/run/load functions, QC template and export function,
  mirroring the FACS structure, rather than one generic "table assay" pipeline.

## Approaches considered

Timepoints:
- **1A (chosen):** collection date on the sample. Matches how the data is
  recorded (titer columns are dates); days since BMT / start come from
  `subjects`; works later for dated FACS. Cost: sample lookups change
  everywhere and existing rows need a date (see Migration).
- 1B: free-text timepoint label on the sample (`'end'`, `'d14'`), date
  optional. FACS untouched, but labels drift between experiments and need a
  lookup table for cross-experiment queries.
- Dropped: timepoint only on ELISA rows -- FACS could never use it.

Pipeline structure:
- 2A: one shared curated-table pipeline, ELISA / histo differ only in reader
  and checks. Less code, one place to fix.
- **2B (chosen):** a copy of the FACS structure per domain. Each is readable on
  its own; fixes may need doing twice. *(My call within 2B: identical
  low-level internals -- exclusion matching, the replace-import transaction,
  the whisker-outlier rule, report folder paths -- are shared internal
  helpers. Everything a user calls or sees -- run/save/load, checks, report
  template, export -- is per domain.)*

## Curated layouts *(my call -- you agreed to curate; the exact layout is mine)*

One `.xlsx` per measurement set, named `<root>.xlsx` (e.g.
`spleen-restim-48h_26-1.xlsx`, `anti-mpo-titer_26-1.xlsx`, `histo_26-1.xlsx`).
The **root** (basename without extension) is the import unit: the database
never holds two imports of one root. No date suffix -- the file is edited in
place; the prep staleness check compares file times as for the `.wsp`.

**ELISA** -- sheet `values`, long by sample, wide by analyte:

| column | required | meaning |
|---|---|---|
| `mouse_ID` | yes | |
| `tissue` | no | absent / empty = script's `tissue` argument (e.g. `"spleen"`, `"serum"`) |
| `date` | no | collection date; absent / empty = the mouse's `end_date` |
| one column per analyte | >= 1 | numeric, or `OOR<` / `OOR>` |

Repeated `mouse_ID x tissue x date` rows are technical replicates, numbered
in row order (confirmed by you for `spleen-restim-48h_26-1.xlsx`: the two
rows per mouse are duplicates, not conditions -- so ELISA has no `condition`
dimension). Numeric values outside LLOQ / ULOQ (the restim sheet has
extrapolated values below LLOQ) are stored as they are with status `ok`; the
report counts them and `lloq` / `uloq` travel with each row *(my call)*. Optional sheet `range`: first column holds the row labels
`LLOQ`, `ULOQ` and optionally `unit`; one column per analyte (this is the
layout `spleen-restim-48h_26-1.xlsx` already has). The titer file needs
reshaping: one row per mouse x date, analyte column e.g. `anti-MPO_OD`; the
end-point bleed (`15/16.6.2026`) gets an empty `date`.

**Histology** -- sheet `values`: `mouse_ID`, optional `tissue` (default from
the script, e.g. `"kidney"`), optional `date`, one numeric column per score
(`gloms`, `crescents`, `necrosis`, `abnormal`, ...). Other columns (e.g. the
`ID` slide label in `histo_26-1.xlsx`) are ignored *(my call)*. An all-empty
row (e.g. `26-2-01`) is kept as missing values so completeness shows it.

Readers: `elisa_read_curated(path, tissue)` and `histo_read_curated(path,
tissue)` -> long tibble `mouse_ID`, `tissue`, `collected_at`, `analyte` /
`metric`, (`replicate`), `value`, (`result_status`, `unit`, `lloq`, `uloq`).
`collected_at` stays `NA` when the sheet has no date; it is resolved to
`end_date` against the database at export and against `mice$death_date` in
the QC report.

## The loop (per domain, mirroring FACS)

```
prep script (final <- FALSE)
  read curated .xlsx -> prepare() (optional derivations)          all samples
  -> annotate group / batch from meta -> save prep (.rds) -> render QC report
  review -> fix the curated file / edit exclusions sheet -> rerun   (repeat)
prep script (final <- TRUE)
  load saved .rds (refuse if older than the .xlsx or report missing)
  -> apply exclusions -> replace root (ask if stored) -> import_log
  + import_exclusions -> snapshot
```

`prepare` is a function like FACS's, e.g. for histology
`\(d) histo_calc_ratio(d, "crescents", "gloms", "pct_crescents")` -- derived
values are calculated before export and stored as extra metric rows, as for
FACS. Which helper functions exist (ratio at least) is left to the builder.

Files: `<domain folder>/processed/<root>_prep.rds` and
`<domain folder>/qc/<root>_qc.html`, where the domain folder is the folder
holding the curated file (`raw/ELISA/`, `raw/histology/`).

## QC report contents

Both domains, from the prep `.rds`:
- Overview: samples per tissue x date, file and prep timestamps.
- Completeness against `expected_mice` (script argument; same `early`
  convention as FACS) per tissue x date; `mouse_ID`s not in the meta sheet.
- Empty values grid (analyte / metric x tissue x date).
- Exclusion flags (rows that would be dropped; non-matching rows as warnings).
- **Outlier plots**: per analyte / metric, values by `group_col`, boxplot with
  points outside the whiskers labelled by `mouse_ID` (the FACS rule: whiskers
  computed on the plot axis from the same points; `ggrepel` if installed). For
  multi-date data, one panel per date plus a per-mouse line plot over time.
  A ranked flag table of the labelled points. Optional `batch_col` plots.

ELISA adds: per analyte counts of `OOR<` / `OOR>` and of values below LLOQ /
above ULOQ; replicate agreement -- samples whose replicate CV exceeds
`max_cv` (script argument, default 20 % *(my call)*), with the replicate
values listed. Outlier plots use the replicate mean; log scale optional.

Histology adds: user-supplied sanity rules as strings evaluated per row, e.g.
`rules = c("crescents <= gloms", "necrosis <= crescents")`; violations listed.

## Database changes

**Samples (1A).** `collected_at DATE NOT NULL`, `UNIQUE(subject_id, tissue,
collected_at)`. `db_write_samples()` fills a missing `collected_at` from the
subject's `end_date` and errors listing subjects without one. Sample lookups
(`lookup_sample_id(s)_`) take `collected_at` with the same default.
`db_export_facs()` resolves its samples at `end_date`.

**`elisa_measurements`** (recreated):
`sample_id`, `assay_id NOT NULL`, `import_unit NOT NULL` (the root),
`analyte` (renamed from `cytokine` -- titers aren't cytokines), `replicate
NOT NULL`, `value`, `unit`, `result_status` (`ok` / `OOR<` / `OOR>`),
`source_file`, `imported_at`; `UNIQUE(sample_id, assay_id, analyte,
replicate)`. `sample_id_raw` is dropped (curated sheets carry `mouse_ID`).
LLOQ / ULOQ are stored per row as `lloq`, `uloq` *(my call: they can change
between plates; storing them lets analysis substitute below-range values
later)*.

**`histo_measurements`** (recreated): `sample_id`, `assay_id NOT NULL`,
`import_unit NOT NULL`, `metric`, `value`, `source_file`, `imported_at`;
`UNIQUE(sample_id, assay_id, metric)`.

**`import_exclusions`**: add `collected_at DATE`. The `panel` column holds the
assay name for every domain (documented as "panel / assay").

**Exclusions sheet**: `panel` matches the ELISA / histology assay name too;
new optional column `collected_at` (empty = all dates). `date` keeps its
meaning (when the exclusion was decided). `condition` / `file_name` are
ignored for these domains, and an applicable row that sets them is an error.

**Writers.** New `db_export_elisa(con, prep, exclusions, ...)` and
`db_export_histo(...)` -- same contract as `db_export_facs()` (checks first,
write nothing on failure; overwrite prompt; one transaction: delete the root,
insert, `import_log` with `supersedes_import_id`, exclusion rows). Export
errors if one `mouse_ID x tissue x date x analyte x replicate` occurs twice.
`db_write_elisa()` / `db_write_histo()` (append-only) are removed *(my call:
their skip-on-conflict behaviour is the anti-goal above)*.

**Queries.** `db_query_elisa()` / `db_query_histo()` / `db_query_facs()` add
`collected_at` and `import_unit` / `wsp_root`; ELISA renames `cytokine` ->
`analyte`. Days since BMT etc. stay in analysis (the joined `subjects`
columns are available).

**Migration** *(my call, same pattern as the `wsp_root` schema change)*:
`db_connect()` detects the old `samples` / measurement layout and errors
asking to recreate the database file. Recreating means rerunning the FACS
panel scripts with `final <- TRUE` (preps are saved, so no re-read) -- for
26-1 that is overview, Treg, MPO-Restim -- then the new ELISA / histo scripts.
Take a snapshot first.

## Package pieces

- `R/elisa_read.R`: `elisa_read_curated()`.
- `R/elisa_prep.R`: `elisa_prep_save()`, `elisa_prep_load()`,
  `elisa_prep_run(final, path, assay, tissue, prepare, mice, exclusions, con,
  group_col, batch_col, expected_mice, max_cv, ..., data_dir, snapshot_root,
  overwrite)`.
- `R/elisa_qc.R`: internal check / plot helpers.
- `R/histo_*.R`: same set (`histo_read_curated()`, `histo_prep_*()`,
  `histo_calc_ratio()`, rules check). **Needs a new `histo_` domain in
  CLAUDE.md** (the naming rule requires declaring it first).
- Templates `inst/rmarkdown/templates/elisa-qc-template/`,
  `histo-qc-template/`; knit wrappers `report_knit_elisa_qc()`,
  `report_knit_histo_qc()`.
- Script templates `inst/scripts/elisa_prep_template.R`,
  `inst/scripts/histo_prep_template.R`, shaped like `facs_prep_template.R`.
- Schema / writers / queries in `R/db_*.R` as above; schema spec and
  CLAUDE.md data-structures section updated.
- Tests mirror sources; small curated fixtures in `tests/fixtures/`.

## Project side (outside the package)

Reshape `anti-mpo-titer_26-1.xlsx` into the curated layout (keep the original
alongside); check `spleen-restim-48h_26-1.xlsx` (add `tissue` or rely on the
script default) and both `histo_*.xlsx`. One script per file in the
experiment's `src/`, e.g. `elisa_prep_26-1_spleen-restim.R`,
`histo_prep_26-1.R`.

## Left to the builder

- Exact argument order and internal helper names.
- Report layout beyond the listed sections; plot styling.
- Message wording; how violations and unmatched rows are printed.
- Which `histo_calc_*()` helpers beyond a ratio.
