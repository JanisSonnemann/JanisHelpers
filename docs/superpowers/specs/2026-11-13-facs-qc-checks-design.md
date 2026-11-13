# FACS QC checks in the per-panel QC report

Stage 1 of `docs/superpowers/plans/2026-11-12-facs-pipeline-overview.md`.
Builds on `2026-11-12-facs-prep-qc-export-design.md`, which left the
"QC checks" section of `facs-qc-template` empty. Source draft: the user's
`qc-report-template.Rmd` (repo root, untracked), checks 5-8 of it.
Approved example: the mockup shown during brainstorming (made-up treg panel).

## Why

The QC report is where the user decides `final <- TRUE` for a panel. Today it
only lists files, keyword clashes and exclusions. It needs to show, at a
glance:

- **single mice or samples that are off** (a bad gate, a failed acquisition,
  low event counts), and
- **day-to-day batch effects** (a staining / acquisition day that is shifted),
  without unbalanced groups faking or hiding a day effect.

A flag means "go look at this in FlowJo", not a verdict. Consistent ranking
matters more than statistical formality.

## Goals

- Four checks in the per-panel report: completeness, gating hierarchy, raw
  counts, MFI.
- Each of hierarchy / counts / MFI gets an **individuals** plot and a
  **batch** plot.
- A final **flagged samples** table ranking mice by how often they are flagged.
- The report template stays experiment-agnostic: everything experiment-specific
  (group column, batch column, expected mice, QC selection) comes from the prep
  script via the prep `.rds`.

## Non-goals

- Experiment-level checks across panels or on the meta workbook (draft checks
  1-4: group x sex overview, missing volumes, dropouts, mice missing a panel).
  Possibly a separate report later; not part of this work.
- Formal statistical tests (Grubbs, Dixon, tests for batch differences).
- Excluding or dropping samples in the report -- it still shows all samples
  and only flags exclusions (unchanged rule).
- Any automatic QC pass/fail. `final <- TRUE` stays the user's decision.

## Anti-goals

- Experiment-specific code in the template (hardcoded group names / colours
  like `MPO409` / `irrel`, hardcoded panel object names).
- An RStudio dependency (`rstudioapi`) -- rendering stays headless via
  `report_knit_facs_qc()`.
- A batch shift showing up as a group difference because groups are
  unbalanced across days.

## Data flow (decided)

The **prep script** prepares everything; the report only reads the `.rds`.

- The prep script joins mouse-level columns (at least the group column and the
  batch column) onto `data` with `meta_annotate()` before `facs_prep_save()`.
- New inputs to `facs_prep_save()`, stored in `info`:
  - `group_col` -- name of the treatment-group column in `data` (e.g. `"group"`).
  - `batch_col` -- name of the batch column (e.g. `"death_date"`; it depends on
    the experiment whether that is the staining / acquisition day, so it is a
    setting, not a fixed column).
  - `expected_mice` -- character vector of mouse IDs the panel should contain:
    meta mice **minus dropouts**. The prep script computes it (the dropout
    column is experiment-specific).
  - `expected_tissues` -- character vector of tissues the panel should contain.
  - `qc_populations` -- populations to show in the QC plots; default all
    populations in `data`.
  - `qc_metrics` -- MFI metrics (markers) to show; default all MFI metrics in
    `data`.
  - `qc_mfi_populations` -- the short list of populations whose MFIs are
    checked (e.g. `"Treg"`). Default: none -- the MFI check is **opt-in per
    panel**, because every population carries a median for every marker and
    scoring all of them would swamp the flag table.
- The template `inst/scripts/facs_prep_template.R` gets the matching settings
  and the `meta_annotate()` join.

## The checks (decided)

All plot checks use the same outlier rule and the same two-plot pattern.

### Outlier rule

- **Robust z-score**: `z = (x - median) / (1.4826 * MAD)`, flag `|z| > 3.5`
  (Iglewicz-Hoaglin).
- Computed **within group x tissue x population x metric** (and condition, for
  multi-well panels).
- **Excluded samples** (exclusions sheet) don't feed a box's median / MAD and
  don't count toward the 4-sample minimum, but still get a z-score and are
  shown and marked.
- MFI values `<= 0` (possible after compensation) are not scored (log10).
- **Minimum events (user's decision):** values from a population with fewer
  than **50 events** (`count < 50`) are not scored -- shown in the plots, but
  never flagged and not counted toward a box's median / MAD. Applies to all
  three checks (fraction, count, MFI of that population). On the Treg data,
  15 of ~100 flags came from gates with 2-50 events (counting noise).
- **Empty gates**: `count = 0` / `fraction_of_parent = 0` occur in rare gates.
  Transforms use a pseudo-count so they stay finite: `log10(count + 1)` and
  logit of `(count + 0.5) / (parent_count + 1)`. (With the 50-event minimum
  these are not scored anyway; the pseudo-count keeps the plots and code
  free of `-Inf`.)
- The **batch column** is treated as a factor (e.g. `death_date` is read as a
  date-time; format to a date label).
- On a transformed scale: **logit** for `fraction_of_parent`, **log10** for
  counts and MFIs.
- Boxes with **fewer than 4 samples** are shown but not scored.
- **Minimum spread**: the MAD is floored so a box whose values happen to sit
  very close together can't inflate z (the mockup flagged ordinary values at
  z = +7 without it).
- Plots stay as in the draft (boxplot + jittered points). Only the label
  choice changes: `mouse_ID` labels go on points with `|z| > 3.5`. Boxplot
  whiskers stay Tukey-style; a labeled point may sit inside the whiskers --
  accepted.

### 5. Completeness

Expected `expected_mice x expected_tissues` vs. mouse x tissue present in
`data`. For multi-well panels the unit is mouse x tissue x condition (expected
conditions = those present in `data`), so a missing well is a missing sample.
Output: a line with expected / present counts and a table of missing
samples, or a "none missing" line.

### Multi-well panels

Each condition gets its own set of plots (stimulated and unstimulated wells
never share a box). Single-well panels (condition `"none"`) look as described
below with no extra level.

### 6. Gating hierarchy -- `fraction_of_parent`

- Individuals: x = group, facets population (rows) x tissue (cols), free y,
  percent axis, logit-based flags.
- Batch: see below.

### 7. Raw counts -- `count`

Same layout as 6, log10 y-axis.

### 8. MFI -- median metrics

One plot per population in `qc_mfi_populations` (skipped with a note if the
list is empty): facets marker (rows) x tissue (cols), log10 y-axis.

### Batch plots (layout C, user's choice)

For each of 6-8, **one plot per tissue**: facets population/marker (rows) x
group (cols), x = batch, boxes and points coloured by batch. Days are compared
**within** each group, so unbalanced groups can't fake an effect. No labels,
no test.

### Flagged samples table (end of the report)

Every flagged value: `mouse_ID`, flag score for that mouse, tissue,
group, **batch**, population, metric, value, z.

**Ranking (user's decision):** mice are ranked by the number of **distinct
sample (mouse x tissue [x condition]) x check** combinations flagged -- check
being hierarchy, counts or MFI -- not by the raw number of flagged values. One
bad parent gate flags itself and its children in fraction and count; that
counts once per sample and check, not ten times. Every flagged value is still
listed. Within a mouse, rows are sorted by `|z|`. The batch column lets the user see when all of a
mouse's flags come from one shifted day (in the mockup, a 50 % brighter day
put batch-driven MFI flags next to a real count problem). Rows for samples on
the exclusions sheet are marked as excluded.

## Report header

Add to the existing header table: group / batch column, and the QC selection
("all (n)" or the list).

## Calls made by Claude (user to check)

- Group and batch column names are settings stored in `info`, not hardcoded.
- Colours are ggplot2 defaults; no fixed group palette.
- Plot / flag helpers live in the package as internal functions in a new
  `R/facs_qc.R` (`facs_` domain), so the template stays short and the scoring
  is unit-testable.
- No new hard dependencies: tables via `knitr::kable` (no `kableExtra`); no
  `ggpubr` (no p-values in QC). `ggrepel` for labels only if installed
  (Suggests), plain `geom_text` otherwise.
- The flagged-samples table marks excluded samples rather than hiding them.
- Multi-well panels: one plot set per condition; completeness per condition.
- Excluded samples left out of median / MAD; MFI `<= 0` not scored.
- Pseudo-counts for empty gates; batch column as a factor; MAD floor starting
  values 0.2 (logit) / 0.1 (log10).

## Real data for tuning

26-1 overview prep result (read-only, OneDrive):
`~/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/Shared Documents/MPO_mRNA/data/26-1_BMT-mRNA-V4/raw/FACS/overview/processed/overview-26-1_20260806_prep.rds`.
41 populations, ~13 mice per tissue, 6 tissues; metrics `count`,
`fraction_of_parent`, `pct_of_live_CD45+`, `count_per_g`. It has **no MFI
metrics** and no group / batch columns yet (prep written before this change)
-- join them from the 26-1 meta workbook in the same `data/` folder for tuning.
Never write into that folder.

26-1 Treg prep result (same folder layout):
`.../raw/FACS/Treg/processed/Treg_26-1_20260828_prep.rds`. 48 files (16 mice x
kidney / lung / spleen), 34 populations, condition `"none"`; MFI metrics
`median_<marker>` for 12 markers (179 values `<= 0`, 77 `NA` in empty gates;
some populations carry MFIs in only 11 samples). Meta `meta_26-1.xlsx`
(sheet `meta`: `group` MPO409 / irrel, `death_date` 2026-06-15 / 06-16,
`early` = dropout). Per group x tissue n = 8; per group x day x tissue n = 3-5.

**Prototype result (2026-11-13, throwaway script):** with no MAD floor 148
flagged values; floor 0.2 (logit) / 0.1 (log10) -> 98, top-ranked mice
stable. Top mouse 26-1-2 (RORgt gates in lung and spleen, z up to 18) -- a
plausible real gating issue, i.e. the method catches what it should.

## Left to the builder

- Exact MAD floor. **Starting values: 0.2 on the logit scale, 0.1 on the
  log10 scale** (from the prototype above); confirm on overview + Treg after
  the 50-event rule and document the final values.
- How MFI metrics are recognised by default: metric names starting with
  `median_` (confirmed on the Treg prep).
- Plot sizing (`fig.height` scaling with number of facets).
- Argument names and defaults for the new `facs_prep_save()` inputs; keeping
  existing calls working where sensible.
- Whether `db_export_facs()` needs to drop the joined meta columns before
  writing (check it tolerates extra columns).
