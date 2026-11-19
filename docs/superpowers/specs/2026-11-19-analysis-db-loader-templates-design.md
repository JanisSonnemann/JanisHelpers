# Experiment loader, analysis templates, DB port of `analysis_26-1.Rmd`

Builds on `2026-11-19-population-order-plot-helpers-design.md` (step 1, merged
at 92f5e3d). This is step 2: the analysis doc reads from the lab DB.

## Intent

The user analyses each experiment in one exploration Rmd. Data now lives in
the lab DuckDB (FACS panels, ELISA, histology, subject meta, exclusions), but
`analysis_26-1.Rmd` still re-reads workspaces/xlsx, recalculates
`count_per_g`, keeps its own exclusion list, and leaves many loose objects in
the workspace (`ov_data`, `treg_data`, `restim_data`, `histo`, `wells`,
`limits`, `meta`, population vectors, ...). The user cares about **workspace
tidiness**: setup should leave as few top-level objects as possible, and the
setup boilerplate should exist once (package + template), not be re-typed per
experiment. A quick workflow straight from `.wsp` files must stay available
for analyses that are not (yet) in the DB.

## Goals

1. `db_load_experiment()` loads one experiment from the DB into **one named
   list**: `meta` + one collected tibble per import unit.
2. Group colours set once per document, used by the `plot_*` box helpers by
   default.
3. Two RStudio templates: `analysis-db` (DB loader setup) and `analysis-wsp`
   (quick `.wsp` workflow), both producing the same `dat` / `pops` shape so
   plot chunks are portable between them.
4. `analysis_26-1_db.Rmd`: a faithful DB port of `analysis_26-1.Rmd` --
   every section and plot reproduced, from the DB.

## Non-goals

- Changing or deleting the existing `analysis_26-1.Rmd` (stays as reference).
- Re-exporting panels or changing prep scripts (user's job; see "Data gaps").
- New analysis content (e.g. the anti-MPO titer series is loaded but not
  plotted -- the port reproduces, it does not extend).
- Moving the doc's ELISA statistics (`compare_groups()`, `run_comparison()`,
  `show_comparison()`) or `qc_*` helpers into the package.
- Changing the `experiment-report` template.

## Anti-goals

- Many top-level objects after setup. Target for the DB template: `dat`,
  `pops` (plus knitr's own). No open connection left behind.
- The loader writing to the DB (must work on a read-only connection, also on
  DB files that predate `population_order`).
- Silent data differences between old and new doc beyond the documented gaps.

## Approaches (data shape)

- **One list, one tibble per import unit (chosen by user).** Per-panel
  population order works (shared labels like `CD4+` keep their own position
  per panel); each tibble has only its relevant columns; chunks change from
  `treg_data$data |>` to `dat$treg |>`.
- One tibble per domain (all FACS together): rejected -- one factor across
  panels breaks per-panel order, every plot needs a `wsp_root` filter.
- One top-level object per assay: rejected -- workspace clutter.

## Decisions

### 1. `db_load_experiment()` (`db_` domain, `R/db_query.R`)

```r
db_load_experiment(con, experiment, names = NULL, group_levels = NULL,
                   subject_cols = "group", keep_excluded = FALSE)
```

- `con`: DBI connection (read-only is the intended use).
- Returns a named list. First element `meta` = collected
  `db_query_subjects()` filtered to `experiment` (one row per mouse;
  excluded mice dropped, see `keep_excluded`). Then
  one element per distinct import unit of that experiment, FACS first
  (alphabetical by `wsp_root`), then ELISA, then histology (alphabetical by
  `import_unit` within domain).
  - FACS import unit = `wsp_root`; ELISA / histology = `import_unit`.
- **Automatic names:** import unit with the experiment code and a trailing
  `_long` removed, plus the separator next to the code (`overview-26-1` ->
  `overview`, `Treg_26-1` -> `Treg`, `MPO-Restim_26-1` -> `MPO-Restim`,
  `MPO409-spleen-restim-48h_26-1` -> `MPO409-spleen-restim-48h`,
  `anti-mpo-titer_26-1_long` -> `anti-mpo-titer`, `histo_26-1` -> `histo`).
  Names that would be empty or collide (with each other or with `meta`) ->
  error listing them.
- **`names`**: named character vector `c(new_name = "import_unit")`
  overriding the automatic name of those units; unknown import units -> error
  listing the available ones.
- **`group_levels`**: if given, `group` becomes a factor with these levels in
  `meta` (before joining, so every tibble inherits it). Groups present but not
  in `group_levels` -> error.
- **`subject_cols`**: columns of `meta` joined onto each tibble via
  `meta_annotate(by = c("experiment_code", "mouse_id"))`. Default `"group"`.
- **`keep_excluded`** (user decision): `FALSE` drops from `meta` every mouse
  that is excluded mouse-wide (exclusion row with `tissue`, `panel`,
  `condition`, `file_name`, `collected_at` all `NULL`) in the *current*
  import of **every** import unit of the experiment (current = latest
  `import_id` per import unit in `import_log`). Partially excluded mice stay.
  Message naming the dropped mice. `TRUE` keeps all subjects. Measurement
  tibbles are unaffected (excluded rows are not in the DB).
- **`exclusions`** element: distinct `mouse_id` x `tissue` x `panel` x
  `condition` x `file_name` x `collected_at` x `reason` x `decided_on` rows of
  the current imports (builder's call: DB wording is the source of truth),
  with the import unit; stored as `dat$exclusions` (name reserved like
  `meta`).
- Each tibble is `dplyr::collect()`ed.
  - FACS: `facs_order_populations()` applied per panel.
  - Columns kept as the `db_query_*()` return them (`mouse_id`, `condition`,
    `collected_at`, ...). No renaming to the old CamelCase.
- **Older DB files without `facs_measurements.population_order`:** the loader
  (or `db_query_facs()`, builder's choice) must not fail on a read-only
  connection: treat the column as all-`NA`, so `facs_order_populations()`
  messages and orders alphabetically. Today's 26-1 DB file is in this state.
- Experiment code not in the DB -> error listing available codes.

### 2. Group colours default

- New package option `JanisHelpers.group_colors` (named character vector).
- `plot_box()`, `plot_population()`, `plot_paired()`:
  `colors = getOption("JanisHelpers.group_colors")` as default. Explicit
  `colors` still wins; option unset = ggplot default fill (as today).
- Documented in each function's `@param colors`.
- `plot_bubble_fc()` unchanged (no group fill).

### 3. Templates

Both under `inst/rmarkdown/templates/`, each with `template.yaml` +
`skeleton/skeleton.Rmd`, same YAML header / knit function
(`report_knit_exp()`) and knitr option block as `experiment-report` (uses
`knitr::current_input()`, no `rstudioapi`).

**`analysis-db`** setup chunk, in this order:
1. knitr options (as experiment-report).
2. `library(tidyverse)`, `library(JanisHelpers)`.
3. User block: `db_path`, `experiment`, `group_levels`,
   `options(JanisHelpers.group_colors = c(...))`.
4. `local({ con <- DBI::dbConnect(duckdb::duckdb(), db_path, read_only = TRUE); on.exit(DBI::dbDisconnect(con, shutdown = TRUE)); db_load_experiment(con, experiment, group_levels = group_levels) })`
   assigned to `dat` -- connection never reaches the global env. The user
   block variables live inside a `local()` / list too (builder's choice how)
   so only `dat` and `pops` remain.
5. Derived data stored back into `dat` (commented example:
   `dat$restim_fc <- dat$restim |> facs_calc_log2fc(..., mouse_col = "mouse_id", restim_col = "condition")`).
6. `pops <- list(...)` for population vectors.
Body: one empty heading per domain (FACS, ELISA, Histology) with one
commented example chunk each using `plot_population()` / `plot_box()`.

**`analysis-wsp`** (quick, no DB):
1. Same knitr / library / colour block.
2. `dat <- local({ ... })` building `list(meta = ..., <panel> = ...)`:
   `meta_read()` + `meta_clean()`, `facs_read_wsp()`, `facs_calc_pct_of()`,
   `facs_calc_count_per_g()`, `facs_select_populations()` for order, join
   `group` with `meta_annotate()`, then rename `mouse_ID` -> `mouse_id` (and
   the restim keyword -> `condition` where present) so chunks match the DB
   template and `plot_paired()`'s default `pair = mouse_id`.
3. `pops <- list(...)`; same body headings.

CLAUDE.md template table gains both rows (knit function `report_knit_exp()`).

### 4. `analysis_26-1_db.Rmd` (OneDrive, next to the old doc)

Path: `.../MPO_mRNA/data/26-1_BMT-mRNA-V4/analysis/analysis_26-1_db.Rmd`.
Built from the `analysis-db` template. Faithful = same sections, same plots
(data selection, metric, tissues, scale, facet, p-values, titles), same
tables, in the same order. Mapping:

| Old | New |
|---|---|
| `ov_data$data`, `treg_data$data`, `restim_data$data` | `dat$overview`, `dat$treg`, `dat$restim` (`names =` override) |
| `histo` (+ `cresc_perc`) | `dat$histo`, metric `pct_crescents` (stored 0-100 -> divide by 100 for the percent axis) |
| `wells` / `limits` (ELISA xlsx) | `dat$restim_48h`: `oor_low = result_status == "OOR<"`, `conc = value`, `LLOQ = lloq`, `well = replicate` (derived in setup, stored back into `dat`) |
| `mouse_ID`, `restimulation` | `mouse_id`, `condition` |
| restim labels `...#` | stored without `#` (e.g. `CD4+`, `CD4+_CD44+_IFNg+`) |
| `plot_cells()` / local `plot_population()` / paired ggplot copies | `plot_box()` / `plot_population()` / `plot_paired(levels = c("unstim", "MPO"))` |
| `normalized_log2FC()` | `facs_calc_log2fc(ref_pop = "CD4+", mouse_col = "mouse_id", restim_col = "condition", ...)` -- builder verifies same numbers as the local function on the same data, else keeps a local function |
| `exclude_reasons` tribble | `dat$exclusions` (DB reason wording, current imports only) |
| `meta$meta` group x sex table | from `dat$meta` |
| dropouts table (`early`, `reason`) | not stored in DB: read those two columns from the meta workbook (`meta_read()`) in setup, stored in `dat` |
| missing volumes / organ weights table | dropped: volumes are prep inputs, checked in the QC reports; replaced by a one-line note |
| `qc_*` helpers (completeness, hierarchy, raw counts, MFI drift) | kept as local functions in the doc, adapted to DB column names; `death_date` -> `end_date` |
| local helper functions generally | stay top-level, all in one `helpers` chunk right after setup (user decision: functions are not counted against the `dat` / `pops` tidiness target) |
| `pop_exclude` / population vectors | in `pops` |
| "without 26-1-16" ELISA sensitivity table | DB has no 26-1-16 (excluded on export) -> chunk kept, replaced by a note that the main table is already without 26-1-16 |

Obvious defects in the old doc are fixed and listed in a short "Changes vs.
the wsp version" section at the top of the new doc:
- `mfis` used before definition -> defined in `pops` in setup.
- Titles that name the wrong tissue (e.g. "Kidney" on lung plots) -> match
  the plotted tissue. **Data selections are not changed**, even where they
  look inconsistent (`lung` vs `lung-whole`, one lung chunk filtering
  `kidney`); these are listed in the section for the user to decide.
- Group colours from the option, not hard-coded per plot.

### Data gaps (documented in the new doc, not fixed here)

Any old-vs-DB difference the builder finds (missing gate, missing metric x
tissue combination, different value) is added to this list in the doc and
the plot is ported with the facet missing -- the port does not stop for it
(builder's call). Known so far: Treg has no spleen `count_per_g` (dropped on
purpose at prep), so the Treg per-g chunks lose the spleen facet.

Gates the old doc plots that the DB does not hold -- their facets are missing
until the user adds them to the prep selection and re-exports:
- overview: `Mac_Siglec-F+`, `Mac_Siglec-F-`, `Mac_kidney-infiltrating`,
  `Mac_kidney-resident`, `Mac_Ly6Chi`, `Mac_Siglec-F+_CD11b+/-`,
  `Mac_Siglec-F-_CD11b+`, `Mac_F4-80+_MHC-II+`, `Alv-Mac_SiglecF+`,
  `Int-Mac_SiglecF-`.
- Treg: `time` (hierarchy QC), `FoxP3-_RORgt+`.
- Population order: alphabetical until the three panels are re-exported with
  master >= 92f5e3d.

## Testing

- `db_load_experiment()` (temporary DuckDB, `helper-db.R`): list names and
  order; `meta` one row per mouse; automatic names incl. `_long` stripping;
  `names` override; collision / unknown unit / unknown experiment errors;
  `group_levels` factor in `meta` and joined tibbles, unknown group error;
  `subject_cols`; FACS population factor per panel; `keep_excluded`
  (mouse-wide exclusion in all units dropped, partial or superseded-import
  exclusion kept, `TRUE` keeps all); `exclusions` element; works on a read-only
  connection to a DB file whose `facs_measurements` lacks
  `population_order`.
- Plot helpers: option used when `colors` missing; explicit `colors` wins;
  unset option -> no manual scale.
- Templates: each skeleton parses (`knitr::purl` + `parse`).
- `devtools::check()`: 0 errors, 0 warnings.
- New doc: knits against the live DB (read-only); builder compares a sample
  of values (one population per panel x tissue x metric, the ELISA table, the
  histology plot data) with the old doc's in-memory data for the same mice.

## Left to the builder

- Internal helper names; error / message wording.
- Whether the missing-column fallback lives in `db_query_facs()` or the
  loader.
- How the template's user block stays out of the global env.
- Chunk labels in the new doc.
