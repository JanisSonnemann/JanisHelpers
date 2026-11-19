# Population display order in the DB, plot helpers, DB-ready calc functions

Builds on `2026-07-17-facs-db-schema-design.md` and
`2026-11-12-facs-prep-qc-export-design.md`. First of two steps: this one
changes JanisHelpers; the second (separate design) rewrites
`analysis_26-1.Rmd` to read from the DB.

## Intent

The user is moving their exploration doc `analysis_26-1.Rmd` (experiment 26-1,
MPO_mRNA) from re-reading workspaces to querying the lab DuckDB. The doc today
is too long and carries quick fixes applied in some chunks but not others
(population order vectors, label fixes, group colours, copy-pasted plots).
Everything that is not specific to 26-1 should live once, in the package, so
the doc's setup chunk only queries, cleans and calculates, and the rest only
plots and text. The next experiment gets the same fixes for free.

## Goals

1. `db_query_facs()` output can be put back into the display order chosen in
   the panel's prep script (`facs_select_populations()`), per panel.
2. Subject meta (including the treatment group) is queried once with a new
   `db_query_subjects()` and joined to each dataset as needed, instead of
   every measurement query carrying meta columns.
3. The two generic plot helpers from the Rmd (`plot_cells()`,
   `plot_population()`) and the copy-pasted paired stim-vs-unstim plot are
   package functions.
4. `facs_calc_log2fc()` / `facs_calc_diff()` run on `db_query_facs()` output.

## Non-goals

- Rewriting `analysis_26-1.Rmd` (next step).
- New group-level treatment table or storing `early` / `reason` / `comment`
  from the meta sheet -- the user confirmed the stored subject meta is enough.
- Moving `compare_groups()` / `run_comparison()`, the ELISA per-mouse
  summary, `show_comparison()`, or the Rmd's `qc_*` helpers into the package.
  `qc_*`, histology ratio calculations and ELISA range parsing are dropped
  from the doc (covered by prep + QC reports + DB); `show_comparison()` stays
  in the doc.
- A separate populations table (approach B below).
- Storing display labels anywhere other than the existing `population` column.
- Renaming `subjects.treatment_group` in the DB. It was renamed away from
  `group` on purpose (SQL reserved word, see the schema spec); only the
  query output calls it `group`.

## Anti-goals

- An existing database must not need recreating: the schema change is
  additive. `db_old_schema_()` must not start rejecting today's DB.
- Stored values must not change on re-export -- only the new column is added.
- No back-filling of the column outside the logged export path (`import_log`
  must describe what is in the DB).

## Approaches (population order)

- **A. `population_order` column on `facs_measurements` (chosen).** Integer
  position of the population in that panel's display order, written by the
  export on every measurement row. Positions are only meaningful within one
  panel (`wsp_root`); the same label can sit at different positions in
  different panels (e.g. PMN early in one panel, after T cells in another).
  Smallest change; repeated values cost nothing in DuckDB's columnar storage.
- **B. New `facs_populations` table** (one row per `wsp_root` x
  `population_full_path`: label + display order), joined in
  `db_query_facs()`. Cleaner model and would let labels live once per panel,
  but a new table, a new join and more rows to replace in the export
  transaction. Rejected: same output for the user, more moving parts.
- Back-filling from the prep `.rds` files instead of re-exporting: rejected
  (bypasses `import_log`, see anti-goals).

## Decisions

### 1. Schema: `facs_measurements.population_order`

- New nullable `INTEGER` column `population_order`. New databases get it in
  the `CREATE TABLE`; existing ones via
  `ALTER TABLE facs_measurements ADD COLUMN IF NOT EXISTS population_order INTEGER`
  in `create_schema_()` (idempotent, runs on every `db_connect()`). Existing
  rows stay `NULL` until their panel is re-exported.
- Meaning: 1-based position of the row's `population` label in the levels of
  the prep's `population` factor. Rows of one stain/root share one numbering.
- `NULL` when the prep's `population` is not a factor (no ordering known).
- Update the schema spec (`2026-07-17-facs-db-schema-design.md`) and
  CLAUDE.md's `db_` entry.

### 2. Export writes the order

- `insert_facs_stain_()` takes an optional `population_order` column in its
  `data` and inserts it. `db_export_facs()` derives it from the prep data's
  `population` factor (`as.integer(population)`) before `as.character()`
  drops it. The factor's levels are the labels (after
  `facs_select_populations(labels = )`), so label and position match.
- Legacy `db_write_facs()` does the same when `data$Population` is a factor,
  else writes `NULL`.

### 3. Queries: `population_order` in FACS, new `db_query_subjects()`

- `db_query_facs()` adds `population_order` to its `select()`. The ELISA and
  histology queries are unchanged. No measurement query gets meta columns.
- **`db_query_subjects(con)`** (new, `R/db_query.R`): lazy `tbl` of
  `subjects` joined to `experiments`, one row per mouse, columns
  `experiment_code`, `mouse_id`, `group` (= `subjects.treatment_group`,
  renamed in the select; dbplyr quotes the reserved word), `mouse_treatment`,
  `sex`, `cage`, `mouse_strain`, `generation`, `dob`, `start_date`,
  `bmt_date`, `end_date`. No surrogate IDs, no `created_at`.
- Intended workflow (the analysis doc's setup chunk, next step): query
  `db_query_subjects()` once, filter to the experiment, `collect()`, set the
  group factor levels once; then per dataset join only the needed columns
  with the existing checked join
  `meta_annotate(data, subjects_cols, by = c("experiment_code", "mouse_id"))`
  (errors on missing keys / name collisions, warns on unmatched rows).
  Connect with `DBI::dbConnect(duckdb::duckdb(), path, read_only = TRUE)`,
  not `db_connect()` (which writes the schema). Document this in the
  `db_query_subjects()` roxygen `@details` / example.

### 4. `facs_order_populations(data)` (new, `facs_` domain, `R/facs_calc.R`)

Turns a collected `db_query_facs()` result's `population` into a factor.

- Requires columns `population`, `population_order`, `wsp_root`.
- Single panel: levels = labels sorted by `population_order`.
- Several panels: panels in alphabetical `wsp_root` order, within each panel
  by `population_order`; a label already placed by an earlier panel keeps that
  position (no error -- different positions across panels are expected).
- Rows with `NA` order (panel not yet re-exported): their labels go after all
  ordered ones, alphabetically, with a message naming the affected `wsp_root`s.
- Returns `data` with `population` as factor; other columns untouched.

### 5. `facs_calc_log2fc()` / `facs_calc_diff()` accept DB column names

- New argument `mouse_col = "mouse_ID"` used everywhere the helper
  `calc_restim_proportions_()` hard-codes `mouse_ID`. DB data is called with
  `mouse_col = "mouse_id", restim_col = "condition"`.
- The passthrough-column detection must not break on the extra DB columns
  (`experiment_code`, `collected_at`, `assay_name`, `wsp_root`,
  `condition_keyword`, `source_file`, `imported_at`, `population_order`,
  joined `group`): columns that are constant per mouse x tissue pass
  through, others drop -- as today. Add a test on a DB-shaped tibble.
- A `population` factor survives (as it does today via `pivot_wider` id cols).

### 6. Plot helpers (`plot_` domain, new file `R/plot_box.R`)

Column arguments are bare names (tidy eval, `{{ }}`), so existing Rmd call
sites keep their shape. Returned objects are `ggplot`s the caller extends with
`+ labs()`, `+ facet_*()`, etc. All use `theme_bw()`.

- **`plot_box(data, x, y, fill, scale = c("linear", "percent", "log"), colors = NULL)`**
  -- the Rmd's `plot_cells()`: dodged boxplot (no outlier points) + jittered
  points (`position_jitterdodge`), x labels dodged over two rows.
  `scale = "percent"` -> percent y labels; `"log"` -> `scale_y_log10()` with
  log ticks (no offset added here). `colors` = named vector for
  `scale_fill_manual()`; `NULL` = ggplot default fill scale.
- **`plot_population(data, populations, metric, tissues = NULL, x = tissue, fill = group, facet = population, scale = "linear", p_val = FALSE, colors = NULL, pseudocount = 1)`**
  -- the Rmd's `plot_population()`: filters `metric`, `population`, optional
  `tissue`; for `scale = "log"` adds `pseudocount` to `value`; calls
  `plot_box()`; `facet_wrap(free_y)` by `facet`; `p_val = TRUE` adds
  `ggpubr::stat_compare_means(label = "p.format")`. `fill` defaults to `group`
  (the column `db_query_subjects()` returns and the Rmd already uses).
- **`plot_paired()`** -- replaces the five copies of the paired plot.
  Signature:
  `plot_paired(data, condition = condition, levels, fill = group, pair = mouse_id, facet = population, y = value, scale = "percent", colors = NULL)`.
  Filters `condition` to `levels` and orders it by them; x =
  `interaction(condition, fill)`; boxplot (alpha 0.6, no outliers) filled by
  `fill`; grey lines joining each `pair` across conditions; points shape 21;
  `facet_wrap(free_y, nrow = 2)`; percent y labels for `scale = "percent"`.
- **Dependency:** add `ggpubr` to `Imports` (DESCRIPTION + CLAUDE.md list).

### 7. Re-export (user, after merge)

The user reruns the three 26-1 FACS panel scripts (`overview-prep_26-1.R`,
`treg_prep_26-1.R`, `mpo-restim-facs_prep_26-1.R`) with `final <- TRUE` and
confirms each overwrite. The preps are reloaded from `processed/`, not
re-read. Check afterwards: row counts and values per root unchanged vs. the
current imports (6, 2, 3 -- see `import_log`), `population_order` non-NULL.

## Testing

`testthat` edition 3, temporary DuckDB via `helper-db.R`.

- `db_connect()` on a DB whose `facs_measurements` lacks the column adds it,
  keeps rows (`NULL` order), and `db_old_schema_()` returns nothing.
- `db_export_facs()` stores `population_order` matching the prep factor;
  character `population` -> `NULL`.
- `db_query_facs()` returns `population_order`; `db_query_subjects()` returns
  one row per mouse with `group` = stored `treatment_group`, filterable lazily.
- `facs_order_populations()`: single panel; two panels with a shared label at
  different positions (first panel wins, no error); `NA` order -> message.
- `facs_calc_log2fc()` / `facs_calc_diff()` with `mouse_col = "mouse_id"`,
  `restim_col = "condition"` on a DB-shaped tibble give the same numbers as
  the CamelCase input.
- Plot helpers return `ggplot` objects; `plot_population()` filters correctly
  and offsets for log; `plot_paired()` drops conditions not in `levels`.
- `devtools::check()`: 0 errors, 0 warnings.

## Left to the builder

- Exact message wording; internal helper names.
- Whether `plot_paired()` hides x-axis text (the Rmd did in some copies) --
  leave it to the caller via `+ theme()`.
- Roxygen examples: `\dontrun{}` where DB/files are needed.
