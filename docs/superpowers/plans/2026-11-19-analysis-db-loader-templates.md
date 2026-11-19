# Experiment loader, analysis templates, 26-1 DB port -- Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Load one experiment from the lab DuckDB into a single named list, give the box plots a document-wide colour default, ship two analysis templates (DB and quick `.wsp`), and port `analysis_26-1.Rmd` to the DB.

**Architecture:** New `db_load_experiment()` in `R/db_load.R` composes the existing lazy `db_query_*()` tables, collects per import unit, joins subject columns with `meta_annotate()`, orders FACS populations per panel with `facs_order_populations()`. `db_query_facs()` gains a fallback for DB files without `population_order`. Plot helpers read `getOption("JanisHelpers.group_colors")`. Templates are plain RStudio skeletons. The 26-1 doc lives on OneDrive, outside the repo.

**Tech Stack:** R, dplyr/dbplyr, DBI/duckdb, ggplot2, rmarkdown, testthat 3 + withr.

**Spec:** `docs/superpowers/specs/2026-11-19-analysis-db-loader-templates-design.md`

## Global Constraints

- CLAUDE.md rules: `|>` only; every external call `pkg::fn()`; `dplyr::select(!col)`; non-ASCII as `\uXXXX` in R strings, ASCII-only roxygen; every export has `@param`, `@returns`, `@export`; `\dontrun{}` for examples needing a DB/file.
- New functions belong to declared domains: `db_load_experiment()` is `db_`.
- `db_load_experiment()` never writes to the DB; must work on `DBI::dbConnect(duckdb::duckdb(), path, read_only = TRUE)`.
- No new `Imports`. `ggrepel`/`kableExtra` are used only in the Rmd docs, not the package.
- Reserved list names: `meta`, `exclusions`.
- `devtools::document()` after roxygen changes; `devtools::check()` 0 errors, 0 warnings before the final commit (expected NOTEs per CLAUDE.md).
- The old `analysis_26-1.Rmd` is never modified.
- R lib with DBI/duckdb for running checks/scripts: `R_LIBS=~/Library/Caches/org.R-project.R/R/renv/library/facs-pipeline-workflow-0ff3068e/macos/R-4.6/aarch64-apple-darwin23`. Live DB: `~/databases/MPO_mRNA/MPO_mRNA.duckdb` (open read-only only).

## Review Focus

1. **Experiment with only some domains** (e.g. FACS only, no ELISA/histo) -> list has only the units that exist, no empty elements, no error. Test in Task 2.
2. **Experiment code stripped literally, at the end only** (`26-1` contains `-`; regex use would misbehave; `overview-26-10` must not become `overview0`). Test in Task 2.
3. **Mouse with `NA` group while `group_levels` is given** -> stays `NA`, not an "unknown group" error. Test in Task 2.
4. **Exclusion only in a superseded import, or only tissue-specific** -> mouse stays in `meta`. Test in Task 3.
5. **Explicit `colors = NULL` while the option is set** -> ggplot default fills (the option is only the default for a missing argument), so a single plot can opt out. Test in Task 4.

---

### Task 1: `db_query_facs()` works on DB files without `population_order`

**Files:**
- Modify: `R/db_query.R` (`db_query_facs()`)
- Test: `tests/testthat/test-db_query.R`

**Interfaces:**
- Produces: `db_query_facs(con)` -- unchanged output columns; `population_order` is all `NA` (integer) when the table lacks it.

- [ ] **Step 1: Write the failing test** `"db_query_facs falls back when population_order is missing (read-only)"`: build a temp DB with `db_connect()`, `seed_minimal_fixture()`, `db_write_facs()` one row (as in the existing test), `ALTER TABLE facs_measurements DROP COLUMN population_order`, disconnect, reconnect with `DBI::dbConnect(duckdb::duckdb(), path, read_only = TRUE)` (defer disconnect). Assert `res <- dplyr::collect(db_query_facs(con))`; `expect_true("population_order" %in% names(res))`; `expect_true(all(is.na(res$population_order)))`; `expect_equal(res$population, "CD4+")`.
- [ ] **Step 2: Run** `devtools::test(filter = "db_query")` -- expect FAIL (binder error on `population_order`).
- [ ] **Step 3: Implement**: if `"population_order" %in% DBI::dbListFields(con, "facs_measurements")` is false, add `dplyr::mutate(population_order = NA_integer_)` before the `select()` (translates to `NULL`; cast so the collected column is integer). Note the fallback in the roxygen.
- [ ] **Step 4: Run** `devtools::test(filter = "db_query")` -- PASS.
- [ ] **Step 5: Commit** `fix(db): db_query_facs() reads DB files without population_order`.

### Task 2: `db_load_experiment()` -- list, names, groups, FACS order

**Files:**
- Create: `R/db_load.R`
- Test: `tests/testthat/test-db_load.R` (fixture helper `seed_load_fixture(con)` local to this file)

**Interfaces:**
- Consumes: `db_query_subjects()`, `db_query_facs()` (Task 1), `db_query_elisa()`, `db_query_histo()`, `meta_annotate()`, `facs_order_populations()`.
- Produces: `db_load_experiment(con, experiment, names = NULL, group_levels = NULL, subject_cols = "group", keep_excluded = FALSE)` -> named list `meta`, `exclusions` (Task 3 fills it; this task returns a zero-row tibble with the Task 3 columns), then units. Internal `unit_names_(units, experiment)` -> character.

Fixture (`seed_load_fixture`): experiment `"26-1"` with subjects `26-1-1` (group `MPO409`), `26-1-2` (`irrel`), `26-1-3` (`NA` group); FACS rows via `db_write_facs()` for two roots `overview-26-1` and `Treg_26-1` with a factor `Population` so `population_order` is stored, both containing `CD4+` at different positions; one ELISA row with `import_unit = "anti-mpo-titer_26-1_long"` and one histo row `histo_26-1` (direct INSERTs as in the existing `db_query_elisa` test). A second experiment `"26-10"` with root `overview-26-10`.

- [ ] **Step 1: Write failing tests:**
  - `"returns meta + exclusions + one tibble per import unit in domain order"`: `names(dat) == c("meta", "exclusions", "overview", "Treg", "anti-mpo-titer", "histo")`; `nrow(dat$meta) == 3`; every unit is a collected `tbl_df`; `dat$overview$experiment_code` all `"26-1"`.
  - `"automatic names strip the code literally and only at the end"`: `unit_names_(c("overview-26-1", "Treg_26-1", "anti-mpo-titer_26-1_long", "MPO409-spleen-restim-48h_26-1"), "26-1") == c("overview", "Treg", "anti-mpo-titer", "MPO409-spleen-restim-48h")`; `unit_names_("overview-26-10", "26-1") == "overview-26-10"`; loading `"26-10"` gives element `"overview"`.
  - `"names overrides selected units; unknown unit errors"`: `names = c(titer = "anti-mpo-titer_26-1_long")` -> element `titer`, others automatic; `names = c(x = "nope")` -> `expect_error(..., "nope")` and the message lists available units.
  - `"colliding or empty names error"`: `names = c(overview = "Treg_26-1")` -> error mentioning `overview`; `names = c(meta = "histo_26-1")` -> error.
  - `"unknown experiment errors listing available codes"`: `expect_error(db_load_experiment(con, "99-9"), "26-1")`.
  - `"group_levels sets the factor in meta and joined tibbles; NA group allowed"`: `levels(dat$meta$group) == c("MPO409", "irrel")`; `levels(dat$overview$group)` same; `26-1-3` group `NA` without error; `group_levels = "MPO409"` -> error naming `irrel`.
  - `"subject_cols joins requested columns only"`: default -> `"group" %in% names(dat$histo)`, `"sex"` not; `subject_cols = c("group", "sex")` -> both present.
  - `"FACS population factor is per panel"`: `levels(dat$overview$population)` and `levels(dat$Treg$population)` equal each panel's own stored order (shared `CD4+` at its own position in each).
  - `"only existing domains appear"`: experiment `"26-10"` (FACS only) -> `names == c("meta", "exclusions", "overview")`.
  - `"works on a read-only connection"`: reopen the fixture file read-only and load.
- [ ] **Step 2: Run** `devtools::test(filter = "db_load")` -- FAIL (function missing).
- [ ] **Step 3: Implement** `db_load_experiment()` in `R/db_load.R`:
  - Validate `experiment` against `experiments.experiment_code`.
  - `meta` = `db_query_subjects()` filtered + collected; apply `group_levels` (error on non-`NA` groups not in it).
  - Units: distinct `wsp_root` (FACS), `import_unit` (ELISA, histo) for the experiment; order FACS, ELISA, histo, alphabetical within.
  - `unit_names_()`: remove trailing `_long`, then a trailing `[-_]` + `experiment` matched as a fixed string (e.g. `endsWith()` + `substr()`), no regex on the code.
  - Apply `names`, check unknown / empty / duplicate / reserved (`meta`, `exclusions`).
  - Per unit: filter lazily, `collect()`, `meta_annotate(dplyr::select(meta, experiment_code, mouse_id, dplyr::all_of(subject_cols)), by = c("experiment_code", "mouse_id"))`; FACS additionally `facs_order_populations()` (its message passes through).
  - Roxygen: params, `@returns` list structure, `\dontrun{}` example with read-only connect + `on.exit` disconnect.
- [ ] **Step 4: Run** `devtools::test(filter = "db_load")` -- PASS.
- [ ] **Step 5: Commit** `feat(db): db_load_experiment() loads one experiment as a named list`.

### Task 3: exclusions element and `keep_excluded`

**Files:**
- Modify: `R/db_load.R`
- Test: `tests/testthat/test-db_load.R`

**Interfaces:**
- Produces: `dat$exclusions` tibble with columns `import_unit`, `mouse_id`, `tissue`, `panel`, `condition`, `file_name`, `collected_at`, `reason`, `decided_on`, distinct rows, from the current import (max `import_id` per `import_unit` in `import_log` for the experiment). `meta` drops mice excluded mouse-wide in the current import of every unit unless `keep_excluded = TRUE`; message names dropped mice.

Fixture additions: insert `import_log` rows via `JanisHelpers:::write_import_log_()` for each of the fixture's 4 units, plus a superseded earlier import of `overview-26-1`; insert `import_exclusions` rows directly:
  - `26-1-1` mouse-wide (all scope columns `NULL`) in the current import of all 4 units -> dropped.
  - `26-1-2` mouse-wide in the superseded overview import only, and tissue-specific (`tissue = "kidney"`) in the current imports -> kept.

- [ ] **Step 1: Write failing tests:**
  - `"exclusions lists current-import rows only"`: `dat$exclusions` has no row from the superseded import; `26-1-1` appears once per unit; columns as above.
  - `"keep_excluded = FALSE drops mice excluded mouse-wide everywhere"`: `expect_message(dat <- db_load_experiment(con, "26-1"), "26-1-1")`; `!"26-1-1" %in% dat$meta$mouse_id`; `"26-1-2" %in% dat$meta$mouse_id`.
  - `"keep_excluded = TRUE keeps all subjects"`: `nrow(dat$meta) == 3`.
  - `"no import_log rows -> empty exclusions, all mice kept"`: on the Task 2 fixture without log rows, `nrow(dat$exclusions) == 0`, `nrow(dat$meta) == 3`.
- [ ] **Step 2: Run** `devtools::test(filter = "db_load")` -- FAIL.
- [ ] **Step 3: Implement**: query `import_log` x `import_exclusions` for the experiment (`dplyr::tbl()` + `dplyr::slice_max(import_id, by = import_unit)` or SQL), build `exclusions`; a mouse is dropped when its mouse-wide rows cover every unit name in the experiment's unit list. Filter `meta` before the per-unit joins.
- [ ] **Step 4: Run** `devtools::test(filter = "db_load")` -- PASS.
- [ ] **Step 5: Commit** `feat(db): db_load_experiment() returns exclusions and drops excluded mice`.

### Task 4: group colour option in plot helpers

**Files:**
- Modify: `R/plot_box.R` (`plot_box()`, `plot_population()`, `plot_paired()` defaults and `@param colors`)
- Test: `tests/testthat/test-plot_box.R`

**Interfaces:**
- Produces: `colors = getOption("JanisHelpers.group_colors")` default on all three; explicit `colors` wins; `NULL` (explicit or unset option) -> no manual fill scale.

- [ ] **Step 1: Write failing tests** (use `withr::local_options(JanisHelpers.group_colors = c(a = "#111111", b = "#222222"))` and a two-group tibble):
  - `"plot helpers use the group colour option when colors is missing"`: for each of the three helpers, `ggplot2::ggplot_build(p)$data[[1]]$fill` contains only `"#111111"`/`"#222222"`.
  - `"explicit colors beat the option"`: `colors = c(a = "#333333", b = "#444444")` -> those fills.
  - `"explicit NULL gives ggplot default fills"`: `colors = NULL` -> fills not in the option values.
- [ ] **Step 2: Run** `devtools::test(filter = "plot_box")` -- FAIL.
- [ ] **Step 3: Implement**: change the three defaults; `plot_population()` passes its `colors` on to `plot_box()` unchanged. Document the option in each `@param colors`.
- [ ] **Step 4: Run** `devtools::test(filter = "plot_box")` -- PASS.
- [ ] **Step 5: Commit** `feat(plot): group colours default from JanisHelpers.group_colors option`.

### Task 5: `analysis-db` and `analysis-wsp` templates

**Files:**
- Create: `inst/rmarkdown/templates/analysis-db/template.yaml`, `.../analysis-db/skeleton/skeleton.Rmd`
- Create: `inst/rmarkdown/templates/analysis-wsp/template.yaml`, `.../analysis-wsp/skeleton/skeleton.Rmd`
- Create: `tests/testthat/test-report_templates.R`
- Modify: `CLAUDE.md` (template table, `db_load_experiment()` entry in Data structures, colour option under `plot_box()`)

**Interfaces:**
- Consumes: `db_load_experiment()` (Tasks 2-3), `options(JanisHelpers.group_colors)` (Task 4), `facs_read_wsp()`, `facs_calc_pct_of()`, `facs_calc_count_per_g()`, `facs_select_populations()`, `meta_read()`, `meta_clean()`, `meta_annotate()`.
- Produces: skeletons whose setup leaves only `dat` and `pops` in the global env.

- [ ] **Step 1: Write failing test** `"analysis templates parse and leave only dat and pops"`: for each template, `system.file("rmarkdown/templates/<t>/skeleton/skeleton.Rmd", package = "JanisHelpers")` exists; `knitr::purl(..., output = tempfile(), quiet = TRUE)` then `parse()` succeeds; the purled setup code contains `dat <-` and `pops <-` and (DB template) `read_only = TRUE` and `db_load_experiment(`; no `rstudioapi`.
- [ ] **Step 2: Run** `devtools::test(filter = "report_templates")` -- FAIL.
- [ ] **Step 3: Write the templates** per spec section 3. `template.yaml`: `name`, one-line `description`, `create_dir: FALSE`. Both: YAML header + knitr block copied from `experiment-report`; user settings (`db_path`/`experiment`/`group_levels`, or meta/wsp paths) held inside the `local()` that builds `dat`, with `options(JanisHelpers.group_colors = ...)` at top level. DB template: connection opened and closed inside `local()` via `on.exit(DBI::dbDisconnect(con, shutdown = TRUE))`; commented `dat$<x>_fc <- ... facs_calc_log2fc(mouse_col = "mouse_id", restim_col = "condition")` example. WSP template: rename `mouse_ID` -> `mouse_id` and the restim keyword column -> `condition` after the calc steps. Body: headings FACS / ELISA / Histology with one commented example chunk each.
- [ ] **Step 4: Update CLAUDE.md** (rows `analysis-db` / `analysis-wsp` -> `report_knit_exp()`; Data structures line for `db_load_experiment()`; colour option note).
- [ ] **Step 5: Run** `devtools::document(); devtools::test()` -- all PASS; then `devtools::check()` -- 0 errors, 0 warnings.
- [ ] **Step 6: Commit** `feat(report): analysis-db and analysis-wsp templates`.

### Task 6: `analysis_26-1_db.Rmd` (DB port)

**Files:**
- Create: `/Users/janis/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/Shared Documents/MPO_mRNA/data/26-1_BMT-mRNA-V4/analysis/analysis_26-1_db.Rmd` (outside the repo; not committed)
- Read only: `analysis_26-1.Rmd` in the same folder; `meta/meta_26-1.xlsx`

**Interfaces:**
- Consumes: the installed branch build of JanisHelpers (Tasks 1-5).

- [ ] **Step 1: Install the branch build to a temp lib**: `R CMD INSTALL -l "$TMPDIR/rlib" .` (with the R_LIBS above on the path).
- [ ] **Step 2: Write the doc** from the `analysis-db` skeleton, following the spec's section 4 mapping table section by section in the old order. Setup: `db_load_experiment(con, "26-1", names = c(overview = "overview-26-1", treg = "Treg_26-1", restim = "MPO-Restim_26-1", restim_48h = "MPO409-spleen-restim-48h_26-1", titer = "anti-mpo-titer_26-1_long"), group_levels = c("MPO409", "irrel"))`; colours `c(MPO409 = "#00BFC4", irrel = "#F8766D")`; derived ELISA columns (`oor_low`, `conc`, `LLOQ`, `well`) and `dropouts` (`early`/`reason` from `meta_read()`) stored into `dat`; all population vectors and `mfis` in `pops`. A `helpers` chunk right after setup holds the adapted `qc_*`, `compare_groups()`, `run_comparison()`, `show_comparison()` (and `normalized_log2FC()` only if Step 4 shows `facs_calc_log2fc()` differs). Top section "Changes vs. the wsp version" + "Data gaps" (spec list + anything found in Step 4).
- [ ] **Step 3: Knit against the live DB**: `R_LIBS="$TMPDIR/rlib:<lib above>" Rscript -e 'rmarkdown::render("<doc>", output_dir = Sys.getenv("TMPDIR"))'`. Expected: renders without error; global env after setup holds only `dat`, `pops` and helper functions.
- [ ] **Step 4: Compare with the old doc's data**: purl the old doc's setup chunk (replace the `rstudioapi` fig.path line), run it, and compare per panel for the same mice: one population x tissue x metric per panel (`value` equal within 1e-9), histology `cresc_perc` vs `pct_crescents / 100`, ELISA per-mouse `value` and `res_all` table, restim log2FC from the old local function vs `facs_calc_log2fc()`. Record every difference in the doc's "Data gaps".
- [ ] **Step 5: Report**: list of data gaps found and the rendered HTML path; the doc is not committed (OneDrive).

### Task 7: Final verification and memory

- [ ] **Step 1: Run** `devtools::check()` on the branch -- 0 errors, 0 warnings, only the CLAUDE.md-listed NOTEs.
- [ ] **Step 2: Update memory** `facs-db-import-review-2026-11-17.md`: 26-1 analysis doc ported to `analysis_26-1_db.Rmd`; population order still needs the three FACS re-exports.
- [ ] **Step 3: Commit** any remaining doc changes (`docs: db_load_experiment in CLAUDE.md`), if not already in Task 5.
