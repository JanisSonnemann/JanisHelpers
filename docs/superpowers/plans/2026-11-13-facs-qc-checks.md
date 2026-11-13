# FACS QC checks Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fill the empty "QC checks" section of `facs-qc-template` with completeness, gating-hierarchy, raw-count and MFI checks, robust-z outlier labels, layout-C batch plots and a ranked flagged-samples table.

**Architecture:** The prep script passes QC settings to `facs_prep_save()`, which stores them in `info$qc`. A new `R/facs_qc.R` holds internal (unexported) helpers: settings resolution, scoring, tables, plots. The template calls them and stays experiment-agnostic.

**Tech Stack:** R, dplyr / tidyr / ggplot2 (already in Imports), knitr::kable, rmarkdown; `ggrepel` optional (Suggests); testthat 3.

**Spec:** `docs/superpowers/specs/2026-11-13-facs-qc-checks-design.md`

## Global Constraints

- CLAUDE.md conventions apply: base pipe `|>`, every external call `pkg::fn()`, `dplyr::select(!col)` not `-col`, non-ASCII in R strings as `\uXXXX`, ASCII-only roxygen.
- No new `Imports`. `ggrepel` goes into `Suggests` and is used only behind `requireNamespace("ggrepel", quietly = TRUE)`; otherwise `ggplot2::geom_text()`.
- No `rstudioapi`, `kableExtra`, `ggpubr`; tables via `knitr::kable`.
- Outlier rule: robust z = `(t - median) / max(1.4826 * MAD, floor)`, flag `|z| > 3.5`; within group x tissue x condition x population x metric; boxes with < 4 scorable samples not scored.
- Scales: logit for `fraction_of_parent`, log10 for counts and MFIs. Pseudo-counts: `log10(count + 1)`; fraction as `(count + 0.5) / (count / fraction + 1)`.
- MAD floor starting values: `0.2` (logit), `0.1` (log10).
- Not scored (but plotted): population `count < 50`; MFI `<= 0` or `NA`; excluded files don't feed median / MAD but get a z.
- MFI metrics = names starting with `median_`. MFI check is opt-in via `qc_mfi_populations` (default none).
- Ranking: mice by number of distinct (tissue x condition x check) combinations flagged.
- Real data is read-only: never write into the OneDrive `MPO_mRNA` folder.
- `devtools::check()`: 0 errors, 0 warnings, no new NOTEs beyond CLAUDE.md's known list.

## Review Focus

1. **A prep `.rds` saved before this change** (overview and Treg 26-1 exist, no `info$qc`) must still render; checks needing a group column show a "rerun prep with `group_col`" note instead of erroring. -> Task 1 (`qc_settings_` defaults), Task 5 (render test).
2. **A box whose values are all identical** (MAD = 0) must give finite z (floor), not `Inf` / `NaN`. -> Task 2.
3. **Fraction exactly 0 or 1, count 0** must never produce `Inf` / `NaN` in `t` or `z`. -> Task 2.
4. **Mice with `NA` group** (not in meta) are plotted but not scored, no error. -> Task 2.
5. **Panel without any `median_` metric, or `qc_mfi_populations` empty** -> MFI section shows a note, no empty-plot error. -> Task 4 (plot returns `NULL`), Task 5 (template note).

---

### Task 1: QC settings in `facs_prep_save()`

**Files:**
- Modify: `R/facs_prep.R` (`facs_prep_save()` signature, validation, `info$qc`)
- Create: `R/facs_qc.R` (start with `qc_settings_()`)
- Test: `tests/testthat/test-facs_prep.R`, `tests/testthat/test-facs_qc.R`

**Interfaces:**
- Produces: `facs_prep_save(data, wsp_path, assay, panel_dir, condition_keyword = NA_character_, group_col = NULL, batch_col = NULL, expected_mice = NULL, expected_tissues = NULL, qc_populations = NULL, qc_mfi_populations = character(), qc_metrics = NULL)` -- stores `info$qc`, a list with exactly these names: `group_col`, `batch_col` (character scalar or `NULL`), `expected_mice` (character or `NULL`), `expected_tissues`, `qc_populations`, `qc_mfi_populations`, `qc_metrics` (character vectors, resolved -- never `NULL` except as noted below).
- Produces: `qc_settings_(prep) -> list` (internal, same names as `info$qc`). Returns `prep$info$qc` when present; for an old prep without it: `group_col = NULL`, `batch_col = NULL`, `expected_mice = NULL`, and the defaults below.
- Produces: `qc_resolve_(data, qc = list()) -> list` (internal) -- takes a partial settings list (any of the seven names), validates columns / populations against `data` (errors naming every unknown value), fills defaults, returns the full list. Used by `facs_prep_save()`, `qc_settings_()` and tests.
- Defaults when resolving: `expected_tissues` = sorted unique `data$tissue`; `qc_populations` = all populations in factor-level order (or first-appearance order if not a factor); `qc_metrics` = all metrics starting `median_`; `qc_mfi_populations` = `character()`.

- [ ] **Step 1: Write the failing tests**

In `test-facs_prep.R`:
- `facs_prep_save stores resolved qc settings in info` -- save `fake_prep_data()` with a `grp` and `day` column added, `group_col = "grp"`, `batch_col = "day"`, `expected_mice = c("25-7-1","25-7-2")`; `readRDS()` the result; `expect_equal(info$qc$group_col, "grp")`, `expect_equal(info$qc$qc_populations, c("live_CD45+", "CD4+"))`, `expect_equal(info$qc$expected_tissues, "spleen")`, `expect_equal(info$qc$qc_mfi_populations, character())`.
- `facs_prep_save errors on unknown qc columns or populations` -- `group_col = "nope"` -> `expect_error(..., "nope")`; `qc_populations = "B cells"` -> `expect_error(..., "B cells")`; `qc_mfi_populations = "Treg"` (absent) -> `expect_error(..., "Treg")`.
- `facs_prep_save keeps extra meta columns in data` -- `grp` column survives in `readRDS(path)$data`.

In `test-facs_qc.R`:
- `qc_settings_ fills defaults for a prep saved without qc settings` -- build `prep <- list(data = fake_prep_data(files), info = list())`; `s <- qc_settings_(prep)`; `expect_null(s$group_col)`; `expect_equal(s$qc_populations, c("live_CD45+", "CD4+"))`; `expect_equal(s$qc_metrics, character())`.

- [ ] **Step 2: Run to verify they fail**

Run: `devtools::test(filter = "facs_prep|facs_qc")` -- Expected: FAIL (unused arguments / `qc_settings_` not found).

- [ ] **Step 3: Implement** the new arguments in `facs_prep_save()`, delegating validation and defaults to `qc_resolve_()`. Document the new `@param`s (ASCII roxygen) and extend `@returns` with `info$qc`.

- [ ] **Step 4: Run tests, then the full suite** -- `devtools::test()`; Expected: 0 failures.

- [ ] **Step 5: Document and commit**

```bash
Rscript -e 'devtools::document()'
git add R/facs_prep.R R/facs_qc.R man/ tests/testthat/test-facs_prep.R tests/testthat/test-facs_qc.R
git commit -m "feat(facs): store QC settings in the prep result"
```

---

### Task 2: Scoring -- `qc_score_()`

**Files:**
- Modify: `R/facs_qc.R`
- Modify: `tests/testthat/helper-prep.R` (add `fake_qc_data()`)
- Test: `tests/testthat/test-facs_qc.R`

**Interfaces:**
- Consumes: `qc_settings_()` list (Task 1).
- Produces: `qc_score_(data, settings, excluded_files = character(), min_n = 4, min_events = 50, z_cut = 3.5, mad_floor = c(logit = 0.2, log10 = 0.1)) -> tibble`: the rows of `data` that belong to a check, plus columns `check` (`"hierarchy"` = `fraction_of_parent` in `qc_populations`; `"counts"` = `count` in `qc_populations`; `"mfi"` = `qc_metrics` in `qc_mfi_populations`), `n_events` (the population's `count` in that file), `t` (transformed value, `NA` if unscorable), `excluded` (lgl), `z` (dbl, `NA` when not scored), `flag` (lgl, never `NA`). Requires `settings$group_col`; if it is `NULL`, returns the rows with `z = NA`, `flag = FALSE`.
- Produces (helper): `fake_qc_data(n_per_group = 8, outlier = "m-3")` -- long tibble with `file_name`, `population_full_path`, `population`, `metric`, `value`, `mouse_ID`, `tissue` (`"spleen"`), `condition` (`"none"`), `grp` (`"A"`/`"B"`), `day` (two Dates); populations `live` (count ~1e5) and `live/T` (count ~1e4) with `count` and `fraction_of_parent`, a rare `live/T/rare` (count 10), and `median_CD25` for `live/T`; small deterministic scatter; mouse `outlier` has `live/T` fraction far off.

Transform and scoring (algorithm the signature doesn't fix):

```r
# t per check
hierarchy: p <- (count + 0.5) / (count / value + 1); t <- qlogis(p)   # NA if value == 0
counts:    t <- log10(value + 1)
mfi:       t <- ifelse(value > 0, log10(value), NA)
# then t <- NA where n_events < min_events or group is NA
# per group x tissue x condition x population_full_path x metric,
# reference = rows with !is.na(t) & !excluded:
#   if n(reference) < min_n -> z = NA
#   else z = (t - median(ref)) / max(1.4826 * mad_raw(ref), floor[scale])
```

- [ ] **Step 1: Write the failing tests** (all on `fake_qc_data()` with `settings <- list(group_col = "grp", qc_populations = c("live","T","rare"), qc_metrics = "median_CD25", qc_mfi_populations = "T", ...)` built via `qc_resolve_()`):
- `qc_score_ flags the planted outlier only` -- `flagged <- unique(r$mouse_ID[r$flag])`; `expect_equal(flagged, "m-3")`.
- `qc_score_ does not score populations under 50 events` -- `expect_true(all(is.na(r$z[r$population == "rare"])))`, `expect_false(any(r$flag[r$population == "rare"]))`.
- `qc_score_ gives finite z when a box has zero spread` -- set every `live` count to 1e5; `expect_true(all(is.finite(r$z[r$population == "live" & r$metric == "count"])))`.
- `qc_score_ keeps fraction 0 / 1 and count 0 finite` -- one row `fraction_of_parent = 1`, one `count = 0` with fraction 0; `expect_false(any(is.infinite(r$t) | is.nan(r$t)))`, `expect_false(any(is.nan(r$z)))`.
- `qc_score_ skips boxes under 4 samples` -- `fake_qc_data(n_per_group = 3)`; `expect_true(all(is.na(r$z)))`.
- `qc_score_ leaves excluded files out of the reference but scores them` -- pass the outlier's file as `excluded_files`; `expect_true(all(r$excluded[r$mouse_ID == "m-3"]))`; its `z` is non-`NA` and `flag` still `TRUE`; the other mice's medians unaffected (`r$flag` count for others = 0).
- `qc_score_ ignores mice with NA group` -- set `grp` to `NA` for one mouse; no error; that mouse's `z` all `NA`.
- `qc_score_ only scores MFI for qc_mfi_populations and positive values` -- with `qc_mfi_populations = character()`, `expect_false("mfi" %in% r$check)`; with `"T"` and one `median_CD25 = -5`, that row has `t = NA`.

- [ ] **Step 2: Run to verify they fail** -- `devtools::test(filter = "facs_qc")`; Expected: FAIL, `qc_score_` not found.

- [ ] **Step 3: Implement `qc_score_()`** in `R/facs_qc.R` per the block above; `n_events` via a join of the `count` rows on `file_name` + `population_full_path`.

- [ ] **Step 4: Run tests, then the full suite** -- Expected: 0 failures.

- [ ] **Step 5: Commit**

```bash
git add R/facs_qc.R tests/testthat/helper-prep.R tests/testthat/test-facs_qc.R
git commit -m "feat(facs): robust-z QC scoring for hierarchy, counts and MFI"
```

---

### Task 3: Tables -- completeness and flagged samples

**Files:**
- Modify: `R/facs_qc.R`
- Test: `tests/testthat/test-facs_qc.R`

**Interfaces:**
- Consumes: `qc_score_()` output (Task 2), `qc_settings_()` (Task 1).
- Produces: `qc_completeness_(data, settings) -> list(n_expected = int, n_present = int, missing = tibble(mouse_ID, tissue, condition))`. Unit = `expected_mice x expected_tissues x` conditions present in `data`. Returns `NULL` when `settings$expected_mice` is `NULL`.
- Produces: `qc_flag_table_(scored, settings) -> tibble` with columns `mouse_ID`, `score` (int: distinct `tissue x condition x check` flagged for that mouse), `tissue`, `condition`, `group`, `batch` (character; `batch` formatted as a date label when the column is Date / POSIXct, `NA` when `batch_col` is `NULL`), `population`, `metric`, `value`, `z`, `excluded`. One row per flagged value; sorted by `score` desc, then `mouse_ID`, then `abs(z)` desc. Zero rows when nothing is flagged.

- [ ] **Step 1: Write the failing tests**
- `qc_completeness_ lists missing mouse x tissue` -- `expected_mice = c("m-1","m-99")` with data only for `m-1`; `expect_equal(res$missing$mouse_ID, "m-99")`; `expect_equal(res$n_expected, 2L)`; `expect_equal(res$n_present, 1L)`.
- `qc_completeness_ counts a missing condition well as missing` -- data with conditions `unstim`/`MPO`, drop `m-1`'s `MPO` rows; `missing` has one row `m-1`, `MPO`.
- `qc_completeness_ returns NULL without expected_mice`.
- `qc_flag_table_ counts a cascading gate once per sample and check` -- fake `scored` tibble: mouse `m-1` flagged in 3 populations x 2 checks in spleen, `m-2` flagged in hierarchy in spleen and lung; `expect_equal(tbl$score[tbl$mouse_ID == "m-1"][1], 2L)`; `expect_equal(tbl$score[tbl$mouse_ID == "m-2"][1], 2L)`; `nrow(tbl) == 8`.
- `qc_flag_table_ formats a POSIXct batch as a date` -- `expect_equal(unique(tbl$batch), "2026-06-15")`.

- [ ] **Step 2: Run to verify they fail** -- Expected: FAIL, functions not found.

- [ ] **Step 3: Implement both functions** in `R/facs_qc.R`.

- [ ] **Step 4: Run tests, then the full suite** -- Expected: 0 failures.

- [ ] **Step 5: Commit**

```bash
git add R/facs_qc.R tests/testthat/test-facs_qc.R
git commit -m "feat(facs): QC completeness and ranked flag tables"
```

---

### Task 4: Plots -- individuals and batch (layout C)

**Files:**
- Modify: `R/facs_qc.R`, `DESCRIPTION` (`Suggests: ggrepel`), CLAUDE.md (Suggests line)
- Test: `tests/testthat/test-facs_qc.R`

**Interfaces:**
- Consumes: `qc_score_()` output, `qc_settings_()`.
- Produces: `qc_plot_individuals_(scored, settings, check, condition = NULL, population = NULL) -> ggplot | NULL` -- x = group, boxplot (`outlier.shape = NA`) + jittered points (fixed seed), labels = `mouse_ID` on `flag` rows only; facets `population ~ tissue` (for `check = "mfi"`: rows = `metric`, one `population` -- the argument), `scales = "free_y"`; y percent for hierarchy, log10 for counts / MFI. Returns `NULL` if no rows.
- Produces: `qc_plot_batch_(scored, settings, check, tissue, condition = NULL, population = NULL) -> ggplot | NULL` -- one tissue; x = `.batch`, a column the function adds to the plot data (factor, date label for Date / POSIXct), fill/colour = batch, facets `population ~ group` (MFI: `metric ~ group`), no labels. Returns `NULL` if no rows or `batch_col` is `NULL`.
- Produces: `qc_fig_height_(n_rows) -> numeric` -- inches, `max(3, 1.6 * n_rows + 1)`.
- Produces: `qc_has_ggrepel_() -> logical` -- `requireNamespace("ggrepel", quietly = TRUE)`; the only place that check lives, so tests can mock it.

- [ ] **Step 1: Write the failing tests** (on `qc_score_(fake_qc_data(), settings)`):
- `qc_plot_individuals_ labels only flagged mice` -- `p <- qc_plot_individuals_(r, s, "hierarchy")`; `expect_s3_class(p, "ggplot")`; the text/label layer's data has `mouse_ID` values `== "m-3"` only (find the layer whose `geom` inherits `GeomText` or `GeomTextRepel`).
- `qc_plot_batch_ facets by group and uses batch on x` -- `p <- qc_plot_batch_(r, s, "counts", tissue = "spleen")`; `expect_s3_class(p$facet, "FacetGrid")`; `x <- p$data$.batch`; `expect_true(is.factor(x))`; `expect_equal(nlevels(x), 2L)`; `expect_match(levels(x), "^\\d{4}-\\d{2}-\\d{2}$")`.
- `qc plots return NULL without data` -- `qc_plot_individuals_(r, s, "mfi")` with `qc_mfi_populations = character()` -> `expect_null()`; `qc_plot_batch_()` with `batch_col = NULL` -> `expect_null()`.
- `qc plots work without ggrepel` -- `local_mocked_bindings(qc_has_ggrepel_ = function() FALSE)`; `expect_no_error(ggplot2::ggplot_build(qc_plot_individuals_(r, s, "hierarchy")))`.

- [ ] **Step 2: Run to verify they fail** -- Expected: FAIL, functions not found.

- [ ] **Step 3: Implement** the three plot functions plus `qc_has_ggrepel_()`; colours are ggplot2 defaults. Add `ggrepel` to `Suggests` in DESCRIPTION and to the Suggests line in CLAUDE.md.

- [ ] **Step 4: Run tests, then the full suite** -- Expected: 0 failures.

- [ ] **Step 5: Commit**

```bash
git add R/facs_qc.R DESCRIPTION CLAUDE.md tests/testthat/test-facs_qc.R
git commit -m "feat(facs): QC individual and batch plots"
```

---

### Task 5: Template, prep script, export check

**Files:**
- Modify: `inst/rmarkdown/templates/facs-qc-template/skeleton/skeleton.Rmd`
- Modify: `inst/scripts/facs_prep_template.R`
- Modify: `R/report_knit.R` (roxygen of `report_knit_facs_qc()`: drop "placeholder sections", describe the checks)
- Test: `tests/testthat/test-report_knit.R`, `tests/testthat/test-db_export.R`

**Interfaces:**
- Consumes: `qc_settings_`, `qc_score_`, `qc_completeness_`, `qc_flag_table_`, `qc_plot_individuals_`, `qc_plot_batch_`, `qc_fig_height_` (Tasks 1-4), called as `JanisHelpers:::fn()` from the template.

Template changes:
- Header table: add rows "Group / batch column" and "QC populations / MFI populations" (`"all (n)"` when equal to every population, else the list).
- Replace `## QC checks` placeholders with sections `### Completeness`, `### Gating hierarchy`, `### Raw counts`, `### MFI`, `### Flagged samples`. Each plot section: individuals plot, then `#### Batch -- <tissue>` per tissue. Multi-well panels: loop conditions, one sub-heading per condition (single condition `"none"` -> no sub-heading). Generate plot chunks dynamically (`knitr::knit_child()` or a loop with `results = "asis"` and per-plot `fig.height` from `qc_fig_height_()`).
- Notes instead of errors: no `group_col` -> "Group column not set -- rerun the prep script with `group_col`." (individual scoring and flag table skipped); no `batch_col` -> batch plots skipped with a note; MFI `NULL` -> "No MFI populations selected (`qc_mfi_populations`)."; completeness `NULL` -> "Expected mice not set."; no flags -> "No samples flagged."
- Flag table: `excluded` rendered as "excluded" / "".

Prep script template: add settings `group_col <- "group"`, `batch_col <- "death_date"`, `qc_populations <- NULL`, `qc_mfi_populations <- character()`, `qc_metrics <- NULL`; join `meta$meta` columns with `meta_annotate(data, dplyr::select(meta$meta, mouse_ID, dplyr::all_of(c(group_col, batch_col))), by = "mouse_ID")`; compute `expected_mice <- meta$meta$mouse_ID[!is.na(meta$meta$mouse_ID) & is.na(meta$meta$early)]` (comment: dropout column is experiment-specific); pass all to `facs_prep_save()`.

- [ ] **Step 1: Write the failing tests**

In `test-report_knit.R`:
- `report_knit_facs_qc renders QC checks with flags` -- `skip_if_not(rmarkdown::pandoc_available())`; save `fake_qc_data()` via `facs_prep_save(..., group_col = "grp", batch_col = "day", expected_mice = c(unique(mouse_ID), "m-99"), qc_mfi_populations = "T")` (adapt `file_name` / wsp to `local_prep_dirs()`); render; `expect_match(html, "Flagged samples")`, `expect_match(html, "m-3")`, `expect_match(html, "m-99")` (completeness), `expect_match(html, "Batch")`.
- `report_knit_facs_qc renders an old prep without qc settings` -- `saveRDS()` a prep whose `info` lacks `qc`; render succeeds; `expect_match(html, "rerun the prep script")`.
- `report_knit_facs_qc notes when no MFI populations are selected` -- `expect_match(html, "No MFI populations selected")`.

In `test-db_export.R`:
- `db_export_facs ignores extra meta columns in prep data` -- like the existing happy-path export test but with `grp` / `day` columns in the data; export succeeds and `facs_measurements` row count equals the run without them.

- [ ] **Step 2: Run to verify they fail** -- `devtools::test(filter = "report_knit|db_export")`; Expected: the three report tests FAIL (the export test may already pass -- keep it as a guard).

- [ ] **Step 3: Implement** the template and prep-script changes above; update `report_knit_facs_qc()` roxygen. Confirm `test-facs_prep.R`'s "prep script template ships and parses" still passes.

- [ ] **Step 4: Run tests, then the full suite** -- Expected: 0 failures.

- [ ] **Step 5: Document and commit**

```bash
Rscript -e 'devtools::document()'
git add inst/ R/report_knit.R man/ tests/testthat/test-report_knit.R tests/testthat/test-db_export.R
git commit -m "feat(report): QC checks section in the facs QC report"
```

---

### Task 6: Tune on real data, docs, check

**Files:**
- Modify: `R/facs_qc.R` (only if the floor defaults change)
- Modify: `docs/superpowers/specs/2026-11-13-facs-qc-checks-design.md` (final floor values), `docs/superpowers/plans/2026-11-12-facs-pipeline-overview.md` (tick Stage 1 QC checks), CLAUDE.md (data-structures FACS prep bullet: `info$qc`, QC checks)

- [ ] **Step 1: Score the real 26-1 data read-only.** In a script under `$TMPDIR` (never in the OneDrive folder), for `overview-26-1_20260806_prep.rds` and `Treg_26-1_20260828_prep.rds` (paths in the spec): join `group` and `death_date` from `meta_26-1.xlsx` sheet `meta`, build settings with `qc_resolve_()` (`qc_mfi_populations = "Treg_FoxP3+_CD25hi"` for Treg), run `qc_score_()` with floors `c(0.2, 0.1)`, `c(0.3, 0.15)`, and `qc_flag_table_()`. Expected: Treg top mouse 26-1-2 (RORgt gates); no flags in populations under 50 events; total flagged values below the prototype's 98.

- [ ] **Step 2: Render one real report into `$TMPDIR`** -- copy the Treg prep `.rds` with added `info$qc` into a temp `processed/` folder and call `report_knit_facs_qc()`; open the HTML path for the user to look at. Expected: renders without error, all five QC sections present.

- [ ] **Step 3: Record the chosen floor** in the spec's "Left to the builder" section (replace "Starting values" with the final values and one line of evidence); change the `mad_floor` default only if it moved.

- [ ] **Step 4: Update CLAUDE.md and the overview plan** -- FACS prep bullet: `facs_prep_save()` also stores `info$qc` (group / batch column, expected mice, QC selections) and the QC report runs completeness, hierarchy, counts, MFI checks with robust-z flags. Tick "QC checks in the report" in Stage 1.

- [ ] **Step 5: Full check** -- `Rscript -e 'devtools::check()'`; Expected: 0 errors, 0 warnings, only known NOTEs.

- [ ] **Step 6: Commit**

```bash
git add R/ docs/ CLAUDE.md
git commit -m "docs: QC checks tuned on 26-1, CLAUDE.md and pipeline overview"
```
