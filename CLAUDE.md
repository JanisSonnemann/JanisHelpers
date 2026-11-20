# CLAUDE.md

Personal R utility package for biomedical research workflows at Charité Berlin: FlowJo/FCS import and analysis, experiment metadata, multiplex ELISA import, statistical tables, plots, RMarkdown knitting helpers, and a DuckDB lab database. Installed from GitHub via `devtools::install_github("JanisSonnemann/JanisHelpers", auth_token = gh::gh_token())`.

---

## How to work in this package

1. Read this file fully before writing any code
2. Check existing functions in the relevant domain file before writing new ones
3. Run `devtools::load_all()` before testing interactively
4. Run `devtools::check()` before every commit — target 0 errors, 0 warnings
5. Follow naming convention strictly — new domains require updating this file first

---

## Development commands

```r
devtools::document()   # regenerate man/*.Rd after any @param/@returns/@export change
devtools::check()      # full R CMD check — target: 0 errors, 0 warnings
devtools::load_all()   # fast iteration without reinstalling
devtools::install()    # full install
```

---

## File and function naming

| Layer | Pattern | Examples |
|---|---|---|
| File | `domain_verb.R` | `facs_import.R`, `analysis_stats.R`, `report_knit.R` |
| Function | `domain_verb()` or `domain_verb_qualifier()` | `facs_import_wsp()`, `report_knit_dated()`, `report_knit_wide()` |

**Domains** (use exactly these prefixes):
- `facs_` — FlowJo / flow cytometry
- `analysis_` — statistical summaries and tests
- `plot_` — ggplot2 visualization of FACS/analysis data
- `report_` — RMarkdown rendering
- `meta_` — experiment/subject metadata import and annotation
- `elisa_` — ELISA: multiplex bead (Luminex/Bio-Plex) results import; curated sheets import, QC and export
- `histo_` — histology scoring import (curated sheets), QC and export
- `wrangle_` — general data wrangling [stub — no functions yet]
- `db_` — database access (DuckDB: connect/write/query the FACS/ELISA/histology schema)

New functions must follow the `domain_verb` pattern. Never add a function that doesn't belong to a declared domain.

---

## Code style

- **Pipe**: use `|>` (base pipe) everywhere. `%>%` is only tolerated inside `gtsummary` chains where method dispatch requires it; annotate with a comment if used.
- **Namespacing**: always call external functions as `pkg::fn()` inside exported function bodies. No bare `map()`, `filter()`, `select()`, `set_names()` etc. — these cause R CMD check NOTEs.
  - `dplyr::filter()`, `dplyr::select()`, `dplyr::mutate()`, `dplyr::left_join()`, etc.
  - `purrr::map()`, `purrr::map_chr()`
  - `tidyr::pivot_longer()`, `tidyr::pivot_wider()`, `tidyr::unnest()`
  - `stringr::str_extract()`
  - `glue::glue()`
  - `gt::tab_style()`, `gt::tab_header()`, etc.
  - `gtsummary::tbl_summary()`, etc.
  - **Known debt**: `analysis_stats.R` currently has bare calls — fix when touching that file.
- **Column references in tidy eval**: bare column names inside `dplyr` verbs are fine. For `.data$col` pronoun use, prefix with `.data` only when the column name is stored in a variable.
- **Deselection**: `dplyr::select(!col)`, not `-col`.
- **`invisible()`**: functions that return data frames/tibbles primarily for side-effects (e.g., import + message) should return `invisible(result)`.
- **Non-ASCII**: escape all non-ASCII characters in R source with `\uXXXX` (e.g., em-dash = `—`). Only inside R string literals -- in roxygen `#'` text write plain ASCII (`--`), since `\uXXXX` there breaks the `.Rd` build.

---

## Data structures

Roxygen docs (`?fn`) are the reference for arguments and return shapes; this is the map.

- **`facs_read_wsp()`** (current) -- named list `data` / `meta` / `panel`. `data` is long, one row per `file_name x population_full_path x metric`; `population` is a factor ordered by gating-hierarchy traversal, not alphabetically; requested keywords (e.g. `mouse_ID`, `tissue`, `panel`) appended as columns. `facs_import_wsp()` is the older CamelCase variant (`FileName`, `PopulationFullPath`, ...).
- **`facs_calc_*()`** -- derived metrics on the long tibble (`pct_of`, `count_per_g`, `log2fc`, `diff`, cluster frequencies/medians). `facs_calc_count_per_g()` takes `meta` from `meta_clean()` and joins on `mouse_ID` + `tissue`. `facs_calc_log2fc()` / `facs_calc_diff()` take `mouse_col` (default `"mouse_ID"`); on `db_query_facs()` output use `mouse_col = "mouse_id", restim_col = "condition"`. `facs_order_populations()` turns collected `db_query_facs()` output's `population` back into the prep's display-order factor (`population_order`; panels alphabetical by `wsp_root`, a shared label keeps its first position, `NA` orders last with a message).
- **`facs_read_fcs_gated()` / `facs_cluster_flowsom()` / `facs_reduce_umap()` / `facs_test_cluster_abundance()` / `facs_plot_*()`** -- event-level unsupervised pipeline (gated events -> FlowSOM -> UMAP -> diffcyt).
- **`meta_read()`** -- Excel workbook -> named list of cleaned tibbles (one per sheet). **`meta_clean()`** -- joins `organ_weights` (pivoted long) with `facs_volumes` by `mouse_ID` + base tissue. **`meta_annotate(data, meta, by)`** -- checked left join; errors on missing keys / name collisions, warns on unmatched rows.
- **`meta_summary_table(workbook, dat)` / `meta_mouse_table(dat, workbook)` / `meta_label_weeks(data, meta)`** (`R/meta_tables.R`, `R/meta_weeks.R`, spec: `docs/superpowers/specs/2026-11-20-cohort-tables-weeks-timecourse-design.md`) -- cohort tables for analysis docs (`gt`; data in `tbl[["_data"]]`). Cohort = `dat$cohort` (all mice incl. facility-ID dropouts, from the DB; `workbook = NULL` uses it) or, as before, the meta workbook's `meta` sheet (`meta_read()$meta`, `dat$workbook`); DB names (`mouse_id`, `end_date`) are accepted (spec: `docs/superpowers/specs/2026-11-20-exclusions-in-db-design.md`). Excluded = `dat$cohort$excluded_mouse_wide` (else mouse-wide rules in `dat$exclusions`); analysed = end-point data in >= 1 unit and not excluded (early-ending mice count). Summary: per group + `not assigned` + `total`. Mouse grid: one column per unit -- `n/N` tissues (missing named), `n/N` serial dates + 1 end-point, or a check mark; `status` `ok` / `early` / `excluded` / `not in DB` (+ partial exclusions). `meta_label_weeks()` adds `days_post_bmt`, `week_num` (days / 7 rounded, halves up) and ordered `week`; end-points of one experiment in the same ISO calendar week share the week of their median day.
- **`analysis_summary_table()` / `analysis_posthoc_tables()`** -- wide tibble with `group` (+ optional `tissue`) -> `gt` table(s); non-significant post-hoc variables return a character string.
- **`plot_bubble_fc()`** -- long tibble (one metric, one tissue) -> `ggplot` bubble plot of log2FC vs `control` with significance stars; `x_col` switches the x-axis to another column. Colour limits symmetric around 0, legend ends say `higher than <control>` (red) / `lower than <control>` (blue).
- **Excluded mice / labels on per-mouse plots** (`plot_box()`, `plot_population()`, `plot_paired()`, `plot_timecourse()`; spec: `docs/superpowers/specs/2026-11-20-plot-excluded-labels-design.md`) -- rows with `excluded == TRUE` (data loaded with `keep_excluded = TRUE`) are dropped with a message by default (`plot_timecourse_summary()` too); `show_excluded = TRUE` draws them as red crosses with a red mouse-ID label (timecourse: red dashed line, label at the last time point), never in boxes / medians / p-values. `label = TRUE` labels every point (`id` column in `plot_box()`/`plot_population()`, `pair` in `plot_paired()`, `subject` in `plot_timecourse()`); internals `plot_split_excluded_()` / `plot_add_labels_()` / `plot_label_layer_()` in `R/plot_box.R`. Point and label layers share rows and jitter seed, with `group` pinned to the fill so the `.excluded` aesthetics don't change the dodge.
- **`plot_timecourse()` / `plot_timecourse_summary()`** (`R/plot_timecourse.R`) -- long data (ELISA after `meta_label_weeks()`) over `x = week_num` (numeric = real time, labelled `week <n>`; a factor = categories); `NA` values dropped with a message, replicates averaged per mouse x time (error if several analyte / tissue / condition / population / metric values would be averaged -- filter or facet). Full: each mouse + thick group median; summary: median + IQR error bars, groups dodged, n per group at the bottom. Time points with fewer than `min_n` (2) mice are bridged by the median line (summary: hollow point, no bar). Optional `lloq` dashed line, `facet`, `scale = "log"`, `p_val` (Wilcoxon / Kruskal per time point, BH per facet); colours as `plot_box()` (`add_fill_scale_(aesthetics = c("colour", "fill"))`).
- **`plot_box()` / `plot_population()` / `plot_paired()`** (`R/plot_box.R`) -- group-comparison plots for the analysis docs; column arguments are bare names (tidy eval), `fill` defaults to `group` (as from `db_query_subjects()`), `scale = "linear"/"percent"/"log"`, optional named `colors` (default `getOption("JanisHelpers.group_colors")`, so a doc sets group colours once; explicit `NULL` = ggplot defaults). `plot_population()` filters populations / metric(s) / tissues, adds `pseudocount` for log, facets by population, `p_val = TRUE` adds `ggpubr::stat_compare_means()`. `plot_paired()` keeps the given condition `levels` (in order) and joins each mouse across conditions: one cluster per group (numeric x), axis labels = condition names, group name above each cluster, conditions as shades of the group colour (first level lightest, `alpha` legend).
- **`elisa_read_results()`** -- `Results_<cytokine>.xlsx` (`*Unknowns` sheet) -> long tibble `cytokine`, `sample_id`, `replicate`, `value`, `unit`, `result_status`. `sample_id` is not assumed to equal `mouse_ID`. Outside the prep -> QC -> export pipeline below.
- **`elisa_read_curated()` / `histo_read_curated()`** -- one curated `.xlsx` per assay (sheet `values`: `mouse_ID`, optional `tissue`, optional `date`, one column per analyte / metric; ELISA optional sheet `range` with `LLOQ`/`ULOQ`/`unit` rows) -> long tibble with `mouse_ID`, `tissue`, `collected_at` (`NA` = end-point), `analyte` (ELISA, + `replicate`, `value`, `result_status` `"ok"`/`"OOR<"`/`"OOR>"`, `unit`, `lloq`, `uloq`) or `metric` (histo). Repeated mouse x tissue x date rows are ELISA replicates. Histology: a `comment`/`comments` column (any case) is free text -> long column `comment` (only when present; QC report "Comments" section, never exported). `histo_calc_ratio()` appends ratio metrics -- one (`name`) or several via a named `numerator` vector (names = new metrics), copying other columns from the numerator rows. The histology QC report has no batch / time outlier plots. ELISA / histology outlier and CV plots draw excluded samples as crosses with red labels (`excluded_keys`); outlier / high-CV labels are black. Blank templates (with a README sheet): `inst/templates/{elisa,histo}_curated_template.xlsx`, rebuilt by `data-raw/make_curated_templates.R`.
- **`elisa_reshape_wide()`** -- optional step before ELISA prep: wide sheet (one row per mouse [x analyte x tissue], optional `analyte`/`tissue` columns, one column per date, `final` = end-point) -> curated `<stem>_long.xlsx` next to it (`_wide` suffix replaced); other columns dropped; empty cells give no row; rewrites only when the content changes (keeps the prep stale check); returns the path for `elisa_prep_run(path = )`.
- **ELISA / histology prep -> QC -> export loop** (spec: `docs/superpowers/specs/2026-11-18-elisa-histo-prep-qc-export-design.md`) -- same shape as the FACS loop, one script per curated file from `inst/scripts/elisa_prep_template.R` / `histo_prep_template.R` with a `final` switch; the one call is `elisa_prep_run()` / `histo_prep_run(final, path, assay, tissue, ...)`. Import unit (root) = file name without extension; `*_prep_save()` writes `list(data, info)` to `<dir>/processed/<root>_prep.rds`, `report_knit_elisa_qc()` / `report_knit_histo_qc()` render to `<dir>/qc/<root>_qc.html`, `*_prep_load()` refuses if stale or the report is missing. `final = TRUE` calls `db_export_elisa()` / `db_export_histo()`: writes every row (excluded ones are flagged on read, not dropped), replaces the import unit's earlier rows in one transaction, writes `import_log` (not `import_exclusions`, kept as history); `*_prep_run(final = TRUE, meta_path)` also records the `exclusions` sheet via `db_write_exclusions()`. The sheet may carry `collected_at` (empty = all dates; empty sample dates resolve to the mouse's `end_date` before matching, spec: `docs/superpowers/specs/2026-11-20-exclusions-in-db-design.md`). `*_prep_save(end_dates)` (filled from `mice$death_date` by `*_prep_run()`) places end-point samples in the QC report (completeness: an end-point-only date expects only the mice ending that day; a serial date expects mice still alive); unknown exclusion panels are a message from `db_write_exclusions()`.
- **`report_knit_*()`** -- render an `.Rmd` into dated output; return the output path. No RStudio dependency.
- **FACS prep / QC loop** (spec: `docs/superpowers/specs/2026-11-12-facs-prep-qc-export-design.md`) -- per-panel prep script from `inst/scripts/facs_prep_template.R` with a `final` switch; its one call is `facs_prep_run(final, wsp_path, panel_dir, assay, prepare, mice, exclusions, con, ...)`, which chains the steps below (QC settings pass through `...` to `facs_prep_save()`; `con` is opened and closed by the panel script, which registers experiment + mice with `db_write_meta(con, meta, experiment_code, ...)` from the standard `meta` sheet; one self-contained script per panel, no shared project script). `facs_prep_save()` writes `list(data, info)` (all samples, plus a `condition` column) to `<panel folder>/processed/<wsp>_prep.rds`; `report_knit_facs_qc()` renders `facs-qc-template` to `<panel folder>/qc/<wsp>_qc.html`; `facs_prep_load()` reloads it for the export (adds `paths`; refuses if stale or the report is missing). `facs_prep_save()` also stores `info$qc` (group / batch column, expected mice / tissues, QC population and MFI selections, `min_events`; older preps get the default via `qc_settings_()`); it also records `info$dropped_files` (files without `mouse_ID`/`tissue`). The QC report runs completeness, empty-value / missing-stat / keyword checks (metric x tissue grid of `NA` values -- empties in 0-event populations count as expected; gates or stats present in only some files of a tissue; dropped files; empty group/batch), gating-hierarchy, raw-count and MFI checks (outlier labels = points outside the boxplot whiskers, computed on the plot axis from the same points, per-tissue batch plots, ranked flag table; samples on the exclusions sheet are drawn as crosses on every one-point-per-sample plot, with a red `mouse_ID` label where plots label (outlier labels black; an excluded outlier gets only the red one), stats unchanged, spec: `docs/superpowers/specs/2026-11-20-qc-excluded-markers-design.md`) via internal helpers in `R/facs_qc.R` (spec: `docs/superpowers/specs/2026-11-13-facs-qc-checks-design.md`). The `prepare` function selects populations with `facs_select_populations(data, populations, labels)` (filter + display-order factor; errors on names not in the workspace, suggesting the closest). `facs_match_exclusions()` resolves the meta workbook's `exclusions` sheet (`mouse_ID` + `reason` required; empty/absent `tissue`/`panel`/`condition`/`file_name` = all) -> `applicable` / `excluded_files` / `unmatched`; `facs_prep_run(final = TRUE, meta_path)` also writes the rule set via `db_write_exclusions()` (experiment code inferred; unmatched rules are messages).
- **`db_load_experiment(con, experiment, names, group_levels, subject_cols = "group", keep_excluded = FALSE)`** (`R/db_load.R`, spec: `docs/superpowers/specs/2026-11-19-analysis-db-loader-templates-design.md`) -- one experiment as one named list for analysis docs: `meta` (`db_query_subjects()` rows; default drops facility-ID dropouts and mouse-wide-excluded mice), `cohort` (all mice, flags kept), `exclusions` (current rule set; message if none), then one collected tibble per import unit (default drops excluded rows; `keep_excluded = TRUE` keeps all with `excluded` / `exclusion_reason`; spec: `docs/superpowers/specs/2026-11-20-exclusions-in-db-design.md`) (FACS `wsp_root`, then ELISA, then histology `import_unit`, each alphabetical, case-insensitive) with `subject_cols` joined via `meta_annotate()` and FACS populations ordered per panel. Names = import unit minus trailing `_long` and `-`/`_` + experiment code; `names = c(new = "import_unit")` overrides; `meta`/`cohort`/`exclusions` reserved. Read-only connection is the intended use; `db_query_facs()` returns all-`NA` `population_order` on files without the column.
- **`db_connect()` / `db_write_*()` / `db_query_*()`** -- DuckDB lab database (schema: `docs/superpowers/specs/2026-07-17-facs-db-schema-design.md`). `db_connect()` creates the schema idempotently. `db_write_*()` are idempotent; dimension writers error listing every unresolved key. Samples are keyed by `mouse_id` x `tissue` x `collected_at` (`DATE NOT NULL`); a missing `collected_at` means the subject's `end_date`. `db_write_elisa()`/`db_write_histo()` no longer exist: `db_export_elisa()`/`db_export_histo()` (like `db_export_facs()`) are the import paths, writing `elisa_measurements` / `histo_measurements` (keyed by `import_unit`, replaced on reimport; ELISA stores `analyte`, `replicate`, `result_status`, `lloq`/`uloq`, never substituted values). `db_write_facs()` tags stains with `wsp_root` (one root = one experiment x one panel) and refuses a second root or workspace version. `db_export_facs(con, prep)` is the FACS import path: writes all files (stains carry `file_name`), asks before overwriting a stored root, and in one transaction replaces the root and writes `import_log` (`qc_passed = TRUE`). Wells differing in one variable (e.g. restimulation) are separate stains via `condition` + `condition_keyword` (default `"none"`); ratios between wells are computed in analysis, not stored. No volumes/weights are stored -- they are calculation inputs. `source_file` is a path relative to the project's `data/` folder. Snapshots: `db_snapshot(con, root)` -> `<root>/<YYYY-MM-DD>/` (same-day rerun replaces it; con stays open), via DuckDB `EXPORT DATABASE`, which rejects self-referencing foreign keys -- don't add any. `db_query_*()` return lazy `dbplyr` tables keyed by `experiment_code`, `mouse_id`, `tissue`, `collected_at`, `assay_name` (+ `import_unit` / `wsp_root`) -- `dplyr::collect()` to materialize. `db_query_facs/elisa/histo()` add `excluded` (never `NA`) + `exclusion_reason` from views over the rule set. Exclusions (spec: `docs/superpowers/specs/2026-11-20-exclusions-in-db-design.md`): `db_write_exclusions(con, exclusions, experiment_code, source_file, note)` stores the meta workbook's `exclusions` sheet as a versioned per-experiment rule set (`exclusion_sets` / `exclusion_rules`; new set only when changed; every written set is reported in a message, replacing a non-empty set by an empty one warns; a `*_prep_run(final = TRUE)` without `exclusions` leaves the rule set unchanged); `facs_/elisa_/histo_exclusion_flags` views match rules to rows on read. `db_write_meta()` registers every mouse -- dropouts without `mouse_ID` by facility `id` (`mouse_id` = facility ID, `id_source = 'facility'`), plus `early` / `early_reason`; it refuses `end_date` changes and re-keying a facility ID. Migration is additive on `db_connect()`; re-run `db_write_meta()` + each unit's final script to backfill. They carry no other subject meta: `db_query_subjects()` returns one row per mouse (`experiment_code`, `mouse_id`, `group` = stored `subjects.treatment_group` -- not renamed in the DB, `group` is an SQL reserved word -- plus sex, cage, strain, dates, `mdc_id`, `early`, `early_reason`, `id_source`, `excluded_mouse_wide`, `exclusion_reason`); analysis docs query it once (read-only connection), set group levels once and join needed columns per dataset with `meta_annotate(by = c("experiment_code", "mouse_id"))`. `facs_measurements.population_order` (spec: `docs/superpowers/specs/2026-11-19-population-order-plot-helpers-design.md`) stores each population's position in its panel's display order (from the prep's factor, written by the export; added on open to older files, `NULL` until the panel is re-exported).

---

## Dependency philosophy

- **`Imports`**: every package called inside an exported function body must be listed here. All calls must be namespaced `pkg::fn()`.
- **`Suggests`**: packages only needed in examples, vignettes, or tests -- or used optionally behind `requireNamespace()`. Currently `ggrepel` (QC report outlier labels; plain `geom_text` without it), `knitr` (template parse test), `testthat` and `withr` (used for `withr::defer()` cleanup in `db_` domain tests).
- **Non-CRAN packages** (`fcexpr`): list in `Imports` as normal; document the install requirement in README. Do not add `remotes::` calls inside function bodies.

Current `Imports`: `CytoML`, `DBI`, `dbplyr`, `diffcyt`, `dplyr`, `duckdb`, `fcexpr`, `flowCore`, `FlowSOM`, `flowWorkspace`, `ggplot2`, `ggpubr`, `glue`, `gt`, `gtsummary`, `janitor`, `parallel`, `purrr`, `readxl`, `rmarkdown`, `rstatix`, `scales`, `stats`, `stringr`, `SummarizedExperiment`, `tibble`, `tidyr`, `tools`, `utils`, `uwot`, `writexl`, `xml2`, `xfun` (order matches `DESCRIPTION`'s `Imports:` block).

---

## Testing conventions

- `testthat` (edition 3). Test files mirror source: `R/facs_read.R` -> `tests/testthat/test-facs_read.R`.
- Workspace-reading functions use real fixtures in `tests/fixtures/`. Do not mock `fcexpr` calls.
- `db_` tests use a temporary DuckDB file cleaned up with `withr::defer()` (see `helper-db.R`).

---

## RMarkdown templates

`inst/rmarkdown/templates/` — available in RStudio *New File → R Markdown*:

| Template dir | Knit function |
|---|---|
| `auto-knit-report` | `report_knit_dated()` via YAML `knit:` field |
| `experiment-report` | `report_knit_exp()` |
| `wide-html` | `report_knit_wide()` |
| `facs-qc-template` | `report_knit_facs_qc()` |
| `histo-qc-template` | `report_knit_histo_qc()` |
| `elisa-qc-template` | `report_knit_elisa_qc()` |
| `analysis-db` | `report_knit_exp()` -- setup: read-only `db_load_experiment()` -> `dat` (`dat$cohort` feeds the cohort tables), population vectors in `pops`; overview = cohort tables, ELISA = timecourse example |
| `analysis-wsp` | `report_knit_exp()` -- quick route from `.wsp`, same `dat` / `pops` shape (`mouse_id`, `condition`) |
| `protocol-template` | standalone, no custom knit function |

`report_knit_wide()` injects `inst/resources/wide-output.css` via `system.file()`.

---

## Roxygen conventions

- Every exported function needs `@param`, `@returns`, and `@export`.
- `@returns` must have a value — empty tag causes a roxygen warning and blocks `devtools::document()`.
- `@examples` — omit the tag entirely rather than leaving it empty. Use `\dontrun{}` when: (a) the example requires an external file (`.wsp`, `.Rmd`), or (b) the function lives in `analysis_stats.R` (bare-call debt causes examples to fail under `R CMD check`'s clean environment).
- After any roxygen change: run `devtools::document()` before committing.

---

## Known check output

Target: **0 errors, 0 warnings**. Expected NOTEs (do not chase):

- "no visible binding for global variable" / `.data` / `:=` / `%>%` -- tidy-eval column names, `ggplot2::aes()`, `gt` `rows =` predicates. False positives. Column names that collide with functions (`median`, `stars`) also trigger an `importFrom` suggestion -- ignore it.
- `all_pops_fun` multiple definitions in `analysis_posthoc_tables()` -- known, low priority.
- Environment-only: future file timestamps; `.git` hidden file (only when checking from a git worktree); top-level `docs` directory.

A new NOTE that is *not* one of these is a real issue -- fix it.
