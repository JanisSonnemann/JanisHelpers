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
- `elisa_` — multiplex bead ELISA (Luminex/Bio-Plex) import
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
- **`facs_calc_*()`** -- derived metrics on the long tibble (`pct_of`, `count_per_g`, `log2fc`, `diff`, cluster frequencies/medians). `facs_calc_count_per_g()` takes `meta` from `meta_clean()` and joins on `mouse_ID` + `tissue`.
- **`facs_read_fcs_gated()` / `facs_cluster_flowsom()` / `facs_reduce_umap()` / `facs_test_cluster_abundance()` / `facs_plot_*()`** -- event-level unsupervised pipeline (gated events -> FlowSOM -> UMAP -> diffcyt).
- **`meta_read()`** -- Excel workbook -> named list of cleaned tibbles (one per sheet). **`meta_clean()`** -- joins `organ_weights` (pivoted long) with `facs_volumes` by `mouse_ID` + base tissue. **`meta_annotate(data, meta, by)`** -- checked left join; errors on missing keys / name collisions, warns on unmatched rows.
- **`analysis_summary_table()` / `analysis_posthoc_tables()`** -- wide tibble with `group` (+ optional `tissue`) -> `gt` table(s); non-significant post-hoc variables return a character string.
- **`plot_bubble_fc()`** -- long tibble (one metric, one tissue) -> `ggplot` bubble plot of log2FC vs `control` with significance stars; `x_col` switches the x-axis to another column.
- **`elisa_read_results()`** -- `Results_<cytokine>.xlsx` (`*Unknowns` sheet) -> long tibble `cytokine`, `sample_id`, `replicate`, `value`, `unit`, `result_status`. `sample_id` is not assumed to equal `mouse_ID`.
- **`report_knit_*()`** -- render an `.Rmd` into dated output; return the output path. No RStudio dependency.
- **FACS prep / QC loop** (spec: `docs/superpowers/specs/2026-11-12-facs-prep-qc-export-design.md`) -- per-panel prep script from `inst/scripts/facs_prep_template.R` with a `final` switch; its one call is `facs_prep_run(final, wsp_path, panel_dir, assay, prepare, mice, exclusions, con, ...)`, which chains the steps below (QC settings pass through `...` to `facs_prep_save()`; `con` is opened and closed by the project script, which also registers experiment + subjects). `facs_prep_save()` writes `list(data, info)` (all samples, plus a `condition` column) to `<panel folder>/processed/<wsp>_prep.rds`; `report_knit_facs_qc()` renders `facs-qc-template` to `<panel folder>/qc/<wsp>_qc.html`; `facs_prep_load()` reloads it for the export (adds `paths`; refuses if stale or the report is missing). `facs_prep_save()` also stores `info$qc` (group / batch column, expected mice / tissues, QC population and MFI selections); the QC report runs completeness, gating-hierarchy, raw-count and MFI checks (robust-z outlier labels, per-tissue batch plots, ranked flag table) via internal helpers in `R/facs_qc.R` (spec: `docs/superpowers/specs/2026-11-13-facs-qc-checks-design.md`). `facs_match_exclusions()` resolves the meta workbook's `exclusions` sheet (`mouse_ID` + `reason` required; empty/absent `tissue`/`panel`/`condition`/`file_name` = all) -> `applicable` / `excluded_files` / `unmatched`.
- **`db_connect()` / `db_write_*()` / `db_query_*()`** -- DuckDB lab database (schema: `docs/superpowers/specs/2026-07-17-facs-db-schema-design.md`). `db_connect()` creates the schema idempotently. `db_write_*()` are idempotent; dimension writers error listing every unresolved key. `db_write_elisa()`/`db_write_histo()` take whole tables with `mouse_id`/`tissue` columns and write one `import_log` row per call (FACS logging comes from `db_export_facs()`). `db_write_facs()` tags stains with `wsp_root` (one root = one experiment x one panel) and refuses a second root or workspace version. `db_export_facs(con, prep, exclusions)` is the FACS import path: drops excluded files, asks before overwriting a stored root, and in one transaction replaces the root, writes `import_log` (`qc_passed = TRUE`) and `import_exclusions`. Wells differing in one variable (e.g. restimulation) are separate stains via `condition` + `condition_keyword` (default `"none"`); ratios between wells are computed in analysis, not stored. No volumes/weights are stored -- they are calculation inputs. `source_file` is a path relative to the project's `data/` folder. Snapshots: `db_snapshot(con, root)` -> `<root>/<YYYY-MM-DD>/` (same-day rerun replaces it; con stays open), via DuckDB `EXPORT DATABASE`, which rejects self-referencing foreign keys -- don't add any. `db_query_*()` return lazy `dbplyr` tables keyed by `experiment_code`, `mouse_id`, `tissue`, `assay_name` -- `dplyr::collect()` to materialize.

---

## Dependency philosophy

- **`Imports`**: every package called inside an exported function body must be listed here. All calls must be namespaced `pkg::fn()`.
- **`Suggests`**: packages only needed in examples, vignettes, or tests -- or used optionally behind `requireNamespace()`. Currently `ggrepel` (QC report outlier labels; plain `geom_text` without it), `testthat` and `withr` (the latter used for `withr::defer()` cleanup in `db_` domain tests).
- **Non-CRAN packages** (`fcexpr`): list in `Imports` as normal; document the install requirement in README. Do not add `remotes::` calls inside function bodies.

Current `Imports`: `CytoML`, `DBI`, `dbplyr`, `diffcyt`, `dplyr`, `duckdb`, `fcexpr`, `flowCore`, `FlowSOM`, `flowWorkspace`, `ggplot2`, `glue`, `gt`, `gtsummary`, `janitor`, `parallel`, `purrr`, `readxl`, `rmarkdown`, `rstatix`, `stats`, `stringr`, `SummarizedExperiment`, `tibble`, `tidyr`, `tools`, `utils`, `uwot`, `xml2`, `xfun` (order matches `DESCRIPTION`'s `Imports:` block).

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
