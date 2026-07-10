# CLAUDE.md

Personal R utility package for biomedical research workflows at Charité Berlin. Wraps three domains: (1) FlowJo `.wsp` cytometry data import via `fcexpr`, producing tidy long-format tibbles; (2) statistical summary and post-hoc tables via `gtsummary`/`gt`/`rstatix`; (3) RMarkdown knitting helpers that render dated HTML/PDF output into structured experiment directories. Installed from GitHub via `devtools::install_github("JanisSonnemann/JanisHelpers", auth_token = gh::gh_token())`.

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
- `wrangle_` — general data wrangling [stub — no functions yet]
- `db_` — database access [stub — no functions yet]

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
- **Non-ASCII**: escape all non-ASCII characters in R source with `\uXXXX` (e.g., em-dash = `—`).

---

## Data structures

### `facs_import_wsp()` — FlowJo workspace
- **Input**: `path` = path to `.wsp` file
- **Output**: long tibble, one row per `FileName × PopulationFullPath × metric`

| Column | Type | Notes |
|---|---|---|
| `FileName` | chr | FCS filename as stored in workspace |
| `PopulationFullPath` | chr | Full gating hierarchy path |
| `Population` | chr | Leaf gate name (`basename(PopulationFullPath)`) |
| `metric` | chr | `"Count"`, `"FractionOfParent"`, or `"<Stat>_<Label>"` e.g. `"Median_CD4"` |
| `value` | dbl | Numeric measurement |
| `<keyword>` | chr | One column per requested keyword (e.g. `mouse_ID`, `group`) |

- Stats label priority: stain label (`$PnS`) over channel name (`$PnN`) via `dplyr::coalesce(stain, channel)`.
- Keyword columns appended on the right; warnings issued for missing keywords (filled `NA_character_`).

### `analysis_summary_table()` / `analysis_posthoc_tables()` — stats
- **Input**: wide tibble with at minimum a `group` column and ≥1 numeric column. Optional `tissue` column triggers per-tissue subsetting.
- **Output**:
  - `analysis_summary_table()`: single `gt` table, or named list `tissue → gt` when `tissue_col` is set.
  - `analysis_posthoc_tables()`: named list `tissue → population → gt` (or flat `population → gt` without tissue). Non-significant variables return a plain character string instead of a table.

### `report_knit_*()` — knitting helpers
- **Input**: `input` = path to `.Rmd` file (usually the active document).
- **Output**: rendered HTML/PDF written to disk; returns the output path (from `rmarkdown::render()`).
- `report_knit_exp()` and `report_knit_wide()` use `normalizePath(input)` to resolve the output directory — no RStudio dependency, works from Positron or the terminal.

### `meta_read()` / `meta_clean()` / `meta_annotate()` — experiment metadata
- **`meta_read(path)`**: reads every sheet of an Excel workbook into a named
  list of cleaned tibbles (one per sheet, named after the sheet) -- blank
  rows/columns dropped, column names snake_cased (except `mouse_ID`,
  preserved verbatim), character columns trimmed, date columns coerced to
  `Date`, `group` coerced to a factor if present (each rule applied
  independently per sheet, only where applicable columns exist).
- **`meta_clean(meta_list)`**: takes the list from `meta_read()`, pivots its
  `organ_weights` element from wide-by-tissue (e.g. `kidney_total_weight`,
  `kidney_facs_weight`) to long (`mouse_ID`, `tissue`, `total_weight`,
  `facs_weight`), then joins it with `facs_volumes` (already long by
  `tissue`) via `meta_annotate(by = c("mouse_ID", "tissue"))`. Returns one
  tibble ready to pass as `meta` to `facs_calc_count_per_g()`. Errors if
  `organ_weights` or `facs_volumes` is missing from `meta_list`.
- **`meta_annotate(data, meta, by = "mouse_ID")`**: left-joins `meta` onto
  `data` by `by`, a character vector of one or more shared column names
  (e.g. `c("mouse_ID", "tissue")`). Errors if any `by` column is missing
  from either side or if non-`by` column names collide. Warns (rows kept,
  `NA` filled) if a `by` combination in `data` has no match in `meta`.

---

## Dependency philosophy

- **`Imports`**: every package called inside an exported function body must be listed here. All calls must be namespaced `pkg::fn()`.
- **`Suggests`**: packages only needed in examples or vignettes. Currently none.
- **Non-CRAN packages** (`fcexpr`): list in `Imports` as normal; document the install requirement in README. Do not add `remotes::` calls inside function bodies.

Current `Imports`: `CytoML`, `diffcyt`, `dplyr`, `fcexpr`, `flowCore`, `FlowSOM`, `flowWorkspace`, `ggplot2`, `glue`, `gt`, `gtsummary`, `janitor`, `parallel`, `purrr`, `readxl`, `rmarkdown`, `rstatix`, `stats`, `stringr`, `SummarizedExperiment`, `tibble`, `tidyr`, `tools`, `uwot`, `xml2`, `xfun` (order matches `DESCRIPTION`'s `Imports:` block).

---

## Testing conventions

- No tests exist yet. When adding tests, use `testthat` (`usethis::use_testthat()`).
- Test files mirror source: `tests/testthat/test-facs_import.R`, `test-analysis_stats.R`, etc.
- `facs_import_wsp()` requires a real `.wsp` file — use a minimal fixture in `tests/fixtures/`. Do not mock `fcexpr` calls.
- [TODO: define fixture `.wsp` file and add at least one smoke test per domain]

---

## RMarkdown templates

`inst/rmarkdown/templates/` — available in RStudio *New File → R Markdown*:

| Template dir | Knit function |
|---|---|
| `auto-knit-report` | `report_knit_dated()` via YAML `knit:` field |
| `experiment-report` | `report_knit_exp()` |
| `wide-html` | `report_knit_wide()` |
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

`devtools::check()` currently produces **0 errors, 0 warnings**, plus a small set of pre-existing notes:

- `future file timestamps` — network/environment issue, not code.
- `checking for hidden files and directories` (`.git`) — only observed when running `check()` from inside a git worktree, where `.git` is a file (pointing at the real gitdir) rather than a directory; R's default VCS-ignore logic doesn't recognize the file form. Does not appear when checking from a normal clone.
- `checking top-level files` (`docs`) — the `docs/superpowers/` directory (planning/spec notes from the `sdd` workflow, not package documentation) isn't covered by `.Rbuildignore`. Pre-existing since the `facs_read_wsp` work; not introduced or addressed by the `meta_` domain work.
- `R code for possible problems` — remaining items are all tidy-eval / gt false positives, not bare calls:
  - `%>%` in `analysis_stats.R` — kept intentionally in the gtsummary method-dispatch chain (see pipe rule above).
  - `.data` in `analysis_stats.R` — the rlang tidy-eval pronoun; required for dynamic column selection, cannot be namespaced.
  - `p.adj`, `p.value` in `analysis_stats.R` — bare column names inside gt `rows =` predicates; this is how gt row selection works.
  - `all_pops_fun` multiple definitions in `analysis_posthoc_tables` — two branches define the function with different signatures; structural issue, low priority.
  - All `facs_import_wsp` variable-binding notes (`FileName`, `PopulationFullPath`, etc.) — bare column names inside dplyr verbs; valid tidy eval, flagged as a static-analysis false positive.
  - All `meta_read`/`meta_annotate` variable-binding notes (`mouse_id`, `mouse_ID`, `group`) — bare column names inside dplyr verbs; valid tidy eval, flagged as a static-analysis false positive.
  - `meta_clean`/`pivot_organ_weights_long_` variable-binding notes (`tissue`, `tissue_variant_`, `mouse_ID`) — bare column names inside dplyr verbs, same false-positive pattern as above. `meta_clean()` derives a base organ from each `facs_volumes` tissue value (substring before the first hyphen, via the unexported `extract_base_tissue_()`) and joins `organ_weights` on that base organ instead of the literal tissue string, so a tissue variant of the same organ with no dedicated `organ_weights` columns (e.g. `kidney-whole`, an aliquot taken before further processing) reuses that organ's weight instead of requiring duplicate `organ_weights` columns.
  - All `facs_calc_pct_of`/`facs_calc_count_per_g` variable-binding notes (`population`, `metric`, `value`, `file_name`, `population_full_path`, `ref_count`, `mouse_ID`, `tissue`, `method`, `bead_count`, `.env`) — bare column names inside dplyr verbs; valid tidy eval, flagged as a static-analysis false positive. `tissue` is a bare column reference (not a formal argument, unlike in the old single-tissue signature) — `facs_calc_count_per_g()` now joins `meta` on `mouse_ID` + `tissue` and processes every tissue present in `data` in one call. Note: the function's own formal arguments (`vol_total`, `vol_stained`, `vol_resuspended`, `vol_measured`, `organ_piece_weight`, `method_col`, `bead_pop`, `bead_concentration`) are also referenced bare inside `dplyr` verbs but are not flagged, since codetools resolves them against the matching parameter name in scope.
  - `pivot_organ_weights_long_` (the internal `meta_clean()` helper that pivots `organ_weights` to long format) variable-binding note (`mouse_ID`) — bare column name inside a `tidyr::pivot_longer()` deselection (`cols = !mouse_ID`); valid tidy eval, flagged as a static-analysis false positive, same pattern as the existing `meta_read`/`meta_annotate` note above. (`meta_clean_sheet_`, the internal `meta_read()` per-sheet helper, also produces `mouse_id`/`group` notes -- already covered by the `meta_read`/`meta_annotate` bullet above.)
  - All `calc_restim_proportions_` variable-binding notes (`population`, `metric`, `value`, `n`, `ref_count`, `proportion`, `.data`) — bare column names inside dplyr verbs, same false-positive pattern as `facs_calc_pct_of`/`facs_calc_count_per_g` above. `calc_restim_proportions_` is unexported (an internal helper shared by both `facs_calc_log2fc()` and `facs_calc_diff()`), so its notes are attributed to `R/facs_calc.R` as a whole rather than one function. Its nested `is_constant_` closure (used to detect passthrough columns constant within each `mouse_ID` x `tissue` group) adds its own notes (`mouse_ID`, `tissue`, `.data`, `n_distinct`) for the same bare-column-reference reason. `facs_calc_log2fc`'s and `facs_calc_diff`'s own top-level notes are limited to `.data`, from the anonymous `function(x) ...` closure each passes to `dplyr::mutate(dplyr::across(...))` to compute the log2-fold-change/difference column.
  - `facs_read_fcs_gated` variable-binding notes (`key`, `value`) — bare column names inside `dplyr::filter()`/`tidyr::pivot_wider()`, same false-positive pattern as `facs_read_wsp`/`parse_keywords_` above. `facs_read_fcs_gated()` reads its `keywords` argument from the `.wsp` XML via the shared `parse_keywords_()` helper (not from `flowWorkspace::keyword()`, which only sees keywords physically embedded in the raw `.fcs` file's TEXT segment and misses keywords typed into FlowJo's UI but never written back to the `.fcs` file) — this is the same long-tibble (`file_name`, `key`, `value`) shape `facs_read_wsp()` already pivots, hence the identical notes.
  - `facs_calc_cluster_freq` variable-binding notes (`file_name`, `.data`, `n`, `total`) — bare column names inside `dplyr`/`tidyr` verbs, same false-positive pattern as `facs_calc_pct_of`/`facs_calc_count_per_g` above. Its nested `is_constant_` closure (used to detect passthrough columns constant within each `file_name` group -- the single-tissue-per-file analog of `calc_restim_proportions_`'s own `is_constant_` helper) adds its own notes (`file_name`, `.data`, `n_distinct`) for the same bare-column-reference reason. `facs_cluster_flowsom()`, the FlowSOM clustering function whose output this helper summarizes, produces no `R CMD check` notes of its own.
  - `facs_plot_umap` variable-binding notes (`UMAP1`, `UMAP2`, `.data`) — bare column names and the rlang tidy-eval pronoun inside `ggplot2::aes()`, same false-positive pattern as `.data` usage elsewhere in the package (see the `analysis_stats.R` bullet above). `facs_reduce_umap()` (Stage 3's other new function, which computes the `UMAP1`/`UMAP2` columns `facs_plot_umap()` consumes) produces no `R CMD check` notes of its own -- it accesses columns via `$`/`[[` subsetting rather than bare tidy-eval names.
  - `facs_test_cluster_abundance` variable-binding notes (`file_name`, `:=`, `cluster_id`) — `file_name` is a bare column name inside `dplyr::rename(sample_id = file_name)`/`dplyr::select(dplyr::all_of(...))`/`tidyr::pivot_wider()`, same false-positive pattern as elsewhere. `:=` is the rlang walrus operator used in the closing `dplyr::rename(results, !!cluster_col := cluster_id)` call (renames the `diffcyt` result column back to the user-supplied `cluster_col` name) -- a standard dplyr/rlang non-standard-evaluation idiom, not an actually-undefined function. `cluster_id` is the bare right-hand-side reference in that same `rename()` call, matching the `cluster_id` column name `diffcyt::testDA_*()` returns in its `rowData`.
  - `facs_plot_cluster_abundance` variable-binding notes (`p_adj`, `.data`, `fraction`) — `p_adj` is a bare column name inside `dplyr::filter(p_adj <= p_adj_threshold)`; `.data` is the rlang tidy-eval pronoun used in `.data[[cluster_col]]`/`.data[[group_col]]` to subset by the user-supplied cluster/grouping column names; `fraction` is a bare column name inside `ggplot2::aes(y = fraction)` -- same false-positive pattern as `facs_plot_umap` above.
