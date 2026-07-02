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

---

## Dependency philosophy

- **`Imports`**: every package called inside an exported function body must be listed here. All calls must be namespaced `pkg::fn()`.
- **`Suggests`**: packages only needed in examples or vignettes. Currently none.
- **Non-CRAN packages** (`fcexpr`): list in `Imports` as normal; document the install requirement in README. Do not add `remotes::` calls inside function bodies.

Current `Imports`: `fcexpr`, `dplyr`, `tidyr`, `stringr`, `tibble`, `purrr`, `glue`, `rmarkdown`, `xfun`, `gt`, `gtsummary`, `rstatix`.

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

`devtools::check()` currently produces **0 errors, 0 warnings, 2 notes**. Both notes are pre-existing:

- `future file timestamps` — network/environment issue, not code.
- `R code for possible problems` — remaining items are all tidy-eval / gt false positives, not bare calls:
  - `%>%` in `analysis_stats.R` — kept intentionally in the gtsummary method-dispatch chain (see pipe rule above).
  - `.data` in `analysis_stats.R` — the rlang tidy-eval pronoun; required for dynamic column selection, cannot be namespaced.
  - `p.adj`, `p.value` in `analysis_stats.R` — bare column names inside gt `rows =` predicates; this is how gt row selection works.
  - `all_pops_fun` multiple definitions in `analysis_posthoc_tables` — two branches define the function with different signatures; structural issue, low priority.
  - All `facs_import_wsp` variable-binding notes (`FileName`, `PopulationFullPath`, etc.) — bare column names inside dplyr verbs; valid tidy eval, flagged as a static-analysis false positive.
