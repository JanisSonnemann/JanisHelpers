# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package overview

`JanisHelpers` is a personal R utility package for biomedical research workflows at Charité. It wraps FlowJo cytometry data import, statistical analysis, and RMarkdown knitting into reusable functions. Installed from GitHub via `devtools::install_github("JanisSonnemann/JanisHelpers", auth_token = gh::gh_token())`.

## Development commands

```r
# Install/reinstall the package locally
devtools::install()

# Rebuild documentation (run after editing roxygen2 comments)
devtools::document()

# Check package
devtools::check()

# Load without installing (faster iteration)
devtools::load_all()
```

After editing any `@param`, `@returns`, or `@export` tags, always run `devtools::document()` to regenerate `man/*.Rd` and `NAMESPACE`.

## Code architecture

All exported functions live in `R/` across four files:

- **`flow-functions.R`** — FlowJo `.wsp` import pipeline. `import_workspace()` is the core function: it calls `fcexpr::wsx_get_popstats()` for cell counts/stats and `fcexpr::wsx_get_keywords()` for FCS header keywords, then joins them. `import_fcs()` and `import_fcs_clean()` wrap it with Excel caching (`xlsx::write.xlsx`). `import_workspace_long()` is the newer long-format variant that keeps data tidy (one row per file × population × metric) and is better suited for downstream `dplyr`/`ggplot2` workflows.

- **`analysis-functions.R`** — Statistical summary tables. `create_descriptive_table()` produces `gtsummary`/`gt` tables with Kruskal-Wallis p-values; `create_posthoc_tables()` performs Dunn's post-hoc on all significant numeric columns. Both accept an optional `tissue_col` that causes them to loop over tissues and return named lists.

- **`knitting_functions.R`** — RMarkdown rendering helpers. All functions wrap `rmarkdown::render()` and write dated output files (`basename-YYYY-MM-DD.html`). `knit_exp_structure()` and `knit_wide_html()` use `rstudioapi::getActiveDocumentContext()` to resolve output paths relative to the active script, so they require an RStudio session. `knit_wide_html()` injects `inst/resources/wide-output.css`.

- **`housekeeping.R`** — Single function `update_JanisHelpers_git()` that reinstalls the package from GitHub.

## Key dependencies

- **`fcexpr`**: non-CRAN package that parses `.wsp` files. All flow data entry points depend on it.
- `dplyr`, `tidyr`, `stringr`: used directly inside function bodies without `library()` — they must be in `Imports`.
- `gtsummary` + `gt`: table rendering in analysis functions.
- `rmarkdown` + `xfun` + `rstudioapi`: knitting helpers.

## RMarkdown templates

`inst/rmarkdown/templates/` provides four templates available in RStudio's *New File → R Markdown* dialog:
- `auto-knit-report` — uses `knit_multiple_dated` in the YAML `knit:` field
- `experiment-report` — uses `knit_exp_structure`
- `wide-html` — uses `knit_wide_html`
- `protocol-template` — standalone protocol skeleton

## Conventions

- Functions use base-pipe `|>` and tidy evaluation (`dplyr`/`tidyr` verbs) throughout.
- The `flow-functions.R` channel-label extraction relies on FCS keyword patterns `$P[0-9]+N` (channel name) and `$P[0-9]+S` (stain label); changes to that regex affect all stat column naming.
- `import_workspace_long()` prioritises stain label over channel name via `coalesce(stain, channel)` when building the `metric` column — keep this consistent if extending stat extraction.
