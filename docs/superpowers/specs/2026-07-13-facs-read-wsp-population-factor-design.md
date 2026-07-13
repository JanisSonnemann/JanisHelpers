# Design: `population` as gating-order factor in `facs_read_wsp()`

Date: 2026-07-13

## Problem

`facs_read_wsp()`'s `data$population` column (the leaf gate name, e.g.
`"CD4+"`) is currently a plain character vector. Downstream consumers
(plots, `gt`/`gtsummary` tables) sort it alphabetically, which does not
match the biologically meaningful order of the FlowJo gating hierarchy
(e.g. `Singlets` -> `non-debris-` -> `CD45+` -> `CD4+`).

## Goal

Make `data$population` a plain `factor` (not R's `ordered` class) whose
level order follows the depth-first traversal order of the gating tree,
as encoded by `population_full_path`.

## Why plain factor, not `ordered`

`population` feeds into `gtsummary`/`rstatix` calls in `analysis_stats.R`.
An `ordered` factor signals ordinal data to R and changes default
behavior there (e.g. polynomial contrasts). `population` is a nominal
label, not an ordinal one — a plain `factor()` with a fixed level order
gives correct display/sort order in plots and tables without that side
effect.

## Why leaf name (`population`), not `population_full_path`

User requirement: sort the leaf name by the order of its
`population_full_path`. `population_full_path` remains a plain character
column, available for callers who need unambiguous identity when leaf
names collide across branches (see Edge Cases below).

## Implementation

`walk_pops_()` already performs a depth-first pre-order traversal of the
gating tree: a population's own `count`/`fraction_of_parent` rows
(`base_rows`, `R/facs_read.R:111-117`) are always emitted before its
descendants' rows (the recursive call at `R/facs_read.R:119-122`). This
means the row order of the assembled `pops` tibble already encodes
gating-hierarchy traversal order — no changes are needed to
`walk_pops_()` or `parse_populations_()`.

Add one small internal helper in `R/facs_read.R`, following the existing
`_`-suffix naming convention for internal helpers in that file
(`filter_samples_`, `parse_keywords_`, `parse_panel_`, `walk_pops_`,
`parse_populations_`):

```r
population_gating_order_ <- function(pops) {
  full_path_order <- unique(pops$population_full_path)
  unique(basename(full_path_order))
}
```

Apply it in `facs_read_wsp()` immediately after `data <- pops`
(`R/facs_read.R:295`), before the optional keyword join (ordering does
not depend on keywords):

```r
data$population <- factor(data$population, levels = population_gating_order_(pops))
```

## Edge cases

- **Empty result** (e.g. a `group` filter that matches zero samples):
  `pops` has 0 rows, so `levels` is `character(0)` and
  `factor(character(0), levels = character(0))` succeeds with no error.
- **Leaf-name collisions**: if the same leaf name (e.g. `"CD4+"`) occurs
  under two different parents in the tree, it collapses to a single
  factor level, positioned at its *first* traversal occurrence. This is
  an inherent consequence of ordering by leaf name rather than full path.
  Documented explicitly in roxygen so it isn't a silent surprise;
  `population_full_path` remains available for disambiguation.
- **Group filtering**: order is derived from whatever `pops` rows survive
  the `group`/`sample_ids` filter already applied earlier in
  `facs_read_wsp()`, so the factor levels are always self-consistent with
  the returned `data` — never based on unfiltered data.

## Downstream impact

Checked all uses of `population` in `R/facs_calc.R` (the only other file
referencing it, besides docs):

- `population == ref_pop`, `population != bead_pop`, etc. — `==`/`!=`
  work transparently between a factor and a character scalar.
- `paste(...)`, `purrr::pmap_chr(...)` over rows including a `population`
  factor column — `paste()` calls `as.character()` internally, which for
  factors returns the level label, not the integer code.
- `dplyr::group_by()`, `tidyr::pivot_wider(id_cols = ...)` with
  `population` as a grouping/id column — both work normally with factor
  columns.

No breakage expected. `facs_read_fcs.R` (`facs_read_fcs_gated()`) does
not produce a `population` column, so it is out of scope.

## Documentation updates

- `facs_read_wsp()` roxygen `@returns` (`R/facs_read.R:207-221`): update
  the `data` slot's `population` column description to state it is a
  `factor`, levels ordered by first-encountered depth-first traversal of
  the gating hierarchy (matching `population_full_path` order), including
  the leaf-collision caveat.
- `CLAUDE.md` "Data structures" section: add one line noting `population`
  is a gating-ordered factor, using the correct current column name
  (`population`, snake_case). Note: the existing table in that section
  already has stale camelCase column names (`PopulationFullPath`,
  `Population`) predating this change — that drift is pre-existing and
  out of scope here.

## Testing plan

New tests in `tests/testthat/test-facs_read.R`:

1. `result$data$population` is `is.factor() == TRUE` and
   `is.ordered() == FALSE`.
2. Level order matches
   `unique(basename(unique(result$data$population_full_path)))` computed
   independently from the same fixture result.
3. Using the fixture's known nested boolean-gate chain (`Singlets` ->
   `Singlets/non-debris-` -> `Singlets/non-debris-/CD45+` ->
   `.../CD4+`, already exercised by the tests at
   `tests/testthat/test-facs_read.R:116-148`): assert each ancestor's
   factor level index is less than its descendant's.
4. Run the full `testthat` suite plus `devtools::check()` to confirm no
   regressions in `facs_calc.R`/`analysis_stats.R` consumers (0 errors, 0
   warnings target per `CLAUDE.md`).
