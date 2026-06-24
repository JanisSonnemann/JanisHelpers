# Design: `facs_read_wsp()` — direct XML parser replacing `fcexpr`

**Date:** 2026-06-24
**Status:** Approved

---

## Goal

Replace the `fcexpr` dependency in the FACS import layer with a direct `xml2`-based parser. The new function `facs_read_wsp()` returns a structured list object containing all population data plus two metadata tables extracted from embedded FCS keywords. The existing `facs_import_wsp()` is preserved and marked deprecated.

---

## Architecture

Two files in `R/`:

| File | Contents |
|---|---|
| `facs_import.R` | `facs_import_wsp()` preserved as-is; `@deprecated` roxygen tag added; `.Deprecated()` call inserted at top of function body |
| `facs_read.R` | `facs_read_wsp()` (exported) + three unexported helpers: `parse_populations_()`, `parse_panel_()`, `parse_keywords_()` |

`xml2` added to `Imports` in `DESCRIPTION`. `fcexpr` remains in `Imports` (still required by the deprecated function). No other files change.

---

## Public API

```r
facs_read_wsp(path, group = NULL, keywords = NULL)
```

| Argument | Type | Default | Notes |
|---|---|---|---|
| `path` | chr | — | Path to `.wsp` file |
| `group` | chr | `NULL` | Group name to extract; `NULL` = all groups |
| `keywords` | chr vector | `NULL` | FCS keywords to join into `data`; missing ones filled `NA_character_` with a warning |

`r_stats` is removed — statistics are always extracted.

Returns a named list (visibly) with three slots: `data`, `meta`, `panel`.

---

## Output structure

### `data` — long tibble, one row per `FileName × PopulationFullPath × metric`

| Column | Type | Notes |
|---|---|---|
| `FileName` | chr | FCS filename as stored in workspace |
| `PopulationFullPath` | chr | Full gating hierarchy path |
| `Population` | chr | Leaf gate name (`basename(PopulationFullPath)`) |
| `metric` | chr | `"Count"`, `"FractionOfParent"`, or `"<Stat>_<Label>"` e.g. `"Median_CD19"` |
| `value` | dbl | Numeric measurement |
| `<keyword>` | chr | One column per entry in `keywords` arg, joined by `FileName` |

Stat labels use stain (`$PnS`) over channel (`$PnN`) via `dplyr::coalesce()`, matching current behaviour.

### `meta` — wide tibble, one row per file

Fixed set of system keywords always extracted, columns named without `$` prefix:

`DATE`, `BTIM`, `ETIM`, `CYT`, `INST`, `OP`, `TOT`

Missing keywords filled `NA_character_` silently (system keywords are optional in the FCS standard).

### `panel` — wide tibble, one row per file

One column per parameter, named by channel (`$PnN`), value = stain label (`$PnS`). Empty stain strings coerced to `NA_character_` (scatter/time channels typically have no stain).

---

## Internal parsing logic

`facs_read_wsp()` opens the XML once:

```r
doc <- xml2::read_xml(path)
```

The document object is passed to all three helpers. The file is read exactly once.

### `parse_populations_(doc, group)`

1. If `group` is not `NULL`, locate the matching `<GroupNode>` in `<Groups>`, collect its sample IDs from `<SampleRefs>`.
2. Iterate over `<Sample>` nodes in `<SampleList>`, filtering to the relevant sample IDs.
3. For each sample, extract `FileName` from the `<DataSet>` URI attribute.
4. Traverse the `<Population>` / `<Subpopulations>` tree recursively, accumulating the full path string (`"root/Lymphocytes/CD4"`) as it descends.
5. From each population node extract: count and fraction-of-parent (attributes), and any `<Statistic>` child nodes (Median, Mean, etc.).
6. Pivot all metrics to long format; bind rows across samples.

Returns a long tibble with columns: `FileName`, `PopulationFullPath`, `Population`, `metric`, `value`.

### `parse_panel_(doc)`

1. For each sample, find all `<Keyword>` nodes whose `name` attribute matches `$P[0-9]+[NS]`.
2. Extract parameter number, type (N or S), and value.
3. Pivot wide: channel name (`N`) becomes column name, stain label (`S`) becomes value.
4. Coerce empty stain strings to `NA_character_`.

Returns a wide tibble with columns: `FileName`, then one column per channel.

### `parse_keywords_(doc)`

1. For each sample, collect all `<Keyword>` nodes excluding `$P[0-9]+[NS]` (panel) and excluding the fixed meta set (`$DATE`, `$BTIM`, `$ETIM`, `$CYT`, `$INST`, `$OP`, `$TOT`).
2. Return as a long tibble (`FileName`, `key`, `value`).

Used two ways by the main function:
- Filtered to entries in the `keywords` arg → pivoted wide → joined into `data` by `FileName`
- The fixed meta set extracted separately (same XML pass) → `meta`

---

## Error handling

| Condition | Behaviour |
|---|---|
| Population node missing count attribute | Fill `NA_real_` |
| Requested `keywords` entry absent for a file | Fill `NA_character_`, issue `warning()` listing missing keys (current behaviour preserved) |
| `group` arg doesn't match any group in workspace | `stop()` with informative message listing available group names |
| Malformed XML | `xml2::read_xml()` raises its own error; not caught |

The legacy fallback (`fcexpr::wsx_get_popstats_legacy`) is not replicated — it handled a specific `fcexpr` bug with renamed FCS files, not a general XML parsing issue.

---

## Testing

A real `.wsp` fixture is required in `tests/fixtures/`. No mocking of XML parsing.

| File | Covers |
|---|---|
| `tests/testthat/test-facs_read.R` | `facs_read_wsp()` returns a named list with slots `data`, `meta`, `panel`; each is a tibble; `data` has required columns; `meta` has one row per file; `panel` has one row per file; requesting a missing keyword issues a warning |
| `tests/testthat/test-facs_import.R` | `facs_import_wsp()` issues a deprecation warning; output structure unchanged |

The fixture should include: at least one stat (to validate stat parsing), at least one custom keyword, and at least one parameter without a stain label (to validate `NA_character_` coercion in `panel`).

Internal helpers are unexported and tested indirectly through `facs_read_wsp()`.

---

## Dependencies

| Package | Change |
|---|---|
| `xml2` | Add to `Imports` |
| `fcexpr` | Keep in `Imports` (deprecated function still uses it) |
| `dplyr`, `tidyr`, `stringr`, `tibble`, `purrr`, `glue` | No change |
