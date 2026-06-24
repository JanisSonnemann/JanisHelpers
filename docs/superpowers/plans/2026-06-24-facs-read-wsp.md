# facs_read_wsp() Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the `fcexpr` dependency in the FACS import layer with a direct `xml2`-based parser in a new `facs_read_wsp()` function, while preserving `facs_import_wsp()` as deprecated.

**Architecture:** Three unexported helpers (`parse_keywords_()`, `parse_panel_()`, `parse_populations_()`) each receive the parsed `xml_document` and return tidy tibbles. `facs_read_wsp()` reads the XML once, passes the document to all three, assembles `list(data, meta, panel)`, and handles group filtering and keyword joining. The old `facs_import_wsp()` is kept untouched in `facs_import.R` with only a `.Deprecated()` call added.

**Tech Stack:** R ≥ 4.1, `xml2` (new), `dplyr`, `tidyr`, `stringr`, `tibble`, `purrr`, `glue` (all existing in Imports)

## Global Constraints

- Base pipe `|>` everywhere; no `%>%` in new code
- All external calls namespaced as `pkg::fn()` — no bare `filter()`, `map()`, `select()`, etc.
- Output column names are snake_case: `file_name`, `population_full_path`, `population`, `metric`, `value`
- `metric` values are lowercase: `"count"`, `"fraction_of_parent"`, `"median_CD19"` (stat name lowercase, label case-preserved)
- R CMD check target: 0 errors, 0 warnings (2 pre-existing notes are acceptable)
- `devtools::document()` after any roxygen change; `devtools::check()` before final commit
- Fixture path: `tests/fixtures/minimal.wsp` (must be present before Task 2)
- Working directory for all commands: the package root

---

## File Map

| Action | File | Purpose |
|---|---|---|
| Modify | `DESCRIPTION` | Add `xml2` to `Imports` |
| Modify | `R/facs_import.R` | Add `.Deprecated()` call at top of function body |
| Create | `R/facs_read.R` | `facs_read_wsp()` + `filter_samples_()`, `parse_keywords_()`, `parse_panel_()`, `parse_populations_()`, `walk_pops_()` |
| Modify | `tests/testthat/test-facs_import.R` | Add deprecation test |
| Create | `tests/testthat/test-facs_read.R` | Integration tests for `facs_read_wsp()` |

---

### Task 1: Add xml2 dependency + discover WSP XML structure

**Files:**
- Modify: `DESCRIPTION`
- (Read-only) `tests/fixtures/minimal.wsp`

**Interfaces:**
- Produces: `xml2` in Imports; a comment block in `R/facs_read.R` documenting the exact element/attribute names found in the fixture, used in Tasks 4–5.

- [ ] **Step 1: Add xml2 to DESCRIPTION Imports**

Open `DESCRIPTION` and add `xml2` to the `Imports:` block in alphabetical position:

```
Imports:
    dplyr,
    fcexpr,
    glue,
    gt,
    gtsummary,
    purrr,
    rmarkdown,
    rstatix,
    stringr,
    tibble,
    tidyr,
    xml2,
    xfun
```

- [ ] **Step 2: Verify fixture exists**

```bash
ls tests/fixtures/minimal.wsp
```

Expected: file listed. If absent, create it from FlowJo 10 following `tests/fixtures/README.md` before continuing.

- [ ] **Step 3: Run XML discovery script**

In an R session from the package root:

```r
devtools::load_all()
library(xml2)

doc <- xml2::read_xml("tests/fixtures/minimal.wsp")

# 1. Check for XML namespaces — note any listed
xml2::xml_ns(doc)

# 2. Top-level children
xml2::xml_name(xml2::xml_children(doc))

# 3. Samples
samples <- xml2::xml_find_all(doc, ".//SampleList/Sample")
cat("Samples found:", length(samples), "\n")

# 4. Direct children of first sample (shows node names)
xml2::xml_name(xml2::xml_children(samples[[1]]))

# 5. Keywords (expect name/value attributes on Keyword nodes)
kw <- xml2::xml_find_all(samples[[1]], ".//keywords/Keyword")
head(data.frame(name  = xml2::xml_attr(kw, "name"),
                value = xml2::xml_attr(kw, "value")), 20)

# 6. Population tree: check root node name and count attribute
for (child in xml2::xml_children(samples[[1]])) {
  cat("Node:", xml2::xml_name(child),
      "| name attr:", xml2::xml_attr(child, "name"),
      "| count attr:", xml2::xml_attr(child, "count"), "\n")
}

# 7. Population sub-nodes: what are child populations called?
root_node <- xml2::xml_children(samples[[1]])[[3]]  # adjust index if needed
xml2::xml_structure(root_node)

# 8. Statistic nodes
stat_nodes <- xml2::xml_find_all(samples[[1]], ".//Statistic")
if (length(stat_nodes) > 0) {
  cat("Statistic node count:", length(stat_nodes), "\n")
  cat("Attributes of first stat:\n")
  print(xml2::xml_attrs(stat_nodes[[1]]))
} else {
  cat("No <Statistic> nodes found — try .//stat or .//Statistics/Statistic\n")
}

# 9. Group structure
groups <- xml2::xml_find_all(doc, ".//Groups/GroupNode")
cat("Groups:", xml2::xml_attr(groups, "name"), "\n")
sample_refs <- xml2::xml_find_all(groups[[1]], ".//SampleRefs/SampleRef")
cat("SampleRef ID attribute:", xml2::xml_attr(sample_refs[[1]], "sampleID"), "\n")
cat("DataSet ID attribute:", xml2::xml_attr(
  xml2::xml_find_first(samples[[1]], "DataSet"), "sampleID"), "\n")
```

- [ ] **Step 4: Create facs_read.R with discovery notes**

Create `R/facs_read.R` containing only a comment block filled in from Step 3 output:

```r
# FlowJo WSP XML structure — verified against tests/fixtures/minimal.wsp
#
# Root sample tree node name (child of <Sample>):          e.g. "SampleNode"
# Population container element:                            e.g. "Subpopulations"
# Individual population element:                           e.g. "PopulationNode"
# Population count attribute:                              e.g. "count"
# Population name attribute:                               e.g. "name"
# Statistic element name:                                  e.g. "Statistic"
# Statistic type attribute (Median/Mean/etc.):             e.g. "name"
# Statistic channel attribute:                             e.g. "channel"
# Statistic value attribute:                               e.g. "value"
# SampleRef ID attribute:                                  e.g. "sampleID"
# DataSet sample ID attribute:                             e.g. "sampleID"
# Namespace prefix required for XPath (none if ns empty):  e.g. none
```

Fill every "e.g." with what Step 3 actually showed.

- [ ] **Step 5: Commit**

```bash
git add DESCRIPTION R/facs_read.R
git commit -m "chore: add xml2 to Imports and document WSP XML structure"
```

---

### Task 2: Deprecate facs_import_wsp() and add deprecation test

**Files:**
- Modify: `R/facs_import.R`
- Modify: `tests/testthat/test-facs_import.R`

**Interfaces:**
- Produces: `facs_import_wsp()` that emits a deprecation warning before executing its existing logic (no other behaviour change)

- [ ] **Step 1: Write the failing test**

Add to `tests/testthat/test-facs_import.R`:

```r
test_that("facs_import_wsp() emits a deprecation warning", {
  skip_if_not(
    file.exists(testthat::test_path("../fixtures/minimal.wsp")),
    "WSP fixture not available"
  )
  expect_warning(
    suppressMessages(
      facs_import_wsp(testthat::test_path("../fixtures/minimal.wsp"))
    ),
    regexp = "deprecated",
    ignore.case = TRUE
  )
})
```

- [ ] **Step 2: Run the test to verify it fails**

```r
devtools::load_all()
testthat::test_file("tests/testthat/test-facs_import.R")
```

Expected: FAIL — no deprecation warning found.

- [ ] **Step 3: Add .Deprecated() to facs_import_wsp()**

In `R/facs_import.R`, insert `.Deprecated("facs_read_wsp")` as the very first line of the function body, and add a deprecation note to the roxygen block:

```r
#' Import FlowJo workspace data in long format
#'
#' @description
#' **Deprecated.** Use \code{\link{facs_read_wsp}} instead.
#'
#' @param path path to .wsp file
#' @param group group to extract from workspace, default = NULL (all groups)
#' @param r_stats logical, whether to extract statistics such as MFI, default = FALSE
#' @param keywords character vector of FCS keywords to attach to each row, default = NULL
#'
#' @returns tibble in long format with one row per file x population x metric,
#'   returned invisibly -- assign the result explicitly
#' @export
#'
#' @examples
#' \dontrun{
#'   dat <- facs_import_wsp(
#'     path     = "experiment.wsp",
#'     group    = "Spleen",
#'     r_stats  = TRUE,
#'     keywords = c("mouse_ID", "group")
#'   )
#' }
facs_import_wsp <- function(path, group = NULL, r_stats = FALSE, keywords = NULL) {
  .Deprecated("facs_read_wsp")

  # Import raw workspace; fall back to legacy importer if FCS files were renamed after export
  ps_raw <- tryCatch(
  # ... rest of function body unchanged ...
```

- [ ] **Step 4: Regenerate docs and run test**

```r
devtools::document()
devtools::load_all()
testthat::test_file("tests/testthat/test-facs_import.R")
```

Expected: new test PASS; existing tests unaffected (they use `suppressWarnings()` if needed).

- [ ] **Step 5: Commit**

```bash
git add R/facs_import.R man/facs_import_wsp.Rd tests/testthat/test-facs_import.R
git commit -m "feat: deprecate facs_import_wsp() in favour of facs_read_wsp()"
```

---

### Task 3: Scaffold facs_read.R and write all integration tests (failing)

**Files:**
- Modify: `R/facs_read.R`
- Create: `tests/testthat/test-facs_read.R`

**Interfaces:**
- Produces: stub implementations (all `stop("not yet implemented")`); complete test suite that all fail

- [ ] **Step 1: Add stubs and roxygen to facs_read.R**

Append to `R/facs_read.R` (after the comment block from Task 1):

```r
# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

filter_samples_ <- function(samples, ids) {
  stop("not yet implemented")
}

parse_keywords_ <- function(doc, sample_ids = NULL) {
  stop("not yet implemented")
}

parse_panel_ <- function(doc, sample_ids = NULL) {
  stop("not yet implemented")
}

walk_pops_ <- function(node, file_name, path, parent_count, stain_lookup) {
  stop("not yet implemented")
}

parse_populations_ <- function(doc, sample_ids = NULL) {
  stop("not yet implemented")
}

# ---------------------------------------------------------------------------
# Exported function
# ---------------------------------------------------------------------------

#' Read a FlowJo workspace into a structured list
#'
#' Parses a FlowJo \code{.wsp} file directly via \code{xml2} with no dependency
#' on \code{fcexpr}. All population statistics are always extracted.
#'
#' @param path path to \code{.wsp} file
#' @param group character; group name to extract. \code{NULL} (default) extracts
#'   all groups.
#' @param keywords character vector of FCS keyword names to join into
#'   \code{data}. Keywords absent from the workspace are filled
#'   \code{NA_character_} with a warning.
#'
#' @returns A named list with three elements:
#'   \describe{
#'     \item{\code{data}}{Long-format tibble, one row per
#'       \code{file_name \times population_full_path \times metric}.
#'       Columns: \code{file_name}, \code{population_full_path},
#'       \code{population}, \code{metric}, \code{value}, plus any
#'       requested \code{keywords}.}
#'     \item{\code{meta}}{Wide-format tibble, one row per file.
#'       Columns: \code{file_name}, \code{DATE}, \code{BTIM},
#'       \code{ETIM}, \code{CYT}, \code{INST}, \code{OP}, \code{TOT}.
#'       Missing system keywords filled \code{NA_character_}.}
#'     \item{\code{panel}}{Wide-format tibble, one row per file.
#'       One column per cytometer parameter (channel name); value is
#'       the stain label, or \code{NA} for unlabelled channels.}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#'   res <- facs_read_wsp(
#'     path     = "experiment.wsp",
#'     group    = "Spleen",
#'     keywords = c("mouse_ID", "group")
#'   )
#'   res$data
#'   res$meta
#'   res$panel
#' }
facs_read_wsp <- function(path, group = NULL, keywords = NULL) {
  stop("not yet implemented")
}
```

- [ ] **Step 2: Create test-facs_read.R**

Create `tests/testthat/test-facs_read.R`:

```r
library(testthat)
library(JanisHelpers)

wsp_path <- testthat::test_path("../fixtures/minimal.wsp")
skip_msg  <- "WSP fixture not available"

test_that("facs_read_wsp() returns a named list with slots data, meta, panel", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))
  expect_type(result, "list")
  expect_named(result, c("data", "meta", "panel"))
})

test_that("data slot is a tibble with required columns", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))
  expect_s3_class(result$data, "tbl_df")
  expect_true(all(
    c("file_name", "population_full_path", "population", "metric", "value")
    %in% names(result$data)
  ))
})

test_that("data metric contains count and fraction_of_parent", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))
  expect_true("count" %in% result$data$metric)
  expect_true("fraction_of_parent" %in% result$data$metric)
})

test_that("data value column is numeric", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))
  expect_type(result$data$value, "double")
})

test_that("meta slot is a wide tibble with one row per file", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))
  expect_s3_class(result$meta, "tbl_df")
  n_files <- dplyr::n_distinct(result$data$file_name)
  expect_equal(nrow(result$meta), n_files)
  expect_true("file_name" %in% names(result$meta))
})

test_that("panel slot is a wide tibble with one row per file", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))
  expect_s3_class(result$panel, "tbl_df")
  n_files <- dplyr::n_distinct(result$data$file_name)
  expect_equal(nrow(result$panel), n_files)
  expect_true("file_name" %in% names(result$panel))
})

test_that("requested keywords are joined into data", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  doc <- xml2::read_xml(wsp_path)
  kw_nodes <- xml2::xml_find_all(doc, ".//SampleList/Sample[1]//keywords/Keyword")
  all_kws  <- xml2::xml_attr(kw_nodes, "name")
  user_kw  <- all_kws[!grepl("^\\$", all_kws)][1]
  skip_if(is.na(user_kw), "No user-level keywords in fixture")

  result <- suppressMessages(facs_read_wsp(wsp_path, keywords = user_kw))
  expect_true(user_kw %in% names(result$data))
})

test_that("missing keyword warns and fills NA", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  expect_warning(
    result <- suppressMessages(
      facs_read_wsp(wsp_path, keywords = "KW_DOES_NOT_EXIST_XYZ")
    ),
    regexp = "not found",
    ignore.case = TRUE
  )
  expect_true("KW_DOES_NOT_EXIST_XYZ" %in% names(result$data))
  expect_true(all(is.na(result$data$KW_DOES_NOT_EXIST_XYZ)))
})

test_that("group filtering restricts to that group's samples", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  doc    <- xml2::read_xml(wsp_path)
  groups <- xml2::xml_attr(xml2::xml_find_all(doc, ".//Groups/GroupNode"), "name")
  skip_if(length(groups) == 0L, "No groups in fixture")

  result_all   <- suppressMessages(facs_read_wsp(wsp_path))
  result_group <- suppressMessages(facs_read_wsp(wsp_path, group = groups[[1]]))
  expect_lte(
    dplyr::n_distinct(result_group$data$file_name),
    dplyr::n_distinct(result_all$data$file_name)
  )
})

test_that("invalid group name stops with informative error", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  expect_error(
    facs_read_wsp(wsp_path, group = "GROUP_DOES_NOT_EXIST_XYZ"),
    regexp = "not found",
    ignore.case = TRUE
  )
})
```

- [ ] **Step 3: Run tests to confirm they all fail**

```r
devtools::document()
devtools::load_all()
testthat::test_file("tests/testthat/test-facs_read.R")
```

Expected: all tests FAIL with "not yet implemented".

- [ ] **Step 4: Commit**

```bash
git add R/facs_read.R man/facs_read_wsp.Rd tests/testthat/test-facs_read.R
git commit -m "test: scaffold facs_read_wsp integration tests (all failing)"
```

---

### Task 4: Implement filter_samples_(), parse_keywords_(), parse_panel_()

**Files:**
- Modify: `R/facs_read.R`

**Interfaces:**
- Consumes: XML node/attribute names from Task 1; `xml2::xml_find_all`, `xml2::xml_find_first`, `xml2::xml_attr`
- Produces:
  - `filter_samples_(samples, ids)` → `xml_nodeset` (subset of samples matching IDs, or all if `ids` is NULL)
  - `parse_keywords_(doc, sample_ids)` → long tibble `(file_name chr, key chr, value chr)`, panel keys excluded
  - `parse_panel_(doc, sample_ids)` → wide tibble `(file_name chr, <channel_name> chr, ...)` one row per file

- [ ] **Step 1: Replace filter_samples_() stub**

```r
filter_samples_ <- function(samples, ids) {
  if (is.null(ids)) return(samples)
  found_ids <- purrr::map_chr(samples, function(s) {
    xml2::xml_attr(xml2::xml_find_first(s, "DataSet"), "sampleID")
  })
  samples[found_ids %in% ids]
}
```

- [ ] **Step 2: Replace parse_keywords_() stub**

```r
parse_keywords_ <- function(doc, sample_ids = NULL) {
  PANEL_PAT <- "^\\$P[0-9]+[NS]$"
  samples   <- filter_samples_(
    xml2::xml_find_all(doc, ".//SampleList/Sample"),
    sample_ids
  )

  purrr::map(samples, function(sample) {
    file_name <- basename(
      xml2::xml_attr(xml2::xml_find_first(sample, "DataSet"), "uri")
    )
    kw_nodes <- xml2::xml_find_all(sample, ".//keywords/Keyword")
    tibble::tibble(
      file_name = file_name,
      key       = xml2::xml_attr(kw_nodes, "name"),
      value     = xml2::xml_attr(kw_nodes, "value")
    ) |>
      dplyr::filter(!grepl(PANEL_PAT, key))
  }) |>
    dplyr::bind_rows()
}
```

- [ ] **Step 3: Replace parse_panel_() stub**

```r
parse_panel_ <- function(doc, sample_ids = NULL) {
  samples <- filter_samples_(
    xml2::xml_find_all(doc, ".//SampleList/Sample"),
    sample_ids
  )

  purrr::map(samples, function(sample) {
    file_name <- basename(
      xml2::xml_attr(xml2::xml_find_first(sample, "DataSet"), "uri")
    )
    kw_nodes <- xml2::xml_find_all(sample, ".//keywords/Keyword")

    panel_wide <- tibble::tibble(
      key   = xml2::xml_attr(kw_nodes, "name"),
      value = xml2::xml_attr(kw_nodes, "value")
    ) |>
      dplyr::filter(grepl("^\\$P[0-9]+[NS]$", key)) |>
      dplyr::mutate(
        number = stringr::str_extract(key, "(?<=\\$P)[0-9]+"),
        type   = stringr::str_extract(key, "[NS]$")
      ) |>
      dplyr::select(!key) |>
      tidyr::pivot_wider(
        names_from  = type,
        values_from = value,
        values_fill = NA_character_
      ) |>
      dplyr::rename(channel = N, stain = S) |>
      dplyr::mutate(stain = dplyr::na_if(stain, "")) |>
      dplyr::select(!number) |>
      tidyr::pivot_wider(names_from = channel, values_from = stain)

    tibble::add_column(panel_wide, file_name = file_name, .before = 1L)
  }) |>
    dplyr::bind_rows()
}
```

- [ ] **Step 4: Sanity-check interactively**

```r
devtools::load_all()
doc <- xml2::read_xml("tests/fixtures/minimal.wsp")

kws <- JanisHelpers:::parse_keywords_(doc)
dplyr::glimpse(kws)
# Expect: columns file_name, key, value; no rows with key matching $P[0-9]+[NS]

pnl <- JanisHelpers:::parse_panel_(doc)
dplyr::glimpse(pnl)
# Expect: one row per sample; columns = file_name + one per channel; stain labels or NA
```

If `parse_panel_()` errors with "Column `S` doesn't exist", the fixture has no `$PnS` keywords. Add `values_fill = NA_character_` already present handles this — but also add a guard:

```r
# After the first pivot_wider, check if S column exists
if (!"S" %in% names(...)) { add stain = NA_character_ column }
```

Alternatively restructure as:

```r
panel_n <- ... |> dplyr::filter(type == "N") |> dplyr::select(number, channel = value)
panel_s <- ... |> dplyr::filter(type == "S") |> dplyr::select(number, stain  = value)
panel_long <- dplyr::left_join(panel_n, panel_s, by = "number") |>
  dplyr::mutate(stain = dplyr::na_if(stain, "")) |>
  dplyr::select(!number) |>
  tidyr::pivot_wider(names_from = channel, values_from = stain)
```

Use whichever form works with your fixture.

- [ ] **Step 5: Commit**

```bash
git add R/facs_read.R
git commit -m "feat: implement filter_samples_(), parse_keywords_(), parse_panel_()"
```

---

### Task 5: Implement walk_pops_() and parse_populations_()

**Files:**
- Modify: `R/facs_read.R`

**Interfaces:**
- Consumes: node/attribute names from Task 1; `filter_samples_()` from Task 4
- Produces:
  - `walk_pops_(node, file_name, path, parent_count, stain_lookup)` → long tibble `(file_name, population_full_path, population, metric, value)` for all descendant populations of `node`
  - `parse_populations_(doc, sample_ids)` → same schema across all samples

**Note:** The XPath strings `"Subpopulations/PopulationNode"` and `"Statistic"` below match standard FlowJo 10. If Task 1 showed different names, substitute them here. Similarly for attribute names `"count"`, `"name"`, `"channel"`, `"value"` on stat nodes.

- [ ] **Step 1: Replace walk_pops_() stub**

```r
walk_pops_ <- function(node, file_name, path, parent_count, stain_lookup) {
  pop_nodes <- xml2::xml_find_all(node, "Subpopulations/PopulationNode")
  if (length(pop_nodes) == 0L) return(tibble::tibble())

  purrr::map(pop_nodes, function(pop) {
    pop_name  <- xml2::xml_attr(pop, "name")
    pop_count <- suppressWarnings(as.numeric(xml2::xml_attr(pop, "count")))
    pop_path  <- if (nzchar(path)) paste0(path, "/", pop_name) else pop_name

    fop <- if (!is.na(parent_count) && parent_count > 0) {
      pop_count / parent_count
    } else {
      NA_real_
    }

    base_rows <- tibble::tibble(
      file_name            = file_name,
      population_full_path = pop_path,
      population           = pop_name,
      metric               = c("count", "fraction_of_parent"),
      value                = c(pop_count, fop)
    )

    # Statistic child nodes — skip count/frequency duplicates
    stat_nodes <- xml2::xml_find_all(pop, "Statistic")
    skip_stats <- c("count", "freq. of parent", "freq of parent",
                    "frequency of parent")

    stat_rows <- if (length(stat_nodes) > 0L) {
      purrr::map(stat_nodes, function(stat) {
        stat_type    <- xml2::xml_attr(stat, "name")
        stat_channel <- xml2::xml_attr(stat, "channel")
        stat_value   <- suppressWarnings(as.numeric(xml2::xml_attr(stat, "value")))

        if (is.na(stat_type) || tolower(stat_type) %in% skip_stats) {
          return(tibble::tibble())
        }
        if (is.na(stat_channel) || !nzchar(stat_channel)) {
          return(tibble::tibble())
        }

        matched <- stain_lookup$label[stain_lookup$channel == stat_channel]
        label   <- if (length(matched) > 0L && !is.na(matched[[1L]])) {
          matched[[1L]]
        } else {
          stat_channel
        }

        tibble::tibble(
          file_name            = file_name,
          population_full_path = pop_path,
          population           = pop_name,
          metric               = paste0(tolower(stat_type), "_", label),
          value                = stat_value
        )
      }) |>
        dplyr::bind_rows()
    } else {
      tibble::tibble()
    }

    dplyr::bind_rows(
      base_rows,
      stat_rows,
      walk_pops_(pop, file_name, pop_path, pop_count, stain_lookup)
    )
  }) |>
    dplyr::bind_rows()
}
```

- [ ] **Step 2: Replace parse_populations_() stub**

```r
parse_populations_ <- function(doc, sample_ids = NULL) {
  samples <- filter_samples_(
    xml2::xml_find_all(doc, ".//SampleList/Sample"),
    sample_ids
  )

  purrr::map(samples, function(sample) {
    file_name <- basename(
      xml2::xml_attr(xml2::xml_find_first(sample, "DataSet"), "uri")
    )

    # Build stain lookup: tibble(channel, label) for this sample
    kw_nodes <- xml2::xml_find_all(sample, ".//keywords/Keyword")
    stain_lookup <- tibble::tibble(
      key   = xml2::xml_attr(kw_nodes, "name"),
      value = xml2::xml_attr(kw_nodes, "value")
    ) |>
      dplyr::filter(grepl("^\\$P[0-9]+[NS]$", key)) |>
      dplyr::mutate(
        number = stringr::str_extract(key, "(?<=\\$P)[0-9]+"),
        type   = stringr::str_extract(key, "[NS]$")
      ) |>
      dplyr::select(!key) |>
      tidyr::pivot_wider(
        names_from  = type,
        values_from = value,
        values_fill = NA_character_
      ) |>
      dplyr::rename(channel = N, stain = S) |>
      dplyr::mutate(
        stain = dplyr::na_if(stain, ""),
        label = dplyr::coalesce(stain, channel)
      ) |>
      dplyr::select(channel, label)

    # Root tree node — name verified in Task 1 (standard: "SampleNode")
    root_node  <- xml2::xml_find_first(sample, "SampleNode")
    root_count <- suppressWarnings(as.numeric(xml2::xml_attr(root_node, "count")))

    walk_pops_(root_node, file_name, "", root_count, stain_lookup)
  }) |>
    dplyr::bind_rows()
}
```

- [ ] **Step 3: Sanity-check interactively**

```r
devtools::load_all()
doc  <- xml2::read_xml("tests/fixtures/minimal.wsp")
pops <- JanisHelpers:::parse_populations_(doc)
dplyr::glimpse(pops)
dplyr::count(pops, metric)
```

Expected: `metric` column contains `"count"`, `"fraction_of_parent"`, and at least one stat (e.g. `"median_CD4"`) if the fixture has exported statistics. `population_full_path` values should look like `"Lymphocytes/CD4+"`.

Troubleshooting:
- If `pops` is empty: the root node name is wrong. Try `xml2::xml_find_first(sample, "SampleNode")` vs other names from Task 1.
- If namespace errors occur: add `xml2::xml_ns_strip(doc)` immediately after `xml2::read_xml(path)` in `parse_populations_()` and in `facs_read_wsp()`.

- [ ] **Step 4: Commit**

```bash
git add R/facs_read.R
git commit -m "feat: implement parse_populations_() and walk_pops_() recursive traversal"
```

---

### Task 6: Implement facs_read_wsp() and pass all tests

**Files:**
- Modify: `R/facs_read.R`

**Interfaces:**
- Consumes: all helpers from Tasks 4–5
- Produces: exported `facs_read_wsp(path, group, keywords)` returning `list(data, meta, panel)`; all integration tests green

- [ ] **Step 1: Replace facs_read_wsp() stub body**

Keep the roxygen block unchanged; replace only the function body:

```r
facs_read_wsp <- function(path, group = NULL, keywords = NULL) {
  META_KEYS <- c("$DATE", "$BTIM", "$ETIM", "$CYT", "$INST", "$OP", "$TOT")

  doc <- xml2::read_xml(path)

  # Resolve group -> sample IDs for filtering
  sample_ids <- NULL
  if (!is.null(group)) {
    group_node <- xml2::xml_find_first(
      doc,
      glue::glue(".//Groups/GroupNode[@name='{group}']")
    )
    if (inherits(group_node, "xml_missing")) {
      available <- xml2::xml_attr(
        xml2::xml_find_all(doc, ".//Groups/GroupNode"),
        "name"
      )
      stop(glue::glue(
        "Group '{group}' not found in workspace. ",
        "Available groups: {paste(available, collapse = ', ')}"
      ))
    }
    sample_ids <- xml2::xml_attr(
      xml2::xml_find_all(group_node, ".//SampleRefs/SampleRef"),
      "sampleID"
    )
  }

  # Parse all components from the single open document
  kws  <- parse_keywords_(doc, sample_ids)
  pnl  <- parse_panel_(doc, sample_ids)
  pops <- parse_populations_(doc, sample_ids)

  # Build meta: system keywords, strip $ prefix, one row per file
  meta_long <- kws |>
    dplyr::filter(key %in% META_KEYS) |>
    dplyr::mutate(key = stringr::str_remove(key, "^\\$"))

  meta <- tidyr::pivot_wider(meta_long, names_from = key, values_from = value)

  missing_meta <- setdiff(stringr::str_remove(META_KEYS, "^\\$"), names(meta))
  if (length(missing_meta) > 0L) meta[missing_meta] <- NA_character_

  # Build data: population rows + optional keyword join
  data <- pops

  if (!is.null(keywords) && length(keywords) > 0L) {
    user_kws <- kws |>
      dplyr::filter(key %in% keywords) |>
      tidyr::pivot_wider(names_from = key, values_from = value)

    missing_kws <- setdiff(keywords, names(user_kws))
    if (length(missing_kws) > 0L) {
      warning(
        "The following requested keywords were not found in the workspace ",
        "and were filled with NA: ",
        paste(missing_kws, collapse = ", ")
      )
      user_kws[missing_kws] <- NA_character_
    }

    data <- dplyr::left_join(data, user_kws, by = "file_name")
  }

  n_files <- dplyr::n_distinct(data$file_name)
  message(glue::glue(
    "\nExtraction Summary",
    "\n----------------------------------------------",
    "\nExtracted groups:  {if (is.null(group)) 'all' else group}",
    "\nNumber of samples: {n_files}",
    "\n----------------------------------------------\n"
  ))

  list(data = data, meta = meta, panel = pnl)
}
```

- [ ] **Step 2: Run full test suite**

```r
devtools::load_all()
testthat::test_file("tests/testthat/test-facs_read.R")
```

Expected: all tests PASS (or SKIP if fixture absent).

- [ ] **Step 3: If tests fail, diagnose**

Common failure modes and fixes:

| Symptom | Fix |
|---|---|
| `xml_find_all` returns 0 nodes | Document has a default namespace. Add `xml2::xml_ns_strip(doc)` right after `xml2::read_xml(path)` in `facs_read_wsp()`. |
| `meta` has 0 rows | Fixture has no `$DATE` etc. keywords — expected; `missing_meta` guard fills columns with NA. Check `meta` has `file_name` column from `pivot_wider`. If `meta_long` is empty, `pivot_wider` produces a 0-row tibble — join `pops |> dplyr::distinct(file_name)` to it: `dplyr::left_join(dplyr::distinct(pops, file_name), meta, by = "file_name")`. |
| `panel` has duplicate channel columns | Two `$PnN` entries with same number. De-duplicate: add `dplyr::distinct(number, type, .keep_all = TRUE)` before `pivot_wider`. |
| `data` missing `file_name` after keyword join | `user_kws` is empty tibble (no matching keyword rows). The `missing_kws` guard adds the column manually. Verify `user_kws` has `file_name` before join. |

- [ ] **Step 4: Run devtools::check()**

```r
devtools::document()
devtools::check()
```

Expected: 0 errors, 0 warnings, ≤ 2 notes (pre-existing).

If a new NOTE appears for `facs_read.R`: likely a bare name inside a `dplyr` verb. Namespace it as `pkg::fn()` or use `.data$col` as appropriate.

- [ ] **Step 5: Commit**

```bash
git add R/facs_read.R man/facs_read_wsp.Rd
git commit -m "feat: implement facs_read_wsp() — xml2 WSP parser returning list(data, meta, panel)"
```

---

### Task 7: Final verification and commit

**Files:** None — verification only

**Interfaces:**
- Consumes: all previous tasks complete
- Produces: clean `devtools::check()` confirming 0 errors, 0 warnings

- [ ] **Step 1: Run full test suite across all files**

```r
devtools::load_all()
testthat::test_dir("tests/testthat/")
```

Expected: all tests PASS or SKIP; 0 failures.

- [ ] **Step 2: Run devtools::check() one final time**

```r
devtools::check()
```

Expected last lines:
```
0 errors ✔ | 0 warnings ✔ | 2 notes ✖
```

Any new note beyond the 2 pre-existing ones is a bug — do not commit until resolved.

- [ ] **Step 3: Final commit**

```bash
git add -p  # stage any remaining unstaged changes selectively
git commit -m "chore: final cleanup after facs_read_wsp implementation"
```
