# Export Boolean Gate Populations in facs_read_wsp() Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix `facs_read_wsp()` so that FlowJo boolean-gate populations (`NotNode`/`AndNode`/`OrNode` elements in the `.wsp` XML) — and every population gated further downstream of them — are included in the exported `data` tibble instead of being silently dropped.

**Architecture:** `walk_pops_()` in `R/facs_read.R` recurses through `<Subpopulations>` children and only recognizes the XML element name `"Population"` as an extractable/recursable population node; every other element name (including boolean-gate nodes) falls into a catch-all `else` branch that returns an empty tibble **and never recurses into it**, so the boolean population itself and its entire descendant subtree disappear. The fix generalizes the node-type check to include FlowJo's boolean-gate element names, which carry the same `name`/`count` attributes and the same optional `<Subpopulations>` child as a regular `<Population>` node, so no other logic needs to change.

**Tech Stack:** R, `xml2`, `testthat`, existing fixture `tests/fixtures/minimal.wsp` (already contains a real `NotNode` — no new fixture needed).

## Global Constraints

- Follow `CLAUDE.md` naming/style rules: base pipe `|>`, `pkg::fn()` namespacing, no bare tidyverse calls.
- `devtools::check()` must stay at 0 errors, 0 warnings after the change.
- No test mocking of `fcexpr`/`xml2` — use the real fixture file (per existing testing convention in `CLAUDE.md`).
- Minimal fix only: do not add new output columns (e.g. a "gate type" column) or parse `<Dependents>` — out of scope for this bug fix.

---

## Root cause (confirmed)

- `R/facs_read.R`, `walk_pops_()` (around lines 91–139), branches on `xml2::xml_name(child)`:
  - `"Population"` → extract `count`/`fraction_of_parent` rows and recurse.
  - `"Statistic"` → extract a stat row attributed to the current population.
  - anything else → `tibble::tibble()` (dropped, **no recursion**).
- FlowJo represents boolean-gate populations (AND/OR/NOT combinations) as sibling elements to `<Population>` inside the same `<Subpopulations>` container, using element names `<NotNode>`, `<AndNode>`, `<OrNode>` instead of `<Population>`. This was confirmed two ways:
  1. The `fcexpr` package (already an `Imports` dependency) hardcodes exactly this set of node names in its own wsp-walking helpers, e.g. `wsx_get_poppaths()`: `xml2::xml_name(prnts) %in% c("AndNode", "OrNode", "NotNode", "Population")`, and `recursive_walk_xml()` treats all of `SampleNode`/`Population`/`OrNode`/`AndNode`/`NotNode` uniformly when reading `name`/`count` attributes and recursing into `Subpopulations`.
  2. The existing fixture `tests/fixtures/minimal.wsp` already contains a real example at (approx.) line 2272, inside `.//SampleList/Sample[2]/SampleNode`:
     ```xml
     <NotNode name="non-debris-" ... count="378685" >
       ...
       <Subpopulations>
         <Statistic name="Median" ... id="Comp-PE-Cy5-A" ... value="15.401223183" />
         <Population name="CD45+" ... count="984" >...</Population>
       </Subpopulations>
       <Dependents>
         <Dependent name="Singlets/non-debris" />
       </Dependents>
     </NotNode>
     ```
     `NotNode` sits as a sibling of `<Population name="non-debris">` inside `Singlets`'s `<Subpopulations>`, has `name`/`count` attributes directly on it exactly like `Population`, and has its own `<Subpopulations>` with a further child population (`CD45+`) and a `<Statistic>`.
- Empirically verified by running `facs_read_wsp("tests/fixtures/minimal.wsp")`: `"Singlets/non-debris-"` and `"Singlets/non-debris-/CD45+"` are completely absent from `result$data$population_full_path`, even though the fixture contains real count data for both.

## Fix

Generalize the type check in `walk_pops_()` from `nm == "Population"` to `nm %in% POPULATION_NODE_TYPES` where `POPULATION_NODE_TYPES <- c("Population", "OrNode", "AndNode", "NotNode")`. No other code in that branch needs to change: `xml_attr(child, "name")`/`xml_attr(child, "count")` and the recursive `walk_pops_(child, ...)` call already work identically for these node types (confirmed above), because the recursion itself doesn't care what tag name the population lives in — it just calls `xml2::xml_find_first(node, "Subpopulations")` on whatever node it's given.

---

## Task 1: Export boolean-gate populations in `walk_pops_()`

**Files:**
- Modify: `R/facs_read.R:85-142` (`walk_pops_()`)
- Modify: `R/facs_read.R:1-18` (top-of-file XML structure notes comment block)
- Test: `tests/testthat/test-facs_read.R`

**Interfaces:**
- Consumes: existing fixture `tests/fixtures/minimal.wsp` (already checked in, already contains a `NotNode`).
- Produces: no new exported symbols; `facs_read_wsp()`'s return shape is unchanged, it just contains more rows.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-facs_read.R`:

```r
test_that("boolean gate populations (NotNode/AndNode/OrNode) are exported", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))

  boolean_pop <- result$data |>
    dplyr::filter(
      file_name == "26-1-17_percoll-kidney_E06.fcs",
      population_full_path == "Singlets/non-debris-"
    )

  expect_true(nrow(boolean_pop) > 0)
  expect_equal(
    boolean_pop$value[boolean_pop$metric == "count"],
    378685
  )
})

test_that("populations gated downstream of a boolean gate are exported", {
  skip_if_not(file.exists(wsp_path), skip_msg)
  result <- suppressMessages(facs_read_wsp(wsp_path))

  child_pop <- result$data |>
    dplyr::filter(
      file_name == "26-1-17_percoll-kidney_E06.fcs",
      population_full_path == "Singlets/non-debris-/CD45+"
    )

  expect_true(nrow(child_pop) > 0)
  expect_equal(
    child_pop$value[child_pop$metric == "count"],
    984
  )
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-facs_read.R")'`
Expected: FAIL — both new tests report `nrow(boolean_pop) > 0` / `nrow(child_pop) > 0` as `FALSE` (0 rows found), because the boolean population and its child are currently dropped entirely.

- [ ] **Step 3: Implement the minimal fix**

In `R/facs_read.R`, inside `walk_pops_()` (currently starts at line 85), add the node-type constant next to the existing `skip_stats` constant and change the branch condition:

```r
walk_pops_ <- function(node, file_name, path, parent_count, stain_lookup) {
  subpops <- xml2::xml_find_first(node, "Subpopulations")
  if (inherits(subpops, "xml_missing")) return(tibble::tibble())

  skip_stats <- c("count", "freq. of parent", "freq of parent", "frequency of parent")
  # FlowJo represents boolean-gate populations (AND/OR/NOT of other gates) as
  # sibling element types to <Population>, not as <Population> itself, but
  # they carry the same name/count attributes and Subpopulations structure.
  pop_node_types <- c("Population", "OrNode", "AndNode", "NotNode")

  purrr::map(xml2::xml_children(subpops), function(child) {
    nm <- xml2::xml_name(child)

    if (nm %in% pop_node_types) {
      pop_name  <- xml2::xml_attr(child, "name")
      pop_count <- suppressWarnings(as.numeric(xml2::xml_attr(child, "count")))
      pop_path  <- if (nzchar(path)) paste0(path, "/", pop_name) else pop_name

      fop <- if (!is.na(parent_count) && parent_count > 0L) {
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

      dplyr::bind_rows(
        base_rows,
        walk_pops_(child, file_name, pop_path, pop_count, stain_lookup)
      )

    } else if (nm == "Statistic") {
      stat_type    <- xml2::xml_attr(child, "name")
      stat_channel <- xml2::xml_attr(child, "id")
      stat_value   <- suppressWarnings(as.numeric(xml2::xml_attr(child, "value")))

      if (is.na(stat_type) || tolower(stat_type) %in% skip_stats) return(tibble::tibble())
      if (is.na(stat_channel) || !nzchar(stat_channel))           return(tibble::tibble())

      matched <- stain_lookup$label[stain_lookup$channel == stat_channel]
      label   <- if (length(matched) > 0L && !is.na(matched[[1L]])) matched[[1L]] else stat_channel

      tibble::tibble(
        file_name            = file_name,
        population_full_path = path,
        population           = basename(path),
        metric               = paste0(tolower(stat_type), "_", label),
        value                = stat_value
      )

    } else {
      tibble::tibble()
    }
  }) |>
    dplyr::bind_rows()
}
```

Also update the top-of-file XML structure notes comment block (`R/facs_read.R:1-18`) — add one line documenting the boolean-gate element names, right after the existing "Individual population element" line:

```
# Individual population element:                            Population
# Boolean-gate population elements (same attrs/children
# as Population; AND/OR/NOT combinations of other gates):   OrNode, AndNode, NotNode
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-facs_read.R")'`
Expected: PASS — all tests in the file green, including the two new ones.

- [ ] **Step 5: Run full check**

Run: `Rscript -e 'devtools::check(quiet = TRUE)'`
Expected: 0 errors, 0 warnings (same pre-existing notes as documented in `CLAUDE.md`, nothing new).

- [ ] **Step 6: Commit**

```bash
git add R/facs_read.R tests/testthat/test-facs_read.R
git commit -m "fix: export boolean-gate populations (NotNode/AndNode/OrNode) in facs_read_wsp()"
```

---

## Self-review notes

- Spec coverage: the only requirement ("boolean gates are not exported — find out why, plan a proper fix") is covered by Task 1: root cause documented above, fix generalizes node-type matching, regression tests cover both the boolean population itself and a population gated downstream of it (the two ways the old code silently lost data).
- No placeholders: all steps include literal code/commands.
- Type consistency: no new functions/signatures introduced; `walk_pops_()`'s signature and return shape (a tibble with `file_name`, `population_full_path`, `population`, `metric`, `value`) are unchanged.
- Explicitly out of scope (call out if the user wants it later): recording which populations are ORed/ANDed/NOTed together (`<Dependents>`) as extra columns, and distinguishing boolean-derived rows from regular gate rows via a `gate_type` column. Neither was requested; both would be additive, non-breaking follow-ups if wanted.
